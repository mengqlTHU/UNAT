!-------------------------------------------------------------------------------
!   set the volume and face for FV.
!-------------------------------------------------------------------------------
    subroutine set_vol_face(lev)
    use var_kind_def
    use var_cgns
    use var_global, only: n_dim,R13,translation_axis,rotation_axis,err_mem
    use var_mesh
    use var_parallel
    implicit none
    integer(dpI),intent(in):: lev
    integer(dpI):: isec,ele,npe,iele,i,v(8),im,sL,eL,sR,eR,sL_g,sR_g
    real   (dpR):: xyz(3,8),xyz_L(3,8),xyz_R(3,8),c(3),d(3),cL(3),cR(3),n(3),area, &
                &  vL,ome(3),rtmp

    call mesh_set_face_neighbour(lev)

!   ----------------------------------------------------------------------------
!   get the volume of the element.
    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(.not. allocated(sec(isec)%vol))  allocate(sec(isec)%vol(sec(isec)%n_ele))
        if(sec(isec)%is_bnd)    cycle
        ele =  sec(isec)%ele_type
        npe =  sec(isec)%npe

        do iele=1,sec(isec)%n_ele
            v(1:npe)=  sec(isec)%n2e(1:npe,iele)
            forall(i=1:npe) xyz(1:n_dim,i)  =  mesh(lev)%xyz(1:n_dim,v(i))

            if(ele .eq. TRI_3) then
                call tri_get_area(xyz, rtmp)
            elseif(ele .eq. QUAD_4) then
                call quad_get_area(xyz, rtmp)
            elseif(ele .eq. TETRA_4) then
                call tetra_get_vol(xyz(1,1), xyz(1,2), xyz(1,3), xyz(1,4), rtmp)
            elseif(ele .eq. PYRA_5) then
                call pyra_get_vol(xyz(1,1), xyz(1,2), xyz(1,3), xyz(1,4), xyz(1,5), rtmp)
            elseif(ele .eq. PENTA_6) then
                call penta_get_vol(xyz, rtmp)
            elseif(ele .eq. HEXA_8) then
                call hexa_get_vol(xyz, rtmp)
            else
                stop 'Error: element type not supported.'
            end if
            if(rtmp .le. 0.0d0) print*,isec,iele,rtmp
            if(rtmp .le. 0.0d0) stop 'Error: non-positive volume.'
            sec(isec)%vol(iele) =  rtmp
        end do
    end do
!   get the volume of the element.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   set the mortar.
    allocate(mesh(lev)%mortar_cen (3, mesh(lev)%n_mortar), stat=err_mem)
    allocate(mesh(lev)%mortar_n_vg(8, mesh(lev)%n_mortar), stat=err_mem)
    mesh(lev)%mortar_cen    =  0.0d0
    mesh(lev)%mortar_n_vg   =  0.0d0

    do im=1,mesh(lev)%n_mortar
        sL  =  mesh(lev)%mortar_LR(1,im)
        eL  =  mesh(lev)%mortar_LR(2,im)
        sR  =  mesh(lev)%mortar_LR(3,im)
        eR  =  mesh(lev)%mortar_LR(4,im)
        ele =  mesh(lev)%mortar_ele_type(im)
        npe =  npe_ele(ele)
        v(1:npe)=  mesh(lev)%mortar_n2e(1:npe,im)
        forall(i=1:npe) xyz(1:n_dim,i)  =  mesh(lev)%xyz(1:n_dim,v(i))
        n   =  0.0d0

        if(ele .eq. BAR_2) then
            c(1:2)  = (xyz(1:2,1)+xyz(1:2,2))*0.5d0
            n(1:2)  = (/xyz(2,2)-xyz(2,1),-xyz(1,2)+xyz(1,1)/)
        elseif(ele .eq. TRI_3) then
            c(1:3)  = (xyz(1:3,1)+xyz(1:3,2)+xyz(1:3,3))*R13
            call crs_prd(xyz(1:3,2)-xyz(1:3,1), xyz(1:3,3)-xyz(1:3,1), n)
            n       =  0.5d0*n
        elseif(ele .eq. QUAD_4) then
            c(1:3)  = (xyz(1:3,1)+xyz(1:3,2)+xyz(1:3,3)+xyz(1:3,4))*0.25d0
            call crs_prd(xyz(1:3,3)-xyz(1:3,1), xyz(1:3,4)-xyz(1:3,2), n)
            n       =  0.5d0*n
        else
            stop 'Error: mortar type not supported.'
        end if
        call norm_vec(n_dim, n, area)

        mesh(lev)%mortar_cen(1:3    ,im)    =  0.0d0
        mesh(lev)%mortar_cen(1:n_dim,im)    =  c(1:n_dim)

!       ------------------------------------------------------------------------
!       specify the grid velocity.
        sL_g=  sec(sL)%ID_sec_g
        sR_g=  sec(sR)%ID_sec_g
        if(im .le. mesh(lev)%n_mortar_b) then
        else
            if((sec_motion_type(sL_g) .ne. sec_motion_type(sR_g)) .or. &
              &(abs(sec_motion_speed(sL_g)-sec_motion_speed(sR_g)) .ge. 1.0d-10)) &
                &  stop 'Error: speed of L&R section not match.'
        end if

        if(sec_motion_type(sL_g) .eq. 0) then
            mesh(lev)%mortar_n_vg(6:8,im)   =  0.0d0
        elseif(sec_motion_type(sL_g) .eq. 1) then
            mesh(lev)%mortar_n_vg(6:8,im)   =  sec_motion_speed(sL_g)*translation_axis
        elseif(sec_motion_type(sL_g) .eq. 2) then
            ome =  sec_motion_speed(sL_g)*rotation_axis
            mesh(lev)%mortar_n_vg(6,im) =  ome(2)*c(3)-ome(3)*c(2)
            mesh(lev)%mortar_n_vg(7,im) =  ome(3)*c(1)-ome(1)*c(3)
            mesh(lev)%mortar_n_vg(8,im) =  ome(1)*c(2)-ome(2)*c(1)
        end if
!       specify the grid velocity.
!       ------------------------------------------------------------------------

        npe =  npe_ele(sec(sL)%ele_type)
        forall(i=1:npe) xyz_L(1:n_dim,i)=  mesh(lev)%xyz(1:n_dim,sec(sL)%n2e(i,eL))
        cL(1:n_dim) =  xyz_L(1:n_dim,1)
        do i=2,npe
            cL(1:n_dim) =  cL(1:n_dim)+xyz_L(1:n_dim,i)
        end do
        cL  =  cL/real(npe, dpR)
        vL  =  sec(sL)%vol(eL)

        npe =  npe_ele(sec(sR)%ele_type)
        forall(i=1:npe) xyz_R(1:n_dim,i)=  mesh(lev)%xyz(1:n_dim,sec(sR)%n2e(i,eR))
        cR(1:n_dim) =  xyz_R(1:n_dim,1)
        do i=2,npe
            cR(1:n_dim) =  cR(1:n_dim)+xyz_R(1:n_dim,i)
        end do
        cR  =  cR/real(npe, dpR)

        d   =  cR-cL
        call norm_vec(n_dim, d, rtmp)
        rtmp=  n(1)*d(1)
        do i=2,n_dim
            rtmp=  rtmp+n(i)*d(i)
        end do
        if(rtmp .gt. 0.0d0) then
!           do nothing.
        elseif(rtmp .lt. 0.0d0) then
            if(ele .eq. BAR_2) then
                mesh(lev)%mortar_n2e(1:2,im)= (/v(2), v(1)/)
            elseif(ele .eq. TRI_3) then
                mesh(lev)%mortar_n2e(1:3,im)= (/v(1), v(3), v(2)/)
            else
                mesh(lev)%mortar_n2e(1:4,im)= (/v(1), v(4), v(3), v(2)/)
            end if
            n   = -n
        else
            stop 'Error: normal of face perpendicular to cL-->cR.'
        end if

        mesh(lev)%mortar_n_vg(1:3,im)   =  n(1:3)
        mesh(lev)%mortar_n_vg(4  ,im)   =  area
        rtmp= (cR(1)-cL(1))*n(1)+(cR(2)-cL(2))*n(2)+(cR(3)-cL(3))*n(3)
        if(.true.) then
            if(im .le. mesh(lev)%n_mortar_b) then
                mesh(lev)%mortar_n_vg(5,im) =  0.5d0/abs(rtmp)
            else
                mesh(lev)%mortar_n_vg(5,im) =  1.0d0/abs(rtmp)
            end if
        else
            if(im .le. mesh(lev)%n_mortar_b) then
                mesh(lev)%mortar_n_vg(5,im) =        area/ vL
            else
                mesh(lev)%mortar_n_vg(5,im) =  2.0d0*area/(vL+sec(sR)%vol(eR))
            end if
        end if
    end do
!   set the mortar.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   check the closeness of internal element.
    call mesh_check(lev)
!   check the closeness of internal element.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   compute the center.
    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        allocate(sec(isec)%cen(3, sec(isec)%n_ele), stat=err_mem)
        npe             =  sec(isec)%npe
        sec(isec)%cen   =  0.0d0

        do iele=1,sec(isec)%n_ele
            v(1:npe)=  sec(isec)%n2e(1:npe,iele)
            sec(isec)%cen(1:n_dim,iele) =  mesh(lev)%xyz(1:n_dim,v(1))
            do i=2,npe
                sec(isec)%cen(1:n_dim,iele) =  sec(isec)%cen(1:n_dim,iele) &
                                            & +mesh(lev)%xyz(1:n_dim,v(i))
            end do
            sec(isec)%cen(1:n_dim,iele) =  sec(isec)%cen(1:n_dim,iele)/real(npe, dpR)
        end do
    end do
!   compute the center.
!   ----------------------------------------------------------------------------

    call set_ls_gra(lev)

!   ----------------------------------------------------------------------------
!   set the velocity of the boundary section.
    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(.not. sec(isec)%is_bnd)  cycle
        allocate(sec(isec)%vg(3, sec(isec)%n_ele), stat=err_mem)
        sL_g=  sec(sL)%ID_sec_g

        if(sec_motion_type(sL_g) .eq. 0) then
            sec(isec)%vg=  0.0d0
        elseif(sec_motion_type(sL_g) .eq. 1) then
            do iele=1,sec(isec)%n_ele
                sec(isec)%vg(1:3,iele)  =  sec_motion_speed(sL_g)*translation_axis
            end do
        elseif(sec_motion_type(sL_g) .eq. 1) then
            ome =  sec_motion_speed(sL_g)*rotation_axis
            do iele=1,sec(isec)%n_ele
                c(1:3)  =  sec(isec)%cen(1:3,iele)
                sec(isec)%vg(1,iele)=  ome(2)*c(3)-ome(3)*c(2)
                sec(isec)%vg(2,iele)=  ome(3)*c(1)-ome(1)*c(3)
                sec(isec)%vg(3,iele)=  ome(1)*c(2)-ome(2)*c(1)
            end do
        else
            stop 'Error: motion type not supported.'
        end if
    end do
!   set the velocity of the boundary section.
!   ----------------------------------------------------------------------------

    return
    end subroutine set_vol_face
!-------------------------------------------------------------------------------
!   get the area of triangle.
!-------------------------------------------------------------------------------
    subroutine tri_get_area(xyz,area)
    use var_kind_def
    implicit none
    real(dpR),intent(in):: xyz(3,*)
    real(dpR):: area,d1(2),d2(2)

    d1(1:2) =  xyz(1:2,2)-xyz(1:2,1)
    d2(1:2) =  xyz(1:2,3)-xyz(1:2,1)
    area    =  0.5d0*(d1(1)*d2(2)-d2(1)*d1(2))

    return
    end subroutine tri_get_area
!-------------------------------------------------------------------------------
!   get the area of quad.
!-------------------------------------------------------------------------------
    subroutine quad_get_area(xyz,area)
    use var_kind_def
    implicit none
    real(dpR),intent(in):: xyz(3,*)
    real(dpR):: area,d1(2),d2(2)

    d1(1:2) =  xyz(1:2,3)-xyz(1:2,1)
    d2(1:2) =  xyz(1:2,4)-xyz(1:2,2)
    area    =  0.5d0*(d1(1)*d2(2)-d2(1)*d1(2))

    return
    end subroutine quad_get_area
!-------------------------------------------------------------------------------
!   get the vol of tetra.
!-------------------------------------------------------------------------------
    subroutine tetra_get_vol(v1,v2,v3,v4,vol)
    use var_kind_def
    implicit none
    real(dpR),intent(in):: v1(*),v2(*),v3(*),v4(*)
    real(dpR):: vol,d1(3),d2(3),n(3)

    d1(1:3) =  v2(1:3)-v1(1:3)
    d2(1:3) =  v3(1:3)-v1(1:3)
    call crs_prd(d1, d2, n)
    vol = (n(1)*(v4(1)-v1(1))+n(2)*(v4(2)-v1(2))+n(3)*(v4(3)-v1(3)))/6.0d0

    return
    end subroutine tetra_get_vol
!-------------------------------------------------------------------------------
!   get the vol of pyra.
!-------------------------------------------------------------------------------
    subroutine pyra_get_vol(v1,v2,v3,v4,v5,vol)
    use var_kind_def
    implicit none
    real(dpR),intent(in):: v1(*),v2(*),v3(*),v4(*),v5(*)
    real(dpR):: vol,d1(3),d2(3),n(3)

    d1(1:3) =  v3(1:3)-v1(1:3)
    d2(1:3) =  v4(1:3)-v2(1:3)
    call crs_prd(d1, d2, n)
    vol = (n(1)*(v5(1)-v1(1))+n(2)*(v5(2)-v1(2))+n(3)*(v5(3)-v1(3)))/6.0d0

    return
    end subroutine pyra_get_vol
!-------------------------------------------------------------------------------
!   get the vol of penta.
!-------------------------------------------------------------------------------
    subroutine penta_get_vol(xyz,vol)
    use var_kind_def
    implicit none
    real(dpR),intent(in):: xyz(3,*)
    real(dpR):: vol,c(3),rtmp

    c(1:3)  = (xyz(1:3,1)+xyz(1:3,2)+xyz(1:3,3)+xyz(1:3,4)+xyz(1:3,5)+xyz(1:3,6))/6.0d0

    call tetra_get_vol(xyz(1,1), xyz(1,2), xyz(1,3), c, rtmp)
    vol     =  rtmp
    call tetra_get_vol(xyz(1,4), xyz(1,6), xyz(1,5), c, rtmp)
    vol     =  vol+rtmp

    call pyra_get_vol(xyz(1,1), xyz(1,4), xyz(1,5), xyz(1,2), c, rtmp)
    vol     =  vol+rtmp
    call pyra_get_vol(xyz(1,2), xyz(1,5), xyz(1,6), xyz(1,3), c, rtmp)
    vol     =  vol+rtmp
    call pyra_get_vol(xyz(1,1), xyz(1,3), xyz(1,6), xyz(1,4), c, rtmp)
    vol     =  vol+rtmp

    return
    end subroutine penta_get_vol
!-------------------------------------------------------------------------------
!   get the vol of hexa.
!-------------------------------------------------------------------------------
    subroutine hexa_get_vol(xyz,vol)
    use var_kind_def
    implicit none
    real(dpR),intent(in):: xyz(3,*)
    real(dpR):: vol,c(3),rtmp

    c(1:3)  = (xyz(1:3,1)+xyz(1:3,2)+xyz(1:3,3)+xyz(1:3,4) &
            & +xyz(1:3,5)+xyz(1:3,6)+xyz(1:3,7)+xyz(1:3,8))*0.125d0

    call pyra_get_vol(xyz(1,1), xyz(1,2), xyz(1,3), xyz(1,4), c, rtmp)
    vol     =  rtmp
    call pyra_get_vol(xyz(1,5), xyz(1,8), xyz(1,7), xyz(1,6), c, rtmp)
    vol     =  vol+rtmp
    call pyra_get_vol(xyz(1,1), xyz(1,4), xyz(1,8), xyz(1,5), c, rtmp)
    vol     =  vol+rtmp
    call pyra_get_vol(xyz(1,2), xyz(1,6), xyz(1,7), xyz(1,3), c, rtmp)
    vol     =  vol+rtmp
    call pyra_get_vol(xyz(1,1), xyz(1,5), xyz(1,6), xyz(1,2), c, rtmp)
    vol     =  vol+rtmp
    call pyra_get_vol(xyz(1,3), xyz(1,7), xyz(1,8), xyz(1,4), c, rtmp)
    vol     =  vol+rtmp

    return
    end subroutine hexa_get_vol
!-------------------------------------------------------------------------------
!   get the ID of ghost element, if its donor beglongs to myid.
!-------------------------------------------------------------------------------
    subroutine get_ID_ghost_ele(isec,iele,s,e,per)
    use var_kind_def
    use var_mesh
    use var_parallel, only: myid
    implicit none
    integer(dpI),intent(in):: isec,iele
    integer(dpI):: s,e,per(*),isr,i

    if(.not. sec(isec)%is_ghost)    stop 'Error: fails to get ID of ghost element.'
    s   =  0
    e   =  0
    isr =  sec(isec)%id_recv(1,iele)
    i   =  sec(isec)%id_recv(2,iele)
    if(p2p(isr)%ip_remote .ne. myid)    return
    per(1:3)=  sec(isec)%per_path(1:3,iele)
    s   =  p2p(isr)%id_ele_send(1,i)
    e   =  p2p(isr)%id_ele_send(2,i)

    return
    end subroutine get_ID_ghost_ele
!-------------------------------------------------------------------------------
!   check the mesh.
!-------------------------------------------------------------------------------
    subroutine mesh_check(lev)
    use var_kind_def
    use var_global, only: err_mem,L_ref
    use var_mesh
    implicit none
    integer(dpI),intent(in):: lev
    integer(dpI):: nele,isec,iele,im,sL,eL,sR,eR
    real   (dpR):: n(3),rtmp
    integer(dpI),allocatable:: iA(:)
    real   (dpR),allocatable:: r(:,:)

    allocate(iA(mesh(lev)%sec_1:mesh(lev)%sec_0+1), stat=err_mem)
    iA(mesh(lev)%sec_1) =  1
    nele                =  0
    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(sec(isec)%is_int)    nele=  nele+sec(isec)%n_ele
        iA(isec+1)  =  nele+1
    end do
    allocate(r(3,nele), stat=err_mem)
    r   =  0.0d0

    do im=1,mesh(lev)%n_mortar
        sL  =  mesh(lev)%mortar_LR(1,im)
        eL  =  mesh(lev)%mortar_LR(2,im)
        sR  =  mesh(lev)%mortar_LR(3,im)
        eR  =  mesh(lev)%mortar_LR(4,im)
        n(1:3)  =  mesh(lev)%mortar_n_vg(1:3,im)*mesh(lev)%mortar_n_vg(4,im)

        r(1:3,eL+iA(sL)-1)  =  r(1:3,eL+iA(sL)-1)+n(1:3)
        if(sec(sR)%is_int) then
        r(1:3,eR+iA(sR)-1)  =  r(1:3,eR+iA(sR)-1)-n(1:3)
        end if
    end do
    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(.not. sec(isec)%is_int)  cycle
        do iele=1,sec(isec)%n_ele
            rtmp=  maxval(abs(r(1:3,iele+iA(isec)-1)))
            if(rtmp .ge. 1.0d-13*L_ref) stop 'Error: element is not closed.'
        end do
    end do

    if(allocated(iA))   deallocate(iA)
    if(allocated(r ))   deallocate(r )

    return
    end subroutine mesh_check
!-------------------------------------------------------------------------------
!   set face neighbour.
!-------------------------------------------------------------------------------
    subroutine mesh_set_face_neighbour(lev)
    use var_kind_def
    use var_global, only: err_mem
    use var_mesh
    implicit none
    integer(dpI),intent(in):: lev
    logical(dpL):: ltmp
    integer(dpI):: isec,iele,i,j,im,sL,eL,sR,eR,ifac

    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(.not. sec(isec)%is_int)  cycle

        allocate(sec(isec)%iA_face_neighbour(sec(isec)%n_ele+1), stat=err_mem)
        sec(isec)%iA_face_neighbour =  0
    end do
    do im=1,mesh(lev)%n_mortar
        sL  =  mesh(lev)%mortar_LR(1,im)
        eL  =  mesh(lev)%mortar_LR(2,im)
        sR  =  mesh(lev)%mortar_LR(3,im)
        eR  =  mesh(lev)%mortar_LR(4,im)

        sec(sL)%iA_face_neighbour(eL+1) =  sec(sL)%iA_face_neighbour(eL+1)+1
        if(sec(sR)%is_int) then
        sec(sR)%iA_face_neighbour(eR+1) =  sec(sR)%iA_face_neighbour(eR+1)+1
        end if
    end do
    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(.not. sec(isec)%is_int)  cycle
        i   =  0
        do iele=1,sec(isec)%n_ele
            i   =  i+sec(isec)%iA_face_neighbour(iele+1)
        end do
        allocate(sec(isec)%jA_face_neighbour(i), stat=err_mem)
        sec(isec)%iA_face_neighbour(1)  =  1
        sec(isec)%jA_face_neighbour     =  0
        do iele=1,sec(isec)%n_ele
            sec(isec)%iA_face_neighbour(iele+1) =  sec(isec)%iA_face_neighbour(iele+1) &
                & +sec(isec)%iA_face_neighbour(iele)
        end do
        do iele=1,sec(isec)%n_ele
            i   =  sec(isec)%iA_face_neighbour(iele+1)-1
            sec(isec)%jA_face_neighbour(i)  =  sec(isec)%iA_face_neighbour(iele)
        end do
    end do
    do im=1,mesh(lev)%n_mortar
        sL  =  mesh(lev)%mortar_LR(1,im)
        eL  =  mesh(lev)%mortar_LR(2,im)
        sR  =  mesh(lev)%mortar_LR(3,im)
        eR  =  mesh(lev)%mortar_LR(4,im)

        i   =  sec(sL)%iA_face_neighbour(eL+1)-1
        j   =  sec(sL)%jA_face_neighbour(i)
        sec(sL)%jA_face_neighbour(j)=  im
        if(j .lt. i)    sec(sL)%jA_face_neighbour(i)=  sec(sL)%jA_face_neighbour(i)+1

        if(sec(sR)%is_int) then
        i   =  sec(sR)%iA_face_neighbour(eR+1)-1
        j   =  sec(sR)%jA_face_neighbour(i)
        sec(sR)%jA_face_neighbour(j)= -im
        if(j .lt. i)    sec(sR)%jA_face_neighbour(i)=  sec(sR)%jA_face_neighbour(i)+1
        end if
    end do
    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(.not. sec(isec)%is_int)  cycle
        do iele=1,sec(isec)%n_ele
        do ifac=sec(isec)%iA_face_neighbour(iele),sec(isec)%iA_face_neighbour(iele+1)-1
            im  =  sec(isec)%jA_face_neighbour(ifac)
            if(im .gt. 0) then
                ltmp= (mesh(lev)%mortar_LR(1, im) .eq. isec) .and. &
                    & (mesh(lev)%mortar_LR(2, im) .eq. iele)
            else
                ltmp= (mesh(lev)%mortar_LR(3,-im) .eq. isec) .and. &
                    & (mesh(lev)%mortar_LR(4,-im) .eq. iele)
            end if
            if(.not. ltmp)  stop 'Error: fails to set face_neighbour.'
        end do
        end do
    end do

    return
    end subroutine mesh_set_face_neighbour
!-------------------------------------------------------------------------------
!   interpolate cell-center variables to vertex.
!-------------------------------------------------------------------------------
    subroutine mesh_get_vc_weights
    use var_kind_def
    use var_global, only: is_hanging_node,err_mem,n_dim
    use var_mesh
    implicit none
    integer(dpI):: isec,iele,i,j,k,ivtx,im,imh,sL,eL,sR,eR,idx,st(2,100),n_st,iA(101)
    real   (dpR):: d(3),w(100),rtmp
    integer(dpI),allocatable:: e2n(:,:)

    allocate(mesh(0)%iA_e2n(mesh(0)%n_vtx+1), stat=err_mem)
    if(err_mem .ne. 0)  stop 'Error: fv_get_vc fails to allocate memory.'
    mesh(0)%iA_e2n  =  0

!   ----------------------------------------------------------------------------
!   cell-to-vertex defined by n2e.
    do isec=mesh(0)%sec_1,mesh(0)%sec_0
        if(sec(isec)%is_bnd)    cycle

        do iele=1,sec(isec)%n_ele
        do i=1,sec(isec)%npe
            ivtx=  sec(isec)%n2e(i,iele)
            if(ivtx .le. mesh(0)%n_vtx) mesh(0)%iA_e2n(ivtx)=  mesh(0)%iA_e2n(ivtx)+1
        end do
        end do
    end do

!   cell-to-vertex defined by interface, hanging node.
    if(is_hanging_node) then
        do imh=1,mesh(0)%n_mortar_hanging
            im  =  mesh(0)%mortar_hanging(imh)
            sL  =  mesh(0)%mortar_LR(1,im)
            eL  =  mesh(0)%mortar_LR(2,im)
            sR  =  mesh(0)%mortar_LR(3,im)
            eR  =  mesh(0)%mortar_LR(4,im)
            do i=1,2*n_dim-2
                ivtx=  mesh(0)%mortar_n2e(i,im)
                if(ivtx .le. 0) exit
                if(ivtx .le. mesh(0)%n_vtx) mesh(0)%iA_e2n(ivtx)=  mesh(0)%iA_e2n(ivtx)+2
            end do
        end do
    end if
!   ----------------------------------------------------------------------------

    i   =  0
    do ivtx=1,mesh(0)%n_vtx
        i   =  i+mesh(0)%iA_e2n(ivtx)
    end do
    allocate(e2n(3,i), stat=err_mem)
    if(err_mem .ne. 0)  stop 'Error: me_get_vc_weights fails to allocate memory.'
    e2n =  0

!   ----------------------------------------------------------------------------
!   record the cell-to-vertex defined by n2e.
    idx =  0
    do isec=mesh(0)%sec_1,mesh(0)%sec_0
        if(sec(isec)%is_bnd)    cycle

        do iele=1,sec(isec)%n_ele
        do i=1,sec(isec)%npe
            ivtx=  sec(isec)%n2e(i,iele)
            if(ivtx .gt. mesh(0)%n_vtx) cycle
            idx =  idx+1
            e2n(1,idx)  =  ivtx
            e2n(2,idx)  =  isec
            e2n(3,idx)  =  iele
        end do
        end do
    end do

!   record the cell-to-vertex defined by interface, hanging node.
    if(is_hanging_node) then
        do imh=1,mesh(0)%n_mortar_hanging
            im  =  mesh(0)%mortar_hanging(imh)
            sL  =  mesh(0)%mortar_LR(1,im)
            eL  =  mesh(0)%mortar_LR(2,im)
            sR  =  mesh(0)%mortar_LR(3,im)
            eR  =  mesh(0)%mortar_LR(4,im)
            do i=1,2*n_dim-2
                ivtx=  mesh(0)%mortar_n2e(i,im)
                if(ivtx .le. 0) exit
                if(ivtx .gt. mesh(0)%n_vtx) cycle
                idx =  idx+1
                e2n(1,idx)  =  ivtx
                e2n(2,idx)  =  sL
                e2n(3,idx)  =  eL
                idx =  idx+1
                e2n(1,idx)  =  ivtx
                e2n(2,idx)  =  sR
                e2n(3,idx)  =  eR
            end do
        end do
    end if
!   ----------------------------------------------------------------------------

    call iqsortcols(.true., 1, idx, 1, 3, e2n)
    i   =  1
    sL  =  0
    eL  =  0
    do while(i .le. idx)
        k   =  i
        ivtx=  e2n(1,i)
        eL  =  eL+1
        if(eL .ne. ivtx)    stop 'Error: fails to get vc_weights.'
        do j=i+1,idx
            if(e2n(1,j) .ne. ivtx) then
                exit
            else
                k   =  j
            end if
        end do

        do j=i,k
            st(1:2,j-i+1) =  e2n(2:3,j)
        end do
        n_st=  k-i+1
        if(is_hanging_node) call simplify_matrix(2, n_st, 1, 2, st)
        do j=1,n_st
            sL  =  sL+1
            e2n(1  ,sL) =  ivtx
            e2n(2:3,sL) =  st(1:2,j)
        end do

        i   =  k+1
    end do
    allocate(mesh(0)%jA_e2n(sL), stat=err_mem)
    allocate(mesh(0)%w_e2n (sL), stat=err_mem)
    mesh(0)%iA_e2n      =  0
    mesh(0)%iA_e2n(1)   =  1
    iA(1)=  1
    do isec=mesh(0)%sec_1,mesh(0)%sec_0
        iA(isec+1)  =  iA(isec)+sec(isec)%n_ele
    end do

    i   =  1
    do while(i .le. sL)
        k   =  i
        ivtx=  e2n(1,i)
        do j=i+1,sL
            if(e2n(1,j) .ne. ivtx) then
                exit
            else
                k   =  j
            end if
        end do

        rtmp=  0.0d0
        do j=i,k
            isec                =  e2n(2,j)
            iele                =  e2n(3,j)
            mesh(0)%jA_e2n(j)   =  iA(isec)+iele-1

            if(n_dim .eq. 2) then
                d(1:2)  =  sec(isec)%cen(1:2,iele)-mesh(0)%xyz(1:2,ivtx)
                w(j-i+1)=  1.0d0/sqrt(d(1)**2+d(2)**2        )
            else
                d(1:3)  =  sec(isec)%cen(1:3,iele)-mesh(0)%xyz(1:3,ivtx)
                w(j-i+1)=  1.0d0/sqrt(d(1)**2+d(2)**2+d(3)**2)
            end if
            rtmp    =  rtmp+w(j-i+1)
        end do
        w(1:k-i+1)  =  w(1:k-i+1)/rtmp
        call DCOPY(k-i+1, w, 1, mesh(0)%w_e2n(i), 1)

        i   =  k+1
        mesh(0)%iA_e2n(ivtx+1)  =  i
    end do

    if(allocated(e2n))  deallocate(e2n)

    return
    end subroutine mesh_get_vc_weights
