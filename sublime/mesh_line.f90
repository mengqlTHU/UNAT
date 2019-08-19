!-------------------------------------------------------------------------------
!   module 1D element lines.
!-------------------------------------------------------------------------------
    module var_eline
        use var_kind_def
        implicit none

        type type_eline
            integer(dpI):: n_ele    =  0
            integer(dpI):: mortar_L =  0
            integer(dpI):: mortar_R =  0
            integer(dpI),allocatable:: ele(:,:)
            real   (dpR),allocatable:: d  (:  )

            real   (dpR):: BL_thickness =  0.0d0
            real   (dpR):: intermittency=  0.0d0
        end type type_eline
        type(type_eline),allocatable:: eline(:)
        integer(dpI):: n_eline=  0
        contains
!       ------------------------------------------------------------------------
!       output the 1D lines.
!       ------------------------------------------------------------------------
        subroutine eline_output
        use var_global, only: is_2d_cal,n_dim
        use var_mesh
        use var_parallel
        implicit none
        character(len=80):: str
        integer(dpI):: i,j,sL,sR,isec,iele

        write(str,*),myid
        str =  './data/implicit_lines_'//trim(adjustl(str))//'.dat'
        open(unit=10,file=trim(str))
        if(is_2d_cal) then
            write(unit=10,fmt=*),'variables="CoordinateX","CoordinateY"'
        else
            write(unit=10,fmt=*),'variables="CoordinateX","CoordinateY","CoordinateZ"'
        end if

        do i=1,n_eline
            if(eline(i)%intermittency .le. 0.0d0)   cycle
            write(str,*),i
            write(unit=10,fmt=*),'zone T="'//trim(adjustl(str))//'"'
            sL  =  eline(i)%mortar_L
            sR  =  eline(i)%mortar_R
            if(sL .gt. 0)   write(unit=10,fmt=*),mesh(0)%mortar_cen(1:n_dim,sL)
            do j=1,eline(i)%n_ele
                isec=  eline(i)%ele(1,j)
                iele=  eline(i)%ele(2,j)
                write(unit=10,fmt=*),sec(isec)%cen(1:n_dim,iele)
            end do
            if(sR .gt. 0)   write(unit=10,fmt=*),mesh(0)%mortar_cen(1:n_dim,sR)
        end do
        close(10)

        return
        end subroutine eline_output
    end module var_eline
!-------------------------------------------------------------------------------
!   get implicit lines.
!-------------------------------------------------------------------------------
    subroutine mesh_get_1d_lines
    use var_kind_def
    use var_cgns
    use var_eline, n_line=>n_eline
    use var_global, only: err_mem,BL_thick
    use var_mesh
    implicit none
    logical(dpL):: is_output_line,ltmp
    integer(dpI):: iA(1000),n_ele,isec,iele,i,j,k,im,n_ele_line,sL,eL,sR,face,line(2,1000)
    real   (dpR):: max_thickness,A_max,A_min,A(6),min_AR
    integer(dpI),allocatable:: ID(:),lines(:,:),iA_line(:),jA_line(:),LR_line(:,:)
    real   (dpR),allocatable:: AR(:),AR_face(:)

    is_output_line  =  .false.
    min_AR          =  2.0d1
    max_thickness   =  2.0d0*BL_thick

    if(allocated(eline))    return

!   ----------------------------------------------------------------------------
!   memory allocation.
    n_ele   =  0
    iA(1)   =  1
    do isec=mesh(0)%sec_1,mesh(0)%sec_0
        if(.not. sec(isec)%is_int)  cycle
        n_ele       =  n_ele+sec(isec)%n_ele
        iA(isec+1)  =  n_ele+1
    end do
    allocate(ID(n_ele), stat=err_mem)
    allocate(AR(n_ele), stat=err_mem)
    allocate(AR_face(mesh(0)%n_mortar), stat=err_mem)
!   memory allocation.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   get the AR of element and mortar.
    AR      =  0.0d0
    AR_face =  0.0d0
    do isec=mesh(0)%sec_1,mesh(0)%sec_0
        if(.not. is_quad_ele(isec)) cycle
        do iele=1,sec(isec)%n_ele
            A_max   =  0.0d0
            A_min   =  huge(1.0d0)
            do i=sec(isec)%iA_face_neighbour(iele),sec(isec)%iA_face_neighbour(iele+1)-1
                im      =  sec(isec)%jA_face_neighbour(i)
                k       =  i-sec(isec)%iA_face_neighbour(iele)+1
                A(k)    =  mesh(0)%mortar_n_vg(4,abs(im))
                A_max   =  max(A_max, A(k))
                A_min   =  min(A_min, A(k))
            end do

            AR(iA(isec)+iele-1) =  A_max/A_min
            do i=sec(isec)%iA_face_neighbour(iele),sec(isec)%iA_face_neighbour(iele+1)-1
                im  =  abs(sec(isec)%jA_face_neighbour(i))
                AR_face(im) =  max(AR_face(im), mesh(0)%mortar_n_vg(4,abs(im))/A_min)
            end do
        end do
    end do
!   get the AR of element and mortar.
!   ----------------------------------------------------------------------------

    ID      =  0
    n_line  =  0

!   ----------------------------------------------------------------------------
!   construct lines from the wall elements.
    do im=1,mesh(0)%n_mortar_b
        sL  =  mesh(0)%mortar_LR(1,im)
        eL  =  mesh(0)%mortar_LR(2,im)
        sR  =  mesh(0)%mortar_LR(3,im)
        if(sec(sR)%bct .ne. BCWallViscous)  cycle

!       only QUAD, PENTA and HEXA considered.
        if(.not. is_quad_ele(sL))   cycle

!       this element has already been indexed.
        if(ID(iA(sL)+eL-1) .gt. 0)  cycle

        n_line  =  n_line+1
        ID(iA(sL)+eL-1) =  n_line
        face            =  im
        n_ele_line      =  1
        line(1:2,1)     = (/sL, eL/)

        call get_line(.true., iA, n_line, 5, ID, n_ele_line, face, sL, eL)
    end do
!   construct lines from the wall elements.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   construct lines from the internal cells.
    do isec=mesh(0)%sec_1,mesh(0)%sec_0
        if(.not. sec(isec)%is_int)  cycle
        if(.not. is_quad_ele(isec)) cycle

        do iele=1,sec(isec)%n_ele
!           this element has been indexed.
            if(ID(iA(isec)+iele-1) .gt. 0)  cycle

!           do we have available face?
            ltmp=  .false.
            do i=sec(isec)%iA_face_neighbour(iele),sec(isec)%iA_face_neighbour(iele+1)-1
                im  =  abs(sec(isec)%jA_face_neighbour(i))
                ltmp=  ltmp .or. (AR_face(im) .ge. min_AR)
                if(ltmp)    exit
            end do
            if(.not. ltmp)  cycle

            n_line  =  n_line+1
            ID(iA(isec)+iele-1) =  n_line
            n_ele_line          =  1
            line(1:2,1)         = (/isec, iele/)

            do i=sec(isec)%iA_face_neighbour(iele),sec(isec)%iA_face_neighbour(iele+1)-1
                im  =  abs(sec(isec)%jA_face_neighbour(i))
                if(AR_face(im) .lt. min_AR) cycle

                face=  abs(sec(isec)%jA_face_neighbour(i))
                sL  =  isec
                eL  =  iele
                call get_line(.false., iA, n_line, 0, ID, n_ele_line, face, sL, eL)
            end do
        end do
    end do
!   construct lines from the internal cells.
!   ----------------------------------------------------------------------------

    if(allocated(AR     ))  deallocate(AR)
    if(allocated(AR_face))  deallocate(AR_face)

    i   =  0
    do iele=1,n_ele
        if(ID(iele) .gt. 0) i   =  i+1
    end do
    allocate(lines(3,i), stat=err_mem)
    allocate(iA_line(n_line+1), stat=err_mem)
    allocate(jA_line(i       ), stat=err_mem)
    allocate(LR_line(2,n_line), stat=err_mem)

    lines       =  0
    n_ele_line  =  0
    do iele=1,n_ele
        if(ID(iele) .le. 0) cycle

        n_ele_line  =  n_ele_line+1
        do isec=mesh(0)%sec_1,mesh(0)%sec_0
            if((iele .ge. iA(isec)) .and. (iele .lt. iA(isec+1))) then
                lines(1,n_ele_line) =  ID(iele)
                lines(2,n_ele_line) =  isec
                lines(3,n_ele_line) =  iele-iA(isec)+1
                exit
            end if
        end do
        if(lines(1,n_ele_line) .le. 0)  &
            &  stop 'Error: fails to get the section of element, 1d line.'
    end do
    call iqsortcols(.true., 1, n_ele_line, 1, 3, lines)

    n_line      =  0
    i           =  1
    iA_line(1)  =  1
    do while(i .le. n_ele_line)
        k   =  i
        do j=i+1,n_ele_line
            if(lines(1,j) .eq. lines(1,i)) then
                k   =  j
            else
                exit
            end if
        end do

!       this line is too short.
        if(k-i .le. 2) then
            i   =  k+1
            cycle
        end if

        n_line  =  n_line+1
        forall(j=i:k)   lines(1,j)  =  n_line
        call get_order(k-i+1, lines(1,i), LR_line(1,n_line), LR_line(2,n_line))
!       print*,n_line,':',LR_line(1:2,n_line)
        iA_line(n_line+1)   =  iA_line(n_line)+k-i+1
        do j=i,k
            im  =  iA_line(n_line)+j-i
            isec=  lines(2,j)
            iele=  lines(3,j)
            jA_line(im) =  iA(isec)+iele-1
        end do

        i   =  k+1
    end do

    allocate(eline(n_line), stat=err_mem)

    do i=1,n_line
        eline(i)%n_ele  =  iA_line(i+1)-iA_line(i)
        if(LR_line(1,i) .gt. 0) eline(i)%mortar_L   =  LR_line(1,i)
        if(LR_line(2,i) .gt. 0) eline(i)%mortar_R   =  LR_line(2,i)
        allocate(eline(i)%ele(2, eline(i)%n_ele), stat=err_mem)
        eline(i)%ele    =  0
        do j=iA_line(i),iA_line(i+1)-1
            k   =  j-iA_line(i)+1
            do isec=mesh(0)%sec_1,mesh(0)%sec_0
                if(.not. sec(isec)%is_int)  cycle
                if(jA_line(j) .lt. iA(isec  ))  cycle
                if(jA_line(j) .ge. iA(isec+1))  cycle

                eline(i)%ele(1,k)   =  isec
                eline(i)%ele(2,k)   =  jA_line(j)-iA(isec)+1
                exit
            end do
            if(eline(i)%ele(1,k) .le. 0)    stop 'Error: mesh_get_1d_line fails.'
        end do
    end do

    if(is_output_line)  call eline_output

    if(allocated(ID     ))  deallocate(ID     )
    if(allocated(lines  ))  deallocate(lines  )
    if(allocated(iA_line))  deallocate(iA_line)
    if(allocated(jA_line))  deallocate(jA_line)
    if(allocated(LR_line))  deallocate(LR_line)
    if(allocated(AR     ))  deallocate(AR     )
    if(allocated(AR_face))  deallocate(AR_face)

    return
    contains
!       ------------------------------------------------------------------------
!       element n2e minus face n2e.
!       ------------------------------------------------------------------------
        subroutine get_opposite_face(n_e,e,n_f,f,n_r,r)
        use var_kind_def
        implicit none
        integer(dpI),intent(in):: n_e,e(*),n_f,f(*)
        logical(dpl):: ltmp
        integer(dpI):: n_r,r(*),i,LL,RR

        call iqsort(.true., 1, n_e, e)
        r(1:n_e)=  e(1:n_e)
        do i=1,n_f
            call ib_search(1, n_e, 1, 1, e, f(i), ltmp, LL, RR)
            if((.not. ltmp) .or. (LL .ne. RR))  stop 'Error: face not found in element.'
            r(LL)   = -1
        end do
        n_r =  0
        do i=1,n_e
            if(r(i) .le. 0) cycle
            n_r     =  n_r+1
            r(n_r)  =  r(i)
        end do

        return
        end subroutine get_opposite_face
!       ------------------------------------------------------------------------
!       get mortar form vertex.
!       ------------------------------------------------------------------------
        subroutine get_mortar(isec,iele,npe,v,im,s_nb,e_nb)
        use var_kind_def
        use var_mesh
        implicit none
        integer(dpI),intent(in):: isec,iele,npe,v(*)
        logical(dpL):: ltmp
        integer(dpI):: im,s_nb,e_nb,i,j,f(4)

        do i=sec(isec)%iA_face_neighbour(iele),sec(isec)%iA_face_neighbour(iele+1)-1
            im      =  sec(isec)%jA_face_neighbour(i)
            f(1:npe)=  mesh(0)%mortar_n2e(1:npe,abs(im))
            call iqsort(.true., 1, npe, f)
            ltmp=  .true.
            do j=1,npe
                ltmp=  v(j) .eq. f(j)
                if(.not. ltmp)  exit
            end do
            if(ltmp) then
                if(im .gt. 0) then
                    s_nb=  mesh(0)%mortar_LR(3,im)
                    e_nb=  mesh(0)%mortar_LR(4,im)
                else
                    s_nb=  mesh(0)%mortar_LR(1,-im)
                    e_nb=  mesh(0)%mortar_LR(2,-im)
                end if
                im  =  abs(im)
                return
            end if
        end do
        stop 'Error: fails to get mortar from n2e.'

        return
        end subroutine get_mortar
!       ------------------------------------------------------------------------
!       get the line.
!       ------------------------------------------------------------------------
        subroutine get_line(is_from_wall,iA,n_line,min_n_ele,ID,n_ele_line,face,sL,eL)
        use var_kind_def
        use var_mesh
        implicit none
        logical(dpL),intent(in):: is_from_wall
        integer(dpI),intent(in):: iA(*),n_line,min_n_ele
        logical(dpL):: ltmp
        integer(dpI):: ID(*),face,sL,eL,n_ele_line,n_f,v_f(4),n_e,v_e(4),n_r, &
                    &  v_r(4),im_nb,s_nb,e_nb,id_nb

        do while(.true.)
            n_f         =  npe_ele(mesh(0)%mortar_ele_type(face))
            v_f(1:n_f)  =  mesh(0)%mortar_n2e(1:n_f,face)
            n_e         =  sec(sL)%npe
            v_e(1:n_e)  =  sec(sL)%n2e(1:n_e,eL)
            call get_opposite_face(n_e, v_e, n_f, v_f, n_r, v_r)

            call get_mortar(sL, eL, n_r, v_r, im_nb, s_nb, e_nb)

!           this element is not QUAD.
            if(.not. is_quad_ele(s_nb)) exit

            id_nb   =  iA(s_nb)+e_nb-1

!           this element has already been indexed.
            if(ID(id_nb) .gt. 0)    exit

            if(is_from_wall) then
                ltmp=  sec(s_nb)%dnw(e_nb) .gt. max_thickness
            else
!               AR of this element is small; or AR of this face is small.
                ltmp= ((n_ele_line .ge. min_n_ele) .and. (AR(id_nb) .le. min_AR)) .or. &
                    & ((n_ele_line .ge. min_n_ele) .and. (AR_face(im_nb) .le. min_AR))
            end if
            if(ltmp)    exit

!           add this element to the line.
            n_ele_line  =  n_ele_line+1
            ID(id_nb)   =  n_line
            line(1:2,n_ele_line)= (/s_nb, e_nb/)

!           to the next one.
            face=  im_nb
            sL  =  s_nb
            eL  =  e_nb
        end do

        return
        end subroutine get_line
!       ------------------------------------------------------------------------
!       QUAD, PENTA, HEXA.
!       ------------------------------------------------------------------------
        function is_quad_ele(isec) result(is_qph)
        use var_kind_def
        use var_cgns
        use var_mesh, only: sec
        implicit none
        integer(dpI),intent(in):: isec
        logical(dpL):: is_qph

        is_qph  =  sec(isec)%is_quad .or. sec(isec)%is_hexa .or. &
                & (sec(isec)%ele_type .eq. PENTA_6)

        return
        end function is_quad_ele
!       ------------------------------------------------------------------------
!       order the line.
!       ------------------------------------------------------------------------
        subroutine get_order(n_ele,line,L,R)
        use var_kind_def
        use var_mesh, only: mesh,sec
        implicit none
        integer(dpI),intent(in):: n_ele
        logical(dpL):: ltmv(n_ele)
        integer(dpI):: line(3,*),ID,nb(3,n_ele),im,i,j,k,isec,iele,order(n_ele), &
                    &  N,seed,L,R,s_nb,e_nb,v(3),v1(8),v2(8)
        logical(dpL),external:: list_eq_list

        ID  =  line(1,1)
        nb  =  0
        do i=1,n_ele
            isec=  line(2,i)
            iele=  line(3,i)

!           find neighbours.
            do k=sec(isec)%iA_face_neighbour(iele),sec(isec)%iA_face_neighbour(iele+1)-1
                im  =  sec(isec)%jA_face_neighbour(k)
                if(im .gt. 0) then
                    s_nb=  mesh(0)%mortar_LR(3, im)
                    e_nb=  mesh(0)%mortar_LR(4, im)
                else
                    s_nb=  mesh(0)%mortar_LR(1,-im)
                    e_nb=  mesh(0)%mortar_LR(2,-im)
                end if

                do j=1,n_ele
                    if((s_nb .eq. line(2,j)) .and. (e_nb .eq. line(3,j))) then
                        nb(1,i) =  nb(1,i)+1
                        nb(nb(1,i)+1,i) =  j
                        exit
                    end if
                end do
            end do
        end do

!       check whether the line is valid.
        j   =  0
        k   =  0
        do i=1,n_ele
            if(nb(1,i) .eq. 1) then
                j   =  j+1
            elseif(nb(1,i) .eq. 2) then
                k   =  k+1
            end if
        end do
        if((j .ne. 2) .or. (j+k .ne. n_ele))    stop 'Error: fails to order the line.'

        ltmv=  .true.

!       add the first two elements to the line.
        do i=1,n_ele
            if(nb(1,i) .eq. 1) then
                order(1)=  i
                order(2)=  nb(2,i)

                ltmv(i)         =  .false.
                ltmv(nb(2,i))   =  .false.
                exit
            end if
        end do

        seed=  order(2)
        N   =  2
        do while(.true.)
            L   =  nb(2,seed)
            R   =  nb(3,seed)

!           seed is the end element.
            if(R .le. 0)    exit

            if(ltmv(L) .and. ltmv(R))   stop 'Error: fails to order the line.'

            if(.not. (ltmv(L) .or. ltmv(R))) then
!               L and R are not available.
                stop 'Error: fails to order the line.'
            elseif(ltmv(L) .and. ltmv(R)) then
!               L and R are both available.
                stop 'Error: fails to order the line.'
            elseif(ltmv(L)) then
                N   =  N+1
                order(N)=  L
                ltmv (L)=  .false.
            else
                N   =  N+1
                order(N)    =  R
                ltmv (R)=  .false.
            end if
            if(N .eq. n_ele)    exit
            seed=  order(N)
        end do

        do i=1,n_ele
            line(1,order(i))=  i
        end do
        call iqsortcols(.true., 1, n_ele, 1, 3, line)
        line(1,1:n_ele) =  ID

!       get left wall surface.
        L       =  0
        isec    =  line(2,1)
        iele    =  line(3,1)
        j       =  line(2,2)
        k       =  line(3,2)
        i       =  sec(isec)%npe
        v1(1:i) =  sec(isec)%n2e(1:i,iele)
        call list_minus_list(i, v1, sec(j)%npe, sec(j)%n2e(1,k))
        do k=sec(isec)%iA_face_neighbour(iele),sec(isec)%iA_face_neighbour(iele+1)-1
            im  =  sec(isec)%jA_face_neighbour(k)
            if(im .gt. 0) then
                s_nb=  mesh(0)%mortar_LR(3, im)
                e_nb=  mesh(0)%mortar_LR(4, im) 
            else
                s_nb=  mesh(0)%mortar_LR(1,-im)
                e_nb=  mesh(0)%mortar_LR(2,-im)
            end if
            if((sec(s_nb)%bct .ne. BCWallViscous) .or. (sec(s_nb)%npe .ne. i))  cycle
            v2(1:i) =  sec(s_nb)%n2e(1:i,e_nb)
            call iqsort(.true., 1, i, v2)
            if(list_eq_list(i, v1, v2)) then
                L   =  abs(im)
                exit
            end if
        end do

!       get right wall surface.
        R       =  0
        isec    =  line(2,n_ele  )
        iele    =  line(3,n_ele  )
        j       =  line(2,n_ele-1)
        k       =  line(3,n_ele-1)
        i       =  sec(isec)%npe
        v1(1:i) =  sec(isec)%n2e(1:i,iele)
        call list_minus_list(i, v1, sec(j)%npe, sec(j)%n2e(1,k))
        do k=sec(isec)%iA_face_neighbour(iele),sec(isec)%iA_face_neighbour(iele+1)-1
            im  =  sec(isec)%jA_face_neighbour(k)
            if(im .gt. 0) then
                s_nb=  mesh(0)%mortar_LR(3, im)
                e_nb=  mesh(0)%mortar_LR(4, im) 
            else
                s_nb=  mesh(0)%mortar_LR(1,-im)
                e_nb=  mesh(0)%mortar_LR(2,-im)
            end if
            if((sec(s_nb)%bct .ne. BCWallViscous) .or. (sec(s_nb)%npe .ne. i))  cycle
            v2(1:i) =  sec(s_nb)%n2e(1:i,e_nb)
            call iqsort(.true., 1, i, v2)
            if(list_eq_list(i, v1, v2)) then
                R   =  abs(im)
                exit
            end if
        end do

!       at last, reverse the ordering if the element-beyond-wall is not the first one.
        if((L .le. 0) .and. (R .gt. 0)) then
            do k=1,n_ele/2
                v   (1:3  ) =  line(1:3,k)
                line(1:3,k) =  line(1:3,n_ele+1-k)
                line(1:3,n_ele+1-k) =  v(1:3)
            end do
            L   =  R
            R   =  0
        end if

        return
        end subroutine get_order
    end subroutine mesh_get_1d_lines
