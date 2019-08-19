    include './include/rtree.f90'
!-------------------------------------------------------------------------------
!   definition of the solid wall.
!-------------------------------------------------------------------------------
    module var_wall
        use var_kind_def
        use var_rtree
        implicit none

        integer(dpI):: n_vtx_wall   =  0
        integer(dpI):: n_ele_wall   =  0
        type(type_rtree),save:: tree
        integer(dpI),allocatable:: iA_wall(:),jA_wall(:),ID_wall(:,:)
        real   (dpR),allocatable:: xyz_wall(:,:)
        contains
!       ------------------------------------------------------------------------
!       output the solid wall if necessary.
!       ------------------------------------------------------------------------
        subroutine output_wall
        implicit none
        character(len=100):: str,tec
        integer(dpI):: i,j,v(3),isec,n_sec

        open(unit=10,file='./data/solid_wall.dat')
        write(unit=10,fmt='(A)'),'variables="x","y","z"'

        v   =  0
        do i=1,n_ele_wall
            j   =  iA_wall(i+1)-iA_wall(i)
            if(j .eq. 2) then
                v(1)=  v(1)+1
            elseif(j .eq. 3) then
                v(2)=  v(2)+1
            elseif(j .eq. 4) then
                v(3)=  v(3)+1
            else
                stop 'Error: element type not supported for wall.'
            end if
        end do

        n_sec   =  0
        do isec=1,3
            if(v(isec) .le. 0)  cycle
            n_sec   =  n_sec+1

            tec =  'zone N='
            write(str,*),size(xyz_wall, dim=2)
            tec =  trim(adjustl(tec))//trim(adjustl(str))
            write(str,*),v(isec)
            tec =  trim(adjustl(tec))//',E='//trim(adjustl(str))
            if(isec .eq. 1) then
                str =  'FELINESEG'
            elseif(isec .eq. 2) then
                str =  'FETRIANGLE'
            else
                str =  'FEQUADRILATERAL'
            end if
            tec =  trim(adjustl(tec))//','//'zonetype='//trim(adjustl(str))
            if(n_sec .gt. 1) then
                write(str,*),3
                tec =  trim(adjustl(tec))//',varsharelist=([1-'//trim(adjustl(str)) &
                    &  //']=1)'
            end if
            tec =  trim(adjustl(tec))//',datapacking=point'
            write(unit=10,fmt='(A)'),trim(adjustl(tec))
            if(isec .eq. 1) then
                do i=1,size(xyz_wall, dim=2)
                    write(unit=10,fmt='(3ES20.12)'),xyz_wall(1:3,i)
                end do
            end if

            do i=1,n_ele_wall
                j   =  iA_wall(i+1)-iA_wall(i)
                if((isec .eq. 1) .and. (j .eq. 2)) then
                    write(unit=10,fmt='(2I8)'),jA_wall(iA_wall(i):iA_wall(i+1)-1)
                elseif((isec .eq. 2) .and. (j .eq. 3)) then
                    write(unit=10,fmt='(3I8)'),jA_wall(iA_wall(i):iA_wall(i+1)-1)
                elseif((isec .eq. 3) .and. (j .eq. 4)) then
                    write(unit=10,fmt='(4I8)'),jA_wall(iA_wall(i):iA_wall(i+1)-1)
                end if
            end do
        end do

        close(10)

        return
        end subroutine output_wall
    end module var_wall
!-------------------------------------------------------------------------------
!   output vertex based dnw.
!-------------------------------------------------------------------------------
    subroutine wr_vertex_dnw_local(dnw)
    use var_kind_def
    use var_global, only: mesh_name
    use var_mesh
    use var_parallel, only: myid
    implicit none
    real   (dpR),intent(in):: dnw(*)
    character(len=1000):: str,tec
    integer(dpI):: npe,n_vtx,isec,n_ele,i,iele

    n_vtx   =  mesh(0)%n_vtx
    n_ele   =  0
    do isec=mesh(0)%sec_1,mesh(0)%sec_0
        if(sec(isec)%is_int)    n_ele   =  n_ele+sec(isec)%n_ele
    end do

    write(str,*),myid
    str =  trim(adjustl(mesh_name))//'_sol_'//trim(adjustl(str))//'.dat'
    open(unit=10,file=trim(adjustl(str)))
    write(unit=10,fmt='(A)'),'variables="x","y","d"'

    do isec=mesh(0)%sec_1,mesh(0)%sec_0
        if(sec(isec)%is_ghost .or. sec(isec)%is_bnd)    cycle

        tec =  'zone N='
        write(str,*),n_vtx
        tec =  trim(adjustl(tec))//trim(adjustl(str))
        write(str,*),sec(isec)%n_ele
        tec =  trim(adjustl(tec))//',E='//trim(adjustl(str))
        if(sec(isec)%is_quad) then
            str =  'FEQUADRILATERAL'
        else
            str =  'FETRIANGLE'
        end if
        tec =  trim(adjustl(tec))//','//'zonetype='//trim(adjustl(str))
        if(isec .gt. 1) then
            write(str,*),3
            tec =  trim(adjustl(tec))//',varsharelist=([1-'//trim(adjustl(str)) &
                &  //']=1)'
        end if
        tec =  trim(adjustl(tec))//',datapacking=point'
        write(unit=10,fmt='(A)'),trim(adjustl(tec))
        if(isec .eq. 1) then
            do i=1,n_vtx
                write(unit=10,fmt='(3ES20.12)'),mesh(0)%xyz(1:2,i),dnw(i)
            end do
        end if
        npe =  sec(isec)%npe
        do iele=1,sec(isec)%n_ele
            write(unit=10,fmt='(8I8)'),sec(isec)%n2e(1:npe,iele)
        end do
    end do

    return
    end subroutine wr_vertex_dnw_local
!-------------------------------------------------------------------------------
!   get distance to a BAR element.
!-------------------------------------------------------------------------------
    subroutine get_dis_to_bar(p,p1,p2,is_inner,d)
    use var_kind_def
    implicit none
    real   (dpR),intent(in):: p(*),p1(*),p2(*)
    logical(dpL):: is_inner
    real   (dpR):: d,a,b,t

    a   = (p2(1)-p1(1))*(p2(1)-p1(1))+(p2(2)-p1(2))*(p2(2)-p1(2)) &
        &+(p2(3)-p1(3))*(p2(3)-p1(3))
    b   = (p (1)-p1(1))*(p2(1)-p1(1))+(p (2)-p1(2))*(p2(2)-p1(2)) &
        &+(p (3)-p1(3))*(p2(3)-p1(3))
    t   =  b/a
    is_inner= (t .ge. 0.0d0) .and. (t .le. 1.0d0)
    call norm_vec(3, (1.0d0-t)*p1(1:3)+t*p2(1:3)-p(1:3), d)

    return
    end subroutine get_dis_to_bar
!-------------------------------------------------------------------------------
!   get distance to a QUAD element.
!-------------------------------------------------------------------------------
    subroutine get_dis_to_quad(p,p1,p2,p3,p4,is_inner,d,n)
    use var_kind_def
    implicit none
    real   (dpR),intent(in):: p(*),p1(*),p2(*),p3(*),p4(*)
    logical(dpL):: is_inner
    real   (dpR):: n(*),x(3),y(3),c(3),LHS(3,3),d,q(2,4),O(2),a,b

    x(1:3)  =  p3(1:3)-p1(1:3)
    y(1:3)  =  p4(1:3)-p2(1:3)
    c(1:3)  = (p1(1:3)+p2(1:3)+p3(1:3)+p4(1:3))*0.25d0
    n(1)    =  x(2)*y(3)-x(3)*y(2)
    n(2)    = -x(1)*y(3)+x(3)*y(1)
    n(3)    =  x(1)*y(2)-x(2)*y(1)
    n(1:3)  =  n(1:3)/sqrt(n(1)*n(1)+n(2)*n(2)+n(3)*n(3))
    LHS(1:3,1)  =  x(1:3)
    LHS(1:3,2)  =  y(1:3)
    LHS(1:3,3)  =  n(1:3)
    call mat_inv(3, LHS)
    d   =  abs(n(1)*(p(1)-c(1))+n(2)*(p(2)-c(2))+n(3)*(p(3)-c(3)))

    q(1:2,1)=  LHS(1:2,1)*(p1(1)-c(1))+LHS(1:2,2)*(p1(2)-c(2))+LHS(1:2,3)*(p1(3)-c(3))
    q(1:2,2)=  LHS(1:2,1)*(p2(1)-c(1))+LHS(1:2,2)*(p2(2)-c(2))+LHS(1:2,3)*(p2(3)-c(3))
    q(1:2,3)=  LHS(1:2,1)*(p3(1)-c(1))+LHS(1:2,2)*(p3(2)-c(2))+LHS(1:2,3)*(p3(3)-c(3))
    q(1:2,4)=  LHS(1:2,1)*(p4(1)-c(1))+LHS(1:2,2)*(p4(2)-c(2))+LHS(1:2,3)*(p4(3)-c(3))
    O(1:2)  =  (/q(1,2), q(2,1)/)
    n(1:2)  =  LHS(1:2,1)*(p(1)-c(1))+LHS(1:2,2)*(p(2)-c(2)) &
            & +LHS(1:2,3)*(p(3)-c(3))-O(1:2)
    a       =  q(1,3)-O(1)
    b       =  q(2,4)-O(2)
    is_inner=  .false.
    if(a*n(2)+ b       *n(1) .ge. a* b       )  return
    if(a*n(2)+(b-1.0d0)*n(1) .le. a*(b-1.0d0))  return
    if((a-1.0d0)*n(2)+ b       *n(1) .le. (a-1.0d0)* b       )  return
    if((a-1.0d0)*n(2)+(b-1.0d0)*n(1) .ge. (a-1.0d0)*(b-1.0d0))  return
    is_inner=  .true.

    return
    end subroutine get_dis_to_quad
!-------------------------------------------------------------------------------
!   get distance to a TRI element.
!-------------------------------------------------------------------------------
    subroutine get_dis_to_tri(p,p1,p2,p3,is_inner,d,n)
    use var_kind_def
    implicit none
    real   (dpR),intent(in):: p(*),p1(*),p2(*),p3(*)
    logical(dpL):: is_inner
    real   (dpR):: d,n(*),x(3),y(3),LHS(3,3),rst(3)

    x(1:3)  =  p2(1:3)-p1(1:3)
    y(1:3)  =  p3(1:3)-p1(1:3)
    n(1)    =  x(2)*y(3)-x(3)*y(2)
    n(2)    = -x(1)*y(3)+x(3)*y(1)
    n(3)    =  x(1)*y(2)-x(2)*y(1)
    n(1:3)  =  n(1:3)/sqrt(n(1)*n(1)+n(2)*n(2)+n(3)*n(3))
    LHS(1:3,1)  =  x(1:3)
    LHS(1:3,2)  =  y(1:3)
    LHS(1:3,3)  =  n(1:3)
    rst(1:3  )  =  p(1:3)
    call slv_AxB(3, 1, LHS, rst)
    d       =  abs(rst(3))
    is_inner= (rst(1) .ge. 0.0d0) .and. (rst(2) .ge. 0.0d0) .and. &
            & (rst(1)+rst(2) .le. 1.0d0)

    return
    end subroutine get_dis_to_tri
!-------------------------------------------------------------------------------
!   get wall element definition and setup rtree.
!-------------------------------------------------------------------------------
    subroutine mesh_get_wall
    use var_kind_def
    use var_cgns, only: BCWallViscous
    use var_global, only: is_2d_cal,err_mem,n_dim,pi
    use var_mesh
    use var_parallel
    use var_slv, only: is_vis_cal
    use var_wall
    implicit none
    logical(dpL):: ltmp
    integer(dpI):: n_vtx_L,n_ele_L,im,sR,eR,i,j,k,v(4),iper,nper,isec_g,L,R, &
                &  LR(2,4),npe,pivot
    real   (dpR):: d(3),rm(3,3),rtmp
    integer(dpI),allocatable:: iA_L(:),jA_L(:),vtx_L(:,:),ID_L(:,:)
    real   (dpR),allocatable:: xyz_L(:,:)

    if(.not. is_vis_cal)    return

!   ----------------------------------------------------------------------------
!   number the vertex used to define the wall elements.
    n_vtx_L =  0
    n_ele_L =  0
    do im=1,mesh(0)%n_mortar_b
        sR  =  mesh(0)%mortar_LR(3,im)
        if(sec(sR)%bct .ne. BCWallViscous)  cycle
        isec_g  =  sec(sR)%ID_sec_g

        nper=  0
        do iper=1,n_sec_per_info
            if(nint(sec_per_info(1,iper), dpI) .ne. isec_g) cycle
            nper=  nper+1
        end do
        n_ele_L =  n_ele_L+1+2*nper
        n_vtx_L =  n_vtx_L+(1+2*nper)*sec(sR)%npe
    end do
    allocate(iA_L (  n_ele_L), stat=err_mem)
    allocate(ID_L (2,n_ele_L), stat=err_mem)
    allocate(jA_L (  n_vtx_L), stat=err_mem)
    allocate(vtx_L(2,n_vtx_L), stat=err_mem)

    n_vtx_L =  0
    do im=1,mesh(0)%n_mortar_b
        sR  =  mesh(0)%mortar_LR(3,im)
        eR  =  mesh(0)%mortar_LR(4,im)
        if(sec(sR)%bct .ne. BCWallViscous)  cycle
        isec_g  =  sec(sR)%ID_sec_g
        npe     =  sec(sR)%npe

        do i=1,npe
            n_vtx_L =  n_vtx_L+1
            vtx_L(1,n_vtx_L)=  sec(sR)%n2e(i,eR)
            vtx_L(2,n_vtx_L)=  0
        end do

        do iper=1,n_sec_per_info
            if(nint(sec_per_info(1,iper), dpI) .ne. isec_g) cycle

            do i=1,npe
                n_vtx_L         =  n_vtx_L+1
                vtx_L(1,n_vtx_L)=  sec(sR)%n2e(i,eR)
                vtx_L(2,n_vtx_L)=  iper

                n_vtx_L         =  n_vtx_L+1
                vtx_L(1,n_vtx_L)=  sec(sR)%n2e(i,eR)
                vtx_L(2,n_vtx_L)= -iper
            end do
        end do
    end do
    call iqsortcols(.true., 1, n_vtx_L, 1, 2, vtx_L)
    i       =  1
    sR      =  n_vtx_L
    n_vtx_L =  0
    do while(i .le. sR)
        k   =  i
        do j=i+1,sR
            if(vtx_L(1,j) .ne. vtx_L(1,i)) then
                exit
            else
                k   =  j
            end if
        end do

        call iqsortcols(.true., i, k, 2, 2, vtx_L)
        n_vtx_L             =  n_vtx_L+1
        vtx_L(1:2,n_vtx_L)  =  vtx_L(1:2,i)
        do j=i+1,k
            if(vtx_L(2,j) .eq. vtx_L(2,j-1))    cycle
            n_vtx_L             =  n_vtx_L+1
            vtx_L(1:2,n_vtx_L)  =  vtx_L(1:2,j)
        end do

        i   =  k+1
    end do
!   number the vertex used to define the wall elements.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   set the coordinates of the wall elements.
    allocate(xyz_L(3,n_vtx_L), stat=err_mem)
    do i=1,n_vtx_L
        xyz_L(1:n_dim,i)=  mesh(0)%xyz(1:n_dim, vtx_L(1,i))
        if(is_2d_cal)   xyz_L(3,i)=  0.0d0

        iper=  vtx_L(2,i)
        if(iper .eq. 0) then
            cycle
        elseif(iper .gt. 0) then
            if(nint(sec_per_info(2,iper), dpI) .eq. 1) then
                xyz_L(1:3,i)=  xyz_L(1:3,i)+sec_per_info(3:5,iper)
            else
                d(1:3)  =  sec_per_info(3:5,iper)
                call norm_vec(3, d, rtmp)
                rtmp    =  2.0d0*pi/real(nint(2.0d0*pi/rtmp, dpI), dpR)
                call per_mat_rot(d, rtmp, rm)
                call rotate_vector(1, rm, xyz_L(1,i))
            end if
        else
            if(nint(sec_per_info(2,-iper), dpI) .eq. 1) then
                xyz_L(1:3,i)=  xyz_L(1:3,i)-sec_per_info(3:5,-iper)
            else
                d(1:3)  =  sec_per_info(3:5,-iper)
                call norm_vec(3, d, rtmp)
                rtmp    =  2.0d0*pi/real(nint(2.0d0*pi/rtmp, dpI), dpR)
                call per_mat_rot(d,-rtmp, rm)
                call rotate_vector(1, rm, xyz_L(1,i))
            end if
        end if
    end do
!   set the coordinates of the wall elements.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   allgather the coordinates of the wall elements.
    call mpi_allgather(n_vtx_L, 1, mpi_dpI, prc_info(0,1), 1, mpi_dpI, &
        &  mpi_comm_world, mpi_err)
    pivot   =  0
    do i=0,myid-1
        pivot   =  pivot+prc_info(i,1)
    end do

    n_vtx_wall  =  prc_info(0,1)
    do i=1,nprc-1
        n_vtx_wall  =  n_vtx_wall+prc_info(i,1)
    end do
    if(n_vtx_wall .le. 0)   return
    allocate(xyz_wall(3,n_vtx_wall), stat=err_mem)
    i   =  3*n_vtx_L
    call r_allgather(mpi_comm_world, nprc, i, xyz_L, xyz_wall)
    if(allocated(xyz_L))    deallocate(xyz_L)
!   allgather the coordinates of the wall elements.
!   ----------------------------------------------------------------------------

    n_ele_L =  0
    k       =  0
    do im=1,mesh(0)%n_mortar_b
        sR  =  mesh(0)%mortar_LR(3,im)
        eR  =  mesh(0)%mortar_LR(4,im)
        if(sec(sR)%bct .ne. BCWallViscous)  cycle
        isec_g  =  sec(sR)%ID_sec_g
        npe     =  sec(sR)%npe
        v(1:npe)=  sec(sR)%n2e(1:npe,eR)

!       the original.
        n_ele_L =  n_ele_L+1
        do i=1,npe
            call ib_search(1, n_vtx_L, 1, 2, vtx_L, v(i), ltmp, LR(1,i), LR(2,i))
            if(.not. ltmp)  stop 'Error: mesh_get_dnw fails.'
            call ib_search(LR(1,i), LR(2,i), 2, 2, vtx_L, 0, ltmp, L, R)
            if((.not. ltmp) .or. (L .ne. R))    stop 'Error: mesh_get_dnw fails.'
            k       =  k+1
            jA_L(k) =  L
        end do
        iA_L(  n_ele_L) =  npe
        ID_L(1,n_ele_L) =  myid
        ID_L(2,n_ele_L) =  im

        do iper=1,n_sec_per_info
            if(nint(sec_per_info(1,iper), dpI) .ne. isec_g) cycle

!           rotate in the positive direction.
            n_ele_L =  n_ele_L+1
            do i=1,npe
                call ib_search(LR(1,i), LR(2,i), 2, 2, vtx_L, iper, ltmp, L, R)
                if((.not. ltmp) .or. (L .ne. R))    stop 'Error: mesh_get_dnw fails.'
                k       =  k+1
                jA_L(k) =  L
            end do
            iA_L(  n_ele_L) =  npe
            ID_L(1,n_ele_L) =  myid
            ID_L(2,n_ele_L) =  im

!           rotate in the negative direction.
            n_ele_L =  n_ele_L+1
            do i=1,npe
                call ib_search(LR(1,i), LR(2,i), 2, 2, vtx_L, -iper, ltmp, L, R)
                if((.not. ltmp) .or. (L .ne. R))    stop 'Error: mesh_get_dnw fails.'
                k       =  k+1
                jA_L(k) =  L
            end do
            iA_L(  n_ele_L) =  npe
            ID_L(1,n_ele_L) =  myid
            ID_L(2,n_ele_L) =  im
        end do
    end do
    jA_L=  jA_L+pivot

    call mpi_allreduce((/n_ele_L, size(jA_L)/), v, 2, mpi_dpI, mpi_sum, &
        &  mpi_comm_world, mpi_err)
    n_ele_wall  =  v(1)
    allocate(iA_wall(v(1)+1), stat=err_mem)
    allocate(jA_wall(v(2)  ), stat=err_mem)
    i   =  n_ele_L
    call i_allgather(mpi_comm_world, nprc, i, iA_L, iA_wall(2))
    iA_wall(1)  =  1
    do i=1,n_ele_wall
        j           =  iA_wall(i)+iA_wall(i+1)
        iA_wall(i+1)=  j
    end do

    allocate(ID_wall(2,n_ele_wall), stat=err_mem)
    i   =  2*n_ele_L
    call i_allgather(mpi_comm_world, nprc, i, ID_L, ID_wall)
    if(allocated(ID_L)) deallocate(ID_L)

    i   =  size(jA_L)
    call i_allgather(mpi_comm_world, nprc, i, jA_L, jA_wall)
    if(allocated(iA_L)) deallocate(iA_L)
    if(allocated(jA_L)) deallocate(jA_L)

!   ----------------------------------------------------------------------------
!   setup rtree.
    call rtree_setup(3, n_vtx_wall, xyz_wall, tree)
!   setup rtree.
!   ----------------------------------------------------------------------------

!   call output_wall

    return
    end subroutine mesh_get_wall
