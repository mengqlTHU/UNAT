!-------------------------------------------------------------------------------
!   module GMG.
!-------------------------------------------------------------------------------
    module var_mg
        use var_kind_def
        implicit none

        logical(dpL):: is_amg               =  .true.
        logical(dpL):: is_mesh_agglomerated =  .false.
        integer(dpI):: mlev                 =  1
        type type_gmg
!           number of element and number of boundary element.
            integer(dpI):: n_ele    =  0
            integer(dpI):: n_ele_b  =  0
!           number of element on the coarse level.
            integer(dpI):: n_ele_c  =  0
            integer(dpI),allocatable:: ID(:)
        end type type_gmg
        type(type_gmg),allocatable:: mg(:)
    end module var_mg
!-------------------------------------------------------------------------------
!   sort the element in the wave-front pattern.
!-------------------------------------------------------------------------------
    subroutine sort_graph_wavefront(nele,iA,jA,n_wf,ID_wf)
    use var_kind_def
    implicit none
    integer(dpI),intent(in):: nele,iA(*),jA(*)
    integer(dpI):: n_wf,ID_wf(*),seed,i,j,k,iele,n_current,n_next,current(10000), &
                &  next(10000),itmp

    ID_wf(1:nele)   =  0
    n_wf            =  0
    do while(.true.)
!       first, we find the seed element.
        seed=  0
        i   =  huge(1)
        do iele=1,nele
            if(ID_wf(iele) .gt. 0)  cycle
            itmp=  iA(iele+1)-iA(iele)
            if(itmp .lt. i) then
                seed=  iele
                i   =  itmp
            end if
        end do
        if(seed .le. 0) exit

        n_current   =  1
        current(1)  =  seed
        do while(.true.)
            n_wf    =  n_wf+1
            n_next  =  0
            do i=1,n_current
                iele        =  current(i)
                ID_wf(iele) =  n_wf
                do j=iA(iele),iA(iele+1)-1
                    k           =  jA(j)
                    if(ID_wf(k) .gt. 0) cycle
                    n_next      =  n_next+1
                    next(n_next)=  k
                end do
            end do
            if(n_next .le. 0)   exit

            call simplify_series(n_next, 1, 1 ,next)
            n_current           =  n_next
            current(1:n_next)   =  next(1:n_next)
        end do
    end do

    return
    end subroutine sort_graph_wavefront
!-------------------------------------------------------------------------------
!   output 2D solution, tecplot.
!-------------------------------------------------------------------------------
    subroutine gmg_wr_ID_tec(ID)
    use var_kind_def
    use var_cgns
    use var_global, only: is_2d_cal,n_dim
    use var_mesh, only: sec,mesh
    use var_parallel
    implicit none
    integer(dpI),intent(in):: ID(*)
    character(len=80):: tec,str
    integer(dpI):: isec,iele,nvtx,i,v(8),ele1,ele0
    integer(dpI),allocatable:: ibuf(:)
    real   (dpR),allocatable:: rbuf(:,:)

    write(str,*),myid
    str =  './data/mesh_local_ID_'//trim(adjustl(str))//'.dat'
    open(unit=10,file=trim(adjustl(str)))
    if(is_2d_cal) then
        write(unit=10,fmt='(A)'),'variables="x","y","ID"'
    else
        write(unit=10,fmt='(A)'),'variables="x","y","z","ID"'
    end if

    allocate(ibuf(mesh(0)%n_vtx   ))
    allocate(rbuf(mesh(0)%n_vtx, 3))
    ele1=  1
    do isec=mesh(0)%sec_1,mesh(0)%sec_0
        if(.not. sec(isec)%is_int)  cycle
        ele0=  ele1+sec(isec)%n_ele-1

        ibuf=  0
        do iele=1,sec(isec)%n_ele
        do i=1,sec(isec)%npe
            ibuf(sec(isec)%n2e(i,iele)) =  1
        end do
        end do
        nvtx=  0
        do i=1,mesh(0)%n_vtx
            if(ibuf(i) .gt. 0) then
                nvtx    =  nvtx+1
                ibuf(i) =  nvtx
                rbuf(nvtx,1:n_dim)  =  mesh(0)%xyz(1:n_dim,i)
            end if
        end do

        tec =  'zone N='
        write(str,*),nvtx
        tec =  trim(adjustl(tec))//trim(adjustl(str))//',E='
        write(str,*),sec(isec)%n_ele
        if(is_2d_cal) then
        tec =  trim(adjustl(tec))//trim(adjustl(str))//',varlocation=([3]=cellcentered)'
        else
        tec =  trim(adjustl(tec))//trim(adjustl(str))//',varlocation=([4]=cellcentered)'
        end if
        if(sec(isec)%ele_type .eq. TRI_3) then
            tec =  trim(adjustl(tec))//',zonetype=FETRIANGLE'
        elseif(sec(isec)%ele_type .eq. QUAD_4) then
            tec =  trim(adjustl(tec))//',zonetype=FEQUADRILATERAL'
        elseif(sec(isec)%ele_type .eq. HEXA_8) then
            tec =  trim(adjustl(tec))//',zonetype=FEBRICK'
        elseif(sec(isec)%ele_type .eq. PENTA_6) then
            tec =  trim(adjustl(tec))//',zonetype=FEBRICK'
        else
            stop 'Error: element type not supported.'
        end if
        write(unit=10,fmt='(A)'),trim(adjustl(tec))
        write(unit=10,fmt=*),rbuf(1:nvtx,1)
        write(unit=10,fmt=*),rbuf(1:nvtx,2)
        if(.not. is_2d_cal) write(unit=10,fmt=*),rbuf(1:nvtx,3)
        write(unit=10,fmt=*),ID(ele1:ele0)
        do iele=1,sec(isec)%n_ele
            do i=1,sec(isec)%npe
                v(i)=  ibuf(sec(isec)%n2e(i,iele))
            end do

            if(sec(isec)%ele_type .eq. TRI_3) then
                write(unit=10,fmt='(3I8)'),v(1:3)
            elseif(sec(isec)%ele_type .eq. QUAD_4) then
                write(unit=10,fmt='(4I8)'),v(1:4)
            elseif(sec(isec)%ele_type .eq. HEXA_8) then
                write(unit=10,fmt='(8I8)'),v(1:8)
            elseif(sec(isec)%ele_type .eq. PENTA_6) then
                write(unit=10,fmt='(8I8)'),v(1:3),v(3),v(4:6),v(6)
            else
                stop 'Error: element type not supported.'
            end if
        end do
        ele1=  ele0+1
    end do
    close(10)

    if(allocated(ibuf)) deallocate(ibuf)
    if(allocated(rbuf)) deallocate(rbuf)

    return
    end subroutine gmg_wr_ID_tec
!-------------------------------------------------------------------------------
!   define the coarse level mesh.
!-------------------------------------------------------------------------------
    subroutine mesh_get_coarse_level(L0)
    use var_kind_def
    use var_global, only: err_mem
    use var_mesh
    use var_mg
    use var_parallel
    implicit none
    integer(dpI),intent(in):: L0
    logical(dpL):: ltmp
    integer(dpI):: L1,isr,isr0,isec,iele,nele,i,j,L,R,s,e,ip_remote,sL,eL,sR,eR, &
                &  im,v(4)
    integer(dpI),allocatable:: ibuf(:),iA(:),s_coarse(:)

    if(L0 .ge. mlev-1)  return
    L1  =  L0+1

    allocate(s_coarse(mesh(L0)%sec_1:mesh(L0)%sec_0), stat=err_mem)

!   ----------------------------------------------------------------------------
!   set the internal section of the coarse mesh.
    s               =  mesh(L0)%sec_0+1
    mesh(L1)%sec_1  =  s
    mesh(L1)%sec_0  =  s
    sec(s)%n_ele    =  mg(L0)%n_ele_c
    sec(s)%n_ele_b  =  mg(L1)%n_ele_b
    sec(s)%is_int   =  .true.
    sec(s)%is_ghost =  .false.
    sec(s)%is_bnd   =  .false.
    allocate(sec(s)%vol(sec(s)%n_ele), stat=err_mem)
    sec(s)%vol      =  0.0d0

    i   =  0
    do isec=mesh(L0)%sec_1,mesh(L0)%sec_0
        if(.not. sec(isec)%is_int)  cycle
        allocate(sec(isec)%ID_coarse_mesh(2,sec(isec)%n_ele), stat=err_mem)
        s_coarse(isec)  =  mesh(L1)%sec_1
        do iele=1,sec(isec)%n_ele
            i   =  i+1
            s   =  mesh(L1)%sec_1
            e   =  mg(L0)%ID(i)
            sec(isec)%ID_coarse_mesh(1,iele)=  s
            sec(isec)%ID_coarse_mesh(2,iele)=  e
            sec(s)%vol(e)   =  sec(s)%vol(e)+sec(isec)%vol(iele)
        end do
    end do

    allocate(iA(mesh(L1)%sec_1:mesh(L1)%sec_0+1), stat=err_mem)
    iA(mesh(L1)%sec_1)  =  1
    do isec=mesh(L1)%sec_1,mesh(L1)%sec_0
        iA(isec+1)=  iA(isec)+sec(isec)%n_ele
    end do
!   set the internal section of the coarse mesh.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   set p2p of the coarse mesh.
    mesh(L1)%p2p_1  =  mesh(L0)%p2p_0+1
    mesh(L1)%p2p_0  =  mesh(L1)%p2p_1+mesh(L0)%p2p_0-mesh(L0)%p2p_1
    nele            =  0
    do isr=mesh(L1)%p2p_1,mesh(L1)%p2p_0
        isr0                =  isr-(mesh(L0)%p2p_0-mesh(L0)%p2p_1+1)
        p2p(isr)%ip_remote  =  p2p(isr0)%ip_remote
        if(p2p(isr0)%n_ele_send .gt. 0) then
            nele=  max(nele, p2p(isr0)%n_ele_send)
            if(allocated(p2p(isr0)%isend))  deallocate(p2p(isr0)%isend)
            allocate(p2p(isr0)%isend(p2p(isr0)%n_ele_send), stat=err_mem)
        end if
        if(p2p(isr0)%n_ele_recv .gt. 0) then
            if(allocated(p2p(isr0)%irecv))  deallocate(p2p(isr0)%irecv)
            allocate(p2p(isr0)%irecv(p2p(isr0)%n_ele_recv), stat=err_mem)
        end if
    end do
    allocate(ibuf(max(3*nele, 5*mesh(L0)%n_mortar)), stat=err_mem)
!   set p2p of the coarse mesh.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   record the element to be sent out on the coarse mesh.
    do isr0=mesh(L0)%p2p_1,mesh(L0)%p2p_0
        if(p2p(isr0)%n_ele_send .le. 0) cycle
        isr =  isr0+(mesh(L0)%p2p_0-mesh(L0)%p2p_1+1)

        do i=1,p2p(isr0)%n_ele_send
            isec=  p2p(isr0)%id_ele_send(1,i)
            iele=  p2p(isr0)%id_ele_send(2,i)

            s   =  sec(isec)%ID_coarse_mesh(1,iele)
            e   =  sec(isec)%ID_coarse_mesh(2,iele)
            ibuf(3*i-2) =  s
            ibuf(3*i-1) =  e
            ibuf(3*i  ) =  e+iA(s)-1
        end do
        p2p(isr)%n_ele_send =  p2p(isr0)%n_ele_send
        call simplify_series(p2p(isr)%n_ele_send, 3, 3, ibuf)
        allocate(p2p(isr)%id_ele_send(2,p2p(isr)%n_ele_send), stat=err_mem)
        do i=1,p2p(isr)%n_ele_send
            p2p(isr)%id_ele_send(1,i)   =  ibuf(3*i-2)
            p2p(isr)%id_ele_send(2,i)   =  ibuf(3*i-1)
        end do

        do i=1,p2p(isr0)%n_ele_send
            isec=  p2p(isr0)%id_ele_send(1,i)
            iele=  p2p(isr0)%id_ele_send(2,i)

            s   =  sec(isec)%ID_coarse_mesh(1,iele)
            e   =  sec(isec)%ID_coarse_mesh(2,iele)
            j   =  e+iA(s)-1
            call ib_search(1, p2p(isr)%n_ele_send, 3, 3, ibuf, j, ltmp, L, R)
            if((.not. ltmp) .or. (L .ne. R))    stop 'Error: mesh_coarsen fails.'
            p2p(isr0)%isend(i)  =  L
        end do
    end do
!   record the element to be sent out on the coarse mesh.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   send the coarse element index to the other part.
    mpi_nreq=  0
    do isr=mesh(L0)%p2p_1,mesh(L0)%p2p_0
        ip_remote   =  p2p(isr)%ip_remote
        if(myid .eq. ip_remote) then
            if(p2p(isr)%n_ele_send .ne. p2p(isr)%n_ele_recv)    stop 'Error: S!=R.'
            if(p2p(isr)%n_ele_send .gt. 0)  call ICOPY(p2p(isr)%n_ele_send, &
                &  p2p(isr)%isend, 1, p2p(isr)%irecv, 1)
            cycle
        end if

        if(p2p(isr)%n_ele_send .gt. 0) then
            mpi_nreq=  mpi_nreq+1
            if(mpi_nreq .gt. max_p2p)   stop 'Error: you are working with a bad LB.'
            call mpi_isend(p2p(isr)%isend, p2p(isr)%n_ele_send, mpi_dpI, ip_remote, &
                &  myid, mpi_comm_world, mpi_req(mpi_nreq), mpi_err)
        end if

        if(p2p(isr)%n_ele_recv .gt. 0) then
            mpi_nreq=  mpi_nreq+1
            if(mpi_nreq .gt. max_p2p)   stop 'Error: you are working with a bad LB.'
            call mpi_irecv(p2p(isr)%irecv, p2p(isr)%n_ele_recv, mpi_dpI, ip_remote, &
                &  ip_remote, mpi_comm_world, mpi_req(mpi_nreq), mpi_err)
        end if
    end do
    if(mpi_nreq .gt. 0) call mpi_waitall(mpi_nreq, mpi_req, mpi_sta, mpi_err)
!   send the coarse element index to the other part.
!   ----------------------------------------------------------------------------

    do isr=mesh(L1)%p2p_1,mesh(L1)%p2p_0
        isr0=  isr-(mesh(L0)%p2p_0-mesh(L0)%p2p_1+1)
        if(p2p(isr0)%n_ele_recv .le. 0) cycle

        p2p(isr)%n_ele_recv =  0
        do i=1,p2p(isr0)%n_ele_recv
            p2p(isr)%n_ele_recv =  max(p2p(isr)%n_ele_recv, p2p(isr0)%irecv(i))
        end do
    end do

!   ----------------------------------------------------------------------------
!   record the ghost element of the coarse level.
    nele=  0
    do isec=mesh(L0)%sec_1,mesh(L0)%sec_0
        if(.not. sec(isec)%is_ghost)    nele=  max(nele, sec(isec)%n_ele)
    end do
    if(allocated(ibuf)) then
        if(size(ibuf) .lt. 6*nele) then
            deallocate(ibuf)
            allocate(ibuf(5*nele), stat=err_mem)
        end if
    else
        allocate(ibuf(5*nele), stat=err_mem)
    end if
    do isec=mesh(L0)%sec_1,mesh(L0)%sec_0
        if(.not. sec(isec)%is_ghost)    cycle
        do i=1,sec(isec)%n_ele
            isr0=  sec(isec)%id_recv(1,i)
            iele=  sec(isec)%id_recv(2,i)

            ibuf(5*i-4) =  p2p(isr0)%irecv(iele)
            ibuf(5*i-3) =  isr0+(mesh(L0)%p2p_0-mesh(L0)%p2p_1+1)
            ibuf(5*i-2) =  sec(isec)%per_path(1,i)
            ibuf(5*i-1) =  sec(isec)%per_path(2,i)
            ibuf(5*i  ) =  sec(isec)%per_path(3,i)
        end do

        s               =  mesh(L1)%sec_0+1
        mesh(L1)%sec_0  =  s
        sec(s)%is_int   =  .false.
        sec(s)%is_ghost =  .true.
        sec(s)%is_bnd   =  .false.
        sec(s)%n_ele    =  sec(isec)%n_ele
        s_coarse(isec)  =  s
        call simplify_matrix(5, sec(s)%n_ele, 1, 5, ibuf)
        allocate(sec(s)%ID_recv (2,sec(s)%n_ele), stat=err_mem)
        allocate(sec(s)%per_path(3,sec(s)%n_ele), stat=err_mem)
        do i=1,sec(s)%n_ele
            sec(s)%ID_recv (1  ,i)  =  ibuf(2+5*(i-1))
            sec(s)%ID_recv (2  ,i)  =  ibuf(1+5*(i-1))
            sec(s)%per_path(1:3,i)  =  ibuf(3+5*(i-1):5+5*(i-1))
        end do

        allocate(sec(isec)%ID_coarse_mesh(2,sec(isec)%n_ele), stat=err_mem)
        do i=1,sec(isec)%n_ele
            isr0=  sec(isec)%id_recv(1,i)
            iele=  sec(isec)%id_recv(2,i)
            call ib_search(1, sec(s)%n_ele, 1, 5, ibuf, p2p(isr0)%irecv(iele), ltmp, L, R)
            if(.not. ltmp)  stop 'Error: mesh_coarsen fails.'
            ltmp=  .true.
            do j=L,R
                if((ibuf(5*j-3) .eq. (isr0+(mesh(L0)%p2p_0-mesh(L0)%p2p_1+1))) .and. &
                 & (ibuf(5*j-2) .eq. sec(isec)%per_path(1,i)) .and. &
                 & (ibuf(5*j-1) .eq. sec(isec)%per_path(2,i)) .and. &
                 & (ibuf(5*j  ) .eq. sec(isec)%per_path(3,i))) then
                    sec(isec)%ID_coarse_mesh(1,i)   =  s
                    sec(isec)%ID_coarse_mesh(2,i)   =  j
                    ltmp                            =  .false.
                end if
            end do
            if(ltmp)    stop 'Error: mesh_coarsen fails.'
        end do
    end do
!   record the ghost element of the coarse level.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   set boundary section of the coarse mesh.
    do isec=mesh(L0)%sec_1,mesh(L0)%sec_0
        if(.not. sec(isec)%is_bnd)  cycle

        s               =  mesh(L1)%sec_0+1
        mesh(L1)%sec_0  =  s
        sec(s)%is_int   =  .false.
        sec(s)%is_ghost =  .false.
        sec(s)%is_bnd   =  .true.
        sec(s)%n_ele    =  sec(isec)%n_ele
        sec(s)%bct      =  sec(isec)%bct
        sec(s)%ID_group =  sec(isec)%ID_group
        sec(s)%ID_bnd   =  sec(isec)%ID_bnd
        s_coarse(isec)  =  s
    end do
!   set boundary section of the coarse mesh.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   set mortar of the coarse mesh.
    mesh(L1)%n_mortar   =  0
    do im=1,mesh(L0)%n_mortar
        sL  =  mesh(L0)%mortar_LR(1,im)
        eL  =  mesh(L0)%mortar_LR(2,im)
        sR  =  mesh(L0)%mortar_LR(3,im)
        eR  =  mesh(L0)%mortar_LR(4,im)

        if(sec(sL)%is_int .or. sec(sL)%is_ghost) then
            v(1)=  sec(sL)%ID_coarse_mesh(1,eL)
            v(2)=  sec(sL)%ID_coarse_mesh(2,eL)
        elseif(sec(sL)%is_bnd) then
            v(1)=  s_coarse(sL)
            v(2)=  eL
        else
            stop 'Error: mesh_coarsen fails.'
        end if

        if(sec(sR)%is_int .or. sec(sR)%is_ghost) then
            v(3)=  sec(sR)%ID_coarse_mesh(1,eR)
            v(4)=  sec(sR)%ID_coarse_mesh(2,eR)
        elseif(sec(sR)%is_bnd) then
            v(3)=  s_coarse(sR)
            v(4)=  eR
        else
            stop 'Error: mesh_coarsen fails.'
        end if

        if(sec(v(1))%is_int .and. sec(v(3))%is_int) then
            ltmp= (v(1) .eq. v(3)) .and. (v(2) .eq. v(4))
        else
            ltmp=  .false.
        end if
        if(ltmp)    cycle
        mesh(L1)%n_mortar   =  mesh(L1)%n_mortar+1
        ibuf(5*mesh(L1)%n_mortar-4:5*mesh(L1)%n_mortar-1)   =  v(1:4)
        ibuf(5*mesh(L1)%n_mortar                        )   =  im
    end do
    allocate(mesh(L1)%mortar_LR  (4,mesh(L1)%n_mortar), stat=err_mem)
    allocate(mesh(L1)%mortar_n_vg(8,mesh(L1)%n_mortar), stat=err_mem)
    mesh(L1)%n_mortar_b =  0
    do im=1,mesh(L1)%n_mortar
        sL  =  ibuf(5*im-4)
        eL  =  ibuf(5*im-3)
        sR  =  ibuf(5*im-2)
        eR  =  ibuf(5*im-1)
        if(sec(sL)%is_bnd)  stop 'Error: left side of mortar is boundary.'
        if(sec(sR)%is_bnd)  mesh(L1)%n_mortar_b =  im
        mesh(L1)%mortar_LR  (1:4,im)=  ibuf(5*im-4:5*im-1)
        mesh(L1)%mortar_n_vg(1:8,im)=  mesh(L0)%mortar_n_vg(1:8,ibuf(5*im))
    end do
!   set mortar of the coarse mesh.
!   ----------------------------------------------------------------------------

    call set_p2p_element_order(L1)

    call mesh_check(L1)

    call mesh_set_face_neighbour(L1)

    do isr=mesh(L0)%p2p_1,mesh(L0)%p2p_0
        if(allocated(p2p(isr)%isend))   deallocate(p2p(isr)%isend)
        if(allocated(p2p(isr)%irecv))   deallocate(p2p(isr)%irecv)
    end do
    if(allocated(ibuf    )) deallocate(ibuf    )
    if(allocated(iA      )) deallocate(iA      )
    if(allocated(s_coarse)) deallocate(s_coarse)

    return
    end subroutine mesh_get_coarse_level
!-------------------------------------------------------------------------------
!   pairwise agglomeration MG.
!-------------------------------------------------------------------------------
    module var_agglomeration_pairwise
        use var_kind_def
        implicit none
        private
        public:: mesh_agglomeration_pairwise,mesh_agglomeration_pairwise_destroy
        public:: mg

        type type_mg
            integer(dpI):: nele     =  0
            integer(dpI):: nele_c   =  0
            integer(dpI),allocatable:: iA(:),jA(:)
            integer(dpI),allocatable:: perm(:)
            integer(dpI),allocatable:: ID_c(:)
            real   (dpR),allocatable:: A(:)
        end type type_mg
        type(type_mg),save:: mg(0:10)

        contains
!       ------------------------------------------------------------------------
!       get the LHS for diffusion equation.
!       ------------------------------------------------------------------------
        subroutine get_LHS_diffusion
        use var_mesh
        implicit none
        integer(dpI):: v(100),isec,iele,j,k,m,im,s,e,err_mem
        real   (dpR):: c(2,20),rtmp

        mg(0)%nele  =  0
        v           =  0
        j           =  0
        do isec=mesh(0)%sec_1,mesh(0)%sec_0
            if(.not. sec(isec)%is_int)  cycle
            v(isec)     =  mg(0)%nele+1
            mg(0)%nele  =  mg(0)%nele+sec(isec)%n_ele
            do iele=1,sec(isec)%n_ele
                j   =  j+sec(isec)%iA_face_neighbour(iele+1) &
                    &   -sec(isec)%iA_face_neighbour(iele)+1
            end do
        end do
        allocate(mg(0)%iA(mg(0)%nele+1), stat=err_mem)
        allocate(mg(0)%jA(j           ), stat=err_mem)
        allocate(mg(0)%A (j           ), stat=err_mem)
        mg(0)%iA=  0
        mg(0)%jA=  0
        mg(0)%A =  0.0d0

        mg(0)%iA(1) =  1
        do isec=mesh(0)%sec_1,mesh(0)%sec_0
        if(.not. sec(isec)%is_int)  cycle
        do iele=1,sec(isec)%n_ele
            m       =  1
            c(1,1)  =  real(iele+v(isec)-1, dpR)
            c(2,1)  =  0.0d0
            do k=sec(isec)%iA_face_neighbour(iele),sec(isec)%iA_face_neighbour(iele+1)-1
                im  =  sec(isec)%jA_face_neighbour(k)
                if(im .gt. 0) then
                    s   =  mesh(0)%mortar_LR(3, im)
                    e   =  mesh(0)%mortar_LR(4, im)
                else
                    s   =  mesh(0)%mortar_LR(1,-im)
                    e   =  mesh(0)%mortar_LR(2,-im)
                end if
                if(.not. sec(s)%is_int) cycle
                m   =  m+1
                rtmp=  mesh(0)%mortar_n_vg(4,abs(im))*mesh(0)%mortar_n_vg(5,abs(im))
                c(2,1)  =  c(2,1)+rtmp
                c(1,m)  =  real(e+v(s)-1, dpR)
                c(2,m)  = -rtmp
            end do

            call dqsortcols(.true., 1, m, 1, 2, c)
            do j=1,m
                mg(0)%jA(j+mg(0)%iA(iele+v(isec)-1)-1)  =  nint(c(1,j), dpI)
                mg(0)%A     (j+mg(0)%iA(iele+v(isec)-1)-1)  =  c(2,j)
            end do
            mg(0)%iA(iele+v(isec))  =  mg(0)%iA(iele+v(isec)-1)+m
        end do
        end do

        return
        end subroutine get_LHS_diffusion
!       ------------------------------------------------------------------------
!       pairwise aggregation.
!       ------------------------------------------------------------------------
        subroutine pairwise_aggregation(N,iA,jA,A,n_c,ID_c)
        use var_kind_def
        implicit none
        integer(dpI),intent(in):: N,iA(*),jA(*)
        real   (dpR),intent(in):: A(*)
        real   (dpR),parameter:: threshold  =  1.0d1
        logical(dpL):: ltmp
        integer(dpI):: n_c,ID_c(*),err_mem,i,j,k,ii,jj,ij,ji,next,seed
        real   (dpR):: v1,v2,mu,rtmp
        integer(dpI),allocatable:: perm(:),ordered_list(:)
        real   (dpR),allocatable:: s(:)

        if(N .le. 0)    stop 'Error: wrong input, pairwise aggregation.'

        allocate(perm        (N), stat=err_mem)
        allocate(ordered_list(N), stat=err_mem)
        call matrix_ordering(N, iA, jA, perm)
        do i=1,N
            ordered_list(perm(i))   =  i
        end do

        allocate(s(n), stat=err_mem)
        s   =  0.0d0
        do i=1,N
        do k=iA(i),iA(i+1)-1
            j   =  jA(k)
            if(i .eq. j)    cycle
            s(i)=  s(i)-0.5d0*A(k)
            s(j)=  s(j)-0.5d0*A(k)
        end do
        end do

        ID_c(1:N)   =  0
        n_c         =  0
        seed        =  1
        do while(.true.)
!           find the seed element.
            i   = -1
            do j=seed,N
                if(ID_c(ordered_list(j)) .le. 0) then
                    i   =  ordered_list(j)
                    seed=  j
                    exit
                end if
            end do
            if(i .le. 0)    exit

            call existelm(iA, jA, i, i, ltmp, ii)
            if(.not. ltmp)  stop 'Error: fails to get A(i,i), pairwise aggregation.'

            mu  =  huge(1.0d0)
            next=  0
            do k=iA(i),iA(i+1)-1
                j   =  jA(k)
                if((i .eq. j) .or. (ID_c(j) .gt. 0))    cycle

                call existelm(iA, jA, j, j, ltmp, jj)
                if(.not. ltmp)  stop 'Error: fails to get A(j,j), pairwise aggregation.'
!               if(A(ii)-s(i)+A(jj)-s(j) .lt. 0.0d0)    cycle

                call existelm(iA, jA, i, j, ltmp, ij)
                if(.not. ltmp)  stop 'Error: fails to get A(i,j), pairwise aggregation.'
                call existelm(iA, jA, j, i, ltmp, ji)
                if(.not. ltmp)  stop 'Error: fails to get A(j,i), pairwise aggregation.'
                v1  =  2.0d0/(1.0d0/A(ii)+1.0d0/A(jj))
                v2  = -0.5d0*(A(ij)+A(ji))+1.0d0/(1.0d0/(A(ii)-s(i))+1.0d0/(A(jj)-s(j)))
                rtmp=  v1/v2
                if(rtmp .le. 0.0d0) cycle
                if(rtmp .lt. mu) then
                    mu  =  rtmp
                    next=  j
                end if
            end do

            n_c =  n_c+1
            if(mu .le. threshold) then
                ID_c(i   )  =  n_c
                ID_c(next)  =  n_c
            else
                ID_c(i   )  =  n_c
            end if
        end do

        if(allocated(s))    deallocate(s)

        return
        end subroutine pairwise_aggregation
!       ------------------------------------------------------------------------
!       construct iA and jA from the finer mesh.
!       ------------------------------------------------------------------------
        subroutine get_coarse_topology(L1)
        use var_kind_def
        use var_global, only: err_mem
        implicit none
        integer(dpI),intent(in):: L1
        logical(dpL):: ltmp
        integer(dpI):: L0,ne2e,iele,i,j,k,L,R,nb,v(1000),ii,jj
        integer(dpI),allocatable:: e2e(:,:)

        L0          =  L1-1
        mg(L1)%nele =  mg(L0)%nele_c
        allocate(mg(L1)%iA(mg(L1)%nele+1), stat=err_mem)
        allocate(e2e(2,size(mg(L0)%jA)), stat=err_mem)
        ne2e=  0
        do iele=1,mg(L0)%nele
            L   =  mg(L0)%ID_C(iele)
            do j=mg(L0)%iA(iele),mg(L0)%iA(iele+1)-1
                R               =  mg(L0)%ID_C(mg(L0)%jA(j))
                ne2e            =  ne2e+1
                e2e(1:2,ne2e)   = (/L, R/)
            end do
        end do
        call iqsortcols(.true., 1, ne2e, 1, 2, e2e)

!       setup iA.
        mg(L1)%iA   =  0
        i   =  1
        do while(i .le. ne2e)
            k   =  i
            do j=i+1,ne2e
                if(e2e(1,j) .eq. e2e(1,i)) then
                    k   =  j
                else
                    exit
                end if
            end do

            do j=i,k
                v(j-i+1)=  e2e(2,j)
            end do
            nb  =  k-i+1
            call simplify_series(nb, 1, 1, v)
            mg(L1)%iA(e2e(1,i)+1)   =  nb

            i   =  k+1
        end do

!       allocate memory for jA.
        j           =  0
        mg(L1)%iA(1)=  1
        do iele=1,mg(L1)%nele
            j   =  j+mg(L1)%iA(iele+1)
            mg(L1)%iA(iele+1)   =  j+1
        end do
        if(j .le. 0)    stop 'Error: fails to get the topology of the coarse level.'
        allocate(mg(L1)%jA(j), stat=err_mem)
        allocate(mg(L1)% A(j), stat=err_mem)
        mg(L1)%jA   =  0
        mg(L1)% A   =  0.0d0

        i   =  1
        do while(i .le. ne2e)
            k   =  i
            do j=i+1,ne2e
                if(e2e(1,j) .eq. e2e(1,i)) then
                    k   =  j
                else
                    exit
                end if
            end do

            do j=i,k
                v(j-i+1)=  e2e(2,j)
            end do
            nb  =  k-i+1
            call simplify_series(nb, 1, 1, v)

            iele=  e2e(1,i)
            do j=1,nb
                mg(L1)%jA(mg(L1)%iA(iele)+j-1)  =  v(j)
            end do

            i   =  k+1
        end do

        do i=1,mg(L0)%nele
            ii  =  mg(L0)%ID_c(i)
            do k=mg(L0)%iA(i),mg(L0)%iA(i+1)-1
                jj  =  mg(L0)%ID_c(mg(L0)%jA(k))
                ltmp=  .true.
                do L=mg(L1)%iA(ii),mg(L1)%iA(ii+1)-1
                    if(mg(L1)%jA(L) .eq. jj) then
                        ltmp        =  .false.
                        mg(L1)%A(L) =  mg(L1)%A(L)+mg(L0)%A(k)
                        exit
                    end if
                end do
                if(ltmp)    stop 'Error: fails to find coarse (i,j).'
            end do
        end do

        if(allocated(e2e))  deallocate(e2e)

        return
        end subroutine get_coarse_topology
!       ------------------------------------------------------------------------
!       mesh agglomeration.
!       ------------------------------------------------------------------------
        subroutine mesh_agglomeration_pairwise(mlev)
        use var_global, only: n_dim
        implicit none
        integer(dpI),intent(in):: mlev
        integer(dpI):: lev,i
        integer(dpI),allocatable:: ID(:)

        call get_LHS_diffusion
        allocate(ID(mg(0)%nele))
        forall(lev=1:mg(0)%nele)    ID(lev) =  lev

        do lev=0,n_dim*(mlev-1)-1
            allocate(mg(lev)%ID_c(mg(lev)%nele))
            call pairwise_aggregation(mg(lev)%nele,mg(lev)%iA,mg(lev)%jA,mg(lev)%A, &
                &  mg(lev)%nele_c, mg(lev)%ID_c)

            call get_coarse_topology(lev+1)

            do i=1,mg(0)%nele
                ID(i)   =  mg(lev)%ID_c(ID(i))
            end do
        end do
!       do lev=0,n_dim*(mlev-1)
!           print*,lev,'|',mg(lev)%nele
!       end do

!       call gmg_wr_ID_tec(ID)
        if(allocated(ID))   deallocate(ID)

        return
        end subroutine mesh_agglomeration_pairwise
!       ------------------------------------------------------------------------
!       mesh agglomeration destroy.
!       ------------------------------------------------------------------------
        subroutine mesh_agglomeration_pairwise_destroy
        implicit none
        integer(dpI):: lev

        do lev=0,10
            if(allocated(mg(lev)%iA  )) deallocate(mg(lev)%iA)
            if(allocated(mg(lev)%jA  )) deallocate(mg(lev)%jA)
            if(allocated(mg(lev)%perm)) deallocate(mg(lev)%perm)
            if(allocated(mg(lev)%ID_c)) deallocate(mg(lev)%ID_c)
            if(allocated(mg(lev)%A   )) deallocate(mg(lev)%A   )
            mg(lev)%nele= -1
        end do

        return
        end subroutine mesh_agglomeration_pairwise_destroy
    end module
!-------------------------------------------------------------------------------
!   module GMG.
!-------------------------------------------------------------------------------
    module var_agglomeration_circle
        use var_kind_def
        implicit none
        private
        public:: mg,mesh_agglomeration_circle,mesh_agglomeration_circle_destroy

        type type_gmg
!           number of element and number of boundary element.
            integer(dpI):: n_ele    =  0
            integer(dpI):: n_ele_b  =  0
!           number of element on the coarse level.
            integer(dpI):: n_ele_c  =  0
            integer(dpI),allocatable:: ID(:)
            integer(dpI),allocatable:: iA_e2e(:),jA_e2e(:)
            real   (dpR),allocatable:: A(:),vol(:)

            integer(dpI):: n_wf =  0
            integer(dpI),allocatable:: ID_wf(:)
        end type type_gmg
        type(type_gmg),save:: mg(0:9)

        contains
!       ------------------------------------------------------------------------
!       construct coarse mesh.
!       ------------------------------------------------------------------------
        subroutine mesh_agglomeration_circle(mlev)
        use var_kind_def
        use var_global, only: err_mem
        use var_parallel
        implicit none
        integer(dpI),intent(in):: mlev
        integer(dpI):: L1,L0,i
        integer(dpI),allocatable:: ID(:)

!       get the e2e information of the finest level.
        call gmg_get_e2e

        if(mlev .gt. 1) then
            allocate(ID(mg(0)%n_ele), stat=err_mem)
            forall(L1=1:mg(0)%n_ele)    ID(L1)  =  L1
        end if

        do L1=1,mlev-1
            L0  =  L1-1
            allocate(mg(L0)%ID_wf(mg(L0)%n_ele), stat=err_mem)
            allocate(mg(L0)%ID   (mg(L0)%n_ele), stat=err_mem)

            call sort_graph_wavefront(mg(L0)%n_ele, mg(L0)%iA_e2e, mg(L0)%jA_e2e, &
                &  mg(L0)%n_wf, mg(L0)%ID_wf)
            call gmg_step_ag(mg(L0)%n_ele, mg(L0)%iA_e2e, mg(L0)%jA_e2e, mg(L0)%n_wf, &
                &  mg(L0)%ID_wf, mg(L0)%n_ele_c, mg(L0)%ID)
            call gmg_push_boundary_element(L0)
            if(L1 .ne. mlev-1)  call gmg_get_coarse_topology(L1)

            do i=1,mg(0)%n_ele
                ID(i)   =  mg(L0)%ID   (ID(i))
!               ID(i)   =  mg(L0)%ID_wf(ID(i))
            end do
        end do

!       call gmg_wr_ID_tec(ID)
        if(allocated(ID))   deallocate(ID)

        return
        end subroutine mesh_agglomeration_circle
!       ------------------------------------------------------------------------
!       get e2e of the original mesh.
!       ------------------------------------------------------------------------
        subroutine gmg_get_e2e
        use var_kind_def
        use var_global, only: err_mem
        use var_mesh
        implicit none
        integer(dpI):: nele,isec,iele,nb,s,e,im,i,j,k,nfac,itmv(1000)
        real   (dpR):: A

        itmv(1) =  1
        do isec=2,mesh(0)%sec_0
            if(.not. sec(isec)%is_int)  exit
            itmv(isec)  =  itmv(isec-1)+sec(isec-1)%n_ele
        end do

!       ------------------------------------------------------------------------
!       record volume.
        nele=  0
        do isec=1,mesh(0)%sec_0
            if(.not. sec(isec)%is_int)  exit
            nele=  nele+sec(isec)%n_ele
        end do
        allocate(mg(0)%vol(nele), stat=err_mem)
        nele=  0
        do isec=1,mesh(0)%sec_0
            if(.not. sec(isec)%is_int)  exit
            do iele=1,sec(isec)%n_ele
                nele            =  nele+1
                mg(0)%vol(nele) =  sec(isec)%vol(iele)
            end do
        end do
!       record volume.
!       ------------------------------------------------------------------------

        allocate(mg(0)%iA_e2e(nele+1), stat=err_mem)
        mg(0)%n_ele     =  nele
        mg(0)%iA_e2e    =  0
        mg(0)%iA_e2e(1) =  1
        do isec=1,mesh(0)%sec_0
            if(.not. sec(isec)%is_int)  cycle
            do iele=1,sec(isec)%n_ele
                nb  =  0
                do j=sec(isec)%iA_face_neighbour(iele),sec(isec)%iA_face_neighbour(iele+1)-1
                    im  =  sec(isec)%jA_face_neighbour(j)
                    if(im .gt. 0) then
                        s   =  mesh(0)%mortar_LR(3, im)
                        e   =  mesh(0)%mortar_LR(4, im)
                    else
                        s   =  mesh(0)%mortar_LR(1,-im)
                        e   =  mesh(0)%mortar_LR(2,-im)
                    end if
                    if(.not. sec(s)%is_int) cycle
                    nb  =  nb+1
                end do
                j   =  iele+itmv(isec)-1
                mg(0)%iA_e2e(j+1)   =  mg(0)%iA_e2e(j)+nb
            end do
        end do
        nfac=  mg(0)%iA_e2e(nele+1)-1
        if(nfac .le. 0) stop 'Error: fails to get_e2e.'
        allocate(mg(0)%jA_e2e(nfac), stat=err_mem)
        allocate(mg(0)%A     (nfac), stat=err_mem)

        do isec=1,mesh(0)%sec_0
            if(.not. sec(isec)%is_int)  cycle
            do iele=1,sec(isec)%n_ele
                i   =  iele+itmv(isec)-1

                nb  =  0
                do j=sec(isec)%iA_face_neighbour(iele),sec(isec)%iA_face_neighbour(iele+1)-1
                    im  =  sec(isec)%jA_face_neighbour(j)
                    if(im .gt. 0) then
                        s   =  mesh(0)%mortar_LR  (3, im)
                        e   =  mesh(0)%mortar_LR  (4, im)
                        A   =  mesh(0)%mortar_n_vg(4, im)
                    else
                        s   =  mesh(0)%mortar_LR  (1,-im)
                        e   =  mesh(0)%mortar_LR  (2,-im)
                        A   =  mesh(0)%mortar_n_vg(4,-im)
                    end if
                    if(.not. sec(s)%is_int) cycle
                    nb  =  nb+1
                    k   =  mg(0)%iA_e2e(i)+nb-1
                    mg(0)%jA_e2e(k) =  e+itmv(s)-1
                    mg(0)%A     (k) =  A
                end do
            end do
        end do

        return
        end subroutine gmg_get_e2e
!       ------------------------------------------------------------------------
!       coarsen the mesh with the aggregation method.
!       ------------------------------------------------------------------------
        subroutine gmg_step_ag(nele,iA_e2e,jA_e2e,n_wf,ID_wf,nele_c,ID)
        use var_kind_def
        use var_global, only: is_2d_cal,err_mem
        use var_parallel
        implicit none
        integer(dpI),intent(in):: nele,iA_e2e(*),jA_e2e(*),n_wf,ID_wf(*)
        logical(dpL):: ltmp
        integer(dpI):: nele_c,ID(*),iele,i,j,k,n_cd,cd(2,1000),path(8,2000),n_c,c(8), &
                    &  max_ele_circle,max_distance_circle
        integer(dpI),allocatable:: iA(:),jA(:),ag(:)

        if(is_2d_cal) then
            max_ele_circle      =  4
            max_distance_circle =  2
        else
            max_ele_circle      =  8
            max_distance_circle =  3
        end if

        allocate(iA(n_wf+1), stat=err_mem)
        allocate(jA(nele  ), stat=err_mem)
        call itrans_spy(.true., 1, nele, cd, ID_wf, n_wf, iA, jA)

        ID(1:nele)  =  0
        nele_c      =  0

        do i=1,nele
            iele=  jA(i)
            if(ID(iele) .gt. 0) cycle

!           only aggregates the element whose neighbours have not been aggregated. 
            ltmp=  .false.
            do j=iA_e2e(iele),iA_e2e(iele+1)-1
                ltmp=  ID(jA_e2e(j)) .gt. 0
                if(ltmp)    exit
            end do
            if(ltmp)    cycle

            nele_c  =  nele_c+1

!           try to find a circle in the graph.
            call get_circle(max_ele_circle, max_distance_circle, iele, n_c, c)
            if(n_c .ge. 3)  ID(c(1:n_c))=  nele_c
            if(n_c .ge. 3)  cycle

            ID(iele)=  nele_c
            do j=iA_e2e(iele),iA_e2e(iele+1)-1
                ID(jA_e2e(j))   =  nele_c
            end do
        end do

!       we try to find circle again with less restrict condition.
        do i=1,nele
            iele=  jA(i)
            if(ID(iele) .gt. 0) cycle

!           try to find a circle in the graph.
            call get_circle(max_ele_circle, max_distance_circle, iele, n_c, c)
            if(n_c .ge. 3) then
                nele_c      =  nele_c+1
                ID(c(1:n_c))=  nele_c
                cycle
            end if
        end do

!       record the number of elements in every agglomeration.
        allocate(ag(nele_c), stat=err_mem)
        ag  =  0
        do iele=1,nele
            j   =  ID(iele)
            if(j .gt. 0)    ag(j)   =  ag(j)+1
        end do

!       attach the non-agglomerated element to its agglomerated neighbour.
        do i=1,nele
            iele=  jA(i)
            if(ID(iele) .gt. 0) cycle

            n_cd=  0
            do j=iA_e2e(iele),iA_e2e(iele+1)-1
                k   =  jA_e2e(j)
                if(ID(k) .le. 0)    cycle
                n_cd        =  n_cd+1
                cd(1,n_cd)  =  ID(k)
            end do
            if(n_cd .le. 0) stop 'Error: fails to find a neighbour that is coarsened.'
            call simplify_series(n_cd, 1, 2, cd)
            do j=1,n_cd
                cd(2,j) =  ag(cd(1,j))
            end do
            call iqsortcols(.true., 1, n_cd, 2, 2, cd)
            ag(cd(1,1)) =  ag(cd(1,1))+1
            ID(iele)= -cd(1,1)
        end do
        ID(1:nele)  =  abs(ID(1:nele))
        do i=1,nele
            if(ID(i) .le. 0)    stop 'Error: gmg_step_ag fails.'
        end do

        if(allocated(iA))   deallocate(iA)
        if(allocated(jA))   deallocate(jA)
        if(allocated(ag))   deallocate(ag)

        return
        contains
!           ------------------------------------------------------------------------
!           try to find a circle in the graph.
!           ------------------------------------------------------------------------
            subroutine get_circle(max_ele,max_d,seed,n_c,c)
            implicit none
            integer(dpI),intent(in):: max_ele,max_d,seed
            logical(dpL):: ltmp
            integer(dpI):: n_c,c(*),n_path,n_circle,R,iter,iele,i,j,ii,nb,circle(8,1000)

            path(1,1)   =  seed
            n_path      =  1
            n_circle    =  0
            circle      =  0
            do iter=2,max_ele+1
                if(n_path .le. 0)   exit

                R   =  n_path
                do i=1,n_path
                    iele=  path(iter-1,i)
                    do j=iA_e2e(iele),iA_e2e(iele+1)-1
                        nb  =  jA_e2e(j)

!                       iteself is not allowed.
                        if(nb .eq. iele)    cycle

!                       this element is already numbered.
                        if(ID(nb) .gt. 0)   cycle

!                       this element is too far away.
                        if(abs(ID_wf(nb)-ID_wf(path(1,i))) .gt. max_d)  cycle

!                       this element is already recorded.
                        ltmp=  .false.
                        do ii=2,iter-1
                            ltmp=  path(ii,i) .eq. nb
                            if(ltmp)    exit
                        end do
                        if(ltmp)    cycle

                        if((nb .eq. seed) .and. (iter .gt. 3)) then
!                           find a circle.
                            n_circle    =  n_circle+1
                            circle(1:iter-1,n_circle)   =  path(1:iter-1,i)
                        elseif((nb .eq. seed) .and. (iter .le. 3)) then
!                           do nothing.
                        else
                            R   =  R+1
                            path(1:iter-1,R)=  path(1:iter-1,i)
                            path(  iter  ,R)=  nb
                        end if
                    end do
                end do

                i       =  n_path
                n_path  =  0
                do j=i+1,R
                    n_path          =  n_path+1
                    path(:,n_path)  =  path(:,j)
                end do
            end do

            n_c =  0
            do j=1,1000
                ii  =  0
                do i=1,8
                    if(circle(i,j) .gt. 0)  ii  =  ii+1
                end do
                if(ii .gt. n_c) then
                    n_c     =  ii
                    c(1:n_c)=  circle(1:n_c,j)
                end if
            end do

            return
            end subroutine get_circle
        end subroutine gmg_step_ag
!       ------------------------------------------------------------------------
!       construct iA and jA from the finer mesh.
!       ------------------------------------------------------------------------
        subroutine gmg_get_coarse_topology(L1)
        use var_kind_def
        use var_global, only: err_mem
        implicit none
        integer(dpI),intent(in):: L1
        integer(dpI):: L0,ne2e,iele,i,j,k,L,R,nb,v(1000)
        integer(dpI),allocatable:: e2e(:,:)

        L0          =  L1-1
        mg(L1)%n_ele=  mg(L0)%n_ele_c
        allocate(mg(L1)%iA_e2e(mg(L1)%n_ele+1), stat=err_mem)
        allocate(e2e(2,size(mg(L0)%jA_e2e   )), stat=err_mem)
        ne2e=  0
        do iele=1,mg(L0)%n_ele
            L   =  mg(L0)%ID(iele)
            do j=mg(L0)%iA_e2e(iele),mg(L0)%iA_e2e(iele+1)-1
                R               =  mg(L0)%ID(mg(L0)%jA_e2e(j))
                ne2e            =  ne2e+1
                e2e(1:2,ne2e)   = (/L, R/)
            end do
        end do

        call iqsortcols(.true., 1, ne2e, 1, 2, e2e)

!       delete the self-connection.
        j   =  ne2e
        ne2e=  0
        do i=1,j
            if(e2e(1,i) .eq. e2e(2,i))  cycle
            ne2e            =  ne2e+1
            e2e(1:2,ne2e)   =  e2e(1:2,i)
        end do

!       setup iA.
        mg(L1)%iA_e2e   =  0
        i   =  1
        do while(i .le. ne2e)
            k   =  i
            do j=i+1,ne2e
                if(e2e(1,j) .eq. e2e(1,i)) then
                    k   =  j
                else
                    exit
                end if
            end do

            do j=i,k
                v(j-i+1)=  e2e(2,j)
            end do
            nb  =  k-i+1
            call simplify_series(nb, 1, 1, v)
            mg(L1)%iA_e2e(e2e(1,i)+1)   =  nb

            i   =  k+1
        end do

!       allocate memory for jA.
        j               =  0
        mg(L1)%iA_e2e(1)=  1
        do iele=1,mg(L1)%n_ele
            j   =  j+mg(L1)%iA_e2e(iele+1)
            mg(L1)%iA_e2e(iele+1)   =  j+1
        end do
        if(j .le. 0)    stop 'Error: fails to get the topology of the coarse level.'
        allocate(mg(L1)%jA_e2e(j), stat=err_mem)

        i   =  1
        do while(i .le. ne2e)
            k   =  i
            do j=i+1,ne2e
                if(e2e(1,j) .eq. e2e(1,i)) then
                    k   =  j
                else
                    exit
                end if
            end do

            do j=i,k
                v(j-i+1)=  e2e(2,j)
            end do
            nb  =  k-i+1
            call simplify_series(nb, 1, 1, v)

            iele=  e2e(1,i)
            do j=1,nb
                mg(L1)%jA_e2e(mg(L1)%iA_e2e(iele)+j-1)  =  v(j)
            end do

            i   =  k+1
        end do

        if(allocated(e2e))  deallocate(e2e)

        return
        end subroutine gmg_get_coarse_topology
!       ------------------------------------------------------------------------
!       adjust the element ordering, push the boundary element to the beginning.
!       ------------------------------------------------------------------------
        subroutine gmg_push_boundary_element(L0)
        use var_kind_def
        use var_global, only: err_mem
        use var_mesh
        implicit none
        integer(dpI),intent(in):: L0
        integer(dpI):: L1,i,j,isec,iele
        integer(dpI),allocatable:: ID(:)

        L1  =  L0+1
        allocate(ID(mg(L0)%n_ele_c), stat=err_mem)
        ID  =  0

        if(L0 .eq. 0) then
            i   =  0
            do isec=mesh(0)%sec_1,mesh(0)%sec_0
                if(.not. sec(isec)%is_int)  cycle

                do iele=1,sec(isec)%n_ele_b
                    j   =  iele+i
                    ID(mg(0)%ID(j)) = -1
                end do
                i   =  i+sec(isec)%n_ele
            end do
        else
            do iele=1,mg(L0)%n_ele_b
                iD(mg(L0)%ID(iele)) = -1
            end do
        end if

        j   =  0
        do i=1,mg(L0)%n_ele_c
            if(ID(i) .ge. 0)    cycle
            j       =  j+1
            ID(i)   =  j
        end do
        mg(L1)%n_ele_b  =  j
        do i=1,mg(L0)%n_ele_c
            if(ID(i) .ne. 0)    cycle
            j       =  j+1
            ID(i)   =  j
        end do
        if(j .ne. mg(L0)%n_ele_c)   stop 'Error: fails to push b_element to the beginning.'

        do i=1,mg(L0)%n_ele
            j           =  ID(mg(L0)%ID(i))
            if((j .le. 0) .or. (j .gt. mg(L0)%n_ele_c)) &
                &  stop 'Error: fails to push b_element to the beginning.'
            mg(L0)%ID(i)=  j
        end do

        if(allocated(ID))   deallocate(ID)

        return
        end subroutine gmg_push_boundary_element
!       ------------------------------------------------------------------------
!       mesh agglomeration destroy.
!       ------------------------------------------------------------------------
        subroutine mesh_agglomeration_circle_destroy
        implicit none
        integer(dpI):: lev

        do lev=0,10
            if(allocated(mg(lev)%ID  )) deallocate(mg(lev)%ID)
            if(allocated(mg(lev)%iA_e2e))   deallocate(mg(lev)%iA_e2e)
            if(allocated(mg(lev)%jA_e2e))   deallocate(mg(lev)%jA_e2e)
            if(allocated(mg(lev)%A     ))   deallocate(mg(lev)%A)
            if(allocated(mg(lev)%vol   ))   deallocate(mg(lev)%vol)
            if(allocated(mg(lev)%ID_wf ))   deallocate(mg(lev)%ID_wf)
            mg(lev)%n_ele   = -1
        end do

        return
        end subroutine mesh_agglomeration_circle_destroy
    end module var_agglomeration_circle
!-------------------------------------------------------------------------------
!   point-based agglomeration MG.
!-------------------------------------------------------------------------------
    module var_agglomeration_point
        use var_kind_def
        use var_global, only: err_mem,n_dim
        implicit none
        private
        public:: mesh_agglomeration_point,mesh_agglomeration_point_destroy
        public:: mg

        integer(dpI),allocatable:: perm(:),ordered_list(:)
        integer(dpI),allocatable:: e2e(:,:),id_v(:)

        type type_mg
            integer(dpI):: n_ele    =  0
            integer(dpI):: n_ele_c  =  0
            integer(dpI):: n_vtx    =  0
            integer(dpI),allocatable:: iA_e2e(:),jA_e2e(:)
            integer(dpI),allocatable:: iA_n2e(:),jA_n2e(:)
            integer(dpI),allocatable:: iA_e2n(:),jA_e2n(:)
            integer(dpI),allocatable:: ID_c(:)
        end type type_mg
        type(type_mg),save:: mg(0:9)

        contains
!       ------------------------------------------------------------------------
!       mesh agglomeration.
!       ------------------------------------------------------------------------
        subroutine mesh_agglomeration_point(mlev)
        implicit none
        integer(dpI),intent(in):: mlev
        integer(dpI):: lev,i
        integer(dpI),allocatable:: ID(:)

        call get_e2e_L0
        call get_n2e_e2n_L0
        allocate(ID(mg(0)%n_ele), stat=err_mem)
        forall(lev=1:mg(0)%n_ele)   ID(lev) =  lev

        do lev=0,mlev-2
            allocate(mg(lev)%ID_c(mg(lev)%n_ele))
            call point_aggregation(mg(lev)%n_ele, mg(lev)%iA_e2e, mg(lev)%jA_e2e, &
                &  mg(lev)%iA_e2n, mg(lev)%jA_e2n, mg(lev)%iA_n2e, mg(lev)%jA_n2e, &
                &  mg(lev)%n_ele_c, mg(lev)%ID_c)

            mg(lev+1)%n_ele =  mg(lev)%n_ele_c
            if(mg(lev)%n_ele_c .le. 25) exit

            call get_coarse_topology(lev+1)

            do i=1,mg(0)%n_ele
                ID(i)   =  mg(lev)%ID_c(ID(i))
            end do
        end do

!       do lev=0,mlev-1
!           print*,lev,'|',mg(lev)%n_ele
!       end do

!       call gmg_wr_ID_tec(ID)
!       stop

        if(allocated(ID))   deallocate(ID)

        return
        end subroutine mesh_agglomeration_point
!       ------------------------------------------------------------------------
!       get e2e of the finest level.
!       ------------------------------------------------------------------------
        subroutine get_e2e_L0
        use var_mesh
        implicit none
        integer(dpI):: v(200),isec,iele,j,k,m,im,s,e,c(200)

        mg(0)%n_ele =  0
        v           =  0
        j           =  0
        do isec=mesh(0)%sec_1,mesh(0)%sec_0
            if(.not. sec(isec)%is_int)  cycle
            v(isec)     =  mg(0)%n_ele+1
            mg(0)%n_ele =  mg(0)%n_ele+sec(isec)%n_ele
            do iele=1,sec(isec)%n_ele
                j   =  j+sec(isec)%iA_face_neighbour(iele+1) &
                    &   -sec(isec)%iA_face_neighbour(iele)+1
            end do
        end do
        allocate(mg(0)%iA_e2e(mg(0)%n_ele+1), stat=err_mem)
        allocate(mg(0)%jA_e2e(j            ), stat=err_mem)
        mg(0)%iA_e2e=  0
        mg(0)%jA_e2e=  0

        mg(0)%iA_e2e(1) =  1
        do isec=mesh(0)%sec_1,mesh(0)%sec_0
        if(.not. sec(isec)%is_int)  cycle
        do iele=1,sec(isec)%n_ele
            m   =  0
            do k=sec(isec)%iA_face_neighbour(iele),sec(isec)%iA_face_neighbour(iele+1)-1
                im  =  sec(isec)%jA_face_neighbour(k)
                if(im .gt. 0) then
                    s   =  mesh(0)%mortar_LR(3, im)
                    e   =  mesh(0)%mortar_LR(4, im)
                else
                    s   =  mesh(0)%mortar_LR(1,-im)
                    e   =  mesh(0)%mortar_LR(2,-im)
                end if
                if(.not. sec(s)%is_int) cycle
                m   =  m+1
                c(m)=  e+v(s)-1
            end do

            do j=1,m
                mg(0)%jA_e2e(j+mg(0)%iA_e2e(iele+v(isec)-1)-1)  =  c(j)
            end do
            mg(0)%iA_e2e(iele+v(isec))  =  mg(0)%iA_e2e(iele+v(isec)-1)+m
        end do
        end do

        return
        end subroutine get_e2e_L0
!       ------------------------------------------------------------------------
!       get n2e and e2n of the finest level.
!       ------------------------------------------------------------------------
        subroutine get_n2e_e2n_L0
        use var_mesh
        implicit none
        integer(dpI):: isec,iele,i,j,idx,v(8)

        i   =  0
        do isec=mesh(0)%sec_1,mesh(0)%sec_0
            if(.not. sec(isec)%is_int)  cycle

            i   =  i+sec(isec)%n_ele*sec(isec)%npe
        end do
        allocate(mg(0)%iA_n2e(mg(0)%n_ele+1), stat=err_mem)
        allocate(mg(0)%jA_n2e(i            ), stat=err_mem)
        allocate(mg(0)%jA_e2n(i            ), stat=err_mem)

        if(.not. allocated(id_v))   allocate(id_v(mesh(0)%n_vtx), stat=err_mem)
        id_v=  0
        do isec=mesh(0)%sec_1,mesh(0)%sec_0
            if(.not. sec(isec)%is_int)  cycle

            do iele=1,sec(isec)%n_ele
            do i=1,sec(isec)%npe
                id_v(sec(isec)%n2e(i,iele)) =  1
            end do
            end do
        end do
        mg(0)%n_vtx =  0
        do i=1,mesh(0)%n_vtx
            if(id_v(i) .le. 0)  cycle
            mg(0)%n_vtx =  mg(0)%n_vtx+1
            id_v(i)     =  mg(0)%n_vtx
        end do

        mg(0)%iA_n2e(1) =  1
        idx             =  1
        j               =  1
        do isec=mesh(0)%sec_1,mesh(0)%sec_0
            if(.not. sec(isec)%is_int)  cycle

            do iele=1,sec(isec)%n_ele
                do i=1,sec(isec)%npe
                    v(i)=  id_v(sec(isec)%n2e(i,iele))
                end do

                call ICOPY(sec(isec)%npe, v, 1, mg(0)%jA_n2e(idx), 1)
                idx             =  idx+sec(isec)%npe
                j               =  j+1
                mg(0)%iA_n2e(j) =  idx
            end do
        end do

        allocate(mg(0)%iA_e2n(mg(0)%n_vtx+1), stat=err_mem)
        call itrans_spy(.false., 1, mg(0)%n_ele, mg(0)%iA_n2e, mg(0)%jA_n2e, &
            &  mg(0)%n_vtx, mg(0)%iA_e2n, mg(0)%jA_e2n)

        return
        end subroutine get_n2e_e2n_L0
!       ------------------------------------------------------------------------
!       point-based aggregation.
!       ------------------------------------------------------------------------
        subroutine point_aggregation(N,iA_e2e,jA_e2e,iA_e2n,jA_e2n,iA_n2e,jA_n2e,n_c,ID_c)
        use var_kind_def
        implicit none
        integer(dpI),intent(in):: N,iA_e2e(*),jA_e2e(*),iA_e2n(*),jA_e2n(*), &
                            &     iA_n2e(*),jA_n2e(*)
        integer(dpI):: n_c,ID_c(*),i,j,k,L,M,seed,v(2,1000)

        if(N .le. 0)    stop 'Error: wrong input, pairwise aggregation.'

        if(.not. allocated(perm)) then
            allocate(perm        (N), stat=err_mem)
            allocate(ordered_list(N), stat=err_mem)
        end if
        call matrix_ordering(N, iA_e2e, jA_e2e, perm)
        do i=1,N
            ordered_list(perm(i))   =  i
        end do

        ID_c(1:N)   =  0
        n_c         =  0
        seed        =  1
        do while(.true.)
!           find the seed element.
            i   = -1
            do j=seed,N
                if(ID_c(ordered_list(j)) .le. 0) then
                    i   =  ordered_list(j)
                    seed=  j
                    exit
                end if
            end do
            if(i .le. 0)    exit

!           get the number of the surrouding element for each vertex of element-i.
            do k=iA_n2e(i),iA_n2e(i+1)-1
                j   =  jA_n2e(k)
                v(1,k-iA_n2e(i)+1)  =  j
                v(2,k-iA_n2e(i)+1)  =  0
                do L=iA_e2n(j),iA_e2n(j+1)-1
                    M   =  jA_e2n(L)
                    if(ID_c(M) .gt. 0)  cycle
                    v(2,k-iA_n2e(i)+1)  =  v(2,k-iA_n2e(i)+1)+1
                end do
            end do

            call iqsortcols(.false., 1, iA_n2e(i+1)-iA_n2e(i), 2, 2, v)
            n_c =  n_c+1
            j   =  v(1,1)
            do L=iA_e2n(j),iA_e2n(j+1)-1
                M   =  jA_e2n(L)
                if(ID_c(M) .gt. 0)  cycle
                ID_c(M) =  n_c
            end do
        end do

        return
        end subroutine point_aggregation
!       ------------------------------------------------------------------------
!       construct iA and jA from the finer mesh.
!       ------------------------------------------------------------------------
        subroutine get_coarse_topology(L1)
        use var_kind_def
        implicit none
        integer(dpI),intent(in):: L1
        logical(dpL):: ltmp
        integer(dpI):: L0,ne2e,iele,i,j,k,L,R,nb,v(1000),ivtx,c,idx

        L0          =  L1-1
        mg(L1)%n_ele=  mg(L0)%n_ele_c
        allocate(mg(L1)%iA_e2e(mg(L1)%n_ele+1), stat=err_mem)
        if(.not. allocated(e2e))    allocate(e2e(2,size(mg(L0)%jA_e2e)), stat=err_mem)
        ne2e=  0
        do iele=1,mg(L0)%n_ele
            L   =  mg(L0)%ID_C(iele)
            do j=mg(L0)%iA_e2e(iele),mg(L0)%iA_e2e(iele+1)-1
                R               =  mg(L0)%ID_C(mg(L0)%jA_e2e(j))
                if(L .eq. R)    cycle
                ne2e            =  ne2e+1
                e2e(1:2,ne2e)   = (/L, R/)
            end do
        end do
        call iqsortcols(.true., 1, ne2e, 1, 2, e2e)

!       setup iA_e2e.
        mg(L1)%iA_e2e   =  0
        i   =  1
        do while(i .le. ne2e)
            k   =  i
            do j=i+1,ne2e
                if(e2e(1,j) .eq. e2e(1,i)) then
                    k   =  j
                else
                    exit
                end if
            end do

            do j=i,k
                v(j-i+1)=  e2e(2,j)
            end do
            nb  =  k-i+1
            call simplify_series(nb, 1, 1, v)
            mg(L1)%iA_e2e(e2e(1,i)+1)   =  nb

            i   =  k+1
        end do

!       allocate memory for jA_e2e.
        j           =  0
        mg(L1)%iA_e2e(1)=  1
        do iele=1,mg(L1)%n_ele
            j   =  j+mg(L1)%iA_e2e(iele+1)
            mg(L1)%iA_e2e(iele+1)   =  j+1
        end do
        if(j .le. 0)    stop 'Error: fails to get the topology of the coarse level.'
        allocate(mg(L1)%jA_e2e(j), stat=err_mem)
        mg(L1)%jA_e2e   =  0

        i   =  1
        do while(i .le. ne2e)
            k   =  i
            do j=i+1,ne2e
                if(e2e(1,j) .eq. e2e(1,i)) then
                    k   =  j
                else
                    exit
                end if
            end do

            do j=i,k
                v(j-i+1)=  e2e(2,j)
            end do
            nb  =  k-i+1
            call simplify_series(nb, 1, 1, v)

            iele=  e2e(1,i)
            do j=1,nb
                mg(L1)%jA_e2e(mg(L1)%iA_e2e(iele)+j-1)  =  v(j)
            end do

            i   =  k+1
        end do

!       get the ID of vertex on the coarse level.
        mg(L1)%n_vtx=  0
        do ivtx=1,mg(L0)%n_vtx
            do j=mg(L0)%iA_e2n(ivtx),mg(L0)%iA_e2n(ivtx+1)-1
                v(j-mg(L0)%iA_e2n(ivtx)+1)  =  mg(L0)%ID_c(mg(L0)%jA_e2n(j))
            end do
            i   =  mg(L0)%iA_e2n(ivtx+1)-mg(L0)%iA_e2n(ivtx)
            j   =  i
            call simplify_series(j, 1, 1, v)
            if((i .gt. 1) .and. (j .le. 1)) then
                id_v(ivtx)  =  0
            else
                mg(L1)%n_vtx=  mg(L1)%n_vtx+1
                id_v(ivtx)  =  mg(L1)%n_vtx
            end if
        end do

        allocate(mg(L1)%iA_n2e(mg(L1)%n_ele+1), stat=err_mem)
        mg(L1)%iA_n2e   =  0
        i               =  0
        do iele=1,mg(L0)%n_ele
            c   =  mg(L0)%ID_c(iele)
            do j=mg(L0)%iA_n2e(iele),mg(L0)%iA_n2e(iele+1)-1
                if(id_v(mg(L0)%jA_n2e(j)) .le. 0)   cycle
                mg(L1)%iA_n2e(c+1)  =  mg(L1)%iA_n2e(c+1)+1
                i                   =  i+1
            end do
        end do
        allocate(mg(L1)%jA_n2e(i), stat=err_mem)
        mg(L1)%jA_n2e   =  0

        mg(L1)%iA_n2e(1)=  1
        do iele=1,mg(L1)%n_ele
            mg(L1)%iA_n2e(iele+1)   =  mg(L1)%iA_n2e(iele)+mg(L1)%iA_n2e(iele+1)
        end do

!       get iA_n2e and jA_n2e of coarse level.
        do iele=1,mg(L0)%n_ele
            c   =  mg(L0)%ID_c(iele)
            do j=mg(L0)%iA_n2e(iele),mg(L0)%iA_n2e(iele+1)-1
                k   =  id_v(mg(L0)%jA_n2e(j))
                if(k .le. 0)    cycle

                ltmp=  .true.
                do L=mg(L1)%iA_n2e(c),mg(L1)%iA_n2e(c+1)-1
                    if(mg(L1)%jA_n2e(L) .le. 0) then
                        mg(L1)%jA_n2e(L)=  k
                        ltmp            =  .false.
                        exit
                    end if
                end do
                if(ltmp)    stop 'Error: fails to get coarse iA_n2e.'
            end do
        end do

!       simplify the coarse level iA_n2e jA_n2e.
        idx =  1
        do iele=1,mg(L1)%n_ele
            L   =  mg(L1)%iA_n2e(iele)
            R   =  mg(L1)%iA_n2e(iele+1)
            i   =  R-L
            call ICOPY(i, mg(L1)%jA_n2e(L), 1, v, 1)
            call simplify_series(i, 1, 1, v)
            call ICOPY(i, v, 1, mg(L1)%jA_n2e(idx), 1)
            mg(L1)%iA_n2e(iele) =  i
            idx                 =  idx+i
        end do
        do iele=mg(L1)%n_ele,1,-1
            mg(L1)%iA_n2e(iele+1)   =  mg(L1)%iA_n2e(iele)
        end do
        mg(L1)%iA_n2e(1)=  1
        do iele=1,mg(L1)%n_ele
            mg(L1)%iA_n2e(iele+1)   =  mg(L1)%iA_n2e(iele)+mg(L1)%iA_n2e(iele+1)
        end do
        L   =  mg(L1)%iA_n2e(mg(L1)%n_ele+1)-1

        if(allocated(mg(L0)%iA_e2e))    deallocate(mg(L0)%iA_e2e)
        if(allocated(mg(L0)%jA_e2e))    deallocate(mg(L0)%jA_e2e)
        if(allocated(mg(L0)%iA_n2e))    deallocate(mg(L0)%iA_n2e)
        if(allocated(mg(L0)%jA_n2e))    deallocate(mg(L0)%jA_n2e)
        if(allocated(mg(L0)%iA_e2n))    deallocate(mg(L0)%iA_e2n)
        if(allocated(mg(L0)%jA_e2n))    deallocate(mg(L0)%jA_e2n)

        allocate(mg(L1)%iA_e2n(mg(L1)%n_vtx+1), stat=err_mem)
        allocate(mg(L1)%jA_e2n(L             ), stat=err_mem)
        call itrans_spy(.false., 1, mg(L1)%n_ele, mg(L1)%iA_n2e, mg(L1)%jA_n2e, &
            &  mg(L1)%n_vtx, mg(L1)%iA_e2n, mg(L1)%jA_e2n)
!       do i=1,mg(L1)%n_ele
!           print*,i,'|',mg(L1)%jA_n2e(mg(L1)%iA_n2e(i) : mg(L1)%iA_n2e(i+1)-1)
!       end do
!       do i=1,mg(L1)%n_vtx
!           print*,i,mg(L1)%jA_e2n(mg(L1)%iA_e2n(i) : mg(L1)%iA_e2n(i+1)-1)
!       end do

        return
        end subroutine get_coarse_topology
!       ------------------------------------------------------------------------
!       mesh agglomeration destroy.
!       ------------------------------------------------------------------------
        subroutine mesh_agglomeration_point_destroy
        implicit none
        integer(dpI):: lev

        do lev=0,9
            if(allocated(mg(lev)%iA_e2e))   deallocate(mg(lev)%iA_e2e)
            if(allocated(mg(lev)%jA_e2e))   deallocate(mg(lev)%jA_e2e)
            if(allocated(mg(lev)%iA_n2e))   deallocate(mg(lev)%iA_n2e)
            if(allocated(mg(lev)%jA_n2e))   deallocate(mg(lev)%jA_n2e)
            if(allocated(mg(lev)%iA_e2n))   deallocate(mg(lev)%iA_e2n)
            if(allocated(mg(lev)%jA_e2n))   deallocate(mg(lev)%jA_e2n)
            if(allocated(mg(lev)%ID_c  ))   deallocate(mg(lev)%ID_c)
            mg(lev)%n_ele   = -1
        end do
        if(allocated(perm        )) deallocate(perm)
        if(allocated(ordered_list)) deallocate(ordered_list)
        if(allocated(e2e         )) deallocate(e2e)
        if(allocated(id_v        )) deallocate(id_v)

        return
        end subroutine mesh_agglomeration_point_destroy
    end module var_agglomeration_point
!-------------------------------------------------------------------------------
!   mesh agglomeration.
!-------------------------------------------------------------------------------
    subroutine mesh_agglomeration
    use var_kind_def
    use var_global, only: is_has_cfg,cfg_file,err_mem
    use var_mg
    use var_agglomeration_circle, c=>mg
!   use var_agglomeration_pairwise, p2=>mg
    use var_agglomeration_point, p=>mg
    use var_parallel
    implicit none
    logical(dpL):: lbuf(1)
    integer(dpI):: ibuf(2),io_err,lev

!   ----------------------------------------------------------------------------
!   read and synchronize the GMG parameters.
    namelist /multigrid/    mlev,is_amg
    if(is_mesh_agglomerated)    return
    if(is_has_cfg .and. (myid .eq. 0)) then
        open(unit=10,file=trim(adjustl(cfg_file)))
        read(unit=10, nml=multigrid, iostat=io_err)
        if(io_err .gt. 0)   stop 'Error: fails to read namelist:multigrid.'
        close(10)
    end if
    lbuf(1) =  is_amg
    ibuf(1) =  mlev
    call mpi_bcast(lbuf, 2, mpi_dpL, 0, mpi_comm_world, mpi_err)
    call mpi_bcast(ibuf, 1, mpi_dpI, 0, mpi_comm_world, mpi_err)
    if(myid .ne. 0) then
        is_amg  =  lbuf(1)
        mlev    =  ibuf(1)
    end if
!   read and synchronize the GMG parameters.
!   ----------------------------------------------------------------------------

    mlev=  max(mlev, 1)
    if(mlev .le. 1) return

    allocate(mg(0:mlev-1))
    if(.true.) then
        call mesh_agglomeration_point(mlev)
        do lev=0,mlev-1
            mg(lev)%n_ele   =  p(lev)%n_ele
!           print*,lev,mg(lev)%n_ele

            if(lev .eq. mlev-1) cycle
            allocate(mg(lev)%ID(mg(lev)%n_ele), stat=err_mem)
            call ICOPY(mg(lev)%n_ele, p(lev)%ID_c, 1, mg(lev)%ID, 1)
        end do
        call mesh_agglomeration_point_destroy
    elseif(.false.) then
        call mesh_agglomeration_circle(mlev)
        do lev=0,mlev-1
            mg(lev)%n_ele   =  c(lev)%n_ele

            if(lev .eq. mlev-1) cycle
            allocate(mg(lev)%ID(mg(lev)%n_ele), stat=err_mem)
            call ICOPY(mg(lev)%n_ele, c(lev)%ID, 1, mg(lev)%ID, 1)
        end do
        call mesh_agglomeration_circle_destroy
    elseif(.false.) then
!       call mesh_agglomeration_pairwise(mlev)
!       do lev=0,mlev-1
!           mg(lev)%n_ele   =  ap(n_dim*lev)%n_ele

!           if(lev .eq. mlev-1) cycle
!           allocate(mg(lev)%ID(mg(lev)%n_ele), stat=err_mem)
!           forall(i=1:mg(lev)%n_ele)   mg(lev)%ID(i)   =  i
!           do j=0,n_dim-1
!           do i=1,mg(lev)%n_ele
!               mg(lev)%ID(i)   =  ap(j+n_dim*lev)%ID_c(mg(lev)%ID(i))
!           end do
!           end do
!       end do
!       call mesh_agglomeration_destroy
    end if
    is_mesh_agglomerated=  .true.

    return
    end subroutine mesh_agglomeration
