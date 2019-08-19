!-------------------------------------------------------------------------------
!   read in topology, CGNS format.
!-------------------------------------------------------------------------------
    subroutine read_mesh_cgns
    use var_kind_def
    use var_global
    use var_cgns
    use var_load_balance, xyz=>xyz_parallel
    use var_mesh
    use var_parallel
    use var_temp, only: ibf1
    implicit none
    character(len=80):: str
    logical(kind=4):: ltmp,ltm1,ltm2
    integer(kind=4):: isec,isize(6),itype,isec_ele,ele1,ele0,ip,ele1_sec,ele0_sec, &
                    & ib,i,j,k,nnz,L,M
    real   (dpR),allocatable:: x(:,:)

    if(myid .eq. 0) then
        inquire(file=mesh_file, exist=ltmp)
        if(.not. ltmp)  stop 'Error: cgns mesh file not found.'

        call cg_open_f(mesh_file, CG_MODE_READ, ifile, cg_err)
        call cg_nbases_f(ifile, n_base, cg_err)
        if(n_base .ne. 1)   stop 'Error: only 1 base in the mesh supported.'
        call cg_base_read_f(ifile, 1, str, n_dim, itype, cg_err)

        call cg_nzones_f(ifile, 1, n_zone, cg_err)
        if(n_zone .ne. 1)   stop 'Error: only 1 zone in the mesh supported.'
        call cg_zone_read_f(ifile, 1, 1, zone_name, isize, cg_err)
        n_vtx_g =  isize(1)

        call cg_nsections_f(ifile, 1, 1, isize, cg_err)
        n_sec_g =  isize(1)
        allocate(sec_g(5,n_sec_g))

        allocate(sec_name_g(n_sec_g))
        do isec=1,n_sec_g
            call cg_section_read_f(ifile, 1, 1, isec, sec_name_g(isec), itype, &
                &  isize(2), isize(3), isize(4), isize(5), cg_err)

            sec_g(1,isec)   =  itype
            sec_g(3,isec)   =  isize(2)
            sec_g(4,isec)   =  isize(3)
            sec_g(5,isec)   =  isize(5)
            call cg_npe_f(itype, sec_g(2,isec), cg_err)
        end do
    end if

    isize(1)=  n_dim
    isize(2)=  n_sec_g
    isize(3)=  n_vtx_g
    call mpi_bcast(isize, 3, mpi_dpI, 0, mpi_comm_world, mpi_err)
    n_dim   =  isize(1)
    n_sec_g =  isize(2)
    n_vtx_g =  isize(3)
    is_2d_cal   =  n_dim .eq. 2
    if(myid .ne. 0) allocate(sec_g(5,n_sec_g))
    call mpi_bcast(sec_g, 5*n_sec_g, mpi_dpI, 0, mpi_comm_world, mpi_err)

!   ----------------------------------------------------------------------------
!   cal number and nnz of boundary elements.
    n_eleb_g    =  0
    nnz_n2e_b   =  0
    n_seci_g    =  0
    n_elei_g    =  1
    n_ele_g     =  0
    allocate(seci_g(7,n_sec_g))
    if(myid .eq. 0) then
        allocate(sec_name (n_sec_g))
        allocate(bsec_name(n_sec_g))
    end if
    i   =  0
    do isec=1,n_sec_g
        itype   =  sec_g(1,isec)
        if(is_2d_cal) then
            ltmp= (itype .eq. BAR_2) .or. (itype .eq. BAR_3)
        else
            ltmp= (itype .eq. TRI_3) .or. (itype .eq. QUAD_4) .or. &
                & (itype .eq. QUAD_8) .or. (itype .eq. QUAD_9)
        end if
        n_ele_g =  n_ele_g+sec_g(4,isec)-sec_g(3,isec)+1
        if(ltmp) then
            n_eleb_g    =  n_eleb_g+sec_g(4,isec)-sec_g(3,isec)+1
            nnz_n2e_b   =  nnz_n2e_b+sec_g(2,isec)*(sec_g(4,isec)-sec_g(3,isec)+1)
            i           =  i+1
            if(myid .eq. 0) bsec_name(i)=  sec_name_g(isec)
        else
            n_seci_g            =  n_seci_g+1
            seci_g(1:4,n_seci_g)=  sec_g(1:4,isec)
            seci_g(5  ,n_seci_g)=  isec
            seci_g(6  ,n_seci_g)=  n_elei_g
            seci_g(7  ,n_seci_g)=  n_elei_g+sec_g(4,isec)-sec_g(3,isec)
            n_elei_g            =  seci_g(7,n_seci_g)+1
            if(myid .eq. 0) sec_name(n_seci_g)=  sec_name_g(isec)
        end if
    end do
    n_elei_g=  n_elei_g-1
!   cal number and nnz of boundary elements.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   distribute cell elements to all processes.
    allocate(elei_prc(2,0:nprc-1))

!   here we have to determine the size of jA_n2c.
    ele1=  1
    k   =  0
    do ip=0,nprc-1
        ele0=  ele1+n_elei_g/nprc-1
        if(ip .eq. nprc-1)  ele0=  n_elei_g

!       myid should save the information of cell(ele1:ele0)
        nnz =  0
        do isec=1,n_seci_g
            if(seci_g(6,isec) .gt. ele0)    cycle
            if(seci_g(7,isec) .lt. ele1)    cycle

            i   =  max(ele1, seci_g(6,isec))
            j   =  min(ele0, seci_g(7,isec))
            nnz =  nnz+(j-i+1)*seci_g(2,isec)
        end do
        if(myid .eq. 0 )    k   =  max(k, nnz)
        if(myid .eq. ip) then
            allocate(iA_n2e  (ele0-ele1+2), stat=err_mem)
            allocate(ele_type(ele0-ele1+1), stat=err_mem)
            allocate(jA_n2e  (nnz        ), stat=err_mem)
        end if

        elei_prc(1,ip)  =  ele1
        elei_prc(2,ip)  =  ele0
        ele1=  ele0+1
    end do
    if(myid .eq. 0) allocate(ibf1(k), stat=err_mem)

    ele1        =  1
    n_ele_l     =  0
    iA_n2e(1)   =  1
    do ip=0,nprc-1
        ele0=  ele1+n_elei_g/nprc-1
        if(ip .eq. nprc-1)  ele0=  n_elei_g

!       myid should save the information of cell(ele1:ele0)
        nnz =  0
        do isec_ele=1,n_seci_g
            if(seci_g(6,isec_ele) .gt. ele0)    cycle
            if(seci_g(7,isec_ele) .lt. ele1)    cycle

            ele1_sec=  max(ele1, seci_g(6,isec_ele))
            ele0_sec=  min(ele0, seci_g(7,isec_ele))
            isec    =  seci_g(5,isec_ele)

            if(myid .eq. 0) then
                i   =  ele1_sec-seci_g(6,isec_ele)+sec_g(3,isec)
                j   =  ele0_sec-seci_g(7,isec_ele)+sec_g(4,isec)
                call cg_elements_partial_read_f(ifile, 1, 1, isec, i, j, &
                    &  ibf1(nnz+1), CG_Null, cg_err)
            end if
            if(myid .eq. ip) then
                j   =  n_ele_l+1
                k   =  n_ele_l+ele0_sec-ele1_sec+1
                do i=j,k
                    iA_n2e(i+1) =  iA_n2e(i)+seci_g(2,isec_ele)
                end do
                forall(i=j:k)   ele_type(i) =  seci_g(1,isec_ele)
                n_ele_l =  n_ele_l+ele0_sec-ele1_sec+1
            end if
            nnz =  nnz+seci_g(2,isec_ele)*(ele0_sec+1-ele1_sec)
        end do

        if(myid .eq. 0) then
            if(ip .eq. 0) then
                jA_n2e(1:nnz)   =  ibf1(1:nnz)
            else
                call mpi_send(ibf1, nnz, mpi_dpI, ip, ip, mpi_comm_world, mpi_err)
            end if
        else
            if(myid .eq. ip) then
                call mpi_recv(jA_n2e, nnz, mpi_dpI, 0 , ip, mpi_comm_world, &
                    &  mpi_status, mpi_err)
            end if
        end if

        ele1=  ele0+1
    end do
!   distribute cell elements to all processes.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   distribute xyz to all processes.
    allocate(vtx_prc(2,0:nprc-1))
    j   =  n_vtx_g/nprc
    do ip=0,nprc-1
        vtx_prc(1:2,ip) = (/1+ip*j, min((ip+1)*j, n_vtx_g)/)
    end do
    vtx_prc(2,nprc-1)   =  n_vtx_g
    j   =  0
    do ip=0,nprc-1
        j   =  max(j, vtx_prc(2,ip)-vtx_prc(1,ip)+1)
    end do
    allocate(xyz(n_dim*j), stat=err_mem)

    do ip=nprc-1,0,-1
        i   =  vtx_prc(1,ip)
        j   =  vtx_prc(2,ip)
        L   =  j-i+1

        if(myid .eq. 0) then
            call cg_coord_read_f(ifile, 1, 1, 'CoordinateX', RealDouble, &
                &  i, j, xyz       , cg_err)
            call cg_coord_read_f(ifile, 1, 1, 'CoordinateY', RealDouble, &
                &  i, j, xyz(  L+1), cg_err)
            if(.not. is_2d_cal) then
            call cg_coord_read_f(ifile, 1, 1, 'CoordinateZ', RealDouble, &
                &  i, j, xyz(2*L+1), cg_err)
            end if

!           offset the Coordinates if needed.
            do m=1,n_dim
            do k=1,L
                xyz(k+L*(m-1))  =  zoom*xyz(k+L*(m-1))+offset(m)
            end do
            end do

            if(ip .ne. myid) then
                call mpi_send(xyz, n_dim*L, mpi_dpR, ip, ip, mpi_comm_world, mpi_err)
            end if
        else
            if(myid .eq. ip)    call mpi_recv(xyz, n_dim*L, mpi_dpR, 0, myid, &
                &  mpi_comm_world, mpi_status, mpi_err)
        end if
    end do
    k   =  vtx_prc(2,myid)-vtx_prc(1,myid)+1
    allocate(x(n_dim,k), stat=err_mem)
    do i=1,k
    do j=1,n_dim
        x(j,i)  =  xyz((j-1)*k+i)
    end do
    end do
    call DCOPY(n_dim*k, x, 1, xyz, 1)
    if(allocated(x))    deallocate(x)
!   distribute xyz to all processes.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   def the n2c info of boundary elements.
    allocate(iA_n2e_b(  n_eleb_g+1), stat=err_mem)
    allocate(bct     (2,n_eleb_g  ), stat=err_mem)
    allocate(jA_n2e_b(  nnz_n2e_b ), stat=err_mem)
    iA_n2e_b(1) =  1
    nnz =  1
    ele1=  1
    do isec=1,n_sec_g
        itype   =  sec_g(1,isec)
        if(is_2d_cal) then
            ltmp= (itype .eq. BAR_2) .or. (itype .eq. BAR_3)
        else
            ltmp= (itype .eq. TRI_3) .or. (itype .eq. QUAD_4) .or. &
                & (itype .eq. QUAD_8) .or. (itype .eq. QUAD_9)
        end if
        if(.not. ltmp)  cycle

        ele0=  ele1+sec_g(4,isec)-sec_g(3,isec)
        do i=ele1,ele0
            iA_n2e_b(i+1)   =  iA_n2e_b(i)+sec_g(2,isec)
        end do
        if(myid .eq. 0) then
            if(sec_g(5,isec) .eq. 0) then
                call cg_elements_read_f(ifile, 1, 1, isec, jA_n2e_b(nnz), cg_null, cg_err)
            else
                call cg_elements_read_f(ifile, 1, 1, isec, jA_n2e_b(nnz), ibf1, cg_err)
            end if
        end if

        nnz =  nnz+sec_g(2,isec)*(sec_g(4,isec)-sec_g(3,isec)+1)
        ele1=  ele0+1
    end do
    call mpi_bcast(jA_n2e_b, nnz_n2e_b, mpi_dpI, 0, mpi_comm_world, mpi_err)
!   def the n2c info of boundary elements.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   read the boundary conditions.
    if(myid .eq. 0) then
    call cg_nbocos_f(ifile, 1, 1, n_bocos, cg_err)
    j   =  0
    if(n_bocos .gt. 0)  allocate(boconame(  n_bocos), stat=err_mem)
    if(n_bocos .gt. 0)  allocate(bocoinfo(4,n_bocos), stat=err_mem)
    do ib=1,n_bocos
        call cg_goto_f(ifile, 1, cg_err, 'Zone_t', 1, 'ZoneBC_t', 1, 'BC_t', ib, 'end')
        call cg_gridlocation_read_f(igr, cg_err)

        call cg_boco_info_f(ifile, 1, 1, ib, boconame(ib), ibocotype, iptset, &
            &  npts, normalindex, normallistflag, normaldatatype, ndataset, cg_err)
        if(iptset .ne. ElementRange)    stop 'Error: bc def with ElementRange prefered.'
        call cg_boco_read_f(ifile, 1, 1, ib, isize, normallist, cg_err)

        ltm1= (ibocotype .eq. BCWall) .or. &
            & (ibocotype .eq. BCWallInviscid) .or. &
            & (ibocotype .eq. BCSymmetryPlane) .or. &
            & (ibocotype .eq. BCExtrapolate)
        ltm2= (ibocotype .eq. BCWallViscous) .or. &
            & (ibocotype .eq. BCWallViscousHeatFlux) .or. &
            & (ibocotype .eq. BCWallViscousIsothermal)
        if(ltm1)    ibocotype   =  BCWallInviscid
!       if((.not. is_vis_cal) .and. ltm2)   ibocotype   =  BCWallInviscid

        bocoinfo(1  ,ib)=  ibocotype
        bocoinfo(2:3,ib)=  isize(1:2)
        do i=isize(1),isize(2)
            j           =  j+1
            bct(1:2,j)  = (/ib, ibocotype/)
        end do

        ltmp=  .false.
        do isec=1,n_sec_g
            ltmp= (bocoinfo(2,ib) .eq. sec_g(3,isec)) .and. &
                & (bocoinfo(3,ib) .eq. sec_g(4,isec))
            if(ltmp)    exit
        end do
        if(.not. ltmp)  stop 'Error: fails to map boundary to section.'
        bocoinfo(4,ib)  =  isec
    end do
    call cg_close_f(ifile, cg_err)
    end if
    call mpi_bcast(n_bocos, 1, mpi_dpI, 0, mpi_comm_world, mpi_err)
    if(n_bocos .gt. 0) then
        if(.not. allocated(bocoinfo))   allocate(bocoinfo(4,n_bocos), stat=err_mem)
        call mpi_bcast(bocoinfo, 4*n_bocos, mpi_dpI, 0, mpi_comm_world, mpi_err)
    end if
    call mpi_bcast(bct, 2*n_eleb_g, mpi_dpI, 0, mpi_comm_world, mpi_err)
!   read the boundary conditions.
!   ----------------------------------------------------------------------------

    return
    end subroutine read_mesh_cgns
!-------------------------------------------------------------------------------
!   read in periodic boundary, CGNS format.
!-------------------------------------------------------------------------------
    subroutine get_vtx_pair_per
    use var_kind_def
    use var_cgns
    use var_global
    use var_load_balance
    use var_parallel
    use var_per_bnd
    implicit none
    integer(kind=4),parameter:: max_npnts=120000
    character(len=80):: con_name,donorname
    integer(kind=4):: n_cons,location,connect_type,ptset_type,npnts,donor_zonetype, &
                    & donor_ptset_type,donor_datatype,ndata_donor,pnts(max_npnts), &
                    & donor_data(max_npnts),iper,io_err
    logical(dpL):: ltmp
    integer(dpI):: icon,i,j,k,m,ivtx,n_pair,pair(0:7),n_n2n_per,per_path(3,7),L,R,n_vtx
    real   (kind=4):: c(3),a(3),t(3)
    real   (dpR):: z(3),Translation(3,100)
    integer(dpI),allocatable:: n2n_per(:,:),vtx(:)

    namelist /per_bnd/  Translation
    c           =  0.0E0
    a           =  0.0E0
    t           =  0.0E0
    Translation = -1.0d99
    n_vtx       =  size(jA_n2e)
    allocate(vtx(n_vtx))
    call ICOPY(n_vtx, jA_n2e, 1, vtx, 1)
    call simplify_series(n_vtx, 1, 1, vtx)

!   ----------------------------------------------------------------------------
!   read and synchronize the periodic boundary condition.
    if(myid .eq. 0) then
        call cg_open_f(mesh_file, CG_MODE_READ, ifile, cg_err)
        call cg_nconns_f(ifile, 1, 1, n_cons, cg_err)
        if(mod(n_cons,2) .ne. 0) stop 'Error: failed to read per_bnd.'
        n_per   =  0
        do icon=1,n_cons
            call cg_conn_periodic_read_f(ifile, 1, 1, icon, c, a, t, cg_err)
            if(cg_err .eq. 0)   n_per   =  n_per+1
        end do
        if(mod(n_per,2) .ne. 0) stop 'Error: failed to read per_bnd.'
        allocate(per(n_per))

        n_per   =  0
        do icon=1,n_cons
            call cg_conn_periodic_read_f(ifile, 1, 1, icon, c, a, t, cg_err)
            if(cg_err .ne. 0)   cycle
            n_per   =  n_per+1

            call cg_conn_info_f(ifile, 1, 1, icon, con_name, location, connect_type, &
                &  ptset_type, npnts, donorname, donor_zonetype, donor_ptset_type, &
                &  donor_datatype, ndata_donor, cg_err)
            call cg_conn_read_f(ifile,1,1,icon, pnts, donor_datatype, donor_data, cg_err)

            if(ptset_type .eq. PointRange) then
                npnts   =  pnts(2)-pnts(1)
                per(n_per)%nvtx =  npnts
            else
                per(n_per)%nvtx =  npnts
            end if
            if(ndata_donor .ne. npnts)  stop 'Error: failed to read per_bnd.'
            if(npnts .gt. max_npnts)    stop 'Error: max_npnts too small.'

            allocate(per(n_per)%vtx(2,npnts))
            if(ptset_type .eq. PointRange) then
                forall(i=1:npnts)  per(n_per)%vtx(1,i)  =  pnts(1)+i-1
            else
                forall(i=1:npnts)  per(n_per)%vtx(1,i)  =  pnts(i)
            end if
            per(n_per)%vtx(2,1:npnts)   =  donor_data(1:npnts)
            call iqsortcols(.true., 1, npnts, 1, 2, per(n_per)%vtx)

            if(maxval(abs(real(a, kind=dpR))) .ge. 1.0d-10) then
                per(n_per)%per_type         =  2
                per(n_per)%rotationcenter   =  real(c, kind=dpR)
                z   =  real(a, kind=dpR)
                do i=1,3
                    if(abs(z(i)) .ge. 1.0d-5) then
                        if(z(i) .gt.  pi)   z(i)=  z(i)-2.0d0*pi
                        if(z(i) .lt. -pi)   z(i)=  z(i)+2.0d0*pi
                        per(n_per)%rotationangle(i) =  nint(2.0d0*pi/z(i), dpI)
                    else
                        per(n_per)%rotationangle(i) =  0
                    end if
                end do
            elseif(maxval(abs(real(t, kind=dpR))) .gt. 1.0d-10) then
                per(n_per)%per_type     =  1
                per(n_per)%translation  =  real(t, kind=dpR)
            else
                stop 'Error: failed to read per_bnd.'
            end if
        end do
        call cg_close_f(ifile, cg_err)

        if(.false.) then
            do iper=1,n_per
                write(unit=6,fmt='(2I4,3ES20.12,3I6)'),iper,per(iper)%per_type, &
                    &  per(iper)%translation(1:3),per(iper)%rotationangle(1:3)
            end do
            stop
        end if

        if(is_has_cfg) then
            open(unit=10,file=trim(adjustl(cfg_file)))
            read(unit=10,nml=per_bnd,iostat=io_err)
            if(io_err .gt. 0)   stop 'Error: fails to read namelist:per_bnd.'
            close(10)

            do iper=1,n_per
                if(per(iper)%per_type .ne. 1)   cycle
                if(Translation(1,iper) .le. -1.0d20)    cycle
                per(iper)%translation(1:3)  =  Translation(1:3,iper)
            end do
        end if
    end if
    call mpi_bcast(n_per, 1, mpi_dpI, 0, mpi_comm_world, mpi_err)
    if(n_per .le. 0)    return

    if(myid .eq. 0) then
        n_vtx_pair_per  =  0
        do i=1,n_per
            n_vtx_pair_per  =  n_vtx_pair_per+per(i)%nvtx
        end do
        allocate(vtx_pair_per(3,n_vtx_pair_per))

!       here we try to find the same periodic definition.
        n_per_info  =  0
        pnts        =  0
        do j=1,n_per
        ltmp=  .false.
        do i=1,j-1
            if(per(i)%per_type .ne. per(j)%per_type)    cycle
            if(maxval(abs(per(i)%rotationcenter-per(j)%rotationcenter)) .gt. 0.0d0) cycle
            if(maxval(abs(per(i)%rotationangle -per(j)%rotationangle )) .gt. 0.0d0) cycle
            if(maxval(abs(per(i)%translation   -per(j)%translation   )) .gt. 0.0d0) cycle
            ltmp=  .true.
            exit
        end do

        if(ltmp) then
            pnts(j) =  pnts(i)
        else
            n_per_info  =  n_per_info+1
            pnts(j)     =  n_per_info
        end if
        end do

        allocate(per_info(10,n_per_info))
        per_info=  0.0d0
        do i=1,n_per
            per_info(1   ,pnts(i))  =  real(per(i)%per_type, kind=dpR)
            per_info(2:4 ,pnts(i))  =  per(i)%rotationcenter
            if(per(i)%rotationangle(1) .ne. 0)  per_info(5,pnts(i)) =  &
                &  2.0d0*pi/real(per(i)%rotationangle(1), kind=dpR)
            if(per(i)%rotationangle(2) .ne. 0)  per_info(6,pnts(i)) =  &
                &  2.0d0*pi/real(per(i)%rotationangle(2), kind=dpR)
            if(per(i)%rotationangle(3) .ne. 0)  per_info(7,pnts(i)) =  &
                &  2.0d0*pi/real(per(i)%rotationangle(3), kind=dpR)
            per_info(8:10,pnts(i))  =  per(i)%translation
        end do

        j   =  0
        do i=1,n_per
        do k=1,per(i)%nvtx
            j   =  j+1
            vtx_pair_per(1:3,j) = (/pnts(i), per(i)%vtx(1:2,k)/)
        end do
        end do
    end if

    call mpi_bcast(n_per_info, 1, mpi_dpI, 0, mpi_comm_world, mpi_err)
    if(myid .ne. 0) allocate(per_info(10,n_per_info))
    call mpi_bcast(per_info, 10*n_per_info, mpi_dpR, 0, mpi_comm_world, mpi_err)
    allocate(per_mat(9,n_per_info))
    do i=1,n_per_info
        if(nint(per_info(1,i)) .eq. 2) then
            call cal_rot_mat(per_info(5,i), per_mat(1,i))
        else
            per_mat(1:9,i)  = (/1.0d0,0.0d0,0.0d0, 0.0d0,1.0d0,0.0d0, 0.0d0,0.0d0,1.0d0/)
        end if
    end do

    call mpi_bcast(n_vtx_pair_per, 1, mpi_dpI, 0, mpi_comm_world, mpi_err)
    if(myid .ne. 0) allocate(vtx_pair_per(3,n_vtx_pair_per))
    call mpi_bcast(vtx_pair_per,3*n_vtx_pair_per, mpi_dpI, 0, mpi_comm_world, mpi_err)
!   read and synchronize the periodic boundary condition.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   find out all the periodic vertex pair.
    call iqsortcols(.true., 1, n_vtx_pair_per, 2, 3, vtx_pair_per)

    i   =  1
    do while(i .le. n_vtx_pair_per)
        k   =  i
        do j=i+1,n_vtx_pair_per
            if(vtx_pair_per(2,j) .ne. vtx_pair_per(2,i)) then
                exit
            else
                k   =  j
            end if
        end do
        if(k .gt. i)    call iqsortcols(.true., i, k, 3, 3, vtx_pair_per)

        do j=i+1,k
            if(abs(vtx_pair_per(3,j)) .eq. abs(vtx_pair_per(3,j-1))) then
!               between two vertex, only one periodic mapping is allowed.
                vtx_pair_per(1:3,j) = -abs(vtx_pair_per(1:3,j))
            end if
        end do

        i   =  k+1
    end do

!   between two vertex, only one periodic mapping is allowed. So delete redundant info.
    j   =  n_vtx_pair_per
    n_vtx_pair_per  =  0
    do i=1,j
        if(vtx_pair_per(1,i) .le. 0)    cycle
        n_vtx_pair_per  =  n_vtx_pair_per+1
        vtx_pair_per(1:3,n_vtx_pair_per)=  vtx_pair_per(1:3,i)
    end do

    i   =  0
    do j=1,n_vtx_pair_per
        call ib_search(1, n_vtx, 1, 1, vtx, vtx_pair_per(2,j), ltmp, L, R)
        if(.not. ltmp)  cycle
        i   =  i+1
    end do
    allocate(n2n_per(5,i*2+100))

!   my vertex is periodically connected to another vertex.
    n_n2n_per   =  0
    do k=1,n_vtx_pair_per
        ivtx=  vtx_pair_per(2,k)
        call ib_search(1, n_vtx, 1, 1, vtx, ivtx, ltmp, L, R)
        if(.not. ltmp)  cycle
        n_pair  =  0
        per_path=  0
        pair(0) =  ivtx
        call get_pair(ivtx, 0, 0, 0)

        do j=1,n_pair
            n_n2n_per   =  n_n2n_per+1
            if(n_n2n_per .gt. size(n2n_per,2))  stop 'Error: n2n_per too small.'
            n2n_per(1:2,n_n2n_per)  = (/ivtx, pair(j)/)
            n2n_per(3:5,n_n2n_per)  =  per_path(1:3,j)
        end do
    end do
    if(allocated(vtx_pair_per)) deallocate(vtx_pair_per)
!   find out all the periodic vertex pair.
!   ----------------------------------------------------------------------------

    call mpi_allreduce(n_n2n_per,n_vtx_pair_per,1,mpi_dpI,mpi_sum,mpi_comm_world,mpi_err)
    allocate(vtx_pair_per(5,n_vtx_pair_per), stat=err_mem)
    call i_allgather(mpi_comm_world, nprc, 5*n_n2n_per, n2n_per, vtx_pair_per)
    if(allocated(n2n_per))  deallocate(n2n_per)

    call iqsortcols(.true., 1, n_vtx_pair_per, 1, 5, vtx_pair_per)
    i   =  1
    do while(i .le. n_vtx_pair_per)
        k   =  i
        do j=i+1,n_vtx_pair_per
            if(vtx_pair_per(1,j) .eq. vtx_pair_per(1,i)) then
                k   =  j
            else
                exit
            end if
        end do
        call iqsortcols(.true., i, k, 2, 5, vtx_pair_per)
        do j=i+1,k
            ltmp=  .false.
            do m=i,j-1
                if(vtx_pair_per(1,m) .lt. 0)    cycle
                ltmp= (vtx_pair_per(3,j) .eq. vtx_pair_per(3,m)) .and. &
                    & (vtx_pair_per(4,j) .eq. vtx_pair_per(4,m)) .and. &
                    & (vtx_pair_per(5,j) .eq. vtx_pair_per(5,m))
                if(ltmp)    exit
            end do
            if(ltmp)    vtx_pair_per(1:2,j) = -abs(vtx_pair_per(1:2,j))
        end do

        i   =  k+1
    end do
    m               =  n_vtx_pair_per
    n_vtx_pair_per  =  0
    do i=1,m
        if(vtx_pair_per(1,i) .le. 0)    cycle
        n_vtx_pair_per  =  n_vtx_pair_per+1
        vtx_pair_per(1:5,n_vtx_pair_per)=  vtx_pair_per(1:5,i)
    end do

    if(.not. allocated(per))    return
    do i=1,n_per
        if(allocated(per(i)%vtx))   deallocate(per(i)%vtx)
    end do
    if(allocated(per))  deallocate(per)
    if(allocated(vtx))  deallocate(vtx)

    return
    contains
!       ------------------------------------------------------------------------
!       recursively find all the R part of vertex pair.
!       ------------------------------------------------------------------------
        recursive subroutine get_pair(root,dep,i1,i0)
        implicit none
        integer(dpI),intent(in):: root,dep,i1,i0
        logical(dpL):: ltmp
        integer(dpI):: RR,L,R,i,j,k,dep_new,j1,j0

        if(i1 .gt. i0)  return
        dep_new =  dep+1

        j1  =  i0+1
        j0  =  i0
        do i=i1,i0
            call ib_search(1, n_vtx_pair_per, 2, 3, vtx_pair_per, pair(i), ltmp, L, R)
            if(.not. ltmp)  cycle

            do k=L,R
                RR  =  vtx_pair_per(3,k)
                if(RR .eq. root)    cycle

                ltmp=  .false.
                do j=1,n_pair
                    ltmp=  pair(j) .eq. RR
                    if(ltmp)    exit
                end do
!               this vtx_pair_per has already been recorded.
                if(ltmp)    cycle

                if(dep_new .gt. n_dim)  stop 'Error: periodic path depth larger than 3.'
                n_pair  =  n_pair+1
                j0      =  j0+1
                pair(j0)=  RR
                if(i .gt. 0)    per_path(1:dep,n_pair) =  per_path(1:dep,i)
                per_path(1+dep,n_pair) =  vtx_pair_per(1,k)
                call iqsort(.true., 1, dep_new, per_path(1,n_pair))
            end do
        end do
        if((dep_new .eq. 1) .and. (.not. is_include_corner))    return
        if(j0 .ge. j1)  call get_pair(root, dep_new, j1, j0)

        return
        end subroutine get_pair
    end subroutine get_vtx_pair_per
!-------------------------------------------------------------------------------
!   output local mesh.
!-------------------------------------------------------------------------------
    subroutine wr_local_mesh
    use var_kind_def
    use var_cgns
    use var_global, only: is_2d_cal,n_dim,err_mem
    use var_mesh, only: sec,mesh
    use var_parallel
    implicit none
    character(len=80):: str
    character(len=1000):: tec
    logical(dpL):: ltmp
    integer(dpI):: isec,iele,npe,n_vtx,i,j,L,R,ele,v(8)
    integer(dpI),allocatable:: n2e(:),vtx_ID(:)
    real   (dpR),allocatable:: xyz(:,:)

    j   =  0
    do isec=mesh(0)%sec_1,mesh(0)%sec_0
        npe =  npe_ele_1st(sec(isec)%ele_type)
        j   =  max(j, npe*sec(isec)%n_ele)
    end do
    allocate(xyz (4,j), stat=err_mem)
    allocate(n2e (  j), stat=err_mem)
    allocate(vtx_ID(j), stat=err_mem)

    write(str,*),myid
    str =  './data/mesh_local_'//trim(adjustl(str))//'.dat'
    open(unit=10,file=trim(adjustl(str)))
    if(is_2d_cal) then
        write(unit=10,fmt='(A)'),'variables="CoordinateX","CoordinateY","vertex_ID","element_ID"'
    else
        write(unit=10,fmt='(A)'),'variables="CoordinateX","CoordinateY","CoordinateZ","vertex_ID","element_ID"'
    end if

    do isec=mesh(0)%sec_1,mesh(0)%sec_0
!       if(sec(isec)%is_bnd)    cycle
        ele     =  sec(isec)%ele_type
        npe     =  npe_ele_1st(ele)
        n_vtx   =  npe*sec(isec)%n_ele
        do iele=1,sec(isec)%n_ele
            call ICOPY(npe, sec(isec)%n2e(1,iele), 1, vtx_ID(1+npe*(iele-1)), 1)
        end do

        call simplify_series(n_vtx, 1, 1, vtx_ID)
        do i=1,n_vtx
            xyz(1:n_dim,i)  =  mesh(0)%xyz(1:n_dim,vtx_ID(i))
            if(vtx_ID(i) .le. mesh(0)%n_vtx) then
                xyz(4,i)=  real(mesh(0)%ID_vtx(vtx_ID(i)), dpR)
            else
                xyz(4,i)=  real(mesh(0)%ID_vtx_per(vtx_ID(i)-mesh(0)%n_vtx), dpR)
            end if
        end do

        j   =  0
        do iele=1,sec(isec)%n_ele
        do i=1,npe
            call ib_search(1, n_vtx, 1, 1, vtx_ID, sec(isec)%n2e(i,iele), ltmp, L, R)
            j       =  j+1
            n2e(j)  =  L
        end do
        end do

        if(sec(isec)%is_ghost .and. sec(isec)%is_bnd) then
            tec =  'zone T="bnd&ghost",N='
        elseif(sec(isec)%is_ghost) then
            tec =  'zone T="ghost",N='
        elseif(sec(isec)%is_bnd) then
            tec =  'zone T="bnd",N='
        else
            tec =  'zone N='
        end if
        write(str,*),n_vtx
        tec =  trim(adjustl(tec))//trim(adjustl(str))
        write(str,*),sec(isec)%n_ele
        tec =  trim(adjustl(tec))//',E='//trim(adjustl(str))
        if(ele_type_1st(ele) .eq. BAR_2) then
            str =  'FELINESEG'
        elseif(ele_type_1st(ele) .eq. QUAD_4) then
            str =  'FEQUADRILATERAL'
        elseif(ele_type_1st(ele) .eq. TRI_3) then
            str =  'FETRIANGLE'
        elseif(ele_type_1st(ele) .eq. HEXA_8) then
            str =  'FEBRICK'
        elseif(ele_type_1st(ele) .eq. PENTA_6) then
            str =  'FEBRICK'
        end if
        tec =  trim(adjustl(tec))//','//'zonetype='//trim(adjustl(str))
        tec =  trim(adjustl(tec))//',datapacking=block'
        if(n_dim .eq. 2) then
            tec =  trim(adjustl(tec))//',varlocation=([4]=cellcentered)'
        else
            tec =  trim(adjustl(tec))//',varlocation=([5]=cellcentered)'
        end if
        write(unit=10,fmt='(A)'),trim(adjustl(tec))
        do L=1,n_vtx
            write(unit=10,fmt=*),xyz(1,L)
        end do
        do L=1,n_vtx
            write(unit=10,fmt=*),xyz(2,L)
        end do
        if(n_dim .eq. 3) then
        do L=1,n_vtx
            write(unit=10,fmt=*),xyz(3,L)
        end do
        end if
        do L=1,n_vtx
            write(unit=10,fmt=*),xyz(4,L)
        end do

        if(sec(isec)%is_ghost) then
            do iele=1,sec(isec)%n_ele
                if(sec(isec)%per_path(1,iele) .gt. 0) then
                    write(unit=10,fmt=*),-sec(isec)%ID_ele_i(iele)
                else
                    write(unit=10,fmt=*), sec(isec)%ID_ele_i(iele)
                end if
            end do
        elseif(sec(isec)%is_int) then
!           do i=1,sec(isec)%n_ele
!               sec(isec)%ID_ele_i(i)   =  i
!           end do
            write(unit=10,fmt=*),sec(isec)%ID_ele_i
        else
            write(unit=10,fmt=*),(0,i=1,sec(isec)%n_ele)
        end if

        if(ele_type_1st(ele) .eq. PENTA_6) then
            do iele=1,sec(isec)%n_ele
                v(1:6)  =  n2e(1+6*(iele-1):6*iele)
                write(unit=10,fmt='(8I9)'),v(1:3),v(3),v(4:6),v(6)
            end do
        else
            do iele=1,sec(isec)%n_ele
                write(unit=10,fmt=*),n2e(1+npe*(iele-1):npe*iele)
            end do
        end if
    end do
    close(10)

    if(allocated(n2e   ))   deallocate(n2e)
    if(allocated(vtx_ID))   deallocate(vtx_ID)
    if(allocated(xyz   ))   deallocate(xyz)

    return
    end subroutine wr_local_mesh
!-------------------------------------------------------------------------------
!   setup monitors.
!-------------------------------------------------------------------------------
    subroutine set_monitor
    use var_kind_def
    use var_global
    use var_load_balance
    use var_mesh
    use var_parallel
    use var_slv
    implicit none
    logical(dpL):: is_bnd
    integer(dpI):: i,io_err,isec,iele,ID,idx(2,max_monitor)

    namelist /monitor/  monitor_idx
    if(.not. is_has_cfg)    return

    if(myid .eq. 0) then
        open(unit=10,file=trim(adjustl(cfg_file)))
        read(unit=10, nml=monitor, iostat=io_err)
        if(io_err .gt. 0)   stop 'Error: fails to read namelist:monitor.'
        close(10)
    end if
    call mpi_bcast(monitor_idx, 2*max_monitor, mpi_dpI, 0, mpi_comm_world, mpi_err)
    n_monitor_L =  0

    if(.false.) then
        do i=1,max_monitor
            if((monitor_idx(1,i) .le. 0) .or. (monitor_idx(1,i) .gt. n_seci_g)) cycle
            if((monitor_idx(2,i) .le. 0) .or. (monitor_idx(2,i) .gt. n_elei_g)) cycle
            ID  =  seci_g(6,monitor_idx(1,i))+monitor_idx(2,i)-1
            call get_isec_iele(ID, isec, iele)
            if(isec .le. 0) cycle

            n_monitor_L =  n_monitor_L+1
            monitor_idx_L(1,n_monitor_L)=  isec
            monitor_idx_L(2,n_monitor_L)=  iele
            idx        (1:2,n_monitor_L)=  monitor_idx(1:2,i)
        end do
    else
        do i=1,max_monitor
            if(monitor_idx(1,i) .le. 0) cycle
            call get_isec_iele_2(monitor_idx(1,i), monitor_idx(2,i), is_bnd, isec, iele)
            if(isec .le. 0) cycle

            n_monitor_L =  n_monitor_L+1
            monitor_idx_L(1,n_monitor_L)=  isec
            monitor_idx_L(2,n_monitor_L)=  iele
            idx        (1:2,n_monitor_L)=  monitor_idx(1:2,i)
        end do
    end if

    allocate(n_monitor_prc(0:nprc-1))
    call mpi_gather(n_monitor_L, 1, mpi_dpI, n_monitor_prc, 1, mpi_dpI, 0, &
        &  mpi_comm_world, mpi_err)
    if(myid .eq. 0) then
        n_monitor   =  n_monitor_prc(0)
        do i=1,nprc-1
            n_monitor   =  n_monitor+n_monitor_prc(i)
        end do
    end if
    i   =  2*n_monitor_L
    call i_gather(mpi_comm_world, nprc, i, idx, monitor_idx)

    return
    contains
!   ----------------------------------------------------------------------------
!   get local isec and iele.
!   ----------------------------------------------------------------------------
    subroutine get_isec_iele_2(s,e,is_bnd,isec,iele)
    use var_kind_def
    use var_cgns, only: BAR_2,TRI_3,QUAD_4,ele_type_1st
    use var_global, only: is_2d_cal
    use var_load_balance, only: n_seci_g
    use var_mesh
    implicit none
    integer(dpI),intent(in):: s,e
    logical(dpL):: is_bnd,ltmp
    integer(dpI):: isec,iele,ele_type,ID,i

    isec    =  0
    iele    =  0
    is_bnd  =  .false.
    if((s .le. 0) .or. (s .gt. n_sec_g))    return
    if((e .le. 0) .or. (e .gt. sec_g(4,n_sec_g)))   return
    ele_type=  ele_type_1st(sec_g(1,s))
    if(is_2d_cal) then
        is_bnd  =  ele_type .eq. BAR_2
    else
        is_bnd  = (ele_type .eq. TRI_3) .or. (ele_type .eq. QUAD_4)
    end if

    if(is_bnd) then
        ID  =  0
        do i=1,s
            ele_type=  ele_type_1st(sec_g(1,i))
            if(is_2d_cal) then
                ltmp=  ele_type .eq. BAR_2
            else
                ltmp= (ele_type .eq. TRI_3) .or. (ele_type .eq. QUAD_4)
            end if
            if(.not. ltmp)  cycle
            if(i .lt. s) then
                ID  =  ID+sec_g(4,i)-sec_g(3,i)+1
            else
                ID  =  ID+e
            end if
        end do
        do isec=mesh(0)%sec_1,mesh(0)%sec_0
            if((.not. sec(isec)%is_bnd) .or. (sec(isec)%ID_sec_g .ne. s))   cycle
            do iele=1,sec(isec)%n_ele
                if(sec(isec)%ID_ele_b(iele) .eq. ID)    return
            end do
        end do
    else
        do i=1,n_seci_g
            if(seci_g(5,i) .eq. s) then
                call get_isec_iele(e+seci_g(6,i)-1, isec, iele)
                return
            end if
        end do
    end if
    isec    =  0
    iele    =  0
    is_bnd  =  .false.

    return
    end subroutine get_isec_iele_2
    end subroutine set_monitor
