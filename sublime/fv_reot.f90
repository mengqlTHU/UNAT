!-------------------------------------------------------------------------------
!   set the REOT boundary condition.
!-------------------------------------------------------------------------------
    subroutine fv_set_reot
    use var_kind_def
    use var_avg
    use var_cgns, only: BCOutflow
    use var_fv
    use var_global, only: is_has_cfg,cfg_file,err_mem
    use var_mesh
    use var_parallel
    implicit none
    logical(dpL):: ltmp
    integer(dpI):: isec,iele,ID,io_err,n_vtx,L,R,i,ipe,idx,iavg,ifac
    real   (dpR):: r0_reot(10),p0_reot(10),rbuf(20)
    integer(dpI),allocatable:: n2e(:),vtx(:),IDv(:)
    real   (dpR),allocatable:: xyz(:,:)

    namelist /avg_def/  r0_reot,p0_reot
    if(.not. is_has_cfg)    return

    r0_reot = -1.0d10
    p0_reot = -1.0d10
    if(myid .eq. 0) then
        open(unit=10,file=trim(adjustl(cfg_file)))
        read(unit=10, nml=avg_def, iostat=io_err)
        if(io_err .gt. 0)   stop 'Error: fails to read namelist:avg_def.'
        close(10)

        rbuf(1 :10) =  r0_reot(1:10)
        rbuf(11:20) =  p0_reot(1:10)
    end if
    call mpi_bcast(rbuf, 20, mpi_dpR, 0, mpi_comm_world, mpi_err)
    if(myid .ne. 0) then
        r0_reot(1:10)   =  rbuf(1 :10)
        p0_reot(1:10)   =  rbuf(11:20)
    end if

    idx =  0
    do isec=mesh(0)%sec_1,mesh(0)%sec_0
        if(sec(isec)%bct .ne. BCOutflow)    cycle
        ID  =  sec(isec)%ID_group
        if((r0_reot(ID) .le. 0.0d0) .or. (p0_reot(ID) .le. 0.0d0))  cycle
        fv(isec)%is_reot=  .true.
        idx             =  max(idx, sec(isec)%npe*sec(isec)%n_ele)
    end do
    if(idx .gt. 0) then
        allocate(n2e(  idx), stat=err_mem)
        allocate(vtx(  idx), stat=err_mem)
        allocate(IDv(  idx), stat=err_mem)
        allocate(xyz(3,idx), stat=err_mem)
    end if

    do isec=mesh(0)%sec_1,mesh(0)%sec_0
        if(.not. fv(isec)%is_reot)  cycle

        n_vtx   =  sec(isec)%npe*sec(isec)%n_ele
        call ICOPY(n_vtx, sec(isec)%n2e, 1, vtx, 1)
        call simplify_series(n_vtx, 1, 1, vtx)
        do i=1,n_vtx
            xyz(1:3,i)  =  mesh(0)%xyz(1:3,vtx(i))
            IDv(    i)  =  mesh(0)%ID_vtx( vtx(i))
        end do

        idx =  1
        do iele=1,sec(isec)%n_ele
        do ipe =1,sec(isec)%npe
            call ib_search(1, n_vtx, 1, 1, vtx, sec(isec)%n2e(ipe,iele), ltmp, L, R)
            if((.not. ltmp) .or. (L .ne. R))    stop 'Error: fv_set_reot fails.'
            n2e(idx)=  L
            idx     =  idx+1
        end do
        end do
        call avg_add(0, sec(isec)%ID_group, BCoutflow, iavg)

        ID  =  sec(isec)%ID_group
        avg(iavg)%r_reot=  r0_reot(ID)
        avg(iavg)%p_reot=  p0_reot(ID)

        call avg_add_patch(avg(iavg), .true., myid, isec, sec(isec)%ele_type, &
            &  sec(isec)%n_ele, n2e, n_vtx, IDv, xyz, ifac)

!       record the mortar information.
        allocate(avg(iavg)%facL(ifac)%mortar(sec(isec)%n_ele), stat=err_mem)
        avg(iavg)%facL(ifac)%mortar =  0
        do i=1,mesh(0)%n_mortar_b
            if(mesh(0)%mortar_LR(3,i) .ne. isec)    cycle
            iele=  mesh(0)%mortar_LR(4,i)
            avg(iavg)%facL(ifac)%mortar(iele)   =  i
        end do
        do iele=1,sec(isec)%n_ele
            if(avg(iavg)%facL(ifac)%mortar(iele) .le. 0)    stop 'Error: fv_reot fails.'
        end do
    end do

    call avg_synchronize

    if(allocated(n2e))  deallocate(n2e)
    if(allocated(vtx))  deallocate(vtx)
    if(allocated(IDv))  deallocate(IDv)
    if(allocated(xyz))  deallocate(xyz)

    return
    end subroutine fv_set_reot
!-------------------------------------------------------------------------------
!   fv REOT boundary condition.
!-------------------------------------------------------------------------------
    subroutine fv_reot
    use var_kind_def
    use var_avg
    use var_fv
    use var_global, only: z=>rotation_axis
    use var_mesh
    use var_parallel
    implicit none
    logical(dpL):: ltmp
    integer(dpI):: iavg,ifac,iele,im,sL,eL,sR,eR,i,j,M,idx,lev,ID,bct,idmL,idmR,prc
    real   (dpR):: n(4),u(5),r,s1,s2,r1,r2,p0,r0,a,b,c,dr

    if(myid_avg .lt. 0) return

!   ----------------------------------------------------------------------------
!   pack and allgather data.
    idx =  0
    do iavg=1,n_avg
!       lev, ID, bct, idmL, idmR, myid
        abf1(idx+1) =  real(avg(iavg)%lev , dpR)
        abf1(idx+2) =  real(avg(iavg)%ID  , dpR)
        abf1(idx+3) =  real(avg(iavg)%bct , dpR)
        abf1(idx+4) =  real(avg(iavg)%idmL, dpR)
        abf1(idx+5) =  real(avg(iavg)%idmR, dpR)
        abf1(idx+6) =  real(myid          , dpR)
        idx         =  idx+6

        avg(iavg)%u =  0.0d0

        do ifac=1,avg(iavg)%nfacL
            if(avg(iavg)%facL(ifac)%myid .ne. myid) cycle
            do iele=1,avg(iavg)%facL(ifac)%n_ele
                im      =  avg(iavg)%facL(ifac)%mortar(iele)
                i       =  avg(iavg)%facL(ifac)%ID_r  (iele)
                sL      =  mesh(0)%mortar_LR  (1  ,im)
                eL      =  mesh(0)%mortar_LR  (2  ,im)
                sR      =  mesh(0)%mortar_LR  (3  ,im)
                eR      =  mesh(0)%mortar_LR  (4  ,im)
                n(1:4)  =  mesh(0)%mortar_n_vg(1:4,im)

                u(1:3)  =  sec(sR)%cen(1:3,eR)
                a       =  z(1)*u(1)+z(2)*u(2)+z(3)*u(3)
                u       =  u-a*u
                call norm_vec(3, u, a)
                call crs_prd(z, u, n)

                u(1:5)  =  fv(sL)%u(1:5,eL)
                avg(iavg)%u(1,i)=  avg(iavg)%u(1,i)+n(4)*u(1)
                avg(iavg)%u(2,i)=  avg(iavg)%u(2,i)+n(4)*(n(1)*u(2)+n(2)*u(3)+n(3)*u(4))
                avg(iavg)%u(3,i)=  avg(iavg)%u(3,i)+n(4)*u(5)
            end do
        end do
        call DCOPY(3*(avg(iavg)%idmL-1), avg(iavg)%u, 1, abf1(idx+1), 1)
        idx =  idx+3*(avg(iavg)%idmL-1)
    end do
    M   =  idx
    call r_allgather(mpi_comm_avg, nprc_avg, M, abf1, abf2)
!   pack and allgather data.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   update the REOT.
    idx =  0
    do while(idx .lt. M)
        lev =  nint(abf2(idx+1), dpI)
        ID  =  nint(abf2(idx+2), dpI)
        bct =  nint(abf2(idx+3), dpI)
        idmL=  nint(abf2(idx+4), dpI)
        idmR=  nint(abf2(idx+5), dpI)
        prc =  nint(abf2(idx+6), dpI)
        idx =  idx+6

        ltmp=  .false.
        do iavg=1,n_avg
            ltmp= (lev .eq. avg(iavg)%lev) .and. (ID .eq. avg(iavg)%ID) .and. &
                & (bct .eq. avg(iavg)%bct) .and. (prc .ne. myid)
            if(ltmp)    exit
        end do

        if(.not. ltmp) then
            idx =  idx+3*(idmL-1)
            cycle
        end if

        call DAXPY(3*(idmL-1), 1.0d0, abf2(idx+1), 1, avg(iavg)%u, 1)
    end do
!   update the REOT.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   solving radial equilibrium equation.
    do iavg=1,n_avg
        do i=1,avg(iavg)%idmL-1
            u          (1:3  )  =  avg(iavg)%u(1:3,i)/avg(iavg)%A_L(i)
            avg(iavg)%u(1  ,i)  =  u(1)*u(2)**2/avg(iavg)%r_L(i)
        end do

!       ------------------------------------------------------------------------
!       solving radial equilibrium equation.
        r   =  avg(iavg)%r_reot
        do j=1,avg(iavg)%idmL-2
            if(avg(iavg)%r_L(j) .le. r .and. avg(iavg)%r_L(j+1) .ge. r) then
                if(.true.) then
                    s1  =  avg(iavg)%u(1,j  )
                    s2  =  avg(iavg)%u(1,j+1)
                    r1  =  avg(iavg)%r_L(j  )
                    r2  =  avg(iavg)%r_L(j+1)
                    p0  =  avg(iavg)%p_reot
                    r0  =  r
                    a   =  0.5d0*(s1-s2)/(r1-r2)
                    b   =  s1-2.0d0*a*r1
                    c   =  p0-a*r0*r0-b*r0
                    avg(iavg)%u(3,j  )  =  a*r1*r1+b*r1+c
                    avg(iavg)%u(3,j+1)  =  a*r2*r2+b*r2+c
                else
                    dr                  =  r-avg(iavg)%r_L(j  )
                    avg(iavg)%u(3,j  )  =  avg(iavg)%p_reot-dr*avg(iavg)%u(1,j  )
                    dr                  =  r-avg(iavg)%r_L(j+1)
                    avg(iavg)%u(3,j+1)  =  avg(iavg)%p_reot-dr*avg(iavg)%u(1,j+1)
                end if

                exit
            end if
        end do
        do i=j-1,1,-1
            avg(iavg)%u(3,i)=  avg(iavg)%u(3,i+1) &
                & -0.5d0*(avg(iavg)%u(1,i)+avg(iavg)%u(1,i+1)) &
                & *(avg(iavg)%r_L(i+1)-avg(iavg)%r_L(i))
        end do
        do i=j+2,avg(iavg)%idmL-1
            avg(iavg)%u(3,i)=  avg(iavg)%u(3,i-1) &
                & +0.5d0*(avg(iavg)%u(1,i-1)+avg(iavg)%u(1,i)) &
                & *(avg(iavg)%r_L(i)-avg(iavg)%r_L(i-1))
        end do
!       solving radial equilibrium equation.
!       ------------------------------------------------------------------------

        do ifac=1,avg(iavg)%nfacL
            if(avg(iavg)%facL(ifac)%myid .ne. myid) cycle
            do iele=1,avg(iavg)%facL(ifac)%n_ele
                im      =  avg(iavg)%facL(ifac)%mortar(iele)
                i       =  avg(iavg)%facL(ifac)%ID_r  (iele)
                sR      =  mesh(0)%mortar_LR(3,im)
                eR      =  mesh(0)%mortar_LR(4,im)
                fv(sR)%bc(1,eR) =  avg(iavg)%u(3,i)
            end do
        end do
    end do
!   solving radial equilibrium equation.
!   ----------------------------------------------------------------------------

    return
    end subroutine fv_reot

!-------------------------------------------------------------------------------
!   fv REOT boundary condition.
!-------------------------------------------------------------------------------
    subroutine fv_reot_sw
    use var_kind_def
    use var_avg
    use var_fv
    use var_global, only: z=>rotation_axis
    use var_mesh
    use var_parallel
    implicit none
    logical(dpL):: ltmp
    integer(dpI):: iavg,ifac,iele,im,sL,eL,sR,eR,i,j,M,idx,lev,ID,bct,idmL,idmR,prc
    real   (dpR):: n(4),u(5),r,s1,s2,r1,r2,p0,r0,a,b,c,dr

    if(myid_avg .lt. 0) return
    
    stop 'fv_reot has not been implemented on Sunway platform'

    end subroutine fv_reot_sw