!-------------------------------------------------------------------------------
!   read bnd information.
!-------------------------------------------------------------------------------
    subroutine read_bnd
    use var_kind_def
    use var_air
    use var_bndv
    use var_global
    use var_parallel
    use var_slv
    use var_turb
    implicit none
    logical(dpL):: lbuf(20)
    integer(dpI):: io_err,ibuf(300),i
    real   (dpR):: rbuf(200),rtmp

    namelist /bnd_ID/   inl_ID,out_ID,sol_ID
    namelist /boundary/ L_ref,Area_ref,Ma_fs,Re_fs,t_fs,AoA_fs,tt_inl,pt_inl,vdir_inl, &
                &  pb_out,is_inl_cyl,tu_fs,BL_thick_inl,InflowSupersonic,t_wall, &
                &  flux_wall,nut_nu
    if(.not. is_has_cfg)    return

    if(myid .eq. 0) then
        open(unit=10,file=trim(adjustl(cfg_file)))

        read(unit=10,nml=boundary,iostat=io_err)
        if(io_err .gt. 0)   stop 'Error: fails to read namelist:boundary.'

        rewind(10)
        read(unit=10,nml=bnd_ID,iostat=io_err)
        if(io_err .gt. 0)   stop 'Error: fails to read namelist:bnd_ID.'

        close(10)

        lbuf(1 :10 )=  is_inl_cyl(1:10)
        ibuf(1 :300)= (/inl_ID, out_ID, sol_ID/)
        rbuf(1 :8 ) = (/L_ref, Area_ref, Ma_fs, Re_fs, t_fs, AoA_fs/)
        rbuf(9 :18) =  tt_inl(1:10)
        rbuf(19:28) =  pt_inl(1:10)
        call DCOPY(30, vdir_inl, 1, rbuf(29), 1)
        rbuf(59:68) =  pb_out(1:10)
        rbuf(69   ) =  tu_fs
        call DCOPY(50, InflowSupersonic, 1, rbuf(70), 1)
        rbuf(120:129)   =  BL_thick_inl(1:10)
        rbuf(130:139)   =  t_wall(1:10)
        rbuf(140:149)   =  flux_wall(1:10)
        rbuf(150)       =  nut_nu
    end if
    call mpi_bcast(lbuf, 10 , mpi_dpL, 0, mpi_comm_world, mpi_err)
    call mpi_bcast(ibuf, 300, mpi_dpI, 0, mpi_comm_world, mpi_err)
    call mpi_bcast(rbuf, 150, mpi_dpR, 0, mpi_comm_world, mpi_err)

    if(myid .ne. 0) then
        is_inl_cyl(1:10)=  lbuf(1 :10)
        inl_ID(1:100)   =  ibuf(1  :100)
        out_ID(1:100)   =  ibuf(101:200)
        sol_ID(1:100)   =  ibuf(201:300)
        L_ref           =  rbuf(1)
        Area_ref        =  rbuf(2)
        Ma_fs           =  rbuf(3)
        Re_fs           =  rbuf(4)
        t_fs            =  rbuf(5)
        AoA_fs(1:3)     =  rbuf(6 :8 )
        tt_inl(1:10)    =  rbuf(9 :18)
        pt_inl(1:10)    =  rbuf(19:28)
        call DCOPY(30, rbuf(29), 1, vdir_inl, 1)
        pb_out(1:10)    =  rbuf(59:68)
        tu_fs           =  rbuf(69   )
        call DCOPY(50, rbuf(70), 1, InflowSupersonic, 1)
        BL_thick_inl(1:10)  =  rbuf(120:129)
        t_wall(1:10)        =  rbuf(130:139)
        flux_wall(1:10)     =  rbuf(140:149)
        nut_nu              =  rbuf(150)
    end if

    is_inl_cyl  =  is_inl_cyl .and. (.not. is_2d_cal)
    do i=1,10
        is_inl_vis(i)   =  is_vis_cal .and. (BL_thick_inl(i) .gt. 0.0d0)
    end do
    u_fs(2:4)   =  cos(AoA_fs(1:3)*pi/1.8d2)
    call norm_vec(3, u_fs(2), rtmp)
    a_fs        =  sqrt(gk*rr*t_fs)
    ua_fs       =  Ma_fs*a_fs
    u_fs(2:4)   =  u_fs(2:4)*ua_fs
    mu_fs       =  1.461d-6*sqrt(t_fs**3)/(t_fs+t_Suth)
    u_fs(1)     =  Re_fs*mu_fs/(ua_fs*L_ref)
    u_fs(5)     =  u_fs(1)*rr*t_fs

    if(is_flow_int) then
        tref    =  pref/(rr*rref)
        mu_ref  =  1.461d-6*sqrt(tref**3)/(tref+t_suth)
        Re_ref  =  rref*uref*L_ref/mu_ref
        call nut2nu(nut_nu*mu_ref/rref, mu_ref/rref, nut_ref)
    else
        rref=  u_fs(1)
        uref=  ua_fs
        pref=  u_fs(5)
        dir_drag(1:3)   =  u_fs(2:4)
        call norm_vec(3, dir_drag, rtmp)
        rtmp=  dir_lift(1)*dir_drag(1)+dir_lift(2)*dir_drag(2)+dir_lift(3)*dir_drag(3)
        dir_lift=  dir_lift-rtmp*dir_drag
        call norm_vec(3, dir_lift, rtmp)
        if(RANS_model .eq. RANS_SA) then
            call nut2nu(nut_nu*mu_fs/rref, mu_fs/rref, nut_ref)
        elseif(RANS_model .eq. RANS_WA) then
            nut_ref =  4.0d0*mu_fs/rref
        end if
    end if

    tm_ref  =  L_ref/uref
    if(is_2d_cal) then
        res_ref(1)  =  tm_ref/(rref*L_ref**3)
    else
        res_ref(1)  =  tm_ref/(rref*L_ref**3)
    end if
    res_ref(2:4)=  res_ref(1)/uref
    res_ref(5  )=  res_ref(2)/uref

    return
    end subroutine read_bnd
!-------------------------------------------------------------------------------
!   setup bnd ID.
!-------------------------------------------------------------------------------
    subroutine set_bnd_ID(lev)
    use var_kind_def
    use var_bndv
    use var_cgns
    use var_load_balance, only: n_bocos,bocoinfo
    use var_mesh
    use var_slv, only: is_vis_cal
    implicit none
    integer(dpI),intent(in):: lev
    logical(dpL):: ltmp
    integer(dpI):: isec,bct,ID,ID_group,v(1000),n_in,n_ot,n_wall,i
    real   (dpR):: rtmp

    if(n_bocos .le. 0)  return
    v       =  0
    n_in    =  0
    n_ot    =  0
    n_wall  =  0
    do i=1,n_bocos
        bct =  bocoinfo(1,i)
        if(bct .eq. BCInflow) then
            n_in    =  n_in+1
            v(i)    =  n_in
        elseif(bct .eq. BCInflowSupersonic) then
            n_in    =  n_in+1
            v(i)    =  n_in
        elseif(bct .eq. BCOutflow) then
            n_ot    =  n_ot+1
            v(i)    =  n_ot
        elseif(bct .eq. BCOutflowSupersonic) then
            n_ot    =  n_ot+1
            v(i)    =  n_ot
        elseif(bct .eq. BCWallViscous) then
            n_wall  =  n_wall+1
            v(i)    =  n_wall
        elseif(bct .eq. BCWallInviscid) then
            n_wall  =  n_wall+1
            v(i)    =  n_wall
        else
!           do nothing.
        end if
    end do

    n_inl   =  0
    n_out   =  0
    n_sol   =  0
    do i=1,n_in
        n_inl   =  max(n_inl, inl_ID(i))
    end do
    do i=1,n_ot
        n_out   =  max(n_out, out_ID(i))
    end do
    do i=1,n_wall
        n_sol   =  max(n_sol, sol_ID(i))
    end do

    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(.not. sec(isec)%is_bnd)  cycle

        bct =  sec(isec)%bct
        ID  =  v(sec(isec)%ID_bnd)
        if(bct .eq. BCInflow) then
            if(ID .le. 0)   stop 'Error: set bnd_ID fails.'
            ID_group=  inl_ID(ID)
            ltmp= (tt_inl(ID_group) .le. 0.0d0) .or. (pt_inl(ID_group) .le. 0.0d0) .or. &
                & (vdir_inl(1,ID_group) .le. -1.0d3)
            if(ltmp)    stop 'Error: inflow boundary condition not defined.'
            call norm_vec(3, vdir_inl(1,ID_group), rtmp)
        elseif(bct .eq. BCInflowSupersonic) then
            if(ID .le. 0)   stop 'Error: set bnd_ID fails.'
            ID_group=  inl_ID(ID)
            ltmp=  InflowSupersonic(1,ID) .le. 0.0d0
            if(ltmp)    stop 'Error: InflowSupersonic boundary condition not defined.'
        elseif(bct .eq. BCOutflow) then
            if(ID .le. 0)   stop 'Error: set bnd_ID fails.'
            ID_group=  out_ID(ID)
            ltmp=  pb_out(ID_group) .le. 0.0d0
            if(ltmp)    stop 'Error: outflow boundary condition not defined.'
        elseif(bct .eq. BCOutflowSupersonic) then
            if(ID .le. 0)   stop 'Error: set bnd_ID fails.'
            ID_group=  out_ID(ID)
        elseif(bct .eq. BCWallViscous) then
            if(ID .le. 0)   stop 'Error: set bnd_ID fails.'
            ID_group=  sol_ID(ID)
            sec(isec)%is_t      =  is_vis_cal .and. (t_wall(ID_group) .gt. 0.0d0)
            sec(isec)%is_flux   =  is_vis_cal .and. (flux_wall(ID_group) .gt. -1.0d-20)
        elseif(bct .eq. BCWallInviscid) then
            if(ID .le. 0)   stop 'Error: set bnd_ID fails.'
            ID_group=  sol_ID(ID)
        else
            ID_group=  0
        end if
        sec(isec)%ID_group  =  ID_group
    end do

    return
    end subroutine set_bnd_ID
!-------------------------------------------------------------------------------
!   Inflow boundary.
!-------------------------------------------------------------------------------
    subroutine bnd_inflow(n,v_dir,tt_inl,pt_inl,ue,ub)
    use var_kind_def
    use var_air
    use var_global
    implicit none
    real   (dpR),intent(in):: n(*),v_dir(*),tt_inl,pt_inl,ue(*)
    logical(dpL):: lava(2),ltmp
    real   (dpR):: ub(*),dn,Jn,H,coe1,coe2,coe3,rt(2),t,p,rh,a,ke,uaa

    dn  =  n(1)*v_dir(1)+n(2)*v_dir(2)+n(3)*v_dir(3)
    Jn  =  ue(2)*n(1)+ue(3)*n(2)+ue(4)*n(3)+2.0d0*dsqrt(gk*ue(5)/ue(1))/gk1
    H   =  gk*rr*tt_inl
    coe1=  H*dn*dn-0.5d0*gk1*Jn*Jn
    coe2=  4.0d0*H*dn/gk1
    coe3=  4.0d0*H/(gk1*gk1)-Jn*Jn
    lava=  .false.
    if(coe2*coe2-4.0d0*coe1*coe3 .gt. 0.0d0) then
        rt(1)   = (-coe2+dsqrt(coe2*coe2-4.0d0*coe1*coe3))*0.5d0/coe1
        rt(2)   = (-coe2-dsqrt(coe2*coe2-4.0d0*coe1*coe3))*0.5d0/coe1
        if((rt(1) .gt. 0.0d0) .and. (rt(1) .lt. 1.0d0)) lava(1) =  .true.
        if((rt(2) .gt. 0.0d0) .and. (rt(2) .lt. 1.0d0)) lava(2) =  .true.
    end if
    ltmp=  lava(1) .or. lava(2)
    if(lava(1) .and. lava(2))   ltmp=  .false.
    if(lava(1)) uaa =  rt(1)
    if(lava(2)) uaa =  rt(2)

    if(ltmp) then
        t       =  tt_inl/(1.0d0+0.5d0*gk1*uaa*uaa)
        p       =  pt_inl*(tt_inl/t)**(-gk/gk1)
        rh      =  p/(rr*t)
        rt(1)   =  uaa*sqrt(gk*p/rh)
        ub(1)   =  rh
        ub(2)   =  rt(1)*v_dir(1)
        ub(3)   =  rt(1)*v_dir(2)
        ub(4)   =  rt(1)*v_dir(3)
        ub(5)   =  p
    else
        a       =  sqrt(ue(2)*ue(2)+ue(3)*ue(3)+ue(4)*ue(4))
        ub(2)   =  a*v_dir(1)
        ub(3)   =  a*v_dir(2)
        ub(4)   =  a*v_dir(3)
        ke      =  0.5d0*(ub(2)**2+ub(3)**2+ub(4)**2)
        t       =  tt_inl-ke/cp
        p       =  pt_inl*(tt_inl/t)**(-gk/gk1)
        ub(1)   =  p/(rr*t)
        ub(5)   =  p
    end if

    return
    end subroutine bnd_inflow
!-------------------------------------------------------------------------------
!   Inflow boundary, viscous.
!-------------------------------------------------------------------------------
    subroutine bnd_inflow_vis(v_dir,tt_inl,pt_inl,ue,dnw,BL_thick,ub)
    use var_kind_def
    use var_air
    implicit none
    real(dpR),intent(in):: v_dir(*),tt_inl,pt_inl,ue(*),dnw,BL_thick
    real(dpR):: ub(*),p,t,coe2

    p       =  ue(5)
    t       = (pt_inl/p)**(-gk1/gk)*tt_inl
    coe2    =  dsqrt(2.0d0*cp*dabs(tt_inl-t))*(dnw/BL_thick)**(1.0d0/7.0d0)
    ub(2)   =  coe2*v_dir(1)
    ub(3)   =  coe2*v_dir(2)
    ub(4)   =  coe2*v_dir(3)
!   t       =  Tt-coe2*coe2/cp*0.5d0
    ub(1)   =  p/(rr*t)
    ub(5)   =  p

    return
    end subroutine bnd_inflow_vis
!-------------------------------------------------------------------------------
!   Outflow boundary.
!-------------------------------------------------------------------------------
    subroutine bnd_outflow(n,ue,pb_out,ub)
    use var_kind_def
    use var_air, only: gk,gk1
    implicit none
    real(kind=8):: n(*),ue(*),pb_out,ub(*),a,b,rh,un

    a       =  ue(5)/ue(1)**gk
    rh      = (pb_out/a)**(1.0d0/gk)
    b       =  ue(2)*n(1)+ue(3)*n(2)+ue(4)*n(3)
    a       =  b+2.0d0*dsqrt(gk*ue(5)/ue(1))/gk1
    un      =  a-2.0d0*dsqrt(gk*pb_out/rh)/gk1
    if(un .lt. 0.0d0)   un  =  0.0d0

    ub(1)   =  rh
    ub(2)   =  ue(2)+(un-b)*n(1)
    ub(3)   =  ue(3)+(un-b)*n(2)
    ub(4)   =  ue(4)+(un-b)*n(3)
    ub(5)   =  pb_out

    return
    end subroutine bnd_outflow
!-------------------------------------------------------------------------------
!   Farfield boundary.
!-------------------------------------------------------------------------------
    subroutine bnd_ffd(n,ue,ub)
    use var_kind_def
    use var_air
    use var_bndv
    use var_global
    implicit none
    real(dpR):: n(*),ue(*),ub(*),un_s,un_e,coe1,coe2,un,a,t

    un_s=  u_fs(2)*n(1)+u_fs(3)*n(2)+u_fs(4)*n(3)
    un_e=  ue  (2)*n(1)+ue  (3)*n(2)+ue  (4)*n(3)

    coe1=  un_s-2.0d0*a_fs/gk1
    coe2=  un_e+2.0d0*sqrt(gk*ue(5)/ue(1))/gk1
    un  =  0.5d0*(coe1+coe2)
    a   =  0.25d0*gk1*(coe2-coe1)
    if(un .le. 0.0d0) then
        ub(1)   = (a*a*u_fs(1)**gk/(gk*u_fs(5)))**(1.0d0/gk1)
        t       =  un-un_s
        ub(2)   =  u_fs(2)+t*n(1)
        ub(3)   =  u_fs(3)+t*n(2)
        ub(4)   =  u_fs(4)+t*n(3)
    else
        ub(1)   = (a*a*ue(1)**gk/(gk*ue(5)))**(1.0d0/gk1)
        t       =  un-un_e
        ub(2)   =  ue(2)+t*n(1)
        ub(3)   =  ue(3)+t*n(2)
        ub(4)   =  ue(4)+t*n(3)
    end if
    ub(5)   =  ub(1)*a*a/gk

    return
    end subroutine bnd_ffd
!-------------------------------------------------------------------------------
!   read in the boundary profiles.
!-------------------------------------------------------------------------------
    subroutine read_bnd_profile
    use var_kind_def
    use var_bndv
    use var_cgns
    use var_global, only: is_has_cfg,cfg_file
    use var_mesh
    use var_parallel
    implicit none
    character(len=80):: str,s
    integer(dpI):: ID,npnt,idx,io_err,i,bct,isec
    real   (dpR):: rbf1(10000),rbf2(10000)

    if(.not. is_has_cfg)    return

!   ----------------------------------------------------------------------------
!   read profile.
    idx =  0
    if(myid .eq. 0) then
        open(unit=10,file=cfg_file)
        do while(.true.)
            read(unit=10,fmt=*,iostat=io_err),str
            if(io_err .ne. 0)   exit

            s   =  str
            call upper_string(s)
            if(trim(s) .eq. '&INL_PROFILE') then
                read(unit=10,fmt='(A)'),str
                call upper_string(str)
                if(index(str, 'ID'  ) .gt. 0)   call get_ID_string(str, ID)

                read(unit=10,fmt='(A)'),str
                call upper_string(str)
                if(index(str, 'NPNT') .gt. 0)   call get_ID_string(str, npnt)

                rbf1(idx+1) =  real(BCInflow, dpR)
                rbf1(idx+2) =  real(ID, dpR)
                rbf1(idx+3) =  real(npnt, dpR)
                idx =  idx+3
                do i=1,npnt
                    read(unit=10,fmt=*),rbf1(idx+1:idx+6)
                    idx =  idx+6
                end do
            elseif(trim(s) .eq. '&OUT_PROFILE') then
                read(unit=10,fmt='(A)'),str
                call upper_string(str)
                if(index(str, 'ID'  ) .gt. 0)   call get_ID_string(str, ID)

                read(unit=10,fmt='(A)'),str
                call upper_string(str)
                if(index(str, 'NPNT') .gt. 0)   call get_ID_string(str, npnt)

                rbf1(idx+1) =  real(BCOutflow, dpR)
                rbf1(idx+2) =  real(ID, dpR)
                rbf1(idx+3) =  real(npnt, dpR)
                idx =  idx+3
                do i=1,npnt
                    read(unit=10,fmt=*),rbf1(idx+1:idx+2)
                    idx =  idx+2
                end do
            end if
        end do
        close(10)
    end if
!   read profile.
!   ----------------------------------------------------------------------------

    call mpi_bcast(idx, 1, mpi_dpI, 0, mpi_comm_world, mpi_err)
    if(idx .le. 0)  return
    call mpi_bcast(rbf1, idx, mpi_dpR, 0, mpi_comm_world, mpi_err)

    io_err  =  idx
    idx     =  0
    do while(idx .lt. io_err)
        bct =  nint(rbf1(idx+1), dpI)
        ID  =  nint(rbf1(idx+2), dpI)
        npnt=  nint(rbf1(idx+3), dpI)
        idx =  idx+3

        if(bct .eq. BCInflow) then
            rbf2(1:6*npnt)  =  rbf1(idx+1:idx+6*npnt)
            idx =  idx+6*npnt
        elseif(bct .eq. BCOutflow) then
            rbf2(1:2*npnt)  =  rbf1(idx+1:idx+2*npnt)
            idx =  idx+2*npnt
        end if

        do isec=mesh(0)%sec_1,mesh(0)%sec_0
            if(.not. sec(isec)%is_bnd)  cycle
            if((sec(isec)%bct .ne. bct) .or. (sec(isec)%ID_group .ne. ID))  cycle

            sec(isec)%is_profiled   =  .true.
            sec(isec)%n_vtx_profile =  npnt
            if(bct .eq. BCInflow) then
                allocate(sec(isec)%bc_profile(6,npnt))
                call DCOPY(6*npnt, rbf2, 1, sec(isec)%bc_profile, 1)
            elseif(bct .eq. BCOutflow) then
                allocate(sec(isec)%bc_profile(2,npnt))
                call DCOPY(2*npnt, rbf2, 1, sec(isec)%bc_profile, 1)
            end if
        end do
    end do

    return
    end subroutine read_bnd_profile
!-------------------------------------------------------------------------------
!   get boundary ID from string.
!-------------------------------------------------------------------------------
    subroutine get_ID_string(s,ID)
    use var_kind_def
    implicit none
    character(len=*),intent(in):: s
    character(len=80):: str
    integer(dpI):: ID,i

    do i=1,len_trim(s)
        if(s(i:i) .eq. '=') then
            str =  trim(adjustl(s(i+1:len_trim(s))))
            read(str,*),ID
            return
        end if
    end do
    stop 'Error: failed to get ID from string.'

    return
    end subroutine get_ID_string
!-------------------------------------------------------------------------------
!   interpolate the boundary profile, 1D.
!-------------------------------------------------------------------------------
    subroutine get_bnd_profile_1d(ibc,x,y)
    use var_kind_def
    use var_bndv
    use var_cgns
    use var_mesh
    implicit none
    integer(dpI):: ibc,i,N,LDA
    real   (dpR):: x,y(*),x1,x0,c

    if(sec(ibc)%bct .eq. BCInflow) then
        LDA =  5
    elseif(sec(ibc)%bct .eq. BCOutflow) then
        LDA =  1
    else
        stop 'Error: wrong bc type for profile definition.'
    end if

    N   =  sec(ibc)%n_vtx_profile
    if(x .ge. sec(ibc)%bc_profile(1,N)) then
        y(1:LDA)=   sec(ibc)%bc_profile(2:LDA+1,N)
    elseif(x .le. sec(ibc)%bc_profile(1,1)) then
        y(1:LDA)=   sec(ibc)%bc_profile(2:LDA+1,1)
    else
        do i=1,N-1
            x1  =  sec(ibc)%bc_profile(1,i  )
            x0  =  sec(ibc)%bc_profile(1,i+1)
            if((x .ge. x1) .and. (x .le. x0)) then
                c   = (x-x1)/(x0-x1)
                y(1:LDA)= (1.0d0-c)*sec(ibc)%bc_profile(2:LDA+1,i  ) &
                        &       +c *sec(ibc)%bc_profile(2:LDA+1,i+1)
                return
            end if
        end do
    end if

    return
    end subroutine get_bnd_profile_1d
!-------------------------------------------------------------------------------
!   cal nu for given nut, SA.
!-------------------------------------------------------------------------------
    subroutine nut2nu(nut,nul,nu)
    use var_kind_def
    implicit none
    integer(dpI):: iter
    real   (dpR):: nut,nul,nu,f,fx,df

    nu  =  nut
    do iter=1,10
        f   =  nu**4-nut*nu**3-(7.1d0*nul)**3*nut
        fx  =  4.0d0*nu**3-3.0d0*nut*nu*nu
        df  = -f/fx
        nu  =  nu+df
        if(abs(df) .ge. 1.0d2*nut) then
            nu  =  nut
            return
        end if
        if(abs(df) .le. 1.0d-6*nut) return
    end do

    return
    end subroutine nut2nu
