!-------------------------------------------------------------------------------
!   bnd.
!-------------------------------------------------------------------------------
    subroutine fv_bnd(lev)
    use var_kind_def
    use var_bndv, only: is_inl_vis,InflowSupersonic,BL_thick_inl
    use var_cgns
    use var_fv, only: fv
    use var_mesh
    implicit none
    integer(dpI),intent(in):: lev
    integer(dpI):: im,sL,eL,sR,eR,bct,ID
    real   (dpR):: un,vg(3),n(3),uL(5),um(5)

    call fv_reot
    do im=1,mesh(lev)%n_mortar_b
        sL      =  mesh(lev)%mortar_LR(1,im)
        eL      =  mesh(lev)%mortar_LR(2,im)
        sR      =  mesh(lev)%mortar_LR(3,im)
        eR      =  mesh(lev)%mortar_LR(4,im)
        bct     =  sec(sR)%bct
        ID      =  sec(sR)%ID_group
        vg(1:3) =  sec(sR)%vg(1:3,eR)
        n (1:3) =  mesh(lev)%mortar_n_vg(1:3,im)
        uL(1:5) =  fv  (sR )%uL         (1:5,eR)

        if(bct .eq. BCFarfield) then
            call bnd_ffd(n, uL, um)
        elseif(bct .eq. BCWallInviscid) then
            un      = (vg(1)-uL(2))*n(1)+(vg(2)-uL(3))*n(2)+(vg(3)-uL(4))*n(3)
            um(1)   =  uL(1)
            um(2:4) =  uL(2:4)+un*n(1:3)
            um(5)   =  uL(5)
        elseif(bct .eq. BCWallViscous) then
            um(1)   =  uL(1)
            um(2:4) =  vg(1:3)
            um(5)   =  uL(5)
        elseif(bct .eq. BCInflow) then
            if(.not. is_inl_vis(ID)) then
                call bnd_inflow(n,fv(sR)%bc(3,eR),fv(sR)%bc(1,eR),fv(sR)%bc(2,eR),uL,um)
            elseif(sec(sL)%dnw(eL) .gt. BL_thick_inl(ID)) then
                call bnd_inflow(n,fv(sR)%bc(3,eR),fv(sR)%bc(1,eR),fv(sR)%bc(2,eR),uL,um)
            else
                call bnd_inflow_vis(fv(sR)%bc(3,eR),fv(sR)%bc(1,eR),fv(sR)%bc(2,eR), &
                    &  uL,sec(sL)%dnw(eL),BL_thick_inl(ID),um)
            end if
        elseif(bct .eq. BCInflowSupersonic) then
            um(1:5) =  InflowSupersonic(1:5,ID)
        elseif(bct .eq. BCOutflow) then
            call bnd_outflow(n, uL, fv(sR)%bc(1,eR), um)
        elseif(bct .eq. BCOutflowSupersonic) then
            um(1:5) =  uL(1:5)
        else
            stop 'Error: boundary type not supported.'
        end if
        fv(sR)%u (1:5,eR)   =  um(1:5)
        fv(sR)%uR(1:5,eR)   =  2.0d0*um(1:5)-uL(1:5)
    end do

    return
    end subroutine fv_bnd
!-------------------------------------------------------------------------------
!   bnd. Sunway version.
!-------------------------------------------------------------------------------
    subroutine fv_bnd_sw(lev)
    use var_kind_def
    use var_bndv, only: is_inl_vis,InflowSupersonic,BL_thick_inl
    use var_cgns
    use var_fv, only: fv
    use var_fv_array
    use var_mesh
    implicit none
    integer(dpI),intent(in):: lev
    integer(dpI):: im,sL,eL,sR,eR,bct,ID,idxL,idxR
    real   (dpR):: un,vg(3),n(3),uL(5),um(5)

    call fv_reot_sw

! call fv_struct_to_array(lev)

    do im=1,mesh(lev)%n_mortar_b
        sL      =  mesh(lev)%mortar_LR(1,im)
        eL      =  mesh(lev)%mortar_LR(2,im)
        sR      =  mesh(lev)%mortar_LR(3,im)
        eR      =  mesh(lev)%mortar_LR(4,im)
        idxL    =  sec(sL)%ID_ele_g(eL)
        idxR    =  sec(sR)%ID_ele_g(eR)
        bct     =  sec(sR)%bct
        ID      =  sec(sR)%ID_group
        vg(1:3) =  sec(sR)%vg(1:3,eR)
        n (1:3) =  mesh(lev)%mortar_n_vg(1:3,im)
        uL(1:5) =  fv_uL(1:5,idxR)

        if(bct .eq. BCFarfield) then
            call bnd_ffd(n, uL, um)
        elseif(bct .eq. BCWallInviscid) then
            un      = (vg(1)-uL(2))*n(1)+(vg(2)-uL(3))*n(2)+(vg(3)-uL(4))*n(3)
            um(1)   =  uL(1)
            um(2:4) =  uL(2:4)+un*n(1:3)
            um(5)   =  uL(5)
        elseif(bct .eq. BCWallViscous) then
            um(1)   =  uL(1)
            um(2:4) =  vg(1:3)
            um(5)   =  uL(5)
        elseif(bct .eq. BCInflow) then
            write(*,*),'BCInflow boundary condition has not been verified'
            ! if(.not. is_inl_vis(ID)) then
                ! call bnd_inflow(n,fv(sR)%bc(3,eR),fv(sR)%bc(1,eR),fv(sR)%bc(2,eR),uL,um)
            ! elseif(sec(sL)%dnw(eL) .gt. BL_thick_inl(ID)) then
                ! call bnd_inflow(n,fv(sR)%bc(3,eR),fv(sR)%bc(1,eR),fv(sR)%bc(2,eR),uL,um)
            ! else
                ! call bnd_inflow_vis(fv(sR)%bc(3,eR),fv(sR)%bc(1,eR),fv(sR)%bc(2,eR), &
                    ! &  uL,sec(sL)%dnw(eL),BL_thick_inl(ID),um)
            ! end if
        elseif(bct .eq. BCInflowSupersonic) then
            um(1:5) =  InflowSupersonic(1:5,ID)
        elseif(bct .eq. BCOutflow) then
            write(*,*),'BCOutflow boundary condition has not been verified'
            ! call bnd_outflow(n, uL, fv(sR)%bc(1,eR), um)
        elseif(bct .eq. BCOutflowSupersonic) then
            um(1:5) =  uL(1:5)
        else
            stop 'Error: boundary type not supported.'
        end if
        fv_u (1:5,idxR)   =  um(1:5)
        fv_uR(1:5,idxR)   =  2.0d0*um(1:5)-uL(1:5)
    end do

! call fv_array_to_struct(lev)

    return
    end subroutine fv_bnd_sw
!-------------------------------------------------------------------------------
!   boundary treatment for FV in the parallel environment..
!-------------------------------------------------------------------------------
    subroutine fv_bnd_parallel(lev,is_uc_to_u)
    use var_kind_def
    use var_fv
    use var_fv_array
    use var_mesh
    use var_parallel
    use var_global, only: sw_slave
    implicit none
    logical(dpL),intent(in):: is_uc_to_u
    integer(dpI),intent(in):: lev
    integer(dpI):: isr,isec,iele,i,idx,ip_remote,ID
    ! character(20):: get_gra='get_gra'
    ! character(20):: mpi_waitall_t = 'mpi_waitall_t'

!   ----------------------------------------------------------------------------
!   get u from uc for the elements used in the data exchange.
    if(is_uc_to_u)  call fv_lev_uc2u(lev, 1)
!   get u from uc for the elements used in the data exchange.
!   ----------------------------------------------------------------------------
! call starttime(mpi_waitall_t)
  ! ----------------------------------------------------------------------------
!   prepare data.
    do isr=mesh(lev)%p2p_1,mesh(lev)%p2p_0
        idx =  1
        do i=1,p2p(isr)%n_ele_send
            isec=  p2p(isr)%id_ele_send(1,i)
            iele=  p2p(isr)%id_ele_send(2,i)
            p2p(isr)%rsend(idx:idx+4)   =  fv(isec)%u(1:5,iele)
            idx =  idx+5
        end do
        p2p(isr)%n_send =  idx-1
        p2p(isr)%n_recv =  5*p2p(isr)%n_ele_recv
    end do
!   prepare data.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   data exchange.
    mpi_nreq=  0
    do isr=mesh(lev)%p2p_1,mesh(lev)%p2p_0
        ip_remote   =  p2p(isr)%ip_remote
        if(myid .eq. ip_remote) cycle

        if(p2p(isr)%n_send .gt. 0) then
            mpi_nreq=  mpi_nreq+1
            if(mpi_nreq .gt. max_p2p)   stop 'Error: you are working with a bad LB.'
            call mpi_isend(p2p(isr)%rsend, p2p(isr)%n_send, mpi_dpR, ip_remote, &
                &  myid, mpi_comm_world, mpi_req(mpi_nreq), mpi_err)
        end if

        if(p2p(isr)%n_recv .gt. 0) then
            mpi_nreq=  mpi_nreq+1
            if(mpi_nreq .gt. max_p2p)   stop 'Error: you are working with a bad LB.'
            call mpi_irecv(p2p(isr)%rrecv, p2p(isr)%n_recv, mpi_dpR, ip_remote, &
                &  ip_remote, mpi_comm_world, mpi_req(mpi_nreq), mpi_err)
        end if
    end do
!   data exchange.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   get u from uc for the elements not used in the data exchange.
    if(is_uc_to_u)  call fv_lev_uc2u(lev, 2)
!   get u from uc for the elements not used in the data exchange.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   physical boundary.
    call fv_mortar_get_u(lev, 1)
    call fv_bnd(lev)
!   physical boundary.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   data exchange on myid.
    do isr=mesh(lev)%p2p_1,mesh(lev)%p2p_0
        ip_remote   =  p2p(isr)%ip_remote
        if(myid .ne. ip_remote) cycle
        if(p2p(isr)%n_send .ne. p2p(isr)%n_recv)    stop 'Error: s&r not match on myid.'

        if(p2p(isr)%n_send .gt. 0)  call DCOPY(p2p(isr)%n_send, p2p(isr)%rsend, 1, &
            &  p2p(isr)%rrecv, 1)
    end do
!   data exchange on myid.
!   ----------------------------------------------------------------------------

    call mpi_waitall(mpi_nreq, mpi_req, mpi_sta, mpi_err)

!   ----------------------------------------------------------------------------
!   unpack data received.
    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(.not. sec(isec)%is_ghost)    cycle
        do i=1,sec(isec)%n_ele
            isr =  sec(isec)%id_recv(1,i)
            iele=  sec(isec)%id_recv(2,i)
            call per_rot_vec(1, sec(isec)%per_path(1,i), p2p(isr)%rrecv(5*iele-4), &
                &  fv(isec)%u(1,i))
        end do
    end do
!   unpack data received.
!   ----------------------------------------------------------------------------
! call endtime(mpi_waitall_t)

    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(sec(isec)%is_int)    cycle

        call u_to_uc(sec(isec)%n_ele, fv(isec)%u, fv(isec)%uc, fv(isec)%t, fv(isec)%mu)
    end do

! call starttime(get_gra)
    if(lev .eq. 0)  call fv_get_gra(lev)
! call endtime(get_gra)

    return
    end subroutine fv_bnd_parallel
!-------------------------------------------------------------------------------
!   boundary treatment for FV in the parallel environment. Sunway version.
!-------------------------------------------------------------------------------
    subroutine fv_bnd_parallel_sw(lev,is_uc_to_u)
    use var_kind_def
    use var_fv
    use var_fv_array
    use var_mesh
    use var_parallel
    use var_global, only: sw_slave
    implicit none
    logical(dpL),intent(in):: is_uc_to_u
    integer(dpI),intent(in):: lev
    integer(dpI):: isr,isec,iele,i,idx,ip_remote,ID
    character(20):: get_gra='get_gra'
    character(20):: mpi_waitall_t = 'mpi_waitall_t'
    character(20):: mpi_rot_vec_t = 'mpi_rot_vec_t'
    character(20):: fv_bnd_parallel_t = 'fv_bnd_parallel_t'
    character(20):: fv_lev_uc2u_t = 'fv_lev_uc2u_t'
    character(20):: fv_lev_u2uc_t = 'fv_lev_u2uc_t'

! call fv_struct_to_array(lev)
call starttime(fv_bnd_parallel_t)

call starttime(fv_lev_uc2u_t)
    if(is_uc_to_u) call fv_lev_uc2u_sw(lev)
call endtime(fv_lev_uc2u_t)

call starttime(mpi_waitall_t)
!   ----------------------------------------------------------------------------
!   prepare data.
    do isr=mesh(lev)%p2p_1,mesh(lev)%p2p_0
        idx =  1
        do i=1,p2p(isr)%n_ele_send
            isec=  p2p(isr)%id_ele_send(1,i)
            iele=  p2p(isr)%id_ele_send(2,i)
            ID = sec(isec)%ID_ele_g(iele)
            if(ID .le. cellNum) ID = perm(ID)
            ! p2p(isr)%rsend(idx:idx+4)   =  fv(isec)%u(1:5,iele)
            p2p(isr)%rsend(idx:idx+4)   =  fv_u(1:5,ID)
            idx =  idx+5
        end do
        p2p(isr)%n_send =  idx-1
        p2p(isr)%n_recv =  5*p2p(isr)%n_ele_recv
    end do
!   prepare data.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   data exchange.
    mpi_nreq=  0
    do isr=mesh(lev)%p2p_1,mesh(lev)%p2p_0
        ip_remote   =  p2p(isr)%ip_remote
        if(myid .eq. ip_remote) cycle

        if(p2p(isr)%n_send .gt. 0) then
            mpi_nreq=  mpi_nreq+1
            if(mpi_nreq .gt. max_p2p)   stop 'Error: you are working with a bad LB.'
            call mpi_isend(p2p(isr)%rsend, p2p(isr)%n_send, mpi_dpR, ip_remote, &
                &  myid, mpi_comm_world, mpi_req(mpi_nreq), mpi_err)
        end if

        if(p2p(isr)%n_recv .gt. 0) then
            mpi_nreq=  mpi_nreq+1
            if(mpi_nreq .gt. max_p2p)   stop 'Error: you are working with a bad LB.'
            call mpi_irecv(p2p(isr)%rrecv, p2p(isr)%n_recv, mpi_dpR, ip_remote, &
                &  ip_remote, mpi_comm_world, mpi_req(mpi_nreq), mpi_err)
        end if
    end do
!   data exchange.

!   ----------------------------------------------------------------------------
!   data exchange on myid.
    do isr=mesh(lev)%p2p_1,mesh(lev)%p2p_0
        ip_remote   =  p2p(isr)%ip_remote
        if(myid .ne. ip_remote) cycle
        if(p2p(isr)%n_send .ne. p2p(isr)%n_recv)    stop 'Error: s&r not match on myid.'

        if(p2p(isr)%n_send .gt. 0)  call DCOPY(p2p(isr)%n_send, p2p(isr)%rsend, 1, &
            &  p2p(isr)%rrecv, 1)
    end do
!   data exchange on myid.
!   ----------------------------------------------------------------------------
!   ----------------------------------------------------------------------------
!   physical boundary.
    ! call fv_mortar_get_u(lev, 1)
    call fv_mortar_get_u_sw(lev, 1)
    call fv_bnd_sw(lev)
!   physical boundary.
!   ----------------------------------------------------------------------------

    call mpi_waitall(mpi_nreq, mpi_req, mpi_sta, mpi_err)
!   ----------------------------------------------------------------------------
!   unpack data received.
    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(.not. sec(isec)%is_ghost)    cycle
        do i=1,sec(isec)%n_ele
            isr =  sec(isec)%id_recv(1,i)
            iele=  sec(isec)%id_recv(2,i)
            ID = sec(isec)%ID_ele_g(i)
            if(ID .le. cellNum) ID = perm(ID)
            call per_rot_vec(1, sec(isec)%per_path(1,i), p2p(isr)%rrecv(5*iele-4), &
                &  fv_u(1,ID))
        end do
    end do
!   unpack data received.
!   ----------------------------------------------------------------------------
call endtime(mpi_waitall_t)
! call fv_array_to_struct(lev)

    ! do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
    !     if(sec(isec)%is_int)    cycle

    !     call u_to_uc(sec(isec)%n_ele, fv(isec)%u, fv(isec)%uc, fv(isec)%t, fv(isec)%mu)
    ! end do

call starttime(fv_lev_u2uc_t)
    call fv_lev_u2uc_sw(lev)
call endtime(fv_lev_u2uc_t)

call starttime(get_gra)
    if(lev .eq. 0)  call fv_get_gra_sw(lev)
call endtime(get_gra)

call endtime(fv_bnd_parallel_t)
! call fv_array_to_struct(lev)

    return
    end subroutine fv_bnd_parallel_sw
!-------------------------------------------------------------------------------
!   get L and R flow variables for the mortar element.
!-------------------------------------------------------------------------------
    subroutine fv_mortar_get_u(lev,mode)
    use var_kind_def
    use var_cgns
    use var_fv
    use var_mesh
    implicit none
    integer(dpI),intent(in):: lev,mode
    integer(dpI):: im,sL,eL,sR,eR,im1,im0

    if(mode .le. 0) then
        im1 =  1
        im0 =  mesh(lev)%n_mortar
    elseif(mode .eq. 1) then
        im1 =  1
        im0 =  mesh(lev)%n_mortar_b
    else
        im1 =  mesh(lev)%n_mortar_b+1
        im0 =  mesh(lev)%n_mortar
    end if
    do im=im1,im0
        sL  =  mesh(lev)%mortar_LR(1,im)
        eL  =  mesh(lev)%mortar_LR(2,im)
        sR  =  mesh(lev)%mortar_LR(3,im)
        eR  =  mesh(lev)%mortar_LR(4,im)
        if(.not. sec(sL)%is_bnd)    fv(sR)%uL(1:5,eR)   =  fv(sL)%u(1:5,eL)
        if(.not. sec(sR)%is_bnd)    fv(sR)%uR(1:5,eR)   =  fv(sR)%u(1:5,eR)
    end do

    return
    end subroutine fv_mortar_get_u
!-------------------------------------------------------------------------------
!   get L and R flow variables for the mortar element. Sunway version
!-------------------------------------------------------------------------------
    subroutine fv_mortar_get_u_sw(lev,mode)
    use var_kind_def
    use var_cgns
    use var_fv
    use var_fv_array
    use var_mesh
    use var_sec_array
    implicit none
    integer(dpI),intent(in):: lev,mode
    integer(dpI):: im,sL,eL,sR,eR,im1,im0,idxL,idxR

! call fv_struct_to_array(lev)

    if(mode .le. 0) then
        im1 =  1
        im0 =  mesh(lev)%n_mortar
    elseif(mode .eq. 1) then
        im1 =  1
        im0 =  mesh(lev)%n_mortar_b
    else
        im1 =  mesh(lev)%n_mortar_b+1
        im0 =  mesh(lev)%n_mortar
    end if
    do im=im1,im0
        sL      =  mesh(lev)%mortar_LR(1,im)
        eL      =  mesh(lev)%mortar_LR(2,im)
        sR      =  mesh(lev)%mortar_LR(3,im)
        eR      =  mesh(lev)%mortar_LR(4,im)
        idxL =  sec(sL)%ID_ele_g(eL)
        idxR =  sec(sR)%ID_ele_g(eR)
        if(idxL .lt. cellNum) idxL = perm(idxL)
        if(idxR .lt. cellNum) idxR = perm(idxR)
        if(sec_is_bnd(idxL) .ne. 1)    fv_uL(1:5,idxR)   =  fv_u(1:5,idxL)
        if(sec_is_bnd(idxR) .ne. 1)    fv_uR(1:5,idxR)   =  fv_u(1:5,idxR)
        ! if(.not. sec(sL)%is_bnd)    fv(sR)%uL(1:5,eR)   =  fv(sL)%u(1:5,eL)
        ! if(.not. sec(sR)%is_bnd)    fv(sR)%uR(1:5,eR)   =  fv(sR)%u(1:5,eR)
    end do

! call fv_array_to_struct(lev)

    return
    end subroutine fv_mortar_get_u_sw
!-------------------------------------------------------------------------------
!   set the boundary condition.
!-------------------------------------------------------------------------------
    subroutine fv_set_boundary_condition
    use var_kind_def
    use var_bndv
    use var_cgns, only: BCInflow,BCOutflow
    use var_fv
    use var_global, only: z=>rotation_axis
    use var_mesh
    use var_slv, only: lev0_cal
    implicit none
    integer(dpI):: lev,isec,bct,ID,iele
    real   (dpR):: r(3),rc,c(3),a(3),t(3),rtmp

    do lev=lev0_cal,mlev_fv-1
    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(.not. sec(isec)%is_bnd)  cycle
        bct =  sec(isec)%bct
        ID  =  sec(isec)%ID_group
        if(bct .eq. BCInflow) then
            if(.not. allocated(fv(isec)%bc))    allocate(fv(isec)%bc(5,sec(isec)%n_ele))
        elseif(bct .eq. BCOutflow) then
            if(.not. allocated(fv(isec)%bc))    allocate(fv(isec)%bc(1,sec(isec)%n_ele))
        else
            cycle
        end if

        do iele=1,sec(isec)%n_ele
            c(1:3)  =  sec(isec)%cen(1:3,iele)
            r       =  c-(c(1)*z(1)+c(2)*z(2)+c(3)*z(3))*z
            call norm_vec(3, r, rc)

            if(bct .eq. BCInflow) then
                if(is_inl_cyl(ID)) then
                    call crs_prd(z, r, t)
                    if(sec(isec)%is_profiled) then
                        call get_bnd_profile_1d(isec, rc, fv(isec)%bc(1,iele))
                        a(1:3)  =  fv(isec)%bc(3:5,iele)
                    else
                        a(1:3)  =  vdir_inl(1:3,ID)
                    end if
                    fv(isec)%bc(3:5,iele)   =  a(1)*r+a(2)*t+a(3)*z
                else
                    if(sec(isec)%is_profiled) then
                        stop 'Error: INL boundary condition profile not supported.'
                    else
                        fv(isec)%bc(1  ,iele)   =  Tt_inl(ID)
                        fv(isec)%bc(2  ,iele)   =  Pt_inl(ID)
                        fv(isec)%bc(3:5,iele)   =  vdir_inl(1:3,ID)
                    end if
                end if
                call norm_vec(3, fv(isec)%bc(3,iele), rtmp)
            else
                if(sec(isec)%is_profiled) then
                    call get_bnd_profile_1d(isec, rc, fv(isec)%bc(1,iele))
                else
                    fv(isec)%bc(1,iele) =  pb_out(ID)
                end if
            end if
        end do
    end do
    end do

    call fv_set_reot

    return
    end subroutine fv_set_boundary_condition
