!-------------------------------------------------------------------------------
!   setup the explicit solving.
!-------------------------------------------------------------------------------
    subroutine fv_set_exp_slv(lev)
    use var_kind_def
    use var_fv
    use var_mesh
    implicit none
    integer(dpI),intent(in):: lev
    integer(dpI):: isec

    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(sec(isec)%is_ghost .or. sec(isec)%is_bnd)    cycle
        if(.not. allocated(fv(isec)%LHS_s)) allocate(fv(isec)%LHS_s(sec(isec)%n_ele))
    end do

    return
    end subroutine fv_set_exp_slv
!-------------------------------------------------------------------------------
!   compute the LHS, scalar.
!-------------------------------------------------------------------------------
    subroutine fv_cal_LHS_scalar(lev,is_cal_dt,is_conv_only)
    use var_kind_def
    use var_air, only: gk
    use var_fv
    use var_fv_array
    use var_mesh
    use var_slv, only: is_vis_cal,CFL_lev
    use var_turb,only: is_tur_cal
    use var_sec_array
    implicit none
    logical(dpI),intent(in):: is_cal_dt,is_conv_only
    integer(dpI),intent(in):: lev
    integer(dpI):: isec,i,im,sL,eL,sR,eR,idxL,idxR,imi
    real   (dpR):: CFL,u(5),a,n(5),s,eps,CFLu,uc(5),R(5),mu,mut
    character(20):: fv_cal_LHS_scalar_t = 'fv_cal_LHS_scalar_t'


    eps =  1.0d-3
    CFL =  CFL_lev(lev)

    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(sec(isec)%is_int)    fv(isec)%LHS_s  =  0.0d0
    end do

    do im=1,mesh(lev)%n_mortar
        sL  =  mesh(lev)%mortar_LR(1,im)
        eL  =  mesh(lev)%mortar_LR(2,im)
        sR  =  mesh(lev)%mortar_LR(3,im)
        eR  =  mesh(lev)%mortar_LR(4,im)

        if(im .le. mesh(lev)%n_mortar_b) then
            u(1:5)  =  fv(sR)%u(1:5,eR)
            if(is_vis_cal)  mu  =  fv(sR)%mu(1,eR)
            if(is_tur_cal)  mut =  fv(sR)%mu(2,eR)
        else
            u(1:5)  = (fv(sL)%u(1:5,eL)+fv(sR)%u(1:5,eR))*0.5d0
            if(is_vis_cal)  mu  = (fv(sL)%mu(1,eL)+fv(sR)%mu(1,eR))*0.5d0
            if(is_tur_cal)  mut = (fv(sL)%mu(2,eL)+fv(sR)%mu(2,eR))*0.5d0            
        end if
        a       =  sqrt(gk*u(5)/u(1))
        n(1:5)  =  mesh(lev)%mortar_n_vg(1:5,im)
        u(2:4)  =  u(2:4)-mesh(lev)%mortar_n_vg(6:8,im)
        s       =  0.5d0*n(4)*(abs(u(2)*n(1)+u(3)*n(2)+u(4)*n(3))+a)

        if(.not. is_conv_only) then
            if(is_vis_cal)  s   =  s+gk*mu *n(4)*n(5)/u(1)
            if(is_tur_cal)  s   =  s+gk*mut*n(4)*n(5)/u(1)
        end if

        fv(sL)%LHS_s(eL)=  fv(sL)%LHS_s(eL)+s
        if(sec(sR)%is_int)  fv(sR)%LHS_s(eR)=  fv(sR)%LHS_s(eR)+s
    end do

    if(.not. is_cal_dt) return
    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(sec(isec)%is_ghost .or. sec(isec)%is_bnd)    cycle
        do i=1,sec(isec)%n_ele
            s       =  fv(isec)%LHS_s(  i)
            uc(1:5) =  fv(isec)%uc (1:5,i)
            R (1:5) =  fv(isec)%rhs(1:5,i)
            CFLu    =  CFL
            CFLu    =  min(CFLu, eps*uc(1)*s/(abs(R(1))+1.0d-9))
            CFLu    =  min(CFLu, eps*uc(5)*s/(abs(R(5))+1.0d-9))
            CFLu    =  max(CFLu, 1.0d-1*CFL)

            fv(isec)%LHS_s(i)   =  CFLu/s
        end do
    end do

    return
    end subroutine fv_cal_LHS_scalar
!-------------------------------------------------------------------------------
!   compute the LHS, scalar.
!-------------------------------------------------------------------------------
    subroutine fv_cal_LHS_scalar_sw(lev,is_cal_dt,is_conv_only)
    use var_kind_def
    use var_air, only: gk
    use var_fv
    use var_fv_array
    use var_mesh
    use var_slv, only: is_vis_cal,CFL_lev
    use var_turb,only: is_tur_cal
    use var_sec_array
    use var_global_real
    implicit none
    logical(dpI),intent(in):: is_cal_dt,is_conv_only
    integer(dpI),intent(in):: lev
    integer(dpI):: isec,i,im,sL,eL,sR,eR,idxL,idxR,imi,iele
    real   (dpR):: CFL,u(5),a,n(5),s,eps,CFLu,uc(5),R(5),mu,mut
    character(20):: fv_cal_LHS_scalar_t = 'fv_cal_LHS_scalar_t'
    real   (dpR):: is_conv_only_r

call starttime(fv_cal_LHS_scalar_t)

    eps =  1.0d-3
    CFL =  CFL_lev(lev)

    ! do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
    !     if(sec(isec)%is_int)    fv(isec)%LHS_s  =  0.0d0
    ! end do
    ! fv_LHS_s = 0.0d0
    call lhs_scalar_zero_host(tot_ele, fv_LHS_s)

! call fv_struct_to_array(lev)
    do im=1,mesh(lev)%n_mortar_b
        sL  =  mesh(lev)%mortar_LR(1,im)
        eL  =  mesh(lev)%mortar_LR(2,im)
        sR  =  mesh(lev)%mortar_LR(3,im)
        eR  =  mesh(lev)%mortar_LR(4,im)

        idxL = sec(sL)%ID_ele_g(eL)
        idxR = sec(sR)%ID_ele_g(eR)
        if(idxL .le. cellNum) idxL = perm(idxL)
        if(idxR .le. cellNum) idxR = perm(idxR)

        u(1:5)  =  fv_u(1:5,idxR)
        if(is_vis_cal)  mu  =  fv_mu(1,idxR)
        if(is_tur_cal)  mut =  fv_mu(2,idxR)
        a       =  sqrt(gk*u(5)/u(1))
        n(1:5)  =  mesh(lev)%mortar_n_vg(1:5,im)
        u(2:4)  =  u(2:4)-mesh(lev)%mortar_n_vg(6:8,im)
        s       =  0.5d0*n(4)*(abs(u(2)*n(1)+u(3)*n(2)+u(4)*n(3))+a)

        if(.not. is_conv_only) then
            if(is_vis_cal)  s   =  s+gk*mu *n(4)*n(5)/u(1)
            if(is_tur_cal)  s   =  s+gk*mut*n(4)*n(5)/u(1)
        end if

        fv_LHS_s(idxL)=  fv_LHS_s(idxL)+s
        if(sec(sR)%is_int)  fv_LHS_s(idxR)=  fv_LHS_s(idxR)+s
    end do

    if(is_conv_only) then
        is_conv_only_r = 1
    else
        is_conv_only_r = -1
    end if
    call lhs_scalar_host(mesh_reordered(lev)%owner, mesh_reordered(lev)%neighbor, &
        & faceNum, cellNum, mesh_reordered(lev)%mortar_transform, &
        & mesh_reordered(lev)%mortar_n_vg, fv_u, fv_mu, fv_LHS_s, sec_is_int, &
        is_vis_cal_r, is_tur_cal_r, is_conv_only_r, gk)

    ! imi = 0
    ! do im=1+mesh(lev)%n_mortar_b,mesh(lev)%n_mortar
    !     imi = imi +1
    !     if(mesh_reordered(lev)%mortar_transform(imi) .eq. 0) then
    !         idxL = mesh_reordered(lev)%mortar_own_ID(imi)
    !         idxR = mesh_reordered(lev)%mortar_nei_ID(imi)
    !     else
    !         idxL = mesh_reordered(lev)%mortar_nei_ID(imi)
    !         idxR = mesh_reordered(lev)%mortar_own_ID(imi)
    !     end if
    !     u(1:5)  = (fv_u(1:5,idxL)+fv_u(1:5,idxR))*0.5d0
    !     if(is_vis_cal)  mu  = (fv_mu(1,idxL)+fv_mu(1,idxR))*0.5d0
    !     if(is_tur_cal)  mut = (fv_mu(2,idxL)+fv_mu(2,idxR))*0.5d0
    !     a       =  sqrt(gk*u(5)/u(1))
    !     n(1:5)  =  mesh_reordered(lev)%mortar_n_vg(1:5,imi)
    !     u(2:4)  =  u(2:4)-mesh_reordered(lev)%mortar_n_vg(6:8,imi)
    !     s       =  0.5d0*n(4)*(abs(u(2)*n(1)+u(3)*n(2)+u(4)*n(3))+a)

    !     if(.not. is_conv_only) then
    !         if(is_vis_cal)  s   =  s+gk*mu *n(4)*n(5)/u(1)
    !         if(is_tur_cal)  s   =  s+gk*mut*n(4)*n(5)/u(1)
    !     end if

    !     fv_LHS_s(idxL)=  fv_LHS_s(idxL)+s
    !     if(sec_is_int(idxR) .eq. 1)  fv_LHS_s(idxR)=  fv_LHS_s(idxR)+s
    ! end do

    if(.not. is_cal_dt) return
    ! do iele = 1,tot_ele
    !     if(sec_is_int(iele) .eq. 1) then
    !         s       = fv_LHS_s(iele)
    !         uc(1:5) = fv_uc(1:5,iele)
    !         R (1:5) = fv_rhs(1:5,iele)
    !         CFLu    = CFL
    !         CFLu    = min(CFLu, eps*uc(1)*s/(abs(R(1))+1.0d-9))
    !         CFLu    =  min(CFLu, eps*uc(5)*s/(abs(R(5))+1.0d-9))
    !         CFLu    =  max(CFLu, 1.0d-1*CFL)

    !         ! fv_LHS_s(iele) = CFLu/s
    !     end if
    ! end do

    call lhs_scalar_cal_dt_host(tot_ele, fv_LHS_s, fv_uc, fv_rhs, sec_is_int, &
        & CFL, eps)
call endtime(fv_cal_LHS_scalar_t)
! call fv_array_to_struct(lev)
    ! do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
    !     if(sec(isec)%is_ghost .or. sec(isec)%is_bnd)    cycle
    !     do i=1,sec(isec)%n_ele
    !         s       =  fv(isec)%LHS_s(  i)
    !         uc(1:5) =  fv(isec)%uc (1:5,i)
    !         R (1:5) =  fv(isec)%rhs(1:5,i)
    !         CFLu    =  CFL
    !         CFLu    =  min(CFLu, eps*uc(1)*s/(abs(R(1))+1.0d-9))
    !         CFLu    =  min(CFLu, eps*uc(5)*s/(abs(R(5))+1.0d-9))
    !         CFLu    =  max(CFLu, 1.0d-1*CFL)

    !         fv(isec)%LHS_s(i)   =  CFLu/s
    !     end do
    ! end do

    return
    end subroutine fv_cal_LHS_scalar_sw
!-------------------------------------------------------------------------------
!   RK(4,3) time marching, SSP.
!-------------------------------------------------------------------------------
    subroutine fv_slv_exp_rk43_ssp(lev)
    use var_kind_def
    use var_fv
    use var_global, only: R16,R13,R23,err_mem,sw_slave
    use var_mesh
    use var_slv, only: is_local_dt,is_st_cal_now,lev_out,RK
    implicit none
    integer(dpI),intent(in):: lev
    integer(dpI):: isec,ip
    character(20):: fv_slv_exp_rk43_ssp_t = 'fv_slv_exp_rk43_ssp_t'

call starttime(fv_slv_exp_rk43_ssp_t)
    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(.not. allocated(fv(isec)%uc0))  then
            allocate(fv(isec)%uc0(5, sec(isec)%n_ele), stat=err_mem)
        endif
        if(sec(isec)%is_int)    fv(isec)%uc0=  fv(isec)%uc
    end do

!   ----------------------------------------------------------------------------
!   the first sub-iteration.
    RK  =  1
    call fv_get_rhs(lev, is_st_cal_now)
    if(lev .eq. lev_out)    call fv_get_residual(lev, res_NS)
    if(is_local_dt) call fv_cal_LHS_scalar(lev, .true., .false.)
    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(sec(isec)%is_ghost .or. sec(isec)%is_bnd)    cycle
        do ip=1,sec(isec)%n_ele
            fv(isec)%uc(:,ip)   =  fv(isec)%uc0(:,ip) &
                                & -0.5d0*fv(isec)%rhs(:,ip)*fv(isec)%LHS_s(ip)
        end do
    end do
    call fv_bnd_parallel(lev, .true.)
!   the first sub-iteration.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   the second sub-iteration.
    RK  =  2
    call fv_get_rhs(lev, is_st_cal_now)
    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(sec(isec)%is_ghost .or. sec(isec)%is_bnd)    cycle
        do ip=1,sec(isec)%n_ele
            fv(isec)%uc(:,ip)   =  fv(isec)%uc(:,ip) &
                                & -0.5d0*fv(isec)%rhs(:,ip)*fv(isec)%LHS_s(ip)
        end do
    end do
    call fv_bnd_parallel(lev, .true.)
!   the second sub-iteration.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   the third sub-iteration.
    RK  =  3
    call fv_get_rhs(lev, is_st_cal_now)
    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(sec(isec)%is_ghost .or. sec(isec)%is_bnd)    cycle
        do ip=1,sec(isec)%n_ele
            fv(isec)%uc(:,ip)   =  R23*fv(isec)%uc0(:,ip)+R13*fv(isec)%uc(:,ip) &
                                & -R16*fv(isec)%rhs(:,ip)*fv(isec)%LHS_s(ip)
        end do
    end do
    call fv_bnd_parallel(lev, .true.)
!   the third sub-iteration.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   the fourth sub-iteration.
    RK  =  4
    call fv_get_rhs(lev, is_st_cal_now)
    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(sec(isec)%is_ghost .or. sec(isec)%is_bnd)    cycle
        do ip=1,sec(isec)%n_ele
            fv(isec)%uc(:,ip)   =  fv(isec)%uc(:,ip) &
                                & -0.5d0*fv(isec)%rhs(:,ip)*fv(isec)%LHS_s(ip)
        end do
    end do
    call fv_bnd_parallel(lev, .true.)
!   the fourth sub-iteration.
!   ----------------------------------------------------------------------------
call endtime(fv_slv_exp_rk43_ssp_t)

    return
    end subroutine fv_slv_exp_rk43_ssp
!-------------------------------------------------------------------------------
!   RK(4,3) time marching, SSP. Sunway version
!-------------------------------------------------------------------------------
    subroutine fv_slv_exp_rk43_ssp_sw(lev)
    use var_kind_def
    use var_fv
    use var_fv_array
    use var_sec_array
    use var_global, only: R16,R13,R23,err_mem,sw_slave
    use var_mesh
    use var_slv, only: is_local_dt,is_st_cal_now,lev_out,RK
    implicit none
    integer(dpI),intent(in):: lev
    integer(dpI):: isec,ip
    character(20):: fv_slv_exp_rk43_ssp_t = 'fv_slv_exp_rk43_ssp_t'

call fv_struct_to_array(lev)
call starttime(fv_slv_exp_rk43_ssp_t)

    ! fv_uc0 = fv_uc
    call exp_rk43_copy_host(tot_ele, fv_uc0, fv_uc)

!   ----------------------------------------------------------------------------
!   the first sub-iteration.
    RK  =  1
    call fv_get_rhs_sw(lev, is_st_cal_now)
    if(lev .eq. lev_out)    call fv_get_residual_sw(lev, res_NS)
    if(is_local_dt) call fv_cal_LHS_scalar_sw(lev, .true., .false.)

    call exp_rk43_1st_host(tot_ele, fv_uc0, fv_uc, fv_rhs, fv_LHS_s, sec_is_int)

    call fv_bnd_parallel_sw(lev, .true.)
!   the first sub-iteration.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   the second sub-iteration.
    RK  =  2
    call fv_get_rhs_sw(lev, is_st_cal_now)
    ! do ip=1,tot_ele
    !     if(sec_is_int(ip) .eq. 1) then
    !         fv_uc(:,ip) = fv_uc(:,ip)-0.5d0*fv_rhs(:,ip)*fv_LHS_s(ip)
    !     end if
    ! end do
    call exp_rk43_2nd_host(tot_ele, fv_uc, fv_rhs, fv_LHS_s, sec_is_int)

    call fv_bnd_parallel_sw(lev, .true.)
!   the second sub-iteration.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   the third sub-iteration.
    RK  =  3
    call fv_get_rhs_sw(lev, is_st_cal_now)
    ! do ip=1,tot_ele
    !     if(sec_is_int(ip) .eq. 1) then
    !         fv_uc(:,ip) = R23*fv_uc0(:,ip)+R13*fv_uc(:,ip)-R16*fv_rhs(:,ip)*fv_LHS_s(ip)
    !     end if
    ! end do
    call exp_rk43_3rd_host(tot_ele, fv_uc, fv_uc0, fv_rhs, fv_LHS_s, sec_is_int, &
        & R23, R13, R16)
    call fv_bnd_parallel_sw(lev, .true.)
!   the third sub-iteration.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   the fourth sub-iteration.
    RK  =  4
    call fv_get_rhs_sw(lev, is_st_cal_now)
    ! do ip=1,tot_ele
    !     if(sec_is_int(ip) .eq. 1) then
    !         fv_uc(:,ip) = fv_uc(:,ip)-0.5d0*fv_rhs(:,ip)*fv_LHS_s(ip)
    !     end if
    ! end do
    call exp_rk43_2nd_host(tot_ele, fv_uc, fv_rhs, fv_LHS_s, sec_is_int)
    
    call fv_bnd_parallel_sw(lev, .true.)
!   the fourth sub-iteration.
!   ----------------------------------------------------------------------------
call endtime(fv_slv_exp_rk43_ssp_t)
call fv_array_to_struct(lev)

    return
    end subroutine fv_slv_exp_rk43_ssp_sw
