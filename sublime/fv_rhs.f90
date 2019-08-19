!-------------------------------------------------------------------------------
!   cal Rhs.
!-------------------------------------------------------------------------------
    subroutine fv_get_rhs(lev,is_ff)
    use var_kind_def
    use var_fv
    use var_mesh
    use var_slv, only: is_vis_cal,lev_out
    use var_turb, only: is_LES_now,LES_model
    use var_uns_cal, only: is_bdf_now
    use var_global, only: sw_slave
    implicit none
    logical(dpL),intent(in):: is_ff
    integer(dpI),intent(in):: lev
    integer(dpI):: isec
    character(20):: fv_get_rhs_t = 'fv_get_rhs'
    character(20):: fv_get_rhs_upw_t = 'fv_get_rhs_upw_t'
    character(20):: fv_get_rhs_vis_t = 'fv_get_rhs_vis_t'

call starttime(fv_get_rhs_t)
    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(sec(isec)%is_int)    fv(isec)%rhs=  0.0d0
    end do

    if(is_LES_now .and. (LES_model .ge. 1)) call fv_les(lev)
    if(is_shock_sensor) call fv_get_shock_sensor(lev)

    if(is_kep) then
        call fv_get_rhs_kep(lev)
    else
        if(rhs_lev(lev) .eq. 1) then
            if(is_jst_su2) then
                call fv_get_rhs_jst_su2(lev, .true., .true.)
            else
                call fv_get_rhs_jst    (lev, .true., .true.)
            end if
        else
call starttime(fv_get_rhs_upw_t)
            call fv_get_rhs_upw(lev, .true., .true.)
call endtime(fv_get_rhs_upw_t)
        end if
    end if
    call fv_get_rhs_src_rotation(lev)
    if(is_bdf_now)  call fv_get_src_bdf(lev)

call starttime(fv_get_rhs_vis_t)
    if(is_vis_cal)  call fv_get_rhs_vis(lev, .true.)
call endtime(fv_get_rhs_vis_t)

    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(sec(isec)%is_ghost .or. sec(isec)%is_bnd)    cycle
        if((lev .gt. lev_out) .and. is_ff)  fv(isec)%rhs=  fv(isec)%rhs+fv(isec)%ff
    end do
call endtime(fv_get_rhs_t)

    return
    end subroutine fv_get_rhs
!-------------------------------------------------------------------------------
!   cal Rhs. Sunway version.
!-------------------------------------------------------------------------------
    subroutine fv_get_rhs_sw(lev,is_ff)
    use var_kind_def
    use var_fv
    use var_fv_array
    use var_mesh
    use var_slv, only: is_vis_cal,lev_out
    use var_turb, only: is_LES_now,LES_model
    use var_uns_cal, only: is_bdf_now
    use var_global, only: sw_slave
    implicit none
    logical(dpL),intent(in):: is_ff
    integer(dpI),intent(in):: lev
    integer(dpI):: isec
    character(20):: fv_get_rhs_t = 'fv_get_rhs'
    character(20):: fv_get_rhs_upw_t = 'fv_get_rhs_upw_t'
    character(20):: fv_get_rhs_vis_t = 'fv_get_rhs_vis_t'

! call fv_struct_to_array(lev)
call starttime(fv_get_rhs_t)
    fv_rhs = 0.0d0

    if(is_LES_now .and. (LES_model .ge. 1)) stop 'fv_les has not been implemented on Sunway platform'
    if(is_shock_sensor) stop 'fv_get_shock_sensor has not been implemented on Sunway platform'

    if(is_kep) then
        stop 'fv_get_rhs_kep has not been implemented on Sunway platform'
        call fv_get_rhs_kep(lev)
    else
        if(rhs_lev(lev) .eq. 1) then
            if(is_jst_su2) then
                stop 'fv_get_rhs_jst_su2 has not been implemented on Sunway platform'
                call fv_get_rhs_jst_su2(lev, .true., .true.)
            else
                stop 'fv_get_rhs_jst has not been implemented on Sunway platform'
                call fv_get_rhs_jst    (lev, .true., .true.)
            end if
        else
call starttime(fv_get_rhs_upw_t)
            call fv_get_rhs_upw_sw(lev, .true., .true.)
call endtime(fv_get_rhs_upw_t)
        end if
    end if
    ! call fv_get_rhs_src_rotation(lev)
    if(is_bdf_now)  stop 'fv_get_src_bdf has not been implemented on Sunway platform'

call starttime(fv_get_rhs_vis_t)
    if(is_vis_cal)  call fv_get_rhs_vis_sw(lev, .true.)
call endtime(fv_get_rhs_vis_t)

    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(sec(isec)%is_ghost .or. sec(isec)%is_bnd)    cycle
        if((lev .gt. lev_out) .and. is_ff)  stop 'fv_ff has not been implemented on Sunway platform'
    end do
call endtime(fv_get_rhs_t)

! call fv_array_to_struct(lev)

    return
    end subroutine fv_get_rhs_sw
!-------------------------------------------------------------------------------
!   cal residual.
!-------------------------------------------------------------------------------
    subroutine fv_get_residual(lev,res)
    use var_kind_def
    use var_fv
    use var_mesh
    use var_slv, only: res_ref
    implicit none
    integer(dpI),intent(in):: lev
    integer(dpI):: isec,i
    real   (dpR):: res,r(5)

    res =  0.0d0
    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(sec(isec)%is_ghost .or. sec(isec)%is_bnd)    cycle
        do i=1,sec(isec)%n_ele
            r(1:5)  =  fv(isec)%rhs(1:5,i)*res_ref(1:5)
            res     =  res+r(1)**2+r(2)**2+r(3)**2+r(4)**2+r(5)**2
            if(res .ne. res)    stop 'Error: solution diverges.'
        end do
    end do

    return
    end subroutine fv_get_residual
!-------------------------------------------------------------------------------
!   cal residual. Sunway version.
!-------------------------------------------------------------------------------
    subroutine fv_get_residual_sw(lev,res)
    use var_kind_def
    use var_fv
    use var_fv_array
    use var_mesh
    use var_slv, only: res_ref
    implicit none
    integer(dpI),intent(in):: lev
    integer(dpI):: isec,i,ID
    real   (dpR):: res,r(5)
    character(20):: fv_get_residual_t = 'fv_get_residual_t'

! call fv_struct_to_array(lev)
call starttime(fv_get_residual_t)
    res =  0.0d0
    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(sec(isec)%is_ghost .or. sec(isec)%is_bnd)    cycle
        do i=1,sec(isec)%n_ele
            ID = sec(isec)%ID_ele_g(i)
            if(ID .le. cellNum) ID = perm(ID)
            r(1:5)  =  fv_rhs(1:5,ID)*res_ref(1:5)
            res     =  res+r(1)**2+r(2)**2+r(3)**2+r(4)**2+r(5)**2
            if(res .ne. res)    stop 'Error: solution diverges.'
        end do
    end do
call endtime(fv_get_residual_t)
! call fv_array_to_struct(lev)

    return
    end subroutine fv_get_residual_sw
!-------------------------------------------------------------------------------
!   cal Rhs, viscous part, fv.
!-------------------------------------------------------------------------------
    subroutine fv_get_rhs_vis(lev,is_addD)
    use var_kind_def
    use var_fv
    use var_mesh
    use var_temp, only: uLf,uRf,fac_1d,rhsD
    use var_turb, only: is_tur_cal,is_KO
    implicit none
    logical(dpL),intent(in):: is_addD
    integer(dpI),intent(in):: lev
    integer(dpI):: im,sL,eL,sR,eR,j,k
    real   (dpR):: g(12),u(5),d(3),du(4),tL,tR,nne(3),dr,mu(2),tke(1),rtmp
    character(20):: vis_kernel='vis_kernel'

    d       =  0.0d0
    mu      =  0.0d0
    do im=1,mesh(lev)%n_mortar_b
        sL  =  mesh(lev)%mortar_LR(1,im)
        eL  =  mesh(lev)%mortar_LR(2,im)
        sR  =  mesh(lev)%mortar_LR(3,im)
        eR  =  mesh(lev)%mortar_LR(4,im)

        u  (1 :5 )  =  fv(sR)%u(1:5,eR)
        mu (1    )  =  fv(sR)%mu(1,eR)
        if(is_tur_cal)  mu(2)   =  fv(sR)%mu  (2,eR)
        if(is_KO     )  tke     =  fv(sR)%turb(1,eR)
        g  (1 :9 )  =  fv(sR)%gra(4 :12,eR)
        g  (10:12)  =  fv(sR)%gra(16:18,eR)

        call get_vis_fac(1, mesh(lev)%mortar_n_vg(1,im), u, mu, g, tke, rhsD)

        if(is_addD) then
            fv(sL)%rhs(1:5,eL)  =  fv(sL)%rhs(1:5,eL)-rhsD(1:5,1)
        else
            fv(sL)%duc(1:5,eL)  =  fv(sL)%duc(1:5,eL)+rhsD(1:5,1)
        end if
    end do
call starttime(vis_kernel)
    do im=1+mesh(lev)%n_mortar_b,mesh(lev)%n_mortar
        sL  =  mesh(lev)%mortar_LR(1,im)
        eL  =  mesh(lev)%mortar_LR(2,im)
        sR  =  mesh(lev)%mortar_LR(3,im)
        eR  =  mesh(lev)%mortar_LR(4,im)
        fac_1d(1:4,1)   =  mesh(lev)%mortar_n_vg(1:4,im)

        uLf(1:5,1)  =  fv(sL)%u(1:5,eL)
        uRf(1:5,1)  =  fv(sR)%u(1:5,eR)
        u  (1:5  )  =  0.5d0*(uLf(1:5,1)+uRf(1:5,1))
        mu (1    )  = (fv(sL)%mu(1,eL)+fv(sR)%mu(1,eR))*0.5d0
        if(is_tur_cal)  mu(2)   = (fv(sL)%mu  (2,eL)+fv(sR)%mu  (2,eR))*0.5d0
        if(is_KO     )  tke     = (fv(sL)%turb(1,eL)+fv(sR)%turb(1,eR))*0.5d0
        tL          =  fv(sL)%t(    eL)
        tR          =  fv(sR)%t(    eR)
        g  (1 :9 )  = (fv(sL)%gra(4 :12,eL)+fv(sR)%gra(4 :12,eR))*0.5d0
        g  (10:12)  = (fv(sL)%gra(16:18,eL)+fv(sR)%gra(16:18,eR))*0.5d0

        d(1:3)      =  sec(sR)%cen(1:3,eR)-sec(sL)%cen(1:3,eL)
        rtmp        =  1.0d0/sqrt(d(1)**2+d(2)**2+d(3)**2)
        d           =  d*rtmp
        du(1:3)     =  rtmp*(uRf(2:4,1)-uLf(2:4,1))
        du(4  )     =  rtmp*(tR-tL)

        nne(1:3)=  fac_1d(1:3,1)/(fac_1d(1,1)*d(1)+fac_1d(2,1)*d(2)+fac_1d(3,1)*d(3))
        do j=1,4
            k       =  3*j-2
            dr      =  du(j)-g(k)*d(1)-g(k+1)*d(2)-g(k+2)*d(3)
            g(k:k+2)=  g(k:k+2)+dr*nne(1:3)
        end do

        call get_vis_fac(1, fac_1d, u, mu, g, tke, rhsD)

        if(is_addD) then
            fv(sL)%rhs(1:5,eL)  =  fv(sL)%rhs(1:5,eL)-rhsD(1:5,1)
            if(sec(sR)%is_int)  fv(sR)%rhs(1:5,eR)  =  fv(sR)%rhs(1:5,eR)+rhsD(1:5,1)
        else
            fv(sL)%duc(1:5,eL)  =  fv(sL)%duc(1:5,eL)+rhsD(1:5,1)
            if(sec(sR)%is_int)  fv(sR)%duc(1:5,eR)  =  fv(sR)%duc(1:5,eR)-rhsD(1:5,1)
        end if
    end do
call endtime(vis_kernel)

    if(.not. is_cht)    return
    do im=1,mesh(lev)%n_mortar_b
        sL  =  mesh(lev)%mortar_LR(1,im)
        eL  =  mesh(lev)%mortar_LR(2,im)
        sR  =  mesh(lev)%mortar_LR(3,im)
        eR  =  mesh(lev)%mortar_LR(4,im)

        if(sec(sR)%is_cht) then
            fv(sL)%rhs(5,eL)=  fv(sL)%rhs(5,eL)-fv(sR)%heat_flux(eR)
        end if
    end do

    return
    end subroutine fv_get_rhs_vis
!-------------------------------------------------------------------------------
!   cal Rhs, viscous part, fv. sunway edition
!-------------------------------------------------------------------------------
    subroutine fv_get_rhs_vis_sw(lev,is_addD)
    use var_kind_def
    use var_fv
    use var_fv_array
    use var_mesh
    use var_sec_array
    use var_temp, only: uLf,uRf,fac_1d,rhsD
    use var_turb, only: is_tur_cal,is_KO
    use var_air, only: cp,prL,prT
    use var_global_real
    use var_global, only: sw_time, sw_slave
    use var_parallel, only: myid
    implicit none
    logical(dpL),intent(in):: is_addD
    integer(dpI),intent(in):: lev
    integer(dpI):: im,sL,eL,sR,eR,j,k,imi,idxL,idxR
    real   (dpR):: g(12),u(5),d(3),du(4),tL,tR,nne(3),dr,mu(2),tke(1),rtmp
    integer(kind=8):: t1,t2,t3,t4
    real   (dpR):: is_addD_r
    character(20):: vis_slave = 'vis_slave'

    ! call fv_struct_to_array(lev)
    d       =  0.0d0
    mu      =  0.0d0
    do im=1,mesh(lev)%n_mortar_b
        sL  =  mesh(lev)%mortar_LR(1,im)
        eL  =  mesh(lev)%mortar_LR(2,im)
        sR  =  mesh(lev)%mortar_LR(3,im)
        eR  =  mesh(lev)%mortar_LR(4,im)

        idxL = sec(sL)%ID_ele_g(eL)
        idxR = sec(sR)%ID_ele_g(eR)
        if(idxL .le. cellNum) idxL = perm(idxL)
        if(idxR .le. cellNum) idxR = perm(idxR)

        u  (1 :5 )  =  fv_u(1:5,idxR)
        mu (1    )  =  fv_mu(1,idxR)
        if(is_tur_cal)  mu(2)   =  fv_mu  (2,idxR)
        if(is_KO     )  tke     =  fv_turb(1,idxR)
        g  (1 :9 )  =  fv_gra(4 :12,idxR)
        g  (10:12)  =  fv_gra(16:18,idxR)

        call get_vis_fac(1, mesh(lev)%mortar_n_vg(1,im), u, mu, g, tke, rhsD)

        if(is_addD) then
            fv_rhs(1:5,idxL)  =  fv_rhs(1:5,idxL)-rhsD(1:5,1)
        else
            fv_duc(1:5,idxL)  =  fv_duc(1:5,idxL)+rhsD(1:5,1)
        end if
    end do

    ! imi = 1
    ! call sec_struct_to_array(lev)
    ! call transform_parameter_to_real(lev)

    if(is_addD) then
        is_addD_r = 1
    else
        is_addD_r = -1
    end if

if(sw_time) call system_clock(t1)
call starttime(vis_slave)
    call rhs_vis_host(mesh_reordered(lev)%owner, mesh_reordered(lev)%neighbor, &
        & faceNum, cellNum, mesh_reordered(lev)%mortar_transform, &
        & mesh_reordered(lev)%mortar_n_vg, fv_u, sec_cen, fv_gra, fv_rhs, &
        & fv_duc, sec_is_int, fv_mu, fv_turb, fv_t, is_addD_r, is_tur_cal_r, &
        & cp, prL, prT, is_KO_r, is_RANS_r, RANS_model_r)
call endtime(vis_slave)
if(sw_time) call system_clock(t2)

if(sw_time) then
    call system_clock(t3)

    do im=1+mesh(lev)%n_mortar_b,mesh(lev)%n_mortar
        sL  =  mesh(lev)%mortar_LR(1,im)
        eL  =  mesh(lev)%mortar_LR(2,im)
        sR  =  mesh(lev)%mortar_LR(3,im)
        eR  =  mesh(lev)%mortar_LR(4,im)
        fac_1d(1:4,1)   =  mesh(lev)%mortar_n_vg(1:4,im)

        uLf(1:5,1)  =  fv(sL)%u(1:5,eL)
        uRf(1:5,1)  =  fv(sR)%u(1:5,eR)
        u  (1:5  )  =  0.5d0*(uLf(1:5,1)+uRf(1:5,1))
        mu (1    )  = (fv(sL)%mu(1,eL)+fv(sR)%mu(1,eR))*0.5d0
        if(is_tur_cal)  mu(2)   = (fv(sL)%mu  (2,eL)+fv(sR)%mu  (2,eR))*0.5d0
        if(is_KO     )  tke     = (fv(sL)%turb(1,eL)+fv(sR)%turb(1,eR))*0.5d0
        tL          =  fv(sL)%t(    eL)
        tR          =  fv(sR)%t(    eR)
        g  (1 :9 )  = (fv(sL)%gra(4 :12,eL)+fv(sR)%gra(4 :12,eR))*0.5d0
        g  (10:12)  = (fv(sL)%gra(16:18,eL)+fv(sR)%gra(16:18,eR))*0.5d0

        d(1:3)      =  sec(sR)%cen(1:3,eR)-sec(sL)%cen(1:3,eL)
        rtmp        =  1.0d0/sqrt(d(1)**2+d(2)**2+d(3)**2)
        d           =  d*rtmp
        du(1:3)     =  rtmp*(uRf(2:4,1)-uLf(2:4,1))
        du(4  )     =  rtmp*(tR-tL)

        nne(1:3)=  fac_1d(1:3,1)/(fac_1d(1,1)*d(1)+fac_1d(2,1)*d(2)+fac_1d(3,1)*d(3))
        do j=1,4
            k       =  3*j-2
            dr      =  du(j)-g(k)*d(1)-g(k+1)*d(2)-g(k+2)*d(3)
            g(k:k+2)=  g(k:k+2)+dr*nne(1:3)
        end do

        call get_vis_fac(1, fac_1d, u, mu, g, tke, rhsD)

        if(is_addD) then
            fv(sL)%rhs(1:5,eL)  =  fv(sL)%rhs(1:5,eL)-rhsD(1:5,1)
            if(sec(sR)%is_int)  fv(sR)%rhs(1:5,eR)  =  fv(sR)%rhs(1:5,eR)+rhsD(1:5,1)
        else
            fv(sL)%duc(1:5,eL)  =  fv(sL)%duc(1:5,eL)+rhsD(1:5,1)
            if(sec(sR)%is_int)  fv(sR)%duc(1:5,eR)  =  fv(sR)%duc(1:5,eR)-rhsD(1:5,1)
        end if
    end do

    call system_clock(t4)
    write(*,*),'Processor ID: ', myid, ', Speed-up of viscous term: ',real(t4-t3)/(t2-t1)
! write(*,*),'cpu_time: ', end-start
end if
    ! imi = 0
    ! do im=1+mesh(lev)%n_mortar_b,mesh(lev)%n_mortar
    !     imi = imi + 1
    !     if(mesh_reordered(lev)%mortar_transform(imi) .eq. 0) then
    !         idxL = mesh_reordered(lev)%mortar_own_ID(imi)
    !         idxR = mesh_reordered(lev)%mortar_nei_ID(imi)
    !     else 
    !         idxL = mesh_reordered(lev)%mortar_nei_ID(imi)
    !         idxR = mesh_reordered(lev)%mortar_own_ID(imi)
    !     end if
    !     fac_1d(1:4,1) = mesh_reordered(lev)%mortar_n_vg(1:4,imi)

    !     uLf(1:5,1) = fv_u(1:5,idxL)
    !     uRf(1:5,1) = fv_u(1:5,idxR)

    !     u(1:5    ) = 0.5d0*(uLf(1:5,1)+uRf(1:5,1))
    !     mu(1     ) = 0.5d0*(fv_mu(1,idxL)+fv_mu(1,idxR))
    !     if(is_tur_cal) mu(2)  = 0.5d0*(fv_mu(2,idxL)+fv_mu(2,idxR))
    !     if(is_KO)      tke    = 0.5d0*(fv_turb(1,idxL)+fv_turb(1,idxR))
    !     tL         = fv_t(idxL)
    !     tR         = fv_t(idxR)
    !     g(1:9)     = 0.5d0*(fv_gra(4:12,idxL)+fv_gra(4:12,idxR))
    !     g(10:12)   = 0.5d0*(fv_gra(16:18,idxL)+fv_gra(16:18,idxR))

    !     d(1:3)     = sec_cen(1:3,idxR)-sec_cen(1:3,idxL)
    !     rtmp       = 1.0d0/sqrt(d(1)**2+d(2)**2+d(3)**2)
    !     d          = d*rtmp
    !     du(1:3)    = rtmp*(uRf(2:4,1)-uLf(2:4,1))
    !     du(4)      = rtmp*(tR-tL)

    !     nne(1:3)   = fac_1d(1:3,1)/(fac_1d(1,1)*d(1)+fac_1d(2,1)*d(2)+fac_1d(3,1)*d(3))
    !     do j=1,4
    !         k        = 3*j-2
    !         dr       = du(j)-g(k)*d(1)-g(k+1)*d(2)-g(k+2)*d(3)
    !         g(k:k+2) = g(k:k+2)+dr*nne(1:3)
    !     end do
    !     call get_vis_fac_test(1, fac_1d, u, mu, g, tke, rhsD)

    !     if(is_addD) then
    !         fv_rhs(1:5,idxL)  =  fv_rhs(1:5,idxL)-rhsD(1:5,1)
    !         if(sec_is_int(idxR) .eq. 1)  fv_rhs(1:5,idxR)  =  fv_rhs(1:5,idxR)+rhsD(1:5,1)
    !     else
    !         fv_duc(1:5,idxL)  =  fv_duc(1:5,idxL)+rhsD(1:5,1)
    !         if(sec_is_int(idxR) .eq. 1)  fv_duc(1:5,idxR)  =  fv_duc(1:5,idxR)-rhsD(1:5,1)
    !     end if

    ! end do

    ! imi = 0
    ! do im=1+mesh(lev)%n_mortar_b,mesh(lev)%n_mortar
    !     imi = imi + 1
    !     if(mesh_reordered(lev)%mortar_transform(imi) .eq. 0) then
    !         idxL = mesh_reordered(lev)%mortar_own_ID(imi)
    !         idxR = mesh_reordered(lev)%mortar_nei_ID(imi)
    !     else 
    !         idxL = mesh_reordered(lev)%mortar_nei_ID(imi)
    !         idxR = mesh_reordered(lev)%mortar_own_ID(imi)
    !     end if
    !     fac_1d(1:4,1) = mesh_reordered(lev)%mortar_n_vg(1:4,imi)

    !     uLf(1:5,1) = fv_u(1:5,idxL)
    !     uRf(1:5,1) = fv_u(1:5,idxR)

    !     u(1:5    ) = 0.5d0*(uLf(1:5,1)+uRf(1:5,1))
    !     mu(1     ) = 0.5d0*(fv_mu(1,idxL)+fv_mu(1,idxR))
    !     if(is_tur_cal) mu(2)  = 0.5d0*(fv_mu(2,idxL)+fv_mu(2,idxR))
    !     if(is_KO)      tke    = 0.5d0*(fv_turb(1,idxL)+fv_turb(1,idxR))
    !     tL         = fv_t(idxL)
    !     tR         = fv_t(idxR)
    !     g(1:9)     = 0.5d0*(fv_gra(4:12,idxL)+fv_gra(4:12,idxR))
    !     g(10:12)   = 0.5d0*(fv_gra(16:18,idxL)+fv_gra(16:18,idxR))

    !     d(1:3)     = sec_cen(1:3,idxR)-sec_cen(1:3,idxL)
    !     rtmp       = 1.0d0/sqrt(d(1)**2+d(2)**2+d(3)**2)
    !     d          = d*rtmp
    !     du(1:3)    = rtmp*(uRf(2:4,1)-uLf(2:4,1))
    !     du(4)      = rtmp*(tR-tL)

    !     nne(1:3)   = fac_1d(1:3,1)/(fac_1d(1,1)*d(1)+fac_1d(2,1)*d(2)+fac_1d(3,1)*d(3))
    !     do j=1,4
    !         k        = 3*j-2
    !         dr       = du(j)-g(k)*d(1)-g(k+1)*d(2)-g(k+2)*d(3)
    !         g(k:k+2) = g(k:k+2)+dr*nne(1:3)
    !     end do
    !     call get_vis_fac(1, fac_1d, u, mu, g, tke, rhsD)

    !     if(is_addD) then
    !         fv_rhs(1:5,idxL)  =  fv_rhs(1:5,idxL)-rhsD(1:5,1)
    !         if(sec_is_int(idxR) .eq. 1)  fv_rhs(1:5,idxR)  =  fv_rhs(1:5,idxR)+rhsD(1:5,1)
    !     else
    !         fv_duc(1:5,idxL)  =  fv_duc(1:5,idxL)+rhsD(1:5,1)
    !         if(sec_is_int(idxR) .eq. 1)  fv_duc(1:5,idxR)  =  fv_duc(1:5,idxR)-rhsD(1:5,1)
    !     end if

    ! end do
    ! call fv_array_to_struct(lev)

    ! do im=1+mesh(lev)%n_mortar_b,mesh(lev)%n_mortar
    !     sL  =  mesh(lev)%mortar_LR(1,im)
    !     eL  =  mesh(lev)%mortar_LR(2,im)
    !     sR  =  mesh(lev)%mortar_LR(3,im)
    !     eR  =  mesh(lev)%mortar_LR(4,im)
    !     fac_1d(1:4,1)   =  mesh(lev)%mortar_n_vg(1:4,im)

    !     uLf(1:5,1)  =  fv(sL)%u(1:5,eL)
    !     uRf(1:5,1)  =  fv(sR)%u(1:5,eR)
    !     u  (1:5  )  =  0.5d0*(uLf(1:5,1)+uRf(1:5,1))
    !     mu (1    )  = (fv(sL)%mu(1,eL)+fv(sR)%mu(1,eR))*0.5d0
    !     if(is_tur_cal)  mu(2)   = (fv(sL)%mu  (2,eL)+fv(sR)%mu  (2,eR))*0.5d0
    !     if(is_KO     )  tke     = (fv(sL)%turb(1,eL)+fv(sR)%turb(1,eR))*0.5d0
    !     tL          =  fv(sL)%t(    eL)
    !     tR          =  fv(sR)%t(    eR)
    !     g  (1 :9 )  = (fv(sL)%gra(4 :12,eL)+fv(sR)%gra(4 :12,eR))*0.5d0
    !     g  (10:12)  = (fv(sL)%gra(16:18,eL)+fv(sR)%gra(16:18,eR))*0.5d0

    !     d(1:3)      =  sec(sR)%cen(1:3,eR)-sec(sL)%cen(1:3,eL)
    !     rtmp        =  1.0d0/sqrt(d(1)**2+d(2)**2+d(3)**2)
    !     d           =  d*rtmp
    !     du(1:3)     =  rtmp*(uRf(2:4,1)-uLf(2:4,1))
    !     du(4  )     =  rtmp*(tR-tL)

    !     nne(1:3)=  fac_1d(1:3,1)/(fac_1d(1,1)*d(1)+fac_1d(2,1)*d(2)+fac_1d(3,1)*d(3))
    !     do j=1,4
    !         k       =  3*j-2
    !         dr      =  du(j)-g(k)*d(1)-g(k+1)*d(2)-g(k+2)*d(3)
    !         g(k:k+2)=  g(k:k+2)+dr*nne(1:3)
    !     end do

    !     call get_vis_fac(1, fac_1d, u, mu, g, tke, rhsD)

    !     if(is_addD) then
    !         fv(sL)%rhs(1:5,eL)  =  fv(sL)%rhs(1:5,eL)+rhsD(1:5,1)
    !         if(sec(sR)%is_int)  fv(sR)%rhs(1:5,eR)  =  fv(sR)%rhs(1:5,eR)-rhsD(1:5,1)
    !     else
    !         fv(sL)%duc(1:5,eL)  =  fv(sL)%duc(1:5,eL)-rhsD(1:5,1)
    !         if(sec(sR)%is_int)  fv(sR)%duc(1:5,eR)  =  fv(sR)%duc(1:5,eR)+rhsD(1:5,1)
    !     end if
    ! end do

    ! do im=1+mesh(lev)%n_mortar_b,mesh(lev)%n_mortar
    !     sL  =  mesh(lev)%mortar_LR(1,im)
    !     eL  =  mesh(lev)%mortar_LR(2,im)
    !     sR  =  mesh(lev)%mortar_LR(3,im)
    !     eR  =  mesh(lev)%mortar_LR(4,im)
    !     fac_1d(1:4,1)   =  mesh(lev)%mortar_n_vg(1:4,im)

    !     uLf(1:5,1)  =  fv(sL)%u(1:5,eL)
    !     uRf(1:5,1)  =  fv(sR)%u(1:5,eR)
    !     u  (1:5  )  =  0.5d0*(uLf(1:5,1)+uRf(1:5,1))
    !     mu (1    )  = (fv(sL)%mu(1,eL)+fv(sR)%mu(1,eR))*0.5d0
    !     if(is_tur_cal)  mu(2)   = (fv(sL)%mu  (2,eL)+fv(sR)%mu  (2,eR))*0.5d0
    !     if(is_KO     )  tke     = (fv(sL)%turb(1,eL)+fv(sR)%turb(1,eR))*0.5d0
    !     tL          =  fv(sL)%t(    eL)
    !     tR          =  fv(sR)%t(    eR)
    !     g  (1 :9 )  = (fv(sL)%gra(4 :12,eL)+fv(sR)%gra(4 :12,eR))*0.5d0
    !     g  (10:12)  = (fv(sL)%gra(16:18,eL)+fv(sR)%gra(16:18,eR))*0.5d0

    !     d(1:3)      =  sec(sR)%cen(1:3,eR)-sec(sL)%cen(1:3,eL)
    !     rtmp        =  1.0d0/sqrt(d(1)**2+d(2)**2+d(3)**2)
    !     d           =  d*rtmp
    !     du(1:3)     =  rtmp*(uRf(2:4,1)-uLf(2:4,1))
    !     du(4  )     =  rtmp*(tR-tL)

    !     nne(1:3)=  fac_1d(1:3,1)/(fac_1d(1,1)*d(1)+fac_1d(2,1)*d(2)+fac_1d(3,1)*d(3))
    !     do j=1,4
    !         k       =  3*j-2
    !         dr      =  du(j)-g(k)*d(1)-g(k+1)*d(2)-g(k+2)*d(3)
    !         g(k:k+2)=  g(k:k+2)+dr*nne(1:3)
    !     end do

    !     call get_vis_fac(1, fac_1d, u, mu, g, tke, rhsD)

    !     if(is_addD) then
    !         fv(sL)%rhs(1:5,eL)  =  fv(sL)%rhs(1:5,eL)-rhsD(1:5,1)
    !         if(sec(sR)%is_int)  fv(sR)%rhs(1:5,eR)  =  fv(sR)%rhs(1:5,eR)+rhsD(1:5,1)
    !     else
    !         fv(sL)%duc(1:5,eL)  =  fv(sL)%duc(1:5,eL)+rhsD(1:5,1)
    !         if(sec(sR)%is_int)  fv(sR)%duc(1:5,eR)  =  fv(sR)%duc(1:5,eR)-rhsD(1:5,1)
    !     end if
    ! end do

    if(.not. is_cht)    return
    stop 'is_cht has not been implemented on Sunway platform'
    do im=1,mesh(lev)%n_mortar_b
        sL  =  mesh(lev)%mortar_LR(1,im)
        eL  =  mesh(lev)%mortar_LR(2,im)
        sR  =  mesh(lev)%mortar_LR(3,im)
        eR  =  mesh(lev)%mortar_LR(4,im)

        if(sec(sR)%is_cht) then
            fv(sL)%rhs(5,eL)=  fv(sL)%rhs(5,eL)-fv(sR)%heat_flux(eR)
        end if
    end do

    return
    end subroutine fv_get_rhs_vis_sw
!-------------------------------------------------------------------------------
!   compute the source term due to rotation.
!-------------------------------------------------------------------------------
    subroutine fv_get_rhs_src_rotation(lev)
    use var_kind_def
    use var_fv
    use var_global, only: rotation_axis,body_force
    use var_mesh
    implicit none
    integer(dpI),intent(in):: lev
    integer(dpI):: isec,ID,iele
    real   (dpR):: x(3),y(3),r(3),u(3),f(5)

    if(lev .ne. 0)  return
    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(.not. sec(isec)%is_int)  cycle
        ID  =  sec(isec)%ID_sec_g
        if(sec_motion_type(ID) .ne. 2)  cycle

        x   =  sec_motion_speed(ID)*rotation_axis
        do iele=1,sec(isec)%n_ele
            y(1:3)  =  fv(isec)%u(1,iele)*fv(isec)%u(2:4,iele)
            r(1)    =  x(2)*y(3)-x(3)*y(2)
            r(2)    = -x(1)*y(3)+x(3)*y(1)
            r(3)    =  x(1)*y(2)-x(2)*y(1)
            fv(isec)%rhs(2:4,iele)  =  fv(isec)%rhs(2:4,iele)+sec(isec)%vol(iele)*r(1:3)
        end do
    end do

    if(body_force(1) .le. -1.0d10)  return

    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(.not. sec(isec)%is_int)  cycle

        do iele=1,sec(isec)%n_ele
            u(1:3)  =  fv(isec)%u(2:4,iele)
            f(2:4)  =  body_force(1:3)
            f(5  )  =  f(2)*u(1)+f(3)*u(2)+f(4)*u(3)
            fv(isec)%rhs(2:5,iele)  =  fv(isec)%rhs(2:5,iele)-sec(isec)%vol(iele)*f(2:5)
        end do
    end do

    return
    end subroutine fv_get_rhs_src_rotation
!-------------------------------------------------------------------------------
!   compute the source term due to rotation.
!-------------------------------------------------------------------------------
    subroutine fv_get_rhs_src_rotation_sw(lev)
    use var_kind_def
    use var_fv
    use var_global, only: rotation_axis,body_force
    use var_mesh
    implicit none
    integer(dpI),intent(in):: lev
    integer(dpI):: isec,ID,iele
    real   (dpR):: x(3),y(3),r(3),u(3),f(5)

    if(lev .ne. 0)  return
    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(.not. sec(isec)%is_int)  cycle
        ID  =  sec(isec)%ID_sec_g
        if(sec_motion_type(ID) .ne. 2)  cycle

        stop 'src_rotation_sw has not been implemented on Sunway platform'
    end do

    if(body_force(1) .le. -1.0d10)  return

    stop 'body_force has not been implemented on Sunway platform'

    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(.not. sec(isec)%is_int)  cycle

        do iele=1,sec(isec)%n_ele
            u(1:3)  =  fv(isec)%u(2:4,iele)
            f(2:4)  =  body_force(1:3)
            f(5  )  =  f(2)*u(1)+f(3)*u(2)+f(4)*u(3)
            fv(isec)%rhs(2:5,iele)  =  fv(isec)%rhs(2:5,iele)-sec(isec)%vol(iele)*f(2:5)
        end do
    end do

    return
    end subroutine fv_get_rhs_src_rotation_sw
!-------------------------------------------------------------------------------
!   get gradient for the mortar element. Sunway version
!-------------------------------------------------------------------------------
    subroutine fv_mortar_get_gra_sw(lev)
    use var_kind_def
    use var_cgns, only: BCWallViscous
    use var_fv
    use var_fv_array
    use var_mesh
    implicit none
    integer(dpI),intent(in):: lev
    integer(dpI):: im,sL,eL,sR,eR,j,idxL,idxR
    real   (dpR):: n(3),gn,du(6),g(3,6),d

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
        n(1:3)  =  sec(sL)%cen(1:3,eL)-mesh(lev)%mortar_cen(1:3,im)
        call norm_vec(3, n, d)

        call DCOPY(18, fv_gra(1,idxL), 1, g, 1)
        du(1:5) =  fv_u(1:5,idxL)-fv_u(1:5,idxR)
        du(6  ) =  fv_t(    idxL)-fv_t(  idxR)
        do j=1,6
            gn      =  g(1,j)*n(1)+g(2,j)*n(2)+g(3,j)*n(3)
            g(1:3,j)=  g(1:3,j)+(du(j)/d-gn)*n(1:3)
        end do
        call DCOPY(18, g, 1, fv_gra(1,idxR), 1)

        if(sec(sR)%bct .ne. BCWallViscous)  cycle
        gn  =  fv_gra(16,idxR)*mesh(lev)%mortar_n_vg(1,im) &
            & +fv_gra(17,idxR)*mesh(lev)%mortar_n_vg(2,im) &
            & +fv_gra(18,idxR)*mesh(lev)%mortar_n_vg(3,im)
        fv_gra(16:18,idxR) =  fv_gra(16:18,idxR)-gn*mesh(lev)%mortar_n_vg(1:3,im)
    end do
! call fv_array_to_struct(lev)

    return
    end subroutine fv_mortar_get_gra_sw
!-------------------------------------------------------------------------------
!   get gradient for the mortar element.
!-------------------------------------------------------------------------------
    subroutine fv_mortar_get_gra(lev)
    use var_kind_def
    use var_cgns, only: BCWallViscous
    use var_fv
    use var_mesh
    implicit none
    integer(dpI),intent(in):: lev
    integer(dpI):: im,sL,eL,sR,eR,j
    real   (dpR):: n(3),gn,du(6),g(3,6),d

    do im=1,mesh(lev)%n_mortar_b
        sL  =  mesh(lev)%mortar_LR(1,im)
        eL  =  mesh(lev)%mortar_LR(2,im)
        sR  =  mesh(lev)%mortar_LR(3,im)
        eR  =  mesh(lev)%mortar_LR(4,im)
        n(1:3)  =  sec(sL)%cen(1:3,eL)-mesh(lev)%mortar_cen(1:3,im)
        call norm_vec(3, n, d)

        call DCOPY(18, fv(sL)%gra(1,eL), 1, g, 1)
        du(1:5) =  fv(sL)%u(1:5,eL)-fv(sR)%u(1:5,eR)
        du(6  ) =  fv(sL)%t(    eL)-fv(sR)%t(  eR)
        do j=1,6
            gn      =  g(1,j)*n(1)+g(2,j)*n(2)+g(3,j)*n(3)
            g(1:3,j)=  g(1:3,j)+(du(j)/d-gn)*n(1:3)
        end do
        call DCOPY(18, g, 1, fv(sR)%gra(1,eR), 1)

        if(sec(sR)%bct .ne. BCWallViscous)  cycle
        gn  =  fv(sR)%gra(16,eR)*mesh(lev)%mortar_n_vg(1,im) &
            & +fv(sR)%gra(17,eR)*mesh(lev)%mortar_n_vg(2,im) &
            & +fv(sR)%gra(18,eR)*mesh(lev)%mortar_n_vg(3,im)
        fv(sR)%gra(16:18,eR) =  fv(sR)%gra(16:18,eR)-gn*mesh(lev)%mortar_n_vg(1:3,im)
    end do

    return
    end subroutine fv_mortar_get_gra
!-------------------------------------------------------------------------------
!   get the Venkatakrishnan limiter.
!-------------------------------------------------------------------------------
    subroutine fv_get_Venka_limiter(isec,mode)
    use var_kind_def
    use var_air, only: gk
    use var_fv
    use var_global, only: n_dim,pi,rref,pref,L_ref
    use var_mesh
    implicit none
    integer(dpI),intent(in):: isec,mode
    real   (dpR),parameter:: beta=1.0d0,K_venka=1.0d0
    integer(dpI):: e1,e0,iele,ss,ee,i,j,im
    real   (dpR):: phi,u(6),rh_min,rh_max,p_min,p_max,alp,d_rh,d_p, &
                &  v1,v2,d(3),L_c,eps_rh,eps_p,ua,aa,Ma(6),um(5)

    if(mode .eq. 1) then
        e1  =  1
        e0  =  sec(isec)%n_ele_b
    elseif(mode .eq. 2) then
        e1  =  1+sec(isec)%n_ele_b
        e0  =  sec(isec)%n_ele
    else
        e1  =  1
        e0  =  sec(isec)%n_ele
    end if

    do iele=e1,e0
        u(1:5)  =  fv(isec)%u(1:5,iele)
        rh_min  =  u(1)
        rh_max  =  u(1)
        p_min   =  u(5)
        p_max   =  u(5)
        do i=sec(isec)%iA_face_neighbour(iele),sec(isec)%iA_face_neighbour(iele+1)-1
            im  =  sec(isec)%jA_face_neighbour(i)
            if(im .gt. 0) then
                ss  =  mesh(0)%mortar_LR(3, im)
                ee  =  mesh(0)%mortar_LR(4, im)
            else
                ss  =  mesh(0)%mortar_LR(1,-im)
                ee  =  mesh(0)%mortar_LR(2,-im)
            end if

            rh_min  =  min(rh_min, fv(ss)%u(1,ee))
            rh_max  =  max(rh_max, fv(ss)%u(1,ee))
            p_min   =  min(p_min , fv(ss)%u(5,ee))
            p_max   =  max(p_max , fv(ss)%u(5,ee))

            im      =  abs(sec(isec)%jA_face_neighbour(i))
            um(1:5) =  0.5d0*(u(1:5)+fv(ss)%u(1:5,ee))
            ua      = (u(2)-mesh(0)%mortar_n_vg(6,im))**2 &
                    &+(u(3)-mesh(0)%mortar_n_vg(7,im))**2 &
                    &+(u(4)-mesh(0)%mortar_n_vg(8,im))**2
            aa      =  gk*um(5)/um(1)
            Ma(i-sec(isec)%iA_face_neighbour(iele)+1)   =  sqrt(ua/aa)
        end do

        d       =  0.0d0
        phi     =  1.0d0
        L_c     =  sec(isec)%vol(iele)**(1.0d0/real(n_dim, dpR))
        alp     = (K_venka*L_c/L_ref)**3
        eps_rh  =  alp*rref**2
        eps_p   =  alp*pref**2
        do i=sec(isec)%iA_face_neighbour(iele),sec(isec)%iA_face_neighbour(iele+1)-1
            im          =  abs(sec(isec)%jA_face_neighbour(i))
            d(1:n_dim)  =  mesh(0)%mortar_cen(1:n_dim,im)-sec(isec)%cen(1:n_dim,iele)

            d_rh=  fv(isec)%gra(1,iele)*d(1)+fv(isec)%gra(2,iele)*d(2) &
                & +fv(isec)%gra(3,iele)*d(3)
            if(d_rh .ge. 0.0d0) then
                alp =  beta*(rh_max-u(1))
            else
                alp =  beta*(rh_min-u(1))
            end if
            v1  = (alp**2+2.0d0*alp*d_rh+eps_rh)/(alp**2+alp*d_rh+2.0d0*d_rh**2+eps_rh)
            phi =  min(phi, v1)

            d_p =  fv(isec)%gra(13,iele)*d(1)+fv(isec)%gra(14,iele)*d(2) &
                & +fv(isec)%gra(15,iele)*d(3)
            if(d_p .ge. 0.0d0) then
                alp =  beta*(p_max-u(5))
            else
                alp =  beta*(p_min-u(5))
            end if
            v2  = (alp**2+2.0d0*alp*d_p +eps_p )/(alp**2+alp*d_p +2.0d0*d_p **2+eps_p )
            phi =  min(phi, v2)

            j   =  i-sec(isec)%iA_face_neighbour(iele)+1
            phi =  max(phi, 0.5d0*(1.0d0-tanh(5.0d0*pi*(Ma(j)-1.0d0))))
        end do
        fv(isec)%gra(19,iele)   =  phi
    end do

    return
    end subroutine fv_get_Venka_limiter
!-------------------------------------------------------------------------------
!   cal gradient.
!-------------------------------------------------------------------------------
    subroutine fv_get_gra(lev)
    use var_kind_def
    use var_fv
    use var_mesh
    use var_parallel
    use var_slv, only: gradient_method,gradient_GG,gradient_LS
    use var_turb, only: RANS_model,RANS_WA
    use var_global, only: sw_slave
    implicit none
    integer(dpI),intent(in):: lev
    logical(dpL):: is_WA
    integer(dpI):: isec,iele,i,ip_remote,isr,LDA,idx
    ! character(20):: gra_gg_pre = 'gra_gg_pre'
    ! character(20):: gra_gg = 'gra_gg'
    ! character(20):: gra_gg_mortar = 'gra_gg_mortar'

    is_WA   = (RANS_model .eq. RANS_WA) .and. (gradient_method .ne. gradient_GG)
! call starttime(gra_gg_pre)
    if(gradient_method .eq. gradient_GG)  then
        ! if(sw_slave) then
            ! call fv_get_gra_gg_pre_sw(lev)
        ! else
            call fv_get_gra_gg_pre(lev)
        ! end if
    end if
! call endtime(gra_gg_pre)

!   ----------------------------------------------------------------------------
!   cal gradient for the elements involved in the data exchange.
    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(sec(isec)%is_int)    call get_gra(isec, 1)
    end do
!   cal gradient for the elements involved in the data exchange.
!   ----------------------------------------------------------------------------

    LDA =  19

!   ----------------------------------------------------------------------------
!   prepare data.
    if(is_WA) then
        do isr=mesh(lev)%p2p_1,mesh(lev)%p2p_0
            idx =  1
            do i=1,p2p(isr)%n_ele_send
                isec=  p2p(isr)%id_ele_send(1,i)
                iele=  p2p(isr)%id_ele_send(2,i)
                call DCOPY(LDA, fv(isec)%gra(1,iele), 1, p2p(isr)%rsend(idx), 1)
                idx =  idx+LDA
            end do
            p2p(isr)%n_send =  idx-1
            p2p(isr)%n_recv =  LDA*p2p(isr)%n_ele_recv
        end do
    else
        do isr=mesh(lev)%p2p_1,mesh(lev)%p2p_0
            idx =  1
            do i=1,p2p(isr)%n_ele_send_face
                isec=  p2p(isr)%id_ele_send(1,i)
                iele=  p2p(isr)%id_ele_send(2,i)
                call DCOPY(LDA, fv(isec)%gra(1,iele), 1, p2p(isr)%rsend(idx), 1)
                idx =  idx+LDA
            end do
            p2p(isr)%n_send =  idx-1
            p2p(isr)%n_recv =  LDA*p2p(isr)%n_ele_recv_face
        end do
    end if
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
            call mpi_isend(p2p(isr)%rsend, p2p(isr)%n_send, mpi_dpR, ip_remote, &
                &  myid, mpi_comm_world, mpi_req(mpi_nreq), mpi_err)
        end if

        if(p2p(isr)%n_recv .gt. 0) then
            mpi_nreq=  mpi_nreq+1
            call mpi_irecv(p2p(isr)%rrecv, p2p(isr)%n_recv, mpi_dpR, ip_remote, &
                &  ip_remote, mpi_comm_world, mpi_req(mpi_nreq), mpi_err)
        end if
    end do
!   data exchange.
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

! call starttime(gra_gg)
!   ----------------------------------------------------------------------------
!   cal gradient for the elements not involved in the data exchange.
    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(sec(isec)%is_int)    call get_gra(isec, 2)
    end do
!   cal gradient for the elements not involved in the data exchange.
!   ----------------------------------------------------------------------------
! call endtime(gra_gg)

    call mpi_waitall(mpi_nreq, mpi_req, mpi_sta, mpi_err)

!   ----------------------------------------------------------------------------
!   unpack data received.
    if(is_WA) then
        do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
            if(.not. sec(isec)%is_ghost)    cycle
            do i=1,sec(isec)%n_ele
                isr =  sec(isec)%id_recv(1,i)
                iele=  sec(isec)%id_recv(2,i)
                call per_rot_vec(11,sec(isec)%per_path(1,i), &
                    &  p2p(isr)%rrecv(1+LDA*(iele-1)), fv(isec)%gra(1,i))
            end do
        end do
    else
        do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
            if(.not. sec(isec)%is_ghost)    cycle
            do i=1,sec(isec)%n_ele
                isr =  sec(isec)%id_recv(1,i)
                iele=  sec(isec)%id_recv(2,i)
                if(iele .gt. p2p(isr)%n_ele_recv_face)  cycle
                call per_rot_vec(11,sec(isec)%per_path(1,i), &
                    &  p2p(isr)%rrecv(1+LDA*(iele-1)), fv(isec)%gra(1,i))
            end do
        end do
    end if
!   unpack data received.
!   ----------------------------------------------------------------------------

! call starttime(gra_gg_mortar)
    call fv_mortar_get_gra(lev)
! call endtime(gra_gg_mortar)

    return
    contains
!       ------------------------------------------------------------------------
!       get gra.
!       ------------------------------------------------------------------------
        subroutine get_gra(isec,mode)
        use var_kind_def
        implicit none
        integer(dpI),intent(in):: isec,mode

        if(gradient_method .eq. gradient_GG) then
            call fv_get_gra_gg (isec, mode)
        elseif(gradient_method .eq. gradient_LS) then
            call fv_get_gra_LS (isec, mode)
        else
            call fv_get_gra_LSm(isec, mode)
        end if
        if(is_limiter_on)   call  fv_get_Venka_limiter(isec, mode)

        return
        end subroutine get_gra
    end subroutine fv_get_gra

!-------------------------------------------------------------------------------
!   cal gradient. Sunway version.
!-------------------------------------------------------------------------------
    subroutine fv_get_gra_sw(lev)
    use var_kind_def
    use var_fv
    use var_fv_array
    use var_sec_array
    use var_mesh
    use var_parallel
    use var_slv, only: gradient_method,gradient_GG,gradient_LS
    use var_turb, only: RANS_model,RANS_WA
    use var_global, only: sw_slave
    implicit none
    integer(dpI),intent(in):: lev
    logical(dpL):: is_WA
    integer(dpI):: isec,iele,i,ip_remote,isr,LDA,idx,ID
    character(20):: gra_gg_pre = 'gra_gg_pre'
    character(20):: gra_gg = 'gra_gg'
    character(20):: gra_gg_mortar = 'gra_gg_mortar'
    character(20):: gra_gg_mpi = 'gra_gg_mpi'

! call fv_struct_to_array(lev)
call starttime(gra_gg_pre)
    is_WA   = (RANS_model .eq. RANS_WA) .and. (gradient_method .ne. gradient_GG)
    if(gradient_method .eq. gradient_GG)  then
        call fv_get_gra_gg_pre_sw(lev)
    end if
call endtime(gra_gg_pre)
! call fv_array_to_struct(lev)
! call fv_struct_to_array(lev)
!   ----------------------------------------------------------------------------
!   cal gradient for the elements involved in the data exchange.
call starttime(gra_gg)
    call fv_get_gra_gg_sw(lev, 0)
call endtime(gra_gg)
!   cal gradient for the elements involved in the data exchange.
!   ----------------------------------------------------------------------------

call starttime(gra_gg_mpi)
    LDA =  19
!   ----------------------------------------------------------------------------
!   prepare data.
    if(is_WA) then
        do isr=mesh(lev)%p2p_1,mesh(lev)%p2p_0
            idx =  1
            do i=1,p2p(isr)%n_ele_send
                isec=  p2p(isr)%id_ele_send(1,i)
                iele=  p2p(isr)%id_ele_send(2,i)
                ID = sec(isec)%ID_ele_g(iele)
                if(ID .le. cellNum) ID = perm(ID)
                call DCOPY(LDA, fv_gra(1,ID), 1, p2p(isr)%rsend(idx), 1)
                idx =  idx+LDA
            end do
            p2p(isr)%n_send =  idx-1
            p2p(isr)%n_recv =  LDA*p2p(isr)%n_ele_recv
        end do
    else
        do isr=mesh(lev)%p2p_1,mesh(lev)%p2p_0
            idx =  1
            do i=1,p2p(isr)%n_ele_send_face
                isec=  p2p(isr)%id_ele_send(1,i)
                iele=  p2p(isr)%id_ele_send(2,i)
                ID = sec(isec)%ID_ele_g(iele)
                if(ID .lt. cellNum) ID = perm(ID)
                call DCOPY(LDA, fv_gra(1,ID), 1, p2p(isr)%rsend(idx), 1)
                idx =  idx+LDA
            end do
            p2p(isr)%n_send =  idx-1
            p2p(isr)%n_recv =  LDA*p2p(isr)%n_ele_recv_face
        end do
    end if
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
            call mpi_isend(p2p(isr)%rsend, p2p(isr)%n_send, mpi_dpR, ip_remote, &
                &  myid, mpi_comm_world, mpi_req(mpi_nreq), mpi_err)
        end if

        if(p2p(isr)%n_recv .gt. 0) then
            mpi_nreq=  mpi_nreq+1
            call mpi_irecv(p2p(isr)%rrecv, p2p(isr)%n_recv, mpi_dpR, ip_remote, &
                &  ip_remote, mpi_comm_world, mpi_req(mpi_nreq), mpi_err)
        end if
    end do
!   data exchange.
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
    if(is_WA) then
        do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
            if(.not. sec(isec)%is_ghost)    cycle
            do i=1,sec(isec)%n_ele
                isr =  sec(isec)%id_recv(1,i)
                iele=  sec(isec)%id_recv(2,i)
                ID = sec(isec)%ID_ele_g(i)
                if(ID .le. cellNum) ID = perm(ID)
                call per_rot_vec(11,sec(isec)%per_path(1,i), &
                    &  p2p(isr)%rrecv(1+LDA*(iele-1)), fv_gra(1,ID))
            end do
        end do
    else
        do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
            if(.not. sec(isec)%is_ghost)    cycle
            do i=1,sec(isec)%n_ele
                isr =  sec(isec)%id_recv(1,i)
                iele=  sec(isec)%id_recv(2,i)
                ID = sec(isec)%ID_ele_g(i)
                if(ID .le. cellNum) ID = perm(ID)
                if(iele .gt. p2p(isr)%n_ele_recv_face)  cycle
                call per_rot_vec(11,sec(isec)%per_path(1,i), &
                    &  p2p(isr)%rrecv(1+LDA*(iele-1)), fv_gra(1,ID))
            end do
        end do
    end if
!   unpack data received.
!   ----------------------------------------------------------------------------
call endtime(gra_gg_mpi)

call starttime(gra_gg_mortar)
    call fv_mortar_get_gra_sw(lev)
call endtime(gra_gg_mortar)

! call fv_array_to_struct(lev)

    return
    end subroutine fv_get_gra_sw
!-------------------------------------------------------------------------------
!   get gradient using the Green-Gauss method.
!-------------------------------------------------------------------------------
    subroutine fv_get_gra_gg_pre(lev)
    use var_kind_def
    use var_fv
    use var_mesh
    implicit none
    integer(dpI),intent(in):: lev
    integer(dpI):: im,sL,eL,sR,eR,isec,i
    real   (dpR):: um(6),vL,vR,n(3)
    character(20):: gg_pre_kernel = 'gg_pre_kernel'

    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(sec(isec)%is_int)    fv(isec)%gra=  0.0d0
    end do

    do im=1,mesh(lev)%n_mortar_b
        sL  =  mesh(lev)%mortar_LR(1,im)
        eL  =  mesh(lev)%mortar_LR(2,im)
        sR  =  mesh(lev)%mortar_LR(3,im)
        eR  =  mesh(lev)%mortar_LR(4,im)
        n(1:3)  =  mesh(lev)%mortar_n_vg(1:3,im)*mesh(lev)%mortar_n_vg(4,im)

        um(1:5) =  fv(sR)%u(1:5,eR)
        um(6  ) =  fv(sR)%t(    eR)
        do i=1,6
            fv(sL)%gra(3*i-2:3*i,eL)=  fv(sL)%gra(3*i-2:3*i,eL)+um(i)*n(1:3)
        end do
    end do

call starttime(gg_pre_kernel)
    do im=1+mesh(lev)%n_mortar_b,mesh(lev)%n_mortar
        sL  =  mesh(lev)%mortar_LR(1,im)
        eL  =  mesh(lev)%mortar_LR(2,im)
        sR  =  mesh(lev)%mortar_LR(3,im)
        eR  =  mesh(lev)%mortar_LR(4,im)
        n(1:3)  =  mesh(lev)%mortar_n_vg(1:3,im)*mesh(lev)%mortar_n_vg(4,im)

        vL      =  1.0d0/sec(sL)%vol(eL)
        vR      =  1.0d0/sec(sR)%vol(eR)
        um(1:5) = (vL*fv(sL)%u(1:5,eL)+vR*fv(sR)%u(1:5,eR))/(vL+vR)
        um(6  ) = (vL*fv(sL)%t(    eL)+vR*fv(sR)%t(    eR))/(vL+vR)
        do i=1,6
            fv(sL)%gra(3*i-2:3*i,eL)=  fv(sL)%gra(3*i-2:3*i,eL)+um(i)*n(1:3)
        end do
        if(sec(sR)%is_int) then
        do i=1,6
            fv(sR)%gra(3*i-2:3*i,eR)=  fv(sR)%gra(3*i-2:3*i,eR)-um(i)*n(1:3)
        end do
        end if
    end do
call endtime(gg_pre_kernel)

    return
    end subroutine fv_get_gra_gg_pre
!-------------------------------------------------------------------------------
!   get gradient using the Green-Gauss method. Sunway version
!-------------------------------------------------------------------------------
    subroutine fv_get_gra_gg_pre_sw(lev)
    use var_kind_def
    use var_fv
    use var_fv_array
    use var_sec_array
    use var_mesh
    implicit none
    integer(dpI),intent(in):: lev
    integer(dpI):: im,sL,eL,sR,eR,isec,i,idxL,idxR,imi,ID
    real   (dpR):: um(6),vL,vR,n(3)
    character(20):: gra_gg_pre_slave = 'gra_gg_pre_slave'
    character(20):: gra_gg_pre_zero = 'gra_gg_pre_zero'
    character(20):: gra_gg_pre_bnd = 'gra_gg_pre_bnd'

    ! call fv_struct_to_array(lev)
    ! call sec_struct_to_array(lev)

call starttime(gra_gg_pre_zero)
    ! fv_gra = 0.0d0
    call gra_gg_zero_host(tot_ele, fv_gra)
call endtime(gra_gg_pre_zero)

call starttime(gra_gg_pre_bnd)
    do im=1,mesh(lev)%n_mortar_b
        sL  =  mesh(lev)%mortar_LR(1,im)
        eL  =  mesh(lev)%mortar_LR(2,im)
        sR  =  mesh(lev)%mortar_LR(3,im)
        eR  =  mesh(lev)%mortar_LR(4,im)
        idxL = sec(sL)%ID_ele_g(eL)
        idxR = sec(sR)%ID_ele_g(eR)
        if(idxL .le. cellNum) idxL = perm(idxL)
        if(idxR .le. cellNum) idxR = perm(idxR)
        n(1:3)  =  mesh(lev)%mortar_n_vg(1:3,im)*mesh(lev)%mortar_n_vg(4,im)

        um(1:5) =  fv_u(1:5,idxR)
        um(6  ) =  fv_t(    idxR)
        do i=1,6
            fv_gra(3*i-2:3*i,idxL)=  fv_gra(3*i-2:3*i,idxL)+um(i)*n(1:3)
        end do
    end do
call endtime(gra_gg_pre_bnd)

call starttime(gra_gg_pre_slave)
    call gra_gg_pre_host(mesh_reordered(lev)%owner, mesh_reordered(lev)%neighbor, &
        & faceNum, cellNum, mesh_reordered(lev)%mortar_transform, &
        & mesh_reordered(lev)%mortar_n_vg, fv_u, fv_gra, sec_vol, &
        & fv_t, sec_is_int)
call endtime(gra_gg_pre_slave)

    ! imi = 0
    ! do im=1+mesh(lev)%n_mortar_b,mesh(lev)%n_mortar
    !     imi = imi + 1
    !     if(mesh_reordered(lev)%mortar_transform(imi) .eq. 0) then
    !         idxL = mesh_reordered(lev)%mortar_own_ID(imi)
    !         idxR = mesh_reordered(lev)%mortar_nei_ID(imi)
    !     else
    !         idxL = mesh_reordered(lev)%mortar_nei_ID(imi)
    !         idxR = mesh_reordered(lev)%mortar_own_ID(imi)
    !     end if
    !     n(1:3)  =  mesh_reordered(lev)%mortar_n_vg(1:3,imi)*mesh_reordered(lev)%mortar_n_vg(4,imi)

    !     vL      =  1.0d0/sec_vol(idxL)
    !     vR      =  1.0d0/sec_vol(idxR)
    !     um(1:5) = (vL*fv_u(1:5,idxL)+vR*fv_u(1:5,idxR))/(vL+vR)
    !     um(6  ) = (vL*fv_t(    idxL)+vR*fv_t(    idxR))/(vL+vR)
    !     ! fv_gra(1,idxL) = fv_gra(1,idxL) - um(6)
    !     ! fv_gra(1,idxR) = fv_gra(1,idxR) - um(1)
    !     do i=1,6
    !         fv_gra(3*i-2:3*i,idxL)=  fv_gra(3*i-2:3*i,idxL)-um(i)*n(1:3)
    !     end do
    !     if(sec_is_int(idxR) .eq. 1) then
    !         do i=1,6
    !             fv_gra(3*i-2:3*i,idxR)=  fv_gra(3*i-2:3*i,idxR)+um(i)*n(1:3)
    !         end do
    !     end if
    ! end do

    ! call fv_array_to_struct(lev)
! 
    
    ! do im=1+mesh(lev)%n_mortar_b,mesh(lev)%n_mortar
    !     sL  =  mesh(lev)%mortar_LR(1,im)
    !     eL  =  mesh(lev)%mortar_LR(2,im)
    !     sR  =  mesh(lev)%mortar_LR(3,im)
    !     eR  =  mesh(lev)%mortar_LR(4,im)
    !     n(1:3)  =  mesh(lev)%mortar_n_vg(1:3,im)*mesh(lev)%mortar_n_vg(4,im)

    !     vL      =  1.0d0/sec(sL)%vol(eL)
    !     vR      =  1.0d0/sec(sR)%vol(eR)
    !     um(1:5) = (vL*fv(sL)%u(1:5,eL)+vR*fv(sR)%u(1:5,eR))/(vL+vR)
    !     um(6  ) = (vL*fv(sL)%t(    eL)+vR*fv(sR)%t(    eR))/(vL+vR)
    !     do i=1,6
    !         fv(sL)%gra(3*i-2:3*i,eL)=  fv(sL)%gra(3*i-2:3*i,eL)+um(i)*n(1:3)
    !     end do
    !     if(sec(sR)%is_int) then
    !     do i=1,6
    !         fv(sR)%gra(3*i-2:3*i,eR)=  fv(sR)%gra(3*i-2:3*i,eR)-um(i)*n(1:3)
    !     end do
    !     end if
    ! end do

    return
    end subroutine fv_get_gra_gg_pre_sw
!-------------------------------------------------------------------------------
!   get gradient using the Green-Gauss method.
!-------------------------------------------------------------------------------
    subroutine fv_get_gra_gg(isec,mode)
    use var_kind_def
    use var_fv
    use var_mesh
    implicit none
    integer(dpI),intent(in):: isec,mode
    integer(dpI):: iele,e1,e0

    if(mode .eq. 1) then
        e1  =  1
        e0  =  sec(isec)%n_ele_b
    elseif(mode .eq. 2) then
        e1  =  1+sec(isec)%n_ele_b
        e0  =  sec(isec)%n_ele
    else
        e1  =  1
        e0  =  sec(isec)%n_ele
    end if

    do iele=e1,e0
        fv(isec)%gra(1:18,iele) =  fv(isec)%gra(1:18,iele)/sec(isec)%vol(iele)
    end do

    return
    end subroutine fv_get_gra_gg
!-------------------------------------------------------------------------------
!   get gradient using the Green-Gauss method. Sunway version
!-------------------------------------------------------------------------------
    subroutine fv_get_gra_gg_sw(lev,mode)
    use var_kind_def
    use var_fv
    use var_fv_array
    use var_sec_array
    use var_mesh
    use var_parallel, only: myid
    implicit none
    integer(dpI),intent(in):: mode,lev
    integer(dpI):: iele,e1,e0
    real   (dpR),allocatable:: tmp(:,:)

    call gra_gg_host(tot_ele, fv_gra, sec_vol)

    return
    end subroutine fv_get_gra_gg_sw
!-------------------------------------------------------------------------------
!   get gradient using the Least-square method.
!-------------------------------------------------------------------------------
    subroutine fv_get_gra_LS(isec,mode)
    use var_kind_def
    use var_fv
    use var_mesh
    implicit none
    integer(dpI),intent(in):: isec,mode
    integer(dpI):: e1,e0,ss,ee,iele,j,k
    real   (dpR):: u(6),du(6)

    if(mode .eq. 1) then
        e1  =  1
        e0  =  sec(isec)%n_ele_b
    elseif(mode .eq. 2) then
        e1  =  1+sec(isec)%n_ele_b
        e0  =  sec(isec)%n_ele
    else
        e1  =  1
        e0  =  sec(isec)%n_ele
    end if

    do iele=e1,e0
        u(1:5)  =  fv(isec)%u(1:5,iele)
        u(6  )  =  fv(isec)%t(    iele)
        fv(isec)%gra(1:18,iele) =  0.0d0
        do j=sec(isec)%iA_ls(iele),sec(isec)%iA_ls(iele+1)-1
            ss  =  sec(isec)%jA_ls(1,j)
            ee  =  sec(isec)%jA_ls(2,j)

            du(1:5) =  fv(ss)%u(1:5,ee)-u(1:5)
            du(6  ) =  fv(ss)%t(    ee)-u(6  )
            do k=1,6
                fv(isec)%gra(3*k-2:3*k,iele)=  fv(isec)%gra(3*k-2:3*k,iele) &
                    & +du(k)*sec(isec)%coe_ls(1:3,j)
            end do
        end do
    end do

    return
    end subroutine fv_get_gra_LS
!-------------------------------------------------------------------------------
!   get gradient using the Least-square method, memory efficient.
!-------------------------------------------------------------------------------
    subroutine fv_get_gra_LSm(isec,mode)
    use var_kind_def
    use var_fv
    use var_mesh
    use var_slv, only: ls_weight_order
    implicit none
    integer(dpI),intent(in):: isec,mode
    integer(dpI):: e1,e0,ss,ee,iele,j,k
    real   (dpR):: u(6),du(6),dp(3),d,c(3),LHS(3,3)

    if(mode .eq. 1) then
        e1  =  1
        e0  =  sec(isec)%n_ele_b
    elseif(mode .eq. 2) then
        e1  =  1+sec(isec)%n_ele_b
        e0  =  sec(isec)%n_ele
    else
        e1  =  1
        e0  =  sec(isec)%n_ele
    end if

    do iele=e1,e0
        LHS(1,1)=  sec(isec)%coe_ls(1,iele)
        LHS(2,1)=  sec(isec)%coe_ls(2,iele)
        LHS(3,1)=  sec(isec)%coe_ls(3,iele)
        LHS(1,2)=  LHS(2,1)
        LHS(2,2)=  sec(isec)%coe_ls(4,iele)
        LHS(3,2)=  sec(isec)%coe_ls(5,iele)
        LHS(1,3)=  LHS(3,1)
        LHS(2,3)=  LHS(3,2)
        LHS(3,3)=  sec(isec)%coe_ls(6,iele)

        u(1:5)  =  fv(isec)%u(1:5,iele)
        u(6  )  =  fv(isec)%t(    iele)
        fv(isec)%gra(1:18,iele) =  0.0d0
        do j=sec(isec)%iA_ls(iele),sec(isec)%iA_ls(iele+1)-1
            ss  =  sec(isec)%jA_ls(1,j)
            ee  =  sec(isec)%jA_ls(2,j)

            dp(1:3) =  sec(ss)%cen(1:3,ee)-sec(isec)%cen(1:3,iele)
            d       =  1.0d0/sqrt(dp(1)*dp(1)+dp(2)*dp(2)+dp(3)*dp(3))**ls_weight_order
            c(1:3)  =  d*(LHS(1:3,1)*dp(1)+LHS(1:3,2)*dp(2)+LHS(1:3,3)*dp(3))

            du(1:5) =  fv(ss)%u(1:5,ee)-u(1:5)
            du(6  ) =  fv(ss)%t(    ee)-u(6  )
            do k=1,6
                fv(isec)%gra(3*k-2:3*k,iele)=  fv(isec)%gra(3*k-2:3*k,iele)+du(k)*c(1:3)
            end do
        end do
    end do

    return
    end subroutine fv_get_gra_LSm
