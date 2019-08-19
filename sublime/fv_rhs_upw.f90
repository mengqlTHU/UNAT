!-------------------------------------------------------------------------------
!   cal Rhs, convective part, boundary.
!-------------------------------------------------------------------------------
    subroutine fv_get_rhs_boundary(lev,calD,addD)
    use var_kind_def
    use var_cgns
    use var_fv
    use var_mesh
    use var_temp, only: uLf,uRf,vgf,fac_1d,rhsl,rhsD
    implicit none
    logical(dpL),intent(in):: calD,addD
    integer(dpI),intent(in):: lev
    logical(dpL):: savD,is_upwind
    integer(dpI):: im,sL,eL,sR,eR,bct

    savD=  calD .and. (.not. addD)
    do im=1,mesh(lev)%n_mortar_b
        sL  =  mesh(lev)%mortar_LR(1,im)
        eL  =  mesh(lev)%mortar_LR(2,im)
        sR  =  mesh(lev)%mortar_LR(3,im)
        eR  =  mesh(lev)%mortar_LR(4,im)
        fac_1d(1:4,1)   =  mesh(lev)%mortar_n_vg(1:4,im)
        vgf   (1:3,1)   =  mesh(lev)%mortar_n_vg(6:8,im)
        bct =  sec(sR)%bct

        is_upwind   = (bct .eq. BCSymmetryPlane)
        if(is_upwind) then
            uLf(1:5,1)  =  fv(sR)%uL(1:5,eR)
            uRf(1:5,1)  =  fv(sR)%uR(1:5,eR)
            if(rhs_lev(lev) .eq. 3) then
                call rhs_conv_hllc(calD, addD, 1)
            elseif(rhs_lev(lev) .eq. 4) then
                call rhs_conv_ausm(calD, addD, 1)
            else
                call rhs_conv_roe (calD, addD, 1)
            end if
        else
            call u_to_F(1, fac_1d(1:3,1), fv(sR)%u(1,eR), vgf, rhsl)
            rhsl(1:5,1) =  rhsl(1:5,1)*fac_1d(4,1)
            rhsD(1:5,1) =  0.0d0
        end if

        fv(sL)%rhs(1:5,eL)  =  fv(sL)%rhs(1:5,eL)+rhsl(1:5,1)
        if(savD)    fv(sL)%duc(1:5,eL)  =  fv(sL)%duc(1:5,eL)+rhsD(1:5,1)
    end do

    return
    end subroutine fv_get_rhs_boundary
!-------------------------------------------------------------------------------
!   3rd order WENO reconstruction, characteristic based.
!-------------------------------------------------------------------------------
    subroutine recons_3rd_c(n,vg,uLf,uRf)
    use var_kind_def
    use var_air, only: gk
    use var_global, only: rref,R16
    use var_temp, only: u_1d=>u_structured
    implicit none
    real   (dpR),parameter:: eps_WENO=1.0d-4
    real   (dpR),parameter:: eps_vis =1.0d-1
    integer(dpI):: j,K
    real   (dpR),intent(in):: n(*),vg(*)
    real   (dpR):: uLf(*),uRf(*),a,R(5,3),L(3,5),c(3,4),F(3),dU(3),dC(3),rh,S, &
                &  ua,Ma,Re,wC(3),wU(3),a_rh,rh_a,E(3),eps(3),Amp,one(3),u(5,4)

    one     =  1.0d0
    E       =  rref*rref
    rh      =  0.5d0*(u_1d(1,1)+u_1d(1,2))
    F(1:3)  =  0.5d0*(u_1d(2:4,1)+u_1d(2:4,2))-vg(1:3)
    ua      =  sqrt(F(1)**2+F(2)**2+F(3)**2)
    a       =  sqrt(gk*(u_1d(5,1)+u_1d(5,2))/(u_1d(1,1)+u_1d(1,2)))
    Ma      =  ua/a
    Re      =  1.0d20

    if((Ma .le. 0.25d0) .or. (Re .le. 5.0d0)) then
        uLf(1:5)= (     -u_1d(1:5,0)+5.0d0*u_1d(1:5,1)+2.0d0*u_1d(1:5,2))*R16
        uRf(1:5)= (2.0d0*u_1d(1:5,1)+5.0d0*u_1d(1:5,2)-      u_1d(1:5,3))*R16
        return
    elseif(Ma .ge. 0.55d0) then
        Amp =  1.0d0
    else
        Amp =  0.5d0*(6.0d0-4.0d0*dtanh(2.0d1*(Ma-0.4d0)))
    end if
    eps =  E*(Amp*Amp*eps_WENO+eps_vis/Re)

    S       =  a*a
    a_rh    =  a/rh
    rh_a    =  rh/a
    Amp     =  1.0d0/S
    R(1:5,1)= (/1.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0/)
    R(1:5,2)= (/1.0d0, a_rh*n(1), a_rh*n(2), a_rh*n(3), S/)
    R(1:5,3)= (/1.0d0,-a_rh*n(1),-a_rh*n(2),-a_rh*n(3), S/)
    L(1,1:5)= (/1.0d0, 0.0d0, 0.0d0, 0.0d0, -Amp/)
    L(2,1:5)= (/0.0d0, rh_a*n(1), rh_a*n(2), rh_a*n(3), Amp/)*0.5d0
    L(3,1:5)= (/0.0d0,-rh_a*n(1),-rh_a*n(2),-rh_a*n(3), Amp/)*0.5d0

    do j=0,3
        K   =  j+1
        c(1:3,K)=  L(1:3,1)*u_1d(1,j)+L(1:3,2)*u_1d(2,j)+L(1:3,3)*u_1d(3,j) &
                & +L(1:3,4)*u_1d(4,j)+L(1:3,5)*u_1d(5,j)
        u(1:5,K)=  u_1d(1:5,j)-R(1:5,1)*c(1,K)-R(1:5,2)*c(2,K)-R(1:5,3)*c(3,K)
    end do
    uLf(1:5)= (     -u(1:5,1)+5.0d0*u(1:5,2)+2.0d0*u(1:5,3))*R16
    uRf(1:5)= (2.0d0*u(1:5,2)+5.0d0*u(1:5,3)-      u(1:5,4))*R16

    dU(1:3) =  c(1:3,2)-c(1:3,1)
    dC(1:3) =  c(1:3,3)-c(1:3,2)
    F       = (dC-dU)**2
    wC      =  2.0d0*(one+F/(dC*dC+eps))
    wU      =         one+F/(dU*dU+eps)
    F(1:3)  =  0.5d0*(wC*(c(1:3,2)+c(1:3,3))+wU*(3.0d0*c(1:3,2)-c(1:3,1)))/(wC+wU)
    uLf(1:5)=  uLf(1:5)+R(1:5,1)*F(1)+R(1:5,2)*F(2)+R(1:5,3)*F(3)

    dC(1:3) =  c(1:3,3)-c(1:3,2)
    dU(1:3) =  c(1:3,4)-c(1:3,3)
    F       = (dC-dU)**2
    wC      =  2.0d0*(one+F/(dC*dC+eps))
    wU      =         one+F/(dU*dU+eps)
    F(1:3)  =  0.5d0*(wC*(c(1:3,2)+c(1:3,3))+wU*(3.0d0*c(1:3,3)-c(1:3,4)))/(wC+wU)
    uRf(1:5)=  uRf(1:5)+R(1:5,1)*F(1)+R(1:5,2)*F(2)+R(1:5,3)*F(3)

    return
    end subroutine recons_3rd_c
!-------------------------------------------------------------------------------
!   cal Rhs, convective part, fv.
!-------------------------------------------------------------------------------
    subroutine fv_get_rhs_upw(lev,calD,addD)
    use var_kind_def
    use var_fv
    use var_global, only: n_dim,R16,R23,sw_slave
    use var_mesh
    use var_parallel
    use var_temp, only: uLf,uRf,vgf,fac_1d,rhsl,rhsD,u=>u_structured
    use var_turb, only: low_dissipation
    use var_uns_cal, only: is_uns_cal_now
    implicit none
    logical(dpL),intent(in):: calD,addD
    integer(dpI),intent(in):: lev
    logical(dpL):: savD
    integer(dpI):: im,s_LL,e_LL,sL,eL,sR,eR,s_RR,e_RR,j
    real   (dpR):: dL(3),dR(3),limiter_L,limiter_R,duL(5),duR(5),du(5)
    character(20):: upw_kernel='upw_kernel'

    if(is_uns_cal_now) then
        if(low_dissipation .eq. 1) then
            call fv_get_rhs_upw_LD    (lev, calD, addD)
        elseif(low_dissipation .eq. 2) then
            call fv_get_rhs_upw_LDmesh(lev, calD, addD)
        end if
        return
    end if

    dL          =  0.0d0
    dR          =  0.0d0
    limiter_L   =  1.0d0
    limiter_R   =  1.0d0
    savD=  calD .and. (.not. addD)

    call fv_get_rhs_boundary(lev, calD, addD)

    do im=1+mesh(lev)%n_mortar_b,mesh(lev)%n_mortar_ss
        s_LL=  mesh(lev)%mortar_structured_stencil(1,im)
        e_LL=  mesh(lev)%mortar_structured_stencil(2,im)
        sL  =  mesh(lev)%mortar_structured_stencil(3,im)
        eL  =  mesh(lev)%mortar_structured_stencil(4,im)
        sR  =  mesh(lev)%mortar_structured_stencil(5,im)
        eR  =  mesh(lev)%mortar_structured_stencil(6,im)
        s_RR=  mesh(lev)%mortar_structured_stencil(7,im)
        e_RR=  mesh(lev)%mortar_structured_stencil(8,im)
        fac_1d(1:4,1)   =  mesh(lev)%mortar_n_vg(1:4,im)
        vgf   (1:3,1)   =  mesh(lev)%mortar_n_vg(6:8,im)

        u(1:5,0)=  fv(s_LL)%u(1:5,e_LL)
        u(1:5,1)=  fv(sL  )%u(1:5,eL  )
        u(1:5,2)=  fv(sR  )%u(1:5,eR  )
        u(1:5,3)=  fv(s_RR)%u(1:5,e_RR)
        if(is_limiter_on) then
            call recons_3rd_c(fac_1d, vgf, uLf, uRf)
        else
            if(fv(sL)%order .le. 1) then
                uLf(1:5,1)  =  u(1:5,1)
                uRf(1:5,1)  =  u(1:5,2)
            else
                uLf(1:5,1)  =(     -u(1:5,0)+5.0d0*u(1:5,1)+2.0d0*u(1:5,2))*R16
                uRf(1:5,1)  =(2.0d0*u(1:5,1)+5.0d0*u(1:5,2)      -u(1:5,3))*R16
            end if
        end if

        if(rhs_lev(lev) .eq. 2) then
            call rhs_conv_roe (calD, addD, 1)
        elseif(rhs_lev(lev) .eq. 3) then
            call rhs_conv_hllc(calD, addD, 1)
        else
            call rhs_conv_ausm(calD, addD, 1)
        end if

        fv(sL)%rhs(1:5,eL)  =  fv(sL)%rhs(1:5,eL)+rhsl(1:5,1)
        if(sec(sR)%is_int)  fv(sR)%rhs(1:5,eR)  =  fv(sR)%rhs(1:5,eR)-rhsl(1:5,1)
        if(savD) then
            fv(sL)%duc(1:5,eL)  =  fv(sL)%duc(1:5,eL)+rhsD(1:5,1)
            if(sec(sR)%is_int)  fv(sR)%duc(1:5,eR)  =  fv(sR)%duc(1:5,eR)-rhsD(1:5,1)
        end if
    end do
call starttime(upw_kernel)
    do im=1+mesh(lev)%n_mortar_ss,mesh(lev)%n_mortar
        sL  =  mesh(lev)%mortar_LR(1,im)
        eL  =  mesh(lev)%mortar_LR(2,im)
        sR  =  mesh(lev)%mortar_LR(3,im)
        eR  =  mesh(lev)%mortar_LR(4,im)
        fac_1d(1:4,1)   =  mesh(lev)%mortar_n_vg(1:4,im)
        vgf   (1:3,1)   =  mesh(lev)%mortar_n_vg(6:8,im)

        uLf(1:5,1)  =  fv(sL)%u(1:5,eL)
        uRf(1:5,1)  =  fv(sR)%u(1:5,eR)
        if(fv(sL)%order .gt. 1) then
            dL(1:n_dim) =  mesh(lev)%mortar_cen(1:n_dim,im)-sec(sL)%cen(1:n_dim,eL)
            dR(1:n_dim) =  mesh(lev)%mortar_cen(1:n_dim,im)-sec(sR)%cen(1:n_dim,eR)
            do j=1,5
                duL(j)  =  fv(sL)%gra(3*j-2,eL)*dL(1)+fv(sL)%gra(3*j-1,eL)*dL(2) &
                        & +fv(sL)%gra(3*j  ,eL)*dL(3)
                duR(j)  =  fv(sR)%gra(3*j-2,eR)*dR(1)+fv(sR)%gra(3*j-1,eR)*dR(2) &
                        & +fv(sR)%gra(3*j  ,eR)*dR(3)
            end do
            if(is_limiter_on) then
                limiter_L   =  fv(sL)%gra(19,eL)
                limiter_R   =  fv(sR)%gra(19,eR)
            end if
            if(fv(sL)%order .eq. 2) then
                uLf(1:5,1)  =  uLf(1:5,1)+limiter_L*duL
                uRf(1:5,1)  =  uRf(1:5,1)+limiter_R*duR
            else
                du (1:5  )  =  uRf(1:5,1)-uLf(1:5,1)
                uLf(1:5,1)  =  uLf(1:5,1)+R23*duL(1:5)+R16*du(1:5)
                uRf(1:5,1)  =  uRf(1:5,1)+R23*duR(1:5)-R16*du(1:5)
            end if

            if((uLf(1,1) .le. 0.0d0) .or. (uLf(5,1) .le. 0.0d0) .or. &
             & (uRf(1,1) .le. 0.0d0) .or. (uRf(5,1) .le. 0.0d0)) then
                uLf(1:5,1)  =  fv(sL)%u(1:5,eL)
                uRf(1:5,1)  =  fv(sR)%u(1:5,eR)
            end if
        end if

        if(rhs_lev(lev) .eq. 2) then
            call rhs_conv_roe (calD, addD, 1)
        elseif(rhs_lev(lev) .eq. 3) then
            call rhs_conv_hllc(calD, addD, 1)
        elseif(rhs_lev(lev) .eq. 4) then
            call rhs_conv_ausm(calD, addD, 1)
        else
            call rhs_conv_lf  (calD, addD, 1)
        end if

        fv(sL)%rhs(1:5,eL)  =  fv(sL)%rhs(1:5,eL)+rhsl(1:5,1)
        if(sec(sR)%is_int)  fv(sR)%rhs(1:5,eR)  =  fv(sR)%rhs(1:5,eR)-rhsl(1:5,1)
        if(savD) then
            fv(sL)%duc(1:5,eL)  =  fv(sL)%duc(1:5,eL)+rhsD(1:5,1)
            if(sec(sR)%is_int)  fv(sR)%duc(1:5,eR)  =  fv(sR)%duc(1:5,eR)-rhsD(1:5,1)
        end if
    end do
call endtime(upw_kernel)

! if(sw_slave) then
    ! call fv_struct_to_array(lev)
    ! call check_struct_to_array(lev)
    ! call fv_array_to_struct(lev)
! end if
    return
    end subroutine fv_get_rhs_upw
!-------------------------------------------------------------------------------
!   cal Rhs, convective part, fv, low-dissipation.
!-------------------------------------------------------------------------------
    subroutine fv_get_rhs_upw_LD(lev,calD,addD)
    use var_kind_def
    use var_fv
    use var_global, only: n_dim,R16,R23
    use var_mesh
    use var_parallel
    use var_temp, only: uLf,uRf,vgf,fac_1d,rhsl,rhsD,u=>u_structured
    use var_turb, only: is_DDES_now,is_LES_now
    implicit none
    logical(dpL),intent(in):: calD,addD
    integer(dpI),intent(in):: lev
    logical(dpL):: savD
    integer(dpI):: im,s_LL,e_LL,sL,eL,sR,eR,s_RR,e_RR,j
    real   (dpR):: dL(3),dR(3),limiter_L,limiter_R,gls,g(9),nu,nut,LD(100),rhL,rhR, &
                &  duL(5),duR(5),du(5),v(3)

    dL          =  0.0d0
    dR          =  0.0d0
    limiter_L   =  1.0d0
    limiter_R   =  1.0d0
    LD          =  1.0d0
    savD=  calD .and. (.not. addD)

    call fv_get_rhs_boundary(lev, calD, addD)

    do im=1+mesh(lev)%n_mortar_b,mesh(lev)%n_mortar_ss
        s_LL=  mesh(lev)%mortar_structured_stencil(1,im)
        e_LL=  mesh(lev)%mortar_structured_stencil(2,im)
        sL  =  mesh(lev)%mortar_structured_stencil(3,im)
        eL  =  mesh(lev)%mortar_structured_stencil(4,im)
        sR  =  mesh(lev)%mortar_structured_stencil(5,im)
        eR  =  mesh(lev)%mortar_structured_stencil(6,im)
        s_RR=  mesh(lev)%mortar_structured_stencil(7,im)
        e_RR=  mesh(lev)%mortar_structured_stencil(8,im)
        fac_1d(1:4,1)   =  mesh(lev)%mortar_n_vg(1:4,im)
        vgf   (1:3,1)   =  mesh(lev)%mortar_n_vg(6:8,im)
        gls             =  1.0d0/mesh(lev)%mortar_n_vg(5,im)

        u(1:5,0)=  fv(s_LL)%u(1:5,e_LL)
        u(1:5,1)=  fv(sL  )%u(1:5,eL  )
        u(1:5,2)=  fv(sR  )%u(1:5,eR  )
        u(1:5,3)=  fv(s_RR)%u(1:5,e_RR)
        v(1:3)  =  0.5d0*(u(2:4,1)+u(2:4,2))
        if(is_limiter_on) then
            call recons_3rd_c(fac_1d, vgf, uLf, uRf)
        else
            if(fv(sL)%order .eq. 1) then
                uLf(1:5,1)  =  u(1:5,1)
                uRf(1:5,1)  =  u(1:5,2)
            else
                uLf(1:5,1)  =(     -u(1:5,0)+5.0d0*u(1:5,1)+2.0d0*u(1:5,2))*R16
                uRf(1:5,1)  =(2.0d0*u(1:5,1)+5.0d0*u(1:5,2)      -u(1:5,3))*R16
            end if
        end if

        rhL =  1.0d0/fv(sL)%u(1,eL)
        rhR =  1.0d0/fv(sR)%u(1,eR)
        nu  =  0.5d0*(fv(sL)%mu(1,eL)*rhL+fv(sR)%mu(1,eR)*rhR)
        if(is_LES_now .or. is_DDES_now) then
            nut =  0.5d0*(fv(sL)%mu(2,eL)*rhL+fv(sR)%mu(2,eR)*rhR)
        else
            nut =  0.0d0
        end if
        g(1:9)  =  0.5d0*(fv(sL)%gra(4:12,eL)+fv(sR)%gra(4:12,eR))
        call fv_get_blending_vortex(gls,mesh(lev)%mortar_n_vg(5,im), v,g,nu,nut, LD(1))
        if(rhs_lev(lev) .eq. 2) then
            call rhs_conv_roe_ld (calD, addD, 1, LD(1))
        else
            call rhs_conv_AUSM_ld(calD, addD, 1, LD(1))
        end if

        fv(sL)%rhs(1:5,eL)  =  fv(sL)%rhs(1:5,eL)+rhsl(1:5,1)
        if(sec(sR)%is_int)  fv(sR)%rhs(1:5,eR)  =  fv(sR)%rhs(1:5,eR)-rhsl(1:5,1)
        if(savD) then
            fv(sL)%duc(1:5,eL)  =  fv(sL)%duc(1:5,eL)+rhsD(1:5,1)
            if(sec(sR)%is_int)  fv(sR)%duc(1:5,eR)  =  fv(sR)%duc(1:5,eR)-rhsD(1:5,1)
        end if
    end do

    do im=1+mesh(lev)%n_mortar_ss,mesh(lev)%n_mortar
        sL  =  mesh(lev)%mortar_LR(1,im)
        eL  =  mesh(lev)%mortar_LR(2,im)
        sR  =  mesh(lev)%mortar_LR(3,im)
        eR  =  mesh(lev)%mortar_LR(4,im)
        fac_1d(1:4,1)   =  mesh(lev)%mortar_n_vg(1:4,im)
        vgf   (1:3,1)   =  mesh(lev)%mortar_n_vg(6:8,im)
        gls             =  1.0d0/mesh(lev)%mortar_n_vg(5,im)

        uLf(1:5,1)  =  fv(sL)%u(1:5,eL)
        uRf(1:5,1)  =  fv(sR)%u(1:5,eR)
        v  (1:3  )  =  0.5d0*(uLf(2:4,1)+uRf(2:4,1))
        if(fv(sL)%order .gt. 1) then
            dL(1:n_dim) =  mesh(lev)%mortar_cen(1:n_dim,im)-sec(sL)%cen(1:n_dim,eL)
            dR(1:n_dim) =  mesh(lev)%mortar_cen(1:n_dim,im)-sec(sR)%cen(1:n_dim,eR)
            do j=1,5
                duL(j)  =  fv(sL)%gra(3*j-2,eL)*dL(1)+fv(sL)%gra(3*j-1,eL)*dL(2) &
                        & +fv(sL)%gra(3*j  ,eL)*dL(3)
                duR(j)  =  fv(sR)%gra(3*j-2,eR)*dR(1)+fv(sR)%gra(3*j-1,eR)*dR(2) &
                        & +fv(sR)%gra(3*j  ,eR)*dR(3)
            end do
            if(is_limiter_on) then
                limiter_L   =  fv(sL)%gra(19,eL)
                limiter_R   =  fv(sR)%gra(19,eR)
            end if
            if(fv(sL)%order .eq. 2) then
                uLf(1:5,1)  =  uLf(1:5,1)+limiter_L*duL
                uRf(1:5,1)  =  uRf(1:5,1)+limiter_R*duR
            else
                du (1:5  )  =  uRf(1:5,1)-uLf(1:5,1)
                uLf(1:5,1)  =  uLf(1:5,1)+R23*duL(1:5)+R16*du(1:5)
                uRf(1:5,1)  =  uRf(1:5,1)+R23*duR(1:5)-R16*du(1:5)
            end if

            if((uLf(1,1) .le. 0.0d0) .or. (uLf(5,1) .le. 0.0d0) .or. &
             & (uRf(1,1) .le. 0.0d0) .or. (uRf(5,1) .le. 0.0d0)) then
                uLf(1:5,1)  =  fv(sL)%u(1:5,eL)
                uRf(1:5,1)  =  fv(sR)%u(1:5,eR)
            end if
        end if

        rhL =  1.0d0/fv(sL)%u(1,eL)
        rhR =  1.0d0/fv(sR)%u(1,eR)
        nu  =  0.5d0*(fv(sL)%mu(1,eL)*rhL+fv(sR)%mu(1,eR)*rhR)
        if(is_LES_now .or. is_DDES_now) then
            nut =  0.5d0*(fv(sL)%mu(2,eL)*rhL+fv(sR)%mu(2,eR)*rhR)
        else
            nut =  0.0d0
        end if
        g(1:9)  =  0.5d0*(fv(sL)%gra(4:12,eL)+fv(sR)%gra(4:12,eR))
        call fv_get_blending_vortex(gls,mesh(lev)%mortar_n_vg(5,im), v,g,nu,nut, LD(1))
        if(rhs_lev(lev) .eq. 2) then
            call rhs_conv_roe_ld (calD, addD, 1, LD(1))
        else
            call rhs_conv_AUSM_ld(calD, addD, 1, LD(1))
        end if

        fv(sL)%rhs(1:5,eL)  =  fv(sL)%rhs(1:5,eL)+rhsl(1:5,1)
        if(sec(sR)%is_int)  fv(sR)%rhs(1:5,eR)  =  fv(sR)%rhs(1:5,eR)-rhsl(1:5,1)
        if(savD) then
            fv(sL)%duc(1:5,eL)  =  fv(sL)%duc(1:5,eL)+rhsD(1:5,1)
            if(sec(sR)%is_int)  fv(sR)%duc(1:5,eR)  =  fv(sR)%duc(1:5,eR)-rhsD(1:5,1)
        end if
    end do

    return
    end subroutine fv_get_rhs_upw_LD
!-------------------------------------------------------------------------------
!   cal Rhs, convective part, fv, low-dissipation based on mesh.
!-------------------------------------------------------------------------------
    subroutine fv_get_rhs_upw_LDmesh(lev,calD,addD)
    use var_kind_def
    use var_fv
    use var_global, only: n_dim,R16,R23,uref
    use var_mesh
    use var_temp, only: uLf,uRf,vgf,fac_1d,rhsl,rhsD,u=>u_structured
    use var_turb, only: is_DDES_now,is_LES_now
    implicit none
    real   (dpR),parameter:: min_sigma =  2.5d-2
    logical(dpL),intent(in):: calD,addD
    integer(dpI),intent(in):: lev
    logical(dpL):: savD
    integer(dpI):: im,sL,eL,sR,eR,j,s_LL,e_LL,s_RR,e_RR
    real   (dpR):: dL(3),dR(3),limiter_L,limiter_R,gls,g(9),nu,nut,LD(100),rhL,rhR, &
                &  duL(5),duR(5),du(5),v(3),div2,Re,vor(3),vor2,e,d(3)

    dL          =  0.0d0
    dR          =  0.0d0
    limiter_L   =  1.0d0
    limiter_R   =  1.0d0
    LD          =  1.0d0
    savD=  calD .and. (.not. addD)

    call fv_get_rhs_boundary(lev, calD, addD)

    do im=1+mesh(lev)%n_mortar_b,mesh(lev)%n_mortar_ss
        s_LL=  mesh(lev)%mortar_structured_stencil(1,im)
        e_LL=  mesh(lev)%mortar_structured_stencil(2,im)
        sL  =  mesh(lev)%mortar_structured_stencil(3,im)
        eL  =  mesh(lev)%mortar_structured_stencil(4,im)
        sR  =  mesh(lev)%mortar_structured_stencil(5,im)
        eR  =  mesh(lev)%mortar_structured_stencil(6,im)
        s_RR=  mesh(lev)%mortar_structured_stencil(7,im)
        e_RR=  mesh(lev)%mortar_structured_stencil(8,im)
        fac_1d(1:4,1)   =  mesh(lev)%mortar_n_vg(1:4,im)
        vgf   (1:3,1)   =  mesh(lev)%mortar_n_vg(6:8,im)
        gls             =  1.0d0/mesh(lev)%mortar_n_vg(5,im)

        u(1:5,0)=  fv(s_LL)%u(1:5,e_LL)
        u(1:5,1)=  fv(sL  )%u(1:5,eL  )
        u(1:5,2)=  fv(sR  )%u(1:5,eR  )
        u(1:5,3)=  fv(s_RR)%u(1:5,e_RR)
        v(1:3)  =  0.5d0*(u(2:4,1)+u(2:4,2))
        if(is_limiter_on) then
            call recons_3rd_c(fac_1d, vgf, uLf, uRf)
        else
            uLf(1:5,1)  =(     -u(1:5,0)+5.0d0*u(1:5,1)+2.0d0*u(1:5,2))*R16
            uRf(1:5,1)  =(2.0d0*u(1:5,1)+5.0d0*u(1:5,2)      -u(1:5,3))*R16
        end if

!       Ducros approach.
        rhL =  1.0d0/fv(sL)%u(1,eL)
        rhR =  1.0d0/fv(sR)%u(1,eR)
        nu  =  0.5d0*(fv(sL)%mu(1,eL)*rhL+fv(sR)%mu(1,eR)*rhR)
        if(is_LES_now .or. is_DDES_now) then
            nut =  0.5d0*(fv(sL)%mu(2,eL)*rhL+fv(sR)%mu(2,eR)*rhR)
        else
            nut =  0.0d0
        end if
        g(1:9)  =  0.5d0*(fv(sL)%gra(4:12,eL)+fv(sR)%gra(4:12,eR))
        div2    = (g(1)+g(5)+g(9))**2
        vor(1)  =  g(8)-g(6)
        vor(2)  =  g(3)-g(7)
        vor(3)  =  g(4)-g(2)
        vor2    =  vor(1)*vor(1)+vor(2)*vor(2)+vor(3)*vor(3)
        Re      =  uref*gls/(nu+nut)
        e       = (1.0d-6+1.0d0/Re)*(uref/gls)**2
        LD(1)   =  max(min_sigma, div2/(div2+vor2+e))

!       mesh based approach.
        v(1:3)  =  sec(sL)%cen(1:3,eL)-sec(sR)%cen(1:3,eR)
        d(1:3)  = (sec(sL)%cen(1:3,eL)+sec(sR)%cen(1:3,eR))*0.5d0 &
                & -mesh(0)%mortar_cen(1:3,im)
        LD(1)   =  max(LD(1), sqrt(d(1)**2+d(2)**2+d(3)**2)/sqrt(v(1)**2+v(2)**2+v(3)**2))

        if(rhs_lev(lev) .eq. 2) then
            call rhs_conv_roe_ld (calD, addD, 1, LD(1))
        else
            call rhs_conv_AUSM_ld(calD, addD, 1, LD(1))
        end if

        fv(sL)%rhs(1:5,eL)  =  fv(sL)%rhs(1:5,eL)+rhsl(1:5,1)
        if(sec(sR)%is_int)  fv(sR)%rhs(1:5,eR)  =  fv(sR)%rhs(1:5,eR)-rhsl(1:5,1)
        if(savD) then
            fv(sL)%duc(1:5,eL)  =  fv(sL)%duc(1:5,eL)+rhsD(1:5,1)
            if(sec(sR)%is_int)  fv(sR)%duc(1:5,eR)  =  fv(sR)%duc(1:5,eR)-rhsD(1:5,1)
        end if
    end do

    do im=1+mesh(lev)%n_mortar_ss+1,mesh(lev)%n_mortar
        sL  =  mesh(lev)%mortar_LR(1,im)
        eL  =  mesh(lev)%mortar_LR(2,im)
        sR  =  mesh(lev)%mortar_LR(3,im)
        eR  =  mesh(lev)%mortar_LR(4,im)
        fac_1d(1:4,1)   =  mesh(lev)%mortar_n_vg(1:4,im)
        vgf   (1:3,1)   =  mesh(lev)%mortar_n_vg(6:8,im)
        gls             =  1.0d0/mesh(lev)%mortar_n_vg(5,im)

        uLf(1:5,1)  =  fv(sL)%u(1:5,eL)
        uRf(1:5,1)  =  fv(sR)%u(1:5,eR)

        dL(1:n_dim) =  mesh(lev)%mortar_cen(1:n_dim,im)-sec(sL)%cen(1:n_dim,eL)
        dR(1:n_dim) =  mesh(lev)%mortar_cen(1:n_dim,im)-sec(sR)%cen(1:n_dim,eR)
        do j=1,5
            duL(j)  =  fv(sL)%gra(3*j-2,eL)*dL(1)+fv(sL)%gra(3*j-1,eL)*dL(2) &
                    & +fv(sL)%gra(3*j  ,eL)*dL(3)
            duR(j)  =  fv(sR)%gra(3*j-2,eR)*dR(1)+fv(sR)%gra(3*j-1,eR)*dR(2) &
                    & +fv(sR)%gra(3*j  ,eR)*dR(3)
        end do
        if(is_limiter_on) then
            limiter_L   =  fv(sL)%gra(19,eL)
            limiter_R   =  fv(sR)%gra(19,eR)
        end if
        if(fv(sL)%order .eq. 2) then
            uLf(1:5,1)  =  uLf(1:5,1)+limiter_L*duL
            uRf(1:5,1)  =  uRf(1:5,1)+limiter_R*duR
        else
            du (1:5  )  =  uRf(1:5,1)-uLf(1:5,1)
            uLf(1:5,1)  =  uLf(1:5,1)+R23*duL(1:5)+R16*du(1:5)
            uRf(1:5,1)  =  uRf(1:5,1)+R23*duR(1:5)-R16*du(1:5)
        end if

        if((uLf(1,1) .le. 0.0d0) .or. (uLf(5,1) .le. 0.0d0) .or. &
         & (uRf(1,1) .le. 0.0d0) .or. (uRf(5,1) .le. 0.0d0)) then
            uLf(1:5,1)  =  fv(sL)%u(1:5,eL)
            uRf(1:5,1)  =  fv(sR)%u(1:5,eR)
        end if

!       Ducros approach.
        rhL =  1.0d0/fv(sL)%u(1,eL)
        rhR =  1.0d0/fv(sR)%u(1,eR)
        nu  =  0.5d0*(fv(sL)%mu(1,eL)*rhL+fv(sR)%mu(1,eR)*rhR)
        if(is_LES_now .or. is_DDES_now) then
            nut =  0.5d0*(fv(sL)%mu(2,eL)*rhL+fv(sR)%mu(2,eR)*rhR)
        else
            nut =  0.0d0
        end if
        g(1:9)  =  0.5d0*(fv(sL)%gra(4:12,eL)+fv(sR)%gra(4:12,eR))
        div2    = (g(1)+g(5)+g(9))**2
        vor(1)  =  g(8)-g(6)
        vor(2)  =  g(3)-g(7)
        vor(3)  =  g(4)-g(2)
        vor2    =  vor(1)*vor(1)+vor(2)*vor(2)+vor(3)*vor(3)
        Re      =  uref*gls/(nu+nut)
        e       = (1.0d-6+1.0d0/Re)*(uref/gls)**2
        LD(1)   =  max(min_sigma, div2/(div2+vor2+e))

!       mesh based approach.
        v(1:3)  =  sec(sL)%cen(1:3,eL)-sec(sR)%cen(1:3,eR)
        d(1:3)  = (sec(sL)%cen(1:3,eL)+sec(sR)%cen(1:3,eR))*0.5d0 &
                & -mesh(0)%mortar_cen(1:3,im)
        LD(1)   =  max(LD(1), sqrt(d(1)**2+d(2)**2+d(3)**2)/sqrt(v(1)**2+v(2)**2+v(3)**2))

        if(rhs_lev(lev) .eq. 2) then
            call rhs_conv_roe_ld (calD, addD, 1, LD(1))
        else
            call rhs_conv_AUSM_ld(calD, addD, 1, LD(1))
        end if

        fv(sL)%rhs(1:5,eL)  =  fv(sL)%rhs(1:5,eL)+rhsl(1:5,1)
        if(sec(sR)%is_int)  fv(sR)%rhs(1:5,eR)  =  fv(sR)%rhs(1:5,eR)-rhsl(1:5,1)
        if(savD) then
            fv(sL)%duc(1:5,eL)  =  fv(sL)%duc(1:5,eL)+rhsD(1:5,1)
            if(sec(sR)%is_int)  fv(sR)%duc(1:5,eR)  =  fv(sR)%duc(1:5,eR)-rhsD(1:5,1)
        end if
    end do

    return
    end subroutine fv_get_rhs_upw_LDmesh

!-------------------------------------------------------------------------------
!   cal Rhs, convective part, fv. Sunway version
!-------------------------------------------------------------------------------
    subroutine fv_get_rhs_upw_sw(lev,calD,addD)
    use var_kind_def
    use var_fv
    use var_fv_array
    use var_global, only: n_dim,R16,R23,sw_time
    use var_mesh
    use var_cgns
    use var_sec_array
    use var_parallel
    use var_temp, only: uLf,uRf,vgf,fac_1d,rhsl,rhsD,u=>u_structured
    use var_turb, only: low_dissipation
    use var_uns_cal, only: is_uns_cal_now
    use var_global_real
    implicit none
    logical(dpL),intent(in):: calD,addD
    integer(dpI),intent(in):: lev
    logical(dpL):: savD
    integer(dpI):: im,s_LL,e_LL,sL,eL,sR,eR,s_RR,e_RR,j
    real   (dpR):: dL(3),dR(3),limiter_L,limiter_R,duL(5),duR(5),du(5),tmp
    integer(dpI):: imi,idxL,idxR,bct
    real   (dpR):: calD_r, addD_r
    logical(dpL):: is_upwind
    integer(kind=8):: t1,t2,t3,t4
    character(20):: upw_slave = 'upw_slave'

    if(is_uns_cal_now) then
        if(low_dissipation .eq. 1) then
            stop 'fv_get_rhs_upw_LD has not implemented on Sunway platform'
            call fv_get_rhs_upw_LD    (lev, calD, addD)
        elseif(low_dissipation .eq. 2) then
            stop 'fv_get_rhs_upw_LDmesh has not implemented on Sunway platform'
            call fv_get_rhs_upw_LDmesh(lev, calD, addD)
        end if
        return
    end if

    dL          =  0.0d0
    dR          =  0.0d0
    limiter_L   =  1.0d0
    limiter_R   =  1.0d0
    savD=  calD .and. (.not. addD)

    ! call fv_get_rhs_boundary(lev, calD, addD)
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
        fac_1d(1:4,1)   =  mesh(lev)%mortar_n_vg(1:4,im)
        vgf   (1:3,1)   =  mesh(lev)%mortar_n_vg(6:8,im)
        bct =  sec(sR)%bct

        is_upwind   = (bct .eq. BCSymmetryPlane)
        if(is_upwind) then
            uLf(1:5,1)  =  fv_uL(1:5,idxR)
            uRf(1:5,1)  =  fv_uR(1:5,idxR)
            if(rhs_lev(lev) .eq. 3) then
                call rhs_conv_hllc(calD, addD, 1)
            elseif(rhs_lev(lev) .eq. 4) then
                call rhs_conv_ausm(calD, addD, 1)
            else
                call rhs_conv_roe (calD, addD, 1)
            end if
        else
            call u_to_F(1, fac_1d(1:3,1), fv_u(1,idxR), vgf, rhsl)
            rhsl(1:5,1) =  rhsl(1:5,1)*fac_1d(4,1)
            rhsD(1:5,1) =  0.0d0
        end if

        fv_rhs(1:5,idxL)  =  fv_rhs(1:5,idxL)+rhsl(1:5,1)
        if(savD)    fv_duc(1:5,idxL)  =  fv_duc(1:5,idxL)+rhsD(1:5,1)
    end do

    do im=1+mesh(lev)%n_mortar_b,mesh(lev)%n_mortar_ss
        stop 'mortar_structured_stencil has not supported on Sunway platform'
    end do

    

    if(calD) then
        calD_r = 1
    else
        calD_r = -1
    end if
    if(addD) then
        addD_r = 1
    else
        addD_r = -1
    end if

call starttime(upw_slave)
    call upw_conv_host(mesh_reordered(lev)%owner, mesh_reordered(lev)%neighbor, &
        & faceNum, cellNum, mesh_reordered(lev)%mortar_transform, &
        & mesh_reordered(lev)%mortar_n_vg, fv_u, mesh_reordered(lev)%mortar_cen, &
        sec_cen, fv_gra, fv_rhs, fv_duc, fv_order, sec_is_int, rhs_lev_r(lev), &
        is_limiter_on_r, calD_r, addD_r)
call endtime(upw_slave)

if(sw_time) then
    call system_clock(t3)

    do im=1+mesh(lev)%n_mortar_ss,mesh(lev)%n_mortar
        sL  =  mesh(lev)%mortar_LR(1,im)
        eL  =  mesh(lev)%mortar_LR(2,im)
        sR  =  mesh(lev)%mortar_LR(3,im)
        eR  =  mesh(lev)%mortar_LR(4,im)
        fac_1d(1:4,1)   =  mesh(lev)%mortar_n_vg(1:4,im)
        vgf   (1:3,1)   =  mesh(lev)%mortar_n_vg(6:8,im)

        uLf(1:5,1)  =  fv(sL)%u(1:5,eL)
        uRf(1:5,1)  =  fv(sR)%u(1:5,eR)
        if(fv(sL)%order .gt. 1) then
            dL(1:n_dim) =  mesh(lev)%mortar_cen(1:n_dim,im)-sec(sL)%cen(1:n_dim,eL)
            dR(1:n_dim) =  mesh(lev)%mortar_cen(1:n_dim,im)-sec(sR)%cen(1:n_dim,eR)
            do j=1,5
                duL(j)  =  fv(sL)%gra(3*j-2,eL)*dL(1)+fv(sL)%gra(3*j-1,eL)*dL(2) &
                        & +fv(sL)%gra(3*j  ,eL)*dL(3)
                duR(j)  =  fv(sR)%gra(3*j-2,eR)*dR(1)+fv(sR)%gra(3*j-1,eR)*dR(2) &
                        & +fv(sR)%gra(3*j  ,eR)*dR(3)
            end do
            if(is_limiter_on) then
                limiter_L   =  fv(sL)%gra(19,eL)
                limiter_R   =  fv(sR)%gra(19,eR)
            end if
            if(fv(sL)%order .eq. 2) then
                uLf(1:5,1)  =  uLf(1:5,1)+limiter_L*duL
                uRf(1:5,1)  =  uRf(1:5,1)+limiter_R*duR
            else
                du (1:5  )  =  uRf(1:5,1)-uLf(1:5,1)
                uLf(1:5,1)  =  uLf(1:5,1)+R23*duL(1:5)+R16*du(1:5)
                uRf(1:5,1)  =  uRf(1:5,1)+R23*duR(1:5)-R16*du(1:5)
            end if

            if((uLf(1,1) .le. 0.0d0) .or. (uLf(5,1) .le. 0.0d0) .or. &
             & (uRf(1,1) .le. 0.0d0) .or. (uRf(5,1) .le. 0.0d0)) then
                uLf(1:5,1)  =  fv(sL)%u(1:5,eL)
                uRf(1:5,1)  =  fv(sR)%u(1:5,eR)
            end if
        end if

        if(rhs_lev(lev) .eq. 2) then
            call rhs_conv_roe (calD, addD, 1)
        elseif(rhs_lev(lev) .eq. 3) then
            call rhs_conv_hllc(calD, addD, 1)
        elseif(rhs_lev(lev) .eq. 4) then
            call rhs_conv_ausm(calD, addD, 1)
        else
            call rhs_conv_lf  (calD, addD, 1)
        end if

        fv(sL)%rhs(1:5,eL)  =  fv(sL)%rhs(1:5,eL)+rhsl(1:5,1)
        if(sec(sR)%is_int)  fv(sR)%rhs(1:5,eR)  =  fv(sR)%rhs(1:5,eR)-rhsl(1:5,1)
        if(savD) then
            fv(sL)%duc(1:5,eL)  =  fv(sL)%duc(1:5,eL)+rhsD(1:5,1)
            if(sec(sR)%is_int)  fv(sR)%duc(1:5,eR)  =  fv(sR)%duc(1:5,eR)-rhsD(1:5,1)
        end if
    end do

    call system_clock(t4)
    write(*,*),'Processor ID: ', myid, ', Speed-up of convective term: ',real(t4-t3)/(t2-t1)
! write(*,*),'cpu_time: ', end-start
end if

!     imi = 0
!     do im=1+mesh(lev)%n_mortar_ss,mesh(lev)%n_mortar
!         imi = imi+1
!         ! sL  =  mesh(lev)%mortar_LR(1,im)
!         ! sR  =  mesh(lev)%mortar_LR(3,im)
!         if(mesh_reordered(lev)%mortar_transform(imi) .eq. 0) then
!             idxL = mesh_reordered(lev)%mortar_own_ID(imi)
!             idxR = mesh_reordered(lev)%mortar_nei_ID(imi)   
!         else 
!             idxL = mesh_reordered(lev)%mortar_nei_ID(imi)
!             idxR = mesh_reordered(lev)%mortar_own_ID(imi)   
!         end if
!         fac_1d(1:4,1) = mesh_reordered(lev)%mortar_n_vg(1:4,imi)
!         vgf   (1:3,1) = mesh_reordered(lev)%mortar_n_vg(6:8,imi)
!         uLf(1:5,1)  =  fv_u(1:5,idxL)
!         uRf(1:5,1)  =  fv_u(1:5,idxR)
!         if(fv_order(idxL) .gt. 1) then
!             ! if(myid==3) write(*,*),idxL,tot_ele,size(sec_cen)
!             ! if(myid==3) sec_cen(1,idxL) = 1.0d0
!             dL(1:n_dim) =  mesh_reordered(lev)%mortar_cen(1:n_dim,imi)-sec_cen(1:n_dim,idxL)
!             dR(1:n_dim) =  mesh_reordered(lev)%mortar_cen(1:n_dim,imi)-sec_cen(1:n_dim,idxR)
!             do j=1,5
!                 duL(j)  =  fv_gra(3*j-2,idxL)*dL(1)+fv_gra(3*j-1,idxL)*dL(2) &
!                         & +fv_gra(3*j  ,idxL)*dL(3)
!                 duR(j)  =  fv_gra(3*j-2,idxR)*dR(1)+fv_gra(3*j-1,idxR)*dR(2) &
!                         & +fv_gra(3*j  ,idxR)*dR(3)
!             end do
!             if(is_limiter_on) then
!                 limiter_L   =  fv_gra(19,idxL)
!                 limiter_R   =  fv_gra(19,idxR)
!             end if
!             if(fv_order(idxL) .eq. 2) then
!                 uLf(1:5,1)  =  uLf(1:5,1)+limiter_L*duL
!                 uRf(1:5,1)  =  uRf(1:5,1)+limiter_R*duR
!             else
!                 du (1:5  )  =  uRf(1:5,1)-uLf(1:5,1)
!                 uLf(1:5,1)  =  uLf(1:5,1)+R23*duL(1:5)+R16*du(1:5)
!                 uRf(1:5,1)  =  uRf(1:5,1)+R23*duR(1:5)-R16*du(1:5)
!             end if
!             if((uLf(1,1) .le. 0.0d0) .or. (uLf(5,1) .le. 0.0d0) .or. &
!              & (uRf(1,1) .le. 0.0d0) .or. (uRf(5,1) .le. 0.0d0)) then
!                 ! uLf(1:5,1)  =  fv_u(1:5,idxL)
!                 ! uRf(1:5,1)  =  fv_u(1:5,idxR)
!                 uLf(1:5,1)  =  fv_u(1:5,idxL)
!                 uRf(1:5,1)  =  fv_u(1:5,idxR)
!             end if
!         end if     
!         if(rhs_lev(lev) .eq. 2) then
!             call rhs_conv_roe (calD, addD, 1)
!         elseif(rhs_lev(lev) .eq. 3) then
!             call rhs_conv_hllc(calD, addD, 1)
!         elseif(rhs_lev(lev) .eq. 4) then
!             call rhs_conv_ausm(calD, addD, 1)
!         else
!             call rhs_conv_lf  (calD, addD, 1)
!         end if

!         fv_rhs(1:5,idxL)  =  fv_rhs(1:5,idxL)+rhsl(1:5,1)
!         if(sec_is_int(idxR)==1)  then
!             fv_rhs(1:5,idxR)  =  fv_rhs(1:5,idxR)-rhsl(1:5,1)
!         end if
!         ! write(*,*),savD
!         if(savD) then
!             fv_duc(1:5,idxL)  =  fv_duc(1:5,idxL)+rhsD(1:5,1)
!             if(sec_is_int(idxR)==1)  then
!                 fv_duc(1:5,idxR)  =  fv_duc(1:5,idxR)-rhsD(1:5,1)
!             end if
!         end if

!         ! fv_rhs(1:5,idxL) = fv_rhs(1:5,idxL) - rhsl(1:5,1)
!         ! fv_rhs(1:5,idxR) = fv_rhs(1:5,idxR) + rhsl(1:5,1)
!     end do
! ! end if

    ! call fv_array_to_struct(lev)

    ! imi = 0
    !  do im=1+mesh(lev)%n_mortar_ss,mesh(lev)%n_mortar
    !     imi = imi+1
    !     sL  =  mesh(lev)%mortar_LR(1,im)
    !     eL  =  mesh(lev)%mortar_LR(2,im)
    !     sR  =  mesh(lev)%mortar_LR(3,im)
    !     eR  =  mesh(lev)%mortar_LR(4,im)
    !     fac_1d(1:4,1)   =  mesh(lev)%mortar_n_vg(1:4,im)
    !     vgf   (1:3,1)   =  mesh(lev)%mortar_n_vg(6:8,im)
    !     uLf(1:5,1)  =  fv(sL)%u(1:5,eL)
    !     uRf(1:5,1)  =  fv(sR)%u(1:5,eR)
    !     ! if(fv(sL)%order .gt. 1) then
    !         dL(1:n_dim) =  mesh(lev)%mortar_cen(1:n_dim,im)-sec(sL)%cen(1:n_dim,eL)
    !         dR(1:n_dim) =  mesh(lev)%mortar_cen(1:n_dim,im)-sec(sR)%cen(1:n_dim,eR)
    !         do j=1,5
    !             duL(j)  =  fv(sL)%gra(3*j-2,eL)*dL(1)+fv(sL)%gra(3*j-1,eL)*dL(2) &
    !                     & +fv(sL)%gra(3*j  ,eL)*dL(3)
    !             duR(j)  =  fv(sR)%gra(3*j-2,eR)*dR(1)+fv(sR)%gra(3*j-1,eR)*dR(2) &
    !                     & +fv(sR)%gra(3*j  ,eR)*dR(3)
    !         end do
    !         if(is_limiter_on) then
    !             ! limiter_L   =  fv_gra(19,idxL)
    !             ! limiter_R   =  fv_gra(19,idxR)
    !             limiter_L   =  fv(sL)%gra(19,eL)
    !             limiter_R   =  fv(sR)%gra(19,eR)
    !         end if
    !         if(fv(sL)%order .eq. 2) then
    !             uLf(1:5,1)  =  uLf(1:5,1)+limiter_L*duL
    !             uRf(1:5,1)  =  uRf(1:5,1)+limiter_R*duR
    !         else
    !             du (1:5  )  =  uRf(1:5,1)-uLf(1:5,1)
    !             uLf(1:5,1)  =  uLf(1:5,1)+R23*duL(1:5)+R16*du(1:5)
    !             uRf(1:5,1)  =  uRf(1:5,1)+R23*duR(1:5)-R16*du(1:5)
    !         end if
    !         if((uLf(1,1) .le. 0.0d0) .or. (uLf(5,1) .le. 0.0d0) .or. &
    !          & (uRf(1,1) .le. 0.0d0) .or. (uRf(5,1) .le. 0.0d0)) then
    !             ! uLf(1:5,1)  =  fv_u(1:5,idxL)
    !             ! uRf(1:5,1)  =  fv_u(1:5,idxR)
    !             uLf(1:5,1)  =  fv(sL)%u(1:5,eL)
    !             uRf(1:5,1)  =  fv(sR)%u(1:5,eR)
    !         end if
    !     ! end if     
    !     if(rhs_lev(lev) .eq. 2) then
    !         call rhs_conv_roe (calD, addD, 1)
    !     elseif(rhs_lev(lev) .eq. 3) then
    !         call rhs_conv_hllc(calD, addD, 1)
    !     elseif(rhs_lev(lev) .eq. 4) then
    !         call rhs_conv_ausm(calD, addD, 1)
    !     else
    !         call rhs_conv_lf  (calD, addD, 1)
    !     end if

    !     fv(sL)%rhs(1:5,eL)  =  fv(sL)%rhs(1:5,eL)-rhsl(1:5,1)
    !     if(sec(sR)%is_int)  fv(sR)%rhs(1:5,eR)  =  fv(sR)%rhs(1:5,eR)+rhsl(1:5,1)
    !     if(savD) then
    !         fv(sL)%duc(1:5,eL)  =  fv(sL)%duc(1:5,eL)-rhsD(1:5,1)
    !         if(sec(sR)%is_int)  fv(sR)%duc(1:5,eR)  =  fv(sR)%duc(1:5,eR)+rhsD(1:5,1)
    !     end if
    !     ! fv(sL)%rhs(1:5,eL) = fv(sL)%rhs(1:5,eL) + rhsl(1:5,1)
    !     ! fv(sR)%rhs(1:5,eR) = fv(sR)%rhs(1:5,eR) + rhsl(1:5,1)
    ! end do
       
    ! imi = 0
    ! do im=1+mesh(lev)%n_mortar_ss,mesh(lev)%n_mortar
    !     imi =  imi+1
    !     sL  =  mesh(lev)%mortar_LR(1,im)
    !     eL  =  mesh(lev)%mortar_LR(2,im)
    !     sR  =  mesh(lev)%mortar_LR(3,im)
    !     eR  =  mesh(lev)%mortar_LR(4,im)
    !     idxL = mesh(lev)%mortar_own_ID(imi)
    !     idxR = mesh(lev)%mortar_nei_ID(imi)
    !     fac_1d(1:4,1)   =  mesh(lev)%mortar_n_vg(1:4,im)
    !     vgf   (1:3,1)   =  mesh(lev)%mortar_n_vg(6:8,im)

    !     uLf(1:5,1)  =  fv(sL)%u(1:5,eL)
    !     uRf(1:5,1)  =  fv(sR)%u(1:5,eR)
    !     if(fv(sL)%order .gt. 1) then
    !         dL(1:n_dim) =  mesh(lev)%mortar_cen(1:n_dim,im)-sec(sL)%cen(1:n_dim,eL)
    !         dR(1:n_dim) =  mesh(lev)%mortar_cen(1:n_dim,im)-sec(sR)%cen(1:n_dim,eR)
    !         do j=1,5
    !             ! duL(j)  =  fv_gra(3*j-2,idxL)*dL(1)+fv_gra(3*j-1,idxL)*dL(2) &
    !                     ! & +fv_gra(3*j  ,idxL)*dL(3)
    !             ! duR(j)  =  fv_gra(3*j-2,idxR)*dR(1)+fv_gra(3*j-1,idxR)*dR(2) &
    !                     ! & +fv_gra(3*j  ,idxR)*dR(3)
    !             duL(j)  =  fv(sL)%gra(3*j-2,eL)*dL(1)+fv(sL)%gra(3*j-1,eL)*dL(2) &
    !                     & +fv(sL)%gra(3*j  ,eL)*dL(3)
    !             duR(j)  =  fv(sR)%gra(3*j-2,eR)*dR(1)+fv(sR)%gra(3*j-1,eR)*dR(2) &
    !                     & +fv(sR)%gra(3*j  ,eR)*dR(3)
    !         end do
    !         if(is_limiter_on) then
    !             ! limiter_L   =  fv_gra(19,idxL)
    !             ! limiter_R   =  fv_gra(19,idxR)
    !             limiter_L   =  fv(sL)%gra(19,eL)
    !             limiter_R   =  fv(sR)%gra(19,eR)
    !         end if
    !         if(fv(sL)%order .eq. 2) then
    !             uLf(1:5,1)  =  uLf(1:5,1)+limiter_L*duL
    !             uRf(1:5,1)  =  uRf(1:5,1)+limiter_R*duR
    !         else
    !             du (1:5  )  =  uRf(1:5,1)-uLf(1:5,1)
    !             uLf(1:5,1)  =  uLf(1:5,1)+R23*duL(1:5)+R16*du(1:5)
    !             uRf(1:5,1)  =  uRf(1:5,1)+R23*duR(1:5)-R16*du(1:5)
    !         end if

    !         if((uLf(1,1) .le. 0.0d0) .or. (uLf(5,1) .le. 0.0d0) .or. &
    !          & (uRf(1,1) .le. 0.0d0) .or. (uRf(5,1) .le. 0.0d0)) then
    !             ! uLf(1:5,1)  =  fv_u(1:5,idxL)
    !             ! uRf(1:5,1)  =  fv_u(1:5,idxR)
    !             uLf(1:5,1)  =  fv(sL)%u(1:5,eL)
    !             uRf(1:5,1)  =  fv(sR)%u(1:5,eR)
    !         end if
    !     end if

    !     if(rhs_lev(lev) .eq. 2) then
    !         call rhs_conv_roe (calD, addD, 1)
    !     elseif(rhs_lev(lev) .eq. 3) then
    !         call rhs_conv_hllc(calD, addD, 1)
    !     elseif(rhs_lev(lev) .eq. 4) then
    !         call rhs_conv_ausm(calD, addD, 1)
    !     else
    !         call rhs_conv_lf  (calD, addD, 1)
    !     end if

    !     ! fv_rhs(1:5,idxL)  =  fv_rhs(1:5,idxL)+rhsl(1:5,1)
    !     ! if(sec(sR)%is_int)  then
    !         ! fv_rhs(1:5,idxR)  =  fv_rhs(1:5,idxR)-rhsl(1:5,1)
    !     ! end if
    !     ! if(savD) then
    !         ! fv_duc(1:5,idxL)  =  fv_duc(1:5,idxL)+rhsD(1:5,1)
    !         ! if(sec(sR)%is_int)  then
    !             ! fv_duc(1:5,idxR)  =  fv_duc(1:5,idxR)-rhsD(1:5,1)
    !         ! end if
    !     ! end if
    !     fv(sL)%rhs(1:5,eL)  =  fv(sL)%rhs(1:5,eL)+rhsl(1:5,1)
    !     if(sec(sR)%is_int)  fv(sR)%rhs(1:5,eR)  =  fv(sR)%rhs(1:5,eR)-rhsl(1:5,1)
    !     if(savD) then
    !         fv(sL)%duc(1:5,eL)  =  fv(sL)%duc(1:5,eL)+rhsD(1:5,1)
    !         if(sec(sR)%is_int)  fv(sR)%duc(1:5,eR)  =  fv(sR)%duc(1:5,eR)-rhsD(1:5,1)
    !     end if
    ! end do

    ! call fv_array_to_struct(lev)
    return
    end subroutine fv_get_rhs_upw_sw

