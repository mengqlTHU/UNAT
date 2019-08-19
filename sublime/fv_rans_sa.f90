!-------------------------------------------------------------------------------
!   get mut from nut.
!-------------------------------------------------------------------------------
    subroutine fv_SA_nut_to_mut(lev,mode)
    use var_kind_def
    use var_fv
    use var_mesh, only: mesh,sec
    implicit none
    integer(dpI),intent(in):: lev,mode
    integer(dpI):: isec,i,ele1,ele0
    real   (dpR):: rh,nut,mul,nu,chi3,fv1

    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(sec(isec)%is_ghost)  cycle
        if(mode .le. 0) then
            ele1=  1
            ele0=  sec(isec)%n_ele
        elseif(mode .eq. 1) then
            ele1=  1
            ele0=  sec(isec)%n_ele_b
        else
            ele1=  sec(isec)%n_ele_b+1
            ele0=  sec(isec)%n_ele
        end if

        do i=ele1,ele0
            rh  =  fv(isec)%u   (1,i)
            nut =  fv(isec)%turb(1,i)
            mul =  fv(isec)%mu  (1,i)
            nu  =  mul/rh
            chi3= (max(nut/nu, -5.0d0))**3
            fv1 =  chi3/(chi3+3.57911d2)
            fv(isec)%mu(2,i)=  rh*nut*fv1
        end do
    end do

    return
    end subroutine fv_SA_nut_to_mut
!-------------------------------------------------------------------------------
!   cal Rhs, turbulence.
!-------------------------------------------------------------------------------
    subroutine fv_SA_get_rhs(lev,is_LHS)
    use var_kind_def
    use var_fv
    use var_global, only: n_dim
    use var_mesh
    use var_parallel
    use var_prec
    use var_turb, transition=>transition_model
    use var_uns_cal, only: is_BDF_now,dt_uns,uns_iter,is_uns_initialized
    implicit none
    logical(dpL),intent(in):: is_LHS
    integer(dpI),intent(in):: lev
    integer(dpI):: isec,iele,im,sL,eL,sR,eR,iA(1000)
    real   (dpR):: n(5),uLf(5),uRf(5),tL,tR,nul,nut,gra(3),rL,rR,rV,d(3),g,dt, &
                &  nne(3),unL,unR,vol,dnw,rh,mul,mut,ux,uy,uz,vx,vy,vz,wx,wy, &
                &  wz,gran,oo,d2_1,chi,fv1,fv2,s,st,r,fw,pp,dd,chi_nut,fv1_nut, &
                &  fv2_nut,r_nut,g_nut,fw_nut,fLtL(1),fLtR(1),fRtL(1),fRtR(1), &
                &  vg(3),t1,t2,rtmp

    if(is_LHS) then
        if(.not. allocated(fv_rans_prec%iA))    call fv_rans_set_slv(lev, fv_rans_prec)
        fv_rans_prec%A  =  0.0d0

        iA  =  1
        iele=  0
        do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
            if(sec(isec)%is_ghost .or. sec(isec)%is_bnd)    cycle
            iA(isec)=  iele+1
            iele    =  iele+sec(isec)%n_ele
        end do
    end if

!   ----------------------------------------------------------------------------
!   source term.
    t1  =  1.0d0
    t2  =  1.0d0
    if(is_DDES_now) then
        call fv_DDES_get_src(lev, .true.)
    else
        do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(sec(isec)%is_ghost .or. sec(isec)%is_bnd)    cycle

        do iele=1,sec(isec)%n_ele
            vol =  sec(isec)%vol( iele)
            dnw =  sec(isec)%dnw( iele)
            rh  =  fv (isec)%u (1,iele)
            mul =  fv (isec)%mu(1,iele)
            mut =  fv (isec)%mu(2,iele)
            nul =  mul/rh
            nut =  fv (isec)%turb(1,iele)

            ux  =  fv(isec)%gra(4 ,iele)
            uy  =  fv(isec)%gra(5 ,iele)
            uz  =  fv(isec)%gra(6 ,iele)
            vx  =  fv(isec)%gra(7 ,iele)
            vy  =  fv(isec)%gra(8 ,iele)
            vz  =  fv(isec)%gra(9 ,iele)
            wx  =  fv(isec)%gra(10,iele)
            wy  =  fv(isec)%gra(11,iele)
            wz  =  fv(isec)%gra(12,iele)
            gran=  fv(isec)%turb_gra(1,iele)**2+fv(isec)%turb_gra(2,iele)**2 &
                & +fv(isec)%turb_gra(3,iele)**2
            oo  =  sqrt((uy-vx)**2+(uz-wx)**2+(vz-wy)**2)
            d2_1=  1.0d0/(dnw*dnw)

            chi =  nut/nul
            rtmp=  chi**3
            fv1 =  rtmp/(rtmp+3.57911d2)

            s   =  oo
            fv2 =  1.0d0-chi/(1.0d0+chi*fv1)
            rtmp=  nut*fv2*k2_1*d2_1
            if(rtmp .ge. -0.7d0*s) then
                st  =  s+rtmp
            else
                st  =  s+s*(0.49d0*s+0.9d0*rtmp)/(-0.5d0*s-rtmp)
            end if
            r   =  min(nut*k2_1*d2_1/st, 1.0d1)
            g   =  r+0.3d0*(r**6-r)
            fw  =  g*(6.5d1/(g**6+6.4d1))**(1.0d0/6.0d0)
            if(transition .eq. transition_PTM) then
                t1  =  fv(isec)%intermittency(iele)
                t2  =  max(0.02d0, t1)
            elseif(transition .eq. transition_BC) then
                t1  =  fv(isec)%intermittency(iele)
            end if
            pp  =  t1*(cb1*st*nut)+0.933d0*gran
            dd  =  t2*cw1*fw*d2_1*nut**2
            fv(isec)%rhs(1,iele)= -vol*(pp-dd)

            if(is_LHS) then
                chi_nut =  1.0d0/nul
                fv1_nut =  1073.733d0*chi_nut*chi*chi/(chi**3+3.57911d2)**2
                fv2_nut = -chi_nut/(1.0d0+chi*fv1)+chi*(fv1*chi_nut+chi*fv1_nut) &
                        & /(1.0d0+chi*fv1)**2
                if(r .lt. 1.0d1) then
                    r_nut   =  r/nut-r*r*(fv2/nut+fv2_nut)
                else
                    r_nut   =  0.0d0
                end if
                g_nut   = (1.0d0-cw2+6.0d0*cw2*r**5)*r_nut
                fw_nut  = (6.5d1/(6.4d1+g**6))**(1.0d0/6.0d0)*g_nut* &
                        & (1.0d0-g**6/(g**6+6.4d1))
                g       = -t1*(cb1*nut*nut*k2_1*d2_1*min(fv2_nut, 0.0d0)) &
                        & +t2*(2.0d0*cw1*fw*nut*d2_1+cw1*nut*nut*d2_1*max(fw_nut, 0.0d0))

                eL  =  iele+iA(isec)-1
                call prec_add_eleR(fv_rans_prec, eL, eL, (/vol*g/))
            end if
        end do
        end do
    end if
!   source term.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   convective and viscous terms.
    d   =  0.0d0
    do im=1,mesh(lev)%n_mortar
        sL      =  mesh(lev)%mortar_LR(1,im)
        eL      =  mesh(lev)%mortar_LR(2,im)
        sR      =  mesh(lev)%mortar_LR(3,im)
        eR      =  mesh(lev)%mortar_LR(4,im)
        n (1:5) =  mesh(lev)%mortar_n_vg(1:5,im)
        vg(1:3) =  mesh(lev)%mortar_n_vg(6:8,im)
        if(im .le. mesh(lev)%n_mortar_b) then
            uLf(1:5)=  fv(sR)%u(1:5,eR)
            uRf(1:5)=  fv(sR)%u(1:5,eR)
            tL      =  fv(sR)%turb(1,eR)
            tR      =  tL
            call DCOPY(3, fv(sR)%turb_gra(1,eR), 1, gra, 1)
            nul     =  fv(sR)%mu(1,eR)/uRf(1)
            nut     =  tL
        else
            uLf(1:5)=  fv(sL)%u(1:5,eL)
            uRf(1:5)=  fv(sR)%u(1:5,eR)
            tL      =  fv(sL)%turb(1,eL)
            tR      =  fv(sR)%turb(1,eR)
            nul     =  0.5d0*(fv(sL)%mu(1,eL)/uLf(1)+fv(sR)%mu(1,eR)/uRf(1))
            nut     =  0.5d0*(tL+tR)

            if(.true.) then
                gra(1:3)    = (tR-tL)*n(1:3)*n(5)
            else
                gra(1:3)    =  0.5d0*(fv(sL)%turb_gra(1:3,eL)+fv(sR)%turb_gra(1:3,eR))
                d(1:n_dim)  =  sec(sR)%cen(1:n_dim,eR)-sec(sL)%cen(1:n_dim,eL)
                rtmp        =  1.0d0/sqrt(d(1)**2+d(2)**2+d(3)**2)
                d           =  d*rtmp
                dt          =  rtmp*(tR-tL)
                rtmp        =  n(1)*d(1)+n(2)*d(2)+n(3)*d(3)
                nne(1:3)    =  n(1:3)/rtmp
                gra(1:3)    =  gra(1:3)+(dt-gra(1)*d(1)-gra(2)*d(2)-gra(3)*d(3))*nne(1:3)
            end if
        end if

        unL = (uLf(2)-vg(1))*n(1)+(uLf(3)-vg(2))*n(2)+(uLf(4)-vg(3))*n(3)
        unR = (uRf(2)-vg(1))*n(1)+(uRf(3)-vg(2))*n(2)+(uRf(4)-vg(3))*n(3)
        rV  = (n(1)*gra(1)+n(2)*gra(2)+n(3)*gra(3))*(nul+nut)*Pr_SA*n(4)

        fLtL=  0.0d0
        fLtR=  0.0d0
        fRtL=  0.0d0
        fRtR=  0.0d0
        if(unL .ge. 0.0d0) then
            rL  =  unL*tL*n(4)
            fLtL=  unL*n(4)
        else
            rL  =  unL*tR*n(4)
            fLtR=  unL*n(4)
        end if
        if(unR .ge. 0.0d0) then
            rR  =  unR*tL*n(4)
            fRtL= -unR*n(4)
        else
            rR  =  unR*tR*n(4)
            fRtR= -unR*n(4)
        end if
        fv(sL)%rhs(1,eL)=  fv(sL)%rhs(1,eL)+(rL-rV)
        if(sec(sR)%is_int)  fv(sR)%rhs(1,eR)=  fv(sR)%rhs(1,eR)-(rR-rV)
        if(is_LHS) then
            eL  =  eL+iA(sL)-1
            eR  =  eR+iA(sR)-1
            fLtL=  fLtL+n(4)*n(5)*Pr_SA*(nul+nut)
            fLtR=  fLtR-n(4)*n(5)*Pr_SA*(nul+nut)
            fRtR=  fRtR+n(4)*n(5)*Pr_SA*(nul+nut)
            fRtL=  fRtL-n(4)*n(5)*Pr_SA*(nul+nut)
            call prec_add_eleR(fv_rans_prec, eL, eL, fLtL)
            if(sec(sR)%is_int) then
                call prec_add_eleR(fv_rans_prec, eL, eR, fLtR)
                call prec_add_eleR(fv_rans_prec, eR, eR, fRtR)
                call prec_add_eleR(fv_rans_prec, eR, eL, fRtL)
            end if
        end if
    end do
!   convective and viscous terms.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   source term due to BDF.
    if(is_BDF_now) then
        dt  =  1.0d0/dt_uns

        do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
            if(sec(isec)%is_ghost .or. sec(isec)%is_bnd)    cycle

            do iele=1,sec(isec)%n_ele
                vol =  sec(isec)%vol( iele)
                if((uns_iter .le. 1) .and. (.not. is_uns_initialized)) then
                    s   =  fv(isec)%turb(1,iele)-fv(isec)%uns_turb(1,iele)
                else
                    s   =  1.5d0*fv(isec)%turb(1,iele)-2.0d0*fv(isec)%uns_turb(1,iele) &
                        & +0.5d0*fv(isec)%uns_turb(2,iele)
                end if
                fv(isec)%rhs(1,iele)=  fv(isec)%rhs(1,iele)+vol*s*dt

                if(is_LHS) then
                    eL  =  iele+iA(isec)-1
                    call prec_add_eleR(fv_rans_prec, eL, eL, (/1.5d0*vol*dt/))
                end if
            end do
        end do
    end if
!   source term due to BDF.
!   ----------------------------------------------------------------------------

    if(is_LHS)  call prec_get_decomposition(fv_rans_prec)

    res_RANS=  0.0d0
    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(sec(isec)%is_ghost .or. sec(isec)%is_bnd)    cycle

        do iele=1,sec(isec)%n_ele
            res_RANS=  res_RANS+fv(isec)%rhs(1,iele)**2
        end do
    end do

    return
    end subroutine fv_SA_get_rhs
!-------------------------------------------------------------------------------
!   DDES source term.
!-------------------------------------------------------------------------------
    subroutine fv_DDES_get_src(lev,is_LHS)
    use var_kind_def
    use var_fv
    use var_mesh
    use var_prec
    use var_turb
    implicit none
    logical(dpL),intent(in):: is_LHS
    integer(dpI),intent(in):: lev
    integer(dpI):: isec,iele,eL,iA(1000)

    if(is_LHS) then
        if(.not. allocated(fv_rans_prec%iA))    call fv_rans_set_slv(lev, fv_rans_prec)
        fv_rans_prec%A  =  0.0d0

        iA  =  1
        iele=  0
        do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
            if(sec(isec)%is_ghost .or. sec(isec)%is_bnd)    cycle
            iA(isec)=  iele+1
            iele    =  iele+sec(isec)%n_ele
        end do
    end if

    if(DDES_model .eq. DDES_Spalart) then
        call fv_DDES_Spalart(lev, is_LHS)
    elseif(DDES_model .eq. DDES_WALE) then
        call fv_DDES_WALE(lev, is_LHS)
    elseif(DDES_model .eq. DDES_sigma) then
        call fv_DDES_sigma(lev, is_LHS)
    elseif(DDES_model .eq. DDES_IDDES) then
        call fv_DDES_IDDES(lev, is_LHS)
    else
        stop 'Error: DDES model not supported.'
    end if

    return
    contains
!   ----------------------------------------------------------------------------
!   DDES source term, Spalart model.
!   ----------------------------------------------------------------------------
    subroutine fv_DDES_Spalart(lev,is_LHS)
    use var_global, only: L_ref,uref
    implicit none
    logical(dpL),intent(in):: is_LHS
    integer(dpI),intent(in):: lev
    integer(dpI):: isec,iele,eL
    real   (dpR):: nul,nut,g,ug(3,3),vol,dnw,rh,mul,mut,gran,oo,gls,rd,fd,d2_1, &
                &  chi,fv1,fv2,st,r,fw,pp,dd,uij,psi,L_RANS,L_LES,chi_nut,fv1_nut, &
                &  fv2_nut,r_nut,g_nut,fw_nut,rtmp

    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(sec(isec)%is_ghost .or. sec(isec)%is_bnd)    cycle

        do iele=1,sec(isec)%n_ele
            vol =  sec(isec)%vol( iele)
            dnw =  sec(isec)%dnw( iele)
            rh  =  fv (isec)%u(1 ,iele)
            mul =  fv (isec)%mu(1,iele)
            mut =  fv (isec)%mu(2,iele)
            nul =  mul/rh
            nut =  fv (isec)%turb(1,iele)

            ug(1,1:3)   =  fv(isec)%gra(4 :6 ,iele)
            ug(2,1:3)   =  fv(isec)%gra(7 :9 ,iele)
            ug(3,1:3)   =  fv(isec)%gra(10:12,iele)
            gran=  fv(isec)%turb_gra(1,iele)**2+fv(isec)%turb_gra(2,iele)**2 &
                & +fv(isec)%turb_gra(3,iele)**2
            oo  =  sqrt((ug(1,2)-ug(2,1))**2+(ug(1,3)-ug(3,1))**2+(ug(2,3)-ug(3,2))**2)
            gls =  sec(isec)%gls(iele)
            uij =  sqrt(ug(1,1)*ug(1,1)+ug(1,2)*ug(1,2)+ug(1,3)*ug(1,3) &
                &      +ug(2,1)*ug(2,1)+ug(2,2)*ug(2,2)+ug(2,3)*ug(2,3) &
                &      +ug(3,1)*ug(3,1)+ug(3,2)*ug(3,2)+ug(3,3)*ug(3,3))
            uij =  max(uij, 1.0d-6*uref/L_ref)
            chi =  nut/nul
            rtmp=  chi**3
            fv1 =  rtmp/(rtmp+3.57911d2)
            fv2 =  1.0d0-chi/(1.0d0+chi*fv1)
!           low Reynolds number correction.
            Psi =  sqrt(min(1.0d2, (1.0d0-cb1*fv2/(0.424d0*cw1*kappa*kappa))/fv1))

            L_RANS  =  dnw
            L_LES   =  Psi*C_DES*gls
            if(.true.) then
                rd  = (nul+nut)*k2_1/(uij*dnw*dnw)
                fd  =  1.0d0-tanh((8.0d0*rd)**3)
            else
                rd  =  nut*k2_1/(uij*dnw*dnw)
                fd  =  1.0d0-tanh((1.4d1*rd)**3)
            end if
            dnw =  L_RANS-fd*max(0.0d0, L_RANS-L_LES)
            d2_1=  1.0d0/(dnw*dnw)
            st  =  max(oo+nut*fv2*k2_1*d2_1, 0.3d0*oo)
            if(is_get_DDES_quality) fv(isec)%DDES_quality(3,iele)   =  1.0d0-fd

            r   =  min(nut*k2_1*d2_1/st, 1.0d1)
            g   =  r+0.3d0*(r**6-r)
            fw  =  g*(6.5d1/(g**6+6.4d1))**(1.0d0/6.0d0)
            pp  =  cb1*st*nut+0.933d0*gran
            dd  =  cw1*fw*d2_1*nut**2
            fv(isec)%rhs(1,iele)= -vol*(pp-dd)

            if(is_LHS) then
                chi_nut =  1.0d0/nul
                fv1_nut =  1073.733d0*chi_nut*chi*chi/(chi**3+3.57911d2)**2
                fv2_nut = -chi_nut/(1.0d0+chi*fv1)+chi*(fv1*chi_nut+chi*fv1_nut) &
                        & /(1.0d0+chi*fv1)**2
                if(r .lt. 1.0d1) then
                    r_nut   =  r/nut-r*r*(fv2/nut+fv2_nut)
                else
                    r_nut   =  0.0d0
                end if
                g_nut   = (1.0d0-cw2+6.0d0*cw2*r**5)*r_nut
                fw_nut  = (6.5d1/(6.4d1+g**6))**(1.0d0/6.0d0)*g_nut* &
                        & (1.0d0-g**6/(g**6+6.4d1))
                g       = -cb1*nut*nut*k2_1*d2_1*min(fv2_nut, 0.0d0) &
                        & +2.0d0*cw1*fw*nut*d2_1+cw1*nut*nut*d2_1*max(fw_nut, 0.0d0)

                eL  =  iele+iA(isec)-1
                call prec_add_eleR(fv_rans_prec, eL, eL, (/vol*g/))
            end if
        end do
    end do

    return
    end subroutine fv_DDES_Spalart
!   ----------------------------------------------------------------------------
!   DDES source term, WALE model.
!   ----------------------------------------------------------------------------
    subroutine fv_DDES_WALE(lev,is_LHS)
    use var_global, only: L_ref,uref,R13
    implicit none
    logical(dpL),intent(in):: is_LHS
    integer(dpI),intent(in):: lev
    integer(dpI):: isec,iele,eL,i,j,k
    real   (dpR):: nul,nut,g,ug(3,3),vol,dnw,rh,mul,mut,gran,oo,gls,rd,fd,d2_1, &
                &  chi,fv1,fv2,st,r,fw,pp,dd,uij,psi,L_RANS,L_LES,S_RANS,S_LES, &
                &  Sd(3,3),S(6),SS,SdSd,chi_nut,fv1_nut,fv2_nut,r_nut,g_nut,fw_nut, &
                &  rtmp

    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(sec(isec)%is_ghost .or. sec(isec)%is_bnd)    cycle

        do iele=1,sec(isec)%n_ele
            vol =  sec(isec)%vol( iele)
            dnw =  sec(isec)%dnw( iele)
            rh  =  fv (isec)%u(1 ,iele)
            mul =  fv (isec)%mu(1,iele)
            mut =  fv (isec)%mu(2,iele)
            nul =  mul/rh
            nut =  fv (isec)%turb(1,iele)

            ug(1,1:3)   =  fv(isec)%gra(4 :6 ,iele)
            ug(2,1:3)   =  fv(isec)%gra(7 :9 ,iele)
            ug(3,1:3)   =  fv(isec)%gra(10:12,iele)
            gran=  fv(isec)%turb_gra(1,iele)**2+fv(isec)%turb_gra(2,iele)**2 &
                & +fv(isec)%turb_gra(3,iele)**2
            oo  =  sqrt((ug(1,2)-ug(2,1))**2+(ug(1,3)-ug(3,1))**2+(ug(2,3)-ug(3,2))**2)
            gls =  sec(isec)%gls(iele)
            uij =  sqrt(ug(1,1)*ug(1,1)+ug(1,2)*ug(1,2)+ug(1,3)*ug(1,3) &
                &      +ug(2,1)*ug(2,1)+ug(2,2)*ug(2,2)+ug(2,3)*ug(2,3) &
                &      +ug(3,1)*ug(3,1)+ug(3,2)*ug(3,2)+ug(3,3)*ug(3,3))
            uij =  max(uij, 1.0d-6*uref/L_ref)
            chi =  nut/nul
            rtmp=  chi**3
            fv1 =  rtmp/(rtmp+3.57911d2)
            fv2 =  1.0d0-chi/(1.0d0+chi*fv1)
!           low Reynolds number correction.
            Psi =  sqrt(min(1.0d2, (1.0d0-cb1*fv2/(0.424d0*cw1*kappa*kappa))/fv1))

            L_RANS  =  dnw
            L_LES   =  Psi*C_DES*gls
            rd  = (nul+nut)*k2_1/(uij*dnw*dnw)
            fd  =  1.0d0-tanh((1.0d1*rd)**3)
            dnw =  L_RANS-fd*max(0.0d0, L_RANS-L_LES)
            d2_1=  1.0d0/(dnw*dnw)

            Sd  =  0.0d0
            rtmp=  0.0d0
            do j=1,3
            do i=1,3
                rtmp=  rtmp+ug(i,j)*ug(j,i)
                do k=1,3
                    Sd(i,j) =  Sd(i,j)+ug(i,k)*ug(k,j)+ug(j,k)*ug(k,i)
                end do
            end do
            end do
            Sd  =  Sd*0.5d0
            Sd(1,1) =  Sd(1,1)-R13*rtmp
            Sd(2,2) =  Sd(2,2)-R13*rtmp
            Sd(3,3) =  Sd(3,3)-R13*rtmp
            SdSd    =  Sd(1,1)*Sd(1,1)+Sd(1,2)*Sd(1,2)+Sd(1,3)*Sd(1,3) &
                    & +Sd(2,1)*Sd(2,1)+Sd(2,2)*Sd(2,2)+Sd(2,3)*Sd(2,3) &
                    & +Sd(3,1)*Sd(3,1)+Sd(3,2)*Sd(3,2)+Sd(3,3)*Sd(3,3)
            S(1)    =  ug(1,1)
            S(2)    =  0.5d0*(ug(2,1)+ug(1,2))
            S(3)    =  0.5d0*(ug(3,1)+ug(1,3))
            S(4)    =  ug(2,2)
            S(5)    =  0.5d0*(ug(3,2)+ug(2,3))
            S(6)    =  ug(3,3)
            SS      =  S(1)*S(1)+S(4)*S(4)+S(6)*S(6) &
                    & +2.0d0*(S(2)*S(2)+S(3)*S(3)+S(5)*S(5))
            S_LES   =  SdSd**1.5d0/(SS**2.5d0+SdSd**1.25d0)

            S_RANS  =  oo
            if(L_RANS .gt. L_LES) then
                rtmp=  1.0d0
            else
                rtmp=  0.0d0
            end if
            oo      =  S_RANS-fd*rtmp*(S_RANS-0.1192024d0*S_LES)
            st      =  oo+nut*fv2*k2_1*d2_1
            if(is_get_DDES_quality) fv(isec)%DDES_quality(3,iele)   =  1.0d0-fd

            r   =  min(nut*k2_1*d2_1/st, 1.0d1)
            g   =  r+0.3d0*(r**6-r)
            fw  =  g*(6.5d1/(g**6+6.4d1))**(1.0d0/6.0d0)
            pp  =  cb1*st*nut+0.933d0*gran
            dd  =  cw1*fw*d2_1*nut**2
            fv(isec)%rhs(1,iele)= -vol*(pp-dd)

            if(is_LHS) then
                chi_nut =  1.0d0/nul
                fv1_nut =  1073.733d0*chi_nut*chi*chi/(chi**3+3.57911d2)**2
                fv2_nut = -chi_nut/(1.0d0+chi*fv1)+chi*(fv1*chi_nut+chi*fv1_nut) &
                        & /(1.0d0+chi*fv1)**2
                if(r .lt. 1.0d1) then
                    r_nut   =  r/nut-r*r*(fv2/nut+fv2_nut)
                else
                    r_nut   =  0.0d0
                end if
                g_nut   = (1.0d0-cw2+6.0d0*cw2*r**5)*r_nut
                fw_nut  = (6.5d1/(6.4d1+g**6))**(1.0d0/6.0d0)*g_nut* &
                        & (1.0d0-g**6/(g**6+6.4d1))
                g       = -cb1*nut*nut*k2_1*d2_1*min(fv2_nut, 0.0d0) &
                        & +2.0d0*cw1*fw*nut*d2_1+cw1*nut*nut*d2_1*max(fw_nut, 0.0d0)

                eL  =  iele+iA(isec)-1
                call prec_add_eleR(fv_rans_prec, eL, eL, (/vol*g/))
            end if
        end do
    end do

    return
    end subroutine fv_DDES_WALE
!   ----------------------------------------------------------------------------
!   DDES source term, sigma model.
!   ----------------------------------------------------------------------------
    subroutine fv_DDES_sigma(lev,is_LHS)
    use var_global, only: L_ref,uref,R13,pi
    implicit none
    logical(dpL),intent(in):: is_LHS
    integer(dpI),intent(in):: lev
    integer(dpI):: isec,iele,eL
    real   (dpR):: nul,nut,g,ug(9),vol,dnw,rh,mul,mut,gran,oo,gls,rd,fd,d2_1, &
                &  chi,fv1,fv2,st,r,fw,pp,dd,uij,psi,L_RANS,L_LES,S_RANS,S_LES, &
                &  S(6),chi_nut,fv1_nut,tt1,tt2,tt3,alp1,alp2,alp3,alp12,sig1, &
                &  sig2,sig3,dsig,fv2_nut,r_nut,g_nut,fw_nut,rtmp

    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(sec(isec)%is_ghost .or. sec(isec)%is_bnd)    cycle

        do iele=1,sec(isec)%n_ele
            vol =  sec(isec)%vol( iele)
            dnw =  sec(isec)%dnw( iele)
            gls =  sec(isec)%gls( iele)
            rh  =  fv (isec)%u (1,iele)
            mul =  fv (isec)%mu(1,iele)
            mut =  fv (isec)%mu(2,iele)
            nul =  mul/rh
            nut =  fv (isec)%turb(1,iele)

            ug(1:9) =  fv(isec)%gra(4:12,iele)
            gran=  fv(isec)%turb_gra(1,iele)**2+fv(isec)%turb_gra(2,iele)**2 &
                & +fv(isec)%turb_gra(3,iele)**2
            oo  =  sqrt((ug(2)-ug(4))**2+(ug(3)-ug(7))**2+(ug(6)-ug(8))**2)
            uij =  sqrt(ug(1)*ug(1)+ug(2)*ug(2)+ug(3)*ug(3) &
                &      +ug(4)*ug(4)+ug(5)*ug(5)+ug(6)*ug(6) &
                &      +ug(7)*ug(7)+ug(8)*ug(8)+ug(9)*ug(9))
            uij =  max(uij, 1.0d-6*uref/L_ref)
            chi =  nut/nul
            rtmp=  chi**3
            fv1 =  rtmp/(rtmp+3.57911d2)
            fv2 =  1.0d0-chi/(1.0d0+chi*fv1)
!           low Reynolds number correction.
            Psi =  sqrt(min(1.0d2, (1.0d0-cb1*fv2/(0.424d0*cw1*kappa*kappa))/fv1))

            L_RANS  =  dnw
            L_LES   =  Psi*C_DES*gls

            rd  = (nul+nut)*k2_1/(uij*dnw*dnw)
            fd  =  1.0d0-tanh((1.0d1*rd)**3)
            dnw =  L_RANS-fd*max(0.0d0, L_RANS-L_LES)
            d2_1=  1.0d0/(dnw*dnw)

            S(1)=  ug(1)*ug(1)+ug(4)*ug(4)+ug(7)*ug(7)
            S(2)=  ug(1)*ug(2)+ug(4)*ug(5)+ug(7)*ug(8)
            S(3)=  ug(1)*ug(3)+ug(4)*ug(6)+ug(7)*ug(9)
            S(4)=  ug(2)*ug(2)+ug(5)*ug(5)+ug(8)*ug(8)
            S(5)=  ug(2)*ug(3)+ug(5)*ug(6)+ug(8)*ug(9)
            S(6)=  ug(3)*ug(3)+ug(6)*ug(6)+ug(9)*ug(9)
            tt1 =  S(1)+S(4)+S(6)
            tt2 =  S(1)*S(4)+S(1)*S(6)+S(4)*S(6)-(S(2)*S(2)+S(3)*S(3)+S(5)*S(5))
            tt3 =  S(1)*S(4)*S(6)+2.0d0*S(2)*S(3)*S(5) &
                & -S(4)*S(3)*S(3)-S(1)*S(5)*S(5)-S(6)*S(2)*S(2)
            alp1=  tt1*tt1/9.0d0-tt2*R13
            alp2=  tt1*tt1*tt1/2.7d1-tt1*tt2/6.0d0+tt3*0.5d0
            rtmp=  alp2/((alp1+1.0d-30)*sqrt(alp1+1.0d-30))
            if(rtmp .gt. 1.0d0) then
                dsig=  0.0d0
            elseif(rtmp .lt. -1.0d0) then
                dsig=  0.0d0
            else
                alp3    =  acos(rtmp)*R13
                alp12   =  2.0d0*sqrt(alp1)
                sig1    =  sqrt(max(tt1*R13+alp12*cos(alp3), 0.0d0))
                sig2    =  sqrt(max(tt1*R13-alp12*cos(pi*R13+alp3), 0.0d0))
                sig3    =  sqrt(max(tt1*R13-alp12*cos(pi*R13-alp3), 0.0d0))
                if((sig1 .lt. sig2) .or. (sig2 .lt. sig3)) then
                    print*, 'Error: wrong sigma for LES sigma:',sig1,sig2,sig3
                    stop
                end if
                if(sig1 .le. 1.0d-10) then
                    dsig=  0.0d0
                else
                    dsig=  sig3*(sig1-sig2)*(sig2-sig3)/(sig1*sig1)
                end if
            end if
            S_LES   =  dsig
            S_RANS  =  oo

            if(L_RANS .gt. L_LES) then
                rtmp=  1.0d0
            else
                rtmp=  0.0d0
            end if
            oo  =  S_RANS-fd*rtmp*(S_RANS-67.8d0*S_LES)
            st  =  oo+nut*fv2*k2_1*d2_1

            if(is_get_DDES_quality) fv(isec)%DDES_quality(3,iele)   =  1.0d0-fd

            r   =  min(nut*k2_1*d2_1/st, 1.0d1)
            g   =  r+0.3d0*(r**6-r)
            fw  =  g*(6.5d1/(g**6+6.4d1))**(1.0d0/6.0d0)
            pp  =  cb1*st*nut+0.933d0*gran
            dd  =  cw1*fw*d2_1*nut**2
            fv(isec)%rhs(1,iele)= -vol*(pp-dd)

            if(is_LHS) then
                chi_nut =  1.0d0/nul
                fv1_nut =  1073.733d0*chi_nut*chi*chi/(chi**3+3.57911d2)**2
                fv2_nut = -chi_nut/(1.0d0+chi*fv1)+chi*(fv1*chi_nut+chi*fv1_nut) &
                        & /(1.0d0+chi*fv1)**2
                if(r .lt. 1.0d1) then
                    r_nut   =  r/nut-r*r*(fv2/nut+fv2_nut)
                else
                    r_nut   =  0.0d0
                end if
                g_nut   = (1.0d0-cw2+6.0d0*cw2*r**5)*r_nut
                fw_nut  = (6.5d1/(6.4d1+g**6))**(1.0d0/6.0d0)*g_nut* &
                        & (1.0d0-g**6/(g**6+6.4d1))
                g       =  2.0d0*cw1*fw*nut*d2_1+cw1*nut*nut*d2_1*max(fw_nut, 0.0d0)

                eL  =  iele+iA(isec)-1
                call prec_add_eleR(fv_rans_prec, eL, eL, (/vol*g/))
            end if
        end do
    end do

    return
    end subroutine fv_DDES_sigma
!   ----------------------------------------------------------------------------
!   DDES source term, IDDES model.
!   ----------------------------------------------------------------------------
    subroutine fv_DDES_IDDES(lev,is_LHS)
    use var_global, only: L_ref,uref
    implicit none
    logical(dpL),intent(in):: is_LHS
    integer(dpI),intent(in):: lev
    integer(dpI):: isec,iele
    real   (dpR):: nul,nut,g,ug(3,3),vol,dnw,rh,mul,mut,gran,oo,gls,fd,d2_1, &
                &  chi,fv1,fv2,st,r,fw,pp,dd,uij,psi,L_RANS,L_LES,L_IDDES, &
                &  a,fb,fe,fe1,fe2,ft,fl,fdt,chi_nut,fv1_nut,fv2_nut,r_nut, &
                &  g_nut,fw_nut,rtmp

    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(sec(isec)%is_ghost .or. sec(isec)%is_bnd)    cycle

        do iele=1,sec(isec)%n_ele
            vol =  sec(isec)%vol( iele)
            dnw =  sec(isec)%dnw( iele)
            rh  =  fv (isec)%u (1,iele)
            mul =  fv (isec)%mu(1,iele)
            mut =  fv (isec)%mu(2,iele)
            nul =  mul/rh
            nut =  fv (isec)%turb(1,iele)

            ug(1,1:3)   =  fv(isec)%gra(4 :6 ,iele)
            ug(2,1:3)   =  fv(isec)%gra(7 :9 ,iele)
            ug(3,1:3)   =  fv(isec)%gra(10:12,iele)
            gran=  fv(isec)%turb_gra(1,iele)**2+fv(isec)%turb_gra(2,iele)**2 &
                & +fv(isec)%turb_gra(3,iele)**2
            oo  =  sqrt((ug(1,2)-ug(2,1))**2+(ug(1,3)-ug(3,1))**2+(ug(2,3)-ug(3,2))**2)
            gls =  sec(isec)%gls(iele)
            uij =  sqrt(ug(1,1)*ug(1,1)+ug(1,2)*ug(1,2)+ug(1,3)*ug(1,3) &
                &      +ug(2,1)*ug(2,1)+ug(2,2)*ug(2,2)+ug(2,3)*ug(2,3) &
                &      +ug(3,1)*ug(3,1)+ug(3,2)*ug(3,2)+ug(3,3)*ug(3,3))
            uij =  max(uij, 1.0d-6*uref/L_ref)
            chi =  nut/nul
            rtmp=  chi**3
            fv1 =  rtmp/(rtmp+3.57911d2)
            fv2 =  1.0d0-chi/(1.0d0+chi*fv1)
            Psi =  sqrt(min(1.0d2, (1.0d0-cb1*fv2/(0.424d0*cw1*kappa*kappa))/fv1))

            L_RANS  =  dnw
            rtmp    =  min(max(sec(isec)%gls_wn(iele), 0.15d0*max(dnw, gls)), gls)
            L_LES   =  Psi*C_DES*rtmp
            a       =  0.25d0-dnw/gls
            fB      =  min(1.0d0, 2.0d0*exp(-9.0d0*a*a))
            if(a .ge. 0.0d0) then
                fe1 =  2.0d0*exp(-11.09d0*a*a)
            else
                fe1 =  2.0d0*exp(-9.0d0*a*a)
            end if
            rtmp    =  nut/(kappa**2*dnw**2*max(uij, 1.0d-10*uref/L_ref))
            ft      =  tanh((1.63d0**2*rtmp)**3)
            fdt     =  1.0d0-tanh((8.0d0*rtmp)**3)
            rtmp    =  nul/(kappa**2*dnw**2*max(uij, 1.0d-10*uref/L_ref))
            fl      =  tanh((3.55d0**2*rtmp)**10)
            fe2     =  1.0d0-max(ft, fl)
            fe      =  max(0.0d0, fe1-1.0d0)*Psi*fe2
!           L_WMLES =  fB*(1.0d0+fe)*L_RANS+(1.0d0-fB)*L_LES
            fd      =  max(fB, 1.0d0-fdt)
            L_IDDES =  fd*(1.0d0+fe)*L_RANS+(1.0d0-fd)*L_LES
            if(is_get_DDES_quality) then
                fv(isec)%DDES_quality(3,iele)   =  fd*(1.0d0+fe)
                fv(isec)%DDES_quality(4,iele)   =  1.0d0-fd
            end if

            d2_1=  1.0d0/(L_IDDES*L_IDDES)
            st  =  oo+nut*fv2*k2_1*d2_1
            r   =  min(nut*k2_1*d2_1/st, 1.0d1)
            g   =  r+0.3d0*(r**6-r)
            fw  =  g*(6.5d1/(g**6+6.4d1))**(1.0d0/6.0d0)
            pp  =  cb1*st*nut+0.933d0*gran
            dd  =  cw1*fw*d2_1*nut**2
            fv(isec)%rhs(1,iele)= -vol*(pp-dd)

            if(is_LHS) then
                chi_nut =  1.0d0/nul
                fv1_nut =  1073.733d0*chi_nut*chi*chi/(chi**3+3.57911d2)**2
                fv2_nut = -chi_nut/(1.0d0+chi*fv1)+chi*(fv1*chi_nut+chi*fv1_nut) &
                        & /(1.0d0+chi*fv1)**2
                if(r .lt. 1.0d1) then
                    r_nut   =  r/nut-r*r*(fv2/nut+fv2_nut)
                else
                    r_nut   =  0.0d0
                end if
                g_nut   = (1.0d0-cw2+6.0d0*cw2*r**5)*r_nut
                fw_nut  = (6.5d1/(6.4d1+g**6))**(1.0d0/6.0d0)*g_nut* &
                        & (1.0d0-g**6/(g**6+6.4d1))
                g       = -cb1*nut*nut*k2_1*d2_1*min(fv2_nut, 0.0d0) &
                        & +2.0d0*cw1*fw*nut*d2_1+cw1*nut*nut*d2_1*max(fw_nut, 0.0d0)

                eL  =  iele+iA(isec)-1
                call prec_add_eleR(fv_rans_prec, eL, eL, (/vol*g/))
            end if
        end do
    end do

    return
    end subroutine fv_DDES_IDDES
    end subroutine fv_DDES_get_src
!-------------------------------------------------------------------------------
!   get the DDES quality.
!-------------------------------------------------------------------------------
    subroutine fv_get_DDES_quality(lev)
    use var_kind_def
    use var_fv
    use var_global, only: mesh_name
    use var_mesh
    use var_turb
    use var_uns_cal, only: uns_iter,uns_ra_begin,max_uns_iter
    implicit none
    integer(dpI),intent(in):: lev
    integer(dpI):: isec,iele
    real   (dpR):: nut,g(9),s12,s13,s23,eps,ss,k_modeled

    if(.not. is_get_DDES_quality)   return
    eps =  1.0d0/real(1+uns_iter-uns_ra_begin, dpR)
    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(.not. sec(isec)%is_int)  cycle

        do iele=1,sec(isec)%n_ele
            nut     =  fv(isec)%mu(2,iele)/fv(isec)%u(1,iele)
            g(1:9)  =  fv(isec)%gra(4:12,iele)
            s12 =  0.5d0*(g(2)+g(4))
            s13 =  0.5d0*(g(3)+g(7))
            s23 =  0.5d0*(g(6)+g(8))
            ss  =  sqrt(2.0d0*(g(1)**2+g(5)**2+g(9)**2)+4.0d0*(s12*s12+s13*s13+s23*s23))

!           kinetic energy modeled.
!           k_modeled   =  nut*ss/0.3d0
            k_modeled   = (nut/(0.09137d0*sec(isec)%gls(iele)))**2

            fv(isec)%DDES_quality(1,iele)   =  fv(isec)%DDES_quality(1,iele) &
                & +(k_modeled-fv(isec)%DDES_quality(1,iele))*eps

            if(uns_iter .eq. max_uns_iter) then
!               kinetic energy resolved.
                fv(isec)%DDES_quality(2,iele)   =  0.5d0*(fv(isec)%stress(1,iele) &
                    & +fv(isec)%stress(4,iele)+fv(isec)%stress(6,iele))
            end if
        end do
    end do
    if(uns_iter .eq. max_uns_iter)  &
        &  call fv_wr_cgns_cc(3, trim(adjustl(mesh_name))//'_DDES_quality.cgns')

    return
    end subroutine fv_get_DDES_quality
