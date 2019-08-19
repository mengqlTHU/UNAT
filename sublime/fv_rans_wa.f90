!-------------------------------------------------------------------------------
!   module for Wray-Agarwal one-equation model.
!-------------------------------------------------------------------------------
    module var_wa
        use var_kind_def
        implicit none

        logical(dpL),parameter:: is_WA_log  =  .false.
        integer(dpI),parameter:: WA_version =  2018
        real   (dpR),parameter:: C_b        =  1.66d0

        real   (dpR),parameter:: C1_ko      =  0.144d0
!       real   (dpR),parameter:: C1_ko      =  0.0833d0
        real   (dpR),parameter:: C1_ke      =  0.1127d0
        real   (dpR),parameter:: s_ko       =  0.5d0
        real   (dpR),parameter:: s_ke       =  1.0d0
        real   (dpR),parameter:: C_w        =  13.0d0
        real   (dpR),parameter:: C_w3       =  C_w**3
        real   (dpR),parameter:: C2_ko      =  C1_ko/0.41d0**2+s_ko
        real   (dpR),parameter:: C2_ke      =  C1_ke/0.41d0**2+s_ke

        real   (dpR),parameter:: c1_ko_2017 =  0.0829d0
        real   (dpR),parameter:: c1_ke_2017 =  0.1127d0
        real   (dpR),parameter:: s_ko_2017  =  0.72d0
        real   (dpR),parameter:: s_ke_2017  =  1.0d0
        real   (dpR),parameter:: C_w_2017   =  8.54d0
        real   (dpR),parameter:: C_w3_2017  =  C_w_2017**3
        real   (dpR),parameter:: C2_ko_2017 =  C1_ko_2017/0.41d0**2+s_ko_2017
        real   (dpR),parameter:: C2_ke_2017 =  C1_ke_2017/0.41d0**2+s_ke_2017
        real   (dpR),parameter:: C_m_2017   =  8.0d0

        real   (dpR),parameter:: c1_ko_2018 =  0.0829d0
        real   (dpR),parameter:: c1_ke_2018 =  0.1284d0
        real   (dpR),parameter:: s_ko_2018  =  0.72d0
        real   (dpR),parameter:: s_ke_2018  =  1.0d0
        real   (dpR),parameter:: C_w_2018   =  8.54d0
        real   (dpR),parameter:: C_w3_2018  =  C_w_2018**3
        real   (dpR),parameter:: C2_ko_2018 =  C1_ko_2018/0.41d0**2+s_ko_2018
        real   (dpR),parameter:: C2_ke_2018 =  C1_ke_2018/0.41d0**2+s_ke_2018
        real   (dpR),parameter:: C_m_2018   =  8.0d0
    end module var_wa
!-------------------------------------------------------------------------------
!   get mut from nut.
!-------------------------------------------------------------------------------
    subroutine fv_WA_nut_to_mut(lev,mode)
    use var_kind_def
    use var_fv
    use var_mesh, only: mesh,sec
    use var_wa
    implicit none
    integer(dpI),intent(in):: lev,mode
    integer(dpI):: isec,i,ele1,ele0
    real   (dpR):: rh,nut,nu,chi3,fv1

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

        if(WA_version .eq. 2015) then
            do i=ele1,ele0
                rh  =  fv(isec)%u   (1,i)
                nut =  fv(isec)%turb(1,i)
                nu  =  fv(isec)%mu  (1,i)/rh
                chi3= (nut/nu)**3
                fv1 =  chi3/(chi3+C_w3)
                fv(isec)%mu(2,i)=  rh*nut*fv1
            end do
        elseif(WA_version .eq. 2017) then
            do i=ele1,ele0
                rh  =  fv(isec)%u   (1,i)
                nut =  fv(isec)%turb(1,i)
                nu  =  fv(isec)%mu  (1,i)/rh
                chi3= (nut/nu)**3
                fv1 =  chi3/(chi3+C_w3_2017)
                fv(isec)%mu(2,i)=  rh*nut*fv1
            end do
        else
            do i=ele1,ele0
                rh  =  fv(isec)%u   (1,i)
                nut =  fv(isec)%turb(1,i)
                nu  =  fv(isec)%mu  (1,i)/rh
                chi3= (nut/nu)**3
                fv1 =  chi3/(chi3+C_w3_2018)
                fv(isec)%mu(2,i)=  rh*nut*fv1
            end do
        end if
    end do

    return
    end subroutine fv_WA_nut_to_mut
!-------------------------------------------------------------------------------
!   get f1, WA model.
!-------------------------------------------------------------------------------
    subroutine fv_wa_get_f1(lev)
    use var_kind_def
    use var_fv
    use var_mesh
    use var_parallel
    use var_slv, p=>ls_weight_order
    use var_WA, only: is_WA_log,C_b,WA_version
    implicit none
    integer(dpI),intent(in):: lev
    integer(dpI):: isec,iele,sL,eL,sR,eR,im,ss,ee,j
    real   (dpR):: g(9),s12,s13,s23,s,dRS,a,nu,n(3),vL,vR,LHS(3,3),dp(3),d,c(3),s1,d1,w

    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        fv(isec)%rhs=  0.0d0

        do iele=1,sec(isec)%n_ele
            g(1:9)  =  fv(isec)%gra(4:12,iele)
            s12 =  0.5d0*(g(2)+g(4))
            s13 =  0.5d0*(g(3)+g(7))
            s23 =  0.5d0*(g(6)+g(8))
            s   =  sqrt(2.0d0*(g(1)**2+g(5)**2+g(9)**2)+4.0d0*(s12*s12+s13*s13+s23*s23))
            fv(isec)%rhs(1,iele)=  s
            if(sec(isec)%is_bnd)    cycle

            nu  =  fv(isec)%mu(1,iele)/fv(isec)%u(1,iele)
            if(WA_version .eq. 2018) then
                w   =  sqrt((g(2)-g(4))**2+(g(3)-g(7))**2+(g(6)-g(8))**2)
                a   =  S*max(1.0d0, abs(w/s))
                vL  =  fv(isec)%mu(2,iele)/fv(isec)%u(1,iele)*s/0.3d0
                vR  =  s/0.3d0
                a   =  0.5d0*(nu+fv(isec)%turb(1,iele))*a*a/(0.3d0*vL*vR)
                fv(isec)%rhs(2,iele)=  tanh(a**4)
            elseif(.false.) then
                d1  =  1.0d0/sec(isec)%dnw(iele)
                a   =  min(C_b*fv(isec)%turb(1,iele)*d1*d1/(S*0.41d0**2), &
                    & (fv(isec)%turb(1,iele)/nu+1.0d0)**2)
!               a   =  C_b*(fv(isec)%turb(1,iele)+nu)*d1*d1/(S*0.41d0**2)
!               a   =  C_b*(fv(isec)%turb(1,iele)+0.0d0)*d1*d1/(S*0.41d0**2)
                fv(isec)%rhs(2,iele)=  min(tanh(a**4), 0.9d0)
            else
                dRS =  sec(isec)%dnw(iele)*sqrt(fv(isec)%turb(1,iele)*s)
                a   = (1.0d0+dRS/nu)/ &
                    & (1.0d0+(max(dRS, 1.5d0*fv(isec)%turb(1,iele))/(2.0d1*nu))**2)
                fv(isec)%rhs(2,iele)=  min(tanh(a**4), 0.9d0)
            end if
        end do
    end do

    if(gradient_method .eq. gradient_GG) then
        do im=1,mesh(lev)%n_mortar_b
            sL  =  mesh(lev)%mortar_LR(1,im)
            eL  =  mesh(lev)%mortar_LR(2,im)
            sR  =  mesh(lev)%mortar_LR(3,im)
            eR  =  mesh(lev)%mortar_LR(4,im)
            n(1:3)  =  mesh(lev)%mortar_n_vg(1:3,im)*mesh(lev)%mortar_n_vg(4,im)

            if(is_WA_log) then
                fv(sL)%rhs(3:5,eL)  =  fv(sL)%rhs(3:5,eL)+log(fv(sR)%rhs(1,eR))*n(1:3)
            else
                fv(sL)%rhs(3:5,eL)  =  fv(sL)%rhs(3:5,eL)+fv(sR)%rhs(1,eR)*n(1:3)
            end if
        end do

        do im=1+mesh(lev)%n_mortar_b,mesh(lev)%n_mortar
            sL  =  mesh(lev)%mortar_LR(1,im)
            eL  =  mesh(lev)%mortar_LR(2,im)
            sR  =  mesh(lev)%mortar_LR(3,im)
            eR  =  mesh(lev)%mortar_LR(4,im)
            n(1:3)  =  mesh(lev)%mortar_n_vg(1:3,im)*mesh(lev)%mortar_n_vg(4,im)

            vL  =  1.0d0/sec(sL)%vol(eL)
            vR  =  1.0d0/sec(sR)%vol(eR)
            s   = (vL*fv(sL)%rhs(1,eL)+vR*fv(sR)%rhs(1,eR))/(vL+vR)
            if(is_WA_log)   s   =  log(s)
            fv(sL)%rhs(3:5,eL)  =  fv(sL)%rhs(3:5,eL)+s*n(1:3)
            if(sec(sR)%is_int) then
            fv(sR)%rhs(3:5,eR)  =  fv(sR)%rhs(3:5,eR)-s*n(1:3)
            end if
        end do

        do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
            if(.not. sec(isec)%is_int)  cycle
            do iele=1,sec(isec)%n_ele
                fv(isec)%rhs(3:5,iele)  =  fv(isec)%rhs(3:5,iele)/sec(isec)%vol(iele)
            end do
        end do
    elseif(gradient_method .eq. gradient_LS) then
        do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
            if(.not. sec(isec)%is_int)  cycle
            do iele=1,sec(isec)%n_ele
                if(is_WA_log) then
                    s   =  log(fv(isec)%rhs(1,iele))
                else
                    s   =  fv(isec)%rhs(1,iele)
                end if
                do j=sec(isec)%iA_ls(iele),sec(isec)%iA_ls(iele+1)-1
                    ss  =  sec(isec)%jA_ls(1,j)
                    ee  =  sec(isec)%jA_ls(2,j)
                    if(is_WA_log) then
                        s1  =  log(fv(ss)%rhs(1,ee))
                    else
                        s1  =  fv(ss)%rhs(1,ee)
                    end if
                    fv(isec)%rhs(3:5,iele)  =  fv(isec)%rhs(3:5,iele) &
                        & +(s1-s)*sec(isec)%coe_ls(1:3,j)
                end do
            end do
        end do
    else
        do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
            if(.not. sec(isec)%is_int)  cycle
            do iele=1,sec(isec)%n_ele
                LHS(1,1)=  sec(isec)%coe_ls(1,iele)
                LHS(2,1)=  sec(isec)%coe_ls(2,iele)
                LHS(3,1)=  sec(isec)%coe_ls(3,iele)
                LHS(1,2)=  LHS(2,1)
                LHS(2,2)=  sec(isec)%coe_ls(4,iele)
                LHS(3,2)=  sec(isec)%coe_ls(5,iele)
                LHS(1,3)=  LHS(3,1)
                LHS(2,3)=  LHS(3,2)
                LHS(3,3)=  sec(isec)%coe_ls(6,iele)

                if(is_WA_log) then
                    s   =  log(fv(isec)%rhs(1,iele))
                else
                    s   =  fv(isec)%rhs(1,iele)
                end if
                do j=sec(isec)%iA_ls(iele),sec(isec)%iA_ls(iele+1)-1
                    ss      =  sec(isec)%jA_ls(1,j)
                    ee      =  sec(isec)%jA_ls(2,j)
                    dp(1:3) =  sec(ss)%cen(1:3,ee)-sec(isec)%cen(1:3,iele)
                    d       =  1.0d0/sqrt(dp(1)*dp(1)+dp(2)*dp(2)+dp(3)*dp(3))**p
                    c(1:3)  =  d*(LHS(1:3,1)*dp(1)+LHS(1:3,2)*dp(2)+LHS(1:3,3)*dp(3))

                    if(is_WA_log) then
                        s1  =  log(fv(ss)%rhs(1,ee))
                    else
                        s1  =  fv(ss)%rhs(1,ee)
                    end if
                    fv(isec)%rhs(3:5,iele)  =  fv(isec)%rhs(3:5,iele)+(s1-s)*c(1:3)
                end do
            end do
        end do
    end if

    return
    end subroutine fv_wa_get_f1
!-------------------------------------------------------------------------------
!   cal Rhs, WA turbulence model.
!-------------------------------------------------------------------------------
    subroutine fv_WA_get_rhs(lev,is_LHS)
    use var_kind_def
    use var_fv
    use var_global, only: n_dim,L_ref,uref
    use var_mesh
    use var_parallel
    use var_prec
    use var_turb
    use var_uns_cal, only: is_BDF_now,dt_uns,uns_iter,is_uns_initialized
    use var_wa
    implicit none
    logical(dpL),intent(in):: is_LHS
    integer(dpI),intent(in):: lev
    integer(dpI):: isec,iele,im,sL,eL,sR,eR,iA(1000),bct,s_donor,e_donor,per(3)
    real   (dpR):: dt,vol,s

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

    call fv_wa_get_f1(lev)

    if(WA_version .eq. 2015) then
        call fv_WA2015_get_rhs(lev, is_LHS)
    elseif(WA_version .eq. 2017) then
        call fv_WA2017_get_rhs(lev, is_LHS)
    else
        call fv_WA2018_get_rhs(lev, is_LHS)
    end if

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
    contains
!   ----------------------------------------------------------------------------
!   cal Rhs, WA turbulence model, version 2017.
!   ----------------------------------------------------------------------------
    subroutine fv_WA2017_get_rhs(lev,is_LHS)
    implicit none
    logical(dpL),intent(in):: is_LHS
    integer(dpI),intent(in):: lev
    real   (dpR):: n(5),uLf(5),uRf(5),tL,tR,nul,nut,gra(3),d(3),g,R,src,f,dt, &
                &  nne(3),unL,unR,vol,RS,fLtL(1),fLtR(1),fRtL(1),fRtR(1),vg(3), &
                &  f1,s,gR(3),gS(3),SS,R_S,C1_R,dnw,nu,RR,C1,rtmp

!   ----------------------------------------------------------------------------
!   source term.
    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(sec(isec)%is_ghost .or. sec(isec)%is_bnd)    cycle

        do iele=1,sec(isec)%n_ele
            vol     =  sec(isec)%vol     (    iele)
            dnw     =  sec(isec)%dnw     (    iele)
            R       =  fv (isec)%turb    (1  ,iele)
            s       =  fv (isec)%rhs     (1  ,iele)
            f1      =  fv (isec)%rhs     (2  ,iele)
            gs(1:3) =  fv (isec)%rhs     (3:5,iele)
            gR(1:3) =  fv (isec)%turb_gra(1:3,iele)
            nu      =  fv (isec)%mu(1,iele)/fv(isec)%u(1,iele)

            C1  =  f1*(C1_ko_2017-C1_ke_2017)+C1_ke_2017
            R_S =  R/max(S, 1.0d-14*uref/L_ref)
            RS  =  gR(1)*gs(1)+gR(2)*gs(2)+gR(3)*gs(3)
            SS  =  gS(1)*gs(1)+gS(2)*gs(2)+gS(3)*gs(3)
            RR  =  gR(1)*gR(1)+gR(2)*gR(2)+gR(3)*gR(3)
            if(is_WA_log) then
                src =  C1*R*S+f1*C2_ko_2017*R  *RS &
                    & -(1.0d0-f1)*min(C2_ke_2017*SS*R  **2, C_m_2017*RR)
            else
                src =  C1*R*S+f1*C2_ko_2017*R_S*RS &
                    & -(1.0d0-f1)*min(C2_ke_2017*SS*R_S**2, C_m_2017*RR)
            end if
            fv(isec)%rhs(1,iele)= -vol*src

            if(is_LHS) then
                C1_R=  0.0d0
                if(is_WA_log) then
                    g   =  R*S*min(C1_R, 0.0d0)+min(f1*C2_ko_2017*RS  , 0.0d0) &
                        & -2.0d0*(1.0d0-f1)*C2_ke_2017*R*SS
                else
                    g   =  R*S*min(C1_R, 0.0d0)+min(f1*C2_ko_2017*RS/S, 0.0d0) &
                        & -2.0d0*(1.0d0-f1)*C2_ke_2017*R*SS/max(S, 1.0d-10*uref/L_ref)**2
                end if

                eL  =  iele+iA(isec)-1
                call prec_add_eleR(fv_rans_prec, eL, eL, (/-vol*g/))
            end if
        end do
    end do
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
            bct     =  sec(sR)%bct
            uLf(1:5)=  fv(sR)%u(1:5,eR)
            uRf(1:5)=  fv(sR)%u(1:5,eR)
            tL      =  fv(sR)%turb(1,eR)
            tR      =  tL
            call DCOPY(3, fv(sR)%turb_gra(1,eR), 1, gra, 1)
            nul     =  fv(sR)%mu(1,eR)/uRf(1)
            nut     =  tL
            f1      =  fv(sL)%rhs(2,eL)
        else
            uLf(1:5)=  fv(sL)%u(1:5,eL)
            uRf(1:5)=  fv(sR)%u(1:5,eR)
            tL      =  fv(sL)%turb(1,eL)
            tR      =  fv(sR)%turb(1,eR)
            nul     =  0.5d0*(fv(sL)%mu(1,eL)/uLf(1)+fv(sR)%mu(1,eR)/uRf(1))
            nut     =  0.5d0*(tL+tR)
            f1      =  0.5d0*(fv(sL)%rhs(2,eL)+fv(sR)%rhs(2,eR))

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
        s   =  f1*(s_ko_2017-s_ke_2017)+s_ke_2017
        f   = (0.5d0*(unL+abs(unL))*tL+0.5d0*(unR-abs(unR))*tR &
            & -(n(1)*gra(1)+n(2)*gra(2)+n(3)*gra(3))*(nul+s*nut))*n(4)

        fv(sL)%rhs(1,eL)=  fv(sL)%rhs(1,eL)+f
        if(sec(sR)%is_int)  fv(sR)%rhs(1,eR)=  fv(sR)%rhs(1,eR)-f
        if(is_LHS) then
            fLtL=  n(4)*( 0.5d0*(unL+abs(unL))+n(5)*(nul+s*nut))
            fLtR=  n(4)*( 0.5d0*(unR-abs(unR))-n(5)*(nul+s*nut))
            fRtL=  n(4)*(-0.5d0*(unL+abs(unL))-n(5)*(nul+s*nut))
            fRtR=  n(4)*(-0.5d0*(unR-abs(unR))+n(5)*(nul+s*nut))

            eL  =  eL+iA(sL)-1
            eR  =  eR+iA(sR)-1
            call prec_add_eleR(fv_rans_prec, eL, eL, fLtL)
            if(sec(sR)%is_int) then
                call prec_add_eleR(fv_rans_prec, eL, eR, fLtR)
                call prec_add_eleR(fv_rans_prec, eR, eR, fRtR)
                call prec_add_eleR(fv_rans_prec, eR, eL, fRtL)
            elseif(sec(sR)%is_ghost) then
                call get_ID_ghost_ele(sR, eR, s_donor, e_donor, per)
                if(s_donor .gt. 0) then
                    eR  =  e_donor+iA(s_donor)-1
                    call prec_add_eleR(fv_rans_prec, eL, eR, fLtR)
                    call prec_add_eleR(fv_rans_prec, eR, eR, fRtR)
                    call prec_add_eleR(fv_rans_prec, eR, eL, fRtL)
                end if
            end if
        end if
    end do
!   convective and viscous terms.
!   ----------------------------------------------------------------------------

    return
    end subroutine fv_WA2017_get_rhs
!   ----------------------------------------------------------------------------
!   cal Rhs, WA turbulence model, version 2015.
!   ----------------------------------------------------------------------------
    subroutine fv_WA2015_get_rhs(lev,is_LHS)
    implicit none
    logical(dpL),intent(in):: is_LHS
    integer(dpI),intent(in):: lev
    real   (dpR):: n(5),uLf(5),uRf(5),tL,tR,nul,nut,gra(3),d(3),g,R,src,f,dt, &
                &  nne(3),unL,unR,vol,RS,fLtL(1),fLtR(1),fRtL(1),fRtR(1),vg(3), &
                &  f1,s,gR(3),gS(3),SS,R_S,C1_R,dnw,nu,RR,rtmp

!   ----------------------------------------------------------------------------
!   source term.
    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(sec(isec)%is_ghost .or. sec(isec)%is_bnd)    cycle

        do iele=1,sec(isec)%n_ele
            vol     =  sec(isec)%vol     (    iele)
            dnw     =  sec(isec)%dnw     (    iele)
            R       =  fv (isec)%turb    (1  ,iele)
            s       =  fv (isec)%rhs     (1  ,iele)
            f1      =  fv (isec)%rhs     (2  ,iele)
            gs(1:3) =  fv (isec)%rhs     (3:5,iele)
            gR(1:3) =  fv (isec)%turb_gra(1:3,iele)
            nu      =  fv (isec)%mu(1,iele)/fv(isec)%u(1,iele)

            R_S =  R/max(S, 1.0d-14*uref/L_ref)
            RS  =  gR(1)*gs(1)+gR(2)*gs(2)+gR(3)*gs(3)
            SS  =  gS(1)*gs(1)+gS(2)*gs(2)+gS(3)*gs(3)
            RR  =  gR(1)*gR(1)+gR(2)*gR(2)+gR(3)*gR(3)
            if(is_WA_log) then
                src =  C1_ko*R*S+f1*C2_ko*R  *RS-(1.0d0-f1)*C2_ke*SS*R  **2
            else
                src =  C1_ko*R*S+f1*C2_ko*R_S*RS-(1.0d0-f1)*C2_ke*SS*R_S**2
            end if
            fv(isec)%rhs(1,iele)= -vol*src

            if(is_LHS) then
                C1_R=  0.0d0
                if(is_WA_log) then
                    g   =  R*S*min(C1_R, 0.0d0)+min(f1*C2_ko*RS, 0.0d0) &
                        & -2.0d0*(1.0d0-f1)*C2_ke*R*SS
                else
                    g   =  R*S*min(C1_R, 0.0d0)+min(f1*C2_ko*RS/S, 0.0d0) &
                        & -2.0d0*(1.0d0-f1)*C2_ke*R*SS/max(S, 1.0d-10*uref/L_ref)**2
                end if

                eL  =  iele+iA(isec)-1
                call prec_add_eleR(fv_rans_prec, eL, eL, (/-vol*g/))
            end if
        end do
    end do
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
            bct     =  sec(sR)%bct
            uLf(1:5)=  fv(sR)%u(1:5,eR)
            uRf(1:5)=  fv(sR)%u(1:5,eR)
            tL      =  fv(sR)%turb(1,eR)
            tR      =  tL
            call DCOPY(3, fv(sR)%turb_gra(1,eR), 1, gra, 1)
            nul     =  fv(sR)%mu(1,eR)/uRf(1)
            nut     =  tL
            f1      =  fv(sL)%rhs(2,eL)
        else
            uLf(1:5)=  fv(sL)%u(1:5,eL)
            uRf(1:5)=  fv(sR)%u(1:5,eR)
            tL      =  fv(sL)%turb(1,eL)
            tR      =  fv(sR)%turb(1,eR)
            nul     =  0.5d0*(fv(sL)%mu(1,eL)/uLf(1)+fv(sR)%mu(1,eR)/uRf(1))
            nut     =  0.5d0*(tL+tR)
            f1      =  0.5d0*(fv(sL)%rhs(2,eL)+fv(sR)%rhs(2,eR))

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
        s   =  f1*(s_ko-s_ke)+s_ke
        f   = (0.5d0*(unL+abs(unL))*tL+0.5d0*(unR-abs(unR))*tR &
            & -(n(1)*gra(1)+n(2)*gra(2)+n(3)*gra(3))*(nul+s*nut))*n(4)

        fv(sL)%rhs(1,eL)=  fv(sL)%rhs(1,eL)+f
        if(sec(sR)%is_int)  fv(sR)%rhs(1,eR)=  fv(sR)%rhs(1,eR)-f
        if(is_LHS) then
            fLtL=  n(4)*( 0.5d0*(unL+abs(unL))+n(5)*(nul+s*nut))
            fLtR=  n(4)*( 0.5d0*(unR-abs(unR))-n(5)*(nul+s*nut))
            fRtL=  n(4)*(-0.5d0*(unL+abs(unL))-n(5)*(nul+s*nut))
            fRtR=  n(4)*(-0.5d0*(unR-abs(unR))+n(5)*(nul+s*nut))

            eL  =  eL+iA(sL)-1
            eR  =  eR+iA(sR)-1
            call prec_add_eleR(fv_rans_prec, eL, eL, fLtL)
            if(sec(sR)%is_int) then
                call prec_add_eleR(fv_rans_prec, eL, eR, fLtR)
                call prec_add_eleR(fv_rans_prec, eR, eR, fRtR)
                call prec_add_eleR(fv_rans_prec, eR, eL, fRtL)
            elseif(sec(sR)%is_ghost) then
                call get_ID_ghost_ele(sR, eR, s_donor, e_donor, per)
                if(s_donor .gt. 0) then
                    eR  =  e_donor+iA(s_donor)-1
                    call prec_add_eleR(fv_rans_prec, eL, eR, fLtR)
                    call prec_add_eleR(fv_rans_prec, eR, eR, fRtR)
                    call prec_add_eleR(fv_rans_prec, eR, eL, fRtL)
                end if
            end if
        end if
    end do
!   convective and viscous terms.
!   ----------------------------------------------------------------------------

    return
    end subroutine fv_WA2015_get_rhs
!   ----------------------------------------------------------------------------
!   cal Rhs, WA turbulence model, version 2018.
!   ----------------------------------------------------------------------------
    subroutine fv_WA2018_get_rhs(lev,is_LHS)
    implicit none
    logical(dpL),intent(in):: is_LHS
    integer(dpI),intent(in):: lev
    real   (dpR):: n(5),uLf(5),uRf(5),tL,tR,nul,nut,gra(3),d(3),g,R,src,f,dt, &
                &  nne(3),unL,unR,vol,RS,fLtL(1),fLtR(1),fRtL(1),fRtR(1),vg(3), &
                &  f1,s,gR(3),gS(3),SS,R_S,C1_R,dnw,nu,RR,C1,rtmp

!   ----------------------------------------------------------------------------
!   source term.
    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(sec(isec)%is_ghost .or. sec(isec)%is_bnd)    cycle

        do iele=1,sec(isec)%n_ele
            vol     =  sec(isec)%vol     (    iele)
            dnw     =  sec(isec)%dnw     (    iele)
            R       =  fv (isec)%turb    (1  ,iele)
            s       =  fv (isec)%rhs     (1  ,iele)
            f1      =  fv (isec)%rhs     (2  ,iele)
            gs(1:3) =  fv (isec)%rhs     (3:5,iele)
            gR(1:3) =  fv (isec)%turb_gra(1:3,iele)
            nu      =  fv (isec)%mu(1,iele)/fv(isec)%u(1,iele)

            C1  =  f1*(C1_ko_2018-C1_ke_2018)+C1_ke_2018
            R_S =  R/max(S, 1.0d-14*uref/L_ref)
            RS  =  gR(1)*gs(1)+gR(2)*gs(2)+gR(3)*gs(3)
            SS  =  gS(1)*gs(1)+gS(2)*gs(2)+gS(3)*gs(3)
            RR  =  gR(1)*gR(1)+gR(2)*gR(2)+gR(3)*gR(3)
            if(is_WA_log) then
                src =  C1*R*S+f1*C2_ko_2018*R  *RS &
                    & -(1.0d0-f1)*min(C2_ke_2018*SS*R  **2, C_m_2018*RR)
            else
                src =  C1*R*S+f1*C2_ko_2018*R_S*RS &
                    & -(1.0d0-f1)*min(C2_ke_2018*SS*R_S**2, C_m_2018*RR)
            end if
            fv(isec)%rhs(1,iele)= -vol*src

            if(is_LHS) then
                C1_R=  0.0d0
                if(is_WA_log) then
                    g   =  R*S*min(C1_R, 0.0d0)+min(f1*C2_ko_2018*RS  , 0.0d0) &
                        & -2.0d0*(1.0d0-f1)*C2_ke_2018*R*SS
                else
                    g   =  R*S*min(C1_R, 0.0d0)+min(f1*C2_ko_2018*RS/S, 0.0d0) &
                        & -2.0d0*(1.0d0-f1)*C2_ke_2018*R*SS/max(S, 1.0d-10*uref/L_ref)**2
                end if

                eL  =  iele+iA(isec)-1
                call prec_add_eleR(fv_rans_prec, eL, eL, (/-vol*g/))
            end if
        end do
    end do
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
            bct     =  sec(sR)%bct
            uLf(1:5)=  fv(sR)%u(1:5,eR)
            uRf(1:5)=  fv(sR)%u(1:5,eR)
            tL      =  fv(sR)%turb(1,eR)
            tR      =  tL
            call DCOPY(3, fv(sR)%turb_gra(1,eR), 1, gra, 1)
            nul     =  fv(sR)%mu(1,eR)/uRf(1)
            nut     =  tL
            f1      =  fv(sL)%rhs(2,eL)
        else
            uLf(1:5)=  fv(sL)%u(1:5,eL)
            uRf(1:5)=  fv(sR)%u(1:5,eR)
            tL      =  fv(sL)%turb(1,eL)
            tR      =  fv(sR)%turb(1,eR)
            nul     =  0.5d0*(fv(sL)%mu(1,eL)/uLf(1)+fv(sR)%mu(1,eR)/uRf(1))
            nut     =  0.5d0*(tL+tR)
            f1      =  0.5d0*(fv(sL)%rhs(2,eL)+fv(sR)%rhs(2,eR))

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
        s   =  f1*(s_ko_2018-s_ke_2018)+s_ke_2018
        f   = (0.5d0*(unL+abs(unL))*tL+0.5d0*(unR-abs(unR))*tR &
            & -(n(1)*gra(1)+n(2)*gra(2)+n(3)*gra(3))*(nul+s*nut))*n(4)

        fv(sL)%rhs(1,eL)=  fv(sL)%rhs(1,eL)+f
        if(sec(sR)%is_int)  fv(sR)%rhs(1,eR)=  fv(sR)%rhs(1,eR)-f
        if(is_LHS) then
            fLtL=  n(4)*( 0.5d0*(unL+abs(unL))+n(5)*(nul+s*nut))
            fLtR=  n(4)*( 0.5d0*(unR-abs(unR))-n(5)*(nul+s*nut))
            fRtL=  n(4)*(-0.5d0*(unL+abs(unL))-n(5)*(nul+s*nut))
            fRtR=  n(4)*(-0.5d0*(unR-abs(unR))+n(5)*(nul+s*nut))

            eL  =  eL+iA(sL)-1
            eR  =  eR+iA(sR)-1
            call prec_add_eleR(fv_rans_prec, eL, eL, fLtL)
            if(sec(sR)%is_int) then
                call prec_add_eleR(fv_rans_prec, eL, eR, fLtR)
                call prec_add_eleR(fv_rans_prec, eR, eR, fRtR)
                call prec_add_eleR(fv_rans_prec, eR, eL, fRtL)
            elseif(sec(sR)%is_ghost) then
                call get_ID_ghost_ele(sR, eR, s_donor, e_donor, per)
                if(s_donor .gt. 0) then
                    eR  =  e_donor+iA(s_donor)-1
                    call prec_add_eleR(fv_rans_prec, eL, eR, fLtR)
                    call prec_add_eleR(fv_rans_prec, eR, eR, fRtR)
                    call prec_add_eleR(fv_rans_prec, eR, eL, fRtL)
                end if
            end if
        end if
    end do
!   convective and viscous terms.
!   ----------------------------------------------------------------------------

    return
    end subroutine fv_WA2018_get_rhs
    end subroutine fv_WA_get_rhs
