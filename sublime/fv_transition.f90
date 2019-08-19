!-------------------------------------------------------------------------------
!   module AGS transition model.
!-------------------------------------------------------------------------------
    module var_ags
        use var_kind_def
        implicit none

        real(dpR),allocatable:: intermittency_line_L(:,:)
        real(dpR),allocatable:: intermittency_line_G(:,:)
        real(dpR),allocatable:: intermittency_wall(:),BL_thickness_wall(:)
    end module var_ags
!-------------------------------------------------------------------------------
!   get the intermittency with AGS model.
!-------------------------------------------------------------------------------
    subroutine fv_get_intermittency_ags
    use var_kind_def
    use var_ags
    use var_eline
    use var_fv
    use var_global, only: err_mem,R23
    use var_mesh
    use var_parallel, only: myid
    use var_wall
    implicit none
    logical(dpL):: ltmp
    integer(dpI):: isec,iele,i,j,L,R,im
    real   (dpR):: ux,uy,uz,vx,vy,vz,wx,wy,wz,k,Tu_L,omega,mu,s12,s13,s23,F,v(3), &
                &  Re_theta,Pohl,dnw,rh,u_norm,u_x,u_y,u_z,u_s,s,mut,Re_theta_s

    call fv_get_bl_parameter

!   ----------------------------------------------------------------------------
!   memory allocation.
    if(.not. allocated(intermittency_wall)) then
        allocate(BL_thickness_wall (n_ele_wall), stat=err_mem)
        allocate(intermittency_wall(n_ele_wall), stat=err_mem)

        j   =  0
        do i=1,n_eline
            if(eline(i)%mortar_L .gt. 0)    j   =  j+1
        end do
        allocate(intermittency_line_L(4,j), stat=err_mem)
        allocate(intermittency_line_G(4,j), stat=err_mem)
    end if
!   memory allocation.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   synchronize the intermittency of 1D line from AGS model.
    j   =  0
    do i=1,n_eline
        if(eline(i)%mortar_L .le. 0)    cycle
        j   =  j+1
        intermittency_line_L(1,j)   =  real(myid, 8)
        intermittency_line_L(2,j)   =  real(eline(i)%mortar_L, 8)
        intermittency_line_L(3,j)   =  eline(i)%BL_thickness
        intermittency_line_L(4,j)   =  eline(i)%intermittency
    end do
    intermittency_line_G=  intermittency_line_L
!   synchronize the intermittency of 1D line from AGS model.
!   ----------------------------------------------------------------------------

    intermittency_wall  =  0.0d0
    do i=1,size(intermittency_line_G, dim=2)
        myid=  nint(intermittency_line_G(1,i), dpI)
        call ib_search(1, n_ele_wall, 1, 2, ID_wall, myid, ltmp, L, R)
        if(.not. ltmp)  stop 'Error: AGS fails to find wall element.'
        im  =  nint(intermittency_line_G(2,i), dpI)
        call ib_search(L, R, 2, 2, ID_wall, im, ltmp, L, R)
        if((.not. ltmp) .or. (L .ne. R))    stop 'Error: AGS fails to find wall element.'
        BL_thickness_wall (L)   =  intermittency_line_G(3,i)
        intermittency_wall(L)   =  intermittency_line_G(4,i)
    end do

    do isec=mesh(0)%sec_1,mesh(0)%sec_0
        if(.not. sec(isec)%is_int)  cycle
        do iele=1,sec(isec)%n_ele
            i   =  fv(isec)%nearest_wall(iele)

            if(i .gt. 0) then
                if(sec(isec)%dnw(iele) .ge. BL_thickness_wall(i)) then
                    fv(isec)%intermittency(iele)=  1.0d0
                else
                    fv(isec)%intermittency(iele)=  intermittency_wall(i)
                end if
            else
                fv(isec)%intermittency(iele)=  1.0d0
            end if
        end do
    end do
!   call fv_wr_tec_cc(5)

    do isec=mesh(0)%sec_1,mesh(0)%sec_0
        if(.not. sec(isec)%is_int)  cycle

        do iele=1,sec(isec)%n_ele
            dnw =  sec(isec)%dnw(  iele)
            ux  =  fv(isec)%gra(4 ,iele)
            uy  =  fv(isec)%gra(5 ,iele)
            uz  =  fv(isec)%gra(6 ,iele)
            vx  =  fv(isec)%gra(7 ,iele)
            vy  =  fv(isec)%gra(8 ,iele)
            vz  =  fv(isec)%gra(9 ,iele)
            wx  =  fv(isec)%gra(10,iele)
            wy  =  fv(isec)%gra(11,iele)
            wz  =  fv(isec)%gra(12,iele)
            rh      =  fv(isec)%u (1  ,iele)
            v(1:3)  =  fv(isec)%u (2:4,iele)
            mu      =  fv(isec)%mu(1  ,iele)
            mut     =  fv(isec)%mu(2  ,iele)
            call norm_vec(3, v, u_norm)

            U_x = (v(1)*ux+v(2)*vx+v(3)*wx)/sqrt(u_norm)
            U_y = (v(1)*uy+v(2)*vy+v(3)*wy)/sqrt(u_norm)
            U_z = (v(1)*uz+v(2)*vz+v(3)*wz)/sqrt(u_norm)
            U_s = (v(1)*U_x+v(2)*U_y+v(3)*U_z)/u_norm
            Pohl=  dnw**2*U_s/(mu/rh)

            s12 =  0.5d0*(uy+vx)
            s13 =  0.5d0*(uz+wx)
            s23 =  0.5d0*(vz+wy)
            S   =  sqrt(2.0d0*(ux*ux+vy*vy+wz*wz)+4.0d0*(s12*s12+s13*s13+s23*s23))
            omega   =  S/0.3d0
            k       =  mut*omega/rh
            Tu_L    =  min(sqrt(R23*k)/max(omega*dnw, 0.0d0), 1.0d0)*1.0d2
            Re_theta=  rh*u_norm*dnw/mu

            if(Pohl .le. 0.0d0) then
                F   =  6.91d0+12.75d0*Pohl+63.34d0*Pohl**2
            else
                F   =  6.91d0+2.48d0*Pohl-12.27d0*Pohl**2
            end if
            Re_theta_s  =  1.63d2+exp(F-F/6.91d0*Tu_L*1.0d2)
            if(Re_theta .ge. Re_theta_s) then
                fv(isec)%intermittency  =  1.0d0
            else
                fv(isec)%intermittency  =  0.0d0
            end if
        end do
    end do

    return
    end subroutine fv_get_intermittency_ags
!-------------------------------------------------------------------------------
!   get the intermittency with PTM model.
!-------------------------------------------------------------------------------
    subroutine fv_get_intermittency_ptm
    use var_kind_def
    use var_air, only: gk
    use var_fv
    use var_global, only: R23,pi
    use var_mesh
    implicit none
    integer(dpI):: isec,iele
    real   (dpR):: ux,uy,uz,vx,vy,vz,wx,wy,wz,s,k,Tu_L,omega,mu,s12,s13,s23,F, &
                &  v(3),U,dnw,rh,mut,Re_v,ptm1,ptm2,F3,R_t,ptm,ps,M2

    do isec=mesh(0)%sec_1,mesh(0)%sec_0
        if(.not. sec(isec)%is_int)  cycle

        do iele=1,sec(isec)%n_ele
            dnw =  sec(isec)%dnw(  iele)
            ux  =  fv(isec)%gra(4 ,iele)
            uy  =  fv(isec)%gra(5 ,iele)
            uz  =  fv(isec)%gra(6 ,iele)
            vx  =  fv(isec)%gra(7 ,iele)
            vy  =  fv(isec)%gra(8 ,iele)
            vz  =  fv(isec)%gra(9 ,iele)
            wx  =  fv(isec)%gra(10,iele)
            wy  =  fv(isec)%gra(11,iele)
            wz  =  fv(isec)%gra(12,iele)
            rh      =  fv(isec)%u (1  ,iele)
            v(1:3)  =  fv(isec)%u (2:4,iele)
            mu      =  fv(isec)%mu(1  ,iele)
            mut     =  fv(isec)%mu(2  ,iele)
            call norm_vec(3, v, U)

            s12 =  0.5d0*(uy+vx)
            s13 =  0.5d0*(uz+wx)
            s23 =  0.5d0*(vz+wy)
            S   =  sqrt(2.0d0*(ux*ux+vy*vy+wz*wz)+4.0d0*(s12*s12+s13*s13+s23*s23))
            omega   =  S/0.3d0
            k       =  mut*omega/rh
            Tu_L    =  min(sqrt(R23*k)/max(omega*dnw, 0.0d0), 1.0d0)*1.0d2
            Re_v    =  rh*dnw*dnw*S/mu
            R_t     =  mut/mu

            if(Re_v .le. 1.0d3) then
                ptm1=  1.0d0-(3.28d-4*Re_v-3.94d-7*Re_v**2+1.43d-10*Re_v**3)
            else
                ptm1=  1.0d0-(0.12d0+1.0d-5*Re_v)
            end if

            ps      =  fv(isec)%gra(13,iele)*v(1)+fv(isec)%gra(14,iele)*v(2) &
                    & +fv(isec)%gra(15,iele)*v(3)
            M2      =  U*U/(gk*fv(isec)%u(5,iele)/rh)
            K       = -mu*abs(1.0d0-M2)*ps/(rh**2*U**3)
            if(K .lt. 0.0d0)  then
                ptm2= -abs(K)**0.4d0*Re_v*0.0125d0
            else
                ptm2=  0.0d0
            end if

            if(.false.) then
                F3  =  exp(-(R_t/6.5d0)**4)
                ptm =  1.0d0-ptm1*F3
            elseif(.false.) then
                F   =  2.5d0*sqrt(2.0d0*pi)*exp(-0.5d0*(R_t-3.0d0)**2)
                F3  =  exp(-(R_t/3.0d0)**2)*(1.0d0-F)+0.5d0*F
                ptm =  1.0d0-0.94d0*ptm1*F3
            else
                F3  =  exp(-(R_t/5.0d0)**4)
                ptm =  max(0.0d0, min(1.0d0, 1.0d0-0.94d0*ptm1*F3))
                ptm =  1.0d0-0.94d0*F3*(ptm1+ptm2)
            end if

            fv(isec)%intermittency(iele)=  max(0.0d0, min(1.0d0, ptm))
        end do
    end do

    return
    end subroutine fv_get_intermittency_ptm
!-------------------------------------------------------------------------------
!   get the intermittency with BC model.
!-------------------------------------------------------------------------------
    subroutine fv_get_intermittency_BC
    use var_kind_def
    use var_bndv, only: tu_fs
    use var_fv
    use var_global, only: uref
    use var_mesh
    implicit none
    integer(dpI):: isec,iele
    real   (dpR):: ux,uy,uz,vx,vy,vz,wx,wy,wz,rh,mu,mut,oo,Re_v,Re_t,Re_c,T1,T2, &
                &  U,nu_b,dnw

    do isec=mesh(0)%sec_1,mesh(0)%sec_0
        if(.not. sec(isec)%is_int)  cycle

        do iele=1,sec(isec)%n_ele
            dnw =  sec(isec)%dnw(  iele)
            ux  =  fv(isec)%gra(4 ,iele)
            uy  =  fv(isec)%gra(5 ,iele)
            uz  =  fv(isec)%gra(6 ,iele)
            vx  =  fv(isec)%gra(7 ,iele)
            vy  =  fv(isec)%gra(8 ,iele)
            vz  =  fv(isec)%gra(9 ,iele)
            wx  =  fv(isec)%gra(10,iele)
            wy  =  fv(isec)%gra(11,iele)
            wz  =  fv(isec)%gra(12,iele)
            rh  =  fv(isec)%u  (1 ,iele)
            mu  =  fv(isec)%mu (1 ,iele)
            mut =  fv(isec)%mu (2 ,iele)

            oo  =  sqrt((uy-vx)**2+(uz-wx)**2+(vz-wy)**2)
            Re_v=  rh*dnw*dnw*oo/mu
            Re_t=  Re_v/2.193d0
            Re_c=  803.73d0*(Tu_fs*1.0d2+0.6067d0)**(-1.027d0)
            T1  =  max(Re_t-Re_c, 0.0d0)/(2.0d-3*Re_c)
            U   =  sqrt(fv(isec)%u(2,iele)**2+fv(isec)%u(3,iele)**2+fv(isec)%u(4,iele)**2)
            nu_b=  mut/(rh*dnw*max(U, 1.0d-6*uref))
            T2  =  max(nu_b-5.0d0, 0.0d0)/5.0d0

            fv(isec)%intermittency(iele)=  1.0d0-exp(-sqrt(T1)-sqrt(T2))
        end do
    end do

    return
    end subroutine fv_get_intermittency_BC
!-------------------------------------------------------------------------------
!   cal Rhs, SA with the one-equation Menter transition model.
!-------------------------------------------------------------------------------
    subroutine fv_SAM_get_rhs(lev,is_LHS)
    use var_kind_def
    use var_bndv, only: Tu_fs
    use var_fv
    use var_global, only: R23
    use var_mesh
    use var_prec
    use var_turb
    use var_uns_cal, only: is_BDF_now,dt_uns,uns_iter,is_uns_initialized
    implicit none
    logical(dpL),intent(in):: is_LHS
    integer(dpI),intent(in):: lev
    integer(dpI):: isec,iele,im,sL,eL,sR,eR,iA(1000)
    real   (dpR):: n(5),uLf(5),uRf(5),tL(2),tR(2),nul,nut,gra(6),fL(2),fR(2),t1,t2, &
                &  d(2,2),g,unL,unR,vol,dnw,rh,mul,mut,ux,uy,uz,vx,vy,vz,wx,wy, &
                &  wz,gran,oo,d2_1,chi,fv1,fv2,s,st,r,fw,pp,dd,chi_nut,fv1_nut, &
                &  fv2_nut,r_nut,g_nut,fw_nut,fLtL(2,2),fLtR(2,2),fRtL(2,2),Pohl, &
                &  fRtR(2,2),vg(3),u(3),spL,spR,rV(2),s12,s13,s23,Re_v,R_T,omega, &
                &  k,Tu_L,u_norm,U_x,U_y,U_z,U_s,lambda_L,F_PG,Re_thetac,F_onset1, &
                &  F_onset2,F_onset3,F_onset,F_turb,P_gamma,intermittency,D_gamma,rtmp

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
    D   =  0.0d0
    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(sec(isec)%is_ghost .or. sec(isec)%is_bnd)    cycle

        do iele=1,sec(isec)%n_ele
            vol             =  sec(isec)%vol( iele)
            dnw             =  sec(isec)%dnw( iele)
            rh              =  fv (isec)%u (1,iele)
            u(1:3)          =  fv (isec)%u (2:4,iele)
            mul             =  fv (isec)%mu(1,iele)
            mut             =  fv (isec)%mu(2,iele)
            nul             =  mul/rh
            nut             =  fv (isec)%turb(1,iele)
            intermittency   =  fv (isec)%turb(2,iele)

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

            s12 =  0.5d0*(uy+vx)
            s13 =  0.5d0*(uz+wx)
            s23 =  0.5d0*(vz+wy)
            S   =  sqrt(2.0d0*(ux*ux+vy*vy+wz*wz)+4.0d0*(s12*s12+s13*s13+s23*s23))
            d2_1=  1.0d0/(dnw*dnw)

            chi =  nut/nul
            rtmp=  chi**3
            fv1 =  rtmp/(rtmp+3.57911d2)
            fv2 =  1.0d0-chi/(1.0d0+chi*fv1)
            st  =  max(oo+nut*k2_1*d2_1*fv2, 0.3d0*oo)
            st  =  oo+nut*k2_1*d2_1*fv2
            r   =  min(nut*k2_1*d2_1/st, 1.0d1)
            g   =  r+0.3d0*(r**6-r)
            fw  =  g*(6.5d1/(g**6+6.4d1))**(1.0d0/6.0d0)
            t1  =  max(0.0d0, min(1.0d0, intermittency))
            t2  =  t1
            pp  =  t1*cb1*st*nut+0.933d0*gran
            dd  =  t2*cw1*fw*d2_1*nut**2

            Re_v    =  dnw*dnw*S/nul
            R_T     =  mut/mul
            F_turb  =  exp(-(0.5d0*R_T)**4)
            omega   =  S/0.3d0
            k       =  mut*omega/rh
!           Tu_L    =  min(sqrt(R23*k)/max(omega*dnw, ua_fs), 1.0d0)*1.0d2
            Tu_L    =  min(sqrt(R23*k)/max(omega*dnw, 0.0d0), 1.0d0)*1.0d2
            u_norm  =  sqrt(u(1)**2+u(2)**2+u(3)**2)
            U_x     = (u(1)*ux+u(2)*vx+u(3)*wx)/sqrt(u_norm)
            U_y     = (u(1)*uy+u(2)*vy+u(3)*wy)/sqrt(u_norm)
            U_z     = (u(1)*uz+u(2)*vz+u(3)*wz)/sqrt(u_norm)
            U_s     = (u(1)*U_x+u(2)*U_y+u(3)*U_z)/u_norm
            Pohl    =  U_s*dnw*dnw/nul

            if(.true.) then
                lambda_L=  min(max(-1.0d0, 7.57d-3*Pohl+1.28d-2), 1.0d0)
                if(lambda_L .ge. 0.0d0) then
                    F_PG=  min(1.0d0+14.68d0*lambda_L, 1.5d0)
                else
                    F_PG=  min(1.0d0-7.34d0*lambda_L, 3.0d0)
                end if
                F_PG        =  max(F_PG, 0.0d0)
                Re_thetac   =  1.0d2+1.0d3*exp(-Tu_L*F_PG)
            else
                omega       =  2.2d2*exp(-Tu_fs/9.0d-3)*(atan(Pohl/2.0d-2+0.842d0)-0.7d0)
                Re_thetac   =  1.045d3*exp(-Tu_fs/1.0d-2)+omega+1.55d2
            end if

            F_onset1    =  Re_v/(2.2d0*Re_thetac)
            F_onset2    =  min(F_onset1, 2.0d0)
            F_onset3    =  max(1.0d0-(R_T/3.5d0)**3, 0.0d0)
            F_onset     =  max(F_onset2-F_onset3, 0.0d0)
            P_gamma     =  1.0d2*rh*S*F_onset*intermittency*(1.0d0-intermittency)
            D_gamma     =  6.0d-2*rh*oo*F_turb*intermittency*(5.0d1*intermittency-1.0d0)

            fv(isec)%rhs(1,iele)= -vol*(pp     -dd     )
            fv(isec)%rhs(2,iele)= -vol*(P_gamma-D_gamma)

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
                D(1,1)  =  vol*g
                D(2,2)  = (1.0d2*S*F_onset*2.0d0*intermittency &
                        & +6.0d-2*oo*F_turb*1.0d2*intermittency)*vol
                call prec_add_eleR(fv_rans_prec, eL, eL, D)
            end if
        end do
    end do
!   source term.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   convective and viscous terms.
    D   =  0.0d0
    do im=1,mesh(lev)%n_mortar
        sL      =  mesh(lev)%mortar_LR(1,im)
        eL      =  mesh(lev)%mortar_LR(2,im)
        sR      =  mesh(lev)%mortar_LR(3,im)
        eR      =  mesh(lev)%mortar_LR(4,im)
        n(1:5)  =  mesh(lev)%mortar_n_vg(1:5,im)
        vg(1:3) =  mesh(lev)%mortar_n_vg(6:8,im)
        if(im .le. mesh(lev)%n_mortar_b) then
            uLf(1:5)=  fv(sR)%u(1:5,eR)
            uRf(1:5)=  fv(sR)%u(1:5,eR)
            tL      =  fv(sR)%turb(1:2,eR)
            tR      =  tL
            nul     =  fv(sR)%mu(1,eR)/uRf(1)
            gra(1:3)=  fv(sR)%turb_gra(1:3,eR)
            gra(4:6)=  0.0d0
        else
            uLf(1:5)=  fv(sL)%u(1:5,eL)
            uRf(1:5)=  fv(sR)%u(1:5,eR)
            tL      =  fv(sL)%turb(1:2,eL)
            tR      =  fv(sR)%turb(1:2,eR)
            gra(1:3)=  0.5d0*(fv(sL)%duc(1:3,eL)+fv(sR)%duc(1:3,eR))
            nul     =  0.5d0*(fv(sL)%mu(1,eL)/uLf(1)+fv(sR)%mu(1,eR)/uRf(1))
            gra(1:3)    = (tR(1)-tL(1))*n(1:3)*n(5)
            gra(4:6)    = (tR(2)-tL(2))*n(1:3)*n(5)
        end if
        nut =  0.5d0*(tL(1)+tR(1))
        unL = (uLf(2)-vg(1))*n(1)+(uLf(3)-vg(2))*n(2)+(uLf(4)-vg(3))*n(3)
        unR = (uRf(2)-vg(1))*n(1)+(uRf(3)-vg(2))*n(2)+(uRf(4)-vg(3))*n(3)

        fLtL=  0.0d0
        fLtR=  0.0d0
        fRtL=  0.0d0
        fRtR=  0.0d0
        if(unL .ge. 0.0d0) then
            fL(1)       =  unL*tL(1)*n(4)
            fLtL(1,1)   =  unL*n(4)
        else
            fL(1)       =  unL*tR(1)*n(4)
            fLtR(1,1)   =  unL*n(4)
        end if
        if(unR .ge. 0.0d0) then
            fR(1)       =  unR*tL(1)*n(4)
            fRtL(1,1)   =  unR*n(4)
        else
            fR(1)       =  unR*tR(1)*n(4)
            fRtR(1,1)   =  unR*n(4)
        end if

        spL         =  0.5d0*uLf(1)*(unL+abs(unL))*n(4)
        spR         =  0.5d0*uRf(1)*(unR-abs(unR))*n(4)
        fL(2)       =  spL*tL(2)+spR*tR(2)
        fR(2)       =  fL(2)
        fLtL(2,2)   =  spL
        fLtR(2,2)   =  spR
        fRtL(2,2)   =  spL
        fRtR(2,2)   =  spR

        rV(1)   = (n(1)*gra(1)+n(2)*gra(2)+n(3)*gra(3))*(nul+nut)*Pr_SA*n(4)
        rV(2)   = (n(1)*gra(4)+n(2)*gra(5)+n(3)*gra(6))*(mul+mut)      *n(4)
        fL      =  fL-rV
        fR      =  fR-rV
        D(1,1)  =  n(4)*n(5)*(nul+nut)*Pr_SA
        D(2,2)  =  n(4)*n(5)*(mul+mut)
        fLtL    =  fLtL+D
        fLtR    =  fLtR-D
        fRtL    =  fRtL+D
        fRtR    =  fRtR-D

        fv(sL)%rhs(1:2,eL)  =  fv(sL)%rhs(1:2,eL)+fL
        if(sec(sR)%is_int)  fv(sR)%rhs(1:2,eR)  =  fv(sR)%rhs(1:2,eR)-fR
        if(is_LHS) then
            eL  =  eL+iA(sL)-1
            eR  =  eR+iA(sR)-1
            call prec_add_eleR(fv_rans_prec, eL, eL, fLtL)
            if(sec(sR)%is_int) then
                call prec_add_eleR(fv_rans_prec, eL, eR,  fLtR)
                call prec_add_eleR(fv_rans_prec, eR, eR, -fRtR)
                call prec_add_eleR(fv_rans_prec, eR, eL, -fRtL)
            end if
        end if
    end do
!   convective and viscous terms.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   source term due to BDF.
    if(is_BDF_now) then
        g   =  1.0d0/dt_uns
        D   =  0.0d0

        do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
            if(sec(isec)%is_ghost .or. sec(isec)%is_bnd)    cycle

            do iele=1,sec(isec)%n_ele
                vol =  sec(isec)%vol( iele)
                n(1)=  fv(isec)%turb(1,iele)
                n(2)=  fv(isec)%turb(2,iele)*fv(isec)%u(1,iele)

                if((uns_iter .le. 1) .and. (.not. is_uns_initialized)) then
                    u(1:2)  =  n(1:2)-fv(isec)%uns_turb(1:2,iele)
                else
                    u(1:2)  =  1.5d0*n(1:2)-2.0d0*fv(isec)%uns_turb(1:2,iele) &
                            & +0.5d0*fv(isec)%uns_turb(3:4,iele)
                end if
                fv(isec)%rhs(1:2,iele)  =  fv(isec)%rhs(1:2,iele)+vol*g*u(1:2)

                if(is_LHS) then
                    eL  =  iele+iA(isec)-1
                    D(1,1)  =  1.5d0*vol*g
                    D(2,2)  =  D(1,1)
                    call prec_add_eleR(fv_rans_prec, eL, eL, D)
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
    end subroutine fv_SAM_get_rhs
!-------------------------------------------------------------------------------
!   get localized shape factor with the Coder approach.
!-------------------------------------------------------------------------------
    subroutine fv_get_shape_factor
    use var_kind_def
    use var_fv
    use var_mesh
    use var_slv, only: is_vis_cal
    implicit none
    integer(dpI):: isec,iele,sL,eL,sR,eR,im
    real   (dpR):: d(3),v(3),n(3),u(3),H

    if(.not. is_vis_cal)    return

    do isec=mesh(0)%sec_1,mesh(0)%sec_0
        if(sec(isec)%is_int)    fv(isec)%rhs=  0.0d0
    end do

    do im=1,mesh(0)%n_mortar
        sL  =  mesh(0)%mortar_LR(1,im)
        eL  =  mesh(0)%mortar_LR(2,im)
        sR  =  mesh(0)%mortar_LR(3,im)
        eR  =  mesh(0)%mortar_LR(4,im)
        n(1:3)  =  mesh(0)%mortar_n_vg(1:3,im)*mesh(0)%mortar_n_vg(4,im)
        v(1:3)  =  mesh(0)%mortar_n_vg(6:8,im)

        if(im .le. mesh(0)%n_mortar_b) then
            u(1:3)  =  fv(sR)%u(1,eR)*(fv(sR)%u(2:4,eR)-v(1:3))
            d(1:3)  =  sec(sL)%wn(1:3,eL)
        else
            u(1:3)  =  0.25d0*(fv(sL)%u(1,eL)+fv(sR)%u(1,eR)) &
                    & *(fv(sL)%u(2:4,eL)+fv(sR)%u(2:4,eR)-2.0d0*v(1:3))
            d(1:3)  =  0.5d0*(sec(sL)%wn(1:3,eL)+sec(sR)%wn(1:3,eR))
        end if
        v(1:3)  = (d(1)*u(1)+d(2)*u(2)+d(3)*u(3))*n(1:3)
        fv(sL)%rhs(1:3,eL)  =  fv(sL)%rhs(1:3,eL)+v(1:3)
        if(sec(sR)%is_int)  fv(sR)%rhs(1:3,eR)  =  fv(sR)%rhs(1:3,eR)-v(1:3)
    end do

    do isec=mesh(0)%sec_1,mesh(0)%sec_0
        if(.not. sec(isec)%is_int)  cycle
        do iele=1,sec(isec)%n_ele
            H   = (fv(isec)%rhs(1,iele)*sec(isec)%wn(1,iele) &
                & +fv(isec)%rhs(2,iele)*sec(isec)%wn(2,iele) &
                & +fv(isec)%rhs(3,iele)*sec(isec)%wn(3,iele))&
                & *sec(isec)%dnw(iele)**2/fv(isec)%mu(1,iele)
            H   =  max(-0.24d0, min(H, 2.0d2))
            fv(isec)%rhs(1,iele)=  0.1818d0+sqrt((H+2.3609d0)/0.5178d0)
        end do
    end do

    return
    end subroutine fv_get_shape_factor
!-------------------------------------------------------------------------------
!   cal Rhs, SA with the one-equation Coder transition model.
!-------------------------------------------------------------------------------
    subroutine fv_SAC_get_rhs(lev,is_LHS)
    use var_kind_def
    use var_bndv, only: Tu_fs
    use var_fv
    use var_mesh
    use var_prec
    use var_turb
    use var_uns_cal, only: is_BDF_now,dt_uns,uns_iter,is_uns_initialized
    implicit none
    logical(dpL),intent(in):: is_LHS
    integer(dpI),intent(in):: lev
    integer(dpI):: isec,iele,im,sL,eL,sR,eR,iA(1000)
    real   (dpR):: n(5),uLf(5),uRf(5),tL(2),tR(2),nul,nut,gra(6),fL(2),fR(2), &
                &  d(2,2),g,unL,unR,vol,dnw,rh,mul,mut,ux,uy,uz,vx,vy,vz,wx,wy, &
                &  wz,gran,oo,d2_1,chi,fv1,fv2,s,st,r,fw,pp,dd,chi_nut,fv1_nut, &
                &  fv2_nut,r_nut,g_nut,fw_nut,fLtL(2,2),fLtR(2,2),fRtL(2,2), &
                &  fRtR(2,2),vg(3),u(3),spL,spR,rV(2),s12,s13,s23,amp,H12,kv,Red2, &
                &  Rev0,Rev,F_crit,l,m,F_growth,n_Re,n_crit,f_t2,DH,rtmp

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

    call fv_get_shape_factor

!   ----------------------------------------------------------------------------
!   source term.
    D   =  0.0d0
    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(sec(isec)%is_ghost .or. sec(isec)%is_bnd)    cycle

        do iele=1,sec(isec)%n_ele
            vol     =  sec(isec)%vol( iele)
            dnw     =  sec(isec)%dnw( iele)
            rh      =  fv (isec)%u (1,iele)
            u(1:3)  =  fv (isec)%u (2:4,iele)
            mul     =  fv (isec)%mu(1,iele)
            mut     =  fv (isec)%mu(2,iele)
            nul     =  mul/rh
            nut     =  fv (isec)%turb(1,iele)
            amp     =  fv (isec)%turb(2,iele)

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

            s12 =  0.5d0*(uy+vx)
            s13 =  0.5d0*(uz+wx)
            s23 =  0.5d0*(vz+wy)
            S   =  sqrt(2.0d0*(ux*ux+vy*vy+wz*wz)+4.0d0*(s12*s12+s13*s13+s23*s23))
            d2_1=  1.0d0/(dnw*dnw)

            H12 =  fv(isec)%rhs(1,iele)
            DH  =  H12/(0.6202d0*H12-0.6387d0)
            kv  =  0.2231d0*H12**2-0.1617d0*H12+0.1121d0
            Red2=  0.7d0*tanh(1.4d1/(H12-1.0d0)-9.24d0)+2.492d0/(H12-1.0d0)**0.43d0+0.62d0
            Rev0=  kv*1.0d1**Red2
            Rev =  rh*S*dnw*dnw/(mul+mut)
            if(Rev .lt. Rev0) then
                F_crit  =  0.0d0
            else
                F_crit  =  1.0d0
            end if
            l       = (6.54d0*H12-14.07d0)/(H12*H12)
            m       = (0.058d0*(H12-4.0d0)/(H12-1.0d0)-0.068d0)/l
            F_growth=  0.5d0*DH*l*(1.0d0+m)
            n_Re    =  0.028d0*(H12-1.0d0)-0.0345d0*exp(-(3.87d0/(H12-1.0d0)-2.52d0)**2)
            N_crit  =  max(9.0d0, min(-8.43d0-2.4d0*log(Tu_fs), 9.0d0))

            chi =  nut/nul
            f_t2=  1.2d0*(1.0d0-exp(2.0d0*(amp-N_crit)))*exp(-5.0d-2*chi**2)
            rtmp=  chi**3
            fv1 =  rtmp/(rtmp+3.57911d2)
            fv2 =  1.0d0-chi/(1.0d0+chi*fv1)
            st  =  max(oo+nut*k2_1*d2_1*fv2, 0.3d0*oo)
            r   =  min(nut*k2_1*d2_1/st, 1.0d1)
            g   =  r+0.3d0*(r**6-r)
            fw  =  g*(6.5d1/(g**6+6.4d1))**(1.0d0/6.0d0)

            pp  =  cb1*st*nut*(1.0d0-f_t2)+0.933d0*gran
            dd  = (cw1*fw-cb1*k2_1*f_t2)*nut**2*d2_1

            fv(isec)%rhs(1,iele)= -vol*(pp-dd)
            fv(isec)%rhs(2,iele)= -vol*rh*oo*F_crit*F_growth*n_Re

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
                D(1,1)  =  vol*g
                D(2,2)  =  0.0d0
                call prec_add_eleR(fv_rans_prec, eL, eL, D)
            end if
        end do
    end do
!   source term.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   convective and viscous terms.
    D   =  0.0d0
    do im=1,mesh(lev)%n_mortar
        sL      =  mesh(lev)%mortar_LR(1,im)
        eL      =  mesh(lev)%mortar_LR(2,im)
        sR      =  mesh(lev)%mortar_LR(3,im)
        eR      =  mesh(lev)%mortar_LR(4,im)
        n(1:5)  =  mesh(lev)%mortar_n_vg(1:5,im)
        vg(1:3) =  mesh(lev)%mortar_n_vg(6:8,im)
        if(im .le. mesh(lev)%n_mortar_b) then
            uLf(1:5)=  fv(sR)%u(1:5,eR)
            uRf(1:5)=  fv(sR)%u(1:5,eR)
            tL      =  fv(sR)%turb(1:2,eR)
            tR      =  tL
            nul     =  fv(sR)%mu(1,eR)/uRf(1)
            gra(1:3)=  fv(sR)%turb_gra(1:3,eR)
            gra(4:6)=  0.0d0
        else
            uLf(1:5)=  fv(sL)%u(1:5,eL)
            uRf(1:5)=  fv(sR)%u(1:5,eR)
            tL      =  fv(sL)%turb(1:2,eL)
            tR      =  fv(sR)%turb(1:2,eR)
            gra(1:3)=  0.5d0*(fv(sL)%duc(1:3,eL)+fv(sR)%duc(1:3,eR))
            nul     =  0.5d0*(fv(sL)%mu(1,eL)/uLf(1)+fv(sR)%mu(1,eR)/uRf(1))
            gra(1:3)    = (tR(1)-tL(1))*n(1:3)*n(5)
            gra(4:6)    = (tR(2)-tL(2))*n(1:3)*n(5)
        end if
        nut =  0.5d0*(tL(1)+tR(1))
        unL = (uLf(2)-vg(1))*n(1)+(uLf(3)-vg(2))*n(2)+(uLf(4)-vg(3))*n(3)
        unR = (uRf(2)-vg(1))*n(1)+(uRf(3)-vg(2))*n(2)+(uRf(4)-vg(3))*n(3)

        fLtL=  0.0d0
        fLtR=  0.0d0
        fRtL=  0.0d0
        fRtR=  0.0d0
        if(unL .ge. 0.0d0) then
            fL(1)       =  unL*tL(1)*n(4)
            fLtL(1,1)   =  unL*n(4)
        else
            fL(1)       =  unL*tR(1)*n(4)
            fLtR(1,1)   =  unL*n(4)
        end if
        if(unR .ge. 0.0d0) then
            fR(1)       =  unR*tL(1)*n(4)
            fRtL(1,1)   =  unR*n(4)
        else
            fR(1)       =  unR*tR(1)*n(4)
            fRtR(1,1)   =  unR*n(4)
        end if

        spL         =  0.5d0*uLf(1)*(unL+abs(unL))*n(4)
        spR         =  0.5d0*uRf(1)*(unR-abs(unR))*n(4)
        fL(2)       =  spL*tL(2)+spR*tR(2)
        fR(2)       =  fL(2)
        fLtL(2,2)   =  spL
        fLtR(2,2)   =  spR
        fRtL(2,2)   =  spL
        fRtR(2,2)   =  spR

        rV(1)   = (n(1)*gra(1)+n(2)*gra(2)+n(3)*gra(3))*(nul+nut)*Pr_SA*n(4)
        rV(2)   = (n(1)*gra(4)+n(2)*gra(5)+n(3)*gra(6))*(mul+mut)      *n(4)
        fL      =  fL-rV
        fR      =  fR-rV
        D(1,1)  =  n(4)*n(5)*(nul+nut)*Pr_SA
        D(2,2)  =  n(4)*n(5)*(mul+mut)
        fLtL    =  fLtL+D
        fLtR    =  fLtR-D
        fRtL    =  fRtL+D
        fRtR    =  fRtR-D

        fv(sL)%rhs(1:2,eL)  =  fv(sL)%rhs(1:2,eL)+fL
        if(sec(sR)%is_int)  fv(sR)%rhs(1:2,eR)  =  fv(sR)%rhs(1:2,eR)-fR
        if(is_LHS) then
            eL  =  eL+iA(sL)-1
            eR  =  eR+iA(sR)-1
            call prec_add_eleR(fv_rans_prec, eL, eL, fLtL)
            if(sec(sR)%is_int) then
                call prec_add_eleR(fv_rans_prec, eL, eR,  fLtR)
                call prec_add_eleR(fv_rans_prec, eR, eR, -fRtR)
                call prec_add_eleR(fv_rans_prec, eR, eL, -fRtL)
            end if
        end if
    end do
!   convective and viscous terms.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   source term due to BDF.
    if(is_BDF_now) then
        g   =  1.0d0/dt_uns
        D   =  0.0d0

        do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
            if(sec(isec)%is_ghost .or. sec(isec)%is_bnd)    cycle

            do iele=1,sec(isec)%n_ele
                vol =  sec(isec)%vol( iele)
                n(1)=  fv(isec)%turb(1,iele)
                n(2)=  fv(isec)%turb(2,iele)*fv(isec)%u(1,iele)

                if((uns_iter .le. 1) .and. (.not. is_uns_initialized)) then
                    u(1:2)  =  n(1:2)-fv(isec)%uns_turb(1:2,iele)
                else
                    u(1:2)  =  1.5d0*n(1:2)-2.0d0*fv(isec)%uns_turb(1:2,iele) &
                            & +0.5d0*fv(isec)%uns_turb(3:4,iele)
                end if
                fv(isec)%rhs(1:2,iele)  =  fv(isec)%rhs(1:2,iele)+vol*g*u(1:2)

                if(is_LHS) then
                    eL  =  iele+iA(isec)-1
                    D(1,1)  =  1.5d0*vol*g
                    D(2,2)  =  D(1,1)
                    call prec_add_eleR(fv_rans_prec, eL, eL, D)
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
    end subroutine fv_SAC_get_rhs
!-------------------------------------------------------------------------------
!   exchange intermittency.
!-------------------------------------------------------------------------------
    subroutine fv_exchange_intermittency(lev)
    use var_kind_def
    use var_fv
    use var_mesh
    use var_parallel
    implicit none
    integer(dpI),intent(in):: lev
    integer(dpI):: isr,isec,iele,i,idx,ip_remote

!   ----------------------------------------------------------------------------
!   prepare data.
    do isr=mesh(lev)%p2p_1,mesh(lev)%p2p_0
        idx =  1
        do i=1,p2p(isr)%n_ele_send
            isec=  p2p(isr)%id_ele_send(1,i)
            iele=  p2p(isr)%id_ele_send(2,i)
            p2p(isr)%rsend(idx) =  fv(isec)%intermittency(iele)
            idx =  idx+1
        end do
        p2p(isr)%n_send =  idx-1
        p2p(isr)%n_recv =  p2p(isr)%n_ele_recv
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
    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(.not. sec(isec)%is_ghost)    cycle
        do i=1,sec(isec)%n_ele
            isr =  sec(isec)%id_recv(1,i)
            iele=  sec(isec)%id_recv(2,i)
            fv(isec)%intermittency(i)   =  p2p(isr)%rrecv(iele)
        end do
    end do
!   unpack data received.
!   ----------------------------------------------------------------------------

    return
    end subroutine fv_exchange_intermittency
