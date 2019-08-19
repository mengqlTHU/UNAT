!-------------------------------------------------------------------------------
!   transformation between entropy and conservative variables.
!-------------------------------------------------------------------------------
    subroutine S_over_U(u,SU)
    use var_kind_def
    use var_air, only: cv,gk,gk1,rr
    implicit none
    real(dpR),intent(in):: u(*)
    real(dpR):: SU(*),s

    s       =  cv*log(u(5)/u(1)**gk)
    SU(1)   =  gk/gk1-s/rr-0.5d0*u(1)*(u(2)*u(2)+u(3)*u(3)+u(4)*u(4))/u(5)
    SU(2:4) =  u(1)*u(2:4)/u(5)
    SU(5)   = -u(1)/u(5)

    return
    end subroutine S_over_U
!-------------------------------------------------------------------------------
!   cal Rhs.
!-------------------------------------------------------------------------------
    subroutine fv_get_entropy_residual(lev)
    use var_kind_def
    use var_fv
    use var_mesh
    use var_slv, only: is_vis_cal
    implicit none
    integer(dpI),intent(in):: lev
    integer(dpI):: isec,iele
    real   (dpR):: mu,S(3,3),gra(9)

    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(sec(isec)%is_ghost .or. sec(isec)%is_bnd)    cycle
        fv(isec)%duc=  0.0d0
    end do

    if(is_kep) then
    else
        if(rhs_lev(lev) .eq. 1) then
            call fv_get_rhs_cen_entropy_residual(lev)
        else
        end if
    end if

    if(is_vis_cal)  call fv_get_rhs_vis_entropy_residual(lev)

    call fv_get_shock_sensor(lev)
    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(sec(isec)%is_ghost .or. sec(isec)%is_bnd)    cycle
        do iele=1,sec(isec)%n_ele
            fv(isec)%duc(1,iele)=  fv(isec)%duc(1,iele)/sec(isec)%vol(iele)

            gra(1:9)=  fv(isec)%gra(4:12,iele)
            S(1,1)  =  gra(1)
            S(2,1)  =  0.5d0*(gra(2)+gra(4))
            S(3,1)  =  0.5d0*(gra(3)+gra(7))
            S(1,2)  =  S(2,1)
            S(2,2)  =  gra(5)
            S(3,2)  =  0.5d0*(gra(6)+gra(8))
            S(1,3)  =  S(3,1)
            S(2,3)  =  S(3,2)
            S(3,3)  =  gra(9)
            mu      =  fv(isec)%mu(1,iele)+fv(isec)%mu(2,iele)
            fv(isec)%duc(2,iele)=  4.0d0*(s(1,2)**2+s(1,3)**2+s(2,3)**2) &
                & +2.0d0*((s(1,1)-s(2,2))**2+(s(1,1)-s(3,3))**2+(s(2,2)-s(3,3))**2)/3.0d0
            fv(isec)%duc(2,iele)=  log10(mu*fv(isec)%duc(2,iele))
        end do
    end do
    call fv_wr_entropy_residual

    return
    end subroutine fv_get_entropy_residual
!-------------------------------------------------------------------------------
!   cal Rhs, convective part, fv.
!-------------------------------------------------------------------------------
    subroutine fv_get_rhs_cen_entropy_residual(lev)
    use var_kind_def
    use var_air, only: gk,gk1
    use var_fv
    use var_mesh
    use var_temp, only: rhsl,rhsD
    implicit none
    integer(dpI),intent(in):: lev
    integer(dpI):: im,sL,eL,sR,eR
    real   (dpR):: u(5),un,rhun,a,H,e2,e4,d1(5),d3(5),spra,n(5), &
                &  ke,d(5),ssL,ssR,akhL(3),akhR(3),vg(3),vn,SUL(5),SUR(5)

    call fv_get_ss(lev)

    call fv_get_rhs_boundary_entropy_residual(lev)
    do im=1+mesh(lev)%n_mortar_b,mesh(lev)%n_mortar
        sL  =  mesh(lev)%mortar_LR(1,im)
        eL  =  mesh(lev)%mortar_LR(2,im)
        sR  =  mesh(lev)%mortar_LR(3,im)
        eR  =  mesh(lev)%mortar_LR(4,im)
        n (1:5) =  mesh(lev)%mortar_n_vg(1:5,im)
        vg(1:3) =  mesh(lev)%mortar_n_vg(6:8,im)

        u(1:5)  =  0.5d0*(fv(sL)%u(1:5,eL)+fv(sR)%u(1:5,eR))
        call u_to_akh(1, fv(sL)%u(1,eL), akhL)
        call u_to_akh(1, fv(sR)%u(1,eR), akhR)
        call S_over_U(fv(sL)%u(1,eL), SUL)
        call S_over_U(fv(sR)%u(1,eR), SUR)
        a           =  sqrt(gk*u(5)/u(1))
        ke          =  0.5d0*(u(2)*u(2)+u(3)*u(3)+u(4)*u(4))
        H           =  a*a/gk1+ke
        un          =  u (2)*n(1)+u (3)*n(2)+u (4)*n(3)
        vn          =  vg(1)*n(1)+vg(2)*n(2)+vg(3)*n(3)
        rhun        =  u(1)*(un-vn)
        rhsl(1  ,1) =  rhun
        rhsl(2:4,1) =  rhun*u(2:4)+n(1:3)*u(5)
        rhsl(5  ,1) =  rhun*0.5d0*(akhL(3)+akhR(3))+u(5)*vn
        spra    =  abs(un)+a
        d1(1:5) =  fv(sR)%uc    (1:5,eR)-fv(sL)%uc    (1:5,eL)
        d3(1:5) =  fv(sR)%ss_lap(1:5,eR)-fv(sL)%ss_lap(1:5,eL)
        ssL     =  fv(sL)%ss_lap(6,eL)
        ssR     =  fv(sR)%ss_lap(6,eR)

        if(.true.) then
            e2  =  1.0d0*max(ssL, ssR)
            e4  =  max(0.0d0, 1.0d0/5.0d1-e2)
            d   =  e2*d1-e4*d3
            rhsD(1:5,1) =  spra*d
        else
        end if
        rhsl(1:5,1) =  n(4)*(rhsl(1:5,1)-rhsD(1:5,1))

        fv(sL)%duc(1,eL)=  fv(sL)%duc(1,eL)+rhsl(1,1)*SUL(1)+rhsl(2,1)*SUL(2) &
            & +rhsl(3,1)*SUL(3)+rhsl(4,1)*SUL(4)+rhsl(5,1)*SUL(5)
        if(sec(sR)%is_int) then
            fv(sR)%duc(1,eR)=  fv(sR)%duc(1,eR)-rhsl(1,1)*SUR(1)-rhsl(2,1)*SUR(2) &
                & -rhsl(3,1)*SUR(3)-rhsl(4,1)*SUR(4)-rhsl(5,1)*SUR(5)
        end if
    end do

    return
    end subroutine fv_get_rhs_cen_entropy_residual
!-------------------------------------------------------------------------------
!   cal Rhs, convective part, boundary.
!-------------------------------------------------------------------------------
    subroutine fv_get_rhs_boundary_entropy_residual(lev)
    use var_kind_def
    use var_cgns
    use var_fv
    use var_mesh
    use var_temp, only: uLf,uRf,vgf,fac_1d,rhsl
    implicit none
    integer(dpI),intent(in):: lev
    logical(dpL):: is_upwind
    integer(dpI):: im,sL,eL,sR,eR,bct
    real   (dpR):: SUL(5)

    do im=1,mesh(lev)%n_mortar_b
        sL  =  mesh(lev)%mortar_LR(1,im)
        eL  =  mesh(lev)%mortar_LR(2,im)
        sR  =  mesh(lev)%mortar_LR(3,im)
        eR  =  mesh(lev)%mortar_LR(4,im)
        fac_1d(1:4,1)   =  mesh(lev)%mortar_n_vg(1:4,im)
        vgf   (1:3,1)   =  mesh(lev)%mortar_n_vg(6:8,im)
        bct     =  sec(sR)%bct

        is_upwind   = (bct .eq. BCSymmetryPlane)
        if(is_upwind) then
            uLf(1:5,1)  =  fv(sR)%uL(1:5,eR)
            uRf(1:5,1)  =  fv(sR)%uR(1:5,eR)
            call rhs_conv_roe(.true., .true., 1)
        else
            call u_to_F(1, fac_1d(1:3,1), fv(sR)%u(1:5,eR), vgf, rhsl)
            rhsl(1:5,1) =  rhsl(1:5,1)*fac_1d(4,1)
        end if

        call S_over_U(fv(sL)%u(1,eL), SUL)
        fv(sL)%duc(1,eL)=  fv(sL)%duc(1,eL)+rhsl(1,1)*SUL(1)+rhsl(2,1)*SUL(2) &
            & +rhsl(3,1)*SUL(3)+rhsl(4,1)*SUL(4)+rhsl(5,1)*SUL(5)
    end do

    return
    end subroutine fv_get_rhs_boundary_entropy_residual
!-------------------------------------------------------------------------------
!   cal Rhs, viscous part, fv.
!-------------------------------------------------------------------------------
    subroutine fv_get_rhs_vis_entropy_residual(lev)
    use var_kind_def
    use var_cgns, only: BCWallViscous
    use var_fv
    use var_global, only: n_dim
    use var_mesh
    use var_temp, only: uLf,uRf,fac_1d,rhsD
    use var_turb, only: is_tur_cal,is_KO
    implicit none
    integer(dpI),intent(in):: lev
    integer(dpI):: im,sL,eL,sR,eR,j,k
    real   (dpR):: g(12),u(5),d(3),du(4),tL,tR,nne(3),dr,mu(2),SUL(5),SUR(5),tke(1),rtmp

    d   =  0.0d0
    mu  =  0.0d0
    do im=1,mesh(lev)%n_mortar
        sL  =  mesh(lev)%mortar_LR(1,im)
        eL  =  mesh(lev)%mortar_LR(2,im)
        sR  =  mesh(lev)%mortar_LR(3,im)
        eR  =  mesh(lev)%mortar_LR(4,im)
        fac_1d(1:4,1)   =  mesh(lev)%mortar_n_vg(1:4,im)

        if(im .le. mesh(lev)%n_mortar_b) then
            u  (1 :5 )  =  fv(sR)%u(1:5,eR)
            mu (1    )  =  fv(sR)%mu(1,eR)
            if(is_tur_cal)  mu(2)   =  fv(sR)%mu  (2,eR)
            if(is_KO     )  tke     =  fv(sR)%turb(1,eR)
            g  (1 :9 )  =  fv(sR)%gra(4 :12,eR)
            g  (10:12)  =  fv(sR)%gra(16:18,eR)
            if(sec(sR)%bct .eq. BCWallViscous)  g(10:12)=  0.0d0
        else
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

            d(1:n_dim)  =  sec(sR)%cen(1:n_dim,eR)-sec(sL)%cen(1:n_dim,eL)
            rtmp        =  1.0d0/sqrt(d(1)**2+d(2)**2+d(3)**2)
            d           =  d*rtmp
            du(1:3)     =  rtmp*(uRf(2:4,1)-uLf(2:4,1))
            du(4  )     =  rtmp*(tR-tL)

            rtmp=  fac_1d(1,1)*d(1)+fac_1d(2,1)*d(2)+fac_1d(3,1)*d(3)
            nne(1:3)=  fac_1d(1:3,1)/rtmp
            do j=1,4
                k       =  3*j-2
                dr      =  du(j)-g(k)*d(1)-g(k+1)*d(2)-g(k+2)*d(3)
                g(k:k+2)=  g(k:k+2)+dr*nne(1:3)
            end do
        end if

        call get_vis_fac(1, fac_1d, u, mu, g, tke, rhsD)
        call S_over_U(fv(sL)%u(1,eL), SUL)
        call S_over_U(fv(sR)%u(1,eR), SUR)

        fv(sL)%duc(1,eL)=  fv(sL)%duc(1,eL)-rhsD(1,1)*SUL(1)-rhsD(2,1)*SUL(2) &
            & -rhsD(3,1)*SUL(3)-rhsD(4,1)*SUL(4)-rhsD(5,1)*SUL(5)
        if(sec(sR)%is_int) then
            fv(sR)%duc(1,eR)=  fv(sR)%duc(1,eR)+rhsD(1,1)*SUR(1)+rhsD(2,1)*SUR(2) &
                & +rhsD(3,1)*SUR(3)+rhsD(4,1)*SUR(4)+rhsD(5,1)*SUR(5)
        end if
    end do

    return
    end subroutine fv_get_rhs_vis_entropy_residual
!-------------------------------------------------------------------------------
!   output cell-centered solution, tecplot.
!-------------------------------------------------------------------------------
    subroutine fv_wr_entropy_residual
    use var_kind_def
    use var_cgns
    use var_fv
    use var_global, only: mesh_name,n_dim
    use var_mesh
    use var_parallel, only: myid
    implicit none
    character(len=1000):: str,tec
    logical(dpL):: ltmp
    integer(dpI):: n_vtx,isec,i,iele,v(8),L,R
    integer(dpI),allocatable:: ibuf(:)

    write(str,*),myid
    str =  trim(adjustl(mesh_name))//'_entropy_'//trim(adjustl(str))//'.dat'
    open(unit=10,file=trim(adjustl(str)))

    tec =  'variables="CoordinateX","CoordinateY"'
    if(n_dim .eq. 3)    tec =  trim(adjustl(tec))//',"CoordinateZ"'
    tec =  trim(adjustl(tec))//',"entropy","dissipation","ss"'
    write(unit=10,fmt='(A)'),trim(tec)

    allocate(ibuf(mesh(0)%n_vtx))
    do isec=mesh(0)%sec_1,mesh(0)%sec_0
        if(sec(isec)%is_ghost .or. sec(isec)%is_bnd)    cycle

        ibuf=  0
        do iele=1,sec(isec)%n_ele
        do i=1,sec(isec)%npe
            ibuf(sec(isec)%n2e(i,iele)) =  1
        end do
        end do

        n_vtx   =  0
        do i=1,mesh(0)%n_vtx
            if(ibuf(i) .le. 0)  cycle
            n_vtx       =  n_vtx+1
            ibuf(n_vtx) =  i
        end do

        write(str,*),n_vtx
        tec =  'zone N='//trim(adjustl(str))//',E='
        write(str,*),sec(isec)%n_ele
        tec =  trim(tec)//trim(adjustl(str))

        if(sec(isec)%is_quad) then
            str =  'FEQUADRILATERAL'
        else
            str =  'FETRIANGLE'
        end if
        tec =  trim(tec)//',zonetype='//trim(adjustl(str))//',datapacking=block'
        tec =  trim(tec)//',varlocation=([3-5]=cellcentered)'
        write(unit=10,fmt='(A)'),trim(adjustl(tec))

        do i=1,n_vtx
            write(unit=10,fmt=*),mesh(0)%xyz(1,ibuf(i))
        end do
        do i=1,n_vtx
            write(unit=10,fmt=*),mesh(0)%xyz(2,ibuf(i))
        end do
        write(unit=10,fmt=*),fv(isec)%duc(1,:)
        write(unit=10,fmt=*),fv(isec)%duc(2,:)
        write(unit=10,fmt=*),fv(isec)%shock_sensor
        do iele=1,sec(isec)%n_ele
            do i=1,sec(isec)%npe
                call ib_search(1, n_vtx, 1, 1, ibuf, sec(isec)%n2e(i,iele), ltmp, L, R)
                if((.not. ltmp) .or. (L .ne. R))    stop 'Error: fv_wr_tec_cc.'
                v(i)=  L
            end do
            write(unit=10,fmt='(8I9)'),v(1:sec(isec)%npe)
        end do
    end do

    return
    end subroutine fv_wr_entropy_residual
!-------------------------------------------------------------------------------
!   get Q criterion.
!-------------------------------------------------------------------------------
    subroutine fv_get_q
    use var_kind_def
    use var_fv
    use var_global, only: R23,uref
    use var_mesh
    implicit none
    integer(dpI):: isec,iele
    real   (dpR):: ux,uy,uz,vx,vy,vz,wx,wy,wz,oo,s12,s13,s23,ss,Re_v,ptm,ptm1,R_t, &
                &  F3,rh,dnw,mul,mut,F_turb,omega,k,TU_L,u_norm,u_x,u_y,u_z, &
                &  u_s,Lambda_L,F_PG,F_onset1,F_onset2,F_onset3,Re_thetac,u(3),nul, &
                &  Pohl,U_edge

    do isec=mesh(0)%sec_1,mesh(0)%sec_0
        if(.not. sec(isec)%is_int)  cycle
        if(.not. allocated(fv(isec)%q)) allocate(fv(isec)%q(sec(isec)%n_ele))

        do iele=1,sec(isec)%n_ele
            ux  =  fv(isec)%gra(4 ,iele)
            uy  =  fv(isec)%gra(5 ,iele)
            uz  =  fv(isec)%gra(6 ,iele)
            vx  =  fv(isec)%gra(7 ,iele)
            vy  =  fv(isec)%gra(8 ,iele)
            vz  =  fv(isec)%gra(9 ,iele)
            wx  =  fv(isec)%gra(10,iele)
            wy  =  fv(isec)%gra(11,iele)
            wz  =  fv(isec)%gra(12,iele)
            oo  =  sqrt((uy-vx)**2+(uz-wx)**2+(vz-wy)**2)
            s12 =  0.5d0*(uy+vx)
            s13 =  0.5d0*(uz+wx)
            s23 =  0.5d0*(vz+wy)
            ss  =  sqrt(2.0d0*(ux*ux+vy*vy+wz*wz)+4.0d0*(s12*s12+s13*s13+s23*s23))
            fv(isec)%q(iele)=  0.5d0*(oo*oo-ss*ss)

            dnw     =  sec(isec)%dnw(iele)
            rh      =  fv(isec)%u(1,iele)
            u(1:3)  =  fv(isec)%u(2:4,iele)
            mul     =  fv(isec)%mu(1,iele)
            mut     =  fv(isec)%mu(2,iele)
            nul     =  mul/rh
            Re_v    =  rh*dnw*dnw*SS/mul
            R_T     =  mut/mul
            F_turb  =  exp(-(0.5d0*R_T)**4)
            omega   =  ss/0.3d0
            k       =  mut*omega/rh
            U_edge  =  max(omega*dnw, 2.0d-1*uref)
            Tu_L    =  min(sqrt(R23*k)/U_edge, 1.0d0)*1.0d2
            u_norm  =  sqrt(u(1)**2+u(2)**2+u(3)**2)
            U_x     = (u(1)*ux+u(2)*vx+u(3)*wx)/sqrt(u_norm)
            U_y     = (u(1)*uy+u(2)*vy+u(3)*wy)/sqrt(u_norm)
            U_z     = (u(1)*uz+u(2)*vz+u(3)*wz)/sqrt(u_norm)
            U_s     = (u(1)*U_x+u(2)*U_y+u(3)*U_z)/u_norm
            Pohl    =  U_s*dnw*dnw/nul
!           Pohl    = -vy *dnw*dnw/nul
            Pohl    =  min(0.1d0, max(Pohl, -0.1d0))
            lambda_L=  min(max(-1.0d0, 7.57d-3*Pohl+1.28d-2), 1.0d0)
            if(lambda_L .ge. 0.0d0) then
                F_PG=  min(1.0d0+14.68d0*lambda_L, 1.5d0)
            else
                F_PG=  min(1.0d0-7.34d0*lambda_L, 3.0d0)
            end if
            F_PG        =  max(F_PG, 0.0d0)
            Re_thetac   =  1.0d2+1.0d3*exp(-Tu_L*F_PG)
!           rtmp        =  2.2d2*exp(-Tu_fs/9.0d-3)*(atan(Pohl/2.0d-2+0.842d0)-0.7d0)
!           Re_thetac   =  1.045d3*exp(-Tu_fs/1.0d-2)+rtmp+1.55d2
            F_onset1    =  Re_v/(2.2d0*Re_thetac)
            F_onset2    =  min(F_onset1, 2.0d0)
            F_onset3    =  max(1.0d0-(R_T/3.5d0)**3, 0.0d0)

            if(Re_v .le. 1.0d3) then
                ptm1=  1.0d0-(3.28d-4*Re_v-3.94d-7*Re_v**2+1.43d-10*Re_v**3)
            else
                ptm1=  1.0d0-(0.12d0+1.0d-5*Re_v)
            end if
            F3  =  exp(-(R_t/5.0d0)**4)
            ptm =  max(0.0d0, min(1.0d0, 1.0d0-0.94d0*ptm1*F3))

            fv(isec)%rhs(1,iele)=  F_onset2
            fv(isec)%rhs(2,iele)=  F_onset3
            fv(isec)%rhs(3,iele)=  max(F_onset2-F_onset3, 0.0d0)
            fv(isec)%rhs(4,iele)=  fv(isec)%intermittency(iele)
            fv(isec)%rhs(5,iele)=  ptm
        end do
    end do

    return
    end subroutine fv_get_q
