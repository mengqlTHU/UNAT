!-------------------------------------------------------------------------------
!   cal gradient.
!-------------------------------------------------------------------------------
    subroutine fv_SST_get_gra(lev)
    use var_kind_def
    use var_fv
    use var_mesh
    use var_parallel
    use var_slv
    implicit none
    integer(dpI),intent(in):: lev
    integer(dpI):: isec,iele,i,ip_remote,isr,LDA,idx,im,sL,eL,sR,eR
    real   (dpR):: g(6),n(3),d,dt,gn

!   ----------------------------------------------------------------------------
!   cal gradient for the elements involved in the data exchange.
    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(sec(isec)%is_int)    call get_gra(isec, 1)
    end do
!   cal gradient for the elements involved in the data exchange.
!   ----------------------------------------------------------------------------

    LDA =  6

!   ----------------------------------------------------------------------------
!   prepare data.
    do isr=mesh(lev)%p2p_1,mesh(lev)%p2p_0
        idx =  1
        do i=1,p2p(isr)%n_ele_send_face
            isec=  p2p(isr)%id_ele_send(1,i)
            iele=  p2p(isr)%id_ele_send(2,i)
            call DCOPY(LDA, fv(isec)%turb_gra(1,iele), 1, p2p(isr)%rsend(idx), 1)
            idx =  idx+LDA
        end do
        p2p(isr)%n_send =  idx-1
        p2p(isr)%n_recv =  LDA*p2p(isr)%n_ele_recv_face
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

!   ----------------------------------------------------------------------------
!   cal gradient for the elements not involved in the data exchange.
    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(sec(isec)%is_int)    call get_gra(isec, 2)
    end do
!   cal gradient for the elements not involved in the data exchange.
!   ----------------------------------------------------------------------------

    call mpi_waitall(mpi_nreq, mpi_req, mpi_sta, mpi_err)

!   ----------------------------------------------------------------------------
!   unpack data received.
    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(.not. sec(isec)%is_ghost)    cycle
        do i=1,sec(isec)%n_ele
            isr =  sec(isec)%id_recv(1,i)
            iele=  sec(isec)%id_recv(2,i)
            if(iele .gt. p2p(isr)%n_ele_recv_face)  cycle
            call per_rot_vec(12, sec(isec)%per_path(1,i), &
                &  p2p(isr)%rrecv(1+LDA*(iele-1)), fv(isec)%turb_gra(1,i))
        end do
    end do
!   unpack data received.
!   ----------------------------------------------------------------------------

    do im=1,mesh(lev)%n_mortar_b
        sL  =  mesh(lev)%mortar_LR(1,im)
        eL  =  mesh(lev)%mortar_LR(2,im)
        sR  =  mesh(lev)%mortar_LR(3,im)
        eR  =  mesh(lev)%mortar_LR(4,im)
        n(1:3)  =  sec(sL)%cen(1:3,eL)-mesh(lev)%mortar_cen(1:3,im)
        call norm_vec(3, n, d)

        g(1:6)  =  fv(sL)%turb_gra(1:6,eL)

        dt      =  fv(sL)%turb(1,eL)-fv(sR)%turb(1,eR)
        gn      =  g(1)*n(1)+g(2)*n(2)+g(3)*n(3)
        fv(sR)%turb_gra(1:3,eR) =  g(1:3)+(dt/d-gn)*n(1:3)

        dt      =  fv(sL)%turb(2,eL)-fv(sR)%turb(2,eR)
        gn      =  g(4)*n(1)+g(5)*n(2)+g(6)*n(3)
        fv(sR)%turb_gra(4:6,eR) =  g(4:6)+(dt/d-gn)*n(1:3)
    end do

    return
    contains
!   ----------------------------------------------------------------------------
!   get gra for turbulence equation.
!   ----------------------------------------------------------------------------
    subroutine get_gra(isec,mode)
    use var_kind_def
    implicit none
    integer(dpI),intent(in):: isec,mode

    if(.not. sec(isec)%is_int)  return
    if(gradient_method .eq. gradient_GG) then
        call fv_sst_get_gra_gg (isec, mode)
    elseif(gradient_method .eq. gradient_LS) then
        call fv_sst_get_gra_LS (isec, mode)
    else
        call fv_sst_get_gra_LSm(isec, mode)
    end if

    return
    end subroutine get_gra
!   ----------------------------------------------------------------------------
!   get gra for turbulence equation, Green-Gauss.
!   ----------------------------------------------------------------------------
    subroutine fv_sst_get_gra_gg(isec,mode)
    use var_kind_def
    implicit none
    integer(dpI),intent(in):: isec,mode
    integer(dpI):: iele,e1,e0,j,im,s,e
    real   (dpR):: tL(2),vL,n(3),tm(2),vR

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
        tL  =  fv(isec)%turb(1:2,iele)
        vL  =  1.0d0/sec(isec)%vol(iele)
        fv(isec)%turb_gra(1:6,iele) =  0.0d0
        do j=sec(isec)%iA_face_neighbour(iele),sec(isec)%iA_face_neighbour(iele+1)-1
            im  =  sec(isec)%jA_face_neighbour(j)
            if(im .gt. 0) then
                s       =  mesh(lev)%mortar_LR(3, im)
                e       =  mesh(lev)%mortar_LR(4, im)
                n(1:3)  =  mesh(lev)%mortar_n_vg(1:3, im)*mesh(lev)%mortar_n_vg(4, im)
            else
                s       =  mesh(lev)%mortar_LR(1,-im)
                e       =  mesh(lev)%mortar_LR(2,-im)
                n(1:3)  = -mesh(lev)%mortar_n_vg(1:3,-im)*mesh(lev)%mortar_n_vg(4,-im)
            end if
            if(sec(s)%is_bnd) then
                tm  =  fv(s)%turb(1:2,e)
            else
                vR  =  1.0d0/sec(s)%vol(e)
                tm  = (vL*tL+vR*fv(s)%turb(1:2,e))/(vL+vR)
            end if
            fv(isec)%turb_gra(1:3,iele) =  fv(isec)%turb_gra(1:3,iele)+tm(1)*n(1:3)
            fv(isec)%turb_gra(4:6,iele) =  fv(isec)%turb_gra(4:6,iele)+tm(2)*n(1:3)
        end do
        fv(isec)%turb_gra(1:6,iele) =  fv(isec)%turb_gra(1:6,iele)*vL
    end do

    return
    end subroutine fv_sst_get_gra_gg
!   ----------------------------------------------------------------------------
!   get gra for turbulence equation, Least-Square.
!   ----------------------------------------------------------------------------
    subroutine fv_sst_get_gra_LS(isec,mode)
    use var_kind_def
    implicit none
    integer(dpI),intent(in):: isec,mode
    integer(dpI):: iele,e1,e0,j,ss,ee
    real   (dpR):: t(2)

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
        t(1:2)                      =  fv(isec)%turb(1,iele)
        fv(isec)%turb_gra(1:6,iele) =  0.0d0
        do j=sec(isec)%iA_ls(iele),sec(isec)%iA_ls(iele+1)-1
            ss  =  sec(isec)%jA_ls(1,j)
            ee  =  sec(isec)%jA_ls(2,j)
            fv(isec)%turb_gra(1:3,iele) =  fv(isec)%turb_gra(1:3,iele) &
                & +sec(isec)%coe_ls(1:3,j)*(fv(ss)%turb(1,ee)-t(1))
            fv(isec)%turb_gra(4:6,iele) =  fv(isec)%turb_gra(4:6,iele) &
                & +sec(isec)%coe_ls(1:3,j)*(fv(ss)%turb(2,ee)-t(2))
        end do
    end do

    return
    end subroutine fv_sst_get_gra_LS
!   ----------------------------------------------------------------------------
!   get gra for turbulence equation, Least-Square modified for memory.
!   ----------------------------------------------------------------------------
    subroutine fv_sst_get_gra_LSm(isec,mode)
    use var_kind_def
    implicit none
    integer(dpI),intent(in):: isec,mode
    integer(dpI):: iele,e1,e0,j,ss,ee
    real   (dpR):: dp(3),d,c(3),LHS(3,3),t(2)

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

        t(1:2)                      =  fv(isec)%turb(1,iele)
        fv(isec)%turb_gra(1:6,iele) =  0.0d0
        do j=sec(isec)%iA_ls(iele),sec(isec)%iA_ls(iele+1)-1
            ss      =  sec(isec)%jA_ls(1,j)
            ee      =  sec(isec)%jA_ls(2,j)
            dp(1:3) =  sec(ss)%cen(1:3,ee)-sec(isec)%cen(1:3,iele)
            d       =  1.0d0/sqrt(dp(1)*dp(1)+dp(2)*dp(2)+dp(3)*dp(3))**ls_weight_order
            c(1:3)  =  d*(LHS(1:3,1)*dp(1)+LHS(1:3,2)*dp(2)+LHS(1:3,3)*dp(3))

            fv(isec)%turb_gra(1:3,iele) =  fv(isec)%turb_gra(1:3,iele) &
                & +(fv(ss)%turb(1,ee)-t(1))*c(1:3)
            fv(isec)%turb_gra(4:6,iele) =  fv(isec)%turb_gra(4:6,iele) &
                & +(fv(ss)%turb(2,ee)-t(2))*c(1:3)
        end do
    end do

    return
    end subroutine fv_sst_get_gra_LSm
    end subroutine fv_SST_get_gra
!-------------------------------------------------------------------------------
!   get mut from k and omega.
!-------------------------------------------------------------------------------
    subroutine fv_SST_ko_to_mut(lev,mode)
    use var_kind_def
    use var_fv
    use var_mesh, only: mesh,sec
    implicit none
    integer(dpI),intent(in):: lev,mode
    integer(dpI):: isec,i,ele1,ele0

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
            fv(isec)%mu(2,i)=  fv(isec)%u(1,i)*fv(isec)%turb(1,i)*0.31d0 &
                & /max(0.31d0*fv(isec)%turb(2,i), fv(isec)%rhs(3,i))
        end do
    end do

    return
    end subroutine fv_SST_ko_to_mut
!-------------------------------------------------------------------------------
!   cal Rhs, turbulence.
!-------------------------------------------------------------------------------
    subroutine fv_SST_get_rhs(lev,is_LHS)
    use var_kind_def
    use var_fv
    use var_global, only: R23,I2
    use var_mesh
    use var_parallel
    use var_prec
    use var_turb
    use var_uns_cal, only: is_BDF_now,dt_uns,uns_iter,is_uns_initialized
    implicit none
    logical(dpL),intent(in):: is_LHS
    integer(dpI),intent(in):: lev
    integer(dpI):: isec,iele,im,sL,eL,sR,eR,iA(1000)
    real   (dpR):: n(5),uLf(5),uRf(5),tL(2),tR(2),mul,mut,nul,nut,gra(6),d(3),dt, &
                &  unL,unR,vol,dnw,rh,ux,uy,uz,vx,vy,vz,wx,wy,wz,d1,v1,v2,crs,cd, &
                &  s12,s13,s23,F1,F11,F2,S,O,Pk,Co,t(6),src(2),turk,turo,ncdu,sk, &
                &  so,f(2),vg(3),ftL(2,2),ftR(2,2),rtmp

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
!   get F1.
    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
    if(sec(isec)%is_bnd)    cycle
    do iele=1,sec(isec)%n_ele
        vol =  sec(isec)%vol( iele)
        dnw =  sec(isec)%dnw( iele)
        rh  =  fv (isec)%u (1,iele)
        mul =  fv (isec)%mu(1,iele)
        nul =  mul/rh
        turk=  fv (isec)%turb(1,iele)
        turo=  fv (isec)%turb(2,iele)
        d1  =  1.0d0/dnw

        v1  =  sqrt(turk)*d1/(turo*9.0d-2)
        v2  =  5.0d2*nul*d1*d1/turo
        crs =  fv(isec)%turb_gra(1,iele)*fv(isec)%turb_gra(4,iele) &
            & +fv(isec)%turb_gra(2,iele)*fv(isec)%turb_gra(5,iele) &
            & +fv(isec)%turb_gra(3,iele)*fv(isec)%turb_gra(6,iele)
        if(SST_version .eq. 1994) then
            CD  =  max(1.712d0*rh*crs/turo, 1.0d-20)
        else
            CD  =  max(1.712d0*rh*crs/turo, 1.0d-10)
        end if
        s12 =  max(v1, v2)
        s13 =  min(s12, 3.424d0*rh*turk*d1*d1/CD)
!       F1.
        fv(isec)%rhs(4,iele)=  tanh(s13**4)
    end do
    end do
!   get F1.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   source term.
    ftL =  0.0d0
    if(is_DDES_now) then
        call fv_DDES_get_src(lev, .true.)
        stop 'Error: SST does not support DDES right now.'
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
            nut =  mut/rh
            turk=  fv (isec)%turb(1,iele)
            turo=  fv (isec)%turb(2,iele)
            d1  =  1.0d0/dnw

            ux  =  fv(isec)%gra(4 ,iele)
            uy  =  fv(isec)%gra(5 ,iele)
            uz  =  fv(isec)%gra(6 ,iele)
            vx  =  fv(isec)%gra(7 ,iele)
            vy  =  fv(isec)%gra(8 ,iele)
            vz  =  fv(isec)%gra(9 ,iele)
            wx  =  fv(isec)%gra(10,iele)
            wy  =  fv(isec)%gra(11,iele)
            wz  =  fv(isec)%gra(12,iele)
            ncdu=  ux+vy+wz

            v1  =  sqrt(turk)*d1/(turo*9.0d-2)
            v2  =  5.0d2*nul*d1*d1/turo
            crs =  fv(isec)%turb_gra(1,iele)*fv(isec)%turb_gra(4,iele) &
                & +fv(isec)%turb_gra(2,iele)*fv(isec)%turb_gra(5,iele) &
                & +fv(isec)%turb_gra(3,iele)*fv(isec)%turb_gra(6,iele)
            F1  =  fv(isec)%rhs(4,iele)
            F11 =  1.0d0-F1

            s12 =  0.5d0*(uy+vx)
            s13 =  0.5d0*(uz+wx)
            s23 =  0.5d0*(vz+wy)
            S   =  sqrt(2.0d0*(ux*ux+vy*vy+wz*wz)+4.0d0*(s12*s12+s13*s13+s23*s23))
            O   =  sqrt((uy-vx)**2+(uz-wx)**2+(vz-wy)**2)

            rtmp= -ncdu/3.0d0
            Pk  = -R23*rh*turk
            Co  =  2.0d0*mut
            t(1)=  Co*(ux+rtmp)+Pk
            t(2)=  Co* s12
            t(3)=  Co* s13
            t(4)=  Co*(vy+rtmp)+Pk
            t(5)=  Co* s23
            t(6)=  Co*(wz+rtmp)+Pk
            if(SST_version .eq. 1994) then
                Pk  =  t(1)*ux+t(4)*vy+t(6)*wz+2.0d0*(t(2)*s12+t(3)*s13+t(5)*s23)
                Pk  =  min(Pk, 1.8d0*rh*turk*turo)
            else
                Pk  =  mut*S*S
                Pk  =  min(Pk, 0.9d0*rh*turk*turo)
            end if
            Co  =  1.712d0*F11*rh*crs/turo

            if(is_DDES_now) then
            else
                src(1)  =  Pk-9.0d-2*rh*turk*turo
            end if
            if(SST_version .eq. 1994) then
                ux  =  F1*0.5531667d0+F11*0.440354667d0
            else
                ux  =  F1*0.5555556d0+F11*0.44d0
            end if
            uy      =  F1*7.5d-2     +F11*8.28d-2
            src(2)  =  Pk*ux/nut-uy*rh*turo*turo+Co
            fv(isec)%rhs(1:2,iele)  =  fv(isec)%rhs(1:2,iele)-vol*src(1:2)

            uz  =  max(2.0d0*v1, v2)
            F2  =  tanh(uz*uz)
            if(SST_version .eq. 1994) then
                fv(isec)%rhs(3,iele)=  O*F2
            else
                fv(isec)%rhs(3,iele)=  S*F2
            end if

            if(is_LHS) then
                if(is_DDES_now) then
                else
!                   linearization based on the SST-V model.
                    ftL(1,1)=  9.0d-2*turo-min(ncdu, 0.0d0)*R23
                end if
                ftL(2,2)=  max(Co, 0.0d0)/turo+2.0d0*uy*turo

                eL  =  iele+iA(isec)-1
                call prec_add_eleR(fv_rans_prec, eL, eL, vol*ftL)
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
            tL(1:2) =  fv(sR)%turb(1:2,eR)
            tR      =  tL
            call DCOPY(6, fv(sR)%turb_gra(1,eR), 1, gra, 1)
            mul     =  fv(sR)%mu(1,eR)
            mut     =  fv(sR)%mu(2,eR)
            F1      =  fv(sL)%rhs(4,eL)
        else
            uLf(1:5)=  fv(sL)%u(1:5,eL)
            uRf(1:5)=  fv(sR)%u(1:5,eR)
            tL(1:2) =  fv(sL)%turb(1:2,eL)
            tR(1:2) =  fv(sR)%turb(1:2,eR)
            mul     =  0.5d0*(fv(sL)%mu(1,eL)+fv(sR)%mu(1,eR))
            mut     =  0.5d0*(fv(sL)%mu(2,eL)+fv(sR)%mu(2,eR))
            F1      =  0.5d0*(fv(sL)%rhs(4,eL)+fv(sR)%rhs(4,eR))

            if(.true.) then
                gra(1:3)    = (tR(1)-tL(1))*n(1:3)*n(5)
                gra(4:6)    = (tR(2)-tL(2))*n(1:3)*n(5)
            else
            end if
        end if

        F11 =  1.0d0-F1
        sk  =  F1*0.85d0+F11
        so  =  F1*0.5d0 +F11*0.856d0
        unL =  uLf(1)*((uLf(2)-vg(1))*n(1)+(uLf(3)-vg(2))*n(2)+(uLf(4)-vg(3))*n(3))*n(4)
        unR =  uRf(1)*((uRf(2)-vg(1))*n(1)+(uRf(3)-vg(2))*n(2)+(uRf(4)-vg(3))*n(3))*n(4)

        unL = (uLf(2)-vg(1))*n(1)+(uLf(3)-vg(2))*n(2)+(uLf(4)-vg(3))*n(3)
        unR = (uRf(2)-vg(1))*n(1)+(uRf(3)-vg(2))*n(2)+(uRf(4)-vg(3))*n(3)
        rh  =  0.5d0*(uLf(1)+uRf(1))

        f(1:2)  =  max(unL, 0.0d0)*uLf(1)*tL(1:2)*n(4) &
                & +min(unR, 0.0d0)*uRf(1)*tR(1:2)*n(4)
        f(1  )  =  f(1)-(mul+sk*mut)*(n(1)*gra(1)+n(2)*gra(2)+n(3)*gra(3))*n(4)
        f(2  )  =  f(2)-(mul+so*mut)*(n(1)*gra(4)+n(2)*gra(5)+n(3)*gra(6))*n(4)

        fv(sL)%rhs(1:2,eL)  =  fv(sL)%rhs(1:2,eL)+f(1:2)
        if(sec(sR)%is_int)  fv(sR)%rhs(1:2,eR)  =  fv(sR)%rhs(1:2,eR)-f(1:2)

        if(is_LHS) then
            eL  =  eL+iA(sL)-1
            eR  =  eR+iA(sR)-1

            ftL(1,1)=  max(unL, 0.0d0)*n(4)+(mul+sk*mut)*n(4)*n(5)/rh
            ftL(2,1)=  0.0d0
            ftL(1,2)=  0.0d0
            ftL(2,2)=  max(unL, 0.0d0)*n(4)+(mul+so*mut)*n(4)*n(5)/rh

            ftR(1,1)=  min(unR, 0.0d0)*n(4)-(mul+sk*mut)*n(4)*n(5)/rh
            ftR(2,1)=  0.0d0
            ftR(1,2)=  0.0d0
            ftR(2,2)=  min(unR, 0.0d0)*n(4)-(mul+so*mut)*n(4)*n(5)/rh

            call prec_add_eleR(fv_rans_prec, eL, eL, ftL)
            if(sec(sR)%is_int) then
                call prec_add_eleR(fv_rans_prec, eL, eR, ftR)
                call prec_add_eleR(fv_rans_prec, eR, eR,-ftR)
                call prec_add_eleR(fv_rans_prec, eR, eL,-ftL)
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
                    f(1:2)  =  fv(isec)%u(1,iele)*fv(isec)%turb(1:2,iele) &
                            & -fv(isec)%uns_turb(1:2,iele)
                else
                    f(1:2)  =  1.5d0*fv(isec)%u(1,iele)*fv(isec)%turb(1:2,iele) &
                            & -2.0d0*fv(isec)%uns_turb(1:2,iele) &
                            & +0.5d0*fv(isec)%uns_turb(3:4,iele)
                end if
                fv(isec)%rhs(1:2,iele)  =  fv(isec)%rhs(1:2,iele)+vol*f*dt

                if(is_LHS) then
                    eL  =  iele+iA(isec)-1
                    call prec_add_eleR(fv_rans_prec, eL, eL, 1.5d0*vol*dt*I2)
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
    end subroutine fv_SST_get_rhs
