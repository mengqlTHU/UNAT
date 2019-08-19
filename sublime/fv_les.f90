!-------------------------------------------------------------------------------
!   get the grid length scale used in DDES and LES.
!-------------------------------------------------------------------------------
    subroutine fv_get_gls
    use var_kind_def
    use var_global, only: err_mem,span,R13,is_2d_cal
    use var_mesh
    use var_turb
    implicit none
    logical(dpL):: is_2d_geo
    integer(dpI):: isec
    real   (dpR):: d3_to_d2(3,3)

    is_2d_geo   =  span(1) .ge. -1.0d2
    if(is_2d_geo) then
        call gen_coord_n3(span, d3_to_d2(1,1), d3_to_d2(1,2))
        d3_to_d2(1:3,3) =  span(1:3)
        call mat_inv(3, d3_to_d2)
    end if

    do isec=mesh(0)%sec_1,mesh(0)%sec_0
        if(.not. sec(isec)%is_int)  cycle
        if(.not. allocated(sec(isec)%gls)) &
            &  allocate(sec(isec)%gls(sec(isec)%n_ele), stat=err_mem)

        if(is_2d_geo) then
            call get_gls_2d(isec)
        elseif(length_scale .eq. length_max) then
            call get_gls_max(isec)
        elseif(length_scale .eq. length_Scotti) then
            call get_gls_Scotti(isec)
        elseif(length_scale .eq. length_Shur) then
            call get_gls_Shur(isec)
        elseif(length_scale .eq. length_SSM) then
            call get_gls_SSM(isec)
        elseif(length_scale .eq. length_LSQ) then
            call get_gls_LSQ(isec)
        else
            stop 'Error: length scale method not supported.'
        end if
    end do

    return
    contains
!       ------------------------------------------------------------------------
!       compute grid length scale, 2D element.
!       ------------------------------------------------------------------------
        subroutine get_gls_2d(isec)
        implicit none
        integer(dpI),intent(in):: isec
        integer(dpI):: i,iele
        real   (dpR):: p(3,8),p2(2,8),c(3),xyz2obb(9),d(3),L(3),R(3)

        do iele=1,sec(isec)%n_ele
            do i=1,sec(isec)%npe
                if(is_2d_cal) then
                    p(1:2,i)=  mesh(0)%xyz(1:2,sec(isec)%n2e(i,iele))
                    p(3  ,i)=  0.0d0
                else
                    p(1:3,i)=  mesh(0)%xyz(1:3,sec(isec)%n2e(i,iele))
                end if
                p2(1:2,i)   =  d3_to_d2(1:2,1)*(p(1,i)-p(1,1)) &
                            & +d3_to_d2(1:2,2)*(p(2,i)-p(2,1)) &
                            & +d3_to_d2(1:2,3)*(p(3,i)-p(3,1))
            end do
            call cal_obb(2, sec(isec)%npe, 2, p2, c, xyz2obb, L, R)
            d(1:2)  =  R(1:2)-L(1:2)
            sec(isec)%gls(iele) =  max(d(1), d(2))
        end do

        return
        end subroutine get_gls_2d
!       ------------------------------------------------------------------------
!       compute grid length scale, max.
!       ------------------------------------------------------------------------
        subroutine get_gls_max(isec)
        implicit none
        integer(dpI),intent(in):: isec
        integer(dpI):: i,iele
        real   (dpR):: p(3,8),c(3),xyz2obb(9),d(3),L(3),R(3)

        do iele=1,sec(isec)%n_ele
            do i=1,sec(isec)%npe
                if(is_2d_cal) then
                    p(1:2,i)=  mesh(0)%xyz(1:2,sec(isec)%n2e(i,iele))
                    p(3  ,i)=  0.0d0
                else
                    p(1:3,i)=  mesh(0)%xyz(1:3,sec(isec)%n2e(i,iele))
                end if
            end do
            call cal_obb(3, sec(isec)%npe, 3, p , c, xyz2obb, L, R)
            d(1:3)  =  R(1:3)-L(1:3)
            sec(isec)%gls(iele) =  maxval(d(1:3))
        end do

        return
        end subroutine get_gls_max
!       ------------------------------------------------------------------------
!       compute grid length scale, Scotti.
!       ------------------------------------------------------------------------
        subroutine get_gls_Scotti(isec)
        implicit none
        integer(dpI),intent(in):: isec
        integer(dpI):: i,iele
        real   (dpR):: p(3,8),d(3),L(3),R(3),c(3),xyz2obb(9),d1,d2,d3,f,a1,a2

        do iele=1,sec(isec)%n_ele
            do i=1,sec(isec)%npe
                if(is_2d_cal) then
                    p(1:2,i)=  mesh(0)%xyz(1:2,sec(isec)%n2e(i,iele))
                    p(3  ,i)=  0.0d0
                else
                    p(1:3,i)=  mesh(0)%xyz(1:3,sec(isec)%n2e(i,iele))
                end if
            end do
            call cal_obb(3, sec(isec)%npe, 3, p , c, xyz2obb, L, R)
            d(1:3)  =  R(1:3)-L(1:3)

!           CTR-1996.
            d1  =  minval(d(1:3))
            d3  =  maxval(d(1:3))
            if((d(1) .ge. d1) .and. (d(1) .le. d3)) then
                d2  =  d(1)
            elseif((d(2) .ge. d1) .and. (d(2) .le. d3)) then
                d2  =  d(2)
            else
                d2  =  d(3)
            end if
            a1  =  d1/d3
            a2  =  d2/d3
            f   =  log(a1)**2-log(a1)*log(a2)+log(a2)**2
            sec(isec)%gls(iele) = (d1*d2*d3)**R13*cosh(sqrt(4.0d0*f/2.7d1))
        end do

        return
        end subroutine get_gls_Scotti
!       ------------------------------------------------------------------------
!       compute grid length scale, Shur.
!       ------------------------------------------------------------------------
        subroutine get_gls_Shur(isec)
        implicit none
        integer(dpI),intent(in):: isec
        integer(dpI):: i,iele
        real   (dpR):: p(3,8),d(3),d3,c(3),L(3),R(3),xyz2obb(9)

        do iele=1,sec(isec)%n_ele
            do i=1,sec(isec)%npe
                if(is_2d_cal) then
                    p(1:2,i)=  mesh(0)%xyz(1:2,sec(isec)%n2e(i,iele))
                    p(3  ,i)=  0.0d0
                else
                    p(1:3,i)=  mesh(0)%xyz(1:3,sec(isec)%n2e(i,iele))
                end if
            end do
            call cal_obb(3, sec(isec)%npe, 3, p , c, xyz2obb, L, R)
            d(1:3)  =  R(1:3)-L(1:3)

            d3  =  maxval(d(1:3))
            sec(isec)%gls(iele) =  min(max(sec(isec)%gls_wn(iele), &
                &  0.15d0*max(sec(isec)%dnw(iele), d3)), d3)
        end do

        return
        end subroutine get_gls_Shur
!       ------------------------------------------------------------------------
!       compute grid length scale, SSM.
!       ------------------------------------------------------------------------
        subroutine get_gls_SSM(isec)
        implicit none
        integer(dpI),intent(in):: isec
        integer(dpI):: i,iele
        real   (dpR):: gls,p(3,8),d(3),L(3),R(3),c(3),xyz2obb(9),d1,d2,d3,a1,a2,f

        do iele=1,sec(isec)%n_ele
            do i=1,sec(isec)%npe
                if(is_2d_cal) then
                    p(1:2,i)=  mesh(0)%xyz(1:2,sec(isec)%n2e(i,iele))
                    p(3  ,i)=  0.0d0
                else
                    p(1:3,i)=  mesh(0)%xyz(1:3,sec(isec)%n2e(i,iele))
                end if
            end do
            call cal_obb(3, sec(isec)%npe, 3, p , c, xyz2obb, L, R)
            d(1:3)  =  R(1:3)-L(1:3)

!           The Scotti length scale.
            d1  =  minval(d(1:3))
            d3  =  maxval(d(1:3))
            if((d(1) .ge. d1) .and. (d(1) .le. d3)) then
                d2  =  d(1)
            elseif((d(2) .ge. d1) .and. (d(2) .le. d3)) then
                d2  =  d(2)
            else
                d2  =  d(3)
            end if
            a1  =  d1/d3
            a2  =  d2/d3
            f   =  log(a1)**2-log(a1)*log(a2)+log(a2)**2
            gls = (d1*d2*d3)**R13*cosh(sqrt(4.0d0*f/2.7d1))

!           The Shur length scale.
            d3  =  maxval(d(1:3))
            a2  =  min(max(sec(isec)%gls_wn(iele), 0.15d0*max(sec(isec)%dnw(iele), &
                &  d3)), d3)

!           The Shur-Scotti-Min length scale, AIAA-2017-4282.
            gls =  min(gls, a2)
            sec(isec)%gls(iele) =  min(gls, a2)
        end do

        return
        end subroutine get_gls_SSM
!       ------------------------------------------------------------------------
!       get the grid length scale used in DDES, Mockett.
!       ------------------------------------------------------------------------
        subroutine get_gls_Mockett(isec)
        use var_cgns
        use var_fv
        implicit none
        integer(dpI),intent(in):: isec
        integer(dpI):: iele,m,n
        real   (dpR):: L(3,8),ug(9),o(3),gls,d(3)

        if(sec(isec)%ele_type .eq. HEXA_8) then
            do iele=1,sec(isec)%n_ele
                forall(m=1:8)   L(1:3,m)=  mesh(0)%xyz(1:3,sec(isec)%n2e(m,iele))
                ug(1:9) =  fv(isec)%gra(4:12,iele)
                o(1)    =  ug(6)-ug(8)
                o(2)    =  ug(7)-ug(3)
                o(3)    =  ug(2)-ug(4)
                o       =  o/sqrt(o(1)**2+o(2)**2+o(3)**2)
                gls     =  0.0d0
                do m=1,8
                do n=1,8
                    if(m .eq. n)    cycle
                    d(1)=  o(2)*(L(3,m)-L(3,n))-o(3)*(L(2,m)-L(2,n))
                    d(2)= -o(1)*(L(3,m)-L(3,n))+o(3)*(L(1,m)-L(1,n))
                    d(3)=  o(1)*(L(2,m)-L(2,n))-o(2)*(L(1,m)-L(1,n))
                    gls =  max(gls, sqrt(d(1)**2+d(2)**2+d(3)**2))
                end do
                end do
                sec(isec)%gls(iele) =  1.025d0*gls/sqrt(3.0d0)
            end do
        else
        end if

        return
        end subroutine get_gls_Mockett
!       ------------------------------------------------------------------------
!       get the grid length scale used in DDES, LSQ.
!       ------------------------------------------------------------------------
        subroutine get_gls_LSQ(isec)
        use var_cgns
        use var_fv
        implicit none
        integer(dpI),intent(in):: isec
        integer(dpI):: iele,k,im
        real   (dpR):: v,J(3),n(3),G(3,3),GT(3,3),GTG(3,3),JGTG(3,3),a,b

        do iele=1,sec(isec)%n_ele
            v   =  1.0d0/sec(isec)%vol(iele)
            J   =  0.0d0
            do k=sec(isec)%iA_face_neighbour(iele),sec(isec)%iA_face_neighbour(iele+1)-1
                im  =  sec(isec)%jA_face_neighbour(k)
                if(im .gt. 0) then
                    n(1:3)  =  mesh(0)%mortar_n_vg(1:3, im)*mesh(0)%mortar_n_vg(4, im)
                else
                    n(1:3)  =  mesh(0)%mortar_n_vg(1:3,-im)*mesh(0)%mortar_n_vg(4,-im)
                end if
                n   =  0.5d0*n*v
                J   =  J+abs(n)
            end do
            J(1)=  1.0d0/J(1)
            J(2)=  1.0d0/J(2)
            if(.not. is_2d_cal) J(3)=  1.0d0/J(3)

            G (1,1:3)   =  fv(isec)%gra(4 :6 ,iele)
            G (2,1:3)   =  fv(isec)%gra(7 :9 ,iele)
            G (3,1:3)   =  fv(isec)%gra(10:12,iele)
            GT(1:3,1)   =  fv(isec)%gra(4 :6 ,iele)
            GT(1:3,2)   =  fv(isec)%gra(7 :9 ,iele)
            GT(1:3,3)   =  fv(isec)%gra(10:12,iele)
            GTG(1:3,1)  =  GT(1:3,1)*G(1,1)+GT(1:3,2)*G(2,1)+GT(1:3,3)*G(3,1)
            GTG(1:3,2)  =  GT(1:3,1)*G(1,2)+GT(1:3,2)*G(2,2)+GT(1:3,3)*G(3,2)
            GTG(1:3,3)  =  GT(1:3,1)*G(1,3)+GT(1:3,2)*G(2,3)+GT(1:3,3)*G(3,3)
            JGTG(1,1:3) =  J(1)*GTG(1,1:3)
            JGTG(2,1:3) =  J(2)*GTG(2,1:3)
            JGTG(3,1:3) =  J(3)*GTG(3,1:3)

            a   =  JGTG(1,1)**2+JGTG(2,1)**2+JGTG(3,1)**2 &
                & +JGTG(1,2)**2+JGTG(2,2)**2+JGTG(3,2)**2 &
                & +JGTG(1,3)**2+JGTG(2,3)**2+JGTG(3,3)**2
            b   =   GTG(1,1)**2+ GTG(2,1)**2+ GTG(3,1)**2 &
                & + GTG(1,2)**2+ GTG(2,2)**2+ GTG(3,2)**2 &
                & + GTG(1,3)**2+ GTG(2,3)**2+ GTG(3,3)**2

            sec(isec)%gls(iele) =  sqrt(a/b)
        end do

        return
        end subroutine get_gls_LSQ
    end subroutine fv_get_gls
!-------------------------------------------------------------------------------
!   get LES mut for FV in the parallel environment on lev(lev).
!-------------------------------------------------------------------------------
    subroutine fv_LES(lev)
    use var_kind_def
    use var_fv
    use var_mesh
    use var_parallel
    use var_turb, only: is_LES_now
    implicit none
    integer(dpI),intent(in):: lev
    integer(dpI):: isr,isec,iele,i,idx,ip_remote

    if(.not. is_LES_now)    return

!   ----------------------------------------------------------------------------
!   get mut from gradient for the elements used in the data exchange.
    call fv_get_LES_mut(lev, 1)
!   get mut from gradient for the elements used in the data exchange.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   prepare data.
    do isr=mesh(lev)%p2p_1,mesh(lev)%p2p_0
        idx =  1
        do i=1,p2p(isr)%n_ele_send_face
            isec=  p2p(isr)%id_ele_send(1,i)
            iele=  p2p(isr)%id_ele_send(2,i)
            p2p(isr)%rsend(idx) =  fv(isec)%mu(2,iele)
            idx =  idx+1
        end do
        p2p(isr)%n_send =  idx-1
        p2p(isr)%n_recv =  p2p(isr)%n_ele_recv_face
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
!   get mut from gradient for the elements not used in the data exchange.
    call fv_get_LES_mut(lev, 2)
!   get mut from gradient for the elements not used in the data exchange.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   turb boundary.
    call fv_bnd_LES_turb(lev)
!   trub boundary.
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
            if(iele .gt. p2p(isr)%n_ele_recv_face)  cycle
            fv(isec)%mu(2,i)=  p2p(isr)%rrecv(iele)
        end do
    end do
!   unpack data received.
!   ----------------------------------------------------------------------------

    return
    end subroutine fv_LES
!-------------------------------------------------------------------------------
!   modify(update) the turbulent variables on lev(lev).
!-------------------------------------------------------------------------------
    subroutine fv_get_LES_mut(lev,mode)
    use var_kind_def
    use var_mesh
    use var_turb
    use var_fv
    implicit none
    integer(dpI),intent(in):: lev,mode
    integer(dpI):: ele1,ele0,i,isec
    real   (dpR):: rh,vol,ug(9)

    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(sec(isec)%is_ghost .or. sec(isec)%is_bnd)    cycle
        if(mode .eq. 0) then
            ele1=  1
            ele0=  sec(isec)%n_ele
        elseif(mode.eq.1) then
            ele1=  1
            ele0=  sec(isec)%n_ele_b
        else
            ele1=  sec(isec)%n_ele_b+1
            ele0=  sec(isec)%n_ele
        end if
        do i=ele1,ele0
            rh      =  fv(isec)%u  (1   ,i)
            ug(1:9) =  fv(isec)%gra(4:12,i)
            vol     =  sec(isec)%vol(    i)
            call get_LES_mut(rh, ug, vol, fv(isec)%mu(2,i))
        end do
    end do

    return
    end subroutine fv_get_LES_mut
!-------------------------------------------------------------------------------
!   LES turb bnd modify.
!-------------------------------------------------------------------------------
    subroutine fv_bnd_LES_turb(lev)
    use var_kind_def
    use var_cgns
    use var_mesh
    use var_fv
    implicit none
    integer(dpI),intent(in):: lev
    integer(dpI):: im,sL,eL,sR,eR,bct
    real   (dpR):: rh,vol,mut,ug(9)

    do im=1,mesh(lev)%n_mortar_b
        sL  =  mesh(lev)%mortar_LR(1,im)
        eL  =  mesh(lev)%mortar_LR(2,im)
        sR  =  mesh(lev)%mortar_LR(3,im)
        eR  =  mesh(lev)%mortar_LR(4,im)
        bct =  sec(sR)%bct

        rh      =  fv(sR)%u  (1   ,eR)
        ug(1:9) =  fv(sR)%gra(4:12,eR)
        vol     =  sec(sL)%vol(    eL)
        call get_LES_mut(rh, ug, vol, mut)

        if(bct .eq. BCFarfield) then
        elseif(bct .eq. BCWallInviscid) then
        elseif(bct .eq. BCWallViscous) then
            mut =  0.0d0
        elseif(bct .eq. BCInflow) then
        elseif(bct .eq. BCOutflow) then
        end if
        fv(sR)%mu(2,eR) =  mut
    end do

    return
    end subroutine fv_bnd_LES_turb
!-------------------------------------------------------------------------------
!   cal the LES mut.
!-------------------------------------------------------------------------------
    subroutine get_LES_mut(rh,ug,vol,mut)
    use var_kind_def
    use var_global, only: pi,R13
    use var_turb
    implicit none
    integer(dpI):: i,j,k
    real   (dpR),parameter :: eps=1.0d-16
    real   (dpR):: rh,vol,mut,rtmp,S(6),S_,flength,Sd(3,3),Sd_,ug(9),tt1,tt2,tt3, &
                &  alp1,alp2,alp3,sig1,sig2,sig3,alp12,dsig

    if(LES_model .eq. LES_sSMARG) then
        S(1)=  ug(1)
        S(2)=  0.5d0*(ug(2)+ug(4))
        S(3)=  0.5d0*(ug(3)+ug(7))
        S(4)=  ug(5)
        S(5)=  0.5d0*(ug(6)+ug(8))
        S(6)=  ug(9)
        S_  =  S(1)*S(1)+2d0*S(2)*S(2)+2d0*S(3)*S(3)+S(4)*S(4)+2d0*S(5)*S(5)+S(6)*S(6)
        S_  =  sqrt(2.0d0*S_)
        flength =  vol**R13
        mut     =  rh*flength*flength*cs_SMARG*cs_SMARG*S_
    elseif(LES_model .eq. LES_WALE) then
        S(1)=  ug(1)
        S(2)=  0.5d0*(ug(2)+ug(4))
        S(3)=  0.5d0*(ug(3)+ug(7))
        S(4)=  ug(5)
        S(5)=  0.5d0*(ug(6)+ug(8))
        S(6)=  ug(9)
        S_  =  S(1)*S(1)+2d0*S(2)*S(2)+2d0*S(3)*S(3)+S(4)*S(4)+2d0*S(5)*S(5)+S(6)*S(6)
        Sd_ =  0.0d0
        rtmp=  0.0d0
        Sd  =  0.0d0
        do i=1,3
        do j=1,3
            rtmp=  rtmp+ug(3*(i-1)+j)*ug(3*(j-1)+i)
            do k=1,3
                Sd(i,j) =  Sd(i,j)+ug(3*(i-1)+k)*ug(3*(k-1)+i)+ug(3*(j-1)+k)*ug(3*(k-1)+j)
            end do
        end do
        end do
        Sd  =  Sd*0.5d0
        rtmp=  rtmp/3d0
        Sd(1,1) =  Sd(1,1)+rtmp
        Sd(2,2) =  Sd(2,2)+rtmp
        Sd(3,3) =  Sd(3,3)+rtmp
        do i=1,3
        do j=1,3
            Sd_ =  Sd_+Sd(i,j)*Sd(i,j)
        end do
        end do
        flength =  vol**R13
        if(abs(S_) .lt. eps) then
            Sd_ =  0.0d0
            S_  =  1.0d0
        endif
        mut =  rh*cw_WALE*cw_WALE*flength*flength*Sd_**1.5d0/(S_**2.5d0+Sd_**1.25d0)
    elseif(LES_model .eq. LES_VSS) then
        S(1)=  ug(1)
        S(2)=  0.5d0*(ug(2)+ug(4))
        S(3)=  0.5d0*(ug(3)+ug(7))
        S(4)=  ug(5)
        S(5)=  0.5d0*(ug(6)+ug(8))
        S(6)=  ug(9)
        S_  =  S(1)*S(1)+2d0*S(2)*S(2)+2d0*S(3)*S(3)+S(4)*S(4)+2d0*S(5)*S(5)+S(6)*S(6)
        Sd_ = (S(2)*S(2)+S(3)*S(3)+S(5)*S(5))*(ug(1)*ug(1)+ug(5)*ug(5)+ug(9)*ug(9))
        flength =  vol**R13
        if(abs(S_).lt.eps) then
            Sd_ =  0.0d0
            S_  =  1.0d0
        endif
        mut =  rh*cr_VSS*cr_VSS*flength*flength*Sd_**1.5d0/S_**2.5d0
    elseif(LES_model .eq. LES_sigma) then
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

        flength =  vol**R13
        mut     =  rh*cr_sigma*cr_sigma*flength*flength*dsig
    else
        stop 'Error: LES model not supported.'
    end if

    return
    end subroutine get_LES_mut
!-------------------------------------------------------------------------------
!   get the shock sensor for FV in the parallel environment.
!-------------------------------------------------------------------------------
    subroutine fv_get_shock_sensor(lev)
    use var_kind_def
    use var_fv
    use var_mesh
    use var_parallel
    implicit none
    integer(dpI),intent(in):: lev
    integer(dpI):: isr,isec,iele,i,idx,ip_remote

    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(sec(isec)%is_bnd)    cycle
        i   =  sec(isec)%n_ele
        if(.not. allocated(fv(isec)%shock_sensor))  allocate(fv(isec)%shock_sensor(i))
    end do

!   ----------------------------------------------------------------------------
!   get shock sensor for the elments used in the data exchange
    call fv_get_Ducros_sensor(lev, 1)
!   get shock sensor for the elments used in the data exchange
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   prepare data.
    do isr=mesh(lev)%p2p_1,mesh(lev)%p2p_0
        idx =  1
        do i=1,p2p(isr)%n_ele_send
            isec=  p2p(isr)%id_ele_send(1,i)
            iele=  p2p(isr)%id_ele_send(2,i)
            call DCOPY(1, fv(isec)%shock_sensor(iele), 1, p2p(isr)%rsend(idx), 1)
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
!   get shock sensor for the elments not used in the data exchange
    call fv_get_Ducros_sensor(lev, 2)
!   get shock sensor for the elments not used in the data exchange
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
            fv(isec)%shock_sensor(i)=  p2p(isr)%rrecv(iele)
        end do
    end do
!   unpack data received.
!   ----------------------------------------------------------------------------

    return
    contains
!       ------------------------------------------------------------------------
!       get the Ducros shock sensor.
!       ------------------------------------------------------------------------
        subroutine fv_get_Ducros_sensor(lev,mode)
        use var_kind_def
        use var_fv
        use var_global, only: n_dim,uref
        use var_mesh
        use var_slv, only: is_vis_cal
        use var_turb, only: is_RANS,is_LES_now
        implicit none
        integer(dpI),intent(in):: lev,mode
        integer(dpI):: ele1,ele0,i,isec
        real   (dpR):: ug(9),div,div2,vor(3),vor2,L,mu,Re,e
        real   (dpR),parameter :: eps=1.0d-6,s0=2.0d-2

        Re  =  1.0d20
        do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
            if(sec(isec)%is_ghost .or. sec(isec)%is_bnd)    cycle
            if(mode .eq. 0) then
                ele1=  1
                ele0=  sec(isec)%n_ele
            elseif(mode.eq.1) then
                ele1=  1
                ele0=  sec(isec)%n_ele_b
            else
                ele1=  sec(isec)%n_ele_b+1
                ele0=  sec(isec)%n_ele
            end if
            do i=ele1,ele0
                L       =  sec(isec)%vol(i)**(1.0d0/real(n_dim, dpR))
                ug(1:9) =  fv(isec)%gra(4:12,i)
                div     =  ug(1)+ug(5)+ug(9)
                div2    =  div*div
                vor(1)  =  ug(8)-ug(6)
                vor(2)  =  ug(3)-ug(7)
                vor(3)  =  ug(4)-ug(2)
                vor2    =  vor(1)*vor(1)+vor(2)*vor(2)+vor(3)*vor(3)
                if(is_vis_cal) then
                    mu  =  fv(isec)%mu(1,i)
                    if(is_RANS .or. is_LES_now) mu  =  mu+fv(isec)%mu(2,i)
                    Re  =  fv(isec)%u(1,i)*uref*L/mu
                end if
                e   = (eps+1.0d0/Re)*(uref/L)**2

!               fv(isec)%shock_sensor(i)=  sign(1.0d0, div)*div2/(div2+vor2+e)
                fv(isec)%shock_sensor(i)=  max(s0, div2/(div2+vor2+e))
            end do
        end do

        return
        end subroutine fv_get_Ducros_sensor
    end subroutine fv_get_shock_sensor
!-------------------------------------------------------------------------------
!   get the blending function.
!-------------------------------------------------------------------------------
    subroutine fv_get_blending_vortex(gls,d,u,gu,nu,nut,sigma)
    use var_kind_def
    use var_global, only: uref,L_ref
    implicit none
    real(dpR),parameter:: min_sigma =  2.5d-2
    real(dpR),intent(in):: gu(*),d,u(*),nu,nut,gls
    real(dpR):: sigma,oo,s12,s13,s23,ss,B,g,K,L_t,A,div,div2,vor(3),vor2,Re,e,vortex

    if(.true.) then
!       DDES approach, Reynolds number.
        oo  =  sqrt((gu(2)-gu(4))**2+(gu(3)-gu(7))**2+(gu(6)-gu(8))**2)
        s12 =  0.5d0*(gu(2)+gu(4))
        s13 =  0.5d0*(gu(3)+gu(7))
        s23 =  0.5d0*(gu(6)+gu(8))
        ss  =  sqrt(2.0d0*(gu(1)**2+gu(5)**2+gu(9)**2)+4.0d0*(s12*s12+s13*s13+s23*s23))

        B       =  2.0d0*oo*max(oo, ss)/max(0.5d0*(ss*ss+oo*oo), 1.0d-10*(uref/L_ref)**2)
        g       =  tanh(B**4)
        K       =  max(sqrt(0.5d0*(ss*ss+oo*oo)), 0.1d0*uref/L_ref)
        L_t     =  sqrt((nu+nut)/(K*0.09d0**1.5d0))
        A       =  max(0.65d0*gls/(L_t*g)-0.5d0, 0.0d0)
        vortex  =  tanh(A**3)
        Re      =  sqrt(u(1)**2+u(2)**2+u(3)**2)*d/(nu+nut)
        Re      =  max(0.0d0, 1.0d0-2.0d0/Re)
        sigma   =  max(min_sigma, min(vortex, Re))
    else
!       Ducros approach.
        div     =  gu(1)+gu(5)+gu(9)
        div2    =  div*div
        vor(1)  =  gu(8)-gu(6)
        vor(2)  =  gu(3)-gu(7)
        vor(3)  =  gu(4)-gu(2)
        vor2    =  vor(1)*vor(1)+vor(2)*vor(2)+vor(3)*vor(3)
        Re      =  uref*gls/(nu+nut)
        e       = (1.0d-6+1.0d0/Re)*(uref/gls)**2
        sigma   =  max(min_sigma, div2/(div2+vor2+e))
    end if

    return
    end subroutine fv_get_blending_vortex
!-------------------------------------------------------------------------------
!   get the geometry based blending function.
!-------------------------------------------------------------------------------
    subroutine fv_get_blending_mesh(im,sigma)
    use var_kind_def
    use var_mesh
    implicit none
    integer(dpI),intent(in):: im
    integer(dpI):: sL,eL,sR,eR
    real   (dpR):: sigma,dc(3),d(3)

    sL      =  mesh(0)%mortar_LR(1,im)
    eL      =  mesh(0)%mortar_LR(2,im)
    sR      =  mesh(0)%mortar_LR(3,im)
    eR      =  mesh(0)%mortar_LR(4,im)
    dc(1:3) =  sec(sL)%cen(1:3,eL)-sec(sR)%cen(1:3,eR)
    d (1:3) = (sec(sL)%cen(1:3,eL)+sec(sR)%cen(1:3,eR))*0.5d0-mesh(0)%mortar_cen(1:3,im)
    sigma   =  sqrt(d(1)**2+d(2)**2+d(3)**2)/sqrt(dc(1)**2+dc(2)**2+dc(3)**2)

    return
    end subroutine fv_get_blending_mesh
!-------------------------------------------------------------------------------
!   output low dissipation parameter.
!-------------------------------------------------------------------------------
    subroutine fv_output_LD(lev)
    use var_kind_def
    use var_fv
    use var_global, only: uref
    use var_mesh
    use var_temp, only: uLf,uRf,vgf,fac_1d
    use var_turb, only: is_DDES_now,is_LES_now
    implicit none
    real   (dpR),parameter:: min_sigma =  2.5d-2
    integer(dpI),intent(in):: lev
    integer(dpI):: isec,im,sL,eL,sR,eR
    real   (dpR):: gls,g(9),nu,nut,sigma,rhL,rhR,v(3),div2,Re,vor(3),vor2,e,d(3)

    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(sec(isec)%is_int)    fv(isec)%rhs=  0.0d0
    end do

    do im=1+mesh(lev)%n_mortar_b,mesh(lev)%n_mortar
        sL  =  mesh(lev)%mortar_LR(1,im)
        eL  =  mesh(lev)%mortar_LR(2,im)
        sR  =  mesh(lev)%mortar_LR(3,im)
        eR  =  mesh(lev)%mortar_LR(4,im)

!       Ducros approach.
        gls =  0.5d0*(sec(sL)%gls(eL)+sec(sR)%gls(eR))
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
        sigma   =  max(min_sigma, div2/(div2+vor2+e))

!       mesh based approach.
        v(1:3)  =  sec(sL)%cen(1:3,eL)-sec(sR)%cen(1:3,eR)
        d(1:3)  = (sec(sL)%cen(1:3,eL)+sec(sR)%cen(1:3,eR))*0.5d0 &
                & -mesh(0)%mortar_cen(1:3,im)
        sigma   =  max(sigma, sqrt(d(1)**2+d(2)**2+d(3)**2)/sqrt(v(1)**2+v(2)**2+v(3)**2))

        fv(sL)%rhs(1,eL)=  max(fv(sL)%rhs(1,eL), sigma)
        if(sec(sR)%is_int)  fv(sR)%rhs(1,eR)=  max(fv(sR)%rhs(1,eR), sigma)
    end do

    do im=1+mesh(lev)%n_mortar_b,mesh(lev)%n_mortar
        sL  =  mesh(lev)%mortar_LR(1,im)
        eL  =  mesh(lev)%mortar_LR(2,im)
        sR  =  mesh(lev)%mortar_LR(3,im)
        eR  =  mesh(lev)%mortar_LR(4,im)
        fac_1d(1:4,1)   =  mesh(lev)%mortar_n_vg(1:4,im)
        vgf   (1:3,1)   =  mesh(lev)%mortar_n_vg(6:8,im)

        uLf(1:5,1)  =  fv(sL)%u(1:5,eL)
        uRf(1:5,1)  =  fv(sR)%u(1:5,eR)
        v  (1:3  )  =  0.5d0*(uLf(2:4,1)+uRf(2:4,1))

        gls =  0.5d0*(sec(sL)%gls(eL)+sec(sR)%gls(eR))
        rhL =  1.0d0/fv(sL)%u(1,eL)
        rhR =  1.0d0/fv(sR)%u(1,eR)
        nu  =  0.5d0*(fv(sL)%mu(1,eL)*rhL+fv(sR)%mu(1,eR)*rhR)
        if(is_LES_now .or. is_DDES_now) then
            nut =  0.5d0*(fv(sL)%mu(2,eL)*rhL+fv(sR)%mu(2,eR)*rhR)
        else
            nut =  0.0d0
        end if
        g(1:9)  =  0.5d0*(fv(sL)%gra(4:12,eL)+fv(sR)%gra(4:12,eR))
        call fv_get_blending_vortex(gls,mesh(lev)%mortar_n_vg(5,im), v,g,nu,nut, sigma)

        fv(sL)%rhs(2,eL)=  max(fv(sL)%rhs(2,eL), sigma)
        if(sec(sR)%is_int)  fv(sR)%rhs(2,eR)=  max(fv(sR)%rhs(2,eR), sigma)
    end do

    call fv_wr_tec_cc(4)

    return
    end subroutine fv_output_LD
