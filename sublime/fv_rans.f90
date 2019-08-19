!-------------------------------------------------------------------------------
!   bnd.
!-------------------------------------------------------------------------------
    subroutine fv_rans_bnd(lev)
    use var_kind_def
    use var_bndv
    use var_cgns
    use var_fv, only: fv
    use var_mesh
    use var_turb
    implicit none
    integer(dpI),intent(in):: lev
    integer(dpI):: im,sL,eL,sR,eR,bct
    real   (dpR):: un,t(2),u(5),transition

    if(RANS_model .le. RANS_WA) then
        do im=1,mesh(lev)%n_mortar_b
            sL      =  mesh(lev)%mortar_LR(1,im)
            eL      =  mesh(lev)%mortar_LR(2,im)
            sR      =  mesh(lev)%mortar_LR(3,im)
            eR      =  mesh(lev)%mortar_LR(4,im)
            bct     =  sec(sR)%bct
            u(1:5)  =  fv(sR)%u(1:5,eR)

            if(bct .eq. BCWallInviscid) then
                t(1)=  fv(sL)%turb(1,eL)
            elseif(bct .eq. BCWallViscous) then
                t(1)=  0.0d0
            else
                un  = (u(2)-mesh(lev)%mortar_n_vg(6,im))*mesh(lev)%mortar_n_vg(1,im) &
                    &+(u(3)-mesh(lev)%mortar_n_vg(7,im))*mesh(lev)%mortar_n_vg(2,im) &
                    &+(u(4)-mesh(lev)%mortar_n_vg(8,im))*mesh(lev)%mortar_n_vg(3,im)
                if(un .gt. 0.0d0) then
                    t(1)=  fv(sL)%turb(1,eL)
                else
                    t(1)=  nut_ref
                end if
            end if
            fv(sR)%turb(1,eR)   =  t(1)
        end do
    elseif(RANS_model .lt. RANS_SST) then
        if(RANS_model .eq. RANS_SAM) then
            transition  =  1.0d0
        else
            transition  =  0.0d0
        end if

        do im=1,mesh(lev)%n_mortar_b
            sL      =  mesh(lev)%mortar_LR(1,im)
            eL      =  mesh(lev)%mortar_LR(2,im)
            sR      =  mesh(lev)%mortar_LR(3,im)
            eR      =  mesh(lev)%mortar_LR(4,im)
            bct     =  sec(sR)%bct
            u(1:5)  =  fv(sR)%u(1:5,eR)

            if(bct .eq. BCWallInviscid) then
                t(1)=  fv(sL)%turb(1,eL)
                t(2)=  fv(sL)%turb(2,eL)
            elseif(bct .eq. BCWallViscous) then
                t(1)=  0.0d0
                t(2)=  fv(sL)%turb(2,eL)
            else
                un  = (u(2)-mesh(lev)%mortar_n_vg(6,im))*mesh(lev)%mortar_n_vg(1,im) &
                    &+(u(3)-mesh(lev)%mortar_n_vg(7,im))*mesh(lev)%mortar_n_vg(2,im) &
                    &+(u(4)-mesh(lev)%mortar_n_vg(8,im))*mesh(lev)%mortar_n_vg(3,im)
                if(un .gt. 0.0d0) then
                    t(1)=  fv(sL)%turb(1,eL)
                    t(2)=  fv(sL)%turb(2,eL)
                else
                    t(1)=  nut_ref
                    t(2)=  transition
                end if
            end if
            fv(sR)%turb(1:2,eR) =  t(1:2)
        end do
    else
        do im=1,mesh(lev)%n_mortar_b
            sL      =  mesh(lev)%mortar_LR(1,im)
            eL      =  mesh(lev)%mortar_LR(2,im)
            sR      =  mesh(lev)%mortar_LR(3,im)
            eR      =  mesh(lev)%mortar_LR(4,im)
            bct     =  sec(sR)%bct
            u(1:5)  =  fv(sR)%u(1:5,eR)

            if(bct .eq. BCWallInviscid) then
                t(1)=  fv(sL)%turb(1,eL)
                t(2)=  fv(sL)%turb(2,eL)
            elseif(bct .eq. BCWallViscous) then
                t(1)=  0.0d0
                t(2)=  8.0d3*fv(sL)%mu(1,eL)/(fv(sL)%u(1,eL)*sec(sL)%dnw(eL)**2)
            else
                un  = (u(2)-mesh(lev)%mortar_n_vg(6,im))*mesh(lev)%mortar_n_vg(1,im) &
                    &+(u(3)-mesh(lev)%mortar_n_vg(7,im))*mesh(lev)%mortar_n_vg(2,im) &
                    &+(u(4)-mesh(lev)%mortar_n_vg(8,im))*mesh(lev)%mortar_n_vg(3,im)
                if(un .gt. 0.0d0) then
                    t(1)=  fv(sL)%turb(1,eL)
                    t(2)=  fv(sL)%turb(2,eL)
                else
                    t(1)=  k_ref
                    t(2)=  o_ref
                end if
            end if
            fv(sR)%turb(1:2,eR) =  t(1:2)
        end do
    end if

    return
    end subroutine fv_rans_bnd
!-------------------------------------------------------------------------------
!   boundary treatment for FV in the parallel environment.
!-------------------------------------------------------------------------------
    subroutine fv_rans_bnd_parallel(lev,is_nut_to_mut)
    use var_kind_def
    use var_fv
    use var_mesh
    use var_parallel
    use var_turb
    implicit none
    logical(dpL),intent(in):: is_nut_to_mut
    integer(dpI),intent(in):: lev
    integer(dpI):: isr,isec,iele,i,idx,ip_remote

!   ----------------------------------------------------------------------------
!   get mut for the elements used in the data exchange.
    if(is_nut_to_mut) then
        if(RANS_model .eq. RANS_WA) then
            call fv_WA_nut_to_mut(lev, 1)
        elseif(RANS_model .lt. RANS_SST) then
            call fv_SA_nut_to_mut(lev, 1)
        else
            call fv_SST_ko_to_mut(lev, 1)
        end if
    end if
!   get mut for the elements used in the data exchange.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   prepare data.
    if(RANS_model .le. RANS_WA) then
        do isr=mesh(lev)%p2p_1,mesh(lev)%p2p_0
            idx =  1
            do i=1,p2p(isr)%n_ele_send
                isec=  p2p(isr)%id_ele_send(1,i)
                iele=  p2p(isr)%id_ele_send(2,i)
                p2p(isr)%rsend(idx  )   =  fv(isec)%turb(1,iele)
                p2p(isr)%rsend(idx+1)   =  fv(isec)%mu  (2,iele)
                idx =  idx+2
            end do
            p2p(isr)%n_send =  idx-1
            p2p(isr)%n_recv =  2*p2p(isr)%n_ele_recv
        end do
    elseif(RANS_model .lt. RANS_SST) then
        do isr=mesh(lev)%p2p_1,mesh(lev)%p2p_0
            idx =  1
            do i=1,p2p(isr)%n_ele_send
                isec=  p2p(isr)%id_ele_send(1,i)
                iele=  p2p(isr)%id_ele_send(2,i)
                p2p(isr)%rsend(idx  )   =  fv(isec)%turb(1,iele)
                p2p(isr)%rsend(idx+1)   =  fv(isec)%turb(2,iele)
                p2p(isr)%rsend(idx+2)   =  fv(isec)%mu  (2,iele)
                idx =  idx+3
            end do
            p2p(isr)%n_send =  idx-1
            p2p(isr)%n_recv =  3*p2p(isr)%n_ele_recv
        end do
    else
        do isr=mesh(lev)%p2p_1,mesh(lev)%p2p_0
            idx =  1
            do i=1,p2p(isr)%n_ele_send
                isec=  p2p(isr)%id_ele_send(1,i)
                iele=  p2p(isr)%id_ele_send(2,i)
                p2p(isr)%rsend(idx  )   =  fv(isec)%turb(1,iele)
                p2p(isr)%rsend(idx+1)   =  fv(isec)%turb(2,iele)
                p2p(isr)%rsend(idx+2)   =  fv(isec)%mu  (2,iele)
                idx =  idx+3
            end do
            p2p(isr)%n_send =  idx-1
            p2p(isr)%n_recv =  3*p2p(isr)%n_ele_recv
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
!   physical boundary.
    call fv_rans_bnd(lev)
!   physical boundary.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   get mut for the elements not used in the data exchange.
    if(is_nut_to_mut) then
        if(RANS_model .eq. RANS_WA) then
            call fv_WA_nut_to_mut(lev, 2)
        elseif(RANS_model .lt. RANS_SST) then
            call fv_SA_nut_to_mut(lev, 2)
        else
            call fv_SST_ko_to_mut(lev, 2)
        end if
    end if
!   get mut for the elements not used in the data exchange.
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
    if(RANS_model .le. RANS_WA) then
        do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
            if(.not. sec(isec)%is_ghost)    cycle
            do i=1,sec(isec)%n_ele
                isr =  sec(isec)%id_recv(1,i)
                iele=  sec(isec)%id_recv(2,i)
                fv(isec)%turb(1,i)  =  p2p(isr)%rrecv(2*iele-1)
                fv(isec)%mu  (2,i)  =  p2p(isr)%rrecv(2*iele  )
            end do
        end do
    elseif(RANS_model .lt. RANS_SST) then
        do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
            if(.not. sec(isec)%is_ghost)    cycle
            do i=1,sec(isec)%n_ele
                isr =  sec(isec)%id_recv(1,i)
                iele=  sec(isec)%id_recv(2,i)
                fv(isec)%turb(1,i)  =  p2p(isr)%rrecv(3*iele-2)
                fv(isec)%turb(2,i)  =  p2p(isr)%rrecv(3*iele-1)
                fv(isec)%mu  (2,i)  =  p2p(isr)%rrecv(3*iele  )
            end do
        end do
    else
        do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
            if(.not. sec(isec)%is_ghost)    cycle
            do i=1,sec(isec)%n_ele
                isr =  sec(isec)%id_recv(1,i)
                iele=  sec(isec)%id_recv(2,i)
                fv(isec)%turb(1,i)  =  p2p(isr)%rrecv(3*iele-2)
                fv(isec)%turb(2,i)  =  p2p(isr)%rrecv(3*iele-1)
                fv(isec)%mu  (2,i)  =  p2p(isr)%rrecv(3*iele  )
            end do
        end do
    end if
!   unpack data received.
!   ----------------------------------------------------------------------------

    if(RANS_model .lt. RANS_SST) then
        call fv_rans_get_gra(lev)
    else
        call fv_SST_get_gra (lev)
    end if

    return
    end subroutine fv_rans_bnd_parallel
!-------------------------------------------------------------------------------
!   setup the implicit solver, FV, RANS.
!-------------------------------------------------------------------------------
    subroutine fv_rans_set_slv(lev,p)
    use var_kind_def
    use var_fv
    use var_global, only: err_mem
    use var_mesh
    use var_prec
    use var_turb
    implicit none
    integer(dpI),intent(in):: lev
    integer(dpI):: isec
    type(type_prec):: p

    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(sec(isec)%is_ghost .or. sec(isec)%is_bnd)    cycle
        if(.not. allocated(fv(isec)%LHS_s)) &
            &  allocate(fv(isec)%LHS_s(sec(isec)%n_ele), stat=err_mem)
    end do
    if(RANS_model .eq. RANS_SA) then
        call fv_set_prec(lev, p, 1, 1, .false.)
    elseif(RANS_model .eq. RANS_WA) then
        call fv_set_prec(lev, p, 1, 1, .true.)
    elseif(RANS_model .eq. RANS_SAM) then
        call fv_set_prec(lev, p, 2, 1, .false.)
    elseif(RANS_model .eq. RANS_SAC) then
        call fv_set_prec(lev, p, 2, 1, .false.)
    elseif(RANS_model .eq. RANS_SST) then
        call fv_set_prec(lev, p, 2, 1, .false.)
    else
        stop 'Error: RANS model not supported.'
    end if
    call prec_set_ordering(p)

    return
    end subroutine fv_rans_set_slv
!-------------------------------------------------------------------------------
!   cal gradient.
!-------------------------------------------------------------------------------
    subroutine fv_rans_get_gra(lev)
    use var_kind_def
    use var_fv
    use var_mesh
    use var_parallel
    use var_slv
    implicit none
    integer(dpI),intent(in):: lev
    integer(dpI):: isec,iele,i,ip_remote,isr,LDA,idx,im,sL,eL,sR,eR
    real   (dpR):: g(3),n(3),d,dt,gn

!   ----------------------------------------------------------------------------
!   cal gradient for the elements involved in the data exchange.
    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(sec(isec)%is_int)    call get_gra(isec, 1)
    end do
!   cal gradient for the elements involved in the data exchange.
!   ----------------------------------------------------------------------------

    LDA =  3

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
            call per_rot_vec(12,sec(isec)%per_path(1,i), &
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

        g(1:3)  =  fv(sL)%turb_gra(1:3,eL)
        dt      =  fv(sL)%turb(1,eL)-fv(sR)%turb(1,eR)
        gn      =  g(1)*n(1)+g(2)*n(2)+g(3)*n(3)
        fv(sR)%turb_gra(1:3,eR) =  g(1:3)+(dt/d-gn)*n(1:3)
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
        call fv_sa_get_gra_gg (isec, mode)
    elseif(gradient_method .eq. gradient_LS) then
        call fv_sa_get_gra_LS (isec, mode)
    else
        call fv_sa_get_gra_LSm(isec, mode)
    end if

    return
    end subroutine get_gra
!   ----------------------------------------------------------------------------
!   get gra for turbulence equation, Green-Gauss.
!   ----------------------------------------------------------------------------
    subroutine fv_sa_get_gra_gg(isec,mode)
    use var_kind_def
    implicit none
    integer(dpI),intent(in):: isec,mode
    integer(dpI):: iele,e1,e0,j,im,s,e
    real   (dpR):: tL,vL,n(3),tm,vR

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
        tL  =  fv(isec)%turb(1,iele)
        vL  =  1.0d0/sec(isec)%vol(iele)
        fv(isec)%turb_gra(1:3,iele) =  0.0d0
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
                tm  =  fv(s)%turb(1,e)
            else
                vR  =  1.0d0/sec(s)%vol(e)
                tm  = (vL*tL+vR*fv(s)%turb(1,e))/(vL+vR)
            end if
            fv(isec)%turb_gra(1:3,iele) =  fv(isec)%turb_gra(1:3,iele)+tm*n(1:3)
        end do
        fv(isec)%turb_gra(1:3,iele) =  fv(isec)%turb_gra(1:3,iele)*vL
    end do

    return
    end subroutine fv_sa_get_gra_gg
!   ----------------------------------------------------------------------------
!   get gra for turbulence equation, Least-Square.
!   ----------------------------------------------------------------------------
    subroutine fv_sa_get_gra_LS(isec,mode)
    use var_kind_def
    implicit none
    integer(dpI),intent(in):: isec,mode
    integer(dpI):: iele,e1,e0,j,ss,ee
    real   (dpR):: t

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
        t                           =  fv(isec)%turb(1,iele)
        fv(isec)%turb_gra(1:3,iele) =  0.0d0
        do j=sec(isec)%iA_ls(iele),sec(isec)%iA_ls(iele+1)-1
            ss  =  sec(isec)%jA_ls(1,j)
            ee  =  sec(isec)%jA_ls(2,j)
            fv(isec)%turb_gra(1:3,iele) =  fv(isec)%turb_gra(1:3,iele) &
                & +sec(isec)%coe_ls(1:3,j)*(fv(ss)%turb(1,ee)-t)
        end do
    end do

    return
    end subroutine fv_sa_get_gra_LS
!   ----------------------------------------------------------------------------
!   get gra for turbulence equation, Least-Square modified for memory.
!   ----------------------------------------------------------------------------
    subroutine fv_sa_get_gra_LSm(isec,mode)
    use var_kind_def
    implicit none
    integer(dpI),intent(in):: isec,mode
    integer(dpI):: iele,e1,e0,j,ss,ee
    real   (dpR):: dp(3),d,c(3),LHS(3,3),t

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

        t                           =  fv(isec)%turb(1,iele)
        fv(isec)%turb_gra(1:3,iele) =  0.0d0
        do j=sec(isec)%iA_ls(iele),sec(isec)%iA_ls(iele+1)-1
            ss      =  sec(isec)%jA_ls(1,j)
            ee      =  sec(isec)%jA_ls(2,j)
            dp(1:3) =  sec(ss)%cen(1:3,ee)-sec(isec)%cen(1:3,iele)
            d       =  1.0d0/sqrt(dp(1)*dp(1)+dp(2)*dp(2)+dp(3)*dp(3))**ls_weight_order
            c(1:3)  =  d*(LHS(1:3,1)*dp(1)+LHS(1:3,2)*dp(2)+LHS(1:3,3)*dp(3))

            fv(isec)%turb_gra(1:3,iele) =  fv(isec)%turb_gra(1:3,iele) &
                & +(fv(ss)%turb(1,ee)-t)*c(1:3)
        end do
    end do

    return
    end subroutine fv_sa_get_gra_LSm
    end subroutine fv_rans_get_gra
!-------------------------------------------------------------------------------
!   implicit solver for RANS model.
!-------------------------------------------------------------------------------
    subroutine fv_rans_slv(lev)
    use var_kind_def
    use var_fv
    use var_mesh
    use var_prec
    use var_turb
    implicit none
    integer(dpI),intent(in):: lev

    if(transition_model .eq. transition_AGS) then
        call fv_get_intermittency_AGS
    elseif(transition_model .eq. transition_PTM) then
        call fv_get_intermittency_PTM
    elseif(transition_model .eq. transition_BC) then
        call fv_get_intermittency_BC
    end if

    if(RANS_model .eq. RANS_SA) then
        call fv_SA_get_rhs (lev, .true.)
    elseif(RANS_model .eq. RANS_WA) then
        call fv_WA_get_rhs (lev, .true.)
    elseif(RANS_model .eq. RANS_SAM) then
        call fv_SAM_get_rhs(lev, .true.)
    elseif(RANS_model .eq. RANS_SAC) then
        call fv_SAC_get_rhs(lev, .true.)
    elseif(RANS_model .eq. RANS_SST) then
        call fv_SST_get_rhs(lev, .true.)
    else
        stop 'Error: RANS model not supported.'
    end if
    call iteration(lev, fv_rans_prec)
    call fv_rans_bnd_parallel(lev, .true.)

    return
    contains
!       ------------------------------------------------------------------------
!       prec iteration.
!       ------------------------------------------------------------------------
        subroutine iteration(lev,p)
        use var_prec
        implicit none
        integer(dpI),intent(in):: lev
        integer(dpI):: isec,i,iele
        real   (dpR):: nut,t(2),dt(2),rtmp
        type(type_prec):: p

        i   =  1
        do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
            if(.not. sec(isec)%is_int)  cycle

            if(RANS_model .le. RANS_WA) then
                do iele=1,sec(isec)%n_ele
                    p%RHS(i)= -fv(isec)%rhs(1,iele)
                    i       =  i+1
                end do
            else
                do iele=1,sec(isec)%n_ele
                    p%RHS(i  )  = -fv(isec)%rhs(1,iele)
                    p%RHS(i+1)  = -fv(isec)%rhs(2,iele)
                    i           =  i+2
                end do
            end if
        end do

        p%x =  0.0d0
        call prec_solve(p, .true., 2)

        i   =  1
        do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
            if(.not. sec(isec)%is_int)  cycle

            if(RANS_model .le. RANS_WA) then
                do iele=1,sec(isec)%n_ele
                    nut =  fv(isec)%turb(1,iele)
                    rtmp=  min(max(p%x(i), -0.1d0*nut), 0.1d0*nut)
                    fv(isec)%turb(1,iele)   =  max(nut+rtmp, 1.0d-4*nut_ref)
                    i   =  i+1
                end do
            elseif(RANS_model .eq. RANS_SAM) then
                do iele=1,sec(isec)%n_ele
                    dt(1  ) =  p%x(i  )
                    dt(2  ) =  p%x(i+1)/fv(isec)%u(1,iele)
                    t (1:2) =  fv(isec)%turb(1:2,iele) 
                    dt      =  min(max(dt, -0.1d0*t), 0.1d0*t)
                    fv(isec)%turb(1,iele)   =  max(t(1)+dt(1), 1.0d-4*nut_ref)
                    fv(isec)%turb(2,iele)   =  max(t(2)+dt(2), 0.0d0         )
                    i       =  i+2
                end do
            elseif(RANS_model .eq. RANS_SAC) then
                do iele=1,sec(isec)%n_ele
                    dt(1  ) =  p%x(i  )
!                   dt(2  ) =  p%x(i+1)/fv(isec)%u(1,iele)
                    dt(2  ) =  p%x(i+1)
                    t (1:2) =  fv(isec)%turb(1:2,iele) 
                    fv(isec)%turb(1,iele)   =  max(t(1)+dt(1), 1.0d-4*nut_ref)
                    fv(isec)%turb(2,iele)   =      t(2)+dt(2)
                    i       =  i+2
                end do
            elseif(RANS_model .eq. RANS_SST) then
                do iele=1,sec(isec)%n_ele
                    dt(1  ) =  p%x(i  )/fv(isec)%u(1,iele)
                    dt(2  ) =  p%x(i+1)/fv(isec)%u(1,iele)
                    t (1:2) =  fv(isec)%turb(1:2,iele) 
                    fv(isec)%turb(1,iele)   =  max(t(1)+dt(1), 1.0d-4*k_ref)
                    fv(isec)%turb(2,iele)   =  max(t(2)+dt(2), 1.0d-4*o_ref)
                    i       =  i+2
                end do
            end if
        end do

        return
        end subroutine iteration
    end subroutine fv_rans_slv
