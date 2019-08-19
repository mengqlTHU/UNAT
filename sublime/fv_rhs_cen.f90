!-------------------------------------------------------------------------------
!   cal shock sensor and Laplace.
!-------------------------------------------------------------------------------
    subroutine fv_get_ss_su2(lev)
    use var_kind_def
    use var_air, only: gk
    use var_fv
    use var_mesh
    use var_parallel
    implicit none
    integer(dpI),intent(in):: lev
    integer(dpI):: nele,isec,iele,iA(2,1000),i,j,im,s,e,isr,LDA,idx,ip_remote
    real   (dpR):: lap(5),ss,ss_2,s_sp,uL(5),uR(5),aL,aR,uc(5),n(4),u(5)

    iele=  0
    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(sec(isec)%is_bnd)    cycle 
        nele=  sec(isec)%n_ele
        if(.not. allocated(fv(isec)%ss_lap))    allocate(fv(isec)%ss_lap(7,nele))
        fv(isec)%ss_lap =  0.0d0
        if(sec(isec)%is_int) then
            iA(1,isec)  =  iele+1
            iele        =  iele+nele
            iA(2,isec)  =  iele
        end if
    end do

!   ----------------------------------------------------------------------------
!   cal ss and lap for the elements involved in the data exchange.
    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(.not. sec(isec)%is_int)  cycle

        do iele=1,sec(isec)%n_ele_b
            uL(1:5) =  fv(isec)%u (1:5,iele)
            uc(1:5) =  fv(isec)%uc(1:5,iele)
            i       =  iele+iA(1,isec)-1

            lap     =  0.0d0
            ss      =  0.0d0
            ss_2    =  0.0d0
            s_sp    =  0.0d0
            do j=sec(isec)%iA_face_neighbour(iele),sec(isec)%iA_face_neighbour(iele+1)-1
                im  =  sec(isec)%jA_face_neighbour(j)
                if(im .gt. 0) then
                    s   =  mesh(lev)%mortar_LR(3, im)
                    e   =  mesh(lev)%mortar_LR(4, im)
                else
                    s   =  mesh(lev)%mortar_LR(1,-im)
                    e   =  mesh(lev)%mortar_LR(2,-im)
                end if
                im      =  abs(im)
                n (1:4) =  mesh(lev)%mortar_n_vg(1:4,im)
                uR(1:5) =  fv(s)%u(1:5,e)

                lap     =  lap +fv(s)%uc(1:5,e)-uc(1:5)
                ss      =  ss  +(uR(5)-uL(5))
                ss_2    =  ss_2+ uR(5)+uL(5)

                aL  =  sqrt(gk*uL(5)/uL(1))
                aR  =  sqrt(gk*uR(5)/uR(1))
                u(1:3)  =  0.5d0*(uL(2:4)+uR(2:4))-mesh(lev)%mortar_n_vg(6:8,im)
                s_sp    =  s_sp+(abs(u(1)*n(1)+u(2)*n(2)+u(3)*n(3))+0.5d0*(aL+aR))*n(4)
            end do
            fv(isec)%ss_lap(1:5,iele)   =  lap(1:5)
            fv(isec)%ss_lap(6  ,iele)   =  abs(ss)/ss_2
            fv(isec)%ss_lap(7  ,iele)   =  s_sp
        end do
    end do
!   cal ss and lap for the elements involved in the data exchange.
!   ----------------------------------------------------------------------------

    LDA =  7

!   ----------------------------------------------------------------------------
!   prepare data.
    do isr=mesh(lev)%p2p_1,mesh(lev)%p2p_0
        idx =  1
        do i=1,p2p(isr)%n_ele_send
            isec=  p2p(isr)%id_ele_send(1,i)
            iele=  p2p(isr)%id_ele_send(2,i)
            call DCOPY(LDA, fv(isec)%ss_lap(1,iele), 1, p2p(isr)%rsend(idx), 1)
            idx =  idx+LDA
        end do
        p2p(isr)%n_send =  idx-1
        p2p(isr)%n_recv =  LDA*p2p(isr)%n_ele_recv
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
!   cal ss and lap for the elements not involved in the data exchange.
    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(.not. sec(isec)%is_int)  cycle

        do iele=1+sec(isec)%n_ele_b,sec(isec)%n_ele
            uL(1:5) =  fv(isec)%u (1:5,iele)
            uc(1:5) =  fv(isec)%uc(1:5,iele)
            i       =  iele+iA(1,isec)-1

            lap     =  0.0d0
            ss      =  0.0d0
            ss_2    =  0.0d0
            s_sp    =  0.0d0
            do j=sec(isec)%iA_face_neighbour(iele),sec(isec)%iA_face_neighbour(iele+1)-1
                im  =  sec(isec)%jA_face_neighbour(j)
                if(im .gt. 0) then
                    s   =  mesh(lev)%mortar_LR(3, im)
                    e   =  mesh(lev)%mortar_LR(4, im)
                else
                    s   =  mesh(lev)%mortar_LR(1,-im)
                    e   =  mesh(lev)%mortar_LR(2,-im)
                end if
                im      =  abs(im)
                n (1:4) =  mesh(lev)%mortar_n_vg(1:4,im)
                uR(1:5) =  fv(s)%u(1:5,e)

                lap     =  lap +fv(s)%uc(1:5,e)-uc(1:5)
                ss      =  ss  +(uR(5)-uL(5))
                ss_2    =  ss_2+ uR(5)+uL(5)

                aL  =  sqrt(gk*uL(5)/uL(1))
                aR  =  sqrt(gk*uR(5)/uR(1))
                u(1:3)  =  0.5d0*(uL(2:4)+uR(2:4))-mesh(lev)%mortar_n_vg(6:8,im)
                s_sp    =  s_sp+(abs(u(1)*n(1)+u(2)*n(2)+u(3)*n(3))+0.5d0*(aL+aR))*n(4)
            end do
            fv(isec)%ss_lap(1:5,iele)   =  lap(1:5)
            fv(isec)%ss_lap(6  ,iele)   =  abs(ss)/ss_2
            fv(isec)%ss_lap(7  ,iele)   =  s_sp
        end do
    end do
!   cal ss and lap for the elements not involved in the data exchange.
!   ----------------------------------------------------------------------------

    call mpi_waitall(mpi_nreq, mpi_req, mpi_sta, mpi_err)

!   ----------------------------------------------------------------------------
!   unpack data received.
    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(.not. sec(isec)%is_ghost)    cycle
        do i=1,sec(isec)%n_ele
            isr =  sec(isec)%id_recv(1,i)
            iele=  sec(isec)%id_recv(2,i)
!           call DCOPY(LDA, p2p(isr)%rrecv(1+LDA*(iele-1)), 1, fv(isec)%ss_lap(1,i), 1)
            call per_rot_vec(1,sec(isec)%per_path(1,i),p2p(isr)%rrecv(1+LDA*(iele-1)), &
                &  fv(isec)%ss_lap(1,i))
            fv(isec)%ss_lap(6,i)=  p2p(isr)%rrecv(6+LDA*(iele-1))
            fv(isec)%ss_lap(7,i)=  p2p(isr)%rrecv(7+LDA*(iele-1))
        end do
    end do
!   unpack data received.
!   ----------------------------------------------------------------------------

    return
    end subroutine fv_get_ss_su2
!-------------------------------------------------------------------------------
!   cal Rhs, convective part, fv.
!-------------------------------------------------------------------------------
    subroutine fv_get_rhs_jst_su2(lev,calD,addD)
    use var_kind_def
    use var_air, only: gk,gk1
    use var_cgns, only: nface_ele
    use var_fv
    use var_mesh
    use var_slv
    use var_temp, only: rhsl,rhsD
    implicit none
    logical(dpL),intent(in):: calD,addD
    integer(dpI),intent(in):: lev
    logical(dpL):: savD
    integer(dpI):: im,sL,eL,sR,eR,nL,nR
    real   (dpR):: u(5),un,rhun,a,ke,H,e2,e4,d1(5),d3(5),spra,n(5),eig(3),v1(5), &
                &  v2(5),d(5),efix(3),un_a,eiga,eigb,akhL(3),akhR(3),vg(3),vn,unL, &
                &  unR,spL,spR,phi_L,phi_R,phi,s2,s4,tmp1,tmp2

    savD=  calD .and. (.not. addD)
    call fv_get_ss_su2(lev)

    call fv_get_rhs_boundary(lev, calD, addD)
    do im=1+mesh(lev)%n_mortar_b,mesh(lev)%n_mortar
        sL      =  mesh(lev)%mortar_LR(1,im)
        eL      =  mesh(lev)%mortar_LR(2,im)
        sR      =  mesh(lev)%mortar_LR(3,im)
        eR      =  mesh(lev)%mortar_LR(4,im)
        n (1:5) =  mesh(lev)%mortar_n_vg(1:5,im)
        vg(1:3) =  mesh(lev)%mortar_n_vg(6:8,im)

        call u_to_akh(1, fv(sL)%u(1,eL), akhL)
        call u_to_akh(1, fv(sR)%u(1,eR), akhR)
        u(1:5)      =  0.5d0*(fv(sL)%u(1:5,eL)+fv(sR)%u(1:5,eR))
        a           =  sqrt(gk*u(5)/u(1))
        ke          =  0.5d0*(u(2)*u(2)+u(3)*u(3)+u(4)*u(4))
        H           =  a*a/gk1+ke
        un_a        =  u (2)*n(1)+u (3)*n(2)+u (4)*n(3)
        vn          =  vg(1)*n(1)+vg(2)*n(2)+vg(3)*n(3)
        un          =  un_a-vn
        rhun        =  u(1)*un
        rhsl(1  ,1) =  rhun
        rhsl(2:4,1) =  rhun*u(2:4)+n(1:3)*u(5)
        rhsl(5  ,1) =  rhun*0.5d0*(akhL(3)+akhR(3))+u(5)*vn

        d1(1:5) =  fv(sR)%uc    (1:5,eR)-fv(sL)%uc    (1:5,eL)
        d3(1:5) =  fv(sR)%ss_lap(1:5,eR)-fv(sL)%ss_lap(1:5,eL)
        unL     =  fv(sL)%u(2,eL)*n(1)+fv(sL)%u(3,eL)*n(2)+fv(sL)%u(4,eL)*n(3)
        unR     =  fv(sR)%u(2,eR)*n(1)+fv(sR)%u(3,eR)*n(2)+fv(sR)%u(4,eR)*n(3)
        spL     =  abs(unL-vn)+akhL(1)
        spR     =  abs(unR-vn)+akhR(1)
        spra    =  0.5d0*(spL+spR)
        spra    =  abs(un)+a

        if(.true.) then
            phi_L   = (0.25d0*fv(sL)%ss_lap(7,eL)/(spra*n(4)))**0.3d0
            phi_R   = (0.25d0*fv(sR)%ss_lap(7,eR)/(spra*n(4)))**0.3d0
            phi     =  4.0d0*phi_L*phi_R/(phi_L+phi_R)
            nL      =  nface_ele(sec(sL)%ele_type)
            nR      =  nface_ele(sec(sR)%ele_type)
            s2      =  real(3*(nL+nR)/(nL*nR), dpR)
            s4      =  s2*s2*0.25d0
            e2      =  JST_k2*s2*0.5d0*(fv(sL)%ss_lap(6,eL)+fv(sR)%ss_lap(6,eR))
            e4      =  s4*max(0.0d0, JST_k4-e2)

            rhsD(1:5,1) =  spra*phi*(e2*d1-e4*d3)
        else
            e2          =  2.0d0*max(fv(sL)%ss_lap(6,eL), fv(sR)%ss_lap(6,eR))
            e4          =  max(0.0d0, 1.0d0/3.2d1-e2)
            d           =  e2*d1-e4*d3
            eig         = (/un, un+a, un-a/)
            efix        = (/0.2d0, 0.2d0, 0.2d0/)*(a+abs(un))
            eig         =  max(abs(eig), efix)
            eiga        =  0.5d0*(eig(2)+eig(3))-eig(1)
            eigb        =  0.5d0*(eig(2)-eig(3))
            tmp1        = -un_a*d(1)+n(1)*d(2)+n(2)*d(3)+n(3)*d(4)
            tmp2        = (ke  *d(1)-u(2)*d(2)-u(3)*d(3)-u(4)*d(4)+d(5))*gk1
            v1          = (/1.0d0, u(2:4), H   /)
            v2          = (/0.0d0, n(1:3), un_a/)
            rhsD(1:5,1) =  eig(1)*d &
                        & +(eiga*tmp2/a+eigb*tmp1)/a*v1+(eigb*tmp2/a+eiga*tmp1)*v2
        end if
        if(addD)    rhsl(1:5,1) =  rhsl(1:5,1)-rhsD(1:5,1)

        fv(sL)%rhs(1:5,eL)  =  fv(sL)%rhs(1:5,eL)+rhsl(1:5,1)*n(4)
        if(sec(sR)%is_int)  fv(sR)%rhs(1:5,eR)  =  fv(sR)%rhs(1:5,eR)-rhsl(1:5,1)*n(4)
        if(savD) then
            fv(sL)%duc(1:5,eL)  =  fv(sL)%duc(1:5,eL)+rhsD(1:5,1)*n(4)
            if(sec(sR)%is_int)  fv(sR)%duc(1:5,eR)  =  fv(sR)%duc(1:5,eR)-rhsD(1:5,1)*n(4)
        end if
    end do

    return
    end subroutine fv_get_rhs_jst_su2
!-------------------------------------------------------------------------------
!   cal shock sensor and Laplace.
!-------------------------------------------------------------------------------
    subroutine fv_get_ss(lev)
    use var_kind_def
    use var_fv
    use var_mesh
    use var_parallel
    implicit none
    integer(dpI),intent(in):: lev
    integer(dpI):: nele,isec,iele,i,j,im,s,e,isr,LDA,idx,ip_remote
    real   (dpR):: lap(5),lap_2,ss,ss_2,uc(5),uc_R(5),t_R(2),mu_R(2),p,Ad,pR

    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(sec(isec)%is_bnd)    cycle 
        nele=  sec(isec)%n_ele
        if(.not. allocated(fv(isec)%ss_lap))    allocate(fv(isec)%ss_lap(6,nele))
        fv(isec)%ss_lap =  0.0d0
    end do

!   ----------------------------------------------------------------------------
!   cal ss and lap for the elements involved in the data exchange.
    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(.not. sec(isec)%is_int)  cycle

        do iele=1,sec(isec)%n_ele_b
            include './include/get_ss.f90'
        end do
    end do
!   cal ss and lap for the elements involved in the data exchange.
!   ----------------------------------------------------------------------------

    LDA =  6

!   ----------------------------------------------------------------------------
!   prepare data.
    do isr=mesh(lev)%p2p_1,mesh(lev)%p2p_0
        idx =  1
        do i=1,p2p(isr)%n_ele_send
            isec=  p2p(isr)%id_ele_send(1,i)
            iele=  p2p(isr)%id_ele_send(2,i)
            call DCOPY(LDA, fv(isec)%ss_lap(1,iele), 1, p2p(isr)%rsend(idx), 1)
            idx =  idx+LDA
        end do
        p2p(isr)%n_send =  idx-1
        p2p(isr)%n_recv =  LDA*p2p(isr)%n_ele_recv
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
!   cal ss and lap for the elements not involved in the data exchange.
    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(.not. sec(isec)%is_int)  cycle

        do iele=1+sec(isec)%n_ele_b,sec(isec)%n_ele
            include './include/get_ss.f90'
        end do
    end do
!   cal ss and lap for the elements not involved in the data exchange.
!   ----------------------------------------------------------------------------

    call mpi_waitall(mpi_nreq, mpi_req, mpi_sta, mpi_err)

!   ----------------------------------------------------------------------------
!   unpack data received.
    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(.not. sec(isec)%is_ghost)    cycle
        do i=1,sec(isec)%n_ele
            isr =  sec(isec)%id_recv(1,i)
            iele=  sec(isec)%id_recv(2,i)
!           call DCOPY(LDA, p2p(isr)%rrecv(1+LDA*(iele-1)), 1, fv(isec)%ss_lap(1,i), 1)
            call per_rot_vec(1,sec(isec)%per_path(1,i),p2p(isr)%rrecv(1+LDA*(iele-1)), &
                &  fv(isec)%ss_lap(1,i))
            fv(isec)%ss_lap(6,i)=  p2p(isr)%rrecv(6*iele)
        end do
    end do
!   unpack data received.
!   ----------------------------------------------------------------------------

    return
    end subroutine fv_get_ss
!-------------------------------------------------------------------------------
!   cal Rhs, convective part, fv.
!-------------------------------------------------------------------------------
    subroutine fv_get_rhs_jst(lev,calD,addD)
    use var_kind_def
    use var_air, only: gk,gk1
    use var_fv
    use var_mesh
    use var_slv
    use var_temp, only: rhsl,rhsD
    implicit none
    real   (dpR),parameter:: JST_ome=5.0d-2
    logical(dpL),intent(in):: calD,addD
    integer(dpI),intent(in):: lev
    logical(dpL):: savD
    integer(dpI):: im,s_LL,e_LL,sL,eL,sR,eR,s_RR,e_RR
    real   (dpR):: u(5),un,rhun,a,ke,H,e2,e4,d1(5),d3(5),spra,n(5),eig(3),v1(5), &
                &  v2(5),d(5),efix(3),un_a,eiga,eigb,akhL(3),akhR(3),vg(3),vn,unL, &
                &  unR,spL,spR,shock,Amp,Ma,w(5,4),uc(5,4),dL,dR,t(4),mu(8),tmp1,tmp2

    savD=  calD .and. (.not. addD)
    call fv_get_ss(lev)

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
        n (1:5) =  mesh(lev)%mortar_n_vg(1:5,im)
        vg(1:3) =  mesh(lev)%mortar_n_vg(6:8,im)

        w(1:5,1)=  fv(s_LL)%u(1:5,e_LL)
        w(1:5,2)=  fv(sL  )%u(1:5,eL  )
        w(1:5,3)=  fv(sR  )%u(1:5,eR  )
        w(1:5,4)=  fv(s_RR)%u(1:5,e_RR)
        call u_to_uc(4, w, uc, t, mu)

        call u_to_akh(1, w(1,2), akhL)
        call u_to_akh(1, w(1,3), akhR)
        u(1:5)      =  0.5d0*(w(1:5,2)+w(1:5,3))
        un_a        =  u (2)*n(1)+u (3)*n(2)+u (4)*n(3)
        vn          =  vg(1)*n(1)+vg(2)*n(2)+vg(3)*n(3)
        un          =  un_a-vn
        rhun        =  u(1)*un
        rhsl(1  ,1) =  rhun
        rhsl(2:4,1) =  rhun*u(2:4)+n(1:3)*u(5)
        rhsl(5  ,1) =  rhun*0.5d0*(akhL(3)+akhR(3))+u(5)*vn

        unL     =  w(2,2)*n(1)+w(3,2)*n(2)+w(4,2)*n(3)
        unR     =  w(2,3)*n(1)+w(3,3)*n(2)+w(4,3)*n(3)
        spL     =  abs(unL-vn)+akhL(1)
        spR     =  abs(unR-vn)+akhR(1)
        spra    =  0.5d0*(spL+spR)

        if(.true.) then
            dL  =  w(5,1)-w(5,2)
            dR  =  w(5,2)-w(5,3)
            rhun= (1.0d0-JST_ome)*(abs(dL)+abs(dR))+JST_ome*(w(5,1)+2.0d0*w(5,2)+w(5,3))
            spL =  abs(dL-dR)/rhun
            dL  =  w(5,2)-w(5,3)
            dR  =  w(5,3)-w(5,4)
            rhun= (1.0d0-JST_ome)*(abs(dL)+abs(dR))+JST_ome*(w(5,2)+2.0d0*w(5,3)+w(5,4))
            spR =  abs(dL-dR)/rhun
        else
            spL =  abs(w(5,1)-2.0d0*w(5,2)+w(5,3))/(w(5,1)+2.0d0*w(5,2)+w(5,3))
            spR =  abs(w(5,2)-2.0d0*w(5,3)+w(5,4))/(w(5,2)+2.0d0*w(5,3)+w(5,4))
        end if
        shock   =  max(spL, spR)
        if(.true.) then
            u(2:4)  =  u(2:4)-vg(1:3)
            spL     =  u(2)*u(2)+u(3)*u(3)+u(4)*u(4)
            spR     =  gk*u(5)/u(1)
            Ma      =  sqrt(spL/spR)
            if(Ma .le. 0.3d0) then
                Amp =  0.0d0
            elseif(Ma .ge. 0.7d0) then
                Amp =  1.0d0
            else
                Amp =  0.5d0+0.5d0*tanh(2.0d1*(Ma-0.5d0))
            end if
            shock   =  shock*Amp
        end if

        e2      =  0.25d0*shock
        e4      =  max(0.0d0, 1.0d0/3.2d1-e2)
        d       =  e2*(uc(1:5,3)-uc(1:5,2)) &
                & -e4*(uc(1:5,4)-3.0d0*uc(1:5,3)+3.0d0*uc(1:5,2)-uc(1:5,1))
        if(is_shock_sensor) then
            rhsD(1:5,1) =  spra*d*max(fv(sL)%shock_sensor(eL), fv(sR)%shock_sensor(eR))
        else
            rhsD(1:5,1) =  spra*d
        end if

        if(addD)    rhsl(1:5,1) =  rhsl(1:5,1)-rhsD(1:5,1)
        fv(sL)%rhs(1:5,eL)  =  fv(sL)%rhs(1:5,eL)+rhsl(1:5,1)*n(4)
        if(sec(sR)%is_int)  fv(sR)%rhs(1:5,eR)  =  fv(sR)%rhs(1:5,eR)-rhsl(1:5,1)*n(4)
        if(savD) then
            fv(sL)%duc(1:5,eL)  =  fv(sL)%duc(1:5,eL)+rhsD(1:5,1)*n(4)
            if(sec(sR)%is_int)  fv(sR)%duc(1:5,eR)  =  fv(sR)%duc(1:5,eR)-rhsD(1:5,1)*n(4)
        end if
    end do

    do im=1+mesh(lev)%n_mortar_ss,mesh(lev)%n_mortar
        sL      =  mesh(lev)%mortar_LR(1,im)
        eL      =  mesh(lev)%mortar_LR(2,im)
        sR      =  mesh(lev)%mortar_LR(3,im)
        eR      =  mesh(lev)%mortar_LR(4,im)
        n (1:5) =  mesh(lev)%mortar_n_vg(1:5,im)
        vg(1:3) =  mesh(lev)%mortar_n_vg(6:8,im)

        call u_to_akh(1, fv(sL)%u(1,eL), akhL)
        call u_to_akh(1, fv(sR)%u(1,eR), akhR)
        u(1:5)      =  0.5d0*(fv(sL)%u(1:5,eL)+fv(sR)%u(1:5,eR))
        un_a        =  u (2)*n(1)+u (3)*n(2)+u (4)*n(3)
        vn          =  vg(1)*n(1)+vg(2)*n(2)+vg(3)*n(3)
        un          =  un_a-vn
        rhun        =  u(1)*un
        rhsl(1  ,1) =  rhun
        rhsl(2:4,1) =  rhun*u(2:4)+n(1:3)*u(5)
        rhsl(5  ,1) =  rhun*0.5d0*(akhL(3)+akhR(3))+u(5)*vn

        d1(1:5) =  fv(sR)%uc    (1:5,eR)-fv(sL)%uc    (1:5,eL)
        d3(1:5) =  fv(sR)%ss_lap(1:5,eR)-fv(sL)%ss_lap(1:5,eL)
        unL     =  fv(sL)%u(2,eL)*n(1)+fv(sL)%u(3,eL)*n(2)+fv(sL)%u(4,eL)*n(3)
        unR     =  fv(sR)%u(2,eR)*n(1)+fv(sR)%u(3,eR)*n(2)+fv(sR)%u(4,eR)*n(3)
        spL     =  abs(unL-vn)+akhL(1)
        spR     =  abs(unR-vn)+akhR(1)
        spra    =  0.5d0*(spL+spR)
        shock   =  max(fv(sL)%ss_lap(6,eL),fv(sR)%ss_lap(6,eR))

        if(.true.) then
            e2          =  JST_k2*shock
            e4          =  max(0.0d0, JST_k4-e2)
            d           =  e2*d1-e4*d3
            rhsD(1:5,1) =  spra*d
        else
            e2          =  JST_k2*shock
            e4          =  max(0.0d0, JST_k4-e2)
            d           =  e2*d1-e4*d3
            a           =  sqrt(gk*u(5)/u(1))
            ke          =  0.5d0*(u(2)*u(2)+u(3)*u(3)+u(4)*u(4))
            H           =  a*a/gk1+ke
            eig         = (/un, un+a, un-a/)
            efix        = (/0.2d0, 0.2d0, 0.2d0/)*spra
            eig         =  max(abs(eig), efix)
            eiga        =  0.5d0*(eig(2)+eig(3))-eig(1)
            eigb        =  0.5d0*(eig(2)-eig(3))
            tmp1        = -un_a*d(1)+n(1)*d(2)+n(2)*d(3)+n(3)*d(4)
            tmp2        = (ke  *d(1)-u(2)*d(2)-u(3)*d(3)-u(4)*d(4)+d(5))*gk1
            v1          = (/1.0d0, u(2:4), H   /)
            v2          = (/0.0d0, n(1:3), un_a/)
            rhsD(1:5,1) =  eig(1)*d &
                        & +(eiga*tmp2/a+eigb*tmp1)/a*v1+(eigb*tmp2/a+eiga*tmp1)*v2
        end if
        if(is_shock_sensor) rhsD(1:5,1) =  rhsD(1:5,1) &
            & *max(fv(sL)%shock_sensor(eL), fv(sR)%shock_sensor(eR))
        if(addD)    rhsl(1:5,1) =  rhsl(1:5,1)-rhsD(1:5,1)

        fv(sL)%rhs(1:5,eL)  =  fv(sL)%rhs(1:5,eL)+rhsl(1:5,1)*n(4)
        if(sec(sR)%is_int)  fv(sR)%rhs(1:5,eR)  =  fv(sR)%rhs(1:5,eR)-rhsl(1:5,1)*n(4)
        if(savD) then
            fv(sL)%duc(1:5,eL)  =  fv(sL)%duc(1:5,eL)+rhsD(1:5,1)*n(4)
            if(sec(sR)%is_int)  fv(sR)%duc(1:5,eR)  =  fv(sR)%duc(1:5,eR)-rhsD(1:5,1)*n(4)
        end if
    end do

    return
    end subroutine fv_get_rhs_jst
!-------------------------------------------------------------------------------
!   cal Rhs, convective part, fv, KEP.
!-------------------------------------------------------------------------------
    subroutine fv_get_rhs_kep(lev)
    use var_kind_def
    use var_air, only: gk,gk1
    use var_fv
    use var_global, only: R16,R23
    use var_mesh
    implicit none
    integer(dpI),intent(in):: lev
    logical(dpL),parameter:: is_hybrid  =  .true.
    integer(dpI):: cen_order,im,sL,eL,sR,eR,j
    real   (dpR):: uL(5),uR(5),rhun,hL,hR,n(4),r(5),duL(5),duR(5),dL(3),dR(3), &
                &  du(5),vn,pm,p,f_upwind(5),unL,unR,aL,aR,asL,asR,fL(5),fR(5), &
                &  ua,am,MaL,MaR,mL,pL,mR,pR,mf,alp,vg(3),u(5),H

    cen_order   =  3
    call fv_get_rhs_boundary(lev, .true., .true.)
    do im=1+mesh(lev)%n_mortar_b,mesh(lev)%n_mortar
        sL  =  mesh(lev)%mortar_LR(1,im)
        eL  =  mesh(lev)%mortar_LR(2,im)
        sR  =  mesh(lev)%mortar_LR(3,im)
        eR  =  mesh(lev)%mortar_LR(4,im)

        n (1:4) =  mesh(lev)%mortar_n_vg(1:4,im)
        vg(1:3) =  mesh(lev)%mortar_n_vg(6:8,im)
        vn      =  vg(1)*n(1)+vg(2)*n(2)+vg(3)*n(3)
        uL(1:5) =  fv(sL)%u(1:5,eL)
        uR(1:5) =  fv(sR)%u(1:5,eR)

        dL(1:3) =  mesh(lev)%mortar_cen(1:3,im)-sec(sL)%cen(1:3,eL)
        dR(1:3) =  mesh(lev)%mortar_cen(1:3,im)-sec(sR)%cen(1:3,eR)
        if(cen_order .gt. 2) then
            do j=1,5
                duL(j)  =  fv(sL)%gra(3*j-2,eL)*dL(1)+fv(sL)%gra(3*j-1,eL)*dL(2) &
                        & +fv(sL)%gra(3*j  ,eL)*dL(3)
            end do
            do j=1,5
                duR(j)  =  fv(sR)%gra(3*j-2,eR)*dR(1)+fv(sR)%gra(3*j-1,eR)*dR(2) &
                        & +fv(sR)%gra(3*j  ,eR)*dR(3)
            end do
            du  =  uR-uL
            uL  =  uL+R23*duL+R16*du
            uR  =  uR+R23*duR-R16*du
!           uL  =  uL+R23*duL
!           uR  =  uR+R23*duR
        end if

        u       =  0.5d0*(uL+uR)
        H       =  gk*u(5)/(gk1*u(1))+0.5d0*(u(2)*u(2)+u(3)*u(3)+u(4)*u(4))
        pm      =  0.5d0*(uL(5)+uR(5))
        unL     =  uL(2)*n(1)+uL(3)*n(2)+uL(4)*n(3)-vn
        unR     =  uR(2)*n(1)+uR(3)*n(2)+uR(4)*n(3)-vn
        aL      =  sqrt(gk*uL(5)/uL(1))
        aR      =  sqrt(gk*uR(5)/uR(1))
        hL      =  aL*aL/gk1+0.5d0*(uL(2)*uL(2)+uL(3)*uL(3)+uL(4)*uL(4))
        hR      =  aR*aR/gk1+0.5d0*(uR(2)*uR(2)+uR(3)*uR(3)+uR(4)*uR(4))

        if(.false.) then
            rhun= (uL(1)*unL+uR(1)*unR)*0.5d0
        else
            rhun=  sqrt(uL(1)*uR(1))*0.5d0*(unL+unR)
        end if
        r(1  )  =  rhun
        r(2:4)  =  rhun*(uL(2:4)+uR(2:4))*0.5d0+pm*n(1:3)
        if(.false.) then
            r(5)=  rhun*H+pm*vn
        else
            r(5)=  rhun*(0.5d0*(uL(2)*uR(2)+uL(3)*uR(3)+uL(4)*uR(4)) &
                & +aL*aR/(gk*gk1))+0.5d0*(uL(5)*unL+uR(5)*unR)
        end if

        dL  =  dL/sqrt(dL(1)*dL(1)+dL(2)*dL(2)+dL(3)*dL(3))
        dR  = -dR/sqrt(dR(1)*dR(1)+dR(2)*dR(2)+dR(3)*dR(3))
        alp =  max(0.0d0, min(1.0d0, 1.0d0-(dL(1)*dR(1)+dL(2)*dR(2)+dL(3)*dR(3))))
        if((alp .ge. 1.0d-5) .and. is_hybrid) then
!           use AUSM upwind flux instead on the skewed mesh.
            asL     =  sqrt(2.0d0*hL*gk1/(gk+1.0d0))
            asR     =  sqrt(2.0d0*hR*gk1/(gk+1.0d0))
            fL(1:5) = (/1.0d0, uL(2:4), hL/)*uL(1)
            fR(1:5) = (/1.0d0, uR(2:4), hR/)*uR(1)
            ua      = (uL(2)-vg(1))**2+(uL(3)-vg(2))**2+(uL(4)-vg(3))**2 &
                    &+(uR(2)-vg(1))**2+(uR(3)-vg(2))**2+(uR(4)-vg(3))**2 
            ua      =  sqrt(0.5d0*ua)

            am  =  min(asL*asL/max(asL, unL), asR*asR/max(asR, -unR))
            maL =  unL/am
            maR =  unR/am
            if    (maL .ge. 1.0d0) then
                mL  =  maL
                pL  =  1.0d0
            elseif(maL .le. -1.0d0) then
                mL  =  0.0d0
                pL  =  0.0d0
            else
                mL  =  0.25d0*(MaL+1.0d0)**2+0.125d0*(MaL*MaL-1.0d0)**2
!               pL  =  0.25d0*(MaL+1.0d0)**2*(2.0d0-MaL)+0.1875d0*MaL*(MaL**2-1.0d0)**2
                pL  =  0.25d0*(MaL+1.0d0)**2*(2.0d0-MaL)
            end if

            if    (maR .le. -1.0d0) then
                mR  =  maR
                pR  =  1.0d0
            elseif(maR .ge. 1.0d0) then
                mR  =  0.0d0
                pR  =  0.0d0
            else
                mR  = -0.25d0*(MaR-1.0d0)**2-0.125d0*(MaR*MaR-1.0d0)**2
!               pR  =  0.25d0*(MaR-1.0d0)**2*(2.0d0+MaR)-0.1875d0*MaR*(MaR**2-1.0d0)**2
                pR  =  0.25d0*(MaR-1.0d0)**2*(2.0d0+MaR)
            end if
            mf  = (mL+mR)*am
!           p   =  pm+0.5d0*(pL-pR)*(uL(5)-uR(5)) &
!               &+(pL+pR-1.0d0)*ua*0.25d0*(uL(1)+uR(1))*(aL+aR)
            p   =  pL*uL(5)+pR*uR(5)-pL*pR*(uL(1)+uR(1))*ua*(unR-unL)
            mL  =  max(mf, 0.0d0)
            mR  =  min(mf ,0.0d0)
            f_upwind(1)     =  mL*fL(1)  +mR*fR(1)
            f_upwind(2:4)   =  mL*fL(2:4)+mR*fR(2:4)+p*n(1:3)
            f_upwind(5)     =  mL*fL(5)  +mR*fR(5)  +p*vn

            r   = (1.0d0-alp)*r+alp*f_upwind
        end if

        r       =  r*n(4)
        fv(sL)%rhs(1:5,eL)  =  fv(sL)%rhs(1:5,eL)+r(1:5)
        if(sec(sR)%is_int)  fv(sR)%rhs(1:5,eR)  =  fv(sR)%rhs(1:5,eR)-r(1:5)
    end do

    return
    end subroutine fv_get_rhs_kep
