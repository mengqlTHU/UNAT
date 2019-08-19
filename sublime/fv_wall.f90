!-------------------------------------------------------------------------------
!   get BL thickness along the BL line.
!-------------------------------------------------------------------------------
    subroutine fv_get_bl_parameter
    use var_kind_def
    use var_air, only: cp,gk,gk1,rr,t_Suth
    use var_bndv, only: Tu_fs
    use var_cgns, only: BCWallViscous
    use var_eline
    use var_fv
    use var_global
    use var_mesh
    use var_slv, only: is_vis_cal
    use var_turb, only: k2_1,is_RANS
    implicit none
    integer(dpI):: m,i,isec,iele,im,id_edge,id_edge_1,id_edge_2,sR,eR
    real   (dpR):: n(3),u_wall(3),rh(0:1000),u(0:1000),p(0:1000),v(3),rh_edge, &
                &  u_edge,t1(1000),t2(1000),fd(1000),d1,d2,g,ux,uy,uz,vx,vy,vz, &
                &  wx,wy,wz,rd,pt(0:1000),t,tt,max_pt,mu,Re_theta,U_x,U_y,U_z, &
                &  U_s,P_s,Pohl,F,Re_theta_s,Re_theta_e,L_theta,f_Lambda,H12,rtmp, &
                &  s12,s13,s23,Re_v

    if((n_eline .le. 0) .or. (.not. is_vis_cal))    return

    do isec=mesh(0)%sec_1,mesh(0)%sec_0
        if(sec(isec)%bct .ne. BCWallViscous)    cycle
        if(.not. allocated(fv(isec)%Re_theta)) &
            &  allocate(fv(isec)%Re_theta(sec(isec)%n_ele))
    end do

    if(.not. allocated(eline(1)%d)) then
        do m=1,n_eline
            allocate(eline(m)%d(eline(m)%n_ele), stat=err_mem)
            do i=1,eline(m)%n_ele
                isec=  eline(m)%ele(1,i)
                iele=  eline(m)%ele(2,i)
                eline(m)%d(i)   =  sec(isec)%dnw(iele)
            end do
        end do
    end if

    do m=1,n_eline
        im  =  eline(m)%mortar_L
        if(im .le. 0)   cycle

        sR  =  mesh(0)%mortar_LR(3,im)
        eR  =  mesh(0)%mortar_LR(4,im)
        n(1:3)      =  mesh(0)%mortar_n_vg(1:3,im)
        u_wall(1:3) =  fv(sR)%u(2:4,eR)
        p(0)        =  fv(sR)%u(5  ,eR)
        rh(0)       =  fv(sR)%u(1  ,eR)
        u(0)        =  0.0d0
        pt(0)       =  p(0)
        t           =  p(0)/(rr*rh(0))
        mu          =  1.461d-6*sqrt(t**3)/(t+t_suth)

        Re_v=  0.0d0
        do i=1,eline(m)%n_ele
            isec    =  eline(m)%ele(1,i)
            iele    =  eline(m)%ele(2,i)
            rh(i)   =  fv(isec)%u(1,iele)
            v(1:3)  =  fv(isec)%u(2:4,iele)-u_wall(1:3)
            p(i)    =  fv(isec)%u(5,iele)
            rtmp    =  v(1)*n(1)+v(2)*n(2)+v(3)*n(3)
            v       =  v-rtmp*n
            call norm_vec(3, v, u(i))

            t       =  fv(isec)%t(iele)
            tt      =  t+0.5d0*(fv(isec)%u(2,iele)**2+fv(isec)%u(3,iele)**2 &
                    & +fv(isec)%u(4,iele)**2)/cp
            pt(i)   =  p(i)*(t/tt)**(-gk/gk1)

            if(is_RANS) then
                ux  =  fv(isec)%gra(4 ,iele)
                uy  =  fv(isec)%gra(5 ,iele)
                uz  =  fv(isec)%gra(6 ,iele)
                vx  =  fv(isec)%gra(7 ,iele)
                vy  =  fv(isec)%gra(8 ,iele)
                vz  =  fv(isec)%gra(9 ,iele)
                wx  =  fv(isec)%gra(10,iele)
                wy  =  fv(isec)%gra(11,iele)
                wz  =  fv(isec)%gra(12,iele)
                g   =  sqrt(ux*ux+uy*uy+uz*uz+vx*vx+vy*vy+vz*vz+wx*wx+wy*wy+wz*wz)
                g   =  max(g, 1.0d-6*uref/L_ref)
                rd  = (fv(isec)%mu(1,iele)+fv(isec)%mu(2,iele))*k2_1 &
                    &/(g*rh(i)*sec(isec)%dnw(iele)**2)
                fd(i)   =  1.0d0-tanh(5.12d2*rd**3)

                s12 =  0.5d0*(uy+vx)
                s13 =  0.5d0*(uz+wx)
                s23 =  0.5d0*(vz+wy)
                g   =  sqrt(2.0d0*(ux*ux+vy*vy+wz*wz)+4.0d0*(s12*s12+s13*s13+s23*s23))
                rd  =  fv(isec)%u(1,iele)*sec(isec)%dnw(iele)**2*g/fv(isec)%mu(1,iele)
                Re_v=  max(Re_v, rd)
            end if
        end do

!       ------------------------------------------------------------------------
!       the method of total pressure defect to detect the boundary layer edge.
        id_edge_1   = -1
        max_pt      =  p(0)
        do i=1,eline(m)%n_ele
            if(pt(i) .gt. max_pt)   max_pt  =  pt(i)
        end do
        do i=1,eline(m)%n_ele
            if((pt(i)-p(0)) .ge. 0.985d0*(max_pt-p(0))) then
                id_edge_1   =  i
                exit
            end if
        end do
!       the method of total pressure defect to detect the boundary layer edge.
!       ------------------------------------------------------------------------

!       ------------------------------------------------------------------------
!       the method of fd.
        id_edge_2   = -1
        if(is_RANS) then
            do i=1,eline(m)%n_ele
                if(fd(i) .ge. 0.999d0) then
                    id_edge_2   =  i
                    exit
                end if
            end do
        end if
!       the method of fd.
!       ------------------------------------------------------------------------

        if((id_edge_1 .gt. 0) .and. (id_edge_2 .gt. 0)) then
            id_edge =  max(id_edge_1, id_edge_2)
        elseif(id_edge_2 .gt. 0) then
            id_edge =  id_edge_2
        elseif(id_edge_1 .gt. 0) then
            id_edge =  id_edge_1
        else
            cycle
        end if

        isec    =  eline(m)%ele(1,id_edge)
        iele    =  eline(m)%ele(2,id_edge)
        rh_edge =  rh(id_edge)
        u_edge  =  u (id_edge)
        v(1:3)  =  fv(isec)%u(2:4,iele)-u_wall(1:3)
        rtmp    =  v(1)*n(1)+v(2)*n(2)+v(3)*n(3)
        v       =  v-rtmp*n
        call norm_vec(3, v, rtmp)

        do i=1,id_edge
            t1(i)   =  1.0d0-rh(i)*u(i)/(rh_edge*u_edge)
            t2(i)   = (1.0d0-u(i)/u_edge)*(rh(i)*u(i))/(rh_edge*u_edge)
        end do
        d1  =  0.5d0*eline(m)%d(1)*(1.0d0+t1(1))
        d2  =  0.5d0*eline(m)%d(1)*(0.0d0+t2(1))
        do i=1,id_edge-1
            d1  =  0.5d0*(eline(m)%d(i+1)-eline(m)%d(i))*(t1(i+1)+t1(i))+d1
            d2  =  0.5d0*(eline(m)%d(i+1)-eline(m)%d(i))*(t2(i+1)+t2(i))+d2
        end do
        H12     =  max(1.01d0, d1/d2)
        Re_theta=  max(2.0d1, rh_edge*u_edge*d2/mu)
!       print*,m,Re_v/2.193,Re_theta
!       print*,eline(m)%n_ele,id_edge_1,id_edge_2
!       write(unit=6,fmt='(I4,4ES20.12)'),m,eline(m)%d(id_edge),d1,d2,H12

        ux  =  fv(isec)%gra(4 ,iele)
        uy  =  fv(isec)%gra(5 ,iele)
        uz  =  fv(isec)%gra(6 ,iele)
        vx  =  fv(isec)%gra(7 ,iele)
        vy  =  fv(isec)%gra(8 ,iele)
        vz  =  fv(isec)%gra(9 ,iele)
        wx  =  fv(isec)%gra(10,iele)
        wy  =  fv(isec)%gra(11,iele)
        wz  =  fv(isec)%gra(12,iele)

        if(.true.) then
            U_x = (v(1)*ux+v(2)*vx+v(3)*wx)/sqrt(u_edge)
            U_y = (v(1)*uy+v(2)*vy+v(3)*wy)/sqrt(u_edge)
            U_z = (v(1)*uz+v(2)*vz+v(3)*wz)/sqrt(u_edge)
            U_s = (v(1)*U_x+v(2)*U_y+v(3)*U_z)/u_edge
            Pohl=  d2**2*U_s/(mu/rh(0))
        else
            P_s = (fv(isec)%gra(13,iele)*v(1)+fv(isec)%gra(14,iele)*v(2) &
                & +fv(isec)%gra(15,iele)*v(3))/u_edge
            Pohl= -d2**2/(mu*u_edge)*P_s
        end if
        Pohl=  min(max(-0.1d0, Pohl), 0.1d0)

        if(.true.) then
            L_theta     =  min(Pohl, 5.8d-2*(H12-4.0d0)**2/(H12-1.0d0)-6.8d-2)
            f_Lambda    =  2.2d2*exp(-Tu_fs/9.0d-3)*(atan(L_theta/2.0d-2+0.842d0)-0.7d0)
            Re_theta_s  =  1.045d3*exp(-Tu_fs/1.0d-2)+f_Lambda+1.55d2
            Re_theta_e  =  Re_theta_s*(1.6d0+1.3d0*exp(-Tu_fs/2.0d-2))
            if(Re_theta .ge. Re_theta_e) then
                eline(m)%intermittency  =  1.0d0
            elseif(Re_theta .le. Re_theta_s) then
                eline(m)%intermittency  =  0.0d0
            else
                rtmp= (Re_theta-Re_theta_s)/(Re_theta_e-Re_theta_s)
                eline(m)%intermittency  =  1.0d0-exp(-5.0d0*rtmp**1.2d0)
            end if
        else
            if(Pohl .le. 0.0d0) then
                F   =  6.91d0+12.75d0*Pohl+63.34d0*Pohl**2
            else
                F   =  6.91d0+2.48d0*Pohl-12.27d0*Pohl**2
            end if
            Re_theta_s  =  1.63d2+exp(F-F/6.91d0*Tu_fs*1.0d2)
            if(Re_theta .ge. Re_theta_s) then
!               print*,m,':',Re_theta,'>',Re_theta_s
                eline(m)%intermittency  =  1.0d0
            else
!               print*,m,':',Re_theta,'<',Re_theta_s
                eline(m)%intermittency  =  0.0d0
            end if
        end if

        eline(m)%BL_thickness   =  eline(m)%d(id_edge)

        isec=  mesh(0)%mortar_LR(3,im)
        iele=  mesh(0)%mortar_LR(4,im)
        fv(isec)%Re_theta(iele) =  Re_theta
!       fv(isec)%Re_theta(iele) =  eline(m)%intermittency
    end do

!   call eline_output

    return
    end subroutine fv_get_bl_parameter
!-------------------------------------------------------------------------------
!   get distance to the nearest wall.
!-------------------------------------------------------------------------------
    subroutine fv_get_dnw
    use var_kind_def
    use var_fv
    use var_global, only: is_2d_cal,err_mem,BL_thick
    use var_mesh
    use var_rtree
    use var_slv, only: is_vis_cal
    use var_turb, only: transition_model,transition_Coder,transition_AGS,is_DDES
    use var_wall, iA=>iA_wall,jA=>jA_wall,xyz=>xyz_wall
    implicit none
    logical(dpL):: is_inner
    integer(dpI):: i,j,isec,iele,n_ele,idx,ip,ivtx,nnp,list(1000),ivtx_sol
    real   (dpR):: d(3),r1(3,3),p(3,4),n(3),un,c(3,1000),v(8),rtmp
    integer(dpI),allocatable:: nearest_vtx(:,:),iA_e2n(:),jA_e2n(:)
    real   (dpR),allocatable:: dnw(:)

    if(.not. is_vis_cal)    return

    call mesh_get_wall

!   ----------------------------------------------------------------------------
!   first pass: get dnw for vertex.
    allocate(nearest_vtx(max_nearest_point, mesh(0)%n_vtx), stat=err_mem)
    allocate(dnw        (                   mesh(0)%n_vtx), stat=err_mem)
    nearest_vtx =  0
    do i=1,mesh(0)%n_vtx
        if(is_2d_cal) then
            d(1:2)  =  mesh(0)%xyz(1:2,i)
            d(3  )  =  0.0d0
        else
            d(1:3)  =  mesh(0)%xyz(1:3,i)
        end if
        call rtree_search(tree, d, huge(1.0d0), nnp, nearest_vtx(1,i), dnw(i))
    end do
!   first pass: get dnw for vertex.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   get the wall-normal direction for Coder transition model.
    if(transition_model .eq. transition_Coder)  call fv_get_wn(dnw)
!   get the wall-normal direction for Coder transition model.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   second pass: get dnw for element.
    allocate(iA_e2n(n_vtx_wall+1), stat=err_mem)
    allocate(jA_e2n(size(jA)    ), stat=err_mem)
    call itrans_spy(.false., 1, n_ele_wall, iA, jA, n_vtx_wall, iA_e2n, jA_e2n)

    do isec=mesh(0)%sec_1,mesh(0)%sec_0
        if((.not. allocated(sec(isec)%dnw)) .and. (.not. sec(isec)%is_bnd)) &
            &  allocate(sec(isec)%dnw(sec(isec)%n_ele), stat=err_mem)

        if(.not. sec(isec)%is_int)  cycle

        if(transition_model .eq. transition_AGS) then
            allocate(fv(isec)%nearest_wall(sec(isec)%n_ele), stat=err_mem)
            fv(isec)%nearest_wall   =  0
        end if
        do iele=1,sec(isec)%n_ele
            rtmp    =  0.0d0
            n_ele   =  0

!           step 1: averaged dnw.
            do i=1,sec(isec)%npe
                ivtx=  sec(isec)%n2e(i,iele)
                rtmp=  rtmp+dnw(ivtx)

                do j=1,max_nearest_point
                    ivtx_sol=  nearest_vtx(j,ivtx)
                    if(ivtx_sol .le. 0) cycle
                    do ip=iA_e2n(ivtx_sol),iA_e2n(ivtx_sol+1)-1
                        n_ele       =  n_ele+1
                        list(n_ele) =  jA_e2n(ip)
                    end do
                end do
            end do
            rtmp=  rtmp/real(sec(isec)%npe, dpR)

!           step 2: correction.
            call simplify_series(n_ele, 1, 1, list)
            d(1:3        )  =  sec(isec)%cen(1:3,iele)
            c(1:3,1:n_ele)  =  0.0d0
            do i=1,n_ele
                idx =  list(i)
                do j=iA(idx),iA(idx+1)-1
                    p(1:3,j-iA(idx)+1)  =  xyz(1:3,jA(j))
                    c(1:3,j-iA(idx)+1)  =  c(1:3,j-iA(idx)+1)+xyz(1:3,jA(j))
                end do
                c(1:3,j-iA(idx)+1)  =  c(1:3,j-iA(idx)+1)/real(iA(idx+1)-iA(idx), dpR)

                if(iA(idx+1)-iA(idx) .eq. 2) then
                    call get_dis_to_bar(d, p(1,1), p(1,2), is_inner, r1(1,1))
                    if(is_inner) then
                        n(1)= -p(2,2)+p(2,1)
                        n(2)=  p(1,2)-p(1,1)
                        n(3)=  0.0d0
                        call norm_vec(3, n, un)
                    end if
                elseif(iA(idx+1)-iA(idx) .eq. 3) then
                    call get_dis_to_tri(d, p(1,1), p(1,2), p(1,3), is_inner, r1(1,1), n)
                elseif(iA(idx+1)-iA(idx) .eq. 4) then
                    call get_dis_to_quad(d,p(1,1),p(1,2),p(1,3),p(1,4),is_inner,r1(1,1),n)
                else
                    stop 'Error: wrong wall element type.'
                end if

                if(is_inner) then
                    rtmp=  min(rtmp, r1(1,1))
                    un  =  fv(isec)%u(2,iele)*n(1)+fv(isec)%u(3,iele)*n(2) &
                        & +fv(isec)%u(4,iele)*n(3)
                    if(rtmp .le. BL_thick) &
                        &  fv(isec)%u(2:4,iele)=  fv(isec)%u(2:4,iele)-un*n(1:3)
                end if
            end do
            sec(isec)%dnw(iele) =  rtmp

            if(transition_model .eq. transition_AGS) then
                un  =  huge(1.0d0)
                do i=1,n_ele
                    rtmp= (d(1)-c(1,i))**2+(d(2)-c(2,i))**2+(d(3)-c(3,i))**2
                    if(rtmp .lt. un) then
                        un  =  rtmp
                        fv(isec)%nearest_wall(iele) =  list(i)
                    end if
                end do
            end if
        end do
    end do
!   second pass: get dnw for element.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   get wall-normal grid leng scale for IDDES.
    if(is_DDES) then
        do isec=mesh(0)%sec_1,mesh(0)%sec_0
            if(.not. sec(isec)%is_int)  cycle
            allocate(sec(isec)%gls_wn(sec(isec)%n_ele), stat=err_mem)
            do iele=1,sec(isec)%n_ele
                do ivtx=1,sec(isec)%npe
                    v(ivtx) =  dnw(sec(isec)%n2e(ivtx,iele))
                end do
                sec(isec)%gls_wn(iele)  =  maxval(v(1:sec(isec)%npe)) &
                                        & -minval(v(1:sec(isec)%npe))
            end do
        end do
    end if
!   get wall-normal grid leng scale for IDDES.
!   ----------------------------------------------------------------------------

    call rtree_destroy(tree)
    if(allocated(iA         ))  deallocate(iA)
    if(allocated(jA         ))  deallocate(jA)
    if(allocated(nearest_vtx))  deallocate(nearest_vtx)
    if(allocated(iA_e2n     ))  deallocate(iA_e2n)
    if(allocated(jA_e2n     ))  deallocate(jA_e2n)
    if(allocated(xyz        ))  deallocate(xyz)
    if(allocated(dnw        ))  deallocate(dnw)

!   call fv_wr_cc_L(1)

    call fv_exchange_dnw(0)

    return
    end subroutine fv_get_dnw
!-------------------------------------------------------------------------------
!   exchange dnw.
!-------------------------------------------------------------------------------
    subroutine fv_exchange_dnw(lev)
    use var_kind_def
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
            p2p(isr)%rsend(idx) =  sec(isec)%dnw(iele)
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
            sec(isec)%dnw(i)=  p2p(isr)%rrecv(iele)
        end do
    end do
!   unpack data received.
!   ----------------------------------------------------------------------------

    return
    end subroutine fv_exchange_dnw
!-------------------------------------------------------------------------------
!   get wall normal direction.
!-------------------------------------------------------------------------------
    subroutine fv_get_wn(dnw)
    use var_kind_def
    use var_cgns
    use var_fv
    use var_global, only: err_mem
    use var_mesh
    use var_parallel
    implicit none
    real   (dpR),intent(in):: dnw(*)
    integer(dpI):: isec,iele,i,npe,im,sL,eL,sR,eR,LDA,isr,ip_remote,idx
    real   (dpR):: n(3),d

!   ----------------------------------------------------------------------------
!   memory allocation.
    do isec=mesh(0)%sec_1,mesh(0)%sec_0
        if(.not. sec(isec)%is_bnd) then
            allocate(sec(isec)%wn(3,sec(isec)%n_ele), stat=err_mem)
            sec(isec)%wn=  0.0d0
        end if
    end do
!   memory allocation.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   surface integration.
    do im=1,mesh(0)%n_mortar
        sL  =  mesh(0)%mortar_LR(1,im)
        eL  =  mesh(0)%mortar_LR(2,im)
        sR  =  mesh(0)%mortar_LR(3,im)
        eR  =  mesh(0)%mortar_LR(4,im)

        npe =  npe_ele(mesh(0)%mortar_ele_type(im))
        d   =  0.0d0
        do i=1,npe
            d   =  d+dnw(mesh(0)%mortar_n2e(i,im))
        end do
        n(1:3)  =  mesh(0)%mortar_n_vg(1:3,im)*mesh(0)%mortar_n_vg(4,im) &
                & *d/real(npe, dpR)
        sec(sL)%wn(1:3,eL)  =  sec(sL)%wn(1:3,eL)+n(1:3)
        if(sec(sR)%is_int)  sec(sR)%wn(1:3,eR)  =  sec(sR)%wn(1:3,eR)-n(1:3)
    end do
!   surface integration.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   get the wall normal direction.
    do isec=mesh(0)%sec_1,mesh(0)%sec_0
        if(.not. sec(isec)%is_int)  cycle
        do iele=1,sec(isec)%n_ele
            n(1:3)  =  sec(isec)%wn(1:3,iele)/sec(isec)%vol(iele)
            sec(isec)%wn(1:3,iele)  =  n(1:3)/sqrt(n(1)**2+n(2)**2+n(3)**2)
        end do
    end do
!   get the wall normal direction.
!   ----------------------------------------------------------------------------

    LDA =  3

!   ----------------------------------------------------------------------------
!   prepare data.
    do isr=mesh(0)%p2p_1,mesh(0)%p2p_0
        idx =  1
        do i=1,p2p(isr)%n_ele_send_face
            isec=  p2p(isr)%id_ele_send(1,i)
            iele=  p2p(isr)%id_ele_send(2,i)
            call DCOPY(LDA, sec(isec)%wn(1,iele), 1, p2p(isr)%rsend(idx), 1)
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
    do isr=mesh(0)%p2p_1,mesh(0)%p2p_0
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
    do isr=mesh(0)%p2p_1,mesh(0)%p2p_0
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
    do isec=mesh(0)%sec_1,mesh(0)%sec_0
        if(.not. sec(isec)%is_ghost)    cycle
        do i=1,sec(isec)%n_ele
            isr =  sec(isec)%id_recv(1,i)
            iele=  sec(isec)%id_recv(2,i)
            if(iele .gt. p2p(isr)%n_ele_recv_face)  cycle
            call per_rot_vec(12,sec(isec)%per_path(1,i),p2p(isr)%rrecv(1+LDA*(iele-1)), &
                &  sec(isec)%wn(1,i))
        end do
    end do
!   unpack data received.
!   ----------------------------------------------------------------------------

    return
    end subroutine fv_get_wn
!-------------------------------------------------------------------------------
!   get y+.
!-------------------------------------------------------------------------------
    subroutine fv_get_yp(lev)
    use var_kind_def
    use var_cgns, only: BCWallViscous
    use var_fv
    use var_global, only: is_show_res
    use var_mesh
    use var_parallel
    use var_slv, only: is_vis_cal
    implicit none
    integer(dpI),intent(in):: lev
    integer(dpI):: im,sL,eL,sR,eR
    real   (dpR):: yp_1,yp_0,uL(5),u_n,tw,u_f,yp

    if(.not. is_vis_cal)    return

    yp_1=  1.0d99
    yp_0= -1.0d99
    do im=1,mesh(lev)%n_mortar_b
        sL      =  mesh(lev)%mortar_LR(1,im)
        eL      =  mesh(lev)%mortar_LR(2,im)
        sR      =  mesh(lev)%mortar_LR(3,im)
        eR      =  mesh(lev)%mortar_LR(4,im)
        if(sec(sR)%bct .ne. BCWallViscous)  cycle
        uL(1:5) =  fv(sR)%uL(1:5,eR)

        u_n = sqrt((fv(sL)%u(2,eL)-mesh(lev)%mortar_n_vg(6,im))**2 &
                & +(fv(sL)%u(3,eL)-mesh(lev)%mortar_n_vg(7,im))**2 &
                & +(fv(sL)%u(4,eL)-mesh(lev)%mortar_n_vg(8,im))**2)/sec(sL)%dnw(eL)
        tw  =  fv(sL)%mu(1,eL)*u_n
        u_f =  sqrt(tw/uL(1))
        yp  =  u_f*sec(sL)%dnw(eL)*uL(1)/fv(sL)%mu(1,eL)
        yp_1=  min(yp, yp_1)
        yp_0=  max(yp, yp_0)
    end do
    call mpi_allreduce(yp_1, min_yp, 1, mpi_dpR, mpi_min, mpi_comm_world, mpi_err)
    call mpi_allreduce(yp_0, max_yp, 1, mpi_dpR, mpi_max, mpi_comm_world, mpi_err)
    if((myid .eq. 0) .and. is_show_res) then
        print*,'----------------------------------------------------------'
        write(unit=6,fmt='(A8,es12.4,A11,es12.4)'),'Min_y+=',min_yp,', max_y+=',max_yp
    end if

    return
    end subroutine fv_get_yp
