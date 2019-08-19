!-------------------------------------------------------------------------------
!   output data for DMD, cc.
!-------------------------------------------------------------------------------
    subroutine fv_wr_slice_cc
    use var_kind_def
    use var_fv
    use var_global, only: err_mem
    use var_load_balance
    use var_mesh
    use var_parallel
    use var_uns_cal, only: uns_iter,uns_ra_begin,uns_write_delta
    implicit none
    character(len=80):: str
    integer(dpI):: LDA,nele,iseg,e1,e0,M,isec,iele,ID,ip,ib

    if((uns_iter .lt. uns_ra_begin) .or. (uns_write_delta .gt. 100000)) return
    if(mod(uns_iter-uns_ra_begin, uns_write_delta) .ne. 0)  return

    if(myid .eq. 0) then
        write(str,*),uns_iter
        str =  './data/dmd/'//trim(adjustl(str))//'.dat'
        open(unit=10,file=trim(adjustl(str)))
    end if

    LDA =  6
    nele=  min(2**20, n_ele_g)

    if(.not. allocated(r_send)) then
        allocate(r_send(LDA*nele), stat=err_mem)
    else
        if(size(r_send) .lt. nele*LDA) then
            deallocate(r_send)
            allocate(r_send(LDA*nele), stat=err_mem)
        end if
    end if
    if(.not. allocated(r_recv)) then
        allocate(r_recv(LDA*nele), stat=err_mem)
    else
        if(size(r_recv) .lt. nele*LDA) then
            deallocate(r_recv)
            allocate(r_recv(LDA*nele), stat=err_mem)
        end if
    end if

    do iseg=1,(n_elei_g-1)/nele+1
        e1  =  1+nele*(iseg-1)
        e0  =  min(iseg*nele, n_elei_g)

!       every processor checks the number of elements to be output in this partial write.
        m   =  0
        do isec=mesh(0)%sec_1,mesh(0)%sec_0
            if(.not. sec(isec)%is_int)  cycle
            do iele=1,sec(isec)%n_ele
                ID  =  sec(isec)%ID_ele_i(iele)
                if((ID .lt. e1) .or. (ID .gt. e0))  cycle
                m   =  m+1
                r_send(1+LDA*(m-1)            ) =  real(ID, dpR)
                r_send(2+LDA*(m-1):6+LDA*(m-1)) =  fv(isec)%u(1:5,iele)
            end do
        end do

        call mpi_gather(m, 1, mpi_dpI, prc_info, 1, mpi_dpI, 0, mpi_comm_world, mpi_err)

        if((myid .ne. 0) .and. (m .gt. 0))  &
            &  call mpi_send(r_send, LDA*m, mpi_dpR, 0, myid, mpi_comm_world, mpi_err)

        if(myid .ne. 0) cycle
        m   =  1
        do ip=0,nprc-1
            if(prc_info(ip,1) .le. 0)   cycle

            ib  =  prc_info(ip,1)*LDA
            if(ip .eq. 0) then
                call DCOPY(ib, r_send, 1, r_recv(m), 1)
            else
                call mpi_recv(r_recv(m), ib, mpi_dpR, ip, ip, mpi_comm_world, &
                    &  mpi_status, mpi_err)
            end if
            m   =  m+ib
        end do
        if(m-1 .ne. LDA*(e0-e1+1))  stop 'Error: fails to gather solution, cc.'
        call dqsortcols(.true., 1, e0-e1+1, 1, LDA, r_recv)

        do iele=1,e0-e1+1
            ID  =  LDA*(iele-1)
            write(unit=10,fmt='(5ES20.12)'),r_recv(ID+2:ID+6)
        end do
    end do

    if(myid .eq. 0) close(10)

    return
    end subroutine fv_wr_slice_cc
!-------------------------------------------------------------------------------
!   output data for DMD, cc.
!-------------------------------------------------------------------------------
    subroutine fv_wr_slice_vc
    use var_kind_def
    use var_fv
    use var_global
    use var_load_balance, only: n_vtx_g
    use var_mesh
    use var_parallel
    use var_uns_cal, only: uns_iter,uns_ra_begin,uns_write_delta
    implicit none
    character(len=80):: str
    integer(dpI):: nvtx,ivtx,ID,v1,v0,iseg,i,j,m,LDA,ip
    real   (dpR):: v(10)

    if((uns_iter .lt. uns_ra_begin) .or. (uns_write_delta .gt. 100000)) return
    if(mod(uns_iter-uns_ra_begin, uns_write_delta) .ne. 0)  return

!   ----------------------------------------------------------------------------
!   every processor prepares the vertex based data to be output.
    call fv_get_vc(1)
!   every processor prepares the vertex based data to be output.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   how many vertex/elements will be output in one call of partial_write.
    nvtx=  min(2**20, n_vtx_g)
!   how many vertex/elements will be output in one call of partial_write.
!   ----------------------------------------------------------------------------

    LDA =  6
    if(.not. allocated(r_send)) then
        allocate(r_send(LDA*nvtx) , stat=err_mem)
    else
        if(size(r_send) .lt. LDA*nvtx) then
            deallocate(r_send)
            allocate(r_send(LDA*nvtx) , stat=err_mem)
        end if
    end if
    if(myid .eq. 0) then
        if(.not. allocated(r_recv)) then
            allocate(r_recv(LDA*nvtx) , stat=err_mem)
        else
            deallocate(r_recv)
            allocate(r_recv(LDA*nvtx) , stat=err_mem)
        end if
        if(.not. allocated(r_recv_all)) then
            allocate(r_recv_all(nvtx*(LDA-1)), stat=err_mem)
        else
            deallocate(r_recv_all)
            allocate(r_recv_all(nvtx*(LDA-1)), stat=err_mem)
        end if
    end if

    if(myid .eq. 0) then
        write(str,*),uns_iter
        str =  './data/dmd/'//trim(adjustl(str))//'.dat'
        open(unit=10,file=trim(adjustl(str)))
    end if

!   ----------------------------------------------------------------------------
!   output vertex based information.
    do iseg=1,(n_vtx_g-1)/nvtx+1
        v1  =  1+nvtx*(iseg-1)
        v0  =  min(nvtx*iseg, n_vtx_g)

        m   =  0
        do ivtx=1,mesh(0)%n_vtx
            ID  =  mesh(0)%id_vtx(ivtx)
            if((ID .lt. v1) .or. (ID .gt. v0))  cycle
            m   =  m+1
            r_send(1+LDA*(m-1)            ) =  real(ID, dpR)
            if(is_2d_cal) then
                r_send(2+LDA*(m-1)) =  fv_vertex_variables(1,ivtx)
                r_send(3+LDA*(m-1)) =  fv_vertex_variables(2,ivtx)
                r_send(4+LDA*(m-1)) =  fv_vertex_variables(3,ivtx)
                r_send(5+LDA*(m-1)) =  0.0d0
                r_send(6+LDA*(m-1)) =  fv_vertex_variables(4,ivtx)
            else
                r_send(2+LDA*(m-1):6+LDA*(m-1)) =  fv_vertex_variables(1:5,ivtx)
            end if
        end do

        call mpi_gather(m, 1, mpi_dpI, prc_info, 1, mpi_dpI, 0, mpi_comm_world, mpi_err)

        if((myid .ne. 0) .and. (m .gt. 0))  &
            &  call mpi_send(r_send, LDA*m, mpi_dpR, 0, myid, mpi_comm_world, mpi_err)

        if(myid .ne. 0) cycle
        do ip=0,nprc-1
            if(prc_info(ip,1) .le. 0)   cycle
            if(ip .eq. 0) then
                call DCOPY(LDA*prc_info(ip,1), r_send, 1, r_recv, 1)
            else
                call mpi_recv(r_recv, LDA*prc_info(ip,1), mpi_dpR, ip, ip, &
                    &  mpi_comm_world, mpi_status, mpi_err)
            end if
!           decode the received information.
            do i=1,prc_info(ip,1)
                ID  =  nint(r_recv(1+LDA*(i-1)), dpI)
                ivtx=  ID-v1+1
                forall(j=1:LDA-1)   r_recv_all(ivtx+nvtx*(j-1)) =  r_recv(j+1+LDA*(i-1))
            end do
        end do

        do i=v1,v0
            forall(j=1:5)   v(j)=  r_recv_all(i+nvtx*(j-1))
            write(unit=10,fmt='(5ES20.12)'),v(1:5)
        end do
    end do
!   output vertex based information.
!   ----------------------------------------------------------------------------

    if(myid .eq. 0) close(10)

    return
    end subroutine fv_wr_slice_vc
!-------------------------------------------------------------------------------
!   output unsteady records.
!-------------------------------------------------------------------------------
    subroutine fv_wr_uns_record
    use var_kind_def
    use var_cgns, only: BCWallViscous
    use var_fv
    use var_global, only: mesh_name
    use var_mesh
    use var_parallel
    use var_slv
    use var_uns_cal
    implicit none
    character(len=80):: str(10)
    character(len=10000):: ot,tec
    integer(dpI):: i,j,ip,N,isec,iele,c
    real   (dpR):: r(10000)

    n_uns_record=  n_uns_record+1
    if(myid .eq. 0) then
        uns_record(n_uns_record,1 ) =  real(uns_iter, dpR)
        uns_record(n_uns_record,2 ) =  tph
        uns_record(n_uns_record,3 ) =  lift
        uns_record(n_uns_record,4 ) =  drag
        uns_record(n_uns_record,5 ) =  lift_p
        uns_record(n_uns_record,6 ) =  lift_v
        uns_record(n_uns_record,7 ) =  drag_p
        uns_record(n_uns_record,8 ) =  drag_v
        uns_record(n_uns_record,9 ) =  mfin
        uns_record(n_uns_record,10) =  mfot
        j   =  10
    else
        j   =  0
    end if

    do i=1,n_monitor_L
        isec=  monitor_idx_L(1,i)
        iele=  monitor_idx_L(2,i)
        if((sec(isec)%bct .eq. BCWallViscous) .and. is_vis_cal) then
            uns_record(n_uns_record,5*i-4+j)=  fv(isec)%u(1,iele)
            uns_record(n_uns_record,5*i-3+j)=  fv(isec)%bnd_solution(1,iele)
            uns_record(n_uns_record,5*i-2+j)=  fv(isec)%bnd_solution(2,iele)
            uns_record(n_uns_record,5*i-1+j)=  0.0d0
            uns_record(n_uns_record,5*i  +j)=  fv(isec)%u(5,iele)
        else
            uns_record(n_uns_record,5*i-4+j)=  fv(isec)%u(1,iele)
            uns_record(n_uns_record,5*i-3+j)=  fv(isec)%u(2,iele)
            uns_record(n_uns_record,5*i-2+j)=  fv(isec)%u(3,iele)
            uns_record(n_uns_record,5*i-1+j)=  fv(isec)%u(4,iele)
            uns_record(n_uns_record,5*i  +j)=  fv(isec)%u(5,iele)
        end if
    end do

    ot  =  trim(adjustl(mesh_name))//'_uns_data.dat'
    if((uns_iter .eq. 1) .and. (myid .eq. 0)) then
        tec     =  'variables='
        str(1 ) =  '"iteration",'
        str(2 ) =  '"time",'
        str(3 ) =  '"lift",'
        str(4 ) =  '"drag",'
        str(5 ) =  '"lift_p",'
        str(6 ) =  '"lift_v",'
        str(7 ) =  '"drag_p",'
        str(8 ) =  '"drag_v",'
        str(9 ) =  '"mfin",'
        str(10) =  '"mfot"'

        do i=1,10
            tec =  trim(adjustl(tec))//trim(adjustl(str(i)))
        end do

        open(unit=10,file=trim(adjustl(ot)))
        write(unit=10,fmt='(A21,I8)'),'#Number_of_monitors= ',n_monitor
        write(unit=10,fmt='(A21,I8)'),'#Number_of_variables=',10+5*n_monitor
        do i=1,n_monitor
            write(unit=10,fmt='(A9,I4,A1,2I8)'),'#Monitor ',i,'=',monitor_idx(1:2,i)
            write(str(1),*),i
            tec =  trim(adjustl(tec))//',"rh_'//trim(adjustl(str(1)))//'"'
            tec =  trim(adjustl(tec))//',"u_' //trim(adjustl(str(1)))//'"'
            tec =  trim(adjustl(tec))//',"v_' //trim(adjustl(str(1)))//'"'
            tec =  trim(adjustl(tec))//',"w_' //trim(adjustl(str(1)))//'"'
            tec =  trim(adjustl(tec))//',"p_' //trim(adjustl(str(1)))//'"'
        end do
        write(unit=10,fmt='(A)'),trim(tec)
        close(10)
    end if

    if((n_uns_record .eq. max_uns_record) .or. (uns_iter .ge. max_uns_iter)) then
        c   =  1
        do ip=0,nprc-1
            if(ip .eq. 0) then
                N   =  max_uns_record*(5*n_monitor_prc(ip)+10)
            else
                if(myid .eq. 0) then
                    N   =  max_uns_record*(5*n_monitor_prc(ip))
                else
                    N   =  max_uns_record*(5*n_monitor_L)
                end if
            end if
            if(N .le. 0)    cycle

            if((myid .eq. ip) .and. (myid .ne. 0)) then
                call mpi_send(uns_record, N, mpi_dpR, 0, myid, mpi_comm_world, mpi_err)
            end if
            if((myid .ne. ip) .and. (myid .eq. 0)) then
                call mpi_recv(uns_record(1,c), N, mpi_dpR, ip, ip, &
                    &  mpi_comm_world, mpi_status, mpi_err)
            end if
            c   =  c+N/max_uns_record
        end do

        if(myid .eq. 0) then
            write(tec,*),'(',10+5*n_monitor,'ES24.16)'
            open(unit=10,file=trim(adjustl(ot)),status='old',position='append')
            do i=1,n_uns_record
                call DCOPY(10+5*n_monitor, uns_record(i,1), max_uns_record, r, 1)
                write(unit=10,fmt=trim(adjustl(tec))),r(1:10+5*n_monitor)
            end do
            close(10)
        end if
        n_uns_record=  0
    end if

    return
    end subroutine fv_wr_uns_record
!-------------------------------------------------------------------------------
!   get CAA source.
!-------------------------------------------------------------------------------
    subroutine fv_get_caa_source
    use var_kind_def
    use var_fv
    use var_global, only: err_mem,mesh_name
    use var_mesh
    use var_parallel
    use var_slv, ord=>ls_weight_order
    use var_turb, only: is_tur_cal
    use var_uns_cal, only: uns_iter
    implicit none
    character(len=100):: str
    integer(dpI):: isec,iele,ss,ee,i,j,isr,idx,ip_remote
    real   (dpR):: T(3,3),T0(3,3),mul,gra(9),S(3,3),dT(3,3),v(3),u(3),dU3, &
                &  Tlam(3,3),LHS(3,3),c(3),d,dp(3)

    do isec=mesh(0)%sec_1,mesh(0)%sec_0
        if(.not. allocated(fv(isec)%caa_temp)) &
            &  allocate(fv(isec)%caa_temp(6,sec(isec)%n_ele), stat=err_mem)
        if((.not. allocated(fv(isec)%caa_source)) .and. sec(isec)%is_int) &
            &  allocate(fv(isec)%caa_source(sec(isec)%n_ele), stat=err_mem)
    end do

!   ----------------------------------------------------------------------------
!   $\rho u u-T_vis$.
    do isec=mesh(0)%sec_1,mesh(0)%sec_0
        do iele=1,sec(isec)%n_ele
            u(1:3)  =  fv(isec)%u(2:4,iele)
            T(1:3,1)=  fv(isec)%u(1,iele)*fv(isec)%u(2,iele)*fv(isec)%u(2:4,iele)
            T(1:3,2)=  fv(isec)%u(1,iele)*fv(isec)%u(3,iele)*fv(isec)%u(2:4,iele)
            T(1:3,3)=  fv(isec)%u(1,iele)*fv(isec)%u(4,iele)*fv(isec)%u(2:4,iele)
!           T(1  ,1)=  T(1  ,1)-fv(isec)%u(5,iele)
!           T(2  ,2)=  T(2  ,2)-fv(isec)%u(5,iele)
!           T(3  ,3)=  T(3  ,3)-fv(isec)%u(5,iele)

            if(is_vis_cal) then
                mul     =  fv(isec)%mu (1   ,iele)
                gra(1:9)=  fv(isec)%gra(4:12,iele)
                dU3     = (gra(1)+gra(5)+gra(9))/3.0d0
                S(1,1)  =  gra(1)-dU3
                S(2,1)  =  0.5d0*(gra(2)+gra(4))
                S(3,1)  =  0.5d0*(gra(3)+gra(7))
                S(1,2)  =  S(2,1)
                S(2,2)  =  gra(5)-dU3
                S(3,2)  =  0.5d0*(gra(6)+gra(8))
                S(1,3)  =  S(3,1)
                S(2,3)  =  S(3,2)
                S(3,3)  =  gra(9)-dU3
                Tlam    =  2.0d0*mul*S
                if(is_tur_cal)  Tlam=  Tlam+2.0d0*fv(isec)%mu(2,iele)*S
                T       =  T-Tlam
            end if
            fv(isec)%caa_temp(1:3,iele) =  T(1:3,1)
            fv(isec)%caa_temp(4:5,iele) =  T(2:3,2)
            fv(isec)%caa_temp(6  ,iele) =  T(3  ,3)
        end do
    end do
!   $\rho u u-T_vis$.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   prepare data.
    do isr=mesh(0)%p2p_1,mesh(0)%p2p_0
        idx =  1
        do i=1,p2p(isr)%n_ele_send
            isec=  p2p(isr)%id_ele_send(1,i)
            iele=  p2p(isr)%id_ele_send(2,i)
            p2p(isr)%rsend(idx:idx+5)   =  fv(isec)%caa_temp(1:6,iele)
            idx =  idx+6
        end do
        p2p(isr)%n_send =  idx-1
        p2p(isr)%n_recv =  6*p2p(isr)%n_ele_recv
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
            if(mpi_nreq .gt. max_p2p)   stop 'Error: you are working with a bad LB.'
            call mpi_isend(p2p(isr)%rsend, p2p(isr)%n_send, mpi_dpR, ip_remote, &
                &  myid, mpi_comm_world, mpi_req(mpi_nreq), mpi_err)
        end if

        if(p2p(isr)%n_recv .gt. 0) then
            mpi_nreq=  mpi_nreq+1
            if(mpi_nreq .gt. max_p2p)   stop 'Error: you are working with a bad LB.'
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
            call DCOPY(6, p2p(isr)%rrecv(6*iele-5), 1, fv(isec)%caa_temp(1,i), 1)
        end do
    end do
!   unpack data received.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   divergence of $\rho u u-T_vis$.
    do isec=mesh(0)%sec_1,mesh(0)%sec_0
        if(.not. sec(isec)%is_int)  cycle
        do iele=1,sec(isec)%n_ele
            T0(1:3,1)   =  fv(isec)%caa_temp(1:3,iele)
            T0(2:3,2)   =  fv(isec)%caa_temp(4:5,iele)
            T0(3  ,3)   =  fv(isec)%caa_temp(6  ,iele)
            T0(1  ,2)   =  T0(2,1)
            T0(1:2,3)   =  T0(3,1:2)
            v           =  0.0d0
            if(gradient_method .eq. gradient_LSm) then
                LHS(1,1)=  sec(isec)%coe_ls(1,iele)
                LHS(2,1)=  sec(isec)%coe_ls(2,iele)
                LHS(3,1)=  sec(isec)%coe_ls(3,iele)
                LHS(1,2)=  LHS(2,1)
                LHS(2,2)=  sec(isec)%coe_ls(4,iele)
                LHS(3,2)=  sec(isec)%coe_ls(5,iele)
                LHS(1,3)=  LHS(3,1)
                LHS(2,3)=  LHS(3,2)
                LHS(3,3)=  sec(isec)%coe_ls(6,iele)
            end if

            do j=sec(isec)%iA_ls(iele),sec(isec)%iA_ls(iele+1)-1
                ss  =  sec(isec)%jA_ls(1,j)
                ee  =  sec(isec)%jA_ls(2,j)

                T(1:3,1)=  fv(ss)%caa_temp(1:3,ee)
                T(2:3,2)=  fv(ss)%caa_temp(4:5,ee)
                T(3  ,3)=  fv(ss)%caa_temp(6  ,ee)
                T(1  ,2)=  T(2,1)
                T(1:2,3)=  T(3,1:2)
                dT      =  T-T0

                if(gradient_method .eq. gradient_LS) then
                    v(1)=  dT(1,1)*sec(isec)%coe_ls(1,j)+dT(2,1)*sec(isec)%coe_ls(2,j) &
                        & +dT(3,1)*sec(isec)%coe_ls(3,j)+v(1)
                    v(2)=  dT(1,2)*sec(isec)%coe_ls(1,j)+dT(2,2)*sec(isec)%coe_ls(2,j) &
                        & +dT(3,2)*sec(isec)%coe_ls(3,j)+v(2)
                    v(3)=  dT(1,3)*sec(isec)%coe_ls(1,j)+dT(2,3)*sec(isec)%coe_ls(2,j) &
                        & +dT(3,3)*sec(isec)%coe_ls(3,j)+v(3)
                elseif(gradient_method .eq. gradient_LSm) then
                    dp(1:3) =  sec(ss)%cen(1:3,ee)-sec(isec)%cen(1:3,iele)
                    d       =  1.0d0/sqrt(dp(1)*dp(1)+dp(2)*dp(2)+dp(3)*dp(3))**ord
                    c(1:3)  =  d*(LHS(1:3,1)*dp(1)+LHS(1:3,2)*dp(2)+LHS(1:3,3)*dp(3))

                    v(1)=  dT(1,1)*c(1)+dT(2,1)*c(2)+dT(3,1)*c(3)+v(1)
                    v(2)=  dT(1,2)*c(1)+dT(2,2)*c(2)+dT(3,2)*c(3)+v(2)
                    v(3)=  dT(1,3)*c(1)+dT(2,3)*c(2)+dT(3,3)*c(3)+v(3)
                end if
            end do
            fv(isec)%duc(1:3,iele)  =  v(1:3)
        end do
    end do
    do isec=mesh(0)%sec_1,mesh(0)%sec_0
        if(.not. sec(isec)%is_int)  cycle
        do iele=1,sec(isec)%n_ele
            fv(isec)%caa_temp(1:3,iele) =  fv(isec)%duc(1:3,iele)
        end do
    end do
    do i=1,mesh(0)%n_mortar_b
        isec=  mesh(0)%mortar_LR(1,i)
        iele=  mesh(0)%mortar_LR(2,i)
        ss  =  mesh(0)%mortar_LR(3,i)
        ee  =  mesh(0)%mortar_LR(4,i)
        fv(ss)%CAA_temp(:,ee)   =  fv(isec)%CAA_temp(:,iele)
    end do
!   divergence of $\rho u u-T_vis$.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   prepare data.
    do isr=mesh(0)%p2p_1,mesh(0)%p2p_0
        idx =  1
        do i=1,p2p(isr)%n_ele_send
            isec=  p2p(isr)%id_ele_send(1,i)
            iele=  p2p(isr)%id_ele_send(2,i)
            p2p(isr)%rsend(idx:idx+2)   =  fv(isec)%CAA_temp(1:3,iele)
            idx =  idx+3
        end do
        p2p(isr)%n_send =  idx-1
        p2p(isr)%n_recv =  3*p2p(isr)%n_ele_recv
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
            if(mpi_nreq .gt. max_p2p)   stop 'Error: you are working with a bad LB.'
            call mpi_isend(p2p(isr)%rsend, p2p(isr)%n_send, mpi_dpR, ip_remote, &
                &  myid, mpi_comm_world, mpi_req(mpi_nreq), mpi_err)
        end if

        if(p2p(isr)%n_recv .gt. 0) then
            mpi_nreq=  mpi_nreq+1
            if(mpi_nreq .gt. max_p2p)   stop 'Error: you are working with a bad LB.'
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
            call DCOPY(3, p2p(isr)%rrecv(3*iele-2), 1, fv(isec)%caa_temp(1,i), 1)
        end do
    end do
!   unpack data received.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   CAA source.
    do isec=mesh(0)%sec_1,mesh(0)%sec_0
        if(.not. sec(isec)%is_int)  cycle
        fv(isec)%caa_source =  0.0d0
        do iele=1,sec(isec)%n_ele
            if(gradient_method .eq. gradient_LSm) then
                LHS(1,1)=  sec(isec)%coe_ls(1,iele)
                LHS(2,1)=  sec(isec)%coe_ls(2,iele)
                LHS(3,1)=  sec(isec)%coe_ls(3,iele)
                LHS(1,2)=  LHS(2,1)
                LHS(2,2)=  sec(isec)%coe_ls(4,iele)
                LHS(3,2)=  sec(isec)%coe_ls(5,iele)
                LHS(1,3)=  LHS(3,1)
                LHS(2,3)=  LHS(3,2)
                LHS(3,3)=  sec(isec)%coe_ls(6,iele)
            end if

            do j=sec(isec)%iA_ls(iele),sec(isec)%iA_ls(iele+1)-1
                ss  =  sec(isec)%jA_ls(1,j)
                ee  =  sec(isec)%jA_ls(2,j)

                v(1:3)  =  fv(ss)%caa_temp(1:3,ee)-fv(isec)%caa_temp(1:3,iele)

                if(gradient_method .eq. gradient_LS) then
                    fv(isec)%caa_source(iele)   =  fv(isec)%caa_source(iele) &
                        & +v(1)*sec(isec)%coe_ls(1,j)+v(2)*sec(isec)%coe_ls(2,j) &
                        & +v(3)*sec(isec)%coe_ls(3,j)
                elseif(gradient_method .eq. gradient_LSm) then
                    dp(1:3) =  sec(ss)%cen(1:3,ee)-sec(isec)%cen(1:3,iele)
                    d       =  1.0d0/sqrt(dp(1)*dp(1)+dp(2)*dp(2)+dp(3)*dp(3))**ord
                    c(1:3)  =  d*(LHS(1:3,1)*dp(1)+LHS(1:3,2)*dp(2)+LHS(1:3,3)*dp(3))
                    fv(isec)%caa_source(iele)   =  fv(isec)%caa_source(iele) &
                        & +v(1)*c(1)+v(2)*c(2)+v(3)*c(3)
                end if
            end do
        end do
    end do
!   CAA source.
!   ----------------------------------------------------------------------------

    write(str,*),uns_iter
    str =  trim(adjustl(mesh_name))//'_caa_'//trim(adjustl(str))//'.cgns'
    call fv_wr_cgns_cc(6, str)

    return
    end subroutine fv_get_caa_source
!-------------------------------------------------------------------------------
!   get entropy generation rate.
!-------------------------------------------------------------------------------
    subroutine fv_get_entropy_generation
    use var_kind_def
    use var_air
    use var_fv
    use var_global, only: mesh_name
    use var_mesh
    use var_slv, only: gradient_method,gradient_LS,gradient_LSm,ord=>LS_weight_order
    use var_turb, only: is_tur_cal,is_Spalart_QCR
    implicit none
    character(len=100):: str
    integer(dpI):: isec,iele,j,ss,ee
    real   (dpR):: kq,d0(3),d(3),LHS(3,3),mul,gra(9),dU3,S(3,3),Tlam(3,3),dp(3),L, &
                &  mut,Ttur(3,3),c(3)

    do isec=mesh(0)%sec_1,mesh(0)%sec_0
        if(.not. sec(isec)%is_int)  cycle
        fv(isec)%rhs=  0.0d0

        if(gradient_method .eq. gradient_LS) then
            do iele=1,sec(isec)%n_ele
                kq  =  cp*fv(isec)%mu(1,iele)/prL
                if(is_tur_cal)  kq  =  kq+cp*fv(isec)%mu(2,iele)/prT
                d0(1:3) =  kq*fv(isec)%gra(16:18,iele)

                do j=sec(isec)%iA_ls(iele),sec(isec)%iA_ls(iele+1)-1
                    ss  =  sec(isec)%jA_ls(1,j)
                    ee  =  sec(isec)%jA_ls(2,j)
                    kq  =  cp*fv(ss)%mu(1,ee)/prL
                    if(is_tur_cal)  kq  =  kq+cp*fv(ss)%mu(2,ee)/prT
                    d   =  kq*fv(ss)%gra(16:18,ee)-d0
                    fv(isec)%rhs(1,iele)=  fv(isec)%rhs(1,iele) &
                        & +d(1)*sec(isec)%coe_ls(1,j)+d(2)*sec(isec)%coe_ls(2,j) &
                        & +d(3)*sec(isec)%coe_ls(3,j)
                end do
            end do
        elseif(gradient_method .eq. gradient_LSm) then
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
                kq      =  cp*fv(isec)%mu(1,iele)/prL
                if(is_tur_cal)  kq  =  kq+cp*fv(isec)%mu(2,iele)/prT
                d0(1:3) =  kq*fv(isec)%gra(16:18,iele)
                do j=sec(isec)%iA_ls(iele),sec(isec)%iA_ls(iele+1)-1
                    ss  =  sec(isec)%jA_ls(1,j)
                    ee  =  sec(isec)%jA_ls(2,j)

                    dp(1:3) =  sec(ss)%cen(1:3,ee)-sec(isec)%cen(1:3,iele)
                    L       =  1.0d0/sqrt(dp(1)*dp(1)+dp(2)*dp(2)+dp(3)*dp(3))**ord
                    c(1:3)  =  L*(LHS(1:3,1)*dp(1)+LHS(1:3,2)*dp(2)+LHS(1:3,3)*dp(3))

                    kq  =  cp*fv(ss)%mu(1,ee)/prL
                    if(is_tur_cal)  kq  =  kq+cp*fv(ss)%mu(2,ee)/prT
                    d   =  kq*fv(ss)%gra(16:18,ee)-d0
                    fv(isec)%rhs(1,iele)=  fv(isec)%rhs(1,iele) &
                                        & +d(1)*c(1)+d(2)*c(2)+d(3)*c(3)
                end do
            end do
        else
            stop 'Error: gradient method not supported, fv_get_entropy_generation.'
        end if
    end do

    do isec=mesh(0)%sec_1,mesh(0)%sec_0
        if(.not. sec(isec)%is_int)  cycle
        do iele=1,sec(isec)%n_ele
            mul     =  fv(isec)%mu(1,iele)
            gra(1:9)=  fv(isec)%gra(4:12,iele)
            dU3     = (gra(1)+gra(5)+gra(9))/3.0d0
            S(1,1)  =  gra(1)-dU3
            S(2,1)  =  0.5d0*(gra(2)+gra(4))
            S(3,1)  =  0.5d0*(gra(3)+gra(7))
            S(1,2)  =  S(2,1)
            S(2,2)  =  gra(5)-dU3
            S(3,2)  =  0.5d0*(gra(6)+gra(8))
            S(1,3)  =  S(3,1)
            S(2,3)  =  S(3,2)
            S(3,3)  =  gra(9)-dU3
            Tlam    =  2.0d0*mul*S
            if(is_tur_cal) then
                mut =  fv(isec)%mu(2,iele)
                Ttur=  2.0d0*mut*S
                if(is_Spalart_QCR)  call get_QCR_stress(mut, gra, S, Ttur)
                Tlam=  Tlam+Ttur
            end if

            fv(isec)%rhs(2,iele)=  Tlam(1,1)*S(1,1)+Tlam(2,1)*S(2,1)+Tlam(3,1)*S(3,1) &
                & +Tlam(1,2)*S(1,2)+Tlam(2,2)*S(2,2)+Tlam(3,2)*S(3,2) &
                & +Tlam(1,3)*S(1,3)+Tlam(2,3)*S(2,3)+Tlam(3,3)*S(3,3)
            mul =  1.0d0/(fv(isec)%u(1,iele)*fv(isec)%t(iele))
            fv(isec)%rhs(1,iele)= -fv(isec)%rhs(1,iele)*mul
            fv(isec)%rhs(2,iele)=  fv(isec)%rhs(2,iele)*mul
        end do
    end do

    str =  trim(adjustl(mesh_name))//'_entropy.cgns'
    call fv_wr_cgns_cc(7, str)

    return
    end subroutine fv_get_entropy_generation
!-------------------------------------------------------------------------------
!   output unsteady solution.
!-------------------------------------------------------------------------------
    subroutine fv_uns_wr_results
    use var_kind_def
    use var_global, only: is_output_cc,is_output_vc,mesh_name
    use var_uns_cal
    implicit none
    character(len=80):: str
    logical(dpL):: ltmp

    call fv_wr_uns_record
    call fv_get_uns_ra(0)

!   save cell-centered solution at some timesteps.
    ltmp=((mod(uns_iter-uns_ra_begin, uns_write_delta) .eq. 0) .and. &
        & (uns_iter+uns_write_delta .le. max_uns_iter) .and. &
        & (uns_iter .ge. uns_ra_begin) .and. is_output_cc) .or. &
        & (uns_iter .gt. max_uns_iter-uns_n_stencil+1)
    if(ltmp) then
        write(str,*),uns_iter
        str =  trim(adjustl(mesh_name))//'_uns_cc_'//trim(adjustl(str))//'.cgns'
        call fv_wr_cgns_cc(2, str)
    end if

!   save vertex-centered solution at some timesteps.
    ltmp=((mod(uns_iter-uns_ra_begin, uns_write_delta) .eq. 0) .and. &
        & (uns_iter+uns_write_delta .le. max_uns_iter) .and. &
        & (uns_iter .ge. uns_ra_begin) .and. is_output_vc) .or. &
        & (uns_iter .eq. max_uns_iter)
    if(ltmp) then
        write(str,*),uns_iter
        str =  trim(adjustl(mesh_name))//'_uns_vc_'//trim(adjustl(str))//'.cgns'
        call fv_wr_cgns(1, str)
    end if

!   if(uns_iter .eq. max_uns_iter)  call fv_get_caa_source

    return
    end subroutine fv_uns_wr_results
