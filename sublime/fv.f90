!-------------------------------------------------------------------------------
!   module FV. 
!-------------------------------------------------------------------------------
    module var_fv
        use var_kind_def
        use var_mesh, only: max_sec_lev
        use var_prec, only: type_prec
        use var_gmres
        implicit none

        real   (dpR),parameter:: JST_k2 =  0.5d0
        real   (dpR),parameter:: JST_k4 =  1.0d0/3.2d1

        logical(dpL):: is_kep               =  .false.
        logical(dpL):: is_shock_sensor      =  .false.
        logical(dpL):: is_jst_su2           =  .false.
        logical(dpL):: is_limiter_on        =  .false.
        logical(dpL):: is_rslv_LMa          =  .true.
        logical(dpL):: is_mg                =  .false.
        logical(dpL):: is_cht               =  .false.
        integer(dpI):: mlev_fv              =  1
        integer(dpI):: order(1000)          =  1
        integer(dpI):: rhs_lev(0:9)         =  1
        integer(dpI):: freeze_limiter       = -1
        real   (dpR):: res_NS               =  0.0d0
        real   (dpR):: res_RANS             =  0.0d0
        real   (dpR):: res_NS_now           =  0.0d0
        real   (dpR):: res_NS_old           =  0.0d0
        real   (dpR):: res_NS_1st           =  0.0d0
        real   (dpR):: LHS_efix             =  5.0d-2
        real   (dpR):: min_yp               =  1.0d10
        real   (dpR):: max_yp               =  0.0d0
        type type_fv_section
            logical(dpL):: is_reot  =  .false.
            integer(dpI):: order    =  1
            real   (dpR),allocatable:: u  (:,:)
            real   (dpR),allocatable:: t  (:  )
            real   (dpR),allocatable:: uc (:,:)
            real   (dpR),allocatable:: uc0(:,:)
            real   (dpR),allocatable:: duc(:,:)
            real   (dpR),allocatable:: rhs(:,:)
            real   (dpR),allocatable:: rv0(:,:)
            real   (dpR),allocatable:: LHS_s(:)
            real   (dpR),allocatable:: mu(:,:),gra(:,:)
            real   (dpR),allocatable:: heat_flux(:)
            real   (dpR),allocatable:: ss_lap(:,:)
            real   (dpR),allocatable:: shock_sensor(:)
            real   (dpR),allocatable:: turb(:,:),turb_gra(:,:),q(:)
            real   (dpR),allocatable:: intermittency(:)
            real   (dpR),allocatable:: DDES_quality(:,:)
            real   (dpR),allocatable:: Re_theta(:),bnd_solution(:,:)
            real   (dpR),allocatable:: uns_uc_1(:,:),uns_uc_2(:,:),uns_uc_3(:,:)
            real   (dpR),allocatable:: uns_turb(:,:)

!           running-average of cell-centered unsteady solutions.
            real   (dpR),allocatable:: uns_ra(:,:)
            real   (dpR),allocatable:: stress(:,:)

            real   (dpR),allocatable:: adj(:,:),adj0(:,:),adj_rhs(:,:),adj_pJu(:,:)
            real   (dpR),allocatable:: adj_rhs_um(:,:),adj_rhs_ug(:,:)
            real   (dpR),allocatable:: adj_JST(:,:),adj_rhs_gradient(:,:)
            real   (dpR),allocatable:: um_ui(:,:),ug_ui(:,:),ub_bc(:,:),J_bc(:,:)
            real   (dpR),allocatable:: J_u(:,:)

            real   (dpR),allocatable:: ff(:,:),uc_f2c(:,:)

            real   (dpR),allocatable:: bc(:,:),uL(:,:),uR(:,:)

            integer(dpI),allocatable:: nearest_wall(:)

!           for IMEX.
            integer(dpI),allocatable:: IMEX_ID(:)
            real   (dpR),allocatable:: IMEX_k  (:,:),IMEX_kt (:,:)
            real   (dpR),allocatable:: IMEX_uc0(:,:),IMEX_uct(:,:)

            real   (dpR),allocatable:: caa_source(:),caa_temp(:,:)
        end type type_fv_section
        type(type_fv_section),save:: fv(max_sec_lev*5)
        type(type_prec      ),save:: fv_prec(0:5),fv_rans_prec,fv_imex_prec
        type(type_gmres),save:: fv_gmres

        integer(dpI),allocatable:: ele_send(:,:),ele_recv(:,:)
        real   (dpR),allocatable:: r_send(:),r_recv(:),r_recv_all(:)
        real   (dpR),allocatable:: fv_vertex_variables(:,:)
        character(len=80):: vc_sol_name(20)
        character(len=80):: ini_file=  'NONE'
        integer(dpI):: n_vc_to_write=  0
    end module var_fv
!-------------------------------------------------------------------------------
!   module FV (array). 
!-------------------------------------------------------------------------------
    module var_fv_array
        use var_kind_def
        use var_mesh, only: max_sec_lev
        use var_prec, only: type_prec
        use var_gmres
        implicit none

        real   (dpR),allocatable:: fv_order(:)
        real   (dpR),allocatable:: fv_u  (:,:)
        real   (dpR),allocatable:: fv_t  (:  )
        real   (dpR),allocatable:: fv_uc (:,:)
        real   (dpR),allocatable:: fv_uc0(:,:)
        real   (dpR),allocatable:: fv_duc(:,:)
        real   (dpR),allocatable:: fv_rhs(:,:)
        real   (dpR),allocatable:: fv_rv0(:,:)
        real   (dpR),allocatable:: fv_LHS_s(:)
        real   (dpR),allocatable:: fv_mu(:,:),fv_gra(:,:)
        real   (dpR),allocatable:: fv_heat_flux(:)
        real   (dpR),allocatable:: fv_ss_lap(:,:)
        real   (dpR),allocatable:: fv_shock_sensor(:)
        real   (dpR),allocatable:: fv_turb(:,:),fv_turb_gra(:,:),fv_q(:)
        !real   (dpR),allocatable:: intermittency(:)
        !real   (dpR),allocatable:: DDES_quality(:,:)
        real   (dpR),allocatable:: fv_Re_theta(:),fv_bnd_solution(:,:)
        !real   (dpR),allocatable:: uns_uc_1(:,:),uns_uc_2(:,:),uns_uc_3(:,:)
        !real   (dpR),allocatable:: uns_turb(:,:)

!       running-average of cell-centered unsteady solutions.
        !real   (dpR),allocatable:: uns_ra(:,:)
        !real   (dpR),allocatable:: stress(:,:)

        !real   (dpR),allocatable:: adj(:,:),adj0(:,:),adj_rhs(:,:),adj_pJu(:,:)
        !real   (dpR),allocatable:: adj_rhs_um(:,:),adj_rhs_ug(:,:)
        !real   (dpR),allocatable:: adj_JST(:,:),adj_rhs_gradient(:,:)
        !real   (dpR),allocatable:: um_ui(:,:),ug_ui(:,:),ub_bc(:,:),J_bc(:,:)
        !real   (dpR),allocatable:: J_u(:,:)

        !real   (dpR),allocatable:: ff(:,:),uc_f2c(:,:)

        real   (dpR),allocatable:: fv_bc(:,:),fv_uL(:,:),fv_uR(:,:)

        !integer(dpI),allocatable:: nearest_wall(:)

!       for IMEX.
        !integer(dpI),allocatable:: IMEX_ID(:)
        !real   (dpR),allocatable:: IMEX_k  (:,:),IMEX_kt (:,:)
        !real   (dpR),allocatable:: IMEX_uc0(:,:),IMEX_uct(:,:)

        !real   (dpR),allocatable:: caa_source(:),caa_temp(:,:)
        !end type type_fv_section
    end module var_fv_array

!-------------------------------------------------------------------------------
!   setup solver.
!-------------------------------------------------------------------------------
    subroutine fv_set_slv
    use var_kind_def
    use var_bndv, only: Re_fs
    use var_fv
    use var_global
    use var_mg, only: mlev
    use var_parallel
    use var_slv
    use var_turb, only: is_DDES,is_LES,is_tur_cal
    use var_uns_cal, only: max_uns_iter
    implicit none
    logical(dpL):: lbuf(100)
    integer(dpI):: lev,io_err,ibuf(10000)

!   ----------------------------------------------------------------------------
!   read and broadcast the FV parameters.
    namelist /fv_parameters/    is_kep,order,rhs_lev,is_limiter_on,freeze_limiter, &
                            &   is_mg,is_shock_sensor

    if((myid .eq. 0) .and. is_has_cfg) then
        open(unit=10,file=trim(adjustl(cfg_file)))
        read(unit=10, nml=fv_parameters, iostat=io_err)
        if(io_err .gt. 0)   stop 'Error: fails to read namelist:fv_parameters.'
        close(10)
    end if

    lbuf(1:4   )= (/is_kep, is_limiter_on, is_mg, is_shock_sensor/)
    ibuf(1:1011)= (/order, rhs_lev, freeze_limiter/)
    call mpi_bcast(lbuf, 4   , mpi_dpL, 0, mpi_comm_world, mpi_err)
    call mpi_bcast(ibuf, 1011, mpi_dpI, 0, mpi_comm_world, mpi_err)
    if(myid .ne. 0) then
        is_kep          =  lbuf(1)
        is_limiter_on   =  lbuf(2)
        is_mg           =  lbuf(3)
        is_shock_sensor =  lbuf(4)
        order(1:1000)   =  ibuf(1:1000)
        rhs_lev(0:9)    =  ibuf(1001:1010)
        freeze_limiter  =  ibuf(1011)
    end if
!   read and broadcast the FV parameters.
!   ----------------------------------------------------------------------------

    solver_lev      =  max(1, min(solver_lev, 3))
    is_mg           =  is_mg .and. (mlev .gt. 1)
    if(is_mg)   call fv_set_coarse_mesh
    max_iter_lev    =  max(0, max_iter_lev)
    lev0_cal        =  max(0, min(mlev_fv-1, lev0_cal))
    is_limiter_on   =  is_limiter_on .and. (rhs_lev(0) .gt. 1)
    if(freeze_limiter .le. 0)   freeze_limiter  =  huge(1)

    do lev=lev0_cal,mlev_fv-1
        call fv_set_mesh_flow(lev)

        if((max_iter_lev(lev_out) .gt. 0) .or. (max_uns_iter .gt. 0)) then
            if(solver_lev(lev) .eq. solver_exp) then
                call fv_set_exp_slv(lev)
            elseif(solver_lev(lev) .le. solver_gmres_d) then
                call fv_set_imp_slv(lev, fv_prec(lev))
if(sw_slave) then
                call fv_set_lhs_imp_sw(lev, fv_prec(lev))
end if
            end if
        end if

        if(lev .gt. lev0_cal) then
            if(rhs_lev(0) .eq. 1) then
                rhs_lev(lev)=  2
            else
                rhs_lev(lev)=  rhs_lev(0)
            end if
        end if
    end do
    tot_ite =  0
    if(is_tur_cal) then
        BL_thick=  0.382d0*L_ref/Re_fs**0.2d0
        BL_thick=  0.1d0/sqrt(Re_fs)
    else
        BL_thick=  4.910d0*L_ref/sqrt(Re_fs)
    end if

    if(is_DDES .or. is_LES) call fv_get_gls
    call mesh_get_vc_weights

    return
    end subroutine fv_set_slv
!-------------------------------------------------------------------------------
!   fv set mesh.
!-------------------------------------------------------------------------------
    subroutine fv_set_mesh_flow(lev)
    use var_kind_def
    use var_fv
    use var_fv_array
    use var_global, only: err_mem,sw_slave
    use var_mesh
    use var_slv, only: lev0_cal,LDA_res_record
    use var_turb
    use var_parallel, only: myid
    implicit none
    integer(dpI),intent(in):: lev
    integer(dpI):: isec,n_ele

    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(sec(isec)%is_bnd)    cycle
        if(lev .eq. 0) then
            fv(isec)%order  =  order(sec(isec)%id_sec_i)
        else
            fv(isec)%order  =  1
        end if
    end do

    tot_ele = 0
    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        n_ele   =  sec(isec)%n_ele
        tot_ele =  tot_ele + n_ele
        allocate(fv(isec)%u  (5 ,n_ele), stat=err_mem)
        allocate(fv(isec)%uc (5 ,n_ele), stat=err_mem)
        allocate(fv(isec)%uc0(5 ,n_ele), stat=err_mem)
        allocate(fv(isec)%rhs(5 ,n_ele), stat=err_mem)
        allocate(fv(isec)%duc(5 ,n_ele), stat=err_mem)
        allocate(fv(isec)%gra(19,n_ele), stat=err_mem)
        allocate(fv(isec)%t  (   n_ele), stat=err_mem)
        if(lev .gt. lev0_cal) then
            allocate(fv(isec)%uc_f2c(5,n_ele), stat=err_mem)
            allocate(fv(isec)%ff    (5,n_ele), stat=err_mem)
        end if

        ! if(sec(isec)%is_bnd) then
            allocate(fv(isec)%uL          (5,n_ele), stat=err_mem)
            allocate(fv(isec)%uR          (5,n_ele), stat=err_mem)
            allocate(fv(isec)%bnd_solution(2,n_ele), stat=err_mem)
            fv(isec)%uL             =  0.0d0
            fv(isec)%uR             =  0.0d0
            fv(isec)%bnd_solution   =  0.0d0
        ! end if

        fv(isec)%u  =  0.0d0
        fv(isec)%uc =  0.0d0
        fv(isec)%rhs=  0.0d0
        fv(isec)%duc=  0.0d0
        fv(isec)%gra=  0.0d0
        fv(isec)%t  =  0.0d0

        if(is_tur_cal) then
            allocate(fv(isec)%mu(2,n_ele), stat=err_mem)
        else
            allocate(fv(isec)%mu(1,n_ele), stat=err_mem)
        end if
        fv(isec)%mu =  0.0d0
        if(is_RANS) then
            if(RANS_model .eq. 1) then
                allocate(fv(isec)%turb    (1,n_ele), stat=err_mem)
                allocate(fv(isec)%turb_gra(3,n_ele), stat=err_mem)
            elseif(RANS_model .eq. 2) then
                allocate(fv(isec)%turb    (1,n_ele), stat=err_mem)
                allocate(fv(isec)%turb_gra(3,n_ele), stat=err_mem)
            elseif(RANS_model .eq. 3) then
                allocate(fv(isec)%turb    (2,n_ele), stat=err_mem)
                allocate(fv(isec)%turb_gra(3,n_ele), stat=err_mem)
            elseif(RANS_model .eq. 4) then
                allocate(fv(isec)%turb    (2,n_ele), stat=err_mem)
                allocate(fv(isec)%turb_gra(3,n_ele), stat=err_mem)
            elseif(RANS_model .eq. RANS_SST) then
                allocate(fv(isec)%turb    (2,n_ele), stat=err_mem)
                allocate(fv(isec)%turb_gra(6,n_ele), stat=err_mem)
            end if
            fv(isec)%turb_gra   =  0.0d0
            LDA_res_record      =  11
        else
            LDA_res_record      =  10
        end if

        if(transition_model .ge. transition_PTM) &
            &  allocate(fv(isec)%intermittency(sec(isec)%n_ele), stat=err_mem)
    end do

if(sw_slave) then
    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        allocate(fv_u    (5, tot_ele), stat=err_mem)
        allocate(fv_uc   (5 ,tot_ele), stat=err_mem)
        allocate(fv_uc0  (5 ,tot_ele), stat=err_mem)
        allocate(fv_rhs  (5 ,tot_ele), stat=err_mem)
        allocate(fv_duc  (5 ,tot_ele), stat=err_mem)
        allocate(fv_gra  (19,tot_ele), stat=err_mem)
        allocate(fv_t    (   tot_ele), stat=err_mem)
        if(is_tur_cal) then
            allocate(fv_mu(2,tot_ele), stat=err_mem)
        else
            allocate(fv_mu(1,tot_ele), stat=err_mem)
        end if
        if(is_RANS) then
            if(RANS_model .eq. 1) then
                allocate(fv_turb    (1,tot_ele), stat=err_mem)
                allocate(fv_turb_gra(3,tot_ele), stat=err_mem)
            elseif(RANS_model .eq. 2) then
                allocate(fv_turb    (1,tot_ele), stat=err_mem)
                allocate(fv_turb_gra(3,tot_ele), stat=err_mem)
            elseif(RANS_model .eq. 3) then
                allocate(fv_turb    (2,tot_ele), stat=err_mem)
                allocate(fv_turb_gra(3,tot_ele), stat=err_mem)
            elseif(RANS_model .eq. 4) then
                allocate(fv_turb    (2,tot_ele), stat=err_mem)
                allocate(fv_turb_gra(3,tot_ele), stat=err_mem)
            elseif(RANS_model .eq. RANS_SST) then
                allocate(fv_turb    (2,tot_ele), stat=err_mem)
                allocate(fv_turb_gra(6,tot_ele), stat=err_mem)
            end if
        end if

        allocate(fv_uL          (5,tot_ele), stat=err_mem)
        allocate(fv_uR          (5,tot_ele), stat=err_mem)
        allocate(fv_bnd_solution(2,tot_ele), stat=err_mem)
    end do
    if(myid .eq. 0) write(*,*), "total_ele", tot_ele
    call sec_struct_to_array(lev)
    call transform_parameter_to_real(lev)
end if

    call fv_set_p2p(lev)

    return
    end subroutine fv_set_mesh_flow
!-------------------------------------------------------------------------------
!   set the p2p, fv.
!-------------------------------------------------------------------------------
    subroutine fv_set_p2p(lev)
    use var_kind_def
    use var_global, only: err_mem
    use var_mesh
    implicit none
    integer(dpI),intent(in):: lev
    integer(dpI):: isr,M,N

    do isr=mesh(lev)%p2p_1,mesh(lev)%p2p_0
        M   =  19*p2p(isr)%n_ele_send
        if(M .gt. 0)    allocate(p2p(isr)%rsend(M), stat=err_mem)

        N   =  19*p2p(isr)%n_ele_recv
        if(N .gt. 0)    allocate(p2p(isr)%rrecv(N), stat=err_mem)
    end do

    return
    end subroutine fv_set_p2p
!-------------------------------------------------------------------------------
!   transfrom u into uc.
!-------------------------------------------------------------------------------
    subroutine fv_lev_u2uc(lev,mode)
    use var_kind_def
    use var_fv
    use var_mesh, only: mesh,sec
    implicit none
    integer(dpI),intent(in):: lev,mode
    integer(dpI):: isec,ele1,ele0

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
        if(ele0 .lt. ele1)  cycle

        call u_to_uc(ele0-ele1+1, fv(isec)%u(1,ele1), fv(isec)%uc(1,ele1), &
            &  fv(isec)%t(ele1), fv(isec)%mu(1,ele1))
    end do

    return
    end subroutine fv_lev_u2uc
!-------------------------------------------------------------------------------
!   transfrom uc into u.
!-------------------------------------------------------------------------------
    subroutine fv_lev_uc2u(lev,mode)
    use var_kind_def
    use var_fv
    use var_mesh, only: mesh,sec
    implicit none
    integer(dpI),intent(in):: lev,mode
    integer(dpI):: isec,ele1,ele0

    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(sec(isec)%is_ghost .or. sec(isec)%is_bnd)    cycle
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
        if(ele0 .lt. ele1)  cycle

        call uc_to_u(ele0-ele1+1, fv(isec)%uc(1,ele1), fv(isec)%u(1,ele1), &
            &  fv(isec)%t(ele1), fv(isec)%mu(1,ele1))
    end do

    return
    end subroutine fv_lev_uc2u
!-------------------------------------------------------------------------------
!   transfrom uc into u. Sunway version
!-------------------------------------------------------------------------------
    subroutine fv_lev_uc2u_sw(lev)
    use var_kind_def
    use var_fv
    use var_mesh, only: mesh,sec,perm,cellNum,tot_ele
    use var_fv_array
    use var_air, only: rr,gk1,t_suth
    use var_turb, only: is_tur_cal
    use var_sec_array
    use var_global_real, only: is_tur_cal_r
    implicit none
    integer(dpI),intent(in):: lev
    integer(dpI):: isec,ele1,ele0,iele,ID
    real   (dpR):: ke,mu

! call fv_struct_to_array(lev)

    ! do iele = 1,tot_ele
    !     if(sec_is_ghost(iele) .eq. 1 .or. sec_is_bnd(iele) .eq. 1) cycle
    !         ID = iele
    !         fv_u(1,ID)   = fv_uc(1,ID)
    !         fv_u(2:4,ID) = fv_uc(2:4,ID)/fv_uc(1,ID)
    !         ke = 0.5d0*(fv_u(2,ID)**2+fv_u(3,ID)**2+fv_u(4,ID)**2)
    !         fv_u(5,ID) = abs(fv_uc(5,ID)-fv_uc(1,ID)*ke)*gk1
    !         fv_t(  ID) = fv_u(5,ID)/(rr*fv_uc(1,ID))
    !         mu = 1.461d-6*sqrt(fv_t(ID)**3)/(fv_t(ID)+t_suth)
    !         fv_mu(1,ID) = mu           
    !     ! end if
    ! end do

    call uc2u_host(tot_ele, fv_u, fv_uc, fv_t, fv_mu, sec_is_int, &
        & rr, gk1, t_suth, is_tur_cal_r)

! call fv_array_to_struct(lev)

    return
    end subroutine fv_lev_uc2u_sw
!-------------------------------------------------------------------------------
!   transfrom u into uc. Sunway version
!-------------------------------------------------------------------------------
    subroutine fv_lev_u2uc_sw(lev)
    use var_kind_def
    use var_fv
    use var_mesh, only: mesh,sec,perm,cellNum,tot_ele
    use var_fv_array
    use var_air, only: rr,gk1,gk,t_suth
    use var_turb, only: is_tur_cal
    use var_sec_array
    implicit none
    integer(dpI),intent(in):: lev
    integer(dpI):: isec,ele1,ele0,iele,ID
    real   (dpR):: ke,mu,a,h

! call fv_struct_to_array(lev)

    do iele = 1,tot_ele
        if(sec_is_int(iele) .eq. 1) cycle
        a = sqrt(abs(gk*fv_u(5,iele)/fv_u(1,iele)))
        ke = 0.5d0*(fv_u(2,iele)**2+fv_u(3,iele)**2+fv_u(4,iele)**2)
        h = a*a/gk1+ke
        fv_uc(1,iele) = fv_u(1,iele)
        fv_uc(2:4,iele) = fv_u(1,iele)*fv_u(2:4,iele)
        fv_uc(5,iele) = fv_u(5,iele)/gk1+ke*fv_u(1,iele)
        fv_t(iele) = fv_u(5,iele)/(rr*fv_u(1,iele))
        fv_mu(1,iele) = 1.461d-6*sqrt(fv_t(iele)**3)/(fv_t(iele)+t_suth)          
    end do

! call fv_array_to_struct(lev)

    return
    end subroutine fv_lev_u2uc_sw
!-------------------------------------------------------------------------------
!   prefld.
!-------------------------------------------------------------------------------
    subroutine fv_prefld
    use var_kind_def
    use var_bndv, only: u_fs,vdir_inl
    use var_fv
    use var_global
    use var_mesh
    use var_slv, only: is_vis_cal
    use var_turb
    implicit none
    integer(dpI):: lev,isec,iele,npe,i,ID
    real   (dpR):: c(3),ome(3),vg(3)

    do isec=mesh(0)%sec_1,mesh(mlev_fv-1)%sec_0
        npe =  sec(isec)%npe
        do iele=1,sec(isec)%n_ele
            c   =  0.0d0
            do i=1,npe
                c(1:n_dim)  =  c(1:n_dim)+mesh(0)%xyz(1:n_dim,sec(isec)%n2e(i,iele))
            end do
            c   =  c/real(npe, dpR)

            if(is_flow_int) then
                fv(isec)%u(1:5,iele)= (/rref, uref*vdir_inl(1:3,1), pref/)
            else
                fv(isec)%u(1:5,iele)=  u_fs(1:5)
            end if

            if(RANS_model .eq. RANS_SA) then
                fv(isec)%turb(1,iele)   =  nut_ref
            elseif(RANS_model .eq. RANS_WA) then
                fv(isec)%turb(1,iele)   =  nut_ref
            elseif(RANS_model .eq. RANS_SAM) then
                fv(isec)%turb(1,iele)   =  nut_ref
                fv(isec)%turb(2,iele)   =  1.0d0
            elseif(RANS_model .eq. RANS_SAC) then
                fv(isec)%turb(1,iele)   =  nut_ref
                fv(isec)%turb(2,iele)   =  0.0d0
            elseif(RANS_model .eq. RANS_SST) then
                fv(isec)%turb(1,iele)   =  k_ref
                fv(isec)%turb(2,iele)   =  o_ref
            end if

!           call iv2d((/0.0d0, 0.0d0/), c, fv(isec)%u(1,iele))
        end do
    end do

    call fv_prefld_along_axis

!   ----------------------------------------------------------------------------
!   change the velocity to relative frame.
    do isec=mesh(0)%sec_1,mesh(0)%sec_0
        if(.not. sec(isec)%is_int)  cycle
        ID  =  sec(isec)%ID_sec_g
        if(sec_motion_type(ID) .eq. 0)  cycle

        if(sec_motion_type(ID) .eq. 1) then
            do iele=1,sec(isec)%n_ele
                fv(isec)%u(2:4,iele)=  fv(isec)%u(2:4,iele) &
                    & -sec_motion_speed(ID)*translation_axis(1:3)
            end do
        elseif(sec_motion_type(ID) .eq. 2) then
            do iele=1,sec(isec)%n_ele
                ome =  sec_motion_speed(isec)*rotation_axis
                call crs_prd(ome, sec(isec)%cen(1,iele), vg)
                fv(isec)%u(2:4,iele)=  fv(isec)%u(2:4,iele)-vg(1:3)
            end do
        end if
    end do
!   change the velocity to relative frame.
!   ----------------------------------------------------------------------------

    if(is_vis_cal) then
        call fv_get_dnw
        do isec=mesh(0)%sec_1,mesh(mlev_fv-1)%sec_0
            if(.not. sec(isec)%is_int)  cycle
            if(.not. allocated(sec(isec)%dnw))  cycle
            do iele=1,sec(isec)%n_ele
                c(1)=  min(sec(isec)%dnw(iele)/BL_thick, 1.0d0)
                fv(isec)%u(2:4,iele)=  fv(isec)%u(2:4,iele)*c(1)
                if(is_RANS) fv(isec)%turb(:,iele)   =  fv(isec)%turb(:,iele)*c(1)
            end do
        end do
    end if

!   ----------------------------------------------------------------------------
!   change the velocity to absolution frame.
    do isec=mesh(0)%sec_1,mesh(0)%sec_0
        if(.not. sec(isec)%is_int)  cycle
        ID  =  sec(isec)%ID_sec_g
        if(sec_motion_type(ID) .eq. 0)  cycle

        if(sec_motion_type(ID) .eq. 1) then
            do iele=1,sec(isec)%n_ele
                fv(isec)%u(2:4,iele)=  fv(isec)%u(2:4,iele) &
                    & +sec_motion_speed(ID)*translation_axis(1:3)
            end do
        elseif(sec_motion_type(ID) .eq. 2) then
            do iele=1,sec(isec)%n_ele
                ome =  sec_motion_speed(isec)*rotation_axis
                call crs_prd(ome, sec(isec)%cen(1,iele), vg)
                fv(isec)%u(2:4,iele)=  fv(isec)%u(2:4,iele)+vg(1:3)
            end do
        end if
    end do
!   change the velocity to absolution frame.
!   ----------------------------------------------------------------------------

    if(ini_file(1:4) .ne. 'NONE')   call fv_prefld_cc(ini_file)

    do lev=0,mlev_fv-1
        call fv_lev_u2uc(lev, 0)
        call fv_bnd_parallel(lev, .false.)
        if(is_RANS)     call fv_rans_bnd_parallel(lev, .true.)
    end do

    return
    end subroutine fv_prefld
!-------------------------------------------------------------------------------
!   initialize the solution field along the axis.
!-------------------------------------------------------------------------------
    subroutine fv_prefld_along_axis
    use var_kind_def
    use var_fv
    use var_global, only: is_has_cfg,cfg_file,rotation_axis
    use var_mesh
    use var_parallel
    implicit none
    integer(dpI):: isec,iele,n_axis,i,io_err
    real   (dpR):: u_axis(6,100),a,r

    namelist /ini/  u_axis,ini_file
    if(.not. is_has_cfg)    return

    u_axis  = -1.0d10
    if(myid .eq. 0) then
        open(unit=10,file=trim(adjustl(cfg_file)))

        read(unit=10,nml=ini,iostat=io_err)
        if(io_err .gt. 0)   stop 'Error: fails to read namelist: ini.'

        close(10)
    end if
    call mpi_bcast(u_axis, 600, mpi_dpR, 0, mpi_comm_world, mpi_err)
    call mpi_bcast(ini_file, 80, mpi_character, 0, mpi_comm_world, mpi_err)
    if(u_axis(1,1) .le. -1.0d8) return
    do n_axis=1,99
        if((u_axis(1,n_axis) .gt. -1.0d8) .and. (u_axis(1,n_axis+1) .le. -1.0d8))   exit
    end do
    call dqsortcols(.true., 1, n_axis, 1, 6, u_axis)

    do isec=mesh(0)%sec_1,mesh(mlev_fv-1)%sec_0
    do iele=1,sec(isec)%n_ele
        a   =  sec(isec)%cen(1,iele)*rotation_axis(1) &
            & +sec(isec)%cen(2,iele)*rotation_axis(2) &
            & +sec(isec)%cen(3,iele)*rotation_axis(3)
        if(a .le. u_axis(1,1)) then
            fv(isec)%u(1:5,iele)=  u_axis(2:6,1)
        elseif(a .ge. u_axis(1,n_axis)) then
            fv(isec)%u(1:5,iele)=  u_axis(2:6,n_axis)
        else
            do i=1,n_axis-1
                if((a .ge. u_axis(1,i)) .and. (a .le. u_axis(1,i+1))) then
                    r   = (a-u_axis(1,i))/(u_axis(1,i+1)-u_axis(1,i))
                    fv(isec)%u(1:5,iele)=  u_axis(2:6,i  )*(1.0d0-r) &
                                        & +u_axis(2:6,i+1)*r
                    exit
                end if
            end do
        end if
    end do
    end do

    return
    end subroutine fv_prefld_along_axis
!-------------------------------------------------------------------------------
!   solution initializtion with existing cell-centered solution.
!-------------------------------------------------------------------------------
    subroutine fv_prefld_cc(file_name)
    use var_kind_def
    use var_cgns
    use var_fv
    use var_global, only: is_show_res,err_mem
    use var_load_balance
    use var_mesh
    use var_parallel
    use var_turb
    implicit none
    character(len=*),intent(in):: file_name
    logical(dpL):: ltmp
    integer(dpI):: isec,iele,e1,e0,nele,iseg,ID,M,LDA

    if(myid .eq. 0) then
        inquire(file=trim(adjustl(file_name)), exist=ltmp)
        if(.not. ltmp)  stop 'Error: cell-centered solution file not found.'

        call cg_open_f(trim(adjustl(file_name)), MODE_READ, ifile, cg_err)
    end if
    ibase   =  1
    izone   =  1
    iflow   =  1

    nele=  min(2**20, n_ele_g)
    LDA =  5
    if(is_RANS) then
        LDA =  6
        if(transition_model .ge. 1) LDA =  7
    end if

    if(.not. allocated(r_recv)) then
        allocate(r_recv(nele*LDA), stat=err_mem)
    else
        if(size(r_recv) .lt. nele*LDA) then
            deallocate(r_recv)
            allocate(r_recv(nele*LDA), stat=err_mem)
        end if
    end if
    do iseg=1,(n_elei_g-1)/nele+1
        e1  =  1+nele*(iseg-1)
        e0  =  min(iseg*nele, n_elei_g)
        M   =  e0-e1+1

        if(myid .eq. 0) then
            call cg_field_read_f(ifile, ibase, izone, 1, 'Density', RealDouble, &
                &  e1, e0, r_recv(    1), cg_err)
            if(cg_err .ne. 0)   stop 'Error: fails to read cell-centered Density.'

            call cg_field_read_f(ifile, ibase, izone, 1, 'VelocityX', RealDouble, &
                &  e1, e0, r_recv(  M+1), cg_err)
            if(cg_err .ne. 0)   stop 'Error: fails to read cell-centered VelocityX.'

            call cg_field_read_f(ifile, ibase, izone, 1, 'VelocityY', RealDouble, &
                &  e1, e0, r_recv(2*M+1), cg_err)
            if(cg_err .ne. 0)   stop 'Error: fails to read cell-centered VelocityY.'

            call cg_field_read_f(ifile, ibase, izone, 1, 'VelocityZ', RealDouble, &
                &  e1, e0, r_recv(3*M+1), cg_err)
            if(cg_err .ne. 0)   stop 'Error: fails to read cell-centered VelocityZ.'

            call cg_field_read_f(ifile, ibase, izone, 1, 'Pressure', RealDouble, &
                &  e1, e0, r_recv(4*M+1), cg_err)
            if(cg_err .ne. 0)   stop 'Error: fails to read cell-centered Pressure.'

            if(is_RANS) then
                call cg_field_read_f(ifile, ibase, izone, 1, 'TurbulentSANuTilde', &
                    &  RealDouble, e1, e0, r_recv(5*M+1), cg_err)
                if(cg_err .ne. 0)   stop 'Error: fails to read cell-centered SANuTilde.'

                if(transition_model .ge. 1) then
                    call cg_field_read_f(ifile, ibase, izone, 1, 'Intermittency', &
                        &  RealDouble, e1, e0, r_recv(6*M+1), cg_err)
                    if(cg_err .ne. 0)   r_recv(6*M+1:7*M)   =  1.0d0
                end if
            end if
        end if

        call mpi_bcast(r_recv, LDA*M, mpi_dpR, 0, mpi_comm_world, mpi_err)

        do isec=mesh(0)%sec_1,mesh(0)%sec_0
            if(.not. sec(isec)%is_int)  cycle

            do iele=1,sec(isec)%n_ele
                ID  =  sec(isec)%id_ele_i(iele)
                if((ID .lt. e1) .or. (ID .gt. e0))  cycle
                ID  =  ID-e1+1

                fv(isec)%u(1,iele)  =  r_recv(ID    )
                fv(isec)%u(2,iele)  =  r_recv(ID+  M)
                fv(isec)%u(3,iele)  =  r_recv(ID+2*M)
                fv(isec)%u(4,iele)  =  r_recv(ID+3*M)
                fv(isec)%u(5,iele)  =  r_recv(ID+4*M)

                if(is_RANS) then
                    fv(isec)%turb(1,iele)   =  r_recv(ID+5*M)
                    if((transition_model .eq. 1) .or. (transition_model .eq. 2)) then
                        fv(isec)%turb(2,iele)   =  r_recv(ID+6*M)
                    elseif(transition_model .gt. 2) then
                        fv(isec)%intermittency(iele)=  r_recv(ID+6*M)
                    end if
                end if
            end do
        end do
    end do

    if(myid .eq. 0) then
        call cg_close_f(ifile, cg_err)
        if(is_show_res) then
            print*,'----------------------------------------------------------'
            print*,'Solution initialized from file: ',trim(adjustl(file_name))
        end if
    end if

    return
    end subroutine fv_prefld_cc
!-------------------------------------------------------------------------------
!   time marching solution.
!-------------------------------------------------------------------------------
    subroutine fv_solution
    use var_slv
    use var_turb, only: is_RANS,transition_model,transition_AGS
    implicit none

    call fv_set_slv
    call fv_set_boundary_condition
    call fv_prefld
    if(transition_model .eq. transition_AGS)    call mesh_get_1d_lines

    lev_out =  0
    do ite_wrk=1,max_iter_lev(lev_out)
        call fv_gmg_cycle(lev_out)
        if(is_RANS) call fv_rans_slv(lev_out)
        call fv_slv_monitor(.false., lev_out)

        if(is_converged)    exit
    end do

    call fv_unsteady_cal
    call fv_slv_monitor(.true., 0)

!   call fv_get_entropy_generation

    return
    end subroutine fv_solution

    subroutine fv_struct_to_array(lev)
    use var_kind_def
    use var_slv
    use var_global
    use var_fv
    use var_fv_array
    use var_mesh
    use var_parallel
    use var_turb
    implicit none

    integer(dpI):: isec,iele,n_ele,ID
    integer(dpI),intent(in):: lev

    if(.not. allocated(fv_order)) allocate(fv_order(tot_ele))
    if(.not. allocated(fv_LHS_s)) allocate(fv_LHS_s(tot_ele))
    if(.not. allocated(fv_rv0))   allocate(fv_rv0(5,tot_ele))
    do isec = mesh(lev)%sec_1,mesh(lev)%sec_0
        do iele = 1,sec(isec)%n_ele
            ID = sec(isec)%ID_ele_g(iele)
            if(ID .le. cellNum) ID = perm(ID)
            fv_u       (:,ID) = fv(isec)%u    (:,iele)
            fv_t       (  ID) = fv(isec)%t    (  iele)
            fv_uc      (:,ID) = fv(isec)%uc   (:,iele)
            fv_uc0     (:,ID) = fv(isec)%uc0  (:,iele)
            fv_duc     (:,ID) = fv(isec)%duc  (:,iele)
            fv_rhs     (:,ID) = fv(isec)%rhs  (:,iele)
            if(allocated(fv(isec)%LHS_s)) then
                fv_LHS_s(ID) = fv(isec)%LHS_s(iele)
            end if
            fv_mu      (:,ID) = fv(isec)%mu   (:,iele)
            fv_gra     (:,ID) = fv(isec)%gra  (:,iele)
            if(is_RANS) fv_turb    (:,ID) = fv(isec)%turb (:,iele)
            if(is_RANS) fv_turb_gra(:,ID) = fv(isec)%turb_gra(:,iele)

            fv_order   (  ID) = fv(isec)%order

            fv_uL          (:,ID) = fv(isec)%uL          (:,iele)
            fv_uR          (:,ID) = fv(isec)%uR          (:,iele)
            fv_bnd_solution(:,ID) = fv(isec)%bnd_solution(:,iele)

            if(allocated(fv(isec)%rv0)) then
                fv_rv0(:,ID) = fv(isec)%rv0(:,iele)
            end if
        end do
    end do

    return
    end subroutine fv_struct_to_array


    subroutine fv_array_to_struct(lev)
    use var_slv
    use var_global
    use var_fv
    use var_fv_array
    use var_mesh
    use var_turb
    implicit none

    integer(dpI):: isec,iele,n_ele,ID
    integer(dpI),intent(in):: lev

    if(.not. allocated(fv_order)) allocate(fv_order(tot_ele))
    if(.not. allocated(fv_LHS_s)) allocate(fv_LHS_s(tot_ele))
    if(.not. allocated(fv_rv0))   allocate(fv_rv0(5,tot_ele))
    do isec = mesh(lev)%sec_1,mesh(lev)%sec_0
        do iele = 1,sec(isec)%n_ele
            ID = sec(isec)%ID_ele_g(iele)
            if(ID .le. cellNum) ID = perm(ID)
            fv(isec)%u       (:,iele) = fv_u    (:,ID) 
            fv(isec)%t       (  iele) = fv_t    (  ID) 
            fv(isec)%uc      (:,iele) = fv_uc   (:,ID) 
            fv(isec)%uc0     (:,iele) = fv_uc0  (:,ID) 
            fv(isec)%duc     (:,iele) = fv_duc  (:,ID) 
            fv(isec)%rhs     (:,iele) = fv_rhs  (:,ID) 
            if(allocated(fv(isec)%LHS_s)) then
                fv(isec)%LHS_s(iele) = fv_LHS_s(ID) 
            end if
            fv(isec)%mu      (:,iele) = fv_mu   (:,ID) 
            fv(isec)%gra     (:,iele) = fv_gra  (:,ID) 
            if(is_RANS) fv(isec)%turb    (:,iele) = fv_turb (:,ID) 
            if(is_RANS) fv(isec)%turb_gra(:,iele) = fv_turb_gra(:,ID) 

            fv(isec)%order            = fv_order(ID)

            fv(isec)%uL          (:,iele) = fv_uL          (:,ID)
            fv(isec)%uR          (:,iele) = fv_uR          (:,ID)
            fv(isec)%bnd_solution(:,iele) = fv_bnd_solution(:,ID)
            if(allocated(fv(isec)%rv0)) then
                fv(isec)%rv0(:,iele) = fv_rv0(:,ID)
            end if
        end do
    end do

    return
    end subroutine fv_array_to_struct

    subroutine check_struct_to_array(lev)
    use var_slv
    use var_global
    use var_fv
    use var_fv_array
    use var_mesh
    use var_parallel
    implicit none

    integer(dpI):: isec,iele,n_ele,i
    integer(dpI),intent(in):: lev

    n_ele = 0
    do isec = mesh(lev)%sec_1,mesh(lev)%sec_0
!    if(myid .ne. 34) cycle
        do iele = 1,sec(isec)%n_ele
!            write(*,*),isec,iele,sec(isec)%ID_ele_g(iele)
            do i = 1,5
                if(fv_u(i,sec(isec)%ID_ele_g(iele)) &
                            & .ne. fv(isec)%u(i,iele)) then
                    write(*,*),fv_u(i,sec(isec)%ID_ele_g(iele)),&
                                & fv(isec)%u(i,iele)
                end if
            end do
        end do
    end do

    return
    end subroutine check_struct_to_array

    subroutine sec_struct_to_array(lev)
    use var_kind_def
    ! use var_slv
    use var_global
    ! use var_fv
    ! use var_fv_array
    use var_mesh
    use var_sec_array
    use var_parallel
    implicit none

    integer(dpI):: isec,iele,n_ele,ID
    integer(dpI),intent(in):: lev

    if(.not. allocated(sec_cen))      allocate(sec_cen(3, tot_ele))
    if(.not. allocated(sec_vol))      allocate(sec_vol(   tot_ele))
    if(.not. allocated(sec_is_int))   allocate(sec_is_int(tot_ele))
    if(.not. allocated(sec_is_ghost)) allocate(sec_is_ghost(tot_ele))
    if(.not. allocated(sec_is_bnd))   allocate(sec_is_bnd(tot_ele))
    if(.not. allocated(sec_bct))      allocate(sec_bct(tot_ele))
    do isec = mesh(lev)%sec_1,mesh(lev)%sec_0
        do iele = 1,sec(isec)%n_ele
            ID = sec(isec)%ID_ele_g(iele)
            if(ID .le. cellNum) ID = perm(ID)
            ! if(myid==3) write(*,*),isec,iele,ID,size(sec_cen)
            sec_cen(:,ID) = sec(isec)%cen(:,iele)
            sec_vol(  ID) = sec(isec)%vol(  iele)
            sec_bct(  ID) = sec(isec)%bct
            if(sec(isec)%is_int) then
                sec_is_int(ID) = 1
            else
                sec_is_int(ID) = -1
            end if
            if(sec(isec)%is_ghost) then
                sec_is_ghost(ID) = 1
            else
                sec_is_ghost(ID) = -1
            end if
            if(sec(isec)%is_bnd) then
                sec_is_bnd(ID) = 1
            else
                sec_is_bnd(ID) = -1
            end if
        end do
    end do

    return
    end subroutine sec_struct_to_array

    subroutine transform_parameter_to_real(lev)
        use var_kind_def
        use var_global_real
        use var_fv
        use var_turb
        use var_slv, only: is_vis_cal
        implicit none
        integer(dpI),intent(IN):: lev
        
        if(is_limiter_on) then
            is_limiter_on_r = 1
        else
            is_limiter_on_r = -1
        end if

        if(is_tur_cal) then
            is_tur_cal_r = 1
        else
            is_tur_cal_r = -1
        end if

        if(is_vis_cal) then
            is_vis_cal_r = 1
        else
            is_vis_cal_r = -1
        end if

        if(is_KO) then
            is_KO_r = 1
        else
            is_KO_r = -1
        end if

        if(is_RANS) then
            is_RANS_r = 1
        else
            is_RANS_r = -1
        end if

        RANS_model_r = real(RANS_model, dpR)
        rhs_lev_r(lev) = rhs_lev(lev)
    return
    end subroutine transform_parameter_to_real