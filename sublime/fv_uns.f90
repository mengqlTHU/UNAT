!-------------------------------------------------------------------------------
!   Unsteady time marching, read data from existing file.
!-------------------------------------------------------------------------------
    subroutine fv_unsteady_prefld
    use var_kind_def
    use var_fv
    use var_global, only: err_mem
    use var_mesh
    use var_parallel
    use var_turb
    use var_uns_cal
    implicit none
    logical(dpL):: ltm1,ltm2,ltm3,alive(3)
    integer(dpI):: isec,isol,iele

    if(is_DDES .and. (uns_ra_begin .le. max_uns_iter)) then
        do isec=mesh(0)%sec_1,mesh(0)%sec_0
            if(.not. sec(isec)%is_int)  cycle
            allocate(fv(isec)%DDES_quality(n_DDES_quality,sec(isec)%n_ele), stat=err_mem)
            fv(isec)%DDES_quality   =  0.0d0
        end do
    end if

    ltm1=  .false.
    ltm2=  .false.
    ltm3=  .false.
    if(uns_solution_1(1:4) .ne. 'NONE') &
        &  inquire(file=trim(adjustl(uns_solution_1)), exist=ltm1)
    if(uns_solution_2(1:4) .ne. 'NONE') &
        &  inquire(file=trim(adjustl(uns_solution_2)), exist=ltm2)
    if(uns_solution_3(1:4) .ne. 'NONE') &
        &  inquire(file=trim(adjustl(uns_solution_3)), exist=ltm3)
    alive(1:3)  = (/ltm1, ltm2, ltm3/)
    call mpi_bcast(alive, 3, mpi_dpL, 0, mpi_comm_world, mpi_err)

    if((uns_method .eq. 1) .or. (uns_method .eq. 4)) then
        ltm1    =  alive(1)
        if(ltm1)    uns_n_ff=  1
    elseif(uns_method .eq. 2) then
        if(alive(1) .and. alive(2)) then
            ltm1    =  .true.
            uns_n_ff=  2
        elseif(alive(1)) then
            ltm1    =  .true.
            uns_n_ff=  1
        end if
    elseif(uns_method .eq. 3) then
        if(alive(1) .and. alive(2) .and. alive(3)) then
            ltm1    =  .true.
            uns_n_ff=  3
        elseif(alive(1) .and. alive(2)) then
            ltm1    =  .true.
            uns_n_ff=  2
        elseif(alive(1)) then
            ltm1    =  .true.
            uns_n_ff=  1
        end if
    else
        stop 'Error: uns_method not supported.'
    end if
    if(.not. ltm1)  return

    do isol=uns_n_ff,1,-1
        if(isol .eq. 3) then
!           read solution for u^{n-3}.
            call fv_prefld_cc(uns_solution_3)
            call fv_lev_u2uc(0, 0)
            do isec=mesh(0)%sec_1,mesh(0)%sec_0
                if(.not. sec(isec)%is_int)  cycle
                call DCOPY(5*sec(isec)%n_ele, fv(isec)%uc, 1, fv(isec)%uns_uc_2, 1)
            end do
        elseif(isol .eq. 2) then
!           read solution for u^{n-2}.
            call fv_prefld_cc(uns_solution_2)
            call fv_lev_u2uc(0, 0)
            do isec=mesh(0)%sec_1,mesh(0)%sec_0
                if(.not. sec(isec)%is_int)  cycle
                call DCOPY(5*sec(isec)%n_ele, fv(isec)%uc, 1, fv(isec)%uns_uc_1, 1)
                if(RANS_model .eq. RANS_SA) then
                    call DCOPY(sec(isec)%n_ele, fv(isec)%turb,1, fv(isec)%uns_turb(1,1),2)
                elseif(RANS_model .eq. RANS_WA) then
                    call DCOPY(sec(isec)%n_ele, fv(isec)%turb,1, fv(isec)%uns_turb(1,1),2)
                elseif(RANS_model .eq. RANS_SST) then
                    do iele=1,sec(isec)%n_ele
                        fv(isec)%uns_turb(1:2,iele) =  fv(isec)%u(1,iele) &
                            & *fv(isec)%turb(1:2,iele)
                    end do
                end if
            end do
        elseif(isol .eq. 1) then
!           read solution for u^{n-1}.
            call fv_prefld_cc(uns_solution_1)
            call fv_lev_u2uc(0, 0)
            call fv_bnd_parallel(0, .false.)
            if(is_RANS) call fv_rans_bnd_parallel(0, .true.)
        end if
    end do
    is_uns_initialized  =  .true.
    if(myid .eq. 0) then
        print*,'----------------------------------------------------------'
        print*,'Unsteady solution initialized from files.'
    end if

    return
    end subroutine fv_unsteady_prefld
!-------------------------------------------------------------------------------
!   add the source term, BDF unsteady.
!-------------------------------------------------------------------------------
    subroutine fv_get_src_bdf(lev)
    use var_kind_def
    use var_fv
    use var_mesh
    use var_uns_cal
    implicit none
    integer(dpI),intent(in):: lev
    integer(dpI):: isec,iele,method
    real   (dpR):: dt

    if(.not. is_bdf_now)    return

    dt  =  1.0d0/dt_uns
    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(.not. sec(isec)%is_int)  cycle

        if(uns_n_ff+uns_iter .ge. uns_n_stencil) then
            method  =  uns_method
        elseif(uns_n_ff+uns_iter .eq. 3) then
            method  =  2
        else
            method  =  1
        end if

        if(method .eq. 1) then
            do iele=1,sec(isec)%n_ele
                fv(isec)%rhs(1:5,iele)  =  fv(isec)%rhs(1:5,iele) &
                    &+(fv(isec)%uc      (1:5,iele) &
                    &- fv(isec)%uns_uc_1(1:5,iele))*sec(isec)%vol(iele)*dt
            end do
        elseif(method .eq. 2) then
            do iele=1,sec(isec)%n_ele
                fv(isec)%rhs(1:5,iele)  =  fv(isec)%rhs(1:5,iele) &
                    &+(1.5d0*fv(isec)%uc      (1:5,iele) &
                    &- 2.0d0*fv(isec)%uns_uc_1(1:5,iele) &
                    &+ 0.5d0*fv(isec)%uns_uc_2(1:5,iele))*sec(isec)%vol(iele)*dt
            end do
        elseif(method .eq. 3) then
            do iele=1,sec(isec)%n_ele
                fv(isec)%rhs(1:5,iele)  =  fv(isec)%rhs(1:5,iele) &
                    &+(BDF2opt_c1*fv(isec)%uc      (1:5,iele) &
                    &+ BDF2opt_c2*fv(isec)%uns_uc_1(1:5,iele) &
                    &+ BDF2opt_c3*fv(isec)%uns_uc_2(1:5,iele) &
                    &+ BDF2opt_c4*fv(isec)%uns_uc_3(1:5,iele))*sec(isec)%vol(iele)*dt
            end do
        else
            stop 'Error: uns_method not supported.'
        end if
    end do

    return
    end subroutine fv_get_src_bdf
!-------------------------------------------------------------------------------
!   get the running-average and Reynolds stress.
!-------------------------------------------------------------------------------
    subroutine fv_get_uns_ra(lev)
    use var_kind_def
    use var_fv
    use var_global, only: is_output_cc,translation_axis,rotation_axis,mesh_name
    use var_mesh, only: sec,mesh,sec_motion_type,sec_motion_speed
    use var_turb, only: is_RANS,is_DDES_now
    use var_uns_cal
    implicit none
    integer(dpI),intent(in):: lev
    integer(dpI):: isec,iele,ID,M
    real   (dpR):: eps,u(6),ome(3),c(3),u_old(6),s_new(7),s_old(7)

    if((.not. is_uns_ra) .or. (uns_iter .lt. uns_ra_begin)) return

    if(is_RANS) then
        M   =  6
    else
        M   =  5
    end if
    eps =  1.0d0/real(1+uns_iter-uns_ra_begin, dpR)
    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
!       if(sec(isec)%is_ghost)  cycle

        ID  =  sec(isec)%ID_sec_g
        do iele=1,sec(isec)%n_ele
            u(1:5)  =  fv(isec)%u(1:5,iele)
            if(sec_motion_type(ID) .eq. 1) then
                u(2:4)  =  u(2:4)-sec_motion_speed(ID)*translation_axis(1:3)
            elseif(sec_motion_type(ID) .eq. 2) then
                ome     =  sec_motion_speed(ID)*rotation_axis
                c(1:3)  =  sec(isec)%cen(1:3,iele)
                u(2)    =  u(2)-(ome(2)*c(3)-ome(3)*c(2))
                u(3)    =  u(3)-(ome(3)*c(1)-ome(1)*c(3))
                u(4)    =  u(4)-(ome(1)*c(2)-ome(2)*c(1))
            end if
            if(is_RANS) u(6)=  fv(isec)%turb(1,iele)

            u_old(1:M)  =  fv(isec)%uns_ra(1:M,iele)
            fv(isec)%uns_ra(1:M,iele)   =  fv(isec)%uns_ra(1:M,iele) &
                & +(u(1:M)-u_old(1:M))*eps
            s_old(1:7)  =  fv(isec)%stress(1:7,iele)
            s_new(1:3)  =  u(2)*u(2:4)
            s_new(4:5)  =  u(3)*u(3:4)
            s_new(6)    =  u(4)*u(4)
            s_new(7)    =  u(5)*u(5)
            fv(isec)%stress(1:7,iele)   =  fv(isec)%stress(1:7,iele) &
                & +(s_new(1:7)-s_old(1:7))*eps

            if(uns_iter .eq. max_uns_iter) then
!               compute the Reynolds stress and RMSpressure.
                s_old(1:7)  =  fv(isec)%stress(1:7,iele)
                u    (1:5)  =  fv(isec)%uns_ra(1:5,iele)
                fv(isec)%stress(1,iele) =  s_old(1)-u(2)*u(2)
                fv(isec)%stress(2,iele) =  s_old(2)-u(2)*u(3)
                fv(isec)%stress(3,iele) =  s_old(3)-u(2)*u(4)
                fv(isec)%stress(4,iele) =  s_old(4)-u(3)*u(3)
                fv(isec)%stress(5,iele) =  s_old(5)-u(3)*u(4)
                fv(isec)%stress(6,iele) =  s_old(6)-u(4)*u(4)
                fv(isec)%stress(7,iele) =  sqrt(abs(s_old(7)-u(5)*u(5)))
            end if
        end do
    end do
    if(uns_iter .eq. max_uns_iter) then
        call fv_wr_cgns(2, trim(adjustl(mesh_name))//'_uns_ra_vc.cgns')
        if(is_output_cc)    call fv_wr_cgns_cc(1, trim(adjustl(mesh_name))//'_uns_ra_cc.cgns')
    end if

    if(is_DDES_now) call fv_get_DDES_quality(lev)

    return
    end subroutine fv_get_uns_ra
!-------------------------------------------------------------------------------
!   Unsteady time marching.
!-------------------------------------------------------------------------------
    subroutine fv_unsteady_cal
    use var_kind_def
    use var_fv
    use var_global, only: err_mem
    use var_mesh
    use var_parallel
    use var_slv
    use var_turb
    use var_uns_cal
    implicit none
    integer(dpI):: isec,i
    real   (dpR):: dt,s

    if(.not. is_uns_cal)    return
    is_st_cal_now       =  .false.
    is_uns_cal_now      =  .true.
    is_local_dt         =  .false.
    lev_out             =  0
    is_bdf_now          = (uns_method .ge. 2) .and. (uns_method .le. 3)
    is_DDES_now         =  is_DDES
    is_LES_now          =  is_LES
    if(is_RANS .and. (.not. is_BDF_now))    stop 'Error: SA not supported by UNS_METHOD.'
    if(myid .eq. 0) then
        allocate(uns_record(max_uns_record,10+5*n_monitor  ), stat=err_mem)
    else
        allocate(uns_record(max_uns_record,10+5*n_monitor_L), stat=err_mem)
    end if
    uns_record  =  0.0d0
    if(uns_method .eq. uns_BDF2opt) LHS_uns_c1  =  BDF2opt_c1

!   ----------------------------------------------------------------------------
!   check the timestep.
    call fv_cal_LHS_scalar(0, .false., uns_method .eq. uns_IMEX)
    dt  =  huge(1.0d0)
    do isec=mesh(0)%sec_1,mesh(0)%sec_0
        if(sec(isec)%is_ghost .or. sec(isec)%is_bnd)    cycle
        do i=1,sec(isec)%n_ele
            dt  =  min(dt, sec(isec)%vol(i)/fv(isec)%LHS_s(i))
            fv(isec)%LHS_s(i)   =  dt_uns/sec(isec)%vol(i)
        end do

        allocate(fv(isec)%uns_uc_1(5,sec(isec)%n_ele))
        if(is_BDF_now)  allocate(fv(isec)%uns_uc_2(5,sec(isec)%n_ele))
        if(is_BDF2opt)  allocate(fv(isec)%uns_uc_3(5,sec(isec)%n_ele))
        if(is_RANS   ) then
            if(RANS_model .lt. RANS_SST) then
                allocate(fv(isec)%uns_turb(2,sec(isec)%n_ele))
            else
                allocate(fv(isec)%uns_turb(4,sec(isec)%n_ele))
            end if
        end if
    end do
    if(is_uns_ra) then
        do isec=mesh(0)%sec_1,mesh(0)%sec_0
            if(is_RANS) then
                allocate(fv(isec)%uns_ra(6,sec(isec)%n_ele), stat=err_mem)
            else
                allocate(fv(isec)%uns_ra(5,sec(isec)%n_ele), stat=err_mem)
            end if
            allocate(fv(isec)%stress(7,sec(isec)%n_ele), stat=err_mem)
            fv(isec)%uns_ra =  0.0d0
            fv(isec)%stress =  0.0d0
        end do
    end if
    s   =  dt
    call mpi_allreduce(s, dt, 1, mpi_dpR, mpi_min, mpi_comm_world, mpi_err)
    if(myid .eq. 0) then
        print*,'----------------------------------------------------------'
        write(unit=6,fmt='(A17,ES20.12)'),'CFL=1D0, dt_uns=',dt
    end if
!   check the timestep.
!   ----------------------------------------------------------------------------

    call fv_unsteady_prefld

    do uns_iter=1,max_uns_iter
        tph =  tph_1+real(uns_iter, dpR)*dt_uns

        if((is_DDES .or. is_LES) .and. (length_scale .eq. length_LSQ))  call fv_get_gls

        if(is_bdf_now) then
            do isec=mesh(0)%sec_1,mesh(0)%sec_0
                if(.not. sec(isec)%is_int)  cycle
                i   =  5*sec(isec)%n_ele
                if(is_BDF2opt)  call DCOPY(i, fv(isec)%uns_uc_2, 1, fv(isec)%uns_uc_3, 1)
                call DCOPY(i, fv(isec)%uns_uc_1, 1, fv(isec)%uns_uc_2, 1)
                call DCOPY(i, fv(isec)%uc      , 1, fv(isec)%uns_uc_1, 1)
                if(RANS_model .eq. RANS_SA) then
                    do i=1,sec(isec)%n_ele
                        fv(isec)%uns_turb(2,i)  =  fv(isec)%uns_turb(1,i)
                        fv(isec)%uns_turb(1,i)  =  fv(isec)%turb    (1,i)
                    end do
                elseif(RANS_model .eq. RANS_WA) then
                    do i=1,sec(isec)%n_ele
                        fv(isec)%uns_turb(2,i)  =  fv(isec)%uns_turb(1,i)
                        fv(isec)%uns_turb(1,i)  =  fv(isec)%turb    (1,i)
                    end do
                elseif(RANS_model .eq. RANS_SAM) then
                    do i=1,sec(isec)%n_ele
                        fv(isec)%uns_turb(3:4,i)=  fv(isec)%uns_turb(1:2,i)
                        fv(isec)%uns_turb(1  ,i)=  fv(isec)%turb(1,i)
                        fv(isec)%uns_turb(2  ,i)=  fv(isec)%turb(2,i)*fv(isec)%u(1,i)
                    end do
                elseif(RANS_model .eq. RANS_SST) then
                    do i=1,sec(isec)%n_ele
                        fv(isec)%uns_turb(3:4,i)=  fv(isec)%uns_turb(1:2,i)
                        fv(isec)%uns_turb(1:2,i)=  fv(isec)%turb(1:2,i)*fv(isec)%u(1,i)
                    end do
                end if
            end do

            is_get_DDES_quality = (uns_iter .ge. uns_ra_begin) .and. is_DDES_now
            do ite_wrk=1,max_subiter_uns
                call fv_gmg_cycle(0)
                if(is_RANS) call fv_rans_slv(0)
                call fv_slv_monitor(.false., 0)
                if(is_converged)    exit
            end do
            if(ite_wrk .gt. max_subiter_uns)    uns_fail=  uns_fail+1
        elseif(uns_method .eq. uns_IMEX) then
            call fv_slv_imex
            call fv_slv_monitor(.false., 0)
        else
            call fv_slv_exp_rk43_ssp(0)
            call fv_slv_monitor(.false., 0)
        end if

        call fv_uns_wr_results
    end do

    if((myid .eq. 0) .and. (uns_fail .gt. 0)) then
        print*,'----------------------------------------------------------'
        write(6,fmt='(A6,i5,A23)'),'Note:',uns_fail,' unsteady steps failed.'
    end if

    return
    end subroutine fv_unsteady_cal
