!-------------------------------------------------------------------------------
!   set the cfg_file.
!-------------------------------------------------------------------------------
    subroutine set_cfg_file(str)
    use var_global, only: cfg_file
    implicit none
    character(len=*),intent(in):: str

    cfg_file=  trim(adjustl(str))

    return
    end subroutine set_cfg_file
!-------------------------------------------------------------------------------
!   time marching solution interface for the HEAT.
!-------------------------------------------------------------------------------
    subroutine heat_solution_interface(residual)
    use var_heat, only: res_heat
    use var_slv
    implicit none
    real(kind=8),intent(out):: residual

    lev_out =  0
    ite_wrk =  ite_wrk+1
    call heat_solution_step(0)
    call heat_slv_monitor(.false., 0)
    residual=  res_heat

    return
    end subroutine heat_solution_interface
!-------------------------------------------------------------------------------
!   time marching solution interface for the FV.
!-------------------------------------------------------------------------------
    subroutine fv_solution_interface(residual)
    use var_fv, only: res_NS
    use var_slv
    use var_turb, only: is_RANS
    implicit none
    real(kind=8),intent(out):: residual
!f2py real(kind=8),intent(out):: residual

    lev_out =  0
    ite_wrk =  ite_wrk+1
    call fv_solution_step(0)
    if(is_RANS) call fv_rans_slv(0)
    call fv_slv_monitor(.false., 0)
    residual=  res_NS

    return
    end subroutine fv_solution_interface
!-------------------------------------------------------------------------------
!   time marching unsteady solution interface for the FV.
!-------------------------------------------------------------------------------
    subroutine fv_unsteady_interface(iter_outer,iter_inner,residual)
    use var_kind_def
    use var_fv
    use var_global, only: mesh_name,err_mem,is_show_res
    use var_mesh
    use var_parallel
    use var_slv, only: is_converged,is_local_dt,is_st_cal_now,lev_out,ite_wrk
    use var_turb, only: is_DDES,is_DDES_now,is_RANS
    use var_uns_cal
    implicit none
    integer(dpI),intent(in):: iter_outer,iter_inner
    character(len=80):: str
    integer(dpI):: isec,i
    real   (dpR):: residual,dt,s

    if(.not. is_uns_cal)    stop 'Error: unsteady solution is not defined.'
    if(uns_method .ne. 2)   stop 'Error: py_unsteady supports BDF only.'

    uns_iter=  iter_outer
    ite_wrk =  iter_inner

    if((uns_iter .le. 1) .and. (ite_wrk .le. 1)) then
    is_st_cal_now   =  .false.
    is_uns_cal_now  =  .true.
    is_local_dt     =  .false.
    lev_out         =  0
    is_bdf_now      =  uns_method .eq. 2
    is_DDES_now     =  is_DDES
    allocate(uns_record(10,max_uns_record), stat=err_mem)

!   ----------------------------------------------------------------------------
!   check the timestep.
    call fv_cal_LHS_scalar(0, .false.)
    dt  =  huge(1.0d0)
    do isec=mesh(0)%sec_1,mesh(0)%sec_0
        if(sec(isec)%is_ghost .or. sec(isec)%is_bnd)    cycle
        do i=1,sec(isec)%n_ele
            dt  =  min(dt, sec(isec)%vol(i)/fv(isec)%LHS_s(i))
            fv(isec)%LHS_s(i)   =  dt_uns/sec(isec)%vol(i)
        end do

        if(is_bdf_now) then
            allocate(fv(isec)%uns_uc_1(5,sec(isec)%n_ele))
            allocate(fv(isec)%uns_uc_2(5,sec(isec)%n_ele))
            if(is_RANS) allocate(fv(isec)%uns_turb(2,sec(isec)%n_ele))
        end if
    end do
    if(is_uns_ra) then
        do isec=mesh(0)%sec_1,mesh(0)%sec_0
            if(is_RANS) then
                allocate(fv(isec)%uns_ra(6,sec(isec)%n_ele))
            else
                allocate(fv(isec)%uns_ra(5,sec(isec)%n_ele))
            end if
            allocate(fv(isec)%stress(7,sec(isec)%n_ele))
            fv(isec)%uns_ra =  0.0d0
            fv(isec)%stress =  0.0d0
        end do
    end if
    s   =  dt
    call mpi_allreduce(s, dt, 1, mpi_dpR, mpi_min, mpi_comm_world, mpi_err)
    if((myid .eq. 0) .and. is_show_res) then
        print*,'----------------------------------------------------------'
        write(unit=6,fmt='(A17,ES20.12)'),'CFL=1D0, dt_uns=',dt
    end if
!   check the timestep.
!   ----------------------------------------------------------------------------
    end if

    if(ite_wrk .le. 1) then
        tph =  tph_1+real(uns_iter, dpR)*dt_uns
        do isec=mesh(0)%sec_1,mesh(0)%sec_0
            if(.not. sec(isec)%is_int)  cycle
            call DCOPY(5*sec(isec)%n_ele, fv(isec)%uns_uc_1, 1, fv(isec)%uns_uc_2, 1)
            call DCOPY(5*sec(isec)%n_ele, fv(isec)%uc      , 1, fv(isec)%uns_uc_1, 1)
            if(is_RANS) then
                do i=1,sec(isec)%n_ele
                    fv(isec)%uns_turb(2,i)  =  fv(isec)%uns_turb(1,i)
                    fv(isec)%uns_turb(1,i)  =  fv(isec)%turb    (  i)
                end do
            end if
        end do
    end if

    call fv_solution_step(0)
    if(is_RANS) call fv_rans_slv(0)
    call fv_slv_monitor(.false., 0)
    residual=  res_NS

    if(ite_wrk .eq. max_subiter_uns) then
        call fv_wr_uns_record
        call fv_get_uns_ra(0)
    end if

    if((uns_iter .ge. max_uns_iter-1) .and. (ite_wrk .eq. max_subiter_uns)) then
!       save the solution fields at last two timesteps.
        write(str,*),uns_iter
        str =  trim(adjustl(mesh_name))//'_uns_'//trim(adjustl(str))//'.cgns'
        call fv_wr_cgns_cc(2, str)
    end if

    return
    end subroutine fv_unsteady_interface
!-------------------------------------------------------------------------------
!   time marching unsteady solution interface for the HEAT.
!-------------------------------------------------------------------------------
    subroutine heat_unsteady_interface(iter_outer,iter_inner,residual)
    use var_kind_def
    use var_global, only: is_show_res,err_mem
    use var_heat
    use var_mesh
    use var_parallel
    use var_slv, only: is_converged,is_local_dt,is_st_cal_now,lev_out,ite_wrk
    use var_uns_cal
    implicit none
    integer(dpI),intent(in):: iter_outer,iter_inner
    integer(dpI):: isec,i
    real   (dpR):: residual,dt,s

    if(.not. is_uns_cal)    stop 'Error: unsteady solution is not defined.'
    if(uns_method .ne. 2)   stop 'Error: py_unsteady supports BDF only.'
    uns_iter=  iter_outer
    ite_wrk =  iter_inner

    if((uns_iter .le. 1) .and. (ite_wrk .le. 1)) then
    is_st_cal_now   =  .false.
    is_uns_cal_now  =  .true.
    is_local_dt     =  .false.
    lev_out         =  0
    is_bdf_now      =  uns_method .eq. 2
    allocate(uns_record(10,max_uns_record), stat=err_mem)

!   ----------------------------------------------------------------------------
!   check the timestep.
    call heat_cal_LHS_scalar(0, .false.)
    dt  =  huge(1.0d0)
    do isec=mesh(0)%sec_1,mesh(0)%sec_0
        if(sec(isec)%is_ghost .or. sec(isec)%is_bnd)    cycle
        do i=1,sec(isec)%n_ele
            dt  =  min(dt, sec(isec)%vol(i)/heat(isec)%LHS_s(i))
            heat(isec)%LHS_s(i) =  dt_uns/sec(isec)%vol(i)
        end do

        if(is_bdf_now) then
            allocate(heat(isec)%uns_t_1(sec(isec)%n_ele))
            allocate(heat(isec)%uns_t_2(sec(isec)%n_ele))
        end if
    end do
    s   =  dt
    call mpi_allreduce(s, dt, 1, mpi_dpR, mpi_min, mpi_comm_world, mpi_err)
    if((myid .eq. 0) .and. is_show_res) then
        print*,'----------------------------------------------------------'
        write(unit=6,fmt='(A17,ES20.12)'),'CFL=1D0, dt_uns=',dt
    end if
!   check the timestep.
!   ----------------------------------------------------------------------------
    end if

    if(ite_wrk .le. 1) then
        tph =  tph_1+real(uns_iter, dpR)*dt_uns
        do isec=mesh(0)%sec_1,mesh(0)%sec_0
            if(.not. sec(isec)%is_int)  cycle
            call DCOPY(sec(isec)%n_ele, heat(isec)%uns_t_1, 1, heat(isec)%uns_t_2, 1)
            call DCOPY(sec(isec)%n_ele, heat(isec)%t      , 1, heat(isec)%uns_t_1, 1)
        end do
    end if
    call heat_solution_step(lev_out)
    call heat_slv_monitor(.false., lev_out)
    residual=  res_heat

    return
    end subroutine heat_unsteady_interface
