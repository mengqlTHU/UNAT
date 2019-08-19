    include './include/var_def.f90'

    program main
    use var_global, only: mbt_t1,mbt_t0,sw_slave
    use var_parallel, only: myid,mpi_err
    use var_slv, only: is_FV,is_heat
    implicit none

if(sw_slave) then
    call CG_init()
end if

    call setup(.true.)

if(sw_slave) then
    ! call swlu_debug_init
!    call swlu_prof_init
!    call swlu_prof_start
end if

    call CPU_time(mbt_t1)

    if(is_FV) then
        call fv_solution
    elseif(is_heat) then
    else
    end if

    call CPU_time(mbt_t0)
    if(myid .eq. 0) print*,'----------------------------------------------------------'
    if(myid .eq. 0) write(unit=6,fmt='(A14,ES20.12)'),'Elapsed time=',mbt_t0-mbt_t1

if(sw_slave) then
!    call swlu_prof_print
end if

if(sw_slave) then
    call CG_halt()
end if
    
    call getmaxtime()

    if(myid .eq. 0) call printtime()

    call mpi_finalize(mpi_err)
    stop
    end program main
