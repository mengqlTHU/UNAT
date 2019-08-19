!-------------------------------------------------------------------------------
!   py call to heat pre.
!-------------------------------------------------------------------------------
    subroutine py_heat_pre(cfg_file,mpi_started)
    implicit none
    character(len=*),intent(in):: cfg_file
    integer(kind=4),intent(inout):: mpi_started
!f2py intent(in):: cfg_file
!f2py intent(in,out):: mpi_started

    call set_cfg_file(cfg_file)
    call setup(mpi_started .ne. 1)
    mpi_started =  1
    call heat_set_slv

    return
    end subroutine py_heat_pre
!-------------------------------------------------------------------------------
!   py call to heat prefld.
!-------------------------------------------------------------------------------
    subroutine py_heat_prefld
    implicit none

    call heat_prefld

    return
    end subroutine py_heat_prefld
!-------------------------------------------------------------------------------
!   py call to heat solution.
!-------------------------------------------------------------------------------
    subroutine py_heat_solution(residual)
    implicit none
    real(kind=8),intent(out):: residual
!f2py real(kind=8),intent(out):: residual

    call heat_solution_interface(residual)

    return
    end subroutine py_heat_solution
!-------------------------------------------------------------------------------
!   py call to FV pre.
!-------------------------------------------------------------------------------
    subroutine py_fv_pre(cfg_file,mpi_started)
    implicit none
    character(len=*),intent(in):: cfg_file
    integer(kind=4),intent(inout):: mpi_started
!f2py intent(in):: cfg_file
!f2py intent(in,out):: mpi_started

    call set_cfg_file(cfg_file)
    call setup(mpi_started .ne. 1)
    mpi_started =  1
    call fv_set_slv

    return
    end subroutine py_fv_pre
!-------------------------------------------------------------------------------
!   py call to FV prefld.
!-------------------------------------------------------------------------------
    subroutine py_fv_prefld
    implicit none

    call fv_prefld

    return
    end subroutine py_fv_prefld
!-------------------------------------------------------------------------------
!   py call to FV solution.
!-------------------------------------------------------------------------------
    subroutine py_fv_solution(residual)
    implicit none
    real(kind=8),intent(out):: residual
!f2py real(kind=8),intent(out):: residual

    call fv_solution_interface(residual)

    return
    end subroutine py_fv_solution
!-------------------------------------------------------------------------------
!   get the number of cht sections.
!-------------------------------------------------------------------------------
    subroutine py_cht_get_n_sections(nsec)
    implicit none
    integer(kind=4):: nsec
!f2py intent(out):: nsec

    call cht_get_n_sections(nsec)

    return
    end subroutine py_cht_get_n_sections
!-------------------------------------------------------------------------------
!   get the information of a cht section.
!-------------------------------------------------------------------------------
    subroutine py_cht_get_section(isec,info,n2e,xyz)
    implicit none
    integer(kind=4),intent(in):: isec
    integer(kind=4):: info(*),n2e(*)
    real   (kind=8):: xyz(*)
!f2py intent(in):: isec
!f2py intent(in,out):: info,n2e(*)
!f2py intent(in,out):: xyz(*)

    call cht_get_section(isec, info, n2e, xyz)

    return
    end subroutine py_cht_get_section
!-------------------------------------------------------------------------------
!   set the ID of a cht section.
!-------------------------------------------------------------------------------
    subroutine py_cht_set_ID(isec,ID_coupler,ID_sec)
    implicit none
    integer(kind=4),intent(in):: isec,ID_coupler,ID_sec
!f2py intent(in):: isec,ID_coupler,ID_sec

    call cht_set_ID(isec, ID_coupler, ID_sec)

    return
    end subroutine py_cht_set_ID
!-------------------------------------------------------------------------------
!   send the geo and kappa to the coupler.
!-------------------------------------------------------------------------------
    subroutine py_cht_send_geo(idx,rbuf)
    implicit none
    integer(kind=4):: idx
    real   (kind=8):: rbuf(*)
!f2py intent(in,out):: idx
!f2py intent(in,out):: rbuf(*)

    call cht_send_geo(idx, rbuf)

    return
    end subroutine py_cht_send_geo
!-------------------------------------------------------------------------------
!   recv the geo and kappa from the coupler.
!-------------------------------------------------------------------------------
    subroutine py_cht_recv_geo(size_rbuf,rbuf)
    implicit none
    integer(kind=4),intent(in):: size_rbuf
    real   (kind=8),intent(in):: rbuf(*)
!f2py intent(in):: size_rbuf
!f2py intent(in):: rbuf(*)

    call cht_recv_geo(size_rbuf, rbuf)

    return
    end subroutine py_cht_recv_geo
!-------------------------------------------------------------------------------
!   send the temperature to the coupler.
!-------------------------------------------------------------------------------
    subroutine py_cht_send_t(idx,rbuf)
    implicit none
    integer(kind=4):: idx
    real   (kind=8):: rbuf(*)
!f2py intent(in,out):: idx
!f2py intent(in,out):: rbuf(*)

    call cht_send_t(idx, rbuf)

    return
    end subroutine py_cht_send_t
!-------------------------------------------------------------------------------
!   recv the temperature from the coupler.
!-------------------------------------------------------------------------------
    subroutine py_cht_recv_t(size_rbuf,rbuf)
    implicit none
    integer(kind=4),intent(in):: size_rbuf
    real   (kind=8),intent(in):: rbuf(*)
!f2py intent(in):: size_rbuf
!f2py intent(in):: rbuf(*)

    call cht_recv_t(size_rbuf, rbuf)

    return
    end subroutine py_cht_recv_t
!-------------------------------------------------------------------------------
!   output the cht.
!-------------------------------------------------------------------------------
    subroutine py_cht_wr
    implicit none

    call cht_wr

    return
    end subroutine py_cht_wr
!-------------------------------------------------------------------------------
!   py call to FV unsteady solution.
!-------------------------------------------------------------------------------
    subroutine py_fv_unsteady(uns_iter,iter,residual)
    implicit none
    integer(kind=4),intent(in ):: uns_iter,iter
    real   (kind=8),intent(out):: residual
!f2py integer(kind=4),intent(in ):: uns_iter,iter
!f2py real   (kind=8),intent(out):: residual

    call fv_unsteady_interface(uns_iter,iter,residual)

    return
    end subroutine py_fv_unsteady
!-------------------------------------------------------------------------------
!   py call to HEAT unsteady solution.
!-------------------------------------------------------------------------------
    subroutine py_heat_unsteady(uns_iter,iter,residual)
    implicit none
    integer(kind=4),intent(in ):: uns_iter,iter
    real   (kind=8),intent(out):: residual
!f2py integer(kind=4),intent(in ):: uns_iter,iter
!f2py real   (kind=8),intent(out):: residual

    call heat_unsteady_interface(uns_iter,iter,residual)

    return
    end subroutine py_heat_unsteady
!-------------------------------------------------------------------------------
!   py call to FV finish.
!-------------------------------------------------------------------------------
    subroutine py_fv_finish
    implicit none

    call fv_slv_monitor(.true., 0)

    return
    end subroutine py_fv_finish
!-------------------------------------------------------------------------------
!   py call to HEAT finish.
!-------------------------------------------------------------------------------
    subroutine py_heat_finish
    implicit none

    call heat_slv_monitor(.true., 0)

    return
    end subroutine py_heat_finish
