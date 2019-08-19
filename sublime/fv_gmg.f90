!-------------------------------------------------------------------------------
!   set coarse mesh.
!-------------------------------------------------------------------------------
    subroutine fv_set_coarse_mesh
    use var_kind_def
    use var_fv, only: mlev_fv
    use var_mg, only: is_amg,mlev
    use var_slv
    implicit none
    integer(dpI):: lev

    if(is_amg) then
        if(is_matrix_free .or. (solver_lev(0) .ne. solver_rk_sgs))  &
            &  stop 'Error: current solver does not support AMG.'
        mlev_fv =  1
    else
        mlev_fv =  max(1, min(mlev, 5))
    end if
    if(mlev_fv .le. 1)  return
    do lev=0,mlev_fv-2
        call mesh_get_coarse_level(lev)
    end do

    return
    end subroutine fv_set_coarse_mesh
!-------------------------------------------------------------------------------
!   injection the solution to the coarse level.
!-------------------------------------------------------------------------------
    subroutine fv_gmg_injection(L0)
    use var_kind_def
    use var_fv
    use var_mesh
    use var_slv
    implicit none
    integer(dpI),intent(in):: L0
    integer(dpI):: L1,isec,iele,s,e

    if(L0 .ge. mlev_fv-1)   return
    L1  =  L0+1

    do isec=mesh(L1)%sec_1,mesh(L1)%sec_0
        if(sec(isec)%is_int)    fv(isec)%uc =  0.0d0
    end do

    do isec=mesh(L0)%sec_1,mesh(L0)%sec_0
        if(.not. sec(isec)%is_int)  cycle
        do iele=1,sec(isec)%n_ele
            s   =  sec(isec)%ID_coarse_mesh(1,iele)
            e   =  sec(isec)%ID_coarse_mesh(2,iele)
            fv(s)%uc(1:5,e) =  fv(s)%uc(1:5,e)+fv(isec)%uc(1:5,iele)*sec(isec)%vol(iele)
        end do
    end do

    do isec=mesh(L1)%sec_1,mesh(L1)%sec_0
        if(.not. sec(isec)%is_int)  cycle
        do iele=1,sec(isec)%n_ele
            fv(isec)%uc(1:5,iele)   =  fv(isec)%uc(1:5,iele)/sec(isec)%vol(iele)
        end do
        call DCOPY(5*sec(isec)%n_ele, fv(isec)%uc, 1, fv(isec)%uc_f2c, 1)
    end do

    call fv_bnd_parallel(L1, .true.)

    return
    end subroutine fv_gmg_injection
!-------------------------------------------------------------------------------
!   get forcing function of the coarse level.
!-------------------------------------------------------------------------------
    subroutine fv_gmg_get_ff(L0)
    use var_kind_def
    use var_fv
    use var_mesh
    use var_slv
    implicit none
    integer(dpI),intent(in):: L0
    integer(dpI):: L1,isec,iele,s,e

    if(L0 .ge. mlev_fv-1)   return
    L1  =  L0+1

    call fv_get_rhs(L0, .true.)
    do isec=mesh(L1)%sec_1,mesh(L1)%sec_0
        if(sec(isec)%is_int)    fv(isec)%ff =  0.0d0
    end do

    do isec=mesh(L0)%sec_1,mesh(L0)%sec_0
        if(.not. sec(isec)%is_int)  cycle
        do iele=1,sec(isec)%n_ele
            s   =  sec(isec)%ID_coarse_mesh(1,iele)
            e   =  sec(isec)%ID_coarse_mesh(2,iele)
            fv(s)%ff(1:5,e) =  fv(s)%ff(1:5,e)+fv(isec)%rhs(1:5,iele)
        end do
    end do

    call fv_get_rhs(L1, .false.)
    do isec=mesh(L1)%sec_1,mesh(L1)%sec_0
        if(.not. sec(isec)%is_int)  cycle
        do iele=1,sec(isec)%n_ele
            fv(isec)%ff(1:5,iele)   =  fv(isec)%ff(1:5,iele)-fv(isec)%rhs(1:5,iele)
        end do
    end do

    return
    end subroutine fv_gmg_get_ff
!-------------------------------------------------------------------------------
!   prolongate coarser level correction to the current level.
!-------------------------------------------------------------------------------
    subroutine fv_gmg_prolongation(L0)
    use var_kind_def
    use var_fv
    use var_mesh
    use var_slv
    implicit none
    integer(dpI),intent(in):: L0
    integer(dpI):: L1,isec,iele,s,e

    if(L0 .ge. mlev_fv-1)   return
    L1  =  L0+1

    do isec=mesh(L0)%sec_1,mesh(L0)%sec_0
        if(.not. sec(isec)%is_int)  cycle
        do iele=1,sec(isec)%n_ele
            s   =  sec(isec)%ID_coarse_mesh(1,iele)
            e   =  sec(isec)%ID_coarse_mesh(2,iele)
            fv(isec)%uc(1:5,iele)   =  fv(isec)%uc(1:5,iele) &
                & +(fv(s)%uc(1:5,e)-fv(s)%uc_f2c(1:5,e))
        end do
    end do

    call fv_bnd_parallel(L0, .true.)

    return
    end subroutine fv_gmg_prolongation
!-------------------------------------------------------------------------------
!   FV multigrid cycle.
!-------------------------------------------------------------------------------
    recursive subroutine fv_gmg_cycle(lev)
    use var_kind_def
    use var_fv, only: mlev_fv,fv_gmres
    use var_mg, only: is_amg
    use var_slv
    use var_global, only: sw_slave
    implicit none
    integer(dpI):: lev,iter

!   pre-relaxation: time marching on the current level.
    if(mg_cycle_pre .or. (lev .eq. mlev_fv-1)) then
        if(solver_lev(lev) .eq. solver_exp) then
if(sw_slave) then
            call fv_slv_exp_rk43_ssp_sw(lev)
else
            call fv_slv_exp_rk43_ssp(lev)
end if
        elseif(solver_lev(lev) .eq. solver_rk_sgs) then
            if(is_matrix_free) then
                call fv_slv_rk3_sgs(lev)
            else
if(sw_slave) then
                call fv_slv_imp_sw(lev)
else
                call fv_slv_imp(lev)
end if
            end if
        elseif(solver_lev(lev) .eq. solver_gmres_d) then
            call fv_slv_gmres_decoupled(lev, fv_gmres)
        else
            stop 'Error: solver not supported, FV.'
        end if
    end if
    if((lev .eq. 0) .and. (solver_lev(lev) .gt. solver_exp) .and. is_amg)   return
    if(lev .eq. mlev_fv-1)  return

!   solution injection to the coarser level.
    call fv_gmg_injection(lev)

!   compute rhs on the current level and add forcing function.
    call fv_gmg_get_ff(lev)

    do iter=1,mg_cycle_vw
        call fv_gmg_cycle(lev+1)
    end do

!   prolongate coarser level correction to the current level.
    call fv_gmg_prolongation(lev)

!   post-relaxation: time marching on the current level.
    if(mg_cycle_post) then
        if(solver_lev(lev) .eq. solver_exp) then
            call fv_slv_exp_rk43_ssp(lev)
        elseif(solver_lev(lev) .eq. solver_rk_sgs) then
            if(is_matrix_free) then
                call fv_slv_rk3_sgs(lev)
            else
                call fv_slv_imp(lev)
            end if
        else
            stop 'Error: solver not supported, FV.'
        end if
    end if

    return
    end subroutine fv_gmg_cycle
