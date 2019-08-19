!-------------------------------------------------------------------------------
!   GMRES solver, decoupled.
!-------------------------------------------------------------------------------
    subroutine fv_slv_gmres_decoupled(lev,s)
    use var_kind_def
    use var_air, only: gk1
    use var_global, only: rref,uref,pref
    use var_gmres
    use var_slv, only: CFL_lev
    implicit none
    integer(dpI),intent(in):: lev
    real   (dpR):: uc_ref(5),norm_u
    type(type_gmres):: s

    uc_ref(1:5) = (/rref, rref*uref, rref*uref, rref*uref, pref/gk1+0.5d0*rref*uref*uref/)
    if(.not. s%is_allocated)    call kry_setup(lev, s)
    call CFL_evolution(lev)

    call gmres_solution(s, fv_slv_gmres_set_b, fv_slv_gmres_matvec, fv_slv_gmres_precond)
    call fv_slv_gmres_get_x(s%is_scaled, s%ref, s%v_k, 1.0d0)

    return
    contains
!   ----------------------------------------------------------------------------
!   memory allocation for KRY.
!   ----------------------------------------------------------------------------
    subroutine kry_setup(lev,s)
    use var_kind_def
    use var_gmres
    use var_mesh, only: mesh,sec
    implicit none
    integer(dpI),intent(in):: lev
    integer(dpI):: isec,n_vtx
    type(type_gmres):: s

    n_vtx   =  0
    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(sec(isec)%is_int)    n_vtx   =  n_vtx+sec(isec)%n_ele
    end do
    call gmres_set(s, n_vtx, 5, 15, 1, uc_ref, 1.0d-1)
    s%TOL_gmres     =  1.0d-1
    s%is_output_res =  .false.

    return
    end subroutine kry_setup
!   ----------------------------------------------------------------------------
!   CFL evolution strategy.
!   ----------------------------------------------------------------------------
    subroutine CFL_evolution(lev)
    use var_kind_def
    use var_fv
    use var_slv
    implicit none
    integer(dpI),intent(in):: lev
    real   (dpR):: eps,beta,CFL,CFL_min,CFL_max,r

    eps     =  1.0d-2
    beta    =  1.5d0
    CFL_min =  1.0d1
    CFL_max =  1.0d4

    if(ite_wrk .le. 2)  return
    r   =  1.0d0-res_NS_now/res_NS_old-eps
    CFL =  CFL_lev(lev)*beta**r
    CFL_lev(lev)=  max(CFL_min, min(CFL_max, CFL))

    return
    end subroutine CFL_evolution
!   ----------------------------------------------------------------------------
!   set b.
!   ----------------------------------------------------------------------------
    subroutine fv_slv_gmres_set_b(LDA,b)
    use var_kind_def
    use var_fv
    use var_mesh
    use var_parallel
    implicit none
    integer(dpI),intent(in):: LDA
    integer(dpI):: isec,iele,idx
    real   (dpR):: b(*),local(2),global(2),v(5)

    call fv_get_rhs     (lev, .false.)
    call fv_get_residual(lev, res_NS)
    call fv_get_imp_LHS (lev, fv_prec(lev))
    idx     =  1
    local(1)=  0.0d0
    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(.not. sec(isec)%is_int)  cycle
        do iele=1,sec(isec)%n_ele
            b(idx:idx+4)= -fv(isec)%rhs(1:5,iele)
            idx         =  idx+5
        end do

        do iele=1,sec(isec)%n_ele
            fv(isec)%uc0(1:5,iele)  =  fv(isec)%uc(1:5,iele)
            v(1:5)  =  fv(isec)%uc(1:5,iele)/uc_ref(1:5)
            local(1)=  local(1)+v(1)**2+v(2)**2+v(3)**2+v(4)**2+v(5)**2
        end do
    end do
    if(idx-1 .ne. LDA)  stop 'Error: GMRES set_b fails.'

    local(2)=  real(idx-1, dpR)
    if(nprc .le. 1) then
        global  =  local
    else
        call mpi_allreduce(local, global, 2, mpi_dpR, mpi_sum, mpi_comm_world, mpi_err)
    end if
    norm_u  =  sqrt(global(1)/global(2))

    return
    end subroutine fv_slv_gmres_set_b
!   ----------------------------------------------------------------------------
!   get x.
!   ----------------------------------------------------------------------------
    subroutine fv_slv_gmres_get_x(is_scale,ref,x,c)
    use var_kind_def
    use var_fv
    use var_mesh
    implicit none
    logical(dpL),intent(in):: is_scale
    real   (dpR),intent(in):: ref(*),x(*),c
    integer(dpI):: isec,iele,idx

    idx =  1
    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(.not. sec(isec)%is_int)  cycle
        if(is_scale) then
            do iele=1,sec(isec)%n_ele
                fv(isec)%uc(1:5,iele)   =  fv(isec)%uc0(1:5,iele)+c*ref(1:5)*x(idx:idx+4)
                idx                     =  idx+5
            end do
        else
            call DAXPBYZ(5*sec(isec)%n_ele, c, x(idx), 1.0d0, fv(isec)%uc0, fv(isec)%uc)
            idx =  idx+5*sec(isec)%n_ele
        end if
    end do
    call fv_bnd_parallel(lev, .true.)

    return
    end subroutine fv_slv_gmres_get_x
!   ----------------------------------------------------------------------------
!   precondition.
!   ----------------------------------------------------------------------------
    subroutine fv_slv_gmres_precond(LDA,b)
    use var_kind_def
    use var_fv
    use var_prec, only: prec_solve
    implicit none
    integer(dpI),intent(in):: LDA
    real   (dpR):: b(*)

    call DCOPY(LDA, b, 1, fv_prec(lev)%RHS, 1)
    fv_prec(lev)%x  =  0.0d0
    call prec_solve(fv_prec(lev), .true., 1)
    call DCOPY(LDA, fv_prec(lev)%x, 1, b, 1)

    return
    end subroutine fv_slv_gmres_precond
!   ----------------------------------------------------------------------------
!   A*x.
!   ----------------------------------------------------------------------------
    subroutine fv_slv_gmres_matvec(is_scale,n_vtx,bsize,ref,x,Ax)
    use var_kind_def
    use var_fv
    use var_gmres
    use var_mesh
    implicit none
    logical(dpL),intent(in):: is_scale
    integer(dpI),intent(in):: n_vtx,bsize
    real   (dpR),intent(in):: ref(*),x(*)
    integer(dpI):: isec,iele,idx
    real   (dpR):: Ax(*),v(5),eps,norm_x
    real   (dpR),external:: DNRM_mpi

    norm_x  =  DNRM_mpi(n_vtx*bsize, x)
    if(norm_x .le. 0.0d0) then
        Ax(1:n_vtx*bsize)   =  0.0d0
        return
    end if
    eps     =  sqrt(epsilon(1.0d0))*norm_u/norm_x
    call fv_slv_gmres_get_x(is_scale, ref, x, eps)
    call fv_get_rhs(lev, .false.)
    idx     =  1
    v(1:5)  =  ref(1:5)/CFL_lev(lev)
    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(.not. sec(isec)%is_int)  cycle
        do iele=1,sec(isec)%n_ele
            Ax(idx:idx+4)   = (fv(isec)%rhs(1:5,iele)+s%b(idx:idx+4))/eps
!                           & +fv(isec)%LHS_s(iele)*v(1:5)*v_k(idx:idx+4)
            idx             =  idx+5
        end do
    end do

    return
    end subroutine fv_slv_gmres_matvec
    end subroutine fv_slv_gmres_decoupled
!-------------------------------------------------------------------------------
!   GMRES solver, coupled.
!-------------------------------------------------------------------------------
    subroutine fv_slv_gmres(lev,s)
    use var_kind_def
    use var_air, only: gk1
    use var_fv
    use var_global, only: rref,uref,pref
    use var_gmres
    use var_slv
    use var_turb, only: RANS_model,RANS_SA,nut_ref
    implicit none
    integer(dpI),intent(in):: lev
    real   (dpR):: uc_ref(6),norm_u
    type(type_gmres):: s

    if((solver_lev(lev) .ne. solver_gmres) .or. (RANS_model .ne. RANS_SA))  return

    uc_ref(1:6) = (/rref, rref*uref, rref*uref, rref*uref, &
                &  pref/gk1+0.5d0*rref*uref*uref, 1.0d3*nut_ref/)
    if(.not. s%is_allocated)    call kry_setup(lev, s)
    call fv_gmres_set_prec(lev, fv_prec(lev))

!   call fv_gmres_get_prec(lev, fv_prec(lev))

    call CFL_evolution(lev)

    call gmres_solution(s, fv_slv_gmres_set_b, fv_slv_gmres_matvec, fv_slv_gmres_precond)
    call fv_slv_gmres_get_x(s%is_scaled, s%ref, s%v_k, 1.0d0)

    return
    contains
!   ----------------------------------------------------------------------------
!   memory allocation for KRY.
!   ----------------------------------------------------------------------------
    subroutine kry_setup(lev,s)
    use var_mesh, only: mesh,sec
    implicit none
    integer(dpI),intent(in):: lev
    integer(dpI):: isec,n_vtx
    type(type_gmres):: s

    n_vtx   =  0
    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(sec(isec)%is_int)    n_vtx   =  n_vtx+sec(isec)%n_ele
    end do
    call gmres_set(s, n_vtx, 6, 15, 1, uc_ref, 1.0d-1)
    s%TOL_gmres     =  1.0d-1
    s%is_output_res =  .false.

    return
    end subroutine kry_setup
!   ----------------------------------------------------------------------------
!   set prec.
!   ----------------------------------------------------------------------------
    subroutine fv_gmres_set_prec(lev,p)
    use var_global, only: err_mem
    use var_mesh
    use var_mg, only: is_amg,mlev,mg
    use var_prec
    use var_slv
    implicit none
    integer(dpI),intent(in):: lev
    integer(dpI):: isec
    type(type_prec):: p

    if(p%is_spy)    return

    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(.not. sec(isec)%is_int)  cycle
        allocate(fv(isec)%uc0(6,sec(isec)%n_ele), stat=err_mem)
    end do

    if(prec_type .eq. 1) then
        call fv_set_prec(lev, p, 6, 1, .true.)
    elseif(prec_type .eq. 2) then
        call fv_set_prec(lev, p, 6, 3, .true.)
    end if

!   ----------------------------------------------------------------------------
!   agglomeration multigrid.
    if((lev .eq. 0) .and. is_amg) then
        do isec=0,mlev-2
            call prec_set_coarse(p, isec, mg(isec)%ID)
        end do
    end if
!   agglomeration multigrid.
!   ----------------------------------------------------------------------------

    call prec_set_ordering(p)

    return
    end subroutine fv_gmres_set_prec
!   ----------------------------------------------------------------------------
!   CFL evolution strategy.
!   ----------------------------------------------------------------------------
    subroutine CFL_evolution(lev)
    use var_slv
    implicit none
    integer(dpI),intent(in):: lev
    real   (dpR):: eps,beta,CFL,CFL_min,CFL_max,r

    eps     =  1.0d-2
    beta    =  1.5d0
    CFL_min =  1.0d1
    CFL_max =  1.0d4

    if(ite_wrk .le. 2)  return
    r   =  1.0d0-res_NS_now/res_NS_old-eps
    CFL =  CFL_lev(lev)*beta**r
    CFL_lev(lev)=  max(CFL_min, min(CFL_max, CFL))

    return
    end subroutine CFL_evolution
!   ----------------------------------------------------------------------------
!   set b.
!   ----------------------------------------------------------------------------
    subroutine fv_slv_gmres_set_b(LDA,b)
    use var_mesh
    use var_parallel
    implicit none
    integer(dpI),intent(in):: LDA
    integer(dpI):: isec,iele,idx
    real   (dpR):: b(*),local(2),global(2),v(5)

!   ----------------------------------------------------------------------------
!   residual for the Navier-Stokes part.
    call fv_get_rhs     (lev, .false.)
    call fv_get_residual(lev, res_NS)
    idx     =  1
    local(1)=  0.0d0
    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(.not. sec(isec)%is_int)  cycle
        do iele=1,sec(isec)%n_ele
            b(idx:idx+4)= -fv(isec)%rhs(1:5,iele)
            idx         =  idx+6
        end do

        do iele=1,sec(isec)%n_ele
            fv(isec)%uc0(1:5,iele)  =  fv(isec)%uc(1:5,iele)
            v(1:5)  =  fv(isec)%uc(1:5,iele)/uc_ref(1:5)
            local(1)=  local(1)+v(1)**2+v(2)**2+v(3)**2+v(4)**2+v(5)**2
        end do
    end do
!   residual for the Navier-Stokes part.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   residual for the SA part.
    call fv_SA_get_rhs(lev, .false.)
    idx =  6
    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(.not. sec(isec)%is_int)  cycle
        do iele=1,sec(isec)%n_ele
            b(idx)  = -fv(isec)%rhs(1,iele)
            idx     =  idx+6
        end do

        do iele=1,sec(isec)%n_ele
            fv(isec)%uc0(6,iele)=  fv(isec)%turb(1,iele)
            local(1)=  local(1)+(fv(isec)%rhs(1,iele)/uc_ref(6))**2
        end do
    end do
!   residual for the SA part.
!   ----------------------------------------------------------------------------

    if(idx-1 .ne. LDA)  stop 'Error: GMRES set_b fails.'
    local(2)=  real(idx-1, dpR)
    if(nprc .le. 1) then
        global  =  local
    else
        call mpi_allreduce(local, global, 2, mpi_dpR, mpi_sum, mpi_comm_world, mpi_err)
    end if
    norm_u  =  sqrt(global(1)/global(2))

    return
    end subroutine fv_slv_gmres_set_b
!   ----------------------------------------------------------------------------
!   get x.
!   ----------------------------------------------------------------------------
    subroutine fv_slv_gmres_get_x(is_scale,ref,x,c)
    use var_mesh
    implicit none
    logical(dpL),intent(in):: is_scale
    real   (dpR),intent(in):: ref(*),x(*),c
    integer(dpI):: isec,iele,idx

    idx =  1
    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(.not. sec(isec)%is_int)  cycle
        if(is_scale) then
            do iele=1,sec(isec)%n_ele
                fv(isec)%uc(1:5,iele)   =  fv(isec)%uc0(1:5,iele)+c*ref(1:5)*x(idx:idx+4)
                fv(isec)%turb(1,iele)   =  fv(isec)%uc0(6  ,iele)+c*ref(6  )*x(    idx+5)
                idx                     =  idx+6
            end do
        else
            do iele=1,sec(isec)%n_ele
                fv(isec)%uc(1:5,iele)   =  fv(isec)%uc0(1:5,iele)+c*x(idx:idx+4)
                fv(isec)%turb(1,iele)   =  fv(isec)%uc0(6  ,iele)+c*x(    idx+5)
                idx                     =  idx+6
            end do
        end if
    end do
    call fv_bnd_parallel(lev, .true.)
    call fv_rans_bnd_parallel(lev, .true.)

    return
    end subroutine fv_slv_gmres_get_x
!   ----------------------------------------------------------------------------
!   precondition.
!   ----------------------------------------------------------------------------
    subroutine fv_slv_gmres_precond(LDA,b)
    use var_prec, only: prec_solve
    implicit none
    integer(dpI),intent(in):: LDA
    real   (dpR):: b(*)

    call DCOPY(LDA, b, 1, fv_prec(lev)%RHS, 1)
    fv_prec(lev)%x  =  0.0d0
    call prec_solve(fv_prec(lev), .true., 1)
    call DCOPY(LDA, fv_prec(lev)%x, 1, b, 1)

    return
    end subroutine fv_slv_gmres_precond
!   ----------------------------------------------------------------------------
!   A*x.
!   ----------------------------------------------------------------------------
    subroutine fv_slv_gmres_matvec(is_scale,n_vtx,bsize,ref,x,Ax)
    use var_mesh
    implicit none
    logical(dpL),intent(in):: is_scale
    integer(dpI),intent(in):: n_vtx,bsize
    real   (dpR),intent(in):: ref(*),x(*)
    integer(dpI):: isec,iele,idx
    real   (dpR):: Ax(*),v(6),eps,norm_x
    real   (dpR),external:: DNRM_mpi

    norm_x  =  DNRM_mpi(n_vtx*bsize, x)
    if(norm_x .le. 0.0d0) then
        Ax(1:n_vtx*bsize)   =  0.0d0
        return
    end if
    eps     =  sqrt(epsilon(1.0d0))*norm_u/norm_x
    call fv_slv_gmres_get_x(is_scale, ref, x, eps)
    v(1:6)  =  ref(1:6)/CFL_lev(lev)

    call fv_get_rhs(lev, .false.)
    idx     =  1
    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(.not. sec(isec)%is_int)  cycle
        do iele=1,sec(isec)%n_ele
            Ax(idx:idx+4)   = (fv(isec)%rhs(1:5,iele)+s%b(idx:idx+4))/eps
!                           & +fv(isec)%LHS_s(iele)*v(1:5)*v_k(idx:idx+4)
            idx             =  idx+5
        end do
    end do

    call fv_SA_get_rhs(lev, .false.)
    idx     =  6
    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(.not. sec(isec)%is_int)  cycle
        do iele=1,sec(isec)%n_ele
            Ax(idx) = (fv(isec)%rhs(1,iele)+s%b(idx))/eps
!                   & +fv(isec)%LHS_s(iele)*v(1:5)*v_k(idx:idx+4)
            idx     =  idx+6
        end do
    end do

    return
    end subroutine fv_slv_gmres_matvec
    end subroutine fv_slv_gmres
