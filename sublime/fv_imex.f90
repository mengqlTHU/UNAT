!-------------------------------------------------------------------------------
!   IMEX for FV, initialization.
!-------------------------------------------------------------------------------
    subroutine fv_imex_ini
    use var_kind_def
    use var_fv
    use var_global, only: is_has_cfg,cfg_file,err_mem
    use var_imex
    use var_mesh
    use var_parallel
    implicit none
    integer(dpI):: isec,i,iele,ibuf(10),io_err
    real   (dpR):: max_vol,min_vol

    namelist /imex/ imex_method
    if(is_IMEX_initialized) return

    if(is_has_cfg) then
        if(myid .eq. 0) then
            open(unit=10,file=trim(adjustl(cfg_file)))
            read(unit=10, nml=imex, iostat=io_err)
            if(io_err .gt. 0)   stop 'Error: fails to read namelist:IMEX.'
            close(10)
        end if
        ibuf(1) =  imex_method
        call mpi_bcast(ibuf, 1, mpi_dpI, 0, mpi_comm_world, mpi_err)
        if(myid .ne. 0) then
            imex_method =  ibuf(1)
        end if
    end if
    imex_method =  min(2, max(imex_method, 1))

    call get_ars_232
    do isec=mesh(0)%sec_1,mesh(0)%sec_0
        if(.not. sec(isec)%is_int)  cycle
        allocate(fv(isec)%IMEX_k  (5,n_stage*sec(isec)%n_ele), stat=err_mem)
        allocate(fv(isec)%IMEX_kt (5,n_stage*sec(isec)%n_ele), stat=err_mem)
        allocate(fv(isec)%IMEX_uc0(5,        sec(isec)%n_ele), stat=err_mem)
        allocate(fv(isec)%IMEX_uct(5,        sec(isec)%n_ele), stat=err_mem)
        fv(isec)%IMEX_k     =  0.0d0
        fv(isec)%IMEX_kt    =  0.0d0
        fv(isec)%IMEX_uc0   =  0.0d0
        fv(isec)%IMEX_uct   =  0.0d0
    end do

    if(imex_method .eq. 1) then
        max_vol = sec(1)%vol(1)
        min_vol = sec(1)%vol(1)
        do isec=mesh(0)%sec_1,mesh(0)%sec_0
            if(.not. sec(isec)%is_int)  cycle
            do iele = 1,sec(isec)%n_ele
                if(sec(isec)%vol(iele) .lt. min_vol) then
                    min_vol =  sec(isec)%vol(iele)
                else if(sec(isec)%vol(iele) .gt. max_vol) then
                    max_vol =  sec(isec)%vol(iele)
                end if
            end do
        end do

        i   =  0
        do isec=mesh(0)%sec_1,mesh(0)%sec_0
            if(.not. sec(isec)%is_int)  cycle

            if(.not. allocated(fv(isec)%IMEX_id)) then
                allocate(fv(isec)%IMEX_id(sec(isec)%n_ele), stat=err_mem)
                fv(isec)%IMEX_id=  0
            end if

            do iele=1,sec(isec)%n_ele
                if(sec(isec)%vol(iele) .gt. 5.0d1*min_vol)  cycle
                i   =  i+1
                fv(isec)%IMEX_id(iele)  =  i
            end do
        end do

        call fv_imex_set_prec(0, fv_imex_prec)
    end if

    is_IMEX_initialized =  .true.

    return
    end subroutine fv_imex_ini
!-------------------------------------------------------------------------------
!   IMEX for FV.
!-------------------------------------------------------------------------------
    subroutine fv_slv_imex
    use var_imex
    use var_uns_cal
    implicit none

    if(.not. (is_uns_cal_now .and. (uns_method .eq. uns_IMEX))) return
    if(.not. is_IMEX_initialized)   call fv_imex_ini

    if(IMEX_method .eq. 1) then
        call fv_slv_imex_region
    else
        call fv_slv_imex_term
    end if

    return
    end subroutine fv_slv_imex
!-------------------------------------------------------------------------------
!   IMEX for FV, splitting by convective/viscous term.
!-------------------------------------------------------------------------------
    subroutine fv_slv_imex_term
    use var_kind_def
    use var_fv
    use var_imex
    use var_mesh
    use var_uns_cal
    implicit none
    integer(dpI):: isec,iele,j

    if(.not. (is_uns_cal_now .and. (uns_method .eq. uns_IMEX))) return
    if(.not. is_IMEX_initialized)   call fv_imex_ini

!   kt  =  0; k =  0; u0=  u
    do isec=mesh(0)%sec_1,mesh(0)%sec_0
        if(.not. sec(isec)%is_int)  cycle
        fv(isec)%IMEX_k     =  0.0d0
        fv(isec)%IMEX_kt    =  0.0d0
        call DCOPY(5*sec(isec)%n_ele, fv(isec)%uc, 1, fv(isec)%IMEX_uc0, 1)
    end do

    do i_stage=1,n_stage
!       u   =  u0 if(i>1)
        if(i_stage .gt. 1) then
            do isec=mesh(0)%sec_1,mesh(0)%sec_0
                if(.not. sec(isec)%is_int)  cycle
                call DCOPY(5*sec(isec)%n_ele, fv(isec)%IMEX_uc0, 1, fv(isec)%uc, 1)
            end do
        end if

!       the contribution from the explicit part.
        do j=1,i_stage-1
!           u   =  u+dt*A_ex(i,j)*kt(j)
            if(abs(A_ex(i_stage,j)) .le. 0.0d0) cycle
            do isec=mesh(0)%sec_1,mesh(0)%sec_0
                if(.not. sec(isec)%is_int)  cycle
                call DAXPY(5*sec(isec)%n_ele, A_ex(i_stage,j)*dt_uns, &
                    &  fv(isec)%IMEX_kt(1,1+(j-1)*sec(isec)%n_ele), 1, &
                    &  fv(isec)%uc, 1)
            end do
        end do

!       the contribution from the implicit part, without the last one.
        do j=1,i_stage-1
!           u   =  u+dt*A_im(i,j)*k(j)
            if(abs(A_im(i_stage,j)) .le. 0.0d0) cycle
            do isec=mesh(0)%sec_1,mesh(0)%sec_0
                if(.not. sec(isec)%is_int)  cycle
                call DAXPY(5*sec(isec)%n_ele, A_im(i_stage,j)*dt_uns, &
                    &  fv(isec)%IMEX_k (1,1+(j-1)*sec(isec)%n_ele), 1, &
                    &  fv(isec)%uc, 1)
            end do
        end do
        call fv_bnd_parallel(0, .true.)

        if(abs(A_im(i_stage,i_stage)) .le. 0.0d0) then
!           $g(u+\delta t*A_im(i,i)*k_i)=k_i$
            call imex_get_implicit_rhs(0)
            do isec=mesh(0)%sec_1,mesh(0)%sec_0
                if(.not. sec(isec)%is_int)  cycle
                j   =  1+(i_stage-1)*sec(isec)%n_ele
                call DCOPY(5*sec(isec)%n_ele, fv(isec)%rhs, 1, fv(isec)%IMEX_k(1,j), 1)
            end do
        else
!           then we have to solve $g(u+\delta t*A_im(i,i)*k_i)=k_i$ for $k_i$.
            call imex_get_implicit_k(0)

!           u   =  u+dt*A_im(i,i)*k(i)
        end if

!       kt(i)   =  f(u)
        call imex_get_explicit_rhs(0)
        do isec=mesh(0)%sec_1,mesh(0)%sec_0
            if(.not. sec(isec)%is_int)  cycle
            do iele=1,sec(isec)%n_ele
                j   =  iele+(i_stage-1)*sec(isec)%n_ele
                fv(isec)%IMEX_kt(1:5,j) = -fv(isec)%rhs(1:5,iele)/sec(isec)%vol(iele)
            end do
        end do
    end do

!   u   =  u0
    do isec=mesh(0)%sec_1,mesh(0)%sec_0
        if(.not. sec(isec)%is_int)  cycle
        call DCOPY(5*sec(isec)%n_ele, fv(isec)%IMEX_uc0, 1, fv(isec)%uc, 1)
    end do

    do i_stage=1,n_stage
!       u   =  u+dt*b_im(i)*k (i)
        if(abs(b_im(i_stage)) .le. 0.0d0)   cycle
        do isec=mesh(0)%sec_1,mesh(0)%sec_0
            if(.not. sec(isec)%is_int)  cycle
            call DAXPY(5*sec(isec)%n_ele, b_im(i_stage)*dt_uns, &
                &  fv(isec)%IMEX_k (1,1+(i_stage-1)*sec(isec)%n_ele), 1, fv(isec)%uc, 1)
        end do
    end do
    do i_stage=1,n_stage
!       u   =  u+dt*b_ex(i)*kt(i)
        if(abs(b_ex(i_stage)) .le. 0.0d0)   cycle
        do isec=mesh(0)%sec_1,mesh(0)%sec_0
            if(.not. sec(isec)%is_int)  cycle
            call DAXPY(5*sec(isec)%n_ele, b_ex(i_stage)*dt_uns, &
                &  fv(isec)%IMEX_kt(1,1+(i_stage-1)*sec(isec)%n_ele), 1, fv(isec)%uc, 1)
        end do
    end do
    call fv_bnd_parallel(0, .true.)

    return
    contains
!   ----------------------------------------------------------------------------
!   get the explicit part of RHS.
!   ----------------------------------------------------------------------------
    subroutine imex_get_explicit_rhs(lev)
    use var_turb, only: is_LES_now,LES_model
    implicit none
    integer(dpI),intent(in):: lev
    integer(dpI):: isec

!   only the convective term.
    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(sec(isec)%is_int)    fv(isec)%rhs=  0.0d0
    end do

    if(is_LES_now .and. (LES_model .ge. 1)) call fv_les(lev)
    if(is_shock_sensor) call fv_get_shock_sensor(lev)

    if(is_kep) then
        call fv_get_rhs_kep(lev)
    else
        if(rhs_lev(lev) .eq. 1) then
            if(is_jst_su2) then
                call fv_get_rhs_jst_su2(lev, .true., .true.)
            else
                call fv_get_rhs_jst    (lev, .true., .true.)
            end if
        else
            call fv_get_rhs_upw(lev, .true., .true.)
        end if
    end if

    return
    end subroutine imex_get_explicit_rhs
!   ----------------------------------------------------------------------------
!   get the implicit part of RHS.
!   ----------------------------------------------------------------------------
    subroutine imex_get_implicit_rhs(lev)
    use var_slv, only: is_vis_cal
    implicit none
    integer(dpI),intent(in):: lev
    integer(dpI):: isec

!   only the viscous term.
    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(sec(isec)%is_int)    fv(isec)%rhs=  0.0d0
    end do
    if(is_vis_cal)  call fv_get_rhs_vis(lev, .true.)

    return
    end subroutine imex_get_implicit_rhs
!   ----------------------------------------------------------------------------
!   slv the nonlinear equation for implicit part.
!   ----------------------------------------------------------------------------
    subroutine imex_get_implicit_k(lev)
    use var_slv
    use var_parallel
    implicit none
    integer(dpI),intent(in):: lev
    integer(dpI):: isec,iele
    real   (dpR):: IMEX_a1,res_IMEX,res_IMEX_1

    do isec=mesh(0)%sec_1,mesh(0)%sec_0
        if(.not. sec(isec)%is_int)  cycle
        call DCOPY(5*sec(isec)%n_ele, fv(isec)%uc, 1, fv(isec)%IMEX_uct, 1)
    end do

    LHS_uns_c1  =  1.0d0/ A_im(i_stage,i_stage)
    IMEX_a1     =  1.0d0/(A_im(i_stage,i_stage)*dt_uns)
    do ite_wrk=1,max_subiter_uns
        call imex_get_implicit_rhs(lev)
        do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
            if(.not. sec(isec)%is_int)  cycle
            do iele=1,sec(isec)%n_ele
                fv(isec)%rhs(1:5,iele)  =  fv(isec)%rhs(1:5,iele) &
                    & +IMEX_a1*sec(isec)%vol(iele) &
                    &*(fv(isec)%uc(1:5,iele)-fv(isec)%IMEX_uct(1:5,iele))
            end do
        end do

        call fv_get_residual(lev, res_NS)
        call mpi_allreduce(res_NS**2,res_IMEX,1,mpi_dpR,mpi_sum,mpi_comm_world,mpi_err)
        res_IMEX=  sqrt(res_IMEX)
        if(ite_wrk .eq. 1)  res_IMEX_1  =  res_IMEX

        call fv_get_imp_LHS(lev, fv_prec(lev))

        call iteration(lev, fv_prec(lev))

        do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
            if(sec(isec)%is_int)    fv(isec)%uc =  fv(isec)%uc+fv(isec)%duc
        end do
        call fv_bnd_parallel(lev, .true.)

        is_converged=  res_IMEX .le. rhs_cvg_uns*res_IMEX_1

        if(is_converged)    exit
    end do
!   to be consistent with FV_SLV_MONITOR.
    res_NS  = (res_IMEX/res_IMEX_1)**2/real(nprc, dpR)

    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(.not. sec(isec)%is_int)  cycle
        do iele=1,sec(isec)%n_ele
            j   =  iele+(i_stage-1)*sec(isec)%n_ele
            fv(isec)%IMEX_k(1:5,j)  =  IMEX_a1* &
                & (fv(isec)%uc(1:5,iele)-fv(isec)%IMEX_uct(1:5,iele))
        end do
    end do

    return
    end subroutine imex_get_implicit_k
!   ----------------------------------------------------------------------------
!   prec iteration.
!   ----------------------------------------------------------------------------
    subroutine iteration(lev,p)
    use var_prec
    implicit none
    integer(dpI),intent(in):: lev
    integer(dpI):: isec,i,iele
    type(type_prec):: p

    i   =  1
    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(.not. sec(isec)%is_int)  cycle
        do iele=1,sec(isec)%n_ele
            p%RHS(i:i+4)= -fv(isec)%rhs(1:5,iele)
            i           =  i+5
        end do
    end do

    p%x =  0.0d0
    call prec_solve(p, .true., 1)

    i   =  1
    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(.not. sec(isec)%is_int)  cycle
        call DCOPY(5*sec(isec)%n_ele, p%x(i), 1, fv(isec)%duc, 1)
        i   =  i+5*sec(isec)%n_ele
    end do

    return
    end subroutine iteration
    end subroutine fv_slv_imex_term
!-------------------------------------------------------------------------------
!   set LHS for the implicit part, IMEX.
!-------------------------------------------------------------------------------
    subroutine fv_imex_set_prec(lev,p)
    use var_kind_def
    use var_fv
    use var_global, only: err_mem
    use var_mesh
    use var_prec
    use var_slv, only: prec_type
    implicit none
    integer(dpI),intent(in):: lev
    integer(dpI):: n_ele,isec,iele,n_jA,n_nb,nb(100),sR,eR,s_donor,e_donor,j,im,per(3)
    type(type_prec):: p

    if(p%is_spy)    return

!   ----------------------------------------------------------------------------
!   get the size of iA.
    n_ele   =  0
    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(.not. sec(isec)%is_int)  cycle
        do iele=1,sec(isec)%n_ele
            if(fv(isec)%IMEX_ID(iele) .gt. 0)   n_ele   =  n_ele+1
        end do
    end do
    call prec_add(p, n_ele, 5, prec_type)
    allocate(p%iA(p%nvtx+1), stat=err_mem)
    if(err_mem .ne. 0)  stop 'Error: fails to allocate memory, fv_imex_set_prec.'
!   get the size of iA.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   get the size of jA.
    p%iA=  0
    n_jA=  0
    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(.not. sec(isec)%is_int)  cycle
        do iele=1,sec(isec)%n_ele
            if(fv(isec)%IMEX_ID(iele) .le. 0)   cycle
            n_jA=  n_jA+1

            do j=sec(isec)%iA_face_neighbour(iele),sec(isec)%iA_face_neighbour(iele+1)-1
                im  =  sec(isec)%jA_face_neighbour(j)
                if(im .gt. 0) then
                    sR  =  mesh(lev)%mortar_LR(3, im)
                    eR  =  mesh(lev)%mortar_LR(4, im)
                else
                    sR  =  mesh(lev)%mortar_LR(1,-im)
                    eR  =  mesh(lev)%mortar_LR(2,-im)
                end if
                if(sec(sR)%is_int) then
                    if(fv(sR)%IMEX_id(eR) .gt. 0)   n_jA=  n_jA+1
                elseif(sec(sR)%is_ghost .and. .true.) then
                    call get_ID_ghost_ele(sR, eR, s_donor, e_donor, per)
                    if(s_donor .gt. 0) then
                        if(fv(s_donor)%IMEX_id(e_donor) .gt. 0) n_jA=  n_jA+1
                    end if
                end if
            end do

            p%iA(fv(isec)%IMEX_id(iele)+1)  =  n_jA+1
        end do
    end do
    allocate(p%jA(n_jA), stat=err_mem)
    if(err_mem .ne. 0)  stop 'Error: fails to allocate memory, fv_imex_set_prec.'
!   get the size of jA.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   sparse pattern of A.
    p%iA    =  0
    p%iA(1) =  1
    n_jA    =  1
    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(.not. sec(isec)%is_int)  cycle
        do iele=1,sec(isec)%n_ele
            if(fv(isec)%IMEX_ID(iele) .le. 0)   cycle
            n_nb    =  1
            nb(1)   =  fv(isec)%IMEX_id(iele)

            do j=sec(isec)%iA_face_neighbour(iele),sec(isec)%iA_face_neighbour(iele+1)-1
                im  =  sec(isec)%jA_face_neighbour(j)
                if(im .gt. 0) then
                    sR  =  mesh(lev)%mortar_LR(3, im)
                    eR  =  mesh(lev)%mortar_LR(4, im)
                else
                    sR  =  mesh(lev)%mortar_LR(1,-im)
                    eR  =  mesh(lev)%mortar_LR(2,-im)
                end if
                if(sec(sR)%is_int) then
                    if(fv(sR)%IMEX_id(eR) .gt. 0) then
                        n_nb    =  n_nb+1
                        nb(n_nb)=  fv(sR)%IMEX_id(eR)
                    end if
                elseif(sec(sR)%is_ghost .and. .true.) then
                    call get_ID_ghost_ele(sR, eR, s_donor, e_donor, per)
                    if(s_donor .gt. 0) then
                        if(fv(s_donor)%IMEX_id(e_donor) .gt. 0) then
                            n_nb    =  n_nb+1
                            nb(n_nb)=  fv(s_donor)%IMEX_id(e_donor)
                        end if
                    end if
                end if
            end do
            call simplify_series(n_nb, 1, 1, nb)
            call ICOPY(n_nb, nb, 1, p%jA(n_jA), 1)
            n_jA                            =  n_jA+n_nb
            p%iA(fv(isec)%IMEX_id(iele)+1)  =  n_jA
        end do
    end do
!   sparse pattern of A.
!   ----------------------------------------------------------------------------

    if(.not. allocated(p%A)) then
        allocate(p%A(p%bsize,p%bsize*(n_jA-1)), stat=err_mem)
        if(err_mem .ne. 0)  stop 'Error: fails to allocate memory, fv_set_prec.'
        p%A =  0.0d0
    end if

    p%is_spy=  .true.

    call prec_set_ordering(p)

    return
    end subroutine fv_imex_set_prec
!-------------------------------------------------------------------------------
!   IMEX for FV, splitting by exp/imp region.
!-------------------------------------------------------------------------------
    subroutine fv_slv_imex_region
    use var_kind_def
    use var_fv
    use var_imex
    use var_mesh
    use var_uns_cal
    implicit none
    integer(dpI):: isec,iele,j

    if(.not. (is_uns_cal_now .and. (uns_method .eq. uns_IMEX))) return
    if(.not. is_IMEX_initialized)   call fv_imex_ini

!   kt  =  0; k =  0; u0=  u
    do isec=mesh(0)%sec_1,mesh(0)%sec_0
        if(.not. sec(isec)%is_int)  cycle
        fv(isec)%IMEX_k     =  0.0d0
        fv(isec)%IMEX_kt    =  0.0d0
        call DCOPY(5*sec(isec)%n_ele, fv(isec)%uc, 1, fv(isec)%IMEX_uc0, 1)
    end do

    do i_stage=1,n_stage
!       u   =  u0 if(i>1)
        if(i_stage .gt. 1) then
            do isec=mesh(0)%sec_1,mesh(0)%sec_0
                if(.not. sec(isec)%is_int)  cycle
                call DCOPY(5*sec(isec)%n_ele, fv(isec)%IMEX_uc0, 1, fv(isec)%uc, 1)
            end do
        end if

!       the contribution from the explicit part.
        do j=1,i_stage-1
!           u   =  u+dt*A_ex(i,j)*kt(j)
            if(abs(A_ex(i_stage,j)) .le. 0.0d0) cycle
            do isec=mesh(0)%sec_1,mesh(0)%sec_0
                if(.not. sec(isec)%is_int)  cycle
                call DAXPY(5*sec(isec)%n_ele, A_ex(i_stage,j)*dt_uns, &
                    &  fv(isec)%IMEX_kt(1,1+(j-1)*sec(isec)%n_ele), 1, &
                    &  fv(isec)%uc, 1)
            end do
        end do

!       the contribution from the implicit part, without the last one.
        do j=1,i_stage-1
!           u   =  u+dt*A_im(i,j)*k(j)
            if(abs(A_im(i_stage,j)) .le. 0.0d0) cycle
            do isec=mesh(0)%sec_1,mesh(0)%sec_0
                if(.not. sec(isec)%is_int)  cycle
                call DAXPY(5*sec(isec)%n_ele, A_im(i_stage,j)*dt_uns, &
                    &  fv(isec)%IMEX_k (1,1+(j-1)*sec(isec)%n_ele), 1, &
                    &  fv(isec)%uc, 1)
            end do
        end do
        call fv_bnd_parallel(0, .true.)

        if(abs(A_im(i_stage,i_stage)) .le. 0.0d0) then
!           $g(u+\delta t*A_im(i,i)*k_i)=k_i$
            call imex_get_implicit_rhs(0)
            do isec=mesh(0)%sec_1,mesh(0)%sec_0
                if(.not. sec(isec)%is_int)  cycle
                j   =  1+(i_stage-1)*sec(isec)%n_ele
                call DCOPY(5*sec(isec)%n_ele, fv(isec)%rhs, 1, fv(isec)%IMEX_k(1,j), 1)
            end do
        else
!           then we have to solve $g(u+\delta t*A_im(i,i)*k_i)=k_i$ for $k_i$.
            call imex_get_implicit_k(0)

!           u   =  u+dt*A_im(i,i)*k(i)
        end if

!       kt(i)   =  f(u)
        call imex_get_explicit_rhs(0)
        do isec=mesh(0)%sec_1,mesh(0)%sec_0
            if(.not. sec(isec)%is_int)  cycle
            do iele=1,sec(isec)%n_ele
                j   =  iele+(i_stage-1)*sec(isec)%n_ele
                fv(isec)%IMEX_kt(1:5,j) = -fv(isec)%rhs(1:5,iele)/sec(isec)%vol(iele)
            end do
        end do
    end do

!   u   =  u0
    do isec=mesh(0)%sec_1,mesh(0)%sec_0
        if(.not. sec(isec)%is_int)  cycle
        call DCOPY(5*sec(isec)%n_ele, fv(isec)%IMEX_uc0, 1, fv(isec)%uc, 1)
    end do

    do i_stage=1,n_stage
!       u   =  u+dt*b_im(i)*k (i)
        if(abs(b_im(i_stage)) .le. 0.0d0)   cycle
        do isec=mesh(0)%sec_1,mesh(0)%sec_0
            if(.not. sec(isec)%is_int)  cycle
            call DAXPY(5*sec(isec)%n_ele, b_im(i_stage)*dt_uns, &
                &  fv(isec)%IMEX_k (1,1+(i_stage-1)*sec(isec)%n_ele), 1, fv(isec)%uc, 1)
        end do
    end do
    do i_stage=1,n_stage
!       u   =  u+dt*b_ex(i)*kt(i)
        if(abs(b_ex(i_stage)) .le. 0.0d0)   cycle
        do isec=mesh(0)%sec_1,mesh(0)%sec_0
            if(.not. sec(isec)%is_int)  cycle
            call DAXPY(5*sec(isec)%n_ele, b_ex(i_stage)*dt_uns, &
                &  fv(isec)%IMEX_kt(1,1+(i_stage-1)*sec(isec)%n_ele), 1, fv(isec)%uc, 1)
        end do
    end do
    call fv_bnd_parallel(0, .true.)

    return
    contains
!   ----------------------------------------------------------------------------
!   get the explicit part of RHS.
!   ----------------------------------------------------------------------------
    subroutine imex_get_explicit_rhs(lev)
    implicit none
    integer(dpI),intent(in):: lev

    call fv_get_rhs_upw_imex(lev, .true., .true., .true.)
    call fv_get_rhs_vis_imp (lev, .true.,         .true.)

    return
    end subroutine imex_get_explicit_rhs
!   ----------------------------------------------------------------------------
!   get the implicit part of RHS.
!   ----------------------------------------------------------------------------
    subroutine imex_get_implicit_rhs(lev)
    implicit none
    integer(dpI),intent(in):: lev

    call fv_get_rhs_upw_imex(lev, .true., .true., .false.)
    call fv_get_rhs_vis_imp (lev, .true.,         .false.)

    return
    end subroutine imex_get_implicit_rhs
!   ----------------------------------------------------------------------------
!   slv the nonlinear equation for implicit part.
!   ----------------------------------------------------------------------------
    subroutine imex_get_implicit_k(lev)
    use var_slv
    use var_parallel
    implicit none
    integer(dpI),intent(in):: lev
    integer(dpI):: isec,iele
    real   (dpR):: IMEX_a1,res_IMEX,res_IMEX_1

    do isec=mesh(0)%sec_1,mesh(0)%sec_0
        if(.not. sec(isec)%is_int)  cycle
        call DCOPY(5*sec(isec)%n_ele, fv(isec)%uc, 1, fv(isec)%IMEX_uct, 1)
    end do

    LHS_uns_c1  =  1.0d0/ A_im(i_stage,i_stage)
    IMEX_a1     =  1.0d0/(A_im(i_stage,i_stage)*dt_uns)
    do ite_wrk=1,max_subiter_uns
        call imex_get_implicit_rhs(lev)
        do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
            if(.not. sec(isec)%is_int)  cycle
            do iele=1,sec(isec)%n_ele
                fv(isec)%rhs(1:5,iele)  =  fv(isec)%rhs(1:5,iele) &
                    & +IMEX_a1*sec(isec)%vol(iele) &
                    &*(fv(isec)%uc(1:5,iele)-fv(isec)%IMEX_uct(1:5,iele))
            end do
        end do

        call fv_get_residual(lev, res_NS)
        call mpi_allreduce(res_NS**2,res_IMEX,1,mpi_dpR,mpi_sum,mpi_comm_world,mpi_err)
        res_IMEX=  sqrt(res_IMEX)
        if(ite_wrk .eq. 1)  res_IMEX_1  =  res_IMEX

!       call fv_get_imp_LHS(lev, fv_imex_prec)
        call fv_imex_get_lhs(lev, fv_imex_prec)

        call iteration(lev, fv_imex_prec)

        do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
            if(sec(isec)%is_int)    fv(isec)%uc =  fv(isec)%uc+fv(isec)%duc
        end do
        call fv_bnd_parallel(lev, .true.)

        is_converged=  res_IMEX .le. rhs_cvg_uns*res_IMEX_1

        if(is_converged)    exit
    end do
!   to be consistent with FV_SLV_MONITOR.
    res_NS  = (res_IMEX/res_IMEX_1)**2/real(nprc, dpR)

    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(.not. sec(isec)%is_int)  cycle
        do iele=1,sec(isec)%n_ele
            j   =  iele+(i_stage-1)*sec(isec)%n_ele
            fv(isec)%IMEX_k(1:5,j)  =  IMEX_a1* &
                & (fv(isec)%uc(1:5,iele)-fv(isec)%IMEX_uct(1:5,iele))
        end do
    end do

    return
    end subroutine imex_get_implicit_k
!   ----------------------------------------------------------------------------
!   prec iteration.
!   ----------------------------------------------------------------------------
    subroutine iteration(lev,p)
    use var_prec
    implicit none
    integer(dpI),intent(in):: lev
    integer(dpI):: isec,i,iele
    type(type_prec):: p

    i   =  1
    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(.not. sec(isec)%is_int)  cycle
        do iele=1,sec(isec)%n_ele
            if(fv(isec)%IMEX_id(iele) .le. 0)   cycle
            p%RHS(i:i+4)= -fv(isec)%rhs(1:5,iele)
            i           =  i+5
        end do
    end do

    p%x =  0.0d0
    call prec_solve(p, .true., 1)

    i   =  1
    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(.not. sec(isec)%is_int)  cycle
        do iele=1,sec(isec)%n_ele
            if(fv(isec)%IMEX_id(iele) .le. 0)   cycle
            fv(isec)%duc(1:5,iele)  =  p%x(i:i+4)
            i                       =  i+5
        end do
    end do

    return
    end subroutine iteration
!   ----------------------------------------------------------------------------
!   cal Rhs, convective part, fv, IMEX.
!   ----------------------------------------------------------------------------
    subroutine fv_get_rhs_upw_imex(lev,calD,addD,is_exp)
    use var_kind_def
    use var_fv
    use var_global, only: n_dim,R16,R23
    use var_mesh
    use var_parallel
    use var_temp, only: uLf,uRf,vgf,fac_1d,rhsl,rhsD,u=>u_structured
    implicit none
    logical(dpL),intent(in):: calD,addD,is_exp
    integer(dpI),intent(in):: lev
    logical(dpL):: savD,is_cal
    integer(dpI):: im,s_LL,e_LL,sL,eL,sR,eR,s_RR,e_RR,j
    real   (dpR):: dL(3),dR(3),limiter_L,limiter_R,duL(5),duR(5),du(5)

    dL          =  0.0d0
    dR          =  0.0d0
    limiter_L   =  1.0d0
    limiter_R   =  1.0d0
    savD=  calD .and. (.not. addD)

    call fv_get_rhs_boundary_imex(lev, calD, addD, is_exp)

    do im=1+mesh(lev)%n_mortar_b,mesh(lev)%n_mortar_ss
        s_LL=  mesh(lev)%mortar_structured_stencil(1,im)
        e_LL=  mesh(lev)%mortar_structured_stencil(2,im)
        sL  =  mesh(lev)%mortar_structured_stencil(3,im)
        eL  =  mesh(lev)%mortar_structured_stencil(4,im)
        sR  =  mesh(lev)%mortar_structured_stencil(5,im)
        eR  =  mesh(lev)%mortar_structured_stencil(6,im)
        s_RR=  mesh(lev)%mortar_structured_stencil(7,im)
        e_RR=  mesh(lev)%mortar_structured_stencil(8,im)
        fac_1d(1:4,1)   =  mesh(lev)%mortar_n_vg(1:4,im)
        vgf   (1:3,1)   =  mesh(lev)%mortar_n_vg(6:8,im)

        if(sec(sR)%is_int) then
            is_cal  =((fv(sL)%IMEX_ID(eL) .le. 0) .eqv. is_exp) .or. &
                    &((fv(sR)%IMEX_ID(eR) .le. 0) .eqv. is_exp)
        else
            is_cal  =((fv(sL)%IMEX_ID(eL) .le. 0) .eqv. is_exp)
        end if
        if(.not. is_cal)    cycle

        u(1:5,0)=  fv(s_LL)%u(1:5,e_LL)
        u(1:5,1)=  fv(sL  )%u(1:5,eL  )
        u(1:5,2)=  fv(sR  )%u(1:5,eR  )
        u(1:5,3)=  fv(s_RR)%u(1:5,e_RR)
        if(is_limiter_on) then
            call recons_3rd_c(fac_1d, vgf, uLf, uRf)
        else
            if(fv(sL)%order .le. 1) then
                uLf(1:5,1)  =  u(1:5,1)
                uRf(1:5,1)  =  u(1:5,2)
            else
                uLf(1:5,1)  =(     -u(1:5,0)+5.0d0*u(1:5,1)+2.0d0*u(1:5,2))*R16
                uRf(1:5,1)  =(2.0d0*u(1:5,1)+5.0d0*u(1:5,2)      -u(1:5,3))*R16
            end if
        end if

        if(rhs_lev(lev) .eq. 2) then
            call rhs_conv_roe (calD, addD, 1)
        elseif(rhs_lev(lev) .eq. 3) then
            call rhs_conv_hllc(calD, addD, 1)
        else
            call rhs_conv_ausm(calD, addD, 1)
        end if

        if((fv(sL)%IMEX_ID(eL) .le. 0) .eqv. is_exp) then
            fv(sL)%rhs(1:5,eL)  =  fv(sL)%rhs(1:5,eL)+rhsl(1:5,1)
            if(savD)    fv(sL)%duc(1:5,eL)  =  fv(sL)%duc(1:5,eL)+rhsD(1:5,1)
        end if
        if(((fv(sR)%IMEX_ID(eR) .le. 0) .eqv. is_exp) .and. sec(sR)%is_int) then
            fv(sR)%rhs(1:5,eR)  =  fv(sR)%rhs(1:5,eR)-rhsl(1:5,1)
            if(sec(sR)%is_int)  fv(sR)%duc(1:5,eR)  =  fv(sR)%duc(1:5,eR)-rhsD(1:5,1)
        end if
    end do

    do im=1+mesh(lev)%n_mortar_ss,mesh(lev)%n_mortar
        sL  =  mesh(lev)%mortar_LR(1,im)
        eL  =  mesh(lev)%mortar_LR(2,im)
        sR  =  mesh(lev)%mortar_LR(3,im)
        eR  =  mesh(lev)%mortar_LR(4,im)
        fac_1d(1:4,1)   =  mesh(lev)%mortar_n_vg(1:4,im)
        vgf   (1:3,1)   =  mesh(lev)%mortar_n_vg(6:8,im)

        if(sec(sR)%is_int) then
            is_cal  =((fv(sL)%IMEX_ID(eL) .le. 0) .eqv. is_exp) .or. &
                    &((fv(sR)%IMEX_ID(eR) .le. 0) .eqv. is_exp)
        else
            is_cal  =((fv(sL)%IMEX_ID(eL) .le. 0) .eqv. is_exp)
        end if
        if(.not. is_cal)    cycle

        uLf(1:5,1)  =  fv(sL)%u (1:5,eL)
        uRf(1:5,1)  =  fv(sR)%u (1:5,eR)
        if(fv(sL)%order .gt. 1) then
            dL(1:n_dim) =  mesh(lev)%mortar_cen(1:n_dim,im)-sec(sL)%cen(1:n_dim,eL)
            dR(1:n_dim) =  mesh(lev)%mortar_cen(1:n_dim,im)-sec(sR)%cen(1:n_dim,eR)
            if(is_limiter_on)   limiter_L   =  fv(sL)%gra(19,eL)
            if(is_limiter_on)   limiter_R   =  fv(sR)%gra(19,eR)

            do j=1,5
                duL(j)  =  fv(sL)%gra(3*j-2,eL)*dL(1)+fv(sL)%gra(3*j-1,eL)*dL(2) &
                        & +fv(sL)%gra(3*j  ,eL)*dL(3)
                duR(j)  =  fv(sR)%gra(3*j-2,eR)*dR(1)+fv(sR)%gra(3*j-1,eR)*dR(2) &
                        & +fv(sR)%gra(3*j  ,eR)*dR(3)
            end do
            if(fv(sL)%order .eq. 2) then
                uLf(1:5,1)  =  uLf(1:5,1)+limiter_L*duL
                uRf(1:5,1)  =  uRf(1:5,1)+limiter_R*duR
            else
                du (1:5  )  =  uRf(1:5,1)-uLf(1:5,1)
                uLf(1:5,1)  =  uLf(1:5,1)+R23*duL(1:5)+R16*du(1:5)
                uRf(1:5,1)  =  uRf(1:5,1)+R23*duR(1:5)-R16*du(1:5)
            end if

            if((uLf(1,1) .le. 0.0d0) .or. (uLf(5,1) .le. 0.0d0) .or. &
             & (uRf(1,1) .le. 0.0d0) .or. (uRf(5,1) .le. 0.0d0)) then
                uLf(1:5,1)  =  fv(sL)%u(1:5,eL)
                uRf(1:5,1)  =  fv(sR)%u(1:5,eR)
            end if
        end if

        if(rhs_lev(lev) .eq. 2) then
            call rhs_conv_roe (calD, addD, 1)
        elseif(rhs_lev(lev) .eq. 3) then
            call rhs_conv_hllc(calD, addD, 1)
        elseif(rhs_lev(lev) .eq. 4) then
            call rhs_conv_ausm(calD, addD, 1)
        else
            call rhs_conv_lf  (calD, addD, 1)
        end if

        if((fv(sL)%IMEX_ID(eL) .le. 0) .eqv. is_exp) then
            fv(sL)%rhs(1:5,eL)  =  fv(sL)%rhs(1:5,eL)+rhsl(1:5,1)
            if(savD)    fv(sL)%duc(1:5,eL)  =  fv(sL)%duc(1:5,eL)+rhsD(1:5,1)
        end if
        if(((fv(sR)%IMEX_ID(eR) .le. 0) .eqv. is_exp) .and. sec(sR)%is_int) then
            fv(sR)%rhs(1:5,eR)  =  fv(sR)%rhs(1:5,eR)-rhsl(1:5,1)
            if(sec(sR)%is_int)  fv(sR)%duc(1:5,eR)  =  fv(sR)%duc(1:5,eR)-rhsD(1:5,1)
        end if
    end do

    return
    end subroutine fv_get_rhs_upw_imex
!   ----------------------------------------------------------------------------
!   cal Rhs, convective part, boundary.
!   ----------------------------------------------------------------------------
    subroutine fv_get_rhs_boundary_imex(lev,calD,addD,is_exp)
    use var_kind_def
    use var_cgns
    use var_fv
    use var_mesh
    use var_temp, only: uLf,uRf,vgf,fac_1d,rhsl,rhsD
    implicit none
    logical(dpL),intent(in):: calD,addD,is_exp
    integer(dpI),intent(in):: lev
    logical(dpL):: savD,is_upwind,is_cal
    integer(dpI):: im,sL,eL,sR,eR,bct

    savD=  calD .and. (.not. addD)
    do im=1,mesh(lev)%n_mortar_b
        sL  =  mesh(lev)%mortar_LR(1,im)
        eL  =  mesh(lev)%mortar_LR(2,im)
        sR  =  mesh(lev)%mortar_LR(3,im)
        eR  =  mesh(lev)%mortar_LR(4,im)
        fac_1d(1:4,1)   =  mesh(lev)%mortar_n_vg(1:4,im)
        vgf   (1:3,1)   =  mesh(lev)%mortar_n_vg(6:8,im)
        bct =  sec(sR)%bct

        is_cal  =((fv(sL)%IMEX_ID(eL) .le. 0) .eqv. is_exp)
        if(.not. is_cal)    cycle

        is_upwind   = (bct .eq. BCSymmetryPlane)
        if(is_upwind) then
            uLf(1:5,1)  =  fv(sR)%uL(1:5,eR)
            uRf(1:5,1)  =  fv(sR)%uR(1:5,eR)
            if(rhs_lev(lev) .eq. 3) then
                call rhs_conv_hllc(calD, addD, 1)
            elseif(rhs_lev(lev) .eq. 4) then
                call rhs_conv_ausm(calD, addD, 1)
            else
                call rhs_conv_roe (calD, addD, 1)
            end if
        else
            call u_to_F(1, fac_1d(1:3,1), fv(sR)%u(1,eR), vgf, rhsl)
            rhsl(1:5,1) =  rhsl(1:5,1)*fac_1d(4,1)
            rhsD(1:5,1) =  0.0d0
        end if

        if(is_cal) then
            fv(sL)%rhs(1:5,eL)  =  fv(sL)%rhs(1:5,eL)+rhsl(1:5,1)
            if(savD)    fv(sL)%duc(1:5,eL)  =  fv(sL)%duc(1:5,eL)+rhsD(1:5,1)
        end if
    end do

    return
    end subroutine fv_get_rhs_boundary_imex
!   ----------------------------------------------------------------------------
!   cal Rhs, viscous part, fv.
!   ----------------------------------------------------------------------------
    subroutine fv_get_rhs_vis_imp(lev,is_addD,is_exp)
    use var_kind_def
    use var_fv
    use var_mesh
    use var_temp, only: uLf,uRf,fac_1d,rhsD
    use var_turb, only: is_tur_cal,is_KO
    implicit none
    logical(dpL),intent(in):: is_addD,is_exp
    integer(dpI),intent(in):: lev
    logical(dpL):: is_cal
    integer(dpI):: im,sL,eL,sR,eR,j,k
    real   (dpR):: g(12),u(5),d(3),du(4),tL,tR,nne(3),dr,mu(2),tke(1),rtmp

    d   =  0.0d0
    mu  =  0.0d0
    do im=1,mesh(lev)%n_mortar_b
        sL  =  mesh(lev)%mortar_LR(1,im)
        eL  =  mesh(lev)%mortar_LR(2,im)
        sR  =  mesh(lev)%mortar_LR(3,im)
        eR  =  mesh(lev)%mortar_LR(4,im)
        is_cal  =((fv(sL)%IMEX_ID(eL) .le. 0) .eqv. is_exp)
        if(.not. is_cal)    cycle

        u  (1 :5 )  =  fv(sR)%u(1:5,eR)
        mu (1    )  =  fv(sR)%mu(1,eR)
        if(is_tur_cal)  mu(2)   =  fv(sR)%mu  (2,eR)
        if(is_KO     )  tke     =  fv(sR)%turb(1,eR)
        g  (1 :9 )  =  fv(sR)%gra(4 :12,eR)
        g  (10:12)  =  fv(sR)%gra(16:18,eR)

        call get_vis_fac(1, mesh(lev)%mortar_n_vg(1,im), u, mu, g, tke, rhsD)

        if(is_cal) then
            fv(sL)%rhs(1:5,eL)  =  fv(sL)%rhs(1:5,eL)-rhsD(1:5,1)
         else
            fv(sL)%duc(1:5,eL)  =  fv(sL)%duc(1:5,eL)+rhsD(1:5,1)
        end if
    end do

    do im=1+mesh(lev)%n_mortar_b,mesh(lev)%n_mortar
        sL  =  mesh(lev)%mortar_LR(1,im)
        eL  =  mesh(lev)%mortar_LR(2,im)
        sR  =  mesh(lev)%mortar_LR(3,im)
        eR  =  mesh(lev)%mortar_LR(4,im)
        fac_1d(1:4,1)   =  mesh(lev)%mortar_n_vg(1:4,im)

        if(sec(sR)%is_int) then
            is_cal  =((fv(sL)%IMEX_ID(eL) .le. 0) .eqv. is_exp) .or. &
                    &((fv(sR)%IMEX_ID(eR) .le. 0) .eqv. is_exp)
        else
            is_cal  =((fv(sL)%IMEX_ID(eL) .le. 0) .eqv. is_exp)
        end if
        if(.not. is_cal)    cycle

        uLf(1:5,1)  =  fv(sL)%u(1:5,eL)
        uRf(1:5,1)  =  fv(sR)%u(1:5,eR)
        u  (1:5  )  =  0.5d0*(uLf(1:5,1)+uRf(1:5,1))
        mu (1    )  = (fv(sL)%mu(1,eL)+fv(sR)%mu(1,eR))*0.5d0
        if(is_tur_cal)  mu(2)   = (fv(sL)%mu  (2,eL)+fv(sR)%mu  (2,eR))*0.5d0
        if(is_KO     )  tke     = (fv(sL)%turb(1,eL)+fv(sR)%turb(1,eR))*0.5d0
        tL          =  fv(sL)%t(    eL)
        tR          =  fv(sR)%t(    eR)
        g  (1 :9 )  = (fv(sL)%gra(4 :12,eL)+fv(sR)%gra(4 :12,eR))*0.5d0
        g  (10:12)  = (fv(sL)%gra(16:18,eL)+fv(sR)%gra(16:18,eR))*0.5d0

        d(1:3)      =  sec(sR)%cen(1:3,eR)-sec(sL)%cen(1:3,eL)
        rtmp        =  1.0d0/sqrt(d(1)**2+d(2)**2+d(3)**2)
        d           =  d*rtmp
        du(1:3)     =  rtmp*(uRf(2:4,1)-uLf(2:4,1))
        du(4  )     =  rtmp*(tR-tL)

        nne(1:3)=  fac_1d(1:3,1)/(fac_1d(1,1)*d(1)+fac_1d(2,1)*d(2)+fac_1d(3,1)*d(3))
        do j=1,4
            k       =  3*j-2
            dr      =  du(j)-g(k)*d(1)-g(k+1)*d(2)-g(k+2)*d(3)
            g(k:k+2)=  g(k:k+2)+dr*nne(1:3)
        end do

        call get_vis_fac(1, fac_1d, u, mu, g, tke, rhsD)

        if(is_addD) then
            if((fv(sL)%IMEX_ID(eL) .le. 0) .eqv. is_exp) then
                fv(sL)%rhs(1:5,eL)  =  fv(sL)%rhs(1:5,eL)-rhsD(1:5,1)
            end if
            if(((fv(sR)%IMEX_ID(eR) .le. 0) .eqv. is_exp) .and. sec(sR)%is_int) then
                fv(sR)%rhs(1:5,eR)  =  fv(sR)%rhs(1:5,eR)+rhsD(1:5,1)
            end if
        else
            if((fv(sL)%IMEX_ID(eL) .le. 0) .eqv. is_exp) then
                fv(sL)%duc(1:5,eL)  =  fv(sL)%duc(1:5,eL)+rhsD(1:5,1)
            end if
            if(((fv(sR)%IMEX_ID(eR) .le. 0) .eqv. is_exp) .and. sec(sR)%is_int) then
                fv(sR)%rhs(1:5,eR)  =  fv(sR)%rhs(1:5,eR)-rhsD(1:5,1)
            end if
        end if
    end do

    return
    end subroutine fv_get_rhs_vis_imp
!   ----------------------------------------------------------------------------
!   get LHS for the implicit part, IMEX.
!   ----------------------------------------------------------------------------
    subroutine fv_imex_get_lhs(lev,p)
    use var_kind_def
    use var_cgns
    use var_fv
    use var_mesh
    use var_prec
    use var_slv
    use var_turb, only: is_tur_cal
    use var_uns_cal, only: is_uns_cal_now,dt_uns,uns_iter,LHS_uns_c1
    implicit none
    integer(dpI),intent(in):: lev
    logical(dpL):: is_update_LHS
    integer(dpI):: isec,sL,eL,sR,eR,im,i,j,iele,bct,per(3),s_donor,e_donor
    real   (dpR):: JL(5,5),JR(5,5),JV(5,5),uL(5),uR(5),um(5),vg(3),spra,eig(3), &
                &  sprv,n(5),mul,eps,CFL,CFLu,R(5),uc(5),s,mut,JR_donor(5,5), &
                &  CFL_min
    type(type_prec):: p

    is_update_LHS   =  .true.
    if(.not. is_update_LHS) return

    vg  =  0.0d0
    mut =  0.0d0
    p%A =  0.0d0
    do im=1,mesh(lev)%n_mortar
        sL  =  mesh(lev)%mortar_LR(1,im)
        eL  =  mesh(lev)%mortar_LR(2,im)
        sR  =  mesh(lev)%mortar_LR(3,im)
        eR  =  mesh(lev)%mortar_LR(4,im)
        per =  0
        i   =  fv(sL)%IMEX_id(eL)
        j   =  0

        if(sec(sR)%is_int) then
            j   =  fv(sR)%IMEX_id(eR)
        elseif(sec(sR)%is_ghost) then
            call get_ID_ghost_ele(sR, eR, s_donor, e_donor, per)
            if(s_donor .gt. 0)  j   =  fv(s_donor)%IMEX_id(e_donor)
        else
!           nothing to do here.
        end if
        if((i .le. 0) .and. (j .le. 0)) cycle

        bct     =  sec(sR)%bct
        n (1:5) =  mesh(lev)%mortar_n_vg(1:5,im)
        vg(1:3) =  mesh(lev)%mortar_n_vg(6:8,im)

        uL(1:5) =  fv(sL)%u(1:5,eL)
        uR(1:5) =  fv(sR)%u(1:5,eR)

        call reigl_new(uL, vg, n, 1.0d0, 1.0d0, LHS_efix, JL, spra, eig)
        call reigl_new(uR, vg, n,-1.0d0, 1.0d0, LHS_efix, JR, spra, eig)

        if(is_vis_cal) then
            um  =  0.5d0*(uL+uR)
            mul =  0.5d0*(fv(sL)%mu(1,eL)+fv(sR)%mu(1,eR))
            if(is_tur_cal)  mut =  0.5d0*(fv(sL)%mu(2,eL)+fv(sR)%mu(2,eR))
            call reigl_vis(um, mul, mut, n, JV, sprv)
            JL  =  JL+JV
            JR  =  JR-JV
        end if

        if(i .gt. 0)    call prec_add_eleR(p, i, i, JL)
        if(sec(sR)%is_int .and. (j .gt. 0)) then
            call prec_add_eleR(p, i, j, JR)
            call prec_add_eleR(p, j, j,-JR)
            call prec_add_eleR(p, j, i,-JL)
        elseif(sec(sR)%is_ghost) then
            call get_ID_ghost_ele(sR, eR, s_donor, e_donor, per)
            if((s_donor .gt. 0) .and. (j .gt. 0)) then
                call get_J_donor(per, JR, JR_donor)
                call prec_add_eleR(p, i, j, JR_donor)
                call prec_add_eleR(p, j, j,-JR_donor)
                call prec_add_eleR(p, j, i,-JL      )
            end if
        end if
    end do

    JL  =  0.0d0
    eps =  1.0d-2
    CFL =  CFL_lev(lev)
    if(is_uns_cal_now .and. (uns_iter .ge. 100)) then
        CFL_min =  5.0d0
    else
        CFL_min =  1.0d0
    end if
    call fv_cal_LHS_scalar(lev, .false., .false.)
    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(sec(isec)%is_ghost .or. sec(isec)%is_bnd)    cycle

        do iele=1,sec(isec)%n_ele
            j       =  fv(isec)%IMEX_id(iele)
            if(j .le. 0)    cycle
            s       =  fv(isec)%LHS_s(  iele)
            uc(1:5) =  fv(isec)%uc (1:5,iele)
            R (1:5) =  fv(isec)%rhs(1:5,iele)
            CFLu    =  CFL
            CFLu    =  min(CFLu, eps*uc(1)*s/(abs(R(1))+1.0d-9))
            CFLu    =  min(CFLu, eps*uc(5)*s/(abs(R(5))+1.0d-9))
            CFLu    =  max(CFLu, CFL_min)

            do i=1,5
                JL(i,i) =  s/CFLu
                if(is_uns_cal_now)  JL(i,i) =  JL(i,i)+LHS_uns_c1*sec(isec)%vol(iele)/dt_uns
            end do

            fv(isec)%LHS_s(iele)=  s
            call prec_add_eleR(p, j, j, JL)
        end do
    end do

    call prec_get_decomposition(p)

    return
    end subroutine fv_imex_get_lhs
!   ----------------------------------------------------------------------------
!   transform Jacobian about ghost element to donor element.
!   ----------------------------------------------------------------------------
    subroutine get_J_donor(p,J,JJ)
    use var_kind_def
    use var_per_bnd
    implicit none
    integer(dpI),intent(in):: p(*)
    real   (dpR),intent(in):: J(5,*)
    integer(dpI):: N,iper
    real   (dpR):: JJ(5,*),T(5,5)

    call DCOPY(25, J, 1, JJ, 1)
    do iper=3,1,-1
        N   =  p(iper)
        if((N .le. 0) .or. (N .gt. n_per_info)) return
        if(nint(per_info(1,N)) .ne. 2)  cycle
        T(1:5,1)= (/1.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0/)
        T(1:5,2)= (/1.0d0, per_mat(1:3,N), 0.0d0/)
        T(1:5,3)= (/1.0d0, per_mat(4:6,N), 0.0d0/)
        T(1:5,4)= (/1.0d0, per_mat(7:9,N), 0.0d0/)
        T(1:5,5)= (/0.0d0, 0.0d0, 0.0d0, 0.0d0, 1.0d0/)
        JJ(1:5,1:5) =  matmul(JJ(1:5,1:5), T)
    end do

    return
    end subroutine get_J_donor
    end subroutine fv_slv_imex_region
