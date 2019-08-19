!-------------------------------------------------------------------------------
!   calculate Jacobian for the solid boundary.
!-------------------------------------------------------------------------------
    subroutine jac_bnd_sol(is_vis_wall,is_ghost,n,JL)
    use var_kind_def
    use var_global, only: I5
    implicit none
    logical(dpL),intent(in):: is_vis_wall,is_ghost
    real   (dpR):: n(*),JL(5,*)

    JL(1,1:5)   = (/1.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0/)
    if(is_vis_wall) then
        JL(2:4,1:5) =  0.0d0
    else
        JL(2,1:5)   = (/0.0d0, 1.0d0-n(1)*n(1), -n(1)*n(2), -n(1)*n(3), 0.0d0/)
        JL(3,1:5)   = (/0.0d0, -n(2)*n(1), 1.0d0-n(2)*n(2), -n(2)*n(3), 0.0d0/)
        JL(4,1:5)   = (/0.0d0, -n(3)*n(1), -n(3)*n(2), 1.0d0-n(3)*n(3), 0.0d0/)
    end if
    JL(5,1:5)   = (/0.0d0, 0.0d0, 0.0d0, 0.0d0, 1.0d0/)
    if(is_ghost)    JL(1:5,1:5) =  2.0d0*JL(1:5,1:5)-I5(1:5,1:5)

    return
    end subroutine jac_bnd_sol
!-------------------------------------------------------------------------------
!   setup the implicit solver, FV.
!-------------------------------------------------------------------------------
    subroutine fv_set_imp_slv(lev,p)
    use var_kind_def
    use var_eline
    use var_fv
    use var_global, only: err_mem,rref,uref
    use var_mesh
    use var_mg
    use var_prec
    use var_slv, solver=>solver_lev
    implicit none
    integer(dpI),intent(in):: lev
    integer(dpI):: isec,i,j,iA_sec(1000)
    real   (dpR):: scale_L(5),scale_R(5)
    type(type_prec):: p
    integer(dpI),allocatable:: iA(:),jA(:)

    if((solver(lev) .ne. solver_rk_sgs) .and. (solver(lev) .ne. solver_gmres_d))    return

    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(sec(isec)%is_ghost .or. sec(isec)%is_bnd)    cycle
        allocate(fv(isec)%LHS_s(sec(isec)%n_ele), stat=err_mem)
        allocate(fv(isec)%uc0(5,sec(isec)%n_ele), stat=err_mem)

        if(solver(lev) .eq. solver_rk_sgs) then
            allocate(fv(isec)%rv0(5,sec(isec)%n_ele), stat=err_mem)
        end if
    end do

    if(is_matrix_free .and. (solver(lev) .eq. solver_rk_sgs)) then
        call fv_set_prec(lev, p, 5, slv_SGS_MF, .true.)
    else
        if(prec_type .eq. 1) then
            call fv_set_prec(lev, p, 5, 1, .true.)
        elseif(prec_type .eq. 2) then
            call fv_set_prec(lev, p, 5, 3, .true.)
        end if
        scale_L =  1.0d0/(rref*(/uref, uref**2, uref**2, uref**2, uref**3/))
        scale_R = (/rref, rref*uref, rref*uref, rref*uref, rref*uref**2/)
!       call prec_set_scale(p, scale_L, scale_R)

        if((lev .eq. 0) .and. is_line_implicit .and. (n_eline .gt. 0)) then
            call mesh_get_1d_lines
            i           =  1
            iA_sec(1)   =  1
            do isec=mesh(0)%sec_1,mesh(0)%sec_0
                if(sec(isec)%is_int)    i   =  i+sec(isec)%n_ele
                iA_sec(isec+1)  =  i
            end do

            allocate(iA(n_eline+1), stat=err_mem)
            iA(1)   =  1
            do i=1,n_eline
                iA(i+1) =  iA(i)+eline(i)%n_ele
            end do
            allocate(jA(iA(n_eline+1)-1), stat=err_mem)
            do i=1,n_eline
            do j=1,eline(i)%n_ele
                isec            =  eline(i)%ele(1,i)
                jA(j+iA(i+1)-1) =  iA_sec(isec)+eline(i)%ele(2,i)-1
            end do
            end do
            call prec_set_lines(p, n_eline, iA, jA)
            if(allocated(iA))   deallocate(iA)
            if(allocated(jA))   deallocate(jA)
        end if
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
    end subroutine fv_set_imp_slv
!-------------------------------------------------------------------------------
!   setup the preconditioner, FV.
!-------------------------------------------------------------------------------
    subroutine fv_set_prec(lev,p,bsize,itype,is_include_per)
    use var_kind_def
    use var_fv
    use var_global, only: err_mem
    use var_mesh
    use var_prec
    implicit none
    logical(dpL),intent(in):: is_include_per
    integer(dpI),intent(in):: lev,bsize,itype
    integer(dpI):: isec,iele,nele,sR,eR,iA(1000),im,n_jA,ivtx,j,n_nb,nb(1000), &
                &  s_donor,e_donor,per(3)
    type(type_prec):: p

    if(p%is_spy)    return

    iA  =  1
    nele=  0
    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(.not. sec(isec)%is_int)  cycle
        iA(isec)=  nele+1
        nele    =  nele+sec(isec)%n_ele
    end do

    call prec_add(p, nele, bsize, itype)
    allocate(p%iA(nele+1), stat=err_mem)
    if(err_mem .ne. 0)  stop 'Error: fails to allocate memory, fv_set_prec.'
    p%iA(1) =  1

    n_jA=  0
    ivtx=  0
    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(.not. sec(isec)%is_int)  cycle
        do iele=1,sec(isec)%n_ele
            ivtx=  ivtx+1
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
                    n_jA=  n_jA+1
                elseif(sec(sR)%is_ghost .and. is_include_per) then
                    call get_ID_ghost_ele(sR, eR, s_donor, e_donor, per)
                    if(s_donor .gt. 0)  n_jA=  n_jA+1
                end if
            end do
            p%iA(ivtx+1)=  n_jA+1
        end do
    end do
    allocate(p%jA(n_jA), stat=err_mem)
    if(err_mem .ne. 0)  stop 'Error: fails to allocate memory, fv_set_prec.'

    n_jA    =  1
    ivtx    =  0
    p%iA(1) =  1
    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(.not. sec(isec)%is_int)  cycle
        do iele=1,sec(isec)%n_ele
            ivtx=  ivtx+1

            n_nb    =  1
            nb(1)   =  iele+iA(isec)-1
            do j=sec(isec)%iA_face_neighbour(iele),sec(isec)%iA_face_neighbour(iele+1)-1
                im  =  sec(isec)%jA_face_neighbour(j)
                if(im .gt. 0) then
                    sR  =  mesh(lev)%mortar_LR(3 ,im)
                    eR  =  mesh(lev)%mortar_LR(4 ,im)
                else
                    sR  =  mesh(lev)%mortar_LR(1,-im)
                    eR  =  mesh(lev)%mortar_LR(2,-im)
                end if

                if(sec(sR)%is_int) then
                    n_nb    =  n_nb+1
                    nb(n_nb)=  eR+iA(sR)-1
                elseif(sec(sR)%is_ghost .and. is_include_per) then
                    call get_ID_ghost_ele(sR, eR, s_donor, e_donor, per)
                    if(s_donor .le. 0)  cycle
                    n_nb    =  n_nb+1
                    nb(n_nb)=  e_donor+iA(s_donor)-1
                end if
            end do
            call simplify_series(n_nb, 1, 1, nb)

            call ICOPY(n_nb, nb, 1, p%jA(n_jA), 1)
            n_jA        =  n_jA+n_nb
            p%iA(ivtx+1)=  n_jA
        end do
    end do
    if(.not. allocated(p%A)) then
        allocate(p%A(p%bsize,p%bsize*(n_jA-1)), stat=err_mem)
        if(err_mem .ne. 0)  stop 'Error: fails to allocate memory, fv_set_prec.'
        p%A =  0.0d0
    end if

    p%is_spy=  .true.

    return
    end subroutine fv_set_prec
!-------------------------------------------------------------------------------
!   setup the implicit solver, FV.
!-------------------------------------------------------------------------------
    subroutine fv_get_imp_LHS(lev,p)
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
    integer(dpI):: isec,nele,sL,eL,sR,eR,iA(1000),im,i,iele,bct,per(3),s_donor,e_donor
    real   (dpR):: JL(5,5),JR(5,5),JV(5,5),uL(5),uR(5),um(5),vg(3),spra,eig(3), &
                &  sprv,n(5),mul,eps,CFL,CFLu,R(5),uc(5),s,mut,JR_donor(5,5), &
                &  CFL_min
    type(type_prec):: p
    character(20):: imp_LHS_imp = 'imp_LHS_imp'
    character(20):: imp_LHS_sca = 'imp_LHS_sca'
    character(20):: imp_LHS_decom = 'imp_LHS_decom'

    if(is_uns_cal_now) then
        is_update_LHS   =  mod(ite_wrk, 5) .eq. 1
    else
        is_update_LHS   =  .true.
    end if
    if(.not. is_update_LHS) return

    nele    =  0
    iA      =  1
    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(sec(isec)%is_ghost .or. sec(isec)%is_bnd)    cycle
        iA(isec)=  nele+1
        nele    =  nele+sec(isec)%n_ele
    end do

    vg  =  0.0d0
    mut =  0.0d0
    p%A =  0.0d0
    if((p%solver .eq. slv_SGS_MF) .and. (p%mlev .gt. 1))    p%clev%A=  0.0d0
call starttime(imp_LHS_imp)
    do im=1,mesh(lev)%n_mortar
        sL      =  mesh(lev)%mortar_LR(1  ,im)
        eL      =  mesh(lev)%mortar_LR(2  ,im)
        sR      =  mesh(lev)%mortar_LR(3  ,im)
        eR      =  mesh(lev)%mortar_LR(4  ,im)
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
!       if((sec(sR)%bct .eq. BCWallViscous) .or. (sec(sR)%bct .eq. BCWallInviscid)) then
!           call jac_bnd_sol(sec(sR)%bct .eq. BCWallViscous, .true., n, UR_UL)
!           JL  =  JL+matmul(JR, UR_UL)
!       end if

        eL  =  eL+iA(sL)-1
        eR  =  eR+iA(sR)-1
        call prec_add_eleR(p, eL, eL, JL)
        if(sec(sR)%is_int) then
            call prec_add_eleR(p, eL, eR, JR)
            call prec_add_eleR(p, eR, eR,-JR)
            call prec_add_eleR(p, eR, eL,-JL)
        elseif(sec(sR)%is_ghost) then
            call get_ID_ghost_ele(sR, eR, s_donor, e_donor, per)
            if(s_donor .gt. 0) then
                call get_J_donor(per, JR, JR_donor)
                eR  =  e_donor+iA(s_donor)-1
                call prec_add_eleR(p, eL, eR, JR_donor)
                call prec_add_eleR(p, eR, eR,-JR_donor)
                call prec_add_eleR(p, eR, eL,-JL      )
            end if
        end if
    end do
call endtime(imp_LHS_imp)

    JL  =  0.0d0
    eps =  1.0d-2
    CFL =  CFL_lev(lev)
    if(is_uns_cal_now .and. (uns_iter .ge. 100)) then
        CFL_min =  5.0d0
    else
        CFL_min =  1.0d1
    end if
call starttime(imp_LHS_sca)
    call fv_cal_LHS_scalar(lev, .false., .false.)
call endtime(imp_LHS_sca)
    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(sec(isec)%is_ghost .or. sec(isec)%is_bnd)    cycle

        do iele=1,sec(isec)%n_ele
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
            eL  =  iele+iA(isec)-1
            call prec_add_eleR(p, eL, eL, JL)
        end do
    end do

call starttime(imp_LHS_decom)
    call prec_get_decomposition(p)
call endtime(imp_LHS_decom)

    return
    contains
!       ------------------------------------------------------------------------
!       transform Jacobian about ghost element to donor element.
!       ------------------------------------------------------------------------
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
    end subroutine fv_get_imp_LHS
!-------------------------------------------------------------------------------
!   setup the implicit solver, FV. Sunway version.
!-------------------------------------------------------------------------------
    subroutine fv_get_imp_LHS_sw(lev,p)
    use var_kind_def
    use var_cgns
    use var_fv
    use var_fv_array
    use var_sec_array
    use var_global_real
    use var_air
    use var_mesh
    use var_prec
    use var_slv
    use var_turb, only: is_tur_cal
    use var_lhs_imp
    use var_uns_cal, only: is_uns_cal_now,dt_uns,uns_iter,LHS_uns_c1
    use var_parallel, only: myid
    implicit none
    integer(dpI),intent(in):: lev
    logical(dpL):: is_update_LHS
    integer(dpI):: isec,nele,sL,eL,sR,eR,iA(1000),im,i,iele,bct
    integer(dpI):: per(3),s_donor,e_donor,idxL,idxR,imi
    real   (dpR):: JL(5,5),JR(5,5),JV(5,5),uL(5),uR(5),um(5),vg(3),spra,eig(3), &
                &  sprv,n(5),mul,eps,CFL,CFLu,R(5),uc(5),s,mut,JR_donor(5,5), &
                &  CFL_min
    real   (dpR):: JL_tmp(25),JR_tmp(25),JV_tmp(25)
    type(type_prec):: p
    character(20):: imp_LHS_imp = 'imp_LHS_imp'
    character(20):: imp_LHS_sca = 'imp_LHS_sca'
    character(20):: imp_LHS_decom = 'imp_LHS_decom'

    if(is_uns_cal_now) then
        is_update_LHS   =  mod(ite_wrk, 5) .eq. 1
    else
        is_update_LHS   =  .true.
    end if
    if(.not. is_update_LHS) return

    nele    =  0
    iA      =  1
    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(sec(isec)%is_ghost .or. sec(isec)%is_bnd)    cycle
        iA(isec)=  nele+1
        nele    =  nele+sec(isec)%n_ele
    end do

    vg  =  0.0d0
    mut =  0.0d0
    p%A =  0.0d0
    if((p%solver .eq. slv_SGS_MF) .and. (p%mlev .gt. 1))    p%clev%A=  0.0d0

call fv_struct_to_array(lev)
    imi = 0
    do im=1+mesh(lev)%n_mortar_b,mesh(lev)%n_mortar
        imi = imi+1
        if(mesh_reordered(lev)%mortar_transform(imi) .eq. 0) then
            idxL = mesh_reordered(lev)%mortar_own_ID(imi)
            idxR = mesh_reordered(lev)%mortar_nei_ID(imi)
        else
            idxL = mesh_reordered(lev)%mortar_nei_ID(imi)
            idxR = mesh_reordered(lev)%mortar_own_ID(imi)
        end if
        bct     =  sec_bct(idxR)
        n (1:5) =  mesh_reordered(lev)%mortar_n_vg(1:5,imi)
        vg(1:3) =  mesh_reordered(lev)%mortar_n_vg(6:8,imi)

        uL(1:5) =  fv_u(1:5,idxL)
        uR(1:5) =  fv_u(1:5,idxR)

        call reigl_new_test(uL, vg, n, 1.0d0, 1.0d0, LHS_efix, JL_tmp, spra, eig)
        call reigl_new_test(uR, vg, n,-1.0d0, 1.0d0, LHS_efix, JR_tmp, spra, eig)

        if(is_vis_cal) then
            um  =  0.5d0*(uL+uR)
            mul =  0.5d0*(fv_mu(1,idxL)+fv_mu(1,idxR))
            if(is_tur_cal)  mut =  0.5d0*(fv_mu(2,idxL)+fv_mu(2,idxR))
            call reigl_vis_test(um, mul, mut, n, JV_tmp, sprv)
            JL_tmp  =  JL_tmp+JV_tmp
            JR_tmp  =  JR_tmp-JV_tmp
        end if
if(idxL .eq. 2020 .and. myid .eq. 0) then
    ! write(*,*),'reordered'
    ! write(*,*),fv_u(:,idxL)
    ! call print_matrix(5,JL_tmp)
end if

        diag(:,idxL) = diag(:,idxL)+JL_tmp(:)
        if(sec_is_int(idxR) .eq. 1) then
            ! diag(:,idxR) = diag(:,idxR)-JR_tmp(:)
            if(mesh_reordered(lev)%mortar_transform(imi) .eq. 0) then
                ! upper(:,imi) = upper(:,imi)+JR_tmp(:)
                ! lower(:,imi) = lower(:,imi)-JL_tmp(:)
            else
                ! upper(:,imi) = upper(:,imi)-JL_tmp(:)
                ! lower(:,imi) = lower(:,imi)+JR_tmp(:)
            end if
        end if
        fv_uc(1,idxL) = fv_uc(1,idxL)-JL_tmp(1)-JR_tmp(6)
        fv_uc(2,idxL) = fv_uc(2,idxL)-JL_tmp(2)-JR_tmp(7)
        fv_uc(3,idxL) = fv_uc(3,idxL)-JL_tmp(3)-JR_tmp(8)
        fv_uc(4,idxL) = fv_uc(4,idxL)-JL_tmp(4)-JR_tmp(9)
        fv_uc(5,idxL) = fv_uc(5,idxL)-JL_tmp(5)-JR_tmp(10)
    end do
    call imp_LHS_imp_host(mesh_reordered(lev)%owner, mesh_reordered(lev)%neighbor, &
        & faceNum, cellNum, mesh_reordered(lev)%mortar_transform, &
        & mesh_reordered(lev)%mortar_n_vg, sec_bct, fv_u, fv_mu, is_vis_cal_r, &
        & is_tur_cal_r, gk, gk1, LHS_efix, fv_uc, prl, prt, rr, cp)
call fv_array_to_struct(lev)


call starttime(imp_LHS_imp)
    do im=1+mesh(lev)%n_mortar_b,mesh(lev)%n_mortar
        sL      =  mesh(lev)%mortar_LR(1  ,im)
        eL      =  mesh(lev)%mortar_LR(2  ,im)
        sR      =  mesh(lev)%mortar_LR(3  ,im)
        eR      =  mesh(lev)%mortar_LR(4  ,im)
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
!       if((sec(sR)%bct .eq. BCWallViscous) .or. (sec(sR)%bct .eq. BCWallInviscid)) then
!           call jac_bnd_sol(sec(sR)%bct .eq. BCWallViscous, .true., n, UR_UL)
!           JL  =  JL+matmul(JR, UR_UL)
!       end if

        eL  =  eL+iA(sL)-1
        eR  =  eR+iA(sR)-1
if(eL .eq. 1 .and. myid .eq. 0) then
    ! write(*,*),'no-reorder'
    ! write(*,*),fv(sL)%u(:,eL)
    ! call print_matrix(5,JL)
end if
        call prec_add_eleR(p, eL, eL, JL)
        if(sec(sR)%is_int) then
            ! call prec_add_eleR(p, eL, eR, JR)
            ! call prec_add_eleR(p, eR, eR,-JR)
            ! call prec_add_eleR(p, eR, eL,-JL)
        elseif(sec(sR)%is_ghost) then
            call get_ID_ghost_ele(sR, eR, s_donor, e_donor, per)
            ! s_donor = sec(sR)%ID_ghost_ele(1,eR)
            ! e_donor = sec(sR)%ID_ghost_ele(2,eR)
            if(s_donor .gt. 0) then
                stop 'get_ID_ghost_ele has not been verified on Sunway platform'
                call get_J_donor(per, JR, JR_donor)
                eR  =  e_donor+iA(s_donor)-1
                call prec_add_eleR(p, eL, eR, JR_donor)
                call prec_add_eleR(p, eR, eR,-JR_donor)
                call prec_add_eleR(p, eR, eL,-JL      )
            end if
        end if
    end do
call endtime(imp_LHS_imp)

if(myid .eq. 0) call check_lhs(p,lev)
! stop

    do im=1,mesh(lev)%n_mortar_b
        sL      =  mesh(lev)%mortar_LR(1  ,im)
        eL      =  mesh(lev)%mortar_LR(2  ,im)
        sR      =  mesh(lev)%mortar_LR(3  ,im)
        eR      =  mesh(lev)%mortar_LR(4  ,im)
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
!       if((sec(sR)%bct .eq. BCWallViscous) .or. (sec(sR)%bct .eq. BCWallInviscid)) then
!           call jac_bnd_sol(sec(sR)%bct .eq. BCWallViscous, .true., n, UR_UL)
!           JL  =  JL+matmul(JR, UR_UL)
!       end if

        eL  =  eL+iA(sL)-1
        eR  =  eR+iA(sR)-1
        call prec_add_eleR(p, eL, eL, JL)
        if(sec(sR)%is_int) then
            call prec_add_eleR(p, eL, eR, JR)
            call prec_add_eleR(p, eR, eR,-JR)
            call prec_add_eleR(p, eR, eL,-JL)
        elseif(sec(sR)%is_ghost) then
            call get_ID_ghost_ele(sR, eR, s_donor, e_donor, per)
            ! s_donor = sec(sR)%ID_ghost_ele(1,eR)
            ! e_donor = sec(sR)%ID_ghost_ele(2,eR)
            if(s_donor .gt. 0) then
                stop 'get_ID_ghost_ele has not been verified on Sunway platform'
                call get_J_donor(per, JR, JR_donor)
                eR  =  e_donor+iA(s_donor)-1
                call prec_add_eleR(p, eL, eR, JR_donor)
                call prec_add_eleR(p, eR, eR,-JR_donor)
                call prec_add_eleR(p, eR, eL,-JL      )
            end if
        end if
    end do

    JL  =  0.0d0
    eps =  1.0d-2
    CFL =  CFL_lev(lev)
    if(is_uns_cal_now .and. (uns_iter .ge. 100)) then
        CFL_min =  5.0d0
    else
        CFL_min =  1.0d1
    end if
call starttime(imp_LHS_sca)
    call fv_cal_LHS_scalar(lev, .false., .false.)
call endtime(imp_LHS_sca)
    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(sec(isec)%is_ghost .or. sec(isec)%is_bnd)    cycle

        do iele=1,sec(isec)%n_ele
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
            eL  =  iele+iA(isec)-1
            call prec_add_eleR(p, eL, eL, JL)
        end do
    end do

call starttime(imp_LHS_decom)
    call prec_get_decomposition(p)
call endtime(imp_LHS_decom)

    return
    contains
!       ------------------------------------------------------------------------
!       transform Jacobian about ghost element to donor element.
!       ------------------------------------------------------------------------
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
    end subroutine fv_get_imp_LHS_sw
!-------------------------------------------------------------------------------
!   cal Rhs with dissipation lagging.
!-------------------------------------------------------------------------------
    subroutine fv_get_rhs_rk3(lev,RK)
    use var_kind_def
    use var_fv
    use var_mesh
    use var_slv, only: is_vis_cal,lev_out
    use var_turb, only: is_LES_now,LES_model
    use var_uns_cal, only: is_bdf_now
    use var_global, only: sw_slave,sw_time
    use var_parallel, only: myid
    implicit none
    integer(dpI),intent(in):: lev,RK
    integer(dpI):: isec
    integer(kind=8):: t1,t2,t3,t4
    character(20):: rhs_upw = 'rhs_upw'
    character(20):: rhs_vis = 'rhs_vis'

    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(sec(isec)%is_ghost .or. sec(isec)%is_bnd)    cycle
        fv(isec)%rhs=  0.0d0
        fv(isec)%duc=  0.0d0
    end do

    if(is_LES_now .and. (LES_model .ge. 1)) call fv_les(lev)
    if(is_shock_sensor) call fv_get_shock_sensor(lev)

    if(is_kep) then
        call fv_get_rhs_kep(lev)
    else
        if(rhs_lev(lev) .eq. 1) then
            if(is_jst_su2) then
                call fv_get_rhs_jst_su2(lev, .true., .false.)
            else
                call fv_get_rhs_jst    (lev, .true., .false.)
            end if
        else
call starttime(rhs_upw)
! if(sw_slave) then
    ! call fv_get_rhs_upw_sw(lev, .true., .false.)
! else 
    call fv_get_rhs_upw(lev, .true., .false.)
! end if
call endtime(rhs_upw)
        end if
    end if
    call fv_get_rhs_src_rotation(lev)
    if(is_bdf_now)  call fv_get_src_bdf(lev)

call starttime(rhs_vis)
! if(sw_slave) then
    ! if(is_vis_cal)  call fv_get_rhs_vis_sw(lev, .false.)
! else
    if(is_vis_cal)  call fv_get_rhs_vis(lev, .false.)
! end if
call endtime(rhs_vis)

    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(sec(isec)%is_ghost .or. sec(isec)%is_bnd)    cycle
        if(RK .eq. 1) then
            call DCOPY(5*sec(isec)%n_ele, fv(isec)%duc, 1, fv(isec)%rv0, 1)
        else
            fv(isec)%duc=  0.5d0*fv(isec)%rv0+0.5d0*fv(isec)%duc
        end if
        fv(isec)%rhs=  fv(isec)%rhs-fv(isec)%duc
        if(lev .gt. lev_out)    fv(isec)%rhs=  fv(isec)%rhs+fv(isec)%ff
    end do

    return
    end subroutine fv_get_rhs_rk3
!-------------------------------------------------------------------------------
!   cal Rhs with dissipation lagging. Sunway version.
!-------------------------------------------------------------------------------
    subroutine fv_get_rhs_rk3_sw(lev,RK)
    use var_kind_def
    use var_fv
    use var_fv_array
    use var_sec_array
    use var_mesh
    use var_slv, only: is_vis_cal,lev_out
    use var_turb, only: is_LES_now,LES_model
    use var_uns_cal, only: is_bdf_now
    use var_global, only: sw_slave,sw_time
    use var_parallel, only: myid
    implicit none
    integer(dpI),intent(in):: lev,RK
    integer(dpI):: isec,iele
    integer(kind=8):: t1,t2,t3,t4
    real   (dpR):: RK_r
    character(20):: rhs_upw = 'rhs_upw'
    character(20):: rhs_vis = 'rhs_vis'

    ! do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        ! if(sec(isec)%is_ghost .or. sec(isec)%is_bnd)    cycle
        ! fv(isec)%rhs=  0.0d0
        ! fv(isec)%duc=  0.0d0
    ! end do
    call rhs_rk3_zero_host(tot_ele, fv_duc, fv_rhs)

    if(is_LES_now .and. (LES_model .ge. 1)) stop 'fv_les has not been implemented on Sunway platform'
    if(is_shock_sensor) stop 'fv_get_shock_sensor has not been implemented on Sunway platform'

    if(is_kep) then
        stop 'rhs_kep has not been implemented on Sunway platform'
        call fv_get_rhs_kep(lev)
    else
        if(rhs_lev(lev) .eq. 1) then
            if(is_jst_su2) then
                stop 'rhs_jst_su2 has not been implemented on Sunway platform'
                call fv_get_rhs_jst_su2(lev, .true., .false.)
            else
                stop 'rhs_jst has not been implemented on Sunway platform'
                call fv_get_rhs_jst    (lev, .true., .false.)
            end if
        else
call starttime(rhs_upw)
    call fv_get_rhs_upw_sw(lev, .true., .false.)
call endtime(rhs_upw)
        end if
    end if
    call fv_get_rhs_src_rotation_sw(lev)
    if(is_bdf_now)  stop  'fv_get_src_bdf has not been implemented on Sunway platform'

call starttime(rhs_vis)
    if(is_vis_cal)  call fv_get_rhs_vis_sw(lev, .false.)
call endtime(rhs_vis)

    ! do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
    !     if(sec(isec)%is_ghost .or. sec(isec)%is_bnd)    cycle
    !     if(RK .eq. 1) then
    !         call DCOPY(5*sec(isec)%n_ele, fv(isec)%duc, 1, fv(isec)%rv0, 1)
    !     else
    !         fv(isec)%duc=  0.5d0*fv(isec)%rv0+0.5d0*fv(isec)%duc
    !     end if
    !     fv(isec)%rhs=  fv(isec)%rhs-fv(isec)%duc
    !     if(lev .gt. lev_out)    fv(isec)%rhs=  fv(isec)%rhs+fv(isec)%ff
    ! end do
    ! do iele=1,tot_ele
    !     ! if(sec_is_ghost(iele) .eq. 1 .or. sec_is_bnd(iele) .eq. 1) cycle
    !     if(RK .eq. 1) then
    !         fv_rv0(:,iele) = fv_duc(:,iele)
    !     else
    !         fv_duc(:,iele) = 0.5d0*fv_rv0(:,iele)+0.5d0*fv_duc(:,iele)
    !     end if
    !     fv_rhs(:,iele) = fv_rhs(:,iele)-fv_duc(:,iele)
    !     if(lev .gt. lev_out) stop 'ff has not been supported on Sunway platform'
    ! end do
    RK_r = real(RK, dpR)
    call rhs_rk3_upt_rhs_host(tot_ele, fv_rv0, fv_duc, fv_rhs, RK_r)

    return
    end subroutine fv_get_rhs_rk3_sw
!-------------------------------------------------------------------------------
!   RK3/implicit solver, FV.
!-------------------------------------------------------------------------------
    subroutine fv_slv_imp(lev)
    use var_kind_def
    use var_fv
    use var_jst_rk, only: RKgs3
    use var_mesh
    use var_prec
    use var_slv, only: RK,lev_out
    implicit none
    integer(dpI),intent(in):: lev
    integer(dpI):: isec
    real   (dpR):: CFL_RK3
    character(20):: slv_imp = 'slv_imp'
    character(20):: rhs_rk3 = 'rhs_rk3'
    character(20):: imp_lhs = 'imp_lhs'
    character(20):: iteration_t = 'iteration'
    character(20):: bnd_par = 'ban_par'

call starttime(slv_imp)
    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(sec(isec)%is_int)    fv(isec)%uc0=  fv(isec)%uc
    end do

    CFL_RK3 =  2.0d0
    do RK=1,3
call starttime(rhs_rk3)
        call fv_get_rhs_rk3(lev, RK)
call endtime(rhs_rk3)
        if((lev .eq. lev_out) .and. (RK .eq. 1))    call fv_get_residual(lev, res_NS)
call starttime(imp_lhs)
        if(RK .eq. 1)   call fv_get_imp_LHS(lev, fv_prec(lev))
call endtime(imp_lhs)

call starttime(iteration_t)
        call iteration(lev, fv_prec(lev))
call endtime(iteration_t)
        do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
            if(sec(isec)%is_int)    fv(isec)%duc=  CFL_RK3*RKgs3(RK)*fv(isec)%RHS
        end do
        call positivity_constraint
call starttime(bnd_par)
        call fv_bnd_parallel(lev, .true.)
call endtime(bnd_par)
    end do

call endtime(slv_imp)

    return
    contains
!       ------------------------------------------------------------------------
!       prec iteration.
!       ------------------------------------------------------------------------
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
            call DCOPY(5*sec(isec)%n_ele, p%x(i), 1, fv(isec)%rhs, 1)
            i   =  i+5*sec(isec)%n_ele
        end do

        return
        end subroutine iteration
!       ------------------------------------------------------------------------
!       positivity constraint.
!       ------------------------------------------------------------------------
        subroutine positivity_constraint
        implicit none
        logical(dpL):: is_positive
        integer(dpI):: iter,max_iter,isec,iele
        real   (dpR):: relax,uc(5)

        relax   =  1.0d0
        max_iter=  3
        do iter=1,max_iter
            is_positive =  .true.
            positivity:do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
                if(sec(isec)%is_ghost .or. sec(isec)%is_bnd)    cycle
                do iele=1,sec(isec)%n_ele
                    uc(1:5) =  fv(isec)%uc0(1:5,iele)+relax*fv(isec)%duc(1:5,iele)
                    is_positive = (uc(1) .gt. 0.0d0) .and. &
                                & (2.0d0*uc(1)*uc(5) .gt. uc(2)**2+uc(3)**2+uc(4)**2)
                    if(.not. is_positive)   exit positivity
                end do
            end do positivity

            if(is_positive) then
                do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
                    if(sec(isec)%is_int)    fv(isec)%uc =  fv(isec)%uc0+relax*fv(isec)%duc
                end do
                return
            else
                relax   =  relax*0.67d0
            end if
        end do
        stop 'Error: fails to satisfy the positivity constraint.'

        return
        end subroutine positivity_constraint
    end subroutine fv_slv_imp
!-------------------------------------------------------------------------------
!   RK3/implicit solver, FV. Sunway version
!-------------------------------------------------------------------------------
    subroutine fv_slv_imp_sw(lev)
    use var_kind_def
    use var_fv
    use var_jst_rk, only: RKgs3
    use var_mesh
    use var_prec
    use var_slv, only: RK,lev_out
    implicit none
    integer(dpI),intent(in):: lev
    integer(dpI):: isec
    real   (dpR):: CFL_RK3
    character(20):: slv_imp = 'slv_imp'
    character(20):: rhs_rk3 = 'rhs_rk3'
    character(20):: imp_lhs = 'imp_lhs'
    character(20):: iteration_t = 'iteration'
    character(20):: bnd_par = 'ban_par'

call starttime(slv_imp)
    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(sec(isec)%is_int)    fv(isec)%uc0=  fv(isec)%uc
    end do

    CFL_RK3 =  2.0d0
    do RK=1,3
call fv_struct_to_array(lev)
call starttime(rhs_rk3)
        call fv_get_rhs_rk3_sw(lev, RK)
call endtime(rhs_rk3)
        if((lev .eq. lev_out) .and. (RK .eq. 1))    call fv_get_residual_sw(lev, res_NS)
call fv_array_to_struct(lev)
call starttime(imp_lhs)
        if(RK .eq. 1)   call fv_get_imp_LHS_sw(lev, fv_prec(lev))
call endtime(imp_lhs)

call starttime(iteration_t)
        call iteration(lev, fv_prec(lev))
call endtime(iteration_t)
        do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
            if(sec(isec)%is_int)    fv(isec)%duc=  CFL_RK3*RKgs3(RK)*fv(isec)%RHS
        end do
call fv_struct_to_array(lev)
        call positivity_constraint_sw
call fv_array_to_struct(lev)
call fv_struct_to_array(lev)
call starttime(bnd_par)
        call fv_bnd_parallel_sw(lev, .true.)
call endtime(bnd_par)
call fv_array_to_struct(lev)
    end do

call endtime(slv_imp)

    return
    contains
!       ------------------------------------------------------------------------
!       prec iteration.
!       ------------------------------------------------------------------------
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
            call DCOPY(5*sec(isec)%n_ele, p%x(i), 1, fv(isec)%rhs, 1)
            i   =  i+5*sec(isec)%n_ele
        end do

        return
        end subroutine iteration
!       ------------------------------------------------------------------------
!       positivity constraint. Sunway version
!       ------------------------------------------------------------------------
        subroutine positivity_constraint_sw
        use var_fv_array
        use var_sec_array
        implicit none
        logical(dpL):: is_positive
        integer(dpI):: iter,max_iter,isec,iele
        real   (dpR):: relax,uc(5)

        relax   =  1.0d0
        max_iter=  3
        do iter=1,max_iter
            is_positive =  .true.
            positivity:do iele=1,tot_ele
                if(sec_is_int(iele) .eq. 1) then
                    uc(1:5) = fv_uc0(1:5,iele)+relax*fv_duc(1:5,iele)
                    is_positive = (uc(1) .gt. 0.0d0) .and. &
                                & (2.0d0*uc(1)*uc(5) .gt. uc(2)**2+uc(3)**2+uc(4)**2)
                    if(.not. is_positive) exit positivity
                end if
            end do positivity
            if(is_positive) then
                do iele=1,tot_ele
                    if(sec_is_int(iele) .eq. 1) then
                        fv_uc(:,iele) = fv_uc0(:,iele)+relax*fv_duc(:,iele)
                    end if
                end do
                return
            else
                relax = relax*0.67d0
            end if
        end do
        stop 'Error: fails to satisfy the positivity constraint.'

        return
        end subroutine positivity_constraint_sw
    end subroutine fv_slv_imp_sw
!-------------------------------------------------------------------------------
!   RK3/implicit solver, FV.
!-------------------------------------------------------------------------------
    subroutine fv_slv_rk3_sgs(lev)
    use var_kind_def
    use var_fv
    use var_jst_rk, only: RKgs2,RKgs3
    use var_mesh
    use var_slv, only: is_reorder,is_vis_cal,RK,lev_out
    implicit none
    integer(dpI),parameter:: RK_max =  3
    integer(dpI),intent(in):: lev
    logical(dpL):: is_exchange_duc
    integer(dpI):: isec
    real   (dpR):: CFL_RK,RKgs(3)

    is_exchange_duc =  .true.

    if(RK_max .eq. 3) then
        CFL_RK      =  2.0d0
        RKgs(1:3)   =  RKgs3(1:3)
    else
        CFL_RK      =  1.0d0
        RKgs(1:2)   =  RKgs2(1:2)
    end if

    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(sec(isec)%is_int)    fv(isec)%uc0=  fv(isec)%uc
    end do

    do RK=1,RK_max
        call fv_get_rhs_rk3(lev, RK)

        if((lev .eq. lev_out) .and. (RK .eq. 1))    call fv_get_residual(lev, res_NS)
        if(RK .eq. 1)   call fv_get_imp_LHS(lev, fv_prec(lev))

        if(is_reorder) then
            call iteration_rcm(lev, fv_prec(lev))
        else
            call iteration    (lev, fv_prec(lev))
        end if
        do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
            if(sec(isec)%is_int)    fv(isec)%duc=  CFL_RK*RKgs(RK)*fv(isec)%duc
        end do
        call add_duc
        call fv_bnd_parallel(lev, .true.)
    end do

    return
    contains
!       ------------------------------------------------------------------------
!       prec iteration.
!       ------------------------------------------------------------------------
        subroutine iteration(lev,p)
        use var_prec, only: prec_solve
        use var_turb, only: is_tur_cal
        implicit none
        integer(dpI),intent(in):: lev
        type(type_prec),intent(in):: p
        integer(dpI):: max_iter,iter,isec,iele,n_nb,i_nb,j_nb,s_nb,e_nb,im,iA(1000)
        real   (dpR):: vg(3,100),e(5,100),uR(5,100),d(5,100),mu_R(2,100),akh_R(3,100), &
                    &  sgn(100),efix(100),v(5,100),r(5),uL(5),mu_L(2)

        max_iter=  1
        iA(1)   =  1
        do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
            iA(isec+1)  =  iA(isec)+sec(isec)%n_ele
            fv(isec)%duc=  0.0d0
        end do
        do iter=1,max_iter
!           Forward iteration.
            do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
                if(.not. sec(isec)%is_int)  cycle

                do iele=1,sec(isec)%n_ele
                    include './include/fv_axd.f90'
                end do
            end do

!           exchange duc.
            if(is_exchange_duc) call exchange_duc

!           Backward iteration.
            do isec=mesh(lev)%sec_0,mesh(lev)%sec_1,-1
                if(.not. sec(isec)%is_int)  cycle

                do iele=sec(isec)%n_ele,1,-1
                    include './include/fv_axd.f90'
                end do
            end do

!           exchange duc.
            if(is_exchange_duc .and. (iter .lt. max_iter))  call exchange_duc
        end do

        if(p%mlev .le. 1)   return

        call get_bAx(lev, p)
        call prec_solve(p, .true., 1)

        iele=  1
        do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
            if(.not. sec(isec)%is_int)  cycle
            call DCOPY(5*sec(isec)%n_ele, p%x(iele), 1, fv(isec)%duc, 1)
            iele=  iele+5*sec(isec)%n_ele
        end do

        return
        end subroutine iteration
!       ------------------------------------------------------------------------
!       prec iteration with RCM ordering.
!       ------------------------------------------------------------------------
        subroutine iteration_rcm(lev,p)
        use var_turb, only: is_tur_cal
        implicit none
        integer(dpI),intent(in):: lev
        type(type_prec),intent(in):: p
        integer(dpI):: max_iter,iter,isec,iele,n_nb,i_nb,j_nb,s_nb,e_nb,im,i,iA(1000)
        real   (dpR):: vg(3,100),e(5,100),uR(5,100),d(5,100),mu_R(2,100),akh_R(3,100), &
                    &  sgn(100),efix(100),v(5,100),r(5),uL(5),mu_L(2)

        iA(1)   =  1
        do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
            iA(isec+1)  =  iA(isec)+sec(isec)%n_ele
            fv(isec)%duc=  0.0d0
        end do

        max_iter=  2
        do iter=1,max_iter
!           Forward iteration.
            do i=1,p%nvtx
                isec=  mesh(lev)%rcm_order(1,i)
                iele=  mesh(lev)%rcm_order(2,i)

                include './include/fv_axd.f90'
            end do

!           exchange duc.
            if(is_exchange_duc) call exchange_duc

!           Backward iteration.
            do i=p%nvtx,1,-1
                isec=  mesh(lev)%rcm_order(1,i)
                iele=  mesh(lev)%rcm_order(2,i)

                include './include/fv_axd.f90'
            end do

!           exchange duc.
            if(is_exchange_duc .and. (iter .lt. max_iter))  call exchange_duc
        end do

        return
        end subroutine iteration_rcm
!       ------------------------------------------------------------------------
!       get b-Ax.
!       ------------------------------------------------------------------------
        subroutine get_bAx(lev,p)
        use var_turb, only: is_tur_cal
        implicit none
        integer(dpI),intent(in):: lev
        integer(dpI):: isec,iele,n_nb,i_nb,j_nb,s_nb,e_nb,im,iA(1000)
        real   (dpR):: vg(3,100),e(5,100),uR(5,100),d(5,100),mu_R(2,100),akh_R(3,100), &
                    &  sgn(100),efix(100),v(5,100),r(5),uL(5),mu_L(2)
        type(type_prec):: p

        iA(1)   =  1
        do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
            iA(isec+1)  =  iA(isec)+sec(isec)%n_ele
        end do

        do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(.not. sec(isec)%is_int)  cycle

        do iele=1,sec(isec)%n_ele
        n_nb=  sec(isec)%iA_face_neighbour(iele+1)-sec(isec)%iA_face_neighbour(iele)
        do i_nb=sec(isec)%iA_face_neighbour(iele),sec(isec)%iA_face_neighbour(iele+1)-1
            j_nb=  i_nb-sec(isec)%iA_face_neighbour(iele)+1
            im  =  sec(isec)%jA_face_neighbour(i_nb)
            if(im .gt. 0) then
                s_nb        =  mesh(lev)%mortar_LR(3  ,im)
                e_nb        =  mesh(lev)%mortar_LR(4  ,im)
                e (1:5,j_nb)=  mesh(lev)%mortar_n_vg(1:5,im)
                vg(1:3,j_nb)=  mesh(lev)%mortar_n_vg(6:8,im)
            else
                s_nb        =  mesh(lev)%mortar_LR(1  ,-im)
                e_nb        =  mesh(lev)%mortar_LR(2  ,-im)
                e (1:3,j_nb)= -mesh(lev)%mortar_n_vg(1:3,-im)
                e (4:5,j_nb)=  mesh(lev)%mortar_n_vg(4:5,-im)
                vg(1:3,j_nb)=  mesh(lev)%mortar_n_vg(6:8,-im)
            end if
            uR(1:5,j_nb)=  fv(s_nb)%u  (1:5,e_nb)
            d (1:5,j_nb)=  fv(s_nb)%duc(1:5,e_nb)
            if(is_vis_cal)  mu_R(1,j_nb)=  fv(s_nb)%mu(1,e_nb)
            if(is_tur_cal)  mu_R(2,j_nb)=  fv(s_nb)%mu(2,e_nb)
        end do
        call u_to_akh(n_nb, uR, akh_R)

        sgn (1:n_nb)= -1.0d0
        efix(1:n_nb)=  LHS_efix
        call AprdV(n_nb, e, vg, uR, akh_R, sgn, efix, d, v)
        r(1:5)  = -fv(isec)%rhs(1:5,iele)
        do j_nb=1,n_nb
            r(1:5)  =  r(1:5)-v(1:5,j_nb)
        end do

        if(is_vis_cal) then
            uL(1:5) =  fv(isec)%u (1:5,iele)
            mu_L(1) =  fv(isec)%mu(1  ,iele)
            if(is_tur_cal)  mu_L(2) =  fv(isec)%mu(2,iele)

            call AprdV_vis(n_nb, e, uL, uR, mu_L, mu_R, d, v)
            do j_nb=1,n_nb
                r(2:5)  =  r(2:5)+v(2:5,j_nb)
            end do
        end if

        n_nb=  iele+iA(isec)-1
        p%x  (5*n_nb-4:5*n_nb)  =  fv(isec)%duc(1:5,iele)
        p%RHS(5*n_nb-4:5*n_nb)  =  r(1:5) &
                                & -p%A(1:5,5*n_nb-4)*fv(isec)%duc(1,iele) &
                                & -p%A(1:5,5*n_nb-3)*fv(isec)%duc(2,iele) &
                                & -p%A(1:5,5*n_nb-2)*fv(isec)%duc(3,iele) &
                                & -p%A(1:5,5*n_nb-1)*fv(isec)%duc(4,iele) &
                                & -p%A(1:5,5*n_nb  )*fv(isec)%duc(5,iele)
        end do
        end do

        return
        end subroutine get_bAx
!       ------------------------------------------------------------------------
!       add duc.
!       ------------------------------------------------------------------------
        subroutine add_duc
        use var_fv
        implicit none
        integer(dpI):: isec,iele
        real   (dpR):: r1,r5,r,eps

        eps =  5.0d-2
        do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
            if(.not. sec(isec)%is_int)  cycle
            do iele=1,sec(isec)%n_ele
                r1  =  abs(fv(isec)%duc(1,iele))/fv(isec)%uc0(1,iele)
                r5  =  abs(fv(isec)%duc(5,iele))/fv(isec)%uc0(5,iele)
                r   =  min(1.0d0, eps/max(r1, r5))
                fv(isec)%uc(1:5,iele)   =  fv(isec)%uc0(1:5,iele) &
                                        & +r*fv(isec)%duc(1:5,iele)
            end do
        end do

        return
        end subroutine add_duc
!       ------------------------------------------------------------------------
!       exchange duc in the parallel environment.
!       ------------------------------------------------------------------------
        subroutine exchange_duc
        use var_kind_def
        use var_fv
        use var_mesh
        use var_parallel
        implicit none
        integer(dpI):: isr,isec,iele,i,idx,ip_remote

!       ------------------------------------------------------------------------
!       prepare data.
        do isr=mesh(lev)%p2p_1,mesh(lev)%p2p_0
            idx =  1
            do i=1,p2p(isr)%n_ele_send
                isec=  p2p(isr)%id_ele_send(1,i)
                iele=  p2p(isr)%id_ele_send(2,i)
                call DCOPY(5, fv(isec)%duc(1,iele), 1, p2p(isr)%rsend(idx), 1)
                idx =  idx+5
            end do
            p2p(isr)%n_send =  idx-1
            p2p(isr)%n_recv =  5*p2p(isr)%n_ele_recv
        end do
!       prepare data.
!       ------------------------------------------------------------------------

!       ------------------------------------------------------------------------
!       data exchange.
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
!       data exchange.
!       ------------------------------------------------------------------------

!       ------------------------------------------------------------------------
!       data exchange on myid.
        do isr=mesh(lev)%p2p_1,mesh(lev)%p2p_0
            ip_remote   =  p2p(isr)%ip_remote
            if(myid .ne. ip_remote) cycle
            if(p2p(isr)%n_send .ne. p2p(isr)%n_recv)    stop 'Error: s&r not match on myid.'

            if(p2p(isr)%n_send .gt. 0)  call DCOPY(p2p(isr)%n_send, p2p(isr)%rsend, 1, &
                &  p2p(isr)%rrecv, 1)
        end do
!       data exchange on myid.
!       ------------------------------------------------------------------------

        call mpi_waitall(mpi_nreq, mpi_req, mpi_sta, mpi_err)

!       ------------------------------------------------------------------------
!       unpack data received.
        do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
            if(.not. sec(isec)%is_ghost)    cycle
            do i=1,sec(isec)%n_ele
                isr =  sec(isec)%id_recv(1,i)
                iele=  sec(isec)%id_recv(2,i)
                call per_rot_vec(21, sec(isec)%per_path(1,i), p2p(isr)%rrecv(5*iele-4), &
                    &  fv(isec)%duc(1,i))
            end do
        end do
!       unpack data received.
!       ------------------------------------------------------------------------

        return
        end subroutine exchange_duc
    end subroutine fv_slv_rk3_sgs

    subroutine print_matrix(N,A)
    use var_kind_def
    use var_mesh
    implicit none

    real   (dpR),intent(IN):: A(*)
    integer(dpI),intent(IN):: N
    integer(dpI):: i,j

    do i=1,N
        do j=1,N
            write(*,'(f12.6,$)'),A((i-1)*N+j)
        end do
        write(*,*), ''
    end do
    end subroutine print_matrix

    function compare_float(a,b)
    use var_kind_def
    implicit none
    real   (dpR),intent(IN):: a,b
    real   (dpR),parameter:: eps = 1.0d-14
    logical(dpL):: compare_float

    compare_float = .true.
    if(a .eq. 0) then
        if(b .ne. 0) compare_float = .false.
    else if(abs(a-b) .gt. eps) then
        if(abs(a-b)/a .gt. eps) compare_float = .false.
    end if

    return
    end function

    subroutine check_matrix(N,A,B)
    use var_kind_def
    use var_mesh
    implicit none

    real   (dpR),intent(IN):: A(*),B(*)
    integer(dpI),intent(IN):: N
    integer(dpI):: i,j
    logical(dpL):: ltmp
    logical(dpL),external :: compare_float

    ltmp = .true.
    do i=1,N
        do j=1,N
            if(.not. compare_float(A((i-1)*N+j),B((i-1)*N+j))) ltmp = .false.
            if(.not. compare_float(A((i-1)*N+j),B((i-1)*N+j))) write(*,*),A((i-1)*N+j),B((i-1)*N+j)
        end do
    end do

    if(.not. ltmp) then
        write(*,*),'!!!----Error------------------'
        call print_matrix(5,A)
        call print_matrix(5,B)
        stop
    end if
    end subroutine check_matrix

    subroutine check_lhs(p,lev)
    use var_kind_def
    use var_mesh
    use var_prec
    use var_slv
    use var_fv_array
    use var_lhs_imp
    implicit none

    integer(dpI),intent(IN):: lev
    integer(dpI):: isec,iele,iA(1000),nele,idxL,imi,im,k,idxR,idx
    integer(dpI):: idxL0,idxR0,sL,eL,sR,eR
    type(type_prec):: p

    nele    =  0
    iA      =  1
    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(sec(isec)%is_ghost .or. sec(isec)%is_bnd)    cycle
        iA(isec)=  nele+1
        nele    =  nele+sec(isec)%n_ele
    end do
    imi = 0
    do im=1+mesh(lev)%n_mortar_b,mesh(lev)%n_mortar
        imi = imi+1
        sL      =  mesh(lev)%mortar_LR(1  ,im)
        eL      =  mesh(lev)%mortar_LR(2  ,im)
        sR      =  mesh(lev)%mortar_LR(3  ,im)
        eR      =  mesh(lev)%mortar_LR(4  ,im)

        idxL    =  sec(sL)%ID_ele_g(eL)
        idxR    =  sec(sR)%ID_ele_g(eR)

        if(idxL .le. cellNum) idxL = perm(idxL)
        if(idxR .le. cellNum) idxR = perm(idxR)

        eL  =  eL+iA(sL)-1
        eR  =  eR+iA(sR)-1

        if(p%is_reordered) then
            idxL0  =  p%perm(eL)
        else
            idxL0  =  eL
        end if
        if(p%is_spy) then
            do k=p%iA(idxL0),p%iA(idxL0+1)-1
                if(p%jA(k) .eq. idxL0) then
                    write(*,*),'left',eL,eR
                    call check_matrix(p%bsize,p%A(1,1+p%bsize*(k-1)),diag(1,idxL))
                    ! write(*,*),'right'
                    ! call check_matrix(p%bsize,p%A(1,1+p%bsize*(k-1)),diag(1,idxR))
                    ! call print_matrix(p%bsize,p%A(1,1+p%bsize*(k-1)))
                    ! call print_matrix(p%bsize,diag(1,idxL))
                    ! call print_matrix(p%bsize,diag(1,idxR))
                    ! stop
                end if
            end do
        else
            stop 'is_spy has not been implemented on Sunway platform'
        end if
    end do
        ! diag(:,idxL)
        ! if(sec_is_int(idxR) .eq. 1) then
        !     diag(:,idxR) = diag(:,idxR)-JR_tmp(:)
        !     if(mesh_reordered(lev)%mortar_transform(imi) .eq. 0) then
        !         upper(:,imi) = upper(:,imi)+JR_tmp(:)
        !         lower(:,imi) = lower(:,imi)-JL_tmp(:)
        !     else
        !         upper(:,imi) = upper(:,imi)-JL_tmp(:)
        !         lower(:,imi) = lower(:,imi)+JR_tmp(:)
        !     end if
        ! end if
    end subroutine check_lhs
