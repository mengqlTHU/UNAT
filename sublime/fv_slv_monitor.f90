!-------------------------------------------------------------------------------
!   calculate total quantities, perform integration.
!-------------------------------------------------------------------------------
    subroutine fv_integrate_force_mass
    use var_kind_def
    use var_bndv
    use var_parallel
    implicit none
    integer(dpI):: i
    real   (dpR):: inflow(100),outflow(100)

    if(n_inl .gt. 0) then
        call mpi_allreduce(mass_inl, inflow , n_inl, mpi_dpR, mpi_sum, &
            &  mpi_comm_world, mpi_err)
        mass_inl(1:n_inl)   =  inflow(1:n_inl)
    end if

    if(n_out .gt. 0) then
        call mpi_allreduce(mass_out, outflow, n_out, mpi_dpR, mpi_sum, &
            &  mpi_comm_world, mpi_err)
        mass_out(1:n_out)   =  outflow(1:n_out)
    end if

    if(myid .eq. 0) then
        do i=1,n_inl
!           print*,mass_inl(i)
        end do
    end if

    return
    end subroutine fv_integrate_force_mass
!-------------------------------------------------------------------------------
!   calculate total quantities, perform integration.
!-------------------------------------------------------------------------------
    subroutine fv_get_force_mass(lev)
    use var_kind_def
    use var_air, only: rr,t_Suth
    use var_bndv, only: u_fs,ua_fs,Area_ref,mass_inl,mass_out
    use var_cgns
    use var_fv
    use var_mesh, only: sec,mesh
    use var_slv
    implicit none
    integer(dpI),intent(in):: lev
    integer(dpI):: im,bct,sL,eL,sR,eR,ID
    real   (dpR):: A,u(5),n(3),g(18),gra(9),t,mu,nu(6),f(3),dm,ua

    mfin    =  0.0d0
    mfot    =  0.0d0
    drag_p  =  0.0d0
    lift_p  =  0.0d0
    drag_v  =  0.0d0
    lift_v  =  0.0d0
    mass_inl=  0.0d0
    mass_out=  0.0d0
    do im=1,mesh(lev)%n_mortar_b
        sR      =  mesh(lev)%mortar_LR(3,im)
        eR      =  mesh(lev)%mortar_LR(4,im)
        n(1:3)  =  mesh(lev)%mortar_n_vg(1:3,im)
        A       =  mesh(lev)%mortar_n_vg(4  ,im)
        u(1:5)  =  fv(sR)%u(1:5,eR)
        if(is_vis_cal)  g(1:18) =  fv(sR)%gra(1:18,eR)
        bct     =  sec(sR)%bct
        ID      =  sec(sR)%ID_group

        if(bct .eq. BCInflow) then
            dm          = -A*u(1)*(u(2)*n(1)+u(3)*n(2)+u(4)*n(3))
            mfin        =  mfin+dm
            mass_inl(ID)=  mass_inl(ID)+dm
        elseif(bct .eq. BCInflowSupersonic) then
            dm          = -A*u(1)*(u(2)*n(1)+u(3)*n(2)+u(4)*n(3))
            mfin        =  mfin+dm
            mass_inl(ID)=  mass_inl(ID)+dm
        elseif(bct .eq. BCOutflow) then
            dm          =  A*u(1)*(u(2)*n(1)+u(3)*n(2)+u(4)*n(3))
            mfot        =  mfot+dm
            mass_out(ID)=  mass_out(ID)+dm
        elseif(bct .eq. BCOutflowSupersonic) then
            dm          =  A*u(1)*(u(2)*n(1)+u(3)*n(2)+u(4)*n(3))
            mfot        =  mfot+dm
            mass_out(ID)=  mass_out(ID)+dm
        end if

        if((bct .eq. BCWallInviscid) .or. (bct .eq. BCWallViscous)) then
            f       =-(u(5)-u_fs(5))*n*A
            drag_p  =  drag_p+f(1)*dir_drag(1)+f(2)*dir_drag(2)+f(3)*dir_drag(3)
            lift_p  =  lift_p+f(1)*dir_lift(1)+f(2)*dir_lift(2)+f(3)*dir_lift(3)
        end if

        if(bct .eq. BCWallViscous) then
            gra(1:9)=  g(4:12)
            dm      = -2.0d0/3.0d0*(gra(1)+gra(5)+gra(9))
            nu(1)   =  gra(1)*2.0d0+dm
            nu(2)   =  gra(2)+gra(4)
            nu(3)   =  gra(3)+gra(7)
            nu(4)   =  gra(5)*2.0d0+dm
            nu(5)   =  gra(6)+gra(8)
            nu(6)   =  gra(9)*2.0d0+dm
            t       =  u(5)/(rr*u(1))
            mu      =  1.461d-6*sqrt(t**3)/(t+t_suth)
            nu      =  mu*nu
            f(1)    =  A*(nu(1)*n(1)+nu(2)*n(2)+nu(3)*n(3))
            f(2)    =  A*(nu(2)*n(1)+nu(4)*n(2)+nu(5)*n(3))
            f(3)    =  A*(nu(3)*n(1)+nu(5)*n(2)+nu(6)*n(3))
            drag_v  =  drag_v+f(1)*dir_drag(1)+f(2)*dir_drag(2)+f(3)*dir_drag(3)
            lift_v  =  lift_v+f(1)*dir_lift(1)+f(2)*dir_lift(2)+f(3)*dir_lift(3)
        end if
    end do

    A       =  1.0d0/(0.5d0*u_fs(1)*ua_fs*ua_fs*Area_ref)
    lift_p  = -lift_p*A
    lift_v  = -lift_v*A
    drag_p  = -drag_p*A
    drag_v  = -drag_v*A
    lift    =  lift_p+lift_v
    drag    =  drag_p+drag_v

!   get cp and cf for solid wall.
    dm  =  1.0d0/(0.5d0*u_fs(1)*ua_fs**2)
    do im=1,mesh(0)%n_mortar_b
        sL  =  mesh(0)%mortar_LR(1,im)
        eL  =  mesh(0)%mortar_LR(2,im)
        sR  =  mesh(0)%mortar_LR(3,im)
        eR  =  mesh(0)%mortar_LR(4,im)
        if(sec(sR)%bct .ne. BCWallViscous)  cycle

        ua  =  sqrt(fv(sL)%u(2,eL)**2+fv(sL)%u(3,eL)**2+fv(sL)%u(4,eL)**2)
        fv(sR)%bnd_solution(1,eR)   = (fv(sR)%u(5,eR)-u_fs(5))*dm
        fv(sR)%bnd_solution(2,eR)   =  fv(sR)%mu(1,eR)*ua/sec(sL)%dnw(eL)*dm
    end do

    return
    end subroutine fv_get_force_mass
!-------------------------------------------------------------------------------
!   monitor the solution process.
!-------------------------------------------------------------------------------
    subroutine fv_slv_monitor(is_end,lev)
    use var_kind_def
    use var_fv
    use var_global
    use var_parallel
    use var_slv
    use var_turb, only: is_RANS,RANS_model,RANS_WA,RANS_SST
    use var_uns_cal
    implicit none
    logical(dpL),intent(in):: is_end
    integer(dpI),intent(in):: lev
    character(len=1000):: tec
    logical(dpL):: ltmp
    integer(dpI):: i,iter
    real   (dpR):: s(10),r(10),cost(nprc)

    call fv_get_force_mass(lev)
    s(1 )   =  res_NS
    s(2 )   =  res_RANS
    s(3 )   =  lift
    s(4 )   =  drag
    s(5 )   =  lift_p
    s(6 )   =  lift_v
    s(7 )   =  drag_p
    s(8 )   =  drag_v
    s(9 )   =  mfin
    s(10)   =  mfot
    call mpi_allreduce(s, r, 10, mpi_dpR, mpi_sum, mpi_comm_world, mpi_err)
    res_NS  =  sqrt(r(1))
    res_RANS=  sqrt(r(2))
    lift    =  r(3)
    drag    =  r(4)
    lift_p  =  r(5)
    lift_v  =  r(6)
    drag_p  =  r(7)
    drag_v  =  r(8)
    mfin    =  r(9)
    mfot    =  r(10)

    res_NS_old  =  res_NS_now
    res_NS_now  =  res_NS

    if(.not. is_end) then
        n_res_record=  n_res_record+1
        tot_ite     =  tot_ite+1

        res_record(1 ,n_res_record) =  real(tot_ite, dpR)
        res_record(2 ,n_res_record) =  log10(res_NS)
        res_record(3 ,n_res_record) =  log10(res_RANS)
        res_record(4 ,n_res_record) =  lift
        res_record(5 ,n_res_record) =  drag
        res_record(6 ,n_res_record) =  lift_p
        res_record(7 ,n_res_record) =  lift_v
        res_record(8 ,n_res_record) =  drag_p
        res_record(9 ,n_res_record) =  drag_v
        res_record(10,n_res_record) =  mfin
        res_record(11,n_res_record) =  mfot

        if(ite_wrk .eq. 1)  res_NS_1st  =  res_NS

        if(is_st_cal_now) then
            is_converged=  res_NS .le. res_cvg_lev(lev)
            iter        =  ite_wrk
            ltmp=  is_converged .or. (mod(iter,10) .eq. 1) .or. &
                & (iter .eq. max_iter_lev(lev))
        else
            if(is_bdf_now) then
                iter=  ite_wrk
                is_converged= (res_NS/res_NS_1st .le. rhs_cvg_uns) .and. &
                            & (ite_wrk .ge. min_subiter_uns)
                ltmp= (mod(iter,10) .eq. 1) .or. (iter .eq. max_subiter_uns) .or. &
                    &  is_converged
            elseif(uns_method .eq. uns_IMEX) then
                iter=  uns_iter
                ltmp=  .true.
            else
                iter=  uns_iter
                ltmp= (mod(iter,10) .eq. 1) .or. (iter .eq. max_uns_iter)
            end if
        end if
        if(is_show_res .and. (myid .eq. 0) .and. (iter .eq. 1)) then
            print*,'----------------------------------------------------------'
            if(is_bdf_now) then
                print*,'Unsteady iter=',uns_iter
                print*,'Physical Time=',tph
            end if

            if(is_flow_int) then
                if(RANS_model .eq. RANS_WA) then
                    write(6,fmt='(A8,6A12)'),'Iter','Res','Res_WA','Lift','Drag', &
                        &  'mf_inl','mf_out'
                elseif(RANS_model .eq. RANS_SST) then
                    write(6,fmt='(A8,6A12)'),'Iter','Res','Res_SST','Lift','Drag', &
                        &  'mf_inl','mf_out'
                elseif(is_RANS) then
                    write(6,fmt='(A8,6A12)'),'Iter','Res','Res_SA','Lift','Drag', &
                        &  'mf_inl','mf_out'
                else
                    write(6,fmt='(A8,5A12)'),'Iter','Res','Lift','Drag','mf_inl','mf_out'
                end if
            else
                if(RANS_model .eq. RANS_WA) then
                    write(6,fmt='(A8,6A12)'),'Iter','Res','Res_WA','Lift','Drag', &
                        &  'Drag_p','Drag_v'
                elseif(RANS_model .eq. RANS_SST) then
                    write(6,fmt='(A8,6A12)'),'Iter','Res','Res_SST','Lift','Drag', &
                        &  'mf_inl','mf_out'
                elseif(is_RANS) then
                    write(6,fmt='(A8,6A12)'),'Iter','Res','Res_SA','Lift','Drag', &
                        &  'Drag_p','Drag_v'
                else
                    write(6,fmt='(A8,5A12)'),'Iter','Res','Lift','Drag','Drag_p','Drag_v'
                end if
            end if
        end if

        if(is_show_res .and. (myid .eq. 0) .and. ltmp) then
            if(is_flow_int) then
                if(is_RANS) then
                    write(unit=6,fmt='(I8, 6ES12.4)'),iter,res_NS,res_RANS,lift,drag,mfin,mfot
                else
                    write(unit=6,fmt='(I8, 5ES12.4)'),iter,res_NS,lift,drag,mfin,mfot
                end if
            else
                if(is_RANS) then
                    write(unit=6,fmt='(I8, 6ES12.4)'),iter,res_NS,res_RANS,lift,drag,drag_p,drag_v
                else
                    write(unit=6,fmt='(I8, 5ES12.4)'),iter,res_NS,lift,drag,drag_p,drag_v
                end if
            end if
        end if

        if((tot_ite .eq. 1) .and. (myid .eq. 0)) then
            open(unit=10,file=trim(adjustl(res_file)))
            tec =  'variables="Iteration","log<sub>10</sub>(residual)"'
            if(is_RANS) tec =  trim(adjustl(tec))//',"log<sub>10</sub>(residual_SA)"'
            tec =  trim(adjustl(tec))//',"lift","drag","lift_p","lift_v"'
            tec =  trim(adjustl(tec))//',"drag_p","drag_v","mf_inl","mf_out"'
            write(unit=10,fmt='(A)'),trim(adjustl(tec))
            write(unit=10,fmt='(A)'),'zone T="Navier-Stokes"'
            close(10)
        end if
    end if
    ltmp=((.not. is_end) .and. (n_res_record .eq. max_res_record)) .or. &
        & (is_end .and. (n_res_record .gt. 0))
    if(ltmp) then
        if(myid .eq. 0) then
            open(unit=10,file=trim(adjustl(res_file)),status='old',position='append')
            if(is_RANS) then
                do i=1,n_res_record
                    write(unit=10,fmt='(F9.1,1x,10ES20.12)'),res_record(1:11,i)
                end do
            else
                do i=1,n_res_record
                    write(unit=10,fmt='(F9.1,1x, 9ES20.12)'),res_record(1:2,i), &
                        &  res_record(4:11,i)
                end do
            end if
        end if
        n_res_record=  0
    end if
    if(is_end) then
        if(is_output_cgns) then
            if(is_output_vc .and. (.not. is_uns_cal)) &
                &  call fv_wr_cgns   (1, trim(adjustl(mesh_name))//'_slv.cgns')
            if(is_output_cc .and. (.not. is_uns_cal)) &
                &  call fv_wr_cgns_cc(2, trim(adjustl(mesh_name))//'_slv_cc.cgns')
        else
            call fv_wr_tec(1)
        end if
        if(is_output_cpcf)  call fv_wr_boundary_cgns
        call fv_get_yp(lev)

        call CPU_time(s(1))
        call mpi_gather(s(1)-mbt_t1, 1, mpi_dpR, cost, 1, mpi_dpR, 0, &
            &  mpi_comm_world, mpi_err)
        if((myid .eq. 0) .and. (tot_ite .gt. 0)) then
            open(unit=10,file=trim(adjustl(res_file)),status='old',position='append')
            do i=0,nprc-1
                write(unit=10,fmt='(A18,ES20.12,A9,I8)'),'#Elapsed time    =', &
                    &  cost(i+1),'s on prc',i
            end do
            if(is_vis_cal) then
                write(10, fmt='(A18,ES20.12)'),'#Minimum y+      =',min_yp
                write(10, fmt='(A18,ES20.12)'),'#Maximum y+      =',max_yp
            end if
            close(10)
        end if

        call fv_integrate_force_mass
    end if

    return
    end subroutine fv_slv_monitor
