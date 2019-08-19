!-------------------------------------------------------------------------------
!   output vertex-centered solution, tecplot.
!-------------------------------------------------------------------------------
    subroutine fv_wr_tec(content)
    use var_kind_def
    use var_cgns
    use var_fv
    use var_global, only: mesh_name,n_dim
    use var_mesh
    use var_parallel, only: myid
    implicit none
    integer(dpL),intent(in):: content
    character(len=1000):: str,tec
    integer(dpI):: npe,n_vtx,isec,n_ele,i,iele

    if(content .eq. 1) then
        call fv_get_vc(content)
    elseif(content .eq. 2) then
        n_vc_to_write   =  3
        vc_sol_name(1)  =  'J_CoordinateX'
        vc_sol_name(2)  =  'J_CoordinateY'
        vc_sol_name(3)  =  'J_CoordinateZ'
    else
        stop 'Error: wrong input, fv_wr_tec.'
    end if

    n_vtx   =  mesh(0)%n_vtx
    n_ele   =  0
    do isec=mesh(0)%sec_1,mesh(0)%sec_0
        if(sec(isec)%is_int)    n_ele   =  n_ele+sec(isec)%n_ele
    end do

    write(str,*),myid
    str =  trim(adjustl(mesh_name))//'_sol_'//trim(adjustl(str))//'.dat'
    open(unit=10,file=trim(adjustl(str)))

    tec =  'variables="CoordinateX","CoordinateY"'
    if(n_dim .eq. 3)    tec =  trim(adjustl(tec))//',"CoordinateZ"'
    do i=1,n_vc_to_write
        tec =  trim(adjustl(tec))//',"'//trim(adjustl(vc_sol_name(i)))//'"'
    end do
    write(unit=10,fmt='(A)'),trim(tec)

    do isec=mesh(0)%sec_1,mesh(0)%sec_0
        if(sec(isec)%is_ghost .or. sec(isec)%is_bnd)    cycle

        tec =  'zone N='
        write(str,*),n_vtx
        tec =  trim(adjustl(tec))//trim(adjustl(str))
        write(str,*),sec(isec)%n_ele
        tec =  trim(adjustl(tec))//',E='//trim(adjustl(str))
        if(sec(isec)%is_quad) then
            str =  'FEQUADRILATERAL'
        elseif(sec(isec)%is_tri) then
            str =  'FETRIANGLE'
        elseif(sec(isec)%is_hexa) then
            str =  'FEBRICK'
        else
            stop 'Error: ele_type not supported.'
        end if
        tec =  trim(adjustl(tec))//','//'zonetype='//trim(adjustl(str))
        if(isec .gt. 1) then
            write(str,*),7
            tec =  trim(adjustl(tec))//',varsharelist=([1-'//trim(adjustl(str)) &
                &  //']=1)'
        end if
        tec =  trim(adjustl(tec))//',datapacking=point'
        write(unit=10,fmt='(A)'),trim(adjustl(tec))
        if(isec .eq. 1) then
            if(content .eq. 1) then
                do i=1,n_vtx
                    write(unit=10,fmt='(10ES20.12)'),mesh(0)%xyz(1:n_dim,i), &
                        &  fv_vertex_variables(1:n_vc_to_write,i)
                end do
            else
                do i=1,n_vtx
                    write(unit=10,fmt='(10ES20.12)'),mesh(0)%xyz(1:n_dim,i), &
                        &  mesh(0)%J_xyz(1:3,i)
                end do
            end if
        end if
        npe =  sec(isec)%npe
        do iele=1,sec(isec)%n_ele
            write(unit=10,fmt='(8I8)'),sec(isec)%n2e(1:npe,iele)
        end do
    end do

    return
    end subroutine fv_wr_tec
!-------------------------------------------------------------------------------
!   output something for debugging.
!-------------------------------------------------------------------------------
    subroutine fv_debug(lev)
    use var_kind_def
    use var_fv
    use var_mesh
    use var_parallel
    implicit none
    integer(dpI),intent(in):: lev
    character(len=100):: str
    integer(dpI):: LDA,i,isec,iele,ip
    real   (dpR),allocatable:: send(:,:),recv(:,:)

    LDA =  1+1
    i   =  0
    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(sec(isec)%is_int)    i   =  i+sec(isec)%n_ele
    end do
    allocate(send(LDA,i))

    call mpi_gather(i, 1, mpi_dpI, prc_info, 1, mpi_dpI, 0, mpi_comm_world, mpi_err)
    if(myid .eq. 0) then
        i   =  prc_info(0,1)
        do ip=1,nprc-1
            i   =  i+prc_info(ip,1)
        end do
        allocate(recv(LDA,i))
    else
        allocate(recv(LDA,1))
    end if

    i   =  0
    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(.not. sec(isec)%is_int)  cycle
        do iele=1,sec(isec)%n_ele
            i   =  i+1
            send(1  ,i) =  real(sec(isec)%ID_ele_i(iele), dpR)
!           send(2:6,i) =  fv(isec)%u  (1:5,iele)
!           send(2:9,i) =  fv(isec)%gra(1:8,iele)
!           send(2:6,i) =  fv(isec)%rhs(1:5,iele)
!           send(2:6,i) =  fv(isec)%adj(1:5,iele)
!           send(2:6,i) =  fv(isec)%adj_pJu(1:5,iele)
!           send(2:6,i) =  fv(isec)%adj_rhs(1:5,iele)
!           send(2  ,i) =  fv(isec)%LHS_s(iele)
!           send(2  ,i) =  fv(isec)%turb(iele)
            send(2  ,i) =  fv(isec)%rhs (1,iele)
!           send(3  ,i) =  fv(isec)%mut (iele)
!           send(4  ,i) =  fv(isec)%rhs_turb(iele)
        end do
    end do
    i   =  i*LDA
    call r_gather(mpi_comm_world, nprc, i, send, recv)
    i   =  i/LDA

    if(myid .eq. 0) then
        call dqsortcols(.true., 1, i, 1, LDA, recv)
        write(str,*),nprc
        str =  './data/debug_'//trim(adjustl(str))//'.dat'
        open(unit=10,file=trim(adjustl(str)))

        if(.false.) then
            write(unit=10,fmt=*),LDA
            do iele=1,i
                write(unit=10,fmt='(I5,10ES20.12)'),nint(recv(1,iele)),recv(2:LDA,iele)
            end do
        else
            do iele=1,i
                write(unit=10,fmt='(10ES20.12)'),recv(2:LDA,iele)
            end do
        end if
        close(10)
    end if
    call mpi_barrier(mpi_comm_world, mpi_err)
    call mpi_finalize(mpi_err)
    stop

    return
    end subroutine fv_debug
!-------------------------------------------------------------------------------
!   interpolate cell-center variables to vertex.
!-------------------------------------------------------------------------------
    subroutine fv_get_vc(cnt)
    use var_kind_def
    use var_cgns, only: BCWallViscous
    use var_fv, v=>fv_vertex_variables,LDA=>n_vc_to_write
    use var_global
    use var_mesh
    use var_turb, transition=>transition_model
    implicit none
    integer(dpI),intent(in):: cnt
    logical(dpL):: ltmp
    integer(dpI):: isec,iele,i,j,ivtx,ID,iA(101)
    real   (dpR):: u(12),ome(3),vg(3),c(3)

    if(cnt .eq. 1) then
        vc_sol_name(1)  = 'Density'
        vc_sol_name(2)  = 'VelocityX'
        vc_sol_name(3)  = 'VelocityY'
        if(is_output_rel_uvw) then
            vc_sol_name(2)  = 'RotatingVelocityX'
            vc_sol_name(3)  = 'RotatingVelocityY'
        end if
        if(is_2d_cal) then
            LDA             =  4
            vc_sol_name(4)  = 'Pressure'
        else
            LDA             =  5
            vc_sol_name(4)  = 'VelocityZ'
            vc_sol_name(5)  = 'Pressure'
            if(is_output_rel_uvw)   vc_sol_name(4)  = 'RotatingVelocityZ'
        end if
        if(RANS_model .eq. RANS_SA) then
            LDA             =  LDA+1
            vc_sol_name(LDA)= 'TurbulentSANuTilde'
        elseif(RANS_model .eq. RANS_WA) then
            LDA             =  LDA+1
            vc_sol_name(LDA)= 'TurbulentSANuTilde'
        elseif(RANS_model .eq. RANS_SAM) then
            LDA             =  LDA+1
            vc_sol_name(LDA)= 'TurbulentSANuTilde'
        elseif(RANS_model .eq. RANS_SAC) then
            LDA             =  LDA+1
            vc_sol_name(LDA)= 'TurbulentSANuTilde'
        elseif(RANS_model .eq. RANS_SST) then
            LDA             =  LDA+1
            vc_sol_name(LDA)= 'TurbulentEnergyKinetic'
            LDA             =  LDA+1
            vc_sol_name(LDA)= 'TurbulentDissipationRate'
        end if

        if((transition .ge. transition_Menter) .and. (transition .le. transition_BC)) then
            LDA             =  LDA+1
            vc_sol_name(LDA)= 'Intermittency'
        end if
        if((transition .ge. transition_PTM) .and. (transition .le. transition_BC)) then
            call fv_exchange_intermittency(0)
        end if
    elseif(cnt .eq. 2) then
        LDA             =  12
        vc_sol_name(1 ) =  'Density'
        vc_sol_name(2 ) =  'VelocityX'
        vc_sol_name(3 ) =  'VelocityY'
        vc_sol_name(4 ) =  'VelocityZ'
        vc_sol_name(5 ) =  'Pressure'
        vc_sol_name(6 ) =  'ReynoldsStressXX'
        vc_sol_name(7 ) =  'ReynoldsStressXY'
        vc_sol_name(8 ) =  'ReynoldsStressXZ'
        vc_sol_name(9 ) =  'ReynoldsStressYY'
        vc_sol_name(10) =  'ReynoldsStressYZ'
        vc_sol_name(11) =  'ReynoldsStressZZ'
        vc_sol_name(12) =  'PressureRMS'
    else
        stop 'Error: wrong input, fv_get_vc.'
    end if
    if(.not. allocated(v)) then
        allocate(v(LDA, mesh(0)%n_vtx), stat=err_mem)
        if(err_mem .ne. 0)  stop 'Error: fv_get_vc fails to allocate memory.'
    else
        if(size(v) .lt. LDA*mesh(0)%n_vtx) then
            deallocate(v)
            allocate(v(LDA, mesh(0)%n_vtx), stat=err_mem)
            if(err_mem .ne. 0)  stop 'Error: fv_get_vc fails to allocate memory.'
        end if
    end if

!   ----------------------------------------------------------------------------
!   cell centered value to vertex centered value.
    iA(1)   =  1
    do isec=mesh(0)%sec_1,mesh(0)%sec_0
        iA(isec+1)  =  iA(isec)+sec(isec)%n_ele
    end do

    v   =  0.0d0
    do ivtx=1,mesh(0)%n_vtx
    do j=mesh(0)%iA_e2n(ivtx),mesh(0)%iA_e2n(ivtx+1)-1
        i   =  mesh(0)%jA_e2n(j)
        ltmp=  .true.
        do isec=mesh(0)%sec_1,mesh(0)%sec_0
            if((i .ge. iA(isec)) .and. (i .lt. iA(isec+1))) then
                ltmp=  .false.
                iele=  i-iA(isec)+1
                exit
            end if
        end do
        if(ltmp)    stop 'Error: fv_get_vc fails to get isec.'

        if(cnt .eq. 1) then
            u(1:3)  =  fv(isec)%u(1:3,iele)
            if(is_2d_cal) then
                u(4  )  =  fv(isec)%u(  5,iele)
                if(RANS_model .eq. RANS_SA) then
                    u(5)=  fv(isec)%turb(1,iele)
                elseif(RANS_model .eq. RANS_WA) then
                    u(5)=  fv(isec)%turb(1,iele)
                elseif(RANS_model .eq. RANS_SAM) then
                    u(5)=  fv(isec)%turb(1,iele)
                elseif(RANS_model .eq. RANS_SAC) then
                    u(5)=  fv(isec)%turb(1,iele)
                elseif(RANS_model .eq. RANS_SST) then
                    u(5)=  fv(isec)%turb(1,iele)
                    u(6)=  fv(isec)%turb(2,iele)
                end if

                if(transition .eq. transition_Menter) then
                    u(6)=  fv(isec)%turb(2,iele)
                elseif(transition .eq. transition_Coder) then
                    u(6)=  fv(isec)%turb(2,iele)
                elseif(transition .eq. transition_PTM) then
                    u(6)=  fv(isec)%intermittency(iele)
                elseif(transition .eq. transition_AGS) then
                    u(6)=  fv(isec)%intermittency(iele)
                elseif(transition .eq. transition_BC) then
                    u(6)=  fv(isec)%intermittency(iele)
                end if
            else
                u(4:5)  =  fv(isec)%u(4:5,iele)
                if(RANS_model .eq. RANS_SA) then
                    u(6)=  fv(isec)%turb(1,iele)
                elseif(RANS_model .eq. RANS_WA) then
                    u(6)=  fv(isec)%turb(1,iele)
                elseif(RANS_model .eq. RANS_SAM) then
                    u(6)=  fv(isec)%turb(1,iele)
                elseif(RANS_model .eq. RANS_SAC) then
                    u(6)=  fv(isec)%turb(1,iele)
                elseif(RANS_model .eq. RANS_SST) then
                    u(6)=  fv(isec)%turb(1,iele)
                    u(7)=  fv(isec)%turb(2,iele)
                end if

                if(transition .eq. transition_Menter) then
                    u(7)=  fv(isec)%turb(2,iele)
                elseif(transition .eq. transition_Coder) then
                    u(7)=  fv(isec)%turb(2,iele)
                elseif(transition .eq. transition_PTM) then
                    u(7)=  fv(isec)%intermittency(iele)
                elseif(transition .eq. transition_AGS) then
                    u(7)=  fv(isec)%intermittency(iele)
                elseif(transition .eq. transition_BC) then
                    u(7)=  fv(isec)%intermittency(iele)
                end if
            end if
        else
            u(1:5 ) =  fv(isec)%uns_ra(1:5,iele)
            u(6:12) =  fv(isec)%stress(1:7,iele)
        end if
        v(1:LDA,ivtx)   =  v(1:LDA,ivtx)+u(1:LDA)*mesh(0)%w_e2n(j)
    end do
    end do
!   cell centered value to vertex centered value.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   enforce the boundary condition.
    do isec=mesh(0)%sec_1,mesh(0)%sec_0
        if(.not. sec(isec)%is_bnd)  cycle
        if(sec(isec)%bct .ne. BCWallViscous)    cycle
        ID  =  sec(isec)%ID_sec_g

        do iele=1,sec(isec)%n_ele
        do i=1,sec(isec)%npe
            ivtx=  sec(isec)%n2e(i,iele)
            if(ivtx .gt. mesh(0)%n_vtx) cycle

            if(cnt .eq. 1) then
                if(sec_motion_type(ID) .eq. 0) then
                    vg  =  0.0d0
                elseif(sec_motion_type(ID) .eq. 1) then
                    vg  =  sec_motion_speed(ID)*translation_axis
                else
                    ome     =  sec_motion_speed(ID)*rotation_axis
                    c(1:3)  =  sec(isec)%cen(1:3,iele)
                    vg(1)   =  ome(2)*c(3)-ome(3)*c(2)
                    vg(2)   =  ome(3)*c(1)-ome(1)*c(3)
                    vg(3)   =  ome(1)*c(2)-ome(2)*c(1)
                end if

                if(is_2d_cal) then
                    v(2:3,ivtx) =  vg(1:2)
                else
                    v(2:4,ivtx) =  vg(1:3)
                end if
                if((RANS_model .eq. 1) .or. (RANS_model .eq. 2)) then
                    v(LDA,ivtx) =  0.0d0
                elseif(RANS_model .eq. RANS_SAM) then
                    v(LDA-1,ivtx)   =  0.0d0
                elseif(RANS_model .eq. RANS_SAC) then
                    v(LDA-1,ivtx)   =  0.0d0
                elseif(RANS_model .eq. RANS_SST) then
                    v(LDA-1,ivtx)   =  0.0d0
                end if
            else
                v(6:11,ivtx)=  0.0d0
            end if
        end do
        end do
    end do
!   enforce the boundary condition.
!   ----------------------------------------------------------------------------

    if((.not. is_output_rel_uvw) .or. (cnt .ne. 1)) return
    do isec=mesh(0)%sec_1,mesh(0)%sec_0
        if(.not. sec(isec)%is_int)  cycle
        ID  =  sec(isec)%ID_sec_g

        do iele=1,sec(isec)%n_ele
        do i=1,sec(isec)%npe
            ivtx=  sec(isec)%n2e(i,iele)
            if(ivtx .gt. mesh(0)%n_vtx) cycle

            if(sec_motion_type(ID) .eq. 0) then
                vg  =  0.0d0
            elseif(sec_motion_type(ID) .eq. 1) then
                vg  =  sec_motion_speed(ID)*translation_axis
            else
                ome     =  sec_motion_speed(ID)*rotation_axis
                c(1:3)  =  mesh(0)%xyz(1:3,ivtx)
                vg(1)   =  ome(2)*c(3)-ome(3)*c(2)
                vg(2)   =  ome(3)*c(1)-ome(1)*c(3)
                vg(3)   =  ome(1)*c(2)-ome(2)*c(1)
            end if

            if(is_2d_cal) then
                v(2:3,ivtx) =  v(2:3,ivtx)-vg(1:2)
            else
                v(2:4,ivtx) =  v(2:4,ivtx)-vg(1:3)
            end if
        end do
        end do
    end do

    return
    end subroutine fv_get_vc
!-------------------------------------------------------------------------------
!   output vertex variables for FV.
!-------------------------------------------------------------------------------
    subroutine fv_wr_cgns(content,file_name)
    use var_kind_def
    use var_cgns
    use var_fv
    use var_global
    use var_load_balance
    use var_mesh
    use var_parallel
    implicit none
    character(len=*),intent(in):: file_name
    integer(dpI),intent(in):: content
    logical(dpL):: ltmp
    integer(dpI):: nvtx,nele,isize(3),ivtx,ID,v1,v0,iseg,i,j,m,n,LDA,ip,e1,e0, &
                &  isec,iele,nnz,itype,ib,iboco
    real   (dpR):: v(100)

!   ----------------------------------------------------------------------------
!   every processor prepares the vertex based data to be output.
    call fv_get_vc(content)
!   every processor prepares the vertex based data to be output.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   how many vertex/elements will be output in one call of partial_write.
    nvtx=  min(2**20, n_vtx_g)
    nele=  min(2**20, n_ele_g)
!   how many vertex/elements will be output in one call of partial_write.
!   ----------------------------------------------------------------------------

    n   =  n_vc_to_write
    LDA =  1+n_dim+n
    if(.not. allocated(r_send)) then
        allocate(r_send    (LDA * nvtx  ), stat=err_mem)
    else
        if(size(r_send) .lt. LDA*nvtx) then
            deallocate(r_send)
            allocate(r_send(LDA*nvtx), stat=err_mem)
        end if
    end if

    if(.not. allocated(r_recv)) then
        allocate(r_recv    (LDA * nvtx  ), stat=err_mem)
    else
        if(size(r_recv) .lt. LDA*nvtx) then
            deallocate(r_recv)
            allocate(r_recv(LDA*nvtx), stat=err_mem)
        end if
    end if

    if(.not. allocated(r_recv_all)) then
        allocate(r_recv_all(nvtx*(LDA-1)), stat=err_mem)
    else
        if(size(r_recv_all) .lt. (LDA-1)*nvtx) then
            deallocate(r_recv_all)
            allocate(r_recv_all((LDA-1)*nvtx), stat=err_mem)
        end if
    end if

    if(.not. allocated(ele_send)) then
        if(is_2d_cal) then
            ID  =  5
        else
            ID  =  9
        end if
        allocate(ele_send(ID, max(nele, mesh(0)%n_vtx)), stat=err_mem)
        if(myid .eq. 0) allocate(ele_recv(ID, nele), stat=err_mem)
    end if

    if(myid .eq. 0) then
        call cg_open_f(trim(adjustl(file_name)), MODE_WRITE, ifile, cg_err)
        if(is_2d_cal) then
            call cg_base_write_f(ifile,'Base',2,2,ibase,cg_err)
        else
            call cg_base_write_f(ifile,'Base',3,3,ibase,cg_err)
        end if
        if(is_output_bnd) then
            isize   = (/n_vtx_g, n_ele_g , 0/)
        else
            isize   = (/n_vtx_g, n_elei_g, 0/)
        end if
        call cg_zone_write_f(ifile, ibase, 'zone_1', isize, unstructured, &
            &  izone, cg_err)
        call cg_sol_write_f(ifile, ibase, izone, 'FlowSolution', Vertex, &
            &  iflow, cg_err)
    end if

!   ----------------------------------------------------------------------------
!   output vertex based information.
    do iseg=1,(n_vtx_g-1)/nvtx+1
        v1  =  1+nvtx*(iseg-1)
        v0  =  min(nvtx*iseg, n_vtx_g)

!       every processor checks the number of vertex to be output in this partial write.
        ele_send=  0
        do isec=mesh(0)%sec_1,mesh(0)%sec_0
            if(.not. sec(isec)%is_int)  cycle
            do iele=1,sec(isec)%n_ele
            do i=1,sec(isec)%npe
                ele_send(1,sec(isec)%n2e(i,iele))   =  1
            end do
            end do
        end do

        m   =  0
        do ivtx=1,mesh(0)%n_vtx
            ID  =  mesh(0)%id_vtx(ivtx)
            if((ID .lt. v1) .or. (ID .gt. v0))  cycle
            if(ele_send(1,ivtx) .ne. 1) cycle
            m   =  m+1
            v(1          )  =  real(ID, dpR)
            v(2:1+n_dim  )  =  mesh(0)%xyz(1:n_dim,ivtx)
            v(2+n_dim:LDA)  =  fv_vertex_variables(1:n,ivtx)
            r_send(1+LDA*(m-1):LDA*m)   =  v(1:LDA)
        end do

        call mpi_gather(m, 1, mpi_dpI, prc_info, 1, mpi_dpI, 0, mpi_comm_world, mpi_err)

        if((myid .ne. 0) .and. (m .gt. 0))  &
            &  call mpi_send(r_send, LDA*m, mpi_dpR, 0, myid, mpi_comm_world, mpi_err)

        if(myid .ne. 0) cycle
        do ip=0,nprc-1
            if(prc_info(ip,1) .le. 0)   cycle
            if(ip .eq. 0) then
                call DCOPY(LDA*prc_info(ip,1), r_send, 1, r_recv, 1)
            else
                call mpi_recv(r_recv, LDA*prc_info(ip,1), mpi_dpR, ip, ip, &
                    &  mpi_comm_world, mpi_status, mpi_err)
            end if
!           decode the received information.
            do i=1,prc_info(ip,1)
                ID  =  nint(r_recv(1+LDA*(i-1)), dpI)
                ivtx=  ID-v1+1
                forall(j=1:LDA-1)   r_recv_all(ivtx+nvtx*(j-1)) =  r_recv(j+1+LDA*(i-1))
            end do
        end do

        call cg_coord_partial_write_f(ifile, ibase, izone, RealDouble,'CoordinateX', &
            &  v1, v0, r_recv_all(1       ), icoord, cg_err)
        call cg_coord_partial_write_f(ifile, ibase, izone, RealDouble,'CoordinateY', &
            &  v1, v0, r_recv_all(1+nvtx  ), icoord, cg_err)
        if(.not. is_2d_cal) then
        call cg_coord_partial_write_f(ifile, ibase, izone, RealDouble,'CoordinateZ', &
            &  v1, v0, r_recv_all(1+nvtx*2), icoord, cg_err)
        end if
        do i=n_dim+1,LDA-1
            m   =  i-n_dim
            call cg_field_partial_write_f(ifile, ibase, izone, iflow, RealDouble, &
                &  vc_sol_name(m), v1, v0, r_recv_all(1+nvtx*(i-1)), ifield, cg_err)
            if(cg_err .ne. 0)   stop 'Error: fails to output vertex variable.'
        end do
    end do
!   output vertex based information.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   output element based information.
    if(is_2d_cal) then
        LDA =  5
    else
        LDA =  9
    end if
    do iseg=1,(n_elei_g-1)/nele+1
        e1  =  1+nele*(iseg-1)
        e0  =  min(iseg*nele, n_elei_g)

!       every processor checks the number of elements to be output in this partial write.
        m   =  0
        do isec=mesh(0)%sec_1,mesh(0)%sec_0
            if(.not. sec(isec)%is_int)  cycle
            do iele=1,sec(isec)%n_ele
                ID  =  sec(isec)%ID_ele_i(iele)
                if((ID .lt. e1) .or. (ID .gt. e0))  cycle
                m   =  m+1
                ele_send(1,m)   =  ID
                do i=1,sec(isec)%npe
                    ele_send(i+1,m) =  mesh(0)%id_vtx(sec(isec)%n2e(i,iele))
                end do
            end do
        end do

        call mpi_gather(m, 1, mpi_dpI, prc_info, 1, mpi_dpI, 0, mpi_comm_world, mpi_err)

        if((myid .ne. 0) .and. (m .gt. 0))  &
            &  call mpi_send(ele_send, LDA*m, mpi_dpI, 0, myid, mpi_comm_world, mpi_err)

        if(myid .ne. 0) cycle
        m   =  1
        do ip=0,nprc-1
            if(prc_info(ip,1) .le. 0)   cycle
            if(ip .eq. 0) then
                call ICOPY(LDA*prc_info(ip,1), ele_send, 1, ele_recv(1,m), 1)
            else
                call mpi_recv(ele_recv(1,m), LDA*prc_info(ip,1), mpi_dpI, ip, ip, &
                    &  mpi_comm_world, mpi_status, mpi_err)
            end if
            m   =  m+prc_info(ip,1)
        end do
        call iqsortcols(.true., 1, m-1, 1, LDA, ele_recv)

        do isec=1,n_seci_g
            if((e1 .gt. seci_g(4,isec)) .or. (e0 .lt. seci_g(3,isec)))  cycle
            i   =  max(e1, seci_g(3,isec))
            j   =  min(e0, seci_g(4,isec))
            call prepare_n2e(LDA, seci_g(2,isec), i-e1+1, j-e1+1, ele_recv)
            call cg_section_partial_write_f(ifile, ibase, izone, sec_name(isec), &
                &  seci_g(1,isec), i, j, 0, ele_recv, isection, cg_err)
            if(cg_err .ne. 0)   stop 'Error: fails to output n2e.'
        end do
    end do
!   output element based information.
!   ----------------------------------------------------------------------------

    if(is_output_bnd .and. (myid .eq. 0)) then
        ID  =  0
        nnz =  1
        e1  =  n_elei_g+1
        do isec=1,n_sec_g
            itype   =  sec_g(1,isec)
            if(is_2d_cal) then
                ltmp= (itype .eq. BAR_2) .or. (itype .eq. BAR_3)
            else
                ltmp= (itype .eq. TRI_3) .or. (itype .eq. QUAD_4) .or. &
                    & (itype .eq. QUAD_8) .or. (itype .eq. QUAD_9)
            end if
            if(.not. ltmp)  cycle
            ID  =  ID+1
            e0  =  e1+sec_g(4,isec)-sec_g(3,isec)

            ltmp=  .false.
            do ib=1,n_bocos
                ltmp= (bocoinfo(2,ib) .eq. sec_g(3,isec)) .and. &
                    & (bocoinfo(3,ib) .eq. sec_g(4,isec))
                if(ltmp)    exit
            end do
            if(.not. ltmp)  stop 'Error: fails to map boundary to section.'

            call cg_section_write_f(ifile, ibase, izone, bsec_name(ID), &
                &  sec_g(1,isec), e1, e0, 0, jA_n2e_b(nnz), isection, cg_err)
            if(cg_err .ne. 0)   stop 'Error: fails to output boundary n2e.'

            if(.false.) then
            call cg_boco_write_f(ifile, ibase, izone, boconame(ib), bocoinfo(1,ib), &
                &  ElementRange, 2, (/e1, e0/), iboco, cg_err)
            if(cg_err .ne. 0)   stop 'Error: fails to output boundary info.'
            call cg_goto_f(ifile,ibase,cg_err,'Zone_t',1,'ZoneBC_t',1,'BC_t',iboco,'end')
            call cg_gridlocation_write_f(FaceCenter, cg_err)
            if(cg_err .ne. 0)   stop 'Error: fails to output boundary location.'
            end if

            nnz =  nnz+sec_g(2,isec)*(sec_g(4,isec)-sec_g(3,isec)+1)
            e1  =  e0+1
        end do
    end if

    if(myid .eq. 0) call cg_close_f(ifile, cg_err)

    if(allocated(ele_send)) deallocate(ele_send)
    if(allocated(ele_recv)) deallocate(ele_recv)

    return
    contains
!       ------------------------------------------------------------------------
!       prepare n2e.
!       ------------------------------------------------------------------------
        subroutine prepare_n2e(LDA,npe,L,R,b)
        implicit none
        integer(dpI),intent(in):: LDA,npe,L,R
        integer(dpI):: b(*),i,j,v(10)

        j   =  1
        do i=L,R
            call ICOPY(npe, b(2+LDA*(i-1)), 1, v   , 1)
            call ICOPY(npe, v             , 1, b(j), 1)
            j   =  j+npe
        end do

        return
        end subroutine prepare_n2e
    end subroutine fv_wr_cgns
!-------------------------------------------------------------------------------
!   output cell-centered solution of the boundary.
!-------------------------------------------------------------------------------
    subroutine fv_wr_boundary
    use var_kind_def
    use var_bndv
    use var_cgns
    use var_fv
    use var_global, only: is_2d_cal,mesh_name,err_mem
    use var_mesh
    use var_parallel, only: myid
    use var_slv, only: is_vis_cal
    use var_turb, only: transition_model,transition_AGS
    implicit none
    character(len=1000):: sol_name(9),str,tec
    integer(dpI):: n_vtx,n_sec,isec,iele,i,N,v(4)
    integer(dpI),allocatable:: ID_vtx(:)

    allocate(ID_vtx(mesh(0)%n_vtx), stat=err_mem)

    sol_name(1) =  '"CoordinateX"'
    sol_name(2) =  '"CoordinateY"'
    sol_name(3) =  '"CoordinateZ"'
    sol_name(4) =  '"Cp"'
    sol_name(5) =  '"Cf"'
    N           =  5
    if(transition_model .eq. transition_AGS) then
        sol_name(6) =  '"Re_theta"'
        N           =  6
    end if

    if(.not. is_vis_cal)    N   =  N-1

    write(str,*),myid
    str =  trim(adjustl(mesh_name))//'_bnd_'//trim(adjustl(str))//'.dat'
    open(unit=10,file=trim(adjustl(str)))

    str =  'variables='
    do i=1,N
        str =  trim(adjustl(str))//trim(adjustl(sol_name(i)))
    end do
    write(unit=10,fmt='(A)'),trim(adjustl(str))
    n_sec   =  0
    do isec=mesh(0)%sec_1,mesh(0)%sec_0
        if(.not. sec(isec)%is_bnd)  cycle

        n_sec   =  n_sec+1

        ID_vtx  =  0
        do iele=1,sec(isec)%n_ele
        do i=1,sec(isec)%npe
            ID_vtx(sec(isec)%n2e(i,iele))   =  1
        end do
        end do

        n_vtx   =  0
        do i=1,mesh(0)%n_vtx
            if(ID_vtx(i) .le. 0)    cycle
            n_vtx       =  n_vtx+1
            ID_vtx(i)   =  n_vtx
        end do

        tec =  'zone N='
        write(str,*),n_vtx
        tec =  trim(adjustl(tec))//trim(adjustl(str))
        write(str,*),sec(isec)%n_ele
        tec =  trim(adjustl(tec))//',E='//trim(adjustl(str))
        if(sec(isec)%is_quad) then
            str =  'FEQUADRILATERAL'
        elseif(sec(isec)%is_bar) then
            str =  'FELINESEG'
        else
            str =  'FETRIANGLE'
        end if
        tec =  trim(adjustl(tec))//','//'zonetype='//trim(adjustl(str))
        tec =  trim(adjustl(tec))//',datapacking=block'
        write(str,*),N
        if(N .gt. 4) then
            tec =  trim(adjustl(tec))//',varlocation=([4-'//trim(adjustl(str))//']=cellcentered)'
        else
            tec =  trim(adjustl(tec))//',varlocation=([4]=cellcentered)'
        end if
        write(unit=10,fmt='(A)'),trim(adjustl(tec))

!       output Coodinates.
        do i=1,mesh(0)%n_vtx
            if(ID_vtx(i) .gt. 0)    write(unit=10,fmt=*),mesh(0)%xyz(1,i)
        end do
        do i=1,mesh(0)%n_vtx
            if(ID_vtx(i) .gt. 0)    write(unit=10,fmt=*),mesh(0)%xyz(2,i)
        end do
        do i=1,mesh(0)%n_vtx
            if(ID_vtx(i) .le. 0)    cycle
            if(is_2d_cal) then
                write(unit=10,fmt=*),0.0d0
            else
                write(unit=10,fmt=*),mesh(0)%xyz(3,i)
            end if
        end do

!       output cp.
        do iele=1,sec(isec)%n_ele
            write(unit=10,fmt=*),(fv(isec)%u(5,iele)-u_fs(5))/(0.5d0*u_fs(1)*ua_fs**2)
        end do

!       output cf.
        if(is_vis_cal) then
            if(sec(isec)%bct .eq. BCWallViscous) then
                write(unit=10,fmt=*),fv(isec)%bnd_solution(2,:)
            else
                write(unit=10,fmt=*),(0.0d0, iele=1,sec(isec)%n_ele)
            end if
        end if

!       output Re_theta.
        if(transition_model .eq. transition_AGS) then
            if(sec(isec)%bct .eq. BCWallViscous) then
                write(unit=10,fmt=*),fv(isec)%Re_theta
            else
                write(unit=10,fmt=*),(0.0d0, iele=1,sec(isec)%n_ele)
            end if
        end if

        do iele=1,sec(isec)%n_ele
            do i=1,sec(isec)%npe
                v(i)=  ID_vtx(sec(isec)%n2e(i,iele))
            end do
            if(sec(isec)%is_bar) then
                write(unit=10,fmt='(2I8)'),v(1:2)
            else
                write(unit=10,fmt='(4I8)'),v(1:sec(isec)%npe)
            end if
        end do
    end do
    close(10)

    if(allocated(ID_vtx))   deallocate(ID_vtx)

    return
    end subroutine fv_wr_boundary
!-------------------------------------------------------------------------------
!   output local element based solution.
!-------------------------------------------------------------------------------
    subroutine fv_wr_cc_L(solution)
    use var_kind_def
    use var_cgns
    use var_fv
    use var_global, only: is_2d_cal,n_dim,err_mem
    use var_mesh, only: sec,mesh
    use var_parallel
    implicit none
    integer(dpI),intent(in):: solution
    character(len=80):: str,sec_name(2000),solution_name
    integer(dpI):: isec,nvtx,nele,ele1,ele0,i,isize(3),ele
    real   (dpR),allocatable:: xyz(:,:),rbuf(:)

    write(str,*),myid
    if(solution .eq. 1) then
        str =  './data/dnw_local_'//trim(adjustl(str))//'.cgns'
        solution_name   = 'TurbulentDistance'
    elseif(solution .eq. 2) then
        str =  './data/gls_local_'//trim(adjustl(str))//'.cgns'
        solution_name   = 'DDESlength'
    else
        str =  './data/Ducros_local_'//trim(adjustl(str))//'.cgns'
        solution_name   = 'DucrosSensor'
    end if

    nvtx=  mesh(0)%n_vtx
    nele=  0
    do isec=mesh(0)%sec_1,mesh(0)%sec_0
        if(sec(isec)%is_int)    nele=  nele+sec(isec)%n_ele
    end do
    allocate(xyz (nvtx,n_dim), stat=err_mem)
    allocate(rbuf(nele      ), stat=err_mem)

    call DCOPY(mesh(0)%n_vtx, mesh(0)%xyz(1,1), n_dim, xyz(1,1), 1)
    call DCOPY(mesh(0)%n_vtx, mesh(0)%xyz(2,1), n_dim, xyz(1,2), 1)
    if(.not. is_2d_cal) then
    call DCOPY(mesh(0)%n_vtx, mesh(0)%xyz(3,1), n_dim, xyz(1,3), 1)
    end if

    i   =  1
    do isec=mesh(0)%sec_1,mesh(0)%sec_0
        if(.not. sec(isec)%is_int)  cycle

        if(solution .eq. 1) then
            call DCOPY(sec(isec)%n_ele, sec(isec)%dnw, 1, rbuf(i), 1)
        elseif(solution .eq. 2) then
            call DCOPY(sec(isec)%n_ele, sec(isec)%gls, 1, rbuf(i), 1)
        else
            call DCOPY(sec(isec)%n_ele, fv(isec)%shock_sensor, 1, rbuf(i), 1)
        end if
        i   =  i+sec(isec)%n_ele
    end do

    call cg_open_f(trim(adjustl(str)), MODE_WRITE, ifile, cg_err)
    if(is_2d_cal) then
        call cg_base_write_f(ifile,'Base',2,2,ibase,cg_err)
    else
        call cg_base_write_f(ifile,'Base',3,3,ibase,cg_err)
    end if
    isize   = (/mesh(0)%n_vtx, nele, 0/)
    call cg_zone_write_f(ifile, ibase, 'zone_1', isize, unstructured, izone, cg_err)
    call cg_coord_write_f(ifile, ibase, izone, RealDouble, 'CoordinateX', &
        &  xyz(1,1), icoord, cg_err)
    call cg_coord_write_f(ifile, ibase, izone, RealDouble, 'CoordinateY', &
        &  xyz(1,2), icoord, cg_err)
    if(.not. is_2d_cal) then
    call cg_coord_write_f(ifile, ibase, izone, RealDouble, 'CoordinateZ', &
        &  xyz(1,3), icoord, cg_err)
    end if

    str =  'section'
    call num_to_str(mesh(0)%sec_0, str, sec_name)

    ele1=  1
    ele0=  0
    do isec=mesh(0)%sec_1,mesh(0)%sec_0
        if(.not. sec(isec)%is_int)  cycle
        ele0=  ele1+sec(isec)%n_ele-1
        ele =  sec(isec)%ele_type
        call cg_section_write_f(ifile, ibase, izone, sec_name(isec), &
            &  ele_type_1st(ele), ele1, ele0, 0, sec(isec)%n2e, i, cg_err)
        ele1=  ele0+1
    end do

    call cg_sol_write_f(ifile, ibase, izone, 'FlowSolution', CellCenter, iflow, cg_err)
    if(cg_err .ne. 0)   stop 'Error: cg_sol_write_f fails.'
    call cg_field_write_f(ifile, ibase, izone, iflow, RealDouble, &
        &  trim(solution_name), rbuf, ifield, cg_err)
    if(cg_err .ne. 0)   stop 'Error: cg_field_write_f fails.'

    call cg_close_f(ifile, cg_err)

    if(allocated(xyz )) deallocate(xyz)
    if(allocated(rbuf)) deallocate(rbuf)

    return
    end subroutine fv_wr_cc_L
!-------------------------------------------------------------------------------
!   output cell-centered solution, tecplot, local.
!-------------------------------------------------------------------------------
    subroutine fv_wr_tec_cc(solution)
    use var_kind_def
    use var_cgns
    use var_fv
    use var_global, only: is_2d_cal,mesh_name,n_dim
    use var_mesh
    use var_parallel, only: myid
    implicit none
    integer(dpI),intent(in):: solution
    character(len=1000):: str,tec
    logical(dpL):: ltmp
    integer(dpI):: n_vtx,isec,i,iele,v(8),L,R
    integer(dpI),allocatable:: ibuf(:)

    write(str,*),myid
    str =  trim(adjustl(mesh_name))//'_sol_'//trim(adjustl(str))//'.dat'
    open(unit=10,file=trim(adjustl(str)))

    tec =  'variables="CoordinateX","CoordinateY"'
    if(n_dim .eq. 3)    tec =  trim(adjustl(tec))//',"CoordinateZ"'
    if(solution .eq. 1) then
        tec =  trim(adjustl(tec))//',"Density","VelocityX","VelocityY","VelocityZ","Pressure"'
    elseif(solution .eq. 2) then
        tec =  trim(adjustl(tec))//',"shock_sensor"'
    elseif(solution .eq. 3) then
        tec =  trim(adjustl(tec))//',"k_modeled","k_resolved","DDES_fd"'
    elseif(solution .eq. 4) then
        tec =  trim(adjustl(tec))//',"v1","v2","v3","v4","v5"'
    elseif(solution .eq. 5) then
        tec =  trim(adjustl(tec))//',"intermittency"'
    else
        stop 'Error: output cell-centered data not supported.'
    end if
    write(unit=10,fmt='(A)'),trim(tec)

    allocate(ibuf(mesh(0)%n_vtx))
    do isec=mesh(0)%sec_1,mesh(0)%sec_0
        if(.not. sec(isec)%is_int)  cycle

        ibuf=  0
        do iele=1,sec(isec)%n_ele
        do i=1,sec(isec)%npe
            ibuf(sec(isec)%n2e(i,iele)) =  1
        end do
        end do

        n_vtx   =  0
        do i=1,mesh(0)%n_vtx
            if(ibuf(i) .le. 0)  cycle
            n_vtx       =  n_vtx+1
            ibuf(n_vtx) =  i
        end do

        write(str,*),n_vtx
        tec =  'zone N='//trim(adjustl(str))//',E='
        write(str,*),sec(isec)%n_ele
        tec =  trim(tec)//trim(adjustl(str))

        if(sec(isec)%is_quad) then
            str =  'FEQUADRILATERAL'
        elseif(sec(isec)%is_tri) then
            str =  'FETRIANGLE'
        else
            str =  'FEBRICK'
        end if
        tec =  trim(tec)//',zonetype='//trim(adjustl(str))//',datapacking=block'

        if(is_2d_cal) then
            if(solution .eq. 1) then
                tec =  trim(tec)//',varlocation=([3-6]=cellcentered)'
            elseif(solution .eq. 2) then
                tec =  trim(tec)//',varlocation=([3]=cellcentered)'
            elseif(solution .eq. 3) then
                tec =  trim(tec)//',varlocation=([3-5]=cellcentered)'
            elseif(solution .eq. 4) then
                tec =  trim(tec)//',varlocation=([3-7]=cellcentered)'
            else
                tec =  trim(tec)//',varlocation=([3]=cellcentered)'
            end if
        else
            if(solution .eq. 1) then
                tec =  trim(tec)//',varlocation=([4-8]=cellcentered)'
            elseif(solution .eq. 2) then
                tec =  trim(tec)//',varlocation=([4]=cellcentered)'
            else
                tec =  trim(tec)//',varlocation=([4-6]=cellcentered)'
            end if
        end if
        write(unit=10,fmt='(A)'),trim(adjustl(tec))

        do i=1,n_vtx
            write(unit=10,fmt=*),mesh(0)%xyz(1,ibuf(i))
        end do
        do i=1,n_vtx
            write(unit=10,fmt=*),mesh(0)%xyz(2,ibuf(i))
        end do
        if(.not. is_2d_cal) then
            do i=1,n_vtx
                write(unit=10,fmt=*),mesh(0)%xyz(3,ibuf(i))
            end do
        end if

        if(solution .eq. 1) then
            write(unit=10,fmt=*),fv(isec)%u(1,:)
            write(unit=10,fmt=*),fv(isec)%u(2,:)
            write(unit=10,fmt=*),fv(isec)%u(3,:)
            write(unit=10,fmt=*),fv(isec)%u(4,:)
            write(unit=10,fmt=*),fv(isec)%u(5,:)
        elseif(solution .eq. 2) then
            write(unit=10,fmt=*),fv(isec)%shock_sensor
        elseif(solution .eq. 3) then
            write(unit=10,fmt=*),fv(isec)%DDES_quality(1,:)
            write(unit=10,fmt=*),fv(isec)%DDES_quality(2,:)
            write(unit=10,fmt=*),fv(isec)%DDES_quality(3,:)
        elseif(solution .eq. 4) then
            write(unit=10,fmt='(4ES20.12)'),fv(isec)%rhs(1,:)
            write(unit=10,fmt='(4ES20.12)'),fv(isec)%rhs(2,:)
            write(unit=10,fmt='(4ES20.12)'),fv(isec)%rhs(3,:)
            write(unit=10,fmt='(4ES20.12)'),fv(isec)%rhs(4,:)
            write(unit=10,fmt='(4ES20.12)'),fv(isec)%rhs(5,:)
        else
            write(unit=10,fmt=*),fv(isec)%intermittency
        end if
        do iele=1,sec(isec)%n_ele
            do i=1,sec(isec)%npe
                call ib_search(1, n_vtx, 1, 1, ibuf, sec(isec)%n2e(i,iele), ltmp, L, R)
                if((.not. ltmp) .or. (L .ne. R))    stop 'Error: fv_wr_tec_cc.'
                v(i)=  L
            end do
            if(sec(isec)%ele_type .eq. PENTA_6) then
                write(unit=10,fmt='(8I9)'),v(1:3),v(3),v(4:6),v(6)
            else
                write(unit=10,fmt='(8I9)'),v(1:sec(isec)%npe)
            end if
        end do
    end do

    return
    end subroutine fv_wr_tec_cc
!-------------------------------------------------------------------------------
!   output cell-centered variables for FV.
!-------------------------------------------------------------------------------
    subroutine fv_wr_cgns_cc(cnt,file_name)
    use var_kind_def
    use var_cgns
    use var_fv
    use var_global
    use var_load_balance
    use var_mesh
    use var_parallel
    use var_turb, transition=>transition_model,D=>n_DDES_quality
    implicit none
    character(len=*),intent(in):: file_name
    integer(dpI),intent(in):: cnt
    character(len=80):: str(20)
    integer(dpI):: nvtx,nele,isize(3),ivtx,ID,v1,v0,iseg,i,j,m,LDA,ip,e1,e0, &
                &  isec,iele,ib
    real   (dpR):: v(100)

!   ----------------------------------------------------------------------------
!   how many vertex/elements will be output in one call of partial_write.
    nvtx=  min(2**20, n_vtx_g)
    nele=  min(2**20, n_ele_g)
!   how many vertex/elements will be output in one call of partial_write.
!   ----------------------------------------------------------------------------

    LDA =  1+n_dim
    if(.not. allocated(r_send)) then
        allocate(r_send    (LDA * nvtx  ), stat=err_mem)
    else
        if(size(r_send) .le. LDA*nvtx) then
            deallocate(r_send)
            allocate(r_send(LDA*nvtx), stat=err_mem)
        end if
    end if
    if(.not. allocated(r_recv)) then
        allocate(r_recv    (LDA * nvtx  ), stat=err_mem)
    else
        if(size(r_recv) .le. LDA*nvtx) then
            deallocate(r_recv)
            allocate(r_recv(LDA*nvtx), stat=err_mem)
        end if
    end if
    if(.not. allocated(r_recv_all)) then
        allocate(r_recv_all(nvtx*(LDA-1)), stat=err_mem)
    else
        if(size(r_recv_all) .le. (LDA-1)*nvtx) then
            deallocate(r_recv_all)
            allocate(r_recv_all((LDA-1)*nvtx), stat=err_mem)
        end if
    end if

    if(.not. allocated(ele_send)) then
        if(is_2d_cal) then
            ID  =  5
        else
            ID  =  9
        end if
        allocate(ele_send(ID, max(nele, mesh(0)%n_vtx)), stat=err_mem)
        if(myid .eq. 0) allocate(ele_recv(ID, nele), stat=err_mem)
    end if

!   ----------------------------------------------------------------------------
!   open the CGNS file.
    if(myid .eq. 0) then
        call cg_open_f(trim(adjustl(file_name)), MODE_WRITE, ifile, cg_err)
        if(is_2d_cal) then
            call cg_base_write_f(ifile,'Base',2,2,ibase,cg_err)
        else
            call cg_base_write_f(ifile,'Base',3,3,ibase,cg_err)
        end if
        isize   = (/n_vtx_g, n_elei_g, 0/)
        call cg_zone_write_f(ifile, ibase, 'zone_1', isize, unstructured, &
            &  izone, cg_err)
        call cg_sol_write_f(ifile,ibase,izone,'FlowSolution',CellCenter,iflow,cg_err)
    end if
!   open the CGNS file.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   output vertex based information.
    do iseg=1,(n_vtx_g-1)/nvtx+1
        v1  =  1+nvtx*(iseg-1)
        v0  =  min(nvtx*iseg, n_vtx_g)

!       every processor checks the number of vertex to be output in this partial write.
        ele_send=  0
        do isec=mesh(0)%sec_1,mesh(0)%sec_0
            if(.not. sec(isec)%is_int)  cycle
            do iele=1,sec(isec)%n_ele
            do i=1,sec(isec)%npe
                ele_send(1,sec(isec)%n2e(i,iele))   =  1
            end do
            end do
        end do

        m   =  0
        do ivtx=1,mesh(0)%n_vtx
            ID  =  mesh(0)%id_vtx(ivtx)
            if((ID .lt. v1) .or. (ID .gt. v0))  cycle
            if(ele_send(1,ivtx) .ne. 1) cycle
            m   =  m+1
            v(1        )=  real(ID, dpR)
            v(2:1+n_dim)=  mesh(0)%xyz(1:n_dim,ivtx)
            r_send(1+LDA*(m-1):LDA*m)   =  v(1:LDA)
        end do

        call mpi_gather(m, 1, mpi_dpI, prc_info, 1, mpi_dpI, 0, mpi_comm_world, mpi_err)

        if((myid .ne. 0) .and. (m .gt. 0))  &
            &  call mpi_send(r_send, LDA*m, mpi_dpR, 0, myid, mpi_comm_world, mpi_err)

        if(myid .ne. 0) cycle
        do ip=0,nprc-1
            if(prc_info(ip,1) .le. 0)   cycle
            if(ip .eq. 0) then
                call DCOPY(LDA*prc_info(ip,1), r_send, 1, r_recv, 1)
            else
                call mpi_recv(r_recv, LDA*prc_info(ip,1), mpi_dpR, ip, ip, &
                    &  mpi_comm_world, mpi_status, mpi_err)
            end if
!           decode the received information.
            do i=1,prc_info(ip,1)
                ID  =  nint(r_recv(1+LDA*(i-1)), dpI)
                ivtx=  ID-v1+1
                forall(j=1:LDA-1)   r_recv_all(ivtx+nvtx*(j-1)) =  r_recv(j+1+LDA*(i-1))
            end do
        end do

        call cg_coord_partial_write_f(ifile, ibase, izone, RealDouble,'CoordinateX', &
            &  v1, v0, r_recv_all(1       ), icoord, cg_err)
        call cg_coord_partial_write_f(ifile, ibase, izone, RealDouble,'CoordinateY', &
            &  v1, v0, r_recv_all(1+nvtx  ), icoord, cg_err)
        if(.not. is_2d_cal) then
        call cg_coord_partial_write_f(ifile, ibase, izone, RealDouble,'CoordinateZ', &
            &  v1, v0, r_recv_all(1+nvtx*2), icoord, cg_err)
        end if
    end do
!   output vertex based information.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   output element connectivity.
    if(is_2d_cal) then
        LDA =  5
    else
        LDA =  9
    end if
    do iseg=1,(n_elei_g-1)/nele+1
        e1  =  1+nele*(iseg-1)
        e0  =  min(iseg*nele, n_elei_g)

!       every processor checks the number of elements to be output in this partial write.
        m   =  0
        do isec=mesh(0)%sec_1,mesh(0)%sec_0
            if(.not. sec(isec)%is_int)  cycle
            do iele=1,sec(isec)%n_ele
                ID  =  sec(isec)%ID_ele_i(iele)
                if((ID .lt. e1) .or. (ID .gt. e0))  cycle
                m   =  m+1
                ele_send(1,m)   =  ID
                do i=1,sec(isec)%npe
                    ele_send(i+1,m) =  mesh(0)%id_vtx(sec(isec)%n2e(i,iele))
                end do
            end do
        end do

        call mpi_gather(m, 1, mpi_dpI, prc_info, 1, mpi_dpI, 0, mpi_comm_world, mpi_err)

        if((myid .ne. 0) .and. (m .gt. 0))  &
            &  call mpi_send(ele_send, LDA*m, mpi_dpI, 0, myid, mpi_comm_world, mpi_err)

        if(myid .ne. 0) cycle
        m   =  1
        do ip=0,nprc-1
            if(prc_info(ip,1) .le. 0)   cycle
            if(ip .eq. 0) then
                call ICOPY(LDA*prc_info(ip,1), ele_send, 1, ele_recv(1,m), 1)
            else
                call mpi_recv(ele_recv(1,m), LDA*prc_info(ip,1), mpi_dpI, ip, ip, &
                    &  mpi_comm_world, mpi_status, mpi_err)
            end if
            m   =  m+prc_info(ip,1)
        end do
        call iqsortcols(.true., 1, m-1, 1, LDA, ele_recv)

        do isec=1,n_seci_g
            if((e1 .gt. seci_g(4,isec)) .or. (e0 .lt. seci_g(3,isec)))  cycle
            i   =  max(e1, seci_g(3,isec))
            j   =  min(e0, seci_g(4,isec))
            call prepare_n2e(LDA, seci_g(2,isec), i-e1+1, j-e1+1, ele_recv)
            call cg_section_partial_write_f(ifile, ibase, izone, sec_name(isec), &
                &  seci_g(1,isec), i, j, 0, ele_recv, isection, cg_err)
            if(cg_err .ne. 0)   stop 'Error: fails to output n2e.'
        end do
    end do
!   output element connectivity.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   output cell-centered solution.
    if(cnt .eq. 1) then
        LDA     =  12
        str(1 ) =  'Density'
        str(2 ) =  'VelocityX'
        str(3 ) =  'VelocityY'
        str(4 ) =  'VelocityZ'
        str(5 ) =  'Pressure'
        str(6 ) =  'ReynoldsStressXX'
        str(7 ) =  'ReynoldsStressXY'
        str(8 ) =  'ReynoldsStressXZ'
        str(9 ) =  'ReynoldsStressYY'
        str(10) =  'ReynoldsStressYZ'
        str(11) =  'ReynoldsStressZZ'
        str(12) =  'PressureRMS'
    elseif(cnt .eq. 2) then
        if(is_RANS) then
            LDA =  6
            if(transition .ge. 1)   LDA =  7
        else
            LDA =  5
        end if
        str(1 ) =  'Density'
        str(2 ) =  'VelocityX'
        str(3 ) =  'VelocityY'
        str(4 ) =  'VelocityZ'
        str(5 ) =  'Pressure'
        str(6 ) =  'TurbulentSANuTilde'
        str(7 ) =  'Intermittency'
    elseif(cnt .eq. 3) then
        LDA     =  D
        str(1 ) =  'k_modeled'
        str(2 ) =  'k_resolved'
        str(3 ) =  'L_RANS'
        str(4 ) =  'L_LES'
    elseif(cnt .eq. 4) then
        LDA     =  1
        str(1 ) =  'Q'
        call fv_get_q
    elseif(cnt .eq. 5) then
        LDA     =  5
        str(1 ) =  'Density'
        str(2 ) =  'VelocityX'
        str(3 ) =  'VelocityY'
        str(4 ) =  'VelocityZ'
        str(5 ) =  'Pressure'
    elseif(cnt .eq. 6) then
        LDA     =  1
        str(1 ) =  'source'
    elseif(cnt .eq. 7) then
        LDA     =  2
        str(1 ) =  'Entropy_thermal'
        str(2 ) =  'Entropy_viscous'
    else
        stop 'Error: cell-centered solution not defined.'
    end if
    LDA =  LDA+1

    if(.not. allocated(r_send)) then
        allocate(r_send(LDA*nele), stat=err_mem)
    else
        if(size(r_send) .lt. LDA*nele) then
            deallocate(r_send)
            allocate(r_send(LDA*nele), stat=err_mem)
        end if
    end if
    if(.not. allocated(r_recv)) then
        allocate(r_recv(LDA*nele), stat=err_mem)
    else
        if(size(r_recv) .lt. nele*LDA) then
            deallocate(r_recv)
            allocate(r_recv(LDA*nele), stat=err_mem)
        end if
    end if

    do iseg=1,(n_elei_g-1)/nele+1
        e1  =  1+nele*(iseg-1)
        e0  =  min(iseg*nele, n_elei_g)

!       every processor checks the number of elements to be output in this partial write.
        m   =  0
        do isec=mesh(0)%sec_1,mesh(0)%sec_0
            if(.not. sec(isec)%is_int)  cycle
            do iele=1,sec(isec)%n_ele
                ID  =  sec(isec)%ID_ele_i(iele)
                if((ID .lt. e1) .or. (ID .gt. e0))  cycle
                m   =  m+1
                v(1)=  real(ID, dpR)

                if(cnt .eq. 1) then
                    v(2:6 ) =  fv(isec)%uns_ra(1:5,iele)
                    v(7:13) =  fv(isec)%stress(1:7,iele)
                elseif(cnt .eq. 2) then
                    v(2:6 ) =  fv(isec)%u     (1:5,iele)
                    if(is_RANS) v(7)=  fv(isec)%turb(1,iele)
                    if(transition .eq. transition_Menter) then
                        v(8)=  fv(isec)%turb(2,iele)
                    elseif(transition .eq. transition_Coder) then
                        v(8)=  fv(isec)%turb(2,iele)
                    elseif(transition .gt. transition_Coder) then
                        v(8)=  fv(isec)%intermittency(iele)
                    end if
                elseif(cnt .eq. 3) then
                    v(2:D+1)=  fv(isec)%DDES_quality(1:D,iele)
                elseif(cnt .eq. 4) then
                    v(2   ) =  fv(isec)%q(iele)
                elseif(cnt .eq. 5) then
                    v(2:6 ) =  fv(isec)%adj(1:5,iele)
                elseif(cnt .eq. 6) then
                    v(2   ) =  fv(isec)%caa_source(iele)
                else
                    v(2:3 ) =  fv(isec)%rhs(1:2,iele)
                end if
                r_send(1+LDA*(m-1):LDA*m)   =  v(1:LDA)
            end do
        end do

        call mpi_gather(m, 1, mpi_dpI, prc_info, 1, mpi_dpI, 0, mpi_comm_world, mpi_err)

        if((myid .ne. 0) .and. (m .gt. 0))  &
            &  call mpi_send(r_send, LDA*m, mpi_dpR, 0, myid, mpi_comm_world, mpi_err)

        if(myid .ne. 0) cycle
        m   =  1
        do ip=0,nprc-1
            if(prc_info(ip,1) .le. 0)   cycle

            ib  =  prc_info(ip,1)*LDA
            if(ip .eq. 0) then
                call DCOPY(ib, r_send, 1, r_recv(m), 1)
            else
                call mpi_recv(r_recv(m), ib, mpi_dpR, ip, ip, mpi_comm_world, &
                    &  mpi_status, mpi_err)
            end if
            m   =  m+ib
        end do
        if(m-1 .ne. LDA*(e0-e1+1))  stop 'Error: fails to gather solution, cc.'
        call dqsortcols(.true., 1, e0-e1+1, 1, LDA, r_recv)

        do m=1,LDA-1
            call DCOPY(e0-e1+1, r_recv(m+1), LDA, r_send, 1)
            if(n_elei_g .le. nele) then
                call cg_field_write_f(ifile, ibase, izone, iflow, RealDouble, &
                    &  str(m), r_send, ifield, cg_err)
            else
                call cg_field_partial_write_f(ifile, ibase, izone, iflow, RealDouble, &
                    &  str(m), e1, e0, r_send, ifield, cg_err)
            end if
            if(cg_err .ne. 0)   stop 'Error: fails to output cell-centered solution.'
        end do
    end do
!   output cell-centered solution.
!   ----------------------------------------------------------------------------

    if(myid .eq. 0) call cg_close_f(ifile, cg_err)

    if(allocated(ele_send)) deallocate(ele_send)
    if(allocated(ele_recv)) deallocate(ele_recv)

    return
    contains
!       ------------------------------------------------------------------------
!       prepare n2e.
!       ------------------------------------------------------------------------
        subroutine prepare_n2e(LDA,npe,L,R,b)
        implicit none
        integer(dpI),intent(in):: LDA,npe,L,R
        integer(dpI):: b(*),i,j,v(10)

        j   =  1
        do i=L,R
            call ICOPY(npe, b(2+LDA*(i-1)), 1, v   , 1)
            call ICOPY(npe, v             , 1, b(j), 1)
            j   =  j+npe
        end do

        return
        end subroutine prepare_n2e
    end subroutine fv_wr_cgns_cc
!-------------------------------------------------------------------------------
!   output boundary with cell-centered solutions, CGNS.
!-------------------------------------------------------------------------------
    subroutine fv_wr_boundary_cgns
    use var_kind_def
    use var_bndv, only: u_fs,ua_fs
    use var_cgns
    use var_fv
    use var_global, only: is_2d_cal,mesh_name,err_mem,n_dim
    use var_load_balance
    use var_mesh
    use var_parallel
    use var_slv, only: is_vis_cal
    use var_turb, only: transition_model,transition_AGS
    implicit none
    character(len=80):: str,sol_name(10)
    logical(dpL):: ltmp
    integer(dpI):: isec,iele,ivtx,n_vtx,n_n2e,i,isize(3),ID,e1,e0,ib,itype,nnz, &
                &  L,R,LDA
    integer(dpI),allocatable:: vtx(:),n2e(:)
    real   (dpR),allocatable:: xyz_L(:),xyz_g(:),rbuf(:)

    if(n_eleb_g .le. 0) return

    LDA         =  2
    sol_name(1) =  'Cp'
    if(is_vis_cal) then
        sol_name(2) =  'Cf'
        LDA         =  3
    end if
    if(transition_model .eq. transition_AGS) then
        sol_name(3) =  'Re_theta'
        LDA         =  4
    end if

!   record the vertex used by boundary elements.
    if(myid .eq. 0) then
        n_vtx   =  iA_n2e_b(n_eleb_g+1)-1
        allocate(vtx(n_vtx), stat=err_mem)
        call ICOPY(n_vtx, jA_n2e_b, 1, vtx, 1)
        call simplify_series(n_vtx, 1, 1, vtx)
        allocate(rbuf(max(n_vtx, n_eleb_g)), stat=err_mem)
    end if

    n_n2e   =  0
    i       =  0
    do isec=mesh(0)%sec_1,mesh(0)%sec_0
        if(.not. sec(isec)%is_bnd)  cycle
        i       =  i    +sec(isec)%n_ele
        n_n2e   =  n_n2e+sec(isec)%n_ele*sec(isec)%npe
    end do
    if(n_n2e .gt. 0)    allocate(xyz_L(max(i*LDA, (1+n_dim)*n_n2e)), stat=err_mem)

!   collect the xyz of boudary elements and output.
    n_n2e   =  0
    do isec=mesh(0)%sec_1,mesh(0)%sec_0
        if(.not. sec(isec)%is_bnd)  cycle
        do iele=1,sec(isec)%n_ele
        do i=1,sec(isec)%npe
            n_n2e                   =  n_n2e+1
            ivtx                    =  mesh(0)%ID_vtx(sec(isec)%n2e(i,iele))
            xyz_L(1+(1+n_dim)*(n_n2e-1))=  real(ivtx, dpR)
            call DCOPY(n_dim, mesh(0)%xyz(1,sec(isec)%n2e(i,iele)), 1, &
                &  xyz_L(2+(1+n_dim)*(n_n2e-1)), 1)
        end do
        end do
    end do
    call dsimplify_series(n_n2e, 1, 1+n_dim, xyz_L)
    call mpi_reduce(n_n2e, isec, 1, mpi_dpI, mpi_sum, 0, mpi_comm_world, mpi_err)
    if(myid .eq. 0) then
        allocate(xyz_g(max((1+n_dim)*isec, LDA*n_eleb_g)), stat=err_mem)
    else
        allocate(xyz_g(1))
    end if
    n_n2e   =  n_n2e*(1+n_dim)
    call r_gather(mpi_comm_world, nprc, n_n2e, xyz_L, xyz_g)
    if(myid .eq. 0) then
        call dsimplify_series(isec, 1, 1+n_dim, xyz_g)
        if(isec .ne. n_vtx) stop 'Error: fails to collect xyz_L, fv_wr_boundary_cgns.'
        do i=1,n_vtx
            if(nint(xyz_g(1+(1+n_dim)*(i-1)), dpI) .ne. vtx(i)) &
                &  stop 'Error: fails to collect xyz_L, fv_wr_boundary_cgns.'
        end do

        str =  trim(adjustl(mesh_name))//'_bnd_solution.cgns'
        call cg_open_f(trim(adjustl(str)), MODE_WRITE, ifile, cg_err)
        if(is_2d_cal) then
            call cg_base_write_f(ifile,'Base',2,2,ibase,cg_err)
        else
            call cg_base_write_f(ifile,'Base',3,3,ibase,cg_err)
        end if
        isize   = (/n_vtx, n_eleb_g , 0/)
        call cg_zone_write_f(ifile, ibase, 'zone_1', isize, unstructured, izone, cg_err)

        call DCOPY(n_vtx, xyz_g(2), 1+n_dim, rbuf, 1)
        call cg_coord_write_f(ifile, ibase, izone, RealDouble, 'CoordinateX', &
            &  rbuf, icoord, cg_err)

        call DCOPY(n_vtx, xyz_g(3), 1+n_dim, rbuf, 1)
        call cg_coord_write_f(ifile, ibase, izone, RealDouble, 'CoordinateY', &
            &  rbuf, icoord, cg_err)

        if(.not. is_2d_cal) then
        call DCOPY(n_vtx, xyz_g(4), 1+n_dim, rbuf, 1)
        call cg_coord_write_f(ifile, ibase, izone, RealDouble, 'CoordinateZ', &
            &  rbuf, icoord, cg_err)
        end if

        i   =  0
        do isec=1,n_sec_g
            itype   =  sec_g(1,isec)
            if(is_2d_cal) then
                ltmp= (itype .eq. BAR_2) .or. (itype .eq. BAR_3)
            else
                ltmp= (itype .eq. TRI_3) .or. (itype .eq. QUAD_4) .or. &
                    & (itype .eq. QUAD_8) .or. (itype .eq. QUAD_9)
            end if
            if(.not. ltmp)  cycle
            i   =  max(i, npe_ele_1st(itype)*(sec_g(4,isec)-sec_g(3,isec)+1))
        end do
        if(i .gt. 0)    allocate(n2e(i), stat=err_mem)

        ID  =  0
        nnz =  1
        e1  =  1
        do isec=1,n_sec_g
            itype   =  sec_g(1,isec)
            if(is_2d_cal) then
                ltmp= (itype .eq. BAR_2) .or. (itype .eq. BAR_3)
            else
                ltmp= (itype .eq. TRI_3) .or. (itype .eq. QUAD_4) .or. &
                    & (itype .eq. QUAD_8) .or. (itype .eq. QUAD_9)
            end if
            if(.not. ltmp)  cycle
            ID  =  ID+1
            e0  =  e1+sec_g(4,isec)-sec_g(3,isec)

            ltmp=  .false.
            do ib=1,n_bocos
                ltmp= (bocoinfo(2,ib) .eq. sec_g(3,isec)) .and. &
                    & (bocoinfo(3,ib) .eq. sec_g(4,isec))
                if(ltmp)    exit
            end do
            if(.not. ltmp)  stop 'Error: fails to map boundary to section.'

            call ICOPY(npe_ele_1st(itype)*(e0-e1+1), jA_n2e_b(nnz), 1, n2e, 1)
            do i=1,npe_ele_1st(itype)*(e0-e1+1)
                call ib_search(1, n_vtx, 1, 1, vtx, n2e(i), ltmp, L, R)
                if(.not. ltmp)  stop 'Error: fails to get vtx_b, fv_wr_boundary.'
                n2e(i)  =  L
            end do

            call cg_section_write_f(ifile, ibase, izone, bsec_name(ID), &
                &  sec_g(1,isec), e1, e0, 0, n2e, isection, cg_err)
            if(cg_err .ne. 0)   stop 'Error: fails to output boundary n2e.'

            nnz =  nnz+sec_g(2,isec)*(sec_g(4,isec)-sec_g(3,isec)+1)
            e1  =  e0+1
        end do
    end if

!   collect cell-centered solution.
    i   =  0
    do isec=mesh(0)%sec_1,mesh(0)%sec_0
        if(.not. sec(isec)%is_bnd)  cycle
        do iele=1,sec(isec)%n_ele
            i   =  i+1
            xyz_L(1+LDA*(i-1))  =  real(sec(isec)%ID_ele_b(iele), dpR)
            xyz_L(2+LDA*(i-1))  = (fv(isec)%u(5,iele)-u_fs(5))/(0.5d0*u_fs(1)*ua_fs**2)
            if(is_vis_cal) then
                if(sec(isec)%bct .eq. BCWallViscous) then
                    xyz_L(3+LDA*(i-1))  =  fv(isec)%bnd_solution(2,iele)
                else
                    xyz_L(3+LDA*(i-1))  =  0.0d0
                end if
            end if
            if(transition_model .eq. transition_AGS) then
                if(sec(isec)%bct .eq. BCWallViscous) then
                    xyz_L(4+LDA*(i-1))  =  fv(isec)%Re_theta(iele)
                else
                    xyz_L(4+LDA*(i-1))  =  0.0d0
                end if
            end if
        end do
    end do
    i   =  i*LDA
    call r_gather(mpi_comm_world, nprc, i, xyz_L, xyz_g)

    if(myid .eq. 0) then
        call dqsortcols(.true., 1, n_eleb_g, 1, LDA, xyz_g)
        call cg_sol_write_f(ifile, ibase, izone, 'FlowSolution', CellCenter, &
            &  iflow, cg_err)
        if(cg_err .ne. 0)   stop 'Error: fv_wr_boundary_cgns fails to write sol.'
        do i=1,LDA-1
            call DCOPY(n_eleb_g, xyz_g(i+1), LDA, rbuf, 1)
            call cg_field_write_f(ifile, ibase, izone, iflow, RealDouble, &
                &  trim(sol_name(i)), rbuf, ifield, cg_err)
            if(cg_err .ne. 0)   stop 'Error: fv_wr_boundary_cgns fails to write sol.'
        end do
    end if

    if(myid .eq. 0) call cg_close_f(ifile, cg_err)

    if(allocated(vtx  ))    deallocate(vtx  )
    if(allocated(n2e  ))    deallocate(n2e  )
    if(allocated(xyz_L))    deallocate(xyz_L)
    if(allocated(xyz_g))    deallocate(xyz_g)
    if(allocated(rbuf ))    deallocate(rbuf )

    return
    end subroutine fv_wr_boundary_cgns
!-------------------------------------------------------------------------------
!   module designed for tecplot output, unstructured mesh.
!-------------------------------------------------------------------------------
    module var_tec
        use var_kind_def
        use var_global, only: err_mem,n_dim
        implicit none
        private
        public:: type_tec_section,tec_get_section,tec_get_JBC,tec_output

        type type_tec_section
            integer(dpI):: n_ele    =  0
            integer(dpI):: n_vtx    =  0
            integer(dpI):: ele_type =  0
            integer(dpI):: npe      =  0
            integer(dpI):: n_sol    =  0
            integer(dpI),allocatable:: n2e(:,:)
            real   (dpR),allocatable:: xyz(:,:)
            real   (dpR),allocatable:: sol(:,:)
        end type type_tec_section

        integer(dpI),allocatable:: n2e(:)
        real   (dpR),allocatable:: xyz(:)
    contains
!   ----------------------------------------------------------------------------
!   collect data for a TEC section.
!   ----------------------------------------------------------------------------
    subroutine tec_get_section(ID,s)
    use var_mesh
    use var_parallel
    implicit none
    integer(dpI),intent(in):: ID
    logical(dpL):: ltmp
    integer(dpI):: npe,isec,iele,j,n_ele,idx,ip,n_vtx,L,R
    type(type_tec_section),optional:: s

    if((ID .le. 0) .or. (ID .gt. n_sec_g))  return
    npe =  sec_g(2,ID)

!   memory allocation.
    if(myid .eq. 0) then
        s%ele_type  =  sec_g(1,ID)
        s%npe       =  sec_g(2,ID)
        s%n_ele     =  sec_g(4,ID)-sec_g(3,ID)+1
        allocate(s%n2e(s%npe, s%n_ele), stat=err_mem)
    end if

!   ----------------------------------------------------------------------------
!   collect n2e information, parallel.
    n_ele   =  0
    do isec=mesh(0)%sec_1,mesh(0)%sec_0
        if(sec(isec)%ID_sec_g .ne. ID)  cycle
        n_ele   =  n_ele+sec(isec)%n_ele
    end do
    if(n_ele .gt. 0) then
        if(.not. allocated(n2e)) then
            allocate(n2e(npe*n_ele), stat=err_mem)
        else
            if(size(n2e) .lt. n_ele*npe) then
                deallocate(n2e)
                allocate(n2e(npe*n_ele), stat=err_mem)
            end if
        end if
    end if
    n_ele   =  0
    do isec=mesh(0)%sec_1,mesh(0)%sec_0
        if(sec(isec)%ID_sec_g .ne. ID)  cycle
        do iele=1,sec(isec)%n_ele
            n_ele   =  n_ele+1
            do j=1,npe
                n2e(j+npe*(n_ele-1))=  mesh(0)%id_vtx(sec(isec)%n2e(j,iele))
            end do
        end do
    end do
    call mpi_gather(n_ele*npe,1,mpi_dpI,prc_info(0,1),1,mpi_dpI,0,mpi_comm_world,mpi_err)
    idx =  1
    do ip=0,nprc-1
        if((ip .eq. 0) .and. (myid .eq. 0)) then
            if(prc_info(ip,1) .gt. 0)   call ICOPY(prc_info(ip,1), n2e, 1, s%n2e, 1)
        else
            if((myid .eq. ip) .and. (n_ele .gt. 0)) then
                call mpi_send(n2e, n_ele*npe, mpi_dpI, 0, ip, mpi_comm_world, mpi_err)
            end if
            if((myid .eq. 0) .and. (prc_info(ip,1) .gt. 0)) then
                call mpi_recv(s%n2e(1,idx), prc_info(ip,1), mpi_dpI, ip, ip, &
                    &  mpi_comm_world, mpi_status, mpi_err)
            end if
        end if

        if(myid .eq. 0) idx =  idx+prc_info(ip,1)/npe
    end do
!   collect n2e information, parallel.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   collect xyz, parallel.
    n_vtx   =  n_ele*npe
    if(n_vtx .gt. 0)    call simplify_series(n_vtx, 1, 1, n2e)
    call mpi_gather(n_vtx,1,mpi_dpI,prc_info(0,1),1,mpi_dpI,0,mpi_comm_world,mpi_err)

    if(myid .eq. 0) then
        j   =  prc_info(0,1)
        do ip=1,nprc-1
            j   =  j+prc_info(ip,1)
        end do
    else
        j   =  n_vtx
    end if
    if(j .gt. 0) then
        if(.not. allocated(xyz)) then
            allocate(xyz(4*j), stat=err_mem)
        else
            if(size(xyz) .lt. 4*j) then
                deallocate(xyz)
                allocate(xyz(4*j), stat=err_mem)
            end if
        end if
    end if

    do isec=mesh(0)%sec_1,mesh(0)%sec_0
        if(sec(isec)%ID_sec_g .ne. ID)  cycle
        do iele=1,sec(isec)%n_ele
        do j=1,npe
            call ib_search(1,n_vtx,1,1,n2e,mesh(0)%id_vtx(sec(isec)%n2e(j,iele)),ltmp,L,R)
            if(.not. ltmp)  stop 'Error: fails to collect xyz.'

            xyz(4*L-3)  =  real(mesh(0)%id_vtx(sec(isec)%n2e(j,iele)), dpR)
            xyz(4*L-2)  =  mesh(0)%xyz(1,sec(isec)%n2e(j,iele))
            xyz(4*L-1)  =  mesh(0)%xyz(2,sec(isec)%n2e(j,iele))
            if(n_dim .eq. 2) then
                xyz(4*L)=  0.0d0
            else
                xyz(4*L)=  mesh(0)%xyz(3,sec(isec)%n2e(j,iele))
            end if
        end do
        end do
    end do

    if(myid .eq. 0) idx =  n_vtx+1
    do ip=1,nprc-1
        if((myid .eq. ip) .and. (n_vtx .gt. 0)) then
            call mpi_send(xyz, n_vtx*4, mpi_dpR, 0, ip, mpi_comm_world, mpi_err)
        end if
        if((myid .eq. 0) .and. (prc_info(ip,1) .gt. 0)) then
            call mpi_recv(xyz(4*idx-3), 4*prc_info(ip,1), mpi_dpR, ip, ip, &
                    &  mpi_comm_world, mpi_status, mpi_err)
            idx =  idx+prc_info(ip,1)
        end if
    end do
    if(myid .ne. 0) return

    s%n_vtx =  idx-1
    call dsimplify_series(s%n_vtx, 1, 4, xyz)
    allocate(s%xyz(3, s%n_vtx), stat=err_mem)

    if(.not. allocated(n2e)) then
        allocate(n2e(s%n_vtx), stat=err_mem)
    else
        if(size(n2e) .lt. s%n_vtx) then
            deallocate(n2e)
            allocate(n2e(s%n_vtx), stat=err_mem)
        end if
    end if

    do j=1,s%n_vtx
        n2e  (    j)=  nint(xyz(4*j-3))
        s%xyz(1:3,j)=  xyz(4*j-2:4*j)
    end do
!   collect xyz, parallel.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   reset n2e.
    do iele=1,s%n_ele
        do j=1,npe
            call ib_search(1, s%n_vtx, 1, 1, n2e, s%n2e(j,iele), ltmp, L, R)
            if(.not. ltmp)  stop 'Error: fails to map n2e from global to local.'
            s%n2e(j,iele)   =  L
        end do
    end do
!   reset n2e.
!   ----------------------------------------------------------------------------

    return
    end subroutine tec_get_section
!   ----------------------------------------------------------------------------
!   collect J_BC.
!   ----------------------------------------------------------------------------
    subroutine tec_get_JBC(ID,LDA,s)
    use var_fv
    use var_mesh
    use var_parallel
    implicit none
    integer(dpI),intent(in):: ID,LDA
    integer(dpI):: isec,iele,i,idx,n_ele,ip
    type(type_tec_section),optional:: s

    if(myid .eq. 0) then
        s%n_sol =  LDA
        allocate(s%sol(LDA, s%n_ele), stat=err_mem)
    end if

    n_ele   =  0
    do isec=mesh(0)%sec_1,mesh(0)%sec_0
        if(sec(isec)%ID_sec_g .eq. ID)  n_ele   =  n_ele+sec(isec)%n_ele
    end do
    call mpi_gather(n_ele*LDA,1,mpi_dpI,prc_info,1,mpi_dpI,0,mpi_comm_world,mpi_err)

    if(n_ele .gt. 0) then
        if(.not. allocated(xyz)) then
            allocate(xyz(LDA*n_ele), stat=err_mem)
        else
            if(size(xyz) .lt. LDA*n_ele) then
                deallocate(xyz)
                allocate(xyz(LDA*n_ele), stat=err_mem)
            end if
        end if
    end if
    i   =  0
    do isec=mesh(0)%sec_1,mesh(0)%sec_0
        if(sec(isec)%ID_sec_g .ne. ID)  cycle
        do iele=1,sec(isec)%n_ele
            i   =  i+1
            xyz(1+LDA*(i-1):LDA*i)  =  fv(isec)%J_bc(1:LDA,iele)
        end do
    end do

    idx =  1
    do ip=0,nprc-1
        if((ip .eq. 0) .and. (myid .eq. 0)) then
            if(prc_info(ip,1) .gt. 0)   call DCOPY(prc_info(ip,1), xyz, 1, s%sol, 1)
        else
            if((myid .eq. ip) .and. (n_ele .gt. 0)) then
                call mpi_send(xyz, n_ele*LDA, mpi_dpR, 0, ip, mpi_comm_world, mpi_err)
            end if
            if((myid .eq. 0) .and. (prc_info(ip,1) .gt. 0)) then
                call mpi_recv(s%sol(1,idx), prc_info(ip,1), mpi_dpR, ip, ip, &
                    &  mpi_comm_world, mpi_status, mpi_err)
            end if
        end if

        if(myid .eq. 0) idx =  idx+prc_info(ip,1)/LDA
    end do

    return
    end subroutine tec_get_JBC
!   ----------------------------------------------------------------------------
!   output section.
!   ----------------------------------------------------------------------------
    subroutine tec_output(n_sec,sec,file_name,is_output_sol)
    use var_cgns
    use var_global, only: mesh_name
    use var_parallel
    implicit none
    character(len=*),intent(in):: file_name
    logical(dpL),intent(in):: is_output_sol
    integer(dpI),intent(in):: n_sec
    type(type_tec_section),intent(in):: sec(*)
    character(len=1000):: str,str1
    integer(dpI):: isec,iele,LDA,i,j

    if((myid .ne. 0) .or. (n_sec .le. 0))   return

    LDA =  0
    if(is_output_sol) then
        do isec=1,n_sec
            LDA =  max(LDA, sec(isec)%n_sol)
        end do
    end if

    str =  trim(adjustl(mesh_name))//'_'//trim(adjustl(file_name))//'.dat'
    open(unit=10,file=trim(adjustl(str)))

    str =  'variables="CoordinateX","CoordinateY","CoordinateZ"'
    do i=1,LDA
        write(str1,fmt=*),i
        str =  trim(adjustl(str))//',"JBC_'//trim(adjustl(str1))//'"'
    end do
    write(unit=10,fmt='(A)'),trim(adjustl(str))

    do isec=1,n_sec
        str =  'zone N='
        write(str1,fmt=*),sec(isec)%n_vtx
        str =  trim(adjustl(str))//trim(adjustl(str1))

        write(str1,fmt=*),sec(isec)%n_ele
        str =  trim(adjustl(str))//',E='//trim(adjustl(str1))

        if(sec(isec)%ele_type .eq. QUAD_4) then
            str1=  'FEQUADRILATERAL'
        elseif(sec(isec)%ele_type .eq. BAR_2) then
            str1=  'FELINESEG'
        elseif(sec(isec)%ele_type .eq. TRI_3) then
            str1=  'FETRIANGLE'
        else
            stop 'Error: ele_type not supported.'
        end if
        str =  trim(str)//',zonetype='//trim(adjustl(str1))//',datapacking=block'

        if(LDA .gt. 0) then
            write(str1,fmt=*),LDA+3
            if(LDA .eq. 1) then
                str =  trim(adjustl(str))//',varlocation=([4]=cellcentered)'
            else
                str =  trim(adjustl(str))//',varlocation=([4-'//trim(adjustl(str1))//']=cellcentered)'
            end if
        end if

        write(unit=10,fmt='(A)'),trim(adjustl(str))

        write(unit=10,fmt=*),sec(isec)%xyz(1,:)
        write(unit=10,fmt=*),sec(isec)%xyz(2,:)
        write(unit=10,fmt=*),sec(isec)%xyz(3,:)

        do i=1,LDA
            if(i .le. sec(isec)%n_sol) then
                write(unit=10,fmt=*),sec(isec)%sol(i,:)
            else
                write(unit=10,fmt=*),(0.0d0, j=1,sec(isec)%n_ele)
            end if
        end do

        do iele=1,sec(isec)%n_ele
            write(unit=10,fmt='(4I8)'),sec(isec)%n2e(:,iele)
        end do
    end do
    close(10)

    return
    end subroutine tec_output
    end module var_tec
