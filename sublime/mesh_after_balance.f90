!-------------------------------------------------------------------------------
!   distribute elements, after load balance.
!-------------------------------------------------------------------------------
    subroutine set_section
    use var_kind_def
    use var_cgns, only: npe_ele
    use var_global, only: err_mem
    use var_load_balance
    use var_mesh
    use var_parallel
    use var_per_bnd
    implicit none
    logical(dpL):: ltmp
    integer(dpI):: iele,i,j,k,ip,sec_info(3,100),isec,itype,nsec,npe
    integer(dpI):: v(100),L,R,ID_ele_g

!   ----------------------------------------------------------------------------
!   collect information about idx, ele_type, n_n2e and n2e to their destiny.
    prc_info=  0
    do iele=1,n_ele_L
        ip              =  part(iele)
        prc_info(ip,1)  =  prc_info(ip,1)+3+iA_n2e(iele+1)-iA_n2e(iele)
    end do

    call mpi_alltoall(prc_info(0,1), 1, mpi_dpI, prc_info(0,3), 1, mpi_dpI, &
        &  mpi_comm_world, mpi_err)
    i   =  prc_info(0,1)
    j   =  prc_info(0,3)
    do ip=1,nprc-1
        i   =  i+prc_info(ip,1)
        j   =  j+prc_info(ip,3)
    end do
    i   =  0
    do ip=0,nprc-1
        prc_info(ip,2)  =  i
        i   =  i+prc_info(ip,1)
    end do
    if(allocated(n2e_i_have)) then
        if(i .gt. size(n2e_i_have)) then
            deallocate(n2e_i_have)
            allocate(n2e_i_have(i), stat=err_mem)
        end if
    else
        if(i .gt. 0)    allocate(n2e_i_have(i), stat=err_mem)
    end if
    if(allocated(n2e_i_need)) then
        if(j .gt. size(n2e_i_need)) then
            deallocate(n2e_i_need)
            allocate(n2e_i_need(j), stat=err_mem)
        end if
    else
        if(j .gt. 0)    allocate(n2e_i_need(j), stat=err_mem)
    end if

    do iele=1,n_ele_L
        ip  =  part(iele)
        k   =  prc_info(ip,2)

        k               =  k+1
        n2e_i_have(k)   =  iele+elei_prc(1,myid)-1
        k               =  k+1
        n2e_i_have(k)   =  ele_type(iele)
        k               =  k+1
        n2e_i_have(k)   =  iA_n2e(iele+1)-iA_n2e(iele)
        do j=iA_n2e(iele),iA_n2e(iele+1)-1
            k               =  k+1
            n2e_i_have(k)   =  jA_n2e(j)
        end do
        prc_info(ip,2)  =  k
    end do
    call collect_data_parallel(prc_info(0,1), prc_info(0,2), n2e_i_have, prc_info(0,3), &
        &  prc_info(0,4), n2e_i_need)
    nnz_n2e_i_need  =  prc_info(nprc-1,3)+prc_info(nprc-1,4)-1

    i       =  0
    nsec    =  0
    sec_info=  0
    do while(i .lt. nnz_n2e_i_need)
        i       =  i+1
        iele    =  n2e_i_need(i)
        i       =  i+1
        itype   =  n2e_i_need(i)
        i       =  i+1
        npe     =  n2e_i_need(i)
        i       =  i+npe

        if((iele .lt. seci_g(6,1)) .or. (iele .gt. seci_g(7,n_seci_g))) &
            &  stop 'Error: ele_idx out of range.'
        do j=1,n_seci_g
            if((iele .ge. seci_g(6,j)) .and. (iele .le. seci_g(7,j)))   exit
        end do
        ltmp=  .true.
        do isec=1,nsec
            if((j .eq. sec_info(1,isec)) .and. (itype .eq. sec_info(2,isec))) then
                ltmp=  .false.
                sec_info(3,isec)=  sec_info(3,isec)+1
                exit
            end if
        end do
        if(ltmp) then
            nsec=  nsec+1
            sec_info(1:3,nsec)  = (/j, itype, 1/)
        end if
    end do
    if(nsec .gt. max_sec_lev)   stop 'Error: max_sec_lev too small.'

    mesh(0)%sec_1   =  1
    mesh(0)%sec_0   =  nsec
    ID_ele_g = 1
    do isec=1,nsec
        sec(isec)%ID_sec_i  =  sec_info(1,isec)
        sec(isec)%ele_type  =  sec_info(2,isec)
        sec(isec)%n_ele     =  sec_info(3,isec)
        sec(isec)%npe       =  npe_ele(sec_info(2,isec))
        allocate(sec(isec)%n2e     (sec(isec)%npe, sec(isec)%n_ele))
        allocate(sec(isec)%id_ele_i(               sec(isec)%n_ele))
!#ifdef SW_slave
!        allocate(sec(isec)%ID_ele_g(               sec(isec)%n_ele))
!        do iele=1,sec(isec)%n_ele
!            sec(isec)%ID_ele_g(iele) = ID_ele_g
!            ID_ele_g = ID_ele_g+1
!        end do
!#endif
    end do

    sec_info(3,1:nsec)  =  0
    i   =  0
    do while(i .lt. nnz_n2e_i_need)
        i       =  i+1
        iele    =  n2e_i_need(i)
        i       =  i+1
        itype   =  n2e_i_need(i)
        i       =  i+1
        npe     =  n2e_i_need(i)
        do j=1,npe
            i   =  i+1
            v(j)=  n2e_i_need(i)
        end do

        do j=1,n_seci_g
            if((iele .ge. seci_g(6,j)) .and. (iele .le. seci_g(7,j)))   exit
        end do
        do isec=1,nsec
            if((j .eq. sec_info(1,isec)) .and. (itype .eq. sec_info(2,isec))) then
                sec_info(3,isec)=  sec_info(3,isec)+1
                npe             =  sec(isec)%npe
                sec(isec)%id_ele_i( sec_info(3,isec))   =  iele
                sec(isec)%n2e(1:npe,sec_info(3,isec))   =  v(1:npe)
                exit
            end if
        end do
    end do
!   collect information about idx, ele_type, n_n2e and n2e to their destiny.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   sort out which elements have already been recorded.
    n_ele_i_have=  0
    do isec=mesh(0)%sec_1,mesh(0)%sec_0
        n_ele_i_have=  n_ele_i_have+sec(isec)%n_ele
    end do
    if(allocated(ele_i_have)) then
        if(n_ele_i_have .gt. size(ele_i_have)) then
            deallocate(ele_i_have)
            allocate(ele_i_have(n_ele_i_have), stat=err_mem)
        end if
    else
        if(n_ele_i_have .gt. 0) allocate(ele_i_have(n_ele_i_have), stat=err_mem)
    end if
    n_ele_i_have=  0
    do isec=mesh(0)%sec_1,mesh(0)%sec_0
    do iele=1,sec(isec)%n_ele
        n_ele_i_have=  n_ele_i_have+1
        ele_i_have(n_ele_i_have)=  sec(isec)%id_ele_i(iele)
    end do
    end do
    call iqsort(.true., 1, n_ele_i_have, ele_i_have)
!   sort out which elements have already been recorded.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   get to know which vertex is used by the internal element.
    k   =  0
    i   =  0
    do while(i .lt. nnz_n2e_i_need)
!       ID, ele_type, n_n2e, n2e.
        i   =  i+3
        do j=1,n2e_i_need(i)
            i   =  i+1
            k   =  k+1
            n2e_i_need(k)   =  n2e_i_need(i)
        end do
    end do
    call simplify_series(k, 1, 1, n2e_i_need)
    if(allocated(vtx_i_need)) then
        if(k .gt. size(vtx_i_need)) then
            deallocate(vtx_i_need)
            allocate(vtx_i_need(k))
        end if
    else
        if(k .gt. 0)    allocate(vtx_i_need(k), stat=err_mem)
    end if
    call ICOPY(k, n2e_i_need, 1, vtx_i_need, 1)
    n_vtx_i_need=  k
!   get to know which vertex is used by the internal element.
!   ----------------------------------------------------------------------------

    call get_ele_ghost(.false., .true.)

!   ----------------------------------------------------------------------------
!   record the donor vertex due to periodic matching.
    if(n_vtx_pair_per .gt. 0) then
        n_vtx_i_have=  n_vtx_i_need
        if(allocated(vtx_i_have)) then
            if(size(vtx_i_have) .lt. n_vtx_i_have) then
                deallocate(vtx_i_have)
                allocate(vtx_i_have(n_vtx_i_have), stat=err_mem)
            end if
        else
            allocate(vtx_i_have(n_vtx_i_have), stat=err_mem)
        end if
        call ICOPY(n_vtx_i_have, vtx_i_need, 1, vtx_i_have, 1)
        n_vtx_i_need=  0

        k   =  size(vtx_i_need)
        call iqsortcols(.true., 1, n_vtx_pair_per, 2, 5, vtx_pair_per)
        do i=1,n_vtx_i_have
            call ib_search(1,n_vtx_pair_per,2,5,vtx_pair_per,vtx_i_have(i),ltmp,L,R)
            if(.not. ltmp)  cycle
            do j=L,R
                n_vtx_i_need=  n_vtx_i_need+1
                if(n_vtx_i_need .gt. k) stop 'Error: vtx_i_need is small.'
                vtx_i_need(n_vtx_i_need)=  vtx_pair_per(1,j)
            end do
        end do
        call simplify_series(n_vtx_i_need, 1, 1, vtx_i_need)
        call get_ele_per(.false., .false.)
    end if
!   record the donor vertex due to periodic matching.
!   ----------------------------------------------------------------------------

    call check_ghost_element

    call get_structured_stencil

    return
    end subroutine set_section
!-------------------------------------------------------------------------------
!   delete some ghost elements not needed.
!-------------------------------------------------------------------------------
    subroutine check_ghost_element
    use var_kind_def
    use var_cgns
    use var_global, only: n_dim,err_mem
    use var_load_balance
    use var_mesh
    use var_parallel, only: is_include_corner
    use var_per_bnd
    implicit none
    logical(dpL):: ltmp
    integer(dpI):: isec,iele,i,j,k,M,L,R,ivtx,n_ele,npe

!   ----------------------------------------------------------------------------
!   record the vertex used by the internal element.
    nnz_n2e_i_have  =  0
    do isec=mesh(0)%sec_1,mesh(0)%sec_0
        if(.not. sec(isec)%is_int)  cycle
        nnz_n2e_i_have  =  nnz_n2e_i_have+sec(isec)%n_ele*npe_ele_1st(sec(isec)%ele_type)
    end do
    if(allocated(n2e_i_have)) then
        if(nnz_n2e_i_have .gt. size(n2e_i_have)) then
            deallocate(n2e_i_have)
            allocate(n2e_i_have(nnz_n2e_i_have), stat=err_mem)
        end if
    else
        allocate(n2e_i_have(max(nnz_n2e_i_have,1)), stat=err_mem)
    end if
    nnz_n2e_i_have  =  0
    do isec=mesh(0)%sec_1,mesh(0)%sec_0
        if(.not. sec(isec)%is_int)  cycle
        do iele=1,sec(isec)%n_ele
            call ICOPY(npe_ele_1st(sec(isec)%ele_type), sec(isec)%n2e(1,iele), 1, &
                &  n2e_i_have(nnz_n2e_i_have+1), 1)
            nnz_n2e_i_have  =  nnz_n2e_i_have+npe_ele_1st(sec(isec)%ele_type)
        end do
    end do
    call simplify_series(nnz_n2e_i_have, 1, 1, n2e_i_have)
!   record the vertex used by the internal element.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   check the ghost element.
    n_ele       =  n_ele_ghost
    n_ele_ghost =  0
    do i=1,n_ele
        npe =  npe_ele_1st(ele_type_ghost(i))
        M   =  0
        do j=iA_ghost(i),iA_ghost(i)+npe-1
            call ib_search(1, nnz_n2e_i_have, 1, 1, n2e_i_have, jA_ghost(j), ltmp, L, R)
            if(ltmp)    M   =  M+1
        end do
        if(is_include_corner) then
            ltmp=  M .gt. 0
        else
            if(n_dim .eq. 2) then
                ltmp=  M .ge. 2
            else
                ltmp=  M .ge. 3
            end if
        end if
        if(.not. ltmp)  cycle
        n_ele_ghost =  n_ele_ghost+1
        ele_type_ghost(n_ele_ghost) =  ele_type_ghost(i)
        ele_idx_ghost (n_ele_ghost) =  ele_idx_ghost (i)
        j   =  iA_ghost(i)
        k   =  iA_ghost(n_ele_ghost)
        M   =  iA_ghost(i+1)-iA_ghost(i)
        call ICOPY(M, jA_ghost(j), 1, jA_ghost(k), 1)
        iA_ghost(n_ele_ghost+1) =  iA_ghost(n_ele_ghost)+M
    end do
!   check the ghost element.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   check the periodic ghost element.
    if(n_vtx_pair_per .le. 0)   return
    call iqsortcols(.true., 1, n_vtx_pair_per, 1, 5, vtx_pair_per)
    if(.not. allocated(is_corner_per))  allocate(is_corner_per(n_ele_per))
    n_ele       =  n_ele_per
    n_ele_per   =  0
    do i=1,n_ele
        npe =  npe_ele_1st(ele_type_per(i))
        M   =  0
        do j=iA_per(i),iA_per(i)+npe-1
!           try to get the receiver vertex_ID.
            call ib_search(1,n_vtx_pair_per,1,5,vtx_pair_per,abs(jA_per(j)),ltmp,L,R)
            if(.not. ltmp)  cycle
            ivtx= -1
            do k=L,R
                if((ele_per(2,i) .eq. vtx_pair_per(3,k)) .and. &
                 & (ele_per(3,i) .eq. vtx_pair_per(4,k)) .and. &
                 & (ele_per(4,i) .eq. vtx_pair_per(5,k))) then
                    ivtx=  vtx_pair_per(2,k)
                    exit
                end if
            end do
            if(ivtx .le. 0) cycle

!           to check whether the receiver vertex is already recorded.
            call ib_search(1, nnz_n2e_i_have, 1, 1, n2e_i_have, ivtx, ltmp, L, R)
            if(ltmp)    M   =  M+1
        end do

        if(is_include_corner) then
            ltmp=  M .gt. 0
        else
            if(n_dim .eq. 2) then
                ltmp=  M .ge. 2
            else
                ltmp=  M .ge. 3
            end if
        end if
        if(.not. ltmp)  cycle

        n_ele_per   =  n_ele_per+1
        if(n_dim .eq. 2) then
            is_corner_per(n_ele_per)=  M .eq. 1
        else
            is_corner_per(n_ele_per)=  M .le. 2
        end if
        if(n_ele_per .eq. i)    cycle
        ele_per   (:,n_ele_per) =  ele_per   (:,i)
        ele_type_per(n_ele_per) =  ele_type_per(i)

        j   =  iA_per(i)
        k   =  iA_per(n_ele_per)
        M   =  iA_per(i+1)-iA_per(i)
        call ICOPY(M, jA_per(j), 1, jA_per(k), 1)
        iA_per(n_ele_per+1) =  iA_per(n_ele_per)+M
    end do
!   check the periodic ghost element.
!   ----------------------------------------------------------------------------

    return
    end subroutine check_ghost_element
!-------------------------------------------------------------------------------
!   set the ghost sections after load balance.
!-------------------------------------------------------------------------------
    subroutine set_ghost_section
    use var_kind_def
    use var_cgns, only: npe_ele
    use var_global, only: err_mem
    use var_mesh
    use var_load_balance
    use var_parallel
    implicit none
    logical(dpL):: ltmp
    integer(dpI):: i,j,k,L,R,M,isec,iele,isec_g,sec_info(5,1000),n_sec,npe
    integer(dpI),allocatable:: ghost(:,:)

    M   =  n_ele_ghost+n_ele_per+n_ss
    if(M .le. 0)    return
    allocate(ghost(7,M), stat=err_mem)
    if(err_mem .ne. 0)  stop 'Error: fails to allocate memory, set_ghost_section.'
    ghost   =  0

!   ----------------------------------------------------------------------------
!   record the ID and per_path of ghost elements.
    M   =  0
    do i=1,n_ele_ghost
        M   =  M+1
        ghost(1  ,M)=  ele_idx_ghost(i)
        ghost(3:5,M)=  0
        ghost(6  ,M)=  1
        ghost(7  ,M)=  i
    end do
    do i=1,n_ele_per
        M   =  M+1
        ghost(1  ,M)=  ele_per(1  ,i)
        ghost(3:5,M)=  ele_per(2:4,i)
        ghost(6  ,M)=  2
        ghost(7  ,M)=  i
    end do
    do i=1,n_ss
        M   =  M+1
        ghost(1  ,M)=  ele_ID_ss( i)
        ghost(3:5,M)=  per_ss(1:3,i)
        ghost(6  ,M)=  3
        ghost(7  ,M)=  i
    end do

    do i=1,M
        j   =  ghost(1,i)
        if((j .le. 0) .or. (j .gt. n_elei_g))   stop 'Error: ID of ghost ele is wrong.'
        do isec=1,n_seci_g
            if((j .ge. seci_g(6,isec)) .and. (j .le. seci_g(7,isec)))   exit
        end do
        ghost(2,i)  =  isec

        if(ghost(3,i) .gt. 0)   cycle
        call get_isec_iele(ghost(1,i), isec, iele)
        if(isec .gt. 0) stop 'Error: ghost element should be internal element.'
    end do
!   record the ID and per_path of ghost elements.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   delete redundant record.
    call iqsortcols(.true., 1, M, 1, 7, ghost)
    i   =  1
    do while(i .le. M)
        k   =  i
        do j=i+1,M
            if(ghost(1,j) .eq. ghost(1,i)) then
                k   =  j
            else
                exit
            end if
        end do
        do j=i+1,k
        do L=i,j-1
            if(ghost(1,L) .le. 0)   cycle
            ltmp= (ghost(3,L) .eq. ghost(3,j)) .and. (ghost(4,L) .eq. ghost(4,j)) .and. &
                & (ghost(5,L) .eq. ghost(5,j))
            if(ltmp) then
                ghost(1:7,j)= -1
                exit
            end if
        end do
        end do

        i   =  k+1
    end do
    j   =  M
    M   =  0
    do i=1,j
        if(ghost(1,i) .le. 0)   cycle
        M           =  M+1
        ghost(1:7,M)=  ghost(1:7,i)
    end do
!   delete redundant record.
!   ----------------------------------------------------------------------------

    call iqsortcols(.true., 1, M, 2, 7, ghost)

    i       =  1
    n_sec   =  0
    do while(i .le. M)
        k   =  i
        do j=i+1,M
            if(ghost(2,j) .ne. ghost(2,i)) then
                exit
            else
                k   =  j
            end if
        end do
        n_sec   =  n_sec+1
        sec_info(1,n_sec)   =  ghost(2,i)
        sec_info(2,n_sec)   =  k-i+1
        sec_info(3,n_sec)   =  i
        sec_info(4,n_sec)   =  k
        call iqsortcols(.true., i, k, 1, 7, ghost)

        i   =  k+1
    end do

    do i=1,n_sec
        isec=  mesh(0)%sec_0+1
        mesh(0)%sec_0   =  isec

        isec_g  =  sec_info(1,i)
        sec(isec)%ID_sec_i  =  isec_g
        sec(isec)%ele_type  =  seci_g(1,isec_g)
        sec(isec)%n_ele     =  sec_info(2,i)
        sec(isec)%npe       =  npe_ele(seci_g(1,isec_g))
        sec(isec)%is_int    =  .false.
        sec(isec)%is_bnd    =  .false.
        sec(isec)%is_ghost  =  .true.
        allocate(sec(isec)%n2e     (sec(isec)%npe, sec(isec)%n_ele))
        allocate(sec(isec)%id_ele_i(               sec(isec)%n_ele))
        allocate(sec(isec)%per_path(            3, sec(isec)%n_ele))
        sec(isec)%per_path  =  0
        npe =  sec(isec)%npe
        do j=sec_info(3,i),sec_info(4,i)
            iele=  j-sec_info(3,i)+1
            sec(isec)%ID_ele_i(    iele)=  ghost(1  ,j)
            sec(isec)%per_path(1:3,iele)=  ghost(3:5,j)
            L   =  ghost(7,j)

            if(ghost(6,j) .eq. 1) then
                R   =  iA_ghost(L)
                sec(isec)%n2e(1:npe,iele)   =  jA_ghost(R:R+npe-1)
            elseif(ghost(6,j) .eq. 2) then
                R   =  iA_per(L)
                sec(isec)%n2e(1:npe,iele)   =  jA_per  (R:R+npe-1)
            else
                R   =  iA_ss(L)
                sec(isec)%n2e(1:npe,iele)   =  jA_ss   (R:R+npe-1)
            end if
        end do
    end do

    if(allocated(iA_ghost      ))   deallocate(iA_ghost)
    if(allocated(jA_ghost      ))   deallocate(jA_ghost)
    if(allocated(ele_type_ghost))   deallocate(ele_type_ghost)
    if(allocated(ele_idx_ghost ))   deallocate(ele_idx_ghost)
    if(allocated(ele_per       ))   deallocate(ele_per)
    if(allocated(iA_per        ))   deallocate(iA_per)
    if(allocated(jA_per        ))   deallocate(jA_per)
    if(allocated(ele_type_per  ))   deallocate(ele_type_per)
    if(allocated(iA_ss         ))   deallocate(iA_ss)
    if(allocated(jA_ss         ))   deallocate(jA_ss)
    if(allocated(ele_ID_ss     ))   deallocate(ele_ID_ss)
    if(allocated(ele_type_ss   ))   deallocate(ele_type_ss)
    if(allocated(per_ss        ))   deallocate(per_ss)
    if(allocated(ghost         ))   deallocate(ghost)

    return
    end subroutine set_ghost_section
!-------------------------------------------------------------------------------
!   set p2p.
!-------------------------------------------------------------------------------
    subroutine set_p2p
    use var_kind_def
    use var_global, only: err_mem
    use var_load_balance
    use var_mesh
    use var_parallel
    implicit none
    integer(dpI):: M,isec,iele,i,j,k,N,ip,ip2p,L
    integer(dpI),allocatable:: ghost(:,:)

!   ----------------------------------------------------------------------------
!   record the ghost element myid needs.
    M   =  0
    do isec=mesh(0)%sec_1,mesh(0)%sec_0
        if(sec(isec)%is_ghost)  M   =  M+sec(isec)%n_ele
    end do
    if(allocated(ele_i_need)) then
        if(M .gt. size(ele_i_need)) then
            deallocate(ele_i_need)
            allocate(ele_i_need(M), stat=err_mem)
        end if
    else
        allocate(ele_i_need(max(M,1)), stat=err_mem)
    end if

    M   =  0
    do isec=mesh(0)%sec_1,mesh(0)%sec_0
        if(.not. sec(isec)%is_ghost)    cycle
        do i=1,sec(isec)%n_ele
            M               =  M+1
            ele_i_need(M)   =  sec(isec)%ID_ele_i(i)
        end do
    end do
    call simplify_series(M, 1, 1, ele_i_need)
    if(M .gt. 0) then
        allocate(ghost(2,M), stat=err_mem)
        do i=1,M
            ghost(1,i)  =  ele_i_need(i)
            ghost(2,i)  =  0
        end do
    end if
!   record the ghost element myid needs.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   get the ip of the ghost element.
    prc_info=  0
    do i=1,M
        iele=  ele_i_need(i)
        if(iele .ge. elei_prc(1,nprc-1)) then
            ip  =  nprc-1
        else
            ip  = (iele-1)/(n_elei_g/nprc)
        end if
        prc_info(ip,1)  =  prc_info(ip,1)+1
    end do
    call mpi_alltoall(prc_info(0,1), 1, mpi_dpI, prc_info(0,3), 1, mpi_dpI, &
        &  mpi_comm_world, mpi_err)
    n   =  prc_info(0,3)
    do ip=1,nprc-1
        n   =  n+prc_info(ip,3)
    end do
    if(allocated(ele_i_have)) then
        if(n .gt. size(ele_i_have)) then
            deallocate(ele_i_have)
            allocate(ele_i_have(n), stat=err_mem)
        end if
    else
        allocate(ele_i_have(max(n,1)), stat=err_mem)
    end if
    call collect_data_parallel(prc_info(0,1), prc_info(0,2), ele_i_need, prc_info(0,3), &
        &  prc_info(0,4), ele_i_have)
    do i=1,prc_info(nprc-1,3)+prc_info(nprc-1,4)-1
        iele=  ele_i_have(i)
        if((iele .lt. elei_prc(1,myid)) .or. (iele .gt. elei_prc(2,myid))) &
            &  stop 'Error: fails to set p2p.'
        ele_i_have(i)   =  part(iele-elei_prc(1,myid)+1)
    end do
    prc_info(0:nprc-1,1)=  prc_info(0:nprc-1,3)
    call collect_data_parallel(prc_info(0,1), prc_info(0,2), ele_i_have, prc_info(0,3), &
        &  prc_info(0,4), ele_i_need)
!   get the ip of the ghost element.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   set p2p_r.
    do i=1,M
        ghost(2,i)  =  ele_i_need(i)
    end do
    if(M .gt. 0)    call iqsortcols(.true., 1, M, 2, 2, ghost)

    i   =  1
    do while(i .le. M)
        k   =  i
        do j=i+1,M
            if(ghost(2,j) .eq. ghost(2,i)) then
                k   =  j
            else
                exit
            end if
        end do
        call iqsortcols(.true., i, k, 1, 2, ghost)

        i   =  k+1
    end do

    prc_info=  0
    do i=1,M
        prc_info(ghost(2,i),1)  =  prc_info(ghost(2,i),1)+1
    end do
    call mpi_alltoall(prc_info(0,1), 1, mpi_dpI, prc_info(0,3), 1, mpi_dpI, &
        &  mpi_comm_world, mpi_err)
    ip2p=  0
    k   =  1
    do ip=0,nprc-1
!       myid receives data from ip.
        i   =  prc_info(ip,1)
!       myid sends    data to   ip.
        j   =  prc_info(ip,3)
        if((i .le. 0) .and. (j .le. 0)) cycle
        if((i .gt. 0) .neqv. (j .gt. 0))    stop 'Error: send&recv not match.'
        if((myid .eq. ip) .and. (i .ne. j)) stop 'Error: send&recv on myid not match.'

        ip2p=  ip2p+1
        mesh(0)%p2p_0       =  ip2p
        p2p(ip2p)%ip_remote =  ip

        if(j .gt. 0)    allocate(p2p(ip2p)%ID_ele_send(2,j), stat=err_mem)

        if(i .gt. 0) then
            p2p(ip2p)%n_ele_recv=  i
            allocate(p2p(ip2p)%ID_ele_recv(i), stat=err_mem)
            do j=1,i
                p2p(ip2p)%ID_ele_recv(j)=  ghost(1,j+k-1)
            end do
            k   =  k+i
        end if
    end do
!   set p2p_r.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   set p2p_s.
    do i=1,M
        ele_i_need(i)   =  ghost(1,i)
    end do
    call mpi_alltoall(prc_info(0,1), 1, mpi_dpI, prc_info(0,3), 1, mpi_dpI, &
        &  mpi_comm_world, mpi_err)
    j   =  prc_info(0,3)
    do ip=1,nprc-1
        j   =  j+prc_info(ip,3)
    end do
    if(allocated(n2e_i_have)) then
        if(size(n2e_i_have) .lt. j) then
            deallocate(n2e_i_have)
            allocate(n2e_i_have(max(j, 1)), stat=err_mem)
        end if
    else
        allocate(n2e_i_have(max(j, 1)), stat=err_mem)
    end if

    call collect_data_parallel(prc_info(0,1), prc_info(0,2), ele_i_need, prc_info(0,3), &
        &  prc_info(0,4), n2e_i_have)

    do ip=0,nprc-1
        if(prc_info(ip,3) .le. 0)   cycle
        ip2p=  0
        do j=mesh(0)%p2p_1,mesh(0)%p2p_0
            if(p2p(j)%ip_remote .eq. ip) then
                ip2p=  j
                exit
            end if
        end do
        if(ip2p .le. 0) stop 'Error: fails to set p2p.'

        do i=prc_info(ip,4),prc_info(ip,3)+prc_info(ip,4)-1
            call get_isec_iele(n2e_i_have(i), isec, L)
            if(isec .lt. 0) stop 'Error: fails to find local element for p2p_s.'
            p2p(ip2p)%n_ele_send=  p2p(ip2p)%n_ele_send+1
            p2p(ip2p)%ID_ele_send(1:2,p2p(ip2p)%n_ele_send) = (/isec, L/)
        end do

        if(p2p(ip2p)%n_ele_send .ne. prc_info(ip,3))    stop 'Error: fails to set p2p_s.'
    end do
!   set p2p_s.
!   ----------------------------------------------------------------------------

    if(allocated(ghost))    deallocate(ghost)

    return
    end subroutine set_p2p
!-------------------------------------------------------------------------------
!   set the vertex.
!-------------------------------------------------------------------------------
    subroutine set_vertex
    use var_kind_def
    use var_global, only: err_mem,n_dim
    use var_load_balance
    use var_mesh
    use var_parallel
    use var_per_bnd
    implicit none
    logical(dpL):: ltmp
    integer(dpI):: isec,iele,L,R,M,N,i,j,k
    real   (dpR):: x(3)
    integer(dpI),allocatable:: ghost(:,:)

!   ----------------------------------------------------------------------------
!   record the vertex used and delete duplicate record.
    M   =  0
    do isec=mesh(0)%sec_1,mesh(0)%sec_0
        if(.not. sec(isec)%is_bnd)  M   =  M+sec(isec)%n_ele*sec(isec)%npe
    end do
    allocate(ghost(5,M), stat=err_mem)
    ghost   =  0

    M   =  0
    do isec=mesh(0)%sec_1,mesh(0)%sec_0
        if(sec(isec)%is_bnd)    cycle
        do iele=1,sec(isec)%n_ele
        do i=1,sec(isec)%npe
            M   =  M+1
            ghost(1  ,M)=  abs(sec(isec)%n2e(i,iele))
            if(sec(isec)%is_int) then
                ghost(2:4,M)=  0
            else
                ghost(2:4,M)=  sec(isec)%per_path(1:3,iele)
            end if
        end do
!       write(6,fmt='(4I6,A2,3I6)'),sec(isec)%n2e(:,iele),':',sec(isec)%per_path(:,iele)
        end do
    end do

    call iqsortcols(.true., 1, M, 1, 5, ghost)
    i   =  1
    do while(i .le. M)
        k   =  i
        do j=i+1,M
            if(ghost(1,j) .ne. ghost(1,i)) then
                exit
            else
                k   =  j
            end if
        end do

        do j=i+1,k
        do n=i,j-1
            if(ghost(1,n) .le. 0)   cycle
            ltmp= (ghost(2,n) .eq. ghost(2,j)) .and. (ghost(3,n) .eq. ghost(3,j)) .and. &
                & (ghost(4,n) .eq. ghost(4,j))
            if(ltmp) then
                ghost(1:5,j)= -1
                exit
            end if
        end do
        end do

        i   =  k+1
    end do
    j   =  M
    M   =  0
    do i=1,j
        if(ghost(1,i) .lt. 0)   cycle
        M           =  M+1
        ghost(1:5,M)=  ghost(1:5,i)
    end do
!   record the vertex used and delete duplicate record.
!   ----------------------------------------------------------------------------

    call iqsortcols(.true., 1, M, 2, 5, ghost)
    mesh(0)%n_vtx   =  0
    do i=1,M
        if(ghost(2,i) .le. 0) then
            mesh(0)%n_vtx   =  mesh(0)%n_vtx+1
        else
            exit
        end if
    end do
    call iqsortcols(.true., 1, mesh(0)%n_vtx, 1, 5, ghost)
    allocate(mesh(0)%ID_vtx(mesh(0)%n_vtx), stat=err_mem)
    do i=1,mesh(0)%n_vtx
        ghost(5,i)  =  i
        mesh(0)%ID_vtx(i)   =  ghost(1,i)
    end do

!   ----------------------------------------------------------------------------
!   to check whether the periodic vertex has its mapping recorded by myid.
    if(M .ge. mesh(0)%n_vtx)    call iqsortcols(.true., mesh(0)%n_vtx+1, M, 1, 5, ghost)
    if(n_vtx_pair_per .gt. 0)   call iqsortcols(.true.,1,n_vtx_pair_per,1,5,vtx_pair_per)
    do i=mesh(0)%n_vtx+1,M
        k   = -1
        call ib_search(1, n_vtx_pair_per, 1, 5, vtx_pair_per, ghost(1,i), ltmp, L, R)
        if(ltmp) then
            do j=L,R
                if((ghost(2,i) .eq. vtx_pair_per(3,j)) .and. &
                  &(ghost(3,i) .eq. vtx_pair_per(4,j)) .and. &
                  &(ghost(4,i) .eq. vtx_pair_per(5,j))) then
                    k   =  vtx_pair_per(2,j)
                    exit
                end if
            end do
        end if
        if(k .le. 0)    cycle
!       if(k .le. 0)    stop 'Error: fails to get the receiver of periodic vertex.'

!       check whether myid already has the receiver vertex.
        call ib_search(1, mesh(0)%n_vtx, 1, 1, mesh(0)%ID_vtx, k, ltmp, L, R)
        if(ltmp)    ghost(5,i)  =  L
    end do
!   to check whether the periodic vertex has its mapping recorded by myid.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   reocord the ID and per_path of ghost vertex.
    mesh(0)%n_vtx_per   =  0
    do i=1,M
        if(ghost(5,i) .le. 0)   mesh(0)%n_vtx_per   =  mesh(0)%n_vtx_per+1
    end do
    if(mesh(0)%n_vtx_per .gt. 0) then
        allocate(mesh(0)%ID_vtx_per(mesh(0)%n_vtx_per), stat=err_mem)
        allocate(mesh(0)%per_path(3,mesh(0)%n_vtx_per), stat=err_mem)
        mesh(0)%per_path=  0
    end if
    mesh(0)%n_vtx_per   =  0
    do i=1,M
        if(ghost(5,i) .gt. 0)   cycle
        mesh(0)%n_vtx_per   =  mesh(0)%n_vtx_per+1
        mesh(0)%ID_vtx_per(  mesh(0)%n_vtx_per) =  ghost(1  ,i)
        mesh(0)%per_path(1:3,mesh(0)%n_vtx_per) =  ghost(2:4,i)
        ghost(5,i)  =  mesh(0)%n_vtx_per+mesh(0)%n_vtx
    end do
!   reocord the ID and per_path of ghost vertex.
!   ----------------------------------------------------------------------------

    call iqsortcols(.true., 1, M, 1, 5, ghost)

!   ----------------------------------------------------------------------------
!   map vertex of ghost element from global to local.
    do isec=mesh(0)%sec_1,mesh(0)%sec_0
        if(sec(isec)%is_bnd)    cycle
        do iele=1,sec(isec)%n_ele
        do i=1,sec(isec)%npe
            call ib_search(1, M, 1, 5, ghost, abs(sec(isec)%n2e(i,iele)), ltmp, L, R)
            if(.not. ltmp)  stop 'Error: set_vertex fails.'
            sec(isec)%n2e(i,iele)   =  0
            if(sec(isec)%is_int) then
                do j=L,R
                    if((ghost(2,j) .eq. 0) .and. (ghost(3,j) .eq. 0) .and. &
                      &(ghost(4,j) .eq. 0)) then
                        sec(isec)%n2e(i,iele)   =  ghost(5,j)
                        exit
                    end if
                end do
            else
                do j=L,R
                    if((ghost(2,j) .eq. sec(isec)%per_path(1,iele)) .and. &
                      &(ghost(3,j) .eq. sec(isec)%per_path(2,iele)) .and. &
                      &(ghost(4,j) .eq. sec(isec)%per_path(3,iele))) then
                        sec(isec)%n2e(i,iele)   =  ghost(5,j)
                        exit
                    end if
                end do
            end if
            if(sec(isec)%n2e(i,iele) .le. 0)    stop 'Error: fails to set the vertex ID.'
        end do
        end do
    end do
!   map vertex of ghost element from global to local.
!   ----------------------------------------------------------------------------

    if(allocated(ghost))    deallocate(ghost)

!   ----------------------------------------------------------------------------
!   get vertex coordinates.
    n_xyz_i_need=  mesh(0)%n_vtx+mesh(0)%n_vtx_per
    if(allocated(xyz_i_need)) then
        if(size(xyz_i_need) .lt. n_xyz_i_need) then
            deallocate(xyz_i_need)
            allocate(xyz_i_need(n_xyz_i_need), stat=err_mem)
        end if
    else
        allocate(xyz_i_need(n_xyz_i_need), stat=err_mem)
    end if

    call ICOPY(mesh(0)%n_vtx    , mesh(0)%ID_vtx    , 1, xyz_i_need                 , 1)
    if(mesh(0)%n_vtx_per .gt. 0) then
    call ICOPY(mesh(0)%n_vtx_per, mesh(0)%ID_vtx_per, 1, xyz_i_need(1+mesh(0)%n_vtx), 1)
    end if

    call simplify_series(n_xyz_i_need, 1, 1, xyz_i_need)

    call get_xyz_i_need
    allocate(mesh(0)%xyz(n_dim,mesh(0)%n_vtx+mesh(0)%n_vtx_per), stat=err_mem)
    do i=1,mesh(0)%n_vtx
        call ib_search(1, n_xyz_i_need, 1, 1, xyz_i_need, mesh(0)%ID_vtx(i), ltmp, L, R)
        if((.not. ltmp) .or. (L .ne. R))    stop 'Error: set_vertex fails.'
        mesh(0)%xyz(1:n_dim,i)  =  xyz_need(1+n_dim*(L-1):n_dim*L)
    end do
    do i=1,mesh(0)%n_vtx_per
        call ib_search(1,n_xyz_i_need,1,1,xyz_i_need,mesh(0)%ID_vtx_per(i),ltmp,L,R)
        if((.not. ltmp) .or. (L .ne. R))    stop 'Error: set_vertex fails.'

        x(1:n_dim)  =  xyz_need(1+n_dim*(L-1):n_dim*L)
        j           =  i+mesh(0)%n_vtx
        if(mesh(0)%per_path(1,i) .le. 0) then
            mesh(0)%xyz(1:n_dim,j)  =  x(1:n_dim)
        else
            call per_rot_vec(0, mesh(0)%per_path(1,i), x, mesh(0)%xyz(1,j))
        end if
    end do
!   get vertex coordinates.
!   ----------------------------------------------------------------------------

    if(allocated(xyz_i_need))   deallocate(xyz_i_need)
    if(allocated(xyz_i_have))   deallocate(xyz_i_have)
    if(allocated(xyz_need  ))   deallocate(xyz_need  )
    if(allocated(xyz_have  ))   deallocate(xyz_have  )
    if(allocated(ghost     ))   deallocate(ghost     )

    return
    end subroutine set_vertex
!-------------------------------------------------------------------------------
!   distribute boundary.
!-------------------------------------------------------------------------------
    subroutine set_bnd_section
    use var_kind_def
    use var_cgns
    use var_global, only: is_2d_cal,err_mem
    use var_load_balance
    use var_mesh
    use var_parallel
    implicit none
    logical(dpL):: ltmp
    integer(dpI):: isec,iele,i,j,k,nfac,ele,v(99),vr(4),npe,ifac,n_matched,matched(10), &
                &  nsec,sec_info(3,100),L,R,iele_all
    integer(dpI),allocatable:: fac(:,:),b_info(:,:)

    if(n_eleb_g .le. 0) return

!   ----------------------------------------------------------------------------
!   record the faces of the internal elements.
    nfac=  0
    do isec=mesh(0)%sec_1,mesh(0)%sec_0
        if(.not. sec(isec)%is_int)  cycle
        ele =  sec(isec)%ele_type
        nfac=  nfac+sec(isec)%n_ele*nface_ele(ele)
    end do
    allocate(fac(5,nfac), stat=err_mem)
    k   =  huge(1)-1
    fac =  k

    nfac    =  0
    iele_all=  0
    do isec=mesh(0)%sec_1,mesh(0)%sec_0
        if(.not. sec(isec)%is_int)  cycle
        ele =  sec(isec)%ele_type
        npe =  npe_ele_1st(ele)
        do iele=1,sec(isec)%n_ele
            iele_all=  iele_all+1
            do i=1,npe
                v(i)=  mesh(0)%id_vtx(sec(isec)%n2e(i,iele))
            end do

            if((ele .eq. TRI_3)) then
                fac(1:5,nfac+1) = (/v(1), v(2), k, k, iele_all/)
                fac(1:5,nfac+2) = (/v(2), v(3), k, k, iele_all/)
                fac(1:5,nfac+3) = (/v(3), v(1), k, k, iele_all/)
                nfac    =  nfac+3
            elseif((ele .eq. QUAD_4) .or. (ele .eq. QUAD_8) .or. (ele .eq. QUAD_9)) then
                fac(1:5,nfac+1) = (/v(1), v(2), k, k, iele_all/)
                fac(1:5,nfac+2) = (/v(2), v(3), k, k, iele_all/)
                fac(1:5,nfac+3) = (/v(3), v(4), k, k, iele_all/)
                fac(1:5,nfac+4) = (/v(4), v(1), k, k, iele_all/)
                nfac    =  nfac+4
            elseif(ele .eq. TETRA_4) then
                fac(1:5,nfac+1) = (/v(1), v(2), v(3), k, iele_all/)
                fac(1:5,nfac+2) = (/v(1), v(2), v(4), k, iele_all/)
                fac(1:5,nfac+3) = (/v(2), v(3), v(4), k, iele_all/)
                fac(1:5,nfac+4) = (/v(1), v(3), v(4), k, iele_all/)
                nfac    =  nfac+4
            elseif(ele .eq. PYRA_5) then
                fac(1:5,nfac+1) = (/v(1), v(2), v(3), v(4), iele_all/)
                fac(1:5,nfac+2) = (/v(1), v(2), v(5), k, iele_all/)
                fac(1:5,nfac+3) = (/v(2), v(3), v(5), k, iele_all/)
                fac(1:5,nfac+4) = (/v(3), v(4), v(5), k, iele_all/)
                fac(1:5,nfac+5) = (/v(4), v(1), v(5), k, iele_all/)
                nfac    =  nfac+5
            elseif(ele .eq. PENTA_6) then
                fac(1:5,nfac+1) = (/v(1), v(2), v(3), k, iele_all/)
                fac(1:5,nfac+2) = (/v(4), v(5), v(6), k, iele_all/)
                fac(1:5,nfac+3) = (/v(1), v(2), v(5), v(4), iele_all/)
                fac(1:5,nfac+4) = (/v(2), v(3), v(6), v(5), iele_all/)
                fac(1:5,nfac+5) = (/v(1), v(3), v(6), v(4), iele_all/)
                nfac    =  nfac+5
            elseif((ele .eq. HEXA_8) .or. (ele .eq. HEXA_20) .or. (ele .eq. HEXA_27)) then
                fac(1:5,nfac+1) = (/v(1), v(2), v(3), v(4), iele_all/)
                fac(1:5,nfac+2) = (/v(5), v(6), v(7), v(8), iele_all/)
                fac(1:5,nfac+3) = (/v(1), v(5), v(8), v(4), iele_all/)
                fac(1:5,nfac+4) = (/v(2), v(6), v(7), v(3), iele_all/)
                fac(1:5,nfac+5) = (/v(1), v(2), v(6), v(5), iele_all/)
                fac(1:5,nfac+6) = (/v(4), v(3), v(7), v(8), iele_all/)
                nfac    =  nfac+6
            else
                stop 'Error: element type not supported.'
            end if
        end do
    end do
    do ifac=1,nfac
        call iqsort(.true., 1, 4, fac(1,ifac))
    end do
    call iqsortcols(.true., 1, nfac, 1, 5, fac)
!   record the faces of the internal elements.
!   ----------------------------------------------------------------------------

    if(allocated(n2e_i_have)) then
        if(size(n2e_i_have) .lt. iele_all) then
            deallocate(n2e_i_have)
            allocate(n2e_i_have(iele_all), stat=err_mem)
        end if
    else
        allocate(n2e_i_have(iele_all), stat=err_mem)
    end if
    n2e_i_have  =  0

!   ----------------------------------------------------------------------------
!   record the boundary patch related to internal element.
    if(n_eleb_g .gt. 0) allocate(b_info(3,n_eleb_g), stat=err_mem)
    b_info  =  0
    do iele=1,n_eleb_g
        i   =  iA_n2e_b(iele)
        j   =  iA_n2e_b(iele+1)-1
        npe =  j-i+1
        if(is_2d_cal) then
            if(npe .eq. 2) then
                ele =  BAR_2
            elseif(npe .eq. 3) then
                ele =  BAR_3
            else
                stop 'Error: boundary element type not supported.'
            end if
        else
            if(npe .eq. 3) then
                ele =  TRI_3
            elseif(npe .eq. 4) then
                ele =  QUAD_4
            elseif(npe .eq. 8) then
                ele =  QUAD_8
            elseif(npe .eq. 9) then
                ele =  QUAD_9
            else
                stop 'Error: boundary element type not supported.'
            end if
        end if
        v(1:npe)=  jA_n2e_b(i:j)
        vr      =  k
        call ICOPY(npe_ele_1st(ele), v, 1, vr, 1)

        call iqsort(.true., 1, 4, vr)
        call find_col_in_matrix(5, 1, 4, 1, nfac, fac, vr, n_matched, matched)
        if(n_matched .le. 0)    cycle
        if(n_matched .gt. 1)    stop 'Error: multiple matching of ele_b&ele.'
        b_info(1:2,iele)= (/ele, bct(1,iele)/)

        n2e_i_have(fac(5,matched(1)))   =  1
    end do
!   record the boundary patch related to internal element.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   record the sections formed by the above elements.
    nsec=  0
    do iele=1,n_eleb_g
        if(b_info(1,iele) .le. 0)   cycle
        ltmp=  .false.
        do isec=1,nsec
            ltmp= (b_info(1,iele) .eq. sec_info(1,isec)) .and. &
                & (b_info(2,iele) .eq. sec_info(2,isec))
            if(ltmp) then
                sec_info(3,isec)=  sec_info(3,isec)+1
                exit
            end if
        end do
        if(.not. ltmp) then
            nsec=  nsec+1
            sec_info(1:3,nsec)  = (/b_info(1:2,iele), 1/)
        end if
    end do
    do i=1,nsec
        isec=  mesh(0)%sec_0+1
        mesh(0)%sec_0   =  isec

        sec(isec)%is_int    =  .false.
        sec(isec)%is_ghost  =  .false.
        sec(isec)%is_bnd    =  .true.
        sec(isec)%ID_bnd    =  sec_info(2,i)
        sec(isec)%ele_type  =  sec_info(1,i)
        sec(isec)%n_ele     =  sec_info(3,i)
        sec(isec)%npe       =  npe_ele(sec_info(1,i))
        allocate(sec(isec)%n2e     (sec(isec)%npe, sec(isec)%n_ele))
        allocate(sec(isec)%ID_ele_b(               sec(isec)%n_ele))
        sec(isec)%n_ele     =  0
        sec(isec)%ID_ele_b  =  0
    end do

    do iele=1,n_eleb_g
        if(b_info(1,iele) .le. 0)   cycle
        ltmp=  .false.
        do isec=mesh(0)%sec_1,mesh(0)%sec_0
            if(.not. sec(isec)%is_bnd)  cycle
            ltmp= (b_info(1,iele) .eq. sec(isec)%ele_type) .and. &
                & (b_info(2,iele) .eq. sec(isec)%ID_bnd)
            if(ltmp) then
                sec_info(3,isec)=  sec_info(3,isec)+1
                exit
            end if
        end do
        if(.not. ltmp)  stop 'Error: fails to map ele_b to sec.'
        sec(isec)%bct   =  bct(2,iele)
        npe =  sec(isec)%npe
        i   =  iA_n2e_b(iele)
        call ICOPY(npe, jA_n2e_b(i), 1, v, 1)
        do i=1,npe
            call ib_search(1, mesh(0)%n_vtx, 1, 1, mesh(0)%id_vtx, v(i), ltmp, L, R)
            if(.not. ltmp)  stop 'Error: fails to find local id_vtx for ele_b.'
            v(i)=  L
        end do
        sec(isec)%n_ele =  sec(isec)%n_ele+1
        call ICOPY(npe, v, 1, sec(isec)%n2e(1,sec(isec)%n_ele), 1)
        sec(isec)%ID_ele_b(sec(isec)%n_ele) =  iele
    end do
    call set_bnd_ID(0)
!   record the sections formed by the above elements.
!   ----------------------------------------------------------------------------

    if(allocated(fac   ))   deallocate(fac)
    if(allocated(b_info))   deallocate(b_info)

    return
    end subroutine set_bnd_section
!-------------------------------------------------------------------------------
!   reorder the internal elements.
!-------------------------------------------------------------------------------
    subroutine reorder_internal_element
    use var_kind_def
    use var_global, only: err_mem
    use var_load_balance
    use var_mesh
    use var_parallel, only: myid
    implicit none
    logical(dpL):: ltmp
    integer(dpI):: isec,iele,i,ip2p,j,L,R,N,iele_all,isr,v(2)

    i   =  0
    do isec=mesh(0)%sec_1,mesh(0)%sec_0
        if(sec(isec)%is_int)    i   =  max(i, sec(isec)%n_ele)
    end do
    if(allocated(ele_i_have)) then
        if(size(ele_i_have) .lt. i) then
            deallocate(ele_i_have)
            allocate(ele_i_have(max(i, 1)), stat=err_mem)
        end if
    else
        allocate(ele_i_have(max(i, 1)), stat=err_mem)
    end if

    iele_all=  0
    do isec=mesh(0)%sec_1,mesh(0)%sec_0
        if(.not. sec(isec)%is_int)  cycle
        N   =  3+sec(isec)%npe
        if(N*sec(isec)%n_ele .gt. size(n2e_i_need)) stop 'Error: n2e_i_need has wrong size.'
        ele_i_have  =  0

!       element used in the parallel communication.
        do ip2p=mesh(0)%p2p_1,mesh(0)%p2p_0
        do i=1,p2p(ip2p)%n_ele_send
            if(p2p(ip2p)%ID_ele_send(1,i) .ne. isec)    cycle
            ele_i_have(p2p(ip2p)%ID_ele_send(2,i))  =  1
        end do
        end do

!       element used in the boundary condition.
        do i=1,sec(isec)%n_ele
            iele_all=  iele_all+1
            if(n2e_i_have(iele_all) .gt. 0) ele_i_have(i)   =  1
        end do

        j   =  0
        do i=1,sec(isec)%n_ele
            if(ele_i_have(i) .le. 0)    cycle
            j   =  j+1
            n2e_i_need(N*(j-1)+1)       =  i
            n2e_i_need(N*(j-1)+2)       =  j
            n2e_i_need(N*(j-1)+3)       =  sec(isec)%ID_ele_i(i)
            n2e_i_need(N*(j-1)+4:N*j)   =  sec(isec)%n2e(1:sec(isec)%npe,i)
        end do
        sec(isec)%n_ele_b   =  j
        if(j .le. 0)    cycle

        do i=1,sec(isec)%n_ele
            if(ele_i_have(i) .gt. 0)    cycle
            j   =  j+1
            n2e_i_need(N*(j-1)+1)       =  i
            n2e_i_need(N*(j-1)+2)       =  j
            n2e_i_need(N*(j-1)+3)       =  sec(isec)%ID_ele_i(i)
            n2e_i_need(N*(j-1)+4:N*j)   =  sec(isec)%n2e(1:sec(isec)%npe,i)
        end do
        if(j .ne. sec(isec)%n_ele)  stop 'Error: fails to move element_b to the top.'

        call ICOPY(sec(isec)%n_ele, n2e_i_need(3), N, sec(isec)%ID_ele_i, 1)
        do i=1,sec(isec)%n_ele
            sec(isec)%n2e(1:sec(isec)%npe,i)=  n2e_i_need(N*(i-1)+4:N*i)
        end do

        call iqsortcols(.true., 1, sec(isec)%n_ele, 1, N, n2e_i_need)
        do ip2p=mesh(0)%p2p_1,mesh(0)%p2p_0
        do i=1,p2p(ip2p)%n_ele_send
            if(p2p(ip2p)%ID_ele_send(1,i) .ne. isec)    cycle
            call ib_search(1, sec(isec)%n_ele, 1, N, n2e_i_need, &
                &  p2p(ip2p)%ID_ele_send(2,i), ltmp, L, R)
            if((.not. ltmp) .or. (L .ne. R))    stop 'Error: fails to move eb to the top.'
            p2p(ip2p)%ID_ele_send(2,i)  =  n2e_i_need(N*(L-1)+2)
        end do
        end do
    end do

    do isec=mesh(0)%sec_1,mesh(0)%sec_0
        if(.not. sec(isec)%is_ghost)    cycle
        allocate(sec(isec)%ID_recv(2,sec(isec)%n_ele))
        sec(isec)%ID_recv   =  0

        do iele=1,sec(isec)%n_ele
            i   =  sec(isec)%id_ele_i(iele)

            do isr=mesh(0)%p2p_1,mesh(0)%p2p_0
                call ib_search(1, p2p(isr)%n_ele_recv, 1, 1, p2p(isr)%id_ele_recv, &
                    &  i, ltmp, L, R)
                if(ltmp) then
                    sec(isec)%ID_recv(1:2,iele) = (/isr, L/)
                    exit
                end if
            end do
            if(sec(isec)%ID_recv(1,iele) .le. 0)    stop 'Error: fails to map ele to recv.'
            if(p2p(isr)%ip_remote .eq. myid) then
                v(1)=  p2p(isr)%id_ele_send(1,L)
                v(2)=  p2p(isr)%id_ele_send(2,L)
                if(.not. sec(v(1))%is_int)  stop 'Error: ID of ghost and donor not match.'
                if(i .ne. sec(v(1))%id_ele_i(v(2))) stop 'Error: ID of ghost and donor not match.'
            end if
        end do
    end do

    return
    end subroutine reorder_internal_element
!-------------------------------------------------------------------------------
!   clean memory after load balance.
!-------------------------------------------------------------------------------
    subroutine load_balance_clean
    use var_load_balance
    implicit none

    if(allocated(iA_n2e        ))   deallocate(iA_n2e)
    if(allocated(jA_n2e        ))   deallocate(jA_n2e)
    if(allocated(ele_type      ))   deallocate(ele_type)
    if(allocated(elei_prc      ))   deallocate(elei_prc)
    if(allocated(iA_e2n        ))   deallocate(iA_e2n)
    if(allocated(jA_e2n        ))   deallocate(jA_e2n)
    if(allocated(vtx_prc       ))   deallocate(vtx_prc)
!   if(allocated(iA_n2e_b      ))   deallocate(iA_n2e_b)
!   if(allocated(jA_n2e_b      ))   deallocate(jA_n2e_b)
    if(allocated(bct           ))   deallocate(bct)
    if(allocated(iA_ghost      ))   deallocate(iA_ghost)
    if(allocated(jA_ghost      ))   deallocate(jA_ghost)
    if(allocated(ele_type_ghost))   deallocate(ele_type_ghost)
    if(allocated(ele_idx_ghost ))   deallocate(ele_idx_ghost)
    if(allocated(vtx_i_need    ))   deallocate(vtx_i_need)
    if(allocated(vtx_i_have    ))   deallocate(vtx_i_have)
    if(allocated(e2n_to_send   ))   deallocate(e2n_to_send)
    if(allocated(e2n_i_need    ))   deallocate(e2n_i_need)
    if(allocated(ele_i_need    ))   deallocate(ele_i_need)
    if(allocated(ele_i_have    ))   deallocate(ele_i_have)
    if(allocated(n2e_i_need    ))   deallocate(n2e_i_need)
    if(allocated(n2e_i_have    ))   deallocate(n2e_i_have)
    if(allocated(is_corner_per ))   deallocate(is_corner_per)
    if(allocated(ele_per       ))   deallocate(ele_per)
    if(allocated(iA_per        ))   deallocate(iA_per)
    if(allocated(jA_per        ))   deallocate(jA_per)
    if(allocated(ele_type_per  ))   deallocate(ele_type_per)
    if(allocated(vtxdist       ))   deallocate(vtxdist)
    if(allocated(vwgt          ))   deallocate(vwgt)
    if(allocated(part          ))   deallocate(part)
    if(allocated(ele_prc       ))   deallocate(ele_prc)
    if(allocated(tpwgts        ))   deallocate(tpwgts)
    if(allocated(iA_ss         ))   deallocate(iA_ss)
    if(allocated(jA_ss         ))   deallocate(jA_ss)
    if(allocated(ele_ID_ss     ))   deallocate(ele_ID_ss)
    if(allocated(ele_type_ss   ))   deallocate(ele_type_ss)
    if(allocated(per_ss        ))   deallocate(per_ss)
    if(allocated(xyz_i_need    ))   deallocate(xyz_i_need)
    if(allocated(xyz_i_have    ))   deallocate(xyz_i_have)
    if(allocated(xyz_need      ))   deallocate(xyz_need)
    if(allocated(xyz_have      ))   deallocate(xyz_have)

    return
    end subroutine load_balance_clean
!-------------------------------------------------------------------------------
!   move the face-ghost element to the top in the p2p send/recv list.
!-------------------------------------------------------------------------------
    subroutine set_p2p_element_order(lev)
    use var_kind_def
    use var_global, only: err_mem
    use var_mesh, only: mesh,sec,p2p,max_p2p
    use var_parallel
    implicit none
    integer(dpI),intent(in):: lev
    integer(dpI):: isr,sR,eR,iele,im,ip_remote,isec

    if(.false.) then
        do isr=mesh(lev)%p2p_1,mesh(lev)%p2p_0
            p2p(isr)%n_ele_send_face=  p2p(isr)%n_ele_send
            p2p(isr)%n_ele_recv_face=  p2p(isr)%n_ele_recv
        end do
        return
    end if

    do isr=mesh(lev)%p2p_1,mesh(lev)%p2p_0
        if(p2p(isr)%n_ele_send .gt. 0) then
            allocate(p2p(isr)%isend(p2p(isr)%n_ele_send), stat=err_mem)
            p2p(isr)%isend  =  0
        end if
        if(p2p(isr)%n_ele_recv .gt. 0) then
            allocate(p2p(isr)%irecv(p2p(isr)%n_ele_recv), stat=err_mem)
            p2p(isr)%irecv  =  0
        end if
    end do

!   structured stencil.
    do im=1+mesh(lev)%n_mortar_b,mesh(lev)%n_mortar_ss
        sR  =  mesh(lev)%mortar_structured_stencil(5,im)
        eR  =  mesh(lev)%mortar_structured_stencil(6,im)
        if(.not. sec(sR)%is_ghost)  cycle
        isr =  sec(sR)%ID_recv(1,eR)
        iele=  sec(sR)%ID_recv(2,eR)
        p2p(isr)%irecv(iele)=  1
    end do

!   mortar stencil.
    do im=1+mesh(lev)%n_mortar_ss,mesh(lev)%n_mortar
        sR  =  mesh(lev)%mortar_LR(3,im)
        eR  =  mesh(lev)%mortar_LR(4,im)
        if(.not. sec(sR)%is_ghost)  cycle

        isr =  sec(sR)%ID_recv(1,eR)
        iele=  sec(sR)%ID_recv(2,eR)
        p2p(isr)%irecv(iele)=  1
    end do

    do isr=mesh(lev)%p2p_1,mesh(lev)%p2p_0
        if(p2p(isr)%n_ele_recv .le. 0)  cycle
        im  =  0
        do iele=1,p2p(isr)%n_ele_recv
            if(p2p(isr)%irecv(iele) .le. 0) cycle
            im                  =  im+1
            p2p(isr)%irecv(iele)=  im
        end do
    end do

!   ----------------------------------------------------------------------------
!   data exchange.
    mpi_nreq=  0
    do isr=mesh(lev)%p2p_1,mesh(lev)%p2p_0
        ip_remote   =  p2p(isr)%ip_remote
        if(myid .eq. ip_remote) cycle

        if(p2p(isr)%n_ele_recv .gt. 0) then
            mpi_nreq=  mpi_nreq+1
            if(mpi_nreq .gt. max_p2p)   stop 'Error: you are working with a bad LB.'
            call mpi_isend(p2p(isr)%irecv, p2p(isr)%n_ele_recv, mpi_dpI, ip_remote, &
                &  myid, mpi_comm_world, mpi_req(mpi_nreq), mpi_err)
        end if

        if(p2p(isr)%n_ele_send .gt. 0) then
            mpi_nreq=  mpi_nreq+1
            if(mpi_nreq .gt. max_p2p)   stop 'Error: you are working with a bad LB.'
            call mpi_irecv(p2p(isr)%isend, p2p(isr)%n_ele_send, mpi_dpI, ip_remote, &
                &  ip_remote, mpi_comm_world, mpi_req(mpi_nreq), mpi_err)
        end if
    end do
!   data exchange.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   data exchange on myid.
    do isr=mesh(lev)%p2p_1,mesh(lev)%p2p_0
        ip_remote   =  p2p(isr)%ip_remote
        if(myid .ne. ip_remote) cycle
        if(p2p(isr)%n_ele_send .ne. p2p(isr)%n_ele_recv) &
            &  stop 'Error: s&r not match on myid.'

        if(p2p(isr)%n_ele_send .gt. 0)  call ICOPY(p2p(isr)%n_ele_send, p2p(isr)%irecv, &
            &  1, p2p(isr)%isend, 1)
    end do
!   data exchange on myid.
!   ----------------------------------------------------------------------------

    call mpi_waitall(mpi_nreq, mpi_req, mpi_sta, mpi_err)

!   ----------------------------------------------------------------------------
!   reorder the ele_send and ele_recv.
    do isr=mesh(lev)%p2p_1,mesh(lev)%p2p_0
        if(p2p(isr)%n_ele_recv .gt. 0) then
            im  =  0
            do iele=1,p2p(isr)%n_ele_recv
                if(p2p(isr)%irecv(iele) .le. 0) cycle
                im  =  im+1
            end do
            p2p(isr)%n_ele_recv_face=  im
            do iele=1,p2p(isr)%n_ele_recv
                if(p2p(isr)%irecv(iele) .gt. 0) cycle
                im                  =  im+1
                p2p(isr)%irecv(iele)=  im
            end do
            if(im .ne. p2p(isr)%n_ele_recv) stop 'Error: set_p2p_element_order fails.'

            call permute_list(1,p2p(isr)%n_ele_recv,p2p(isr)%irecv,p2p(isr)%id_ele_recv)
!           print*,'RECV:',p2p(isr)%n_ele_recv,'|',p2p(isr)%n_ele_recv_face
        end if

        if(p2p(isr)%n_ele_send .gt. 0) then
            im  =  0
            do iele=1,p2p(isr)%n_ele_send
                if(p2p(isr)%isend(iele) .le. 0) cycle
                im  =  im+1
            end do
            p2p(isr)%n_ele_send_face=  im
            do iele=1,p2p(isr)%n_ele_send
                if(p2p(isr)%isend(iele) .gt. 0) cycle
                im                  =  im+1
                p2p(isr)%isend(iele)=  im
            end do
            if(im .ne. p2p(isr)%n_ele_send) stop 'Error: set_p2p_element_order fails.'

            call permute_list(2,p2p(isr)%n_ele_send,p2p(isr)%isend,p2p(isr)%id_ele_send)
!           print*,'SEND:',p2p(isr)%n_ele_send,'|',p2p(isr)%n_ele_send_face
        end if
    end do
!   reorder the ele_send and ele_recv.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   update the information of ghost section.
    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(.not. sec(isec)%is_ghost)    cycle
        do im=1,sec(isec)%n_ele
            isr =  sec(isec)%id_recv(1,im)
            iele=  sec(isec)%id_recv(2,im)
            sec(isec)%id_recv(2,im) =  p2p(isr)%irecv(iele)
        end do
    end do
!   update the information of ghost section.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   check.
    do im=1+mesh(lev)%n_mortar_b,mesh(lev)%n_mortar_ss
        sR  =  mesh(lev)%mortar_structured_stencil(5,im)
        eR  =  mesh(lev)%mortar_structured_stencil(6,im)
        if(.not. sec(sR)%is_ghost)  cycle
        isr =  sec(sR)%ID_recv(1,eR)
        iele=  sec(sR)%ID_recv(2,eR)
        if(iele .gt. p2p(isr)%n_ele_recv_face)  stop 'Error: set_p2p_ele_order fails.'
    end do

    do im=1+mesh(lev)%n_mortar_ss,mesh(lev)%n_mortar
        sR  =  mesh(lev)%mortar_LR(3,im)
        eR  =  mesh(lev)%mortar_LR(4,im)
        if(.not. sec(sR)%is_ghost)  cycle

        isr =  sec(sR)%ID_recv(1,eR)
        iele=  sec(sR)%ID_recv(2,eR)
        if(iele .gt. p2p(isr)%n_ele_recv_face)  stop 'Error: set_p2p_ele_order fails.'
    end do
!   check.
!   ----------------------------------------------------------------------------

    do isr=mesh(lev)%p2p_1,mesh(lev)%p2p_0
        if(allocated(p2p(isr)%isend))   deallocate(p2p(isr)%isend)
        if(allocated(p2p(isr)%irecv))   deallocate(p2p(isr)%irecv)
    end do

    return
    end subroutine set_p2p_element_order
