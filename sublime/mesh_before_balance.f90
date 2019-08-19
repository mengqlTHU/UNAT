!-------------------------------------------------------------------------------
!   module used in the load balance.
!-------------------------------------------------------------------------------
    module var_load_balance
        use var_kind_def
        implicit none

        character(len=80):: zone_name
        character(len=80),allocatable:: boconame(:)
        integer(dpI):: n_bocos  =  0
!       bocotype, ele1, ele0, ID_sec_g
        integer(dpI),allocatable:: bocoinfo(:,:)

        integer(dpI):: n_base   =  0
        integer(dpI):: n_zone   =  0
        integer(dpI):: n_ele_g  =  0
        integer(dpI):: n_elei_g =  0
        integer(dpI):: n_seci_g =  0

        integer(dpI):: n_ele_L  =  0
        integer(dpI),allocatable:: iA_n2e(:),jA_n2e(:),ele_type(:),elei_prc(:,:)

        integer(dpI):: n_vtx_g  =  0
        integer(dpI),allocatable:: iA_e2n(:),jA_e2n(:),vtx_prc(:,:)

        integer(dpI):: nnz_n2e_b=  0
        integer(dpI),allocatable:: iA_n2e_b(:),jA_n2e_b(:)
        integer(dpI),allocatable:: bct(:,:)

        integer(dpI):: n_ele_ghost  =  0
        integer(dpI),allocatable:: iA_ghost(:),jA_ghost(:)
        integer(dpI),allocatable:: ele_type_ghost(:),ele_idx_ghost(:)

        integer(dpI):: n_vtx_i_have     =  0
        integer(dpI):: n_vtx_i_need     =  0
        integer(dpI):: n_ele_i_need     =  0
        integer(dpI):: n_ele_i_have     =  0
        integer(dpI):: nnz_e2n_i_need   =  0
        integer(dpI):: nnz_n2e_i_need   =  0
        integer(dpI):: nnz_n2e_i_have   =  0
        integer(dpI),allocatable:: vtx_i_need(:),vtx_i_have(:)
        integer(dpI),allocatable:: e2n_to_send(:),e2n_i_need(:)
        integer(dpI),allocatable:: ele_i_need(:),ele_i_have(:)
        integer(dpI),allocatable:: n2e_i_need(:),n2e_i_have(:)

        integer(dpI):: n_ele_per    =  0
        logical(dpL),allocatable:: is_corner_per(:)
        integer(dpI),allocatable:: ele_per(:,:),iA_per(:),jA_per(:),ele_type_per(:)

        integer(dpI):: edgecut,options(3)
        integer(dpI),allocatable:: vtxdist(:),vwgt(:),part(:),ele_prc(:)
        real   (kind=8),allocatable:: tpwgts(:)

        integer(dpI),allocatable:: iA_ss(:),jA_ss(:)
        integer(dpI),allocatable:: ele_ID_ss(:),ele_type_ss(:),per_ss(:,:)

        real   (dpR),allocatable:: xyz_parallel(:)
        integer(dpI):: n_xyz_i_need =  0
        integer(dpI):: n_xyz_i_have =  0
        integer(dpI),allocatable:: xyz_i_need(:),xyz_i_have(:)
        real   (dpR),allocatable:: xyz_need(:),xyz_have(:)
    end module var_load_balance
!-------------------------------------------------------------------------------
!   get the coordinates of vertex.
!-------------------------------------------------------------------------------
    subroutine get_xyz_i_need
    use var_kind_def
    use var_global, only: n_dim,err_mem
    use var_load_balance
    use var_parallel
    implicit none
    integer(dpI):: i,j,ivtx,ip

    prc_info=  0
    j       =  n_vtx_g/nprc
    do i=1,n_xyz_i_need
        ivtx=  xyz_i_need(i)
        if(ivtx .ge. vtx_prc(1,nprc-1)) then
            ip  =  nprc-1
        else
            ip  = (ivtx-1)/j
        end if
        prc_info(ip,1)  =  prc_info(ip,1)+1
    end do
    call mpi_alltoall(prc_info(0,1), 1, mpi_dpI, prc_info(0,3), 1, mpi_dpI, &
        &  mpi_comm_world, mpi_err)

    n_xyz_i_have=  prc_info(0,3)
    do ip=1,nprc-1
        n_xyz_i_have=  n_xyz_i_have+prc_info(ip,3)
    end do
    if(allocated(xyz_i_have)) then
        if(n_xyz_i_have .gt. size(xyz_i_have)) then
            deallocate(xyz_i_have)
            allocate(xyz_i_have(n_xyz_i_have), stat=err_mem)
        end if
    else
        allocate(xyz_i_have(max(1, n_xyz_i_have)), stat=err_mem)
    end if
    call collect_data_parallel(prc_info(0,1), prc_info(0,2), xyz_i_need, prc_info(0,3), &
        &  prc_info(0,4), xyz_i_have)

    if(allocated(xyz_need)) then
        if(n_dim*n_xyz_i_need .gt. size(xyz_need)) then
            deallocate(xyz_need)
            allocate(xyz_need(n_dim*n_xyz_i_need), stat=err_mem)
        end if
    else
        allocate(xyz_need(n_dim*max(1, n_xyz_i_need)), stat=err_mem)
    end if
    if(allocated(xyz_have)) then
        if(n_dim*n_xyz_i_have .gt. size(xyz_have)) then
            deallocate(xyz_have)
            allocate(xyz_have(n_dim*n_xyz_i_have), stat=err_mem)
        end if
    else
        allocate(xyz_have(n_dim*max(1, n_xyz_i_have)), stat=err_mem)
    end if
    do i=1,n_xyz_i_have
        j   =  xyz_i_have(i)
        if((j .lt. vtx_prc(1,myid)) .or. (j .gt. vtx_prc(2,myid))) &
            &  stop 'Error: fails to get xyz.'
        j   =  j-vtx_prc(1,myid)+1
        xyz_have(1+n_dim*(i-1):n_dim*i) =  xyz_parallel(1+n_dim*(j-1):n_dim*j)
    end do
    prc_info(0:nprc-1,1)=  prc_info(0:nprc-1,3)*n_dim
    call dcollect_data_parallel(prc_info(0,1), prc_info(0,2), xyz_have, prc_info(0,3), &
        &  prc_info(0,4), xyz_need)

    return
    end subroutine get_xyz_i_need
!-------------------------------------------------------------------------------
!   get the e2n information, parallel.
!-------------------------------------------------------------------------------
    subroutine get_e2n_parallel
    use var_kind_def
    use var_global, only: err_mem
    use var_load_balance
    use var_parallel
    implicit none
    integer(dpI):: i,j,k,iele,ivtx,ip
    integer(dpI),allocatable:: e2n(:,:),ibuf(:,:)

!   ----------------------------------------------------------------------------
!   record the element-to-node information.
    k   =  iA_n2e(n_ele_L+1)-1
    allocate(e2n(2,k), stat=err_mem)
    i   =  0
    do iele=1,n_ele_L
    do j=iA_n2e(iele),iA_n2e(iele+1)-1
        ivtx=  jA_n2e(j)
        i   =  i+1
        if(i .gt. k)    stop 'Error: e2n is of wrong size.'
        e2n(1:2,i)  = (/ivtx, iele+elei_prc(1,myid)-1/)
    end do
    end do
    call iqsortcols(.true., 1, k, 1, 2, e2n)
!   record the element-to-node information.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   memory allocation and collect e2n infomation.
    prc_info=  0
    do i=1,k
        if(e2n(1,i) .ge. vtx_prc(1,nprc-1)) then
            ip  =  nprc-1
        else
            ip  = (e2n(1,i)-1)/(n_vtx_g/nprc)
        end if
        if((e2n(1,i) .lt. vtx_prc(1,ip)) .or. (e2n(1,i) .gt. vtx_prc(2,ip))) &
            &  stop 'Error: fails to find the ip of node.'
        prc_info(ip,1)  =  prc_info(ip,1)+1
    end do
    call mpi_alltoall(prc_info(0,1), 1, mpi_dpI, prc_info(0,3), 1, mpi_dpI, &
        &  mpi_comm_world, mpi_err)
    i   =  prc_info(0,3)
    do ip=1,nprc-1
        i   =  i+prc_info(ip,3)
    end do
    allocate(ibuf(2,i))
    prc_info(0:nprc-1,1)=  prc_info(0:nprc-1,1)*2
    call collect_data_parallel(prc_info(0,1), prc_info(0,2), e2n, prc_info(0,3), &
        &  prc_info(0,4), ibuf)
!   memory allocation and collect e2n infomation.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   save the received information in the CSR format.
    if(allocated(e2n))  deallocate(e2n)
    ip  = (prc_info(nprc-1,4)+prc_info(nprc-1,3)-1)/2
    call iqsortcols(.true., 1, ip, 1, 2, ibuf)
    allocate(iA_e2n(vtx_prc(2,myid)-vtx_prc(1,myid)+2))
    allocate(jA_e2n(ip))
    i   =  1
    ivtx=  0
    do while(i .le. ip)
        k   =  i
        do j=i+1,ip
            if(ibuf(1,j) .eq. ibuf(1,i)) then
                k   =  j
            else
                exit
            end if
        end do
        call iqsortcols(.true., i, k, 2, 2, ibuf)
        ivtx=  ivtx+1
        if(ibuf(1,i)-vtx_prc(1,myid)+1 .ne. ivtx)   stop 'Error: fails to cal e2n.'
        iA_e2n(ivtx  )  =  i
        iA_e2n(ivtx+1)  =  k+1
        jA_e2n(i:k)     =  ibuf(2,i:k)

        i   =  k+1
    end do
!   save the received information in the CSR format.
!   ----------------------------------------------------------------------------

    if(allocated(ibuf)) deallocate(ibuf)

    return
    end subroutine get_e2n_parallel
!-------------------------------------------------------------------------------
!   get the e2n information that myid needs.
!-------------------------------------------------------------------------------
    subroutine get_e2n_i_need
    use var_kind_def
    use var_global
    use var_load_balance
    use var_parallel
    implicit none
    integer(dpI):: i,j,ivtx,ip,n

    prc_info=  0
    do i=1,n_vtx_i_need
        ivtx=  vtx_i_need(i)
        if(ivtx .gt. vtx_prc(1,nprc-1)) then
            ip  =  nprc-1
        else
            ip  = (ivtx-1)/(n_vtx_g/nprc)
        end if
        prc_info(ip,1)  =  prc_info(ip,1)+1
    end do
    call mpi_alltoall(prc_info(0,1), 1, mpi_dpI, prc_info(0,3), 1, mpi_dpI, &
        &  mpi_comm_world, mpi_err)
    n   =  prc_info(0,3)
    do ip=1,nprc-1
        n   =  n+prc_info(ip,3)
    end do
    if(allocated(vtx_i_have)) then
        if(n .gt. size(vtx_i_have)) then
            deallocate(vtx_i_have)
            allocate(vtx_i_have(n), stat=err_mem)
        end if
    else
        allocate(vtx_i_have(max(n,1)), stat=err_mem)
    end if
    call collect_data_parallel(prc_info(0,1), prc_info(0,2), vtx_i_need, prc_info(0,3), &
        &  prc_info(0,4), vtx_i_have)

!   collect the e2n information to be sent.
    n   =  0
    do ip=0,nprc-1
    do i=prc_info(ip,4),prc_info(ip,3)+prc_info(ip,4)-1
        ivtx=  vtx_i_have(i)
        if((ivtx .lt. vtx_prc(1,myid)) .or. (ivtx .gt. vtx_prc(2,myid))) &
            &  stop 'Error: fails to collect the e2n.'
        ivtx=  ivtx-vtx_prc(1,myid)+1
        n   =  n+2+iA_e2n(ivtx+1)-iA_e2n(ivtx)
    end do
    end do
    if(allocated(e2n_to_send)) then
        if(n .gt. size(e2n_to_send)) then
            deallocate(e2n_to_send)
            allocate(e2n_to_send(n), stat=err_mem)
        end if
    else
        allocate(e2n_to_send(max(n,1)), stat=err_mem)
    end if

    n   =  0
    prc_info(0:nprc-1,1)=  0
    do ip=0,nprc-1
        prc_info(ip,2)  =  n+1
        do i=prc_info(ip,4),prc_info(ip,3)+prc_info(ip,4)-1
            ivtx=  vtx_i_have(i)-vtx_prc(1,myid)+1
            n   =  n+1
            e2n_to_send(n)  =  vtx_i_have(i)
            n   =  n+1
            e2n_to_send(n)  =  iA_e2n(ivtx+1)-iA_e2n(ivtx)
            do j=iA_e2n(ivtx),iA_e2n(ivtx+1)-1
                n   =  n+1
                e2n_to_send(n)  =  jA_e2n(j)
            end do
        end do
        prc_info(ip,1)  =  n-prc_info(ip,2)+1
    end do

    call mpi_alltoall(prc_info(0,1), 1, mpi_dpI, prc_info(0,3), 1, mpi_dpI, &
        &  mpi_comm_world, mpi_err)
    n   =  prc_info(0,3)
    do ip=1,nprc-1
        n   =  n+prc_info(ip,3)
    end do
    if(allocated(e2n_i_need)) then
        if(n .gt. size(e2n_i_need)) then
            deallocate(e2n_i_need)
            allocate(e2n_i_need(n), stat=err_mem)
        end if
    else
        if(n .gt. 0)    allocate(e2n_i_need(n), stat=err_mem)
    end if
    nnz_e2n_i_need  =  n
    call collect_data_parallel(prc_info(0,1), prc_info(0,2), e2n_to_send, prc_info(0,3), &
        &  prc_info(0,4), e2n_i_need)

    return
    end subroutine get_e2n_i_need
!-------------------------------------------------------------------------------
!   get the n2e information that myid needs.
!-------------------------------------------------------------------------------
    subroutine get_n2e_i_need
    use var_kind_def
    use var_global
    use var_load_balance
    use var_parallel
    implicit none
    integer(dpI):: i,j,iele,ip,n

    prc_info=  0
    do i=1,n_ele_i_need
        if(ele_i_need(i) .ge. elei_prc(1,nprc-1)) then
            ip  =  nprc-1
        else
            ip  = (ele_i_need(i)-1)/(n_elei_g/nprc)
        end if
        prc_info(ip,1)  =  prc_info(ip,1)+1
    end do

!   record the element that myid has.
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

!   collect the n2e information to be sent out.
    n   =  0
    do i=1,prc_info(nprc-1,3)+prc_info(nprc-1,4)-1
        iele=  ele_i_have(i)
        if((iele .lt. elei_prc(1,myid)) .or. (iele .gt. elei_prc(2,myid))) &
            &  stop 'Error: fails to collect ele_i_have.'
        iele=  iele-elei_prc(1,myid)+1
        n   =  n+3+iA_n2e(iele+1)-iA_n2e(iele)
    end do
    if(allocated(n2e_i_have)) then
        if(n .gt. size(n2e_i_have)) then
            deallocate(n2e_i_have)
            allocate(n2e_i_have(n), stat=err_mem)
        end if
    else
        allocate(n2e_i_have(max(n,10000)), stat=err_mem)
    end if
    n   =  0
    do ip=0,nprc-1
        prc_info(ip,2)  =  n+1
        do i=prc_info(ip,4),prc_info(ip,3)+prc_info(ip,4)-1
            iele=  ele_i_have(i)-elei_prc(1,myid)+1
!           ID, ele_type, n_n2e, n2e.
            n   =  n+1
            n2e_i_have(n)   =  ele_i_have(i)
            n   =  n+1
            n2e_i_have(n)   =  ele_type(iele)
            n   =  n+1
            n2e_i_have(n)   =  iA_n2e(iele+1)-iA_n2e(iele)
            do j=iA_n2e(iele),iA_n2e(iele+1)-1
                n   =  n+1
                n2e_i_have(n)   =  jA_n2e(j)
            end do
        end do
        prc_info(ip,1)  =  n-prc_info(ip,2)+1
    end do

    call mpi_alltoall(prc_info(0,1), 1, mpi_dpI, prc_info(0,3), 1, mpi_dpI, &
        &  mpi_comm_world, mpi_err)
    n   =  prc_info(0,3)
    do ip=1,nprc-1
        n   =  n+prc_info(ip,3)
    end do
    nnz_n2e_i_need  =  n
    if(allocated(n2e_i_need)) then
        if(n .gt. size(n2e_i_need)) then
            deallocate(n2e_i_need)
            allocate(n2e_i_need(n), stat=err_mem)
        end if
    else
        allocate(n2e_i_need(max(n,1)), stat=err_mem)
    end if
    nnz_n2e_i_need  =  n
    call collect_data_parallel(prc_info(0,1), prc_info(0,2), n2e_i_have, prc_info(0,3), &
        &  prc_info(0,4), n2e_i_need)

    return
    end subroutine get_n2e_i_need
!-------------------------------------------------------------------------------
!   get the element information, ghost.
!-------------------------------------------------------------------------------
    subroutine get_ele_ghost(is_cal_vertex,is_balanced)
    use var_kind_def
    use var_cgns, only: npe_ele_1st
    use var_global, only: err_mem
    use var_load_balance
    use var_parallel
    use var_slv, only: is_FV,is_heat
    implicit none
    logical(dpL),intent(in):: is_cal_vertex,is_balanced
    logical(dpL):: ltmp
    integer(dpI):: i,j,n,ip,iele,ivtx,npe,L,R

    if(is_cal_vertex) then
!       ------------------------------------------------------------------------
!       memory allocation.
        n_vtx_i_need=  0
        if(is_FV .or. is_heat) then
            n_vtx_i_need=  iA_n2e(n_ele_L+1)-1
        else
            do iele=1,n_ele_L
                n_vtx_i_need=  n_vtx_i_need+npe_ele_1st(ele_type(iele))
            end do
        end if
        if(allocated(vtx_i_need)) then
            if(n_vtx_i_need .gt. size(vtx_i_need)) then
                deallocate(vtx_i_need)
                allocate(vtx_i_need(n_vtx_i_need), stat=err_mem)
            end if
        else
            allocate(vtx_i_need(n_vtx_i_need), stat=err_mem)
        end if
!       memory allocation.
!       ------------------------------------------------------------------------

!       ------------------------------------------------------------------------
!       record all vertex to be queried.
        if(is_FV .or. is_heat) then
            n_vtx_i_need=  iA_n2e(n_ele_L+1)-1
            call ICOPY(n_vtx_i_need, jA_n2e, 1, vtx_i_need, 1)
        else
            n_vtx_i_need=  0
            do iele=1,n_ele_L
                do j=iA_n2e(iele),iA_n2e(iele)+npe_ele_1st(ele_type(iele))-1
                    n_vtx_i_need            =  n_vtx_i_need+1
                    vtx_i_need(n_vtx_i_need)=  jA_n2e(j)
                end do
            end do
        end if
        call simplify_series(n_vtx_i_need, 1, 1, vtx_i_need)
!       record all vertex to be queried.
!       ------------------------------------------------------------------------
    end if

    call get_e2n_i_need
    call mpi_barrier(mpi_comm_world, mpi_err)

!   ----------------------------------------------------------------------------
!   record the element that myid needs.
    i   =  0
    n   =  0
    do while(i .lt. nnz_e2n_i_need)
        i   =  i+1
        ivtx=  e2n_i_need(i)
        i   =  i+1
        ip  =  e2n_i_need(i)
        do j=1,ip
            i   =  i+1
            iele=  e2n_i_need(i)
            if(is_balanced) then
                call ib_search(1, n_ele_i_have, 1, 1, ele_i_have, iele, ltmp, L, R)
            else
                ltmp= (iele .ge. elei_prc(1,myid)) .and. (iele .le. elei_prc(2,myid))
            end if
            if(ltmp)    cycle
            n   =  n+1
        end do
    end do

    if(allocated(ele_i_need)) then
        if(n .gt. size(ele_i_need)) then
            deallocate(ele_i_need)
            allocate(ele_i_need(n), stat=err_mem)
        end if
    else
        allocate(ele_i_need(max(n,1)), stat=err_mem)
    end if

    i   =  0
    n   =  0
    do while(i .lt. nnz_e2n_i_need)
        i   =  i+1
        ivtx=  e2n_i_need(i)
        i   =  i+1
        ip  =  e2n_i_need(i)
        do j=1,ip
            i   =  i+1
            iele=  e2n_i_need(i)

            ltmp=  .false.
            if(is_balanced) then
                call ib_search(1, n_ele_i_have, 1, 1, ele_i_have, iele, ltmp, L, R)
            else
                ltmp= (iele .ge. elei_prc(1,myid)) .and. (iele .le. elei_prc(2,myid))
            end if
            if(ltmp)    cycle

            n   =  n+1
            ele_i_need(n)   =  iele
        end do
    end do
    n_ele_i_need=  n
    if(n_ele_i_need .gt. 0) call simplify_series(n_ele_i_need, 1, 1, ele_i_need)
!   record the element that myid needs.
!   ----------------------------------------------------------------------------

    call get_n2e_i_need

    if(allocated(iA_ghost  ))   deallocate(iA_ghost  )
    if(allocated(jA_ghost  ))   deallocate(jA_ghost  )
    if(allocated(ele_type_ghost))   deallocate(ele_type_ghost)
    if(allocated(ele_idx_ghost ))   deallocate(ele_idx_ghost )
    n   =  0
    i   =  0
    ip  =  0
    do while(i .lt. nnz_n2e_i_need)
        n   =  n+1

!       ID, ele_type, n_n2e, n2e.
        i   =  i+3
        npe =  n2e_i_need(i)
        ip  =  ip+npe
        i   =  i +npe
    end do
    n_ele_ghost =  n
    allocate(iA_ghost(n+1))
    allocate(jA_ghost(ip))
    allocate(ele_type_ghost(n))
    allocate(ele_idx_ghost(n))
    n   =  0
    i   =  0
    j   =  0
    iA_ghost(1) =  1
    do while(i .lt. nnz_n2e_i_need)
        n   =  n+1

!       ID, ele_type, n_n2e, n2e.
        i   =  i+1
        ele_idx_ghost (n)   =  n2e_i_need(i)
        i   =  i+1
        ele_type_ghost(n)   =  n2e_i_need(i)
        i   =  i+1
        iA_ghost(n+1)   =  iA_ghost(n)+n2e_i_need(i)
        do ip=1,iA_ghost(n+1)-iA_ghost(n)
            i   =  i+1
            j   =  j+1
            jA_ghost(j) =  n2e_i_need(i)
        end do
    end do

    return
    end subroutine get_ele_ghost
!-------------------------------------------------------------------------------
!   get the element information, ghost, periodic.
!-------------------------------------------------------------------------------
    subroutine get_ele_per(is_cal_vertex,is_receiver_vertex)
    use var_kind_def
    use var_cgns, only: npe_ele
    use var_global, only: err_mem
    use var_load_balance
    use var_parallel
    use var_per_bnd
    implicit none
    logical(dpL),intent(in):: is_cal_vertex,is_receiver_vertex
    logical(dpL):: ltmp
    integer(dpI):: i,j,k,iele,ele,ivtx,L,R,npe,LL,RR,L2,R2
    integer(dpI),allocatable:: idx_v(:),iA_v(:),jA_v(:),idx_e(:),type_e(:),iA_e(:), &
                &  jA_e(:)

    if(n_vtx_pair_per .le. 0)   return
    call iqsortcols(.true., 1, n_vtx_pair_per, 2, 5, vtx_pair_per)

    if(is_cal_vertex) then
!       collect the donor vertex myid needs.
        if(allocated(vtx_i_need)) then
            if(n_vtx_pair_per .gt. size(vtx_i_need)) then
                deallocate(vtx_i_need)
                allocate(vtx_i_need(n_vtx_pair_per), stat=err_mem)
            end if
        else
            allocate(vtx_i_need(n_vtx_pair_per), stat=err_mem)
        end if
        vtx_i_need  =  0
        do iele=1,n_ele_L
            ele =  iele+elei_prc(1,myid)-1
            do j=iA_n2e(iele),iA_n2e(iele+1)-1
                ivtx=  jA_n2e(j)
                call ib_search(1, n_vtx_pair_per, 2, 5, vtx_pair_per, ivtx, ltmp, L, R)
                if(.not. ltmp)  cycle

                do i=L,R
                    vtx_i_need(i)   =  1
                end do
            end do
        end do
        n_vtx_i_need=  0
        do i=1,n_vtx_pair_per
            if(vtx_i_need(i) .ne. 1)    cycle
            n_vtx_i_need=  n_vtx_i_need+1
            vtx_i_need(n_vtx_i_need)=  vtx_pair_per(1,i)
        end do
        call simplify_series(n_vtx_i_need, 1, 1, vtx_i_need)
    end if

!   ----------------------------------------------------------------------------
!   collect the e2n information of the donor vertex.
    call get_e2n_i_need
    i   =  0
    j   =  0
    do while(i .lt. nnz_e2n_i_need)
        i   =  i+2
        j   =  j+e2n_i_need(i)
        i   =  i+e2n_i_need(i)
    end do
    allocate(idx_v(n_vtx_i_need))
    allocate(iA_v(n_vtx_i_need+1))
    allocate(jA_v(j))

    i   =  0
    j   =  0
    R   =  0
    iA_v(1) =  1
    do while(i .lt. nnz_e2n_i_need)
        R   =  R+1

        i       =  i+1
        idx_v(R)=  e2n_i_need(i)
        i       =  i+1
        ele     =  e2n_i_need(i)
        do L=1,ele
            i   =  i+1
            j   =  j+1
            jA_v(j) =  e2n_i_need(i)
        end do
        iA_v(R+1)   =  j+1
    end do
!   collect the e2n information of the donor vertex.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   collect the n2e information of the donor element.
    n_ele_i_need=  size(jA_v)
    if(allocated(ele_i_need)) then
        if(n_ele_i_need .gt. size(ele_i_need)) then
            deallocate(ele_i_need)
            allocate(ele_i_need(n_ele_i_need))
        end if
    else
        if(n_ele_i_need .gt. 0) allocate(ele_i_need(n_ele_i_need))
    end if
    call ICOPY(n_ele_i_need, jA_v, 1, ele_i_need, 1)
    call simplify_series(n_ele_i_need, 1, 1, ele_i_need)
    call get_n2e_i_need
!   collect the n2e information of the donor element.
!   ----------------------------------------------------------------------------

    allocate(idx_e(n_ele_i_need))
    allocate(type_e(n_ele_i_need))
    allocate(iA_e(n_ele_i_need+1))
    i   =  0
    j   =  0
    do while(i .lt. nnz_n2e_i_need)
!       ID, ele_type, n_n2e, n2e.
        i   =  i+3
        npe =  n2e_i_need(i)
        j   =  j+npe
        i   =  i+npe
    end do
    allocate(jA_e(j))

    i       =  0
    j       =  0
    iA_e(1) =  1
    L       =  0
    do while(i .lt. nnz_n2e_i_need)
        L   =  L+1

!       ID, ele_type, n_n2e, n2e.
        i   =  i+1
        idx_e(L)=  n2e_i_need(i)
        i   =  i+1
        type_e(L)   =  n2e_i_need(i)
        i   =  i+1
        npe =  n2e_i_need(i)
        do R=1,npe
            i   =  i+1
            j   =  j+1
            jA_e(j) =  n2e_i_need(i)
        end do
        iA_e(L+1)   =  j+1
    end do

!   ----------------------------------------------------------------------------
!   record the element and periodic path information.
    call iqsortcols(.true., 1, n_vtx_pair_per, 1, 5, vtx_pair_per)
    ele =  0
    do i=1,n_vtx_i_need
        ivtx=  vtx_i_need(i)
        call ib_search(1, n_vtx_pair_per, 1, 5, vtx_pair_per, ivtx, ltmp, LL, RR)
        if(.not. ltmp)  stop 'Error: fails to locate vtx_per_donor.'

        do j=LL,RR
            call ib_search(1, n_vtx_i_need, 1, 1, idx_v, ivtx, ltmp, L, R)
            if((.not. ltmp) .or. (L .ne. R))  stop 'Error: fails to get id_vtx_per_donor.'
            ele =  ele+iA_v(L+1)-iA_v(L)
        end do
    end do

    if(allocated(ele_per)) then
        if(size(ele_per, dim=2) .lt. ele) then
            deallocate(ele_per)
            allocate(ele_per(4,ele), stat=err_mem)
        end if
    else
        if(ele .gt. 0)  allocate(ele_per(4,ele), stat=err_mem)
    end if

    ele =  0
    do i=1,n_vtx_i_need
        ivtx=  vtx_i_need(i)
        call ib_search(1, n_vtx_pair_per, 1, 5, vtx_pair_per, ivtx, ltmp, LL, RR)

        do j=LL,RR
            call ib_search(1, n_vtx_i_need, 1, 1, idx_v, ivtx, ltmp, L, R)
            do R=iA_v(L),iA_v(L+1)-1
                ele =  ele+1
                ele_per(1:4,ele)= (/jA_v(R), vtx_pair_per(3:5,j)/)
            end do
        end do
    end do

    do i=1,ele
        iele=  ele_per(1,i)
        call ib_search(1, n_ele_i_need, 1, 1, idx_e, iele, ltmp, L, R)
        if((.not. ltmp) .or. (L .ne. R))    stop 'Error: fails to find doner ele, per.'
        ele_per(1,i)=  L
    end do
    if(ele .gt. 0)  call iqsortcols(.true., 1, ele, 1, 4, ele_per)
    i   =  1
    do while(i .le. ele)
        k   =  i
        do j=i+1,ele
            if(ele_per(1,j) .eq. ele_per(1,i)) then
                k   =  j
            else
                exit
            end if
        end do
        do j=i+1,k
            ltmp=  .false.
            do L=i,j-1
                if(ele_per(1,L) .lt. 0) cycle
                ltmp= (ele_per(2,L) .eq. ele_per(2,j)) .and. &
                    & (ele_per(3,L) .eq. ele_per(3,j)) .and. &
                    & (ele_per(4,L) .eq. ele_per(4,j))
                if(ltmp)    exit
            end do
            if(ltmp)    ele_per(1:4,j)  = -1
        end do

        i   =  k+1
    end do
    n_ele_per   =  0
    j           =  0
    do i=1,ele
        if(ele_per(1,i) .lt. 0) cycle
        n_ele_per   =  n_ele_per+1
        ele_per(1:4,n_ele_per)  =  ele_per(1:4,i)
        npe =  npe_ele(type_e(ele_per(1,n_ele_per)))
        j   =  j+npe
    end do
!   record the element and periodic path information.
!   ----------------------------------------------------------------------------

    if(allocated(iA_per)) then
        if(size(iA_per) .lt. n_ele_per+1) then
            deallocate(iA_per)
            allocate(iA_per(n_ele_per+1), stat=err_mem)
        end if
    else
        allocate(iA_per(n_ele_per+1), stat=err_mem)
    end if
    if(allocated(jA_per)) then
        if(size(jA_per) .lt. j) then
            deallocate(jA_per)
            allocate(jA_per(j), stat=err_mem)
        end if
    else
        if(j .gt. 0)    allocate(jA_per(j), stat=err_mem)
    end if
    if(allocated(ele_type_per)) then
        if(size(ele_type_per) .lt. n_ele_per) then
            deallocate(ele_type_per)
            allocate(ele_type_per(n_ele_per), stat=err_mem)
        end if
    else
        if(n_ele_per .gt. 0)    allocate(ele_type_per(n_ele_per), stat=err_mem)
    end if
    iA_per(1)   =  1
    j           =  0
    do i=1,n_ele_per
        k   =  ele_per(1,i)
        ele_type_per(i) =  type_e(k)
        ele_per   (1,i) =  idx_e (k)

        do L=iA_e(k),iA_e(k)+npe_ele(type_e(k))-1
            j   =  j+1
            jA_per(j)   = -jA_e(L)
            if(.not. is_receiver_vertex)    cycle

!           try to express the n2e with the final vertex ID.
            call ib_search(1, n_vtx_pair_per, 1, 5, vtx_pair_per, jA_e(L), ltmp, L2, R2)
            if(ltmp) then
                do k=L2,R2
                    if((vtx_pair_per(3,k) .eq. ele_per(2,i)) .and. &
                      &(vtx_pair_per(4,k) .eq. ele_per(3,i)) .and. &
                      &(vtx_pair_per(5,k) .eq. ele_per(4,i))) then
                        jA_per(j)   =  vtx_pair_per(2,k)
                        exit
                    end if
                end do
            end if
        end do
        iA_per(i+1) =  j+1
    end do

    if(allocated(idx_v ))   deallocate(idx_v)
    if(allocated(iA_v  ))   deallocate(iA_v)
    if(allocated(jA_v  ))   deallocate(jA_v)
    if(allocated(idx_e ))   deallocate(idx_e)
    if(allocated(type_e))   deallocate(type_e)
    if(allocated(iA_e  ))   deallocate(iA_e)
    if(allocated(jA_e  ))   deallocate(jA_e)

    return
    end subroutine get_ele_per
!-------------------------------------------------------------------------------
!   load balance.
!-------------------------------------------------------------------------------
    subroutine load_balance
    use var_kind_def
    use var_cgns
    use var_global, only: err_mem
    use var_load_balance
    use var_mesh, only: n_LR_hanging,LR_hanging
    use var_parallel
    implicit none
    integer(dpI):: nfac,iele,ele,i,j,k,v(8),nfac_i,n_matched,ip
    integer(dpI),allocatable:: fac(:,:),iA(:),jA(:)

!   ----------------------------------------------------------------------------
!   memory allocation.
    nfac=  0
    do iele=1,n_ele_L
        nfac=  nfac+nface_ele(ele_type(iele))
    end do
    do iele=1,n_ele_ghost
        nfac=  nfac+nface_ele(ele_type_ghost(iele))
    end do
    do iele=1,n_ele_per
        nfac=  nfac+nface_ele(ele_type_per(iele))
    end do

    do i=1,n_LR_hanging
        k   =  0
        do j=1,5
            if(LR_hanging(j,i) .le. 0) then
                exit
            else
                k   =  j
            end if
        end do
        ip  =  LR_hanging(1,i)
        if((ip .ge. elei_prc(1,myid)) .and. (ip .le. elei_prc(2,myid))) nfac=  nfac+k-1
        do j=2,k
            ip  =  LR_hanging(j,i)
            if((ip .ge. elei_prc(1,myid)) .and. (ip .le. elei_prc(2,myid))) nfac=  nfac+1
        end do
    end do

    allocate(fac(6,nfac), stat=err_mem)
    nfac=  0
    fac =  0
!   memory allocation.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   record the faces of elements of myid.
    do iele=1,n_ele_L
        ele =  ele_type(iele)
        i   =  iA_n2e(iele)

        if((ele .eq. TRI_3)) then
            v(1:3)  =  jA_n2e(i:i+2)
            fac(1:5,nfac+1) = (/v(1), v(2), 0, 0, iele+elei_prc(1,myid)-1/)
            fac(1:5,nfac+2) = (/v(2), v(3), 0, 0, iele+elei_prc(1,myid)-1/)
            fac(1:5,nfac+3) = (/v(3), v(1), 0, 0, iele+elei_prc(1,myid)-1/)
            nfac    =  nfac+3
        elseif((ele .eq. QUAD_4) .or. (ele .eq. QUAD_8) .or. (ele .eq. QUAD_9)) then
            v(1:4)  =  jA_n2e(i:i+3)
            fac(1:5,nfac+1) = (/v(1), v(2), 0, 0, iele+elei_prc(1,myid)-1/)
            fac(1:5,nfac+2) = (/v(2), v(3), 0, 0, iele+elei_prc(1,myid)-1/)
            fac(1:5,nfac+3) = (/v(3), v(4), 0, 0, iele+elei_prc(1,myid)-1/)
            fac(1:5,nfac+4) = (/v(4), v(1), 0, 0, iele+elei_prc(1,myid)-1/)
            nfac    =  nfac+4
        elseif(ele .eq. TETRA_4) then
            v(1:4)  =  jA_n2e(i:i+3)
            fac(1:5,nfac+1) = (/v(1), v(2), v(3), 0, iele+elei_prc(1,myid)-1/)
            fac(1:5,nfac+2) = (/v(1), v(2), v(4), 0, iele+elei_prc(1,myid)-1/)
            fac(1:5,nfac+3) = (/v(2), v(3), v(4), 0, iele+elei_prc(1,myid)-1/)
            fac(1:5,nfac+4) = (/v(1), v(3), v(4), 0, iele+elei_prc(1,myid)-1/)
            nfac    =  nfac+4
        elseif(ele .eq. PYRA_5) then
            v(1:5)  =  jA_n2e(i:i+4)
            fac(1:5,nfac+1) = (/v(1), v(2), v(3), v(4), iele+elei_prc(1,myid)-1/)
            fac(1:5,nfac+2) = (/v(1), v(2), v(5), 0, iele+elei_prc(1,myid)-1/)
            fac(1:5,nfac+3) = (/v(2), v(3), v(5), 0, iele+elei_prc(1,myid)-1/)
            fac(1:5,nfac+4) = (/v(3), v(4), v(5), 0, iele+elei_prc(1,myid)-1/)
            fac(1:5,nfac+5) = (/v(4), v(1), v(5), 0, iele+elei_prc(1,myid)-1/)
            nfac    =  nfac+5
        elseif(ele .eq. PENTA_6) then
            v(1:6)  =  jA_n2e(i:i+5)
            fac(1:5,nfac+1) = (/v(1), v(2), v(3), 0, iele+elei_prc(1,myid)-1/)
            fac(1:5,nfac+2) = (/v(4), v(5), v(6), 0, iele+elei_prc(1,myid)-1/)
            fac(1:5,nfac+3) = (/v(1), v(2), v(5), v(4), iele+elei_prc(1,myid)-1/)
            fac(1:5,nfac+4) = (/v(2), v(3), v(6), v(5), iele+elei_prc(1,myid)-1/)
            fac(1:5,nfac+5) = (/v(1), v(3), v(6), v(4), iele+elei_prc(1,myid)-1/)
            nfac    =  nfac+5
        elseif((ele .eq. HEXA_8) .or. (ele .eq. HEXA_20) .or. (ele .eq. HEXA_27)) then
            v(1:8)  =  jA_n2e(i:i+7)
            fac(1:5,nfac+1) = (/v(1), v(2), v(3), v(4), iele+elei_prc(1,myid)-1/)
            fac(1:5,nfac+2) = (/v(5), v(6), v(7), v(8), iele+elei_prc(1,myid)-1/)
            fac(1:5,nfac+3) = (/v(1), v(5), v(8), v(4), iele+elei_prc(1,myid)-1/)
            fac(1:5,nfac+4) = (/v(2), v(6), v(7), v(3), iele+elei_prc(1,myid)-1/)
            fac(1:5,nfac+5) = (/v(1), v(2), v(6), v(5), iele+elei_prc(1,myid)-1/)
            fac(1:5,nfac+6) = (/v(4), v(3), v(7), v(8), iele+elei_prc(1,myid)-1/)
            nfac    =  nfac+6
        else
            stop 'Error: element type not supported.'
        end if
    end do
    nfac_i  =  nfac
    do i=1,nfac
        call iqsort(.true., 1, 4, fac(1,i))
    end do
    call iqsortcols(.true., 1, nfac_i, 4, 6, fac)
!   record the faces of elements of myid.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   record the faces of ghost elements.
    do iele=1,n_ele_ghost
        ele =  ele_type_ghost(iele)
        i   =  iA_ghost(iele)

        if((ele .eq. TRI_3)) then
            v(1:3)  =  jA_ghost(i:i+2)
            fac(1:5,nfac+1) = (/v(1), v(2), 0, 0, ele_idx_ghost(iele)/)
            fac(1:5,nfac+2) = (/v(2), v(3), 0, 0, ele_idx_ghost(iele)/)
            fac(1:5,nfac+3) = (/v(3), v(1), 0, 0, ele_idx_ghost(iele)/)
            nfac    =  nfac+3
        elseif((ele .eq. QUAD_4) .or. (ele .eq. QUAD_8) .or. (ele .eq. QUAD_9)) then
            v(1:4)  =  jA_ghost(i:i+3)
            fac(1:5,nfac+1) = (/v(1), v(2), 0, 0, ele_idx_ghost(iele)/)
            fac(1:5,nfac+2) = (/v(2), v(3), 0, 0, ele_idx_ghost(iele)/)
            fac(1:5,nfac+3) = (/v(3), v(4), 0, 0, ele_idx_ghost(iele)/)
            fac(1:5,nfac+4) = (/v(4), v(1), 0, 0, ele_idx_ghost(iele)/)
            nfac    =  nfac+4
        elseif(ele .eq. TETRA_4) then
            v(1:4)  =  jA_ghost(i:i+3)
            fac(1:5,nfac+1) = (/v(1), v(2), v(3), 0, ele_idx_ghost(iele)/)
            fac(1:5,nfac+2) = (/v(1), v(2), v(4), 0, ele_idx_ghost(iele)/)
            fac(1:5,nfac+3) = (/v(2), v(3), v(4), 0, ele_idx_ghost(iele)/)
            fac(1:5,nfac+4) = (/v(1), v(3), v(4), 0, ele_idx_ghost(iele)/)
            nfac    =  nfac+4
        elseif(ele .eq. PYRA_5) then
            v(1:5)  =  jA_ghost(i:i+4)
            fac(1:5,nfac+1) = (/v(1), v(2), v(3), v(4), ele_idx_ghost(iele)/)
            fac(1:5,nfac+2) = (/v(1), v(2), v(5), 0, ele_idx_ghost(iele)/)
            fac(1:5,nfac+3) = (/v(2), v(3), v(5), 0, ele_idx_ghost(iele)/)
            fac(1:5,nfac+4) = (/v(3), v(4), v(5), 0, ele_idx_ghost(iele)/)
            fac(1:5,nfac+5) = (/v(4), v(1), v(5), 0, ele_idx_ghost(iele)/)
            nfac    =  nfac+5
        elseif(ele .eq. PENTA_6) then
            v(1:6)  =  jA_ghost(i:i+5)
            fac(1:5,nfac+1) = (/v(1), v(2), v(3), 0, ele_idx_ghost(iele)/)
            fac(1:5,nfac+2) = (/v(4), v(5), v(6), 0, ele_idx_ghost(iele)/)
            fac(1:5,nfac+3) = (/v(1), v(2), v(5), v(4), ele_idx_ghost(iele)/)
            fac(1:5,nfac+4) = (/v(2), v(3), v(6), v(5), ele_idx_ghost(iele)/)
            fac(1:5,nfac+5) = (/v(1), v(3), v(6), v(4), ele_idx_ghost(iele)/)
            nfac    =  nfac+5
        elseif((ele .eq. HEXA_8) .or. (ele .eq. HEXA_20) .or. (ele .eq. HEXA_27)) then
            v(1:8)  =  jA_ghost(i:i+7)
            fac(1:5,nfac+1) = (/v(1), v(2), v(3), v(4), ele_idx_ghost(iele)/)
            fac(1:5,nfac+2) = (/v(5), v(6), v(7), v(8), ele_idx_ghost(iele)/)
            fac(1:5,nfac+3) = (/v(1), v(5), v(8), v(4), ele_idx_ghost(iele)/)
            fac(1:5,nfac+4) = (/v(2), v(6), v(7), v(3), ele_idx_ghost(iele)/)
            fac(1:5,nfac+5) = (/v(1), v(2), v(6), v(5), ele_idx_ghost(iele)/)
            fac(1:5,nfac+6) = (/v(4), v(3), v(7), v(8), ele_idx_ghost(iele)/)
            nfac    =  nfac+6
        else
            stop 'Error: element type not supported.'
        end if
    end do
!   record the faces of ghost elements.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   record the faces of ghost_per elements.
    do iele=1,n_ele_per
        ele =  ele_type_per(iele)
        i   =  iA_per(iele)

        if((ele .eq. TRI_3)) then
            v(1:3)  =  jA_per(i:i+2)
            fac(1:5,nfac+1) = (/v(1), v(2), 0, 0, ele_per(1,iele)/)
            fac(1:5,nfac+2) = (/v(2), v(3), 0, 0, ele_per(1,iele)/)
            fac(1:5,nfac+3) = (/v(3), v(1), 0, 0, ele_per(1,iele)/)
            nfac    =  nfac+3
        elseif((ele .eq. QUAD_4) .or. (ele .eq. QUAD_8) .or. (ele .eq. QUAD_9)) then
            v(1:4)  =  jA_per(i:i+3)
            fac(1:5,nfac+1) = (/v(1), v(2), 0, 0, ele_per(1,iele)/)
            fac(1:5,nfac+2) = (/v(2), v(3), 0, 0, ele_per(1,iele)/)
            fac(1:5,nfac+3) = (/v(3), v(4), 0, 0, ele_per(1,iele)/)
            fac(1:5,nfac+4) = (/v(4), v(1), 0, 0, ele_per(1,iele)/)
            nfac    =  nfac+4
        elseif(ele .eq. TETRA_4) then
            v(1:4)  =  jA_per(i:i+3)
            fac(1:5,nfac+1) = (/v(1), v(2), v(3), 0, ele_per(1,iele)/)
            fac(1:5,nfac+2) = (/v(1), v(2), v(4), 0, ele_per(1,iele)/)
            fac(1:5,nfac+3) = (/v(2), v(3), v(4), 0, ele_per(1,iele)/)
            fac(1:5,nfac+4) = (/v(1), v(3), v(4), 0, ele_per(1,iele)/)
            nfac    =  nfac+4
        elseif(ele .eq. PYRA_5) then
            v(1:5)  =  jA_per(i:i+4)
            fac(1:5,nfac+1) = (/v(1), v(2), v(3), v(4), ele_per(1,iele)/)
            fac(1:5,nfac+2) = (/v(1), v(2), v(5), 0, ele_per(1,iele)/)
            fac(1:5,nfac+3) = (/v(2), v(3), v(5), 0, ele_per(1,iele)/)
            fac(1:5,nfac+4) = (/v(3), v(4), v(5), 0, ele_per(1,iele)/)
            fac(1:5,nfac+5) = (/v(4), v(1), v(5), 0, ele_per(1,iele)/)
            nfac    =  nfac+5
        elseif(ele .eq. PENTA_6) then
            v(1:6)  =  jA_per(i:i+5)
            fac(1:5,nfac+1) = (/v(1), v(2), v(3), 0, ele_per(1,iele)/)
            fac(1:5,nfac+2) = (/v(4), v(5), v(6), 0, ele_per(1,iele)/)
            fac(1:5,nfac+3) = (/v(1), v(2), v(5), v(4), ele_per(1,iele)/)
            fac(1:5,nfac+4) = (/v(2), v(3), v(6), v(5), ele_per(1,iele)/)
            fac(1:5,nfac+5) = (/v(1), v(3), v(6), v(4), ele_per(1,iele)/)
            nfac    =  nfac+5
        elseif((ele .eq. HEXA_8) .or. (ele .eq. HEXA_20) .or. (ele .eq. HEXA_27)) then
            v(1:8)  =  jA_per(i:i+7)
            fac(1:5,nfac+1) = (/v(1), v(2), v(3), v(4), ele_per(1,iele)/)
            fac(1:5,nfac+2) = (/v(5), v(6), v(7), v(8), ele_per(1,iele)/)
            fac(1:5,nfac+3) = (/v(1), v(5), v(8), v(4), ele_per(1,iele)/)
            fac(1:5,nfac+4) = (/v(2), v(6), v(7), v(3), ele_per(1,iele)/)
            fac(1:5,nfac+5) = (/v(1), v(2), v(6), v(5), ele_per(1,iele)/)
            fac(1:5,nfac+6) = (/v(4), v(3), v(7), v(8), ele_per(1,iele)/)
            nfac    =  nfac+6
        else
            stop 'Error: element type not supported.'
        end if
    end do
    do i=nfac_i+1,nfac
        call iqsort(.true., 1, 4, fac(1,i))
    end do
    call iqsortcols(.true., nfac_i+1, nfac, 4, 6, fac)
!   record the faces of ghost_per elements.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   find the element-to-element connectivity, hanging-face is not considered yet.
    n_matched   =  0
    do i=1,nfac_i
        if(fac(6,i) .gt. 0) cycle

!       matching with internal element.
        call find_face(nfac_i, fac, i, fac(1,i), j, v)
        if(j .gt. 0) then
            if(j .ne. 1)    stop 'Error: more than 2 faces share the same vertex.'
            n_matched   =  n_matched+1
            fac(6,i)    =  fac(5,v(1))
            fac(6,v(1)) =  fac(5,i)
            cycle
        end if

        if(nfac_i .eq. nfac)    cycle
!       matching with ghost element.
        call find_face(nfac-nfac_i, fac(1,nfac_i+1), 0, fac(1,i), j, v)
        if(j .gt. 0) then
            if(j .ne. 1)    stop 'Error: more than 2 faces share the same vertex.'
            n_matched           =  n_matched+1
            fac(6,i)            =  fac(5,v(1)+nfac_i)
            fac(6,v(1)+nfac_i)  =  fac(5,i)
        end if
    end do

    n_matched   =  0
    do i=1,nfac
        if(fac(6,i) .le. 0) cycle
        iele=  fac(5,i)
        ele =  fac(6,i)
        if((iele .lt. elei_prc(1,myid)) .or. (iele .gt. elei_prc(2,myid)))  cycle
        n_matched   =  n_matched+1
        fac(5:6,n_matched)  =  fac(5:6,i)
    end do
!   find the element-to-element connectivity, hanging-face is not considered yet.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   record the element-to-element connectivity due to hanging-face.
    do i=1,n_LR_hanging
        k   =  0
        do j=1,5
            if(LR_hanging(j,i) .le. 0) then
                exit
            else
                k   =  j
            end if
        end do
        ip  =  LR_hanging(1,i)
        if((ip .ge. elei_prc(1,myid)) .and. (ip .le. elei_prc(2,myid))) then
            do j=2,k
                n_matched       =  n_matched+1
                fac(5,n_matched)=  ip
                fac(6,n_matched)=  LR_hanging(j,i)
            end do
        end if

        do j=2,k
            ip  =  LR_hanging(j,i)
            if((ip .ge. elei_prc(1,myid)) .and. (ip .le. elei_prc(2,myid))) then
                n_matched       =  n_matched+1
                fac(5,n_matched)=  ip
                fac(6,n_matched)=  LR_hanging(1,i)
            end if
        end do
    end do
!   record the element-to-element connectivity due to hanging-face.
!   ----------------------------------------------------------------------------

    call iqsortcols(.true., 1, n_matched, 5, 6, fac)

    allocate(iA(n_ele_L+1))
    allocate(jA(n_matched))
    allocate(vwgt(n_ele_L))
    allocate(part(n_ele_L))
    i       =  1
    ele     =  0
    iA(1)   =  1
    do while(i .le. n_matched)
        k   =  i
        do j=i+1,n_matched
            if(fac(5,j) .ne. fac(5,i)) then
                exit
            else
                k   =  j
            end if
        end do
        ele =  ele+1
        if(fac(5,i) .ne. ele+elei_prc(1,myid)-1)    stop 'Error: fails to get iA.'
        iele=  k-i+1
        call simplify_series(iele, 6, 6, fac(1,i))
        iA(ele+1)   =  iA(ele)+iele
        vwgt(ele)   =  iele
        do j=1,iele
            jA(iA(ele)+j-1) =  fac(6,i+j-1)
        end do
        i           =  k+1
    end do

!   ----------------------------------------------------------------------------
!   prc(ip) stores ele(ip):ele(ip+1)-1.
    allocate(ele_prc(0:nprc))
    ele_prc(0)  =  1
    do ip=1,nprc
        ele_prc(ip) =  elei_prc(2,ip-1)+1
    end do
!   prc(ip) stores ele(ip):ele(ip+1)-1.
!   ----------------------------------------------------------------------------

    allocate(tpwgts(nprc), stat=err_mem)
    tpwgts  =  real(1, kind=8)/real(nprc, kind=8)
    if(nprc .gt. 1) then
        options =  0
        call parmetis_v3_partkway(ele_prc, iA, jA, vwgt, 0, 2, 1, &
            &  1, nprc, tpwgts, 1.02d0, options, edgecut, part, mpi_comm_world)
    else
        part=  1
    end if
    part=  part-1
    do i=1,n_ele_L
        if((part(i) .lt. 0) .or. (part(i) .ge. nprc))   stop 'Error: wrong ParMETIS out.'
    end do

    if(allocated(iA     ))  deallocate(iA     )
    if(allocated(jA     ))  deallocate(jA     )
    if(allocated(vtxdist))  deallocate(vtxdist)
    if(allocated(vwgt   ))  deallocate(vwgt   )
    if(allocated(ele_prc))  deallocate(ele_prc)
    if(allocated(fac    ))  deallocate(fac    )

    return
    contains
!       ------------------------------------------------------------------------
!       find matched face.
!       ------------------------------------------------------------------------
        subroutine find_face(nfac,fac,idx_ff,ff,nele,ele)
        implicit none
        integer(dpI),intent(in):: nfac,fac(6,*),idx_ff,ff(*)
        logical(dpL):: ltmp
        integer(dpI):: nele,ele(*),f(5),L,R,i,j

        nele    =  0
        f(1:4)  =  ff(1:4)
        call ib_search(1, nfac, 4, 6, fac, f(4), ltmp, L, R)
        if(.not. ltmp)  return

        do i=L,R
            if(i .eq. idx_ff)   cycle
            ltmp=  .false.
            do j=1,3
                ltmp=  fac(j,i) .ne. f(j)
                if(ltmp)    exit
            end do
            if(ltmp)    cycle
            nele=  nele+1
            ele(nele)   =  i
        end do

        return
        end subroutine find_face
    end subroutine load_balance
