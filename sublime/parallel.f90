!-------------------------------------------------------------------------------
!   module mpi setup.
!-------------------------------------------------------------------------------
    subroutine set_parallel(start_mpi)
    use var_kind_def
    use var_parallel
    implicit none
    logical(dpL),intent(in):: start_mpi

!   ----------------------------------------------------------------------------
!   mpi initialization.
    if(start_mpi)   call mpi_init(mpi_err)
    call mpi_comm_size(mpi_comm_world, nprc, mpi_err)
    call mpi_comm_rank(mpi_comm_world, myid, mpi_err)
    call mpi_get_processor_name(mpi_prc_name, mpi_prc_name_len, mpi_err)
    call mpi_comm_group(mpi_comm_world, mpi_group_world, mpi_err)
!   write(unit=6,fmt='(a8,i8,a3,i8,a4,a12,a9)'),'process',myid,'of',nprc,'on ', &
!       & mpi_prc_name,'started.'
    call mpi_barrier(mpi_comm_world, mpi_err)
    if(myid .eq. 0) then
        print*,'----------------------------------------------------------'
        print*,'Note: MPI initialization finished.'
    end if
!   mpi initialization.
!   ----------------------------------------------------------------------------

    allocate(mpi_req(                2*nprc))
    allocate(mpi_sta(mpi_status_size,2*nprc))
    allocate(prc_info(0:nprc-1,4))
    prc_info=  0
    if(dpL .eq. 1) then
        mpi_dpL =  mpi_logical
    elseif(dpL .eq. 2) then
        mpi_dpL =  mpi_logical
    elseif(dpL .eq. 4) then
        mpi_dpL =  mpi_logical
    elseif(dpL .eq. 8) then
        mpi_dpL =  mpi_logical
    else
        stop 'Error: logical type not supported by MPI.'
    end if

    if(dpI .eq. 1) then
        mpi_dpI =  mpi_integer1
    elseif(dpI .eq. 2) then
        mpi_dpI =  mpi_integer2
    elseif(dpI .eq. 4) then
        mpi_dpI =  mpi_integer4
    elseif(dpI .eq. 8) then
        mpi_dpI =  mpi_integer8
    elseif(dpI .eq. 16) then
        mpi_dpI =  mpi_integer16
    else
        stop 'Error: integer type not supported by MPI.'
    end if

    if(dpR .eq. 4) then
        mpi_dpR =  mpi_real4
    elseif(dpR .eq. 8) then
        mpi_dpR =  mpi_real8
    elseif(dpR .eq. 16) then
        mpi_dpR =  mpi_real16
    else
        stop 'Error: real type not supported by MPI.'
    end if

    return
    end subroutine set_parallel
!-------------------------------------------------------------------------------
!   collect data distributely.
!-------------------------------------------------------------------------------
    subroutine collect_data_parallel(count_S,ind_S,S,count_R,ind_R,R)
    use var_kind_def
    use var_parallel
    implicit none
    integer(dpI):: count_S(*),ind_S(*),S(*),count_R(*),ind_R(*),R(*),i,j,ip

    ind_S  (1:nprc) =  0
    count_R(1:nprc) =  0
    ind_R  (1:nprc) =  0
    call mpi_alltoall(count_S, 1, mpi_dpI, count_R, 1, mpi_dpI, mpi_comm_world, mpi_err)

    i   =  1
    do ip=0,nprc-1
        ind_S(ip+1) =  i
        i   =  i+count_S(ip+1)
    end do
    i   =  1
    do ip=0,nprc-1
        ind_R(ip+1) =  i
        i   =  i+count_R(ip+1)
    end do

    mpi_nreq=  0
    do ip=0,nprc-1
        i   =  ind_S(ip+1)
        j   =  ind_R(ip+1)

        if(myid .eq. ip) then
            call ICOPY(count_S(ip+1), S(i), 1, R(j), 1)
            cycle
        end if

        if(count_S(ip+1) .gt. 0) then
            mpi_nreq=  mpi_nreq+1
            call mpi_isend(S(i), count_S(ip+1), mpi_dpI, ip, myid, &
                &  mpi_comm_world, mpi_req(mpi_nreq), mpi_err)
        end if

        if(count_R(ip+1) .gt. 0) then
            mpi_nreq=  mpi_nreq+1
            call mpi_irecv(R(j), count_R(ip+1), mpi_dpI, ip, ip  , &
                &  mpi_comm_world, mpi_req(mpi_nreq), mpi_err)
        end if
    end do
    call mpi_waitall(mpi_nreq, mpi_req, mpi_sta, mpi_err)

    return
    end subroutine collect_data_parallel
!-------------------------------------------------------------------------------
!   collect data distributely, double precision.
!-------------------------------------------------------------------------------
    subroutine dcollect_data_parallel(count_S,ind_S,S,count_R,ind_R,R)
    use var_kind_def
    use var_parallel
    implicit none
    integer(dpI):: count_S(*),ind_S(*),count_R(*),ind_R(*),i,j,ip
    real   (dpR):: S(*),R(*)

    ind_S  (1:nprc) =  0
    count_R(1:nprc) =  0
    ind_R  (1:nprc) =  0
    call mpi_alltoall(count_S, 1, mpi_dpI, count_R, 1, mpi_dpI, mpi_comm_world, mpi_err)

    i   =  1
    do ip=0,nprc-1
        ind_S(ip+1) =  i
        i   =  i+count_S(ip+1)
    end do
    i   =  1
    do ip=0,nprc-1
        ind_R(ip+1) =  i
        i   =  i+count_R(ip+1)
    end do

    mpi_nreq=  0
    do ip=0,nprc-1
        i   =  ind_S(ip+1)
        j   =  ind_R(ip+1)

        if(myid .eq. ip) then
            call DCOPY(count_S(ip+1), S(i), 1, R(j), 1)
            cycle
        end if

        if(count_S(ip+1) .gt. 0) then
            mpi_nreq=  mpi_nreq+1
            call mpi_isend(S(i), count_S(ip+1), mpi_dpR, ip, myid, &
                &  mpi_comm_world, mpi_req(mpi_nreq), mpi_err)
        end if

        if(count_R(ip+1) .gt. 0) then
            mpi_nreq=  mpi_nreq+1
            call mpi_irecv(R(j), count_R(ip+1), mpi_dpR, ip, ip  , &
                &  mpi_comm_world, mpi_req(mpi_nreq), mpi_err)
        end if
    end do
    call mpi_waitall(mpi_nreq, mpi_req, mpi_sta, mpi_err)

    return
    end subroutine dcollect_data_parallel
!-------------------------------------------------------------------------------
!   allgather info from all processes.
!-------------------------------------------------------------------------------
    subroutine i_allgather(comm,nprc,N,b1,b2)
    use var_kind_def
    use var_parallel, only: mpi_dpI,mpi_err,prc_info
    implicit none
    integer(dpI):: comm,nprc,N,i,b1(*),b2(*)

    call mpi_allgather(N, 1, mpi_dpI, prc_info(0,2), 1, mpi_dpI, comm, mpi_err)
    prc_info(:,3)   =  0
    do i=1,nprc-1
        prc_info(i,3)   =  prc_info(i-1,3)+prc_info(i-1,2)
    end do
    call mpi_allgatherv(b1, N, mpi_dpI, b2, prc_info(0,2), &
        &  prc_info(0,3), mpi_dpI, comm, mpi_err)
    N   =  prc_info(nprc-1,3)+prc_info(nprc-1,2)

    return
    end subroutine i_allgather
!-------------------------------------------------------------------------------
!   allgather info from all processes.
!-------------------------------------------------------------------------------
    subroutine r_allgather(comm,nprc,N,b1,b2)
    use var_kind_def
    use var_parallel, only: mpi_dpI,mpi_dpR,mpi_err,prc_info
    implicit none
    integer(dpI),intent(in):: comm,nprc
    real   (dpR),intent(in):: b1(*)
    integer(dpI):: N,i
    real   (dpR):: b2(*)

    call mpi_allgather(N, 1, mpi_dpI, prc_info(0,2), 1, mpi_dpI, comm, mpi_err)
    prc_info(:,3)   =  0
    do i=1,nprc-1
        prc_info(i,3)   =  prc_info(i-1,3)+prc_info(i-1,2)
    end do
    call mpi_allgatherv(b1, N, mpi_dpR, b2, prc_info(0,2), prc_info(0,3), &
        &  mpi_dpR, comm, mpi_err)
    N   =  prc_info(nprc-1,3)+prc_info(nprc-1,2)

    return
    end subroutine r_allgather
!-------------------------------------------------------------------------------
!   gather info from all processes.
!-------------------------------------------------------------------------------
    subroutine r_gather(comm,nprc,N,b1,b2)
    use var_kind_def
    use var_parallel, only: mpi_dpI,mpi_dpR,mpi_err
    implicit none
    integer(dpI):: comm,nprc,N,i,prc_info(0:nprc-1,2:3)
    real   (dpR):: b1(*),b2(*)

    call mpi_gather(N, 1, mpi_dpI, prc_info(0,2), 1, mpi_dpI, 0, comm, mpi_err)
    prc_info(0,3)   =  0
    do i=1,nprc-1
        prc_info(i,3)   =  prc_info(i-1,3)+prc_info(i-1,2)
    end do
    call mpi_gatherv(b1, N, mpi_dpR, b2, prc_info(0,2), prc_info(0,3), mpi_dpR, &
        &  0, comm, mpi_err)
    N   =  prc_info(nprc-1,3)+prc_info(nprc-1,2)

    return
    end subroutine r_gather
!-------------------------------------------------------------------------------
!   gather info from all processes, integer.
!-------------------------------------------------------------------------------
    subroutine i_gather(comm,nprc,N,b1,b2)
    use var_kind_def
    use var_parallel, only: mpi_dpI,mpi_err,prc_info
    implicit none
    integer(dpI),intent(in):: comm,nprc
    integer(dpI):: N,i,b1(*),b2(*)

    call mpi_gather(N, 1, mpi_dpI, prc_info(0,2), 1, mpi_dpI, 0, comm, mpi_err)
    prc_info(0,3)   =  0
    do i=1,nprc-1
        prc_info(i,3)   =  prc_info(i-1,3)+prc_info(i-1,2)
    end do
    call mpi_gatherv(b1, N, mpi_dpI, b2, prc_info(0,2), prc_info(0,3), mpi_dpI, &
        &  0, comm, mpi_err)
    N   =  prc_info(nprc-1,3)+prc_info(nprc-1,2)

    return
    end subroutine i_gather
!-------------------------------------------------------------------------------
!   compute the norm of vector, parallel.
!-------------------------------------------------------------------------------
    function DNRM_mpi(N,v) result(nrm)
    use var_kind_def
    use var_parallel
    implicit none
    integer(dpI),intent(in):: N
    real   (dpR),intent(in):: v(*)
    real   (dpR):: nrm,rtmp
    real   (kind=8),external:: DNRM2

    rtmp=  DNRM2(N, v, 1)
    if(nprc .le. 1) then
        nrm =  rtmp
        return
    end if
    call mpi_allreduce(rtmp*rtmp, nrm, 1, mpi_dpR, mpi_sum, mpi_comm_world, mpi_err)
    nrm =  sqrt(nrm)

    return
    end function DNRM_mpi
!-------------------------------------------------------------------------------
!   compute the dot product of two vectors, parallel.
!-------------------------------------------------------------------------------
    function DDOT_mpi(N,v1,v2) result(dot)
    use var_kind_def
    use var_parallel
    implicit none
    integer(dpI),intent(in):: N
    real   (dpR),intent(in):: v1(*),v2(*)
    real   (dpR):: dot,rtmp
    real   (kind=8),external:: DDOT

    rtmp=  DDOT(N, v1, 1, v2, 1)
    if(nprc .le. 1) then
        dot =  rtmp
        return
    end if
    call mpi_allreduce(rtmp, dot, 1, mpi_dpR, mpi_sum, mpi_comm_world, mpi_err)

    return
    end function DDOT_mpi
!-------------------------------------------------------------------------------
!   compute the root-mean-square of a vector, parallel.
!-------------------------------------------------------------------------------
    subroutine get_rms(N,v,rms)
    use var_kind_def
    use var_parallel
    implicit none
    integer(dpI),intent(in):: N
    real   (dpR),intent(in):: v(*)
    real   (dpR):: rms,local(2),global(2)
    real   (kind=8),external:: DNRM2

    local(1)=  real(N, dpR)
    local(2)=  DNRM2(N, v, 1)**2
    if(nprc .le. 1) then
        global  =  local
    else
        call mpi_allreduce(local, global, 2, mpi_dpR, mpi_sum, mpi_comm_world, mpi_err)
    end if
    rms =  sqrt(global(2)/global(1))

    return
    end subroutine get_rms
!-------------------------------------------------------------------------------
!   stop the program.
!-------------------------------------------------------------------------------
    subroutine mbt_stop
    use var_parallel
    implicit none

    call mpi_barrier(mpi_comm_world, mpi_err)
    call mpi_finalize(mpi_err)
    stop

    return
    end subroutine mbt_stop
