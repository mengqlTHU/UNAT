!-------------------------------------------------------------------------------
!   QR decomposition with Modified Gram-Schmidt method.
!-------------------------------------------------------------------------------
    subroutine QR_mgs(is_pivot,LDA,n_col,A,Q,R_diag)
    use var_kind_def
    implicit none
    logical(dpL),intent(in):: is_pivot
    integer(dpI),intent(in):: LDA,n_col
    real   (dpR),intent(in):: A(LDA,*)
    integer(dpI):: piv(n_col),i,j,k,m,itmp
    real   (dpR):: Q(LDA,*),R_diag(*),rnorm(n_col),R(n_col,n_col),dtmp
    integer(dpI),external:: IDAMAX
    real   (dpR),external:: DNRM2,DDOT

    call DCOPY(LDA*n_col, A, 1, Q, 1)
    do j=1,n_col
        rnorm(j)=  DNRM2(LDA, Q(1,j), 1)**2
        piv  (j)=  j
    end do
    if(is_pivot) then
        k   =  maxloc(rnorm, dim=1)
    else
        k   =  1
    end if

    R   =  0.0d0
    do m=1,n_col
        if(k .ne. m) then
            call DSWAP(LDA, Q(1,m), 1, Q(1,k), 1)

            dtmp    =  rnorm(m)
            rnorm(m)=  rnorm(k)
            rnorm(k)=  dtmp

            itmp    =  piv(m)
            piv(m)  =  piv(k)
            piv(k)  =  itmp
        end if

        R(m    ,m)  =  sqrt(rnorm(m))
        Q(1:LDA,m)  =  Q(1:LDA,m)/R(m,m)
        do i=m+1,n_col
            R(m    ,i)  =  DDOT(LDA, Q(1,m), 1, Q(1,i), 1)
            Q(1:LDA,i)  =  Q(1:LDA,i)-Q(1:LDA,m)*R(m,i)
        end do

        do i=m+1,n_col
!           rnorm(i)=  rnorm(i)-R(m,i)*R(m,i)
            rnorm(i)=  DNRM2(LDA, Q(1,i), 1)**2
        end do
        if(m .lt. n_col) then
            if(is_pivot) then
                k   =  IDAMAX(n_col-m, rnorm(m+1:n_col), 1)
                k   =  k+m
            else
                k   =  m+1
            end if
        end if
    end do
    forall(i=1:n_col)   R_diag(i)   =  R(i,i)

    return
    end subroutine QR_mgs
!-------------------------------------------------------------------------------
!   QR decomposition with Modified Gram-Schmidt method, parallel.
!-------------------------------------------------------------------------------
    subroutine QR_MGS_parallel(is_pivot,LDA,n_col,A,Q,R_diag)
    use var_kind_def
    use var_parallel
    implicit none
    logical(dpL),intent(in):: is_pivot
    integer(dpI),intent(in):: LDA,n_col
    real   (dpR),intent(in):: A(LDA,*)
    integer(dpI):: piv(n_col),i,j,k,m,itmp
    real   (dpR):: Q(LDA,*),R_diag(*),rnorm(n_col),R(n_col,n_col),v(n_col),dtmp
    integer(dpI),external:: IDAMAX
    real   (dpR),external:: DNRM2,DDOT,DDOT_MPI

    call DCOPY(LDA*n_col, A, 1, Q, 1)
    do j=1,n_col
        rnorm(j)=  DNRM2(LDA, Q(1,j), 1)**2
        piv  (j)=  j
    end do
    if(nprc .gt. 1) then
        call mpi_allreduce(rnorm, v, n_col, mpi_dpR, mpi_sum, mpi_comm_world, mpi_err)
        rnorm   =  v
    end if

    if(is_pivot) then
        k   =  maxloc(rnorm, dim=1)
    else
        k   =  1
    end if

    R   =  0.0d0
    do m=1,n_col
        if(k .ne. m) then
            call DSWAP(LDA, Q(1,m), 1, Q(1,k), 1)

            dtmp    =  rnorm(m)
            rnorm(m)=  rnorm(k)
            rnorm(k)=  dtmp

            itmp    =  piv(m)
            piv(m)  =  piv(k)
            piv(k)  =  itmp
        end if

        R(m    ,m)  =  sqrt(rnorm(m))
        Q(1:LDA,m)  =  Q(1:LDA,m)/R(m,m)

        do i=m+1,n_col
            v(i)=  DDOT(LDA, Q(1,m), 1, Q(1,i), 1)
        end do
        if((nprc .gt. 1) .and. (m .lt. n_col)) then
            call mpi_allreduce(v(m+1), R(m,m+1:n_col), n_col-m, mpi_dpR, mpi_sum, &
                &  mpi_comm_world, mpi_err)
        else
            R(m,m+1:n_col)  =  v(m+1:n_col)
        end if
        do i=m+1,n_col
            Q(1:LDA,i)  =  Q(1:LDA,i)-Q(1:LDA,m)*R(m,i)
!           rnorm(i)    =  rnorm(i)-R(m,i)*R(m,i)
            rnorm(i)    =  DNRM2(LDA, Q(1,i), 1)**2
        end do

        if((nprc .gt. 1) .and. (m .lt. n_col)) then
            call mpi_allreduce(rnorm(m+1), v(m+1), n_col-m, mpi_dpR, mpi_sum, &
                &  mpi_comm_world, mpi_err)
            rnorm(m+1:n_col)=  v(m+1:n_col)
        end if

        if(m .lt. n_col) then
            if(is_pivot) then
                k   =  IDAMAX(n_col-m, rnorm(m+1:n_col), 1)
                k   =  k+m
            else
                k   =  m+1
            end if
        end if
    end do
    forall(i=1:n_col)   R_diag(i)   =  R(i,i)

    return
    end subroutine QR_MGS_parallel
!-------------------------------------------------------------------------------
!   module for Recursive Projection Method.
!-------------------------------------------------------------------------------
    module var_rpm
        use var_kind_def
        implicit none
        private
        public:: type_rpm,rpm_setup,rpm_stabilize_solution

        real(dpR),parameter:: eig_ratio     =  1.0d2
        real(dpR),parameter:: eig_amplitude =  1.0d4
        real(dpR),parameter:: delta         =  1.0d-1

        type type_rpm
            logical(dpL):: is_scaled        =  .false.
            integer(dpI):: max_unstable     =  0
            integer(dpI):: bsize            =  1
            integer(dpI):: LDA              =  0
            integer(dpI):: LDA_global       =  0
            integer(dpI):: n_dq             =  3
            integer(dpI):: n_unstable       =  0
            integer(dpI):: iter             =  0
            integer(dpI):: step             =  0
            real   (dpR):: ref(100)         =  1.0d0

            real   (dpR),allocatable:: H(:,:),inv_H(:,:)
            real   (dpR),allocatable:: q(:),z_new(:),z_old(:)
            real   (dpR),allocatable:: x(:)
            real   (dpR),allocatable:: dq(:,:),unstable_mode(:,:)
            real   (dpR),allocatable:: QR_q(:,:),QR_eig(:)
            real   (dpR),allocatable:: H_q(:,:),H_eig(:)
        end type type_rpm
        contains
!       ------------------------------------------------------------------------
!       setup the RPM.
!       ------------------------------------------------------------------------
        subroutine rpm_setup(p,LDA,max_unstable,bsize,ref)
        use var_parallel
        implicit none
        integer(dpI),intent(in):: LDA,max_unstable,bsize
        real   (dpR),intent(in):: ref(*)
        real   (dpR):: a,b
        type(type_rpm):: p

        if(p%max_unstable .gt. 0)   stop 'Error: this RPM has been used.'
        p%max_unstable  =  max(max_unstable, 2)
        p%LDA           =  LDA
        p%bsize         =  bsize
        allocate(p%H    (max_unstable, max_unstable ))
        allocate(p%inv_H(max_unstable, max_unstable ))
        allocate(p%H_q  (max_unstable, max_unstable ))
        allocate(p%H_eig(max_unstable               ))
        allocate(p%z_old(p%max_unstable             ))
        allocate(p%z_new(p%max_unstable             ))
        allocate(p%q            (LDA                ))
        allocate(p%x            (LDA                ))
        allocate(p%unstable_mode(LDA, p%max_unstable))
        allocate(p%dq           (LDA, p%n_dq        ))
        allocate(p%QR_q         (LDA, p%n_dq        ))
        allocate(p%QR_eig       (     p%n_dq        ))

        a   =  maxval(ref(1:p%bsize))
        b   =  minval(ref(1:p%bsize))
        if((abs(a-b) .gt. 1.0d-10*a) .and. (p%bsize .gt. 1)) then
            p%is_scaled     =  .true.
            p%ref(1:p%bsize)=  ref(1:p%bsize)
        end if

        if(nprc .gt. 1) then
        call mpi_allreduce(p%LDA,p%LDA_global,1,mpi_dpI,mpi_sum,mpi_comm_world,mpi_err)
        else
            p%LDA_global=  p%LDA
        end if

        return
        end subroutine rpm_setup
!       ------------------------------------------------------------------------
!       stabilize the solution from external solver.
!       ------------------------------------------------------------------------
        subroutine rpm_stabilize_solution(p,is_pre,exchange_solution,iteration)
        implicit none
        logical(dpL),intent(in):: is_pre
        type(type_rpm),intent(in):: p
        external:: exchange_solution,iteration

        call exchange_solution(p%is_scaled, p%ref, 1, p%x)
        call stabilize_solution(p, is_pre, p%x, exchange_solution, iteration)
        call exchange_solution(p%is_scaled, p%ref, 2, p%x)

        return
        end subroutine rpm_stabilize_solution
!       ------------------------------------------------------------------------
!       stabilize the solution from external solver.
!       ------------------------------------------------------------------------
        subroutine stabilize_solution(p,is_pre,x,exchange_solution,iteration)
        use var_parallel
        implicit none
        logical(dpL),intent(in):: is_pre
        logical(dpL):: ltmp
        integer(dpI):: i,M,iz,itmp
        real   (dpR):: x(*),y(100),z(100),v(100)
        type(type_rpm):: p
        real   (dpR),external:: DDOT
        external:: exchange_solution,iteration

        M   =  p%n_unstable

!       compute Z^Tx
        do i=1,M
            v(i)=  DDOT(p%LDA, p%unstable_mode(1,i), 1, x, 1)
        end do
        if(nprc .gt. 1) then
            call mpi_allreduce(v, z, M, mpi_dpR, mpi_sum, mpi_comm_world, mpi_err)
        else
            z(1:M)  =  v(1:M)
        end if

        if(is_pre) then
            p%z_old(1:M)=  z(1:M)
            do i=1,p%n_dq-1
                call DCOPY(p%LDA, p%dq(1,i+1), 1, p%dq(1,i), 1)
            end do
            if(M .le. 0)    call DCOPY(p%LDA, x, 1, p%dq(1,p%n_dq), 1)

            return
        end if

        p%iter  =  p%iter+1
        p%step  =  p%step+1
        if(M .le. 0) then
            p% q(1:p%LDA       )=  x(1:p%LDA)
            p%dq(1:p%LDA,p%n_dq)=  x(1:p%LDA)-p%dq(1:p%LDA,p%n_dq)
        else
!           stabilize the solution.
            p%z_new(1:M)=  p%z_old(1:M)+matmul(p%inv_H(1:M,1:M), z(1:M)-p%z_old(1:M))
            p%dq(1:p%LDA,p%n_dq)=  p%q(1:p%LDA)
            p%q (1:p%LDA       )=  x  (1:p%LDA)
            do i=1,M
                p%q(1:p%LDA)=  p%q(1:p%LDA)-z(i)*p%unstable_mode(1:p%LDA,i)
            end do
            p%dq(1:p%LDA,p%n_dq)=  p%q(1:p%LDA)-p%dq(1:p%LDA,p%n_dq)
            x(1:p%LDA)  =  p%q(1:p%LDA)
            do i=1,M
                x(1:p%LDA)  =  x(1:p%LDA)+p%z_new(i)*p%unstable_mode(1:p%LDA,i)
            end do
        end if

!       ------------------------------------------------------------------------
!       increase the unstable basis.
        if(p%step .le. p%n_dq+1)    return
!       if((p%iter .le. p%n_dq+1) .or. (mod(p%iter, p%n_dq*3) .ne. 0))  return

        call QR_mgs_parallel(.true., p%LDA, p%n_dq, p%dq, p%QR_q, p%QR_eig)
        iz  =  0
        do i=1,p%n_dq
            if(p%QR_eig(i) .ge. eig_amplitude) then
                iz  =  max(iz, i)
            elseif(i .lt. p%n_dq) then
                if(p%QR_eig(i) .ge. p%QR_eig(i+1)*eig_ratio)    iz  =  max(iz, i)
            end if
        end do
        if(iz .le. 0)   return

        itmp=  p%n_unstable+iz-p%max_unstable
        if(itmp .gt. 0) then
            do i=itmp+1,p%n_unstable
                call DCOPY(p%LDA, p%unstable_mode(1,i), 1, p%unstable_mode(1,i-itmp), 1)
            end do
            call DCOPY(p%LDA*iz, p%QR_q, 1, p%unstable_mode(1,p%n_unstable-itmp+1), 1)
            p%n_unstable=  p%max_unstable
        else
            call DCOPY(p%LDA*iz, p%QR_q, 1, p%unstable_mode(1,p%n_unstable+1), 1)
            p%n_unstable=  p%n_unstable+iz
        end if
        call get_H(p, x, exchange_solution, iteration)

        if(.false.) then
            iz  =  p%n_unstable
            call QR_mgs(.true., iz, iz, p%H(1:iz,1:iz), p%H_q(1:iz,1:iz), p%H_eig)
            ltmp=  .false.
            do i=1,iz
                ltmp=  abs(p%H_eig(i)) .lt. 1.0d0-delta
                if(ltmp)    exit
            end do
            if(ltmp) then
                p%n_unstable=  i-1
                if(p%n_unstable .le. 0) return
                do i=1,p%LDA
                    z(1:iz) =  p%unstable_mode(i,1:iz)
                    call DGEMM('N', 'N', 1, p%n_unstable, iz, 1.0d0, z, 1, &
                        &  p%H_q(1:iz,1:iz), iz, 0.0d0, y, 1)
                    p%unstable_mode(i,1:p%n_unstable)   =  y(1:p%n_unstable)
                end do
            end if
            call get_H(p, x, exchange_solution, iteration)
        end if

        do i =1,p%n_unstable
        do iz=1,p%n_unstable
            if(i .eq. iz) then
                p%inv_H(i,iz)   =  1.0d0-p%H(i,i)
            else
                p%inv_H(i,iz)   =       -p%H(i,iz)
            end if
        end do
        end do
        call mat_inv(p%n_unstable, p%inv_H(1:p%n_unstable, 1:p%n_unstable))

        p%step  =  0
!       increase the unstable basis.
!       ------------------------------------------------------------------------

        return
        end subroutine stabilize_solution
!       ------------------------------------------------------------------------
!       cal Z^T\fpp{F}{u}Z.
!       ------------------------------------------------------------------------
        subroutine get_H(p,x,exchange_solution,iteration)
        use var_parallel
        implicit none
        real   (dpR),intent(in):: x(*)
        integer(dpI):: iz,i
        real   (dpR):: e,v(2,100),w(2,100)
        type(type_rpm):: p
        real   (dpR),external:: DDOT,DNRM_MPI
        external:: exchange_solution,iteration

!       conduct one iteration and save F(u) in p%QR_q(:,1)
        call exchange_solution(p%is_scaled, p%ref, 2, x)
        call iteration
        call exchange_solution(p%is_scaled, p%ref, 1, p%QR_q(1,1))

!       e   =  sqrt(epsilon(1.0d0))*DNRM2(p%LDA, x, 1)/sqrt(real(p%LDA, dpR))
        e   =  sqrt(epsilon(1.0d0))*DNRM_MPI(p%LDA, x)/sqrt(real(p%LDA_global, dpR))
        do iz=1,p%n_unstable
!           set u<--u+e*unstable_mode
            p%QR_q(1:p%LDA,2)   =  x(1:p%LDA)+e*p%unstable_mode(1:p%LDA,iz)
!           conduct one iteration and save F(u) in p%QR_q(:,2)
            call exchange_solution(p%is_scaled, p%ref, 2, p%QR_q(1,2))
            call iteration
            call exchange_solution(p%is_scaled, p%ref, 1, p%QR_q(1,2))

            do i=1,p%n_unstable
                v(1,i)  =  DDOT(p%LDA, p%unstable_mode(1,i), 1, p%QR_q(1,2), 1)
                v(2,i)  =  DDOT(p%LDA, p%unstable_mode(1,i), 1, p%QR_q(1,1), 1)
            end do
            if(nprc .gt. 1) then
                call mpi_allreduce(v, w, 2*p%n_unstable, mpi_dpR, mpi_sum, &
                    &  mpi_comm_world, mpi_err)
                call DCOPY(2*p%n_unstable, w, 1, v, 1)
            end if
            do i=1,p%n_unstable
                p%H(i,iz)   = (v(1,i)-v(2,i))/e
            end do
        end do

        return
        end subroutine get_H
    end module var_rpm
