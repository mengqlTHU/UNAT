!-------------------------------------------------------------------------------
!   linked list.
!-------------------------------------------------------------------------------
    module var_list
        implicit none
        private
        public:: type_list,list_get_nele,list_destroy
        public:: list_add_eleR,list_get_eleR,list_add_ele,list_get_ele

        integer(kind=4),parameter:: dpL =  kind(.true.)
        integer(kind=4),parameter:: dpI =  kind(1)
        integer(kind=4),parameter:: dpR =  kind(1.0d0)

        type type_list
            integer(dpI):: j=  0
            real   (dpR),allocatable:: A(:)
            type(type_list),pointer:: next  => null()
        end type type_list
        contains
!       ------------------------------------------------------------------------
!       add element to the list.
!       ------------------------------------------------------------------------
        subroutine list_add_ele(p,j)
        implicit none
        integer(dpI),intent(in):: j
        type(type_list),pointer:: p
        type(type_list),pointer:: now=>null()

        if(.not. associated(p)) then
            allocate(p)
            p%j =  j
            return
        end if

        now => p
        do while(.true.)
            if(now%j .eq. j) then
                return
            elseif(associated(now%next)) then
                now => now%next
            else
                allocate(now%next)
                now%next%j  =  j
                return
            end if
        end do

        return
        end subroutine list_add_ele
!       ------------------------------------------------------------------------
!       get the elements.
!       ------------------------------------------------------------------------
        subroutine list_get_ele(p,ibuf)
        implicit none
        type(type_list),pointer,intent(in):: p
        integer   (dpI):: ibuf(*),i
        type(type_list),pointer:: now

        now => p
        i   =  0
        do while(associated(now))
            i   =  i+1
            ibuf(i) =  now%j
            now => now%next
        end do

        return
        end subroutine list_get_ele
!       ------------------------------------------------------------------------
!       add element to the list, together with the real entry.
!       ------------------------------------------------------------------------
        subroutine list_add_eleR(p,j,LDA,A)
        implicit none
        integer(dpI),intent(in):: j,LDA
        real   (dpR),intent(in):: A(*)
        type(type_list),pointer:: p
        type(type_list),pointer:: now=>null()

        if(.not. associated(p)) then
            allocate(p)
            allocate(p%A(LDA))
            p%j =  j
            call DCOPY(LDA, A, 1, p%A, 1)
            return
        end if

        now => p
        do while(.true.)
            if(now%j .eq. j) then
                call DAXPY(LDA, 1.0d0, A, 1, now%A, 1)
                return
            elseif(associated(now%next)) then
                now => now%next
            else
                allocate(now%next)
                allocate(now%next%A(LDA))
                now%next%j  =  j
                call DCOPY(LDA, A, 1, now%next%A, 1)
                return
            end if
        end do

        return
        end subroutine list_add_eleR
!       ------------------------------------------------------------------------
!       get the elements, together with the real entries.
!       ------------------------------------------------------------------------
        subroutine list_get_eleR(p,ibuf,LDA,rbuf)
        implicit none
        type(type_list),pointer,intent(in):: p
        integer   (dpI):: ibuf(*),i,LDA
        real      (dpR):: rbuf(*)
        type(type_list),pointer:: now

        now => p
        i   =  0
        do while(associated(now))
            i   =  i+1
            ibuf(i) =  now%j
            if(.not. allocated(now%A))  stop 'Error: element of list not defined.'
            rbuf(1+LDA*(i-1):LDA*i) =  now%A(1:LDA)
            now => now%next
        end do

        return
        end subroutine list_get_eleR
!       ------------------------------------------------------------------------
!       get the number of elements.
!       ------------------------------------------------------------------------
        subroutine list_get_nele(p,nele)
        implicit none
        type(type_list),pointer,intent(in):: p
        integer   (dpI):: nele
        type(type_list),pointer:: now=>null()

        now => p
        nele=  0
        do while(associated(now))
            nele=  nele+1
            now => now%next
        end do

        return
        end subroutine list_get_nele
!       ------------------------------------------------------------------------
!       destroy the list.
!       ------------------------------------------------------------------------
        recursive subroutine list_destroy(p)
        implicit none
        type(type_list),pointer:: p

        if(associated(p%next))  call list_destroy(p%next)
        if(allocated (p%A   ))  deallocate(p%A)
        deallocate(p)

        end subroutine list_destroy
    end module var_list
!-------------------------------------------------------------------------------
!   module for the preconditioner.
!-------------------------------------------------------------------------------
    module var_prec
        use var_list
        implicit none
        private
        public:: type_prec,slv_SGS_MF
        public:: prec_add,prec_add_ele,prec_add_eleR,prec_add_eleR_partial
        public:: prec_set_base,prec_set_ordering,prec_set_coarse
        public:: prec_set_lines
        public:: prec_get_decomposition,prec_solve
        public:: prec_delete

        integer(kind=4),parameter:: dpL =  kind(.true.)
        integer(kind=4),parameter:: dpI =  kind(1)
        integer(kind=4),parameter:: dpR =  kind(1.0d0)

        logical(dpL),parameter:: is_reorder     =  .true.
        integer(dpI),parameter:: max_level      =  9
        integer(dpI),parameter:: slv_SGS        =  1
        integer(dpI),parameter:: slv_SGS_MF     =  2
        integer(dpI),parameter:: slv_ILU0       =  3
        integer(dpI),parameter:: slv_LINE       =  4
        integer(dpI),parameter:: slv_DIRECT     =  9
        integer(dpI),parameter:: nvtx_direct    =  0
        integer(dpI):: err_mem                  =  0

        type type_pointer_list
            type(type_list),pointer:: p => null()
        end type type_pointer_list

        type type_prec
            logical(dpL):: is_spy       =  .false.
            logical(dpL):: is_reordered =  .false.
            logical(dpL):: is_coarsened =  .false.
            logical(dpL):: is_damped    =  .false.
            logical(dpL):: is_LU        =  .true.
            integer(dpI):: nvtx         =  0
            integer(dpI):: bsize        =  1
            integer(dpI):: LDA          =  0
            integer(dpI):: solver       =  0
            real   (dpR):: relx_max     =  1.0d0
            real   (dpR):: relx_min     =  0.5d0
            real   (dpR):: relx_min_mg  =  1.0d0
            integer(dpI),allocatable:: iA(:),jA(:)
            integer(dpI),allocatable:: perm(:)
            real   (dpR),allocatable:: A(:,:),D(:,:),LU(:,:)

            type(type_pointer_list),allocatable:: row(:)

            real   (dpR),allocatable:: RHS(:),x(:)
            real   (dpR),allocatable:: buf1(:),buf2(:),buf3(:)

            integer(dpI):: lev  =  0
            integer(dpI):: mlev =  1
            integer(dpI),allocatable:: ci(:)
            type(type_prec),pointer:: clev  => null()

            integer(dpI):: n_line       =  0
            integer(dpI):: n_vtx_line   =  0
            integer(dpI),allocatable:: iA_line(:),jA_line(:),ID_LU_line(:,:)
            real   (dpR),allocatable:: L_line(:,:),D_line(:,:),U_line(:,:),R_line(:,:)
        end type type_prec

        contains
!       ------------------------------------------------------------------------
!       add a new prec.
!       ------------------------------------------------------------------------
        subroutine prec_add(p,nvtx,bsize,solver)
        implicit none
        integer(dpI),intent(in):: nvtx,bsize,solver
        type(type_prec):: p

        p%nvtx  =  nvtx
        p%bsize =  bsize
        p%LDA   =  nvtx*bsize
        if(p%LDA .le. nvtx_direct) then
            p%solver=  slv_direct
        else
            p%solver=  solver
        end if
        allocate(p%RHS (p%LDA), stat=err_mem)
        allocate(p%x   (p%LDA), stat=err_mem)
        allocate(p%buf1(p%LDA), stat=err_mem)
        allocate(p%buf2(p%LDA), stat=err_mem)
        allocate(p%buf3(p%LDA), stat=err_mem)
        if(solver .eq. slv_SGS) then
            allocate(p%D(p%bsize, p%LDA), stat=err_mem)
        elseif(solver .eq. slv_SGS_MF) then
            allocate(p%A(p%bsize, p%LDA), stat=err_mem)
            allocate(p%D(p%bsize, p%LDA), stat=err_mem)
        elseif(solver .eq. slv_LINE) then
            allocate(p%D(p%bsize, p%LDA), stat=err_mem)
        end if
        p%RHS   =  0.0d0
        p%x     =  0.0d0
        p%buf1  =  0.0d0
        p%buf2  =  0.0d0
        p%buf3  =  0.0d0

        return
        end subroutine prec_add
!       ------------------------------------------------------------------------
!       add an matrix element to the prec, only the (i,j) coordinates.
!       ------------------------------------------------------------------------
        subroutine prec_add_ele(p,i,j)
        implicit none
        integer(dpI),intent(in):: i,j
        type(type_prec):: p

        if((i .le. 0) .or. (i .gt. p%nvtx)) stop 'Error: wrong index for prec_add_ele.'
        if((j .le. 0) .or. (j .gt. p%nvtx)) stop 'Error: wrong index for prec_add_ele.'
        if(p%is_spy) then
            stop 'Error: the sparse pattern of LHS is already fixed.'
        else
            if(.not. allocated(p%row))  allocate(p%row(p%nvtx))
            call list_add_ele(p%row(i)%p, j)
        end if

        return
        end subroutine prec_add_ele
!       ------------------------------------------------------------------------
!       add an matrix element to the prec, real.
!       ------------------------------------------------------------------------
        recursive subroutine prec_add_eleR(p,i0,j0,LHS)
        implicit none
        integer(dpI),intent(in):: i0,j0
        real   (dpR),intent(in):: LHS(*)
        integer(dpI):: i,j,k,M
        type(type_prec):: p

        if(p%is_reordered) then
            i   =  p%perm(i0)
            j   =  p%perm(j0)
        else
            i   =  i0
            j   =  j0
        end if

        if((i .le. 0) .or. (i .gt. p%nvtx)) stop 'Error: wrong index for prec_add_eleR.'
        if((j .le. 0) .or. (j .gt. p%nvtx)) stop 'Error: wrong index for prec_add_eleR.'
        M   =  p%bsize**2
        if((p%lev .eq. 0) .and. (p%solver .eq. slv_SGS_MF)) then
!           only the diagonal element is saved.
            if(i .eq. j)    call DAXPY(M, 1.0d0, LHS, 1, p%A(1,1+p%bsize*(i-1)), 1)
!           also we have to compute LHS on level(1).
            if(p%is_coarsened)  call prec_add_eleR(p%clev, p%ci(i), p%ci(j), LHS)
        else
            if(p%is_spy) then
                do k=p%iA(i),p%iA(i+1)-1
                    if(p%jA(k) .ne. j)  cycle
                    if(p%bsize .eq. 1) then
                        p%A(1,k)=  p%A(1,k)+LHS(1)
                    else
                        call DAXPY(M, 1.0d0, LHS, 1, p%A(1,1+p%bsize*(k-1)), 1)
                    end if
                    return
                end do
                stop 'Error: fails to find (i,j) in the SPY.'
            else
                if(.not. allocated(p%row))  allocate(p%row(p%nvtx), stat=err_mem)
                call list_add_eleR(p%row(i)%p, j, M, LHS)
            end if
        end if

        return
        end subroutine prec_add_eleR
!       ------------------------------------------------------------------------
!       add an matrix element to the prec, real, only a part of the entry.
!       ------------------------------------------------------------------------
        subroutine prec_add_eleR_partial(p,i0,j0,L,R,LHS)
        implicit none
        integer(dpI),intent(in):: i0,j0,L(*),R(*)
        real   (dpR),intent(in):: LHS(*)
        integer(dpI):: i,j,k,m,ii,jj
        type(type_prec):: p

        if(p%is_reordered) then
            i   =  p%perm(i0)
            j   =  p%perm(j0)
        else
            i   =  i0
            j   =  j0
        end if

        if((i .le. 0) .or. (i .gt. p%nvtx)) stop 'Error: wrong index for prec_add_eleR.'
        if((j .le. 0) .or. (j .gt. p%nvtx)) stop 'Error: wrong index for prec_add_eleR.'
        if(.not. p%is_spy)  stop 'Error: the SPY of LHS must be ready.'

        do k=p%iA(i),p%iA(i+1)-1
            if(p%jA(k) .eq. j) then
                m   =  0
                do jj=L(2),R(2)
                do ii=L(1),R(1)
                    m   =  m+1
                    p%A(ii,jj+p%bsize*(k-1))=  p%A(ii,jj+p%bsize*(k-1))+LHS(m)
                end do
                end do
                return
            end if
        end do
        stop 'Error: fails to find (i,j) in the SPY.'

        return
        end subroutine prec_add_eleR_partial
!       ------------------------------------------------------------------------
!       set the base level.
!       ------------------------------------------------------------------------
        subroutine prec_set_base(p,is_set_LHS)
        implicit none
        logical(dpL),intent(in):: is_set_LHS
        integer(dpI):: i,j,LDA
        type(type_prec):: p

        if(p%is_spy)    return
        allocate(p%iA(p%nvtx+1), stat=err_mem)
        do i=1,p%nvtx
            call list_get_nele(p%row(i)%p, p%iA(i+1))
        end do
        j   =  0
        do i=1,p%nvtx
            j   =  j+p%iA(i+1)
        end do
        allocate(p%jA(j), stat=err_mem)
        allocate(p%A(p%bsize,p%bsize*j), stat=err_mem)
        p%A =  0.0d0

        p%iA(1) =  1
        LDA     =  p%bsize**2
        do i=1,p%nvtx
            j           =  p%iA(i)
            p%iA(i+1)   =  p%iA(i+1)+j
            if(is_set_LHS) then
                call list_get_eleR(p%row(i)%p, p%jA(j), LDA, p%A(1,1+p%bsize*(j-1)))
                call iqsortcols_mat(.true., j, p%iA(i+1)-1, 1, 1, p%jA, LDA, p%A)
            else
                call list_get_ele (p%row(i)%p, p%jA(j))
            end if
            call list_destroy(p%row(i)%p)
        end do
        if(allocated(p%row))    deallocate(p%row)
        p%is_spy=  .true.

        return
        end subroutine prec_set_base
!       ------------------------------------------------------------------------
!       set the SPY of coarse level.
!       ------------------------------------------------------------------------
        subroutine prec_set_coarse(p,lev,ci)
        implicit none
        integer(dpI),intent(in):: lev,ci(*)
        logical(dpL):: is_coarsened
        integer(dpI):: iter
        type(type_prec):: p

        if(lev .lt. 0)  stop 'Error: wrong input for prec_set_coarse.'
        iter=  0
        call prec_coarsen_node(is_coarsened, p)
        if(is_coarsened) then
            p%mlev  =  p%mlev+1
            call set_mlev(p%clev, p%mlev)
        end if

        return
        contains
!           --------------------------------------------------------------------
!           setup the coarse level.
!           --------------------------------------------------------------------
            recursive subroutine prec_coarsen_node(is_coarsened,p)
            implicit none
            logical(dpL):: is_coarsened
            type(type_prec):: p

            is_coarsened=  .false.
            if(p%nvtx .le. nvtx_direct) return

            iter=  iter+1
            if(iter .lt. lev+1) then
                if(.not. associated(p%clev))    stop "Error: fine lev of prec isn't OK."
                call prec_coarsen_node(is_coarsened, p%clev)
            elseif(iter .eq. lev+1) then
                if(associated(p%clev)) then
                    stop 'Error: coarse level of prec was OK.'
                else
                    allocate(p%clev)
                    allocate(p%ci(p%nvtx))
                    call ICOPY(p%nvtx, ci, 1, p%ci, 1)
                end if
                call prec_set_clev_spy(p, p%ci, p%clev)
                is_coarsened=  .true.
            end if

            return
            end subroutine prec_coarsen_node
!           --------------------------------------------------------------------
!           setup the SPY of coarse level.
!           --------------------------------------------------------------------
            subroutine prec_set_clev_spy(p0,ci,p1)
            implicit none
            integer(dpI),intent(in):: ci(*)
            integer(dpI):: i,j,k,m,nvtx
            type(type_prec):: p0,p1
            integer(dpI),allocatable:: jA(:,:)

            if(p0%nvtx .le. nvtx_direct)    stop 'Error: lowest threshold of MG reached.'

            nvtx=  0
            do i=1,p0%nvtx
                if(ci(i) .le. 0)    stop 'Error: wrong coarse level index.'
                nvtx=  max(nvtx, ci(i))
            end do
!           call prec_add(p1, nvtx, p0%bsize, p0%solver)
            call prec_add(p1, nvtx, p0%bsize, slv_ILU0)
            p1%lev          =  p0%lev+1
            p1%relx_min     =  p0%relx_min
            p1%relx_min_mg  =  p0%relx_min_mg
            p1%is_damped    =  p0%is_damped

            allocate(p1%iA(p1%nvtx+1), stat=err_mem)
            allocate(jA(2,p0%iA(p0%nvtx+1)-1), stat=err_mem)
            do i=1,p0%nvtx
                do j=p0%iA(i),p0%iA(i+1)-1
                    jA(1:2,j)   = (/ci(i), ci(p0%jA(j))/)
                end do
            end do
            call iqsortcols(.true., 1, p0%iA(p0%nvtx+1)-1, 1, 2, jA)
            p1%iA   =  0
            i       =  1
            do while(i .le. p0%iA(p0%nvtx+1)-1)
                k   =  i
                do j=i,p0%iA(p0%nvtx+1)-1
                    if(jA(1,j) .ne. jA(1,i)) then
                        exit
                    else
                        k   =  j
                    end if
                end do
                call iqsortcols(.true., i, k, 2, 2, jA)
                do j=i+1,k
!                   delete duplicate adjacency entry.
                    if(abs(jA(2,j)) .eq. abs(jA(2,j-1)))    jA(2,j) = -jA(2,j)
                end do
                do j=i,k
                    if(jA(2,j) .gt. 0)  p1%iA(jA(1,i)+1)=  p1%iA(jA(1,i)+1)+1
                end do

                i   =  k+1
            end do
            p1%iA(1)=  1
            do i=1,p1%nvtx
                if(p1%iA(i+1) .le. 0)   stop 'Error: fails to get coarse level SPY.'
                p1%iA(i+1)  =  p1%iA(i+1)+p1%iA(i)
            end do
            allocate(p1%jA(p1%iA(p1%nvtx+1)-1), stat=err_mem)
            allocate(p1% A(p1%bsize, p1%bsize*(p1%iA(p1%nvtx+1)-1)), stat=err_mem)
            i   =  1
            do while(i .le. p0%iA(p0%nvtx+1)-1)
                k   =  i
                do j=i,p0%iA(p0%nvtx+1)-1
                    if(jA(1,j) .ne. jA(1,i)) then
                        exit
                    else
                        k   =  j
                    end if
                end do
                m   =  0
                do j=i,k
                    if(jA(2,j) .le. 0)  cycle
                    m   =  m+1
                    p1%jA(m-1+p1%iA(jA(1,i)))   =  jA(2,j)
                end do

                i   =  k+1
            end do
            if(allocated(jA))   deallocate(jA)
            p1%is_spy       =  .true.
            p0%is_coarsened =  .true.

            return
            end subroutine prec_set_clev_spy
!           --------------------------------------------------------------------
!           set mlev.
!           --------------------------------------------------------------------
            recursive subroutine set_mlev(p,mlev)
            implicit none
            integer(dpI),intent(in):: mlev
            type(type_prec):: p

            p%mlev  =  mlev
            if(associated(p%clev))  call set_mlev(p%clev, mlev)

            return
            end subroutine set_mlev
        end subroutine prec_set_coarse
!       ------------------------------------------------------------------------
!       reordering the linear system.
!       ------------------------------------------------------------------------
        subroutine prec_set_ordering(p)
        implicit none
        type(type_prec):: p
        integer(dpI),allocatable:: ci(:)

        if(p%lev .ne. 0)    stop 'Error: wrong input for prec_set_ordering.'
        if((.not. is_reorder) .or. p%is_reordered)  return
        call reorder(p)
        if(p%mlev .gt. 1) then
            allocate(ci(p%nvtx), stat=err_mem)
            call reset_ci(p, p%clev)
            if(allocated(ci))   deallocate(ci)
        end if
!       delete the perm infomation of coarse level to avoid confusion.
        call delete_coarse_perm(p)
        if((p%lev .eq. 0) .and. (p%solver .eq. slv_SGS_MF)) deallocate(p%perm)

        return
        contains
!       ------------------------------------------------------------------------
!       reorder the matrix.
!       ------------------------------------------------------------------------
        recursive subroutine reorder(p)
        implicit none
        integer(dpI):: i
        type(type_prec):: p

        if(p%is_reordered)  return

        if((p%lev .eq. 0) .and. (p%solver .eq. slv_SGS_MF)) then
            allocate(p%perm(p%nvtx), stat=err_mem)
            forall(i=1:p%nvtx)  p%perm(i)   =  i
        else
            if(.not. allocated(p%iA))   stop 'Error: SPY of PREC is not available.'
            allocate(p%perm(p%nvtx), stat=err_mem)
            call matrix_ordering(p%nvtx, p%iA, p%jA, p%perm)
            call matrix_reorder(p%nvtx, p%bsize, p%perm, p%iA, p%jA, p%A)
            p%is_reordered  =  p%lev .eq. 0
        end if

        if(associated(p%clev))  call reorder(p%clev)

        return
        end subroutine reorder
!       ------------------------------------------------------------------------
!       adjust the coarse cell index if the mesh is reordered.
!       ------------------------------------------------------------------------
        recursive subroutine reset_ci(p1,p2)
        implicit none
        integer(dpI):: i
        type(type_prec):: p1,p2

        do i=1,p1%nvtx
            ci(p1%perm(i))  =  p2%perm(p1%ci(i))
        end do
        call ICOPY(p1%nvtx, ci, 1, p1%ci, 1)
        if(associated(p2%clev)) call reset_ci(p2, p2%clev)

        return
        end subroutine reset_ci
!       ------------------------------------------------------------------------
!       delete the perm of coarse level.
!       ------------------------------------------------------------------------
        recursive subroutine delete_coarse_perm(p)
        implicit none
        type(type_prec):: p

        if(.not. associated(p%clev))        return
        if(.not. allocated(p%clev%perm))    return
        deallocate(p%clev%perm)
        p%clev%is_reordered =  .false.
        call delete_coarse_perm(p%clev)

        return
        end subroutine delete_coarse_perm
        end subroutine prec_set_ordering
!       ------------------------------------------------------------------------
!       decompose the matrix before the solution.
!       ------------------------------------------------------------------------
        recursive subroutine prec_get_decomposition(p)
        implicit none
        logical(dpL):: ltmp
        integer(dpI):: i,k,M
        type(type_prec):: p

        M   =  p%bsize**2
        if(p%solver .eq. slv_SGS) then
            do i=1,p%nvtx
            ltmp=  .true.
            do k=p%iA(i),p%iA(i+1)-1
                if(p%jA(k) .ne. i)  cycle
                ltmp=  .false.
                call DCOPY(M, p%A(1,1+p%bsize*(k-1)), 1, p%D(1,1+p%bsize*(i-1)), 1)
                call mat_inv(p%bsize, p%D(1,1+p%bsize*(i-1)))
                exit
            end do
            if(ltmp)    stop 'Error: fails to find the diagonal element.'
            end do
        elseif(p%solver .eq. slv_SGS_MF) then
            if(p%lev .eq. 0) then
                if(allocated(p%iA)) deallocate(p%iA)
                if(allocated(p%jA)) deallocate(p%jA)
            else
                stop 'Error: slv_SGS_MF is supported only on the BASE level.'
            end if
            do i=1,p%nvtx
                call DCOPY(M, p%A(1,1+p%bsize*(i-1)), 1, p%D(1,1+p%bsize*(i-1)), 1)
                call mat_inv(p%bsize, p%D(1,1+p%bsize*(i-1)))
            end do
        elseif(p%solver .eq. slv_ILU0) then
            i   =  p%iA(p%nvtx+1)-1
            if((p%lev .eq. 0) .and. (p%mlev .eq. 1) .and. (.not. p%is_damped)) then
                p%is_LU =  .false.
                if(.not. allocated(p%LU))   allocate(p%LU(p%bsize, 1), stat=err_mem)
                call ilu0_decomposition(p%nvtx, p%iA, p%jA, p%bsize, p%A)
            else
                p%is_LU =  .true.
                if(.not. allocated(p%LU))   allocate(p%LU(p%bsize, p%bsize*i), stat=err_mem)
                call DCOPY(M*i, p%A, 1, p%LU, 1)
                call ilu0_decomposition(p%nvtx, p%iA, p%jA, p%bsize, p%LU)
            end if
        elseif(p%solver .eq. slv_LINE) then
            do i=1,p%nvtx
            ltmp=  .true.
            do k=p%iA(i),p%iA(i+1)-1
                if(p%jA(k) .ne. i)  cycle
                ltmp=  .false.
                call DCOPY(M, p%A(1,1+p%bsize*(k-1)), 1, p%D(1,1+p%bsize*(i-1)), 1)
                call mat_inv(p%bsize, p%D(1,1+p%bsize*(i-1)))
                exit
            end do
            if(ltmp)    stop 'Error: fails to find the diagonal element.'
            end do
            call prec_set_line_LHS(p)
        elseif(p%solver .eq. slv_DIRECT) then
            call direct_decomposition(p)
        else
            stop 'Error: solver not supported, PREC.'
        end if
        if(associated(p%clev)) then
            ltmp= (p%lev .eq. 0) .and. (p%solver .eq. slv_SGS_MF)
!           if the solver of base is slv_SGS_MF, LHS of level(1) should be available.
            if(.not. ltmp)  call get_coarse_LHS(p, p%clev)
            call prec_get_decomposition(p%clev)
        end if

        return
        end subroutine prec_get_decomposition
!       ------------------------------------------------------------------------
!       compute the LU decomposition for small size sparse matrix.
!       ------------------------------------------------------------------------
        subroutine direct_decomposition(p)
        implicit none
        integer(dpI):: M,N,i,j,k,ij,err
        type(type_prec):: p

        M   =  p%bsize
        N   =  M*M
        if(.not. allocated(p%LU)) then
            allocate(p%LU(M,M*p%nvtx*p%nvtx), stat=err)
            if(err .ne. 0)  stop 'Error: fails to allocate LU.'
        end if

        p%LU=  0.0d0
        do i=1,p%nvtx
        do k=p%iA(i),p%iA(i+1)-1
            j   =  p%jA(k)
            ij  =  j+p%nvtx*(i-1)
            call DCOPY(N, p%A(1,1+M*(k-1)), 1, p%LU(1,1+M*(ij-1)), 1)
        end do
        end do
        call dgetrf(p%nvtx, M, p%LU)

        return
        end subroutine direct_decomposition
!       ------------------------------------------------------------------------
!       cal the LHS of coarse level.
!       ------------------------------------------------------------------------
        subroutine get_coarse_LHS(p0,p1)
        implicit none
        integer(dpI):: i,j,k,ic,jc,M
        type(type_prec):: p0,p1

        p1%A=  0.0d0
        M   =  p0%bsize**2
        do i=1,p0%nvtx
            ic  =  p0%ci(i)
            do j=p0%iA(i),p0%iA(i+1)-1
                jc  =  p0%ci(p0%jA(j))
                do k=p1%iA(ic),p1%iA(ic+1)-1
                    if(p1%jA(k) .ne. jc)    cycle
                    call DAXPY(M, 1.0d0, p0%A(1,1+p0%bsize*(j-1)), 1, &
                        &  p1%A(1,1+p1%bsize*(k-1)), 1)
                    exit
                end do
            end do
        end do

        return
        end subroutine get_coarse_LHS
!       ------------------------------------------------------------------------
!       solve Ax=RHS using the prec.
!       ------------------------------------------------------------------------
        subroutine prec_solve(p,is_zero_x,max_iter)
        implicit none
        logical(dpL),intent(in):: is_zero_x
        integer(dpI),intent(in):: max_iter
        integer(dpI):: mlev,i,j1,j0
        type(type_prec):: p

        if(p%lev .ne. 0)    stop 'Error: wrong input for prec_solve.'

!       RHS should already be set by the user.
        if(p%is_reordered) then
            do i=1,p%nvtx
                j1  =  1+p%bsize*(p%perm(i)-1)
                j0  =    p%bsize* p%perm(i)
                p%buf1(j1:j0)   =  p%RHS(1+p%bsize*(i-1):p%bsize*i)
            end do
            call DCOPY(p%bsize*p%nvtx, p%buf1, 1, p%RHS, 1)

            if(.not. is_zero_x) then
                do i=1,p%nvtx
                    j1  =  1+p%bsize*(p%perm(i)-1)
                    j0  =    p%bsize* p%perm(i)
                    p%buf1(j1:j0)   =  p%x(1+p%bsize*(i-1):p%bsize*i)
                end do
                call DCOPY(p%bsize*p%nvtx, p%buf1, 1, p%x, 1)
            end if
        end if

        mlev=  p%mlev
!       ------------------------------------------------------------------------
!       initialize x of all the coarser levels.
        if(mlev .gt. 1) call ini(p%clev)
!       initialize x of all the coarser levels.
!       ------------------------------------------------------------------------

        call prec_cycle(p, p%buf1, p%buf2, p%buf3)

!       the user should copy p%x to x.
        if(p%is_reordered) then
            do i=1,p%nvtx
                j1  =  1+p%bsize*(p%perm(i)-1)
                j0  =    p%bsize* p%perm(i)
                p%buf1(1+p%bsize*(i-1):p%bsize*i)   =  p%x(j1:j0)
            end do
            call DCOPY(p%bsize*p%nvtx, p%buf1, 1, p%x, 1)
        end if

        return
        contains
!           --------------------------------------------------------------------
!           ini for amg cycle.
!           --------------------------------------------------------------------
            recursive subroutine ini(p)
            implicit none
            type(type_prec):: p

            p%x =  0.0d0
            if(associated(p%clev))  call ini(p%clev)

            return
            end subroutine ini
!           --------------------------------------------------------------------
!           amg cycle.
!           --------------------------------------------------------------------
            recursive subroutine prec_cycle(p,b1,b2,b3)
            implicit none
            real(dpR):: b1(*),b2(*),b3(*)
            type(type_prec):: p

            if(p%solver .eq. slv_SGS) then
                call prec_iteration_damped(p,max_iter,p%RHS,p%x,b1,b2,b3)
            elseif(p%solver .eq. slv_SGS_MF) then
                if(p%lev .gt. 0)    stop 'Error: slv_SGS_MF not work with coarse level.'
            elseif(p%solver .eq. slv_ILU0) then
                call prec_iteration_damped(p,max_iter,p%RHS,p%x,b1,b2,b3)
            elseif(p%solver .eq. slv_LINE) then
                call prec_iteration_damped(p,max_iter,p%RHS,p%x,b1,b2,b3)
            elseif(p%solver .eq. slv_DIRECT) then
                call DCOPY(p%LDA, p%RHS, 1, p%x, 1)
                call dgetrs(p%nvtx, p%bsize, p%LU, p%x)
                return
            else
                stop 'Error: solver not supported, PREC.'
            end if
            if(p%lev .eq. mlev-1)   return

!           restriction of the residual to coarser level.
            call restriction(p, p%RHS, p%x, p%clev%RHS)

!           iteration on the next level.
            call prec_cycle(p%clev, b1, b2, b3)

!           prolongation of the increment from the coarser level.
            if(abs(p%relx_min_mg-1.0d0) .le. 1.0d-14) then
                call prolongation       (p, p%x)
            else
                call prolongation_damped(p, p%x)
            end if

            return
            end subroutine prec_cycle
!           --------------------------------------------------------------------
!           amg restriction of RHS to the coarser level.
!           --------------------------------------------------------------------
            subroutine restriction(p,RHS,x,RHS_c)
            implicit none
            type(type_prec),intent(in):: p
            real      (dpR),intent(in):: RHS(*),x(*)
            integer   (dpI):: i,j,k,ii,M
            real      (dpR):: RHS_c(*)

            M   =  p%bsize
            RHS_c(1:p%clev%LDA) =  0.0d0
            if((p%lev .eq. 0) .and. (p%solver .eq. slv_SGS_MF)) then
!               RHS=b-Ax should be done by the USER.
                do i=1,p%nvtx
                    ii  =  p%ci(i)
                    RHS_c(1+M*(ii-1):M*ii)  =  RHS_c(1+M*(ii-1):M*ii)+RHS(1+M*(i-1):M*i)
                end do
            else
                do i=1,p%nvtx
                    ii  =  p%ci(i)
                    RHS_c(1+M*(ii-1):M*ii)  =  RHS_c(1+M*(ii-1):M*ii)+RHS(1+M*(i-1):M*i)
                    do k=p%iA(i),p%iA(i+1)-1
                        j   =  p%jA(k)
                        call DGEMV('N', M, M,-1.0d0, p%A(1,1+M*(k-1)), M, x(1+M*(j-1)), &
                            &  1, 1.0d0, RHS_c(1+M*(ii-1)), 1)
                    end do
                end do
            end if

            return
            end subroutine restriction
!           --------------------------------------------------------------------
!           amg prolongation of x from the coarser level.
!           --------------------------------------------------------------------
            subroutine prolongation(p,x)
            implicit none
            type(type_prec),intent(in):: p
            real      (dpR):: x(*)
            integer   (dpI):: i,ii,M

            M   =  p%bsize
            do i=1,p%nvtx
                ii  =  p%ci(i)
                x(1+M*(i-1):M*i)=  x(1+M*(i-1):M*i)+p%clev%x(1+M*(ii-1):M*ii)
            end do

            return
            end subroutine prolongation
!           --------------------------------------------------------------------
!           amg prolongation of x from the coarser level, with damping.
!           --------------------------------------------------------------------
            subroutine prolongation_damped(p,x)
            implicit none
            type(type_prec),intent(in):: p
            integer   (dpI):: i,j,k,ii,jc,M
            real      (dpR):: x(*),a,b,relx,r(100),Adx(100)
            real      (dpR),external:: DDOT,DNRM2

            M   =  p%bsize
            a   =  0.0d0
            b   =  0.0d0
            do i=1,p%nvtx
                r  (1:M)=  p%RHS(1+M*(i-1):M*i)
                Adx(1:M)=  0.0d0
                do k=p%iA(i),p%iA(i+1)-1
                    j   =  p%jA(k)
                    jc  =  p%ci(j)
                    call DGEMV('N', M, M,-1.0d0, p%A(1,1+M*(k-1)), M, x(1+M*(j-1)), &
                        &  1, 1.0d0, r, 1)
                    call DGEMV('N', M, M, 1.0d0, p%A(1,1+M*(k-1)), M, &
                        &  p%clev%x(1+M*(jc-1)), 1, 1.0d0, Adx, 1)
                end do
                a   =  a+DDOT (M, Adx, 1, r, 1)
                b   =  b+DNRM2(M, Adx, 1)**2
            end do
            relx=  min(max(a/b, p%relx_min_mg), 1.0d0)
            do i=1,p%nvtx
                ii  =  p%ci(i)
                x(1+M*(i-1):M*i)=  x(1+M*(i-1):M*i)+relx*p%clev%x(1+M*(ii-1):M*ii)
            end do

            return
            end subroutine prolongation_damped
        end subroutine prec_solve
!       ------------------------------------------------------------------------
!       iteration, with damping.
!       ------------------------------------------------------------------------
        subroutine prec_iteration_damped(p,max_iter,RHS,x,r,dx,Adx)
        implicit none
        real   (dpR),parameter:: eps=2.0d-1
        integer(dpI),intent(in):: max_iter
        real   (dpR),intent(in):: RHS(*)
        integer(dpI):: i,j,k,M,iter
        real   (dpR):: x(*),r(*),dx(*),Adx(*),v(100),residual(0:100),a,b,relx
        type(type_prec):: p
        real   (dpR),external:: DDOT,DNRM2

        if(p%solver .eq. slv_SGS) then
            if(.not. p%is_damped) then
                call gs_iteration(p%nvtx, p%iA, p%jA, p%bsize, max_iter, p%A, p%D, RHS, x)
            elseif(p%bsize .eq. 1) then
                call prec_gs_damped (p, max_iter, RHS, x, r, dx, Adx)
            else
                call prec_bgs_damped(p, max_iter, RHS, x, r, dx, Adx)
            end if
            return
        elseif(p%solver .eq. slv_ILU0) then
            if((max_iter .eq. 1) .and. (.not. p%is_damped)) then
                if(p%is_LU) then
                    call ilu0_solution(p%nvtx, p%iA, p%jA, p%bsize, p%LU, RHS, x)
                else
                    call ilu0_solution(p%nvtx, p%iA, p%jA, p%bsize, p%A , RHS, x)
                end if
                return
            end if
        end if

        M   =  p%bsize
!       r   =  b-A*x
        call DCOPY(p%LDA, RHS, 1, r, 1)
        do i=1,p%nvtx
        do k=p%iA(i),p%iA(i+1)-1
            j   =  p%jA(k)
            call DGEMV('N', M, M,-1.0d0, p%A(1,1+M*(k-1)), M, x(1+M*(j-1)), &
                &  1, 1.0d0, r(1+M*(i-1)), 1)
        end do
        end do
        residual(0) =  DNRM2(p%LDA, r, 1)

        do iter=1,max_iter
!           step 1: dx=A^(-1)r
            dx(1:p%LDA) =  0.0d0
            if(p%solver .eq. slv_SGS) then
                stop 'Error: prec_iteration_damped fails.'
            elseif(p%solver .eq. slv_ILU0) then
                if(p%is_LU) then
                    call ilu0_solution(p%nvtx, p%iA, p%jA, p%bsize, p%LU, r, dx)
                else
                    call ilu0_solution(p%nvtx, p%iA, p%jA, p%bsize, p%A , r, dx)
                end if
            elseif(p%solver .eq. slv_LINE) then
                call line_iteration(p, r, dx)
            end if

!           step 2: get damp coefficient.
            a   =  0.0d0
            b   =  0.0d0
            do i=1,p%nvtx
                v(1:M)  =  0.0d0
                do k=p%iA(i),p%iA(i+1)-1
                    j   =  p%jA(k)
                    call DGEMV('N', M, M, 1.0d0, p%A(1,1+M*(k-1)), M, dx(1+M*(j-1)), &
                        &  1, 1.0d0, v, 1)
                end do
                a   =  a+DDOT (M, v, 1, r(1+M*(i-1)), 1)
                b   =  b+DNRM2(M, v, 1)**2
                Adx(1+M*(i-1):M*i)  =  v(1:M)
            end do
            relx=  min(max(a/b, p%relx_min), p%relx_max)
            x(1:p%LDA)  =  x(1:p%LDA)+relx*dx(1:p%LDA)
            call DAXPY(p%LDA, -relx, Adx, 1, r, 1)
            residual(iter)  =  DNRM2(p%LDA, r, 1)
!           print*,iter,relx,residual(iter)/residual(0)
            if(residual(iter) .le. eps*residual(0)) exit
        end do

        return
        end subroutine prec_iteration_damped
!       ------------------------------------------------------------------------
!       GS iteration, with damping.
!       ------------------------------------------------------------------------
        subroutine prec_gs_damped(p,max_iter,RHS,x,r,dx,Adx)
        implicit none
        real   (dpR),parameter:: eps=2.0d-1
        integer(dpI),intent(in):: max_iter
        real   (dpR),intent(in):: RHS(*)
        integer(dpI):: i,j,k,iter
        real   (dpR):: x(*),r(*),dx(*),Adx(*),v,w,residual(0:100),a,b,relx
        type(type_prec):: p
        real   (dpR),external:: DDOT,DNRM2

!       r   =  b-A*x
        call DCOPY(p%LDA, RHS, 1, r, 1)
        do i=1,p%nvtx
        do k=p%iA(i),p%iA(i+1)-1
            j   =  p%jA(k)
            r(i)=  r(i)-p%A(1,k)*x(j)
        end do
        end do
        residual(0) =  DNRM2(p%LDA, r, 1)

        do iter=1,max_iter
!           step 1: forward iteration.
            do i=1,p%nvtx
                v   =  r(i)
                do k=p%iA(i),p%iA(i+1)-1
                    j   =  p%jA(k)
                    if(j .ge. i)    cycle
                    v   =  v-p%A(1,k)*dx(j)
                end do
                dx(i)   =  p%D(1,i)*v
            end do

!           step 2: backward iteration.
            Adx(1:p%LDA)=  0.0d0
            do i=p%nvtx,1,-1
                v   =  r(i)
                do k=p%iA(i+1)-1,p%iA(i),-1
                    j   =  p%jA(k)
                    w   =  p%A(1,k)*dx(j)
                    if(i .ne. j)    v       =  v-w
                    if(i .lt. j)    Adx(i)  =  Adx(i)+w
                end do
                dx(i)   =  p%D(1,i)*v
            end do

!           step 3: A*dx.
            do i=1,p%nvtx
                do k=p%iA(i),p%iA(i+1)-1
                    j   =  p%jA(k)
                    if(i .lt. j)    cycle
                    Adx(i)  =  Adx(i)+p%A(1,k)*dx(j)
                end do
            end do

!           step 4: get damping coefficient.
            a   =  DDOT (p%LDA, Adx, 1, r, 1)
            b   =  DNRM2(p%LDA, Adx, 1)**2
            relx=  min(max(a/b, p%relx_min), p%relx_max)

            x(1:p%LDA)  =  x(1:p%LDA)+relx*dx(1:p%LDA)
            call DAXPY(p%LDA, -relx, Adx, 1, r, 1)
            residual(iter)  =  DNRM2(p%LDA, r, 1)
!           print*,iter,residual(iter)/residual(0)
            if(residual(iter) .le. eps*residual(0)) exit
        end do

        return
        end subroutine prec_gs_damped
!       ------------------------------------------------------------------------
!       Block-GS iteration, with damping.
!       ------------------------------------------------------------------------
        subroutine prec_bgs_damped(p,max_iter,RHS,x,r,dx,Adx)
        implicit none
        real   (dpR),parameter:: eps=2.0d-1
        integer(dpI),intent(in):: max_iter
        real   (dpR),intent(in):: RHS(*)
        integer(dpI):: i,j,k,M,iter
        real   (dpR):: x(*),r(*),dx(*),Adx(*),v(100),w(100),residual(0:100),a,b,relx
        type(type_prec):: p
        real   (dpR),external:: DDOT,DNRM2

        M   =  p%bsize
!       r   =  b-A*x
        call DCOPY(p%LDA, RHS, 1, r, 1)
        do i=1,p%nvtx
        do k=p%iA(i),p%iA(i+1)-1
            j   =  p%jA(k)
            call DGEMV('N', M, M,-1.0d0, p%A(1,1+M*(k-1)), M, x(1+M*(j-1)), &
                &  1, 1.0d0, r(1+M*(i-1)), 1)
        end do
        end do
        residual(0) =  DNRM2(p%LDA, r, 1)

        do iter=1,max_iter
!           step 1: forward iteration.
            do i=1,p%nvtx
                v(1:M)  =  r(1+M*(i-1):M*i)
                do k=p%iA(i),p%iA(i+1)-1
                    j   =  p%jA(k)
                    if(j .ge. i)    cycle
                    call DGEMV('N', M, M,-1.0d0, p%A(1,1+M*(k-1)), M, dx(1+M*(j-1)), &
                        &  1, 1.0d0, v, 1)
                end do
                call DGEMV('N', M, M, 1.0d0, p%D(1,1+M*(i-1)), M, v, 1, 0.0d0, &
                    &  dx(1+M*(i-1)), 1)
            end do

!           step 2: backward iteration.
            Adx(1:p%LDA)=  0.0d0
            do i=p%nvtx,1,-1
                v(1:M)  =  r(1+M*(i-1):M*i)
                do k=p%iA(i+1)-1,p%iA(i),-1
                    j   =  p%jA(k)
                    call DGEMV('N', M, M, 1.0d0, p%A(1,1+M*(k-1)), M, dx(1+M*(j-1)), &
                        &  1, 0.0d0, w, 1)
                    if(i .ne. j)    v(1:M)  =  v(1:M)-w(1:M)
                    if(i .lt. j)    Adx(1+M*(i-1):M*i)  =  Adx(1+M*(i-1):M*i)+w(1:M)
                end do
                call DGEMV('N', M, M, 1.0d0, p%D(1,1+M*(i-1)), M, v, 1, 0.0d0, &
                    &  dx(1+M*(i-1)), 1)
            end do

!           step 3: A*dx.
            do i=1,p%nvtx
                do k=p%iA(i),p%iA(i+1)-1
                    j   =  p%jA(k)
                    if(i .lt. j)    cycle
                    call DGEMV('N', M, M, 1.0d0, p%A(1,1+M*(k-1)), M, dx(1+M*(j-1)), &
                        &  1, 1.0d0, Adx(1+M*(i-1)), 1)
                end do
            end do

!           step 4: get damping coefficient.
            a   =  DDOT (p%LDA, Adx, 1, r, 1)+2.0d0*tiny(1.0d0)
            b   =  DNRM2(p%LDA, Adx, 1)**2   +2.0d0*tiny(1.0d0)
            relx=  min(max(a/b, p%relx_min), p%relx_max)

            x(1:p%LDA)  =  x(1:p%LDA)+relx*dx(1:p%LDA)
            call DAXPY(p%LDA, -relx, Adx, 1, r, 1)
            residual(iter)  =  DNRM2(p%LDA, r, 1)
!           print*,iter,relx,residual(iter)/residual(0)
            if(residual(iter) .le. eps*residual(0)) exit
        end do

        return
        end subroutine prec_bgs_damped
!       ------------------------------------------------------------------------
!       line implicit iteration.
!       ------------------------------------------------------------------------
        subroutine line_iteration(p,RHS,x)
        implicit none
        real   (dpR),intent(in):: RHS(*)
        integer(dpI):: i,j,k,M,N,L,D,U,L1,L0,S
        real   (dpR):: x(*),v(100)
        type(type_prec):: p

        M   =  p%bsize

!       ------------------------------------------------------------------------
!       forward sweep: line solver and then point Gauss-Seidel.
        do i=1,p%n_line
            L1  =  p%iA_line(i)
            L0  =  p%iA_line(i+1)-1
            do j=L1,L0
!               the lower neighbour.
                if(j .eq. L1) then
                    L   =  0
                else
                    L   =  p%jA_line(j-1)
                end if

                D   =  p%jA_line(j)

!               the upper neighbour.
                if(j .eq. L0) then
                    U   =  0
                else
                    U   =  p%jA_line(j+1)
                end if

!               b-A*x'.
                p%R_line(1:M,j) =  RHS(1+M*(D-1):M*D)
                do k=p%iA(D),p%iA(D+1)-1
                    n   =  p%jA(k)
                    if((n .eq. L) .or. (n .eq. D) .or. (n .eq. U))  cycle
                    call DGEMV('N', M, M,-1.0d0, p%A(1,1+M*(k-1)), M, x(1+M*(n-1)), &
                        &  1, 1.0d0, p%R_line(1,j), 1)
                end do
            end do

            j   =  1+M*(L1-1)
            call B_dgttrs(L0-L1+1, p%L_line(1,j), p%D_line(1,j), p%U_line(1,j), &
                &  p%R_line(1,L1))
            do j=L1,L0
                D   =  p%jA_line(j)
                x(1+M*(D-1):M*D)=  p%R_line(1:M,j)
            end do
        end do

        do S=p%n_line+1,p%n_vtx_line
            i       =  p%jA_line(p%iA_line(S))
            v(1:M)  =  RHS(1+M*(i-1):M*i)
            do k=p%iA(i),p%iA(i+1)-1
                j   =  p%jA(k)
                if(i .eq. j)    cycle
                call DGEMV('N', M, M,-1.0d0, p%A(1,1+M*(k-1)), M, x(1+M*(j-1)), &
                    &  1, 1.0d0, v, 1)
            end do
            call DGEMV('N',M,M,1.0d0,p%D(1,1+M*(i-1)),M,v,1,0.0d0,x(1+M*(i-1)),1)
        end do
!       forward sweep: line solver and then point Gauss-Seidel.
!       ------------------------------------------------------------------------

!       ------------------------------------------------------------------------
!       backward sweep: point Gauss-Seidel and then line solver.
        do S=p%n_line+1,p%n_vtx_line
            i       =  p%jA_line(p%iA_line(S))
            v(1:M)  =  RHS(1+M*(i-1):M*i)
            do k=p%iA(i),p%iA(i+1)-1
                j   =  p%jA(k)
                if(i .eq. j)    cycle
                call DGEMV('N', M, M,-1.0d0, p%A(1,1+M*(k-1)), M, x(1+M*(j-1)), &
                    &  1, 1.0d0, v, 1)
            end do
            call DGEMV('N',M,M,1.0d0,p%D(1,1+M*(i-1)),M,v,1,0.0d0,x(1+M*(i-1)),1)
        end do

        do i=p%n_line,1,-1
            L1  =  p%iA_line(i)
            L0  =  p%iA_line(i+1)-1
            do j=L1,L0
!               the lower neighbour.
                if(j .eq. L1) then
                    L   =  0
                else
                    L   =  p%jA_line(j-1)
                end if

                D   =  p%jA_line(j)

!               the upper neighbour.
                if(j .eq. L0) then
                    U   =  0
                else
                    U   =  p%jA_line(j+1)
                end if

!               b-A*x'.
                p%R_line(1:M,j) =  RHS(1+M*(D-1):M*D)
                do k=p%iA(D),p%iA(D+1)-1
                    n   =  p%jA(k)
                    if((n .eq. L) .or. (n .eq. D) .or. (n .eq. U))  cycle
                    call DGEMV('N', M, M,-1.0d0, p%A(1,1+M*(k-1)), M, x(1+M*(n-1)), &
                        &  1, 1.0d0, p%R_line(1,j), 1)
                end do
            end do

            j   =  1+M*(L1-1)
            call B_dgttrs(L0-L1+1, p%L_line(1,j), p%D_line(1,j), p%U_line(1,j), &
                &  p%R_line(1,L1))
            do j=L1,L0
                D   =  p%jA_line(j)
                x(1+M*(D-1):M*D)=  p%R_line(1:M,j)
            end do
        end do
!       backward sweep: point Gauss-Seidel and then line solver.
!       ------------------------------------------------------------------------

        return
        end subroutine line_iteration
!       ------------------------------------------------------------------------
!       delete the prec.
!       ------------------------------------------------------------------------
        subroutine prec_delete(p)
        implicit none
        type(type_prec):: p

        if(allocated(p%iA  ))   deallocate(p%iA)
        if(allocated(p%jA  ))   deallocate(p%jA)
        if(allocated(p%perm))   deallocate(p%perm)
        if(allocated(p%A   ))   deallocate(p%A)
        if(allocated(p%D   ))   deallocate(p%D)
        if(allocated(p%LU  ))   deallocate(p%LU)
        if(allocated(p%RHS ))   deallocate(p%RHS)
        if(allocated(p%x   ))   deallocate(p%x)
        if(allocated(p%ci  ))   deallocate(p%ci)
        p%LDA   =  0
        p%bsize =  0

        return
        end subroutine prec_delete
!       ------------------------------------------------------------------------
!       set the implicit lines.
!       ------------------------------------------------------------------------
        subroutine prec_set_lines(p,n_line,iA_line,jA_line)
        implicit none
        integer(dpI),intent(in):: n_line,iA_line(*),jA_line(*)
        logical(dpL):: ltmp
        integer(dpI):: i,j,k,N,L,R,jj
        type(type_prec):: p
        integer(dpI),allocatable:: ID(:)

        if(p%lev .ne. 0)    stop 'Error: wrong input for set_line.'
        if(.not. p%is_spy)  stop 'Error: SPY of prec is not available.'
        if(n_line .le. 0)   return
        allocate(ID(p%nvtx))
        ID  =  0

        if(iA_line(1) .ne. 1)   stop 'Error: wrong input for set_line.'

        allocate(p%ID_LU_line(3,iA_line(n_line+1)-1), stat=err_mem)
        p%ID_LU_line=  0
        do i=1,n_line
        do j=iA_line(i),iA_line(i+1)-1
            k   =  jA_line(j)
            if((k .le. 0) .or. (k .gt. p%nvtx)) stop 'Error: wrong input for set_line.'
            if(ID(k) .gt. 0)    stop 'Error: wrong input for set_line.'

!           check the lower element.
            if(j .gt. iA_line(i)) then
                L   =  jA_line(j-1)
                ltmp=  .false.
                do jj=p%iA(k),p%iA(k+1)-1
                    ltmp=  p%jA(jj) .eq. L
                    if(ltmp)    p%ID_LU_line(1,j)   =  jj
                    if(ltmp)    exit
                end do
                if(.not. ltmp)  stop 'Error: wrong input for set_line.'
            end if

            do jj=p%iA(k),p%iA(k+1)-1
                if(p%jA(jj) .eq. k) then
                    p%ID_LU_line(2,j)   =  jj
                    exit
                end if
            end do

!           check the upper element.
            if(j .lt. iA_line(i+1)-1) then
                R   =  jA_line(j+1)
                ltmp=  .false.
                do jj=p%iA(k),p%iA(k+1)-1
                    ltmp=  p%jA(jj) .eq. R
                    if(ltmp)    p%ID_LU_line(3,j)   =  jj
                    if(ltmp)    exit
                end do
                if(.not. ltmp)  stop 'Error: wrong input for set_line.'
            end if

            ID(k)   =  i
        end do
        end do

        p%n_vtx_line=  n_line
        p%n_line    =  n_line
        do i=1,p%nvtx
            if(ID(i) .le. 0)    p%n_vtx_line=  p%n_vtx_line+1
        end do
        allocate(p%iA_line(p%n_vtx_line+1), stat=err_mem)
        allocate(p%jA_line(p%nvtx        ), stat=err_mem)

        N   =  iA_line(n_line+1)-1
        call ICOPY(n_line+1, iA_line, 1, p%iA_line, 1)
        call ICOPY(N       , jA_line, 1, p%jA_line, 1)

        j   =  n_line
        k   =  N
        do i=1,p%nvtx
            if(ID(i) .gt. 0)    cycle
            j   =  j+1
            k   =  k+1
            p%iA_line(j+1)  =  p%iA_line(j)+1
            p%jA_line(k  )  =  i
        end do
        if((j .ne. p%n_vtx_line) .or. (k .ne. p%nvtx))  stop 'Error: fails to set line.'
        allocate(p%L_line(p%bsize, p%bsize*N), stat=err_mem)
        allocate(p%D_line(p%bsize, p%bsize*N), stat=err_mem)
        allocate(p%U_line(p%bsize, p%bsize*N), stat=err_mem)
        allocate(p%R_line(p%bsize,         N), stat=err_mem)
        p%L_line=  0.0d0
        p%D_line=  0.0d0
        p%U_line=  0.0d0
        p%R_line=  0.0d0
        if(p%solver .eq. slv_SGS)   p%solver=  slv_LINE

        if(allocated(ID))   deallocate(ID)

        return
        end subroutine prec_set_lines
!       ------------------------------------------------------------------------
!       set the LHS of lines.
!       ------------------------------------------------------------------------
        subroutine prec_set_line_LHS(p)
        implicit none
        integer(dpI):: i,iele,idx,L,D,U,M,L1,L0
        type(type_prec):: p

        if(p%n_line .le. 0) return
        if(.not. allocated(p%iA_line))  stop 'Error: iA_line not available.'

        M   =  p%bsize**2
        do i=1,p%n_line
            L1  =  p%iA_line(i)
            L0  =  p%iA_line(i+1)-1
            do iele=L1,L0
                idx =  1+p%bsize*(iele-1)

                if(iele .gt. p%iA_line(i)) then
                    L   =  p%ID_LU_line(1,iele)
                    call DCOPY(M, p%A(1,1+p%bsize*(L-1)), 1, p%L_line(1,idx), 1)
                end if

                D   =  p%ID_LU_line(2,iele)
                call DCOPY(M, p%A(1,1+p%bsize*(D-1)), 1, p%D_line(1,idx), 1)

                if(iele .lt. p%iA_line(i+1)-1) then
                    U   =  p%ID_LU_line(3,iele)
                    call DCOPY(M, p%A(1,1+p%bsize*(U-1)), 1, p%U_line(1,idx), 1)
                end if
            end do

            idx =  1+p%bsize*(L1-1)
            call B_dgttrf(L0-L1+1, p%L_line(1,idx), p%D_line(1,idx), p%U_line(1,idx))
        end do

        return
        end subroutine prec_set_line_LHS
    end module var_prec

    module var_lhs_imp
    use var_kind_def
    implicit none

    integer(dpI):: n_ele_i,n_face_i
    integer(dpI),allocatable:: owner(:),neighbor(:)
    integer(dpI),allocatable:: owner_int(:),neighbor_int(:)
    real   (dpR),allocatable:: upper(:,:),lower(:,:),diag(:,:),x(:,:),b(:,:)
    real   (dpR),allocatable:: upper_int(:,:),lower_int(:,:),diag_int(:,:)
    real   (dpR),allocatable:: x_int(:,:),b_int(:,:)

    end module var_lhs_imp