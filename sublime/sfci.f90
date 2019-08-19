!-------------------------------------------------------------------------------
!   add elelemts to list.
!-------------------------------------------------------------------------------
    subroutine add_ele_to_list(nele,ele,nlst,lst)
    use var_kind_def
    implicit none
    logical(dpL):: ltmp
    integer(dpI):: Nlst,Nele,lst(*),ele(*),i,j

    if(nele .le. 0) return
    do i=1,nele
        ltmp=  .false.
        do j=1,nlst
            if(lst(j) .eq. ele(i)) then
                ltmp=  .true.
                exit
            end if
        end do
        if(ltmp)    cycle
        nlst        =  nlst+1
        lst(nlst)   =  ele(i)
    end do

    return
    end subroutine add_ele_to_list
!-------------------------------------------------------------------------------
!   find a column in the matrix.
!-------------------------------------------------------------------------------
    subroutine find_col_in_matrix(LDA,R1,R0,c1,c0,matrix,col,nele,ele)
    use var_kind_def
    implicit none
    integer(dpI),intent(in):: LDA,R1,R0,c1,c0,matrix(LDA,*),col(*)
    logical(dpL):: ltmp
    integer(dpI):: nele,ele(*),c(LDA),L,R,i,j

    nele    =  0
    c       =  0
    c(R1:R0)=  col(R1:R0)
    call iqsort(.true., R1, R0, c)
    call ib_search(c1, c0, R1, LDA, matrix, c(1), ltmp, L, R)
    if(.not. ltmp)  return
    do j=L,R
        ltmp=  .false.
        do i=R1+1,R0
            ltmp=  c(i-R1+1) .ne. matrix(i,j)
            if(ltmp)    exit
        end do
        if(ltmp)    cycle
        nele        =  nele+1
        ele(nele)   =  j
    end do

    return
    end subroutine find_col_in_matrix
!-------------------------------------------------------------------------------
!   copies a vector, x, to a vector, y, integer.
!-------------------------------------------------------------------------------
    subroutine ICOPY(N,dx,incx,dy,incy)
    use var_kind_def
    implicit none
    integer(dpI),intent(in):: N,dx(*),incx,incy
    integer(dpI):: dy(*),i,M,ix,iy

    if(N .le. 0)    return
    if((incx .eq. 1) .and. (incy .eq. 1)) then
        M   =  mod(N, 7)
        do i=1,M
            dy(i)   =  dx(i)
        end do
        if(N .lt. 7)    return
        do i=M+1,N,7
            dy(i)   =  dx(i)
            dy(i+1) =  dx(i+1)
            dy(i+2) =  dx(i+2)
            dy(i+3) =  dx(i+3)
            dy(i+4) =  dx(i+4)
            dy(i+5) =  dx(i+5)
            dy(i+6) =  dx(i+6)
        end do
    else
        ix  =  1
        iy  =  1
        if(incx .lt. 0) ix  =  incx*(-N+1)+1
        if(incy .lt. 0) iy  =  incy*(-N+1)+1
        do i=1,N
            dy(iy)  =  dx(ix)
            ix      =  ix+incx
            iy      =  iy+incy
        end do
    end if

    return
    end subroutine ICOPY
!-------------------------------------------------------------------------------
!   compute the tanspose of CSR matrix, integer, only the sparse pattern.
!-------------------------------------------------------------------------------
    subroutine itrans_spy(is_constant_num,LDA,dimA,iA,jA,dimT,iT,jT)
    use var_kind_def
    implicit none
    logical(dpL),intent(in):: is_constant_num
    integer(dpI),intent(in):: LDA,dimA,iA(*),jA(*),dimT
    integer(dpI):: iT(*),jT(*),i,j,k,M,N

    jT(1:dimT)  =  0
    if(is_constant_num) then
        do i=1,dimA
        do k=1+(i-1)*LDA,i*LDA
!           A(i,j) is an element of A.
            j       =  jA(k)
            jT(j)   =  jT(j)+1
        end do
        end do
    else
        do i=1,dimA
        do k=iA(i),iA(i+1)-1
!           A(i,j) is an element of A.
            j       =  jA(k)
            jT(j)   =  jT(j)+1
        end do
        end do
    end if

    iT(1)   =  1
    do i=1,dimT
        iT(i+1) =  iT(i)+jT(i)
    end do
    do i=1,dimT
        jT(iT(i+1)-1)   =  iT(i)
    end do

    if(is_constant_num) then
        do i=1,dimA
        do k=1+(i-1)*LDA,i*LDA
            j   =  jA(k)

            M   =  iT(j+1)-1
            N   =  jT(M)
            if(N .gt. M)    stop 'Error: ITRANS_SPY fails.'
            jT(N)   =  i
            if(N .lt. M)    jT(M)   =  N+1
        end do
        end do
    else
        do i=1,dimA
        do k=iA(i),iA(i+1)-1
            j   =  jA(k)

            M   =  iT(j+1)-1
            N   =  jT(M)
            if(N .gt. M)    stop 'Error: ITRANS_SPY fails.'
            jT(N)   =  i
            if(N .lt. M)    jT(M)   =  N+1
        end do
        end do
    end if

    return
    end subroutine itrans_spy
!-------------------------------------------------------------------------------
!   sort and simplify a series of numbers.
!-------------------------------------------------------------------------------
    subroutine simplify_series(N,col,LDA,v)
    use var_kind_def
    implicit none
    integer(dpI),intent(in):: col,LDA
    integer(dpI):: N,v(*),N0,i

    if(N .le. 0)    return
    N0  =  N
    N   =  1
    if(LDA .eq. 1) then
        call iqsort    (.true., 1, N0,          v)
        do i=2,N0
            if(v(i) .eq. v(i-1))    cycle
            N   =  N+1
            v(N)=  v(i)
        end do
    else
        call iqsortcols(.true., 1, N0, col, LDA, v)
        do i=2,N0
            if(v(col+LDA*(i-1)) .eq. v(col+LDA*(i-2)))  cycle
            N   =  N+1
            v(1+LDA*(N-1):LDA*N)=  v(1+LDA*(i-1):LDA*i)
        end do
    end if

    return
    end subroutine simplify_series
!-------------------------------------------------------------------------------
!   sort and simplify a list, also return the frequency.
!-------------------------------------------------------------------------------
    subroutine simplify_list_frequency(N,col,LDA,v,freq)
    use var_kind_def
    implicit none
    integer(dpI),intent(in):: col,LDA
    integer(dpI):: N,v(*),freq(*),N0,i,j,k

    if(N .le. 0)    return
    call iqsortcols(.true., 1, N, col, LDA, v)
    N0  =  N
    i   =  1
    N   =  0
    do while(i .le. N0)
        k   =  i
        do j=i+1,N0
            if(v(col+LDA*(j-1)) .eq. v(col+LDA*(i-1))) then
                k   =  j
            else
                exit
            end if
        end do
        N       =  N+1
        freq(N) =  k-i+1
        v(1+LDA*(N-1):LDA*N)=  v(1+LDA*(i-1):LDA*i)
        i       =  k+1
    end do

    return
    end subroutine simplify_list_frequency
!-------------------------------------------------------------------------------
!   subtract a list from another list.
!-------------------------------------------------------------------------------
    subroutine list_minus_list(n_a,a,n_b,b)
    use var_kind_def
    implicit none
    integer(dpI),intent(in):: n_b,b(*)
    logical(dpL):: ltmv(n_a),ltmp
    integer(dpI):: n_a,a(*),i,j,L,R

    call simplify_series(n_a, 1, 1, a)
    ltmv=  .true.
    do j=1,n_b
        call ib_search(1, n_a, 1, 1, a, b(j), ltmp, L, R)
        if(ltmp)    ltmv(L:R)   =  .false.
    end do
    j   =  n_a
    n_a =  0
    do i=1,j
        if(ltmv(i)) then
            n_a     =  n_a+1
            a(n_a)  =  a(i)
        end if
    end do

    return
    end subroutine list_minus_list
!-------------------------------------------------------------------------------
!   sort and simplify a matrix, column-by-column.
!-------------------------------------------------------------------------------
    subroutine simplify_matrix(LDA,N,r1,r0,A)
    use var_kind_def
    implicit none
    integer(dpI),intent(in):: LDA,r1,r0
    logical(dpL):: ii_oo,ltmp
    integer(dpI):: N,A(*),i,j,k,ii,oo,NA,r,T

    call iqsortcols(.true., 1, N, r1, LDA, A)
    i   =  1
    NA  =  0
    do while(i .le. N)
        k   =  i
        do j=i+1,N
            if(A(r1+LDA*(j-1)) .eq. A(r1+LDA*(i-1))) then
                k   =  j
            else
                exit
            end if
        end do

        NA  =  NA+1
        T   =  NA
        if(NA .ne. i)   call ICOPY(LDA, A(1+LDA*(i-1)), 1, A(1+LDA*(NA-1)), 1)

        do oo=i+1,k
            ltmp=  .false.
            do ii=T,NA
                ii_oo   =  .true.
                do r=r1+1,r0
                    if(A(r+LDA*(ii-1)) .ne. A(r+LDA*(oo-1))) then
                        ii_oo   =  .false.
                        exit
                    end if
                end do
                if(ii_oo) then
                    ltmp=  .true.
                    exit
                end if
            end do
            if(ltmp)    cycle
            NA  =  NA+1
            if(NA .ne. oo)  call ICOPY(LDA, A(1+LDA*(oo-1)), 1, A(1+LDA*(NA-1)), 1)
        end do

        i   =  k+1
    end do
    N   =  NA

    return
    end subroutine simplify_matrix
!-------------------------------------------------------------------------------
!   flip an integer matrix, left-and-right.
!-------------------------------------------------------------------------------
    subroutine ifliplr(LDA,L,R,M)
    use var_kind_def
    implicit none
    integer(dpI),intent(in):: LDA,L,R
    integer(dpI):: i,j
    integer(dpI):: M(LDA,*),v(LDA)

    do i=L,R
        j   =  R+1-i
        if(i .ge. j)    return
        v(1:LDA)    =  M(1:LDA,i)
        M(1:LDA,i)  =  M(1:LDA,j)
        M(1:LDA,j)  =  v(1:LDA)
    end do

    return
    end subroutine ifliplr
!-------------------------------------------------------------------------------
!   whether list A equals to list B.
!-------------------------------------------------------------------------------
    function list_eq_list(N,a,b) result(is_equal)
    use var_kind_def
    implicit none
    integer(dpI),intent(in):: N,a(*),b(*)
    logical(dpL):: is_equal
    integer(dpI):: i

    is_equal=  .false.
    do i=1,N
        if(a(i) .ne. b(i))  return
    end do
    is_equal=  .true.

    return
    end function list_eq_list
!-------------------------------------------------------------------------------
!   permute a list according to the given ordering.
!-------------------------------------------------------------------------------
    subroutine permute_list(LDA,N,order,L)
    use var_kind_def
    implicit none
    integer(dpI),intent(in):: LDA,N,order(*)
    logical(dpL):: ltmp
    integer(dpI):: L(LDA,*),p(N),v(LDA),w(LDA),id,id_seed,i,id_n

    p(1:N)  =  order(1:N)
    id_seed =  1
    do while(.true.)
!       get the seed.
        ltmp=  .true.
        do i=id_seed,N
            if(p(i) .gt. 0) then
                ltmp    =  .false.
                id_seed =  i
                exit
            end if
        end do
        if(ltmp)    exit

        id      =  id_seed
        v(1:LDA)=  L(1:LDA,id)
        do while(.true.)
            id_n    =  p(id)
            p(id)   = -abs(id_n)
            if(id_n .eq. id) then
                exit
            elseif(id_n .gt. 0) then
                w(1:LDA     )   =  L(1:LDA,id_n)
                L(1:LDA,id_n)   =  v(1:LDA)
                v               =  w
                id              =  id_n
            elseif(id_n .lt. 0) then
                L(1:LDA,-id_n)  =  v(1:LDA)
                exit
            else
                stop 'Error: permute_list fails, id_n=0.'
            end if
        end do
    end do

    return
    end subroutine permute_list
