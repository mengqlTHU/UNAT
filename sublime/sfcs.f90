!-------------------------------------------------------------------------------
!   binary search, integer.
!-------------------------------------------------------------------------------
    subroutine ib_search(ibeg,iend,col,LDA,A,x,find,L,R)
    use var_kind_def
    implicit none
    integer(dpI),intent(in):: ibeg,iend,col,LDA,A(*),x
    logical(dpL):: find
    integer(dpI):: M,L,R,P,d

    find=  .false.
    if(iend .lt. ibeg)  return

    L   =  ibeg-1
    R   =  iend+1
    do while(.true.)
        p   = (R-L)/2
        if(P .le. 0)    exit
        M   =  L+p
        d   =  x-A(col+LDA*(M-1))
        if(d .lt. 0) then
            R   =  M
        elseif(d .eq. 0) then
            find=  .true.
            exit
        else
            L   =  M
        end if
    end do
    if(.not. find)  return

    L   =  M
    R   =  M
    M   =  L-1
    do while(M .ge. ibeg)
        if(x .eq. A(col+LDA*(M-1))) then
            L   =  M
            M   =  M-1
        else
            exit
        end if
    end do
    M   =  R+1
    do while(M .le. iend)
        if(x .eq. A(col+LDA*(M-1))) then
            R   =  M
            M   =  M+1
        else
            exit
        end if
    end do

    return
    end subroutine ib_search
!-------------------------------------------------------------------------------
!   quicksort, integer.
!-------------------------------------------------------------------------------
    subroutine iqsort(asc,lef,rig,A)
    use var_kind_def
    implicit none
    integer(dpI),parameter:: mlen=9,mstk=30
    logical(dpL),intent(in):: asc
    integer(dpI),intent(in):: lef,rig
    integer(dpI):: A(*),L,R,M,i,j,nstk,pvt,Lstk(mstk),Rstk(mstk),itmp

    nstk=  0
    L   =  lef
    R   =  rig

    do
!       use insertion sort to sort A(L:R) and pop new sequence from stack.
        if(R-L+1 .le. mlen) then
            do j=L+1,R
                itmp=  A(j)
                do i=j-1,L,-1
                    if(A(i) .le. itmp)  exit
                    A(i+1)  =  A(i)
                end do
                A(i+1)  =  itmp
            end do

            if(nstk .le. 0) exit
            L   =  Lstk(nstk)
            R   =  Rstk(nstk)
            nstk=  nstk-1
        end if
        if(R .le. L) then
            if(nstk .le. 0) exit
            L   =  Lstk(nstk)
            R   =  Rstk(nstk)
            nstk=  nstk-1
        end if

        M   =  L+(R-L)/2
        if(A(L) .gt. A(M)) then
            itmp=  A(L)
            A(L)=  A(M)
            A(M)=  itmp
        end if
        if(A(M) .gt. A(R)) then
            itmp=  A(M)
            A(M)=  A(R)
            A(R)=  itmp
        end if
        if(A(L) .gt. A(M)) then
            itmp=  A(L)
            A(L)=  A(M)
            A(M)=  itmp
        end if
        pvt =  A(M)

!       now A(L)<=pvt<=A(R), we have to find the final location of pvt.
        i   =  L
        j   =  R
        do
            do
                i   =  i+1
                if(A(i) .ge. pvt)   exit
            end do
            do
                j   =  j-1
                if(A(j) .le. pvt)   exit
            end do
            if(i .lt. j) then
!               swap A(i) and A(j).
                itmp=  A(i)
                A(i)=  A(j)
                A(j)=  itmp
            else
                exit
            end if
        end do

!       now i>=j and A(L)<=A(i)=pvt<=A(R). For the two sequences A(L:i-1) and
!       A(j+1:R), if both are more than mlen long, push the longer one to the
!       stack and go back to sort the other one; if only one is more than mlen
!       long, go back to sort this one; or pop a sequence from the stack and 
!       sort it.

        if(R-j .gt. i-L) then
!           store A(j+1:R) and sort A(L:i-1).
            if(R .gt. j+1) then
                nstk=  nstk+1
                if(nstk .gt. mstk)  stop 'Error: stack size too small.'
                Lstk(nstk)  =  j+1
                Rstk(nstk)  =  R
            end if
            R   =  i-1
        else
!           store A(L:i-1) and sort A(j+1:R).
            if(i .gt. L+1) then
                nstk=  nstk+1
                if(nstk .gt. mstk)  stop 'Error: stack size too small.'
                Lstk(nstk)  =  L
                Rstk(nstk)  =  i-1
            end if
            L   =  j+1
        end if
    end do

    if(.not. asc) then
        M   =  lef+rig
        do i=lef,lef+(rig-lef-1)/2
            j   =  M-i
            itmp=  A(i)
            A(i)=  A(j)
            A(j)=  itmp
        end do
    end if

    return
    end subroutine iqsort
!-------------------------------------------------------------------------------
!   quicksort for double precision real.
!-------------------------------------------------------------------------------
    subroutine dqsort(asc,lef,rig,A)
    use var_kind_def
    implicit none
    integer(dpI),parameter:: mlen=9,mstk=30
    logical(dpL),intent(in):: asc
    integer(dpI),intent(in):: lef,rig
    integer(dpI):: L,R,M,i,j,nstk,Lstk(mstk),Rstk(mstk)
    real   (dpR):: A(*),pvt,rtmp

    nstk=  0
    L   =  lef
    R   =  rig

    do
!       use insertion sort to sort A(L:R) and pop new sequence from stack.
        if(R-L+1 .le. mlen) then
            do j=L+1,R
                rtmp=  A(j)
                do i=j-1,L,-1
                    if(A(i) .le. rtmp)  exit
                    A(i+1)  =  A(i)
                end do
                A(i+1)  =  rtmp
            end do

            if(nstk .le. 0) exit
            L   =  Lstk(nstk)
            R   =  Rstk(nstk)
            nstk=  nstk-1
        end if
        if(R .le. L) then
            if(nstk .le. 0) exit
            L   =  Lstk(nstk)
            R   =  Rstk(nstk)
            nstk=  nstk-1
        end if

        M   =  L+(R-L)/2
        if(A(L) .gt. A(M)) then
            rtmp=  A(L)
            A(L)=  A(M)
            A(M)=  rtmp
        end if
        if(A(M) .gt. A(R)) then
            rtmp=  A(M)
            A(M)=  A(R)
            A(R)=  rtmp
        end if
        if(A(L) .gt. A(M)) then
            rtmp=  A(L)
            A(L)=  A(M)
            A(M)=  rtmp
        end if
        pvt =  A(M)

!       now A(L)<=pvt<=A(R), we have to find the final location of pvt.
        i   =  L
        j   =  R
        do
            do
                i   =  i+1
                if(A(i) .ge. pvt)   exit
            end do
            do
                j   =  j-1
                if(A(j) .le. pvt)   exit
            end do
            if(i .lt. j) then
!               swap A(i) and A(j).
                rtmp=  A(i)
                A(i)=  A(j)
                A(j)=  rtmp
            else
                exit
            end if
        end do

!       now i>=j and A(L)<=A(i)=pvt<=A(R). For the two sequences A(L:i-1) and
!       A(j+1:R), if both are more than mlen long, push the longer one to the
!       stack and go back to sort the other one; if only one is more than mlen
!       long, go back to sort this one; or pop a sequence from the stack and 
!       sort it.

        if(R-j .gt. i-L) then
!           store A(j+1:R) and sort A(L:i-1).
            if(R .gt. j+1) then
                nstk=  nstk+1
                if(nstk .gt. mstk)  stop 'Error: stack size too small.'
                Lstk(nstk)  =  j+1
                Rstk(nstk)  =  R
            end if
            R   =  i-1
        else
!           store A(L:i-1) and sort A(j+1:R).
            if(i .gt. L+1) then
                nstk=  nstk+1
                if(nstk .gt. mstk)  stop 'Error: stack size too small.'
                Lstk(nstk)  =  L
                Rstk(nstk)  =  i-1
            end if
            L   =  j+1
        end if
    end do

    if(.not. asc) then
        M   =  lef+rig
        do i=lef,lef+(rig-lef-1)/2
            j   =  M-i
            rtmp=  A(i)
            A(i)=  A(j)
            A(j)=  rtmp
        end do
    end if

    return
    end subroutine dqsort
!-------------------------------------------------------------------------------
!   quicksort, integer.
!-------------------------------------------------------------------------------
    subroutine iqsortcols(asc,lef,rig,row,bsze,A)
    use var_kind_def
    implicit none
    integer(dpI),parameter:: mlen=9,mstk=30
    logical(dpL),intent(in):: asc
    integer(dpI),intent(in):: lef,rig,row,bsze
    integer(dpI):: A(*),L,R,M,i,j,nstk,pvt,Lstk(mstk),Rstk(mstk),v(bsze),itmp

    nstk=  0
    L   =  lef
    R   =  rig

    do
!       use insertion sort to sort A(L:R) and pop new sequence from stack.
        if(R-L+1 .le. mlen) then
            do j=L+1,R
                itmp=  A(row+bsze*(j-1))
                v(1:bsze)   =  A(1+bsze*(j-1):bsze*j)
                do i=j-1,L,-1
                    if(A(row+bsze*(i-1)) .le. itmp) exit
                    A(1+bsze*i:bsze*(i+1))  =  A(1+bsze*(i-1):bsze*i)
                end do
                A(1+bsze*i:bsze*(i+1))  =  v(1:bsze)
            end do

            if(nstk .le. 0) exit
            L   =  Lstk(nstk)
            R   =  Rstk(nstk)
            nstk=  nstk-1
        end if
        if(R .le. L) then
            if(nstk .le. 0) exit
            L   =  Lstk(nstk)
            R   =  Rstk(nstk)
            nstk=  nstk-1
        end if

        M   =  L+(R-L)/2
        if(A(row+bsze*(L-1)) .gt. A(row+bsze*(M-1))) then
            V(1:bsze)   =  A(1+bsze*(L-1):bsze*L)
            A(1+bsze*(L-1):bsze*L)  =  A(1+bsze*(M-1):bsze*M)
            A(1+bsze*(M-1):bsze*M)  =  V(1:bsze)
        end if
        if(A(row+bsze*(M-1)) .gt. A(row+bsze*(R-1))) then
            V(1:bsze)   =  A(1+bsze*(M-1):bsze*M)
            A(1+bsze*(M-1):bsze*M)  =  A(1+bsze*(R-1):bsze*R)
            A(1+bsze*(R-1):bsze*R)  =  V(1:bsze)
        end if
        if(A(row+bsze*(L-1)) .gt. A(row+bsze*(M-1))) then
            V(1:bsze)   =  A(1+bsze*(L-1):bsze*L)
            A(1+bsze*(L-1):bsze*L)  =  A(1+bsze*(M-1):bsze*M)
            A(1+bsze*(M-1):bsze*M)  =  V(1:bsze)
        end if
        pvt =  A(row+bsze*(M-1))

!       now A(L)<=pvt<=A(R), we have to find the final location of pvt.
        i   =  L
        j   =  R
        do
            do
                i   =  i+1
                if(A(row+bsze*(i-1)) .ge. pvt)  exit
            end do
            do
                j   =  j-1
                if(A(row+bsze*(j-1)) .le. pvt)  exit
            end do
            if(i .lt. j) then
!               swap A(i) and A(j).
                V(1:bsze)   =  A(1+bsze*(i-1):bsze*i)
                A(1+bsze*(i-1):bsze*i)  =  A(1+bsze*(j-1):bsze*j)
                A(1+bsze*(j-1):bsze*j)  =  V(1:bsze)
            else
                exit
            end if
        end do

!       now i>=j and A(L)<=A(i)=pvt<=A(R). For the two sequences A(L:i-1) and
!       A(j+1:R), if both are more than mlen long, push the longer one to the
!       stack and go back to sort the other one; if only one is more than mlen
!       long, go back to sort this one; or pop a sequence from the stack and 
!       sort it.

        if(R-j .gt. i-L) then
!           store A(j+1:R) and sort A(L:i-1).
            if(R .gt. j+1) then
                nstk=  nstk+1
                if(nstk .gt. mstk)  stop 'Error: stack size too small.'
                Lstk(nstk)  =  j+1
                Rstk(nstk)  =  R
            end if
            R   =  i-1
        else
!           store A(L:i-1) and sort A(j+1:R).
            if(i .gt. L+1) then
                nstk=  nstk+1
                if(nstk .gt. mstk)  stop 'Error: stack size too small.'
                Lstk(nstk)  =  L
                Rstk(nstk)  =  i-1
            end if
            L   =  j+1
        end if
    end do

    if(.not. asc) then
        M   =  lef+rig
        do i=lef,lef+(rig-lef-1)/2
            j   =  M-i
            V(1:bsze)   =  A(1+bsze*(i-1):bsze*i)
            A(1+bsze*(i-1):bsze*i)  =  A(1+bsze*(j-1):bsze*j)
            A(1+bsze*(j-1):bsze*j)  =  V(1:bsze)
        end do
    end if

    return
    end subroutine iqsortcols
!-------------------------------------------------------------------------------
!   quicksort of integer array, plus a real array.
!-------------------------------------------------------------------------------
    subroutine iqsortcols_mat(asc,lef,rig,row,bsze,A,bsze_d,D)
    use var_kind_def
    implicit none
    integer(dpI),parameter:: mlen=9,mstk=30
    logical(dpL),intent(in):: asc
    integer(dpI),intent(in):: lef,rig,row,bsze,bsze_d
    integer(dpI):: A(*),L,R,M,i,j,nstk,pvt,Lstk(mstk),Rstk(mstk),v(bsze),itmp
    real   (dpR):: D(*),C(bsze_d)

    nstk=  0
    L   =  lef
    R   =  rig

    do
!       use insertion sort to sort A(L:R) and pop new sequence from stack.
        if(R-L+1 .le. mlen) then
            do j=L+1,R
                itmp=  A(row+bsze*(j-1))
                v(1:bsze  ) =  A(1+bsze  *(j-1):bsze  *j)
                C(1:bsze_d) =  D(1+bsze_d*(j-1):bsze_d*j)
                do i=j-1,L,-1
                    if(A(row+bsze*(i-1)) .le. itmp) exit
                    A(1+bsze  *i:bsze  *(i+1))  =  A(1+bsze  *(i-1):bsze  *i)
                    D(1+bsze_d*i:bsze_d*(i+1))  =  D(1+bsze_d*(i-1):bsze_d*i)
                end do
                A(1+bsze  *i:bsze  *(i+1))  =  v(1:bsze  )
                D(1+bsze_d*i:bsze_d*(i+1))  =  C(1:bsze_d)
            end do

            if(nstk .le. 0) exit
            L   =  Lstk(nstk)
            R   =  Rstk(nstk)
            nstk=  nstk-1
        end if
        if(R .le. L) then
            if(nstk .le. 0) exit
            L   =  Lstk(nstk)
            R   =  Rstk(nstk)
            nstk=  nstk-1
        end if

        M   =  L+(R-L)/2
        if(A(row+bsze*(L-1)) .gt. A(row+bsze*(M-1))) then
            V(1:bsze)   =  A(1+bsze*(L-1):bsze*L)
            A(1+bsze*(L-1):bsze*L)  =  A(1+bsze*(M-1):bsze*M)
            A(1+bsze*(M-1):bsze*M)  =  V(1:bsze)

            C(1:bsze_d) =  D(1+bsze_d*(L-1):bsze_d*L)
            D(1+bsze_d*(L-1):bsze_d*L)  =  D(1+bsze_d*(M-1):bsze_d*M)
            D(1+bsze_d*(M-1):bsze_d*M)  =  C(1:bsze_d)
        end if
        if(A(row+bsze*(M-1)) .gt. A(row+bsze*(R-1))) then
            V(1:bsze)   =  A(1+bsze*(M-1):bsze*M)
            A(1+bsze*(M-1):bsze*M)  =  A(1+bsze*(R-1):bsze*R)
            A(1+bsze*(R-1):bsze*R)  =  V(1:bsze)

            C(1:bsze_d) =  D(1+bsze_d*(M-1):bsze_d*M)
            D(1+bsze_d*(M-1):bsze_d*M)  =  D(1+bsze_d*(R-1):bsze_d*R)
            D(1+bsze_d*(R-1):bsze_d*R)  =  C(1:bsze_d)
        end if
        if(A(row+bsze*(L-1)) .gt. A(row+bsze*(M-1))) then
            V(1:bsze)   =  A(1+bsze*(L-1):bsze*L)
            A(1+bsze*(L-1):bsze*L)  =  A(1+bsze*(M-1):bsze*M)
            A(1+bsze*(M-1):bsze*M)  =  V(1:bsze)

            C(1:bsze_d) =  D(1+bsze_d*(L-1):bsze_d*L)
            D(1+bsze_d*(L-1):bsze_d*L)  =  D(1+bsze_d*(M-1):bsze_d*M)
            D(1+bsze_d*(M-1):bsze_d*M)  =  C(1:bsze_d)
        end if
        pvt =  A(row+bsze*(M-1))

!       now A(L)<=pvt<=A(R), we have to find the final location of pvt.
        i   =  L
        j   =  R
        do
            do
                i   =  i+1
                if(A(row+bsze*(i-1)) .ge. pvt)  exit
            end do
            do
                j   =  j-1
                if(A(row+bsze*(j-1)) .le. pvt)  exit
            end do
            if(i .lt. j) then
!               swap A(i) and A(j).
                V(1:bsze)   =  A(1+bsze*(i-1):bsze*i)
                A(1+bsze*(i-1):bsze*i)  =  A(1+bsze*(j-1):bsze*j)
                A(1+bsze*(j-1):bsze*j)  =  V(1:bsze)

                C(1:bsze_d) =  D(1+bsze_d*(i-1):bsze_d*i)
                D(1+bsze_d*(i-1):bsze_d*i)  =  D(1+bsze_d*(j-1):bsze_d*j)
                D(1+bsze_d*(j-1):bsze_d*j)  =  C(1:bsze_d)
            else
                exit
            end if
        end do

!       now i>=j and A(L)<=A(i)=pvt<=A(R). For the two sequences A(L:i-1) and
!       A(j+1:R), if both are more than mlen long, push the longer one to the
!       stack and go back to sort the other one; if only one is more than mlen
!       long, go back to sort this one; or pop a sequence from the stack and 
!       sort it.

        if(R-j .gt. i-L) then
!           store A(j+1:R) and sort A(L:i-1).
            if(R .gt. j+1) then
                nstk=  nstk+1
                if(nstk .gt. mstk)  stop 'Error: stack size too small.'
                Lstk(nstk)  =  j+1
                Rstk(nstk)  =  R
            end if
            R   =  i-1
        else
!           store A(L:i-1) and sort A(j+1:R).
            if(i .gt. L+1) then
                nstk=  nstk+1
                if(nstk .gt. mstk)  stop 'Error: stack size too small.'
                Lstk(nstk)  =  L
                Rstk(nstk)  =  i-1
            end if
            L   =  j+1
        end if
    end do

    if(.not. asc) then
        M   =  lef+rig
        do i=lef,lef+(rig-lef-1)/2
            j   =  M-i
            V(1:bsze)   =  A(1+bsze*(i-1):bsze*i)
            A(1+bsze*(i-1):bsze*i)  =  A(1+bsze*(j-1):bsze*j)
            A(1+bsze*(j-1):bsze*j)  =  V(1:bsze)

            C(1:bsze_d) =  D(1+bsze_d*(i-1):bsze_d*i)
            D(1+bsze_d*(i-1):bsze_d*i)  =  D(1+bsze_d*(j-1):bsze_d*j)
            D(1+bsze_d*(j-1):bsze_d*j)  =  C(1:bsze_d)
        end do
    end if

    return
    end subroutine iqsortcols_mat
!-------------------------------------------------------------------------------
!   quicksort for double precision real, column.
!-------------------------------------------------------------------------------
    subroutine dqsortcols(asc,lef,rig,row,bsze,A)
    use var_kind_def
    implicit none
    integer(dpI),parameter:: mlen=9,mstk=30
    logical(dpL),intent(in):: asc
    integer(dpI),intent(in):: lef,rig,row,bsze
    integer(dpI):: L,R,M,i,j,nstk,Lstk(mstk),Rstk(mstk)
    real   (dpR):: A(*),pvt,v(bsze),rtmp

    nstk=  0
    L   =  lef
    R   =  rig

    do
!       use insertion sort to sort A(L:R) and pop new sequence from stack.
        if(R-L+1 .le. mlen) then
            do j=L+1,R
                rtmp=  A(row+bsze*(j-1))
                v(1:bsze)   =  A(1+bsze*(j-1):bsze*j)
                do i=j-1,L,-1
                    if(A(row+bsze*(i-1)) .le. rtmp) exit
                    A(1+bsze*i:bsze*(i+1))  =  A(1+bsze*(i-1):bsze*i)
                end do
                A(1+bsze*i:bsze*(i+1))  =  v(1:bsze)
            end do

            if(nstk .le. 0) exit
            L   =  Lstk(nstk)
            R   =  Rstk(nstk)
            nstk=  nstk-1
        end if
        if(R .le. L) then
            if(nstk .le. 0) exit
            L   =  Lstk(nstk)
            R   =  Rstk(nstk)
            nstk=  nstk-1
        end if

        M   =  L+(R-L)/2
        if(A(row+bsze*(L-1)) .gt. A(row+bsze*(M-1))) then
            V(1:bsze)   =  A(1+bsze*(L-1):bsze*L)
            A(1+bsze*(L-1):bsze*L)  =  A(1+bsze*(M-1):bsze*M)
            A(1+bsze*(M-1):bsze*M)  =  V(1:bsze)
        end if
        if(A(row+bsze*(M-1)) .gt. A(row+bsze*(R-1))) then
            V(1:bsze)   =  A(1+bsze*(M-1):bsze*M)
            A(1+bsze*(M-1):bsze*M)  =  A(1+bsze*(R-1):bsze*R)
            A(1+bsze*(R-1):bsze*R)  =  V(1:bsze)
        end if
        if(A(row+bsze*(L-1)) .gt. A(row+bsze*(M-1))) then
            V(1:bsze)   =  A(1+bsze*(L-1):bsze*L)
            A(1+bsze*(L-1):bsze*L)  =  A(1+bsze*(M-1):bsze*M)
            A(1+bsze*(M-1):bsze*M)  =  V(1:bsze)
        end if
        pvt =  A(row+bsze*(M-1))

!       now A(L)<=pvt<=A(R), we have to find the final location of pvt.
        i   =  L
        j   =  R
        do
            do
                i   =  i+1
                if(A(row+bsze*(i-1)) .ge. pvt)  exit
            end do
            do
                j   =  j-1
                if(A(row+bsze*(j-1)) .le. pvt)  exit
            end do
            if(i .lt. j) then
!               swap A(i) and A(j).
                V(1:bsze)   =  A(1+bsze*(i-1):bsze*i)
                A(1+bsze*(i-1):bsze*i)  =  A(1+bsze*(j-1):bsze*j)
                A(1+bsze*(j-1):bsze*j)  =  V(1:bsze)
            else
                exit
            end if
        end do

!       now i>=j and A(L)<=A(i)=pvt<=A(R). For the two sequences A(L:i-1) and
!       A(j+1:R), if both are more than mlen long, push the longer one to the
!       stack and go back to sort the other one; if only one is more than mlen
!       long, go back to sort this one; or pop a sequence from the stack and 
!       sort it.

        if(R-j .gt. i-L) then
!           store A(j+1:R) and sort A(L:i-1).
            if(R .gt. j+1) then
                nstk=  nstk+1
                if(nstk .gt. mstk)  stop 'Error: stack size too small.'
                Lstk(nstk)  =  j+1
                Rstk(nstk)  =  R
            end if
            R   =  i-1
        else
!           store A(L:i-1) and sort A(j+1:R).
            if(i .gt. L+1) then
                nstk=  nstk+1
                if(nstk .gt. mstk)  stop 'Error: stack size too small.'
                Lstk(nstk)  =  L
                Rstk(nstk)  =  i-1
            end if
            L   =  j+1
        end if
    end do

    if(.not. asc) then
        M   =  lef+rig
        do i=lef,lef+(rig-lef-1)/2
            j   =  M-i
            V(1:bsze)   =  A(1+bsze*(i-1):bsze*i)
            A(1+bsze*(i-1):bsze*i)  =  A(1+bsze*(j-1):bsze*j)
            A(1+bsze*(j-1):bsze*j)  =  V(1:bsze)
        end do
    end if

    return
    end subroutine dqsortcols
!-------------------------------------------------------------------------------
!   sort matrix, column by column, two tags.
!-------------------------------------------------------------------------------
    subroutine sort_matrix_col(LDA,L,R,tag1,tag2,D)
    use var_kind_def
    implicit none
    integer(dpI):: LDA,L,R,tag1,tag2,D(LDA,*),i,j,k

    if(L .ge. R)    return
    call iqsortcols(.true., L, R, tag1, LDA, D)
    i   =  L
    do while(i .le. R)
        k   =  i
        do j=i+1,R
            if(D(tag1,j) .ne. D(tag1,i)) then
                exit
            else
                k   =  j
            end if
        end do
        if(k .gt. i)    call iqsortcols(.true., i, k, tag2, LDA, D)

        i   =  k+1
    end do

    return
    end subroutine sort_matrix_col
!-------------------------------------------------------------------------------
!   sort matrix, first row-by-row(in column), then column-by-column.
!-------------------------------------------------------------------------------
    subroutine sort_matrix_row_col(npe,L,R,tag,LDA,mat)
    use var_kind_def
    implicit none
    integer(dpI):: npe,L,R,tag,LDA,mat(LDA,*),ifac,i,j,itmp

    do ifac=L,R
        do i=npe-1,1,-1
        do j=1,i
            if(mat(j,ifac) .gt. mat(j+1,ifac)) then
                itmp=  mat(j,ifac)
                mat(j  ,ifac)   =  mat(j+1,ifac)
                mat(j+1,ifac)   =  itmp
            end if
        end do
        end do
    end do
    call iqsortcols(.true., L, R, tag, LDA, mat)

    return
    end subroutine sort_matrix_row_col
