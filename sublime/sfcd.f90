!-------------------------------------------------------------------------------
!   cal the Lagrange interpolation coefficients.
!-------------------------------------------------------------------------------
    subroutine cal_Lagrange_coe(M,x0,x,coe)
    use var_kind_def
    implicit none
    integer(dpI),intent(in):: M
    real   (dpR),intent(in):: x0(*),x
    integer(dpI):: j,k
    real   (dpR):: coe(*)

    do j=1,M
        coe(j)  =  1.0d0
        do k=1,M
            if(k .ne. j)    coe(j)  =  coe(j)*(x-x0(k))/(x0(j)-x0(k))
        end do
    end do

    return
    end subroutine cal_Lagrange_coe
!-------------------------------------------------------------------------------
!   cal the differential coefficients with the Lagrange interpolation.
!-------------------------------------------------------------------------------
    subroutine cal_Lagrange_dcoe(M,x0,x,coe)
    use var_kind_def
    implicit none
    integer(dpI),intent(in):: M
    real   (dpR),intent(in):: x0(*),x
    integer(dpI):: j,k,r
    real   (dpR):: coe(*),c

    do j=1,M
        coe(j)  =  0.0d0
        do r=1,M
            if(r .eq. j)    cycle

            c   =  1.0d0
            do k=1,M
                if((k .eq. j) .or. (k .eq. r))  cycle
                c   =  c*(x-x0(k))/(x0(j)-x0(k))
            end do
            coe(j)  =  coe(j)+c/(x0(j)-x0(r))
        end do
    end do

    return
    end subroutine cal_Lagrange_dcoe
!-------------------------------------------------------------------------------
!   compute the oriented boundary box for 2D/3D points set.
!-------------------------------------------------------------------------------
    subroutine cal_obb(LDA,npnt,dim,xyz,c,xyz2obb,L,R)
    use var_kind_def
    implicit none
    integer(dpI),intent(in):: LDA,npnt,dim
    real   (dpR),intent(in):: xyz(LDA,*)
    integer(dpI):: i,j
    real   (dpR):: c(*),xyz2obb(*),L(*),R(*),LHS(dim,dim),d(dim)

    if((dim .lt. 2) .or. (dim .gt. 3))  stop 'Error: cal_obb works for 2D/3D only.'
    c(1:dim)=  0.0d0
    do i=1,npnt
        c(1:dim)=  c(1:dim)+xyz(1:dim,i)
    end do
    c(1:dim)=  c(1:dim)/real(npnt, dpR)
    LHS =  0.0d0
    do i=1,npnt
        d(1:dim)=  xyz(1:dim,i)-c(1:dim)
        do j=1,dim
            LHS(1:dim,j)=  LHS(1:dim,j)+d(1:dim)*d(j)
        end do
    end do
    call DSYEV(dim, LHS, L, xyz2obb)

    if(abs(L(dim)) .gt. abs(L(1))) then
        do i=1,dim/2
            j   =  dim+1-i
            L(1:dim)    =  xyz2obb(1+(i-1)*dim:i*dim)
            xyz2obb(1+(i-1)*dim:i*dim)  =  xyz2obb(1+(j-1)*dim:j*dim)
            xyz2obb(1+(j-1)*dim:j*dim)  =  L(1:dim)
        end do
    end if
    do i=1,dim
        call norm_vec(dim, xyz2obb(1+(i-1)*dim), L(1))
    end do
    call mat_inv(dim, xyz2obb)

    L(1:dim)=  huge(1.0d0)
    R(1:dim)= -L(1:dim)
    do i=1,npnt
        d(1:dim)=  xyz(1:dim,i)-c(1:dim)
        LHS(1:dim,1)=  xyz2obb(1:dim)*d(1)
        do j=2,dim
            LHS(1:dim,1)=  LHS(1:dim,1)+xyz2obb(1+(j-1)*dim:j*dim)*d(j)
        end do
        L(1:dim)=  min(L(1:dim), LHS(1:dim,1))
        R(1:dim)=  max(R(1:dim), LHS(1:dim,1))
    end do

    return
    end subroutine cal_obb
!-------------------------------------------------------------------------------
!   cal the rotation matrix.
!-------------------------------------------------------------------------------
    subroutine cal_rot_mat(a,M)
    use var_kind_def
    implicit none
    real(dpR):: a(3),M(3,*),Mx(3,3),My(3,3),Mz(3,3)

    Mx(1:3,1)   = (/1.0d0, 0.0d0    , 0.0d0    /)
    Mx(1:3,2)   = (/0.0d0, cos(a(1)), sin(a(1))/)
    Mx(1:3,3)   = (/0.0d0,-sin(a(1)), cos(a(1))/)

    My(1:3,1)   = (/cos(a(2)), 0.0d0,-sin(a(2))/)
    My(1:3,2)   = (/0.0d0    , 1.0d0, 0.0d0    /)
    My(1:3,3)   = (/sin(a(2)), 0.0d0, cos(a(2))/)

    Mz(1:3,1)   = (/cos(a(3)), sin(a(3)), 0.0d0/)
    Mz(1:3,2)   = (/-sin(a(3)),cos(a(3)), 0.0d0/)
    Mz(1:3,3)   = (/0.0d0     , 0.0d0   , 1.0d0/)
    M(1:3,1:3)  =  matmul(Mz, matmul(My, Mx))

    return
    end subroutine cal_rot_mat
!-------------------------------------------------------------------------------
!   Calculate the cross product of two vectors.
!-------------------------------------------------------------------------------
    subroutine crs_prd(x,y,z)
    use var_kind_def
    implicit none
    real(dpR):: x(*),y(*),z(*)

    z(1)=  x(2)*y(3)-x(3)*y(2)
    z(2)= -x(1)*y(3)+x(3)*y(1)
    z(3)=  x(1)*y(2)-x(2)*y(1)

    return
    end subroutine crs_prd
!-------------------------------------------------------------------------------
!   y=a*x, double precision.
!-------------------------------------------------------------------------------
    subroutine DAXY(N,a,x,y)
    use var_kind_def
    implicit none
    integer(dpI),intent(in):: N
    real   (dpR),intent(in):: a,x(*)
    integer(dpI):: i,M
    real   (dpR):: y(*)

    M   =  mod(N, 7)
    forall(i=1:M)   y(i)=  a*x(i)

    do i=M+1,N,7
        y(i  )  =  a*x(i  )
        y(i+1)  =  a*x(i+1)
        y(i+2)  =  a*x(i+2)
        y(i+3)  =  a*x(i+3)
        y(i+4)  =  a*x(i+4)
        y(i+5)  =  a*x(i+5)
        y(i+6)  =  a*x(i+6)
    end do

    return
    end subroutine DAXY
!-------------------------------------------------------------------------------
!   z=b*y+a*x, double precision.
!-------------------------------------------------------------------------------
    subroutine DAXPBYZ(N,a,x,b,y,z)
    implicit none
    integer(kind=4),intent(in):: N
    real   (kind=8),intent(in):: a,b,x(*),y(*)
    integer(kind=4):: i,M
    real   (kind=8):: z(*)

    M   =  mod(N, 7)
    forall(i=1:M)   z(i)=  b*y(i)+a*x(i)

    do i=M+1,N,7
        z(i  )  =  b*y(i  )+a*x(i  )
        z(i+1)  =  b*y(i+1)+a*x(i+1)
        z(i+2)  =  b*y(i+2)+a*x(i+2)
        z(i+3)  =  b*y(i+3)+a*x(i+3)
        z(i+4)  =  b*y(i+4)+a*x(i+4)
        z(i+5)  =  b*y(i+5)+a*x(i+5)
        z(i+6)  =  b*y(i+6)+a*x(i+6)
    end do

    return
    end subroutine DAXPBYZ
!-------------------------------------------------------------------------------
!   generating orthogonal coordinate system.
!-------------------------------------------------------------------------------
    subroutine gen_coord_n3(n1,n2,n3)
    use var_kind_def
    implicit none
    integer(dpI):: K
    real   (dpR):: n1(*),n2(*),n3(*),rtmp

    K       =  minloc(abs(n1(1:3)), dim=1)
    rtmp    =  n1(K)
    n2(1:3) =  rtmp*n1(1:3)
    n2(K)   =  rtmp*rtmp-1.0d0
    call norm_vec(3, n2, rtmp)
    n3(1)   =  n1(2)*n2(3)-n1(3)*n2(2)
    n3(2)   = -n1(1)*n2(3)+n1(3)*n2(1)
    n3(3)   =  n1(1)*n2(2)-n1(2)*n2(1)

    return
    end subroutine gen_coord_n3
!-------------------------------------------------------------------------------
!   Def the rotation matrix for periodic bnds.
!-------------------------------------------------------------------------------
    subroutine per_mat_rot(unt,ang,rmat)
    use var_kind_def
    implicit none
    real(dpR):: unt(3),ang,rmat(3,3),c,t,s

!   p_new=  per_mat*(p_old-axis_start_point)+axis_start_point

    c   =  cos(ang)
    t   =  1.0d0-c
    s   =  sin(ang)

    rmat(1,1)   =  t*unt(1)*unt(1)+c
    rmat(2,1)   =  t*unt(2)*unt(1)+s*unt(3)
    rmat(3,1)   =  t*unt(3)*unt(1)-s*unt(2)
    rmat(1,2)   =  t*unt(1)*unt(2)-s*unt(3)
    rmat(2,2)   =  t*unt(2)*unt(2)+c
    rmat(3,2)   =  t*unt(3)*unt(2)+s*unt(1)
    rmat(1,3)   =  t*unt(1)*unt(3)+s*unt(2)
    rmat(2,3)   =  t*unt(2)*unt(3)-s*unt(1)
    rmat(3,3)   =  t*unt(3)*unt(3)+c

    return
    end subroutine per_mat_rot
!-------------------------------------------------------------------------------
!   evaluate polynomial values.
!-------------------------------------------------------------------------------
    subroutine polyval(order,LDA,N,p,x,y)
    use var_kind_def
    implicit none
    integer(dpI),intent(in):: order,LDA,N
    real   (dpR),intent(in):: x(*)
    integer(dpI):: i,j,j1,j0
    real   (dpR):: p(*),y(*)

    do j=1,N
        j1  =  1+LDA*(j-1)
        j0  =    LDA* j
        y(j1:j0)=  p(1:LDA)
        do i=1,order-1
            y(j1:j0)=  y(j1:j0)*x(j)+p(1+LDA*i:LDA*(i+1))
        end do
    end do

    return
    end subroutine polyval
!-------------------------------------------------------------------------------
!   cal the differential of polynominal.
!-------------------------------------------------------------------------------
    subroutine polyder(N,p,dp)
    use var_kind_def
    implicit none
    integer(dpI),intent(in):: N
    real   (dpR),intent(in):: p(*)
    integer(dpI):: i
    real   (dpR):: dp(*)

    do i=1,N-1
        dp(i)   =  p(i)*real(N-i, dpR)
    end do

    return
    end subroutine polyder
!-------------------------------------------------------------------------------
!   transpose matrix.
!-------------------------------------------------------------------------------
    subroutine matrix_transpose(LDA,A)
    use var_kind_def
    implicit none
    integer(dpI):: LDA,i,j
    real   (dpR):: A(*),rtmp

    do i=1,LDA
    do j=1,i-1
        rtmp=  A(i+LDA*(j-1))
        A(i+LDA*(j-1))  =  A(j+LDA*(i-1))
        A(j+LDA*(i-1))  =  rtmp
    end do
    end do

    return
    end subroutine matrix_transpose
!-------------------------------------------------------------------------------
!   Normalize given vec with error control.
!-------------------------------------------------------------------------------
    subroutine norm_vec(N,v,nrm)
    use var_kind_def
    implicit none
    integer(dpI):: N,i
    real   (dpR):: v(*),nrm

    nrm =  v(1)*v(1)
    do i=2,N
        nrm =  nrm+v(i)*v(i)
    end do
    if(nrm .le. 0.0d0)  return
    nrm =  sqrt(nrm)
    v(1:N)  =  v(1:N)/nrm

    return
    end subroutine norm_vec
!-------------------------------------------------------------------------------
!   convert string to uppercase.
!-------------------------------------------------------------------------------
    subroutine upper_string(str)
    implicit none
    character(len=*):: str
    integer(kind=4):: i,j

    do i=1,len_trim(str)
        j   =  ichar(str(i:i))
        if(j .ge. 97 .and. j .le. 122)  str(i:i)=  char(j-32)
    end do
    str =  trim(str)

    return
    end subroutine upper_string
!----------------------------------------------------------------------------------
!   Zonename.
!----------------------------------------------------------------------------------
    subroutine num_to_str(mblk,str,zonename)
    implicit none
    integer  (kind=4):: blk,blkn,mblk,i,j,k,m
    character(len=* ):: str
    character(len=80):: zonename(mblk)

    if(mblk .le. 9) then
        do blk=1,mblk
            zonename(blk)   =  str(1:len_trim(str))//'_'//char(blk+48)
        end do
    elseif(mblk .le. 99) then
        do blk=1,mblk
            i   =  blk/10
            j   =  blk-10*i
            zonename(blk)   =  str(1:len_trim(str))//'_'//char(i+48)//char(j+48)
        end do
    elseif(mblk .le. 999) then
        do blk=1,mblk
            blkn=  blk
            i   =  blkn/100
            blkn=  blkn-i*100
            j   =  blkn/10
            k   =  blkn-j*10
            zonename(blk)   =  str(1:len_trim(str))//'_'//char(i+48)//char(j+48)//char(k+48)
        end do
    else
        do blk=1,mblk
            blkn=  mod(blk,10000)
            i   =  blkn/1000
            blkn=  blkn-i*1000
            j   =  blkn/100
            blkn=  blkn-j*100
            k   =  blkn/10
            m   =  blkn-10*k
            zonename(blk)   =  str(1:len_trim(str))//'_'//char(i+48)//char(j+48)//char(k+48)//char(m+48)
        end do
    end if

    return
    end subroutine num_to_str
!-------------------------------------------------------------------------------
!   Rotate and offset vec.
!-------------------------------------------------------------------------------
    subroutine rotate_vector(N,rotm,vec)
    use var_kind_def
    implicit none
    integer(dpI),intent(in):: N
    real   (dpR),intent(in):: rotm(*)
    integer(dpI):: i
    real   (dpR):: vec(*),rtm(3)

    do i=1,N
        rtm(1:3)=  rotm(1:3)*vec(3*i-2)+rotm(4:6)*vec(3*i-1)+rotm(7:9)*vec(3*i)
        vec(3*i-2:3*i)  =  rtm(1:3)
    end do

    return
    end subroutine rotate_vector
!-------------------------------------------------------------------------------
!   sort and simplify a series of numbers, real.
!-------------------------------------------------------------------------------
    subroutine dsimplify_series(N,col,LDA,v)
    use var_kind_def
    implicit none
    integer(dpI),intent(in):: col,LDA
    integer(dpI):: N,N0,i,j,k
    real   (dpR):: v(*)

    if(N .le. 0)    return
    call dqsortcols(.true., 1, N, col, LDA, v)
    i   =  1
    N0  =  N
    N   =  0
    do while(i .le. N0)
        k   =  i
        do j=i+1,N0
            if(nint(v(col+LDA*(i-1)), dpI) .ne. nint(v(col+LDA*(j-1)), dpI)) then
                exit
            else
                k   =  j
            end if
        end do
        N   =  N+1
        v(1+LDA*(N-1):LDA*N)=  v(1+LDA*(i-1):LDA*i)

        i   =  k+1
    end do

    return
    end subroutine dsimplify_series
!-------------------------------------------------------------------------------
!   xyz to RTZ.
!-------------------------------------------------------------------------------
    subroutine xyz_to_rtz(is_rot,n,xyz,rtz)
    use var_kind_def
    use var_global, only: pi,rot_axs=>rotation_axis,trs_axs=>translation_axis
    implicit none
    logical(dpI),intent(in):: is_rot
    integer(dpI),intent(in):: n
    real   (dpR),intent(in):: xyz(3,*)
    integer(dpI):: i
    real   (dpR):: rtz(*),d(3),rr,tt,a,b,rad_axs(3)

    if(is_rot) then
        call crs_prd(trs_axs, rot_axs, rad_axs)
        do i=1,N
            rtz(3*i)=  xyz(1,i)*rot_axs(1)+xyz(2,i)*rot_axs(2)+xyz(3,i)*rot_axs(3)
            d(1:3)  =  xyz(1:3,i)-rtz(3*i)*rot_axs(1:3)
            rr      =  d(1)*rad_axs(1)+d(2)*rad_axs(2)+d(3)*rad_axs(3)
            tt      =  d(1)*trs_axs(1)+d(2)*trs_axs(2)+d(3)*trs_axs(3)
            b           =  atan2(tt, rr)
            rtz(3*i-2)  =  sqrt(rr*rr+tt*tt)
            if(i .gt. 1) then
                if(b .ge. a+0.5d0*pi) then
                    b   =  b-2.0d0*pi
                    if(b .ge. a+0.5d0*pi)   stop 'Error: xyz_to_rtz fails, check.'
                elseif(b .le. a-0.5d0*pi) then
                    b   =  b+2.0d0*pi
                    if(b .ge. a-0.5d0*pi)   stop 'Error: xyz_to_rtz fails, check.'
                end if
            end if
            rtz(3*i-1)  =  b
            a           =  b
        end do
    else
        do i=1,N
            rtz(3*i-2)  =  xyz(1,i)*rad_axs(1)+xyz(2,i)*rad_axs(2)+xyz(3,i)*rad_axs(3)
            rtz(3*i-1)  =  xyz(1,i)*trs_axs(1)+xyz(2,i)*trs_axs(2)+xyz(3,i)*trs_axs(3)
            rtz(3*i  )  =  xyz(1,i)*rot_axs(1)+xyz(2,i)*rot_axs(2)+xyz(3,i)*rot_axs(3)
        end do
    end if

    return
    end subroutine xyz_to_rtz
