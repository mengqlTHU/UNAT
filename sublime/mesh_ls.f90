!-------------------------------------------------------------------------------
!   set least-square gradient for elements.
!-------------------------------------------------------------------------------
    subroutine set_ls_gra(lev)
    use var_kind_def
    use var_global, only: is_2d_cal,n_dim,err_mem
    use var_mesh
    use var_parallel
    use var_slv, only: gradient_method,gradient_GG,gradient_LS,gradient_LSm
    implicit none
    integer(dpI),intent(in):: lev
    integer(dpI):: isec,iele,i,j,k,m,n_e2n,stencil(2,200)
    real   (dpR):: c0(3),c(n_dim,100),ls_coe(n_dim,100)
    integer(dpI),allocatable:: e2n(:,:)

    if(gradient_method .le. gradient_GG )   return
    if(gradient_method .gt. gradient_LSm)   stop 'Error: gradient method not supported.'

!   ----------------------------------------------------------------------------
!   get the element-to-node information.
    n_e2n   =  0
    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        n_e2n   =  n_e2n+sec(isec)%n_ele*sec(isec)%npe
    end do
    allocate(e2n(3,n_e2n), stat=err_mem)
    n_e2n   =  0
    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
    do iele=1,sec(isec)%n_ele
    do j=1,sec(isec)%npe
        n_e2n   =  n_e2n+1
        e2n(1:3,n_e2n)  = (/sec(isec)%n2e(j,iele), isec, iele/)
    end do
    end do
    end do
    call iqsortcols(.true., 1, n_e2n, 1, 3, e2n)
    i   =  1
    do while(i .le. n_e2n)
        k   =  i
        do j=i+1,n_e2n
            if(e2n(1,j) .eq. e2n(1,i)) then
                k   =  j
            else
                exit
            end if
        end do

        do j=i+1,k
        do m=i,j-1
            if(e2n(1,m) .lt. 0) cycle
            if((e2n(2,j) .eq. e2n(2,m)) .and. (e2n(3,j) .eq. e2n(3,m))) then
                e2n(1:3,j)  = -abs(e2n(1:3,j))
                exit
            end if
        end do
        end do

        i   =  k+1
    end do
    j       =  n_e2n
    n_e2n   =  0
    do i=1,j
        if(e2n(1,i) .lt. 0) cycle
        n_e2n   =  n_e2n+1
        e2n(1:3,n_e2n)  =  e2n(1:3,i)
    end do
!   get the element-to-node information.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   estimate the size of jA_ls.
    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(.not. sec(isec)%is_int)  cycle

        allocate(sec(isec)%iA_ls(sec(isec)%n_ele+1), stat=err_mem)
        sec(isec)%iA_ls     =  0
        sec(isec)%iA_ls(1)  =  1
        do iele=1,sec(isec)%n_ele
            call set_LS_stencil(isec, iele, .true., i, stencil)
            sec(isec)%iA_ls(iele+1) =  sec(isec)%iA_ls(iele)+i
        end do
        i   =  sec(isec)%iA_ls(sec(isec)%n_ele+1)-1
        allocate(sec(isec)%jA_ls (2,i), stat=err_mem)
        if(err_mem .ne. 0)  stop 'Error: fails to allocate jA_ls.'

        if(gradient_method .eq. gradient_LS) then
!           store the coefficients.
            allocate(sec(isec)%coe_ls(3,i), stat=err_mem)
            if(err_mem .ne. 0)  stop 'Error: fails to allocate coe_ls.'
        elseif(gradient_method .eq. gradient_LSm) then
!           store the upper part of the LS matrix only.
            allocate(sec(isec)%coe_ls(6,sec(isec)%n_ele), stat=err_mem)
            if(err_mem .ne. 0)  stop 'Error: fails to allocate coe_ls.'
        end if
        sec(isec)%coe_ls=  0.0d0
    end do
!   estimate the size of jA_ls.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   compute coefficients.
    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(.not. sec(isec)%is_int)  cycle

        do iele=1,sec(isec)%n_ele
            call set_LS_stencil(isec, iele, .false., i, stencil)

            c0(1:n_dim) =  sec(isec)%cen(1:n_dim,iele)
            do j=1,i
                sec(isec)%jA_ls(1:2,j+sec(isec)%iA_ls(iele)-1)  =  stencil(1:2,j)
                c(1:n_dim,j)=  sec(stencil(1,j))%cen(1:n_dim,stencil(2,j))
            end do

            if(gradient_method .eq. gradient_LS) then
                if(is_2d_cal) then
                    call ls_coe_2d(.false., i, c0, c, ls_coe)
                else
                    call ls_coe_3d(.false., i, c0, c, ls_coe)
                end if
                k   =  sec(isec)%iA_ls(iele)-1
                do j=1,i
                    sec(isec)%coe_ls(1:n_dim,j+k)   =  ls_coe(1:n_dim,j)
                end do
            else
                if(is_2d_cal) then
                    call ls_coe_2d(.true., i, c0, c, ls_coe)
                    sec(isec)%coe_ls(1,iele)=  ls_coe(1,1)
                    sec(isec)%coe_ls(2,iele)=  ls_coe(2,1)
                    sec(isec)%coe_ls(3,iele)=  0.0d0
                    sec(isec)%coe_ls(4,iele)=  ls_coe(2,2)
                    sec(isec)%coe_ls(5,iele)=  0.0d0
                    sec(isec)%coe_ls(6,iele)=  0.0d0
                else
                    call ls_coe_3d(.true., i, c0, c, ls_coe)
                    sec(isec)%coe_ls(1,iele)=  ls_coe(1,1)
                    sec(isec)%coe_ls(2,iele)=  ls_coe(2,1)
                    sec(isec)%coe_ls(3,iele)=  ls_coe(3,1)
                    sec(isec)%coe_ls(4,iele)=  ls_coe(2,2)
                    sec(isec)%coe_ls(5,iele)=  ls_coe(3,2)
                    sec(isec)%coe_ls(6,iele)=  ls_coe(3,3)
                end if
            end if
        end do
    end do
!   compute coefficients.
!   ----------------------------------------------------------------------------

    if(allocated(e2n))  deallocate(e2n)

    return
    contains
!   ----------------------------------------------------------------------------
!   choose the least-square stecil.
!   ----------------------------------------------------------------------------
    subroutine set_LS_stencil(isec,iele,is_ini,n_stencil,stencil)
    use var_kind_def
    use var_mesh, only: mesh,sec
    use var_slv
    implicit none
    logical(dpL),intent(in):: is_ini
    integer(dpI),intent(in):: isec,iele
    logical(dpL):: ltmp
    integer(dpI):: n_stencil,stencil(2,*),c(2,300),n_c,i,j,k,s,e,im,L,R,v(8)

!   ----------------------------------------------------------------------------
!   the face-based stencil must be used.
    n_stencil   =  0
    do j=sec(isec)%iA_face_neighbour(iele),sec(isec)%iA_face_neighbour(iele+1)-1
        im  =  sec(isec)%jA_face_neighbour(j)
        if(im .gt. 0) then
            s   =  mesh(0)%mortar_LR(3, im)
            e   =  mesh(0)%mortar_LR(4, im)
        else
            s   =  mesh(0)%mortar_LR(1,-im)
            e   =  mesh(0)%mortar_LR(2,-im)
        end if
        n_stencil   =  n_stencil+1
        stencil(1,n_stencil)=  s
        stencil(2,n_stencil)=  e
    end do
    if(LS_stencil .eq. LS_stencil_face) return
!   the face-based stencil must be used.
!   ----------------------------------------------------------------------------

    v(1:sec(isec)%npe)  =  sec(isec)%n2e(1:sec(isec)%npe,iele)
    call iqsort(.true., 1, sec(isec)%npe, v)

    n_c =  0
    do i=1,sec(isec)%npe
        call ib_search(1, n_e2n, 1, 3, e2n, sec(isec)%n2e(i,iele), ltmp, L, R)
        if(.not. ltmp)  stop 'Error: fails to find LS stencil.'
        do j=L,R
            if((e2n(2,j) .eq. isec) .and. (e2n(3,j) .eq. iele)) cycle
            n_c         =  n_c+1
            c(1:2,n_c)  =  e2n(2:3,j)
        end do
    end do
    call simplify_matrix(2, n_c, 1, 2, c)

    if(LS_stencil .eq. LS_stencil_hybrid) then
        k   =  n_c
        n_c =  0
        do i=1,k
            s   =  c(1,i)
            e   =  c(2,i)
            m   =  0
            do j=1,sec(s)%npe
                if(sec(s)%n2e(j,e) .lt. v(1            ))   cycle
                if(sec(s)%n2e(j,e) .gt. v(sec(isec)%npe))   cycle
                call ib_search(1, sec(isec)%npe, 1, 1, v, sec(s)%n2e(j,e), ltmp, L, R)
                if(ltmp)    m   =  m+1
            end do
            if(m .ge. n_dim-1) then
                n_c         =  n_c+1
                c(1:2,n_c)  =  c(1:2,i)
            end if
        end do
    end if

    do i=1,n_c
        if((c(1,i) .eq. isec) .and. (c(2,i) .eq. iele)) cycle
        if(sec(c(1,i))%is_bnd)  cycle
        ltmp=  .false.
        do j=1,sec(isec)%iA_face_neighbour(iele+1)-sec(isec)%iA_face_neighbour(iele)
            if((c(1,i) .eq. stencil(1,j)) .and. (c(2,i) .eq. stencil(2,j))) then
                ltmp=  .true.
                exit
            end if
        end do
        if(ltmp)    cycle

        n_stencil               =  n_stencil+1
        stencil(1:2,n_stencil)  =  c(1:2,i)
    end do
    if(.not. is_ini)    call sort_matrix_col(2, 1, n_stencil, 1, 2, stencil)

    return
    end subroutine set_LS_stencil
    end subroutine set_ls_gra
!-------------------------------------------------------------------------------
!   cal the coe for least-square, 2D.
!-------------------------------------------------------------------------------
    subroutine ls_coe_2d(is_return_LHS,n,p0,p,rhs)
    use var_kind_def
    use var_slv, only: ls_weight_order
    implicit none
    logical(dpL),intent(in):: is_return_LHS
    integer(dpI),intent(in):: n
    real   (dpR),intent(in):: p0(*),p(2,*)
    logical(dpL):: weighted,ltmp(40)
    integer(dpI):: nn,i,j
    real   (dpR):: rhs(2,*),lhs(2,2),dp(2),w(42),pn(2,40),c(2,40)

!   ----------------------------------------------------------------------------
!   first delete redundant points.
    ltmp=  .false.
    do i=2,n
    do j=1,i-1
        if(ltmp(j)) cycle
        dp(1:2) =  abs(p(1:2,i)-p(1:2,j))
        if(dp(1) .le. 1.0d-12 .and. dp(2) .le. 1.0d-12) then
            ltmp(i) =  .true.
            exit
        end if
    end do
    end do
    nn  =  0
    do i=1,n
        if(.not. ltmp(i)) then
            nn  =  nn+1
            pn(1:2,nn)  =  p(1:2,i)
        else
            stop 'Error: duplicate LS stencil.'
        end if
    end do
!   first delete redundant points.
!   ----------------------------------------------------------------------------

    lhs =  0.0d0
!   call cal_ang_weights(nn, p0, pn, w)
    weighted=  ls_weight_order .ne. 0
    do i=1,nn
        dp(1:2) =  pn(1:2,i)-p0(1:2)
        if(weighted) then
            w(i)=  1.0d0/sqrt(dp(1)*dp(1)+dp(2)*dp(2))**ls_weight_order
        else
            w(i)=  1.0d0
        end if
        c(1:2,i)=  dp(1:2)*w(i)
        lhs(:,1)=  lhs(:,1)+dp(1)*c(:,i)
        lhs(:,2)=  lhs(:,2)+dp(2)*c(:,i)
    end do
    if(is_return_LHS) then
        call mat_inv(2, LHS)
        rhs(1:2,1)  =  LHS(1:2,1)
        rhs(1:2,2)  =  LHS(1:2,2)
        return
    end if
    call slv_axb(2, nn, lhs, c)

    nn  =  0
    do i=1,n
        if(ltmp(i)) then
            rhs(1:2,i)  =  0.0d0
        else
            nn  =  nn+1
            rhs(1:2,i)  =  c(1:2,nn)
        end if
    end do

    return
    end subroutine ls_coe_2d
!-------------------------------------------------------------------------------
!   cal the coe for least-square, 3D.
!-------------------------------------------------------------------------------
    subroutine ls_coe_3d(is_return_LHS,n,p0,p,rhs)
    use var_kind_def
    use var_slv, only: ls_weight_order
    implicit none
    logical(dpL),intent(in):: is_return_LHS
    integer(dpI),intent(in):: n
    real   (dpR),intent(in):: p0(*),p(3,*)
    logical(dpL):: weighted
    integer(dpI):: i
    real   (dpR):: rhs(3,*),dp(3),lhs(3,3),d(n)

    weighted=  ls_weight_order .ne. 0
    do i=1,n
        dp(1:3) =  p(1:3,i)-p0(1:3)
        if(weighted) then
            d(i)=  1.0d0/sqrt(dp(1)*dp(1)+dp(2)*dp(2)+dp(3)*dp(3))**ls_weight_order
        else
            d(i)=  1.0d0
        end if
    end do

    lhs =  0.0d0
    do i=1,n
        dp (1:3)=  p(1:3,i)-p0(1:3)
        rhs(:,i)=  dp(:)*d(i)

        lhs(:,1)=  lhs(:,1)+dp(1)*rhs(:,i)
        lhs(:,2)=  lhs(:,2)+dp(2)*rhs(:,i)
        lhs(:,3)=  lhs(:,3)+dp(3)*rhs(:,i)
    end do
    if(is_return_LHS) then
        call mat_inv(3, LHS)
        rhs(1:3,1)  =  LHS(1:3,1)
        rhs(1:3,2)  =  LHS(1:3,2)
        rhs(1:3,3)  =  LHS(1:3,3)
        return
    end if
    call slv_AxB(3, n, lhs, rhs)

    return
    end subroutine ls_coe_3d
