!-------------------------------------------------------------------------------
!   module GCROT.
!-------------------------------------------------------------------------------
    module var_gcrot
        use var_kind_def
        implicit none
        private
        public:: type_gcrot,gcrot_setup,gcrot_solution

        type type_gcrot
            logical(dpL):: is_preconditioned=  .true.
            integer(dpI):: max_iter         =  100
            integer(dpI):: n_vtx            =  0
            integer(dpI):: LDA              =  0
            integer(dpI):: bsize            =  1
            integer(dpI):: max_subspace     =  10
            integer(dpI):: max_augment      =  2
            real   (dpR):: restol           =  1.0d-10
            real   (dpR):: gcrot_ref(100)   =  1.0d0

            real   (dpR),allocatable:: C(:,:),U(:,:),H(:,:)
            real   (dpR),allocatable:: g(:),x(:),sn(:),cs(:),r(:),cnew(:),unew(:)
            real   (dpR),allocatable:: sol(:)
        end type type_gcrot

        contains
!       ------------------------------------------------------------------------
!       setup.
!       ------------------------------------------------------------------------
        subroutine gcrot_setup(g,bsize,n_vtx,max_iter,max_subspace,max_augment)
        implicit none
        integer(dpI),intent(in):: bsize,n_vtx,max_iter,max_subspace,max_augment
        integer(dpI):: err_mem
        type(type_gcrot):: g

        g%bsize     =  bsize
        g%n_vtx     =  n_vtx
        g%LDA       =  n_vtx*bsize
        if(max_iter     .gt. 0) g%max_iter      =  max_iter
        if(max_subspace .gt. 0) g%max_subspace  =  max_subspace
        if(max_augment  .gt. 0) g%max_augment   =  max_augment

        allocate(g%C   (g%LDA           , g%max_subspace+1), stat=err_mem)
        allocate(g%U   (g%LDA           , g%max_subspace  ), stat=err_mem)
        allocate(g%H   (g%max_subspace+1, g%max_subspace  ), stat=err_mem)
        allocate(g%g   (g%max_subspace+1                  ), stat=err_mem)
        allocate(g%x   (g%max_subspace+1                  ), stat=err_mem)
        allocate(g%sn  (g%max_subspace+1                  ), stat=err_mem)
        allocate(g%cs  (g%max_subspace+1                  ), stat=err_mem)
        allocate(g%r   (g%LDA                             ), stat=err_mem)
        allocate(g%cnew(g%LDA                             ), stat=err_mem)
        allocate(g%unew(g%LDA                             ), stat=err_mem)
        allocate(g%sol (g%LDA                             ), stat=err_mem)
        g%C     =  0.0d0
        g%U     =  0.0d0
        g%H     =  0.0d0
        g%g     =  0.0d0
        g%x     =  0.0d0
        g%sn    =  0.0d0
        g%cs    =  0.0d0
        g%r     =  0.0d0
        g%cnew  =  0.0d0
        g%unew  =  0.0d0
        g%sol   =  0.0d0

        return
        end subroutine gcrot_setup
!       ------------------------------------------------------------------------
!       gcrot solution.
!       ------------------------------------------------------------------------
        subroutine gcrot_solution(g,set_b,matvec,precond)
        use var_kind_def
        implicit none
        integer(dpI):: nkryv,jmax,kmax,outdone,indone,iter,i,j,miniter,i2,ptr
        real   (dpR):: beta,norm0,alpha,tmp
        type(type_gcrot):: g
        real   (dpR),external:: DNRM2,DDOT
        external:: set_b,matvec,precond

        miniter =  1
        nkryv   =  0
        call set_b(g%r)

        beta    =  DNRM2(g%LDA, g%r, 1)
        norm0   =  beta
        jmax    =  g%max_subspace
        kmax    =  0
        outdone =  0
        do iter=1,g%max_iter
            tmp =  1.0d0/beta
            g%C(:,kmax+1)   =  g%r(:)*tmp

            g%H     =  0.0d0
            g%g     =  0.0d0
            g%g(1)  =  beta
            indone  =  0

            do j=1,jmax
                nkryv   =  nkryv+1
                if(g%is_preconditioned) then
                    call precond(g%LDA, g%C(1,kmax+j), g%U(1,kmax+j))
                else
                    call DCOPY(g%LDA, g%C(1,kmax+j), 1, g%U(1,kmax+j), 1)
                end if

                call matvec(g%U(1,kmax+j), g%C(1,kmax+j+1))
                call mgs(kmax+j+1, g%LDA, g%H(1,j), g%C)

                do i=1,j-1
                    call applyGivens(g%sn(i), g%cs(i), g%H(kmax+i,j), g%H(kmax+i+1,j))
                end do
                call generateGivens(g%H(kmax+j,j), g%H(kmax+j+1,j), g%sn(j), g%cs(j))
                call applyGivens(g%sn(j), g%cs(j), g%g(j), g%g(j+1))
                beta=  abs(g%g(j+1))

                if(((beta < norm0*g%restol) .and. (j>=miniter)) .or. (j .eq. jmax) .or. &
                   &(nkryv .ge. g%max_iter))    indone  =  1
                if(indone .eq. 1)   exit
            end do

            g%x(1:j)=  g%g(1:j)
            do i=j,1,-1
                g%x(i)  =  g%x(i)/g%H(kmax+i, i)
                do i2=i-1,1,-1
                    g%x(i2) =  g%x(i2)-g%H(kmax+i2,i)*g%x(i)
                end do
            end do
            g%unew  =  0.0d0
            do i=1,j
                g%unew  =  g%unew+g%x(i)*g%U(:,kmax+i)
            end do
            do i2=1,kmax
                tmp =  0.0d0
                do i=1,j
                    tmp =  tmp+g%H(i2,i)*g%x(i)
                end do
                g%unew  =  g%unew-tmp*g%U(:,i2)
            end do

            g%x     =  g%g
            g%x(j+1)=  0.0d0
            do i=j,1,-1
                call applyGivens(-g%sn(i), g%cs(i), g%x(i), g%x(i+1))
            end do
            g%cnew  =  0.0d0
            do i=1,j+1
                g%cnew  =  g%cnew+g%x(i)*g%C(:,kmax+i)
            end do

            alpha   =  1.0d0/DNRM2(g%LDA, g%cnew, 1)
            g%cnew  =  alpha*g%cnew
            g%unew  =  alpha*g%unew
            alpha   =  DDOT(g%LDA, g%cnew, 1, g%r, 1)
            g%r     =  g%r-alpha*g%cnew
            g%sol   =  g%sol+alpha*g%unew
            beta    =  DNRM2(g%LDA, g%r, 1)

            if(((beta .lt. norm0*g%restol) .and. (iter .ge. miniter)) .or. &
               &(nkryv .ge. g%max_iter))    outdone  =  1
            print*,iter,beta/norm0
            if(outdone .eq. 1)  exit

            if(kmax .lt. g%max_augment) then
                kmax=  kmax+1
                jmax=  jmax-1
                ptr =  kmax
            else
                if(g%max_augment .gt. 0)    ptr =  mod(ptr, kmax)+1
            end if
            if(g%max_augment .eq. 0)    cycle
            call DCOPY(g%LDA, g%cnew(1), 1, g%C(1,ptr), 1)
            call DCOPY(g%LDA, g%unew(1), 1, g%U(1,ptr), 1)
        end do

        return
        end subroutine gcrot_solution
!       ------------------------------------------------------------------------
!       modified Gram-Schmidt procedure.
!       ------------------------------------------------------------------------
        subroutine mgs(i,LDA,H,w)
        implicit none
        integer(dpI),intent(in):: i,LDA
        integer(dpI):: k
        real   (dpR):: H(*),w(LDA,*),reorth,nrm0,thr,fct,nrm1
        real   (dpR),external:: DDOT,DNRM2

        reorth  =  0.98d0
        nrm0    =  DNRM2(LDA, w(1,i), 1)**2
        thr     =  nrm0*reorth

        do k=1,i-1
            fct     =  DDOT(LDA, w(1,i), 1, w(1,k), 1)
            H(k)    =  fct
            w(:,i)  =  w(:,i)-fct*w(:,k)

            if(fct*fct .gt. thr) then
                fct =  DDOT(LDA, w(1,i), 1, w(1,k), 1)
                H(k)=  H(k)+fct
                w(:,i)  =  w(:,i)-fct*w(:,k)
            end if

            nrm0=  nrm0-H(k)**2
            if(nrm0 .lt. 0.0d0) nrm0=  0.0d0

            thr =  nrm0*reorth
        end do

        nrm1    =  DNRM2(LDA, w(1,i), 1)
        H(i)    =  nrm1
        fct     =  1.0d0/nrm1
        w(:,i)  =  w(:,i)*fct

        return
        end subroutine mgs
!       ------------------------------------------------------------------------
!       Givens rotation.
!       ------------------------------------------------------------------------
        subroutine applyGivens(s,c,h1,h2)
        implicit none
        real(dpR),intent(in):: s,c
        real(dpR):: h1,h2,temp

        temp=  c*h1+s*h2
        h2  =  c*h2-s*h1
        h1  =  temp

        return
        end subroutine applyGivens
!       ------------------------------------------------------------------------
!       generate Givens rotation.
!       ------------------------------------------------------------------------
        subroutine generateGivens(dx,dy,s,c)
        implicit none
        real(dpR):: dx,dy,s,c,t

        if((dx .eq. 0.0d0) .and. (dy .eq. 0.0d0)) then
            c   =  1.0d0
            s   =  1.0d0
        elseif(abs(dy) .gt. abs(dx)) then
            t   =  dx/dy
            dx  =  sqrt(1.0d0+t*t)
            s   = (1.0d0/dx)*(dy/abs(dy))
            c   =  t*s
        elseif(abs(dy) .le. abs(dx)) then
            t   =  dy/dx
            dy  =  sqrt(1.0d0+ t*t)
            c   = (1.0d0/dy)*(dx/abs(dx))
            s   =  t*c
        else
            dx  = 0.0d0
            dy  = 0.0d0
            c   = 1.0d0
            s   = 0.0d0
        end if
        dx  =  abs(dx*dy)
        dy  =  0.0d0

        return
        end subroutine generateGivens
    end module var_gcrot
!-------------------------------------------------------------------------------
!   module GMRES.
!-------------------------------------------------------------------------------
    module var_gmres
        use var_kind_def
        implicit none
        private
        public:: type_gmres,gmres_set,machine_tiny,gmres_solution

        real   (dpR),parameter:: machine_tiny   =  tiny(1.0d0)

        type type_gmres
            logical(dpL):: is_preconditioned=  .true.
            logical(dpL):: is_scaled        =  .false.
            logical(dpL):: is_allocated     =  .false.
            logical(dpL):: is_output_res    =  .true.
            integer(dpI):: max_kry          =  10
            integer(dpI):: max_kry_rst      =  1
            real   (dpR):: TOL_gmres        =  1.0d-10
            real   (dpR):: ref(100)         =  1.0d0

            integer(dpI):: n_vtx=  0
            integer(dpI):: LDA  =  0
            integer(dpI):: bsize=  1
            real   (dpR):: bnorm,relerr,rnorm
            real   (dpR),pointer:: H_k(:,:),s_k(:)
            real   (dpR),allocatable:: sn_k(:),cs_k(:),y_k(:)
            real   (dpR),allocatable:: Rhs(:),b(:)

            real   (dpR),pointer:: window(:)
            real   (dpR),pointer:: v_k(:,:),R_k(:)
        end type type_gmres
        contains
!       ------------------------------------------------------------------------
!       compute sin and cos.
!       ------------------------------------------------------------------------
        subroutine Gmresm_Get_rotation(vector_in,cos_theta,sin_theta)
        implicit none
        real(dpR),intent(in)  :: vector_in(2)
        real(dpR),intent(out) :: cos_theta,sin_theta
        real(dpR):: temp

        if ( abs(vector_in (2)) .lt. machine_tiny ) then
            cos_theta = 1.0d0
            sin_theta = 0.0d0
        elseif ( abs ( vector_in(2) ) .gt. abs ( vector_in (1) ) ) then
            temp        = -vector_in(1)/vector_in(2)
            sin_theta   =  1.0d0/sqrt(1.0d0+temp*temp)
            cos_theta   =  temp*sin_theta
        else
            temp        = -vector_in(2)/vector_in(1)
            cos_theta   =  1.0d0/sqrt(1.0d0+temp*temp)
            sin_theta   =  temp*cos_theta
        end if

        return
        end subroutine Gmresm_Get_rotation
!       ------------------------------------------------------------------------
!       rotate vector.
!       ------------------------------------------------------------------------
        function Gmresm_Rotate_vector ( vec_in, cos_theta, sin_theta )
        implicit none
        real(dpR) , intent (in) :: vec_in (2), cos_theta, sin_theta
        real(dpR)               :: Gmresm_Rotate_vector (2)
        real(dpR) :: temp

        temp                    = cos_theta * vec_in (1) - sin_theta * vec_in (2)
        Gmresm_Rotate_vector(2) = sin_theta * vec_in (1) + cos_theta * vec_in (2)
        Gmresm_Rotate_vector(1) = temp

        return
        end function Gmresm_Rotate_vector
!       ------------------------------------------------------------------------
!       solve the upper triangle matrix.
!       ------------------------------------------------------------------------
        function Gmresm_Solve_upper_triang (a, b_rhs, n)
        implicit none
        integer(dpI),intent(in):: n
        real   (dpR),intent(in):: a(:,:),b_rhs(:)
        integer(dpI):: j
        real   (dpR):: Gmresm_Solve_upper_triang (n)

        if (ubound(a,DIM=1) /= ubound (a,DIM=2) + 1) stop 'PROBLEM 1 in SOLVETR'
        if (n > ubound (a,DIM = 2) .OR. n < 1) stop 'PROBLEM 2 in SOLVETR'

        Gmresm_Solve_upper_triang = b_rhs(1:n)
        do j = N, 1, -1
            Gmresm_Solve_upper_triang (j)     = Gmresm_Solve_upper_triang (j) / a (j,j)
            Gmresm_Solve_upper_triang (1:j-1) = Gmresm_Solve_upper_triang (1:j-1) - &
                                      Gmresm_Solve_upper_triang (j) * a(1:j-1,j)
        end do

        return
        end function Gmresm_Solve_upper_triang
!       ------------------------------------------------------------------------
!       set GMRES solver.
!       ------------------------------------------------------------------------
        subroutine gmres_set(s,n_vtx,bsize,max_subspace,max_restart,ref,cvg)
        implicit none
        integer(dpI),intent(in):: n_vtx,bsize,max_subspace,max_restart
        real   (dpR),intent(in):: ref(*),cvg
        integer(dpI):: err_mem
        real   (dpR):: a,b
        type(type_gmres):: s

        if(s%is_allocated)  return
        s%n_vtx         =  n_vtx
        s%bsize         =  bsize
        s%LDA           =  s%n_vtx*s%bsize
        s%max_kry       =  min(s%LDA, max_subspace)
        s%max_kry_rst   =  max_restart
        if(cvg .gt. 0.0d0)  s%TOL_gmres =  cvg

        allocate(s%v_k(s%LDA,0:s%max_kry), stat=err_mem)
        if(err_mem .ne. 0)  stop 'Error: GMRES_set fails to allocate memory.'
        allocate(s%R_k(s%LDA), stat=err_mem)
        allocate(s%b  (s%LDA), stat=err_mem)
        s%v_k   =  0.0d0
        allocate(s%H_k (s%max_kry+1,s%max_kry), stat=err_mem)
        allocate(s%s_k (s%max_kry+1), stat=err_mem)
        allocate(s%sn_k(s%max_kry), stat=err_mem)
        allocate(s%cs_k(s%max_kry), stat=err_mem)
        allocate(s%y_k (s%max_kry), stat=err_mem)
        a   =  maxval(ref(1:s%bsize))
        b   =  minval(ref(1:s%bsize))
        if((abs(a-b) .gt. 1.0d-10*a) .and. (s%bsize .gt. 1)) then
            s%is_scaled     =  .true.
            s%ref(1:s%bsize)=  ref(1:s%bsize)
        end if
        s%is_allocated  =  .true.

        return
        end subroutine gmres_set
!       ------------------------------------------------------------------------
!       GMRES solution.
!       ------------------------------------------------------------------------
        subroutine gmres_solution(s,set_b,matvec,precond)
        use var_parallel
        implicit none
        integer(dpI):: ikry,i,L,R,iter
        real   (dpR):: x1(100)
        real   (dpR),external:: DNRM_mpi,DDOT_mpi
        type   (type_gmres):: s
        external:: set_b,matvec,precond

        x1(1:s%bsize)   =  1.0d0/s%ref(1:s%bsize)
        call set_b(s%LDA, s%b)
        s%v_k(:,0)  =  0.0d0
        s%H_k       =  0.0d0
        rst_ite:do iter=1,s%max_kry_rst
!           R_k =  b-A*x0
            if(iter .gt. 1) then
                call matvec(s%is_scaled, s%n_vtx, s%bsize, s%ref, s%v_k, s%R_k)
                s%R_k   =  s%b-s%R_k
            else
                s%R_k   =  s%b
            end if

!           R_k =  preconditioner(R_k)
            if(s%is_preconditioned) call precond(s%LDA, s%R_k)
            if(s%is_scaled) then
                do i=1,s%n_vtx
                    L   =  s%bsize*(i-1)+1
                    R   =  s%bsize* i
                    s%R_k(L:R)  =  s%R_k(L:R)*x1(1:s%bsize)
                end do
            end if

            s%rnorm =  DNRM_mpi(s%LDA, s%R_k)
            if(iter .eq. 1) s%bnorm =  s%rnorm

!           s%v_k(1:s%LDA,1)    =  s%R_k(1:s%LDA)/s%rnorm
            call DAXY(s%LDA, 1.0d0/s%rnorm, s%R_k, s%v_k(:,1))
            s%s_k(1)            =  s%rnorm
            s%s_k(2:s%max_kry+1)=  0.0d0

            inner_ite:do ikry=1,s%max_kry
!               R_k =  A*v_k(:,ikry)
                call matvec(s%is_scaled, s%n_vtx, s%bsize, s%ref, s%v_k(1,ikry), s%R_k)

!               R_k =  preconditioner(R_k)
                if(s%is_preconditioned) call precond(s%LDA, s%R_k)
                if(s%is_scaled) then
                    do i=1,s%n_vtx
                        L   =  s%bsize*(i-1)+1
                        R   =  s%bsize* i
                        s%R_k(L:R)  =  s%R_k(L:R)*x1(1:s%bsize)
                    end do
                end if

                do i=1,ikry
                    s%H_k(i,ikry)   =  DDOT_mpi(s%LDA, s%R_k, s%v_k(:,i))
                    call DAXPY(s%LDA, -s%H_k(i,ikry), s%v_k(:,i), 1, s%R_k, 1)
                end do

                s%H_k(ikry+1,ikry)  =  DNRM_mpi(s%LDA, s%R_k)
                if(ikry .lt. s%max_kry) call DAXY(s%LDA, 1.0d0/s%H_k(ikry+1,ikry), &
                    &  s%R_k, s%v_k(:,ikry+1))

                do i=1,ikry-1
                    s%window=> s%H_k(i:i+1,ikry)
                    s%window=  Gmresm_Rotate_vector(s%window,s%cs_k(i),s%sn_k(i))
                end do

                s%window=> s%H_k(ikry:ikry+1,ikry)
                call Gmresm_Get_rotation(s%window, s%cs_k(ikry), s%sn_k(ikry))

                s%window=> s%H_k(ikry:ikry+1,ikry)
                s%window=  Gmresm_Rotate_vector(s%window, s%cs_k(ikry), s%sn_k(ikry))
                s%H_k(ikry+1,ikry)= 0.0d0
     
                s%window=> s%s_k(ikry:ikry+1)
                s%window=  Gmresm_Rotate_vector(s%window, s%cs_k(ikry), s%sn_k(ikry))

                s%relerr=  abs(s%s_k(ikry+1))/s%bnorm

                if(s%relerr .le. s%TOL_gmres) then
                    s%y_k(1:ikry)   =  Gmresm_Solve_upper_triang(s%H_k,s%s_k,ikry)
                    do i=1,ikry
                        call DAXPY(s%LDA, s%y_k(i), s%v_k(:,i), 1, s%v_k, 1)
                    end do
                    if((myid .eq. 0) .and. s%is_output_res) &
                        &  write(unit=6,fmt='(I8,ES12.4)'),iter,s%relerr
                    exit rst_ite
                end if
            end do inner_ite
            if((myid .eq. 0) .and. s%is_output_res) &
                &  write(unit=6,fmt='(I8,ES12.4)'),iter,s%relerr

            s%y_k(1:s%max_kry)  =  Gmresm_Solve_upper_triang(s%H_k,s%s_k,s%max_kry)
            do i=1,s%max_kry
                call DAXPY(s%LDA, s%y_k(i), s%v_k(:,i), 1, s%v_k, 1)
            end do
        end do rst_ite

        return
        end subroutine gmres_solution
    end module var_gmres
!-------------------------------------------------------------------------------
!   LOOSE GMRES.
!-------------------------------------------------------------------------------
    module var_lgmres
        use var_kind_def
        implicit none
        private
        public:: type_lgmres,lgmres_set,lgmres_solution

        type type_lgmres
            logical(dpL):: is_preconditioned=  .true.
            logical(dpL):: is_scaled        =  .false.
            integer(dpI):: n_vtx            =  0
            integer(dpI):: bsize            =  1
            integer(dpI):: LDA              =  0
            integer(dpI):: m                =  10
            integer(dpI):: k                =  1
            integer(dpI):: max_restart      =  100
            real   (dpR):: restol           =  1.0d-10
            real   (dpR):: ref(100)         =  1.0d0
            real   (dpR),allocatable:: b(:),x(:),v(:,:),z(:,:),H(:,:)
            real   (dpR),allocatable:: residual(:)
        end type type_lgmres
        contains
!       ------------------------------------------------------------------------
!       set LGMRES solver.
!       ------------------------------------------------------------------------
        subroutine lgmres_set(g,n_vtx,bsize,max_subspace,max_augment,max_restart,ref)
        implicit none
        integer(dpI),intent(in):: n_vtx,bsize,max_subspace,max_augment,max_restart
        real   (dpR),intent(in):: ref(*)
        integer(dpI):: err_mem
        real   (dpR):: a,b
        type(type_lgmres):: g

        g%n_vtx =  n_vtx
        g%bsize =  bsize
        g%LDA   =  n_vtx*bsize
        if(max_subspace .gt. 0) g%m =  min(max_subspace, g%LDA)
        if(max_augment  .ge. 0) g%k =  max_augment
        if(g%m .eq. g%LDA)      g%k =  0
        if(max_restart  .gt. 0) g%max_restart   =  max_restart
        allocate(g%b(g%LDA               ), stat=err_mem)
        allocate(g%x(g%LDA               ), stat=err_mem)
        allocate(g%v(g%LDA    , g%m+g%k+1), stat=err_mem)
        allocate(g%H(g%m+g%k+1, g%m+g%k  ), stat=err_mem)
        allocate(g%residual(g%max_restart), stat=err_mem)
        if(g%k .gt. 0)  allocate(g%z(g%LDA, g%k), stat=err_mem)
        g%b =  0.0d0
        g%x =  0.0d0
        g%v =  0.0d0
        g%H =  0.0d0
        g%residual  =  0.0d0

        a   =  maxval(ref(1:g%bsize))
        b   =  minval(ref(1:g%bsize))
        if((abs(a-b) .gt. 1.0d-10*a) .and. (g%bsize .gt. 1)) then
            g%is_scaled     =  .true.
            g%ref(1:g%bsize)=  ref(1:g%bsize)
        end if

        return
        end subroutine lgmres_set
!       ------------------------------------------------------------------------
!       LGMRES solution.
!       ------------------------------------------------------------------------
        subroutine lgmres_solution(g,set_b,matvec,precond)
        use var_parallel, only: myid
        implicit none
        integer(dpI):: iter,i,j,l,s,mL,r,err_mem,M,N
        real   (dpR):: beta,res,x1(100)
        type   (type_lgmres):: g
        real   (dpR),external:: DNRM2,DDOT,DNRM_mpi,DDOT_mpi
        real   (dpR),allocatable:: u(:),y(:)
        external:: set_b,matvec,precond

        allocate(u(g%LDA      ), stat=err_mem)
        allocate(y(g%m  +g%k+1), stat=err_mem)
        u   =  0.0d0
        y   =  0.0d0
        x1(1:g%bsize)   =  1.0d0/g%ref(1:g%bsize)
        call set_b(g%LDA, g%b)

        iter=  0
        do r=0,g%max_restart-1
!           b-A*x0
            call matvec(g%is_scaled, g%n_vtx, g%bsize, g%ref, g%x, u)
            u(1:g%LDA)      =  g%b(1:g%LDA)-u(1:g%LDA)
            if(g%is_preconditioned)     call precond(g%LDA, u)
            if(g%is_scaled) then
                do i=1,g%n_vtx
                    M   =  g%bsize*(i-1)+1
                    N   =  g%bsize* i
                    u(M:N)  =  u(M:N)*x1(1:g%bsize)
                end do
            end if

!           beta            =  DNRM2(g%LDA, u, 1)
            beta            =  DNRM_mpi(g%LDA, u)
            if(beta .le. g%restol*g%residual(1))    exit
            iter            =  iter+1
            g%residual(iter)=  beta
            if(myid .eq. 0) write(unit=6,fmt='(I8,ES12.4)'),iter,beta

            g%v(1:g%LDA,1)  =  u(1:g%LDA)/beta
            s               =  g%m+g%k
            mL              =  g%m+max(g%k-r, 0)

            do j=1,s
                if(j .le. mL) then
                    call matvec(g%is_scaled, g%n_vtx, g%bsize, g%ref, g%v(1,j   ), u)
                else
                    call matvec(g%is_scaled, g%n_vtx, g%bsize, g%ref, g%z(1,j-mL), u)
                end if
                if(g%is_preconditioned) call precond(g%LDA, u)
                if(g%is_scaled) then
                    do i=1,g%n_vtx
                        M   =  g%bsize*(i-1)+1
                        N   =  g%bsize* i
                        u(M:N)  =  u(M:N)*x1(1:g%bsize)
                    end do
                end if

                do l=1,j
!                   g%H(l,j)=  DDOT(g%LDA, u, 1, g%v(1,l), 1)
                    g%H(l,j)=  DDOT_mpi(g%LDA, u, g%v(1,l))
                    call DAXPY(g%LDA, -g%H(l,j), g%v(1,l), 1, u, 1)
                end do
!               g%H(j+1,j)  =  DNRM2(g%LDA, u, 1)
                g%H(j+1,j)  =  DNRM_mpi(g%LDA, u)
                g%v(1:g%LDA, j+1)   =  u(1:g%LDA)/g%H(j+1,j)
            end do
            call least_square(s+1, s, g%H, beta, y, res)

            u(1:g%LDA)  =  g%v(1:g%LDA,1)*y(1)
            do i=2,g%m
                call DAXPY(g%LDA, y(i    ), g%v(1,i), 1, u, 1)
            end do
            do i=1,g%k
                call DAXPY(g%LDA, y(i+g%M), g%z(1,i), 1, u, 1)
            end do
            do i=g%k-1,1,-1
                call DCOPY(g%LDA, g%z(1,i), 1, g%z(1,i+1), 1)
            end do
            if(g%k .gt. 0)  call DCOPY(g%LDA, u, 1, g%z, 1)
            call DAXPY(g%LDA, 1.0d0, u, 1, g%x, 1)
        end do
!       print*,iter,beta/g%residual(1),res
        if(allocated(u))    deallocate(u)
        if(allocated(y))    deallocate(y)

        return
        end subroutine lgmres_solution
!       ------------------------------------------------------------------------
!       generate Givens rotation.
!       ------------------------------------------------------------------------
        subroutine generateGivens(dx,dy,c,s)
        implicit none
        real(dpR),intent(in):: dx,dy
        real(dpR):: s,c,t

        if((dx .eq. 0.0d0) .and. (dy .eq. 0.0d0)) then
            c   =  1.0d0
            s   =  1.0d0
        elseif(abs(dy) .gt. abs(dx)) then
            t   =  dx/dy
            s   = (1.0d0/sqrt(1.0d0+t*t))*(dy/abs(dy))
            c   =  t*s
        elseif(abs(dy) .le. abs(dx)) then
            t   =  dy/dx
            c   = (1.0d0/sqrt(1.0d0+t*t))*(dx/abs(dx))
            s   =  t*c
        else
            c   = 1.0d0
            s   = 0.0d0
        end if

        return
        end subroutine generateGivens
!       ------------------------------------------------------------------------
!       solve the least-square problem, H(M,N)*x-beta*e1.
!       ------------------------------------------------------------------------
        subroutine least_square(M,N,H,beta,x,res)
        use var_kind_def
        implicit none
        integer(dpI),intent(in):: M,N
        real   (dpR),intent(in):: beta
        integer(dpI):: i,j
        real   (dpR):: H(M,*),x(*),res,g(M),a,b,c,s

        g(1  )  =  beta
        g(2:M)  =  0.0d0
        do j=1,N
            call generateGivens(H(j,j), H(j+1,j), c, s)
            do i=j,N
                a   =  c*H(j,i)+s*H(j+1,i)
                b   = -s*H(j,i)+c*H(j+1,i)
                H(j  ,i)=  a
                H(j+1,i)=  b
            end do
            a   =  c*g(j)+s*g(j+1)
            b   = -s*g(j)+c*g(j+1)
            g(j  )  =  a
            g(j+1)  =  b
        end do
        x(1:N)  =  g(1:N)
        do j=N,1,-1
            x(  j  )=  x(j)/H(j,j)
            x(1:j-1)=  x(1:j-1)-x(j)*H(1:j-1,j)
        end do
        res =  abs(g(M))

        return
        end subroutine least_square
    end module var_lgmres
!-------------------------------------------------------------------------------
!   Invert matrix by Gauss method with pivot.
!-------------------------------------------------------------------------------
    subroutine mat_inv(n,a)
    use var_kind_def
    implicit none
    integer(dpI),intent(in):: n
    integer(dpI):: j
    real   (dpR):: a(n,*),I(n,n)

    I   =  0.0d0
    forall(j=1:N)   I(j,j)  =  1.0d0
    call slv_AxB(N, N, A, I)
    A(1:N,1:N)  =  I(1:N,1:N)

    return
    end subroutine mat_inv
!-------------------------------------------------------------------------------
!   slv Ax=b using Gauss Elimination with scaling and pivoting.
!-------------------------------------------------------------------------------
    subroutine slv_AxB(LDA,NRHS,a,b)
    use var_kind_def
    implicit none 
    integer(dpI),intent(in):: LDA,NRHS
    integer(dpI):: i,j,k,l
    real   (dpR):: a(LDA,*),b(LDA,*)
    real   (dpR):: s(LDA),c(NRHS),pivot,rtmp

!   step 1: begin forward elimination
    do k=1,LDA-1
!       step 2: "scaling", s(i) will have the largest element from row i 
        do i=k,LDA
            s(i)=  0.0d0
            do j=k,LDA
                s(i)=  max(s(i),abs(a(i,j)))
            end do
        end do

!       step 3: "pivoting 1", find a row with the largest pivoting element
        pivot   =  abs(a(k,k)/s(k))
        l       =  k
        do j=k+1,LDA
            if(abs(a(j,k)/s(j)) .gt. pivot) then
                pivot   =  abs(a(j,k)/s(j))
                l       =  j
            end if
        end do
        if(pivot .eq. 0.0d0)    stop 'Error: the matrix is singular.'

!       step 4: "pivoting 2" interchange rows k and l (if needed)
        if(l .ne. k) then
            do j=k,LDA
                rtmp    =  a(k,j)
                a(k,j)  =  a(l,j)
                a(l,j)  =  rtmp
            end do
            c(  1:NRHS) =  b(k,1:NRHS)
            b(k,1:NRHS) =  b(l,1:NRHS)
            b(l,1:NRHS) =  c(  1:NRHS)
        end if

!       step 5: the elimination (after scaling and pivoting)
        do i=k+1,LDA
            rtmp    =  a(i,k)/a(k,k)
            a(i,k)  =  0.0d0
            b(i,1:NRHS) =  b(i,1:NRHS)-rtmp*b(k,1:NRHS)
            do j=k+1,LDA
                a(i,j)  =  a(i,j)-rtmp*a(k,j)
            end do
        end do
    end do

!   step 6: back substiturion 
    b(LDA,1:NRHS)   =  b(LDA,1:NRHS)/a(LDA,LDA)
    do i=LDA-1,1,-1
        c   =  a(i,i+1)*b(i+1,1:NRHS)
        do j=i+2,LDA
            c   =  c+a(i,j)*b(j,1:NRHS)
        end do 
        b(i,1:NRHS)  = (b(i,1:NRHS)-c)/a(i,i)
    end do

    return
    end subroutine slv_AxB
!-------------------------------------------------------------------------------
!   add an element to sparse matrix in CSR format.
!-------------------------------------------------------------------------------
    subroutine add_ele_to_csr(is_trans,LDA,row,col,ele,iA,jA,A)
    use var_kind_def
    implicit none
    logical(dpL),intent(in):: is_trans
    integer(dpI),intent(in):: LDA,row,col,iA(*),jA(*)
    real   (dpR),intent(in):: ele(LDA,*)
    integer(dpI):: i,j,k
    real   (dpR):: A(LDA,*)

    do k=iA(row),iA(row+1)-1
        if(jA(k) .eq. col) then
            if(is_trans) then
                do j=1,LDA
                do i=1,LDA
                    A(i,LDA*(k-1)+j)=  A(i,LDA*(k-1)+j)+ele(j,i)
                end do
                end do
            else
                do j=1,LDA
                do i=1,LDA
                    A(i,LDA*(k-1)+j)=  A(i,LDA*(k-1)+j)+ele(i,j)
                end do
                end do
            end if
            return
        end if
    end do
    print*,'Error: failed to find A(',row,',',col,') in the CSR.'
    stop

    return
    end subroutine add_ele_to_csr
!-------------------------------------------------------------------------------
!   Existence tesing of A(i,j).
!-------------------------------------------------------------------------------
    subroutine existelm(iA,jA,i,j,exist_ij,ind_ij)
    use var_kind_def
    implicit none
    logical(dpL):: exist_ij
    integer(dpI),intent(in):: iA(*),jA(*),i,j
    integer(dpI):: ind_ij,ind

    exist_ij=  .false.
    do ind=iA(i),iA(i+1)-1
        if(jA(ind) .eq. j) then
            exist_ij    =  .true.
            ind_ij      =  ind
            return
        end if
    end do

    return
    end subroutine existelm
!-------------------------------------------------------------------------------
!   Solving LUx=b for BILU0.
!-------------------------------------------------------------------------------
    subroutine slv_x_BILU0(dimA,iA,jA,DA,bsize,A,x)
    use var_kind_def
    implicit none
    integer(dpI),intent(in):: dimA,bsize,iA(*),jA(*),DA(*)
    real   (dpR),intent(in):: A(bsize,*)
    real   (dpR):: x(bsize,*)
    integer(dpI):: i,ind,col,inc
    real   (dpR):: Rhs(bsize)

    inc =  bsize-1

!   ----------------------------------------------------------------------------
!   Forward iteration.
    do i=1,dimA
        Rhs(:)  =  x(:,i)
        do ind=iA(i),DA(i)-1
            col =  bsize*(ind-1)+1
            call DGEMV('N', bsize, bsize, -1.0d0, A(1,col), bsize, x(1,jA(ind)), &
                &  1, 1.0d0, Rhs, 1)
        end do
        x(:,i)  =  Rhs(:)
    end do
!   Forward iteration.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   Backward iteration.
    do i=dimA,1,-1
        Rhs(:)  =  x(:,i)
        do ind=DA(i)+1,iA(i+1)-1
            col =  bsize*(ind-1)+1
            call DGEMV('N', bsize, bsize, -1.0d0, A(1,col), bsize, x(1,jA(ind)), &
                &  1, 1.0d0, Rhs, 1)
        end do
        col =  bsize*(DA(i)-1)+1
        call DGEMV('N', bsize, bsize, 1.0d0, A(1,col), bsize, Rhs, 1, 0.0d0, x(1,i), 1)
    end do
!   Backward iteration.
!   ----------------------------------------------------------------------------

    return
    end subroutine slv_x_BILU0
!-------------------------------------------------------------------------------
!   get the eigenvalues and eigenvectors of a real-symmetric matrix.
!-------------------------------------------------------------------------------
    subroutine dsyev(LDA,A,e,v)
    use var_kind_def
    implicit none
    integer(dpI),intent(in):: LDA
    real   (dpR),intent(in):: A(*)
    integer(dpI):: ierr
    real   (dpR):: e(*),v(*),rtmv(LDA)

    call TRED2(LDA, LDA, A, e, rtmv, v)
    call TQL2 (LDA, LDA, e, rtmv, v, ierr)
    if(ierr .ne. 0) stop 'Error: DSYEV failed.'

    return
    end subroutine dsyev
!-------------------------------------------------------------------------------
!   TRED2.
!-------------------------------------------------------------------------------
    SUBROUTINE TRED2(NM,N,A,D,E,Z)
!---------------------------------------------------------------------------
!     TRIDIAGONALIZATION OF A SYMMETRIC MATRIX BY ORTHOGONAL TRANSFORMATIONS
!     (ALGORITHM OF HOUSEHOLDER)
!     CALLING MODE:
!               CALL TRED2(NM,N,A,D,E,Z)
!     INPUTS:
!     NM  (I4)  1ST DIMENSION OF MATRICES A AND Z IN CALLING PROGRAM
!     N   (I4)  SIZE OF A
!     A  (R*8)  TABLE(NM,N) STORING THE COEFFICIENTS OF SYMMETRIC A MATRIX
!               (LOWER HALF), A IS NOT DESTROYED DURING THE PROCESS
!               IF Z MATRIX HAS NOT THE SAME ADDRESS.
!     OUTPUTS:
!     D  (R*8)  MAIN DIAGONAL (N) OF REDUCED TRIDIAGONAL MATRIX
!     E  (R*8)  SUB-DIAGONAL (N) OF REDUCED TRIDIAGONAL MATRIX
!     Z  (R*8)  TABLE (NM,N) STORING THE ELEMENTS OF THE ORTHOGONAL 
!               TRANSFORMATION MATRIX.
!     REFERENCE:
!     J.H.WILKINSON,-C.REINSCH,R.S.MARTIN
!     HANDBOOK FOR AUTOMATIC COMPUTATION, VOL.2, LINEAR ALGEBRA
!     SPRINGER-VERLAG 1971.
!-----------------------------------------------------------------------
      use var_kind_def
      implicit none
      INTEGER(dpI):: I,J,K,L,N,NM
      REAL   (dpR):: A(NM,N),D(N),E(N),Z(NM,N),F,G,H,HH,SCALE

!     LOWER HALF OF A PUT INTO Z

      DO 10 I = 1,N
      DO 10 J = 1,I
   10 Z(I,J) = A(I,J)
      IF (N.EQ.1) GO TO 32

!     N-2 STAGE OF TRANSFORMATION

      DO 30 I = N,2,-1
      L = I-1
      H = 0.

!     CONDITIONNING BY NORM OF A

      SCALE = 0.
      IF (L.LT.2) GO TO 14
      DO 12 K = 1,L
   12 SCALE = SCALE+ABS(Z(I,K))
      IF (SCALE.NE.0.) GO TO 16

   14 E(I) = Z(I,L)
      GO TO 28

   16 DO 18 K = 1,L
      Z(I,K) = Z(I,K)/SCALE
      H = H+Z(I,K)*Z(I,K)
   18 CONTINUE

      F = Z(I,L)
      G = -SIGN(SQRT(H),F)
      E(I) = SCALE*G
      H = H-F*G
      Z(I,L) = F-G
      F = 0.
      DO 24 J = 1,L
      Z(J,I) = Z(I,J)/H
      G = 0.

!     ELEMENT OF A*U
      DO 20 K = 1,J
   20 G = G+Z(J,K)*Z(I,K)
      IF (L.GE.J+1) THEN
      DO 22 K = J+1,L
   22 G = G+Z(K,J)*Z(I,K)

!     ELEMENT OF P = A*U/H

      END IF
      E(J) = G/H
      F = F+E(J)*Z(I,J)
   24 CONTINUE

!     ELEMENT OF K

      HH = F/(H+H)

!     REDUCED FORM OF A

      DO 26 J = 1,L
      F = Z(I,J)
      G = E(J)-HH*F
      E(J) = G
      DO 26 K = 1,J
      Z(J,K) = Z(J,K)-F*E(K)-G*Z(I,K)
   26 CONTINUE
!
   28 D(I) = H
   30 CONTINUE

!     END OF TRANSFORMATION

   32 D(1) = 0.
      E(1) = 0.

!     ACCUMULATE TRANSFORMATION MATRICES IN Z

      DO 40 I = 1,N
      L = I-1
      IF (D(I).NE.0.) THEN
      DO 36 J = 1,L
      G = 0.
      DO 34 K = 1,L
   34 G = G+Z(I,K)*Z(K,J)
      DO 36 K = 1,L
      Z(K,J) = Z(K,J)-G*Z(K,I)
   36 CONTINUE
      END IF
      D(I) = Z(I,I)
      Z(I,I) = 1.
      IF (L.LT.1) GO TO 40
      DO 38 J = 1,L
      Z(I,J) = 0.
      Z(J,I) = 0.
   38 CONTINUE
   40 CONTINUE

      RETURN
      END
!-------------------------------------------------------------------------------
!   TQL2.
!-------------------------------------------------------------------------------
    SUBROUTINE TQL2(NM,N,D,E,Z,IER)
      use var_kind_def
      implicit none
!-------------------------------------------------------------------------
!     QL METHOD TO DETERMINE THE EIGENVALUES AND EIGENVECTORS OF:
!
!       1)  A SYMMETRIC TRIDIAGONAL MATRIX.
!       2)  A FULL SYMMETRIC MATRIX AFTER A PREVIOUS CALL TO TRED2.
!
!     CALLING MODE:
!               CALL TQL2(NM,N,D,E,Z,IER)
!     INPUTSS:
!     NM  (I4)  1ST DIMENSION OF MATRICES A AND Z IN CALLING PROGRAM
!     N   (I4)  SIZE OF Z
!     D  (R*8)  MAIN DIAGONAL (N) OF THE TRIDIAGONAL MATRIX
!     E  (R*8)  SUB-DIAGONAL (N) OF THE TRIDIAGONAL MATRIX
!     Z  (R*8)  TABLE (NM,N) STORING THE UNITY MATRIX IF THE TRIDIAGONAL
!               MATRIX IS DEFINED BY D AND E, CASE #1.
!               FOR CASE #2, IT CONTAINS THE ELEMENTS OF THE TRANSFORMATION
!               MATRIX AFTER A CALL TO TRED2.
!     OUTPUTS:
!     D  (R*8)  EIGENVALUES
!     Z  (R*8)  EIGENVECTORS
!     IER (I4)  ERROR CODE = 0,  CONVERGENCE OK.
!                          = L,  NO CONVERGENCE FOR THE Lth EIGENVALUE
!
!     REFERENCE:
!     J.H.WILKINSON,-C.REINSCH,R.S.MARTIN
!     HANDBOOK FOR AUTOMATIC COMPUTATION, VOL.2, LINEAR ALGEBRA
!     SPRINGER-VERLAG 1971.
!-------------------------------------------------------------------------
      INTEGER(dpI):: I,J,K,L,M,N,NM,JM,ier
      REAL   (dpR):: D(N),E(N),Z(NM,N),B,C,F,G,H,P,R,S,EPS,EPS1
      DATA EPS /0.D0/,JM /30/
      IER = 0
      IF (N.EQ.1) GO TO 38
!
!     MACHINE EPSILON
!
      IF (EPS.NE.0.D0) GO TO 12
      EPS = 1.D0
   10 EPS = EPS/2.D0
      EPS1 = 1.D0+EPS
      IF (EPS1.GT.1.D0) GO TO 10
!
   12 DO 14 I = 2,N
   14 E(I-1) = E(I)
      E(N) = 0.D0
      F = 0.D0
      B = 0.D0
!
      DO 28 L = 1,N
      J = 0
      H = EPS*(ABS(D(L))+ABS(E(L)))
      IF (B.LT.H) B = H
!
!     SEEK SMALLEST ELEMENT OF SUBDIAGONAL
!
      DO 16 M = L,N
      IF (ABS(E(M)).LE.B) GO TO 18
   16 CONTINUE
   18 IF (M.EQ.L) GO TO 26

!     START ITERATION

   20 IF (J.EQ.JM) GO TO 36
      J = J+1

!     SHIFT

      G = D(L)
      P = (D(L+1)-G)/(2.D0*E(L))
      R = SQRT(P*P+1.D0)
      D(L) = E(L)/(P+SIGN(R,P))
      H = G-D(L)
      DO 22 I = L+1,N
   22 D(I) = D(I)-H
      F = F+H

!     QL TRANSFORMATION

      P = D(M)
      C = 1.D0
      S = 0.D0
      DO 24 I = M-1,L,-1
      G = C*E(I)
      H = C*P
      IF (ABS(P).GE.ABS(E(I))) THEN
      C = E(I)/P
      R = SQRT(C*C+1.D0)
      E(I+1) = S*P*R
      S = C/R
      C = 1.D0/R
      ELSE
      C = P/E(I)
      R = SQRT(C*C+1.D0)
      E(I+1) = S*E(I)*R
      S = 1.D0/R
      C = C*S
      ENDIF
      P = C*D(I)-S*G
      D(I+1) = H+S*(C*G+S*D(I))

!     ELEMENTS OF EIGENVECTORS

      DO 24 K = 1,N
      H = Z(K,I+1)
      Z(K,I+1) = S*Z(K,I)+C*H
      Z(K,I) = Z(K,I)*C-S*H
   24 CONTINUE
      E(L) = S*P
      D(L) = C*P
      IF (ABS(E(L)).GT.B) GO TO 20

!     CONVERGENCE

   26 D(L) = D(L)+F
   28 CONTINUE

!     SORT EIGENVALUES AND EIGENVECTORS
!     IN ASVENDING ORDER

      DO 34 L = 2,N
      I = L-1
      K = I
      P = D(I)
      DO 30 J = L,N
      IF (D(J).GE.P) GO TO 30
      K = J
      P = D(J)
   30 CONTINUE
      IF (K.EQ.I) GO TO 34
      D(K) = D(I)
      D(I) = P
      DO 32 J = 1,N
      P = Z(J,I)
      Z(J,I) = Z(J,K)
   32 Z(J,K) = P
   34 CONTINUE
      GO TO 38

!     NO CONVERGENCE

   36 IER = L
   38 RETURN
      END
!-------------------------------------------------------------------------------
!   Symmetric Gauss-Seidel iteration.
!-------------------------------------------------------------------------------
    subroutine gs_iteration(N,iA,jA,M,max_iter,A,D,RHS,x)
    use var_kind_def
    implicit none
    integer(dpI),intent(in):: N,iA(*),jA(*),M,max_iter
    real   (dpR),intent(in):: A(M,*),D(M,*),RHS(*)
    integer(dpI):: i,j,k,iter,L,S
    real   (dpR):: x(*),v(100)

    if(M .eq. 5) then
        do iter=1,max_iter
            do i=1,N
                v(1:5)  =  RHS(5*i-4:5*i)
                do k=iA(i),iA(i+1)-1
                    j   =  jA(k)
                    if(i .eq. j)    cycle
                    L   =  5*k-4
                    S   =  5*j-4
                    v(1:5)  =  v(1:5)-A(1:5,L)*x(S)-A(1:5,L+1)*x(S+1)-A(1:5,L+2)*x(S+2) &
                            & -A(1:5,L+3)*x(S+3)-A(1:5,L+4)*x(S+4)
                end do
                S   =  5*i-4
                x(S:S+4)=  D(1:5,S  )*v(1)+D(1:5,S+1)*v(2)+D(1:5,S+2)*v(3) &
                        & +D(1:5,S+3)*v(4)+D(1:5,S+4)*v(5)
            end do

            do i=N,1,-1
                v(1:5)  =  RHS(5*i-4:5*i)
                do k=iA(i+1)-1,iA(i),-1
                    j   =  jA(k)
                    if(i .eq. j)    cycle
                    L   =  5*k-4
                    S   =  5*j-4
                    v(1:5)  =  v(1:5)-A(1:5,L)*x(S)-A(1:5,L+1)*x(S+1)-A(1:5,L+2)*x(S+2) &
                            & -A(1:5,L+3)*x(S+3)-A(1:5,L+4)*x(S+4)
                end do
                S   =  5*i-4
                x(S:S+4)=  D(1:5,S  )*v(1)+D(1:5,S+1)*v(2)+D(1:5,S+2)*v(3) &
                        & +D(1:5,S+3)*v(4)+D(1:5,S+4)*v(5)
            end do
        end do
        return
    elseif(M .eq. 1) then
        do iter=1,max_iter
            do i=1,N
                v(1)=  RHS(i)
                do k=iA(i),iA(i+1)-1
                    j   =  jA(k)
                    if(i .eq. j)    cycle
                    v(1)=  v(1)-A(1,k)*x(j)
                end do
                x(i)=  D(1,i)*v(1)
            end do

            do i=N,1,-1
                v(1)=  RHS(i)
                do k=iA(i+1)-1,iA(i),-1
                    j   =  jA(k)
                    if(i .eq. j)    cycle
                    L   =  5*k-4
                    S   =  5*j-4
                    v(1)=  v(1)-A(1,k)*x(j)
                end do
                x(i)=  D(1,i)*v(1)
            end do
        end do
        return
    end if

    do iter=1,max_iter
        do i=1,N
            v(1:M)  =  RHS(1+M*(i-1):M*i)
            do k=iA(i),iA(i+1)-1
                j   =  jA(k)
                if(i .eq. j)    cycle
                call DGEMV('N', M, M,-1.0d0, A(1,1+M*(k-1)), M, x(1+M*(j-1)), &
                    &  1, 1.0d0, v, 1)
            end do
            call DGEMV('N', M, M, 1.0d0, D(1,1+M*(i-1)), M, v, 1, 0.0d0, x(1+M*(i-1)), 1)
        end do

        do i=N,1,-1
            v(1:M)  =  RHS(1+M*(i-1):M*i)
            do k=iA(i+1)-1,iA(i),-1
                j   =  jA(k)
                if(i .eq. j)    cycle
                call DGEMV('N', M, M,-1.0d0, A(1,1+M*(k-1)), M, x(1+M*(j-1)), &
                    &  1, 1.0d0, v, 1)
            end do
            call DGEMV('N', M, M, 1.0d0, D(1,1+M*(i-1)), M, v, 1, 0.0d0, x(1+M*(i-1)), 1)
        end do
    end do

    return
    end subroutine gs_iteration
!-------------------------------------------------------------------------------
!   compute the LU decomposition for small size sparse matrix.
!-------------------------------------------------------------------------------
    subroutine dgetrf(N,M,A)
    use var_kind_def
    implicit none
    integer(dpI),intent(in):: N,M
    real   (dpR),intent(in):: A(M*M,*)
    integer(dpI):: i,j,k,ii,ik,ij,kj,kk,S
    real   (dpR):: X(M*M)

    S   =  M*M
    do i=1,N
        do k=1,i-1
            ik  =  k+(i-1)*N
            kk  =  k+(k-1)*N
            call DGEMM('N','N', M, M, M, 1.0d0, A(1,ik), M, A(1,kk), M, 0.0d0, X, M)
            call DCOPY(S, X, 1, A(1,ik), 1)

            do j=k+1,N
                ij  =  j+(i-1)*N
                kj  =  j+(k-1)*N
                call DGEMM('N','N', M, M, M,-1.0d0, X, M, A(1,kj), M, 1.0d0, A(1,ij), M)
            end do
        end do
        ii  =  i+(i-1)*N
        call mat_inv(M, A(1,ii))
    end do

    return
    end subroutine dgetrf
!-------------------------------------------------------------------------------
!   solution use the LU decomposition for small size sparse matrix.
!-------------------------------------------------------------------------------
    subroutine dgetrs(N,M,A,x)
    use var_kind_def
    implicit none
    integer(dpI),intent(in):: N,M
    real   (dpR),intent(in):: A(M*M,*)
    integer(dpI):: i,j,k,ii
    real   (dpR):: x(M,*),R(M)

!   ----------------------------------------------------------------------------
!   forward iteration.
    do i=2,N
        R(1:M)  =  x(1:M,i)
        do j=1,i-1
            k   =  j+N*(i-1)
            call dgemv('N', M, M,-1.0d0, A(1,k), M, x(1,j), 1, 1.0d0, R, 1)
        end do
        x(1:M,i)    =  R(1:M)
    end do
!   forward iteration.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   backward iteration.
    do i=N,1,-1
        R(1:M)  =  x(1:M,i)
        do j=N,i+1,-1
            k   =  j+N*(i-1)
            call dgemv('N', M, M,-1.0d0, A(1,k), M, x(1,j), 1, 1.0d0, R, 1)
        end do
        ii  =  i+N*(i-1)
        call dgemv('N', M, M, 1.0d0, A(1,ii), M, R, 1, 0.0d0, x(1,i), 1)
    end do
!   backward iteration.
!   ----------------------------------------------------------------------------

    return
    end subroutine dgetrs
!-------------------------------------------------------------------------------
!   BILU(0) decomposition, storage of LU is re-used.
!-------------------------------------------------------------------------------
    subroutine ilu0_decomposition(N,iA,jA,M,LU)
    use var_kind_def
    implicit none
    integer(dpI),intent(in):: N,iA(*),jA(*),M
    integer(dpI):: i,j,k,ii,ik,ij,kj,kk
    real   (dpR):: LU(M*M,*),X(M*M)

    if(M .eq. 5) then
    do i=1,N
    do ik=iA(i),iA(i+1)-1
        k   =  jA(ik)
        if(k .eq. i)    ii  =  ik
        if(k .ge. i)    cycle

        call is_exist(k, k, kk)
        x(1 :5 )=  LU(1 :5 ,ik)*LU(1 ,kk)+LU(6 :10,ik)*LU(2 ,kk)+LU(11:15,ik)*LU(3 ,kk) &
                & +LU(16:20,ik)*LU(4 ,kk)+LU(21:25,ik)*LU(5 ,kk)
        x(6 :10)=  LU(1 :5 ,ik)*LU(6 ,kk)+LU(6 :10,ik)*LU(7 ,kk)+LU(11:15,ik)*LU(8 ,kk) &
                & +LU(16:20,ik)*LU(9 ,kk)+LU(21:25,ik)*LU(10,kk)
        x(11:15)=  LU(1 :5 ,ik)*LU(11,kk)+LU(6 :10,ik)*LU(12,kk)+LU(11:15,ik)*LU(13,kk) &
                & +LU(16:20,ik)*LU(14,kk)+LU(21:25,ik)*LU(15,kk)
        x(16:20)=  LU(1 :5 ,ik)*LU(16,kk)+LU(6 :10,ik)*LU(17,kk)+LU(11:15,ik)*LU(18,kk) &
                & +LU(16:20,ik)*LU(19,kk)+LU(21:25,ik)*LU(20,kk)
        x(21:25)=  LU(1 :5 ,ik)*LU(21,kk)+LU(6 :10,ik)*LU(22,kk)+LU(11:15,ik)*LU(23,kk) &
                & +LU(16:20,ik)*LU(24,kk)+LU(21:25,ik)*LU(25,kk)
        call DCOPY(25, X, 1, LU(1,ik), 1)

        do ij=iA(i),iA(i+1)-1
            j   =  jA(ij)
            if(j .le. k)    cycle

            call is_exist(k, j, kj)
            if(kj .le. 0)   cycle
            LU(1 :5 ,ij)=  LU(1 :5 ,ij)-x(1:5)*LU(1 ,kj)-x(6 :10)*LU(2 ,kj) &
                        & -x(11:15)*LU(3 ,kj)-x(16:20)*LU(4 ,kj)-x(21:25)*LU(5 ,kj)
            LU(6 :10,ij)=  LU(6 :10,ij)-x(1:5)*LU(6 ,kj)-x(6 :10)*LU(7 ,kj) &
                        & -x(11:15)*LU(8 ,kj)-x(16:20)*LU(9 ,kj)-x(21:25)*LU(10,kj)
            LU(11:15,ij)=  LU(11:15,ij)-x(1:5)*LU(11,kj)-x(6 :10)*LU(12,kj) &
                        & -x(11:15)*LU(13,kj)-x(16:20)*LU(14,kj)-x(21:25)*LU(15,kj)
            LU(16:20,ij)=  LU(16:20,ij)-x(1:5)*LU(16,kj)-x(6 :10)*LU(17,kj) &
                        & -x(11:15)*LU(18,kj)-x(16:20)*LU(19,kj)-x(21:25)*LU(20,kj)
            LU(21:25,ij)=  LU(21:25,ij)-x(1:5)*LU(21,kj)-x(6 :10)*LU(22,kj) &
                        & -x(11:15)*LU(23,kj)-x(16:20)*LU(24,kj)-x(21:25)*LU(25,kj)
        end do
    end do
    call mat_inv(5, LU(1,ii))
    end do
    return
    end if

    do i=1,N
        do ik=iA(i),iA(i+1)-1
            k   =  jA(ik)
            if(k .eq. i)    ii  =  ik
            if(k .ge. i)    cycle

            call is_exist(k, k, kk)
            call DGEMM('N','N', M, M, M, 1.0d0, LU(1,ik), M, LU(1,kk), M, 0.0d0, X, M)
            call DCOPY(M*M, X, 1, LU(1,ik), 1)

            do ij=iA(i),iA(i+1)-1
                j   =  jA(ij)
                if(j .le. k)    cycle

                call is_exist(k, j, kj)
                if(kj .le. 0)   cycle
                call DGEMM('N','N', M, M, M,-1.0d0, X, M, LU(1,kj), M, 1.0d0, LU(1,ij), M)
            end do
        end do
        call mat_inv(M, LU(1,ii))
    end do

    return
    contains
!       ------------------------------------------------------------------------
!       check the existence of element(i,j).
!       ------------------------------------------------------------------------
        subroutine is_exist(i,j,ij)
        implicit none
        integer(dpI),intent(in):: i,j
        integer(dpI):: k,ij

        ij  = -1
        do k=iA(i),iA(i+1)-1
            if(jA(k) .eq. j) then
                ij  =  k
                return
            end if
        end do

        return
        end subroutine is_exist
    end subroutine ilu0_decomposition
!-------------------------------------------------------------------------------
!   BILU(0) solution.
!-------------------------------------------------------------------------------
    subroutine ilu0_solution(N,iA,jA,M,A,RHS,x)
    use var_kind_def
    implicit none
    integer(dpI),intent(in):: N,iA(*),jA(*),M
    real   (dpR),intent(in):: A(M*M,*),RHS(M,*)
    integer(dpI):: i,j,k,ii
    real   (dpR):: x(M,*),R(M)

    if(M .eq. 5) then
!   ----------------------------------------------------------------------------
!   forward iteration.
    x(1:M,1)    =  RHS(1:M,1)
    do i=2,N
        R(1:M)  =  RHS(1:M,i)
        do k=iA(i),iA(i+1)-1
            j   =  jA(k)
            if(j .ge. i)    cycle
!           call dgemv('N', M, M,-1.0d0, A(1,k), M, x(1,j), 1, 1.0d0, R, 1)
            R(1:5)  =  R(1:5)-A(1:5,k)*x(1,j)-A(6:10,k)*x(2,j)-A(11:15,k)*x(3,j) &
                    & -A(16:20,k)*x(4,j)-A(21:25,k)*x(5,j)
        end do
        x(1:M,i)=  R(1:M)
    end do
!   forward iteration.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   backward iteration.
    do i=N,1,-1
        R(1:M)  =  x(1:M,i)
        do k=iA(i+1)-1,iA(i),-1
            j   =  jA(k)
            if(i .eq. j)    ii  =  k
            if(i .ge. j)    cycle
!           call dgemv('N', M, M,-1.0d0, A(1,k), M, x(1,j), 1, 1.0d0, R, 1)
            R(1:5)  =  R(1:5)-A(1:5,k)*x(1,j)-A(6:10,k)*x(2,j)-A(11:15,k)*x(3,j) &
                    & -A(16:20,k)*x(4,j)-A(21:25,k)*x(5,j)
        end do
!       call dgemv('N', M, M, 1.0d0, A(1,ii), M, R, 1, 0.0d0, x(1,i), 1)
        x(1:5,i)=  A(1 :5 ,ii)*R(1)+A(6 :10,ii)*R(2)+A(11:15,ii)*R(3) &
                & +A(16:20,ii)*R(4)+A(21:25,ii)*R(5)
    end do
!   backward iteration.
!   ----------------------------------------------------------------------------
    return
    end if

!   ----------------------------------------------------------------------------
!   forward iteration.
    x(1:M,1)    =  RHS(1:M,1)
    do i=2,N
        R(1:M)  =  RHS(1:M,i)
        do k=iA(i),iA(i+1)-1
            j   =  jA(k)
            if(j .ge. i)    cycle
            call dgemv('N', M, M,-1.0d0, A(1,k), M, x(1,j), 1, 1.0d0, R, 1)
        end do
        x(1:M,i)=  R(1:M)
    end do
!   forward iteration.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   backward iteration.
    do i=N,1,-1
        R(1:M)  =  x(1:M,i)
        do k=iA(i+1)-1,iA(i),-1
            j   =  jA(k)
            if(i .eq. j)    ii  =  k
            if(i .ge. j)    cycle
            call dgemv('N', M, M,-1.0d0, A(1,k), M, x(1,j), 1, 1.0d0, R, 1)
        end do
        call dgemv('N', M, M, 1.0d0, A(1,ii), M, R, 1, 0.0d0, x(1,i), 1)
    end do
!   backward iteration.
!   ----------------------------------------------------------------------------

    return
    end subroutine ilu0_solution
!-------------------------------------------------------------------------------
!   LU decompostion of block-tridiagonal system.
!-------------------------------------------------------------------------------
    subroutine B_dgttrf(N,L,D,U)
    use var_kind_def
    implicit none
    integer(dpI),parameter:: M=5
    integer(dpI):: N,i,i1,i0
    real   (dpR):: L(M,*),D(M,*),U(M,*),tm1(M,M)

!   trid(L,D,U)=[LL,I]*[DD,U]
!   Note that D(i,1) is replaced by inv(DD(i,i  )) to save memory.
!   Note that L(i,1) is replaced by     LL(i,i-1)  to save memory.
!   Note that U(i,1) equals to UU(i,i+1).

    call mat_inv(M, D)
    do i=2,N
        i1  = (i-1)*M+1
        i0  =  i   *M
        tm1(:,1:M)  =  L(:,i1:i0)
        call DGEMM('N', 'N', M, M, M, 1.0d0, tm1, M, D(1,i1-M), M, 0.0d0, L(1,i1), M)

        call DGEMM('N', 'N', M, M, M,-1.0d0, L(1,i1), M, U(1,i1-M), M, 1.0d0, D(1,i1), M)
        call mat_inv(M, D(1,i1))
    end do

    return
    end subroutine B_dgttrf
!-------------------------------------------------------------------------------
!   Block tridiagonal solver using LU decomposition.
!-------------------------------------------------------------------------------
    subroutine B_dgttrs(N,L,D,U,Rhs)
    use var_kind_def
    implicit none
    integer(dpI),parameter:: M=5
    integer(dpI):: N,i,i1
    real   (dpR):: L(M,*),D(M,*),U(M,*),Rhs(M,*),tmv1(M)

    if(N .lt. 2)    stop 'Error: wrong input for B_dgttrs.'
    do i=2,N
        i1  = (i-1)*M+1
        call dgemv('N', M, M,-1.0d0, L(1,i1), M, Rhs(:,i-1), 1, 1.0d0, Rhs(:,i), 1)
    end do

    tmv1(:) =  Rhs(:,N)
    call dgemv('N', M, M, 1.0d0, D(1,i1), M, tmv1, 1, 0.0d0, Rhs(:,N), 1)
    do i=N-1,1,-1
        i1  = (i-1)*M+1

        call dgemv('N', M, M,-1.0d0, U(1,i1), M, Rhs(1,i+1), 1, 1.0d0, Rhs(1,i), 1)
        tmv1(:) =  Rhs(:,i)
        call dgemv('N', M, M, 1.0d0, D(1,i1), M, tmv1, 1, 0.0d0, Rhs(1,i), 1)
    end do

    return
    end subroutine B_dgttrs
!-------------------------------------------------------------------------------
!   matrix re-ordering using the HSL-MC60.
!-------------------------------------------------------------------------------
    subroutine matrix_ordering(N,iA,jA,perm)
    use var_kind_def
    implicit none
    integer(dpI),intent(in):: N,iA(*),jA(*)
    logical(dpL):: is_supervariables
    integer(dpI):: perm(*),info(4),options(2),n_sup,i
    real   (dpR):: weight(2),rinfo(4)
    integer(dpI),allocatable:: iw(:),vars(:),permsv(:),pair(:,:),svar(:)
    real   (dpR),allocatable:: w(:)

    is_supervariables   =  .false.
    options             = (/1, 0/)
    weight              = (/2.0d0, 1.0d0/)

    allocate(vars  (N))
    allocate(permsv(N))
    allocate(svar  (N))
    allocate(pair(2,N/2))
    allocate(iw(3*N+1))
    allocate(w(N))
    permsv  =  1

    if(is_supervariables) then
        call mc60bd(N, iA(N+1)-1, jA, iA, n_sup, svar, vars, iw)
        call mc60cd(N, n_sup, iA(N+1)-1, jA, iA, vars, options, permsv, weight, &
            &  pair, info, iw, w)
        call mc60fd(N, n_sup, iA(N+1)-1, jA, iA, vars, permsv, iw, rinfo)
        call mc60dd(N, n_sup, svar, vars, permsv, perm, iw)
    else
        n_sup   =  N
        vars    =  1
        call mc60cd(N, n_sup, iA(N+1)-1, jA, iA, vars, options, perm, weight, &
            &  pair, info, iw, w)
        call mc60fd(N, n_sup, iA(N+1)-1, jA, iA, vars, perm, iw, rinfo)
    end if
    do i=1,N
        if((perm(i) .le. 0) .or. (perm(i) .gt. N))  stop 'Error: HSL fails to reorder.'
    end do

    if(allocated(iw    ))   deallocate(iw    )
    if(allocated(vars  ))   deallocate(vars  )
    if(allocated(permsv))   deallocate(permsv)
    if(allocated(pair  ))   deallocate(pair  )
    if(allocated(svar  ))   deallocate(svar  )
    if(allocated(w     ))   deallocate(w     )

    return
    end subroutine matrix_ordering
!-------------------------------------------------------------------------------
!   reorder the matrix according to permutation.
!-------------------------------------------------------------------------------
    subroutine matrix_reorder(nvtx,bsize,perm,iA,jA,A)
    use var_kind_def
    implicit none
    integer(dpI),intent(in):: nvtx,bsize,perm(*)
    logical(dpL):: ltmp
    integer(dpI):: iA(*),jA(*),i,j,k,i_new,j_new,M,j_min,j_max,band_1,band_2
    real   (dpR):: A(*),band_all_1,band_all_2
    integer(dpI),allocatable:: iA_new(:),jA_new(:)
    real   (dpR),allocatable:: A_new(:)

!   ----------------------------------------------------------------------------
!   get the original bandwidth.
    band_1      =  0
    band_all_1  =  0.0d0
    do i=1,nvtx
        j_min   =  huge(1)
        j_max   = -1
        do j=iA(i),iA(i+1)-1
            j_min   =  min(j_min, jA(j))
            j_max   =  max(j_max, jA(j))
        end do
        band_1      =  max(band_1, j_max-j_min+1)
        band_all_1  =  band_all_1+real(j_max-j_min+1, dpR)
    end do
!   get the original bandwidth.
!   ----------------------------------------------------------------------------

    M   =  bsize*bsize
    allocate(iA_new(nvtx+1))
    allocate(jA_new(iA(nvtx+1)-1))
    allocate( A_new(M*(iA(nvtx+1)-1)))
    do i=1,nvtx
        iA_new(perm(i)+1)   =  iA(i+1)-iA(i)
    end do
    iA_new(1)   =  1
    do i=1,nvtx
        iA_new(i+1) =  iA_new(i)+iA_new(i+1)
    end do
    do i=1,nvtx
        i_new   =  perm(i)
        do j=iA(i),iA(i+1)-1
            jA_new(iA_new(i_new)+j-iA(i))   =  perm(jA(j))
        end do
    end do
    do i=1,nvtx
        i_new   =  perm(i)
        call iqsort(.true., iA_new(i_new), iA_new(i_new+1)-1, jA_new)
        do j=iA(i),iA(i+1)-1
            j_new   =  perm(jA(j))
            ltmp    =  .false.
            do k=iA_new(i_new),iA_new(i_new+1)-1
                ltmp=  jA_new(k) .eq. j_new
                if(ltmp)    exit
            end do
            if(.not. ltmp)  stop 'Error: failed to reorder the LHS.'
            call DCOPY(M, A(1+M*(j-1)), 1, A_new(1+M*(k-1)), 1)
        end do
    end do

    iA(1:nvtx+1)=  iA_new(1:nvtx+1)
    jA(1:iA(nvtx+1)-1)  =  jA_new(1:iA(nvtx+1)-1)
    call DCOPY(M*(iA(nvtx+1)-1), A_new, 1, A, 1)

!   ----------------------------------------------------------------------------
!   get the new bandwidth.
    band_2      =  0
    band_all_2  =  0.0d0
    do i=1,nvtx
        j_min   =  huge(1)
        j_max   = -1
        do j=iA(i),iA(i+1)-1
            j_min   =  min(j_min, jA(j))
            j_max   =  max(j_max, jA(j))
        end do
        band_2      =  max(band_2, j_max-j_min+1)
        band_all_2  =  band_all_2+real(j_max-j_min+1, dpR)
    end do
!   if(band_all_2 .gt. band_all_1)  stop 'Error: fails to improve the bandwidth.'
!   get the new bandwidth.
!   ----------------------------------------------------------------------------

    if(allocated(iA_new))   deallocate(iA_new)
    if(allocated(jA_new))   deallocate(jA_new)
    if(allocated( A_new))   deallocate( A_new)

    return
    end subroutine matrix_reorder

!-------------------------------------------------------------------------------
!   reorder the symmetric matrix with LDU format according to permutation.
!-------------------------------------------------------------------------------
    subroutine matrix_reorder_ldu(uptTopo,bsize,nvtx,perm,iA,jA,A)
    use var_kind_def
    use var_parallel, only: myid
    implicit none
    integer(dpI),intent(in):: nvtx,perm(*),bsize
    logical(dpL),intent(in):: uptTopo
    logical(dpL):: ltmp
    integer(dpI):: iA(*),jA(*),i,j,k,i_new,j_new,j_min,j_max,band_1,band_2,i_ldu,l
    real   (dpR):: A(*),band_all_1,band_all_2
    integer(dpI),allocatable:: iA_new(:),jA_new(:),iA_ldu(:),jA_ldu(:)
    real   (dpR),allocatable:: A_new(:)

!   ----------------------------------------------------------------------------
!   get the original bandwidth.
    band_1      =  0
    band_all_1  =  0.0d0
    do i=1,nvtx
        j_min   =  huge(1)
        j_max   = -1
        do j=iA(i),iA(i+1)-1
            j_min   =  min(j_min, jA(j))
            j_max   =  max(j_max, jA(j))
        end do
        band_1      =  max(band_1, j_max-j_min+1)
        if(j_max .ne. -1) band_all_1  =  band_all_1+real(j_max-j_min+1, dpR)
    end do
    ! write(*,*),"bandwidth", band_all_1
!   get the original bandwidth.
!   ----------------------------------------------------------------------------

    allocate(iA_new(nvtx+1))
    allocate(jA_new(iA(nvtx+1)-1))
    allocate( A_new(bsize*(iA(nvtx+1)-1)))
    A_new = -1
    do i=1,nvtx
        iA_new(perm(i)+1)   =  iA(i+1)-iA(i)
    end do
    iA_new(1)   =  1
    do i=1,nvtx
        iA_new(i+1) =  iA_new(i)+iA_new(i+1)
    end do
    do i=1,nvtx
        i_new   =  perm(i)
        do j=iA(i),iA(i+1)-1
            jA_new(iA_new(i_new)+j-iA(i))   =  perm(jA(j))
        end do
    end do

    allocate(iA_ldu(nvtx+1))
    allocate(jA_ldu((iA(nvtx+1)-1)/2))
    iA_ldu = 0
    i_ldu = 1
    do i=1,nvtx
        do j=iA(i),iA(i+1)-1
            if(i .gt. jA(j)) cycle
            iA_ldu(i+1) = iA_ldu(i+1) + 1
            jA_ldu(i_ldu) = jA(j)
            i_ldu = i_ldu + 1
        end do
    end do
    ! write(*,*),i_ldu,(iA(nvtx+1)-1)/2

    iA_ldu(1) = 1
    do i=1,nvtx
        iA_ldu(i+1) = iA_ldu(i+1) + iA_ldu(i)
    end do

    ! allocate(edge_map_csr((iA(nvtx+1)-1)))
    ! allocate(edge_map_ldu((iA(nvtx+1)-1)/2))
    i_ldu = 0
    do i=1,nvtx
        i_new   =  perm(i)
        call iqsort(.true., iA_new(i_new), iA_new(i_new+1)-1, jA_new)
        do j=iA(i),iA(i+1)-1
            j_new   =  perm(jA(j))
            ltmp    =  .false.
            do k=iA_new(i_new),iA_new(i_new+1)-1
                ltmp=  jA_new(k) .eq. j_new
                if(ltmp)    exit
            end do
            if(.not. ltmp)  stop 'Error: failed to reorder the LHS.'
            if(i .gt. jA(j)) then
                ltmp    =  .false.
                do l=iA_ldu(jA(j)),iA_ldu(jA(j)+1)-1
                    ltmp = jA_ldu(l) .eq. i
                    if(ltmp) exit
                end do
                if(.not. ltmp)  stop 'Error: failed to find the lower triangle.'
                ! A_new(k) = A(l)
                call DCOPY(bsize, A(bsize*(l-1)+1), 1, A_new(bsize*(k-1)+1), 1)
                ! edge_map_csr(j) = k
            else 
                i_ldu = i_ldu+1
                ! edge_map_ldu(j) = k
                ! if(k==1 .and. myid==0) write(*,*),i_ldu,A_new(k),A(i_ldu)
                call DCOPY(bsize, A(bsize*(i_ldu-1)+1), 1, A_new(bsize*(k-1)+1), 1)
            end if
            ! A_new(k) = A(i_ldu)
        end do
    end do

    if(uptTopo) then
        iA(1:nvtx+1)=  iA_new(1:nvtx+1)
        jA(1:iA(nvtx+1)-1)  =  jA_new(1:iA(nvtx+1)-1)
    end if 

    i_ldu = 0
    do i=1,nvtx
        do j=iA_new(i),iA_new(i+1)-1
            if(i .gt. jA_new(j)) cycle
            i_ldu = i_ldu+1
            call DCOPY(bsize, A_new(bsize*(j-1)+1), 1, A(bsize*(i_ldu-1)+1), 1)
            ! A(i_ldu) = A_new(j)
        end do
    end do
    ! write(*,*),i_ldu,(iA(nvtx+1)-1)/2
    ! call DCOPY(iA_new(nvtx+1)-1, A_new, 1, A, 1)

!   ----------------------------------------------------------------------------
!   get the new bandwidth.
    band_2      =  0
    band_all_2  =  0.0d0
    do i=1,nvtx
        j_min   =  huge(1)
        j_max   = -1
        do j=iA(i),iA(i+1)-1
            j_min   =  min(j_min, jA(j))
            j_max   =  max(j_max, jA(j))
        end do
        band_2      =  max(band_2, j_max-j_min+1)
        if(j_max .ne. -1) band_all_2  =  band_all_2+real(j_max-j_min+1, dpR)
    end do
    ! write(*,*),"bandwidth", band_all_2
  ! if(band_all_2 .gt. band_all_1)  stop 'Error: fails to improve the bandwidth.'
!   get the new bandwidth.
!   ----------------------------------------------------------------------------

    if(allocated(iA_new))   deallocate(iA_new)
    if(allocated(jA_new))   deallocate(jA_new)
    if(allocated( A_new))   deallocate( A_new)

    return
    end subroutine matrix_reorder_ldu

