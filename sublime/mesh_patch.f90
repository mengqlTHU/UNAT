!-------------------------------------------------------------------------------
!   module AVG.
!-------------------------------------------------------------------------------
    module var_avg
        use var_kind_def
        implicit none
        private
        public:: avg,n_avg
        public:: nprc_avg,myid_avg,mpi_comm_avg
        public:: avg_add,avg_add_patch,avg_synchronize
        public:: abf1,abf2

        integer(dpI),parameter:: max_patch  =  10

        type type_patch
            integer(dpI):: myid     =  0
            integer(dpI):: ID       =  0
            integer(dpI):: ele_type =  0
            integer(dpI):: n_ele    =  0
            integer(dpI):: n_vtx    =  0
            integer(dpI):: npe      =  0
            integer(dpI),allocatable:: n2e(:,:),ID_vtx(:)
            integer(dpI),allocatable:: mortar(:),ID_r(:)
            real   (dpR),allocatable:: xyz(:,:)
        end type type_patch

        type type_avg
            integer(dpI):: lev  =  0
            integer(dpI):: ID   =  0
            integer(dpI):: bct  =  0
            integer(dpI):: nfacL=  0
            integer(dpI):: nfacR=  0
            integer(dpI):: idmL =  0
            integer(dpI):: idmR =  0
            type(type_patch):: facL(max_patch),facR(max_patch)

            integer(dpI):: n_vtx_L  =  0
            real   (dpR):: r_reot       =  1.0d0
            real   (dpR):: p_reot       =  1.0d5
            integer(dpI),allocatable:: ID_vtx(:),ID_R(:)
            real   (dpR),allocatable:: xyz(:,:)
            real   (dpR),allocatable:: r_L(:),A_L(:)
            real   (dpR),allocatable:: u  (:,:)
        end type type_avg
        type(type_avg),save:: avg(10)
        integer(dpI):: n_avg    =  0
        integer(dpI):: mpi_group_avg,mpi_comm_avg
        integer(dpI):: nprc_avg = -1
        integer(dpI):: myid_avg = -1
        real   (dpR),allocatable:: abf1(:),abf2(:)

        contains
!       ------------------------------------------------------------------------
!       add a patch.
!       ------------------------------------------------------------------------
        subroutine patch_add(p,myid,ID,ele_type,n_ele,n2e,n_vtx,ID_vtx,xyz)
        use var_kind_def
        use var_cgns, only: npe_ele
        use var_global, only: err_mem
        implicit none
        integer(dpI),intent(in):: myid,ID,ele_type,n_ele,n2e(*),n_vtx,ID_vtx(*)
        real   (dpR),intent(in):: xyz(*)
        type(type_patch):: p

        p%myid      =  myid
        p%ID        =  ID
        p%ele_type  =  ele_type
        p%npe       =  npe_ele(ele_type)
        p%n_ele     =  n_ele
        p%n_vtx     =  n_vtx
        allocate(p%n2e   (p%npe, p%n_ele), stat=err_mem)
        allocate(p%ID_vtx(       p%n_vtx), stat=err_mem)
        allocate(p%xyz   (3    , p%n_vtx), stat=err_mem)
        call ICOPY(p%npe*p%n_ele, n2e   , 1, p%n2e   , 1)
        call ICOPY(      p%n_vtx, ID_vtx, 1, p%ID_vtx, 1)
        call DCOPY(3    *p%n_vtx, xyz   , 1, p%xyz   , 1)

        return
        end subroutine patch_add
!       ------------------------------------------------------------------------
!       sort the elements.
!       ------------------------------------------------------------------------
        subroutine patch_get_r(iavg,sec,n_sec)
        use var_kind_def
        use var_global, only: err_mem
        implicit none
        integer(dpI),intent(in):: iavg,n_sec
        logical(dpL):: ltmp
        integer(dpI):: isec,n_vtx,iele,ipe,L,R,i,n_R
        real   (dpR):: a,b,c,eps,rtmp
        type(type_patch):: sec(*)
        integer(dpI),allocatable:: v(:)
        real   (dpR),allocatable:: rz(:,:)

!       ------------------------------------------------------------------------
!       record the vertex used by the quilt.
        n_vtx   =  0
        do isec=1,n_sec
            n_vtx   =  n_vtx+sec(isec)%n_vtx
        end do
        if(n_vtx .gt. 0)    allocate(v(n_vtx), stat=err_mem)
        n_vtx   =  0
        do isec=1,n_sec
            call ICOPY(sec(isec)%n_vtx, sec(isec)%ID_vtx, 1, v(n_vtx+1), 1)
            n_vtx   =  n_vtx+sec(isec)%n_vtx
        end do
        call simplify_series(n_vtx, 1, 1, v)
        allocate(avg(iavg)%xyz(3 ,n_vtx), stat=err_mem)
        allocate(avg(iavg)%ID_vtx(n_vtx), stat=err_mem)
        allocate(avg(iavg)%ID_r  (n_vtx), stat=err_mem)
        call ICOPY(n_vtx, v, 1, avg(iavg)%ID_vtx, 1)
        avg(iavg)%n_vtx_L   =  n_vtx
!       record the vertex used by the quilt.
!       ------------------------------------------------------------------------

!       ------------------------------------------------------------------------
!       merge the patches into quilt.
        do isec=1,n_sec
        do iele=1,sec(isec)%n_ele
        do ipe=1,sec(isec)%npe
            call ib_search(1, n_vtx, 1, 1, v, sec(isec)%ID_vtx(sec(isec)%n2e(ipe,iele)), &
                &  ltmp, L, R)
            if((.not. ltmp) .or. (L .ne. R))    stop 'Error: patch_get_r fails.'
            avg(iavg)%xyz(1:3,L)    =  sec(isec)%xyz(1:3,sec(isec)%n2e(ipe,iele))
            sec(isec)%n2e(ipe,iele) =  L
        end do
        end do
        end do
!       merge the patches into quilt.
!       ------------------------------------------------------------------------

!       ------------------------------------------------------------------------
!       try to find the R and azimuthal directions.
        allocate(rz(3,n_vtx), stat=err_mem)
        do i=1,n_vtx
            rz(1,i) =  real(i, dpR)
            rz(2,i) =  sqrt(avg(iavg)%xyz(1,i)**2+avg(iavg)%xyz(2,i)**2 &
                    &      +avg(iavg)%xyz(3,i)**2)
            rz(3,i) =  1.0d0
        end do
        call dqsortcols(.true., 1, n_vtx, 2, 3, rz)
        rtmp=  1.0d-2*(rz(2,n_vtx)-rz(2,1))
        eps =  1.0d-8*(rz(2,n_vtx)-rz(2,1))

        n_r =  1
        do i=2,n_vtx-2
            a   =  abs(rz(2,i  )-rz(2,i-1))
            b   =  abs(rz(2,i+1)-rz(2,i  ))
            c   =  abs(rz(2,i+2)-rz(2,i+1))
            if((b .gt. 1.0d1*(a+eps)) .and. (b .gt. 1.0d1*(c+eps))) then
                n_r             =  n_r+1
                rz(3,i+1:n_vtx) =  real(n_r, dpR)
            elseif(b .ge. rtmp) then
                n_r             =  n_r+1
                rz(3,i+1:n_vtx) =  real(n_r, dpR)
            end if
        end do
        call dqsortcols(.true., 1, n_vtx, 1, 3, rz)
!       try to find the R and azimuthal directions.
!       ------------------------------------------------------------------------

!       ------------------------------------------------------------------------
!       check.
        do isec=1,n_sec
        do iele=1,sec(isec)%n_ele
            do ipe=1,sec(isec)%npe
                v(ipe)  =  nint(rz(3,sec(isec)%n2e(ipe,iele)), dpI)
            end do
            L   =  minval(v(1:sec(isec)%npe))
            R   =  maxval(v(1:sec(isec)%npe))
            if(R-L .ne. 1)  stop 'Error: patch_get_r fails.'
        end do
        end do
        avg(iavg)%idmL  =  0
        do i=1,avg(iavg)%n_vtx_L
            avg(iavg)%ID_r(i)   =  nint(rz(3,i), dpI)
            avg(iavg)%idmL      =  max(avg(iavg)%idmL, avg(iavg)%ID_r(i))
        end do
!       check.
!       ------------------------------------------------------------------------

        if(allocated(v ))   deallocate(v )
        if(allocated(rz))   deallocate(rz)

        return
        end subroutine patch_get_r
!       ------------------------------------------------------------------------
!       add an avg interface.
!       ------------------------------------------------------------------------
        subroutine avg_add(lev,ID,bct,i_avg)
        use var_kind_def
        implicit none
        integer(dpI),intent(in):: lev,ID,bct
        integer(dpI):: i_avg

        do i_avg=1,n_avg
            if((avg(i_avg)%lev .eq. lev) .and. (avg(i_avg)%ID .eq. ID) .and. &
              &(avg(i_avg)%bct .eq. bct))   return
        end do
        n_avg           =  n_avg+1
        avg(n_avg)%lev  =  lev
        avg(n_avg)%ID   =  ID
        avg(n_avg)%bct  =  bct
        i_avg           =  n_avg

        return
        end subroutine avg_add
!       ------------------------------------------------------------------------
!       add a patch to the avg.
!       ------------------------------------------------------------------------
        subroutine avg_add_patch(a,is_L,prc,ID,ele_type,n_ele,n2e,n_vtx,ID_vtx,xyz,ifac)
        use var_kind_def
        implicit none
        logical(dpL),intent(in):: is_L
        integer(dpI),intent(in):: prc,ID,ele_type,n_ele,n2e(*),n_vtx,ID_vtx(*)
        real   (dpR),intent(in):: xyz(*)
        logical(dpL):: ltmp
        integer(dpI):: ifac
        type(type_avg):: a

        ltmp=  .true.
        if(is_L) then
            do ifac=1,a%nfacL
                if((a%facL(ifac)%ID .eq. ID) .and. (a%facL(ifac)%myid .eq. prc)) then
                    ltmp=  .false.
                    exit
                end if
            end do
            if(ltmp)    a%nfacL =  a%nfacL+1
            call patch_add(a%facL(a%nfacL),prc,ID,ele_type,n_ele,n2e,n_vtx,ID_vtx,xyz)
        else
            do ifac=1,a%nfacR
                if(a%facR(ifac)%ID .eq. ID) then
                    ltmp=  .false.
                    exit
                end if
            end do
            if(ltmp)    a%nfacR =  a%nfacR+1
            call patch_add(a%facR(a%nfacR),prc,ID,ele_type,n_ele,n2e,n_vtx,ID_vtx,xyz)
        end if

        return
        end subroutine avg_add_patch
!       ------------------------------------------------------------------------
!       set mpi environment for AVG.
!       ------------------------------------------------------------------------
        subroutine avg_set_mpi
        use var_kind_def
        use var_parallel
        implicit none
        integer(dpI):: i,gp_old(0:nprc-1),gp_new(0:nprc-1)

        call mpi_allgather(n_avg, 1, mpi_dpI, gp_old, 1, mpi_dpI, mpi_comm_world, mpi_err)
        nprc_avg=  0
        do i=0,nprc-1
            if(gp_old(i) .le. 0)    cycle
            nprc_avg            =  nprc_avg+1
            gp_new(nprc_avg-1)  =  i
        end do
        call mpi_group_incl(mpi_group_world, nprc_avg, gp_new, mpi_group_avg, mpi_err)
        call mpi_comm_create(mpi_comm_world, mpi_group_avg, mpi_comm_avg, mpi_err)
        if(n_avg .le. 0)    return
        call mpi_comm_size(mpi_comm_avg, nprc_avg, mpi_err)
        call mpi_comm_rank(mpi_comm_avg, myid_avg, mpi_err)

        return
        end subroutine avg_set_mpi
!       ------------------------------------------------------------------------
!       synchronize the AVG.
!       ------------------------------------------------------------------------
        subroutine avg_synchronize
        use var_kind_def
        use var_cgns, only: TRI_3
        use var_global, only: err_mem,z=>rotation_axis
        use var_parallel, only: mpi_dpI,mpi_err,mpi_sum,myid
        implicit none
        logical(dpL):: ltmp
        integer(dpI):: i,idx,isec,iavg,iele,ipe,lev,ID,bct,nfacL,nfacR,prc,ele_type, &
                    &  n_ele,n_vtx,npe,ifac
        real   (dpR):: c(3),n(3),xyz(3,4),rtmp
        integer(dpI),allocatable:: v1(:),v2(:)
        real   (dpR),allocatable:: S(:),R(:)

        call avg_set_mpi
        if(myid_avg .lt. 0) return

!       ------------------------------------------------------------------------
!       allocate memory.
        idx =  0
        do iavg=1,n_avg
!           lev, ID, bct, nfacL, nfacR
            idx =  idx+5

            do isec=1,avg(iavg)%nfacL
!               myid, ID, ele_type, n_ele, n_vtx, npe
                idx =  idx+6
!               n2e, ID_vtx, xyz
                idx =  idx+avg(iavg)%facL(isec)%npe*avg(iavg)%facL(isec)%n_ele &
                    & +4*avg(iavg)%facL(isec)%n_vtx
            end do

            do isec=1,avg(iavg)%nfacR
!               myid, ID, ele_type, n_ele, n_vtx, npe
                idx =  idx+6
!               n2e, ID_vtx, xyz
                idx =  idx+avg(iavg)%facL(isec)%npe*avg(iavg)%facL(isec)%n_ele &
                    & +4*avg(iavg)%facL(isec)%n_vtx
            end do
        end do
        call mpi_allreduce(idx, i, 1, mpi_dpI, mpi_sum, mpi_comm_avg, mpi_err)
        allocate(S (max(idx, 1)), stat=err_mem)
        allocate(R (max(i  , 1)), stat=err_mem)
        if(i .gt. 0)    allocate(v1(i), stat=err_mem)
        if(i .gt. 0)    allocate(v2(i), stat=err_mem)
!       allocate memory.
!       ------------------------------------------------------------------------

!       ------------------------------------------------------------------------
!       allgather and synchronize.
        idx =  0
        do iavg=1,n_avg
!           lev, ID, bct, nfacL, nfacR
            S(idx+1)=  real(avg(iavg)%lev  , dpR)
            S(idx+2)=  real(avg(iavg)%ID   , dpR)
            S(idx+3)=  real(avg(iavg)%bct  , dpR)
            S(idx+4)=  real(avg(iavg)%nfacL, dpR)
            S(idx+5)=  real(avg(iavg)%nfacR, dpR)
            idx     =  idx+5

            do isec=1,avg(iavg)%nfacL
!               myid, ID, ele_type, n_ele, n_vtx, npe
                S(idx+1)=  real(avg(iavg)%facL(isec)%myid    , dpR)
                S(idx+2)=  real(avg(iavg)%facL(isec)%ID      , dpR)
                S(idx+3)=  real(avg(iavg)%facL(isec)%ele_type, dpR)
                S(idx+4)=  real(avg(iavg)%facL(isec)%n_ele   , dpR)
                S(idx+5)=  real(avg(iavg)%facL(isec)%n_vtx   , dpR)
                S(idx+6)=  real(avg(iavg)%facL(isec)%npe     , dpR)
                idx     =  idx+6

!               n2e.
                do iele=1,avg(iavg)%facL(isec)%n_ele
                    do ipe=1,avg(iavg)%facL(isec)%npe
                        S(idx+ipe)  =  real(avg(iavg)%facL(isec)%n2e(ipe,iele), dpR)
                    end do
                    idx =  idx+avg(iavg)%facL(isec)%npe
                end do

!               ID_vtx.
                do i=1,avg(iavg)%facL(isec)%n_vtx
                    S(idx+i)=  real(avg(iavg)%facL(isec)%ID_vtx(i), dpR)
                end do
                idx =  idx+avg(iavg)%facL(isec)%n_vtx

!               xyz.
                do i=1,avg(iavg)%facL(isec)%n_vtx
                    S(idx+1)=  avg(iavg)%facL(isec)%xyz(1,i)
                    S(idx+2)=  avg(iavg)%facL(isec)%xyz(2,i)
                    S(idx+3)=  avg(iavg)%facL(isec)%xyz(3,i)
                    idx     =  idx+3
                end do
            end do
        end do
        call r_allgather(mpi_comm_avg, nprc_avg, idx, S, R)
        if(allocated(S))    deallocate(S)

        idx =  0
        do while(idx .lt. size(R))
!           lev, ID, bct, nfacL, nfacR
            lev     =  nint(R(idx+1), dpI)
            ID      =  nint(R(idx+2), dpI)
            bct     =  nint(R(idx+3), dpI)
            nfacL   =  nint(R(idx+4), dpI)
            nfacR   =  nint(R(idx+5), dpI)
            idx     =  idx+5

            ltmp=  .false.
            do iavg=1,n_avg
                ltmp= (lev .eq. avg(iavg)%lev) .and. (ID .eq. avg(iavg)%ID) .and. &
                    & (bct .eq. avg(iavg)%bct)
                if(ltmp)    exit
            end do

            do isec=1,nfacL
!               myid, ID, ele_type, n_ele, n_vtx, npe
                prc     =  nint(R(idx+1), dpI)
                ID      =  nint(R(idx+2), dpI)
                ele_type=  nint(R(idx+3), dpI)
                n_ele   =  nint(R(idx+4), dpI)
                n_vtx   =  nint(R(idx+5), dpI)
                npe     =  nint(R(idx+6), dpI)
                idx     =  idx+6

                if((.not. ltmp) .or. (prc .eq. myid)) then
                    idx =  idx+npe*n_ele+4*n_vtx
                else
!                   n2e.
                    do i=1,npe*n_ele
                        v1(i)   =  nint(R(idx+i), dpI)
                    end do
                    idx =  idx+npe*n_ele

!                   ID_vtx.
                    do i=1,n_vtx
                        v2(i)   =  nint(R(idx+i), dpI)
                    end do
                    idx =  idx+n_vtx

                    call avg_add_patch(avg(iavg), .true., prc, ID, ele_type, n_ele, &
                        &  v1, n_vtx, v2, R(idx+1), ifac)
                    idx =  idx+3*n_vtx
                end if
            end do
        end do
!       allgather and synchronize.
!       ------------------------------------------------------------------------

        do i=1,n_avg
            if(avg(i)%nfacL .le. 0) cycle
            call patch_get_r(i, avg(i)%facL, avg(i)%nfacL)
!           call avg_output(i)

            idx =  avg(i)%idmL
            allocate(avg(i)%r_L(  idx-1), stat=err_mem)
            allocate(avg(i)%A_L(  idx-1), stat=err_mem)
            allocate(avg(i)%u  (3,idx-1), stat=err_mem)

            avg(i)%r_L  =  0.0d0
            avg(i)%A_L  =  0.0d0
            v2          =  0
            do isec=1,avg(i)%nfacL
                ltmp=  avg(i)%facL(isec)%myid .eq. myid
                if(ltmp)    allocate(avg(i)%facL(isec)%ID_r(avg(i)%facL(isec)%n_ele))

                do iele=1,avg(i)%facL(isec)%n_ele
                    do ipe=1,avg(i)%facL(isec)%npe
                        xyz(1:3,ipe)=  avg(i)%xyz(1:3,avg(iavg)%facL(isec)%n2e(ipe,iele))
                        v1 (    ipe)=  avg(i)%ID_r(   avg(iavg)%facL(isec)%n2e(ipe,iele))
                    end do
                    idx =  minval(v1(1:avg(i)%facL(isec)%npe))
                    if(ltmp)    avg(i)%facL(isec)%ID_r(iele)=  idx

                    if(avg(i)%facL(isec)%ele_type .eq. TRI_3) then
                        call crs_prd(xyz(1:3,2)-xyz(1:3,1), xyz(1:3,3)-xyz(1:3,1), n)
                        c   = (xyz(1:3,1)+xyz(1:3,2)+xyz(1:3,3))/3.0d0
                    else
                        call crs_prd(xyz(1:3,3)-xyz(1:3,1), xyz(1:3,4)-xyz(1:3,2), n)
                        c   = (xyz(1:3,1)+xyz(1:3,2)+xyz(1:3,3)+xyz(1:3,4))*0.25d0
                    end if
                    avg(i)%A_L(idx) =  avg(i)%A_L(idx)+0.5d0*sqrt(n(1)**2+n(2)**2+n(3)**2)

                    rtmp=  c(1)*z(1)+c(2)*z(2)+c(3)*z(3)
                    c   =  c-rtmp*z
                    avg(i)%r_L(idx) =  avg(i)%r_L(idx)+sqrt(c(1)**2+c(2)**2+c(3)**2)
                    v2        (idx) =  v2        (idx)+1
                end do
            end do

            do iele=1,avg(i)%idmL-1
                avg(i)%r_L(iele)=  avg(i)%r_L(iele)/real(v2(iele), dpR)
            end do
        end do

        lev =  0
        do i=1,n_avg
            do ID=1,avg(i)%nfacL
                if(allocated(avg(i)%facL(ID)%xyz   ))   deallocate(avg(i)%facL(ID)%xyz   )
                if(allocated(avg(i)%facL(ID)%ID_vtx))   deallocate(avg(i)%facL(ID)%ID_vtx)
                if(avg(i)%facL(ID)%myid .eq. myid)  cycle

                if(allocated(avg(i)%facL(ID)%n2e   ))   deallocate(avg(i)%facL(ID)%n2e   )
            end do
            lev =  lev+6+3*(avg(i)%idmL-1)
        end do
        call mpi_allreduce(lev, bct, 1, mpi_dpI, mpi_sum, mpi_comm_avg, mpi_err)
        allocate(abf1(max(lev, 1)), stat=err_mem)
        allocate(abf2(max(bct, 1)), stat=err_mem)

        if(allocated(v1))   deallocate(v1)
        if(allocated(v2))   deallocate(v2)
        if(allocated(S ))   deallocate(S )
        if(allocated(R ))   deallocate(R )

        return
        end subroutine avg_synchronize
!       ------------------------------------------------------------------------
!       output AVG.
!       ------------------------------------------------------------------------
        subroutine avg_output(iavg)
        use var_kind_def
        use var_cgns
        use var_parallel
        implicit none
        integer(dpI),intent(in):: iavg
        character(len=80):: str
        character(len=1000):: tec
        integer(dpI):: isec,iele,i,npe

        write(str,*),myid
        str =  './data/'//trim(adjustl(str))//'_avg_'
        write(tec,*),iavg
        str =  trim(adjustl(str))//trim(adjustl(tec))//'.dat'
        open(unit=10,file=trim(adjustl(str)))

        str =  'variables="CoordinateX","CoordinateY","CoordinateZ","R"'
        write(unit=10,fmt='(A)'),trim(str)

        do isec=1,avg(iavg)%nfacL
            tec =  'zone N='

            write(str,*),avg(iavg)%n_vtx_L
            tec =  trim(adjustl(tec))//trim(adjustl(str))

            write(str,*),avg(iavg)%facL(isec)%n_ele
            tec =  trim(adjustl(tec))//',E='//trim(adjustl(str))

            if(avg(iavg)%facL(isec)%ele_type .eq. QUAD_4) then
                str =  'FEQUADRILATERAL'
            elseif(avg(iavg)%facL(isec)%ele_type .eq. TRI_3) then
                str =  'FETRIANGLE'
            else
                stop 'Error: ele_type not supported.'
            end if
            tec =  trim(adjustl(tec))//','//'zonetype='//trim(adjustl(str))

            if(isec .gt. 1) then
                tec =  trim(adjustl(tec))//',varsharelist=([1-4]=1)'
            end if
            tec =  trim(adjustl(tec))//',datapacking=point'

            write(unit=10,fmt='(A)'),trim(adjustl(tec))
            if(isec .eq. 1) then
                do i=1,avg(iavg)%n_vtx_L
                    write(unit=10,fmt='(3ES20.12,I8)'),avg(iavg)%xyz(1:3,i), &
                        &  avg(iavg)%ID_r(i)
                end do
            end if
            npe =  avg(iavg)%facL(isec)%npe
            do iele=1,avg(iavg)%facL(isec)%n_ele
                write(unit=10,fmt='(8I8)'),avg(iavg)%facL(isec)%n2e(1:npe,iele)
            end do
        end do

        close(10)

        return
        end subroutine avg_output
    end module var_avg
