!-------------------------------------------------------------------------------
!   module rtree.
!-------------------------------------------------------------------------------
    module var_rtree
        implicit none
        private
        public:: type_rtree,rtree_setup,rtree_search,rtree_destroy,max_nearest_point

        integer, parameter:: dpL=  kind(.true.)
        integer, parameter:: dpI=  kind(1)
        integer, parameter:: dpR=  kind(1.0d0)
        integer(dpI),parameter:: max_nearest_point  =  5
        integer(dpI),parameter:: max_npnt   =  12
        real   (dpR),parameter:: max_size   =  1.0d-8

        type type_rtree
            integer(dpI):: dim  =  2
            integer(dpI):: npnt =  0
            real   (dpR),allocatable:: xyz(:,:)
            integer(dpI),allocatable:: idx(:)

            integer(dpI):: max_node =  0
            integer(dpI),allocatable:: node_LR(:,:),node_cdim(:)
            real   (dpR),allocatable:: node_xyz1(:,:),node_xyz0(:,:)
            real   (dpR),allocatable:: node_cL(:),node_cR(:),node_cut(:)
        end type type_rtree
        integer(dpI):: inode=  0

        contains
!       ------------------------------------------------------------------------
!       setup rtree.
!       ------------------------------------------------------------------------
        subroutine rtree_setup(dim,npnt,xyz,p)
        implicit none
        integer(dpI),intent(in):: dim,npnt
        real   (dpR),intent(in):: xyz(*)
        logical(dpL):: wrte
        integer(dpI):: i
        type(type_rtree):: p

        if((dim .le. 0) .or. (npnt .le. 0)) stop 'Error: wrong input for rtree_setup.'
        p%dim   =  dim
        p%npnt  =  npnt
        if(.not. allocated(p%xyz))  allocate(p%xyz(dim,npnt))
        if(.not. allocated(p%idx))  allocate(p%idx(    npnt))
        if(size(p%idx) .lt. npnt)   stop 'Error: size(idx) too small, rtree.'
        do i=1,npnt
            p%xyz(1:dim,i)  =  xyz(1+dim*(i-1):dim*i)
        end do
        forall(i=1:npnt)    p%idx(i)=  i

        p%max_node  = (3*p%npnt)/max_npnt+1000
        if(.not. allocated(p%node_LR  ))    allocate(p%node_LR  (2    ,p%max_node))
        if(.not. allocated(p%node_cdim))    allocate(p%node_cdim(      p%max_node))
        if(.not. allocated(p%node_xyz1))    allocate(p%node_xyz1(p%dim,p%max_node))
        if(.not. allocated(p%node_xyz0))    allocate(p%node_xyz0(p%dim,p%max_node))
        if(.not. allocated(p%node_cL  ))    allocate(p%node_cL  (      p%max_node))
        if(.not. allocated(p%node_cR  ))    allocate(p%node_cR  (      p%max_node))
        if(.not. allocated(p%node_cut ))    allocate(p%node_cut (      p%max_node))
        if(size(p%node_cL) .lt. p%max_node) stop 'Error: size(node_cL) too small, rtree.'
        p%node_LR   =  0
        p%node_cdim =  0
        p%node_xyz1 =  0.0d0
        p%node_xyz0 =  0.0d0
        p%node_cL   =  0.0d0
        p%node_cR   =  0.0d0
        p%node_cut  =  0.0d0

        wrte    =  .false.
        inode   =  1
        if(wrte) then
            open(unit=10,file='./lowest_level.dat')
            write(unit=10,fmt=*),'variables="x","y"'
            close(10)
        end if
        call rtree_build(1, p%npnt, p, wrte)
!       print*,inode,p%max_node

        return
        end subroutine rtree_setup
!       ------------------------------------------------------------------------
!       recursively build the tree.
!       ------------------------------------------------------------------------
        recursive subroutine rtree_build(L,R,p,wrte)
        implicit none
        logical(dpL),intent(in):: wrte
        integer(dpI),intent(in):: L,R
        type(type_rtree):: p
        integer(dpI):: LL,RR,i,j,k,m,cdim,node
        real   (dpR):: xyz1(p%dim),xyz0(p%dim),cut,wid

        node=  inode
        xyz1(1:p%dim)   =  p%xyz(1:p%dim,p%idx(L))
        xyz0(1:p%dim)   =  p%xyz(1:p%dim,p%idx(L))
        do i=L+1,R
            xyz1(1:p%dim)  =  min(xyz1(1:p%dim), p%xyz(1:p%dim,p%idx(i)))
            xyz0(1:p%dim)  =  max(xyz0(1:p%dim), p%xyz(1:p%dim,p%idx(i)))
        end do
        p%node_xyz1(1:p%dim,node)   =  xyz1(1:p%dim)
        p%node_xyz0(1:p%dim,node)   =  xyz0(1:p%dim)

        cdim=  maxloc(xyz0(1:p%dim)-xyz1(1:p%dim), 1)
        wid =  xyz0(cdim)-xyz1(cdim)
        if((R-L+1 .le. max_npnt) .or. (wid .le. max_size)) then
            p%node_LR(1,node)   =  L
            p%node_LR(2,node)   =  R
            if(wrte) then
                open(unit=10,file='./lowest_level.dat',status='old',position='append')
                write(unit=10,fmt=*),'zone I=',2,',J=',2
                write(unit=10,fmt=*),p%node_xyz1(1:2,node)
                write(unit=10,fmt=*),p%node_xyz0(1,node),p%node_xyz1(2,node)
                write(unit=10,fmt=*),p%node_xyz1(1,node),p%node_xyz0(2,node)
                write(unit=10,fmt=*),p%node_xyz0(1:2,node)
                close(10)
            end if
            return
        end if

!       ------------------------------------------------------------------------
!       cal cdim and sort.
        cut =  p%xyz(cdim,p%idx(L))
        do i=L+1,R
            cut =  cut+p%xyz(cdim,p%idx(i))
        end do
        cut =  cut/dble(R-L+1)

        i   =  L
        j   =  R
        do while(i .lt. j)
            if(p%xyz(cdim,p%idx(i)) .le. cut) then
                i   =  i+1
            else
                k   =  p%idx(j)
                p%idx(j)=  p%idx(i)
                p%idx(i)=  k
                j   =  j-1
            end if
        end do
        if(p%xyz(cdim,p%idx(i)) .le. cut) then
            M   =  i
        else
            M   =  i-1
        end if
!       cal cdim and sort.
!       ------------------------------------------------------------------------

        p%node_cdim(node)   =  cdim
        p%node_cut (node)   =  cut
        LL  =  0
        RR  =  0
        if(L .le. M) then
            inode   =  inode+1
            if(inode .gt. p%max_node)   stop 'Error: max_node not correct.'
            LL      =  inode
            call rtree_build(L  , M, p, wrte)
        end if
        if(M .lt. R) then
            inode   =  inode+1
            if(inode .gt. p%max_node)   stop 'Error: max_node not correct.'
            RR      =  inode
            call rtree_build(M+1, R, p, wrte)
        end if
        p%node_LR(1:2,node) =-(/LL, RR/)
        if(RR .eq. 0) then
            p%node_xyz1(:,node) =  p%node_xyz1(:,LL)
            p%node_xyz0(:,node) =  p%node_xyz0(:,LL)
            p%node_cL  (  node) =  p%node_xyz0(cdim,node)
            p%node_cut (  node) =  p%node_xyz0(cdim,node)
        elseif(LL .eq. 0) then
            p%node_xyz1(:,node) =  p%node_xyz1(:,RR)
            p%node_xyz0(:,node) =  p%node_xyz0(:,RR)
            p%node_cL  (  node) =  p%node_xyz1(cdim,node)
            p%node_cut (  node) =  p%node_xyz1(cdim,node)
        else
            p%node_cL  (  node) =  p%node_xyz0(cdim,LL)
            p%node_cR  (  node) =  p%node_xyz1(cdim,RR)
            p%node_cut (  node) =  0.5d0*(p%node_cL(node)+p%node_cR(node))
        end if

        return
        end subroutine rtree_build
!       ------------------------------------------------------------------------
!       find the nearest point.
!       ------------------------------------------------------------------------
        subroutine rtree_search(tree,p,dis0,nnp,np,dis)
        implicit none
        type(type_rtree),intent(in):: tree
        logical(dpL):: user_give_d0
        integer(dpI):: nnp,np(*)
        real   (dpR):: p(*),dis0,dis,rtmp

        if(tree%npnt .le. 0)    stop 'Error: empty tree.'
        rtmp        =  sum((p(1:tree%dim)-tree%xyz(1:tree%dim,1))**2)
        user_give_d0=  dis0 .lt. sqrt(rtmp)
        if(user_give_d0) then
            nnp     =  0
            dis     =  min(dis0*dis0, huge(1.0d0))
        else
            nnp     =  1
            np(1)   =  1
            dis     =  rtmp
        end if

        call search(1, nnp, np, p, dis, tree)
        dis =  sqrt(dis)
        np(nnp+1:max_nearest_point) = -1

        return
        end subroutine rtree_search
!       ------------------------------------------------------------------------
!       recursively find the nearest point using rtree.
!       ------------------------------------------------------------------------
        recursive subroutine search(inode,nvtx,vtx,p0,dis,p)
        implicit none
        integer(dpI),intent(in):: inode
        real   (dpR),intent(in):: p0(*)
        type(type_rtree),intent(in):: p
        logical(dpL):: ltmp
        integer(dpI):: nvtx,vtx(*),i,near,far
        real   (dpR):: dis,d(p%dim),rtmp

        ltmp= (p%node_LR(1,inode) .lt. 0) .or. (p%node_LR(2,inode) .lt. 0)
        if(.not. ltmp) then
            do i=p%node_LR(1,inode),p%node_LR(2,inode)
                if(i .eq. 0)    cycle
                d(1:p%dim)  =  p%xyz(1:p%dim,p%idx(i))-p0(1:p%dim)
                rtmp        =  sum(d*d)
                if(rtmp .lt. dis) then
                    nvtx    =  1
                    vtx(1)  =  p%idx(i)
                    dis     =  rtmp
                elseif(rtmp .eq. dis) then
                    if(nvtx .lt. max_nearest_point) then
                        nvtx        =  nvtx+1
                        vtx(nvtx)   =  p%idx(i)
                    end if
                end if
            end do
            return
        end if

        if(p0(p%node_cdim(inode)) .le. p%node_cut(inode)) then
            near=  p%node_LR(1,inode)
            far =  p%node_LR(2,inode)
            rtmp= (p0(p%node_cdim(inode))-p%node_cR(inode))**2
        else
            near=  p%node_LR(2,inode)
            far =  p%node_LR(1,inode)
            rtmp= (p0(p%node_cdim(inode))-p%node_cL(inode))**2
        end if
        if(near .lt. 0) call search(-near, nvtx, vtx, p0, dis, p)
        if(rtmp .gt. dis)   return
        if(far .ge. 0)  return
        do i=1,p%dim
            if(i .eq. p%node_cdim(inode))   cycle
            if(p0(i) .le. p%node_xyz1(i,inode)) then
                rtmp=  rtmp+(p%node_xyz1(i,inode)-p0(i))**2
            elseif(p0(i) .ge. p%node_xyz0(i,inode)) then
                rtmp=  rtmp+(p%node_xyz0(i,inode)-p0(i))**2
            else
!               do nothing.
            end if
            if(rtmp .gt. dis)   return
        end do
        if(far .lt. 0)  call search(-far , nvtx, vtx, p0, dis, p)

        return
        end subroutine search
!       ------------------------------------------------------------------------
!       delete the tree.
!       ------------------------------------------------------------------------
        subroutine rtree_destroy(p)
        implicit none
        type(type_rtree):: p

        if(allocated(p%xyz))    deallocate(p%xyz)
        if(allocated(p%idx))    deallocate(p%idx)
        if(allocated(p%node_LR))    deallocate(p%node_LR)
        if(allocated(p%node_cdim))    deallocate(p%node_cdim)
        if(allocated(p%node_xyz1))    deallocate(p%node_xyz1)
        if(allocated(p%node_xyz0))    deallocate(p%node_xyz0)
        if(allocated(p%node_cL))    deallocate(p%node_cL)
        if(allocated(p%node_cR))    deallocate(p%node_cR)
        if(allocated(p%node_cut))    deallocate(p%node_cut)
        p%dim       =  0
        p%npnt      = -1
        p%max_node  = -1

        return
        end subroutine rtree_destroy
    end module var_rtree
