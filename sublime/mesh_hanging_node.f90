!-------------------------------------------------------------------------------
!   get all faces related to hanging node.
!-------------------------------------------------------------------------------
    subroutine get_faces_hanging_node
    use var_kind_def
    use var_cgns
    use var_global, only: is_2d_cal,is_hanging_node,err_mem,n_dim
    use var_load_balance
    use var_mesh, match=>LR_hanging,n_match=>n_LR_hanging
    use var_parallel
    use var_per_bnd
    implicit none
    logical(dpL):: is_per,ltmp
    integer(dpI):: huge_1,iele,i,j,k,nfac,ele,v(8),npe,nfac_all,w(16),n_f2n, &
                &  iA_f(0:nprc),n_n2e
    integer(dpI),allocatable:: fac(:,:),fac_b(:,:),fac_all(:,:),f2n(:,:),n2e(:,:)

    if(.not. is_hanging_node)   return

    huge_1  =  huge(1)-1

!   ----------------------------------------------------------------------------
!   record the faces of all elements.
    nfac=  0
    do iele=1,n_ele_L
        nfac=  nfac+nface_ele(ele_type(iele))
    end do
    allocate(fac(5,nfac), stat=err_mem)

    nfac=  0
    do iele=1,n_ele_L
        ele =  ele_type(iele)
        i   =  iA_n2e(iele)

        if((ele .eq. TRI_3)) then
            v(1:3)  =  jA_n2e(i:i+2)
            fac(1:5,nfac+1) = (/v(1), v(2), huge_1, huge_1, iele+elei_prc(1,myid)-1/)
            fac(1:5,nfac+2) = (/v(2), v(3), huge_1, huge_1, iele+elei_prc(1,myid)-1/)
            fac(1:5,nfac+3) = (/v(3), v(1), huge_1, huge_1, iele+elei_prc(1,myid)-1/)
            nfac    =  nfac+3
        elseif((ele .eq. QUAD_4) .or. (ele .eq. QUAD_8) .or. (ele .eq. QUAD_9)) then
            v(1:4)  =  jA_n2e(i:i+3)
            fac(1:5,nfac+1) = (/v(1), v(2), huge_1, huge_1, iele+elei_prc(1,myid)-1/)
            fac(1:5,nfac+2) = (/v(2), v(3), huge_1, huge_1, iele+elei_prc(1,myid)-1/)
            fac(1:5,nfac+3) = (/v(3), v(4), huge_1, huge_1, iele+elei_prc(1,myid)-1/)
            fac(1:5,nfac+4) = (/v(4), v(1), huge_1, huge_1, iele+elei_prc(1,myid)-1/)
            nfac    =  nfac+4
        elseif(ele .eq. TETRA_4) then
            v(1:4)  =  jA_n2e(i:i+3)
            fac(1:5,nfac+1) = (/v(1), v(2), v(3), huge_1, iele+elei_prc(1,myid)-1/)
            fac(1:5,nfac+2) = (/v(1), v(2), v(4), huge_1, iele+elei_prc(1,myid)-1/)
            fac(1:5,nfac+3) = (/v(2), v(3), v(4), huge_1, iele+elei_prc(1,myid)-1/)
            fac(1:5,nfac+4) = (/v(1), v(3), v(4), huge_1, iele+elei_prc(1,myid)-1/)
            nfac    =  nfac+4
        elseif(ele .eq. PYRA_5) then
            v(1:5)  =  jA_n2e(i:i+4)
            fac(1:5,nfac+1) = (/v(1), v(2), v(3), v(4)  , iele+elei_prc(1,myid)-1/)
            fac(1:5,nfac+2) = (/v(1), v(2), v(5), huge_1, iele+elei_prc(1,myid)-1/)
            fac(1:5,nfac+3) = (/v(2), v(3), v(5), huge_1, iele+elei_prc(1,myid)-1/)
            fac(1:5,nfac+4) = (/v(3), v(4), v(5), huge_1, iele+elei_prc(1,myid)-1/)
            fac(1:5,nfac+5) = (/v(4), v(1), v(5), huge_1, iele+elei_prc(1,myid)-1/)
            nfac    =  nfac+5
        elseif(ele .eq. PENTA_6) then
            v(1:6)  =  jA_n2e(i:i+5)
            fac(1:5,nfac+1) = (/v(1), v(2), v(3), huge_1, iele+elei_prc(1,myid)-1/)
            fac(1:5,nfac+2) = (/v(4), v(5), v(6), huge_1, iele+elei_prc(1,myid)-1/)
            fac(1:5,nfac+3) = (/v(1), v(2), v(5), v(4), iele+elei_prc(1,myid)-1/)
            fac(1:5,nfac+4) = (/v(2), v(3), v(6), v(5), iele+elei_prc(1,myid)-1/)
            fac(1:5,nfac+5) = (/v(1), v(3), v(6), v(4), iele+elei_prc(1,myid)-1/)
            nfac    =  nfac+5
        elseif((ele .eq. HEXA_8) .or. (ele .eq. HEXA_20) .or. (ele .eq. HEXA_27)) then
            v(1:8)  =  jA_n2e(i:i+7)
            fac(1:5,nfac+1) = (/v(1), v(2), v(3), v(4), iele+elei_prc(1,myid)-1/)
            fac(1:5,nfac+2) = (/v(5), v(6), v(7), v(8), iele+elei_prc(1,myid)-1/)
            fac(1:5,nfac+3) = (/v(1), v(5), v(8), v(4), iele+elei_prc(1,myid)-1/)
            fac(1:5,nfac+4) = (/v(2), v(6), v(7), v(3), iele+elei_prc(1,myid)-1/)
            fac(1:5,nfac+5) = (/v(1), v(2), v(6), v(5), iele+elei_prc(1,myid)-1/)
            fac(1:5,nfac+6) = (/v(4), v(3), v(7), v(8), iele+elei_prc(1,myid)-1/)
            nfac    =  nfac+6
        else
            stop 'Error: element type not supported.'
        end if
    end do
    do i=1,nfac
        call iqsort(.true., 1, 4, fac(1,i))
    end do
    call iqsortcols(.true., 1, nfac, 1, 5, fac)
!   record the faces of all elements.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   delete the internal faces.
    do i=1,nfac
        if(fac(5,i) .le. 0) cycle
        v(1:4)  =  fac(1:4,i)
        if(is_2d_cal) then
            call find_internal_face(v, 2, nfac, fac, j, w)
        else
            call find_internal_face(v, 4, nfac, fac, j, w)
        end if
        if(j .le. 1) then
!           do nothing.
        elseif(j .gt. 2) then
            stop 'Error: more than 2 faces share the same vertex.'
        else
            fac(5,w(1)) = -abs(fac(5,w(1)))
            fac(5,w(2)) = -abs(fac(5,w(2)))
        end if
    end do
    j   =  nfac
    nfac=  0
    do i=1,j
        if(fac(5,i) .le. 0) cycle
        nfac            =  nfac+1
        fac(1:5,nfac)   =  fac(1:5,i)
    end do
!   delete the internal faces.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   delete the faces covered by physical boundary.
    allocate(fac_b(4,n_eleb_g), stat=err_mem)
    fac_b   =  huge_1
    do iele=1,n_eleb_g
        i   =  iA_n2e_b(iele)
        j   =  iA_n2e_b(iele+1)-1
        npe =  j-i+1
        if(is_2d_cal) then
            fac_b(1:2,iele) = (/jA_n2e_b(i:i+1)/)
        else
            if(npe .eq. 3) then
                fac_b(1:3,iele) = (/jA_n2e_b(i:i+2)/)
            elseif(npe .eq. 4) then
                fac_b(1:4,iele) = (/jA_n2e_b(i:i+3)/)
            elseif(npe .eq. 8) then
                fac_b(1:4,iele) = (/jA_n2e_b(i:i+3)/)
            elseif(npe .eq. 9) then
                fac_b(1:4,iele) = (/jA_n2e_b(i:i+3)/)
            else
                stop 'Error: boundary element type not supported.'
            end if
        end if
        call iqsort(.true., 1, 4, fac_b(1,iele))
    end do
    call iqsortcols(.true., 1, n_eleb_g, 1, 4, fac_b)

    j   =  nfac
    nfac=  0
    do iele=1,j
        call ib_search(1, n_eleb_g, 1, 4, fac_b, fac(1,iele), ltmp, v(1), v(2))
        i   =  0
        if(ltmp) then
        do k=v(1),v(2)
            if((fac(2,iele) .eq. fac_b(2,k)) .and. (fac(3,iele) .eq. fac_b(3,k)) .and. &
              &(fac(4,iele) .eq. fac_b(4,k)))   i   =  i+1
        end do
        end if

        if(i .gt. 1) then
            stop 'Error: more than 1 boundary face match the internal element.'
        elseif(i .eq. 1) then
!           do nothing.
        else
            nfac            =  nfac+1
            fac(1:5,nfac)   =  fac(1:5,iele)
        end if
    end do
    if(allocated(fac_b))    deallocate(fac_b)
!   delete the faces covered by physical boundary.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   delete the faces covered by periodic boundary.
    if(n_vtx_pair_per .gt. 0) then
    k   =  nfac
    nfac=  0
    do i=1,k
        is_per  =  .true.
        do j=1,4
            if(fac(j,i) .lt. huge_1) then
                call ib_search(1,n_vtx_pair_per,1,5,vtx_pair_per,fac(j,i),ltmp,v(1),v(2))
                if(.not. ltmp) then
                    is_per  =  .false.
                    exit
                end if
            else
                exit
            end if
        end do
        if(is_per) then
!           do nothing.
        else
            nfac            =  nfac+1
            fac(1:5,nfac)   =  fac(1:5,i)
        end if
    end do
    end if
!   delete the faces covered by periodic boundary.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   synchronize all the remaining faces.
    if(allocated(fac_b))    deallocate(fac_b)
    call mpi_allgather(nfac, 1, mpi_dpI, prc_info, 1, mpi_dpI, mpi_comm_world, mpi_err)
    iA_f(0) =  1
    do i=0,nprc-1
        iA_f(i+1)   =  iA_f(i)+prc_info(i,1)
    end do
    nfac_all=  iA_f(nprc)-1
    if(nfac_all .le. 0) then
        is_hanging_node =  .false.
        return
    end if
    allocate(fac_all(5, nfac_all), stat=err_mem)

    nfac_all=  nfac*5
    call i_allgather(mpi_comm_world, nprc, nfac_all, fac, fac_all)
    nfac_all=  nfac_all/5
    call iqsortcols(.true., 1, nfac_all, 1, 5, fac_all)
    if(nfac_all .le. 0) return
!   synchronize all the remaining faces.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   delete duplicate faces in case of parallel computation.
    if(nprc .gt. 1) then
        do i=1,nfac_all
            if(fac_all(5,i) .lt. 0) cycle
            v(1:4)  =  fac_all(1:4,i)
            if(is_2d_cal) then
                call find_internal_face(v, 2, nfac_all, fac_all, j, w)
            else
                call find_internal_face(v, 4, nfac_all, fac_all, j, w)
            end if
            if(j .ne. 2)    cycle
            fac_all(5,w(1)) = -1
            fac_all(5,w(2)) = -1
        end do
        j       =  nfac_all
        nfac_all=  0
        do i=1,j
            if(fac_all(5,i) .le. 0) cycle
            nfac_all=  nfac_all+1
            fac_all(1:5,nfac_all)   =  fac_all(1:5,i)
        end do
    end if
!   delete duplicate faces in case of parallel computation.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   for 2D case, the coordinate info is required to assist the procedure.
    if(n_dim .eq. 2) then
        n_xyz_i_need=  0
        do i=1,nfac_all
        do j=1,4
            if(fac_all(j,i) .lt. huge_1) then
                n_xyz_i_need=  n_xyz_i_need+1
            end if
        end do
        end do
        if(allocated(xyz_i_need)) then
            if(n_xyz_i_need .gt. size(xyz_i_need)) then
                deallocate(xyz_i_need)
                allocate(xyz_i_need(n_xyz_i_need), stat=err_mem)
            end if
        else
            allocate(xyz_i_need(n_xyz_i_need), stat=err_mem)
        end if
        n_xyz_i_need=  0
        do i=1,nfac_all
        do j=1,4
            if(fac_all(j,i) .lt. huge_1) then
                n_xyz_i_need=  n_xyz_i_need+1
                xyz_i_need(n_xyz_i_need)=  fac_all(j,i)
            end if
        end do
        end do
        call simplify_series(n_xyz_i_need, 1, 1, xyz_i_need)
        call get_xyz_i_need
    end if
!   for 2D case, the coordinate info is required to assist the procedure.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   we need the face-to-node.
    n_f2n   =  0
    do i=1,nfac_all
        k   =  0
        do j=1,4
            if(fac_all(j,i) .eq. huge_1) then
                exit
            else
                k   =  k+1
            end if
        end do

        if((k .ge. 2) .and. (k .le. 4)) then
            n_f2n   =  n_f2n+k
        else
            stop 'Error: fails to get face_all.'
        end if
    end do
    allocate(f2n(2,n_f2n), stat=err_mem)
    n_f2n   =  0
    do i=1,nfac_all
        k   =  0
        do j=1,4
            if(fac_all(j,i) .eq. huge_1) then
                exit
            else
                k   =  k+1
            end if
        end do

        do j=1,k
            n_f2n   =  n_f2n+1
            f2n(1:2,n_f2n)  = (/fac_all(j,i), i/)
        end do
    end do
    call iqsortcols(.true., 1, n_f2n, 1, 2, f2n)
!   we need the face-to-node.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   find the right-part.
    do i=1,nfac
        k   =  0
        do j=1,4
            if(fac(j,i) .eq. huge_1) then
                exit
            else
                k   =  k+1
            end if
        end do

        call find_internal_face(fac(1,i), k, nfac_all, fac_all, j, w)
        if(j .ne. 1)    cycle

        if(k .eq. 2) then
            call find_bar(i, w(1), w(2), fac, fac_all, n_f2n, f2n)
            if(w(2) .gt. 0) then
                fac_all(5,w(1)) = -abs(fac_all(5,w(1)))
                fac_all(5,w(2)) = -abs(fac_all(5,w(2)))
                fac_all(5,w(3)) = -abs(fac_all(5,w(3)))
                fac(1:3,i)      = -w(1:3)
                fac(4:5,i)      =  0
            end if
        elseif(k .eq. 3) then
            call find_tri(i, w(1), w(2), fac, fac_all, n_f2n, f2n)
            if(w(2) .gt. 0) then
                fac_all(5,w(1)) = -abs(fac_all(5,w(1)))
                fac_all(5,w(2)) = -abs(fac_all(5,w(2)))
                fac_all(5,w(3)) = -abs(fac_all(5,w(3)))
                fac_all(5,w(4)) = -abs(fac_all(5,w(4)))
                fac_all(5,w(5)) = -abs(fac_all(5,w(5)))
                fac(1:5,i)      = -w(1:5)
            end if
        elseif(k .eq. 4) then
            call find_quad(i, w(1), w(2), fac, fac_all, n_f2n, f2n)
            if(w(2) .gt. 0) then
                fac_all(5,w(1)) = -abs(fac_all(5,w(1)))
                fac_all(5,w(2)) = -abs(fac_all(5,w(2)))
                fac_all(5,w(3)) = -abs(fac_all(5,w(3)))
                fac_all(5,w(4)) = -abs(fac_all(5,w(4)))
                fac_all(5,w(5)) = -abs(fac_all(5,w(5)))
                fac(1:5,i)      = -w(1:5)
            end if
        else
            stop 'Error: fails to get faces with hanging nodes.'
        end if
    end do

    if(allocated(f2n))  deallocate(f2n)

    k   =  nfac
    nfac=  0
    do i=1,k
        if(fac(1,i) .ge. 0) cycle
        nfac            =  nfac+1
        fac(1:5,nfac)   =  abs(fac(1:5,i))
    end do
!   find the right-part.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   synchronize all the matched records and check.
    call mpi_allreduce(nfac, n_match, 1, mpi_dpI, mpi_sum, mpi_comm_world, mpi_err)
    allocate(match(5,n_match), stat=err_mem)
    match   =  0

    n_match =  5*nfac
    call i_allgather(mpi_comm_world, nprc, n_match, fac, match)
    n_match =  n_match/5

    do i=1,nfac_all
        fac_all(5,i)= -abs(fac_all(5,i))
    end do
    do i=1,n_match
    do k=1,5
        if(match(k,i) .gt. 0) then
            j           =  match(k,i)
            fac_all(5,j)=  abs(fac_all(5,j))
        else
            exit
        end if
    end do
    end do
    do i=1,nfac_all
        if(fac_all(5,i) .le. 0) then
            print*,'Error: a face has not find the matching part:',myid,i
            stop
        end if
    end do
!   synchronize all the matched records and check.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   get the hanging_node-to-left_element connectivity. 
    n_n2e   =  0
    do i=1,n_match
        if(fac_all(5,match(1,i)) .lt. elei_prc(1,myid)) cycle
        if(fac_all(5,match(1,i)) .gt. elei_prc(2,myid)) cycle
        npe =  0
        do j=1,4
            if(fac_all(j,match(1,i)) .lt. huge_1) then
                npe     =  npe+1
            else
                exit
            end if
        end do
        if(npe .eq. 2) then
            n_n2e   =  n_n2e+2
        elseif(npe .eq. 3) then
            n_n2e   =  n_n2e+3
        else
            n_n2e   =  n_n2e+5
        end if
    end do
    allocate(n2e(2,max(n_n2e, 1)))
    n_n2e   =  0
    do i=1,n_match
        if(fac_all(5,match(1,i)) .lt. elei_prc(1,myid)) cycle
        if(fac_all(5,match(1,i)) .gt. elei_prc(2,myid)) cycle
        npe =  0
        do j=1,4
            if(fac_all(j,match(1,i)) .lt. huge_1) then
                npe     =  npe+1
            else
                exit
            end if
        end do

        if(npe .eq. 2) then
            w(1 :2 )=  fac_all(1:2,match(2,i))
            w(3 :4 )=  fac_all(1:2,match(3,i))
            k       =  4
        elseif(npe .eq. 3) then
            w(1 :3 )=  fac_all(1:3,match(2,i))
            w(4 :6 )=  fac_all(1:3,match(3,i))
            w(7 :9 )=  fac_all(1:3,match(4,i))
            w(10:12)=  fac_all(1:3,match(5,i))
            k       =  12
        else
            w(1 :4 )=  fac_all(1:4,match(2,i))
            w(5 :8 )=  fac_all(1:4,match(3,i))
            w(9 :12)=  fac_all(1:4,match(4,i))
            w(13:16)=  fac_all(1:4,match(5,i))
            k       =  16
        end if
        call list_minus_list(k, w, npe, fac_all(1,match(1,i)))
        if(npe .eq. 2) then
            if(k .ne. 1)    stop 'Error: fails to find the hanging node.'
        elseif(npe .eq. 3) then
            if(k .ne. 3)    stop 'Error: fails to find the hanging node.'
        else
            if(k .ne. 5)    stop 'Error: fails to find the hanging node.'
        end if

        do j=1,k
            n_n2e           =  n_n2e+1
            n2e(1:2,n_n2e)  = (/w(j), fac_all(5,match(1,i))/)
        end do
    end do
!   get the hanging_node-to-left_element connectivity. 
!   ----------------------------------------------------------------------------

    do i=1,n_match
    do k=1,5
        if(match(k,i) .gt. 0) then
            j           =  match(k,i)
            match(k,i)  =  abs(fac_all(5,j))
        else
            exit
        end if
    end do
    end do
    call iqsortcols(.true., 1, n_match, 1, 5, match)

    if(allocated(fac    ))  deallocate(fac    )
    if(allocated(fac_b  ))  deallocate(fac_b  )
    if(allocated(fac_all))  deallocate(fac_all)
    if(allocated(f2n    ))  deallocate(f2n    )

    if(n_n2e .gt. 0)    call add_hanging_n2e(n_n2e, n2e)
    if(allocated(n2e))  deallocate(n2e)

    return
    contains
!       ------------------------------------------------------------------------
!       find the internal face.
!       ------------------------------------------------------------------------
        subroutine find_internal_face(LL,LDA,nfac,fac,n_f,f)
        implicit none
        integer(dpI),intent(in):: LL(*),LDA,nfac,fac(5,*)
        logical(dpL):: ltmp
        integer(dpI):: n_f,f(*),L,R,i,j

        n_f =  0
        call ib_search(1, nfac, 1, 5, fac, LL(1), ltmp, L, R)
        if(.not. ltmp)  return
        do i=L,R
            ltmp=  .false.
            do j=2,LDA
                ltmp=  fac(j,i) .ne. LL(j)
                if(ltmp)    exit
            end do
            if(ltmp)    cycle
            n_f     =  n_f+1
            f(n_f)  =  i
        end do

        return
        end subroutine find_internal_face
!       ------------------------------------------------------------------------
!       find the fine-mesh part of non-matching face, QUAD_4.
!       ------------------------------------------------------------------------
        subroutine find_quad(iele,L_face,LR,fac,f,n_f2n,f2n)
        use var_kind_def
        implicit none
        integer(dpI),intent(in):: iele,L_face
        logical(dpL):: ltmp
        integer(dpI):: LR(*),vL(4),v(16),freq(16),w(2,16),i,j,k,L,M,R,fR(0:100,4)
        integer(dpI):: n_f2n,fac(5,*),f(5,*),f2n(2,*)

        LR(1:4) =  0
        if(fac(5,iele) .le. 0)  return
        vL(1:4) =  fac(1:4,iele)

!       ------------------------------------------------------------------------
!       record the candidate faces.
        fR  =  0
        do i=1,4
            call ib_search(1, n_f2n, 1, 2, f2n, vL(i), ltmp, L, R)
            if(ltmp) then
                do j=L,R
                    k   =  f2n(2,j)
                    if((f(5,k) .le. 0) .or. (k .eq. L_face))    cycle
                    fR(0      ,i)   =  fR(0,i)+1
                    fR(fR(0,i),i)   =  k
                end do
            end if
        end do
!       record the candidate faces.
!       ------------------------------------------------------------------------

        do L=1,fR(0,4)
        if(f(5,fR(L,4)) .le. 0) cycle
        do k=1,fR(0,3)
        if(f(5,fR(k,3)) .le. 0) cycle
        do j=1,fR(0,2)
        if(f(5,fR(j,2)) .le. 0) cycle
        do i=1,fR(0,1)
            if(f(5,fR(i,1)) .le. 0) cycle
            v(1 :4 )=  f(1:4,fR(i,1))
            v(5 :8 )=  f(1:4,fR(j,2))
            v(9 :12)=  f(1:4,fR(k,3))
            v(13:16)=  f(1:4,fR(L,4))
            M       =  16
            call simplify_list_frequency(M, 1, 1, v, freq)
            if(M .ne. 9)    cycle
            do M=1,9
                w(1:2,M)= (/v(M), freq(M)/)
            end do
            call iqsortcols(.true., 1, 9, 2, 2, w)
            if((w(2,4) .ne. 1) .or. (w(2,8) .ne. 2) .or. (w(2,9) .ne. 4))   cycle
            call iqsortcols(.true., 1, 4, 1, 2, w)
            if(w(1,1) .ne. vL(1))   cycle
            if(w(1,2) .ne. vL(2))   cycle
            if(w(1,3) .ne. vL(3))   cycle
            if(w(1,4) .ne. vL(4))   cycle

            LR(1:4) = (/fR(i,1), fR(j,2), fR(k,3), fR(L,4)/)
            return
        end do
        end do
        end do
        end do

        return
        end subroutine find_quad
!       ------------------------------------------------------------------------
!       find the fine-mesh part of non-matching face, TRI_3.
!       ------------------------------------------------------------------------
        subroutine find_tri(iele,L_face,LR,fac,f,n_f2n,f2n)
        use var_kind_def
        implicit none
        integer(dpI),intent(in):: iele,L_face
        logical(dpL):: ltmp
        integer(dpI):: LR(*),vL(3),v(9),freq(9),w(2,9),i,j,k,L,M,R,fR(0:100,4)
        integer(dpI):: n_f2n,fac(5,*),f(5,*),f2n(2,*)

        LR(1:4) =  0
        if(fac(5,iele) .le. 0)  return
        vL(1:3) =  fac(1:3,iele)

!       ------------------------------------------------------------------------
!       record the candidate faces.
        fR  =  0
        do i=1,3
            call ib_search(1, n_f2n, 1, 2, f2n, vL(i), ltmp, L, R)
            if(ltmp) then
                do j=L,R
                    k   =  f2n(2,j)
                    if((f(5,k) .le. 0) .or. (k .eq. L_face))    cycle
                    fR(0      ,i)   =  fR(0,i)+1
                    fR(fR(0,i),i)   =  k
                end do
            end if
        end do
!       record the candidate faces.
!       ------------------------------------------------------------------------

        do k=1,fR(0,3)
        if(f(5,fR(k,3)) .le. 0) cycle
        do j=1,fR(0,2)
        if(f(5,fR(j,2)) .le. 0) cycle
        do i=1,fR(0,1)
            if(f(5,fR(i,1)) .le. 0) cycle
            v(1:3)  =  f(1:3,fR(i,1))
            v(4:6)  =  f(1:3,fR(j,2))
            v(7:9)  =  f(1:3,fR(k,3))
            M       =  9
            call simplify_list_frequency(M, 1, 1, v, freq)
            if(M .ne. 6)    cycle
            do M=1,6
                w(1:2,M)= (/v(M), freq(M)/)
            end do
            call iqsortcols(.true., 1, 6, 2, 2, w)
            if((w(2,3) .ne. 1) .or. (w(2,6) .ne. 2))    cycle

            call iqsortcols(.true., 1, 3, 1, 2, w)
            if(w(1,1) .ne. vL(1))   cycle
            if(w(1,2) .ne. vL(2))   cycle
            if(w(1,3) .ne. vL(3))   cycle

            call ib_search(1, nfac_all, 1, 5, fac_all, w(1,4), ltmp, v(1), v(2))
            if(.not. ltmp)  cycle
            do M=v(1),v(2)
                if(fac_all(5,M) .le. 0) cycle
                if(fac_all(5,M) .eq. fac_all(5,fR(i,1)))  cycle
                if(fac_all(5,M) .eq. fac_all(5,fR(j,2)))  cycle
                if(fac_all(5,M) .eq. fac_all(5,fR(k,3)))  cycle
                ltmp= (w(1,5) .eq. fac_all(2,M)) .and. (w(1,6) .eq. fac_all(3,M)) .and. &
                    & (fac_all(4,M) .eq. huge_1)
                if(.not. ltmp)  cycle
                LR(1:4) = (/fR(i,1), fR(j,2), fR(k,3), M/)
                return
            end do
        end do
        end do
        end do

        return
        end subroutine find_tri
!       ------------------------------------------------------------------------
!       find the fine-mesh part of non-matching face, BAR_2.
!       ------------------------------------------------------------------------
        subroutine find_bar(iele,L_face,LR,fac,f,n_f2n,f2n)
        use var_kind_def
        implicit none
        integer(dpI),intent(in):: iele,L_face
        real   (dpR),parameter:: eps=1.0d-6
        logical(dpL):: ltmp
        integer(dpI):: LR(*),vL(2),v(6),i,j,k,L,R,fR(0:100,2)
        integer(dpI):: n_f2n,fac(5,*),f(5,*),f2n(2,*)
        real   (dpR):: xyz(2,6),b(2),LL,rtmp

        LR(1:2) =  0
        if(fac(5,iele) .le. 0)  return
        vL(1:2) =  fac(1:2,iele)

!       ------------------------------------------------------------------------
!       record the candidate faces.
        fR  =  0
        do i=1,2
            call ib_search(1, n_f2n, 1, 2, f2n, vL(i), ltmp, L, R)
            if(ltmp) then
                do j=L,R
                    k   =  f2n(2,j)
                    if((f(5,k) .le. 0) .or. (k .eq. L_face))    cycle
                    fR(0      ,i)   =  fR(0,i)+1
                    fR(fR(0,i),i)   =  k
                end do
            end if
        end do
!       record the candidate faces.
!       ------------------------------------------------------------------------

        loop:do j=1,fR(0,2)
        do i=1,fR(0,1)
            v(1:2)  =  f(1:2,fR(i,1))
            v(3:4)  =  f(1:2,fR(j,2))
            k   =  4
            call simplify_series(k, 1, 1, v)
            if(k .ne. 3)    cycle
            do k=1,3
                if((v(k) .eq. vL(1)) .or. (v(k) .eq. vL(2)))    v(k)=  0
            end do
            L   =  0
            do k=1,3
                if(v(k) .gt. 0) L   =  L+1
            end do
            if(L .eq. 1) then
                LR(1:2) = (/fR(i,1), fR(j,2)/)
                exit loop
            end if
        end do
        end do loop
        if(LR(1) .le. 0)    return
        v(1:2)  =  vL(1:2)
        v(3:4)  =  f(1:2,LR(1))
        v(5:6)  =  f(1:2,LR(2))
        do i=1,6
            call ib_search(1, n_xyz_i_need, 1, 1, xyz_i_need, v(i), ltmp, L, R)
            if((.not. ltmp) .or. (L .ne. R))    stop 'Error: find_bar fails.'
            xyz(1:2,i)  =  xyz_need(2*L-1:2*L)
        end do
        b(1:2)  =  xyz(1:2,2)-xyz(1:2,1)
        call norm_vec(2, b, LL)
        do i=3,6
            rtmp=((xyz(1,i)-xyz(1,1))*b(1)+(xyz(2,i)-xyz(2,1))*b(2))/LL
            if((rtmp .le. -eps) .or. (rtmp .ge. 1.0d0+eps)) then
                LR(1:2) =  0
                return
            end if
        end do

        return
        end subroutine find_bar
!       ------------------------------------------------------------------------
!       add the hanging_node-to-left_face connectivity to the n2e information.
!       ------------------------------------------------------------------------
        subroutine add_hanging_n2e(n_n2e,n2e)
        use var_kind_def
        implicit none
        integer(dpI),intent(in):: n_n2e,n2e(2,*)
        integer(dpI):: i,j,k,iele,n_ele,nnz,ivtx
        integer(dpI),allocatable:: iA(:),jA(:),npe(:)

        n_ele   =  elei_prc(2,myid)-elei_prc(1,myid)+1
        allocate(npe(n_ele))
        npe     =  0
        do i=1,n_ele
            npe(i)  =  iA_n2e(i+1)-iA_n2e(i)
        end do

        if(n_n2e .gt. 0)    call iqsortcols(.true., 1, n_n2e, 2, 2, n2e)
        do i=1,n_n2e
            iele        =  n2e(2,i)
            if((iele .lt. elei_prc(1,myid)) .or. (iele .gt. elei_prc(2,myid))) &
                &  stop 'Error: fails to add n2e of hanging node.'
            iele        =  iele-elei_prc(1,myid)+1
            npe(iele)   =  npe(iele)+1
        end do
        nnz =  0
        do i=1,n_ele
            nnz =  nnz+npe(i)
        end do
        allocate(iA(n_ele+1), stat=err_mem)
        allocate(jA(nnz    ), stat=err_mem)
        iA(1)   =  1
        jA      =  0
        do i=1,n_ele
            iA(i+1) =  iA(i)+npe(i)
        end do
        do i=1,n_ele
            call ICOPY(iA_n2e(i+1)-iA_n2e(i), jA_n2e(iA_n2e(i)), 1, jA(iA(i)), 1)
        end do
        do j=1,n_n2e
            ivtx=  n2e(1,j)
            i   =  n2e(2,j)-elei_prc(1,myid)+1
            ltmp=  .false.
            do k=iA(i),iA(i+1)-1
                ltmp=  jA(k) .le. 0
                if(ltmp) then
                    jA(k)   =  ivtx
                    exit
                end if
            end do
            if(.not. ltmp)  stop 'Error: fails to add hanging n2e to n2e.'
        end do
        call ICOPY(n_ele+1, iA, 1, iA_n2e, 1)
        if(allocated(jA_n2e))   deallocate(jA_n2e)
        allocate(jA_n2e(iA_n2e(n_ele+1)-1), stat=err_mem)
        call ICOPY(iA_n2e(n_ele+1)-1, jA, 1, jA_n2e, 1)

        if(allocated(iA ))  deallocate(iA )
        if(allocated(jA ))  deallocate(jA )
        if(allocated(npe))  deallocate(npe)

        return
        end subroutine add_hanging_n2e
    end subroutine get_faces_hanging_node
!-------------------------------------------------------------------------------
!   reset the volume of with hanging node
!-------------------------------------------------------------------------------
    subroutine reset_vol_hanging_node(lev)
    use var_kind_def
    use var_cgns
    use var_global,only: n_dim,R13
    use var_mesh,only: mesh,sec
    implicit none
    integer(dpI),intent(in):: lev
    integer(dpI):: im,isec,sL,eL,sR,eR,npe,v(8),i,ele,m2f(4)
    real   (dpR):: xyz(3,8),n(4),cen0(3),cenL(3),cenR(3),dL(3),dR(3),vL,vR,rtm
    logical(dpL):: ltm

    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(sec(isec)%is_bnd)    cycle
        if(.not. allocated(sec(isec)%vol))  allocate(sec(isec)%vol(sec(isec)%n_ele))
        sec(isec)%vol   =  0.0d0
    enddo

    if(n_dim .eq. 3) then
        do im=1,mesh(lev)%n_mortar
            ele =  mesh(lev)%mortar_ele_type(im)
            npe =  npe_ele(ele)
!           m2f(1:npe)  =  mesh(lev)%mortar_m2f(1:npe,im)
            ltm =  .false.
            do i=1,npe
                if(m2f(i).le.0) cycle
                ltm =  .true.
                exit
            end do
            if(.not.ltm) cycle
            v(1:npe)=  mesh(lev)%mortar_n2e(1:npe,im)
            forall(i=1:npe) xyz(1:n_dim,i)  =  mesh(lev)%xyz(1:n_dim,v(i))
            do i=1,npe
                if(m2f(i).le.0) cycle
                xyz(1:n_dim,i+4)    =  mesh(lev)%xyz(1:n_dim,m2f(i))
            end do
            if(ele .eq. TRI_3) then
                cen0(1:3)   = (xyz(1:3,1)+xyz(1:3,2)+xyz(1:3,3))*R13
            elseif(ele .eq. QUAD_4) then
                cen0(1:3)   = (xyz(1:3,1)+xyz(1:3,2)+xyz(1:3,3)+xyz(1:3,4))*0.25d0
            else
                stop 'only 3D mortar_n need to be modified'
            end if
            n   =  0d0
            do i=1,npe
                call cal_normal_m2f(i,n)
            end do
            n   =  n*0.5d0
            call norm_vec(3, n, n(4))
            rtm =  n(1)*mesh(lev)%mortar_n_vg(1,im)+n(2)*mesh(lev)%mortar_n_vg(2,im) &
                & +n(3)*mesh(lev)%mortar_n_vg(3,im)
            if(rtm .ge. 0d0) then
                mesh(lev)%mortar_n_vg(1:4,im)  =  n(1:4)
            else
                mesh(lev)%mortar_n_vg(1:4,im)  = -n(1:4)
            end if
        end do
!       if(allocated(mesh(lev)%mortar_m2f)) deallocate(mesh(lev)%mortar_m2f)
    end if
    call mesh_check(lev)

    do im=1,mesh(lev)%n_mortar
        sL  =  mesh(lev)%mortar_LR(1,im)
        eL  =  mesh(lev)%mortar_LR(2,im)
        sR  =  mesh(lev)%mortar_LR(3,im)
        eR  =  mesh(lev)%mortar_LR(4,im)

        cen0(1:3)   =  mesh(lev)%mortar_cen(1:n_dim,im)
        n(1:3)      =  mesh(lev)%mortar_n_vg(1:3,im)*mesh(lev)%mortar_n_vg(4,im)
        if(sec(sL)%is_int) then
            cenL(1:3)   =  sec(sL)%cen(1:3,eL)
            dL(1:3)     =  cen0(1:3)-cenL(1:3)
            vL          =  (dL(1)*n(1)+dL(2)*n(2)+dL(3)*n(3))/dble(n_dim)
            sec(sL)%vol(eL) =  sec(sL)%vol(eL)+vL
        end if
        if(sec(sR)%is_int) then
            cenR(1:3)   =  sec(sR)%cen(1:3,eR)
            dR(1:3)     =  cenR(1:3)-cen0(1:3)
            vR          =  (dR(1)*n(1)+dR(2)*n(2)+dR(3)*n(3))/dble(n_dim)
            sec(sR)%vol(eR) =  sec(sR)%vol(eR)+vR
        end if
    end do
    call reset_vol_hanging_node_parallel(lev)
    do im=1,mesh(lev)%n_mortar
        sL  =  mesh(lev)%mortar_LR(1,im)
        eL  =  mesh(lev)%mortar_LR(2,im)
        sR  =  mesh(lev)%mortar_LR(3,im)
        eR  =  mesh(lev)%mortar_LR(4,im)

        vL  =  sec(sL)%vol(eL)
        if(im .le. mesh(lev)%n_mortar_b) then
            vR  =  vL
        else
            vR  =  sec(sR)%vol(eR)
        end if
        mesh(lev)%mortar_n_vg(5,im) =  mesh(lev)%mortar_n_vg(4,im)/(0.5d0*(vL+vR))
    end do

    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(sec(isec)%is_bnd)    cycle
        do eL=1,sec(isec)%n_ele
            vL  =  sec(isec)%vol(eL)
            if(vL.le.0d0) then
                npe =  sec(isec)%npe
                v(1:npe)=  sec(isec)%n2e(1:npe,eL)
                forall(i=1:npe) xyz(1:n_dim,i)  =  mesh(lev)%xyz(1:n_dim,v(i))
                print*, vL
                print*, xyz(1:n_dim,1:npe)
                stop 'Error: non-posotive volume.'
            end if
        end do
    end do

    return
    contains
!   ------------------------------------------------------------------------
!   cal normal with middle vertex of face
!   ------------------------------------------------------------------------
    subroutine cal_normal_m2f(i,n)
    implicit none
    integer(dpI):: i,j
    real   (dpR):: n(*),d1(1:3),d2(1:3),d3(1:3),n1(3)

    if(i .eq. npe) then
        j   =  1
    else
        j   =  i+1
    end if
    d1(1:3) =  xyz(1:3,i)-cen0(1:3)
    d3(1:3) =  xyz(1:3,j)-cen0(1:3)
    if(m2f(i).le.0) then
        call crs_prd(d1, d3, n1)
        n(1:3)  =  n(1:3)+n1(1:3)
    else
        d2(1:3) =  xyz(1:3,i+4)-cen0(1:3)
        call crs_prd(d1, d2, n1)
        n(1:3)  =  n(1:3)+n1(1:3)
        call crs_prd(d2, d3, n1)
        n(1:3)  =  n(1:3)+n1(1:3)
    end if

    end subroutine cal_normal_m2f
    end subroutine reset_vol_hanging_node
!-------------------------------------------------------------------------------
!   get vol for FV in the parallel environment on lev(lev) with hanging node
!-------------------------------------------------------------------------------
    subroutine reset_vol_hanging_node_parallel(lev)
    use var_kind_def
    use var_mesh
    use var_parallel
    implicit none
    integer(dpI),intent(in):: lev
    integer(dpI):: isr,isec,iele,i,idx,ip_remote

!   ----------------------------------------------------------------------------
!   prepare data.
    do isr=mesh(lev)%p2p_1,mesh(lev)%p2p_0
        idx =  1
        do i=1,p2p(isr)%n_ele_send
            isec=  p2p(isr)%id_ele_send(1,i)
            iele=  p2p(isr)%id_ele_send(2,i)
            p2p(isr)%rsend(idx) =  sec(isec)%vol(iele)
            idx =  idx+1
        end do
        p2p(isr)%n_send =  idx-1
        p2p(isr)%n_recv =  p2p(isr)%n_ele_recv
    end do
!   prepare data.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   data exchange.
    mpi_nreq=  0
    do isr=mesh(lev)%p2p_1,mesh(lev)%p2p_0
        ip_remote   =  p2p(isr)%ip_remote
        if(myid .eq. ip_remote) cycle

        if(p2p(isr)%n_send .gt. 0) then
            mpi_nreq=  mpi_nreq+1
            call mpi_isend(p2p(isr)%rsend, p2p(isr)%n_send, mpi_dpR, ip_remote, &
                &  myid, mpi_comm_world, mpi_req(mpi_nreq), mpi_err)
        end if

        if(p2p(isr)%n_recv .gt. 0) then
            mpi_nreq=  mpi_nreq+1
            call mpi_irecv(p2p(isr)%rrecv, p2p(isr)%n_recv, mpi_dpR, ip_remote, &
                &  ip_remote, mpi_comm_world, mpi_req(mpi_nreq), mpi_err)
        end if
    end do
!   data exchange.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   data exchange on myid.
    do isr=mesh(lev)%p2p_1,mesh(lev)%p2p_0
        ip_remote   =  p2p(isr)%ip_remote
        if(myid .ne. ip_remote) cycle
        if(p2p(isr)%n_send .ne. p2p(isr)%n_recv)    stop 'Error: s&r not match on myid.'

        if(p2p(isr)%n_send .gt. 0)  call DCOPY(p2p(isr)%n_send, p2p(isr)%rsend, 1, &
            &  p2p(isr)%rrecv, 1)
    end do
!   data exchange on myid.
!   ----------------------------------------------------------------------------

    call mpi_waitall(mpi_nreq, mpi_req, mpi_sta, mpi_err)

!   ----------------------------------------------------------------------------
!   unpack data received.
    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(.not. sec(isec)%is_ghost)    cycle
        do i=1,sec(isec)%n_ele
            isr =  sec(isec)%id_recv(1,i)
            iele=  sec(isec)%id_recv(2,i)
            sec(isec)%vol(i)=  p2p(isr)%rrecv(iele)
        end do
    end do
!   unpack data received.
!   ----------------------------------------------------------------------------

    return
    end subroutine reset_vol_hanging_node_parallel
