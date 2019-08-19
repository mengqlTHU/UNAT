!-------------------------------------------------------------------------------
!   get the structured stecil.
!-------------------------------------------------------------------------------
    subroutine get_structured_stencil
    use var_kind_def
    use var_cgns
    use var_global, only: is_2d_cal,is_structured_stencil,err_mem
    use var_load_balance
    use var_mesh
    use var_parallel
    use var_per_bnd, p=>per
    implicit none
    logical(dpL):: ltmp
    integer(dpI):: nfac,ifac,i,j,k,M,N,isec,iele,ele,v(8),w(4),LL,RR,L,R,npe,ip, &
                &  f(4,6),v_L1(4),v_L2(4),v_R1(4),v_R2(4),n_fuse,LDA
    integer(dpI),allocatable:: fac(:,:),fac_b(:,:),per(:,:)

    if(.not. is_structured_stencil) return

!   ----------------------------------------------------------------------------
!   record the vertex and element of face.
    nfac=  0
    do isec=mesh(0)%sec_1,mesh(0)%sec_0
        nfac=  nfac+sec(isec)%n_ele*nface_ele(sec(isec)%ele_type)
    end do
    allocate(fac(7,nfac), stat=err_mem)
    fac =  0

    nfac=  0
    k   =  huge(1)-1
    do isec=mesh(0)%sec_1,mesh(0)%sec_0
    do iele=1,sec(isec)%n_ele
        ele =  sec(isec)%ele_type

        if((ele .eq. TRI_3)) then
            v(1:3)  =  sec(isec)%n2e(1:3,iele)
            fac(1:5,nfac+1) = (/v(1), v(2), k, k, sec(isec)%ID_ele_i(iele)/)
            fac(1:5,nfac+2) = (/v(2), v(3), k, k, sec(isec)%ID_ele_i(iele)/)
            fac(1:5,nfac+3) = (/v(3), v(1), k, k, sec(isec)%ID_ele_i(iele)/)
            nfac=  nfac+3
        elseif((ele .eq. QUAD_4) .or. (ele .eq. QUAD_8) .or. (ele .eq. QUAD_9)) then
            v(1:4)  =  sec(isec)%n2e(1:4,iele)
            fac(1:5,nfac+1) = (/v(1), v(2), k, k, sec(isec)%ID_ele_i(iele)/)
            fac(1:5,nfac+2) = (/v(2), v(3), k, k, sec(isec)%ID_ele_i(iele)/)
            fac(1:5,nfac+3) = (/v(3), v(4), k, k, sec(isec)%ID_ele_i(iele)/)
            fac(1:5,nfac+4) = (/v(4), v(1), k, k, sec(isec)%ID_ele_i(iele)/)
            nfac=  nfac+4
        elseif(ele .eq. TETRA_4) then
            v(1:4)  =  sec(isec)%n2e(1:4,iele)
            fac(1:5,nfac+1) = (/v(1), v(2), v(3), k, sec(isec)%ID_ele_i(iele)/)
            fac(1:5,nfac+2) = (/v(1), v(2), v(4), k, sec(isec)%ID_ele_i(iele)/)
            fac(1:5,nfac+3) = (/v(2), v(3), v(4), k, sec(isec)%ID_ele_i(iele)/)
            fac(1:5,nfac+4) = (/v(1), v(3), v(4), k, sec(isec)%ID_ele_i(iele)/)
            nfac=  nfac+4
        elseif(ele .eq. PYRA_5) then
            v(1:5)  =  sec(isec)%n2e(1:5,iele)
            fac(1:5,nfac+1) = (/v(1), v(2), v(3), v(4), sec(isec)%ID_ele_i(iele)/)
            fac(1:5,nfac+2) = (/v(1), v(2), v(5), k   , sec(isec)%ID_ele_i(iele)/)
            fac(1:5,nfac+3) = (/v(2), v(3), v(5), k   , sec(isec)%ID_ele_i(iele)/)
            fac(1:5,nfac+4) = (/v(3), v(4), v(5), k   , sec(isec)%ID_ele_i(iele)/)
            fac(1:5,nfac+5) = (/v(4), v(1), v(5), k   , sec(isec)%ID_ele_i(iele)/)
            nfac=  nfac+5
        elseif(ele .eq. PENTA_6) then
            v(1:6)  =  sec(isec)%n2e(1:6,iele)
            fac(1:5,nfac+1) = (/v(1), v(2), v(3), k   , sec(isec)%ID_ele_i(iele)/)
            fac(1:5,nfac+2) = (/v(4), v(5), v(6), k   , sec(isec)%ID_ele_i(iele)/)
            fac(1:5,nfac+3) = (/v(1), v(2), v(5), v(4), sec(isec)%ID_ele_i(iele)/)
            fac(1:5,nfac+4) = (/v(2), v(3), v(6), v(5), sec(isec)%ID_ele_i(iele)/)
            fac(1:5,nfac+5) = (/v(1), v(3), v(6), v(4), sec(isec)%ID_ele_i(iele)/)
            nfac=  nfac+5
        elseif((ele .eq. HEXA_8) .or. (ele .eq. HEXA_20) .or. (ele .eq. HEXA_27)) then
            v(1:8)  =  sec(isec)%n2e(1:8,iele)
            fac(1:5,nfac+1) = (/v(1), v(2), v(3), v(4), sec(isec)%ID_ele_i(iele)/)
            fac(1:5,nfac+2) = (/v(5), v(6), v(7), v(8), sec(isec)%ID_ele_i(iele)/)
            fac(1:5,nfac+3) = (/v(1), v(5), v(8), v(4), sec(isec)%ID_ele_i(iele)/)
            fac(1:5,nfac+4) = (/v(2), v(6), v(7), v(3), sec(isec)%ID_ele_i(iele)/)
            fac(1:5,nfac+5) = (/v(1), v(2), v(6), v(5), sec(isec)%ID_ele_i(iele)/)
            fac(1:5,nfac+6) = (/v(4), v(3), v(7), v(8), sec(isec)%ID_ele_i(iele)/)
            nfac=  nfac+6
        else
            stop 'Error: element type not supported.'
        end if
    end do
    end do

    do i=1,nfac
        call iqsort(.true., 1, 4, fac(1,i))
    end do
    call iqsortcols(.true., 1, nfac, 1, 7, fac)
!   record the vertex and element of face.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   delete the faces covered by physical boundary.
    allocate(fac_b(4,n_eleb_g), stat=err_mem)
    fac_b   =  huge(1)-1
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
            fac(1:7,nfac)   =  fac(1:7,iele)
        end if
    end do
    if(allocated(fac_b))    deallocate(fac_b)
!   delete the faces covered by physical boundary.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   find L&R elements of face on myid.
    do i=1,nfac
        if(fac(6,i) .ne. 0) cycle
        v(1:4)  =  fac(1:4,i)

        call ib_search(1, nfac, 1, 7, fac, v(1), ltmp, L, R)
        if(.not. ltmp)  cycle
        do j=L,R
            if((fac(6,j) .gt. 0) .or. (i .eq. j))   cycle
            if((fac(2,j) .eq. v(2)) .and. (fac(3,j) .eq. v(3)) .and. &
              &(fac(4,j) .eq. v(4))) then
                LL  =  min(fac(5,i), fac(5,j))
                RR  =  max(fac(5,i), fac(5,j))
                if(fac(5,i) .eq. LL) then
                    fac(6,i)=  RR
                else
                    fac(6,i)=  LL
                end if
                fac(2:6,j)  = -1
                exit
            end if
        end do
    end do
    j   =  nfac
    nfac=  0
    do i=1,j
        if(fac(2,i) .lt. 0) cycle
        nfac            =  nfac+1
        fac(1:7,nfac)   =  fac(1:7,i)
    end do
!   find L&R elements of face on myid.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   find L&R elements of face using n2e of ghost element.
    k   =  huge(1)-1
    do i=1,n_ele_ghost
        ele =  ele_type_ghost(i)
        npe =  npe_ele_1st(ele)
        v(1:npe)=  jA_ghost(iA_ghost(i):iA_ghost(i)+npe-1)

        if((ele .eq. TRI_3)) then
            f(1:4,1)= (/v(1), v(2), k, k/)
            f(1:4,2)= (/v(2), v(3), k, k/)
            f(1:4,3)= (/v(3), v(1), k, k/)
        elseif((ele .eq. QUAD_4) .or. (ele .eq. QUAD_8) .or. (ele .eq. QUAD_9)) then
            f(1:4,1)= (/v(1), v(2), k, k/)
            f(1:4,2)= (/v(2), v(3), k, k/)
            f(1:4,3)= (/v(3), v(4), k, k/)
            f(1:4,4)= (/v(4), v(1), k, k/)
        elseif(ele .eq. TETRA_4) then
            f(1:4,1)= (/v(1), v(2), v(3), k/)
            f(1:4,2)= (/v(1), v(2), v(4), k/)
            f(1:4,3)= (/v(2), v(3), v(4), k/)
            f(1:4,4)= (/v(1), v(3), v(4), k/)
        elseif(ele .eq. PYRA_5) then
            f(1:4,1)= (/v(1), v(2), v(3), v(4)/)
            f(1:4,2)= (/v(1), v(2), v(5), k   /)
            f(1:4,3)= (/v(2), v(3), v(5), k   /)
            f(1:4,4)= (/v(3), v(4), v(5), k   /)
            f(1:4,5)= (/v(4), v(1), v(5), k   /)
        elseif(ele .eq. PENTA_6) then
            f(1:4,1)= (/v(1), v(2), v(3), k   /)
            f(1:4,2)= (/v(4), v(5), v(6), k   /)
            f(1:4,3)= (/v(1), v(2), v(5), v(4)/)
            f(1:4,4)= (/v(2), v(3), v(6), v(5)/)
            f(1:4,5)= (/v(1), v(3), v(6), v(4)/)
        elseif((ele .eq. HEXA_8) .or. (ele .eq. HEXA_20) .or. (ele .eq. HEXA_27)) then
            f(1:4,1)= (/v(1), v(2), v(3), v(4)/)
            f(1:4,2)= (/v(5), v(6), v(7), v(8)/)
            f(1:4,3)= (/v(1), v(5), v(8), v(4)/)
            f(1:4,4)= (/v(2), v(6), v(7), v(3)/)
            f(1:4,5)= (/v(1), v(2), v(6), v(5)/)
            f(1:4,6)= (/v(4), v(3), v(7), v(8)/)
        else
            stop 'Error: element type not supported.'
        end if

        do j=1,nface_ele(ele)
            call find_col_in_matrix(7, 1, 4, 1, nfac, fac, f(1,j), ip, w)
            if(ip .le. 0) then
!               do nothing.
            elseif(ip .eq. 1) then
                fac(6,w(1)) =  ele_idx_ghost(i)
            else
                stop 'Error: find 2 faces between element&ghost_element.'
            end if
        end do
    end do
!   find L&R elements of face using n2e of ghost element.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   find L&R elements of face using n2e of periodic ghost element.
    allocate(per(4,nfac), stat=err_mem)
    if(err_mem .ne. 0)  stop 'Error: fails to allocate memory, get_structured_stencil.'
    k   =  huge(1)-1
    per =  k

    if(n_vtx_pair_per .gt. 0) then
    do i=1,n_ele_per
        if(is_corner_per(i) .or. (ele_per(3,i) .gt. 0)) cycle
        ele =  ele_type_per(i)
        npe =  npe_ele_1st(ele)
        v(1:npe)=  jA_per(iA_per(i):iA_per(i)+npe-1)

        if((ele .eq. TRI_3)) then
            f(1:4,1)= (/v(1), v(2), k, k/)
            f(1:4,2)= (/v(2), v(3), k, k/)
            f(1:4,3)= (/v(3), v(1), k, k/)
        elseif((ele .eq. QUAD_4) .or. (ele .eq. QUAD_8) .or. (ele .eq. QUAD_9)) then
            f(1:4,1)= (/v(1), v(2), k, k/)
            f(1:4,2)= (/v(2), v(3), k, k/)
            f(1:4,3)= (/v(3), v(4), k, k/)
            f(1:4,4)= (/v(4), v(1), k, k/)
        elseif(ele .eq. TETRA_4) then
            f(1:4,1)= (/v(1), v(2), v(3), k/)
            f(1:4,2)= (/v(1), v(2), v(4), k/)
            f(1:4,3)= (/v(2), v(3), v(4), k/)
            f(1:4,4)= (/v(1), v(3), v(4), k/)
        elseif(ele .eq. PYRA_5) then
            f(1:4,1)= (/v(1), v(2), v(3), v(4)/)
            f(1:4,2)= (/v(1), v(2), v(5), k   /)
            f(1:4,3)= (/v(2), v(3), v(5), k   /)
            f(1:4,4)= (/v(3), v(4), v(5), k   /)
            f(1:4,5)= (/v(4), v(1), v(5), k   /)
        elseif(ele .eq. PENTA_6) then
            f(1:4,1)= (/v(1), v(2), v(3), k   /)
            f(1:4,2)= (/v(4), v(5), v(6), k   /)
            f(1:4,3)= (/v(1), v(2), v(5), v(4)/)
            f(1:4,4)= (/v(2), v(3), v(6), v(5)/)
            f(1:4,5)= (/v(1), v(3), v(6), v(4)/)
        elseif((ele .eq. HEXA_8) .or. (ele .eq. HEXA_20) .or. (ele .eq. HEXA_27)) then
            f(1:4,1)= (/v(1), v(2), v(3), v(4)/)
            f(1:4,2)= (/v(5), v(6), v(7), v(8)/)
            f(1:4,3)= (/v(1), v(5), v(8), v(4)/)
            f(1:4,4)= (/v(2), v(6), v(7), v(3)/)
            f(1:4,5)= (/v(1), v(2), v(6), v(5)/)
            f(1:4,6)= (/v(4), v(3), v(7), v(8)/)
        else
            stop 'Error: element type not supported.'
        end if
        f(1:4,1:nface_ele(ele)) =  abs(f(1:4,1:nface_ele(ele)))

        do j=1,nface_ele(ele)
            npe =  0
            do M=1,4
                if(f(M,j) .eq. k) then
                    exit
                else
                    npe =  npe+1
                end if
            end do

            LL  =  0
            do M=1,npe
                call ib_search(1, n_vtx_pair_per, 1, 5, vtx_pair_per, f(M,j), ltmp, L, R)
                if(.not. ltmp)  cycle
                do ip=L,R
                    if((ele_per(2,i) .eq. vtx_pair_per(3,ip)) .and. &
                      &(ele_per(3,i) .eq. vtx_pair_per(4,ip)) .and. &
                      &(ele_per(4,i) .eq. vtx_pair_per(5,ip))) then
                        LL  =  LL+1
                        v(M)=  vtx_pair_per(2,ip)
                        exit
                    end if
                end do
            end do
            if(LL .ne. npe) cycle

            call find_col_in_matrix(7, 1, npe, 1, nfac, fac, v, ip, w)
            if(ip .le. 0) then
!               do nothing.
            elseif(ip .eq. 1) then
                fac(6:7  ,w(1)) =  ele_per(1:2,i)
                per(1:npe,w(1)) =  f(1:npe,j)
            else
                stop 'Error: find 2 faces between element&ghost_element.'
            end if
        end do
    end do
    end if
!   find L&R elements of face using n2e of periodic ghost element.
!   ----------------------------------------------------------------------------

    j   =  nfac
    nfac=  0
    do i=1,j
        if(fac(6,i) .le. 0) cycle
        nfac            =  nfac+1
        fac(1:7,nfac)   =  fac(1:7,i)
        per(1:4,nfac)   =  per(1:4,i)
    end do

!   ----------------------------------------------------------------------------
!   record the stencil for which RR can not be found in the internal elements.
    n_fuse      =  0
    n2e_i_have  =  0
    LDA         =  9
    do ifac=1,nfac
        L   =  fac(5,ifac)
        R   =  fac(6,ifac)

        call get_isec_iele(R, isec, iele)
!       R is internal element of myid, then no need to worry about the structured stencil.
        if((isec .gt. 0) .and. (fac(7,ifac) .le. 0)) cycle

!       check the element type, L side.
        do isec=1,n_seci_g
            if((L .ge. seci_g(6,isec)) .and. (L .le. seci_g(7,isec)))   exit
        end do
        ele =  seci_g(1,isec)
        ltmp= (ele .eq. QUAD_4) .or. ((ele .eq. PENTA_6) .and. (npe .eq. 3)) .or. &
            & (ele .eq. HEXA_8)
        if(.not. ltmp)  cycle

!       check the element type, R side.
        do isec=1,n_seci_g
            if((R .ge. seci_g(6,isec)) .and. (R .le. seci_g(7,isec)))   exit
        end do
        ele =  seci_g(1,isec)
        ltmp= (ele .eq. QUAD_4) .or. ((ele .eq. PENTA_6) .and. (npe .eq. 3)) .or. &
            & (ele .eq. HEXA_8)
        if(.not. ltmp)  cycle

!       how many vertex this face has?
        npe =  0
        do k=1,4
            if(fac(k,ifac) .lt. huge(1)-1)  npe =  npe+1
        end do

        v_L1(1:npe) =  fac(1:npe,ifac)
        if(fac(7,ifac) .le. 0) then
            v_R1(1:npe) =  v_L1(1:npe)
        else
            v_R1(1:npe) =  per(1:npe,ifac)
        end if

        call get_vtx_opposite(L, npe, v_L1, j, v_L2)
        call get_ID_face(L, j, v_L2, LL)
        if(LL .le. 0)   cycle

        call get_vtx_opposite(R, npe, v_R1, j, v_R2)
        n_fuse  =  n_fuse+1
        if(LDA*n_fuse .ge. size(n2e_i_have))    stop 'Error: size(n2e_i_have) too small.'
!       n2e_R2, R, R_per, RR, RR_per, myid.
        n2e_i_have(LDA*(n_fuse-1)+1:LDA*(n_fuse-1)+j)   =  abs(v_R2(1:j))
        n2e_i_have(LDA*(n_fuse-1)+5                 )   =  R
        n2e_i_have(LDA*(n_fuse-1)+6                 )   =  fac(7,ifac)
!       n2e_i_have(LDA*(n_fuse-1)+7                 )   =  RR
!       n2e_i_have(LDA*(n_fuse-1)+8                 )   =  per
        n2e_i_have(LDA*(n_fuse-1)+9                 )   =  myid
    end do
!   record the stencil for which RR can not be found in the internal elements.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   get the ip of R.
    call iqsortcols(.true., 1, n_fuse, 5, 9, n2e_i_have)
    prc_info=  0
    do i=1,n_fuse
        iele=  n2e_i_have(LDA*(i-1)+5)
        if(iele .ge. elei_prc(1,nprc-1)) then
            ip  =  nprc-1
        else
            ip  = (iele-1)/(n_elei_g/nprc)
        end if
        prc_info(ip,1)  =  prc_info(ip,1)+1
        ele_i_need(i)   =  iele
    end do
    call mpi_alltoall(prc_info(0,1), 1, mpi_dpI, prc_info(0,3), 1, mpi_dpI, &
        &  mpi_comm_world, mpi_err)
    n   =  prc_info(0,3)
    do ip=1,nprc-1
        n   =  n+prc_info(ip,3)
    end do
    if(allocated(ele_i_have)) then
        if(n .gt. size(ele_i_have)) then
            deallocate(ele_i_have)
            allocate(ele_i_have(n), stat=err_mem)
        end if
    else
        allocate(ele_i_have(max(n,1)), stat=err_mem)
    end if
    call collect_data_parallel(prc_info(0,1), prc_info(0,2), ele_i_need, prc_info(0,3), &
        &  prc_info(0,4), ele_i_have)
    do i=1,prc_info(nprc-1,3)+prc_info(nprc-1,4)-1
        iele=  ele_i_have(i)
        if((iele .lt. elei_prc(1,myid)) .or. (iele .gt. elei_prc(2,myid))) &
            &  stop 'Error: fails to get ip of element.'
        ele_i_have(i)   =  part(iele-elei_prc(1,myid)+1)
    end do
    prc_info(0:nprc-1,1)=  prc_info(0:nprc-1,3)
    call collect_data_parallel(prc_info(0,1), prc_info(0,2), ele_i_have, prc_info(0,3), &
        &  prc_info(0,4), ele_i_need)
    do i=1,n_fuse
        n2e_i_have(LDA*(i-1)+7) =  ele_i_need(i)
    end do
!   get the ip of R.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   send the above information to the host and find RR.
    call iqsortcols(.true., 1, n_fuse, 7, 9, n2e_i_have)
    prc_info=  0
    do i=1,n_fuse
        ip              =  n2e_i_have(LDA*(i-1)+7)
        prc_info(ip,1)  =  prc_info(ip,1)+LDA
    end do
    call mpi_alltoall(prc_info(0,1), 1, mpi_dpI, prc_info(0,3), 1, mpi_dpI, &
        &  mpi_comm_world, mpi_err)
    n   =  prc_info(0,3)
    do ip=1,nprc-1
        n   =  n+prc_info(ip,3)
    end do
    if(allocated(n2e_i_need)) then
        if(n .gt. size(n2e_i_need)) then
            deallocate(n2e_i_need)
            allocate(n2e_i_need(n), stat=err_mem)
        end if
    else
        allocate(n2e_i_need(max(n,1)), stat=err_mem)
    end if
    call collect_data_parallel(prc_info(0,1), prc_info(0,2), n2e_i_have, prc_info(0,3), &
        &  prc_info(0,4), n2e_i_need)
    do i=1,(prc_info(nprc-1,3)+prc_info(nprc-1,4)-1)/LDA
!       n2e_R2, R, R_per, RR, RR_per, myid.
        n2e_i_need(LDA*(i-1)+7) =  0
        R   =  n2e_i_need(LDA*(i-1)+5)
        call get_isec_iele(R, isec, iele)
        if(isec .le. 0) stop 'Error: fails to find R in the internal elements.'

        v_R2(1:4)   =  n2e_i_need(LDA*(i-1)+1:LDA*(i-1)+4)
        npe =  0
        do j=1,4
            if(v_R2(j) .gt. 0) then
                npe =  npe+1
            else
                exit
            end if
        end do
        call find_col_in_matrix(7, 1, npe, 1, nfac, fac, v_R2, k, w)
        if(k .le. 0) then
            cycle
        elseif(k .ge. 2) then
            stop 'Error: find 2 faces between elements.'
        end if
        if(fac(6,w(1)) .le. 0)  cycle

        if(fac(7,w(1)) .gt. 0) then
            if(fac(5,w(1)) .eq. R) then
                n2e_i_need(LDA*(i-1)+7) =  fac(6,w(1))
                n2e_i_need(LDA*(i-1)+8) =  fac(7,w(1))
            end if
        else
            if(fac(5,w(1)) .eq. R) then
                n2e_i_need(LDA*(i-1)+7) =  fac(6,w(1))
            elseif(fac(6,w(1)) .eq. R) then
                n2e_i_need(LDA*(i-1)+7) =  fac(5,w(1))
            end if
        end if
    end do
    prc_info(0:nprc-1,1)=  prc_info(0:nprc-1,3)
    call collect_data_parallel(prc_info(0,1), prc_info(0,2), n2e_i_need, prc_info(0,3), &
        &  prc_info(0,4), n2e_i_have)
!   send the above information to the host and find RR.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   record RR and delete duplicate information.
    n_ss=  0
    if(allocated(per))  per =  0
    do i=1,(prc_info(nprc-1,3)+prc_info(nprc-1,4)-1)/LDA
!       n2e_R2, R, R_per, RR, RR_per, myid.
        if(n2e_i_have(LDA*(i-1)+9) .ne. myid)   stop 'Error: get_structured_stencil fail.'
        if(n2e_i_have(LDA*(i-1)+7) .le. 0   )   cycle
        n_ss=  n_ss+1
        per(1,n_ss) =  n2e_i_have(LDA*(i-1)+7)
        j           =  n2e_i_have(LDA*(i-1)+6)
        k           =  n2e_i_have(LDA*(i-1)+8)
        if((j .le. 0) .and. (k .le. 0)) then
!           do nothing.
        elseif((j .gt. 0) .and. (k .gt. 0)) then
            per(2,n_ss) =  k
            per(3,n_ss) =  j
        else
            per(2,n_ss) =  max(j,k)
        end if
    end do
    if(n_ss .gt. 0) call iqsortcols(.true., 1, n_ss, 1, 4, per)
    i   =  1
    do while(i .le. n_ss)
        k   =  i
        do j=i+1,n_ss
            if(per(1,j) .eq. per(1,i)) then
                k   =  j
            else
                exit
            end if
        end do
        do j=i+1,k
        do M=i,j-1
            if(per(1,M) .le. 0) cycle
            if((per(2,j) .eq. per(2,M)) .and. (per(3,j) .eq. per(3,M)) .and. &
              &(per(4,j) .eq. per(4,M))) then
                per(1:4,j)  = -1
                exit
            end if
        end do
        end do

        i   =  k+1
    end do

    do i=1,n_ss
        if((per(1,i) .le. 0) .or. (per(2,i) .gt. 0))    cycle
        call get_isec_iele(per(1,i), isec, iele)
!       this structured stencil belongs to myid.
        if(isec .gt. 0) per(1:4,i)  = -1
    end do

    j   =  n_ss
    n_ss=  0
    do i=1,j
        if(per(1,i) .le. 0) cycle
        n_ss            =  n_ss+1
        per(1:4,n_ss)   =  per(1:4,i)
    end do
!   record RR and delete duplicate information.
!   ----------------------------------------------------------------------------

    if(allocated(fac_b))    deallocate(fac_b)
    if(allocated(fac  ))    deallocate(fac  )

!   ----------------------------------------------------------------------------
!   get n2e of RR.
    if(n_ss .gt. 0) then
        call iqsortcols(.true., 1, n_ss, 1, 4, n2e_i_have)
        allocate(ele_ID_ss  (n_ss), stat=err_mem)
        allocate(ele_type_ss(n_ss), stat=err_mem)
        allocate(per_ss   (3,n_ss), stat=err_mem)
        allocate(iA_ss    (n_ss+1), stat=err_mem)
        per_ss  =  0
        iA_ss(1)=  1
    end if
    if(n_ss .gt. 0) call ICOPY(n_ss, per, 4, ele_i_need, 1)
    n_ele_i_need=  n_ss
    call get_n2e_i_need
    i   =  0
    M   =  0
    do while(i .lt. nnz_n2e_i_need)
!       ID, ele_type, n_n2e, n2e.
        i   =  i+1
        RR  =  n2e_i_need(i)
        i   =  i+1
        ele =  n2e_i_need(i)
        i   =  i+1
        npe =  n2e_i_need(i)
        i   =  i+npe

        M   =  M+1
        if(RR .ne. per(1,M))    stop 'Error: fails to collect n2e of RR.'
        iA_ss(M+1)  =  iA_ss(M)+npe_ele(ele)
    end do
    if(M .ne. n_ss) stop 'Error: fails to collect n2e of RR.'
    if(n_ss .gt. 0) allocate(jA_ss(iA_ss(n_ss+1)-1), stat=err_mem)
    i   =  0
    M   =  0
    do while(i .lt. nnz_n2e_i_need)
!       ID, ele_type, n_n2e, n2e.
        i   =  i+1
        RR  =  n2e_i_need(i)
        i   =  i+1
        ele =  n2e_i_need(i)
        i   =  i+1
        npe =  n2e_i_need(i)

        M   =  M+1
        ele_ID_ss  (M)  =  RR
        ele_type_ss(M)  =  ele
        per_ss (1:2,M)  =  per(2:3,M)
        do j=1,npe_ele(ele)
            jA_ss(iA_ss(M)+j-1) =  n2e_i_need(i+j)
        end do
        i   =  i+npe
    end do
!   get n2e of RR.
!   ----------------------------------------------------------------------------

    if(allocated(per))  deallocate(per)

    return
    contains
!       ------------------------------------------------------------------------
!       get element ID of the other side of a face.
!       ------------------------------------------------------------------------
        subroutine get_ID_face(L,npe,f,R)
        use var_kind_def
        implicit none
        integer(dpI),intent(in):: L,npe,f(*)
        integer(dpI):: R,i,v(10)

        R   = -1
        call find_col_in_matrix(7, 1, npe, 1, nfac, fac, f, i, v)
        if(i .ne. 1)    return
        if(fac(5,v(1)) .eq. L) then
            R   =  fac(6,v(1))
        elseif(fac(6,v(1)) .eq. L) then
            R   =  fac(5,v(1))
        else
            stop 'Error: get_ID_face fails.'
        end if

        return
        end subroutine get_ID_face
!       ------------------------------------------------------------------------
!       get vertex of the opposite face.
!       ------------------------------------------------------------------------
        subroutine get_vtx_opposite(ID,n1,f1,n2,f2)
        use var_kind_def
        use var_cgns, only: npe_ele_1st
        use var_load_balance
        use var_mesh, only: sec
        implicit none
        integer(dpI),intent(in):: ID,n1,f1(*)
        logical(dpL):: ltmp
        integer(dpI):: n2,f2(*),ff(8),L,R,isec,iele

        ff(1)   = -1
        do while(.true.)
            call get_isec_iele(ID, isec, iele)
            if(isec .gt. 0) then
                n2      =  npe_ele_1st(sec(isec)%ele_type)
                ff(1:n2)=  sec(isec)%n2e(1:n2,iele)
                exit
            end if

            call ib_search(1, n_ele_ghost, 1, 1, ele_idx_ghost, ID, ltmp, L, R)
            if(ltmp) then
                n2      =  npe_ele_1st(ele_type_ghost(L))
                ff(1:n2)=  jA_ghost(iA_ghost(L):iA_ghost(L)+n2-1)
                exit
            end if

            call ib_search(1, n_ele_per, 1, 4, ele_per, ID, ltmp, L, R)
            if(ltmp) then
                n2      =  npe_ele_1st(ele_type_per(L))
                ff(1:n2)=  abs(jA_per(iA_per(L):iA_per(L)+n2-1))
                exit
            end if
            stop 'Error: fails to find the element.'
        end do
        call list_minus_list(n2, ff, n1, f1)
        f2(1:n2)=  ff(1:n2)

        return
        end subroutine get_vtx_opposite
    end subroutine get_structured_stencil
