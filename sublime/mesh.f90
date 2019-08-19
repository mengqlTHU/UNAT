!-------------------------------------------------------------------------------
!   definition of mesh.
!-------------------------------------------------------------------------------
    module var_mesh
        use var_kind_def
        implicit none

!       ------------------------------------------------------------------------
!       information of the original mesh file.
        character(len=80),allocatable:: sec_name_g(:),sec_name(:),bsec_name(:)
        integer(dpI):: n_eleb_g =  0
        integer(dpI):: n_sec_g  =  0
!       ele_type, npe, ele1, ele0, parent_flag
        integer(dpI),allocatable:: sec_g(:,:)
!       ele_type, npe, eleg_1, eleg_0, isec_g, elei_1, elei_0
        integer(dpI),allocatable:: seci_g(:,:)

        integer(dpI):: n_LR_hanging =  0
        integer(dpI),allocatable:: LR_hanging(:,:)

        integer(dpI):: n_ss =  0
!       information of the original mesh file.
!       ------------------------------------------------------------------------

        integer(dpI),parameter:: max_sec_lev=  50
        type type_section
            integer(dpI):: ID_sec_g =  0
            integer(dpI):: ID_sec_i =  0
            integer(dpI):: ele_type =  0
            integer(dpI):: npe      =  0
            integer(dpI):: n_ele    =  0

            integer(dpI):: n_ele_b  =  0
            integer(dpI),allocatable:: ID_ele_i(:),ID_ele_b(:)
            ! SW_slave
            integer(dpI),allocatable:: ID_ele_g(:)

!           For ghost element only, ip2p and ele_p2p.
            integer(dpI),allocatable:: ID_recv (:,:)
            integer(dpI),allocatable:: per_path(:,:)

            integer(dpI),allocatable:: ID_ghost_ele(:,:)

!           For mesh coarsen only.
            integer(dpI),allocatable:: ID_coarse_mesh(:,:)

            integer(dpI),allocatable:: n2e(:,:)

            logical(dpL):: is_int   =  .true.
            logical(dpL):: is_ghost =  .false.
            logical(dpL):: is_bnd   =  .false.
            logical(dpL):: is_cht   =  .false.
            logical(dpL):: is_t     =  .false.
            logical(dpL):: is_flux  =  .false.
            integer(dpI):: ID_bnd   =  0
            integer(dpI):: ID_group =  0
            integer(dpI):: bct      =  0

            logical(dpL):: is_bar   =  .false.
            logical(dpL):: is_tri   =  .false.
            logical(dpL):: is_quad  =  .false.
            logical(dpL):: is_hexa  =  .false.

            integer(dpI):: lev  =  0

            integer(dpI),allocatable:: iA_face_neighbour(:),jA_face_neighbour(:)
            integer(dpI),allocatable:: iA_ls(:),jA_ls(:,:)
            real   (dpR),allocatable:: cen(:,:),vol(:),gls(:),dnw(:),gls_wn(:),wn(:,:)
            real   (dpR),allocatable:: coe_ls(:,:),vg(:,:)

!           for boundary section ONLY.
            logical(dpL):: is_profiled  =  .false.
            integer(dpI):: n_vtx_profile=  0
            real   (dpR),allocatable:: bc_profile(:,:)
        end type type_section
        type(type_section),save:: sec(max_sec_lev*5)

        type type_mesh
            integer(dpI):: sec_1    =  1
            integer(dpI):: sec_0    =  0

            integer(dpI):: n_vtx    =  0
            integer(dpI),allocatable:: id_vtx(:)

            integer(dpI):: n_vtx_per=  0
            integer(dpI),allocatable:: id_vtx_per(:),per_path(:,:)

            real   (dpR),allocatable:: xyz(:,:)
            real   (dpR),allocatable:: J_xyz(:,:)

            integer(dpI):: p2p_1=  1
            integer(dpI):: p2p_0=  0

            integer(dpI):: n_mortar_b       =  0
            integer(dpI):: n_mortar_ss      =  0
            integer(dpI):: n_mortar         =  0
            integer(dpI):: n_mortar_hanging =  0
            integer(dpI),allocatable:: mortar_ele_type(:)
            integer(dpI),allocatable:: mortar_LR      (:,:)
            ! SW_slave
            integer(dpI),allocatable:: mortar_own_ID    (:)
            integer(dpI),allocatable:: mortar_nei_ID    (:)
            integer(dpI),allocatable:: owner            (:)
            integer(dpI),allocatable:: neighbor         (:)
            real   (dpR),allocatable:: mortar_transform (:)
            integer(dpI),allocatable:: mortar_order     (:,:)

            integer(dpI),allocatable:: mortar_n2e     (:,:)
            integer(dpI),allocatable:: mortar_structured_stencil(:,:)
            integer(dpI),allocatable:: mortar_hanging (:  )
            real   (dpR),allocatable:: mortar_cen     (:,:)
            real   (dpR),allocatable:: mortar_n_vg    (:,:)

            integer(dpI):: rcm_b1   =  huge(1)
            integer(dpI):: rcm_b0   =  0
            integer(dpI),allocatable:: rcm_order(:,:)

            integer(dpI),allocatable:: iA_e2n(:),jA_e2n(:)
            real   (dpR),allocatable:: w_e2n(:)
        end type type_mesh
        type(type_mesh),save:: mesh(0:4)
        type(type_mesh),save:: mesh_reordered(0:4)

        integer(dpI),parameter:: max_p2p=  50
        type type_p2p_sr
            integer(dpI):: lev      =  0
            integer(dpI):: ip_remote=  0

            integer(dpI):: n_ele_send       =  0
            integer(dpI):: n_ele_recv       =  0
            integer(dpI):: n_ele_send_face  =  0
            integer(dpI):: n_ele_recv_face  =  0
            integer(dpI):: n_send           =  0
            integer(dpI):: n_recv           =  0
!           local ID of the element to be sent, isec and iele.
            integer(dpI),allocatable:: id_ele_send(:,:)
!           global ID of the internal element to be received.
            integer(dpI),allocatable:: id_ele_recv(:),iA_n_sp(:),iA_n_fp(:)
            integer(dpI),allocatable:: isend(:),irecv(:)
            real   (dpR),allocatable:: rsend(:),rrecv(:)
            real   (dpR),allocatable:: adj_send(:),adj_recv(:)
        end type type_p2p_sr
        type(type_p2p_sr),save:: p2p(max_p2p*5)

        integer(dpI):: n_sec_per_info   =  0
        real   (dpR),allocatable:: sec_per_info(:,:)

        integer(dpI),allocatable:: sec_motion_type (:)
        real   (dpR),allocatable:: sec_motion_speed(:)

        ! SW_slave
        integer(dpI),allocatable:: perm(:)
        integer(dpI):: cellNum = 0
        integer(dpI):: faceNum = 0
        integer(dpI):: tot_ele = 0
    end module var_mesh

!-------------------------------------------------------------------------------
!   SW_slave: section array
!-------------------------------------------------------------------------------
    module var_sec_array
    use var_kind_def
    implicit none
    real   (dpR),allocatable:: sec_cen(:,:),sec_vol(:),sec_coe_ls(:,:)
    real   (dpR),allocatable:: sec_is_int(:),sec_is_ghost(:),sec_is_bnd(:)
    real   (dpR),allocatable:: sec_bct(:),sec_ID_ghost_ele(:)
        
    end module var_sec_array
!-------------------------------------------------------------------------------
!   module to define periodic boundary.
!-------------------------------------------------------------------------------
    module var_per_bnd
        use var_kind_def
        implicit none

        integer(dpI):: n_per        =  0
        integer(dpI):: n_per_info   =  0
        real   (dpR),allocatable:: per_info(:,:),per_mat(:,:)

        type type_per_bnd
            integer(dpI):: per_type =  1
            integer(dpI):: nvtx     =  0
            integer(dpI),allocatable:: vtx(:,:)

            real   (dpR):: translation   (3)=  0.0d0
            real   (dpR):: rotationcenter(3)=  0.0d0
            integer(dpI):: rotationangle (3)=  10
        end type type_per_bnd
        type(type_per_bnd),allocatable:: per(:)
        integer(dpI):: n_vtx_pair_per   =  0
        integer(dpI),allocatable:: vtx_pair_per(:,:)
    end module var_per_bnd
!-------------------------------------------------------------------------------
!   read mesh.
!-------------------------------------------------------------------------------
    subroutine read_mesh
    use var_kind_def
    use var_cgns
    use var_global
    use var_mesh
    use var_parallel
    implicit none
    integer(dpI):: isec

    if(mesh_file(len_trim(mesh_file)-3:len_trim(mesh_file)) .eq. 'cgns') then
        call read_mesh_cgns
    elseif(mesh_file(len_trim(mesh_file)-2:len_trim(mesh_file)) .eq. 'msh') then
!       call read_mesh_gmsh
    end if

    call get_vtx_pair_per
    call get_faces_hanging_node
    call get_e2n_parallel
    call get_ele_ghost(.true., .false.)
    call get_ele_per(.true., .true.)

    call load_balance

    call set_section
    call set_ghost_section
    call set_p2p
    call set_vertex
    call set_bnd_section
    call reorder_internal_element
    call load_balance_clean
!   call wr_local_mesh

    if(mesh(0)%sec_0 .gt. max_sec_lev)  stop 'Error: max_sec_lev too small.'
    do isec=mesh(0)%sec_1,mesh(0)%sec_0
        sec(isec)%is_bar    = (sec(isec)%ele_type .eq. BAR_2) .or. &
                            & (sec(isec)%ele_type .eq. BAR_3) .or. &
                            & (sec(isec)%ele_type .eq. BAR_5)
        sec(isec)%is_tri    = (sec(isec)%ele_type .eq. TRI_3)
        sec(isec)%is_quad   = (sec(isec)%ele_type .eq. QUAD_4) .or. &
                            & (sec(isec)%ele_type .eq. QUAD_8) .or. &
                            & (sec(isec)%ele_type .eq. QUAD_9) .or. &
                            & (sec(isec)%ele_type .eq. QUAD_25)
        sec(isec)%is_hexa   = (sec(isec)%ele_type .eq. HEXA_8 ) .or. &
                            & (sec(isec)%ele_type .eq. HEXA_20) .or. &
                            & (sec(isec)%ele_type .eq. HEXA_27)
        sec(isec)%lev       =  0
    end do
    call get_section_ID
    call get_sec_per_info
    call get_sec_motion_info

if(sw_slave) then
    call set_ID_ele_g(0)
end if

    call set_mortar
    call set_p2p_element_order(0)
    call set_monitor

if(sw_slave) then
    call set_ID_ghost_ele_sw(0)
end if

    return
    end subroutine read_mesh
!-------------------------------------------------------------------------------
!   rotate a vector.
!-------------------------------------------------------------------------------
    subroutine per_rot_vec(cnt,p,v,w)
    use var_kind_def
    use var_global, only: n_dim
    use var_per_bnd
    implicit none
    integer(dpL),intent(in):: cnt,p(*)
    real   (dpR),intent(in):: v(*)
    integer(dpI):: i,N,per_type
    real   (dpR):: w(*),t(9)

!   cnt=0 , for xyz.
!   cnt=1 , for u.
!   cnt=2 , for turb.
!   cnt=11, for gradient of u.
!   cnt=12, for gradient of turb.
!   cnt=21, for du.
!   cnt=31, gradient of (u, v, w, t) in 2D.

    if(cnt .eq. 0) then
        w(1:n_dim)  =  v(1:n_dim)
    elseif(cnt .eq. 1) then
        w(1:5)  =  v(1:5)
    elseif(cnt .eq. 2) then
        w(1)    =  v(1)
        return
    elseif(cnt .eq. 11) then
        w(1:19) =  v(1:19)
    elseif(cnt .eq. 12) then
        w(1:3)  =  v(1:3)
    elseif(cnt .eq. 21) then
        w(1:5)  =  v(1:5)
    elseif(cnt .eq. 31) then
        w(1:8)  =  v(1:8)
        return
    else
        stop 'Error: wrong parameter for per_rot_vec.'
    end if

    do i=1,3
        N   =  p(i)
        if((N .le. 0) .or. (N .gt. n_per_info)) return
        per_type=  nint(per_info(1,N), dpI)

        if(per_type .eq. 1) then
            if(cnt .eq. 0)  w(1:n_dim)  =  w(1:n_dim)+per_info(8:7+n_dim,N)
        elseif(per_type .eq. 2) then
            if(cnt .eq. 0) then
                w(1:3)  =  per_mat(1:3,N)*w(1)+per_mat(4:6,N)*w(2)+per_mat(7:9,N)*w(3)
            elseif(cnt .eq. 1) then
                w(2:4)  =  per_mat(1:3,N)*w(2)+per_mat(4:6,N)*w(3)+per_mat(7:9,N)*w(4)
            elseif(cnt .eq. 2) then
                return
            elseif(cnt .eq. 11) then
                w(1 :3 )=  per_mat(1:3,N)*w(1 )+per_mat(4:6,N)*w(2 )+per_mat(7:9,N)*w(3 )
                call dgemm('N','N',3,3,3,1.0d0,per_mat(1,N),3,w(4),3,0.0d0,t,3)
                call dgemm('N','T',3,3,3,1.0d0,t,3,per_mat(1,N),3,0.0d0,w(4),3)
                w(13:15)=  per_mat(1:3,N)*w(13)+per_mat(4:6,N)*w(14)+per_mat(7:9,N)*w(15)
                w(16:18)=  per_mat(1:3,N)*w(16)+per_mat(4:6,N)*w(17)+per_mat(7:9,N)*w(18)
            elseif(cnt .eq. 12) then
                w(1 :3 )=  per_mat(1:3,N)*w(1 )+per_mat(4:6,N)*w(2 )+per_mat(7:9,N)*w(3 )
            elseif(cnt .eq. 21) then
                w(2:4)  =  per_mat(1:3,N)*w(2)+per_mat(4:6,N)*w(3)+per_mat(7:9,N)*w(4)
            else
                stop 'Error: wrong parameter for per_rot_vec.'
            end if
        end if
    end do

    return
    end subroutine per_rot_vec
!-------------------------------------------------------------------------------
!   get local isec and iele of ID_ele_i.
!-------------------------------------------------------------------------------
    subroutine get_isec_iele(ID,isec,iele)
    use var_kind_def
    use var_mesh
    implicit none
    integer(dpI),intent(in):: ID
    logical(dpL):: ltmp
    integer(dpI):: isec,iele,i,L,R

    isec= -1
    do i=mesh(0)%sec_1,mesh(0)%sec_0
        if(.not. sec(i)%is_int)  cycle
!       ID_ele_i should be sorted.
        call ib_search(1, sec(i)%n_ele_b, 1, 1, sec(i)%ID_ele_i, ID, ltmp, L, R)
        if(ltmp) then
            isec=  i
            iele=  L
            return
        end if

!       ID_ele_i should be sorted.
        call ib_search(1+sec(i)%n_ele_b,sec(i)%n_ele,1,1,sec(i)%ID_ele_i,ID,ltmp,L,R)
        if(ltmp) then
            isec=  i
            iele=  L
            return
        end if
    end do

    return
    end subroutine get_isec_iele
!-------------------------------------------------------------------------------
!   get the global section index.
!-------------------------------------------------------------------------------
    subroutine get_section_ID
    use var_kind_def
    use var_cgns
    use var_global, only: is_2d_cal
!   use var_load_balance
    use var_mesh
    implicit none
    logical(dpL):: ltmp
    integer(dpI):: isec,i,j,ele

    do isec=mesh(0)%sec_1,mesh(0)%sec_0
        sec(isec)%id_sec_g  =  0
        if(sec(isec)%is_bnd) then
            j   =  0
            do i=1,n_sec_g
                ele =  ele_type_1st(sec_g(1,i))
                if(is_2d_cal) then
                    ltmp=  ele .eq. BAR_2
                else
                    ltmp= (ele .eq. TRI_3) .or. (ele .eq. QUAD_4)
                end if
                if(ltmp)    j   =  j+1
                if(sec(isec)%id_bnd .eq. j) then
                    sec(isec)%id_sec_g  =  i
                    exit
                end if
            end do
        else
            if(sec(isec)%ID_sec_I .le. 0)   stop 'Error: wrong id of seci.'
            sec(isec)%id_sec_g  =  seci_g(5,sec(isec)%ID_sec_I)
        end if
        if(sec(isec)%id_sec_g .le. 0)   stop 'Error: fails to get id_sec_g.'
    end do

    return
    end subroutine get_section_ID
!-------------------------------------------------------------------------------
!   get the periodic info of all sections.
!-------------------------------------------------------------------------------
    subroutine get_sec_per_info
    use var_kind_def
    use var_global, only: is_has_cfg,cfg_file,pi
    use var_mesh, only: n_sec_per_info,sec_per_info
    use var_parallel
    implicit none
    character(len=80):: str,s
    integer(dpI):: io_err,i
    real   (dpR):: rtmp

    if(.not. is_has_cfg)    return
    if(myid .eq. 0) then
        open(unit=10,file=trim(cfg_file))
        do while(.true.)
            read(unit=10,fmt=*,iostat=io_err),str
            if(io_err .ne. 0)   exit

            s   =  str
            call upper_string(s)
            if(trim(s) .eq. '&SEC_PER_INFO') then
                read(unit=10,fmt=*),n_sec_per_info
                if(n_sec_per_info .gt. 0)   allocate(sec_per_info(5,n_sec_per_info))
                do i=1,n_sec_per_info
                    read(unit=10,fmt=*),sec_per_info(1:5,i)

                    if(nint(sec_per_info(2,i), dpI) .eq. 1) then
                        sec_per_info(2,i)   =  1.0d0
                    elseif(nint(sec_per_info(2,i), dpI) .eq. 2) then
                        sec_per_info(2,i)   =  2.0d0
                        call norm_vec(3, sec_per_info(3,i), rtmp)
                        rtmp=  2.0d0*pi/real(nint(2.0d0*pi/rtmp, dpI), dpR)
                        sec_per_info(3:5,i) =  sec_per_info(3:5,i)*rtmp
                    else
                        stop 'Error: wrong sec_per_info.'
                    end if
                end do

                exit
            end if
        end do
        close(10)
    end if
    call mpi_bcast(n_sec_per_info, 1, mpi_dpI, 0, mpi_comm_world, mpi_err)
    if(n_sec_per_info .le. 0)   return
    if(myid .ne. 0) allocate(sec_per_info(5,n_sec_per_info))
    call mpi_bcast(sec_per_info, 5*n_sec_per_info, mpi_dpR, 0, mpi_comm_world, mpi_err)

    return
    end subroutine get_sec_per_info
!-------------------------------------------------------------------------------
!   get the motion info of all sections.
!-------------------------------------------------------------------------------
    subroutine get_sec_motion_info
    use var_kind_def
    use var_global, only: is_has_cfg,cfg_file,pi
    use var_mesh
    use var_parallel
    implicit none
    character(len=80):: str,s
    integer(dpI):: io_err,i,n,isec,itype
    real   (dpR):: speed

    allocate(sec_motion_type (n_sec_g))
    allocate(sec_motion_speed(n_sec_g))
    sec_motion_type =  0
    sec_motion_speed=  0.0d0

    if(.not. is_has_cfg)    return
    if(myid .eq. 0) then
        open(unit=10,file=trim(cfg_file))
        do while(.true.)
            read(unit=10,fmt=*,iostat=io_err),str
            if(io_err .ne. 0)   exit

            s   =  str
            call upper_string(s)
            if(trim(s) .eq. '&SEC_MOTION_INFO') then
                read(unit=10,fmt=*),n
                do i=1,n
                    read(unit=10,fmt=*),isec,itype,speed

                    if((isec .le. 0) .or. (isec .gt. n_sec_g))  cycle
                    if(itype .eq. 1) then
                        sec_motion_type (isec)  =  1
                        sec_motion_speed(isec)  =  speed
                    elseif(itype .eq. 2) then
                        sec_motion_type (isec)  =  2
                        sec_motion_speed(isec)  =  speed*pi/3.0d1
                    end if
                end do

                exit
            end if
        end do
        close(10)
    end if
    call mpi_bcast(sec_motion_type , n_sec_g, mpi_dpI, 0, mpi_comm_world, mpi_err)
    call mpi_bcast(sec_motion_speed, n_sec_g, mpi_dpR, 0, mpi_comm_world, mpi_err)

    return
    end subroutine get_sec_motion_info
!-------------------------------------------------------------------------------
!   set the mortar.
!-------------------------------------------------------------------------------
    subroutine set_mortar
    use var_kind_def
    use var_cgns
    use var_global, only: is_structured_stencil,is_hanging_node,err_mem,n_dim,sw_slave
    use var_mesh
    use var_parallel
    use var_slv, only: 
    implicit none
    logical(dpL):: ltmp,ltmv(10)
    integer(dpI):: isec,iele,i,j,k,M,v(8),im,nfac,ifac,npe,L,R,LR(2,5),vL(2,8), &
                &  n_ele_R,n_f2n,s_LL,e_LL,sL,eL,sR,eR,s_RR,e_RR,n_ID_g2L,f(4,6), &
                &  imh,itmv(10)
    integer(dpI),allocatable:: fac(:,:),f2n(:,:),ss(:,:),ID_g2L(:)
    integer(dpI):: istart,iend,idx,n_ele_i,ref_val,owner

!   ----------------------------------------------------------------------------
!   record the faces.
    i   =  0
    do isec=mesh(0)%sec_1,mesh(0)%sec_0
        if(sec(isec)%is_bnd)    cycle
        i   =  i+sec(isec)%n_ele*nface_ele(sec(isec)%ele_type)
    end do
    allocate(fac(9,i), stat=err_mem)
    if(err_mem .ne. 0)  stop 'Error: fails to allocate memory, set mortar.'
    fac =  0

    nfac=  0
    k   =  huge(1)-1
    do isec=mesh(0)%sec_1,mesh(0)%sec_0
        if(sec(isec)%is_bnd)    cycle
        if(sec(isec)%is_quad) then
            do iele=1,sec(isec)%n_ele
                v(1:4)          =  sec(isec)%n2e(1:4,iele)
                fac(1:6,nfac+1) = (/v(1), v(2), k, k, isec, iele/)
                fac(1:6,nfac+2) = (/v(2), v(3), k, k, isec, iele/)
                fac(1:6,nfac+3) = (/v(3), v(4), k, k, isec, iele/)
                fac(1:6,nfac+4) = (/v(4), v(1), k, k, isec, iele/)
                nfac            =  nfac+4
            end do
        elseif(sec(isec)%is_tri) then
            do iele=1,sec(isec)%n_ele
                v(1:3)          =  sec(isec)%n2e(1:3,iele)
                fac(1:6,nfac+1) = (/v(1), v(2), k, k, isec, iele/)
                fac(1:6,nfac+2) = (/v(2), v(3), k, k, isec, iele/)
                fac(1:6,nfac+3) = (/v(3), v(1), k, k, isec, iele/)
                nfac            =  nfac+3
            end do
        elseif(sec(isec)%is_hexa) then
            do iele=1,sec(isec)%n_ele
                v(1:8)          =  sec(isec)%n2e(1:8,iele)
                fac(1:6,nfac+1) = (/v(1), v(4), v(3), v(2), isec, iele/)
                fac(1:6,nfac+2) = (/v(5), v(6), v(7), v(8), isec, iele/)
                fac(1:6,nfac+3) = (/v(1), v(5), v(8), v(4), isec, iele/)
                fac(1:6,nfac+4) = (/v(2), v(3), v(7), v(6), isec, iele/)
                fac(1:6,nfac+5) = (/v(1), v(2), v(6), v(5), isec, iele/)
                fac(1:6,nfac+6) = (/v(4), v(8), v(7), v(3), isec, iele/)
                nfac            =  nfac+6
            end do
        elseif(sec(isec)%ele_type .eq. TETRA_4) then
            do iele=1,sec(isec)%n_ele
                v(1:4)          =  sec(isec)%n2e(1:4,iele)
                fac(1:6,nfac+1) = (/v(1), v(3), v(2), k, isec, iele/)
                fac(1:6,nfac+2) = (/v(1), v(2), v(4), k, isec, iele/)
                fac(1:6,nfac+3) = (/v(2), v(3), v(4), k, isec, iele/)
                fac(1:6,nfac+4) = (/v(3), v(1), v(4), k, isec, iele/)
                nfac            =  nfac+4
            end do
        elseif(sec(isec)%ele_type .eq. PYRA_5) then
            do iele=1,sec(isec)%n_ele
                v(1:5)          =  sec(isec)%n2e(1:5,iele)
                fac(1:6,nfac+1) = (/v(1), v(4), v(3), v(2), isec, iele/)
                fac(1:6,nfac+2) = (/v(1), v(2), v(5), k   , isec, iele/)
                fac(1:6,nfac+3) = (/v(2), v(3), v(5), k   , isec, iele/)
                fac(1:6,nfac+4) = (/v(3), v(4), v(5), k   , isec, iele/)
                fac(1:6,nfac+5) = (/v(4), v(1), v(5), k   , isec, iele/)
                nfac            =  nfac+5
            end do
        elseif(sec(isec)%ele_type .eq. PENTA_6) then
            do iele=1,sec(isec)%n_ele
                v(1:6)          =  sec(isec)%n2e(1:6,iele)
                fac(1:6,nfac+1) = (/v(1), v(2), v(5), v(4), isec, iele/)
                fac(1:6,nfac+2) = (/v(2), v(3), v(6), v(5), isec, iele/)
                fac(1:6,nfac+3) = (/v(3), v(1), v(4), v(6), isec, iele/)
                fac(1:6,nfac+4) = (/v(1), v(3), v(2), k   , isec, iele/)
                fac(1:6,nfac+5) = (/v(4), v(5), v(6), k   , isec, iele/)
                nfac            =  nfac+5
            end do
        else
            stop 'Error: element type not supported.'
        end if
    end do
    do i=1,nfac
        call iqsort(.true., 1, 4, fac(1,i))
    end do
    call iqsortcols(.true., 1, nfac, 1, 9, fac)
!   record the faces.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   find the faces covered by boundary condition.
    do isec=mesh(0)%sec_1,mesh(0)%sec_0
        if(.not. sec(isec)%is_bnd)  cycle
        npe =  npe_ele_1st(sec(isec)%ele_type)
        do iele=1,sec(isec)%n_ele
            v       =  huge(1)-1
            v(1:npe)=  sec(isec)%n2e(1:npe,iele)
            call find_col_in_matrix(9, 1, 4, 1, nfac, fac, v, i, itmv)
            if((i .le. 0) .or. (i .gt. 1))  stop 'Error: fails to map face to element.'

            fac(7,itmv(1))  =  1
            fac(8,itmv(1))  =  isec
            fac(9,itmv(1))  =  iele
        end do
    end do
!   find the faces covered by boundary condition.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   find internal faces.
    do ifac=1,nfac
        if(fac(7,ifac) .gt. 0)  cycle
        v(1:4)  =  fac(1:4,ifac)
        call find_col_in_matrix(9, 1, 4, 1, nfac, fac, v, i, itmv)
        if(i .le. 1) then
!           do nothing.
        elseif(i .gt. 2) then
            stop 'Error: more than 2 faces matched, set_mortar.'
        else
            if(itmv(1) .eq. ifac) then
                L   =  itmv(2)
            else
                L   =  itmv(1)
            end if
            fac(8:9,ifac)   =  fac(5:6,L)
            if(sec(fac(5,ifac))%is_ghost .and. sec(fac(5,L))%is_ghost) then
!               do nothing.
            else
                fac(7,ifac) =  2
                fac(7,L   ) = -1
            end if
        end if
    end do
    do ifac=1,nfac
        if(fac(7,ifac) .ne. 2)  cycle
        sL  =  fac(5,ifac)
        eL  =  fac(6,ifac)
        sR  =  fac(8,ifac)
        eR  =  fac(9,ifac)
        if(.not. (sec(sL)%is_int .and. sec(sR)%is_int)) cycle
        if(sec(sL)%ID_ele_i(eL) .le. sec(sR)%ID_ele_i(eR))  cycle
        fac(5:6,ifac)   = (/sR, eR/)
        fac(8:9,ifac)   = (/sL, eL/)
    end do
!   find internal faces.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   record the elements used by myid.
    i   =  0
    do isec=mesh(0)%sec_1,mesh(0)%sec_0
        if(sec(isec)%is_int) then
            i   =  i+sec(isec)%n_ele
        elseif(sec(isec)%is_ghost) then
            do j=1,sec(isec)%n_ele
                if(sec(isec)%per_path(1,j) .eq. 0)  i   =  i+1
            end do
        end if
    end do
    i   =  i*3
    if(i .gt. 0)    allocate(ID_g2L(i), stat=err_mem)
    i   =  0
    do isec=mesh(0)%sec_1,mesh(0)%sec_0
        if(sec(isec)%is_int) then
            do j=1,sec(isec)%n_ele
                i   =  i+1
                ID_g2L(3*i-2:3*i)   = (/sec(isec)%ID_ele_i(j), isec, j/)
            end do
        elseif(sec(isec)%is_ghost) then
            do j=1,sec(isec)%n_ele
                if(sec(isec)%per_path(1,j) .gt. 0)  cycle
                i   =  i+1
                ID_g2L(3*i-2:3*i)   = (/sec(isec)%ID_ele_i(j), isec, j/)
            end do
        end if
    end do
    n_ID_g2L=  i
    call iqsortcols(.true., 1, n_ID_g2L, 1, 3, ID_g2L)
!   record the elements used by myid.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   record the face-to-vertex information.
    n_f2n   =  0
    do i=1,nfac
        if(fac(7,i) .ne. 0) cycle
        npe =  0
        do j=1,4
            if(fac(j,i) .eq. huge(1)-1) then
                exit
            else
                npe =  npe+1
            end if
        end do
        n_f2n   =  n_f2n+npe
    end do
    if(n_f2n .gt. 0)    allocate(f2n(2,n_f2n), stat=err_mem)
    n_f2n   =  0
    do i=1,nfac
        if(fac(7,i) .ne. 0) cycle
        do j=1,4
            if(fac(j,i) .eq. huge(1)-1) then
                exit
            else
                n_f2n           =  n_f2n+1
                f2n(1:2,n_f2n)  = (/fac(j,i), i/)
            end if
        end do
    end do
    if(n_f2n .gt. 0)    call iqsortcols(.true., 1, n_f2n, 1, 2, f2n)
!   record the face-to-vertex information.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   find face with hanging nodes.
    do i=1,n_LR_hanging
        k   =  0
        M   =  0
        do j=1,5
            if(LR_hanging(j,i) .le. 0) then
                exit
            else
                k   =  k+1
            end if
            call ib_search(1, n_ID_g2L, 1, 3, ID_g2L, LR_hanging(j,i), ltmv(j), L, R)
            if(.not. ltmv(j))   cycle
            if(L .ne. R)    stop 'Error: set_mortar fails.'
            LR(1:2,j)   =  ID_g2L(3*L-1:3*L)
            if(sec(LR(1,j))%is_int) M   =  M+1
        end do
        if(k .le. 0)    stop 'Error: LR_hanging record is not correct.'
!       if none of the L/R element belongs to the internal element, no need to proceed.
        if(M .le. 0)    cycle

        do j=1,k
            if(.not. ltmv(j))   stop 'Error: L/R element of hanging face is not found.'
        end do
        n_ele_R =  k-1

        npe =  sec(LR(1,1))%npe
        vL(1,1:npe) =  sec(LR(1,1))%n2e(1:npe,LR(2,1))
        vL(2,1:npe) =  0
        call iqsortcols(.true., 1, npe, 1, 2, vL)
        do j=2,k
        do M=1,sec(LR(1,j))%npe
            call ib_search(1, npe, 1, 2, vL, sec(LR(1,j))%n2e(M,LR(2,j)), ltmp, L, R)
            if(ltmp)    vL(2,L) =  vL(2,L)+1
        end do
        end do
        call iqsortcols(.false., 1, npe, 2, 2, vL)
        j   =  0
        do k=1,npe
            if(vL(2,k) .ne. 1) then
                exit
            else
                j   =  j+1
            end if
        end do
        if(n_dim .eq. 2) then
            ltmp=  j .ne. 2
        else
            ltmp= (j .ne. 3) .and. (j .ne. 4)
        end if
        if(ltmp)    stop 'Error: fails to get the L-vertex of hanging face.'

        v       =  huge(1)-1
        v(1:j)  =  vL(1,1:j)
        call find_col_in_matrix(9, 1, 4, 1, nfac, fac, v, M, itmv)
        if(M .ne. 1)    stop 'Error: fails to find L/R hanging face.'
!       de-activate the Left side hanging face.
        fac(7,itmv(1))  = -1

        call hanging_face(j, v, n_ele_R, LR(1,2), M, itmv)
        do j=1,M
            k   =  itmv(j)
            fac(8:9,itmv(j))=  LR(1:2,1)
            if(.not. (sec(fac(5,k))%is_int .or. sec(fac(8,k))%is_int)) then
                fac(7,itmv(j))  = -1
            else
                fac(7,itmv(j))  =  3
            end if
        end do
    end do
    if(allocated(f2n       ))   deallocate(f2n       )
    if(allocated(ID_g2L    ))   deallocate(ID_g2L    )
    if(allocated(LR_hanging))   deallocate(LR_hanging)
!   find face with hanging nodes.
!   ----------------------------------------------------------------------------

    mesh(0)%n_mortar_b  =  0
    mesh(0)%n_mortar    =  0
    do ifac=1,nfac
        i   =  fac(7,ifac)
        if((i .lt. 1) .or. (i .gt. 3))  cycle
        mesh(0)%n_mortar=  mesh(0)%n_mortar+1
        if(i .eq. 1)    mesh(0)%n_mortar_b  =  mesh(0)%n_mortar_b+1
    end do
    allocate(mesh(0)%mortar_ele_type(mesh(0)%n_mortar), stat=err_mem)
    allocate(mesh(0)%mortar_LR    (4,mesh(0)%n_mortar), stat=err_mem)
if(sw_slave) then
    allocate(mesh(0)%mortar_order (10,mesh(0)%n_mortar), stat=err_mem)
    mesh(0)%mortar_order = huge(1)-1
end if
    mesh(0)%mortar_LR   =  0
    if(n_dim .eq. 2) then
        allocate(mesh(0)%mortar_n2e(2,mesh(0)%n_mortar), stat=err_mem)
    else
        allocate(mesh(0)%mortar_n2e(4,mesh(0)%n_mortar), stat=err_mem)
    end if
    mesh(0)%mortar_n2e  =  0

!   ----------------------------------------------------------------------------
!   record the mortar covered by boundary condition.
    im  =  0
    do ifac=1,nfac
        if(fac(7,ifac) .ne. 1)  cycle
        im  =  im+1
        mesh(0)%mortar_LR(1:2,im)   =  fac(5:6,ifac)
        mesh(0)%mortar_LR(3:4,im)   =  fac(8:9,ifac)
        if(n_dim .eq. 2) then
            npe                         =  2
            mesh(0)%mortar_ele_type(im) =  BAR_2
        else
            npe =  0
            do i=1,4
                if(fac(i,ifac) .lt. huge(1)-1) then
                    npe =  npe+1
                else
                    exit
                end if
            end do
            if(npe .eq. 3) then
                mesh(0)%mortar_ele_type(im) =  TRI_3
            else
                mesh(0)%mortar_ele_type(im) =  QUAD_4
            end if
        end if
!       mesh(0)%mortar_n2e(1:npe,im)=  fac(1:npe,ifac)
        mesh(0)%mortar_n2e(1:npe,im)=  sec(fac(8,ifac))%n2e(1:npe,fac(9,ifac))
    end do
!   record the mortar covered by boundary condition.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   record the structured stencil.
    if(is_structured_stencil) then
        i   =  0
        do ifac=1,nfac
            if(fac(7,ifac) .eq. 2)  i   =  i+1
        end do
        if(i .gt. 0)    allocate(ss(9,i), stat=err_mem)

        n_ss=  0
        do ifac=1,nfac
            if(fac(7,ifac) .ne. 2)  cycle
            if(fac(5,ifac) .le. fac(8,ifac)) then
                sL  =  fac(5,ifac)
                eL  =  fac(6,ifac)
                sR  =  fac(8,ifac)
                eR  =  fac(9,ifac)
            else
                sL  =  fac(8,ifac)
                eL  =  fac(9,ifac)
                sR  =  fac(5,ifac)
                eR  =  fac(6,ifac)
            end if
            npe =  0
            do i=1,4
                if(fac(i,ifac) .lt. huge(1)-1) then
                    npe =  npe+1
                else
                    exit
                end if
            end do

            if(n_dim .eq. 2) then
                ltmp=  sec(sL)%ele_type .eq. QUAD_4
            else
                ltmp=((sec(sL)%ele_type .eq. PENTA_6) .and. (npe .eq. 3)) .or. &
                    & (sec(sL)%ele_type .eq. HEXA_8)
            end if
            if(.not. ltmp)  cycle

            if(n_dim .eq. 2) then
                ltmp=  sec(sR)%ele_type .eq. QUAD_4
            else
                ltmp=((sec(sR)%ele_type .eq. PENTA_6) .and. (npe .eq. 3)) .or. &
                    & (sec(sR)%ele_type .eq. HEXA_8)
            end if
            if(.not. ltmp)  cycle

            s_LL=  0
            e_LL=  0
            s_RR=  0
            e_RR=  0

            L       =  sec(sL)%npe
            v(1:L)  =  sec(sL)%n2e(1:L, eL)
            call list_minus_list(L, v, npe, fac(1,ifac))
            v(L+1:4)=  huge(1)-1
            call find_col_in_matrix(9, 1, 4, 1, nfac, fac, v, k, itmv)
            if(k .le. 0)    stop 'Error: fails to set structured stencil, set_mortar, 1.'
            do i=1,k
                if((fac(7,itmv(i)) .eq. 1) .or. (fac(7,itmv(i)) .eq. 3))    cycle
                if(fac(8,itmv(i)) .le. 0)   cycle
                if((fac(5,itmv(i)) .eq. sL) .and. (fac(6,itmv(i)) .eq. eL)) then
                    s_LL=  fac(8,itmv(i))
                    e_LL=  fac(9,itmv(i))
                    exit
                elseif((fac(8,itmv(i)) .eq. sL) .and. (fac(9,itmv(i)) .eq. eL)) then
                    s_LL=  fac(5,itmv(i))
                    e_LL=  fac(6,itmv(i))
                    exit
                end if
            end do
            if((s_LL .le. 0) .and. (k .eq. 1))  cycle
            if(s_LL .le. 0) stop 'Error: fails to set structured stencil, set_mortar, 2.'

            R       =  sec(sR)%npe
            v(1:R)  =  sec(sR)%n2e(1:R, eR)
            call list_minus_list(R, v, npe, fac(1,ifac))
            v(R+1:4)=  huge(1)-1
            call find_col_in_matrix(9, 1, 4, 1, nfac, fac, v, k, itmv)
            if(k .le. 0)    stop 'Error: fails to set structured stencil, set_mortar, 3.'
            do i=1,k
                if((fac(7,itmv(i)) .eq. 1) .or. (fac(7,itmv(i)) .eq. 3))    cycle
                if(fac(8,itmv(i)) .le. 0)   cycle
                if((fac(5,itmv(i)) .eq. sR) .and. (fac(6,itmv(i)) .eq. eR)) then
                    s_RR=  fac(8,itmv(i))
                    e_RR=  fac(9,itmv(i))
                    exit
                elseif((fac(8,itmv(i)) .eq. sR) .and. (fac(9,itmv(i)) .eq. eR)) then
                    s_RR=  fac(5,itmv(i))
                    e_RR=  fac(6,itmv(i))
                    exit
                end if
            end do
            if((s_RR .le. 0) .and. (k .eq. 1))  cycle
            if(s_RR .le. 0) stop 'Error: fails to set structured stencil, set_mortar, 4.'

            n_ss        =  n_ss+1
            ss(1  ,n_ss)=  ifac
            ss(2:3,n_ss)= (/s_LL, e_LL/)
            ss(4:5,n_ss)= (/sL  , eL  /)
            ss(6:7,n_ss)= (/sR  , eR  /)
            ss(8:9,n_ss)= (/s_RR, e_RR/)
            fac(7,ifac) = -1
        end do

        if(n_ss .gt. 0) then
            mesh(0)%n_mortar_ss =  im+n_ss
            allocate(mesh(0)%mortar_structured_stencil(8,im+1:im+n_ss), stat=err_mem)
            do i=1,n_ss
                im  =  im+1
                mesh(0)%mortar_LR(1:4,im)   =  ss(4:7,i)
                mesh(0)%mortar_structured_stencil(1:8,im)   =  ss(2:9,i)

                ifac=  ss(1,i)
                if(n_dim .eq. 2) then
                    npe                         =  2
                    mesh(0)%mortar_ele_type(im) =  BAR_2
                else
                    npe =  0
                    do j=1,4
                        if(fac(j,ifac) .lt. huge(1)-1) then
                            npe =  npe+1
                        else
                            exit
                        end if
                    end do
                    if(npe .eq. 3) then
                        mesh(0)%mortar_ele_type(im) =  TRI_3
                    else
                        mesh(0)%mortar_ele_type(im) =  QUAD_4
                    end if
                end if
                mesh(0)%mortar_n2e(1:npe,im)=  fac(1:npe,ifac)
            end do
        else
            mesh(0)%n_mortar_ss =  mesh(0)%n_mortar_b
        end if

        if(allocated(ss))   deallocate(ss)
    else
        mesh(0)%n_mortar_ss =  mesh(0)%n_mortar_b
    end if

    do i=1+mesh(0)%n_mortar_b,mesh(0)%n_mortar_ss
        ltmp=  sec(mesh(0)%mortar_structured_stencil(1,i))%is_bnd .or. &
            &  sec(mesh(0)%mortar_structured_stencil(3,i))%is_bnd .or. &
            &  sec(mesh(0)%mortar_structured_stencil(5,i))%is_bnd .or. &
            &  sec(mesh(0)%mortar_structured_stencil(7,i))%is_bnd
        if(ltmp)    stop 'Error: structured stencil is not correct.'
    end do
!   record the structured stencil.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   record the mortar defining the interface with hanging node.
    if(is_hanging_node) then
        i   =  0
        do ifac=1,nfac
            if(fac(7,ifac) .eq. 3)  i   =  i+1
        end do
        mesh(0)%n_mortar_hanging=  i
        if(i .gt. 0)    allocate(mesh(0)%mortar_hanging(i), stat=err_mem)
    end if
!   record the mortar defining the interface with hanging node.
!   ----------------------------------------------------------------------------

!   ----------------------------------------------------------------------------
!   record the mortar between element and element/ghost.
if(sw_slave) then
    ! record the face information
    n_ele_i = 0
    do ifac = 1,nfac
        if((fac(7,ifac) .ne. 2) .and. (fac(7,ifac) .ne. 3)) cycle
        n_ele_i  =  n_ele_i+1

        if(fac(5,ifac) .lt. fac(8,ifac)) then
            sL  =  fac(5,ifac)
            eL  =  fac(6,ifac)
            sR  =  fac(8,ifac)
            eR  =  fac(9,ifac)
        else if(fac(5,ifac) .eq. fac(8,ifac) &
                    & .and. fac(6,ifac) .le. fac(9,ifac)) then
            sL  =  fac(5,ifac)
            eL  =  fac(6,ifac)
            sR  =  fac(8,ifac)
            eR  =  fac(9,ifac)
        else
            sL  =  fac(8,ifac)
            eL  =  fac(9,ifac)
            sR  =  fac(5,ifac)
            eR  =  fac(6,ifac)
        end if
        mesh(0)%mortar_order(1 ,n_ele_i) = fac(1,ifac)
        mesh(0)%mortar_order(2 ,n_ele_i) = fac(2,ifac)
        mesh(0)%mortar_order(3 ,n_ele_i) = fac(3,ifac)
        mesh(0)%mortar_order(4 ,n_ele_i) = fac(4,ifac)
        mesh(0)%mortar_order(5 ,n_ele_i) = sL
        mesh(0)%mortar_order(6 ,n_ele_i) = eL
        mesh(0)%mortar_order(7 ,n_ele_i) = sR
        mesh(0)%mortar_order(8 ,n_ele_i) = eR
        mesh(0)%mortar_order(9 ,n_ele_i) = sec(sL)%ID_ele_g(eL)
        mesh(0)%mortar_order(10,n_ele_i) = sec(sR)%ID_ele_g(eR)
    end do
    
    ! sort faces according to owner
    call iqsortcols(.true., 1, n_ele_i, 9, 10, mesh(0)%mortar_order)

    ! sort faces according to neighbor
    istart = 1
    iend   = 1
    ref_val = mesh(0)%mortar_order(9,1)
    do idx=2,n_ele_i
        owner = mesh(0)%mortar_order(9,idx)
        if(owner .eq. ref_val) then
            iend = idx
            cycle
        else
            call iqsortcols(.true., istart, iend, 10, 10, &
                        & mesh(0)%mortar_order)
            istart  = idx
            iend    = idx
            ref_val = owner
        end if
    end do
            
    idx = 0
    do ifac=1,nfac
        if((fac(7,ifac) .ne. 2) .and. (fac(7,ifac) .ne. 3)) cycle
        im  =  im+1
        idx =  idx+1

        if(fac(7,ifac) .eq. 3) then
            imh =  imh+1
            mesh(0)%mortar_hanging(imh) =  im
        end if

        sL = mesh(0)%mortar_order(5,idx)
        eL = mesh(0)%mortar_order(6,idx)
        sR = mesh(0)%mortar_order(7,idx)
        eR = mesh(0)%mortar_order(8,idx)

        if(n_dim .eq. 2) then
            npe                         =  2
            mesh(0)%mortar_ele_type(im) =  BAR_2
        else
            npe =  0
            do i=1,4
                if(mesh(0)%mortar_order(i,idx) .lt. huge(1)-1) then
                    npe =  npe+1
                else
                    exit
                end if
            end do
            if(npe .eq. 3) then
                mesh(0)%mortar_ele_type(im) =  TRI_3
            else
                mesh(0)%mortar_ele_type(im) =  QUAD_4
            end if
        end if
        mesh(0)%mortar_n2e(1:npe,im)=  mesh(0)%mortar_order(1:npe,idx)
        mesh(0)%mortar_LR (1:4  ,im)= (/sL, eL, sR, eR/)
    end do
else
    imh =  0
    do ifac=1,nfac
        if((fac(7,ifac) .ne. 2) .and. (fac(7,ifac) .ne. 3)) cycle
        im  =  im+1

        if(fac(7,ifac) .eq. 3) then
            imh =  imh+1
            mesh(0)%mortar_hanging(imh) =  im
        end if

        if(fac(5,ifac) .lt. fac(8,ifac)) then
            sL  =  fac(5,ifac)
            eL  =  fac(6,ifac)
            sR  =  fac(8,ifac)
            eR  =  fac(9,ifac)
        else if(fac(5,ifac) .eq. fac(8,ifac) &
                    & .and. fac(6,ifac) .le. fac(9,ifac)) then
            sL  =  fac(5,ifac)
            eL  =  fac(6,ifac)
            sR  =  fac(8,ifac)
            eR  =  fac(9,ifac)
        else
            sL  =  fac(8,ifac)
            eL  =  fac(9,ifac)
            sR  =  fac(5,ifac)
            eR  =  fac(6,ifac)
        end if

        if(n_dim .eq. 2) then
            npe                         =  2
            mesh(0)%mortar_ele_type(im) =  BAR_2
        else
            npe =  0
            do i=1,4
                if(fac(i,ifac) .lt. huge(1)-1) then
                    npe =  npe+1
                else
                    exit
                end if
            end do
            if(npe .eq. 3) then
                mesh(0)%mortar_ele_type(im) =  TRI_3
            else
                mesh(0)%mortar_ele_type(im) =  QUAD_4
            end if
        end if
        mesh(0)%mortar_n2e(1:npe,im)=  fac(1:npe,ifac)
        mesh(0)%mortar_LR (1:4  ,im)= (/sL, eL, sR, eR/)
    end do
end if

!   record the mortar between element and element/ghost.
!   ----------------------------------------------------------------------------
!   print*,myid,mesh(0)%n_mortar_b,mesh(0)%n_mortar_ss,mesh(0)%n_mortar

    if(im .ne. mesh(0)%n_mortar)    stop 'Error: set_mortar fails.'

    do im=1,mesh(0)%n_mortar
        sL  =  mesh(0)%mortar_LR(1,im)
!       sR  =  mesh(0)%mortar_LR(3,im)
        if(.not. sec(sL)%is_int)    stop 'Error: left side of mortar is not INTERNAL.'
    end do

!   ----------------------------------------------------------------------------
!   recover the vertex order for QUAD_4 type mortar.
    k   =  huge(1)-1
    do im=1+mesh(0)%n_mortar_b,mesh(0)%n_mortar
        if(mesh(0)%mortar_ele_type(im) .ne. QUAD_4) cycle
        sL  =  mesh(0)%mortar_LR(1,im)
        eL  =  mesh(0)%mortar_LR(2,im)
        sR  =  mesh(0)%mortar_LR(3,im)
        eR  =  mesh(0)%mortar_LR(4,im)
        M   =  ele_type_1st(sec(sL)%ele_type)
        npe =  npe_ele_1st (sec(sL)%ele_type)
        v(1:npe)=  sec(sL)%n2e(1:npe,eL)

        if((M .eq. TRI_3)) then
            nfac    =  0
        elseif((M .eq. QUAD_4) .or. (M .eq. QUAD_8) .or. (M .eq. QUAD_9)) then
            nfac    =  0
        elseif(M .eq. TETRA_4) then
            nfac    =  0
        elseif(M .eq. PYRA_5) then
            nfac    =  1
            f(1:4,1)= (/v(1), v(2), v(3), v(4)/)
        elseif(M .eq. PENTA_6) then
            nfac    =  3
            f(1:4,1)= (/v(1), v(2), v(5), v(4)/)
            f(1:4,2)= (/v(2), v(3), v(6), v(5)/)
            f(1:4,3)= (/v(1), v(3), v(6), v(4)/)
        elseif((M .eq. HEXA_8) .or. (M .eq. HEXA_20) .or. (M .eq. HEXA_27)) then
            nfac    =  6
            f(1:4,1)= (/v(1), v(2), v(3), v(4)/)
            f(1:4,2)= (/v(5), v(6), v(7), v(8)/)
            f(1:4,3)= (/v(1), v(5), v(8), v(4)/)
            f(1:4,4)= (/v(2), v(6), v(7), v(3)/)
            f(1:4,5)= (/v(1), v(2), v(6), v(5)/)
            f(1:4,6)= (/v(4), v(3), v(7), v(8)/)
        else
            stop 'Error: element type not supported.'
        end if
        if(nfac .le. 0) cycle

        ltmv(1) =  .true.
        do i=1,nfac
            v(1:4)  =  f(1:4,i)
            call iqsort(.true., 1, 4, v)
            ltmp=  .true.
            do j=1,4
                ltmp=  ltmp .and. (mesh(0)%mortar_n2e(j,im) .eq. v(j))
                if(.not. ltmp)  exit
            end do
            if(ltmp) then
                mesh(0)%mortar_n2e(1:4,im)  =  f(1:4,i)
                ltmv(1)                     =  .false.
                exit
            end if
        end do
        if(ltmv(1)) stop 'Error: fails to map mortar(QUAD_4) to internal element.'
    end do
!   recover the vertex order for QUAD_4 type mortar.
!   ----------------------------------------------------------------------------

    if(allocated(fac))  deallocate(fac)

    return
    contains
!   ----------------------------------------------------------------------------
!   find the R-part face.
!   ----------------------------------------------------------------------------
        subroutine hanging_face(npe,vL,n_ele_R,eR,n_fac_R,R_face)
        use var_kind_def
        implicit none
        integer(dpI),intent(in):: npe,vL(*),n_ele_R,eR(2,*)
        logical(dpL):: ltmp
        integer(dpI):: n_fac_R,R_face(*),i,j,k,L,R,f(0:100,4),M,v(16),freq(16), &
                    &  w(2,16),isec,iele,ifac,itmv(10)

        f   =  0
        do i=1,npe
            call ib_search(1, n_f2n, 1, 2, f2n, vL(i), ltmp, L, R)
            if(.not. ltmp)  stop 'Error: fails to find faces connected to hanging vertex.'
            do j=L,R
                ifac=  f2n(2,j)
                if(fac(7,ifac) .ne. 0)  cycle
                isec=  fac(5,ifac)
                iele=  fac(6,ifac)
                do k=1,n_ele_R
                    if((isec .eq. eR(1,k)) .and. (iele .eq. eR(2,k))) then
                        f(0,i)      =  f(0,i)+1
                        f(f(0,i),i) =  ifac
                        exit
                    end if
                end do
            end do
        end do

        if(npe .eq. 2) then
            n_fac_R =  2
        else
            n_fac_R =  4
        end if
        R_face(1:n_fac_R)   =  0

        if(npe .eq. 2) then
            loop_1:do j=1,f(0,2)
            do i=1,f(0,1)
                v(1:2)  =  fac(1:2,f(j,2))
                v(3:4)  =  fac(1:2,f(i,1))
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
                    R_face(1:2) = (/f(i,1), f(j,2)/)
                    exit loop_1
                end if
            end do
            end do loop_1
        elseif(npe .eq. 3) then
            loop_2:do k=1,f(0,3)
            do j=1,f(0,2)
            do i=1,f(0,1)
                v(1:3)  =  fac(1:3,f(k,3))
                v(4:6)  =  fac(1:3,f(j,2))
                v(7:9)  =  fac(1:3,f(i,1))
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

                call find_col_in_matrix(9, 1, 3, 1, nfac, fac, w(1,3:6), L, itmv)
                if(L .ne. 1)    cycle

                R_face(1:4) = (/f(i,1), f(j,2), f(k,3), itmv(1)/)
                exit loop_2
            end do
            end do
            end do loop_2
        else
            loop_3:do L=1,f(0,4)
            do k=1,f(0,3)
            do j=1,f(0,2)
            do i=1,f(0,1)
                v(1 :4 )=  fac(1:4,f(L,4))
                v(5 :8 )=  fac(1:4,f(k,3))
                v(9 :12)=  fac(1:4,f(j,2))
                v(13:16)=  fac(1:4,f(i,1))
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

                R_face(1:4) = (/f(i,1), f(j,2), f(k,3), f(L,4)/)
                exit loop_3
            end do
            end do
            end do
            end do loop_3
        end if
        do i=1,n_fac_R
            if(R_face(i) .le. 0)    stop 'Error: fails to find R_face of hanging face.'
        end do

        return
        end subroutine hanging_face
    end subroutine set_mortar
!-------------------------------------------------------------------------------
!   memory usage.
!-------------------------------------------------------------------------------
    subroutine get_memory_usage
    use var_kind_def
    use var_mesh
    use var_global, only: sw_slave
    implicit none
    character(len=32):: str(100)
    integer(dpI):: isec
    real   (dpR):: v(10)

    str(1)  = 'n_ele'
    str(2)  = 'iA_ls'
    str(3)  = 'jA_ls'
    str(4)  = 'coe_ls'
    str(5)  = 'iA_nb'
    str(6)  = 'jA_nb'
    open(unit=10,file='./data/memory.dat')
    write(10,fmt='(5A12)'),str(1:5)
    do isec=mesh(0)%sec_1,mesh(0)%sec_0
        v   =  0.0d0
        if(allocated(sec(isec)%iA_ls))  v(1)=  size(sec(isec)%iA_ls)*4.0d0
        if(allocated(sec(isec)%jA_ls))  v(2)=  size(sec(isec)%jA_ls)*4.0d0
        if(allocated(sec(isec)%coe_ls)) v(3)=  size(sec(isec)%coe_ls)*8.0d0
        if(allocated(sec(isec)%iA_face_neighbour))  v(4)=  size(sec(isec)%iA_face_neighbour)*4.0d0
        if(allocated(sec(isec)%jA_face_neighbour))  v(5)=  size(sec(isec)%jA_face_neighbour)*4.0d0
        write(unit=10,fmt='(I8,5ES12.4)'),sec(isec)%n_ele,v(1:5)/1.0d6
    end do
    v(1)=  size(mesh(0)%mortar_LR  )*4.0d0
    v(2)=  size(mesh(0)%mortar_cen )*8.0d0
    v(3)=  size(mesh(0)%mortar_n_vg)*8.0d0
    write(unit=10,fmt='(5ES12.4)'),v(1:4)/1.0d6
    close(10)

    return
    end subroutine get_memory_usage

    subroutine set_ID_ele_g(lev)
    use var_mesh
    implicit none

    integer(dpI):: isec,iele,nele
    integer(dpI),intent(in):: lev

    nele = 0
    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(.not. allocated(sec(isec)%ID_ele_g)) then
            allocate(sec(isec)%ID_ele_g(sec(isec)%n_ele))
        end if
        do iele=1,sec(isec)%n_ele
            sec(isec)%ID_ele_g(iele) = nele+iele
        end do
        nele = nele+sec(isec)%n_ele
    end do
    
    return
    end subroutine set_ID_ele_g

    subroutine set_mortar_ID(lev)
    use var_mesh
    use var_parallel
    use var_global, only: n_dim
    implicit none

    integer(dpI):: isec,iele,sL,eL,sR,eR,im,imi,icell,s_donor,e_donor,per(3)
    integer(dpI),intent(in):: lev
    integer(dpI),allocatable:: cellWeights(:), faceWeights(:),iA(:),jA(:)
    real   (dpR),allocatable:: A(:),transform(:)
    integer(dpI),allocatable:: owner(:),neighbor(:)

    faceNum = mesh(lev)%n_mortar-mesh(lev)%n_mortar_ss

    if(.not. allocated(mesh(lev)%mortar_own_ID)) then
        allocate(mesh(lev)%mortar_own_ID &
                    & (mesh(lev)%n_mortar-mesh(lev)%n_mortar_ss))
        allocate(mesh(lev)%mortar_nei_ID &
                    & (mesh(lev)%n_mortar-mesh(lev)%n_mortar_ss))
    end if
    imi = 1
    cellNum = 1
! if(myid==0) write(*,*),mesh(lev)%n_mortar_b,mesh(lev)%n_mortar_ss,mesh(lev)%n_mortar
    do im=1+mesh(lev)%n_mortar_ss,mesh(lev)%n_mortar
        sL = mesh(lev)%mortar_LR(1,im)
        eL = mesh(lev)%mortar_LR(2,im)
        sR = mesh(lev)%mortar_LR(3,im)
        eR = mesh(lev)%mortar_LR(4,im)
        ! if(myid==0 .and. sL==2 .and. eL==3757) write(*,*),"left",im
        ! if(myid==0 .and. sR==2 .and. eR==3757) write(*,*),"right",im
!        if(sec(sL)%ID_ele_g(eL) .gt. sec(sR)%ID_ele_g(eR)) &
!& write(*,*),im,sL,sR,sec(sL)%ID_ele_g(eL),sec(sR)%ID_ele_g(eR)
        mesh(lev)%mortar_own_ID(imi) = sec(sL)%ID_ele_g(eL)
        mesh(lev)%mortar_nei_ID(imi) = sec(sR)%ID_ele_g(eR)
! if(myid==0) write(*,*),mesh(lev)%mortar_own_ID(imi)
        cellNum = max(cellNum, mesh(lev)%mortar_nei_ID(imi))
        cellNum = max(cellNum, mesh(lev)%mortar_own_ID(imi))
        imi = imi+1
! if(sec(sR)%is_ghost) then
    ! call get_ID_ghost_ele(sR, eR, s_donor, e_donor, per)
    ! write(*,*),myid,sR,im
! end if
    end do
if(myid==0) then
    write(*,*),cellNum
    do im = mesh(lev)%sec_1,mesh(lev)%sec_0
        write(*,*),im,sec(im)%n_ele,sec(im)%is_int,sec(im)%is_ghost,sec(im)%is_bnd,sec(im)%ele_type
    end do
end if

    allocate(iA(cellNum+2))
    allocate(jA(faceNum*2))
    iA = 0
    do im=1,faceNum
        iA(mesh(lev)%mortar_nei_ID(im)+2) = iA(mesh(lev)%mortar_nei_ID(im)+2)+1
        iA(mesh(lev)%mortar_own_ID(im)+2) = iA(mesh(lev)%mortar_own_ID(im)+2)+1
    end do
    iA(1) = 1
    iA(2) = 1
    do icell=1,cellNum
        iA(icell+2) = iA(icell+2)+iA(icell+1)
    end do
    do im=1,faceNum
        jA(iA(mesh(lev)%mortar_own_ID(im)+1)) = mesh(lev)%mortar_nei_ID(im)
        jA(iA(mesh(lev)%mortar_nei_ID(im)+1)) = mesh(lev)%mortar_own_ID(im)
        iA(mesh(lev)%mortar_own_ID(im)+1) = iA(mesh(lev)%mortar_own_ID(im)+1)+1
        iA(mesh(lev)%mortar_nei_ID(im)+1) = iA(mesh(lev)%mortar_nei_ID(im)+1)+1
    end do
    if(.not. allocated(perm)) allocate(perm(cellNum))
    if(.not. allocated(A))    allocate(A(faceNum))
    ! A = 2
    call matrix_ordering(cellNum, iA, jA, perm)

    ! The faces located in the upper triangle may be in the lower triangle 
    ! after reordering. Transform records the situation with value -1.
    if(.not. allocated(transform)) allocate(transform(faceNum))
    imi = 1
    transform = 0
    do icell=1,cellNum
        do im=iA(icell),iA(icell+1)-1
            if(icell .lt. jA(im)) then
                if(perm(icell) .gt. perm(jA(im))) then
                    transform(imi) = 1
                end if
                imi = imi + 1
            end if
        end do
    end do
    call matrix_reorder_ldu(.false., 1, cellNum, perm, iA, jA, transform)
    
    ! call topo_reorder_ldu(cellNum, perm, iA, jA, &
        ! & mesh(lev)%mortar_own_ID, mesh(lev)%mortar_nei_ID)

    ! reorder mesh
    call init_mesh_reorder(lev)
    ! call matrix_reorder_ldu(.false., 1, cellNum, perm, iA, jA, &
        ! & transform)
    call matrix_reorder_ldu(.false., 3, cellNum, perm, iA, jA, &
        & mesh_reordered(lev)%mortar_cen)
    call matrix_reorder_ldu(.true., 8, cellNum, perm, iA, jA, &
        & mesh_reordered(lev)%mortar_n_vg)
    ! if(myid==0) write(*,*),mesh_reordered(lev)%mortar_n_vg(1,1)

    ! if(myid==0) write(*,*),mesh_reordered(lev)%mortar_n_vg(1,1)
    ! call matrix_reorder_ldu(.true., cellNum, perm, iA, jA, &
        ! & mesh_reordered(lev)%mortar_nei_ID)
    if(.not. allocated(mesh_reordered(lev)%mortar_own_ID)) then
        allocate(mesh_reordered(lev)%mortar_own_ID(faceNum))
        allocate(mesh_reordered(lev)%mortar_nei_ID(faceNum))
        allocate(mesh_reordered(lev)%mortar_transform(faceNum))
        allocate(mesh_reordered(lev)%owner(faceNum))
        allocate(mesh_reordered(lev)%neighbor(faceNum))
    end if
    imi = 1
    mesh_reordered(lev)%mortar_transform(:) = transform(:)
    do icell=1,cellNum
        do im=iA(icell),iA(icell+1)-1
            ! write(*,*),iA(icell),icell,jA(im)
            if(icell .gt. jA(im)) cycle
            mesh_reordered(lev)%mortar_own_ID(imi) = icell
            mesh_reordered(lev)%mortar_nei_ID(imi) = jA(im)
            imi = imi+1
        end do
    end do

    mesh_reordered(lev)%owner(:)    = mesh_reordered(lev)%mortar_own_ID(:) - 1
    mesh_reordered(lev)%neighbor(:) = mesh_reordered(lev)%mortar_nei_ID(:) - 1
    ! do im=1,faceNum
        ! if(myid==0) write(*,*),mesh_reordered(lev)%neighbor(im)
    ! end do

    return 
    end subroutine set_mortar_ID

    subroutine init_mesh_reorder(lev)
    use var_kind_def
    use var_mesh
    use var_parallel
    use var_global, only: n_dim
    implicit none

    integer(dpI),intent(in):: lev
    integer(dpI):: im,imi

    if(.not. allocated(mesh_reordered(lev)%mortar_n_vg)) then
        allocate(mesh_reordered(lev)%mortar_n_vg(8, faceNum))
    end if
    if(.not. allocated(mesh_reordered(lev)%mortar_cen)) then
        allocate(mesh_reordered(lev)%mortar_cen(3, faceNum))
    end if

! if(myid .ne. 0) return
    imi = 0
    do im=1+mesh(lev)%n_mortar_ss,mesh(lev)%n_mortar
        imi = imi + 1
        ! write(*,*),im,mesh(lev)%n_mortar
        mesh_reordered(lev)%mortar_n_vg(:,imi) = mesh(lev)%mortar_n_vg(:,im)
        mesh_reordered(lev)%mortar_cen (:,imi) = mesh(lev)%mortar_cen (:,im)
    end do

    return
    end subroutine init_mesh_reorder

    subroutine mesh_reorder(lev)
    use var_kind_def
    use var_mesh
    use var_parallel
    use var_global, only: n_dim
    implicit none

    integer(dpI),intent(in):: lev
    integer(dpI):: im,imi

    ! call matrix_reorder_ldu(.false., n_dim, cellNum, perm, iA, jA, &
    !     & mesh_reordered(lev)%mortar_cen)
    ! call matrix_reorder_ldu(.true., 8, cellNum, perm, iA, jA, &
    !     & mesh_reordered(lev)%mortar_n_vg)

    return
    end subroutine mesh_reorder

!-------------------------------------------------------------------------------
!   get the ID of ghost element, if its donor beglongs to myid. Sunway version.
!-------------------------------------------------------------------------------
    subroutine set_ID_ghost_ele_sw(lev)
    use var_kind_def
    use var_mesh
    use var_sec_array
    use var_parallel, only: myid
    implicit none
    integer(dpI),intent(IN):: lev
    integer(dpI):: s,e,isr,i,ID,isec,iele

    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(.not. allocated(sec(isec)%ID_ghost_ele)) then
            allocate(sec(isec)%ID_ghost_ele(5,sec(isec)%n_ele))
        end if
        sec(isec)%ID_ghost_ele = 0
        if(.not. sec(isec)%is_ghost)    cycle
        do iele=1,sec(isec)%n_ele
            isr =  sec(isec)%id_recv(1,iele)
            i   =  sec(isec)%id_recv(2,iele)
            if(p2p(isr)%ip_remote .ne. myid)    return
            s   =  p2p(isr)%id_ele_send(1,i)
            e   =  p2p(isr)%id_ele_send(2,i)
            sec(isec)%ID_ghost_ele(1,iele) = s
            sec(isec)%ID_ghost_ele(2,iele) = e
        end do
    end do

    return
    end subroutine set_ID_ghost_ele_sw

    subroutine fv_set_lhs_imp_sw(lev,p)
    use var_mesh
    use var_prec
    use var_lhs_imp
    implicit none

    integer(dpI),intent(IN):: lev
    type(type_prec):: p
    integer(dpI):: isec,LDA,ivtx,im,j,iele,sR,eR,s_donor,e_donor,per(3)
    logical(dpL):: is_include_per

    n_ele_i = 0
    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(sec(isec)%is_int) then
            n_ele_i = n_ele_i+sec(isec)%n_ele
        end if
    end do

    is_include_per = .true.
    n_face_i =  0
    do isec=mesh(lev)%sec_1,mesh(lev)%sec_0
        if(.not. sec(isec)%is_int)  cycle
        do iele=1,sec(isec)%n_ele
            do j=sec(isec)%iA_face_neighbour(iele),sec(isec)%iA_face_neighbour(iele+1)-1
                im  =  sec(isec)%jA_face_neighbour(j)
                if(im .gt. 0) then
                    sR  =  mesh(lev)%mortar_LR(3, im)
                    eR  =  mesh(lev)%mortar_LR(4, im)
                else
                    cycle
                end if
                if(sec(sR)%is_int) then
                    n_face_i=  n_face_i+1
                elseif(sec(sR)%is_ghost .and. is_include_per) then
                    call get_ID_ghost_ele(sR, eR, s_donor, e_donor, per)
                    if(s_donor .gt. 0)  n_face_i =  n_face_i+1
                end if
            end do
        end do
    end do

    ! write(*,*),p%nvtx,n_ele_i,cellNum,size(p%A),n_face_i,faceNum

    if(.not. allocated(owner))        allocate(owner(faceNum))
    if(.not. allocated(neighbor))     allocate(neighbor(faceNum))
    if(.not. allocated(owner_int))    allocate(owner_int(n_face_i))
    if(.not. allocated(neighbor_int)) allocate(neighbor_int(n_face_i))

    LDA = p%bsize*p%bsize

    if(.not. allocated(upper))        allocate(upper(LDA,faceNum))
    if(.not. allocated(lower))        allocate(lower(LDA,faceNum))
    if(.not. allocated(diag))         allocate(diag (LDA,cellNum))
    if(.not. allocated(x))            allocate(x    (p%bsize,cellNum))
    if(.not. allocated(b))            allocate(b    (p%bsize,cellNum))

    if(.not. allocated(upper_int))    allocate(upper_int(LDA,n_face_i))
    if(.not. allocated(lower_int))    allocate(lower_int(LDA,n_face_i))
    if(.not. allocated(diag_int))     allocate(diag_int (LDA,n_ele_i))
    if(.not. allocated(x_int))        allocate(x_int    (p%bsize,n_ele_i))
    if(.not. allocated(b_int))        allocate(b_int    (p%bsize,n_ele_i))

    return
    end subroutine fv_set_lhs_imp_sw