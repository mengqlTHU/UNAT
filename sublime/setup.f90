!-------------------------------------------------------------------------------
!   setup.
!-------------------------------------------------------------------------------
    subroutine setup(is_start_mpi)
    use var_kind_def
    use var_slv, only: is_FV,is_heat
    use var_global, only: sw_slave
    implicit none
    logical(dpL),intent(in):: is_start_mpi

    call set_parallel(is_start_mpi)
    call set_cgns
    call set_global
    call set_unsteady
    call set_turbulence
    call read_bnd
    call read_mesh
    call read_bnd_profile

    if(is_fv .or. is_heat) then
        call set_vol_face(0)
        call mesh_agglomeration
    end if

if(sw_slave) then
    call set_mortar_ID(0)
end if

    return
    end subroutine setup
!-------------------------------------------------------------------------------
!   setup cgns.
!-------------------------------------------------------------------------------
    subroutine set_cgns
    use var_cgns
    implicit none

    nface_ele(BAR_2  )  =  0
    nface_ele(BAR_3  )  =  0
    nface_ele(TRI_3  )  =  3
    nface_ele(QUAD_4 )  =  4
    nface_ele(QUAD_8 )  =  4
    nface_ele(QUAD_9 )  =  4
    nface_ele(TETRA_4)  =  4
    nface_ele(PYRA_5 )  =  5
    nface_ele(PENTA_6)  =  5
    nface_ele(HEXA_8 )  =  6
    nface_ele(HEXA_20)  =  6
    nface_ele(HEXA_27)  =  6

    npe_ele(BAR_2  )=  2
    npe_ele(BAR_3  )=  3
    npe_ele(TRI_3  )=  3
    npe_ele(QUAD_4 )=  4
    npe_ele(QUAD_8 )=  8
    npe_ele(QUAD_9 )=  9
    npe_ele(TETRA_4)=  4
    npe_ele(PYRA_5 )=  5
    npe_ele(PENTA_6)=  6
    npe_ele(HEXA_8 )=  8
    npe_ele(HEXA_20)=  20
    npe_ele(HEXA_27)=  27

    npe_ele_1st(BAR_2  )=  2
    npe_ele_1st(BAR_3  )=  2
    npe_ele_1st(TRI_3  )=  3
    npe_ele_1st(QUAD_4 )=  4
    npe_ele_1st(QUAD_8 )=  4
    npe_ele_1st(QUAD_9 )=  4
    npe_ele_1st(TETRA_4)=  4
    npe_ele_1st(PYRA_5 )=  5
    npe_ele_1st(PENTA_6)=  6
    npe_ele_1st(HEXA_8 )=  8
    npe_ele_1st(HEXA_20)=  8
    npe_ele_1st(HEXA_27)=  8

    ele_type_1st(BAR_2  )   =  BAR_2
    ele_type_1st(BAR_3  )   =  BAR_2
    ele_type_1st(TRI_3  )   =  TRI_3
    ele_type_1st(QUAD_4 )   =  QUAD_4
    ele_type_1st(QUAD_8 )   =  QUAD_4
    ele_type_1st(QUAD_9 )   =  QUAD_4
    ele_type_1st(TETRA_4)   =  TETRA_4
    ele_type_1st(PYRA_5 )   =  PYRA_5
    ele_type_1st(PENTA_6)   =  PENTA_6
    ele_type_1st(HEXA_8 )   =  HEXA_8
    ele_type_1st(HEXA_20)   =  HEXA_8
    ele_type_1st(HEXA_27)   =  HEXA_8

    return
    end subroutine set_cgns
!-------------------------------------------------------------------------------
!   setup parameters.
!-------------------------------------------------------------------------------
    subroutine set_global
    use var_kind_def
    use var_bndv
    use var_global
    use var_parallel
    use var_slv
    implicit none
    logical(dpL):: lbuf(100)
    integer(dpI):: io_err,ibuf(1000)
    real   (dpR):: rbuf(100),rtmp

    namelist /model/    mesh_file,is_vis_cal,offset,is_show_res,is_output_cpcf, &
                &       translation_axis,rotation_axis,is_output_rel_uvw,zoom, &
                &       is_output_cgns,is_structured_stencil,is_hanging_node, &
                &       is_output_cc,span,body_force
    namelist /slv/      is_flow_int,is_FV,is_heat,max_iter_lev,CFL_lev,lev0_cal, &
                &       solver_lev,res_cvg_lev,is_line_implicit,prec_type,is_matrix_free

    if(myid .eq. 0) inquire(file=trim(adjustl(cfg_file)),exist=is_has_cfg)
    call mpi_bcast(is_has_cfg, 1, mpi_dpL, 0, mpi_comm_world, mpi_err)
    if(.not. is_has_cfg)    return

    if(myid .eq. 0) then
        open(unit=10,file=trim(adjustl(cfg_file)))

        read(unit=10, nml=model, iostat=io_err)
        if(io_err .gt. 0)   stop 'Error: fails to read namelist:model.'

        rewind(10)
        read(unit=10, nml=slv, iostat=io_err)
        if(io_err .gt. 0)   stop 'Error: fails to read namelist:slv.'

        close(10)
    end if

    lbuf(1:13)  = (/is_vis_cal, is_FV, is_flow_int, is_heat, is_show_res,is_output_cpcf, &
                &   is_output_rel_uvw, is_output_cgns, is_structured_stencil, &
                &   is_line_implicit, is_hanging_node, is_output_cc, is_matrix_free/)
    ibuf(1:22)  = (/lev0_cal, max_iter_lev, solver_lev, prec_type/)
    rbuf(1:36)  = (/res_cvg_lev, CFL_lev, offset, translation_axis, rotation_axis, &
                &   zoom, span, body_force/)
    call mpi_bcast(lbuf, 13, mpi_dpL, 0, mpi_comm_world, mpi_err)
    call mpi_bcast(ibuf, 22, mpi_dpI, 0, mpi_comm_world, mpi_err)
    call mpi_bcast(rbuf, 36, mpi_dpR, 0, mpi_comm_world, mpi_err)
    call mpi_bcast(mesh_file, 80, mpi_character, 0, mpi_comm_world, mpi_err)

    if(mesh_file(len_trim(mesh_file)-3:len_trim(mesh_file)) .eq. 'cgns') then
        mesh_name   =  mesh_file(1:len_trim(mesh_file)-5)
        slv_file    =  mesh_file(1:len_trim(mesh_file)-5)//'_slv.cgns'
        res_file    =  mesh_file(1:len_trim(mesh_file)-5)//'_res.dat'
    elseif(mesh_file(len_trim(mesh_file)-2:len_trim(mesh_file)) .eq. 'msh') then
        mesh_name   =  mesh_file(1:len_trim(mesh_file)-4)
        slv_file    =  mesh_file(1:len_trim(mesh_file)-4)//'_slv.cgns'
        res_file    =  mesh_file(1:len_trim(mesh_file)-4)//'_res.dat'
    end if

    if(myid .ne. 0) then
        is_vis_cal              =  lbuf(1)
        is_FV                   =  lbuf(2)
        is_flow_int             =  lbuf(3)
        is_heat                 =  lbuf(4)
        is_show_res             =  lbuf(5)
        is_output_cpcf          =  lbuf(6)
        is_output_rel_uvw       =  lbuf(7)
        is_output_cgns          =  lbuf(8)
        is_structured_stencil   =  lbuf(9)
        is_line_implicit        =  lbuf(10)
        is_hanging_node         =  lbuf(11)
        is_output_cc            =  lbuf(12)
        is_matrix_free          =  lbuf(13)
        lev0_cal                =  ibuf( 1   )
        max_iter_lev(0:9)       =  ibuf( 2:11)
        solver_lev  (0:9)       =  ibuf(12:21)
        prec_type               =  ibuf(22)
        res_cvg_lev (0:9)       =  rbuf( 1:10)
        CFL_lev     (0:9)       =  rbuf(11:20)
        offset      (1:3)       =  rbuf(21:23)
        translation_axis(1:3)   =  rbuf(24:26)
        rotation_axis   (1:3)   =  rbuf(27:29)
        zoom                    =  rbuf(30)
        span                    =  rbuf(31:33)
        body_force(1:3)         =  rbuf(34:36)
    end if
    call norm_vec(3, translation_axis, rtmp)
    call norm_vec(3, rotation_axis   , rtmp)
    if(is_FV .or. is_heat)  is_include_corner   =  .true.
    is_hanging_node =  is_hanging_node .and. (is_FV .or. is_heat)
    is_structured_stencil   =  is_structured_stencil .and. is_FV
    if(span(1) .ge. -1.0d10)    call norm_vec(3, span, rtmp)
    prec_type   =  max(1, min(prec_type, 2))

    return
    end subroutine set_global
!-------------------------------------------------------------------------------
!   setup turbulence.
!-------------------------------------------------------------------------------
    subroutine set_turbulence
    use var_kind_def
    use var_global
    use var_parallel
    use var_slv, only: is_vis_cal
    use var_turb
    implicit none
    logical(dpL):: lbuf(1)
    integer(dpI):: io_err,ibuf(6)

    namelist /turbulence/   is_Spalart_QCR,LES_model,RANS_model,low_dissipation, &
        &  transition_model,DDES_model,length_scale
    if(.not. is_has_cfg)    return

    if(myid .eq. 0) then
        open(unit=10,file=trim(adjustl(cfg_file)))
        read(unit=10, nml=turbulence, iostat=io_err)
        if(io_err .gt. 0)   stop 'Error: fails to read namelist:turbulence.'
        close(10)

        lbuf(1:1)   = (/is_Spalart_QCR/)
        ibuf(1:6)   = (/RANS_model, LES_model, transition_model, DDES_model, &
                    &   low_dissipation, length_scale/)
    end if
    call mpi_bcast(lbuf, 1, mpi_dpL, 0, mpi_comm_world, mpi_err)
    call mpi_bcast(ibuf, 6, mpi_dpI, 0, mpi_comm_world, mpi_err)
    if(myid .ne. 0) then
        is_Spalart_QCR  =  lbuf(1)
        RANS_model      =  ibuf(1)
        LES_model       =  ibuf(2)
        transition_model=  ibuf(3)
        DDES_model      =  ibuf(4)
        low_dissipation =  ibuf(5)
        length_scale    =  ibuf(6)
    end if

    if(.not. is_vis_cal) then
        is_RANS         =  .false.
        is_DDES         =  .false.
        is_LES          =  .false.
        is_Spalart_QCR  =  .false.
        RANS_model      =  0
        DDES_model      =  0
        LES_model       =  0
        transition_model= -1
    elseif((RANS_model .eq. 1) .or. (RANS_model .eq. 2) .or. (RANS_model .eq. 11)) then
        is_RANS         =  .true.
        is_DDES         = (DDES_model .ge. 1) .and. (DDES_model .le. 3) .and. &
                        & (RANS_model .eq. 1)
        is_LES          =  .false.
        is_Spalart_QCR  =  is_Spalart_QCR .and. is_RANS .and. (.not. is_DDES)
        LES_model       =  0
        if(.not. is_DDES) then
            if(RANS_model .eq. 1) then
                if(transition_model .eq. transition_Menter) then
                    RANS_model  =  RANS_SAM
                elseif(transition_model .eq. transition_Coder) then
                    RANS_model  =  RANS_SAC
                end if
            elseif(RANS_model .eq. 2) then
                RANS_model  =  RANS_WA
            elseif(RANS_model .eq. RANS_SST) then
                RANS_model  =  RANS_SST
                is_KO       =  .true.
            end if
        end if
        if(.not. is_DDES)   DDES_model  =  0
        if(transition_model .lt. Transition_Menter) transition_model= -1
        if(transition_model .gt. Transition_BC    ) transition_model= -1
        if((transition_model .eq. transition_AGS) .and. (nprc .ne. 1)) &
            &  stop 'Error: AGS does not support parallel computation.'
    elseif((LES_model .ge. LES_WALE) .and. (LES_model .le. LES_sSMARG)) then
        is_RANS         =  .false.
        is_DDES         =  .false.
        is_LES          =  .true.
        is_Spalart_QCR  =  .false.
        RANS_model      =  0
        DDES_model      =  0
        transition_model= -1
    else
        is_RANS         =  .false.
        is_Spalart_QCR  =  .false.
        is_DDES         =  .false.
        is_LES          =  .false.
        RANS_model      =  0
        DDES_model      =  0
        LES_model       =  0
        transition_model= -1
    end if
    is_tur_cal      =  is_RANS .or. is_LES
    if(is_vis_cal) then
        low_dissipation =  min(2, max(0, low_dissipation))
    else
        low_dissipation =  0
    end if
    if(DDES_model .eq. 1) then
        n_DDES_quality  =  3
    elseif(DDES_model .eq. 2) then
        n_DDES_quality  =  3
    elseif(DDES_model .eq. 3) then
        n_DDES_quality  =  4
    end if

    return
    end subroutine set_turbulence
!-------------------------------------------------------------------------------
!   setup unsteady.
!-------------------------------------------------------------------------------
    subroutine set_unsteady
    use var_kind_def
    use var_global
    use var_parallel
    use var_uns_cal
    implicit none
    character(len=80):: str(3)
    integer(dpI):: io_err,ibuf(10)
    real   (dpR):: rbuf(10)

    namelist /unsteady/ max_uns_iter,dt_uns,tph_1,CFL_uns,uns_ra_begin,uns_method, &
                    &   max_subiter_uns,rhs_cvg_uns,uns_solution_1,uns_solution_2, &
                    &   uns_solution_3,min_subiter_uns,uns_write_delta
    if(.not. is_has_cfg)    return

    if(myid .eq. 0) then
        open(unit=10,file=trim(adjustl(cfg_file)))
        read(unit=10, nml=unsteady, iostat=io_err)
        if(io_err .gt. 0)   stop 'Error: fails to read namelist:unsteady.'
        close(10)
    end if
    str(1)      =  uns_solution_1
    str(2)      =  uns_solution_2
    str(3)      =  uns_solution_3
    ibuf(1:6)   = (/max_uns_iter, uns_ra_begin, uns_method, max_subiter_uns, &
                &   min_subiter_uns, uns_write_delta/)
    rbuf(1:4)   = (/dt_uns, tph_1, CFL_uns, rhs_cvg_uns/)
    call mpi_bcast(str , 240, mpi_character, 0, mpi_comm_world, mpi_err)
    call mpi_bcast(ibuf, 6  , mpi_dpI, 0, mpi_comm_world, mpi_err)
    call mpi_bcast(rbuf, 4  , mpi_dpR, 0, mpi_comm_world, mpi_err)
    if(myid .ne. 0) then
        uns_solution_1  =  str(1)
        uns_solution_2  =  str(2)
        uns_solution_3  =  str(3)
        max_uns_iter    =  ibuf(1)
        uns_ra_begin    =  ibuf(2)
        uns_method      =  ibuf(3)
        max_subiter_uns =  ibuf(4)
        min_subiter_uns =  ibuf(5)
        uns_write_delta =  ibuf(6)
        dt_uns          =  rbuf(1)
        tph_1           =  rbuf(2)
        CFL_uns         =  rbuf(3)
        rhs_cvg_uns     =  rbuf(4)
    end if

    max_uns_iter=  max(0, max_uns_iter)
    is_uns_cal  =  max_uns_iter .gt. 0
    is_uns_ra   =  is_uns_cal .and. (uns_ra_begin .le. max_uns_iter) .and. &
                & (uns_ra_begin .ge. 1)
    uns_method  =  max(1, min(4, uns_method))
    is_BDF2opt  =  uns_method .eq. uns_BDF2opt
    min_subiter_uns =  max(0, min(min_subiter_uns, max_subiter_uns))
    max_subiter_uns =  max(min_subiter_uns, max_subiter_uns)
    if(uns_method .eq. uns_RK4) then
        uns_n_stencil   =  2
    elseif(uns_method .eq. uns_BDF2) then
        uns_n_stencil   =  3
    elseif(uns_method .eq. uns_BDF2opt) then
        uns_n_stencil   =  4
    elseif(uns_method .eq. uns_IMEX) then
        uns_n_stencil   =  2
    else
        stop 'Error: unsteady method not supported.'
    end if

    return
    end subroutine set_unsteady
