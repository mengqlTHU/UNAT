!   THIS IS ONLY A SUITE OF RESEARCH CODE AND NEVER WORK AS COMMERCIAL SOFTWARE.
!   IT IS DEVELOPED BY DR. XINRONG SU AND IF YOU WERE ENTITLED TO USE THIS CODE,
!   REMEMBER TO CITE HIS CONTRIBUTIONS IN YOUR WORK.

!-------------------------------------------------------------------------------
!   module defines precisions of logical, integer and real.
!-------------------------------------------------------------------------------
    module var_kind_def
        implicit none

        integer(kind=4),parameter:: dpL =  kind(.true.)
        integer(kind=4),parameter:: dpI =  kind(1)
        integer(kind=4),parameter:: dpR =  kind(1.0d0)
    end module var_kind_def
!-------------------------------------------------------------------------------
!   module boundary values.
!-------------------------------------------------------------------------------
    module var_bndv
        use var_kind_def
        implicit none

        real   (dpR):: u_fs(5)  = (/1.0d0, 1.0d2, 0.0d0, 0.0d0, 1.0d5/)
        real   (dpR):: ua_fs    =  1.0d2
        real   (dpR):: a_fs     =  374.16574d0
        real   (dpR):: t_fs     =  2.8815d2
        real   (dpR):: mu_fs    =  1.5d-5
        real   (dpR):: Tu_fs    =  1.0d-2
        real   (dpR):: Re_fs    =  1.0d6
        real   (dpR):: Ma_fs    =  0.2d0
        real   (dpR):: Area_ref =  1.0d0
        real   (dpR):: AoA_fs(3)= (/0.0d0, 9.0d1, 9.0d1/)
        real   (dpR):: nut_nu   =  4.0d0

        real   (dpR):: flux_wall(10)= -1.0d99
        real   (dpR):: t_wall   (10)=  0.0d0

!       number of BND groups.
        integer(dpI):: n_inl            =  0
        integer(dpI):: n_out            =  0
        integer(dpI):: n_sol            =  0
        integer(dpI):: inl_ID(100)      =  1
        integer(dpI):: out_ID(100)      =  1
        integer(dpI):: sol_ID(100)      =  1

!       Inflow conditions.
        logical(dpL):: is_inl_cyl(10)   =  .false.
        logical(dpL):: is_inl_vis(10)   =  .false.
        real   (dpR):: tt_inl(10)       = -3.0d2
        real   (dpR):: pt_inl(10)       = -1.2d5
        real   (dpR):: vdir_inl(3,10)   = -1.0d99
        real   (dpR):: BL_thick_inl(10) = -1.0d0

!       Outflow condition.
        real   (dpR):: pb_out(10)   = -1.1d5

!       for BCInflowSupersonic.
        real   (dpR):: InflowSupersonic(5,10)   = -1.0d99

        real   (dpR):: mass_inl(10)             =  0.0d0
        real   (dpR):: mass_out(10)             =  0.0d0
    end module var_bndv
!-------------------------------------------------------------------------------
!   module CGNS.
!-------------------------------------------------------------------------------
    module var_cgns
        use var_kind_def
        implicit none

        include 'cgnslib_f.h'

        integer(kind=4),parameter:: BAR_5   =  27
        integer(kind=4),parameter:: QUAD_25 =  37

        integer(kind=4):: ifile,ibase,izone,iflow,ifield,icoord,nbases,isection, &
                        & ibocotype,iptset,npts,normalindex(3),normallistflag, &
                        & normaldatatype,ndataset,normallist,igr,cg_err
        integer(dpI):: nface_ele   (100)=  0
        integer(dpI):: npe_ele     (100)=  0
        integer(dpI):: npe_ele_1st (100)=  0
        integer(dpI):: ele_type_1st(100)=  0
    end module var_cgns
!-------------------------------------------------------------------------------
!   module constant.
!-------------------------------------------------------------------------------
    module var_global
        use var_kind_def
        implicit none

        real   (dpR),parameter:: pi     =  4.0d0*datan(1.0d0)
        real   (dpR),parameter:: R16    =  1.0d0/6.0d0
        real   (dpR),parameter:: R13    =  1.0d0/3.0d0
        real   (dpR),parameter:: R12    =  0.5d0
        real   (dpR),parameter:: R23    =  2.0d0/3.0d0
        real   (dpR),parameter:: R43    =  4.0d0/3.0d0
        real   (dpR),parameter:: I2(2,2)=reshape((/1.0d0,0.0d0,0.0d0,1.0d0/), (/2,2/))
        real   (dpR),parameter:: I5(5,5)= reshape( &
            & (/1.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, &
            &   0.0d0, 1.0d0, 0.0d0, 0.0d0, 0.0d0, &
            &   0.0d0, 0.0d0, 1.0d0, 0.0d0, 0.0d0, &
            &   0.0d0, 0.0d0, 0.0d0, 1.0d0, 0.0d0, &
            &   0.0d0, 0.0d0, 0.0d0, 0.0d0, 1.0d0/), (/5,5/))

        character(len=80):: cfg_file        = './data/mbt.cfg'
        character(len=80):: mesh_name       = './data/mesh'
        character(len=80):: mesh_file       = './data/mesh.cgns'
        character(len=80):: slv_file        = './data/mesh_slv.cgns'
        character(len=80):: res_file        = './data/mesh_res.dat'
        logical(dpL):: is_hanging_node      =  .false.
        logical(dpL):: is_structured_stencil=  .false.
        logical(dpL):: is_show_res          =  .true.
        logical(dpL):: is_has_cfg           =  .true.
        logical(dpL):: is_2d_cal            =  .true.
        logical(dpL):: is_flow_int          =  .false.
        logical(dpL):: is_debug             =  .false.
        logical(dpL):: is_output_rel_uvw    =  .false.
        logical(dpL):: is_output_cc         =  .false.
        logical(dpL):: is_output_vc         =  .true.
        logical(dpL):: is_output_cgns       =  .true.
        logical(dpL):: is_output_bnd        =  .true.
        logical(dpL):: is_output_cpcf       =  .false.
        integer(dpI):: n_dim                =  2
        integer(dpI):: err_mem              =  0

        real   (dpR):: offset(3)=  0.0d0
        real   (dpR):: zoom     =  1.0d0
        real   (dpR):: L_ref    =  1.0d0
        real   (dpR):: rref     =  1.0d0
        real   (dpR):: uref     =  1.0d2
        real   (dpR):: pref     =  1.0d5
        real   (dpR):: tm_ref   =  1.0d-2
        real   (dpR):: tref,mu_ref,Re_ref
        real   (dpR):: BL_thick =  1.0d-3

        real   (dpR):: mbt_t1   =  0.0d0
        real   (dpR):: mbt_t0   =  0.0d0
        real   (dpR):: translation_axis(3)  = (/0.0d0, 1.0d0, 0.0d0/)
        real   (dpR):: rotation_axis   (3)  = (/1.0d0, 0.0d0, 0.0d0/)
        real   (dpR):: span(1:3)            = -1.0d99
        real   (dpR):: body_force(1:3)      = -1.0d99

        logical(dpL):: sw_slave = .true.
        logical(dpL):: sw_time = .false.
    end module var_global
!-------------------------------------------------------------------------------
!   module global parameter with real type
!-------------------------------------------------------------------------------
    module var_global_real
        use var_kind_def
        implicit none

        real(dpR):: rhs_lev_r(9)
        real(dpR):: is_limiter_on_r,is_tur_cal_r,is_KO_r,is_vis_cal_r
        real(dpR):: is_RANS_r,RANS_model_r
    end module var_global_real
!-------------------------------------------------------------------------------
!   module jst_rk.
!-------------------------------------------------------------------------------
    module var_jst_rk
        use var_kind_def
        implicit none

        real(dpR),parameter:: RKa(5)    = (/0.25d0,1.0d0/6.0d0,0.375d0,0.5d0,1.0d0/)
        real(dpR),parameter:: RKgs2(2)  = (/0.5d0, 1.0d0/)
        real(dpR),parameter:: RKgs3(3)  = (/0.1048d0, 0.4d0, 1.0d0/)
        real(dpR),parameter:: RKgs5(5)  = (/6.95d-2, 1.602d-1, 2.898d-1, 5.06d-1, 1.0d0/)
    end module var_jst_rk
!-------------------------------------------------------------------------------
!   module parallel, both MPI and OPENMP.
!-------------------------------------------------------------------------------
    module var_parallel
        use var_kind_def
        implicit none

!       include 'mpif.f90'
        include 'mpif.h'

!       integer(dpI),parameter:: max_mpi_req=  30

        character(len=mpi_max_processor_name):: mpi_prc_name
        integer(kind=4):: mpi_prc_name_len
        integer(kind=4):: myid  =  0
        integer(kind=4):: nprc  =  1
        integer(kind=4):: mpi_group_world
        integer(kind=4):: psrc,pdst,mpi_err,mpi_errcode,mpi_nreq
        integer(kind=4):: mpi_status(mpi_status_size)
        integer(kind=4),allocatable:: mpi_req(:),mpi_sta(:,:)

        logical(dpL):: is_no_data_exchange  =  .true.
        logical(dpL):: is_include_corner    =  .false.
        integer(dpI),allocatable:: prc_info(:,:)
        integer:: mpi_dpL   =  31
        integer:: mpi_dpI   =  10
        integer:: mpi_dpR   =  15
    end module var_parallel
!-------------------------------------------------------------------------------
!   module solution. 
!-------------------------------------------------------------------------------
    module var_slv
        use var_kind_def
        implicit none

        integer(dpI),parameter:: solver_exp     =  1
        integer(dpI),parameter:: solver_rk_sgs  =  2
        integer(dpI),parameter:: solver_gmres_d =  3
        integer(dpI),parameter:: solver_gmres   =  4
        integer(dpI),parameter:: max_res_record =  1000

        logical(dpL):: mg_cycle_pre     =  .true.
        logical(dpL):: mg_cycle_post    =  .false.
        integer(dpI):: mg_cycle_vw      =  1

        logical(dpL):: is_reorder       =  .false.
        logical(dpL):: is_line_implicit =  .false.
        logical(dpL):: is_converged     =  .false.
        logical(dpL):: is_FV            =  .false.
        logical(dpL):: is_heat          =  .false.
        logical(dpL):: is_vis_cal       =  .false.
        logical(dpL):: is_local_dt      =  .true.
        logical(dpL):: is_st_cal_now    =  .true.
        logical(dpL):: is_matrix_free   =  .false.
        integer(dpI):: lev0_cal         =  0
        integer(dpI):: lev_out          =  0
        integer(dpI):: ite_wrk          =  0
        integer(dpI):: RK               =  0
        integer(dpI):: prec_type        =  1
        integer(dpI):: max_iter_lev(0:9)=  100
        integer(dpI):: solver_lev(0:9)  =  1
        integer(dpI):: tot_ite          =  0
        integer(dpI):: LDA_res_record   =  0
        integer(dpI):: n_res_record     =  0
        integer(dpI):: LS_weight_order  =  1
        integer(dpI):: gradient_method  =  1
        integer(dpI),parameter:: gradient_GG        =  1
        integer(dpI),parameter:: gradient_LS        =  2
        integer(dpI),parameter:: gradient_LSm       =  3
        integer(dpI),parameter:: LS_stencil         =  2
        integer(dpI),parameter:: LS_stencil_face    =  1
        integer(dpI),parameter:: LS_stencil_hybrid  =  2
        integer(dpI),parameter:: LS_stencil_vertex  =  3

        real   (dpR):: CFL_lev(0:9)     =  1.0d-1
        real   (dpR):: res_cvg_lev(0:9) =  1.0d-15
        real   (dpR):: dir_lift(3)      = (/0.0d0, 1.0d0, 0.0d0/)
        real   (dpR):: dir_drag(3)      = (/1.0d0, 0.0d0, 0.0d0/)
        real   (dpR):: mfin,mfot,drag_p,lift_p,drag_v,lift_v,drag,lift
        real   (dpR):: res_record(11,max_res_record)
        real   (dpR):: res_ref(5)   =  1.0d0

        integer(dpI),parameter:: max_monitor=  1000
        integer(dpI):: n_monitor    =  0
        integer(dpI):: n_monitor_L  =  0
        integer(dpI),allocatable:: n_monitor_prc(:)
        integer(dpI):: monitor_idx  (2,max_monitor) = -1
        integer(dpI):: monitor_idx_L(2,max_monitor) = -1
    end module var_slv
!-------------------------------------------------------------------------------
!   moduel temp.
!-------------------------------------------------------------------------------
    module var_temp
        use var_kind_def
        implicit none

        real   (dpR):: uLf(5,6),uRf(5,6),vgf(3,6)
        real   (dpR):: rhsl(5,6),rhsD(5,6),fac_1d(4,6),g_1d(12,6)
        real   (dpR):: u_structured(5,0:3)
        integer(dpI),allocatable:: ibf1(:)
    end module var_temp
!-------------------------------------------------------------------------------
!   moduel turbulence.
!-------------------------------------------------------------------------------
    module var_turb
        use var_kind_def
        implicit none

        logical(dpL):: is_DDES_now          =  .false.
        logical(dpL):: is_LES_now           =  .false.
        logical(dpL):: is_tur_cal           =  .false.
        logical(dpL):: is_RANS              =  .false.
        logical(dpL):: is_Spalart_QCR       =  .false.
        logical(dpL):: is_DDES              =  .false.
        logical(dpL):: is_get_DDES_quality  =  .false.
        logical(dpL):: is_KO                =  .false.
        integer(dpI):: RANS_model           =  0
        integer(dpI):: DDES_model           =  0
        integer(dpI):: n_DDES_quality       =  0
        integer(dpI):: low_dissipation      =  0
        integer(dpI),parameter:: RANS_SA    =  1
        integer(dpI),parameter:: RANS_WA    =  2
        integer(dpI),parameter:: RANS_SAM   =  3
        integer(dpI),parameter:: RANS_SAC   =  4
        integer(dpI),parameter:: RANS_SST   =  11
        integer(dpI),parameter:: SST_version=  1994

        integer(dpI),parameter:: DDES_Spalart   =  1
        integer(dpI),parameter:: DDES_WALE      =  2
        integer(dpI),parameter:: DDES_sigma     =  3
        integer(dpI),parameter:: DDES_IDDES     =  4

!       length scale.
        integer(dpI),parameter:: length_max     =  1
        integer(dpI),parameter:: length_Scotti  =  2
        integer(dpI),parameter:: length_Shur    =  3
        integer(dpI),parameter:: length_SSM     =  4
        integer(dpI),parameter:: length_LSQ     =  5
        integer(dpI),parameter:: length_Mockett =  6
        integer(dpI):: length_scale             =  1

!       LES variables.
        logical(dpL):: is_LES               =  .false.
        integer(dpI):: LES_model            =  0
        integer(dpI),parameter:: LES_WALE   =  1
        integer(dpI),parameter:: LES_sigma  =  2
        integer(dpI),parameter:: LES_VSS    =  3
        integer(dpI),parameter:: LES_sSMARG =  4
        real   (dpR),parameter:: cs_SMARG   =  0.18d0
        real   (dpR),parameter:: cw_WALE    =  sqrt(10.6d0)*cs_SMARG
        real   (dpR),parameter:: cr_VSS     =  1.3d0
        real   (dpR),parameter:: cr_sigma   =  1.5d0

        integer(dpI),parameter:: SA_version =  1
        integer(dpI),parameter:: QCR_version=  1
        real   (dpR),parameter:: kappa      =  0.41d0
        real   (dpR),parameter:: k2_1       =  5.94884d0
        real   (dpR),parameter:: cb1        =  0.1355d0
        real   (dpR),parameter:: cv1        =  7.1d0
        real   (dpR),parameter:: cw1        =  cb1*k2_1+2.433d0
        real   (dpR),parameter:: cw2        =  0.3d0
        real   (dpR),parameter:: Pr_SA      =  1.5d0
        real   (dpR),parameter:: C_DES      =  0.65d0
!       nut_ref =  sqrt(1.5)*U*Tu_fs*turbulent_length_scale
        real   (dpR):: nut_ref              =  5.0d-5
        real   (dpR):: k_ref                =  1.0d-2
        real   (dpR):: o_ref                =  5.0d4

        integer(dpI),parameter:: transition_Menter  =  1
        integer(dpI),parameter:: transition_Coder   =  2
        integer(dpI),parameter:: transition_PTM     =  3
        integer(dpI),parameter:: transition_AGS     =  4
        integer(dpI),parameter:: transition_BC      =  5
        integer(dpI):: transition_model             = -1
    end module var_turb
!-------------------------------------------------------------------------------
!   definition of unsteady computation.
!-------------------------------------------------------------------------------
    module var_uns_cal
        use var_kind_def
        implicit none

        integer(dpI),parameter:: uns_RK4    =  1
        integer(dpI),parameter:: uns_BDF2   =  2
        integer(dpI),parameter:: uns_BDF2opt=  3
        integer(dpI),parameter:: uns_IMEX   =  4

        real   (dpR),parameter:: BDF2opt_c1 =  1.66d0
        real   (dpR),parameter:: BDF2opt_c2 = -2.48d0
        real   (dpR),parameter:: BDF2opt_c3 =  0.98d0
        real   (dpR),parameter:: BDF2opt_c4 = -0.16d0

        logical(dpL):: is_uns_cal       =  .false.
        logical(dpL):: is_uns_cal_now   =  .false.
        logical(dpL):: is_bdf_now       =  .false.
        logical(dpL):: is_BDF2opt       =  .false.
        logical(dpL):: is_uns_ra        =  .false.
        integer(dpI):: max_uns_iter     =  0
        integer(dpI):: uns_iter         =  0
        integer(dpI):: uns_fail         =  0
        integer(dpI):: uns_ra_begin     =  huge(1)
        integer(dpI):: uns_method       =  3
        integer(dpI):: max_subiter_uns  =  20
        integer(dpI):: min_subiter_uns  =  15
        integer(dpI):: uns_write_delta  =  huge(1)/2
        real   (dpR):: CFL_uns          =  1.0d0
        real   (dpR):: dt_uns           =  0.0d0
        real   (dpR):: tph              =  0.0d0
        real   (dpR):: tph_1            =  0.0d0
        real   (dpR):: tph_0            =  0.0d0
        real   (dpR):: rhs_cvg_uns      =  1.0d-4
        real   (dpR):: LHS_uns_c1       =  1.5d0

        integer(dpI),parameter:: max_uns_record =  100
        integer(dpI):: n_uns_record             =  0
        real   (dpR),allocatable:: uns_record(:,:)

        character(len=80):: uns_solution_1  =  'NONE'
        character(len=80):: uns_solution_2  =  'NONE'
        character(len=80):: uns_solution_3  =  'NONE'
        logical(dpL):: is_uns_initialized   =  .false.
        integer(dpI):: uns_n_stencil        =  4
        integer(dpI):: uns_n_ff             =  1
    end module var_uns_cal
!-------------------------------------------------------------------------------
!   definition of IMEX method.
!-------------------------------------------------------------------------------
    module var_imex
        use var_kind_def
        implicit none

        logical(dpL):: is_IMEX_initialized  =  .false.
        integer(dpI):: IMEX_method  =  1
        integer(dpI):: n_stage      =  0
        integer(dpI):: i_stage      =  0
        real   (dpR):: A_im(4,4)    =  0.0d0
        real   (dpR):: b_im(4  )    =  0.0d0
        real   (dpR):: A_ex(4,4)    =  0.0d0
        real   (dpR):: b_ex(4  )    =  0.0d0

        contains
!       ------------------------------------------------------------------------
!       ARS(2,3,2) scheme.
!       ------------------------------------------------------------------------
        subroutine get_ars_232
        implicit none

        n_stage     =  3
        A_im(1,1)   =  0.0d0
        A_im(2,1)   =  0.0d0
        A_im(3,1)   =  0.0d0
        A_im(1,2)   =  0.0d0
        A_im(2,2)   =  1.0d0-0.5d0*sqrt(2.0d0)
        A_im(3,2)   =  1.0d0-A_im(2,2)
        A_im(1,3)   =  0.0d0
        A_im(2,3)   =  0.0d0
        A_im(3,3)   =  A_im(2,2)
        b_im(1)     =  0.0d0
        b_im(2)     =  A_im(3,2)
        b_im(3)     =  A_im(3,3)

        A_ex(1,1)   =  0.0d0
        A_ex(2,1)   =  1.0d0-0.5d0*sqrt(2.0d0)
        A_ex(3,1)   = -2.0d0*sqrt(2.0d0)/3.0d0
        A_ex(1,2)   =  0.0d0
        A_ex(2,2)   =  0.0d0
        A_ex(3,2)   =  1.0d0-A_ex(3,1)
        A_ex(1,3)   =  0.0d0
        A_ex(2,3)   =  0.0d0
        A_ex(3,3)   =  0.0d0
        b_ex(1)     =  0.0d0
        b_ex(2)     =  1.0d0-A_ex(2,1)
        b_ex(3)     =  A_ex(2,1)

        return
        end subroutine get_ars_232
!       ------------------------------------------------------------------------
!       ARS(2,2,2) scheme.
!       ------------------------------------------------------------------------
        subroutine get_ars_222
        implicit none
        real(kind=8):: gk,delta

        gk      =  1.0d0-0.5d0*sqrt(2.0d0)
        delta   =  1.0d0-0.5d0/gk

        n_stage     =  3
        A_im(1,1)   =  0.0d0
        A_im(2,1)   =  0.0d0
        A_im(3,1)   =  0.0d0
        A_im(1,2)   =  0.0d0
        A_im(2,2)   =  gk
        A_im(3,2)   =  1.0d0-gk
        A_im(1,3)   =  0.0d0
        A_im(2,3)   =  0.0d0
        A_im(3,3)   =  gk
        b_im(1)     =  0.0d0
        b_im(2)     =  1.0d0-gk
        b_im(3)     =  gk

        A_ex(1,1)   =  0.0d0
        A_ex(2,1)   =  gk
        A_ex(3,1)   =  delta
        A_ex(1,2)   =  0.0d0
        A_ex(2,2)   =  0.0d0
        A_ex(3,2)   =  1.0d0-delta
        A_ex(1,3)   =  0.0d0
        A_ex(2,3)   =  0.0d0
        A_ex(3,3)   =  0.0d0
        b_ex(1)     =  delta
        b_ex(2)     =  1.0d0-delta
        b_ex(3)     =  0.0d0

        return
        end subroutine get_ars_222
!       ------------------------------------------------------------------------
!       SSP2(2,2,2) scheme.
!       ------------------------------------------------------------------------
        subroutine get_ssp2_222
        implicit none
        real(kind=8):: gk

        gk      =  1.0d0-0.5d0*sqrt(2.0d0)

        n_stage     =  2
        A_im(1,1)   =  gk
        A_im(2,1)   =  1.0d0-2.0d0*gk
        A_im(1,2)   =  0.0d0
        A_im(2,2)   =  gk
        b_im(1)     =  0.5d0
        b_im(2)     =  0.5d0

        A_ex(1,1)   =  0.0d0
        A_ex(2,1)   =  1.0d0
        A_ex(1,2)   =  0.0d0
        A_ex(2,2)   =  0.0d0
        b_ex(1)     =  0.5d0
        b_ex(2)     =  0.5d0

        return
        end subroutine get_ssp2_222
    end module var_imex
