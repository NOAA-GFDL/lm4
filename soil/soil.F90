! ============================================================================
! soil model module
! ============================================================================
module soil_mod

#include "../shared/debug.inc"

#ifdef INTERNAL_FILE_NML
use mpp_mod, only: input_nml_file
#else
use fms_mod, only: open_namelist_file
#endif

use fms_mod, only: error_mesg, file_exist, check_nml_error, &
     stdlog, write_version_number, close_file, mpp_pe, mpp_root_pe, FATAL, WARNING, NOTE
use fms_io_mod, only: read_compressed, restart_file_type, free_restart_type, &
      field_exist, save_restart, &
      register_restart_field, set_domain, nullify_domain, get_field_size
use time_manager_mod,   only: time_type, increment_time, time_type_to_real
use diag_manager_mod,   only: diag_axis_init
use constants_mod,      only: tfreeze, hlv, hlf, dens_h2o
use tracer_manager_mod, only: NO_TRACER

use land_constants_mod, only : NBANDS, BAND_VIS, BAND_NIR, seconds_per_year
use soil_tile_mod, only : GW_LM2, GW_LINEAR, GW_HILL_AR5, GW_HILL, GW_TILED, &
     soil_tile_type, soil_pars_type, read_soil_data_namelist, &
     soil_data_radiation, soil_data_diffusion, soil_data_thermodynamics, &
     soil_data_hydraulic_properties, soil_data_psi_for_rh, &
     soil_data_gw_hydraulics, soil_data_gw_hydraulics_ar5, &
     soil_data_vwc_for_init_only, &
     soil_data_init_derive_subsurf_pars, &
     soil_data_init_derive_subsurf_pars_ar5, soil_data_init_derive_subsurf_pars_tiled, &
     max_lev, psi_wilt, cpw, clw, csw, g_iso, g_vol, g_geo, g_RT, aspect, &
     gw_scale_length, gw_scale_relief, gw_scale_soil_depth, &
     slope_exp, gw_scale_perm, k0_macro_x, retro_a0n1, &
     soil_type_file, &
     soil_tile_stock_pe, initval, comp, soil_theta, soil_ice_porosity

use soil_carbon_mod, only: poolTotals, get_pool_data_accessors, soilMaxCohorts, &
     update_pool, add_litter,add_C_N_to_rhizosphere, &
     tracer_leaching_with_litter,transfer_pool_fraction, n_c_types, &
     soil_carbon_option, SOILC_CENTURY, SOILC_CENTURY_BY_LAYER, SOILC_CORPSE, SOILC_CORPSE_N, &
     A_function, debug_pool, mycorrhizal_mineral_N_uptake_rate, mycorrhizal_decomposition, hypothetical_mycorrhizal_decomposition


use land_tile_mod, only : land_tile_type, land_tile_enum_type, &
     first_elmt, tail_elmt, next_elmt, prev_elmt, current_tile, get_elmt_indices, &
     operator(/=)
use land_utils_mod, only : put_to_tiles_r0d_fptr, put_to_tiles_r1d_fptr
use land_tile_diag_mod, only : diag_buff_type, &
     register_tiled_static_field, register_tiled_diag_field, &
     send_tile_data, send_tile_data_r0d_fptr, send_tile_data_r1d_fptr, &
     send_tile_data_i0d_fptr, &
     add_tiled_diag_field_alias, add_tiled_static_field_alias, &
     set_default_diag_filter
use land_data_mod, only : land_state_type, lnd, land_time
use land_io_mod, only : read_field
use land_tile_io_mod, only : create_tile_out_file, write_tile_data_r0d_fptr,&
     write_tile_data_r1d_fptr, read_tile_data_r0d_fptr, read_tile_data_r1d_fptr,&
     read_tile_data_layered_cohort_fptr, write_tile_data_layered_cohort_fptr,&
     write_tile_data_i1d_fptr_all,read_tile_data_i1d_fptr_all,&
     print_netcdf_error, get_input_restart_name, sync_nc_files, gather_tile_data, assemble_tiles
use nf_utils_mod, only : nfu_def_dim, nfu_put_att, nfu_inq_var
use vegn_data_mod, only: K1, K2, root_NH4_uptake_rate, root_NO3_uptake_rate
use vegn_tile_mod, only : vegn_tile_type, vegn_uptake_profile, vegn_hydraulic_properties
use land_debug_mod, only : is_watch_point, get_current_point, set_current_point, &
    check_var_range, check_conservation, is_watch_cell
use uptake_mod, only : UPTAKE_LINEAR, UPTAKE_DARCY2D, UPTAKE_DARCY2D_LIN, &
     uptake_init, darcy2d_uptake, darcy2d_uptake_solver

use hillslope_mod, only : do_hillslope_model, max_num_topo_hlsps, &
     num_vertclusters, hlsp_coldfracs, use_geohydrodata, & !pond, &
     horiz_wt_depth_to_init, calculate_wt_init, simple_inundation
use land_io_mod, only : &
     init_cover_field
use soil_tile_mod, only : n_dim_soil_types, soil_to_use, &
     soil_index_constant, input_cover_types
use hillslope_hydrology_mod, only: hlsp_hydro_lev_init, hlsp_hydrology_2, &
     stiff_explicit_gwupdate

! Test tridiagonal solution for advection
use land_numerics_mod, only : tridiag
implicit none
private

! ==== public interfaces =====================================================
public :: read_soil_namelist
public :: soil_cover_cold_start
public :: retrieve_soil_tags
public :: soil_init
public :: soil_end
public :: save_soil_restart
public :: save_soil_restart_new

public :: soil_get_sfc_temp
public :: soil_radiation
public :: soil_diffusion
public :: soil_step_1
public :: soil_step_2
public :: soil_step_3
public :: Dsdt
public :: add_root_litter
public :: add_root_exudates
public :: myc_scavenger_N_uptake
public :: hypothetical_myc_scavenger_N_uptake
public :: myc_miner_N_uptake
public :: hypothetical_myc_miner_N_uptake
public :: root_N_uptake
public :: redistribute_peat_carbon
! =====end of public interfaces ==============================================



! ==== module constants ======================================================
character(len=*), parameter, private   :: &
    module_name = 'soil',&
    version     = '$Id$',&
    tagname     = '$Name$'

! ==== module variables ======================================================

!---- namelist ---------------------------------------------------------------
logical :: lm2                  = .false.
logical :: use_E_min            = .false.     ! prohibit condensation
logical :: use_E_max            = .true.      ! theoretical effiltration capacity flag
real    :: init_temp            = 288.        ! cold-start soil T
real    :: init_w               = 150.        ! cold-start w(l)/dz(l)
real    :: init_wtdep           = -1.         ! positive value activates hydrostatic IC,
                                              ! overriding init_w
real    :: init_groundwater     =   0.        ! cold-start gw storage
real    :: lrunf_ie_min         = -1.0e-4     ! trigger for clip and runoff
real    :: lrunf_ie_tol         =  1.e-12
character(len=16) :: albedo_to_use = ''       ! or 'albedo-map' or 'brdf-maps'
character(len=24) :: uptake_to_use = 'linear' ! or 'darcy2d', or 'darcy2d-linearized'
logical :: uptake_oneway        = .false.     ! if true, roots cannot loose water to soil
logical :: uptake_from_sat      = .true.      ! if false, the uptake from saturated soil is prohibited
logical :: allow_negative_rie   = .false.
logical :: baseflow_where_frozen = .false.
logical :: write_when_flagged   = .false.
logical :: corrected_lm2_gw     = .true.
logical :: use_fringe           = .false.
logical :: push_down_sfc_excess = .true.
logical :: lrunf_from_div       = .true.
logical :: use_tall_plants      = .false.
logical :: cold_infilt          = .false.
logical :: bottom_up_cold_infilt= .false.
logical :: use_depth_to_wt_4    = .false.
logical :: always_use_bsw       = .false.
logical :: harmonic_mean_K      = .false.
logical :: verbose              = .false.
logical :: no_min_Dpsi          = .false.
logical :: div_bug_fix          = .false.
logical :: require_pos_hcap     = .false.
logical :: use_richards_clean   = .false.
real    :: active_layer_drainage_acceleration = 0.
real    :: hlf_factor           = 1.
real    :: gw_flux_max          = 1.e10
real    :: aquifer_heat_cap     = 0.         ! in equivalent liquid water amount, kg/m2
logical :: use_tridiag_foradvec = .false.    ! use tridiagonal solution for advection
                        ! information for soil carbon acceleration
logical :: write_soil_carbon_restart = .FALSE. ! indicates whether to write
                        ! information for soil carbon acceleration
logical :: horiz_init_wt        = .false.   ! initialize horizontal water table, if gw_option == GW_TILED
logical :: use_coldstart_wtt_data = .false. ! read additional data for soil initialization
character(len=256)  :: coldstart_datafile = 'INPUT/soil_wtt.nc'
logical :: allow_neg_rnu        = .false.   ! Refill from stream if wl < 0 with warning, i.e. during spinup.
logical :: allow_neg_wl         = .false.   ! Warn rather than abort if wl < 0, even if .not. allow_neg_rnu
logical :: prohibit_negative_water_div = .false. ! if TRUE, div_bf abd dif_if are
  ! set to zero in case water content of *any* layer is negative
real    :: zeta_bar_override    = -1.
real    :: cold_depth           = 0.
real    :: Wl_min               = -1.e20
real    :: bwood_macinf         = -1.
integer :: max_iter_trans = 100 ! max number of iterations for psi_crown_min
integer :: layer_for_gw_switch = 1000000 ! to accelerate permafrost gw shutoff
real    :: eps_trans     = 1.e-7 ! convergence crit for psi_crown_min
logical :: supercooled_rnu = .true. ! Excess ice converted to supercooled water for runoff.
real    :: wet_depth = 0.6 ! [m] water table depth threshold for diagnosing wetland fraction
real    :: thetathresh = -0.01 ! [-] threshold for negative soil liquid water volume
  ! before warning or abort
real    :: negrnuthresh = -0.1 ! [mm/s] threshold for negative lrunf_nu
  ! before warning or abort

real :: max_soil_C_density = 50.0   !(kgC/m3) -- for redistribution of peat
real :: max_litter_thickness = 0.05 ! m of litter layer thickness before it gets redistributed
real :: r_rhiz = 0.001              !Radius of rhizosphere around root (m)

namelist /soil_nml/ lm2, use_E_min, use_E_max,           &
                    init_temp,      &
                    init_w,   init_wtdep,    &
                    init_groundwater, lrunf_ie_min, lrunf_ie_tol, &
                    cpw, clw, csw, &
                    albedo_to_use, &
                    uptake_to_use, uptake_oneway, uptake_from_sat, &
                    allow_negative_rie, &
                    baseflow_where_frozen, &
                    write_when_flagged, &
                    corrected_lm2_gw, &
                    use_fringe, &
                    use_tall_plants, cold_infilt, bottom_up_cold_infilt, &
                    use_depth_to_wt_4, always_use_bsw, &
                    harmonic_mean_K, verbose, no_min_Dpsi, div_bug_fix, &
                    require_pos_hcap, use_richards_clean, &
                    push_down_sfc_excess, lrunf_from_div, &
                    active_layer_drainage_acceleration, hlf_factor, &
                    gw_flux_max, aquifer_heat_cap, use_tridiag_foradvec, &
                    horiz_init_wt, use_coldstart_wtt_data, coldstart_datafile, &
                    allow_neg_rnu, allow_neg_wl, prohibit_negative_water_div, &
                    zeta_bar_override, &
                    cold_depth, Wl_min, &
                    bwood_macinf, &
                    max_iter_trans, layer_for_gw_switch, eps_trans, &
                    supercooled_rnu, wet_depth, thetathresh, negrnuthresh, &
                    write_soil_carbon_restart, &
                    max_soil_C_density, max_litter_thickness, r_rhiz
!---- end of namelist --------------------------------------------------------

logical         :: module_is_initialized =.FALSE.
logical         :: use_brdf = .false.
real            :: delta_time ! fast (physical) time step, s
real            :: dt_fast_yr ! fast (physical) time step, yr (year is defined as 365 days)
logical         :: use_single_geo
integer         :: num_l              ! # of water layers
real            :: dz    (max_lev)    ! thicknesses of layers
real            :: zfull (max_lev)
real            :: zhalf (max_lev+1)
real            :: Eg_min

integer         :: uptake_option = -1
integer         :: gw_option = -1


real, allocatable :: soilCCohort_data(:)

! ---- diagnostic field IDs
integer :: id_fast_soil_C, id_slow_soil_C, id_protected_C, id_fsc, id_ssc,&
    id_fast_soil_N, id_slow_soil_N, id_protected_N, id_fsn, id_ssn,&
    id_leaflitter_deadmic_C, id_leaflitter_livemic_C, id_leaflitter_fast_C, id_leaflitter_slow_C, id_nleaflittercohorts, &
    id_leaflitter_deadmic_N, id_leaflitter_livemic_N, id_leaflitter_fast_N, id_leaflitter_slow_N, &
    id_finewoodlitter_deadmic_C, id_finewoodlitter_livemic_C, id_finewoodlitter_fast_C, id_finewoodlitter_slow_C, id_nfinewoodlittercohorts, &
    id_finewoodlitter_deadmic_N, id_finewoodlitter_livemic_N, id_finewoodlitter_fast_N, id_finewoodlitter_slow_N, &
    id_coarsewoodlitter_deadmic_C, id_coarsewoodlitter_livemic_C, id_coarsewoodlitter_fast_C, id_coarsewoodlitter_slow_C, id_ncoarsewoodlittercohorts, &
    id_coarsewoodlitter_deadmic_N, id_coarsewoodlitter_livemic_N, id_coarsewoodlitter_fast_N, id_coarsewoodlitter_slow_N, &
    id_deadmic_C, id_deadmic_N, id_livemic_C, id_livemic_N, id_nsoilcohorts, id_Qmax, id_protectedC, id_protectedN,  id_deadmic_C_total, id_deadmic_N_total, &
    id_livemic_C_total,id_livemic_N_total, &
    id_total_soil_C,id_dissolved_C_total,id_coarseWoodlitter_total_C,id_fineWoodlitter_total_C,id_leaflitter_total_C,&
    id_total_soil_N,id_dissolved_N_total,id_coarseWoodlitter_total_N,id_fineWoodlitter_total_N,id_leaflitter_total_N,&
    id_total_carbon_layered,id_total_nitrogen_layered,&
    id_lwc, id_swc, id_psi, id_temp, &
    id_ie, id_sn, id_bf, id_if, id_al, id_nu, id_sc, &
    id_hie, id_hsn, id_hbf, id_hif, id_hal, id_hnu, id_hsc, &
    id_heat_cap, id_thermal_cond, id_type, id_tau_gw, id_slope_l, &
    id_slope_Z, id_zeta_bar, id_e_depth, id_vwc_sat, id_vwc_fc, &
    id_vwc_wilt, id_K_sat, id_K_gw, id_w_fc, id_alpha, &
    id_refl_dry_dif, id_refl_dry_dir, id_refl_sat_dif, id_refl_sat_dir, &
    id_f_iso_dry, id_f_vol_dry, id_f_geo_dry, &
    id_f_iso_sat, id_f_vol_sat, id_f_geo_sat, &
    id_evap, id_uptk_n_iter, id_uptk, id_psi_x0, id_uptk_residual, &
    id_sws_n_iter, id_psi_x0_sws, &
    id_excess, id_deficit, id_deficit_2, id_deficit_3, id_deficit_4, id_zeta, id_tau, &
    id_psi_bot, id_sat_frac, id_stor_frac, id_sat_depth, id_sat_dept2, &
    id_cf_1, id_cf_3, id_wt_1, id_wt_2, id_wt_2a, id_wt_2b, id_wt_3, id_wt2_3, id_wt_4, &
    id_div_bf, id_div_if, id_div_al, &
    id_z_cap, id_active_layer, id_surface_water, id_inun_frac, id_rsn_frac, id_flow, id_reflux, &
    id_fast_C_leaching,id_slow_C_leaching,id_livemic_C_leaching,id_deadmic_C_leaching,&
    id_fast_N_leaching,id_slow_N_leaching,id_livemic_N_leaching,id_deadmic_N_leaching,&
    id_protected_C_leaching,id_protected_C_total,id_protected_N_leaching,id_protected_N_total,&
    id_nitrification_rate,id_denitrification_rate,id_leaflitter_nitrification_rate,&
    id_finewoodlitter_nitrification_rate, id_coarsewoodlitter_nitrification_rate,&
    id_leaflitter_denitrification_rate,&
    id_finewoodlitter_denitrification_rate, id_coarsewoodlitter_denitrification_rate,&
    id_N_immobilization_rate,id_N_mineralization_rate,id_leaflitter_N_immobilization_rate,id_leaflitter_N_mineralization_rate,&
    id_finewoodlitter_N_immobilization_rate, id_coarsewoodlitter_N_immobilization_rate,&
    id_finewoodlitter_N_mineralization_rate, id_coarsewoodlitter_N_mineralization_rate,&
    id_total_N_mineralization_rate, id_total_N_immobilization_rate, id_total_nitrification_rate, id_total_denitrification_rate,&
    id_fast_dissolved_C,id_slow_dissolved_C,id_deadmic_dissolved_C,&
    id_fast_dissolved_N,id_slow_dissolved_N,id_deadmic_dissolved_N,id_NO3,id_NH4,id_nitrif,id_denitrif,&
    id_leaflitter_fast_dissolved_C,id_leaflitter_slow_dissolved_C,id_leaflitter_deadmic_dissolved_C,&
    id_leaflitter_fast_dissolved_N,id_leaflitter_slow_dissolved_N,id_leaflitter_deadmic_dissolved_N,id_leaflitter_NO3,id_leaflitter_NH4,&
    id_finewoodlitter_fast_dissolved_C,id_finewoodlitter_slow_dissolved_C,id_finewoodlitter_deadmic_dissolved_C,&
    id_finewoodlitter_fast_dissolved_N,id_finewoodlitter_slow_dissolved_N,id_finewoodlitter_deadmic_dissolved_N,id_finewoodlitter_NO3,id_finewoodlitter_NH4,&
    id_coarsewoodlitter_fast_dissolved_C,id_coarsewoodlitter_slow_dissolved_C,id_coarsewoodlitter_deadmic_dissolved_C,&
    id_coarsewoodlitter_fast_dissolved_N,id_coarsewoodlitter_slow_dissolved_N,id_coarsewoodlitter_deadmic_dissolved_N,id_coarsewoodlitter_NO3,id_coarsewoodlitter_NH4,&
    id_leaflitter_nitrif,id_leaflitter_denitrif,id_finewoodlitter_nitrif,id_finewoodlitter_denitrif,id_coarsewoodlitter_nitrif,id_coarsewoodlitter_denitrif,&
    id_rsoil_leaflitter_deadmic, id_rsoil_leaflitter_fast, id_rsoil_leaflitter_slow, &
    id_rsoil_finewoodlitter_deadmic, id_rsoil_finewoodlitter_fast, id_rsoil_finewoodlitter_slow, &
    id_rsoil_coarsewoodlitter_deadmic, id_rsoil_coarsewoodlitter_fast, id_rsoil_coarsewoodlitter_slow, &
    id_rsoil_leaflitter_deadmic_N, id_rsoil_leaflitter_fast_N, id_rsoil_leaflitter_slow_N, &
    id_rsoil_finewoodlitter_deadmic_N, id_rsoil_finewoodlitter_fast_N, id_rsoil_finewoodlitter_slow_N, &
    id_rsoil_coarsewoodlitter_deadmic_N, id_rsoil_coarsewoodlitter_fast_N, id_rsoil_coarsewoodlitter_slow_N, &
    id_rsoil_fast, id_rsoil_slow, id_resp,id_rsoil_deadmic,id_asoil,id_rsoil,&
    id_rsoil_fast_N, id_rsoil_slow_N, id_rsoil_deadmic_N,id_rsoil_N,&
    id_leaflitter_C_dissolve_rate_fast,id_leaflitter_C_dissolve_rate_slow,id_leaflitter_C_dissolve_rate_deadmic,&
    id_leaflitter_C_deposition_fast,id_leaflitter_C_deposition_slow,id_leaflitter_C_deposition_deadmic,&
    id_leaflitter_N_dissolve_rate_fast,id_leaflitter_N_dissolve_rate_slow,id_leaflitter_N_dissolve_rate_deadmic,&
    id_leaflitter_N_deposition_fast,id_leaflitter_N_deposition_slow,id_leaflitter_N_deposition_deadmic,&
    id_finewoodlitter_C_dissolve_rate_fast,id_finewoodlitter_C_dissolve_rate_slow,id_finewoodlitter_C_dissolve_rate_deadmic,&
    id_finewoodlitter_N_deposition_fast,id_finewoodlitter_N_deposition_slow,id_finewoodlitter_N_deposition_deadmic,&
    id_finewoodlitter_N_dissolve_rate_fast,id_finewoodlitter_N_dissolve_rate_slow,id_finewoodlitter_N_dissolve_rate_deadmic,&
    id_finewoodlitter_C_deposition_fast,id_finewoodlitter_C_deposition_slow,id_finewoodlitter_C_deposition_deadmic,&
    id_coarsewoodlitter_C_dissolve_rate_fast,id_coarsewoodlitter_C_dissolve_rate_slow,id_coarsewoodlitter_C_dissolve_rate_deadmic,&
    id_coarsewoodlitter_C_deposition_fast,id_coarsewoodlitter_C_deposition_slow,id_coarsewoodlitter_C_deposition_deadmic,&
    id_coarsewoodlitter_N_dissolve_rate_fast,id_coarsewoodlitter_N_dissolve_rate_slow,id_coarsewoodlitter_N_dissolve_rate_deadmic,&
    id_coarsewoodlitter_N_deposition_fast,id_coarsewoodlitter_N_deposition_slow,id_coarsewoodlitter_N_deposition_deadmic,&
    id_C_dissolve_rate_fast,id_C_dissolve_rate_slow,id_C_dissolve_rate_deadmic,&
    id_N_dissolve_rate_fast,id_N_dissolve_rate_slow,id_N_dissolve_rate_deadmic,&
    id_C_deposition_fast,id_C_deposition_slow,id_C_deposition_deadmic, &
    id_N_deposition_fast,id_N_deposition_slow,id_N_deposition_deadmic, &
    id_leaflitter_fast_C_leaching,id_leaflitter_slow_C_leaching,id_leaflitter_deadmic_C_leaching,&
    id_finewoodlitter_fast_C_leaching,id_finewoodlitter_slow_C_leaching,id_finewoodlitter_deadmic_C_leaching,&
    id_coarsewoodlitter_fast_C_leaching,id_coarsewoodlitter_slow_C_leaching,id_coarsewoodlitter_deadmic_C_leaching,&
    id_fast_DOC_div_loss,id_slow_DOC_div_loss,id_deadmic_DOC_div_loss, &
    id_leaflitter_fast_N_leaching,id_leaflitter_slow_N_leaching,id_leaflitter_deadmic_N_leaching,&
    id_finewoodlitter_fast_N_leaching,id_finewoodlitter_slow_N_leaching,id_finewoodlitter_deadmic_N_leaching,&
    id_coarsewoodlitter_fast_N_leaching,id_coarsewoodlitter_slow_N_leaching,id_coarsewoodlitter_deadmic_N_leaching,&
    id_fast_DON_div_loss,id_slow_DON_div_loss,id_deadmic_DON_div_loss, &
    id_NO3_leaching, id_leaflitter_NO3_leaching, id_finewoodlitter_NO3_leaching, id_coarsewoodlitter_NO3_leaching,&
    id_NH4_leaching, id_leaflitter_NH4_leaching, id_finewoodlitter_NH4_leaching, id_coarsewoodlitter_NH4_leaching,&
    id_NO3_div_loss, &
    id_NH4_div_loss, &
    id_total_soil_carbon,id_total_soil_organic_N, id_total_soil_ammonium, id_total_soil_nitrate, id_wet_frac, id_macro_infilt, &
    id_surf_DOC_loss, id_total_C_leaching, id_total_DOC_div_loss, &
    id_surf_DON_loss, id_total_ON_leaching, id_total_DON_div_loss, &
    id_surf_NO3_loss, id_surf_NH4_loss, id_total_NO3_leaching, id_total_NH4_leaching, id_total_NO3_div_loss, id_total_NH4_div_loss,&
    id_rhizosphere_frac, id_soil_nitrate, id_soil_ammonium, id_leaflitter_nitrate, id_leaflitter_ammonium, id_finewoodlitter_nitrate, id_finewoodlitter_ammonium,&
    id_coarsewoodlitter_nitrate,id_coarsewoodlitter_ammonium, id_root_profile

    integer :: id_leaflitter_CUE_CO2,id_leaflitter_NUE_MINERAL,id_leaflitter_CUE_mic,id_leaflitter_NUE_mic,&
    			id_finewoodlitter_CUE_CO2,id_finewoodlitter_NUE_MINERAL,id_finewoodlitter_CUE_mic,id_finewoodlitter_NUE_mic,&
    			id_coarsewoodlitter_CUE_CO2,id_coarsewoodlitter_NUE_MINERAL,id_coarsewoodlitter_CUE_mic,id_coarsewoodlitter_NUE_mic,&
    			id_soil_CUE_CO2,id_soil_NUE_MINERAL,id_soil_CUE_mic,id_soil_NUE_mic

! test tridiagonal solver for advection
integer :: id_st_diff

integer, allocatable, dimension(:,:,:), private :: soil_tags ! module copy of soil tags for cold start

! ==== end of module variables ===============================================

! ==== NetCDF declarations ===================================================
include 'netcdf.inc'
#define __NF_ASRT__(x) call print_netcdf_error((x),module_name,__LINE__)

contains

! ============================================================================
subroutine read_soil_namelist()
  ! ---- local vars
  integer :: unit         ! unit for namelist i/o
  integer :: io           ! i/o status for the namelist
  integer :: ierr         ! error code, returned by i/o routines
  integer :: l

  call read_soil_data_namelist(num_l,dz,use_single_geo,gw_option)

  call write_version_number(version, tagname)
#ifdef INTERNAL_FILE_NML
  read (input_nml_file, nml=soil_nml, iostat=io)
  ierr = check_nml_error(io, 'soil_nml')
#else
  if (file_exist('input.nml')) then
     unit = open_namelist_file()
     ierr = 1;
     do while (ierr /= 0)
        read (unit, nml=soil_nml, iostat=io, end=10)
        ierr = check_nml_error (io, 'soil_nml')
     enddo
10   continue
     call close_file (unit)
  endif
#endif
  if (mpp_pe() == mpp_root_pe()) then
     unit=stdlog()
     write(unit, nml=soil_nml)
  endif

  ! ---- set up vertical discretization
  zhalf(1) = 0
  do l = 1, num_l;
     zhalf(l+1) = zhalf(l) + dz(l)
     zfull(l) = 0.5*(zhalf(l+1) + zhalf(l))
  enddo

  ! ---- convert symbolic names of options into numeric IDs to speed up
  ! the selection during run-time
  if (trim(uptake_to_use)=='linear') then
     uptake_option = UPTAKE_LINEAR
  else if (trim(uptake_to_use)=='darcy2d') then
     uptake_option = UPTAKE_DARCY2D
  else if (trim(uptake_to_use)=='darcy2d-linearized') then
     uptake_option = UPTAKE_DARCY2D_LIN
  else
     call error_mesg('soil_init',&
          'soil uptake option uptake_to_use="'//&
          trim(uptake_to_use)//'" is invalid, use "linear", "darcy2d" or "darcy2d-linearized"',&
          FATAL)
  endif

  if (use_E_min) then
     Eg_min = 0.
  else
     Eg_min = -HUGE(Eg_min)
  endif

  ! Configuration checking
  if (gw_option == GW_TILED .and. .not. use_tridiag_foradvec) then
     call error_mesg(module_name, 'read_soil_namelist: over-riding "use_tridiag_foradvec" value set '// &
                     'in namelist or module default. Tiled groundwater model / full hillslope model requires '// &
                     '"use_tridiag_foradvec" == .true.', NOTE)
     use_tridiag_foradvec = .true.
  end if

end subroutine read_soil_namelist


! ============================================================================
! initialize soil model
subroutine soil_init ( id_lon, id_lat, id_band, id_zfull, new_land_io)
  integer, intent(in)  :: id_lon  ! ID of land longitude (X) axis
  integer, intent(in)  :: id_lat  ! ID of land latitude (Y) axis
  integer, intent(in)  :: id_band ! ID of spectral band axis
  integer, intent(out) :: id_zfull ! ID of vertical soil axis
  logical, intent(in) :: new_land_io

  ! ---- local vars
  integer :: unit, unit1  ! unit numbers for various i/o
  type(land_tile_enum_type)     :: te,ce  ! tail and current tile list elements
  type(land_tile_type), pointer :: tile   ! pointer to current tile
  ! input data buffers for respective variables:
  real, allocatable :: gw_param(:,:), gw_param2(:,:), gw_param3(:,:), albedo(:,:,:)
  real, allocatable :: f_iso(:,:,:), f_vol(:,:,:), f_geo(:,:,:), refl_dif(:,:,:)

  real :: local_wt_depth ! [m] water table depth for tile (+ for below surface)
  real, allocatable :: ref_soil_t(:,:) ! reference soil temperature (based on 5 m or surface air temperature)
                                       ! for cold-start initialization
  real, allocatable :: wetmask(:,:)    ! input mask for zones with high water table
  logical :: drypoint                  ! This point is predicted to have a falling water table.

  integer :: i, li, lj, k ! indices
  real :: psi(num_l), mwc(num_l)
  character(len=256) :: restart_file_name
  character(len=17)  :: restart_base_name='INPUT/soil.res.nc'
  integer :: siz(4), isize,  nz, ncc, nccoh
  integer, allocatable :: idx(:)  ! I/O domain vector of compressed indices
  real,    allocatable :: r2d(:,:,:), r1d(:,:), r1dc(:,:), r0d(:) ! I/O domain level dependent vector of real data
  integer, allocatable :: i1d(:,:)  ! I/O domain level dependent vector of integer data
  logical :: restart_exists, found

  module_is_initialized = .TRUE.
  delta_time = time_type_to_real(lnd%dt_fast)
  dt_fast_yr = delta_time/seconds_per_year

  allocate(soilCCohort_data(soilMaxCohorts))
  soilCCohort_data = (/(nccoh,nccoh=1,soilMaxCohorts)/)

  call uptake_init(num_l,dz,zfull)
  call hlsp_hydro_lev_init(num_l,dz,zfull)



  ! -------- initialize soil model diagnostic fields
  call soil_diag_init ( id_lon, id_lat, id_band, id_zfull)

  ! -------- read spatially distributed fields for groundwater parameters, if requested
  if (.not.use_single_geo) then
     select case (gw_option)
     case (GW_LINEAR,GW_LM2)
        allocate(gw_param(lnd%is:lnd%ie,lnd%js:lnd%je))
        call read_field( 'INPUT/groundwater_residence.nc','tau', lnd%lon, lnd%lat, &
             gw_param, interp='bilinear' )
        call put_to_tiles_r0d_fptr( gw_param, lnd%tile_map, soil_tau_groundwater_ptr )
        deallocate(gw_param)
     case (GW_HILL, GW_HILL_AR5)
        allocate(gw_param (lnd%is:lnd%ie,lnd%js:lnd%je))
        allocate(gw_param2(lnd%is:lnd%ie,lnd%js:lnd%je))
        allocate(gw_param3(lnd%is:lnd%ie,lnd%js:lnd%je))
        call read_field( 'INPUT/geohydrology.nc','hillslope_length',  lnd%lon, lnd%lat, &
          gw_param, interp='bilinear' )
        call put_to_tiles_r0d_fptr( gw_param*gw_scale_length, lnd%tile_map, soil_hillslope_length_ptr )
        call read_field( 'INPUT/geohydrology.nc','slope', lnd%lon, lnd%lat, &
          gw_param2, interp='bilinear' )
        gw_param = gw_param*gw_param2
        call put_to_tiles_r0d_fptr( gw_param*gw_scale_relief, lnd%tile_map, soil_hillslope_relief_ptr )

        if (retro_a0n1 .or. gw_option.eq.GW_HILL_AR5) then
            gw_param = 0.
            call put_to_tiles_r0d_fptr( gw_param, lnd%tile_map, soil_hillslope_a_ptr )
            gw_param = 1.
            call put_to_tiles_r0d_fptr( gw_param, lnd%tile_map, soil_hillslope_n_ptr )
!            call read_field( 'INPUT/geohydrology.nc','hillslope_zeta_bar', &
!              lnd%lon, lnd%lat, gw_param, interp='bilinear' )
            gw_param = 0.5
            call put_to_tiles_r0d_fptr( gw_param, lnd%tile_map, soil_hillslope_zeta_bar_ptr )
        else
            call read_field( 'INPUT/geohydrology.nc','hillslope_a', &
              lnd%lon, lnd%lat, gw_param, interp='bilinear' )
            call put_to_tiles_r0d_fptr( gw_param, lnd%tile_map, soil_hillslope_a_ptr )
            call read_field( 'INPUT/geohydrology.nc','hillslope_n', &
              lnd%lon, lnd%lat, gw_param2, interp='bilinear' )
            call put_to_tiles_r0d_fptr( gw_param2, lnd%tile_map, soil_hillslope_n_ptr )
            gw_param3 = (1./(gw_param2+1.)+gw_param/(gw_param2+2.))/(1.+gw_param/2.)
            call put_to_tiles_r0d_fptr( gw_param3, lnd%tile_map, soil_hillslope_zeta_bar_ptr )
        endif

        call read_field( 'INPUT/geohydrology.nc','soil_e_depth', &
          lnd%lon, lnd%lat, gw_param, interp='bilinear' )
        if (slope_exp.gt.0.01) then
            call put_to_tiles_r0d_fptr( gw_param*gw_scale_soil_depth*(0.08/gw_param2)**slope_exp, &
                                                  lnd%tile_map, soil_soil_e_depth_ptr )
        else
            call put_to_tiles_r0d_fptr( gw_param*gw_scale_soil_depth, lnd%tile_map, soil_soil_e_depth_ptr )
        endif
        if (gw_option /= GW_HILL_AR5) then
            call read_field( 'INPUT/geohydrology.nc','perm', lnd%lon, lnd%lat, &
                 gw_param, interp='bilinear' )
            call put_to_tiles_r0d_fptr(9.8e9*gw_scale_perm*gw_param, lnd%tile_map, &
                                            soil_k_sat_gw_ptr )
        endif
        deallocate(gw_param, gw_param2, gw_param3)
        te = tail_elmt (lnd%tile_map)
        ce = first_elmt(lnd%tile_map)
        do while(ce /= te)
            tile=>current_tile(ce)  ! get pointer to current tile
            ce=next_elmt(ce)        ! advance position to the next tile
            if (.not.associated(tile%soil)) cycle
            select case (gw_option)
            case (GW_HILL)
                call soil_data_init_derive_subsurf_pars(tile%soil)
            case (GW_HILL_AR5)
                call soil_data_init_derive_subsurf_pars_ar5(tile%soil)
            end select
        enddo
     case (GW_TILED)
        if (use_geohydrodata) then
           allocate(gw_param (lnd%is:lnd%ie,lnd%js:lnd%je))
           allocate(gw_param2(lnd%is:lnd%ie,lnd%js:lnd%je))
           call read_field( 'INPUT/geohydrology.nc','hillslope_length',  lnd%lon, lnd%lat, &
             gw_param, interp='bilinear' )
           call put_to_tiles_r0d_fptr( gw_param*gw_scale_length, lnd%tile_map, soil_hillslope_length_ptr )
           call read_field( 'INPUT/geohydrology.nc','slope', lnd%lon, lnd%lat, &
             gw_param2, interp='bilinear' )
           gw_param = gw_param*gw_param2
           call put_to_tiles_r0d_fptr( gw_param*gw_scale_relief, lnd%tile_map, soil_hillslope_relief_ptr )
           call read_field( 'INPUT/geohydrology.nc','hillslope_zeta_bar', &
             lnd%lon, lnd%lat, gw_param, interp='bilinear' )
           if (zeta_bar_override.gt.0.) gw_param=zeta_bar_override
           call put_to_tiles_r0d_fptr( gw_param, lnd%tile_map, soil_hillslope_zeta_bar_ptr )
           call read_field( 'INPUT/geohydrology.nc','soil_e_depth', &
             lnd%lon, lnd%lat, gw_param, interp='bilinear' )

           if (slope_exp.gt.0.01) then
           ! ZMS It's probably inconsistent to leave in this if statement.
               call error_mesg(module_name, 'soil_init: "slope_exp" > 0.0 requested even though '// &
                               'running with tiled hillslope model.  This may be inconsistent.', WARNING)
               call put_to_tiles_r0d_fptr( gw_param*gw_scale_soil_depth*(0.08/gw_param2)**slope_exp, &
                                                     lnd%tile_map, soil_soil_e_depth_ptr )
           else
               call put_to_tiles_r0d_fptr( gw_param*gw_scale_soil_depth, lnd%tile_map, soil_soil_e_depth_ptr )
           endif
           call read_field( 'INPUT/geohydrology.nc','perm', lnd%lon, lnd%lat, &
                  gw_param, interp='bilinear' )
           call put_to_tiles_r0d_fptr(9.8e9*gw_scale_perm*gw_param, lnd%tile_map, &
                                          soil_k_sat_gw_ptr )
           deallocate(gw_param, gw_param2)
        end if
        te = tail_elmt (lnd%tile_map)
        ce = first_elmt(lnd%tile_map)
        do while(ce /= te)
            tile=>current_tile(ce)  ! get pointer to current tile
            ce=next_elmt(ce)        ! advance position to the next tile
            if (.not.associated(tile%soil)) cycle
            call soil_data_init_derive_subsurf_pars_tiled(tile%soil, use_geohydrodata)
        end do
     end select ! gw_option
  else if (gw_option == GW_TILED) then ! and use_single_geo
     ! Error checking
     if (.not. use_geohydrodata) then
        call error_mesg(module_name, 'soil_init: incompatible namelist options selected. gw_option =='// &
                        ' tiled, use_geohydrodata == .false., and use_single_geo == .true.', FATAL)
     else
        call error_mesg(module_name, 'soil_init: Warning: using tiled hillslope groundwater model '// &
                        'with single global values for soil hydrological properties (i.e. "use_single_geo).', &
                        WARNING)
     end if
  endif ! single geo

  ! -------- set dry soil albedo values, if requested
  if (trim(albedo_to_use)=='albedo-map') then
     allocate(albedo(lnd%is:lnd%ie,lnd%js:lnd%je,NBANDS))
     call read_field( 'INPUT/soil_albedo.nc','SOIL_ALBEDO_VIS',&
          lnd%lon, lnd%lat, albedo(:,:,BAND_VIS),'bilinear')
     call read_field( 'INPUT/soil_albedo.nc','SOIL_ALBEDO_NIR',&
          lnd%lon, lnd%lat, albedo(:,:,BAND_NIR),'bilinear')
     call put_to_tiles_r1d_fptr( albedo, lnd%tile_map, soil_refl_dry_dir_ptr )
     call put_to_tiles_r1d_fptr( albedo, lnd%tile_map, soil_refl_dry_dif_ptr )
     ! for now, put the same value into the saturated soil albedo, so that
     ! the albedo does not depend on soil wetness
     call put_to_tiles_r1d_fptr( albedo, lnd%tile_map, soil_refl_sat_dir_ptr )
     call put_to_tiles_r1d_fptr( albedo, lnd%tile_map, soil_refl_sat_dif_ptr )
     deallocate(albedo)
  else if (trim(albedo_to_use)=='brdf-maps') then
     use_brdf = .true.
     allocate(   f_iso(lnd%is:lnd%ie,lnd%js:lnd%je,NBANDS))
     allocate(   f_vol(lnd%is:lnd%ie,lnd%js:lnd%je,NBANDS))
     allocate(   f_geo(lnd%is:lnd%ie,lnd%js:lnd%je,NBANDS))
     allocate(refl_dif(lnd%is:lnd%ie,lnd%js:lnd%je,NBANDS))
     call read_field( 'INPUT/soil_brdf.nc','f_iso_vis',&
          lnd%lon, lnd%lat, f_iso(:,:,BAND_VIS),'bilinear')
     call read_field( 'INPUT/soil_brdf.nc','f_vol_vis',&
          lnd%lon, lnd%lat, f_vol(:,:,BAND_VIS),'bilinear')
     call read_field( 'INPUT/soil_brdf.nc','f_geo_vis',&
          lnd%lon, lnd%lat, f_geo(:,:,BAND_VIS),'bilinear')
     call read_field( 'INPUT/soil_brdf.nc','f_iso_nir',&
          lnd%lon, lnd%lat, f_iso(:,:,BAND_NIR),'bilinear')
     call read_field( 'INPUT/soil_brdf.nc','f_vol_nir',&
          lnd%lon, lnd%lat, f_vol(:,:,BAND_NIR),'bilinear')
     call read_field( 'INPUT/soil_brdf.nc','f_geo_nir',&
          lnd%lon, lnd%lat, f_geo(:,:,BAND_NIR),'bilinear')
     refl_dif = g_iso*f_iso + g_vol*f_vol + g_geo*f_geo
     call put_to_tiles_r1d_fptr( f_iso,    lnd%tile_map, soil_f_iso_dry_ptr )
     call put_to_tiles_r1d_fptr( f_vol,    lnd%tile_map, soil_f_vol_dry_ptr )
     call put_to_tiles_r1d_fptr( f_geo,    lnd%tile_map, soil_f_geo_dry_ptr )
     call put_to_tiles_r1d_fptr( refl_dif, lnd%tile_map, soil_refl_dry_dif_ptr )
     ! for now, put the same value into the saturated soil albedo, so that
     ! the albedo does not depend on soil wetness
     call put_to_tiles_r1d_fptr( f_iso,    lnd%tile_map, soil_f_iso_sat_ptr )
     call put_to_tiles_r1d_fptr( f_vol,    lnd%tile_map, soil_f_vol_sat_ptr )
     call put_to_tiles_r1d_fptr( f_geo,    lnd%tile_map, soil_f_geo_sat_ptr )
     call put_to_tiles_r1d_fptr( refl_dif, lnd%tile_map, soil_refl_sat_dif_ptr )
     deallocate(f_iso, f_vol, f_geo, refl_dif)
  else if (trim(albedo_to_use)=='') then
     ! do nothing, that is leave soil albedo parameters as defined based on the data table
  else
     call error_mesg('soil_init',&
          'option albedo_to_use="'// trim(albedo_to_use)//&
          '" is invalid, use "albedo-map", "brdf-maps", or empty line ("")',&
          FATAL)
  endif

  ! Call calculate_wt_init outside tile loop so that it is done once per hillslope
  if (init_wtdep .gt. 0. .and. gw_option == GW_TILED) then
     call calculate_wt_init(init_wtdep)
  end if

  if (use_coldstart_wtt_data) then
     allocate(ref_soil_t(lnd%is:lnd%ie,lnd%js:lnd%je), wetmask(lnd%is:lnd%ie,lnd%js:lnd%je))
     call read_field( coldstart_datafile, 'REFSOILT', &
             lnd%lon, lnd%lat, ref_soil_t, interp='bilinear' )
     call read_field( coldstart_datafile, 'WETMASK', &
             lnd%lon, lnd%lat, wetmask, interp='bilinear' )
  end if

  ! -------- initialize soil state --------
  te = tail_elmt (lnd%tile_map)
  ce = first_elmt(lnd%tile_map, is=lnd%is, js=lnd%js) ! Use global indices here because element indices
                                                      ! needed.
  do while(ce /= te)
     tile=>current_tile(ce)  ! get pointer to current tile
     ce=next_elmt(ce)        ! advance position to the next tile
     if (.not.associated(tile%soil)) cycle
     ! Retrieve indices
     call get_elmt_indices(prev_elmt(ce),i=li,j=lj,k=k)
     call set_current_point(li,lj,k)
     if (init_wtdep .gt. 0.) then
        if (.not. use_coldstart_wtt_data) then
           if (horiz_init_wt .and. gw_option == GW_TILED) then
              call horiz_wt_depth_to_init(tile%soil, prev_elmt(ce), local_wt_depth)
              ! Note: if restart_exists, then this function returns dummy local_wt_depth == 0.
              ! prev_elmt(ce) passed because indices will be needed.
              psi = zfull(1:num_l) - local_wt_depth
           else
              psi = zfull(1:num_l) - init_wtdep
           end if
        else
           if (wetmask(li, lj) > 0.5) then ! wet point
              drypoint = .false.
           else
              drypoint = .true.
           end if
           if (gw_option == GW_TILED) then
              call horiz_wt_depth_to_init(tile%soil, prev_elmt(ce), local_wt_depth, dry=drypoint)
              psi = zfull(1:num_l) - local_wt_depth
           else if (drypoint) then
              psi = zfull(1:num_l) - tile%soil%pars%hillslope_relief*tile%soil%pars%hillslope_zeta_bar
           else
              psi = zfull(1:num_l) - init_wtdep
           end if
        end if
        call soil_data_vwc_for_init_only(tile%soil, psi, mwc)
        mwc = mwc * dens_h2o
     else if (init_w .ge. 0.) then
        mwc = init_w
     else ! negative init_w is to be intrepreted as prescribed saturation
        mwc = -init_w*tile%soil%pars%vwc_sat*dens_h2o
     endif
     if (.not. use_coldstart_wtt_data) then
        if (init_temp.ge.tile%soil%pars%tfreeze) then
           tile%soil%wl = mwc*dz(1:num_l)
           tile%soil%ws = 0
        else
           tile%soil%wl = 0
           tile%soil%ws = mwc*dz(1:num_l)
        endif
        tile%soil%T             = init_temp
        tile%soil%groundwater   = init_groundwater
        tile%soil%groundwater_T = init_temp
        tile%soil%uptake_T           = init_temp
     else
        call init_soil_twc(tile%soil, ref_soil_t(li, lj), mwc)
     end if
  end do

  if (use_coldstart_wtt_data) then
     deallocate(ref_soil_t, wetmask)
  end if

  call get_input_restart_name(restart_base_name,restart_exists,restart_file_name,new_land_io)
  if (restart_exists) then
     if(new_land_io) then
        restart_file_name = restart_base_name

        call error_mesg('soil_init', 'Using new soil restart read', NOTE)
        call get_field_size(restart_file_name, 'tile_index', siz, field_found=found, domain=lnd%domain)
        if ( .not.found ) call error_mesg(trim(module_name), &
             'tile_index axis not found in '//trim(restart_file_name), FATAL)
        isize = siz(1)

        call get_field_size(restart_file_name, 'zfull', siz, field_found=found, domain=lnd%domain)
        if ( .not.found ) call error_mesg(trim(module_name), 'Z axis not found in '//trim(restart_file_name), FATAL)
        nz = siz(1)

        allocate(idx(isize),r1d(isize,nz),r0d(isize),i1d(isize,nz))
        call read_compressed(restart_file_name,'tile_index',idx, domain=lnd%domain, timelevel=1)

        call read_compressed(restart_file_name,'temp',r1d, domain=lnd%domain, timelevel=1)
        call assemble_tiles(soil_T_ptr,idx,r1d)

        call read_compressed(restart_file_name,'wl',r1d, domain=lnd%domain, timelevel=1)
        call assemble_tiles(soil_wl_ptr,idx,r1d)

        call read_compressed(restart_file_name,'ws',r1d, domain=lnd%domain, timelevel=1)
        call assemble_tiles(soil_ws_ptr,idx,r1d)

        call read_compressed(restart_file_name,'groundwater',r1d, domain=lnd%domain, timelevel=1)
        call assemble_tiles(soil_groundwater_ptr,idx,r1d)

        call read_compressed(restart_file_name,'groundwater_T',r1d, domain=lnd%domain, timelevel=1)
        call assemble_tiles(soil_groundwater_T_ptr,idx,r1d)

        if ( field_exist(restart_file_name,'uptake_T', domain=lnd%domain) ) then
           call read_compressed(restart_file_name,'uptake_T',r0d, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soil_uptake_T_ptr,idx,r0d)
        endif

        if (soil_carbon_option==SOILC_CENTURY.or.soil_carbon_option==SOILC_CENTURY_BY_LAYER) then
           if ( field_exist(restart_file_name,'fsc', domain=lnd%domain) ) then
              call read_compressed(restart_file_name,'fsc',r1d, domain=lnd%domain, timelevel=1)
              call assemble_tiles(soil_fast_soil_C_ptr,idx,r1d)
              call read_compressed(restart_file_name,'ssc',r1d, domain=lnd%domain, timelevel=1)
              call assemble_tiles(soil_slow_soil_C_ptr,idx,r1d)
           else
              ! try to read fsc and ssc from vegetation restart
              call get_input_restart_name('INPUT/vegn2.res.nc',restart_exists,restart_file_name,new_land_io)
              if (restart_exists) then
                restart_file_name = 'INPUT/vegn2.res.nc'
                call read_compressed(restart_file_name,'fsc',r0d, domain=lnd%domain, timelevel=1)
                call assemble_tiles(soil_fast_soil_C_ptr,idx,r0d,1)
                call read_compressed(restart_file_name,'ssc',r0d, domain=lnd%domain, timelevel=1)
                call assemble_tiles(soil_slow_soil_C_ptr,idx,r0d,1)
              endif
           endif
        endif

        if ( field_exist(restart_file_name,'fast_soil_C', domain=lnd%domain) ) then
           call get_field_size(restart_file_name, 'soilCCohort', siz, field_found=found, domain=lnd%domain)
           if ( .not.found ) call error_mesg(trim(module_name), &
                'soil carbon cohort axis not found in '//trim(restart_file_name), FATAL)
           ncc = siz(1)
           allocate(r2d(isize,nz,ncc),r1dc(isize,ncc))
           call read_compressed(restart_file_name,'fast_soil_C',r2d, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soilc_fast_soil_C_ptr,idx,r2d)
           call read_compressed(restart_file_name,'slow_soil_C',r2d, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soilc_slow_soil_C_ptr,idx,r2d)
           call read_compressed(restart_file_name,'deadMic',r2d, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soilc_deadMicrobeC_ptr,idx,r2d)
           call read_compressed(restart_file_name,'fastProtectedC',r2d, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soilc_fast_protected_C_ptr,idx,r2d)
           call read_compressed(restart_file_name,'slowProtectedC',r2d, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soilc_slow_protected_C_ptr,idx,r2d)
           call read_compressed(restart_file_name,'deadMicrobeProtectedC',r2d, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soilc_deadMicrobe_protected_C_ptr,idx,r2d)
           call read_compressed(restart_file_name,'liveMic',r2d, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soilc_livingMicrobeC_ptr,idx,r2d)
           call read_compressed(restart_file_name,'CO2',r2d, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soilc_CO2_ptr,idx,r2d)
           call read_compressed(restart_file_name,'Rtot',r2d, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soilc_Rtot_ptr,idx,r2d)
           call read_compressed(restart_file_name,'originalCohortC',r2d, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soilc_originalLitterC_ptr,idx,r2d)

if(soil_carbon_option == SOILC_CORPSE_N) then
           call read_compressed(restart_file_name,'fast_soil_N',r2d, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soilc_fast_soil_N_ptr,idx,r2d)
           call read_compressed(restart_file_name,'slow_soil_N',r2d, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soilc_slow_soil_N_ptr,idx,r2d)
           call read_compressed(restart_file_name,'deadMic_N',r2d, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soilc_deadMicrobeN_ptr,idx,r2d)
           call read_compressed(restart_file_name,'fastProtectedN',r2d, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soilc_fast_protected_N_ptr,idx,r2d)
           call read_compressed(restart_file_name,'slowProtectedN',r2d, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soilc_slow_protected_N_ptr,idx,r2d)
           call read_compressed(restart_file_name,'deadMicrobeProtectedN',r2d, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soilc_deadMicrobe_protected_N_ptr,idx,r2d)
           call read_compressed(restart_file_name,'liveMic_N',r2d, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soilc_livingMicrobeN_ptr,idx,r2d)
           call read_compressed(restart_file_name,'originalCohortN',r2d, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soilc_originalLitterC_ptr,idx,r2d)
           ! Leaving out IMM_N_max, IMM_N_gross, MINER_gross, etc for now -- BNS

           call read_compressed(restart_file_name,'soil_DON_fast',r1d, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soil_fast_DON_ptr,idx,r1d)
           call read_compressed(restart_file_name,'soil_DON_slow',r1d, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soil_slow_DON_ptr,idx,r1d)
           call read_compressed(restart_file_name,'soil_DON_deadmic',r1d, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soil_deadmicrobe_DON_ptr,idx,r1d)
           call read_compressed(restart_file_name,'soil_NO3',r1d, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soil_nitrate_ptr,idx,r1d)
           call read_compressed(restart_file_name,'soil_NH4',r1d, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soil_ammonium_ptr,idx,r1d)
           call read_compressed(restart_file_name,'soil_nitrif',r1d, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soil_nitrif_ptr,idx,r1d)
           call read_compressed(restart_file_name,'soil_denitrif',r1d, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soil_denitrif_ptr,idx,r1d)

           call read_compressed(restart_file_name,'fast_DON_leached',r0d, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soil_fast_DON_leached_ptr,idx,r0d)
           call read_compressed(restart_file_name,'slow_DON_leached',r0d, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soil_slow_DON_leached_ptr,idx,r0d)
           call read_compressed(restart_file_name,'deadmic_DON_leached',r0d, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soil_deadmic_DON_leached_ptr,idx,r0d)

           call read_compressed(restart_file_name,'leaf_litter_fast_N',r1dc, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soilc_leafLitter_fast_soil_N_ptr,idx,r1dc)
           call read_compressed(restart_file_name,'leaf_litter_slow_N',r1dc, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soilc_leafLitter_slow_soil_N_ptr,idx,r1dc)
           call read_compressed(restart_file_name,'leaf_litter_deadMic_N',r1dc, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soilc_leafLitter_deadMicrobeN_ptr,idx,r1dc)
           call read_compressed(restart_file_name,'leaf_litter_liveMic_N',r1dc, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soilc_leafLitter_livingMicrobeN_ptr,idx,r1dc)
           call read_compressed(restart_file_name,'leaf_litter_originalCohortN',r1dc, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soilc_leafLitter_originalLitterN_ptr,idx,r1dc)
           call read_compressed(restart_file_name,'leaf_litter_fastProtectedN',r1dc, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soilc_leafLitter_fast_protected_N_ptr,idx,r1dc)
           call read_compressed(restart_file_name,'leaf_litter_slowProtectedN',r1dc, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soilc_leafLitter_slow_protected_N_ptr,idx,r1dc)
           call read_compressed(restart_file_name,'leaf_litter_deadMicrobeProtectedN',r1dc, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soilc_leafLitter_deadMicrobe_protected_N_ptr,idx,r1dc)

           call read_compressed(restart_file_name,'leaf_litter_DON_fast',r0d, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soilc_leaflitter_fast_DON_ptr,idx,r0d)
           call read_compressed(restart_file_name,'leaf_litter_DON_slow',r0d, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soilc_leaflitter_slow_DON_ptr,idx,r0d)
           call read_compressed(restart_file_name,'leaf_litter_DON_deadmic',r0d, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soilc_leaflitter_deadMicrobe_DON_ptr,idx,r0d)
           call read_compressed(restart_file_name,'leaf_litter_NO3',r1d, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soilc_leaflitter_nitrate_ptr,idx,r1d)
           call read_compressed(restart_file_name,'leaf_litter_NH4',r1d, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soilc_leaflitter_ammonium_ptr,idx,r1d)
           call read_compressed(restart_file_name,'leaf_litter_nitrif',r1d, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soilc_leaflitter_nitrif_ptr,idx,r1d)
           call read_compressed(restart_file_name,'leaf_litter_denitrif',r1d, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soilc_leaflitter_denitrif_ptr,idx,r1d)

           call read_compressed(restart_file_name,'fineWood_litter_fast_N',r1dc, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soilc_fineWoodLitter_fast_soil_N_ptr,idx,r1dc)
           call read_compressed(restart_file_name,'fineWood_litter_slow_N',r1dc, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soilc_fineWoodLitter_slow_soil_N_ptr,idx,r1dc)
           call read_compressed(restart_file_name,'fineWood_litter_deadMic_N',r1dc, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soilc_fineWoodLitter_deadMicrobeN_ptr,idx,r1dc)
           call read_compressed(restart_file_name,'fineWood_litter_liveMic_N',r1dc, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soilc_fineWoodLitter_livingMicrobeN_ptr,idx,r1dc)
           call read_compressed(restart_file_name,'fineWood_litter_originalCohortN',r1dc, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soilc_fineWoodLitter_originalLitterN_ptr,idx,r1dc)
           call read_compressed(restart_file_name,'fineWood_litter_fastProtectedN',r1dc, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soilc_fineWoodLitter_fast_protected_N_ptr,idx,r1dc)
           call read_compressed(restart_file_name,'fineWood_litter_slowProtectedN',r1dc, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soilc_fineWoodLitter_slow_protected_N_ptr,idx,r1dc)
           call read_compressed(restart_file_name,'fineWood_litter_deadMicrobeProtectedN',r1dc, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soilc_fineWoodLitter_deadMicrobe_protected_N_ptr,idx,r1dc)

           call read_compressed(restart_file_name,'coarseWood_litter_fast_N',r1dc, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soilc_coarseWoodLitter_fast_soil_N_ptr,idx,r1dc)
           call read_compressed(restart_file_name,'coarseWood_litter_slow_N',r1dc, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soilc_coarseWoodLitter_slow_soil_N_ptr,idx,r1dc)
           call read_compressed(restart_file_name,'coarseWood_litter_deadMic_N',r1dc, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soilc_coarseWoodLitter_deadMicrobeN_ptr,idx,r1dc)
           call read_compressed(restart_file_name,'coarseWood_litter_liveMic_N',r1dc, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soilc_coarseWoodLitter_livingMicrobeN_ptr,idx,r1dc)
           call read_compressed(restart_file_name,'coarseWood_litter_originalCohortN',r1dc, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soilc_coarseWoodLitter_originalLitterN_ptr,idx,r1dc)
           call read_compressed(restart_file_name,'coarseWood_litter_fastProtectedN',r1dc, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soilc_coarseWoodLitter_fast_protected_N_ptr,idx,r1dc)
           call read_compressed(restart_file_name,'coarseWood_litter_slowProtectedN',r1dc, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soilc_coarseWoodLitter_slow_protected_N_ptr,idx,r1dc)
           call read_compressed(restart_file_name,'coarseWood_litter_deadMicrobeProtectedN',r1dc, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soilc_coarseWoodLitter_deadMicrobe_protected_N_ptr,idx,r1dc)
           call read_compressed(restart_file_name,'coarseWood_litter_DON_fast',r0d, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soilc_coarseWoodlitter_fast_DON_ptr,idx,r0d)
           call read_compressed(restart_file_name,'coarseWood_litter_DON_slow',r0d, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soilc_coarseWoodlitter_slow_DON_ptr,idx,r0d)
           call read_compressed(restart_file_name,'coarseWood_litter_DON_deadmic',r0d, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soilc_coarseWoodlitter_deadMicrobe_DON_ptr,idx,r0d)
           call read_compressed(restart_file_name,'coarseWood_litter_NO3',r1d, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soilc_coarseWoodlitter_nitrate_ptr,idx,r1d)
           call read_compressed(restart_file_name,'coarseWood_litter_NH4',r1d, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soilc_coarseWoodlitter_ammonium_ptr,idx,r1d)
           call read_compressed(restart_file_name,'coarseWood_litter_nitrif',r1d, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soilc_coarseWoodlitter_nitrif_ptr,idx,r1d)
           call read_compressed(restart_file_name,'coarseWood_litter_denitrif',r1d, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soilc_coarseWoodlitter_denitrif_ptr,idx,r1d)

           call read_compressed(restart_file_name,'fineWood_litter_DON_fast',r0d, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soilc_fineWoodlitter_fast_DON_ptr,idx,r0d)
           call read_compressed(restart_file_name,'fineWood_litter_DON_slow',r0d, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soilc_fineWoodlitter_slow_DON_ptr,idx,r0d)
           call read_compressed(restart_file_name,'fineWood_litter_DON_deadmic',r0d, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soilc_fineWoodlitter_deadMicrobe_DON_ptr,idx,r0d)
           call read_compressed(restart_file_name,'fineWood_litter_NO3',r1d, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soilc_fineWoodlitter_nitrate_ptr,idx,r1d)
           call read_compressed(restart_file_name,'fineWood_litter_NH4',r1d, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soilc_fineWoodlitter_ammonium_ptr,idx,r1d)
           call read_compressed(restart_file_name,'fineWood_litter_nitrif',r1d, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soilc_fineWoodlitter_nitrif_ptr,idx,r1d)
           call read_compressed(restart_file_name,'fineWood_litter_denitrif',r1d, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soilc_fineWoodlitter_denitrif_ptr,idx,r1d)
       endif

           call read_compressed(restart_file_name,'soil_DOC_fast',r1d, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soil_fast_DOC_ptr,idx,r1d)
           call read_compressed(restart_file_name,'soil_DOC_slow',r1d, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soil_slow_DOC_ptr,idx,r1d)
           call read_compressed(restart_file_name,'soil_DOC_deadmic',r1d, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soil_deadmicrobe_DOC_ptr,idx,r1d)



           if(field_exist(restart_file_name,'fast_DOC_leached',domain=lnd%domain)) then
              call read_compressed(restart_file_name,'fast_DOC_leached',r0d, domain=lnd%domain, timelevel=1)
              call assemble_tiles(soil_fast_DOC_leached_ptr,idx,r0d)
              call read_compressed(restart_file_name,'slow_DOC_leached',r0d, domain=lnd%domain, timelevel=1)
              call assemble_tiles(soil_slow_DOC_leached_ptr,idx,r0d)
              call read_compressed(restart_file_name,'deadmic_DOC_leached',r0d, domain=lnd%domain, timelevel=1)
              call assemble_tiles(soil_deadmic_DOC_leached_ptr,idx,r0d)
           endif


           call read_compressed(restart_file_name,'leaf_litter_fast_C',r1dc, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soilc_leafLitter_fast_soil_C_ptr,idx,r1dc)
           call read_compressed(restart_file_name,'leaf_litter_slow_C',r1dc, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soilc_leafLitter_slow_soil_C_ptr,idx,r1dc)
           call read_compressed(restart_file_name,'leaf_litter_deadMic_C',r1dc, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soilc_leafLitter_deadMicrobeC_ptr,idx,r1dc)
           call read_compressed(restart_file_name,'leaf_litter_liveMic_C',r1dc, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soilc_leafLitter_livingMicrobeC_ptr,idx,r1dc)
           call read_compressed(restart_file_name,'leaf_litter_CO2',r1dc, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soilc_leafLitter_CO2_ptr,idx,r1dc)
           call read_compressed(restart_file_name,'leaf_litter_Rtot',r1dc, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soilc_leafLitter_Rtot_ptr,idx,r1dc)
           call read_compressed(restart_file_name,'leaf_litter_originalCohortC',r1dc, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soilc_leafLitter_originalLitterC_ptr,idx,r1dc)
           call read_compressed(restart_file_name,'leaf_litter_fastProtectedC',r1dc, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soilc_leafLitter_fast_protected_C_ptr,idx,r1dc)
           call read_compressed(restart_file_name,'leaf_litter_slowProtectedC',r1dc, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soilc_leafLitter_slow_protected_C_ptr,idx,r1dc)
           call read_compressed(restart_file_name,'leaf_litter_deadMicrobeProtectedC',r1dc, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soilc_leafLitter_deadMicrobe_protected_C_ptr,idx,r1dc)



           call read_compressed(restart_file_name,'leaf_litter_DOC_fast',r0d, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soilc_leaflitter_fast_DOC_ptr,idx,r0d)
           call read_compressed(restart_file_name,'leaf_litter_DOC_slow',r0d, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soilc_leaflitter_slow_DOC_ptr,idx,r0d)
           call read_compressed(restart_file_name,'leaf_litter_DOC_deadmic',r0d, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soilc_leaflitter_deadMicrobe_DOC_ptr,idx,r0d)




           call read_compressed(restart_file_name,'fineWood_litter_fast_C',r1dc, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soilc_fineWoodLitter_fast_soil_C_ptr,idx,r1dc)
           call read_compressed(restart_file_name,'fineWood_litter_slow_C',r1dc, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soilc_fineWoodLitter_slow_soil_C_ptr,idx,r1dc)
           call read_compressed(restart_file_name,'fineWood_litter_deadMic_C',r1dc, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soilc_fineWoodLitter_deadMicrobeC_ptr,idx,r1dc)
           call read_compressed(restart_file_name,'fineWood_litter_liveMic_C',r1dc, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soilc_fineWoodLitter_livingMicrobeC_ptr,idx,r1dc)
           call read_compressed(restart_file_name,'fineWood_litter_CO2',r1dc, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soilc_fineWoodLitter_CO2_ptr,idx,r1dc)
           call read_compressed(restart_file_name,'fineWood_litter_Rtot',r1dc, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soilc_fineWoodLitter_Rtot_ptr,idx,r1dc)
           call read_compressed(restart_file_name,'fineWood_litter_originalCohortC',r1dc, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soilc_fineWoodLitter_originalLitterC_ptr,idx,r1dc)
           call read_compressed(restart_file_name,'fineWood_litter_fastProtectedC',r1dc, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soilc_fineWoodLitter_fast_protected_C_ptr,idx,r1dc)
           call read_compressed(restart_file_name,'fineWood_litter_slowProtectedC',r1dc, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soilc_fineWoodLitter_slow_protected_C_ptr,idx,r1dc)
           call read_compressed(restart_file_name,'fineWood_litter_deadMicrobeProtectedC',r1dc, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soilc_fineWoodLitter_deadMicrobe_protected_C_ptr,idx,r1dc)



           call read_compressed(restart_file_name,'fineWood_litter_DOC_fast',r0d, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soilc_fineWoodlitter_fast_DOC_ptr,idx,r0d)
           call read_compressed(restart_file_name,'fineWood_litter_DOC_slow',r0d, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soilc_fineWoodlitter_slow_DOC_ptr,idx,r0d)
           call read_compressed(restart_file_name,'fineWood_litter_DOC_deadmic',r0d, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soilc_fineWoodlitter_deadMicrobe_DOC_ptr,idx,r0d)



           call read_compressed(restart_file_name,'coarseWood_litter_fast_C',r1dc, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soilc_coarseWoodLitter_fast_soil_C_ptr,idx,r1dc)
           call read_compressed(restart_file_name,'coarseWood_litter_slow_C',r1dc, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soilc_coarseWoodLitter_slow_soil_C_ptr,idx,r1dc)
           call read_compressed(restart_file_name,'coarseWood_litter_deadMic_C',r1dc, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soilc_coarseWoodLitter_deadMicrobeC_ptr,idx,r1dc)
           call read_compressed(restart_file_name,'coarseWood_litter_liveMic_C',r1dc, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soilc_coarseWoodLitter_livingMicrobeC_ptr,idx,r1dc)
           call read_compressed(restart_file_name,'coarseWood_litter_CO2',r1dc, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soilc_coarseWoodLitter_CO2_ptr,idx,r1dc)
           call read_compressed(restart_file_name,'coarseWood_litter_Rtot',r1dc, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soilc_coarseWoodLitter_Rtot_ptr,idx,r1dc)
           call read_compressed(restart_file_name,'coarseWood_litter_originalCohortC',r1dc, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soilc_coarseWoodLitter_originalLitterC_ptr,idx,r1dc)
           call read_compressed(restart_file_name,'coarseWood_litter_fastProtectedC',r1dc, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soilc_coarseWoodLitter_fast_protected_C_ptr,idx,r1dc)
           call read_compressed(restart_file_name,'coarseWood_litter_slowProtectedC',r1dc, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soilc_coarseWoodLitter_slow_protected_C_ptr,idx,r1dc)
           call read_compressed(restart_file_name,'coarseWood_litter_deadMicrobeProtectedC',r1dc, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soilc_coarseWoodLitter_deadMicrobe_protected_C_ptr,idx,r1dc)



           call read_compressed(restart_file_name,'coarseWood_litter_DOC_fast',r0d, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soilc_coarseWoodlitter_fast_DOC_ptr,idx,r0d)
           call read_compressed(restart_file_name,'coarseWood_litter_DOC_slow',r0d, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soilc_coarseWoodlitter_slow_DOC_ptr,idx,r0d)
           call read_compressed(restart_file_name,'coarseWood_litter_DOC_deadmic',r0d, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soilc_coarseWoodlitter_deadMicrobe_DOC_ptr,idx,r0d)



           if(field_exist(restart_file_name,'is_peat',domain=lnd%domain)) then
              call read_compressed(restart_file_name,'is_peat',i1d, domain=lnd%domain, timelevel=1)
              call assemble_tiles(soil_is_peat_ptr,idx,i1d)
           endif
           deallocate(r2d, r1dc)
        endif
        deallocate(idx,r1d,r0d,i1d)
     else
        call error_mesg('soil_init', 'reading NetCDF restart "'//trim(restart_file_name)//'"', NOTE)
        __NF_ASRT__(nf_open(restart_file_name,NF_NOWRITE,unit))
        call read_tile_data_r1d_fptr(unit, 'temp'         , soil_T_ptr  )
        call read_tile_data_r1d_fptr(unit, 'wl'           , soil_wl_ptr )
        call read_tile_data_r1d_fptr(unit, 'ws'           , soil_ws_ptr )
        call read_tile_data_r1d_fptr(unit, 'groundwater'  , soil_groundwater_ptr )
        call read_tile_data_r1d_fptr(unit, 'groundwater_T', soil_groundwater_T_ptr)
        if(nfu_inq_var(unit, 'uptake_T')==NF_NOERR) call read_tile_data_r0d_fptr(unit, 'uptake_T', soil_uptake_T_ptr)

        if (soil_carbon_option==SOILC_CENTURY.or.soil_carbon_option==SOILC_CENTURY_BY_LAYER) then
           if(nfu_inq_var(unit, 'fsc')==NF_NOERR) then
              call read_tile_data_r1d_fptr(unit,'fsc',soil_fast_soil_C_ptr)
              call read_tile_data_r1d_fptr(unit,'ssc',soil_slow_soil_C_ptr)
           else
              ! try to read fsc and ssc from vegetation restart
              call get_input_restart_name('INPUT/vegn2.res.nc',restart_exists,restart_file_name)
              if (restart_exists) then
                 __NF_ASRT__(nf_open(restart_file_name,NF_NOWRITE,unit1))
                 ! read old (scalar) fsc and ssc into the first element of the fast_soil_C
                 ! and slow_soil_C arrays
                 call read_tile_data_r1d_fptr(unit1,'fsc',soil_fast_soil_C_ptr,1)
                 call read_tile_data_r1d_fptr(unit1,'ssc',soil_slow_soil_C_ptr,1)
              endif
           endif
        endif

        if (nfu_inq_var(unit,'fast_soil_C')==NF_NOERR) then
           call read_tile_data_layered_cohort_fptr(unit, 'fast_soil_C', soilc_fast_soil_C_ptr)
           call read_tile_data_layered_cohort_fptr(unit, 'slow_soil_C', soilc_slow_soil_C_ptr)
           call read_tile_data_layered_cohort_fptr(unit, 'deadMic'    , soilc_deadMicrobeC_ptr)
           call read_tile_data_layered_cohort_fptr(unit, 'fastProtectedC'    , soilc_fast_protected_C_ptr)
           call read_tile_data_layered_cohort_fptr(unit, 'slowProtectedC'    , soilc_slow_protected_C_ptr)
           call read_tile_data_layered_cohort_fptr(unit, 'deadMicrobeProtectedC'    , soilc_deadMicrobe_protected_C_ptr)
           call read_tile_data_layered_cohort_fptr(unit, 'liveMic'    , soilc_livingMicrobeC_ptr)
           call read_tile_data_layered_cohort_fptr(unit, 'CO2'        , soilc_CO2_ptr)
           call read_tile_data_layered_cohort_fptr(unit, 'Rtot'       , soilc_Rtot_ptr)
           call read_tile_data_layered_cohort_fptr(unit, 'originalCohortC',soilc_originalLitterC_ptr)

if (soil_carbon_option == SOILC_CORPSE_N) then

           call read_tile_data_layered_cohort_fptr(unit, 'fast_soil_N', soilc_fast_soil_N_ptr)
           call read_tile_data_layered_cohort_fptr(unit, 'slow_soil_N', soilc_slow_soil_N_ptr)
           call read_tile_data_layered_cohort_fptr(unit, 'deadMic_N'    , soilc_deadMicrobeN_ptr)
           call read_tile_data_layered_cohort_fptr(unit, 'fastProtectedN'    , soilc_fast_protected_N_ptr)
           call read_tile_data_layered_cohort_fptr(unit, 'slowProtectedN'    , soilc_slow_protected_N_ptr)
           call read_tile_data_layered_cohort_fptr(unit, 'deadMicrobeProtectedN'    , soilc_deadMicrobe_protected_N_ptr)
           call read_tile_data_layered_cohort_fptr(unit, 'liveMicN'    , soilc_livingMicrobeN_ptr)
           call read_tile_data_layered_cohort_fptr(unit, 'originalCohortN',soilc_originalLitterN_ptr)
           ! Leaving out cohort-level immobilization and mineralization fields for now -- BNS

           call read_tile_data_r1d_fptr(unit, 'soil_DON_fast',soil_fast_DON_ptr)
           call read_tile_data_r1d_fptr(unit, 'soil_DON_slow',soil_slow_DON_ptr)
           call read_tile_data_r1d_fptr(unit, 'soil_DON_deadmic',soil_deadmicrobe_DON_ptr)
           call read_tile_data_r1d_fptr(unit, 'soil_NO3',soil_nitrate_ptr)
           call read_tile_data_r1d_fptr(unit, 'soil_NH4',soil_ammonium_ptr)
           call read_tile_data_r1d_fptr(unit, 'soil_nitrif',soil_nitrif_ptr)
           call read_tile_data_r1d_fptr(unit, 'soil_denitrif',soil_denitrif_ptr)

           call read_tile_data_r1d_fptr(unit, 'leaf_litter_fast_N'    , soilc_leafLitter_fast_soil_N_ptr)
           call read_tile_data_r1d_fptr(unit, 'leaf_litter_slow_N'    , soilc_leafLitter_slow_soil_N_ptr)
           call read_tile_data_r1d_fptr(unit, 'leaf_litter_deadMic_N'    , soilc_leafLitter_deadMicrobeN_ptr)
           call read_tile_data_r1d_fptr(unit, 'leaf_litter_liveMic_N'    , soilc_leafLitter_livingMicrobeN_ptr)
           call read_tile_data_r1d_fptr(unit, 'leaf_litter_originalCohortN',soilc_leafLitter_originalLitterN_ptr)
           call read_tile_data_r1d_fptr(unit, 'leaf_litter_fastProtectedN',soilc_leafLitter_fast_protected_N_ptr)
           call read_tile_data_r1d_fptr(unit, 'leaf_litter_slowProtectedN',soilc_leafLitter_slow_protected_N_ptr)
           call read_tile_data_r1d_fptr(unit, 'leaf_litter_deadMicrobeProtectedN',soilc_leafLitter_deadMicrobe_protected_N_ptr)

           call read_tile_data_r1d_fptr(unit, 'fineWood_litter_fast_N'    , soilc_fineWoodLitter_fast_soil_N_ptr)
           call read_tile_data_r1d_fptr(unit, 'fineWood_litter_slow_N'    , soilc_fineWoodLitter_slow_soil_N_ptr)
           call read_tile_data_r1d_fptr(unit, 'fineWood_litter_deadMic_N'    , soilc_fineWoodLitter_deadMicrobeN_ptr)
           call read_tile_data_r1d_fptr(unit, 'fineWood_litter_liveMic_N'    , soilc_fineWoodLitter_livingMicrobeN_ptr)
           call read_tile_data_r1d_fptr(unit, 'fineWood_litter_originalCohortN',soilc_fineWoodLitter_originalLitterN_ptr)
           call read_tile_data_r1d_fptr(unit, 'fineWood_litter_fastProtectedN',soilc_fineWoodLitter_fast_protected_N_ptr)
           call read_tile_data_r1d_fptr(unit, 'fineWood_litter_slowProtectedN',soilc_fineWoodLitter_slow_protected_N_ptr)
           call read_tile_data_r1d_fptr(unit, 'fineWood_litter_deadMicrobeProtectedN',soilc_fineWoodLitter_deadMicrobe_protected_N_ptr)

           call read_tile_data_r0d_fptr(unit, 'fineWood_litter_DON_fast',soilc_fineWoodLitter_fast_DON_ptr)
           call read_tile_data_r0d_fptr(unit, 'fineWood_litter_DON_slow',soilc_fineWoodLitter_slow_DON_ptr)
           call read_tile_data_r0d_fptr(unit, 'fineWood_litter_DON_deadmic',soilc_fineWoodLitter_deadmicrobe_DON_ptr)
           call read_tile_data_r0d_fptr(unit, 'fineWood_litter_NO3',soilc_fineWoodLitter_nitrate_ptr)
           call read_tile_data_r0d_fptr(unit, 'fineWood_litter_NH4',soilc_fineWoodLitter_ammonium_ptr)
           call read_tile_data_r0d_fptr(unit, 'fineWood_litter_nitrif',soilc_fineWoodLitter_nitrif_ptr)
           call read_tile_data_r0d_fptr(unit, 'fineWood_litter_denitrif',soilc_fineWoodLitter_denitrif_ptr)

           call read_tile_data_r0d_fptr(unit, 'leaf_litter_DON_fast',soilc_leafLitter_fast_DON_ptr)
           call read_tile_data_r0d_fptr(unit, 'leaf_litter_DON_slow',soilc_leafLitter_slow_DON_ptr)
           call read_tile_data_r0d_fptr(unit, 'leaf_litter_DON_deadmic',soilc_leafLitter_deadmicrobe_DON_ptr)
           call read_tile_data_r0d_fptr(unit, 'leaf_litter_NO3',soilc_leafLitter_nitrate_ptr)
           call read_tile_data_r0d_fptr(unit, 'leaf_litter_NH4',soilc_leafLitter_ammonium_ptr)
           call read_tile_data_r0d_fptr(unit, 'leaf_litter_nitrif',soilc_leafLitter_nitrif_ptr)
           call read_tile_data_r0d_fptr(unit, 'leaf_litter_denitrif',soilc_leafLitter_denitrif_ptr)

           call read_tile_data_r1d_fptr(unit, 'coarseWood_litter_fast_N'    , soilc_coarseWoodLitter_fast_soil_N_ptr)
           call read_tile_data_r1d_fptr(unit, 'coarseWood_litter_slow_N'    , soilc_coarseWoodLitter_slow_soil_N_ptr)
           call read_tile_data_r1d_fptr(unit, 'coarseWood_litter_deadMic_N'    , soilc_coarseWoodLitter_deadMicrobeN_ptr)
           call read_tile_data_r1d_fptr(unit, 'coarseWood_litter_liveMic_N'    , soilc_coarseWoodLitter_livingMicrobeN_ptr)
           call read_tile_data_r1d_fptr(unit, 'coarseWood_litter_originalCohortN',soilc_coarseWoodLitter_originalLitterN_ptr)
           call read_tile_data_r1d_fptr(unit, 'coarseWood_litter_fastProtectedN',soilc_coarseWoodLitter_fast_protected_N_ptr)
           call read_tile_data_r1d_fptr(unit, 'coarseWood_litter_slowProtectedN',soilc_coarseWoodLitter_slow_protected_N_ptr)
           call read_tile_data_r1d_fptr(unit, 'coarseWood_litter_deadMicrobeProtectedN',soilc_coarseWoodLitter_deadMicrobe_protected_N_ptr)

           call read_tile_data_r0d_fptr(unit, 'coarseWood_litter_DON_fast',soilc_coarseWoodLitter_fast_DON_ptr)
           call read_tile_data_r0d_fptr(unit, 'coarseWood_litter_DON_slow',soilc_coarseWoodLitter_slow_DON_ptr)
           call read_tile_data_r0d_fptr(unit, 'coarseWood_litter_DON_deadmic',soilc_coarseWoodLitter_deadmicrobe_DON_ptr)
           call read_tile_data_r0d_fptr(unit, 'coarseWood_litter_NO3',soilc_coarseWoodLitter_nitrate_ptr)
           call read_tile_data_r0d_fptr(unit, 'coarseWood_litter_NH4',soilc_coarseWoodLitter_ammonium_ptr)
           call read_tile_data_r0d_fptr(unit, 'coarseWood_litter_nitrif',soilc_coarseWoodLitter_nitrif_ptr)
           call read_tile_data_r0d_fptr(unit, 'coarseWood_litter_denitrif',soilc_coarseWoodLitter_denitrif_ptr)


       endif

           call read_tile_data_r1d_fptr(unit, 'soil_DOC_fast',soil_fast_DOC_ptr)
           call read_tile_data_r1d_fptr(unit, 'soil_DOC_slow',soil_slow_DOC_ptr)
           call read_tile_data_r1d_fptr(unit, 'soil_DOC_deadmic',soil_deadmicrobe_DOC_ptr)



           if(nfu_inq_var(unit, 'fast_DOC_leached')==NF_NOERR) then
              call read_tile_data_r0d_fptr(unit,'fast_DOC_leached',     soil_fast_DOC_leached_ptr)
              call read_tile_data_r0d_fptr(unit,'slow_DOC_leached',     soil_slow_DOC_leached_ptr)
              call read_tile_data_r0d_fptr(unit,'deadmic_DOC_leached',     soil_deadmic_DOC_leached_ptr)
           endif
           if(nfu_inq_var(unit, 'fast_DON_leached')==NF_NOERR) then
              call read_tile_data_r0d_fptr(unit,'fast_DON_leached',     soil_fast_DON_leached_ptr)
              call read_tile_data_r0d_fptr(unit,'slow_DON_leached',     soil_slow_DON_leached_ptr)
              call read_tile_data_r0d_fptr(unit,'deadmic_DON_leached',     soil_deadmic_DON_leached_ptr)
              call read_tile_data_r0d_fptr(unit,'NO3_leached',     soil_NO3_leached_ptr)
              call read_tile_data_r0d_fptr(unit,'NH4_leached',     soil_NH4_leached_ptr)
           endif

           call read_tile_data_r1d_fptr(unit, 'leaf_litter_fast_C'    , soilc_leafLitter_fast_soil_C_ptr)
           call read_tile_data_r1d_fptr(unit, 'leaf_litter_slow_C'    , soilc_leafLitter_slow_soil_C_ptr)
           call read_tile_data_r1d_fptr(unit, 'leaf_litter_deadMic_C'    , soilc_leafLitter_deadMicrobeC_ptr)
           call read_tile_data_r1d_fptr(unit, 'leaf_litter_liveMic_C'    , soilc_leafLitter_livingMicrobeC_ptr)
           call read_tile_data_r1d_fptr(unit, 'leaf_litter_CO2'        , soilc_leafLitter_CO2_ptr)
           call read_tile_data_r1d_fptr(unit, 'leaf_litter_Rtot'       , soilc_leafLitter_Rtot_ptr)
           call read_tile_data_r1d_fptr(unit, 'leaf_litter_originalCohortC',soilc_leafLitter_originalLitterC_ptr)
           call read_tile_data_r1d_fptr(unit, 'leaf_litter_fastProtectedC',soilc_leafLitter_fast_protected_C_ptr)
           call read_tile_data_r1d_fptr(unit, 'leaf_litter_slowProtectedC',soilc_leafLitter_slow_protected_C_ptr)
           call read_tile_data_r1d_fptr(unit, 'leaf_litter_deadMicrobeProtectedC',soilc_leafLitter_deadMicrobe_protected_C_ptr)

           call read_tile_data_r0d_fptr(unit, 'leaf_litter_DOC_fast',soilc_leaflitter_fast_DOC_ptr)
           call read_tile_data_r0d_fptr(unit, 'leaf_litter_DOC_slow',soilc_leaflitter_slow_DOC_ptr)
           call read_tile_data_r0d_fptr(unit, 'leaf_litter_DOC_deadmic',soilc_leaflitter_deadMicrobe_DOC_ptr)




           call read_tile_data_r1d_fptr(unit, 'fineWood_litter_fast_C'    , soilc_fineWoodLitter_fast_soil_C_ptr)
           call read_tile_data_r1d_fptr(unit, 'fineWood_litter_slow_C'    , soilc_fineWoodLitter_slow_soil_C_ptr)
           call read_tile_data_r1d_fptr(unit, 'fineWood_litter_deadMic_C'    , soilc_fineWoodLitter_deadMicrobeC_ptr)
           call read_tile_data_r1d_fptr(unit, 'fineWood_litter_liveMic_C'    , soilc_fineWoodLitter_livingMicrobeC_ptr)
           call read_tile_data_r1d_fptr(unit, 'fineWood_litter_CO2'        , soilc_fineWoodLitter_CO2_ptr)
           call read_tile_data_r1d_fptr(unit, 'fineWood_litter_Rtot'       , soilc_fineWoodLitter_Rtot_ptr)
           call read_tile_data_r1d_fptr(unit, 'fineWood_litter_originalCohortC',soilc_fineWoodLitter_originalLitterC_ptr)
           call read_tile_data_r1d_fptr(unit, 'fineWood_litter_fastProtectedC',soilc_fineWoodLitter_fast_protected_C_ptr)
           call read_tile_data_r1d_fptr(unit, 'fineWood_litter_slowProtectedC',soilc_fineWoodLitter_slow_protected_C_ptr)
           call read_tile_data_r1d_fptr(unit, 'fineWood_litter_deadMicrobeProtectedC',soilc_fineWoodLitter_deadMicrobe_protected_C_ptr)

           call read_tile_data_r0d_fptr(unit, 'fineWood_litter_DOC_fast',soilc_fineWoodlitter_fast_DOC_ptr)
           call read_tile_data_r0d_fptr(unit, 'fineWood_litter_DOC_slow',soilc_fineWoodlitter_slow_DOC_ptr)
           call read_tile_data_r0d_fptr(unit, 'fineWood_litter_DOC_deadmic',soilc_fineWoodlitter_deadMicrobe_DOC_ptr)


           call read_tile_data_r1d_fptr(unit, 'coarseWood_litter_fast_C'    , soilc_coarseWoodLitter_fast_soil_C_ptr)
           call read_tile_data_r1d_fptr(unit, 'coarseWood_litter_slow_C'    , soilc_coarseWoodLitter_slow_soil_C_ptr)
           call read_tile_data_r1d_fptr(unit, 'coarseWood_litter_deadMic_C'    , soilc_coarseWoodLitter_deadMicrobeC_ptr)
           call read_tile_data_r1d_fptr(unit, 'coarseWood_litter_liveMic_C'    , soilc_coarseWoodLitter_livingMicrobeC_ptr)
           call read_tile_data_r1d_fptr(unit, 'coarseWood_litter_CO2'        , soilc_coarseWoodLitter_CO2_ptr)
           call read_tile_data_r1d_fptr(unit, 'coarseWood_litter_Rtot'       , soilc_coarseWoodLitter_Rtot_ptr)
           call read_tile_data_r1d_fptr(unit, 'coarseWood_litter_originalCohortC',soilc_coarseWoodLitter_originalLitterC_ptr)
           call read_tile_data_r1d_fptr(unit, 'coarseWood_litter_fastProtectedC',soilc_coarseWoodLitter_fast_protected_C_ptr)
           call read_tile_data_r1d_fptr(unit, 'coarseWood_litter_slowProtectedC',soilc_coarseWoodLitter_slow_protected_C_ptr)
           call read_tile_data_r1d_fptr(unit, 'coarseWood_litter_deadMicrobeProtectedC',soilc_coarseWoodLitter_deadMicrobe_protected_C_ptr)

           call read_tile_data_r0d_fptr(unit, 'coarseWood_litter_DOC_fast',soilc_coarseWoodlitter_fast_DOC_ptr)
           call read_tile_data_r0d_fptr(unit, 'coarseWood_litter_DOC_slow',soilc_coarseWoodlitter_slow_DOC_ptr)
           call read_tile_data_r0d_fptr(unit, 'coarseWood_litter_DOC_deadmic',soilc_coarseWoodlitter_deadMicrobe_DOC_ptr)



           if(nfu_inq_var(unit, 'is_peat')==NF_NOERR) then
              call read_tile_data_i1d_fptr_all(unit, 'is_peat',soil_is_peat_ptr)
           endif
        endif
        __NF_ASRT__(nf_close(unit))
     endif
  else
     call error_mesg('soil_init', 'cold-starting soil', NOTE)
  endif

  ! read soil carbon restart, if present
  call get_input_restart_name('INPUT/soil_carbon.res.nc',restart_exists,restart_file_name)
  if (restart_exists) then
     call error_mesg('veg_data_init','reading soil_carbon restart',NOTE)
     if(new_land_io) then
        restart_file_name = 'INPUT/soil_carbon.res.nc'
        call get_field_size(restart_file_name, 'tile_index', siz, field_found=found, domain=lnd%domain)
        if ( .not.found ) call error_mesg(trim(module_name), &
             'tile_index axis not found in '//trim(restart_file_name), FATAL)
        isize = siz(1)

        call get_field_size(restart_file_name, 'zfull', siz, field_found=found, domain=lnd%domain)
        if ( .not.found ) call error_mesg(trim(module_name), &
             'Z axis not found in '//trim(restart_file_name), FATAL)
        nz = siz(1)

        allocate(idx(isize),r0d(isize),r1d(isize,nz))
        call read_compressed(restart_file_name,'tile_index',idx, domain=lnd%domain, timelevel=1)

        call read_compressed(restart_file_name,'asoil_in',r1d, domain=lnd%domain, timelevel=1)
        call assemble_tiles(soil_asoil_in_ptr,idx,r1d)

        call read_compressed(restart_file_name,'fsc_in',r1d, domain=lnd%domain, timelevel=1)
        call assemble_tiles(soil_fsc_in_ptr,idx,r1d)

        call read_compressed(restart_file_name,'ssc_in',r1d, domain=lnd%domain, timelevel=1)
        call assemble_tiles(soil_ssc_in_ptr,idx,r1d)

        if ((soil_carbon_option==SOILC_CORPSE .or. soil_carbon_option==SOILC_CORPSE_N) .and. field_exist(restart_file_name,'deadmic_C_in', domain=lnd%domain)) then
              ! Note: Probably want to make backward compatible with restarts that have deadmic_in
           call read_compressed(restart_file_name,'deadmic_C_in',r1d, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soil_deadmic_C_in_ptr,idx,r1d)
           call read_compressed(restart_file_name,'fast_protected_C_in',r1d, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soil_fast_protected_C_in_ptr,idx,r1d)
           call read_compressed(restart_file_name,'slow_protected_C_in',r1d, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soil_slow_protected_C_in_ptr,idx,r1d)
           call read_compressed(restart_file_name,'deadmic_protected_C_in',r1d, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soil_deadmic_protected_C_in_ptr,idx,r1d)


           call read_compressed(restart_file_name,'leaflitter_fsc_in',r0d, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soil_leaflitter_fsc_in_ptr,idx,r0d)
           call read_compressed(restart_file_name,'leaflitter_ssc_in',r0d, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soil_leaflitter_ssc_in_ptr,idx,r0d)
           call read_compressed(restart_file_name,'leaflitter_deadmic_C_in',r0d, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soil_leaflitter_deadmic_C_in_ptr,idx,r0d)


           call read_compressed(restart_file_name,'finewoodlitter_fsc_in',r0d, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soil_finewoodlitter_fsc_in_ptr,idx,r0d)
           call read_compressed(restart_file_name,'finewoodlitter_ssc_in',r0d, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soil_finewoodlitter_ssc_in_ptr,idx,r0d)
           call read_compressed(restart_file_name,'finewoodlitter_deadmic_C_in',r0d, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soil_finewoodlitter_deadmic_C_in_ptr,idx,r0d)


           call read_compressed(restart_file_name,'coarsewoodlitter_fsc_in',r0d, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soil_coarsewoodlitter_fsc_in_ptr,idx,r0d)
           call read_compressed(restart_file_name,'coarsewoodlitter_ssc_in',r0d, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soil_coarsewoodlitter_ssc_in_ptr,idx,r0d)
           call read_compressed(restart_file_name,'coarsewoodlitter_deadmic_in',r0d, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soil_coarsewoodlitter_deadmic_C_in_ptr,idx,r0d)


           call read_compressed(restart_file_name,'fast_C_turnover_accumulated',r1d, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soil_fast_C_turnover_accumulated_ptr,idx,r1d)
           call read_compressed(restart_file_name,'slow_C_turnover_accumulated',r1d, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soil_slow_C_turnover_accumulated_ptr,idx,r1d)
           call read_compressed(restart_file_name,'deadmic_C_turnover_accumulated',r1d, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soil_deadmic_C_turnover_accumulated_ptr,idx,r1d)


           call read_compressed(restart_file_name,'fast_protected_C_turnover_accumulated',r1d, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soil_fast_protected_C_turnover_accumulated_ptr,idx,r1d)
           call read_compressed(restart_file_name,'slow_C_protected_C_turnover_accumulated',r1d, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soil_slow_protected_C_turnover_accumulated_ptr,idx,r1d)
           call read_compressed(restart_file_name,'deadmic_protected_C_turnover_accumulated',r1d, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soil_deadmic_protected_C_turnover_accumulated_ptr,idx,r1d)


           call read_compressed(restart_file_name,'leaflitter_fast_C_turnover_accumulated',r0d, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soil_leaflitter_fast_C_turnover_accumulated_ptr,idx,r0d)
           call read_compressed(restart_file_name,'leaflitter_slow_C_turnover_accumulated',r0d, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soil_leaflitter_slow_C_turnover_accumulated_ptr,idx,r0d)
           call read_compressed(restart_file_name,'leaflitter_deadmic_C_turnover_accumulated',r0d, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soil_leaflitter_deadmic_C_turnover_accumulated_ptr,idx,r0d)
           call read_compressed(restart_file_name,'finewoodlitter_fast_C_turnover_accumulated',r0d, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soil_finewoodlitter_fast_C_turnover_accumulated_ptr,idx,r0d)
           call read_compressed(restart_file_name,'finewoodlitter_slow_C_turnover_accumulated',r0d, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soil_finewoodlitter_slow_C_turnover_accumulated_ptr,idx,r0d)
           call read_compressed(restart_file_name,'finewoodlitter_deadmic_C_turnover_accumulated',r0d, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soil_finewoodlitter_deadmic_C_turnover_accumulated_ptr,idx,r0d)
           call read_compressed(restart_file_name,'coarsewoodlitter_fast_C_turnover_accumulated',r0d, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soil_coarsewoodlitter_fast_C_turnover_accumulated_ptr,idx,r0d)
           call read_compressed(restart_file_name,'coarsewoodlitter_slow_C_turnover_accumulated',r0d, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soil_coarsewoodlitter_slow_C_turnover_accumulated_ptr,idx,r0d)
           call read_compressed(restart_file_name,'coarsewoodlitter_deadmic_C_turnover_accumulated',r0d, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soil_coarsewoodlitter_deadmic_C_turnover_accumulated_ptr,idx,r0d)

           IF(soil_carbon_option == SOILC_CORPSE_N) THEN
               call read_compressed(restart_file_name,'fsn_in',r1d, domain=lnd%domain, timelevel=1)
               call assemble_tiles(soil_fsn_in_ptr,idx,r1d)
               call read_compressed(restart_file_name,'ssn_in',r1d, domain=lnd%domain, timelevel=1)
               call assemble_tiles(soil_ssn_in_ptr,idx,r1d)

               call read_compressed(restart_file_name,'deadmic_N_in',r1d, domain=lnd%domain, timelevel=1)
               call assemble_tiles(soil_deadmic_N_in_ptr,idx,r1d)
               call read_compressed(restart_file_name,'fast_protected_N_in',r1d, domain=lnd%domain, timelevel=1)
               call assemble_tiles(soil_fast_protected_N_in_ptr,idx,r1d)
               call read_compressed(restart_file_name,'slow_protected_N_in',r1d, domain=lnd%domain, timelevel=1)
               call assemble_tiles(soil_slow_protected_N_in_ptr,idx,r1d)
               call read_compressed(restart_file_name,'deadmic_protected_N_in',r1d, domain=lnd%domain, timelevel=1)
               call assemble_tiles(soil_deadmic_protected_N_in_ptr,idx,r1d)

               call read_compressed(restart_file_name,'leaflitter_fsn_in',r0d, domain=lnd%domain, timelevel=1)
               call assemble_tiles(soil_leaflitter_fsn_in_ptr,idx,r0d)
               call read_compressed(restart_file_name,'leaflitter_ssn_in',r0d, domain=lnd%domain, timelevel=1)
               call assemble_tiles(soil_leaflitter_ssn_in_ptr,idx,r0d)
               call read_compressed(restart_file_name,'leaflitter_deadmic_N_in',r0d, domain=lnd%domain, timelevel=1)
               call assemble_tiles(soil_leaflitter_deadmic_N_in_ptr,idx,r0d)

               call read_compressed(restart_file_name,'finewoodlitter_fsn_in',r0d, domain=lnd%domain, timelevel=1)
               call assemble_tiles(soil_finewoodlitter_fsn_in_ptr,idx,r0d)
               call read_compressed(restart_file_name,'finewoodlitter_ssn_in',r0d, domain=lnd%domain, timelevel=1)
               call assemble_tiles(soil_finewoodlitter_ssn_in_ptr,idx,r0d)
               call read_compressed(restart_file_name,'finewoodlitter_deadmic_N_in',r0d, domain=lnd%domain, timelevel=1)
               call assemble_tiles(soil_finewoodlitter_deadmic_N_in_ptr,idx,r0d)

               call read_compressed(restart_file_name,'coarsewoodlitter_fsn_in',r0d, domain=lnd%domain, timelevel=1)
               call assemble_tiles(soil_coarsewoodlitter_fsn_in_ptr,idx,r0d)
               call read_compressed(restart_file_name,'coarsewoodlitter_ssn_in',r0d, domain=lnd%domain, timelevel=1)
               call assemble_tiles(soil_coarsewoodlitter_ssn_in_ptr,idx,r0d)
               call read_compressed(restart_file_name,'coarsewoodlitter_deadmic_N_in',r0d, domain=lnd%domain, timelevel=1)
               call assemble_tiles(soil_coarsewoodlitter_deadmic_N_in_ptr,idx,r0d)

               call read_compressed(restart_file_name,'fast_N_turnover_accumulated',r1d, domain=lnd%domain, timelevel=1)
               call assemble_tiles(soil_fast_N_turnover_accumulated_ptr,idx,r1d)
               call read_compressed(restart_file_name,'slow_N_turnover_accumulated',r1d, domain=lnd%domain, timelevel=1)
               call assemble_tiles(soil_slow_N_turnover_accumulated_ptr,idx,r1d)
               call read_compressed(restart_file_name,'deadmic_N_turnover_accumulated',r1d, domain=lnd%domain, timelevel=1)
               call assemble_tiles(soil_deadmic_N_turnover_accumulated_ptr,idx,r1d)

               call read_compressed(restart_file_name,'fast_protected_N_turnover_accumulated',r1d, domain=lnd%domain, timelevel=1)
               call assemble_tiles(soil_fast_protected_N_turnover_accumulated_ptr,idx,r1d)
               call read_compressed(restart_file_name,'slow_protected_N_turnover_accumulated',r1d, domain=lnd%domain, timelevel=1)
               call assemble_tiles(soil_slow_protected_N_turnover_accumulated_ptr,idx,r1d)
               call read_compressed(restart_file_name,'deadmic_protected_N_turnover_accumulated',r1d, domain=lnd%domain, timelevel=1)
               call assemble_tiles(soil_deadmic_protected_N_turnover_accumulated_ptr,idx,r1d)

           call read_compressed(restart_file_name,'leaflitter_fast_N_turnover_accumulated',r0d, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soil_leaflitter_fast_N_turnover_accumulated_ptr,idx,r0d)
           call read_compressed(restart_file_name,'leaflitter_slow_N_turnover_accumulated',r0d, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soil_leaflitter_slow_N_turnover_accumulated_ptr,idx,r0d)
           call read_compressed(restart_file_name,'leaflitter_deadmic_N_turnover_accumulated',r0d, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soil_leaflitter_deadmic_N_turnover_accumulated_ptr,idx,r0d)
           call read_compressed(restart_file_name,'finewoodlitter_fast_N_turnover_accumulated',r0d, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soil_finewoodlitter_fast_N_turnover_accumulated_ptr,idx,r0d)
           call read_compressed(restart_file_name,'finewoodlitter_slow_N_turnover_accumulated',r0d, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soil_finewoodlitter_slow_N_turnover_accumulated_ptr,idx,r0d)
           call read_compressed(restart_file_name,'finewoodlitter_deadmic_N_turnover_accumulated',r0d, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soil_finewoodlitter_deadmic_N_turnover_accumulated_ptr,idx,r0d)
           call read_compressed(restart_file_name,'coarsewoodlitter_fast_N_turnover_accumulated',r0d, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soil_coarsewoodlitter_fast_N_turnover_accumulated_ptr,idx,r0d)
           call read_compressed(restart_file_name,'coarsewoodlitter_slow_N_turnover_accumulated',r0d, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soil_coarsewoodlitter_slow_N_turnover_accumulated_ptr,idx,r0d)
           call read_compressed(restart_file_name,'coarsewoodlitter_deadmic_N_turnover_accumulated',r0d, domain=lnd%domain, timelevel=1)
           call assemble_tiles(soil_coarsewoodlitter_deadmic_N_turnover_accumulated_ptr,idx,r0d)

       ENDIF  ! End of SOILC_CORPSE_N option
        endif
        deallocate(idx,r0d,r1d)
     else
        __NF_ASRT__(nf_open(restart_file_name,NF_NOWRITE,unit))
        call error_mesg('veg_data_init','reading soil_carbon restart',NOTE)
        call read_tile_data_r1d_fptr(unit,'asoil_in',soil_asoil_in_ptr)
        call read_tile_data_r1d_fptr(unit,'fsc_in',soil_fsc_in_ptr)
        call read_tile_data_r1d_fptr(unit,'ssc_in',soil_ssc_in_ptr)
        if ((soil_carbon_option==SOILC_CORPSE .or. soil_carbon_option==SOILC_CORPSE_N) .and. nfu_inq_var(unit, 'deadmic_C_in')==NF_NOERR) then
           call read_tile_data_r1d_fptr(unit,'deadmic_C_in',soil_deadmic_C_in_ptr)

           call read_tile_data_r1d_fptr(unit,'fast_protected_C_in',soil_fast_protected_C_in_ptr)
           call read_tile_data_r1d_fptr(unit,'slow_protected_C_in', soil_slow_protected_C_in_ptr)
           call read_tile_data_r1d_fptr(unit,'deadmic_protected_C_in', soil_deadmic_protected_C_in_ptr)


           call read_tile_data_r0d_fptr(unit,'leaflitter_fsc_in', soil_leaflitter_fsc_in_ptr)
           call read_tile_data_r0d_fptr(unit,'leaflitter_ssc_in', soil_leaflitter_ssc_in_ptr)
           call read_tile_data_r0d_fptr(unit,'leaflitter_deadmic_C_in', soil_leaflitter_deadmic_C_in_ptr)


           call read_tile_data_r0d_fptr(unit,'finewoodlitter_fsc_in', soil_finewoodlitter_fsc_in_ptr)
           call read_tile_data_r0d_fptr(unit,'finewoodlitter_ssc_in', soil_finewoodlitter_ssc_in_ptr)
           call read_tile_data_r0d_fptr(unit,'finewoodlitter_deadmic_C_in', soil_finewoodlitter_deadmic_C_in_ptr)


           call read_tile_data_r0d_fptr(unit,'coarsewoodlitter_fsc_in', soil_coarsewoodlitter_fsc_in_ptr)
           call read_tile_data_r0d_fptr(unit,'coarsewoodlitter_ssc_in', soil_coarsewoodlitter_ssc_in_ptr)
           call read_tile_data_r0d_fptr(unit,'coarsewoodlitter_deadmic_C_in', soil_coarsewoodlitter_deadmic_C_in_ptr)


           call read_tile_data_r1d_fptr(unit,'fast_C_turnover_accumulated', soil_fast_C_turnover_accumulated_ptr)
           call read_tile_data_r1d_fptr(unit,'slow_C_turnover_accumulated', soil_slow_C_turnover_accumulated_ptr)
           call read_tile_data_r1d_fptr(unit,'deadmic_C_turnover_accumulated', soil_deadmic_C_turnover_accumulated_ptr)
           call read_tile_data_r1d_fptr(unit,'fast_protected_C_turnover_accumulated', soil_fast_protected_C_turnover_accumulated_ptr)
           call read_tile_data_r1d_fptr(unit,'slow_protected_C_turnover_accumulated', soil_slow_protected_C_turnover_accumulated_ptr)
           call read_tile_data_r1d_fptr(unit,'deadmic_protected_C_turnover_accumulated', soil_deadmic_protected_C_turnover_accumulated_ptr)


           call read_tile_data_r0d_fptr(unit,'leaflitter_fast_C_turnover_accumulated', soil_leaflitter_fast_C_turnover_accumulated_ptr)
           call read_tile_data_r0d_fptr(unit,'leaflitter_slow_C_turnover_accumulated', soil_leaflitter_slow_C_turnover_accumulated_ptr)
           call read_tile_data_r0d_fptr(unit,'leaflitter_deadmic_C_turnover_accumulated', soil_leaflitter_deadmic_C_turnover_accumulated_ptr)
           call read_tile_data_r0d_fptr(unit,'finewoodlitter_fast_C_turnover_accumulated', soil_finewoodlitter_fast_C_turnover_accumulated_ptr)
           call read_tile_data_r0d_fptr(unit,'finewoodlitter_slow_C_turnover_accumulated', soil_finewoodlitter_slow_C_turnover_accumulated_ptr)
           call read_tile_data_r0d_fptr(unit,'finewoodlitter_deadmic_C_turnover_accumulated', soil_finewoodlitter_deadmic_C_turnover_accumulated_ptr)
           call read_tile_data_r0d_fptr(unit,'coarsewoodlitter_fast_C_turnover_accumulated', soil_coarsewoodlitter_fast_C_turnover_accumulated_ptr)
           call read_tile_data_r0d_fptr(unit,'coarsewoodlitter_slow_C_turnover_accumulated', soil_coarsewoodlitter_slow_C_turnover_accumulated_ptr)
           call read_tile_data_r0d_fptr(unit,'coarsewoodlitter_deadmic_C_turnover_accumulated', soil_coarsewoodlitter_deadmic_C_turnover_accumulated_ptr)

IF(soil_carbon_option == SOILC_CORPSE_N) THEN
    call read_tile_data_r1d_fptr(unit,'fast_protected_N_turnover_accumulated', soil_fast_protected_N_turnover_accumulated_ptr)
    call read_tile_data_r1d_fptr(unit,'slow_protected_N_turnover_accumulated', soil_slow_protected_N_turnover_accumulated_ptr)
    call read_tile_data_r1d_fptr(unit,'deadmic_protected_N_turnover_accumulated', soil_deadmic_protected_N_turnover_accumulated_ptr)

    call read_tile_data_r1d_fptr(unit,'fast_N_turnover_accumulated', soil_fast_N_turnover_accumulated_ptr)
    call read_tile_data_r1d_fptr(unit,'slow_N_turnover_accumulated', soil_slow_N_turnover_accumulated_ptr)
    call read_tile_data_r1d_fptr(unit,'deadmic_N_turnover_accumulated', soil_deadmic_N_turnover_accumulated_ptr)

    call read_tile_data_r0d_fptr(unit,'coarsewoodlitter_fsn_in', soil_coarsewoodlitter_fsn_in_ptr)
    call read_tile_data_r0d_fptr(unit,'coarsewoodlitter_ssn_in', soil_coarsewoodlitter_ssn_in_ptr)
    call read_tile_data_r0d_fptr(unit,'coarsewoodlitter_deadmic_N_in', soil_coarsewoodlitter_deadmic_N_in_ptr)

    call read_tile_data_r0d_fptr(unit,'leaflitter_fsn_in', soil_leaflitter_fsn_in_ptr)
    call read_tile_data_r0d_fptr(unit,'leaflitter_ssn_in', soil_leaflitter_ssn_in_ptr)
    call read_tile_data_r0d_fptr(unit,'leaflitter_deadmic_N_in', soil_leaflitter_deadmic_N_in_ptr)
    call read_tile_data_r0d_fptr(unit,'finewoodlitter_fsn_in', soil_finewoodlitter_fsn_in_ptr)
    call read_tile_data_r0d_fptr(unit,'finewoodlitter_ssn_in', soil_finewoodlitter_ssn_in_ptr)
    call read_tile_data_r0d_fptr(unit,'finewoodlitter_deadmic_N_in', soil_finewoodlitter_deadmic_N_in_ptr)

    call read_tile_data_r1d_fptr(unit,'fsn_in',soil_fsn_in_ptr)
    call read_tile_data_r1d_fptr(unit,'ssn_in',soil_ssn_in_ptr)
    call read_tile_data_r1d_fptr(unit,'deadmic_N_in',soil_deadmic_N_in_ptr)
    call read_tile_data_r1d_fptr(unit,'fast_protected_N_in',soil_fast_protected_N_in_ptr)
    call read_tile_data_r1d_fptr(unit,'slow_protected_N_in', soil_slow_protected_N_in_ptr)
    call read_tile_data_r1d_fptr(unit,'deadmic_protected_N_in', soil_deadmic_protected_N_in_ptr)

           call read_tile_data_r0d_fptr(unit,'leaflitter_fast_N_turnover_accumulated', soil_leaflitter_fast_N_turnover_accumulated_ptr)
           call read_tile_data_r0d_fptr(unit,'leaflitter_slow_N_turnover_accumulated', soil_leaflitter_slow_N_turnover_accumulated_ptr)
           call read_tile_data_r0d_fptr(unit,'leaflitter_deadmic_N_turnover_accumulated', soil_leaflitter_deadmic_N_turnover_accumulated_ptr)
           call read_tile_data_r0d_fptr(unit,'finewoodlitter_fast_N_turnover_accumulated', soil_finewoodlitter_fast_N_turnover_accumulated_ptr)
           call read_tile_data_r0d_fptr(unit,'finewoodlitter_slow_N_turnover_accumulated', soil_finewoodlitter_slow_N_turnover_accumulated_ptr)
           call read_tile_data_r0d_fptr(unit,'finewoodlitter_deadmic_N_turnover_accumulated', soil_finewoodlitter_deadmic_N_turnover_accumulated_ptr)
           call read_tile_data_r0d_fptr(unit,'coarsewoodlitter_fast_N_turnover_accumulated', soil_coarsewoodlitter_fast_N_turnover_accumulated_ptr)
           call read_tile_data_r0d_fptr(unit,'coarsewoodlitter_slow_N_turnover_accumulated', soil_coarsewoodlitter_slow_N_turnover_accumulated_ptr)
           call read_tile_data_r0d_fptr(unit,'coarsewoodlitter_deadmic_N_turnover_accumulated', soil_coarsewoodlitter_deadmic_N_turnover_accumulated_ptr)
ENDIF

        endif
        __NF_ASRT__(nf_close(unit))
     endif
  endif

  ! ---- static diagnostic section
  call send_tile_data_r0d_fptr(id_tau_gw,       lnd%tile_map, soil_tau_groundwater_ptr)
  call send_tile_data_r0d_fptr(id_slope_l,      lnd%tile_map, soil_hillslope_length_ptr)
  call send_tile_data_r0d_fptr(id_slope_Z,      lnd%tile_map, soil_hillslope_relief_ptr)
  call send_tile_data_r0d_fptr(id_zeta_bar,     lnd%tile_map, soil_hillslope_zeta_bar_ptr)
  call send_tile_data_r0d_fptr(id_e_depth,      lnd%tile_map, soil_soil_e_depth_ptr)
  call send_tile_data_r0d_fptr(id_zeta,         lnd%tile_map, soil_zeta_ptr)
  call send_tile_data_r0d_fptr(id_tau,          lnd%tile_map, soil_tau_ptr)
  call send_tile_data_r0d_fptr(id_vwc_wilt,     lnd%tile_map, soil_vwc_wilt_ptr)
  call send_tile_data_r0d_fptr(id_vwc_fc,       lnd%tile_map, soil_vwc_fc_ptr)
  call send_tile_data_r0d_fptr(id_vwc_sat,      lnd%tile_map, soil_vwc_sat_ptr)
  call send_tile_data_r0d_fptr(id_K_sat,        lnd%tile_map, soil_k_sat_ref_ptr)
  call send_tile_data_r0d_fptr(id_K_gw,         lnd%tile_map, soil_k_sat_gw_ptr)
  call send_tile_data_r1d_fptr(id_w_fc,         lnd%tile_map, soil_w_fc_ptr)
  call send_tile_data_r1d_fptr(id_alpha,        lnd%tile_map, soil_alpha_ptr)
  call send_tile_data_r1d_fptr(id_refl_dry_dir, lnd%tile_map, soil_refl_dry_dir_ptr)
  call send_tile_data_r1d_fptr(id_refl_dry_dif, lnd%tile_map, soil_refl_dry_dif_ptr)
  call send_tile_data_r1d_fptr(id_refl_sat_dir, lnd%tile_map, soil_refl_sat_dir_ptr)
  call send_tile_data_r1d_fptr(id_refl_sat_dif, lnd%tile_map, soil_refl_sat_dif_ptr)
  call send_tile_data_r1d_fptr(id_f_iso_dry, lnd%tile_map, soil_f_iso_dry_ptr)
  call send_tile_data_r1d_fptr(id_f_vol_dry, lnd%tile_map, soil_f_vol_dry_ptr)
  call send_tile_data_r1d_fptr(id_f_geo_dry, lnd%tile_map, soil_f_geo_dry_ptr)
  call send_tile_data_r1d_fptr(id_f_iso_sat, lnd%tile_map, soil_f_iso_sat_ptr)
  call send_tile_data_r1d_fptr(id_f_vol_sat, lnd%tile_map, soil_f_vol_sat_ptr)
  call send_tile_data_r1d_fptr(id_f_geo_sat, lnd%tile_map, soil_f_geo_sat_ptr)
  call send_tile_data_i0d_fptr(id_type,         lnd%tile_map, soil_tag_ptr)
end subroutine soil_init


! ============================================================================
subroutine soil_diag_init ( id_lon, id_lat, id_band, id_zfull)
  integer, intent(in) :: id_lon  ! ID of land longitude (X) axis
  integer, intent(in) :: id_lat  ! ID of land latitude (Y) axis
  integer, intent(in) :: id_band ! ID of spectral band axis
  integer, intent(out) :: id_zfull ! ID of vertical soil axis

  ! ---- local vars
  integer :: axes(3)
  integer :: id_zhalf

  ! define vertical axis and its edges
  id_zhalf = diag_axis_init ( &
       'zhalf_soil', zhalf(1:num_l+1), 'meters', 'z', 'half level',  -1, set_name='soil' )
  id_zfull = diag_axis_init ( &
       'zfull_soil', zfull(1:num_l),   'meters', 'z', 'full level',  -1, set_name='soil', &
       edges=id_zhalf )

  ! define array of axis indices
  axes = (/ id_lon, id_lat, id_zfull /)

  ! set the default sub-sampling filter for the fields below
  call set_default_diag_filter('soil')

  ! define diagnostic fields
  id_rhizosphere_frac = register_tiled_diag_field ( module_name, 'rhizosphere_frac', axes, &
       land_time, 'rhizosphere fraction of soil', 'm3/m3', missing_value=-100.0 )
  id_fast_soil_C = register_tiled_diag_field ( module_name, 'fast_soil_C', axes,  &
       land_time, 'fast soil carbon', 'kg C/m3', missing_value=-100.0 )
  id_slow_soil_C = register_tiled_diag_field ( module_name, 'slow_soil_C', axes,  &
       land_time, 'slow soil carbon', 'kg C/m3', missing_value=-100.0 )
  id_deadmic_C = register_tiled_diag_field ( module_name, 'dead_microbe_C', axes,  &
       land_time, 'dead microbe soil carbon per layer', 'kg C/m3', missing_value=-100.0 )
  id_protectedC = register_tiled_diag_field ( module_name, 'protected_soil_C', axes,  &
       land_time, 'protected soil carbon per layer', 'kg C/m3', missing_value=-100.0 )
  id_livemic_C = register_tiled_diag_field ( module_name, 'livemic',  &
       (/id_lon,id_lat,id_zfull/), land_time, 'live microbe soil carbon', 'kg C/m3', &
       missing_value=-100.0 )
   id_fast_soil_N = register_tiled_diag_field ( module_name, 'fast_soil_N', axes,  &
        land_time, 'fast soil nitrogen', 'kg N/m3', missing_value=-100.0 )
   id_slow_soil_N = register_tiled_diag_field ( module_name, 'slow_soil_N', axes,  &
        land_time, 'slow soil nitrogen', 'kg N/m3', missing_value=-100.0 )
   id_deadmic_N = register_tiled_diag_field ( module_name, 'dead_microbe_N', axes,  &
        land_time, 'dead microbe soil nitrogen per layer', 'kg N/m3', missing_value=-100.0 )
   id_protectedN = register_tiled_diag_field ( module_name, 'protected_soil_N', axes,  &
        land_time, 'protected soil nitrogen per layer', 'kg N/m3', missing_value=-100.0 )
   id_livemic_N = register_tiled_diag_field ( module_name, 'live_microbe_N',  &
        (/id_lon,id_lat,id_zfull/), land_time, 'live microbe soil nitrogen', 'kg N/m3', &
        missing_value=-100.0 )
  id_fast_dissolved_C = register_tiled_diag_field ( module_name, 'fast_dissolved_C', axes,  &
       land_time, 'fast dissolved carbon per layer', 'kg C/m3', missing_value=-100.0 )
  id_slow_dissolved_C = register_tiled_diag_field ( module_name, 'slow_dissolved_C', axes,  &
       land_time, 'slow dissolved carbon per layer', 'kg C/m3', missing_value=-100.0 )
  id_deadmic_dissolved_C = register_tiled_diag_field ( module_name, 'dead_microbe_dissolved_C', axes,  &
       land_time, 'dead microbe dissolved carbon per layer', 'kg C/m3', missing_value=-100.0 )
   id_fast_dissolved_N = register_tiled_diag_field ( module_name, 'fast_dissolved_N', axes,  &
        land_time, 'fast dissolved nitrogen per layer', 'kg N/m3', missing_value=-100.0 )
   id_slow_dissolved_N = register_tiled_diag_field ( module_name, 'slow_dissolved_N', axes,  &
        land_time, 'slow dissolved nitrogen per layer', 'kg N/m3', missing_value=-100.0 )
   id_deadmic_dissolved_N = register_tiled_diag_field ( module_name, 'dead_microbe_dissolved_N', axes,  &
        land_time, 'dead microbe dissolved nitrogen per layer', 'kg N/m3', missing_value=-100.0 )
    id_soil_nitrate = register_tiled_diag_field ( module_name, 'soil_NO3', axes,  &
         land_time, 'NO3 per layer', 'kg N/m3', missing_value=-100.0 )
     id_soil_ammonium = register_tiled_diag_field ( module_name, 'soil_NH4', axes,  &
          land_time, 'NO3 per layer', 'kg N/m3', missing_value=-100.0 )
      id_nitrif = register_tiled_diag_field ( module_name, 'soil_nitrif', axes,  &
           land_time, 'Cumulative nitrification per layer', 'kg N/m3', missing_value=-100.0 )
   id_denitrif = register_tiled_diag_field ( module_name, 'soil_denitrif', axes,  &
        land_time, 'Cumulative denitrification per layer', 'kg N/m3', missing_value=-100.0 )
      id_nitrification_rate = register_tiled_diag_field ( module_name, 'soil_nitrification_rate', axes,  &
           land_time, 'Nitrification rate per layer', 'kg N/m3/year', missing_value=-100.0 )
   id_denitrification_rate = register_tiled_diag_field ( module_name, 'soil_denitrification_rate', axes,  &
        land_time, 'Denitrification rate per layer', 'kg N/m3/year', missing_value=-100.0 )
  id_total_carbon_layered = register_tiled_diag_field ( module_name, 'total_soil_carbon_layered', axes,  &
       land_time, 'total soil carbon per layer', 'kg C/m3', missing_value=-100.0 )
   id_total_nitrogen_layered = register_tiled_diag_field ( module_name, 'total_soil_nitrogen_layered', axes,  &
        land_time, 'total soil nitrogen per layer', 'kg N/m3', missing_value=-100.0 )
  id_fast_DOC_div_loss = register_tiled_diag_field ( module_name, 'fast_DOC_div_loss', (/id_lon,id_lat/),  &
       land_time, 'total fast DOC divergence loss', 'kg C/m2', missing_value=-100.0 )
  id_slow_DOC_div_loss = register_tiled_diag_field ( module_name, 'slow_DOC_div_loss', (/id_lon,id_lat/),  &
       land_time, 'total slow DOC divergence loss', 'kg C/m2', missing_value=-100.0 )
  id_deadmic_DOC_div_loss = register_tiled_diag_field ( module_name, 'deadmic_DOC_div_loss', (/id_lon,id_lat/),  &
       land_time, 'total dead microbe DOC divergence loss', 'kg C/m2', missing_value=-100.0 )
  id_total_DOC_div_loss = register_tiled_diag_field ( module_name, 'total_DOC_div', axes(1:2), &
       land_time, 'total rate of DOC divergence loss', 'kg C/m^2/s', missing_value=initval)

       id_fast_DON_div_loss = register_tiled_diag_field ( module_name, 'fast_DON_div_loss', (/id_lon,id_lat/),  &
            land_time, 'total fast DON divergence loss', 'kg N/m2', missing_value=-100.0 )
       id_slow_DON_div_loss = register_tiled_diag_field ( module_name, 'slow_DON_div_loss', (/id_lon,id_lat/),  &
            land_time, 'total slow DON divergence loss', 'kg N/m2', missing_value=-100.0 )
       id_deadmic_DON_div_loss = register_tiled_diag_field ( module_name, 'deadmic_DON_div_loss', (/id_lon,id_lat/),  &
            land_time, 'total dead microbe DON divergence loss', 'kg N/m2', missing_value=-100.0 )
       id_total_DON_div_loss = register_tiled_diag_field ( module_name, 'total_DON_div', axes(1:2), &
            land_time, 'total rate of DON divergence loss', 'kg N/m^2/s', missing_value=initval)
    id_NO3_div_loss = register_tiled_diag_field ( module_name, 'total_NO3_div', axes(1:2), &
         land_time, 'total rate of NO3 divergence loss', 'kg N/m^2/s', missing_value=initval)
     id_NH4_div_loss = register_tiled_diag_field ( module_name, 'total_NH4_div', axes(1:2), &
          land_time, 'total rate of NH4 divergence loss', 'kg N/m^2/s', missing_value=initval)

  id_rsoil = register_tiled_diag_field ( module_name, 'rsoil',  &
       (/id_lon,id_lat/), land_time, 'soil respiration', 'kg C/(m2 year)', &
       missing_value=-100.0 )
  id_rsoil_fast = register_tiled_diag_field ( module_name, 'rsoil_fast',  &
       axes, land_time, 'fast soil carbon respiration', 'kg C/(m3 year)', &
       missing_value=-100.0 )
  id_rsoil_fast_N = register_tiled_diag_field ( module_name, 'rsoil_fast_N',  &
       axes, land_time, 'fast soil nitrogen degradation', 'kg N/(m3 year)', &
       missing_value=-100.0 )
  id_rsoil_slow = register_tiled_diag_field ( module_name, 'rsoil_slow',  &
       axes, land_time, 'slow soil carbon respiration', 'kg C/(m3 year)', &
       missing_value=-100.0 )
  id_rsoil_slow_N = register_tiled_diag_field ( module_name, 'rsoil_slow_N',  &
       axes, land_time, 'slow soil nitrogen degradation', 'kg N/(m3 year)', &
       missing_value=-100.0 )
  id_rsoil_deadmic = register_tiled_diag_field ( module_name, 'rsoil_deadmic',  &
       axes, land_time, 'dead microbe soil carbon respiration', 'kg C/(m3 year)', &
       missing_value=-100.0 )
  id_rsoil_deadmic_N = register_tiled_diag_field ( module_name, 'rsoil_deadmic_N',  &
       axes, land_time, 'dead microbe soil nitrogen degradation', 'kg N/(m3 year)', &
       missing_value=-100.0 )
  id_N_mineralization_rate = register_tiled_diag_field ( module_name, 'N_mineralization_rate',  &
       axes, land_time, 'Soil N mineralization rate', 'kg N/(m3 year)', &
       missing_value=-100.0 )
  id_N_immobilization_rate = register_tiled_diag_field ( module_name, 'N_immobilization_rate',  &
       axes, land_time, 'Soil N immobilization rate', 'kg N/(m3 year)', &
       missing_value=-100.0 )
  id_rsoil_leaflitter_fast = register_tiled_diag_field ( module_name, 'rsoil_leaflitter_fast',  &
       (/id_lon,id_lat/), land_time, 'surface leaf litter fast C respiration', 'kg C/(m2 year)', &
       missing_value=-100.0 )
  id_rsoil_leaflitter_slow = register_tiled_diag_field ( module_name, 'rsoil_leaflitter_slow',  &
       (/id_lon,id_lat/), land_time, 'surface leaf litter slow C respiration', 'kg C/(m2 year)', &
       missing_value=-100.0 )
  id_rsoil_leaflitter_deadmic = register_tiled_diag_field ( module_name, 'rsoil_leaflitter_deadmic',  &
       (/id_lon,id_lat/), land_time, 'surface leaf litter dead microbe C respiration', 'kg C/(m2 year)', &
       missing_value=-100.0 )
  id_rsoil_leaflitter_fast_N = register_tiled_diag_field ( module_name, 'rsoil_leaflitter_fast_N',  &
       (/id_lon,id_lat/), land_time, 'surface leaf litter fast N degradation', 'kg N/(m2 year)', &
       missing_value=-100.0 )
  id_rsoil_leaflitter_slow_N = register_tiled_diag_field ( module_name, 'rsoil_leaflitter_slow_N',  &
       (/id_lon,id_lat/), land_time, 'surface leaf litter slow N degradation', 'kg N/(m2 year)', &
       missing_value=-100.0 )
  id_rsoil_leaflitter_deadmic_N = register_tiled_diag_field ( module_name, 'rsoil_leaflitter_deadmic_N',  &
       (/id_lon,id_lat/), land_time, 'surface leaf litter dead microbe N degradation', 'kg N/(m2 year)', &
       missing_value=-100.0 )
  id_leaflitter_N_mineralization_rate = register_tiled_diag_field ( module_name, 'leaflitter_N_mineralization_rate',  &
       (/id_lon,id_lat/), land_time, 'surface leaf litter N mineralization', 'kg N/(m2 year)', &
       missing_value=-100.0 )
  id_leaflitter_N_immobilization_rate = register_tiled_diag_field ( module_name, 'leaflitter_N_immobilization_rate',  &
       (/id_lon,id_lat/), land_time, 'surface leaf litter N immobilization', 'kg N/(m2 year)', &
       missing_value=-100.0 )
  id_leaflitter_nitrification_rate = register_tiled_diag_field ( module_name, 'leaflitter_nitrification_rate',  &
       (/id_lon,id_lat/), land_time, 'surface leaf litter nitrification', 'kg N/(m2 year)', &
       missing_value=-100.0 )
  id_leaflitter_denitrification_rate = register_tiled_diag_field ( module_name, 'leaflitter_denitrification_rate',  &
       (/id_lon,id_lat/), land_time, 'surface leaf litter denitrification', 'kg N/(m2 year)', &
       missing_value=-100.0 )
  id_rsoil_coarsewoodlitter_fast = register_tiled_diag_field ( module_name, 'rsoil_coarsewoodlitter_fast',  &
       (/id_lon,id_lat/), land_time, 'surface coarse wood litter fast C respiration', 'kg C/(m2 year)', &
       missing_value=-100.0 )
  id_rsoil_coarsewoodlitter_slow = register_tiled_diag_field ( module_name, 'rsoil_coarsewoodlitter_slow',  &
       (/id_lon,id_lat/), land_time, 'surface coarse wood litter slow C respiration', 'kg C/(m2 year)', &
       missing_value=-100.0 )
  id_rsoil_coarsewoodlitter_deadmic = register_tiled_diag_field ( module_name, 'rsoil_coarsewoodlitter_deadmic',  &
       (/id_lon,id_lat/), land_time, 'surface coarse wood litter dead microbe C respiration', 'kg C/(m2 year)', &
       missing_value=-100.0 )
  id_rsoil_coarsewoodlitter_fast_N = register_tiled_diag_field ( module_name, 'rsoil_coarsewoodlitter_fast_N',  &
       (/id_lon,id_lat/), land_time, 'surface coarse wood litter fast N degradation', 'kg N/(m2 year)', &
       missing_value=-100.0 )
  id_rsoil_coarsewoodlitter_slow_N = register_tiled_diag_field ( module_name, 'rsoil_coarsewoodlitter_slow_N',  &
       (/id_lon,id_lat/), land_time, 'surface coarse wood litter slow N degradation', 'kg N/(m2 year)', &
       missing_value=-100.0 )
  id_rsoil_coarsewoodlitter_deadmic_N = register_tiled_diag_field ( module_name, 'rsoil_coarsewoodlitter_deadmic_N',  &
       (/id_lon,id_lat/), land_time, 'surface coarse wood litter dead microbe N degradation', 'kg N/(m2 year)', &
       missing_value=-100.0 )
  id_coarsewoodlitter_N_mineralization_rate = register_tiled_diag_field ( module_name, 'coarsewoodlitter_N_mineralization_rate',  &
       (/id_lon,id_lat/), land_time, 'surface coarse wood litter N mineralization', 'kg N/(m2 year)', &
       missing_value=-100.0 )
  id_coarsewoodlitter_N_immobilization_rate = register_tiled_diag_field ( module_name, 'coarsewoodlitter_N_immobilization_rate',  &
       (/id_lon,id_lat/), land_time, 'surface coarse wood litter N immobilization', 'kg N/(m2 year)', &
       missing_value=-100.0 )
  id_coarsewoodlitter_nitrification_rate = register_tiled_diag_field ( module_name, 'coarsewoodlitter_nitrification_rate',  &
       (/id_lon,id_lat/), land_time, 'surface coarse wood litter nitrification', 'kg N/(m2 year)', &
       missing_value=-100.0 )
  id_coarsewoodlitter_denitrification_rate = register_tiled_diag_field ( module_name, 'coarsewoodlitter_denitrification_rate',  &
       (/id_lon,id_lat/), land_time, 'surface coarse wood litter denitrification', 'kg N/(m2 year)', &
       missing_value=-100.0 )
  id_rsoil_finewoodlitter_fast = register_tiled_diag_field ( module_name, 'rsoil_finewoodlitter_fast',  &
       (/id_lon,id_lat/), land_time, 'surface fine wood litter fast C respiration', 'kg C/(m2 year)', &
       missing_value=-100.0 )
  id_rsoil_finewoodlitter_slow = register_tiled_diag_field ( module_name, 'rsoil_finewoodlitter_slow',  &
       (/id_lon,id_lat/), land_time, 'surface fine wood litter slow C respiration', 'kg C/(m2 year)', &
       missing_value=-100.0 )
  id_rsoil_finewoodlitter_deadmic = register_tiled_diag_field ( module_name, 'rsoil_finewoodlitter_deadmic',  &
       (/id_lon,id_lat/), land_time, 'surface fine wood litter dead microbe C respiration', 'kg C/(m2 year)', &
       missing_value=-100.0 )
  id_rsoil_finewoodlitter_fast_N = register_tiled_diag_field ( module_name, 'rsoil_finewoodlitter_fast_N',  &
       (/id_lon,id_lat/), land_time, 'surface fine wood litter fast N degradation', 'kg N/(m2 year)', &
       missing_value=-100.0 )
  id_rsoil_finewoodlitter_slow_N = register_tiled_diag_field ( module_name, 'rsoil_finewoodlitter_slow_N',  &
       (/id_lon,id_lat/), land_time, 'surface fine wood litter slow N degradation', 'kg N/(m2 year)', &
       missing_value=-100.0 )
  id_rsoil_finewoodlitter_deadmic_N = register_tiled_diag_field ( module_name, 'rsoil_finewoodlitter_deadmic_N',  &
       (/id_lon,id_lat/), land_time, 'surface fine wood litter dead microbe N degradation', 'kg N/(m2 year)', &
       missing_value=-100.0 )
  id_finewoodlitter_N_mineralization_rate = register_tiled_diag_field ( module_name, 'finewoodlitter_N_mineralization_rate',  &
       (/id_lon,id_lat/), land_time, 'surface fine wood litter N mineralization', 'kg N/(m2 year)', &
       missing_value=-100.0 )
  id_finewoodlitter_N_immobilization_rate = register_tiled_diag_field ( module_name, 'finewoodlitter_N_immobilization_rate',  &
       (/id_lon,id_lat/), land_time, 'surface fine wood litter N immobilization', 'kg N/(m2 year)', &
       missing_value=-100.0 )
  id_finewoodlitter_nitrification_rate = register_tiled_diag_field ( module_name, 'finewoodlitter_nitrification_rate',  &
       (/id_lon,id_lat/), land_time, 'surface fine wood litter nitrification', 'kg N/(m2 year)', &
       missing_value=-100.0 )
  id_finewoodlitter_denitrification_rate = register_tiled_diag_field ( module_name, 'finewoodlitter_denitrification_rate',  &
       (/id_lon,id_lat/), land_time, 'surface fine wood litter denitrification', 'kg N/(m2 year)', &
       missing_value=-100.0 )
  id_total_N_mineralization_rate = register_tiled_diag_field ( module_name, 'total_N_mineralization_rate',  &
       (/id_lon,id_lat/), land_time, 'Total N mineralization', 'kg N/(m2 year)', &
       missing_value=-100.0 )
  id_total_N_immobilization_rate = register_tiled_diag_field ( module_name, 'total_N_immobilization_rate',  &
       (/id_lon,id_lat/), land_time, 'Total N immobilization', 'kg N/(m2 year)', &
       missing_value=-100.0 )
  id_total_nitrification_rate = register_tiled_diag_field ( module_name, 'total_nitrification_rate',  &
       (/id_lon,id_lat/), land_time, 'Total nitrification', 'kg N/(m2 year)', &
       missing_value=-100.0 )
  id_total_denitrification_rate = register_tiled_diag_field ( module_name, 'total_denitrification_rate',  &
       (/id_lon,id_lat/), land_time, 'Total denitrification', 'kg N/(m2 year)', &
       missing_value=-100.0 )
  id_C_dissolve_rate_fast = register_tiled_diag_field ( module_name, 'fast_C_dissolve_rate',  &
       axes, land_time, 'fast soil carbon dissolving rate', 'kg C/(m3 year)', &
       missing_value=-100.0 )
  id_C_dissolve_rate_slow = register_tiled_diag_field ( module_name, 'slow_C_dissolve_rate',  &
       axes, land_time, 'slow soil carbon dissolving rate', 'kg C/(m3 year)', &
       missing_value=-100.0 )
  id_C_dissolve_rate_deadmic = register_tiled_diag_field ( module_name, 'deadmic_C_dissolve_rate',  &
       axes, land_time, 'dead microbe soil carbon dissolving rate', 'kg C/(m3 year)', &
       missing_value=-100.0 )
  id_C_deposition_fast = register_tiled_diag_field ( module_name, 'fast_C_deposition_rate',  &
       axes, land_time, 'fast soil carbon deposition from DOC rate', 'kg C/(m3 year)', &
       missing_value=-100.0 )
  id_C_deposition_slow = register_tiled_diag_field ( module_name, 'slow_C_deposition_rate',  &
       axes, land_time, 'slow soil carbon deposition from DOC rate', 'kg C/(m3 year)', &
       missing_value=-100.0 )
  id_C_deposition_deadmic = register_tiled_diag_field ( module_name, 'deadmic_C_deposition_rate',  &
       axes, land_time, 'dead microbe soil carbon deposition from DOC rate', 'kg C/(m3 year)', &
       missing_value=-100.0 )
  id_leaflitter_C_dissolve_rate_fast = register_tiled_diag_field ( module_name, 'leaflitter_fast_C_dissolve_rate',  &
       axes(1:2), land_time, 'fast leaf litter carbon dissolving rate', 'kg C/(m2 year)', &
       missing_value=-100.0 )
  id_leaflitter_C_dissolve_rate_slow = register_tiled_diag_field ( module_name, 'leaflitter_slow_C_dissolve_rate',  &
       axes(1:2), land_time, 'slow leaf litter carbon dissolving rate', 'kg C/(m2 year)', &
       missing_value=-100.0 )
  id_leaflitter_C_dissolve_rate_deadmic = register_tiled_diag_field ( module_name, 'leaflitter_deadmic_C_dissolve_rate',  &
       axes(1:2), land_time, 'dead microbe leaf litter carbon dissolving rate', 'kg C/(m2 year)', &
       missing_value=-100.0 )
  id_leaflitter_C_deposition_fast = register_tiled_diag_field ( module_name, 'leaflitter_fast_C_deposition_rate',  &
       axes(1:2), land_time, 'fast leaf litter carbon deposition from DOC rate', 'kg C/(m2 year)', &
       missing_value=-100.0 )
  id_leaflitter_C_deposition_slow = register_tiled_diag_field ( module_name, 'leaflitter_slow_C_deposition_rate',  &
       axes(1:2), land_time, 'slow leaf litter carbon deposition from DOC rate', 'kg C/(m2 year)', &
       missing_value=-100.0 )
  id_leaflitter_C_deposition_deadmic = register_tiled_diag_field ( module_name, 'leaflitter_deadmic_C_deposition_rate',  &
       axes(1:2), land_time, 'dead microbe leaf litter carbon deposition from DOC rate', 'kg C/(m2 year)', &
       missing_value=-100.0 )
  id_finewoodlitter_C_dissolve_rate_fast = register_tiled_diag_field ( module_name, 'finewoodlitter_fast_C_dissolve_rate',  &
       axes(1:2), land_time, 'fast fine wood litter carbon dissolving rate', 'kg C/(m2 year)', &
       missing_value=-100.0 )
  id_finewoodlitter_C_dissolve_rate_slow = register_tiled_diag_field ( module_name, 'finewoodlitter_slow_C_dissolve_rate',  &
       axes(1:2), land_time, 'slow fine wood litter carbon dissolving rate', 'kg C/(m2 year)', &
       missing_value=-100.0 )
  id_finewoodlitter_C_dissolve_rate_deadmic = register_tiled_diag_field ( module_name, 'finewoodlitter_deadmic_C_dissolve_rate',  &
       axes(1:2), land_time, 'dead microbe fine wood litter carbon dissolving rate', 'kg C/(m2 year)', &
       missing_value=-100.0 )
  id_finewoodlitter_C_deposition_fast = register_tiled_diag_field ( module_name, 'finewoodlitter_fast_C_deposition_rate',  &
       axes(1:2), land_time, 'fast fine wood litter carbon deposition from DOC rate', 'kg C/(m2 year)', &
       missing_value=-100.0 )
  id_finewoodlitter_C_deposition_slow = register_tiled_diag_field ( module_name, 'finewoodlitter_slow_C_deposition_rate',  &
       axes(1:2), land_time, 'slow fine wood litter carbon deposition from DOC rate', 'kg C/(m2 year)', &
       missing_value=-100.0 )
  id_finewoodlitter_C_deposition_deadmic = register_tiled_diag_field ( module_name, 'finewoodlitter_deadmic_C_deposition_rate',  &
       axes(1:2), land_time, 'dead microbe fine wood litter carbon deposition from DOC rate', 'kg C/(m2 year)', &
       missing_value=-100.0 )
    id_coarsewoodlitter_C_dissolve_rate_fast = register_tiled_diag_field ( module_name, 'coarsewoodlitter_fast_C_dissolve_rate',  &
       axes(1:2), land_time, 'fast coarse wood litter carbon dissolving rate', 'kg C/(m2 year)', &
       missing_value=-100.0 )
  id_coarsewoodlitter_C_dissolve_rate_slow = register_tiled_diag_field ( module_name, 'coarsewoodlitter_slow_C_dissolve_rate',  &
       axes(1:2), land_time, 'slow coarse wood litter carbon dissolving rate', 'kg C/(m2 year)', &
       missing_value=-100.0 )
  id_coarsewoodlitter_C_dissolve_rate_deadmic = register_tiled_diag_field ( module_name, 'coarsewoodlitter_deadmic_C_dissolve_rate',  &
       axes(1:2), land_time, 'dead microbe coarse wood litter carbon dissolving rate', 'kg C/(m2 year)', &
       missing_value=-100.0 )
  id_coarsewoodlitter_C_deposition_fast = register_tiled_diag_field ( module_name, 'coarsewoodlitter_fast_C_deposition_rate',  &
       axes(1:2), land_time, 'fast coarse wood litter carbon deposition from DOC rate', 'kg C/(m2 year)', &
       missing_value=-100.0 )
  id_coarsewoodlitter_C_deposition_slow = register_tiled_diag_field ( module_name, 'coarsewoodlitter_slow_C_deposition_rate',  &
       axes(1:2), land_time, 'slow coarse wood litter carbon deposition from DOC rate', 'kg C/(m2 year)', &
       missing_value=-100.0 )
  id_coarsewoodlitter_C_deposition_deadmic = register_tiled_diag_field ( module_name, 'coarsewoodlitter_deadmic_C_deposition_rate',  &
       axes(1:2), land_time, 'dead microbe coarse wood litter carbon deposition from DOC rate', 'kg C/(m2 year)', &
       missing_value=-100.0 )
  id_resp = register_tiled_diag_field ( module_name, 'resp', (/id_lon,id_lat/), &
       land_time, 'Total soil respiration', 'kg C/(m2 year)', missing_value=-100.0 )
  id_Qmax = register_tiled_diag_field ( module_name, 'Qmax', axes(1:2),  &
       land_time, 'Maximum clay sorptive capacity', 'kg C/m3', missing_value=-100.0 )

  id_leaflitter_fast_C = register_tiled_diag_field ( module_name, 'fast_leaflitter_C', axes(1:2),  &
       land_time, 'fast leaf litter carbon', 'kg C/m2', missing_value=-100.0 )
  id_leaflitter_slow_C = register_tiled_diag_field ( module_name, 'slow_leaflitter_C', axes(1:2),  &
       land_time, 'slow leaf litter carbon', 'kg C/m2', missing_value=-100.0 )
  id_leaflitter_deadmic_C = register_tiled_diag_field ( module_name, 'leaflitter_dead_microbe_C', axes(1:2),  &
       land_time, 'dead microbe leaf litter carbon', 'kg C/m2', missing_value=-100.0 )
  id_leaflitter_livemic_C = register_tiled_diag_field ( module_name, 'leaflitter_live_microbe_C', axes(1:2),  &
       land_time, 'live microbe leaf litter carbon', 'kg C/m2', missing_value=-100.0 )
       id_leaflitter_total_C = register_tiled_diag_field ( module_name, 'leaflitter_total_C', axes(1:2),  &
            land_time, 'leaf litter total carbon', 'kg C/m2', missing_value=-100.0 )
            id_leaflitter_total_N = register_tiled_diag_field ( module_name, 'leaflitter_total_N', axes(1:2),  &
                 land_time, 'leaf litter total nitrogen', 'kg N/m2', missing_value=-100.0 )
       id_leaflitter_fast_N = register_tiled_diag_field ( module_name, 'fast_leaflitter_N', axes(1:2),  &
            land_time, 'fast leaf litter nitrogen', 'kg N/m2', missing_value=-100.0 )
       id_leaflitter_slow_N = register_tiled_diag_field ( module_name, 'slow_leaflitter_N', axes(1:2),  &
            land_time, 'slow leaf litter nitrogen', 'kg N/m2', missing_value=-100.0 )
       id_leaflitter_deadmic_N = register_tiled_diag_field ( module_name, 'leaflitter_dead_microbe_N', axes(1:2),  &
            land_time, 'dead microbe leaf litter nitrogen', 'kg N/m2', missing_value=-100.0 )
       id_leaflitter_livemic_N = register_tiled_diag_field ( module_name, 'leaflitter_live_microbe_N', axes(1:2),  &
            land_time, 'live microbe leaf litter nitrogen', 'kg N/m2', missing_value=-100.0 )
        id_leaflitter_nitrate = register_tiled_diag_field ( module_name, 'leaflitter_NO3', axes(1:2),  &
             land_time, 'Leaf litter nitrate', 'kg N/m2', missing_value=-100.0 )
     id_leaflitter_ammonium = register_tiled_diag_field ( module_name, 'leaflitter_NH4', axes(1:2),  &
          land_time, 'Leaf litter ammonium', 'kg N/m2', missing_value=-100.0 )

  id_coarsewoodlitter_fast_C = register_tiled_diag_field ( module_name, 'fast_coarsewoodlitter_C', axes(1:2),  &
       land_time, 'fast coarse wood litter carbon', 'kg C/m2', missing_value=-100.0 )
  id_coarsewoodlitter_slow_C = register_tiled_diag_field ( module_name, 'slow_coarsewoodlitter_C', axes(1:2),  &
       land_time, 'slow coarse wood litter carbon', 'kg C/m2', missing_value=-100.0 )
  id_coarsewoodlitter_deadmic_C = register_tiled_diag_field ( module_name, 'coarsewoodlitter_dead_microbe_C', axes(1:2),  &
       land_time, 'dead microbe coarse wood litter carbon', 'kg C/m2', missing_value=-100.0 )
  id_coarsewoodlitter_livemic_C = register_tiled_diag_field ( module_name, 'coarsewoodlitter_live_microbe_C', axes(1:2),  &
       land_time, 'live microbe coarse wood litter carbon', 'kg C/m2', missing_value=-100.0 )
  id_coarseWoodlitter_total_C = register_tiled_diag_field ( module_name, 'coarsewoodlitter_total_C', axes(1:2),  &
       land_time, 'coarse wood litter total carbon', 'kg C/m2', missing_value=-100.0 )

       id_coarseWoodlitter_total_N = register_tiled_diag_field ( module_name, 'coarsewoodlitter_total_N', axes(1:2),  &
            land_time, 'coarse wood litter total nitrogen', 'kg N/m2', missing_value=-100.0 )

       id_coarsewoodlitter_fast_N = register_tiled_diag_field ( module_name, 'fast_coarsewoodlitter_N', axes(1:2),  &
            land_time, 'fast coarse wood litter nitrogen', 'kg N/m2', missing_value=-100.0 )
       id_coarsewoodlitter_slow_N = register_tiled_diag_field ( module_name, 'slow_coarsewoodlitter_N', axes(1:2),  &
            land_time, 'slow coarse wood litter nitrogen', 'kg N/m2', missing_value=-100.0 )
       id_coarsewoodlitter_deadmic_N = register_tiled_diag_field ( module_name, 'coarsewoodlitter_dead_microbe_N', axes(1:2),  &
            land_time, 'dead microbe coarse wood litter nitrogen', 'kg N/m2', missing_value=-100.0 )
       id_coarsewoodlitter_livemic_N = register_tiled_diag_field ( module_name, 'coarsewoodlitter_live_microbe_N', axes(1:2),  &
            land_time, 'live microbe coarse wood litter nitrogen', 'kg N/m2', missing_value=-100.0 )
            id_coarsewoodlitter_nitrate = register_tiled_diag_field ( module_name, 'coarsewoodlitter_NO3', axes(1:2),  &
                 land_time, 'Coarse wood litter nitrate', 'kg N/m2', missing_value=-100.0 )
         id_coarsewoodlitter_ammonium = register_tiled_diag_field ( module_name, 'coarsewoodlitter_NH4', axes(1:2),  &
              land_time, 'Coarse wood litter ammonium', 'kg N/m2', missing_value=-100.0 )

  id_finewoodlitter_fast_C = register_tiled_diag_field ( module_name, 'fast_finewoodlitter_C', axes(1:2),  &
       land_time, 'fast fine wood litter carbon', 'kg C/m2', missing_value=-100.0 )
  id_finewoodlitter_slow_C = register_tiled_diag_field ( module_name, 'slow_finewoodlitter_C', axes(1:2),  &
       land_time, 'slow fine wood litter carbon', 'kg C/m2', missing_value=-100.0 )
  id_finewoodlitter_deadmic_C = register_tiled_diag_field ( module_name, 'finewoodlitter_dead_microbe_C', axes(1:2),  &
       land_time, 'dead microbe fine wood litter carbon', 'kg C/m2', missing_value=-100.0 )
  id_finewoodlitter_livemic_C = register_tiled_diag_field ( module_name, 'finewoodlitter_live_microbe_C', axes(1:2),  &
       land_time, 'live microbe fine wood litter carbon', 'kg C/m2', missing_value=-100.0 )
       id_fineWoodlitter_total_C = register_tiled_diag_field ( module_name, 'finewoodlitter_total_C', axes(1:2),  &
            land_time, 'fine wood litter total carbon', 'kg C/m2', missing_value=-100.0 )
      id_fineWoodlitter_total_N = register_tiled_diag_field ( module_name, 'finewoodlitter_total_N', axes(1:2),  &
           land_time, 'fine wood litter total nitrogen', 'kg N/m2', missing_value=-100.0 )
       id_finewoodlitter_fast_N = register_tiled_diag_field ( module_name, 'fast_finewoodlitter_N', axes(1:2),  &
            land_time, 'fast fine wood litter nitrogen', 'kg N/m2', missing_value=-100.0 )
       id_finewoodlitter_slow_N = register_tiled_diag_field ( module_name, 'slow_finewoodlitter_N', axes(1:2),  &
            land_time, 'slow fine wood litter nitrogen', 'kg N/m2', missing_value=-100.0 )
       id_finewoodlitter_deadmic_N = register_tiled_diag_field ( module_name, 'finewoodlitter_dead_microbe_N', axes(1:2),  &
            land_time, 'dead microbe fine wood litter nitrogen', 'kg N/m2', missing_value=-100.0 )
       id_finewoodlitter_livemic_N = register_tiled_diag_field ( module_name, 'finewoodlitter_live_microbe_N', axes(1:2),  &
            land_time, 'live microbe fine wood litter nitrogen', 'kg N/m2', missing_value=-100.0 )
            id_finewoodlitter_nitrate = register_tiled_diag_field ( module_name, 'finewoodlitter_NO3', axes(1:2),  &
                 land_time, 'Fine wood litter nitrate', 'kg N/m2', missing_value=-100.0 )
         id_finewoodlitter_ammonium = register_tiled_diag_field ( module_name, 'finewoodlitter_NH4', axes(1:2),  &
              land_time, 'Fine wood litter ammonium', 'kg N/m2', missing_value=-100.0 )


  id_leaflitter_fast_dissolved_C = register_tiled_diag_field ( module_name, 'fast_leaflitter_dissolved_C', axes(1:2),  &
       land_time, 'fast leaf litter dissolved carbon', 'kg C/m2', missing_value=-100.0 )
  id_leaflitter_slow_dissolved_C = register_tiled_diag_field ( module_name, 'slow_leaflitter_dissolved_C', axes(1:2),  &
       land_time, 'slow leaf litter dissolved carbon', 'kg C/m2', missing_value=-100.0 )
  id_leaflitter_deadmic_dissolved_C = register_tiled_diag_field ( module_name, 'leaflitter_dead_microbe_dissolved_C', axes(1:2),  &
       land_time, 'dead microbe leaf litter dissolved carbon', 'kg C/m2', missing_value=-100.0 )
  id_finewoodlitter_fast_dissolved_C = register_tiled_diag_field ( module_name, 'fast_finewoodlitter_dissolved_C', axes(1:2),  &
       land_time, 'fast fine wood litter dissolved carbon', 'kg C/m2', missing_value=-100.0 )
  id_finewoodlitter_slow_dissolved_C = register_tiled_diag_field ( module_name, 'slow_finewoodlitter_dissolved_C', axes(1:2),  &
       land_time, 'slow fine wood litter dissolved carbon', 'kg C/m2', missing_value=-100.0 )
  id_finewoodlitter_deadmic_dissolved_C = register_tiled_diag_field ( module_name, 'finewoodlitter_dead_microbe_dissolved_C', axes(1:2),  &
       land_time, 'dead microbe fine wood litter dissolved carbon', 'kg C/m2', missing_value=-100.0 )
  id_coarsewoodlitter_fast_dissolved_C = register_tiled_diag_field ( module_name, 'fast_coarsewoodlitter_dissolved_C', axes(1:2),  &
       land_time, 'fast coarse woodlitter dissolved carbon', 'kg C/m2', missing_value=-100.0 )
  id_coarsewoodlitter_slow_dissolved_C = register_tiled_diag_field ( module_name, 'slow_coarsewoodlitter_dissolved_C', axes(1:2),  &
       land_time, 'slow coarse wood litter dissolved carbon', 'kg C/m2', missing_value=-100.0 )
  id_coarsewoodlitter_deadmic_dissolved_C = register_tiled_diag_field ( module_name, 'coarsewoodlitter_dead_microbe_dissolved_C', axes(1:2),  &
       land_time, 'dead microbe coarse wood litter dissolved carbon', 'kg C/m2', missing_value=-100.0 )

       id_leaflitter_fast_dissolved_N = register_tiled_diag_field ( module_name, 'fast_leaflitter_dissolved_N', axes(1:2),  &
            land_time, 'fast leaf litter dissolved nitrogen', 'kg N/m2', missing_value=-100.0 )
       id_leaflitter_slow_dissolved_N = register_tiled_diag_field ( module_name, 'slow_leaflitter_dissolved_N', axes(1:2),  &
            land_time, 'slow leaf litter dissolved nitrogen', 'kg N/m2', missing_value=-100.0 )
       id_leaflitter_deadmic_dissolved_N = register_tiled_diag_field ( module_name, 'leaflitter_dead_microbe_dissolved_N', axes(1:2),  &
            land_time, 'dead microbe leaf litter dissolved nitrogen', 'kg N/m2', missing_value=-100.0 )
       id_finewoodlitter_fast_dissolved_N = register_tiled_diag_field ( module_name, 'fast_finewoodlitter_dissolved_N', axes(1:2),  &
            land_time, 'fast fine wood litter dissolved nitrogen', 'kg N/m2', missing_value=-100.0 )
       id_finewoodlitter_slow_dissolved_N = register_tiled_diag_field ( module_name, 'slow_finewoodlitter_dissolved_N', axes(1:2),  &
            land_time, 'slow fine wood litter dissolved nitrogen', 'kg N/m2', missing_value=-100.0 )
       id_finewoodlitter_deadmic_dissolved_N = register_tiled_diag_field ( module_name, 'finewoodlitter_dead_microbe_dissolved_N', axes(1:2),  &
            land_time, 'dead microbe fine wood litter dissolved nitrogen', 'kg N/m2', missing_value=-100.0 )
       id_coarsewoodlitter_fast_dissolved_N = register_tiled_diag_field ( module_name, 'fast_coarsewoodlitter_dissolved_N', axes(1:2),  &
            land_time, 'fast coarse woodlitter dissolved nitrogen', 'kg N/m2', missing_value=-100.0 )
       id_coarsewoodlitter_slow_dissolved_N = register_tiled_diag_field ( module_name, 'slow_coarsewoodlitter_dissolved_N', axes(1:2),  &
            land_time, 'slow coarse wood litter dissolved nitrogen', 'kg N/m2', missing_value=-100.0 )
       id_coarsewoodlitter_deadmic_dissolved_N = register_tiled_diag_field ( module_name, 'coarsewoodlitter_dead_microbe_dissolved_N', axes(1:2),  &
            land_time, 'dead microbe coarse wood litter dissolved nitrogen', 'kg N/m2', missing_value=-100.0 )


  id_livemic_C = register_tiled_diag_field ( module_name, 'live_microbe_C', axes,  &
       land_time, 'Live microbe soil carbon', 'kg C/m3', missing_value=-100.0 )
  id_nsoilcohorts = register_tiled_diag_field ( module_name, 'n_soil_cohorts', axes,  &
       land_time, 'number of soil cohorts', missing_value=-100.0 )
  id_nleaflittercohorts = register_tiled_diag_field ( module_name, 'n_leaflitter_cohorts', axes(1:2),  &
       land_time, 'number of leaf litter cohorts', missing_value=-100.0 )
  id_nfinewoodlittercohorts = register_tiled_diag_field ( module_name, 'n_finewoodlitter_cohorts', axes(1:2),  &
       land_time, 'number of fine wood litter cohorts', missing_value=-100.0 )
  id_ncoarsewoodlittercohorts = register_tiled_diag_field ( module_name, 'n_coarsewoodlitter_cohorts', axes(1:2),  &
       land_time, 'number of coarse wood litter cohorts', missing_value=-100.0 )
  id_deadmic_C_total = register_tiled_diag_field ( module_name, 'deadmic_total_C', axes(1:2),  &
       land_time, 'total dead microbe soil carbon including litter', 'kg C/m2', missing_value=-100.0 )
  id_livemic_C_total = register_tiled_diag_field ( module_name, 'livemic_total_C', axes(1:2),  &
       land_time, 'total live microbe soil carbon including litter', 'kg C/m2', missing_value=-100.0 )
  id_protected_C_total = register_tiled_diag_field ( module_name, 'protected_total_C', axes(1:2),  &
       land_time, 'total protected soil carbon', 'kg C/m2', missing_value=-100.0 )
  id_dissolved_C_total = register_tiled_diag_field ( module_name, 'dissolved_total_C', axes(1:2),  &
       land_time, 'total dissolved soil carbon including litter', 'kg C/m2', missing_value=-100.0 )
  id_total_soil_C = register_tiled_diag_field ( module_name, 'total_soil_C', axes(1:2),  &
       land_time, 'total soil carbon including litter', 'kg C/m2', missing_value=-100.0 )

       id_deadmic_N_total = register_tiled_diag_field ( module_name, 'deadmic_total_N', axes(1:2),  &
            land_time, 'total dead microbe soil nitrogen including litter', 'kg N/m2', missing_value=-100.0 )
       id_livemic_N_total = register_tiled_diag_field ( module_name, 'livemic_total_N', axes(1:2),  &
            land_time, 'total live microbe soil nitrogen including litter', 'kg N/m2', missing_value=-100.0 )
       id_protected_N_total = register_tiled_diag_field ( module_name, 'protected_total_N', axes(1:2),  &
            land_time, 'total protected soil nitrogen including litter', 'kg N/m2', missing_value=-100.0 )
       id_dissolved_N_total = register_tiled_diag_field ( module_name, 'dissolved_total_N', axes(1:2),  &
            land_time, 'total dissolved soil nitrogen including litter', 'kg N/m2', missing_value=-100.0 )
       id_total_soil_N = register_tiled_diag_field ( module_name, 'total_soil_N', axes(1:2),  &
            land_time, 'total soil nitrogen including litter', 'kg N/m2', missing_value=-100.0 )
        id_total_soil_organic_N = register_tiled_diag_field ( module_name, 'total_soil_organic_N', axes(1:2),  &
             land_time, 'total soil organic nitrogen', 'kg N/m2', missing_value=-100.0 )
         id_total_soil_ammonium = register_tiled_diag_field ( module_name, 'total_soil_ammonium', axes(1:2),  &
              land_time, 'total soil ammonium including litter', 'kg N/m2', missing_value=-100.0 )
          id_total_soil_nitrate = register_tiled_diag_field ( module_name, 'total_soil_nitrate', axes(1:2),  &
               land_time, 'total soil nitrate including litter', 'kg N/m2', missing_value=-100.0 )

  id_fast_C_leaching = register_tiled_diag_field ( module_name, 'fast_C_leaching', axes, &
       land_time, 'net layer fast soil C leaching',  'kg/(m2 s)', missing_value=-100.0)
  id_slow_C_leaching = register_tiled_diag_field ( module_name, 'slow_C_leaching', axes, &
       land_time, 'net layer slow soil C leaching',  'kg/(m2 s)', missing_value=-100.0)
  id_deadmic_C_leaching = register_tiled_diag_field ( module_name, 'deadmic_C_leaching', axes, &
       land_time, 'net layer dead microbe soil C leaching',  'kg/(m2 s)', missing_value=-100.0)
  id_total_C_leaching = register_tiled_diag_field ( module_name, 'total_C_leaching', axes, &
       land_time, 'total vertical soil C leaching', 'kg/(m2 s)', missing_value=initval)

       id_fast_N_leaching = register_tiled_diag_field ( module_name, 'fast_N_leaching', axes, &
            land_time, 'net layer fast soil N leaching',  'kg/(m2 s)', missing_value=-100.0)
       id_slow_N_leaching = register_tiled_diag_field ( module_name, 'slow_N_leaching', axes, &
            land_time, 'net layer slow soil N leaching',  'kg/(m2 s)', missing_value=-100.0)
       id_deadmic_N_leaching = register_tiled_diag_field ( module_name, 'deadmic_N_leaching', axes, &
            land_time, 'net layer dead microbe soil N leaching',  'kg/(m2 s)', missing_value=-100.0)
       id_total_ON_leaching = register_tiled_diag_field ( module_name, 'total_ON_leaching', axes, &
            land_time, 'total vertical soil organic N leaching', 'kg/(m2 s)', missing_value=initval)
    id_NO3_leaching = register_tiled_diag_field ( module_name, 'NO3_leaching', axes, &
         land_time, 'net layer NO3 leaching',  'kgN/(m2 s)', missing_value=-100.0)
     id_NH4_leaching = register_tiled_diag_field ( module_name, 'NH4_leaching', axes, &
          land_time, 'net layer NH4 leaching',  'kgN/(m2 s)', missing_value=-100.0)
      id_total_NO3_leaching = register_tiled_diag_field ( module_name, 'total_NO3_leaching', axes, &
           land_time, 'Total NO3 leaching',  'kgN/(m2 s)', missing_value=-100.0)
       id_total_NH4_leaching = register_tiled_diag_field ( module_name, 'total_NH4_leaching', axes, &
            land_time, 'Total NH4 leaching',  'kgN/(m2 s)', missing_value=-100.0)

  id_livemic_C_leaching = register_tiled_diag_field ( module_name, 'livemic_C_leaching', axes, &
       land_time, 'net layer live microbe C leaching',  'kg/(m2 s)', missing_value=-100.0)
   id_livemic_N_leaching = register_tiled_diag_field ( module_name, 'livemic_N_leaching', axes, &
        land_time, 'net layer live microbe N leaching',  'kg/(m2 s)', missing_value=-100.0)
  id_protected_C_leaching = register_tiled_diag_field ( module_name, 'protected_C_leaching', axes, &
      land_time, 'net layer protected soil C leaching',  'kg/(m2 s)', missing_value=-100.0)
  id_protected_N_leaching = register_tiled_diag_field ( module_name, 'protected_N_leaching', axes, &
      land_time, 'net layer protected soil N leaching',  'kg/(m2 s)', missing_value=-100.0)
  id_leaflitter_fast_C_leaching = register_tiled_diag_field ( module_name, 'fast_leaflitter_C_leaching', axes(1:2), &
        land_time, 'Leaf litter fast C leaching','kg/(m2 s)', missing_value=-100.0)
  id_leaflitter_slow_C_leaching = register_tiled_diag_field ( module_name, 'slow_leaflitter_C_leaching', axes(1:2), &
        land_time, 'Leaf litter slow C leaching','kg/(m2 s)', missing_value=-100.0)
  id_leaflitter_deadmic_C_leaching = register_tiled_diag_field ( module_name, 'deadmic_leaflitter_C_leaching', axes(1:2), &
        land_time, 'Leaf litter dead microbe C leaching','kg/(m2 s)', missing_value=-100.0)
  id_coarsewoodlitter_fast_C_leaching = register_tiled_diag_field ( module_name, 'fast_coarsewoodlitter_C_leaching', axes(1:2), &
        land_time, 'Coarse wood litter fast C leaching','kg/(m2 s)', missing_value=-100.0)
  id_coarsewoodlitter_slow_C_leaching = register_tiled_diag_field ( module_name, 'slow_coarsewoodlitter_C_leaching', axes(1:2), &
        land_time, 'Coarse wood litter slow C leaching','kg/(m2 s)', missing_value=-100.0)
  id_coarsewoodlitter_deadmic_C_leaching = register_tiled_diag_field ( module_name, 'deadmic_coarsewoodlitter_C_leaching', axes(1:2), &
        land_time, 'Coarse wood litter dead microbe C leaching','kg/(m2 s)', missing_value=-100.0)
  id_finewoodlitter_fast_C_leaching = register_tiled_diag_field ( module_name, 'fast_finewoodlitter_C_leaching', axes(1:2), &
        land_time, 'Fine wood litter fast C leaching','kg/(m2 s)', missing_value=-100.0)
  id_finewoodlitter_slow_C_leaching = register_tiled_diag_field ( module_name, 'slow_finewoodlitter_C_leaching', axes(1:2), &
        land_time, 'Fine wood litter slow C leaching','kg/(m2 s)', missing_value=-100.0)
  id_finewoodlitter_deadmic_C_leaching = register_tiled_diag_field ( module_name, 'deadmic_finewoodlitter_C_leaching', axes(1:2), &
        land_time, 'Fine wood litter dead microbe C leaching','kg/(m2 s)', missing_value=-100.0)

        id_leaflitter_fast_N_leaching = register_tiled_diag_field ( module_name, 'fast_leaflitter_N_leaching', axes(1:2), &
              land_time, 'Leaf litter fast N leaching','kg/(m2 s)', missing_value=-100.0)
        id_leaflitter_slow_N_leaching = register_tiled_diag_field ( module_name, 'slow_leaflitter_N_leaching', axes(1:2), &
              land_time, 'Leaf litter slow N leaching','kg/(m2 s)', missing_value=-100.0)
        id_leaflitter_deadmic_N_leaching = register_tiled_diag_field ( module_name, 'deadmic_leaflitter_N_leaching', axes(1:2), &
              land_time, 'Leaf litter dead microbe N leaching','kg/(m2 s)', missing_value=-100.0)
        id_coarsewoodlitter_fast_N_leaching = register_tiled_diag_field ( module_name, 'fast_coarsewoodlitter_N_leaching', axes(1:2), &
              land_time, 'Coarse wood litter fast N leaching','kg/(m2 s)', missing_value=-100.0)
        id_coarsewoodlitter_slow_N_leaching = register_tiled_diag_field ( module_name, 'slow_coarsewoodlitter_N_leaching', axes(1:2), &
              land_time, 'Coarse wood litter slow N leaching','kg/(m2 s)', missing_value=-100.0)
        id_coarsewoodlitter_deadmic_C_leaching = register_tiled_diag_field ( module_name, 'deadmic_coarsewoodlitter_C_leaching', axes(1:2), &
              land_time, 'Coarse wood litter dead microbe N leaching','kg/(m2 s)', missing_value=-100.0)
        id_finewoodlitter_fast_N_leaching = register_tiled_diag_field ( module_name, 'fast_finewoodlitter_N_leaching', axes(1:2), &
              land_time, 'Fine wood litter fast N leaching','kg/(m2 s)', missing_value=-100.0)
        id_finewoodlitter_slow_N_leaching = register_tiled_diag_field ( module_name, 'slow_finewoodlitter_N_leaching', axes(1:2), &
              land_time, 'Fine wood litter slow N leaching','kg/(m2 s)', missing_value=-100.0)
        id_finewoodlitter_deadmic_N_leaching = register_tiled_diag_field ( module_name, 'deadmic_finewoodlitter_N_leaching', axes(1:2), &
              land_time, 'Fine wood litter dead microbe N leaching','kg/(m2 s)', missing_value=-100.0)
      id_finewoodlitter_NO3_leaching = register_tiled_diag_field ( module_name, 'finewoodlitter_NO3_leaching', axes(1:2), &
            land_time, 'Fine wood litter NO3 leaching','kg/(m2 s)', missing_value=-100.0)
    id_coarsewoodlitter_NO3_leaching = register_tiled_diag_field ( module_name, 'coarsewoodlitter_NO3_leaching', axes(1:2), &
          land_time, 'Coarse wood litter NO3 leaching','kg/(m2 s)', missing_value=-100.0)
      id_leaflitter_NO3_leaching = register_tiled_diag_field ( module_name, 'leaflitter_NO3_leaching', axes(1:2), &
            land_time, 'Leaf litter NO3 leaching','kg/(m2 s)', missing_value=-100.0)
            id_finewoodlitter_NH4_leaching = register_tiled_diag_field ( module_name, 'finewoodlitter_NH4_leaching', axes(1:2), &
                  land_time, 'Fine wood litter NH4 leaching','kg/(m2 s)', missing_value=-100.0)
          id_coarsewoodlitter_NH4_leaching = register_tiled_diag_field ( module_name, 'coarsewoodlitter_NH4_leaching', axes(1:2), &
                land_time, 'Coarse wood litter NH4 leaching','kg/(m2 s)', missing_value=-100.0)
            id_leaflitter_NH4_leaching = register_tiled_diag_field ( module_name, 'leaflitter_NH4_leaching', axes(1:2), &
                  land_time, 'Leaf litter NH4 leaching','kg/(m2 s)', missing_value=-100.0)


        id_root_profile = register_tiled_diag_field ( module_name, 'root_profile', axes, &
              land_time, 'Vertical root profile in soil','Fraction', missing_value=-100.0)

  ! ZMS
  id_total_soil_carbon = register_tiled_diag_field ( module_name, 'total_lit_SOM_C', axes(1:2), &
       land_time, 'vertical sum of all litter and soil carbon pools', 'kg C/m^2', missing_value=-100.0)
  id_surf_DOC_loss = register_tiled_diag_field ( module_name, 'surf_DOC_loss', axes(1:2), &
       land_time, 'loss of top layer DOC to surface runoff due to efflux', 'kg C/m^2/s', &
       missing_value=initval)

       id_surf_DON_loss = register_tiled_diag_field ( module_name, 'surf_DON_loss', axes(1:2), &
            land_time, 'loss of top layer DON to surface runoff due to efflux', 'kg N/m^2/s', &
            missing_value=initval)
        id_surf_NO3_loss = register_tiled_diag_field ( module_name, 'surf_NO3_loss', axes(1:2), &
             land_time, 'loss of top layer NO3 to surface runoff due to efflux', 'kg N/m^2/s', &
             missing_value=initval)
         id_surf_NH4_loss = register_tiled_diag_field ( module_name, 'surf_NH4_loss', axes(1:2), &
              land_time, 'loss of top layer NH4 to surface runoff due to efflux', 'kg N/m^2/s', &
              missing_value=initval)

  id_fsc = register_tiled_diag_field ( module_name, 'fsc', axes(1:2),  &
       land_time, 'total fast soil carbon', 'kg C/m2', missing_value=-100.0 )
  id_ssc = register_tiled_diag_field ( module_name, 'ssc', axes(1:2),  &
       land_time, 'total slow soil carbon', 'kg C/m2', missing_value=-100.0 )
   id_fsn = register_tiled_diag_field ( module_name, 'fsn', axes(1:2),  &
        land_time, 'total fast soil nitrogen', 'kg N/m2', missing_value=-100.0 )
   id_ssn = register_tiled_diag_field ( module_name, 'ssn', axes(1:2),  &
        land_time, 'total slow soil nitrogen', 'kg N/m2', missing_value=-100.0 )
  id_lwc = register_tiled_diag_field ( module_name, 'soil_liq', axes,  &
       land_time, 'bulk density of liquid water', 'kg/m3', missing_value=-100.0 )
  id_swc  = register_tiled_diag_field ( module_name, 'soil_ice',  axes,  &
       land_time, 'bulk density of solid water', 'kg/m3',  missing_value=-100.0 )
  id_psi = register_tiled_diag_field ( module_name, 'soil_psi', axes,  &
       land_time, 'soil-water matric head', 'm', missing_value=-100.0 )
  id_temp  = register_tiled_diag_field ( module_name, 'soil_T',  axes,       &
       land_time, 'temperature',            'degK',  missing_value=-100.0 )
  id_ie  = register_tiled_diag_field ( module_name, 'soil_rie',  axes(1:2),  &
       land_time, 'inf exc runf',            'kg/(m2 s)',  missing_value=-100.0 )
  id_sn  = register_tiled_diag_field ( module_name, 'soil_rsn',  axes(1:2),  &
       land_time, 'satn runf',            'kg/(m2 s)',  missing_value=-100.0 )
  id_bf  = register_tiled_diag_field ( module_name, 'soil_rbf',  axes(1:2),  &
       land_time, 'baseflow',            'kg/(m2 s)',  missing_value=-100.0 )
  id_if  = register_tiled_diag_field ( module_name, 'soil_rif',  axes(1:2),  &
       land_time, 'interflow',            'kg/(m2 s)',  missing_value=-100.0 )
  id_al  = register_tiled_diag_field ( module_name, 'soil_ral',  axes(1:2),  &
       land_time, 'active layer flow',    'kg/(m2 s)',  missing_value=-100.0 )
  id_nu  = register_tiled_diag_field ( module_name, 'soil_rnu',  axes(1:2),  &
       land_time, 'numerical runoff',    'kg/(m2 s)',  missing_value=-100.0 )
  id_sc  = register_tiled_diag_field ( module_name, 'soil_rsc',  axes(1:2),  &
       land_time, 'lm2 groundwater runoff',    'kg/(m2 s)',  missing_value=-100.0 )
  id_hie  = register_tiled_diag_field ( module_name, 'soil_hie',  axes(1:2), &
       land_time, 'heat ie runf',            'W/m2',  missing_value=-100.0 )
  id_hsn  = register_tiled_diag_field ( module_name, 'soil_hsn',  axes(1:2), &
       land_time, 'heat sn runf',            'W/m2',  missing_value=-100.0 )
  id_hbf  = register_tiled_diag_field ( module_name, 'soil_hbf',  axes(1:2), &
       land_time, 'heat bf runf',            'W/m2',  missing_value=-100.0 )
  id_hif  = register_tiled_diag_field ( module_name, 'soil_hif',  axes(1:2), &
       land_time, 'heat if runf',            'W/m2',  missing_value=-100.0 )
  id_hal  = register_tiled_diag_field ( module_name, 'soil_hal',  axes(1:2), &
       land_time, 'heat al runf',            'W/m2',  missing_value=-100.0 )
  id_hnu  = register_tiled_diag_field ( module_name, 'soil_hnu',  axes(1:2), &
       land_time, 'heat nu runoff',          'W/m2',  missing_value=-100.0 )
  id_hsc  = register_tiled_diag_field ( module_name, 'soil_hsc',  axes(1:2), &
       land_time, 'heat sc runoff',          'W/m2',  missing_value=-100.0 )
  id_evap  = register_tiled_diag_field ( module_name, 'soil_evap',  axes(1:2), &
       land_time, 'soil evap',            'kg/(m2 s)',  missing_value=-100.0 )
  id_excess  = register_tiled_diag_field ( module_name, 'sfc_excess',  axes(1:2),  &
       land_time, 'sfc excess pushed down',    'kg/(m2 s)',  missing_value=-100.0 )

  id_uptk_n_iter  = register_tiled_diag_field ( module_name, 'uptake_n_iter',  axes(1:2), &
       land_time, 'number of iterations for soil uptake',  missing_value=-100.0 )
  id_uptk = register_tiled_diag_field ( module_name, 'soil_uptk', axes, &
       land_time, 'uptake of water by roots', 'kg/(m2 s)',  missing_value=-100.0 )
  id_psi_x0 = register_tiled_diag_field ( module_name, 'soil_psix0', axes(1:2), &
       land_time, 'xylem potential at z=0', 'm',  missing_value=-100.0 )
  id_sws_n_iter  = register_tiled_diag_field ( module_name, 'sws_n_iter',  axes(1:2), &
       land_time, 'number of iterations for soil water supply',  missing_value=-100.0 )
  id_psi_x0_sws = register_tiled_diag_field ( module_name, 'soil_psix0_sws', axes(1:2), &
       land_time, 'xylem potential at z=0 for max transpiration', 'm',  missing_value=-100.0 )
  id_deficit = register_tiled_diag_field ( module_name, 'soil_def', axes(1:2), &
       land_time, 'groundwater storage deficit', '-',  missing_value=-100.0 )
  id_deficit_2 = register_tiled_diag_field ( module_name, 'soil_def2', axes(1:2), &
       land_time, 'groundwater storage deficit2', '-',  missing_value=-100.0 )
  id_deficit_3 = register_tiled_diag_field ( module_name, 'soil_def3', axes(1:2), &
       land_time, 'groundwater storage deficit3', '-',  missing_value=-100.0 )
  id_deficit_4 = register_tiled_diag_field ( module_name, 'soil_def4', axes(1:2), &
       land_time, 'groundwater storage deficit4', '-',  missing_value=-100.0 )
  id_psi_bot = register_tiled_diag_field ( module_name, 'soil_psi_n', axes(1:2), &
       land_time, 'psi at bottom of soil column', 'm',  missing_value=-100.0 )
  id_sat_frac = register_tiled_diag_field ( module_name, 'soil_fsat', axes(1:2), &
       land_time, 'fraction of soil area saturated at surface', '-',  missing_value=-100.0 )
  id_stor_frac = register_tiled_diag_field ( module_name, 'soil_fgw', axes(1:2), &
       land_time, 'groundwater storage frac above base elev', '-',  missing_value=-100.0 )
  id_sat_depth = register_tiled_diag_field ( module_name, 'soil_wtdep', axes(1:2), &
       land_time, 'depth below sfc to saturated soil', 'm',  missing_value=-100.0 )
  id_sat_dept2 = register_tiled_diag_field ( module_name, 'soil_wtdp2', axes(1:2), &
       land_time, 'alt depth below sfc to saturated soil', 'm',  missing_value=-100.0 )
  id_z_cap = register_tiled_diag_field ( module_name, 'soil_zcap', axes(1:2), &
       land_time, 'depth below sfc to capillary fringe', 'm',  missing_value=-100.0 )

  id_div_bf = register_tiled_diag_field ( module_name, 'soil_dvbf', axes, &
       land_time, 'baseflow by layer', 'kg/(m2 s)',  missing_value=-100.0 )
  id_div_if = register_tiled_diag_field ( module_name, 'soil_dvif', axes, &
       land_time, 'interflow by layer', 'kg/(m2 s)',  missing_value=-100.0 )
  id_div_al = register_tiled_diag_field ( module_name, 'soil_dval', axes, &
       land_time, 'active-layer flow by layer', 'kg/(m2 s)',  missing_value=-100.0 )

  id_cf_1 = register_tiled_diag_field ( module_name, 'soil_cf_1', axes(1:2), &
       land_time, 'soil_cf_1', 'm',  missing_value=-100.0 )
  id_cf_3 = register_tiled_diag_field ( module_name, 'soil_cf_3', axes(1:2), &
       land_time, 'soil_cf_3', 'm',  missing_value=-100.0 )
  id_wt_1 = register_tiled_diag_field ( module_name, 'soil_wt_1', axes(1:2), &
       land_time, 'soil_wt_1', 'm',  missing_value=-100.0 )
  id_wt_2 = register_tiled_diag_field ( module_name, 'soil_wt_2', axes(1:2), &
       land_time, 'soil_wt_2', 'm',  missing_value=-100.0 )
  id_wt_2a = register_tiled_diag_field ( module_name, 'soil_wt_2a', axes(1:2), &
       land_time, 'Water Table Depth from Surface to Saturation', 'm',  missing_value=-100.0 )
  id_wt_2b = register_tiled_diag_field ( module_name, 'soil_wt_2b', axes(1:2), &
       land_time, 'Water Table Depth from Surface to Liquid Saturation', 'm',  missing_value=-100.0 )
  id_wt_3 = register_tiled_diag_field ( module_name, 'soil_wt_3', axes(1:2), &
       land_time, 'soil_wt_3', 'm',  missing_value=-100.0 )
  id_wt2_3 = register_tiled_diag_field ( module_name, 'soil_wt2_3', axes(1:2), &
       land_time, 'soil_wt2_3', 'm',  missing_value=-100.0 )
  id_wt_4 = register_tiled_diag_field ( module_name, 'soil_wt_4', axes(1:2), &
       land_time, 'Interpolated psi = 0 from Bottom Up', 'm',  missing_value=-100.0 )

  id_active_layer = register_tiled_diag_field ( module_name, 'soil_alt', axes(1:2), &
       land_time, 'active-layer thickness', 'm',  missing_value=-100.0 )
  id_heat_cap = register_tiled_diag_field ( module_name, 'soil_heat_cap',  &
       axes, land_time, 'heat capacity of dry soil','J/(m3 K)', missing_value=-100.0 )
  id_thermal_cond =  register_tiled_diag_field ( module_name, 'soil_tcon', &
       axes, land_time, 'soil thermal conductivity', 'W/(m K)',  missing_value=-100.0 )

  id_surface_water = register_tiled_diag_field (module_name, 'surface_water', &
       axes(1:2), land_time, 'surface water storage', 'm', missing_value=-100.0 )
  id_inun_frac = register_tiled_diag_field (module_name, 'inun_fraction', &
       axes(1:2), land_time, 'inundated area fraction', '-', missing_value=-100.0 )
  if (gw_option == GW_TILED) then
     id_wet_frac = register_tiled_diag_field (module_name, 'wet_fraction', &
       axes(1:2), land_time, 'diagnostic wetland fraction', '-', missing_value=-100.0 )
  end if
  if (gw_option == GW_TILED .and. simple_inundation) then
      id_rsn_frac = register_tiled_diag_field (module_name, 'surface_runoff_frac', &
         axes(1:2), land_time, 'effective fraction of throughfall converted to sat-excess surface runoff', '-', missing_value=-100.0 )
  end if
  id_flow = register_tiled_diag_field (module_name, 'flow', axes, &
       land_time, 'vertical soil water flow at interface above (+ downward)', 'mm/s', missing_value=initval )
  id_reflux = register_tiled_diag_field (module_name, 'reflux', axes(1:2), &
       land_time, 'upwards flow of soil water at surface; zero if flow into surface', 'mm/s', missing_value=-100.0 )
  id_macro_infilt = register_tiled_diag_field (module_name, 'macro_inf', axes(1:2), &
       land_time, 'infiltration (decrease to IE runoff) at soil surface due to vertical macroporosity', 'mm/s', missing_value=-100.0 )

  id_type = register_tiled_static_field ( module_name, 'soil_type',  &
       axes(1:2), 'soil type', missing_value=-1.0 )
  id_tau_gw = register_tiled_static_field ( module_name, 'tau_gw',  &
       axes(1:2), 'groundwater residence time', 's', missing_value=-100.0 )
  id_slope_l = register_tiled_static_field ( module_name, 'slope_l',  &
       axes(1:2), 'hillslope length', 'm', missing_value=-100.0 )
  id_slope_Z = register_tiled_static_field ( module_name, 'soil_rlief',  &
       axes(1:2), 'hillslope relief', 'm', missing_value=-100.0 )
  id_zeta_bar = register_tiled_static_field ( module_name, 'zeta_bar',  &
       axes(1:2), 'hillslope zeta bar', '-', missing_value=-100.0 )
  id_e_depth = register_tiled_static_field ( module_name, 'soil_depth',  &
       axes(1:2), 'soil hydraulic e-folding depth', 'm', missing_value=-100.0 )
  id_zeta = register_tiled_static_field ( module_name, 'soil_zeta',      &
       axes(1:2), 'soil depth/topo relief', '-',  missing_value=-100.0 )
  id_tau = register_tiled_static_field ( module_name, 'soil_tau',        &
       axes(1:2), 'gw transmissivity/soil transmissivity', '-',  missing_value=-100.0 )
  id_vwc_wilt = register_tiled_static_field ( module_name, 'soil_wilt',  &
       axes(1:2), 'wilting water content', '-', missing_value=-100.0 )
  id_vwc_fc = register_tiled_static_field ( module_name, 'soil_fc',  &
       axes(1:2), 'field capacity', '-', missing_value=-100.0 )
  id_vwc_sat = register_tiled_static_field ( module_name, 'soil_sat',  &
       axes(1:2), 'soil porosity', '-', missing_value=-100.0 )
  id_K_sat = register_tiled_static_field ( module_name, 'soil_Ksat',  &
       axes(1:2), 'soil sat. hydraulic conductivity', 'kg /(m2 s)', missing_value=-100.0 )
  id_K_gw  = register_tiled_static_field ( module_name, 'soil_K_gw',  &
       axes(1:2), 'deep hydraulic conductivity', 'kg /(m2 s)', missing_value=-100.0 )
  id_w_fc = register_tiled_static_field ( module_name, 'w_fc',  &
       axes, 'soil field capacity', missing_value=-1.0 )
  id_alpha = register_tiled_static_field ( module_name, 'soil_alpha',  &
       axes, 'soil microscopic length scale', missing_value=-1.0 )
  id_refl_dry_dir = register_tiled_static_field ( module_name, 'refl_dry_dir',  &
       (/id_lon, id_lat, id_band/), 'reflectance of dry soil for direct light', &
       missing_value=-1.0 )
  id_refl_dry_dif = register_tiled_static_field ( module_name, 'refl_dry_dif',  &
       (/id_lon, id_lat, id_band/), 'reflectance of dry soil for diffuse light', &
       missing_value=-1.0 )
  id_refl_sat_dir = register_tiled_static_field ( module_name, 'refl_sat_dir',  &
       (/id_lon, id_lat, id_band/), 'reflectance of saturated soil for direct light', &
       missing_value=-1.0 )
  id_refl_sat_dif = register_tiled_static_field ( module_name, 'refl_sat_dif',  &
       (/id_lon, id_lat, id_band/), 'reflectance of saturated soil for diffuse light', &
       missing_value=-1.0 )
  id_f_iso_dry = register_tiled_static_field ( module_name, 'f_iso_dry',  &
       (/id_lon, id_lat, id_band/), 'isotropic brdf weight, dry soil', &
       missing_value=-1.0 )
  id_f_vol_dry = register_tiled_static_field ( module_name, 'f_vol_dry',  &
       (/id_lon, id_lat, id_band/), 'volumetric brdf weight, dry soil', &
       missing_value=-1.0 )
  id_f_geo_dry = register_tiled_static_field ( module_name, 'f_geo_dry',  &
       (/id_lon, id_lat, id_band/), 'geometric brdf weight, dry soil', &
       missing_value=-1.0 )
  id_f_iso_sat = register_tiled_static_field ( module_name, 'f_iso_sat',  &
       (/id_lon, id_lat, id_band/), 'isotropic brdf weight, saturated soil', &
       missing_value=-1.0 )
  id_f_vol_sat = register_tiled_static_field ( module_name, 'f_vol_sat',  &
       (/id_lon, id_lat, id_band/), 'volumetric brdf weight, saturated soil', &
       missing_value=-1.0 )
  id_f_geo_sat = register_tiled_static_field ( module_name, 'f_geo_sat',  &
       (/id_lon, id_lat, id_band/), 'geometric brdf weight, saturated soil', &
       missing_value=-1.0 )

  id_asoil = register_tiled_diag_field ( module_name, 'asoil', &
       (/id_lon,id_lat/), land_time, 'aerobic activity modifier', &
       missing_value=-100.0 )

  ! the following fields are for compatibility with older diag tables only
  call add_tiled_static_field_alias ( id_slope_Z, module_name, 'slope_Z',  &
       axes(1:2), 'hillslope relief (obsolete, use "soil_rlief" instead)',&
       'm', missing_value=-100.0 )
  call add_tiled_static_field_alias ( id_e_depth, module_name, 'e_depth',  &
       axes(1:2), 'soil e-folding depth (obsolete, use "soil_depth" instead)', &
       'm', missing_value=-100.0 )
  call add_tiled_static_field_alias ( id_vwc_wilt, module_name, 'vwc_wilt',  &
       axes(1:2), 'wilting water content (obsolete, use "soil_wilt" instead)', &
       '-', missing_value=-100.0 )
  call add_tiled_static_field_alias ( id_vwc_fc, module_name, 'vwc_fc',  &
       axes(1:2), 'field capacity (obsolete, use "soil_fc" instead)', &
       '-', missing_value=-100.0 )
  call add_tiled_static_field_alias ( id_vwc_sat, module_name, 'vwc_sat',  &
       axes(1:2), 'soil porosity (obsolete, use "soil_sat")', &
       '-', missing_value=-100.0 )
  call add_tiled_static_field_alias ( id_K_sat, module_name, 'K_sat',  &
       axes(1:2), 'soil sat. hydraulic conductivity (obsolte, use "soil_Ksat" instead)', &
       'kg /(m2 s)', missing_value=-100.0 )

#ifdef ZMSDEBUG_TRIDIAGTEST
  ! For testing tridiagonal solution for advection
  id_st_diff = register_tiled_diag_field ( module_name, 'soil_T_diff', axes,  &
      land_time, 'soil Temperature difference after advection with tridiagonal solution', &
      'K', missing_value=-100.0 )
#endif

end subroutine soil_diag_init


! ============================================================================
subroutine soil_end ()

  module_is_initialized =.FALSE.

end subroutine soil_end


! ============================================================================
subroutine save_soil_restart (tile_dim_length, timestamp)
  integer, intent(in) :: tile_dim_length ! length of tile dim. in the output file
  character(*), intent(in) :: timestamp ! timestamp to add to the file name

  ! ---- local vars ----------------------------------------------------------
  integer :: unit            ! restart file i/o unit

  call error_mesg('soil_end','writing NetCDF restart',NOTE)
  ! create output file, including internal structure necessary for output
  call create_tile_out_file(unit,'RESTART/'//trim(timestamp)//'soil.res.nc', &
          lnd%coord_glon, lnd%coord_glat, soil_tile_exists, tile_dim_length )
  ! in addition, define vertical coordinate
  if (mpp_pe()==lnd%io_pelist(1)) then
     __NF_ASRT__(nfu_def_dim(unit,'zfull',zfull(1:num_l),'full level','m'))
     __NF_ASRT__(nfu_put_att(unit,'zfull','positive','down'))
     if (soil_carbon_option==SOILC_CORPSE .or. soil_carbon_option==SOILC_CORPSE_N) then
        __NF_ASRT__(nfu_def_dim(unit,name='soilCCohort',size=soilMaxCohorts,xtype=NF_INT,long_name='Soil carbon cohort'))
     endif
  endif
  call sync_nc_files(unit)

  ! write out fields
  call write_tile_data_r1d_fptr(unit,'temp'         ,soil_T_ptr   ,'zfull','soil temperature','degrees_K')
  call write_tile_data_r1d_fptr(unit,'wl'           ,soil_wl_ptr  ,'zfull','liquid water content','kg/m2')
  call write_tile_data_r1d_fptr(unit,'ws'           ,soil_ws_ptr  ,'zfull','solid water content','kg/m2')
  call write_tile_data_r1d_fptr(unit,'groundwater'  ,soil_groundwater_ptr  ,'zfull')
  call write_tile_data_r1d_fptr(unit,'groundwater_T',soil_groundwater_T_ptr ,'zfull')
  call write_tile_data_r0d_fptr(unit,'uptake_T',     soil_uptake_T_ptr, 'temperature of transpiring water', 'degrees_K')
  select case(soil_carbon_option)
  case (SOILC_CENTURY, SOILC_CENTURY_BY_LAYER)
     call write_tile_data_r1d_fptr(unit,'fsc',          soil_fast_soil_C_ptr,'zfull','fast soil carbon', 'kg C/m2')
     call write_tile_data_r1d_fptr(unit,'ssc',          soil_slow_soil_C_ptr,'zfull','slow soil carbon', 'kg C/m2')
 case (SOILC_CORPSE, SOILC_CORPSE_N)
     call write_tile_data_layered_cohort_fptr(unit,'fast_soil_C',soilc_fast_soil_C_ptr ,'zfull','soilCCohort','Fast soil carbon','kg/m2')
     call write_tile_data_layered_cohort_fptr(unit,'slow_soil_C',soilc_slow_soil_C_ptr ,'zfull','soilCCohort','Slow soil carbon','kg/m2')
     call write_tile_data_layered_cohort_fptr(unit,'deadMic',soilc_deadMicrobeC_ptr ,'zfull','soilCCohort','Dead microbe carbon','kg/m2')
     call write_tile_data_layered_cohort_fptr(unit,'fastProtectedC',soilc_fast_protected_C_ptr ,'zfull','soilCCohort','Protected fast carbon','kg/m2')
     call write_tile_data_layered_cohort_fptr(unit,'slowProtectedC',soilc_slow_protected_C_ptr ,'zfull','soilCCohort','Protected slow carbon','kg/m2')
     call write_tile_data_layered_cohort_fptr(unit,'deadMicrobeProtectedC',soilc_deadMicrobe_protected_C_ptr ,'zfull','soilCCohort','Protected dead microbe carbon','kg/m2')
     call write_tile_data_layered_cohort_fptr(unit,'liveMic',soilc_livingMicrobeC_ptr ,'zfull','soilCCohort','Living microbial carbon','kg/m2')
     call write_tile_data_layered_cohort_fptr(unit,'CO2',soilc_CO2_ptr ,'zfull','soilCCohort','Cohort CO2 generated','kg/m2')
     call write_tile_data_layered_cohort_fptr(unit,'Rtot',soilc_Rtot_ptr ,'zfull','soilCCohort','Total degradation','kg/m2')
     call write_tile_data_layered_cohort_fptr(unit,'originalCohortC',soilc_originalLitterC_ptr ,'zfull','soilCCohort','Cohort original carbon','g/m2')

! Do we want to write N stuff if not using SOILC_CORPSE_N? It should be valid, but all zeros
if (soil_carbon_option == SOILC_CORPSE_N) then
     call write_tile_data_layered_cohort_fptr(unit,'fast_soil_N',soilc_fast_soil_N_ptr ,'zfull','soilCCohort','Fast soil nitrogen','kg/m2')
     call write_tile_data_layered_cohort_fptr(unit,'slow_soil_N',soilc_slow_soil_N_ptr ,'zfull','soilCCohort','Slow soil nitrogen','kg/m2')
     call write_tile_data_layered_cohort_fptr(unit,'deadMic_N',soilc_deadMicrobeN_ptr ,'zfull','soilCCohort','Dead microbe nitrogen','kg/m2')
     call write_tile_data_layered_cohort_fptr(unit,'fastProtectedN',soilc_fast_protected_N_ptr ,'zfull','soilCCohort','Protected fast nitrogen','kg/m2')
     call write_tile_data_layered_cohort_fptr(unit,'slowProtectedN',soilc_slow_protected_N_ptr ,'zfull','soilCCohort','Protected slow nitrogen','kg/m2')
     call write_tile_data_layered_cohort_fptr(unit,'deadMicrobeProtectedN',soilc_deadMicrobe_protected_N_ptr ,'zfull','soilCCohort','Protected dead microbe nitrogen','kg/m2')
     call write_tile_data_layered_cohort_fptr(unit,'liveMicN',soilc_livingMicrobeN_ptr ,'zfull','soilCCohort','Living microbial nitrogen','kg/m2')
     call write_tile_data_layered_cohort_fptr(unit,'originalCohortN',soilc_originalLitterN_ptr ,'zfull','soilCCohort','Cohort original nitrogen','g/m2')
     ! Leaving out cohort-level immobilization and mineralization fields for now -- BNS

     call write_tile_data_r1d_fptr(unit,'soil_DON_fast',soil_fast_DON_ptr,'zfull','Dissolved fast nitrogen','kg/m2')
     call write_tile_data_r1d_fptr(unit,'soil_DON_slow',soil_slow_DON_ptr,'zfull','Dissolved slow nitrogen','kg/m2')
     call write_tile_data_r1d_fptr(unit,'soil_DON_deadmic',soil_deadmicrobe_DON_ptr,'zfull','Dissolved dead microbe nitrogen','kg/m2')
     call write_tile_data_r1d_fptr(unit,'soil_NO3',soil_nitrate_ptr,'zfull','Soil nitrate content','kg/m2')
     call write_tile_data_r1d_fptr(unit,'soil_NH4',soil_ammonium_ptr,'zfull','Soil ammonium content','kg/m2')
     call write_tile_data_r1d_fptr(unit,'soil_nitrif',soil_nitrif_ptr,'zfull','Soil cumulative nitrification','kg/m2')
     call write_tile_data_r1d_fptr(unit,'soil_denitrif',soil_denitrif_ptr,'zfull','Soil cumulative denitrification','kg/m2')

     call write_tile_data_r0d_fptr(unit,'fast_DON_leached',     soil_fast_DON_leached_ptr, 'Cumulative fast DON leached out of the column', 'kg/m2')
     call write_tile_data_r0d_fptr(unit,'slow_DON_leached',     soil_slow_DON_leached_ptr, 'Cumulative slow DON leached out of the column', 'kg/m2')
     call write_tile_data_r0d_fptr(unit,'deadmic_DON_leached',     soil_deadmic_DON_leached_ptr, 'Cumulative dead microbe DON leached out of the column', 'kg/m2')
     call write_tile_data_r0d_fptr(unit,'NO3_leached',     soil_NO3_leached_ptr, 'Cumulative NO3 leached out of the column', 'kgN/m2')
     call write_tile_data_r0d_fptr(unit,'NH4_leached',     soil_NH4_leached_ptr, 'Cumulative NH4 leached out of the column', 'kgN/m2')

     call write_tile_data_r1d_fptr(unit,'leaf_litter_fast_N',soilc_leafLitter_fast_soil_N_ptr,'soilCCohort','Leaf litter fast N','kg/m2')
     call write_tile_data_r1d_fptr(unit,'leaf_litter_slow_N',soilc_leafLitter_slow_soil_N_ptr,'soilCCohort','Leaf litter slow N','kg/m2')
     call write_tile_data_r1d_fptr(unit,'leaf_litter_deadMic_N',soilc_leafLitter_deadMicrobeN_ptr,'soilCCohort','Leaf litter dead microbe N','kg/m2')
     call write_tile_data_r1d_fptr(unit,'leaf_litter_liveMic_N',soilc_leafLitter_livingMicrobeN_ptr,'soilCCohort','Leaf litter live microbe N','kg/m2')
     call write_tile_data_r1d_fptr(unit,'leaf_litter_originalCohortN',soilc_leafLitter_originalLitterN_ptr,'soilCCohort','Leaf litter cohort original N','kg/m2')
     call write_tile_data_r1d_fptr(unit,'leaf_litter_fastProtectedN',soilc_leafLitter_fast_protected_N_ptr,'soilCCohort','Leaf litter fast protected N','kg/m2')
     call write_tile_data_r1d_fptr(unit,'leaf_litter_slowProtectedN',soilc_leafLitter_slow_protected_N_ptr,'soilCCohort','Leaf litter slow protected N','kg/m2')
     call write_tile_data_r1d_fptr(unit,'leaf_litter_deadMicrobeProtectedN',soilc_leafLitter_deadMicrobe_protected_N_ptr,'soilCCohort','Leaf litter dead microbe protected N','kg/m2')

     call write_tile_data_r0d_fptr(unit,'leaf_litter_DON_fast',soilc_leafLitter_fast_DON_ptr,'Leaf litter dissolved fast nitrogen','kg/m2')
     call write_tile_data_r0d_fptr(unit,'leaf_litter_DON_slow',soilc_leafLitter_slow_DON_ptr,'Leaf litter dissolved slow nitrogen','kg/m2')
     call write_tile_data_r0d_fptr(unit,'leaf_litter_DON_deadmic',soilc_leafLitter_deadmicrobe_DON_ptr,'Leaf litter dissolved dead microbe nitrogen','kg/m2')
     call write_tile_data_r0d_fptr(unit,'leaf_litter_NO3',soilc_leafLitter_nitrate_ptr,'Leaf litter nitrate content','kg/m2')
     call write_tile_data_r0d_fptr(unit,'leaf_litter_NH4',soilc_leafLitter_ammonium_ptr,'Leaf litter ammonium content','kg/m2')
     call write_tile_data_r0d_fptr(unit,'leaf_litter_nitrif',soilc_leafLitter_nitrif_ptr,'Leaf litter cumulative nitrification','kg/m2')
     call write_tile_data_r0d_fptr(unit,'leaf_litter_denitrif',soilc_leafLitter_denitrif_ptr,'Leaf litter cumulative denitrification','kg/m2')

     call write_tile_data_r1d_fptr(unit,'fineWood_litter_fast_N',soilc_fineWoodLitter_fast_soil_N_ptr,'soilCCohort','fineWood litter fast N','kg/m2')
     call write_tile_data_r1d_fptr(unit,'fineWood_litter_slow_N',soilc_fineWoodLitter_slow_soil_N_ptr,'soilCCohort','fineWood litter slow N','kg/m2')
     call write_tile_data_r1d_fptr(unit,'fineWood_litter_deadMic_N',soilc_fineWoodLitter_deadMicrobeN_ptr,'soilCCohort','fineWood litter dead microbe N','kg/m2')
     call write_tile_data_r1d_fptr(unit,'fineWood_litter_liveMic_N',soilc_fineWoodLitter_livingMicrobeN_ptr,'soilCCohort','fineWood litter live microbe N','kg/m2')
     call write_tile_data_r1d_fptr(unit,'fineWood_litter_originalCohortN',soilc_fineWoodLitter_originalLitterN_ptr,'soilCCohort','fineWood litter cohort original N','kg/m2')
     call write_tile_data_r1d_fptr(unit,'fineWood_litter_fastProtectedN',soilc_fineWoodLitter_fast_protected_N_ptr,'soilCCohort','fineWood litter fast protected N','kg/m2')
     call write_tile_data_r1d_fptr(unit,'fineWood_litter_slowProtectedN',soilc_fineWoodLitter_slow_protected_N_ptr,'soilCCohort','fineWood litter slow protected N','kg/m2')
     call write_tile_data_r1d_fptr(unit,'fineWood_litter_deadMicrobeProtectedN',soilc_fineWoodLitter_deadMicrobe_protected_N_ptr,'soilCCohort','fineWood litter dead microbe protected N','kg/m2')

     call write_tile_data_r0d_fptr(unit,'fineWood_litter_DON_fast',soilc_fineWoodLitter_fast_DON_ptr,'Fine wood litter dissolved fast nitrogen','kg/m2')
     call write_tile_data_r0d_fptr(unit,'fineWood_litter_DON_slow',soilc_fineWoodLitter_slow_DON_ptr,'Fine wood litter dissolved slow nitrogen','kg/m2')
     call write_tile_data_r0d_fptr(unit,'fineWood_litter_DON_deadmic',soilc_fineWoodLitter_deadmicrobe_DON_ptr,'Fine wood litter dissolved dead microbe nitrogen','kg/m2')
     call write_tile_data_r0d_fptr(unit,'fineWood_litter_NO3',soilc_fineWoodLitter_nitrate_ptr,'Fine wood litter nitrate content','kg/m2')
     call write_tile_data_r0d_fptr(unit,'fineWood_litter_NH4',soilc_fineWoodLitter_ammonium_ptr,'Fine wood litter ammonium content','kg/m2')
     call write_tile_data_r0d_fptr(unit,'fineWood_litter_nitrif',soilc_fineWoodLitter_nitrif_ptr,'Fine wood litter cumulative nitrification','kg/m2')
     call write_tile_data_r0d_fptr(unit,'fineWood_litter_denitrif',soilc_fineWoodLitter_denitrif_ptr,'Fine wood litter cumulative denitrification','kg/m2')

     call write_tile_data_r1d_fptr(unit,'coarseWood_litter_fast_N',soilc_coarseWoodLitter_fast_soil_N_ptr,'soilCCohort','coarseWood litter fast N','kg/m2')
     call write_tile_data_r1d_fptr(unit,'coarseWood_litter_slow_N',soilc_coarseWoodLitter_slow_soil_N_ptr,'soilCCohort','coarseWood litter slow N','kg/m2')
     call write_tile_data_r1d_fptr(unit,'coarseWood_litter_deadMic_N',soilc_coarseWoodLitter_deadMicrobeN_ptr,'soilCCohort','coarseWood litter dead microbe N','kg/m2')
     call write_tile_data_r1d_fptr(unit,'coarseWood_litter_liveMic_N',soilc_coarseWoodLitter_livingMicrobeN_ptr,'soilCCohort','coarseWood litter live microbe N','kg/m2')
     call write_tile_data_r1d_fptr(unit,'coarseWood_litter_originalCohortN',soilc_coarseWoodLitter_originalLitterN_ptr,'soilCCohort','coarseWood litter cohort original N','kg/m2')
     call write_tile_data_r1d_fptr(unit,'coarseWood_litter_fastProtectedN',soilc_coarseWoodLitter_fast_protected_N_ptr,'soilCCohort','coarseWood litter fast protected N','kg/m2')
     call write_tile_data_r1d_fptr(unit,'coarseWood_litter_slowProtectedN',soilc_coarseWoodLitter_slow_protected_N_ptr,'soilCCohort','coarseWood litter slow protected N','kg/m2')
     call write_tile_data_r1d_fptr(unit,'coarseWood_litter_deadMicrobeProtectedN',soilc_coarseWoodLitter_deadMicrobe_protected_N_ptr,'soilCCohort','coarseWood litter dead microbe protected N','kg/m2')


     call write_tile_data_r0d_fptr(unit,'coarseWood_litter_DON_fast',soilc_coarseWoodLitter_fast_DON_ptr,'Coarse wood litter dissolved fast nitrogen','kg/m2')
     call write_tile_data_r0d_fptr(unit,'coarseWood_litter_DON_slow',soilc_coarseWoodLitter_slow_DON_ptr,'Coarse wood litter dissolved slow nitrogen','kg/m2')
     call write_tile_data_r0d_fptr(unit,'coarseWood_litter_DON_deadmic',soilc_coarseWoodLitter_deadmicrobe_DON_ptr,'Coarse wood litter dissolved dead microbe nitrogen','kg/m2')
     call write_tile_data_r0d_fptr(unit,'coarseWood_litter_NO3',soilc_coarseWoodLitter_nitrate_ptr,'Coarse wood litter nitrate content','kg/m2')
     call write_tile_data_r0d_fptr(unit,'coarseWood_litter_NH4',soilc_coarseWoodLitter_ammonium_ptr,'Coarse wood litter ammonium content','kg/m2')
     call write_tile_data_r0d_fptr(unit,'coarseWood_litter_nitrif',soilc_coarseWoodLitter_nitrif_ptr,'Coarse wood litter cumulative nitrification','kg/m2')
     call write_tile_data_r0d_fptr(unit,'coarseWood_litter_denitrif',soilc_coarseWoodLitter_denitrif_ptr,'Coarse wood litter cumulative denitrification','kg/m2')

endif

     call write_tile_data_r1d_fptr(unit,'soil_DOC_fast',soil_fast_DOC_ptr,'zfull','Dissolved fast carbon','kg/m2')
     call write_tile_data_r1d_fptr(unit,'soil_DOC_slow',soil_slow_DOC_ptr,'zfull','Dissolved slow carbon','kg/m2')
     call write_tile_data_r1d_fptr(unit,'soil_DOC_deadmic',soil_deadmicrobe_DOC_ptr,'zfull','Dissolved dead microbe carbon','kg/m2')


     call write_tile_data_r0d_fptr(unit,'fast_DOC_leached',     soil_fast_DOC_leached_ptr, 'Cumulative fast DOC leached out of the column', 'kg/m2')
     call write_tile_data_r0d_fptr(unit,'slow_DOC_leached',     soil_slow_DOC_leached_ptr, 'Cumulative slow DOC leached out of the column', 'kg/m2')
     call write_tile_data_r0d_fptr(unit,'deadmic_DOC_leached',     soil_deadmic_DOC_leached_ptr, 'Cumulative dead microbe DOC leached out of the column', 'kg/m2')


     call write_tile_data_r1d_fptr(unit,'leaf_litter_fast_C',soilc_leafLitter_fast_soil_C_ptr,'soilCCohort','Leaf litter fast C','kg/m2')
     call write_tile_data_r1d_fptr(unit,'leaf_litter_slow_C',soilc_leafLitter_slow_soil_C_ptr,'soilCCohort','Leaf litter slow C','kg/m2')
     call write_tile_data_r1d_fptr(unit,'leaf_litter_deadMic_C',soilc_leafLitter_deadMicrobeC_ptr,'soilCCohort','Leaf litter dead microbe C','kg/m2')
     call write_tile_data_r1d_fptr(unit,'leaf_litter_liveMic_C',soilc_leafLitter_livingMicrobeC_ptr,'soilCCohort','Leaf litter live microbe C','kg/m2')
     call write_tile_data_r1d_fptr(unit,'leaf_litter_CO2',soilc_leafLitter_CO2_ptr,'soilCCohort','Leaf litter CO2 generated','kg/m2')
     call write_tile_data_r1d_fptr(unit,'leaf_litter_Rtot',soilc_leafLitter_Rtot_ptr,'soilCCohort','Leaf litter total degradation','kg/m2')
     call write_tile_data_r1d_fptr(unit,'leaf_litter_originalCohortC',soilc_leafLitter_originalLitterC_ptr,'soilCCohort','Leaf litter cohort original carbon','kg/m2')
     call write_tile_data_r1d_fptr(unit,'leaf_litter_fastProtectedC',soilc_leafLitter_fast_protected_C_ptr,'soilCCohort','Leaf litter fast protected C','kg/m2')
     call write_tile_data_r1d_fptr(unit,'leaf_litter_slowProtectedC',soilc_leafLitter_slow_protected_C_ptr,'soilCCohort','Leaf litter slow protected C','kg/m2')
     call write_tile_data_r1d_fptr(unit,'leaf_litter_deadMicrobeProtectedC',soilc_leafLitter_deadMicrobe_protected_C_ptr,'soilCCohort','Leaf litter dead microbe protected C','kg/m2')


     call write_tile_data_r0d_fptr(unit,'leaf_litter_DOC_fast',soilc_leafLitter_fast_DOC_ptr,'Dissolved leaf litter fast carbon','kg/m2')
     call write_tile_data_r0d_fptr(unit,'leaf_litter_DOC_slow',soilc_leafLitter_slow_DOC_ptr,'Dissolved leaf litter slow carbon','kg/m2')
     call write_tile_data_r0d_fptr(unit,'leaf_litter_DOC_deadmic',soilc_leafLitter_deadmicrobe_DOC_ptr,'Dissolved leaf litter dead microbe carbon','kg/m2')



     call write_tile_data_r1d_fptr(unit,'fineWood_litter_fast_C',soilc_fineWoodLitter_fast_soil_C_ptr,'soilCCohort','Fine wood litter fast C','kg/m2')
     call write_tile_data_r1d_fptr(unit,'fineWood_litter_slow_C',soilc_fineWoodLitter_slow_soil_C_ptr,'soilCCohort','Fine wood litter slow C','kg/m2')
     call write_tile_data_r1d_fptr(unit,'fineWood_litter_deadMic_C',soilc_fineWoodLitter_deadMicrobeC_ptr,'soilCCohort','Fine wood litter dead microbe C','kg/m2')
     call write_tile_data_r1d_fptr(unit,'fineWood_litter_liveMic_C',soilc_fineWoodLitter_livingMicrobeC_ptr,'soilCCohort','Fine wood litter live microbe C','kg/m2')
     call write_tile_data_r1d_fptr(unit,'fineWood_litter_CO2',soilc_fineWoodLitter_CO2_ptr,'soilCCohort','Fine wood litter CO2 generated','kg/m2')
     call write_tile_data_r1d_fptr(unit,'fineWood_litter_Rtot',soilc_fineWoodLitter_Rtot_ptr,'soilCCohort','Fine wood litter total degradation','kg/m2')
     call write_tile_data_r1d_fptr(unit,'fineWood_litter_originalCohortC',soilc_fineWoodLitter_originalLitterC_ptr,'soilCCohort','Fine wood litter cohort original carbon','kg/m2')
     call write_tile_data_r1d_fptr(unit,'fineWood_litter_fastProtectedC',soilc_fineWoodLitter_fast_protected_C_ptr,'soilCCohort','Fine wood litter fast protected C','kg/m2')
     call write_tile_data_r1d_fptr(unit,'fineWood_litter_slowProtectedC',soilc_fineWoodLitter_slow_protected_C_ptr,'soilCCohort','Fine wood litter slow protected C','kg/m2')
     call write_tile_data_r1d_fptr(unit,'fineWood_litter_deadMicrobeProtectedC',soilc_fineWoodLitter_deadMicrobe_protected_C_ptr,'soilCCohort','Fine wood litter dead microbe protected C','kg/m2')



     call write_tile_data_r0d_fptr(unit,'fineWood_litter_DOC_fast',soilc_fineWoodLitter_fast_DOC_ptr,'Dissolved fine wood litter fast carbon','kg/m2')
     call write_tile_data_r0d_fptr(unit,'fineWood_litter_DOC_slow',soilc_fineWoodLitter_slow_DOC_ptr,'Dissolved fine wood litter slow carbon','kg/m2')
     call write_tile_data_r0d_fptr(unit,'fineWood_litter_DOC_deadmic',soilc_fineWoodLitter_deadmicrobe_DOC_ptr,'Dissolved fine wood litter dead microbe carbon','kg/m2')


     call write_tile_data_r1d_fptr(unit,'coarseWood_litter_fast_C',soilc_coarseWoodLitter_fast_soil_C_ptr,'soilCCohort','Coarse wood litter fast C','kg/m2')
     call write_tile_data_r1d_fptr(unit,'coarseWood_litter_slow_C',soilc_coarseWoodLitter_slow_soil_C_ptr,'soilCCohort','Coarse wood litter slow C','kg/m2')
     call write_tile_data_r1d_fptr(unit,'coarseWood_litter_deadMic_C',soilc_coarseWoodLitter_deadMicrobeC_ptr,'soilCCohort','Coarse wood litter dead microbe C','kg/m2')
     call write_tile_data_r1d_fptr(unit,'coarseWood_litter_liveMic_C',soilc_coarseWoodLitter_livingMicrobeC_ptr,'soilCCohort','Coarse wood litter live microbe C','kg/m2')
     call write_tile_data_r1d_fptr(unit,'coarseWood_litter_CO2',soilc_coarseWoodLitter_CO2_ptr,'soilCCohort','Coarse wood litter CO2 generated','kg/m2')
     call write_tile_data_r1d_fptr(unit,'coarseWood_litter_Rtot',soilc_coarseWoodLitter_Rtot_ptr,'soilCCohort','Coarse wood litter total degradation','kg/m2')
     call write_tile_data_r1d_fptr(unit,'coarseWood_litter_originalCohortC',soilc_coarseWoodLitter_originalLitterC_ptr,'soilCCohort','Coarse wood litter cohort original carbon','kg/m2')
     call write_tile_data_r1d_fptr(unit,'coarseWood_litter_fastProtectedC',soilc_coarseWoodLitter_fast_protected_C_ptr,'soilCCohort','Coarse wood litter fast protected C','kg/m2')
     call write_tile_data_r1d_fptr(unit,'coarseWood_litter_slowProtectedC',soilc_coarseWoodLitter_slow_protected_C_ptr,'soilCCohort','Coarse wood litter slow protected C','kg/m2')
     call write_tile_data_r1d_fptr(unit,'coarseWood_litter_deadMicrobeProtectedC',soilc_coarseWoodLitter_deadMicrobe_protected_C_ptr,'soilCCohort','Coarse wood litter dead microbe protected C','kg/m2')


     call write_tile_data_r0d_fptr(unit,'coarseWood_litter_DOC_fast',soilc_coarseWoodLitter_fast_DOC_ptr,'Dissolved coarse wood litter fast carbon','kg/m2')
     call write_tile_data_r0d_fptr(unit,'coarseWood_litter_DOC_slow',soilc_coarseWoodLitter_slow_DOC_ptr,'Dissolved coarse wood litter slow carbon','kg/m2')
     call write_tile_data_r0d_fptr(unit,'coarseWood_litter_DOC_deadmic',soilc_coarseWoodLitter_deadmicrobe_DOC_ptr,'Dissolved coarse wood litter dead microbe carbon','kg/m2')



     call write_tile_data_i1d_fptr_all(unit,'is_peat',soil_is_peat_ptr,'zfull','Is layer peat?','Boolean')
  case default
     call error_mesg('save_soil_restart','soil_carbon_option is invalid. This should never happen. Contact developer', FATAL)
  end select

  ! close file
  __NF_ASRT__(nf_close(unit))

  if (write_soil_carbon_restart) then
     call create_tile_out_file(unit,'RESTART/'//trim(timestamp)//'soil_carbon.res.nc', &
          lnd%coord_glon, lnd%coord_glat, soil_tile_exists, tile_dim_length )
     ! in addition, define vertical coordinate
     if (mpp_pe()==lnd%io_pelist(1)) then
        __NF_ASRT__(nfu_def_dim(unit,'zfull',zfull(1:num_l),'full level','m'))
        __NF_ASRT__(nfu_put_att(unit,'zfull','positive','down'))
     endif
     call sync_nc_files(unit)

     call write_tile_data_r1d_fptr(unit,'asoil_in',soil_asoil_in_ptr,'zfull','aerobic activity modifier', 'unitless')
     call write_tile_data_r1d_fptr(unit,'fsc_in',soil_fsc_in_ptr,'zfull','fast soil carbon input', 'kg C/m2')
     call write_tile_data_r1d_fptr(unit,'ssc_in',soil_ssc_in_ptr,'zfull','slow soil carbon input', 'kg C/m2')
     if (soil_carbon_option == SOILC_CORPSE .or. soil_carbon_option == SOILC_CORPSE_N) then
         call write_tile_data_r1d_fptr(unit,'fsn_in',soil_fsn_in_ptr,'zfull','fast soil nitrogen input', 'kg N/m2')
         call write_tile_data_r1d_fptr(unit,'ssn_in',soil_ssn_in_ptr,'zfull','slow soil nitrogen input', 'kg N/m2')
        call write_tile_data_r1d_fptr(unit,'deadmic_C_in',soil_deadmic_C_in_ptr,'zfull','dead microbe soil carbon input', 'kg C/m2')
        call write_tile_data_r1d_fptr(unit,'deadmic_N_in',soil_deadmic_N_in_ptr,'zfull','dead microbe soil nitrogen input', 'kg N/m2')
        call write_tile_data_r1d_fptr(unit,'fast_protected_C_in',soil_fast_protected_C_in_ptr,'zfull','protected fast soil carbon input', 'kg C/m2')
        call write_tile_data_r1d_fptr(unit,'slow_protected_C_in',soil_slow_protected_C_in_ptr,'zfull','protected slow soil carbon input', 'kg C/m2')
        call write_tile_data_r1d_fptr(unit,'deadmic_protected_C_in',soil_deadmic_protected_C_in_ptr,'zfull','protected dead microbe soil carbon input', 'kg C/m2')
        call write_tile_data_r1d_fptr(unit,'fast_protected_N_in',soil_fast_protected_N_in_ptr,'zfull','protected fast soil nitrogen input', 'kg N/m2')
        call write_tile_data_r1d_fptr(unit,'slow_protected_N_in',soil_slow_protected_N_in_ptr,'zfull','protected slow soil nitrogen input', 'kg N/m2')
        call write_tile_data_r1d_fptr(unit,'deadmic_protected_N_in',soil_deadmic_protected_N_in_ptr,'zfull','protected dead microbe soil nitrogen input', 'kg N/m2')

        call write_tile_data_r1d_fptr(unit,'fast_C_turnover_accumulated',soil_fast_C_turnover_accumulated_ptr,'zfull','fast soil carbon turnover', 'year-1')
        call write_tile_data_r1d_fptr(unit,'slow_C_turnover_accumulated',soil_slow_C_turnover_accumulated_ptr,'zfull','slow soil carbon turnover', 'year-1')
        call write_tile_data_r1d_fptr(unit,'deadmic_C_turnover_accumulated',soil_deadmic_C_turnover_accumulated_ptr,'zfull','dead microbe soil carbon turnover', 'year-1')
        call write_tile_data_r1d_fptr(unit,'fast_protected_C_turnover_accumulated',soil_fast_protected_C_turnover_accumulated_ptr,'zfull','fast protected soil carbon turnover', 'year-1')
        call write_tile_data_r1d_fptr(unit,'slow_protected_C_turnover_accumulated',soil_slow_protected_C_turnover_accumulated_ptr,'zfull','slow protectedsoil carbon turnover', 'year-1')
        call write_tile_data_r1d_fptr(unit,'deadmic_protected_C_turnover_accumulated',soil_deadmic_protected_C_turnover_accumulated_ptr,'zfull','dead microbe protected soil carbon turnover', 'year-1')

        call write_tile_data_r1d_fptr(unit,'fast_N_turnover_accumulated',soil_fast_N_turnover_accumulated_ptr,'zfull','fast soil carbon turnover', 'year-1')
        call write_tile_data_r1d_fptr(unit,'slow_N_turnover_accumulated',soil_slow_N_turnover_accumulated_ptr,'zfull','slow soil carbon turnover', 'year-1')
        call write_tile_data_r1d_fptr(unit,'deadmic_N_turnover_accumulated',soil_deadmic_N_turnover_accumulated_ptr,'zfull','dead microbe soil carbon turnover', 'year-1')
        call write_tile_data_r1d_fptr(unit,'fast_protected_N_turnover_accumulated',soil_fast_protected_N_turnover_accumulated_ptr,'zfull','fast protected soil carbon turnover', 'year-1')
        call write_tile_data_r1d_fptr(unit,'slow_protected_N_turnover_accumulated',soil_slow_protected_N_turnover_accumulated_ptr,'zfull','slow protectedsoil carbon turnover', 'year-1')
        call write_tile_data_r1d_fptr(unit,'deadmic_protected_N_turnover_accumulated',soil_deadmic_protected_N_turnover_accumulated_ptr,'zfull','dead microbe protected soil carbon turnover', 'year-1')


        call write_tile_data_r0d_fptr(unit,'leaflitter_fast_C_turnover_accumulated',soil_leaflitter_fast_C_turnover_accumulated_ptr,'fast leaf litter carbon turnover', 'year-1')
        call write_tile_data_r0d_fptr(unit,'leaflitter_slow_C_turnover_accumulated',soil_leaflitter_slow_C_turnover_accumulated_ptr,'slow leaf litter carbon turnover', 'year-1')
        call write_tile_data_r0d_fptr(unit,'leaflitter_deadmic_C_turnover_accumulated',soil_leaflitter_deadmic_C_turnover_accumulated_ptr,'dead microbe leaf litter carbon turnover', 'year-1')
        call write_tile_data_r0d_fptr(unit,'leaflitter_fsc_in',soil_leaflitter_fsc_in_ptr,'fast leaf litter carbon input', 'kg C/m2')
        call write_tile_data_r0d_fptr(unit,'leaflitter_ssc_in',soil_leaflitter_ssc_in_ptr,'slow leaf litter carbon input', 'kg C/m2')
        call write_tile_data_r0d_fptr(unit,'leaflitter_deadmic_C_in',soil_leaflitter_deadmic_C_in_ptr,'dead microbe leaf litter carbon input', 'kg C/m2')

        call write_tile_data_r0d_fptr(unit,'leaflitter_fast_N_turnover_accumulated',soil_leaflitter_fast_N_turnover_accumulated_ptr,'fast leaf litter nitrogen turnover', 'year-1')
        call write_tile_data_r0d_fptr(unit,'leaflitter_slow_N_turnover_accumulated',soil_leaflitter_slow_N_turnover_accumulated_ptr,'slow leaf litter nitrogen turnover', 'year-1')
        call write_tile_data_r0d_fptr(unit,'leaflitter_deadmic_N_turnover_accumulated',soil_leaflitter_deadmic_N_turnover_accumulated_ptr,'dead microbe leaf litter nitrogen turnover', 'year-1')
        call write_tile_data_r0d_fptr(unit,'leaflitter_fsn_in',soil_leaflitter_fsn_in_ptr,'fast leaf litter nitrogen input', 'kg N/m2')
        call write_tile_data_r0d_fptr(unit,'leaflitter_ssn_in',soil_leaflitter_ssn_in_ptr,'slow leaf litter nitrogen input', 'kg N/m2')
        call write_tile_data_r0d_fptr(unit,'leaflitter_deadmic_N_in',soil_leaflitter_deadmic_N_in_ptr,'dead microbe leaf litter nitrogen input', 'kg N/m2')

        call write_tile_data_r0d_fptr(unit,'finewoodlitter_fast_C_turnover_accumulated',soil_finewoodlitter_fast_C_turnover_accumulated_ptr,'fast fine wood litter carbon turnover', 'year-1')
        call write_tile_data_r0d_fptr(unit,'finewoodlitter_slow_C_turnover_accumulated',soil_finewoodlitter_slow_C_turnover_accumulated_ptr,'slow fine wood litter carbon turnover', 'year-1')
        call write_tile_data_r0d_fptr(unit,'finewoodlitter_deadmic_C_turnover_accumulated',soil_finewoodlitter_deadmic_C_turnover_accumulated_ptr,'dead microbe fine wood litter carbon turnover', 'year-1')
        call write_tile_data_r0d_fptr(unit,'finewoodlitter_fsc_in',soil_finewoodlitter_fsc_in_ptr,'fast fine wood litter carbon input', 'kg C/m2')
        call write_tile_data_r0d_fptr(unit,'finewoodlitter_ssc_in',soil_finewoodlitter_ssc_in_ptr,'slow fine wood litter carbon input', 'kg C/m2')
        call write_tile_data_r0d_fptr(unit,'finewoodlitter_deadmic_C_in',soil_finewoodlitter_deadmic_C_in_ptr,'dead microbe fine wood litter carbon input', 'kg C/m2')

        call write_tile_data_r0d_fptr(unit,'finewoodlitter_fast_N_turnover_accumulated',soil_finewoodlitter_fast_N_turnover_accumulated_ptr,'fast fine wood litter nitrogen turnover', 'year-1')
        call write_tile_data_r0d_fptr(unit,'finewoodlitter_slow_N_turnover_accumulated',soil_finewoodlitter_slow_N_turnover_accumulated_ptr,'slow fine wood litter nitrogen turnover', 'year-1')
        call write_tile_data_r0d_fptr(unit,'finewoodlitter_deadmic_N_turnover_accumulated',soil_finewoodlitter_deadmic_N_turnover_accumulated_ptr,'dead microbe fine wood litter nitrogen turnover', 'year-1')
        call write_tile_data_r0d_fptr(unit,'finewoodlitter_fsn_in',soil_finewoodlitter_fsn_in_ptr,'fast fine wood litter nitrogen input', 'kg N/m2')
        call write_tile_data_r0d_fptr(unit,'finewoodlitter_ssn_in',soil_finewoodlitter_ssn_in_ptr,'slow fine wood litter nitrogen input', 'kg N/m2')
        call write_tile_data_r0d_fptr(unit,'finewoodlitter_deadmic_N_in',soil_finewoodlitter_deadmic_N_in_ptr,'dead microbe fine wood litter nitrogen input', 'kg N/m2')

        call write_tile_data_r0d_fptr(unit,'coarsewoodlitter_fast_C_turnover_accumulated',soil_coarsewoodlitter_fast_C_turnover_accumulated_ptr,'fast coarse wood litter carbon turnover', 'year-1')
        call write_tile_data_r0d_fptr(unit,'coarsewoodlitter_slow_C_turnover_accumulated',soil_coarsewoodlitter_slow_C_turnover_accumulated_ptr,'slow coarse wood litter carbon turnover', 'year-1')
        call write_tile_data_r0d_fptr(unit,'coarsewoodlitter_deadmic_C_turnover_accumulated',soil_coarsewoodlitter_deadmic_C_turnover_accumulated_ptr,'dead microbe coarse wood litter carbon turnover', 'year-1')
        call write_tile_data_r0d_fptr(unit,'coarsewoodlitter_fsc_in',soil_coarsewoodlitter_fsc_in_ptr,'fast coarse wood litter carbon input', 'kg C/m2')
        call write_tile_data_r0d_fptr(unit,'coarsewoodlitter_ssc_in',soil_coarsewoodlitter_ssc_in_ptr,'slow coarse wood litter carbon input', 'kg C/m2')
        call write_tile_data_r0d_fptr(unit,'coarsewoodlitter_deadmic_C_in',soil_coarsewoodlitter_deadmic_C_in_ptr,'dead microbe coarse wood litter carbon input', 'kg C/m2')

        call write_tile_data_r0d_fptr(unit,'coarsewoodlitter_fast_N_turnover_accumulated',soil_coarsewoodlitter_fast_N_turnover_accumulated_ptr,'fast coarse wood litter nitrogen turnover', 'year-1')
        call write_tile_data_r0d_fptr(unit,'coarsewoodlitter_slow_N_turnover_accumulated',soil_coarsewoodlitter_slow_N_turnover_accumulated_ptr,'slow coarse wood litter nitrogen turnover', 'year-1')
        call write_tile_data_r0d_fptr(unit,'coarsewoodlitter_deadmic_N_turnover_accumulated',soil_coarsewoodlitter_deadmic_N_turnover_accumulated_ptr,'dead microbe coarse wood litter nitrogen turnover', 'year-1')
        call write_tile_data_r0d_fptr(unit,'coarsewoodlitter_fsn_in',soil_coarsewoodlitter_fsn_in_ptr,'fast coarse wood litter nitrogen input', 'kg N/m2')
        call write_tile_data_r0d_fptr(unit,'coarsewoodlitter_ssn_in',soil_coarsewoodlitter_ssn_in_ptr,'slow coarse wood litter nitrogen input', 'kg N/m2')
        call write_tile_data_r0d_fptr(unit,'coarsewoodlitter_deadmic_N_in',soil_coarsewoodlitter_deadmic_N_in_ptr,'dead microbe coarse wood litter nitrogen input', 'kg N/m2')


     endif
     __NF_ASRT__(nf_close(unit))
  endif


end subroutine save_soil_restart

! ============================================================================
subroutine save_soil_restart_new (tile_dim_length, timestamp)
  integer, intent(in) :: tile_dim_length ! length of tile dim. in the output file
  character(*), intent(in) :: timestamp ! timestamp to add to the file name

  ! ---- local vars ----------------------------------------------------------
  integer :: unit            ! restart file i/o unit
  character(267) :: fname
  type(restart_file_type) :: soil_restart, soil_carbon_restart ! restart files i/o object
  integer, allocatable :: idx(:)
  real, allocatable, dimension(:,:) :: temp, wl, ws, groundwater, groundwater_T, fsc, ssc! (tile_index, zfull)
  real, allocatable, dimension(:)   :: uptake_T
  real, allocatable, dimension(:,:) :: asoil_in, fsc_in, ssc_in

  real, allocatable, dimension(:,:) :: fsn_in, ssn_in

  real, allocatable, dimension(:,:,:) :: & ! (tile_index, zfull, soilCCohort)
                  fast_soil_C,           &
                  slow_soil_C,           &
                  deadMic,               &
                  fastProtectedC,        &
                  slowProtectedC,        &
                  deadMicrobeProtectedC, &
                  liveMic,               &
                  CO2,                   &
                  Rtot,                  &
                  originalCohortC

  real, allocatable, dimension(:,:,:) :: & ! (tile_index, zfull, soilCCohort)
                  fast_soil_N,           &
                  slow_soil_N,           &
                  deadMic_N,               &
                  fastProtectedN,        &
                  slowProtectedN,        &
                  deadMicrobeProtectedN, &
                  liveMic_N,               &
                  originalCohortN

  real, allocatable, dimension(:,:) :: & ! (tile_index, zfull)
                  soil_DOC_fast,       &
                  soil_DOC_slow,       &
                  soil_DOC_deadmic,    &
                  deadmic_in ,          &
                  fast_protected_in ,   &
                  slow_protected_in ,   &
                  deadmic_protected_in , &
                  fast_turnover_accumulated , &
                  slow_turnover_accumulated , &
                  deadmic_turnover_accumulated , &
                  fast_protected_turnover_accumulated , &
                  slow_protected_turnover_accumulated , &
                  deadmic_protected_turnover_accumulated

  real, allocatable, dimension(:,:) :: & ! (tile_index, zfull)
                  soil_DON_fast,       &
                  soil_DON_slow,       &
                  soil_DON_deadmic,    &
                  soil_NO3,            &
                  soil_NH4,            &
                  soil_nitrif,         &
                  soil_denitrif,       &
                  deadmic_N_in ,          &
                  fast_protected_N_in ,   &
                  slow_protected_N_in ,   &
                  deadmic_protected_N_in , &
                  fast_N_turnover_accumulated , &
                  slow_N_turnover_accumulated , &
                  deadmic_N_turnover_accumulated , &
                  fast_protected_N_turnover_accumulated , &
                  slow_protected_N_turnover_accumulated , &
                  deadmic_protected_N_turnover_accumulated

  real, allocatable, dimension(:,:) ::                   & ! (tile_index, soilCCohort)
                  leaf_litter_fast_C,                    &
                  leaf_litter_slow_C,                    &
                  leaf_litter_deadMic_C,                 &
                  leaf_litter_liveMic_C,                 &
                  leaf_litter_CO2,                       &
                  leaf_litter_Rtot,                      &
                  leaf_litter_originalCohortC,           &
                  leaf_litter_fastProtectedC,            &
                  leaf_litter_slowProtectedC,            &
                  leaf_litter_deadMicrobeProtectedC,     &
                  fineWood_litter_fast_C,                &
                  fineWood_litter_slow_C,                &
                  fineWood_litter_deadMic_C,             &
                  fineWood_litter_liveMic_C,             &
                  fineWood_litter_CO2,                   &
                  fineWood_litter_Rtot,                  &
                  fineWood_litter_originalCohortC,       &
                  fineWood_litter_fastProtectedC,        &
                  fineWood_litter_slowProtectedC,        &
                  fineWood_litter_deadMicrobeProtectedC, &
                  coarseWood_litter_fast_C,              &
                  coarseWood_litter_slow_C,              &
                  coarseWood_litter_deadMic_C,           &
                  coarseWood_litter_liveMic_C,           &
                  coarseWood_litter_CO2,                 &
                  coarseWood_litter_Rtot,                &
                  coarseWood_litter_originalCohortC,     &
                  coarseWood_litter_fastProtectedC,      &
                  coarseWood_litter_slowProtectedC,      &
                  coarseWood_litter_deadMicrobeProtectedC


  real, allocatable, dimension(:,:) ::                   & ! (tile_index, soilCCohort)
                  leaf_litter_fast_N,                    &
                  leaf_litter_slow_N,                    &
                  leaf_litter_deadMic_N,                 &
                  leaf_litter_liveMic_N,                 &
                  leaf_litter_originalCohortN,           &
                  leaf_litter_fastProtectedN,            &
                  leaf_litter_slowProtectedN,            &
                  leaf_litter_deadMicrobeProtectedN,     &
                  fineWood_litter_fast_N,                &
                  fineWood_litter_slow_N,                &
                  fineWood_litter_deadMic_N,             &
                  fineWood_litter_liveMic_N,             &
                  fineWood_litter_originalCohortN,       &
                  fineWood_litter_fastProtectedN,        &
                  fineWood_litter_slowProtectedN,        &
                  fineWood_litter_deadMicrobeProtectedN, &
                  coarseWood_litter_fast_N,              &
                  coarseWood_litter_slow_N,              &
                  coarseWood_litter_deadMic_N,           &
                  coarseWood_litter_liveMic_N,           &
                  coarseWood_litter_originalCohortN,     &
                  coarseWood_litter_fastProtectedN,      &
                  coarseWood_litter_slowProtectedN,      &
                  coarseWood_litter_deadMicrobeProtectedN

  real, allocatable, dimension(:) :: & ! (tile_index)
        fast_DOC_leached ,            &
        slow_DOC_leached ,            &
        deadmic_DOC_leached ,         &
        leaf_litter_DOC_fast ,        &
        leaf_litter_DOC_slow ,        &
        leaf_litter_DOC_deadmic ,     &
        fineWood_litter_DOC_fast ,    &
        fineWood_litter_DOC_slow ,    &
        fineWood_litter_DOC_deadmic , &
        coarseWood_litter_DOC_fast ,  &
        coarseWood_litter_DOC_slow ,  &
        coarseWood_litter_DOC_deadmic ,                 &
        leaflitter_fast_turnover_accumulated ,          &
        leaflitter_slow_turnover_accumulated ,          &
        leaflitter_deadmic_turnover_accumulated ,       &
        leaflitter_fsc_in ,                             &
        leaflitter_ssc_in ,                             &
        leaflitter_deadmic_in ,                         &
        finewoodlitter_fast_turnover_accumulated ,      &
        finewoodlitter_slow_turnover_accumulated ,      &
        finewoodlitter_deadmic_turnover_accumulated ,   &
        finewoodlitter_fsc_in ,                         &
        finewoodlitter_ssc_in ,                         &
        finewoodlitter_deadmic_in ,                     &
        coarsewoodlitter_fast_turnover_accumulated ,    &
        coarsewoodlitter_slow_turnover_accumulated ,    &
        coarsewoodlitter_deadmic_turnover_accumulated , &
        coarsewoodlitter_fsc_in ,                       &
        coarsewoodlitter_ssc_in ,                       &
        coarsewoodlitter_deadmic_in


real, allocatable, dimension(:) :: & ! (tile_index)
      fast_DON_leached ,            &
      slow_DON_leached ,            &
      deadmic_DON_leached ,         &
      NO3_leached ,                 &
      NH4_leached ,                 &
      leaf_litter_DON_fast ,        &
      leaf_litter_DON_slow ,        &
      leaf_litter_DON_deadmic ,     &
      leaf_litter_NO3,              &
      leaf_litter_NH4,              &
      leaf_litter_nitrif,           &
      leaf_litter_denitrif,         &
      fineWood_litter_DON_fast ,    &
      fineWood_litter_DON_slow ,    &
      fineWood_litter_DON_deadmic , &
      fineWood_litter_NO3,              &
      fineWood_litter_NH4,              &
      fineWood_litter_nitrif,           &
      fineWood_litter_denitrif,         &
      coarseWood_litter_DON_fast ,  &
      coarseWood_litter_DON_slow ,  &
      coarseWood_litter_DON_deadmic ,                 &
      coarseWood_litter_NO3,              &
      coarseWood_litter_NH4,              &
      coarseWood_litter_nitrif,           &
      coarseWood_litter_denitrif,         &
      leaflitter_fast_N_turnover_accumulated ,          &
      leaflitter_slow_N_turnover_accumulated ,          &
      leaflitter_deadmic_N_turnover_accumulated ,       &
      leaflitter_fsn_in ,                             &
      leaflitter_ssn_in ,                             &
      leaflitter_deadmic_N_in ,                         &
      finewoodlitter_fast_N_turnover_accumulated ,      &
      finewoodlitter_slow_N_turnover_accumulated ,      &
      finewoodlitter_deadmic_N_turnover_accumulated ,   &
      finewoodlitter_fsn_in ,                         &
      finewoodlitter_ssn_in ,                         &
      finewoodlitter_deadmic_N_in ,                     &
      coarsewoodlitter_fast_N_turnover_accumulated ,    &
      coarsewoodlitter_slow_N_turnover_accumulated ,    &
      coarsewoodlitter_deadmic_N_turnover_accumulated , &
      coarsewoodlitter_fsn_in ,                       &
      coarsewoodlitter_ssn_in ,                       &
      coarsewoodlitter_deadmic_N_in

  integer, allocatable, dimension(:,:) :: is_peat ! (tile_index, zfull)
  integer :: id_restart, isize, nccoh

  call error_mesg('soil_end','writing new format NetCDF restart',NOTE)

  ! must set domain so that io_domain is available
  call set_domain(lnd%domain)
! Note that fname is updated for tile & rank numbers during file creation
  fname = trim(timestamp)//'soil.res.nc'

  select case(soil_carbon_option)
  case (SOILC_CENTURY, SOILC_CENTURY_BY_LAYER)
     call create_tile_out_file(soil_restart,idx,fname,soil_tile_exists,tile_dim_length,zfull(1:num_l))
 case (SOILC_CORPSE, SOILC_CORPSE_N)
     call create_tile_out_file(soil_restart,idx,fname,soil_tile_exists,tile_dim_length,zfull(1:num_l),soilCCohort_data=soilCCohort_data)
  case default
     call error_mesg('save_soil_restart','soil_carbon_option is invalid. This should never happen. Contact developer', FATAL)
  end select

  isize = size(idx)
  allocate(temp(isize,num_l),wl(isize,num_l),ws(isize,num_l),groundwater(isize,num_l),groundwater_T(isize,num_l),uptake_T(isize))

 ! Output data provides signature
  call gather_tile_data(soil_T_ptr,idx,temp)
  id_restart = register_restart_field(soil_restart,fname,'temp',temp,compressed=.true., &
                                      longname='soil temperature',units='degrees_K')

  call gather_tile_data(soil_wl_ptr,idx,wl)
  id_restart = register_restart_field(soil_restart,fname,'wl',wl,compressed=.true., &
                                      longname='liquid water content',units='kg/m2')

  call gather_tile_data(soil_ws_ptr,idx,ws)
  id_restart = register_restart_field(soil_restart,fname,'ws',ws,compressed=.true., &
                                      longname='solid water content',units='kg/m2')

  call gather_tile_data(soil_groundwater_ptr,idx,groundwater)
  id_restart = register_restart_field(soil_restart,fname,'groundwater',groundwater,compressed=.true., &
                                      longname='groundwater',units='kg/m2')

  call gather_tile_data(soil_groundwater_T_ptr,idx,groundwater_T)
  id_restart = register_restart_field(soil_restart,fname,'groundwater_T',groundwater_T,compressed=.true., &
                                      longname='groundwater temperature',units='degrees_K')

  call gather_tile_data(soil_uptake_T_ptr,idx,uptake_T)
  id_restart = register_restart_field(soil_restart,fname,'uptake_T',uptake_T, &
                                      longname='temperature of transpiring water',units='degrees_K')

  select case(soil_carbon_option)
  case (SOILC_CENTURY, SOILC_CENTURY_BY_LAYER)
     allocate(fsc(isize,num_l), ssc(isize,num_l))

     call gather_tile_data(soil_fast_soil_C_ptr,idx,fsc)
     id_restart = register_restart_field(soil_restart,fname,'fsc',fsc,compressed=.true., &
                                         longname='fast soil carbon',units='kg C/m2')

     call gather_tile_data(soil_slow_soil_C_ptr,idx,ssc)
     id_restart = register_restart_field(soil_restart,fname,'ssc',ssc,compressed=.true., &
                                         longname='slow soil carbon',units='kg C/m2')
  case (SOILC_CORPSE, SOILC_CORPSE_N)

     ! Register the fields that have the soil carbon cohort axis
     allocate(fast_soil_C(isize,num_l,soilMaxCohorts))
     call gather_tile_data(soilc_fast_soil_C_ptr,idx,fast_soil_C)
     id_restart = register_restart_field(soil_restart,fname,'fast_soil_C',fast_soil_C,compressed=.true., &
                                         longname='Fast soil carbon',units='kg/m2')
     allocate(slow_soil_C(isize,num_l,soilMaxCohorts))
     call gather_tile_data(soilc_slow_soil_C_ptr,idx,slow_soil_C)
     id_restart = register_restart_field(soil_restart,fname,'slow_soil_C',slow_soil_C,compressed=.true., &
                                         longname='Slow soil carbon',units='kg/m2')
     allocate(deadMic(isize,num_l,soilMaxCohorts))
     call gather_tile_data(soilc_deadMicrobeC_ptr,idx,deadMic)
     id_restart = register_restart_field(soil_restart,fname,'deadMic',deadMic,compressed=.true., &
                                         longname='Dead microbe carbon',units='kg/m2')



     allocate(fastProtectedC(isize,num_l,soilMaxCohorts))
     call gather_tile_data(soilc_fast_protected_C_ptr,idx,fastProtectedC)
     id_restart = register_restart_field(soil_restart,fname,'fastProtectedC',fastProtectedC,compressed=.true., &
                                         longname='Protected fast carbon',units='kg/m2')
     allocate(slowProtectedC(isize,num_l,soilMaxCohorts))
     call gather_tile_data(soilc_slow_protected_C_ptr,idx,slowProtectedC)
     id_restart = register_restart_field(soil_restart,fname,'slowProtectedC',slowProtectedC,compressed=.true., &
                                         longname='Protected slow carbon',units='kg/m2')
     allocate(deadMicrobeProtectedC(isize,num_l,soilMaxCohorts))
     call gather_tile_data(soilc_deadMicrobe_protected_C_ptr,idx,deadMicrobeProtectedC)
     id_restart = register_restart_field(soil_restart,fname,'deadMicrobeProtectedC',deadMicrobeProtectedC,compressed=.true., &
                                         longname='Protected dead microbe carbon',units='kg/m2')



     allocate(liveMic(isize,num_l,soilMaxCohorts))
     call gather_tile_data(soilc_livingMicrobeC_ptr,idx,liveMic)
     id_restart = register_restart_field(soil_restart,fname,'liveMic',liveMic,compressed=.true., &
                                         longname='Living microbial carbon',units='kg/m2')


     allocate(CO2(isize,num_l,soilMaxCohorts))
     call gather_tile_data(soilc_CO2_ptr,idx,CO2)
     id_restart = register_restart_field(soil_restart,fname,'CO2',CO2,compressed=.true., &
                                         longname='Cohort CO2 generated',units='kg/m2')
     allocate(Rtot(isize,num_l,soilMaxCohorts))
     call gather_tile_data(soilc_Rtot_ptr,idx,Rtot)
     id_restart = register_restart_field(soil_restart,fname,'Rtot',Rtot,compressed=.true., &
                                         longname='Total degradation',units='kg/m2')
     allocate(originalCohortC(isize,num_l,soilMaxCohorts))
     call gather_tile_data(soilc_originalLitterC_ptr,idx,originalCohortC)
     id_restart = register_restart_field(soil_restart,fname,'originalCohortC',originalCohortC,compressed=.true., &
                                         longname='Cohort original carbon',units='kg/m2')




     allocate(soil_DOC_fast(isize,num_l))
     call gather_tile_data(soil_fast_DOC_ptr,idx,soil_DOC_fast)
     id_restart = register_restart_field(soil_restart,fname,'soil_DOC_fast',soil_DOC_fast,compressed=.true., &
                                         longname='Dissolved fast carbon',units='kg/m2')
     allocate(soil_DOC_slow(isize,num_l))
     call gather_tile_data(soil_slow_DOC_ptr,idx,soil_DOC_slow)
     id_restart = register_restart_field(soil_restart,fname,'soil_DOC_slow',soil_DOC_slow,compressed=.true., &
                                         longname='Dissolved slow carbon',units='kg/m2')
     allocate(soil_DOC_deadmic(isize,num_l))
     call gather_tile_data(soil_deadmicrobe_DOC_ptr,idx,soil_DOC_deadmic)
     id_restart = register_restart_field(soil_restart,fname,'soil_DOC_deadmic',soil_DOC_deadmic,compressed=.true., &
                                         longname='Dissolved dead microbe carbon',units='kg/m2')



     allocate(fast_DOC_leached(isize))
     call gather_tile_data(soil_fast_DOC_leached_ptr,idx,fast_DOC_leached)
     id_restart = register_restart_field(soil_restart,fname,'fast_DOC_leached',fast_DOC_leached, &
                                         longname='Cumulative fast DOC leached out of the column',units='kg/m2')
     allocate(slow_DOC_leached(isize))
     call gather_tile_data(soil_slow_DOC_leached_ptr,idx,slow_DOC_leached)
     id_restart = register_restart_field(soil_restart,fname,'slow_DOC_leached',slow_DOC_leached, &
                                         longname='Cumulative slow DOC leached out of the column',units='kg/m2')
     allocate(deadmic_DOC_leached(isize))
     call gather_tile_data(soil_deadmic_DOC_leached_ptr,idx,deadmic_DOC_leached)
     id_restart = register_restart_field(soil_restart,fname,'deadmic_DOC_leached',deadmic_DOC_leached, &
                                         longname='Cumulative dead microbe DOC leached out of the column',units='kg/m2')



     allocate(leaf_litter_fast_C(isize,soilMaxCohorts))
     call gather_tile_data(soilc_leafLitter_fast_soil_C_ptr,idx,leaf_litter_fast_C)
     id_restart = register_restart_field(soil_restart,fname,'leaf_litter_fast_C',leaf_litter_fast_C,compressed=.true., &
                                         longname='Leaf litter fast C',units='kg/m2', compressed_axis='C_CC')
     allocate(leaf_litter_slow_C(isize,soilMaxCohorts))
     call gather_tile_data(soilc_leafLitter_slow_soil_C_ptr,idx,leaf_litter_slow_C)
     id_restart = register_restart_field(soil_restart,fname,'leaf_litter_slow_C',leaf_litter_slow_C,compressed=.true., &
                                         longname='Leaf litter slow C',units='kg/m2', compressed_axis='C_CC')
     allocate(leaf_litter_deadMic_C(isize,soilMaxCohorts))
     call gather_tile_data(soilc_leafLitter_deadMicrobeC_ptr,idx,leaf_litter_deadMic_C)
     id_restart = register_restart_field(soil_restart,fname,'leaf_litter_deadMic_C',leaf_litter_deadMic_C,compressed=.true., &
                                         longname='Leaf litter dead microbe C',units='kg/m2', compressed_axis='C_CC')
     allocate(leaf_litter_liveMic_C(isize,soilMaxCohorts))
     call gather_tile_data(soilc_leafLitter_livingMicrobeC_ptr,idx,leaf_litter_liveMic_C)
     id_restart = register_restart_field(soil_restart,fname,'leaf_litter_liveMic_C',leaf_litter_liveMic_C,compressed=.true., &
                                         longname='Leaf litter live microbe C',units='kg/m2', compressed_axis='C_CC')
     allocate(leaf_litter_CO2(isize,soilMaxCohorts))
     call gather_tile_data(soilc_leafLitter_CO2_ptr,idx,leaf_litter_CO2)
     id_restart = register_restart_field(soil_restart,fname,'leaf_litter_CO2',leaf_litter_CO2,compressed=.true., &
                                         longname='Leaf litter CO2 generated',units='kg/m2', compressed_axis='C_CC')
     allocate(leaf_litter_Rtot(isize,soilMaxCohorts))
     call gather_tile_data(soilc_leafLitter_Rtot_ptr,idx,leaf_litter_Rtot)
     id_restart = register_restart_field(soil_restart,fname,'leaf_litter_Rtot',leaf_litter_Rtot,compressed=.true., &
                                         longname='Leaf litter total degradation',units='kg/m2', compressed_axis='C_CC')
     allocate(leaf_litter_originalCohortC(isize,soilMaxCohorts))
     call gather_tile_data(soilc_leafLitter_originalLitterC_ptr,idx,leaf_litter_originalCohortC)
     id_restart = register_restart_field(soil_restart,fname,'leaf_litter_originalCohortC',leaf_litter_originalCohortC,compressed=.true., &
                                         longname='Leaf litter cohort original carbon',units='kg/m2', compressed_axis='C_CC')
     allocate(leaf_litter_fastProtectedC(isize,soilMaxCohorts))
     call gather_tile_data(soilc_leafLitter_fast_protected_C_ptr,idx,leaf_litter_fastProtectedC)
     id_restart = register_restart_field(soil_restart,fname,'leaf_litter_fastProtectedC',leaf_litter_fastProtectedC,compressed=.true., &
                                         longname='Leaf litter fast protected C',units='kg/m2', compressed_axis='C_CC')
     allocate(leaf_litter_slowProtectedC(isize,soilMaxCohorts))
     call gather_tile_data(soilc_leafLitter_slow_protected_C_ptr,idx,leaf_litter_slowProtectedC)
     id_restart = register_restart_field(soil_restart,fname,'leaf_litter_slowProtectedC',leaf_litter_slowProtectedC,compressed=.true., &
                                         longname='Leaf litter slow protected C',units='kg/m2', compressed_axis='C_CC')
     allocate(leaf_litter_deadMicrobeProtectedC(isize,soilMaxCohorts))
     call gather_tile_data(soilc_leafLitter_deadMicrobe_protected_C_ptr,idx,leaf_litter_deadMicrobeProtectedC)
     id_restart = register_restart_field(soil_restart,fname,'leaf_litter_deadMicrobeProtectedC',leaf_litter_deadMicrobeProtectedC,compressed=.true., &
                                         longname='Leaf litter dead microbe protected C',units='kg/m2', compressed_axis='C_CC')
     allocate(leaf_litter_DOC_fast(isize))
     call gather_tile_data(soilc_leafLitter_fast_DOC_ptr,idx,leaf_litter_DOC_fast)
     id_restart = register_restart_field(soil_restart,fname,'leaf_litter_DOC_fast',leaf_litter_DOC_fast, &
                                         longname='Dissolved leaf litter fast carbon',units='kg/m2')
     allocate(leaf_litter_DOC_slow(isize))
     call gather_tile_data(soilc_leafLitter_slow_DOC_ptr,idx,leaf_litter_DOC_slow)
     id_restart = register_restart_field(soil_restart,fname,'leaf_litter_DOC_slow',leaf_litter_DOC_slow, &
                                         longname='Dissolved leaf litter slow carbon',units='kg/m2')
     allocate(leaf_litter_DOC_deadmic(isize))
     call gather_tile_data(soilc_leafLitter_deadmicrobe_DOC_ptr,idx,leaf_litter_DOC_deadmic)
     id_restart = register_restart_field(soil_restart,fname,'leaf_litter_DOC_deadmic',leaf_litter_DOC_deadmic, &
                                         longname='Dissolved leaf litter dead microbe carbon',units='kg/m2')



     allocate(fineWood_litter_fast_C(isize,soilMaxCohorts))
     call gather_tile_data(soilc_fineWoodLitter_fast_soil_C_ptr,idx,fineWood_litter_fast_C)
     id_restart = register_restart_field(soil_restart,fname,'fineWood_litter_fast_C',fineWood_litter_fast_C,compressed=.true., &
                                         longname='Fine wood litter fast C',units='kg/m2', compressed_axis='C_CC')
     allocate(fineWood_litter_slow_C(isize,soilMaxCohorts))
     call gather_tile_data(soilc_fineWoodLitter_slow_soil_C_ptr,idx,fineWood_litter_slow_C)
     id_restart = register_restart_field(soil_restart,fname,'fineWood_litter_slow_C',fineWood_litter_slow_C,compressed=.true., &
                                         longname='Fine wood litter slow C',units='kg/m2', compressed_axis='C_CC')
     allocate(fineWood_litter_deadMic_C(isize,soilMaxCohorts))
     call gather_tile_data(soilc_fineWoodLitter_deadMicrobeC_ptr,idx,fineWood_litter_deadMic_C)
     id_restart = register_restart_field(soil_restart,fname,'fineWood_litter_deadMic_C',fineWood_litter_deadMic_C,compressed=.true., &
                                         longname='Fine wood litter dead microbe C',units='kg/m2', compressed_axis='C_CC')
     allocate(fineWood_litter_liveMic_C(isize,soilMaxCohorts))
     call gather_tile_data(soilc_fineWoodLitter_livingMicrobeC_ptr,idx,fineWood_litter_liveMic_C)
     id_restart = register_restart_field(soil_restart,fname,'fineWood_litter_liveMic_C',fineWood_litter_liveMic_C,compressed=.true., &
                                         longname='Fine wood litter live microbe C',units='kg/m2', compressed_axis='C_CC')
     allocate(fineWood_litter_CO2(isize,soilMaxCohorts))
     call gather_tile_data(soilc_fineWoodLitter_CO2_ptr,idx,fineWood_litter_CO2)
     id_restart = register_restart_field(soil_restart,fname,'fineWood_litter_CO2',fineWood_litter_CO2,compressed=.true., &
                                         longname='Fine wood litter CO2 generated',units='kg/m2', compressed_axis='C_CC')
     allocate(fineWood_litter_Rtot(isize,soilMaxCohorts))
     call gather_tile_data(soilc_fineWoodLitter_Rtot_ptr,idx,fineWood_litter_Rtot)
     id_restart = register_restart_field(soil_restart,fname,'fineWood_litter_Rtot',fineWood_litter_Rtot,compressed=.true., &
                                         longname='Fine wood litter total degradation',units='kg/m2', compressed_axis='C_CC')
     allocate(fineWood_litter_originalCohortC(isize,soilMaxCohorts))
     call gather_tile_data(soilc_fineWoodLitter_originalLitterC_ptr,idx,fineWood_litter_originalCohortC)
     id_restart = register_restart_field(soil_restart,fname,'fineWood_litter_originalCohortC',fineWood_litter_originalCohortC,compressed=.true., &
                                         longname='Fine wood litter cohort original carbon',units='kg/m2', compressed_axis='C_CC')
     allocate(fineWood_litter_fastProtectedC(isize,soilMaxCohorts))
     call gather_tile_data(soilc_fineWoodLitter_fast_protected_C_ptr,idx,fineWood_litter_fastProtectedC)
     id_restart = register_restart_field(soil_restart,fname,'fineWood_litter_fastProtectedC',fineWood_litter_fastProtectedC,compressed=.true., &
                                         longname='Fine wood litter fast protected C',units='kg/m2', compressed_axis='C_CC')
     allocate(fineWood_litter_slowProtectedC(isize,soilMaxCohorts))
     call gather_tile_data(soilc_fineWoodLitter_slow_protected_C_ptr,idx,fineWood_litter_slowProtectedC)
     id_restart = register_restart_field(soil_restart,fname,'fineWood_litter_slowProtectedC',fineWood_litter_slowProtectedC,compressed=.true., &
                                         longname='Fine wood litter slow protected C',units='kg/m2', compressed_axis='C_CC')
     allocate(fineWood_litter_deadMicrobeProtectedC(isize,soilMaxCohorts))
     call gather_tile_data(soilc_fineWoodLitter_deadMicrobe_protected_C_ptr,idx,fineWood_litter_deadMicrobeProtectedC)
     id_restart = register_restart_field(soil_restart,fname,'fineWood_litter_deadMicrobeProtectedC',fineWood_litter_deadMicrobeProtectedC,compressed=.true., &
                                         longname='Fine wood litter dead microbe protected C',units='kg/m2', compressed_axis='C_CC')
     allocate(fineWood_litter_DOC_fast(isize))
     call gather_tile_data(soilc_fineWoodLitter_fast_DOC_ptr,idx,fineWood_litter_DOC_fast)
     id_restart = register_restart_field(soil_restart,fname,'fineWood_litter_DOC_fast',fineWood_litter_DOC_fast, &
                                         longname='Dissolved fine wood litter fast carbon',units='kg/m2')
     allocate(fineWood_litter_DOC_slow(isize))
     call gather_tile_data(soilc_fineWoodLitter_slow_DOC_ptr,idx,fineWood_litter_DOC_slow)
     id_restart = register_restart_field(soil_restart,fname,'fineWood_litter_DOC_slow',fineWood_litter_DOC_slow, &
                                         longname='Dissolved fine wood litter slow carbon',units='kg/m2')
     allocate(fineWood_litter_DOC_deadmic(isize))
     call gather_tile_data(soilc_fineWoodLitter_deadmicrobe_DOC_ptr,idx,fineWood_litter_DOC_deadmic)
     id_restart = register_restart_field(soil_restart,fname,'fineWood_litter_DOC_deadmic',fineWood_litter_DOC_deadmic, &
                                         longname='Dissolved fine wood litter dead microbe carbon',units='kg/m2')



     allocate(coarseWood_litter_fast_C(isize,soilMaxCohorts))
     call gather_tile_data(soilc_coarseWoodLitter_fast_soil_C_ptr,idx,coarseWood_litter_fast_C)
     id_restart = register_restart_field(soil_restart,fname,'coarseWood_litter_fast_C',coarseWood_litter_fast_C,compressed=.true., &
                                         longname='Coarse wood litter fast C',units='kg/m2', compressed_axis='C_CC')
     allocate(coarseWood_litter_slow_C(isize,soilMaxCohorts))
     call gather_tile_data(soilc_coarseWoodLitter_slow_soil_C_ptr,idx,coarseWood_litter_slow_C)
     id_restart = register_restart_field(soil_restart,fname,'coarseWood_litter_slow_C',coarseWood_litter_slow_C,compressed=.true., &
                                         longname='Coarse wood litter slow C',units='kg/m2', compressed_axis='C_CC')
     allocate(coarseWood_litter_deadMic_C(isize,soilMaxCohorts))
     call gather_tile_data(soilc_coarseWoodLitter_deadMicrobeC_ptr,idx,coarseWood_litter_deadMic_C)
     id_restart = register_restart_field(soil_restart,fname,'coarseWood_litter_deadMic_C',coarseWood_litter_deadMic_C,compressed=.true., &
                                         longname='Coarse wood litter dead microbe C',units='kg/m2', compressed_axis='C_CC')
     allocate(coarseWood_litter_liveMic_C(isize,soilMaxCohorts))
     call gather_tile_data(soilc_coarseWoodLitter_livingMicrobeC_ptr,idx,coarseWood_litter_liveMic_C)
     id_restart = register_restart_field(soil_restart,fname,'coarseWood_litter_liveMic_C',coarseWood_litter_liveMic_C,compressed=.true., &
                                         longname='Coarse wood litter live microbe C',units='kg/m2', compressed_axis='C_CC')
     allocate(coarseWood_litter_CO2(isize,soilMaxCohorts))
     call gather_tile_data(soilc_coarseWoodLitter_CO2_ptr,idx,coarseWood_litter_CO2)
     id_restart = register_restart_field(soil_restart,fname,'coarseWood_litter_CO2',coarseWood_litter_CO2,compressed=.true., &
                                         longname='Coarse wood litter CO2 generated',units='kg/m2', compressed_axis='C_CC')
     allocate(coarseWood_litter_Rtot(isize,soilMaxCohorts))
     call gather_tile_data(soilc_coarseWoodLitter_Rtot_ptr,idx,coarseWood_litter_Rtot)
     id_restart = register_restart_field(soil_restart,fname,'coarseWood_litter_Rtot',coarseWood_litter_Rtot,compressed=.true., &
                                         longname='Coarse wood litter total degradation',units='kg/m2', compressed_axis='C_CC')
     allocate(coarseWood_litter_originalCohortC(isize,soilMaxCohorts))
     call gather_tile_data(soilc_coarseWoodLitter_originalLitterC_ptr,idx,coarseWood_litter_originalCohortC)
     id_restart = register_restart_field(soil_restart,fname,'coarseWood_litter_originalCohortC',coarseWood_litter_originalCohortC,compressed=.true., &
                                         longname='Coarse wood litter cohort original carbon',units='kg/m2', compressed_axis='C_CC')
     allocate(coarseWood_litter_fastProtectedC(isize,soilMaxCohorts))
     call gather_tile_data(soilc_coarseWoodLitter_fast_protected_C_ptr,idx,coarseWood_litter_fastProtectedC)
     id_restart = register_restart_field(soil_restart,fname,'coarseWood_litter_fastProtectedC',coarseWood_litter_fastProtectedC,compressed=.true., &
                                         longname='Coarse wood litter fast protected C',units='kg/m2', compressed_axis='C_CC')
     allocate(coarseWood_litter_slowProtectedC(isize,soilMaxCohorts))
     call gather_tile_data(soilc_coarseWoodLitter_slow_protected_C_ptr,idx,coarseWood_litter_slowProtectedC)
     id_restart = register_restart_field(soil_restart,fname,'coarseWood_litter_slowProtectedC',coarseWood_litter_slowProtectedC,compressed=.true., &
                                         longname='Coarse wood litter slow protected C',units='kg/m2', compressed_axis='C_CC')
     allocate(coarseWood_litter_deadMicrobeProtectedC(isize,soilMaxCohorts))
     call gather_tile_data(soilc_coarseWoodLitter_deadMicrobe_protected_C_ptr,idx,coarseWood_litter_deadMicrobeProtectedC)
     id_restart = register_restart_field(soil_restart,fname,'coarseWood_litter_deadMicrobeProtectedC',coarseWood_litter_deadMicrobeProtectedC,compressed=.true., &
                                         longname='Coarse wood litter dead microbe protected C',units='kg/m2', compressed_axis='C_CC')
     allocate(coarseWood_litter_DOC_fast(isize))
     call gather_tile_data(soilc_coarseWoodLitter_fast_DOC_ptr,idx,coarseWood_litter_DOC_fast)
     id_restart = register_restart_field(soil_restart,fname,'coarseWood_litter_DOC_fast',coarseWood_litter_DOC_fast, &
                                         longname='Dissolved coarse wood litter fast carbon',units='kg/m2')
     allocate(coarseWood_litter_DOC_slow(isize))
     call gather_tile_data(soilc_coarseWoodLitter_slow_DOC_ptr,idx,coarseWood_litter_DOC_slow)
     id_restart = register_restart_field(soil_restart,fname,'coarseWood_litter_DOC_slow',coarseWood_litter_DOC_slow, &
                                         longname='Dissolved coarse wood litter slow carbon',units='kg/m2')
     allocate(coarseWood_litter_DOC_deadmic(isize))
     call gather_tile_data(soilc_coarseWoodLitter_deadmicrobe_DOC_ptr,idx,coarseWood_litter_DOC_deadmic)
     id_restart = register_restart_field(soil_restart,fname,'coarseWood_litter_DOC_deadmic',coarseWood_litter_DOC_deadmic, &
                                         longname='Dissolved coarse wood litter dead microbe carbon',units='kg/m2')

     allocate(is_peat(isize,num_l))
     call gather_tile_data(soil_is_peat_ptr,idx,is_peat)
     id_restart = register_restart_field(soil_restart,fname,'is_peat',is_peat,compressed=.true., &
                                         longname='Is layer peat?',units='Boolean')


     if (soil_carbon_option == SOILC_CORPSE_N) then
          allocate(fast_soil_N(isize,num_l,soilMaxCohorts))
          call gather_tile_data(soilc_fast_soil_N_ptr,idx,fast_soil_N)
          id_restart = register_restart_field(soil_restart,fname,'fast_soil_N',fast_soil_N,compressed=.true., &
                                              longname='Fast soil nitrogen',units='kg/m2')
          allocate(slow_soil_N(isize,num_l,soilMaxCohorts))
          call gather_tile_data(soilc_slow_soil_N_ptr,idx,slow_soil_N)
          id_restart = register_restart_field(soil_restart,fname,'slow_soil_N',slow_soil_N,compressed=.true., &
                                              longname='Slow soil nitrogen',units='kg/m2')
          allocate(deadMic_N(isize,num_l,soilMaxCohorts))
          call gather_tile_data(soilc_deadMicrobeN_ptr,idx,deadMic_N)
          id_restart = register_restart_field(soil_restart,fname,'deadMic_N',deadMic_N,compressed=.true., &
                                              longname='Dead microbe nitrogen',units='kg/m2')
         allocate(fastProtectedN(isize,num_l,soilMaxCohorts))
         call gather_tile_data(soilc_fast_protected_N_ptr,idx,fastProtectedN)
         id_restart = register_restart_field(soil_restart,fname,'fastProtectedN',fastProtectedN,compressed=.true., &
                                             longname='Protected fast nitrogen',units='kg/m2')
         allocate(slowProtectedN(isize,num_l,soilMaxCohorts))
         call gather_tile_data(soilc_slow_protected_N_ptr,idx,slowProtectedN)
         id_restart = register_restart_field(soil_restart,fname,'slowProtectedN',slowProtectedN,compressed=.true., &
                                             longname='Protected slow nitrogen',units='kg/m2')
         allocate(deadMicrobeProtectedN(isize,num_l,soilMaxCohorts))
         call gather_tile_data(soilc_deadMicrobe_protected_N_ptr,idx,deadMicrobeProtectedN)
         id_restart = register_restart_field(soil_restart,fname,'deadMicrobeProtectedN',deadMicrobeProtectedN,compressed=.true., &
                                             longname='Protected dead microbe nitrogen',units='kg/m2')
         allocate(originalCohortN(isize,num_l,soilMaxCohorts))
         call gather_tile_data(soilc_originalLitterN_ptr,idx,originalCohortN)
         id_restart = register_restart_field(soil_restart,fname,'originalCohortN',originalCohortN,compressed=.true., &
                                             longname='Cohort original nitrogen',units='kg/m2')

         allocate(liveMic_N(isize,num_l,soilMaxCohorts))
         call gather_tile_data(soilc_livingMicrobeN_ptr,idx,liveMic_N)
         id_restart = register_restart_field(soil_restart,fname,'liveMic_N',liveMic_N,compressed=.true., &
                                             longname='Living microbial nitrogen',units='kg/m2')

         allocate(soil_DON_fast(isize,num_l))
         call gather_tile_data(soil_fast_DON_ptr,idx,soil_DON_fast)
         id_restart = register_restart_field(soil_restart,fname,'soil_DON_fast',soil_DON_fast,compressed=.true., &
                                             longname='Dissolved fast nitrogen',units='kg/m2')
         allocate(soil_DON_slow(isize,num_l))
         call gather_tile_data(soil_slow_DON_ptr,idx,soil_DON_slow)
         id_restart = register_restart_field(soil_restart,fname,'soil_DON_slow',soil_DON_slow,compressed=.true., &
                                             longname='Dissolved slow nitrogen',units='kg/m2')
         allocate(soil_DON_deadmic(isize,num_l))
         call gather_tile_data(soil_deadmicrobe_DON_ptr,idx,soil_DON_deadmic)
         id_restart = register_restart_field(soil_restart,fname,'soil_DON_deadmic',soil_DON_deadmic,compressed=.true., &
                                             longname='Dissolved dead microbe nitrogen',units='kg/m2')
         allocate(soil_NO3(isize,num_l))
         call gather_tile_data(soil_nitrate_ptr,idx,soil_NO3)
         id_restart = register_restart_field(soil_restart,fname,'soil_NO3',soil_NO3,compressed=.true., &
                                             longname='NO3 content',units='kgN/m2')
         allocate(soil_NH4(isize,num_l))
         call gather_tile_data(soil_ammonium_ptr,idx,soil_NH4)
         id_restart = register_restart_field(soil_restart,fname,'soil_NH4',soil_NH4,compressed=.true., &
                                             longname='NH4 content',units='kgN/m2')
         allocate(soil_nitrif(isize,num_l))
         call gather_tile_data(soil_nitrif_ptr,idx,soil_nitrif)
         id_restart = register_restart_field(soil_restart,fname,'soil_nitrif',soil_nitrif,compressed=.true., &
                                             longname='Cumulative nitrification',units='kgN/m2')
         allocate(soil_denitrif(isize,num_l))
         call gather_tile_data(soil_denitrif_ptr,idx,soil_denitrif)
         id_restart = register_restart_field(soil_restart,fname,'soil_denitrif',soil_denitrif,compressed=.true., &
                                             longname='Cumulative denitrification',units='kgN/m2')

         allocate(fast_DON_leached(isize))
         call gather_tile_data(soil_fast_DON_leached_ptr,idx,fast_DON_leached)
         id_restart = register_restart_field(soil_restart,fname,'fast_DON_leached',fast_DON_leached, &
                                             longname='Cumulative fast DON leached out of the column',units='kg/m2')
         allocate(slow_DON_leached(isize))
         call gather_tile_data(soil_slow_DOC_leached_ptr,idx,slow_DON_leached)
         id_restart = register_restart_field(soil_restart,fname,'slow_DON_leached',slow_DON_leached, &
                                             longname='Cumulative slow DON leached out of the column',units='kg/m2')
         allocate(deadmic_DON_leached(isize))
         call gather_tile_data(soil_deadmic_DON_leached_ptr,idx,deadmic_DON_leached)
         id_restart = register_restart_field(soil_restart,fname,'deadmic_DON_leached',deadmic_DON_leached, &
                                             longname='Cumulative dead microbe DON leached out of the column',units='kg/m2')

         allocate(NO3_leached(isize))
         call gather_tile_data(soil_NO3_leached_ptr,idx,NO3_leached)
         id_restart = register_restart_field(soil_restart,fname,'NO3_leached',NO3_leached, &
                                             longname='Cumulative NO3 leached out of the column',units='kgN/m2')

         allocate(NH4_leached(isize))
         call gather_tile_data(soil_NH4_leached_ptr,idx,NH4_leached)
         id_restart = register_restart_field(soil_restart,fname,'NH4_leached',NH4_leached, &
                                             longname='Cumulative NH4 leached out of the column',units='kgN/m2')


         allocate(leaf_litter_fast_N(isize,soilMaxCohorts))
         call gather_tile_data(soilc_leafLitter_fast_soil_N_ptr,idx,leaf_litter_fast_N)
         id_restart = register_restart_field(soil_restart,fname,'leaf_litter_fast_N',leaf_litter_fast_N,compressed=.true., &
                                             longname='Leaf litter fast N',units='kg/m2', compressed_axis='C_CC')
         allocate(leaf_litter_slow_N(isize,soilMaxCohorts))
         call gather_tile_data(soilc_leafLitter_slow_soil_N_ptr,idx,leaf_litter_slow_N)
         id_restart = register_restart_field(soil_restart,fname,'leaf_litter_slow_N',leaf_litter_slow_N,compressed=.true., &
                                             longname='Leaf litter slow N',units='kg/m2', compressed_axis='C_CC')
         allocate(leaf_litter_deadMic_N(isize,soilMaxCohorts))
         call gather_tile_data(soilc_leafLitter_deadMicrobeN_ptr,idx,leaf_litter_deadMic_N)
         id_restart = register_restart_field(soil_restart,fname,'leaf_litter_deadMic_N',leaf_litter_deadMic_N,compressed=.true., &
                                             longname='Leaf litter dead microbe N',units='kg/m2', compressed_axis='C_CC')
         allocate(leaf_litter_liveMic_N(isize,soilMaxCohorts))
         call gather_tile_data(soilc_leafLitter_livingMicrobeN_ptr,idx,leaf_litter_liveMic_N)
         id_restart = register_restart_field(soil_restart,fname,'leaf_litter_liveMic_N',leaf_litter_liveMic_N,compressed=.true., &
                                             longname='Leaf litter live microbe N',units='kg/m2', compressed_axis='C_CC')
         allocate(leaf_litter_originalCohortN(isize,soilMaxCohorts))
         call gather_tile_data(soilc_leafLitter_originalLitterN_ptr,idx,leaf_litter_originalCohortN)
         id_restart = register_restart_field(soil_restart,fname,'leaf_litter_originalCohortN',leaf_litter_originalCohortN,compressed=.true., &
                                             longname='Leaf litter cohort original nitrogen',units='kg/m2', compressed_axis='C_CC')
         allocate(leaf_litter_fastProtectedN(isize,soilMaxCohorts))
         call gather_tile_data(soilc_leafLitter_fast_protected_N_ptr,idx,leaf_litter_fastProtectedN)
         id_restart = register_restart_field(soil_restart,fname,'leaf_litter_fastProtectedN',leaf_litter_fastProtectedN,compressed=.true., &
                                             longname='Leaf litter fast protected N',units='kg/m2', compressed_axis='C_CC')
         allocate(leaf_litter_slowProtectedN(isize,soilMaxCohorts))
         call gather_tile_data(soilc_leafLitter_slow_protected_N_ptr,idx,leaf_litter_slowProtectedN)
         id_restart = register_restart_field(soil_restart,fname,'leaf_litter_slowProtectedN',leaf_litter_slowProtectedN,compressed=.true., &
                                             longname='Leaf litter slow protected N',units='kg/m2', compressed_axis='C_CC')
         allocate(leaf_litter_deadMicrobeProtectedN(isize,soilMaxCohorts))
         call gather_tile_data(soilc_leafLitter_deadMicrobe_protected_N_ptr,idx,leaf_litter_deadMicrobeProtectedN)
         id_restart = register_restart_field(soil_restart,fname,'leaf_litter_deadMicrobeProtectedN',leaf_litter_deadMicrobeProtectedN,compressed=.true., &
                                             longname='Leaf litter dead microbe protected N',units='kg/m2', compressed_axis='C_CC')
         allocate(leaf_litter_DON_fast(isize))
         call gather_tile_data(soilc_leafLitter_fast_DON_ptr,idx,leaf_litter_DON_fast)
         id_restart = register_restart_field(soil_restart,fname,'leaf_litter_DON_fast',leaf_litter_DON_fast, &
                                             longname='Dissolved leaf litter fast nitrogen',units='kg/m2')
         allocate(leaf_litter_DON_slow(isize))
         call gather_tile_data(soilc_leafLitter_slow_DON_ptr,idx,leaf_litter_DON_slow)
         id_restart = register_restart_field(soil_restart,fname,'leaf_litter_DON_slow',leaf_litter_DON_slow, &
                                             longname='Dissolved leaf litter slow nitrogen',units='kg/m2')
         allocate(leaf_litter_DON_deadmic(isize))
         call gather_tile_data(soilc_leafLitter_deadmicrobe_DON_ptr,idx,leaf_litter_DON_deadmic)
         id_restart = register_restart_field(soil_restart,fname,'leaf_litter_DON_deadmic',leaf_litter_DON_deadmic, &
                                             longname='Dissolved leaf litter dead microbe nitrogen',units='kg/m2')
         allocate(leaf_litter_NO3(isize))
         call gather_tile_data(soilc_leafLitter_nitrate_ptr,idx,leaf_litter_NO3)
         id_restart = register_restart_field(soil_restart,fname,'leaf_litter_NO3',leaf_litter_NO3, &
                                             longname='Leaf litter NO3 content',units='kgN/m2')
         allocate(leaf_litter_NH4(isize))
         call gather_tile_data(soilc_leafLitter_ammonium_ptr,idx,leaf_litter_NH4)
         id_restart = register_restart_field(soil_restart,fname,'leaf_litter_NH4',leaf_litter_NH4, &
                                             longname='Leaf litter NH4 content',units='kgN/m2')
         allocate(leaf_litter_nitrif(isize))
         call gather_tile_data(soilc_leafLitter_nitrif_ptr,idx,leaf_litter_nitrif)
         id_restart = register_restart_field(soil_restart,fname,'leaf_litter_nitrif',leaf_litter_nitrif, &
                                             longname='Leaf litter cumulative nitrification',units='kgN/m2')
         allocate(leaf_litter_denitrif(isize))
         call gather_tile_data(soilc_leafLitter_denitrif_ptr,idx,leaf_litter_denitrif)
         id_restart = register_restart_field(soil_restart,fname,'leaf_litter_denitrif',leaf_litter_denitrif, &
                                             longname='Leaf litter cumulative denitrification',units='kgN/m2')

         allocate(fineWood_litter_fast_N(isize,soilMaxCohorts))
         call gather_tile_data(soilc_fineWoodLitter_fast_soil_N_ptr,idx,fineWood_litter_fast_N)
         id_restart = register_restart_field(soil_restart,fname,'fineWood_litter_fast_N',fineWood_litter_fast_N,compressed=.true., &
                                             longname='Fine wood litter fast N',units='kg/m2', compressed_axis='C_CC')
         allocate(fineWood_litter_slow_N(isize,soilMaxCohorts))
         call gather_tile_data(soilc_fineWoodLitter_slow_soil_N_ptr,idx,fineWood_litter_slow_N)
         id_restart = register_restart_field(soil_restart,fname,'fineWood_litter_slow_N',fineWood_litter_slow_N,compressed=.true., &
                                             longname='Fine wood litter slow N',units='kg/m2', compressed_axis='C_CC')
         allocate(fineWood_litter_deadMic_N(isize,soilMaxCohorts))
         call gather_tile_data(soilc_fineWoodLitter_deadMicrobeN_ptr,idx,fineWood_litter_deadMic_N)
         id_restart = register_restart_field(soil_restart,fname,'fineWood_litter_deadMic_N',fineWood_litter_deadMic_N,compressed=.true., &
                                             longname='Fine wood litter dead microbe N',units='kg/m2', compressed_axis='C_CC')
         allocate(fineWood_litter_liveMic_N(isize,soilMaxCohorts))
         call gather_tile_data(soilc_fineWoodLitter_livingMicrobeN_ptr,idx,fineWood_litter_liveMic_N)
         id_restart = register_restart_field(soil_restart,fname,'fineWood_litter_liveMic_N',fineWood_litter_liveMic_N,compressed=.true., &
                                             longname='Fine wood litter live microbe N',units='kg/m2', compressed_axis='C_CC')
         allocate(fineWood_litter_originalCohortN(isize,soilMaxCohorts))
         call gather_tile_data(soilc_fineWoodLitter_originalLitterN_ptr,idx,fineWood_litter_originalCohortN)
         id_restart = register_restart_field(soil_restart,fname,'fineWood_litter_originalCohortN',fineWood_litter_originalCohortN,compressed=.true., &
                                             longname='Fine wood litter cohort original nitrogen',units='kg/m2', compressed_axis='C_CC')
         allocate(fineWood_litter_fastProtectedN(isize,soilMaxCohorts))
         call gather_tile_data(soilc_fineWoodLitter_fast_protected_N_ptr,idx,fineWood_litter_fastProtectedN)
         id_restart = register_restart_field(soil_restart,fname,'fineWood_litter_fastProtectedN',fineWood_litter_fastProtectedN,compressed=.true., &
                                             longname='Fine wood litter fast protected N',units='kg/m2', compressed_axis='C_CC')
         allocate(fineWood_litter_slowProtectedN(isize,soilMaxCohorts))
         call gather_tile_data(soilc_fineWoodLitter_slow_protected_N_ptr,idx,fineWood_litter_slowProtectedN)
         id_restart = register_restart_field(soil_restart,fname,'fineWood_litter_slowProtectedN',fineWood_litter_slowProtectedN,compressed=.true., &
                                             longname='Fine wood litter slow protected N',units='kg/m2', compressed_axis='C_CC')
         allocate(fineWood_litter_deadMicrobeProtectedN(isize,soilMaxCohorts))
         call gather_tile_data(soilc_fineWoodLitter_deadMicrobe_protected_N_ptr,idx,fineWood_litter_deadMicrobeProtectedN)
         id_restart = register_restart_field(soil_restart,fname,'fineWood_litter_deadMicrobeProtectedN',fineWood_litter_deadMicrobeProtectedN,compressed=.true., &
                                             longname='Fine wood litter dead microbe protected N',units='kg/m2', compressed_axis='C_CC')
         allocate(fineWood_litter_DON_fast(isize))
         call gather_tile_data(soilc_fineWoodLitter_fast_DON_ptr,idx,fineWood_litter_DON_fast)
         id_restart = register_restart_field(soil_restart,fname,'fineWood_litter_DON_fast',fineWood_litter_DON_fast, &
                                             longname='Dissolved fine wood litter fast nitrogen',units='kg/m2')
         allocate(fineWood_litter_DON_slow(isize))
         call gather_tile_data(soilc_fineWoodLitter_slow_DON_ptr,idx,fineWood_litter_DON_slow)
         id_restart = register_restart_field(soil_restart,fname,'fineWood_litter_DON_slow',fineWood_litter_DON_slow, &
                                             longname='Dissolved fine wood litter slow nitrogen',units='kg/m2')
         allocate(fineWood_litter_DON_deadmic(isize))
         call gather_tile_data(soilc_fineWoodLitter_deadmicrobe_DON_ptr,idx,fineWood_litter_DON_deadmic)
         id_restart = register_restart_field(soil_restart,fname,'fineWood_litter_DON_deadmic',fineWood_litter_DON_deadmic, &
                                             longname='Dissolved fine wood litter dead microbe nitrogen',units='kg/m2')
         allocate(fineWood_litter_NO3(isize))
         call gather_tile_data(soilc_fineWoodLitter_nitrate_ptr,idx,fineWood_litter_NO3)
         id_restart = register_restart_field(soil_restart,fname,'fineWood_litter_NO3',fineWood_litter_NO3, &
                                             longname='Fine wood litter NO3 content',units='kgN/m2')
         allocate(fineWood_litter_NH4(isize))
         call gather_tile_data(soilc_fineWoodLitter_ammonium_ptr,idx,fineWood_litter_NH4)
         id_restart = register_restart_field(soil_restart,fname,'fineWood_litter_NH4',fineWood_litter_NH4, &
                                             longname='Fine wood litter NH4 content',units='kgN/m2')
         allocate(fineWood_litter_nitrif(isize))
         call gather_tile_data(soilc_fineWoodLitter_nitrif_ptr,idx,fineWood_litter_nitrif)
         id_restart = register_restart_field(soil_restart,fname,'fineWood_litter_nitrif',fineWood_litter_nitrif, &
                                             longname='Fine wood litter cumulative nitrification',units='kgN/m2')
         allocate(fineWood_litter_denitrif(isize))
         call gather_tile_data(soilc_fineWoodLitter_denitrif_ptr,idx,fineWood_litter_denitrif)
         id_restart = register_restart_field(soil_restart,fname,'fineWood_litter_denitrif',fineWood_litter_denitrif, &
                                             longname='Fine wood litter cumulative denitrification',units='kgN/m2')

          allocate(coarseWood_litter_fast_N(isize,soilMaxCohorts))
          call gather_tile_data(soilc_coarseWoodLitter_fast_soil_N_ptr,idx,coarseWood_litter_fast_N)
          id_restart = register_restart_field(soil_restart,fname,'coarseWood_litter_fast_N',coarseWood_litter_fast_N,compressed=.true., &
                                              longname='coarse wood litter fast N',units='kg/m2', compressed_axis='C_CC')
          allocate(coarseWood_litter_slow_N(isize,soilMaxCohorts))
          call gather_tile_data(soilc_coarseWoodLitter_slow_soil_N_ptr,idx,coarseWood_litter_slow_N)
          id_restart = register_restart_field(soil_restart,fname,'coarseWood_litter_slow_N',coarseWood_litter_slow_N,compressed=.true., &
                                              longname='coarse wood litter slow N',units='kg/m2', compressed_axis='C_CC')
          allocate(coarseWood_litter_deadMic_N(isize,soilMaxCohorts))
          call gather_tile_data(soilc_coarseWoodLitter_deadMicrobeN_ptr,idx,coarseWood_litter_deadMic_N)
          id_restart = register_restart_field(soil_restart,fname,'coarseWood_litter_deadMic_N',coarseWood_litter_deadMic_N,compressed=.true., &
                                              longname='coarse wood litter dead microbe N',units='kg/m2', compressed_axis='C_CC')
          allocate(coarseWood_litter_liveMic_N(isize,soilMaxCohorts))
          call gather_tile_data(soilc_coarseWoodLitter_livingMicrobeN_ptr,idx,coarseWood_litter_liveMic_N)
          id_restart = register_restart_field(soil_restart,fname,'coarseWood_litter_liveMic_N',coarseWood_litter_liveMic_N,compressed=.true., &
                                              longname='coarse wood litter live microbe N',units='kg/m2', compressed_axis='C_CC')
          allocate(coarseWood_litter_originalCohortN(isize,soilMaxCohorts))
          call gather_tile_data(soilc_coarseWoodLitter_originalLitterN_ptr,idx,coarseWood_litter_originalCohortN)
          id_restart = register_restart_field(soil_restart,fname,'coarseWood_litter_originalCohortN',coarseWood_litter_originalCohortN,compressed=.true., &
                                              longname='coarse wood litter cohort original nitrogen',units='kg/m2', compressed_axis='C_CC')
          allocate(coarseWood_litter_fastProtectedN(isize,soilMaxCohorts))
          call gather_tile_data(soilc_coarseWoodLitter_fast_protected_N_ptr,idx,coarseWood_litter_fastProtectedN)
          id_restart = register_restart_field(soil_restart,fname,'coarseWood_litter_fastProtectedN',coarseWood_litter_fastProtectedN,compressed=.true., &
                                              longname='coarse wood litter fast protected N',units='kg/m2', compressed_axis='C_CC')
          allocate(coarseWood_litter_slowProtectedN(isize,soilMaxCohorts))
          call gather_tile_data(soilc_coarseWoodLitter_slow_protected_N_ptr,idx,coarseWood_litter_slowProtectedN)
          id_restart = register_restart_field(soil_restart,fname,'coarseWood_litter_slowProtectedN',coarseWood_litter_slowProtectedN,compressed=.true., &
                                              longname='coarse wood litter slow protected N',units='kg/m2', compressed_axis='C_CC')
          allocate(coarseWood_litter_deadMicrobeProtectedN(isize,soilMaxCohorts))
          call gather_tile_data(soilc_coarseWoodLitter_deadMicrobe_protected_N_ptr,idx,coarseWood_litter_deadMicrobeProtectedN)
          id_restart = register_restart_field(soil_restart,fname,'coarseWood_litter_deadMicrobeProtectedN',coarseWood_litter_deadMicrobeProtectedN,compressed=.true., &
                                              longname='coarse wood litter dead microbe protected N',units='kg/m2', compressed_axis='C_CC')
          allocate(coarseWood_litter_DON_fast(isize))
          call gather_tile_data(soilc_coarseWoodLitter_fast_DON_ptr,idx,coarseWood_litter_DON_fast)
          id_restart = register_restart_field(soil_restart,fname,'coarseWood_litter_DON_fast',coarseWood_litter_DON_fast, &
                                              longname='Dissolved coarse wood litter fast nitrogen',units='kg/m2')
          allocate(coarseWood_litter_DON_slow(isize))
          call gather_tile_data(soilc_coarseWoodLitter_slow_DON_ptr,idx,coarseWood_litter_DON_slow)
          id_restart = register_restart_field(soil_restart,fname,'coarseWood_litter_DON_slow',coarseWood_litter_DON_slow, &
                                              longname='Dissolved coarse wood litter slow nitrogen',units='kg/m2')
          allocate(coarseWood_litter_DON_deadmic(isize))
          call gather_tile_data(soilc_coarseWoodLitter_deadmicrobe_DON_ptr,idx,coarseWood_litter_DON_deadmic)
          id_restart = register_restart_field(soil_restart,fname,'coarseWood_litter_DON_deadmic',coarseWood_litter_DON_deadmic, &
                                              longname='Dissolved coarse wood litter dead microbe nitrogen',units='kg/m2')
          allocate(coarseWood_litter_NO3(isize))
          call gather_tile_data(soilc_coarseWoodLitter_nitrate_ptr,idx,coarseWood_litter_NO3)
          id_restart = register_restart_field(soil_restart,fname,'coarseWood_litter_NO3',coarseWood_litter_NO3, &
                                              longname='coarse wood litter NO3 content',units='kgN/m2')
          allocate(coarseWood_litter_NH4(isize))
          call gather_tile_data(soilc_coarseWoodLitter_ammonium_ptr,idx,coarseWood_litter_NH4)
          id_restart = register_restart_field(soil_restart,fname,'coarseWood_litter_NH4',coarseWood_litter_NH4, &
                                              longname='coarse wood litter NH4 content',units='kgN/m2')
          allocate(coarseWood_litter_nitrif(isize))
          call gather_tile_data(soilc_coarseWoodLitter_nitrif_ptr,idx,coarseWood_litter_nitrif)
          id_restart = register_restart_field(soil_restart,fname,'coarseWood_litter_nitrif',coarseWood_litter_nitrif, &
                                              longname='coarse wood litter cumulative nitrification',units='kgN/m2')
          allocate(coarseWood_litter_denitrif(isize))
          call gather_tile_data(soilc_coarseWoodLitter_denitrif_ptr,idx,coarseWood_litter_denitrif)
          id_restart = register_restart_field(soil_restart,fname,'coarseWood_litter_denitrif',coarseWood_litter_denitrif, &
                                              longname='Coarse wood litter cumulative denitrification',units='kgN/m2')


     endif

  end select

  ! save performs io domain aggregation through mpp_io as with regular domain data
  call save_restart(soil_restart)
  call free_restart_type(soil_restart)

  deallocate(temp, wl, ws, groundwater, groundwater_T, uptake_T)
  select case(soil_carbon_option)
  case (SOILC_CENTURY, SOILC_CENTURY_BY_LAYER)
    deallocate(fsc, ssc)
  case (SOILC_CORPSE, SOILC_CORPSE_N)
    deallocate(fast_soil_C, slow_soil_C, deadMic, fastProtectedC, slowProtectedC, deadMicrobeProtectedC, liveMic, CO2, Rtot,     &
      originalCohortC, soil_DOC_fast, soil_DOC_slow, soil_DOC_deadmic, fast_DOC_leached, slow_DOC_leached, deadmic_DOC_leached,  &
      leaf_litter_fast_C, leaf_litter_slow_C, leaf_litter_deadMic_C, leaf_litter_liveMic_C, leaf_litter_CO2,                     &
      leaf_litter_Rtot, leaf_litter_originalCohortC, leaf_litter_fastProtectedC, leaf_litter_slowProtectedC,                     &
      leaf_litter_deadMicrobeProtectedC, leaf_litter_DOC_fast, leaf_litter_DOC_slow, leaf_litter_DOC_deadmic,                    &
      fineWood_litter_fast_C, fineWood_litter_slow_C, fineWood_litter_deadMic_C, fineWood_litter_liveMic_C, fineWood_litter_CO2, &
      fineWood_litter_Rtot, fineWood_litter_originalCohortC, fineWood_litter_fastProtectedC, fineWood_litter_slowProtectedC,     &
      fineWood_litter_deadMicrobeProtectedC, fineWood_litter_DOC_fast, fineWood_litter_DOC_slow, fineWood_litter_DOC_deadmic,    &
      coarseWood_litter_fast_C, coarseWood_litter_slow_C, coarseWood_litter_deadMic_C, coarseWood_litter_liveMic_C,              &
      coarseWood_litter_CO2, coarseWood_litter_Rtot, coarseWood_litter_originalCohortC, coarseWood_litter_fastProtectedC,        &
      coarseWood_litter_slowProtectedC, coarseWood_litter_deadMicrobeProtectedC, coarseWood_litter_DOC_fast,                     &
      coarseWood_litter_DOC_slow, coarseWood_litter_DOC_deadmic, is_peat)

if(soil_carbon_option == SOILC_CORPSE_N) then

      deallocate(fast_soil_N, slow_soil_N, deadMic_N, fastProtectedN, slowProtectedN, deadMicrobeProtectedN, liveMic_N,     &
        originalCohortN, soil_DON_fast, soil_DON_slow, soil_DON_deadmic, fast_DON_leached, slow_DON_leached, deadmic_DON_leached,  &
        NO3_leached, NH4_leached, &
        leaf_litter_fast_N, leaf_litter_slow_N, leaf_litter_deadMic_N, leaf_litter_liveMic_N,                     &
        leaf_litter_originalCohortN, leaf_litter_fastProtectedN, leaf_litter_slowProtectedN,                     &
        leaf_litter_deadMicrobeProtectedN, leaf_litter_DON_fast, leaf_litter_DON_slow, leaf_litter_DON_deadmic,                    &
        fineWood_litter_fast_N, fineWood_litter_slow_N, fineWood_litter_deadMic_N, fineWood_litter_liveMic_N, &
        fineWood_litter_originalCohortN, fineWood_litter_fastProtectedN, fineWood_litter_slowProtectedN,     &
        fineWood_litter_deadMicrobeProtectedN, fineWood_litter_DON_fast, fineWood_litter_DON_slow, fineWood_litter_DON_deadmic,    &
        coarseWood_litter_fast_N, coarseWood_litter_slow_N, coarseWood_litter_deadMic_N, coarseWood_litter_liveMic_N,              &
        coarseWood_litter_originalCohortN, coarseWood_litter_fastProtectedN,        &
        coarseWood_litter_slowProtectedN, coarseWood_litter_deadMicrobeProtectedN, coarseWood_litter_DON_fast,                     &
        coarseWood_litter_DON_slow, coarseWood_litter_DON_deadmic,                                                   &
        soil_NO3, soil_NH4, soil_nitrif, soil_denitrif, leaf_litter_NO3, leaf_litter_NH4, leaf_litter_nitrif, leaf_litter_denitrif, &
        fineWood_litter_NO3, fineWood_litter_NH4, fineWood_litter_nitrif, fineWood_litter_denitrif, &
        coarseWood_litter_NO3, coarseWood_litter_NH4, coarseWood_litter_nitrif, coarseWood_litter_denitrif    )
    endif

  end select

  if (write_soil_carbon_restart) then
    fname = trim(timestamp)//'soil_carbon.res.nc'
    call create_tile_out_file(soil_carbon_restart,idx,fname,soil_tile_exists,tile_dim_length,zaxis_data=zfull(1:num_l))
    allocate(asoil_in(isize,num_l), fsc_in(isize,num_l),  ssc_in(isize,num_l))

    call gather_tile_data(soil_asoil_in_ptr,idx,asoil_in)
    id_restart = register_restart_field(soil_carbon_restart,fname,'asoil_in',asoil_in,compressed=.true., &
                                        longname='aerobic activity modifier',units='unitless')

    call gather_tile_data(soil_fsc_in_ptr,idx,fsc_in)
    id_restart = register_restart_field(soil_carbon_restart,fname,'fsc_in',fsc_in,compressed=.true., &
                                        longname='fast soil carbon input',units='kg C/m2')

    call gather_tile_data(soil_ssc_in_ptr,idx,ssc_in)
    id_restart = register_restart_field(soil_carbon_restart,fname,'ssc_in',ssc_in,compressed=.true., &
                                        longname='slow soil carbon input',units='kg C/m2')

    if (soil_carbon_option == SOILC_CORPSE .or. soil_carbon_option == SOILC_CORPSE_N) then

       allocate(deadmic_in(isize,num_l))
       call gather_tile_data(soil_deadmic_C_in_ptr ,idx,deadmic_in)
       id_restart = register_restart_field(soil_carbon_restart,fname,'deadmic_C_in',deadmic_in,compressed=.true., longname='dead microbe soil carbon input',units='kg C/m2')
       allocate(fast_protected_in(isize,num_l))
       call gather_tile_data(soil_fast_protected_C_in_ptr ,idx,fast_protected_in)
       id_restart = register_restart_field(soil_carbon_restart,fname,'fast_protected_C_in',fast_protected_in,compressed=.true., longname='protected fast soil carbon input',units='kg C/m2')
       allocate(slow_protected_in(isize,num_l))
       call gather_tile_data(soil_slow_protected_C_in_ptr ,idx,slow_protected_in)
       id_restart = register_restart_field(soil_carbon_restart,fname,'slow_protected_C_in',slow_protected_in,compressed=.true., longname='protected slow soil carbon input',units='kg C/m2')
       allocate(deadmic_protected_in(isize,num_l))
       call gather_tile_data(soil_deadmic_protected_C_in_ptr ,idx,deadmic_protected_in)
       id_restart = register_restart_field(soil_carbon_restart,fname,'deadmic_protected_C_in',deadmic_protected_in,compressed=.true., longname='protected dead microbe soil carbon input',units='kg C/m2')
       allocate(fast_turnover_accumulated(isize,num_l))
       call gather_tile_data(soil_fast_C_turnover_accumulated_ptr ,idx,fast_turnover_accumulated)
       id_restart = register_restart_field(soil_carbon_restart,fname,'fast_C_turnover_accumulated',fast_turnover_accumulated,compressed=.true., longname='fast soil carbon turnover',units='year-1')
       allocate(slow_turnover_accumulated(isize,num_l))
       call gather_tile_data(soil_slow_C_turnover_accumulated_ptr ,idx,slow_turnover_accumulated)
       id_restart = register_restart_field(soil_carbon_restart,fname,'slow_C_turnover_accumulated',slow_turnover_accumulated,compressed=.true., longname='slow soil carbon turnover',units='year-1')
       allocate(deadmic_turnover_accumulated(isize,num_l))
       call gather_tile_data(soil_deadmic_C_turnover_accumulated_ptr ,idx,deadmic_turnover_accumulated)
       id_restart = register_restart_field(soil_carbon_restart,fname,'deadmic_C_turnover_accumulated',deadmic_turnover_accumulated,compressed=.true., longname='dead microbe soil carbon turnover',units='year-1')
       allocate(fast_protected_turnover_accumulated(isize,num_l))
       call gather_tile_data(soil_fast_protected_C_turnover_accumulated_ptr ,idx,fast_protected_turnover_accumulated)
       id_restart = register_restart_field(soil_carbon_restart,fname,'fast_protected_C_turnover_accumulated',fast_protected_turnover_accumulated,compressed=.true., longname='fast protected soil carbon turnover',units='year-1')
       allocate(slow_protected_turnover_accumulated(isize,num_l))
       call gather_tile_data(soil_slow_protected_C_turnover_accumulated_ptr ,idx,slow_protected_turnover_accumulated)
       id_restart = register_restart_field(soil_carbon_restart,fname,'slow_protected_C_turnover_accumulated',slow_protected_turnover_accumulated,compressed=.true., longname='slow protected soil carbon turnover',units='year-1')
       allocate(deadmic_protected_turnover_accumulated(isize,num_l))
       call gather_tile_data(soil_deadmic_protected_C_turnover_accumulated_ptr ,idx,deadmic_protected_turnover_accumulated)
       id_restart = register_restart_field(soil_carbon_restart,fname,'deadmic_protected_C_turnover_accumulated',deadmic_protected_turnover_accumulated,compressed=.true., longname='dead microbe protected soil carbon turnover',units='year-1')



       allocate(leaflitter_fast_turnover_accumulated(isize))
       call gather_tile_data(soil_leaflitter_fast_C_turnover_accumulated_ptr ,idx,leaflitter_fast_turnover_accumulated)
       id_restart = register_restart_field(soil_carbon_restart,fname,'leaflitter_fast_C_turnover_accumulated',leaflitter_fast_turnover_accumulated, longname='fast leaf litter carbon turnover',units='year-1')
       allocate(leaflitter_slow_turnover_accumulated(isize))
       call gather_tile_data(soil_leaflitter_slow_C_turnover_accumulated_ptr ,idx,leaflitter_slow_turnover_accumulated)
       id_restart = register_restart_field(soil_carbon_restart,fname,'leaflitter_slow_C_turnover_accumulated',leaflitter_slow_turnover_accumulated, longname='slow leaf litter carbon turnover',units='year-1')
       allocate(leaflitter_deadmic_turnover_accumulated(isize))
       call gather_tile_data(soil_leaflitter_deadmic_C_turnover_accumulated_ptr ,idx,leaflitter_deadmic_turnover_accumulated)
       id_restart = register_restart_field(soil_carbon_restart,fname,'leaflitter_deadmic_C_turnover_accumulated',leaflitter_deadmic_turnover_accumulated, longname='dead microbe leaf litter carbon turnover',units='year-1')
       allocate(leaflitter_fsc_in(isize))
       call gather_tile_data(soil_leaflitter_fsc_in_ptr ,idx,leaflitter_fsc_in)
       id_restart = register_restart_field(soil_carbon_restart,fname,'leaflitter_fsc_in',leaflitter_fsc_in, longname='fast leaf litter carbon input',units='kg C/m2')
       allocate(leaflitter_ssc_in(isize))
       call gather_tile_data(soil_leaflitter_ssc_in_ptr ,idx,leaflitter_ssc_in)
       id_restart = register_restart_field(soil_carbon_restart,fname,'leaflitter_ssc_in',leaflitter_ssc_in, longname='slow leaf litter carbon input',units='kg C/m2')
       allocate(leaflitter_deadmic_in(isize))
       call gather_tile_data(soil_leaflitter_deadmic_C_in_ptr ,idx,leaflitter_deadmic_in)
       id_restart = register_restart_field(soil_carbon_restart,fname,'leaflitter_deadmic_C_in',leaflitter_deadmic_in, longname='dead microbe leaf litter carbon input',units='kg C/m2')



       allocate(finewoodlitter_fast_turnover_accumulated(isize))
       call gather_tile_data(soil_finewoodlitter_fast_C_turnover_accumulated_ptr ,idx,finewoodlitter_fast_turnover_accumulated)
       id_restart = register_restart_field(soil_carbon_restart,fname,'finewoodlitter_fast_C_turnover_accumulated',finewoodlitter_fast_turnover_accumulated, longname='fast fine wood litter carbon turnover',units='year-1')
       allocate(finewoodlitter_slow_turnover_accumulated(isize))
       call gather_tile_data(soil_finewoodlitter_slow_C_turnover_accumulated_ptr ,idx,finewoodlitter_slow_turnover_accumulated)
       id_restart = register_restart_field(soil_carbon_restart,fname,'finewoodlitter_slow_C_turnover_accumulated',finewoodlitter_slow_turnover_accumulated, longname='slow fine wood litter carbon turnover',units='year-1')
       allocate(finewoodlitter_deadmic_turnover_accumulated(isize))
       call gather_tile_data(soil_finewoodlitter_deadmic_C_turnover_accumulated_ptr ,idx,finewoodlitter_deadmic_turnover_accumulated)
       id_restart = register_restart_field(soil_carbon_restart,fname,'finewoodlitter_deadmic_C_turnover_accumulated',finewoodlitter_deadmic_turnover_accumulated, longname='dead microbe fine wood litter carbon turnover',units='year-1')
       allocate(finewoodlitter_fsc_in(isize))
       call gather_tile_data(soil_finewoodlitter_fsc_in_ptr ,idx,finewoodlitter_fsc_in)
       id_restart = register_restart_field(soil_carbon_restart,fname,'finewoodlitter_fsc_in',finewoodlitter_fsc_in, longname='fast fine wood litter carbon input',units='kg C/m2')
       allocate(finewoodlitter_ssc_in(isize))
       call gather_tile_data(soil_finewoodlitter_ssc_in_ptr ,idx,finewoodlitter_ssc_in)
       id_restart = register_restart_field(soil_carbon_restart,fname,'finewoodlitter_ssc_in',finewoodlitter_ssc_in, longname='slow fine wood litter carbon input',units='kg C/m2')
       allocate(finewoodlitter_deadmic_in(isize))
       call gather_tile_data(soil_finewoodlitter_deadmic_C_in_ptr ,idx,finewoodlitter_deadmic_in)
       id_restart = register_restart_field(soil_carbon_restart,fname,'finewoodlitter_deadmic_C_in',finewoodlitter_deadmic_in, longname='dead microbe fine wood litter carbon input',units='kg C/m2')



       allocate(coarsewoodlitter_fast_turnover_accumulated(isize))
       call gather_tile_data(soil_coarsewoodlitter_fast_C_turnover_accumulated_ptr ,idx,coarsewoodlitter_fast_turnover_accumulated)
       id_restart = register_restart_field(soil_carbon_restart,fname,'coarsewoodlitter_fast_C_turnover_accumulated',coarsewoodlitter_fast_turnover_accumulated, longname='fast coarse wood litter carbon turnover',units='year-1')
       allocate(coarsewoodlitter_slow_turnover_accumulated(isize))
       call gather_tile_data(soil_coarsewoodlitter_slow_C_turnover_accumulated_ptr ,idx,coarsewoodlitter_slow_turnover_accumulated)
       id_restart = register_restart_field(soil_carbon_restart,fname,'coarsewoodlitter_slow_C_turnover_accumulated',coarsewoodlitter_slow_turnover_accumulated, longname='slow coarse wood litter carbon turnover',units='year-1')
       allocate(coarsewoodlitter_deadmic_turnover_accumulated(isize))
       call gather_tile_data(soil_coarsewoodlitter_deadmic_C_turnover_accumulated_ptr ,idx,coarsewoodlitter_deadmic_turnover_accumulated)
       id_restart = register_restart_field(soil_carbon_restart,fname,'coarsewoodlitter_deadmic_C_turnover_accumulated',coarsewoodlitter_deadmic_turnover_accumulated, longname='dead microbe coarse wood litter carbon turnover',units='year-1')
       allocate(coarsewoodlitter_fsc_in(isize))
       call gather_tile_data(soil_coarsewoodlitter_fsc_in_ptr ,idx,coarsewoodlitter_fsc_in)
       id_restart = register_restart_field(soil_carbon_restart,fname,'coarsewoodlitter_fsc_in',coarsewoodlitter_fsc_in, longname='fast coarse wood litter carbon input',units='kg C/m2')
       allocate(coarsewoodlitter_ssc_in(isize))
       call gather_tile_data(soil_coarsewoodlitter_ssc_in_ptr ,idx,coarsewoodlitter_ssc_in)
       id_restart = register_restart_field(soil_carbon_restart,fname,'coarsewoodlitter_ssc_in',coarsewoodlitter_ssc_in, longname='slow coarse wood litter carbon input',units='kg C/m2')
       allocate(coarsewoodlitter_deadmic_in(isize))
       call gather_tile_data(soil_coarsewoodlitter_deadmic_C_in_ptr ,idx,coarsewoodlitter_deadmic_in)
       id_restart = register_restart_field(soil_carbon_restart,fname,'coarsewoodlitter_deadmic_C_in',coarsewoodlitter_deadmic_in, longname='dead microbe coarse wood litter carbon input',units='kg C/m2')


if(soil_carbon_option == SOILC_CORPSE_N) then
    allocate(fsn_in(isize,num_l),  ssn_in(isize,num_l))
    call gather_tile_data(soil_fsn_in_ptr,idx,fsn_in)
    id_restart = register_restart_field(soil_carbon_restart,fname,'fsn_in',fsn_in,compressed=.true., &
                                        longname='fast soil nitrogen input',units='kg N/m2')

    call gather_tile_data(soil_ssn_in_ptr,idx,ssn_in)
    id_restart = register_restart_field(soil_carbon_restart,fname,'ssn_in',ssn_in,compressed=.true., &
                                        longname='slow soil nitrogen input',units='kg N/m2')


    allocate(coarsewoodlitter_fast_N_turnover_accumulated(isize))
    call gather_tile_data(soil_coarsewoodlitter_fast_N_turnover_accumulated_ptr ,idx,coarsewoodlitter_fast_N_turnover_accumulated)
    id_restart = register_restart_field(soil_carbon_restart,fname,'coarsewoodlitter_fast_N_turnover_accumulated',coarsewoodlitter_fast_N_turnover_accumulated, longname='fast coarse wood litter nitrogen turnover',units='year-1')
    allocate(coarsewoodlitter_slow_N_turnover_accumulated(isize))
    call gather_tile_data(soil_coarsewoodlitter_slow_N_turnover_accumulated_ptr ,idx,coarsewoodlitter_slow_N_turnover_accumulated)
    id_restart = register_restart_field(soil_carbon_restart,fname,'coarsewoodlitter_slow_N_turnover_accumulated',coarsewoodlitter_slow_N_turnover_accumulated, longname='slow coarse wood litter nitrogen turnover',units='year-1')
    allocate(coarsewoodlitter_deadmic_N_turnover_accumulated(isize))
    call gather_tile_data(soil_coarsewoodlitter_deadmic_N_turnover_accumulated_ptr ,idx,coarsewoodlitter_deadmic_N_turnover_accumulated)
    id_restart = register_restart_field(soil_carbon_restart,fname,'coarsewoodlitter_deadmic_N_turnover_accumulated',coarsewoodlitter_deadmic_N_turnover_accumulated, longname='dead microbe coarse wood litter nitrogen turnover',units='year-1')
    allocate(coarsewoodlitter_fsn_in(isize))
    call gather_tile_data(soil_coarsewoodlitter_fsn_in_ptr ,idx,coarsewoodlitter_fsn_in)
    id_restart = register_restart_field(soil_carbon_restart,fname,'coarsewoodlitter_fsn_in',coarsewoodlitter_fsn_in, longname='fast coarse wood litter nitrogen input',units='kg N/m2')
    allocate(coarsewoodlitter_ssn_in(isize))
    call gather_tile_data(soil_coarsewoodlitter_ssn_in_ptr ,idx,coarsewoodlitter_ssn_in)
    id_restart = register_restart_field(soil_carbon_restart,fname,'coarsewoodlitter_ssn_in',coarsewoodlitter_ssn_in, longname='slow coarse wood litter nitrogen input',units='kg N/m2')
    allocate(coarsewoodlitter_deadmic_N_in(isize))
    call gather_tile_data(soil_coarsewoodlitter_deadmic_N_in_ptr ,idx,coarsewoodlitter_deadmic_N_in)
    id_restart = register_restart_field(soil_carbon_restart,fname,'coarsewoodlitter_deadmic_N_in',coarsewoodlitter_deadmic_N_in, longname='dead microbe coarse wood litter nitrogen input',units='kg N/m2')

    allocate(finewoodlitter_fast_N_turnover_accumulated(isize))
    call gather_tile_data(soil_finewoodlitter_fast_N_turnover_accumulated_ptr ,idx,finewoodlitter_fast_N_turnover_accumulated)
    id_restart = register_restart_field(soil_carbon_restart,fname,'finewoodlitter_fast_N_turnover_accumulated',finewoodlitter_fast_N_turnover_accumulated, longname='fast fine wood litter carbon turnover',units='year-1')
    allocate(finewoodlitter_slow_N_turnover_accumulated(isize))
    call gather_tile_data(soil_finewoodlitter_slow_N_turnover_accumulated_ptr ,idx,finewoodlitter_slow_N_turnover_accumulated)
    id_restart = register_restart_field(soil_carbon_restart,fname,'finewoodlitter_slow_N_turnover_accumulated',finewoodlitter_slow_N_turnover_accumulated, longname='slow fine wood litter carbon turnover',units='year-1')
    allocate(finewoodlitter_deadmic_N_turnover_accumulated(isize))
    call gather_tile_data(soil_finewoodlitter_deadmic_N_turnover_accumulated_ptr ,idx,finewoodlitter_deadmic_N_turnover_accumulated)
    id_restart = register_restart_field(soil_carbon_restart,fname,'finewoodlitter_deadmic_N_turnover_accumulated',finewoodlitter_deadmic_N_turnover_accumulated, longname='dead microbe fine wood litter carbon turnover',units='year-1')
    allocate(finewoodlitter_fsn_in(isize))
    call gather_tile_data(soil_finewoodlitter_fsn_in_ptr ,idx,finewoodlitter_fsn_in)
    id_restart = register_restart_field(soil_carbon_restart,fname,'finewoodlitter_fsn_in',finewoodlitter_fsn_in, longname='fast fine wood litter nitrogen input',units='kg N/m2')
    allocate(finewoodlitter_ssn_in(isize))
    call gather_tile_data(soil_finewoodlitter_ssn_in_ptr ,idx,finewoodlitter_ssn_in)
    id_restart = register_restart_field(soil_carbon_restart,fname,'finewoodlitter_ssn_in',finewoodlitter_ssn_in, longname='slow fine wood litter nitrogen input',units='kg N/m2')
    allocate(finewoodlitter_deadmic_N_in(isize))
    call gather_tile_data(soil_finewoodlitter_deadmic_N_in_ptr ,idx,finewoodlitter_deadmic_N_in)
    id_restart = register_restart_field(soil_carbon_restart,fname,'finewoodlitter_deadmic_N_in',finewoodlitter_deadmic_N_in, longname='dead microbe fine wood litter nitrogen input',units='kg N/m2')

    allocate(leaflitter_fast_N_turnover_accumulated(isize))
    call gather_tile_data(soil_leaflitter_fast_N_turnover_accumulated_ptr ,idx,leaflitter_fast_N_turnover_accumulated)
    id_restart = register_restart_field(soil_carbon_restart,fname,'leaflitter_fast_N_turnover_accumulated',leaflitter_fast_N_turnover_accumulated, longname='fast leaf litter nitrogen turnover',units='year-1')
    allocate(leaflitter_slow_N_turnover_accumulated(isize))
    call gather_tile_data(soil_leaflitter_slow_N_turnover_accumulated_ptr ,idx,leaflitter_slow_N_turnover_accumulated)
    id_restart = register_restart_field(soil_carbon_restart,fname,'leaflitter_slow_N_turnover_accumulated',leaflitter_slow_N_turnover_accumulated, longname='slow leaf litter nitrogen turnover',units='year-1')
    allocate(leaflitter_deadmic_N_turnover_accumulated(isize))
    call gather_tile_data(soil_leaflitter_deadmic_N_turnover_accumulated_ptr ,idx,leaflitter_deadmic_N_turnover_accumulated)
    id_restart = register_restart_field(soil_carbon_restart,fname,'leaflitter_deadmic_N_turnover_accumulated',leaflitter_deadmic_N_turnover_accumulated, longname='dead microbe leaf litter nitrogen turnover',units='year-1')
    allocate(leaflitter_fsn_in(isize))
    call gather_tile_data(soil_leaflitter_fsn_in_ptr ,idx,leaflitter_fsn_in)
    id_restart = register_restart_field(soil_carbon_restart,fname,'leaflitter_fsn_in',leaflitter_fsn_in, longname='fast leaf litter nitrogen input',units='kg N/m2')
    allocate(leaflitter_ssn_in(isize))
    call gather_tile_data(soil_leaflitter_ssn_in_ptr ,idx,leaflitter_ssn_in)
    id_restart = register_restart_field(soil_carbon_restart,fname,'leaflitter_ssn_in',leaflitter_ssn_in, longname='slow leaf litter nitrogen input',units='kg N/m2')
    allocate(leaflitter_deadmic_N_in(isize))
    call gather_tile_data(soil_leaflitter_deadmic_N_in_ptr ,idx,leaflitter_deadmic_N_in)
    id_restart = register_restart_field(soil_carbon_restart,fname,'leaflitter_deadmic_N_in',leaflitter_deadmic_N_in, longname='dead microbe leaf litter nitrogen input',units='kg N/m2')

    allocate(deadmic_N_in(isize,num_l))
    call gather_tile_data(soil_deadmic_N_in_ptr ,idx,deadmic_N_in)
    id_restart = register_restart_field(soil_carbon_restart,fname,'deadmic_N_in',deadmic_N_in,compressed=.true., longname='dead microbe soil nitrogen input',units='kg N/m2')
    allocate(fast_protected_N_in(isize,num_l))
    call gather_tile_data(soil_fast_protected_N_in_ptr ,idx,fast_protected_N_in)
    id_restart = register_restart_field(soil_carbon_restart,fname,'fast_protected_N_in',fast_protected_N_in,compressed=.true., longname='protected fast soil nitrogen input',units='kg N/m2')
    allocate(slow_protected_N_in(isize,num_l))
    call gather_tile_data(soil_slow_protected_N_in_ptr ,idx,slow_protected_N_in)
    id_restart = register_restart_field(soil_carbon_restart,fname,'slow_protected_N_in',slow_protected_N_in,compressed=.true., longname='protected slow soil nitrogen input',units='kg N/m2')
    allocate(deadmic_protected_N_in(isize,num_l))
    call gather_tile_data(soil_deadmic_protected_N_in_ptr ,idx,deadmic_protected_N_in)
    id_restart = register_restart_field(soil_carbon_restart,fname,'deadmic_protected_N_in',deadmic_protected_N_in,compressed=.true., longname='protected dead microbe soil nitrogen input',units='kg N/m2')
    allocate(fast_N_turnover_accumulated(isize,num_l))
    call gather_tile_data(soil_fast_N_turnover_accumulated_ptr ,idx,fast_N_turnover_accumulated)
    id_restart = register_restart_field(soil_carbon_restart,fname,'fast_N_turnover_accumulated',fast_N_turnover_accumulated,compressed=.true., longname='fast soil nitrogen turnover',units='year-1')
    allocate(slow_turnover_accumulated(isize,num_l))
    call gather_tile_data(soil_slow_N_turnover_accumulated_ptr ,idx,slow_N_turnover_accumulated)
    id_restart = register_restart_field(soil_carbon_restart,fname,'slow_N_turnover_accumulated',slow_N_turnover_accumulated,compressed=.true., longname='slow soil nitrogen turnover',units='year-1')
    allocate(deadmic_N_turnover_accumulated(isize,num_l))
    call gather_tile_data(soil_deadmic_N_turnover_accumulated_ptr ,idx,deadmic_N_turnover_accumulated)
    id_restart = register_restart_field(soil_carbon_restart,fname,'deadmic_N_turnover_accumulated',deadmic_N_turnover_accumulated,compressed=.true., longname='dead microbe soil nitrogen turnover',units='year-1')
    allocate(fast_protected_N_turnover_accumulated(isize,num_l))
    call gather_tile_data(soil_fast_protected_N_turnover_accumulated_ptr ,idx,fast_protected_N_turnover_accumulated)
    id_restart = register_restart_field(soil_carbon_restart,fname,'fast_protected_N_turnover_accumulated',fast_protected_N_turnover_accumulated,compressed=.true., longname='fast protected soil nitrogen turnover',units='year-1')
    allocate(slow_protected_turnover_accumulated(isize,num_l))
    call gather_tile_data(soil_slow_protected_N_turnover_accumulated_ptr ,idx,slow_protected_N_turnover_accumulated)
    id_restart = register_restart_field(soil_carbon_restart,fname,'slow_protected_N_turnover_accumulated',slow_protected_N_turnover_accumulated,compressed=.true., longname='slow protected soil nitrogen turnover',units='year-1')
    allocate(deadmic_protected_N_turnover_accumulated(isize,num_l))
    call gather_tile_data(soil_deadmic_protected_N_turnover_accumulated_ptr ,idx,deadmic_protected_N_turnover_accumulated)
    id_restart = register_restart_field(soil_carbon_restart,fname,'deadmic_protected_N_turnover_accumulated',deadmic_protected_N_turnover_accumulated,compressed=.true., longname='dead microbe protected soil nitrogen turnover',units='year-1')


endif

    endif

  call save_restart(soil_carbon_restart)
  call free_restart_type(soil_carbon_restart)
  deallocate(idx, asoil_in, fsc_in, ssc_in)
  if (soil_carbon_option == SOILC_CORPSE .or. soil_carbon_option==SOILC_CORPSE_N) then
    deallocate(deadmic_in,fast_protected_in,slow_protected_in,deadmic_protected_in,fast_turnover_accumulated,                   &
       slow_turnover_accumulated,deadmic_turnover_accumulated,fast_protected_turnover_accumulated,                              &
       slow_protected_turnover_accumulated,deadmic_protected_turnover_accumulated,leaflitter_fast_turnover_accumulated,         &
       leaflitter_slow_turnover_accumulated,leaflitter_deadmic_turnover_accumulated,leaflitter_fsc_in,leaflitter_ssc_in,        &
       leaflitter_deadmic_in,finewoodlitter_fast_turnover_accumulated,finewoodlitter_slow_turnover_accumulated,                 &
       finewoodlitter_deadmic_turnover_accumulated,finewoodlitter_fsc_in,finewoodlitter_ssc_in,finewoodlitter_deadmic_in,       &
       coarsewoodlitter_fast_turnover_accumulated,coarsewoodlitter_slow_turnover_accumulated,                                   &
       coarsewoodlitter_deadmic_turnover_accumulated,coarsewoodlitter_fsc_in,coarsewoodlitter_ssc_in,coarsewoodlitter_deadmic_in)

       deallocate(deadmic_N_in,fast_protected_N_in,slow_protected_N_in,deadmic_protected_N_in,fast_N_turnover_accumulated,                   &
          slow_N_turnover_accumulated,deadmic_N_turnover_accumulated,fast_protected_N_turnover_accumulated,fsn_in,ssn_in ,                &
          slow_protected_N_turnover_accumulated,deadmic_protected_N_turnover_accumulated,leaflitter_fast_N_turnover_accumulated,         &
          leaflitter_slow_N_turnover_accumulated,leaflitter_deadmic_N_turnover_accumulated,leaflitter_fsn_in,leaflitter_ssn_in,        &
          leaflitter_deadmic_N_in,finewoodlitter_fast_N_turnover_accumulated,finewoodlitter_slow_N_turnover_accumulated,                 &
          finewoodlitter_deadmic_N_turnover_accumulated,finewoodlitter_fsn_in,finewoodlitter_ssn_in,finewoodlitter_deadmic_N_in,       &
          coarsewoodlitter_fast_N_turnover_accumulated,coarsewoodlitter_slow_N_turnover_accumulated,                                   &
          coarsewoodlitter_deadmic_N_turnover_accumulated,coarsewoodlitter_fsn_in,coarsewoodlitter_ssn_in,coarsewoodlitter_deadmic_N_in)
  endif
  endif

end subroutine save_soil_restart_new


! ============================================================================
subroutine soil_get_sfc_temp ( soil, soil_T )
  type(soil_tile_type), intent(in) :: soil
  real, intent(out) :: soil_T

  soil_T= soil%T(1)
end subroutine soil_get_sfc_temp


! ============================================================================
! compute soil radiative properties
subroutine soil_radiation ( soil, cosz, &
     soil_refl_dir, soil_refl_dif, soil_refl_lw, soil_emis )
  type(soil_tile_type), intent(in) :: soil
  real, intent(in)  :: cosz
  real, intent(out) :: soil_refl_dir(NBANDS), soil_refl_dif(NBANDS), soil_refl_lw, soil_emis

  call soil_data_radiation ( soil, cosz, use_brdf, soil_refl_dir, soil_refl_dif, soil_emis )
  soil_refl_lw = 1 - soil_emis
  call check_var_range(soil_refl_dir(BAND_VIS), 0.0, 1.0, 'soil_radiation', 'soil_refl_dir(VIS)', FATAL)
  call check_var_range(soil_refl_dir(BAND_NIR), 0.0, 1.0, 'soil_radiation', 'soil_refl_dir(NIR)', FATAL)
  call check_var_range(soil_refl_dif(BAND_VIS), 0.0, 1.0, 'soil_radiation', 'soil_refl_dif(VIS)', FATAL)
  call check_var_range(soil_refl_dif(BAND_NIR), 0.0, 1.0, 'soil_radiation', 'soil_refl_dif(NIR)', FATAL)
end subroutine soil_radiation


! ============================================================================
! compute soil roughness
subroutine soil_diffusion ( soil, soil_z0s, soil_z0m )
  type(soil_tile_type), intent(in) :: soil
  real, intent(out) :: soil_z0s, soil_z0m

  call soil_data_diffusion ( soil, soil_z0s, soil_z0m )
end subroutine soil_diffusion


! ============================================================================
! compute beta function
! after Manabe (1969), but distributed vertically.
subroutine soil_data_beta ( soil, vegn, diag, soil_beta, soil_water_supply, &
                            soil_uptake_T, soil_rh, soil_rh_psi )
  type(soil_tile_type), intent(inout) :: soil
  type(vegn_tile_type), intent(in)    :: vegn
  type(diag_buff_type), intent(inout) :: diag
  real, intent(out) :: soil_beta
  real, intent(out) :: soil_water_supply ! max rate of water supply to roots, kg/(m2 s)
  real, intent(out) :: soil_uptake_T ! an estimate of temperature of the water
             ! taken up by transpiration. In case of 'linear' uptake it is an exact
             ! value; in case of 'darcy*' treatments the actual uptake profile
             ! is calculated only in step 2, so the value returned is an estimate
  real, intent(out) :: soil_rh
  real, intent(out) :: soil_rh_psi

  ! ---- local vars
  integer :: l, iter
  real, dimension(num_l) :: &
       uptake_frac_max, & ! root distribution
       vegn_uptake_term, &
       vlc, vsc, & ! volumetric fractions of water and ice in the layer
       VRL, & ! vertical distribution of volumetric root length, m/m3
       u, du  ! uptake and its derivative (the latter is not used)
  real :: psi_for_rh
  real :: K_r, r_r ! root properties
  real :: psi_crown_min, grav_head, plant_height, xylem_resist, xylem_area_frac, sws, dum4
  real :: DsapDpsi, psi_left, psi_right, psi_mid, f_left, f_right, f_last, f_mid
  real :: z  !  soil depth

  call vegn_uptake_profile (vegn, dz(1:num_l), uptake_frac_max, vegn_uptake_term )

  vlc=0;vsc=0
  do l = 1, num_l
    vlc(l) = max(0., soil%wl(l) / (dens_h2o*dz(l)))
    vsc(l) = max(0., soil%ws(l) / (dens_h2o*dz(l)))
  enddo

  soil%uptake_frac = 0
  do l = 1, num_l
     soil%uptake_frac(l) = uptake_frac_max(l) &
          * max(0.0, min(1.0,(vlc(l)-soil%w_wilt(l))/&
               (0.75*(soil%w_fc(l)-soil%w_wilt(l)))))
  enddo
  soil_beta = sum(soil%uptake_frac)
  do l = 1, num_l
     if (soil_beta /= 0) then
          soil%uptake_frac(l) = soil%uptake_frac(l) / soil_beta
     else
          soil%uptake_frac(l) = uptake_frac_max(l)
     endif
  enddo
  if (lm2) soil%uptake_frac = uptake_frac_max

  ! calculate relative humidity at soil surface
  call soil_data_psi_for_rh ( soil, vlc, vsc, soil%psi, psi_for_rh )
  soil_rh = exp(psi_for_rh*g_RT)
  soil_rh_psi = g_RT*soil_rh

  ! calculate total water supply
  select case (uptake_option)
  case(UPTAKE_LINEAR)
     soil_water_supply = 0
     z = 0
     do l = 1, num_l
        soil_water_supply = soil_water_supply + &
          vegn_uptake_term(l)*max(0.0,soil%wl(l)/dz(l)-soil%w_wilt(l)*dens_h2o)
        z = z + dz(l)
     enddo
     soil_water_supply = z * soil_water_supply
     soil_water_supply = soil_water_supply/delta_time
     soil_uptake_T = sum(soil%uptake_frac*soil%T)
  case(UPTAKE_DARCY2D, UPTAKE_DARCY2D_LIN)
     call vegn_hydraulic_properties (vegn, dz(1:num_l), always_use_bsw, &
                                     VRL, K_r, r_r, &
                                     plant_height, xylem_area_frac, xylem_resist, &
                                     dum4 )
     grav_head = 0
     iter = 0
     psi_mid = 0
     if (use_tall_plants) grav_head = plant_height
     if (xylem_area_frac.gt.0. .and. xylem_resist.gt.0. .and. plant_height.gt.0. ) then
        DsapDpsi = DENS_H2O*xylem_area_frac/(plant_height*xylem_resist)
        psi_left = psi_wilt + grav_head
        call darcy2d_uptake ( soil, psi_left, VRL, &
               K_r, r_r, uptake_option, uptake_oneway, uptake_from_sat, u, du )
        f_left = max(0.0,sum(u)) - (psi_left-psi_wilt-grav_head)*DsapDpsi
        psi_right = 0.
        call darcy2d_uptake ( soil, psi_right, VRL, &
               K_r, r_r, uptake_option, uptake_oneway, uptake_from_sat, u, du )
        f_right = max(0.0,sum(u)) - (psi_right-psi_wilt-grav_head)*DsapDpsi
        f_last = f_left
        do iter = 1, max_iter_trans
           if (is_watch_point()) then
               write(*,*)'##### sws iteration iter=',iter
               __DEBUG5__(psi_left,f_left,psi_right,f_right,f_last)
           endif
           psi_mid = 0.5*(psi_left+psi_right)
           call darcy2d_uptake ( soil, psi_mid, VRL, &
                  K_r, r_r, uptake_option, uptake_oneway, uptake_from_sat, u, du )
           f_mid = max(0.0,sum(u)) - (psi_mid-psi_wilt-grav_head)*DsapDpsi
           if (abs(f_mid-f_last).lt.eps_trans) exit
           f_last = f_mid
           if (f_mid.lt.0.) then
               psi_right = psi_mid
               f_right = f_mid
           else
               psi_left = psi_mid
               f_left = f_mid
           endif
        enddo
     else
        call darcy2d_uptake ( soil, psi_wilt+grav_head, VRL, &
              K_r, r_r, uptake_option, uptake_oneway, uptake_from_sat, u, du )
     endif
     soil_water_supply = max(0.0,sum(u))
     soil_uptake_T = soil%uptake_T
     call send_tile_data(id_sws_n_iter, real(iter), diag)
     call send_tile_data(id_psi_x0_sws, psi_mid, diag)
  end select
end subroutine soil_data_beta


! ============================================================================
! update soil properties explicitly for time step.
! MAY WISH TO INTRODUCE 'UNIVERSAL' SENSITIVITIES FOR SIMPLICITY.
! T-DEPENDENCE OF HYDRAULIC PROPERTIES COULD BE DONE LESS FREQUENTLY.
! integrate soil-heat conduction equation upward from bottom of soil
! to surface, delivering linearization of surface ground heat flux.
subroutine soil_step_1 ( soil, vegn, diag, &
                         soil_T, soil_uptake_T, soil_beta, soil_water_supply, &
                         soil_E_min, soil_E_max, &
                         soil_rh, soil_rh_psi, soil_liq, soil_ice, soil_subl, soil_tf, &
                         soil_G0, soil_DGDT )
  type(soil_tile_type), intent(inout) :: soil
  type(vegn_tile_type), intent(in)    :: vegn
  type(diag_buff_type), intent(inout) :: diag
  real, intent(out) :: &
       soil_T, &    ! temperature of the upper layer of the soil, degK
       soil_uptake_T, & ! estimate of the temperature of the water taken up by transpiration
       soil_beta, &
       soil_water_supply, & ! supply of water to vegetation per unit total active root biomass, kg/m2
       soil_E_min, &
       soil_E_max, &
       soil_rh,   & ! soil surface relative humidity
       soil_rh_psi,& ! derivative of soil_rh w.r.t. soil surface matric head
       soil_liq,  & ! amount of liquid water available for implicit freeze (=0)
       soil_ice,  & ! amount of ice available for implicit melt (=0)
       soil_subl, & ! part of sublimation in water vapor flux, dimensionless [0,1]
       soil_tf,   & ! soil freezing temperature, degK
       soil_G0, soil_DGDT ! linearization of ground heat flux
  ! ---- local vars
  real :: bbb, denom, dt_e
  real, dimension(num_l) :: aaa, ccc, thermal_cond, heat_capacity, vlc, vsc
  integer :: l

  if(is_watch_point()) then
     write(*,*) 'soil%tag', soil%tag
     write(*,*) 'soil%pars%k_sat_ref', soil%pars%k_sat_ref
     write(*,*) 'soil%pars%psi_sat_ref', soil%pars%psi_sat_ref
     write(*,*) 'soil%pars%chb', soil%pars%chb
     write(*,*) 'soil%pars%w_sa', soil%pars%vwc_sat
  endif
! ----------------------------------------------------------------------------
! in preparation for implicit energy balance, determine various measures
! of water availability, so that vapor fluxes will not exceed mass limits
! ----------------------------------------------------------------------------

  soil_T = soil%T(1)
  call soil_data_beta ( soil, vegn, diag, soil_beta, soil_water_supply, soil_uptake_T, &
                        soil_rh, soil_rh_psi )

  do l = 1, num_l
     vlc(l) = max(0.0, soil%wl(l) / (dens_h2o * dz(l)))
     vsc(l) = max(0.0, soil%ws(l) / (dens_h2o * dz(l)))
  enddo
  call soil_data_thermodynamics ( soil, vlc, vsc,  &
                                  soil_E_max, thermal_cond )
  if (.not.use_E_max) soil_E_max =  HUGE(soil_E_max)
  soil_E_min = Eg_min

  do l = 1, num_l
     heat_capacity(l) = soil%heat_capacity_dry(l) *dz(l) &
          + clw*soil%wl(l) + csw*soil%ws(l)
  enddo

  soil_liq  = max(soil%wl(1), 0.)
  soil_ice  = max(soil%ws(1), 0.)
  if (soil_liq + soil_ice > 0) then
     soil_subl = soil_ice / (soil_liq + soil_ice)
  else
     soil_subl = 0
  endif
  soil_liq = 0
  soil_ice = 0

  if(num_l > 1) then
     do l = 1, num_l-1
        dt_e = 2 / ( dz(l+1)/thermal_cond(l+1) &
                     + dz(l)/thermal_cond(l)   )
        aaa(l+1) = - dt_e * delta_time / heat_capacity(l+1)
        ccc(l)   = - dt_e * delta_time / heat_capacity(l)
     enddo

     bbb = 1.0 - aaa(num_l)
     denom = bbb
     dt_e = aaa(num_l)*(soil%T(num_l) - soil%T(num_l-1)) &
               + soil%geothermal_heat_flux * delta_time / heat_capacity(num_l)
     soil%e(num_l-1) = -aaa(num_l)/denom
     soil%f(num_l-1) = dt_e/denom
     do l = num_l-1, 2, -1
        bbb = 1.0 - aaa(l) - ccc(l)
        denom = bbb + ccc(l)*soil%e(l)
        dt_e = - ( ccc(l)*(soil%T(l+1) - soil%T(l)  ) &
                  -aaa(l)*(soil%T(l)   - soil%T(l-1)) )
        soil%e(l-1) = -aaa(l)/denom
        soil%f(l-1) = (dt_e - ccc(l)*soil%f(l))/denom
     end do
     denom = delta_time/(heat_capacity(1) )
     soil_G0   = ccc(1)*(soil%T(2)- soil%T(1) + soil%f(1)) / denom
     soil_DGDT = (1 - ccc(1)*(1-soil%e(1))) / denom
  else  ! one-level case
     denom = delta_time/heat_capacity(1)
     soil_G0    = 0.
     soil_DGDT  = 1. / denom
  end if

  ! set soil freezing temperature
  soil_tf = soil%pars%tfreeze

  if(is_watch_point()) then
     write(*,*) '#### soil_step_1 checkpoint 1 ####'
     write(*,*) 'mask    ', .true.
     write(*,*) 'T       ', soil_T
     write(*,*) 'uptake_T', soil_uptake_T
     write(*,*) 'beta    ', soil_beta
     write(*,*) 'E_max   ', soil_E_max
     write(*,*) 'rh      ', soil_rh
     write(*,*) 'liq     ', soil_liq
     write(*,*) 'ice     ', soil_ice
     write(*,*) 'subl    ', soil_subl
     write(*,*) 'G0      ', soil_G0
     write(*,*) 'DGDT    ', soil_DGDT
     __DEBUG1__(soil_water_supply)
  endif

  call send_tile_data(id_thermal_cond, thermal_cond, diag)

end subroutine soil_step_1


! ============================================================================
! apply boundary flows to soil water and move soil water vertically.
  subroutine soil_step_2 ( soil, vegn, diag, soil_subl, snow_lprec, snow_hlprec,  &
                           vegn_uptk, &
                           subs_DT, subs_M_imp, subs_evap, &
                           use_tfreeze_in_grnd_latent, &
                           soil_levap, soil_fevap, soil_melt, &
                           soil_lrunf, soil_hlrunf, soil_Ttop, soil_Ctop, &
                           soil_frunf, soil_hfrunf, soil_doc_runf, soil_don_runf, soil_NO3_runf, soil_NH4_runf)
  type(soil_tile_type), intent(inout) :: soil
  type(vegn_tile_type), intent(in)    :: vegn
  type(diag_buff_type), intent(inout) :: diag
  real, intent(in) :: & ! ZMS assign tentative annotations below with "??"
       soil_subl     ! ?? solution for soil surface sublimation [mm/s]
  real, intent(in) :: &
       snow_lprec, & ! ?? solid / liquid throughfall infiltrating the snow [mm/s]
       snow_hlprec, & ! ?? heat associated with snow_lprec [W/m^2]
       vegn_uptk, &  ! vegetation soil water uptake flux [mm/s]
       subs_DT,       & ! ?? soil surface layer temperature tendency [K]
       subs_M_imp,       &! rate of phase change of non-evaporated soil water ?? [mm/s]
       subs_evap         ! ?? solution for soil surface evaporation [mm/s]
  logical, intent(in) :: use_tfreeze_in_grnd_latent
  real, intent(out) :: &
       soil_levap, & ! ?? liquid soil surface evaporation [mm/s]
       soil_fevap, & ! ?? solid soil surface sublimation [mm/s]
       soil_melt, &  ! ?? net melt of non-evaporated soil ice over timestep [mm]
       soil_lrunf, & ! ?? total liquid runoff from soil [mm/s]
       soil_hlrunf, & ! ?? heat associated with runoff from soil [W/m^2]
       soil_Ttop, & ! ?? soil surface layer temperature [K]
       soil_Ctop, & ! ?? soil surface layer heat capacity [J/m^2.K]
       soil_frunf, & ! ?? frozen runoff from soil [mm/s]
       soil_hfrunf, & ! ?? heat associated with frozen runoff from soil [W/m^2]
       soil_doc_runf, & ! dissolved carbon tracer runoff from soil [kgC/m^2/s]
       soil_don_runf, & ! dissolved organic N runoff from soil [kgN/m2/s]
       soil_NO3_runf, soil_NH4_runf ! dissolved nitrate and ammonium runoff from soil [kgN/m2/s]


  ! ---- local vars ----------------------------------------------------------
  real, dimension(num_l)   :: del_t, &! ?? temperature tendency [K]
       psi, &! soil moisture potential wrt local elevation [m]
       DThDP, & ! ?? deriv. of vol. liq. cont. wrt psi [1/m]
       K_z, K_x, &! soil horiz. and vert. hydraulic conductivity, respectively [mm/s]
       DKDP, & ! ?? deriv. of hyd. cond. wrt psi [kg/m^3]
       vlc, & ! volumetric liquid water content [-]
       vsc, & ! volumetric solid water content [-]
       dW_l, & ! tendency of soil water mass [mm]
       DPsi
  real, dimension(num_l+1) :: flow, & ! downwards flow at layer interface above [mm/timestep]
       infilt
  real, dimension(num_l  ) :: div, & ! total divergence of soil water [mm/s]
       div_it    ! divergence of water due to inter-tile flow (incl. to stream)
  ! set in hlsp_hydrology_1 [mm/s]
  real, dimension(num_l  ) :: hdiv_it, &! divergence of heat due to inter-tile water flow [W/m^2]
       div_bf, & ! baseflow [mm/s]
       div_if, & ! interlow [mm/s]
       div_al, & ! div from active layer [mm/s]
       dq, div_active, &
       air_depth, macro_frac, extra

  real  :: &
       lprec_eff, & ! infiltrating throughfall (less saturated runoff) [mm/s], and
       hlprec_eff, & !  associated heat [W/m^2]
       tflow, & ! assumed temperature of liquid flow entering top layer [K]
       hcap, & ! heat capacity [J/m^2.K]
       dheat, &
       melt_per_deg, melt, macro_inf, &
       extra_cum, ddW, sum_air, denom, h1, h2, &
       flow_macro, & ! infiltration-excess runoff penetrating the soil due to macroporosity [mm/s]
       lrunf_sn, & ! runoff from saturated ground [mm/s]
       lrunf_ie,lrunf_bf,lrunf_if,lrunf_al, & ! ??, baseflow (incl inter-tile), interflow, active-layer runoff [mm/s]
       lrunf_nu, & ! runoff due to numerical saturation excess [mm/s]
       lrunf_sc, & ! sub-column runoff [mm/s]
       frunf, & ! frozen runoff (due to sat excess) [mm/s]
       hfrunf, & ! heat from frozen runoff (due to sat excess) [W/m^2]
       d_GW, &
       hlrunf_sn,hlrunf_ie,hlrunf_bf,hlrunf_if,hlrunf_al,hlrunf_nu,hlrunf_sc, & ! runoff heat fluxes [W/m^2]
       c0, c1, c2, Dpsi_min, Dpsi_max, &
       sat_area_frac, sat_thick, sum_trans, &
       gw_flux, depth_to_wt2_3, &
       depth_to_gw_flow, depth_to_gw_flow_3, &
       active_layer_thickness, d_psi, d_psi_s, psi_star, &
       depth_to_cf_1, depth_to_cf_3, &
       depth_to_wt_1, depth_to_wt_2, depth_to_wt_2a, depth_to_wt_2b, depth_to_wt_3, depth_to_wt_4, &
       storage_2, deficit_2, deficit_3, deficit_4, deficit, dum1, dum2, dum3
  logical :: stiff
  logical :: hlsp_stiff ! were all the horizontal conductivities zero at the call to hlsp_hydrology_1?
  real :: zimh, ziph, dTr_g(num_l), dTr_s(num_l)
  integer :: n_iter, l, l_max_active_layer
  real :: &
       VRL(num_l), & ! volumetric root length, m/m3
       K_r, & ! root membrame permeability, kg/(m3 s)
       r_r, & ! root radius, m
       bwood, & ! heartwood biomass kg C/m2
       uptake(num_l),   & ! uptake by roots per layer, kg/(m2 s)
       uptake_pos,      & ! sum of the positive uptake, kg/(m2 s)
       uptake_T_new, & ! updated average temperature of uptaken water, deg K
       uptake_T_corr,& ! correction for uptake temperature, deg K
       Tu,           & ! temperature of water taken up from (or added to) a layer, deg K
       psi_x0          ! water potential inside roots (in xylem) at zero depth, m
  ! For testing tridiagonal solution
  real, dimension(num_l)   :: t_soil_tridiag ! soil temperature based on generic tridiagonal solution [K]
  real, dimension(num_l)   :: t_diff ! difference from original advection subroutine [K]

  real :: DOC_leached(n_c_types,num_l), div_DOC_loss(n_c_types,num_l),  &     ! C leaching
         leaflitter_DOC_loss(n_c_types),woodlitter_DOC_loss(n_c_types)        ! Surface litter C leaching loss
 real :: DON_leached(n_c_types,num_l), div_DON_loss(n_c_types,num_l),  &     ! N leaching
        leaflitter_DON_loss(n_c_types),woodlitter_DON_loss(n_c_types)        ! Surface litter N leaching loss
real :: NO3_leached(num_l), div_NO3_loss(num_l),  &     ! NO3 leaching
       leaflitter_NO3_loss,woodlitter_NO3_loss        ! Surface litter NO3 leaching loss
real :: NH4_leached(num_l), div_NH4_loss(num_l),  &     ! NH4 leaching
      leaflitter_NH4_loss,woodlitter_NH4_loss        ! Surface litter NH4 leaching loss

  real :: surface_water ! diagnostic surface water storage [m]
  real :: inundated_frac ! diagnostic inundated area fraction [-]
  real :: wet_frac ! diagnostic wetland area fraction
  real :: sat_runf_frac ! [-] effective saturated fraction used for calculating surface runoff, used
  ! when simple inundation scheme is applied
  real :: flow_s(num_l) ! flow downwards at layer interfaces [mm/s] (flow/delta_time excluding bottom interface)
  real :: reflux        ! [mm/s] upwards flow at soil surface in excess of rejected throughfall
#ifdef ZMSDEBUG
  ! Checking diagnostics from Richards that are used in advection solution.
  real :: w1, w2
#endif
  real, dimension(num_l)   :: wl_before ! water content before call to Richards [mm]
  real, parameter :: wthresh = 1.e-9   ! [mm] tolerance for roundoff error for water balance
  real :: wsum1, wsum2   ! total water stocks for balance check [mm, or kg/m^2]
  real :: sliq, sice     ! for call to soil_tile_stock_pe
  character(len=512) :: mesg ! for error message
  integer :: ipt, jpt, kpt, fpt ! for debug
  real :: surf_DOC_loss(n_c_types)! [kg C/m^2] DOC loss from top soil layer to surface runoff due
                                  ! to efflux
  real :: surf_DON_loss(n_c_types)! [kg N/m^2] DON loss from top soil layer to surface runoff due
                              ! to efflux
  real :: surf_NO3_loss, surf_NH4_loss ! [kg N/m^2] NH4 and NO3 loss from top soil layer to surface runoff due to efflux
  real :: total_C_leaching(num_l) ! [kg C/m^2/s] net total vertical DOC leaching by layer
  real :: total_DOC_div           ! [kg C/m^2/s] net total DOC divergence loss rate
  real :: total_DON_leaching(num_l) ! [kg N/m^2/s] net total vertical DON leaching by layer
  real :: total_DON_div           ! [kg N/m^2/s] net total DON divergence loss rate
  real :: total_NO3_leaching(num_l) ! [kg N/m^2/s] net total vertical NO3 leaching by layer
  real :: total_NO3_div           ! [kg N/m^2/s] net total NO3 divergence loss rate
  real :: total_NH4_leaching(num_l) ! [kg N/m^2/s] net total vertical NH4 leaching by layer
  real :: total_NH4_div           ! [kg N/m^2/s] net total NH4 divergence loss rate
  real :: rhizosphere_frac(num_l)
  real, dimension(num_l) :: passive_ammonium_uptake, passive_nitrate_uptake ! Uptake of dissolved mineral N by roots through water uptake
  ! --------------------------------------------------------------------------
  div_active(:) = 0.0

  !.........................................................................
  if(is_watch_point()) then
     write(*,*) ' ##### soil_step_2 checkpoint 1 #####'
     write(*,*) 'mask    ', .true.
     __DEBUG1__(subs_evap)
     __DEBUG1__(snow_lprec)
     __DEBUG1__(vegn_uptk)
     __DEBUG1__(subs_M_imp)
     call dpri('theta_s ',soil%pars%vwc_sat); write(*,*)
     do l = 1, num_l
        write(*,'(a,i2.2)',advance='NO') 'level=', l
        call dpri(' T =', soil%T(l))
        call dpri(' Th=', (soil%ws(l)+soil%wl(l))/(dens_h2o*dz(l)))
        call dpri(' wl=', soil%wl(l))
        call dpri(' ws=', soil%ws(l))
        call dpri(' gw=', soil%groundwater(l))
        write(*,*)
     enddo
  endif
  !.........................................................................

!  if (do_component_balchecks) then
  ! Sum total water mass at beginning of soil_step_2
  call soil_tile_stock_pe (soil, sliq, sice )
  wsum1 = sliq + sice
!  end if

  ! ---- record fluxes -----------------------------------------------------
  soil_levap  = subs_evap*(1-soil_subl)
  soil_fevap  = subs_evap*   soil_subl
  soil_melt   = subs_M_imp / delta_time

  ! ---- load surface temp change and perform back substitution ------------
  del_t(1) = subs_DT
  soil%T(1) = soil%T(1) + del_t(1)
  if ( num_l > 1) then
     do l = 1, num_l-1
        del_t(l+1) = soil%e(l) * del_t(l) + soil%f(l)
        soil%T(l+1) = soil%T(l+1) + del_t(l+1)
     enddo
  endif

  !.........................................................................
  if(is_watch_point()) then
     write(*,*) ' ##### soil_step_2 checkpoint 2 #####'
     do l = 1, num_l
        write(*,'(a,i2.2)',advance='NO') 'level=',l
        call dpri('T=', soil%T(l))
        call dpri('del_t=', del_t(l))
        call dpri('e=', soil%e(l))
        call dpri('f=', soil%f(l))
        write(*,*)
     enddo
  endif
  !.........................................................................

  ! ---- extract evap from soil, adjusting T, and do implicit melt ---------
  IF (LM2) THEN  ! (extract surface E--is there any?--uniformly from bucket)
     do l = 1, num_l
        soil%wl(l) = soil%wl(l) &
                      - soil%uptake_frac(l)*soil_levap*delta_time
     enddo
  ELSE
     soil%wl(1) = soil%wl(1) - soil_levap*delta_time
     soil%ws(1) = soil%ws(1) - soil_fevap*delta_time
  ENDIF
  hcap = soil%heat_capacity_dry(1)*dz(1) &
                       + clw*soil%wl(1) + csw*soil%ws(1)
  ! T adjustment for nonlinear terms (del_T)*(del_W)
  dheat = delta_time*(clw*soil_levap+csw*soil_fevap)*del_T(1)
  ! take out extra heat not claimed in advance for evaporation
  if (use_tfreeze_in_grnd_latent) dheat = dheat &
          - delta_time*((cpw-clw)*soil_levap+(cpw-csw)*soil_fevap) &
                                 *(soil%T(1)-del_T(1)-tfreeze)
  soil%T(1)  = soil%T(1)  + dheat/hcap
  soil%wl(1) = soil%wl(1) + subs_M_imp
  soil%ws(1) = soil%ws(1) - subs_M_imp
  soil%T(1)  = tfreeze + (hcap*(soil%T(1)-tfreeze) ) &
                              / ( hcap + (clw-csw)*subs_M_imp )

  ! ---- calculate actual uptake and update its T --------------------------
  select case(uptake_option)
  case ( UPTAKE_LINEAR )
     uptake_T_corr = 0
     n_iter = 0
     uptake = soil%uptake_frac*vegn_uptk
     soil%psi_x0 = 0.
     bwood = 0.
  case ( UPTAKE_DARCY2D, UPTAKE_DARCY2D_LIN )
     ! for Darcy-flow uptake, find the root water potential to satify actual
     ! transpiration by the vegetation
     ! **** introduce optional arguments for dum1 etc ? ********
     call vegn_hydraulic_properties (vegn, dz(1:num_l), always_use_bsw, &
                               VRL, K_r, r_r, dum1, dum2, dum3, bwood)

     call darcy2d_uptake_solver (soil, vegn_uptk, VRL, K_r, r_r, &
             uptake_option, uptake_oneway, uptake_from_sat, uptake, psi_x0, n_iter)
     soil%psi_x0 = psi_x0

     uptake_pos = sum(uptake(:),mask=uptake(:)>0)
     if (uptake_pos > 0) then
        ! calculate actual temperature of uptake
        uptake_T_new  = sum(uptake*soil%T,mask=uptake>0)/uptake_pos
        ! and temperature correction
        uptake_T_corr = soil%uptake_T - uptake_T_new
        if(is_watch_point()) then
           __DEBUG3__(soil%uptake_T, uptake_T_new, uptake_T_corr)
        endif
        ! save new uptake for the next time step to serve as an estimate of uptake
        ! temperature
        soil%uptake_T = uptake_T_new
     else
        uptake_T_corr = 0.0
        ! and do not change the soil%uptake_T
     endif
  case default
     call error_mesg('soil_step_2', 'invalid soil uptake option', FATAL)
  end select

  !.........................................................................
  if (is_watch_point())then
     write(*,*) ' ##### soil_step_2 checkpoint 2.1 #####'
     __DEBUG2__(vegn_uptk,sum(uptake))
     do l = 1,num_l
        write(*,'(a,i2.2)',advance='NO')'level=',l
        call dpri('uptake=',uptake(l))
        call dpri('dwl=',-uptake(l)*delta_time)
        call dpri('wl=',soil%wl(l))
        call dpri('new wl=',soil%wl(l) - uptake(l)*delta_time)
        write(*,*)
     enddo
  endif
  !.........................................................................

  call send_tile_data(id_uptk_n_iter, real(n_iter), diag)
  call send_tile_data(id_uptk, uptake, diag)
  call send_tile_data(id_psi_x0, psi_x0, diag)

  ! ---- perform the uptake ------------------------------------------------
  do l = 1, num_l
     ! calculate the temperature of water that is taken from the layer (or added
     ! to the layer), including energy balance correction
     if (uptake(l) > 0) then
        Tu = soil%T(l) + uptake_T_corr
     else
        Tu = soil%uptake_T + uptake_T_corr
     endif
     hcap = soil%heat_capacity_dry(l)*dz(l) &
           + clw*soil%wl(l) + csw*soil%ws(l)
     soil%T(l) = soil%T(l) - &
           uptake(l)*delta_time*clw*( Tu-soil%T(l) ) / &
           ( hcap - uptake(l)*delta_time*clw )
     soil%wl(l) = soil%wl(l) - uptake(l)*delta_time
  enddo

  !.........................................................................
  if(is_watch_point()) then
     write(*,*) ' ##### soil_step_2 checkpoint 3 #####'
     do l = 1, num_l
        write(*,'(x,a,x,i2.2)',advance='NO')' level=', l
        call dpri(' T =',soil%T(l))
        call dpri(' Th=', (soil%ws(l)+soil%wl(l))/(dens_h2o*dz(l)))
        call dpri(' wl=', soil%wl(l))
        call dpri(' ws=', soil%ws(l))
        write(*,*)
     enddo
  endif
  !.........................................................................

  ! ---- push down any excess surface water, with heat ---------------------
  IF (PUSH_DOWN_SFC_EXCESS) THEN
     CALL SOIL_PUSH_DOWN_EXCESS ( soil, diag, lrunf_nu, hlrunf_nu, frunf, hfrunf)
  ELSE
     lrunf_nu=0; hlrunf_nu=0; frunf=0; hfrunf = 0
  ENDIF

  ! ---- fetch soil hydraulic properties -----------------------------------
  do l = 1, num_l
    vlc(l) = max(0., soil%wl(l) / (dens_h2o*dz(l)))
    vsc(l) = max(0., soil%ws(l) / (dens_h2o*dz(l)))
  enddo
  call soil_data_hydraulic_properties (soil, vlc, vsc, &
                   psi, DThDP, K_z, K_x, DKDP, Dpsi_min, Dpsi_max )

  ! ---- compute various measures of water table depth ---------------------
  sat_thick = 0.
  do l=num_l,1,-1
     if(vsc(l)+vlc(l).le.soil%pars%vwc_sat) exit
     sat_thick = sat_thick + dz(l)
  enddo
  depth_to_cf_1 = zhalf(num_l+1) - sat_thick
  depth_to_wt_1 = depth_to_cf_1 - soil%pars%psi_sat_ref/soil%alpha(max(l,1))

  depth_to_wt_2 = zfull(num_l)-psi(num_l)

  depth_to_wt_2a = 0.
  depth_to_wt_3  = 0.0
  do l=1,num_l
     if (soil%wl(l)+soil%ws(l) .lt. &
                      soil%pars%vwc_sat*dens_h2o*dz(l)) then
        depth_to_wt_2a = depth_to_wt_2a + dz(l)
        if (gw_option /= GW_TILED .and. l.eq.num_l) depth_to_wt_2a = -1.
        ! ZMS changed definition here so time-average will remain deep water table
        ! and for use in calculating inundated fraction in hlsp_hydrology_2
     else
        exit
     endif
  enddo
  depth_to_wt_2b = 0. ! This version requires an entirely unfrozen saturated layer
  do l=1,num_l
     if (soil%wl(l) .lt. soil%pars%vwc_sat*dens_h2o*dz(l)) then
        depth_to_wt_2b = depth_to_wt_2b + dz(l)
     else
        exit
     end if
  end do

  storage_2 = 1 - depth_to_wt_2  &
           /(soil%pars%hillslope_zeta_bar*soil%pars%hillslope_relief)
  storage_2 = min( max( 0., storage_2 ) , 1.)
  deficit_2 = 1 - storage_2

  if (vsc(num_l).gt.0.) then   ! permafrost
     if (gw_option /= GW_TILED) then
        depth_to_wt2_3 = 0.
        depth_to_cf_3 = 0.
        depth_to_wt_3 = 0.
     else
        depth_to_wt2_3 = zhalf(l+1)
        depth_to_cf_3 = zhalf(l+1)
        depth_to_wt_3 = zhalf(l+1)
     end if
  else                       ! liquid water at depth
     depth_to_cf_3 = 0.
     if (use_fringe) then
        do l = num_l, 1, -1
           if ( l.eq.num_l .and. psi(l).le.soil%pars%psi_sat_ref/soil%alpha(l) ) then
              depth_to_cf_3 = zfull(l) + soil%pars%psi_sat_ref/soil%alpha(l) - psi(l)
              exit
           else if (psi(l).le.soil%pars%psi_sat_ref/soil%alpha(l) ) then
              d_psi = psi(l+1) - psi(l)
              d_psi_s = (soil%pars%psi_sat_ref/soil%alpha(l+1)) &
                       -(soil%pars%psi_sat_ref/soil%alpha(l))
              psi_star = (psi(l)*d_psi_s - &
                          d_psi*(soil%pars%psi_sat_ref/soil%alpha(l)))&
                          / (d_psi_s - d_psi)
              depth_to_cf_3 = zfull(l) + (zfull(l+1)-zfull(l)) &
                                    * (psi_star-psi(l)) / d_psi
              exit
           else if (l.eq.1) then
              d_psi = psi(l+1) - psi(l)
              d_psi_s = (soil%pars%psi_sat_ref/soil%alpha(l+1)) &
                       -(soil%pars%psi_sat_ref/soil%alpha(l))
              psi_star = (psi(l)*d_psi_s - &
                          d_psi*(soil%pars%psi_sat_ref/soil%alpha(l)))&
                          / (d_psi_s - d_psi)
              depth_to_cf_3 = zfull(l) + (zfull(l+1)-zfull(l)) &
                                    * (psi_star-psi(l)) / d_psi
              depth_to_cf_3 = max(0.,depth_to_cf_3)
           endif
        enddo
     endif
     depth_to_wt_3 = max(0., zhalf(num_l+1)-(psi(num_l)+dz(num_l)/2.))
     depth_to_wt2_3 = depth_to_wt_3
  endif

  if (use_fringe) then
     depth_to_gw_flow_3 = depth_to_cf_3
  else
     depth_to_gw_flow_3 = depth_to_wt_3
  endif
  deficit_3 = depth_to_gw_flow_3  &
           /(soil%pars%hillslope_zeta_bar*soil%pars%hillslope_relief)

  if (vsc(min(num_l,layer_for_gw_switch)).gt.0.) then   ! permafrost
     if (gw_option /= GW_TILED) then
        depth_to_wt_4 = 0.
     else
        depth_to_wt_4 = zhalf(num_l+1)
     end if
  else ! liquid water at depth
     depth_to_wt_4 = 0.
     do l = num_l, 1, -1
        if ( l.eq.num_l .and. psi(l).le.0. ) then
           depth_to_wt_4 = zfull(l) - psi(l)
           exit
        else if (psi(l).le.0. ) then
           d_psi = psi(l+1) - psi(l)
           depth_to_wt_4 = zfull(l) - (zfull(l+1)-zfull(l)) &
                                 * (psi(l)) / d_psi
           exit
        else if (l.eq.1) then
           d_psi = psi(l+1) - psi(l)
           depth_to_wt_4 = zfull(l) - (zfull(l+1)-zfull(l)) &
                                 * (psi(l)) / d_psi
           depth_to_wt_4 = max(0.,depth_to_wt_4)
        endif
     enddo
  endif
  deficit_4 = depth_to_wt_4  &
           /(soil%pars%hillslope_zeta_bar*soil%pars%hillslope_relief)

! ---- get saturated area and column flow divergences --------------------
  SELECT CASE(gw_option)

  CASE(GW_LM2)

   div_bf=0; div_if=0; div_al=0; sat_area_frac = 0; div_it = 0.; hdiv_it = 0.

  CASE(GW_LINEAR)

      IF (CORRECTED_LM2_GW) THEN
         do l = 1, num_l
            if (vlc(l) .ge. soil%pars%vwc_sat .and. vsc(l).le.0.) &
               div_bf(l) = 0.15*dens_h2o*dz(l)/soil%pars%tau_groundwater
         enddo
      ELSE
         do l = 1, num_l
            if ((vsc(l)+vlc(l)) .ge. soil%pars%vwc_sat) &
                div_bf(l) = 0.15*dens_h2o*dz(l)*(vlc(l)/(vsc(l)+vlc(l)))  &
                                    /soil%pars%tau_groundwater
         enddo
      ENDIF
      div_if = 0
      div_al = 0; div_it = 0.; hdiv_it = 0.
      sat_thick = zhalf(num_l+1) - depth_to_cf_1
      sat_area_frac = min((sat_thick/zhalf(num_l+1))**soil%pars%rsa_exp,1.)

  CASE(GW_HILL_AR5)

      call soil_data_gw_hydraulics_ar5(soil, storage_2, &
                                               gw_flux, sat_area_frac)
      dq = 0.
      sat_thick = 0.
      do l=num_l,1,-1
         if(psi(l).le.0.) exit
         if (vsc(l).le.0.) dq(l) = dz(l)
         sat_thick = sat_thick + dz(l)
      enddo
      div_bf = 0.
      if (sat_thick.gt.0.) div_bf = (dq/sat_thick)*gw_flux

      div_active = 0.
      l_max_active_layer = 0
      do l=1,num_l
         if(vsc(l).gt.0.) exit
         l_max_active_layer = l
      enddo
      if (l_max_active_layer.lt.num_l .and. l_max_active_layer.gt.0) then
          do l = 1, l_max_active_layer
             if(vlc(l).gt.0) &
             div_active(l) = K_x(l) * soil%pars%hillslope_relief*dz(l) &
              / (soil%pars%hillslope_length*soil%pars%hillslope_length)
          enddo
      endif

      div_al = 0
      where (div_bf.eq.0.) div_al = div_active*active_layer_drainage_acceleration
   div_if = 0; div_it = 0.; hdiv_it = 0.

  CASE(GW_HILL)
     if (any(soil%wl+soil%ws.lt.0.) .and. prohibit_negative_water_div) then
        sat_area_frac = 0.
        div_bf = 0.
        div_if = 0.
     else if (vsc(min(num_l,layer_for_gw_switch)).gt.0.) then   ! permafrost
        sat_area_frac = 0.
        div_bf = 0.
        div_if = 0.
     else                       ! liquid water at depth
        if (use_depth_to_wt_4) then
           depth_to_gw_flow = depth_to_wt_4
           deficit = deficit_4
        else
           depth_to_gw_flow = depth_to_gw_flow_3
           deficit = deficit_3
        endif
        call soil_data_gw_hydraulics(soil, deficit, gw_flux, sat_area_frac)
        gw_flux = min(gw_flux, gw_flux_max)
        dTr_g = 0.
        dTr_s = 0.
        dTr_g(num_l) = 1.
        l = num_l
        ziph = sum(dz)
        zimh = ziph - dz(num_l)
        if (depth_to_gw_flow .lt. zimh) then
           dTR_g(l) = dz(l)
           dTr_s(l) = (exp(-zimh/soil%pars%soil_e_depth))
           do l = num_l-1, 1, -1
              if (vsc(l).gt.0.) exit
              ziph = zimh
              zimh = ziph - dz(l)
              if (depth_to_gw_flow .lt. zimh) then
                 dTR_g(l) = dz(l)
                 dTr_s(l) = exp(-zimh/soil%pars%soil_e_depth)-exp(-ziph/soil%pars%soil_e_depth)
              else if (depth_to_gw_flow .lt. ziph) then
                 dTR_g(l) =(ziph-depth_to_gw_flow)
                 dTr_s(l) = exp(-depth_to_gw_flow/soil%pars%soil_e_depth)-exp(-ziph/soil%pars%soil_e_depth)
              else
                 exit
              endif
           enddo
        endif
        sum_trans = sum(dTr_g)
        if (sum_trans.ne.0.) then
           dTR_g = dTR_g / sum_trans
           dTR_g = dTR_g * soil%pars%k_sat_gw*aspect*soil%pars%hillslope_length
        endif
        dTR_s = dTR_s * (soil%pars%k_sat_sfc+k0_macro_x)*soil%pars%soil_e_depth
        sum_trans = sum(dTR_g) + sum(dTr_s)
        if (sum_trans.ne.0.) then
           div_bf = gw_flux * dTR_g /sum_trans
           div_if = gw_flux * dTR_s /sum_trans
        else
           div_bf = 0.
           div_if = 0.
        endif
     endif

     div_al = 0; div_it = 0.; hdiv_it = 0.
     l_max_active_layer = 0   ! "active layer" either over permafrost or perched
     do l=1,num_l
        if(vsc(l).gt.0.) exit
        l_max_active_layer = l
     enddo
     if (l_max_active_layer.lt.num_l .and. l_max_active_layer.gt.0) then
        do l = 1, l_max_active_layer
           div_al(l) = K_x(l) * soil%pars%hillslope_relief*dz(l) &
                 / (soil%pars%hillslope_length*soil%pars%hillslope_length)
        enddo
     endif

  CASE (GW_TILED)
     call hlsp_hydrology_2(soil, psi, vlc, vsc, div_it, hdiv_it, &
        sat_area_frac, inundated_frac, storage_2, depth_to_wt_2b, surface_water, sat_runf_frac)
     if (depth_to_wt_2b >= 0. .and. depth_to_wt_2b <= wet_depth) then
        wet_frac = 1.
     else
        wet_frac = 0.
     end if
     div_if=0.
     div_al=0.
     div_bf=0.

  END SELECT

  div = div_bf + div_if + div_al + div_it ! div includes inter-tile flow
  lrunf_bf = sum(div_bf + div_it) ! baseflow runoff includes inter-tile flow
  lrunf_if = sum(div_if)
  lrunf_al = sum(div_al)

  if (snow_lprec.ne.0.) then
     if (gw_option == GW_TILED .and. simple_inundation) then
        ! use sat_runf_frac calculated in hlsp_hydrology_2. Includes microtopography and finite runoff rate.
        lrunf_sn = sat_runf_frac * snow_lprec
     else
        lrunf_sn = sat_area_frac * snow_lprec
        ! if gw_option == GW_TILED, sat_area_frac = 0 or 1 depending on saturation status of top layer
     end if
     hlrunf_sn = lrunf_sn*snow_hlprec/snow_lprec
  else
     lrunf_sn = 0.
     hlrunf_sn = 0.
  endif
  hlrunf_ie=0
  lprec_eff = snow_lprec - lrunf_sn
  hlprec_eff = snow_hlprec - hlrunf_sn

  if(is_watch_point()) then
     do l = 1, num_l
        write(*,'(a,1x,i2.2,100(2x,g23.16))')'div_ac,div_bf,div_if,div_al,div', &
                       l,div_active(l),div_bf(l),div_if(l),div_al(l),div(l)
     enddo
     do l = 1, num_l
        write(*,'(a,1x,i2.2,100(2x,g23.16))')'vsc,psi,dz',l,vsc(l),psi(l),dz(l)
     enddo
     write(*,*)'lrunf_bf',lrunf_bf
     write(*,*)'tau_gw',soil%pars%tau_groundwater
     write(*,*)'dens_h2o',dens_h2o
  endif

  wl_before(1:num_l) = soil%wl(1:num_l)

  ! ---- soil-water flow ----------------------------------------------------
  IF (LM2) THEN
     flow(1) = 0
     do l = 1, num_l
        infilt(l) = soil%uptake_frac(l)*lprec_eff *delta_time
        flow(l+1) = max(0., soil%wl(l) + flow(l) &
              + infilt(l) - soil%w_fc(l)*dz(l)*dens_h2o)
        dW_l(l) = flow(l) - flow(l+1) + infilt(l)
        soil%wl(l) = soil%wl(l) + dW_l(l)
     enddo
     do l = 1, num_l
        flow(l) = flow(l) + infilt(l)
     enddo
     dW_l=0
     dpsi=0
     c0 = delta_time/soil%pars%tau_groundwater
     c1 = exp(-c0)
     c2 = (1-c1)/c0
     l = 1
     d_GW = c1 * soil%groundwater(l) + c2 * flow(num_l+1) &
                           - soil%groundwater(l)
     soil%groundwater(l) = soil%groundwater(l) + d_GW
     lrunf_sc  = (1-c1)*soil%groundwater(l)/delta_time &
                         + (1-c2)*flow(num_l+1)/delta_time
     lrunf_ie=0
  ELSE
     lrunf_sc = 0
     d_GW = 0
     stiff = all(DThDP.eq.0)
     IF(stiff) THEN
        ! Note: in order to maintain energy and water conservation with between tile-flows
        ! in case of the tiled groundwater / full hillslope model, we require stiff ==>
        ! soil%T(l) < tfreeze for all l.
        ! ZMS: since this could change during timestep, impose some explicit updates in case this
        ! is now the case but the flows are nonzero.
        hlsp_stiff = all(div_it==0.) .and. all(hdiv_it==0.)
        if (is_watch_cell()) then
           write(*,*)'Stiff point in watch-cell at hidx_j=', soil%hidx_j
           __DEBUG2__(div_it, hdiv_it)
           if (hlsp_stiff) write(*,*) 'hlsp_stiff'
        end if
        if (.not. hlsp_stiff) then
           call stiff_explicit_gwupdate(soil, div_it, hdiv_it, div, lrunf_bf)
        end if
        flow = 0.
        dW_l = 0.
        if (.not. hlsp_stiff) then
           div  = 0.
           lrunf_bf = 0.
        end if
        div_bf=0.; div_if=0; div_al=0
        lrunf_if = 0; lrunf_al = 0
        lrunf_ie = lprec_eff
        hlrunf_ie = hlprec_eff
        psi=zfull(1:num_l)
        dpsi=0.
     ELSE
        if (USE_RICHARDS_CLEAN) then
           call RICHARDS_clean(soil, psi, DThDP, K_z, DKDP, div, &
                  lprec_eff, Dpsi_min, Dpsi_max, delta_time, &
                  dPsi, dW_l, flow, lrunf_ie)
        else
           call RICHARDS(soil, psi, DThDP, K_z, DKDP, div, &
                  lprec_eff, Dpsi_min, Dpsi_max, delta_time, &
                  dPsi, dW_l, flow, lrunf_ie)
        endif
        do l = 1, num_l
           soil%wl(l) = soil%wl(l) + dW_l(l)
        enddo
#ifdef ZMSDEBUG
        do l = 1, num_l
           w1 = wl_before(l)
           w2 = soil%wl(l) - dW_l(l)
           call check_conservation('soil_step_2: Richards Eqn. Diagnostics, wl_before', 'Water', w1, w2, &
                wthresh, WARNING)
           w1 = dW_l(l)
           if (l < num_l) then
              w2 = flow(l) - flow(l+1) - div(l)*delta_time
           else
              w2 = flow(l) - div(l)*delta_time
           end if
           call check_conservation('soil_step_2: Richards Eqn. Diagnostics, dW_l', 'Water', w1, w2, &
                wthresh, WARNING)
        end do
#endif
      ENDIF
  ENDIF

  ! Check for negative wl
  do l = 1, num_l
     if (soil%wl(l) < 0.) then
        if (soil%wl(l) / (dens_h2o*dz(l)*soil%pars%vwc_sat) < thetathresh) then
           call get_current_point(ipt, jpt, kpt, fpt)
           write(mesg,*) 'soil%wl(l) < 0! l,i,j,k,face:', l, ipt, jpt, kpt, fpt, '. degree of saturation = ', &
                soil%wl(l) / (dens_h2o*dz(l)*soil%pars%vwc_sat), '. If ".not. allow_neg_wl", '// &
                'model will abort.'
           if (.not. allow_neg_wl) then
              call error_mesg(module_name, mesg, FATAL)
           else
              if (verbose) call error_mesg(module_name, mesg, WARNING)
           end if
        end if
     end if
  end do

  ! Check for negative hcap
  if (require_pos_hcap) then
     do l = 1, num_l
        if (soil%wl(l) < 0.) then
            ! Make sure bulk heat capacity stays above zero
            hcap = soil%heat_capacity_dry(l)*dz(l) &
                   + clw*soil%wl(l) + csw*soil%ws(l)
            if (hcap .le. 0.) then
               call get_current_point(ipt, jpt, kpt, fpt)
               write(mesg,*) 'soil%wl(l) < 0! l,i,j,k,face:', l, ipt, jpt, kpt, fpt, '. degree of saturation = ', &
                      soil%wl(l) / (dens_h2o*dz(l)*soil%pars%vwc_sat), '. This makes hcap = ', hcap, &
                      ', which is < 0! Model Aborting!'
               call error_mesg(module_name, mesg, FATAL)
            end if
        end if
     end do
  endif
  ! Check total
  if (lrunf_nu < negrnuthresh) then
     call get_current_point(ipt, jpt, kpt, fpt)
     write(mesg,*) 'soil%wl(l) < 0 at one or more l at i,j,k,face:', ipt, jpt, kpt, fpt, '.', &
          ' Total lrunf_nu required to set to zero = ', lrunf_nu, ' mm/s.'
     call error_mesg(module_name, mesg, WARNING)
  end if


  ! ---- heat advection by water flow ---------------------------------------
  if  (snow_lprec.ne.0.) then
    tflow = tfreeze + snow_hlprec/(clw*snow_lprec)
  else
    tflow = tfreeze
  endif

  if(is_watch_point()) then
     write(*,*) ' ***** soil_step_2 checkpoint 3.4 ***** '
     write(*,*) ' tfreeze', tfreeze
     write(*,*) '  tflow ', tflow
     write(*,*) ' snow_hlprec', snow_hlprec
  endif

#ifndef ZMSDEBUG_TRIDIAGTEST
  if (use_tridiag_foradvec) then
     call advection_tri(soil, flow, dW_l, tflow, d_GW, div, delta_time, t_soil_tridiag, hdiv_it)
  else
     call advection(soil, flow, dW_l, tflow, d_GW, div, delta_time)
  end if
#else
  call advection_tri(soil, flow, dW_l, tflow, d_GW, div, delta_time, t_soil_tridiag, hdiv_it)
  call advection(soil, flow, dW_l, tflow, d_GW, div, delta_time)
  ! Calculate difference between two solutions
  t_diff(1:num_l) = t_soil_tridiag(1:num_l) - soil%T(1:num_l)
  call send_tile_data(id_st_diff, t_diff, diag)
#endif

  if (lprec_eff.ne.0. .and. flow(1).ge.0. ) then
     hlrunf_ie = lrunf_ie*hlprec_eff/lprec_eff
  else if (flow(1).lt.0. ) then
     hlrunf_ie = hlprec_eff - (flow(1)/delta_time)*clw &
                         *(soil%T(1)-tfreeze)
  else
     hlrunf_ie = 0.
  endif

  ! Initialize for use in output below
  macro_inf = 0.
  extra_cum = 0.
  ! ---- allow infiltration-excess runoff to enter soil via macroporosity
  IF (BOTTOM_UP_COLD_INFILT.and.bwood.gt.bwood_macinf.and.lrunf_ie.gt.0.) then
     do l = 1, num_l
        air_depth(l) = soil%pars%vwc_sat*dz(l)-(soil%wl(l)+soil%ws(l))/dens_h2o
        air_depth(l) = min(air_depth(l), soil%pars%vwc_sat*dz(l))
        air_depth(l) = max(air_depth(l), 0.)
        if (zfull(l).gt.COLD_DEPTH) air_depth(l) = 0.
     enddo
     sum_air = sum(air_depth)
     macro_inf = min(lrunf_ie*delta_time,sum_air*dens_h2o)
     tflow = tfreeze + hlrunf_ie/(clw*lrunf_ie)
     ! compute the fill amounts from bottom up
     extra_cum = macro_inf
     dW_l = 0.
     do l = num_l, 1, -1
        if (extra_cum.gt.0.) then
           dW_l(l) = min(extra_cum, air_depth(l)*dens_h2o)
           extra_cum = extra_cum - dW_l(l)
        else
           exit
        endif
     enddo
     ! place the water and heat
     do l= 1, num_l
        h1 = soil%heat_capacity_dry(l)*dz(l) &
             + csw*soil%ws(l) + clw*soil%wl(l)
        h2 = clw*dW_l(l)
        soil%T(l) = (h1 * soil%T(l) &
                      + h2 * tflow )  / (h1+h2)
        soil%wl(l) = soil%wl(l) + dW_l(l)
     enddo
     ! extra_cum should now be 0, but included below for safety
     lrunf_ie = lrunf_ie - (macro_inf-extra_cum)/delta_time
     hlrunf_ie = hlrunf_ie - clw*(tflow-tfreeze)*(macro_inf-extra_cum)/delta_time
  ELSE if (cold_infilt.and.lrunf_ie.gt.0.) then
     do l = 1, num_l
        air_depth(l) = soil%pars%vwc_sat*dz(l)-(soil%wl(l)+soil%ws(l))/dens_h2o
        air_depth(l) = min(air_depth(l), soil%pars%vwc_sat*dz(l))
        air_depth(l) = max(air_depth(l), 0.)
        macro_frac(l) = exp(-zhalf(l)  /soil%pars%soil_e_depth) &
                       -exp(-zhalf(l+1)/soil%pars%soil_e_depth)
     enddo
     sum_air = sum(air_depth)
     macro_inf = min(lrunf_ie*delta_time,sum_air*dens_h2o)
     tflow = tfreeze + hlrunf_ie/(clw*lrunf_ie)
     if (sum_air.gt.0.) then
        denom = sum(air_depth*macro_frac/dz)
        ! tentatively distribute macropore infiltration in proportion
        ! to available pore volume and macropore density
        do l= 1, num_l
           dW_l(l) = macro_inf*(air_depth(l)*macro_frac(l)/dz(l))/denom
        enddo
        extra = 0.
        ! find excess infiltration by layer
        do l = 1, num_l
           extra(l) = max(dW_l(l) - air_depth(l)*dens_h2o, 0.)
           dW_l(l) = dW_l(l) - extra(l)
        enddo
        ! sweep any excess downward
        extra_cum = 0.
        do l = 1, num_l
           extra_cum = extra_cum + extra(l)
           if (air_depth(l)*dens_h2o.gt.dW_l(l)) then
              ddW = min(extra_cum, air_depth(l)*dens_h2o-dW_l(l))
              extra_cum = extra_cum - ddw
              dW_l(l) = dW_l(l) + ddW
           endif
        enddo
        ! sweep any remaining excess upward
        do l = num_l, 1, -1
           extra_cum = extra_cum + extra(l)
           if (air_depth(l)*dens_h2o.gt.dW_l(l)) then
              ddW = min(extra_cum, air_depth(l)*dens_h2o-dW_l(l))
              extra_cum = extra_cum - ddw
              dW_l(l) = dW_l(l) + ddW
           endif
        enddo
        ! place the water and heat
        do l= 1, num_l
           h1 = soil%heat_capacity_dry(l)*dz(l) &
                + csw*soil%ws(l) + clw*soil%wl(l)
           h2 = clw*dW_l(l)
           soil%T(l) = (h1 * soil%T(l) &
                               + h2 * tflow )  / (h1+h2)
           soil%wl(l) = soil%wl(l) + dW_l(l)
        enddo
     endif
     ! extra_cum should now be 0, but included below for safety
     lrunf_ie = lrunf_ie - (macro_inf-extra_cum)/delta_time
     hlrunf_ie = hlrunf_ie - clw*(tflow-tfreeze)*(macro_inf-extra_cum)/delta_time
  ENDIF

  flow_macro = (macro_inf-extra_cum)/delta_time

  hlrunf_bf = clw*sum(div_bf*(soil%T-tfreeze)) + sum(hdiv_it)
  ! div_bf is 0. for GW_TILED, else hdiv_it == zero
  hlrunf_if = clw*sum(div_if*(soil%T-tfreeze))
  hlrunf_al = clw*sum(div_al*(soil%T-tfreeze))
  hlrunf_sc = clw*lrunf_sc  *(soil%groundwater_T(1)-tfreeze)
  if (lrunf_from_div) then
     soil_lrunf  =  lrunf_sn +  lrunf_ie +  sum(div) +  lrunf_nu +  lrunf_sc
     if (gw_option /= GW_TILED) then
        soil_hlrunf = hlrunf_sn + hlrunf_ie +  clw*sum(div*(soil%T-tfreeze)) &
                                                      + hlrunf_nu + hlrunf_sc
     else
        soil_hlrunf = hlrunf_sn + hlrunf_ie +  sum(hdiv_it) &
             + hlrunf_nu + hlrunf_sc
     end if
  else
     soil_lrunf  =  lrunf_sn +  lrunf_ie +  lrunf_bf +  lrunf_if &
                             +  lrunf_al +  lrunf_nu +  lrunf_sc
     soil_hlrunf = hlrunf_sn + hlrunf_ie + hlrunf_bf + hlrunf_if &
                             + hlrunf_al + hlrunf_nu + hlrunf_sc
  endif
  soil_frunf = frunf
  soil_hfrunf = hfrunf

  if (is_watch_cell()) then
     write(*,*)'Soil runoff from point in watch_cell'
     __DEBUG3__(soil%hidx_j, soil_lrunf, lrunf_bf)
     __DEBUG5__(lrunf_sn, lrunf_ie, lrunf_if, lrunf_al, lrunf_sc)
  end if


  do l = 1, num_l
     ! ---- compute explicit melt/freeze --------------------------------------
     hcap = soil%heat_capacity_dry(l)*dz(l) &
              + clw*soil%wl(l) + csw*soil%ws(l)
     melt_per_deg = hcap/(hlf_factor*hlf)
     if       (soil%ws(l)>0 .and. soil%T(l)>soil%pars%tfreeze) then
       melt =  min(soil%ws(l), (soil%T(l)-soil%pars%tfreeze)*melt_per_deg)
     else if (soil%wl(l)>0 .and. soil%T(l)<soil%pars%tfreeze) then
       melt = -min(soil%wl(l), (soil%pars%tfreeze-soil%T(l))*melt_per_deg)
     else
       melt = 0
     endif
     soil%wl(l) = soil%wl(l) + melt
     soil%ws(l) = soil%ws(l) - melt
     soil%T(l) = tfreeze &
        + (hcap*(soil%T(l)-tfreeze) - hlf_factor*hlf*melt) &
                             / ( hcap + (clw-csw)*melt )
     soil_melt = soil_melt + melt / delta_time
  enddo

  if(is_watch_point()) then
     write(*,*) ' ##### soil_step_2 checkpoint 5 #####'
     do l = 1, num_l
        write(*,'(a,i2.2,100(2x,a,g23.16))') ' level=', l,&
             ' T =', soil%T(l),&
             ' Th=', (soil%ws(l) +soil%wl(l))/(dens_h2o*dz(l)),&
             ' wl=', soil%wl(l),&
             ' ws=', soil%ws(l),&
             ' gw=', soil%groundwater(l)
     enddo
     call debug_pool(soil%leafLitter, 'leafLitter')
  endif

  active_layer_thickness = 0.
  do l = 1, num_l
     if (soil%ws(l).gt.0.) then
        active_layer_thickness = active_layer_thickness &
          + dz(l)*max(0.,soil%wl(l))/(soil%wl(l)+soil%ws(l))
        exit
     endif
     active_layer_thickness = active_layer_thickness + dz(l)
  enddo

  soil_Ttop = soil%T(1)
  soil_Ctop = soil%heat_capacity_dry(1)*dz(1) &
    + clw*soil%wl(1) + csw*soil%ws(1)

  soil%psi=psi+dPsi
!  if (do_component_balchecks) then
     ! Sum total water mass at end of soil_step_2
     call soil_tile_stock_pe (soil, sliq, sice )
     wsum2 = sliq + sice

     ! Subtract influxes and add outfluxes
     wsum2 = wsum2 - delta_time *  snow_lprec &
          + delta_time * ( subs_evap + vegn_uptk + soil_lrunf + soil_frunf)

     call check_conservation('soil_mod: soil_step_2', 'Water', wsum1, wsum2, wthresh, FATAL)
! endif


   if (is_watch_point()) then
      write(*,*)'##### soil_step_2 checkpoint 6 #####'
      __DEBUG1__(flow)
      __DEBUG1__(div)
      __DEBUG1__(wl_before)
      __DEBUG1__(gw_option)
      call debug_pool(soil%leafLitter,       'leafLitter')
      call debug_pool(soil%fineWoodLitter,   'fineWoodLitter')
      call debug_pool(soil%coarseWoodLitter, 'coarseWoodLitter')
      do l = 1, num_l
         call debug_pool(soil%soil_organic_matter(l), 'soil_organic_matter(l)')
      enddo
      do l = 1, size(soil%div_hlsp_DOC,2)
         __DEBUG1__(soil%div_hlsp_DOC(:,l))
      enddo
      do l = 1, size(soil%div_hlsp_DON,2)
         __DEBUG1__(soil%div_hlsp_DON(:,l))
      enddo
      do l=1, size(soil%div_hlsp_NO3)
          __DEBUG1__(soil%div_hlsp_NO3(l))
          __DEBUG1__(soil%div_hlsp_NH4(l))
      enddo
   endif


! Do plant nitrogen uptake
! Passive uptake by roots, equal to N concentration times root water uptake
! units of uptake: kg/m2/s
! units of wl_before: mm = kg/m2
! units of N: kg/m2
! N/wl_before -> kg N/kg H2O
where(wl_before>1.0e-4)
  passive_ammonium_uptake = min(soil%soil_organic_matter(:)%ammonium,max(0.0,uptake*soil%soil_organic_matter(:)%ammonium/wl_before))
  passive_nitrate_uptake = min(soil%soil_organic_matter(:)%nitrate,max(0.0,uptake*soil%soil_organic_matter(:)%nitrate/wl_before))
elsewhere
  passive_ammonium_uptake=0.0
  passive_nitrate_uptake=0.0
endwhere

  soil%soil_organic_matter(:)%ammonium=soil%soil_organic_matter(:)%ammonium-passive_ammonium_uptake
  soil%soil_organic_matter(:)%nitrate=soil%soil_organic_matter(:)%nitrate-passive_nitrate_uptake



soil%passive_N_uptake=sum(passive_ammonium_uptake+passive_nitrate_uptake)

!New version that combines the two leaching steps and should do a better job of moving DOC from litter layer
!Note: fine wood litter currently not included, just because we haven't implemented anything with it anywhere
!ZMS Edited to allow for tiled fluxes. Also pass in water content before Richards.
   if (gw_option == GW_TILED) then
      call tracer_leaching_with_litter(soil=soil%soil_organic_matter(:),&
                                        wl=wl_before,&
                                        leaflitter=soil%leafLitter,&
                                        woodlitter=soil%coarsewoodLitter,&
                                        flow=flow, &
                                        litterflow=max(0.0,flow(1)),&
                                        div=div,&
                                        dz=dz(1:num_l),&
                                        dt=delta_time,&
                                        del_soil_DOC=DOC_leached,del_soil_DON=DON_leached,&
                                        del_leaflitter_DOC=leaflitter_DOC_loss,del_leaflitter_DON=leaflitter_DON_loss,&
                                        del_woodlitter_DOC=woodlitter_DOC_loss,del_woodlitter_DON=woodlitter_DON_loss,&
                                        div_DOC_loss=div_DOC_loss,div_DON_loss=div_DON_loss,&
                                        del_soil_NH4=NH4_leached,del_soil_NO3=NO3_leached,&
                                        del_leaflitter_NO3=leaflitter_NO3_loss,del_woodlitter_NO3=woodlitter_NO3_loss,&
                                        del_leaflitter_NH4=leaflitter_NH4_loss,del_woodlitter_NH4=woodlitter_NH4_loss,&
                                        div_NO3_loss=div_NO3_loss,div_NH4_loss=div_NH4_loss,&
                                        tiled=.TRUE.,&
                                        div_hlsp_DOC=soil%div_hlsp_DOC,div_hlsp_DON=soil%div_hlsp_DON,&
                                        div_hlsp_NO3=soil%div_hlsp_NO3,div_hlsp_NH4=soil%div_hlsp_NH4,&
                                        surf_DOC_loss=surf_DOC_loss,surf_DON_loss=surf_DON_loss,&
                                        surf_NO3_loss=surf_NO3_loss,surf_NH4_loss=surf_NH4_loss)
   else
    call tracer_leaching_with_litter(soil=soil%soil_organic_matter(:),&
                                      wl=wl_before,&
                                      leaflitter=soil%leafLitter,&
                                      woodlitter=soil%coarsewoodLitter,&
                                      flow=flow, &
                                      litterflow=max(0.0,flow(1)),&
                                      div=div,&
                                      dz=dz(1:num_l),&
                                      dt=delta_time,&
                                      del_soil_DOC=DOC_leached,del_soil_DON=DON_leached,&
                                      del_leaflitter_DOC=leaflitter_DOC_loss,del_leaflitter_DON=leaflitter_DON_loss,&
                                      del_woodlitter_DOC=woodlitter_DOC_loss,del_woodlitter_DON=woodlitter_DON_loss,&
                                      div_DOC_loss=div_DOC_loss,div_DON_loss=div_DON_loss,&
                                      del_soil_NH4=NH4_leached,del_soil_NO3=NO3_leached,&
                                      del_leaflitter_NO3=leaflitter_NO3_loss,del_woodlitter_NO3=woodlitter_NO3_loss,&
                                      del_leaflitter_NH4=leaflitter_NH4_loss,del_woodlitter_NH4=woodlitter_NH4_loss,&
                                      div_NO3_loss=div_NO3_loss,div_NH4_loss=div_NH4_loss,&
                                      tiled=.FALSE. )

      surf_DOC_loss(:) = 0.
      surf_DON_loss(:) = 0.
      surf_NO3_loss = 0.
      surf_NH4_loss = 0.
   end if

   soil%fast_DOC_leached=soil%fast_DOC_leached+sum(div_DOC_loss(1,:)) + surf_DOC_loss(1)
   soil%slow_DOC_leached=soil%slow_DOC_leached+sum(div_DOC_loss(2,:)) + surf_DOC_loss(2)
   soil%deadmic_DOC_leached=soil%deadmic_DOC_leached+sum(div_DOC_loss(3,:)) + surf_DOC_loss(3)

   soil%fast_DON_leached=soil%fast_DON_leached+sum(div_DON_loss(1,:)) + surf_DON_loss(1)
   soil%slow_DON_leached=soil%slow_DON_leached+sum(div_DON_loss(2,:)) + surf_DON_loss(2)
   soil%deadmic_DON_leached=soil%deadmic_DON_leached+sum(div_DON_loss(3,:)) + surf_DON_loss(3)

   soil%NO3_leached=soil%NO3_leached+sum(div_NO3_loss(:)) + surf_NO3_loss
   soil%NH4_leached=soil%NH4_leached+sum(div_NH4_loss(:)) + surf_NH4_loss

   ! Diagnostic. Later pass this back to land_model for transfer to rivers.
   total_DOC_div = sum(surf_DOC_loss(:))
   total_DON_div = sum(surf_DON_loss(:))
   total_NO3_div = surf_NO3_loss + sum(div_NO3_loss(1:num_l))
   total_NH4_div = surf_NH4_loss + sum(div_NH4_loss(1:num_l))

   do l=1,num_l
      total_DOC_div = total_DOC_div + sum(div_DOC_loss(:,l))
      total_DON_div = total_DON_div + sum(div_DON_loss(:,l))
   end do
   total_DOC_div = total_DOC_div/delta_time
   total_DON_div = total_DON_div/delta_time
   total_NO3_div = total_NO3_div/delta_time
   total_NH4_div = total_NH4_div/delta_time
   soil_doc_runf = total_DOC_div
   soil_don_runf = total_DON_div
   soil_NO3_runf = total_NO3_div
   soil_NH4_runf = total_NH4_div

   if (is_watch_point()) then
      write(*,*)'##### soil_step_2 checkpoint 7 #####'
      call debug_pool(soil%leafLitter,       'leafLitter')
      call debug_pool(soil%fineWoodLitter,   'fineWoodLitter')
      call debug_pool(soil%coarseWoodLitter, 'coarseWoodLitter')
      __DEBUG3__(leaflitter_DOC_loss,woodlitter_DOC_loss,total_DOC_div)
      __DEBUG3__(leaflitter_DON_loss,woodlitter_DON_loss,total_DON_div)
      __DEBUG3__(leaflitter_NO3_loss,woodlitter_NO3_loss,total_NO3_div)
      __DEBUG3__(leaflitter_NH4_loss,woodlitter_NH4_loss,total_NH4_div)
      do l = 1, num_l
         call debug_pool(soil%soil_organic_matter(l), 'soil_organic_matter(l)')
      enddo
   endif

! ----------------------------------------------------------------------------
! given solution for surface energy balance, write diagnostic output.
!

  ! ---- diagnostic section
  call send_tile_data(id_temp, soil%T, diag)
  if (id_lwc > 0) call send_tile_data(id_lwc,  soil%wl/dz(1:num_l), diag)
  if (id_swc > 0) call send_tile_data(id_swc,  soil%ws/dz(1:num_l), diag)
  if (id_psi > 0) call send_tile_data(id_psi,  psi+dPsi, diag)
  ! ZMS uncomment for back-compatibility with diag tables
  if (gw_option == GW_TILED) then
     call send_tile_data(id_deficit, deficit, diag)
     call send_tile_data(id_sat_depth, depth_to_wt_3, diag)
     call send_tile_data(id_sat_dept2, depth_to_wt2_3, diag)
     call send_tile_data(id_z_cap, depth_to_cf_3, diag)
     if (depth_to_wt_2a .ge. -0.5) &
          call send_tile_data(id_sat_depth, depth_to_wt_2a, diag)

  end if
  !f1p: save sat_area_frac for use in tracer deposition calculations
  soil%sat_area_frac = sat_area_frac

  call send_tile_data(id_cf_1, depth_to_cf_1, diag)
  call send_tile_data(id_cf_3, depth_to_cf_3, diag)
  call send_tile_data(id_wt_1, depth_to_wt_1, diag)
  call send_tile_data(id_wt_2, depth_to_wt_2, diag)
  call send_tile_data(id_wt_2a, depth_to_wt_2a, diag)
  call send_tile_data(id_wt_2b, depth_to_wt_2b, diag)
  call send_tile_data(id_wt_3, depth_to_wt_3, diag)
  call send_tile_data(id_wt2_3, depth_to_wt2_3, diag)
  call send_tile_data(id_wt_4, depth_to_wt_4, diag)
  call send_tile_data(id_deficit_2, deficit_2, diag)
  call send_tile_data(id_deficit_3, deficit_3, diag)
  call send_tile_data(id_deficit_4, deficit_4, diag)
  call send_tile_data(id_sat_frac, sat_area_frac, diag)
  call send_tile_data(id_div_bf, div_bf, diag)
  call send_tile_data(id_div_if, div_if, diag)
  call send_tile_data(id_div_al, div_al, diag)

  call send_tile_data(id_ie,   lrunf_ie, diag)
  call send_tile_data(id_sn,   lrunf_sn, diag)
  call send_tile_data(id_bf,   lrunf_bf, diag)
  call send_tile_data(id_if,   lrunf_if, diag)
  call send_tile_data(id_al,   lrunf_al, diag)
  call send_tile_data(id_nu,   lrunf_nu, diag)
  call send_tile_data(id_sc,   lrunf_sc, diag)
  call send_tile_data(id_hie,  hlrunf_ie, diag)
  call send_tile_data(id_hsn,  hlrunf_sn, diag)
  call send_tile_data(id_hbf,  hlrunf_bf, diag)
  call send_tile_data(id_hif,  hlrunf_if, diag)
  call send_tile_data(id_hal,  hlrunf_al, diag)
  call send_tile_data(id_hnu,  hlrunf_nu, diag)
  call send_tile_data(id_hsc,  hlrunf_sc, diag)
  if (id_evap > 0) call send_tile_data(id_evap,  soil_levap+soil_fevap, diag)

   if (id_leaflitter_fast_C_leaching > 0) call send_tile_data(id_leaflitter_fast_C_leaching,leaflitter_DOC_loss(1)/delta_time,diag)
   if (id_leaflitter_slow_C_leaching > 0) call send_tile_data(id_leaflitter_slow_C_leaching,leaflitter_DOC_loss(2)/delta_time,diag)
   if (id_leaflitter_deadmic_C_leaching > 0) call send_tile_data(id_leaflitter_deadmic_C_leaching,leaflitter_DOC_loss(3)/delta_time,diag)

   if (id_leaflitter_fast_N_leaching > 0) call send_tile_data(id_leaflitter_fast_N_leaching,leaflitter_DON_loss(1)/delta_time,diag)
   if (id_leaflitter_slow_N_leaching > 0) call send_tile_data(id_leaflitter_slow_N_leaching,leaflitter_DON_loss(2)/delta_time,diag)
   if (id_leaflitter_deadmic_N_leaching > 0) call send_tile_data(id_leaflitter_deadmic_N_leaching,leaflitter_DON_loss(3)/delta_time,diag)
   if (id_leaflitter_NO3_leaching > 0) call send_tile_data(id_leaflitter_NO3_leaching,leaflitter_NO3_loss/delta_time,diag)
   if (id_leaflitter_NH4_leaching > 0) call send_tile_data(id_leaflitter_NH4_leaching,leaflitter_NH4_loss/delta_time,diag)

   if (id_coarsewoodlitter_fast_C_leaching > 0) call send_tile_data(id_coarsewoodlitter_fast_C_leaching,woodlitter_DOC_loss(1)/delta_time,diag)
   if (id_coarsewoodlitter_slow_C_leaching > 0) call send_tile_data(id_coarsewoodlitter_slow_C_leaching,woodlitter_DOC_loss(2)/delta_time,diag)
   if (id_coarsewoodlitter_deadmic_C_leaching > 0) call send_tile_data(id_coarsewoodlitter_deadmic_C_leaching,woodlitter_DOC_loss(3)/delta_time,diag)

   if (id_coarsewoodlitter_fast_N_leaching > 0) call send_tile_data(id_coarsewoodlitter_fast_N_leaching,woodlitter_DON_loss(1)/delta_time,diag)
   if (id_coarsewoodlitter_slow_N_leaching > 0) call send_tile_data(id_coarsewoodlitter_slow_N_leaching,woodlitter_DON_loss(2)/delta_time,diag)
   if (id_coarsewoodlitter_deadmic_N_leaching > 0) call send_tile_data(id_coarsewoodlitter_deadmic_N_leaching,woodlitter_DON_loss(3)/delta_time,diag)
   if (id_coarsewoodlitter_NO3_leaching > 0) call send_tile_data(id_coarsewoodlitter_NO3_leaching,woodlitter_NO3_loss/delta_time,diag)
   if (id_coarsewoodlitter_NH4_leaching > 0) call send_tile_data(id_coarsewoodlitter_NH4_leaching,woodlitter_NH4_loss/delta_time,diag)

   !Skipping fine wood litter since it's not set up

   do l=1,num_l
        rhizosphere_frac(l)=min(3.141592*((r_rhiz+r_r)**2-r_r**2)*vrl(l),1.0)
   enddo
  if (id_rhizosphere_frac > 0) call send_tile_data(id_rhizosphere_frac, rhizosphere_frac, diag)

  call send_tile_data(id_heat_cap, soil%heat_capacity_dry, diag)
  call send_tile_data(id_active_layer, active_layer_thickness, diag)
  if (gw_option == GW_TILED) then
     call send_tile_data(id_surface_water, surface_water, diag)
     call send_tile_data(id_inun_frac, inundated_frac, diag)
     call send_tile_data(id_wet_frac, wet_frac, diag)
     if (simple_inundation) call send_tile_data(id_rsn_frac, sat_runf_frac, diag)
  end if
  if (id_flow > 0) then
     flow_s(:) = flow(1:num_l) / delta_time
     call send_tile_data(id_flow, flow_s, diag)
  end if
  if (id_reflux > 0) then
     reflux = max(lrunf_ie - lprec_eff, 0.)
     call send_tile_data(id_reflux, reflux, diag)
  end if
  call send_tile_data(id_macro_infilt, flow_macro, diag)
  call send_tile_data(id_fast_C_leaching, DOC_leached(1,:)/delta_time,diag)
  call send_tile_data(id_slow_C_leaching, DOC_leached(2,:)/delta_time,diag)
  call send_tile_data(id_deadmic_C_leaching, DOC_leached(3,:)/delta_time,diag)
  call send_tile_data(id_fast_N_leaching, DON_leached(1,:)/delta_time,diag)
  call send_tile_data(id_slow_N_leaching, DON_leached(2,:)/delta_time,diag)
  call send_tile_data(id_deadmic_N_leaching, DON_leached(3,:)/delta_time,diag)
  call send_tile_data(id_NO3_leaching, NO3_leached(:)/delta_time,diag)
  call send_tile_data(id_NH4_leaching, NH4_leached(:)/delta_time,diag)
  do l=1,num_l
     total_C_leaching(l) = sum(DOC_leached(:,l))/delta_time
     total_DON_leaching(l) = sum(DON_leached(:,l))/delta_time
  end do
  call send_tile_data(id_total_C_leaching, total_C_leaching, diag)
  call send_tile_data(id_total_ON_leaching, total_DON_leaching, diag)
  call send_tile_data(id_total_NO3_leaching, sum(NO3_leached(1:num_l)),diag)
  call send_tile_data(id_total_NH4_leaching, sum(NH4_leached(1:num_l)),diag)
  if (gw_option == GW_TILED) then
     call send_tile_data(id_surf_DOC_loss, sum(surf_DOC_loss(:))/delta_time,diag)
     call send_tile_data(id_surf_DON_loss, sum(surf_DON_loss(:))/delta_time,diag)
     call send_tile_data(id_surf_NO3_loss, surf_NO3_loss/delta_time,diag)
     call send_tile_data(id_surf_NH4_loss, surf_NH4_loss/delta_time,diag)
  end if
  call send_tile_data(id_total_DOC_div_loss, total_DOC_div, diag)
  call send_tile_data(id_total_DON_div_loss, total_DON_div, diag)
  call send_tile_data(id_total_NO3_div_loss, total_NO3_div, diag)
  call send_tile_data(id_total_NH4_div_loss, total_NH4_div, diag)

  if (.not. LM2) call send_tile_data(id_psi_bot, soil%psi(num_l), diag)

end subroutine soil_step_2

! ============================================================================
subroutine soil_step_3(soil, diag)
  type(soil_tile_type), intent(in) :: soil
  type(diag_buff_type), intent(inout) :: diag

  real :: sum_fsc, sum_ssc, sum_deadmic_C, sum_livemic_C, sum_protectedC !, slomtot
  real :: sum_fsn, sum_ssn, sum_deadmic_N, sum_livemic_N, sum_protectedN
  real :: fast_C(num_l), slow_C(num_l), deadMicrobeC(num_l), liveMicrobeC(num_l), protectedC(num_l)
  real :: fast_N(num_l), slow_N(num_l), deadMicrobeN(num_l), liveMicrobeN(num_l), protectedN(num_l)
  real :: litter_fast_C, litter_slow_C, litter_deadmic_C, litter_livemic_C, litter_total_C
  real :: litter_fast_N, litter_slow_N, litter_deadmic_N, litter_livemic_N, litter_total_N
  integer :: layer, ncohorts(num_l), litter_ncohorts
  real, dimension(num_l) :: fast_dissolved_C,slow_dissolved_C,deadmic_dissolved_C
  real, dimension(num_l) :: fast_dissolved_N,slow_dissolved_N,deadmic_dissolved_N
  real :: total_fast_C, total_slow_C, total_deadmic_C, total_livemic_C, total_protected_C, total_dissolved_C, total_carbon
  real :: total_fast_N, total_slow_N, total_deadmic_N, total_livemic_N, total_protected_N, total_dissolved_N, total_nitrogen
  real :: total_carbon_layered(num_l),total_nitrogen_layered(num_l)
  real :: total_nitrate, total_ammonium

  total_carbon_layered=0.0
  total_fast_C=0.0
  total_slow_C=0.0
  total_deadmic_C=0.0
  total_livemic_C=0.0
  total_protected_C=0.0
  total_dissolved_C=0.0
  total_carbon=0.0

  total_nitrogen_layered=0.0
  total_fast_N=0.0
  total_slow_N=0.0
  total_deadmic_N=0.0
  total_livemic_N=0.0
  total_protected_N=0.0
  total_dissolved_N=0.0
  total_nitrogen=0.0

  DO layer=1,num_l
    call poolTotals(soil%soil_organic_matter(layer),fastC=fast_C(layer),fastN=fast_N(layer),slowC=slow_C(layer),slowN=slow_N(layer),&
    deadMicrobeC=deadMicrobeC(layer),liveMicrobeC=liveMicrobeC(layer),protectedC=protectedC(layer),&
    deadMicrobeN=deadMicrobeN(layer),liveMicrobeN=liveMicrobeN(layer),protectedN=protectedN(layer),&
    fast_dissolvedC=fast_dissolved_C(layer),slow_dissolvedC=slow_dissolved_C(layer),&
    fast_dissolvedN=fast_dissolved_N(layer),slow_dissolvedN=slow_dissolved_N(layer),&
    deadmic_dissolvedC=deadmic_dissolved_C(layer),ncohorts=ncohorts(layer),totalCarbon=total_carbon_layered(layer),&
    deadmic_dissolvedN=deadmic_dissolved_N(layer),totalNitrogen=total_nitrogen_layered(layer)  )
  ENDDO

  total_fast_C=sum(fast_C)
  total_slow_C=sum(slow_C)
  total_deadmic_C=sum(deadMicrobeC)
  total_livemic_C=sum(liveMicrobeC)
  total_protected_C=sum(protectedC)
  total_dissolved_C=sum(fast_dissolved_C+slow_dissolved_C+deadmic_dissolved_C)

  total_fast_N=sum(fast_N)
  total_slow_N=sum(slow_N)
  total_deadmic_N=sum(deadMicrobeN)
  total_livemic_N=sum(liveMicrobeN)
  total_protected_N=sum(protectedN)
  total_dissolved_N=sum(fast_dissolved_N+slow_dissolved_N+deadmic_dissolved_N)
  total_nitrate = sum(soil%soil_organic_matter(:)%nitrate)
  total_ammonium = sum(soil%soil_organic_matter(:)%ammonium)

  if (id_fast_soil_C > 0) call send_tile_data(id_fast_soil_C, fast_C/dz, diag)
  if (id_slow_soil_C > 0) call send_tile_data(id_slow_soil_C, slow_C/dz, diag)
  if (id_deadmic_C > 0) call send_tile_data(id_deadmic_C, deadMicrobeC/dz, diag)
  if (id_livemic_C > 0) call send_tile_data(id_livemic_C, liveMicrobeC/dz, diag)
  if (id_protectedC > 0) call send_tile_data(id_protectedC, protectedC/dz, diag)
  if (id_nsoilcohorts > 0) call send_tile_data(id_nsoilcohorts, real(ncohorts), diag)
  if (id_fast_dissolved_C > 0) call send_tile_data(id_fast_dissolved_C, fast_dissolved_C/dz, diag)
  if (id_slow_dissolved_C > 0) call send_tile_data(id_slow_dissolved_C, slow_dissolved_C/dz, diag)
  if (id_deadmic_dissolved_C > 0) call send_tile_data(id_deadmic_dissolved_C, deadmic_dissolved_C/dz, diag)
  if (id_total_carbon_layered > 0) call send_tile_data(id_total_carbon_layered, total_carbon_layered/dz,diag)

  if (id_fast_DOC_div_loss > 0) call send_tile_data(id_fast_DOC_div_loss, soil%fast_DOC_leached,diag)
  if (id_slow_DOC_div_loss > 0) call send_tile_data(id_slow_DOC_div_loss, soil%slow_DOC_leached,diag)
  if (id_deadmic_DOC_div_loss > 0) call send_tile_data(id_deadmic_DOC_div_loss, soil%deadmic_DOC_leached,diag)


  if (id_fast_soil_N > 0) call send_tile_data(id_fast_soil_N, fast_N/dz, diag)
  if (id_slow_soil_N > 0) call send_tile_data(id_slow_soil_N, slow_N/dz, diag)
  if (id_deadmic_N > 0) call send_tile_data(id_deadmic_N, deadMicrobeN/dz, diag)
  if (id_livemic_N > 0) call send_tile_data(id_livemic_N, liveMicrobeN/dz, diag)
  if (id_protectedN > 0) call send_tile_data(id_protectedN, protectedN/dz, diag)
  if (id_fast_dissolved_N > 0) call send_tile_data(id_fast_dissolved_N, fast_dissolved_N/dz, diag)
  if (id_slow_dissolved_N > 0) call send_tile_data(id_slow_dissolved_N, slow_dissolved_N/dz, diag)
  if (id_deadmic_dissolved_N > 0) call send_tile_data(id_deadmic_dissolved_N, deadmic_dissolved_N/dz, diag)
  if (id_total_nitrogen_layered > 0) call send_tile_data(id_total_nitrogen_layered, total_nitrogen_layered/dz,diag)
  if (id_soil_ammonium > 0) call send_tile_data(id_soil_ammonium, soil%soil_organic_matter(1:num_l)%ammonium/dz,diag)
  if (id_soil_nitrate > 0) call send_tile_data(id_soil_nitrate, soil%soil_organic_matter(1:num_l)%nitrate/dz,diag)

  if (id_fast_DON_div_loss > 0) call send_tile_data(id_fast_DON_div_loss, soil%fast_DON_leached,diag)
  if (id_slow_DON_div_loss > 0) call send_tile_data(id_slow_DON_div_loss, soil%slow_DON_leached,diag)
  if (id_deadmic_DON_div_loss > 0) call send_tile_data(id_deadmic_DON_div_loss, soil%deadmic_DON_leached,diag)
  if (id_NO3_div_loss > 0) call send_tile_data(id_NO3_div_loss, soil%NO3_leached,diag)
  if (id_NH4_div_loss > 0) call send_tile_data(id_NH4_div_loss, soil%NH4_leached,diag)

  call poolTotals(soil%leafLitter,fastC=litter_fast_C,slowC=litter_slow_C,deadMicrobeC=litter_deadmic_C,&
                    liveMicrobeC=litter_livemic_C,ncohorts=litter_ncohorts,&
                    fastN=litter_fast_N,slowN=litter_slow_N,deadMicrobeN=litter_deadmic_N,liveMicrobeN=litter_livemic_N,&
                    totalCarbon=litter_total_C, totalNitrogen=litter_total_N)
  total_fast_C=total_fast_C+litter_fast_C
  total_slow_C=total_slow_C+litter_slow_C
  total_deadmic_C=total_deadmic_C+litter_deadmic_C
  total_livemic_C=total_livemic_C+litter_livemic_C
  total_dissolved_C=total_dissolved_C+sum(soil%leafLitter%dissolved_carbon(:))

  total_fast_N=total_fast_N+litter_fast_N
  total_slow_N=total_slow_N+litter_slow_N
  total_deadmic_N=total_deadmic_N+litter_deadmic_N
  total_livemic_N=total_livemic_N+litter_livemic_N
  total_dissolved_N=total_dissolved_N+sum(soil%leafLitter%dissolved_nitrogen(:))
  total_nitrate=total_nitrate+soil%leafLitter%nitrate
  total_ammonium=total_ammonium+soil%leafLitter%ammonium


  if (id_nleaflittercohorts > 0) call send_tile_data(id_nleaflittercohorts, real(litter_ncohorts), diag)
  if (id_leaflitter_fast_C > 0) call send_tile_data(id_leaflitter_fast_C, litter_fast_C, diag)
  if (id_leaflitter_slow_C > 0) call send_tile_data(id_leaflitter_slow_C, litter_slow_C, diag)
  if (id_leaflitter_deadmic_C > 0) call send_tile_data(id_leaflitter_deadmic_C, litter_deadmic_C, diag)
  if (id_leaflitter_livemic_C > 0) call send_tile_data(id_leaflitter_livemic_C, litter_livemic_C, diag)
  if (id_leaflitter_fast_dissolved_C > 0) call send_tile_data(id_leaflitter_fast_dissolved_C, soil%leafLitter%dissolved_carbon(1), diag)
  if (id_leaflitter_slow_dissolved_C > 0) call send_tile_data(id_leaflitter_slow_dissolved_C, soil%leafLitter%dissolved_carbon(2), diag)
  if (id_leaflitter_deadmic_dissolved_C > 0) call send_tile_data(id_leaflitter_deadmic_dissolved_C, soil%leafLitter%dissolved_carbon(3), diag)
  if (id_leaflitter_total_C > 0) call send_tile_data(id_leaflitter_total_C, litter_total_C,diag)

  if (id_leaflitter_fast_N > 0) call send_tile_data(id_leaflitter_fast_N, litter_fast_N, diag)
  if (id_leaflitter_slow_N > 0) call send_tile_data(id_leaflitter_slow_N, litter_slow_N, diag)
  if (id_leaflitter_deadmic_N > 0) call send_tile_data(id_leaflitter_deadmic_N, litter_deadmic_N, diag)
  if (id_leaflitter_livemic_N > 0) call send_tile_data(id_leaflitter_livemic_N, litter_livemic_N, diag)
  if (id_leaflitter_fast_dissolved_N > 0) call send_tile_data(id_leaflitter_fast_dissolved_N, soil%leafLitter%dissolved_nitrogen(1), diag)
  if (id_leaflitter_slow_dissolved_N > 0) call send_tile_data(id_leaflitter_slow_dissolved_N, soil%leafLitter%dissolved_nitrogen(2), diag)
  if (id_leaflitter_deadmic_dissolved_N > 0) call send_tile_data(id_leaflitter_deadmic_dissolved_N, soil%leafLitter%dissolved_nitrogen(3), diag)
  if (id_leaflitter_total_N > 0) call send_tile_data(id_leaflitter_total_N, litter_total_N,diag)


  call poolTotals(soil%fineWoodLitter,fastC=litter_fast_C,slowC=litter_slow_C,deadMicrobeC=litter_deadmic_C,&
                    liveMicrobeC=litter_livemic_C,ncohorts=litter_ncohorts,&
                    fastN=litter_fast_N,slowN=litter_slow_N,deadMicrobeN=litter_deadmic_N,liveMicrobeN=litter_livemic_N,&
                    totalCarbon=litter_total_C, totalNitrogen=litter_total_N)
  total_fast_C=total_fast_C+litter_fast_C
  total_slow_C=total_slow_C+litter_slow_C
  total_deadmic_C=total_deadmic_C+litter_deadmic_C
  total_livemic_C=total_livemic_C+litter_livemic_C
  total_dissolved_C=total_dissolved_C+sum(soil%fineWoodLitter%dissolved_carbon(:))

  total_fast_N=total_fast_N+litter_fast_N
  total_slow_N=total_slow_N+litter_slow_N
  total_deadmic_N=total_deadmic_N+litter_deadmic_N
  total_livemic_N=total_livemic_N+litter_livemic_N
  total_dissolved_N=total_dissolved_N+sum(soil%fineWoodLitter%dissolved_nitrogen(:))
  total_nitrate=total_nitrate+soil%fineWoodLitter%nitrate
  total_ammonium=total_ammonium+soil%fineWoodLitter%ammonium


  if (id_nfineWoodlittercohorts > 0) call send_tile_data(id_nfineWoodlittercohorts, real(litter_ncohorts), diag)
  if (id_fineWoodlitter_fast_C > 0) call send_tile_data(id_fineWoodlitter_fast_C, litter_fast_C, diag)
  if (id_fineWoodlitter_slow_C > 0) call send_tile_data(id_fineWoodlitter_slow_C, litter_slow_C, diag)
  if (id_fineWoodlitter_deadmic_C > 0) call send_tile_data(id_fineWoodlitter_deadmic_C, litter_deadmic_C, diag)
  if (id_fineWoodlitter_livemic_C > 0) call send_tile_data(id_fineWoodlitter_livemic_C, litter_livemic_C, diag)
  if (id_fineWoodlitter_fast_dissolved_C > 0) call send_tile_data(id_fineWoodlitter_fast_dissolved_C, soil%fineWoodLitter%dissolved_carbon(1), diag)
  if (id_fineWoodlitter_slow_dissolved_C > 0) call send_tile_data(id_fineWoodlitter_slow_dissolved_C, soil%fineWoodLitter%dissolved_carbon(2), diag)
  if (id_fineWoodlitter_deadmic_dissolved_C > 0) call send_tile_data(id_fineWoodlitter_deadmic_dissolved_C, soil%fineWoodLitter%dissolved_carbon(3), diag)
  if (id_fineWoodlitter_total_C > 0) call send_tile_data(id_fineWoodlitter_total_C, litter_total_C,diag)

  if (id_fineWoodlitter_fast_N > 0) call send_tile_data(id_fineWoodlitter_fast_N, litter_fast_N, diag)
  if (id_fineWoodlitter_slow_N > 0) call send_tile_data(id_fineWoodlitter_slow_N, litter_slow_N, diag)
  if (id_fineWoodlitter_deadmic_N > 0) call send_tile_data(id_fineWoodlitter_deadmic_N, litter_deadmic_N, diag)
  if (id_fineWoodlitter_livemic_N > 0) call send_tile_data(id_fineWoodlitter_livemic_N, litter_livemic_N, diag)
  if (id_fineWoodlitter_fast_dissolved_N > 0) call send_tile_data(id_fineWoodlitter_fast_dissolved_N, soil%fineWoodLitter%dissolved_nitrogen(1), diag)
  if (id_fineWoodlitter_slow_dissolved_N > 0) call send_tile_data(id_fineWoodlitter_slow_dissolved_N, soil%fineWoodLitter%dissolved_nitrogen(2), diag)
  if (id_fineWoodlitter_deadmic_dissolved_N > 0) call send_tile_data(id_fineWoodlitter_deadmic_dissolved_N, soil%fineWoodLitter%dissolved_nitrogen(3), diag)
  if (id_fineWoodlitter_total_N > 0) call send_tile_data(id_fineWoodlitter_total_N, litter_total_N,diag)


  call poolTotals(soil%coarseWoodLitter,fastC=litter_fast_C,slowC=litter_slow_C,deadMicrobeC=litter_deadmic_C,&
                    liveMicrobeC=litter_livemic_C,ncohorts=litter_ncohorts,&
                    fastN=litter_fast_N,slowN=litter_slow_N,deadMicrobeN=litter_deadmic_N,liveMicrobeN=litter_livemic_N,&
                    totalCarbon=litter_total_C, totalNitrogen=litter_total_N)
  total_fast_C=total_fast_C+litter_fast_C
  total_slow_C=total_slow_C+litter_slow_C
  total_deadmic_C=total_deadmic_C+litter_deadmic_C
  total_livemic_C=total_livemic_C+litter_livemic_C
  total_dissolved_C=total_dissolved_C+sum(soil%coarseWoodLitter%dissolved_carbon(:))

  total_fast_N=total_fast_N+litter_fast_N
  total_slow_N=total_slow_N+litter_slow_N
  total_deadmic_N=total_deadmic_N+litter_deadmic_N
  total_livemic_N=total_livemic_N+litter_livemic_N
  total_dissolved_N=total_dissolved_N+sum(soil%coarseWoodLitter%dissolved_nitrogen(:))
  total_nitrate=total_nitrate+soil%coarseWoodLitter%nitrate
  total_ammonium=total_ammonium+soil%coarseWoodLitter%ammonium


  if (id_ncoarseWoodlittercohorts > 0) call send_tile_data(id_ncoarseWoodlittercohorts, real(litter_ncohorts), diag)
  if (id_coarseWoodlitter_fast_C > 0) call send_tile_data(id_coarseWoodlitter_fast_C, litter_fast_C, diag)
  if (id_coarseWoodlitter_slow_C > 0) call send_tile_data(id_coarseWoodlitter_slow_C, litter_slow_C, diag)
  if (id_coarseWoodlitter_deadmic_C > 0) call send_tile_data(id_coarseWoodlitter_deadmic_C, litter_deadmic_C, diag)
  if (id_coarseWoodlitter_livemic_C > 0) call send_tile_data(id_coarseWoodlitter_livemic_C, litter_livemic_C, diag)
  if (id_coarseWoodlitter_fast_dissolved_C > 0) call send_tile_data(id_coarseWoodlitter_fast_dissolved_C, soil%coarseWoodLitter%dissolved_carbon(1), diag)
  if (id_coarseWoodlitter_slow_dissolved_C > 0) call send_tile_data(id_coarseWoodlitter_slow_dissolved_C, soil%coarseWoodLitter%dissolved_carbon(2), diag)
  if (id_coarseWoodlitter_deadmic_dissolved_C > 0) call send_tile_data(id_coarseWoodlitter_deadmic_dissolved_C, soil%coarseWoodLitter%dissolved_carbon(3), diag)
  if (id_coarseWoodlitter_total_C > 0) call send_tile_data(id_coarseWoodlitter_total_C, litter_total_C,diag)

  if (id_coarseWoodlitter_fast_N > 0) call send_tile_data(id_coarseWoodlitter_fast_N, litter_fast_N, diag)
  if (id_coarseWoodlitter_slow_N > 0) call send_tile_data(id_coarseWoodlitter_slow_N, litter_slow_N, diag)
  if (id_coarseWoodlitter_deadmic_N > 0) call send_tile_data(id_coarseWoodlitter_deadmic_N, litter_deadmic_N, diag)
  if (id_coarseWoodlitter_livemic_N > 0) call send_tile_data(id_coarseWoodlitter_livemic_N, litter_livemic_N, diag)
  if (id_coarseWoodlitter_fast_dissolved_N > 0) call send_tile_data(id_coarseWoodlitter_fast_dissolved_N, soil%coarseWoodLitter%dissolved_nitrogen(1), diag)
  if (id_coarseWoodlitter_slow_dissolved_N > 0) call send_tile_data(id_coarseWoodlitter_slow_dissolved_N, soil%coarseWoodLitter%dissolved_nitrogen(2), diag)
  if (id_coarseWoodlitter_deadmic_dissolved_N > 0) call send_tile_data(id_coarseWoodlitter_deadmic_dissolved_N, soil%coarseWoodLitter%dissolved_nitrogen(3), diag)
  if (id_coarseWoodlitter_total_N > 0) call send_tile_data(id_coarseWoodlitter_total_N, litter_total_N,diag)


  sum_fsc = total_fast_C
  sum_ssc = total_slow_C
  sum_deadmic_C = total_deadmic_C
  sum_livemic_C = total_livemic_C
  sum_protectedC = total_protected_C
  total_carbon=total_fast_C+total_slow_C+total_deadmic_C+total_livemic_C+total_dissolved_C+total_protected_C
  if (id_fsc > 0) call send_tile_data(id_fsc, sum_fsc, diag)
  if (id_ssc > 0) call send_tile_data(id_ssc, sum_ssc, diag)
  if (id_deadmic_C_total > 0) call send_tile_data(id_deadmic_C_total, sum_deadmic_C, diag)
  if (id_livemic_C_total > 0) call send_tile_data(id_livemic_C_total, sum_livemic_C, diag)

  if (id_protected_C_total > 0) call send_tile_data(id_protected_C_total, sum_protectedC, diag)
  if (id_dissolved_C_total > 0) call send_tile_data(id_dissolved_C_total, total_dissolved_C, diag)
  if (id_total_soil_C > 0) call send_tile_data(id_total_soil_C, total_carbon, diag)


  sum_fsn = total_fast_N
  sum_ssn = total_slow_N
  sum_deadmic_N = total_deadmic_N
  sum_livemic_N = total_livemic_N
  sum_protectedN = total_protected_N
  total_nitrogen=total_fast_N+total_slow_N+total_deadmic_N+total_livemic_N+total_dissolved_N+total_protected_N
  if (id_fsn > 0) call send_tile_data(id_fsn, sum_fsn, diag)
  if (id_ssn > 0) call send_tile_data(id_ssn, sum_ssn, diag)
  if (id_deadmic_N_total > 0) call send_tile_data(id_deadmic_N_total, sum_deadmic_N, diag)
  if (id_livemic_N_total > 0) call send_tile_data(id_livemic_N_total, sum_livemic_N, diag)

  if (id_protected_N_total > 0) call send_tile_data(id_protected_N_total, sum_protectedN, diag)
  if (id_dissolved_N_total > 0) call send_tile_data(id_dissolved_N_total, total_dissolved_N, diag)
  if (id_total_soil_organic_N > 0) call send_tile_data(id_total_soil_N, total_nitrogen, diag)
  if (id_total_soil_N > 0) call send_tile_data(id_total_soil_N, total_nitrogen+total_ammonium+total_nitrate, diag)
  if (id_total_soil_ammonium > 0) call send_tile_data(id_total_soil_ammonium, total_ammonium, diag)
  if (id_total_soil_nitrate > 0) call send_tile_data(id_total_soil_nitrate, total_nitrate, diag)


end subroutine soil_step_3


! ============================================================================
subroutine Dsdt(vegn, soil, diag, soilt, theta)
  type(vegn_tile_type), intent(inout) :: vegn
  type(soil_tile_type), intent(inout) :: soil
  type(diag_buff_type), intent(inout) :: diag
  real                , intent(in)    :: soilt ! average soil temperature, deg K
  real                , intent(in)    :: theta ! average soil moisture

  select case (soil_carbon_option)
  case (SOILC_CENTURY, SOILC_CENTURY_BY_LAYER)
     call Dsdt_CENTURY(vegn, soil, diag, soilt, theta)
 case (SOILC_CORPSE, SOILC_CORPSE_N)
     call Dsdt_CORPSE(vegn, soil, diag)
  case default
     call error_mesg('Dsdt','soil_carbon_option is invalid. This should never happen. Contact developer', FATAL)
  end select
end subroutine Dsdt


! ============================================================================
subroutine Dsdt_CORPSE(vegn, soil, diag)
  type(vegn_tile_type), intent(inout) :: vegn
  type(soil_tile_type), intent(inout) :: soil
  type(diag_buff_type), intent(inout) :: diag

  real :: leaflitter_fast_C_loss_rate, leaflitter_slow_C_loss_rate, leaflitter_deadmic_C_loss_rate
  real :: finewoodlitter_fast_C_loss_rate, finewoodlitter_slow_C_loss_rate, finewoodlitter_deadmic_C_loss_rate
  real :: coarsewoodlitter_fast_C_loss_rate, coarsewoodlitter_slow_C_loss_rate, coarsewoodlitter_deadmic_C_loss_rate
  real :: leaflitter_fast_N_loss_rate, leaflitter_slow_N_loss_rate, leaflitter_deadmic_N_loss_rate
  real :: finewoodlitter_fast_N_loss_rate, finewoodlitter_slow_N_loss_rate, finewoodlitter_deadmic_N_loss_rate
  real :: coarsewoodlitter_fast_N_loss_rate, coarsewoodlitter_slow_N_loss_rate, coarsewoodlitter_deadmic_N_loss_rate
  real :: fast_C_loss_rate(size(soil%soil_organic_matter))
  real :: slow_C_loss_rate(size(soil%soil_organic_matter))
  real :: dead_microbe_C_loss_rate(size(soil%soil_organic_matter))
  real :: fast_N_loss_rate(size(soil%soil_organic_matter))
  real :: slow_N_loss_rate(size(soil%soil_organic_matter))
  real :: dead_microbe_N_loss_rate(size(soil%soil_organic_matter))
  real, dimension(size(soil%soil_organic_matter)) :: decomp_T,decomp_theta,ice_porosity
  real :: A          (size(soil%soil_organic_matter)) ! decomp rate reduction due to moisture and temperature

  integer :: badCohort   ! For soil carbon pool carbon balance and invalid number check
  integer :: k
  real :: CO2prod,protected_C_produced(3,size(soil%soil_organic_matter)),protected_C_turnover_rate(3,size(soil%soil_organic_matter))
  real :: leaflitter_protected_C_produced(3),leaflitter_protected_C_turnover_rate(3)
  real :: finewoodlitter_protected_C_produced(3),finewoodlitter_protected_C_turnover_rate(3)
  real :: coarsewoodlitter_protected_C_produced(3),coarsewoodlitter_protected_C_turnover_rate(3)
  real :: protected_N_produced(3,size(soil%soil_organic_matter)),protected_N_turnover_rate(3,size(soil%soil_organic_matter))
  real :: soil_nitrif(size(soil%soil_organic_matter)),soil_denitrif(size(soil%soil_organic_matter))
  real :: leaflitter_protected_N_produced(3),leaflitter_protected_N_turnover_rate(3)
  real :: finewoodlitter_protected_N_produced(3),finewoodlitter_protected_N_turnover_rate(3)
  real :: coarsewoodlitter_protected_N_produced(3),coarsewoodlitter_protected_N_turnover_rate(3)
  real :: leaflitter_C_dissolved(3),leaflitter_C_deposited(3),C_dissolved(3,num_l),C_deposited(3,num_l)
  real :: finewoodlitter_C_dissolved(3),finewoodlitter_C_deposited(3)
  real :: coarsewoodlitter_C_dissolved(3),coarsewoodlitter_C_deposited(3)
  real :: leaflitter_N_dissolved(3),leaflitter_N_deposited(3),N_dissolved(3,num_l),N_deposited(3,num_l)
  real :: finewoodlitter_N_dissolved(3),finewoodlitter_N_deposited(3)
  real :: coarsewoodlitter_N_dissolved(3),coarsewoodlitter_N_deposited(3)
  real :: leaflitter_nitrif, leaflitter_denitrif, finewoodlitter_nitrif, finewoodlitter_denitrif, coarsewoodlitter_nitrif, coarsewoodlitter_denitrif
  real :: deadmic_C_produced(size(soil%soil_organic_matter)), leaflitter_deadmic_C_produced, finewoodlitter_deadmic_C_produced, coarsewoodlitter_deadmic_C_produced
  real :: deadmic_N_produced(size(soil%soil_organic_matter)), leaflitter_deadmic_N_produced, finewoodlitter_deadmic_N_produced, coarsewoodlitter_deadmic_N_produced

  real :: total_fast_C, total_slow_C, total_deadmic_C,total_livemic_C,temp_fast_C,temp_slow_C,temp_livemic_C,temp_deadmic_C,temp_protected_C
  real :: total_fast_N, total_slow_N, total_deadmic_N,total_livemic_N,temp_fast_N,temp_slow_N,temp_livemic_N,temp_deadmic_N,temp_protected_N
  real :: temp_protected_fast_C, temp_protected_slow_C, temp_protected_deadmic_C
  real :: temp_protected_fast_N, temp_protected_slow_N, temp_protected_deadmic_N
  real :: leaflitter_N_mineralization, leaflitter_N_immobilization, finewoodlitter_N_mineralization, finewoodlitter_N_immobilization
  real :: coarsewoodlitter_N_mineralization, coarsewoodlitter_N_immobilization, soil_N_mineralization(size(soil%soil_organic_matter)), soil_N_immobilization(size(soil%soil_organic_matter))
  real,dimension(num_l) :: uptake_frac_max, vegn_uptake_term
  integer :: point_i,point_j,point_k,point_face

  A(:) = A_function(soil%T(:), soil_theta(soil))
  decomp_T = soil%T(:)
  decomp_theta = soil_theta(soil)
  ice_porosity = soil_ice_porosity(soil)

  vegn%rh=0.0  ! Is this the logical place to add total denitrification?
  total_fast_C=0.0
  total_slow_C=0.0
  total_deadmic_C=0.0
  total_livemic_C=0.0
  total_fast_N=0.0
  total_slow_N=0.0
  total_deadmic_N=0.0
  total_livemic_N=0.0

  !  First surface litter is decomposed
  !  update_pool with SOILC_CORPSE should ignore N and set all N outputs to 0.0
  call update_pool(pool=soil%leafLitter,T=decomp_T(1),theta=decomp_theta(1),air_filled_porosity=1.0-(decomp_theta(1)+ice_porosity(1)),&
            liquid_water=soil%wl(1),frozen_water=soil%ws(1),dt=dt_fast_yr,layerThickness=dz(1),&
            fast_C_loss_rate=leaflitter_fast_C_loss_rate,slow_C_loss_rate=leaflitter_slow_C_loss_rate, deadmic_C_loss_rate=leaflitter_deadmic_C_loss_rate, CO2prod=CO2prod, &
            fast_N_loss_rate=leaflitter_fast_N_loss_rate,slow_N_loss_rate=leaflitter_slow_N_loss_rate, deadmic_N_loss_rate=leaflitter_deadmic_N_loss_rate, &
            deadmic_C_produced=leaflitter_deadmic_C_produced, protected_C_produced=leaflitter_protected_C_produced, protected_turnover_rate=leaflitter_protected_C_turnover_rate, &
            deadmic_N_produced=leaflitter_deadmic_N_produced, protected_N_produced=leaflitter_protected_N_produced, protected_N_turnover_rate=leaflitter_protected_N_turnover_rate, &
            C_dissolved=leaflitter_C_dissolved, deposited_C=leaflitter_C_deposited, badCohort=badCohort,&
            N_dissolved=leaflitter_N_dissolved, deposited_N=leaflitter_N_deposited,nitrification=leaflitter_nitrif,denitrification=leaflitter_denitrif,&
            N_mineralization=leaflitter_N_mineralization,N_immobilization=leafLitter_N_immobilization)
  IF (badCohort.ne.0) THEN
        call get_current_point(point_i,point_j,point_k,point_face)
        WRITE (*,*), 'Found bad cohort in leaf litter.  Point i,j,k,face:',point_i,point_j,point_k,point_face
        WRITE (*,*), 'T=',decomp_T(1),'theta=',decomp_theta(1),'dt=',dt_fast_yr
        call error_mesg('Dsdt','Found bad cohort in leaf litter',FATAL)
  ENDIF

  ! loss of C to atmosphere
  vegn%rh=vegn%rh + CO2prod/dt_fast_yr

  call update_pool(pool=soil%fineWoodLitter,T=decomp_T(1),theta=decomp_theta(1),air_filled_porosity=1.0-(decomp_theta(1)+ice_porosity(1)),&
            liquid_water=soil%wl(1),frozen_water=soil%ws(1),dt=dt_fast_yr,layerThickness=dz(1),&
            fast_C_loss_rate=fineWoodlitter_fast_C_loss_rate,slow_C_loss_rate=fineWoodlitter_slow_C_loss_rate, deadmic_C_loss_rate=fineWoodlitter_deadmic_C_loss_rate, CO2prod=CO2prod, &
            fast_N_loss_rate=fineWoodlitter_fast_N_loss_rate,slow_N_loss_rate=fineWoodlitter_slow_N_loss_rate, deadmic_N_loss_rate=fineWoodlitter_deadmic_N_loss_rate, &
            deadmic_C_produced=fineWoodlitter_deadmic_C_produced, protected_C_produced=fineWoodlitter_protected_C_produced, protected_turnover_rate=fineWoodlitter_protected_C_turnover_rate, &
            deadmic_N_produced=fineWoodlitter_deadmic_N_produced, protected_N_produced=fineWoodlitter_protected_N_produced, protected_N_turnover_rate=fineWoodlitter_protected_N_turnover_rate, &
            C_dissolved=fineWoodlitter_C_dissolved, deposited_C=fineWoodlitter_C_deposited, badCohort=badCohort,&
            N_dissolved=fineWoodlitter_N_dissolved, deposited_N=fineWoodlitter_N_deposited,nitrification=fineWoodlitter_nitrif,denitrification=fineWoodlitter_denitrif,&
            N_mineralization=fineWoodlitter_N_mineralization,N_immobilization=fineWoodLitter_N_immobilization)

  IF (badCohort.ne.0) THEN
        call get_current_point(point_i,point_j,point_k,point_face)
        WRITE (*,*), 'Found bad cohort in fineWood litter.  Point i,j,k,face:',point_i,point_j,point_k,point_face
        WRITE (*,*), 'T=',decomp_T(1),'theta=',decomp_theta(1),'dt=',dt_fast_yr
        call error_mesg('Dsdt','Found bad cohort in fineWood litter',FATAL)
  ENDIF

  ! loss of C to atmosphere
  vegn%rh=vegn%rh + CO2prod/dt_fast_yr

  call update_pool(pool=soil%coarseWoodLitter,T=decomp_T(1),theta=decomp_theta(1),air_filled_porosity=1.0-(decomp_theta(1)+ice_porosity(1)),&
            liquid_water=soil%wl(1),frozen_water=soil%ws(1),dt=dt_fast_yr,layerThickness=dz(1),&
            fast_C_loss_rate=coarseWoodlitter_fast_C_loss_rate,slow_C_loss_rate=coarseWoodlitter_slow_C_loss_rate, &
            deadmic_C_loss_rate=coarseWoodlitter_deadmic_C_loss_rate, CO2prod=CO2prod, &
            fast_N_loss_rate=coarseWoodlitter_fast_N_loss_rate,slow_N_loss_rate=coarseWoodlitter_slow_N_loss_rate, &
            deadmic_N_loss_rate=coarseWoodlitter_deadmic_N_loss_rate, &
            deadmic_C_produced=coarseWoodlitter_deadmic_C_produced, protected_C_produced=coarseWoodlitter_protected_C_produced, &
            protected_turnover_rate=coarseWoodlitter_protected_C_turnover_rate, &
            deadmic_N_produced=coarseWoodlitter_deadmic_N_produced, protected_N_produced=coarseWoodlitter_protected_N_produced, &
            protected_N_turnover_rate=coarseWoodlitter_protected_N_turnover_rate, &
            C_dissolved=coarseWoodlitter_C_dissolved, deposited_C=coarseWoodlitter_C_deposited, badCohort=badCohort,&
            N_dissolved=coarseWoodlitter_N_dissolved, deposited_N=coarseWoodlitter_N_deposited,&
            nitrification=coarseWoodlitter_nitrif,denitrification=coarseWoodlitter_denitrif,&
            N_mineralization=coarseWoodlitter_N_mineralization,N_immobilization=coarseWoodLitter_N_immobilization)

  IF (badCohort.ne.0) THEN
        call get_current_point(point_i,point_j,point_k,point_face)
        WRITE (*,*), 'Found bad cohort in coarseWood litter.  Point i,j,k,face:',point_i,point_j,point_k,point_face
        WRITE (*,*), 'T=',decomp_T(1),'theta=',decomp_theta(1),'dt=',dt_fast_yr
        call error_mesg('Dsdt','Found bad cohort in coarseWood litter',FATAL)
  ENDIF

  ! loss of C to atmosphere
  vegn%rh=vegn%rh + CO2prod/dt_fast_yr

  call poolTotals(soil%leafLitter,fastC=temp_fast_C,slowC=temp_slow_C,deadMicrobeC=temp_deadmic_C,liveMicrobeC=temp_livemic_C, &
                                fastN=temp_fast_N,slowN=temp_slow_N,deadMicrobeN=temp_deadmic_N,liveMicrobeN=temp_livemic_N    )
  total_fast_C=total_fast_C+temp_fast_N
  total_slow_C=total_slow_C+temp_slow_N
  total_deadmic_C=total_deadmic_C+temp_deadmic_N
  total_livemic_C=total_livemic_C+temp_livemic_N

  !Accumulate turnover rates for determining steady state pools
  if(temp_fast_C>0)soil%leaflitter_fast_C_turnover_accumulated=soil%leaflitter_fast_C_turnover_accumulated+leaflitter_fast_C_loss_rate/temp_fast_C
  if(temp_slow_C>0)soil%leaflitter_slow_C_turnover_accumulated=soil%leaflitter_slow_C_turnover_accumulated+leaflitter_slow_C_loss_rate/temp_slow_C
  if(temp_deadmic_C>0)soil%leaflitter_deadmic_C_turnover_accumulated=soil%leaflitter_deadmic_C_turnover_accumulated+leaflitter_deadmic_C_loss_rate/temp_deadmic_C
  soil%leaflitter_deadmic_C_in=soil%leaflitter_deadmic_C_in+leaflitter_deadmic_C_produced
  soil%leaflitter_fsc_in=soil%leaflitter_fsc_in+leaflitter_C_deposited(1)-leaflitter_C_dissolved(1)
  soil%leaflitter_ssc_in=soil%leaflitter_ssc_in+leaflitter_C_deposited(2)-leaflitter_C_dissolved(2)
  soil%leaflitter_deadmic_C_in=soil%leaflitter_deadmic_C_in+leaflitter_C_deposited(3)-leaflitter_C_dissolved(3)


  if(temp_fast_N>0)soil%leaflitter_fast_N_turnover_accumulated=soil%leaflitter_fast_N_turnover_accumulated+leaflitter_fast_N_loss_rate/temp_fast_N
  if(temp_slow_N>0)soil%leaflitter_slow_N_turnover_accumulated=soil%leaflitter_slow_N_turnover_accumulated+leaflitter_slow_N_loss_rate/temp_slow_N
  if(temp_deadmic_N>0)soil%leaflitter_deadmic_N_turnover_accumulated=soil%leaflitter_deadmic_N_turnover_accumulated+leaflitter_deadmic_N_loss_rate/temp_deadmic_N
  soil%leaflitter_deadmic_N_in=soil%leaflitter_deadmic_N_in+leaflitter_deadmic_N_produced
  soil%leaflitter_fsn_in=soil%leaflitter_fsn_in+leaflitter_N_deposited(1)-leaflitter_N_dissolved(1)
  soil%leaflitter_ssn_in=soil%leaflitter_ssn_in+leaflitter_N_deposited(2)-leaflitter_N_dissolved(2)
  soil%leaflitter_deadmic_N_in=soil%leaflitter_deadmic_N_in+leaflitter_N_deposited(3)-leaflitter_N_dissolved(3)



  call poolTotals(soil%finewoodLitter,fastC=temp_fast_C,slowC=temp_slow_C,deadMicrobeC=temp_deadmic_C,liveMicrobeC=temp_livemic_C,&
                                        fastN=temp_fast_N,slowN=temp_slow_N,deadMicrobeN=temp_deadmic_N,liveMicrobeN=temp_livemic_N    )

    ! Are these even used for anything?
  total_fast_C=total_fast_C+temp_fast_C
  total_slow_C=total_slow_C+temp_slow_C
  total_deadmic_C=total_deadmic_C+temp_deadmic_C
  total_livemic_C=total_livemic_C+temp_livemic_C

  !Accumulate turnover rates for determining steady state pools
  if(temp_fast_C>0)soil%finewoodlitter_fast_C_turnover_accumulated=soil%finewoodlitter_fast_C_turnover_accumulated+finewoodlitter_fast_C_loss_rate/temp_fast_C
  if(temp_slow_C>0)soil%finewoodlitter_slow_C_turnover_accumulated=soil%finewoodlitter_slow_C_turnover_accumulated+finewoodlitter_slow_C_loss_rate/temp_slow_C
  if(temp_deadmic_C>0)soil%finewoodlitter_deadmic_C_turnover_accumulated=soil%finewoodlitter_deadmic_C_turnover_accumulated+finewoodlitter_deadmic_C_loss_rate/temp_deadmic_C
  soil%finewoodlitter_deadmic_C_in=soil%finewoodlitter_deadmic_C_in+finewoodlitter_deadmic_C_produced
  soil%finewoodlitter_fsc_in=soil%finewoodlitter_fsc_in+finewoodlitter_C_deposited(1)-finewoodlitter_C_dissolved(1)
  soil%finewoodlitter_ssc_in=soil%finewoodlitter_ssc_in+finewoodlitter_C_deposited(2)-finewoodlitter_C_dissolved(2)
  soil%finewoodlitter_deadmic_C_in=soil%finewoodlitter_deadmic_C_in+finewoodlitter_C_deposited(3)-finewoodlitter_C_dissolved(3)

  if(temp_fast_N>0)soil%finewoodlitter_fast_N_turnover_accumulated=soil%finewoodlitter_fast_N_turnover_accumulated+finewoodlitter_fast_N_loss_rate/temp_fast_N
  if(temp_slow_N>0)soil%finewoodlitter_slow_N_turnover_accumulated=soil%finewoodlitter_slow_N_turnover_accumulated+finewoodlitter_slow_N_loss_rate/temp_slow_N
  if(temp_deadmic_N>0)soil%finewoodlitter_deadmic_N_turnover_accumulated=soil%finewoodlitter_deadmic_N_turnover_accumulated+finewoodlitter_deadmic_N_loss_rate/temp_deadmic_N
  soil%finewoodlitter_deadmic_N_in=soil%finewoodlitter_deadmic_N_in+finewoodlitter_deadmic_N_produced
  soil%finewoodlitter_fsn_in=soil%finewoodlitter_fsn_in+finewoodlitter_N_deposited(1)-finewoodlitter_N_dissolved(1)
  soil%finewoodlitter_ssn_in=soil%finewoodlitter_ssn_in+finewoodlitter_N_deposited(2)-finewoodlitter_N_dissolved(2)
  soil%finewoodlitter_deadmic_N_in=soil%finewoodlitter_deadmic_N_in+finewoodlitter_N_deposited(3)-finewoodlitter_N_dissolved(3)

  call poolTotals(soil%coarsewoodLitter,fastC=temp_fast_C,slowC=temp_slow_C,deadMicrobeC=temp_deadmic_C,liveMicrobeC=temp_livemic_C,&
                                        fastN=temp_fast_N,slowN=temp_slow_N,deadMicrobeN=temp_deadmic_N,liveMicrobeN=temp_livemic_N    )

    ! Are these even used for anything?
  total_fast_C=total_fast_C+temp_fast_C
  total_slow_C=total_slow_C+temp_slow_C
  total_deadmic_C=total_deadmic_C+temp_deadmic_C
  total_livemic_C=total_livemic_C+temp_livemic_C

  !Accumulate turnover rates for determining steady state pools
  if(temp_fast_C>0)soil%coarsewoodlitter_fast_C_turnover_accumulated=soil%coarsewoodlitter_fast_C_turnover_accumulated+coarsewoodlitter_fast_C_loss_rate/temp_fast_C
  if(temp_slow_C>0)soil%coarsewoodlitter_slow_C_turnover_accumulated=soil%coarsewoodlitter_slow_C_turnover_accumulated+coarsewoodlitter_slow_C_loss_rate/temp_slow_C
  if(temp_deadmic_C>0)soil%coarsewoodlitter_deadmic_C_turnover_accumulated=soil%coarsewoodlitter_deadmic_C_turnover_accumulated+coarsewoodlitter_deadmic_C_loss_rate/temp_deadmic_C
  soil%coarsewoodlitter_deadmic_C_in=soil%coarsewoodlitter_deadmic_C_in+coarsewoodlitter_deadmic_C_produced
  soil%coarsewoodlitter_fsc_in=soil%coarsewoodlitter_fsc_in+coarsewoodlitter_C_deposited(1)-coarsewoodlitter_C_dissolved(1)
  soil%coarsewoodlitter_ssc_in=soil%coarsewoodlitter_ssc_in+coarsewoodlitter_C_deposited(2)-coarsewoodlitter_C_dissolved(2)
  soil%coarsewoodlitter_deadmic_C_in=soil%coarsewoodlitter_deadmic_C_in+coarsewoodlitter_C_deposited(3)-coarsewoodlitter_C_dissolved(3)

  if(temp_fast_N>0)soil%coarsewoodlitter_fast_N_turnover_accumulated=soil%coarsewoodlitter_fast_N_turnover_accumulated+coarsewoodlitter_fast_N_loss_rate/temp_fast_N
  if(temp_slow_N>0)soil%coarsewoodlitter_slow_N_turnover_accumulated=soil%coarsewoodlitter_slow_N_turnover_accumulated+coarsewoodlitter_slow_N_loss_rate/temp_slow_N
  if(temp_deadmic_N>0)soil%coarsewoodlitter_deadmic_N_turnover_accumulated=soil%coarsewoodlitter_deadmic_N_turnover_accumulated+coarsewoodlitter_deadmic_N_loss_rate/temp_deadmic_N
  soil%coarsewoodlitter_deadmic_N_in=soil%coarsewoodlitter_deadmic_N_in+coarsewoodlitter_deadmic_N_produced
  soil%coarsewoodlitter_fsn_in=soil%coarsewoodlitter_fsn_in+coarsewoodlitter_N_deposited(1)-coarsewoodlitter_N_dissolved(1)
  soil%coarsewoodlitter_ssn_in=soil%coarsewoodlitter_ssn_in+coarsewoodlitter_N_deposited(2)-coarsewoodlitter_N_dissolved(2)
  soil%coarsewoodlitter_deadmic_N_in=soil%coarsewoodlitter_deadmic_N_in+coarsewoodlitter_N_deposited(3)-coarsewoodlitter_N_dissolved(3)


  ! Next we have to go through layers and decompose the soil carbon pools
  do k=1,size(soil%soil_organic_matter)

      call update_pool(pool=soil%soil_organic_matter(k),T=decomp_T(k),theta=decomp_theta(k),air_filled_porosity=1.0-(decomp_theta(k)+ice_porosity(k)),&
                liquid_water=soil%wl(k),frozen_water=soil%ws(k),dt=dt_fast_yr,layerThickness=dz(k),&
                fast_C_loss_rate=fast_C_loss_rate(k),slow_C_loss_rate=slow_C_loss_rate(k), &
                deadmic_C_loss_rate=dead_microbe_C_loss_rate(k), CO2prod=CO2prod, &
                fast_N_loss_rate=fast_N_loss_rate(k),slow_N_loss_rate=slow_N_loss_rate(k), &
                deadmic_N_loss_rate=dead_microbe_N_loss_rate(k), &
                deadmic_C_produced=deadmic_C_produced(k), protected_C_produced=protected_C_produced(:,k), &
                protected_turnover_rate=protected_C_turnover_rate(:,k), &
                deadmic_N_produced=deadmic_N_produced(k), protected_N_produced=protected_N_produced(:,k), &
                protected_N_turnover_rate=protected_N_turnover_rate(:,k), &
                C_dissolved=C_dissolved(:,k), deposited_C=C_deposited(:,k), badCohort=badCohort,&
                N_dissolved=N_dissolved(:,k), deposited_N=N_deposited(:,k),&
                nitrification=soil_nitrif(k),denitrification=soil_denitrif(k),&
                N_mineralization=soil_N_mineralization(k),N_immobilization=soil_N_immobilization(k))


    IF (badCohort.ne.0) THEN
        call get_current_point(point_i,point_j,point_k,point_face)
        WRITE (*,*), 'Found bad cohort in layer',k,'Point i,j,k,face:',point_i,point_j,point_k,point_face
        WRITE (*,*), 'T=',decomp_T(k),'theta=',decomp_theta(k),'dt=',dt_fast_yr
        call error_mesg('Dsdt','Found bad cohort',FATAL)
    ENDIF

    vegn%rh=vegn%rh + CO2prod/dt_fast_yr
    call poolTotals(soil%soil_organic_matter(k),fastC=temp_fast_C,slowC=temp_slow_C,deadMicrobeC=temp_deadmic_C,&
            liveMicrobeC=temp_livemic_C,protectedC=temp_protected_C,&
            fast_protectedC=temp_protected_fast_C,slow_protectedC=temp_protected_slow_C,deadmic_protectedC=temp_protected_deadmic_C,&
            fastN=temp_fast_N,slowN=temp_slow_N,deadMicrobeN=temp_deadmic_N,&
                    liveMicrobeN=temp_livemic_N,protectedN=temp_protected_N,&
                    fast_protectedN=temp_protected_fast_N,slow_protectedN=temp_protected_slow_N,deadmic_protectedN=temp_protected_deadmic_N)
    total_fast_C=total_fast_C+temp_fast_C
    total_slow_C=total_slow_C+temp_slow_C
    total_deadmic_C=total_deadmic_C+temp_deadmic_C
    total_livemic_C=total_livemic_C+temp_livemic_C

    !Accumulate turnover rates for determining steady state pools
    if(temp_fast_C>0)soil%fast_C_turnover_accumulated(k)=soil%fast_C_turnover_accumulated(k)+fast_C_loss_rate(k)/temp_fast_C
    if(temp_slow_C>0)soil%slow_C_turnover_accumulated(k)=soil%slow_C_turnover_accumulated(k)+slow_C_loss_rate(k)/temp_slow_C
    if(temp_deadmic_C>0)soil%deadmic_C_turnover_accumulated(k)=soil%deadmic_C_turnover_accumulated(k)+dead_microbe_C_loss_rate(k)/temp_deadmic_C
    soil%fast_protected_C_in(k)=soil%fast_protected_C_in(k)+protected_C_produced(1,k)
    soil%slow_protected_C_in(k)=soil%slow_protected_C_in(k)+protected_C_produced(2,k)
    soil%deadmic_protected_C_in(k)=soil%deadmic_protected_C_in(k)+protected_C_produced(3,k)
    if(temp_protected_fast_C>0) soil%fast_protected_C_turnover_accumulated(k)=soil%fast_protected_C_turnover_accumulated(k)+protected_C_turnover_rate(1,k)/temp_protected_fast_C
    if(temp_protected_slow_C>0) soil%slow_protected_C_turnover_accumulated(k)=soil%slow_protected_C_turnover_accumulated(k)+protected_C_turnover_rate(2,k)/temp_protected_slow_C
    if(temp_protected_deadmic_C>0) soil%deadmic_protected_C_turnover_accumulated(k)=soil%deadmic_protected_C_turnover_accumulated(k)+protected_C_turnover_rate(3,k)/temp_protected_deadmic_C
    soil%deadmic_C_in(k)=soil%deadmic_C_in(k)+deadmic_C_produced(k)
    soil%fsc_in(k)=soil%fsc_in(k)+C_deposited(1,k)-C_dissolved(1,k)
    soil%ssc_in(k)=soil%ssc_in(k)+C_deposited(2,k)-C_dissolved(2,k)
    soil%deadmic_C_in(k)=soil%deadmic_C_in(k)+C_deposited(3,k)-C_dissolved(3,k)

    if(temp_fast_N>0)soil%fast_N_turnover_accumulated(k)=soil%fast_N_turnover_accumulated(k)+fast_N_loss_rate(k)/temp_fast_N
    if(temp_slow_N>0)soil%slow_N_turnover_accumulated(k)=soil%slow_N_turnover_accumulated(k)+slow_N_loss_rate(k)/temp_slow_N
    if(temp_deadmic_N>0)soil%deadmic_N_turnover_accumulated(k)=soil%deadmic_N_turnover_accumulated(k)+dead_microbe_N_loss_rate(k)/temp_deadmic_N
    soil%fast_protected_N_in(k)=soil%fast_protected_N_in(k)+protected_N_produced(1,k)
    soil%slow_protected_N_in(k)=soil%slow_protected_N_in(k)+protected_N_produced(2,k)
    soil%deadmic_protected_N_in(k)=soil%deadmic_protected_N_in(k)+protected_N_produced(3,k)
    if(temp_protected_fast_N>0) soil%fast_protected_N_turnover_accumulated(k)=soil%fast_protected_N_turnover_accumulated(k)+protected_N_turnover_rate(1,k)/temp_protected_fast_N
    if(temp_protected_slow_N>0) soil%slow_protected_N_turnover_accumulated(k)=soil%slow_protected_N_turnover_accumulated(k)+protected_N_turnover_rate(2,k)/temp_protected_slow_N
    if(temp_protected_deadmic_N>0) soil%deadmic_protected_N_turnover_accumulated(k)=soil%deadmic_protected_N_turnover_accumulated(k)+protected_N_turnover_rate(3,k)/temp_protected_deadmic_N
    soil%deadmic_N_in(k)=soil%deadmic_N_in(k)+deadmic_N_produced(k)
    soil%fsn_in(k)=soil%fsn_in(k)+N_deposited(1,k)-N_dissolved(1,k)
    soil%ssn_in(k)=soil%ssn_in(k)+N_deposited(2,k)-N_dissolved(2,k)
    soil%deadmic_N_in(k)=soil%deadmic_N_in(k)+N_deposited(3,k)-N_dissolved(3,k)
  enddo


  ! for budget check
  ! Still need to set these up for N -- BNS
  vegn%fsc_out = vegn%fsc_out + (sum(fast_C_loss_rate(:)) + leaflitter_fast_C_loss_rate + finewoodlitter_fast_C_loss_rate + coarsewoodlitter_fast_C_loss_rate)*dt_fast_yr
  vegn%ssc_out = vegn%ssc_out + (sum(slow_C_loss_rate(:)) + leaflitter_slow_C_loss_rate + finewoodlitter_slow_C_loss_rate + coarsewoodlitter_slow_C_loss_rate)*dt_fast_yr;


  ! accumulate decomposition rate reduction for the soil carbon restart output
  soil%asoil_in(:) = soil%asoil_in(:) + A(:)



  ! TODO: arithmetic averaging of A doesn't seem correct; we need to invent something better,
  !       e.g. weight it with the carbon loss, or something like that

  ! ---- diagnostic section
  if (id_rsoil_fast>0)  call send_tile_data(id_rsoil_fast, fast_C_loss_rate(:)/dz, diag)
  if (id_rsoil_slow>0)  call send_tile_data(id_rsoil_slow, slow_C_loss_rate(:)/dz, diag)
  if (id_rsoil_deadmic>0) call send_tile_data(id_rsoil_deadmic, dead_microbe_C_loss_rate(:)/dz, diag)
  if (id_rsoil_leaflitter_fast>0) call send_tile_data(id_rsoil_leaflitter_fast, leaflitter_fast_C_loss_rate, diag)
  if (id_rsoil_leaflitter_slow>0) call send_tile_data(id_rsoil_leaflitter_slow, leaflitter_slow_C_loss_rate, diag)
  if (id_rsoil_leaflitter_deadmic>0) call send_tile_data(id_rsoil_leaflitter_deadmic, leaflitter_deadmic_C_loss_rate, diag)
  if (id_rsoil_finewoodlitter_fast>0) call send_tile_data(id_rsoil_finewoodlitter_fast, finewoodlitter_fast_C_loss_rate, diag)
  if (id_rsoil_finewoodlitter_slow>0) call send_tile_data(id_rsoil_finewoodlitter_slow, finewoodlitter_slow_C_loss_rate, diag)
  if (id_rsoil_finewoodlitter_deadmic>0) call send_tile_data(id_rsoil_finewoodlitter_deadmic, finewoodlitter_deadmic_C_loss_rate, diag)
  if (id_rsoil_coarsewoodlitter_fast>0) call send_tile_data(id_rsoil_coarsewoodlitter_fast, coarsewoodlitter_fast_C_loss_rate, diag)
  if (id_rsoil_coarsewoodlitter_slow>0) call send_tile_data(id_rsoil_coarsewoodlitter_slow, coarsewoodlitter_slow_C_loss_rate, diag)
  if (id_rsoil_coarsewoodlitter_deadmic>0) call send_tile_data(id_rsoil_coarsewoodlitter_deadmic, coarsewoodlitter_deadmic_C_loss_rate, diag)
  call send_tile_data(id_rsoil, vegn%rh, diag)
  ! TODO: arithmetic averaging of A doesn't seem correct; we need to invent something better,
  !       e.g. weight it with the carbon loss, or something like that
  if (id_asoil>0) call send_tile_data(id_asoil, sum(A(:))/size(A(:)), diag)


! Should we not send these if N turned off? Currently they should all be zero
  if (id_rsoil_fast_N>0)  call send_tile_data(id_rsoil_fast_N, fast_N_loss_rate(:)/dz, diag)
  if (id_rsoil_slow_N>0)  call send_tile_data(id_rsoil_slow_N, slow_N_loss_rate(:)/dz, diag)
  if (id_rsoil_deadmic_N>0) call send_tile_data(id_rsoil_deadmic_N, dead_microbe_N_loss_rate(:)/dz, diag)
  if (id_rsoil_leaflitter_fast_N>0) call send_tile_data(id_rsoil_leaflitter_fast_N, leaflitter_fast_N_loss_rate, diag)
  if (id_rsoil_leaflitter_slow_N>0) call send_tile_data(id_rsoil_leaflitter_slow_N, leaflitter_slow_N_loss_rate, diag)
  if (id_rsoil_leaflitter_deadmic_N>0) call send_tile_data(id_rsoil_leaflitter_deadmic_N, leaflitter_deadmic_N_loss_rate, diag)
  if (id_rsoil_finewoodlitter_fast_N>0) call send_tile_data(id_rsoil_finewoodlitter_fast_N, finewoodlitter_fast_N_loss_rate, diag)
  if (id_rsoil_finewoodlitter_slow_N>0) call send_tile_data(id_rsoil_finewoodlitter_slow_N, finewoodlitter_slow_N_loss_rate, diag)
  if (id_rsoil_finewoodlitter_deadmic_N>0) call send_tile_data(id_rsoil_finewoodlitter_deadmic_N, finewoodlitter_deadmic_N_loss_rate, diag)
  if (id_rsoil_coarsewoodlitter_fast_N>0) call send_tile_data(id_rsoil_coarsewoodlitter_fast_N, coarsewoodlitter_fast_N_loss_rate, diag)
  if (id_rsoil_coarsewoodlitter_slow_N>0) call send_tile_data(id_rsoil_coarsewoodlitter_slow_N, coarsewoodlitter_slow_N_loss_rate, diag)
  if (id_rsoil_coarsewoodlitter_deadmic_N>0) call send_tile_data(id_rsoil_coarsewoodlitter_deadmic_N, coarsewoodlitter_deadmic_N_loss_rate, diag)

  if (id_nitrification_rate>0) call send_tile_data(id_nitrification_rate, soil_nitrif/dt_fast_yr/dz, diag)
  if (id_denitrification_rate>0) call send_tile_data(id_denitrification_rate, soil_denitrif/dt_fast_yr/dz, diag)
  if (id_N_mineralization_rate>0) call send_tile_data(id_N_mineralization_rate, soil_N_mineralization/dt_fast_yr/dz,diag)
  if (id_N_immobilization_rate>0) call send_tile_data(id_N_immobilization_rate, soil_N_immobilization/dt_fast_yr/dz,diag)

  if (id_leaflitter_nitrification_rate>0) call send_tile_data(id_leaflitter_nitrification_rate, leaflitter_nitrif/dt_fast_yr, diag)
  if (id_leaflitter_denitrification_rate>0) call send_tile_data(id_leaflitter_denitrification_rate, leaflitter_denitrif/dt_fast_yr, diag)
  if (id_leaflitter_N_mineralization_rate>0) call send_tile_data(id_leaflitter_N_mineralization_rate, leaflitter_N_mineralization/dt_fast_yr,diag)
  if (id_leaflitter_N_immobilization_rate>0) call send_tile_data(id_leaflitter_N_immobilization_rate, leaflitter_N_immobilization/dt_fast_yr,diag)
  if (id_finewoodlitter_nitrification_rate>0) call send_tile_data(id_finewoodlitter_nitrification_rate, finewoodlitter_nitrif/dt_fast_yr, diag)
  if (id_finewoodlitter_denitrification_rate>0) call send_tile_data(id_finewoodlitter_denitrification_rate, finewoodlitter_denitrif/dt_fast_yr, diag)
  if (id_finewoodlitter_N_mineralization_rate>0) call send_tile_data(id_finewoodlitter_N_mineralization_rate, finewoodlitter_N_mineralization/dt_fast_yr,diag)
  if (id_finewoodlitter_N_immobilization_rate>0) call send_tile_data(id_finewoodlitter_N_immobilization_rate, finewoodlitter_N_immobilization/dt_fast_yr,diag)
  if (id_coarsewoodlitter_nitrification_rate>0) call send_tile_data(id_coarsewoodlitter_nitrification_rate, coarsewoodlitter_nitrif/dt_fast_yr, diag)
  if (id_coarsewoodlitter_denitrification_rate>0) call send_tile_data(id_coarsewoodlitter_denitrification_rate, coarsewoodlitter_denitrif/dt_fast_yr, diag)
  if (id_coarsewoodlitter_N_mineralization_rate>0) call send_tile_data(id_coarsewoodlitter_N_mineralization_rate, coarsewoodlitter_N_mineralization/dt_fast_yr,diag)
  if (id_coarsewoodlitter_N_immobilization_rate>0) call send_tile_data(id_coarsewoodlitter_N_immobilization_rate, coarsewoodlitter_N_immobilization/dt_fast_yr,diag)

  if (id_total_nitrification_rate>0) call send_tile_data(id_total_nitrification_rate, &
                (sum(soil_nitrif)+leaflitter_nitrif+finewoodlitter_nitrif+coarsewoodlitter_nitrif)/dt_fast_yr,diag)
  if (id_total_denitrification_rate>0) call send_tile_data(id_total_denitrification_rate, &
               (sum(soil_denitrif)+leaflitter_denitrif+finewoodlitter_denitrif+coarsewoodlitter_denitrif)/dt_fast_yr,diag)
  if (id_total_N_immobilization_rate>0) call send_tile_data(id_total_N_immobilization_rate, &
               (sum(soil_N_immobilization)+leaflitter_N_immobilization+finewoodlitter_N_immobilization+coarsewoodlitter_N_immobilization)/dt_fast_yr,diag)
  if (id_total_N_mineralization_rate>0) call send_tile_data(id_total_N_mineralization_rate, &
               (sum(soil_N_mineralization)+leaflitter_N_mineralization+finewoodlitter_N_mineralization+coarsewoodlitter_N_mineralization)/dt_fast_yr,diag)


  if (id_leaflitter_C_dissolve_rate_fast>0) call send_tile_data(id_leaflitter_C_dissolve_rate_fast,leaflitter_C_dissolved(1)/dt_fast_yr,diag)
  if (id_leaflitter_C_dissolve_rate_slow>0) call send_tile_data(id_leaflitter_C_dissolve_rate_slow,leaflitter_C_dissolved(2)/dt_fast_yr,diag)
  if (id_leaflitter_C_dissolve_rate_deadmic>0) call send_tile_data(id_leaflitter_C_dissolve_rate_deadmic,leaflitter_C_dissolved(3)/dt_fast_yr,diag)
  if (id_finewoodlitter_C_dissolve_rate_fast>0) call send_tile_data(id_finewoodlitter_C_dissolve_rate_fast,finewoodlitter_C_dissolved(1)/dt_fast_yr,diag)
  if (id_finewoodlitter_C_dissolve_rate_slow>0) call send_tile_data(id_finewoodlitter_C_dissolve_rate_slow,finewoodlitter_C_dissolved(2)/dt_fast_yr,diag)
  if (id_finewoodlitter_C_dissolve_rate_deadmic>0) call send_tile_data(id_finewoodlitter_C_dissolve_rate_deadmic,finewoodlitter_C_dissolved(3)/dt_fast_yr,diag)
  if (id_coarsewoodlitter_C_dissolve_rate_fast>0) call send_tile_data(id_coarsewoodlitter_C_dissolve_rate_fast,coarsewoodlitter_C_dissolved(1)/dt_fast_yr,diag)
  if (id_coarsewoodlitter_C_dissolve_rate_slow>0) call send_tile_data(id_coarsewoodlitter_C_dissolve_rate_slow,coarsewoodlitter_C_dissolved(2)/dt_fast_yr,diag)
  if (id_coarsewoodlitter_C_dissolve_rate_deadmic>0) call send_tile_data(id_coarsewoodlitter_C_dissolve_rate_deadmic,coarsewoodlitter_C_dissolved(3)/dt_fast_yr,diag)
  if (id_C_dissolve_rate_fast>0) call send_tile_data(id_C_dissolve_rate_fast,C_dissolved(1,:)/dt_fast_yr/dz,diag)
  if (id_C_dissolve_rate_slow>0) call send_tile_data(id_C_dissolve_rate_slow,C_dissolved(2,:)/dt_fast_yr/dz,diag)
  if (id_C_dissolve_rate_deadmic>0) call send_tile_data(id_C_dissolve_rate_deadmic,C_dissolved(3,:)/dt_fast_yr/dz,diag)

  if (id_leaflitter_C_deposition_fast>0) call send_tile_data(id_leaflitter_C_deposition_fast,leaflitter_C_deposited(1)/dt_fast_yr,diag)
  if (id_leaflitter_C_deposition_slow>0) call send_tile_data(id_leaflitter_C_deposition_slow,leaflitter_C_deposited(2)/dt_fast_yr,diag)
  if (id_leaflitter_C_deposition_deadmic>0) call send_tile_data(id_leaflitter_C_deposition_deadmic,leaflitter_C_deposited(3)/dt_fast_yr,diag)
  if (id_finewoodlitter_C_deposition_fast>0) call send_tile_data(id_finewoodlitter_C_deposition_fast,finewoodlitter_C_deposited(1)/dt_fast_yr,diag)
  if (id_finewoodlitter_C_deposition_slow>0) call send_tile_data(id_finewoodlitter_C_deposition_slow,finewoodlitter_C_deposited(2)/dt_fast_yr,diag)
  if (id_finewoodlitter_C_deposition_deadmic>0) call send_tile_data(id_finewoodlitter_C_deposition_deadmic,finewoodlitter_C_deposited(3)/dt_fast_yr,diag)
  if (id_coarsewoodlitter_C_deposition_fast>0) call send_tile_data(id_coarsewoodlitter_C_deposition_fast,coarsewoodlitter_C_deposited(1)/dt_fast_yr,diag)
  if (id_coarsewoodlitter_C_deposition_slow>0) call send_tile_data(id_coarsewoodlitter_C_deposition_slow,coarsewoodlitter_C_deposited(2)/dt_fast_yr,diag)
  if (id_coarsewoodlitter_C_deposition_deadmic>0) call send_tile_data(id_coarsewoodlitter_C_deposition_deadmic,coarsewoodlitter_C_deposited(3)/dt_fast_yr,diag)

  if (id_C_deposition_fast>0) call send_tile_data(id_C_deposition_fast,C_deposited(1,:)/dt_fast_yr/dz,diag)
  if (id_C_deposition_slow>0) call send_tile_data(id_C_deposition_slow,C_deposited(2,:)/dt_fast_yr/dz,diag)
  if (id_C_deposition_deadmic>0) call send_tile_data(id_C_deposition_deadmic,C_deposited(3,:)/dt_fast_yr/dz,diag)


  if (id_leaflitter_N_dissolve_rate_fast>0) call send_tile_data(id_leaflitter_N_dissolve_rate_fast,leaflitter_N_dissolved(1)/dt_fast_yr,diag)
  if (id_leaflitter_N_dissolve_rate_slow>0) call send_tile_data(id_leaflitter_N_dissolve_rate_slow,leaflitter_N_dissolved(2)/dt_fast_yr,diag)
  if (id_leaflitter_N_dissolve_rate_deadmic>0) call send_tile_data(id_leaflitter_N_dissolve_rate_deadmic,leaflitter_N_dissolved(3)/dt_fast_yr,diag)
  if (id_finewoodlitter_N_dissolve_rate_fast>0) call send_tile_data(id_finewoodlitter_N_dissolve_rate_fast,finewoodlitter_N_dissolved(1)/dt_fast_yr,diag)
  if (id_finewoodlitter_N_dissolve_rate_slow>0) call send_tile_data(id_finewoodlitter_N_dissolve_rate_slow,finewoodlitter_N_dissolved(2)/dt_fast_yr,diag)
  if (id_finewoodlitter_N_dissolve_rate_deadmic>0) call send_tile_data(id_finewoodlitter_N_dissolve_rate_deadmic,finewoodlitter_N_dissolved(3)/dt_fast_yr,diag)
  if (id_coarsewoodlitter_N_dissolve_rate_fast>0) call send_tile_data(id_coarsewoodlitter_N_dissolve_rate_fast,coarsewoodlitter_N_dissolved(1)/dt_fast_yr,diag)
  if (id_coarsewoodlitter_N_dissolve_rate_slow>0) call send_tile_data(id_coarsewoodlitter_N_dissolve_rate_slow,coarsewoodlitter_N_dissolved(2)/dt_fast_yr,diag)
  if (id_coarsewoodlitter_N_dissolve_rate_deadmic>0) call send_tile_data(id_coarsewoodlitter_N_dissolve_rate_deadmic,coarsewoodlitter_N_dissolved(3)/dt_fast_yr,diag)
  if (id_N_dissolve_rate_fast>0) call send_tile_data(id_N_dissolve_rate_fast,C_dissolved(1,:)/dt_fast_yr/dz,diag)
  if (id_N_dissolve_rate_slow>0) call send_tile_data(id_C_dissolve_rate_slow,C_dissolved(2,:)/dt_fast_yr/dz,diag)
  if (id_N_dissolve_rate_deadmic>0) call send_tile_data(id_N_dissolve_rate_deadmic,C_dissolved(3,:)/dt_fast_yr/dz,diag)

  if (id_leaflitter_N_deposition_fast>0) call send_tile_data(id_leaflitter_N_deposition_fast,leaflitter_N_deposited(1)/dt_fast_yr,diag)
  if (id_leaflitter_N_deposition_slow>0) call send_tile_data(id_leaflitter_N_deposition_slow,leaflitter_N_deposited(2)/dt_fast_yr,diag)
  if (id_leaflitter_N_deposition_deadmic>0) call send_tile_data(id_leaflitter_N_deposition_deadmic,leaflitter_N_deposited(3)/dt_fast_yr,diag)
  if (id_finewoodlitter_N_deposition_fast>0) call send_tile_data(id_finewoodlitter_N_deposition_fast,finewoodlitter_N_deposited(1)/dt_fast_yr,diag)
  if (id_finewoodlitter_N_deposition_slow>0) call send_tile_data(id_finewoodlitter_N_deposition_slow,finewoodlitter_N_deposited(2)/dt_fast_yr,diag)
  if (id_finewoodlitter_N_deposition_deadmic>0) call send_tile_data(id_finewoodlitter_N_deposition_deadmic,finewoodlitter_N_deposited(3)/dt_fast_yr,diag)
  if (id_coarsewoodlitter_N_deposition_fast>0) call send_tile_data(id_coarsewoodlitter_N_deposition_fast,coarsewoodlitter_N_deposited(1)/dt_fast_yr,diag)
  if (id_coarsewoodlitter_N_deposition_slow>0) call send_tile_data(id_coarsewoodlitter_N_deposition_slow,coarsewoodlitter_N_deposited(2)/dt_fast_yr,diag)
  if (id_coarsewoodlitter_N_deposition_deadmic>0) call send_tile_data(id_coarsewoodlitter_N_deposition_deadmic,coarsewoodlitter_N_deposited(3)/dt_fast_yr,diag)

  if (id_N_deposition_fast>0) call send_tile_data(id_N_deposition_fast,N_deposited(1,:)/dt_fast_yr/dz,diag)
  if (id_N_deposition_slow>0) call send_tile_data(id_N_deposition_slow,N_deposited(2,:)/dt_fast_yr/dz,diag)
  if (id_N_deposition_deadmic>0) call send_tile_data(id_N_deposition_deadmic,N_deposited(3,:)/dt_fast_yr/dz,diag)

  if (id_root_profile > 0) then
    call vegn_uptake_profile (vegn, dz(1:num_l), uptake_frac_max, vegn_uptake_term )
    call send_tile_data(id_root_profile, uptake_frac_max, diag)
  endif
end subroutine Dsdt_CORPSE


! ============================================================================
subroutine Dsdt_CENTURY(vegn, soil, diag, soilt, theta)
  type(vegn_tile_type), intent(inout) :: vegn
  type(soil_tile_type), intent(inout) :: soil
  type(diag_buff_type), intent(inout) :: diag
  real                , intent(in)    :: soilt ! average soil temperature, deg K
  real                , intent(in)    :: theta ! average soil moisture

  real :: fast_C_loss(size(soil%fast_soil_C))
  real :: slow_C_loss(size(soil%slow_soil_C))
  real :: A          (size(soil%slow_soil_C)) ! decomp rate reduction due to moisture and temperature

  select case (soil_carbon_option)
  case(SOILC_CENTURY)
      A(:) = A_function(soilt, theta)
  case(SOILC_CENTURY_BY_LAYER)
      A(:) = A_function(soil%T, soil_theta(soil))
  case default
    call error_mesg('Dsdt_CENTURY','The value of soil_carbon_option is invalid. This should never happen. See developer.',FATAL)
  end select

  fast_C_loss = soil%fast_soil_C(:)*A*K1*dt_fast_yr;
  slow_C_loss = soil%slow_soil_C(:)*A*K2*dt_fast_yr;

  soil%fast_soil_C = soil%fast_soil_C - fast_C_loss;
  soil%slow_soil_C = soil%slow_soil_C - slow_C_loss;

  ! for budget check
  vegn%fsc_out = vegn%fsc_out + sum(fast_C_loss(:));
  vegn%ssc_out = vegn%ssc_out + sum(slow_C_loss(:));

  ! loss of C to atmosphere and leaching
  vegn%rh = sum(fast_C_loss(:)+slow_C_loss(:))/dt_fast_yr;

  ! accumulate decomposition rate reduction for the soil carbon restart output
  soil%asoil_in(:) = soil%asoil_in(:) + A(:)
  ! TODO: arithmetic averaging of A doesn't seem correct; we need to invent something better,
  !       e.g. weight it with the carbon loss, or something like that

  ! ---- diagnostic section
  if (id_fsc>0)         call send_tile_data(id_fsc, sum(soil%fast_soil_C(:)), diag)
  if (id_ssc>0)         call send_tile_data(id_ssc, sum(soil%slow_soil_C(:)), diag)
  if (id_rsoil_fast>0)  call send_tile_data(id_rsoil_fast, sum(fast_C_loss(:))/(dz(:)*dt_fast_yr), diag)
  if (id_rsoil_slow>0)  call send_tile_data(id_rsoil_slow, sum(slow_C_loss(:))/(dz(:)*dt_fast_yr), diag)

  call send_tile_data(id_fast_soil_C, soil%fast_soil_C(:)/dz(1:num_l), diag)
  call send_tile_data(id_slow_soil_C, soil%slow_soil_C(:)/dz(1:num_l), diag)
  call send_tile_data(id_rsoil, vegn%rh, diag)
  ! TODO: arithmetic averaging of A doesn't seem correct; we need to invent something better,
  !       e.g. weight it with the carbon loss, or something like that
  if (id_asoil>0) call send_tile_data(id_asoil, sum(A(:))/size(A(:)), diag)

end subroutine Dsdt_CENTURY


! ============================================================================
subroutine soil_push_down_excess ( soil, diag, lrunf_nu, hlrunf_nu, frunf, hfrunf)
  type(soil_tile_type), intent(inout) :: soil
  type(diag_buff_type), intent(inout) :: diag
  real, intent(out) :: lrunf_nu, hlrunf_nu
  real, intent(out) :: frunf, hfrunf ! frozen runoff (mm/s); frozen runoff heat (W/m^2)

  ! ---- local vars ----------------------------------------------------------
  real      :: &
     liq_frac, excess_wat, excess_liq, excess_ice, excess_t, &
     h1, h2, summax, space_avail, liq_placed, ice_placed
  integer :: k,l

  liq_frac=0;excess_wat=0;excess_liq=0;excess_ice=0;h1=0;h2=0
  l = 1
  summax = max(0.,soil%wl(l))+max(0.,soil%ws(l))
  if (summax > 0) then
     liq_frac = max(0.,soil%wl(l)) / summax
  else
     liq_frac = 1
  endif
  excess_wat = max(0., soil%wl(l) + soil%ws(l) &
       - dens_h2o*dz(l)*soil%vwc_max(l) )
  excess_liq = excess_wat*liq_frac
  excess_ice = excess_wat-excess_liq
  excess_t   = soil%T(l)
  soil%wl(l) = soil%wl(l) - excess_liq
  soil%ws(l) = soil%ws(l) - excess_ice
  call send_tile_data(id_excess, excess_wat/delta_time, diag)

!  if(is_watch_cell()) then
  if(is_watch_point()) then
     write(*,*) ' ##### push_down_excess input #####'
      call get_current_point(k=k)
      write(*,*) 'For watch_cell, tile=',k
     __DEBUG3__(l,summax,liq_frac)
     __DEBUG3__(soil%vwc_max(l),excess_liq,excess_ice)
     __DEBUG2__(dens_h2o,dz(l))
  endif

  do l = 2, num_l
     if (excess_liq+excess_ice>0) then
        space_avail = dens_h2o*dz(l)*soil%vwc_max(l) &
             - (soil%wl(l) + soil%ws(l))
        liq_placed = max(min(space_avail, excess_liq), 0.)
        ice_placed = max(min(space_avail-liq_placed, excess_ice), 0.)
        h1 = (soil%heat_capacity_dry(l)*dz(l) &
             + csw*soil%ws(l) + clw*soil%wl(l))
        h2 = liq_placed*clw+ice_placed*csw
        soil%T(l) = (h1 * soil%T(l) &
             + h2 * excess_T )  / (h1+h2)
        soil%wl(l) = soil%wl(l) + liq_placed
        soil%ws(l) = soil%ws(l) + ice_placed
        excess_liq = excess_liq - liq_placed
        excess_ice = excess_ice - ice_placed
     endif
  enddo

! to avoid adding frozen runoff to soil interface, melt all remaining
! excess ice, even if it results in supercooled liquid runoff
  if (supercooled_rnu) then
     lrunf_nu  = (excess_liq+excess_ice) / delta_time
     hlrunf_nu = (  excess_liq*clw*(excess_T-tfreeze)  &
                  + excess_ice*csw*(excess_T-tfreeze)  &
                  - hlf_factor*hlf*excess_ice                   ) / delta_time
     frunf = 0.
     hfrunf = 0.
  else
!! ZMS: Change this to propagate frozen runoff to avoid lake crashes from supercooled water.
     lrunf_nu = excess_liq / delta_time
     hlrunf_nu = excess_liq*clw*(excess_T-tfreeze) / delta_time
     frunf = excess_ice / delta_time
     hfrunf = excess_ice*csw*(excess_T-tfreeze) / delta_time
  end if

!  if(is_watch_cell()) then
  if(is_watch_point()) then
     write(*,*) ' ##### push_down_excess output #####'
     __DEBUG2__(lrunf_nu,hlrunf_nu)
     write(*,*) 'For watch_cell'
     __DEBUG3__(soil%hidx_k, frunf, hfrunf)
     do l = 1, num_l
        write(*,'(x,a,x,i2.2)',advance='NO')' level=', l
        call dpri(' T =',soil%T(l))
        call dpri(' Th=', (soil%ws(l)+soil%wl(l))/(dens_h2o*dz(l)))
        call dpri(' wl=', soil%wl(l))
        call dpri(' ws=', soil%ws(l))
        write(*,*)
     enddo
  endif
end subroutine soil_push_down_excess

! ============================================================================
  subroutine RICHARDS_clean(soil, psi, DThDP, hyd_cond, DKDP, div, &
                lprec_eff, Dpsi_min, Dpsi_max, dt_richards, &
                 dPsi, dW_l, flow, lrunf_ie)
! ZMS How to set prevent lrunf_ie when surface soil water < pond?

  type(soil_tile_type), intent(inout)   :: soil
  real, intent(in),  dimension(num_l)   :: psi, DThDP, hyd_cond, DKDP, div
  real, intent(in)                      :: lprec_eff, Dpsi_min, Dpsi_max
  real, intent(in)                      :: dt_richards
  real, intent(out), dimension(num_l)   :: dPsi, dW_l
  real, intent(out), dimension(num_l+1) :: flow
  real, intent(out)                     :: lrunf_ie
  ! ---- local vars ----------------------------------------------------------
  integer l, ipt, jpt, kpt, fpt, l_internal
  real, dimension(num_l-1) :: del_z, K, DKDPm, DKDPp, grad, eee, fff
  real aaa, bbb, ccc, ddd, xxx, dpsi_alt, dW_l_internal, w_to_move_up, adj
  logical flag

  flag = .false.
  flow(1) = dt_richards*lprec_eff
  do l = 1, num_l-1
     del_z(l) = zfull(l+1)-zfull(l)
     if (harmonic_mean_K) then
         K(l) = (dz(l)+dz(l+1))/(dz(l)/hyd_cond(l)+dz(l+1)/hyd_cond(l+1))
     else
         K(l) = 0.5*(hyd_cond(l)+hyd_cond(l+1))
     endif
     DKDPm(l) = 0. ! 0.5*DKDP(l)
     DKDPp(l) = 0. ! 0.5*DKDP(l+1)
!        K(l) = hyd_cond(l)
!        DKDPm(l) = DKDP(l)
!        DKDPp(l) = 0
      grad(l)  = (psi(l+1)-psi(l))/del_z(l) - 1
  enddo

  if (is_watch_point()) then
     write(*,*) '##### soil_step_2 checkpoint 3.1 #####'
     do l = 1, num_l
        write(*,'(x,a,x,i2.2,x,a,100(x,g23.16))') 'level=', l, 'DThDP,hyd_cond,psi,DKDP', &
             DThDP(l),&
             hyd_cond(l),&
             psi(l),&
             DKDP(l)
     enddo
     do l = 1, num_l-1
        write(*,'(a,i2.2,1x,a,100(2x,g23.16))') 'interface=', l, 'K,DKDPm,DKDPp,grad,del_z', &
             K(l),&
             DKDPm(l),&
             DKDPp(l),&
             grad(l)
     enddo
  endif


  l = num_l
  xxx = dens_h2o*dz(l)*DThDP(l)/dt_richards
  aaa =     - ( K(l-1)/del_z(l-1) - DKDPm(l-1)*grad(l-1))
  bbb = xxx - (-K(l-1)/del_z(l-1) - DKDPp(l-1)*grad(l-1) )
  ddd = - K(l-1) *grad(l-1) - div(l)
  eee(l-1) = -aaa/bbb
  fff(l-1) =  ddd/bbb

  if(is_watch_point()) then
       write(*,'(a,i2.2,100(2x,g23.16))') 'l,a,b, ,d', l,aaa, bbb,ddd
  endif

  do l = num_l-1, 2, -1
    xxx = dens_h2o*dz(l)*DThDP(l)/dt_richards
    aaa = - ( K(l-1)/del_z(l-1) - DKDPm(l-1)*grad(l-1))
    bbb = xxx-( -K(l-1)/del_z(l-1) - DKDPp(l-1)*grad(l-1)&
                -K(l  )/del_z(l  ) + DKDPm(l  )*grad(l  ))
    ccc =   - (  K(l  )/del_z(l  ) + DKDPp(l  )*grad(l  ))
    ddd =       K(l)*grad(l) - K(l-1)*grad(l-1) &
                          - div(l)
    eee(l-1) =                    -aaa/(bbb+ccc*eee(l))
    fff(l-1) =  (ddd-ccc*fff(l))/(bbb+ccc*eee(l))
    if(is_watch_point()) then
       write(*,'(a,i2.2,100(2x,g23.16))') 'l,a,b,c,d', l,aaa, bbb,ccc,ddd
    endif
  enddo

  l = 1
  xxx = dens_h2o*dz(l)*DThDP(l)/dt_richards
  bbb = xxx - ( -K(l  )/del_z(l  ) + DKDPm(l  )*grad(l  ))
  ccc =     - (  K(l  )/del_z(l  ) + DKDPp(l  )*grad(l  ))
  ddd =          flow(1)/dt_richards +    K(l)     *grad(l) &
                          - div(l)

  if (Dpsi_min.ge.Dpsi_max) call error_mesg(module_name, '=== Dpsi_min.ge.Dpsi_max', FATAL)

  IF (bbb+ccc*eee(l) .EQ. 0.) THEN
      call get_current_point(ipt,jpt,kpt,fpt)
      write(*,*) '===richards b+ce=0 ===','at point ',ipt,jpt,kpt,fpt
        write(*,*) 'bbb+ccc*eee(l) .EQ. 0.'
        write(*,*) 'bbb', bbb
        write(*,*) 'ccc', ccc
        write(*,*) 'ddd', ddd
        write(*,*) 'eee(l)', eee(l)
        write(*,*) 'fff(l)', fff(l)
        write(*,*) 'dPsi(l)', dPsi(l)
        write(*,*) 'dPsi(l)', dPsi(l)
        write(*,*) 'Dpsi_min', Dpsi_min
        write(*,*) 'Dpsi_max', Dpsi_max
    call error_mesg(module_name, 'b+ce=0 in soil-water equations', FATAL)
  ENDIF

     dPsi(l) = (ddd-ccc*fff(l))/(bbb+ccc*eee(l))
     if(is_watch_point()) then
        write(*,*) 'bbb+ccc*eee(l) .NE. 0.'
        write(*,*) 'bbb', bbb
        write(*,*) 'ccc', ccc
        write(*,*) 'ddd', ddd
        write(*,*) 'eee(l)', eee(l)
        write(*,*) 'fff(l)', fff(l)
        write(*,*) 'dPsi(l)', dPsi(l)
        write(*,*) 'dPsi(l)', dPsi(l)
        write(*,*) 'Dpsi_min', Dpsi_min
        write(*,*) 'Dpsi_max', Dpsi_max
     endif
     if (dPsi(l).lt.Dpsi_max) then
        lrunf_ie = 0.
     else
        dPsi(l) = Dpsi_max
        if (div_bug_fix) then
            flow(l) = (dPsi(l)*(bbb+ccc*eee(l))+ccc*fff(l)+div(l) &
                     - K(l)*grad(l))*dt_richards
        else
            flow(l) = (dPsi(l)*(bbb+ccc*eee(l))+ccc*fff(l) &
                     - K(l)*grad(l))*dt_richards
        endif
        lrunf_ie = lprec_eff - flow(l)/dt_richards
     endif

  if(is_watch_point().or.(flag.and.write_when_flagged)) then
     write(*,'(a,i2.2,100(2x,g23.16))') 'l,  b,c,d', l, bbb,ccc,ddd
     write(*,*) ' ##### soil_step_2 checkpoint 3.2 #####'
     write(*,*) 'ie:', lrunf_ie
     do l = 1, num_l-1
        write(*,'(a,i2.2,100(2x,g23.16))') 'l,eee(l),fff(l)',l,eee(l),fff(l)
     enddo
     write(*,*) 'DThDP(1)', DThDP(1)
     write(*,*) 'K(1)', K(1)
     write(*,*) 'grad(1)', grad(1)
     write(*,*) 'ddd(1)', ddd
     write(*,*) 'ccc(1)', ccc
     write(*,*) 'bbb(1)', bbb
     write(*,*) 'dPsi(1)', dPsi(1)
     write(*,*) 'Psi(1)', Psi(1)
     write(*,*) 'div(1)', div(1)
  endif

  do l = 1, num_l-1
     dPsi(l+1) = eee(l)*dPsi(l) + fff(l)
     flow(l+1) = dt_richards*( &
         -K(l)*(grad(l)&
         +(DPsi(l+1)-DPsi(l))/ del_z(l)) &
         -grad(l)*(DKDPp(l)*Dpsi(l+1)+ &
                         DKDPm(l)*Dpsi(l) )  )
     dW_l(l) = flow(l) - flow(l+1) - div(l)*dt_richards
  enddo
  flow(num_l+1) = 0.
  dW_l(num_l) = flow(num_l) - flow(num_l+1) &
                          - div(num_l)*dt_richards

  if(is_watch_point().or.(flag.and.write_when_flagged)) then
     write(*,*) ' ##### soil_step_2 checkpoint 3.21 #####'
     do l = 1, num_l
        write(*,'(i2.2,100(2x,a,g23.16))') l,&
             ' dW_l=', dW_l(l),&
             ' flow=', flow(l),&
             ' div=', div(l)
     enddo
  endif

  if (dPsi(1).lt.Dpsi_min) then
     if (verbose) then
         call get_current_point(ipt,jpt,kpt,fpt)
         write(*,*) '=== warning: dPsi=',dPsi(1),'<min=',dPsi_min,'at',ipt,jpt,kpt,fpt
       endif
     l_internal = 1
     dW_l_internal = -1.e20
     do l = 2, num_l
        if (dW_l(l).gt.dW_l_internal) then
           l_internal = l
           dW_l_internal = dW_l(l)
        endif
     enddo
     w_to_move_up = min(dW_l_internal, -(soil%wl(1)+dW_l(1)))
     w_to_move_up = max(w_to_move_up, 0.)
     write(*,*) 'l_internal=',l_internal
     write(*,*) 'dW_l(l_internal)=',dW_l(l_internal)
     write(*,*) 'soil%wl(1)+dW_l(1)',soil%wl(1)+dW_l(1)
     write(*,*) 'w_to_move_up=',w_to_move_up
     if (l_internal.gt.1) then
        dW_l(1) = dW_l(1) + w_to_move_up
        dW_l(l_internal) = dW_l(l_internal) - w_to_move_up
        do l = 2, l_internal
           flow(l) = flow(l) - w_to_move_up
        enddo
     endif
  endif

  if(is_watch_point().or.(flag.and.write_when_flagged)) then
     write(*,*) ' ##### soil_step_2 checkpoint 3.22 #####'
     do l = 1, num_l
        write(*,'(i2.2,100(2x,a,g23.16))') l,&
             ' dW_l=', dW_l(l),&
             ' flow=', flow(l),&
             ' div=', div(l)
     enddo
  endif

! Check for negative runoff:
  IF (lrunf_ie < lrunf_ie_min) THEN
     call get_current_point(ipt,jpt,kpt,fpt)
     write(*,*) 'note: at point ',ipt,jpt,kpt,fpt,'lrunf_ie=',lrunf_ie,' < lrunf_ie_min=',lrunf_ie_min
     call error_mesg(module_name, 'lrunf_ie < lrunf_ie_min', FATAL)
  ENDIF

  if(is_watch_point().or.(flag.and.write_when_flagged)) then
     write(*,*) ' ***** soil_step_2 checkpoint 3.3 ***** '
     write(*,*) 'psi_sat',soil%pars%psi_sat_ref
     write(*,*) 'Dpsi_max',Dpsi_max
     do l = 1, num_l
        write(*,'(i2.2,100(2x,a,g23.16))') l, &
             'Th=', (soil%ws(l) +soil%wl(l)+dW_l(l))/(dens_h2o*dz(l)), &
             'wl=', soil%wl(l)+dW_l(l), &
             'ws=', soil%ws(l), &
             'dW_l=', dW_l(l), &
             'dPsi=', dPsi(l), &
             'flow=', flow(l)
     enddo
  endif

end subroutine richards_clean

! ============================================================================
  subroutine RICHARDS(soil, psi, DThDP, hyd_cond, DKDP, div, &
                lprec_eff, Dpsi_min, Dpsi_max, dt_richards, &
                 dPsi, dW_l, flow, lrunf_ie)
! ZMS How to set prevent lrunf_ie when surface soil water < pond?

  type(soil_tile_type), intent(inout)   :: soil
  real, intent(in),  dimension(num_l)   :: psi, DThDP, hyd_cond, DKDP, div
  real, intent(in)                      :: lprec_eff, Dpsi_min, Dpsi_max
  real, intent(in)                      :: dt_richards
  real, intent(out), dimension(num_l)   :: dPsi, dW_l
  real, intent(out), dimension(num_l+1) :: flow
  real, intent(out)                     :: lrunf_ie
  ! ---- local vars ----------------------------------------------------------
  integer l, ipt, jpt, kpt, fpt, l_internal
  real, dimension(num_l-1) :: del_z, K, DKDPm, DKDPp, grad, eee, fff
  real aaa, bbb, ccc, ddd, xxx, dpsi_alt, dW_l_internal, w_to_move_up, adj
  logical flag

  flag = .false.
  flow(1) = dt_richards*lprec_eff
  do l = 1, num_l-1
     del_z(l) = zfull(l+1)-zfull(l)
     if (harmonic_mean_K) then
         K(l) = (dz(l)+dz(l+1))/(dz(l)/hyd_cond(l)+dz(l+1)/hyd_cond(l+1))
     else
         K(l) = 0.5*(hyd_cond(l)+hyd_cond(l+1))
     endif
     DKDPm(l) = 0. ! 0.5*DKDP(l)
     DKDPp(l) = 0. ! 0.5*DKDP(l+1)
!        K(l) = hyd_cond(l)
!        DKDPm(l) = DKDP(l)
!        DKDPp(l) = 0
      grad(l)  = (psi(l+1)-psi(l))/del_z(l) - 1
  enddo

  if (is_watch_point()) then
     write(*,*) '##### soil_step_2 checkpoint 3.1 #####'
     do l = 1, num_l
        write(*,'(x,a,x,i2.2,x,a,100(x,g23.16))') 'level=', l, 'DThDP,hyd_cond,psi,DKDP', &
             DThDP(l),&
             hyd_cond(l),&
             psi(l),&
             DKDP(l)
     enddo
     do l = 1, num_l-1
        write(*,'(a,i2.2,1x,a,100(2x,g23.16))') 'interface=', l, 'K,DKDPm,DKDPp,grad,del_z', &
             K(l),&
             DKDPm(l),&
             DKDPp(l),&
             grad(l)
     enddo
  endif


  l = num_l
  xxx = dens_h2o*dz(l)*DThDP(l)/dt_richards
  aaa =     - ( K(l-1)/del_z(l-1) - DKDPm(l-1)*grad(l-1))
  bbb = xxx - (-K(l-1)/del_z(l-1) - DKDPp(l-1)*grad(l-1) )
  ddd = - K(l-1) *grad(l-1) - div(l)
  eee(l-1) = -aaa/bbb
  fff(l-1) =  ddd/bbb

  if(is_watch_point()) then
       write(*,'(a,i2.2,100(2x,g23.16))') 'l,a,b, ,d', l,aaa, bbb,ddd
  endif

  do l = num_l-1, 2, -1
    xxx = dens_h2o*dz(l)*DThDP(l)/dt_richards
    aaa = - ( K(l-1)/del_z(l-1) - DKDPm(l-1)*grad(l-1))
    bbb = xxx-( -K(l-1)/del_z(l-1) - DKDPp(l-1)*grad(l-1)&
                -K(l  )/del_z(l  ) + DKDPm(l  )*grad(l  ))
    ccc =   - (  K(l  )/del_z(l  ) + DKDPp(l  )*grad(l  ))
    ddd =       K(l)*grad(l) - K(l-1)*grad(l-1) &
                          - div(l)
    eee(l-1) =                    -aaa/(bbb+ccc*eee(l))
    fff(l-1) =  (ddd-ccc*fff(l))/(bbb+ccc*eee(l))
    if(is_watch_point()) then
       write(*,'(a,i2.2,100(2x,g23.16))') 'l,a,b,c,d', l,aaa, bbb,ccc,ddd
    endif
  enddo

  l = 1
  xxx = dens_h2o*dz(l)*DThDP(l)/dt_richards
  bbb = xxx - ( -K(l  )/del_z(l  ) + DKDPm(l  )*grad(l  ))
  ccc =     - (  K(l  )/del_z(l  ) + DKDPp(l  )*grad(l  ))
  ddd =          flow(1)/dt_richards +    K(l)     *grad(l) &
                          - div(l)

  if (Dpsi_min.ge.Dpsi_max) call error_mesg(module_name, '=== Dpsi_min.ge.Dpsi_max', FATAL)

  IF (bbb+ccc*eee(l) .NE. 0.) THEN
     dPsi(l) = (ddd-ccc*fff(l))/(bbb+ccc*eee(l))
     if(is_watch_point()) then
        write(*,*) 'bbb+ccc*eee(l) .NE. 0.'
        write(*,*) 'bbb', bbb
        write(*,*) 'ccc', ccc
        write(*,*) 'ddd', ddd
        write(*,*) 'eee(l)', eee(l)
        write(*,*) 'fff(l)', fff(l)
        write(*,*) 'dPsi(l)', dPsi(l)
        write(*,*) 'dPsi(l)', dPsi(l)
        write(*,*) 'Dpsi_min', Dpsi_min
        write(*,*) 'Dpsi_max', Dpsi_max
     endif
     if (verbose.and.dPsi(l).lt.Dpsi_min) then
         call get_current_point(ipt,jpt,kpt,fpt)
         write(*,*) '=== warning: dPsi=',dPsi(l),'<min=',dPsi_min,'at',ipt,jpt,kpt,fpt
     endif
     if ((dPsi(l).gt.Dpsi_min.or.no_min_Dpsi) .and. dPsi(l).lt.Dpsi_max) then
        lrunf_ie = 0.
     else
        dPsi(l) = min (dPsi(l), Dpsi_max)
        if (dPsi(l).lt.Dpsi_min.and.(.not.no_min_Dpsi)) then
           flag = .true.
           dPsi(l) = Dpsi_min
        endif
        if (div_bug_fix) then
           flow(l) = (dPsi(l)*(bbb+ccc*eee(l))+ccc*fff(l)+div(l) &
                   - K(l)*grad(l))*dt_richards
        else
           flow(l) = (dPsi(l)*(bbb+ccc*eee(l))+ccc*fff(l) &
                   - K(l)*grad(l))*dt_richards
        endif
        lrunf_ie = lprec_eff - flow(l)/dt_richards
        if (lrunf_ie.lt.-lrunf_ie_tol) then
           if (verbose) then
              call get_current_point(ipt,jpt,kpt,fpt)
              write(*,*) '=== warning: rie= ',lrunf_ie,'<0 at',ipt,jpt,kpt,fpt
           endif
           if (.not.allow_negative_rie) then
              flag = .true.
              dpsi_alt = (ddd-ccc*fff(l))/(bbb+ccc*eee(l))
              if (verbose) write(*,*) '=== rie= ',lrunf_ie,'<0 reset to 0 and dPsi=',dPsi(l),' reset to ',dpsi_alt,'at ',ipt,jpt,kpt,fpt
              dPsi(l) = dpsi_alt
              lrunf_ie = 0.
              flow(l) = lprec_eff*dt_richards
           endif
        endif
     endif
  ELSE
     call error_mesg(module_name, 'b+ce=0 in soil-water equations', FATAL)
     if (verbose) then
       call get_current_point(ipt,jpt,kpt,fpt)
       write(*,*) '===richards b+ce=0 ===','at point ',ipt,jpt,kpt,fpt
       endif
     if(is_watch_point()) then
        write(*,*) 'bbb+ccc*eee(l) .EQ. 0.'
        write(*,*) 'bbb', bbb
        write(*,*) 'ccc', ccc
        write(*,*) 'ddd', ddd
        write(*,*) 'eee(l)', eee(l)
        write(*,*) 'fff(l)', fff(l)
        write(*,*) 'dPsi(l)', dPsi(l)
        write(*,*) 'dPsi(l)', dPsi(l)
        write(*,*) 'Dpsi_min', Dpsi_min
        write(*,*) 'Dpsi_max', Dpsi_max
     endif
     dPsi(l) = Dpsi_max
     flow(l) = (dPsi(l)*(bbb+ccc*eee(l))+ccc*fff(l) &
                   - K(l)*grad(l))*dt_richards
     lrunf_ie = lprec_eff - flow(l)/dt_richards
     if (lrunf_ie.lt.-lrunf_ie_tol) then
        if (verbose) then
           call get_current_point(ipt,jpt,kpt,fpt)
           write(*,*) '===richards b+ce=0 AND lrunf_ie.lt.-lrunf_ie_tol at point ',ipt,jpt,kpt,fpt
           write(*,*) '===richards b+ce=0 AND lrunf_ie.lt.-lrunf_ie_tol at point ',ipt,jpt,kpt,fpt,'rie=',lrunf_ie
        endif
        if(.not.allow_negative_rie) then
           flag = .true.
           dpsi_alt = 0.
           if (verbose) then
             write(*,*) '===richards b+ce=0 AND lrunf_ie.lt.-lrunf_ie_tol at point ',ipt,jpt,kpt,fpt,'rie=',lrunf_ie,' reset to 0'
             write(*,*) '===richards b+ce=0 AND lrunf_ie.lt.-lrunf_ie_tol at point ',ipt,jpt,kpt,fpt,'dPsi=',dPsi(l),' reset to',dpsi_alt
             write(*,*) '===richards b+ce=0 AND lrunf_ie.lt.-lrunf_ie_tol at point ',ipt,jpt,kpt,fpt,'dPsi_min/max=',dPsi_min,dPsi_max
           endif
           dPsi(l) = dpsi_alt
           lrunf_ie = 0.
           flow(l) = lprec_eff*dt_richards
        endif
     endif
  ENDIF

  if(is_watch_point().or.(flag.and.write_when_flagged)) then
     write(*,'(a,i2.2,100(2x,g23.16))') 'l,  b,c,d', l, bbb,ccc,ddd
     write(*,*) ' ##### soil_step_2 checkpoint 3.2 #####'
     write(*,*) 'ie:', lrunf_ie
     do l = 1, num_l-1
        write(*,'(a,i2.2,100(2x,g23.16))') 'l,eee(l),fff(l)',l,eee(l),fff(l)
     enddo
     write(*,*) 'DThDP(1)', DThDP(1)
     write(*,*) 'K(1)', K(1)
     write(*,*) 'grad(1)', grad(1)
     write(*,*) 'ddd(1)', ddd
     write(*,*) 'ccc(1)', ccc
     write(*,*) 'bbb(1)', bbb
     write(*,*) 'dPsi(1)', dPsi(1)
     write(*,*) 'Psi(1)', Psi(1)
     write(*,*) 'div(1)', div(1)
  endif

  do l = 1, num_l-1
     dPsi(l+1) = eee(l)*dPsi(l) + fff(l)
     flow(l+1) = dt_richards*( &
         -K(l)*(grad(l)&
         +(DPsi(l+1)-DPsi(l))/ del_z(l)) &
         -grad(l)*(DKDPp(l)*Dpsi(l+1)+ &
                         DKDPm(l)*Dpsi(l) )  )
     dW_l(l) = flow(l) - flow(l+1) - div(l)*dt_richards
  enddo
  flow(num_l+1) = 0.
  dW_l(num_l) = flow(num_l) - flow(num_l+1) &
                          - div(num_l)*dt_richards

  if(is_watch_point().or.(flag.and.write_when_flagged)) then
     write(*,*) ' ##### soil_step_2 checkpoint 3.21 #####'
     do l = 1, num_l
        write(*,'(i2.2,100(2x,a,g23.16))') l,&
             ' dW_l=', dW_l(l),&
             ' flow=', flow(l),&
             ' div=', div(l)
     enddo
  endif

  if (flag) then
     l_internal = 1
     dW_l_internal = -1.e20
     do l = 2, num_l
        if (dW_l(l).gt.dW_l_internal) then
           l_internal = l
           dW_l_internal = dW_l(l)
        endif
     enddo
     w_to_move_up = min(dW_l_internal, -(soil%wl(1)+dW_l(1)))
     w_to_move_up = max(w_to_move_up, 0.)
     write(*,*) 'l_internal=',l_internal
     write(*,*) 'dW_l(l_internal)=',dW_l(l_internal)
     write(*,*) 'soil%wl(1)+dW_l(1)',soil%wl(1)+dW_l(1)
     write(*,*) 'w_to_move_up=',w_to_move_up
     if (l_internal.gt.1) then
        dW_l(1) = dW_l(1) + w_to_move_up
        dW_l(l_internal) = dW_l(l_internal) - w_to_move_up
        do l = 2, l_internal
           flow(l) = flow(l) - w_to_move_up
        enddo
     endif
  endif

  if(is_watch_point().or.(flag.and.write_when_flagged)) then
     write(*,*) ' ##### soil_step_2 checkpoint 3.22 #####'
     do l = 1, num_l
        write(*,'(i2.2,100(2x,a,g23.16))') l,&
             ' dW_l=', dW_l(l),&
             ' flow=', flow(l),&
             ' div=', div(l)
     enddo
  endif

! In rare situations where lrunf_ie is large and negative, clip any liquid supersaturation
! layer by layer and recompute lrunf_ie (this is not good, since it ignores 'comp'):
  IF (lrunf_ie < lrunf_ie_min) THEN
     call get_current_point(ipt,jpt,kpt,fpt)
     write(*,*) 'note: at point ',ipt,jpt,kpt,fpt,' clip triggered by lrunf_ie=',lrunf_ie
     call error_mesg(module_name, 'lrunf_ie < lrunf_ie_min', FATAL)
     do l = num_l, 1, -1
        adj = max(dW_l(l)+soil%ws(l)+soil%wl(l) &
             - soil%pars%vwc_sat*dz(l)*dens_h2o, 0. )

        if(is_watch_point()) then
           write(*,*) '3.22 l=', l,&
                ' soil%wl=',soil%wl(l),  &
                ' soil%ws=',soil%ws(l) , &
                ' soil%pars%vwc_sat=', soil%pars%vwc_sat, &
                ' dz=', dz(l), &
                ' adj=', adj
        endif

        adj = min(adj, max(0.,soil%wl(l)))

        if(is_watch_point()) then
           write(*,*) '3.23 l=', l, ' adj=', adj
        endif

        dW_l(l) = dW_l(l) - adj
        flow(l) = flow(l+1) + dW_l(l) + div(l)*dt_richards
     enddo
     lrunf_ie = lprec_eff - flow(1)/dt_richards

  ENDIF

  if(is_watch_point().or.(flag.and.write_when_flagged)) then
     write(*,*) ' ***** soil_step_2 checkpoint 3.3 ***** '
     write(*,*) 'psi_sat',soil%pars%psi_sat_ref
     write(*,*) 'Dpsi_max',Dpsi_max
     do l = 1, num_l
        write(*,'(i2.2,100(2x,a,g23.16))') l, &
             'Th=', (soil%ws(l) +soil%wl(l)+dW_l(l))/(dens_h2o*dz(l)), &
             'wl=', soil%wl(l)+dW_l(l), &
             'ws=', soil%ws(l), &
             'dW_l=', dW_l(l), &
             'dPsi=', dPsi(l), &
             'flow=', flow(l)
     enddo
  endif

end subroutine richards

! ============================================================================
  subroutine advection(soil, flow, dW_l, tflow, d_GW, div, delta_time)
  type(soil_tile_type), intent(inout) :: soil
  real, intent(in), dimension(:) :: flow  ! water tendency downwards into layer [mm]
  real, intent(in), dimension(:) :: dW_l  ! net water tendency in layer [mm]
  real, intent(in), dimension(:) :: div   ! horizontal water flux divergence [mm/s]
  real, intent(in)               :: delta_time ! land model timestep [s]
  real, intent(in) :: &
     tflow, & ! temperature of surface downwards flow [K]
     d_GW
  ! ---- local vars ----------------------------------------------------------
  real, dimension(num_l)   :: u_minus, u_plus, del_t
  real, dimension(num_l-1) :: eee, fff
  real hcap, aaa, bbb, ccc
  integer l
  ! For energy conservation
  real :: esum1, esum2 ! [W/m^2] heat content of soil before and after solution
  real, parameter :: ethresh = 1.e-4 ! [W/m^2] Allowable error in energy solution for roundoff

!  if (do_component_balchecks .and. .not. LM2) then

     ! Sum energy content in soil before solution
     esum1 = clw*max(flow(1), 0.)*(tflow-tfreeze) ! initialize to incoming surface energy tendency
     do l = 1, num_l
        esum1 = esum1 + (soil%heat_capacity_dry(l)*dz(l) + csw*soil%ws(l) + &
                           clw*(soil%wl(l) - dW_l(l)) ) * (soil%T(l)-tfreeze)
                           ! use water content before Richards eq. update to be consistent with
                           ! heat advection solution
     end do

!  end if

! Upstream weighting of advection. Preserving u_plus here for now.
  u_minus = 1.
  where (flow(1:num_l).lt.0.) u_minus = 0.
  do l = 1, num_l-1
     u_plus(l) = 1. - u_minus(l+1)
  enddo
  hcap = (soil%heat_capacity_dry(num_l)*dz(num_l) &
                              + csw*soil%ws(num_l))/clw
  aaa = -flow(num_l) * u_minus(num_l)
  bbb =  hcap + soil%wl(num_l) - dW_l(num_l) - aaa
  eee(num_l-1) = -aaa/bbb
  fff(num_l-1) = aaa*(soil%T(num_l)-soil%T(num_l-1)) / bbb

  do l = num_l-1, 2, -1
    hcap = (soil%heat_capacity_dry(l)*dz(l) &
                              + csw*soil%ws(l))/clw
    aaa = -flow(l)   * u_minus(l)
    ccc =  flow(l+1) * u_plus (l)
    bbb =  hcap + soil%wl(l) - dW_l(l) - aaa - ccc
    eee(l-1) = -aaa / ( bbb +ccc*eee(l) )
    fff(l-1) = (   aaa*(soil%T(l)-soil%T(l-1))    &
                       + ccc*(soil%T(l)-soil%T(l+1))    &
                       - ccc*fff(l) ) / ( bbb +ccc*eee(l) )
  enddo

  hcap = (soil%heat_capacity_dry(1)*dz(1) + csw*soil%ws(1))/clw
  aaa = -flow(1) * u_minus(1)
  ccc =  flow(2) * u_plus (1)
  bbb =  hcap + soil%wl(1) - dW_l(1) - aaa - ccc

  del_t(1) =  (  aaa*(soil%T(1)-tflow          ) &
                     + ccc*(soil%T(1)-soil%T(2)) &
                     - ccc*fff(1) ) / (bbb+ccc*eee(1))
  soil%T(1) = soil%T(1) + del_t(1)

  if(is_watch_point()) then
     write(*,*) ' ***** soil_step_2 checkpoint 3.4.1 ***** '
     write(*,*) 'hcap', hcap
     write(*,*) 'aaa', aaa
     write(*,*) 'bbb', bbb
     write(*,*) 'ccc', ccc
     write(*,*) 'del_t(1)', del_t(1)
     write(*,*) ' T(1)', soil%T(1)
  endif

  do l = 1, num_l-1
     del_t(l+1) = eee(l)*del_t(l) + fff(l)
     soil%T(l+1) = soil%T(l+1) + del_t(l+1)
  enddo

!  if (do_component_balchecks .and. .not. LM2) then
     ! Sum energy content in soil after solution
     esum2 = -clw * min(flow(1), 0.) * (soil%T(1) - tfreeze)
     if (flow(1) < 0. .and. is_watch_point()) then
        write(*,*) 'Apparent reflux in advection'
        write(*,*) esum2, ' J/m^2 energy refluxed.'
     end if
     do l = 1, num_l
        esum2 = esum2 + (soil%heat_capacity_dry(l)*dz(l) + csw*soil%ws(l) + clw*soil%wl(l)) * &
                      (soil%T(l)-tfreeze) ! no phase change here
        ! Add energy lost to stream via divergence
        esum2 = esum2 + clw*div(l) * delta_time * (soil%T(l)-tfreeze)
     end do

     call check_conservation(module_name // ': Advection subroutine.', 'Energy', esum1, esum2, ethresh, &
                             FATAL)

!  end if

  ! (lumped=lm2 groundwater stored in l=1 prog variable, liquid only)
  if (soil%groundwater(1).ne. 0.) soil%groundwater_T(1) =    &
       + ((aquifer_heat_cap+soil%groundwater(1)-d_GW)  &
                                 *soil%groundwater_T(1) &
        + flow(num_l+1)*soil%T(num_l)) &
         /((aquifer_heat_cap+soil%groundwater(1)-d_GW) + flow(num_l+1))

end subroutine advection


! ============================================================================
! Thermal solution for advection of heat by soil water flow,
! using a call to the generic tridiagonal solver.
! Tested in diagnostic mode; required for gw_option == GW_TILED.
subroutine advection_tri(soil, flow, dW_l, tflow, d_GW, div, delta_time, t_soil_tridiag, hdiv_it)
   type(soil_tile_type), intent(inout) :: soil
   real, intent(in), dimension(:) :: flow  ! water tendency downwards into layer [mm]
   real, intent(in), dimension(:) :: dW_l  ! net water tendency in layer [mm]
   real, intent(in), dimension(:) :: div   ! horizontal water flux divergence [mm/s]
                                           ! used with simple or integrated groundwater model
   real, intent(in)               :: delta_time ! land model timestep [s]
   real, intent(out), dimension(:) :: t_soil_tridiag ! soil temperature solution [K]
   real, intent(in) :: &
      tflow, & ! temperature of surface downwards flow [K]
      d_GW     ! groundwater tendency (for LM2) [mm ??]
   real, intent(in), dimension(:) :: hdiv_it ! divergence of heat due to inter-tile water flow [W/m^2]
   ! ---- local vars ----------------------------------------------------------
   real, dimension(num_l)   :: u_minus, &! coefficient for flow into layer from above [-]
                              u_plus, & ! coefficient for flow into layer from below [-] (=1-u_minus for layer below)
                              del_t, &  ! layer T tendency [K]
                              aaa, bbb, ccc, ddd  ! Tridiagonal matrix coefficients on del_t
   real :: hcapr ! heat capacity ratio of layer, c_solid / c_l * Delta z [m]
   integer :: l ! layer index
   ! For energy conservation
   real :: esum1, esum2 ! [W/m^2] heat content of soil before and after solution
   real, parameter :: ethresh = 1.e-4 ! [W/m^2] Allowable error in energy solution for roundoff

!   if (do_component_balchecks) then
      esum1 = clw*max(flow(1), 0.)*(tflow-tfreeze) ! initialize to incoming surface energy tendency
      do l = 1, num_l
         esum1 = esum1 + (soil%heat_capacity_dry(l)*dz(l) + csw*soil%ws(l) + &
                           clw*(soil%wl(l) - dW_l(l)) ) * (soil%T(l)-tfreeze)
                           ! use water content before Richards eq. update to be consistent with
                           ! heat advection solution
         if (is_watch_point()) write(*,*)'l, esum1:', l, esum1
      end do
!   end if

   ! Upstream weighting of advection. Preserving u_plus here for now.
   u_minus = 1.
   where (flow.lt.0.) u_minus = 0.
   do l = 1, num_l-1
      u_plus(l) = 1. - u_minus(l+1)
   enddo

   ! Set Tridiagonal coefficients
   ! Solving equations aaa(l)*delta_t(l-1) + bbb(l)*delta_t(l) + ccc(l)*delta_t(l+1) = ddd(l)
   do l = 1, num_l
      hcapr = (soil%heat_capacity_dry(l)*dz(l) &
                              + csw*soil%ws(l))/clw
      if (l == 1) then
         aaa(l) = 0. ! T_0 is prescribed to tflow
      else
         aaa(l) = - u_minus(l) * flow(l)
      end if

      bbb(l) = hcapr + soil%wl(l) - (1.-u_minus(l)) * flow(l) + &
            (1.-u_plus(l)) * flow(l+1) + delta_time * div(l)

      if (l < num_l) then
         ccc(l) = u_plus(l) * flow(l+1)
      else
         ccc(l) = 0. ! zero flux bottom boundary; flow(l+1) = 0.
      end if

      if (l == 1) then
         ddd(l) = u_minus(l) * flow(l) * tflow  &
                  + soil%T(l) * ( (1.-u_minus(l)) * flow(l)  &
                     - (1-u_plus(l)) * flow(l+1) - dW_l(l) - delta_time * div(l) ) &
                  - u_plus(l) * flow(l+1) * soil%T(l+1)
      else if (l < num_l) then
         ddd(l) = u_minus(l) * flow(l) * soil%T(l-1)  &
                  + soil%T(l) * ( (1.-u_minus(l)) * flow(l)  &
                     - (1-u_plus(l)) * flow(l+1) - dW_l(l) - delta_time * div(l) ) &
                  - u_plus(l) * flow(l+1) * soil%T(l+1)
      else ! l == num_l
         ddd(l) = u_minus(l) * flow(l) * soil%T(l-1)  &
                  + soil%T(l) * ( (1.-u_minus(l)) * flow(l)  &
                      - dW_l(l) - delta_time * div(l) )
      end if

      ! update coefficients if full tiled hillslope model is used
      if (gw_option == GW_TILED) then
         ! Remove div terms. Energy flux will be prescribed explicitly based on inter-tile flows
         ! calculated in hlsp_hydrology_1.
         bbb(l) = bbb(l) - delta_time * div(l)
         ddd(l) = ddd(l) + soil%T(l) * delta_time * div(l)
         ! Add explicit energy sink term to ddd(l)
         ddd(l) = ddd(l) - delta_time / clw * hdiv_it(l) - delta_time * div(l) * tfreeze
         ! kg/m^2 K          s        /(J/kg/K) * J/m^2/s      s     *  kg/m^2/s   *K
      end if
   end do

   ! Update temperature
   call tridiag(aaa, bbb, ccc, ddd, del_t)
   t_soil_tridiag(1:num_l) = soil%T(1:num_l) + del_t(1:num_l)

   if (use_tridiag_foradvec) then
      soil%T(1:num_l) = t_soil_tridiag(1:num_l)
   end if

   if(is_watch_point()) then
      write(*,*) ' ***** soil_step_2 checkpoint 3.4.1 ***** '
      write(*,*) 'hcapr(num_l)', hcapr
      write(*,*) 'aaa', aaa(:)
      write(*,*) 'bbb', bbb(:)
      write(*,*) 'ccc', ccc(:)
      write(*,*) 'ddd', ddd(:)
      write(*,*) ' T', soil%T(:)
      write(*,*) 'hdiv_it', hdiv_it(:)
   endif

!   if (do_component_balchecks) then
      ! Sum energy content in soil after solution
      esum2 = -clw * min(flow(1), 0.) * (t_soil_tridiag(1) - tfreeze) ! if upwards flow at surface
                                      ! will be energy lost to sat runoff ?
      if (flow(1) < 0. .and. is_watch_point()) then
         write(*,*) 'Apparent reflux in advection'
         write(*,*) esum2, ' J/m^2 energy refluxed.'
      end if
      do l = 1, num_l
         esum2 = esum2 + (soil%heat_capacity_dry(l)*dz(l) + csw*soil%ws(l) + clw*soil%wl(l)) * &
                        (t_soil_tridiag(l)-tfreeze) ! no phase change here
         ! Add energy lost via divergence
         if (gw_option /= GW_TILED) then
            esum2 = esum2 + clw*div(l) * delta_time * (t_soil_tridiag(l)-tfreeze)
         else
            esum2 = esum2 + hdiv_it(l) * delta_time
         end if
         if (is_watch_point()) write(*,*)'l, esum2:', l, esum2
      end do

      call check_conservation(module_name // ': Advection_tri subroutine.', 'Energy', esum1, esum2, ethresh, &
                              FATAL)

!   end if


   ! (lumped=lm2 groundwater stored in l=1 prog variable, liquid only)
   if (soil%groundwater(1).ne. 0.) soil%groundwater_T(1) =    &
         + ((aquifer_heat_cap+soil%groundwater(1)-d_GW)  &
                            *soil%groundwater_T(1) &
         + flow(num_l+1)*soil%T(num_l)) &
         /((aquifer_heat_cap+soil%groundwater(1)-d_GW) + flow(num_l+1))

end subroutine advection_tri


! ============================================================================
! Spread new root C through profile, using vertical root profile from vegn_uptake_profile
subroutine add_root_litter(soil,vegn,newlitterC,newlitterN)
    type(soil_tile_type), intent(inout)  :: soil
    type(vegn_tile_type), intent(in)     :: vegn
    real,intent(in) :: newlitterC(:), newlitterN(:)

    real,dimension(num_l) :: uptake_frac_max, vegn_uptake_term, vrl, rhizosphere_frac
    real :: K_r,r_r
    real :: plant_height, xylem_area_frac,xylem_resist,dum4
    integer :: nn

    ! This is just to set the root profile using uptake_frac_max
    call vegn_uptake_profile (vegn, dz(1:num_l), uptake_frac_max, vegn_uptake_term )
    ! r_r and vpl are to set rhizosphere size
    call vegn_hydraulic_properties(vegn, dz(1:num_l),always_use_bsw, vrl, K_r, r_r, plant_height, xylem_area_frac,xylem_resist,dum4)
    if(abs(sum(uptake_frac_max(1:num_l)) - 1.0) > 1e-10 ) then
        print *,'1 - sum(vegn_uptake_frac_max)',1.0-sum(uptake_frac_max(1:num_l))
        call error_mesg('add_root_litter','total of vegn_uptake_frac_max not 1',FATAL)
    endif
    do nn=1,num_l
        !Note, rhizosphere_frac will be ignored if use_rhizosphere_cohort is FALSE in soil_carbon
        rhizosphere_frac(nn)=min(3.141592*((r_rhiz+r_r)**2-r_r**2)*vrl(nn),1.0)
        call add_litter(soil%soil_organic_matter(nn),newLitterC*uptake_frac_max(nn),newLitterN*uptake_frac_max(nn),rhizosphere_frac=rhizosphere_frac(nn))
        soil%fsc_in(nn)=soil%fsc_in(nn)+newLitterC(1)*uptake_frac_max(nn)
        soil%ssc_in(nn)=soil%ssc_in(nn)+newLitterC(2)*uptake_frac_max(nn)
        soil%fsn_in(nn)=soil%fsn_in(nn)+newLitterN(1)*uptake_frac_max(nn)
        soil%ssn_in(nn)=soil%ssn_in(nn)+newLitterN(2)*uptake_frac_max(nn)
    enddo
end subroutine add_root_litter



! ============================================================================
! Spread root exudate C through profile, using vertical root profile from vegn_uptake_profile
! Differs from add_root_litter -- C is distributed through existing cohorts, not deposited as new cohort
subroutine add_root_exudates(soil,vegn,exudateC,exudateN)
    type(soil_tile_type), intent(inout)  :: soil
    type(vegn_tile_type), intent(in)     :: vegn
    real,intent(in) :: exudateC,exudateN

    real,dimension(num_l) :: uptake_frac_max, vegn_uptake_term
    integer :: nn

    call vegn_uptake_profile (vegn, dz(1:num_l), uptake_frac_max, vegn_uptake_term )
    if(abs(sum(uptake_frac_max(1:num_l)) - 1.0) > 1e-10 ) then
        print *,'1 - sum(vegn_uptake_frac_max)',1.0-sum(uptake_frac_max(1:num_l))
        call error_mesg('add_root_litter','total of vegn_uptake_frac_max not 1',FATAL)
    endif

    if (is_watch_point()) then
       write (*,*) '##### add_root_exudates #####'
       __DEBUG1__(exudateC)
       __DEBUG1__(uptake_frac_max)
    endif

    do nn=1,num_l
        if (is_watch_point()) then
           call debug_pool(soil%soil_organic_matter(nn),'soil_organic_matter(nn) before')
        endif
        call add_C_N_to_rhizosphere(soil%soil_organic_matter(nn),(/exudateC*uptake_frac_max(nn),0.0,0.0/),(/exudateN*uptake_frac_max(nn),0.0,0.0/))
        soil%fsc_in(nn)=soil%fsc_in(nn)+exudateC*uptake_frac_max(nn)
        soil%fsn_in(nn)=soil%fsn_in(nn)+exudateN*uptake_frac_max(nn)
        if (is_watch_point()) then
           call debug_pool(soil%soil_organic_matter(nn),'soil_organic_matter(nn) after ')
        endif
    enddo


end subroutine add_root_exudates


! Nitrogen uptake from the rhizosphere by roots (active transport across root-soil interface)
! Mineral nitrogen is taken up from the rhizosphere only
subroutine root_N_uptake(soil,vegn,total_N_uptake,dt)
  real,intent(out)::total_N_uptake
  type(vegn_tile_type),intent(in)::vegn
  type(soil_tile_type),intent(inout)::soil
  real,intent(in)::dt

  real,dimension(num_l) :: uptake_frac_max, vegn_uptake_term, vrl
  real::nitrate_uptake,ammonium_uptake
  integer::nn

  real :: K_r,r_r
  real :: plant_height, xylem_area_frac,xylem_resist,dum4
  real :: rhizosphere_frac

  call vegn_uptake_profile (vegn, dz(1:num_l), uptake_frac_max, vegn_uptake_term )
  ! r_r and vpl are to set rhizosphere size
  call vegn_hydraulic_properties(vegn, dz(1:num_l),always_use_bsw, vrl, K_r, r_r, plant_height, xylem_area_frac,xylem_resist,dum4)

  total_N_uptake=0.0
  do nn=1,num_l
    rhizosphere_frac=min(3.141592*((r_rhiz+r_r)**2-r_r**2)*vrl(nn),1.0)
    ammonium_uptake = soil%soil_organic_matter(nn)%ammonium*rhizosphere_frac*root_NH4_uptake_rate*dt
    nitrate_uptake = soil%soil_organic_matter(nn)%nitrate*rhizosphere_frac*root_NO3_uptake_rate*dt
    soil%soil_organic_matter(nn)%ammonium=soil%soil_organic_matter(nn)%ammonium-ammonium_uptake
    soil%soil_organic_matter(nn)%nitrate=soil%soil_organic_matter(nn)%nitrate-nitrate_uptake
    total_N_uptake = total_N_uptake + ammonium_uptake + nitrate_uptake
  enddo


end subroutine root_N_uptake


! Uptake of mineral N by mycorrhizal "scavengers" -- Should correspond to Arbuscular mycorrhizae
subroutine myc_scavenger_N_uptake(soil,vegn,myc_biomass,total_N_uptake,dt)
  real,intent(in)::myc_biomass,dt  ! dt in years, myc_biomass in kgC/m2
  real,intent(out)::total_N_uptake
  type(vegn_tile_type),intent(in)::vegn
  type(soil_tile_type),intent(inout)::soil

  real,dimension(num_l) :: uptake_frac_max, vegn_uptake_term
  real::nitrate_uptake,ammonium_uptake
  integer::nn

  call vegn_uptake_profile (vegn, dz(1:num_l), uptake_frac_max, vegn_uptake_term )

  total_N_uptake=0.0
  do nn=1,num_l
    call mycorrhizal_mineral_N_uptake_rate(soil%soil_organic_matter(nn),myc_biomass*uptake_frac_max(nn),dz(nn),nitrate_uptake,ammonium_uptake)
    total_N_uptake=total_N_uptake+(ammonium_uptake+nitrate_uptake)*dt
    soil%soil_organic_matter(nn)%ammonium=soil%soil_organic_matter(nn)%ammonium-ammonium_uptake*dt
    soil%soil_organic_matter(nn)%nitrate=soil%soil_organic_matter(nn)%nitrate-nitrate_uptake*dt
  enddo

  ! Mycorrhizae should have access to litter layer too
  ! Might want to update this so it calculates actual layer thickness?
  call mycorrhizal_mineral_N_uptake_rate(soil%leafLitter,myc_biomass*uptake_frac_max(1),dz(1),nitrate_uptake,ammonium_uptake)
  total_N_uptake=total_N_uptake+(nitrate_uptake+ammonium_uptake)*dt
  soil%leafLitter%ammonium=soil%leafLitter%ammonium-ammonium_uptake*dt
  soil%leafLitter%nitrate=soil%leafLitter%nitrate-nitrate_uptake*dt

end subroutine myc_scavenger_N_uptake

! Uptake of mineral N by mycorrhizal "scavengers" -- Should correspond to Arbuscular mycorrhizae
! This version is for calculating marginal N acquisition benefit of change in myc biomass, and does not affect pools
pure subroutine hypothetical_myc_scavenger_N_uptake(soil,vegn,myc_biomass,total_N_uptake,dt)
  real,intent(in)::myc_biomass,dt  ! dt in years, myc_biomass in kgC/m2
  real,intent(out)::total_N_uptake
  type(vegn_tile_type),intent(in)::vegn
  type(soil_tile_type),intent(in)::soil

  real,dimension(num_l) :: uptake_frac_max, vegn_uptake_term
  real::nitrate_uptake,ammonium_uptake
  integer::nn

  call vegn_uptake_profile (vegn, dz(1:num_l), uptake_frac_max, vegn_uptake_term )

  total_N_uptake=0.0
  do nn=1,num_l
    call mycorrhizal_mineral_N_uptake_rate(soil%soil_organic_matter(nn),myc_biomass*uptake_frac_max(nn),dz(nn),nitrate_uptake,ammonium_uptake)
    total_N_uptake=total_N_uptake+(ammonium_uptake+nitrate_uptake)*dt
    !soil%soil_organic_matter(nn)%ammonium=soil%soil_organic_matter(nn)%ammonium-ammonium_uptake*dt
    !soil%soil_organic_matter(nn)%nitrate=soil%soil_organic_matter(nn)%nitrate-nitrate_uptake*dt
  enddo

  ! Mycorrhizae should have access to litter layer too
  ! Might want to update this so it calculates actual layer thickness?
  call mycorrhizal_mineral_N_uptake_rate(soil%leafLitter,myc_biomass*uptake_frac_max(1),dz(1),nitrate_uptake,ammonium_uptake)
  total_N_uptake=total_N_uptake+(nitrate_uptake+ammonium_uptake)*dt
  !soil%leafLitter%ammonium=soil%leafLitter%ammonium-ammonium_uptake*dt
  !soil%leafLitter%nitrate=soil%leafLitter%nitrate-nitrate_uptake*dt

end subroutine hypothetical_myc_scavenger_N_uptake


! Uptake of mineral N by mycorrhizal "miners" -- Should correspond to Ecto mycorrhizae
subroutine myc_miner_N_uptake(soil,vegn,myc_biomass,total_N_uptake,total_C_uptake,total_CO2prod,dt)
  real,intent(in)::myc_biomass,dt  ! dt in years, myc_biomass in kgC/m2
  real,intent(out)::total_N_uptake,total_C_uptake,total_CO2prod
  type(vegn_tile_type),intent(in)::vegn
  type(soil_tile_type),intent(inout)::soil

  real,dimension(num_l) :: uptake_frac_max, vegn_uptake_term
  real::N_uptake,C_uptake,CO2prod
  real,dimension(num_l)::T,theta,air_filled_porosity
  integer::nn

  call vegn_uptake_profile (vegn, dz(1:num_l), uptake_frac_max, vegn_uptake_term )
  T = soil%T(:)
  theta = max(min(soil_theta(soil),1.0),0.0)
  air_filled_porosity=max(min(1.0-theta-soil_ice_porosity(soil),1.0),0.0)

  total_N_uptake=0.0
  total_C_uptake=0.0
  total_CO2prod=0.0
  do nn=1,num_l
    call mycorrhizal_decomposition(soil%soil_organic_matter(nn),myc_biomass*uptake_frac_max(nn),T(nn),theta(nn),air_filled_porosity(nn),N_uptake,C_uptake,CO2prod,dt)
    total_N_uptake=total_N_uptake+N_uptake
    total_C_uptake = total_C_uptake+C_uptake
    total_CO2prod=total_CO2prod+CO2prod
  enddo

  ! Mycorrhizae should have access to litter layer too
  ! Might want to update this so it calculates actual layer thickness?
  call mycorrhizal_decomposition(soil%leafLitter,myc_biomass*uptake_frac_max(1),T(1),theta(1),air_filled_porosity(1),N_uptake,C_uptake,CO2prod,dt)
  total_N_uptake=total_N_uptake+N_uptake
  total_C_uptake = total_C_uptake+C_uptake
  total_CO2prod = total_CO2prod + CO2prod

end subroutine myc_miner_N_uptake


! Uptake of mineral N by mycorrhizal "miners" -- Should correspond to Ecto mycorrhizae
! Pure version, for use in marginal gain calculations
pure subroutine hypothetical_myc_miner_N_uptake(soil,vegn,myc_biomass,total_N_uptake,total_C_uptake,total_CO2prod,dt)
  real,intent(in)::myc_biomass,dt  ! dt in years, myc_biomass in kgC/m2
  real,intent(out)::total_N_uptake,total_C_uptake,total_CO2prod
  type(vegn_tile_type),intent(in)::vegn
  type(soil_tile_type),intent(inout)::soil

  real,dimension(num_l) :: uptake_frac_max, vegn_uptake_term,T,theta,air_filled_porosity
  real::N_uptake,C_uptake,CO2prod
  integer::nn

  call vegn_uptake_profile (vegn, dz(1:num_l), uptake_frac_max, vegn_uptake_term )
  T = soil%T(:)
  theta = max(min(soil_theta(soil),1.0),0.0)
  air_filled_porosity=max(min(1.0-theta-soil_ice_porosity(soil),1.0),0.0)

  total_CO2prod=0.0
  total_N_uptake=0.0
  total_C_uptake=0.0
  do nn=1,num_l
    call hypothetical_mycorrhizal_decomposition(soil%soil_organic_matter(nn),myc_biomass*uptake_frac_max(nn),T(nn),theta(nn),air_filled_porosity(nn),N_uptake,C_uptake,CO2prod,dt)
    total_N_uptake=total_N_uptake+N_uptake
    total_C_uptake = total_C_uptake+C_uptake
    total_CO2prod=total_CO2prod+CO2prod
  enddo

  ! Mycorrhizae should have access to litter layer too
  ! Might want to update this so it calculates actual layer thickness?
  call hypothetical_mycorrhizal_decomposition(soil%leafLitter,myc_biomass*uptake_frac_max(1),T(1),theta(1),air_filled_porosity(1),N_uptake,C_uptake,CO2prod,dt)
  total_N_uptake=total_N_uptake+N_uptake
  total_C_uptake = total_C_uptake+C_uptake
  total_CO2prod=total_CO2prod+CO2prod

end subroutine hypothetical_myc_miner_N_uptake


! ============================================================================
subroutine redistribute_peat_carbon(soil)
    type(soil_tile_type), intent(inout) :: soil

    integer :: nn
    real :: layer_total_C,layer_total_C_2,layer_max_C,layer_extra_C,fraction_to_remove
    real :: layer_total_N,layer_total_N_2,layer_max_N,layer_extra_N
    real :: total_C_before,total_C_after
    real :: leaflitter_total_C, woodlitter_total_C
    real :: total_N_before,total_N_after
    real :: leaflitter_total_N, woodlitter_total_N

    !For conservation check.
    total_C_before=0.0
    total_N_before=0.0
    do nn=1,num_l
    call poolTotals(soil%soil_organic_matter(num_l),totalCarbon=layer_total_C,totalNitrogen=layer_total_N)
    total_C_before=total_C_before+layer_total_C
    total_N_before=total_N_before+layer_total_N
    enddo

    call poolTotals(soil%leaflitter,totalCarbon=leaflitter_total_C,totalNitrogen=leaflitter_total_N)
    call poolTotals(soil%coarseWoodLitter,totalCarbon=woodlitter_total_C,totalNitrogen=woodlitter_total_N)
    layer_total_C=leaflitter_total_C+woodlitter_total_C
    layer_total_N=leaflitter_total_N+woodlitter_total_N

    total_C_before=total_C_before+layer_total_C
    total_N_before=total_N_before+layer_total_N

    layer_max_C=max_litter_thickness*max_soil_C_density
    layer_extra_C = layer_total_C-layer_max_C
    if(layer_extra_C>0) then
        fraction_to_remove=1.0-layer_max_C/layer_total_C
        call transfer_pool_fraction(soil%leaflitter,soil%soil_organic_matter(1),fraction_to_remove)
        call transfer_pool_fraction(soil%coarsewoodlitter,soil%soil_organic_matter(1),fraction_to_remove)
    endif

    !Move carbon down if it exceeds layer_max_C
    do nn=1,num_l-1
        call poolTotals(soil%soil_organic_matter(nn),totalCarbon=layer_total_C)
        layer_max_C=dz(nn)*max_soil_C_density
        layer_extra_C=layer_total_C-layer_max_C
        if (layer_extra_C>0) then
            fraction_to_remove=1.0-layer_max_C/layer_total_C
            call transfer_pool_fraction(soil%soil_organic_matter(nn),soil%soil_organic_matter(nn+1),fraction_to_remove)
            soil%is_peat(nn)=1
        endif

        ! If layer is below limit, carbon is moved back up only if this layer and the layer below are all organic
        if (layer_extra_C < 0 .and. (soil%is_peat(nn).ne.0) .and. (soil%is_peat(nn+1).ne.0)) then
             call poolTotals(soil%soil_organic_matter(nn+1),totalCarbon=layer_total_C_2)
             fraction_to_remove = -layer_extra_C/layer_total_C_2
             if (fraction_to_remove > 0.5) then
                soil%is_peat(nn+1)=0
             else
                call transfer_pool_fraction(soil%soil_organic_matter(nn+1),soil%soil_organic_matter(nn),fraction_to_remove)
             endif
        endif
    enddo

    total_C_after=0.0
    total_N_after=0.0
    do nn=1,num_l
    call poolTotals(soil%soil_organic_matter(num_l),totalCarbon=layer_total_C,totalNitrogen=layer_total_N)
    total_C_after=total_C_after+layer_total_C
    total_N_after=total_N_after+layer_total_N
    enddo

    call poolTotals(soil%leaflitter,totalCarbon=leaflitter_total_C,totalNitrogen=leaflitter_total_N)
    call poolTotals(soil%coarseWoodLitter,totalCarbon=woodlitter_total_C,totalNitrogen=woodlitter_total_N)
    layer_total_C=leaflitter_total_C+woodlitter_total_C
    layer_total_N=leaflitter_total_N+woodlitter_total_N

    total_C_after=total_C_after+layer_total_C
    total_N_after=total_N_after+layer_total_N

    if (abs(total_C_before-total_C_after)>1e-10) then
            print *,'Carbon before:',total_C_before
            print *,'Carbon after:',total_C_after
            call error_mesg('redistribute_peat_carbon','Carbon not conserved after downward move',FATAL)
    endif

    if (abs(total_N_before-total_N_after)>1e-10) then
            print *,'Nitrogen before:',total_N_before
            print *,'Nitrogen after:',total_N_after
            call error_mesg('redistribute_peat_carbon','Nitrogen not conserved after downward move',FATAL)
    endif
end subroutine redistribute_peat_carbon

! ============================================================================
! tile existence detector: returns a logical value indicating whether component
! model tile exists or not
logical function soil_tile_exists(tile)
   type(land_tile_type), pointer :: tile
   soil_tile_exists = associated(tile%soil)
end function soil_tile_exists


!! ZMS Functions moved here from soil_tile_mod to avoid circular dependencies
!! with hillslope_mod and hillslope_tile_mod.
!============================================================================
function soil_cover_cold_start(land_mask, lonb, latb) result(soil_frac)
! creates and initializes a field of fractional soil coverage
  logical, intent(in) :: land_mask(:,:)    ! land mask
  real,    intent(in) :: lonb(:,:), latb(:,:) ! boundaries of the grid cells
  real,    pointer :: soil_frac (:,:,:) ! output: map of soil fractional coverage
  integer :: k

  if (.not. do_hillslope_model) then
     allocate( soil_frac(size(land_mask,1),size(land_mask,2),n_dim_soil_types))
     allocate( soil_tags(size(land_mask,1),size(land_mask,2),n_dim_soil_types))
     do k = 1, n_dim_soil_types
        soil_tags(:,:, k) = k
     end do
  else
     allocate( soil_frac(size(land_mask,1),size(land_mask,2), &
               n_dim_soil_types * num_vertclusters * max_num_topo_hlsps))
     allocate( soil_tags(size(land_mask,1),size(land_mask,2), &
               n_dim_soil_types * num_vertclusters * max_num_topo_hlsps))
     do k = 1,size(soil_tags,3)
        if (mod(k, n_dim_soil_types) == 0) then
           soil_tags(:,:,k) = n_dim_soil_types
        else
           soil_tags(:,:,k) = mod(k, n_dim_soil_types)
        end if
     end do
  end if

  call init_cover_field(soil_to_use, trim(soil_type_file), 'cover','frac', &
       lonb, latb, soil_index_constant, input_cover_types, soil_frac)

  if (do_hillslope_model) then
     call hlsp_coldfracs(soil_frac, n_dim_soil_types)
  end if

end function soil_cover_cold_start


! ============================================================================
subroutine retrieve_soil_tags(soiltags)
! passes back field of soil cover types from soil_cover_cold_start
  integer, pointer :: soiltags(:,:,:)

  if (size(soiltags,1) /= size(soil_tags,1) .or. size(soiltags,2) /= size(soil_tags,2) &
         .or. size(soiltags,3) /= size(soil_tags,3)) &
     call error_mesg(module_name,'Wrong dimension size in "soil_tile_mod:retrieve_soil_tags.',FATAL)

  soiltags(:,:,:) = soil_tags(:,:,:)

  deallocate(soil_tags)
end subroutine retrieve_soil_tags

! ============================================================================
! Calculates soil temperature profile for soil tile, if use_coldstart_wtt_data
subroutine init_soil_twc(soil, ref_soil_t, mwc)
   type(soil_tile_type), intent(inout) :: soil
   real, intent(in)   :: ref_soil_t ! local reference temperature from input data [K]
   real, intent(in), dimension(num_l)  :: mwc  ! total water content to init [kg/m^3]
   integer :: l ! soil level
   real, dimension(num_l) :: vlc, vsc ! volumetric liquid and solid water contents [-]
   real  :: soil_E_max ! not used
   real, dimension(num_l)  :: thermal_cond ! soil thermal conductivity [W/m/K]
   real  :: tres       ! thermal resistance [K / (W/m^2)]

   ! First tentatively initialize soil water / ice content, assuming isothermal profile.
   if (ref_soil_t.ge.soil%pars%tfreeze) then
      soil%wl(1:num_l) = mwc(1:num_l)*dz(1:num_l)
      soil%ws(1:num_l) = 0.
   else
      soil%wl(1:num_l) = 0.
      soil%ws(1:num_l) = mwc(1:num_l)*dz(1:num_l)
   endif
   ! Initialize T and groundwater
   soil%T             = ref_soil_t
   soil%groundwater   = init_groundwater
   soil%groundwater_T = ref_soil_t
   soil%uptake_T           = ref_soil_t
   ! Initialize thermal conductivity
   vlc(:) = soil%wl(1:num_l)/ (dz(1:num_l)*dens_h2o*soil%pars%vwc_sat)
   vsc(:) = soil%ws(1:num_l)/ (dz(1:num_l)*dens_h2o*soil%pars%vwc_sat)
   call soil_data_thermodynamics ( soil, vlc, vsc, soil_E_max, thermal_cond)

   ! Walk down through soil and maintain geothermal heat flux profile.
   do l=2,num_l
      tres = 0.5*(dz(l-1)/thermal_cond(l-1) + dz(l)/thermal_cond(l)) ! [K / (W/m^2)] =  [m / (W/m/K)]
      soil%T(l) = soil%T(l-1) + soil%geothermal_heat_flux * tres

      ! Check for switch to above-freezing
      if (soil%T(l) >= soil%pars%tfreeze .and. soil%T(l-1) < soil%pars%tfreeze) then
         ! Reset wl & ws, and recalculate thermal conductivity, from here down.
         soil%wl(l:num_l) = mwc(l:num_l)*dz(l:num_l)
         soil%ws(l:num_l) = 0.
         vlc(l:num_l) = soil%wl(l:num_l)/ (dz(l:num_l)*dens_h2o*soil%pars%vwc_sat)
         vsc(l:num_l) = soil%ws(l:num_l)/ (dz(l:num_l)*dens_h2o*soil%pars%vwc_sat)
         call soil_data_thermodynamics ( soil, vlc, vsc, soil_E_max, thermal_cond)
      end if
   end do

   ! Debug
   ! current point set above call in soil_init
   if (is_watch_point()) then
      write(*,*)'Initializing soil temperature and water for watch_point.'
      __DEBUG1__(ref_soil_t)
      do l=1,num_l
         __DEBUG4__(l, vlc(l), vsc(l), soil%T(l))
      end do
   end if


end subroutine init_soil_twc


! ============================================================================
! cohort accessor functions: given a pointer to cohort, return a pointer to a
! specific member of the cohort structure
#define DEFINE_SOIL_ACCESSOR_0D(xtype,x) subroutine soil_ ## x ## _ptr(t,p);\
type(land_tile_type),pointer::t;xtype,pointer::p;p=>NULL();if(associated(t))then;if(associated(t%soil))p=>t%soil%x;endif;end subroutine
#define DEFINE_SOIL_ACCESSOR_1D(xtype,x) subroutine soil_ ## x ## _ptr(t,p);\
type(land_tile_type),pointer::t;xtype,pointer::p(:);p=>NULL();if(associated(t))then;if(associated(t%soil))p=>t%soil%x;endif;end subroutine
#define DEFINE_SOIL_COMPONENT_ACCESSOR_0D(xtype,component,x) subroutine soil_ ## x ## _ptr(t,p);\
type(land_tile_type),pointer::t;xtype,pointer::p;p=>NULL();if(associated(t))then;if(associated(t%soil))p=>t%soil%component%x;endif;end subroutine
#define DEFINE_SOIL_COMPONENT_ACCESSOR_1D(xtype,component,x) subroutine soil_ ## x ## _ptr(t,p);\
type(land_tile_type),pointer::t;xtype,pointer::p(:);p=>NULL();if(associated(t))then;if(associated(t%soil))p=>t%soil%component%x;endif;end subroutine

DEFINE_SOIL_ACCESSOR_1D(real,w_fc)
DEFINE_SOIL_ACCESSOR_1D(real,alpha)
DEFINE_SOIL_ACCESSOR_0D(real,uptake_T)
DEFINE_SOIL_ACCESSOR_0D(integer,tag)
DEFINE_SOIL_ACCESSOR_1D(real,fast_soil_C)
DEFINE_SOIL_ACCESSOR_1D(real,slow_soil_C)
DEFINE_SOIL_ACCESSOR_1D(real,T)
DEFINE_SOIL_ACCESSOR_1D(real,wl)
DEFINE_SOIL_ACCESSOR_1D(real,ws)
DEFINE_SOIL_ACCESSOR_1D(real,groundwater)
DEFINE_SOIL_ACCESSOR_1D(real,groundwater_T)
DEFINE_SOIL_ACCESSOR_1D(integer,is_peat)

DEFINE_SOIL_COMPONENT_ACCESSOR_0D(real,pars,tau_groundwater)
DEFINE_SOIL_COMPONENT_ACCESSOR_0D(real,pars,hillslope_length)
DEFINE_SOIL_COMPONENT_ACCESSOR_0D(real,pars,hillslope_relief)
DEFINE_SOIL_COMPONENT_ACCESSOR_0D(real,pars,hillslope_a)
DEFINE_SOIL_COMPONENT_ACCESSOR_0D(real,pars,hillslope_n)
DEFINE_SOIL_COMPONENT_ACCESSOR_0D(real,pars,hillslope_zeta_bar)
DEFINE_SOIL_COMPONENT_ACCESSOR_0D(real,pars,soil_e_depth)
DEFINE_SOIL_COMPONENT_ACCESSOR_0D(real,pars,zeta)
DEFINE_SOIL_COMPONENT_ACCESSOR_0D(real,pars,tau)
DEFINE_SOIL_COMPONENT_ACCESSOR_0D(real,pars,k_sat_gw)
DEFINE_SOIL_COMPONENT_ACCESSOR_0D(real,pars,vwc_wilt)
DEFINE_SOIL_COMPONENT_ACCESSOR_0D(real,pars,vwc_fc)
DEFINE_SOIL_COMPONENT_ACCESSOR_0D(real,pars,vwc_sat)
DEFINE_SOIL_COMPONENT_ACCESSOR_0D(real,pars,k_sat_ref)
DEFINE_SOIL_COMPONENT_ACCESSOR_0D(real,pars,Qmax)

DEFINE_SOIL_COMPONENT_ACCESSOR_1D(real,pars,refl_dry_dir)
DEFINE_SOIL_COMPONENT_ACCESSOR_1D(real,pars,refl_dry_dif)
DEFINE_SOIL_COMPONENT_ACCESSOR_1D(real,pars,refl_sat_dir)
DEFINE_SOIL_COMPONENT_ACCESSOR_1D(real,pars,refl_sat_dif)
DEFINE_SOIL_COMPONENT_ACCESSOR_1D(real,pars,f_iso_dry)
DEFINE_SOIL_COMPONENT_ACCESSOR_1D(real,pars,f_vol_dry)
DEFINE_SOIL_COMPONENT_ACCESSOR_1D(real,pars,f_geo_dry)
DEFINE_SOIL_COMPONENT_ACCESSOR_1D(real,pars,f_iso_sat)
DEFINE_SOIL_COMPONENT_ACCESSOR_1D(real,pars,f_vol_sat)
DEFINE_SOIL_COMPONENT_ACCESSOR_1D(real,pars,f_geo_sat)

DEFINE_SOIL_ACCESSOR_0D(real,fast_DOC_leached)
DEFINE_SOIL_ACCESSOR_0D(real,slow_DOC_leached)
DEFINE_SOIL_ACCESSOR_0D(real,deadmic_DOC_leached)
DEFINE_SOIL_ACCESSOR_0D(real,fast_DON_leached)
DEFINE_SOIL_ACCESSOR_0D(real,slow_DON_leached)
DEFINE_SOIL_ACCESSOR_0D(real,deadmic_DON_leached)
DEFINE_SOIL_ACCESSOR_0D(real,NO3_leached)
DEFINE_SOIL_ACCESSOR_0D(real,NH4_leached)

DEFINE_SOIL_ACCESSOR_1D(real,asoil_in)
DEFINE_SOIL_ACCESSOR_1D(real,fsc_in)
DEFINE_SOIL_ACCESSOR_1D(real,ssc_in)
DEFINE_SOIL_ACCESSOR_1D(real,deadmic_C_in)
DEFINE_SOIL_ACCESSOR_1D(real,fsn_in)
DEFINE_SOIL_ACCESSOR_1D(real,ssn_in)
DEFINE_SOIL_ACCESSOR_1D(real,deadmic_N_in)
DEFINE_SOIL_ACCESSOR_1D(real,fast_protected_C_in)
DEFINE_SOIL_ACCESSOR_1D(real,slow_protected_C_in)
DEFINE_SOIL_ACCESSOR_1D(real,deadmic_protected_C_in)
DEFINE_SOIL_ACCESSOR_1D(real,fast_protected_N_in)
DEFINE_SOIL_ACCESSOR_1D(real,slow_protected_N_in)
DEFINE_SOIL_ACCESSOR_1D(real,deadmic_protected_N_in)

DEFINE_SOIL_ACCESSOR_0D(real,leaflitter_fsc_in)
DEFINE_SOIL_ACCESSOR_0D(real,leaflitter_ssc_in)
DEFINE_SOIL_ACCESSOR_0D(real,leaflitter_deadmic_C_in)
DEFINE_SOIL_ACCESSOR_0D(real,leaflitter_fsn_in)
DEFINE_SOIL_ACCESSOR_0D(real,leaflitter_ssn_in)
DEFINE_SOIL_ACCESSOR_0D(real,leaflitter_deadmic_N_in)

DEFINE_SOIL_ACCESSOR_0D(real,finewoodlitter_fsc_in)
DEFINE_SOIL_ACCESSOR_0D(real,finewoodlitter_ssc_in)
DEFINE_SOIL_ACCESSOR_0D(real,finewoodlitter_deadmic_C_in)
DEFINE_SOIL_ACCESSOR_0D(real,finewoodlitter_fsn_in)
DEFINE_SOIL_ACCESSOR_0D(real,finewoodlitter_ssn_in)
DEFINE_SOIL_ACCESSOR_0D(real,finewoodlitter_deadmic_N_in)

DEFINE_SOIL_ACCESSOR_0D(real,coarsewoodlitter_fsc_in)
DEFINE_SOIL_ACCESSOR_0D(real,coarsewoodlitter_ssc_in)
DEFINE_SOIL_ACCESSOR_0D(real,coarsewoodlitter_deadmic_C_in)
DEFINE_SOIL_ACCESSOR_0D(real,coarsewoodlitter_fsn_in)
DEFINE_SOIL_ACCESSOR_0D(real,coarsewoodlitter_ssn_in)
DEFINE_SOIL_ACCESSOR_0D(real,coarsewoodlitter_deadmic_N_in)

DEFINE_SOIL_ACCESSOR_1D(real,fast_C_turnover_accumulated)
DEFINE_SOIL_ACCESSOR_1D(real,slow_C_turnover_accumulated)
DEFINE_SOIL_ACCESSOR_1D(real,deadmic_C_turnover_accumulated)
DEFINE_SOIL_ACCESSOR_1D(real,fast_protected_C_turnover_accumulated)
DEFINE_SOIL_ACCESSOR_1D(real,slow_protected_C_turnover_accumulated)
DEFINE_SOIL_ACCESSOR_1D(real,deadmic_protected_C_turnover_accumulated)
DEFINE_SOIL_ACCESSOR_0D(real,leaflitter_fast_C_turnover_accumulated)
DEFINE_SOIL_ACCESSOR_0D(real,leaflitter_slow_C_turnover_accumulated)
DEFINE_SOIL_ACCESSOR_0D(real,leaflitter_deadmic_C_turnover_accumulated)
DEFINE_SOIL_ACCESSOR_0D(real,finewoodlitter_fast_C_turnover_accumulated)
DEFINE_SOIL_ACCESSOR_0D(real,finewoodlitter_slow_C_turnover_accumulated)
DEFINE_SOIL_ACCESSOR_0D(real,finewoodlitter_deadmic_C_turnover_accumulated)
DEFINE_SOIL_ACCESSOR_0D(real,coarsewoodlitter_fast_C_turnover_accumulated)
DEFINE_SOIL_ACCESSOR_0D(real,coarsewoodlitter_slow_C_turnover_accumulated)
DEFINE_SOIL_ACCESSOR_0D(real,coarsewoodlitter_deadmic_C_turnover_accumulated)

DEFINE_SOIL_ACCESSOR_1D(real,fast_N_turnover_accumulated)
DEFINE_SOIL_ACCESSOR_1D(real,slow_N_turnover_accumulated)
DEFINE_SOIL_ACCESSOR_1D(real,deadmic_N_turnover_accumulated)
DEFINE_SOIL_ACCESSOR_1D(real,fast_protected_N_turnover_accumulated)
DEFINE_SOIL_ACCESSOR_1D(real,slow_protected_N_turnover_accumulated)
DEFINE_SOIL_ACCESSOR_1D(real,deadmic_protected_N_turnover_accumulated)
DEFINE_SOIL_ACCESSOR_0D(real,leaflitter_fast_N_turnover_accumulated)
DEFINE_SOIL_ACCESSOR_0D(real,leaflitter_slow_N_turnover_accumulated)
DEFINE_SOIL_ACCESSOR_0D(real,leaflitter_deadmic_N_turnover_accumulated)
DEFINE_SOIL_ACCESSOR_0D(real,finewoodlitter_fast_N_turnover_accumulated)
DEFINE_SOIL_ACCESSOR_0D(real,finewoodlitter_slow_N_turnover_accumulated)
DEFINE_SOIL_ACCESSOR_0D(real,finewoodlitter_deadmic_N_turnover_accumulated)
DEFINE_SOIL_ACCESSOR_0D(real,coarsewoodlitter_fast_N_turnover_accumulated)
DEFINE_SOIL_ACCESSOR_0D(real,coarsewoodlitter_slow_N_turnover_accumulated)
DEFINE_SOIL_ACCESSOR_0D(real,coarsewoodlitter_deadmic_N_turnover_accumulated)

! stuff below is for CORPSE
#define DEFINE_SOIL_LAYER_COHORT_COMPONENT_ACCESSOR(xtype,x) subroutine soilc_ ## x ## _ptr(t,p,layer);\
type(land_tile_type),pointer::t;xtype,pointer::p(:);integer,intent(in)::layer;p=>NULL();if(associated(t))then;if(associated(t%soil))call get_pool_data_accessors(t%soil%soil_organic_matter(layer),x=p);endif;end subroutine
#define DEFINE_SOIL_C_POOL_COMPONENT_ACCESSOR(xtype,pool,x) subroutine soilc_ ## pool ## _ ## x ## _ptr(t,p);\
type(land_tile_type),pointer::t;xtype,pointer::p(:);p=>NULL();if(associated(t))then;if(associated(t%soil))call get_pool_data_accessors(t%soil%pool,x=p);endif;end subroutine
#define DEFINE_SOIL_C_POOL_NONCOHORT_COMPONENT_ACCESSOR(xtype,pool,x) subroutine soilc_ ## pool ## _ ## x ## _ptr(t,p);\
type(land_tile_type),pointer::t;xtype,pointer::p;p=>NULL();if(associated(t))then;if(associated(t%soil))call get_pool_data_accessors(t%soil%pool,x=p);endif;end subroutine

DEFINE_SOIL_LAYER_COHORT_COMPONENT_ACCESSOR(real,fast_soil_C)
DEFINE_SOIL_LAYER_COHORT_COMPONENT_ACCESSOR(real,slow_soil_C)
DEFINE_SOIL_LAYER_COHORT_COMPONENT_ACCESSOR(real,deadMicrobeC)
DEFINE_SOIL_LAYER_COHORT_COMPONENT_ACCESSOR(real,fast_protected_C)
DEFINE_SOIL_LAYER_COHORT_COMPONENT_ACCESSOR(real,slow_protected_C)
DEFINE_SOIL_LAYER_COHORT_COMPONENT_ACCESSOR(real,deadMicrobe_protected_C)
DEFINE_SOIL_LAYER_COHORT_COMPONENT_ACCESSOR(real,livingMicrobeC)
DEFINE_SOIL_LAYER_COHORT_COMPONENT_ACCESSOR(real,Rtot)
DEFINE_SOIL_LAYER_COHORT_COMPONENT_ACCESSOR(real,CO2)
DEFINE_SOIL_LAYER_COHORT_COMPONENT_ACCESSOR(real,originalLitterC)

DEFINE_SOIL_LAYER_COHORT_COMPONENT_ACCESSOR(real,fast_soil_N)
DEFINE_SOIL_LAYER_COHORT_COMPONENT_ACCESSOR(real,slow_soil_N)
DEFINE_SOIL_LAYER_COHORT_COMPONENT_ACCESSOR(real,deadMicrobeN)
DEFINE_SOIL_LAYER_COHORT_COMPONENT_ACCESSOR(real,fast_protected_N)
DEFINE_SOIL_LAYER_COHORT_COMPONENT_ACCESSOR(real,slow_protected_N)
DEFINE_SOIL_LAYER_COHORT_COMPONENT_ACCESSOR(real,deadMicrobe_protected_N)
DEFINE_SOIL_LAYER_COHORT_COMPONENT_ACCESSOR(real,livingMicrobeN)
DEFINE_SOIL_LAYER_COHORT_COMPONENT_ACCESSOR(real,originalLitterN)

DEFINE_SOIL_C_POOL_COMPONENT_ACCESSOR(real,leafLitter,fast_soil_C)
DEFINE_SOIL_C_POOL_COMPONENT_ACCESSOR(real,leafLitter,slow_soil_C)
DEFINE_SOIL_C_POOL_COMPONENT_ACCESSOR(real,leafLitter,deadMicrobeC)
DEFINE_SOIL_C_POOL_COMPONENT_ACCESSOR(real,leafLitter,livingMicrobeC)
DEFINE_SOIL_C_POOL_COMPONENT_ACCESSOR(real,leafLitter,Rtot)
DEFINE_SOIL_C_POOL_COMPONENT_ACCESSOR(real,leafLitter,CO2)
DEFINE_SOIL_C_POOL_COMPONENT_ACCESSOR(real,leafLitter,originalLitterC)
DEFINE_SOIL_C_POOL_COMPONENT_ACCESSOR(real,leafLitter,fast_protected_C)
DEFINE_SOIL_C_POOL_COMPONENT_ACCESSOR(real,leafLitter,slow_protected_C)
DEFINE_SOIL_C_POOL_COMPONENT_ACCESSOR(real,leafLitter,deadMicrobe_protected_C)

DEFINE_SOIL_C_POOL_COMPONENT_ACCESSOR(real,leafLitter,fast_soil_N)
DEFINE_SOIL_C_POOL_COMPONENT_ACCESSOR(real,leafLitter,slow_soil_N)
DEFINE_SOIL_C_POOL_COMPONENT_ACCESSOR(real,leafLitter,deadMicrobeN)
DEFINE_SOIL_C_POOL_COMPONENT_ACCESSOR(real,leafLitter,livingMicrobeN)
DEFINE_SOIL_C_POOL_COMPONENT_ACCESSOR(real,leafLitter,originalLitterN)
DEFINE_SOIL_C_POOL_COMPONENT_ACCESSOR(real,leafLitter,fast_protected_N)
DEFINE_SOIL_C_POOL_COMPONENT_ACCESSOR(real,leafLitter,slow_protected_N)
DEFINE_SOIL_C_POOL_COMPONENT_ACCESSOR(real,leafLitter,deadMicrobe_protected_N)

DEFINE_SOIL_C_POOL_NONCOHORT_COMPONENT_ACCESSOR(real,leafLitter,fast_DOC)
DEFINE_SOIL_C_POOL_NONCOHORT_COMPONENT_ACCESSOR(real,leafLitter,slow_DOC)
DEFINE_SOIL_C_POOL_NONCOHORT_COMPONENT_ACCESSOR(real,leafLitter,deadMicrobe_DOC)
DEFINE_SOIL_C_POOL_NONCOHORT_COMPONENT_ACCESSOR(real,leafLitter,fast_DON)
DEFINE_SOIL_C_POOL_NONCOHORT_COMPONENT_ACCESSOR(real,leafLitter,slow_DON)
DEFINE_SOIL_C_POOL_NONCOHORT_COMPONENT_ACCESSOR(real,leafLitter,deadMicrobe_DON)
DEFINE_SOIL_C_POOL_NONCOHORT_COMPONENT_ACCESSOR(real,leafLitter,nitrate)
DEFINE_SOIL_C_POOL_NONCOHORT_COMPONENT_ACCESSOR(real,leafLitter,ammonium)
DEFINE_SOIL_C_POOL_NONCOHORT_COMPONENT_ACCESSOR(real,leafLitter,nitrif)
DEFINE_SOIL_C_POOL_NONCOHORT_COMPONENT_ACCESSOR(real,leafLitter,denitrif)


DEFINE_SOIL_C_POOL_COMPONENT_ACCESSOR(real,fineWoodLitter,fast_soil_C)
DEFINE_SOIL_C_POOL_COMPONENT_ACCESSOR(real,fineWoodLitter,slow_soil_C)
DEFINE_SOIL_C_POOL_COMPONENT_ACCESSOR(real,fineWoodLitter,deadMicrobeC)
DEFINE_SOIL_C_POOL_COMPONENT_ACCESSOR(real,fineWoodLitter,livingMicrobeC)
DEFINE_SOIL_C_POOL_COMPONENT_ACCESSOR(real,fineWoodLitter,Rtot)
DEFINE_SOIL_C_POOL_COMPONENT_ACCESSOR(real,fineWoodLitter,CO2)
DEFINE_SOIL_C_POOL_COMPONENT_ACCESSOR(real,fineWoodLitter,originalLitterC)
DEFINE_SOIL_C_POOL_COMPONENT_ACCESSOR(real,fineWoodLitter,fast_protected_C)
DEFINE_SOIL_C_POOL_COMPONENT_ACCESSOR(real,fineWoodLitter,slow_protected_C)
DEFINE_SOIL_C_POOL_COMPONENT_ACCESSOR(real,fineWoodLitter,deadMicrobe_protected_C)


DEFINE_SOIL_C_POOL_COMPONENT_ACCESSOR(real,fineWoodLitter,fast_soil_N)
DEFINE_SOIL_C_POOL_COMPONENT_ACCESSOR(real,fineWoodLitter,slow_soil_N)
DEFINE_SOIL_C_POOL_COMPONENT_ACCESSOR(real,fineWoodLitter,deadMicrobeN)
DEFINE_SOIL_C_POOL_COMPONENT_ACCESSOR(real,fineWoodLitter,livingMicrobeN)
DEFINE_SOIL_C_POOL_COMPONENT_ACCESSOR(real,fineWoodLitter,originalLitterN)
DEFINE_SOIL_C_POOL_COMPONENT_ACCESSOR(real,fineWoodLitter,fast_protected_N)
DEFINE_SOIL_C_POOL_COMPONENT_ACCESSOR(real,fineWoodLitter,slow_protected_N)
DEFINE_SOIL_C_POOL_COMPONENT_ACCESSOR(real,fineWoodLitter,deadMicrobe_protected_N)

DEFINE_SOIL_C_POOL_NONCOHORT_COMPONENT_ACCESSOR(real,fineWoodLitter,fast_DOC)
DEFINE_SOIL_C_POOL_NONCOHORT_COMPONENT_ACCESSOR(real,fineWoodLitter,slow_DOC)
DEFINE_SOIL_C_POOL_NONCOHORT_COMPONENT_ACCESSOR(real,fineWoodLitter,deadMicrobe_DOC)
DEFINE_SOIL_C_POOL_NONCOHORT_COMPONENT_ACCESSOR(real,fineWoodLitter,fast_DON)
DEFINE_SOIL_C_POOL_NONCOHORT_COMPONENT_ACCESSOR(real,fineWoodLitter,slow_DON)
DEFINE_SOIL_C_POOL_NONCOHORT_COMPONENT_ACCESSOR(real,fineWoodLitter,deadMicrobe_DON)
DEFINE_SOIL_C_POOL_NONCOHORT_COMPONENT_ACCESSOR(real,fineWoodLitter,nitrate)
DEFINE_SOIL_C_POOL_NONCOHORT_COMPONENT_ACCESSOR(real,fineWoodLitter,ammonium)
DEFINE_SOIL_C_POOL_NONCOHORT_COMPONENT_ACCESSOR(real,fineWoodLitter,nitrif)
DEFINE_SOIL_C_POOL_NONCOHORT_COMPONENT_ACCESSOR(real,fineWoodLitter,denitrif)

DEFINE_SOIL_C_POOL_COMPONENT_ACCESSOR(real,coarseWoodLitter,fast_soil_C)
DEFINE_SOIL_C_POOL_COMPONENT_ACCESSOR(real,coarseWoodLitter,slow_soil_C)
DEFINE_SOIL_C_POOL_COMPONENT_ACCESSOR(real,coarseWoodLitter,deadMicrobeC)
DEFINE_SOIL_C_POOL_COMPONENT_ACCESSOR(real,coarseWoodLitter,livingMicrobeC)
DEFINE_SOIL_C_POOL_COMPONENT_ACCESSOR(real,coarseWoodLitter,Rtot)
DEFINE_SOIL_C_POOL_COMPONENT_ACCESSOR(real,coarseWoodLitter,CO2)
DEFINE_SOIL_C_POOL_COMPONENT_ACCESSOR(real,coarseWoodLitter,originalLitterC)
DEFINE_SOIL_C_POOL_COMPONENT_ACCESSOR(real,coarseWoodLitter,fast_protected_C)
DEFINE_SOIL_C_POOL_COMPONENT_ACCESSOR(real,coarseWoodLitter,slow_protected_C)
DEFINE_SOIL_C_POOL_COMPONENT_ACCESSOR(real,coarseWoodLitter,deadMicrobe_protected_C)


DEFINE_SOIL_C_POOL_COMPONENT_ACCESSOR(real,coarseWoodLitter,fast_soil_N)
DEFINE_SOIL_C_POOL_COMPONENT_ACCESSOR(real,coarseWoodLitter,slow_soil_N)
DEFINE_SOIL_C_POOL_COMPONENT_ACCESSOR(real,coarseWoodLitter,deadMicrobeN)
DEFINE_SOIL_C_POOL_COMPONENT_ACCESSOR(real,coarseWoodLitter,livingMicrobeN)
DEFINE_SOIL_C_POOL_COMPONENT_ACCESSOR(real,coarseWoodLitter,originalLitterN)
DEFINE_SOIL_C_POOL_COMPONENT_ACCESSOR(real,coarseWoodLitter,fast_protected_N)
DEFINE_SOIL_C_POOL_COMPONENT_ACCESSOR(real,coarseWoodLitter,slow_protected_N)
DEFINE_SOIL_C_POOL_COMPONENT_ACCESSOR(real,coarseWoodLitter,deadMicrobe_protected_N)


DEFINE_SOIL_C_POOL_NONCOHORT_COMPONENT_ACCESSOR(real,coarseWoodLitter,fast_DOC)
DEFINE_SOIL_C_POOL_NONCOHORT_COMPONENT_ACCESSOR(real,coarseWoodLitter,slow_DOC)
DEFINE_SOIL_C_POOL_NONCOHORT_COMPONENT_ACCESSOR(real,coarseWoodLitter,deadMicrobe_DOC)
DEFINE_SOIL_C_POOL_NONCOHORT_COMPONENT_ACCESSOR(real,coarseWoodLitter,fast_DON)
DEFINE_SOIL_C_POOL_NONCOHORT_COMPONENT_ACCESSOR(real,coarseWoodLitter,slow_DON)
DEFINE_SOIL_C_POOL_NONCOHORT_COMPONENT_ACCESSOR(real,coarseWoodLitter,deadMicrobe_DON)
DEFINE_SOIL_C_POOL_NONCOHORT_COMPONENT_ACCESSOR(real,coarseWoodLitter,nitrate)
DEFINE_SOIL_C_POOL_NONCOHORT_COMPONENT_ACCESSOR(real,coarseWoodLitter,ammonium)
DEFINE_SOIL_C_POOL_NONCOHORT_COMPONENT_ACCESSOR(real,coarseWoodLitter,nitrif)
DEFINE_SOIL_C_POOL_NONCOHORT_COMPONENT_ACCESSOR(real,coarseWoodLitter,denitrif)

subroutine soil_fast_DOC_ptr(t,p)
type(land_tile_type),pointer::t
real,pointer::p(:)
p=>NULL()
if(associated(t))then
if(associated(t%soil))p=>t%soil%soil_organic_matter(:)%dissolved_carbon(1)
endif
end subroutine

subroutine soil_slow_DOC_ptr(t,p)
type(land_tile_type),pointer::t
real,pointer::p(:)
p=>NULL()
if(associated(t))then
if(associated(t%soil))p=>t%soil%soil_organic_matter(:)%dissolved_carbon(2)
endif
end subroutine

subroutine soil_deadMicrobe_DOC_ptr(t,p)
type(land_tile_type),pointer::t
real,pointer::p(:)
p=>NULL()
if(associated(t))then
if(associated(t%soil))p=>t%soil%soil_organic_matter(:)%dissolved_carbon(3)
endif
end subroutine

subroutine soil_fast_DON_ptr(t,p)
type(land_tile_type),pointer::t
real,pointer::p(:)
p=>NULL()
if(associated(t))then
if(associated(t%soil))p=>t%soil%soil_organic_matter(:)%dissolved_nitrogen(1)
endif
end subroutine

subroutine soil_slow_DON_ptr(t,p)
type(land_tile_type),pointer::t
real,pointer::p(:)
p=>NULL()
if(associated(t))then
if(associated(t%soil))p=>t%soil%soil_organic_matter(:)%dissolved_nitrogen(2)
endif
end subroutine

subroutine soil_deadMicrobe_DON_ptr(t,p)
type(land_tile_type),pointer::t
real,pointer::p(:)
p=>NULL()
if(associated(t))then
if(associated(t%soil))p=>t%soil%soil_organic_matter(:)%dissolved_nitrogen(3)
endif
end subroutine

subroutine soil_nitrate_ptr(t,p)
type(land_tile_type),pointer::t
real,pointer::p(:)
p=>NULL()
if(associated(t))then
if(associated(t%soil))p=>t%soil%soil_organic_matter(:)%nitrate
endif
end subroutine

subroutine soil_ammonium_ptr(t,p)
type(land_tile_type),pointer::t
real,pointer::p(:)
p=>NULL()
if(associated(t))then
if(associated(t%soil))p=>t%soil%soil_organic_matter(:)%ammonium
endif
end subroutine

subroutine soil_nitrif_ptr(t,p)
type(land_tile_type),pointer::t
real,pointer::p(:)
p=>NULL()
if(associated(t))then
if(associated(t%soil))p=>t%soil%soil_organic_matter(:)%nitrif
endif
end subroutine

subroutine soil_denitrif_ptr(t,p)
type(land_tile_type),pointer::t
real,pointer::p(:)
p=>NULL()
if(associated(t))then
if(associated(t%soil))p=>t%soil%soil_organic_matter(:)%denitrif
endif
end subroutine

end module soil_mod
