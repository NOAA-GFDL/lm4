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

use fms_mod, only: error_mesg, string, file_exist, check_nml_error, &
     stdlog, close_file, mpp_pe, mpp_root_pe, FATAL, WARNING, NOTE
use time_manager_mod,   only: time_type, time_type_to_real
use diag_manager_mod,   only: diag_axis_init
use constants_mod,      only: pi, tfreeze, hlv, hlf, dens_h2o
use tracer_manager_mod, only: NO_TRACER

use land_constants_mod, only : NBANDS, BAND_VIS, BAND_NIR, seconds_per_year
use land_numerics_mod, only : tridiag
use soil_tile_mod, only : num_l, dz, zfull, zhalf, &
     GW_LM2, GW_LINEAR, GW_HILL_AR5, GW_HILL, GW_TILED, &
     soil_tile_type, read_soil_data_namelist, &
     soil_data_radiation, soil_data_thermodynamics, &
     soil_data_hydraulic_properties, soil_data_psi_for_rh, &
     soil_data_gw_hydraulics, soil_data_gw_hydraulics_ar5, &
     soil_data_vwc_for_init_only, &
     soil_data_init_derive_subsurf_pars, &
     soil_data_init_derive_subsurf_pars_ar5, &
     soil_data_init_derive_subsurf_pars_tiled, &
     psi_wilt, cpw, clw, csw, g_iso, g_vol, g_geo, g_RT, aspect,&
     gw_scale_length, gw_scale_relief, gw_scale_soil_depth, &
     slope_exp, gw_scale_perm, k0_macro_x, retro_a0n1, &
     soil_type_file, &
     soil_tile_stock_pe, initval, comp, soil_theta, soil_ice_porosity, &
     N_LITTER_POOLS,LEAF,CWOOD,l_shortname,l_longname
use soil_util_mod, only: soil_util_init, rhizosphere_frac

use soil_carbon_mod, only: poolTotals, poolTotals1, soilMaxCohorts, litterDensity,&
     update_pool, tracer_leaching_with_litter,transfer_pool_fraction, N_C_TYPES, &
     soil_carbon_option, SOILC_CENTURY, SOILC_CENTURY_BY_LAYER, SOILC_CORPSE, SOILC_CORPSE_N, &
     C_FAST, C_SLOW, C_MIC, A_function, debug_pool, adjust_pool_ncohorts, c_shortname, c_longname, &
     mycorrhizal_mineral_N_uptake_rate, mycorrhizal_decomposition, ammonium_solubility, nitrate_solubility

use land_tile_mod, only : land_tile_map, land_tile_type, land_tile_enum_type, &
     first_elmt, prev_elmt, loop_over_tiles
use land_utils_mod, only : put_to_tiles_r0d_fptr, put_to_tiles_r1d_fptr
use land_tile_diag_mod, only : diag_buff_type, set_default_diag_filter, &
     register_tiled_static_field, register_tiled_diag_field, &
     send_tile_data, send_tile_data_r0d_fptr, send_tile_data_r1d_fptr, &
     send_tile_data_i0d_fptr, cmor_name, cmor_mrsos_depth, &
     add_tiled_diag_field_alias, add_tiled_static_field_alias, &
     send_cohort_data, register_cohort_diag_field, OP_SUM
use land_data_mod, only : lnd, log_version
use land_io_mod, only : read_field
use land_tile_io_mod, only: land_restart_type, &
     init_land_restart, open_land_restart, save_land_restart, free_land_restart, &
     add_tile_data, add_int_tile_data, get_tile_data, get_int_tile_data, &
     add_restart_axis, field_exists
use vegn_data_mod, only: K1, K2, spdata
use vegn_cohort_mod, only : vegn_cohort_type, &
     cohort_uptake_profile, cohort_root_litter_profile

use vegn_tile_mod, only : vegn_tile_type, vegn_tile_bwood
use land_debug_mod, only : is_watch_point, is_watch_cell, get_current_point, &
     set_current_point, check_var_range, check_conservation, land_error_message
use uptake_mod, only : UPTAKE_LINEAR, UPTAKE_DARCY2D, UPTAKE_DARCY2D_LIN, &
     uptake_init, uptake_option, darcy2d_uptake, darcy2d_uptake_solver

use hillslope_mod, only : do_hillslope_model, max_num_topo_hlsps, &
     num_vertclusters, hlsp_coldfracs, use_geohydrodata, & !pond, &
     horiz_wt_depth_to_init, calculate_wt_init, simple_inundation
use land_io_mod, only : &
     init_cover_field
use soil_tile_mod, only : n_dim_soil_types, soil_to_use, &
     soil_index_constant, input_cover_types
use hillslope_hydrology_mod, only: hlsp_hydro_lev_init, hlsp_hydrology_2, &
     stiff_explicit_gwupdate
use river_mod, only : river_tracer_index

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

public :: soil_get_sfc_temp
public :: soil_radiation
public :: soil_step_1
public :: soil_step_2
public :: soil_step_3
public :: soil_data_beta

public :: Dsdt
public :: get_soil_litter_C
public :: root_N_uptake
public :: myc_scavenger_N_uptake
public :: myc_miner_N_uptake
public :: redistribute_peat_carbon

! helper functions that may be better moved elsewhere:
public :: register_litter_soilc_diag_fields

! ==== module constants ======================================================
character(len=*), parameter :: module_name = 'soil'
#include "../shared/version_variable.inc"

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
logical :: allow_negative_rie   = .false.
logical :: baseflow_where_frozen = .false.
logical :: write_when_flagged   = .false.
logical :: corrected_lm2_gw     = .true.
logical :: use_fringe           = .false.
logical :: push_down_sfc_excess = .true.
logical :: lrunf_from_div       = .true.
logical :: cold_infilt          = .false.
logical :: bottom_up_cold_infilt= .false.
logical :: use_depth_to_wt_4    = .false.
logical :: always_use_bsw       = .false.
logical :: harmonic_mean_K      = .false.
logical :: verbose              = .false.
logical :: no_min_Dpsi          = .false.
logical :: div_bug_fix          = .false.
logical :: require_pos_hcap     = .false.
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
logical :: prohibit_negative_water_div = .false. ! if TRUE, div_bf abd dif_if are set to zero
                                            ! in case water content of *any* layer is negative
real    :: zeta_bar_override    = -1.
real    :: cold_depth           = 0.
real    :: Wl_min               = -1.e20
real    :: bwood_macinf         = -1.
integer :: max_iter_trans = 100 ! max number of iterations for psi_crown_min
integer :: layer_for_gw_switch = 1000000 ! to accelerate permafrost gw shutoff
real    :: eps_trans     = 1.e-7 ! convergence crit for psi_crown_min
logical :: supercooled_rnu = .true. ! Excess ice converted to supercooled water for runoff.
real    :: wet_depth = 0.6 ! [m] water table depth threshold for diagnosing wetland fraction
real    :: thetathresh = -0.01 ! [-] threshold for negative soil liquid water saturation
  ! before warning or abort
real    :: negrnuthresh = -0.1 ! [mm/s] threshold for negative lrunf_nu
  ! before warning or abort

real :: max_soil_C_density = 50.0   !(kgC/m3) -- for redistribution of peat
real :: max_litter_thickness = 0.05 ! m of litter layer thickness before it gets redistributed

real :: r_rhiz = 0.001              ! Radius of rhizosphere around root (m)

namelist /soil_nml/ lm2, use_E_min, use_E_max,           &
                    init_temp,      &
                    init_w,   init_wtdep,    &
                    init_groundwater, lrunf_ie_min, lrunf_ie_tol, &
                    cpw, clw, csw, &
                    albedo_to_use, &
                    allow_negative_rie, &
                    baseflow_where_frozen, &
                    write_when_flagged, &
                    corrected_lm2_gw, &
                    use_fringe, &
                    cold_infilt, bottom_up_cold_infilt, &
                    use_depth_to_wt_4, always_use_bsw, &
                    harmonic_mean_K, verbose, no_min_Dpsi, div_bug_fix, &
                    require_pos_hcap, &
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
real            :: Eg_min

integer         :: gw_option = -1

integer :: n_river_tracers = 0
integer :: i_river_DOC     = NO_TRACER
integer :: i_river_DON     = NO_TRACER
integer :: i_river_NO3     = NO_TRACER
integer :: i_river_NH4     = NO_TRACER

! ---- diagnostic field IDs
! unused:
integer ::  &
    id_nsoilcohorts, &
    id_lwc, id_swc, id_psi, id_temp, &
    id_ie, id_sn, id_bf, id_if, id_al, id_nu, id_sc, &
    id_hie, id_hsn, id_hbf, id_hif, id_hal, id_hnu, id_hsc, &
    id_heat_cap, id_thermal_cond, id_type, id_tau_gw, id_slope_l, &
    id_slope_Z, id_zeta_bar, id_e_depth, id_vwc_sat, id_vwc_fc, &
    id_vwc_wilt, id_K_sat, id_K_gw, id_w_fc, id_alpha, &
    id_refl_dry_dif, id_refl_dry_dir, id_refl_sat_dif, id_refl_sat_dir, &
    id_f_iso_dry, id_f_vol_dry, id_f_geo_dry, &
    id_f_iso_sat, id_f_vol_sat, id_f_geo_sat, &
    id_evap, id_uptk_n_iter, id_uptk, id_psi_x0, &
    id_excess, id_deficit, id_deficit_2, id_deficit_3, id_deficit_4, id_zeta, id_tau, &
    id_psi_bot, id_sat_frac, id_sat_depth, id_sat_dept2, &
    id_cf_1, id_cf_3, id_wt_1, id_wt_2, id_wt_2a, id_wt_2b, id_wt_3, id_wt2_3, id_wt_4, &
    id_div_bf, id_div_if, id_div_al, id_div, &
    id_z_cap, id_active_layer, id_surface_water, id_inun_frac, id_rsn_frac, id_flow, id_reflux, &
    id_protected_C_leaching, id_livemic_C_leaching, &
    id_protected_total_C, id_protected_total_N, &
    id_asoil,id_rsoil,&
    id_fast_DOC_div_loss,id_slow_DOC_div_loss,id_deadmic_DOC_div_loss, &
    id_fast_DON_div_loss,id_slow_DON_div_loss,id_deadmic_DON_div_loss, &
    id_wet_frac, id_macro_infilt, &
    id_surf_DOC_loss, id_total_C_leaching, id_total_DOC_div_loss, id_total_ON_leaching, id_NO3_leaching, id_NH4_leaching, &
    id_total_DON_div_loss, id_total_NO3_div_loss, id_total_NH4_div_loss, id_passive_N_uptake,&
    id_Qmax

integer :: &
    id_protected_C, id_livemic_total_C, id_deadmic_total_C, id_fsc, id_ssc, &
    id_protected_N, id_livemic_total_N, id_deadmic_total_N, id_fsn, id_ssn, &
    id_livemic_C, id_total_soil_C, id_dissolved_total_C, id_total_C_layered, &
    id_livemic_N, id_total_soil_N, id_dissolved_total_N, id_total_N_layered, &
    id_tile_N_gain, id_tile_N_loss, &
    id_negative_litter_C(N_C_TYPES), id_tot_negative_litter_C, &
    id_negative_litter_N(N_C_TYPES), id_tot_negative_litter_N

integer, dimension(N_LITTER_POOLS) :: id_nlittercohorts, &
    id_litter_livemic_C, id_litter_total_C, id_litter_total_C_leaching, id_litter_total_ON_leaching, id_litter_NO3_leaching, id_litter_NH4_leaching,&
    id_litter_livemic_N, id_litter_total_N, id_litter_nitrate, id_litter_ammonium
integer, dimension(N_C_TYPES) :: &
    id_soil_C,           id_soil_N, &
    id_soil_dissolved_C, id_soil_dissolved_N, &
    id_soil_protected_C, id_soil_protected_N, &
    id_rsoil_C,          id_rsoil_N, &
    id_C_leaching, id_DON_leaching
integer, dimension(N_LITTER_POOLS,N_C_TYPES) :: &
    id_litter_C, id_litter_N, id_litter_dissolved_C, id_litter_dissolved_N, &
    id_litter_protected_C, id_litter_protected_N, &
    id_litter_rsoil_C,     id_litter_rsoil_N, &
    id_litter_C_leaching, id_litter_DON_leaching

! FIXME: add N leaching terms to diagnostics?

integer :: &
    id_total_NH4,id_total_NO3,&
    id_soil_NO3,id_soil_NH4,&
    id_total_denitrification_rate,id_soil_denitrification_rate,id_NO3_div_loss,&
    id_NH4_div_loss,id_total_N_mineralization_rate,id_total_N_immobilization_rate,&
    id_total_nitrification_rate

! test tridiagonal solver for advection
integer :: id_st_diff

! diag IDs of CMOR variables
integer :: id_mrlsl, id_mrsfl, id_mrsll, id_mrsol, id_mrso, id_mrsos, id_mrlso, id_mrfso, &
    id_mrsofc, id_mrs1mLut, id_mrro, id_mrros, id_csoil, id_rh, &
    id_csoilfast, id_csoilmedium, id_csoilslow, id_cLitter, id_cLitterCwd

! variables for CMOR/CMIP diagnostic calculations
real, allocatable :: mrsos_weight(:) ! weights for mrsos averaging
real, allocatable :: mrs1m_weight(:) ! weights for mrs1m averaging

integer, allocatable :: soil_tags(:,:) ! module copy of soil tags for cold start

! ==== end of module variables ===============================================

contains

! ============================================================================
subroutine read_soil_namelist()
  ! ---- local vars
  integer :: unit         ! unit for namelist i/o
  integer :: io           ! i/o status for the namelist
  integer :: ierr         ! error code, returned by i/o routines
  integer :: l

  call read_soil_data_namelist(use_single_geo,gw_option)

  call log_version(version, module_name, &
  __FILE__)
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

  call soil_util_init(r_rhiz)
end subroutine read_soil_namelist


! ============================================================================
! initialize soil model
subroutine soil_init ( id_ug, id_band, id_zfull )
  integer,intent(in)  :: id_ug    !<Unstructured axis id.
  integer,intent(in)  :: id_band  !<ID of spectral band axis
  integer,intent(out) :: id_zfull !<ID of vertical soil axis

  ! ---- local vars
  type(land_tile_enum_type)     :: ce   ! tile list enumerator
  type(land_tile_type), pointer :: tile ! pointer to current tile
  ! input data buffers for respective variables:
  real, allocatable :: gw_param(:), gw_param2(:), gw_param3(:), albedo(:,:)
  real, allocatable :: f_iso(:,:), f_vol(:,:), f_geo(:,:), refl_dif(:,:)

  real :: total_soil_depth ! for diagnostics only, m
  real :: local_wt_depth ! [m] water table depth for tile (+ for below surface)
  real, allocatable :: ref_soil_t(:) ! reference soil temperature (based on 5 m or surface air temperature)
                                     ! for cold-start initialization
  real, allocatable :: wetmask(:)    ! input mask for zones with high water table
  logical :: drypoint                ! This point is predicted to have a falling water table.

  integer :: i, k, ll ! indices
  real :: psi(num_l), mwc(num_l)

  type(land_restart_type) :: restart, restart1
  logical :: restart_exists
  character(*), parameter :: restart_file_name = 'INPUT/soil.res.nc'

  module_is_initialized = .TRUE.
  delta_time = time_type_to_real(lnd%dt_fast)
  dt_fast_yr = delta_time/seconds_per_year

  call uptake_init(num_l,dz,zfull)
  call hlsp_hydro_lev_init(num_l,dz,zfull)

  ! initialize river tracer indices
  i_river_DOC  = river_tracer_index('doc')
  i_river_DON  = river_tracer_index('don')
  i_river_NO3  = river_tracer_index('no3')
  i_river_NH4  = river_tracer_index('nh4')

  if (i_river_DOC == NO_TRACER .and. (soil_carbon_option==SOILC_CORPSE .or. soil_carbon_option==SOILC_CORPSE_N)) &
      call error_mesg ('soil_init','River tracer for DOC not found: leached DOC goes directly to the atmosphere as CO2 to maintain carbon conservation.', NOTE)

  ! -------- initialize soil model diagnostic fields
  call soil_diag_init(id_ug,id_band,id_zfull)

  ! -------- read spatially distributed fields for groundwater parameters, if requested
  if (.not.use_single_geo) then
     select case (gw_option)
     case (GW_LINEAR,GW_LM2)
        allocate(gw_param(lnd%ls:lnd%le))
        call read_field( 'INPUT/groundwater_residence.nc','tau', gw_param, interp='bilinear' )
        call put_to_tiles_r0d_fptr( gw_param, land_tile_map, soil_tau_groundwater_ptr )
        deallocate(gw_param)
     case (GW_HILL, GW_HILL_AR5)
        allocate(gw_param (lnd%ls:lnd%le))
        allocate(gw_param2(lnd%ls:lnd%le))
        allocate(gw_param3(lnd%ls:lnd%le))
        call read_field( 'INPUT/geohydrology.nc','hillslope_length', gw_param, interp='bilinear' )
        call put_to_tiles_r0d_fptr( gw_param*gw_scale_length, land_tile_map, soil_hillslope_length_ptr )
        call read_field( 'INPUT/geohydrology.nc','slope', gw_param2, interp='bilinear' )
        gw_param = gw_param*gw_param2
        call put_to_tiles_r0d_fptr( gw_param*gw_scale_relief, land_tile_map, soil_hillslope_relief_ptr )

        if (retro_a0n1 .or. gw_option.eq.GW_HILL_AR5) then
            gw_param = 0.
            call put_to_tiles_r0d_fptr( gw_param, land_tile_map, soil_hillslope_a_ptr )
            gw_param = 1.
            call put_to_tiles_r0d_fptr( gw_param, land_tile_map, soil_hillslope_n_ptr )
!            call read_field( 'INPUT/geohydrology.nc','hillslope_zeta_bar', &
!              lnd%sg_lon, lnd%sg_lat, gw_param, interp='bilinear' )
            gw_param = 0.5
            call put_to_tiles_r0d_fptr( gw_param, land_tile_map, soil_hillslope_zeta_bar_ptr )
        else
            call read_field( 'INPUT/geohydrology.nc','hillslope_a', gw_param, interp='bilinear' )
            call put_to_tiles_r0d_fptr( gw_param, land_tile_map, soil_hillslope_a_ptr )
            call read_field( 'INPUT/geohydrology.nc','hillslope_n', gw_param2, interp='bilinear' )
            call put_to_tiles_r0d_fptr( gw_param2, land_tile_map, soil_hillslope_n_ptr )
            gw_param3 = (1./(gw_param2+1.)+gw_param/(gw_param2+2.))/(1.+gw_param/2.)
            call put_to_tiles_r0d_fptr( gw_param3, land_tile_map, soil_hillslope_zeta_bar_ptr )
        endif

        call read_field( 'INPUT/geohydrology.nc','soil_e_depth', gw_param, interp='bilinear' )
        if (slope_exp.gt.0.01) then
            call put_to_tiles_r0d_fptr( gw_param*gw_scale_soil_depth*(0.08/gw_param2)**slope_exp, &
                                                  land_tile_map, soil_soil_e_depth_ptr )
        else
            call put_to_tiles_r0d_fptr( gw_param*gw_scale_soil_depth, land_tile_map, soil_soil_e_depth_ptr )
        endif
        if (gw_option /= GW_HILL_AR5) then
            call read_field( 'INPUT/geohydrology.nc','perm', gw_param, interp='bilinear' )
            call put_to_tiles_r0d_fptr(9.8e9*gw_scale_perm*gw_param, land_tile_map, &
                                            soil_k_sat_gw_ptr )
        endif
        deallocate(gw_param, gw_param2, gw_param3)
        ce = first_elmt(land_tile_map)
        do while(loop_over_tiles(ce,tile))
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
           allocate(gw_param (lnd%ls:lnd%le), gw_param2(lnd%ls:lnd%le))
           call read_field( 'INPUT/geohydrology.nc','hillslope_length', gw_param, interp='bilinear' )
           call put_to_tiles_r0d_fptr( gw_param*gw_scale_length, land_tile_map, soil_hillslope_length_ptr )
           call read_field( 'INPUT/geohydrology.nc','slope', gw_param2, interp='bilinear' )
           gw_param = gw_param*gw_param2
           call put_to_tiles_r0d_fptr( gw_param*gw_scale_relief, land_tile_map, soil_hillslope_relief_ptr )
           call read_field( 'INPUT/geohydrology.nc','hillslope_zeta_bar', gw_param, interp='bilinear' )
           if (zeta_bar_override.gt.0.) gw_param=zeta_bar_override
           call put_to_tiles_r0d_fptr( gw_param, land_tile_map, soil_hillslope_zeta_bar_ptr )
           call read_field( 'INPUT/geohydrology.nc','soil_e_depth', gw_param, interp='bilinear' )

           if (slope_exp.gt.0.01) then
           ! ZMS It's probably inconsistent to leave in this if statement.
               call error_mesg(module_name, 'soil_init: "slope_exp" > 0.0 requested even though '// &
                               'running with tiled hillslope model.  This may be inconsistent.', WARNING)
               call put_to_tiles_r0d_fptr( gw_param*gw_scale_soil_depth*(0.08/gw_param2)**slope_exp, &
                                                     land_tile_map, soil_soil_e_depth_ptr )
           else
               call put_to_tiles_r0d_fptr( gw_param*gw_scale_soil_depth, land_tile_map, soil_soil_e_depth_ptr )
           endif
           call read_field( 'INPUT/geohydrology.nc','perm', gw_param, interp='bilinear' )
           call put_to_tiles_r0d_fptr(9.8e9*gw_scale_perm*gw_param, land_tile_map, &
                                          soil_k_sat_gw_ptr )
           deallocate(gw_param, gw_param2)
        end if
        ce = first_elmt(land_tile_map)
        do while(loop_over_tiles(ce,tile))
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
     allocate(albedo(lnd%ls:lnd%le,NBANDS))
     call read_field( 'INPUT/soil_albedo.nc','SOIL_ALBEDO_VIS', albedo(:,BAND_VIS),'bilinear')
     call read_field( 'INPUT/soil_albedo.nc','SOIL_ALBEDO_NIR', albedo(:,BAND_NIR),'bilinear')
     call put_to_tiles_r1d_fptr( albedo, land_tile_map, soil_refl_dry_dir_ptr )
     call put_to_tiles_r1d_fptr( albedo, land_tile_map, soil_refl_dry_dif_ptr )
     ! for now, put the same value into the saturated soil albedo, so that
     ! the albedo does not depend on soil wetness
     call put_to_tiles_r1d_fptr( albedo, land_tile_map, soil_refl_sat_dir_ptr )
     call put_to_tiles_r1d_fptr( albedo, land_tile_map, soil_refl_sat_dif_ptr )
     deallocate(albedo)
  else if (trim(albedo_to_use)=='brdf-maps') then
     use_brdf = .true.
     allocate(   f_iso(lnd%ls:lnd%le,NBANDS))
     allocate(   f_vol(lnd%ls:lnd%le,NBANDS))
     allocate(   f_geo(lnd%ls:lnd%le,NBANDS))
     allocate(refl_dif(lnd%ls:lnd%le,NBANDS))
     call read_field( 'INPUT/soil_brdf.nc','f_iso_vis', f_iso(:,BAND_VIS),'bilinear')
     call read_field( 'INPUT/soil_brdf.nc','f_vol_vis', f_vol(:,BAND_VIS),'bilinear')
     call read_field( 'INPUT/soil_brdf.nc','f_geo_vis', f_geo(:,BAND_VIS),'bilinear')
     call read_field( 'INPUT/soil_brdf.nc','f_iso_nir', f_iso(:,BAND_NIR),'bilinear')
     call read_field( 'INPUT/soil_brdf.nc','f_vol_nir', f_vol(:,BAND_NIR),'bilinear')
     call read_field( 'INPUT/soil_brdf.nc','f_geo_nir', f_geo(:,BAND_NIR),'bilinear')
     refl_dif = g_iso*f_iso + g_vol*f_vol + g_geo*f_geo
     call put_to_tiles_r1d_fptr( f_iso,    land_tile_map, soil_f_iso_dry_ptr )
     call put_to_tiles_r1d_fptr( f_vol,    land_tile_map, soil_f_vol_dry_ptr )
     call put_to_tiles_r1d_fptr( f_geo,    land_tile_map, soil_f_geo_dry_ptr )
     call put_to_tiles_r1d_fptr( refl_dif, land_tile_map, soil_refl_dry_dif_ptr )
     ! for now, put the same value into the saturated soil albedo, so that
     ! the albedo does not depend on soil wetness
     call put_to_tiles_r1d_fptr( f_iso,    land_tile_map, soil_f_iso_sat_ptr )
     call put_to_tiles_r1d_fptr( f_vol,    land_tile_map, soil_f_vol_sat_ptr )
     call put_to_tiles_r1d_fptr( f_geo,    land_tile_map, soil_f_geo_sat_ptr )
     call put_to_tiles_r1d_fptr( refl_dif, land_tile_map, soil_refl_sat_dif_ptr )
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
     allocate(ref_soil_t(lnd%ls:lnd%le), wetmask(lnd%ls:lnd%le))
     call read_field( coldstart_datafile, 'REFSOILT', ref_soil_t, interp='bilinear' )
     call read_field( coldstart_datafile, 'WETMASK', wetmask, interp='bilinear' )
  end if

  ! -------- initialize soil state --------
  ce = first_elmt(land_tile_map, ls=lnd%ls) ! Use global indices here because element indices
                                            ! needed.
  do while(loop_over_tiles(ce,tile,k=k,l=ll))
     if (.not.associated(tile%soil)) cycle
     call set_current_point(ll,k)
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
           if (wetmask(ll) > 0.5) then ! wet point
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
        tile%soil%uptake_T      = init_temp
     else
        call init_soil_twc(tile%soil, ref_soil_t(ll), mwc)
     end if
  end do

  if (use_coldstart_wtt_data) then
     deallocate(ref_soil_t, wetmask)
  end if

  call open_land_restart(restart,restart_file_name,restart_exists)
  if (restart_exists) then
     call error_mesg('soil_init', 'reading NetCDF restart "'//trim(restart_file_name)//'"', NOTE)
     call get_tile_data(restart, 'temp', 'zfull', soil_T_ptr  )
     call get_tile_data(restart, 'wl', 'zfull', soil_wl_ptr )
     call get_tile_data(restart, 'ws', 'zfull', soil_ws_ptr )
     call get_tile_data(restart, 'groundwater', 'zfull', soil_groundwater_ptr )
     call get_tile_data(restart, 'groundwater_T', 'zfull', soil_groundwater_T_ptr)
     if(field_exists(restart, 'uptake_T')) &
          call get_tile_data(restart, 'uptake_T', soil_uptake_T_ptr)

     if (soil_carbon_option==SOILC_CENTURY.or.soil_carbon_option==SOILC_CENTURY_BY_LAYER) then
        if(field_exists(restart, 'fsc')) then
           call get_tile_data(restart,'fsc','zfull',soil_fast_soil_C_ptr)
           call get_tile_data(restart,'ssc','zfull',soil_slow_soil_C_ptr)
        else
           ! try to read fsc and ssc from vegetation restart
           call open_land_restart(restart1,'INPUT/vegn2.res.nc',restart_exists)
           if (restart_exists) then
              ! read old (scalar) fsc and ssc into the first element of the fast_soil_C
              ! and slow_soil_C arrays
              call get_tile_data(restart1,'fsc',soil_fast_soil_C_ptr,1)
              call get_tile_data(restart1,'ssc',soil_slow_soil_C_ptr,1)
           endif
           call free_land_restart(restart1)
        endif
     endif

     if (field_exists(restart,'fast_soil_C')) then
        ! we are dealing with CORPSE restart
        ce = first_elmt(land_tile_map)
        do while(loop_over_tiles(ce,tile))
            if (.not.associated(tile%soil)) cycle
            do i = 1,N_LITTER_POOLS
               call adjust_pool_ncohorts(tile%soil%litter(i))
            enddo
            do i = 1,num_l
               call adjust_pool_ncohorts(tile%soil%soil_organic_matter(i))
            enddo
        end do
        do i = 1, N_C_TYPES
           call get_tile_data(restart,trim(c_shortname(i))//'_soil_C', 'zfull','soilCCohort', sc_soil_C_ptr,i)
           call get_tile_data(restart,trim(c_shortname(i))//'ProtectedC', 'zfull','soilCCohort', sc_protected_C_ptr,i)
           call get_tile_data(restart,'soil_DOC_'//trim(c_shortname(i)), 'zfull', sc_DOC_ptr,i)

           do k = 1, N_LITTER_POOLS
              call get_tile_data(restart,trim(l_shortname(k))//'_litter_'//trim(c_shortname(i))//'_C','litterCCohort',sc_litter_litterC_ptr,i,k)
              call get_tile_data(restart,trim(l_shortname(k))//'_litter_'//trim(c_shortname(i))//'ProtectedC','litterCCohort',sc_litter_protectedC_ptr,i,k)
              call get_tile_data(restart,trim(l_shortname(k))//'_litter_DOC_'//trim(c_shortname(i)),sc_litter_dissolved_carbon_ptr,i,k)
           enddo
        enddo
        call get_tile_data(restart,'liveMic', 'zfull','soilCCohort',sc_livingMicrobeC_ptr)
        call get_tile_data(restart,'CO2', 'zfull','soilCCohort',sc_CO2_ptr)
        call get_tile_data(restart,'Rtot', 'zfull','soilCCohort',sc_Rtot_ptr)
        call get_tile_data(restart,'originalCohortC', 'zfull','soilCCohort', sc_originalLitterC_ptr)

        do i = 1,N_LITTER_POOLS
           call get_tile_data(restart, trim(l_shortname(i))//'_litter_liveMic_C', 'litterCCohort', sc_litter_livingMicrobeC_ptr, i)
           call get_tile_data(restart, trim(l_shortname(i))//'_litter_CO2',       'litterCCohort', sc_litter_CO2_ptr, i)
           call get_tile_data(restart, trim(l_shortname(i))//'_litter_Rtot',      'litterCCohort', sc_litter_Rtot_ptr, i)
           call get_tile_data(restart, trim(l_shortname(i))//'_litter_originalCohortC', 'litterCCohort',sc_litter_originalLitterC_ptr, i)
        enddo

        if(field_exists(restart, 'fast_DOC_leached')) then
           call get_tile_data(restart,'fast_DOC_leached', soil_fast_DOC_leached_ptr)
           call get_tile_data(restart,'slow_DOC_leached', soil_slow_DOC_leached_ptr)
           call get_tile_data(restart,'deadmic_DOC_leached', soil_deadmic_DOC_leached_ptr)
        endif

        if(field_exists(restart, 'gross_nitrogen_flux_into_tile')) then
           call get_tile_data(restart,'gross_nitrogen_flux_into_tile', soil_gross_nitrogen_flux_into_tile_ptr)
           call get_tile_data(restart,'gross_nitrogen_flux_out_of_tile', soil_gross_nitrogen_flux_out_of_tile_ptr)
         endif

        if(field_exists(restart, 'is_peat')) then
           call get_int_tile_data(restart, 'is_peat','zfull', soil_is_peat_ptr)
        endif
     endif
     if (field_exists(restart,'fast_soil_N')) then
        do i = 1, N_C_TYPES
           call get_tile_data(restart,trim(c_shortname(i))//'_soil_N', 'zfull','soilCCohort', sc_soil_N_ptr,i)
           call get_tile_data(restart,trim(c_shortname(i))//'ProtectedN', 'zfull','soilCCohort', sc_protected_N_ptr,i)
           call get_tile_data(restart,'soil_DON_'//trim(c_shortname(i)), 'zfull', sc_DON_ptr,i)

           do k = 1, N_LITTER_POOLS
              call get_tile_data(restart,trim(l_shortname(k))//'_litter_'//trim(c_shortname(i))//'_N','litterCCohort',sc_litter_litterN_ptr,i,k)
              call get_tile_data(restart,trim(l_shortname(k))//'_litter_'//trim(c_shortname(i))//'ProtectedN','litterCCohort',sc_litter_protectedN_ptr,i,k)
              call get_tile_data(restart,trim(l_shortname(k))//'_litter_DON_'//trim(c_shortname(i)),sc_litter_dissolved_nitrogen_ptr,i,k)
           enddo
        enddo
        call get_tile_data(restart,'liveMicN', 'zfull','soilCCohort', sc_livingMicrobeN_ptr)
        call get_tile_data(restart,'originalCohortN', 'zfull','soilCCohort', sc_originalLitterN_ptr)
        call get_tile_data(restart,'soil_NO3', 'zfull', sc_nitrate_ptr)
        call get_tile_data(restart,'soil_NH4', 'zfull', sc_ammonium_ptr)
        call get_tile_data(restart,'soil_nitrif', 'zfull', sc_nitrif_ptr)
        call get_tile_data(restart,'soil_denitrif', 'zfull', sc_denitrif_ptr)

        ! Leaving out cohort-level immobilization and mineralization fields for now -- BNS

        do k = 1,N_LITTER_POOLS
           call get_tile_data(restart, trim(l_shortname(k))//'_litter_liveMic_N', 'litterCCohort', sc_litter_livingMicrobeN_ptr,k)
           call get_tile_data(restart, trim(l_shortname(k))//'_litter_originalCohortN', 'litterCCohort', sc_litter_originalLitterN_ptr,k)
           call get_tile_data(restart, trim(l_shortname(k))//'_litter_NO3', sc_litter_nitrate_ptr,k)
           call get_tile_data(restart, trim(l_shortname(k))//'_litter_NH4', sc_litter_ammonium_ptr,k)
           call get_tile_data(restart, trim(l_shortname(k))//'_litter_nitrif', sc_litter_nitrif_ptr,k)
           call get_tile_data(restart, trim(l_shortname(k))//'_litter_denitrif', sc_litter_denitrif_ptr,k)
        enddo
     endif
     do i = 1, N_C_TYPES
        if(field_exists(restart, 'negative_litter_C_'//trim(c_shortname(i)))) then
           call get_tile_data(restart,'negative_litter_C_'//trim(c_shortname(i)),sc_negative_litter_C_ptr,i)
        endif
        if(field_exists(restart, 'negative_litter_N_'//trim(c_shortname(i)))) then
           call get_tile_data(restart,'negative_litter_N_'//trim(c_shortname(i)),sc_negative_litter_N_ptr,i)
        endif
     enddo
  else
     call error_mesg('soil_init', 'cold-starting soil', NOTE)
  endif
  call free_land_restart(restart)

  ! read soil carbon restart, if present
  call open_land_restart(restart,'INPUT/soil_carbon.res.nc',restart_exists)
  if (restart_exists) then
     call error_mesg('veg_data_init','reading soil_carbon restart',NOTE)
     select case (soil_carbon_option)
     case (SOILC_CENTURY, SOILC_CENTURY_BY_LAYER)
        call get_tile_data(restart,'asoil_in','zfull',soil_asoil_in_ptr)
        call get_tile_data(restart,'fsc_in','zfull',soil_fsc_in_ptr)
        call get_tile_data(restart,'ssc_in','zfull',soil_ssc_in_ptr)
     case (SOILC_CORPSE, SOILC_CORPSE_N)
        call get_tile_data(restart,'asoil_in','zfull',soil_asoil_in_ptr)
        do i = 1,N_C_TYPES
           ! C inputs
           call get_tile_data(restart,trim(c_shortname(i))//'_protected_C_in','zfull',sc_protected_C_in_ptr, i)
           call get_tile_data(restart,trim(c_shortname(i))//'_C_in','zfull',sc_C_in_ptr, i)
           do k = 1,N_LITTER_POOLS
              call get_tile_data(restart,trim(l_shortname(k))//'litter_'//trim(c_shortname(i))//'_C_in',sc_litter_C_in_ptr,i,k)
              call get_tile_data(restart,trim(l_shortname(k))//'litter_'//trim(c_shortname(i))//'_N_in',sc_litter_N_in_ptr,i,k)
              call get_tile_data(restart,trim(l_shortname(k))//'litter_'//trim(c_shortname(i))//'_C_turnover_accumulated', sc_litter_C_turnover_ptr, i,k)
              call get_tile_data(restart,trim(l_shortname(k))//'litter_'//trim(c_shortname(i))//'_N_turnover_accumulated', sc_litter_N_turnover_ptr, i,k)
           enddo
           ! N inputs
           call get_tile_data(restart,trim(c_shortname(i))//'_protected_N_in','zfull',sc_protected_N_in_ptr, i)
           call get_tile_data(restart,trim(c_shortname(i))//'_N_in','zfull',sc_N_in_ptr, i)
           ! C turnover rates
           call get_tile_data(restart,trim(c_shortname(i))//'_C_turnover_accumulated','zfull',sc_C_turnover_ptr, i)
           call get_tile_data(restart,trim(c_shortname(i))//'_protected_C_turnover_accumulated','zfull',sc_protected_C_turnover_ptr, i)
           ! N turnover rates
           call get_tile_data(restart,trim(c_shortname(i))//'_N_turnover_accumulated','zfull',sc_N_turnover_ptr, i)
           call get_tile_data(restart,trim(c_shortname(i))//'_protected_N_turnover_accumulated','zfull',sc_protected_N_turnover_ptr, i)
        enddo
     case default
        call error_mesg('save_soil_restart', 'unrecognized soil carbon option -- this should never happen', FATAL)
     end select
  endif
  call free_land_restart(restart)

  ! ---- static diagnostic section
  call send_tile_data_r0d_fptr(id_tau_gw,       soil_tau_groundwater_ptr)
  call send_tile_data_r0d_fptr(id_slope_l,      soil_hillslope_length_ptr)
  call send_tile_data_r0d_fptr(id_slope_Z,      soil_hillslope_relief_ptr)
  call send_tile_data_r0d_fptr(id_zeta_bar,     soil_hillslope_zeta_bar_ptr)
  call send_tile_data_r0d_fptr(id_e_depth,      soil_soil_e_depth_ptr)
  call send_tile_data_r0d_fptr(id_zeta,         soil_zeta_ptr)
  call send_tile_data_r0d_fptr(id_tau,          soil_tau_ptr)
  call send_tile_data_r0d_fptr(id_vwc_wilt,     soil_vwc_wilt_ptr)
  call send_tile_data_r0d_fptr(id_vwc_fc,       soil_vwc_fc_ptr)
  call send_tile_data_r0d_fptr(id_vwc_sat,      soil_vwc_sat_ptr)
  call send_tile_data_r0d_fptr(id_K_sat,        soil_k_sat_ref_ptr)
  call send_tile_data_r0d_fptr(id_K_gw,         soil_k_sat_gw_ptr)
  call send_tile_data_r0d_fptr(id_Qmax,         soil_Qmax_ptr)
  call send_tile_data_r1d_fptr(id_w_fc,         soil_w_fc_ptr)
  call send_tile_data_r1d_fptr(id_alpha,        soil_alpha_ptr)
  call send_tile_data_r1d_fptr(id_refl_dry_dir, soil_refl_dry_dir_ptr)
  call send_tile_data_r1d_fptr(id_refl_dry_dif, soil_refl_dry_dif_ptr)
  call send_tile_data_r1d_fptr(id_refl_sat_dir, soil_refl_sat_dir_ptr)
  call send_tile_data_r1d_fptr(id_refl_sat_dif, soil_refl_sat_dif_ptr)
  call send_tile_data_r1d_fptr(id_f_iso_dry,    soil_f_iso_dry_ptr)
  call send_tile_data_r1d_fptr(id_f_vol_dry,    soil_f_vol_dry_ptr)
  call send_tile_data_r1d_fptr(id_f_geo_dry,    soil_f_geo_dry_ptr)
  call send_tile_data_r1d_fptr(id_f_iso_sat,    soil_f_iso_sat_ptr)
  call send_tile_data_r1d_fptr(id_f_vol_sat,    soil_f_vol_sat_ptr)
  call send_tile_data_r1d_fptr(id_f_geo_sat,    soil_f_geo_sat_ptr)
  call send_tile_data_i0d_fptr(id_type,         soil_tag_ptr)
  if (id_mrsofc>0) then
     total_soil_depth = sum(dz(1:num_l))
     ce = first_elmt(land_tile_map)
     do while(loop_over_tiles(ce,tile))
        if (associated(tile%soil)) &
           call send_tile_data(id_mrsofc, tile%soil%pars%vwc_sat*total_soil_depth*dens_h2o, tile%diag)
     end do
  endif
end subroutine soil_init


! ============================================================================
function replace_text (s,text,rep)  result(outs)
character(*), intent(in) :: s,text,rep
character(len(s)+100) :: outs     ! provide outs with extra 100 char len

integer      :: i, nt, nr

outs = s ; nt = len_trim(text) ; nr = len_trim(rep)
do
   i = index(outs,text(:nt)) ; if (i == 0) exit
   outs = outs(:i-1) // rep(:nr) // outs(i+nt:)
end do
end function replace_text

! ============================================================================
function register_soilc_diag_fields(module_name, field_name, axes, init_time, &
     long_name, units, missing_value, range, op, standard_name) result (id)

  integer :: id(N_C_TYPES)

  character(len=*), intent(in) :: module_name
  character(len=*), intent(in) :: field_name
  integer,          intent(in) :: axes(:)
  type(time_type),  intent(in) :: init_time
  character(len=*), intent(in), optional :: long_name
  character(len=*), intent(in), optional :: units
  real,             intent(in), optional :: missing_value
  real,             intent(in), optional :: range(2)
  character(len=*), intent(in), optional :: op ! aggregation operation
  character(len=*), intent(in), optional :: standard_name

  integer :: i

  do i = 1, N_C_TYPES
     id(i) = register_tiled_diag_field(module_name, &
             trim(replace_text(field_name,'<ctype>',trim(c_shortname(i)))), &
             axes, init_time, &
             trim(replace_text(long_name,'<ctype>',trim(c_longname(i)))), &
             units, missing_value, range, op, standard_name)
  enddo
end function register_soilc_diag_fields

! ============================================================================
! registered an array of diag fields, one per litter pool
function register_litter_diag_fields(module_name, field_name, axes, init_time, &
     long_name, units, missing_value, range, op, standard_name) result (id)

  integer :: id(N_LITTER_POOLS)

  character(len=*), intent(in) :: module_name
  character(len=*), intent(in) :: field_name
  integer,          intent(in) :: axes(:)
  type(time_type),  intent(in) :: init_time
  character(len=*), intent(in), optional :: long_name
  character(len=*), intent(in), optional :: units
  real,             intent(in), optional :: missing_value
  real,             intent(in), optional :: range(2)
  character(len=*), intent(in), optional :: op ! aggregation operation
  character(len=*), intent(in), optional :: standard_name

  integer :: i

  do i = 1, N_LITTER_POOLS
     id(i) = register_tiled_diag_field(module_name, &
             trim(replace_text(field_name,'<ltype>',trim(l_shortname(i)))), &
             axes, init_time, &
             trim(replace_text(long_name,'<ltype>',trim(l_longname(i)))), &
             units, missing_value, range, op, standard_name)
  enddo
end function register_litter_diag_fields

! ============================================================================
! registered a 2D array of diag fields, one per litter pool per carbon type
function register_litter_soilc_diag_fields(module_name, field_name, axes, init_time, &
     long_name, units, missing_value, range, op, standard_name) result (id)

  integer :: id(N_LITTER_POOLS, N_C_TYPES)

  character(len=*), intent(in) :: module_name
  character(len=*), intent(in) :: field_name
  integer,          intent(in) :: axes(:)
  type(time_type),  intent(in) :: init_time
  character(len=*), intent(in), optional :: long_name
  character(len=*), intent(in), optional :: units
  real,             intent(in), optional :: missing_value
  real,             intent(in), optional :: range(2)
  character(len=*), intent(in), optional :: op ! aggregation operation
  character(len=*), intent(in), optional :: standard_name

  integer :: i, k
  character(128) :: name
  character(512) :: lname

  do i = 1, N_C_TYPES
     do k = 1, N_LITTER_POOLS
        name = replace_text(field_name,'<ctype>',trim(c_shortname(i)))
        name = replace_text(name,      '<ltype>',trim(l_shortname(k)))
        lname = replace_text(long_name,'<ctype>',trim(c_longname(i)))
        lname = replace_text(lname,    '<ltype>',trim(l_longname(k)))
        id(k,i) = register_tiled_diag_field(module_name, trim(name), axes, init_time, trim(lname), &
             units, missing_value, range, op, standard_name)
     enddo
  enddo
end function register_litter_soilc_diag_fields

! ============================================================================
subroutine soil_diag_init(id_ug,id_band,id_zfull)
  integer,intent(in)  :: id_ug    !<Unstructured axis id.
  integer,intent(in)  :: id_band  !<ID of spectral band axis
  integer,intent(out) :: id_zfull !<ID of vertical soil axis

  ! ---- local vars
  integer :: axes(2)
  integer :: id_zhalf
  integer :: i, l

  ! define vertical axis and its edges
  id_zhalf = diag_axis_init ( &
       'zhalf_soil', zhalf(1:num_l+1), 'meters', 'z', 'half level',  -1, set_name='soil' )
  id_zfull = diag_axis_init ( &
       'zfull_soil', zfull(1:num_l),   'meters', 'z', 'full level',  -1, set_name='soil', &
       edges=id_zhalf )

  ! define array of axis indices
  axes = (/id_ug,id_zfull/)

  ! set the default sub-sampling filter for the fields below
  call set_default_diag_filter('soil')

  ! define diagnostic fields

  ! FIXME slm: generalize DOC fields by carbon type
  id_fast_DOC_div_loss = register_tiled_diag_field ( module_name, 'fast_DOC_div_loss', (/id_ug/),  &
       lnd%time, 'total fast DOC divergence loss', 'kg C/m2', missing_value=-100.0 )
  id_slow_DOC_div_loss = register_tiled_diag_field ( module_name, 'slow_DOC_div_loss', (/id_ug/),  &
       lnd%time, 'total slow DOC divergence loss', 'kg C/m2', missing_value=-100.0 )
  id_deadmic_DOC_div_loss = register_tiled_diag_field ( module_name, 'deadmic_DOC_div_loss', (/id_ug/),  &
       lnd%time, 'total dead microbe DOC divergence loss', 'kg C/m2', missing_value=-100.0 )
  id_total_DOC_div_loss = register_tiled_diag_field ( module_name, 'total_DOC_div', axes(1:1), &
       lnd%time, 'total rate of DOC divergence loss', 'kg C/m^2/s', missing_value=initval)
  id_total_DON_div_loss = register_tiled_diag_field ( module_name, 'total_DON_div', axes(1:1), &
       lnd%time, 'total rate of DON divergence loss', 'kg N/m^2/s', missing_value=initval)
  id_total_NO3_div_loss = register_tiled_diag_field ( module_name, 'total_NO3_div', axes(1:1), &
       lnd%time, 'total rate of NO3 divergence loss', 'kg N/m^2/s', missing_value=initval)
  id_total_NH4_div_loss = register_tiled_diag_field ( module_name, 'total_NH4_div', axes(1:1), &
       lnd%time, 'total rate of NH4 divergence loss', 'kg N/m^2/s', missing_value=initval)
   id_passive_N_uptake = register_cohort_diag_field ( 'vegn', 'passive_N_uptake',  &
        (/id_ug/), lnd%time, 'Plant N uptake by root water flow', 'kg N/m2/year', missing_value=-1.0 )
id_div = register_tiled_diag_field(module_name, 'div',axes,lnd%time,'Water divergence rate by layer','kg/m2/s',missing_value=-100.0)
       id_fast_DON_div_loss = register_tiled_diag_field ( module_name, 'fast_DON_div_loss', (/id_ug/),  &
            lnd%time, 'total fast DON divergence loss', 'kg N/m2', missing_value=-100.0 )
       id_slow_DON_div_loss = register_tiled_diag_field ( module_name, 'slow_DON_div_loss', (/id_ug/),  &
            lnd%time, 'total slow DON divergence loss', 'kg N/m2', missing_value=-100.0 )
       id_deadmic_DON_div_loss = register_tiled_diag_field ( module_name, 'deadmic_DON_div_loss', (/id_ug/),  &
            lnd%time, 'total dead microbe DON divergence loss', 'kg N/m2', missing_value=-100.0 )

  id_rsoil = register_tiled_diag_field ( module_name, 'rsoil',  &
       axes(1:1), lnd%time, 'soil respiration', 'kg C/(m2 year)', missing_value=-100.0 )

  id_protected_C = register_tiled_diag_field ( module_name, 'protected_soil_C', axes,  &
       lnd%time, 'protected soil carbon per layer', 'kg C/m3', missing_value=-100.0 )
  id_protected_N = register_tiled_diag_field ( module_name, 'protected_soil_N', axes,  &
       lnd%time, 'protected soil nitrogen per layer', 'kg N/m3', missing_value=-100.0 )
  id_total_C_layered = register_tiled_diag_field ( module_name, 'total_soil_C_layered', &
       axes, lnd%time, 'total soil carbon', 'kg C/m3', missing_value=-100.0 )
  id_total_N_layered = register_tiled_diag_field ( module_name, 'total_soil_N_layered', axes,  &
       lnd%time, 'total soil nitrogen', 'kg N/m3', missing_value=-100.0 )

  ! by-carbon-species diag fields
  id_soil_C(:) = register_soilc_diag_fields(module_name, '<ctype>_soil_C', &
       axes, lnd%time, '<ctype> soil carbon', 'kg C/m3', missing_value=-100.0 )
  id_soil_dissolved_C(:) = register_soilc_diag_fields(module_name, '<ctype>_dissolved_C', &
       axes, lnd%time, '<ctype> dissolved soil carbon', 'kg C/m3', missing_value=-100.0 )
  id_soil_protected_C(:) = register_soilc_diag_fields(module_name, '<ctype>_protected_C', &
       axes, lnd%time, '<ctype> protected soil carbon', 'kg C/m3', missing_value=-100.0 )
  id_rsoil_C(:) = register_soilc_diag_fields(module_name, 'rsoil_<ctype>', &
       axes, lnd%time, '<ctype> soil carbon respiration', 'kg C/(m3 year)', missing_value=-100.0 )

  id_C_leaching(:) = register_soilc_diag_fields ( module_name, '<ctype>_C_leaching', axes, &
       lnd%time, 'net layer <ctype> soil C leaching',  'kg/(m2 s)', missing_value=-100.0)
  id_DON_leaching(:) = register_soilc_diag_fields ( module_name, '<ctype>_DON_leaching', axes, &
       lnd%time, 'net layer <ctype> soil DON leaching',  'kg/(m2 s)', missing_value=-100.0)

  id_soil_N(:) = register_soilc_diag_fields(module_name, '<ctype>_soil_N', &
       axes, lnd%time, '<ctype> soil nitrogen', 'kg N/m3', missing_value=-100.0 )
  id_soil_dissolved_N(:) = register_soilc_diag_fields(module_name, '<ctype>_dissolved_N', &
       axes, lnd%time, '<ctype> dissolved soil nitrogen', 'kg N/m3', missing_value=-100.0 )
  id_soil_protected_N(:) = register_soilc_diag_fields(module_name, '<ctype>_protected_N', &
       axes, lnd%time, '<ctype> protected soil nitrogen', 'kg N/m3', missing_value=-100.0 )
  id_rsoil_N(:) = register_soilc_diag_fields(module_name, 'rsoil_N_<ctype>', &
       axes, lnd%time, '<ctype> soil nitrogen respiration', 'kg N/(m3 year)', missing_value=-100.0 )

  ! litter fields
  id_litter_C(:,:) = register_litter_soilc_diag_fields ( module_name, '<ctype>_<ltype>litter_C', &
       axes(1:1), lnd%time, '<ctype> <ltype> litter carbon', 'kg C/m2', missing_value=-100.0 )
  id_litter_dissolved_C(:,:) = register_litter_soilc_diag_fields ( module_name, '<ctype>_<ltype>litter_dissolved_C', &
       axes(1:1), lnd%time, '<ctype> <ltype> litter dissolved carbon', 'kg C/m2', missing_value=-100.0 )
  id_litter_protected_C(:,:) = register_litter_soilc_diag_fields ( module_name, '<ctype>_<ltype>litter_protected_C', &
       axes(1:1), lnd%time, '<ctype> <ltype> litter protected carbon', 'kg C/m2', missing_value=-100.0 )
  id_litter_C_leaching(:,:) = register_litter_soilc_diag_fields ( module_name, '<ctype>_<ltype>litter_C_leaching', &
       axes(1:1), lnd%time, '<ltype> litter <ctype> C leaching','kg/(m2 s)', missing_value=-100.0)
  id_litter_rsoil_C(:,:) = register_litter_soilc_diag_fields ( module_name, 'rsoil_<ltype>litter_<ctype>',  &
       axes(1:1), lnd%time, 'surface <ltype> litter <ctype> carbon degradation', 'kg C/(m2 year)', missing_value=-100.0 )
  id_litter_livemic_C(:) = register_litter_diag_fields ( module_name, '<ltype>litter_live_microbe_C', &
       axes(1:1),  lnd%time, 'live microbe <ltype> litter carbon', 'kg C/m2', missing_value=-100.0 )
  id_litter_total_C(:) = register_litter_diag_fields ( module_name, '<ltype>litter_total_C', &
       axes(1:1),  lnd%time, '<ltype> litter total carbon', 'kg C/m2', missing_value=-100.0 )
  id_nlittercohorts(:) = register_litter_diag_fields ( module_name, 'n_<ltype>litter_cohorts', axes(1:1),  &
       lnd%time, 'number of <ltype> litter cohorts', missing_value=-100.0 )

  id_litter_DON_leaching(:,:) = register_litter_soilc_diag_fields ( module_name, '<ctype>_<ltype>litter_DON_leaching', &
       axes(1:1), lnd%time, '<ltype> litter <ctype> DON leaching','kg/(m2 s)', missing_value=-100.0)
  id_litter_total_C_leaching(:) = register_litter_diag_fields ( module_name, '<ltype>litter_total_C_leaching', &
       axes(1:1),  lnd%time, '<ltype> litter total carbon leaching', 'kg C/m2/s', missing_value=-100.0 )
  id_litter_total_ON_leaching(:) = register_litter_diag_fields ( module_name, '<ltype>litter_total_ON_leaching', &
       axes(1:1),  lnd%time, '<ltype> litter total organic N leaching', 'kg N/m2/s', missing_value=-100.0 )
  id_litter_NO3_leaching(:) = register_litter_diag_fields ( module_name, '<ltype>litter_total_NO3_leaching', &
       axes(1:1),  lnd%time, '<ltype> litter NO3 leaching', 'kg N/m2/s', missing_value=-100.0 )
  id_litter_NH4_leaching(:) = register_litter_diag_fields ( module_name, '<ltype>litter_total_NH4_leaching', &
       axes(1:1),  lnd%time, '<ltype> litter NH4 leaching', 'kg N/m2/s', missing_value=-100.0 )

  id_litter_N(:,:) = register_litter_soilc_diag_fields ( module_name, '<ctype>_<ltype>litter_N', &
       axes(1:1), lnd%time, '<ctype> <ltype> litter nitrogen', 'kg N/m2', missing_value=-100.0 )
  id_litter_dissolved_N(:,:) = register_litter_soilc_diag_fields ( module_name, '<ctype>_<ltype>litter_dissolved_N', &
       axes(1:1), lnd%time, '<ctype> <ltype> litter dissolved nitrogen', 'kg N/m2', missing_value=-100.0 )
  id_litter_protected_N(:,:) = register_litter_soilc_diag_fields ( module_name, '<ctype>_<ltype>litter_protected_N', &
       axes(1:1), lnd%time, '<ctype> <ltype> litter protected nitrogen', 'kg N/m2', missing_value=-100.0 )
  id_litter_rsoil_N(:,:) = register_litter_soilc_diag_fields ( module_name, 'rsoil_<ltype>litter_<ctype>_N',  &
       axes(1:1), lnd%time, 'surface <ltype> litter <ctype> nitrogen degradation', 'kg N/(m2 year)', missing_value=-100.0 )
!  id_litter_N_leaching(:,:) = register_litter_soilc_diag_fields ( module_name, '<ctype>_<ltype>litter_N_leaching', &
!       axes(1:1), lnd%time, '<ltype> litter <ctype> N leaching','kg/(m2 s)', missing_value=-100.0)
  id_litter_livemic_N(:) = register_litter_diag_fields ( module_name, '<ltype>litter_live_microbe_N', &
       axes(1:1),  lnd%time, 'live microbe <ltype> litter nitrogen', 'kg N/m2', missing_value=-100.0 )
  id_litter_total_N(:) = register_litter_diag_fields ( module_name, '<ltype>litter_total_N', &
       axes(1:1),  lnd%time, '<ltype> litter total nitrogen', 'kg N/m2', missing_value=-100.0 )
 id_litter_ammonium(:) = register_litter_diag_fields ( module_name, '<ltype>litter_ammonium', &
      axes(1:1),  lnd%time, '<ltype> litter ammonium', 'kg N/m2', missing_value=-100.0 )
  id_litter_nitrate(:) = register_litter_diag_fields ( module_name, '<ltype>litter_nitrate', &
       axes(1:1),  lnd%time, '<ltype> litter nitrate', 'kg N/m2', missing_value=-100.0 )


  id_total_NH4 = register_tiled_diag_field ( module_name, 'total_soil_NH4', axes(1:1),  &
      lnd%time, 'total NH4 including litter', 'kg N/m2', missing_value=-100.0 )
  id_total_NO3 = register_tiled_diag_field ( module_name, 'total_soil_NO3', axes(1:1),  &
       lnd%time, 'total NO3 including litter', 'kg N/m2', missing_value=-100.0 )
  id_soil_NO3 = register_tiled_diag_field ( module_name, 'soil_NO3', axes,  &
       lnd%time, 'NO3 per layer', 'kg N/m3', missing_value=-100.0 )
  id_soil_NH4 = register_tiled_diag_field ( module_name, 'soil_NH4', axes,  &
       lnd%time, 'NH4 per layer', 'kg N/m3', missing_value=-100.0 )
  id_total_denitrification_rate = register_tiled_diag_field ( module_name, 'total_denitrification_rate',  &
       (/id_ug/), lnd%time, 'Total denitrification', 'kg N/(m2 year)', &
       missing_value=-100.0 )
 id_soil_denitrification_rate = register_tiled_diag_field ( module_name, 'soil_denitrification_rate', axes,  &
       lnd%time, 'Denitrification rate per layer', 'kg N/m3/year', missing_value=-100.0 )
  id_NO3_div_loss = register_tiled_diag_field ( module_name, 'NO3_div_loss', axes(1:1), &
       lnd%time, 'total rate of NO3 divergence loss', 'kg N/m^2/s', missing_value=initval)
  id_NH4_div_loss = register_tiled_diag_field ( module_name, 'NH4_div_loss', axes(1:1), &
        lnd%time, 'total rate of NH4 divergence loss', 'kg N/m^2/s', missing_value=initval)
  id_total_N_mineralization_rate = register_tiled_diag_field ( module_name, 'total_N_mineralization_rate',  &
       (/id_ug/), lnd%time, 'Total N mineralization', 'kg N/(m2 year)', &
       missing_value=-100.0 )
 id_total_N_immobilization_rate = register_tiled_diag_field ( module_name, 'total_N_immobilization_rate',  &
      (/id_ug/), lnd%time, 'Total N immobilization', 'kg N/(m2 year)', &
      missing_value=-100.0 )
  id_total_nitrification_rate = register_tiled_diag_field ( module_name, 'total_nitrification_rate',  &
       (/id_ug/), lnd%time, 'Total nitrification', 'kg N/(m2 year)', &
       missing_value=-100.0 )

  id_Qmax = register_tiled_diag_field ( module_name, 'soil_Qmax', axes(1:1),  &
       lnd%time, 'Maximum clay sorptive capacity', 'kg C/m3', missing_value=-100.0 )

  id_nsoilcohorts = register_tiled_diag_field ( module_name, 'n_soil_cohorts', axes,  &
       lnd%time, 'number of soil cohorts', missing_value=-100.0 )
  id_deadmic_total_C = register_tiled_diag_field ( module_name, 'deadmic_total_C', axes(1:1),  &
       lnd%time, 'total dead microbe soil carbon', 'kg C/m2', missing_value=-100.0 )
  id_deadmic_total_N = register_tiled_diag_field ( module_name, 'deadmic_total_N', axes(1:1),  &
       lnd%time, 'total dead microbe soil nitrogen', 'kg N/m2', missing_value=-100.0 )
  id_livemic_C = register_tiled_diag_field ( module_name, 'live_microbe_C', axes,  &
       lnd%time, 'Total live microbe soil carbon', 'kg C/m3', missing_value=-100.0 )
  id_livemic_N = register_tiled_diag_field ( module_name, 'live_microbe_N', axes,  &
       lnd%time, 'Total live microbe soil nitrogen', 'kg N/m3', missing_value=-100.0 )
  id_livemic_total_C = register_tiled_diag_field ( module_name, 'livemic_total_C', axes(1:1),  &
       lnd%time, 'total live microbe soil carbon', 'kg C/m2', missing_value=-100.0 )
  id_livemic_total_N = register_tiled_diag_field ( module_name, 'livemic_total_N', axes(1:1),  &
       lnd%time, 'total live microbe soil nitrogen', 'kg N/m2', missing_value=-100.0 )
  id_protected_total_C = register_tiled_diag_field ( module_name, 'protected_total_C', axes(1:1),  &
       lnd%time, 'total protected soil carbon', 'kg C/m2', missing_value=-100.0 )
  id_protected_total_N = register_tiled_diag_field ( module_name, 'protected_total_N', axes(1:1),  &
       lnd%time, 'total protected soil nitrogen', 'kg N/m2', missing_value=-100.0 )
  id_dissolved_total_C = register_tiled_diag_field ( module_name, 'dissolved_total_C', axes(1:1),  &
       lnd%time, 'total dissolved soil carbon', 'kg C/m2', missing_value=-100.0 )
  id_dissolved_total_N = register_tiled_diag_field ( module_name, 'dissolved_total_N', axes(1:1),  &
       lnd%time, 'total dissolved soil nitrogen', 'kg N/m2', missing_value=-100.0 )
  id_total_soil_C = register_tiled_diag_field ( module_name, 'total_soil_C', axes(1:1),  &
       lnd%time, 'total soil carbon', 'kg C/m2', missing_value=-100.0 )
  id_total_soil_N = register_tiled_diag_field ( module_name, 'total_soil_N', axes(1:1),  &
       lnd%time, 'total soil nitrogen', 'kg N/m2', missing_value=-100.0 )
  id_total_C_leaching = register_tiled_diag_field ( module_name, 'total_C_leaching', axes, &
       lnd%time, 'net layer total vertical soil C leaching', 'kg/(m2 s)', missing_value=initval)
  id_total_ON_leaching = register_tiled_diag_field ( module_name, 'total_ON_leaching', axes, &
       lnd%time, 'net layer total vertical soil organic N leaching', 'kg/(m2 s)', missing_value=initval)
  id_NO3_leaching = register_tiled_diag_field ( module_name, 'NO3_leaching', axes, &
       lnd%time, 'net layer vertical soil NO3 leaching', 'kg/(m2 s)', missing_value=initval)
  id_NH4_leaching = register_tiled_diag_field ( module_name, 'NH4_leaching', axes, &
       lnd%time, 'net layer vertical soil NH4 leaching', 'kg/(m2 s)', missing_value=initval)
  ! id_livemic_C_leaching = register_tiled_diag_field ( module_name, 'livemic_C_leaching', axes, &
  !      lnd%time, 'net layer live microbe C leaching',  'kg/(m2 s)', missing_value=-100.0)
  !id_protected_C_leaching = register_tiled_diag_field ( module_name, 'protected_C_leaching', axes, &
  !     lnd%time, 'net layer protected soil C leaching',  'kg/(m2 s)', missing_value=-100.0)
  id_surf_DOC_loss = register_tiled_diag_field ( module_name, 'surf_DOC_loss', axes(1:1), &
       lnd%time, 'loss of top layer DOC to surface runoff due to efflux', 'kg C/m^2/s', &
       missing_value=initval)
  id_fsc = register_tiled_diag_field ( module_name, 'fsc', axes(1:1),  &
       lnd%time, 'total fast soil carbon', 'kg C/m2', missing_value=-100.0 )
  id_fsn = register_tiled_diag_field ( module_name, 'fsn', axes(1:1),  &
       lnd%time, 'total fast soil nitrogen', 'kg N/m2', missing_value=-100.0 )
  id_ssc = register_tiled_diag_field ( module_name, 'ssc', axes(1:1),  &
       lnd%time, 'total slow soil carbon', 'kg C/m2', missing_value=-100.0 )
  id_ssn = register_tiled_diag_field ( module_name, 'ssn', axes(1:1),  &
       lnd%time, 'total slow soil nitrogen', 'kg N/m2', missing_value=-100.0 )
  id_lwc = register_tiled_diag_field ( module_name, 'soil_liq', axes,  &
       lnd%time, 'bulk density of liquid water', 'kg/m3', missing_value=-100.0 )
  id_swc  = register_tiled_diag_field ( module_name, 'soil_ice',  axes,  &
       lnd%time, 'bulk density of solid water', 'kg/m3',  missing_value=-100.0 )
  id_psi = register_tiled_diag_field ( module_name, 'soil_psi', axes,  &
       lnd%time, 'soil-water matric head', 'm', missing_value=-100.0 )
  id_temp  = register_tiled_diag_field ( module_name, 'soil_T',  axes,       &
       lnd%time, 'temperature',            'degK',  missing_value=-100.0 )
  id_ie  = register_tiled_diag_field ( module_name, 'soil_rie',  axes(1:1),  &
       lnd%time, 'inf exc runf',            'kg/(m2 s)',  missing_value=-100.0 )
  id_sn  = register_tiled_diag_field ( module_name, 'soil_rsn',  axes(1:1),  &
       lnd%time, 'satn runf',            'kg/(m2 s)',  missing_value=-100.0 )
  id_bf  = register_tiled_diag_field ( module_name, 'soil_rbf',  axes(1:1),  &
       lnd%time, 'baseflow',            'kg/(m2 s)',  missing_value=-100.0 )
  id_if  = register_tiled_diag_field ( module_name, 'soil_rif',  axes(1:1),  &
       lnd%time, 'interflow',            'kg/(m2 s)',  missing_value=-100.0 )
  id_al  = register_tiled_diag_field ( module_name, 'soil_ral',  axes(1:1),  &
       lnd%time, 'active layer flow',    'kg/(m2 s)',  missing_value=-100.0 )
  id_nu  = register_tiled_diag_field ( module_name, 'soil_rnu',  axes(1:1),  &
       lnd%time, 'numerical runoff',    'kg/(m2 s)',  missing_value=-100.0 )
  id_sc  = register_tiled_diag_field ( module_name, 'soil_rsc',  axes(1:1),  &
       lnd%time, 'lm2 groundwater runoff',    'kg/(m2 s)',  missing_value=-100.0 )
  id_hie  = register_tiled_diag_field ( module_name, 'soil_hie',  axes(1:1), &
       lnd%time, 'heat ie runf',            'W/m2',  missing_value=-100.0 )
  id_hsn  = register_tiled_diag_field ( module_name, 'soil_hsn',  axes(1:1), &
       lnd%time, 'heat sn runf',            'W/m2',  missing_value=-100.0 )
  id_hbf  = register_tiled_diag_field ( module_name, 'soil_hbf',  axes(1:1), &
       lnd%time, 'heat bf runf',            'W/m2',  missing_value=-100.0 )
  id_hif  = register_tiled_diag_field ( module_name, 'soil_hif',  axes(1:1), &
       lnd%time, 'heat if runf',            'W/m2',  missing_value=-100.0 )
  id_hal  = register_tiled_diag_field ( module_name, 'soil_hal',  axes(1:1), &
       lnd%time, 'heat al runf',            'W/m2',  missing_value=-100.0 )
  id_hnu  = register_tiled_diag_field ( module_name, 'soil_hnu',  axes(1:1), &
       lnd%time, 'heat nu runoff',          'W/m2',  missing_value=-100.0 )
  id_hsc  = register_tiled_diag_field ( module_name, 'soil_hsc',  axes(1:1), &
       lnd%time, 'heat sc runoff',          'W/m2',  missing_value=-100.0 )
  id_evap  = register_tiled_diag_field ( module_name, 'soil_evap',  axes(1:1), &
       lnd%time, 'soil evap',            'kg/(m2 s)',  missing_value=-100.0 )
  id_excess  = register_tiled_diag_field ( module_name, 'sfc_excess',  axes(1:1),  &
       lnd%time, 'sfc excess pushed down',    'kg/(m2 s)',  missing_value=-100.0 )

  id_uptk_n_iter  = register_tiled_diag_field ( module_name, 'uptake_n_iter',  axes(1:1), &
       lnd%time, 'number of iterations for soil uptake',  missing_value=-100.0 )
  id_uptk = register_tiled_diag_field ( module_name, 'soil_uptk', axes, &
       lnd%time, 'uptake of water by roots', 'kg/(m2 s)',  missing_value=-100.0 )
  id_psi_x0 = register_tiled_diag_field ( module_name, 'soil_psix0', axes(1:1), &
       lnd%time, 'xylem potential at z=0', 'm',  missing_value=-100.0 )
  id_deficit = register_tiled_diag_field ( module_name, 'soil_def', axes(1:1), &
       lnd%time, 'groundwater storage deficit', '-',  missing_value=-100.0 )
  id_deficit_2 = register_tiled_diag_field ( module_name, 'soil_def2', axes(1:1), &
       lnd%time, 'groundwater storage deficit2', '-',  missing_value=-100.0 )
  id_deficit_3 = register_tiled_diag_field ( module_name, 'soil_def3', axes(1:1), &
       lnd%time, 'groundwater storage deficit3', '-',  missing_value=-100.0 )
  id_deficit_4 = register_tiled_diag_field ( module_name, 'soil_def4', axes(1:1), &
       lnd%time, 'groundwater storage deficit4', '-',  missing_value=-100.0 )
  id_psi_bot = register_tiled_diag_field ( module_name, 'soil_psi_n', axes(1:1), &
       lnd%time, 'psi at bottom of soil column', 'm',  missing_value=-100.0 )
  id_sat_frac = register_tiled_diag_field ( module_name, 'soil_fsat', axes(1:1), &
       lnd%time, 'fraction of soil area saturated at surface', '-',  missing_value=-100.0 )
  id_sat_depth = register_tiled_diag_field ( module_name, 'soil_wtdep', axes(1:1), &
       lnd%time, 'depth below sfc to saturated soil', 'm',  missing_value=-100.0 )
  id_sat_dept2 = register_tiled_diag_field ( module_name, 'soil_wtdp2', axes(1:1), &
       lnd%time, 'alt depth below sfc to saturated soil', 'm',  missing_value=-100.0 )
  id_z_cap = register_tiled_diag_field ( module_name, 'soil_zcap', axes(1:1), &
       lnd%time, 'depth below sfc to capillary fringe', 'm',  missing_value=-100.0 )

  id_div_bf = register_tiled_diag_field ( module_name, 'soil_dvbf', axes, &
       lnd%time, 'baseflow by layer', 'kg/(m2 s)',  missing_value=-100.0 )
  id_div_if = register_tiled_diag_field ( module_name, 'soil_dvif', axes, &
       lnd%time, 'interflow by layer', 'kg/(m2 s)',  missing_value=-100.0 )
  id_div_al = register_tiled_diag_field ( module_name, 'soil_dval', axes, &
       lnd%time, 'active-layer flow by layer', 'kg/(m2 s)',  missing_value=-100.0 )

  id_cf_1 = register_tiled_diag_field ( module_name, 'soil_cf_1', axes(1:1), &
       lnd%time, 'soil_cf_1', 'm',  missing_value=-100.0 )
  id_cf_3 = register_tiled_diag_field ( module_name, 'soil_cf_3', axes(1:1), &
       lnd%time, 'soil_cf_3', 'm',  missing_value=-100.0 )
  id_wt_1 = register_tiled_diag_field ( module_name, 'soil_wt_1', axes(1:1), &
       lnd%time, 'soil_wt_1', 'm',  missing_value=-100.0 )
  id_wt_2 = register_tiled_diag_field ( module_name, 'soil_wt_2', axes(1:1), &
       lnd%time, 'soil_wt_2', 'm',  missing_value=-100.0 )
  id_wt_2a = register_tiled_diag_field ( module_name, 'soil_wt_2a', axes(1:1), &
       lnd%time, 'Water Table Depth from Surface to Saturation', 'm',  missing_value=-100.0 )
  id_wt_2b = register_tiled_diag_field ( module_name, 'soil_wt_2b', axes(1:1), &
       lnd%time, 'Water Table Depth from Surface to Liquid Saturation', 'm',  missing_value=-100.0 )
  id_wt_3 = register_tiled_diag_field ( module_name, 'soil_wt_3', axes(1:1), &
       lnd%time, 'soil_wt_3', 'm',  missing_value=-100.0 )
  id_wt2_3 = register_tiled_diag_field ( module_name, 'soil_wt2_3', axes(1:1), &
       lnd%time, 'soil_wt2_3', 'm',  missing_value=-100.0 )
  id_wt_4 = register_tiled_diag_field ( module_name, 'soil_wt_4', axes(1:1), &
       lnd%time, 'Interpolated psi = 0 from Bottom Up', 'm',  missing_value=-100.0 )

  id_active_layer = register_tiled_diag_field ( module_name, 'soil_alt', axes(1:1), &
       lnd%time, 'active-layer thickness', 'm',  missing_value=-100.0 )
  id_heat_cap = register_tiled_diag_field ( module_name, 'soil_heat_cap',  &
       axes, lnd%time, 'heat capacity of dry soil','J/(m3 K)', missing_value=-100.0 )
  id_thermal_cond =  register_tiled_diag_field ( module_name, 'soil_tcon', &
       axes, lnd%time, 'soil thermal conductivity', 'W/(m K)',  missing_value=-100.0 )

  id_surface_water = register_tiled_diag_field (module_name, 'surface_water', &
       axes(1:1), lnd%time, 'surface water storage', 'm', missing_value=-100.0 )
  id_inun_frac = register_tiled_diag_field (module_name, 'inun_fraction', &
       axes(1:1), lnd%time, 'inundated area fraction', '-', missing_value=-100.0 )
  if (gw_option == GW_TILED) then
     id_wet_frac = register_tiled_diag_field (module_name, 'wet_fraction', &
       axes(1:1), lnd%time, 'diagnostic wetland fraction', '-', missing_value=-100.0 )
  end if
  if (gw_option == GW_TILED .and. simple_inundation) then
      id_rsn_frac = register_tiled_diag_field (module_name, 'surface_runoff_frac', &
         axes(1:1), lnd%time, 'effective fraction of throughfall converted to sat-excess surface runoff', '-', missing_value=-100.0 )
  end if
  id_flow = register_tiled_diag_field (module_name, 'flow', axes, &
       lnd%time, 'vertical soil water flow at interface above (+ downward)', 'mm/s', missing_value=initval )
  id_reflux = register_tiled_diag_field (module_name, 'reflux', axes(1:1), &
       lnd%time, 'upwards flow of soil water at surface; zero if flow into surface', 'mm/s', missing_value=-100.0 )
  id_macro_infilt = register_tiled_diag_field (module_name, 'macro_inf', axes(1:1), &
       lnd%time, 'infiltration (decrease to IE runoff) at soil surface due to vertical macroporosity', 'mm/s', missing_value=-100.0 )

  id_type = register_tiled_static_field ( module_name, 'soil_type',  &
       axes(1:1), 'soil type', missing_value=-1.0 )
  id_tau_gw = register_tiled_static_field ( module_name, 'tau_gw',  &
       axes(1:1), 'groundwater residence time', 's', missing_value=-100.0 )
  id_slope_l = register_tiled_static_field ( module_name, 'slope_l',  &
       axes(1:1), 'hillslope length', 'm', missing_value=-100.0 )
  id_slope_Z = register_tiled_static_field ( module_name, 'soil_rlief',  &
       axes(1:1), 'hillslope relief', 'm', missing_value=-100.0 )
  id_zeta_bar = register_tiled_static_field ( module_name, 'zeta_bar',  &
       axes(1:1), 'hillslope zeta bar', '-', missing_value=-100.0 )
  id_e_depth = register_tiled_static_field ( module_name, 'soil_depth',  &
       axes(1:1), 'soil hydraulic e-folding depth', 'm', missing_value=-100.0 )
  id_zeta = register_tiled_static_field ( module_name, 'soil_zeta',      &
       axes(1:1), 'soil depth/topo relief', '-',  missing_value=-100.0 )
  id_tau = register_tiled_static_field ( module_name, 'soil_tau',        &
       axes(1:1), 'gw transmissivity/soil transmissivity', '-',  missing_value=-100.0 )
  id_vwc_wilt = register_tiled_static_field ( module_name, 'soil_wilt',  &
       axes(1:1), 'wilting water content', '-', missing_value=-100.0 )
  id_vwc_fc = register_tiled_static_field ( module_name, 'soil_fc',  &
       axes(1:1), 'field capacity', '-', missing_value=-100.0 )
  id_vwc_sat = register_tiled_static_field ( module_name, 'soil_sat',  &
       axes(1:1), 'soil porosity', '-', missing_value=-100.0 )
  id_K_sat = register_tiled_static_field ( module_name, 'soil_Ksat',  &
       axes(1:1), 'soil sat. hydraulic conductivity', 'kg /(m2 s)', missing_value=-100.0 )
  id_K_gw  = register_tiled_static_field ( module_name, 'soil_K_gw',  &
       axes(1:1), 'deep hydraulic conductivity', 'kg /(m2 s)', missing_value=-100.0 )
  id_w_fc = register_tiled_static_field ( module_name, 'w_fc',  &
       axes, 'soil field capacity', missing_value=-1.0 )
  id_alpha = register_tiled_static_field ( module_name, 'soil_alpha',  &
       axes, 'soil microscopic length scale', missing_value=-1.0 )
  id_refl_dry_dir = register_tiled_static_field ( module_name, 'refl_dry_dir',  &
       (/id_ug, id_band/), 'reflectance of dry soil for direct light', &
       missing_value=-1.0 )
  id_refl_dry_dif = register_tiled_static_field ( module_name, 'refl_dry_dif',  &
       (/id_ug, id_band/), 'reflectance of dry soil for diffuse light', &
       missing_value=-1.0 )
  id_refl_sat_dir = register_tiled_static_field ( module_name, 'refl_sat_dir',  &
       (/id_ug, id_band/), 'reflectance of saturated soil for direct light', &
       missing_value=-1.0 )
  id_refl_sat_dif = register_tiled_static_field ( module_name, 'refl_sat_dif',  &
       (/id_ug, id_band/), 'reflectance of saturated soil for diffuse light', &
       missing_value=-1.0 )
  id_f_iso_dry = register_tiled_static_field ( module_name, 'f_iso_dry',  &
       (/id_ug, id_band/), 'isotropic brdf weight, dry soil', &
       missing_value=-1.0 )
  id_f_vol_dry = register_tiled_static_field ( module_name, 'f_vol_dry',  &
       (/id_ug, id_band/), 'volumetric brdf weight, dry soil', &
       missing_value=-1.0 )
  id_f_geo_dry = register_tiled_static_field ( module_name, 'f_geo_dry',  &
       (/id_ug, id_band/), 'geometric brdf weight, dry soil', &
       missing_value=-1.0 )
  id_f_iso_sat = register_tiled_static_field ( module_name, 'f_iso_sat',  &
       (/id_ug, id_band/), 'isotropic brdf weight, saturated soil', &
       missing_value=-1.0 )
  id_f_vol_sat = register_tiled_static_field ( module_name, 'f_vol_sat',  &
       (/id_ug, id_band/), 'volumetric brdf weight, saturated soil', &
       missing_value=-1.0 )
  id_f_geo_sat = register_tiled_static_field ( module_name, 'f_geo_sat',  &
       (/id_ug, id_band/), 'geometric brdf weight, saturated soil', &
       missing_value=-1.0 )

  id_asoil = register_tiled_diag_field ( module_name, 'asoil', &
       (/id_ug/), lnd%time, 'aerobic activity modifier', &
       missing_value=-100.0 )

  ! the following fields are for compatibility with older diag tables only
  call add_tiled_static_field_alias ( id_slope_Z, module_name, 'slope_Z',  &
       axes(1:1), 'hillslope relief (obsolete, use "soil_rlief" instead)',&
       'm', missing_value=-100.0 )
  call add_tiled_static_field_alias ( id_e_depth, module_name, 'e_depth',  &
       axes(1:1), 'soil e-folding depth (obsolete, use "soil_depth" instead)', &
       'm', missing_value=-100.0 )
  call add_tiled_static_field_alias ( id_vwc_wilt, module_name, 'vwc_wilt',  &
       axes(1:1), 'wilting water content (obsolete, use "soil_wilt" instead)', &
       '-', missing_value=-100.0 )
  call add_tiled_static_field_alias ( id_vwc_fc, module_name, 'vwc_fc',  &
       axes(1:1), 'field capacity (obsolete, use "soil_fc" instead)', &
       '-', missing_value=-100.0 )
  call add_tiled_static_field_alias ( id_vwc_sat, module_name, 'vwc_sat',  &
       axes(1:1), 'soil porosity (obsolete, use "soil_sat")', &
       '-', missing_value=-100.0 )
  call add_tiled_static_field_alias ( id_K_sat, module_name, 'K_sat',  &
       axes(1:1), 'soil sat. hydraulic conductivity (obsolte, use "soil_Ksat" instead)', &
       'kg /(m2 s)', missing_value=-100.0 )

#ifdef ZMSDEBUG_TRIDIAGTEST
  ! For testing tridiagonal solution for advection
  id_st_diff = register_tiled_diag_field ( module_name, 'soil_T_diff', axes,  &
      lnd%time, 'soil Temperature difference after advection with tridiagonal solution', &
      'K', missing_value=-100.0 )
#endif

  ! CMOR variables
  ! tsl (soil temperature) should be reported as missing for the non-soil grid cells;
  ! also averaging over non-soil tiles does not make sense and therefore is not done.
  ! The only difference with soil_T is metadata (units and standard_name).
  call add_tiled_diag_field_alias ( id_temp, cmor_name, 'tsl', axes(1:2),  &
       lnd%time, 'Temperature of Soil', 'K', missing_value=-100.0, &
       standard_name='soil_temperature', fill_missing=.FALSE.)
  ! set up weights for mrsos averaging
  allocate(mrsos_weight(num_l), mrs1m_weight(num_l))
  do l = 1,num_l
     mrsos_weight(l) = min(1.0,max(0.0,(cmor_mrsos_depth-zhalf(l))/dz(l)))
     mrs1m_weight(l) = min(1.0,max(0.0,(1.0             -zhalf(l))/dz(l)))
  enddo
  ! set the default sub-sampling filter for the fields below
  call set_default_diag_filter('land')
  id_mrsofc = register_tiled_static_field ( cmor_name, 'mrsofc', axes(1:1),  &
       'Capacity of Soil to Store Water', 'kg m-2', missing_value=-100.0, &
       standard_name='soil_moisture_content_at_field_capacity', fill_missing=.TRUE.)
  id_mrlsl = register_tiled_diag_field ( cmor_name, 'mrlsl', axes,  &
       lnd%time, 'Water Content of Soil Layer', 'kg m-2', missing_value=-100.0, &
       standard_name='moisture_content_of_soil_layer', fill_missing=.TRUE.)
  id_mrsfl = register_tiled_diag_field ( cmor_name, 'mrsfl', axes,  &
       lnd%time, 'Frozen Water Content of Soil Layer', 'kg m-2', missing_value=-100.0, &
       standard_name='frozen_moisture_content_of_soil_layer', fill_missing=.TRUE.)
  id_mrsll = register_tiled_diag_field ( cmor_name, 'mrsll', axes,  &
       lnd%time, 'Liquid Water Content of Soil Layer', 'kg m-2', missing_value=-100.0, &
       standard_name='liquid_moisture_content_of_soil_layer', fill_missing=.TRUE.)
  id_mrsol = register_tiled_diag_field ( cmor_name, 'mrsol', axes,  &
       lnd%time, 'Total Water Content of Soil Layer', 'kg m-2', missing_value=-100.0, &
       standard_name='moisture_content_of_soil_layer', fill_missing=.TRUE.)
  id_mrso  = register_tiled_diag_field ( cmor_name, 'mrso', axes(1:1),  &
       lnd%time, 'Total Soil Moisture Content', 'kg m-2', missing_value=-100.0, &
       standard_name='soil_moisture_content', fill_missing=.TRUE.)
  call add_tiled_diag_field_alias ( id_mrso, cmor_name, 'mrsoLut', axes(1:1),  &
       lnd%time, 'Total Soil Moisture Content', 'kg m-2', missing_value=-100.0, &
       standard_name='soil_moisture_content', fill_missing=.FALSE.)
  id_mrsos  = register_tiled_diag_field ( cmor_name, 'mrsos', axes(1:1),  &
       lnd%time, 'Moisture in Upper Portion of Soil Column', &
       'kg m-2', missing_value=-100.0, standard_name='moisture_content_of_soil_layer', &
       fill_missing=.TRUE.)
  call  add_tiled_diag_field_alias ( id_mrsos, cmor_name, 'mrsosLut', axes(1:1),  &
       lnd%time, 'Moisture in Upper Portion of Soil Column of Land Use Tile', &
       'kg m-2', missing_value=-100.0, standard_name='moisture_content_of_soil_layer', &
       fill_missing=.FALSE.)
  id_mrfso = register_tiled_diag_field ( cmor_name, 'mrfso', axes(1:1),  &
       lnd%time, 'Soil Frozen Water Content', 'kg m-2', missing_value=-100.0, &
       standard_name='soil_frozen_water_content', fill_missing=.TRUE.)
  id_mrlso = register_tiled_diag_field ( cmor_name, 'mrlso', axes(1:1),  &
       lnd%time, 'Soil Liquid Water Content', 'kg m-2', missing_value=-100.0, &
       standard_name='soil_liquid_water_content', fill_missing=.TRUE.)
  id_mrros = register_tiled_diag_field ( cmor_name, 'mrros',  axes(1:1),  &
       lnd%time, 'Surface Runoff', 'kg m-2 s-1',  missing_value=-100.0, &
       standard_name='surface_runoff_flux', fill_missing=.TRUE.)
  id_mrro = register_tiled_diag_field ( cmor_name, 'mrro',  axes(1:1),  &
       lnd%time, 'Total Runoff', 'kg m-2 s-1',  missing_value=-100.0, &
       standard_name='runoff_flux', fill_missing=.TRUE.)
  call add_tiled_diag_field_alias ( id_mrro, cmor_name, 'mrroLut',  axes(1:1),  &
       lnd%time, 'Total Runoff From Land Use Tile', 'kg m-2 s-1',  missing_value=-100.0, &
       standard_name='runoff_flux', fill_missing=.FALSE.)

  id_csoil = register_tiled_diag_field ( cmor_name, 'cSoil', axes(1:1),  &
       lnd%time, 'Carbon in Soil Pool', 'kg m-2', missing_value=-100.0, &
       standard_name='soil_carbon_content', fill_missing=.TRUE.)
  call add_tiled_diag_field_alias ( id_csoil, cmor_name, 'cSoilLut', axes(1:1),  &
       lnd%time, 'Carbon  In Soil Pool On Land Use Tiles', 'kg m-2', missing_value=-100.0, &
       standard_name='soil_carbon_content', fill_missing=.FALSE.)

  id_csoilfast = register_tiled_diag_field ( cmor_name, 'cSoilFast', axes(1:1),  &
       lnd%time, 'Carbon Mass in Fast Soil Pool', 'kg m-2', missing_value=-100.0, &
       standard_name='fast_soil_pool_carbon_content', fill_missing=.TRUE.)
  id_csoilmedium = register_tiled_diag_field ( cmor_name, 'cSoilMedium', axes(1:1),  &
       lnd%time, 'Carbon Mass in Medium Soil Pool', 'kg m-2', missing_value=-100.0, &
       standard_name='medium_soil_pool_carbon_content', fill_missing=.TRUE.)
  id_csoilslow = register_tiled_diag_field ( cmor_name, 'cSoilSlow', axes(1:1),  &
       lnd%time, 'Carbon Mass in Slow Soil Pool', 'kg m-2', missing_value=-100.0, &
       standard_name='slow_soil_pool_carbon_content', fill_missing=.TRUE.)
  id_cLitter = register_tiled_diag_field ( cmor_name, 'cLitter', axes(1:1), &
       lnd%time, 'Carbon Mass in Litter Pool', 'kg m-2', &
       missing_value=-100.0, standard_name='litter_carbon_content', &
       fill_missing=.TRUE.)
  call add_tiled_diag_field_alias ( id_cLitter, cmor_name, 'cLitterLut', axes(1:1),  &
       lnd%time, 'carbon in above and belowground litter pools on land use tiles', &
       'kg m-2', missing_value=-100.0, &
       standard_name='litter_carbon_content', fill_missing=.FALSE.)
  id_cLitterCwd = register_tiled_diag_field ( cmor_name, 'cLitterCwd', axes(1:1), &
       lnd%time, 'Carbon Mass in Coarse Woody Debris', 'kg m-2', &
       missing_value=-100.0, standard_name='litter_wood_debris_carbon_content', &
       fill_missing=.TRUE.)
  id_rh = register_tiled_diag_field ( cmor_name, 'rh', (/id_ug/), &
       lnd%time, 'Heterotrophic Respiration', 'kg m-2 s-1', missing_value=-1.0, &
       standard_name='heterotrophic_respiration_carbon_flux', fill_missing=.TRUE.)
  call add_tiled_diag_field_alias ( id_rh, cmor_name, 'rhLut', axes(1:1),  &
       lnd%time, 'Soil Heterotrophic Respiration On Land Use Tile', 'kg m-2 s-1', &
       standard_name='heterotrophic_respiration_carbon_flux', fill_missing=.FALSE., &
       missing_value=-100.0)
  id_mrs1mLut = register_tiled_diag_field ( cmor_name, 'mrs1mLut', axes(1:1), &
       lnd%time, 'Moisture in Top 1 Meter of Land Use Tile Soil Column', 'kg m-2', &
       missing_value=-100.0, standard_name='moisture_content_of_soil_layer', &
       fill_missing=.FALSE.)

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
  character(267) :: filename
  type(land_restart_type) :: restart ! restart file i/o object
  type(land_tile_enum_type)     :: ce   ! tile list enumerator
  type(land_tile_type), pointer :: tile   ! pointer to current tile
  integer :: i,k

  call error_mesg('soil_end','writing NetCDF restart',NOTE)
! Note that filename is updated for tile & rank numbers during file creation
  filename = trim(timestamp)//'soil.res.nc'
  call init_land_restart(restart, filename, soil_tile_exists, tile_dim_length)
  call add_restart_axis(restart,'zfull',zfull(1:num_l),'Z','m','full level',sense=-1)
  if (soil_carbon_option==SOILC_CORPSE.or.soil_carbon_option==SOILC_CORPSE_N) then
     call add_restart_axis(restart,'soilCCohort',(/(float(i),i=1,soilMaxCohorts)/),'CC')
     call add_restart_axis(restart,'litterCCohort',(/1.0/),'CC')
  endif

  ! write out fields
  call add_tile_data(restart,'temp'         , 'zfull', soil_T_ptr,  'soil temperature','degrees_K')
  call add_tile_data(restart,'wl'           , 'zfull', soil_wl_ptr, 'liquid water content','kg/m2')
  call add_tile_data(restart,'ws'           , 'zfull', soil_ws_ptr, 'solid water content','kg/m2')
  call add_tile_data(restart,'groundwater'  , 'zfull', soil_groundwater_ptr, units='kg/m2' )
  call add_tile_data(restart,'groundwater_T', 'zfull', soil_groundwater_T_ptr, 'groundwater temperature','degrees_K' )
  call add_tile_data(restart,'uptake_T', soil_uptake_T_ptr, 'temperature of transpiring water', 'degrees_K')
  select case(soil_carbon_option)
  case (SOILC_CENTURY, SOILC_CENTURY_BY_LAYER)
     call add_tile_data(restart,'fsc', 'zfull', soil_fast_soil_C_ptr ,'fast soil carbon', 'kg C/m2')
     call add_tile_data(restart,'ssc', 'zfull', soil_slow_soil_C_ptr ,'slow soil carbon', 'kg C/m2')
  case (SOILC_CORPSE, SOILC_CORPSE_N)
     ! make sure all arrays of carbon cohorts are of the same length
     ce = first_elmt(land_tile_map)
     do while (loop_over_tiles(ce,tile))
         if (.not.associated(tile%soil)) cycle
         do i = 1,N_LITTER_POOLS
            call adjust_pool_ncohorts(tile%soil%litter(i))
         enddo
         do i = 1,num_l
            call adjust_pool_ncohorts(tile%soil%soil_organic_matter(i))
         enddo
     end do
     do i = 1, N_C_TYPES
        call add_tile_data(restart,trim(c_shortname(i))//'_soil_C','zfull','soilCCohort',sc_soil_C_ptr,i,trim(c_longname(i))//' soil carbon','kg/m2')
        call add_tile_data(restart,trim(c_shortname(i))//'ProtectedC','zfull','soilCCohort',sc_protected_C_ptr,i,'Protected '//trim(c_longname(i))//' carbon','kg/m2')
        call add_tile_data(restart,'soil_DOC_'//trim(c_shortname(i)),'zfull',sc_DOC_ptr,i,'Dissolved '//trim(c_longname(i))//' carbon','kg/m2')
        do k = 1,N_LITTER_POOLS
           call add_tile_data(restart,trim(l_shortname(k))//'_litter_'//trim(c_shortname(i))//'_C','litterCCohort',sc_litter_litterC_ptr,i,k,trim(l_longname(k))//' litter '//trim(c_longname(i))//' C','kg/m2')
           call add_tile_data(restart,trim(l_shortname(k))//'_litter_'//trim(c_shortname(i))//'ProtectedC','litterCCohort',sc_litter_protectedC_ptr,i,k,trim(l_longname(k))//' litter '//trim(c_longname(i))//' protected C','kg/m2')
           call add_tile_data(restart,trim(l_shortname(k))//'_litter_DOC_'//trim(c_shortname(i)),sc_litter_dissolved_carbon_ptr,i,k,'Dissolved '//trim(l_longname(k))//' litter '//trim(c_longname(i))//' carbon','kg/m2')
        enddo
     enddo
     call add_tile_data(restart,'liveMic' ,'zfull','soilCCohort',sc_livingMicrobeC_ptr,'Living microbial carbon','kg/m2')
     call add_tile_data(restart,'CO2', 'zfull','soilCCohort',sc_CO2_ptr,'Cohort CO2 generated','kg/m2')
     call add_tile_data(restart,'Rtot','zfull','soilCCohort',sc_Rtot_ptr,'Total degradation','kg/m2')
     call add_tile_data(restart,'originalCohortC','zfull','soilCCohort',sc_originalLitterC_ptr,'Cohort original carbon','g/m2')

     do k = 1,N_LITTER_POOLS
        call add_tile_data(restart,trim(l_shortname(k))//'_litter_liveMic_C','litterCCohort',sc_litter_livingMicrobeC_ptr,k,trim(l_longname(k))//' litter live microbe C','kg/m2')
        call add_tile_data(restart,trim(l_shortname(k))//'_litter_CO2','litterCCohort',sc_litter_CO2_ptr,k,trim(l_longname(k))//' litter CO2 generated','kg/m2')
        call add_tile_data(restart,trim(l_shortname(k))//'_litter_Rtot','litterCCohort',sc_litter_Rtot_ptr,k,trim(l_longname(k))//' litter total degradation','kg/m2')
        call add_tile_data(restart,trim(l_shortname(k))//'_litter_originalCohortC','litterCCohort',sc_litter_originalLitterC_ptr,k,trim(l_longname(k))//' litter cohort original carbon','kg/m2')
     enddo

     call add_int_tile_data(restart,'is_peat','zfull',soil_is_peat_ptr,'Is layer peat?','Boolean')

     call add_tile_data(restart,'fast_DOC_leached',     soil_fast_DOC_leached_ptr, 'Cumulative fast DOC leached out of the column', 'kg/m2')
     call add_tile_data(restart,'slow_DOC_leached',     soil_slow_DOC_leached_ptr, 'Cumulative slow DOC leached out of the column', 'kg/m2')
     call add_tile_data(restart,'deadmic_DOC_leached',  soil_deadmic_DOC_leached_ptr, 'Cumulative dead microbe DOC leached out of the column', 'kg/m2')
     if (soil_carbon_option == SOILC_CORPSE_N) then
        do i = 1, N_C_TYPES
           call add_tile_data(restart,trim(c_shortname(i))//'_soil_N', 'zfull','soilCCohort', sc_soil_N_ptr,i,trim(c_longname(i))//' soil nitrogen','kg/m2')
           call add_tile_data(restart,trim(c_shortname(i))//'ProtectedN', 'zfull','soilCCohort', sc_protected_N_ptr,i,'Protected '//trim(c_longname(i))//' soil nitrogen','kg/m2')
           call add_tile_data(restart,'soil_DON_'//trim(c_shortname(i)), 'zfull', sc_DON_ptr,i,'Dissolved '//trim(c_longname(i))//' nitrogen','kg/m2')

           do k = 1,N_LITTER_POOLS
              call add_tile_data(restart,trim(l_shortname(k))//'_litter_'//trim(c_shortname(i))//'_N','litterCCohort',sc_litter_litterN_ptr,i,k,trim(l_longname(k))//' litter '//trim(c_longname(i))//' N','kg/m2')
              call add_tile_data(restart,trim(l_shortname(k))//'_litter_'//trim(c_shortname(i))//'ProtectedN','litterCCohort',sc_litter_protectedN_ptr,i,k,trim(l_longname(k))//' litter '//trim(c_longname(i))//' protected N','kg/m2')
              call add_tile_data(restart,trim(l_shortname(k))//'_litter_DON_'//trim(c_shortname(i)),sc_litter_dissolved_nitrogen_ptr,i,k,'Dissolved '//trim(l_longname(k))//' litter '//trim(c_longname(i))//' nitrogen','kg/m2')
           enddo
        enddo
        call add_tile_data(restart,'liveMicN', 'zfull','soilCCohort', sc_livingMicrobeN_ptr,'Living microbial nitrogen','kg/m2')
        ! FIXME slm: why "original" carbon and nitrogen are is in g/m2
        call add_tile_data(restart,'originalCohortN', 'zfull','soilCCohort', sc_originalLitterN_ptr,'Cohort original nitrogen','g/m2')
        call add_tile_data(restart,'soil_NO3', 'zfull', sc_nitrate_ptr,'Soil nitrate content','kg/m2')
        call add_tile_data(restart,'soil_NH4', 'zfull', sc_ammonium_ptr,'Soil ammonium content','kg/m2')
        call add_tile_data(restart,'soil_nitrif', 'zfull', sc_nitrif_ptr,'Soil cumulative nitrification','kg/m2')
        call add_tile_data(restart,'soil_denitrif', 'zfull', sc_denitrif_ptr,'Soil cumulative denitrification','kg/m2')

        do k = 1,N_LITTER_POOLS
           call add_tile_data(restart,trim(l_shortname(k))//'_litter_liveMic_N', 'litterCCohort', sc_litter_livingMicrobeN_ptr, k, trim(l_longname(k))//' litter live microbe N','kg/m2')
           call add_tile_data(restart,trim(l_shortname(k))//'_litter_originalCohortN', 'litterCCohort', sc_litter_originalLitterN_ptr, k, trim(l_longname(k))//' litter cohort original N','kg/m2')
           call add_tile_data(restart,trim(l_shortname(k))//'_litter_NO3', sc_litter_nitrate_ptr, k, trim(l_longname(k))//' litter nitrate content','kg/m2')
           call add_tile_data(restart,trim(l_shortname(k))//'_litter_NH4', sc_litter_ammonium_ptr, k, trim(l_longname(k))//' litter ammonium content','kg/m2')
           call add_tile_data(restart,trim(l_shortname(k))//'_litter_nitrif', sc_litter_nitrif_ptr, k, trim(l_longname(k))//' litter cumulative nitrification','kg/m2')
           call add_tile_data(restart,trim(l_shortname(k))//'_litter_denitrif', sc_litter_denitrif_ptr, k, trim(l_longname(k))//' litter cumulative denitrification','kg/m2')
        enddo

        call add_tile_data(restart,'gross_nitrogen_flux_into_tile',soil_gross_nitrogen_flux_into_tile_ptr,'Cumulative nitrogen flux into tile','kg/m2')
        call add_tile_data(restart,'gross_nitrogen_flux_out_of_tile',soil_gross_nitrogen_flux_out_of_tile_ptr,'Cumulative nitrogen flux out of tile','kg/m2')

     endif
  case default
     call error_mesg('save_soil_restart','unrecognized soil carbon option -- this should never happen', FATAL)
  end select
  do i = 1, N_C_TYPES
     call add_tile_data(restart,'negative_litter_C_'//trim(c_shortname(i)),sc_negative_litter_C_ptr,i,'accumulated negative '//trim(c_longname(i))//' C litter input','kg/m2')
     call add_tile_data(restart,'negative_litter_N_'//trim(c_shortname(i)),sc_negative_litter_N_ptr,i,'accumulated negative '//trim(c_longname(i))//' N litter input','kg/m2')
  enddo

  call save_land_restart(restart)
  call free_land_restart(restart)

  if (write_soil_carbon_restart) then
     filename = trim(timestamp)//'soil_carbon.res.nc'
     call init_land_restart(restart, filename, soil_tile_exists, tile_dim_length)
     call add_restart_axis(restart,'zfull',zfull(1:num_l),'Z','m','full level',sense=-1)

     select case (soil_carbon_option)
     case (SOILC_CENTURY, SOILC_CENTURY_BY_LAYER)
        call add_tile_data(restart,'asoil_in','zfull',soil_asoil_in_ptr,'aerobic activity modifier', 'unitless')
        call add_tile_data(restart,'fsc_in','zfull',soil_fsc_in_ptr,'fast soil carbon input', 'kg C/m2')
        call add_tile_data(restart,'ssc_in','zfull',soil_ssc_in_ptr,'slow soil carbon input', 'kg C/m2')

     case (SOILC_CORPSE, SOILC_CORPSE_N)
        call add_tile_data(restart,'asoil_in','zfull',soil_asoil_in_ptr,'aerobic activity modifier', 'unitless')

        do i = 1,N_C_TYPES
           ! C inputs
           call add_tile_data(restart,trim(c_shortname(i))//'_protected_C_in','zfull',&
                    sc_protected_C_in_ptr, i,'protected '//trim(c_longname(i))//' soil carbon input', 'kg C/m2')
           call add_tile_data(restart,trim(c_shortname(i))//'_C_in','zfull',&
                    sc_C_in_ptr, i, trim(c_longname(i))//' soil carbon input', 'kg C/m2')
           do k = 1,N_LITTER_POOLS
              call add_tile_data(restart,trim(l_shortname(k))//'litter_'//trim(c_shortname(i))//'_C_in', &
                    sc_litter_C_in_ptr, i, k, trim(c_longname(i))//' '//trim(l_longname(k))//' litter carbon input','kg C/m2')
              call add_tile_data(restart,trim(l_shortname(k))//'litter_'//trim(c_shortname(i))//'_N_in', &
                    sc_litter_N_in_ptr, i, k, trim(c_longname(i))//' '//trim(l_longname(k))//' litter nitrogen input','kg C/m2')
              call add_tile_data(restart,trim(l_shortname(k))//'litter_'//trim(c_shortname(i))//'_C_turnover_accumulated', &
                    sc_litter_C_turnover_ptr, i, k, trim(c_longname(i))//' '//trim(l_longname(k))//' litter carbon turnover', 'year-1')
              call add_tile_data(restart,trim(l_shortname(k))//'litter_'//trim(c_shortname(i))//'_N_turnover_accumulated', &
                    sc_litter_C_turnover_ptr, i, k, trim(c_longname(i))//' '//trim(l_longname(k))//' litter nitrogen turnover', 'year-1')
           enddo
           ! N_inputs
           call add_tile_data(restart,trim(c_shortname(i))//'_protected_N_in','zfull',&
                    sc_protected_N_in_ptr, i,'protected '//trim(c_longname(i))//' soil nitrogen input', 'kg C/m2')
           call add_tile_data(restart,trim(c_shortname(i))//'_N_in','zfull',&
                    sc_N_in_ptr, i, trim(c_longname(i))//' soil nitrogen input', 'kg C/m2')
           ! C turnover rates
           call add_tile_data(restart,trim(c_shortname(i))//'_protected_C_turnover_accumulated','zfull', &
                    sc_protected_C_turnover_ptr, i, trim(c_longname(i))//' protected soil carbon turnover', 'year-1')
           call add_tile_data(restart,trim(c_shortname(i))//'_C_turnover_accumulated','zfull', &
                    sc_C_turnover_ptr, i, trim(c_longname(i))//' soil carbon turnover','year-1')
           ! N turnover rates
           call add_tile_data(restart,trim(c_shortname(i))//'_protected_N_turnover_accumulated','zfull', &
                    sc_protected_N_turnover_ptr, i, trim(c_longname(i))//' protected soil carbon turnover', 'year-1')
           call add_tile_data(restart,trim(c_shortname(i))//'_N_turnover_accumulated','zfull', &
                    sc_N_turnover_ptr, i, trim(c_longname(i))//' soil carbon turnover','year-1')
        enddo
     case default
        call error_mesg('save_soil_restart', 'unrecognized soil carbon option -- this should never happen', FATAL)
     end select
     call save_land_restart(restart)
     call free_land_restart(restart)
  endif
end subroutine save_soil_restart


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
! compute beta function
! after Manabe (1969), but distributed vertically.
subroutine soil_data_beta ( soil, vegn, soil_beta, soil_water_supply, &
                            soil_uptake_T )
  type(soil_tile_type), intent(in)    :: soil
  type(vegn_tile_type), intent(inout) :: vegn ! inout because cc%uptake_frac is updated
  real, intent(out) :: soil_beta(:) ! relative water availability, used only in VEGN_PHOT_SIMPLE treatment
  real, intent(out) :: soil_water_supply(:) ! max rate of water supply to roots, kg/(indiv s)
  real, intent(out) :: soil_uptake_T(:) ! an estimate of temperature of the water
             ! taken up by transpiration. In case of 'linear' uptake it is an exact
             ! value; in case of 'darcy*' treatments the actual uptake profile
             ! is calculated only in step 2, so the value returned is an estimate

  ! ---- local vars
  integer :: k, l
  real, dimension(num_l) :: &
       uptake_frac_max, & ! normalized root distribution
       vegn_uptake_term, &
       vlc, vsc, & ! volumetric fractions of water and ice in the layer
       u, du ! uptake and its derivative (the latter is not used)
  real :: z  !  soil depth
  type (vegn_cohort_type), pointer :: cc

  do l = 1, num_l
     vlc(l) = max(0., soil%wl(l) / (dens_h2o*dz(l)))
     vsc(l) = max(0., soil%ws(l) / (dens_h2o*dz(l)))
  enddo

  do k = 1, vegn%n_cohorts
     cc=>vegn%cohorts(k)
     call cohort_uptake_profile (cc, dz(1:num_l), uptake_frac_max, vegn_uptake_term )

     do l = 1, num_l
        cc%uptake_frac(l) = uptake_frac_max(l) &
             * max(0.0, min(1.0,(vlc(l)-soil%w_wilt(l))/&
                  (0.75*(soil%w_fc(l)-soil%w_wilt(l)))))
     enddo
     soil_beta(k) = sum(cc%uptake_frac(1:num_l))
     if (soil_beta(k) /= 0) then
        cc%uptake_frac(1:num_l) = cc%uptake_frac(1:num_l) / soil_beta(k)
     else
        cc%uptake_frac(1:num_l) = uptake_frac_max(:)
     endif

     if (lm2) cc%uptake_frac(1:num_l) = uptake_frac_max(:)

     ! calculate total water supply
     select case (uptake_option)
     case(UPTAKE_LINEAR)
        z = 0; soil_water_supply(k) = 0
        do l = 1, num_l
           soil_water_supply(k) = soil_water_supply(k) + &
             vegn_uptake_term(k)*max(0.0,soil%wl(l)/dz(l)-soil%w_wilt(l)*dens_h2o)
           z = z + dz(l)
        enddo
        soil_water_supply(k) = z * soil_water_supply(k)/delta_time
        soil_uptake_T(k) = sum(cc%uptake_frac(1:num_l)*soil%T(1:num_l))
     case(UPTAKE_DARCY2D,UPTAKE_DARCY2D_LIN)
        call darcy2d_uptake ( soil, psi_wilt, vegn%root_distance, cc%root_length, &
             cc%K_r, cc%r_r, u, du )
        soil_water_supply(k) = max(0.0,sum(u))
        soil_uptake_T(k) = soil%uptake_T
     end select
  enddo
end subroutine soil_data_beta


! ============================================================================
! update soil properties explicitly for time step.
! MAY WISH TO INTRODUCE 'UNIVERSAL' SENSITIVITIES FOR SIMPLICITY.
! T-DEPENDENCE OF HYDRAULIC PROPERTIES COULD BE DONE LESS FREQUENTLY.
! integrate soil-heat conduction equation upward from bottom of soil
! to surface, delivering linearization of surface ground heat flux.
subroutine soil_step_1 ( soil, vegn, diag, &
                         soil_T, &
                         soil_E_min, soil_E_max, &
                         soil_rh, soil_rh_psi, soil_liq, soil_ice, soil_subl, soil_tf, &
                         soil_G0, soil_DGDT )
  type(soil_tile_type), intent(inout) :: soil
  type(vegn_tile_type), intent(in)    :: vegn
  type(diag_buff_type), intent(inout) :: diag
  real, intent(out) :: &
       soil_T, &    ! temperature of the upper layer of the soil, degK
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
  real :: psi_for_rh

  if(is_watch_point()) then
     write(*,*) '########### soil_step_1 input ###########'
     __DEBUG1__(soil%tag)
     __DEBUG1__(soil%pars%k_sat_ref)
     __DEBUG1__(soil%pars%psi_sat_ref)
     __DEBUG1__(soil%pars%chb)
     __DEBUG1__(soil%pars%vwc_sat)
!     do l = 1,N_LITTER_POOLS
!        call debug_pool(soil%litter(l), trim(l_shortname(l))//'_litter')
!     enddo
!      do l = 1, num_l
!         write(*,'(i2.2,x)',advance='NO') l
!         call debug_pool(soil%soil_organic_matter(l), '')
!      enddo
     write(*,*) '########### end of soil_step_1 input ###########'
  endif
! ----------------------------------------------------------------------------
! in preparation for implicit energy balance, determine various measures
! of water availability, so that vapor fluxes will not exceed mass limits
! ----------------------------------------------------------------------------

  soil_T = soil%T(1)

  do l = 1, num_l
     vlc(l) = max(0.0, soil%wl(l) / (dens_h2o * dz(l)))
     vsc(l) = max(0.0, soil%ws(l) / (dens_h2o * dz(l)))
  enddo
  ! calculate relative humidity at soil surface
  call soil_data_psi_for_rh ( soil, vlc, vsc, soil%psi, psi_for_rh )
  soil_rh = exp(psi_for_rh*g_RT)
  soil_rh_psi = g_RT*soil_rh

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
     __DEBUG1__(soil_T)
     __DEBUG1__(soil_E_max)
     __DEBUG1__(soil_rh)
     __DEBUG1__(soil_liq)
     __DEBUG1__(soil_ice)
     __DEBUG1__(soil_subl)
     __DEBUG1__(soil_G0)
     __DEBUG1__(soil_DGDT)
!     do l = 1,num_l
!        write(*,'(i2.2,x)',advance='NO') l
!        __DEBUG2__(soil%f(l),soil%e(l))
!     enddo
  endif

  call send_tile_data(id_thermal_cond, thermal_cond, diag)
end subroutine soil_step_1


! ============================================================================
! apply boundary flows to soil water and move soil water vertically.
  subroutine soil_step_2 ( soil, vegn, diag, soil_subl, snow_lprec, snow_hlprec,  &
                           vegn_uptk, &
                           subs_DT, subs_M_imp, subs_evap, &
                           use_tfreeze_in_grnd_latent, &
                           ! output
                           soil_levap, soil_fevap, soil_melt, &
                           soil_lrunf, soil_hlrunf, soil_Ttop, soil_Ctop, &
                           soil_frunf, soil_hfrunf, soil_tr_runf, &
                           DOC_to_atmos)
  type(soil_tile_type), intent(inout) :: soil
  type(vegn_tile_type), intent(in)    :: vegn
  type(diag_buff_type), intent(inout) :: diag
  real, intent(in) :: & ! ZMS assign tentative annotations below with "??"
       soil_subl     ! ?? solution for soil surface sublimation [mm/s]
  real, intent(in) :: &
       snow_lprec, & ! ?? solid / liquid throughfall infiltrating the snow [kg/m2/s]
       snow_hlprec, & ! ?? heat associated with snow_lprec [W/m^2]
       vegn_uptk(:), &  ! vegetation soil water uptake flux [kg/m2 of cohort/s]
       subs_DT,       & ! ?? soil surface layer temperature tendency [K]
       subs_M_imp,       &! rate of phase change of non-evaporated soil water ?? [kg/m2/s]
       subs_evap         ! ?? solution for soil surface evaporation [kg/m2/s]
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
       soil_tr_runf(:), & ! tracer runoff from soil [kgX/m^2/s]
       DOC_to_atmos ! if there is no river DOC tracer, this is the loss of DOC
                    ! through leaching [kgC/(m2 s)], otherwise zero.

  ! ---- local vars ----------------------------------------------------------
  real, dimension(num_l) :: del_t, &! ?? temperature tendency [K]
       psi, &! soil moisture potential wrt local elevation [m]
       DThDP, & ! ?? deriv. of vol. liq. cont. wrt psi [1/m]
       K_z, K_x, &! soil horiz. and vert. hydraulic conductivity, respectively [mm/s]
       DKDP, & ! ?? deriv. of hyd. cond. wrt psi [kg/m^3]
       vlc, & ! volumetric liquid water content [-]
       vsc, & ! volumetric solid water content [-]
       dW_l, & ! tendency of soil water mass [mm]
       DPsi
  real, dimension(num_l+1) :: flow, & ! downwards flow at layer interface above [mm/timestep]
       infilt, flow_accum
  real, dimension(num_l  ) :: div, & ! total divergence of soil water [mm/s]
       div_it    ! divergence of water due to inter-tile flow (incl. to stream)
  ! set in hlsp_hydrology_1 [mm/s]
  real, dimension(num_l  ) :: hdiv_it, &! divergence of heat due to inter-tile water flow [W/m^2]
       div_bf, & ! baseflow [mm/s]
       div_if, & ! interlow [mm/s]
       div_al, & ! div from active layer [mm/s]
       dq, div_active, air_depth, macro_frac, extra

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
  integer :: n_iter, l, l_max_active_layer, i
  real :: &
       uptake(num_l),   & ! uptake by roots per layer, kg/(m2 s)
       uptake1(num_l),  & ! uptake by roots per layer per individual, kg/(indiv s)
       soil_uptake_frac(num_l), & ! fraction of uptake for each layer, for LM2 mode only
       tot_nindivs,     & ! total number of individuals, for soil_uptake_frac normalization
       transp1,         & ! transpiration per individual, kg/(indiv s)
       uptake_tot,      & ! total uptake, kg/(m2 s)
       uptake_pos,      & ! sum of the positive uptake, kg/(m2 s)
       uptake_T_new, & ! updated average temperature of uptaken water, deg K
       uptake_T_corr,& ! correction for uptake temperature, deg K
       Tu,           & ! temperature of water taken up from (or added to) a layer, deg K
       psi_x0          ! water potential inside roots (in xylem) at zero depth, m
  type(vegn_cohort_type), pointer :: cc
  real :: Theta ! for debug printout only
  integer :: ic ! cohort iterator
  integer :: severity ! for negative wl checking

  ! For testing tridiagonal solution
  real, dimension(num_l)   :: t_soil_tridiag ! soil temperature based on generic tridiagonal solution [K]
  real, dimension(num_l)   :: t_diff ! difference from original advection subroutine [K]

  real :: DOC_leached(N_C_TYPES,num_l), div_DOC_loss(N_C_TYPES,num_l),  &     ! C leaching
         leaflitter_DOC_loss(N_C_TYPES),woodlitter_DOC_loss(N_C_TYPES)        ! Surface litter C leaching loss
  real :: DON_leached(N_C_TYPES,num_l), div_DON_loss(N_C_TYPES,num_l),  &     ! N leaching
         leaflitter_DON_loss(N_C_TYPES), woodlitter_DON_loss(N_C_TYPES)       ! Surface litter N leaching loss
  real :: NO3_leached(num_l), div_NO3_loss(num_l), &     ! NO3 leaching
         leaflitter_NO3_loss, woodlitter_NO3_loss        ! Surface litter NO3 leaching loss
  real :: NH4_leached(num_l), div_NH4_loss(num_l), &     ! NH4 leaching
         leaflitter_NH4_loss, woodlitter_NH4_loss        ! Surface litter NH4 leaching loss

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
  real :: bwood  ! woody biomass, [kgC/m2], used only as macroporosity flag

  real :: surf_DOC_loss(N_C_TYPES)! [kg C/m^2] DOC loss from top soil layer to surface runoff due
                                  ! to efflux
  real :: surf_DON_loss(N_C_TYPES)! [kg N/m^2] DON loss from top soil layer to surface runoff due
                                  ! to efflux
  real :: surf_NO3_loss, surf_NH4_loss ! [kg N/m^2] NH4 and NO3 loss from top soil layer to surface runoff due to efflux
  real :: total_C_leaching(num_l),total_DON_leaching(num_l) ! [kg C/m^2/s] net total vertical DOC leaching by layer
  real :: total_DOC_div           ! [kg C/m^2/s] net total DOC divergence loss rate
  real :: total_DON_div           ! [kg N/m^2/s] net total DON divergence loss rate
  real :: total_NO3_div,total_NH4_div ! [kg N/m^2/s] net total inorganic N divergence loss rates

  real, dimension(num_l) :: passive_ammonium_uptake, passive_nitrate_uptake ! Uptake of dissolved mineral N by roots through water uptake
  real, dimension(vegn%n_cohorts) :: passive_N_uptake
  ! --------------------------------------------------------------------------
  div_active(:) = 0.0

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
!     do l = 1,N_LITTER_POOLS
!        call debug_pool(soil%litter(l), trim(l_shortname(l))//'_litter')
!     enddo
  endif

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

  ! ---- extract evap from soil and do implicit melt --------------------
  IF(LM2) THEN
    soil_uptake_frac(:) = 0.0 ; tot_nindivs = 0.0
    do ic = 1,vegn%n_cohorts
       cc=>vegn%cohorts(ic)
       soil_uptake_frac(:) = soil_uptake_frac(:)+cc%uptake_frac(:)*cc%nindivs
       tot_nindivs = tot_nindivs+cc%nindivs
    enddo
    if (tot_nindivs>0) &
         soil_uptake_frac(:) = soil_uptake_frac(:)/tot_nindivs
    do l = 1, num_l
       soil%wl(l) = soil%wl(l) - soil_uptake_frac(l)*soil_levap*delta_time
    enddo
  ELSE
    soil%wl(1) = soil%wl(1) - soil_levap*delta_time
    soil%ws(1) = soil%ws(1) - soil_fevap*delta_time
  ENDIF
  hcap = soil%heat_capacity_dry(1)*dz(1) + clw*soil%wl(1) + csw*soil%ws(1)
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

  ! calculate actual vertical distribution of uptake
  uptake(:) = 0.0
  do ic = 1, vegn%n_cohorts
     cc=>vegn%cohorts(ic)
     if (cc%nindivs>0) then
        transp1 = vegn_uptk(ic)*cc%layerfrac/cc%nindivs ! transpiration per individual
     else
        transp1 = 0.0
     endif
     select case(uptake_option)
     case ( UPTAKE_LINEAR )
        n_iter = 0
        uptake1(:) = cc%uptake_frac(:)*transp1
     case ( UPTAKE_DARCY2D, UPTAKE_DARCY2D_LIN )
        ! for Darcy-flow uptake, find the root water potential to satify actual
        ! transpiration by the vegetation
        call darcy2d_uptake_solver     (soil, transp1, vegn%root_distance, &
                cc%root_length, cc%K_r, cc%r_r, &
                uptake1, psi_x0, n_iter)
        ! Solution provides psi inside the skin, given uptake and K_r for each level
        ! This calculates effective psi outside the skin (root-soil interface)
        ! across all levels using single Kri for use in cavitation calculations.
     end select
     if (is_watch_point()) then
        __DEBUG3__(transp1,sum(uptake1(1:num_l)),n_iter)
     endif
     uptake = uptake + uptake1*cc%nindivs

     ! Passive nitrogen uptake by roots, equal to N concentration times root water uptake
     ! units of uptake1: kg/indiv/s
     ! units of wl: mm = kg/m2
     ! units of N: kg/m2
     ! N/wl -> kg N/kg H2O
     ! units of delta_time: s
     ! units of passive_ammonium_uptake, passive_nitrate_uptake, passive_N_uptake: kgN/m2/timestep
     where(soil%wl(1:num_l)>1.0e-4)
        passive_ammonium_uptake(1:num_l) = min(soil%soil_organic_matter(1:num_l)%ammonium,max(0.0,uptake1(1:num_l)*soil%soil_organic_matter(1:num_l)%ammonium*ammonium_solubility/soil%wl(1:num_l)*cc%nindivs*delta_time))
        passive_nitrate_uptake(1:num_l) = min(soil%soil_organic_matter(1:num_l)%nitrate,max(0.0,uptake1(1:num_l)*soil%soil_organic_matter(1:num_l)%nitrate*nitrate_solubility/soil%wl(1:num_l)*cc%nindivs*delta_time))
     elsewhere
        passive_ammonium_uptake(1:num_l)=0.0
        passive_nitrate_uptake(1:num_l)=0.0
     end where
     soil%soil_organic_matter(1:num_l)%ammonium=soil%soil_organic_matter(:)%ammonium-passive_ammonium_uptake(1:num_l)
     soil%soil_organic_matter(1:num_l)%nitrate=soil%soil_organic_matter(:)%nitrate-passive_nitrate_uptake(1:num_l)
     passive_N_uptake(ic) = sum(passive_ammonium_uptake + passive_nitrate_uptake)
     if (cc%nindivs>0) &
           cc%stored_N = cc%stored_N + passive_N_uptake(ic)/cc%nindivs
  enddo

  call send_cohort_data(id_passive_N_uptake,diag,vegn%cohorts(1:vegn%n_cohorts),passive_N_uptake/dt_fast_yr,weight=vegn%cohorts(1:vegn%n_cohorts)%nindivs, op=OP_SUM)

  uptake_pos = sum(uptake(:),mask=uptake(:)>0)
  if (uptake_option/=UPTAKE_LINEAR.and.uptake_pos > 0) then
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

  if (is_watch_point())then
     write(*,*) ' ##### soil_step_2 checkpoint 2.1 #####'
     __DEBUG1__(vegn_uptk)
     __DEBUG2__(sum(vegn_uptk),sum(uptake))
     do l = 1,num_l
        write(*,'(i2.2,x)',advance='NO') l
        call dpri('uptake=',uptake(l))
        call dpri('dwl=',-uptake(l)*delta_time)
        call dpri('wl=',soil%wl(l))
        call dpri('new wl=',soil%wl(l) - uptake(l)*delta_time)
        write(*,*)
     enddo
  endif

  call send_tile_data(id_uptk_n_iter, real(n_iter), diag)
  call send_tile_data(id_uptk, uptake, diag)
  call send_tile_data(id_psi_x0, psi_x0, diag)

  ! update temperature and water content of soil due to root uptake processes
  do l = 1, num_l
     ! calculate the temperature of water that is taken from the layer (or added
     ! to the layer), including energy balance correction
     if (uptake(l) > 0) then
        Tu = soil%T(l) + uptake_T_corr
     else
        Tu = soil%uptake_T + uptake_T_corr
     endif
     ! heat capacity of the layer
     hcap = soil%heat_capacity_dry(l)*dz(l) + clw*soil%wl(l) + csw*soil%ws(l)

     soil%T(l) = soil%T(l) - &
          uptake(l)*delta_time*clw*( Tu-soil%T(l) ) / &
          ( hcap - uptake(l)*delta_time*clw )
     soil%wl(l) = soil%wl(l) - uptake(l)*delta_time
  enddo

  if(is_watch_point()) then
     write(*,*) ' ##### soil_step_2 checkpoint 3 #####'
     do l = 1, num_l
        write(*,'(i2.2,x)',advance='NO')l
        call dpri('T =',soil%T(l))
        call dpri('Th=',(soil%ws(l)+soil%wl(l))/(dens_h2o*dz(l)))
        call dpri('wl=',soil%wl(l))
        call dpri('ws=',soil%ws(l))
        write(*,*)
     enddo
  endif

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
  enddo

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
     call soil_data_gw_hydraulics_ar5(soil, storage_2, gw_flux, sat_area_frac)
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
        ziph = sum(dz(1:num_l))
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

  if (is_watch_point()) then
     do l = 1, num_l
        write(*,'(i2.2,x)',advance='NO')l
        __DEBUG5__(div_active(l),div_bf(l),div_if(l),div_al(l),div(l))
     enddo
     do l = 1, num_l
        write(*,'(i2.2,x)',advance='NO')l
        __DEBUG3__(vsc(l),psi(l),dz(l))
     enddo
     __DEBUG1__(lrunf_bf)
     call dpri('tau_gw',soil%pars%tau_groundwater); write(*,*)
     __DEBUG1__(dens_h2o)
  endif

  wl_before(1:num_l) = soil%wl(1:num_l)

  ! ---- soil-water flow ----------------------------------------------------
  IF (LM2) THEN
     flow(1) = 0
     do l = 1, num_l
        infilt(l) = soil_uptake_frac(l)*lprec_eff *delta_time
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
        call RICHARDS_clean(soil, psi, DThDP, K_z, DKDP, div, &
               lprec_eff, Dpsi_min, Dpsi_max, delta_time, &
               dPsi, dW_l, flow, lrunf_ie)
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
           endif
           call check_conservation('soil_step_2: Richards Eqn. Diagnostics, dW_l', 'Water', w1, w2, &
                wthresh, WARNING)
        enddo
#endif
     ENDIF
  ENDIF

  ! Check for negative wl
  do l = 1, num_l
     if (soil%wl(l) .ge. 0.) cycle
     if (soil%wl(l) / (dens_h2o*dz(l)*soil%pars%vwc_sat) .ge. thetathresh) cycle
     if (.not. allow_neg_wl) then
        call check_var_range(soil%wl(l)/(dens_h2o*dz(l)*soil%pars%vwc_sat), 0.0, HUGE(1.0), 'soil_step_2', 'degree of saturation', FATAL)
     else
        if (verbose) &
            call check_var_range(soil%wl(l)/(dens_h2o*dz(l)*soil%pars%vwc_sat), 0.0, HUGE(1.0), 'soil_step_2', 'degree of saturation', WARNING)
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

  if (is_watch_point()) then
     write(*,*) ' ***** soil_step_2 checkpoint 3.4 ***** '
     __DEBUG1__(tfreeze)
     __DEBUG1__(tflow)
     __DEBUG1__(snow_hlprec)
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
     hlrunf_ie = hlprec_eff - (flow(1)/delta_time)*clw*(soil%T(1)-tfreeze)
  else
     hlrunf_ie = 0.
  endif

  ! Initialize for use in output below
  macro_inf = 0.
  extra_cum = 0.
  ! ---- allow infiltration-excess runoff to enter soil via macroporosity
  bwood = vegn_tile_bwood(vegn)
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
        write(*,'(i2.2,x)',advance='NO')l
        call dpri('T=',soil%T(l))
        call dpri('Th=',(soil%ws(l)+soil%wl(l))/(dens_h2o*dz(l)))
        call dpri('wl=',soil%wl(l))
        call dpri('ws=',soil%ws(l))
        call dpri('gw=',soil%groundwater(l))
        write(*,*)
     enddo
     call debug_pool(soil%litter(LEAF), 'leaf_litter')
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
  soil_Ctop = soil%heat_capacity_dry(1)*dz(1) + clw*soil%wl(1) + csw*soil%ws(1)

  soil%psi=psi+dPsi
!  if (do_component_balchecks) then
     ! Sum total water mass at end of soil_step_2
     call soil_tile_stock_pe (soil, sliq, sice )
     wsum2 = sliq + sice

     ! Subtract influxes and add outfluxes
     wsum2 = wsum2 - delta_time *  snow_lprec &
          + delta_time * ( subs_evap + sum(uptake) + soil_lrunf + soil_frunf)

! slm: call check_conservation('soil_mod: soil_step_2', 'Water', wsum1, wsum2, wthresh, FATAL)
! endif


   if (is_watch_point()) then
      write(*,*)'##### soil_step_2 checkpoint 6 #####'
      __DEBUG1__(flow)
      __DEBUG1__(div)
      __DEBUG1__(wl_before)
      __DEBUG1__(gw_option)
      do l = 1,N_LITTER_POOLS
         call debug_pool(soil%litter(l), trim(l_shortname(l))//'_litter')
      enddo
      do l = 1, num_l
         write(*,'(i2.2,x)',advance='NO') l
         call debug_pool(soil%soil_organic_matter(l), '')
      enddo
      do l = 1, size(soil%div_hlsp_DOC,2)
         __DEBUG1__(soil%div_hlsp_DOC(:,l))
      enddo
   endif

!New version that combines the two leaching steps and should do a better job of moving DOC from litter layer
!Note: fine wood litter currently not included, just because we haven't implemented anything with it anywhere
!ZMS Edited to allow for tiled fluxes. Also pass in water content before Richards.
   if (gw_option == GW_TILED) then
      call tracer_leaching_with_litter(soil=soil%soil_organic_matter(:),&
                                        wl=wl_before,&
                                        leaflitter=soil%litter(LEAF),&
                                        woodlitter=soil%litter(CWOOD),&
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
                                      leaflitter=soil%litter(LEAF),&
                                      woodlitter=soil%litter(CWOOD),&
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
   ! Diagnostic. Later pass this back to land_model for transfer to rivers.
   total_DOC_div = sum(surf_DOC_loss(:))
   total_DON_div = sum(surf_DON_loss(:))
   total_NO3_div = surf_NO3_loss + sum(div_NO3_loss(1:num_l))
   total_NH4_div = surf_NH4_loss + sum(div_NH4_loss(1:num_l))
   do l=1,num_l
      total_DOC_div = total_DOC_div + sum(div_DOC_loss(:,l))
      total_DON_div = total_DON_div + sum(div_DON_loss(:,l))
   end do

   !FIXME BNS: What if there is net flow of nitrogen into tile from other hillslope tiles?
   soil%gross_nitrogen_flux_out_of_tile = soil%gross_nitrogen_flux_out_of_tile + total_DON_div+total_NO3_div+total_NH4_div

   total_DOC_div = total_DOC_div/delta_time
   total_DON_div = total_DON_div/delta_time
   total_NO3_div = total_NO3_div/delta_time
   total_NH4_div = total_NH4_div/delta_time
   if (i_river_DOC/=NO_TRACER) then
       soil_tr_runf(i_river_DOC) = total_DOC_div
       DOC_to_atmos = 0.0
   else
       ! if we don't have DOC river tracer, DOC loss goes directly to the atmosphere
       DOC_to_atmos = total_DOC_div
   endif
   if (i_river_DON/=NO_TRACER) &
       soil_tr_runf(i_river_DON) = total_DON_div
   if (i_river_NO3/=NO_TRACER) &
       soil_tr_runf(i_river_NO3) = total_NO3_div
   if (i_river_NH4/=NO_TRACER) &
       soil_tr_runf(i_river_NH4) = total_NH4_div

   if (is_watch_point()) then
      write(*,*)'##### soil_step_2 checkpoint 7 #####'
      do l = 1,N_LITTER_POOLS
         call debug_pool(soil%litter(l),trim(l_shortname(l))//'Litter')
      enddo
      __DEBUG3__(leaflitter_DOC_loss,woodlitter_DOC_loss,total_DOC_div)
      do l = 1, num_l
         write(*,'(i2.2,x)',advance='NO') l
         call debug_pool(soil%soil_organic_matter(l), '')
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

  ! CMOR variables
  if (id_mrlsl > 0) call send_tile_data(id_mrlsl, soil%wl+soil%ws, diag)
  if (id_mrsfl > 0) call send_tile_data(id_mrsfl, soil%ws, diag)
  if (id_mrsll > 0) call send_tile_data(id_mrsll, soil%wl, diag)
  if (id_mrsol > 0) call send_tile_data(id_mrsol, soil%wl+soil%ws, diag)
  if (id_mrso > 0)  call send_tile_data(id_mrso,  sum(soil%wl+soil%ws), diag)
  if (id_mrsos > 0) call send_tile_data(id_mrsos, sum((soil%wl+soil%ws)*mrsos_weight), diag)
  if (id_mrs1mLut > 0) call send_tile_data(id_mrs1mLut, sum((soil%wl+soil%ws)*mrs1m_weight), diag)
  if (id_mrfso > 0) call send_tile_data(id_mrfso, sum(soil%ws), diag)
  if (id_mrlso > 0) call send_tile_data(id_mrlso, sum(soil%wl), diag)
  call send_tile_data(id_mrros, lrunf_ie+lrunf_sn, diag)
  call send_tile_data(id_mrro,  lrunf_ie+lrunf_sn+lrunf_bf+lrunf_nu, diag)

  ! ZMS uncomment for back-compatibility with diag tables
  if (gw_option == GW_TILED) then
     call send_tile_data(id_deficit, deficit, diag)
     call send_tile_data(id_sat_depth, depth_to_wt_3, diag)
     call send_tile_data(id_sat_dept2, depth_to_wt2_3, diag)
     call send_tile_data(id_z_cap, depth_to_cf_3, diag)
     if (depth_to_wt_2a .ge. -0.5) &
          call send_tile_data(id_sat_depth, depth_to_wt_2a, diag)
  endif
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
  call send_tile_data(id_div,div,diag)

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

  do i = 1, N_C_TYPES
    call send_tile_data(id_litter_C_leaching(LEAF,i),leaflitter_DOC_loss(i)/delta_time,diag)
    call send_tile_data(id_litter_C_leaching(CWOOD,i),woodlitter_DOC_loss(i)/delta_time,diag)
    call send_tile_data(id_C_leaching(i), DOC_leached(i,:)/delta_time,diag)
    call send_tile_data(id_litter_DON_leaching(LEAF,i),leaflitter_DON_loss(i)/delta_time,diag)
    call send_tile_data(id_litter_DON_leaching(CWOOD,i),woodlitter_DON_loss(i)/delta_time,diag)
    call send_tile_data(id_DON_leaching(i), DON_leached(i,:)/delta_time,diag)
  enddo

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
  do l=1,num_l
     total_C_leaching(l) = sum(DOC_leached(:,l))/delta_time
     total_DON_leaching(l) = sum(DON_leached(:,l))/delta_time
  end do
  call send_tile_data(id_total_C_leaching, total_C_leaching, diag)
  call send_tile_data(id_total_ON_leaching, total_DON_leaching, diag)
  call send_tile_data(id_NO3_leaching, NO3_leached/delta_time, diag)
  call send_tile_data(id_NH4_leaching, NH4_leached/delta_time, diag)
  call send_tile_data(id_litter_total_C_leaching(LEAF),sum(leaflitter_DOC_loss)/delta_time,diag)
  call send_tile_data(id_litter_total_ON_leaching(LEAF),sum(leaflitter_DON_loss)/delta_time,diag)
  call send_tile_data(id_litter_NO3_leaching(LEAF),leaflitter_NO3_loss/delta_time,diag)
  call send_tile_data(id_litter_NH4_leaching(LEAF),leaflitter_NH4_loss/delta_time,diag)
  call send_tile_data(id_litter_total_C_leaching(CWOOD),sum(woodlitter_DOC_loss)/delta_time,diag)
  call send_tile_data(id_litter_total_ON_leaching(CWOOD),sum(woodlitter_DON_loss)/delta_time,diag)
  call send_tile_data(id_litter_NO3_leaching(CWOOD),woodlitter_NO3_loss/delta_time,diag)
  call send_tile_data(id_litter_NH4_leaching(CWOOD),woodlitter_NH4_loss/delta_time,diag)
  if (gw_option == GW_TILED) then
     call send_tile_data(id_surf_DOC_loss, sum(surf_DOC_loss(:))/delta_time,diag)
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

  real :: soil_C(N_C_TYPES, num_l),      soil_N(N_C_TYPES, num_l), &
          dissolved_C(N_C_TYPES, num_l), dissolved_N(N_C_TYPES, num_l), &
          protected_C(N_C_TYPES,num_l),  protected_N(N_C_TYPES,num_l), &
          livemic_C(num_l), livemic_N(num_l), &
          layer_C(num_l),   layer_N(num_l)
  real :: litter_C(N_C_TYPES), litter_dissolved_C(N_C_TYPES), litter_protected_C(N_C_TYPES), litter_livemic_C, litter_total_C, &
          litter_N(N_C_TYPES), litter_dissolved_N(N_C_TYPES), litter_protected_N(N_C_TYPES), litter_livemic_N, litter_total_N
  integer :: i, k, l, ncohorts(num_l), litter_ncohorts
  real :: total_C(N_C_TYPES), total_livemic_C, total_prot_C(N_C_TYPES), total_diss_C(N_C_TYPES), &
          total_N(N_C_TYPES), total_livemic_N, total_prot_N(N_C_TYPES), total_diss_N(N_C_TYPES)
  real :: total_NO3, total_NH4
  real :: total_litter_C ! total C in all litter, for diagnostics

  select case (soil_carbon_option)
  case (SOILC_CENTURY, SOILC_CENTURY_BY_LAYER)
     call send_tile_data(id_fsc, sum(soil%fast_soil_C(:)), diag)
     call send_tile_data(id_ssc, sum(soil%slow_soil_C(:)), diag)
     call send_tile_data(id_soil_C(C_FAST), soil%fast_soil_C(:)/dz(1:num_l), diag)
     call send_tile_data(id_soil_C(C_SLOW), soil%slow_soil_C(:)/dz(1:num_l), diag)
     ! --- CMOR vars
     if (id_csoilfast>0)   call send_tile_data(id_csoilfast,   sum(soil%fast_soil_C(:)), diag)
     if (id_csoilmedium>0) call send_tile_data(id_csoilmedium, sum(soil%slow_soil_C(:)), diag)
     call send_tile_data(id_csoilslow, 0.0, diag)
     if (id_csoil>0)       call send_tile_data(id_csoil, sum(soil%fast_soil_C(:))+sum(soil%slow_soil_C(:)), diag)
     ! --- end of CMOR vars
  case (SOILC_CORPSE, SOILC_CORPSE_N)
!     total_carbon=0.0

     do l = 1,num_l
        call poolTotals1 ( soil%soil_organic_matter(l), ncohorts=ncohorts(l), &
            litterC=soil_C(:,l), livemicC=livemic_C(l), protectedC=protected_C(:,l), dissolvedC=dissolved_C(:,l), totalC=layer_C(l), &
            litterN=soil_N(:,l), livemicN=livemic_N(l), protectedN=protected_N(:,l), dissolvedN=dissolved_N(:,l), totalN=layer_N(l)  )
     enddo
     total_C(:)      = sum(soil_C(:,:),2)
     total_livemic_C = sum(livemic_C)
     total_prot_C    = sum(protected_C(:,:),2)
     total_diss_C    = sum(dissolved_C(:,:),2)
     total_N(:)      = sum(soil_N(:,:),2)
     total_livemic_N = sum(livemic_N)
     total_prot_N    = sum(protected_N(:,:),2)
     total_diss_N    = sum(dissolved_N(:,:),2)
     total_NO3       = sum(soil%soil_organic_matter(1:num_l)%nitrate)
     total_NH4       = sum(soil%soil_organic_matter(1:num_l)%ammonium)


     call send_tile_data(id_nsoilcohorts, real(ncohorts), diag)
     do i = 1, N_C_TYPES
        call send_tile_data(id_soil_C(i),           soil_C(i,:)/dz(1:num_l),      diag)
        call send_tile_data(id_soil_N(i),           soil_N(i,:)/dz(1:num_l),      diag)
        call send_tile_data(id_soil_protected_C(i), protected_C(i,:)/dz(1:num_l), diag)
        call send_tile_data(id_soil_protected_N(i), protected_N(i,:)/dz(1:num_l), diag)
        call send_tile_data(id_soil_dissolved_C(i), dissolved_C(i,:)/dz(1:num_l), diag)
        call send_tile_data(id_soil_dissolved_N(i), dissolved_N(i,:)/dz(1:num_l), diag)
     enddo
     call send_tile_data(id_protected_C, sum(protected_C,1)/dz(1:num_l), diag)
     call send_tile_data(id_protected_N, sum(protected_N,1)/dz(1:num_l), diag)
     call send_tile_data(id_livemic_C, livemic_C/dz(1:num_l), diag)
     call send_tile_data(id_livemic_N, livemic_N/dz(1:num_l), diag)
     call send_tile_data(id_total_C_layered, layer_C(:)/dz(1:num_l), diag)
     call send_tile_data(id_total_N_layered, layer_N(:)/dz(1:num_l), diag)

     call send_tile_data(id_fast_DOC_div_loss,    soil%fast_DOC_leached,    diag)
     call send_tile_data(id_slow_DOC_div_loss,    soil%slow_DOC_leached,    diag)
     call send_tile_data(id_deadmic_DOC_div_loss, soil%deadmic_DOC_leached, diag)
     call send_tile_data(id_fast_DON_div_loss,    soil%fast_DON_leached,    diag)
     call send_tile_data(id_slow_DON_div_loss,    soil%slow_DON_leached,    diag)
     call send_tile_data(id_deadmic_DON_div_loss, soil%deadmic_DON_leached, diag)

     call send_tile_data(id_soil_NO3, soil%soil_organic_matter(1:num_l)%nitrate/dz(1:num_l),diag)
     call send_tile_data(id_soil_NH4, soil%soil_organic_matter(1:num_l)%ammonium/dz(1:num_l),diag)
     call send_tile_data(id_total_NO3, total_NO3, diag)
     call send_tile_data(id_total_NH4, total_NH4, diag)

     ! leaf litter diagnostics
     total_litter_C = 0.0
     do k = 1, N_LITTER_POOLS
        call poolTotals1 (soil%litter(k), ncohorts=litter_ncohorts, &
            litterC=litter_C(:), livemicC=litter_livemic_C, protectedC=litter_protected_C(:), dissolvedC=litter_dissolved_C(:), totalC=litter_total_C, &
            litterN=litter_N(:), livemicN=litter_livemic_N, protectedN=litter_protected_N(:), dissolvedN=litter_dissolved_N(:), totalN=litter_total_N  )
        total_C(:)      = total_C(:) + litter_C(:)
        total_livemic_C = total_livemic_C + litter_livemic_C
        total_diss_C    = total_diss_C + sum(soil%litter(k)%dissolved_carbon(:))
        total_prot_C    = total_prot_C + sum(litter_protected_C(:))
        total_N(:)      = total_N(:) + litter_N(:)
        total_livemic_N = total_livemic_N + litter_livemic_N
        total_diss_N    = total_diss_N + sum(soil%litter(k)%dissolved_nitrogen(:))
        total_prot_N    = total_prot_N + sum(litter_protected_N(:))
        total_NO3       = total_NO3 + soil%litter(k)%nitrate
        total_NH4       = total_NO3 + soil%litter(k)%ammonium
        total_litter_C  = total_litter_C + litter_total_C

        call send_tile_data(id_nlittercohorts(k), real(litter_ncohorts), diag)
        call send_tile_data(id_litter_livemic_C(k), litter_livemic_C, diag)
        call send_tile_data(id_litter_livemic_N(k), litter_livemic_N, diag)
        call send_tile_data(id_litter_total_C(k), litter_total_C, diag)
        call send_tile_data(id_litter_total_N(k), litter_total_N, diag)
        call send_tile_data(id_litter_nitrate(k), soil%litter(k)%nitrate, diag)
        call send_tile_data(id_litter_ammonium(k), soil%litter(k)%ammonium, diag)
        do i = 1, N_C_TYPES
           call send_tile_data(id_litter_C(k,i), litter_C(i), diag)
           call send_tile_data(id_litter_N(k,i), litter_N(i), diag)
           call send_tile_data(id_litter_protected_C(k,i), litter_protected_C(i), diag)
           call send_tile_data(id_litter_protected_N(k,i), litter_protected_N(i), diag)
           call send_tile_data(id_litter_dissolved_C(k,i), soil%litter(k)%dissolved_carbon(i), diag)
           call send_tile_data(id_litter_dissolved_N(k,i), soil%litter(k)%dissolved_nitrogen(i), diag)
        enddo
        ! CMOR diagnostics
        if (k==CWOOD) call send_tile_data(id_cLitterCwd, litter_total_C, diag)
     enddo

     ! diagnostic of totals
     call send_tile_data(id_fsc, total_C(C_FAST), diag)
     call send_tile_data(id_fsN, total_N(C_FAST), diag)
     call send_tile_data(id_ssc, total_C(C_SLOW), diag)
     call send_tile_data(id_ssN, total_N(C_SLOW), diag)
     call send_tile_data(id_deadmic_total_C, total_C(C_MIC), diag)
     call send_tile_data(id_deadmic_total_N, total_N(C_MIC), diag)
     call send_tile_data(id_livemic_total_C, total_livemic_C, diag)
     call send_tile_data(id_livemic_total_N, total_livemic_N, diag)
     call send_tile_data(id_protected_total_C, sum(total_prot_C), diag)
     call send_tile_data(id_protected_total_N, sum(total_prot_N), diag)
     call send_tile_data(id_dissolved_total_C, sum(total_diss_C), diag)
     call send_tile_data(id_dissolved_total_N, sum(total_diss_N), diag)
     call send_tile_data(id_total_soil_C, sum(total_C+total_diss_C+total_prot_C)+total_livemic_C, diag)
     call send_tile_data(id_total_soil_N, sum(total_N+total_diss_N+total_prot_N)+total_livemic_N, diag)
     ! --- CMOR vars
     call send_tile_data(id_csoilfast,   total_C(C_FAST), diag)
     call send_tile_data(id_csoilmedium, total_C(C_SLOW), diag)
     call send_tile_data(id_csoilslow,   0.0,             diag)
     call send_tile_data(id_csoil, total_C(C_FAST)+total_C(C_SLOW), diag)
     call send_tile_data(id_cLitter, total_litter_C, diag)
     ! --- end of CMOR vars
  case default
     call error_mesg('soil_step_3','unrecognized soil carbon option -- this should never happen', FATAL)
  end select

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
     call error_mesg('Dsdt','unrecognized soil carbon option -- this should never happen', FATAL)
  end select
  ! CMOR diag
  call send_tile_data(id_rh, vegn%rh/seconds_per_year, diag)
end subroutine Dsdt


! ============================================================================
subroutine Dsdt_CORPSE(vegn, soil, diag)
  type(vegn_tile_type), intent(inout) :: vegn
  type(soil_tile_type), intent(inout) :: soil
  type(diag_buff_type), intent(inout) :: diag

  real :: leaflitter_deadmic_C_produced, leaflitter_deadmic_N_produced
  real :: fineWoodlitter_deadmic_C_produced, fineWoodlitter_deadmic_N_produced
  real :: coarseWoodlitter_deadmic_C_produced, coarseWoodlitter_deadmic_N_produced
  real, dimension(N_C_TYPES) :: &
     leaflitter_protected_C_produced, leaflitter_protected_C_turnover_rate, &
     leaflitter_protected_N_produced, leaflitter_protected_N_turnover_rate, &
     fineWoodlitter_protected_C_produced, fineWoodlitter_protected_C_turnover_rate, &
     fineWoodlitter_protected_N_produced, fineWoodlitter_protected_N_turnover_rate, &
     coarseWoodlitter_protected_C_produced, coarseWoodlitter_protected_C_turnover_rate, &
     coarseWoodlitter_protected_N_produced, coarseWoodlitter_protected_N_turnover_rate
  real :: &
     protected_C_produced(N_C_TYPES,num_l), protected_C_turnover_rate(N_C_TYPES,num_l), &
     protected_N_produced(N_C_TYPES,num_l), protected_N_turnover_rate(N_C_TYPES,num_l)
  real,dimension(N_LITTER_POOLS) :: litter_nitrif, litter_denitrif, litter_N_mineralization, litter_N_immobilization
  real :: deadmic_C_produced(num_l)
  real :: deadmic_N_produced(num_l)
  real :: soil_nitrif(num_l), soil_denitrif(num_l), soil_N_mineralization(num_l), soil_N_immobilization(num_l)
  real :: litter_C_loss_rate(N_C_TYPES)
  real :: litter_N_loss_rate(N_C_TYPES)
  real :: C_loss_rate(size(soil%soil_organic_matter),N_C_TYPES)
  real :: N_loss_rate(size(soil%soil_organic_matter),N_C_TYPES)
  real, dimension(size(soil%soil_organic_matter)) :: decomp_T,decomp_theta,ice_porosity
  real :: A          (size(soil%soil_organic_matter)) ! decomp rate reduction due to moisture and temperature

  integer :: badCohort   ! For soil carbon pool carbon balance and invalid number check
  integer :: i,k
  real :: CO2prod
  integer :: point_i,point_j,point_k,point_face

  A(:) = A_function(soil%T(:), soil_theta(soil))
  decomp_T = soil%T(:)
  decomp_theta = soil_theta(soil)
  ice_porosity = soil_ice_porosity(soil)

  vegn%rh=0.0

  !  First surface litter is decomposed
  do k = 1,N_LITTER_POOLS
     call update_pool(pool=soil%litter(k),T=decomp_T(1),theta=decomp_theta(1),air_filled_porosity=1.0-(decomp_theta(1)+ice_porosity(1)),&
            liquid_water=soil%wl(1),frozen_water=soil%ws(1),dt=dt_fast_yr,layerThickness=dz(1),&
            C_loss_rate=litter_C_loss_rate, CO2prod=CO2prod, &
            N_loss_rate=litter_N_loss_rate, &
            deadmic_C_produced=leaflitter_deadmic_C_produced, protected_C_produced=leaflitter_protected_C_produced, protected_turnover_rate=leaflitter_protected_C_turnover_rate, &
            deadmic_N_produced=leaflitter_deadmic_N_produced, protected_N_produced=leaflitter_protected_N_produced, protected_N_turnover_rate=leaflitter_protected_N_turnover_rate, &
            badCohort=badCohort,&
            nitrification=litter_nitrif(k), denitrification=litter_denitrif(k),&
            N_mineralization=litter_N_mineralization(k), N_immobilization=Litter_N_immobilization(k))
     IF (badCohort.ne.0) THEN
        WRITE (*,*) 'T=',decomp_T(1),'theta=',decomp_theta(1),'dt=',dt_fast_yr
        call land_error_message('Dsdt: Found bad cohort in '//trim(l_longname(k))//' litter.',FATAL)
     ENDIF
     vegn%rh=vegn%rh + CO2prod/dt_fast_yr ! accumulate loss of C to atmosphere
     ! NOTE that the first layer of C_loss_rate and N_loss_rate are used as buffers
     ! for litter diagnostic output.
     do i = 1, N_C_TYPES
        call send_tile_data(id_litter_rsoil_C(k,i), litter_C_loss_rate(i), diag)
        call send_tile_data(id_litter_rsoil_N(k,i), litter_N_loss_rate(i), diag)
     enddo
     ! for budget check
     vegn%fsc_out     = vegn%fsc_out     + litter_C_loss_rate(C_FAST)*dt_fast_yr
     vegn%ssc_out     = vegn%ssc_out     + litter_C_loss_rate(C_SLOW)*dt_fast_yr
     vegn%deadmic_out = vegn%deadmic_out + litter_C_loss_rate(C_MIC) *dt_fast_yr
  enddo


  ! Next we have to go through layers and decompose the soil carbon pools
  do k=1,num_l
      call update_pool(pool=soil%soil_organic_matter(k),T=decomp_T(k),theta=decomp_theta(k),air_filled_porosity=1.0-(decomp_theta(k)+ice_porosity(k)),&
                liquid_water=soil%wl(k),frozen_water=soil%ws(k),dt=dt_fast_yr,layerThickness=dz(k),&
                C_loss_rate=C_loss_rate(k,:), &
                CO2prod=CO2prod, &
                N_loss_rate=N_loss_rate(k,:), &
                deadmic_C_produced=deadmic_C_produced(k), protected_C_produced=protected_C_produced(:,k), &
                protected_turnover_rate=protected_C_turnover_rate(:,k), &
                deadmic_N_produced=deadmic_N_produced(k), protected_N_produced=protected_N_produced(:,k), &
                protected_N_turnover_rate=protected_N_turnover_rate(:,k), &
                badCohort=badCohort,&
                nitrification=soil_nitrif(k),denitrification=soil_denitrif(k),&
                N_mineralization=soil_N_mineralization(k),N_immobilization=soil_N_immobilization(k))
    IF (badCohort.ne.0) THEN
        WRITE (*,*) 'T=',decomp_T(k),'theta=',decomp_theta(k),'dt=',dt_fast_yr
        call land_error_message('Dsdt: Found bad cohort in layer'//trim(string(k))//' of soil carbon.',FATAL)
    ENDIF

    vegn%rh=vegn%rh + CO2prod/dt_fast_yr ! accumulate loss of C to atmosphere
  enddo
  do i = 1, N_C_TYPES
     call send_tile_data(id_rsoil_C(i), C_loss_rate(:,i)/dz(1:num_l), diag)
     call send_tile_data(id_rsoil_N(i), N_loss_rate(:,i)/dz(1:num_l), diag)
  enddo
  ! for budget check
  vegn%fsc_out     = vegn%fsc_out     + sum(C_loss_rate(C_FAST,:))*dt_fast_yr
  vegn%ssc_out     = vegn%ssc_out     + sum(C_loss_rate(C_SLOW,:))*dt_fast_yr
  vegn%deadmic_out = vegn%deadmic_out + sum(C_loss_rate(C_MIC,:)) *dt_fast_yr


  ! accumulate decomposition rate reduction for the soil carbon restart output
  soil%asoil_in(:) = soil%asoil_in(:) + A(:)

  soil%gross_nitrogen_flux_out_of_tile = soil%gross_nitrogen_flux_out_of_tile + (sum(soil_denitrif)+sum(litter_denitrif))

  ! TODO: arithmetic averaging of A does not seem correct; we need to invent something better,
  !       e.g. weight it with the carbon loss, or something like that

  ! ---- diagnostic section
  call send_tile_data(id_rsoil, vegn%rh, diag)
  ! TODO: arithmetic averaging of A does not seem correct; we need to invent something better,
  !       e.g. weight it with the carbon loss, or something like that
  if (id_asoil>0) call send_tile_data(id_asoil, sum(A(:))/size(A(:)), diag)

  if (id_total_denitrification_rate>0) call send_tile_data(id_total_denitrification_rate, &
             (sum(soil_denitrif)+sum(litter_denitrif))/dt_fast_yr,diag)
  if (id_soil_denitrification_rate>0) call send_tile_data(id_soil_denitrification_rate, soil_denitrif/dt_fast_yr/dz, diag)
  if (id_total_N_mineralization_rate>0) call send_tile_data(id_total_N_mineralization_rate, &
             (sum(soil_N_mineralization)+sum(litter_N_mineralization))/dt_fast_yr,diag)
  if (id_total_N_immobilization_rate>0) call send_tile_data(id_total_N_immobilization_rate, &
                (sum(soil_N_immobilization)+sum(litter_N_immobilization))/dt_fast_yr,diag)
  if (id_total_nitrification_rate>0) call send_tile_data(id_total_nitrification_rate, &
          (sum(soil_nitrif)+sum(litter_nitrif))/dt_fast_yr,diag)

  do i = 1, N_C_TYPES
     if (id_negative_litter_C(i)>0) call send_tile_data(id_negative_litter_C(i),soil%neg_litt_C(i),diag)
     if (id_negative_litter_N(i)>0) call send_tile_data(id_negative_litter_N(i),soil%neg_litt_N(i),diag)
  enddo
  if (id_tot_negative_litter_C>0) call send_tile_data(id_tot_negative_litter_C,sum(soil%neg_litt_C),diag)
  if (id_tot_negative_litter_N>0) call send_tile_data(id_tot_negative_litter_N,sum(soil%neg_litt_N),diag)
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

  ! ---- diagnostic section
  call send_tile_data(id_rsoil_C(C_FAST), fast_C_loss(:)/(dz(1:num_l)*dt_fast_yr), diag)
  call send_tile_data(id_rsoil_C(C_SLOW), slow_C_loss(:)/(dz(1:num_l)*dt_fast_yr), diag)
  call send_tile_data(id_rsoil, vegn%rh, diag)

  ! TODO: arithmetic averaging of A does not seem correct; we need to invent something better,
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
  integer :: l, k

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
        soil%T(l) = (h1 * soil%T(l) + h2 * excess_T) / (h1+h2)
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
                  - hlf_factor*hlf*excess_ice          ) / delta_time
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
        write(*,'(i2.2,x)',advance='NO') l
        call dpri('T=',soil%T(l))
        call dpri('Th=',(soil%ws(l)+soil%wl(l))/(dens_h2o*dz(l)))
        call dpri('wl=',soil%wl(l))
        call dpri('ws=',soil%ws(l))
        write(*,*)
     enddo
  endif
end subroutine soil_push_down_excess


! ============================================================================
  subroutine richards_clean(soil, psi, DThDP, hyd_cond, DKDP, div, &
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
        write(*,'(a,i2.2)',advance='NO')'level=', l
        __DEBUG4__(DThDP(l),hyd_cond(l),psi(l),DKDP(l))
     enddo
     do l = 1, num_l-1
        write(*,'(a,i2.2)',advance='NO')'interface=', l
        __DEBUG4__(K(l),DKDPm(l),DKDPp(l),grad(l))
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
     write(*,'(i2.2)',advance='NO') l
     __DEBUG3__(aaa,bbb,ddd)
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
        write(*,'(i2.2)',advance='NO') l
        __DEBUG4__(aaa,bbb,ccc,ddd)
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
      __DEBUG1__(bbb)
      __DEBUG1__(ccc)
      __DEBUG1__(ddd)
      __DEBUG1__(eee(l))
      __DEBUG1__(fff(l))
      __DEBUG1__(dPsi(l))
      __DEBUG1__(dPsi(l))
      __DEBUG1__(Dpsi_min)
      __DEBUG1__(Dpsi_max)
      call error_mesg(module_name, 'b+ce=0 in soil-water equations', FATAL)
  ENDIF

  dPsi(l) = (ddd-ccc*fff(l))/(bbb+ccc*eee(l))
  if(is_watch_point()) then
     __DEBUG1__(bbb)
     __DEBUG1__(ccc)
     __DEBUG1__(ddd)
     __DEBUG1__(eee(l))
     __DEBUG1__(fff(l))
     __DEBUG1__(dPsi(l))
     __DEBUG1__(dPsi(l))
     __DEBUG1__(Dpsi_min)
     __DEBUG1__(Dpsi_max)
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
        write(*,'(i2.2)', advance='NO') l
        call dpri('eee=',eee(l))
        call dpri('fff=',fff(l))
        call dpri('dpsi=',eee(l))
        write(*,*)
     enddo
     __DEBUG1__(DThDP(1))
     __DEBUG1__(K(1))
     __DEBUG1__(grad(1))
     __DEBUG1__(ddd)
     __DEBUG1__(ccc)
     __DEBUG1__(bbb)
     __DEBUG1__(dPsi(1))
     __DEBUG1__(Psi(1))
     __DEBUG1__(div(1))
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
        write(*,'(i2.2)', advance='NO') l
        call dpri('dW_l=',dW_l(l))
        call dpri('flow=',flow(l))
        call dpri('div=', div(l))
        write(*,*)
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
        write(*,'(i2.2)', advance='NO') l
        call dpri('dW_l=',dW_l(l))
        call dpri('flow=',flow(l))
        call dpri('div=', div(l))
        write(*,*)
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
        write(*,'(i2.2)', advance='NO') l
        call dpri('Th=', (soil%ws(l)+soil%wl(l)+dW_l(l))/(dens_h2o*dz(l)))
        call dpri('wl=', soil%wl(l)+dW_l(l))
        call dpri('ws=', soil%ws(l))
        call dpri('dW_l=', dW_l(l))
        call dpri('dPsi=', dPsi(l))
        call dpri('flow=', flow(l))
        write(*,*)
     enddo
  endif
end subroutine richards_clean


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
   u_minus = 1.0; u_plus = 0.0
   where (flow(1:num_l).lt.0.) u_minus = 0.
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
! given soil tile, returns carbon content of various components of litter
subroutine get_soil_litter_C(soil, litter_fast_C, litter_slow_C, litter_deadmic_C)
  type(soil_tile_type), intent(in)  :: soil
  real, intent(out) :: &
     litter_fast_C,    & ! fast litter carbon, [kgC/m2]
     litter_slow_C,    & ! slow litter carbon, [kgC/m2]
     litter_deadmic_C    ! mass of dead microbes in litter, [kgC/m2]

  select case(soil_carbon_option)
  case(SOILC_CENTURY, SOILC_CENTURY_BY_LAYER)
     litter_fast_C    = soil%fast_soil_C(1)
     litter_slow_C    = soil%slow_soil_C(1)
     litter_deadmic_C = 0.0
  case(SOILC_CORPSE, SOILC_CORPSE_N)
     call poolTotals(soil%litter(LEAF),fastC=litter_fast_C,slowC=litter_slow_C,deadMicrobeC=litter_deadmic_C)
  case default
     call error_mesg('get_soil_litter_C','The value of soil_carbon_option is invalid. This should never happen. Contact developer.',FATAL)
  end select
end subroutine get_soil_litter_C

! ============================================================================
! Nitrogen uptake from the rhizosphere by roots (active transport across root-soil interface)
! Mineral nitrogen is taken up from the rhizosphere only
subroutine root_N_uptake(soil,vegn,N_uptake_cohorts,dt,update_pools)
  real,intent(out),dimension(:)::N_uptake_cohorts  ! Units of per individual
  type(vegn_tile_type),intent(in)::vegn
  type(soil_tile_type),intent(inout)::soil
  real,intent(in)::dt
  logical, intent(in) :: update_pools

  real :: nitrate_uptake, ammonium_uptake, ammonium_concentration, nitrate_concentration
  real :: rhiz_frac(num_l)
  real,dimension(num_l) :: profile, vegn_uptake_term, total_ammonium_uptake,total_nitrate_uptake
  real::cohort_root_biomass(num_l,vegn%n_cohorts),total_root_biomass(num_l)
  integer :: k,i

  call rhizosphere_frac(vegn, rhiz_frac)

  do i=1,vegn%n_cohorts
    call cohort_root_litter_profile (vegn%cohorts(i), dz(1:num_l), profile )
    cohort_root_biomass(:,i) =  vegn%cohorts(i)%br*vegn%cohorts(i)%nindivs*profile
  enddo

  total_root_biomass = sum(cohort_root_biomass,dim=2)

  N_uptake_cohorts=0.0
  total_ammonium_uptake=0.0
  total_nitrate_uptake=0.0

  ! If there is no root biomass, skip the rest since there's no uptake in that case
  if(sum(total_root_biomass)==0.0) return

  do k=1,num_l
    ammonium_concentration=soil%soil_organic_matter(k)%ammonium*rhiz_frac(k)/dz(k)
    nitrate_concentration=soil%soil_organic_matter(k)%nitrate*rhiz_frac(k)/dz(k)
    do i=1,vegn%n_cohorts
      ammonium_uptake = ammonium_concentration/(ammonium_concentration+spdata(vegn%cohorts(i)%species)%k_ammonium_root_uptake)*&
                          spdata(vegn%cohorts(i)%species)%root_NH4_uptake_rate*dt*dz(k)*cohort_root_biomass(k,i)/total_root_biomass(k)
      nitrate_uptake = nitrate_concentration/(nitrate_concentration+spdata(vegn%cohorts(i)%species)%k_nitrate_root_uptake)*&
                          spdata(vegn%cohorts(i)%species)%root_NO3_uptake_rate*dt*dz(k)*cohort_root_biomass(k,i)/total_root_biomass(k)
      N_uptake_cohorts(i)=N_uptake_cohorts(i)+(ammonium_uptake+nitrate_uptake)/vegn%cohorts(i)%nindivs
      total_ammonium_uptake(k)=total_ammonium_uptake(k)+ammonium_uptake
      total_nitrate_uptake(k)=total_nitrate_uptake(k)+nitrate_uptake
    enddo
  enddo
  if(any(total_ammonium_uptake(:)>soil%soil_organic_matter(:)%ammonium)) then
     __DEBUG4__(rhiz_frac,ammonium_concentration,total_ammonium_uptake,soil%soil_organic_matter(:)%ammonium)
  endif
  if(any(total_nitrate_uptake(:)>soil%soil_organic_matter(:)%nitrate)) then
     __DEBUG3__(nitrate_concentration,total_nitrate_uptake,soil%soil_organic_matter(:)%nitrate)
  endif
  if (update_pools) then
     soil%soil_organic_matter(:)%ammonium=soil%soil_organic_matter(:)%ammonium-total_ammonium_uptake(:)
     soil%soil_organic_matter(:)%nitrate=soil%soil_organic_matter(:)%nitrate-total_nitrate_uptake(:)
  endif
end subroutine root_N_uptake


! ============================================================================
! Uptake of mineral N by mycorrhizal "scavengers" -- Should correspond to Arbuscular mycorrhizae
subroutine myc_scavenger_N_uptake(soil,vegn,N_uptake_cohorts,myc_efficiency,dt,update_pools)
  type(soil_tile_type),   intent(inout) :: soil
  type(vegn_tile_type),intent(in)::vegn
  real,intent(out),dimension(:) :: N_uptake_cohorts ! Units: kgN/m2 per individual
  real, intent(in) :: dt  ! dt in years
  logical, intent(in) :: update_pools
  real, intent(out) :: myc_efficiency ! units: kgN/kg myc biomass C. Should give N uptake efficiency even when myc biomass is zero

  real,dimension(num_l) :: profile, vegn_uptake_term
  real::nitrate_uptake,ammonium_uptake
  real::cohort_myc_scav_biomass(num_l,vegn%n_cohorts),total_myc_scav_biomass(num_l)
  real,dimension(num_l):: myc_biomass_tiny, N_uptake_tiny ! For calculating return on investment when myc biomass is zero
  real::litterThickness,totalC
  integer::k,i
  logical :: myc_biomass_is_zero

  associate(cc=>vegn%cohorts)
  myc_biomass_is_zero = (sum(cc(:)%myc_scavenger_biomass_C*cc(:)%nindivs)<=0)
  if(.not.myc_biomass_is_zero) then
    do i=1,vegn%n_cohorts
      call cohort_root_litter_profile (cc(i), dz(1:num_l), profile )
      cohort_myc_scav_biomass(:,i) =  cc(i)%myc_scavenger_biomass_C*cc(i)%nindivs*profile
    enddo
  else ! myc_biomass_is_zero
    ! If there is no mycorrhizal biomass, still do calculation so efficiency can be estimated
    ! This allows plants to grow some biomass from zero
    if(sum(cc(:)%bliving*cc(:)%nindivs)==0) then ! No live biomass at all. Assume uptake is zero and skip the rest.
      N_uptake_cohorts(:)=0.0
      myc_efficiency=0.0
      return
    endif

    do i=1,vegn%n_cohorts
      call cohort_root_litter_profile (cc(i), dz(1:num_l), profile )
      cohort_myc_scav_biomass(:,i) =  cc(i)%bliving*cc(i)%nindivs*0.0001*profile
    enddo
  endif


  total_myc_scav_biomass = sum(cohort_myc_scav_biomass,dim=2)
  N_uptake_cohorts=0.0


  do k=1,num_l
    ! Total uptake rate is calculated based on total mycorrhizal biomass
    call mycorrhizal_mineral_N_uptake_rate(soil%soil_organic_matter(k),total_myc_scav_biomass(k),dz(k),nitrate_uptake,ammonium_uptake)
    ammonium_uptake = min(ammonium_uptake,soil%soil_organic_matter(k)%ammonium/dt)
    nitrate_uptake = min(nitrate_uptake,soil%soil_organic_matter(k)%nitrate/dt)

    ! Uptake by each cohort is scaled by the fraction of total mycorrhizal biomass in that layer
    ! that is owned by the cohort
    do i=1,vegn%n_cohorts
      if(cc(i)%nindivs>0) then
          N_uptake_cohorts(i) = N_uptake_cohorts(i)+&
              (ammonium_uptake+nitrate_uptake)*dt*cohort_myc_scav_biomass(k,i)/(total_myc_scav_biomass(k)*cc(i)%nindivs)
      endif
    enddo

    ! These debug macros are too long for PGI! Have to break them up if we're going to leave them commented out anyway...
    !_!_DEBUG3_!_(k,ammonium_uptake/soil%soil_organic_matter(k)%ammonium,nitrate_uptake/soil%soil_organic_matter(k)%nitrate)
    !_!_DEBUG4_!_(k,total_myc_scav_biomass(k)/dz(k),soil%soil_organic_matter(k)%ammonium/dz(k),soil%soil_organic_matter(k)%nitrate/dz(k))

    if(ammonium_uptake*dt>soil%soil_organic_matter(k)%ammonium) then
      __DEBUG2__(ammonium_uptake*dt,soil%soil_organic_matter(k)%ammonium)
    endif
    if(nitrate_uptake*dt>soil%soil_organic_matter(k)%nitrate) then
      __DEBUG2__(nitrate_uptake*dt,soil%soil_organic_matter(k)%nitrate)
    endif

    if (update_pools .and. .not. myc_biomass_is_zero) then  ! If myc biomass is 0, then no actual N is taken up. It is just being used for efficiency calculation
       soil%soil_organic_matter(k)%ammonium=soil%soil_organic_matter(k)%ammonium-ammonium_uptake*dt
       soil%soil_organic_matter(k)%nitrate=soil%soil_organic_matter(k)%nitrate-nitrate_uptake*dt
    endif
  enddo

  ! Mycorrhizae should have access to litter layer too
  ! Assuming volumetric concentration in litter layer is the same as top soil layer
  do k = 1,N_LITTER_POOLS
    call poolTotals(soil%litter(k),totalCarbon=totalC)
    litterThickness=max(totalC/litterDensity,1e-2)
     call mycorrhizal_mineral_N_uptake_rate(soil%litter(k),total_myc_scav_biomass(1)/dz(1)*litterThickness,litterThickness,&
             nitrate_uptake, ammonium_uptake)
     ammonium_uptake = min(ammonium_uptake,soil%litter(k)%ammonium/dt)
     nitrate_uptake  = min(nitrate_uptake,soil%litter(k)%nitrate/dt)

     do i=1,vegn%n_cohorts
       if(cc(i)%nindivs>0) &
           N_uptake_cohorts(i)=N_uptake_cohorts(i)+(ammonium_uptake+nitrate_uptake)*dt*cohort_myc_scav_biomass(1,i)/(total_myc_scav_biomass(1)*cc(i)%nindivs)
     enddo

     if (update_pools .and. .not. myc_biomass_is_zero) then
        soil%litter(k)%ammonium=soil%litter(k)%ammonium-ammonium_uptake*dt
        soil%litter(k)%nitrate=soil%litter(k)%nitrate-nitrate_uptake*dt
     endif
  enddo

  myc_efficiency = sum(N_uptake_cohorts*cc(:)%nindivs)/sum(total_myc_scav_biomass)
  ! If myc biomass was zero, then the calculation used a "virtual" biomass and there is no real uptake
  if (myc_biomass_is_zero) then
     N_uptake_cohorts(:) = 0.0
  endif
  end associate ! cc
end subroutine myc_scavenger_N_uptake


! ============================================================================
! Uptake of mineral N by mycorrhizal "miners" -- Should correspond to Ecto mycorrhizae
subroutine myc_miner_N_uptake(soil,vegn,N_uptake_cohorts,C_uptake_cohorts,total_CO2prod,myc_efficiency,dt,update_pools)
  type(soil_tile_type), intent(inout) :: soil
  type(vegn_tile_type),intent(in)::vegn
  real,    intent(out) :: N_uptake_cohorts(:), C_uptake_cohorts(:)  ! Units kg/m2 of per individual
  real,    intent(out) :: total_CO2prod ! Units of kgC/m2 (not per individual)
  real,    intent(in)  :: dt  ! dt in years
  logical, intent(in)  :: update_pools
  real, intent(out)     :: myc_efficiency  ! units: kgN/kg myc biomass C. Should give N uptake efficiency even when myc biomass is zero

  real, dimension(num_l) :: profile, vegn_uptake_term
  real :: N_uptake,C_uptake,CO2prod
  real, dimension(num_l) :: T, theta, air_filled_porosity
  real::cohort_myc_mine_biomass(num_l,vegn%n_cohorts),total_myc_mine_biomass(num_l)
  real::litterThickness,totalC,cohort_frac
  integer :: k,i
  logical :: myc_biomass_is_zero

  associate(cc=>vegn%cohorts)
  if(sum(cc(:)%myc_miner_biomass_C*cc(:)%nindivs)>0) then
    myc_biomass_is_zero = .FALSE.
    do i=1,vegn%n_cohorts
      call cohort_root_litter_profile (cc(i), dz(1:num_l), profile )
      cohort_myc_mine_biomass(:,i) =  cc(i)%myc_miner_biomass_C*cc(i)%nindivs*profile
    enddo
  else
    ! If there is no mycorrhizal biomass, still do calculation so efficiency can be estimated
    ! This allows plants to grow some biomass from zero
    if(sum(cc(:)%bliving*cc(:)%nindivs)==0) then ! No live biomass at all. Assume uptake is zero and skip the rest.
      N_uptake_cohorts(:)=0.0
      C_uptake_cohorts(:)=0.0
      total_CO2prod=0.0
      myc_efficiency=0.0
      return
    endif

    myc_biomass_is_zero = .TRUE.
    do i=1,vegn%n_cohorts
      call cohort_root_litter_profile (cc(i), dz(1:num_l), profile)
      cohort_myc_mine_biomass(:,i) =  cc(i)%bliving*cc(i)%nindivs*0.0001*profile
    enddo
  endif

  T = soil%T(:)
  theta = max(min(soil_theta(soil),1.0),0.0)
  air_filled_porosity=max(min(1.0-theta-soil_ice_porosity(soil),1.0),0.0)

  total_myc_mine_biomass = sum(cohort_myc_mine_biomass,dim=2)

  N_uptake_cohorts=0.0
  C_uptake_cohorts=0.0
  total_CO2prod=0.0

  do k=1,num_l
     call mycorrhizal_decomposition(soil%soil_organic_matter(k),total_myc_mine_biomass(k),&
          T(k),theta(k),air_filled_porosity(k),N_uptake,C_uptake,CO2prod,dt,&
          update_pools .and. .not. myc_biomass_is_zero)
     total_CO2prod=total_CO2prod+CO2prod

     do i=1,vegn%n_cohorts
       cohort_frac = 0.0
       if(cc(i)%nindivs>0) &
          cohort_frac = cohort_myc_mine_biomass(k,i)/total_myc_mine_biomass(k)/cc(i)%nindivs
       N_uptake_cohorts(i)=N_uptake_cohorts(i)+N_uptake*cohort_frac
       C_uptake_cohorts(i)=C_uptake_cohorts(i)+C_uptake*cohort_frac
     enddo
  enddo

  do k = 1, N_LITTER_POOLS
     call poolTotals(soil%litter(k),totalCarbon=totalC)
     litterThickness=max(totalC/litterDensity,1e-2)
     call mycorrhizal_decomposition(soil%litter(k),total_myc_mine_biomass(1)/dz(1)*litterThickness,&
          T(1),theta(1),air_filled_porosity(1),N_uptake,C_uptake,CO2prod,dt,&
          update_pools .and. .not. myc_biomass_is_zero)
     total_CO2prod  = total_CO2prod + CO2prod

     do i=1,vegn%n_cohorts
       cohort_frac = 0.0
       if(cc(i)%nindivs>0) &
           cohort_frac = cohort_myc_mine_biomass(k,i)/total_myc_mine_biomass(k)/cc(i)%nindivs
       N_uptake_cohorts(i)=N_uptake_cohorts(i)+N_uptake*cohort_frac
       C_uptake_cohorts(i)=C_uptake_cohorts(i)+C_uptake*cohort_frac
     enddo
  enddo

  myc_efficiency = sum(N_uptake_cohorts*cc(:)%nindivs)/sum(total_myc_mine_biomass)

  ! If myc biomass was zero, then the calculation used a tiny "virtual" biomass and there is no real uptake or CO2 production
  if (myc_biomass_is_zero) then
    N_uptake_cohorts(:)=0.0
    C_uptake_cohorts(:)=0.0
    total_CO2prod      =0.0
  endif
  end associate ! cc
end subroutine myc_miner_N_uptake


! ============================================================================
subroutine redistribute_peat_carbon(soil)
    type(soil_tile_type), intent(inout) :: soil

    integer :: nn
    real :: layer_total_C,layer_total_C_2,layer_max_C,layer_extra_C,fraction_to_remove
    real :: total_C_before,total_C_after
    real :: leaflitter_total_C, woodlitter_total_C

    !For conservation check.
    total_C_before=0.0
    do nn=1,num_l
    call poolTotals(soil%soil_organic_matter(num_l),layer_total_C)
    total_C_before=total_C_before+layer_total_C
    enddo

    call poolTotals(soil%litter(LEAF),totalCarbon=leaflitter_total_C)
    call poolTotals(soil%litter(CWOOD),totalCarbon=woodlitter_total_C)
    layer_total_C=leaflitter_total_C+woodlitter_total_C

    layer_max_C=max_litter_thickness*max_soil_C_density
    layer_extra_C = layer_total_C-layer_max_C
    if(layer_extra_C>0) then
        fraction_to_remove=1.0-layer_max_C/layer_total_C
        call transfer_pool_fraction(soil%litter(LEAF),soil%soil_organic_matter(1),fraction_to_remove)
        call transfer_pool_fraction(soil%litter(CWOOD),soil%soil_organic_matter(1),fraction_to_remove)
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
    do nn=1,num_l
    call poolTotals(soil%soil_organic_matter(num_l),layer_total_C)
    total_C_after=total_C_after+layer_total_C
    enddo

    if (abs(total_C_before-total_C_after)>1e-10) then
            print *,'Carbon before:',total_C_before
            print *,'Carbon after:',total_C_after
            call error_mesg('redistribute_peat_carbon','Carbon not conserved after downward move',FATAL)
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
  logical, intent(in) :: land_mask(:)    ! land mask
  real,    intent(in) :: lonb(:,:), latb(:,:) ! boundaries of the grid cells
  real,    pointer :: soil_frac (:,:) ! output: map of soil fractional coverage
  integer :: k

  if (.not. do_hillslope_model) then
     allocate( soil_frac(size(land_mask(:)),n_dim_soil_types))
     allocate( soil_tags(size(land_mask(:)),n_dim_soil_types))
     do k = 1, n_dim_soil_types
        soil_tags(:, k) = k
     end do
  else
     allocate( soil_frac(size(land_mask(:)), &
               n_dim_soil_types * num_vertclusters * max_num_topo_hlsps))
     allocate( soil_tags(size(land_mask(:)), &
               n_dim_soil_types * num_vertclusters * max_num_topo_hlsps))
     do k = 1,size(soil_tags,2)
        if (mod(k, n_dim_soil_types) == 0) then
           soil_tags(:,k) = n_dim_soil_types
        else
           soil_tags(:,k) = mod(k, n_dim_soil_types)
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
  integer, pointer :: soiltags(:,:)

  if (size(soiltags,1) /= size(soil_tags,1) .or. size(soiltags,2) /= size(soil_tags,2)) &
     call error_mesg(module_name,'Wrong dimension size in "soil_tile_mod:retrieve_soil_tags.',FATAL)

  soiltags(:,:) = soil_tags(:,:)

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

!/* #define DEFINE_SOIL_ACCESSOR_0D(xtype,x) subroutine soil_ ## x ## _ptr(t,p);\ */
!/* type(land_tile_type),pointer::t;xtype,pointer::p;p=>NULL();\ */
!/* if(associated(t))then;if(associated(t%soil))p=>t%soil%x;endif;\ */
!/* end subroutine */
!/* #define DEFINE_SOIL_ACCESSOR_1D(xtype,x) subroutine soil_ ## x ## _ptr(t,i,p);\ */
!/* type(land_tile_type),pointer::t;integer,intent(in)::i;xtype,pointer::p;p=>NULL();\ */
!/* if(associated(t))then;if(associated(t%soil))p=>t%soil%x(i);endif;\ */
!/* end subroutine */
!/* #define DEFINE_SOIL_COMPONENT_ACCESSOR_0D(xtype,component,x) subroutine soil_ ## x ## _ptr(t,p);\ */
!/* type(land_tile_type),pointer::t;xtype,pointer::p;p=>NULL();\ */
!/* if(associated(t))then;if(associated(t%soil))p=>t%soil%component%x;endif;\ */
!/* end subroutine */
!/* #define DEFINE_SOIL_COMPONENT_ACCESSOR_1D(xtype,component,x) subroutine soil_ ## x ## _ptr(t,i,p);\ */
!/* type(land_tile_type),pointer::t;integer,intent(in)::i;xtype,pointer::p;p=>NULL();\ */
!/* if(associated(t))then;if(associated(t%soil))p=>t%soil%component%x(i);endif;\ */
!/* end subroutine */

!DEFINE_SOIL_ACCESSOR_1D(real,T)
subroutine soil_T_ptr(t,i,p)
    type(land_tile_type),pointer::t
    integer,intent(in)::i
    real,pointer::p
    p=>NULL()
    if(associated(t))then
        if(associated(t%soil))p=>t%soil%T(i)
    endif
end subroutine

!DEFINE_SOIL_ACCESSOR_1D(real,wl)
subroutine soil_wl_ptr(t,i,p)
    type(land_tile_type),pointer::t
    integer,intent(in)::i
    real,pointer::p
    p=>NULL()
    if(associated(t))then
        if(associated(t%soil))p=>t%soil%wl(i)
    endif
end subroutine

!DEFINE_SOIL_ACCESSOR_1D(real,ws)
subroutine soil_ws_ptr(t,i,p)
    type(land_tile_type),pointer::t
    integer,intent(in)::i
    real,pointer::p;p=>NULL()
    if(associated(t))then
        if(associated(t%soil))p=>t%soil%ws(i)
    endif
end subroutine

!DEFINE_SOIL_ACCESSOR_1D(real,groundwater)
subroutine soil_groundwater_ptr(t,i,p)
    type(land_tile_type),pointer::t
    integer,intent(in)::i
    real,pointer::p
    p=>NULL()
    if(associated(t))then
        if(associated(t%soil))p=>t%soil%groundwater(i)
    endif
end subroutine

!DEFINE_SOIL_ACCESSOR_1D(real,groundwater_T)
subroutine soil_groundwater_T_ptr(t,i,p)
    type(land_tile_type),pointer::t
    integer,intent(in)::i
    real,pointer::p
    p=>NULL()
    if(associated(t))then
        if(associated(t%soil))p=>t%soil%groundwater_T(i)
    endif
end subroutine

!DEFINE_SOIL_ACCESSOR_1D(real,w_fc)
subroutine soil_w_fc_ptr(t,i,p)
    type(land_tile_type),pointer::t
    integer,intent(in)::i
    real,pointer::p
    p=>NULL()
    if(associated(t))then
        if(associated(t%soil))p=>t%soil%w_fc(i)
    endif
end subroutine

!DEFINE_SOIL_ACCESSOR_1D(real,alpha)
subroutine soil_alpha_ptr(t,i,p)
    type(land_tile_type),pointer::t
    integer,intent(in)::i
    real,pointer::p
    p=>NULL()
    if(associated(t))then
        if(associated(t%soil))p=>t%soil%alpha(i)
    endif
end subroutine

!DEFINE_SOIL_ACCESSOR_0D(real,uptake_T)
subroutine soil_uptake_T_ptr(t,p)
    type(land_tile_type),pointer::t
    real,pointer::p
    p=>NULL()
    if(associated(t))then
        if(associated(t%soil))p=>t%soil%uptake_T
    endif
end subroutine

!DEFINE_SOIL_ACCESSOR_0D(integer,tag)
subroutine soil_tag_ptr(t,p)
    type(land_tile_type),pointer::t
    integer,pointer::p
    p=>NULL()
    if(associated(t))then
        if(associated(t%soil))p=>t%soil%tag
        endif
end subroutine

!DEFINE_SOIL_ACCESSOR_1D(real,fast_soil_C)
subroutine soil_fast_soil_C_ptr(t,i,p)
    type(land_tile_type),pointer::t
    integer,intent(in)::i
    real,pointer::p
    p=>NULL()
    if(associated(t))then
        if(associated(t%soil))p=>t%soil%fast_soil_C(i)
    endif
end subroutine

!DEFINE_SOIL_ACCESSOR_1D(real,slow_soil_C)
subroutine soil_slow_soil_C_ptr(t,i,p)
    type(land_tile_type),pointer::t
    integer,intent(in)::i
    real,pointer::p
    p=>NULL()
    if(associated(t))then
        if(associated(t%soil))p=>t%soil%slow_soil_C(i)
    endif
end subroutine

!DEFINE_SOIL_ACCESSOR_1D(real,asoil_in)
subroutine soil_asoil_in_ptr(t,i,p)
    type(land_tile_type),pointer::t
    integer,intent(in)::i
    real,pointer::p
    p=>NULL()
    if(associated(t))then
        if(associated(t%soil))p=>t%soil%asoil_in(i)
    endif
end subroutine

!DEFINE_SOIL_ACCESSOR_1D(real,fsc_in)
subroutine soil_fsc_in_ptr(t,i,p)
    type(land_tile_type),pointer::t
    integer,intent(in)::i
    real,pointer::p
    p=>NULL()
    if(associated(t))then
        if(associated(t%soil))p=>t%soil%fsc_in(i)
    endif
end subroutine

!DEFINE_SOIL_ACCESSOR_1D(real,ssc_in)
subroutine soil_ssc_in_ptr(t,i,p)
    type(land_tile_type),pointer::t
    integer,intent(in)::i
    real,pointer::p
    p=>NULL()
    if(associated(t))then
        if(associated(t%soil))p=>t%soil%ssc_in(i)
    endif
end subroutine

!DEFINE_SOIL_ACCESSOR_1D(integer,is_peat)
subroutine soil_is_peat_ptr(t,i,p)
    type(land_tile_type),pointer::t
    integer,intent(in)::i
    integer,pointer::p
    p=>NULL()
    if(associated(t))then
        if(associated(t%soil))p=>t%soil%is_peat(i)
    endif
end subroutine


!DEFINE_SOIL_COMPONENT_ACCESSOR_0D(real,pars,tau_groundwater)
subroutine soil_tau_groundwater_ptr(t,p)
    type(land_tile_type),pointer::t
    real,pointer::p
    p=>NULL()
    if(associated(t))then
        if(associated(t%soil))p=>t%soil%pars%tau_groundwater
    endif
end subroutine

!DEFINE_SOIL_COMPONENT_ACCESSOR_0D(real,pars,hillslope_length)
subroutine soil_hillslope_length_ptr(t,p)
    type(land_tile_type),pointer::t
    real,pointer::p
    p=>NULL()
    if(associated(t))then
        if(associated(t%soil))p=>t%soil%pars%hillslope_length
    endif
end subroutine

!DEFINE_SOIL_COMPONENT_ACCESSOR_0D(real,pars,hillslope_relief)
subroutine soil_hillslope_relief_ptr(t,p)
    type(land_tile_type),pointer::t
    real,pointer::p
    p=>NULL()
    if(associated(t))then
        if(associated(t%soil))p=>t%soil%pars%hillslope_relief
    endif
end subroutine

!DEFINE_SOIL_COMPONENT_ACCESSOR_0D(real,pars,hillslope_a)
subroutine soil_hillslope_a_ptr(t,p)
    type(land_tile_type),pointer::t
    real,pointer::p
    p=>NULL()
    if(associated(t))then
        if(associated(t%soil))p=>t%soil%pars%hillslope_a
    endif
end subroutine

!DEFINE_SOIL_COMPONENT_ACCESSOR_0D(real,pars,hillslope_n)
subroutine soil_hillslope_n_ptr(t,p)
    type(land_tile_type),pointer::t
    real,pointer::p
    p=>NULL()
    if(associated(t))then
        if(associated(t%soil))p=>t%soil%pars%hillslope_n
    endif
end subroutine

!DEFINE_SOIL_COMPONENT_ACCESSOR_0D(real,pars,hillslope_zeta_bar)
subroutine soil_hillslope_zeta_bar_ptr(t,p)
    type(land_tile_type),pointer::t
    real,pointer::p
    p=>NULL()
    if(associated(t))then
        if(associated(t%soil))p=>t%soil%pars%hillslope_zeta_bar
    endif
end subroutine

!DEFINE_SOIL_COMPONENT_ACCESSOR_0D(real,pars,soil_e_depth)
subroutine soil_soil_e_depth_ptr(t,p)
    type(land_tile_type),pointer::t
    real,pointer::p
    p=>NULL()
    if(associated(t))then
        if(associated(t%soil))p=>t%soil%pars%soil_e_depth
    endif
end subroutine

!DEFINE_SOIL_COMPONENT_ACCESSOR_0D(real,pars,zeta)
subroutine soil_zeta_ptr(t,p)
    type(land_tile_type),pointer::t
    real,pointer::p
    p=>NULL()
    if(associated(t))then
        if(associated(t%soil))p=>t%soil%pars%zeta
    endif
end subroutine

!DEFINE_SOIL_COMPONENT_ACCESSOR_0D(real,pars,tau)
subroutine soil_tau_ptr(t,p)
    type(land_tile_type),pointer::t
    real,pointer::p
    p=>NULL()
    if(associated(t))then
        if(associated(t%soil))p=>t%soil%pars%tau
    endif
end subroutine

!DEFINE_SOIL_COMPONENT_ACCESSOR_0D(real,pars,k_sat_gw)
subroutine soil_k_sat_gw_ptr(t,p)
    type(land_tile_type),pointer::t
    real,pointer::p
    p=>NULL()
    if(associated(t))then
        if(associated(t%soil))p=>t%soil%pars%k_sat_gw
    endif
end subroutine

!DEFINE_SOIL_COMPONENT_ACCESSOR_0D(real,pars,vwc_wilt)
subroutine soil_vwc_wilt_ptr(t,p)
    type(land_tile_type),pointer::t
    real,pointer::p
    p=>NULL()
    if(associated(t))then
        if(associated(t%soil))p=>t%soil%pars%vwc_wilt
    endif
end subroutine

!DEFINE_SOIL_COMPONENT_ACCESSOR_0D(real,pars,vwc_fc)
subroutine soil_vwc_fc_ptr(t,p)
    type(land_tile_type),pointer::t
    real,pointer::p
    p=>NULL()
    if(associated(t))then
        if(associated(t%soil))p=>t%soil%pars%vwc_fc
    endif
end subroutine

!DEFINE_SOIL_COMPONENT_ACCESSOR_0D(real,pars,vwc_sat)
subroutine soil_vwc_sat_ptr(t,p)
    type(land_tile_type),pointer::t
    real,pointer::p
    p=>NULL()
    if(associated(t))then
        if(associated(t%soil))p=>t%soil%pars%vwc_sat
    endif
end subroutine

!DEFINE_SOIL_COMPONENT_ACCESSOR_0D(real,pars,k_sat_ref)
subroutine soil_k_sat_ref_ptr(t,p)
    type(land_tile_type),pointer::t
    real,pointer::p
    p=>NULL()
    if(associated(t))then
        if(associated(t%soil))p=>t%soil%pars%k_sat_ref
    endif
end subroutine

!DEFINE_SOIL_COMPONENT_ACCESSOR_0D(real,pars,Qmax)
subroutine soil_Qmax_ptr(t,p)
    type(land_tile_type),pointer::t
    real,pointer::p
    p=>NULL()
    if(associated(t))then
        if(associated(t%soil))p=>t%soil%pars%Qmax
    endif
end subroutine


!DEFINE_SOIL_COMPONENT_ACCESSOR_1D(real,pars,refl_dry_dir)
subroutine soil_refl_dry_dir_ptr(t,i,p)
    type(land_tile_type),pointer::t
    integer,intent(in)::i
    real,pointer::p
    p=>NULL()
    if(associated(t))then
        if(associated(t%soil))p=>t%soil%pars%refl_dry_dir(i)
    endif
end subroutine

!DEFINE_SOIL_COMPONENT_ACCESSOR_1D(real,pars,refl_dry_dif)
subroutine soil_refl_dry_dif_ptr(t,i,p)
    type(land_tile_type),pointer::t
    integer,intent(in)::i
    real,pointer::p
    p=>NULL()
    if(associated(t))then
        if(associated(t%soil))p=>t%soil%pars%refl_dry_dif(i)
    endif
end subroutine

!DEFINE_SOIL_COMPONENT_ACCESSOR_1D(real,pars,refl_sat_dir)
subroutine soil_refl_sat_dir_ptr(t,i,p)
    type(land_tile_type),pointer::t
    integer,intent(in)::i
    real,pointer::p
    p=>NULL()
    if(associated(t))then
        if(associated(t%soil))p=>t%soil%pars%refl_sat_dir(i)
    endif
end subroutine

!DEFINE_SOIL_COMPONENT_ACCESSOR_1D(real,pars,refl_sat_dif)
subroutine soil_refl_sat_dif_ptr(t,i,p)
    type(land_tile_type),pointer::t
    integer,intent(in)::i
    real,pointer::p
    p=>NULL()
    if(associated(t))then
        if(associated(t%soil))p=>t%soil%pars%refl_sat_dif(i)
    endif
end subroutine

!DEFINE_SOIL_COMPONENT_ACCESSOR_1D(real,pars,f_iso_dry)
subroutine soil_f_iso_dry_ptr(t,i,p)
    type(land_tile_type),pointer::t
    integer,intent(in)::i
    real,pointer::p
    p=>NULL()
    if(associated(t))then
        if(associated(t%soil))p=>t%soil%pars%f_iso_dry(i)
    endif
end subroutine

!DEFINE_SOIL_COMPONENT_ACCESSOR_1D(real,pars,f_vol_dry)
subroutine soil_f_vol_dry_ptr(t,i,p)
    type(land_tile_type),pointer::t
    integer,intent(in)::i
    real,pointer::p
    p=>NULL()
    if(associated(t))then
        if(associated(t%soil))p=>t%soil%pars%f_vol_dry(i)
    endif
end subroutine

!DEFINE_SOIL_COMPONENT_ACCESSOR_1D(real,pars,f_geo_dry)
subroutine soil_f_geo_dry_ptr(t,i,p)
    type(land_tile_type),pointer::t
    integer,intent(in)::i
    real,pointer::p
    p=>NULL()
    if(associated(t))then
        if(associated(t%soil))p=>t%soil%pars%f_geo_dry(i)
    endif
end subroutine

!DEFINE_SOIL_COMPONENT_ACCESSOR_1D(real,pars,f_iso_sat)
subroutine soil_f_iso_sat_ptr(t,i,p)
    type(land_tile_type),pointer::t
    integer,intent(in)::i
    real,pointer::p
    p=>NULL()
    if(associated(t))then
        if(associated(t%soil))p=>t%soil%pars%f_iso_sat(i)
    endif
end subroutine

!DEFINE_SOIL_COMPONENT_ACCESSOR_1D(real,pars,f_vol_sat)
subroutine soil_f_vol_sat_ptr(t,i,p)
    type(land_tile_type),pointer::t
    integer,intent(in)::i
    real,pointer::p
    p=>NULL()
    if(associated(t))then
        if(associated(t%soil))p=>t%soil%pars%f_vol_sat(i)
    endif
end subroutine

!DEFINE_SOIL_COMPONENT_ACCESSOR_1D(real,pars,f_geo_sat)
subroutine soil_f_geo_sat_ptr(t,i,p)
    type(land_tile_type),pointer::t
    integer,intent(in)::i
    real,pointer::p
    p=>NULL();
    if(associated(t))then
        if(associated(t%soil))p=>t%soil%pars%f_geo_sat(i)
    endif
end subroutine


!DEFINE_SOIL_ACCESSOR_0D(real,fast_DOC_leached)
subroutine soil_fast_DOC_leached_ptr(t,p)
    type(land_tile_type),pointer::t
    real,pointer::p
    p=>NULL()
    if(associated(t))then
        if(associated(t%soil))p=>t%soil%fast_DOC_leached
    endif
end subroutine

!DEFINE_SOIL_ACCESSOR_0D(real,slow_DOC_leached)
subroutine soil_slow_DOC_leached_ptr(t,p)
    type(land_tile_type),pointer::t
    real,pointer::p
    p=>NULL()
    if(associated(t))then
        if(associated(t%soil))p=>t%soil%slow_DOC_leached
    endif
end subroutine

!DEFINE_SOIL_ACCESSOR_0D(real,deadmic_DOC_leached)
subroutine soil_deadmic_DOC_leached_ptr(t,p)
    type(land_tile_type),pointer::t
    real,pointer::p
    p=>NULL()
    if(associated(t))then
        if(associated(t%soil))p=>t%soil%deadmic_DOC_leached
    endif
end subroutine


!DEFINE_SOIL_ACCESSOR_0D(real,gross_nitrogen_flux_into_tile)
subroutine soil_gross_nitrogen_flux_into_tile_ptr(t,p)
    type(land_tile_type),pointer::t
    real,pointer::p
    p=>NULL()
    if(associated(t))then
        if(associated(t%soil))p=>t%soil%gross_nitrogen_flux_into_tile
    endif
end subroutine

!DEFINE_SOIL_ACCESSOR_0D(real,gross_nitrogen_flux_out_of_tile)
subroutine soil_gross_nitrogen_flux_out_of_tile_ptr(t,p)
    type(land_tile_type),pointer::t
    real,pointer::p
    p=>NULL()
    if(associated(t))then
        if(associated(t%soil))p=>t%soil%gross_nitrogen_flux_out_of_tile
    endif
end subroutine

! stuff below is for CORPSE
subroutine sc_soil_C_ptr(t,i,j,k,p)
  type(land_tile_type),pointer::t; integer,intent(in)::i,j,k;real,pointer::p
  p=>NULL()
  if(associated(t)) then
     if(associated(t%soil))p=>t%soil%soil_organic_matter(i)%litterCohorts(j)%litterC(k)
  endif
end subroutine

subroutine sc_negative_litter_C_ptr(t,i,p)
  type(land_tile_type),pointer::t; integer,intent(in)::i; real,pointer::p; p=>NULL()
  if(associated(t))then
     if(associated(t%soil))p=>t%soil%neg_litt_C(i)
  endif
end subroutine

subroutine sc_negative_litter_N_ptr(t,i,p)
  type(land_tile_type),pointer::t; integer,intent(in)::i; real,pointer::p; p=>NULL()
  if(associated(t))then
     if(associated(t%soil))p=>t%soil%neg_litt_N(i)
  endif
end subroutine

subroutine sc_soil_N_ptr(t,i,j,k,p)
  type(land_tile_type),pointer::t; integer,intent(in)::i,j,k;real,pointer::p
  p=>NULL()
  if(associated(t)) then
     if(associated(t%soil))p=>t%soil%soil_organic_matter(i)%litterCohorts(j)%litterN(k)
  endif
end subroutine

subroutine sc_protected_C_ptr(t,i,j,k,p)
  type(land_tile_type),pointer::t; integer,intent(in)::i,j,k;real,pointer::p
  p=>NULL()
  if(associated(t)) then
     if(associated(t%soil))p=>t%soil%soil_organic_matter(i)%litterCohorts(j)%protectedC(k)
  endif
end subroutine

subroutine sc_protected_N_ptr(t,i,j,k,p)
  type(land_tile_type),pointer::t; integer,intent(in)::i,j,k;real,pointer::p
  p=>NULL()
  if(associated(t)) then
     if(associated(t%soil))p=>t%soil%soil_organic_matter(i)%litterCohorts(j)%protectedN(k)
  endif
end subroutine

subroutine sc_DOC_ptr(t,i,j,p)
  type(land_tile_type),pointer::t; integer,intent(in)::i,j;real,pointer::p
  p=>NULL()
  if(associated(t)) then
     if(associated(t%soil))p=>t%soil%soil_organic_matter(i)%dissolved_carbon(j)
  endif
end subroutine

subroutine sc_DON_ptr(t,i,j,p)
  type(land_tile_type),pointer::t; integer,intent(in)::i,j;real,pointer::p
  p=>NULL()
  if(associated(t)) then
     if(associated(t%soil))p=>t%soil%soil_organic_matter(i)%dissolved_nitrogen(j)
  endif
end subroutine

subroutine sc_C_in_ptr(t,i,j,p)
  type(land_tile_type),pointer::t; integer,intent(in)::i,j;real,pointer::p
  p=>NULL()
  if(associated(t)) then
     if(associated(t%soil))p=>t%soil%soil_organic_matter(i)%C_in(j)
  endif
end subroutine

subroutine sc_N_in_ptr(t,i,j,p)
  type(land_tile_type),pointer::t; integer,intent(in)::i,j;real,pointer::p
  p=>NULL()
  if(associated(t)) then
     if(associated(t%soil))p=>t%soil%soil_organic_matter(i)%N_in(j)
  endif
end subroutine

subroutine sc_protected_C_in_ptr(t,i,j,p)
  type(land_tile_type),pointer::t; integer,intent(in)::i,j;real,pointer::p
  p=>NULL()
  if(associated(t)) then
     if(associated(t%soil))p=>t%soil%soil_organic_matter(i)%protected_C_in(j)
  endif
end subroutine

subroutine sc_protected_N_in_ptr(t,i,j,p)
  type(land_tile_type),pointer::t; integer,intent(in)::i,j;real,pointer::p
  p=>NULL()
  if(associated(t)) then
     if(associated(t%soil))p=>t%soil%soil_organic_matter(i)%protected_N_in(j)
  endif
end subroutine

subroutine sc_C_turnover_ptr(t,i,j,p)
  type(land_tile_type),pointer::t; integer,intent(in)::i,j;real,pointer::p
  p=>NULL()
  if(associated(t)) then
     if(associated(t%soil))p=>t%soil%soil_organic_matter(i)%C_turnover(j)
  endif
end subroutine

subroutine sc_protected_C_turnover_ptr(t,i,j,p)
  type(land_tile_type),pointer::t; integer,intent(in)::i,j;real,pointer::p
  p=>NULL()
  if(associated(t)) then
     if(associated(t%soil))p=>t%soil%soil_organic_matter(i)%protected_C_turnover(j)
  endif
end subroutine

subroutine sc_N_turnover_ptr(t,i,j,p)
  type(land_tile_type),pointer::t; integer,intent(in)::i,j;real,pointer::p
  p=>NULL()
  if(associated(t)) then
     if(associated(t%soil))p=>t%soil%soil_organic_matter(i)%N_turnover(j)
  endif
end subroutine

subroutine sc_protected_N_turnover_ptr(t,i,j,p)
  type(land_tile_type),pointer::t; integer,intent(in)::i,j;real,pointer::p
  p=>NULL()
  if(associated(t)) then
     if(associated(t%soil))p=>t%soil%soil_organic_matter(i)%protected_N_turnover(j)
  endif
end subroutine

subroutine sc_nitrate_ptr(t,i,p)
  type(land_tile_type),pointer::t; integer,intent(in)::i;real,pointer::p
  p=>NULL()
  if(associated(t)) then
     if(associated(t%soil))p=>t%soil%soil_organic_matter(i)%nitrate
  endif
end subroutine

subroutine sc_ammonium_ptr(t,i,p)
  type(land_tile_type),pointer::t; integer,intent(in)::i;real,pointer::p
  p=>NULL()
  if(associated(t)) then
     if(associated(t%soil))p=>t%soil%soil_organic_matter(i)%ammonium
  endif
end subroutine

subroutine sc_nitrif_ptr(t,i,p)
  type(land_tile_type),pointer::t; integer,intent(in)::i;real,pointer::p
  p=>NULL()
  if(associated(t)) then
     if(associated(t%soil))p=>t%soil%soil_organic_matter(i)%nitrif
  endif
end subroutine

subroutine sc_denitrif_ptr(t,i,p)
  type(land_tile_type),pointer::t; integer,intent(in)::i;real,pointer::p
  p=>NULL()
  if(associated(t)) then
     if(associated(t%soil))p=>t%soil%soil_organic_matter(i)%denitrif
  endif
end subroutine

!/* #define DEFINE_LITTER_COMPONENT_ACCESSOR0(xtype,x) subroutine sc_litter_ ## x ## _ptr(t,i,p);\ */
!/* type(land_tile_type),pointer::t;integer,intent(in)::i;xtype,pointer::p;p=>NULL();if(associated(t))then;if(associated(t%soil))p=>t%soil%litter(i)%x;endif;\ */
!/* end subroutine */

!DEFINE_LITTER_COMPONENT_ACCESSOR0(real,nitrate)
subroutine sc_litter_nitrate_ptr(t,i,p)
    type(land_tile_type),pointer::t
    integer,intent(in)::i
    real,pointer::p
    p=>NULL()
    if(associated(t))then
        if(associated(t%soil))p=>t%soil%litter(i)%nitrate
    endif
end subroutine

!DEFINE_LITTER_COMPONENT_ACCESSOR0(real,ammonium)
subroutine sc_litter_ammonium_ptr(t,i,p)
    type(land_tile_type),pointer::t
    integer,intent(in)::i
    real,pointer::p
    p=>NULL()
    if(associated(t))then
        if(associated(t%soil))p=>t%soil%litter(i)%ammonium
    endif
end subroutine

!DEFINE_LITTER_COMPONENT_ACCESSOR0(real,nitrif)
subroutine sc_litter_nitrif_ptr(t,i,p)
    type(land_tile_type),pointer::t
    integer,intent(in)::i
    real,pointer::p
    p=>NULL()
    if(associated(t))then
        if(associated(t%soil))p=>t%soil%litter(i)%nitrif
    endif
end subroutine

!DEFINE_LITTER_COMPONENT_ACCESSOR0(real,denitrif)
subroutine sc_litter_denitrif_ptr(t,i,p)
    type(land_tile_type),pointer::t
    integer,intent(in)::i
    real,pointer::p
    p=>NULL()
    if(associated(t))then
        if(associated(t%soil))p=>t%soil%litter(i)%denitrif
    endif
end subroutine

!/* #define DEFINE_LITTER_COMPONENT_ACCESSOR1(xtype,x) subroutine sc_litter_ ## x ## _ptr(t,i,j,p);\ */
!/* type(land_tile_type),pointer::t;integer,intent(in)::i,j;xtype,pointer::p;p=>NULL();\ */
!/* if(associated(t))then;if(associated(t%soil))p=>t%soil%litter(j)%x(i);endif;\ */
!/* end subroutine */

!DEFINE_LITTER_COMPONENT_ACCESSOR1(real,C_in)
subroutine sc_litter_C_in_ptr(t,i,j,p)
    type(land_tile_type),pointer::t
    integer,intent(in)::i,j
    real,pointer::p
    p=>NULL()
    if(associated(t)) then
        if(associated(t%soil))p=>t%soil%litter(j)%C_in(i)
    endif
end subroutine

!DEFINE_LITTER_COMPONENT_ACCESSOR1(real,N_in)
subroutine sc_litter_N_in_ptr(t,i,j,p)
    type(land_tile_type),pointer::t
    integer,intent(in)::i,j
    real,pointer::p
    p=>NULL()
    if(associated(t)) then
        if(associated(t%soil))p=>t%soil%litter(j)%N_in(i)
    endif
end subroutine

!DEFINE_LITTER_COMPONENT_ACCESSOR1(real,C_turnover)
subroutine sc_litter_C_turnover_ptr(t,i,j,p)
    type(land_tile_type),pointer::t
    integer,intent(in)::i,j
    real,pointer::p
    p=>NULL()
    if(associated(t)) then
        if(associated(t%soil))p=>t%soil%litter(j)%C_turnover(i)
    endif
end subroutine

!DEFINE_LITTER_COMPONENT_ACCESSOR1(real,N_turnover)
subroutine sc_litter_N_turnover_ptr(t,i,j,p)
    type(land_tile_type),pointer::t
    integer,intent(in)::i,j
    real,pointer::p
    p=>NULL()
    if(associated(t)) then
        if(associated(t%soil))p=>t%soil%litter(j)%N_turnover(i)
    endif
end subroutine

!DEFINE_LITTER_COMPONENT_ACCESSOR1(real,dissolved_carbon)
subroutine sc_litter_dissolved_carbon_ptr(t,i,j,p)
    type(land_tile_type),pointer::t
    integer,intent(in)::i,j
    real,pointer::p
    p=>NULL()
    if(associated(t)) then
        if(associated(t%soil))p=>t%soil%litter(j)%dissolved_carbon(i)
    endif
end subroutine

!DEFINE_LITTER_COMPONENT_ACCESSOR1(real,dissolved_nitrogen)
subroutine sc_litter_dissolved_nitrogen_ptr(t,i,j,p)
    type(land_tile_type),pointer::t
    integer,intent(in)::i,j
    real,pointer::p
    p=>NULL()
    if(associated(t)) then
        if(associated(t%soil))p=>t%soil%litter(j)%dissolved_nitrogen(i)
    endif
end subroutine

!/* #define DEFINE_LITTER_COHORT_COMPONENT_ACCESSOR0(xtype,x) subroutine sc_litter_ ## x ## _ptr(t,i,j,p);\ */
!/* type(land_tile_type),pointer::t;integer,intent(in)::i,j;xtype,pointer::p;p=>NULL();if(associated(t))then;\ */
!/* if(associated(t%soil))p=>t%soil%litter(j)%litterCohorts(i)%x;endif;\ */
!/* end subroutine */

!DEFINE_LITTER_COHORT_COMPONENT_ACCESSOR0(real,livingMicrobeC)
subroutine sc_litter_livingMicrobeC_ptr(t,i,j,p)
    type(land_tile_type),pointer::t
    integer,intent(in)::i,j;xtype,pointer::p
    p=>NULL()
    if(associated(t)) then
        if(associated(t%soil))p=>t%soil%litter(j)%litterCohorts(i)%livingMicrobeC
    endif
end subroutine

!DEFINE_LITTER_COHORT_COMPONENT_ACCESSOR0(real,livingMicrobeN)
subroutine sc_litter_livingMicrobeN_ptr(t,i,j,p)
    type(land_tile_type),pointer::t
    integer,intent(in)::i,j;xtype,pointer::p
    p=>NULL()
    if(associated(t)) then
        if(associated(t%soil))p=>t%soil%litter(j)%litterCohorts(i)%livingMicrobeN
    endif
end subroutine

!DEFINE_LITTER_COHORT_COMPONENT_ACCESSOR0(real,CO2)
subroutine sc_litter_CO2_ptr(t,i,j,p)
    type(land_tile_type),pointer::t
    integer,intent(in)::i,j;xtype,pointer::p
    p=>NULL()
    if(associated(t)) then
        if(associated(t%soil))p=>t%soil%litter(j)%litterCohorts(i)%CO2
    endif
end subroutine

!DEFINE_LITTER_COHORT_COMPONENT_ACCESSOR0(real,Rtot)
subroutine sc_litter_Rtot_ptr(t,i,j,p)
    type(land_tile_type),pointer::t
    integer,intent(in)::i,j;xtype,pointer::p
    p=>NULL()
    if(associated(t)) then
        if(associated(t%soil))p=>t%soil%litter(j)%litterCohorts(i)%Rtot
    endif
end subroutine

!DEFINE_LITTER_COHORT_COMPONENT_ACCESSOR0(real,originalLitterC)
subroutine sc_litter_originalLitterC_ptr(t,i,j,p)
    type(land_tile_type),pointer::t
    integer,intent(in)::i,j;xtype,pointer::p
    p=>NULL()
    if(associated(t)) then
        if(associated(t%soil))p=>t%soil%litter(j)%litterCohorts(i)%originalLitterC
    endif
end subroutine

!DEFINE_LITTER_COHORT_COMPONENT_ACCESSOR0(real,originalLitterN)
subroutine sc_litter_originalLitterN_ptr(t,i,j,p)
    type(land_tile_type),pointer::t
    integer,intent(in)::i,j;xtype,pointer::p
    p=>NULL()
    if(associated(t)) then
        if(associated(t%soil))p=>t%soil%litter(j)%litterCohorts(i)%originalLitterN
    endif
end subroutine

#define DEFINE_LITTER_COHORT_COMPONENT_ACCESSOR1(xtype,x) subroutine sc_litter_ ## x ## _ptr(t,i,j,k,p);\
type(land_tile_type),pointer::t;integer,intent(in)::i,j,k;xtype,pointer::p;p=>NULL();if(associated(t))then;\
if(associated(t%soil))p=>t%soil%litter(k)%litterCohorts(i)%x(j);endif;\
end subroutine
DEFINE_LITTER_COHORT_COMPONENT_ACCESSOR1(real,litterC)
DEFINE_LITTER_COHORT_COMPONENT_ACCESSOR1(real,protectedC)
DEFINE_LITTER_COHORT_COMPONENT_ACCESSOR1(real,litterN)
DEFINE_LITTER_COHORT_COMPONENT_ACCESSOR1(real,protectedN)

#define DEFINE_SOIL_LAYER_COHORT_COMPONENT_ACCESSOR1(xtype,x) subroutine sc_ ## x ## _ptr(t,i,j,p);\
type(land_tile_type),pointer::t;xtype,pointer::p;integer,intent(in)::i,j;p=>NULL();if(associated(t))then;\
if(associated(t%soil))p=>t%soil%soil_organic_matter(i)%litterCohorts(j)%x;endif;\
end subroutine
DEFINE_SOIL_LAYER_COHORT_COMPONENT_ACCESSOR1(real,livingMicrobeC)
DEFINE_SOIL_LAYER_COHORT_COMPONENT_ACCESSOR1(real,Rtot)
DEFINE_SOIL_LAYER_COHORT_COMPONENT_ACCESSOR1(real,CO2)
DEFINE_SOIL_LAYER_COHORT_COMPONENT_ACCESSOR1(real,originalLitterC)

DEFINE_SOIL_LAYER_COHORT_COMPONENT_ACCESSOR1(real,livingMicrobeN)
DEFINE_SOIL_LAYER_COHORT_COMPONENT_ACCESSOR1(real,originalLitterN)

end module soil_mod
