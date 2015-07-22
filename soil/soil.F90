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
use time_manager_mod,   only: time_type, time_type_to_real
use diag_manager_mod,   only: diag_axis_init
use constants_mod,      only: tfreeze, hlv, hlf, dens_h2o

use land_constants_mod, only : NBANDS, BAND_VIS, BAND_NIR
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
     num_zeta_pts, num_tau_pts, soil_type_file, &
     log_rho_table, log_zeta_s, log_tau, gw_scale_perm, &
     z_ref, use_tau_fix, &
     soil_tile_stock_pe, initval, comp

use land_tile_mod, only : land_tile_type, land_tile_enum_type, &
     first_elmt, tail_elmt, next_elmt, prev_elmt, current_tile, get_elmt_indices, &
     operator(/=)
use land_utils_mod, only : put_to_tiles_r0d_fptr, put_to_tiles_r1d_fptr
use land_tile_diag_mod, only : diag_buff_type, &
     register_tiled_static_field, register_tiled_diag_field, &
     send_tile_data, send_tile_data_r0d_fptr, send_tile_data_r1d_fptr, &
     send_tile_data_i0d_fptr, &
     add_tiled_diag_field_alias, add_tiled_static_field_alias
use land_data_mod, only : land_state_type, lnd, land_time
use land_io_mod, only : read_field
use land_tile_io_mod, only : create_tile_out_file, write_tile_data_r0d_fptr,& 
     write_tile_data_r1d_fptr, read_tile_data_r0d_fptr, read_tile_data_r1d_fptr,&
     print_netcdf_error, get_input_restart_name, sync_nc_files
use nf_utils_mod, only : nfu_def_dim, nfu_put_att, nfu_inq_var
use vegn_cohort_mod, only : vegn_cohort_type, & 
     cohort_uptake_profile, cohort_root_properties 

use vegn_tile_mod, only : vegn_tile_type, vegn_tile_bwood
use land_debug_mod, only : is_watch_point, is_watch_cell, get_current_point, &
     set_current_point, check_var_range, check_conservation
use uptake_mod, only : uptake_init, &
     uptake_option, UPTAKE_LINEAR, UPTAKE_DARCY2D, UPTAKE_DARCY2D_LIN, &
     darcy2d_uptake, darcy2d_uptake_solver

use hillslope_mod, only : do_hillslope_model, max_num_topo_hlsps, &
     num_vertclusters, hlsp_coldfracs, use_geohydrodata, & !pond, &
     horiz_wt_depth_to_init, calculate_wt_init, simple_inundation
use land_io_mod, only : &
     init_cover_field
use soil_tile_mod, only : n_dim_soil_types, soil_to_use, &
     soil_index_constant, input_cover_types
use hillslope_hydrology_mod, only: hlsp_hydro_lev_init, hlsp_hydrology_2, &
     stiff_explicit_gwupdate

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
logical :: allow_negative_rie   = .false.
logical :: baseflow_where_frozen = .false.
logical :: write_when_flagged   = .false.
logical :: bypass_richards_when_stiff = .true.
logical :: corrected_lm2_gw     = .true.
logical :: use_fringe           = .false.
logical :: push_down_sfc_excess = .true.
logical :: lrunf_from_div       = .true.
logical :: cold_infilt          = .false.
logical :: bottom_up_cold_infilt= .false.
logical :: use_depth_to_wt_4    = .false.
logical :: always_use_bsw       = .false.
logical :: split_delta_t        = .false.
logical :: harmonic_mean_K      = .false.
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
logical :: allow_neg_wl         = .false.   ! Warn rather than abort if wl < 0
real    :: zeta_bar_override    = -1.
real    :: soil_root_separation_coef = 0.    ! loss of contact with drying
real    :: cold_depth           = 0.
real    :: dW_W_max             = 1.
real    :: Wl_min               = -1.e20
real    :: bwood_macinf         = -1.
integer :: max_iter_trans = 100 ! max number of iterations for psi_crown_min
integer :: layer_for_gw_switch = 1000000 ! to accelerate permafrost gw shutoff
real    :: eps_trans     = 1.e-7 ! convergence crit for psi_crown_min
logical :: supercooled_rnu = .true. ! Excess ice converted to supercooled water for runoff.
real    :: wet_depth = 0.6 ! [m] water table depth threshold for diagnosing wetland fraction

namelist /soil_nml/ lm2, use_E_min, use_E_max,           &
                    init_temp,      &
                    init_w,   init_wtdep,    &
                    init_groundwater, lrunf_ie_min, lrunf_ie_tol, &
                    cpw, clw, csw, &
                    albedo_to_use, &
                    allow_negative_rie, &
                    baseflow_where_frozen, &
                    write_when_flagged, &
                    bypass_richards_when_stiff, corrected_lm2_gw, &
                    use_fringe, &
                    cold_infilt, bottom_up_cold_infilt, &
                    use_depth_to_wt_4, always_use_bsw, split_delta_t, &
                    harmonic_mean_K, &
                    push_down_sfc_excess, lrunf_from_div, &
                    active_layer_drainage_acceleration, hlf_factor, &
                    gw_flux_max, aquifer_heat_cap, use_tridiag_foradvec, &
                    horiz_init_wt, use_coldstart_wtt_data, coldstart_datafile, &
                    allow_neg_wl, &
                    zeta_bar_override, &
                    soil_root_separation_coef, &
                    cold_depth, dW_W_max, Wl_min, &
                    bwood_macinf, &
                    max_iter_trans, layer_for_gw_switch, eps_trans, &
                    supercooled_rnu, wet_depth, &
                    write_soil_carbon_restart
!---- end of namelist --------------------------------------------------------

logical         :: module_is_initialized =.FALSE.
logical         :: use_brdf = .false.
real            :: delta_time
logical         :: use_single_geo
real            :: Eg_min

integer         :: gw_option = -1 

! ---- diagnostic field IDs
integer :: id_fast_soil_C, id_slow_soil_C, id_fsc, id_ssc, &
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
    id_excess, id_deficit, id_deficit_2, id_deficit_3, id_deficit_4, id_zeta, id_tau, &
    id_psi_bot, id_sat_frac, id_stor_frac, id_sat_depth, id_sat_dept2, &
    id_cf_1, id_cf_3, id_wt_1, id_wt_2, id_wt_2a, id_wt_2b, id_wt_3, id_wt2_3, id_wt_4, &
    id_div_bf, id_div_if, id_div_al, &
    id_z_cap, id_active_layer, id_surface_water, id_inun_frac, id_rsn_frac, id_flow, id_reflux, &
    id_wet_frac, &
    id_macro_infilt

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

  call read_soil_data_namelist(use_single_geo,gw_option)

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

  if (use_E_min) then
     Eg_min = 0.
  else
     Eg_min = -HUGE(Eg_min)
  endif
end subroutine read_soil_namelist


! ============================================================================
! initialize soil model
subroutine soil_init ( id_lon, id_lat, id_band, id_zfull )
  integer, intent(in)  :: id_lon   ! ID of land longitude (X) axis  
  integer, intent(in)  :: id_lat   ! ID of land latitude (Y) axis
  integer, intent(in)  :: id_band  ! ID of spectral band axis
  integer, intent(out) :: id_zfull ! ID of vertical soil axis

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
  logical :: restart_exists

  module_is_initialized = .TRUE.
  delta_time = time_type_to_real(lnd%dt_fast)

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
     ! the albedo doesn't depend on soil wetness
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
     ! the albedo doesn't depend on soil wetness
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
              psi = zfull - local_wt_depth
           else
              psi = zfull - init_wtdep
           end if
        else
           if (wetmask(li, lj) > 0.5) then ! wet point
              drypoint = .false.
           else
              drypoint = .true.
           end if
           if (gw_option == GW_TILED) then
              call horiz_wt_depth_to_init(tile%soil, prev_elmt(ce), local_wt_depth, dry=drypoint)
              psi = zfull - local_wt_depth
           else if (drypoint) then
              psi = zfull - tile%soil%pars%hillslope_relief*tile%soil%pars%hillslope_zeta_bar
           else
              psi = zfull - init_wtdep
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
        call init_soil_twc(tile%soil, ref_soil_t(li, lj), mwc)
     end if
  end do

  if (use_coldstart_wtt_data) then
     deallocate(ref_soil_t, wetmask)
  end if

  call get_input_restart_name('INPUT/soil.res.nc',restart_exists,restart_file_name)
  if (restart_exists) then
     call error_mesg('soil_init',&
          'reading NetCDF restart "'//trim(restart_file_name)//'"',&
          NOTE)
     __NF_ASRT__(nf_open(restart_file_name,NF_NOWRITE,unit))
     call read_tile_data_r1d_fptr(unit, 'temp'         , soil_T_ptr  )
     call read_tile_data_r1d_fptr(unit, 'wl'           , soil_wl_ptr )
     call read_tile_data_r1d_fptr(unit, 'ws'           , soil_ws_ptr )
     call read_tile_data_r1d_fptr(unit, 'groundwater'  , soil_groundwater_ptr )
     call read_tile_data_r1d_fptr(unit, 'groundwater_T', soil_groundwater_T_ptr)
     if(nfu_inq_var(unit, 'uptake_T')==NF_NOERR) &
          call read_tile_data_r0d_fptr(unit, 'uptake_T', soil_uptake_T_ptr)
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
          
     __NF_ASRT__(nf_close(unit))     
  else
     call error_mesg('soil_init', 'cold-starting soil', NOTE)
  endif

  ! read soil carbon restart, if present
  call get_input_restart_name('INPUT/soil_carbon.res.nc',restart_exists,restart_file_name)
  if (restart_exists) then
     __NF_ASRT__(nf_open(restart_file_name,NF_NOWRITE,unit))
     call error_mesg('veg_data_init','reading soil_carbon restart',NOTE)
     call read_tile_data_r1d_fptr(unit,'asoil_in',soil_asoil_in_ptr)
     call read_tile_data_r1d_fptr(unit,'fsc_in',soil_fsc_in_ptr)
     call read_tile_data_r1d_fptr(unit,'ssc_in',soil_ssc_in_ptr)
     __NF_ASRT__(nf_close(unit))     
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

  ! define vertical axis and its' edges
  id_zhalf = diag_axis_init ( &
       'zhalf_soil', zhalf(1:num_l+1), 'meters', 'z', 'half level',  -1, set_name='soil' )
  id_zfull = diag_axis_init ( &
       'zfull_soil', zfull(1:num_l),   'meters', 'z', 'full level',  -1, set_name='soil', &
       edges=id_zhalf )

  ! define array of axis indices
  axes = (/ id_lon, id_lat, id_zfull /)

  ! define diagnostic fields
  id_fast_soil_C = register_tiled_diag_field ( module_name, 'fast_soil_C', axes,  &
       land_time, 'fast soil carbon', 'kg C/m3', missing_value=-100.0 )
  id_slow_soil_C = register_tiled_diag_field ( module_name, 'slow_soil_C', axes,  &
       land_time, 'slow soil carbon', 'kg C/m3', missing_value=-100.0 )
  id_fsc = register_tiled_diag_field ( module_name, 'fsc', axes(1:2),  &
       land_time, 'total fast soil carbon', 'kg C/m2', missing_value=-100.0 )
  id_ssc = register_tiled_diag_field ( module_name, 'ssc', axes(1:2),  &
       land_time, 'total slow soil carbon', 'kg C/m2', missing_value=-100.0 )
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
       land_time, 'soil_cf_1', 'm',  missing_value=-100.0 )
  id_wt_1 = register_tiled_diag_field ( module_name, 'soil_wt_1', axes(1:2), &
       land_time, 'soil_wt_1', 'm',  missing_value=-100.0 )
  id_wt_2 = register_tiled_diag_field ( module_name, 'soil_wt_2', axes(1:2), &
       land_time, 'soil_wt_2', 'm',  missing_value=-100.0 )
  id_wt_2a = register_tiled_diag_field ( module_name, 'soil_wt_2a', axes(1:2), &
       land_time, 'soil_wt_2a', 'm',  missing_value=-100.0 )
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
  endif
  call sync_nc_files(unit)
        
  ! write out fields
  call write_tile_data_r1d_fptr(unit,'temp'         ,soil_T_ptr   ,'zfull','soil temperature','degrees_K')
  call write_tile_data_r1d_fptr(unit,'wl'           ,soil_wl_ptr  ,'zfull','liquid water content','kg/m2')
  call write_tile_data_r1d_fptr(unit,'ws'           ,soil_ws_ptr  ,'zfull','solid water content','kg/m2')
  call write_tile_data_r1d_fptr(unit,'groundwater'  ,soil_groundwater_ptr  ,'zfull')
  call write_tile_data_r1d_fptr(unit,'groundwater_T',soil_groundwater_T_ptr ,'zfull')
  call write_tile_data_r0d_fptr(unit,'uptake_T',     soil_uptake_T_ptr, 'temperature of transpiring water', 'degrees_K')
  call write_tile_data_r1d_fptr(unit,'fsc',          soil_fast_soil_C_ptr,'zfull','fast soil carbon', 'kg C/m2')
  call write_tile_data_r1d_fptr(unit,'ssc',          soil_slow_soil_C_ptr,'zfull','slow soil carbon', 'kg C/m2')
  
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
     __NF_ASRT__(nf_close(unit))
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
     __DEBUG1__(soil%tag)
     __DEBUG1__(soil%pars%k_sat_ref)
     __DEBUG1__(soil%pars%psi_sat_ref)
     __DEBUG1__(soil%pars%chb)
     __DEBUG1__(soil%pars%vwc_sat)
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
                           soil_levap, soil_fevap, soil_melt, &
                           soil_lrunf, soil_hlrunf, soil_Ttop, soil_Ctop, &
                           soil_frunf, soil_hfrunf)
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
       soil_hfrunf   ! ?? heat associated with frozen runoff from soil [W/m^2]

  ! ---- local vars ----------------------------------------------------------
  real, dimension(num_l)   :: del_t, &! ?? temperature tendency [K]
       psi, &! soil moisture potential wrt local elevation [m]
       DThDP, & ! ?? deriv. of vol. liq. cont. wrt psi [1/m]
       K_z, K_x, &! soil horiz. and vert. hydraulic conductivity, respectively [mm/s]
       DKDP, & ! ?? deriv. of hyd. cond. wrt psi [kg/m^3]
       vlc, & ! volumetric liquid water content [-]
       vsc, & ! volumetric solid water content [-]
       dW_l, & ! tendency of soil water mass [mm]
       DPsi, dW_l_accum
  real, dimension(num_l+1) :: flow, & ! downwards flow at layer interface above [mm/timestep]
       infilt, flow_accum
  real, dimension(num_l  ) :: div, & ! total divergence of soil water [mm/s]
       div_it    ! divergence of water due to inter-tile flow (incl. to stream)
  ! set in hlsp_hydrology_1 [mm/s]
  real, dimension(num_l  ) :: hdiv_it, &! divergence of heat due to inter-tile water flow [W/m^2]
       div_bf, & ! baseflow [mm/s]
       div_if, & ! interlow [mm/s]
       div_al, & ! div from active layer [mm/s]
       dq, div_active, &
                              air_depth, macro_frac, extra

    real      :: &
       lprec_eff, & ! infiltrating throughfall (less saturated runoff) [mm/s], and
       hlprec_eff, & !  associated heat [W/m^2]
       tflow, & ! assumed temperature of liquid flow entering top layer [K]
       hcap, & ! heat capacity [J/m^2.K]
       cap_flow, dheat, &
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
       c0, c1, c2, Dpsi_min, Dpsi_max, dW_W, dt_richards, lrunf_ie_accum, &
     sat_area_frac, sat_thick, sum_trans, &
       gw_flux, depth_to_wt2_3, depth_to_wt_apparent, &
       depth_to_gw_flow, depth_to_gw_flow_3, z_bot, &
     active_layer_thickness, depth_to_cf, d_psi, d_psi_s, psi_star, &
     depth_to_cf_1, depth_to_cf_3, &
       depth_to_wt_1, depth_to_wt_2, depth_to_wt_2a, depth_to_wt_2b, depth_to_wt_3, depth_to_wt_4, &
       storage_2, deficit_2, deficit_3, deficit_4, deficit, dum1, dum2, dum3
  logical :: stiff
  logical :: hlsp_stiff ! were all the horizontal conductivities zero at the call to hlsp_hydrology_1?
  real :: zimh, ziph, dTr_g(num_l), dTr_s(num_l)
  integer :: n_iter, l, l_max_active_layer, istep, nstep
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
  real :: Theta ! for debug primtout only
  integer :: ic ! cohort iterator
  
  ! For testing tridiagonal solution
  real, dimension(num_l)   :: t_soil_tridiag ! soil temperature based on generic tridiagonal solution [K]
  real, dimension(num_l)   :: t_diff ! difference from original advection subroutine [K]

  real :: surface_water ! diagnostic surface water storage [m]
  real :: inundated_frac ! diagnostic inundated area fraction [-]
  real :: wet_frac ! diagnostic wetland area fraction
  real :: sat_runf_frac ! [-] effective saturated fraction used for calculating surface runoff, used
  ! when simple inundation scheme is applied
  real :: flow_s(num_l) ! flow downwards at layer interfaces [mm/s] (flow/delta_time excluding bottom interface)
  real :: reflux        ! [mm/s] upwards flow at soil surface in excess of rejected throughfall
#ifdef ZMSDEBUG
  ! Checking diagnostics from Richards that are used in advection solution.
  real, dimension(num_l)   :: wl_before ! water content before call to Richards [mm]
  real :: w1, w2
#endif
  real, parameter :: wthresh = 1.e-9   ! [mm] tolerance for roundoff error for water balance
  real :: wsum1, wsum2   ! total water stocks for balance check [mm, or kg/m^2]
  real :: sliq, sice     ! for call to soil_tile_stock_pe
  real, parameter :: thetathresh = -0.01 ! [-] threshold for negative soil liquid water volume
  ! before warning or abort
  real, parameter :: negrnuthresh = -0.1 ! [mm/s] threshold for negative lrunf_nu
  ! before warning or abort
  character(len=512) :: mesg ! for error message
  integer :: ipt, jpt, kpt, fpt ! for debug
  real :: neg_wl ! [mm] negative wl for zeroing
  real :: wl_ssat ! [mm] excess wl to push up for supersupersat
  real :: psimax ! [m] maximum feasible physical psi, above which push up excess to smooth numerics
  real :: Xmax   ! [-] theta associated with psimax
  real :: wl_max ! [mm] water associated with Xmax
  real :: bwood  ! woody biomass, [kgC/m2], used only as macroporosity flag
  ! --------------------------------------------------------------------------
  div_active(:) = 0.0

  
  if(is_watch_point()) then
     write(*,*) ' ##### soil_step_2 checkpoint 1 #####'
     write(*,*) 'mask    ', .true.
     __DEBUG1__(subs_evap)
     __DEBUG1__(snow_lprec)
     __DEBUG1__(vegn_uptk)
     __DEBUG1__(subs_M_imp)
     __DEBUG1__(soil%pars%vwc_sat)
     do l = 1, num_l
        Theta = (soil%ws(l)+soil%wl(l))/(dens_h2o*dz(l))
        write(*,'(x,i2.2,x)',advance='NO')l
        __DEBUG5__(soil%T(l),Theta,soil%wl(l),soil%ws(l),soil%groundwater(l))
     enddo
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
     end do
  end if

  if(is_watch_point()) then
     write(*,*) ' ##### soil_step_2 checkpoint 2 #####'
     do l = 1, num_l
        write(*,'(x,i2.2,x)',advance='NO') l
        __DEBUG4__(soil%T(l),del_t(l),soil%e(l),soil%f(l))
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
  enddo

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
     ! and don't change the soil%uptake_T
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
     hcap = soil%heat_capacity_dry(l)*dz(l) &
          + clw*soil%wl(l) + csw*soil%ws(l)

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
     else
        depth_to_wt2_3 = zhalf(l+1)
        depth_to_cf_3 = zhalf(l+1)
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
     if (vsc(min(num_l,layer_for_gw_switch)).gt.0.) then   ! permafrost
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
     IF(stiff .AND. BYPASS_RICHARDS_WHEN_STIFF) THEN
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
#ifdef ZMSDEBUG
        wl_before(1:num_l) = soil%wl(1:num_l)
#endif
        CALL RICHARDS(soil, psi, DThDP, K_z, DKDP, div, &
                      lprec_eff, Dpsi_min, Dpsi_max, delta_time, stiff, &
                      dPsi, dW_l, flow, lrunf_ie)

        IF (SPLIT_DELTA_T) THEN
           dW_W = 0.
           do l = 1, num_l
              if (soil%wl(l)/dz(l).gt.Wl_min) &
                       dW_W = max(dW_W, -dW_l(l)/soil%wl(l))
           enddo
           IF (dW_W.gt.dW_W_max) THEN
              write(*,*) 'splitting time step'
              call get_current_point(ipt,jpt,kpt,fpt)
              __DEBUG4__(ipt,jpt,kpt,fpt)
              __DEBUG1__(dW_W)
              do l=1,num_l
                 __DEBUG2__(soil%wl(l),dW_l(l))
              enddo
              nstep = 1+int(dW_W/dW_W_max)
              dt_richards = delta_time / nstep
              lrunf_ie_accum = 0.
              flow_accum = 0.
              do istep = 1, nstep
                 do l = 1, num_l
                    vlc(l) = max(0., soil%wl(l) / (dens_h2o*dz(l)))
                    vsc(l) = max(0., soil%ws(l) / (dens_h2o*dz(l)))
                 enddo
                 call soil_data_hydraulic_properties (soil, vlc, vsc, &
                            psi, DThDP, K_z, K_x, DKDP, Dpsi_min, Dpsi_max )
                 CALL RICHARDS(soil, psi, DThDP, K_z, DKDP, div, &
                            lprec_eff, Dpsi_min, Dpsi_max, dt_richards, stiff, &
                            dPsi, dW_l, flow, lrunf_ie)
                 dW_W = 0.
                 do l = 1, num_l
                    if (soil%wl(l)/dz(l).gt.Wl_min) &
                            dW_W = max(dW_W, -dW_l(l)/soil%wl(l))
                 enddo
                 if (dW_W.gt.dW_W_max) &
                    call error_mesg('RICHARDS EQUATION SOLVER', &
                       'EXCESSIVE DEWATERING DURING SPLIT TIME STEP', FATAL)
                 do l = 1, num_l
                    soil%wl(l) = soil%wl(l) + dW_l(l)
                 enddo
                 lrunf_ie_accum = lrunf_ie_accum + lrunf_ie
                 flow_accum = flow_accum + flow
              enddo
              flow = flow_accum
              lrunf_ie = lrunf_ie_accum / float(nstep)
           ELSE
              do l = 1, num_l
                 soil%wl(l) = soil%wl(l) + dW_l(l)
              enddo
           ENDIF
        ELSE
           ! do not split delta_t
           do l = 1, num_l
              soil%wl(l) = soil%wl(l) + dW_l(l)
           enddo
        ENDIF
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
     if (soil%wl(l) < 0.) then
        if (soil%wl(l) / (dens_h2o*dz(l)*soil%pars%vwc_sat) < thetathresh) then
           call get_current_point(ipt, jpt, kpt, fpt)
           write(mesg,*) 'soil%wl(l) < 0! l,i,j,k,face:', l, ipt, jpt, kpt, fpt, '. theta = ', &
                soil%wl(l) / (dens_h2o*dz(l)*soil%pars%vwc_sat), '. If ".not. allow_neg_wl", '// &
                'model will abort.'
           if (.not. allow_neg_wl) then
              call error_mesg(module_name, mesg, FATAL)
           else
              ! Make sure bulk heat capacity stays above zero
              hcap = soil%heat_capacity_dry(l)*dz(l) &
                     + clw*soil%wl(l) + csw*soil%ws(l)
              if (hcap > 0.) then
                 call error_mesg(module_name, mesg, WARNING)
              else
                 write(mesg,*) 'soil%wl(l) < 0! l,i,j,k,face:', l, ipt, jpt, kpt, fpt, '. theta = ', &
                        soil%wl(l) / (dens_h2o*dz(l)*soil%pars%vwc_sat), '. This makes hcap = ', hcap, &
                        ', which is < 0! Model Aborting!'
                 call error_mesg(module_name, mesg, FATAL)
              end if         
           end if
        end if
     end if
  end do
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
        write(*,'(i2.2,x)')l
        call dpri('T=',soil%T(l))
        call dpri('Th=',(soil%ws(l)+soil%wl(l))/(dens_h2o*dz(l)))
        call dpri('wl=',soil%wl(l))
        call dpri('ws=',soil%ws(l))
        call dpri('gw=',soil%groundwater(l))
     enddo
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

  if (.not. LM2) call send_tile_data(id_psi_bot, soil%psi(num_l), diag)
end subroutine soil_step_2

! ============================================================================

subroutine soil_step_3(soil, diag)
  type(soil_tile_type), intent(in) :: soil
  type(diag_buff_type), intent(inout) :: diag

  if(id_fast_soil_C>0) call send_tile_data(id_fast_soil_C, soil%fast_soil_C(:)/dz(1:num_l), diag)
  if(id_slow_soil_C>0) call send_tile_data(id_slow_soil_C, soil%slow_soil_C(:)/dz(1:num_l), diag)
  if (id_fsc > 0)      call send_tile_data(id_fsc, sum(soil%fast_soil_C(:)), diag)
  if (id_ssc > 0)      call send_tile_data(id_ssc, sum(soil%slow_soil_C(:)), diag)
end subroutine soil_step_3


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
  integer :: l

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

  if(is_watch_cell()) then
      write(*,*) ' ##### soil_step_2 checkpoint 3.001 #####'
      write(*,*) 'For watch_cell'
      __DEBUG4__(l,summax,liq_frac,soil%vwc_max(l))
      __DEBUG4__(excess_liq,excess_ice,dens_h2o,dz(l))
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


  if(is_watch_cell()) then
     write(*,*) ' ##### soil_step_2 checkpoint 3.01 #####'
     __DEBUG2__(lrunf_nu, hlrunf_nu)
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
subroutine RICHARDS(soil, psi, DThDP, hyd_cond, DKDP, div, &
                lprec_eff, Dpsi_min, Dpsi_max, dt_richards, stiff, &
                 dPsi, dW_l, flow, lrunf_ie)
! ZMS How to set prevent lrunf_ie when surface soil water < pond?

  type(soil_tile_type), intent(inout)   :: soil
  real, intent(in),  dimension(num_l)   :: psi, DThDP, hyd_cond, DKDP, div
  real, intent(in)                      :: lprec_eff, Dpsi_min, Dpsi_max
  real, intent(in)                      :: dt_richards
  logical, intent(in)                   :: stiff
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
!      where (stiff)
  bbb = xxx - (-K(l-1)/del_z(l-1) - DKDPp(l-1)*grad(l-1) )
  ddd = - K(l-1) *grad(l-1) - div(l)
!        elsewhere
!          Qout = hyd_cond(l) ! gravity drainage
!          DQoutDP = DKDP(l)  ! gravity drainage
!          Qout = 0.                ! no drainage
!          DQoutDP = 0.             ! no drainage
!          where (psi(l).gt.0.) ! linear baseflow from gw
!              Qout = 0.15*psi(l)/soil%pars%tau_groundwater
!              DQoutDP = 0.15/soil%pars%tau_groundwater
!            elsewhere
!              Qout = 0.
!              DQoutDP = 0.
!            endwhere
!          bbb = xxx - (-K(l-1)/del_z(l-1) - DKDPp(l-1)*grad(l-1)&
!                      -DQoutDP )
!          ddd = -Qout - K(l-1) *grad(l-1)
!        endwhere
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
  IF (bbb+ccc*eee(l) .NE. 0.) THEN
     dPsi(l) = (ddd-ccc*fff(l))/(bbb+ccc*eee(l))
     if(is_watch_point()) then
         write(*,*) 'bbb+ccc*eee(l) .NE. 0.'
         __DEBUG1__(bbb)
         __DEBUG1__(ccc)
         __DEBUG1__(ddd)
         __DEBUG1__(eee(l))
         __DEBUG1__(fff(l))
         __DEBUG1__(dPsi(l))
         __DEBUG1__(stiff)
         __DEBUG1__(dPsi(l))
         __DEBUG1__(Dpsi_min)
         __DEBUG1__(Dpsi_max)
     endif
     if (.not.stiff .and. dPsi(l).gt.Dpsi_min .and. dPsi(l).lt.Dpsi_max) then
        lrunf_ie = 0.
     else
      	if (stiff) then
           dPsi(l) = - psi(l)
        else
	   if (dPsi(l).lt.Dpsi_min) then
               flag = .true.
               call get_current_point(ipt,jpt,kpt,fpt)
               write(*,*) 'note 1: at point ',ipt,jpt,kpt,fpt,'computed dPsi(1) too negative'
               write(*,*) 'note 1: at point ',ipt,jpt,kpt,fpt,'dPsi(1)=',dPsi(l)
               write(*,*) 'note 1: at point ',ipt,jpt,kpt,fpt,'dPsi_min =',dPsi_min
           endif
           dPsi(l) = min (dPsi(l), Dpsi_max)
           dPsi(l) = max (dPsi(l), Dpsi_min)
        endif
        flow(l) = (dPsi(l)*(bbb+ccc*eee(l))+ccc*fff(l) &
                      - K(l)*grad(l))*dt_richards
        lrunf_ie = lprec_eff - flow(l)/dt_richards
        if (.not.allow_negative_rie.and.lrunf_ie.lt.-lrunf_ie_tol) then
           flag = .true.
           call get_current_point(ipt,jpt,kpt,fpt)
           dpsi_alt = (ddd-ccc*fff(l))/(bbb+ccc*eee(l))
           write(*,*) 'note 1: at point ',ipt,jpt,kpt,fpt,'rie= ',lrunf_ie,' reset to 0'
           write(*,*) 'note 1: at point ',ipt,jpt,kpt,fpt,'dPsi=',dPsi(l), ' reset to ',dpsi_alt
           write(*,*) 'note 1: at point ',ipt,jpt,kpt,fpt,'dPsi_min/max=',dPsi_min,dPsi_max
           dPsi(l) = dpsi_alt
           lrunf_ie = 0.
           flow(l) = lprec_eff*dt_richards
        endif
     endif
  ELSE
     if(is_watch_point()) then
        write(*,*) 'bbb+ccc*eee(l) .EQ. 0.'
        __DEBUG1__(bbb)
        __DEBUG1__(ccc)
        __DEBUG1__(ddd)
        __DEBUG1__(eee(l))
        __DEBUG1__(fff(l))
        __DEBUG1__(dPsi(l))
        __DEBUG1__(stiff)
        __DEBUG1__(dPsi(l))
        __DEBUG1__(Dpsi_min)
        __DEBUG1__(Dpsi_max)
     endif
     if (stiff) then
        dPsi(l) = - psi(l)
     else
         dPsi(l) = Dpsi_max
     endif
     flow(l) = (dPsi(l)*(bbb+ccc*eee(l))+ccc*fff(l) &
                - K(l)*grad(l))*dt_richards
     lrunf_ie = lprec_eff - flow(l)/dt_richards
     if (.not.allow_negative_rie.and.lrunf_ie.lt.-lrunf_ie_tol) then
        flag = .true.
        call get_current_point(ipt,jpt,kpt,fpt)
        ! next change will not change answers in previous runs, since old version would crash
        ! the only time this point was reached was when DThDP was zero everywhere.
        !     dpsi_alt = (ddd-ccc*fff(l))/(bbb+ccc*eee(l))
        dpsi_alt = 0.
        write(*,*) 'note 2: at point ',ipt,jpt,kpt,fpt,'rie=',lrunf_ie,' reset to 0'
        write(*,*) 'note 2: at point ',ipt,jpt,kpt,fpt,'dPsi=',dPsi(l),' reset to',dpsi_alt
        write(*,*) 'note 2: at point ',ipt,jpt,kpt,fpt,'dPsi_min/max=',dPsi_min,dPsi_max
        dPsi(l) = dpsi_alt
        lrunf_ie = 0.
        flow(l) = lprec_eff*dt_richards
     endif
  ENDIF
      
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

  do l = 2, num_l
     dPsi(l) = eee(l-1)*dPsi(l-1) + fff(l-1)
  enddo
  
  l_internal = 1
  dW_l_internal = -1.e20
  do l = 1, num_l-1
     flow(l+1) = dt_richards*( &
         -K(l)*(grad(l)&
         +(DPsi(l+1)-DPsi(l))/ del_z(l)) &
         -grad(l)*(DKDPp(l)*Dpsi(l+1)+ &
                         DKDPm(l)*Dpsi(l) )  )
     dW_l(l) = flow(l) - flow(l+1) - div(l)*dt_richards
     if (flag .and. l.gt.1. .and. dW_l(l).gt.dW_l_internal) then
        l_internal = l
        dW_l_internal = dW_l(l)
     endif
  enddo
  flow(num_l+1) = 0.
  dW_l(num_l) = flow(num_l) - flow(num_l+1) &
                        - div(num_l)*dt_richards
  if (flag .and. dW_l(num_l).gt.dW_l_internal) then
     l_internal = num_l
     dW_l_internal = dW_l(num_l)
  endif

  if(is_watch_point().or.(flag.and.write_when_flagged)) then
     write(*,*) ' ##### soil_step_2 checkpoint 3.21 #####'
     do l = 1, num_l
        write(*,'(i2.3)', advance='NO') l
        call dpri(' dW_l=',dW_l(l))
        call dpri(' flow=',flow(l))
        call dpri(' div=', div(l))
        write(*,*)
     enddo
  endif

  if (flag) then
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
! tile existence detector: returns a logical value indicating wether component
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

DEFINE_SOIL_ACCESSOR_1D(real,T)
DEFINE_SOIL_ACCESSOR_1D(real,wl)
DEFINE_SOIL_ACCESSOR_1D(real,ws)
DEFINE_SOIL_ACCESSOR_1D(real,groundwater)
DEFINE_SOIL_ACCESSOR_1D(real,groundwater_T)
DEFINE_SOIL_ACCESSOR_1D(real,w_fc)
DEFINE_SOIL_ACCESSOR_1D(real,alpha)
DEFINE_SOIL_ACCESSOR_0D(real,uptake_T)
DEFINE_SOIL_ACCESSOR_0D(integer,tag)
DEFINE_SOIL_ACCESSOR_1D(real,fast_soil_C)
DEFINE_SOIL_ACCESSOR_1D(real,slow_soil_C)
DEFINE_SOIL_ACCESSOR_1D(real,asoil_in)
DEFINE_SOIL_ACCESSOR_1D(real,fsc_in)
DEFINE_SOIL_ACCESSOR_1D(real,ssc_in)

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

end module soil_mod
