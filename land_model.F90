! ============================================================================
! top-level core of the Land Dynamics (LaD) model code
! ============================================================================
module land_model_mod

#include "shared/debug.inc"

use time_manager_mod, only : time_type, get_time, increment_time, time_type_to_real, &
     get_date, operator(+), operator(-)
use mpp_domains_mod, only : domain2d, domainUG, mpp_get_ntile_count, &
     mpp_pass_SG_to_UG, mpp_pass_ug_to_sg, &
     mpp_get_UG_domain_tile_pe_inf, mpp_get_UG_domain_ntiles, &
     mpp_get_UG_compute_domain, mpp_get_UG_domain_grid_index

#ifdef INTERNAL_FILE_NML
use mpp_mod, only: input_nml_file
#else
use fms_mod, only: open_namelist_file
#endif

use mpp_mod, only : mpp_max, mpp_sum, mpp_chksum
use fms_io_mod, only : read_compressed, fms_io_unstructured_read
use fms_mod, only : error_mesg, FATAL, WARNING, NOTE, mpp_pe, &
     mpp_root_pe, file_exist, check_nml_error, close_file, &
     stdlog, stderr, mpp_clock_id, mpp_clock_begin, mpp_clock_end, string, &
     stdout, CLOCK_FLAG_DEFAULT, CLOCK_COMPONENT, CLOCK_ROUTINE
use data_override_mod, only : data_override_ug
use diag_manager_mod, only : diag_axis_init, register_static_field, &
     register_diag_field, send_data, diag_field_add_attribute
use diag_axis_mod, only: diag_axis_add_attribute
use constants_mod, only : radius, hlf, hlv, hls, tfreeze, pi, rdgas, rvgas, cp_air, &
     stefan
use astronomy_mod, only : astronomy_init, diurnal_solar
use sphum_mod, only : qscomp
use tracer_manager_mod, only : NO_TRACER

use land_constants_mod, only : NBANDS, BAND_VIS, BAND_NIR, mol_air, mol_C, mol_co2, d608
use land_tracers_mod, only : land_tracers_init, land_tracers_end, ntcana, isphum, ico2
use land_tracer_driver_mod, only: land_tracer_driver_init, land_tracer_driver_end, &
     update_cana_tracers
use glacier_mod, only : read_glac_namelist, glac_init, glac_end, glac_get_sfc_temp, &
     glac_radiation, glac_step_1, glac_step_2, save_glac_restart
use lake_mod, only : read_lake_namelist, lake_init, lake_end, lake_get_sfc_temp, &
     lake_radiation, lake_step_1, lake_step_2, save_lake_restart
use soil_mod, only : read_soil_namelist, soil_init, soil_end, soil_get_sfc_temp, &
     soil_radiation, soil_step_1, soil_step_2, soil_step_3, save_soil_restart, &
     ! moved here to eliminate circular dependencies with hillslope mods:
     soil_cover_cold_start, retrieve_soil_tags
use soil_carbon_mod, only : read_soil_carbon_namelist, N_C_TYPES, soil_carbon_option, &
    SOILC_CORPSE_N
use snow_mod, only : read_snow_namelist, snow_init, snow_end, snow_get_sfc_temp, &
     snow_get_depth_area, snow_step_1, snow_step_2, &
     save_snow_restart, sweep_tiny_snow
use vegn_data_mod, only : LU_PAST, LU_CROP, LU_NTRL, LU_SCND, LU_RANGE, LU_URBN
use vegetation_mod, only : read_vegn_namelist, vegn_init, vegn_end, &
     vegn_radiation, vegn_diffusion, vegn_step_1, vegn_step_2, vegn_step_3, &
     update_derived_vegn_data, update_vegn_slow, save_vegn_restart, &
     cohort_test_func, cohort_area_frac, any_vegn, is_tree, is_grass, is_c3, is_c4, &
     is_c3grass, is_c4grass
use vegn_disturbance_mod, only : vegn_nat_mortality_ppa
use vegn_fire_mod, only : update_fire_fast, fire_transitions, save_fire_restart
use cana_tile_mod, only : canopy_air_mass, canopy_air_mass_for_tracers, cana_tile_heat, cana_tile_carbon
use canopy_air_mod, only : read_cana_namelist, cana_init, cana_end,&
     cana_roughness, &
     save_cana_restart
use river_mod, only : river_init, river_end, update_river, river_stock_pe, &
     save_river_restart, river_tracers_init, num_river_tracers, river_tracer_index, &
     river_tracer_names, get_river_water
use topo_rough_mod, only : topo_rough_init, topo_rough_end, update_topo_rough
use soil_tile_mod, only : soil_tile_stock_pe, soil_tile_heat, soil_roughness
use vegn_cohort_mod, only : plant_C
use vegn_tile_mod, only : vegn_cover_cold_start, &
                          vegn_tile_stock_pe, vegn_tile_heat, vegn_tile_carbon
use lake_tile_mod, only : lake_cover_cold_start, lake_tile_stock_pe, &
                          lake_tile_heat, lake_roughness
use glac_tile_mod, only : glac_pars_type, glac_cover_cold_start, &
                          glac_tile_stock_pe, glac_tile_heat, glac_roughness
use snow_tile_mod, only : snow_tile_stock_pe, snow_tile_heat, snow_roughness, snow_radiation
use land_numerics_mod, only : ludcmp, lubksb, lubksb_and_improve, nearest, &
     horiz_remap_type, horiz_remap_new, horiz_remap, horiz_remap_del, &
     horiz_remap_print
use land_io_mod, only : read_land_io_namelist, input_buf_size, new_land_io
use land_tile_mod, only : land_tile_map, land_tile_type, land_tile_list_type, &
     land_tile_enum_type, new_land_tile, insert, remove, empty, nitems, &
     first_elmt, tail_elmt, next_elmt, operator(==), current_tile, &
     get_tile_tags, get_tile_water, land_tile_heat, land_tile_nitrogen, &
     land_tile_carbon, max_n_tiles, init_tile_map, free_tile_map, &
     loop_over_tiles, land_tile_list_init, land_tile_list_end, &
     merge_land_tile_into_list, tile_test_func
use land_data_mod, only : land_data_type, atmos_land_boundary_type, &
     land_state_type, land_data_init, land_data_end, lnd, log_version
use nf_utils_mod,  only : nfu_inq_var, nfu_inq_dim, nfu_get_var
use land_utils_mod, only : put_to_tiles_r0d_fptr
use land_tile_io_mod, only: land_restart_type, &
     init_land_restart, open_land_restart, save_land_restart, free_land_restart, &
     add_tile_data, add_int_tile_data, get_tile_data, &
     field_exists, print_netcdf_error
use land_tile_diag_mod, only : OP_SUM, cmor_name, tile_diag_init, tile_diag_end, &
     register_tiled_diag_field, send_tile_data, dump_tile_diag_fields, &
     add_tiled_diag_field_alias, register_cohort_diag_field, send_cohort_data, &
     set_default_diag_filter, register_tiled_area_fields, send_global_land_diag, &
     get_area_id
use land_debug_mod, only : land_debug_init, land_debug_end, set_current_point, &
     is_watch_point, is_watch_cell, is_watch_time, get_watch_point, get_current_point, &
     check_conservation, do_check_conservation, water_cons_tol, carbon_cons_tol, nitrogen_cons_tol, &
     check_var_range, check_temp_range, current_face, log_date, land_error_message
use static_vegn_mod, only : write_static_vegn
use land_transitions_mod, only : &
     land_transitions_init, land_transitions_end, land_transitions, &
     save_land_transitions_restart
use stock_constants_mod, only: ISTOCK_WATER, ISTOCK_HEAT, ISTOCK_SALT
use nitrogen_sources_mod, only : nitrogen_sources_init, nitrogen_sources_end, &
     update_nitrogen_sources, nitrogen_sources
use hillslope_mod, only: retrieve_hlsp_indices, save_hlsp_restart, hlsp_end, &
                         read_hlsp_namelist, hlsp_init, hlsp_config_check
use hillslope_hydrology_mod, only: hlsp_hydrology_1, hlsp_hydro_init

implicit none
private

! ==== public interfaces =====================================================
public land_model_init          ! initialize the land model
public land_model_end           ! finish land model calculations
public land_model_restart       ! saves the land model restart(s)
public update_land_model_fast   ! time-step integration
public update_land_model_slow   ! time-step integration
public atmos_land_boundary_type ! data from coupler to land
public land_data_type           ! data from land to coupler
public land_data_type_chksum    ! routine to print checksums for land_data_type
public atm_lnd_bnd_type_chksum  ! routine to print checksums for atmos_land_boundary_type

public :: Lnd_stock_pe          ! return stocks of conservative quantities

! re-export land diagnostic subroutines for tiled diag in flux exchange
public set_default_diag_filter, register_tiled_diag_field, send_tile_data, dump_tile_diag_fields
public send_global_land_diag
! ==== end of public interfaces ==============================================

! ==== module constants ======================================================
character(len=*), parameter :: module_name = 'land'
#include "shared/version_variable.inc"

! ==== module variables ======================================================

! ---- namelist --------------------------------------------------------------
logical :: use_old_conservation_equations  = .false.
logical :: lm2                             = .false.
logical :: give_stock_details              = .false.
logical :: use_tfreeze_in_grnd_latent      = .false.
logical :: use_atmos_T_for_precip_T        = .false.
logical :: use_atmos_T_for_evap_T          = .false.
real    :: cpw = 1952.  ! specific heat of water vapor at constant pressure
real    :: clw = 4218.  ! specific heat of water (liquid)
real    :: csw = 2106.  ! specific heat of water (ice)
real    :: min_sum_lake_frac = 1.e-8
real    :: min_frac = 0.0 ! minimum fraction of soil, lake, and glacier that is not discarded on cold start
real    :: gfrac_tol         = 1.e-6
real    :: discharge_tol = -1.e20
logical :: improve_solution = .FALSE. ! if true, the solution is improved to reduce
                                      ! numerical errors (important for PPA with meny cohorts)
integer :: max_improv_steps = 5       ! max number of solution improvement steps
real    :: solution_tol     = 1e-10   ! tolerance for solution improvement
real    :: con_fac_large = 1.e6
real    :: con_fac_small = 1.e-6
real    :: tau_snow_T_adj = -1.0 ! time scale of snow temperature adjustment
              ! for the snow-free surface (s); negative means no adjustment
logical :: prohibit_negative_canopy_water = .TRUE. ! if true, the solution of energy/water
              ! balance is iterated (at most max_canopy_water_steps) to ensure
              ! water and snow on leaves do not go negative
integer :: max_canopy_water_steps = 2
real    :: lw_delta_T_thresh = 1.0e36 ! temperature change during time step that
              ! triggers recalculation of longwave radiation derivatives for improvement
              ! of its linearization. Note that linearization around temperature at the start
              ! if time step leads to LW emission underestimate of about 20 W/m2 at
              ! delta_T = 25K. Huge default value turns LW linearization improvement off
character(16) :: nearest_point_search = 'global' ! specifies where to look for
              ! nearest points for missing data, "global" or "face"
logical :: print_remapping = .FALSE. ! if true, full land cover remapping
              ! information is printed on the cold start
integer :: layout(2) = (/0,0/)
integer :: io_layout(2) = (/0,0/)
integer :: npes_io_group = 0
  ! mask_table contains information for masking domain ( n_mask, layout and mask_list).
  !   A text file to specify n_mask, layout and mask_list to reduce number of processor
  !   usage by masking out some domain regions which contain all ocean points.
  !   The default file name of mask_table is "INPUT/land_mask_table". Please note that
  !   the file name must begin with "INPUT/". The first
  !   line of mask_table will be number of region to be masked out. The second line
  !   of the mask_table will be the layout of the model. User need to set land_model_nml
  !   variable layout to be the same as the second line of the mask table.
  !   The following n_mask line will be the position of the processor and face number
  !   to be masked out (First entry is i-direction postion, second entry is
  !   j-direction position, third entry is face number ).
  !   The mask_table could be created by tools check_mask.
  !   For example the mask_table will be as following if n_mask=4, layout=4,6 and
  !   the processor (1,2) and (3,6) will be masked out.
  !     2
  !     4,6
  !     1,2,1
  !     3,6,2
  !     4,5,5
  !     2,4,6
character(len=128) :: mask_table = "INPUT/land_mask_table"
logical :: reset_to_ntrl = .FALSE. ! if true, on model initialization any vegetation tiles
  ! that are not natural are removed and fraction of natural vegetation tiles is scaled up
  ! to occupy the entire soil. If a grid cell does not have any natural tiles, then the tile
  ! with highest biomass becomes natural.
  ! This is intended to be used in one-shot runs to convert existing IC with land use into
  ! IC suitable for potential vegetation runs. Typically one would set this parameter to
  ! TRUE, run for zero time steps, and collect resulting restarts to be an IC for potential
  ! vegetation run.

namelist /land_model_nml/ use_old_conservation_equations, &
                          lm2, give_stock_details, &
                          use_tfreeze_in_grnd_latent, &
                          use_atmos_T_for_precip_T, &
                          use_atmos_T_for_evap_T, &
                          cpw, clw, csw, min_sum_lake_frac, min_frac, &
                          gfrac_tol, discharge_tol, &
                          improve_solution, solution_tol, max_improv_steps, &
                          lw_delta_T_thresh, &
                          con_fac_large, con_fac_small, &
                          tau_snow_T_adj, prohibit_negative_canopy_water, max_canopy_water_steps, &
                          nearest_point_search, print_remapping, &
                          layout, io_layout, npes_io_group, mask_table, reset_to_ntrl
! ---- end of namelist -------------------------------------------------------

logical  :: module_is_initialized = .FALSE.
logical  :: stock_warning_issued  = .FALSE.
logical  :: update_cana_co2 ! if false, cana_co2 is not updated during the model run.
real     :: delta_time ! duration of main land time step (s)
real     :: steps_per_day ! number of fast time steps per day

! ---- indices of river tracers
integer  :: n_river_tracers
integer  :: i_river_ice, i_river_heat, i_river_DOC

! ---- diag field IDs --------------------------------------------------------
integer :: &
 ! COLUMN        VEGN        SNOW      GLAC/LAKE/SOIL  CANOPY-AIR  RIVER
  id_VWS,                                               id_VWSc,           &
  id_LWS,      id_LWSv,     id_LWSs,     id_LWSg,                          &
  id_FWS,      id_FWSv,     id_FWSs,     id_FWSg,                          &
  id_HS,       id_HSv,      id_HSs,      id_HSg,        id_HSc,            &
!
  id_precip,                                                               &
  id_hprec,                                                                &
  id_lprec,    id_lprecv,   id_lprecs,   id_lprecg,                        &
  id_hlprec,   id_hlprecv,  id_hlprecs,  id_hlprecg,                       &
  id_fprec,    id_fprecv,   id_fprecs,                                     &
  id_hfprec,   id_hfprecv,  id_hfprecs,                                    &
  id_evap,                                                                 &
  id_hevap,                                                                &
  id_levap,    id_levapv,   id_levaps,   id_levapg,                        &
  id_hlevap,   id_hlevapv,  id_hlevaps,  id_hlevapg,                       &
  id_fevap,    id_fevapv,   id_fevaps,   id_fevapg,                        &
  id_hfevap,   id_hfevapv,  id_hfevaps,  id_hfevapg,                       &
  id_runf,                                                                 &
  id_hrunf,                                                                &
  id_lrunf,                 id_lrunfs,   id_lrunfg,                        &
  id_hlrunf,                id_hlrunfs,  id_hlrunfg,                       &
  id_frunf,                 id_frunfs,                                     &
  id_hfrunf,                id_hfrunfs,                                    &
  id_melt,     id_meltv,    id_melts,    id_meltg,                         &
  id_fsw,      id_fswv,     id_fsws,     id_fswg,                          &
  id_flw,      id_flwv,     id_flws,     id_flwg,                          &
  id_sens,     id_sensv,    id_senss,    id_sensg,                         &
!
  id_e_res_1,  id_e_res_2,  id_cd_m,     id_cd_t,                          &
  id_cellarea, id_landfrac,                                                &
  id_geolon_t, id_geolat_t,                                                &
  id_frac,     id_area,     id_ntiles,                                     &
  id_z0m,      id_z0s,      id_con_g_h,  id_con_g_v,                       &
  id_transp,                id_wroff,    id_sroff,                         &
  id_htransp,  id_huptake,  id_hroff,    id_gsnow,    id_gequil,           &
  id_grnd_flux,                                                            &
  id_levapg_max,                                                           &
  id_water,    id_snow,     id_snow_frac,                                  &
  id_Trad,     id_Tca,      id_qca,      id_qco2,     id_qco2_dvmr,        &
  id_swdn_dir, id_swdn_dif, id_swup_dir, id_swup_dif, id_lwdn,             &
  id_fco2,     id_co2_mol_flux,                                            &
  id_vegn_cover,    id_vegn_cover_1,  id_vegn_cover_U, id_cosz,            &
  id_albedo_dir,    id_albedo_dif,                                         &
  id_vegn_refl_dir, id_vegn_refl_dif, id_vegn_refl_lw,                     &
  id_vegn_tran_dir, id_vegn_tran_dif, id_vegn_tran_lw,                     &
  id_vegn_sctr_dir,                                                        &
  id_subs_refl_dir, id_subs_refl_dif, id_subs_emis, id_grnd_T, id_total_C, id_total_N, &
  id_water_cons, id_carbon_cons, id_nitrogen_cons, id_grnd_rh, id_cana_rh, id_cTot1
! diagnostic ids for canopy air tracers (moist mass ratio)
integer, allocatable :: id_runf_tr(:), id_dis_tr(:)

! IDs of CMOR/CMIP variables
integer :: id_sftlf, id_sftgif ! static fractions
integer :: id_pcp, id_prra, id_prveg, id_evspsblsoi, id_evspsblveg, &
  id_snw, id_snd, id_snc, id_snm, id_sweLut, id_lwsnl, id_hfdsn, id_tws, &
  id_hflsLut, id_rlusLut, id_rsusLut, id_tslsiLut, id_cLand, id_nbp, &
  id_ec, id_eow, id_esn, id_et, id_nLand, &
! various fractions
  id_vegFrac, id_pastureFrac, id_residualFrac, &
  id_cropFrac, id_cropFracC3, id_cropFracC4, &
  id_grassFrac, id_grassFracC3, id_grassFracC4, &
  id_treeFrac, id_c3pftFrac, id_c4pftFrac, id_nwdFracLut, &
  id_fracLut_psl, id_fracLut_crp, id_fracLut_pst, id_fracLut_urb


! init_value is used to fill most of the allocated boundary condition arrays.
! It is supposed to be double-precision signaling NaN, to trigger a trap when
! the program is compiled with trapping non-initialized values.
! See http://ftp.uniovi.es/~antonio/uned/ieee754/IEEE-754references.html
! real, parameter :: init_value = Z'FFF0000000000001'
real, parameter :: init_value = 0.0

! ---- global clock IDs
integer :: landClock, landFastClock, landSlowClock


! ==== NetCDF declarations ===================================================
include 'netcdf.inc'
#define __NF_ASRT__(x) call print_netcdf_error((x),__FILE__,__LINE__)

contains


! ============================================================================
subroutine land_model_init &
     (cplr2land, land2cplr, time_init, time, dt_fast, dt_slow)
! initialize land model using grid description file as an input. This routine
! reads land grid boundaries and area of land from a grid description file

! NOTES: theoretically, the grid description file can specify any regular
! rectangular grid for land, not just lon/lat grid. Therefore the variables
! "xbl" and "ybl" in NetCDF grid spec file are not necessarily lon and lat
! boundaries of the grid.
!   However, at this time the module land_properties assumes that grid _is_
! lon/lat and therefore the entire module also have to assume that the land
! grid is lon/lat.
!   lon/lat grid is also assumed for the diagnostics, but this is probably not
! so critical.
  type(atmos_land_boundary_type), intent(inout) :: cplr2land ! boundary data
  type(land_data_type)          , intent(inout) :: land2cplr ! boundary data
  type(time_type), intent(in) :: time_init ! initial time of simulation (?)
  type(time_type), intent(in) :: time      ! current time
  type(time_type), intent(in) :: dt_fast   ! fast time step
  type(time_type), intent(in) :: dt_slow   ! slow time step

  ! ---- local vars ----------------------------------------------------------
  integer :: unit, ierr, io
  integer :: id_band, id_zfull ! IDs of land diagnostic axes
  integer :: id_ug !<Unstructured axis id.
  logical :: used                        ! return value of send_data diagnostics routine
  integer :: k,l
  integer :: n_cohorts  ! number of cohorts in the current tile (1 if no vegetation)
  type(land_tile_type), pointer :: tile
  type(land_tile_enum_type) :: ce

  type(land_restart_type) :: restart
  character(*), parameter :: restart_file_name='INPUT/land.res.nc'
  logical :: restart_exists

  ! IDs of local clocks
  integer :: landInitClock

  module_is_initialized = .TRUE.

  ! [1] print out version number
  call log_version(version, module_name, &
  __FILE__)

  ! initialize land model clocks
  landClock      = mpp_clock_id('Land'               ,CLOCK_FLAG_DEFAULT,CLOCK_COMPONENT)
  landFastClock  = mpp_clock_id('Update-Land-Fast'   ,CLOCK_FLAG_DEFAULT,CLOCK_ROUTINE)
  landSlowClock  = mpp_clock_id('Update-Land-Slow'   ,CLOCK_FLAG_DEFAULT,CLOCK_ROUTINE)
  landInitClock  = mpp_clock_id('Land init'          ,CLOCK_FLAG_DEFAULT,CLOCK_ROUTINE)

  call mpp_clock_begin(landInitClock)

  ! [2] read land model namelist
#ifdef INTERNAL_FILE_NML
     read (input_nml_file, nml=land_model_nml, iostat=io)
     ierr = check_nml_error(io, 'land_model_nml')
#else
  if (file_exist('input.nml')) then
     unit = open_namelist_file ( )
     ierr = 1;
     do while (ierr /= 0)
        read (unit, nml=land_model_nml, iostat=io, end=10)
        ierr = check_nml_error (io, 'land_model_nml')
     enddo
10   continue
     call close_file (unit)
  endif
#endif
  if (mpp_pe() == mpp_root_pe()) then
     unit = stdlog()
     write (unit, nml=land_model_nml)
     call close_file (unit)
  endif
  ! initialize astronomy, in case it is not initialized, e.g. when using atmos_null
  call astronomy_init()

  ! initialize land state data, including grid geometry and processor decomposition
  call land_data_init(layout, io_layout, time, dt_fast, dt_slow, mask_table,npes_io_group)

  ! initialize land debug output
  call land_debug_init()

  ! initialize tile-specific diagnostics internals
  call tile_diag_init()

  ! initialize some tracer indices
  call land_tracers_init()
  call river_tracers_init()

  ! read sub-model namelists: then need to be read before component initialization
  ! because they can affect the way cover and tiling is initialized on cold start.
  ! Also, some of them register diagnostic sub-sampling selectors, so they better
  ! be after land_tile_diag_init
  call read_land_io_namelist()
  call read_soil_namelist()
  call read_hlsp_namelist() ! Must be called after read_soil_namelist
  call read_vegn_namelist()
  call read_soil_carbon_namelist()
  call read_lake_namelist()
  call read_glac_namelist()
  call read_snow_namelist()
  call read_cana_namelist()

  delta_time  = time_type_to_real(lnd%dt_fast) ! store in a module variable for convenience
  steps_per_day = 86400.0/delta_time

  call init_tile_map()

  ! initialize subgrid tile distribution
  call open_land_restart(restart,restart_file_name,restart_exists)
  if(restart_exists) then
     call error_mesg('land_model_init',&
          'reading NetCDF restart "'//trim(restart_file_name)//'"',&
          NOTE)
     ! read map of tiles -- retrieve information from
     call land_cover_warm_start(restart)
     ! initialize land model data
     if (field_exists(restart, 'lwup'   )) call get_tile_data(restart,'lwup',   land_lwup_ptr)
     if (field_exists(restart, 'e_res_1')) call get_tile_data(restart,'e_res_1',land_e_res_1_ptr)
     if (field_exists(restart, 'e_res_2')) call get_tile_data(restart,'e_res_2',land_e_res_2_ptr)
  else
     ! initialize map of tiles -- construct it by combining tiles
     ! from component models
     call error_mesg('land_model_init','cold-starting land cover map',NOTE)
     call land_cover_cold_start()
  endif
  call free_land_restart(restart)

  ! initialize land model diagnostics -- must be before *_data_init so that
  ! *_data_init can write static fields if necessary
  call land_sg_diag_init(id_cellarea)
  call land_diag_init( lnd%coord_glonb, lnd%coord_glatb, lnd%coord_glon, lnd%coord_glat, &
      time, id_ug, id_band )

  ! set the land diagnostic axes ids for the flux exchange
  land2cplr%axes = (/id_ug/)
  ! send some static diagnostic fields to output
  if ( id_landfrac > 0 ) used = send_data ( id_landfrac, lnd%ug_landfrac,     lnd%time )
  if ( id_geolon_t > 0 ) used = send_data ( id_geolon_t, lnd%ug_lon*180.0/PI, lnd%time )
  if ( id_geolat_t > 0 ) used = send_data ( id_geolat_t, lnd%ug_lat*180.0/PI, lnd%time )

  ! CMOR variables
  if ( id_sftgif > 0 ) call send_cellfrac_data(id_sftgif,is_glacier)
  if ( id_sftlf > 0 )  used = send_data(id_sftlf,lnd%ug_landfrac*100, lnd%time)

  ! [7] initialize individual sub-models
  call hlsp_init ( id_ug ) ! Must be called before soil_init
  call soil_init ( id_ug, id_band, id_zfull)
  call hlsp_hydro_init (id_ug, id_zfull) ! Must be called after soil_init
  call vegn_init ( id_ug, id_band, id_cellarea )
  call lake_init ( id_ug )
  call glac_init ( id_ug )
  call snow_init ()
  call cana_init ()
  call nitrogen_sources_init ( lnd%time, id_ug )
  call topo_rough_init( lnd%time, lnd%sg_lonb, lnd%sg_latb, lnd%sg_domain, lnd%ug_domain, id_ug)

  call river_init( lnd%sg_lon, lnd%sg_lat, lnd%time, lnd%dt_fast, &
          lnd%sg_domain, lnd%ug_domain, lnd%sg_landfrac, discharge_tol, clw, csw )
  ! initialize river tracer indices
  n_river_tracers = num_river_tracers()
  i_river_ice  = river_tracer_index('ice')
  i_river_heat = river_tracer_index('het')
  i_river_DOC  = river_tracer_index('doc')
  if (i_river_ice  == NO_TRACER) call error_mesg ('land_model_init','required river tracer for ice not found', FATAL)
  if (i_river_heat == NO_TRACER) call error_mesg ('land_model_init','required river tracer for heat not found', FATAL)

  call land_transitions_init (id_ug, id_cellarea)
  ! [8] initialize boundary data
  ! [8.1] allocate storage for the boundary data
  call hlsp_config_check () ! Needs to be done after land_transitions_init and vegn_init
  call land_tracer_driver_init(id_ug)
  call realloc_land2cplr ( land2cplr )
  call realloc_cplr2land ( cplr2land )
  ! [8.2] set the land mask to FALSE everywhere -- update_land_bc_fast
  ! will set it to true where necessary
  land2cplr%mask = .FALSE.
  land2cplr%tile_size = 0.0
  ! [8.3] get the current state of the land boundary for the coupler
  ce = first_elmt(land_tile_map, ls=lnd%ls )
  do while(loop_over_tiles(ce,tile, l,k))
     call set_current_point(l,k)

     ! n_cohorts is calculated and passed down to update_land_bc_fast
     ! for convenience, so that there is no need to make a lot of by-cohort arrays
     ! allocatable -- instead they are created on stack with the size passed
     ! as an argument.
     n_cohorts = 1; if (associated(tile%vegn)) n_cohorts = tile%vegn%n_cohorts
     call update_land_bc_fast (tile, n_cohorts, l,k, land2cplr, is_init=.true.)
  enddo

  ! [8.4] update topographic roughness scaling
  call update_land_bc_slow( land2cplr )

  ! mask error checking
  do l=lnd%ls,lnd%le
     if(lnd%ug_landfrac(l)>0.neqv.ANY(land2cplr%mask(l,:))) then
        call error_mesg('land_model_init','land masks from grid spec and from land restart do not match',FATAL)
     endif
  enddo

  ! [9] check the properties of co2 exchange with the atmosphere and set appropriate
  ! flags
  if (canopy_air_mass_for_tracers==0.and.ico2==NO_TRACER) then
     call error_mesg('land_model_init', &
          'canopy_air_mass_for_tracers is set to zero, and CO2 exchange with the atmosphere is not set up: '// &
          'canopy air CO2 concentration will not be updated',NOTE)
     update_cana_co2 = .FALSE.
  else
     update_cana_co2 = .TRUE.
  end if

  if (reset_to_ntrl) then
     call error_mesg('land_model_init','removing non-primary vegetation',NOTE)
     call remove_non_primary()
  endif

  call mpp_clock_end(landInitClock)

end subroutine land_model_init

! ============================================================================
! removes non-primary vegetation tiles, rescaling the fractions
subroutine remove_non_primary()
  type(land_tile_enum_type):: ce,te,next
  integer :: l, k, n_ntrl
  type(land_tile_type), pointer :: tile, tile_max
  real :: f_ntrl, f_soil, btot, btot_max

  do l = lnd%ls, lnd%le
     call set_current_point(l,1)
     f_ntrl = 0.0; f_soil = 0.0; n_ntrl=0
     ce = first_elmt(land_tile_map(l))
     do while (loop_over_tiles(ce, tile, k=k))
        if(is_watch_cell()) then
           __DEBUG2__(k,tile%frac)
        endif
        if (.not.associated(tile%vegn)) cycle
        if (tile%vegn%landuse==LU_NTRL) then
            f_ntrl = f_ntrl+tile%frac; n_ntrl = n_ntrl+1
        endif
        f_soil = f_soil + tile%frac
        if (is_watch_cell()) then
           __DEBUG3__(k,tile%vegn%landuse,tile%frac)
        endif
     enddo

     if (is_watch_cell()) then
        __DEBUG2__(f_soil, f_ntrl)
     endif
     if (f_soil==0) cycle ! do nothing for grid cells that do not have soil/vegetation tiles

     if (f_ntrl==0) then
        ! if natural area in grid cell is zero, find tile with highest biomass and designate
        ! it as natural
        call land_error_message('remove_non_primary: fraction of NTRL is zero; looking for tile with highest biomass. f_soil='//&
                                             trim(string(f_soil))//' n_ntrl='//string(n_ntrl), WARNING)
        ce = first_elmt(land_tile_map(l)); btot_max = -HUGE(1.0); tile_max => NULL()
        do while (loop_over_tiles(ce, tile, k=k))
           if (.not.associated(tile%vegn)) cycle
           btot = sum(plant_C(tile%vegn%cohorts(1:tile%vegn%n_cohorts)))
           __DEBUG4__(k,tile%frac,tile%vegn%landuse,btot)
           if (btot>btot_max .and. tile%frac>0) then
              btot_max = btot
              tile_max => tile
           endif
        enddo
        if(.not.associated(tile_max)) &
            call land_error_message('remove_non_primary: could not find vegetation tile with max biomass', FATAL)
        __DEBUG3__(tile_max%frac,tile_max%vegn%landuse,btot_max)
        ! change land use type in tile with max biomass to natural
        tile_max%vegn%landuse = LU_NTRL
        ! reconcile total fraction of NTRL in grid cell with land use type changes
        f_ntrl = tile_max%frac
     endif

     ! remove all non-natural tiles, and scale fractions of natural tiles
     ce = first_elmt(land_tile_map(l))
     te = tail_elmt(land_tile_map(l))
     do
        if(ce==te) exit ! reached the end of the list
        tile=>current_tile(ce)
        next = next_elmt(ce)
        if (associated(tile%vegn)) then
           if (tile%vegn%landuse==LU_NTRL) then
              tile%frac = tile%frac*f_soil/f_ntrl
           else
              call remove(ce)
           endif
        endif
        ce = next
     enddo
  enddo
end subroutine remove_non_primary

! ============================================================================
subroutine land_model_end (cplr2land, land2cplr)
  type(atmos_land_boundary_type), intent(inout) :: cplr2land
  type(land_data_type)          , intent(inout) :: land2cplr

  ! ---- local vars

  module_is_initialized = .FALSE.

  call error_mesg('land_model_end','writing NetCDF restart',NOTE)
  call land_model_restart()

  ! we still want to call the *_end procedures for component models, even
  ! if the number of tiles in this domain is zero, in case they are doing
  ! something else besides saving the restart, of if they want to save
  ! restart anyway
  call land_tracer_driver_end()
  call land_transitions_end()
  call glac_end ()
  call lake_end ()
  call soil_end ()
  call hlsp_end ()
  call snow_end ()
  call vegn_end ()
  call cana_end ()
  call nitrogen_sources_end ()
  call topo_rough_end()
  call river_end()

  call dealloc_land2cplr(land2cplr, dealloc_discharges=.TRUE.)
  call dealloc_cplr2land(cplr2land)

  call land_tracers_end()

  call tile_diag_end()

  ! deallocate tiles
  call free_tile_map()
  call land_data_end()

  ! finish up the land debugging diagnostics
  call land_debug_end()

end subroutine land_model_end


! ============================================================================
subroutine land_model_restart(timestamp)
  character(*), intent(in), optional :: timestamp ! timestamp to add to the file name

  ! ---- local vars
  character(267) :: filename
  type(land_restart_type) :: restart ! restart file i/o object
  integer :: tile_dim_length ! length of tile dimension in output files
                             ! global max of number of tiles per gridcell
  integer :: l
  character(256) :: timestamp_

  ! [1] count all land tiles and determine the length of tile dimension
  ! sufficient for the current domain
  tile_dim_length = 0
  do l = lnd%ls, lnd%le
     tile_dim_length = max(tile_dim_length,nitems(land_tile_map(l)))
  enddo

  ! [2] calculate the tile dimension length by taking the max across all domains
  call mpp_max(tile_dim_length)
  if (tile_dim_length==0) then
     call error_mesg('land_model_restart',&
       'No land points exist (tile_dim_length=0), therefore no land restarts will be saved',&
       WARNING)
     return
  endif

  ! [3] create tile output file
  timestamp_=''
  if (present(timestamp)) then
     if(trim(timestamp)/='') timestamp_=trim(timestamp)//'.'
  endif
  ! Note that filename is updated for tile & rank numbers during file creation
  filename = trim(timestamp)//'land.res.nc'
  call init_land_restart(restart, filename, land_tile_exists, tile_dim_length)

  ! [4] write data fields
  ! write fractions and tile tags
  call add_tile_data(restart,'frac',land_frac_ptr,'fractional area of tile')
  call add_int_tile_data(restart,'glac',glac_tag_ptr,'tag of glacier tiles')
  call add_int_tile_data(restart,'lake',lake_tag_ptr,'tag of lake tiles')
  call add_int_tile_data(restart,'soil',soil_tag_ptr,'tag of soil tiles')
  call add_int_tile_data(restart,'vegn',vegn_tag_ptr,'tag of vegetation tiles')
  ! write the upward long-wave flux
  call add_tile_data(restart,'lwup',land_lwup_ptr,'upward long-wave flux')
  ! write energy residuals
  call add_tile_data(restart,'e_res_1',land_e_res_1_ptr,&
       'energy residual in canopy air energy balance equation', 'W/m2')
  call add_tile_data(restart,'e_res_2',land_e_res_2_ptr,&
       'energy residual in canopy energy balance equation', 'W/m2')

  ! [5] close file
  call save_land_restart(restart)
  call free_land_restart(restart)

  ! [6] save component model restarts
  call save_land_transitions_restart(timestamp_)
  call save_glac_restart(tile_dim_length,timestamp_)
  call save_lake_restart(tile_dim_length,timestamp_)
  call save_soil_restart(tile_dim_length,timestamp_)
  call save_hlsp_restart(tile_dim_length,timestamp_)
  call save_snow_restart(tile_dim_length,timestamp_)
  call save_vegn_restart(tile_dim_length,timestamp_)
  call save_cana_restart(tile_dim_length,timestamp_)
  call save_fire_restart(tile_dim_length,timestamp_)
  call save_river_restart(timestamp_)

end subroutine land_model_restart

! ============================================================================
subroutine land_cover_cold_start()
  real, dimension(:,:), pointer :: &
       glac, soil, lake, vegn ! arrays of fractions for respective sub-models
  integer, pointer, dimension(:,:) :: soiltags ! array of soil type tags
  integer, pointer, dimension(:,:) :: hlsp_pos ! hillslope position index
  integer, pointer, dimension(:,:) :: hlsp_par ! hillslope parent index
  real   , allocatable :: rbuffer(:,:) ! real buffer for remap
  logical, dimension(lnd%le-lnd%ls+1) :: &
       land_mask, valid_data, invalid_data
  integer :: iwatch,jwatch,kwatch,lwatch,face
  integer :: l,p,lll
  integer :: tile_root_pe,tile_npes
  integer :: ps,pe ! boundaries of PE list for remapping
  type(horiz_remap_type) :: map

  ! calculate the global land mask
  land_mask = lnd%ug_area > 0

  ! get the global maps of fractional covers for each of the sub-models
  glac=>glac_cover_cold_start(land_mask,lnd%sg_lonb,lnd%sg_latb)
  lake=>lake_cover_cold_start(land_mask,lnd%sg_lonb,lnd%sg_latb,lnd%sg_domain)
  soil=>soil_cover_cold_start(land_mask,lnd%sg_lonb,lnd%sg_latb)
  vegn=>vegn_cover_cold_start(land_mask,lnd%sg_lonb,lnd%sg_latb)

  ! Because of hillslope model, soil tiles may not be returned in order of soil type.
  allocate(soiltags(size(soil,1), size(soil,2)))
  call retrieve_soil_tags(soiltags)
  ! Tiles will be constructed with hillslope data.
  allocate(hlsp_pos(size(soil,1), size(soil,2)))
  allocate(hlsp_par(size(soil,1), size(soil,2)))
  call retrieve_hlsp_indices(hlsp_pos, hlsp_par)

  ! remove any input lake fraction in coastal cells
  where (lnd%ug_landfrac.lt. 1.-gfrac_tol) lake(:,1) = 0.
  ! NOTE that the lake area in the coastal cells can be set to non-zero
  ! again by the "ground fraction reconciliation code" below. Strictly
  ! speaking the above line of code should be replaced with the section
  ! commented out with "!-zero" below, but we preserve the old way to avoid
  ! backward incompatibility with older runs. This needs updating in the
  ! future when the decision about what to do with lakes in coastal cells is
  ! made.

  ! reconcile ground fractions with the land mask within compute domain
  valid_data = land_mask.and.(sum(glac,2)+sum(lake,2)+sum(soil,2)>0)
  invalid_data = land_mask.and..not.valid_data

  call get_watch_point(iwatch,jwatch,kwatch,face,lwatch)
  if (face==lnd%ug_face.and.(lnd%ls<=lwatch.and.lwatch<=lnd%le) ) then
     write(*,*)'###### land_cover_cold_start: input data #####'
     write(*,'(99(a,i4.2,x))')'iwatch=',iwatch,'jwatch=',jwatch,'face=',face
     write(*,'(99(a,g23.16,x))')'lon=',lnd%ug_lon(lwatch)*180/PI,'lat=',lnd%ug_lat(lwatch)*180/PI
     ! calculate local compute domain indices; we assume glac,lake,soil,vegn all
     ! have the same lbounds
     l = lwatch-lnd%ls+lbound(glac,1)
     __DEBUG1__(lnd%ls)
     write(*,'(a,99(a,i4.2,x))')'local indices:','l=',l
     __DEBUG3__(lnd%ug_landfrac(lwatch),land_mask(l),valid_data(l))
     __DEBUG1__(glac(l,:))
     __DEBUG1__(lake(l,:))
     __DEBUG1__(soil(l,:))
     __DEBUG1__(vegn(l,:))
  endif

  if (trim(nearest_point_search)=='global') then
     ps=0 ; pe=size(lnd%pelist)-1
  else if (trim(nearest_point_search)=='face') then
     call mpp_get_UG_domain_tile_pe_inf(lnd%ug_domain, tile_root_pe, npes=tile_npes)
     ps = -1
     do p = 0, size(lnd%pelist(:))-1
        if(lnd%pelist(p) == tile_root_pe) then
           ps = p
           exit
        endif
     enddo
     if( ps == -1) call error_mesg('land_cover_cold_start',&
           'tile_root_pe is not in the lnd%pelist', FATAL)
     pe = ps + tile_npes-1
  else
     call error_mesg('land_cover_cold_start',&
          'option nearest_point_search="'//trim(nearest_point_search)//&
          '" is illegal, use "global" or "face"',&
          FATAL)
  endif
  call horiz_remap_new(invalid_data,valid_data,lnd%ug_lon,lnd%ug_lat,lnd%ug_domain,&
          lnd%pelist(ps:pe),map)
  if (print_remapping) call horiz_remap_print(map,'land cover remap:')
  call horiz_remap(map,lnd%ug_domain,glac)
  call horiz_remap(map,lnd%ug_domain,lake)
  call horiz_remap(map,lnd%ug_domain,soil)
  allocate(rbuffer(size(soil,1), size(soil,2)))
  ! ZMS This is awkward: perhaps horiz_remap should handle integers?
  rbuffer = real(soiltags); call horiz_remap(map,lnd%ug_domain,rbuffer); soiltags = nint(rbuffer)
  rbuffer = real(hlsp_pos); call horiz_remap(map,lnd%ug_domain,rbuffer); hlsp_pos = nint(rbuffer)
  rbuffer = real(hlsp_par); call horiz_remap(map,lnd%ug_domain,rbuffer); hlsp_par = nint(rbuffer)
  deallocate(rbuffer)
  call horiz_remap_del(map)

!-zero  ! remove any input lake fraction in coastal cells
!-zero  do j = lnd%js,lnd%je
!-zero  do i = lnd%is,lnd%ie
!-zero     call set_current_point(i,j,1)
!-zero     if (lnd%landfrac(i,j) < 1-gfrac_tol) then
!-zero        lake(i,j,:) = 0.0
!-zero        if(is_watch_point())then
!-zero           write(*,*)'###### land_cover_cold_start: lake fraction is set to zero #####'
!-zero        endif
!-zero     endif
!-zero  enddo
!-zero  enddo

  ! reconcile vegetation fractions with the land mask within compute domain
  valid_data = sum(vegn,2) > 0
  invalid_data = .FALSE.
  do l = 1,size(land_mask(:))
     if(.not.land_mask(l)) cycle ! skip ocean points
     if(valid_data(l)) cycle ! do not need to do anything with valid points
     if(sum(glac(l,:))+sum(lake(l,:))>=1) &
          cycle                ! skip points fully covered by glaciers or lakes
     invalid_data(l)=.TRUE.
  enddo
  call horiz_remap_new(invalid_data,valid_data,lnd%ug_lon,lnd%ug_lat,lnd%ug_domain,&
       lnd%pelist(ps:pe),map)
  if (print_remapping) call horiz_remap_print(map,'vegetation cover remap:')
  call horiz_remap(map,lnd%ug_domain,vegn)
  call horiz_remap_del(map)

  ! create tiles
  do l = 1,size(land_mask(:))
     lll = l+lnd%ls-1
     if(.not.land_mask(l)) cycle ! skip ocean points
     call set_current_point(lll,1)
     call land_cover_cold_start_0d &
          (land_tile_map(lll),glac(l,:),lake(l,:),soil(l,:),soiltags(l,:),&
               hlsp_pos(l,:), hlsp_par(l,:), vegn(l,:))
     if(nitems(land_tile_map(lll))==0) then
        call error_mesg('land_cover_cold_start',&
             'No tiles were created for a valid land point at i='&
             //trim(string(lnd%i_index(lll)))//' j='//trim(string(lnd%j_index(lll)))&
             //' face='//trim(string(lnd%ug_face)), FATAL)
     endif
  enddo

  deallocate(glac,lake,soil,soiltags,hlsp_pos,hlsp_par,vegn)

end subroutine land_cover_cold_start

! ============================================================================
subroutine land_cover_cold_start_0d (set,glac0,lake0,soil0,soiltags0,&
                                     hlsp_pos0,hlsp_par0,vegn0)
  type(land_tile_list_type), intent(inout) :: set
  real, dimension(:)       , intent(in) :: &
       glac0,lake0,soil0,vegn0 ! fractions of area
  integer, dimension(:)    , intent(in) :: &
       soiltags0, hlsp_pos0, hlsp_par0 ! soil and hillslope tags

  ! ---- local vars
  real :: glac(size(glac0(:))), lake(size(lake0(:))), &
          soil(size(soil0(:))), vegn(size(vegn0(:)))
  type(land_tile_type), pointer :: tile
  integer :: i,j
  real :: factor ! normalizing factor for the tile areas
  real :: frac
  type(land_tile_enum_type) :: first_non_vegn ! position of first non-vegetated tile in the list
  type(land_tile_enum_type) :: ce
  real :: athresh = 1.e-10 ! area threshold allowable for deviation from 1
  real :: area ! area sum for gridcell

  glac = glac0; lake = lake0; soil = soil0; vegn = vegn0
  if (sum(glac)>1) &
       glac=glac/sum(glac)
  if (sum(lake)+sum(glac)>1)&
       lake = lake*(1-sum(glac))/sum(lake)
  if (sum(lake)<min_sum_lake_frac) lake=0
  if (sum(soil)+sum(glac)+sum(lake)>1)&
       soil = soil*(1-sum(lake)-sum(glac))/sum(soil)
  ! make sure that the sum of the fractions of the soil, lake, and glaciers are
  ! either one or zero
  factor = sum(soil)+sum(glac)+sum(lake)
  if(factor>0)then
     glac = glac/factor
     lake = lake/factor
     soil = soil/factor
  endif

  ! remove soil/glac/lake fractions that are too small
  if (min_frac>0) then
     where (glac<min_frac) glac = 0
     where (lake<min_frac) lake = 0
     where (soil<min_frac) soil = 0
     ! do the renormalization again
     factor = sum(soil)+sum(glac)+sum(lake)
     if(factor>0)then
        glac = glac/factor
        lake = lake/factor
        soil = soil/factor
     endif
  endif

  if(is_watch_point()) then
     write(*,*)'#### land_cover_cold_start_0d input data ####'
     __DEBUG1__(glac0)
     __DEBUG1__(lake0)
     __DEBUG1__(soil0)
     __DEBUG1__(vegn0)
     __DEBUG1__(factor)
     write(*,*)'#### land_cover_cold_start_0d renormalized fractions ####'
     __DEBUG1__(glac)
     __DEBUG1__(lake)
     __DEBUG1__(soil)
     __DEBUG1__(vegn)
  endif

  do i = 1,size(glac)
     if (glac(i)>0) then
        tile => new_land_tile(frac=glac(i),glac=i)
        call insert(tile,set)
        if(is_watch_point()) then
           write(*,*)'created glac tile: frac=',glac(i),' tag=',i
        endif
     endif
  enddo
  do i = 1,size(lake)
     if (lake(i)>0) then
        tile => new_land_tile(frac=lake(i),lake=i)
        call insert(tile,set)
        if(is_watch_point()) then
           write(*,*)'created lake tile: frac=',lake(i),' tag=',i
        endif
     endif
  enddo

  factor = sum(soil)*sum(vegn)
  if (factor/=0) factor = 1/factor
  factor = factor*(1-sum(glac)-sum(lake))
  ! vegetation tiles, if any, are inserted in front of non-vegetated tiles;
  ! this really does not matter except for the static vegetation override
  ! case with the data saved by lm3v -- there the vegetation tiles are
  ! in front, so it works more consistently where lad2 has more than
  ! one tile (e.g. glac/soil or lake/soil), if lad2 vegetation tiles are
  ! also in front of the list.
  first_non_vegn=first_elmt(set)
  do i = 1,size(soil)
  do j = 1,size(vegn)
     frac = soil(i)*vegn(j)*factor
     if(frac>0) then
        tile  => new_land_tile(frac=frac,soil=soiltags0(i),vegn=j,&
                               htag_j=hlsp_pos0(i),htag_k=hlsp_par0(i))
        call insert(tile,first_non_vegn)
        if(is_watch_point()) then
           write(*,*)'created soil tile: frac=', frac, ' soil tag=',soiltags0(i), ' veg tag=',j
        endif
     endif
  enddo
  enddo

  ! Check that fractions are set correctly.
  ce = first_elmt(set)
  area = 0.
  do while (loop_over_tiles(ce,tile))
     area = area + tile%frac
  end do

  if (abs(area - 1.) > athresh) then
     call error_mesg(module_name, 'Area fractions do not add up to 1 for tile list in land_cover_cold_start_0d!', &
                     FATAL)
  end if

end subroutine land_cover_cold_start_0d

! ============================================================================
subroutine land_cover_warm_start(restart)
  type(land_restart_type), intent(in) :: restart
  if (new_land_io) then
     call land_cover_warm_start_new(restart)
  else
     call land_cover_warm_start_orig(restart)
  endif
end subroutine land_cover_warm_start

! ============================================================================
! reads the land restart file and restores the tiling structure from this file
subroutine land_cover_warm_start_new (restart)
  type(land_restart_type), intent(in) :: restart

  ! ---- local vars
  integer, allocatable :: glac(:), lake(:), soil(:), vegn(:) ! tile tags
  real,    allocatable :: frac(:) ! fraction of land covered by tile
  integer :: ntiles    ! total number of land tiles in the input file
  integer :: k,it,l,g, npts
  type(land_tile_type), pointer :: tile

  ntiles = size(restart%tidx)
  allocate(glac(ntiles), lake(ntiles), soil(ntiles), vegn(ntiles), frac(ntiles))

  call fms_io_unstructured_read(restart%basename, "frac", frac, lnd%ug_domain, timelevel=1)
  call fms_io_unstructured_read(restart%basename, "glac", glac, lnd%ug_domain, timelevel=1)
  call fms_io_unstructured_read(restart%basename, "lake", lake, lnd%ug_domain, timelevel=1)
  call fms_io_unstructured_read(restart%basename, "soil", soil, lnd%ug_domain, timelevel=1)
  call fms_io_unstructured_read(restart%basename, "vegn", vegn, lnd%ug_domain, timelevel=1)

  npts = lnd%nlon*lnd%nlat
  ! create tiles
  do it = 1,ntiles
     k = restart%tidx(it)
     if (k<0) cycle ! skip negative indices
     g = modulo(k,npts)+1
     if (g<lnd%gs.or.g>lnd%ge) cycle ! skip points outside of domain
     l = lnd%l_index(g)
     ! the size of the tile set at the point (i,j) must be equal to k
     tile=>new_land_tile(frac=frac(it),&
              glac=glac(it),lake=lake(it),soil=soil(it),vegn=vegn(it))
     call insert(tile,land_tile_map(l))
  enddo
  deallocate(glac, lake, soil, vegn, frac)
end subroutine land_cover_warm_start_new


! ============================================================================
! reads the land restart file and restores the tiling structure from this file
subroutine land_cover_warm_start_orig (restart)
  type(land_restart_type), intent(in) :: restart

  ! ---- local vars
  integer, allocatable :: idx(:) ! compressed tile index
  integer, allocatable :: glac(:), lake(:), soil(:), snow(:), cana(:), vegn(:) ! tile tags
  real,    allocatable :: frac(:) ! fraction of land covered by tile
  integer :: ncid ! unit number of the input file
  integer :: ntiles    ! total number of land tiles in the input file
  integer :: bufsize   ! size of the input buffer
  integer :: dimids(1) ! id of tile dimension
  character(NF_MAX_NAME) :: tile_dim_name ! name of the tile dimension and respective variable
  integer :: k,it,npts,g,l
  type(land_tile_type), pointer :: tile;
  integer :: start, count ! slab for reading
  ! netcdf variable IDs
  integer :: id_idx, id_frac, id_glac, id_lake, id_soil, id_vegn

  __NF_ASRT__(nf_open(restart%filename,NF_NOWRITE,ncid))
  ! allocate the input data
  __NF_ASRT__(nfu_inq_var(ncid,'frac',id=id_frac,varsize=ntiles,dimids=dimids))
   ! allocate input buffers for compression index and the variable
  bufsize=min(input_buf_size,ntiles)
  allocate(idx (bufsize), glac(bufsize), lake(bufsize), soil(bufsize), &
           snow(bufsize), cana(bufsize), vegn(bufsize), frac(bufsize)  )
  ! get the name of the fist (and only) dimension of the variable 'frac' -- this
  ! is supposed to be the compressed dimension, and associated variable will
  ! hold the compressed indices
  __NF_ASRT__(nfu_inq_dim(ncid,dimids(1),name=tile_dim_name))
  __NF_ASRT__(nfu_inq_var(ncid,tile_dim_name,id=id_idx))
  ! get the IDs of the variables to read
  __NF_ASRT__(nfu_inq_var(ncid,'glac',id=id_glac))
  __NF_ASRT__(nfu_inq_var(ncid,'lake',id=id_lake))
  __NF_ASRT__(nfu_inq_var(ncid,'soil',id=id_soil))
  __NF_ASRT__(nfu_inq_var(ncid,'vegn',id=id_vegn))

  npts = lnd%nlon*lnd%nlat
  do start = 1,ntiles,bufsize
    count = min(bufsize,ntiles-start+1)
    ! read the compressed tile indices
    __NF_ASRT__(nf_get_vara_int(ncid,id_idx,(/start/),(/count/),idx))
    ! read input data -- fractions and tags
    __NF_ASRT__(nf_get_vara_double(ncid,id_frac,(/start,1/),(/count,1/),frac))
    __NF_ASRT__(nf_get_vara_int(ncid,id_glac,(/start,1/),(/count,1/),glac))
    __NF_ASRT__(nf_get_vara_int(ncid,id_lake,(/start,1/),(/count,1/),lake))
    __NF_ASRT__(nf_get_vara_int(ncid,id_soil,(/start,1/),(/count,1/),soil))
    __NF_ASRT__(nf_get_vara_int(ncid,id_vegn,(/start,1/),(/count,1/),vegn))
    ! create tiles
    do it = 1,count
       k = idx(it)
       if (k<0) cycle ! skip negative indices
       g = modulo(k,npts)+1
       if (g<lnd%gs.or.g>lnd%ge) cycle ! skip points outside of domain
       l = lnd%l_index(g)
       ! the size of the tile set at the point (i,j) must be equal to k
       tile=>new_land_tile(frac=frac(it),&
                glac=glac(it),lake=lake(it),soil=soil(it),vegn=vegn(it))
       call insert(tile,land_tile_map(l))
    enddo
  enddo
  __NF_ASRT__(nf_close(ncid))
  deallocate(idx, glac, lake, soil, snow, cana, vegn, frac)
end subroutine land_cover_warm_start_orig


! ============================================================================
subroutine update_land_model_fast ( cplr2land, land2cplr )
  type(atmos_land_boundary_type), intent(in)    :: cplr2land
  type(land_data_type)          , intent(inout) :: land2cplr

  ! ---- local vars
  real :: &
       ISa_dn_dir(NBANDS), & ! downward direct sw radiation at the top of the canopy
       ISa_dn_dif(NBANDS)    ! downward diffuse sw radiation at the top of the canopy

  ! variables for stock calculations
  real :: &
     cana_VMASS, cana_HEAT,             &
     vegn_LMASS, vegn_FMASS, vegn_HEAT, &
     snow_LMASS, snow_FMASS, snow_HEAT, &
     subs_LMASS, subs_FMASS, subs_HEAT, &
     glac_LMASS, glac_FMASS, glac_HEAT, &
     lake_LMASS, lake_FMASS, lake_HEAT, &
     soil_LMASS, soil_FMASS, soil_HEAT

  real :: runoff(lnd%ls:lnd%le) ! total (liquid+snow) runoff accumulated over tiles in cell
  real :: runoff_c(lnd%ls:lnd%le,n_river_tracers)
                 ! runoff of tracers accumulated over tiles in cell (including ice and heat)
  real :: runoff_sg(lnd%is:lnd%ie,lnd%js:lnd%je)
  real :: runoff_c_sg(lnd%is:lnd%ie,lnd%js:lnd%je,n_river_tracers)

  real :: snc(lnd%ls:lnd%le), snow_depth, snow_area ! snow cover, for CMIP6 diagnostics

  logical :: used          ! return value of send_data diagnostics routine
  integer :: i,j,k,l   ! lon, lat, and tile indices
  integer :: is,ie,js,je ! horizontal bounds of the override buffer
  type(land_tile_enum_type) :: ce ! tile enumerator
  type(land_tile_type), pointer :: tile ! pointer to current tile
  integer :: n_cohorts ! number of cohorts per tile
  integer :: iwatch,jwatch,kwatch,face

  ! variables for data override
  real, allocatable :: phot_co2_data(:)  ! buffer for data
  logical           :: phot_co2_overridden ! flag indicating successful override

  ! variables for total water storage diagnostics
  real :: twsr_sg(lnd%is:lnd%ie,lnd%js:lnd%je), tws(lnd%ls:lnd%le)

  ! start clocks
  call mpp_clock_begin(landClock)
  call mpp_clock_begin(landFastClock)

  ! to avoid output of static vegetation after the transitions worked and
  ! changed the tiling structure, static vegetation output is done here.
  call write_static_vegn()

  ! override data at the beginning of the time step
  is= lnd%is ; ie = lnd%ie
  js= lnd%js ; je = lnd%je
  allocate(phot_co2_data(lnd%ls:lnd%le))
  call data_override_ug('LND','phot_co2',phot_co2_data,lnd%time, &
       override=phot_co2_overridden)

  ! get the fertilization data
  call update_nitrogen_sources(lnd%time, lnd%time+lnd%dt_fast)

  ! clear the runoff values, for accumulation over the tiles
  runoff = 0 ; runoff_c = 0
  ! clear the snc (snow cover diag) values, for accumulation over the tiles
  snc = 0

  ! Calculate groundwater and associated heat fluxes between tiles within each gridcell.
  call hlsp_hydrology_1(n_c_types)
  ! ZMS: Eventually pass these args into river or main tile loop.

  ! main tile loop
!$OMP parallel do default(none) shared(lnd,land_tile_map,cplr2land,land2cplr,phot_co2_overridden, &
!$OMP                                  phot_co2_data,runoff,runoff_c,snc,id_area,id_z0m,id_z0s,       &
!$OMP                                  id_Trad,id_Tca,id_qca,isphum,id_cd_m,id_cd_t,id_snc) &
!$OMP                                  private(i,j,k,ce,tile,ISa_dn_dir,ISa_dn_dif,n_cohorts,snow_depth,snow_area)
  do l = lnd%ls, lnd%le
     i = lnd%i_index(l)
     j = lnd%j_index(l)
!     __DEBUG4__(is,js,i-is+lnd%is,j-js+lnd%js)
     ce = first_elmt(land_tile_map(l))
     do while (loop_over_tiles(ce,tile,k=k))
        ! set this point coordinates as current for debug output
        call set_current_point(i,j,k,l)

        ISa_dn_dir(BAND_VIS) = cplr2land%sw_flux_down_vis_dir(l,k)
        ISa_dn_dir(BAND_NIR) = cplr2land%sw_flux_down_total_dir(l,k)&
                              -cplr2land%sw_flux_down_vis_dir(l,k)
        ISa_dn_dif(BAND_VIS) = cplr2land%sw_flux_down_vis_dif(l,k)
        ISa_dn_dif(BAND_NIR) = cplr2land%sw_flux_down_total_dif(l,k)&
                              -cplr2land%sw_flux_down_vis_dif(l,k)

        ! n_cohorts is calculated and passed down to update_land_model_fast_0d
        ! for convenience, so that there is no need to make a lot of by-cohort arrays
        ! allocatable -- instead they are created on stack with the size passed
        ! as an argument.
        if (associated(tile%vegn)) then
           n_cohorts = tile%vegn%n_cohorts
        else
           n_cohorts = 1
        endif

        call update_land_model_fast_0d(tile, l,k, n_cohorts, land2cplr, &
           cplr2land%lprec(l,k),  cplr2land%fprec(l,k), cplr2land%tprec(l,k), &
           cplr2land%wind(l,k), &
           cplr2land%t_flux(l,k), cplr2land%dhdt(l,k), &
           cplr2land%tr_flux(l,k,:), cplr2land%dfdtr(l,k,:), &
           ISa_dn_dir, ISa_dn_dif, cplr2land%lwdn_flux(l,k), &
           cplr2land%ustar(l,k), cplr2land%p_surf(l,k), cplr2land%drag_q(l,k), &
           phot_co2_overridden, phot_co2_data(l),&
           runoff(l), runoff_c(l,:) &
        )
        ! some of the diagnostic variables are sent from here, purely for coding
        ! convenience: the compute domain-level 2d and 3d vars are generally not
        ! available inside update_land_model_fast_0d, so the diagnostics for those
        ! was left here.
        call send_tile_data(id_area, tile%frac*lnd%ug_area(l),     tile%diag)
        call send_tile_data(id_z0m,  land2cplr%rough_mom(l,k),     tile%diag)
        call send_tile_data(id_z0s,  land2cplr%rough_heat(l,k),    tile%diag)
        call send_tile_data(id_Trad, land2cplr%t_surf(l,k),        tile%diag)
        call send_tile_data(id_Tca,  land2cplr%t_ca(l,k),          tile%diag)
        call send_tile_data(id_qca,  land2cplr%tr(l,k,isphum),     tile%diag)
        call send_tile_data(id_cd_m, cplr2land%cd_m(l,k),          tile%diag)
        call send_tile_data(id_cd_t, cplr2land%cd_t(l,k),          tile%diag)

        if (id_snc>0) then
           call snow_get_depth_area ( tile%snow, snow_depth, snow_area )
           snc(l) = snc(l) + snow_area*tile%frac
        endif
     enddo
  enddo

  !--- pass runoff from unstructured grid to structured grid.
  runoff_sg = 0 ; runoff_c_sg = 0
  call mpp_pass_UG_to_SG(lnd%ug_domain, runoff,   runoff_sg  )
  call mpp_pass_UG_to_SG(lnd%ug_domain, runoff_c, runoff_c_sg)

  call get_watch_point(iwatch,jwatch,kwatch,face)
  if (face==lnd%sg_face.and.(lnd%is<=iwatch.and.iwatch<=lnd%ie).and.&
                            (lnd%js<=jwatch.and.jwatch<=lnd%je).and.&
                            is_watch_time()) then
     __DEBUG1__(runoff_sg(iwatch,jwatch))
     __DEBUG1__(runoff_c_sg(iwatch,jwatch,:))
  endif

  !--- update river state
  call update_river(runoff_sg, runoff_c_sg, land2cplr)

  if(id_tws>0) then
     call get_river_water(twsr_sg)
     call mpp_pass_SG_to_UG(lnd%ug_domain, twsr_sg, tws)
  endif

  ce = first_elmt(land_tile_map, ls=lbound(cplr2land%t_flux,1) )
  do while(loop_over_tiles(ce,tile,l,k))
     cana_VMASS = 0. ;                   cana_HEAT = 0.
     vegn_LMASS = 0. ; vegn_FMASS = 0. ; vegn_HEAT = 0.
     snow_LMASS = 0. ; snow_FMASS = 0. ; snow_HEAT = 0.
     subs_LMASS = 0. ; subs_FMASS = 0. ; subs_HEAT = 0.
     glac_LMASS = 0. ; glac_FMASS = 0. ; glac_HEAT = 0.
     lake_LMASS = 0. ; lake_FMASS = 0. ; lake_HEAT = 0.
     soil_LMASS = 0. ; soil_FMASS = 0. ; soil_HEAT = 0.
     if (associated(tile%cana)) then
         cana_VMASS = canopy_air_mass*tile%cana%tr(isphum)
         cana_HEAT  = cana_tile_heat(tile%cana)
     endif
     if (associated(tile%vegn)) then
         call vegn_tile_stock_pe(tile%vegn, vegn_LMASS, vegn_FMASS)
         vegn_HEAT = vegn_tile_heat(tile%vegn)
     endif
     if(associated(tile%snow)) then
         call snow_tile_stock_pe(tile%snow, snow_LMASS, snow_FMASS)
         snow_HEAT = snow_tile_heat(tile%snow)
     endif
     if (associated(tile%glac)) then
         call glac_tile_stock_pe(tile%glac, subs_LMASS, subs_FMASS)
         subs_HEAT  = glac_tile_heat(tile%glac)
         glac_LMASS = subs_LMASS
         glac_FMASS = subs_FMASS
         glac_HEAT  = subs_HEAT
     else if (associated(tile%lake)) then
         call lake_tile_stock_pe(tile%lake, subs_LMASS, subs_FMASS)
         subs_HEAT  = lake_tile_heat(tile%lake)
         lake_LMASS = subs_LMASS
         lake_FMASS = subs_FMASS
         lake_HEAT  = subs_HEAT
     else if (associated(tile%soil)) then
         call soil_tile_stock_pe(tile%soil, subs_LMASS, subs_FMASS)
         subs_HEAT  = soil_tile_heat(tile%soil)
         soil_LMASS = subs_LMASS
         soil_FMASS = subs_FMASS
         soil_HEAT  = subs_HEAT
     endif

     call send_tile_data(id_VWS,  cana_VMASS, tile%diag)
     call send_tile_data(id_VWSc, cana_VMASS, tile%diag)
     call send_tile_data(id_LWS,  vegn_LMASS+snow_LMASS+subs_LMASS, tile%diag)
     call send_tile_data(id_LWSv, vegn_LMASS, tile%diag)
     call send_tile_data(id_LWSs, snow_LMASS, tile%diag)
     call send_tile_data(id_LWSg, subs_LMASS, tile%diag)
     call send_tile_data(id_FWS,  vegn_FMASS+snow_FMASS+subs_FMASS, tile%diag)
     call send_tile_data(id_FWSv, vegn_FMASS, tile%diag)
     call send_tile_data(id_FWSs, snow_FMASS, tile%diag)
     call send_tile_data(id_FWSg, subs_FMASS, tile%diag)
     call send_tile_data(id_HS,  vegn_HEAT+snow_HEAT+subs_HEAT+cana_HEAT, tile%diag)
     call send_tile_data(id_HSv, vegn_HEAT, tile%diag)
     call send_tile_data(id_HSs, snow_HEAT, tile%diag)
     call send_tile_data(id_HSg, subs_HEAT, tile%diag)
     call send_tile_data(id_HSc, cana_HEAT, tile%diag)
     call send_tile_data(id_water, subs_LMASS+subs_FMASS, tile%diag)
     call send_tile_data(id_snow,  snow_LMASS+snow_FMASS, tile%diag)

     ! CMOR variables
     call send_tile_data(id_snw, snow_FMASS, tile%diag)
     call send_tile_data(id_lwsnl, snow_LMASS, tile%diag)
     ! factor 1000.0 kg/m3 is the liquid water density; it converts mass of water into depth
     call send_tile_data(id_sweLut, max(snow_FMASS+snow_LMASS,0.0)/1000.0, tile%diag)
     if (id_tws>0) then
         ! note that subs_LMASS and subs_FMASS are reused here to hold total water masses
         ! alos note that reported value does not include river storage
         call get_tile_water(tile, subs_LMASS, subs_FMASS)
         tws(l) = tws(l) + (subs_LMASS+subs_FMASS)*tile%frac
     endif
  enddo
  if (id_tws>0) used = send_data(id_tws, tws, lnd%time)

  ! advance land model time
  lnd%time = lnd%time + lnd%dt_fast

  ! send the accumulated diagnostics to the output
  call dump_tile_diag_fields(lnd%time)
  if (id_snc>0) used = send_data(id_snc, snc(:)*lnd%ug_landfrac(:)*100, lnd%time)

  ! deallocate override buffer
  deallocate(phot_co2_data)

  call mpp_clock_end(landFastClock)
  call mpp_clock_end(landClock)
end subroutine update_land_model_fast


! ============================================================================
subroutine update_land_model_fast_0d ( tile, l,itile, N, land2cplr, &
   precip_l, precip_s, atmos_T, atmos_wind, &
   Ha0, DHaDTc, tr_flux, dfdtr, &
   ISa_dn_dir, ISa_dn_dif, ILa_dn, &
   ustar, p_surf, drag_q, &
   phot_co2_overridden, phot_co2_data, &
   runoff, runoff_c &
   )
  type (land_tile_type), pointer :: tile
  integer, intent(in) :: l ! position in unstructured grid
  integer, intent(in) :: itile ! tile number
  integer, intent(in) :: N ! number of cohorts in the tile
  type(land_data_type), intent(inout) :: land2cplr

  real, intent(in) :: &
       precip_l, precip_s, & ! liquid and solid precipitation, kg/(m2 s)
       atmos_T, &        ! incoming precipitation temperature (despite its name), deg K
       atmos_wind, &     ! magnitude of wind at the bottom of the atmosphere, m/s
       tr_flux(:), dfdtr(:), &  ! tracer flux from canopy air to the atmosphere
       Ha0,   DHaDTc, &  ! sensible heat flux from the canopy air to the atmosphere
       ISa_dn_dir(NBANDS), & ! downward direct sw radiation at the top of the canopy
       ISa_dn_dif(NBANDS), & ! downward diffuse sw radiation at the top of the canopy
       ILa_dn,             & ! downward lw radiation at the top of the canopy
       ustar,              & ! friction velocity, m/s
       p_surf,             & ! surface pressure, Pa
       drag_q,             & ! product of atmos_wind*CD_q, m/s
       phot_co2_data         ! data input for the CO2 for photosynthesis

  logical, intent(in):: phot_co2_overridden
  real, intent(inout) :: &
       runoff, &   ! total runoff of H2O, kg/m2
       runoff_c(:) ! runoff of tracers (including ice/snow and heat)

  ! ---- local vars
  real :: A(3*N+2,3*N+2),B0(3*N+2),B1(3*N+2),B2(3*N+2) ! implicit equation matrix and right-hand side vectors
  real :: ALUD(3*N+2,3*N+2) ! LU-decomposed A
  real :: X0(3*N+2),X1(3*N+2),X2(3*N+2) ! copy of the above, only for debugging
  integer :: indx(3*N+2) ! permutation vector
  ! indices of variables and equations for implicit time stepping solution :
  integer :: iqc, iTc, iTv, iwl, iwf
  ! linearization coefficients of various fluxes between components of land
  ! surface scheme
  real :: &
       Ea0,   DEaDqc, &  ! water vapor flux from canopy air to the atmosphere
       fco2_0,Dfco2Dq,&  ! co2 flux from canopy air to the atmosphere
       G0,    DGDTg,  &  ! ground heat flux
       Hg0,   DHgDTg,   DHgDTc, & ! linearization of the sensible heat flux from ground
       Eg0,   DEgDTg,   DEgDqc, DEgDpsig, & ! linearization of evaporation from ground
       flwg0, DflwgDTg, DflwgDTv(N), &  ! linearization of net LW radiation on the ground
       DqsatDTg   ! derivative of sat spec. humidity w.r.t. ground T, kg/kg/K
  real, dimension(N) :: & ! by cohort
       Hv0,   DHvDTv,   DHvDTc, & ! sens heat flux from vegetation
       Et0,   DEtDTv,   DEtDqc,   DEtDwl,   DEtDwf,  & ! transpiration
       Eli0,  DEliDTv,  DEliDqc,  DEliDwl,  DEliDwf, & ! evaporation of intercepted water
       Esi0,  DEsiDTv,  DEsiDqc,  DEsiDwl,  DEsiDwf, & ! sublimation of intercepted snow
       flwv0, DflwvDTg  ! linearization of net LW radiation to the canopy
  real :: &
       DflwvDTv(N,N) ! linearization of net LW radiation to the canopy:
         ! (i,j) == derivative of i-th cohort energy balance w.r.t. j-th cohort temperature

  real :: &
       vegn_prec_l(N), vegn_prec_s(N), & ! precipitation rate on top of each cohort, kg/(m2 s);
       ! the same as precip_{l,s} for the top layer, different for the cohorts below
       vegn_drip_l(N), vegn_drip_s(N), & ! drip rate of water and snow, respectively, kg/(m2 s)
       vegn_lai(N), vegn_lwnet(N)

  ! increments of respective variables over time step, results of the implicit
  ! time step:
  real :: delta_qc, delta_Tc, delta_Tv(N), delta_wl(N), delta_ws(N), delta_Tg, delta_psig, delta_co2
  real :: flwg ! updated value of long-wave ground energy balance
  real :: sum0, sum1, sum2

  real :: ndep_nit, ndep_amm, ndep_org  ! rates of nitrate, ammonium,
                ! and organic nitrogen input to the soil, kg N/(m2 yr)

  integer :: vegn_layer(N) ! layer number of each cohort
  real :: &
       f(N), & ! fraction of each cohort canopy in its layer
       vegn_T(N), vT(N), & ! vegetation (canopy) temperature
       cana_T, cT, & ! canopy air temperature
       evap_T, eT, & ! temperature assigned to vapor going between land and atmosphere
       soil_uptake_T(N), & ! average temperature of water taken up by the vegetation
       vegn_Wl(N),  vegn_Ws(N), & ! water and snow mass of the canopy
       vegn_ifrac(N), & ! intercepted fraction of liquid or frozen precipitation
       vegn_hcap(N),      & ! vegetation heat capacity, including intercepted water and snow
       vegn_fco2, & ! co2 flux from the vegetation, kg CO2/(m2 s)
       hlv_Tv(N), hlv_Tu(N), & ! latent heat of vaporization at vegn and uptake temperatures, respectively
       hls_Tv(N), &         ! latent heat of sublimation at vegn temperature
       con_v_v(N), con_st_v(N), & ! aerodynamic and stomatal conductance, respectively
       grnd_T, gT, & ! ground temperature and its value used for sensible heat advection
       grnd_q,         & ! specific humidity at ground surface
       grnd_rh,        & ! explicit relative humidity at ground surface
       grnd_rh_psi,    & ! psi derivative of relative humidity at ground surface
       grnd_liq, grnd_ice, grnd_subl, &
       grnd_tf, &  ! temperature of freezing on the ground
       grnd_qsat, & ! saturated water vapor on the ground
       grnd_latent, &
       grnd_flux, &
       grnd_E_min, &
       grnd_E_max, &
       soil_E_min, &
       soil_E_max, &
       swdn(N,NBANDS),  & ! downward short-wave radiation on top of the each cohort canopy, W/m2
       swnet(N,NBANDS), & ! net short-wave radiation balance of each cohort canopy, W/m2
       con_g_h, con_g_v, & ! turbulent cond. between ground and canopy air, for heat and vapor respectively
       snow_area, &
       cana_dens, & ! density of canopy air, kg/m3
       cana_q, & ! specific humidity of canopy air, kg/kg
       cana_qsat, & ! saturated specific humidity of canopy air, kg/kg
       cana_co2, & ! co2 moist mixing ratio in canopy air, kg CO2/kg wet air
       cana_co2_mol, & ! co2 dry mixing ratio in canopy air, mol CO2/mol dry air
       fswg, evapg, sensg, &
       subs_G, subs_G2, Mg_imp, snow_G_Z, snow_G_TZ, &
       snow_avrg_T, delta_T_snow,  & ! vertically-averaged snow temperature and its change
       vegn_ovfl_l,  vegn_ovfl_s,  & ! overflow of liquid and solid water from the canopy
       vegn_ovfl_Hl, vegn_ovfl_Hs, & ! heat flux from canopy due to overflow
       delta_fprec, & ! correction of below-canopy solid precip in case its average T > tfreeze

       hprec,              & ! sensible heat flux carried by precipitation
       hevap,              & ! sensible heat flux carried by total evapotranspiration
       land_evap,          & ! total vapor flux from land to atmosphere
       land_sens,          & ! turbulent sensible heat flux from land to atmosphere
       vegn_flw,vegn_sens(N),snow_sens,snow_levap,snow_fevap,snow_melt,&
       snow_lprec, snow_hlprec,snow_lrunf,vegn_levap(N),vegn_fevap(N),vegn_uptk(N),&
       vegn_fsw, vegn_melt, &
       vegn_lprec,  vegn_fprec,  & ! liquid and frozen precip under canopy, kg/(m2 s)
       vegn_hlprec, vegn_hfprec, & ! heat carried by liquid and frozen precip under canopy, J/(m2 s)
       prveg, & ! precip intercepted by vegetation, kg/(m2 s), for CMOR/CMIP output
       precip_T,pT,snow_fsw,snow_flw,snow_frunf,snow_hlrunf,&
       snow_hfrunf,subs_fsw,subs_flw,subs_sens,&
       subs_DT, subs_M_imp, subs_evap, snow_Tbot, snow_Cbot, snow_C, subs_levap,&
       subs_fevap,subs_melt,subs_lrunf,subs_hlrunf, subs_frunf, subs_hfrunf, &
       subs_tr_runf(n_river_tracers), & ! runoff of tracers from soil
       subs_Ttop, subs_Ctop, subs_subl, new_T

  real :: snow_T, snow_rh, snow_liq, snow_ice, snow_subl
  integer :: k, k1 ! cohort indices
  integer :: i, j, ii, jj ! indices for debug output
  integer :: ierr
  integer :: tr ! tracer index
  logical :: conserve_glacier_mass, snow_active, redo_leaf_water
  integer :: canopy_water_step, lw_step
  real :: subs_z0m, subs_z0s, snow_z0m, snow_z0s, grnd_z0s
  real :: lmass0, fmass0, heat0, cmass0, nmass0, nflux0, v0
  real :: lmass1, fmass1, heat1, cmass1, nmass1, nflux1
  real :: DOC_to_atmos
  logical :: calc_water_cons, calc_carbon_cons, calc_nitrogen_cons
  character(*), parameter :: tag = 'update_land_model_fast_0d'
  real :: lswept, fswept, hlswept, hfswept ! amounts of liquid and frozen snow, and corresponding
                                           ! heat swept with tiny snow

  calc_water_cons  = do_check_conservation.or.(id_water_cons>0)
  calc_carbon_cons = do_check_conservation.or.(id_carbon_cons>0)
  calc_nitrogen_cons = do_check_conservation.or.(id_nitrogen_cons>0)
  ! our coordinates in structured grid
  i = lnd%i_index(l); j = lnd%j_index(l)
  if(is_watch_point()) then
     write(*,*)
     call log_date('#### update_land_model_fast_0d begins:',lnd%time)
  endif
  ! send residuals to diag at the beginning of the time step so that the diagnostics captures
  ! adjustments made to residuals in other parts of the model, e.g. during cohort merge.
  ! NOTE that it is not entirely consistent with other energy balance components, which
  ! are sent to diag after the timestep, but it is still better than missing e_res_2 changes
  ! that come from outside of this subroutine.
  call send_tile_data(id_e_res_1, tile%e_res_1, tile%diag)
  call send_tile_data(id_e_res_2, tile%e_res_2, tile%diag)

  ! sanity checks of some input values
  call check_var_range(precip_l,  0.0, 1.0,        'land model input', 'precip_l',    WARNING)
  call check_var_range(precip_s,  0.0, 1.0,        'land model input', 'precip_s',    WARNING)
  call check_temp_range(atmos_T,                   'land model input', 'atmos_T')
  call check_var_range(ISa_dn_dir, 0.0, 1360.0,    'land model input', 'sw.down.dir', WARNING)
  call check_var_range(ISa_dn_dif, 0.0, 1360.0,    'land model input', 'sw.down.dif', WARNING)
  call check_var_range(ILa_dn,     0.0, 1360.0,    'land model input', 'lw.down',     WARNING)
  call check_var_range(ustar,      0.0, HUGE(1.0), 'land model input', 'ustar',       WARNING)
  call check_var_range(drag_q,     0.0, HUGE(1.0), 'land model input', 'drag_q',      WARNING)
  call check_var_range(p_surf,     0.0, HUGE(1.0), 'land model input', 'p_surf',      WARNING)
  ! not checking fluxes and their derivatives, since they can be either positive
  ! or negative, and it is hard to determine valid ranges for them.

  Ea0    = tr_flux(isphum) ; DEaDqc  = dfdtr(isphum)
  fco2_0 = tr_flux(ico2)   ; Dfco2Dq = dfdtr(ico2)

  ! + conservation check, part 1: calculate the pre-transition totals
  if (calc_water_cons)  call get_tile_water(tile,lmass0,fmass0)
  if (calc_carbon_cons) cmass0 = land_tile_carbon(tile)
  if (calc_nitrogen_cons) nmass0 = land_tile_nitrogen(tile)
  if(associated(tile%soil)) then
    nflux0=tile%soil%gross_nitrogen_flux_into_tile - tile%soil%gross_nitrogen_flux_out_of_tile
  else
    nflux0=0.0
  endif
  ! - end of conservation check, part 1

  ! if requested (in snow_nml), sweep tiny snow before calling step_1 subroutines to
  ! avoid numerical issues.
  call sweep_tiny_snow(tile%snow, lswept, fswept, hlswept, hfswept)

  soil_uptake_T(:) = tfreeze ! just to avoid using un-initialized values
  if (associated(tile%glac)) then
     call glac_step_1 ( tile%glac, &
          grnd_T, grnd_rh, grnd_liq, grnd_ice, grnd_subl, grnd_tf, &
          snow_G_Z, snow_G_TZ, conserve_glacier_mass  )
     grnd_E_min = -HUGE(grnd_E_min)
     grnd_E_max =  HUGE(grnd_E_max)
     grnd_rh_psi = 0
  else if (associated(tile%lake)) then
     call lake_step_1 ( ustar, p_surf, &
          lnd%ug_lat(l), tile%lake, &
          grnd_T, grnd_rh, grnd_liq, grnd_ice, grnd_subl, grnd_tf, &
          snow_G_Z, snow_G_TZ)
     grnd_E_min = -HUGE(grnd_E_min)
     grnd_E_max =  HUGE(grnd_E_max)
     grnd_rh_psi = 0
  else if (associated(tile%soil)) then
     call soil_step_1 ( tile%soil, tile%vegn, tile%diag, &
          grnd_T, soil_E_min, soil_E_max, &
          grnd_rh, grnd_rh_psi, grnd_liq, grnd_ice, grnd_subl, grnd_tf, &
          snow_G_Z, snow_G_TZ)
     grnd_E_min = soil_E_min
     grnd_E_max = soil_E_max
     grnd_liq = 0 ! sorry, but solver cannot handle implicit melt anymore
     grnd_ice = 0 ! sorry, but solver cannot handle implicit melt anymore
                  ! no big loss, it is just the surface layer anyway
  else
     call get_current_point(face=ii)
     call error_mesg('update_land_model_fast','none of the surface tiles exist at ('//&
          trim(string(i))//','//trim(string(j))//','//trim(string(itile))//&
          ', face='//trim(string(ii))//')',FATAL)
  endif

  ! + heat conservation check, part 1; land_tile_heat has to be called after
  !   soil_step_1, because soil dry heat capacity is initialized there
  ! heat0  = land_tile_heat(tile)

  subs_subl = grnd_subl

  call snow_step_1 ( tile%snow, snow_G_Z, snow_G_TZ, &
       snow_active, snow_T, snow_rh, snow_liq, snow_ice, &
       snow_subl, snow_area, G0, DGDTg )
  if (snow_active) then
     grnd_T    = snow_T;   grnd_rh   = snow_rh;   grnd_liq  = snow_liq
     grnd_rh_psi = 0
     grnd_ice  = snow_ice; grnd_subl = snow_subl; grnd_tf   = tfreeze
     grnd_E_min = -HUGE(grnd_E_min)
     grnd_E_max =  HUGE(grnd_E_max)
  endif

  cana_T   = tile%cana%T
  cana_q   = tile%cana%tr(isphum)
  cana_co2 = tile%cana%tr(ico2)

  if (associated(tile%vegn)) then
     ! calculate net short-wave radiation input to the vegetation
     do k = 1,N
        swnet(k,:) = tile%Sv_dir (k,:)*ISa_dn_dir + tile%Sv_dif (k,:)*ISa_dn_dif
        swdn (k,:) = tile%Sdn_dir(k,:)*ISa_dn_dir + tile%Sdn_dif(k,:)*ISa_dn_dif
     enddo
     ! calculate roughness of sub-canopy surface
     call soil_roughness(tile%soil, subs_z0s, subs_z0m)
     call snow_roughness(tile%snow, snow_z0s, snow_z0m)
     grnd_z0s = exp( (1-snow_area)*log(subs_z0s) + snow_area*log(snow_z0s))

     ! cana_co2 is moist mass mixing ratio [kg CO2/kg wet air], convert it to dry
     ! volumetric mixing ratio [mol CO2/mol dry air]
     cana_co2_mol = cana_co2*mol_air/mol_CO2/(1-cana_q)
     if (phot_co2_overridden) cana_co2_mol = phot_co2_data

     call vegn_step_1 ( tile%vegn, tile%soil, tile%diag, &
        p_surf, ustar, drag_q, &
        swdn, swnet, precip_l, precip_s, &
        tile%land_d, tile%land_z0s, tile%land_z0m, grnd_z0s, &
        cana_T, cana_q, cana_co2_mol, &
        ! output
        con_g_h, con_g_v, &
        con_v_v, con_st_v, & ! aerodyn. and stomatal conductances
        vegn_T, vegn_Wl, vegn_Ws, & ! temperature, water and snow mass on the canopy
        vegn_ifrac, vegn_lai, &
        vegn_drip_l, vegn_drip_s, &
        vegn_prec_l, vegn_prec_s, &
        vegn_lprec,  vegn_fprec,  &
        vegn_hcap, & ! total vegetation heat capacity (including intercepted water/snow)
        Hv0,   DHvDTv,   DHvDTc,            &
        Et0,   DEtDTv,   DEtDqc,   DEtDwl,   DEtDwf,  &
        Eli0,  DEliDTv,  DEliDqc,  DEliDwl,  DEliDwf, &
        Esi0,  DEsiDTv,  DEsiDqc,  DEsiDwl,  DEsiDwf, &
        soil_uptake_T )
     ! assign cohort layer area fractions (calculated in update_derived_vegn_properties)
     f(:) = tile%vegn%cohorts(1:N)%layerfrac
     vegn_layer(:) = tile%vegn%cohorts(1:N)%layer
     ! calculate precipitation intercepted by vegetation; need to be calculated here
     ! since vegn_lprec and vegn_fprec get modified with drip and overflow later
     prveg = precip_l + precip_s - vegn_lprec - vegn_fprec
  else ! i.e., no vegetation
     swnet    = 0
     con_g_h = con_fac_large ; con_g_v = con_fac_large
     con_v_v = 0.0; con_st_v = 0.0 ! does it matter?
     if(associated(tile%glac).and.conserve_glacier_mass.and..not.snow_active) &
          con_g_v = con_fac_small
     vegn_T  = cana_T ; vegn_Wl = 0 ; vegn_Ws = 0
     vegn_ifrac  = 0 ; vegn_lai    = 0
     vegn_drip_l = 0 ; vegn_drip_s = 0
     vegn_prec_l = 0 ; vegn_prec_s = 0
     vegn_hcap = 1.0
     Hv0 =0;  DHvDTv =0;  DHvDTc=0;
     Et0 =0;  DEtDTv =0;  DEtDqc=0;   DEtDwl=0;   DEtDwf=0
     Eli0=0;  DEliDTv=0;  DEliDqc=0;  DEliDwl=0;  DEliDwf=0
     Esi0=0;  DEsiDTv=0;  DEsiDqc=0;  DEsiDwl=0;  DEsiDwf=0
     f(:)=1;  vegn_layer(:) = 1
     prveg = 0.0
  endif

  ! calculate fluxes between canopy and ground surface, and their derivatives
  ! this used to be in cana_step_1
  call check_temp_range(grnd_T,'updata_land_model_fast_0d','grnd_T')
  call qscomp(grnd_T,p_surf,grnd_qsat,DqsatDTg)
  grnd_q    =  grnd_rh * grnd_qsat
  cana_dens =  p_surf/(rdgas*cana_T*(1+d608*cana_q))
  Hg0       =  cana_dens*cp_air*con_g_h*(grnd_T - cana_T)
  DHgDTg    =  cana_dens*cp_air*con_g_h
  DHgDTc    = -cana_dens*cp_air*con_g_h
  Eg0       =  cana_dens*con_g_v*(grnd_q - cana_q)
  DEgDTg    =  cana_dens*con_g_v*DqsatDTg*grnd_rh
  DEgDqc    = -cana_dens*con_g_v
  DEgDpsig  =  cana_dens*con_g_v*grnd_qsat*grnd_rh_psi

  ! calculate net shortwave for ground and canopy
  fswg     = SUM(tile%Sg_dir*ISa_dn_dir + tile%Sg_dif*ISa_dn_dif)
  vegn_fsw = 0
  do k = 1,N
     vegn_fsw = vegn_fsw+f(k)*SUM(swnet(k,:))
  enddo

! [X.0] calculate the latent heats of vaporization at appropriate temperatures
  if (use_tfreeze_in_grnd_latent) then
    grnd_latent = hlv + hlf*grnd_subl
  else
    grnd_latent = hlv + (cpw-clw)*(grnd_T-tfreeze) &
               + (hlf + (clw-csw)*(grnd_T-tfreeze)) * grnd_subl
  endif
  if (use_atmos_T_for_precip_T) then
    precip_T = atmos_T
  else
    precip_T = cana_T
  endif
  if (use_atmos_T_for_evap_T) then
    evap_T = atmos_T
  else
    evap_T = cana_T
  endif
  if (use_old_conservation_equations) then
    hlv_Tv = hlv       - (cpw-clw)*tfreeze + cpw*vegn_T
    hls_Tv = hlv + hlf - (cpw-csw)*tfreeze + cpw*vegn_T
    hlv_Tu = hlv       - (cpw-clw)*tfreeze + cpw*vegn_T - clw*soil_uptake_T
    pT = precip_T
    cT = cana_T
    eT = evap_T
    gT = grnd_T
    vT = vegn_T
  else
    hlv_Tv = hlv    + cpw*(vegn_T-tfreeze)
    hls_Tv = hlf    + hlv_Tv
    hlv_Tu = hlv_Tv - clw*(soil_uptake_T-tfreeze)
    pT = precip_T-tfreeze
    cT = cana_T-tfreeze
    eT = evap_T-tfreeze
    gT = grnd_T-tfreeze
    vT = vegn_T-tfreeze
  endif

! [X.X] using long-wave optical properties, calculate the explicit long-wave
!       radiative balances and their derivatives w.r.t. temperatures
  delta_Tv(:) = 0.0; delta_Tg = 0.0
  do lw_step = 1,2
     call land_lw_balance(ILa_dn, vegn_layer, f, vegn_T, delta_Tv, grnd_T, delta_Tg, &
        tile%vegn_tran_lw,tile%vegn_refl_lw,tile%surf_refl_lw, &
        flwv0, flwg0, DflwvDTv, DflwvDTg, DflwgDTv, DflwgDTg)

     do canopy_water_step = 1,max_canopy_water_steps
        if(is_watch_point()) then
           write(*,*)'#### input data for the matrix ####'
           __DEBUG1__(canopy_water_step)
           __DEBUG1__(delta_time)
           __DEBUG1__(canopy_air_mass)
           __DEBUG1__(vegn_T)
           __DEBUG1__(vT)
           __DEBUG1__(vegn_Wl)
           __DEBUG1__(vegn_Ws)
           __DEBUG2__(grnd_T,gT)
           __DEBUG1__(grnd_rh)
           __DEBUG2__(cana_T,cT)
           __DEBUG1__(cana_q)
           __DEBUG2__(evap_T,eT)
           __DEBUG2__(precip_T,pT)
           __DEBUG2__(precip_l, precip_s)
           __DEBUG1__(vegn_prec_l)
           __DEBUG1__(vegn_prec_s)
           __DEBUG1__(vegn_drip_l)
           __DEBUG1__(vegn_drip_s)
           __DEBUG1__(vegn_ifrac)
           __DEBUG1__(vegn_lai)
           __DEBUG1__(ILa_dn)
           __DEBUG2__(ISa_dn_dir(1),ISa_dn_dir(2))
           __DEBUG2__(ISa_dn_dif(1),ISa_dn_dif(2))
           __DEBUG1__(sum(swnet(:,:),2))
           __DEBUG2__(fswg, vegn_fsw)
           __DEBUG1__(vegn_hcap)
           __DEBUG1__(hlv_Tv)
           __DEBUG1__(hlv_Tu)
           __DEBUG1__(hls_Tv)
           __DEBUG2__(G0, DGDTg)
           __DEBUG2__(Ha0, DHaDTc)
           __DEBUG2__(Ea0, DEaDqc)
           __DEBUG1__(Hv0)
           __DEBUG1__(DHvDTv)
           __DEBUG1__(DHvDTc)
           __DEBUG1__(Et0)
           __DEBUG1__(DEtDTv)
           __DEBUG1__(DEtDqc)
           __DEBUG1__(DEtDwl)
           __DEBUG1__(DEtDwf)
           __DEBUG1__(Eli0)
           __DEBUG1__(DEliDTv)
           __DEBUG1__(DEliDqc)
           __DEBUG1__(DEliDwl)
           __DEBUG1__(DEliDwf)
           __DEBUG1__(Esi0)
           __DEBUG1__(DEsiDTv)
           __DEBUG1__(DEsiDqc)
           __DEBUG1__(DEsiDwl)
           __DEBUG1__(DEsiDwf)
           __DEBUG3__(Hg0, DHgDTg, DHgDTc)
           __DEBUG3__(Eg0, DEgDTg, DEgDqc)
           __DEBUG1__(flwv0)
           write(*,*)'DflwvDTv:'
           do k = 1,N
              write(*,'(i12.3,99(x,g23.16))') k, DflwvDTv(k,:)
           enddo
           __DEBUG2__(flwg0, DflwgDTg)
           __DEBUG1__(DflwgDTv)
      !     __DEBUG3__(flwv0(1), DflwvDTg(1), DflwvDTv(1,1))
      !     __DEBUG3__(flwg0, DflwgDTg, DflwgDTv(1))
           __DEBUG2__(tile%e_res_1,tile%e_res_2)
           __DEBUG1__(f)
        endif

        ! calculate indices for equation system:
        iqc=1;  iTc=2;  iTv=3;  iwl=iTv+N;  iwf=iwl+N

        A(:,:) = 0
      ! [X.1] form the system of equations for implicit scheme, such that A*X = B1*delta_Tg+B2*delta_psig+B0
      ! [X.1.1] equation of canopy air mass balance
        A(iqc,iqc) = canopy_air_mass/delta_time &
           -sum((DEtDqc(:)+DEliDqc(:)+DEsiDqc(:))*f(:))-DEgDqc+DEaDqc
        A(iqc,iTc) = 0
        do k = 1,N
           A(iqc,iTv+k-1) = -f(k)*(DEtDTv(k)+DEliDTv(k)+DEsiDTv(k))
           A(iqc,iwl+k-1) = -f(k)*(DEtDwl(k)+DEliDwl(k)+DEsiDwl(k))
           A(iqc,iwf+k-1) = -f(k)*(DEtDwf(k)+DEliDwf(k)+DEsiDwf(k))
        enddo
        B0(iqc)  = sum(f(:)*(Esi0(:)+Eli0(:)+Et0(:)))+Eg0-Ea0
        B1(iqc)  = DEgDTg
        B2(iqc)  = DEgDpsig
      ! [X.1.2] equation of canopy air energy balance
#ifdef USE_DRY_CANA_MASS
        A(iTc,iqc) = canopy_air_mass*cpw*cT/delta_time &
#else
        A(iTc,iqc) = canopy_air_mass*(cpw*cT-cp_air*cana_T)/delta_time &
#endif
             - cpw*sum(f(:)*vT(:)*(DEtDqc(:)+DEliDqc(:)+DEsiDqc(:))) &
             - cpw*gT*DEgDqc + cpw*eT*DEaDqc
#ifdef USE_DRY_CANA_MASS
        A(iTc,iTc) = canopy_air_mass*cp_air/delta_time &
#else
        A(iTc,iTc) = canopy_air_mass*(cp_air+cana_q*(cpw-cp_air))/delta_time &
#endif
             - sum(f(:)*DHvDTc(:)) - DHgDTc + DHaDTc
        do k = 1,N
           A(iTc,iTv+k-1) = -f(k)*(DHvDTv(k)+cpw*vT(k)*(DEtDTv(k)+DEliDTv(k)+DEsiDTv(k)))
           A(iTc,iwl+k-1) =            -f(k)*cpw*vT(k)*(DEtDwl(k)+DEliDwl(k)+DEsiDwl(k))
           A(iTc,iwf+k-1) =            -f(k)*cpw*vT(k)*(DEtDwf(k)+DEliDwf(k)+DEsiDwf(k))
        enddo
        B0(iTc)  = sum(f(:)*Hv0(:)) + Hg0 - Ha0 &
           + cpw*sum(f(:)*vT(:)*(Et0(:)+Eli0(:)+Esi0(:)))+cpw*gT*Eg0-cpw*eT*Ea0 &
           - tile%e_res_1 - tile%e_res_2
        B1(iTc)  = DHgDTg + cpw*gT*DEgDTg
        B2(iTc)  =          cpw*gT*DEgDpsig
      ! [X.1.3] equation of canopy energy balance
        do k = 1,N
           A(iTv+k-1,iqc) = hlv_Tu(k)*DEtDqc(k)+hlv_Tv(k)*DEliDqc(k)+hls_Tv(k)*DEsiDqc(k)
           A(iTv+k-1,iTc) = DHvDTc(k)
           A(iTv+k-1,iTv+k-1) = vegn_hcap(k)/delta_time &
             +DHvDTv(k) &
             +hlv_Tu(k)*DEtDTv(k) + hlv_Tv(k)*DEliDTv(k) + hls_Tv(k)*DEsiDTv(k) &
             +clw*vegn_drip_l(k) + csw*vegn_drip_s(k)
           ! add matrix of long-wave derivatives
           do k1=1,N
              A(iTv+k-1,iTv+k1-1) = A(iTv+k-1,iTv+k1-1)-DflwvDTv(k,k1)
           enddo
           A(iTv+k-1,iwl+k-1) = clw*vT(k)/delta_time &
             +hlv_Tu(k)*DEtDwl(k) + hlv_Tv(k)*DEliDwl(k) + hls_Tv(k)*DEsiDwl(k)
           A(iTv+k-1,iwf+k-1) = csw*vT(k)/delta_time &
             +hlv_Tu(k)*DEtDwf(k) + hlv_Tv(k)*DEliDwf(k) + hls_Tv(k)*DEsiDwf(k)
           B0(iTv+k-1) = sum(swnet(k,:)) &
             + flwv0(k) - Hv0(k) - hlv_Tu(k)*Et0(k) - Hlv_Tv(k)*Eli0(k) - hls_Tv(k)*Esi0(k) &
             + clw*vegn_prec_l(k)*vegn_ifrac(k)*pT + csw*vegn_prec_s(k)*vegn_ifrac(k)*pT & ! this is incorrect, needs to be modified. Is it?
             - clw*vegn_drip_l(k)*vT(k) - csw*vegn_drip_s(k)*vT(k)
           B1(iTv+k-1) = DflwvDTg(k)
           B2(iTv+k-1) = 0
        enddo
      ! [X.1.4] equation of intercepted liquid water mass balance
        do k = 1,N
           A(iwl+k-1,iqc) = DEliDqc(k)
           A(iwl+k-1,iTc) = 0
           A(iwl+k-1,iTv+k-1) = DEliDTv(k)
           A(iwl+k-1,iwl+k-1) = 1.0/delta_time + DEliDwl(k)
           A(iwl+k-1,iwf+k-1) = DEliDwf(k)
           B0(iwl+k-1)  = -Eli0(k) + vegn_prec_l(k)*vegn_ifrac(k) - vegn_drip_l(k)
           B1(iwl+k-1)  = 0
           B2(iwl+k-1)  = 0
        enddo
      ! [X.1.5] equation of intercepted frozen water mass balance
        do k = 1,N
           A(iwf+k-1,iqc) = DEsiDqc(k)
           A(iwf+k-1,iTc) = 0
           A(iwf+k-1,iTv+k-1) = DEsiDTv(k)
           A(iwf+k-1,iwl+k-1) = DEsiDwl(k)
           A(iwf+k-1,iwf+k-1) = 1.0/delta_time + DEsiDwf(k)
           B0(iwf+k-1)  = -Esi0(k) + vegn_prec_s(k)*vegn_ifrac(k) - vegn_drip_s(k)
           B1(iwf+k-1)  = 0
           B2(iwf+k-1)  = 0
        enddo
      ! [X.1.6] if LAI becomes zero (and, therefore, all fluxes from vegetation and
      ! their derivatives must be zero too) and heat capacity of the vegetation is
      ! zero, we get a degenerate case. Still, the drip may be non-zero because some
      ! water may remain from before leaf drop, and non-zero energy residual can be
      ! carried over from the previous time step.
      !
      ! To prevent temperature from going haywire in those cases, we simply replace the
      ! equations of canopy energy and mass balance with the following:
      ! vegn_T + delta_Tv = cana_T + delta_Tc
      ! delta_Wl = -vegn_drip_l*delta_time
      ! delta_Ws = -vegn_drip_s*delta_time
      ! The residual vegn_Wl and vegn_Ws, if any, are taken care of by the overflow
      ! calculations.
      !
      ! NOTE: currently vegn_hcap cannot be zero if mcv_min namelist parameter is not
      ! zero (and it is not by default, it is actually pretty big). Also, in non-vegetated
      ! tiles vegn_hcap is set to 1. So this degenerate case never happens in typical
      ! configurations. The only way for this to happen is to set mcv_min=0 and drop
      ! leaves
        do k = 1,N
           if(DHvDTv(k)==0.and.vegn_hcap(k)==0.and.vegn_lai(k)==0) then
             ! vegn_T + delta_Tv = cana_T + delta_Tc
             A(iTv+k-1,:)   = 0
             A(iTv+k-1,iTc) = -1
             A(iTv+k-1,iTv+k-1) = +1
             B0(iTv+k-1) = cana_T - vegn_T(k)
             B1(iTv+k-1) = 0
             ! delta_Wl = -vegn_drip_l*delta_time
             A(iwl+k-1,:)   = 0
             A(iwl+k-1,iwl+k-1) = 1
             B0(iwl+k-1) = -vegn_drip_l(k)*delta_time
             B1(iwl+k-1) = 0
             ! delta_Ws = -vegn_drip_s*delta_time
             A(iwf+k-1,:)   = 0
             A(iwf+k-1,iwf+k-1) = 1
             B0(iwf+k-1) = -vegn_drip_s(k)*delta_time
             B1(iwf+k-1) = 0
           endif
        enddo

        if(is_watch_point()) then
           write(*,*)'#### A ####'
           do ii = 1, size(A,1)
   !           write(*,'(99g23.16)')(A(ii,jj),jj=1,size(A,2))
               do jj = 1,size(A,2)
                  if (A(ii,jj) == 0) then
                     write(*,'(a3)',advance='no') '0'
                  else if (A(ii,jj) > 0) then
                     write(*,'(a3)',advance='no') '+'
                  else if (A(ii,jj) < 0) then
                     write(*,'(a3)',advance='no') '-'
                  endif
               enddo
               write(*,*)
           enddo
           write(*,*)'#### B0, B1, B2 ####'
           do ii = 1, size(A,1)
              write(*,'(i3.3, 99g23.16)')ii, B0(ii),B1(ii),B2(ii)
           enddo
        endif

      ! [X.2] solve the system for free terms and delta_Tg and delta_psig terms, getting
      !       linear equation for delta_Tg and delta_psig
        ALUD = A
        call ludcmp(ALUD,indx,ierr)
        if (ierr/=0)&
             call land_error_message('update_land_model_fast_0D: Matrix is singular',WARNING)
        if (improve_solution) then
           call lubksb_and_improve(A,ALUD,indx,B0,max_improv_steps,solution_tol,X0)
           call lubksb_and_improve(A,ALUD,indx,B1,max_improv_steps,solution_tol,X1)
           call lubksb_and_improve(A,ALUD,indx,B2,max_improv_steps,solution_tol,X2)
        else
           X0 = B0; X1=B1; X2=B2
           call lubksb(ALUD,indx,X0)
           call lubksb(ALUD,indx,X1)
           call lubksb(ALUD,indx,X2)
        endif

      !  if(is_watch_point()) then
      !     write(*,*)'#### solution: X0, X1, X2 ####'
      !     do ii = 1, size(A,1)
      !        __DEBUG3__(X0(ii),X1(ii),X2(ii))
      !     enddo
      !     write(*,*)'#### solution check ####'
      !     do ii = 1, size(A,1)
      !        sum0 = 0; sum1 = 0; sum2=0
      !        do jj = 1, size(A,2)
      !           sum0 = sum0 + A(ii,jj)*X0(jj)
      !           sum1 = sum1 + A(ii,jj)*X1(jj)
      !           sum2 = sum2 + A(ii,jj)*X2(jj)
      !        enddo
      !        __DEBUG3__(sum0-B0(ii),sum1-B1(ii),sum2-B2(ii))
      !     enddo
      !  endif
      ! the result of this solution is a set of expressions for delta_xx in terms
      ! of delta_Tg and delta_psig:
      ! delta_xx(i) = X0(i) + X1(i)*delta_Tg + X2(i)*delta_psig.

        ! solve the non-linear equation for energy balance at the surface.

        call land_surface_energy_balance( &
             grnd_T, grnd_liq, grnd_ice, grnd_latent, grnd_Tf, grnd_E_min, &
             grnd_E_max, fswg, &
             flwg0 + sum(X0(iTv:iTv+N-1)*DflwgDTv(:)), &
             DflwgDTg + sum(X1(iTv:iTv+N-1)*DflwgDTv(:)),&
             sum(X2(iTv:iTv+N-1)*DflwgDTv(:)), &
             Hg0 + X0(iTc)*DHgDTc, DHgDTg + X1(iTc)*DHgDTc, X2(iTc)*DHgDTc,   &
             Eg0 + X0(iqc)*DEgDqc, DEgDTg + X1(iqc)*DEgDqc, DEgDpsig + X2(iqc)*DEgDqc,   &
             G0,                       DGDTg, &
             ! output
             delta_Tg, delta_psig, Mg_imp )

      ! [X.5] calculate final value of other tendencies
        delta_qc = X0(iqc) + X1(iqc)*delta_Tg + X2(iqc)*delta_psig
        delta_Tc = X0(iTc) + X1(iTc)*delta_Tg + X2(iTc)*delta_psig
        delta_Tv(:) = X0(iTv:iTv+N-1) + X1(iTv:iTv+N-1)*delta_Tg + X2(iTv:iTv+N-1)*delta_psig
        delta_wl(:) = X0(iwl:iwl+N-1) + X1(iwl:iwl+N-1)*delta_Tg + X2(iwl:iwl+N-1)*delta_psig
        delta_ws(:) = X0(iwf:iwf+N-1) + X1(iwf:iwf+N-1)*delta_Tg + X2(iwf:iwf+N-1)*delta_psig

      ! [X.6] calculate updated values of energy balance components used in further
      !       calculations
        flwg       = flwg0 + DflwgDTg*delta_Tg + sum(DflwgDTv(:)*delta_Tv(:))
        evapg      = Eg0   + DEgDTg*delta_Tg   + DEgDpsig*delta_psig + DEgDqc*delta_qc
        sensg      = Hg0   + DHgDTg*delta_Tg   + DHgDTc*delta_Tc
        grnd_flux  = G0    + DGDTg*delta_Tg
        vegn_sens  = Hv0   + DHvDTv*delta_Tv   + DHvDTc*delta_Tc
        vegn_flw   = 0
        do k = 1,N
           vegn_levap(k) = Eli0(k)  + DEliDTv(k)*delta_Tv(k)  + DEliDqc(k)*delta_qc + DEliDwl(k)*delta_wl(k) + DEliDwf(k)*delta_ws(k)
           vegn_fevap(k) = Esi0(k)  + DEsiDTv(k)*delta_Tv(k)  + DEsiDqc(k)*delta_qc + DEsiDwl(k)*delta_wl(k) + DEsiDwf(k)*delta_ws(k)
           vegn_uptk (k) = Et0 (k)  + DEtDTv (k)*delta_Tv(k)  + DEtDqc (k)*delta_qc + DEtDwl (k)*delta_wl(k) + DEtDwf (k)*delta_ws(k)
           vegn_lwnet(k) = flwv0(k) + sum(DflwvDTv(k,:)*delta_Tv(:)) + DflwvDTg(k)*delta_Tg
           vegn_flw = vegn_flw + f(k)*vegn_lwnet(k)
        enddo
        land_evap  = Ea0   + DEaDqc*delta_qc
        land_sens  = Ha0   + DHaDTc*delta_Tc
        ! calculate the upward long-wave radiation flux from the land, to be returned to
        ! the flux exchange.
        tile%lwup = ILa_dn - vegn_flw - flwg

        if(is_watch_point())then
           write(*,*)'#### ground balance'
           __DEBUG2__(fswg,flwg)
           __DEBUG2__(sensg,evapg*grnd_latent)
           __DEBUG1__(grnd_flux)
           __DEBUG1__(Mg_imp)
           write(*,*)'#### implicit time steps'
           __DEBUG3__(delta_Tg, grnd_T,  grnd_T+delta_Tg )
           __DEBUG1__(delta_psig)
           __DEBUG2__(delta_qc, cana_q+delta_qc )
           __DEBUG2__(delta_Tc, cana_T+delta_Tc )
           __DEBUG1__(delta_Tv)
           __DEBUG1__(delta_wl)
           __DEBUG1__(delta_ws)
           __DEBUG1__(vegn_T+delta_Tv )
           __DEBUG1__(vegn_Wl+delta_wl)
           __DEBUG1__(vegn_Ws+delta_ws)
           write(*,*)'#### resulting fluxes'
           __DEBUG4__(flwg, evapg, sensg, grnd_flux)
           __DEBUG1__(vegn_levap)
           __DEBUG1__(vegn_fevap)
           __DEBUG1__(vegn_uptk)
           __DEBUG1__(vegn_sens)
           __DEBUG1__(vegn_lwnet)
           __DEBUG1__(vegn_flw)
           __DEBUG1__(land_evap)
        endif

        redo_leaf_water = .FALSE.
        if (prohibit_negative_canopy_water) then
           do k = 1,N
             if (vegn_Wl(k)+delta_wl(k)<0) then
                redo_leaf_water = .TRUE.
                Eli0(k) = vegn_Wl(k)/delta_time + vegn_prec_l(k)*vegn_ifrac(k) - vegn_drip_l(k)
                DEliDTv(k) = 0.0;  DEliDqc(k) = 0.0
                DEliDwl(k) = 0.0;  DEliDwf(k) = 0.0
             endif
             if (vegn_Ws(k)+delta_ws(k)<0) then
                redo_leaf_water = .TRUE.
                Esi0(k) = vegn_Ws(k)/delta_time + vegn_prec_s(k)*vegn_ifrac(k) - vegn_drip_s(k)
                DEsiDTv(k) = 0.0;  DEsiDqc(k) = 0.0
                DEsiDwl(k) = 0.0;  DEsiDwf(k) = 0.0
             endif
           enddo
        endif
        if (.not.redo_leaf_water) exit ! from loop
     enddo ! canopy_water_step
     ! call check_var_range(delta_Tv,  -HUGE(1.0), lw_delta_T_thresh, 'for LW derivative improvement', 'delta_Tv', WARNING)
     ! call check_var_range(delta_Tg,  -HUGE(1.0), lw_delta_T_thresh, 'for LW derivative improvement', 'delta_Tg', WARNING)
     if (all(delta_Tv<lw_delta_T_thresh).and.delta_Tg<lw_delta_T_thresh) exit ! from loop
     ! otherwise redo longwave radiation calculations with new values of delta_Tv, delta_Tg
     ! to get better approximation of long-wave radiation derivatives.
  enddo ! lw_step

  ! calculate energy residuals due to cross-product of time tendencies
#ifdef USE_DRY_CANA_MASS
  tile%e_res_1 = canopy_air_mass*cpw*delta_qc*delta_Tc/delta_time
#else
  tile%e_res_1 = canopy_air_mass*(cpw-cp_air)*delta_qc*delta_Tc/delta_time
#endif
  tile%e_res_2 = sum(f(:)*delta_Tv(:)*(clw*delta_Wl(:)+csw*delta_Ws(:)))/delta_time
  if (is_watch_point()) then
     __DEBUG2__(tile%e_res_1, tile%e_res_2)
  endif

  ! [*] start of step_2 updates
  ! update canopy air temperature and specific humidity
  tile%cana%T = tile%cana%T + delta_Tc
  tile%cana%tr(isphum) = tile%cana%tr(isphum) + delta_qc

  if(associated(tile%vegn)) then
     call vegn_step_2 ( tile%vegn, tile%diag, &
          delta_Tv, delta_wl, delta_ws, &
          vegn_melt,  &
          vegn_ovfl_l,   vegn_ovfl_s, &
          vegn_ovfl_Hl, vegn_ovfl_Hs )
     ! calculate heat carried by liquid and solid precipitation below the canopy
     vegn_hlprec = clw*(vegn_lprec*(precip_T-tfreeze) &
                      + sum(f(:)*vegn_drip_l(:)*(vegn_T(:)+delta_Tv(:)-tfreeze)) &
                      ) + vegn_ovfl_Hl
     vegn_hfprec = csw*(vegn_fprec*(precip_T-tfreeze) &
                      + sum(f(:)*vegn_drip_s(:)*(vegn_T(:)+delta_Tv(:)-tfreeze)) &
                      ) + vegn_ovfl_Hs
     ! calculate total amount of liquid and solid precipitation below the canopy
     vegn_lprec  = vegn_lprec + sum(f(:)*vegn_drip_l(:)) + vegn_ovfl_l
     vegn_fprec  = vegn_fprec + sum(f(:)*vegn_drip_s(:)) + vegn_ovfl_s
     ! make sure the temperature of the snow falling below canopy is below freezing
     ! this correction was introduced in an attempt to fix the problem with fictitious
     ! heat accumulating in near-zero-mass snow; however it does not seem to make a
     ! difference.
     if(vegn_hfprec>0)then
        ! solid precipitation from vegetation carries positive energy -- we cannot have
        ! that, because that would bring snow T above tfreeze, so convert excess to
        ! liquid
        delta_fprec = min(vegn_fprec,vegn_hfprec/hlf)
        vegn_fprec = vegn_fprec - delta_fprec
        vegn_lprec = vegn_lprec + delta_fprec
        vegn_hfprec = vegn_hfprec - hlf*delta_fprec
        ! we do not need to correct the vegn_hlprec since the temperature of additional
        ! liquid precip is tfreeze, and therefore its contribution to vegn_hlprec is
        ! exactly zero
     endif
     ! possibly we need to correct for the opposite situation: negative energy carried
     ! by liquid precipitation.
  else
     vegn_lprec  = precip_l
     vegn_fprec  = precip_s
     vegn_hlprec = precip_l*clw*(precip_T-tfreeze)
     vegn_hfprec = precip_s*csw*(precip_T-tfreeze)
     ! the fields below are only used in diagnostics
     vegn_melt   = 0
     vegn_fsw    = 0
  endif

  call snow_step_2 ( tile%snow, &
       snow_subl, vegn_lprec, vegn_fprec, vegn_hlprec, vegn_hfprec, &
       delta_Tg, Mg_imp, evapg, fswg, flwg, sensg, &
       use_tfreeze_in_grnd_latent, &
       ! output:
       subs_DT, subs_M_imp, subs_evap, subs_fsw, subs_flw, subs_sens, &
       snow_fsw, snow_flw, snow_sens, &
       snow_levap, snow_fevap, snow_melt, &
       snow_lprec, snow_hlprec, snow_lrunf, snow_frunf, &
       snow_hlrunf, snow_hfrunf, snow_Tbot, snow_Cbot, snow_C, snow_avrg_T )
  snow_lrunf  = snow_lrunf  + lswept/delta_time
  snow_frunf  = snow_frunf  + fswept/delta_time
  snow_hlrunf = snow_hlrunf + hlswept/delta_time
  snow_hfrunf = snow_hfrunf + hfswept/delta_time
  if(is_watch_point()) then
     write(*,*) 'subs_M_imp', subs_M_imp
  endif

  if (snow_active) then
     subs_G = snow_G_Z+snow_G_TZ*subs_DT
  else
     subs_G = 0
  endif

  if (associated(tile%glac)) then
     call glac_step_2 &
          ( tile%glac, tile%diag, subs_subl, snow_lprec, snow_hlprec, &
          subs_DT, subs_M_imp, subs_evap, &
          subs_levap, subs_fevap, &
          subs_melt, subs_lrunf, subs_hlrunf, subs_Ttop, subs_Ctop )
     subs_frunf = 0.
     subs_hfrunf = 0.
     subs_tr_runf(:) = 0.
     DOC_to_atmos = 0.
  else if (associated(tile%lake)) then
     call lake_step_2 &
          ( tile%lake, tile%diag, subs_subl, snow_lprec, snow_hlprec, &
          subs_DT, subs_M_imp, subs_evap, &
          use_tfreeze_in_grnd_latent, subs_levap, subs_fevap, &
          subs_melt, subs_Ttop, subs_Ctop )
     subs_lrunf = 0.
     subs_hlrunf = 0.
     subs_frunf = 0.
     subs_hfrunf = 0.
     subs_tr_runf = 0.
     DOC_to_atmos = 0.
  else if (associated(tile%soil)) then
     call soil_step_2 &
          ( tile%soil, tile%vegn, tile%diag, subs_subl, snow_lprec, snow_hlprec, &
          vegn_uptk, subs_DT, subs_M_imp, subs_evap, &
          use_tfreeze_in_grnd_latent, &
          ! output:
          subs_levap, subs_fevap, &
          subs_melt, subs_lrunf, subs_hlrunf, subs_Ttop, subs_Ctop, &
          subs_frunf, subs_hfrunf, subs_tr_runf, DOC_to_atmos)
  endif
  if (is_watch_point()) then
     __DEBUG2__(subs_levap, subs_fevap)
     __DEBUG5__(subs_melt, subs_lrunf, subs_hlrunf, subs_Ttop, subs_Ctop)
  endif

! TEMP FIX: MAIN PROG SHOULD NOT TOUCH CONTENTS OF PROG VARS. ******
! ALSO, DIAGNOSTICS IN COMPONENT MODULES SHOULD _FOLLOW_ THIS ADJUSTMENT******
  if (LM2) then
     tile%snow%T = subs_Ttop
     subs_G2 = 0.
  else
     if (sum(tile%snow%ws(:))>0)then
        new_T = (subs_Ctop*subs_Ttop +snow_Cbot*snow_Tbot) &
                        / (subs_Ctop+snow_Cbot)
        tile%snow%T(size(tile%snow%T)) = new_T
        if(associated(tile%glac)) tile%glac%T(1) = new_T
        if(associated(tile%lake)) tile%lake%T(1) = new_T
        if(associated(tile%soil)) tile%soil%T(1) = new_T
        subs_G2 = subs_Ctop*(new_T-subs_Ttop)/delta_time
     else
        if(tau_snow_T_adj>=0) then
           delta_T_snow = subs_Ctop*(subs_Ttop-snow_avrg_T)/&
                (subs_Ctop*tau_snow_T_adj/delta_time+subs_Ctop+snow_C)
           tile%snow%T(:) = snow_avrg_T + delta_T_snow

           new_T = subs_Ttop-snow_C/subs_Ctop*delta_T_snow
           if(associated(tile%glac)) tile%glac%T(1) = new_T
           if(associated(tile%lake)) tile%lake%T(1) = new_T
           if(associated(tile%soil)) tile%soil%T(1) = new_T
           subs_G2 = subs_Ctop*(new_T-subs_Ttop)/delta_time
        else
           subs_G2 = 0.
        endif
     endif
  endif

  vegn_fco2 = 0
  if (associated(tile%vegn)) then
     ! do the calculations that require updated land surface prognostic variables
     call nitrogen_sources(lnd%time, l, tile%vegn%p_ann, precip_l+precip_s, &
             tile%vegn%landuse, ndep_nit, ndep_amm, ndep_org, tile%diag)
     call vegn_step_3 (tile%vegn, tile%soil, tile%cana%T, precip_l+precip_s, &
          ndep_nit, ndep_amm, ndep_org, vegn_fco2, tile%diag)
     ! if vegn is present, then soil must be too
     call soil_step_3(tile%soil, tile%diag)

     call update_fire_fast(tile, p_surf, atmos_wind, l)
  endif

  ! update co2 concentration in the canopy air. It would be more consistent to do that
  ! in the same place and fashion as the rest of prognostic variables: that is, have the
  ! vegn_step_1 (and perhaps other *_step_1 procedures) calculate fluxes and their
  ! derivatives, then solve the linear equation(s), and finally update
  ! the concentration.
  if(update_cana_co2) then
     delta_co2 = (vegn_fco2 + DOC_to_atmos*mol_CO2/mol_C - fco2_0)/(canopy_air_mass_for_tracers/delta_time+Dfco2Dq)
     tile%cana%tr(ico2) = tile%cana%tr(ico2) + delta_co2
  else
     delta_co2 = 0
  endif
  if(is_watch_point())then
     __DEBUG1__(tile%cana%tr(ico2))
     __DEBUG4__(fco2_0,Dfco2Dq,vegn_fco2, DOC_to_atmos)
  endif

  call update_cana_tracers(tile, l, tr_flux, dfdtr, &
           precip_l, precip_s, p_surf, ustar, con_g_v, con_v_v, con_st_v )

  call update_land_bc_fast (tile, N, l, itile, land2cplr)

  ! accumulate runoff variables over the tiles
  runoff = runoff + (snow_frunf + subs_lrunf + snow_lrunf + subs_frunf)*tile%frac
  do tr = 1,n_river_tracers
     if (tr==i_river_heat) then
        runoff_c(tr) = runoff_c(tr) + (snow_hfrunf + subs_hlrunf + snow_hlrunf + subs_hfrunf)*tile%frac
     else if (tr==i_river_ice) then
        runoff_c(tr) = runoff_c(tr) + (snow_frunf + subs_frunf)*tile%frac
     else
        runoff_c(tr) = runoff_c(tr) + subs_tr_runf(tr) * tile%frac
     endif
  enddo
  hprec = (clw*precip_l+csw*precip_s)*(precip_T-tfreeze)
  hevap = cpw*land_evap*(evap_T-tfreeze)

  if (is_watch_cell()) then
     write(*,*)'Accumulated runoff for watch_cell'
     __DEBUG3__(itile, runoff, runoff_c)
  end if

  ! + conservation check, part 2: calculate totals in final state, and compare
  ! with previous totals
  if (calc_water_cons) then
     call get_tile_water(tile,lmass1,fmass1)
     if (do_check_conservation) call check_conservation (tag,'water', &
         lmass0+fmass0+(precip_l+precip_s-land_evap-(snow_frunf+subs_lrunf+snow_lrunf))*delta_time, &
         lmass1+fmass1, water_cons_tol)
     v0=lmass0+fmass0+(precip_l+precip_s-land_evap-(snow_frunf+subs_lrunf+snow_lrunf))*delta_time
     call send_tile_data(id_water_cons, (lmass1+fmass1-v0)/delta_time, tile%diag)
  endif
  if(calc_carbon_cons) then
     ! BNS: This will likely fail if DOC river tracer is not present and running in CORPSE mode with C_leaching_solubility>0
     v0 = cmass0-(fco2_0+Dfco2Dq*delta_co2)*mol_C/mol_CO2*delta_time
     if (i_river_DOC/=NO_TRACER) &
         v0 = v0 - subs_tr_runf(i_river_DOC)*delta_time
     cmass1 = land_tile_carbon(tile)
     if (do_check_conservation) &
           call check_conservation (tag,'carbon', v0, cmass1, carbon_cons_tol)
     call send_tile_data(id_carbon_cons, (cmass1-v0)/delta_time, tile%diag)
  endif
  if(calc_nitrogen_cons) then
     nmass1 = land_tile_nitrogen(tile)
     if(associated(tile%soil)) then
       nflux1=tile%soil%gross_nitrogen_flux_into_tile - tile%soil%gross_nitrogen_flux_out_of_tile
     else
       nflux1=0.0
     endif

     if (do_check_conservation.and.(soil_carbon_option==SOILC_CORPSE_N)) &
        call check_conservation (tag,'nitrogen', nmass0, nmass1 + (nflux0 - nflux1), nitrogen_cons_tol)
     call send_tile_data(id_nitrogen_cons, (nmass1-nflux1-nmass0+nflux0)/delta_time, tile%diag)
  endif
  ! heat1  = land_tile_heat(tile)
  ! latent heat is missing below, and it is not trivial to add, because there are
  ! multiple components with their own vaporization heat
  !  call check_conservation (tag,'heat content', &
  !      heat0+(hprec-land_sens-hevap &
  !           +sum(ISa_dn_dir*(1-tile%land_refl_dir)+ISa_dn_dif*(1-tile%land_refl_dif)) &
  !           +ILa_dn-tile%lwup &
  !           -(snow_hfrunf + subs_hlrunf + snow_hlrunf) &
  !           )*delta_time, &
  !      heat1, 1e-16, lnd%time)
  ! - end of conservation check, part 2

  ! TODO: go through the diagnostics and verify that they do the right thing in PPA case
  ! ---- diagnostic section ----------------------------------------------
  call send_tile_data(id_total_C, cmass1,                             tile%diag)
  call send_tile_data(id_total_N, nmass1,                             tile%diag)
  call send_tile_data(id_frac,    tile%frac,                          tile%diag)
  call send_tile_data(id_ntiles,  1.0,                                tile%diag)
  call send_tile_data(id_precip,  precip_l+precip_s,                  tile%diag)
  call send_tile_data(id_hprec,   hprec,                              tile%diag)
  call send_tile_data(id_lprec,   precip_l,                           tile%diag)
  call send_tile_data(id_lprecv,  precip_l-vegn_lprec,                tile%diag)
  call send_tile_data(id_lprecs,  vegn_lprec-snow_lprec,              tile%diag)
  call send_tile_data(id_lprecg,  snow_lprec,                         tile%diag)
  call send_tile_data(id_hlprec,  clw*precip_l*(precip_T-tfreeze),    tile%diag)
  call send_tile_data(id_hlprecv, clw*precip_l*(precip_T-tfreeze)-vegn_hlprec, &
                                                                      tile%diag)
  call send_tile_data(id_hlprecs, vegn_hlprec-snow_hlprec,            tile%diag)
  call send_tile_data(id_hlprecg, snow_hlprec,                        tile%diag)
  call send_tile_data(id_fprec,   precip_s,                           tile%diag)
  call send_tile_data(id_fprecv,  precip_s-vegn_fprec,                tile%diag)
  call send_tile_data(id_fprecs,  vegn_fprec,                         tile%diag)
  call send_tile_data(id_hfprec,  csw*precip_s*(precip_T-tfreeze),    tile%diag)
  call send_tile_data(id_hfprecv, csw*precip_s*(precip_T-tfreeze)-vegn_hfprec, &
                                                                      tile%diag)
  call send_tile_data(id_hfprecs, vegn_hfprec,                        tile%diag)
  call send_tile_data(id_evap,    land_evap,                          tile%diag)
  call send_tile_data(id_hevap,   hevap,                              tile%diag)
  call send_tile_data(id_levap,   sum(f(:)*(vegn_levap+vegn_uptk))+snow_levap+subs_levap, &
                                                                      tile%diag)
  call send_tile_data(id_levapv,  sum(f(:)*vegn_levap),               tile%diag)
  call send_tile_data(id_levaps,  snow_levap,                         tile%diag)
  call send_tile_data(id_levapg,  subs_levap,                         tile%diag)
  call send_tile_data(id_hlevap,  sum(f(:)*cpw*vegn_levap*(vegn_T-tfreeze)) &
                                    +cpw*snow_levap*(snow_T-tfreeze) &
                                    +cpw*subs_levap*(grnd_T-tfreeze), tile%diag)
  call send_tile_data(id_hlevapv, sum(f(:)*cpw*vegn_levap*(vegn_T-tfreeze)), tile%diag)
  call send_tile_data(id_hlevaps, cpw*snow_levap*(snow_T-tfreeze),    tile%diag)
  call send_tile_data(id_hlevapg, cpw*subs_levap*(grnd_T-tfreeze),    tile%diag)
  call send_tile_data(id_fevap,   sum(f(:)*vegn_fevap)+snow_fevap+subs_fevap, tile%diag)
  call send_tile_data(id_fevapv,  sum(f(:)*vegn_fevap),                 tile%diag)
  call send_tile_data(id_fevaps,  snow_fevap,                         tile%diag)
  call send_tile_data(id_fevapg,  subs_fevap,                         tile%diag)
  call send_tile_data(id_hfevap,  sum(f(:)*cpw*vegn_fevap*(vegn_T-tfreeze)) &
                                    +cpw*snow_fevap*(snow_T-tfreeze) &
                                    +cpw*subs_fevap*(grnd_T-tfreeze), tile%diag)
  call send_tile_data(id_hfevapv, sum(f(:)*cpw*vegn_fevap*(vegn_T-tfreeze)), tile%diag)
  call send_tile_data(id_hfevaps, cpw*snow_fevap*(snow_T-tfreeze),    tile%diag)
  call send_tile_data(id_hfevapg, cpw*subs_fevap*(grnd_T-tfreeze),    tile%diag)
  call send_tile_data(id_runf,    snow_lrunf+snow_frunf+subs_lrunf,   tile%diag)
  call send_tile_data(id_hrunf,   snow_hlrunf+snow_hfrunf+subs_hlrunf,tile%diag)
  call send_tile_data(id_lrunf,   snow_lrunf+subs_lrunf,              tile%diag)
  call send_tile_data(id_lrunfs,  snow_lrunf,                         tile%diag)
  call send_tile_data(id_lrunfg,  subs_lrunf,                         tile%diag)
  call send_tile_data(id_hlrunf,  snow_hlrunf+subs_hlrunf,            tile%diag)
  call send_tile_data(id_hlrunfs, snow_hlrunf,                        tile%diag)
  call send_tile_data(id_hlrunfg, subs_hlrunf,                        tile%diag)
  call send_tile_data(id_frunf,   snow_frunf + subs_frunf,            tile%diag)
  call send_tile_data(id_frunfs,  snow_frunf,                         tile%diag)
  call send_tile_data(id_hfrunf,  snow_hfrunf + subs_hfrunf,          tile%diag)
  call send_tile_data(id_hfrunfs, snow_hfrunf,                        tile%diag)
  ! TODO: generalize diagnostic for runoff of tracers
  ! BNS: This fails because id_runf_tr is allocated and set in land_diag_init before n_river_tracers is set in land_model_init
  ! do tr = 3, n_river_tracers
  !    if (id_runf_tr(tr)>0) &
  !       call send_tile_data(id_runf_tr(tr), subs_tr_runf(tr),         tile%diag)
  ! enddo
  call send_tile_data(id_melt,    vegn_melt+snow_melt+subs_melt,      tile%diag)
  call send_tile_data(id_meltv,   vegn_melt,                          tile%diag)
  call send_tile_data(id_melts,   snow_melt,                          tile%diag)
  call send_tile_data(id_meltg,   subs_melt,                          tile%diag)
  call send_tile_data(id_fsw,     vegn_fsw+snow_fsw+subs_fsw,         tile%diag)
  call send_tile_data(id_fswv,    vegn_fsw,                           tile%diag)
  call send_tile_data(id_fsws,    snow_fsw,                           tile%diag)
  call send_tile_data(id_fswg,    subs_fsw,                           tile%diag)
  call send_tile_data(id_flw,     vegn_flw+snow_flw+subs_flw,         tile%diag)
  call send_tile_data(id_flwv,    vegn_flw,                           tile%diag)
  call send_tile_data(id_flws,    snow_flw,                           tile%diag)
  call send_tile_data(id_flwg,    subs_flw,                           tile%diag)
  call send_tile_data(id_sens,    land_sens,                          tile%diag)
  call send_tile_data(id_sensv,   vegn_sens,                          tile%diag)
  call send_tile_data(id_senss,   snow_sens,                          tile%diag)
  call send_tile_data(id_sensg,   subs_sens,                          tile%diag)
  call send_tile_data(id_con_g_h, con_g_h,                            tile%diag)
  call send_tile_data(id_con_g_v, con_g_v,                            tile%diag)
  call send_tile_data(id_transp,  sum(f(:)*vegn_uptk),                tile%diag)
  call send_tile_data(id_wroff,   snow_lrunf+subs_lrunf,              tile%diag)
  call send_tile_data(id_sroff,   snow_frunf,                         tile%diag)
  call send_tile_data(id_htransp, sum(f(:)*cpw*vegn_uptk*(vegn_T-tfreeze)), tile%diag)
  call send_tile_data(id_huptake, sum(f(:)*clw*vegn_uptk*(soil_uptake_T-tfreeze)), &
                                                                      tile%diag)
  call send_tile_data(id_hroff,   snow_hlrunf+subs_hlrunf+snow_hfrunf, &
                                                                      tile%diag)
  call send_tile_data(id_gsnow,   subs_G,                             tile%diag)
  call send_tile_data(id_gequil,  subs_G2,                            tile%diag)
  call send_tile_data(id_grnd_flux, grnd_flux,                        tile%diag)
  ! an arbitrary constant 1000.0 is used to avoid numbers that cannot be represented
  ! in floating point history output
  call send_tile_data(id_levapg_max, min(grnd_E_max,1000.0),          tile%diag)
  call send_tile_data(id_qco2,    tile%cana%tr(ico2),                 tile%diag)
  call send_tile_data(id_qco2_dvmr,&
       tile%cana%tr(ico2)*mol_air/mol_co2/(1-tile%cana%tr(isphum)),   tile%diag)
  call send_tile_data(id_fco2,  vegn_fco2*mol_C/mol_CO2 + DOC_to_atmos, tile%diag)
  call send_tile_data(id_co2_mol_flux, fco2_0 + Dfco2Dq*delta_co2,      tile%diag)
  call send_tile_data(id_swdn_dir, ISa_dn_dir,                          tile%diag)
  call send_tile_data(id_swdn_dif, ISa_dn_dif,                          tile%diag)
  call send_tile_data(id_swup_dir, ISa_dn_dir*tile%land_refl_dir,       tile%diag)
  call send_tile_data(id_swup_dif, ISa_dn_dif*tile%land_refl_dif,       tile%diag)
  call send_tile_data(id_lwdn,     ILa_dn,                            tile%diag)
  call send_tile_data(id_subs_emis,1-tile%surf_refl_lw,               tile%diag)

  call send_tile_data(id_grnd_rh, grnd_rh, tile%diag)
  if (id_cana_rh > 0) then
     call qscomp(tile%cana%T, p_surf, cana_qsat)
     call send_tile_data(id_cana_rh, tile%cana%tr(isphum)/cana_qsat, tile%diag)
  endif

  ! CMOR/CMIP variables
  call send_tile_data(id_pcp,    precip_l+precip_s,                   tile%diag)
  call send_tile_data(id_prra,   precip_l,                            tile%diag)
  ! need to send prveg from here (rather than from vegn_step_1) to cover places where
  ! vegetation does not exist
  call send_tile_data(id_prveg,  prveg, tile%diag)
  ! note that the expression below gives only approximate value of latent heat
  ! flux, since in the model specific heat of vaporization is not equal to hlv;
  ! it depends on temperature and phase state.
  call send_tile_data(id_hflsLut, land_evap*hlv,                      tile%diag)
  call send_tile_data(id_rlusLut, tile%lwup,                          tile%diag)
  if (id_rsusLut>0) call send_tile_data(id_rsusLut, &
    sum(ISa_dn_dir*tile%land_refl_dir+ISa_dn_dif*tile%land_refl_dif), tile%diag)
  ! evspsblsoi is evaporation from *soil*, so we send zero from glaciers and lakes;
  ! the result is averaged over the entire land surface, as required by CMIP. evspsblveg
  ! does not need this distinction because it is already zero over glaciers and lakes.
  if (associated(tile%soil)) then
     call send_tile_data(id_evspsblsoi, subs_levap+subs_fevap,        tile%diag)
  else
     call send_tile_data(id_evspsblsoi, 0.0,                          tile%diag)
  endif
  if(id_evspsblveg > 0) call send_tile_data(id_evspsblveg, sum(f(:)*(vegn_levap+vegn_fevap)), tile%diag)
  if(id_ec         > 0) call send_tile_data(id_ec,         sum(f(:)*(vegn_levap+vegn_fevap)), tile%diag)
  if(associated(tile%lake)) then
     call send_tile_data(id_eow, subs_levap,                          tile%diag)
  else
     call send_tile_data(id_eow, 0.0,                                 tile%diag)
  endif
  call send_tile_data(id_esn, snow_fevap+snow_levap,                  tile%diag)
  call send_tile_data(id_et,  land_evap,                              tile%diag)
  call send_tile_data(id_snm, snow_melt,                              tile%diag)
  if (id_hfdsn>0) then
     if (snow_active) then
        call send_tile_data(id_hfdsn, &
            snow_flw + snow_fsw & ! net radiation
            - snow_sens & ! turbilent sensible with canopy air
            + vegn_hfprec + vegn_hlprec & ! sensible heat coming with precipitation
            - cpw*(snow_fevap+snow_levap)*(snow_T-tfreeze) & ! sensible heat carried away by water vapor
            - (snow_fevap+snow_levap)*grnd_latent, & ! latent heat
            tile%diag)
     else
        call send_tile_data(id_hfdsn, 0.0, tile%diag)
     endif
  endif
  if (id_tslsiLut>0) &
      call send_tile_data(id_tslsiLut, (tile%lwup/stefan)**0.25,      tile%diag)
  if (id_cLand > 0) &
      call send_tile_data(id_cLand, land_tile_carbon(tile),           tile%diag)
  if (id_cTot1 > 0) &
      call send_tile_data(id_cTot1, land_tile_carbon(tile)-cana_tile_carbon(tile%cana), tile%diag)
  if (id_nLand > 0) &
      call send_tile_data(id_nLand, land_tile_nitrogen(tile),         tile%diag)
  if (id_nbp>0) call send_tile_data(id_nbp, -vegn_fco2*mol_C/mol_CO2-DOC_to_atmos, tile%diag)
end subroutine update_land_model_fast_0d


! ============================================================================
subroutine update_land_model_slow ( cplr2land, land2cplr )
  type(atmos_land_boundary_type), intent(inout) :: cplr2land
  type(land_data_type)          , intent(inout) :: land2cplr

  ! ---- local vars
  integer :: l,k
  integer :: second, minute, hour, day0, day1, month0, month1, year0, year1
  integer :: n_cohorts
  type(land_tile_type), pointer :: tile
  type(land_tile_enum_type) :: ce
  type(land_tile_list_type) :: tmp

  call mpp_clock_begin(landClock)
  call mpp_clock_begin(landSlowClock)

  call send_cellfrac_data(id_cropFrac, is_crop)
  call send_cellfrac_cohort_data(id_cropFracC3, is_crop, is_c3)
  call send_cellfrac_cohort_data(id_cropFracC4, is_crop, is_c4)
  call send_cellfrac_data(id_pastureFrac,  is_pasture)
  call send_cellfrac_data(id_residualFrac, is_residual)

  call send_cellfrac_cohort_data(id_vegFrac,   any_tile, any_vegn)
  call send_cellfrac_cohort_data(id_treeFrac,  any_tile, is_tree)
  call send_cellfrac_cohort_data(id_c3pftFrac, any_tile, is_c3)
  call send_cellfrac_cohort_data(id_c4pftFrac, any_tile, is_c4)

  call send_cellfrac_cohort_data(id_grassFrac,   is_psl, is_grass)
  call send_cellfrac_cohort_data(id_grassFracC3, is_psl, is_c3grass)
  call send_cellfrac_cohort_data(id_grassFracC4, is_psl, is_c4grass)

  ! LUMIP land use fractions
  call send_cellfrac_data(id_fracLut_psl,  is_psl)
  call send_cellfrac_data(id_fracLut_crp,  is_crop)
  call send_cellfrac_data(id_fracLut_pst,  is_pst)
  call send_cellfrac_data(id_fracLut_urb,  is_urban)

  ! get components of calendar dates for this and previous time step
  call get_date(lnd%time-lnd%dt_slow, year1,month1,day1,hour,minute,second)
  call get_date(lnd%time,             year0,month0,day0,hour,minute,second)

  if (day0/=day1) then
     ! calculate daily average canopy air temperature
     ce = first_elmt(land_tile_map)
     do while(loop_over_tiles(ce,tile))
         if (associated(tile%vegn)) &
            tile%vegn%tc_daily = tile%vegn%tc_daily/steps_per_day
     enddo
  endif

  ! invoke any processes that potentially change tiling
  call vegn_nat_mortality_ppa( )
  call fire_transitions(lnd%time)
  call land_transitions(lnd%time)

  ! try to minimize the number of tiles by merging similar ones
  if (year0/=year1) then
     call land_tile_list_init(tmp)
     do l = lnd%ls,lnd%le
        ! merge all tiles into temporary list
        do while (.not.empty(land_tile_map(l)))
           ce=first_elmt(land_tile_map(l))
           tile=>current_tile(ce)
           call remove(ce)
           call merge_land_tile_into_list(tile,tmp)
        enddo
        ! move all tiles from temporary list to tile map
        do while (.not.empty(tmp))
           ce=first_elmt(tmp)
           tile=>current_tile(ce)
           call remove(ce)
           call insert(tile,land_tile_map(l))
        enddo
     enddo
     call land_tile_list_end(tmp)
  endif

  call update_vegn_slow( )
  ! send the accumulated diagnostics to the output
  call dump_tile_diag_fields(lnd%time)

  ! land_transitions may have changed the number of tiles per grid cell: reallocate
  ! boundary conditions, if necessary
  call realloc_cplr2land( cplr2land )
  call realloc_land2cplr( land2cplr )
  ! set the land mask to FALSE everywhere -- update_land_bc_fast will set it to
  ! true where necessary
  land2cplr%mask = .FALSE.
  land2cplr%tile_size = 0.0

  ! get the current state of the land boundary for the coupler
  ce = first_elmt(land_tile_map, ls=lnd%ls)
  do while(loop_over_tiles(ce,tile,l,k))
     call set_current_point(l,k)

     n_cohorts = 1; if (associated(tile%vegn)) n_cohorts = tile%vegn%n_cohorts
     call update_land_bc_fast (tile, n_cohorts, l,k, land2cplr, is_init=.true.)
  enddo

  call update_land_bc_slow( land2cplr )

  ! reset daily accumulated values
  if (day0/=day1) then
     ! reset daily air temperatures accumulators
     ce = first_elmt(land_tile_map)
     do while(loop_over_tiles(ce,tile))
         if(associated(tile%vegn)) then
            ! reset the accumulated values for the next day
            tile%vegn%tc_daily    =  0.0
            tile%vegn%daily_T_max = -HUGE(1.0)
            tile%vegn%daily_T_min =  HUGE(1.0)
         endif
     enddo
  endif

  call mpp_clock_end(landClock)
  call mpp_clock_end(landSlowClock)
end subroutine update_land_model_slow

! ============================================================================
! solve for surface temperature. ensure that melt does not exceed available
! snow or soil ice (which would create energy-balance problems). also ensure
! that melt and temperature are consistent and that evaporation from liquid
! soil water does not exceed exfiltration rate limit, if one is supplied.
! because the possible combinations of active constraints has multiplied
! greatly, we do not allow phase change for any surface (i.e., soil) at which
! we might apply a constraint on Eg.
subroutine land_surface_energy_balance ( &
     ! surface parameters
     grnd_T,          & ! ground temperature
     grnd_liq, grnd_ice, & ! amount of water available for freeze or melt on the surface, kg/m2
     grnd_latent,     & ! specific heat of vaporization for the ground
     grnd_Tf,         & ! ground freezing temperature
     grnd_E_min,      & ! Eg floor of 0 if condensation is prohibited
     grnd_E_max,      & ! exfiltration rate limit, kg/(m2 s)
     ! components of the ground energy balance linearization. Note that those
     ! are full derivatives, which include the response of the other surface
     ! scheme parameters to the change in ground temperature.
     fswg,            & ! net short-wave
     flwg0, DflwgDTg, DflwgDpsig, & ! net long-wave
     Hg0,   DHgDTg,   DHgDpsig,   & ! sensible heat
     Eg0,   DEgDTg,   DEgDpsig,   & ! latent heat
     G0,    DGDTg,    & ! sub-surface heat
     delta_Tg,        & ! surface temperature change for the time step
     delta_psig,      &
     Mg_imp          )  ! implicit melt, kg/m2

  real, intent(in) :: &
     grnd_T,          & ! ground temperature
     grnd_liq, grnd_ice, & ! amount of water available for freeze or melt on the surface
     grnd_latent,     & ! specific heat of vaporization for the ground
     grnd_Tf,         & ! ground freezing temperature
     grnd_E_min,      & ! Eg floor of 0 if condensation is prohibited
     grnd_E_max,      & ! exfiltration rate limit, kg/(m2 s)
     fswg,            & ! net short-wave
     flwg0, DflwgDTg, DflwgDpsig, & ! net long-wave
     Hg0,   DHgDTg,   DHgDpsig,   & ! sensible heat
     Eg0,   DEgDTg,   DEgDpsig,   & ! latent heat
     G0,    DGDTg                   ! sub-surface heat
  real, intent(out) :: &
     delta_Tg,        & ! change in surface temperature
     delta_psig,      & ! change in surface soil-water matric head
     Mg_imp             ! mass of surface ice melted (or water frozen) during the
                        ! time step, kg/m2

  real :: grnd_B     ! surface energy balance
  real :: grnd_DBDTg ! full derivative of grnd_B w.r.t. surface temperature
  real :: grnd_DBDpsig ! full derivative of grnd_B w.r.t. surface soil-water matric head
  real :: grnd_E_force, Eg_trial, Eg_check, determinant

  grnd_B      = fswg + flwg0      - Hg0      - grnd_latent*Eg0      - G0
  grnd_DBDTg  =        DflwgDTg   - DHgDTg   - grnd_latent*DEgDTg   - DGDTg
  grnd_DBDpsig =       DflwgDpsig - DHgDpsig - grnd_latent*DEgDpsig

  if(is_watch_point())then
     write(*,*)'#### ground balance input'
     __DEBUG1__(grnd_T)
     __DEBUG2__(grnd_liq, grnd_ice)
     __DEBUG1__(grnd_latent)
     __DEBUG1__(grnd_Tf)
     __DEBUG1__(grnd_E_min)
     __DEBUG1__(grnd_E_max)
     __DEBUG1__(fswg)
     __DEBUG3__(flwg0, DflwgDTg,DflwgDpsig)
     __DEBUG3__(Hg0,   DHgDTg,  DHgDpsig  )
     __DEBUG3__(Eg0,   DEgDTg,  DEgDpsig  )
     __DEBUG2__(G0,    DGDTg)
     write(*,*)'#### end of ground balance input'
     __DEBUG3__(grnd_B, grnd_DBDTg, grnd_DBDpsig)
  endif

  ! determine the ground temperature change under the assumptions that
  ! (1) no phase change occurs at surface (always true for soil now), and
  ! (2) delta_psig is zero (always true for snow, lake, glacier now)
  delta_Tg = - grnd_B/grnd_DBDTg
  delta_psig = 0
  ! calculate phase change on the ground, if necessary
  if     (grnd_ice>0.and.grnd_T+delta_Tg>grnd_Tf) then ! melt > 0
     Mg_imp =  min(grnd_ice,  grnd_DBDTg*(grnd_Tf-grnd_T-delta_Tg)*delta_time/hlf)
  elseif (grnd_liq>0.and.grnd_T+delta_Tg<grnd_Tf) then ! melt < 0
     Mg_imp = -min(grnd_liq, -grnd_DBDTg*(grnd_Tf-grnd_T-delta_Tg)*delta_time/hlf)
  else
     Mg_imp = 0
  endif
  ! adjust temperature change for the phase change
  delta_Tg = -(grnd_B - Mg_imp*hlf/delta_time)/grnd_DBDTg
  Eg_trial = Eg0 + DEgDTg*delta_Tg

  if(is_watch_point())then
     write(*,*)'#### ground balance solution with psig constant:'
     __DEBUG2__(grnd_B, grnd_DBDTg)
     __DEBUG3__(Mg_imp, delta_Tg, delta_psig)
     __DEBUG3__(grnd_E_min, Eg_trial, grnd_E_max)
  endif

  ! Solution above (assuming no change in psig) is acceptable if
  ! it does not imply unacceptable value of Eg. this is always
  ! true for lake, glacier, snow. If implied Eg is outside prescribed
  ! bounds (a possibility only for soil), then Eg is set at bound (grnd_E_force)
  ! and the required Tg and psig are found. To accomplish this,
  ! we solve the system
  ! grnd_B + grnd_DBDTg*delta_Tg + grnd_DBDpsig*delta_psig = 0
  ! Eg0    +     DEgDTg*delta_Tg +     DEgDpsig*delta_psig = grnd_E_force
  ! There is no need to revisit the solution for phase change, which
  ! is only done explicitly for soil.

  if (Eg_trial.lt.grnd_E_min .or. Eg_trial.gt.grnd_E_max) then
      grnd_E_force = max(grnd_E_min, Eg_trial)
      grnd_E_force = min(grnd_E_max, grnd_E_force)
      determinant = grnd_DBDTg*DEgDpsig-grnd_DBDpsig*DEgDTg
      delta_Tg   = - (DEgDpsig*grnd_B + grnd_DBDpsig*(grnd_E_force-Eg0))/determinant
      delta_psig =   (DEgDTg  *grnd_B + grnd_DBDTg  *(grnd_E_force-Eg0))/determinant
      Eg_check = Eg0 + DEgDTg*delta_Tg + DEgDpsig*delta_psig
      Mg_imp = 0
         if(is_watch_point())then
            write(*,*)'#### trial solution violated Eg limit, new soln:'
            __DEBUG2__(grnd_B,grnd_DBDTg)
            __DEBUG3__(Mg_imp,delta_Tg,delta_psig)
            __DEBUG3__(grnd_E_min, Eg_check, grnd_E_max)
         endif
    endif
end subroutine land_surface_energy_balance


! ============================================================================
! given downward long-wave flux from the atmosphere, and optical properties
! of vegetation and ground surface, calculates radiative balances of canopy
! layers and ground surface, and their derivatives w.r.t. temperatures
subroutine land_lw_balance(lwdn_atm, layer, frac, vegn_T, vegn_dT, surf_T, surf_dT, &
  vegn_tran_lw, vegn_refl_lw, surf_refl_lw, &
  flwv, flwg, DflwvDTv, DflwvDTg, DflwgDTv, DflwgDTg )

  real, intent(in) :: lwdn_atm ! downward long-wave radiation on top of the canopy, W/m2
  integer, intent(in) :: layer(:) ! layer number for cohort, top-down
  real, intent(in) :: frac(:)   ! fractional crown area of cohorts
  real, intent(in) :: vegn_T(:) ! canopy temperatures, deg K
  real, intent(in) :: vegn_dT(:)! expected change of canopy temperatures during this time step, deg K
  real, intent(in) :: surf_T    ! ground surface temperature, deg K
  real, intent(in) :: surf_dT   ! expected change of ground surface temperature during this time step, deg K
  real, intent(in) :: vegn_tran_lw(:) ! transmittance of cohort canopies to long-wave radiation
  real, intent(in) :: vegn_refl_lw(:) ! reflectance of cohort canopies to long-wave radiation
  real, intent(in) :: surf_refl_lw ! reflectance of ground surface to long-wave

  real, intent(out) :: flwv(:) ! long-wave balances of canopy layers, W/m2
  real, intent(out) :: flwg ! ground surface long-wave balance, W/m2
  real, intent(out) :: DflwvDTv(:,:) ! derivatives of canopy balances w.r.t. canopy temperatures, W/(m2 K)
  ! DflwvDTv(i,j) is the derivative of i-th cohort LW balance w.r.t j-th cohort temperature
  real, intent(out) :: DflwvDTg(:) ! derivatives of canopy balances w.r.t. ground surface temperature, W/(m2 K)
  real, intent(out) :: DflwgDTv(:) ! derivatives of ground surface balance w.r.t. canopy temperatures, W/(m2 K)
  real, intent(out) :: DflwgDTg    ! derivative of ground surface balance w.r.t. ground surface temperature, W/(m2 K)
  ! the radiative balances and their derivatives are per unit area of
  ! respective canopy

  ! ---- local vars
  integer :: N ! number of canopy layers
  integer :: M ! total number of cohorts
  real :: surf_emis_lw   ! surface emissivity = 1-surf_refl_lw
  real, allocatable :: &
     vegn_emis_lw(:),  & ! absorptivity of cohorts
     bbrad(:),         & ! emission from cohort
     layer_tran_lw(:), & ! average transmittance of layers
     layer_refl_lw(:), & ! average reflectance of layers
     layer_emis_lw(:), & ! average absorptivity of layers
     layer_area(:),    & ! total area of the crowns in the layer
     scale(:),         & ! scaling factor due to multiple reflections, 1/(1-refl*refl(i))
     B(:),             & ! black-body emission of individual layers, W/m2
     refl(:),          & ! integral reflectance below layer N, with multiple scattering
     emis(:),          & ! integral emission from all layers below N, with multiple scattering
     lwdn(:),          & ! downward long-wave flux below canopy layer, W/m2
     CL(:)               ! coefficients for radiative balance calculations
  integer :: i,j,k

  surf_emis_lw = 1-surf_refl_lw

  N = maxval(layer) ! number of layers
  M = size(vegn_T)  ! number of cohorts
  ! check argument shapes
!#define __CHECK_SIZE__(x)if(size(x)/=M) call error_mesg('land_lw_balance','Size of '//#x//' is incorrect',FATAL)
!  __CHECK_SIZE__(layer)
!  __CHECK_SIZE__(frac)
!  __CHECK_SIZE__(vegn_tran_lw)
!  __CHECK_SIZE__(vegn_refl_lw)
!  __CHECK_SIZE__(DflwvDTg)
!  __CHECK_SIZE__(DflwgDTv)
!#undef __CHECK_SIZE__
  if (size(DflwvDTv,1)/=M.or.size(DflwvDTv,2)/=M) &
    call error_mesg('land_lw_balance','Shape of DflwvDTv is incorrect',FATAL)

  ! allocate local variables
  allocate(layer_emis_lw(N),layer_tran_lw(N),layer_refl_lw(N),layer_area(N), &
           scale(N),B(N),refl(0:N),emis(0:N),lwdn(0:N),CL(N))
  allocate(vegn_emis_lw(M),bbrad(M))
  layer_emis_lw = 0 ; layer_tran_lw = 0 ; layer_refl_lw = 0 ; layer_area = 0 ;
  B = 0;
  ! [1] calculate some cohort and average layer properties
  do k = 1, M
    ! canopy radiative properties
    vegn_emis_lw(k) = 1-vegn_refl_lw(k)-vegn_tran_lw(k)
    bbrad(k) = stefan*vegn_emis_lw(k)*vegn_T(k)**4
    ! average layer properties
    i = layer(k)
    layer_tran_lw(i) = layer_tran_lw(i)+vegn_tran_lw(k)*frac(k)
    layer_refl_lw(i) = layer_refl_lw(i)+vegn_refl_lw(k)*frac(k)
    ! gray-body radiation emitted by the i-th layer of canopy
    B(i) = B(i) + bbrad(k)*frac(k)
    layer_area(i) = layer_area(i) + frac(k)
  enddo
  ! take into account possible gaps in the canopy layers
  do i = 1,N
    if (layer_area(i) < 1) then
       layer_tran_lw(i) = layer_tran_lw(i) + (1-layer_area(i))
    endif
  enddo

  ! [2] go upward through the canopy and calculate integral reflectances and emissions
  emis(N) = stefan*surf_emis_lw*surf_T**4
  refl(N) = surf_refl_lw
  do i = N,1,-1
     ! calculate absorptivity according to Kircchoff law
     layer_emis_lw(i) = 1-layer_refl_lw(i)-layer_tran_lw(i)
     ! common scale factor to account for multiple reflections
     scale(i) = 1.0/(1-layer_refl_lw(i)*refl(i))
     ! reflectance of all layers below i, including ground
     refl(i-1) = layer_refl_lw(i) + refl(i)*layer_tran_lw(i)**2*scale(i)
     ! emission from all layers below i, including ground
     emis(i-1) = B(i) + layer_tran_lw(i)*(emis(i)+B(i)*refl(i))*scale(i)
     ! coefficients for radiative balance calculations
     CL(i) = (1+layer_tran_lw(i)*refl(i)*scale(i))
  enddo
  ! downward lw flux below each layer
  lwdn(0) = lwdn_atm
  do i = 1, N
    lwdn(i)= (B(i) + layer_refl_lw(i)*emis(i) + layer_tran_lw(i)*lwdn(i-1))*scale(i)
  enddo

  ! [3] go down through the canopy calculates the long-waves balances
  ! TODO: check that cohorts are sorted in the order of layers
  do k = 1, M ! loop over cohorts
    i = layer(k)
    flwv(k) = -2*bbrad(k) + &
      vegn_emis_lw(k)*(B(i)*refl(i)*scale(i) + CL(i)*lwdn(i-1) + emis(i)*scale(i))
  enddo
  flwg = -emis(N) + lwdn(N)*(1-refl(N))

  ! [4] calculate the derivatives of vegetation and ground radiative balances
  ! w.r.t. vegetation temperatures
  do k = 1,M
     ! surrogate "emission" from the k-th cohort -- this is actually
     ! a derivative of gray-body radiation which is going to give us derivatives
     ! for other layers w.r.t k-th temperature due to linearity of the system.
     ! See radiation notes (end of multi-layer long-wave flux section) for
     ! details.
     bbrad(:) = 0 ; bbrad(k)    = 4*stefan*vegn_emis_lw(k)*(vegn_T(k)+0.5*vegn_dT(k))**3
     B(:)     = 0 ; B(layer(k)) = bbrad(k)*frac(k)
     ! recalculate integral emission for "surrogate emissions"
     emis(N) = 0.0
     do i = N,1,-1
        emis(i-1) = B(i) + layer_tran_lw(i)*(emis(i)+B(i)*refl(i))*scale(i)
     enddo
     ! calculate downward longwave below layer i; lwdn(0) is the flux from atmos
     lwdn(0) = 0.0
     do i = 1,N
        lwdn(i) = (B(i) + layer_refl_lw(i)*emis(i) + layer_tran_lw(i)*lwdn(i-1))*scale(i)
     enddo
     do j = 1,M
        ! calculate derivative of j-th cohort energy balance w.r.t. vegn_T(k)
        i = layer(j)
        DflwvDTv(j,k) = -2*bbrad(j) + &
            vegn_emis_lw(j)*(B(i)*refl(i)*scale(i) + CL(i)*lwdn(i-1) + emis(i)*scale(i))
     enddo
     DflwgDTv(k) = -emis(N) + lwdn(N)*(1-refl(N))
  enddo

  ! [5] calculate the derivatives of vegetation and ground radiative balances
  ! w.r.t. ground temperature
  ! surrogate "emission" for the ground
  bbrad(:)= 0 ; B(:) = 0
  emis(N) = 4*stefan*surf_emis_lw*(surf_T+0.5*surf_dT)**3
  do i = N,1,-1
     emis(i-1) = B(i) + layer_tran_lw(i)*(emis(i)+B(i)*refl(i))*scale(i)
  enddo
  lwdn(0) = 0.0
  do i = 1,N
     lwdn(i) = (B(i) + layer_refl_lw(i)*emis(i) + layer_tran_lw(i)*lwdn(i-1))*scale(i)
  enddo
  do k = 1,M
     ! calculate derivative of i-th cohort energy balance w.r.t. vegn_T(k)
     i = layer(k)
     DflwvDTg(k) = -2*bbrad(k) + &
        vegn_emis_lw(k)*(B(i)*refl(i)*scale(i) + CL(i)*lwdn(i-1) + emis(i)*scale(i))
  enddo
  DflwgDTg = -emis(N) + lwdn(N)*(1-refl(N))

  ! delete temporary arrays
  deallocate(vegn_emis_lw, bbrad, layer_tran_lw, layer_refl_lw, layer_emis_lw, &
     layer_area, scale, B, refl, emis, lwdn, CL)
end subroutine land_lw_balance

! ===========================================================================
! given direct and diffuse light on top of the canopy, and canopy optical
! properties, calculates net short-wave radiative balance for each canopy layer
! and underlying surface, and land albedo for direct and diffuse light
subroutine land_sw_balance ( &
  swdn_dif, swdn_dir, &
  vegn_layer, vegn_frac, &
  vegn_refl_dif, vegn_tran_dif, &
  vegn_refl_dir, vegn_sctr_dir, vegn_tran_dir, &
  surf_refl_dif, surf_refl_dir, &
  ! output:
  fswv, fswg, fswdn, &
  land_albedo_dif, land_albedo_dir )
  real, intent(in) :: swdn_dir ! downward direct radiation from atmos, W/m2
  real, intent(in) :: swdn_dif ! downward diffuse radiation from atmos, W/m2
  integer, intent(in) :: vegn_layer(:) ! layer number for each cohort, top-down
  real, intent(in) :: vegn_frac(:)     ! fractional crown area of cohorts
  real, intent(in) :: vegn_tran_dif(:) ! transmittances for diffuse beam
  real, intent(in) :: vegn_refl_dif(:) ! black-background reflectances for diffuse light
  real, intent(in) :: vegn_tran_dir(:) ! transmittances for direct beam
  real, intent(in) :: vegn_refl_dir(:) ! black-background reflectances for direct light
  real, intent(in) :: vegn_sctr_dir(:) ! downward scattering coefficients for direct beam
  real, intent(in) :: surf_refl_dir    ! ground surface albedo for direct light
  real, intent(in) :: surf_refl_dif    ! ground surface albedo for diffuse light

  real, intent(out) :: fswv(:) ! resulting radiative balances of canopy layers, W/m2
  real, intent(out) :: fswg    ! resulting radiative balance of ground surface, W/m2
  real, intent(out) :: fswdn(:) ! downward total flux on top of each cohort, W/m2
  real, intent(out),optional :: land_albedo_dir ! land albedo for direct light
  real, intent(out),optional :: land_albedo_dif ! land albedo for diffuse light

  ! ---- local vars
  integer :: N ! number of canopy layers
  integer :: M ! number of cohorts
  real, allocatable :: &
     scale(:), &    ! scaling factor due to multiple reflections
     refl_dir(:), & ! integral refl. below layer N for direct beam, with multiple scattering
     refl_dif(:)    ! integral refl. below layer N for diffuse light, with multiple scattering
  real, allocatable ::  & ! average optical properties of layers
     layer_tran_dif(:), & ! transmittances for diffuse beam
     layer_refl_dif(:), & ! black-background reflectances for diffuse light
     layer_refl_dir(:), & ! black-background reflectances for direct light
     layer_tran_dir(:), & ! transmittances for direct beam
     layer_sctr_dir(:), & ! downward scattering coefficients for direct beam
     layer_area(:)        ! sum of cohort fractional areas, must be close to 1 for all layers

  real :: dir, dif ! downward direct and diffuse light on top of the current layer, W/m2
  integer :: i,k

  ! verify argument shapes
  M = size(vegn_tran_dir)
!#define __CHECK_SIZE__(x)if(size(x)/=M) call error_mesg('land_sw_balance','Size of '//#x//' is incorrect',FATAL)
!  __CHECK_SIZE__(vegn_tran_dif)
!  __CHECK_SIZE__(vegn_refl_dir)
!  __CHECK_SIZE__(vegn_sctr_dir)
!  __CHECK_SIZE__(vegn_tran_dir)
!  __CHECK_SIZE__(fswv)
!  __CHECK_SIZE__(fswdn)
!#undef __CHECK_SIZE__

  ! allocate local variables
  N = maxval(vegn_layer)
  allocate(scale(N),refl_dir(0:N),refl_dif(0:N))
  allocate(layer_tran_dif(N), layer_refl_dif(N), layer_refl_dir(N), &
           layer_tran_dir(N), layer_sctr_dir(N), layer_area    (N)  )
  ! calculate optical properties of layers
  ! layer_area(i) must be 1 (to the numerical precision) for all i
  ! optical properties of cohorts must be pre-normalized
  layer_tran_dif = 0 ; layer_refl_dif = 0 ; layer_refl_dir = 0
  layer_tran_dir = 0 ; layer_sctr_dir = 0 ; layer_area     = 0
  do k = 1,M ! loop over cohorts
     i = vegn_layer(k)
     layer_tran_dif(i) = layer_tran_dif(i)+vegn_tran_dif(k)*vegn_frac(k)
     layer_refl_dif(i) = layer_refl_dif(i)+vegn_refl_dif(k)*vegn_frac(k)
     layer_refl_dir(i) = layer_refl_dir(i)+vegn_refl_dir(k)*vegn_frac(k)
     layer_tran_dir(i) = layer_tran_dir(i)+vegn_tran_dir(k)*vegn_frac(k)
     layer_sctr_dir(i) = layer_sctr_dir(i)+vegn_sctr_dir(k)*vegn_frac(k)
     layer_area(i)     = layer_area(i)+vegn_frac(k)
  enddo
  ! take into account possible gaps in the canopy
  do i = 1,N
     if (layer_area(i) < 1) then
        layer_tran_dif(i) = layer_tran_dif(i) + (1-layer_area(i))
        layer_tran_dir(i) = layer_tran_dir(i) + (1-layer_area(i))
     endif
  enddo

  ! [1] go upward through the canopy and calculate integral reflectances
  refl_dir(N) = surf_refl_dir
  refl_dif(N) = surf_refl_dif
  do i = N,1,-1
    scale(i) = 1.0/(1 - refl_dif(i)*layer_refl_dif(i))
    refl_dir(i-1) = layer_refl_dir(i) &
      + layer_tran_dif(i)*(refl_dif(i)*layer_sctr_dir(i)+refl_dir(i)*layer_tran_dir(i))&
      * scale(i)
    refl_dif(i-1) = layer_refl_dif(i) &
      + refl_dif(i)*layer_tran_dif(i)**2*scale(i)
  enddo
  ! assign land albedo values, if necessary
  if (present(land_albedo_dir)) land_albedo_dir = refl_dir(0)
  if (present(land_albedo_dif)) land_albedo_dif = refl_dif(0)

  ! [2] go down through the canopy and calculate radiative balances
  ! TODO : assume that the cohorts are sorted top-down and optimize the loop
  ! using this property. Also check that the cohorts are sorted.
  dir = swdn_dir
  dif = swdn_dif
  do i = 1,N
     do k = 1,M ! loop over cohorts
        if (vegn_layer(k)/=i) cycle
        ! assign downward fluxes on top of this cohort
        fswdn(k) = dif + dir
        ! calculate radiative balance of the current cohort
        fswv(k) = &
          dif * (1-vegn_refl_dif(k)-vegn_tran_dif(k)) &
              * (1+layer_tran_dif(i)*refl_dif(i) * scale(i)) &
        + dir * (1-vegn_tran_dir(k)-vegn_refl_dir(k)-vegn_sctr_dir(k)) &
        + dir * (1-vegn_tran_dif(k)-vegn_refl_dif(k)) &
              * (refl_dif(i)*layer_sctr_dir(i)+refl_dir(i)*layer_tran_dir(i))*scale(i)
     enddo
     ! recalculate the fluxes for the lower layer
     dif = (layer_sctr_dir(i)+layer_refl_dif(i)*refl_dir(i)*layer_tran_dir(i))*scale(i)*dir &
         + layer_tran_dif(i)*scale(i)*dif
     dir = layer_tran_dir(i)*dir
  enddo
  fswg = (1-surf_refl_dif)*dif + (1-surf_refl_dir)*dir

  ! deallocate local variables
  deallocate(scale,refl_dir,refl_dif)
  deallocate(layer_tran_dif, layer_refl_dif, layer_refl_dir, &
             layer_tran_dir, layer_sctr_dir, layer_area      )
end subroutine land_sw_balance


! ============================================================================
! calculate fractions of downward short-wave radiation absorbed by vegetation
! and surface, and land albedo values.
subroutine land_sw_radiation (     &
     subs_refl_dir, subs_refl_dif, &
     snow_refl_dir, snow_refl_dif, &
     snow_area, &
     vegn_layer,    vegn_frac, &
     vegn_refl_dif, vegn_tran_dif, &
     vegn_refl_dir, vegn_sctr_dir, vegn_tran_dir, &
     ! output:
     Sg_dir, Sg_dif, Sv_dir, Sv_dif, Sdn_dir, Sdn_dif, &
     land_albedo_dir, land_albedo_dif )

  real, intent(in) :: &
     subs_refl_dir(NBANDS), subs_refl_dif(NBANDS), & ! sub-snow reflectances for direct and diffuse light
     snow_refl_dir(NBANDS), snow_refl_dif(NBANDS), & ! snow reflectances for direct and diffuse light
     snow_area, &
     ! black-background vegetation radiative properties, per cohort and band.
     ! dimensions of the arrays below are (NCOHORTS,NBANDS)
     vegn_refl_dif(:,:), & ! reflectance for diffuse light
     vegn_tran_dif(:,:), & ! transmittance for diffuse light
     vegn_refl_dir(:,:), & ! reflectance (scattering upwards) for direct light
     vegn_sctr_dir(:,:), & ! downward scattering coefficient direct beam
     vegn_tran_dir(:,:)    ! transmittance for direct beam
  integer, intent(in) :: vegn_layer(:) ! layer number for each cohort, top-down
  real, intent(in)    :: vegn_frac(:)  ! fractional crown area of each cohort

  real, intent(out) :: &
     Sg_dir(NBANDS), Sg_dif(NBANDS), & ! fraction of downward short-wave absorbed by ground and snow
     Sv_dir(:,:),    Sv_dif(:,:),    & ! fraction of downward short-wave absorbed by vegetation (NCOHORTS,NBANDS)
     Sdn_dir(:,:),   Sdn_dif(:,:),   & ! fraction of downward short-wave on top of each cohort (NCOHORTS,NBANDS)
     land_albedo_dir(NBANDS), land_albedo_dif(NBANDS) ! land albedo for direct and diffuse light

  ! ---- local vars
  real :: &
     grnd_refl_dir(NBANDS), & ! SW reflectances of ground surface (by spectral band)
     grnd_refl_dif(NBANDS)    ! SW reflectances of ground surface (by spectral band)
  integer :: band, m
  real, allocatable :: layerfrac(:)

  grnd_refl_dir = subs_refl_dir + (snow_refl_dir - subs_refl_dir) * snow_area
  grnd_refl_dif = subs_refl_dif + (snow_refl_dif - subs_refl_dif) * snow_area

  do band = 1, NBANDS
     ! diffuse radiation
     call land_sw_balance(1.0, 0.0, &
        vegn_layer, vegn_frac, &
        vegn_refl_dif(:,band), vegn_tran_dif(:,band), &
        vegn_refl_dir(:,band), vegn_sctr_dir(:,band), vegn_tran_dir(:,band), &
        grnd_refl_dif(band), grnd_refl_dir(band),&
        Sv_dif(:,band), Sg_dif(band), Sdn_dif(:,band), &
        land_albedo_dif=land_albedo_dif(band) )
     ! direct radiation
     call land_sw_balance(0.0, 1.0, &
        vegn_layer, vegn_frac, &
        vegn_refl_dif(:,band), vegn_tran_dif(:,band), &
        vegn_refl_dir(:,band), vegn_sctr_dir(:,band), vegn_tran_dir(:,band), &
        grnd_refl_dif(band), grnd_refl_dir(band),&
        Sv_dir(:,band), Sg_dir(band), Sdn_dir(:,band), &
        land_albedo_dir=land_albedo_dir(band) )
  enddo

  if(is_watch_point()) then
     write(*,*)' #### land_sw_radiation ####'
     allocate(layerfrac(maxval(vegn_layer)))
     layerfrac = 0.0
     do m = 1,size(vegn_layer)
        layerfrac(vegn_layer(m)) = layerfrac(vegn_layer(m)) + vegn_frac(m)
     enddo
     do m = 1, size(layerfrac)
        write(*,'("layer ",i2.2)',advance='NO') m
        call dpri('layerfrac',layerfrac(m))
        write(*,*)
     enddo
     deallocate(layerfrac)

     call debug_rad_properties(BAND_VIS,'(VIS)')
     call debug_rad_properties(BAND_NIR,'(NIR)')
  endif

contains
  subroutine debug_rad_properties(b,tag)
     integer, intent(in) :: b ! band number
     character(*), intent(in) :: tag

     integer :: m ! cohort iterator

     ! write(*,*)
     write(*,*)' #### surface radiative properties '//tag//' ####'
     call dpri('subs_refl_dif',subs_refl_dif(b))
     call dpri('subs_refl_dir',subs_refl_dir(b))
     write(*,*)
     call dpri('snow_refl_dif',snow_refl_dif(b))
     call dpri('snow_refl_dir',snow_refl_dir(b))
     __DEBUG1__(snow_area)
     call dpri('grnd_refl_dif',grnd_refl_dif(b))
     call dpri('grnd_refl_dir',grnd_refl_dir(b))
     write(*,*)
     write(*,*)' #### black-background vegn radiative properties '//tag//' ####'
     do m = 1,size(vegn_refl_dif,1)
        write(*,'(i2.2," : layer ",i2.2)',advance='NO') m, vegn_layer(m)
        call dpri('refl_dif',vegn_refl_dif(m,b))
        call dpri('tran_dif',vegn_tran_dif(m,b))
        call dpri('refl_dir',vegn_refl_dir(m,b))
        call dpri('tran_dir',vegn_tran_dir(m,b))
        call dpri('sctr_dir',vegn_sctr_dir(m,b))
        call dpri('sctr_dir',vegn_sctr_dir(m,b))
        write(*,*)
     enddo

     write(*,*)' #### grnd radiation coefficients '//tag//' ####'
     call dpri('Sg_dif',Sg_dif(b))
     call dpri('Sg_dir',Sg_dir(b))
     write(*,*)
     write(*,*)' #### vegn radiation coefficients '//tag//' ####'
     do m = 1, size(vegn_refl_dif,1)
        write(*,'(i2.2," : layer ",i2.2)',advance='NO') m, vegn_layer(m)
        call dpri('Sv_dif',Sv_dif(m,b))
        call dpri('Sv_dir',Sv_dir(m,b))
        call dpri('Sdn_dif',Sdn_dif(m,b))
        call dpri('Sdn_dir',Sdn_dir(m,b))
        write(*,*)
     enddo
     write(*,*)' #### overall land reflectances '//tag//' ####'
     call dpri('land_refl_dif',land_albedo_dif(b))
     call dpri('land_refl_dir',land_albedo_dir(b))
     write(*,*)
  end subroutine debug_rad_properties
end subroutine land_sw_radiation

subroutine realloc1(x,N)
  real, allocatable, intent(inout) :: x(:)
  integer, intent(in) :: N

  if(allocated(x)) then
     if (size(x)==N) then
        return
     else
        deallocate(x)
     endif
  endif

  allocate(x(N))
end subroutine

subroutine realloc2(x,N)
  real, allocatable, intent(inout) :: x(:,:)
  integer, intent(in) :: N

  if(allocated(x)) then
     if (size(x,1)==N) then
        return
     else
        deallocate(x)
     endif
  endif

  allocate(x(N,NBANDS))
end subroutine

! ============================================================================
subroutine update_land_bc_fast (tile, N, l,k, land2cplr, is_init)
  type(land_tile_type), intent(inout) :: tile
  integer             , intent(in) :: N ! number of cohorts, 1 if no vegetation
  integer             , intent(in) :: l,k
  type(land_data_type), intent(inout) :: land2cplr
  logical, optional :: is_init
  ! Note that N is calculated and passed down to update_land_bc_fast simply
  ! for convenience, so that there is no need to make a lot of by-cohort arrays
  ! allocatable -- instead they are created on stack with the size passed
  ! as an argument.

  ! ---- local vars
  real :: &
         grnd_T, subs_z0m, subs_z0s, &
                 snow_z0s, snow_z0m, &
         snow_area, snow_depth

  real :: subs_refl_dir(NBANDS), subs_refl_dif(NBANDS) ! direct and diffuse albedos
  real :: subs_refl_lw ! reflectance for thermal radiation
  real :: snow_refl_dir(NBANDS), snow_refl_dif(NBANDS) ! direct and diffuse albedos of snow
  real :: snow_refl_lw ! snow reflectance for thermal radiation
  real :: snow_emis ! snow emissivity
  real :: grnd_emis ! ground emissivity
  ! NOTE :  grnd_emis is used only to satisfy xxxx_radiation interfaces; its value is ignored, but
  ! 1-refl is used instead. snow_emis is used in the the vegn_radiation, but it should not be since
  ! properties of intercepted snowpack are, in general, different from the snow on the ground
  real :: snow_area_rad ! "snow area for radiation calculations" -- introduced
                        ! to reproduce lm2 behavior
  real :: &
     ! long-wave black-background radiative properties, per cohort
     vegn_refl_lw(N), vegn_tran_lw(N), & ! reflectance and transmittance, respectively
     ! short-wave black-background vegetation radiative properties, per cohort and band.
     ! dimensions of the arrays below are (NCOHORTS,NBANDS)
     vegn_refl_dif(N,NBANDS), & ! reflectance for diffuse light
     vegn_tran_dif(N,NBANDS), & ! transmittance for diffuse light
     vegn_refl_dir(N,NBANDS), & ! reflectance (scattering upwards) for direct light
     vegn_sctr_dir(N,NBANDS), & ! downward scattering coefficient direct beam
     vegn_tran_dir(N,NBANDS)    ! transmittance for direct beam
  integer :: &
     vegn_layer(N) ! number of layer that respective cohort belongs to
  real :: &
     vegn_frac(N) ! fraction of layer covered by respective cohort canopy, unitless
  real :: &
     vegn_Tv,     &
     vegn_cover,  &
     vegn_height, vegn_lai, vegn_sai
  logical :: do_update

  real :: cosz    ! cosine of solar zenith angle
  real :: fracday ! daytime fraction of time interval
  real :: rrsun   ! earth-sun distance (r) relative to semi-major axis
                  ! of orbital ellipse (a) : (a/r)**2
  integer :: face ! for debugging
  integer :: i, j, m, tr

  i = lnd%i_index(l) ; j = lnd%j_index(l)

  vegn_Tv = 0

  do_update = .not.present(is_init)

  ! on initialization the albedos are calculated for the current time step ( that is, interval
  ! lnd%time, lnd%time+lnd%dt_fast); in the course of the run this subroutine is called
  ! at the end of time step (but before time is advanced) to calculate the radiative properties
  ! for the _next_ time step
  if (do_update) then
     call diurnal_solar(lnd%ug_lat(l), lnd%ug_lon(l), lnd%time+lnd%dt_fast, &
          cosz, fracday, rrsun, lnd%dt_fast)
  else
     call diurnal_solar(lnd%ug_lat(l), lnd%ug_lon(l), lnd%time, &
          cosz, fracday, rrsun, lnd%dt_fast)
  endif

  if (associated(tile%glac)) then
     call glac_radiation(tile%glac, cosz, subs_refl_dir, subs_refl_dif, subs_refl_lw, grnd_emis)
     call glac_roughness(tile%glac, subs_z0s, subs_z0m )
  else if (associated(tile%lake)) then
     call lake_radiation(tile%lake, cosz, subs_refl_dir, subs_refl_dif, subs_refl_lw, grnd_emis)
     call lake_roughness(tile%lake, subs_z0s, subs_z0m )
  else if (associated(tile%soil)) then
     call soil_radiation(tile%soil, cosz, subs_refl_dir, subs_refl_dif, subs_refl_lw, grnd_emis)
     call soil_roughness(tile%soil, subs_z0s, subs_z0m )
  else
     call get_current_point(face=face)
     call error_mesg('update_land_bc_fast','none of the surface tiles exist at ('//&
             trim(string(i))//','//trim(string(j))//','//trim(string(k))//&
             ', face='//trim(string(face))//')',FATAL)
  endif

  call snow_radiation ( tile%snow%T(1), cosz, associated(tile%glac), snow_refl_dir, snow_refl_dif, snow_refl_lw, snow_emis)
  call snow_get_depth_area ( tile%snow, snow_depth, snow_area )
  call snow_roughness ( tile%snow, snow_z0s, snow_z0m )

  ! allocate storage in the land tile, to carry the values calculated here to
  ! update_land_bc_fast
  call realloc2(tile%Sv_dir,N)
  call realloc2(tile%Sv_dif,N)
  call realloc2(tile%Sdn_dir,N)
  call realloc2(tile%Sdn_dif,N)
  call realloc1(tile%vegn_refl_lw,N)
  call realloc1(tile%vegn_tran_lw,N)

  if (associated(tile%vegn)) then
     call update_derived_vegn_data(tile%vegn, tile%soil)
     ! USE OF SNOWPACK RAD PROPERTIES FOR INTERCEPTED SNOW IS ERRONEOUS,
     ! NEEDS TO BE CHANGED. TEMPORARY.
     call vegn_radiation ( tile%vegn, cosz, snow_depth, snow_refl_dif, snow_emis, &
                   vegn_refl_dif, vegn_tran_dif, &
                   vegn_refl_dir, vegn_sctr_dir, vegn_tran_dir, &
                   vegn_refl_lw, vegn_tran_lw)
     ! (later see if we can remove vegn_cover from c-a-radiation...) TEMPORARY
     ! vegn_diffusion returns integral properties of the canopy, relevant for the
     ! calculations of the land roughness and displacement
     call vegn_diffusion ( tile%vegn, snow_depth, &
                   vegn_cover, vegn_height, vegn_lai, vegn_sai)
     ! assign layers and fractions
     vegn_layer(:) = tile%vegn%cohorts(1:N)%layer
     vegn_frac (:) = tile%vegn%cohorts(1:N)%layerfrac
  else
     ! set radiative properties for null vegetation
     vegn_refl_dif = 0
     vegn_tran_dif = 1
     vegn_refl_dir = 0
     vegn_sctr_dir = 0
     vegn_tran_dir = 1
     vegn_refl_lw  = 0
     vegn_tran_lw  = 1
     ! set cover for null vegetation
     vegn_cover    = 0
     ! set other parameters for null vegetation
     vegn_height   = 0
     vegn_lai      = 0
     vegn_sai      = 0

     vegn_layer    = 1
     vegn_frac     = 1.0
  endif

  ! store the values of long-wave optical properties to be used in the update_land_model_fast
  tile%surf_refl_lw = subs_refl_lw  + (snow_refl_lw  - subs_refl_lw ) * snow_area
  tile%vegn_refl_lw = vegn_refl_lw
  tile%vegn_tran_lw = vegn_tran_lw


  if(is_watch_point()) then
     write(*,*) '#### update_land_bc_fast ### checkpoint 1 ####'
     __DEBUG3__(cosz, fracday, rrsun)
     __DEBUG2__(vegn_lai,vegn_sai)
     write(*,*)' #### vegn radiative properties (LW) ####'
     do m = 1,N
        write(*,'(i2.2,2x)',advance='NO') m
        call dpri('refl_lw',vegn_refl_lw(m))
        call dpri('tran_lw',vegn_tran_lw(m))
        write(*,*)
     enddo
  endif

  snow_area_rad = snow_area
  if (lm2) then
     if(associated(tile%glac)                    ) snow_area_rad = 0
     if(associated(tile%soil).and.vegn_cover>0.01) snow_area_rad = 1
  endif

  call land_sw_radiation( &
       subs_refl_dir, subs_refl_dif, &
       snow_refl_dir, snow_refl_dif, &
       snow_area_rad,  &
       vegn_layer, vegn_frac, &
       vegn_refl_dif, vegn_tran_dif, &
       vegn_refl_dir, vegn_sctr_dir, vegn_tran_dir,  &
       ! output:
       tile%Sg_dir, tile%Sg_dif, tile%Sv_dir, tile%Sv_dif, tile%Sdn_dir, tile%Sdn_dif, &
       tile%land_refl_dir, tile%land_refl_dif )

  call cana_roughness( lm2, &
     subs_z0m, subs_z0s, &
     snow_z0m, snow_z0s, snow_area, &
     vegn_cover,  vegn_height, vegn_lai, vegn_sai, &
     tile%land_d, tile%land_z0m, tile%land_z0s)

  if(is_watch_point()) then
     __DEBUG1__(tile%land_z0m)
  endif

  land2cplr%t_surf         (l,k) = tfreeze
  land2cplr%t_ca           (l,k) = tfreeze
  land2cplr%tr             (l,k, isphum) = 0.0
  land2cplr%albedo         (l,k) = 0.0
  land2cplr%albedo_vis_dir (l,k) = 0.0
  land2cplr%albedo_nir_dir (l,k) = 0.0
  land2cplr%albedo_vis_dif (l,k) = 0.0
  land2cplr%albedo_nir_dif (l,k) = 0.0
  land2cplr%rough_mom      (l,k) = 0.1
  land2cplr%rough_heat     (l,k) = 0.1

  ! Calculate radiative surface temperature. lwup cannot be calculated here
  ! based on the available temperatures because it is a result of the implicit
  ! time step: lwup = lwup0 + DlwupDTg*delta_Tg + ..., so we have to carry it
  ! from the update_land_fast
  ! Consequence: since update_landbc_fast is called once at the very beginning of
  ! every run (before update_land_fast is called) lwup from the previous step
  ! must be stored in the in the restart for reproducibility
  land2cplr%t_surf(l,k) = ( tile%lwup/stefan ) ** 0.25

  if (associated(tile%glac)) call glac_get_sfc_temp(tile%glac, grnd_T)
  if (associated(tile%lake)) call lake_get_sfc_temp(tile%lake, grnd_T)
  if (associated(tile%soil)) call soil_get_sfc_temp(tile%soil, grnd_T)
  if (snow_area > 0)         call snow_get_sfc_temp(tile%snow, grnd_T)

  ! set the boundary conditions for the flux exchange
  land2cplr%mask           (l,k) = .TRUE.
  land2cplr%tile_size      (l,k) = tile%frac

  land2cplr%t_ca(l,k) = tile%cana%T
  do tr = 1,ntcana
     land2cplr%tr(l,k,tr) = tile%cana%tr(tr)
  enddo

  land2cplr%albedo_vis_dir (l,k) = tile%land_refl_dir(BAND_VIS)
  land2cplr%albedo_nir_dir (l,k) = tile%land_refl_dir(BAND_NIR)
  land2cplr%albedo_vis_dif (l,k) = tile%land_refl_dif(BAND_VIS)
  land2cplr%albedo_nir_dif (l,k) = tile%land_refl_dif(BAND_NIR)
  land2cplr%albedo         (l,k) = SUM(tile%land_refl_dir + tile%land_refl_dif)/4 ! incorrect, replace with proper weighting later
  land2cplr%rough_mom      (l,k) = tile%land_z0m
  land2cplr%rough_heat     (l,k) = tile%land_z0s

  if(is_watch_point()) then
     write(*,*)'#### update_land_bc_fast ### output ####'
     write(*,*)'land2cplr%mask',land2cplr%mask(l,k)
     write(*,*)'land2cplr%tile_size',land2cplr%tile_size(l,k)
     write(*,*)'land2cplr%t_surf',land2cplr%t_surf(l,k)
     write(*,*)'land2cplr%t_ca',land2cplr%t_ca(l,k)
     write(*,*)'land2cplr%albedo',land2cplr%albedo(l,k)
     write(*,*)'land2cplr%rough_mom',land2cplr%rough_mom(l,k)
     write(*,*)'land2cplr%rough_heat',land2cplr%rough_heat(l,k)
     write(*,*)'land2cplr%tr',land2cplr%tr(l,k,:)
     write(*,*)'#### update_land_bc_fast ### end of output ####'
  endif

  ! ---- diagnostic section
  if (id_vegn_cover_1 > 0) &
     call send_tile_data(id_vegn_cover_1, sum(vegn_frac(:),mask=(vegn_layer==1)), tile%diag)
  if (id_vegn_cover_U > 0) &
     call send_tile_data(id_vegn_cover_U, sum(vegn_frac(:),mask=(vegn_layer>1)), tile%diag)
  call send_tile_data(id_vegn_cover, vegn_cover, tile%diag)
  call send_tile_data(id_cosz, cosz, tile%diag)
  call send_tile_data(id_albedo_dir, tile%land_refl_dir, tile%diag)
  call send_tile_data(id_albedo_dif, tile%land_refl_dif, tile%diag)
  call send_tile_data(id_snow_frac,  snow_area,          tile%diag)
!!$  call send_tile_data(id_vegn_refl_dir, vegn_refl_dir,     tile%diag)
!!$  call send_tile_data(id_vegn_refl_dif, vegn_refl_dif, tile%diag)
!!$  call send_tile_data(id_vegn_refl_lw,  vegn_refl_lw, tile%diag)
!!$  call send_tile_data(id_vegn_tran_dir, vegn_tran_dir, tile%diag)
!!$  call send_tile_data(id_vegn_tran_dif, vegn_tran_dif, tile%diag)
!!$  call send_tile_data(id_vegn_tran_lw,  vegn_tran_lw, tile%diag)
!!$  call send_tile_data(id_vegn_sctr_dir, vegn_sctr_dir,     tile%diag)
  call send_tile_data(id_subs_refl_dir, subs_refl_dir, tile%diag)
  call send_tile_data(id_subs_refl_dif, subs_refl_dif, tile%diag)
  call send_tile_data(id_grnd_T,     grnd_T,     tile%diag)

  ! CMOR variables
  call send_tile_data(id_snd, max(snow_depth,0.0),     tile%diag)

  ! --- debug section
  call check_temp_range(land2cplr%t_ca(l,k),'update_land_bc_fast','T_ca')

end subroutine update_land_bc_fast


! ============================================================================
subroutine update_land_bc_slow (land2cplr)
  type(land_data_type), intent(inout) :: land2cplr

  ! ---- local vars
  integer :: i,j,k,l,face ! coordinates of the watch point, for debug printout

  call update_topo_rough(land2cplr%rough_scale)
  where (land2cplr%mask) &
       land2cplr%rough_scale = max(land2cplr%rough_mom,land2cplr%rough_scale)
  call get_watch_point(i,j,k,face,l)
  if ( lnd%ug_face==face.and.                &
       lnd%ls<=l.and.l<=lnd%le.and.          &
       k<=size(land2cplr%rough_scale,2).and. &
       is_watch_time()) then
     write(*,*)'#### update_land_bc_slow ### output ####'
     write(*,*)'land2cplr%rough_scale',land2cplr%rough_scale(l,k)
     write(*,*)'#### update_land_bc_slow ### end of output ####'
  endif

end subroutine update_land_bc_slow


! ============================================================================
subroutine Lnd_stock_pe(bnd,index,value)

type(land_data_type), intent(in)  :: bnd
integer             , intent(in)  :: index
real                , intent(out) :: value ! Domain water (Kg) or heat (Joules)

type(land_tile_enum_type)     :: ce
type(land_tile_type), pointer :: tile
character(len=128) :: message
integer :: l
real :: area_factor, river_value
! *_twd are tile water densities (kg water per m2 of tile)
real twd_gas_cana,               twd_liq_glac, twd_sol_glac, &
     twd_liq_lake, twd_sol_lake, twd_liq_soil, twd_sol_soil, &
     twd_liq_snow, twd_sol_snow, twd_liq_vegn, twd_sol_vegn
! *_gcwd are grid-cell water densities (kg water per m2 of land in grid cell)
real gcwd_cana, gcwd_glac, gcwd_lake, gcwd_soil, gcwd_snow, gcwd_vegn
! v_* are global masses of water
real v_cana, v_glac, v_lake, v_soil, v_snow, v_vegn
real a_globe

value = 0.0
v_cana = 0.
v_glac = 0.
v_lake = 0.
v_soil = 0.
v_snow = 0.
v_vegn = 0.
if(.not.bnd%pe) return

! The following is a dirty getaround
if(lnd%ug_cellarea(lnd%ls) < 1.0) then
  area_factor = 4*pi*radius**2 ! lnd%area is fraction of globe
else
  area_factor = 1.0 ! lnd%area is actual area (m**2)
endif

select case(index)
case(ISTOCK_WATER)
  do l = lnd%ls, lnd%le
    ce = first_elmt(land_tile_map(l))
    gcwd_cana = 0.0; gcwd_glac = 0.0; gcwd_lake = 0.0
    gcwd_soil = 0.0; gcwd_snow = 0.0; gcwd_vegn = 0.0
    do while(loop_over_tiles(ce,tile))
      twd_gas_cana = 0.0
      twd_liq_glac = 0.0 ; twd_sol_glac = 0.0
      twd_liq_lake = 0.0 ; twd_sol_lake = 0.0
      twd_liq_soil = 0.0 ; twd_sol_soil = 0.0
      twd_liq_snow = 0.0 ; twd_sol_snow = 0.0
      twd_liq_vegn = 0.0 ; twd_sol_vegn = 0.0
      if(associated(tile%cana)) then
         twd_gas_cana = canopy_air_mass*tile%cana%tr(isphum)
      endif
      if(associated(tile%glac)) &
         call glac_tile_stock_pe(tile%glac, twd_liq_glac, twd_sol_glac)
      if(associated(tile%lake)) &
         call lake_tile_stock_pe(tile%lake, twd_liq_lake, twd_sol_lake)
      if(associated(tile%soil)) &
         call soil_tile_stock_pe(tile%soil, twd_liq_soil, twd_sol_soil)
      if(associated(tile%snow)) &
         call snow_tile_stock_pe(tile%snow, twd_liq_snow, twd_sol_snow)
      if(associated(tile%vegn)) &
         call vegn_tile_stock_pe(tile%vegn, twd_liq_vegn, twd_sol_vegn)
      gcwd_cana = gcwd_cana +  twd_gas_cana                 * tile%frac
      gcwd_glac = gcwd_glac + (twd_liq_glac + twd_sol_glac) * tile%frac
      gcwd_lake = gcwd_lake + (twd_liq_lake + twd_sol_lake) * tile%frac
      gcwd_soil = gcwd_soil + (twd_liq_soil + twd_sol_soil) * tile%frac
      gcwd_snow = gcwd_snow + (twd_liq_snow + twd_sol_snow) * tile%frac
      gcwd_vegn = gcwd_vegn + (twd_liq_vegn + twd_sol_vegn) * tile%frac
    enddo
    v_cana = v_cana + gcwd_cana * lnd%ug_area(l)*area_factor
    v_glac = v_glac + gcwd_glac * lnd%ug_area(l)*area_factor
    v_lake = v_lake + gcwd_lake * lnd%ug_area(l)*area_factor
    v_soil = v_soil + gcwd_soil * lnd%ug_area(l)*area_factor
    v_snow = v_snow + gcwd_snow * lnd%ug_area(l)*area_factor
    v_vegn = v_vegn + gcwd_vegn * lnd%ug_area(l)*area_factor
  enddo
  value  = v_cana + v_glac + v_lake + v_soil + v_snow + v_vegn
a_globe = 4. * pi * radius**2
case(ISTOCK_HEAT)
  if(.not.stock_warning_issued) then
    call error_mesg('Lnd_stock_pe','Heat stock not yet implemented',NOTE)
    stock_warning_issued = .true.
  endif
case(ISTOCK_SALT)
  if(.not.stock_warning_issued) then
    call error_mesg('Lnd_stock_pe','Salt stock not yet implemented',NOTE)
    stock_warning_issued = .true.
  endif
case default
  write(message,'(i2,a,i2,a,i2,a,i2,a)') &
  index,' is an invalid stock index. Must be ISTOCK_WATER or ISTOCK_HEAT or ISTOCK_SALT (',ISTOCK_WATER,' or ',ISTOCK_HEAT,' or ', ISTOCK_SALT,')'
  call error_mesg('Lnd_stock_pe',message,FATAL)
end select

call river_stock_pe(index, river_value)
value = value + river_value

if (index.eq.ISTOCK_WATER.and.give_stock_details) then
    call mpp_sum(river_value, pelist=lnd%pelist)
    call mpp_sum(v_cana, pelist=lnd%pelist)
    call mpp_sum(v_glac, pelist=lnd%pelist)
    call mpp_sum(v_lake, pelist=lnd%pelist)
    call mpp_sum(v_soil, pelist=lnd%pelist)
    call mpp_sum(v_snow, pelist=lnd%pelist)
    call mpp_sum(v_vegn, pelist=lnd%pelist)
    write (message,'(a,f10.5)') 'total land storage:',v_cana/a_globe+v_glac/a_globe+ &
        v_lake/a_globe+v_soil/a_globe+v_snow/a_globe+v_vegn/a_globe+river_value/a_globe
    call error_mesg('Lnd_stock_pe',message,NOTE)
    write (message,'(a,7f10.5)') '...canopy air:',v_cana/a_globe
    call error_mesg('Lnd_stock_pe',message,NOTE)
    write (message,'(a,7f10.5)') '......glacier:',v_glac/a_globe
    call error_mesg('Lnd_stock_pe',message,NOTE)
    write (message,'(a,7f10.5)') '.........lake:',v_lake/a_globe
    call error_mesg('Lnd_stock_pe',message,NOTE)
    write (message,'(a,7f10.5)') '.........soil:',v_soil/a_globe
    call error_mesg('Lnd_stock_pe',message,NOTE)
    write (message,'(a,7f10.5)') '.........snow:',v_snow/a_globe
    call error_mesg('Lnd_stock_pe',message,NOTE)
    write (message,'(a,7f10.5)') '.........vegn:',v_vegn/a_globe
    call error_mesg('Lnd_stock_pe',message,NOTE)
    write (message,'(a,7f10.5)') '.......rivers:',river_value/a_globe
    call error_mesg('Lnd_stock_pe',message,NOTE)
endif

end subroutine Lnd_stock_pe

! ============================================================================
! initialize land model output on structured grid
! currently it is only cell_area, since for interpolation to work correctly it
! must have correct values over both land and ocean.
subroutine land_sg_diag_init(id_cellarea)
  integer, intent(out) :: id_cellarea

  ! local vars
  integer :: id_lon, id_lat, id_lonb, id_latb
  integer :: i
  logical :: used

  if(mpp_get_ntile_count(lnd%sg_domain)==1) then
     ! grid has just one tile, so we assume that the grid is regular lat-lon
     ! define longitude axes and its edges
     id_lonb = diag_axis_init ('lonb', lnd%coord_glonb, 'degrees_E', 'X', 'longitude edges', &
          set_name='land_sg', domain2=lnd%sg_domain )
     id_lon  = diag_axis_init ('lon',  lnd%coord_glon, 'degrees_E', 'X',   'longitude', &
          set_name='land_sg',  edges=id_lonb, domain2=lnd%sg_domain )

     ! define latitude axes and its edges
     id_latb = diag_axis_init ('latb', lnd%coord_glatb, 'degrees_N', 'Y', 'latitude edges',  &
          set_name='land_sg',  domain2=lnd%sg_domain   )
     id_lat = diag_axis_init ('lat',  lnd%coord_glat, 'degrees_N', 'Y', 'latitude', &
          set_name='land_sg', edges=id_latb, domain2=lnd%sg_domain)
  else
     id_lon = diag_axis_init ( 'grid_xt', [(real(i),i=1,size(lnd%coord_glon))], 'degrees_E', 'X', &
          'T-cell longitude', set_name='land_sg',  domain2=lnd%sg_domain)
     id_lat = diag_axis_init ( 'grid_yt', [(real(i),i=1,size(lnd%coord_glat))], 'degrees_N', 'Y', &
          'T-cell latitude', set_name='land_sg',  domain2=lnd%sg_domain)
  endif

  ! register land area on structured grid
  id_cellarea = register_static_field ( 'land_sg', 'cell_area', (/id_lon, id_lat/), &
       'total area in grid cell', 'm2', missing_value=-1.0 )
  call diag_field_add_attribute(id_cellarea,'cell_methods','area: sum')
  if ( id_cellarea > 0 ) used = send_data ( id_cellarea, lnd%sg_cellarea,     lnd%time )

end subroutine land_sg_diag_init

! ============================================================================
! initialize horizontal axes for land grid so that all sub-modules can use them,
! instead of creating their own
subroutine land_diag_init(clonb, clatb, clon, clat, time, &
     id_ug, id_band)
  real, intent(in) :: &
       clonb(:), clatb(:), & ! longitudes and latitudes of grid cells vertices,
                             ! specified for the global grid
       clon(:),  clat(:)     ! lon and lat of respective grid cell centers
  type(time_type), intent(in) :: time ! initial time for diagnostic fields
  integer, intent(out) :: &
       id_ug, id_band        ! IDs of respective diag. manager axes

  ! ---- local vars ----------------------------------------------------------
  integer :: nlon, nlat       ! sizes of respective axes
  integer :: axes(1)          ! array of horizontal axes for diag fields
  integer :: ug_dim_size      ! Size of the unstructured axis
  integer :: id_lon, id_lonb, id_lat, id_latb
  integer :: i
  integer, allocatable :: ug_dim_data(:) ! Unstructured axis data.
  character(128) :: long_name, flux_units
  character(32)  :: name

  ! Register the unstructured axis for the unstructured domain.
  call mpp_get_UG_compute_domain(lnd%ug_domain, size=ug_dim_size)
  if (.not. allocated(ug_dim_data)) then
      allocate(ug_dim_data(ug_dim_size))
  endif
  call mpp_get_UG_domain_grid_index(lnd%ug_domain, ug_dim_data)
  !--- grid_index needs to be starting from 0.
  ug_dim_data = ug_dim_data - 1
  id_ug = diag_axis_init("grid_index",  real(ug_dim_data), "none", "U", long_name="grid indices", &
                         set_name="land", DomainU=lnd%ug_domain, aux="geolon_t geolat_t")
  if (allocated(ug_dim_data)) then
      deallocate(ug_dim_data)
  endif

 ! Register horizontal axes that are required by the post-processing so that the output
 ! files can be "decompressed": converted from unstructured back to lon-lat or cubic sphere.
 ! The "grid_xt" and "grid_yt" axes should run from 1 to the total number of x- and
 ! y-points on cubic sphere face. It is assumed that all faces tiles contain the same
 ! number of x- and y-points.

  nlon = size(clon)
  nlat = size(clat)
  if(mpp_get_UG_domain_ntiles(lnd%ug_domain)==1) then
     ! grid has just one tile, so we assume that the grid is regular lat-lon
     ! define geographic axes and its edges
     id_lonb = diag_axis_init ('lonb', clonb, 'degrees_E', 'X', 'longitude edges', set_name='land')
     id_lon  = diag_axis_init ('lon',  clon,  'degrees_E', 'X', 'longitude', set_name='land',  edges=id_lonb)
     id_latb = diag_axis_init ('latb', clatb, 'degrees_N', 'Y', 'latitude edges', set_name='land')
     id_lat  = diag_axis_init ('lat',  clat,  'degrees_N', 'Y', 'latitude', set_name='land', edges=id_latb)
     ! add "compress" attribute to the unstructured grid axis
     call diag_axis_add_attribute(id_ug, "compress", "lat lon")
  else
     id_lon = diag_axis_init ( 'grid_xt', (/(real(i),i=1,nlon)/), 'degrees_E', 'X', &
          'T-cell longitude', set_name='land' )
     id_lat = diag_axis_init ( 'grid_yt', (/(real(i),i=1,nlat)/), 'degrees_N', 'Y', &
          'T-cell latitude', set_name='land' )
     ! add "compress" attribute to the unstructured grid axis
     call diag_axis_add_attribute(id_ug, "compress", "grid_yt grid_xt")
  endif

  id_band = diag_axis_init ('band',  (/1.0,2.0/), 'unitless', 'Z', 'spectral band', set_name='land' )

  ! set up an array of axes, for convenience
  axes = (/id_ug/)

  ! register auxilary coordinate variables
  id_geolon_t = register_static_field ( module_name, 'geolon_t', axes, &
       'longitude of grid cell centers', 'degrees_E', missing_value = -1.0e+20 )
  id_geolat_t = register_static_field ( module_name, 'geolat_t', axes, &
       'latitude of grid cell centers', 'degrees_N', missing_value = -1.0e+20 )

  ! register static diagnostic fields
  id_landfrac = register_static_field ( module_name, 'land_frac', axes, &
       'fraction of land in grid cell','unitless', missing_value=-1.0, area=id_cellarea)
  call diag_field_add_attribute(id_landfrac,'ocean_fillvalue',0.0)

  ! register areas and fractions for the rest of the diagnostic fields
  call register_tiled_area_fields(module_name, axes, time, id_area, id_frac)

  ! set the default filter (for area and subsampling) for consequent calls to
  ! register_tiled_diag_field
  call set_default_diag_filter('land')

  ! register regular (dynamic) diagnostic fields

  id_ntiles = register_tiled_diag_field(module_name,'ntiles', axes,  &
       time, 'number of tiles', 'unitless', missing_value=-1.0, op='sum', &
       cell_methods='area: mean')
  ! override cell methods to use poper method of interpolation from cubic
  ! sphere to lat-lon in post-processing


  id_VWS = register_tiled_diag_field ( module_name, 'VWS', axes, time, &
             'vapor storage on land', 'kg/m2', missing_value=-1.0e+20 )
  id_VWSc    = register_tiled_diag_field ( module_name, 'VWSc', axes, time, &
             'vapor mass in canopy air', 'kg/m2', missing_value=-1.0e+20 )
  id_LWS = register_tiled_diag_field ( module_name, 'LWS', axes, time, &
             'liquid storage on land', 'kg/m2', missing_value=-1.0e+20 )
  id_LWSv    = register_tiled_diag_field ( module_name, 'LWSv', axes, time, &
             'liquid interception storage', 'kg/m2', missing_value=-1.0e+20 )
  id_LWSs    = register_tiled_diag_field ( module_name, 'LWSs', axes, time, &
             'liquid storage in snowpack', 'kg/m2', missing_value=-1.0e+20 )
  id_LWSg    = register_tiled_diag_field ( module_name, 'LWSg', axes, time, &
             'liquid ground storage', 'kg/m2', missing_value=-1.0e+20 )
  id_FWS = register_tiled_diag_field ( module_name, 'FWS', axes, time, &
             'frozen storage on land', 'kg/m2', missing_value=-1.0e+20 )
  id_FWSv    = register_tiled_diag_field ( module_name, 'FWSv', axes, time, &
             'frozen interception storage', 'kg/m2', missing_value=-1.0e+20 )
  id_FWSs    = register_tiled_diag_field ( module_name, 'FWSs', axes, time, &
             'frozen storage in snowpack', 'kg/m2', missing_value=-1.0e+20 )
  id_FWSg    = register_tiled_diag_field ( module_name, 'FWSg', axes, time, &
             'frozen ground storage', 'kg/m2', missing_value=-1.0e+20 )
  id_HS = register_tiled_diag_field ( module_name, 'HS', axes, time, &
             'land heat storage', 'J/m2', missing_value=-1.0e+20 )
  id_HSv     = register_tiled_diag_field ( module_name, 'HSv', axes, time, &
             'interception heat storage', 'J/m2', missing_value=-1.0e+20 )
  id_HSs     = register_tiled_diag_field ( module_name, 'HSs', axes, time, &
             'heat storage in snowpack', 'J/m2', missing_value=-1.0e+20 )
  id_HSg     = register_tiled_diag_field ( module_name, 'HSg', axes, time, &
             'ground heat storage', 'J/m2', missing_value=-1.0e+20 )
  id_HSc     = register_tiled_diag_field ( module_name, 'HSc', axes, time, &
             'canopy-air heat storage', 'J/m2', missing_value=-1.0e+20 )

  ! register diag fields for river tracers
  ! FIXME: This FAILS because n_river_tracers is not set until after land_diag_init --BNS
  allocate(id_runf_tr(n_river_tracers))
  allocate(id_dis_tr(n_river_tracers))
  id_runf_tr(:) = -1
  id_dis_tr(:) = -1
  do i = 3, n_river_tracers
     call river_tracer_names(i,name=name,long_name=long_name,flux_units=flux_units)
     id_runf_tr(i) = register_tiled_diag_field ( module_name, 'lrunf_'//trim(name), axes, time, &
             'total rate of '//trim(long_name)//' runoff', trim(flux_units), missing_value=-1.0e+20 )
     id_dis_tr(i)  = register_diag_field ( module_name, 'dis_'//trim(name), axes, time, &
             trim(long_name)//' discharge to ocean', trim(flux_units), missing_value=-1.0e+20, &
             area = id_cellarea)
  enddo

  id_precip = register_tiled_diag_field ( module_name, 'precip', axes, time, &
             'precipitation rate', 'kg/(m2 s)', missing_value=-1.0e+20 )
  id_hprec = register_tiled_diag_field ( module_name, 'hprec', axes, time, &
             'sensible heat of precipitation', 'W/m2', missing_value=-1.0e+20 )
  id_lprec = register_tiled_diag_field ( module_name, 'lprec_l', axes, time, &
             'rainfall to land', 'kg/(m2 s)', missing_value=-1.0e+20 )
  id_lprecv = register_tiled_diag_field ( module_name, 'lprecv', axes, time, &
             'net rainfall to vegetation', 'kg/(m2 s)', missing_value=-1.0e+20)
  id_lprecs = register_tiled_diag_field ( module_name, 'lprecs', axes, time, &
             'rainfall to snow, minus drainage', 'kg/(m2 s)', missing_value=-1.0e+20)
  id_lprecg = register_tiled_diag_field ( module_name, 'lprecg', axes, time, &
             'effective rainfall to ground sfc', 'kg/(m2 s)', missing_value=-1.0e+20)
  id_hlprec  = register_tiled_diag_field ( module_name, 'hlprec', axes, time, &
             'total liq precipitation heat', 'W/m2', missing_value=-1.0e+20)
  id_hlprecv = register_tiled_diag_field ( module_name, 'hlprecv', axes, time, &
             'net liq heat to vegetation', 'W/m2', missing_value=-1.0e+20)
  id_hlprecs = register_tiled_diag_field ( module_name, 'hlprecs', axes, time, &
             'net liq heat to snow', 'W/m2', missing_value=-1.0e+20)
  id_hlprecg = register_tiled_diag_field ( module_name, 'hlprecg', axes, time, &
             'net liq heat to ground sfc', 'W/m2', missing_value=-1.0e+20)
  id_fprec = register_tiled_diag_field ( module_name, 'fprec_l', axes, time, &
             'snowfall to land', 'kg/(m2 s)', missing_value=-1.0e+20 )
  id_fprecv = register_tiled_diag_field ( module_name, 'fprecv', axes, time, &
             'net snowfall to vegetation', 'kg/(m2 s)', missing_value=-1.0e+20)
  id_fprecs = register_tiled_diag_field ( module_name, 'fprecs', axes, time, &
             'effective snowfall to snowpack', 'kg/(m2 s)', missing_value=-1.0e+20)
  id_hfprec = register_tiled_diag_field ( module_name, 'hfprec', axes, time, &
             'sens heat of snowfall', 'W/m2', missing_value=-1.0e+20)
  id_hfprecv = register_tiled_diag_field ( module_name, 'hfprecv', axes, time, &
             'net sol heat to vegetation', 'W/m2', missing_value=-1.0e+20)
  id_hfprecs = register_tiled_diag_field ( module_name, 'hfprecs', axes, time, &
             'net sol heat to snow', 'W/m2', missing_value=-1.0e+20)
  id_evap = register_tiled_diag_field ( module_name, 'evap', axes, time, &
             'vapor flux up from land', 'kg/(m2 s)', missing_value=-1.0e+20 )
  id_hevap = register_tiled_diag_field ( module_name, 'hevap', axes, time, &
             'sensible heat of evap', 'W/m2', missing_value=-1.0e+20 )
  id_levap   = register_tiled_diag_field ( module_name, 'levap', axes, time, &
             'vapor flux from all liq (inc Tr)', 'kg/(m2 s)', missing_value=-1.0e+20 )
  id_levapv = register_tiled_diag_field ( module_name, 'levapv', axes, time, &
             'vapor flux leaving intercepted liquid', 'kg/(m2 s)', missing_value=-1.0e+20)
  id_levaps = register_tiled_diag_field ( module_name, 'levaps', axes, time, &
             'vapor flux leaving snow liquid', 'kg/(m2 s)', missing_value=-1.0e+20)
  id_levapg = register_tiled_diag_field ( module_name, 'levapg', axes, time, &
             'vapor flux from ground liquid', 'kg/(m2 s)', missing_value=-1.0e+20)
  id_hlevap = register_tiled_diag_field ( module_name, 'hlevap', axes, time, &
             'vapor flux heat from liq source', 'W/m2', missing_value=-1.0e+20)
  id_hlevapv = register_tiled_diag_field ( module_name, 'hlevapv', axes, time, &
             'vapor heat from liq interc', 'W/m2', missing_value=-1.0e+20)
  id_hlevaps = register_tiled_diag_field ( module_name, 'hlevaps', axes, time, &
             'vapor heat from snow liq', 'W/m2', missing_value=-1.0e+20)
  id_hlevapg = register_tiled_diag_field ( module_name, 'hlevapg', axes, time, &
             'vapor heat from ground liq', 'W/m2', missing_value=-1.0e+20)
  id_fevap   = register_tiled_diag_field ( module_name, 'fevap', axes, time, &
             'vapor flux from all ice', 'kg/(m2 s)', missing_value=-1.0e+20 )
  id_fevapv = register_tiled_diag_field ( module_name, 'fevapv', axes, time, &
             'vapor flux leaving vegn ice', 'kg/(m2 s)', missing_value=-1.0e+20)
  id_fevaps = register_tiled_diag_field ( module_name, 'fevaps', axes, time, &
             'vapor flux leaving snow ice', 'kg/(m2 s)', missing_value=-1.0e+20)
  id_fevapg = register_tiled_diag_field ( module_name, 'fevapg', axes, time, &
             'vapor flux leaving ground ice', 'kg/(m2 s)', missing_value=-1.0e+20)
  id_hfevap = register_tiled_diag_field ( module_name, 'hfevap', axes, time, &
             'vapor flux heat from solid source', 'W/m2', missing_value=-1.0e+20)
  id_hfevapv = register_tiled_diag_field ( module_name, 'hfevapv', axes, time, &
             'vapor heat from sol interc', 'W/m2', missing_value=-1.0e+20)
  id_hfevaps = register_tiled_diag_field ( module_name, 'hfevaps', axes, time, &
             'vapor heat from snow sol', 'W/m2', missing_value=-1.0e+20)
  id_hfevapg = register_tiled_diag_field ( module_name, 'hfevapg', axes, time, &
             'vapor heat from ground sol', 'W/m2', missing_value=-1.0e+20)
  id_runf   = register_tiled_diag_field ( module_name, 'runf', axes, time, &
             'total runoff', 'kg/(m2 s)', missing_value=-1.0e+20 )
  id_hrunf   = register_tiled_diag_field ( module_name, 'hrunf', axes, time, &
             'sensible heat of total runoff', 'W/m2', missing_value=-1.0e+20 )
  id_lrunf   = register_tiled_diag_field ( module_name, 'lrunf', axes, time, &
             'total rate of liq runoff', 'kg/(m2 s)', missing_value=-1.0e+20 )
  id_lrunfs  = register_tiled_diag_field ( module_name, 'lrunfs', axes, time, &
             'rate of liq runoff via calving', 'kg/(m2 s)', missing_value=-1.0e+20 )
  id_lrunfg  = register_tiled_diag_field ( module_name, 'lrunfg', axes, time, &
             'rate of liq runoff, ground', 'kg/(m2 s)', missing_value=-1.0e+20 )
  id_hlrunf  = register_tiled_diag_field ( module_name, 'hlrunf', axes, time, &
             'heat of liq runoff', 'W/m2', missing_value=-1.0e+20 )
  id_hlrunfs  = register_tiled_diag_field ( module_name, 'hlrunfs', axes, time, &
             'heat of liq runoff from snow pack', 'W/m2', missing_value=-1.0e+20 )
  id_hlrunfg  = register_tiled_diag_field ( module_name, 'hlrunfg', axes, time, &
             'heat of liq surface runoff', 'W/m2', missing_value=-1.0e+20 )
  id_frunf   = register_tiled_diag_field ( module_name, 'frunf', axes, time, &
             'total rate of solid runoff', 'kg/(m2 s)', missing_value=-1.0e+20 )
  id_frunfs  = register_tiled_diag_field ( module_name, 'frunfs', axes, time, &
             'rate of solid calving', 'kg/(m2 s)', missing_value=-1.0e+20 )
  id_hfrunf  = register_tiled_diag_field ( module_name, 'hfrunf', axes, time, &
             'heat of total ice runoff', 'W/m2', missing_value=-1.0e+20 )
  id_hfrunfs  = register_tiled_diag_field ( module_name, 'hfrunfs', axes, time, &
             'heat of sol snow runoff', 'W/m2', missing_value=-1.0e+20 )
  id_melt    = register_tiled_diag_field ( module_name, 'melt', axes, time, &
             'total rate of melt', 'kg/(m2 s)', missing_value=-1.0e+20 )
  id_meltv   = register_tiled_diag_field ( module_name, 'meltv', axes, time, &
             'rate of melt, interception', 'kg/(m2 s)', missing_value=-1.0e+20 )
  id_melts   = register_tiled_diag_field ( module_name, 'melts', axes, time, &
             'rate of snow melt', 'kg/(m2 s)', missing_value=-1.0e+20 )
  id_meltg   = register_tiled_diag_field ( module_name, 'meltg', axes, time, &
             'rate of substrate thaw', 'kg/(m2 s)', missing_value=-1.0e+20 )
  id_fsw     = register_tiled_diag_field ( module_name, 'fsw', axes, time, &
             'net sw rad to land', 'W/m2', missing_value=-1.0e+20 )
  id_fswv    = register_tiled_diag_field ( module_name, 'fswv', axes, time, &
             'net sw rad to vegetation', 'W/m2', missing_value=-1.0e+20 )
  id_fsws    = register_tiled_diag_field ( module_name, 'fsws', axes, time, &
             'net sw rad to snow', 'W/m2', missing_value=-1.0e+20 )
  id_fswg    = register_tiled_diag_field ( module_name, 'fswg', axes, time, &
             'net sw rad to ground', 'W/m2', missing_value=-1.0e+20 )
  id_flw     = register_tiled_diag_field ( module_name, 'flw', axes, time, &
             'net lw rad to land', 'W/m2', missing_value=-1.0e+20 )
  id_flwv    = register_tiled_diag_field ( module_name, 'flwv', axes, time, &
             'net lw rad to vegetation', 'W/m2', missing_value=-1.0e+20 )
  id_flws    = register_tiled_diag_field ( module_name, 'flws', axes, time, &
             'net lw rad to snow', 'W/m2', missing_value=-1.0e+20 )
  id_flwg    = register_tiled_diag_field ( module_name, 'flwg', axes, time, &
             'net lw rad to ground', 'W/m2', missing_value=-1.0e+20 )
  id_cd_m    = register_tiled_diag_field ( module_name, 'cd_m', axes, time, &
       'drag coefficient for momentum', missing_value=-1e20)
  id_cd_t    = register_tiled_diag_field ( module_name, 'cd_t', axes, time, &
       'drag coefficient for heat and tracers', missing_value=-1e20)
  id_sens    = register_tiled_diag_field ( module_name, 'sens', axes, time, &
             'sens heat flux from land', 'W/m2', missing_value=-1.0e+20 )
  id_sensv   = register_tiled_diag_field ( module_name, 'sensv', axes, time, &
             'sens heat flux from vegn', 'W/m2', missing_value=-1.0e+20 )
  id_senss   = register_tiled_diag_field ( module_name, 'senss', axes, time, &
             'sens heat flux from snow', 'W/m2', missing_value=-1.0e+20 )
  id_sensg   = register_tiled_diag_field ( module_name, 'sensg', axes, time, &
             'sens heat flux from ground', 'W/m2', missing_value=-1.0e+20 )
  id_e_res_1 = register_tiled_diag_field ( module_name, 'e_res_1', axes, time, &
       'canopy air energy residual due to nonlinearities', 'W/m2', missing_value=-1e20)
  id_e_res_2 = register_tiled_diag_field ( module_name, 'e_res_2', axes, time, &
       'canopy energy residual due to nonlinearities', 'W/m2', missing_value=-1e20)
  id_z0m     = register_tiled_diag_field ( module_name, 'z0m', axes, time, &
             'momentum roughness of land', 'm', missing_value=-1.0e+20 )
  id_z0s     = register_tiled_diag_field ( module_name, 'z0s', axes, time, &
             'scalar roughness of land', 'm', missing_value=-1.0e+20 )
  id_con_g_h = register_tiled_diag_field ( module_name, 'con_g_h', axes, time, &
       'conductance for sensible heat between ground surface and canopy air', &
       'm/s', missing_value=-1.0 )
  id_con_g_v = register_tiled_diag_field ( module_name, 'con_g_v', axes, time, &
       'conductance for water vapor between ground surface and canopy air', &
       'm/s', missing_value=-1.0 )
  id_transp  = register_tiled_diag_field ( module_name, 'transp', axes, time, &
             'Transpiration', 'kg/(m2 s)', missing_value=-1.0e+20 )
  id_wroff = register_tiled_diag_field ( module_name, 'wroff', axes, time, &
             'rate of liquid runoff to rivers', 'kg/(m2 s)', missing_value=-1.0e+20 )
  id_sroff = register_tiled_diag_field ( module_name, 'sroff', axes, time, &
             'rate of solid runoff to rivers', 'kg/(m2 s)', missing_value=-1.0e+20 )
  id_htransp = register_tiled_diag_field ( module_name, 'htransp', axes, time, &
             'heat of transpired vapior', 'W/m2', missing_value=-1.0e+20 )
  id_huptake = register_tiled_diag_field ( module_name, 'huptk', axes, time, &
             'heat of soil water uptake', 'W/m2', missing_value=-1.0e+20 )
  id_hroff = register_tiled_diag_field ( module_name, 'hroff', axes, time, &
             'sensible heat of runoff', 'W/m2', missing_value=-1.0e+20 )
  id_gsnow   = register_tiled_diag_field ( module_name, 'gsnow', axes, time, &
             'sens heat into ground from snow', 'W/m2', missing_value=-1.0e+20 )
  call add_tiled_diag_field_alias ( id_gsnow, module_name, 'gflux', axes, time, &
             'obsolete, please use "gsnow" instead', 'W/m2', missing_value=-1.0e+20 )
  id_gequil   = register_tiled_diag_field ( module_name, 'gequil', axes, time, &
             'snow-subs equilibration flux', 'W/m2', missing_value=-1.0e+20 )
  id_grnd_flux = register_tiled_diag_field ( module_name, 'grnd_flux', axes, time, &
             'sensible heat into ground from surface', 'W/m2', missing_value=-1.0e+20 )
  id_levapg_max = register_tiled_diag_field ( module_name, 'Eg_max', axes, time, &
             'soil_water limit on vapor flux from ground liquid', 'kg/(m2 s)', missing_value=-1.0e+20)
  id_water = register_tiled_diag_field ( module_name, 'water', axes, time, &
             'column-integrated soil water', 'kg/m2', missing_value=-1.0e+20 )
  id_snow = register_tiled_diag_field ( module_name, 'snow', axes, time, &
             'column-integrated snow water', 'kg/m2', missing_value=-1.0e+20 )
  id_Trad    = register_tiled_diag_field ( module_name, 'Trad', axes, time, &
             'radiative sfc temperature', 'degK', missing_value=-1.0e+20 )
  id_Tca     = register_tiled_diag_field ( module_name, 'Tca', axes, time, &
             'canopy-air temperature', 'degK', missing_value=-1.0e+20 )
  id_qca     = register_tiled_diag_field ( module_name, 'qca', axes, time, &
             'canopy-air specific humidity', 'kg/kg', missing_value=-1.0 )
  id_cana_rh = register_tiled_diag_field ( module_name, 'rhca', axes, time, &
             'canopy-air relative humidity', 'kg/kg', missing_value=-1.0 )
  id_qco2    = register_tiled_diag_field ( module_name, 'qco2', axes, time, &
             'canopy-air CO2 moist mass mixing ratio', 'kg/kg', missing_value=-1.0 )
  id_qco2_dvmr = register_tiled_diag_field ( module_name, 'qco2_dvmr', axes, time, &
             'canopy-air CO2 dry volumetric mixing ratio', 'mol CO2/mol air', missing_value=-1.0 )
  id_fco2    = register_tiled_diag_field ( module_name, 'fco2', axes, time, &
             'flux of CO2 to canopy air', 'kg C/(m2 s)', missing_value=-1.0 )
  id_co2_mol_flux = register_tiled_diag_field ( module_name, 'co2_mol_flux', axes, time, &
             'flux of CO2 to the atmosphere', 'mol/(m2 s)', missing_value=-1.0 )
  id_swdn_dir = register_tiled_diag_field ( module_name, 'swdn_dir', (/id_ug,id_band/), time, &
       'downward direct short-wave radiation flux to the land surface', 'W/m2', missing_value=-999.0)
  id_swdn_dif = register_tiled_diag_field ( module_name, 'swdn_dif', (/id_ug,id_band/), time, &
       'downward diffuse short-wave radiation flux to the land surface', 'W/m2', missing_value=-999.0)
  id_swup_dir = register_tiled_diag_field ( module_name, 'swup_dir', (/id_ug,id_band/), time, &
       'direct short-wave radiation flux reflected by the land surface', 'W/m2', missing_value=-999.0)
  id_swup_dif = register_tiled_diag_field ( module_name, 'swup_dif', (/id_ug,id_band/), time, &
       'diffuse short-wave radiation flux reflected by the land surface', 'W/m2', missing_value=-999.0)
  id_lwdn = register_tiled_diag_field ( module_name, 'lwdn', axes, time, &
       'downward long-wave radiation flux to the land surface', 'W/m2', missing_value=-999.0)
  id_vegn_cover = register_tiled_diag_field ( module_name, 'vegn_cover', axes, time, &
             'fraction covered by vegetation', missing_value=-1.0 )
  id_vegn_cover_1 = register_tiled_diag_field ( module_name, 'vegn_cover:C', axes, time, &
             'fraction of vegetation in the top layer', missing_value=-1.0 )
  id_vegn_cover_U = register_tiled_diag_field ( module_name, 'vegn_cover:U', axes, time, &
             'fraction of vegetation in the understory', missing_value=-1.0 )
  id_cosz = register_tiled_diag_field ( module_name, 'coszen', axes, time, &
       'cosine of zenith angle', missing_value=-2.0 )
  id_albedo_dir = register_tiled_diag_field ( module_name, 'albedo_dir', &
       (/id_ug,id_band/), time, &
       'land surface albedo for direct light', missing_value=-1.0 )
  id_albedo_dif = register_tiled_diag_field ( module_name, 'albedo_dif', &
       (/id_ug,id_band/), time, &
       'land surface albedo for diffuse light', missing_value=-1.0 )
  id_vegn_refl_dir = register_tiled_diag_field(module_name, 'vegn_refl_dir', &
       (/id_ug, id_band/), time, &
       'black-background canopy reflectivity for direct light',missing_value=-1.0)
  id_vegn_refl_dif = register_tiled_diag_field(module_name, 'vegn_refl_dif', &
       (/id_ug, id_band/), time, &
       'black-background canopy reflectivity for diffuse light',missing_value=-1.0)
  id_vegn_refl_lw = register_tiled_diag_field ( module_name, 'vegn_refl_lw', axes, time, &
       'canopy reflectivity for thermal radiation', missing_value=-1.0)
  id_vegn_tran_dir = register_tiled_diag_field(module_name, 'vegn_tran_dir', &
       (/id_ug, id_band/), time, &
       'part of direct light that passes through canopy unscattered',missing_value=-1.0)
  id_vegn_tran_dif = register_tiled_diag_field(module_name, 'vegn_tran_dif', &
       (/id_ug, id_band/), time, &
       'black-background canopy transmittance for diffuse light',missing_value=-1.0)
  id_vegn_tran_lw = register_tiled_diag_field ( module_name, 'vegn_tran_lw', axes, time, &
       'canopy transmittance for thermal radiation', missing_value=-1.0)
  id_vegn_sctr_dir = register_tiled_diag_field(module_name, 'vegn_sctr_dir', &
       (/id_ug, id_band/), time, &
       'part of direct light scattered downward by canopy',missing_value=-1.0)
  id_subs_refl_dir = register_tiled_diag_field(module_name, 'subs_refl_dir', &
       (/id_ug, id_band/), time, &
       'substrate reflectivity for direct light',missing_value=-1.0)
  id_subs_refl_dif = register_tiled_diag_field(module_name, 'subs_refl_dif', &
       (/id_ug, id_band/), time, &
       'substrate reflectivity for diffuse light',missing_value=-1.0)
  id_subs_emis = register_tiled_diag_field(module_name, 'subs_emis', axes, time, &
       'substrate emissivity for long-wave radiation',missing_value=-1.0)
  id_grnd_T = register_tiled_diag_field ( module_name, 'Tgrnd', axes, time, &
       'ground surface temperature', 'degK', missing_value=-1.0 )
  id_grnd_rh = register_tiled_diag_field ( module_name, 'rh_grnd', axes, time, &
       'explicit ground relative humidity', missing_value=-1.0 )
  id_total_N = register_tiled_diag_field ( module_name, 'Ntot', axes, time, &
       'total land nitrogen', 'kg N/m2', missing_value=-1.0 )


  id_water_cons = register_tiled_diag_field ( module_name, 'water_cons', axes, time, &
       'water non-conservation in update_land_model_fast_0d', 'kg/(m2 s)', missing_value=-1.0 )
  id_carbon_cons = register_tiled_diag_field ( module_name, 'carbon_cons', axes, time, &
       'carbon non-conservation in update_land_model_fast_0d', 'kgC/(m2 s)', missing_value=-1.0 )
  id_nitrogen_cons = register_tiled_diag_field ( module_name, 'nitrogen_cons', axes, time, &
       'nitrogen non-conservation in update_land_model_fast_0d', 'kgN/(m2 s)', missing_value=-1.0 )

  id_snow_frac = register_tiled_diag_field ( module_name, 'snow_frac', axes, time, &
             'fraction of area that is covered by snow','1', missing_value=-1.0e+20)

  ! CMOR/CMIP variables
  id_pcp = register_tiled_diag_field ( cmor_name, 'pcp', axes, time, &
             'Total Precipitation', 'kg m-2 s-1', missing_value=-1.0e+20, &
             standard_name='total_precipitation_flux', fill_missing=.TRUE.)
  id_prra = register_tiled_diag_field ( cmor_name, 'prra', axes, time, &
             'Rainfall Rate', 'kg m-2 s-1', missing_value=-1.0e+20, &
             standard_name='rainfall_flux', fill_missing=.TRUE.)
  id_prveg = register_tiled_diag_field ( cmor_name, 'prveg', axes, time, &
             'Precipitation onto Canopy', 'kg m-2 s-1', missing_value=-1.0e+20, &
             standard_name='precipitation_flux_onto_canopy', fill_missing=.TRUE.)
  id_evspsblveg = register_tiled_diag_field ( cmor_name, 'evspsblveg', axes, time, &
             'Evaporation from Canopy', 'kg m-2 s-1', missing_value=-1.0e+20, &
             standard_name='water_evaporation_flux_from_canopy', fill_missing=.TRUE.)
  id_evspsblsoi = register_tiled_diag_field ( cmor_name, 'evspsblsoi', axes, time, &
             'Water Evaporation from Soil', 'kg m-2 s-1', missing_value=-1.0e+20, &
             standard_name='water_evaporation_flux_from_soil', fill_missing=.TRUE.)
  id_ec = register_tiled_diag_field ( cmor_name, 'ec', axes, time, &
             'Interception evaporation', 'kg m-2 s-1', missing_value=-1.0e+20, &
             standard_name='water_evaporation_flux_from_canopy', fill_missing=.TRUE.)
  id_eow = register_tiled_diag_field ( cmor_name, 'eow', axes, time, &
             'Open Water Evaporation', 'kg m-2 s-1', missing_value=-1.0e+20, &
             standard_name='surface_water_evaporation_flux', fill_missing=.TRUE.)
  id_esn = register_tiled_diag_field ( cmor_name, 'esn', axes, time, &
             'Snow Evaporation', 'kg m-2 s-1', missing_value=-1.0e+20, &
             standard_name='water_evaporation_flux', fill_missing=.TRUE.)
  id_et = register_tiled_diag_field ( cmor_name, 'et', axes, time, &
             'Total Evapotranspiration', 'kg m-2 s-1', missing_value=-1.0e+20, &
             standard_name='surface_evapotranspiration', fill_missing=.TRUE.)

  id_snw = register_tiled_diag_field ( cmor_name, 'snw', axes, time, &
             'Surface Snow Amount','kg m-2', standard_name='surface_snow_amount', &
             missing_value=-1.0e+20, fill_missing=.TRUE.)
  id_snd = register_tiled_diag_field ( cmor_name, 'snd', axes, time, &
             'Snow Depth','m', standard_name='surface_snow_thickness', &
             missing_value=-1.0e+20, fill_missing=.TRUE.)
  id_snc = register_cmor_fraction_field ('snc', 'Snow Area Fraction', axes, &
             standard_name='surface_snow_area_fraction')
  call diag_field_add_attribute (id_snc, 'normalization_corrected', 'per-grid-cell-area')
  id_lwsnl = register_tiled_diag_field ( cmor_name, 'lwsnl', axes, time, &
             'Liquid Water Content of Snow Layer','kg m-2', &
             standard_name='liquid_water_content_of_surface_snow', &
             missing_value=-1.0e+20, fill_missing=.TRUE.)
  id_snm = register_tiled_diag_field ( cmor_name, 'snm', axes, time, &
             'Surface Snow Melt','kg m-2 s-1', standard_name='surface_snow_melt_flux', &
             missing_value=-1.0e+20, fill_missing=.TRUE.)
  id_hfdsn = register_tiled_diag_field ( cmor_name, 'hfdsn', axes, time, &
             'Downward Heat Flux into Snow Where Land over Land','W m-2', standard_name='surface_downward_heat_flux_in_snow', &
             missing_value=-1.0e+20, fill_missing=.TRUE.)

  id_tws = register_diag_field(cmor_name, 'mrtws', axes, time, &
             'Terrestrial Water Storage','kg m-2', &
             standard_name='land_water_amount', &
             area=get_area_id('land'))
  call diag_field_add_attribute(id_tws,'cell_methods','area: mean')
  call diag_field_add_attribute(id_tws,'ocean_fillvalue',0.0)

  id_hflsLut = register_tiled_diag_field ( cmor_name, 'hflsLut', axes, time, &
             'Latent Heat Flux on Land Use Tile', 'W m-2', missing_value=-1.0e+20, &
             standard_name='surface_upward_latent_heat_flux', fill_missing=.FALSE.)
  id_rlusLut = register_tiled_diag_field ( cmor_name, 'rlusLut', axes, time, &
             'Surface Upwelling Longwave on Land Use Tile', 'W m-2', &
             missing_value=-1.0e+20, standard_name='surface_upwelling_longwave_flux_in_air')
  id_rsusLut = register_tiled_diag_field ( cmor_name, 'rsusLut', axes, time, &
             'Surface Upwelling Shortwave on Land Use Tile', 'W m-2', &
             missing_value=-1.0e+20, standard_name='surface_upwelling_shortwave_flux_in_air')
  id_sweLut = register_tiled_diag_field ( cmor_name, 'sweLut', axes, time, &
             'Snow Water Equivalent on Land Use Tile','m', standard_name='lwe_thickness_of_surface_snow_amount', &
             missing_value=-1.0e+20, fill_missing=.FALSE. )
  id_tslsiLut = register_tiled_diag_field ( cmor_name, 'tslsiLut', axes, time, &
             'Surface Skin Temperature on Land Use Tile', 'K', &
             missing_value=-1.0, standard_name='surface_temperature')

  call add_tiled_diag_field_alias(id_transp, cmor_name, 'tran', axes, time, &
      'Transpiration', 'kg m-2 s-1', missing_value=-1.0e+20, &
      standard_name='transpiration_flux', fill_missing=.TRUE.)
  call add_tiled_diag_field_alias(id_sens, cmor_name, 'hfssLut', axes, time, &
      'Sensible Heat Flux on Land Use Tile', 'W m-2', missing_value=-1.0e+20, &
      standard_name='surface_upward_sensible_heat_flux')
  call add_tiled_diag_field_alias(id_fevaps, cmor_name, 'sbl', axes, time, &
      'Surface Snow and Ice Sublimation Flux', 'kg m-2 s-1', missing_value=-1.0e+20, &
      standard_name='tendency_of_atmosphere_mass_content_of_water_vapor_due_to_sublimation_of_surface_snow_and_ice', &
      fill_missing=.TRUE.)

  id_sftlf = register_static_field ( cmor_name, 'sftlf', axes, &
             'Land Area Fraction','%', standard_name='land_area_fraction', &
             area=id_cellarea)
  call diag_field_add_attribute(id_sftlf,'cell_methods','area: mean')
  call diag_field_add_attribute(id_sftlf,'ocean_fillvalue',0.0)

  id_sftgif = register_static_field ( cmor_name, 'sftgif', axes, &
             'Fraction of Grid Cell Covered with Glacier','%', standard_name='land_ice_area_fraction', &
             area=id_cellarea)
  call diag_field_add_attribute(id_sftgif,'cell_methods','area: mean')
  call diag_field_add_attribute(id_sftgif,'ocean_fillvalue',0.0)

  id_residualFrac = register_static_field ( cmor_name, 'residualFrac', axes, &
             'Fraction of Grid Cell that is Land but Neither Vegetation-Covered nor Bare Soil','%', &
             standard_name='area_fraction', &
             area=id_cellarea)
  call diag_field_add_attribute(id_residualFrac,'cell_methods','area: mean')
  call diag_field_add_attribute(id_residualFrac,'ocean_fillvalue',0.0)

  id_cropFrac    = register_cmor_fraction_field ('cropFrac', 'Crop Fraction', axes)
  id_cropFracC3  = register_cmor_fraction_field ('cropFracC3', 'C3 Crop Fraction', axes)
  id_cropFracC4  = register_cmor_fraction_field ('cropFracC4', 'C4 Crop Fraction', axes)

  id_pastureFrac = register_cmor_fraction_field ('pastureFrac', 'Anthropogenic Pasture Fraction', axes)
  id_vegFrac     = register_cmor_fraction_field ('vegFrac',  'Total vegetated fraction', axes)
  id_treeFrac    = register_cmor_fraction_field ('treeFrac', 'Tree Cover Fraction', axes)

  id_grassFrac   = register_cmor_fraction_field ('grassFrac', 'Natural Grass Fraction', axes)
  id_grassFracC3 = register_cmor_fraction_field ('grassFracC3', 'C3 Natural Grass Fraction', axes)
  id_grassFracC4 = register_cmor_fraction_field ('grassFracC4', 'C4 Natural Grass Fraction', axes)
  id_c3pftFrac   = register_cmor_fraction_field ('c3PftFrac', 'Total C3 PFT Cover Fraction', axes)
  id_c4pftFrac   = register_cmor_fraction_field ('c4PftFrac', 'Total C4 PFT Cover Fraction', axes)
  ! Add 'normalization_corrected' attribute to variable that are reported by send_cellfrac_cohort_data.
  ! This attribute works as an indicator that the normalization mistake was fixed in the code
  ! and the variable does not need re-normalization in refineDiag scripts
  call diag_field_add_attribute (id_cropFracC3,  'normalization_corrected', 'per-grid-cell-area')
  call diag_field_add_attribute (id_cropFracC4,  'normalization_corrected', 'per-grid-cell-area')
  call diag_field_add_attribute (id_vegFrac,     'normalization_corrected', 'per-grid-cell-area')
  call diag_field_add_attribute (id_treeFrac,    'normalization_corrected', 'per-grid-cell-area')
  call diag_field_add_attribute (id_c3pftFrac,   'normalization_corrected', 'per-grid-cell-area')
  call diag_field_add_attribute (id_c4pftFrac,   'normalization_corrected', 'per-grid-cell-area')
  call diag_field_add_attribute (id_grassFrac,   'normalization_corrected', 'per-grid-cell-area')
  call diag_field_add_attribute (id_grassFracC3, 'normalization_corrected', 'per-grid-cell-area')
  call diag_field_add_attribute (id_grassFracC4, 'normalization_corrected', 'per-grid-cell-area')

  ! LUMIP land fractions
  id_fracLut_psl = register_diag_field ( cmor_name, 'fracLut_psl', axes, time, &
             'Fraction of Grid Cell for Each Land Use Tile','%', &
             standard_name='area_fraction', area=id_cellarea)
  id_fracLut_crp = register_diag_field ( cmor_name, 'fracLut_crop', axes, time, &
             'Fraction of Grid Cell for Each Land Use Tile','%', &
             standard_name='area_fraction', area=id_cellarea)
  id_fracLut_pst = register_diag_field ( cmor_name, 'fracLut_pst', axes, time, &
             'Fraction of Grid Cell for Each Land Use Tile','%', &
             standard_name='area_fraction', area=id_cellarea)
  id_fracLut_urb = register_diag_field ( cmor_name, 'fracLut_urbn', axes, time, &
             'Fraction of Grid Cell for Each Land Use Tile','%', &
             standard_name='area_fraction', area=id_cellarea)
  call diag_field_add_attribute(id_fracLut_psl,'cell_methods','area: mean')
  call diag_field_add_attribute(id_fracLut_crp,'cell_methods','area: mean')
  call diag_field_add_attribute(id_fracLut_pst,'cell_methods','area: mean')
  call diag_field_add_attribute(id_fracLut_urb,'cell_methods','area: mean')

  id_cLand = register_tiled_diag_field ( cmor_name, 'cLand', axes, time, &
             'Total Carbon in All Terrestrial Carbon Pools', 'kg m-2', &
             standard_name='mass_content_of_carbon_in_vegetation_and_litter_and_soil_and_forestry_and_agricultural_products', &
             missing_value=-1.0, fill_missing=.TRUE. )
  ! add alias for compatibility with older diag tables
  call add_tiled_diag_field_alias(id_cLand, module_name, 'Ctot', axes, time, &
     'total land carbon', 'kg C/m2', missing_value=-1.0)
  id_cTot1 = register_tiled_diag_field ( module_name, 'cTot1', axes, time, &
             'Total Carbon in All Terrestrial Carbon Pools, Except Canopy Air', 'kg m-2', &
             missing_value=-1.0, fill_missing=.TRUE. )

  id_nbp = register_tiled_diag_field ( cmor_name, 'nbp', axes, time, &
             'Carbon Mass Flux out of Atmosphere due to Net Biospheric Production on Land', &
             'kg m-2 s-1', missing_value=-1.0, &
             standard_name='surface_net_downward_mass_flux_of_carbon_dioxide_expressed_as_carbon_due_to_all_land_processes', &
             fill_missing=.TRUE.)
  call add_tiled_diag_field_alias ( id_nbp, cmor_name, 'netAtmosLandCO2Flux', axes, time, &
             'Net flux of CO2 between atmosphere and land (positive into land) as a result of all processes.', &
             'kg m-2 s-1', missing_value=-1.0, &
             standard_name='surface_net_downward_mass_flux_of_carbon_dioxide_expressed_as_carbon_due_to_all_land_processes', &
             fill_missing=.TRUE. )
  call add_tiled_diag_field_alias ( id_nbp, cmor_name, 'necbLut', axes, time, &
             'net rate of C accumulation (or loss) on land use tile', &
             'kg m-2 s-1', missing_value=-1.0, &
             standard_name='surface_net_downward_mass_flux_of_carbon_dioxide_expressed_as_carbon_due_to_all_land_processes', &
             fill_missing=.TRUE. )

  id_nLand = register_tiled_diag_field ( cmor_name, 'nLand', axes, time, &
             'Total nitrogen in all terrestrial nitrogen pools', 'kg m-2', &
             standard_name='mass_content_of_nitrogen_in_vegetation_and_litter_and_soil_and_forestry_and_agricultural_products', &
             missing_value=-1.0, fill_missing=.TRUE. )
end subroutine land_diag_init


! ==============================================================================
function register_cmor_fraction_field(field_name, long_name, axes, standard_name) result(id); integer :: id
  character(*), intent(in) :: field_name, long_name
  integer, intent(in) :: axes(:)
  character(*), intent(in), optional :: standard_name

  if (present(standard_name)) then
     id = register_diag_field ( cmor_name, field_name, axes, lnd%time, &
           long_name,'%', standard_name=standard_name, area=id_cellarea)
  else
     id = register_diag_field ( cmor_name, field_name, axes, lnd%time, &
           long_name,'%', standard_name='area_fraction', area=id_cellarea)
  endif
  call diag_field_add_attribute(id,'cell_methods','area: mean')
  call diag_field_add_attribute(id,'ocean_fillvalue',0.0)
end function register_cmor_fraction_field


! ==============================================================================
subroutine send_cellfrac_data(id, f, scale)
  integer, intent(in) :: id ! id of the diagnostic field
  procedure(tile_test_func) :: f ! existence detector function
  real, intent(in), optional  :: scale ! scaling factor, for unit conversions

  ! ---- local vars
  integer :: l
  logical :: used
  type(land_tile_type), pointer :: tile
  type(land_tile_enum_type) :: ce
  real :: frac(lnd%ls:lnd%le)
  real :: scale_

  if (.not.id>0) return ! do nothing if the field was not registered
  scale_ = 100.0 ! by fractions are in percent
  if (present(scale)) scale_ = scale

  frac(:) = 0.0
  ce = first_elmt(land_tile_map, ls=lnd%ls)
  do while (loop_over_tiles(ce, tile, l))
     ! accumulate fractions
     if (f(tile)) then
        frac(l) = frac(l)+tile%frac*scale_*lnd%ug_landfrac(l)
     endif
  enddo
  used = send_data(id, frac, lnd%time)
end subroutine send_cellfrac_data

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function any_tile(tile) result(answer); logical :: answer
  type(land_tile_type), pointer :: tile

  answer = .TRUE.
end function any_tile

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function is_crop(tile) result(answer); logical :: answer
  type(land_tile_type), pointer :: tile

  answer = .FALSE.
  if (.not.associated(tile)) return
  if (.not.associated(tile%vegn)) return
  answer = (tile%vegn%landuse == LU_CROP)
end function is_crop

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function is_pasture(tile) result(answer); logical :: answer
  type(land_tile_type), pointer :: tile

  answer = .FALSE.
  if (.not.associated(tile)) return
  if (.not.associated(tile%vegn)) return
  answer = (tile%vegn%landuse == LU_PAST)
end function is_pasture

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function is_urban(tile) result(answer); logical :: answer
  type(land_tile_type), pointer :: tile

  answer = .FALSE.
  if (.not.associated(tile)) return
  if (.not.associated(tile%vegn)) return
  answer = (tile%vegn%landuse == LU_URBN)
end function is_urban

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! primary or secondary lamd, for LUMIP
function is_psl(tile) result(answer); logical :: answer
  type(land_tile_type), pointer :: tile

  answer = .FALSE.
  if (.not.associated(tile)) return
  if (.not.associated(tile%vegn)) return
  answer = (tile%vegn%landuse == LU_NTRL.or.tile%vegn%landuse == LU_SCND)
end function is_psl

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! pasture or rangeland, for LUMIP
function is_pst(tile) result(answer); logical :: answer
  type(land_tile_type), pointer :: tile

  answer = .FALSE.
  if (.not.associated(tile)) return
  if (.not.associated(tile%vegn)) return
  answer = (tile%vegn%landuse == LU_PAST.or.tile%vegn%landuse == LU_RANGE)
end function is_pst

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function is_residual(tile) result(answer); logical :: answer
  type(land_tile_type), pointer :: tile

  answer = .FALSE.
  if (.not.associated(tile)) return
  answer = (associated(tile%glac).or.associated(tile%lake))
end function is_residual

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function is_glacier(tile) result(answer); logical :: answer
  type(land_tile_type), pointer :: tile

  answer = .FALSE.
  if (.not.associated(tile)) return
  answer = associated(tile%glac)
end function is_glacier

! ==============================================================================
! given two test functions -- for tiles and for cohorts -- sends to the diagnostics
! the fraction of area covered by the vegetation that matches both criteria.
! For example, tile test might return TRUE for natural vegetation only, and
! cohort test -- for grass, giving fraction of natural grass as a result.
! NOTE that non-vegetated tiles are always skipped by this function.
subroutine send_cellfrac_cohort_data(id, ttest, ctest, scale)
  integer, intent(in) :: id ! id of the diagnostic field
  procedure(tile_test_func)   :: ttest ! returns TRUE for tiles whose fraction is counted
  procedure(cohort_test_func) :: ctest ! returns TRUE for cohorts whose fraction is counted
  real, intent(in), optional  :: scale ! scaling factor, for unit conversions

  ! ---- local vars
  integer :: l
  logical :: used
  type(land_tile_type), pointer :: tile
  type(land_tile_enum_type) :: ce
  real :: frac(lnd%ls:lnd%le)
  real :: scale_

  if (.not.id>0) return ! do nothing if the field was not registered
  scale_ = 100.0 ! by fractions are in percent
  if (present(scale)) scale_ = scale

  frac(:) = 0.0
  ce = first_elmt(land_tile_map, ls=lnd%ls)
  do while (loop_over_tiles(ce, tile, l))
     if (associated(tile%vegn).and.ttest(tile)) then
        ! calculate fraction of area covered by suitable cohorts
        frac(l) = frac(l) + tile%frac*scale_*cohort_area_frac(tile%vegn,ctest)*lnd%ug_landfrac(l)
     endif
  enddo
  used = send_data(id, frac, lnd%time)
end subroutine send_cellfrac_cohort_data


! ============================================================================
! allocates boundary data for land domain and current number of tiles
subroutine realloc_land2cplr ( bnd )
  type(land_data_type), intent(inout) :: bnd     ! data to allocate

  ! ---- local vars
  integer :: n_tiles

  call dealloc_land2cplr(bnd, dealloc_discharges=.FALSE.)

  bnd%domain    = lnd%sg_domain
  bnd%ug_domain = lnd%ug_domain
  n_tiles = max_n_tiles()

  ! allocate data according to the domain boundaries
  allocate( bnd%mask(lnd%ls:lnd%le,n_tiles) )

  allocate( bnd%tile_size(lnd%ls:lnd%le,n_tiles) )
  allocate( bnd%t_surf(lnd%ls:lnd%le,n_tiles) )
  allocate( bnd%t_ca(lnd%ls:lnd%le,n_tiles) )
  allocate( bnd%tr(lnd%ls:lnd%le,n_tiles,ntcana) )
  allocate( bnd%albedo(lnd%ls:lnd%le,n_tiles) )
  allocate( bnd%albedo_vis_dir(lnd%ls:lnd%le,n_tiles) )
  allocate( bnd%albedo_nir_dir(lnd%ls:lnd%le,n_tiles) )
  allocate( bnd%albedo_vis_dif(lnd%ls:lnd%le,n_tiles) )
  allocate( bnd%albedo_nir_dif(lnd%ls:lnd%le,n_tiles) )
  allocate( bnd%rough_mom(lnd%ls:lnd%le,n_tiles) )
  allocate( bnd%rough_heat(lnd%ls:lnd%le,n_tiles) )
  allocate( bnd%rough_scale(lnd%ls:lnd%le,n_tiles) )

  bnd%mask              = .FALSE.
  bnd%tile_size         = init_value
  bnd%t_surf            = init_value
  bnd%t_ca              = init_value
  bnd%tr                = init_value
  bnd%albedo            = init_value
  bnd%albedo_vis_dir    = init_value
  bnd%albedo_nir_dir    = init_value
  bnd%albedo_vis_dif    = init_value
  bnd%albedo_nir_dif    = init_value
  bnd%rough_mom         = init_value
  bnd%rough_heat        = init_value
  bnd%rough_scale       = init_value

  ! in contrast to the rest of the land boundary condition fields, discharges
  ! are specified per grid cell, not per tile; therefore they should not be
  ! re-allocated when the number of tiles changes. In fact, they must not be
  ! changed at all here because their values are assigned in update_land_model_fast,
  ! not in update_land_bc_*, and therefore would be lost if re-allocated.
  if (.not.associated(bnd%discharge)) then
     allocate( bnd%discharge          (lnd%is:lnd%ie,lnd%js:lnd%je) )
     allocate( bnd%discharge_heat     (lnd%is:lnd%ie,lnd%js:lnd%je) )
     allocate( bnd%discharge_snow     (lnd%is:lnd%ie,lnd%js:lnd%je) )
     allocate( bnd%discharge_snow_heat(lnd%is:lnd%ie,lnd%js:lnd%je) )

     ! discharge and discharge_snow must be, in contrast to the rest of the boundary
     ! values, filled with zeroes. The reason is because not all of the usable elements
     ! are updated by the land model (only coastal points are).
     bnd%discharge           = 0.0
     bnd%discharge_heat      = 0.0
     bnd%discharge_snow      = 0.0
     bnd%discharge_snow_heat = 0.0
  endif
end subroutine realloc_land2cplr


! ============================================================================
! deallocates boundary data memory
! NOTE that the discharges should be deallocated only at the final clean-up
! stage; during the model run they should be preserved unchanged even when
! other fields are reallocated.

#define __DEALLOC__(x) if (associated(x)) deallocate(x)

subroutine dealloc_land2cplr ( bnd, dealloc_discharges )
  type(land_data_type), intent(inout) :: bnd  ! data to de-allocate
  logical, intent(in) :: dealloc_discharges

  __DEALLOC__( bnd%tile_size )
  __DEALLOC__( bnd%tile_size )
  __DEALLOC__( bnd%t_surf )
  __DEALLOC__( bnd%t_ca )
  __DEALLOC__( bnd%tr )
  __DEALLOC__( bnd%albedo )
  __DEALLOC__( bnd%albedo_vis_dir )
  __DEALLOC__( bnd%albedo_nir_dir )
  __DEALLOC__( bnd%albedo_vis_dif )
  __DEALLOC__( bnd%albedo_nir_dif )
  __DEALLOC__( bnd%rough_mom )
  __DEALLOC__( bnd%rough_heat )
  __DEALLOC__( bnd%rough_scale )
  __DEALLOC__( bnd%mask )

  if (dealloc_discharges) then
     __DEALLOC__( bnd%discharge           )
     __DEALLOC__( bnd%discharge_heat      )
     __DEALLOC__( bnd%discharge_snow      )
     __DEALLOC__( bnd%discharge_snow_heat )
  end if

end subroutine dealloc_land2cplr


! ============================================================================
! allocates boundary data for land domain and current number of tiles;
! initializes data for data override.
! NOTE: previously the body of the procedure was in the flux_exchange_init,
! currently it is called from land_model_init
subroutine realloc_cplr2land( bnd )
  type(atmos_land_boundary_type), intent(inout) :: bnd

  ! ---- local vars
  integer :: kd

  call dealloc_cplr2land(bnd)

  ! allocate data according to the domain boundaries
  kd = max_n_tiles()

  allocate( bnd%t_flux(lnd%ls:lnd%le,kd) )
  allocate( bnd%lw_flux(lnd%ls:lnd%le,kd) )
  allocate( bnd%sw_flux(lnd%ls:lnd%le,kd) )
  allocate( bnd%lprec(lnd%ls:lnd%le,kd) )
  allocate( bnd%fprec(lnd%ls:lnd%le,kd) )
  allocate( bnd%tprec(lnd%ls:lnd%le,kd) )
  allocate( bnd%dhdt(lnd%ls:lnd%le,kd) )
  allocate( bnd%dhdq(lnd%ls:lnd%le,kd) )
  allocate( bnd%drdt(lnd%ls:lnd%le,kd) )
  allocate( bnd%p_surf(lnd%ls:lnd%le,kd) )
  allocate( bnd%tr_flux(lnd%ls:lnd%le,kd,ntcana) )
  allocate( bnd%dfdtr(lnd%ls:lnd%le,kd,ntcana) )

  allocate( bnd%lwdn_flux(lnd%ls:lnd%le,kd) )
  allocate( bnd%swdn_flux(lnd%ls:lnd%le,kd) )
  allocate( bnd%sw_flux_down_vis_dir(lnd%ls:lnd%le,kd) )
  allocate( bnd%sw_flux_down_total_dir(lnd%ls:lnd%le,kd) )
  allocate( bnd%sw_flux_down_vis_dif(lnd%ls:lnd%le,kd) )
  allocate( bnd%sw_flux_down_total_dif(lnd%ls:lnd%le,kd) )
  allocate( bnd%cd_t(lnd%ls:lnd%le,kd) )
  allocate( bnd%cd_m(lnd%ls:lnd%le,kd) )
  allocate( bnd%bstar(lnd%ls:lnd%le,kd) )
  allocate( bnd%ustar(lnd%ls:lnd%le,kd) )
  allocate( bnd%wind(lnd%ls:lnd%le,kd) )
  allocate( bnd%z_bot(lnd%ls:lnd%le,kd) )

  allocate( bnd%drag_q(lnd%ls:lnd%le,kd) )

  bnd%t_flux                 = init_value
  bnd%lw_flux                = init_value
  bnd%sw_flux                = init_value
  bnd%lprec                  = init_value
  bnd%fprec                  = init_value
  bnd%tprec                  = init_value
  bnd%dhdt                   = init_value
  bnd%dhdq                   = init_value
  bnd%drdt                   = init_value
  bnd%p_surf                 = init_value
  bnd%tr_flux                = init_value
  bnd%dfdtr                  = init_value

  bnd%lwdn_flux              = init_value
  bnd%swdn_flux              = init_value
  bnd%sw_flux_down_vis_dir   = init_value
  bnd%sw_flux_down_total_dir = init_value
  bnd%sw_flux_down_vis_dif   = init_value
  bnd%sw_flux_down_total_dif = init_value
  bnd%cd_t                   = init_value
  bnd%cd_m                   = init_value
  bnd%bstar                  = init_value
  bnd%ustar                  = init_value
  bnd%wind                   = init_value
  bnd%z_bot                  = init_value

  bnd%drag_q                 = init_value

end subroutine realloc_cplr2land


! ============================================================================
subroutine dealloc_cplr2land( bnd )
  type(atmos_land_boundary_type), intent(inout) :: bnd

  __DEALLOC__( bnd%t_flux )
  __DEALLOC__( bnd%lw_flux )
  __DEALLOC__( bnd%sw_flux )
  __DEALLOC__( bnd%lprec )
  __DEALLOC__( bnd%fprec )
  __DEALLOC__( bnd%tprec )
  __DEALLOC__( bnd%dhdt )
  __DEALLOC__( bnd%dhdq )
  __DEALLOC__( bnd%drdt )
  __DEALLOC__( bnd%p_surf )
  __DEALLOC__( bnd%lwdn_flux )
  __DEALLOC__( bnd%swdn_flux )
  __DEALLOC__( bnd%sw_flux_down_vis_dir )
  __DEALLOC__( bnd%sw_flux_down_total_dir )
  __DEALLOC__( bnd%sw_flux_down_vis_dif )
  __DEALLOC__( bnd%sw_flux_down_total_dif )
  __DEALLOC__( bnd%cd_t )
  __DEALLOC__( bnd%cd_m )
  __DEALLOC__( bnd%bstar )
  __DEALLOC__( bnd%ustar )
  __DEALLOC__( bnd%wind )
  __DEALLOC__( bnd%z_bot )
  __DEALLOC__( bnd%drag_q )
  __DEALLOC__( bnd%tr_flux )
  __DEALLOC__( bnd%dfdtr )

end subroutine dealloc_cplr2land


! ===========================================================================
!  Prints checksums of the various fields in the atmos_land_boundary_type.
subroutine atm_lnd_bnd_type_chksum(id, timestep, albt)
    character(len=*), intent(in) :: id  ! Label to differentiate where this
                      ! routine is being called from.
    integer         , intent(in) :: timestep ! An integer to indicate which
                      ! timestep this routine is being called for.
    type(atmos_land_boundary_type), intent(in) :: albt
    integer ::   n, outunit

    outunit = stdout()

    write(outunit,*) 'BEGIN CHECKSUM(atmos_land_boundary_type):: ', id, timestep
    write(outunit,100) 'albt%t_flux                ', mpp_chksum( albt%t_flux)
    write(outunit,100) 'albt%lw_flux               ', mpp_chksum( albt%lw_flux)
    write(outunit,100) 'albt%lwdn_flux             ', mpp_chksum( albt%lwdn_flux)
    write(outunit,100) 'albt%sw_flux               ', mpp_chksum( albt%sw_flux)
    write(outunit,100) 'albt%swdn_flux               ', mpp_chksum( albt%swdn_flux)
    write(outunit,100) 'albt%lprec                 ', mpp_chksum( albt%lprec)
    write(outunit,100) 'albt%fprec                 ', mpp_chksum( albt%fprec)
    write(outunit,100) 'albt%tprec                 ', mpp_chksum( albt%tprec)
    write(outunit,100) 'albt%sw_flux_down_vis_dir  ', mpp_chksum( albt%sw_flux_down_vis_dir)
    write(outunit,100) 'albt%sw_flux_down_total_dir', mpp_chksum( albt%sw_flux_down_total_dir)
    write(outunit,100) 'albt%sw_flux_down_vis_dif  ', mpp_chksum( albt%sw_flux_down_vis_dif)
    write(outunit,100) 'albt%sw_flux_down_total_dif', mpp_chksum( albt%sw_flux_down_total_dif)
    write(outunit,100) 'albt%dhdt                  ', mpp_chksum( albt%dhdt)
    write(outunit,100) 'albt%dhdq                  ', mpp_chksum( albt%dhdq)
    write(outunit,100) 'albt%drdt                  ', mpp_chksum( albt%drdt)
    write(outunit,100) 'albt%cd_m                  ', mpp_chksum( albt%cd_m)
    write(outunit,100) 'albt%cd_t                  ', mpp_chksum( albt%cd_t)
    write(outunit,100) 'albt%ustar                 ', mpp_chksum( albt%ustar)
    write(outunit,100) 'albt%bstar                 ', mpp_chksum( albt%bstar)
    write(outunit,100) 'albt%wind                  ', mpp_chksum( albt%wind)
    write(outunit,100) 'albt%z_bot                 ', mpp_chksum( albt%z_bot)
    write(outunit,100) 'albt%drag_q                ', mpp_chksum( albt%drag_q)
    write(outunit,100) 'albt%p_surf                ', mpp_chksum( albt%p_surf)
    do n = 1,size(albt%tr_flux,3)
    write(outunit,100) 'albt%tr_flux               ', mpp_chksum( albt%tr_flux(:,:,n))
    enddo
    do n = 1,size(albt%dfdtr,3)
    write(outunit,100) 'albt%dfdtr                 ', mpp_chksum( albt%dfdtr(:,:,n))
    enddo

100 FORMAT("CHECKSUM::",A32," = ",Z20)

end subroutine atm_lnd_bnd_type_chksum


! ===========================================================================
!  Prints checksums of the various fields in the land_data_type.
subroutine land_data_type_chksum(id, timestep, land)
    character(len=*), intent(in) :: id ! Label to differentiate where this
        ! routine in being called from.
    integer         , intent(in) :: timestep ! An integer to indicate which
        ! timestep this routine is being called for.
    type(land_data_type), intent(in) :: land
    integer ::   n, outunit

    outunit = stdout()

    write(outunit,*) 'BEGIN CHECKSUM(land_data_type):: ', id, timestep
    write(outunit,100) 'land%tile_size         ',mpp_chksum(land%tile_size)
    write(outunit,100) 'land%t_surf            ',mpp_chksum(land%t_surf)
    write(outunit,100) 'land%t_ca              ',mpp_chksum(land%t_ca)
    write(outunit,100) 'land%albedo            ',mpp_chksum(land%albedo)
    write(outunit,100) 'land%albedo_vis_dir    ',mpp_chksum(land%albedo_vis_dir)
    write(outunit,100) 'land%albedo_nir_dir    ',mpp_chksum(land%albedo_nir_dir)
    write(outunit,100) 'land%albedo_vis_dif    ',mpp_chksum(land%albedo_vis_dif)
    write(outunit,100) 'land%albedo_nir_dif    ',mpp_chksum(land%albedo_nir_dif)
    write(outunit,100) 'land%rough_mom         ',mpp_chksum(land%rough_mom)
    write(outunit,100) 'land%rough_heat        ',mpp_chksum(land%rough_heat)
    write(outunit,100) 'land%rough_scale       ',mpp_chksum(land%rough_scale)

    do n = 1, size(land%tr,3)
    write(outunit,100) 'land%tr                ',mpp_chksum(land%tr(:,:,n))
    enddo
    write(outunit,100) 'land%discharge         ',mpp_chksum(land%discharge)
    write(outunit,100) 'land%discharge_snow    ',mpp_chksum(land%discharge_snow)
    write(outunit,100) 'land%discharge_heat    ',mpp_chksum(land%discharge_heat)


100 FORMAT("CHECKSUM::",A32," = ",Z20)
end subroutine land_data_type_chksum


! the code below defines the accessor routines that are used to access fields of the
! tile data structure in collective operations, like restart i/o. Fore example, a statement
! DEFINE_LAND_ACCESSOR_0D(real,lwup)
! defines a function "land_lwup_ptr" that, given a tile, returns a pointer to the field
! called "lwup" in this tile. The procedure implementing a collective operation would
! enumerate all the tiles within the domain, call accessor routine for each of them, and
! get or set the value pointed to by the accessor routine.
#define DEFINE_LAND_ACCESSOR_0D(xtype,x) subroutine land_ ## x ## _ptr(t,p);\
type(land_tile_type),pointer::t;xtype,pointer::p;p=>NULL();if(associated(t))p=>t%x;end subroutine

DEFINE_LAND_ACCESSOR_0D(real,frac)
DEFINE_LAND_ACCESSOR_0D(real,lwup)
DEFINE_LAND_ACCESSOR_0D(real,e_res_1)
DEFINE_LAND_ACCESSOR_0D(real,e_res_2)

! ============================================================================
! tile existence detector: returns a logical value indicating wether component
! model tile exists or not
logical function land_tile_exists(tile)
  type(land_tile_type), pointer :: tile
  land_tile_exists = associated(tile)
end function land_tile_exists

#define DEFINE_TAG_ACCESSOR(x) subroutine  x ## _tag_ptr(t,p);\
type(land_tile_type),pointer::t;integer,pointer::p;p=>NULL();if(associated(t))\
then;if (associated(t%x)) p=>t%x%tag;endif;end subroutine

DEFINE_TAG_ACCESSOR(glac)
DEFINE_TAG_ACCESSOR(lake)
DEFINE_TAG_ACCESSOR(soil)
DEFINE_TAG_ACCESSOR(vegn)

end module land_model_mod
