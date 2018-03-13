module vegetation_mod

#include "../shared/debug.inc"

#ifdef INTERNAL_FILE_NML
use mpp_mod, only: input_nml_file
#else
use fms_mod, only: open_namelist_file
#endif

use fms_mod, only: error_mesg, NOTE, WARNING, FATAL, file_exist, &
                   close_file, check_nml_error, stdlog, string
use mpp_mod, only: mpp_sum, mpp_max, mpp_pe, mpp_root_pe

use time_manager_mod, only: time_type, time_type_to_real, get_date, day_of_year, operator(-)
use field_manager_mod, only: fm_field_name_len
use constants_mod,    only: tfreeze, rdgas, rvgas, hlv, hlf, cp_air, PI
use sphum_mod, only: qscomp

use vegn_tile_mod, only: vegn_tile_type, &
     vegn_seed_demand, vegn_seed_supply, vegn_seed_N_supply, vegn_add_bliving, &
     vegn_relayer_cohorts_ppa, vegn_mergecohorts_ppa, &
     vegn_tile_LAI, vegn_tile_SAI, &
     cpw, clw, csw
use soil_tile_mod, only: soil_tile_type, num_l, dz, &
     soil_ave_temp, soil_ave_theta0, soil_ave_theta1, soil_psi_stress, &
     N_LITTER_POOLS, LEAF, l_shortname, l_longname
use land_constants_mod, only : NBANDS, BAND_VIS, d608, mol_C, mol_CO2, &
     seconds_per_year
use land_tile_mod, only : land_tile_map, land_tile_type, land_tile_enum_type, &
     first_elmt, loop_over_tiles
use land_tile_diag_mod, only : OP_SUM, OP_AVERAGE, OP_MAX, cmor_name, &
     register_tiled_static_field, register_tiled_diag_field, &
     send_tile_data, diag_buff_type, register_cohort_diag_field, send_cohort_data, &
     set_default_diag_filter, add_tiled_diag_field_alias
use land_utils_mod, only : check_conservation_1, check_conservation_2
use land_data_mod, only : lnd, log_version
use land_io_mod, only : read_field
use land_tile_io_mod, only: land_restart_type, &
     init_land_restart, open_land_restart, save_land_restart, free_land_restart, &
     add_restart_axis, add_tile_data, add_int_tile_data, add_scalar_data, &
     get_scalar_data, get_tile_data, get_int_tile_data, field_exists, &
     add_text_data, get_text_data
use vegn_data_mod, only : read_vegn_data_namelist, FORM_WOODY, FORM_GRASS, PT_C3, PT_C4, &
     LEAF_ON, LU_NTRL, nspecies, C2B, N_HARV_POOLS, HARV_POOL_NAMES, &
     HARV_POOL_CLEARED, HARV_POOL_WOOD_FAST, HARV_POOL_WOOD_MED, HARV_POOL_WOOD_SLOW, &
     HARV_POOL_PAST, HARV_POOL_CROP, &
     spdata, mcv_min, mcv_lai, agf_bs, tau_drip_l, tau_drip_s, T_transp_min, &
     do_ppa, cold_month_threshold, soil_carbon_depth_scale, &
     fsc_pool_spending_time, ssc_pool_spending_time, harvest_spending_time, &
     N_HARV_POOLS, HARV_POOL_NAMES, HARV_POOL_PAST, HARV_POOL_CROP, HARV_POOL_CLEARED, &
     HARV_POOL_WOOD_FAST, HARV_POOL_WOOD_MED, HARV_POOL_WOOD_SLOW, &
     c2n_N_fixer,C2N_SEED, N_limits_live_biomass
use vegn_cohort_mod, only : vegn_cohort_type, &
     init_cohort_allometry_ppa, init_cohort_hydraulics, &
     update_species, update_bio_living_fraction, get_vegn_wet_frac, &
     vegn_data_cover, btotal, height_from_biomass, leaf_area_from_biomass, &
     cohort_root_properties
use canopy_air_mod, only : cana_turbulence
use soil_mod, only : soil_data_beta, get_soil_litter_C, redistribute_peat_carbon, &
     register_litter_soilc_diag_fields

use cohort_io_mod, only :  read_create_cohorts, create_cohort_dimension, &
     add_cohort_data, add_int_cohort_data, get_cohort_data, get_int_cohort_data
use land_debug_mod, only : is_watch_point, set_current_point, check_temp_range, &
     check_var_range
use vegn_radiation_mod, only : vegn_radiation_init, vegn_radiation
use vegn_photosynthesis_mod, only : vegn_photosynthesis_init, vegn_photosynthesis, &
     co2_for_photosynthesis, vegn_phot_co2_option, VEGN_PHOT_CO2_INTERACTIVE
use static_vegn_mod, only : read_static_vegn_namelist, static_vegn_init, static_vegn_end, &
     read_static_vegn
use vegn_dynamics_mod, only : vegn_dynamics_init, vegn_dynamics_end, &
     vegn_carbon_int_lm3, vegn_carbon_int_ppa,    &
     vegn_phenology_lm3,  vegn_phenology_ppa,     &
     vegn_growth, vegn_starvation_ppa, vegn_biogeography, &
     vegn_reproduction_ppa
use vegn_disturbance_mod, only : vegn_disturbance_init, vegn_nat_mortality_lm3, &
     vegn_disturbance, update_fuel
use vegn_harvesting_mod, only : &
     vegn_harvesting_init, vegn_harvesting_end, vegn_harvesting
use vegn_fire_mod, only : vegn_fire_init, vegn_fire_end, update_fire_data, fire_option, FIRE_LM3
use soil_carbon_mod, only : soil_carbon_option, SOILC_CORPSE, SOILC_CORPSE_N, &
     N_C_TYPES, C_FAST, C_SLOW, c_shortname, c_longname, &
     add_litter, soil_NH4_deposition, soil_NO3_deposition, soil_org_N_deposition, &
     cull_cohorts
use vegn_util_mod, only: kill_small_cohorts_ppa

implicit none
private

! ==== public interfaces =====================================================
public :: read_vegn_namelist
public :: vegn_init
public :: vegn_end
public :: save_vegn_restart

public :: vegn_radiation
public :: vegn_diffusion

public :: vegn_step_1
public :: vegn_step_2
public :: vegn_step_3
public :: update_derived_vegn_data  ! given state variables, calculate derived values

public :: update_vegn_slow

public :: cohort_area_frac
public :: cohort_test_func
public :: any_vegn, is_tree, is_grass, is_c3, is_c4, is_c3grass, is_c4grass
! ==== end of public interfaces ==============================================

! ==== module constants ======================================================
character(len=*), parameter :: module_name = 'vegn'
#include "../shared/version_variable.inc"

! size of cohort initial condition array
integer, parameter :: MAX_INIT_COHORTS = 10

abstract interface
  ! the following interface describes the "detector function", which is passed
  ! through the argument list and must return TRUE for any cohort that meets
  ! specific condition (e.g. is tree), FALSE otherwise
  logical function cohort_test_func(cc)
     import vegn_cohort_type
     type(vegn_cohort_type), intent(in) :: cc
  end function cohort_test_func
end interface

! ==== module variables ======================================================

!---- namelist ---------------------------------------------------------------
logical :: lm2               = .false.
real    :: init_Wl           = 0
real    :: init_Ws           = 0
real    :: init_Tv           = 288.0
integer :: init_n_cohorts    = 1
character(fm_field_name_len) :: init_cohort_species(MAX_INIT_COHORTS) = 'tempdec'
real    :: init_cohort_nindivs(MAX_INIT_COHORTS) = 1.0  ! initial individual density, individual/m2
real    :: init_cohort_bl(MAX_INIT_COHORTS)      = 0.05 ! initial biomass of leaves, kg C/individual
real    :: init_cohort_blv(MAX_INIT_COHORTS)     = 0.0  ! initial biomass of labile store, kg C/individual
real    :: init_cohort_br(MAX_INIT_COHORTS)      = 0.05 ! initial biomass of fine roots, kg C/individual
real    :: init_cohort_bsw(MAX_INIT_COHORTS)     = 0.05 ! initial biomass of sapwood, kg C/individual
real    :: init_cohort_bwood(MAX_INIT_COHORTS)   = 0.05 ! initial biomass of heartwood, kg C/individual
real    :: init_cohort_bseed(MAX_INIT_COHORTS)   = 0.05 ! initial biomass of seeds, kg C/individual
real    :: init_cohort_myc_scav(MAX_INIT_COHORTS) = 0.0 ! initial scavenger mycorrhizal biomass, kgC/m2
real    :: init_cohort_myc_mine(MAX_INIT_COHORTS) = 0.0 ! initial miner mycorrhizal biomass, kgC/m2
real    :: init_cohort_n_fixer(MAX_INIT_COHORTS)  = 0.0 ! initial N fixer microbe biomass, kgC/m2
real    :: init_cohort_stored_N_mult(MAX_INIT_COHORTS) = 1.5 ! Multiple of initial leaf N + root N
real    :: init_cohort_age(MAX_INIT_COHORTS)     = 0.0  ! initial cohort age, year
real    :: init_cohort_height(MAX_INIT_COHORTS)  = 0.1  ! initial cohort height, m
real    :: init_cohort_nsc_frac(MAX_INIT_COHORTS)= 3.0  ! initial cohort NSC, as fraction of max. bl
real    :: init_cohort_nsn_frac(MAX_INIT_COHORTS)= 3.0  ! initial cohort NSN, as fraction of max. bl
character(32) :: rad_to_use = 'big-leaf' ! or 'two-stream'
character(32) :: snow_rad_to_use = 'ignore' ! or 'paint-leaves'

logical :: allow_external_gaps = .TRUE. ! if TRUE, there may be gaps between
   ! cohorts of the canopy layers; otherwise canopies are stretched to fill
   ! every layer completely. These gaps are called "external" in contrast to the
   ! "internal" gaps that are created by branch drop processes within cohort
   ! canopies
logical :: do_cohort_dynamics   = .TRUE. ! if true, do vegetation growth
logical :: do_patch_disturbance = .TRUE. !
logical :: do_phenology         = .TRUE.
real :: tau_smooth_theta_phen = 30.0  ! days; time scale of low-band-pass-filter of soil moisture
                                      ! in the drought-deciduous phenology
logical :: xwilt_available      = .TRUE.
logical :: do_biogeography      = .TRUE.
logical :: do_seed_transport    = .TRUE.
real    :: min_Wl=-1.0, min_Ws=-1.0 ! threshold values for condensation numerics, kg/m2:
   ! if water or snow on canopy fall below these values, the derivatives of
   ! condensation are set to zero, thereby prohibiting switching from condensation to
   ! evaporation in one time step.
real    :: tau_smooth_ncm = 0.0 ! Time scale for ncm smoothing (low-pass
   ! filtering), years. 0 retrieves previous behavior (no smoothing)
real :: rav_lit_0         = 0.0 ! constant litter resistance to vapor
real :: rav_lit_vi        = 0.0 ! litter resistance to vapor per LAI+SAI
real :: rav_lit_fsc       = 0.0 ! litter resistance to vapor per fsc
real :: rav_lit_ssc       = 0.0 ! litter resistance to vapor per ssc
real :: rav_lit_deadmic   = 0.0 ! litter resistance to vapor per dead microbe C
real :: rav_lit_bwood     = 0.0 ! litter resistance to vapor per bwood

logical :: do_peat_redistribution = .FALSE.

namelist /vegn_nml/ &
    lm2, init_Wl, init_Ws, init_Tv, cpw, clw, csw, &
    init_n_cohorts, init_cohort_species, init_cohort_nindivs, &
    init_cohort_bl, init_cohort_blv, init_cohort_br, init_cohort_bsw, &
    init_cohort_bwood, init_cohort_bseed, &
    init_cohort_height, &
    init_cohort_myc_scav, init_cohort_myc_mine, init_cohort_n_fixer, &
    init_cohort_stored_N_mult, &
    rad_to_use, snow_rad_to_use, &
    allow_external_gaps, &
    do_cohort_dynamics, do_patch_disturbance, do_phenology, tau_smooth_theta_phen, &
    xwilt_available, &
    do_biogeography, do_seed_transport, &
    min_Wl, min_Ws, tau_smooth_ncm, &
    rav_lit_0, rav_lit_vi, rav_lit_fsc, rav_lit_ssc, rav_lit_deadmic, rav_lit_bwood, &
    do_peat_redistribution

!---- end of namelist --------------------------------------------------------

logical :: module_is_initialized =.FALSE.
real    :: delta_time      ! fast time step
real    :: dt_fast_yr      ! fast time step in years
real    :: steps_per_day   ! number of fast time steps per day
real    :: weight_av_phen  ! weight for low-band-pass soil moisture smoother, for drought-deciduous phenology

! diagnostic field ids
integer :: id_vegn_type, id_height, id_height_ave, &
   id_temp, id_wl, id_ws, &
   id_lai, id_lai_var, id_lai_std, id_sai, id_leafarea, id_leaf_size, id_laii, &
   id_root_density, id_root_zeta, id_rs_min, id_leaf_refl, id_leaf_tran, &
   id_leaf_emis, id_snow_crit, id_stomatal, &
   id_an_op, id_an_cl,&
   id_bl, id_blv, id_br, id_bsw, id_bwood, id_bseed, id_btot, id_nsc, id_bl_max, id_br_max, id_bsw_max, &
   id_leaf_N,id_root_N,id_wood_N,id_sapwood_N,id_seed_N,id_stored_N,id_Ntot,id_Ngain,id_Nloss,&
   id_mrz_scav_C,id_mrz_mine_C,id_Nfix_C,&
   id_mrz_scav_N,id_mrz_mine_N,id_Nfix_N,&
   id_species, id_status, &
   id_con_v_h, id_con_v_v, id_fuel, id_harv_pool_C(N_HARV_POOLS), id_harv_pool_N(N_HARV_POOLS), &
   id_harv_rate_C(N_HARV_POOLS), id_tot_harv_pool_C, id_tot_harv_rate_C, id_tot_harv_pool_N, &
   id_csmoke_pool, id_nsmoke_pool, id_csmoke_rate, id_fsc_in, id_fsc_out, id_ssc_in, &
   id_ssc_out, id_deadmic_in, id_deadmic_out, id_veg_in, id_veg_out, &
   id_tile_nitrogen_gain, id_tile_nitrogen_loss, &
   id_fsc_pool_ag, id_fsc_rate_ag, id_fsc_pool_bg, id_fsc_rate_bg,&
   id_ssc_pool_ag, id_ssc_rate_ag, id_ssc_pool_bg, id_ssc_rate_bg,&
   id_t_ann, id_t_cold, id_p_ann, id_ncm, &
   id_lambda, id_afire, id_atfall, id_closs, id_cgain, id_wdgain, id_leaf_age, &
   id_phot_co2, id_theph, id_psiph, id_evap_demand, &
   id_ncohorts, &
   id_nindivs,  &
   id_nlayers, id_dbh, id_dbh_max, &
   id_crownarea, &
   id_soil_water_supply, id_gdd, id_tc_pheno, id_zstar_1, &
   id_psi_r, id_psi_l, id_psi_x, id_Kxi, id_Kli, id_w_scale, id_RHi, &
   id_brsw, id_topyear, id_growth_prev_day, &
   id_lai_kok, id_DanDlai, id_PAR_dn, id_PAR_net
integer, dimension(N_LITTER_POOLS, N_C_TYPES) :: &
   id_litter_buff_C, id_litter_buff_N, &
   id_litter_rate_C, id_litter_rate_N
! CMOR/CMIP variables
integer :: id_lai_cmor, id_cVeg, id_cLeaf, id_cWood, id_cRoot, id_cMisc, id_cProduct, id_cAnt, &
   id_fFire, id_fGrazing, id_fHarvest, id_fLuc, id_fAnthDisturb, id_fProductDecomp, id_cw, &
   id_nVeg, id_nLeaf, id_nRoot, id_nStem, id_nOther
! ==== end of module variables ===============================================

contains

! ============================================================================
subroutine read_vegn_namelist()
  ! ---- local vars
  integer :: unit         ! unit for namelist i/o
  integer :: io           ! i/o status for the namelist
  integer :: ierr         ! error code, returned by i/o routines
  logical :: use_static_veg ! if true, switch off vegetation dynamics

  call read_vegn_data_namelist()
  call read_static_vegn_namelist(use_static_veg)

  call log_version(version, module_name, &
  __FILE__)
#ifdef INTERNAL_FILE_NML
    read (input_nml_file, nml=vegn_nml, iostat=io)
    ierr = check_nml_error(io, 'vegn_nml')
#else
  if (file_exist('input.nml')) then
     unit = open_namelist_file()
     ierr = 1;
     do while (ierr /= 0)
        read (unit, nml=vegn_nml, iostat=io, end=10)
        ierr = check_nml_error (io, 'vegn_nml')
     enddo
10   continue
     call close_file (unit)
  endif
#endif

  unit=stdlog()

  ! switch off vegetation dynamics if static vegetation is set
  if (use_static_veg) then
     call error_mesg('vegn_init', &
          'use_static_veg=.TRUE., switching off vegetation dynamics', NOTE)
     write(unit,*)'use_static_veg=.TRUE., switching off vegetation dynamics'
     do_cohort_dynamics   = .FALSE.
     do_patch_disturbance = .FALSE.
     do_phenology         = .FALSE.
     do_biogeography      = .FALSE.
     do_seed_transport    = .FALSE.
  endif

  if (mpp_pe() == mpp_root_pe()) then
     write(unit, nml=vegn_nml)
  endif

  ! ---- initialize vegetation radiation options
  call vegn_radiation_init(rad_to_use, snow_rad_to_use)

  ! ---- initialize vegetation photosynthesis options
  call vegn_photosynthesis_init()

end subroutine read_vegn_namelist


! ============================================================================
! initialize vegetation
subroutine vegn_init ( id_ug, id_band, id_cellarea )
  integer, intent(in) :: id_ug   !<Unstructured axis id.
  integer, intent(in) :: id_band ! ID of spectral band axis
  integer, intent(in) :: id_cellarea ! ID of cell area diag field, for cell measures

  ! ---- local vars
  integer :: unit         ! unit for various i/o
  type(land_tile_enum_type)     :: ce    ! current tile list element
  type(land_tile_type), pointer :: tile  ! pointer to current tile
  integer :: n_accum
  integer :: nmn_acm
  type(land_restart_type) :: restart1, restart2
  logical :: restart_1_exists, restart_2_exists
  real, allocatable :: t_ann(:),t_cold(:),p_ann(:),ncm(:) ! buffers for biodata reading
  logical :: did_read_biodata
  integer :: i,j,l,n ! indices of current tile
  integer :: init_cohort_spp(MAX_INIT_COHORTS)

  module_is_initialized = .TRUE.

  ! ---- make module copy of time and calculate time step ------------------
  delta_time = time_type_to_real(lnd%dt_fast)
  steps_per_day = 86400.0/delta_time
  dt_fast_yr = delta_time/seconds_per_year

  ! --- initialize smoothing parameters for phenology; see http://en.wikipedia.org/wiki/Low-pass_filter
  weight_av_phen = delta_time/(delta_time+tau_smooth_theta_phen*86400.0)

  ! ---- initialize vegn state ---------------------------------------------
  n_accum = 0
  nmn_acm = 0
  call open_land_restart(restart1,'INPUT/vegn1.res.nc',restart_1_exists)
  call open_land_restart(restart2,'INPUT/vegn2.res.nc',restart_2_exists)
  if (restart_1_exists) then
     call error_mesg('vegn_init',&
          'reading NetCDF restarts "INPUT/vegn1.res.nc" and "INPUT/vegn2.res.nc"',&
          NOTE)

     ! read the cohort index and generate appropriate number of cohorts
     ! for each vegetation tile
     call read_create_cohorts(restart1)
     call get_cohort_data(restart1, 'tv', cohort_tv_ptr)
     call get_cohort_data(restart1, 'wl', cohort_wl_ptr)
     call get_cohort_data(restart1, 'ws', cohort_ws_ptr)

     ! read global variables
     call get_scalar_data(restart2,'n_accum',n_accum)
     call get_scalar_data(restart2,'nmn_acm',nmn_acm)

     ! read cohort data
     call get_int_cohort_data(restart2, 'species', cohort_species_ptr)
     call get_cohort_data(restart2, 'hite', cohort_height_ptr)
     call get_cohort_data(restart2, 'bl', cohort_bl_ptr)
     call get_cohort_data(restart2, 'blv', cohort_blv_ptr)
     call get_cohort_data(restart2, 'br', cohort_br_ptr)
     call get_cohort_data(restart2, 'bsw', cohort_bsw_ptr)
     call get_cohort_data(restart2, 'bwood', cohort_bwood_ptr)
     if(field_exists(restart2,'nsc')) then
        ! nsc is used as a flag to distinguish between PPA and LM3 restarts
        call get_cohort_data(restart2,'nsc',cohort_nsc_ptr)
        call get_cohort_data(restart2,'bseed',cohort_bseed_ptr)
        call get_cohort_data(restart2,'bl_max',cohort_bl_max_ptr)
        call get_cohort_data(restart2,'br_max',cohort_br_max_ptr)
        ! isa 201707
        if (field_exists(restart2,'bsw_max')) &
           call get_cohort_data(restart2, 'bsw_max', cohort_bsw_max_ptr)
        call get_cohort_data(restart2,'dbh',cohort_dbh_ptr)
        call get_cohort_data(restart2,'crownarea',cohort_crownarea_ptr)
        call get_cohort_data(restart2,'nindivs',cohort_nindivs_ptr)
        call get_int_cohort_data(restart2,'layer',cohort_layer_ptr)
        call get_int_cohort_data(restart2,'firstlayer',cohort_firstlayer_ptr)
        call get_cohort_data(restart2,'cohort_age',cohort_age_ptr)
        ! TODO: possibly initialize cohort age with tile age if the cohort_age is
        !       not present in the restart
        call get_cohort_data(restart2,'BM_ys',cohort_BM_ys_ptr)
        call get_cohort_data(restart2,'DBH_ys',cohort_DBH_ys_ptr)
        call get_cohort_data(restart2,'topyear',cohort_topyear_ptr)
        call get_cohort_data(restart2,'gdd',cohort_gdd_ptr)

        call get_tile_data(restart2,'drop_wl',vegn_drop_wl_ptr)
        call get_tile_data(restart2,'drop_ws',vegn_drop_ws_ptr)
        call get_tile_data(restart2,'drop_hl',vegn_drop_hl_ptr)
        call get_tile_data(restart2,'drop_hs',vegn_drop_hs_ptr)

        ! read hydraulics-related variables
        call get_cohort_data(restart2, 'psi_r', cohort_psi_r_ptr )
        call get_cohort_data(restart2, 'psi_x', cohort_psi_x_ptr )
        call get_cohort_data(restart2, 'psi_l', cohort_psi_l_ptr )
        call get_cohort_data(restart2, 'Kxa',   cohort_Kxa_ptr )
        call get_cohort_data(restart2, 'Kla',   cohort_Kla_ptr )

        call get_cohort_data(restart2, 'growth_prev_day', cohort_growth_previous_day_ptr )
        call get_cohort_data(restart2, 'growth_prev_day_tmp', cohort_growth_previous_day_tmp_ptr )
        call get_cohort_data(restart2, 'brsw', cohort_brsw_ptr)
     endif

     if (field_exists(restart2,'scav_myc_C_reservoir')) then
        call get_cohort_data(restart2, 'scav_myc_C_reservoir', cohort_scav_myc_C_reservoir_ptr)
        call get_cohort_data(restart2, 'scav_myc_N_reservoir', cohort_scav_myc_N_reservoir_ptr)
        call get_cohort_data(restart2, 'mine_myc_C_reservoir', cohort_mine_myc_C_reservoir_ptr)
        call get_cohort_data(restart2, 'mine_myc_N_reservoir', cohort_mine_myc_N_reservoir_ptr)
        call get_cohort_data(restart2, 'N_fixer_C_reservoir', cohort_N_fixer_C_reservoir_ptr)
        call get_cohort_data(restart2, 'N_fixer_N_reservoir', cohort_N_fixer_N_reservoir_ptr)
     endif

     call get_cohort_data(restart2, 'bliving', cohort_bliving_ptr)
     call get_int_cohort_data(restart2, 'status', cohort_status_ptr)
     if(field_exists(restart2,'leaf_age')) &
          call get_cohort_data(restart2,'leaf_age',cohort_leaf_age_ptr)
     call get_cohort_data(restart2, 'npp_prev_day', cohort_npp_previous_day_ptr )

     if (field_exists(restart2,'myc_scavenger_biomass_C')) then
        call get_cohort_data(restart2, 'myc_scavenger_biomass_C', cohort_myc_scavenger_biomass_C_ptr )
        call get_cohort_data(restart2, 'myc_scavenger_biomass_N', cohort_myc_scavenger_biomass_N_ptr )
        call get_cohort_data(restart2, 'myc_miner_biomass_C', cohort_myc_miner_biomass_C_ptr )
        call get_cohort_data(restart2, 'myc_miner_biomass_N', cohort_myc_miner_biomass_N_ptr )
        call get_cohort_data(restart2, 'N_fixer_biomass_C', cohort_N_fixer_biomass_C_ptr )
        call get_cohort_data(restart2, 'N_fixer_biomass_N', cohort_N_fixer_biomass_N_ptr )
     endif
     if (field_exists(restart2,'cohort_leaf_N')) then
        call get_cohort_data(restart2, 'cohort_stored_N', cohort_stored_N_ptr )
        call get_cohort_data(restart2, 'cohort_leaf_N', cohort_leaf_N_ptr )
        call get_cohort_data(restart2, 'cohort_seed_N', cohort_seed_N_ptr )
        call get_cohort_data(restart2, 'cohort_wood_N', cohort_wood_N_ptr )
        call get_cohort_data(restart2, 'cohort_sapwood_N', cohort_sapwood_N_ptr )
        call get_cohort_data(restart2, 'cohort_root_N', cohort_root_N_ptr )
        call get_cohort_data(restart2, 'cohort_total_N', cohort_total_N_ptr )
        call get_cohort_data(restart2, 'nitrogen_stress', cohort_nitrogen_stress_ptr )
     endif
     if(field_exists(restart2,'myc_scav_marginal_gain_smoothed')) then
        call get_cohort_data(restart2, 'myc_scav_marginal_gain_smoothed',cohort_myc_scav_marginal_gain_smoothed_ptr)
        call get_cohort_data(restart2, 'myc_mine_marginal_gain_smoothed',cohort_myc_mine_marginal_gain_smoothed_ptr)
        call get_cohort_data(restart2, 'N_fix_marginal_gain_smoothed',cohort_N_fix_marginal_gain_smoothed_ptr)
        call get_cohort_data(restart2, 'rhiz_exud_marginal_gain_smoothed',cohort_rhiz_exud_marginal_gain_smoothed_ptr)

        call get_cohort_data(restart2, 'max_monthly_scav_alloc',cohort_max_monthly_scav_alloc_ptr)
        call get_cohort_data(restart2, 'max_monthly_mine_alloc',cohort_max_monthly_mine_alloc_ptr)
        call get_cohort_data(restart2, 'max_monthly_Nfix_alloc',cohort_max_monthly_Nfix_alloc_ptr)
        call get_cohort_data(restart2, 'max_scav_allocation',cohort_max_scav_allocation_ptr)
        call get_cohort_data(restart2, 'max_mine_allocation',cohort_max_mine_allocation_ptr)
        call get_cohort_data(restart2, 'max_Nfix_allocation',cohort_max_Nfix_allocation_ptr)
        call get_cohort_data(restart2, 'scav_alloc_accum',cohort_scav_alloc_accum_ptr)
        call get_cohort_data(restart2, 'mine_alloc_accum',cohort_mine_alloc_accum_ptr)
        call get_cohort_data(restart2, 'Nfix_alloc_accum',cohort_Nfix_alloc_accum_ptr)
     endif
     if(field_exists(restart2,'N_stress_smoothed')) &
        call get_cohort_data(restart2, 'N_stress_smoothed',cohort_nitrogen_stress_smoothed_ptr)

     do i = 0,nspecies-1
        if (field_exists(restart2,'drop_seed_C_'//trim(spdata(i)%name))) &
            call get_tile_data(restart2,'drop_seed_C_'//trim(spdata(i)%name),vegn_drop_seed_C_ptr,i)
        if (field_exists(restart2,'drop_seed_N_'//trim(spdata(i)%name))) &
            call get_tile_data(restart2,'drop_seed_N_'//trim(spdata(i)%name),vegn_drop_seed_N_ptr,i)
     enddo

     call get_int_tile_data(restart2,'landuse',vegn_landuse_ptr)

     call get_tile_data(restart2,'age',vegn_age_ptr)

     call get_tile_data(restart2,'fsc_pool_ag',vegn_fsc_pool_ag_ptr)
     call get_tile_data(restart2,'fsc_rate_ag',vegn_fsc_rate_ag_ptr)
     call get_tile_data(restart2,'fsc_pool_bg',vegn_fsc_pool_bg_ptr)
     call get_tile_data(restart2,'fsc_rate_bg',vegn_fsc_rate_bg_ptr)
     call get_tile_data(restart2,'ssc_pool_ag',vegn_ssc_pool_ag_ptr)
     call get_tile_data(restart2,'ssc_rate_ag',vegn_ssc_rate_ag_ptr)
     call get_tile_data(restart2,'ssc_pool_bg',vegn_ssc_pool_bg_ptr)
     call get_tile_data(restart2,'ssc_rate_bg',vegn_ssc_rate_bg_ptr)
     do j = 1,N_LITTER_POOLS
        do i = 1,N_C_TYPES-1 ! "-1" excludes deadmic (which is currently always 0) from restarts
           call get_tile_data(restart2,trim(l_shortname(j))//'litter_buffer_'//c_shortname(i),vegn_litter_buff_C_ptr,i,j)
           call get_tile_data(restart2,trim(l_shortname(j))//'litter_buffer_rate_'//c_shortname(i),vegn_litter_rate_C_ptr,i,j)
        enddo
     enddo
     if (soil_carbon_option==SOILC_CORPSE_N.and.field_exists(restart2,'fsn_pool_bg')) then
        call get_tile_data(restart2,'fsn_pool_bg',vegn_fsn_pool_bg_ptr)
        call get_tile_data(restart2,'fsn_rate_bg',vegn_fsn_rate_bg_ptr)
        call get_tile_data(restart2,'ssn_pool_bg',vegn_ssn_pool_bg_ptr)
        call get_tile_data(restart2,'ssn_rate_bg',vegn_ssn_rate_bg_ptr)
        do j = 1,N_LITTER_POOLS
           do i = 1,N_C_TYPES-1 ! "-1" excludes deadmic (which is currently always 0) from restarts
              call get_tile_data(restart2,trim(l_shortname(j))//'litter_buff_N_'//c_shortname(i),vegn_litter_buff_N_ptr,i,j)
              call get_tile_data(restart2,trim(l_shortname(j))//'litter_rate_N_'//c_shortname(i),vegn_litter_rate_N_ptr,i,j)
           enddo
        enddo
     endif
     ! monthly-mean values
     call get_tile_data(restart2,'tc_av', vegn_tc_av_ptr)
     if(field_exists(restart2,'theta_av_phen')) then
        call get_tile_data(restart2,'theta_av_phen', vegn_theta_av_phen_ptr)
        call get_tile_data(restart2,'theta_av_fire', vegn_theta_av_fire_ptr)
        call get_tile_data(restart2,'psist_av', vegn_psist_av_ptr)
     else
        call get_tile_data(restart2,'theta_av', vegn_theta_av_phen_ptr)
        call get_tile_data(restart2,'theta_av', vegn_theta_av_fire_ptr)
        ! psist_av remains at initial value (equal to 0)
     endif
     call get_tile_data(restart2,'tsoil_av', vegn_tsoil_av_ptr)
     call get_tile_data(restart2,'precip_av', vegn_precip_av_ptr)
     call get_tile_data(restart2,'lambda', vegn_lambda_ptr)
     call get_tile_data(restart2,'fuel', vegn_fuel_ptr)
     ! annual-mean values
     call get_tile_data(restart2,'t_ann', vegn_t_ann_ptr)
     call get_tile_data(restart2,'t_cold', vegn_t_cold_ptr)
     call get_tile_data(restart2,'p_ann', vegn_p_ann_ptr)
     call get_tile_data(restart2,'ncm', vegn_ncm_ptr)
     if (field_exists(restart2,'treeline_T_accum')) then
        call get_tile_data(restart2,'treeline_T_accum',vegn_treeline_T_accum_ptr)
        call get_tile_data(restart2,'treeline_N_accum',vegn_treeline_N_accum_ptr)
     endif
     ! accumulated values for annual averaging
     call get_tile_data(restart2,'t_ann_acm', vegn_t_ann_acm_ptr)
     call get_tile_data(restart2,'t_cold_acm', vegn_t_cold_acm_ptr)
     call get_tile_data(restart2,'p_ann_acm', vegn_p_ann_acm_ptr)
     call get_tile_data(restart2,'ncm_acm', vegn_ncm_acm_ptr)

     if (field_exists(restart2,'tc_pheno')) &
         call get_tile_data(restart2,'tc_pheno', vegn_tc_pheno_ptr)
     ! burned carbon pool and rate
     if(field_exists(restart2,'csmoke_pool')) &
          call get_tile_data(restart2,'csmoke_pool',vegn_csmoke_pool_ptr)
     if(field_exists(restart2,'csmoke_rate')) &
          call get_tile_data(restart2,'csmoke_rate',vegn_csmoke_rate_ptr)
     if(field_exists(restart2,'nsmoke_pool')) &
          call get_tile_data(restart2,'nsmoke_pool',vegn_nsmoke_pool_ptr)
     ! harvesting pools and rates
     do i = 1, N_HARV_POOLS
        if (field_exists(restart2,trim(HARV_POOL_NAMES(i))//'_harv_pool')) &
             call get_tile_data(restart2,trim(HARV_POOL_NAMES(i))//'_harv_pool',vegn_harv_pool_C_ptr,i)
        if (field_exists(restart2,trim(HARV_POOL_NAMES(i))//'_harv_rate')) &
             call get_tile_data(restart2,trim(HARV_POOL_NAMES(i))//'_harv_rate',vegn_harv_rate_C_ptr,i)
        if (field_exists(restart2,trim(HARV_POOL_NAMES(i))//'_harv_pool_nitrogen')) &
             call get_tile_data(restart2,trim(HARV_POOL_NAMES(i))//'_harv_pool_nitrogen',vegn_harv_pool_N_ptr,i)
     enddo
     ! read table of species names, if exists, and remap species as necessary
     call read_remap_species(restart2)
  else
     call error_mesg('vegn_init', 'cold-starting vegetation', NOTE)
  endif
  call free_land_restart(restart1)
  call free_land_restart(restart2)

  ! read climatological fields for initialization of species distribution
  if (file_exist('INPUT/biodata.nc'))then
     allocate(&
          t_ann (lnd%ls:lnd%le),&
          t_cold(lnd%ls:lnd%le),&
          p_ann (lnd%ls:lnd%le),&
          ncm   (lnd%ls:lnd%le) )
     call read_field( 'INPUT/biodata.nc','T_ANN',  t_ann,  interp='nearest')
     call read_field( 'INPUT/biodata.nc','T_COLD', t_cold, interp='nearest')
     call read_field( 'INPUT/biodata.nc','P_ANN',  p_ann,  interp='nearest')
     call read_field( 'INPUT/biodata.nc','NCM',    ncm,    interp='nearest')
     did_read_biodata = .TRUE.
     call error_mesg('vegn_init','did read INPUT/biodata.nc',NOTE)
  else
     did_read_biodata = .FALSE.
     call error_mesg('vegn_init','did NOT read INPUT/biodata.nc',NOTE)
  endif

  ! create a list of species indices for initialization
  init_cohort_spp(:) = -1
  do n = 1,init_n_cohorts
     do i = 0,size(spdata)-1
        if (trim(init_cohort_species(n))==trim(spdata(i)%name)) &
            init_cohort_spp(n) = i
     enddo
  enddo

  ! Go through all tiles and initialize the cohorts that have not been initialized yet --
  ! this allows to read partial restarts. Also initialize accumulation counters to zero
  ! or the values from the restarts.
  ce = first_elmt(land_tile_map, ls=lnd%ls)
  do while(loop_over_tiles(ce,tile,l))
     if (.not.associated(tile%vegn)) cycle

     tile%vegn%n_accum = n_accum
     tile%vegn%nmn_acm = nmn_acm

     if (tile%vegn%n_cohorts>0) cycle ! skip initialized tiles

     ! create and initialize cohorts for this vegetation tile
     tile%vegn%n_cohorts = init_n_cohorts
     tile%vegn%tc_pheno  = init_Tv  ! initial temperature for phenology

     allocate(tile%vegn%cohorts(tile%vegn%n_cohorts))
     do n = 1,tile%vegn%n_cohorts
        associate(cc => tile%vegn%cohorts(n))
        cc%Wl = init_Wl
        cc%Ws = init_Ws
        cc%Tv = init_Tv

        cc%bl      = init_cohort_bl(n)
        cc%blv     = init_cohort_blv(n)
        cc%br      = init_cohort_br(n)
        cc%bsw     = init_cohort_bsw(n)
        cc%bwood   = init_cohort_bwood(n)
        cc%bseed   = init_cohort_bseed(n)
        cc%scav_myc_C_reservoir = 0.0 ; cc%scav_myc_N_reservoir = 0.0
        cc%mine_myc_C_reservoir = 0.0 ; cc%mine_myc_N_reservoir = 0.0
        cc%N_fixer_C_reservoir  = 0.0 ; cc%N_fixer_N_reservoir  = 0.0

        cc%nindivs = init_cohort_nindivs(n)
        cc%age     = init_cohort_age(n)
        cc%bliving = cc%bl+cc%br+cc%blv+cc%bsw
        cc%npp_previous_day = 0.0
        cc%status  = LEAF_ON
        cc%leaf_age = 0.0

        cc%myc_scav_marginal_gain_smoothed = 0.0
        cc%myc_mine_marginal_gain_smoothed = 0.0
        cc%N_fix_marginal_gain_smoothed = 0.0
        cc%rhiz_exud_marginal_gain_smoothed = 0.0
        cc%max_monthly_scav_alloc = 0.0
        cc%max_monthly_mine_alloc = 0.0
        cc%max_monthly_Nfix_alloc = 0.0
        cc%max_scav_allocation = 0.0
        cc%max_mine_allocation = 0.0
        cc%max_Nfix_allocation = 0.0
        cc%scav_alloc_accum = 0.0
        cc%mine_alloc_accum = 0.0
        cc%Nfix_alloc_accum = 0.0
        cc%nitrogen_stress_smoothed = 1.0

        if (do_ppa) then
           cc%species = init_cohort_spp(n)
           if (cc%species < 0) call error_mesg('vegn_init','species "'//trim(init_cohort_species(n))//&
                   '" needed for initialization, but not found in the list of species parameters', FATAL)
           call init_cohort_allometry_ppa(cc, init_cohort_height(n), init_cohort_nsc_frac(n), init_cohort_nsn_frac(n))
           call init_cohort_hydraulics(tile%vegn%cohorts(n), tile%soil%pars%psi_sat_ref) ! adam wolf
           ! initialize DBH_ys
           cc%DBH_ys = cc%dbh
           cc%BM_ys  = cc%bsw + cc%bwood
        else
           if(did_read_biodata.and.do_biogeography) then
              call update_species(cc,t_ann(l),t_cold(l),p_ann(l),ncm(l),LU_NTRL)
              tile%vegn%t_ann  = t_ann (l)
              tile%vegn%t_cold = t_cold(l)
              tile%vegn%p_ann  = p_ann (l)
              tile%vegn%ncm    = ncm   (l)
           else
              cc%species = tile%vegn%tag
           endif
           associate(sp=>spdata(cc%species))
           cc%leaf_N     = (cc%bl+cc%blv)/sp%leaf_live_c2n
           cc%wood_N     = cc%bwood/sp%wood_c2n
           cc%sapwood_N  = cc%bsw/sp%sapwood_c2n
           cc%root_N     = cc%br/sp%froot_live_c2n
           cc%stored_N   = init_cohort_stored_N_mult(n)*(cc%leaf_N+cc%root_N)
           cc%total_N    = cc%stored_N+cc%leaf_N+cc%wood_N+cc%root_N+cc%sapwood_N
           end associate ! sp
        endif

        ! in ppa mode, init_cohort_allometry_ppa sets myc and fixer biomasses to zero, but
        ! we override them with values from namelist here
        associate(sp=>spdata(cc%species))
        cc%myc_scavenger_biomass_C = init_cohort_myc_scav(n)
        cc%myc_scavenger_biomass_N = cc%myc_scavenger_biomass_C/sp%c2n_mycorrhizae
        cc%myc_miner_biomass_C     = init_cohort_myc_mine(n)
        cc%myc_miner_biomass_N     = cc%myc_miner_biomass_C/sp%c2n_mycorrhizae
        cc%N_fixer_biomass_C       = init_cohort_n_fixer(n)
        cc%N_fixer_biomass_N       = cc%N_fixer_biomass_C/c2n_N_fixer

        cc%K_r = spdata(cc%species)%root_perm
        end associate ! sp
        end associate ! cc
     enddo
     if (do_ppa) &
        call vegn_relayer_cohorts_ppa(tile%vegn) ! this can change the number of cohorts
  enddo

  ! initialize carbon integrator
  call vegn_dynamics_init ( id_ug, lnd%time, delta_time )
  call vegn_disturbance_init ( id_ug )

  ! initialize static vegetation
  call static_vegn_init ()
  call read_static_vegn ( lnd%time )

  ! initialize harvesting options
  call vegn_harvesting_init(id_ug)

  ! initialize fire
  call vegn_fire_init(id_ug, id_cellarea, delta_time, lnd%time)

  ! initialize vegetation diagnostic fields
  call vegn_diag_init ( id_ug, id_band, lnd%time )

  ! ---- diagnostic section
  ce = first_elmt(land_tile_map, ls=lnd%ls)
  do while(loop_over_tiles(ce,tile))
     if (.not.associated(tile%vegn)) cycle ! skip non-vegetation tiles
     ! send the data
     call send_tile_data(id_vegn_type,  real(tile%vegn%tag), tile%diag)
  enddo

  if (allocated(t_ann))  deallocate(t_ann)
  if (allocated(t_cold)) deallocate(t_cold)
  if (allocated(p_ann))  deallocate(p_ann)
  if (allocated(ncm))    deallocate(ncm)

end subroutine vegn_init

! ============================================================================
subroutine vegn_diag_init ( id_ug, id_band, time )
  integer        , intent(in) :: id_ug   !<Unstructured axis id.
  integer        , intent(in) :: id_band ! ID of spectral band axis
  type(time_type), intent(in) :: time    ! initial time for diagnostic fields

  ! ---- local vars
  integer :: i

  ! set the default sub-sampling filter for the fields below
  call set_default_diag_filter('soil')

  id_vegn_type = register_tiled_static_field ( module_name, 'vegn_type',  &
       (/id_ug/), 'vegetation type', missing_value=-1.0 )

  id_ncohorts = register_cohort_diag_field( module_name, 'ncohorts', &
       (/id_ug/), time, 'number of cohorts', 'unitless', missing_value=-1.0)
  id_nindivs = register_cohort_diag_field( module_name, 'nindivs', &
       (/id_ug/), time, 'density of individuals', 'individuals/m2', missing_value=-1.0)

  id_nlayers = register_tiled_diag_field( module_name, 'nlayers', &
       (/id_ug/), time, 'number of canopy layers', 'unitless', missing_value=-1.0 )
  id_dbh = register_cohort_diag_field( module_name, 'dbh', &
       (/id_ug/), time, 'diameter at breast height', 'm', missing_value=-1.0)
  id_dbh_max = register_cohort_diag_field( module_name, 'dbh_max', &
       (/id_ug/), time, 'maximum diameter at breast height', 'm', missing_value=-1.0)
  id_crownarea = register_cohort_diag_field( module_name, 'crownarea', &
       (/id_ug/), time, 'mean area of individuals crown', 'm2', missing_value=-1.0)

  id_temp = register_cohort_diag_field ( module_name, 'temp',  &
       (/id_ug/), time, 'canopy temperature', 'degK', missing_value=-1.0)
  id_wl = register_cohort_diag_field ( module_name, 'wl',  &
       (/id_ug/), time, 'canopy liquid water content', 'kg/m2', missing_value=-1.0)
  id_ws = register_cohort_diag_field ( module_name, 'ws',  &
       (/id_ug/), time, 'canopy solid water content', 'kg/m2', missing_value=-1.0)


  id_height = register_tiled_diag_field ( module_name, 'height',  &
       (/id_ug/), time, 'height of tallest vegetation', 'm', missing_value=-1.0)
  id_height_ave = register_cohort_diag_field ( module_name, 'height_ave',  &
       (/id_ug/), time, 'average height of the trees', 'm', missing_value=-1.0)

  id_lai    = register_cohort_diag_field ( module_name, 'lai',  &
       (/id_ug/), time, 'leaf area index', 'm2/m2', missing_value=-1.0)
  id_lai_var = register_cohort_diag_field ( module_name, 'lai_var',  &
       (/id_ug/), time, 'variance of leaf area index across tiles in grid cell', 'm4/m4', &
       missing_value=-1.0 , opt='variance')
  id_lai_std = register_cohort_diag_field ( module_name, 'lai_std',  &
       (/id_ug/), time, 'standard deviation of leaf area index across tiles in grid cell', 'm2/m2', &
       missing_value=-1.0, opt='stdev')
  id_sai    = register_cohort_diag_field ( module_name, 'sai',  &
       (/id_ug/), time, 'stem area index', 'm2/m2', missing_value=-1.0)
  id_leafarea = register_cohort_diag_field ( module_name, 'leafarea',  &
       (/id_ug/), time, 'leaf area per individual', 'm2', missing_value=-1.0)
  id_laii   = register_cohort_diag_field ( module_name, 'laii',  &
       (/id_ug/), time, 'leaf area index per individual', 'm2/m2', missing_value=-1.0)

  id_leaf_size = register_cohort_diag_field ( module_name, 'leaf_size',  &
       (/id_ug/), time, missing_value=-1.0 )
  id_root_density = register_cohort_diag_field ( module_name, 'root_density',  &
       (/id_ug/), time, 'total biomass below ground', 'kg/m2', missing_value=-1.0 )
  id_root_zeta = register_cohort_diag_field ( module_name, 'root_zeta',  &
       (/id_ug/), time, 'e-folding depth of root biomass', 'm',missing_value=-1.0 )
  id_rs_min = register_cohort_diag_field ( module_name, 'rs_min',  &
       (/id_ug/), time, missing_value=-1.0 )
!   id_leaf_refl = register_cohort_diag_field ( module_name, 'leaf_refl',  &
!        (/id_ug,id_band/), time, 'reflectivity of leaf', missing_value=-1.0 )
!   id_leaf_tran = register_cohort_diag_field ( module_name, 'leaf_tran',  &
!        (/id_ug,id_band/), time, 'transmittance of leaf', missing_value=-1.0 )
  id_leaf_emis = register_cohort_diag_field ( module_name, 'leaf_emis',  &
       (/id_ug/), time, 'leaf emissivity', missing_value=-1.0 )
  id_snow_crit = register_cohort_diag_field ( module_name, 'snow_crit',  &
       (/id_ug/), time, missing_value=-1.0 )
  id_stomatal = register_tiled_diag_field ( module_name, 'stomatal_cond',  &
       (/id_ug/), time, 'vegetation stomatal conductance', 'm/s', missing_value=-1.0 )
  id_evap_demand = register_tiled_diag_field ( module_name, 'evap_demand',  &
       (/id_ug/), time, 'plant evaporative water demand',&
       'kg/(m2 s)', missing_value=-1e20 )
  id_an_op = register_cohort_diag_field ( module_name, 'an_op',  &
       (/id_ug/), time, 'net photosynthesis with open stomata', &
       '(mol CO2)(m2 of leaf)^-1 year^-1', missing_value=-1e20)
  id_an_cl = register_cohort_diag_field ( module_name, 'an_cl',  &
       (/id_ug/), time, 'net photosynthesis with closed stomata', &
       '(mol CO2)(m2 of leaf)^-1 year^-1', missing_value=-1e20)
  id_psi_r  = register_cohort_diag_field ( module_name, 'psi_r',  &
       (/id_ug/), time, 'root water potential', 'MPa', missing_value=-1e20)
  id_psi_x  = register_cohort_diag_field ( module_name, 'psi_x',  &
       (/id_ug/), time, 'stem water potential', 'MPa', missing_value=-1e20)
  id_psi_l  = register_cohort_diag_field ( module_name, 'psi_l', &
       (/id_ug/), time, 'leaf water potential', 'MPa', missing_value=-1e20)
  id_Kxi  = register_cohort_diag_field ( module_name, 'Kxi',  &
       (/id_ug/), time, 'stem conductance', 'kg/(indiv s MPa)', missing_value=-1e20)
  id_Kli  = register_cohort_diag_field ( module_name, 'Kli',  &
       (/id_ug/), time, 'leaf conductance', 'kg/(indiv s MPa)', missing_value=-1e20)
  id_w_scale  = register_cohort_diag_field ( module_name, 'w_scale',  &
       (/id_ug/), time, 'reduction of stomatal conductance due to water stress', 'unitless', &
       missing_value=-1e20)
  id_RHi  = register_cohort_diag_field ( module_name, 'RHi',  &
       (/id_ug/), time, 'relative humidity inside leaf', 'percent', missing_value=-1e20)

  id_bl = register_cohort_diag_field ( module_name, 'bl',  &
       (/id_ug/), time, 'biomass of leaves', 'kg C/m2', missing_value=-1.0)
  id_blv = register_cohort_diag_field ( module_name, 'blv',  &
       (/id_ug/), time, 'biomass in labile store', 'kg C/m2', missing_value=-1.0)
  id_br = register_cohort_diag_field ( module_name, 'br',  &
       (/id_ug/), time, 'biomass of fine roots', 'kg C/m2', missing_value=-1.0)
  id_bsw = register_cohort_diag_field ( module_name, 'bsw',  &
       (/id_ug/), time, 'biomass of sapwood', 'kg C/m2', missing_value=-1.0)
  id_bwood = register_cohort_diag_field ( module_name, 'bwood',  &
       (/id_ug/), time, 'biomass of heartwood', 'kg C/m2', missing_value=-1.0)
  id_btot = register_cohort_diag_field ( module_name, 'btot',  &
       (/id_ug/), time, 'total biomass', 'kg C/m2', missing_value=-1.0)
  id_bseed = register_cohort_diag_field ( module_name, 'bseed',  &
       (/id_ug/), time, 'biomass of seed', 'kg C/m2', missing_value=-1.0)
  id_nsc = register_cohort_diag_field ( module_name, 'nsc',  &
       (/id_ug/), time, 'biomass in non-structural pool', 'kg C/m2', missing_value=-1.0)

  id_leaf_N = register_cohort_diag_field ( module_name, 'Nl',  &
       (/id_ug/), time, 'nitrogen content of leaves', 'kg N/m2', missing_value=-1.0 )
  id_root_N = register_cohort_diag_field ( module_name, 'Nr',  &
       (/id_ug/), time, 'nitrogen content of fine roots', 'kg N/m2', missing_value=-1.0 )
  id_wood_N = register_cohort_diag_field ( module_name, 'Nwood',  &
       (/id_ug/), time, 'nitrogen content of wood', 'kg N/m2', missing_value=-1.0 )
  id_sapwood_N = register_cohort_diag_field ( module_name, 'Nsw',  &
       (/id_ug/), time, 'nitrogen content of sapwood', 'kg N/m2', missing_value=-1.0 )
  id_stored_N = register_cohort_diag_field ( module_name, 'Nstore',  &
       (/id_ug/), time, 'vegetation nitrogen storage', 'kg N/m2', missing_value=-1.0 )
  id_Ntot = register_cohort_diag_field ( module_name, 'Nvegtot',  &
       (/id_ug/), time, 'total vegetation nitrogen content', 'kg N/m2', missing_value=-1.0 )
  id_seed_N = register_cohort_diag_field ( module_name, 'Nseed',  &
       (/id_ug/), time, 'nitrogen content of seeds', 'kg N/m2', missing_value=-1.0)


  id_mrz_scav_C = register_cohort_diag_field ( module_name, 'mrz_scav_biomass_C',  &
       (/id_ug/), time, 'scavenger mycorrhizal biomass C', 'kg C/m2', missing_value=-1.0 )
  id_mrz_mine_C = register_cohort_diag_field ( module_name, 'mrz_mine_biomass_C',  &
       (/id_ug/), time, 'miner mycorrhizal biomass C', 'kg C/m2', missing_value=-1.0 )
  id_Nfix_C = register_cohort_diag_field ( module_name, 'Nfix_biomass_C',  &
       (/id_ug/), time, 'symbiotic N fixer biomass C', 'kg C/m2', missing_value=-1.0 )
  id_mrz_scav_N = register_cohort_diag_field ( module_name, 'mrz_scav_biomass_N',  &
       (/id_ug/), time, 'scavenger mycorrhizal biomass N', 'kg N/m2', missing_value=-1.0 )
  id_mrz_mine_N = register_cohort_diag_field ( module_name, 'mrz_mine_biomass_N',  &
       (/id_ug/), time, 'miner mycorrhizal biomass N', 'kg N/m2', missing_value=-1.0 )
  id_Nfix_N = register_cohort_diag_field ( module_name, 'Nfix_biomass_N',  &
       (/id_ug/), time, 'symbiotic N fixer biomass N', 'kg N/m2', missing_value=-1.0 )

  id_bsw_max = register_cohort_diag_field ( module_name, 'bsw_max',  &
       (/id_ug/), time, 'max biomass of sapwood', 'kg C/m2', missing_value=-1.0)
  id_bl_max = register_cohort_diag_field ( module_name, 'bl_max',  &
       (/id_ug/), time, 'max biomass of leaves', 'kg C/m2', missing_value=-1.0)
  id_br_max = register_cohort_diag_field ( module_name, 'br_max',  &
       (/id_ug/), time, 'max biomass of fine roots', 'kg C/m2', missing_value=-1.0)
  ! ens 021617
  id_brsw = register_cohort_diag_field ( module_name, 'brsw',  &
       (/id_ug/), time, 'biomass of branches (only sapwood)', 'kg C/m2', missing_value=-1.0)
  id_growth_prev_day = register_cohort_diag_field ( module_name, 'growth_prev_day',  &
       (/id_ug/), time, 'growth previous day', 'kg C/m2', missing_value=-1.0)
  id_topyear = register_cohort_diag_field ( module_name, 'topyear',  &
       (/id_ug/), time, 'time pants spent in top layer', 'year', missing_value=-1.0)

  id_fuel = register_tiled_diag_field ( module_name, 'fuel',  &
       (/id_ug/), time, 'mass of fuel', 'kg C/m2', missing_value=-1.0 )
  id_lambda = register_tiled_diag_field (module_name, 'lambda',(/id_ug/), &
       time, 'drought', 'months', missing_value=-100.0)

  !ppg 2017-11-03
  id_lai_kok = register_cohort_diag_field ( module_name, 'lai_kok',  &
       (/id_ug/), time, 'leaf area index at Kok effect threshold', 'm2/m2', missing_value=-1.0 )

  id_DanDlai = register_cohort_diag_field ( module_name, 'DanDlai',  &
       (/id_ug/), time, 'derivative of photosynthesis w.r.t. LAI', missing_value=-1.0 )

  id_PAR_dn = register_cohort_diag_field ( module_name, 'PAR_dn',  &
       (/id_ug/), time, 'downward flux of photosynthetically-active radiation on top of the canopy', &
       'W/m2', missing_value=-1.0 )
  id_PAR_net = register_cohort_diag_field ( module_name, 'PAR_net',  &
       (/id_ug/), time, 'net flux of photosynthetically-active radiation to the canopy', &
       'W/m2', missing_value=-1.0 )

  id_species = register_tiled_diag_field ( module_name, 'species',  &
       (/id_ug/), time, 'vegetation species number', missing_value=-1.0 )
  id_status = register_tiled_diag_field ( module_name, 'status',  &
       (/id_ug/), time, 'status of leaves', missing_value=-1.0 )
  id_theph = register_tiled_diag_field ( module_name, 'theph',  &
       (/id_ug/), time, 'theta for phenology', missing_value=-1.0 )
  id_psiph = register_tiled_diag_field ( module_name, 'psiph',  &
       (/id_ug/), time, 'psi stress for phenology', missing_value=-1.0 )
  id_leaf_age = register_cohort_diag_field ( module_name, 'leaf_age',  &
       (/id_ug/), time, 'age of leaves since bud burst', 'days', missing_value=-1.0 )!ens

  id_con_v_h = register_tiled_diag_field ( module_name, 'con_v_h', (/id_ug/), &
       time, 'conductance for sensible heat between canopy and canopy air', &
       'm/s', missing_value=-1.0 )
  id_con_v_v = register_tiled_diag_field ( module_name, 'con_v_v', (/id_ug/), &
       time, 'conductance for water vapor between canopy and canopy air', &
       'm/s', missing_value=-1.0 )
  id_soil_water_supply = register_tiled_diag_field ( module_name, 'soil_water_supply', &
       (/id_ug/), time, 'maximum rate of soil water supply to vegetation', &
       'kg/(m2 s)', missing_value=-1e20)

  id_cgain = register_tiled_diag_field ( module_name, 'cgain', (/id_ug/), &
       time, 'carbon gain', 'kg C', missing_value=-100.0 )
  id_closs = register_tiled_diag_field ( module_name, 'closs', (/id_ug/), &
       time, 'carbon loss', 'kg C', missing_value=-100.0 )
  id_wdgain = register_tiled_diag_field ( module_name, 'wdgain', (/id_ug/), &
       time, 'wood biomass gain', 'kg C', missing_value=-100.0 )

  id_Ngain = register_tiled_diag_field ( module_name, 'Ngain', (/id_ug/), &
       time, 'nitrogen gain', 'kg N', missing_value=-100.0 )
  id_Nloss = register_tiled_diag_field ( module_name, 'Nloss', (/id_ug/), &
       time, 'nitrogen loss', 'kg N', missing_value=-100.0 )

  id_t_ann  = register_tiled_diag_field ( module_name, 't_ann', (/id_ug/), &
       time, 'annual mean temperature', 'degK', missing_value=-999.0 )
  id_t_cold  = register_tiled_diag_field ( module_name, 't_cold', (/id_ug/), &
       time, 'average temperature of the coldest month', 'degK', missing_value=-999.0 )
  id_p_ann  = register_tiled_diag_field ( module_name, 'p_ann', (/id_ug/), &
       time, 'annual mean precipitation', 'kg/(m2 s)', missing_value=-999.0 )
  id_ncm = register_tiled_diag_field ( module_name, 'ncm', (/id_ug/), &
       time, 'number of cold months', 'dimensionless', missing_value=-999.0 )

  id_tot_harv_pool_C = register_tiled_diag_field( module_name, 'harv_pool_C', (/id_ug/), &
       time, 'total harvested carbon', 'kg C/m2', missing_value=-999.0)
  id_tot_harv_rate_C = register_tiled_diag_field( module_name, 'harv_rate_C', (/id_ug/), &
       time, 'total rate of release of harvested carbon to the atmosphere', &
       'kg C/(m2 year)', missing_value=-999.0)
  id_tot_harv_pool_N = register_tiled_diag_field( module_name, 'harv_pool_N', (/id_ug/), &
       time, 'total harvested nitrogen', 'kg N/m2', missing_value=-999.0)
  do i = 1,N_HARV_POOLS
     id_harv_pool_C(i) = register_tiled_diag_field( module_name, &
          trim(HARV_POOL_NAMES(i))//'_harv_pool_C', (/id_ug/), time, &
          'harvested carbon', 'kg C/m2', missing_value=-999.0)
     id_harv_rate_C(i) = register_tiled_diag_field( module_name, &
          trim(HARV_POOL_NAMES(i))//'_harv_rate_C', (/id_ug/), time, &
          'rate of release of harvested carbon to the atmosphere', 'kg C/(m2 year)', &
          missing_value=-999.0)
     id_harv_pool_N(i) = register_tiled_diag_field( module_name, &
           trim(HARV_POOL_NAMES(i))//'_harv_pool_N', (/id_ug/), time, &
           'harvested nitrogen', 'kg N/m2', missing_value=-999.0)
  enddo

  id_litter_buff_C(:,:) = register_litter_soilc_diag_fields(module_name, '<ltype>litter_buff_C_<ctype>', (/id_ug/), &
       time, 'intermediate pool of <ltype> <ctype> litter carbon', 'kg C/m2', missing_value=-999.0)
  id_litter_buff_N(:,:) = register_litter_soilc_diag_fields(module_name, '<ltype>litter_buff_N_<ctype>', (/id_ug/), &
       time, 'intermediate pool of <ltype> <ctype> litter nitrogen', 'kg N/m2', missing_value=-999.0)
  id_litter_rate_C(:,:) = register_litter_soilc_diag_fields(module_name, '<ltype>litter_rate_C_<ctype>', (/id_ug/), &
       time, 'rate of conversion of <ltype> litter buffer to the <ctype> soil carbon', 'kg C/(m2 yr)', missing_value=-999.0)
  id_litter_rate_N(:,:) = register_litter_soilc_diag_fields(module_name, '<ltype>litter_rate_N_<ctype>', (/id_ug/), &
       time, 'rate of conversion of <ltype> litter buffer to the <ctype> soil nitrogen', 'kg N/(m2 yr)', missing_value=-999.0)

  id_fsc_pool_ag = register_tiled_diag_field (module_name, 'fsc_pool_ag', (/id_ug/), &
       time, 'intermediate pool of above-ground fast soil carbon', 'kg C/m2', missing_value=-999.0)
  id_fsc_rate_ag = register_tiled_diag_field (module_name, 'fsc_rate_ag', (/id_ug/), &
       time, 'rate of conversion of above-ground fsc_pool to the fast soil_carbon', 'kg C/(m2 yr)', &
       missing_value=-999.0)
  id_ssc_pool_ag = register_tiled_diag_field (module_name, 'ssc_pool_ag', (/id_ug/), &
       time, 'intermediate pool of above-ground slow soil carbon', 'kg C/m2', missing_value=-999.0)
  id_ssc_rate_ag = register_tiled_diag_field (module_name, 'ssc_rate_ag', (/id_ug/), &
       time, 'rate of conversion of above-ground ssc_pool to the fast soil_carbon', 'kg C/(m2 yr)', &
       missing_value=-999.0)

  id_fsc_pool_bg = register_tiled_diag_field (module_name, 'fsc_pool_bg', (/id_ug/), &
       time, 'intermediate pool of below-ground fast soil carbon', 'kg C/m2', missing_value=-999.0)
  id_fsc_rate_bg = register_tiled_diag_field (module_name, 'fsc_rate_bg', (/id_ug/), &
       time, 'rate of conversion of below-ground fsc_pool to the fast soil_carbon', 'kg C/(m2 yr)', &
       missing_value=-999.0)
  id_ssc_pool_bg = register_tiled_diag_field (module_name, 'ssc_pool_bg', (/id_ug/), &
       time, 'intermediate pool of below-ground slow soil carbon', 'kg C/m2', missing_value=-999.0)
  id_ssc_rate_bg = register_tiled_diag_field (module_name, 'ssc_rate_bg', (/id_ug/), &
       time, 'rate of conversion of below-ground ssc_pool to the fast soil_carbon', 'kg C/(m2 yr)', &
       missing_value=-999.0)

  id_csmoke_pool = register_tiled_diag_field ( module_name, 'csmoke', (/id_ug/), &
       time, 'carbon lost through fire', 'kg C/m2', missing_value=-999.0)
  id_csmoke_rate = register_tiled_diag_field ( module_name, 'csmoke_rate', (/id_ug/), &
       time, 'rate of release of carbon lost through fire to the atmosphere', &
       'kg C/(m2 yr)', missing_value=-999.0)
  id_nsmoke_pool = register_tiled_diag_field ( module_name, 'Nsmoke', (/id_ug/), &
       time, 'nitrogen lost through fire', 'kg N/m2', missing_value=-999.0)

  id_ssc_in = register_tiled_diag_field ( module_name, 'ssc_in',  (/id_ug/), &
     time,  'soil slow carbon in', 'kg C/m2', missing_value=-999.0 )
  id_ssc_out = register_tiled_diag_field ( module_name, 'ssc_out',  (/id_ug/), &
     time,  'soil slow carbon out', 'kg C/m2', missing_value=-999.0 )
  id_deadmic_out = register_tiled_diag_field ( module_name, 'deadmic_out',  (/id_ug/), &
     time,  'daed microbe carbon out', 'kg C/m2', missing_value=-999.0 )
  id_fsc_in = register_tiled_diag_field ( module_name, 'fsc_in',  (/id_ug/), &
     time,  'soil fast carbon in', 'kg C/m2', missing_value=-999.0 )
  id_fsc_out = register_tiled_diag_field ( module_name, 'fsc_out',  (/id_ug/), &
     time,  'soil fast carbon out', 'kg C/m2', missing_value=-999.0 )
  id_veg_in = register_tiled_diag_field ( module_name, 'veg_in',  (/id_ug/), &
     time,  'vegetation carbon in', 'kg C/m2', missing_value=-999.0 )
  id_veg_out = register_tiled_diag_field ( module_name, 'veg_out',  (/id_ug/), &
     time,  'vegetation carbon out', 'kg C/m2', missing_value=-999.0 )

  id_tile_nitrogen_gain = register_tiled_diag_field ( module_name, 'nitrogen_in',  (/id_ug/), &
     time,  'tile nitrogen in', 'kg N/m2', missing_value=-999.0 )
  id_tile_nitrogen_loss = register_tiled_diag_field ( module_name, 'nitrogen_out',  (/id_ug/), &
     time,  'tile nitrogen out', 'kg N/m2', missing_value=-999.0 )

  id_afire = register_tiled_diag_field (module_name, 'afire', (/id_ug/), &
       time, 'area been fired', missing_value=-100.0)
  id_atfall = register_tiled_diag_field (module_name, 'atfall',(/id_ug/), &
       time, 'area been disturbed', missing_value=-100.0)

  id_gdd = register_tiled_diag_field (module_name, 'gdd', (/id_ug/), &
       time, 'growing degree days','degK day', missing_value=-999.0)
  id_tc_pheno = register_tiled_diag_field (module_name, 'tc_pheno', (/id_ug/), &
       time, 'smoothed canopy air temperature for phenology','degK', missing_value=-999.0)

  id_phot_co2 = register_tiled_diag_field (module_name, 'qco2_phot',(/id_ug/), &
       time, 'CO2 mixing ratio for photosynthesis calculations', 'mol CO2/mol dry air', &
       missing_value=-1.0)

  id_zstar_1 = register_tiled_diag_field (module_name, 'zstar:C',(/id_ug/), &
       time, 'critical depth for the top layer', 'm', &
       missing_value=-1.0)

  ! CMOR/CMIP variables
  call set_default_diag_filter('land')

  id_lai_cmor = register_tiled_diag_field(cmor_name, 'lai', (/id_ug/), time, &
       'Leaf Area Index', '1.0', standard_name = 'leaf_area_index', &
       missing_value = -1.0, fill_missing = .TRUE.)
  call add_tiled_diag_field_alias(id_lai_cmor, cmor_name, 'laiLut', (/id_ug/), time, &
       'leaf area index on land use tile', '1.0', standard_name = 'leaf_area_index', &
       missing_value = -1.0, fill_missing = .FALSE.)

  id_cVeg = register_tiled_diag_field (cmor_name, 'cVeg', (/id_ug/), time, &
       'Carbon Mass in Vegetation', 'kg m-2', standard_name='vegetation_carbon_content', &
       missing_value = -1.0, fill_missing=.TRUE.)
  call add_tiled_diag_field_alias (id_cVeg, cmor_name, 'cVegLut', (/id_ug/), time, &
       'Carbon Mass in Vegetation', 'kg m-2', standard_name='vegetation_carbon_content', &
       missing_value = -1.0, fill_missing=.FALSE.)

  id_cLeaf = register_tiled_diag_field ( cmor_name, 'cLeaf',  (/id_ug/), &
       time, 'Carbon Mass in Leaves', 'kg m-2', missing_value=-1.0, &
       standard_name='leaf_carbon_content', fill_missing=.TRUE.)
  id_cWood = register_tiled_diag_field ( cmor_name, 'cWood',  (/id_ug/), &
       time, 'Carbon Mass in Wood', 'kg m-2', missing_value=-1.0, &
       standard_name='wood_carbon_content', fill_missing=.TRUE.)
  id_cRoot = register_tiled_diag_field ( cmor_name, 'cRoot',  (/id_ug/), &
       time, 'Carbon Mass in Roots', 'kg m-2', missing_value=-1.0, &
       standard_name='root_carbon_content', fill_missing=.TRUE.)
  id_cMisc = register_tiled_diag_field ( cmor_name, 'cMisc',  (/id_ug/), &
       time, 'Carbon Mass in Other Living Compartments on Land', 'kg m-2', missing_value=-1.0, &
       standard_name='miscellaneous_living_matter_carbon_content', fill_missing=.TRUE.)

  call add_tiled_diag_field_alias(id_height, cmor_name, 'vegHeight', (/id_ug/), &
       time, 'Vegetation height averaged over all vegetation types and over the vegetated fraction of a grid cell.', &
       'm', missing_value=-1.0, standard_name='canopy_height', fill_missing=.TRUE.)

  id_cProduct = register_tiled_diag_field( cmor_name, 'cProduct', (/id_ug/), &
       time, 'Carbon Mass in Products of Landuse Change', 'kg m-2', missing_value=-999.0, &
       standard_name='carbon_content_of_products_of_anthropogenic_land_use_change', fill_missing=.TRUE.)
  call add_tiled_diag_field_alias(id_cProduct, cmor_name, 'cProductLut', (/id_ug/), &
       time, 'wood and agricultural product pool carbon associated with land use tiles; examples of products include paper, cardboard, timber for construction, and crop harvest for food or fuel.', &
       'kg m-2', missing_value=-999.0, standard_name='carbon_content_in_wood_and_agricultural_products', &
       fill_missing = .FALSE.)

  id_cAnt = register_tiled_diag_field( cmor_name, 'cAnt', (/id_ug/), &
       time, 'Carbon in Anthropogenic Pool', 'kg m-2', missing_value=-999.0, &
       fill_missing=.TRUE.) ! standard_name not known at this time
  call add_tiled_diag_field_alias(id_cAnt, cmor_name, 'cAntLut', (/id_ug/), &
       time, 'Carbon in Anthropogenic Pools Associated with Land Use Tiles', 'kg m-2', &
       missing_value=-999.0, fill_missing = .FALSE.) ! standard_name not known at this time

  id_fGrazing = register_tiled_diag_field( cmor_name, 'fGrazing', (/id_ug/), &
       time, 'Carbon Mass Flux into Atmosphere due to Grazing on Land', 'kg m-2 s-1', missing_value=-1.0, &
       standard_name='surface_upward_mass_flux_of_carbon_dioxide_expressed_as_carbon_due_to_emission_from_grazing', &
       fill_missing=.TRUE.)
  id_fHarvest = register_tiled_diag_field( cmor_name, 'fHarvest', (/id_ug/), &
       time, 'Carbon Mass Flux into Atmosphere due to Crop Harvesting', 'kg m-2 s-1', missing_value=-1.0, &
       standard_name='surface_upward_mass_flux_of_carbon_dioxide_expressed_as_carbon_due_to_emission_from_crop_harvesting', &
       fill_missing=.TRUE.)
  id_fProductDecomp = register_tiled_diag_field( cmor_name, 'fProductDecomp', (/id_ug/), &
       time, 'decomposition out of product pools to CO2 in atmos', 'kg m-2 s-1', missing_value=-1.0, &
       standard_name='Carbon_flux_out_of_storage_product_pools_into_atmos', &
       fill_missing=.TRUE.)
  call add_tiled_diag_field_alias (id_fProductDecomp, cmor_name, 'fProductDecompLut', (/id_ug/), &
       time, 'flux from wood and agricultural product pools on land use tile into atmosphere', 'kg m-2 s-1', missing_value=-1.0, &
       standard_name='tendency_of_atmospheric_mass_content_of_carbon_dioxide_expressed_as_carbon_due_to_emission_from_wood_and_agricultural_product_pool', &
       fill_missing=.FALSE.)
  id_fLuc = register_tiled_diag_field( cmor_name, 'fLuc', (/id_ug/), &
       time, 'Net Carbon Mass Flux into Atmosphere due to Land Use Change', 'kg m-2 s-1', missing_value=-1.0, &
       standard_name='surface_net_upward_mass_flux_of_carbon_dioxide_expressed_as_carbon_due_to_emission_from_anthropogenic_land_use_change', &
       fill_missing=.TRUE.)
  id_fAnthDisturb = register_tiled_diag_field( cmor_name, 'fAnthDisturb', (/id_ug/), &
       time, 'carbon mass flux into atmosphere due to any human activity', 'kg m-2 s-1', missing_value=-1.0, &
       standard_name='surface_upward_mass_flux_of_carbon_dioxide_expressed_as_carbon_due_to_anthrogpogenic_emission', &
       fill_missing=.TRUE.)

  id_cw = register_tiled_diag_field ( cmor_name, 'cw',  &
       (/id_ug/), time, 'Total Canopy Water Storage', 'kg m-2', missing_value=-1.0, &
       standard_name='canopy_water_amount', fill_missing=.TRUE.)

  id_nVeg = register_tiled_diag_field ( cmor_name, 'nVeg', (/id_ug/), &
       time, 'Nitrogen Mass in Vegetation', 'kg m-2', missing_value=-1.0, &
       standard_name='vegetation_nitrogen_content', fill_missing=.TRUE.)
  id_nLeaf = register_tiled_diag_field ( cmor_name, 'nLeaf',  &
       (/id_ug/), time, 'Nitrogen Mass in Leaves', 'kg m-2', missing_value=-1.0, &
       standard_name='leaf_nitrogen_content', fill_missing=.TRUE.)
  id_nStem = register_tiled_diag_field ( cmor_name, 'nStem',  &
       (/id_ug/), time, 'Nitrogen Mass in Stem', 'kg m-2', missing_value=-1.0, &
       standard_name='stem_nitrogen_content', fill_missing=.TRUE.)
  id_nRoot = register_tiled_diag_field ( cmor_name, 'nRoot',  &
       (/id_ug/), time, 'Nitrogen Mass in Roots', 'kg m-2', missing_value=-1.0, &
       standard_name='root_nitrogen_content', fill_missing=.TRUE.)
  id_nOther = register_tiled_diag_field ( cmor_name, 'nOther',  &
       (/id_ug/), time, 'Nitrogen mass in vegetation components other than leaves, stem and root', &
       'kg m-2', missing_value=-1.0, &
       standard_name='other_vegegtation_components_nitrogen_content', fill_missing=.TRUE.)

end subroutine


! ============================================================================
! write restart file and release memory
subroutine vegn_end ()

  module_is_initialized =.FALSE.

  call vegn_fire_end()
  call vegn_harvesting_end()
  call static_vegn_end()
  call vegn_dynamics_end()
end subroutine vegn_end


! ============================================================================
subroutine save_vegn_restart(tile_dim_length,timestamp)
  integer, intent(in) :: tile_dim_length ! length of tile dim. in the output file
  character(*), intent(in) :: timestamp ! timestamp to add to the file name

  ! ---- local vars
  integer ::  i, j
  type(land_tile_enum_type) :: ce
  type(land_tile_type), pointer :: tile
  integer :: n_accum, nmn_acm

  character(267) :: filename
  type(land_restart_type) :: restart1, restart2 ! restart file i/o object
  character:: spnames(fm_field_name_len, nspecies) ! names of the species
  integer :: sp_dim,text_dim,spnames_id

  call error_mesg('vegn_end','writing NetCDF restart',NOTE)

  ! create output file, including internal structure necessary for tile output
  filename = trim(timestamp)//'vegn1.res.nc'
  call init_land_restart(restart1, filename, vegn_tile_exists, tile_dim_length)

  ! create compressed dimension for vegetation cohorts -- must be called even
  ! if restart has not been created, because it calls mpp_max and that should
  ! be called on all PEs to work
  call create_cohort_dimension(restart1)

  call add_cohort_data(restart1,'tv',cohort_tv_ptr,'vegetation temperature','degrees_K')
  call add_cohort_data(restart1,'wl',cohort_wl_ptr,'vegetation liquid water content','kg/individual')
  call add_cohort_data(restart1,'ws',cohort_ws_ptr,'vegetation solid water content','kg/individual')
  call save_land_restart(restart1)
  call free_land_restart(restart1)


  filename = trim(timestamp)//'vegn2.res.nc'
  call init_land_restart(restart2, filename, vegn_tile_exists, tile_dim_length)
  ! create compressed dimension for vegetation cohorts -- see note above
  call create_cohort_dimension(restart2)
  ! store table of species names
  call add_restart_axis(restart2,'nspecies',[(real(i),i=0,nspecies-1)],'Z')
  call add_restart_axis(restart2,'textlen',[(real(i),i=1,fm_field_name_len)],'Z')
  do i = 0, nspecies-1
     do j = 1,size(spnames,1)
        spnames(j,i+1) = ' '
     enddo
     do j = 1,min(len(spdata(i)%name),size(spnames,1))
        spnames(j,i+1) = spdata(i)%name(j:j)
     enddo
  enddo
  call add_text_data(restart2,'species_names','textlen','nspecies',spnames)

  ! store global variables
  ! find first tile and get n_accum and nmn_acm from it
  n_accum = 0; nmn_acm = 0
  ce = first_elmt(land_tile_map)
  do while (loop_over_tiles(ce,tile))
     if(associated(tile%vegn)) then
        n_accum = tile%vegn%n_accum
        nmn_acm = tile%vegn%nmn_acm
     endif
  enddo
  ! n_accum and nmn_acm are currently the same for all tiles; we only call mpp_max
  ! to handle the situation when there are no tiles in the current domain
  call mpp_max(n_accum); call mpp_max(nmn_acm)

  call add_scalar_data(restart2,'n_accum',n_accum,'number of accumulated steps')
  call add_scalar_data(restart2,'nmn_acm',nmn_acm,'number of accumulated months')

  call add_int_cohort_data(restart2,'species', cohort_species_ptr, 'vegetation species')
  call add_cohort_data(restart2,'hite', cohort_height_ptr, 'vegetation height','m')
  call add_cohort_data(restart2,'bl', cohort_bl_ptr, 'biomass of leaves','kg C/individual')
  call add_cohort_data(restart2,'blv', cohort_blv_ptr, 'biomass of virtual leaves (labile store)','kg C/individual')
  call add_cohort_data(restart2,'br', cohort_br_ptr, 'biomass of fine roots','kg C/individual')
  call add_cohort_data(restart2,'bsw', cohort_bsw_ptr, 'biomass of sapwood','kg C/individual')
  call add_cohort_data(restart2,'bwood', cohort_bwood_ptr, 'biomass of heartwood','kg C/individual')
  call add_cohort_data(restart2,'nsc', cohort_nsc_ptr, 'non-structural biomass','kg C/individual')
  call add_cohort_data(restart2,'bseed', cohort_bseed_ptr, 'biomass reserved for future progeny','kg C/individual')
  call add_cohort_data(restart2,'bsw_max', cohort_bsw_max_ptr, 'maximum biomass of sapwood','kg C/individual')
  call add_cohort_data(restart2,'bl_max', cohort_bl_max_ptr, 'maximum biomass of leaves','kg C/individual')
  call add_cohort_data(restart2,'br_max', cohort_br_max_ptr, 'maximum biomass of roots','kg C/individual')
  call add_cohort_data(restart2,'dbh', cohort_dbh_ptr, 'diameter at breast height','m')
  call add_cohort_data(restart2,'crownarea', cohort_crownarea_ptr, 'area of crown','m2/individual')
  call add_cohort_data(restart2,'BM_ys', cohort_BM_ys_ptr, 'bwood+bsw at the end of previous year','kg C/individual')
  call add_cohort_data(restart2,'DBH_ys', cohort_DBH_ys_ptr, 'DBH at the end of previous year','m')
  call add_cohort_data(restart2,'topyear', cohort_topyear_ptr, 'time spent in the top canopy layer','years')
  call add_cohort_data(restart2,'gdd', cohort_gdd_ptr, 'growing degree days','degC day')
  ! wolf restart data - psi, Kxa
  call add_cohort_data(restart2, 'psi_r', cohort_psi_r_ptr, 'psi root', 'm')
  call add_cohort_data(restart2, 'psi_x', cohort_psi_x_ptr, 'psi stem', 'm' )
  call add_cohort_data(restart2, 'psi_l', cohort_psi_l_ptr, 'psi leaf', 'm' )
  call add_cohort_data(restart2, 'Kxa',   cohort_Kxa_ptr, 'K stem per area', 'kg/m2s/(m/m)')
  call add_cohort_data(restart2, 'Kla',   cohort_Kla_ptr, 'K leaf per area', 'kg/m2s/m')

  call add_cohort_data(restart2,'bliving', cohort_bliving_ptr, 'total living biomass','kg C/individual')
  call add_cohort_data(restart2,'nindivs',cohort_nindivs_ptr, 'number of individuals', 'individuals/m2')
  call add_int_cohort_data(restart2,'layer',cohort_layer_ptr, 'canopy layer of cohort', 'unitless')
  call add_int_cohort_data(restart2,'firstlayer',cohort_firstlayer_ptr, &
                     'indicates whether cohort has ever been in the upper layer', 'unitless')
  call add_int_cohort_data(restart2,'status', cohort_status_ptr, 'leaf status')
  call add_cohort_data(restart2,'leaf_age',cohort_leaf_age_ptr, 'age of leaves since bud burst', 'days')
  call add_cohort_data(restart2,'cohort_age',cohort_age_ptr, 'age of cohort', 'years')
  call add_cohort_data(restart2,'npp_prev_day', cohort_npp_previous_day_ptr, 'previous day NPP','kg C/year')

  call add_cohort_data(restart2,'scav_myc_C_reservoir', cohort_scav_myc_C_reservoir_ptr, 'C in scavenger myc reservoir for growth','kg C/m2')
  call add_cohort_data(restart2,'scav_myc_N_reservoir', cohort_scav_myc_N_reservoir_ptr, 'N in scavenger myc reservoir for growth','kg C/m2')
  call add_cohort_data(restart2,'mine_myc_C_reservoir', cohort_mine_myc_C_reservoir_ptr, 'C in miner myc reservoir for growth','kg C/m2')
  call add_cohort_data(restart2,'mine_myc_N_reservoir', cohort_mine_myc_N_reservoir_ptr, 'N in miner myc reservoir for growth','kg C/m2')
  call add_cohort_data(restart2,'N_fixer_C_reservoir', cohort_N_fixer_C_reservoir_ptr, 'C in N fixer reservoir for growth','kg C/m2')
  call add_cohort_data(restart2,'N_fixer_N_reservoir', cohort_N_fixer_N_reservoir_ptr, 'N in N fixer reservoir for growth','kg C/m2')
  call add_cohort_data(restart2,'myc_scavenger_biomass_C', cohort_myc_scavenger_biomass_C_ptr, 'scavenger mycorrhizal biomass C associated with individual','kg C/m2')
  call add_cohort_data(restart2,'myc_scavenger_biomass_N', cohort_myc_scavenger_biomass_N_ptr, 'scavenger mycorrhizal biomass N associated with individual','kg N/m2')
  call add_cohort_data(restart2,'N_fixer_biomass_C', cohort_N_fixer_biomass_C_ptr, 'symbiotic N fixer biomass C associated with individual','kg C/m2')
  call add_cohort_data(restart2,'N_fixer_biomass_N', cohort_N_fixer_biomass_N_ptr, 'symbiotic N fixer biomass N associated with individual','kg N/m2')
  call add_cohort_data(restart2,'myc_miner_biomass_C', cohort_myc_miner_biomass_C_ptr, 'miner mycorrhizal biomass C associated with individual','kg C/m2')
  call add_cohort_data(restart2,'myc_miner_biomass_N', cohort_myc_miner_biomass_N_ptr, 'miner mycorrhizal biomass N associated with individual','kg N/m2')
  call add_cohort_data(restart2,'cohort_stored_N', cohort_stored_N_ptr, 'stored N pool of individual','kg N/m2')
  call add_cohort_data(restart2,'cohort_wood_N', cohort_wood_N_ptr, 'wood N pool of individual','kg N/m2')
  call add_cohort_data(restart2,'cohort_sapwood_N', cohort_sapwood_N_ptr, 'sapwood N pool of individual','kg N/m2')
  call add_cohort_data(restart2,'cohort_leaf_N', cohort_leaf_N_ptr, 'leaf N pool of individual','kg N/m2')
  call add_cohort_data(restart2,'cohort_seed_N', cohort_seed_N_ptr, 'seed N pool of individual','kg N/m2')
  call add_cohort_data(restart2,'cohort_root_N', cohort_root_N_ptr, 'root N pool of individual','kg N/m2')
  call add_cohort_data(restart2,'cohort_total_N', cohort_total_N_ptr, 'total N pool of individual','kg N/m2')
  call add_cohort_data(restart2,'nitrogen_stress', cohort_nitrogen_stress_ptr, 'total N pool of individual','kg N/m2')
  call add_cohort_data(restart2,'myc_scav_marginal_gain_smoothed', cohort_myc_scav_marginal_gain_smoothed_ptr, 'smoothed marginal gain of scavenging','gN/gC')
  call add_cohort_data(restart2,'myc_mine_marginal_gain_smoothed', cohort_myc_mine_marginal_gain_smoothed_ptr, 'smoothed marginal gain of mining','gN/gC')
  call add_cohort_data(restart2,'N_fix_marginal_gain_smoothed', cohort_N_fix_marginal_gain_smoothed_ptr, 'smoothed marginal gain of N fixation','gN/gC')
  call add_cohort_data(restart2,'rhiz_exud_marginal_gain_smoothed', cohort_rhiz_exud_marginal_gain_smoothed_ptr, 'smoothed marginal gain of root exudation','gN/gC')
  call add_cohort_data(restart2,'max_monthly_scav_alloc', cohort_max_monthly_scav_alloc_ptr, 'smoothed allocation to scavenging','kgC/m2/year')
  call add_cohort_data(restart2,'max_monthly_mine_alloc', cohort_max_monthly_mine_alloc_ptr, 'smoothed allocation to mining','kgC/m2/year')
  call add_cohort_data(restart2,'max_monthly_Nfix_alloc', cohort_max_monthly_Nfix_alloc_ptr, 'monthly maximum allocation to N fixation','kgC/m2/year')
  call add_cohort_data(restart2,'scav_alloc_accum', cohort_scav_alloc_accum_ptr, 'scavenging allocation accumulator','kgC/m2')
  call add_cohort_data(restart2,'mine_alloc_accum', cohort_mine_alloc_accum_ptr, 'mining allocation accumulator','kgC/m2')
  call add_cohort_data(restart2,'Nfix_alloc_accum', cohort_Nfix_alloc_accum_ptr, 'N fixation allocation accumulator','kgC/m2')
  call add_cohort_data(restart2,'max_scav_allocation', cohort_max_scav_allocation_ptr, 'max allowed allocation to scavenging','kgC/m2/year')
  call add_cohort_data(restart2,'max_mine_allocation', cohort_max_mine_allocation_ptr, 'max allowed allocation to mining','kgC/m2/year')
  call add_cohort_data(restart2,'max_Nfix_allocation', cohort_max_Nfix_allocation_ptr, 'max allowed allocation to N fixation','kgC/m2/year')
  call add_cohort_data(restart2,'N_stress_smoothed', cohort_nitrogen_stress_smoothed_ptr, 'Smoothed N stress','Dimensionless')

  call add_cohort_data(restart2,'growth_prev_day', cohort_growth_previous_day_ptr, 'pool of growth respiration','kg C')
  call add_cohort_data(restart2,'growth_prev_day_tmp', cohort_growth_previous_day_tmp_ptr, 'rate of growth respiration release to atmos','kg C/year')
  call add_cohort_data(restart2,'brsw', cohort_brsw_ptr, 'biomass of branches (only sapwood)','kg C/individual')

  do i = 0,nspecies-1
     call add_tile_data(restart2,'drop_seed_C_'//trim(spdata(i)%name),vegn_drop_seed_C_ptr,i,&
                        'seed carbon dropped by dying plants', 'kgC/m2')
     call add_tile_data(restart2,'drop_seed_N_'//trim(spdata(i)%name),vegn_drop_seed_C_ptr,i,&
                        'seed nirogen dropped by dying plants', 'kgC/m2')
  enddo

  call add_int_tile_data(restart2,'landuse',vegn_landuse_ptr,'vegetation land use type')
  call add_tile_data(restart2,'age',vegn_age_ptr,'vegetation age', 'yr')

  ! write carbon pools and rates
  call add_tile_data(restart2,'fsc_pool_ag',vegn_fsc_pool_ag_ptr,'intermediate pool for AG fast soil carbon input', 'kg C/m2')
  call add_tile_data(restart2,'fsc_rate_ag',vegn_fsc_rate_ag_ptr,'conversion rate of AG fsc_pool to fast soil carbon', 'kg C/(m2 yr)')
  call add_tile_data(restart2,'ssc_pool_ag',vegn_ssc_pool_ag_ptr,'intermediate pool for AG slow soil carbon input', 'kg C/m2')
  call add_tile_data(restart2,'ssc_rate_ag',vegn_ssc_rate_ag_ptr,'conversion rate of AG ssc_pool to slow soil carbon', 'kg C/(m2 yr)')
  call add_tile_data(restart2,'fsc_pool_bg',vegn_fsc_pool_bg_ptr,'intermediate pool for BG fast soil carbon input', 'kg C/m2')
  call add_tile_data(restart2,'fsc_rate_bg',vegn_fsc_rate_bg_ptr,'conversion rate of BG fsc_pool to fast soil carbon', 'kg C/(m2 yr)')
  call add_tile_data(restart2,'ssc_pool_bg',vegn_ssc_pool_bg_ptr,'intermediate pool for BG slow soil carbon input', 'kg C/m2')
  call add_tile_data(restart2,'ssc_rate_bg',vegn_ssc_rate_bg_ptr,'conversion rate of BG ssc_pool to slow soil carbon', 'kg C/(m2 yr)')

  do j = 1,N_LITTER_POOLS
     do i = 1,N_C_TYPES-1 ! "-1" excludes deadmic from restarts
        call add_tile_data(restart2,trim(l_shortname(j))//'litter_buffer_'//trim(c_shortname(i)),vegn_litter_buff_C_ptr, i, j, &
            'intermediate pool for '//trim(c_longname(i))//' '//trim(l_longname(j))//' litter carbon input', 'kg C/m2')
        call add_tile_data(restart2,trim(l_shortname(j))//'litter_buffer_rate_'//trim(c_shortname(i)),vegn_litter_rate_C_ptr, i, j, &
            'conversion rate of '//trim(c_longname(i))//' '//trim(l_longname(j))//' litter to litter carbon pool', 'kg C/(m2 yr)')
     enddo
  enddo

  if (soil_carbon_option==SOILC_CORPSE_N) then
     call add_tile_data(restart2,'fsn_pool_bg',vegn_fsn_pool_bg_ptr,'intermediate pool for BG fast soil nitrogen input', 'kg N/m2')
     call add_tile_data(restart2,'fsn_rate_bg',vegn_fsn_rate_bg_ptr,'conversion rate of BG fsn_pool to fast soil nitrogen', 'kg N/(m2 yr)')
     call add_tile_data(restart2,'ssn_pool_bg',vegn_ssn_pool_bg_ptr,'intermediate pool for BG slow soil nitrogen input', 'kg N/m2')
     call add_tile_data(restart2,'ssn_rate_bg',vegn_ssn_rate_bg_ptr,'conversion rate of BG ssn_pool to slow soil nitrogen', 'kg N/(m2 yr)')
     do j = 1,N_LITTER_POOLS
        do i = 1,N_C_TYPES-1 ! "-1" excludes deadmic from restarts
           call add_tile_data(restart2,trim(l_shortname(j))//'litter_buff_N_'//trim(c_shortname(i)),vegn_litter_buff_N_ptr, i, j, &
               'intermediate pool for '//trim(c_longname(i))//' '//trim(l_longname(j))//' litter nitrogen input', 'kg N/m2')
           call add_tile_data(restart2,trim(l_shortname(j))//'litter_rate_N_'//trim(c_shortname(i)),vegn_litter_rate_N_ptr, i, j, &
               'conversion rate of '//trim(c_longname(i))//' '//trim(l_longname(j))//' litter to litter nitrogen pool', 'kg N/(m2 yr)')
        enddo
     enddo
  endif

  ! water and heat buffers
  call add_tile_data(restart2,'drop_wl',vegn_drop_wl_ptr,'amount of liquid water dropped by dead trees, etc.', 'kg/m2')
  call add_tile_data(restart2,'drop_ws',vegn_drop_ws_ptr,'amount of solid dropped by dead trees, etc.', 'kg/m2')
  call add_tile_data(restart2,'drop_hl',vegn_drop_hl_ptr,'heat of liquid water dropped by dead trees, etc.', 'J/m2')
  call add_tile_data(restart2,'drop_hs',vegn_drop_hs_ptr,'heat of solid dropped by dead trees, etc.', 'J/m2')

  ! monthly-mean values
  call add_tile_data(restart2,'tc_av', vegn_tc_av_ptr,'average canopy air temperature','degK')
  call add_tile_data(restart2,'theta_av_phen', vegn_theta_av_phen_ptr,'average soil moisture for phenology')
  call add_tile_data(restart2,'theta_av_fire', vegn_theta_av_fire_ptr,'average soil moisture for fire')
  call add_tile_data(restart2,'psist_av', vegn_psist_av_ptr,'average soil-water-stress index')
  call add_tile_data(restart2,'tsoil_av', vegn_tsoil_av_ptr,'average bulk soil temperature for soil carbon','degK')
  call add_tile_data(restart2,'precip_av', vegn_precip_av_ptr,'average total precipitation','kg/(m2 s)')
  call add_tile_data(restart2,'lambda', vegn_lambda_ptr,'dryness parameter')
  call add_tile_data(restart2,'fuel', vegn_fuel_ptr,'fuel density','kg C/m2')
  ! annual-mean values
  call add_tile_data(restart2,'t_ann', vegn_t_ann_ptr,'average annual canopy air temperature','degK')
  call add_tile_data(restart2,'t_cold', vegn_t_cold_ptr,'average canopy air temperature of coldest month','degK')
  call add_tile_data(restart2,'p_ann', vegn_p_ann_ptr,'average annual precipitation','kg/(m2 s)')
  call add_tile_data(restart2,'ncm', vegn_ncm_ptr,'number of cold months')
  call add_tile_data(restart2,'treeline_T_accum', vegn_treeline_T_accum_ptr,'accumulated temperature for tree line calculations','degC')
  call add_tile_data(restart2,'treeline_N_accum', vegn_treeline_N_accum_ptr,'number of temperature samples for tree line calculations')
  ! accumulated values for annual averaging
  call add_tile_data(restart2,'t_ann_acm', vegn_t_ann_acm_ptr,'accumulated annual canopy air temperature','degK')
  call add_tile_data(restart2,'t_cold_acm', vegn_t_cold_acm_ptr,'accumulated temperature of coldest month','degK')
  call add_tile_data(restart2,'p_ann_acm', vegn_p_ann_acm_ptr,'accumulated precipitation','kg/(m2 s)')
  call add_tile_data(restart2,'ncm_acm', vegn_ncm_acm_ptr,'accumulated number of cold months')
  ! accumulated and averaged values for PPA phenology
  call add_tile_data(restart2,'tc_pheno', vegn_tc_pheno_ptr,'smoothed temperature for phenology','degK')
  ! burned carbon pool and rate
  call add_tile_data(restart2,'csmoke_pool',vegn_csmoke_pool_ptr,'carbon lost through fires', 'kg C/m2')
  call add_tile_data(restart2,'csmoke_rate',vegn_csmoke_rate_ptr,'rate of release of carbon lost through fires to the atmosphere', 'kg C/(m2 yr)')
  call add_tile_data(restart2,'nsmoke_pool',vegn_nsmoke_pool_ptr,'nitrogen lost through fires', 'kg N/m2')

  ! harvesting pools and rates
  do i = 1, N_HARV_POOLS
     call add_tile_data(restart2, trim(HARV_POOL_NAMES(i))//'_harv_pool', &
          vegn_harv_pool_C_ptr, i, 'harvested carbon','kg C/m2')
     call add_tile_data(restart2, trim(HARV_POOL_NAMES(i))//'_harv_pool_nitrogen', &
          vegn_harv_pool_N_ptr, i, 'harvested nitrogen','kg N/m2')
     call add_tile_data(restart2, trim(HARV_POOL_NAMES(i))//'_harv_rate', &
          vegn_harv_rate_C_ptr, i, 'rate of release of harvested carbon to the atmosphere','kg C/(m2 yr)')
  enddo

  call save_land_restart(restart2)
  call free_land_restart(restart2)
end subroutine save_vegn_restart


! ============================================================================
! given vegetation state and snow depth, calculate integral diffusion-related
! properties
subroutine vegn_diffusion (vegn, snow_depth, vegn_cover, vegn_height, vegn_lai, vegn_sai)
  type(vegn_tile_type), intent(inout) :: vegn
  real, intent(in) :: snow_depth
  real, intent(out) :: &
       vegn_cover, vegn_height, vegn_lai, vegn_sai

  real :: gaps,      & ! fraction of gaps in the canopy, =1-cover
          layer_gaps   ! fraction of gaps in the canopy in a single layer, accumulator value
  integer :: i, current_layer

  associate(cc=>vegn%cohorts) ! F2003
  ! calculate integral parameters of vegetation
  gaps = 1.0; vegn_lai = 0; vegn_sai = 0
  current_layer = cc(1)%layer ; layer_gaps = 1.0
  vegn_height = 0.0
  do i = 1,vegn%n_cohorts
    call vegn_data_cover(cc(i), snow_depth)
    vegn_lai = vegn_lai + cc(i)%lai*cc(i)%layerfrac
    vegn_sai = vegn_sai + cc(i)%sai*cc(i)%layerfrac
    ! calculate total cover
    if (cc(i)%layer/=current_layer) then
       gaps = gaps*layer_gaps; layer_gaps = 1.0; current_layer = cc(i)%layer
    endif
    layer_gaps = layer_gaps-cc(i)%layerfrac*cc(i)%cover
    vegn_height = max(vegn_height,cc(i)%height)
  enddo
  gaps = gaps*layer_gaps ! take the last layer into account
  vegn_cover  = 1 - gaps
  end associate ! F2003

end subroutine vegn_diffusion


! ============================================================================
subroutine vegn_step_1 ( vegn, soil, diag, &
        p_surf, ustar, drag_q, &
        SWdn, RSv, precip_l, precip_s, &
        land_d, land_z0s, land_z0m, grnd_z0s, &
        cana_T, cana_q, cana_co2_mol, &
        ! output
        con_g_h, con_g_v, & ! aerodynamic conductance between canopy air and ground, for heat and vapor flux
        con_v_v, & ! aerodynamic conductance between canopy air and canopy
        stomatal_cond, & ! integral stomatal conductance of cohort canopy
        vegn_T, vegn_Wl, vegn_Ws, & ! temperature, water and snow mass of the canopy
        vegn_ifrac,               & ! intercepted fraction of liquid and frozen precipitation
        vegn_lai,                 & ! leaf area index
        drip_l, drip_s,           & ! water and snow drip rate from precipitation, kg/(m2 s)
        prec_l, prec_s,           & ! liquid and solid precip on top of each cohort, kg/(m2 s)
        prec_g_l, prec_g_s,       & ! liquid and solid precipitation reaching ground, kg/(m2 s)
        vegn_hcap,                & ! vegetation heat capacity
        Hv0,   DHvDTv,   DHvDTc,                      & ! sens heat flux
        Et0,   DEtDTv,   DEtDqc,   DEtDwl,   DEtDwf,  & ! transpiration
        Eli0,  DEliDTv,  DEliDqc,  DEliDwl,  DEliDwf, & ! evaporation of intercepted water
        Efi0,  DEfiDTv,  DEfiDqc,  DEfiDwl,  DEfiDwf, & ! sublimation of intercepted snow
        soil_uptake_T ) ! average temperature of water lost through transpiration, degK
  type(vegn_tile_type), intent(inout) :: vegn ! vegetation data
  type(soil_tile_type), intent(in)    :: soil ! soil data
  type(diag_buff_type), intent(inout) :: diag ! diagnostic buffer
  real, intent(in) :: &
       p_surf,    & ! surface pressure, N/m2
       ustar,     & ! friction velocity, m/s
       drag_q,    & ! bulk drag coefficient for specific humidity
       SWdn(:,:), & ! downward SW radiation on top of each cohort, W/m2 (NCOHORTS,NBANDS)
       RSv (:,:), & ! net SW radiation balance of the canopy, W/m2, (NCOHORTS,NBANDS)
       precip_l, precip_s, & ! liquid and solid precipitation rates, kg/(m2 s)
       land_d, land_z0s, land_z0m, & ! land displacement height and roughness, m
       grnd_z0s, & ! roughness of ground surface (including snow effect)
       cana_T,    & ! temperature of canopy air, deg K
       cana_q,    & ! specific humidity of canopy air, kg/kg
       cana_co2_mol ! co2 mixing ratio in the canopy air, mol CO2/mol dry air
  ! output -- coefficients of linearized expressions for fluxes
  real, intent(out), dimension(:) ::   &
       con_v_v, & ! aerodyn. conductance between canopy and CAS, for vapor (and tracers)
       stomatal_cond, & ! integral stomatal conductance of cohort canopy
       vegn_T, vegn_Wl, vegn_Ws,& ! temperature, water and snow mass of the canopy
       vegn_ifrac, & ! intercepted fraction of liquid and frozen precipitation
       vegn_lai, & ! vegetation leaf area index
       drip_l, drip_s, & ! water and snow drip rate from precipitation, kg/(m2 s)
       prec_l, prec_s, & ! liquid and solid precip on top of cohort, kg/(m2 s)
       vegn_hcap, & ! total vegetation heat capacity, including intercepted water and snow
       Hv0,   DHvDTv,   DHvDTc, & ! sens heat flux
       Et0,   DEtDTv,   DEtDqc,   DEtDwl,   DEtDwf,  & ! transpiration
       Eli0,  DEliDTv,  DEliDqc,  DEliDwl,  DEliDwf, & ! evaporation of intercepted water
       Efi0,  DEfiDTv,  DEfiDqc,  DEfiDwl,  DEfiDwf, & ! sublimation of intercepted snow
       soil_uptake_T
  real, intent(out) :: &
       prec_g_l, prec_g_s, & ! liquid and solid precipitation reaching ground, kg/(m2 s)
       con_g_h, con_g_v  ! aerodynamic conductance between ground and canopy air

  ! ---- local constants
  real, parameter :: min_reported_psi = -HUGE(1.0_4) ! limit water potential that we send to diagnostics
  ! "_4" in the above indicates single-precision (real*4) literal
  ! ---- local vars
  real :: &
       ft,DftDwl,DftDwf, & ! fraction of canopy not covered by intercepted water/snow, and its
                           ! derivatives w.r.t. intercepted water masses
       fw,DfwDwl,DfwDwf, & ! fraction of canopy covered by intercepted water, and its
                           ! derivatives w.r.t. intercepted water masses
       fs,DfsDwl,DfsDwf, & ! fraction of canopy covered by intercepted snow, and its
                           ! derivatives w.r.t. intercepted water masses
       precip_above_l, precip_above_s, & ! liquid and solid precip on top of the current layer, kg/(m2 s)
       precip_under_l, precip_under_s, & ! liquid and solid precip under the current layer, kg/(m2 s)
       total_stomatal_cond, & ! sum of cohort stomatal conductance values, for diagnostics only
       total_cond, &! overall conductance from inside stomata to canopy air
       qvsat,     & ! sat. specific humidity at the leaf T
       DqvsatDTv, & ! derivative of qvsat w.r.t. leaf T
       rho,       & ! density of canopy air
       gaps,      & ! fraction of gaps in the canopy, used to calculate cover
       layer_gaps,& ! fraction of gaps in the canopy in a single layer, accumulator value
       phot_co2     ! co2 mixing ratio for photosynthesis, mol CO2/mol dry air
  real, dimension(vegn%n_cohorts) :: &
       con_v_h, & ! aerodyn. conductance between canopy and CAS, for heat and vapor
       soil_beta, & ! relative water availability
       soil_water_supply, & ! max rate of water supply to the roots, kg/(indiv s)
       evap_demand, & ! plant evaporative demand, kg/(indiv s)
       RHi, &       ! relative humidity inside the leaf, at the point of vaporization
       lai_kok, &   ! LAI above 40 umoles of light
       An_newleaf   ! derivative of An w.r.t. LAI, for diagnostics only

  type(vegn_cohort_type), pointer :: cc(:)
  integer :: i, current_layer, band, N
  real :: indiv2area ! conversion factor from X per indiv. to X per unit cohort area
  real :: area2indiv ! reciprocal of the indiv2area conversion factor
  real :: rav_lit    ! litter resistance to water vapor
  real :: litter_fast_C, litter_slow_C, litter_deadmic_C ! for rav_lit calculations

  cc => vegn%cohorts(1:vegn%n_cohorts) ! note that the size of cc is always N

  if(is_watch_point()) then
     write(*,*)'#### vegn_step_1 input ####'
     __DEBUG3__(p_surf, ustar, drag_q)
     do band = 1,NBANDS
        __DEBUG1__(SWdn(:,band))
     enddo
     do band = 1,NBANDS
        __DEBUG1__(RSv(:,band))
     enddo
     __DEBUG2__(precip_l, precip_s)
     __DEBUG4__(land_d, land_z0s, land_z0m, grnd_z0s)
     __DEBUG3__(cana_T, cana_q, cana_co2_mol)
     write(*,*)'#### end of vegn_step_1 input ####'
     __DEBUG1__(cc%layer)
     __DEBUG1__(cc%species)
     __DEBUG1__(cc%nindivs)
     __DEBUG1__(cc%crownarea)
     __DEBUG1__(cc%layerfrac)
     __DEBUG1__(cc%height)
     __DEBUG1__(cc%zbot)
     __DEBUG1__(cc%lai)
     __DEBUG1__(cc%sai)
     __DEBUG1__(cc%cover)
     __DEBUG1__(cc%leaf_size)
     __DEBUG1__(cc%Tv)
     __DEBUG1__(cc%Wl)
     __DEBUG1__(cc%Wl_max)
  endif
  ! TODO: check array sizes

  ! TODO: verify cover calculations

  gaps = 1.0 ; current_layer = cc(1)%layer ; layer_gaps = 1.0
  ! check the range of input temperature
  call check_temp_range(cc(1:vegn%n_cohorts)%Tv, 'vegn_step_1','Tv of cohort')
  do i = 1,vegn%n_cohorts
     ! calculate the fractions of intercepted precipitation
     vegn_ifrac(i) = cc(i)%cover
     ! get the lai
     vegn_lai(i)   = cc(i)%lai

     ! calculate total cover
     if (cc(i)%layer/=current_layer) then
        gaps = gaps*layer_gaps; layer_gaps = 1.0; current_layer = cc(i)%layer
     endif
     layer_gaps = layer_gaps-cc(i)%layerfrac*cc(i)%cover
  enddo
  gaps = gaps*layer_gaps ! take the last layer into account

  ! calculate the aerodynamic conductance coefficients
  call cana_turbulence(ustar, 1-gaps, &
     cc(:)%layerfrac, cc(:)%height, cc(:)%zbot, cc(:)%lai, cc(:)%sai, cc(:)%leaf_size, &
     land_d, land_z0m, land_z0s, grnd_z0s, &
     ! output:
     con_v_h, con_v_v, con_g_h, con_g_v)

  ! take into account additional resistance of litter to the water vapor flux.
  ! not a good parameterization, but just using for sensitivity analyses now.
  ! ignores differing biomass and litter turnover rates.
  call get_soil_litter_C(soil, litter_fast_C, litter_slow_C, litter_deadmic_C)
  rav_lit = rav_lit_0 + rav_lit_vi * (vegn_tile_LAI(vegn)+vegn_tile_SAI(vegn)) &
                      + rav_lit_fsc * litter_fast_C &
                      + rav_lit_ssc * litter_slow_C &
                      + rav_lit_deadmic * litter_deadmic_C &
                      + rav_lit_bwood * sum(cc(:)%bwood*cc(:)%nindivs)
  con_g_v = con_g_v/(1.0+rav_lit*con_g_v)

  if(is_watch_point()) then
     __DEBUG1__(con_v_h)
     __DEBUG1__(con_v_v)
     __DEBUG1__(con_g_h)
     __DEBUG1__(con_g_v)
  endif

  call soil_data_beta ( soil, vegn, soil_beta, soil_water_supply, soil_uptake_T )
  if(is_watch_point()) then
     __DEBUG1__(soil_beta)
     __DEBUG1__(soil_water_supply)
  endif

  ! calculate the vegetation photosynthesis and associated stomatal conductance
  if (vegn_phot_co2_option == VEGN_PHOT_CO2_INTERACTIVE) then
     phot_co2 = cana_co2_mol
  else
     phot_co2 = co2_for_photosynthesis
  endif

  total_stomatal_cond = 0
  precip_above_l = precip_l ; precip_under_l = precip_l
  precip_above_s = precip_s ; precip_under_s = precip_s
  current_layer = cc(1)%layer
  do i = 1, vegn%n_cohorts
     call vegn_photosynthesis (soil, vegn, cc(i), &
        SWdn(i,BAND_VIS), RSv(i,BAND_VIS), cana_T, cana_q, phot_co2, p_surf, drag_q, &
        soil_beta(i), soil_water_supply(i), con_v_v(i), &
        ! output
        evap_demand(i), stomatal_cond(i), RHi(i), lai_kok(i), An_newleaf(i))

     ! accumulate total value of stomatal conductance for diagnostics.
     ! stomatal_cond is per unit area of cohort (multiplied by LAI in the
     ! vegn_photosynthesis), so the total_stomatal_cond is per unit area
     ! of tile
     total_stomatal_cond = total_stomatal_cond+stomatal_cond(i)*cc(i)%layerfrac

     ! get the amount of intercepted water and snow per unir area of cohort
     if (cc(i)%nindivs>0) then
        indiv2area = cc(i)%nindivs/cc(i)%layerfrac
        area2indiv = 1/indiv2area
     else
        ! protection from the cohorts whose nindivs and layerfrac are both zero
        ! can happen due to mortality/starvation etc
        indiv2area = 0.0
        area2indiv = 0.0
     endif
     vegn_Wl(i) = cc(i)%Wl*indiv2area
     vegn_Ws(i) = cc(i)%Ws*indiv2area

     call get_vegn_wet_frac ( cc(i), fw, DfwDwl, DfwDwf, fs, DfsDwl, DfsDwf )
     ! derivatives must be renormalized, because the units of canopy water and
     ! snow used in calculations are kg/indiv, and the equations are written for
     ! units of kg/(m2 of stretched cohort)
     DfwDwl = DfwDwl*area2indiv ; DfwDwf = DfwDwf*area2indiv
     DfsDwl = DfsDwl*area2indiv ; DfsDwf = DfsDwf*area2indiv
     ! transpiring fraction and its derivatives
     ft     = 1 - fw - fs
     DftDwl = - DfwDwl - DfsDwl
     DftDwf = - DfwDwf - DfsDwf
     call qscomp(cc(i)%Tv, p_surf, qvsat, DqvsatDTv)

     rho = p_surf/(rdgas*cana_T *(1+d608*cana_q))

     ! get the vegetation temperature
     vegn_T(i)  =  cc(i)%Tv

     if(current_layer/=cc(i)%layer) then
        ! set the precipitation on top of current
        precip_above_l = precip_under_l
        precip_above_s = precip_under_s
        current_layer = cc(i)%layer
     endif
     ! accumulate precipitation under current layer: it is equal to precipitation
     ! above minus the intercepted rainfall
     precip_under_l = precip_under_l - precip_above_l*vegn_ifrac(i)*cc(i)%layerfrac
     precip_under_s = precip_under_s - precip_above_s*vegn_ifrac(i)*cc(i)%layerfrac

     ! set the precipitation on top of the canopy
     prec_l(i)  = precip_above_l; prec_s(i) = precip_above_s
     ! calculate the drip rates, kg/(s m2 of cohort)
     drip_l(i)  = max(vegn_Wl(i),0.0)/tau_drip_l
     drip_s(i)  = max(vegn_Ws(i),0.0)/tau_drip_s
     ! correct the drip rates so that the amount of water and snow accumulated over time step
     ! is no larger then the canopy water-holding capacity
     drip_l(i) = max((vegn_Wl(i)+prec_l(i)*delta_time*vegn_ifrac(i)-cc(i)%Wl_max*indiv2area)/delta_time,drip_l(i))
     drip_s(i) = max((vegn_Ws(i)+prec_s(i)*delta_time*vegn_ifrac(i)-cc(i)%Ws_max*indiv2area)/delta_time,drip_s(i))

     ! calculate the total heat capacity per unit area of cohort
     vegn_hcap(i) = (cc(i)%mcv_dry + clw*cc(i)%Wl + csw*cc(i)%Ws)*indiv2area
     ! calculate the coefficient of sensible heat flux linearization
     Hv0    (i) =  2*rho*cp_air*con_v_h(i)*(cc(i)%Tv - cana_T)
     DHvDTv (i) =  2*rho*cp_air*con_v_h(i)
     DHvDTc (i) = -2*rho*cp_air*con_v_h(i)
     ! calculate the coefficients of the transpiration linearization
     if(con_v_v(i)==0.and.stomatal_cond(i)==0) then
        total_cond = 0.0
     else
        total_cond = stomatal_cond(i)*con_v_v(i)/(stomatal_cond(i)+con_v_v(i))
     endif

     if(qvsat>cana_q*RHi(i))then
        ! Flux is directed FROM the surface: transpiration is possible

        ! prohibit transpiration if leaf temperature below some predefined minimum
        ! typically (268K, but check namelist)
        if(cc(i)%Tv < T_transp_min) total_cond = 0
        ! calculate the transpiration linearization coefficients
        Et0    (i) =  rho*total_cond*ft*(qvsat*RHi(i) - cana_q)
        DEtDTv (i) =  rho*total_cond*ft*DqvsatDTv*RHi(i)
        DEtDqc (i) = -rho*total_cond*ft
        DEtDwl (i) =  rho*total_cond*DftDwl*(qvsat*RHi(i) - cana_q)
        DEtDwf (i) =  rho*total_cond*DftDwf*(qvsat*RHi(i) - cana_q)
     else
        ! Flux is directed TOWARD the surface: no transpiration (assuming plants do not
        ! take water through stomata)
        Et0   (i)  = 0
        DEtDTv(i)  = 0; DEtDwl(i) = 0; DEtDwf(i) = 0;
        DEtDqc(i)  = 0
     endif

     if(qvsat>cana_q)then
        ! Flux is directed FROM the surface:  the evaporation of intercepted water
        ! depends on the fraction of wet/snow covered canopy.

        ! calculate the coefficients of the intercepted liquid evaporation linearization
        Eli0   (i) =  rho*con_v_v(i)*fw*(qvsat - cana_q)
        DEliDTv(i) =  rho*con_v_v(i)*fw*DqvsatDTv
        DEliDqc(i) = -rho*con_v_v(i)*fw
        DEliDwl(i) =  rho*con_v_v(i)*DfwDwl*(qvsat-cana_q)
        DEliDwf(i) =  rho*con_v_v(i)*DfwDwf*(qvsat-cana_q)
        ! calculate the coefficients of the intercepted snow evaporation linearization
        Efi0   (i) =  rho*con_v_v(i)*fs*(qvsat - cana_q)
        DEfiDTv(i) =  rho*con_v_v(i)*fs*DqvsatDTv
        DEfiDqc(i) = -rho*con_v_v(i)*fs
        DEfiDwl(i) =  rho*con_v_v(i)*DfsDwl*(qvsat-cana_q)
        DEfiDwf(i) =  rho*con_v_v(i)*DfsDwf*(qvsat-cana_q)
     else
        ! Flux is directed TOWARD the surface: condensation does not depend on the
        ! fraction of wet canopy -- dew formation occurs on the entire surface

        ! calculate dew or frost formation rates, depending on the temperature
        Eli0   (i) = 0; Efi0   (i) = 0
        DEliDTv(i) = 0; DEfiDTv(i) = 0
        DEliDqc(i) = 0; DEfiDqc(i) = 0
        DEliDwl(i) = 0; DEfiDwl(i) = 0
        DEliDwf(i) = 0; DEfiDwf(i) = 0
        if(vegn_T(i) >= tfreeze) then
           ! calculate the coefficients of the intercepted liquid condensation linearization
           Eli0   (i) =  rho*con_v_v(i)*(qvsat - cana_q)
           DEliDTv(i) =  rho*con_v_v(i)*DqvsatDTv
           DEliDqc(i) = -rho*con_v_v(i)
        else
           ! calculate the coefficients of the intercepted snow condensation linearization
           Efi0   (i) =  rho*con_v_v(i)*(qvsat - cana_q)
           DEfiDTv(i) =  rho*con_v_v(i)*DqvsatDTv
           DEfiDqc(i) = -rho*con_v_v(i)
        endif
        ! prohibit switching from condensation to evaporation if the water content
        ! is below certain (negative) threshold. Default min_Wl == min_Ws = -1 kg/m2
        if (vegn_Wl(i) < min_Wl) then
           Eli0(i) = 0 ; DEliDTv(i) = 0 ; DEliDqc(i) = 0 ; DEliDwl(i) = 0 ; DEliDwf(i) = 0
        endif
        if (vegn_Ws(i) < min_Ws) then
           Efi0(i) = 0 ; DEfiDTv(i) = 0 ; DEfiDqc(i) = 0 ; DEfiDwl(i) = 0 ; DEfiDwf(i) = 0
        endif

     endif
  enddo ! loop by cohorts
  ! assign values of precipitation reaching ground
  prec_g_l = precip_under_l; prec_g_s = precip_under_s

!  if (is_watch_point()) then
!     __DEBUG1__(cc%layer)
!     __DEBUG1__(cc%layerfrac)
!     __DEBUG1__(vegn_ifrac)
!     __DEBUG1__(prec_l)
!     __DEBUG1__(prec_s)
!     __DEBUG2__(prec_g_l, prec_g_s)
!  endif

  ! ---- diagnostic section
  call send_tile_data(id_stomatal, total_stomatal_cond, diag)
  ! An_op and An_cl is per unit area of leaf, so we average over the leaf area
  N = vegn%n_cohorts
  call send_cohort_data(id_an_op, diag, cc(:), cc(:)%An_op, weight=cc(:)%layerfrac*cc(:)%lai, op=OP_AVERAGE)
  call send_cohort_data(id_an_cl, diag, cc(:), cc(:)%An_cl, weight=cc(:)%layerfrac*cc(:)%lai, op=OP_AVERAGE)
  ! con_v_h and con_v_v are per unit area of cohort -- output is per unit tile area
  call send_tile_data(id_con_v_h, sum(con_v_h(:)*cc(:)%layerfrac), diag)
  call send_tile_data(id_con_v_v, sum(con_v_v(:)*cc(:)%layerfrac), diag)
  call send_tile_data(id_phot_co2, phot_co2, diag)
  call send_tile_data(id_soil_water_supply, sum(soil_water_supply(:)*cc(:)%nindivs), diag)
  call send_tile_data(id_evap_demand, sum(evap_demand(:)*cc(:)%nindivs), diag)

  ! plant hydraulics diagnostics
  call send_cohort_data(id_Kxi   , diag, cc(:), cc(:)%Kxi, weight=cc(:)%nindivs, op=OP_AVERAGE)
  call send_cohort_data(id_Kli   , diag, cc(:), cc(:)%Kli, weight=cc(:)%nindivs, op=OP_AVERAGE)
  ! TODO: perhaps use something else for averaging weight
  ! factor 1e-6 converts Pa to MPa
  call send_cohort_data(id_psi_r , diag, cc(:), max(cc(:)%psi_r,min_reported_psi)*1e-6, weight=cc(:)%nindivs, op=OP_AVERAGE)
  call send_cohort_data(id_psi_x , diag, cc(:), max(cc(:)%psi_x,min_reported_psi)*1e-6, weight=cc(:)%nindivs, op=OP_AVERAGE)
  call send_cohort_data(id_psi_l , diag, cc(:), max(cc(:)%psi_l,min_reported_psi)*1e-6, weight=cc(:)%nindivs, op=OP_AVERAGE)
  ! limiting the water potentials that we send to diagnostics is a kludge. From time to
  ! time the solution of hydraulic equation creates very low water potential (~ -1e+109)
  ! to close the stomata if there is no water in the soil. Those unphysical low potentials
  ! do not seem to create any problems in the model itself, but sending them to the
  ! diagnostics causes model to stop to with "out of bounds" check if the diagnostics
  ! debugging is enabled.
  call send_cohort_data(id_w_scale,diag, cc(:), cc(:)%w_scale,    weight=cc(:)%nindivs, op=OP_AVERAGE)
  call send_cohort_data(id_RHi,    diag, cc(:), RHi(:)*100,  weight=cc(:)%layerfrac*cc(:)%lai, op=OP_AVERAGE)
  ! Kok effect ppg 2017-11-03
  call send_cohort_data(id_lai_kok, diag, cc(:), lai_kok(:), weight=cc(:)%layerfrac, op=OP_SUM)

  call send_cohort_data(id_DanDlai, diag, cc(:), An_newleaf(:), weight=cc(:)%layerfrac, op=OP_SUM)
  call send_cohort_data(id_PAR_dn,  diag, cc(:), SWdn(:,BAND_VIS), weight=cc(:)%layerfrac, op=OP_SUM)
  call send_cohort_data(id_PAR_net, diag, cc(:), RSv(:,BAND_VIS), weight=cc(:)%layerfrac, op=OP_SUM)

end subroutine vegn_step_1


! ============================================================================
! Given the surface solution, substitute it back into the vegetation equations
! to determine new vegetation state.
subroutine vegn_step_2 ( vegn, diag, &
     delta_Tv, delta_wl, delta_wf, &
     total_vegn_melt, &
     total_vegn_ovfl_l,  total_vegn_ovfl_s,  & ! overflow of liquid and solid water from the canopy, kg/(m2 s)
     total_vegn_ovfl_Hl, total_vegn_ovfl_Hs  ) ! heat flux carried from canopy by overflow, W/(m2 s)

  ! ---- arguments
  type(vegn_tile_type) , intent(inout) :: vegn
  type(diag_buff_type) , intent(inout) :: diag
  real, intent(in), dimension(:) :: & ! per cohort
       delta_Tv, & ! change in vegetation temperature, degK
       delta_wl, & ! change in intercepted liquid water mass, kg/m2 of (stretched) cohort
       delta_wf    ! change in intercepted frozen water mass, kg/m2 of (stretched) cohort
  real, intent(out) :: &
       total_vegn_melt, &
       total_vegn_ovfl_l,   total_vegn_ovfl_s,   & ! overflow of liquid and solid water from the canopy, kg/(m2 s) (per tike area)
       total_vegn_ovfl_Hl,  total_vegn_ovfl_Hs     ! heat flux from canopy due to overflow, W/m2 of tile

  ! ---- local variables
  real :: &
     cap0, melt_per_deg, &
     Wl, Ws, & ! positively defined amounts of water and snow on canopy
     vegn_melt, &
     vegn_ovfl_l,  vegn_ovfl_s,  & ! overflow of liquid and solid water from the cohort canopy, kg/(individual s)
     vegn_ovfl_Hl, vegn_ovfl_Hs    ! heat flux carried from cohort canopy by overflow, W/(m2 s)
  type(vegn_cohort_type), pointer :: cc
  integer :: i
  integer :: N ! shortcut for number of cohorts
  real :: indiv2area ! conversion factor from values per indiv. to values per unit cohort area
  real :: area2indiv ! reciprocal of the indiv2area conversion factor

  N = vegn%n_cohorts
  if (is_watch_point()) then
     write(*,*)'#### vegn_step_2 input ####'
     __DEBUG1__(vegn%cohorts(1:N)%Tv)
     __DEBUG1__(vegn%cohorts(1:N)%Wl)
     __DEBUG1__(vegn%cohorts(1:N)%Ws)
     __DEBUG1__(delta_Tv)
     __DEBUG1__(delta_wl)
     __DEBUG1__(delta_wf)
  endif

  total_vegn_melt = 0
  total_vegn_ovfl_l  = 0;  total_vegn_ovfl_s  = 0
  total_vegn_ovfl_Hl = 0;  total_vegn_ovfl_Hs = 0
  do i = 1,vegn%n_cohorts
     ! get the pointer to the current cohort
     cc => vegn%cohorts(i)


     ! update vegetation state
     if (cc%nindivs > 0) then
        indiv2area = cc%nindivs/cc%layerfrac
        area2indiv = 1/indiv2area
     else
        indiv2area = 0.0
        area2indiv = 0.0
     endif
     cc%Tv = cc%Tv + delta_Tv(i)
     cc%Wl = cc%Wl + delta_wl(i)*area2indiv
     cc%Ws = cc%Ws + delta_wf(i)*area2indiv

  ! ---- update for evaporation and interception -----------------------------
     cap0 = cc%mcv_dry + clw*cc%Wl + csw*cc%Ws ! J/(K individual)

     ! melt on the vegetation should probably be prohibited altogether, since
     ! the amount of melt or freeze calculated this way is severely underestimated
     ! (depending on the overall vegetation heat capacity) which leads to extended
     ! periods when the canopy temperature is fixed at freezing point.
     if (lm2) then
        vegn_melt = 0
     else
        ! ---- freeze/melt of intercepted water
        ! heat capacity of leaf + intercepted water/snow _can_ go below zero if the
        ! total water content goes below zero as a result of implicit time step.
        ! If it does, we just prohibit melt, setting it to zero.
        if(cap0 > 0)then
           melt_per_deg = cap0 / hlf
              if (cc%Ws>0 .and. cc%Tv>tfreeze) then
                 vegn_melt =  min(cc%Ws, (cc%Tv-tfreeze)*melt_per_deg)
              else if (cc%Wl>0 .and. cc%Tv<tfreeze) then
                 vegn_melt = -min(cc%Wl, (tfreeze-cc%Tv)*melt_per_deg)
           else
              vegn_melt = 0
           endif
              cc%Ws = cc%Ws - vegn_melt
              cc%Wl = cc%Wl + vegn_melt
           if (vegn_melt/=0) &
                   cc%Tv = tfreeze + (cap0*(cc%Tv-tfreeze) - hlf*vegn_melt) &
                / ( cap0 + (clw-csw)*vegn_melt )
           vegn_melt = vegn_melt / delta_time
        else
           vegn_melt = 0
        endif
     endif
     ! vegn_melt is per individual here

     if(is_watch_point()) then
        write (*,*)'#### vegn_step_2 #### 1'
           __DEBUG4__(i,cc%Tv, cc%Wl, cc%Ws)
           __DEBUG4__(cc%nindivs, cc%crownarea, cc%layerfrac, indiv2area)
           __DEBUG1__(vegn_melt)
     endif

     ! ---- update for overflow -------------------------------------------------
     Wl = max(cc%Wl,0.0); Ws = max(cc%Ws,0.0)
     vegn_ovfl_l = max (0.,Wl-cc%Wl_max)/delta_time
     vegn_ovfl_s = max (0.,Ws-cc%Ws_max)/delta_time
     vegn_ovfl_Hl = clw*vegn_ovfl_l*(cc%Tv-tfreeze)
     vegn_ovfl_Hs = csw*vegn_ovfl_s*(cc%Tv-tfreeze)

     cc%Wl = cc%Wl - vegn_ovfl_l*delta_time
     cc%Ws = cc%Ws - vegn_ovfl_s*delta_time

     if(is_watch_point()) then
        write(*,*)'#### vegn_step_2 output #####'
        __DEBUG3__(vegn_melt, vegn_ovfl_l, vegn_ovfl_s)
        __DEBUG2__(vegn_ovfl_Hl,vegn_ovfl_Hs)
     endif
     ! accumulate total values, recalculated per unit area of cc
     total_vegn_melt    = total_vegn_melt    + vegn_melt*cc%nindivs
     total_vegn_ovfl_l  = total_vegn_ovfl_l  + vegn_ovfl_l*cc%nindivs
     total_vegn_ovfl_s  = total_vegn_ovfl_s  + vegn_ovfl_s*cc%nindivs
     total_vegn_ovfl_Hl = total_vegn_ovfl_Hl + vegn_ovfl_Hl*cc%nindivs
     total_vegn_ovfl_Hs = total_vegn_ovfl_Hs + vegn_ovfl_Hs*cc%nindivs
  enddo
  ! take into account water dropped with dead plants, etc
  total_vegn_ovfl_l  = total_vegn_ovfl_l  + vegn%drop_wl/delta_time
  total_vegn_ovfl_s  = total_vegn_ovfl_s  + vegn%drop_ws/delta_time
  total_vegn_ovfl_Hl = total_vegn_ovfl_Hl + vegn%drop_hl/delta_time
  total_vegn_ovfl_Hs = total_vegn_ovfl_Hs + vegn%drop_hs/delta_time
  ! reset buffers for the next time step
  vegn%drop_wl = 0 ; vegn%drop_ws = 0
  vegn%drop_hl = 0 ; vegn%drop_hs = 0

  ! ---- diagnostic section
  ! TODO: invent a way to aggregate diagnostic fields
  ! leaf_size,leaf rad. prop. should probably be averaged with LAI+SAI as weight
  ! root_zeta -- perhaps averaged with root density as weight?
  ! snow_crit???
  associate(c=>vegn%cohorts)
  call send_cohort_data(id_height_ave, diag, c(1:N), c(1:N)%height, weight=c(1:N)%nindivs, op=OP_AVERAGE)
  ! TODO: calculate vegetation temperature as total sensible heat/total heat capacity
  call send_cohort_data(id_temp, diag, c(1:N), c(1:N)%Tv, weight=c(1:N)%nindivs, op=OP_AVERAGE)
  call send_cohort_data(id_wl,   diag, c(1:N), c(1:N)%Wl, weight=c(1:N)%nindivs, op=OP_SUM)
  call send_cohort_data(id_ws,   diag, c(1:N), c(1:N)%Ws, weight=c(1:N)%nindivs, op=OP_SUM)

  call send_tile_data(id_height, maxval(c(1:N)%height), diag) ! tallest
  ! in principle, the first cohort must be the tallest, but since cohorts are
  ! rearranged only once a year, that may not be true for part of the year
  call send_cohort_data(id_lai,     diag, c(1:N), c(1:N)%lai, weight=c(1:N)%layerfrac, op=OP_SUM)
  call send_cohort_data(id_laii,    diag, c(1:N), c(1:N)%lai, weight=c(1:N)%nindivs,   op=OP_AVERAGE)
  call send_cohort_data(id_lai_var, diag, c(1:N), c(1:N)%lai, weight=c(1:N)%layerfrac, op=OP_SUM)
  ! these are LAI variance and standard deviation among *tiles*, not cohorts. So the same data is sent
  ! as for average LAI, but they are aggregated differently by the diagnostics
  call send_cohort_data(id_lai_std, diag, c(1:N), c(1:N)%lai, weight=c(1:N)%layerfrac, op=OP_SUM)
  call send_cohort_data(id_sai,     diag, c(1:N), c(1:N)%sai, weight=c(1:N)%layerfrac, op=OP_SUM)
!  call send_cohort_data(id_leafarea,  diag, c(1:N), c(1:N)%leafarea, weight=c(1:N)%nindivs, op=OP_SUM) -- same as LAI (checked)
  call send_cohort_data(id_leafarea, diag, c(1:N), c(1:N)%leafarea, weight=c(1:N)%nindivs, op=OP_AVERAGE)

  ! leaf size averaging weight is the number of leaves in cohort
  call send_cohort_data(id_leaf_size,    diag, c(1:N), c(1:N)%leaf_size,    weight=c(1:N)%nindivs*c(1:N)%leafarea/c(1:N)%leaf_size, op=OP_AVERAGE)
  call send_cohort_data(id_root_density, diag, c(1:N), c(1:N)%root_density, weight=c(1:N)%nindivs,                 op=OP_SUM)
  call send_cohort_data(id_root_zeta,    diag, c(1:N), c(1:N)%root_zeta,    weight=c(1:N)%nindivs*c(1:N)%br,       op=OP_AVERAGE)
  call send_cohort_data(id_rs_min,       diag, c(1:N), c(1:N)%rs_min,       weight=c(1:N)%nindivs*c(1:N)%leafarea, op=OP_AVERAGE)
!   TODO: implement sending multi-dimensional cohort data and enable leaf_refl, leaf_tran diagnostics
!   call send_cohort_data(id_leaf_refl,    diag, c(1:N), c(1:N)%leaf_refl,    weight=c(1:N)%nindivs*c(1:N)%leafarea, op=OP_AVERAGE)
!   call send_cohort_data(id_leaf_tran,    diag, c(1:N), c(1:N)%leaf_tran,    weight=c(1:N)%nindivs*c(1:N)%leafarea, op=OP_AVERAGE)
  call send_cohort_data(id_leaf_emis,    diag, c(1:N), c(1:N)%leaf_emis,    weight=c(1:N)%nindivs*c(1:N)%leafarea, op=OP_AVERAGE)
  ! snow_crit averaging weight is the cohort leaf area
  call send_cohort_data(id_snow_crit,    diag, c(1:N), c(1:N)%snow_crit,    weight=c(1:N)%nindivs*c(1:N)%leafarea, op=OP_AVERAGE)

  ! CMOR/CMIP variables
  if (id_lai_cmor>0) call send_tile_data(id_lai_cmor, sum(c(1:N)%nindivs*c(1:N)%leafarea), diag)
  if (id_cw>0)       call send_tile_data(id_cw, sum((c(1:N)%Wl+c(1:N)%Ws)*c(1:N)%nindivs), diag)
  end associate
end subroutine vegn_step_2


! ============================================================================
! do the vegetation calculations that require updated (end-of-timestep) values
! of prognostic land variables
subroutine vegn_step_3(vegn, soil, cana_T, precip, ndep_nit, ndep_amm, ndep_org, vegn_fco2, diag)
  type(vegn_tile_type), intent(inout) :: vegn
  type(soil_tile_type), intent(inout) :: soil
  real, intent(in) :: cana_T ! canopy temperature, deg K
  real, intent(in) :: precip ! total (rain+snow) precipitation, kg/(m2 s)
  real, intent(in) :: ndep_nit, ndep_amm, ndep_org ! total nitrate, ammonium,
      ! and organic nitrogen inputs (deposition plus fertilization), kg N/(m2 yr)
  real, intent(out) :: vegn_fco2 ! co2 flux from vegetation, kg CO2/(m2 s)
  type(diag_buff_type), intent(inout) :: diag

  ! ---- local vars
  real :: tsoil ! average temperature of soil for soil carbon decomposition, deg K
  real :: theta ! average soil wetness, unitless
  real :: psist ! psi stress index
  real :: depth_ave! depth for averaging soil moisture based on Jackson function for root distribution
  real :: percentile = 0.95
  real :: harv_pool_nitrogen_loss(N_HARV_POOLS)

  tsoil = soil_ave_temp (soil,soil_carbon_depth_scale)
  ! depth for 95% of root according to Jackson distribution
  depth_ave = -log(1.-percentile)*vegn%cohorts(1)%root_zeta

  theta = soil_ave_theta1(soil, depth_ave)

  if(is_watch_point()) then
     write(*,*)'#### vegn_step_3 drought input ####'
     __DEBUG3__(depth_ave, tsoil, theta)
  endif

  ! Do N deposition first. For now, it all goes to leaf litter
  call check_var_range(ndep_amm, 0.0, HUGE(1.0), 'vegn_step_3', 'ndep_amm', FATAL)
  call check_var_range(ndep_nit, 0.0, HUGE(1.0), 'vegn_step_3', 'ndep_nit', FATAL)
  call check_var_range(ndep_org, 0.0, HUGE(1.0), 'vegn_step_3', 'ndep_org', FATAL)
  call soil_NH4_deposition(ndep_amm*dt_fast_yr,soil%litter(LEAF))
  call soil_NO3_deposition(ndep_nit*dt_fast_yr,soil%litter(LEAF))
  call soil_org_N_deposition(ndep_org*dt_fast_yr,soil%litter(LEAF))

  soil%gross_nitrogen_flux_into_tile = soil%gross_nitrogen_flux_into_tile + (ndep_amm+ndep_nit+ndep_org)*dt_fast_yr

  if (do_ppa) then
     call vegn_carbon_int_ppa(vegn, soil, tsoil, theta, diag)
  else
     call vegn_carbon_int_lm3(vegn, soil, tsoil, theta, diag)
  endif

  ! decrease, if necessary, csmoke spending rate so that csmoke pool
  ! is never depleted below zero
  vegn%csmoke_rate = max( 0.0, &
       min( vegn%csmoke_rate, &
            vegn%csmoke_pool/dt_fast_yr)&
       )
  ! update smoke pool -- stored amount of carbon lost to fire
  vegn%csmoke_pool = vegn%csmoke_pool - &
       vegn%csmoke_rate*dt_fast_yr
  soil%gross_nitrogen_flux_out_of_tile = soil%gross_nitrogen_flux_out_of_tile+vegn%nsmoke_pool*dt_fast_yr
  vegn%nsmoke_pool = vegn%nsmoke_pool - vegn%nsmoke_pool*dt_fast_yr ! Following csmoke_rate, which is set to equal csmoke_pool in vegn_disturbance
  ! decrease harvested rates so that pools are not depleted below zero
  vegn%harv_rate_C(:) = max( 0.0, &
                           min(vegn%harv_rate_C(:), vegn%harv_pool_C(:)/dt_fast_yr) &
                         )
  ! update harvested pools -- amounts of stored harvested carbon by category
  vegn%harv_pool_C(:) = vegn%harv_pool_C(:) - &
       vegn%harv_rate_C(:)*dt_fast_yr
  ! --- calculate total co2 flux from vegetation
  vegn_fco2 = -vegn%nep + vegn%csmoke_rate + sum(vegn%harv_rate_C(:))
  ! --- convert it to kg CO2/(m2 s)
  vegn_fco2 = vegn_fco2*mol_CO2/(mol_C*seconds_per_year)

  ! What to do with harvested nitrogen??
  where(harvest_spending_time(:)>0)
     harv_pool_nitrogen_loss = &
          vegn%harv_pool_N(:)/harvest_spending_time(:)
  elsewhere
     harv_pool_nitrogen_loss(:) = 0.0
  end where
  vegn%harv_pool_N = vegn%harv_pool_N - harv_pool_nitrogen_loss*dt_fast_yr
  soil%gross_nitrogen_flux_out_of_tile = soil%gross_nitrogen_flux_out_of_tile+sum(harv_pool_nitrogen_loss)*dt_fast_yr


  ! --- accumulate values for climatological averages
  vegn%tc_av     = vegn%tc_av + cana_T
  vegn%tsoil_av  = vegn%tsoil_av + tsoil
  vegn%precip_av = vegn%precip_av + precip
  if (xwilt_available) then
     theta = soil_ave_theta1(soil,depth_ave)
  else
     theta = soil_ave_theta0(soil,vegn%cohorts(1)%root_zeta)
  endif
  if (do_ppa) then
     vegn%theta_av_phen = weight_av_phen*theta + (1-weight_av_phen)*vegn%theta_av_phen
  else
     vegn%theta_av_phen = vegn%theta_av_phen + theta
  endif
  vegn%theta_av_fire = vegn%theta_av_fire + soil_ave_theta1(soil,depth_ave)
  psist = soil_psi_stress(soil,vegn%cohorts(1)%root_zeta)
  vegn%psist_av  = vegn%psist_av + psist

  vegn%n_accum   = vegn%n_accum+1

  ! --- accumulate values for daily averaging
  vegn%tc_daily    = vegn%tc_daily + cana_T

  ! ---- update daily min and max temperature, for GDD calculations
  vegn%daily_t_max = max(vegn%daily_t_max, cana_T)
  vegn%daily_t_min = min(vegn%daily_t_min, cana_T)

  call send_tile_data(id_theph, theta, diag)
  call send_tile_data(id_psiph, psist, diag)

end subroutine vegn_step_3


! ============================================================================
! given a vegetation tile with the state variables set up, calculate derived
! parameters to get a consistent state
! NOTE: this subroutine does not call update_biomass_pools, although some
! of the calculations are the same. The reason is because this function may
! be used in the situation when the biomasses are not precisely consistent, for
! example when they come from the data override or from initial conditions.
subroutine update_derived_vegn_data(vegn)
  type(vegn_tile_type), intent(inout) :: vegn

  ! ---- local vars
  type(vegn_cohort_type), pointer :: cc ! pointer to the current cohort
  integer :: k  ! cohort index
  integer :: sp ! shorthand for the vegetation species
  integer :: n_layers ! number of layers in cohort
  real, allocatable :: layer_area(:) ! total area of crowns in the layer
  integer :: current_layer
  real :: zbot ! height of the bottom of the canopy, m (=top of the lower layer)
  real :: stemarea ! individual stem area, for SAI calculations, m2/individual
  real :: VRL(num_l) ! vertical distribution of volumetric root length, m/m3

  ! determine layer boundaries in the array of cohorts
  n_layers = maxval(vegn%cohorts(:)%layer)
  allocate (layer_area(n_layers))

  ! calculate total area of canopies per layer (per unit tile area)
  layer_area(:) = 0
  do k = 1, vegn%n_cohorts
     cc=>vegn%cohorts(k)
     layer_area(cc%layer) = layer_area(cc%layer) + cc%crownarea*cc%nindivs
  enddo

  ! limit layer area so that it only squeezes the cohorts if the total area of
  ! the canopies is higher than tile area
  if (allow_external_gaps) then
     do k = 1,n_layers
        layer_area(k) = max(layer_area(k),1.0)
     enddo
  endif

  ! protect from zero layer area situation: this can happen if all cohorts die
  ! due to mortality or starvation. In this case n_layers is 1, of course.
  do k = 1,n_layers
     if (layer_area(k)<=0) layer_area(k) = 1.0
  enddo


  if (is_watch_point()) then
     write(*,*)'#### update_derived_vegn_data ####'
  endif
  ! given that the cohort state variables are initialized, fill in
  ! the intermediate variables
  do k = 1,vegn%n_cohorts
    cc=>vegn%cohorts(k)

    sp = cc%species
    ! update fractions of the living biomass
    if (.not.do_ppa) then
       cc%height = height_from_biomass(btotal(cc))
    endif
    call update_bio_living_fraction(cc) ! this should not have any effect in PPA,
    ! since it only updates Px fractions of bliving, but I am not sure if this is
    ! implemented consistently right now.
    ! TODO: check that Pl, Pr, Psw, Psw_alphasw are not used in PPA, move the
    ! above call inside "if (.not.do_ppa)" statement

    if(sp<NSPECIES) then ! LM3V species
       ! calculate area fraction that the cohort occupies in its layer
!       if (layer_area(cc%layer)<=0) call error_mesg('update_derived_vegn_data', &
!          'total area of canopy layer '//string(cc%layer)//' is zero', FATAL)
       cc%layerfrac = cc%crownarea*cc%nindivs*(1-spdata(sp)%internal_gap_frac)/layer_area(cc%layer)
       ! calculate the leaf area index based on the biomass of leaves
       ! leaf_area_from_biomass returns the total area of leaves per individual;
       ! convert it to leaf area per m2, and re-normalize to take into account
       ! stretching of canopies
       cc%leafarea = leaf_area_from_biomass(cc%bl, sp, cc%layer, cc%firstlayer)
       cc%lai = cc%leafarea/(cc%crownarea*(1-spdata(sp)%internal_gap_frac))*layer_area(cc%layer)
       ! calculate the root density as the total biomass below ground, in
       ! biomass (not carbon!) units
       cc%root_density = (cc%br + (cc%bsw+cc%bwood+cc%blv)*(1-agf_bs))*C2B
    else
       cc%height        = spdata(sp)%dat_height
       cc%lai           = spdata(sp)%dat_lai
       cc%root_density  = spdata(sp)%dat_root_density
    endif
!     stemarea      = 0.035*cc%height ! Federer and Lash, 1978
!     cc%sai           = stemarea/cc%crownarea*layer_area(cc%layer)
    cc%sai           = spdata(sp)%sai_height_ratio * cc%height ! Federer and Lash, 1978
    cc%leaf_size     = spdata(sp)%leaf_size
    cc%root_zeta     = spdata(sp)%dat_root_zeta
    cc%rs_min        = spdata(sp)%dat_rs_min
    cc%leaf_refl     = spdata(sp)%leaf_refl
    cc%leaf_tran     = spdata(sp)%leaf_tran
    cc%leaf_emis     = spdata(sp)%leaf_emis
    cc%snow_crit     = spdata(sp)%dat_snow_crit

    ! the following variables are per individual
    cc%Wl_max        = spdata(sp)%cmc_lai*cc%leafarea
    cc%Ws_max        = spdata(sp)%csc_lai*cc%leafarea
    cc%mcv_dry       = max(mcv_min, mcv_lai*cc%leafarea)
    if (is_watch_point()) then
       write(*,'(i2.2,2x,":")',advance='NO') k
       call dpri('height',cc%height)
       call dpri('LAI',cc%lai)

       call dpri('bl',cc%bl)
       call dpri('leafarea',cc%leafarea)
       call dpri('crownarea',cc%crownarea)
       call dpri('nindivs',cc%nindivs)
       ! call dpri('gapfrac',spdata(sp)%internal_gap_frac)
       ! call dpri('layerarea',layer_area(cc%layer))
       call dpri('layer',cc%layer)
       call dpri('species',spdata(sp)%name)

       write(*,*)
    endif
  enddo

  ! Calculate height of the canopy bottom: equals to the top of the lower layer.
  ! this code assumes that cohorts are arranged in descending order
  zbot = 0; current_layer = vegn%cohorts(vegn%n_cohorts)%layer
  do k = vegn%n_cohorts, 1, -1
    if (vegn%cohorts(k)%layer/=current_layer) then
       zbot = vegn%cohorts(k+1)%height
       current_layer = vegn%cohorts(k)%layer
    endif
    vegn%cohorts(k)%zbot = zbot
  enddo

  ! calculate volumetric root length for the entire tile
  VRL(:) = 0.0
  do k = 1, vegn%n_cohorts
     cc=>vegn%cohorts(k)
     call cohort_root_properties (cc, dz(1:num_l), cc%root_length(1:num_l), cc%K_r, cc%r_r)
     VRL(:) = VRL(:)+cc%root_length(1:num_l)*cc%nindivs
  enddo

  ! calculate characteristic half-distance between roots, m
  where (VRL(:) > 0)
     vegn%root_distance(1:num_l) = 1.0/sqrt(PI*VRL(:))
  elsewhere
     vegn%root_distance(1:num_l) = 1.0 ! the value does not matter since uptake is 0 anyway
  end where

  deallocate(layer_area)
end subroutine update_derived_vegn_data


! ============================================================================
! update slow components of the vegetation model
subroutine update_vegn_slow( )

  ! ---- local vars ----------------------------------------------------------
  integer :: second, minute, hour, day0, day1, month0, month1, year0, year1, doy
  type(land_tile_enum_type) :: ce
  type(land_tile_type), pointer :: tile
  integer :: i,j,k,l ! current point indices
  integer :: ii ! pool and cohort iterator
  integer :: N ! number of cohorts
  real    :: weight_ncm ! low-pass filter value for the number of cold months
  type(vegn_cohort_type), pointer :: cc(:) ! shorthand for cohort array
  real    :: zstar ! critical height, for diag only
  character(64) :: str
  real, allocatable :: btot(:) ! storage for total biomass
  real :: dheat ! heat residual due to cohort merging
  real :: w ! smoothing weight
  real :: tc_daily ! daly temperature for dormancy detection, degK

  ! variables for conservation checks
  real :: lmass0, fmass0, heat0, cmass0, nmass0
  real :: dbh_max_N ! max dbh for understory; diag only

  ! get components of calendar dates for this and previous time step
  call get_date(lnd%time-lnd%dt_slow, year1,month1,day1,hour,minute,second)
  call get_date(lnd%time,             year0,month0,day0,hour,minute,second)
  ! calculate day of year, to avoid re-calculating it in harvesting
  doy = day_of_year(lnd%time)

  if(month0 /= month1) then
     ! heartbeat
     write(str,'("Current date is ",i4.4,"-",i2.2,"-",i2.2)') year0,month0,day0
     call error_mesg('update_vegn_slow',trim(str),NOTE)
  endif

  call update_fire_data(lnd%time)

  ce = first_elmt(land_tile_map, lnd%ls)
  do while (loop_over_tiles(ce,tile,l,k))
     call set_current_point(l,k) ! this is for debug output only
     if(.not.associated(tile%vegn)) cycle ! skip the rest of the loop body

     ! + conservation check, part 1: calculate the pre-transition totals
     call check_conservation_1(tile,lmass0,fmass0,cmass0,nmass0)

     if (day1 /= day0) then
        do ii = 1, tile%vegn%n_cohorts
           associate (cc=>tile%vegn%cohorts(ii), sp=>spdata(tile%vegn%cohorts(ii)%species)) ! F2003
           cc%npp_previous_day     = cc%npp_previous_day_tmp/steps_per_day
           cc%npp_previous_day_tmp = 0.0
           ! GDD calculated as in http://en.wikipedia.org/wiki/Growing_degree-day
           cc%gdd = cc%gdd + &
               (max(tile%vegn%daily_T_min-sp%gdd_base_T,0.0) + &
                max(tile%vegn%daily_T_max-sp%gdd_base_T,0.0))/2
           end associate ! F2003
        enddo
        tile%vegn%tc_pheno = tile%vegn%tc_pheno * 0.95 + tile%vegn%tc_daily * 0.05
     endif

     call check_conservation_2(tile,'update_vegn_slow 1',lmass0,fmass0,cmass0,nmass0)

     ! monthly averaging
     if (month1 /= month0) then
        ! compute averages from accumulated monthly values
        tile%vegn%tc_av     = tile%vegn%tc_av     / tile%vegn%n_accum
        tile%vegn%tsoil_av  = tile%vegn%tsoil_av  / tile%vegn%n_accum
        if (do_ppa) then
           ! do nothing -- theta_av_phen is already smoothed
        else
           tile%vegn%theta_av_phen  = tile%vegn%theta_av_phen  / tile%vegn%n_accum
        endif
        tile%vegn%theta_av_fire  = tile%vegn%theta_av_fire  / tile%vegn%n_accum
        tile%vegn%psist_av  = tile%vegn%psist_av  / tile%vegn%n_accum
        tile%vegn%precip_av = tile%vegn%precip_av / tile%vegn%n_accum
        ! accumulate annual values
        tile%vegn%p_ann_acm = tile%vegn%p_ann_acm+tile%vegn%precip_av
        tile%vegn%t_ann_acm = tile%vegn%t_ann_acm+tile%vegn%tc_av
        if ( tile%vegn%tc_av < cold_month_threshold ) &
             tile%vegn%ncm_acm = tile%vegn%ncm_acm+1
        tile%vegn%t_cold_acm = min(tile%vegn%t_cold_acm, tile%vegn%tc_av)

        tile%vegn%nmn_acm = tile%vegn%nmn_acm+1 ! increase the number of accumulated months

        ! Calculate monthly average allocation to each N uptake type, then save
        ! the maximum month from each year. So we can then smooth a yearly maximum
        ! estimate for limiting the change over time without messing everything up
        ! because of seasonality
        do i = 1,tile%vegn%n_cohorts
           associate(c=>tile%vegn%cohorts(i))
           c%max_monthly_Nfix_alloc = max(c%max_monthly_Nfix_alloc, c%Nfix_alloc_accum/(tile%vegn%n_accum*dt_fast_yr))
           c%max_monthly_mine_alloc = max(c%max_monthly_mine_alloc, c%mine_alloc_accum/(tile%vegn%n_accum*dt_fast_yr))
           c%max_monthly_scav_alloc = max(c%max_monthly_scav_alloc, c%scav_alloc_accum/(tile%vegn%n_accum*dt_fast_yr))

           c%Nfix_alloc_accum = 0.0
           c%scav_alloc_accum = 0.0
           c%mine_alloc_accum = 0.0
           end associate
        enddo
     endif

     call check_conservation_2(tile,'update_vegn_slow 2',lmass0,fmass0,cmass0,nmass0)

     ! annual averaging
     if (year1 /= year0) then
        ! The ncm smoothing is coded as a low-pass exponential filter. See, for
        ! example http://en.wikipedia.org/wiki/Low-pass_filter
        weight_ncm = 1/(1+tau_smooth_ncm)
        if(tile%vegn%nmn_acm /= 0) then
           ! calculate annual averages from accumulated values
           tile%vegn%p_ann  = tile%vegn%p_ann_acm/tile%vegn%nmn_acm
           tile%vegn%t_ann  = tile%vegn%t_ann_acm/tile%vegn%nmn_acm
           tile%vegn%t_cold = tile%vegn%t_cold_acm
           tile%vegn%ncm    = weight_ncm*tile%vegn%ncm_acm + (1-weight_ncm)*tile%vegn%ncm
           ! reset accumulated values
           tile%vegn%ncm_acm    = 0
           tile%vegn%p_ann_acm  = 0
           tile%vegn%t_ann_acm  = 0
           tile%vegn%t_cold_acm = HUGE(tile%vegn%t_cold_acm)
        endif
        tile%vegn%nmn_acm = 0

        do i = 1,tile%vegn%n_cohorts
           associate(c=>tile%vegn%cohorts(i))
           w = 1.0/(1+spdata(c%species)%tau_smooth_alloc)
           c%max_Nfix_allocation = w*c%max_monthly_Nfix_alloc + (1-w)*c%max_Nfix_allocation
           c%max_mine_allocation = w*c%max_monthly_mine_alloc + (1-w)*c%max_mine_allocation
           c%max_scav_allocation = w*c%max_monthly_scav_alloc + (1-w)*c%max_scav_allocation

           c%max_monthly_Nfix_alloc = 0.0
           c%max_monthly_mine_alloc = 0.0
           c%max_monthly_scav_alloc = 0.0
           end associate
        enddo
      endif

     if (year1 /= year0 .and. do_biogeography) then
        call vegn_biogeography(tile%vegn)
     endif

     call check_conservation_2(tile,'update_vegn_slow 3',lmass0,fmass0,cmass0,nmass0)

     if (year1 /= year0 .and. do_peat_redistribution) then
        call redistribute_peat_carbon(tile%soil)
     endif

     if (month1 /= month0.and.do_patch_disturbance) then
        call update_fuel(tile%vegn,tile%soil%w_wilt(1)/tile%soil%pars%vwc_sat)
        ! assume that all layers are the same soil type and wilting is vertically homogeneous
     endif
     call check_conservation_2(tile,'update_vegn_slow 4',lmass0,fmass0,cmass0,nmass0)

     if (day1 /= day0 .and. do_cohort_dynamics) then
        N = tile%vegn%n_cohorts ; cc=>tile%vegn%cohorts(1:N)
        call send_tile_data(id_cgain,sum(cc(1:N)%carbon_gain*cc(1:N)%nindivs),tile%diag)
        call send_tile_data(id_closs,sum(cc(1:N)%carbon_loss*cc(1:N)%nindivs),tile%diag)
        call send_tile_data(id_wdgain,sum(cc(1:N)%bwood_gain*cc(1:N)%nindivs),tile%diag)

        call send_tile_data(id_Ngain,sum(cc(1:N)%nitrogen_gain*cc(1:N)%nindivs),tile%diag)
        call send_tile_data(id_Nloss,sum(cc(1:N)%nitrogen_loss*cc(1:N)%nindivs),tile%diag)

        call vegn_growth(tile%vegn, tile%vegn%tc_daily, tile%diag) ! selects lm3 or ppa inside
        call check_conservation_2(tile,'update_vegn_slow 4.1',lmass0,fmass0,cmass0)

        if (do_ppa) then
           call vegn_starvation_ppa(tile%vegn, tile%soil)
           call check_conservation_2(tile,'update_vegn_slow 4.2',lmass0,fmass0,cmass0,nmass0)
           call vegn_phenology_ppa (tile)
           call check_conservation_2(tile,'update_vegn_slow 4.3',lmass0,fmass0,cmass0,nmass0)
        else
           call vegn_nat_mortality_lm3(tile%vegn,tile%soil,86400.0)
        endif
     endif
     call check_conservation_2(tile,'update_vegn_slow 5',lmass0,fmass0,cmass0,nmass0)

     if  (month1 /= month0 .and. do_phenology) then
        if (.not.do_ppa) &
            call vegn_phenology_lm3 (tile%vegn,tile%soil)
        ! assume that all layers are the same soil type and wilting is vertically homogeneous
     endif
     call check_conservation_2(tile,'update_vegn_slow 6',lmass0,fmass0,cmass0,nmass0)

     if (year1 /= year0 .AND. fire_option==FIRE_LM3 .AND. do_patch_disturbance) then
        call vegn_disturbance(tile%vegn, tile%soil, seconds_per_year)
     endif
     call check_conservation_2(tile,'update_vegn_slow 7',lmass0,fmass0,cmass0,nmass0)

     call vegn_harvesting(tile, year0/=year1, month0/=month1, day0/=day1, doy, l)

     if (year1 /= year0) then
        tile%vegn%fsc_rate_ag = tile%vegn%fsc_pool_ag/fsc_pool_spending_time
        tile%vegn%ssc_rate_ag = tile%vegn%ssc_pool_ag/ssc_pool_spending_time
        tile%vegn%fsc_rate_bg = tile%vegn%fsc_pool_bg/fsc_pool_spending_time
        tile%vegn%ssc_rate_bg = tile%vegn%ssc_pool_bg/ssc_pool_spending_time
        tile%vegn%fsn_rate_bg = tile%vegn%fsn_pool_bg/fsc_pool_spending_time
        tile%vegn%ssn_rate_bg = tile%vegn%ssn_pool_bg/ssc_pool_spending_time

        tile%vegn%litter_rate_C(C_FAST,:) = tile%vegn%litter_buff_C(C_FAST,:)/fsc_pool_spending_time
        tile%vegn%litter_rate_C(C_SLOW,:) = tile%vegn%litter_buff_C(C_SLOW,:)/ssc_pool_spending_time
        tile%vegn%litter_rate_N(C_FAST,:) = tile%vegn%litter_buff_N(C_FAST,:)/fsc_pool_spending_time
        tile%vegn%litter_rate_N(C_SLOW,:) = tile%vegn%litter_buff_N(C_SLOW,:)/ssc_pool_spending_time
        where(harvest_spending_time(:)>0)
           tile%vegn%harv_rate_C(:) = &
                tile%vegn%harv_pool_C(:)/harvest_spending_time(:)
        elsewhere
           tile%vegn%harv_rate_C(:) = 0.0
        end where
     endif
     call check_conservation_2(tile,'update_vegn_slow 9',lmass0,fmass0,cmass0,nmass0)

     ! + sanity checks
     do ii = 1,tile%vegn%n_cohorts
        associate(cc => tile%vegn%cohorts(ii))
        str = 'cohort('//trim(string(ii))//')%'
        call check_var_range(cc%nindivs, 0.0, HUGE(1.0), 'update_vegn_slow', trim(str)//'nindivs', WARNING)
        call check_var_range(cc%bl,      0.0, HUGE(1.0), 'update_vegn_slow', trim(str)//'bl',      WARNING)
        call check_var_range(cc%blv,     0.0, HUGE(1.0), 'update_vegn_slow', trim(str)//'blv',     WARNING)
        call check_var_range(cc%br,      0.0, HUGE(1.0), 'update_vegn_slow', trim(str)//'br',      WARNING)
        call check_var_range(cc%bsw,     0.0, HUGE(1.0), 'update_vegn_slow', trim(str)//'bsw',     WARNING)
        call check_var_range(cc%bwood,   0.0, HUGE(1.0), 'update_vegn_slow', trim(str)//'bwood',   WARNING)
        !call check_var_range(cc%nsc,     0.0, HUGE(1.0), 'update_vegn_slow', trim(str)//'nsc',     WARNING)
        call check_var_range(cc%bseed,   0.0, HUGE(1.0), 'update_vegn_slow', trim(str)//'bseed',   WARNING)
        end associate
     enddo
     ! - sanity checks

     call check_conservation_2(tile,'update_vegn_slow',lmass0,fmass0,cmass0,nmass0)

     ! perhaps we need to move that either inside vegn_reproduction_ppa, or after that
     if (do_ppa.and.year1 /= year0) then
        call vegn_relayer_cohorts_ppa(tile%vegn)
        call check_conservation_2(tile,'update_vegn_slow 7.2',lmass0,fmass0,cmass0,nmass0)
        call vegn_mergecohorts_ppa(tile%vegn, dheat)
        tile%e_res_2 = tile%e_res_2 - dheat
        call check_conservation_2(tile,'update_vegn_slow 7.3',lmass0,fmass0,cmass0,nmass0)
        ! update DBH_ys
        do ii = 1, tile%vegn%n_cohorts
           tile%vegn%cohorts(ii)%DBH_ys = tile%vegn%cohorts(ii)%dbh
           tile%vegn%cohorts(ii)%BM_ys  = tile%vegn%cohorts(ii)%bsw + &
                                          tile%vegn%cohorts(ii)%bwood
        enddo
     endif

     if (do_ppa.and.day1 /= day0) then
        call kill_small_cohorts_ppa(tile%vegn,tile%soil)
        call check_conservation_2(tile,'update_vegn_slow 8',lmass0,fmass0,cmass0)
     endif

     ! ---- diagnostic section
     call send_tile_data(id_t_ann,   tile%vegn%t_ann,   tile%diag)
     call send_tile_data(id_t_cold,  tile%vegn%t_cold,  tile%diag)
     call send_tile_data(id_lambda,  tile%vegn%lambda,  tile%diag)
     call send_tile_data(id_p_ann,   tile%vegn%p_ann,   tile%diag)
     call send_tile_data(id_ncm,     real(tile%vegn%ncm), tile%diag)
     call send_tile_data(id_afire,   tile%vegn%disturbance_rate(1), tile%diag)
     call send_tile_data(id_atfall,  tile%vegn%disturbance_rate(0), tile%diag)

     call send_tile_data(id_tc_pheno,tile%vegn%tc_pheno,tile%diag)

     do ii = 1,N_HARV_POOLS
        call send_tile_data(id_harv_pool_C(ii),tile%vegn%harv_pool_C(ii),tile%diag)
        call send_tile_data(id_harv_pool_N(ii),tile%vegn%harv_pool_N(ii),tile%diag)
        call send_tile_data(id_harv_rate_C(ii),tile%vegn%harv_rate_C(ii),tile%diag)
     enddo
     call send_tile_data(id_tot_harv_pool_C,sum(tile%vegn%harv_pool_C(:)),tile%diag)
     call send_tile_data(id_tot_harv_rate_C,sum(tile%vegn%harv_rate_C(:)),tile%diag)
     call send_tile_data(id_tot_harv_pool_N,sum(tile%vegn%harv_pool_N(:)),tile%diag)
     call send_tile_data(id_csmoke_pool,tile%vegn%csmoke_pool,tile%diag)
     call send_tile_data(id_csmoke_rate,tile%vegn%csmoke_rate,tile%diag)
     call send_tile_data(id_nsmoke_pool,tile%vegn%nsmoke_pool,tile%diag)

     N=tile%vegn%n_cohorts ; cc=>tile%vegn%cohorts
     call send_cohort_data(id_ncohorts, tile%diag, cc(1:N), (/(1.0,i=1,N)/), op=OP_SUM)
     call send_cohort_data(id_nindivs,  tile%diag, cc(1:N), cc(1:N)%nindivs, op=OP_SUM)

     call send_tile_data(id_nlayers,  real(cc(N)%layer),    tile%diag)

     call send_cohort_data(id_bl,     tile%diag, cc(1:N), cc(1:N)%bl,     weight=cc(1:N)%nindivs, op=OP_SUM)
     call send_cohort_data(id_blv,    tile%diag, cc(1:N), cc(1:N)%blv,    weight=cc(1:N)%nindivs, op=OP_SUM)
     call send_cohort_data(id_br,     tile%diag, cc(1:N), cc(1:N)%br,     weight=cc(1:N)%nindivs, op=OP_SUM)
     call send_cohort_data(id_bsw,    tile%diag, cc(1:N), cc(1:N)%bsw,    weight=cc(1:N)%nindivs, op=OP_SUM)
     call send_cohort_data(id_bwood,  tile%diag, cc(1:N), cc(1:N)%bwood,  weight=cc(1:N)%nindivs, op=OP_SUM)
     call send_cohort_data(id_bseed,  tile%diag, cc(1:N), cc(1:N)%bseed,  weight=cc(1:N)%nindivs, op=OP_SUM)
     call send_cohort_data(id_nsc,    tile%diag, cc(1:N), cc(1:N)%nsc,    weight=cc(1:N)%nindivs, op=OP_SUM)

     call send_cohort_data(id_leaf_N,      tile%diag, cc(1:N), cc(1:N)%leaf_N,    weight=cc(1:N)%nindivs, op=OP_SUM)
     call send_cohort_data(id_stored_N,    tile%diag, cc(1:N), cc(1:N)%stored_N,  weight=cc(1:N)%nindivs, op=OP_SUM)
     call send_cohort_data(id_root_N,      tile%diag, cc(1:N), cc(1:N)%root_N,    weight=cc(1:N)%nindivs, op=OP_SUM)
     call send_cohort_data(id_sapwood_N,   tile%diag, cc(1:N), cc(1:N)%sapwood_N, weight=cc(1:N)%nindivs, op=OP_SUM)
     call send_cohort_data(id_wood_N,      tile%diag, cc(1:N), cc(1:N)%wood_N,    weight=cc(1:N)%nindivs, op=OP_SUM)
     call send_cohort_data(id_seed_N,      tile%diag, cc(1:N), cc(1:N)%seed_N,    weight=cc(1:N)%nindivs, op=OP_SUM)
     call send_cohort_data(id_Ntot,        tile%diag, cc(1:N), cc(1:N)%total_N,   weight=cc(1:N)%nindivs, op=OP_SUM)

     call send_cohort_data(id_mrz_scav_C, tile%diag, cc(1:N), cc(1:N)%myc_scavenger_biomass_C, weight=cc(1:N)%nindivs, op=OP_SUM)
     call send_cohort_data(id_mrz_mine_C, tile%diag, cc(1:N), cc(1:N)%myc_miner_biomass_C, weight=cc(1:N)%nindivs, op=OP_SUM)
     call send_cohort_data(id_Nfix_C,     tile%diag, cc(1:N), cc(1:N)%N_fixer_biomass_C, weight=cc(1:N)%nindivs, op=OP_SUM)
     call send_cohort_data(id_mrz_scav_N, tile%diag, cc(1:N), cc(1:N)%myc_scavenger_biomass_N, weight=cc(1:N)%nindivs, op=OP_SUM)
     call send_cohort_data(id_mrz_mine_N, tile%diag, cc(1:N), cc(1:N)%myc_miner_biomass_N, weight=cc(1:N)%nindivs, op=OP_SUM)
     call send_cohort_data(id_Nfix_N,     tile%diag, cc(1:N), cc(1:N)%N_fixer_biomass_N, weight=cc(1:N)%nindivs, op=OP_SUM)

     ! ens 021517
     call send_cohort_data(id_brsw,   tile%diag, cc(1:N), cc(1:N)%brsw,    weight=cc(1:N)%nindivs, op=OP_SUM)
     call send_cohort_data(id_growth_prev_day, tile%diag, cc(1:N), cc(1:N)%growth_previous_day, weight=cc(1:N)%nindivs, op=OP_SUM)
     allocate(btot(N)) ! has to be allocated since N cohorts changes inside this subroutine
     btot(1:N) = cc(1:N)%bl    + &
                 cc(1:N)%blv   + &
                 cc(1:N)%br    + &
                 cc(1:N)%bsw   + &
                 cc(1:N)%bwood + &
                 cc(1:N)%bseed + &
                 cc(1:N)%nsc
     call send_cohort_data(id_btot,   tile%diag, cc(1:N),  btot, weight=cc(1:N)%nindivs, op=OP_SUM)
     if (id_cVeg>0) call send_tile_data(id_cVeg, sum(cc(1:N)%nindivs*btot(1:N)), tile%diag) ! CMOR/CMIP diagnostics
     deallocate(btot)

     call send_cohort_data(id_bsw_max, tile%diag, cc(1:N), cc(1:N)%bsw_max, weight=cc(1:N)%nindivs, op=OP_SUM)
     call send_cohort_data(id_bl_max, tile%diag, cc(1:N), cc(1:N)%bl_max, weight=cc(1:N)%nindivs, op=OP_SUM)
     call send_cohort_data(id_br_max, tile%diag, cc(1:N), cc(1:N)%br_max, weight=cc(1:N)%nindivs, op=OP_SUM)
     call send_cohort_data(id_topyear, tile%diag, cc(1:N), cc(1:N)%topyear, weight=cc(1:N)%nindivs, op=OP_AVERAGE)

     call send_cohort_data(id_dbh,       tile%diag, cc(1:N), cc(1:N)%dbh,        weight=cc(1:N)%nindivs, op=OP_AVERAGE)
     call send_cohort_data(id_crownarea, tile%diag, cc(1:N), cc(1:N)%crownarea,  weight=cc(1:N)%nindivs, op=OP_AVERAGE)
     call send_cohort_data(id_dbh_max,   tile%diag, cc(1:N), cc(1:N)%dbh, op=OP_MAX)

     call send_tile_data(id_fsc_pool_ag,tile%vegn%fsc_pool_ag,tile%diag)
     call send_tile_data(id_fsc_rate_ag,tile%vegn%fsc_rate_ag,tile%diag)
     call send_tile_data(id_ssc_pool_ag,tile%vegn%ssc_pool_ag,tile%diag)
     call send_tile_data(id_ssc_rate_ag,tile%vegn%ssc_rate_ag,tile%diag)
     call send_tile_data(id_fsc_pool_bg,tile%vegn%fsc_pool_ag,tile%diag)
     call send_tile_data(id_fsc_rate_bg,tile%vegn%fsc_rate_ag,tile%diag)
     call send_tile_data(id_ssc_pool_bg,tile%vegn%ssc_pool_ag,tile%diag)
     call send_tile_data(id_ssc_rate_bg,tile%vegn%ssc_rate_ag,tile%diag)

     do j = 1,N_LITTER_POOLS
     do i = 1,N_C_TYPES
        ! note that order of indices in IDs and buffers/rates is different, due to difference in conventions
        call send_tile_data(id_litter_buff_C(j,i),tile%vegn%litter_buff_C(i,j),tile%diag)
        call send_tile_data(id_litter_rate_C(j,i),tile%vegn%litter_rate_C(i,j),tile%diag)
        call send_tile_data(id_litter_buff_N(j,i),tile%vegn%litter_buff_N(i,j),tile%diag)
        call send_tile_data(id_litter_rate_N(j,i),tile%vegn%litter_rate_N(i,j),tile%diag)
     enddo
     enddo

     call send_tile_data(id_fuel,    tile%vegn%fuel, tile%diag)

     ! TODO: generalize gdd diagnostics
     call send_tile_data(id_gdd, cc(1)%gdd, tile%diag)
     call send_tile_data(id_species, real(cc(1)%species), tile%diag)
     call send_tile_data(id_status,  real(cc(1)%status),  tile%diag)

     call send_cohort_data(id_leaf_age, tile%diag, cc(1:N), cc(1:N)%leaf_age, weight=cc(1:N)%nindivs*cc(1:N)%leafarea, op=OP_AVERAGE)

     ! carbon budget tracking
     call send_tile_data(id_fsc_in,  sum(tile%soil%fsc_in(:)),  tile%diag)
     call send_tile_data(id_fsc_out, tile%vegn%fsc_out, tile%diag)
     call send_tile_data(id_ssc_in,  sum(tile%soil%ssc_in(:)),  tile%diag)
     call send_tile_data(id_ssc_out, tile%vegn%ssc_out, tile%diag)
     call send_tile_data(id_deadmic_out, tile%vegn%deadmic_out, tile%diag)
     call send_tile_data(id_veg_in,  tile%vegn%veg_in,  tile%diag)
     call send_tile_data(id_veg_out, tile%vegn%veg_out, tile%diag)

     call send_tile_data(id_tile_nitrogen_gain,tile%soil%gross_nitrogen_flux_into_tile, tile%diag)
     call send_tile_data(id_tile_nitrogen_loss,tile%soil%gross_nitrogen_flux_out_of_tile, tile%diag)


     ! CMOR/CMIP variables
     if (id_cLeaf>0) call send_tile_data(id_cLeaf, sum(cc(1:N)%bl*cc(1:N)%nindivs), tile%diag)
     if (id_cWood>0) call send_tile_data(id_cWood, &
         sum((cc(1:N)%bwood+cc(1:N)%bsw)*cc(1:N)%nindivs)*agf_bs, tile%diag)
     if (id_cRoot>0) call send_tile_data(id_cRoot, &
         sum(((cc(1:N)%bwood+cc(1:N)%bsw)*(1-agf_bs)+cc(1:N)%br)*cc(1:N)%nindivs), tile%diag)
     if (id_cMisc>0) call send_tile_data(id_cMisc, sum((cc(1:N)%blv+cc(1:N)%nsc+cc(1:N)%bseed)*cc(1:N)%nindivs), tile%diag)

     if (id_nLeaf>0) call send_tile_data(id_nLeaf, sum(cc(1:N)%leaf_N*cc(1:N)%nindivs), tile%diag)
     if (id_nStem>0) call send_tile_data(id_nStem, &
         sum((cc(1:N)%wood_N+cc(1:N)%sapwood_N)*cc(1:N)%nindivs)*agf_bs, tile%diag)
     if (id_nRoot>0) call send_tile_data(id_nRoot, &
         sum(((cc(1:N)%wood_N+cc(1:N)%sapwood_N)*(1-agf_bs)+cc(1:N)%root_N)*cc(1:N)%nindivs), tile%diag)
     if (id_nOther>0) call send_tile_data(id_nOther, &
         sum((cc(1:N)%stored_N+cc(1:N)%seed_N)*cc(1:N)%nindivs), tile%diag)
     if (id_nVeg>0) call send_tile_data(id_nVeg, &
         sum((cc(1:N)%leaf_N+cc(1:N)%root_N+cc(1:N)%stored_N+cc(1:N)%wood_N+cc(1:N)%sapwood_N+cc(1:N)%seed_N)*cc(1:N)%nindivs), &
         tile%diag)

     if (id_cProduct>0) then
        cmass0 = 0.0
        do i = 1, N_HARV_POOLS
           if (i/=HARV_POOL_CLEARED) cmass0 = cmass0 + tile%vegn%harv_pool_C(i)
        enddo
        call send_tile_data(id_cProduct, cmass0, tile%diag)
     endif

     if (id_cAnt>0) then
        call send_tile_data(id_cAnt, sum(tile%vegn%harv_pool_C(:)), tile%diag)
     endif

     if (id_fGrazing>0) call send_tile_data(id_fGrazing, tile%vegn%harv_rate_C(HARV_POOL_PAST)/seconds_per_year, tile%diag)
     if (id_fHarvest>0) call send_tile_data(id_fHarvest, tile%vegn%harv_rate_C(HARV_POOL_CROP)/seconds_per_year, tile%diag)
     if (id_fProductDecomp>0) then
        cmass0 = 0.0
        do i = 1, N_HARV_POOLS
           if (i/=HARV_POOL_CLEARED) cmass0 = cmass0 + tile%vegn%harv_rate_C(i)
        enddo
        call send_tile_data(id_fProductDecomp, cmass0/seconds_per_year, tile%diag)
     endif
     if (id_fLuc>0) call send_tile_data(id_fLuc, &
            (tile%vegn%harv_rate_C(HARV_POOL_CLEARED) &
            +tile%vegn%harv_rate_C(HARV_POOL_WOOD_FAST) &
            +tile%vegn%harv_rate_C(HARV_POOL_WOOD_MED) &
            +tile%vegn%harv_rate_C(HARV_POOL_WOOD_SLOW) &
            )/seconds_per_year, tile%diag)
     if (id_fAnthDisturb>0) &
            call send_tile_data(id_fAnthDisturb, sum(tile%vegn%harv_rate_C(:))/seconds_per_year, tile%diag)

     ! ---- end of diagnostic section

     ! reset averages and number of steps to 0 before the start of new month
     if (month1 /= month0) then
        tile%vegn%n_accum  = 0
        tile%vegn%tc_av    = 0.
        tile%vegn%tsoil_av = 0.
        if (do_ppa) then
           ! do nothing -- theta_av_phen is being smoothed with low-band-pass filter
        else
           tile%vegn%theta_av_phen = 0.
        endif
        tile%vegn%theta_av_fire = 0.
        tile%vegn%psist_av = 0.
        tile%vegn%precip_av= 0.

        ! Reset nitrogen tracking too
        tile%soil%gross_nitrogen_flux_into_tile = 0.0
        tile%soil%gross_nitrogen_flux_out_of_tile = 0.0
     endif

     !reset fuel and drought months before the start of new year
     if (year1 /= year0) then
        tile%vegn%lambda     = 0
        tile%vegn%fuel       = 0
     endif

     ! zstar diagnostics
     if (id_zstar_1 > 0) then
        do ii = 1, tile%vegn%n_cohorts
           if ( tile%vegn%cohorts(ii)%layer > 1 ) exit
        enddo
        if (ii > tile%vegn%n_cohorts) then
           zstar = 0.0
        else
           zstar = tile%vegn%cohorts(ii)%height
        endif
        call send_tile_data(id_zstar_1, zstar, tile%diag)
     endif
  enddo

  if (do_ppa.and.year1 /= year0) then
    if (do_ppa) then
       call vegn_reproduction_ppa(do_seed_transport) ! includes seed transport.
    else
       ! seed transport in lm3
       call vegn_seed_transport_lm3()
    endif
  endif

  if(soil_carbon_option==SOILC_CORPSE.or.soil_carbon_option==SOILC_CORPSE_N) then
     ! Knock soil carbon cohorts down to their maximum number.
     ! For reproducibility across restarts, this must be done after all processes
     ! that can add soil or litter carbon cohorts.
     ce = first_elmt(land_tile_map, lnd%ls)
     do while (loop_over_tiles(ce,tile,l,k))
        call set_current_point(l,k) ! this is for debug output only
        if(.not.associated(tile%vegn)) cycle ! skip the rest of the loop body

        do ii = 1,N_LITTER_POOLS
           call cull_cohorts(tile%soil%litter(ii))
        enddo
        do ii=1,num_l
           call cull_cohorts(tile%soil%soil_organic_matter(ii))
        enddo
     enddo
  endif

  ! override with static vegetation
  if(day1/=day0) &
       call  read_static_vegn(lnd%time)

end subroutine update_vegn_slow



! ============================================================================
subroutine vegn_seed_transport_lm3()

  ! local vars
  type(land_tile_enum_type) :: ce
  type(land_tile_type), pointer :: tile
  integer ::  l ! current point index
  real :: total_seed_supply, total_seed_N_supply
  real :: total_seed_demand, total_seed_N_demand
  real :: f_supply, f_supply_N ! fraction of the supply that gets spent
  real :: f_demand, f_demand_N ! fraction of the demand that gets satisfied

  total_seed_supply = 0.0; total_seed_demand = 0.0; total_seed_N_supply = 0.0
  ce = first_elmt(land_tile_map, lnd%ls)
  do while (loop_over_tiles(ce,tile,l))
     if(.not.associated(tile%vegn)) cycle ! skip the rest of the loop body

     total_seed_supply = total_seed_supply + vegn_seed_supply(tile%vegn)*tile%frac*lnd%ug_area(l)
     total_seed_N_supply = total_seed_N_supply + vegn_seed_N_supply(tile%vegn)*tile%frac*lnd%ug_area(l)
     total_seed_demand = total_seed_demand + vegn_seed_demand(tile%vegn)*tile%frac*lnd%ug_area(l)
  enddo
  total_seed_N_demand = total_seed_demand/C2N_SEED
  ! sum totals globally
  call mpp_sum(total_seed_demand, pelist=lnd%pelist)
  call mpp_sum(total_seed_supply, pelist=lnd%pelist)
  call mpp_sum(total_seed_N_demand, pelist=lnd%pelist)
  call mpp_sum(total_seed_N_supply, pelist=lnd%pelist)
  ! if either demand or supply are zeros we don't need (or can't) transport anything
  if (total_seed_demand==0.or.total_seed_supply==0)then
     return
  end if

  ! calculate the fraction of the supply that is going to be used
  f_supply = MIN(total_seed_demand/total_seed_supply, 1.0)
  f_supply_N = MIN(total_seed_N_demand/total_seed_N_supply,1.0)
  ! calculate the fraction of the demand that is going to be satisfied
  f_demand = MIN(total_seed_supply/total_seed_demand, 1.0)
  f_demand_N = MIN(total_seed_N_supply/total_seed_N_demand,1.0)
  ! note that either f_supply or f_demand is 1; the mass conservation law in the
  ! following calculations is satisfied since
  ! f_demand*total_seed_demand - f_supply*total_seed_supply == 0

  ! redistribute part (or possibly all) of the supply to satisfy part (or possibly all)
  ! of the demand
  ce = first_elmt(land_tile_map)
  do while (loop_over_tiles(ce,tile))
     if(.not.associated(tile%vegn)) cycle ! skip the rest of the loop body

     call vegn_add_bliving(tile%vegn, &
          f_demand*vegn_seed_demand(tile%vegn)-f_supply*vegn_seed_supply(tile%vegn),&
          f_demand_N*vegn_seed_demand(tile%vegn)/C2N_seed-f_supply_N*vegn_seed_N_supply(tile%vegn))
  enddo
end subroutine vegn_seed_transport_lm3


! ============================================================================
! reads species table (if exists) from the input netcdf file and replaces
! species indices with the indices that correspond to the current set of
! species
subroutine read_remap_species(restart)
  type(land_restart_type) :: restart

  ! ---- local vars
  integer :: nsp ! number of input species
  integer :: spnames_id ! id of the species names table in the netcdf
  integer :: spnames_len(2)! sizes of the input species text array
  integer :: i, sp
  character(fm_field_name_len), allocatable :: spnames(:)
  character, allocatable :: text(:,:)
  integer, allocatable :: sptable(:) ! table for remapping
  type(land_tile_enum_type)     :: te,ce ! current and tail tile list elements
  type(land_tile_type), pointer :: tile  ! pointer to current tile

  if (.not.field_exists(restart, 'species_names')) then
     call error_mesg('vegn_init','variable "species_names" is not found in the restart, not remapping species',NOTE)
     return
     ! TODO: perhaps we still need to remap in this case, but using the prescribed
     ! list of LM3 species
  endif

  call get_text_data(restart, 'species_names', text)
  nsp = size(text,2)
  allocate(spnames(0:nsp-1), sptable(0:nsp-1))
  sptable(:) = -1
  do i = 0, nsp-1
     ! convert character array to strings
     call array2str(text(:,i+1),spnames(i))
     ! find corresponding species in the spdata array
     do sp = 0,size(spdata)-1
         if (trim(spdata(sp)%name)==trim(spnames(i))) then
            sptable(i) = sp
         endif
     enddo
  enddo

  ! Go through all tiles and initialize the cohorts that have not been initialized yet --
  ! this allows to read partial restarts. Also initialize accumulation counters to zero
  ! or the values from the restarts.
  ce = first_elmt(land_tile_map, ls=lnd%is)
  do while(loop_over_tiles(ce,tile))
     if (.not.associated(tile%vegn)) cycle
     if (tile%vegn%n_cohorts<=0) cycle ! skip uninitialized tiles
     do i = 1,tile%vegn%n_cohorts
        sp = tile%vegn%cohorts(i)%species
        if (sp<0.or.sp>=nsp) &
             call error_mesg('vegn_init','species index is outside of the bounds', FATAL)
        if (sptable(sp)<0) &
             call error_mesg('vegn_init','species "'//trim(spnames(sp))// &
                            '" from restart are not found in the model species parameter list',&
                            FATAL)
        tile%vegn%cohorts(i)%species = sptable(sp)
     enddo
  enddo
  deallocate(text, spnames, sptable)
end subroutine read_remap_species

! =====================================================================================
! given vegetation tile and cohort test function, returns the fraction of tile area
! occupied by the cohorts selected by the test function, as visible from above.

! This is for CMOR/CMIP diagnostic output, e.g. treeFrac

! Calculations are very similar to the layerfrac calculations in update_derived_vegn_data,
! except that we do not count internal gaps in the canopy.
function cohort_area_frac(vegn,test) result(frac); real :: frac
  type(vegn_tile_type), intent(in) :: vegn
  procedure(cohort_test_func) :: test ! returns TRUE if cohort fraction is counted

  integer :: n_layers ! number of layers in vegetation
  real, allocatable :: &
        layer_area(:), & ! area of _all_ cohorts in the layer
        c_area(:),     & ! area of _all suitable_ cohorts in the layer
        norm_area(:)     ! normalisation layer area
  real :: visible  ! fraction of layer visible from top
  integer :: k, l

  n_layers = maxval(vegn%cohorts(:)%layer)
  allocate(layer_area(n_layers),c_area(n_layers),norm_area(n_layers))

  layer_area(:) = 0; c_area(:) = 0.0
  associate(cc=>vegn%cohorts)
  do k = 1, vegn%n_cohorts
     l = cc(k)%layer
     layer_area(l) = layer_area(l) + cc(k)%crownarea*cc(k)%nindivs
     if (test(cc(k))) &
         c_area(l) = c_area(l)     + cc(k)%crownarea*cc(k)%nindivs
  enddo
  end associate

  if (allow_external_gaps) then
     norm_area(:) = max(1.0,layer_area(:))
  else
     ! if external gaps are not allowed, we stretch cohorts so that the entire layer area
     ! is occupied by the vegetation canopies
     norm_area(:) = layer_area(:)
     layer_area(:) = 1.0 ! to disallow gaps
  endif

  ! protect from zero layer area situation: this can happen if all cohorts die
  ! due to mortality or starvation. In this case n_layers is 1, of course.
  do k = 1,n_layers
     if (layer_area(k)<=0) layer_area(k) = 1.0
  enddo

  visible = 1.0; frac = 0.0
  do k = 1,n_layers
     frac = frac + visible*c_area(k)/norm_area(k)
     visible = visible*max(1.0-layer_area(k), 0.0)
  enddo

  deallocate(layer_area,c_area,norm_area)
end function cohort_area_frac

! ============================================================================
! converts character array to string
subroutine array2str(a,s)
  character,    intent(in)  :: a(:)
  character(*), intent(out) :: s

  integer :: i
  s = ' '
  do i = 1, min(len(s),size(a))
     s(i:i) = a(i)
  enddo
end subroutine array2str


! ============================================================================
! tile existence detector: returns a logical value indicating wether component
! model tile exists or not
logical function vegn_tile_exists(tile)
   type(land_tile_type), pointer :: tile
   vegn_tile_exists = associated(tile%vegn)
end function vegn_tile_exists

function any_vegn(cc) result(answer); logical answer
   type(vegn_cohort_type), intent(in) :: cc
   answer = .TRUE.
end function

function is_tree(cc) result(answer); logical answer
   type(vegn_cohort_type), intent(in) :: cc
   answer = (spdata(cc%species)%lifeform == FORM_WOODY)
end function

function is_grass(cc) result(answer); logical answer
   type(vegn_cohort_type), intent(in) :: cc
   answer = (spdata(cc%species)%lifeform == FORM_GRASS)
end function

function is_c3(cc) result(answer); logical answer
   type(vegn_cohort_type), intent(in) :: cc
   answer = (spdata(cc%species)%pt == PT_C3)
end function

function is_c4(cc) result(answer); logical answer
   type(vegn_cohort_type), intent(in) :: cc
   answer = (spdata(cc%species)%pt == PT_C4)
end function

function is_c3grass(cc) result(answer); logical answer
   type(vegn_cohort_type), intent(in) :: cc
   answer = (spdata(cc%species)%lifeform == FORM_GRASS.and.spdata(cc%species)%pt == PT_C3)
end function

function is_c4grass(cc) result(answer); logical answer
   type(vegn_cohort_type), intent(in) :: cc
   answer = (spdata(cc%species)%lifeform == FORM_GRASS.and.spdata(cc%species)%pt == PT_C4)
end function

! ============================================================================
! cohort accessor functions: given a pointer to cohort, return a pointer to a
! specific member of the cohort structure
#define DEFINE_VEGN_ACCESSOR_0D(xtype,x) subroutine vegn_ ## x ## _ptr(t,p);\
type(land_tile_type),pointer::t;xtype,pointer::p;p=>NULL();if(associated(t))then;if(associated(t%vegn))p=>t%vegn%x;endif;end subroutine

#define DEFINE_VEGN_ACCESSOR_1D(xtype,x) subroutine vegn_ ## x ## _ptr(t,i,p);\
type(land_tile_type),pointer::t;integer,intent(in)::i;xtype,pointer::p;p=>NULL();if(associated(t))then;if(associated(t%vegn))p=>t%vegn%x(i);endif;end subroutine

#define DEFINE_VEGN_ACCESSOR_2D(xtype,x) subroutine vegn_ ## x ## _ptr(t,i,j,p);\
type(land_tile_type),pointer::t;integer,intent(in)::i,j;xtype,pointer::p;p=>NULL();if(associated(t))then;if(associated(t%vegn))p=>t%vegn%x(i,j);endif;end subroutine

#define DEFINE_COHORT_ACCESSOR(xtype,x) subroutine cohort_ ## x ## _ptr(c,p);\
type(vegn_cohort_type),pointer::c;xtype,pointer::p;p=>NULL();if(associated(c))p=>c%x;end subroutine

#define DEFINE_COHORT_COMPONENT_ACCESSOR(xtype,component,x) subroutine cohort_ ## x ## _ptr(c,p);\
type(vegn_cohort_type),pointer::c;xtype,pointer::p;p=>NULL();if(associated(c))p=>c%component%x;end subroutine

DEFINE_VEGN_ACCESSOR_0D(integer,landuse)
DEFINE_VEGN_ACCESSOR_0D(real,age)
DEFINE_VEGN_ACCESSOR_0D(real,fsc_pool_ag)
DEFINE_VEGN_ACCESSOR_0D(real,fsc_rate_ag)
DEFINE_VEGN_ACCESSOR_0D(real,ssc_pool_ag)
DEFINE_VEGN_ACCESSOR_0D(real,ssc_rate_ag)
DEFINE_VEGN_ACCESSOR_0D(real,fsc_pool_bg)
DEFINE_VEGN_ACCESSOR_0D(real,fsc_rate_bg)
DEFINE_VEGN_ACCESSOR_0D(real,ssc_pool_bg)
DEFINE_VEGN_ACCESSOR_0D(real,ssc_rate_bg)

DEFINE_VEGN_ACCESSOR_0D(real,fsn_pool_bg)
DEFINE_VEGN_ACCESSOR_0D(real,fsn_rate_bg)
DEFINE_VEGN_ACCESSOR_0D(real,ssn_pool_bg)
DEFINE_VEGN_ACCESSOR_0D(real,ssn_rate_bg)

DEFINE_VEGN_ACCESSOR_2D(real,litter_buff_C)
DEFINE_VEGN_ACCESSOR_2D(real,litter_buff_N)
DEFINE_VEGN_ACCESSOR_2D(real,litter_rate_C)
DEFINE_VEGN_ACCESSOR_2D(real,litter_rate_N)

DEFINE_VEGN_ACCESSOR_0D(real,tc_av)
DEFINE_VEGN_ACCESSOR_0D(real,theta_av_phen)
DEFINE_VEGN_ACCESSOR_0D(real,theta_av_fire)
DEFINE_VEGN_ACCESSOR_0D(real,psist_av)
DEFINE_VEGN_ACCESSOR_0D(real,tsoil_av)
DEFINE_VEGN_ACCESSOR_0D(real,precip_av)
DEFINE_VEGN_ACCESSOR_0D(real,fuel)
DEFINE_VEGN_ACCESSOR_0D(real,lambda)
DEFINE_VEGN_ACCESSOR_0D(real,t_ann)
DEFINE_VEGN_ACCESSOR_0D(real,p_ann)
DEFINE_VEGN_ACCESSOR_0D(real,t_cold)
DEFINE_VEGN_ACCESSOR_0D(real,ncm)
DEFINE_VEGN_ACCESSOR_0D(real,t_ann_acm)
DEFINE_VEGN_ACCESSOR_0D(real,p_ann_acm)
DEFINE_VEGN_ACCESSOR_0D(real,t_cold_acm)
DEFINE_VEGN_ACCESSOR_0D(real,ncm_acm)
DEFINE_VEGN_ACCESSOR_0D(real,treeline_T_accum)
DEFINE_VEGN_ACCESSOR_0D(real,treeline_N_accum)

DEFINE_VEGN_ACCESSOR_0D(real,tc_pheno)
DEFINE_VEGN_ACCESSOR_0D(real,csmoke_pool)
DEFINE_VEGN_ACCESSOR_0D(real,csmoke_rate)
DEFINE_VEGN_ACCESSOR_0D(real,drop_wl)
DEFINE_VEGN_ACCESSOR_0D(real,drop_ws)
DEFINE_VEGN_ACCESSOR_0D(real,drop_hl)
DEFINE_VEGN_ACCESSOR_0D(real,drop_hs)
DEFINE_VEGN_ACCESSOR_0D(real,nsmoke_pool)

DEFINE_VEGN_ACCESSOR_1D(real,harv_pool_C)
DEFINE_VEGN_ACCESSOR_1D(real,harv_pool_N)
DEFINE_VEGN_ACCESSOR_1D(real,harv_rate_C)
DEFINE_VEGN_ACCESSOR_1D(real,drop_seed_C)
DEFINE_VEGN_ACCESSOR_1D(real,drop_seed_N)

DEFINE_COHORT_ACCESSOR(real,Tv)
DEFINE_COHORT_ACCESSOR(real,Wl)
DEFINE_COHORT_ACCESSOR(real,Ws)

DEFINE_COHORT_ACCESSOR(integer,species)
DEFINE_COHORT_ACCESSOR(real,bl)
DEFINE_COHORT_ACCESSOR(real,br)
DEFINE_COHORT_ACCESSOR(real,blv)
DEFINE_COHORT_ACCESSOR(real,bsw)
DEFINE_COHORT_ACCESSOR(real,bwood)
DEFINE_COHORT_ACCESSOR(real,nsc)
DEFINE_COHORT_ACCESSOR(real,bseed)
DEFINE_COHORT_ACCESSOR(real,bsw_max)
DEFINE_COHORT_ACCESSOR(real,bl_max)
DEFINE_COHORT_ACCESSOR(real,br_max)
DEFINE_COHORT_ACCESSOR(real,dbh)
DEFINE_COHORT_ACCESSOR(real,crownarea)
DEFINE_COHORT_ACCESSOR(real,bliving)
DEFINE_COHORT_ACCESSOR(real,nindivs)
DEFINE_COHORT_ACCESSOR(integer,layer)
DEFINE_COHORT_ACCESSOR(integer,firstlayer)
DEFINE_COHORT_ACCESSOR(integer,status)
DEFINE_COHORT_ACCESSOR(real,leaf_age)
DEFINE_COHORT_ACCESSOR(real,age)
DEFINE_COHORT_ACCESSOR(real,npp_previous_day)
DEFINE_COHORT_ACCESSOR(real,growth_previous_day)
DEFINE_COHORT_ACCESSOR(real,growth_previous_day_tmp)
DEFINE_COHORT_ACCESSOR(real,BM_ys)
DEFINE_COHORT_ACCESSOR(real,DBH_ys)
DEFINE_COHORT_ACCESSOR(real,topyear)
DEFINE_COHORT_ACCESSOR(real,gdd)
DEFINE_COHORT_ACCESSOR(real,height)

DEFINE_COHORT_ACCESSOR(real,scav_myc_C_reservoir)
DEFINE_COHORT_ACCESSOR(real,scav_myc_N_reservoir)
DEFINE_COHORT_ACCESSOR(real,mine_myc_C_reservoir)
DEFINE_COHORT_ACCESSOR(real,mine_myc_N_reservoir)
DEFINE_COHORT_ACCESSOR(real,N_fixer_C_reservoir)
DEFINE_COHORT_ACCESSOR(real,N_fixer_N_reservoir)
DEFINE_COHORT_ACCESSOR(real,myc_scavenger_biomass_C)
DEFINE_COHORT_ACCESSOR(real,myc_scavenger_biomass_N)
DEFINE_COHORT_ACCESSOR(real,myc_miner_biomass_C)
DEFINE_COHORT_ACCESSOR(real,myc_miner_biomass_N)
DEFINE_COHORT_ACCESSOR(real,N_fixer_biomass_C)
DEFINE_COHORT_ACCESSOR(real,N_fixer_biomass_N)
DEFINE_COHORT_ACCESSOR(real,stored_N)
DEFINE_COHORT_ACCESSOR(real,wood_N)
DEFINE_COHORT_ACCESSOR(real,sapwood_N)
DEFINE_COHORT_ACCESSOR(real,leaf_N)
DEFINE_COHORT_ACCESSOR(real,seed_N)
DEFINE_COHORT_ACCESSOR(real,root_N)
DEFINE_COHORT_ACCESSOR(real,total_N)
DEFINE_COHORT_ACCESSOR(real,nitrogen_stress)
DEFINE_COHORT_ACCESSOR(real,myc_scav_marginal_gain_smoothed)
DEFINE_COHORT_ACCESSOR(real,myc_mine_marginal_gain_smoothed)
DEFINE_COHORT_ACCESSOR(real,N_fix_marginal_gain_smoothed)
DEFINE_COHORT_ACCESSOR(real,rhiz_exud_marginal_gain_smoothed)

DEFINE_COHORT_ACCESSOR(real,max_monthly_scav_alloc)
DEFINE_COHORT_ACCESSOR(real,max_monthly_mine_alloc)
DEFINE_COHORT_ACCESSOR(real,max_monthly_Nfix_alloc)
DEFINE_COHORT_ACCESSOR(real,Nfix_alloc_accum)
DEFINE_COHORT_ACCESSOR(real,mine_alloc_accum)
DEFINE_COHORT_ACCESSOR(real,scav_alloc_accum)
DEFINE_COHORT_ACCESSOR(real,max_Nfix_allocation)
DEFINE_COHORT_ACCESSOR(real,max_mine_allocation)
DEFINE_COHORT_ACCESSOR(real,max_scav_allocation)
DEFINE_COHORT_ACCESSOR(real,nitrogen_stress_smoothed)

! wolf
DEFINE_COHORT_ACCESSOR(real,psi_r)
DEFINE_COHORT_ACCESSOR(real,psi_x)
DEFINE_COHORT_ACCESSOR(real,psi_l)
DEFINE_COHORT_ACCESSOR(real,Kxa)
DEFINE_COHORT_ACCESSOR(real,Kla)

! ens for branch sapwood
DEFINE_COHORT_ACCESSOR(real,brsw)

end module vegetation_mod
