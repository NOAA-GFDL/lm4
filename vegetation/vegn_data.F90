module vegn_data_mod

#ifdef INTERNAL_FILE_NML
use mpp_mod, only: input_nml_file
#else
use fms_mod, only: open_namelist_file
#endif

use constants_mod, only : PI, TFREEZE
use fms_mod, only : &
     file_exist, check_nml_error, &
     close_file, stdlog, stdout, string, lowercase, error_mesg, NOTE, FATAL
use field_manager_mod, only: MODEL_LAND, fm_field_name_len, fm_string_len, &
     fm_path_name_len, fm_type_name_len, fm_dump_list, fm_get_length, &
     fm_get_current_list, fm_change_list, fm_list_iter_type, fm_init_loop, fm_loop_over_list
use fm_util_mod, only : fm_util_get_real, fm_util_get_logical, fm_util_get_string

use land_constants_mod, only : NBANDS, BAND_VIS, BAND_NIR, MPa_per_m
use land_data_mod, only : log_version
use land_tile_selectors_mod, only : &
     tile_selector_type, SEL_VEGN, register_tile_selector
use table_printer_mod

implicit none
private

! ==== public interfaces =====================================================
! ---- public constants
integer, public, parameter :: LU_SEL_TAG = 1 ! tag for the land use selectors
integer, public, parameter :: SP_SEL_TAG = 2 ! tag for the species selectors
integer, public, parameter :: NG_SEL_TAG = 3 ! tag for natural grass selector
  ! by "natural" it means non-human-maintained, so secondary vegetation
  ! grassland will be included.

integer, public, parameter :: & ! life form of the plant
 FORM_GRASS = 0, &
 FORM_WOODY = 1
 ! in future, possibly add mosses...

integer, public, parameter :: N_LM3_SPECIES = 5, & ! number of species
 SP_C4GRASS   = 0, & ! c4 grass
 SP_C3GRASS   = 1, & ! c3 grass
 SP_TEMPDEC   = 2, & ! temperate deciduous
 SP_TROPICAL  = 3, & ! non-grass tropical
 SP_EVERGR    = 4    ! non-grass evergreen
character(len=12), parameter :: lm3_species_name(0:N_LM3_SPECIES-1) = &
    (/'c4grass  ',  'c3grass  ' ,  'tempdec  ', 'tropical ','evergreen'/)

integer, public, parameter :: n_dim_vegn_types = 9
integer, public, parameter :: MSPECIES = N_LM3_SPECIES+n_dim_vegn_types-1

integer, public, parameter :: & ! physiology types
 PT_C3        = 0, &
 PT_C4        = 1

integer, public, parameter :: & ! phenology type
 PHEN_DECIDUOUS = 0, &
 PHEN_EVERGREEN = 1

integer, public, parameter :: & ! allometry type
 ALLOM_EW     = 0, & ! Ensheng's original
 ALLOM_EW1    = 1, & ! Ensheng's "alternative"
 ALLOM_HML    = 2    ! Helena's allometry

integer, public, parameter :: & ! status of leaves
 LEAF_ON      = 0, &  ! leaves are displayed
 LEAF_OFF     = 5     ! leaves are dropped

integer, public, parameter :: & ! land use types
 N_LU_TYPES = 6, & ! number of different land use types
 LU_PAST    = 1, & ! pasture
 LU_CROP    = 2, & ! crops
 LU_NTRL    = 3, & ! natural vegetation
 LU_SCND    = 4, & ! secondary vegetation
 LU_URBN    = 5, & ! urban
 LU_RANGE   = 6, & ! rangeland
 LU_PSL     = 1001 ! primary and secondary land, for LUMIP

character(len=5), public, parameter  :: &
     landuse_name (N_LU_TYPES) = (/ 'past ','crop ','ntrl ','scnd ', 'urbn ', 'range'/)
character(len=32), public, parameter :: &
     landuse_longname (N_LU_TYPES) = (/ 'pasture  ', 'crop     ', 'natural  ', 'secondary', 'urban    ','rangeland'/)

integer, public, parameter :: & ! harvesting pools parameters
 N_HARV_POOLS        = 6, & ! number of harvesting pools
 HARV_POOL_PAST      = 1, &
 HARV_POOL_CROP      = 2, &
 HARV_POOL_CLEARED   = 3, &
 HARV_POOL_WOOD_FAST = 4, &
 HARV_POOL_WOOD_MED  = 5, &
 HARV_POOL_WOOD_SLOW = 6
character(len=9), public :: HARV_POOL_NAMES(N_HARV_POOLS)
data HARV_POOL_NAMES &
 / 'past', 'crop', 'cleared', 'wood_fast', 'wood_med', 'wood_slow' /

real, public, parameter :: C2B = 2.0  ! carbon to biomass conversion factor

real, public, parameter :: BSEED = 5e-5 ! seed density for supply/demand calculations, kg C/m2
real, public, parameter :: C2N_SEED = 50 ! seed C:N ratio

! ---- public types
public :: spec_data_type

! ---- public data
integer, public :: nspecies ! total number of species
public :: &
    vegn_to_use,  input_cover_types, &
    mcv_min, mcv_lai, &
    use_bucket, use_mcm_masking, vegn_index_constant, &
    critical_root_density, &
    spdata, &
    min_cosz, &
    agf_bs, K1,K2, &
    tau_drip_l, tau_drip_s, & ! canopy water and snow residence times, for drip calculations
    GR_factor, tg_c3_thresh, T_cold_tropical, tg_c4_thresh, &
    fsc_pool_spending_time, ssc_pool_spending_time, harvest_spending_time, &
    T_transp_min, soil_carbon_depth_scale, &
    cold_month_threshold, scnd_biomass_bins, &
    treeline_thresh_T, treeline_base_T, treeline_season_length, &
    phen_ev1, phen_ev2, cmc_eps, &
    b0_growth, tau_seed, understory_lai_factor, min_lai, min_cohort_nindivs, &
    DBH_mort, A_mort, B_mort, cold_mort, treeline_mort, nsc_starv_frac, &
    DBH_merge_rel, DBH_merge_abs, NSC_merge_rel, &

    mycorrhizal_turnover_time, &
    myc_scav_C_efficiency, myc_mine_C_efficiency, &
    N_fixer_turnover_time, N_fixer_C_efficiency, &
    c2n_N_fixer, N_limits_live_biomass, &
    excess_stored_N_leakage_rate, min_N_stress, calc_SLA_from_lifespan, &
    et_myc, smooth_N_uptake_C_allocation, N_fix_Tdep_Houlton

logical, public :: do_ppa = .FALSE.
logical, public :: nat_mortality_splits_tiles = .FALSE. ! if true, natural mortality
    ! creates disturbed tiles

! ---- public subroutine
public :: read_vegn_data_namelist
! ==== end of public interfaces ==============================================

! ==== constants =============================================================
character(*), parameter :: module_name = 'vegn_data_mod'
character(*), parameter :: data_error_header = 'VEGETATION DATA FATAL ERROR'
#include "../shared/version_variable.inc"
real, parameter :: TWOTHIRDS  = 2.0/3.0


! ==== types ================================================================
type spec_data_type
  character(fm_field_name_len) :: name     = '' ! name of the species
  character(fm_string_len)     :: longname = '' ! longname of the species
  real    :: treefall_disturbance_rate = 0.025
  logical :: mortality_kills_balive    = .false.! if true, then bl, blv, and br are affected by natural mortality
  integer :: pt = PT_C3 ! photosynthetic physiology of species
  integer :: phent = PHEN_DECIDUOUS ! type of phenology
  integer :: allomt = ALLOM_HML ! type of allometry
  integer :: lifeform = FORM_WOODY ! vegetation lifeform

  real    :: c1 = 0.4807692 ! unitless, coefficient for living biomass allocation
  real    :: c2 = 0.4004486 ! 1/m, coefficient for living biomass allocation
  real    :: c3 = 0.4423077 ! unitless, coefficient for calculation of sapwood biomass
                ! fraction times sapwood retirement rate

  ! decay rates of plant carbon pools, 1/yr
  real    :: alpha_leaf=1.0, alpha_root=1.0, alpha_wood=1.2e-2
  real    :: branch_loss_height = 2.0 ! trees do not lose branches till they reach that height, m
  ! respiration rates of plant carbon pools
  real    :: beta_sapwood=0.0, beta_root=1.25, beta_vleaf=0.0

  real    :: dfr         = 5.8   ! fine root diameter ? or parameter relating diameter of fine roots to resistance
  ! the following two parameters are used in the Darcy-law calculations of water supply
  real    :: srl        = 24.4e3 ! specific root length, m/(kg C)
  real    :: root_r     = 2.9e-4 ! radius of the fine roots, m
  real    :: root_perm  = 5.0e-7 ! fine root membrane permeability per unit area, kg/(m3 s)

  real    :: LMA        = 0.036  ! leaf mass per unit area, kg C/m2
  real    :: LMA_understory_factor = 0.5 ! Resp multiplier for understory
  real    :: leaf_size  = 0.04   ! characteristic leaf size

  real    :: alpha_phot = 0.06   ! photosynthesis efficiency
  real    :: m_cond     = 9.0    ! factor of stomatal conductance
  real    :: Vmax       = 2e-5   ! max rubisco rate
  real    :: gamma_resp = 0.02   !
  real    :: wet_leaf_dreg = 0.3 ! wet leaf photosynthesis down-regulation
  real    :: leaf_age_onset = 100.0, leaf_age_tau = 150.0
  real    :: Vmax_understory_factor = 0.5 ! Vmax multiplier for understory
  real    :: Resp_understory_factor = 0.5 ! Resp multiplier for understory

  ! radiation parameters for 2 bands, VIS and NIR
  real    :: leaf_refl (NBANDS) = (/ 0.10, 0.50/) ! reflectance of leaf
  real    :: leaf_tran (NBANDS) = (/ 0.05, 0.25/) ! transmittance of leaf
  real    :: leaf_emis          = 1.00            ! emissivity of leaf

  ! parameters of leaf angle distribution; see also Bonan, NCAR/TN-417+STR (LSM
  ! 1.0 technical description), p.18
  real    :: ksi    =  0.0 ! departure of leaf angles from a random distribution
  real    :: phi1   ! leaf distribution parameter
  real    :: phi2   ! leaf distribution parameter
  real    :: mu_bar ! average inverse diffuse optical depth per unit leaf are

  ! canopy intercepted water parameters
  real    :: cmc_lai = 0.02 ! max amount of liquid water on vegetation, kg/(m2 of leaf)
  real    :: cmc_pow = TWOTHIRDS ! power of wet fraction dependance on amount of canopy water
  real    :: csc_lai = 0.2  ! max amount of snow on vegetation, kg/(m2 of leaf)
  real    :: csc_pow = TWOTHIRDS ! power of snow-covered fraction dependance on amount of canopy snow

  real    :: internal_gap_frac = 0.1 ! fraction of internal gaps in the canopy
  ! "internal" gaps are the gaps that are created within the canopy by the
  ! branch fall processes.

  real    :: fuel_intensity = 0.002


  real    :: tc_crit   = 283.16 ! critical temperature for leaf drop, degK
  real    :: gdd_crit  = 360.0
  real    :: gdd_base_T =  273.16 ! base temperature for GDD calculations, degK
  ! critical soil-water-stress index, used in place of fact_crit_phen and
  ! cnst_crit_phen. It is used if and only if it's value is greater than 0
  real    :: psi_stress_crit_phen = 0.0
  real    :: fact_crit_phen = 0.0, cnst_crit_phen = 0.15 ! wilting factor and offset to
    ! get critical value for leaf drop -- only one is non-zero at any time
  real    :: leaf_C_retrans_frac = 0.5 ! fraction of the leaf biomass re-translocated before leaf drop
  real    :: root_C_retrans_frac = 0.5 ! fraction of the root biomass re-translocated before senescence
  real    :: fact_crit_fire = 0.0, cnst_crit_fire = 0.15 ! wilting factor and offset to
    ! get critical value for fire -- only one is non-zero at the time

  real    :: smoke_fraction   = 0.6 ! fraction of carbon lost as smoke during fires

  ! data from LM3W, temporarily here
  real    :: dat_height       = 6.6
  real    :: dat_lai          = 0.0
  real    :: dat_root_density = 4.2
  real    :: dat_root_zeta    = 0.28909
  real    :: dat_rs_min       = 131.0
  real    :: dat_snow_crit    = 0.0333
  !  for PPA, Weng, 7/25/2011
  real    :: alphaHT = 20.0,  thetaHT = 0.5, gammaHT = 0.6841742 ! height allometry parameters
  real    :: alphaCA = 30.0,  thetaCA = 1.5 ! crown area allometry parameters
  real    :: alphaBM = 0.559, thetaBM = 2.5 ! biomass allometry parameters
  real    :: alphaCSASW    = 2.25e-2, thetaCSASW = 1.5 !
  logical :: limit_tussock_R = .FALSE. ! if true, impose limit on tussock crown radius and area
  real    :: tussock_Ra=0.1, tussock_Rb=0.3 ! max diameter of tussock crown is calculated
     ! as tussock_Ra + tussock_Rb*height
  real    :: layer_height_factor = 1.0  ! scaling of height for cohort relayering.
     ! Grasses bend under the wind, exposing the tree seedlings to light even if the
     ! seedlings are shorter; therefore we can scale the effective height of the grasses
     ! down for the purposes of assigning layers to the vegetation cohorts.

  real    :: maturalage    = 1.0    ! the age that can reproduce
  real    :: v_seed        = 0.1    ! fraction of G_SF to G_F
  logical :: reproduces_in_understory = .FALSE. ! if true, plant can reproduce in the understory,
                                    ! otherwise only in top layer

  ! seed dispersal and transport: obviously the fraction of seed dispersed to other tiles
  ! and grid cells would depend on the size of the tiles grid cells. In future, we should
  ! calculate those fractions, given dispersal radius of seeds for given species, and
  ! geometry of tiles and grid cells. Characteristic size of tiles may need to be assumed.
  real    :: frac_seed_dispersed   = 1e-2 ! fraction of seeds dispersed inside the grid cell (between the tiles)
  real    :: frac_seed_transported = 1e-4 ! fraction of seeds dispersed outside of grid cell (transported to neighbors)
  ! total fraction of seeds leaving tile due to dispersion: frac_seed_dispersed+frac_seed_transported

  real    :: seedling_height   = 0.1 ! height of the seedlings, m
  real    :: seedling_nsc_frac = 3.0 ! initial seedling NSC, as fraction of bl_max (typically > 1)
  real    :: prob_g = 0.45, prob_e = 0.3 ! germination and establishment probabilities
  real    :: mortrate_d_c  = 0.05   ! daily mortality rate in canopy
  real    :: mortrate_d_u  = 0.2    ! daily mortality rate in understory
  real    :: Tmin_mort     = 0.0    ! cold mortality threshold, degK; default zero value turns it off
  logical :: mortality_kills_seeds = .FALSE. ! if true, mortality kills seeds; otherwise seeds survive and used for reproduction later.
  real    :: rho_wood      = 250.0  ! woody density, kg C m-3 wood
  real    :: taperfactor   = 0.9
  real    :: LAImax        = 3.0    ! max. LAI
  real    :: f_NSC         = 0.1    ! max rate of NSC spending on leaf and fine root growth, 1/day, eq.A7 in Weng et al. 2015
  real    :: f_WF          = 0.2    ! max rate of NSC excess spending on sapwood and fruit growth, 1/day, eq. A12 Weng et al 2015
  logical :: use_light_saber = .FALSE. ! if TRUE, then the leaves at the bottom
    ! of the canopy that cannot support themselves by photosynthesis are mercilessly
    ! cut off.
  real    :: light_saber_Hmin = 0.0 ! height below which light_saber is not applied, m.
  real    :: phiRL         = 2.69   ! ratio of fine root to leaf area
  real    :: phiCSA        = 2.5e-4 ! ratio of sapwood CSA to target leaf area
  real    :: SRA           = 44.45982 ! specific fine root area, m2/kg C
  !  for PPA, IMC, 1/8/2017
  real    :: growth_resp   = 0.333  ! fraction of NPP lost as growth respiration
  ! target NSC is calculated as blending between bl_max*NSC2targetbl and bl_max*NSC2targetbl0
  ! between DBH=0 and DBH=NSC2targetbl_dbh, provided the latter is above 0
  real    :: NSC2targetbl  = 4.0    ! ratio of NSC to target biomass of leaves
  real    :: NSC2targetbl0 = 1.5    ! ratio of NSC to target biomass of leaves for zero-size seedlings
  real    :: NSC2targetbl_dbh = -1.0  ! blending diameter for NSC target calculations
  real    :: T_dorm        = TFREEZE  ! dormancy temperature threshold, degK

  real    :: tracer_cuticular_cond = 0.0 ! cuticular conductance for all tracers, m/s

  ! for Kok effect, ppg
  real    :: inib_factor = 0.0 ! default rate of inhibition for leaf respiration
  real    :: light_kok = 0.00004 ! Light intensity above which light inhibition occurs, umol/m2/s

  ! Temperatures for photosynthesis and respiration response, degC
  real    :: ToptP = 35.0
  real    :: TminP = 5.0
  real    :: TmaxP = 45.0
  real    :: ToptR = 42.0
  real    :: TminR = 5.0
  real    :: TmaxR = 45.0

  real    :: Ea_ko = 14.5! Activation Energy for ko
  real    :: Hd_ko = 200.0! Inactivation Energy for ko
  real    :: Ea_kc = 80.5! Activation Energy for kc
  real    :: Hd_kc = 200.0! Inactivation Energy for kc
  real    :: Ea_vm = 48.75! Activation Energy for vm
  real    :: Hd_vm = 200.0! Inactivation Energy for vm
  real    :: Ea_gam = 24.5! Activation Energy for gamma star
  real    :: Hd_gam = 200.0! Inactivation Energy for gamma star
  real    :: Ea_resp = 46.39! Activation Energy for resp
  real    :: Hd_resp = 200.0! Inactivation Energy for resp

  ! for hydraulics, wolf
  real    :: Kxam=0.0, Klam=0.0 ! Conductivity, max, per tissue area: units kg/m2 tissue/s/MPa
  real    :: dx=0.0, dl=0.0     ! Breakpoint of Weibull function, MPa
  real    :: cx=1.0, cl=1.0	    ! Exponent of Weibull function, unitless
  real    :: psi_tlp=0.0        ! psi at turgor loss point

  logical :: do_N_mining_strategy = .TRUE.
  logical :: do_N_scavenging_strategy = .TRUE.
  logical :: do_N_fixation_strategy = .TRUE.
  real    :: N_fixation_rate = 0.1 ! N fixation rate per unit fixer biomass (kgN/kg fixer C/year)
  real    :: kM_myc_growth = 1e-3  ! Half-saturation constant for mycorrhizal growth (on C reservoir)
  real    :: myc_growth_rate = 1.0        ! Max rate at which mycorrhizae gain biomass (kgC/m2/year)
  real    :: myc_N_to_plant_rate = 100.0    ! Rate at which plants send N from reservoir to plant (Fraction per year)
  real    :: root_NH4_uptake_rate = 0.1   ! Max rate of root NH4 uptake (kgN/m3/year)
  real    :: root_NO3_uptake_rate = 0.1   ! Max rate of root NO3 uptake
  real    :: k_ammonium_root_uptake = 3e-2   ! Half-saturation for root NH4 uptake (kgN/m3)
  real    :: k_nitrate_root_uptake = 3e-2    ! Half-saturation for root NO3 uptake (kgN/m3)

  real    :: root_exudate_frac = 0.0 ! fraction of NPP that ends up in root exudates
  real    :: root_exudate_N_frac = 0.0 ! N fraction of root exudates. See e.g. Drake et al 2013
  logical :: dynamic_root_exudation = .FALSE. ! Whether to dynamically determine root exudation rate from plant N limitation

  real    :: fsc_liv    = 0.8 ! Species-specific fsc_liv, separate for backwards compatibility
  real    :: fsc_froot  = 0.3
  real    :: fsc_wood   = 0.2

  real    :: leaf_live_c2n  = 30.0    ! C:N ratio of live leaves.
  real    :: froot_live_c2n = 50.0    ! C:N ratio of live fine roots.
  real    :: sapwood_c2n    = 50.0    ! C:N ratio of sapwood.
  real    :: wood_c2n       = 500.0   ! x2z Wiki  http://en.wikipedia.org/wiki/Carbon-to-nitrogen_ratio
  real    :: c2n_mycorrhizae= 10.0    ! C:N ratio of mycorrhizal biomass
  real    :: seed_c2n       = 30.0    ! C:N ratio of seeds. Note, must contain enough N to support initial sapling

  real    :: leaf_N_retrans_frac = 0.5 ! Fraction of leaf N retranslocated before leaf drop
  real    :: root_N_retrans_frac = 0.0 ! Fraction of fine root N retranslocated before senescence

  real    :: branch_wood_frac = 0.1525 ! fraction of total wood biomass in branches,
                                       ! corresponds to 0.18 of trunk (bole) biomass
                                       ! estimated by Isa from the observations
  real    :: N_stress_root_factor = 0.02 ! Amount of sapwood C that goes to roots instead as N stress increases
  real    :: max_n_stress_for_seed_production = 2.0 ! if N stress is higher than this treshold, seeds are not produced
  real    :: tau_nsc_exudate = 1.0 ! e-folding time for nsc spending to exudates, yr

  real    :: tau_smooth_marginal_gain = 0.0
  real    :: tau_smooth_alloc         = 0.0
  real    :: alloc_allowed_over_limit = 10.0
  real    :: tau_smooth_Nstress       = 0.0

  ! SSR fire-related parameters; default values are for tropical trees in his code
  real    :: ROS_max   = 0.22
  real    :: fire_duration = 86400.0 ! average fire duration, s
  ! combustion completeness; unitless
  real    :: CC_leaf   = 0.70
  real    :: CC_stem   = 0.15
  real    :: CC_litter = 0.15
  ! Post-fire mortality
  real    :: fireMort_leaf = 0.70
  real    :: fireMort_stem = 0.65
  real    :: fireMort_root = 0.10
  ! dsward parameters
  real    :: EF_CO2 = 1643.0, EF_CO = 63.0
  real    :: F_scorch_height = 0.0 ! Fuel bulk density parameter from Thonicke et al. (2010)
             ! Original values:
             ! SP_C4GRASS  : 0.
             ! SP_C3GRASS  : 0.
             ! SP_TEMPDEC  : 0.094
             ! SP_TROPICAL : 0.1487
             ! SP_EVERGR   : 0.11
end type

! ==== module data ===========================================================
integer :: idata ! iterators used in data initialization statements

! ---- namelist --------------------------------------------------------------
type(spec_data_type), allocatable :: spdata(:)

logical :: use_bucket = .false.
logical :: use_mcm_masking = .false.
real    :: mcv_min = 5.   * 4218.
real    :: mcv_lai = 0.15 * 4218.

! ---- remainder are used only for cold start
character(len=16):: vegn_to_use     = 'single-tile'
       ! 'multi-tile' for tiled vegetation
       ! 'single-tile' for geographically varying vegetation with single type per
       !     model grid cell
       ! 'uniform' for global constant vegetation, e.g., to reproduce MCM
integer :: vegn_index_constant   = 1         ! index of global constant vegn,
                                             ! used when vegn_to_use is 'uniform'
real    :: critical_root_density = 0.125

integer, dimension(1:MSPECIES) :: &
 input_cover_types=(/          -1,   -1,   -1,   -1, &
                          1,    2,    3,    4,    5,    6,    7,    8,    9/)

real :: min_cosz = 0.01 ! minimum allowed value of cosz for vegetation radiation
   ! properties calculations.
   ! It probably doesn't make sense to set it any less than the default value, because the angular
   ! diameter of the sun is about 0.01 radian (0.5 degree), so the spread of the direct radiation
   ! zenith angles is about this. Besides, the sub-grid variations of land surface slope are
   ! probably even larger that that.
real :: soil_carbon_depth_scale = 0.2   ! depth of active soil for carbon decomposition
real :: cold_month_threshold    = 283.0 ! monthly temperature threshold for calculations of number of cold months
real :: agf_bs         = 0.8 ! ratio of above ground stem to total stem
real :: K1 = 10.0, K2 = 0.05 ! soil decomposition parameters
real :: tau_drip_l     = 21600.0 ! canopy water residence time, for drip calculations
real :: tau_drip_s     = 86400.0 ! canopy snow residence time, for drip calculations
real :: GR_factor = 0.33 ! growth respiration factor

real :: tg_c3_thresh = 1.5 ! threshold biomass between tree and grass for C3 plants
real :: tg_c4_thresh = 2.0 ! threshold biomass between tree and grass for C4 plants
real :: T_cold_tropical = 278.16 ! Minimum T_cold that makes vegetation tropical (not deciduous)
real :: fsc_pool_spending_time = 1.0 ! time (yrs) during which intermediate pool of
                  ! fast soil carbon is entirely converted to the fast soil carbon
real :: ssc_pool_spending_time = 1.0 ! time (yrs) during which intermediate pool of
                  ! slow soil carbon is entirely converted to the slow soil carbon
real :: harvest_spending_time(N_HARV_POOLS) = &
     (/1.0, 1.0, 1.0, 1.0, 10.0, 100.0/)
     ! time (yrs) during which intermediate pool of harvested carbon is completely
     ! released to the atmosphere.
     ! NOTE: a year in the above *_spending_time definitions is exactly 365*86400 seconds
real :: T_transp_min = 0.0 ! lowest temperature at which transpiration is enabled
                           ! 0 means no limit, lm3v value is 268.0
! Ensheng's growth parameters:
real :: b0_growth   = 0.02   ! min biomass for growth formula, kgC/indiv
real :: tau_seed    = 0.5708 ! characteristic time of nsc spending on seeds, year

! reduction of bl_max and br_max for the understory vegetation, unitless
real, protected :: understory_lai_factor = 0.25

real, protected :: min_lai = 1e-5 ! minimum lai: if leaf fall brings LAI
    ! below this threshold, bl is set to zero

real, protected :: min_cohort_nindivs = 1e-12 ! minimum allowed cohort density, individuals per m2

! boundaries of wood biomass bins for secondary veg. (kg C/m2); used to decide
! whether secondary vegetation tiles can be merged or not. MUST BE IN ASCENDING
! ORDER.
real  :: scnd_biomass_bins(10) &
     = (/ 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 1000.0 /)
real :: phen_ev1 = 0.5, phen_ev2 = 0.9 ! thresholds for evergreen/deciduous
      ! differentiation (see phenology_type in cohort.F90)
real :: cmc_eps = 0.01 ! value of w/w_max for transition to linear function;
                       ! the same value is used for liquid and snow

! tree line parameters from Paulsen and Korner (2014) A climate-based model to predict
! potential treeline position around the globe, Alpine Botany, 124, Issue 1, pp 1–12
real, protected :: treeline_base_T   = Tfreeze + 0.9 ! base temperature for growing season calculations for tree line, degK
real, protected :: treeline_thresh_T = Tfreeze + 6.4 ! threshold temperature for tree line, degK
real, protected :: treeline_season_length = 94.0     ! minimum season length for the trees to survive
! Weng, 7/25/2011
! for understory mortality rate is calculated as:
! deathrate = mortrate_d_u * ( 1 + A * exp(B*(DBH_mort-DBH))/(1 + exp(B*(DBH_mort-DBH)))
real, protected :: DBH_mort   = 0.025 ! characteristic DBH for mortality
real, protected :: A_mort     = 4.0   ! A coefficient in understory mortality rate correction, 1/year
real, protected :: B_mort     = 30.0  ! B coefficient in understory mortality rate correction, 1/m
real, protected :: nsc_starv_frac = 0.01 ! if NSC drops below bl_max multiplied by this value, cohort dies
real, protected :: cold_mort = 2.0 ! mortality rate due to coldest month below Tmin_mort threshold, 1/year
real, protected :: treeline_mort = 2.0 ! mortality rate above treeline, 1/year

real, protected :: DBH_merge_rel = 0.15  ! max relative DBH difference that permits merge of two cohorts
real, protected :: DBH_merge_abs = 0.003 ! max absolute DBH difference (m) that permits merge of two cohorts
real, protected :: NSC_merge_rel = 0.15  ! max relative NSC difference that allows merge of grass cohorts

real :: mycorrhizal_turnover_time = 0.1     ! Mean residence time of live mycorrhizal biomass (yr)
real :: myc_scav_C_efficiency     = 0.8     ! Efficiency of C allocation to scavenger mycorrhizae (remainder goes to CO2)
real :: myc_mine_C_efficiency     = 0.8     ! Efficiency of C allocation to miner mycorrhizae (remainder goes to CO2)
real :: c2n_N_fixer           = 10      ! C:N ratio of N-fixing microbe biomass
real :: N_fixer_turnover_time = 0.1     ! Mean residence time of live N fixer biomass (yr)
real :: N_fixer_C_efficiency  = 0.5     ! Efficiency of C allocation to N fixers (remainder goes to CO2)
logical :: N_limits_live_biomass = .FALSE.  ! Option to have N uptake limit max biomass.  Only relevant with CORPSE_N
real :: excess_stored_N_leakage_rate = 1.0 ! Leaking of excess cohort stored N back to soil (Fraction per year)
real :: min_N_stress = 0.05            ! Minimum value for N stress
real :: et_myc = 0.7                   ! Fraction of mycorrhizal turnover NOT mineralized to CO2 and NH4

logical :: calc_SLA_from_lifespan      ! In LM3, whether to calculate SLA from leaf lifespan or use namelist value

logical :: smooth_N_uptake_C_allocation = .FALSE.
logical :: N_fix_Tdep_Houlton = .FALSE.

namelist /vegn_data_nml/ &
  vegn_to_use,  input_cover_types, &
  mcv_min, mcv_lai, &
  use_bucket, use_mcm_masking, vegn_index_constant, &
  critical_root_density, &
  min_cosz, &
  soil_carbon_depth_scale, cold_month_threshold, &

  agf_bs, K1,K2, &
  tau_drip_l, tau_drip_s, GR_factor, tg_c3_thresh, tg_c4_thresh, T_cold_tropical,&
  fsc_pool_spending_time, ssc_pool_spending_time, harvest_spending_time, &
  T_transp_min, &
  phen_ev1, phen_ev2, &
  treeline_base_T, treeline_thresh_T, treeline_season_length, &
  scnd_biomass_bins, &

  ! PPA-related namelist values
  do_ppa, &
  cmc_eps, &
  DBH_mort, A_mort, B_mort, nsc_starv_frac, cold_mort, treeline_mort, &
  b0_growth, tau_seed, understory_lai_factor, min_lai, min_cohort_nindivs, &
  nat_mortality_splits_tiles, &
  DBH_merge_rel, DBH_merge_abs, NSC_merge_rel, &

  ! N-related namelist values
  mycorrhizal_turnover_time, &
  myc_scav_C_efficiency, myc_mine_C_efficiency, &
  N_fixer_turnover_time, N_fixer_C_efficiency, &
  c2n_N_fixer, N_limits_live_biomass, &
  excess_stored_N_leakage_rate, min_N_stress, calc_SLA_from_lifespan,&
  et_myc, smooth_N_uptake_C_allocation, N_fix_Tdep_Houlton

contains ! ###################################################################


! ============================================================================
subroutine read_vegn_data_namelist()
  ! ---- local vars
  integer :: unit         ! unit for namelist i/o
  integer :: io           ! i/o status for the namelist
  integer :: ierr         ! error code, returned by i/o routines
  character(fm_field_name_len) :: name, name0 ! name of the species
  character(fm_type_name_len)  :: ftype ! type of the field table entry
  type(fm_list_iter_type) :: iter ! iterator over the list of species
  integer :: i, k, n
  integer :: species_errors, total_errors
  logical, allocatable :: spdata_set(:)

  call log_version(version, module_name, &
  __FILE__)
#ifdef INTERNAL_FILE_NML
  read (input_nml_file, nml=vegn_data_nml, iostat=io)
  ierr = check_nml_error(io, 'vegn_data_nml')
#else
  if (file_exist('input.nml')) then
     unit = open_namelist_file()
     ierr = 1;
     do while (ierr /= 0)
        read (unit, nml=vegn_data_nml, iostat=io, end=10)
        ierr = check_nml_error (io, 'vegn_data_nml')
     enddo
10   continue
     call close_file (unit)
  endif
#endif

  unit=stdlog()
  write (unit, nml=vegn_data_nml)

  if(.not.fm_dump_list('/land_mod/species', recursive=.TRUE.)) &
     call error_mesg(module_name,'Cannot dump field list "/land_mod/species"',FATAL)

  ! count species
  nspecies = fm_get_length('/land_mod/species')
  call error_mesg(module_name,'total number of species is '//string(nspecies),NOTE)
  ! check that the number of species is consistent with the model settings
  if((.not.do_ppa).and.nspecies<N_LM3_SPECIES) call error_mesg ( module_name,&
     'Too few species to run in LM3 mode: must be at least '//string(N_LM3_SPECIES), FATAL)
  ! allocate species data
  allocate(spdata(0:nspecies-1), spdata_set(0:nspecies-1))
  spdata_set(:) = .FALSE.

  ! prefill the names of species in LM3 mode: to preserve species indices at they
  ! were in original LM3, since some of the subroutines rely on these indices
  if (.not.do_ppa) then
     do i = 0, N_LM3_SPECIES-1
        spdata(i)%name = lm3_species_name(i)
     enddo
  endif

  ! read data for default species and fill the species table with default parameters.
  ! NOTE that we ignore errors here; they will be caught during the main pass through
  ! the field table anyway.
  call fm_init_loop('/land_mod/species',iter)
  do while (fm_loop_over_list(iter, name, ftype, n))
     if (trim(name) == 'default' ) then
        i = ubound(spdata,1)
        ! write(*,*) 'found default =',i
        call read_species_data(name, spdata(i), species_errors)
        ! write(*,*) 'read default data =',i, species_errors
        do k = lbound(spdata,1),ubound(spdata,1)
           if (i==k) cycle ! skip the default species
           ! set all species data to defaults, except the name
           name0 = spdata(k)%name
           spdata(k) = spdata(i)
           spdata(k)%name = name0
        enddo
     endif
  enddo

  ! go through the species table and read the data
  call fm_init_loop('/land_mod/species',iter)
  total_errors = 0
  do while (fm_loop_over_list(iter, name, ftype, n))
     i = species_slot(name)
     if (i>=nspecies) then
        call error_mesg(data_error_header, &
          'could not place species "'//trim(name)//'" into species table, check if you are running in LM3 mode and misspelled species name', NOTE)
        total_errors = total_errors+1
     endif
     if (spdata_set(i)) then
        call error_mesg(data_error_header, &
          'data for species "'//trim(name)//'" are set more than once in the field table',NOTE)
        total_errors = total_errors+1
     endif

     call read_species_data(name,spdata(i), species_errors)
     if (species_errors==0) then
        call init_derived_species_data(spdata(i))
        spdata_set(i) = .TRUE.
     endif
     total_errors = total_errors + species_errors
  enddo
  if (total_errors>0) &
     call error_mesg(module_name, trim(string(total_errors))//' errors found in species parameters tables, look for "'//&
                                  data_error_header//'" in this output', FATAL)

  if (.not.all(spdata_set)) call error_mesg(module_name,'not all species data were set',FATAL)
  deallocate(spdata_set)

  ! reconcile values of fact_crit_phen and cnst_crit_phen
  spdata(:)%cnst_crit_phen = max(0.0,min(1.0,spdata(:)%cnst_crit_phen))
  spdata(:)%fact_crit_phen = max(0.0,spdata(:)%fact_crit_phen)
  where (spdata(:)%cnst_crit_phen/=0) spdata(:)%fact_crit_phen=0.0
  write(unit,*)'reconciled fact_crit_phen and cnst_crit_phen'

  ! do the same for fire
  spdata(:)%cnst_crit_fire = max(0.0,min(1.0,spdata(:)%cnst_crit_fire))
  spdata(:)%fact_crit_fire = max(0.0,spdata(:)%fact_crit_fire)
  where (spdata(:)%cnst_crit_fire/=0) spdata(:)%fact_crit_fire=0.0
  write(unit,*)'reconciled fact_crit_fire and cnst_crit_fire'

  call print_species_data(stdout(),.FALSE.)
  call print_species_data(stdlog(),.FALSE.)

  ! register selectors for land use type-specific diagnostics
  do i=1, N_LU_TYPES
     call register_tile_selector(landuse_name(i), long_name=landuse_longname(i),&
          tag = SEL_VEGN, idata1 = LU_SEL_TAG, idata2 = i )
  enddo

  ! register selectors for species-specific diagnostics
!   if (.not.do_ppa) then
!      do i=0,nspecies-1
!         call register_tile_selector(spdata(i)%name, long_name=spdata(i)%longname,&
!              tag = SEL_VEGN, idata1 = SP_SEL_TAG, idata2 = i )
!      enddo
!   endif

  ! register selector for natural grass
  call register_tile_selector('ntrlgrass', long_name='natural (non-human-maintained) grass',&
          tag = SEL_VEGN, idata1 = NG_SEL_TAG)

contains

  function species_slot(name) result(i); integer :: i
     character(*), intent(in) :: name

     do i = 0,nspecies-1
        if (trim(spdata(i)%name)==trim(name)) exit ! found slot for species
     enddo
     ! if not found, look for an empty slot in the table
     if (i>=nspecies) then
        do i = 0, nspecies
           if (trim(spdata(i)%name)=='') exit ! found an empty slot
        enddo
     endif
  end function species_slot

end subroutine read_vegn_data_namelist


! ============================================================================
! given a name of the species, and a species data structures, fills the
! the species data with information from the field table
subroutine read_species_data(name, sp, errors_found)
  character(*)        , intent(in)    :: name ! name of the species
  type(spec_data_type), intent(inout) :: sp   ! data to be filled in
  integer             , intent(out)   :: errors_found ! number of errors found in the species parameter list

  character(fm_field_name_len) :: str
  character(fm_type_name_len)  :: ftype ! type of the field table entry
  character(fm_path_name_len)  :: listname  ! name of the field manager list for the vegetation species
  character(fm_path_name_len)  :: current_list ! storage for current location in the fiels manager tree
  character(fm_field_name_len), allocatable :: known_names(:) ! names of all the parameters code attempts to read
  integer :: n_names ! number of parameters code attempts to read
  type(fm_list_iter_type) :: iter ! iterator over the list of species parameters
  integer :: j,n

  ! allocate the table of parameter names for error detection
  ! this code records name every species parameter it tries to read to the table
  ! of known names; then it goes through the list of all parameters in the
  ! input species table, and checks that every name there is in the table of
  ! known names. This triggers an error if one or more of specified parameters
  ! are not know, protecting from user typos.
  n_names = 0
  allocate(known_names(1))

  ! save current position in the field manager tree to restore it on exit
  current_list = fm_get_current_list()
  if (current_list .eq. ' ') &
      call error_mesg(module_name,'Could not get the current list',FATAL)

  listname = '/land_mod/species/'//trim(name)

  if (.not.fm_change_list(listname)) then
     call error_mesg(module_name,'Cannot change fm list to "'//trim(listname)//'"', FATAL)
  endif

  sp%name = name

  call add_known_name('long_name')
  sp%longname = fm_util_get_string('long_name', caller = module_name, default_value = '', scalar = .true.)

  call add_known_name('leaf_refl')
  call add_known_name('leaf_tran')
  do j = 1, NBANDS
     sp%leaf_refl(j) = fm_util_get_real('leaf_refl', &
              caller=module_name, default_value=sp%leaf_refl(j), index=j)
     sp%leaf_tran(j) = fm_util_get_real('leaf_tran', &
              caller=module_name, default_value=sp%leaf_tran(j), index=j)
  enddo

  call add_known_name('physiology_type')
  str = fm_util_get_string('physiology_type', caller = module_name, default_value = 'C3', scalar = .true.)
  select case (trim(lowercase(str)))
  case('c3')
     sp%pt = PT_C3
  case('c4')
     sp%pt = PT_C4
  case default
     call error_mesg(module_name,'Physiology type "'//trim(str)//'" is invalid, use "C3" or "C4"', FATAL)
  end select

  call add_known_name('phenology_type')
  str = fm_util_get_string('phenology_type', caller = module_name, default_value = 'deciduous', scalar = .true.)
  select case (trim(lowercase(str)))
  case('deciduous')
     sp%phent = PHEN_DECIDUOUS
  case('evergreen')
     sp%phent = PHEN_EVERGREEN
  case default
     call error_mesg(module_name,'Phenology type "'//trim(str)//'" is invalid, use "deciduous" or "evergreen"', FATAL)
  end select

  call add_known_name('allometry_type')
  str = fm_util_get_string('allometry_type', caller = module_name, default_value = 'HML', scalar = .true.)
  select case (trim(lowercase(str)))
  case('ew')
     sp%allomt = ALLOM_EW
  case('ew1')
     sp%allomt = ALLOM_EW1
  case('hml')
     sp%allomt = ALLOM_HML
  case default
     call error_mesg(module_name,'Allometry type "'//trim(str)//'" is invalid, use "EW", "EW1", or "HML"', FATAL)
  end select

  call add_known_name('lifeform')
  str = fm_util_get_string('lifeform', caller = module_name, default_value = 'tree', scalar = .true.)
  select case (trim(lowercase(str)))
  case('tree')
     sp%lifeform = FORM_WOODY
  case('grass')
     sp%lifeform = FORM_GRASS
  case default
     call error_mesg(module_name,'Vegetation lifeform "'//trim(str)//'" is invalid, use "tree" or "grass"', FATAL)
  end select

#define __GET_SPDATA_REAL__(v) sp%v = get_spdata_real(#v, sp%v)
#define __GET_SPDATA_LOGICAL__(v) sp%v = get_spdata_logical(#v, sp%v)
  __GET_SPDATA_REAL__(treefall_disturbance_rate)

  __GET_SPDATA_REAL__(c1)
  __GET_SPDATA_REAL__(c2)
  __GET_SPDATA_REAL__(c3)

  __GET_SPDATA_REAL__(alpha_leaf)
  __GET_SPDATA_REAL__(alpha_root)
  __GET_SPDATA_REAL__(alpha_wood)

  __GET_SPDATA_REAL__(branch_loss_height)

  __GET_SPDATA_REAL__(beta_root)
  __GET_SPDATA_REAL__(beta_vleaf)
  __GET_SPDATA_REAL__(beta_sapwood)

  __GET_SPDATA_REAL__(Vmax)
  __GET_SPDATA_REAL__(m_cond)
  __GET_SPDATA_REAL__(alpha_phot)
  __GET_SPDATA_REAL__(gamma_resp)
  __GET_SPDATA_REAL__(wet_leaf_dreg)
  __GET_SPDATA_REAL__(leaf_age_onset)
  __GET_SPDATA_REAL__(leaf_age_tau)
  __GET_SPDATA_REAL__(Vmax_understory_factor)
  __GET_SPDATA_REAL__(Resp_understory_factor)

  __GET_SPDATA_REAL__(dfr)
  __GET_SPDATA_REAL__(srl)
  __GET_SPDATA_REAL__(root_r)
  __GET_SPDATA_REAL__(root_perm)

  __GET_SPDATA_REAL__(cmc_lai)
  __GET_SPDATA_REAL__(cmc_pow)
  __GET_SPDATA_REAL__(csc_lai)
  __GET_SPDATA_REAL__(csc_pow)

  __GET_SPDATA_REAL__(internal_gap_frac)
  __GET_SPDATA_REAL__(leaf_emis)
  __GET_SPDATA_REAL__(ksi)
  __GET_SPDATA_REAL__(leaf_size)

  __GET_SPDATA_REAL__(fuel_intensity)

  __GET_SPDATA_REAL__(tc_crit)
  __GET_SPDATA_REAL__(gdd_crit)
  __GET_SPDATA_REAL__(gdd_base_T)
  __GET_SPDATA_REAL__(psi_stress_crit_phen)
  __GET_SPDATA_REAL__(cnst_crit_phen)
  __GET_SPDATA_REAL__(fact_crit_phen)
  __GET_SPDATA_REAL__(leaf_C_retrans_frac)
  __GET_SPDATA_REAL__(root_C_retrans_frac)
  __GET_SPDATA_REAL__(cnst_crit_fire)
  __GET_SPDATA_REAL__(fact_crit_fire)

  __GET_SPDATA_REAL__(smoke_fraction)
  __GET_SPDATA_REAL__(LMA)
  __GET_SPDATA_REAL__(LMA_understory_factor)

  __GET_SPDATA_REAL__(dat_height)
  __GET_SPDATA_REAL__(dat_lai)
  __GET_SPDATA_REAL__(dat_root_density)
  __GET_SPDATA_REAL__(dat_root_zeta)
  __GET_SPDATA_REAL__(dat_rs_min)
  __GET_SPDATA_REAL__(dat_snow_crit)

  !  for PPA, Weng, 7/25/2011
  __GET_SPDATA_REAL__(alphaHT)
  __GET_SPDATA_REAL__(thetaHT)
  __GET_SPDATA_REAL__(gammaHT)
  __GET_SPDATA_REAL__(alphaCA)
  __GET_SPDATA_REAL__(thetaCA)
  __GET_SPDATA_REAL__(alphaBM)
  __GET_SPDATA_LOGICAL__(limit_tussock_R)
  __GET_SPDATA_REAL__(tussock_Ra)
  __GET_SPDATA_REAL__(tussock_Rb)
  __GET_SPDATA_REAL__(layer_height_factor)
  __GET_SPDATA_REAL__(maturalage)
  __GET_SPDATA_REAL__(v_seed)
  __GET_SPDATA_LOGICAL__(reproduces_in_understory)
  __GET_SPDATA_REAL__(seedling_height)
  __GET_SPDATA_REAL__(seedling_nsc_frac)
  __GET_SPDATA_REAL__(frac_seed_dispersed)
  __GET_SPDATA_REAL__(frac_seed_transported)
  __GET_SPDATA_REAL__(prob_g)
  __GET_SPDATA_REAL__(prob_e)
  __GET_SPDATA_REAL__(mortrate_d_c)
  __GET_SPDATA_REAL__(mortrate_d_u)
  __GET_SPDATA_REAL__(Tmin_mort)
  __GET_SPDATA_LOGICAL__(mortality_kills_balive)
  __GET_SPDATA_LOGICAL__(mortality_kills_seeds)

  __GET_SPDATA_REAL__(rho_wood)
  __GET_SPDATA_REAL__(taperfactor)
  __GET_SPDATA_REAL__(LAImax)
  __GET_SPDATA_REAL__(f_NSC)
  __GET_SPDATA_REAL__(f_WF)
  __GET_SPDATA_LOGICAL__(use_light_saber)
  __GET_SPDATA_REAL__(light_saber_Hmin)
  __GET_SPDATA_REAL__(phiRL)
  __GET_SPDATA_REAL__(phiCSA)
  !  for PPA, IMC, 1/8/2017
  __GET_SPDATA_REAL__(growth_resp)
  __GET_SPDATA_REAL__(NSC2targetbl)
  __GET_SPDATA_REAL__(NSC2targetbl0)
  __GET_SPDATA_REAL__(NSC2targetbl_dbh)
  __GET_SPDATA_REAL__(T_dorm)
  ! for Kok effect, ppg, 17/11/07
  __GET_SPDATA_REAL__(inib_factor)
  __GET_SPDATA_REAL__(light_kok)
  !for Temperature response, ppg, 17/11/08
  __GET_SPDATA_REAL__(ToptP)
  __GET_SPDATA_REAL__(TminP)
  __GET_SPDATA_REAL__(TmaxP)
  __GET_SPDATA_REAL__(ToptR)
  __GET_SPDATA_REAL__(TminR)
  __GET_SPDATA_REAL__(TmaxR)

  __GET_SPDATA_REAL__(Ea_ko)
  __GET_SPDATA_REAL__(Hd_ko)
  __GET_SPDATA_REAL__(Ea_kc)
  __GET_SPDATA_REAL__(Hd_kc)
  __GET_SPDATA_REAL__(Ea_vm)
  __GET_SPDATA_REAL__(Hd_vm)
  __GET_SPDATA_REAL__(Ea_gam)
  __GET_SPDATA_REAL__(Hd_gam)
  __GET_SPDATA_REAL__(Ea_resp)
  __GET_SPDATA_REAL__(Hd_resp)

  ! hydraulics-related parameters
  __GET_SPDATA_REAL__(Kxam)
  __GET_SPDATA_REAL__(dx)
  __GET_SPDATA_REAL__(cx)
  __GET_SPDATA_REAL__(Klam)
  __GET_SPDATA_REAL__(dl)
  __GET_SPDATA_REAL__(cl)
  __GET_SPDATA_REAL__(psi_tlp)

  __GET_SPDATA_LOGICAL__(do_N_mining_strategy)
  __GET_SPDATA_LOGICAL__(do_N_scavenging_strategy)
  __GET_SPDATA_LOGICAL__(do_N_fixation_strategy)
  __GET_SPDATA_REAL__(branch_wood_frac)

  __GET_SPDATA_REAL__(N_fixation_rate)

  __GET_SPDATA_REAL__(root_exudate_frac)
  __GET_SPDATA_REAL__(root_exudate_N_frac)
  __GET_SPDATA_LOGICAL__(dynamic_root_exudation)
  __GET_SPDATA_REAL__(kM_myc_growth)
  __GET_SPDATA_REAL__(myc_growth_rate)
  __GET_SPDATA_REAL__(myc_N_to_plant_rate)
  __GET_SPDATA_REAL__(root_NH4_uptake_rate)
  __GET_SPDATA_REAL__(root_NO3_uptake_rate)
  __GET_SPDATA_REAL__(k_ammonium_root_uptake)
  __GET_SPDATA_REAL__(k_nitrate_root_uptake)
  ! nitrogen-related parameters
  __GET_SPDATA_REAL__(fsc_liv)
  __GET_SPDATA_REAL__(fsc_froot)
  __GET_SPDATA_REAL__(fsc_wood)
  __GET_SPDATA_REAL__(leaf_live_c2n)
  __GET_SPDATA_REAL__(froot_live_c2n)
  __GET_SPDATA_REAL__(wood_c2n)
  __GET_SPDATA_REAL__(sapwood_c2n)
  __GET_SPDATA_REAL__(c2n_mycorrhizae)
  __GET_SPDATA_REAL__(leaf_N_retrans_frac)
  __GET_SPDATA_REAL__(root_N_retrans_frac)
  __GET_SPDATA_REAL__(max_n_stress_for_seed_production)
  __GET_SPDATA_REAL__(N_stress_root_factor)
  __GET_SPDATA_REAL__(tau_nsc_exudate)
  __GET_SPDATA_REAL__(tau_smooth_marginal_gain)
  __GET_SPDATA_REAL__(tau_smooth_alloc)
  __GET_SPDATA_REAL__(alloc_allowed_over_limit)
  __GET_SPDATA_REAL__(tau_smooth_Nstress)
  ! SSR fire parameters
  __GET_SPDATA_REAL__(ROS_max)
  __GET_SPDATA_REAL__(fire_duration)
  __GET_SPDATA_REAL__(CC_leaf)
  __GET_SPDATA_REAL__(CC_stem)
  __GET_SPDATA_REAL__(CC_litter)
  __GET_SPDATA_REAL__(fireMort_leaf)
  __GET_SPDATA_REAL__(fireMort_stem)
  __GET_SPDATA_REAL__(fireMort_root)
  __GET_SPDATA_REAL__(EF_CO2)
  __GET_SPDATA_REAL__(EF_CO)
  __GET_SPDATA_REAL__(F_scorch_height)
#undef __GET_SPDATA_LOGICAL__
#undef __GET_SPDATA_REAL__

  ! check for typos in the namelist: detects parameters that are listed in the
  ! field table, but their name is not one of the parameters model tries to read.
  ! TODO: check types of the parameters
  call fm_init_loop(listname,iter)
  errors_found = 0
  do while (fm_loop_over_list(iter, str, ftype, n))
     ! look for the namelist parameter in the list of names
     do j = 1,n_names
        if (trim(known_names(j))==trim(str)) exit
     enddo
     ! if not found, give an error message
     if (j>n_names) then
        call error_mesg(data_error_header,'Input parameter "'//trim(str)//'" for "'//trim(name)//'" is not known',NOTE)
        errors_found = errors_found + 1
     endif
  enddo

  if (.not.fm_change_list(current_list)) then
     call error_mesg(module_name,'Cannot change fm list back to "'//trim(current_list)//'"', FATAL)
  endif

  deallocate(known_names)

contains

   subroutine add_known_name(name)
     character(*), intent(in) :: name
     integer :: i
     character(fm_field_name_len), allocatable :: tmp(:)

     do i = 1,n_names
        if (trim(name)==trim(known_names(i))) &
            call error_mesg(module_name,'"'//trim(name)//'" is duplicated in species parameters',FATAL)
     enddo
     if (n_names+1>size(known_names)) then
         allocate(tmp(2*size(known_names)))
         tmp(1:n_names) = known_names(1:n_names)
         call move_alloc(tmp,known_names)
     endif
     known_names(n_names+1) = name
     n_names = n_names+1
   end subroutine add_known_name

   real function get_spdata_real(name,dflt) result (v)
      character(*), intent(in) :: name ! name of the field
      real        , intent(in) :: dflt ! default value
      v = fm_util_get_real(name, default_value=dflt, scalar=.true.)
      call add_known_name(name)
   end function get_spdata_real

   logical function get_spdata_logical(name,dflt) result (v)
      character(*), intent(in) :: name ! name of the field
      logical     , intent(in) :: dflt ! default value
      v = fm_util_get_logical(name, default_value=dflt, scalar=.true.)
      call add_known_name(name)
   end function get_spdata_logical
end subroutine read_species_data

! ============================================================================
subroutine init_derived_species_data(sp)
   type(spec_data_type), intent(inout) :: sp

   integer :: j
   real :: specific_leaf_area  ! m2/kgC
   real :: leaf_life_span      ! months

   if (do_ppa .or. .not. calc_SLA_from_lifespan) then
      ! LMA (leaf mass per unit area) comes directly from namelist
   else
      ! calculate specific leaf area (cm2/g(biomass))
      ! Global Raich et al 94 PNAS pp 13730-13734
      leaf_life_span     = 12.0/sp%alpha_leaf ! in months
      specific_leaf_area = 10.0**(2.4 - 0.46*log10(leaf_life_span));
      ! convert to m2/kgC
      specific_leaf_area = C2B*specific_leaf_area*1000.0/10000.0
      sp%LMA = 1.0/specific_leaf_area
      ! rho_wood is not used anywhere?
      ! ! the relationship from Moorcroft, based on Reich
      ! ! units kg C/m^3, hence the factor of 0.001 to convert from g/cm^3
      ! sp%rho_wood = (0.5 + 0.2*(leaf_life_span-1))*0.001;
      ! if (sp%rho_wood > 500.) sp%rho_wood = 0.5*0.001;
   endif
   sp%SRA = 2*PI*sp%root_r*sp%srl

   sp%phi1=0.5-0.633*sp%ksi-0.33*sp%ksi**2;
   sp%phi2=0.877*(1.0-2.0*sp%phi1);
   if(sp%ksi /= 0) then
      sp%mu_bar = &
           (1-sp%phi1/sp%phi2*log(1+sp%phi2/sp%phi1))&
           / sp%phi2
   else
      ! in degenerate case of spherical leaf angular distribution the above
      ! formula for mu_bar gives an undefined value, so we handle it separately
      sp%mu_bar = 1.0
   endif
   select case (sp%allomt)
   case (ALLOM_EW,ALLOM_EW1)
      sp%thetaHT = 0.5
      sp%thetaCA = sp%thetaHT + 1
      sp%thetaBM = sp%thetaHT + 2
      ! calculate alphaBM parameter of allometry
      ! note that rho_wood was re-introduced for this calculation
      ! Isa changed for cross sectional area
      sp%alphaBM = sp%taperfactor * PI/4. * sp%alphaHT
      sp%alphaCSASW = sp%phiCSA*sp%LAImax*sp%alphaCA
      sp%thetaCSASW = sp%thetaCA
   case (ALLOM_HML)
      ! for HML allometry, parameters come directly from the input, so no
      ! calculations are needed here
   end select

  ! Convert units: from MPa to Pa
  sp%Kxam = sp%Kxam * 1e-6
  sp%Klam = sp%Klam * 1e-6
  sp%dx   = sp%dx   * 1e6
  sp%dl   = sp%dl   * 1e6
  sp%psi_tlp = sp%psi_tlp * 1e6  ! pressure in numerator
  ! TLP need not be less than mortality threshold
  ! TODO: make mortality threshold a namelist (or species) parameter?
  sp%psi_tlp = max(sp%psi_tlp, sp%dx*(-log(0.1))**(1./sp%cx))

  ! convert units : from deg C to deg K
  sp%ToptP = sp%ToptP + TFREEZE
  sp%TminP = sp%TminP + TFREEZE
  sp%TmaxP = sp%TmaxP + TFREEZE
  sp%ToptR = sp%ToptR + TFREEZE
  sp%TminR = sp%TminR + TFREEZE
  sp%TmaxR = sp%TmaxR + TFREEZE

  ! TODO: calculate seed C:N ratio

end subroutine init_derived_species_data

! ============================================================================
! prints a table of species parameters to specified unit
subroutine print_species_data(unit, print_default_species)
  integer, intent(in) :: unit ! unit number to print to
  logical, intent(in) :: print_default_species

  type(table_printer_type) :: table
  integer :: i
  integer :: N
  
  if (print_default_species) then
     N = size(spdata)-1
  else
     N = size(spdata)-2
  endif

  call init_with_headers(table, spdata(0:N)%name)
  call add_row(table, 'index', [(i,i=0,N)])

  call add_row(table, 'Treefall dist. rate', spdata(0:N)%treefall_disturbance_rate)
  call add_row(table, 'Mortality kills balive', spdata(0:N)%mortality_kills_balive)
  call add_row(table, 'Physiology Type', spdata(0:N)%pt)
  call add_row(table, 'Phenology Type',  spdata(0:N)%phent)
  call add_row(table, 'Life Form',       spdata(0:N)%lifeform)
  call add_row(table, 'C1',            spdata(0:N)%c1)
  call add_row(table, 'C2',            spdata(0:N)%c2)
  call add_row(table, 'C3',            spdata(0:N)%c3)

  call add_row(table, 'alpha_leaf',    spdata(0:N)%alpha_leaf)
  call add_row(table, 'alpha_root',    spdata(0:N)%alpha_root)
  call add_row(table, 'alpha_wood',    spdata(0:N)%alpha_wood)

  call add_row(table, 'branch_loss_height', spdata(0:N)%branch_loss_height)

  call add_row(table, 'beta_root',     spdata(0:N)%beta_root)
  call add_row(table, 'beta_vleaf',    spdata(0:N)%beta_vleaf)
  call add_row(table, 'beta_sapwood',  spdata(0:N)%beta_sapwood)

  call add_row(table, 'dfr',           spdata(0:N)%dfr)

  call add_row(table, 'srl',           spdata(0:N)%srl)
  call add_row(table, 'root_r',        spdata(0:N)%root_r)
  call add_row(table, 'root_perm',     spdata(0:N)%root_perm)

  call add_row(table, 'leaf_size',     spdata(0:N)%leaf_size)

  call add_row(table, 'alpha_phot',    spdata(0:N)%alpha_phot)
  call add_row(table, 'm_cond',        spdata(0:N)%m_cond)
  call add_row(table, 'Vmax',          spdata(0:N)%Vmax)
  call add_row(table, 'gamma_resp',    spdata(0:N)%gamma_resp)
  call add_row(table, 'wet_leaf_dreg', spdata(0:N)%wet_leaf_dreg)
  call add_row(table, 'leaf_age_onset',spdata(0:N)%leaf_age_onset)
  call add_row(table, 'leaf_age_tau',  spdata(0:N)%leaf_age_tau)
  call add_row(table, 'Vmax_understory_factor', spdata(0:N)%Vmax_understory_factor)
  call add_row(table, 'Resp_understory_factor', spdata(0:N)%Resp_understory_factor)

  ! PPA-related parameters
  call add_row(table, 'Allometry Type',  spdata(0:N)%allomt)
  call add_row(table, 'alphaHT', spdata(0:N)%alphaHT)
  call add_row(table, 'thetaHT', spdata(0:N)%thetaHT)
  call add_row(table, 'gammaHT', spdata(0:N)%gammaHT)
  call add_row(table, 'alphaCA', spdata(0:N)%alphaCA)
  call add_row(table, 'thetaCA', spdata(0:N)%thetaCA)
  call add_row(table, 'alphaBM', spdata(0:N)%alphaBM)
  call add_row(table, 'thetaBM', spdata(0:N)%thetaBM)
  call add_row(table, 'alphaCSASW', spdata(0:N)%alphaCSASW)
  call add_row(table, 'thetaCSASW', spdata(0:N)%thetaCSASW)
  call add_row(table, 'limit_tussock_R', spdata(0:N)%limit_tussock_R)
  call add_row(table, 'tussock_Ra', spdata(0:N)%tussock_Ra)
  call add_row(table, 'tussock_Rb', spdata(0:N)%tussock_Rb)
  call add_row(table, 'layer_height_factor', spdata(0:N)%layer_height_factor)
  call add_row(table, 'maturalage', spdata(0:N)%maturalage)
  call add_row(table, 'v_seed', spdata(0:N)%v_seed)
  call add_row(table, 'reproduces_in_understory', spdata(0:N)%reproduces_in_understory)
  call add_row(table, 'frac_seed_dispersed', spdata(0:N)%frac_seed_dispersed)
  call add_row(table, 'frac_seed_transported', spdata(0:N)%frac_seed_transported)
  call add_row(table, 'seedling_height', spdata(0:N)%seedling_height)
  call add_row(table, 'seedling_nsc_frac', spdata(0:N)%seedling_nsc_frac)
  call add_row(table, 'prob_g', spdata(0:N)%prob_g)
  call add_row(table, 'prob_e', spdata(0:N)%prob_e)
  call add_row(table, 'mortrate_d_c', spdata(0:N)%mortrate_d_c)
  call add_row(table, 'mortrate_d_u', spdata(0:N)%mortrate_d_u)
  call add_row(table, 'Tmin_mort', spdata(0:N)%Tmin_mort)
  call add_row(table, 'mortality_kills_seeds', spdata(0:N)%mortality_kills_seeds)
  call add_row(table, 'LMA', spdata(0:N)%LMA)
  call add_row(table, 'LMA_understory_factor', spdata(0:N)%LMA_understory_factor)
  call add_row(table, 'rho_wood', spdata(0:N)%rho_wood)
  call add_row(table, 'taperfactor', spdata(0:N)%taperfactor)
  call add_row(table, 'LAImax', spdata(0:N)%LAImax)
  call add_row(table, 'f_NSC', spdata(0:N)%f_NSC)
  call add_row(table, 'f_WF', spdata(0:N)%f_WF)
  call add_row(table, 'use_light_saber', spdata(0:N)%use_light_saber)
  call add_row(table, 'light_saber_Hmin', spdata(0:N)%light_saber_Hmin)
  call add_row(table, 'phiRL', spdata(0:N)%phiRL)
  call add_row(table, 'SRA', spdata(0:N)%SRA)
  call add_row(table, 'growth_resp', spdata(0:N)%growth_resp)
  call add_row(table, 'NSC2targetbl', spdata(0:N)%NSC2targetbl)
  call add_row(table, 'NSC2targetbl0', spdata(0:N)%NSC2targetbl0)
  call add_row(table, 'NSC2targetbl_dbh', spdata(0:N)%NSC2targetbl_dbh)
  call add_row(table, 'T_dorm', spdata(0:N)%T_dorm)

  call add_row(table, 'Klam', spdata(0:N)%Klam)
  call add_row(table, 'dl', spdata(0:N)%dl)
  call add_row(table, 'cl', spdata(0:N)%cl)
  call add_row(table, 'Kxam', spdata(0:N)%Kxam)
  call add_row(table, 'dx', spdata(0:N)%dx)
  call add_row(table, 'cx', spdata(0:N)%cx)
  call add_row(table, 'psi_tlp', spdata(0:N)%psi_tlp)

  !for kok effect, ppg, 17/11/07
  call add_row(table, 'inib_factor', spdata(0:N)%inib_factor)
  call add_row(table, 'light_kok', spdata(0:N)%light_kok)

  !for Temperature response, ppg, 17/11/08, reporting in degC, same as input
  call add_row(table, 'ToptP', spdata(0:N)%ToptP-TFREEZE)
  call add_row(table, 'TminP', spdata(0:N)%TminP-TFREEZE)
  call add_row(table, 'TmaxP', spdata(0:N)%TmaxP-TFREEZE)
  call add_row(table, 'ToptR', spdata(0:N)%ToptR-TFREEZE)
  call add_row(table, 'TminR', spdata(0:N)%TminR-TFREEZE)
  call add_row(table, 'TmaxR', spdata(0:N)%TmaxR-TFREEZE)

  call add_row(table, 'Ea_ko', spdata(0:N)%Ea_ko)
  call add_row(table, 'Hd_ko', spdata(0:N)%Hd_ko)
  call add_row(table, 'Ea_kc', spdata(0:N)%Ea_kc)
  call add_row(table, 'Hd_kc', spdata(0:N)%Hd_kc)
  call add_row(table, 'Ea_vm', spdata(0:N)%Ea_vm)
  call add_row(table, 'Hd_vm', spdata(0:N)%Hd_vm)
  call add_row(table, 'Ea_gam', spdata(0:N)%Ea_gam)
  call add_row(table, 'Hd_gam', spdata(0:N)%Hd_gam)
  call add_row(table, 'Ea_resp', spdata(0:N)%Ea_resp)
  call add_row(table, 'Hd_resp', spdata(0:N)%Hd_resp)

  call add_row(table, 'leaf_refl_vis', spdata(0:N)%leaf_refl(BAND_VIS))
  call add_row(table, 'leaf_refl_nir', spdata(0:N)%leaf_refl(BAND_NIR))
  call add_row(table, 'leaf_tran_vis', spdata(0:N)%leaf_tran(BAND_VIS))
  call add_row(table, 'leaf_tran_nir', spdata(0:N)%leaf_tran(BAND_NIR))
  call add_row(table, 'leaf_emis',     spdata(0:N)%leaf_emis)
  call add_row(table, 'ksi',           spdata(0:N)%ksi)
  call add_row(table, 'phi1',          spdata(0:N)%phi1)
  call add_row(table, 'phi2',          spdata(0:N)%phi2)
  call add_row(table, 'mu_bar',        spdata(0:N)%mu_bar)

  call add_row(table, 'cmc_lai',       spdata(0:N)%cmc_lai)
  call add_row(table, 'cmc_pow',       spdata(0:N)%cmc_pow)
  call add_row(table, 'csc_lai',       spdata(0:N)%csc_lai)
  call add_row(table, 'csc_pow',       spdata(0:N)%csc_pow)

  call add_row(table, 'internal_gap_frac', spdata(0:N)%internal_gap_frac)

  call add_row(table, 'fuel_intensity',spdata(0:N)%fuel_intensity)

  call add_row(table, 'tc_crit',       spdata(0:N)%tc_crit)
  call add_row(table, 'gdd_crit',      spdata(0:N)%gdd_crit)
  call add_row(table, 'gdd_base_T',    spdata(0:N)%gdd_base_T)

  call add_row(table, 'psi_stress_crit_phen', spdata(0:N)%psi_stress_crit_phen)
  call add_row(table, 'fact_crit_phen',spdata(0:N)%fact_crit_phen)
  call add_row(table, 'cnst_crit_phen',spdata(0:N)%cnst_crit_phen)
  call add_row(table, 'leaf_C_retrans_frac', spdata(0:N)%leaf_C_retrans_frac)
  call add_row(table, 'root_C_retrans_frac', spdata(0:N)%root_C_retrans_frac)

  call add_row(table, 'fact_crit_fire',spdata(0:N)%fact_crit_fire)
  call add_row(table, 'cnst_crit_fire',spdata(0:N)%cnst_crit_fire)
  call add_row(table, 'smoke_fraction',spdata(0:N)%smoke_fraction)

  call add_row(table, 'dynamic_root_exudation', spdata(0:N)%dynamic_root_exudation)
  call add_row(table, 'root_exudate_frac', spdata(0:N)%root_exudate_frac)
  call add_row(table, 'root_exudate_N_frac', spdata(0:N)%root_exudate_N_frac)
  call add_row(table, 'tau_nsc_exudate', spdata(0:N)%tau_nsc_exudate)

  call add_row(table, 'branch_wood_frac', spdata(0:N)%branch_wood_frac)

  call add_row(table, 'ROS_max',   spdata(0:N)%ROS_max)
  call add_row(table, 'fire_duration', spdata(0:N)%fire_duration)
  call add_row(table, 'CC_leaf',   spdata(0:N)%CC_leaf)
  call add_row(table, 'CC_stem',   spdata(0:N)%CC_stem)
  call add_row(table, 'CC_litter', spdata(0:N)%CC_litter)
  call add_row(table, 'fireMort_leaf', spdata(0:N)%fireMort_leaf)
  call add_row(table, 'fireMort_stem', spdata(0:N)%fireMort_stem)
  call add_row(table, 'fireMort_root', spdata(0:N)%fireMort_root)
  call add_row(table, 'EF_CO2', spdata(0:N)%EF_CO2)
  call add_row(table, 'EF_CO', spdata(0:N)%EF_CO)
  call add_row(table, 'F_scorch_height', spdata(0:N)%F_scorch_height)

  call add_row(table, 'do_N_mining_strategy', spdata(0:N)%do_N_mining_strategy)
  call add_row(table, 'do_N_scavenging_strategy', spdata(0:N)%do_N_scavenging_strategy)
  call add_row(table, 'do_N_fixation_strategy', spdata(0:N)%do_N_fixation_strategy)
  call add_row(table, 'N_fixation_rate', spdata(0:N)%N_fixation_rate)

  call add_row(table, 'kM_myc_growth', spdata(0:N)%kM_myc_growth)
  call add_row(table, 'myc_growth_rate', spdata(0:N)%myc_growth_rate)
  call add_row(table, 'myc_N_to_plant_rate', spdata(0:N)%myc_N_to_plant_rate)
  call add_row(table, 'root_NH4_uptake_rate', spdata(0:N)%root_NH4_uptake_rate)
  call add_row(table, 'root_NO3_uptake_rate', spdata(0:N)%root_NO3_uptake_rate)
  call add_row(table, 'k_ammonium_root_uptake', spdata(0:N)%k_ammonium_root_uptake)
  call add_row(table, 'k_nitrate_root_uptake', spdata(0:N)%k_nitrate_root_uptake)
  ! nitrogen-related parameters
  call add_row(table, 'fsc_liv',       spdata(0:N)%fsc_liv)
  call add_row(table, 'fsc_froot',     spdata(0:N)%fsc_froot)
  call add_row(table, 'fsc_wood',     spdata(0:N)%fsc_wood)
  call add_row(table, 'leaf_live_c2n', spdata(0:N)%leaf_live_c2n)
  call add_row(table, 'froot_live_c2n',spdata(0:N)%froot_live_c2n)
  call add_row(table, 'wood_c2n',      spdata(0:N)%wood_c2n)
  call add_row(table, 'sapwood_c2n',   spdata(0:N)%sapwood_c2n)
  call add_row(table, 'c2n_mycorrhizae', spdata(0:N)%c2n_mycorrhizae)
  call add_row(table, 'leaf_N_retrans_frac', spdata(0:N)%leaf_N_retrans_frac)
  call add_row(table, 'root_N_retrans_frac', spdata(0:N)%root_N_retrans_frac)
  call add_row(table, 'N_stress_root_factor', spdata(0:N)%N_stress_root_factor)
  call add_row(table, 'tau_smooth_marginal_gain', spdata(0:N)%tau_smooth_marginal_gain)
  call add_row(table, 'tau_smooth_alloc', spdata(0:N)%tau_smooth_alloc)
  call add_row(table, 'alloc_allowed_over_limit', spdata(0:N)%alloc_allowed_over_limit)
  call add_row(table, 'tau_smooth_Nstress', spdata(0:N)%tau_smooth_Nstress)
  call add_row(table, 'max_n_stress_for_seed_production', spdata(0:N)%max_n_stress_for_seed_production)

  call add_row(table, 'dat_height',       spdata(0:N)%dat_height)
  call add_row(table, 'dat_lai',          spdata(0:N)%dat_lai)
  call add_row(table, 'dat_root_density', spdata(0:N)%dat_root_density)
  call add_row(table, 'dat_root_zeta',    spdata(0:N)%dat_root_zeta)
  call add_row(table, 'dat_rs_min',       spdata(0:N)%dat_rs_min)
  call add_row(table, 'dat_snow_crit',    spdata(0:N)%dat_snow_crit)

  call print(table,unit,head_width=32)
end subroutine print_species_data

end module
