module vegn_data_mod

#ifdef INTERNAL_FILE_NML
use mpp_mod, only: input_nml_file
#else
use fms_mod, only: open_namelist_file
#endif

use constants_mod, only : PI
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
 N_LU_TYPES = 4, & ! number of different land use types
 LU_PAST    = 1, & ! pasture
 LU_CROP    = 2, & ! crops
 LU_NTRL    = 3, & ! natural vegetation
 LU_SCND    = 4    ! secondary vegetation
character(len=4), public, parameter  :: &
     landuse_name (N_LU_TYPES) = (/ 'past','crop','ntrl','scnd'/)
character(len=32), public, parameter :: &
     landuse_longname (N_LU_TYPES) = (/ 'pasture  ', 'crop     ', 'natural  ', 'secondary' /)

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

! ---- public types
public :: spec_data_type

! ---- public data
integer, public :: nspecies ! total number of species
public :: &
    vegn_to_use,  input_cover_types, &
    mcv_min, mcv_lai, &
    use_bucket, use_mcm_masking, vegn_index_constant, &
    critical_root_density, &
    ! vegetation data, imported from LM3V
    spdata, &
    min_cosz, &
    agf_bs, K1,K2, fsc_liv, fsc_wood, fsc_froot, &
    tau_drip_l, tau_drip_s, & ! canopy water and snow residence times, for drip calculations
    GR_factor, tg_c3_thresh, tg_c4_thresh, &
    fsc_pool_spending_time, ssc_pool_spending_time, harvest_spending_time, &
    l_fract, wood_fract_min, T_transp_min, soil_carbon_depth_scale, &
    cold_month_threshold, scnd_biomass_bins, &
    phen_ev1, phen_ev2, cmc_eps, &
    b0_growth, tau_seed, understory_lai_factor, &
    DBH_mort, A_mort, B_mort, mortrate_s

logical, public :: do_ppa = .FALSE.
logical, public :: do_alt_allometry = .FALSE.
logical, public :: nat_mortality_splits_tiles = .FALSE. ! if true, natural mortality
    ! creates disturbed tiles

! ---- public subroutine
public :: read_vegn_data_namelist
! ==== end of public interfaces ==============================================

! ==== constants =============================================================
character(len=*), parameter :: module_name = 'vegn_data_mod'
#include "../shared/version_variable.inc"
real, parameter :: TWOTHIRDS  = 2.0/3.0


! ==== types ================================================================
type spec_data_type
  character(fm_field_name_len) :: name     = '' ! name of the species
  character(fm_string_len)     :: longname = '' ! longname of the species
  real    :: treefall_disturbance_rate = 0.025;
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
  ! respiration rates of plant carbon pools
  real    :: beta_sapwood=0.0, beta_root=1.25, beta_vleaf=0.0

  real    :: dfr         = 5.8   ! fine root diameter ? or parameter relating diameter of fine roots to resistance
  ! the following two parameters are used in the Darcy-law calculations of water supply
  real    :: srl        = 24.4e3 ! specific root length, m/(kg C)
  real    :: root_r     = 2.9e-4 ! radius of the fine roots, m
  real    :: root_perm  = 5.0e-7 ! fine root membrane permeability per unit area, kg/(m3 s)

  real    :: LMA        = 0.036  ! leaf mass per unit area, kg C/m2
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
  real    :: alphaCA = 30.0,  thetaCA = 1.5 ! crown area alloametry parameters
  real    :: alphaBM = 0.559, thetaBM = 2.5 ! biomass allometry parameters
  real    :: alphaCSASW    = 2.25e-2, thetaCSASW = 1.5 !
  real    :: maturalage    = 1.0    ! the age that can reproduce
  real    :: fecundity     = 0.0    ! max C allocated to next generation per unit canopy area, kg C/m2
  real    :: v_seed        = 0.1    ! fraction of G_SF to G_F
  real    :: seedling_height   = 0.1 ! height of the seedlings, m
  real    :: seedling_nsc_frac = 3.0 ! initial seedling NSC, as fraction of bl_max (typically > 1)
  real    :: prob_g = 0.45, prob_e = 0.3 ! germination and establishment probabilities
  real    :: mortrate_d_c  = 0.05   ! daily mortality rate in canopy
  real    :: mortrate_d_u  = 0.2    ! daily mortality rate in understory
  real    :: rho_wood      = 250.0  ! woody density, kg C m-3 wood
  real    :: taperfactor   = 0.9
  real    :: LAImax        = 3.0    ! max. LAI
  real    :: phiRL         = 2.69   ! ratio of fine root to leaf area
  real    :: phiCSA        = 2.5e-4 ! ratio of sapwood CSA to target leaf area
  real    :: SRA           = 44.45982 ! specific fine root area, m2/kg C
  real    :: tauNSC        = 0.8    ! residence time of C in NSC (to define storage capacity)

  ! for hydraulics, wolf
  real    :: Kxam=0.0, Klam=0.0 ! Conductivity, max, per tissue area: units kg/m2 tissue/s/MPa
  real    :: dx=0.0, dl=0.0       ! Breakpoint of Weibull function, MPa
  real    :: cx=1.0, cl=1.0	  ! Exponent of Weibull function, unitless
  real    :: psi_tlp=0.0                  ! psi at turgor loss point

  real    :: root_exudate_frac = 0.0 ! fraction of NPP that ends up in root exudates
  real    :: branch_wood_frac = 0.1525 ! fraction of total wood biomass in branches,
                                       !  corresponds to 0.18 of trunk ( boule) biomass
                                       ! estimated by Isa from the observations
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
real :: fsc_liv        = 0.8
real :: fsc_wood       = 0.2
real :: fsc_froot      = 0.3
real :: tau_drip_l     = 21600.0 ! canopy water residence time, for drip calculations
real :: tau_drip_s     = 86400.0 ! canopy snow residence time, for drip calculations
real :: GR_factor = 0.33 ! growth respiration factor

real :: tg_c3_thresh = 1.5 ! threshold biomass between tree and grass for C3 plants
real :: tg_c4_thresh = 2.0 ! threshold biomass between tree and grass for C4 plants
real :: fsc_pool_spending_time = 1.0 ! time (yrs) during which intermediate pool of
                  ! fast soil carbon is entirely converted to the fast soil carbon
real :: ssc_pool_spending_time = 1.0 ! time (yrs) during which intermediate pool of
                  ! slow soil carbon is entirely converted to the slow soil carbon
real :: harvest_spending_time(N_HARV_POOLS) = &
     (/1.0, 1.0, 1.0, 1.0, 10.0, 100.0/)
     ! time (yrs) during which intermediate pool of harvested carbon is completely
     ! released to the atmosphere.
     ! NOTE: a year in the above *_spending_time definitions is exactly 365*86400 seconds
real :: l_fract      = 0.5 ! fraction of the leaves retained after leaf drop
real :: T_transp_min = 0.0 ! lowest temperature at which transpiration is enabled
                           ! 0 means no limit, lm3v value is 268.0
! Ensheng's growth parameters:
real :: b0_growth   = 0.02   ! min biomass for growth formula, kgC/indiv
real :: tau_seed    = 0.5708 ! characteristic time of nsc spending on seeds, year
real :: wood_fract_min = 0.33

! reduction of bl_max and br_max for the understory vegetation, unitless
real :: understory_lai_factor = 0.25

! boundaries of wood biomass bins for secondary veg. (kg C/m2); used to decide
! whether secondary vegetation tiles can be merged or not. MUST BE IN ASCENDING
! ORDER.
real  :: scnd_biomass_bins(10) &
     = (/ 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 1000.0 /)
real :: phen_ev1 = 0.5, phen_ev2 = 0.9 ! thresholds for evergreen/deciduous
      ! differentiation (see phenology_type in cohort.F90)
real :: cmc_eps = 0.01 ! value of w/w_max for transition to linear function;
                       ! the same value is used for liquid and snow

! Weng, 7/25/2011
! for understory mortality rate is calculated as:
! deathrate = mortrate_d_u * ( 1 + A * exp(B*(DBH_mort-DBH))/(1 + exp(B*(DBH_mort-DBH)))
real :: DBH_mort   = 0.025 ! characteristic DBH for mortality
real :: A_mort     = 4.0   ! A coefficient in understory mortality rate correction, 1/year
real :: B_mort     = 30.0  ! B coefficient in understory mortality rate correction, 1/m
real :: mortrate_s = 2.3   ! mortality rate of starving plants, 1/year, 2.3 = approx 0.9 plants die in a year


namelist /vegn_data_nml/ &
  vegn_to_use,  input_cover_types, &
  mcv_min, mcv_lai, &
  use_bucket, use_mcm_masking, vegn_index_constant, &
  critical_root_density, &
  min_cosz, &
  soil_carbon_depth_scale, cold_month_threshold, &

  agf_bs, K1,K2, fsc_liv, fsc_wood, fsc_froot, &
  tau_drip_l, tau_drip_s, GR_factor, tg_c3_thresh, tg_c4_thresh, &
  fsc_pool_spending_time, ssc_pool_spending_time, harvest_spending_time, &
  l_fract, wood_fract_min, T_transp_min, &
  phen_ev1, phen_ev2, &
  scnd_biomass_bins, &

  ! PPA-related namelist values
  do_ppa, &
  mortrate_s, cmc_eps, &
  DBH_mort, A_mort, B_mort, &
  b0_growth, tau_seed, understory_lai_factor, &
  do_alt_allometry, nat_mortality_splits_tiles

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

! TODO: possibly create a set of default species parameters, and set them through
!       the namelist
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
        i = species_slot(name)
        write(*,*) 'found default =',i
        call read_species_data(name, spdata(i), species_errors)
        write(*,*) 'read default data =',i, species_errors
        do k = lbound(spdata,1),ubound(spdata,1)
           if (i==k) cycle ! skip the deafult species
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
        call error_mesg(module_name, &
          'could not place species "'//trim(name)//'" into species table, check if you are running in LM3 mode and misspelled species name', NOTE)
        total_errors = total_errors+1
     endif
     if (spdata_set(i)) then
        call error_mesg(module_name, &
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
     call error_mesg(module_name, trim(string(total_errors))//' errors found in species parameters tables, see NOTES above', FATAL)

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

  call print_species_data(stdout())
  call print_species_data(stdlog())

  ! register selectors for land use type-specific diagnostics
  do i=1, N_LU_TYPES
     call register_tile_selector(landuse_name(i), long_name=landuse_longname(i),&
          tag = SEL_VEGN, idata1 = LU_SEL_TAG, idata2 = i )
  enddo

  ! register selectors for species-specific diagnostics
  do i=0,nspecies-1
     call register_tile_selector(spdata(i)%name, long_name=spdata(i)%longname,&
          tag = SEL_VEGN, idata1 = SP_SEL_TAG, idata2 = i )
  enddo

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

  call add_known_name('mortality_kills_balive')
  sp%mortality_kills_balive = fm_util_get_logical('mortality_kills_balive', &
        caller=module_name, default_value=sp%mortality_kills_balive, scalar=.true.)
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
  __GET_SPDATA_REAL__(treefall_disturbance_rate)

  __GET_SPDATA_REAL__(c1)
  __GET_SPDATA_REAL__(c2)
  __GET_SPDATA_REAL__(c3)

  __GET_SPDATA_REAL__(alpha_leaf)
  __GET_SPDATA_REAL__(alpha_root)
  __GET_SPDATA_REAL__(alpha_wood)

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
  __GET_SPDATA_REAL__(cnst_crit_fire)
  __GET_SPDATA_REAL__(fact_crit_fire)

  __GET_SPDATA_REAL__(smoke_fraction)
  __GET_SPDATA_REAL__(LMA)

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
  __GET_SPDATA_REAL__(maturalage)
  __GET_SPDATA_REAL__(v_seed)
  __GET_SPDATA_REAL__(seedling_height)
  __GET_SPDATA_REAL__(seedling_nsc_frac)
  __GET_SPDATA_REAL__(prob_g)
  __GET_SPDATA_REAL__(prob_e)
  __GET_SPDATA_REAL__(mortrate_d_c)
  __GET_SPDATA_REAL__(mortrate_d_u)
  __GET_SPDATA_REAL__(rho_wood)
  __GET_SPDATA_REAL__(taperfactor)
  __GET_SPDATA_REAL__(LAImax)
  __GET_SPDATA_REAL__(tauNSC)
  __GET_SPDATA_REAL__(phiRL)
  __GET_SPDATA_REAL__(phiCSA)
  ! hydraulics-related parameters
  __GET_SPDATA_REAL__(Kxam)
  __GET_SPDATA_REAL__(dx)
  __GET_SPDATA_REAL__(cx)
  __GET_SPDATA_REAL__(Klam)
  __GET_SPDATA_REAL__(dl)
  __GET_SPDATA_REAL__(cl)
  __GET_SPDATA_REAL__(psi_tlp)

  __GET_SPDATA_REAL__(root_exudate_frac)
  __GET_SPDATA_REAL__(branch_wood_frac)
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
        call error_mesg(module_name,'Input parameter "'//trim(str)//'" for "'//trim(name)//'" is not known',NOTE)
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

   function get_spdata_real(name,dflt) result (v) ; real :: v
      character(*), intent(in) :: name ! name of the field
      real        , intent(in) :: dflt ! default value
      v = fm_util_get_real(name, default_value=dflt, scalar=.true.)
      call add_known_name(name)
   end function get_spdata_real

end subroutine read_species_data

! ============================================================================
subroutine init_derived_species_data(sp)
   type(spec_data_type), intent(inout) :: sp

   integer :: j
   real :: specific_leaf_area  ! m2/kgC
   real :: leaf_life_span      ! months

   if (do_ppa) then
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
      !Isa changed for cross sectional area
      sp%alphaBM = sp%taperfactor * PI/4. * sp%alphaHT
      sp%alphaCSASW = sp%phiCSA*sp%LAImax*sp%alphaCA
   !sp%thetaCSASW = sp%thetaCA
   end select
   !sp%alphaCSASW = sp%phiCSA*sp%LAImax*sp%alphaCA
   !sp%thetaCSASW = sp%thetaCA

  ! Convert units: from MPa to Pa
  sp%Kxam = sp%Kxam * 1e-6
  sp%Klam = sp%Klam * 1e-6
  sp%dx   = sp%dx   * 1e6
  sp%dl   = sp%dl   * 1e6
  sp%psi_tlp = sp%psi_tlp * 1e6  ! pressure in numerator
  ! TLP need not be less than mortality threshold
  ! TODO: make mortality threshold a namelist (or species) parameter?
  sp%psi_tlp = max(sp%psi_tlp, sp%dx*(-log(0.1))**(1./sp%cx))

end subroutine init_derived_species_data

! ============================================================================
! prints a table of species parameters to specified unit
subroutine print_species_data(unit)
  integer, intent(in) :: unit ! unit number to print to

  type(table_printer_type) :: table
  integer :: i

  call init_with_headers(table, spdata(:)%name)
  call add_row(table, 'index', [(i,i=0,size(spdata)-1)])

  call add_row(table, 'Treefall dist. rate', spdata(:)%treefall_disturbance_rate)
  call add_row(table, 'Mortality kills balive', spdata(:)%mortality_kills_balive)
  call add_row(table, 'Physiology Type', spdata(:)%pt)
  call add_row(table, 'Phenology Type',  spdata(:)%phent)
  call add_row(table, 'Life Form',       spdata(:)%lifeform)
  call add_row(table, 'C1',            spdata(:)%c1)
  call add_row(table, 'C2',            spdata(:)%c2)
  call add_row(table, 'C3',            spdata(:)%c3)

  call add_row(table, 'alpha_leaf',    spdata(:)%alpha_leaf)
  call add_row(table, 'alpha_root',    spdata(:)%alpha_root)
  call add_row(table, 'alpha_wood',    spdata(:)%alpha_wood)

  call add_row(table, 'beta_root',     spdata(:)%beta_root)
  call add_row(table, 'beta_vleaf',    spdata(:)%beta_vleaf)
  call add_row(table, 'beta_sapwood',  spdata(:)%beta_sapwood)

  call add_row(table, 'dfr',           spdata(:)%dfr)

  call add_row(table, 'srl',           spdata(:)%srl)
  call add_row(table, 'root_r',        spdata(:)%root_r)
  call add_row(table, 'root_perm',     spdata(:)%root_perm)

  call add_row(table, 'leaf_size',     spdata(:)%leaf_size)

  call add_row(table, 'alpha_phot',    spdata(:)%alpha_phot)
  call add_row(table, 'm_cond',        spdata(:)%m_cond)
  call add_row(table, 'Vmax',          spdata(:)%Vmax)
  call add_row(table, 'gamma_resp',    spdata(:)%gamma_resp)
  call add_row(table, 'wet_leaf_dreg', spdata(:)%wet_leaf_dreg)
  call add_row(table, 'leaf_age_onset',spdata(:)%leaf_age_onset)
  call add_row(table, 'leaf_age_tau',  spdata(:)%leaf_age_tau)
  call add_row(table, 'Vmax_understory_factor', spdata(:)%Vmax_understory_factor)
  call add_row(table, 'Resp_understory_factor', spdata(:)%Resp_understory_factor)

  ! PPA-related parameters
  call add_row(table, 'Allometry Type',  spdata(:)%allomt)
  call add_row(table, 'alphaHT', spdata(:)%alphaHT)
  call add_row(table, 'thetaHT', spdata(:)%thetaHT)
  call add_row(table, 'gammaHT', spdata(:)%gammaHT)
  call add_row(table, 'alphaCA', spdata(:)%alphaCA)
  call add_row(table, 'thetaCA', spdata(:)%thetaCA)
  call add_row(table, 'alphaBM', spdata(:)%alphaBM)
  call add_row(table, 'thetaBM', spdata(:)%thetaBM)
  call add_row(table, 'alphaCSASW', spdata(:)%alphaCSASW)
  call add_row(table, 'thetaCSASW', spdata(:)%thetaCSASW)
  call add_row(table, 'maturalage', spdata(:)%maturalage)
  call add_row(table, 'v_seed', spdata(:)%v_seed)
  call add_row(table, 'seedling_height', spdata(:)%seedling_height)
  call add_row(table, 'seedling_nsc_frac', spdata(:)%seedling_nsc_frac)
  call add_row(table, 'prob_g', spdata(:)%prob_g)
  call add_row(table, 'prob_e', spdata(:)%prob_e)
  call add_row(table, 'mortrate_d_c', spdata(:)%mortrate_d_c)
  call add_row(table, 'mortrate_d_u', spdata(:)%mortrate_d_u)
  call add_row(table, 'LMA', spdata(:)%LMA)
  call add_row(table, 'rho_wood', spdata(:)%rho_wood)
  call add_row(table, 'taperfactor', spdata(:)%taperfactor)
  call add_row(table, 'LAImax', spdata(:)%LAImax)
  call add_row(table, 'tauNSC', spdata(:)%tauNSC)
  call add_row(table, 'phiRL', spdata(:)%phiRL)
  call add_row(table, 'SRA', spdata(:)%SRA)

  call add_row(table, 'Klam', spdata(:)%Klam)
  call add_row(table, 'dl', spdata(:)%dl)
  call add_row(table, 'cl', spdata(:)%cl)
  call add_row(table, 'Kxam', spdata(:)%Kxam)
  call add_row(table, 'dx', spdata(:)%dx)
  call add_row(table, 'cx', spdata(:)%cx)
  call add_row(table, 'psi_tlp', spdata(:)%psi_tlp)

  call add_row(table, 'leaf_refl_vis', spdata(:)%leaf_refl(BAND_VIS))
  call add_row(table, 'leaf_refl_nir', spdata(:)%leaf_refl(BAND_NIR))
  call add_row(table, 'leaf_tran_vis', spdata(:)%leaf_tran(BAND_VIS))
  call add_row(table, 'leaf_tran_nir', spdata(:)%leaf_tran(BAND_NIR))
  call add_row(table, 'leaf_emis',     spdata(:)%leaf_emis)
  call add_row(table, 'ksi',           spdata(:)%ksi)
  call add_row(table, 'phi1',          spdata(:)%phi1)
  call add_row(table, 'phi2',          spdata(:)%phi2)
  call add_row(table, 'mu_bar',        spdata(:)%mu_bar)

  call add_row(table, 'cmc_lai',       spdata(:)%cmc_lai)
  call add_row(table, 'cmc_pow',       spdata(:)%cmc_pow)
  call add_row(table, 'csc_lai',       spdata(:)%csc_lai)
  call add_row(table, 'csc_pow',       spdata(:)%csc_pow)

  call add_row(table, 'internal_gap_frac', spdata(:)%internal_gap_frac)

  call add_row(table, 'fuel_intensity',spdata(:)%fuel_intensity)

  call add_row(table, 'tc_crit',       spdata(:)%tc_crit)
  call add_row(table, 'gdd_crit',      spdata(:)%gdd_crit)
  call add_row(table, 'gdd_base_T',    spdata(:)%gdd_base_T)

  call add_row(table, 'psi_stress_crit_phen', spdata(:)%psi_stress_crit_phen)
  call add_row(table, 'fact_crit_phen',spdata(:)%fact_crit_phen)
  call add_row(table, 'cnst_crit_phen',spdata(:)%cnst_crit_phen)
  call add_row(table, 'fact_crit_fire',spdata(:)%fact_crit_fire)
  call add_row(table, 'cnst_crit_fire',spdata(:)%cnst_crit_fire)

  call add_row(table, 'smoke_fraction',spdata(:)%smoke_fraction)

  call add_row(table, 'root_exudate_frac', spdata(:)%root_exudate_frac)
  call add_row(table, 'branch_wood_frac', spdata(:)%branch_wood_frac)

  call add_row(table, 'dat_height',       spdata(:)%dat_height)
  call add_row(table, 'dat_lai',          spdata(:)%dat_lai)
  call add_row(table, 'dat_root_density', spdata(:)%dat_root_density)
  call add_row(table, 'dat_root_zeta',    spdata(:)%dat_root_zeta)
  call add_row(table, 'dat_rs_min',       spdata(:)%dat_rs_min)
  call add_row(table, 'dat_snow_crit',    spdata(:)%dat_snow_crit)

  call print(table,unit)
end subroutine print_species_data

end module
