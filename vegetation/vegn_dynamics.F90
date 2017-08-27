! ============================================================================
! updates carbon pools and rates on the fast time scale
! ============================================================================
module vegn_dynamics_mod

#include "../shared/debug.inc"

use fms_mod, only: error_mesg, FATAL,NOTE
use time_manager_mod, only: time_type

use land_constants_mod, only : seconds_per_year, mol_C
use land_tile_diag_mod, only : register_tiled_diag_field, send_tile_data, &
     set_default_diag_filter, diag_buff_type, cmor_name

use vegn_data_mod, only : spdata, &
     CMPT_VLEAF, CMPT_SAPWOOD, CMPT_ROOT, CMPT_WOOD, CMPT_LEAF, LEAF_ON, LEAF_OFF, &
     fsc_liv, fsc_wood, soil_carbon_depth_scale, C2B, agf_bs, &
     l_fract, &
      root_exudate_N_frac,& !x2z - ens: lets get rid of c2n?
     ! BNS: C2N ratios should be temporary fix, which we can get rid of once N is integrated into vegetation code
     dynamic_root_exudation, c2n_mycorrhizae, mycorrhizal_turnover_time, myc_scav_C_efficiency,myc_mine_C_efficiency,&
     N_fixer_turnover_time, N_fixer_C_efficiency, N_fixation_rate, c2n_N_fixer, excess_stored_N_leakage_rate, N_limits_live_biomass, &
     myc_growth_rate, myc_N_to_plant_rate, et_myc, kM_myc_growth, smooth_N_uptake_C_allocation

use vegn_tile_mod, only: vegn_tile_type,vegn_tile_carbon, vegn_tile_nitrogen
use soil_tile_mod, only: soil_tile_type,soil_tile_carbon, soil_tile_nitrogen
use vegn_cohort_mod, only : vegn_cohort_type, &
     update_biomass_pools, update_bio_living_fraction, update_species
use soil_carbon_mod, only: soil_carbon_option, N_C_TYPES,&
    SOILC_CENTURY, SOILC_CENTURY_BY_LAYER, SOILC_CORPSE, SOILC_CORPSE_N, &
    add_litter, poolTotals, debug_pool,soil_NO3_deposition,soil_NH4_deposition,soil_org_N_deposition,deadmic_slow_frac

use soil_mod, only: add_root_litter, add_root_exudates, Dsdt, myc_scavenger_N_uptake, &
    myc_miner_N_uptake,root_N_uptake

use nitrogen_sources_mod, only: do_nitrogen_deposition

use land_debug_mod, only : is_watch_point
use land_data_mod, only : log_version

implicit none
private

! ==== public interfaces =====================================================
public :: vegn_dynamics_init

public :: vegn_carbon_int   ! fast time-scale integrator of carbon balance
public :: vegn_growth       ! slow time-scale redistributor of accumulated carbon
public :: vegn_daily_npp    ! updates values of daily-average npp
public :: vegn_phenology    !
public :: vegn_biogeography !
! ==== end of public interfaces ==============================================

! ==== module constants ======================================================
character(len=*), parameter :: module_name = 'vegn'
#include "../shared/version_variable.inc"

real, parameter :: GROWTH_RESP=0.333  ! fraction of npp lost as growth respiration

! ==== module data ===========================================================
real    :: dt_fast_yr ! fast (physical) time step, yr (year is defined as 365 days)

! diagnostic field IDs
integer :: id_npp, id_nep, id_gpp, id_resp, id_resl, id_resr, id_resg, &
    id_soilt, id_theta, id_litter, &
    id_mycorrhizal_scav_allocation, id_mycorrhizal_scav_immobilization,&
    id_mycorrhizal_mine_allocation, id_mycorrhizal_mine_immobilization, id_N_fixer_allocation, id_total_plant_N_uptake, &
    id_N_fix_marginal_gain, id_rhiz_exud_marginal_gain, id_myc_scav_marginal_gain, id_myc_mine_marginal_gain, id_rhiz_exudation, id_nitrogen_stress, &
    id_myc_scavenger_N_uptake,id_myc_miner_N_uptake,id_symbiotic_N_fixation,id_active_root_N_uptake,&
    id_scav_plant_N_uptake, id_mine_plant_N_uptake, id_fix_plant_N_uptake,&
    id_mycorrhizal_scav_C_res, id_mycorrhizal_scav_N_res, id_mycorrhizal_mine_C_res, id_mycorrhizal_mine_N_res, &
    id_Nfix_C_res, id_Nfix_N_res

! CMOR diagnostic field IDs
integer :: id_gpp_cmor, id_npp_cmor, id_ra, id_rgrowth

contains

! ============================================================================
subroutine vegn_dynamics_init(id_lon, id_lat, time, delta_time)
  integer        , intent(in) :: id_lon ! ID of land longitude (X) axis
  integer        , intent(in) :: id_lat ! ID of land latitude (Y) axis
  type(time_type), intent(in) :: time       ! initial time for diagnostic fields
  real           , intent(in) :: delta_time ! fast time step, s

  call log_version(version, 'vegn_dynamics_mod', &
  __FILE__)

  ! set up global variables
  dt_fast_yr = delta_time/seconds_per_year

  ! set the default sub-sampling filter for the fields below
  call set_default_diag_filter('soil')

  ! register diagnostic fields
  id_gpp = register_tiled_diag_field ( module_name, 'gpp',  &
       (/id_lon,id_lat/), time, 'gross primary productivity', 'kg C/(m2 year)', &
       missing_value=-100.0 )
  id_npp = register_tiled_diag_field ( module_name, 'npp',  &
       (/id_lon,id_lat/), time, 'net primary productivity', 'kg C/(m2 year)', &
       missing_value=-100.0 )
  id_nep = register_tiled_diag_field ( module_name, 'nep',  &
       (/id_lon,id_lat/), time, 'net ecosystem productivity', 'kg C/(m2 year)', &
       missing_value=-100.0 )
  id_litter = register_tiled_diag_field (module_name, 'litter', (/id_lon,id_lat/), &
       time, 'litter productivity', 'kg C/(m2 year)', missing_value=-100.0)
  id_resp = register_tiled_diag_field ( module_name, 'resp', (/id_lon,id_lat/), &
       time, 'respiration', 'kg C/(m2 year)', missing_value=-100.0 )
  id_resl = register_tiled_diag_field ( module_name, 'resl', (/id_lon,id_lat/), &
       time, 'leaf respiration', 'kg C/(m2 year)', missing_value=-100.0 )
  id_resr = register_tiled_diag_field ( module_name, 'resr', (/id_lon,id_lat/), &
       time, 'root respiration', 'kg C/(m2 year)', missing_value=-100.0 )
  id_resg = register_tiled_diag_field ( module_name, 'resg', (/id_lon,id_lat/), &
       time, 'growth respiration', 'kg C/(m2 year)', missing_value=-100.0 )
  id_soilt = register_tiled_diag_field ( module_name, 'tsoil_av',  &
       (/id_lon,id_lat/), time, 'average soil temperature for carbon decomposition', 'degK', &
       missing_value=-100.0 )
  id_theta = register_tiled_diag_field ( module_name, 'theta',  &
       (/id_lon,id_lat/), time, 'average soil wetness for carbon decomposition', 'm3/m3', &
       missing_value=-100.0 )
  id_mycorrhizal_scav_allocation = register_tiled_diag_field ( module_name, 'mycorrhizal_scav_allocation',  &
       (/id_lon,id_lat/), time, 'C allocation to scavenger mycorrhizae', 'kg C/(m2 year)', &
       missing_value=-100.0 )
   id_mycorrhizal_scav_immobilization = register_tiled_diag_field ( module_name, 'mycorrhizal_scav_immobilization',  &
        (/id_lon,id_lat/), time, 'N immobilization by scavenger mycorrhizae', 'kg N/(m2 year)', &
        missing_value=-100.0 )
        id_mycorrhizal_mine_allocation = register_tiled_diag_field ( module_name, 'mycorrhizal_mine_allocation',  &
             (/id_lon,id_lat/), time, 'C allocation to miner mycorrhizae', 'kg C/(m2 year)', &
             missing_value=-100.0 )
         id_mycorrhizal_mine_immobilization = register_tiled_diag_field ( module_name, 'mycorrhizal_mine_immobilization',  &
              (/id_lon,id_lat/), time, 'N immobilization by miner mycorrhizae', 'kg N/(m2 year)', &
              missing_value=-100.0 )
    id_N_fixer_allocation = register_tiled_diag_field ( module_name, 'N_fixer_allocation',  &
         (/id_lon,id_lat/), time, 'C allocation to N fixers', 'kg C/(m2 year)', &
         missing_value=-100.0 )
     id_N_fix_marginal_gain = register_tiled_diag_field ( module_name, 'N_fix_marginal_gain',  &
          (/id_lon,id_lat/), time, 'Extra N fixation per unit C allocation', 'kg N/(m2 year)/kgC', &
          missing_value=-100.0 )
    id_rhiz_exud_marginal_gain = register_tiled_diag_field ( module_name, 'rhiz_exud_marginal_gain',  &
         (/id_lon,id_lat/), time, 'Extra N acquisition per unit rhiz C exudation', 'kg N/(m2 year)/kgC', &
         missing_value=-100.0 )
    id_myc_scav_marginal_gain = register_tiled_diag_field ( module_name, 'myc_scav_marginal_gain',  &
         (/id_lon,id_lat/), time, 'Extra N acquisition per unit C allocation to scavenger mycorrhizae', 'kg N/(m2 year)/kgC', &
         missing_value=-100.0 )
         id_myc_mine_marginal_gain = register_tiled_diag_field ( module_name, 'myc_mine_marginal_gain',  &
              (/id_lon,id_lat/), time, 'Extra N acquisition per unit C allocation to miner mycorrhizae', 'kg N/(m2 year)/kg C', &
              missing_value=-100.0 )
     id_rhiz_exudation = register_tiled_diag_field ( module_name, 'rhiz_exudation',  &
          (/id_lon,id_lat/), time, 'C allocation to rhizosphere exudation', 'kg C/(m2 year)', &
          missing_value=-100.0 )
      id_nitrogen_stress = register_tiled_diag_field ( module_name, 'nitrogen_stress',  &
           (/id_lon,id_lat/), time, 'Nitrogen stress index', 'Dimensionless', &
           missing_value=-100.0 )
     id_total_plant_N_uptake = register_tiled_diag_field ( module_name, 'plant_N_uptake',  &
            (/id_lon,id_lat/), time, 'Plant N uptake rate', 'kg N/(m2 year)', &
            missing_value=-100.0 )
    id_myc_scavenger_N_uptake = register_tiled_diag_field ( module_name, 'myc_scavenger_N_uptake',  &
         (/id_lon,id_lat/), time, 'N uptake by scavenger mycorrhizae', 'kg N/m2/year', missing_value=-1.0 )
     id_myc_miner_N_uptake = register_tiled_diag_field ( module_name, 'myc_miner_N_uptake',  &
          (/id_lon,id_lat/), time, 'N uptake by miner mycorrhizae', 'kg N/m2/year', missing_value=-1.0 )
     id_symbiotic_N_fixation = register_tiled_diag_field ( module_name, 'symbiotic_N_fixation',  &
          (/id_lon,id_lat/), time, 'Symbiotic N fixation', 'kg N/m2/year', missing_value=-1.0 )
      id_active_root_N_uptake = register_tiled_diag_field ( module_name, 'active_root_N_uptake',  &
           (/id_lon,id_lat/), time, 'N uptake by root active transport', 'kg N/m2/year', missing_value=-1.0 )
      id_scav_plant_N_uptake = register_tiled_diag_field ( module_name, 'plant_scavenger_N_uptake',  &
           (/id_lon,id_lat/), time, 'Plant N uptake from scavenger mycorrhizae', 'kg N/m2/year', missing_value=-1.0 )
     id_mine_plant_N_uptake = register_tiled_diag_field ( module_name, 'plant_miner_N_uptake',  &
          (/id_lon,id_lat/), time, 'Plant N uptake from miner mycorrhizae', 'kg N/m2/year', missing_value=-1.0 )
     id_fix_plant_N_uptake = register_tiled_diag_field ( module_name, 'plant_fixer_N_uptake',  &
          (/id_lon,id_lat/), time, 'Plant N uptake from N fixers', 'kg N/m2/year', missing_value=-1.0 )
    id_mycorrhizal_scav_C_res = register_tiled_diag_field ( module_name, 'myc_scavenger_C_res',  &
         (/id_lon,id_lat/), time, 'Scavenger mycorrhizae C reservoir', 'kg C/m2', missing_value=-1.0 )
    id_mycorrhizal_scav_N_res = register_tiled_diag_field ( module_name, 'myc_scavenger_N_res',  &
         (/id_lon,id_lat/), time, 'Scavenger mycorrhizae N reservoir', 'kg N/m2', missing_value=-1.0 )
    id_mycorrhizal_mine_C_res = register_tiled_diag_field ( module_name, 'myc_miner_C_res',  &
         (/id_lon,id_lat/), time, 'Miner mycorrhizae C reservoir', 'kg C/m2', missing_value=-1.0 )
    id_mycorrhizal_mine_N_res = register_tiled_diag_field ( module_name, 'myc_miner_N_res',  &
         (/id_lon,id_lat/), time, 'Miner mycorrhizae N reservoir', 'kg N/m2', missing_value=-1.0 )
    id_Nfix_C_res = register_tiled_diag_field ( module_name, 'N_fixer_C_res',  &
         (/id_lon,id_lat/), time, 'N fixer C reservoir', 'kg C/m2', missing_value=-1.0 )
    id_Nfix_N_res = register_tiled_diag_field ( module_name, 'N_fixer_N_res',  &
         (/id_lon,id_lat/), time, 'N fixer N reservoir', 'kg N/m2', missing_value=-1.0 )

  ! set the default sub-sampling filter for CMOR variables
  call set_default_diag_filter('land')
  id_gpp_cmor = register_tiled_diag_field ( cmor_name, 'gpp', (/id_lon,id_lat/), &
       time, 'Gross Primary Production', 'kg C m-2 s-1', missing_value=-1.0, &
       standard_name='gross_primary_production', fill_missing=.TRUE.)
  id_npp_cmor = register_tiled_diag_field ( cmor_name, 'npp', (/id_lon,id_lat/), &
       time, 'Net Primary Production', 'kg C m-2 s-1', missing_value=-1.0, &
       standard_name='net_primary_production', fill_missing=.TRUE.)
  id_ra = register_tiled_diag_field ( cmor_name, 'ra', (/id_lon,id_lat/), &
       time, 'Autotrophic (Plant) Respiration', 'kg C m-2 s-1', missing_value=-1.0, &
       standard_name='autotrophic_plant_respiration', fill_missing=.TRUE.)
  id_rgrowth = register_tiled_diag_field ( cmor_name, 'rGrowth', (/id_lon,id_lat/), &
       time, 'Growth Autotrophic Respiration', 'kg C m-2 s-1', missing_value=-1.0, &
       standard_name='growth_autotrophic_respiration', fill_missing=.TRUE.)

end subroutine vegn_dynamics_init


! ============================================================================
subroutine vegn_carbon_int(vegn, soil, soilt, theta, ndep_nit, ndep_amm, ndep_org, diag)
  type(vegn_tile_type), intent(inout) :: vegn
  type(soil_tile_type), intent(inout) :: soil
  real, intent(in) :: soilt ! average temperature of soil for soil carbon decomposition, deg K
  real, intent(in) :: theta ! average soil wetness, unitless
  real, intent(in) :: ndep_nit, ndep_amm, ndep_org ! total nitrate, ammonium,
      ! and organic nitrogen inputs (deposition plus fertilization), kg N/(m2 yr)
  type(diag_buff_type), intent(inout) :: diag

  ! TODO: possibly move soil-related calculations from calling procedure here,
  !       now that we have soil passed as an argument

  type(vegn_cohort_type), pointer :: cc
  real :: resp, resl, resr, resg ! respiration terms accumulated for all cohorts
  real :: cgain, closs ! carbon gain and loss accumulated for entire tile
  real :: md_alive, md_leaf, md_vleaf, md_wood, md_froot, md_sapwood ! component of maintenance demand
  real :: gpp ! gross primary productivity per tile
  real :: root_exudate_C, total_root_exudate_C,root_exudate_N, total_root_exudate_N, root_exudate_frac
  real :: myc_scav_exudate_frac, mycorrhizal_scav_N_immob, myc_mine_exudate_frac, mycorrhizal_mine_N_immob
  real :: myc_scav_N_uptake,total_scavenger_myc_C_allocated,scavenger_myc_C_allocated,total_scav_myc_immob
  real :: scavenger_myc_N_allocated,N_fixer_N_allocated,miner_myc_N_allocated
  real :: myc_mine_N_uptake,myc_mine_C_uptake,total_miner_myc_C_allocated,miner_myc_C_allocated,total_mine_myc_immob,total_myc_mine_C_uptake
  real :: N_fixation, total_N_fixation, total_N_fixer_C_allocated, N_fixer_C_allocated, N_fixer_exudate_frac
  real :: myc_scav_marginal_gain,myc_mine_marginal_gain, N_fix_marginal_gain,rhiz_exud_frac,rhiz_exud_marginal_gain
  real :: weight_marginal_gain,weight_alloc
  real :: root_active_N_uptake, mining_CO2prod,total_mining_CO2prod, excess_mining_C
  real :: N_fixation_2, myc_N_uptake, myc_N_uptake_2, myc_C_uptake, myc_C_uptake_2, dummy1
  real :: total_plant_N_uptake, scavenger_myc_growth, miner_myc_growth, N_fixer_growth
  real :: scav_N_to_plant, mine_N_to_plant, fix_N_to_plant, total_myc_resp, total_myc_Nmin, N_leakage, maint_resp
  real :: d_scav_C_reservoir,d_scav_N_reservoir,d_mine_C_reservoir,d_mine_N_reservoir,d_N_fixer_C_reservoir,d_N_fixer_N_reservoir
  integer :: sp ! shorthand for current cohort specie
  integer :: i
  ! real :: Nbefore,Nafter
  real,dimension(N_C_TYPES) :: leaflitter_C,leaflitter_N,woodlitter_C,woodlitter_N,rootlitter_C,rootlitter_N
  real :: lim_factor
  real :: wood_n2c

  if(is_watch_point()) then
     write(*,*)'#### vegn_carbon_int ####'
     __DEBUG2__(soilt,theta)
     __DEBUG2__(vegn_tile_carbon(vegn),soil_tile_carbon(soil))
  endif

  ! Nbefore=soil_tile_nitrogen(soil)+vegn_tile_nitrogen(vegn)-(soil%gross_nitrogen_flux_into_tile-soil%gross_nitrogen_flux_out_of_tile)

total_scavenger_myc_C_allocated = 0.0
total_miner_myc_C_allocated = 0.0
total_mine_myc_immob = 0.0
total_scav_myc_immob = 0.0
total_N_fixation = 0.0
total_N_fixer_C_allocated = 0.0
total_mining_CO2prod = 0.0
total_myc_mine_C_uptake = 0.0
total_myc_resp = 0.0
total_myc_Nmin = 0.0
N_leakage = 0.0

  ! Do N deposition first. For now, it all goes to leaf litter
  if(ndep_amm > 0.0)then
		call soil_NH4_deposition(ndep_amm*dt_fast_yr,soil%leafLitter)
	endif

	if(ndep_nit > 0.0) then
		call soil_NO3_deposition(ndep_nit*dt_fast_yr,soil%leafLitter)
	endif

  if (ndep_org > 0.0) then
    call soil_org_N_deposition(ndep_org*dt_fast_yr,soil%leafLitter)
  endif

  soil%gross_nitrogen_flux_into_tile = soil%gross_nitrogen_flux_into_tile + (ndep_amm+ndep_nit+ndep_org)*dt_fast_yr

  !  update plant carbon
  vegn%npp = 0
  resp = 0 ; resl = 0 ; resr = 0 ; resg = 0 ; gpp = 0
  cgain = 0 ; closs = 0
  total_root_exudate_C = 0
  total_root_exudate_N = 0
  do i = 1, vegn%n_cohorts
     cc => vegn%cohorts(i)
     sp = cc%species

     call eddy_npp(cc,soilt);
     ! npp2 is for diagnostics and comparison
     cc%npp2 = cc%miami_npp;  ! treat miami npp as above+below npp

     if(dynamic_root_exudation .AND. soil_carbon_option==SOILC_CORPSE_N) then
       ! Initial allocation scheme: root exudation/mycorrhizal allocation depends on ratio of leaf biomass to max (as determined by N uptake)
       ! Root exudation fraction of NPP limited by some maximum value.  Probably need to rename these parameters and not use a hard-coded value
       root_exudate_frac = min(0.9,spdata(sp)%root_exudate_frac*cc%nitrogen_stress)
     else
       root_exudate_frac = spdata(sp)%root_exudate_frac
     endif
     root_exudate_C = max(cc%npp,0.0)*root_exudate_frac
     cc%carbon_gain = cc%carbon_gain + (cc%npp-root_exudate_C)*dt_fast_yr

     if(root_exudate_frac < 0) then
         __DEBUG4__(root_exudate_frac,cc%nitrogen_stress,cc%bl+cc%br,cc%npp)
         call error_mesg('vegn_carbon_int','root_exudate_frac<0',FATAL)
     endif


     ! check if leaves/roots are present and need to be accounted in maintenance
     if(cc%status == LEAF_ON) then
        md_alive = (cc%Pl * spdata(sp)%alpha(CMPT_LEAF) + &
                    cc%Pr * spdata(sp)%alpha(CMPT_ROOT))* &
              cc%bliving*dt_fast_yr
        md_leaf = (spdata(sp)%alpha(CMPT_LEAF)*cc%bl)*dt_fast_yr
        md_vleaf = spdata(sp)%alpha(CMPT_VLEAF)*cc%blv*dt_fast_yr ! This is zero during leaf-on unless there is some extra C in virtual leaves
        md_froot= spdata(sp)%alpha(CMPT_ROOT)*cc%br*dt_fast_yr
        ! NOTE that mathematically. md_alive = md_leaf + md_froot. Unfortunately,
        ! order of operation matters for the bit-wise reproducibility, so all
        ! three need to be calculated separately
     else
        md_alive = 0
        md_leaf  = 0
        md_vleaf = 0
	md_froot = 0
     endif

     ! compute branch and coarse wood losses for tree types
     md_wood =0;
     md_sapwood=0;
     if (sp > 1) then
        md_wood = 0.6 *cc%bwood * spdata(sp)%alpha(CMPT_WOOD)*dt_fast_yr;
        md_sapwood = 0.6 *cc%bsw * spdata(sp)%alpha(CMPT_SAPWOOD)*dt_fast_yr;
     endif
     md_alive = md_alive+md_sapwood

     ! this silliness with the mathematically equivalent formulas is
     ! solely to bit-reproduce old results in both CENTURY and CORPSE
     ! modes: the results are extremely sensitive to the order of operations,
     ! and diverge with time, even in stand-alone land model runs.
     select case (soil_carbon_option)
     case (SOILC_CENTURY, SOILC_CENTURY_BY_LAYER)
        cc%md = md_alive + cc%Psw_alphasw * cc%bliving * dt_fast_yr
    case (SOILC_CORPSE, SOILC_CORPSE_N)
        cc%md = md_leaf + md_vleaf + md_froot + md_sapwood + cc%Psw_alphasw * cc%bliving * dt_fast_yr

        leaflitter_C=(/spdata(sp)%fsc_liv *(md_leaf), (1-spdata(sp)%fsc_liv)*(md_leaf) ,0.0/)
        woodlitter_C=(/fsc_wood *md_wood*agf_bs, (1-fsc_wood)*md_wood*agf_bs, 0.0/)
        woodlitter_C = woodlitter_C + (/spdata(sp)%fsc_liv *md_sapwood*agf_bs, (1-spdata(sp)%fsc_liv)*md_sapwood*agf_bs, 0.0/)
        rootlitter_C=(/spdata(sp)%fsc_froot*md_froot + fsc_wood*md_wood*(1-agf_bs) + spdata(sp)%fsc_liv*md_sapwood*(1-agf_bs)+md_vleaf,&
 					(1-spdata(sp)%fsc_froot)*md_froot + (1-fsc_wood)*md_wood*(1-agf_bs) + (1-spdata(sp)%fsc_liv)*md_sapwood*(1-agf_bs),0.0/)
     end select

     cc%bwood_gain = cc%bwood_gain + cc%Psw_alphasw * cc%bliving * dt_fast_yr;

     ! This is causing nonconservation just below, at a point where md_wood>bwood_gain
     ! Fix: Reduce md_wood to match bwood_gain OR subtract from cc%bwood?
     cc%bwood_gain = cc%bwood_gain - md_wood;


   if(is_watch_point()) then
     __DEBUG2__(cc%bwood_gain,md_wood)
     __DEBUG2__(cc%carbon_gain,cc%md)
     __DEBUG5__(md_leaf , md_vleaf , md_froot , md_sapwood , cc%Psw_alphasw * cc%bliving * dt_fast_yr)
   endif

     if (cc%bwood_gain < 0.0) then
       call error_mesg('vegn_carbon_int','Warning: bwood_gain < 0 after subtracting md_wood',NOTE)
       ! cc%bwood_gain=0.0; BNS: just letting this be <0 for now to solve conservation problem ! potential non-conservation ?
     endif

     cc%carbon_gain = cc%carbon_gain - cc%md;
     cc%carbon_loss = cc%carbon_loss + cc%md; ! used in diagnostics only



     ! add md from leaf and root pools to fast soil carbon
     select case (soil_carbon_option)
     case (SOILC_CENTURY, SOILC_CENTURY_BY_LAYER)
        soil%fast_soil_C(1) = soil%fast_soil_C(1) +    fsc_liv *(md_alive) +    fsc_wood *md_wood;
        soil%slow_soil_C(1) = soil%slow_soil_C(1) + (1-fsc_liv)*(md_alive) + (1-fsc_wood)*md_wood;
     ! for budget tracking
!/*     cp->fsc_in(1)+= data->fsc_liv*md_alive+data->fsc_wood*md_wood; */
!/*     cp->ssc_in(1)+= (1.- data->fsc_liv)*md_alive+(1-data->fsc_wood)*md_wood; */
        soil%fsc_in(1)  = soil%fsc_in(1) + 1*md_alive+0*md_wood;
        soil%ssc_in(1)  = soil%ssc_in(1) + (1.- 1)*md_alive+(1-0)*md_wood;
    case (SOILC_CORPSE_N)
      if(cc%bwood>0) then
        wood_n2c = cc%wood_N/cc%bwood
      else
        wood_n2c=0.0
      endif

      ! Should this N be lost or retranslocated?
      cc%leaf_N = cc%leaf_N - md_leaf/spdata(sp)%leaf_live_c2n*(1.0-spdata(sp)%leaf_retranslocation_frac)
      cc%root_N = cc%root_N - md_froot/spdata(sp)%froot_live_c2n*(1.0-spdata(sp)%froot_retranslocation_frac)
      cc%wood_N = cc%wood_N - md_wood*wood_n2c
      cc%sapwood_N = cc%sapwood_N - md_sapwood/spdata(sp)%sapwood_c2n*(1.0-spdata(sp)%froot_retranslocation_frac)

       leaflitter_N=(/spdata(sp)%fsc_liv *md_leaf/spdata(sp)%leaf_live_c2n*(1.0-spdata(sp)%leaf_retranslocation_frac), &
                  (1-spdata(sp)%fsc_liv)*md_leaf/spdata(sp)%leaf_live_c2n*(1.0-spdata(sp)%leaf_retranslocation_frac) ,0.0/)

       woodlitter_N=(/fsc_wood *md_wood*agf_bs*wood_n2c, (1-fsc_wood)*md_wood*agf_bs*wood_n2c, 0.0/)
       woodlitter_N = woodlitter_N + (/spdata(sp)%fsc_liv *md_sapwood*agf_bs/spdata(sp)%sapwood_c2n*(1.0-spdata(sp)%froot_retranslocation_frac), &
               (1-spdata(sp)%fsc_liv)*md_sapwood*agf_bs/spdata(sp)%sapwood_c2n*(1.0-spdata(sp)%froot_retranslocation_frac), 0.0/)

        rootlitter_N=(/spdata(sp)%fsc_froot*md_froot/spdata(sp)%froot_live_c2n*(1.0-spdata(sp)%froot_retranslocation_frac) &
             + fsc_wood*md_wood*(1-agf_bs)*wood_n2c + spdata(sp)%fsc_liv*md_sapwood*(1-agf_bs)/spdata(sp)%sapwood_c2n*(1.0-spdata(sp)%froot_retranslocation_frac), &
            (1-spdata(sp)%fsc_froot)*md_froot/spdata(sp)%froot_live_c2n*(1.0-spdata(sp)%froot_retranslocation_frac)&
             + (1-fsc_wood)*md_wood*(1-agf_bs)*wood_n2c + (1-spdata(sp)%fsc_liv)*md_sapwood*(1-agf_bs)/spdata(sp)%sapwood_c2n*(1.0-spdata(sp)%froot_retranslocation_frac),0.0/)


        if (is_watch_point()) then
           call debug_pool(soil%leafLitter,      'leafLitter (before)'      )
           call debug_pool(soil%coarseWoodLitter,'coarseWoodLitter (before)')
        endif
        call add_litter(soil%leafLitter,leaflitter_C,leaflitter_N)
        call add_litter(soil%coarseWoodLitter,woodlitter_C,woodlitter_N)
        soil%leaflitter_fsc_in=soil%leaflitter_fsc_in+leaflitter_C(1)
        soil%leaflitter_ssc_in=soil%leaflitter_ssc_in+leaflitter_C(2)
        soil%coarsewoodlitter_fsc_in=soil%coarsewoodlitter_fsc_in + woodlitter_C(1)
        soil%coarsewoodlitter_ssc_in=soil%coarsewoodlitter_ssc_in + woodlitter_C(2)
        soil%leaflitter_fsn_in=soil%leaflitter_fsn_in+leaflitter_N(1)
        soil%leaflitter_ssn_in=soil%leaflitter_ssn_in+leaflitter_N(2)
        soil%coarsewoodlitter_fsn_in=soil%coarsewoodlitter_fsn_in + woodlitter_N(1)
        soil%coarsewoodlitter_ssn_in=soil%coarsewoodlitter_ssn_in + woodlitter_N(2)
	!ssc_in and fsc_in updated in add_root_litter
        call add_root_litter(soil,vegn,rootlitter_C,rootlitter_N)
	if (is_watch_point()) then
           call debug_pool(soil%leafLitter,      'leafLitter (after)'      )
           call debug_pool(soil%coarseWoodLitter,'coarseWoodLitter (after)')
        endif
    case (SOILC_CORPSE)
        if (is_watch_point()) then
            call debug_pool(soil%leafLitter,      'leafLitter (before)'      )
            call debug_pool(soil%coarseWoodLitter,'coarseWoodLitter (before)')
        endif
        call add_litter(soil%leafLitter,leaflitter_C,(/0.0,0.0 ,0.0/))
        call add_litter(soil%coarseWoodLitter,woodlitter_C,(/0.0,0.0, 0.0/))
        soil%leaflitter_fsc_in=soil%leaflitter_fsc_in+leaflitter_C(1)
        soil%leaflitter_ssc_in=soil%leaflitter_ssc_in+leaflitter_C(2)
        soil%coarsewoodlitter_fsc_in=soil%coarsewoodlitter_fsc_in +  woodlitter_C(1)
        soil%coarsewoodlitter_ssc_in=soil%coarsewoodlitter_ssc_in + woodlitter_C(2)

        !ssc_in and fsc_in updated in add_root_litter
        call add_root_litter(soil,vegn,rootlitter_C,(/0.0,0.0,0.0/))
        if (is_watch_point()) then
            call debug_pool(soil%leafLitter,      'leafLitter (after)'      )
            call debug_pool(soil%coarseWoodLitter,'coarseWoodLitter (after)')
        endif

     case default
        call error_mesg('vegn_carbon_int','The value of soil_carbon_option is invalid. This should never happen. Contact developer.',FATAL)
     end select

 if(is_watch_point()) then
   write(*,*) 'After md'
   __DEBUG2__(soil_tile_carbon(soil),vegn_tile_carbon(vegn))
 endif

     vegn%veg_in  = vegn%veg_in  + cc%npp*dt_fast_yr;
     vegn%veg_out = vegn%veg_out + md_alive+md_wood;

     if(is_watch_point()) then
        __DEBUG4__(md_wood,md_leaf,md_froot,md_sapwood)
        __DEBUG4__(cc%bl, cc%br, cc%bsw, cc%bwood)
        __DEBUG3__(cc%An_op, cc%An_cl, cc%lai)
        __DEBUG1__(cc%species)
        __DEBUG2__(cc%npp, cc%gpp)
        __DEBUG4__(cc%resp, cc%resl, cc%resr, cc%resg)
        __DEBUG2__(cc%carbon_gain, cc%bwood_gain)
     endif
     ! accumulate tile-level NPP and GPP
     vegn%npp = vegn%npp + cc%npp
     gpp = gpp + cc%gpp
     ! accumulate respiration terms for tile-level reporting
     resp = resp + cc%resp ; resl = resl + cc%resl
     resr = resr + cc%resr ; resg = resg + cc%resg
     ! accumulate gain/loss terms for tile-level reporting
     cgain = cgain + cc%carbon_gain
     closs = closs + cc%carbon_loss

     ! Mycorrhizal N uptake
     if(soil_carbon_option == SOILC_CORPSE_N) then

       ! Scavenging (AM-style)
       call myc_scavenger_N_uptake(soil,vegn,cc%myc_scavenger_biomass_C,myc_scav_N_uptake,dt_fast_yr,.TRUE.)
       cc%scav_myc_N_reservoir = cc%scav_myc_N_reservoir+myc_scav_N_uptake
      ! Mycorrhizal N mining (ECM-style)
       call myc_miner_N_uptake(soil,vegn,cc%myc_miner_biomass_C,myc_mine_N_uptake,myc_mine_C_uptake,mining_CO2prod,dt_fast_yr,.TRUE.)
       cc%mine_myc_C_reservoir = cc%mine_myc_C_reservoir+myc_mine_C_uptake
       cc%mine_myc_N_reservoir = cc%mine_myc_N_reservoir+myc_mine_N_uptake
      ! Root uptake of nitrogen
       call root_N_uptake(soil,vegn,root_active_N_uptake,dt_fast_yr,.TRUE.)
      ! N fixation
       N_fixation = cc%N_fixer_biomass_C*N_fixation_rate*dt_fast_yr
       cc%N_fixer_N_reservoir = cc%N_fixer_N_reservoir + N_fixation

            if(cc%myc_scavenger_biomass_C<0) call error_mesg('vegn_carbon_int','Mycorrhizal scavenger biomass < 0',FATAL)
            if(cc%myc_miner_biomass_C<0) then
              __DEBUG3__(cc%myc_miner_biomass_C,myc_mine_C_uptake,miner_myc_C_allocated)
              call error_mesg('vegn_carbon_int','Mycorrhizal miner biomass < 0',FATAL)
            endif
            if(cc%N_fixer_biomass_C<0) call error_mesg('vegn_carbon_int','N fixer biomass < 0',FATAL)

            if(isnan(cc%myc_scavenger_biomass_C)) call error_mesg('vegn_carbon_int','Mycorrhizal scavenger biomass is NaN',FATAL)
            if(isnan(cc%myc_miner_biomass_C)) call error_mesg('vegn_carbon_int','Mycorrhizal miner biomass is NaN',FATAL)
            if(isnan(cc%N_fixer_biomass_C)) call error_mesg('vegn_carbon_int','N fixer biomass is NaN',FATAL)

            if(is_watch_point()) then
               __DEBUG3__(cc%myc_scavenger_biomass_C,cc%myc_miner_biomass_C,cc%N_fixer_biomass_C)
            endif


            ! First add mycorrhizal and N fixer turnover to soil C pools
            ! To do: put reservoir turnover into root exudates
            call add_root_litter(soil,vegn,(/0.0,cc%myc_scavenger_biomass_C*(deadmic_slow_frac),cc%myc_scavenger_biomass_C*(1.0-deadmic_slow_frac)/)*dt_fast_yr*et_myc/mycorrhizal_turnover_time,&
                        (/0.0,cc%myc_scavenger_biomass_N*deadmic_slow_frac,cc%myc_scavenger_biomass_N*(1.0-deadmic_slow_frac)/)*dt_fast_yr*et_myc/mycorrhizal_turnover_time)
            call add_root_litter(soil,vegn,(/0.0,cc%myc_miner_biomass_C*deadmic_slow_frac,(cc%myc_miner_biomass_C*(1.0-deadmic_slow_frac))/)*dt_fast_yr*et_myc/mycorrhizal_turnover_time,&
                         (/0.0,cc%myc_miner_biomass_N*deadmic_slow_frac,(cc%myc_miner_biomass_N*(1.0-deadmic_slow_frac))/)*dt_fast_yr*et_myc/mycorrhizal_turnover_time)
            call add_root_litter(soil,vegn,(/0.0,0.0,(cc%N_fixer_biomass_C)/N_fixer_turnover_time/)*dt_fast_yr*et_myc,&
                         (/0.0,0.0,(cc%N_fixer_biomass_N)/N_fixer_turnover_time/)*dt_fast_yr*et_myc)

            total_myc_resp = total_myc_resp + (1.0-et_myc)*dt_fast_yr*&
                 ((cc%myc_scavenger_biomass_C+cc%myc_miner_biomass_C)/mycorrhizal_turnover_time+&
                 (cc%N_fixer_biomass_C)/N_fixer_turnover_time)

            ! total_myc_Nmin = total_myc_Nmin + (1.0-et_myc)*dt_fast_yr*&
                !  ((cc%myc_scavenger_biomass_N+cc%scav_myc_N_reservoir+cc%myc_miner_biomass_N+cc%mine_myc_N_reservoir)/mycorrhizal_turnover_time+&
                !  (cc%N_fixer_biomass_N+cc%N_fixer_N_reservoir)/N_fixer_turnover_time)

            ! d_scav_C_reservoir = - cc%scav_myc_C_reservoir/mycorrhizal_turnover_time*dt_fast_yr
            ! d_scav_N_reservoir = - cc%scav_myc_N_reservoir/mycorrhizal_turnover_time*dt_fast_yr
            ! d_mine_C_reservoir = - cc%mine_myc_C_reservoir/mycorrhizal_turnover_time*dt_fast_yr
            ! d_mine_N_reservoir = - cc%mine_myc_N_reservoir/mycorrhizal_turnover_time*dt_fast_yr
            ! d_N_fixer_C_reservoir = - cc%N_fixer_C_reservoir/N_fixer_turnover_time*dt_fast_yr
            ! d_N_fixer_N_reservoir = - cc%N_fixer_N_reservoir/N_fixer_turnover_time*dt_fast_yr
            d_scav_C_reservoir = 0.0
            d_scav_N_reservoir = 0.0
            d_mine_C_reservoir = 0.0
            d_mine_N_reservoir = 0.0
            d_N_fixer_C_reservoir = 0.0
            d_N_fixer_N_reservoir = 0.0

            ! Then update biomass of scavenger mycorrhizae and N fixers
            scavenger_myc_growth = myc_growth_rate*cc%scav_myc_C_reservoir/(cc%scav_myc_C_reservoir+kM_myc_growth)*myc_scav_C_efficiency*dt_fast_yr
            maint_resp=min(cc%myc_scavenger_biomass_C/mycorrhizal_turnover_time*(1.0-et_myc)*dt_fast_yr,scavenger_myc_growth)
            if (scavenger_myc_growth-maint_resp>c2n_mycorrhizae*cc%scav_myc_N_reservoir*0.9) then
              ! Not enough nitrogen to support growth. Limit to available N, and leave a little bit left over for plant
              scavenger_myc_growth=c2n_mycorrhizae*cc%scav_myc_N_reservoir*0.9+maint_resp
            endif
            total_myc_resp = total_myc_resp + scavenger_myc_growth/myc_scav_C_efficiency*(1.0-myc_scav_C_efficiency)

            cc%scav_myc_N_reservoir=cc%scav_myc_N_reservoir+cc%myc_scavenger_biomass_N*(1-et_myc/mycorrhizal_turnover_time*dt_fast_yr)
            cc%myc_scavenger_biomass_C = cc%myc_scavenger_biomass_C + scavenger_myc_growth - cc%myc_scavenger_biomass_C/mycorrhizal_turnover_time*dt_fast_yr
            ! cc%myc_scavenger_biomass_N = cc%myc_scavenger_biomass_N + (scavenger_myc_growth-maint_resp)/c2n_mycorrhizae - cc%myc_scavenger_biomass_N/mycorrhizal_turnover_time*et_myc*dt_fast_yr
            cc%myc_scavenger_biomass_N = cc%myc_scavenger_biomass_C/c2n_mycorrhizae
            d_scav_C_reservoir = d_scav_C_reservoir - scavenger_myc_growth/myc_scav_C_efficiency
            ! d_scav_N_reservoir = d_scav_N_reservoir - (scavenger_myc_growth-maint_resp)/c2n_mycorrhizae
            cc%scav_myc_N_reservoir=cc%scav_myc_N_reservoir-cc%myc_scavenger_biomass_N

            miner_myc_growth = myc_growth_rate*cc%mine_myc_C_reservoir/(cc%mine_myc_C_reservoir+kM_myc_growth)*myc_mine_C_efficiency*dt_fast_yr
            maint_resp=min(cc%myc_miner_biomass_C/mycorrhizal_turnover_time*(1.0-et_myc)*dt_fast_yr,miner_myc_growth)
            if (miner_myc_growth-maint_resp>c2n_mycorrhizae*cc%mine_myc_N_reservoir*0.9) then
              ! Not enough nitrogen to support growth. Limit to available N, and leave a little bit left over for plant
              miner_myc_growth=c2n_mycorrhizae*cc%mine_myc_N_reservoir*0.9+maint_resp
            endif
            total_myc_resp = total_myc_resp + miner_myc_growth/myc_mine_C_efficiency*(1.0-myc_mine_C_efficiency)

            cc%mine_myc_N_reservoir=cc%mine_myc_N_reservoir+cc%myc_miner_biomass_N*(1-et_myc/mycorrhizal_turnover_time*dt_fast_yr)
            cc%myc_miner_biomass_C = cc%myc_miner_biomass_C + miner_myc_growth - cc%myc_miner_biomass_C/mycorrhizal_turnover_time*dt_fast_yr
            ! cc%myc_miner_biomass_N = cc%myc_miner_biomass_N + (miner_myc_growth-maint_resp)/c2n_mycorrhizae - cc%myc_miner_biomass_N/mycorrhizal_turnover_time*et_myc*dt_fast_yr
            cc%myc_miner_biomass_N = cc%myc_miner_biomass_C/c2n_mycorrhizae
            d_mine_C_reservoir = d_mine_C_reservoir - miner_myc_growth/myc_mine_C_efficiency
            ! d_mine_N_reservoir = d_mine_N_reservoir - (miner_myc_growth-maint_resp)/c2n_mycorrhizae
            cc%mine_myc_N_reservoir=cc%mine_myc_N_reservoir-cc%myc_miner_biomass_N

            N_fixer_growth = myc_growth_rate*cc%N_fixer_C_reservoir/(cc%N_fixer_C_reservoir+kM_myc_growth)*N_fixer_C_efficiency*dt_fast_yr
            maint_resp=min(cc%N_fixer_biomass_C/N_fixer_turnover_time*(1.0-et_myc)*dt_fast_yr,N_fixer_growth)
            ! if (N_fixer_growth>c2n_mycorrhizae*cc%N_fixer_N_reservoir*0.9) then
            !   ! Not enough nitrogen to support growth. Limit to available N, and leave a little bit left over for plant
            !   N_fixer_growth=c2n_mycorrhizae*cc%N_fixer_N_reservoir*0.9
            ! endif

            total_myc_resp = total_myc_resp + N_fixer_growth/N_fixer_C_efficiency*(1.0-N_fixer_C_efficiency)

            N_fixation=N_fixation-cc%N_fixer_biomass_N*(1-et_myc/N_fixer_turnover_time*dt_fast_yr)
            cc%N_fixer_biomass_C = cc%N_fixer_biomass_C + N_fixer_growth - cc%N_fixer_biomass_C/N_fixer_turnover_time*dt_fast_yr
            ! cc%N_fixer_biomass_N = cc%N_fixer_biomass_N + (N_fixer_growth-maint_resp)/c2n_N_fixer - cc%N_fixer_biomass_N/N_fixer_turnover_time*et_myc*dt_fast_yr
            cc%N_fixer_biomass_N = cc%N_fixer_biomass_N/c2n_N_fixer
            d_N_fixer_C_reservoir = d_N_fixer_C_reservoir - (N_fixer_growth)/N_fixer_C_efficiency
            ! d_N_fixer_N_reservoir = d_N_fixer_N_reservoir - N_fixer_growth/c2n_mycorrhizae

            ! N fixers just make all the N they need for their biomass
            ! N_fixation = N_fixation + (N_fixer_growth-maint_resp)/c2n_mycorrhizae
            N_fixation=N_fixation+cc%N_fixer_biomass_N

            if(cc%myc_scavenger_biomass_C<0) call error_mesg('vegn_carbon_int','Mycorrhizal scavenger biomass < 0',FATAL)
            if(cc%myc_miner_biomass_C<0) then
              __DEBUG3__(cc%myc_miner_biomass_C,myc_mine_C_uptake,miner_myc_C_allocated)
              call error_mesg('vegn_carbon_int','Mycorrhizal miner biomass < 0',FATAL)
            endif
            if(cc%N_fixer_biomass_C<0) call error_mesg('vegn_carbon_int','N fixer biomass < 0',FATAL)

            if(cc%myc_scavenger_biomass_N<0) call error_mesg('vegn_carbon_int','Mycorrhizal scavenger biomass N < 0',FATAL)
            if(cc%myc_miner_biomass_N<0) then
              __DEBUG3__(cc%myc_miner_biomass_N,myc_mine_N_uptake,miner_myc_N_allocated)
              call error_mesg('vegn_carbon_int','Mycorrhizal miner biomass N < 0',FATAL)
            endif
            if(cc%N_fixer_biomass_N<0) call error_mesg('vegn_carbon_int','N fixer biomass N < 0',FATAL)



            if(cc%N_fixer_biomass_C<0) call error_mesg('vegn_carbon_int','N fixer biomass < 0',FATAL)
            if(isnan(cc%myc_scavenger_biomass_C)) call error_mesg('vegn_carbon_int','Mycorrhizal scavenger biomass is NaN',FATAL)
            if(isnan(cc%myc_miner_biomass_C)) call error_mesg('vegn_carbon_int','Mycorrhizal miner biomass is NaN',FATAL)
            if(isnan(cc%N_fixer_biomass_C)) call error_mesg('vegn_carbon_int','N fixer biomass is NaN',FATAL)

            if(isnan(cc%myc_scavenger_biomass_N)) call error_mesg('vegn_carbon_int','Mycorrhizal scavenger biomass is NaN',FATAL)
            if(isnan(cc%myc_miner_biomass_N)) call error_mesg('vegn_carbon_int','Mycorrhizal miner biomass is NaN',FATAL)
            if(isnan(cc%N_fixer_biomass_N)) call error_mesg('vegn_carbon_int','N fixer biomass is NaN',FATAL)

            cc%scav_myc_C_reservoir = cc%scav_myc_C_reservoir + d_scav_C_reservoir
            cc%scav_myc_N_reservoir = cc%scav_myc_N_reservoir + d_scav_N_reservoir
            cc%mine_myc_C_reservoir = cc%mine_myc_C_reservoir + d_mine_C_reservoir
            cc%mine_myc_N_reservoir = cc%mine_myc_N_reservoir + d_mine_N_reservoir
            cc%N_fixer_C_reservoir = cc%N_fixer_C_reservoir + d_N_fixer_C_reservoir
            cc%N_fixer_N_reservoir = cc%N_fixer_N_reservoir + d_N_fixer_N_reservoir

            ! Excess C leaks out of reservoir into root exudates at a time scale of one day
            total_root_exudate_C = total_root_exudate_C + (cc%scav_myc_C_reservoir + cc%mine_myc_C_reservoir + cc%N_fixer_C_reservoir)*dt_fast_yr*365
            cc%scav_myc_C_reservoir = cc%scav_myc_C_reservoir - cc%scav_myc_C_reservoir*dt_fast_yr*365
            cc%mine_myc_C_reservoir = cc%mine_myc_C_reservoir - cc%mine_myc_C_reservoir*dt_fast_yr*365
            cc%N_fixer_C_reservoir = cc%N_fixer_C_reservoir - cc%N_fixer_C_reservoir*dt_fast_yr*365

            if(abs(cc%mine_myc_N_reservoir)<1e-20) cc%mine_myc_N_reservoir = 0.0
            if(abs(cc%mine_myc_C_reservoir)<1e-20) cc%mine_myc_C_reservoir = 0.0
            if(abs(cc%scav_myc_N_reservoir)<1e-20) cc%scav_myc_N_reservoir = 0.0
            if(abs(cc%scav_myc_C_reservoir)<1e-20) cc%scav_myc_C_reservoir = 0.0

            if(cc%scav_myc_C_reservoir<-1e-10) call error_mesg('vegn_carbon_int','Mycorrhizal scavenger C reservoir < 0',FATAL)
            if(cc%mine_myc_C_reservoir<-1e-10) call error_mesg('vegn_carbon_int','Mycorrhizal miner C reservoir < 0',FATAL)
            if(cc%N_fixer_C_reservoir<-1e-10) call error_mesg('vegn_carbon_int','N fixer C reservoir < 0',FATAL)
            if(cc%scav_myc_N_reservoir<-1e-10) then
                __DEBUG4__(cc%scav_myc_N_reservoir,cc%scav_myc_C_reservoir,cc%myc_scavenger_biomass_C,cc%myc_scavenger_biomass_N)
                __DEBUG4__(cc%species,cc%total_N,cc%leaf_N,cc%stored_N)
                call error_mesg('vegn_carbon_int','Mycorrhizal scavenger N reservoir < 0',FATAL)
            endif
            if(cc%mine_myc_N_reservoir<-1e-10) then
                __DEBUG4__(cc%mine_myc_N_reservoir,cc%mine_myc_C_reservoir,cc%myc_miner_biomass_C,cc%myc_miner_biomass_N)
                __DEBUG4__(cc%species,cc%total_N,cc%leaf_N,cc%stored_N)
                call error_mesg('vegn_carbon_int','Mycorrhizal miner N reservoir < 0',FATAL)
            endif
            if(cc%N_fixer_N_reservoir<-1e-10) call error_mesg('vegn_carbon_int','N fixer N reservoir < 0',FATAL)

            ! Calculate N released to plant
            scav_N_to_plant = cc%scav_myc_N_reservoir*myc_N_to_plant_rate*dt_fast_yr
            mine_N_to_plant = cc%mine_myc_N_reservoir*myc_N_to_plant_rate*dt_fast_yr
            fix_N_to_plant = cc%N_fixer_N_reservoir*myc_N_to_plant_rate*dt_fast_yr

         if(myc_scav_C_efficiency == 0 .OR. .NOT. spdata(sp)%do_N_scavenging_strategy) then
           myc_scav_marginal_gain = 0
           scav_N_to_plant = 0.0
         elseif (cc%myc_scavenger_biomass_C > 1e-10) then
           myc_scav_marginal_gain = (max(0.0,scav_N_to_plant)/dt_fast_yr)/(cc%myc_scavenger_biomass_C/myc_scav_C_efficiency/mycorrhizal_turnover_time)
         else
           call myc_scavenger_N_uptake(soil,vegn,0.001,myc_N_uptake,dt_fast_yr,.FALSE.)
           myc_scav_marginal_gain = (myc_N_uptake/dt_fast_yr)/(0.001/myc_scav_C_efficiency/mycorrhizal_turnover_time)
         endif
         cc%scav_myc_N_reservoir = cc%scav_myc_N_reservoir - scav_N_to_plant


         if(myc_mine_C_efficiency==0 .OR. .NOT. spdata(sp)%do_N_mining_strategy) then
           myc_mine_marginal_gain = 0.0
           mine_N_to_plant = 0.0
         elseif(cc%myc_miner_biomass_C>1e-10) then
           myc_mine_marginal_gain = (max(0.0,mine_N_to_plant)/dt_fast_yr)/(cc%myc_miner_biomass_C/myc_mine_C_efficiency/mycorrhizal_turnover_time)
         else
           call myc_miner_N_uptake(soil,vegn,0.001,myc_N_uptake,myc_C_uptake,dummy1,dt_fast_yr,.FALSE.)
           myc_mine_marginal_gain = (myc_N_uptake/dt_fast_yr)/(0.001/myc_mine_C_efficiency/mycorrhizal_turnover_time)
         endif
         cc%mine_myc_N_reservoir = cc%mine_myc_N_reservoir - mine_N_to_plant

         if (cc%br>0) then
           rhiz_exud_marginal_gain = (root_active_N_uptake/dt_fast_yr)/(cc%br)!+(myc_mine_marginal_gain+myc_scav_marginal_gain)*0.5
         else
           rhiz_exud_marginal_gain = (myc_mine_marginal_gain+myc_scav_marginal_gain)*0.25
         endif

         if(N_fixer_C_efficiency == 0 .OR. .NOT. spdata(sp)%do_N_fixation_strategy) then
          N_fix_marginal_gain = 0.0
          fix_N_to_plant = 0.0
         elseif(cc%N_fixer_biomass_C>1e-10) then
          !  N_fix_marginal_gain = (N_fixation_2-N_fixation)/dt_fast_yr/(0.05*root_exudate_C*dt_fast_yr)
          N_fix_marginal_gain = (fix_N_to_plant/dt_fast_yr)/(cc%N_fixer_biomass_C/N_fixer_C_efficiency/N_fixer_turnover_time)
         else
          N_fix_marginal_gain = (0.001*N_fixation_rate)/(0.001/N_fixer_C_efficiency/N_fixer_turnover_time)
         endif
         cc%N_fixer_N_reservoir = cc%N_fixer_N_reservoir-fix_N_to_plant

         ! Apply a smoothing filter to marginal gains, so we can control how fast N strategies change at the ecosystem level
         weight_marginal_gain = 1/(1+spdata(sp)%tau_smooth_marginal_gain/dt_fast_yr)
         cc%myc_scav_marginal_gain_smoothed=cc%myc_scav_marginal_gain_smoothed*(1.0-weight_marginal_gain) + myc_scav_marginal_gain*weight_marginal_gain
         cc%myc_mine_marginal_gain_smoothed=cc%myc_mine_marginal_gain_smoothed*(1.0-weight_marginal_gain) + myc_mine_marginal_gain*weight_marginal_gain
         cc%N_fix_marginal_gain_smoothed=cc%N_fix_marginal_gain_smoothed*(1.0-weight_marginal_gain) + N_fix_marginal_gain*weight_marginal_gain
         cc%rhiz_exud_marginal_gain_smoothed=cc%rhiz_exud_marginal_gain_smoothed*(1.0-weight_marginal_gain) + rhiz_exud_marginal_gain*weight_marginal_gain

         ! Calculate relative fractions
          if (cc%myc_scav_marginal_gain_smoothed+cc%N_fix_marginal_gain_smoothed+cc%myc_mine_marginal_gain_smoothed+cc%rhiz_exud_marginal_gain_smoothed>0) then
            myc_scav_exudate_frac = cc%myc_scav_marginal_gain_smoothed/(cc%myc_scav_marginal_gain_smoothed+cc%myc_mine_marginal_gain_smoothed+cc%N_fix_marginal_gain_smoothed+cc%rhiz_exud_marginal_gain_smoothed)
            myc_mine_exudate_frac = cc%myc_mine_marginal_gain_smoothed/(cc%myc_scav_marginal_gain_smoothed+cc%myc_mine_marginal_gain_smoothed+cc%N_fix_marginal_gain_smoothed+cc%rhiz_exud_marginal_gain_smoothed)
            N_fixer_exudate_frac = cc%N_fix_marginal_gain_smoothed/(cc%myc_scav_marginal_gain_smoothed+cc%myc_mine_marginal_gain_smoothed+cc%N_fix_marginal_gain_smoothed+cc%rhiz_exud_marginal_gain_smoothed)
            rhiz_exud_frac = cc%rhiz_exud_marginal_gain_smoothed/(cc%myc_scav_marginal_gain_smoothed+cc%myc_mine_marginal_gain_smoothed+cc%N_fix_marginal_gain_smoothed+cc%rhiz_exud_marginal_gain_smoothed)
          else
            ! Divide evenly if there is no marginal gain.  But this probably only happens if root_exudate_C is zero
            myc_scav_exudate_frac = 0.4*0.7
            myc_mine_exudate_frac = 0.3*0.7
            N_fixer_exudate_frac = 0.3*0.7
            rhiz_exud_frac = 0.3
          endif

          weight_alloc = 1/(1+spdata(sp)%tau_smooth_alloc/dt_fast_yr)
          cc%N_fix_alloc_smoothed = cc%N_fix_alloc_smoothed*(1.0-weight_alloc) + root_exudate_C*N_fixer_exudate_frac*dt_fast_yr*weight_alloc
          cc%myc_mine_alloc_smoothed = cc%myc_mine_alloc_smoothed*(1.0-weight_alloc) + root_exudate_C*myc_mine_exudate_frac*dt_fast_yr*weight_alloc
          cc%myc_scav_alloc_smoothed = cc%myc_scav_alloc_smoothed*(1.0-weight_alloc) + root_exudate_C*myc_scav_exudate_frac*dt_fast_yr*weight_alloc

          if(smooth_N_uptake_C_allocation) then
            N_fixer_C_allocated = min(root_exudate_C*N_fixer_exudate_frac*dt_fast_yr,cc%N_fix_alloc_smoothed)
            miner_myc_C_allocated=min(root_exudate_C*myc_mine_exudate_frac*dt_fast_yr,cc%myc_mine_alloc_smoothed)
          else
            N_fixer_C_allocated = root_exudate_C*N_fixer_exudate_frac*dt_fast_yr
            miner_myc_C_allocated=root_exudate_C*myc_mine_exudate_frac*dt_fast_yr
          endif
          scavenger_myc_C_allocated=root_exudate_C*myc_scav_exudate_frac*dt_fast_yr

          N_fixer_N_allocated = 0.0
          miner_myc_N_allocated = miner_myc_C_allocated*root_exudate_N_frac
          scavenger_myc_N_allocated = scavenger_myc_C_allocated*root_exudate_N_frac

          ! Make sure N allocation doesn't completely deplete stored N
          if (cc%stored_N<=0.0) then
             N_fixer_N_allocated=0.0
             miner_myc_N_allocated=0.0
             scavenger_myc_N_allocated=0.0
          elseif (N_fixer_N_allocated+miner_myc_N_allocated+scavenger_myc_N_allocated > 0.0 .AND. &
                    N_limits_live_biomass .AND. &
                    N_fixer_N_allocated+miner_myc_N_allocated+scavenger_myc_N_allocated > cc%stored_N*0.9) then
             lim_factor=cc%stored_N/(N_fixer_N_allocated+miner_myc_N_allocated+scavenger_myc_N_allocated)*0.9
             N_fixer_N_allocated=N_fixer_N_allocated*lim_factor
             miner_myc_N_allocated=miner_myc_N_allocated*lim_factor
             scavenger_myc_N_allocated=scavenger_myc_N_allocated*lim_factor
          endif

          cc%scav_myc_C_reservoir = cc%scav_myc_C_reservoir + scavenger_myc_C_allocated
          cc%scav_myc_N_reservoir = cc%scav_myc_N_reservoir + scavenger_myc_N_allocated
          cc%mine_myc_C_reservoir = cc%mine_myc_C_reservoir + miner_myc_C_allocated
          cc%mine_myc_N_reservoir = cc%mine_myc_N_reservoir + miner_myc_N_allocated
          cc%N_fixer_C_reservoir = cc%N_fixer_C_reservoir + N_fixer_C_allocated
          cc%N_fixer_N_reservoir = cc%N_fixer_N_reservoir + N_fixer_N_allocated

          if (root_exudate_C>0 .AND. is_watch_point()) then
            print *,'allocation'
            __DEBUG4__(myc_scav_exudate_frac,myc_mine_exudate_frac,N_fixer_exudate_frac,rhiz_exud_frac)
          endif

          root_exudate_N=(root_exudate_C*root_exudate_N_frac*dt_fast_yr - scavenger_myc_N_allocated - N_fixer_N_allocated - miner_myc_N_allocated)/dt_fast_yr
          if(root_exudate_N<0) then
             __DEBUG5__(root_exudate_N*dt_fast_yr,scavenger_myc_N_allocated,N_fixer_N_allocated,miner_myc_N_allocated,miner_myc_C_allocated)
          endif

     else
         scavenger_myc_C_allocated=0.0
         miner_myc_C_allocated=0.0
         scavenger_myc_N_allocated=0.0
         miner_myc_N_allocated=0.0
         myc_scav_N_uptake=0.0
         myc_mine_N_uptake=0.0
         scav_N_to_plant=0.0
         mine_N_to_plant=0.0
         fix_N_to_plant=0.0
         mycorrhizal_scav_N_immob = 0.0
         mycorrhizal_mine_N_immob = 0.0
         N_fixer_C_allocated = 0.0
         root_active_N_uptake = 0.0
         N_fixation = 0.0
         mining_CO2prod = 0.0
         myc_mine_C_uptake = 0.0
         root_exudate_N = 0.0
     endif


 if(is_watch_point()) then
   write(*,*) 'After md, before mycorrhizal turnover'
   __DEBUG2__(soil_tile_carbon(soil),vegn_tile_carbon(vegn))
 endif

     total_plant_N_uptake = scav_N_to_plant + mine_N_to_plant + fix_N_to_plant + root_active_N_uptake
     cc%stored_N = cc%stored_N + total_plant_N_uptake

    !  __DEBUG3__(scavenger_myc_N_allocated,N_fixer_N_allocated,miner_myc_N_allocated)

    !  root_exudate_N = 0.0

     ! To prevent excess stored N buildup under high soil N, leak stored N when stress is zero
     ! NOTE: nitrogen_stress has a lower limit in vegn_cohort, so this logic should be fixed somehow so it doesn't need to be manually kept in sync
     if(cc%nitrogen_stress<=0.05 .AND. cc%stored_N>0 .AND. soil_carbon_option == SOILC_CORPSE_N) then
       N_leakage = N_leakage + cc%stored_N*excess_stored_N_leakage_rate*dt_fast_yr
       cc%stored_N=cc%stored_N - cc%stored_N*excess_stored_N_leakage_rate*dt_fast_yr
     endif

     cc%stored_N=cc%stored_N - root_exudate_N*dt_fast_yr - scavenger_myc_N_allocated - N_fixer_N_allocated - miner_myc_N_allocated


     total_root_exudate_C = total_root_exudate_C + root_exudate_C*dt_fast_yr - scavenger_myc_C_allocated - N_fixer_C_allocated - miner_myc_C_allocated
     total_root_exudate_N = total_root_exudate_N + root_exudate_N*dt_fast_yr

     if(total_root_exudate_C < 0.0) then
         __DEBUG4__(root_exudate_C*dt_fast_yr,scavenger_myc_C_allocated,N_fixer_C_allocated,miner_myc_C_allocated)
         call error_mesg('vegn_carbon_int','Root exudate C < 0',FATAL)
     endif

     total_scavenger_myc_C_allocated = total_scavenger_myc_C_allocated + scavenger_myc_C_allocated
     total_miner_myc_C_allocated = total_miner_myc_C_allocated + miner_myc_C_allocated
     total_scav_myc_immob = total_scav_myc_immob + mycorrhizal_scav_N_immob
     total_mine_myc_immob = total_mine_myc_immob + mycorrhizal_mine_N_immob
     total_mining_CO2prod = total_mining_CO2prod + mining_CO2prod
     total_myc_mine_C_uptake = total_myc_mine_C_uptake + myc_mine_C_uptake

     total_N_fixation = total_N_fixation + N_fixation
     total_N_fixer_C_allocated = total_N_fixer_C_allocated + N_fixer_C_allocated

   enddo

   soil%gross_nitrogen_flux_into_tile = soil%gross_nitrogen_flux_into_tile + total_N_fixation


   ! fsc_in and ssc_in updated in add_root_exudates
   call add_root_exudates(soil,vegn,total_root_exudate_C,total_root_exudate_N,total_myc_Nmin,N_leakage)

 if(is_watch_point()) then
     write(*,*),'Total soil C before Dsdt'
      __DEBUG1__(soil_tile_carbon(soil))
 endif

  ! update soil carbon
  call Dsdt(vegn, soil, diag, soilt, theta)

 if(is_watch_point()) then
     write(*,*),'Total soil C after Dsdt'
     __DEBUG2__(soil_tile_carbon(soil),vegn%rh*dt_fast_yr)
 endif


  ! Add respiration/C waste from mycorrhizae and N fixers
  vegn%rh = vegn%rh + (total_myc_resp + total_mining_CO2prod )/dt_fast_yr

  ! NEP is equal to NPP minus soil respiration
  vegn%nep = vegn%npp - vegn%rh
  if(is_watch_point()) then
     __DEBUG3__(vegn%npp,vegn%rh,vegn%nep)
     __DEBUG5__(total_root_exudate_C,total_scavenger_myc_C_allocated,total_miner_myc_C_allocated,total_N_fixer_C_allocated,total_mining_CO2prod)
  endif

 if(is_watch_point()) then
   write(*,*) 'Before update_soil_pools'
   __DEBUG4__(soil_tile_carbon(soil),vegn_tile_carbon(vegn),vegn%rh*dt_fast_yr,vegn%nep*dt_fast_yr)
 endif


  call update_soil_pools(vegn, soil)
  vegn%age = vegn%age + dt_fast_yr;

 if(is_watch_point()) then
   write(*,*) 'At end of vegn_carbon_int'
   __DEBUG4__(soil_tile_carbon(soil),vegn_tile_carbon(vegn),vegn%rh*dt_fast_yr,vegn%nep*dt_fast_yr)
 endif


  ! ---- diagnostic section
  call send_tile_data(id_gpp,gpp,diag)
  call send_tile_data(id_npp,vegn%npp,diag)
  call send_tile_data(id_nep,vegn%nep,diag)
  call send_tile_data(id_litter,vegn%litter,diag)
  call send_tile_data(id_resp, resp, diag)
  call send_tile_data(id_resl, resl, diag)
  call send_tile_data(id_resr, resr, diag)
  call send_tile_data(id_resg, resg, diag)
  call send_tile_data(id_soilt,soilt,diag)
  call send_tile_data(id_theta,theta,diag)
  call send_tile_data(id_mycorrhizal_scav_allocation,total_scavenger_myc_C_allocated/dt_fast_yr,diag)
  call send_tile_data(id_mycorrhizal_scav_N_res,vegn%cohorts(1)%scav_myc_N_reservoir,diag)
  call send_tile_data(id_mycorrhizal_scav_C_res,vegn%cohorts(1)%scav_myc_C_reservoir,diag)
  call send_tile_data(id_mycorrhizal_mine_N_res,vegn%cohorts(1)%mine_myc_N_reservoir,diag)
  call send_tile_data(id_mycorrhizal_mine_C_res,vegn%cohorts(1)%mine_myc_C_reservoir,diag)
  call send_tile_data(id_Nfix_N_res,vegn%cohorts(1)%N_fixer_N_reservoir,diag)
  call send_tile_data(id_Nfix_C_res,vegn%cohorts(1)%N_fixer_C_reservoir,diag)
  call send_tile_data(id_mycorrhizal_mine_allocation,total_miner_myc_C_allocated/dt_fast_yr,diag)
  call send_tile_data(id_mycorrhizal_mine_immobilization,total_mine_myc_immob/dt_fast_yr,diag)
  call send_tile_data(id_N_fixer_allocation,total_N_fixer_C_allocated/dt_fast_yr,diag)
  call send_tile_data(id_myc_scav_marginal_gain,vegn%cohorts(1)%myc_scav_marginal_gain_smoothed,diag)
  call send_tile_data(id_myc_mine_marginal_gain,vegn%cohorts(1)%myc_mine_marginal_gain_smoothed,diag)
  call send_tile_data(id_N_fix_marginal_gain,vegn%cohorts(1)%N_fix_marginal_gain_smoothed,diag)
  call send_tile_data(id_rhiz_exud_marginal_gain,vegn%cohorts(1)%rhiz_exud_marginal_gain_smoothed,diag)
  call send_tile_data(id_rhiz_exudation,total_root_exudate_C/dt_fast_yr,diag)
  call send_tile_data(id_nitrogen_stress,vegn%cohorts(1)%nitrogen_stress,diag)
  call send_tile_data(id_total_plant_N_uptake,total_plant_N_uptake/dt_fast_yr,diag)
  call send_tile_data(id_myc_scavenger_N_uptake,myc_scav_N_uptake/dt_fast_yr,diag)
  call send_tile_data(id_myc_miner_N_uptake,myc_mine_N_uptake/dt_fast_yr,diag)
  call send_tile_data(id_symbiotic_N_fixation,total_N_fixation/dt_fast_yr,diag)
  call send_tile_data(id_active_root_N_uptake, root_active_N_uptake/dt_fast_yr,diag)
  call send_tile_data(id_scav_plant_N_uptake,scav_N_to_plant/dt_fast_yr,diag)
  call send_tile_data(id_mine_plant_N_uptake,mine_N_to_plant/dt_fast_yr,diag)
  call send_tile_data(id_fix_plant_N_uptake,fix_N_to_plant/dt_fast_yr,diag)

  ! ---- CMOR diagnostics
  call send_tile_data(id_gpp_cmor, gpp/seconds_per_year, diag)
  call send_tile_data(id_npp_cmor, vegn%npp/seconds_per_year, diag)
  call send_tile_data(id_ra, (resp-resg)/seconds_per_year, diag)
  call send_tile_data(id_rgrowth, resg/seconds_per_year, diag)

end subroutine vegn_carbon_int


! ============================================================================
! updates cohort biomass pools, LAI, SAI, and height using accumulated
! carbon_gain and bwood_gain
subroutine vegn_growth (vegn)
  type(vegn_tile_type), intent(inout) :: vegn

  ! ---- local vars
  type(vegn_cohort_type), pointer :: cc    ! current cohort
  integer :: i

  do i = 1, vegn%n_cohorts
     cc => vegn%cohorts(i)

     cc%bwood   = cc%bwood   + cc%bwood_gain
     cc%bliving = cc%bliving + cc%carbon_gain

     if(cc%bliving < 0) then
        cc%bwood    = cc%bwood+cc%bliving
        cc%bliving  = 0
        if (cc%bwood < 0) &
             cc%bwood = 0 ! in principle, that's not conserving carbon
     endif

     call update_biomass_pools(cc)
     cc%root_density = (cc%br + &
            (cc%bsw+cc%bwood+cc%blv)*(1-agf_bs))*C2B
     cc%Wl_max = spdata(cc%species)%cmc_lai*cc%lai
     cc%Ws_max = spdata(cc%species)%csc_lai*cc%lai

     ! reset carbon accumulation terms
     cc%carbon_gain = 0
     cc%carbon_loss = 0
     cc%bwood_gain  = 0
     if (cc%status == LEAF_ON) then
        cc%leaf_age = cc%leaf_age + 1.0
        ! limit the maximum leaf age by the leaf time span (reciprocal of leaf
        ! turnover rate alpha) for given species. alpha is in 1/year, factor of
        ! 365 converts the result to days.
        if (spdata(cc%species)%alpha(CMPT_LEAF) > 0) &
             cc%leaf_age = min(cc%leaf_age,365.0/spdata(cc%species)%alpha(CMPT_LEAF))
     endif
  end do

end subroutine vegn_growth


! ============================================================================
subroutine eddy_npp(cc, tsoil)
  type(vegn_cohort_type), intent(inout) :: cc
  real, intent(in) :: tsoil

  call plant_respiration(cc,tsoil);

  cc%gpp = (cc%An_op - cc%An_cl)*mol_C*cc%lai;
  cc%npp = cc%gpp - cc%resp;

!  if(cc%npp_previous_day > -0.00001/2500.0) then
  if(cc%npp_previous_day > 0) then
     cc%resg = GROWTH_RESP*cc%npp_previous_day;
     cc%npp  = cc%npp  - GROWTH_RESP*cc%npp_previous_day;
     cc%resp = cc%resp + GROWTH_RESP*cc%npp_previous_day;
  else
     cc%resg = 0;
  endif

  ! update, accumulate
  cc%npp_previous_day_tmp = cc%npp_previous_day_tmp + cc%npp;
end subroutine eddy_npp


! ============================================================================
subroutine plant_respiration(cc, tsoil)
  type(vegn_cohort_type), intent(inout) :: cc
  real, intent(in) :: tsoil

  real :: tf,tfs;
  real :: r_leaf, r_vleaf, r_stem, r_root

  integer :: sp ! shorthand for cohort species
  sp = cc%species

  tf = exp(3000.0*(1.0/288.16-1.0/cc%Tv));
  tf = tf / ( &
            (1.0+exp(0.4*(5.0-cc%Tv+273.16)))*&
            (1.0+exp(0.4*(cc%Tv - 273.16-45.0)))&
            )

  tfs = exp(3000.0*(1.0/288.16-1.0/tsoil));
  tfs = tfs / ( &
              (1.0+exp(0.4*(5.0-tsoil+273.16)))* &
              (1.0+exp(0.4*(tsoil - 273.16-45.0)))&
              )

  r_leaf = -mol_C*cc%An_cl*cc%lai;
  r_vleaf = spdata(sp)%beta(CMPT_VLEAF)   * cc%blv*tf;
  r_stem  = spdata(sp)%beta(CMPT_SAPWOOD) * cc%bsw*tf;
  r_root  = spdata(sp)%beta(CMPT_ROOT)    * cc%br*tfs;

  cc%resp = r_leaf + r_vleaf + r_stem + r_root;
  cc%resl = r_leaf;
  cc%resr = r_root;
end subroutine plant_respiration


! ============================================================================
! calculates prev. day average NPP from accumualted values
subroutine vegn_daily_npp(vegn)
  type(vegn_tile_type), intent(inout) :: vegn

  integer :: n_fast_step;
  integer :: i
  type(vegn_cohort_type), pointer :: cc

  n_fast_step = 1.0/365.0/dt_fast_yr;
  do i = 1, vegn%n_cohorts
     cc => vegn%cohorts(i)
     vegn%cohorts(i)%npp_previous_day=vegn%cohorts(i)%npp_previous_day_tmp/n_fast_step;
     vegn%cohorts(i)%npp_previous_day_tmp=0.0
  enddo
end subroutine vegn_daily_npp


! =============================================================================
subroutine vegn_phenology(vegn, soil)
  type(vegn_tile_type), intent(inout) :: vegn
  type(soil_tile_type), intent(inout) :: soil
  ! TODO: possibly move soil-related calculations from calling procedure here,
  !       now that we have soil passed as an argument

  ! ---- local vars
  type(vegn_cohort_type), pointer :: cc
  real :: leaf_litter,root_litter;
  real :: theta_crit; ! critical ratio of average soil water to sat. water
  real :: psi_stress_crit ! critical soil-water-stress index
  real :: wilt ! ratio of wilting to saturated water content
  integer :: i, sp

  wilt = soil%w_wilt(1)/soil%pars%vwc_sat
  vegn%litter = 0

  do i = 1,vegn%n_cohorts
     cc => vegn%cohorts(i)
     sp=cc%species

     if(is_watch_point())then
        write(*,*)'####### vegn_phenology #######'
        __DEBUG4__(vegn%theta_av_phen, wilt, spdata(cc%species)%cnst_crit_phen, spdata(cc%species)%fact_crit_phen)
	__DEBUG2__(vegn%psist_av, spdata(cc%species)%psi_stress_crit_phen)
        __DEBUG1__(cc%species)
        __DEBUG2__(vegn%tc_av,spdata(cc%species)%tc_crit)
     endif
     ! if drought-deciduous or cold-deciduous species
     ! temp=10 degrees C gives good growing season pattern
     ! spp=0 is c3 grass,1 c3 grass,2 deciduous, 3 evergreen
     ! assumption is that no difference between drought and cold deciduous
     cc%status = LEAF_ON; ! set status to indicate no leaf drop

     if(cc%species < 4 )then! deciduous species
        ! actually either fact_crit_phen or cnst_crit_phen is zero, enforced
        ! by logic in the vegn_data.F90
        theta_crit = spdata(cc%species)%cnst_crit_phen &
              + wilt*spdata(cc%species)%fact_crit_phen
        theta_crit = max(0.0,min(1.0, theta_crit))
	psi_stress_crit = spdata(cc%species)%psi_stress_crit_phen
        if (      (psi_stress_crit <= 0. .and. vegn%theta_av_phen < theta_crit) &
	     .or. (psi_stress_crit  > 0. .and. vegn%psist_av > psi_stress_crit) &
             .or. (vegn%tc_av < spdata(cc%species)%tc_crit) ) then
           cc%status = LEAF_OFF; ! set status to indicate leaf drop
           cc%leaf_age = 0;

           leaf_litter = (1.0-l_fract)*cc%bl;
           root_litter = (1.0-l_fract)*cc%br;

           ! add to patch litter flux terms
           vegn%litter = vegn%litter + leaf_litter + root_litter;
           select case (soil_carbon_option)
           case(SOILC_CENTURY, SOILC_CENTURY_BY_LAYER)
              soil%fast_soil_C(1) = soil%fast_soil_C(1) +    fsc_liv *(leaf_litter+root_litter);
              soil%slow_soil_C(1) = soil%slow_soil_C(1) + (1-fsc_liv)*(leaf_litter+root_litter);
              ! soil%fsc_in(1)+=data->fsc_liv*(leaf_litter+root_litter);
              ! soil%ssc_in(1)+=(1.0-data->fsc_liv)*(leaf_litter+root_litter);
              soil%fsc_in(1)  = soil%fsc_in(1)  + leaf_litter+root_litter;
          case(SOILC_CORPSE_N)
            ! Retranslocation: l_fract determines how much leaf and root C is converted to blv !!!
            ! Retranslocation is converting the a potentially different fraction of N back to storage !
              call add_litter(soil%leafLitter,(/spdata(sp)%fsc_liv*leaf_litter,(1-spdata(sp)%fsc_liv)*leaf_litter,0.0/),&
                                                (/spdata(sp)%fsc_liv*cc%leaf_N*(1.0-spdata(sp)%leaf_retranslocation_frac),&
                                                (1-spdata(sp)%fsc_liv)*cc%leaf_N*(1.0-spdata(sp)%leaf_retranslocation_frac),0.0/))
              soil%leaflitter_fsc_in=soil%leaflitter_fsc_in+spdata(sp)%fsc_liv*leaf_litter
              soil%leaflitter_ssc_in=soil%leaflitter_ssc_in+(1-spdata(sp)%fsc_liv)*leaf_litter
              soil%leaflitter_fsn_in=soil%leaflitter_fsn_in+spdata(sp)%fsc_liv*cc%leaf_N*(1.0-spdata(sp)%leaf_retranslocation_frac)
              soil%leaflitter_ssn_in=soil%leaflitter_ssn_in+(1-spdata(sp)%fsc_liv)*cc%leaf_N*(1.0-spdata(sp)%leaf_retranslocation_frac)
              !ssc_in and fsc_in updated in add_root_litter
              call add_root_litter(soil, vegn, (/spdata(sp)%fsc_froot*root_litter,(1-spdata(sp)%fsc_froot)*root_litter,0.0/),&
                                (/spdata(sp)%fsc_froot*cc%root_N*(1.0-spdata(sp)%froot_retranslocation_frac),&
                                (1-spdata(sp)%fsc_froot)*cc%root_N*(1.0-spdata(sp)%froot_retranslocation_frac),0.0/))
            case(SOILC_CORPSE)
                call add_litter(soil%leafLitter,(/spdata(sp)%fsc_liv*leaf_litter,(1-spdata(sp)%fsc_liv)*leaf_litter,0.0/),&
                            (/0.0,0.0,0.0/))
                soil%leaflitter_fsc_in=soil%leaflitter_fsc_in+spdata(sp)%fsc_liv*leaf_litter
                soil%leaflitter_ssc_in=soil%leaflitter_ssc_in+(1-spdata(sp)%fsc_liv)*leaf_litter
                !ssc_in and fsc_in updated in add_root_litter
                call add_root_litter(soil, vegn, (/spdata(sp)%fsc_froot*root_litter,(1-spdata(sp)%fsc_froot)*root_litter,0.0/),&
                            (/0.0,0.0,0.0/))
           case default
              call error_mesg('vegn_phenology','The value of soil_carbon_option is invalid. This should never happen. Contact developer.',FATAL)
           end select
           vegn%veg_out = vegn%veg_out + leaf_litter+root_litter;

           cc%blv = cc%blv + l_fract*(cc%bl+cc%br);
           cc%bl  = 0.0;
           cc%br  = 0.0;
           cc%lai = 0.0;

           cc%stored_N = cc%stored_N + (cc%leaf_N*spdata(sp)%leaf_retranslocation_frac+cc%root_N*spdata(sp)%froot_retranslocation_frac)
           cc%leaf_N = 0.0
           cc%root_N = 0.0


           ! update state
           cc%bliving = cc%blv + cc%br + cc%bl + cc%bsw;
           cc%b = cc%bliving + cc%bwood ;
        endif
     endif
     call update_biomass_pools(cc);
  enddo
end subroutine vegn_phenology


! =============================================================================
subroutine vegn_biogeography(vegn)
  type(vegn_tile_type), intent(inout) :: vegn

  ! ---- local vars
  integer :: i,spp
  real :: excess_stored_N,target_stored_N

  do i = 1, vegn%n_cohorts
    associate (cc=>vegn%cohorts(i))
    spp=cc%species
     call update_species(cc, vegn%t_ann, vegn%t_cold, &
          vegn%p_ann*seconds_per_year, vegn%ncm, vegn%landuse)
    if(spp .ne. cc%species .and. soil_carbon_option .eq. SOILC_CORPSE_N) then
      ! Reset stored nitrogen when species changes
      call update_biomass_pools(cc)
      target_stored_N = 1.5*cc%bliving*(cc%Pl/spdata(cc%species)%leaf_live_c2n+cc%Pr/spdata(cc%species)%froot_live_c2n)
      if(cc%stored_N > target_stored_N) then
        excess_stored_N=cc%stored_N - target_stored_N
        cc%stored_N = cc%stored_N - excess_stored_N
        ! Deposit the excess N as root litter for now, using params from previous species (i.e. assuming some tissue died during transition)
        vegn%fsn_pool_bg = vegn%fsn_pool_bg + excess_stored_N*spdata(spp)%fsc_froot
        vegn%ssn_pool_bg = vegn%ssn_pool_bg + excess_stored_N*(1-spdata(spp)%fsc_froot)
      endif
    endif
    end associate
  enddo
end subroutine

! =============================================================================
! The stuff below comes from she_update.c -- it looks like it belongs here,
! since it is essentially a part of the carbon integration (update_patch_fast
! is only called immediately after carbon_int in lm3v)
! =============================================================================


! =============================================================================
subroutine update_soil_pools(vegn, soil)
  type(vegn_tile_type), intent(inout) :: vegn
  type(soil_tile_type), intent(inout) :: soil

  ! ---- local vars
  real :: delta;
  real :: deltafast, deltaslow, deltafast_N, deltaslow_N;

  !!!!!!!!!!!!!!!x2z /ADD deposition CHECK
  ! This should really be read in from some external data set or something

!  real ::  input_time_fert=130.0


  select case (soil_carbon_option)
  case (SOILC_CENTURY,SOILC_CENTURY_BY_LAYER)
     ! update fsc input rate so that intermediate fsc pool is never
     ! depleted below zero; on the other hand the pool can be only
     ! depleted, never increased
     vegn%fsc_rate_bg = MAX( 0.0, MIN(vegn%fsc_rate_bg, vegn%fsc_pool_bg/dt_fast_yr));
     delta = vegn%fsc_rate_bg*dt_fast_yr;
     soil%fast_soil_C(1) = soil%fast_soil_C(1) + delta;
     vegn%fsc_pool_bg    = vegn%fsc_pool_bg    - delta;

     ! update ssc input rate so that intermediate ssc pool is never
     ! depleted below zero; on the other hand the pool can be only
     ! depleted, never increased
     vegn%ssc_rate_bg = MAX(0.0, MIN(vegn%ssc_rate_bg, vegn%ssc_pool_bg/dt_fast_yr));
     delta = vegn%ssc_rate_bg*dt_fast_yr;
     soil%slow_soil_C(1) = soil%slow_soil_C(1) + delta;
     vegn%ssc_pool_bg    = vegn%ssc_pool_bg    - delta;

  case (SOILC_CORPSE, SOILC_CORPSE_N)


     ! NOTE that this code looks very weird from the point of view of carbon
     ! balance: we are depleting the {fsc,ssc}_pool_ag, but adding to litter
     ! from the pools {leaflitter,coarsewoodLitter}_buffer_ag. The latter two
     ! are not used in the calculations of the total carbon.
     ! -- BNS -- Good catch, I have tried to fix it.
     vegn%leaflitter_buffer_rate_fast = MAX(0.0, MIN(vegn%leaflitter_buffer_rate_fast, vegn%leaflitter_buffer_fast/dt_fast_yr))
     vegn%leaflitter_buffer_rate_slow = MAX(0.0, MIN(vegn%leaflitter_buffer_rate_slow, vegn%leaflitter_buffer_slow/dt_fast_yr))

     deltafast = vegn%leaflitter_buffer_rate_fast*dt_fast_yr
     deltaslow = vegn%leaflitter_buffer_rate_slow*dt_fast_yr

     if(soil_carbon_option == SOILC_CORPSE_N) then
         vegn%leaflitter_buffer_rate_fast_N = MAX(0.0, MIN(vegn%leaflitter_buffer_rate_fast_N, vegn%leaflitter_buffer_fast_N/dt_fast_yr))
         vegn%leaflitter_buffer_rate_slow_N = MAX(0.0, MIN(vegn%leaflitter_buffer_rate_slow_N, vegn%leaflitter_buffer_slow_N/dt_fast_yr))
         deltafast_N = vegn%leaflitter_buffer_rate_fast_N*dt_fast_yr
         deltaslow_N = vegn%leaflitter_buffer_rate_slow_N*dt_fast_yr
     else
         vegn%leaflitter_buffer_rate_fast_N=0.0
         vegn%leaflitter_buffer_rate_slow_N=0.0
         deltafast_N=0.0
         deltaslow_N=0.0
     endif

     call add_litter(soil%leafLitter,(/deltafast,deltaslow,0.0/),(/deltafast_N,deltaslow_N,0.0/))

     vegn%leaflitter_buffer_fast = vegn%leaflitter_buffer_fast - deltafast
     vegn%leaflitter_buffer_slow = vegn%leaflitter_buffer_slow - deltaslow
     vegn%leaflitter_buffer_fast_N = vegn%leaflitter_buffer_fast_N - deltafast_N
     vegn%leaflitter_buffer_slow_N = vegn%leaflitter_buffer_slow_N - deltaslow_N



     vegn%coarsewoodlitter_buffer_rate_fast = MAX(0.0, MIN(vegn%coarsewoodlitter_buffer_rate_fast, vegn%coarsewoodlitter_buffer_fast/dt_fast_yr))
     vegn%coarsewoodlitter_buffer_rate_slow = MAX(0.0, MIN(vegn%coarsewoodlitter_buffer_rate_slow, vegn%coarsewoodlitter_buffer_slow/dt_fast_yr))


     deltafast = vegn%coarsewoodlitter_buffer_rate_fast*dt_fast_yr
     deltaslow = vegn%coarsewoodlitter_buffer_rate_slow*dt_fast_yr

     if (soil_carbon_option == SOILC_CORPSE_N) then
         vegn%coarsewoodlitter_buffer_rate_fast_N = MAX(0.0, MIN(vegn%coarsewoodlitter_buffer_rate_fast_N, vegn%coarsewoodlitter_buffer_fast_N/dt_fast_yr))
         vegn%coarsewoodlitter_buffer_rate_slow_N = MAX(0.0, MIN(vegn%coarsewoodlitter_buffer_rate_slow_N, vegn%coarsewoodlitter_buffer_slow_N/dt_fast_yr))

         deltafast_N = vegn%coarsewoodlitter_buffer_rate_fast_N*dt_fast_yr
         deltaslow_N = vegn%coarsewoodlitter_buffer_rate_slow_N*dt_fast_yr
    else
        vegn%coarsewoodlitter_buffer_rate_fast_N= 0.0
        vegn%coarsewoodlitter_buffer_rate_slow_N=0.0
        deltafast_N=0.0
        deltaslow_N=0.0
    endif

     call add_litter(soil%coarsewoodlitter,(/deltafast,deltaslow,0.0/),(/deltafast_N,deltaslow_N,0.0/))

     vegn%coarsewoodlitter_buffer_fast = vegn%coarsewoodlitter_buffer_fast - deltafast
     vegn%coarsewoodlitter_buffer_slow = vegn%coarsewoodlitter_buffer_slow - deltaslow
     vegn%coarsewoodlitter_buffer_fast_N = vegn%coarsewoodlitter_buffer_fast_N - deltafast_N
     vegn%coarsewoodlitter_buffer_slow_N = vegn%coarsewoodlitter_buffer_slow_N - deltaslow_N



     ! update fsc input rate so that intermediate fsc pool is never
     ! depleted below zero; on the other hand the pool can be only
     ! depleted, never increased
     vegn%fsc_rate_bg = MAX( 0.0, MIN(vegn%fsc_rate_bg, vegn%fsc_pool_bg/dt_fast_yr));
     deltafast = vegn%fsc_rate_bg*dt_fast_yr;
     vegn%fsc_pool_bg    = vegn%fsc_pool_bg    - deltafast;



     ! update ssc input rate so that intermediate ssc pool is never
     ! depleted below zero; on the other hand the pool can be only
     ! depleted, never increased
     vegn%ssc_rate_bg = MAX(0.0, MIN(vegn%ssc_rate_bg, vegn%ssc_pool_bg/dt_fast_yr));
     deltaslow = vegn%ssc_rate_bg*dt_fast_yr;
     vegn%ssc_pool_bg    = vegn%ssc_pool_bg    - deltaslow;

     if(soil_carbon_option == SOILC_CORPSE_N) then
         vegn%fsn_rate_bg = MAX( 0.0, MIN(vegn%fsn_rate_bg, vegn%fsn_pool_bg/dt_fast_yr));
         deltafast_N = vegn%fsn_rate_bg*dt_fast_yr;
         vegn%fsn_pool_bg    = vegn%fsn_pool_bg    - deltafast_N;

         vegn%ssn_rate_bg = MAX(0.0, MIN(vegn%ssn_rate_bg, vegn%ssn_pool_bg/dt_fast_yr));
         deltaslow_N = vegn%ssn_rate_bg*dt_fast_yr;
         vegn%ssn_pool_bg    = vegn%ssn_pool_bg    - deltaslow_N;
     else
         vegn%fsn_rate_bg = 0.0
         deltafast_N = 0.0
         vegn%fsn_pool_bg    = 0.0

         vegn%ssn_rate_bg = 0.0
         deltaslow_N = 0.0
         vegn%ssn_pool_bg    = 0.0
     endif


     call add_root_litter(soil, vegn, (/deltafast,deltaslow,0.0/), (/deltafast_N,deltaslow_N,0.0/))
  case default
     call error_mesg('update_soil_pools','The value of soil_carbon_option is invalid. This should never happen. Contact developer.',FATAL)
  end select

end subroutine update_soil_pools



end module vegn_dynamics_mod
