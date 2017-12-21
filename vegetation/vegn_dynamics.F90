! ============================================================================
! updates carbon pools and rates on the fast time scale
! ============================================================================
module vegn_dynamics_mod

#include "../shared/debug.inc"

use fms_mod, only: error_mesg, FATAL, WARNING
use time_manager_mod, only: time_type
use mpp_mod, only: mpp_sum, mpp_pe, mpp_root_pe
use mpp_domains_mod, only : mpp_pass_UG_to_SG, mpp_pass_SG_to_UG, mpp_update_domains

use constants_mod, only : PI,tfreeze
use land_constants_mod, only : days_per_year, seconds_per_year, mol_C
use land_data_mod, only : lnd,log_version
use land_debug_mod, only : is_watch_point, check_var_range, check_conservation, &
     do_check_conservation, carbon_cons_tol, nitrogen_cons_tol, set_current_point, &
     land_error_message
use land_tile_mod, only : land_tile_map, land_tile_type, land_tile_enum_type, &
     first_elmt, loop_over_tiles, land_tile_carbon, land_tile_nitrogen
use land_tile_diag_mod, only : OP_SUM, OP_AVERAGE, &
     register_tiled_diag_field, send_tile_data, diag_buff_type, &
     register_cohort_diag_field, send_cohort_data, set_default_diag_filter
use vegn_data_mod, only : spdata, nspecies, do_ppa, &
     PHEN_DECIDUOUS, PHEN_EVERGREEN, LEAF_ON, LEAF_OFF, FORM_WOODY, FORM_GRASS, &
     ALLOM_EW, ALLOM_EW1, ALLOM_HML, LU_CROP, &
     agf_bs, l_fract, understory_lai_factor, min_lai, &
     use_light_saber, laimax_ceiling, laimax_floor, &
     myc_scav_C_efficiency, myc_mine_C_efficiency, N_fixer_C_efficiency, N_limits_live_biomass, &
     excess_stored_N_leakage_rate, min_N_stress, &
     c2n_N_fixer, et_myc, smooth_N_uptake_C_allocation, N_fix_Tdep_Houlton, &
     mycorrhizal_turnover_time, N_fixer_turnover_time
use vegn_tile_mod, only: vegn_tile_type, vegn_tile_carbon, vegn_tile_nitrogen, &
     vegn_mergecohorts_ppa, vegn_relayer_cohorts_ppa
use soil_tile_mod, only: num_l, dz, soil_tile_type, N_LITTER_POOLS
use vegn_cohort_mod, only : vegn_cohort_type, &
     update_biomass_pools, update_bio_living_fraction, update_species, &
     leaf_area_from_biomass, cohort_root_litter_profile, cohort_root_exudate_profile, &
     plant_C, plant_N
use vegn_util_mod, only : kill_plants_ppa, add_seedlings_ppa
use vegn_harvesting_mod, only : allow_weeds_on_crops
use soil_carbon_mod, only: N_C_TYPES, soil_carbon_option, &
    SOILC_CENTURY, SOILC_CENTURY_BY_LAYER, SOILC_CORPSE, SOILC_CORPSE_N, &
    add_litter, debug_pool,deadmic_slow_frac
use soil_util_mod, only: add_soil_carbon, add_root_litter, add_root_exudates
use soil_mod, only: Dsdt, root_N_uptake, myc_scavenger_N_uptake, myc_miner_N_uptake

implicit none
private

! ==== public interfaces =====================================================
public :: vegn_dynamics_init, vegn_dynamics_end

public :: vegn_carbon_int_lm3   ! fast time-scale integrator of carbon balance
public :: vegn_carbon_int_ppa   ! fast time-scale integrator of carbon balance
public :: vegn_growth           ! slow time-scale redistributor of accumulated carbon
public :: vegn_phenology_lm3
public :: vegn_phenology_ppa
public :: vegn_biogeography
public :: vegn_reproduction_ppa

public :: vegn_starvation_ppa   !
! ==== end of public interfaces ==============================================

! ==== module constants ======================================================
character(len=*), parameter :: module_name = 'vegn_dynamics_mod'
#include "../shared/version_variable.inc"
character(len=*), parameter :: diag_mod_name = 'vegn'

! ==== module data ===========================================================
real    :: dt_fast_yr ! fast (physical) time step, yr (year is defined as 365 days)
real, allocatable :: sg_soilfrac(:,:) ! fraction of grid cell area occupied by soil, for
  ! normalization in seed transort
real, allocatable :: ug_area_factor(:) ! conversion factor from land area to vegetation area
real :: atot ! global land area for normalization in conservation checks, m2

! diagnostic field IDs
integer :: id_npp, id_nep, id_gpp, id_wood_prod, id_leaf_root_gr, id_sw_seed_gr
integer :: id_rsoil, id_rsoil_fast, id_rsoil_slow
integer :: id_resp, id_resl, id_resr, id_ress, id_resg, id_asoil
integer :: id_soilt, id_theta, id_litter, id_age, id_dbh_growth
integer :: &
    id_mycorrhizal_scav_allocation, id_mycorrhizal_scav_immobilization, &
    id_mycorrhizal_mine_allocation, id_mycorrhizal_mine_immobilization, &
    id_N_fixer_allocation, id_total_plant_N_uptake, &
    id_N_fix_marginal_gain, id_myc_scav_marginal_gain, &
    id_myc_mine_marginal_gain, id_rhiz_exudation, id_nitrogen_stress, &
    id_rhiz_exud_marginal_gain,id_myc_scavenger_N_uptake,&
    id_myc_miner_N_uptake,id_symbiotic_N_fixation,id_active_root_N_uptake,&
    id_scav_plant_N_uptake, id_mine_plant_N_uptake, id_fix_plant_N_uptake,&
    id_exudate, &
    id_mycorrhizal_scav_C_res, id_mycorrhizal_scav_N_res, &
    id_mycorrhizal_mine_C_res, id_mycorrhizal_mine_N_res, &
    id_Nfix_C_res, id_Nfix_N_res, &
    id_N_fix_alloc_smoothed, id_myc_mine_alloc_smoothed, id_myc_scav_alloc_smoothed

contains

! ============================================================================
subroutine vegn_dynamics_init(id_ug, time, delta_time)
  integer        , intent(in) :: id_ug   !<Unstructured axis id.
  type(time_type), intent(in) :: time       ! initial time for diagnostic fields
  real           , intent(in) :: delta_time ! fast time step, s

  type(land_tile_enum_type) :: ce
  type(land_tile_type), pointer :: tile
  integer :: l ! grid cell index on unstructured gris

  call log_version(version, module_name, &
  __FILE__)

  ! set up global variables
  dt_fast_yr = delta_time/seconds_per_year

  ! calculate fraction of the grid cell occupied by soil in ug_area_factor
  allocate(ug_area_factor(lnd%ls:lnd%le))
  ug_area_factor = 0.0
  ce = first_elmt(land_tile_map,lnd%ls)
  do while (loop_over_tiles(ce,tile,l))
     if(associated(tile%vegn)) &
          ug_area_factor(l) = ug_area_factor(l) + tile%frac
  enddo
  ! calculate soil fraction in a grid cell
  allocate (sg_soilfrac(lnd%isd:lnd%ied, lnd%jsd:lnd%jed))
  sg_soilfrac = 0.0
  call mpp_pass_UG_to_SG(lnd%ug_domain,ug_area_factor*lnd%ug_landfrac,sg_soilfrac)
  call mpp_update_domains(sg_soilfrac,lnd%sg_domain)
  ! finish calculating conversion factor from land area to soil area
  where(ug_area_factor>0) &
        ug_area_factor = 1.0/ug_area_factor

  ! calculate total land area for normalization in global conservation checks
  atot = sum(lnd%ug_area)
  call mpp_sum(atot)


  ! set the default sub-sampling filter for the fields below
  call set_default_diag_filter('soil')

  ! register diagnostic fields
  id_gpp = register_cohort_diag_field ( diag_mod_name, 'gpp',  &
       (/id_ug/), time, 'gross primary productivity', 'kg C/(m2 year)', &
       missing_value=-100.0)
  id_npp = register_cohort_diag_field ( diag_mod_name, 'npp',  &
       (/id_ug/), time, 'net primary productivity', 'kg C/(m2 year)', &
       missing_value=-100.0)
  id_nep = register_tiled_diag_field ( diag_mod_name, 'nep',  &
       (/id_ug/), time, 'net ecosystem productivity', 'kg C/(m2 year)', &
       missing_value=-100.0 )
  id_wood_prod = register_cohort_diag_field ( diag_mod_name, 'wood_prod',  &
       (/id_ug/), time, 'total wood (heartwood+sapwood) production', 'kgC/(m2 year)', &
       missing_value=-100.0)
  id_leaf_root_gr = register_cohort_diag_field ( diag_mod_name, 'leaf_root_gr',  &
       (/id_ug/), time, 'C allocated to leaves & fine roots', 'kgC/(m2 year)', &
       missing_value=-100.0)
  id_sw_seed_gr = register_cohort_diag_field ( diag_mod_name, 'sw_seed_gr',  &
       (/id_ug/), time, 'C respired in seeds and sapwood growth', 'kgC/(m2 year)', &
       missing_value=-100.0)
  id_dbh_growth = register_cohort_diag_field ( diag_mod_name, 'dbh_gr',  &
       (/id_ug/), time, 'growth rathe of DBH', 'm/year', &
       missing_value=-100.0)
  id_litter = register_tiled_diag_field (diag_mod_name, 'litter', (/id_ug/), &
       time, 'litter productivity', 'kg C/(m2 year)', missing_value=-100.0)
  id_resp = register_cohort_diag_field ( diag_mod_name, 'resp', (/id_ug/), &
       time, 'respiration', 'kg C/(m2 year)', missing_value=-100.0)
  id_resl = register_cohort_diag_field ( diag_mod_name, 'resl', (/id_ug/), &
       time, 'leaf respiration', 'kg C/(m2 year)', missing_value=-100.0)
  id_resr = register_cohort_diag_field ( diag_mod_name, 'resr', (/id_ug/), &
       time, 'root respiration', 'kg C/(m2 year)', missing_value=-100.0)
  id_ress = register_cohort_diag_field ( diag_mod_name, 'ress', (/id_ug/), &
       time, 'stem respiration', 'kg C/(m2 year)', missing_value=-100.0)
  id_resg = register_cohort_diag_field ( diag_mod_name, 'resg', (/id_ug/), &
       time, 'growth respiration', 'kg C/(m2 year)', missing_value=-100.0)
  id_soilt = register_tiled_diag_field ( diag_mod_name, 'tsoil_av',  &
       (/id_ug/), time, 'average soil temperature for carbon decomposition', 'degK', &
       missing_value=-100.0 )
  id_theta = register_tiled_diag_field ( diag_mod_name, 'theta',  &
       (/id_ug/), time, 'average soil wetness for carbon decomposition', 'm3/m3', &
       missing_value=-100.0 )
  id_age = register_cohort_diag_field ( diag_mod_name, 'age',  &
       (/id_ug/), time, 'average cohort age', 'years', &
       missing_value=-100.0)
  id_exudate = register_cohort_diag_field ( diag_mod_name, 'exudate', (/id_ug/), &
       time, 'carbon root exudates', 'kg C/(m2 year)', missing_value=-100.0)
! FIXME slm: perhaps the the following fields need to be cohort fields?
  id_mycorrhizal_scav_allocation = register_cohort_diag_field ( diag_mod_name, 'mycorrhizal_scav_allocation',  &
       (/id_ug/), time, 'C allocation to scavenger mycorrhizae', 'kg C/(m2 year)', &
       missing_value=-100.0 )
  id_mycorrhizal_scav_immobilization = register_cohort_diag_field ( diag_mod_name, 'mycorrhizal_scav_immobilization',  &
        (/id_ug/), time, 'N immobilization by scavenger mycorrhizae', 'kg N/(m2 year)', &
        missing_value=-100.0 )
  id_mycorrhizal_mine_allocation = register_cohort_diag_field ( diag_mod_name, 'mycorrhizal_mine_allocation',  &
       (/id_ug/), time, 'C allocation to miner mycorrhizae', 'kg C/(m2 year)', &
       missing_value=-100.0 )
  id_mycorrhizal_mine_immobilization = register_cohort_diag_field ( diag_mod_name, 'mycorrhizal_mine_immobilization',  &
       (/id_ug/), time, 'N immobilization by miner mycorrhizae', 'kg N/(m2 year)', &
       missing_value=-100.0 )
  id_N_fixer_allocation = register_cohort_diag_field ( diag_mod_name, 'N_fixer_allocation',  &
       (/id_ug/), time, 'C allocation to N fixers', 'kg C/(m2 year)', &
       missing_value=-100.0 )
  id_N_fix_marginal_gain = register_cohort_diag_field ( diag_mod_name, 'N_fix_marginal_gain',  &
       (/id_ug/), time, 'Extra N fixation per unit C allocation', 'kg N/(m2 year)/kgC', &
       missing_value=-100.0 )
  id_myc_scav_marginal_gain = register_cohort_diag_field ( diag_mod_name, 'myc_scav_marginal_gain',  &
       (/id_ug/), time, 'Extra N acquisition per unit C allocation to scavenger mycorrhizae', 'kg N/(m2 year)/kgC', &
       missing_value=-100.0 )
  id_myc_mine_marginal_gain = register_cohort_diag_field ( diag_mod_name, 'myc_mine_marginal_gain',  &
       (/id_ug/), time, 'Extra N acquisition per unit C allocation to miner mycorrhizae', 'kg N/(m2 year)/kg C', &
       missing_value=-100.0 )
  id_rhiz_exudation = register_cohort_diag_field ( diag_mod_name, 'rhiz_exudation',  &
       (/id_ug/), time, 'C allocation to rhizosphere exudation', 'kg C/(m2 year)', &
       missing_value=-100.0 )
  id_nitrogen_stress = register_cohort_diag_field ( diag_mod_name, 'nitrogen_stress',  &
       (/id_ug/), time, 'Nitrogen stress index', 'Dimensionless', &
       missing_value=-100.0 )
  id_total_plant_N_uptake = register_cohort_diag_field ( diag_mod_name, 'plant_N_uptake',  &
       (/id_ug/), time, 'Plant N uptake rate', 'kg N/(m2 year)', &
       missing_value=-100.0 )
  id_rhiz_exud_marginal_gain = register_cohort_diag_field ( diag_mod_name, 'rhiz_exud_marginal_gain',  &
       (/id_ug/), time, 'Extra N acquisition per unit rhiz C exudation', 'kg N/(m2 year)/kgC', &
       missing_value=-100.0 )
  id_scav_plant_N_uptake = register_cohort_diag_field ( diag_mod_name, 'plant_scavenger_N_uptake',  &
       (/id_ug/), time, 'Plant N uptake from scavenger mycorrhizae', 'kg N/m2/year', missing_value=-1.0 )
  id_mine_plant_N_uptake = register_cohort_diag_field ( diag_mod_name, 'plant_miner_N_uptake',  &
       (/id_ug/), time, 'Plant N uptake from miner mycorrhizae', 'kg N/m2/year', missing_value=-1.0 )
  id_fix_plant_N_uptake = register_cohort_diag_field ( diag_mod_name, 'plant_fixer_N_uptake',  &
        (/id_ug/), time, 'Plant N uptake from N fixers', 'kg N/m2/year', missing_value=-1.0 )

  id_myc_scavenger_N_uptake = register_cohort_diag_field ( diag_mod_name, 'myc_scavenger_N_uptake',  &
       (/id_ug/), time, 'N uptake by scavenger mycorrhizae', 'kg N/m2/year', missing_value=-1.0 )
  id_myc_miner_N_uptake = register_cohort_diag_field ( diag_mod_name, 'myc_miner_N_uptake',  &
       (/id_ug/), time, 'N uptake by miner mycorrhizae', 'kg N/m2/year', missing_value=-1.0 )
  id_symbiotic_N_fixation = register_cohort_diag_field ( diag_mod_name, 'symbiotic_N_fixation',  &
       (/id_ug/), time, 'Symbiotic N fixation', 'kg N/m2/year', missing_value=-1.0 )
  id_active_root_N_uptake = register_cohort_diag_field ( diag_mod_name, 'active_root_N_uptake',  &
       (/id_ug/), time, 'N uptake by root active transport', 'kg N/m2/year', missing_value=-1.0 )

  id_mycorrhizal_scav_C_res = register_tiled_diag_field ( module_name, 'myc_scavenger_C_res',  &
       (/id_ug/), time, 'Scavenger mycorrhizae C reservoir', 'kg C/m2', missing_value=-1.0 )
  id_mycorrhizal_scav_N_res = register_tiled_diag_field ( module_name, 'myc_scavenger_N_res',  &
       (/id_ug/), time, 'Scavenger mycorrhizae N reservoir', 'kg N/m2', missing_value=-1.0 )
  id_mycorrhizal_mine_C_res = register_tiled_diag_field ( module_name, 'myc_miner_C_res',  &
       (/id_ug/), time, 'Miner mycorrhizae C reservoir', 'kg C/m2', missing_value=-1.0 )
  id_mycorrhizal_mine_N_res = register_tiled_diag_field ( module_name, 'myc_miner_N_res',  &
       (/id_ug/), time, 'Miner mycorrhizae N reservoir', 'kg N/m2', missing_value=-1.0 )
  id_Nfix_C_res = register_tiled_diag_field ( module_name, 'N_fixer_C_res',  &
       (/id_ug/), time, 'N fixer C reservoir', 'kg C/m2', missing_value=-1.0 )
  id_Nfix_N_res = register_tiled_diag_field ( module_name, 'N_fixer_N_res',  &
       (/id_ug/), time, 'N fixer N reservoir', 'kg N/m2', missing_value=-1.0 )

  id_N_fix_alloc_smoothed = register_tiled_diag_field ( module_name, 'N_fix_alloc_smoothed',  &
       (/id_ug/), time, 'Plant C allocation to N fixers smoothed', 'kg N/m2/year', missing_value=-1.0 )
  id_myc_mine_alloc_smoothed = register_tiled_diag_field ( module_name, 'myc_mine_alloc_smoothed',  &
       (/id_ug/), time, 'Plant C allocation to N miners smoothed', 'kg N/m2/year', missing_value=-1.0 )
  id_myc_scav_alloc_smoothed = register_tiled_diag_field ( module_name, 'myc_scav_alloc_smoothed',  &
       (/id_ug/), time, 'C allocation to N scavengers smoothed', 'kg N/m2/year', missing_value=-1.0 )
end subroutine vegn_dynamics_init

! ============================================================================
subroutine vegn_dynamics_end()
   deallocate(sg_soilfrac, ug_area_factor)
end subroutine vegn_dynamics_end


subroutine  update_mycorrhizae(cc, soilT, &
  C_allocation_to_N_acq,myc_scav_N_uptake,myc_mine_N_uptake,myc_mine_C_uptake,root_N_uptake,&
  myc_scav_efficiency,myc_mine_efficiency,&
  scavenger_myc_C_allocated,miner_myc_C_allocated,N_fixer_C_allocated,&
  root_exudate_C, root_exudate_N,&
  scav_N_to_plant,mine_N_to_plant,fix_N_to_plant,&
  myc_CO2_prod,N_fixation,total_plant_N_uptake,&
  myc_turnover_C,myc_turnover_N,myc_Nmin)

  type(vegn_cohort_type), intent(inout) :: cc
  real, intent(in) :: soilT ! soil temperature, degK
  real, intent(in) :: myc_scav_N_uptake,myc_mine_N_uptake,myc_mine_C_uptake,root_N_uptake  ! N uptake per individual (kgN/m2)
  real, intent(in) :: C_allocation_to_N_acq
  real, intent(in) :: myc_scav_efficiency,myc_mine_efficiency         ! N uptake/myc biomass -- Only used if myc biomass is zero
  real, intent(out):: scavenger_myc_C_allocated,miner_myc_C_allocated,N_fixer_C_allocated ! Allocation per individual (kgC/m2)
  real, intent(out):: root_exudate_C, root_exudate_N
  real, intent(out):: scav_N_to_plant,mine_N_to_plant,fix_N_to_plant
  real, intent(out):: myc_CO2_prod, myc_Nmin  ! kgC/m2/individual
  real, intent(out):: N_fixation, total_plant_N_uptake
  real, intent(out):: myc_turnover_C,myc_turnover_N

  real :: myc_scav_marginal_gain,myc_mine_marginal_gain,N_fix_marginal_gain,rhiz_exud_marginal_gain ! kgN/kgC allocated
  real :: myc_scav_exudate_frac,myc_mine_exudate_frac,N_fixer_exudate_frac,rhiz_exud_frac
  real :: scavenger_myc_N_allocated,miner_myc_N_allocated,N_fixer_N_allocated
  real :: d_scav_C_reservoir,d_scav_N_reservoir,d_mine_C_reservoir,d_mine_N_reservoir,d_N_fixer_C_reservoir,d_N_fixer_N_reservoir
  real :: scavenger_myc_growth, miner_myc_growth, N_fixer_growth
  real :: lim_factor
  real :: reservoir_C_leakage, maint_resp
  real :: stored_N_loss ! loss of stored N due to export, kg N
  real :: w ! smoothing weight
  real :: mgain ! sum of marginal gains for marginal gain fraction calculations

  if(is_watch_point()) then
     write(*,*)'#### update_mycorrhizae input ####'
     __DEBUG4__(myc_scav_N_uptake,myc_mine_N_uptake,myc_mine_C_uptake,root_N_uptake)
     __DEBUG1__(C_allocation_to_N_acq)
     __DEBUG2__(myc_scav_efficiency,myc_mine_efficiency)
     write(*,*)'#### end of update_mycorrhizae input ####'
  endif

  associate(sp=>spdata(cc%species))
  myc_CO2_prod = 0.0
  myc_Nmin = 0.0
  reservoir_C_leakage = 0.0

  call check_var_range(C_allocation_to_N_acq,0.0,HUGE(0.0),'update_mycorrhizae','C_allocation_to_N_acq',FATAL)
  if(N_limits_live_biomass) &
       call check_var_range(cc%stored_N,0.0,HUGE(1.0),'update_mycorrhizae input','cc%stored_N',FATAL)

  if (soil_carbon_option == SOILC_CORPSE_N) then
     ! Update reservoirs with N uptake for this cohort
     cc%scav_myc_N_reservoir = cc%scav_myc_N_reservoir+myc_scav_N_uptake
     cc%mine_myc_C_reservoir = cc%mine_myc_C_reservoir+myc_mine_C_uptake
     cc%mine_myc_N_reservoir = cc%mine_myc_N_reservoir+myc_mine_N_uptake
     N_fixation = cc%N_fixer_biomass_C*sp%N_fixation_rate*dt_fast_yr

     ! T dependence of N fixation from Houlton et al. (2008) Nature paper (normalized to peak at 1.0)
     ! a=-3.62, b=0.27, c=25.15, T effect = exp(-0.5*b*c+b*Ts*(1-0.5*Ts/c))
     if(N_fix_Tdep_Houlton) N_fixation=N_fixation*exp(-0.5*0.27*25.15 + 0.27*(soilT-273.15)*(1.0-0.5*(soilT-273.15)/25.15))

     cc%N_fixer_N_reservoir = cc%N_fixer_N_reservoir + N_fixation

     call check_var_range(cc%myc_scavenger_biomass_C, 0.0,HUGE(0.0),'update_mycorrhizae','cc%myc_scavenger_biomass_C',FATAL)
     call check_var_range(cc%myc_miner_biomass_C,     0.0,HUGE(0.0),'update_mycorrhizae','cc%myc_miner_biomass_C',FATAL)
     call check_var_range(cc%N_fixer_biomass_C,       0.0,HUGE(0.0),'update_mycorrhizae','cc%N_fixer_biomass_C',FATAL)
     if (is_watch_point()) then
        __DEBUG3__(cc%myc_scavenger_biomass_C,cc%myc_miner_biomass_C,cc%N_fixer_biomass_C)
     endif

     ! Add mycorrhizal and N fixer turnover to soil C pools
     myc_turnover_C = ((cc%myc_scavenger_biomass_C+&
          cc%myc_miner_biomass_C)/mycorrhizal_turnover_time+&
          (cc%N_fixer_biomass_C)/N_fixer_turnover_time)*dt_fast_yr*et_myc
     myc_turnover_N = ((cc%myc_scavenger_biomass_N+&
          cc%myc_miner_biomass_N)/mycorrhizal_turnover_time+&
          (cc%N_fixer_biomass_N)/N_fixer_turnover_time)*dt_fast_yr*et_myc

     myc_CO2_prod = myc_CO2_prod + (1.0-et_myc)*myc_turnover_C/et_myc
     myc_Nmin = myc_Nmin + (1.0-et_myc)*myc_turnover_N/et_myc

     d_scav_C_reservoir = 0.0
     d_scav_N_reservoir = 0.0
     d_mine_C_reservoir = 0.0
     d_mine_N_reservoir = 0.0
     d_N_fixer_C_reservoir = 0.0
     d_N_fixer_N_reservoir = 0.0

     ! Then update biomass of scavenger mycorrhizae and N fixers
     scavenger_myc_growth = sp%myc_growth_rate*cc%scav_myc_C_reservoir/(cc%scav_myc_C_reservoir+sp%kM_myc_growth)*myc_scav_C_efficiency*dt_fast_yr
     maint_resp=min(cc%myc_scavenger_biomass_C/mycorrhizal_turnover_time*(1.0-et_myc)*dt_fast_yr,scavenger_myc_growth)
     if (scavenger_myc_growth-maint_resp>sp%c2n_mycorrhizae*cc%scav_myc_N_reservoir*0.9) then
        ! Not enough nitrogen to support growth. Limit to available N, and leave a little bit left over for plant
        scavenger_myc_growth=sp%c2n_mycorrhizae*cc%scav_myc_N_reservoir*0.9+maint_resp
     endif
     myc_CO2_prod = myc_CO2_prod + scavenger_myc_growth/myc_scav_C_efficiency*(1.0-myc_scav_C_efficiency)

     cc%scav_myc_N_reservoir=cc%scav_myc_N_reservoir+cc%myc_scavenger_biomass_N*(1-et_myc/mycorrhizal_turnover_time*dt_fast_yr)
     cc%myc_scavenger_biomass_C = cc%myc_scavenger_biomass_C + scavenger_myc_growth - cc%myc_scavenger_biomass_C/mycorrhizal_turnover_time*dt_fast_yr
     cc%myc_scavenger_biomass_N = cc%myc_scavenger_biomass_C/sp%c2n_mycorrhizae
     d_scav_C_reservoir = d_scav_C_reservoir - scavenger_myc_growth/myc_scav_C_efficiency
     cc%scav_myc_N_reservoir=cc%scav_myc_N_reservoir-cc%myc_scavenger_biomass_N


     miner_myc_growth = sp%myc_growth_rate*cc%mine_myc_C_reservoir/(cc%mine_myc_C_reservoir+sp%kM_myc_growth)*myc_mine_C_efficiency*dt_fast_yr
     maint_resp=min(cc%myc_miner_biomass_C/mycorrhizal_turnover_time*(1.0-et_myc)*dt_fast_yr,miner_myc_growth)
     if (miner_myc_growth-maint_resp>sp%c2n_mycorrhizae*cc%mine_myc_N_reservoir*0.9) then
       ! Not enough nitrogen to support growth. Limit to available N, and leave a little bit left over for plant
       miner_myc_growth=sp%c2n_mycorrhizae*cc%mine_myc_N_reservoir*0.9+maint_resp
     endif
     myc_CO2_prod = myc_CO2_prod + miner_myc_growth/myc_mine_C_efficiency*(1.0-myc_mine_C_efficiency)

     cc%mine_myc_N_reservoir=cc%mine_myc_N_reservoir+cc%myc_miner_biomass_N*(1-et_myc/mycorrhizal_turnover_time*dt_fast_yr)
     cc%myc_miner_biomass_C = cc%myc_miner_biomass_C + miner_myc_growth - cc%myc_miner_biomass_C/mycorrhizal_turnover_time*dt_fast_yr
     cc%myc_miner_biomass_N = cc%myc_miner_biomass_C/sp%c2n_mycorrhizae
     d_mine_C_reservoir = d_mine_C_reservoir - miner_myc_growth/myc_mine_C_efficiency
     cc%mine_myc_N_reservoir=cc%mine_myc_N_reservoir-cc%myc_miner_biomass_N

     N_fixer_growth = sp%myc_growth_rate*cc%N_fixer_C_reservoir/(cc%N_fixer_C_reservoir+sp%kM_myc_growth)*N_fixer_C_efficiency*dt_fast_yr
     maint_resp=min(cc%N_fixer_biomass_C/N_fixer_turnover_time*(1.0-et_myc)*dt_fast_yr,N_fixer_growth)
     ! if (N_fixer_growth>c2n_mycorrhizae*cc%N_fixer_N_reservoir*0.9) then
     !   ! Not enough nitrogen to support growth. Limit to available N, and leave a little bit left over for plant
     !   N_fixer_growth=c2n_mycorrhizae*cc%N_fixer_N_reservoir*0.9
     ! endif

     myc_CO2_prod = myc_CO2_prod + N_fixer_growth/N_fixer_C_efficiency*(1.0-N_fixer_C_efficiency)

     N_fixation=N_fixation-cc%N_fixer_biomass_N*(1-et_myc/N_fixer_turnover_time*dt_fast_yr)
     cc%N_fixer_biomass_C = cc%N_fixer_biomass_C + N_fixer_growth - cc%N_fixer_biomass_C/N_fixer_turnover_time*dt_fast_yr
     cc%N_fixer_biomass_N = cc%N_fixer_biomass_N/c2n_N_fixer
     d_N_fixer_C_reservoir = d_N_fixer_C_reservoir - N_fixer_growth/N_fixer_C_efficiency
     ! d_N_fixer_N_reservoir = d_N_fixer_N_reservoir - N_fixer_growth/c2n_mycorrhizae

     ! N fixers just make all the N they need for their biomass
     N_fixation=N_fixation+cc%N_fixer_biomass_N


     ! Prevent numerical errors from very small numbers if biomass is decreasing exponentially by killing it all at some point
     if(cc%myc_scavenger_biomass_C < 1e-20) then
        myc_turnover_C = myc_turnover_C + cc%myc_scavenger_biomass_C
        myc_turnover_N = myc_turnover_N + cc%myc_scavenger_biomass_N
        cc%myc_scavenger_biomass_C = 0.0
        cc%myc_scavenger_biomass_N = 0.0
     endif

     if(cc%myc_miner_biomass_C < 1e-20) then
        myc_turnover_C = myc_turnover_C + cc%myc_miner_biomass_C
        myc_turnover_N = myc_turnover_N + cc%myc_miner_biomass_N
        cc%myc_miner_biomass_C = 0.0
        cc%myc_miner_biomass_N = 0.0
     endif

     if(cc%N_fixer_biomass_C < 1e-20) then
        myc_turnover_C = myc_turnover_C + cc%N_fixer_biomass_C
        myc_turnover_N = myc_turnover_N + cc%N_fixer_biomass_N
        cc%N_fixer_biomass_C = 0.0
        cc%N_fixer_biomass_N = 0.0
     endif

     call check_var_range(cc%myc_scavenger_biomass_C, 0.0,HUGE(0.0),'update_mycorrhizae #2','cc%myc_scavenger_biomass_C',FATAL)
     call check_var_range(cc%myc_miner_biomass_C,     0.0,HUGE(0.0),'update_mycorrhizae #2','cc%myc_miner_biomass_C',FATAL)
     call check_var_range(cc%N_fixer_biomass_C,       0.0,HUGE(0.0),'update_mycorrhizae #2','cc%N_fixer_biomass_C',FATAL)
     call check_var_range(cc%myc_scavenger_biomass_N, 0.0,HUGE(0.0),'update_mycorrhizae #2','cc%myc_scavenger_biomass_N',FATAL)
     call check_var_range(cc%myc_miner_biomass_N,     0.0,HUGE(0.0),'update_mycorrhizae #2','cc%myc_miner_biomass_N',FATAL)
     call check_var_range(cc%N_fixer_biomass_N,       0.0,HUGE(0.0),'update_mycorrhizae #2','cc%N_fixer_biomass_N',FATAL)
     if(is_watch_point()) then
        __DEBUG3__(cc%myc_scavenger_biomass_C,cc%myc_miner_biomass_C,cc%N_fixer_biomass_C)
     endif

     cc%scav_myc_C_reservoir = cc%scav_myc_C_reservoir + d_scav_C_reservoir
     cc%scav_myc_N_reservoir = cc%scav_myc_N_reservoir + d_scav_N_reservoir
     cc%mine_myc_C_reservoir = cc%mine_myc_C_reservoir + d_mine_C_reservoir
     cc%mine_myc_N_reservoir = cc%mine_myc_N_reservoir + d_mine_N_reservoir
     cc%N_fixer_C_reservoir  = cc%N_fixer_C_reservoir  + d_N_fixer_C_reservoir
     cc%N_fixer_N_reservoir  = cc%N_fixer_N_reservoir  + d_N_fixer_N_reservoir

     ! Excess C leaks out of reservoir into root exudates at a time scale of one day
     reservoir_C_leakage = reservoir_C_leakage + (cc%scav_myc_C_reservoir + cc%mine_myc_C_reservoir + cc%N_fixer_C_reservoir)*dt_fast_yr*365
     cc%scav_myc_C_reservoir = cc%scav_myc_C_reservoir - cc%scav_myc_C_reservoir*dt_fast_yr*365
     cc%mine_myc_C_reservoir = cc%mine_myc_C_reservoir - cc%mine_myc_C_reservoir*dt_fast_yr*365
     cc%N_fixer_C_reservoir = cc%N_fixer_C_reservoir - cc%N_fixer_C_reservoir*dt_fast_yr*365

     if(abs(cc%mine_myc_N_reservoir)<1e-20) cc%mine_myc_N_reservoir = 0.0
     if(abs(cc%mine_myc_C_reservoir)<1e-20) cc%mine_myc_C_reservoir = 0.0
     if(abs(cc%scav_myc_N_reservoir)<1e-20) cc%scav_myc_N_reservoir = 0.0
     if(abs(cc%scav_myc_C_reservoir)<1e-20) cc%scav_myc_C_reservoir = 0.0
     call check_var_range(cc%N_fixer_C_reservoir, -1e-10,HUGE(0.0),'update_mycorrhizae #3','cc%N_fixer_C_reservoir',FATAL)
     call check_var_range(cc%N_fixer_N_reservoir, -1e-10,HUGE(0.0),'update_mycorrhizae #3','cc%N_fixer_N_reservoir',FATAL)

     ! Calculate N released to plant
     scav_N_to_plant = cc%scav_myc_N_reservoir*sp%myc_N_to_plant_rate*dt_fast_yr
     mine_N_to_plant = cc%mine_myc_N_reservoir*sp%myc_N_to_plant_rate*dt_fast_yr
     fix_N_to_plant = cc%N_fixer_N_reservoir*sp%myc_N_to_plant_rate*dt_fast_yr

     ! Calculate return on investment for each strategy
     ! Scavenging (AM-style)
     ! slm: *_marginal_gain variables have units [kgN/kgC]
     if(myc_scav_C_efficiency == 0 .OR. .NOT. sp%do_N_scavenging_strategy) then
        myc_scav_marginal_gain = 0
        scav_N_to_plant = 0.0
     else
        if (cc%myc_scavenger_biomass_C > 0) then
           myc_scav_marginal_gain = (max(0.0,scav_N_to_plant)/dt_fast_yr)/(cc%myc_scavenger_biomass_C/myc_scav_C_efficiency/mycorrhizal_turnover_time)
        else ! Use the mycorrhizal N uptake efficiency from myc_scavenger_N_uptake, units of (kgN/kg myc biomass C)
           myc_scav_marginal_gain = myc_scav_efficiency/(dt_fast_yr*myc_scav_C_efficiency*mycorrhizal_turnover_time)
        endif
     endif
     cc%scav_myc_N_reservoir = cc%scav_myc_N_reservoir - scav_N_to_plant

     ! Mycorrhizal N mining (ECM-style)
     if(myc_mine_C_efficiency==0 .OR. .NOT. sp%do_N_mining_strategy) then
        myc_mine_marginal_gain = 0.0
        mine_N_to_plant = 0.0
     else
        if(cc%myc_miner_biomass_C>0) then
           myc_mine_marginal_gain = (max(0.0,mine_N_to_plant)/dt_fast_yr)/(cc%myc_miner_biomass_C/myc_mine_C_efficiency/mycorrhizal_turnover_time)
        else
           myc_mine_marginal_gain = myc_mine_efficiency/(dt_fast_yr*myc_mine_C_efficiency*mycorrhizal_turnover_time)
        endif
     endif
     cc%mine_myc_N_reservoir = cc%mine_myc_N_reservoir - mine_N_to_plant

     ! Root uptake of nitrogen
     if (C_allocation_to_N_acq>0) then
        rhiz_exud_marginal_gain = max(0.001,(root_N_uptake/dt_fast_yr)/(C_allocation_to_N_acq))!+(myc_mine_marginal_gain+myc_scav_marginal_gain)*0.5
     else
        rhiz_exud_marginal_gain = (myc_mine_marginal_gain+myc_scav_marginal_gain)*0.25
     endif

     ! N fixer
     if(N_fixer_C_efficiency == 0 .OR. .NOT. sp%do_N_fixation_strategy) then
        N_fix_marginal_gain = 0.0
        fix_N_to_plant = 0.0
     else
        if(cc%N_fixer_biomass_C>0) then
           N_fix_marginal_gain = (fix_N_to_plant/dt_fast_yr)/(cc%N_fixer_biomass_C/N_fixer_C_efficiency/N_fixer_turnover_time)
        else
           N_fix_marginal_gain = sp%N_fixation_rate*N_fixer_C_efficiency*N_fixer_turnover_time
        endif
     endif
     cc%N_fixer_N_reservoir = cc%N_fixer_N_reservoir-fix_N_to_plant

     ! Apply a smoothing filter to marginal gains, so we can control how fast N strategies change at the ecosystem level
     w = 1.0/(1+sp%tau_smooth_marginal_gain/dt_fast_yr)
     cc%myc_scav_marginal_gain_smoothed  = cc%myc_scav_marginal_gain_smoothed*(1-w)  + myc_scav_marginal_gain*w
     cc%myc_mine_marginal_gain_smoothed  = cc%myc_mine_marginal_gain_smoothed*(1-w)  + myc_mine_marginal_gain*w
     cc%N_fix_marginal_gain_smoothed     = cc%N_fix_marginal_gain_smoothed*(1-w)     + N_fix_marginal_gain*w
     cc%rhiz_exud_marginal_gain_smoothed = cc%rhiz_exud_marginal_gain_smoothed*(1-w) + rhiz_exud_marginal_gain*w

     ! Calculate relative fractions
     if (myc_scav_marginal_gain+N_fix_marginal_gain+myc_mine_marginal_gain+rhiz_exud_marginal_gain>0) then
        mgain = cc%myc_scav_marginal_gain_smoothed  &
              + cc%myc_mine_marginal_gain_smoothed  &
              + cc%N_fix_marginal_gain_smoothed     &
              + cc%rhiz_exud_marginal_gain_smoothed
        myc_scav_exudate_frac = cc%myc_scav_marginal_gain_smoothed/mgain
        myc_mine_exudate_frac = cc%myc_mine_marginal_gain_smoothed/mgain
        N_fixer_exudate_frac  = cc%N_fix_marginal_gain_smoothed/mgain
        rhiz_exud_frac        = cc%rhiz_exud_marginal_gain_smoothed/mgain
     else
        ! Divide evenly if there is no marginal gain.  But this probably only happens if C_allocation_to_N_acq is zero
        myc_scav_exudate_frac = 0.4*0.7
        myc_mine_exudate_frac = 0.3*0.7
        N_fixer_exudate_frac = 0.3*0.7
        rhiz_exud_frac = 0.3
     endif

     if (smooth_N_uptake_C_allocation) then
        N_fixer_C_allocated       = min(C_allocation_to_N_acq*N_fixer_exudate_frac*dt_fast_yr, cc%max_Nfix_allocation*dt_fast_yr*(1.0+sp%alloc_allowed_over_limit))
        miner_myc_C_allocated     = min(C_allocation_to_N_acq*myc_mine_exudate_frac*dt_fast_yr,cc%max_mine_allocation*dt_fast_yr*(1.0+sp%alloc_allowed_over_limit))
        scavenger_myc_C_allocated = min(C_allocation_to_N_acq*myc_scav_exudate_frac*dt_fast_yr,cc%max_scav_allocation*dt_fast_yr*(1.0+sp%alloc_allowed_over_limit))
     else
        N_fixer_C_allocated       = C_allocation_to_N_acq*N_fixer_exudate_frac*dt_fast_yr
        miner_myc_C_allocated     = C_allocation_to_N_acq*myc_mine_exudate_frac*dt_fast_yr
        scavenger_myc_C_allocated = C_allocation_to_N_acq*myc_scav_exudate_frac*dt_fast_yr
     endif

     ! We are accumulating potential (not limited/smoothed) allocation so max allocation can increase over time
     cc%Nfix_alloc_accum = cc%Nfix_alloc_accum+C_allocation_to_N_acq*N_fixer_exudate_frac*dt_fast_yr
     cc%scav_alloc_accum = cc%scav_alloc_accum+C_allocation_to_N_acq*myc_scav_exudate_frac*dt_fast_yr
     cc%mine_alloc_accum = cc%mine_alloc_accum+C_allocation_to_N_acq*myc_mine_exudate_frac*dt_fast_yr

     N_fixer_N_allocated = 0.0
     miner_myc_N_allocated = miner_myc_C_allocated*sp%root_exudate_N_frac
     scavenger_myc_N_allocated = scavenger_myc_C_allocated*sp%root_exudate_N_frac

     ! Make sure N allocation doesn't completely deplete stored N
     if (cc%stored_N<=0.0 .and. N_limits_live_biomass) then
        N_fixer_N_allocated=0.0
        miner_myc_N_allocated=0.0
        scavenger_myc_N_allocated=0.0
     elseif (N_fixer_N_allocated+miner_myc_N_allocated+scavenger_myc_N_allocated > 0.0 .AND. &
               N_limits_live_biomass .AND. &
               N_fixer_N_allocated+miner_myc_N_allocated+scavenger_myc_N_allocated > cc%stored_N*0.1) then
        lim_factor=cc%stored_N*0.1/(N_fixer_N_allocated+miner_myc_N_allocated+scavenger_myc_N_allocated)
        N_fixer_N_allocated=N_fixer_N_allocated*lim_factor
        miner_myc_N_allocated=miner_myc_N_allocated*lim_factor
        scavenger_myc_N_allocated=scavenger_myc_N_allocated*lim_factor
     endif

     cc%scav_myc_C_reservoir = cc%scav_myc_C_reservoir + scavenger_myc_C_allocated
     cc%scav_myc_N_reservoir = cc%scav_myc_N_reservoir + scavenger_myc_N_allocated
     cc%mine_myc_C_reservoir = cc%mine_myc_C_reservoir + miner_myc_C_allocated
     cc%mine_myc_N_reservoir = cc%mine_myc_N_reservoir + miner_myc_N_allocated
     cc%N_fixer_C_reservoir  = cc%N_fixer_C_reservoir  + N_fixer_C_allocated
     cc%N_fixer_N_reservoir  = cc%N_fixer_N_reservoir  + N_fixer_N_allocated

     total_plant_N_uptake = scav_N_to_plant + mine_N_to_plant + fix_N_to_plant + root_N_uptake
     cc%stored_N = cc%stored_N + total_plant_N_uptake
     stored_N_loss = min(C_allocation_to_N_acq*dt_fast_yr*sp%root_exudate_N_frac,cc%stored_N)
     cc%stored_N = cc%stored_N - stored_N_loss
     root_exudate_N = stored_N_loss - scavenger_myc_N_allocated - N_fixer_N_allocated - miner_myc_N_allocated
     if(is_watch_point()) then
        __DEBUG4__(scav_N_to_plant, mine_N_to_plant, fix_N_to_plant, root_N_uptake)
        __DEBUG3__(total_plant_N_uptake, cc%stored_N, root_exudate_N)
     endif
     call check_var_range(root_exudate_N, 0.0,HUGE(0.0),'update_mycorrhizae','root_exudate_N',FATAL)
  else  ! If nitrogen turned off, everything is zero
     scavenger_myc_C_allocated = 0.0
     miner_myc_C_allocated = 0.0
     N_fixer_C_allocated = 0.0
     scavenger_myc_N_allocated = 0.0
     miner_myc_N_allocated = 0.0
     N_fixer_N_allocated = 0.0
     myc_CO2_prod = 0.0
     N_fixation = 0.0
     myc_turnover_C = 0.0
     myc_turnover_N = 0.0
     total_plant_N_uptake = 0.0
     root_exudate_N = 0.0
     scav_N_to_plant = 0.0
     mine_N_to_plant = 0.0
     fix_N_to_plant = 0.0
  endif

  root_exudate_C = C_allocation_to_N_acq*dt_fast_yr - scavenger_myc_C_allocated - miner_myc_C_allocated - N_fixer_C_allocated + reservoir_C_leakage

  if (is_watch_point()) then
     ! __DEBUG1__(current_root_exudation*dt_fast_yr)
     __DEBUG5__(cc%stored_N, root_exudate_N, scavenger_myc_N_allocated, N_fixer_N_allocated, miner_myc_N_allocated)
     __DEBUG1__(root_exudate_C)
     __DEBUG1__(scavenger_myc_C_allocated)
     __DEBUG1__(miner_myc_C_allocated)
     __DEBUG1__(myc_mine_C_uptake)
     __DEBUG1__(N_fixer_C_allocated)
     __DEBUG1__(myc_CO2_prod)
  endif
  end associate ! sp
end subroutine update_mycorrhizae

! ============================================================================
subroutine vegn_carbon_int_lm3(vegn, soil, soilt, theta, diag)
  type(vegn_tile_type), intent(inout) :: vegn
  type(soil_tile_type), intent(inout) :: soil
  real, intent(in) :: soilt ! average temperature of soil for soil carbon decomposition, deg K
  real, intent(in) :: theta ! average soil wetness, unitless
  type(diag_buff_type), intent(inout) :: diag

  ! TODO: possibly move soil-related calculations from calling procedure here,
  !       now that we have soil passed as an argument

  type(vegn_cohort_type), pointer :: c(:) ! for debug only
  real :: md_leaf, md_wood, md_froot ! component of maintenance demand
  real :: md ! plant tissue maintenance, kg C/timestep
  real :: total_root_exudate_C(num_l) ! total root exudate per tile, kgC/m2
  real :: total_root_exudate_N(num_l) ! total root exudate per tile, kgN/m2
  real :: total_myc_Nmin(num_l)
  real, dimension(N_C_TYPES) :: &
      leaf_litt_C, leaf_litt_N, & ! fine surface litter per tile, kgC/m2 and kgN/m2
      wood_litt_C, wood_litt_N    ! coarse surface litter per tile, kgC/m2 and kgN/m2
  real, dimension(num_l, N_C_TYPES) :: &
      root_litt_C, root_litt_N ! root litter per soil layer, kgC/m2 and kgN/m2
  real :: profile(num_l) ! storage for vertical profile of exudates and root litter
  real, dimension(vegn%n_cohorts) :: resp, resl, resr, ress, resg, gpp, npp
  real, dimension(vegn%n_cohorts) :: myc_scav_N_uptake, myc_mine_N_uptake, myc_mine_C_uptake, root_active_N_uptake
  integer :: i, l, N

  real :: root_exudate_frac, C_allocation_to_N_acq
  real,dimension(vegn%n_cohorts) :: scavenger_myc_C_allocated,miner_myc_C_allocated, N_fixer_C_allocated
  real,dimension(vegn%n_cohorts) :: N_fixation, root_exudate_C, root_exudate_N, myc_CO2_prod
  real :: mining_CO2prod,myc_turnover_C,myc_turnover_N
  real :: myc_N_uptake, myc_C_uptake
  real,dimension(vegn%n_cohorts) :: total_plant_N_uptake, scav_N_to_plant, mine_N_to_plant, fix_N_to_plant
  real :: excess_C, current_root_exudation, myc_scav_efficiency, myc_mine_efficiency, total_N_leakage(num_l)
  real :: total_myc_CO2_prod, myc_Nmin ! additional heterotrophic respiration from mycorrhizae and N fixers

  c=>vegn%cohorts(1:vegn%n_cohorts)
  N = vegn%n_cohorts
  if(is_watch_point()) then
     write(*,*)'#### vegn_carbon_int_lm3 input ####'
     __DEBUG2__(soilt,theta)
     __DEBUG1__(c%npp_previous_day)
     __DEBUG1__(c%carbon_gain)
     __DEBUG1__(c%carbon_loss)
     __DEBUG1__(c%bwood_gain)
  endif

  !  update plant carbon
  leaf_litt_C = 0 ; wood_litt_C = 0 ; root_litt_C(:,:) = 0
  leaf_litt_N = 0 ; wood_litt_N = 0 ; root_litt_N(:,:) = 0
  total_root_exudate_C(:) = 0 ; total_root_exudate_N(:) = 0
  total_N_leakage = 0

  total_myc_CO2_prod = 0.0
  myc_scav_N_uptake = 0.0
  myc_scav_efficiency = 0.0
  myc_mine_N_uptake = 0.0
  myc_mine_C_uptake = 0.0
  mining_CO2prod = 0.0
  myc_mine_efficiency = 0.0
  root_active_N_uptake = 0.0

  ! First calculate total N uptake for each strategy
  ! The total will then be divided between cohorts based on their relative
  ! root and mycorrhizal biomass

  if (soil_carbon_option == SOILC_CORPSE_N) then
    call myc_scavenger_N_uptake(soil,vegn,myc_scav_N_uptake,myc_scav_efficiency,dt_fast_yr,update_pools=.TRUE.)
    call myc_miner_N_uptake(soil,vegn,myc_mine_N_uptake,myc_mine_C_uptake,mining_CO2prod,myc_mine_efficiency,dt_fast_yr,update_pools=.TRUE.)
    total_myc_CO2_prod = total_myc_CO2_prod + mining_CO2prod
    call root_N_uptake(soil,vegn,root_active_N_uptake,dt_fast_yr, update_pools=.TRUE.)
  endif

  do i = 1, vegn%n_cohorts
     associate(cc=>vegn%cohorts(i), sp=>spdata(vegn%cohorts(i)%species))

     call plant_respiration(cc,soilt, resp(i), resl(i), resr(i), ress(i))

     gpp(i) = (cc%An_op - cc%An_cl)*mol_C*cc%leafarea;
     npp(i) = gpp(i) - resp(i)

     if(cc%npp_previous_day > 0) then
        resg(i) = sp%GROWTH_RESP*cc%npp_previous_day
        npp(i)  = npp(i)  - sp%GROWTH_RESP*cc%npp_previous_day
        resp(i) = resp(i) + sp%GROWTH_RESP*cc%npp_previous_day
     else
        resg(i) = 0
     endif
     ! accumulate npp for the current day
     cc%npp_previous_day_tmp = cc%npp_previous_day_tmp + npp(i);

     if(sp%dynamic_root_exudation .AND. soil_carbon_option==SOILC_CORPSE_N) then
       ! Initial allocation scheme: root exudation/mycorrhizal allocation depends
       ! on ratio of leaf biomass to max (as determined by N uptake)
       ! Root exudation fraction of NPP limited by some maximum value.
       ! Probably need to rename these parameters and not use a hard-coded value
       root_exudate_frac = min(0.9,sp%root_exudate_frac*cc%nitrogen_stress)
     else
       root_exudate_frac = sp%root_exudate_frac
     endif
     C_allocation_to_N_acq = max(npp(i),0.0)*root_exudate_frac
     cc%carbon_gain = cc%carbon_gain + (npp(i)-C_allocation_to_N_acq)*dt_fast_yr

     ! check if leaves/roots are present and need to be accounted in maintenance
     if(cc%status == LEAF_ON) then
        md_leaf = cc%Pl * sp%alpha_leaf*cc%bliving*dt_fast_yr
        md_froot= cc%Pr * sp%alpha_root*cc%bliving*dt_fast_yr
     else
        md_leaf  = 0
        md_froot = 0
     endif

     ! compute branch and coarse wood losses for tree types
     if (sp%lifeform==FORM_WOODY) then
        md_wood    = 0.6 * cc%bwood * sp%alpha_wood    * dt_fast_yr
     else
        md_wood = 0.0
     endif

     md = md_leaf + md_froot + cc%Psw_alphasw * cc%bliving * dt_fast_yr

     cc%bwood_gain = cc%bwood_gain + cc%Psw_alphasw * cc%bliving * dt_fast_yr;
     cc%bwood_gain = cc%bwood_gain - md_wood;
!     if (cc%bwood_gain < 0.0) cc%bwood_gain=0.0; ! potential non-conservation ?
     cc%carbon_gain = cc%carbon_gain - md;
     cc%carbon_loss = cc%carbon_loss + md; ! used in diagnostics only

     ! Should this N be lost or retranslocated?
     cc%leaf_N = cc%leaf_N - md_leaf /sp%leaf_live_c2n*(1.0-sp%leaf_retranslocation_frac)
     cc%root_N = cc%root_N - md_froot/sp%froot_live_c2n*(1.0-sp%froot_retranslocation_frac)
     cc%wood_N = cc%wood_N - md_wood/sp%wood_c2n

     ! add maintenance demand from leaf and root pools to fast soil carbon
     leaf_litt_C(:) = leaf_litt_C(:) + [sp%fsc_liv,  1-sp%fsc_liv,  0.0]*(md_leaf)*cc%nindivs
     leaf_litt_N(:) = leaf_litt_N(:) + [sp%fsc_liv,  1-sp%fsc_liv,  0.0]*md_leaf/sp%leaf_live_c2n*cc%nindivs*(1.0-sp%leaf_retranslocation_frac)
     wood_litt_C(:) = wood_litt_C(:) + cc%nindivs * agf_bs * &
             [sp%fsc_wood, 1-sp%fsc_wood, 0.0]*md_wood
     wood_litt_N(:) = wood_litt_N(:) + cc%nindivs * agf_bs * &
             [sp%fsc_wood, 1-sp%fsc_wood, 0.0]*md_wood/sp%wood_c2n
     call cohort_root_litter_profile(cc, dz, profile)
     do l = 1, num_l
        root_litt_C(l,:) = root_litt_C(l,:) + profile(l)*cc%nindivs * ( &
             [sp%fsc_froot, 1-sp%fsc_froot, 0.0 ]*md_froot + &
             [sp%fsc_wood,  1-sp%fsc_wood,  0.0 ]*md_wood*(1-agf_bs) &
             )
        root_litt_N(l,:) = root_litt_N(l,:) + profile(l)*cc%nindivs*( &
             [sp%fsc_froot, 1-sp%fsc_froot, 0.0 ]*md_froot/sp%froot_live_c2n*(1.0-sp%froot_retranslocation_frac) + &
             [sp%fsc_wood,  1-sp%fsc_wood,  0.0 ]*md_wood*(1-agf_bs)/sp%wood_c2n &
             )
     enddo

     vegn%veg_in  = vegn%veg_in  + npp(i)*cc%nindivs*dt_fast_yr;
     vegn%veg_out = vegn%veg_out + (md_leaf+md_froot+md_wood)*cc%nindivs;

     ! update cohort age
     cc%age = cc%age + dt_fast_yr

     ! Mycorrhizal N uptake
     call update_mycorrhizae(cc, soilt, &
          C_allocation_to_N_acq=C_allocation_to_N_acq,&
          myc_scav_N_uptake=myc_scav_N_uptake(i),myc_mine_C_uptake=myc_mine_C_uptake(i),myc_mine_N_uptake=myc_mine_N_uptake(i),&
          root_N_uptake=root_active_N_uptake(i),&
          myc_scav_efficiency=myc_scav_efficiency,myc_mine_efficiency=myc_mine_efficiency,&
          scavenger_myc_C_allocated=scavenger_myc_C_allocated(i),miner_myc_C_allocated=miner_myc_C_allocated(i),&
          N_fixer_C_allocated=N_fixer_C_allocated(i),&
          root_exudate_C=root_exudate_C(i), root_exudate_N=root_exudate_N(i),&
          myc_CO2_prod=myc_CO2_prod(i),N_fixation=N_fixation(i),&
          scav_N_to_plant=scav_N_to_plant(i), mine_N_to_plant=mine_N_to_plant(i), fix_N_to_plant=fix_N_to_plant(i), &
          total_plant_N_uptake=total_plant_N_uptake(i), &
          myc_turnover_C=myc_turnover_C,myc_turnover_N=myc_turnover_N,myc_Nmin=myc_Nmin)

     ! First add mycorrhizal and N fixer turnover to soil C pools
     do l = 1, num_l
        root_litt_C(l,:) = root_litt_C(l,:) + profile(l)*cc%nindivs*([ 0.0, myc_turnover_C*deadmic_slow_frac, myc_turnover_C*(1-deadmic_slow_frac) ])
        root_litt_N(l,:) = root_litt_N(l,:) + profile(l)*cc%nindivs*([ 0.0, myc_turnover_N*deadmic_slow_frac, myc_turnover_N*(1-deadmic_slow_frac) ])
     enddo

     ! To do: Figure out excess C if mycorrhizae N-limited and add back to root exudate C
     call cohort_root_exudate_profile(cc,dz,profile)
     total_root_exudate_C(:) = total_root_exudate_C(:) + profile(:)*root_exudate_C(i)*cc%nindivs
     total_root_exudate_N(:) = total_root_exudate_N(:) + profile(:)*root_exudate_N(i)*cc%nindivs

     ! To prevent excess stored N buildup under high soil N, leak stored N when stress is zero
     if(cc%nitrogen_stress<=0.05 .AND. cc%stored_N>0 .AND. soil_carbon_option == SOILC_CORPSE_N) then
       total_N_leakage(:) = total_N_leakage(:) + cc%stored_N*excess_stored_N_leakage_rate*profile(:)*cc%nindivs
       cc%stored_N=cc%stored_N - cc%stored_N*excess_stored_N_leakage_rate*dt_fast_yr
     endif

     total_myc_CO2_prod = total_myc_CO2_prod + myc_CO2_prod(i)*cc%nindivs
     total_myc_Nmin(:) = total_myc_Nmin(:) + profile(:)*myc_Nmin*cc%nindivs
     end associate
  enddo

  soil%gross_nitrogen_flux_into_tile = soil%gross_nitrogen_flux_into_tile + sum(N_fixation(1:N)*c(1:N)%nindivs)

  ! fsc_in and ssc_in updated in add_root_exudates
  call add_root_exudates(soil,total_root_exudate_C,total_root_exudate_N,total_myc_Nmin,total_N_leakage*dt_fast_yr)

  ! add litter accumulated over the cohorts
  call add_soil_carbon(soil, vegn, leaf_litt_C, wood_litt_C, root_litt_C, &
                                   leaf_litt_N, wood_litt_N, root_litt_N  )

  if(is_watch_point()) then
     write(*,*)'#### vegn_carbon_int_lm3 output ####'
     __DEBUG1__(c%species)
     __DEBUG1__(c%bl)
     __DEBUG1__(c%br)
     __DEBUG1__(c%bsw)
     __DEBUG1__(c%bwood)
     __DEBUG1__(c%total_N)
     __DEBUG1__(c%stored_N)
     __DEBUG1__(c%leaf_N)
     __DEBUG1__(c%wood_N)
     __DEBUG1__(c%root_N)
     __DEBUG1__(c%sapwood_N)
     __DEBUG1__(c%nitrogen_stress)
     __DEBUG1__(c%An_op)
     __DEBUG1__(c%An_cl)
     __DEBUG1__(c%lai)
     __DEBUG1__(npp)
     __DEBUG1__(gpp)
     __DEBUG1__(resp)
     __DEBUG1__(resl)
     __DEBUG1__(resr)
     __DEBUG1__(resg)
     __DEBUG1__(c%carbon_gain)
     __DEBUG1__(c%carbon_loss)
     __DEBUG1__(c%bwood_gain)
     __DEBUG1__(sum(total_root_exudate_C))
     __DEBUG1__(total_myc_CO2_prod)
  endif

  ! update soil carbon
  call Dsdt(vegn, soil, diag, soilt, theta)

  vegn%rh = vegn%rh + total_myc_CO2_prod/dt_fast_yr

  ! NEP is equal to NPP minus soil respiration
  vegn%nep = sum(npp(1:N)*c(1:N)%nindivs) - vegn%rh

  call update_soil_pools(vegn, soil)
  vegn%age = vegn%age + dt_fast_yr;


  ! ---- diagnostic section
  call send_cohort_data(id_gpp, diag, c(1:N), gpp(1:N), weight=c(1:N)%nindivs, op=OP_SUM)
  call send_cohort_data(id_npp, diag, c(1:N), npp(1:N), weight=c(1:N)%nindivs, op=OP_SUM)
  call send_tile_data(id_nep,vegn%nep,diag)
  call send_tile_data(id_litter,vegn%litter,diag)
  call send_cohort_data(id_resp, diag, c(1:N), resp(1:N), weight=c(1:N)%nindivs, op=OP_SUM)
  call send_cohort_data(id_resl, diag, c(1:N), resl(1:N), weight=c(1:N)%nindivs, op=OP_SUM)
  call send_cohort_data(id_resr, diag, c(1:N), resr(1:N), weight=c(1:N)%nindivs, op=OP_SUM)
  call send_cohort_data(id_ress, diag, c(1:N), ress(1:N), weight=c(1:N)%nindivs, op=OP_SUM)
  call send_cohort_data(id_resg, diag, c(1:N), resg(1:N), weight=c(1:N)%nindivs, op=OP_SUM)
  call send_tile_data(id_soilt,soilt,diag)
  call send_tile_data(id_theta,theta,diag)

  call send_cohort_data(id_mycorrhizal_scav_allocation,diag,c(1:N),scavenger_myc_C_allocated(1:N)/dt_fast_yr,weight=c(1:N)%nindivs, op=OP_SUM)
  call send_cohort_data(id_mycorrhizal_mine_allocation,diag,c(1:N),miner_myc_C_allocated(1:N)/dt_fast_yr,weight=c(1:N)%nindivs, op=OP_SUM)
  call send_cohort_data(id_N_fixer_allocation,diag,c(1:N),N_fixer_C_allocated(1:N)/dt_fast_yr,weight=c(1:N)%nindivs, op=OP_SUM)
  call send_cohort_data(id_myc_scav_marginal_gain,diag,c(1:N),c(1:N)%myc_scav_marginal_gain_smoothed,weight=c(1:N)%nindivs, op=OP_SUM)
  call send_cohort_data(id_myc_mine_marginal_gain,diag,c(1:N),c(1:N)%myc_mine_marginal_gain_smoothed,weight=c(1:N)%nindivs, op=OP_SUM)
  call send_cohort_data(id_N_fix_marginal_gain,diag,c(1:N),c(1:N)%N_fix_marginal_gain_smoothed,weight=c(1:N)%nindivs, op=OP_SUM)
  call send_cohort_data(id_rhiz_exud_marginal_gain,diag,c(1:N),c(1:N)%rhiz_exud_marginal_gain_smoothed,weight=c(1:N)%nindivs, op=OP_SUM)
  call send_cohort_data(id_rhiz_exudation,diag,c(1:N),root_exudate_C(1:N)/dt_fast_yr,weight=c(1:N)%nindivs, op=OP_SUM)
  call send_cohort_data(id_nitrogen_stress,diag,c(1:N),c(1:N)%nitrogen_stress,weight=c(1:N)%nindivs, op=OP_SUM)
  call send_cohort_data(id_total_plant_N_uptake,diag,c(1:N),total_plant_N_uptake(1:N)/dt_fast_yr,weight=c(1:N)%nindivs, op=OP_SUM)

  call send_cohort_data(id_myc_scavenger_N_uptake,diag,c(1:N),myc_scav_N_uptake(1:N)/dt_fast_yr,weight=c(1:N)%nindivs, op=OP_SUM)
  call send_cohort_data(id_myc_miner_N_uptake,diag,c(1:N),myc_mine_N_uptake(1:N)/dt_fast_yr,weight=c(1:N)%nindivs, op=OP_SUM)
  call send_cohort_data(id_symbiotic_N_fixation,diag,c(1:N),N_fixation(1:N)/dt_fast_yr,weight=c(1:N)%nindivs, op=OP_SUM)
  call send_cohort_data(id_active_root_N_uptake,diag,c(1:N), root_active_N_uptake(1:N)/dt_fast_yr,weight=c(1:N)%nindivs, op=OP_SUM)

  call send_cohort_data(id_scav_plant_N_uptake,diag,c(1:N),scav_N_to_plant(1:N)/dt_fast_yr,weight=c(1:N)%nindivs, op=OP_SUM)
  call send_cohort_data(id_mine_plant_N_uptake,diag,c(1:N),mine_N_to_plant(1:N)/dt_fast_yr,weight=c(1:N)%nindivs, op=OP_SUM)
  call send_cohort_data(id_fix_plant_N_uptake,diag,c(1:N),fix_N_to_plant(1:N)/dt_fast_yr,weight=c(1:N)%nindivs, op=OP_SUM)

  call send_cohort_data(id_mycorrhizal_scav_N_res,  diag, c(1:N), c(1:N)%scav_myc_N_reservoir, weight=c(1:N)%nindivs, op=OP_SUM)
  call send_cohort_data(id_mycorrhizal_scav_C_res,  diag, c(1:N), c(1:N)%scav_myc_C_reservoir, weight=c(1:N)%nindivs, op=OP_SUM)
  call send_cohort_data(id_mycorrhizal_mine_N_res,  diag, c(1:N), c(1:N)%mine_myc_N_reservoir, weight=c(1:N)%nindivs, op=OP_SUM)
  call send_cohort_data(id_mycorrhizal_mine_C_res,  diag, c(1:N), c(1:N)%mine_myc_C_reservoir, weight=c(1:N)%nindivs, op=OP_SUM)
  call send_cohort_data(id_Nfix_N_res,              diag, c(1:N), c(1:N)%N_fixer_N_reservoir,  weight=c(1:N)%nindivs, op=OP_SUM)
  call send_cohort_data(id_Nfix_C_res,              diag, c(1:N), c(1:N)%N_fixer_C_reservoir,  weight=c(1:N)%nindivs, op=OP_SUM)

  call send_cohort_data(id_N_fix_alloc_smoothed,    diag, c(1:N), c(1:N)%max_Nfix_allocation,  weight=c(1:N)%nindivs, op=OP_SUM)
  call send_cohort_data(id_myc_mine_alloc_smoothed, diag, c(1:N), c(1:N)%max_mine_allocation,  weight=c(1:N)%nindivs, op=OP_SUM)
  call send_cohort_data(id_myc_scav_alloc_smoothed, diag, c(1:N), c(1:N)%max_scav_allocation,  weight=c(1:N)%nindivs, op=OP_SUM)

end subroutine vegn_carbon_int_lm3


! ============================================================================
subroutine vegn_carbon_int_ppa (vegn, soil, tsoil, theta, diag)
  ! TODO: possibly get rid of tsoil, theta, since they can be calculated here
  type(vegn_tile_type), intent(inout) :: vegn
  type(soil_tile_type), intent(inout) :: soil
  real, intent(in) :: tsoil ! average temperature of soil for soil carbon decomposition, deg K
  real, intent(in) :: theta ! average soil wetness, unitless
  type(diag_buff_type), intent(inout) :: diag

  real :: md_branch_sw ! in ppa we are losing and replacing branchwood
  real :: md_bsw, md_bsw_branch ! we are also overturning sapwood in grasses
  real :: deltaBL, deltaBR, deltaNL, deltaNR ! leaf and fine root carbon tendencies
  integer :: i, l, M
  real :: NSC_supply,LR_demand,LR_deficit
  real :: R_days,fNSC,fLFR,fStem,fSeed
  type(vegn_cohort_type), pointer :: c(:) ! for debug only
  ! accumulators of total input to root litter and soil carbon
  real, dimension(N_C_TYPES) :: &
      leaf_litt_C, leaf_litt_N, & ! fine surface litter per tile, kgC/m2
      wood_litt_C, wood_litt_N    ! coarse surface litter per tile, kgC/m2
  real, dimension (num_l, N_C_TYPES) :: &
      root_litt_C, root_litt_N    ! root litter per soil layer, kgC/m2
  real :: profile(num_l) ! storage for vertical profile of exudates and root litter
  real, dimension(vegn%n_cohorts) :: &
      resp, resl, resr, ress, resg, gpp, npp, &
      scavenger_myc_C_allocated,miner_myc_C_allocated, N_fixer_C_allocated, &
      myc_scav_N_uptake, myc_mine_N_uptake, myc_mine_C_uptake, root_active_N_uptake, &
      myc_scav_marginal_gain,myc_mine_marginal_gain, N_fix_marginal_gain, rhiz_exud_marginal_gain, &
      N_fixation, root_exudate_C, root_exudate_N, myc_CO2_prod, &
      total_plant_N_uptake, scav_N_to_plant, mine_N_to_plant, fix_N_to_plant
  real, dimension(num_l) :: &
      total_root_exudate_C, & ! total root exudate per tile, kgC/m2
      total_root_exudate_N, & ! total root exudate per tile, kgN/m2
      total_N_leakage, total_myc_Nmin
  real :: excess_C, current_root_exudation, myc_scav_efficiency, myc_mine_efficiency
  real :: mining_CO2prod,myc_turnover_C,myc_turnover_N
  real :: root_exudate_frac, C_allocation_to_N_acq
  real :: total_myc_CO2_prod, myc_Nmin ! additional heterotrophic respiration from mycorrhizae and N fixers


  c=>vegn%cohorts(1:vegn%n_cohorts)
  M = vegn%n_cohorts

  if(is_watch_point()) then
     write(*,*)'#### vegn_carbon_int_ppa input ####'
     __DEBUG2__(tsoil,theta)
     __DEBUG1__(c%species)
     __DEBUG1__(c%status)
     __DEBUG1__(c%An_op)
     __DEBUG1__(c%An_cl)
     __DEBUG1__(c%bl_max)
     __DEBUG1__(c%br_max)
     __DEBUG1__(c%dbh)
     __DEBUG1__(c%height)
     __DEBUG1__(c%leafarea)
     __DEBUG1__(c%carbon_gain)
     __DEBUG1__(c%carbon_loss)
  endif
  ! update plant carbon for all cohorts
  leaf_litt_C = 0; wood_litt_C = 0; root_litt_C = 0; total_root_exudate_C = 0
  leaf_litt_N = 0; wood_litt_N = 0; root_litt_N = 0; total_root_exudate_N = 0
  total_N_leakage = 0
  myc_scav_N_uptake = 0.0
  myc_scav_efficiency = 0.0
  myc_mine_N_uptake = 0.0
  myc_mine_C_uptake = 0.0
  mining_CO2prod = 0.0
  myc_mine_efficiency = 0.0
  root_active_N_uptake = 0.0

! 20170617:
  total_myc_CO2_prod = 0.0; total_myc_Nmin = 0.0
  if (soil_carbon_option == SOILC_CORPSE_N) then
    call myc_scavenger_N_uptake(soil,vegn,myc_scav_N_uptake,myc_scav_efficiency,dt_fast_yr,update_pools=.TRUE.)
    call myc_miner_N_uptake(soil,vegn,myc_mine_N_uptake,myc_mine_C_uptake,mining_CO2prod,myc_mine_efficiency,dt_fast_yr,update_pools=.TRUE.)
    total_myc_CO2_prod = total_myc_CO2_prod + mining_CO2prod
    call root_N_uptake(soil,vegn,root_active_N_uptake,dt_fast_yr, update_pools=.TRUE.)
  endif

  do i = 1, vegn%n_cohorts
     associate ( cc => vegn%cohorts(i), sp => spdata(vegn%cohorts(i)%species))

     ! that was in eddy_npp_PPA
     call plant_respiration(cc,tsoil,resp(i),resl(i),resr(i),ress(i))

     gpp(i)  = (cc%An_op - cc%An_cl)*mol_C*cc%leafarea
     ! growth respiration comes from nsc pool not gpp now
     resg(i) = cc%growth_previous_day_tmp ! kg C per individual per year

     cc%growth_previous_day = cc%growth_previous_day - resg(i)*dt_fast_yr

     if (is_watch_point()) then
        write(*,'("####### vegn_carbon_int_ppa #0, cohort ",i2.2)') i
        __DEBUG1__(cc%nindivs)
        __DEBUG5__(cc%bl, cc%br, cc%bsw, cc%bwood, cc%nsc)
        __DEBUG5__(cc%leaf_N, cc%root_N, cc%sapwood_N, cc%wood_N, cc%stored_N)
     endif
     cc%carbon_gain  = cc%carbon_gain + gpp(i)*dt_fast_yr - resp(i)*dt_fast_yr
     if (is_watch_point()) then
        __DEBUG1__(cc%carbon_gain)
     endif

     resp(i) = resp(i) + resg(i) ! add growth respiration and maintenance respiration
     npp(i)  = gpp(i) - resp(i)

     ! Weng, 2013-01-28
     ! Turnover regardless of STATUS
     deltaBL = cc%bl * sp%alpha_leaf * dt_fast_yr
     deltaBR = cc%br * sp%alpha_root * dt_fast_yr
     deltaNL = deltaBL/sp%leaf_live_c2n
     deltaNR = deltaBR/sp%froot_live_c2n

     cc%bl = cc%bl - deltaBL
     cc%br = cc%br - deltaBR

     if(N_limits_live_biomass.and.cc%leaf_N-deltaNL<0) then
         __DEBUG2__(cc%leaf_N,deltaNL)
     endif

     ! 20170617: retranslocate part of N back to storage and reduce nitrogen pools; put
     !           retranslocated N into storage
     cc%leaf_N   = cc%leaf_N - deltaNL
     cc%root_N   = cc%root_N - deltaNR
     if (soil_carbon_option==SOILC_CORPSE_N) &
         cc%stored_N = cc%stored_N + deltaBL/sp%leaf_live_c2n*sp%leaf_retranslocation_frac   &
                                   + deltaBR/sp%froot_live_c2n*sp%froot_retranslocation_frac
     if (is_watch_point()) then
        __DEBUG4__(deltaBL,deltaBR,deltaNL,deltaNR)
        __DEBUG5__(cc%bl, cc%br, cc%leaf_N, cc%root_N, cc%stored_N)
     endif
     if(N_limits_live_biomass) then
          call check_var_range(cc%leaf_N,0.0,HUGE(1.0),'vegn_carbon_int_ppa','cc%leaf_N',FATAL)
          call check_var_range(cc%stored_N,0.0,HUGE(1.0),'vegn_carbon_int_ppa #1','cc%stored_N',FATAL)
     endif
     ! compute branch and coarse wood losses for tree types
     md_branch_sw = 0.0
     if (spdata(cc%species)%lifeform == FORM_WOODY.and.&
         cc%height > sp%branch_loss_height) then
        ! ens 02/07/17
        md_branch_sw = Max(cc%brsw,0.0) * sp%alpha_wood * dt_fast_yr
     endif

     call check_var_range(cc%bsw,  0.0,HUGE(1.0), 'vegn_carbon_int_ppa', 'cc%bsw', FATAL)
     call check_var_range(cc%brsw, 0.0,cc%bsw,    'vegn_carbon_int_ppa', 'cc%brsw',FATAL)

     ! ens, isa, slm 2017-08-03: compute overturning of sapwood in grasses.
     md_bsw = 0.0; md_bsw_branch = 0.0
     if (spdata(cc%species)%lifeform == FORM_GRASS) then
        if (cc%bsw > 0) then
           md_bsw = Max(cc%bsw,0.0) * sp%alpha_wood * dt_fast_yr
           md_bsw_branch = cc%brsw/cc%bsw * md_bsw ! overturning in branches, which are part of bsw. Do they even exist in grass?
        endif
     endif

     ! brsw is a part of bsw, so when we lose branches, we need to reduce both
     cc%brsw = cc%brsw - md_branch_sw - md_bsw_branch
     cc%bsw  = cc%bsw  - md_branch_sw - md_bsw
     ! 20170617: update N pools. For now, assuming that there is no N re-translocation during
     ! either branch loss or turnover of wood and sapwood.
     if (soil_carbon_option==SOILC_CORPSE_N) then
        cc%sapwood_N = cc%sapwood_N - (md_branch_sw+md_bsw)/sp%sapwood_c2n
     endif

     ! accumulate liter and soil carbon inputs across all cohorts
     ! 20170617: deposit lost N into litter and soil
     leaf_litt_C(:) = leaf_litt_C(:) + [sp%fsc_liv,  1-sp%fsc_liv,  0.0]*deltaBL*cc%nindivs
     leaf_litt_N(:) = leaf_litt_N(:) + [sp%fsc_liv,  1-sp%fsc_liv,  0.0]*deltaNL*cc%nindivs*(1.0-sp%leaf_retranslocation_frac)
     wood_litt_C(:) = wood_litt_C(:) + [sp%fsc_wood, 1-sp%fsc_wood, 0.0]*(md_branch_sw+md_bsw)*cc%nindivs
!     wood_litt_N(:) = wood_litt_N(:) + cc%nindivs * agf_bs * & -- FIXME: do we need agf_bs in both C and N wood litter?
     if (soil_carbon_option==SOILC_CORPSE_N) wood_litt_N(:) = wood_litt_N(:) + cc%nindivs * &
             [sp%fsc_wood, 1-sp%fsc_wood, 0.0]*(md_branch_sw+md_bsw)/sp%sapwood_c2n
     call cohort_root_litter_profile(cc, dz, profile)
     do l = 1, num_l
        root_litt_C(l,:) = root_litt_C(l,:) + profile(l)*cc%nindivs* &
             [sp%fsc_froot, 1-sp%fsc_froot, 0.0 ]*deltaBR
        root_litt_N(l,:) = root_litt_N(l,:) + profile(l)*cc%nindivs* &
             [sp%fsc_froot, 1-sp%fsc_froot, 0.0 ]*deltaNR*(1.0-sp%froot_retranslocation_frac)
     enddo

     ! resg(i) = 0 ! slm: that doesn't make much sense to me... why?
! 20170617:
     ! Mycorrhizal N uptake

    ! tau_nsc_exudate is currently a time value, which means it can't be set to zero
    ! Maybe redefine as a rate per year?
    ! Allowing to be zero for now.
     if(sp%tau_nsc_exudate>0) then
        C_allocation_to_N_acq = max(cc%nsc,0.0)/sp%tau_nsc_exudate !This is a rate per year, not per time step
     else
        C_allocation_to_N_acq = 0.0
     endif
     if (sp%dynamic_root_exudation .AND. soil_carbon_option==SOILC_CORPSE_N) then
        ! 20170617: modify frac for exudate depending on N state
         C_allocation_to_N_acq = C_allocation_to_N_acq*cc%nitrogen_stress
         ! we possibly should impose lower and upper limits on exudate stress factor. Currently,
         ! there is a lower limit on nitrogen_stress, so exudates are never 0. Also nitrogen_stress
         ! cannot rise above 2 in current formulation, but we should be careful if the definition
         ! of stress changes.
         ! N exudates are calculated in update_mycorrhizae.
     endif

     cc%nsc = cc%nsc - C_allocation_to_N_acq*dt_fast_yr

     ! this updates N storage in the plant
     ! 20170617: This includes allocation to all N acquisition plus exudates
     call update_mycorrhizae(cc, tsoil, &
           C_allocation_to_N_acq=C_allocation_to_N_acq, &
           myc_scav_N_uptake=myc_scav_N_uptake(i), myc_mine_C_uptake=myc_mine_C_uptake(i), myc_mine_N_uptake=myc_mine_N_uptake(i), &
           root_N_uptake=root_active_N_uptake(i), &
           myc_scav_efficiency=myc_scav_efficiency, myc_mine_efficiency=myc_mine_efficiency, &
           scavenger_myc_C_allocated=scavenger_myc_C_allocated(i), miner_myc_C_allocated=miner_myc_C_allocated(i), &
           N_fixer_C_allocated=N_fixer_C_allocated(i), &
           root_exudate_C=root_exudate_C(i), root_exudate_N=root_exudate_N(i), &
           myc_CO2_prod=myc_CO2_prod(i), N_fixation=N_fixation(i), &
           scav_N_to_plant=scav_N_to_plant(i), mine_N_to_plant=mine_N_to_plant(i), fix_N_to_plant=fix_N_to_plant(i), &
           total_plant_N_uptake=total_plant_N_uptake(i), &
           myc_turnover_C=myc_turnover_C, myc_turnover_N=myc_turnover_N, myc_Nmin=myc_Nmin)
     ! root_exudate_C(i) is updated inside update_mycorrhizae

     ! First add mycorrhizal and N fixer turnover to soil C pools
     do l = 1, num_l
        root_litt_C(l,:) = root_litt_C(l,:) + profile(l)*cc%nindivs*([ 0.0, myc_turnover_C*deadmic_slow_frac, myc_turnover_C*(1-deadmic_slow_frac) ])
        root_litt_N(l,:) = root_litt_N(l,:) + profile(l)*cc%nindivs*([ 0.0, myc_turnover_N*deadmic_slow_frac, myc_turnover_N*(1-deadmic_slow_frac) ])
     enddo
!
! bns: 0.05 must be consistent with 0.05 in stress calculation (NOT by-species)
     ! To do: Figure out excess C if mycorrhizae N-limited and add back to root exudate C
     call cohort_root_exudate_profile(cc,dz,profile)
     total_root_exudate_C(:) = total_root_exudate_C(:) + profile(:)*root_exudate_C(i)*cc%nindivs
     total_root_exudate_N(:) = total_root_exudate_N(:) + profile(:)*root_exudate_N(i)*cc%nindivs

     ! To prevent excess stored N buildup under high soil N, leak stored N when stress is zero
     if(cc%nitrogen_stress<=min_N_stress .AND. cc%stored_N>0 .AND. soil_carbon_option == SOILC_CORPSE_N) then
        total_N_leakage(:) = total_N_leakage(:) + cc%stored_N*excess_stored_N_leakage_rate*profile(:)*cc%nindivs*dt_fast_yr
        cc%stored_N=cc%stored_N - cc%stored_N*excess_stored_N_leakage_rate*dt_fast_yr
     endif

     total_myc_CO2_prod = total_myc_CO2_prod + myc_CO2_prod(i)*cc%nindivs
     total_myc_Nmin(:) = total_myc_Nmin(:) + profile(:)*myc_Nmin*cc%nindivs

     ! increment cohort age
     cc%age = cc%age + dt_fast_yr

     if(N_limits_live_biomass) &
          call check_var_range(cc%stored_N,0.0,HUGE(1.0),'vegn_carbon_int_ppa #2','cc%stored_N',FATAL)

     end associate
  enddo
! 20170617:
  soil%gross_nitrogen_flux_into_tile = soil%gross_nitrogen_flux_into_tile + sum(N_fixation(1:M)*c(1:M)%nindivs)

  ! add litter and exudates accumulated over the cohorts
  ! 20170617: revisit exudates for allocation to different kind of N startegies
  call add_root_exudates(soil, total_root_exudate_C, total_root_exudate_N, total_myc_Nmin, total_N_leakage)
  call add_soil_carbon(soil, vegn, leaf_litt_C, wood_litt_C, root_litt_C, &
                                   leaf_litt_N, wood_litt_N, root_litt_N  )
  ! update soil carbon
  call Dsdt(vegn, soil, diag, tsoil, theta)
  vegn%rh = vegn%rh + total_myc_CO2_prod/dt_fast_yr


  ! NEP is equal to GPP minus autotrophic and soil respiration
!  vegn%nep = sum(npp(1:M)*c(1:M)%nindivs) - vegn%rh
  vegn%nep = sum((gpp(1:M)-resp(1:M))*c(1:M)%nindivs) - vegn%rh

  call update_soil_pools(vegn, soil)
  vegn%age = vegn%age + dt_fast_yr

  if(is_watch_point()) then
     write(*,*)'#### vegn_carbon_int_ppa output ####'
     __DEBUG1__(resl)
     __DEBUG1__(resr)
     __DEBUG1__(resg)
     __DEBUG1__(resp)
     __DEBUG1__(npp)
     __DEBUG1__(gpp)
     __DEBUG1__(c%carbon_gain)
     __DEBUG1__(c%bl)
     __DEBUG1__(c%blv)
     __DEBUG1__(c%br)
     __DEBUG1__(c%bsw)
     __DEBUG1__(c%bwood)
     __DEBUG1__(c%nsc)
     write(*,*)'#### end of vegn_carbon_int_ppa output ####'
  endif

! ------ diagnostic section
  call send_cohort_data(id_gpp,  diag, c(1:M), gpp(1:M),  weight=c(1:M)%nindivs, op=OP_SUM)
  call send_cohort_data(id_npp,  diag, c(1:M), npp(1:M),  weight=c(1:M)%nindivs, op=OP_SUM)
  call send_tile_data(id_nep,vegn%nep,diag)
  call send_tile_data(id_litter,vegn%litter,diag)
  call send_cohort_data(id_resp, diag, c(1:M), resp(1:M), weight=c(1:M)%nindivs, op=OP_SUM)
  call send_cohort_data(id_resl, diag, c(1:M), resl(1:M), weight=c(1:M)%nindivs, op=OP_SUM)
  call send_cohort_data(id_resr, diag, c(1:M), resr(1:M), weight=c(1:M)%nindivs, op=OP_SUM)
  call send_cohort_data(id_ress, diag, c(1:M), ress(1:M), weight=c(1:M)%nindivs, op=OP_SUM)
  call send_cohort_data(id_resg, diag, c(1:M), resg(1:M), weight=c(1:M)%nindivs, op=OP_SUM)
  call send_tile_data(id_soilt,tsoil,diag)
  call send_tile_data(id_theta,theta,diag)
  call send_cohort_data(id_age, diag, c(1:M), c(1:M)%age, weight=c(1:M)%nindivs, op=OP_AVERAGE)

  call send_cohort_data(id_mycorrhizal_scav_allocation, diag, c(1:M), scavenger_myc_C_allocated(1:M)/dt_fast_yr, weight=c(1:M)%nindivs, op=OP_SUM)
  call send_cohort_data(id_mycorrhizal_mine_allocation, diag, c(1:M), miner_myc_C_allocated(1:M)/dt_fast_yr,     weight=c(1:M)%nindivs, op=OP_SUM)
  call send_cohort_data(id_N_fixer_allocation,          diag, c(1:M), N_fixer_C_allocated(1:M)/dt_fast_yr,       weight=c(1:M)%nindivs, op=OP_SUM)
  call send_cohort_data(id_myc_scav_marginal_gain,      diag, c(1:M), c(1:M)%myc_scav_marginal_gain_smoothed,    weight=c(1:M)%nindivs, op=OP_SUM)
  call send_cohort_data(id_myc_mine_marginal_gain,      diag, c(1:M), c(1:M)%myc_mine_marginal_gain_smoothed,    weight=c(1:M)%nindivs, op=OP_SUM)
  call send_cohort_data(id_N_fix_marginal_gain,         diag, c(1:M), c(1:M)%N_fix_marginal_gain_smoothed,       weight=c(1:M)%nindivs, op=OP_SUM)
  call send_cohort_data(id_rhiz_exud_marginal_gain,     diag, c(1:M), c(1:M)%rhiz_exud_marginal_gain_smoothed,   weight=c(1:M)%nindivs, op=OP_SUM)
  call send_cohort_data(id_rhiz_exudation,              diag, c(1:M), root_exudate_C(1:M)/dt_fast_yr,            weight=c(1:M)%nindivs, op=OP_SUM)
  call send_cohort_data(id_nitrogen_stress,             diag, c(1:M), c(1:M)%nitrogen_stress,                    weight=c(1:M)%nindivs, op=OP_AVERAGE)
  call send_cohort_data(id_total_plant_N_uptake,        diag, c(1:M), total_plant_N_uptake(1:M)/dt_fast_yr,      weight=c(1:M)%nindivs, op=OP_SUM)

  call send_cohort_data(id_myc_scavenger_N_uptake,      diag, c(1:M), myc_scav_N_uptake(1:M)/dt_fast_yr,         weight=c(1:M)%nindivs, op=OP_SUM)
  call send_cohort_data(id_myc_miner_N_uptake,          diag, c(1:M), myc_mine_N_uptake(1:M)/dt_fast_yr,         weight=c(1:M)%nindivs, op=OP_SUM)
  call send_cohort_data(id_symbiotic_N_fixation,        diag, c(1:M), N_fixation(1:M)/dt_fast_yr,                weight=c(1:M)%nindivs, op=OP_SUM)
  call send_cohort_data(id_active_root_N_uptake,        diag, c(1:M), root_active_N_uptake(1:M)/dt_fast_yr,      weight=c(1:M)%nindivs, op=OP_SUM)

  call send_cohort_data(id_scav_plant_N_uptake,         diag, c(1:M), scav_N_to_plant(1:M)/dt_fast_yr,           weight=c(1:M)%nindivs, op=OP_SUM)
  call send_cohort_data(id_mine_plant_N_uptake,         diag, c(1:M), mine_N_to_plant(1:M)/dt_fast_yr,           weight=c(1:M)%nindivs, op=OP_SUM)
  call send_cohort_data(id_fix_plant_N_uptake,          diag, c(1:M), fix_N_to_plant(1:M)/dt_fast_yr,            weight=c(1:M)%nindivs, op=OP_SUM)

  call send_cohort_data(id_mycorrhizal_scav_N_res,      diag, c(1:M), c(1:M)%scav_myc_N_reservoir,               weight=c(1:M)%nindivs, op=OP_SUM)
  call send_cohort_data(id_mycorrhizal_scav_C_res,      diag, c(1:M), c(1:M)%scav_myc_C_reservoir,               weight=c(1:M)%nindivs, op=OP_SUM)
  call send_cohort_data(id_mycorrhizal_mine_N_res,      diag, c(1:M), c(1:M)%mine_myc_N_reservoir,               weight=c(1:M)%nindivs, op=OP_SUM)
  call send_cohort_data(id_mycorrhizal_mine_C_res,      diag, c(1:M), c(1:M)%mine_myc_C_reservoir,               weight=c(1:M)%nindivs, op=OP_SUM)
  call send_cohort_data(id_Nfix_N_res,                  diag, c(1:M), c(1:M)%N_fixer_N_reservoir,                weight=c(1:M)%nindivs, op=OP_SUM)
  call send_cohort_data(id_Nfix_C_res,                  diag, c(1:M), c(1:M)%N_fixer_C_reservoir,                weight=c(1:M)%nindivs, op=OP_SUM)

  call send_cohort_data(id_N_fix_alloc_smoothed,        diag, c(1:M), c(1:M)%max_Nfix_allocation,                weight=c(1:M)%nindivs, op=OP_SUM)
  call send_cohort_data(id_myc_mine_alloc_smoothed,     diag, c(1:M), c(1:M)%max_mine_allocation,                weight=c(1:M)%nindivs, op=OP_SUM)
  call send_cohort_data(id_myc_scav_alloc_smoothed,     diag, c(1:M), c(1:M)%max_scav_allocation,                weight=c(1:M)%nindivs, op=OP_SUM)

end subroutine vegn_carbon_int_ppa


! ============================================================================
! updates cohort biomass pools, LAI, SAI, and height using accumulated
! carbon_gain and bwood_gain
subroutine vegn_growth (vegn, diag)
  type(vegn_tile_type), intent(inout) :: vegn
  type(diag_buff_type), intent(inout) :: diag

  ! ---- local vars
  type(vegn_cohort_type), pointer :: cc    ! current cohort
  real :: cmass0,cmass1 ! for debug only
  integer :: i, N
  real, dimension(vegn%n_cohorts) :: &
     wood_prod, leaf_root_gr, sw_seed_gr, deltaDBH
  real :: f, borrowed

  if (is_watch_point()) then
     cmass0 = vegn_tile_carbon(vegn)
  endif

 ! write(*,*) 'counting cohorts: ', vegn%n_cohorts
  do i = 1, vegn%n_cohorts
     cc => vegn%cohorts(i)

     if (do_ppa) then
        if(is_watch_point()) then
           write(*,*)"############## in vegn_growth before allocation PPA #################"
           __DEBUG2__(cc%nsc, cc%carbon_gain)
        endif
        if (cc%nsc+cc%carbon_gain >= 0) then
           cc%nsc=cc%nsc+cc%carbon_gain
        else
           ! slm, ens 2017-08-02
           ! if loss of carbon exceeds current nsc, we borrow the carbon from living tissues
           ! proportionally to the existing amount.
           f = abs(cc%nsc+cc%carbon_gain)/max(abs(cc%nsc+cc%carbon_gain),cc%bl+cc%br+cc%blv+cc%bsw)
           borrowed = f*(cc%bl+cc%br+cc%blv+cc%bsw)
           cc%bl  = (1-f)*cc%bl
           cc%br  = (1-f)*cc%br
           cc%blv = (1-f)*cc%blv
           cc%bsw = (1-f)*cc%bsw
           ! brsw is part of bsw, and should also be reduced to stay within limits [0,bsw]
           cc%brsw = (1-f)*cc%brsw
           ! nsc should become zero, unless there is not enough biomass in all living
           ! tissues to compensate carbon loss
           cc%nsc = cc%nsc+cc%carbon_gain+borrowed
        endif
        call biomass_allocation_ppa(cc,wood_prod(i),leaf_root_gr(i),sw_seed_gr(i),deltaDBH(i))
        if(is_watch_point()) then
           write(*,*)"############## in vegn_growth after allocation PPA #################"
           __DEBUG2__(cc%nsc, cc%carbon_gain)
           __DEBUG5__(cc%bl, cc%blv, cc%br, cc%bsw, cc%bwood)
        endif
        cc%carbon_gain = 0.0
     else
        cc%bwood   = cc%bwood   + cc%bwood_gain
        cc%bliving = cc%bliving + cc%carbon_gain

        if(cc%bliving < 0) then
           cc%bwood    = cc%bwood+cc%bliving
           cc%bliving  = 0
           if (cc%bwood < 0) &
                cc%bwood = 0 ! in principle, that's not conserving carbon
        endif

        call update_biomass_pools(cc)

        ! reset carbon accumulation terms
        cc%carbon_gain = 0
        cc%carbon_loss = 0
        cc%bwood_gain  = 0

        wood_prod(i) = -100.0 ! missing value
     endif

     if (cc%status == LEAF_ON) then
        cc%leaf_age = cc%leaf_age + 1.0
        ! limit the maximum leaf age by the leaf time span (reciprocal of leaf
        ! turnover rate alpha) for given species. alpha is in 1/year, factor of
        ! 365 converts the result to days.
        if (spdata(cc%species)%alpha_leaf > 0) &
             cc%leaf_age = min(cc%leaf_age,365.0/spdata(cc%species)%alpha_leaf)

     endif
  end do
  N = vegn%n_cohorts
  call send_cohort_data(id_wood_prod,diag,vegn%cohorts(1:N),wood_prod(:), weight=vegn%cohorts(1:N)%nindivs, op=OP_SUM)
  call send_cohort_data(id_leaf_root_gr,diag,vegn%cohorts(1:N),leaf_root_gr(:), weight=vegn%cohorts(1:N)%nindivs, op=OP_SUM)
  call send_cohort_data(id_sw_seed_gr,diag,vegn%cohorts(1:N),sw_seed_gr(:), weight=vegn%cohorts(1:N)%nindivs, op=OP_SUM)
  ! conversion of DBH_growth assumes that vegn_growth is called daily
  call send_cohort_data(id_DBH_growth,diag,vegn%cohorts(1:N),deltaDBH(:)*days_per_year, weight=vegn%cohorts(1:N)%nindivs, op=OP_AVERAGE)
end subroutine vegn_growth

!========================================================================
! Starvation due to low NSC
subroutine vegn_starvation_ppa (vegn, soil)
  type(vegn_tile_type), intent(inout) :: vegn
  type(soil_tile_type), intent(inout) :: soil

  ! ---- local vars
  real :: deathrate ! mortality rate, 1/year
  real :: deadtrees ! number of trees that died over the time step
  integer :: i, k
  type(vegn_cohort_type), pointer :: cc(:) ! array to hold new cohorts
  ! accumulators of total input to root litter and soil carbon
  real, dimension(N_C_TYPES) :: &
     leaf_litt_C, leaf_litt_N, & ! fine surface litter per tile, kg/m2
     wood_litt_C, wood_litt_N    ! coarse surface litter per tile, kg/m2
  real :: root_litt_C(num_l, N_C_TYPES) ! root litter per soil layer, kgC/m2
  real :: root_litt_N(num_l, N_C_TYPES) ! root litter per soil layer, kgN/m2

  leaf_litt_C = 0 ; wood_litt_C = 0; root_litt_C = 0
  leaf_litt_N = 0 ; wood_litt_N = 0; root_litt_N = 0
  do i = 1, vegn%n_cohorts
     associate ( cc => vegn%cohorts(i)   , &
                 sp => spdata(vegn%cohorts(i)%species)  )  ! F2003

    ! Mortality due to starvation
    if (cc%bsw<0 .or. cc%nsc < 0.01*cc%bl_max) then
       deathrate = 1.0

       deadtrees = min(cc%nindivs*deathrate,cc%nindivs) ! individuals / m2
       ! kill starved plants and add dead C from leaf and root pools to soil carbon
       call kill_plants_ppa(cc, vegn, deadtrees, 0.0, &
                            leaf_litt_C, wood_litt_C, root_litt_C, &
                            leaf_litt_N, wood_litt_N, root_litt_N  )

       ! for budget tracking - temporary
       vegn%veg_out = deadtrees * (cc%bl+cc%br+cc%bsw+cc%blv+cc%bseed+cc%nsc+cc%bwood)
     endif

     end associate  ! F2003
  enddo

  ! calculate the number of remaining cohorts
  k = 0
  do i = 1, vegn%n_cohorts
     if (vegn%cohorts(i)%nindivs > 0.0) k=k+1
  enddo

!  if (k==0) call error_mesg('vegn_nat_mortality_ppa','All cohorts died',WARNING)
  if (is_watch_point()) then
     write(*,*)'###### vegn_starvation_ppa #######'
     __DEBUG1__(vegn%cohorts%nindivs)
     __DEBUG1__(k)
  endif

  ! remove cohorts that have zero individuals
  if (k < vegn%n_cohorts) then
     allocate(cc(max(k,1)))
     k=0
     do i = 1,vegn%n_cohorts
        if (vegn%cohorts(i)%nindivs > 0.0) then
           k=k+1
           cc(k) = vegn%cohorts(i)
        endif
     enddo

     if (k==0) then
        ! Most of the code assumes that there is at least one cohort present.
        ! So if all cohorts die, preserve single cohort with zero individuals.
        ! It probably doesn't matter which, but let's pick the shortest here.
        cc(1) = vegn%cohorts(vegn%n_cohorts)
        cc(1)%nindivs = 0.0
        vegn%n_cohorts = 1
     else
        vegn%n_cohorts = k
     endif
     deallocate (vegn%cohorts)
     vegn%cohorts=>cc
  endif

  ! add litter accumulated over the cohorts
  call add_soil_carbon(soil, vegn, leaf_litt_C, wood_litt_C, root_litt_C, &
                                   leaf_litt_N, wood_litt_N, root_litt_N  )

!  write(*,*)'vegn_starvation_ppa n_cohorts after: ', vegn%n_cohorts
end subroutine vegn_starvation_ppa



! ==============================================================================
! updates cohort vegetation structure, biomass pools, LAI, SAI, and height spending nsc
subroutine biomass_allocation_ppa(cc, wood_prod,leaf_root_gr,sw_seed_gr,deltaDBH)
  type(vegn_cohort_type), intent(inout) :: cc
  real, intent(out) :: wood_prod ! wood production, kgC/year per individual, diagnostic output
  real, intent(out) :: leaf_root_gr! allocation to leaf and fine root, kgC/year
  real, intent(out) :: sw_seed_gr! allocation to sapwood and seed, kgC/year
  real, intent(out) :: deltaDBH ! tendency of breast height diameter, m

  real :: CSAtot ! total cross section area, m2
  real :: CSAsw  ! Sapwood cross sectional area, m2
  real :: CSAwd  ! Heartwood cross sectional area, m2
  real :: DBHwd  ! diameter of heartwood at breast height, m
  real :: G_LFR  ! amount of carbon spent on leaf and root growth
  real :: G_WF   ! amount of carbon spent on new sapwood growth and seeds
  real :: G_WF_max ! max spending rate
  real :: deltaSeed
  real :: deltaBL, deltaBR ! tendencies of leaf and root biomass, kgC/individual
  real :: deltaNL, deltaNR ! tendencies of leaf and root N, kg/individual
  real :: deltaBSW ! tendency of sapwood biomass, kgC/individual
  real :: deltaBwood ! tendency of wood biomass, kgC/individual
  real :: deltaCA ! tendency of crown area, m2/individual
  real :: deltaHeight ! tendency of vegetation height
  real :: deltaCSAsw ! tendency of sapwood area
  real :: BL_c, BL_u ! canopy and understory target leaf biomasses, kgC/individual
  real :: NSCtarget ! target NSC storage
  real :: N_storage_target
  real :: delta_bsw_branch, delta_nsw_branch
  real :: br_max_Nstress
  real :: b0, b1, n0, n1 ! biomasses for conservation checking


  real, parameter :: DBHtp = 0.8 ! m

  !ens
  real, parameter :: fsf = 0.2 ! in Weng et al 205, page 2677 units are 0.2 day
  real, parameter :: deltaDBH_max = 0.01 ! max growth rate of the trunk, meters per day
  real  :: nsctmp ! temporal nsc budget to avoid biomass loss, KgC/individual

  associate (sp => spdata(cc%species)) ! F2003

  if (do_check_conservation) then
      b0 = plant_C(cc)
      n0 = plant_N(cc)
  endif
  if(is_watch_point()) then
     write(*,*)'#### biomass_allocation_ppa input ####'
     __DEBUG5__(cc%bl, cc%br, cc%bsw, cc%bseed, cc%bwood)
     __DEBUG2__(cc%nsc, cc%npp_previous_day)
  endif

  ! Stress increases as stored N declines relative to total biomass N demand
  ! N stress is calculated based on "potential pools" without N limitation
  ! Elena suggests using 2*(root N + leaf N) as storage target
  ! bliving is already increased after cgain


   NSCtarget = 4.0*cc%bl_max
   N_storage_target = NSCtarget/sp%leaf_live_c2n

   if(N_limits_live_biomass) call check_var_range(cc%stored_N,0.0,HUGE(1.0),'biomass_allocation_ppa','cc%stored_N',FATAL)

  ! 20170617: calculate N stress (function of stored N) here
  !           N traget possibly in relation with NSC target
  if (soil_carbon_option==SOILC_CORPSE_N) then
     if (cc%status == LEAF_ON) then
        cc%nitrogen_stress = max(min_N_stress,(N_storage_target-cc%stored_N)/N_storage_target)
     else
        ! If leaves are off, our N stress would be lower, since we re-translocate some of leaf+root N
        ! but plant should be saving enough to rebuild leaf+root next time.
        ! In this case, we should adjust the calculation:
        cc%nitrogen_stress = max(min_N_stress,(N_storage_target+(cc%bl_max/sp%leaf_live_c2n+cc%br_max/sp%froot_live_c2n)*0.5-cc%stored_N)/N_storage_target)
        ! __DEBUG4__(cc%bl_max/sp%leaf_live_c2n,cc%br_max/sp%froot_live_c2n,N_storage_target,cc%stored_N)
     endif
  else
     cc%nitrogen_stress = 0.0
  endif

  ! TODO: what if carbon_gain is not 0, but leaves are OFF (marginal case? or
  ! typical in lm3?)
  ! set init values for diag output, to be returned when the actual calulations are bypassed:
  wood_prod = 0.0 ; leaf_root_gr = 0.0 ; sw_seed_gr = 0.0 ; deltaDBH = 0.0
  if (cc%status == LEAF_ON) then
     if(is_watch_point()) then
        write(*,*)'########### biomass_allocation_ppa input ###########'
        call dpri('species:',sp%name)
        __DEBUG2__(cc%bl_max, cc%br_max)
        __DEBUG2__(cc%carbon_gain, cc%nsc)
        __DEBUG5__(cc%bl, cc%br, cc%bsw, cc%bseed, cc%bwood)
        __DEBUG3__(cc%dbh, cc%height, cc%crownarea)
        write(*,*)'########### end of biomass_allocation_ppa input ###########'
     endif
     ! update targets based on the nitrogen stress
     ! N_stress_root_factor is something like 0.05 so br_max increases 25% when N stress is 0.5, i.e. stored N is half of target value
     ! N stress or this root factor should not go above some max value, so all biomass doesn't end up in roots in some weird situation
     ! This assumes carbon to build more roots comes from sapwood growth, not leaf growth. But we can make this flexible or have different options
     br_max_Nstress = cc%br_max * (1+cc%nitrogen_stress*sp%N_stress_root_factor)

     !ens first replace lost branch sapwood

     ! in principle re-building wood and sapwood of branches should be taken into account together
     ! since they are the same process and should be limited by the upper limit of
     ! N storage spending rate together.
     call check_var_range(cc%bsw, 0.0,HUGE(1.0),'biomass_allocation_ppa #1', 'cc%bsw', FATAL)
     call check_var_range(cc%brsw,0.0,cc%bsw,   'biomass_allocation_ppa #1', 'cc%brsw',FATAL)

     ! 02/07/17
     ! should this delta_bsw be updated after updating dbh and height?
     ! check if daily branch increase is nt too abrupt
     select case (sp%allomt)
     case (ALLOM_EW,ALLOM_EW1)
        delta_bsw_branch = sp%branch_wood_frac * sp%alphaBM * sp%rho_wood * cc%DBH**sp%thetaBM - cc%brsw
     case (ALLOM_HML)
        delta_bsw_branch = sp%branch_wood_frac * sp%alphaBM * sp%rho_wood * cc%DBH**2 * cc%height - cc%brsw
     end select
     delta_bsw_branch = max(min(delta_bsw_branch,0.1*cc%nsc/(1+sp%GROWTH_RESP)),0.0)
     if (soil_carbon_option==SOILC_CORPSE_N) then
        delta_nsw_branch = delta_bsw_branch/sp%sapwood_c2n
     else
        delta_nsw_branch = 0.0
     endif
     if (delta_nsw_branch>0.1*cc%stored_N.AND.N_limits_live_biomass) then
        delta_bsw_branch = delta_bsw_branch*0.1*cc%stored_N/delta_nsw_branch
        delta_nsw_branch = 0.1*cc%stored_N
     endif
     cc%brsw = cc%brsw + delta_bsw_branch
     cc%bsw  = cc%bsw  + delta_bsw_branch
     if (soil_carbon_option==SOILC_CORPSE_N) then
        cc%sapwood_N = cc%sapwood_N + delta_nsw_branch
        cc%stored_N  = cc%stored_N  - delta_nsw_branch
     endif

     call check_var_range(cc%bsw, 0.0,HUGE(1.0),'biomass_allocation_ppa #2', 'cc%bsw',FATAL)
     call check_var_range(cc%brsw,0.0,cc%bsw,   'biomass_allocation_ppa #2', 'cc%brsw',FATAL)

     ! make sure that if nsc is not spent because N is limited, sapwood is limited too

     ! ens 02/14/17
     ! update seed and sapwood biomass pools with the new growth (branches were updated above)
     select case (sp%lifeform)
     case (FORM_WOODY)
        NSCtarget = sp%NSC2targetbl*cc%bl_max
        G_WF=0.0
        if (cc%nsc > NSCtarget) then ! ens change this
           G_WF = max (0.0, fsf*(cc%nsc - NSCtarget)/(1+sp%GROWTH_RESP))
        endif

        if (cohort_makes_seeds(cc,G_WF)) then
           deltaSeed = sp%v_seed * G_WF
        else
           deltaSeed = 0.0
        endif
        deltaBSW = G_WF - deltaSeed

        if (deltaBSW/sp%sapwood_c2n>0.1*cc%stored_N.and.N_limits_live_biomass) then
            deltaBSW = 0.1*cc%stored_N*sp%sapwood_c2n
            G_WF = deltaBSW + deltaSeed
        endif

        ! calculate carbon spent on growth of leaves and roots
        G_LFR = max(0.0, min(cc%bl_max+br_max_Nstress-cc%bl-cc%br,  &
                             0.1*cc%nsc/(1+sp%GROWTH_RESP))) ! don't allow more than 0.1/(1+GROWTH_RESP) of nsc per day to spend

        ! and distribute it between roots and leaves
        deltaBL  = min(G_LFR, max(0.0, &
             (G_LFR*cc%bl_max + cc%bl_max*cc%br - br_max_Nstress*cc%bl)/(cc%bl_max + br_max_Nstress) &
             ))
        deltaBR  = G_LFR - deltaBL

     case (FORM_GRASS) ! isa 20170705
        ! 20170724 - new scheme
        NSCtarget = sp%NSC2targetbl*cc%bl_max
        G_WF = max(0.0, fsf*(nsctmp - NSCtarget)/(1+sp%GROWTH_RESP))
        ! note that it is only for HML allometry now
        G_WF_max = deltaDBH_max/((sp%gammaHT+cc%DBH**sp%thetaHT)**2/(sp%rho_wood * sp%alphaHT * sp%alphaBM * &
                          (cc%DBH**(1.+sp%thetaHT)*(2.*(sp%gammaHT+cc%DBH**sp%thetaHT)+sp%gammaHT*sp%thetaHT))))
        G_WF = min(G_WF, G_WF_max)
        if (cc%layer == 1 .AND. cc%age > sp%maturalage) then
           deltaSeed=      sp%v_seed * G_WF
           deltaBSW = (1.0-sp%v_seed)* G_WF
        else
           deltaSeed= 0.0
           deltaBSW = G_WF
        endif

        if (deltaBSW/sp%sapwood_c2n>0.1*cc%stored_N.and.N_limits_live_biomass) then
           deltaBSW = 0.1*cc%stored_N*sp%sapwood_c2n
           G_WF = deltaBSW+deltaSeed
        endif

        ! Allocation to leaves and fine root growth
        nsctmp = cc%nsc - G_WF * (1+sp%GROWTH_RESP)
        G_LFR = max(0.0, min( cc%bl_max - cc%bl + max( 0.0, cc%br_max - cc%br ),  &
                              0.1*nsctmp/(1+sp%GROWTH_RESP)))
        ! don't allow more than 0.1/(1+GROWTH_RESP) of nsc per day to spend

        if ( cc%br > cc%br_max ) then      ! roots larger than target
            deltaBL = G_LFR
        else                               ! standard scheme
            deltaBL = min(G_LFR, max(0.0, &
                         (G_LFR*cc%bl_max + cc%bl_max*cc%br - cc%br_max*cc%bl)/(cc%bl_max + cc%br_max) ))
        end if
        deltaBR  = G_LFR - deltaBL
     end select
     ! update leaf and root C tendencies based on N limitations, and calculate
     ! leaf and root N tendencies
     call update_LR_tendencies(cc, deltaBL, deltaBR, deltaNL, deltaNR)
     G_LFR = deltaBL+deltaBR

     ! update biomass pools due to growth
     cc%bl     = cc%bl    + deltaBL;   cc%leaf_N = cc%leaf_N + deltaNL ! updated in vegn_int_ppa
     cc%br     = cc%br    + deltaBR;   cc%root_N = cc%root_N + deltaNR
     cc%bsw    = cc%bsw   + deltaBSW;  if (soil_carbon_option==SOILC_CORPSE_N) cc%sapwood_N = cc%sapwood_N + deltaBSW/sp%sapwood_c2n
     cc%bseed  = cc%bseed + deltaSeed;
     if (soil_carbon_option==SOILC_CORPSE_N) &
           cc%seed_N = cc%seed_N + deltaSeed/sp%seed_c2n ! TODO: calculate seed_c2n as derived quantity in vegn_data initialization
     cc%nsc    = cc%nsc - (deltaBL + deltaBR + deltaSeed + deltaBSW + delta_bsw_branch)*(1+sp%GROWTH_RESP)
     cc%stored_N = cc%stored_N - deltaNL - deltaNR

     if(cc%stored_N - deltaBSW/sp%sapwood_c2n - deltaSeed/sp%seed_c2n<0.and.soil_carbon_option==SOILC_CORPSE_N) then
         __DEBUG3__(cc%stored_N,deltaBSW/sp%sapwood_c2n,deltaSeed/sp%seed_c2n)
     endif
     if (soil_carbon_option==SOILC_CORPSE_N) cc%stored_N = cc%stored_N - deltaBSW/sp%sapwood_c2n - deltaSeed/sp%seed_c2n

     wood_prod = deltaBSW*days_per_year ! conversion from kgC/day to kgC/year
     ! compute daily respiration fluxes
     leaf_root_gr = G_LFR*days_per_year ! conversion from kgC/day to kgC/year
     sw_seed_gr = (G_WF+delta_bsw_branch )*sp%GROWTH_RESP*days_per_year ! conversion from kgC/day to kgC/year

     call check_var_range(cc%nsc,-1e-16,HUGE(1.0),'biomass_allocation_ppa','cc%nsc', WARNING)

     !ens --compute daily growth to compute respiration, apply it next day, use npp_previous day variable, units kg C/(m2 *year)
     cc%growth_previous_day = cc%growth_previous_day+(max(0., G_LFR+G_WF)+delta_bsw_branch)*sp%GROWTH_RESP ! this is for growth respiration to come from nsc

     select case (sp%allomt)
     case (ALLOM_EW, ALLOM_EW1)
        ! calculate tendency of breast height diameter given increase of bsw
        deltaDBH     = deltaBSW / (sp%thetaBM * sp%alphaBM * sp%rho_wood * cc%DBH**(sp%thetaBM-1))
        deltaHeight  = sp%thetaHT * sp%alphaHT * cc%DBH**(sp%thetaHT-1) * deltaDBH
        if(sp%allomt == ALLOM_EW1) then
            deltaCA  = sp%thetaCA * sp%alphaCA *     &
                  (1./(exp(15.*(cc%DBH-DBHtp))+1.))  *     &
                  cc%DBH**(sp%thetaCA-1) * deltaDBH
        else
            deltaCA  = sp%thetaCA * sp%alphaCA * cc%DBH**(sp%thetaCA-1) * deltaDBH
        endif
     case (ALLOM_HML)
        deltaDBH     = deltaBSW*(sp%gammaHT+cc%DBH**sp%thetaHT)**2/(sp%rho_wood * sp%alphaHT * sp%alphaBM * &
                       (cc%DBH**(1.+sp%thetaHT)*(2.*(sp%gammaHT+cc%DBH**sp%thetaHT)+sp%gammaHT*sp%thetaHT)))
!         deltaHeight  = sp%alphaHT*(cc%DBH+deltaDBH)**sp%thetaHT/(sp%gammaHT+(cc%DBH+deltaDBH)**sp%thetaHT) &
!                        - cc%height
        deltaHeight  = deltaDBH* sp%alphaHT *sp%gammaHT*sp%thetaHT/(cc%DBH**(1.-sp%thetaHT) * &
                       (sp%gammaHT + cc%DBH**sp%thetaHT)**2)
        deltaCA      = deltaDBH * sp%alphaCA * sp%thetaCA * cc%DBH**(sp%thetaCA-1.)
     case default
        call error_mesg('biomass_allocation_ppa','Unknown allometry type. This should never happen.', FATAL)
     end select

     cc%DBH       = cc%DBH       + deltaDBH
     cc%height    = cc%height    + deltaHeight
     cc%crownarea = cc%crownarea + deltaCA

     ! calculate DBH, BLmax, BRmax, BSWmax using allometric relationships
     ! Weng 2012-01-31 update_bio_living_fraction
     ! slm 20160523: max biomass of sapwood BSWmax is calculated from allometric
     ! relationships as the B(DBH) - B(DBH of heartwood). DBH of heartwood is
     ! calculated using the allometric relationship for sapwood cross-section
     select case(sp%allomt)
     case (ALLOM_EW,ALLOM_EW1)
        CSAsw    = sp%alphaCSASW * cc%DBH**sp%thetaCSASW ! sapwood cross-section
        CSAtot   = PI * (cc%DBH/2.0)**2 ! total trunk cross-section
        CSAwd    = max(0.0, CSAtot - CSAsw) ! cross-section of heartwood
        DBHwd    = 2*sqrt(CSAwd/PI) ! DBH of heartwood
        ! before only BSWmax, 4 lines above from Isa
        cc%bsw_max = sp%alphaBM * sp%rho_wood * (cc%DBH**sp%thetaBM - DBHwd**sp%thetaBM)
     case (ALLOM_HML)
        CSAsw    = sp%phiCSA * cc%DBH**(sp%thetaCA + sp%thetaHT) / &
                   (sp%gammaHT + cc%DBH** sp%thetaHT)
        CSAtot   = PI * (cc%DBH/2.0)**2 ! total trunk cross-section
        CSAwd    = max(0.0, CSAtot - CSAsw) ! cross-section of heartwood
        DBHwd    = 2*sqrt(CSAwd/PI) ! DBH of heartwood
        cc%bsw_max = sp%alphaBM * sp%rho_wood * cc%height * (cc%DBH**2 - DBHwd**2)
     end select
     ! isa 20170705 - grasses don't form heartwood
     if (sp%lifeform == FORM_GRASS) then
       cc%bsw_max = cc%bsw
       CSAsw = PI * cc%DBH * cc%DBH / 4.0 ! trunk cross-sectional area = sapwood area
     endif

     ! isa and es 201701 - CSAsw different for ALLOM_HML - include into previous case?
     select case(sp%allomt)
     case (ALLOM_EW,ALLOM_EW1)
        deltaCSAsw = CSAsw - sp%alphaCSASW * (cc%DBH - deltaDBH)**sp%thetaCSASW
     case (ALLOM_HML)
        deltaCSAsw = CSAsw - sp%phiCSA * (cc%DBH - deltaDBH)**(sp%thetaCA + sp%thetaHT) / (sp%gammaHT + (cc%DBH - deltaDBH)**sp%thetaHT)
     end select
     ! isa 20170705 - grasses don't form heartwood
     if (sp%lifeform == FORM_GRASS) then
        deltaCSAsw = CSAsw - (PI * (cc%DBH - deltaDBH) * (cc%DBH - deltaDBH) / 4.0 )
     endif
     cc%Kxa = (cc%Kxa*CSAsw + sp%Kxam*deltaCSAsw)/(CSAsw + deltaCSAsw)

     call check_var_range(cc%bsw,  0.0,HUGE(1.0),'biomass_allocation_ppa #3', 'cc%bsw', FATAL)
     call check_var_range(cc%brsw, 0.0,cc%bsw,   'biomass_allocation_ppa #3', 'cc%brsw',FATAL)

     ! ens 02/14/17 replace wood allocation function
     ! slm 06/28/17 Why are we retiring sapwood to wood only if the leaves are displayed?
     !              shouldn't it also happen in leaf-off season?
!     deltaBwood = max(cc%bsw - cc%brsw - (1.0 - ( cc%brsw / cc%bsw ) )*BSWmax, 0.0)
     if (cc%bsw>0) then
        deltaBwood = cc%bsw - cc%brsw - (1.0 - ( cc%brsw / cc%bsw ) )*cc%bsw_max
        ! perhaps it should be min(deltaBwood, cc%bsw-cc%brsw) ?
        deltaBwood = min(deltaBwood, cc%bsw) ! cannot retire more sapwood than we have
        deltaBwood = max(deltaBwood, 0.0)    ! conversion of wood to sapwood is prohibited
     else
        deltaBwood = 0.0
     endif
     cc%bwood   = cc%bwood + deltaBwood
     cc%bsw     = cc%bsw   - deltaBwood
     if (soil_carbon_option==SOILC_CORPSE_N) then
       cc%sapwood_N = cc%sapwood_N - deltaBwood/sp%sapwood_c2n
       cc%wood_N    = cc%wood_N    + deltaBwood/sp%wood_c2n
       cc%stored_N  = cc%stored_N  + deltaBwood/sp%sapwood_c2n - deltaBwood/sp%wood_c2n
     endif

     call check_var_range(cc%bsw,  0.0,HUGE(1.0),'biomass_allocation_ppa #4', 'cc%bsw',FATAL)
     call check_var_range(cc%brsw, 0.0,cc%bsw,   'biomass_allocation_ppa #4', 'cc%brsw',FATAL)

     if (cc%An_newleaf_daily > 0 ) then
         cc%laimax = min(cc%laimax + sp%newleaf_layer,laimax_ceiling)
     else if (cc%An_newleaf_daily < 0) then
         cc%laimax = max(cc%laimax - sp%newleaf_layer,laimax_floor)
     else
         ! do nothing if the derivative is zero
     endif
     cc%An_newleaf_daily = 0.0
     ! update bl_max and br_max daily
     ! slm: why are we updating topyear only when the leaves are displayed? The paper
     !      never mentions this fact (see eq. A6).
     BL_u = sp%LMA * sp%LAImax * cc%crownarea * (1.0-sp%internal_gap_frac) * understory_lai_factor
     BL_c = sp%LMA * sp%LAImax * cc%crownarea * (1.0-sp%internal_gap_frac)
     if (cc%layer > 1 .and. cc%firstlayer == 0) then ! changed back, Weng 2014-01-23
        cc%topyear = 0.0
        cc%bl_max = BL_u
     else
        if(cc%layer == 1)cc%topyear = cc%topyear + 1.0/365.0  ! daily step
        if(cc%layer>1)cc%firstlayer = 0 ! Just for the first year, those who were
                  ! pushed to understory have the characteristics of canopy trees
        if(sp%lifeform == FORM_GRASS) then
           cc%bl_max = BL_c
        else
           cc%bl_max = BL_u + min(cc%topyear/5.0,1.0)*(BL_c - BL_u)
        endif
     endif
     ! in case of light saber override bl_max, but the keep the value of firstlayer
     ! calculated above
     if (use_light_saber) then
        cc%bl_max = sp%LMA * cc%laimax * cc%crownarea * (1.0-sp%internal_gap_frac)
     endif
     cc%br_max = sp%phiRL*cc%bl_max/(sp%LMA*sp%SRA)

     if(is_watch_point()) then
        write(*,*)'########### biomass_allocation_ppa output ###########'
        __DEBUG5__(deltaBL, deltaBR, deltaBSW, deltaSeed, deltaBwood)
        __DEBUG5__(cc%bl,   cc%br,   cc%bsw,   cc%bseed,  cc%bwood)
        __DEBUG3__(deltaDBH, deltaHeight, deltaCA)
        __DEBUG3__(cc%dbh,   cc%height,   cc%crownarea)
        __DEBUG1__(cc%nsc)
        write(*,*)'########### end of biomass_allocation_ppa output ###########'
     endif
  else  ! cc%status == LEAF_OFF
    ! should some nsc go into sapwood and wood or seed?
  endif ! cc%status == LEAF_ON

  cc%total_N    = cc%stored_N+cc%leaf_N+cc%wood_N+cc%root_N+cc%sapwood_N+cc%seed_N

  ! calculate spending rate of growth respiration, to distribute it uniformly
  ! in time over the next day:
  cc%growth_previous_day_tmp = max(0.0,cc%growth_previous_day)*365.0
  ! factor 365.0 converts the rate of growth respiration release to atmosphere
  ! from kgC/day (frequency of this subroutine calls) to kgC/year, the units we
  ! use for other vegetation fluxes

  if(N_limits_live_biomass) call check_var_range(cc%stored_N,0.0,HUGE(1.0),'biomass_allocation_ppa','cc%stored_N',FATAL)

  if (do_check_conservation) then
     b1 = plant_C(cc); n1=plant_N(cc)
     call check_conservation ('biomass_allocation_ppa','carbon', b0, b1, carbon_cons_tol, severity=FATAL)
     call check_conservation ('biomass_allocation_ppa','nitrogen', n0, n1, nitrogen_cons_tol, severity=FATAL)
  endif
  end associate ! F2003
end subroutine biomass_allocation_ppa

! update leaf and root biomass tendencies based on plant nitrogen state, and calculate
! leaf and root nitrogen tendencies
subroutine update_LR_tendencies(cc, deltaBL, deltaBR, deltaNL, deltaNR)
  type(vegn_cohort_type), intent(in) :: cc
  real, intent(inout) :: deltaBL, deltaBR
  real, intent(out)   :: deltaNL, deltaNR

  real :: f

  associate(sp=>spdata(cc%species))
  if (soil_carbon_option==SOILC_CORPSE_N) then
     deltaNL = deltaBL/sp%leaf_live_c2n
     deltaNR = deltaBR/sp%froot_live_c2n
     if (deltaNL+deltaNR > 0.1*cc%stored_N.and.N_limits_live_biomass) then
        ! BNS: make this the minimum of 10% NSC depletion OR 10% stored N depletion
        ! Do this every time we do this limitation.
        ! This does not prioritize spending on leaves and roots
        f = 0.1*cc%stored_N/(deltaNL+deltaNR)
        deltaBR = deltaBR*f; deltaBL = deltaBL*f
        deltaNR = deltaNR*f; deltaNL = deltaNL*f
     endif
  else
     deltaNL=0.0
     deltaNR=0.0
  endif
  end associate
end subroutine

! ============================================================================
! calculated thermal inhibition factor depending on temperature
function thermal_inhibition(T) result(tfs); real tfs
  real, intent(in) :: T ! temperature, degK

  tfs = exp(3000.0*(1.0/288.16-1.0/T));
  tfs = tfs / ( &
              (1.0+exp(0.4*(5.0-T+273.16)))* &
              (1.0+exp(0.4*(T - 273.16-45.0)))&
              )
end function


! ============================================================================
subroutine plant_respiration(cc, tsoil, resp, r_leaf, r_root, r_stem)
  type(vegn_cohort_type), intent(in) :: cc
  real, intent(in) :: tsoil
  real, intent(out) :: resp ! total respiration
  real, intent(out) :: r_leaf ! leaf respiration
  real, intent(out) :: r_root ! fine root respiration
  real, intent(out) :: r_stem ! stem respiration

  real :: tf,tfs ! thermal inhibition factors for above- and below-ground biomass
  real :: r_vleaf ! virtual leaf respiration (perhaps need to be reported, too?)
  real :: Acambium  ! cambium area, m2

  integer :: sp ! shorthand for cohort species
  sp = cc%species

  tf  = thermal_inhibition(cc%Tv)
  tfs = thermal_inhibition(tsoil)

  r_leaf = -mol_C*cc%An_cl*cc%leafarea;
  r_vleaf = spdata(sp)%beta_vleaf * cc%blv*tf;
  if (do_ppa) then
     ! Stem auto-respiration is proportional to cambium area, not sapwood biomass
     Acambium = PI * cc%DBH * cc%height * 1.2
     r_stem   = spdata(sp)%beta_sapwood * Acambium * tf
  else
     r_stem   = spdata(sp)%beta_sapwood * cc%bsw * tf
  endif
  r_root  = spdata(sp)%beta_root * cc%br*tfs;

  resp = r_leaf + r_vleaf + r_stem + r_root
  ! BNS: Should resp be limited so it's not greater than nsc? Otherwise nsc can go less than zero if gpp is zero
end subroutine plant_respiration


! =============================================================================
subroutine vegn_phenology_lm3(vegn, soil)
  type(vegn_tile_type), intent(inout) :: vegn
  type(soil_tile_type), intent(inout) :: soil

  ! ---- local vars
  real :: leaf_litter_C,root_litter_C,leaf_litter_N,root_litter_N
  real :: theta_crit; ! critical ratio of average soil water to sat. water
  real :: psi_stress_crit ! critical soil-water-stress index
  real :: wilt ! ratio of wilting to saturated water content
  real, dimension(N_C_TYPES) :: &
      leaf_litt_C, leaf_litt_N ! fine surface litter per tile, kgC/m2 and kgN/m2
  real, dimension(num_l,n_c_types) :: &
      root_litt_C, root_litt_N ! root litter per soil layer, kgC/m2 and kgN/m2
  real :: profile(num_l) ! storage for vertical profile of root litter
  integer :: i, l

  wilt = soil%w_wilt(1)/soil%pars%vwc_sat
  vegn%litter = 0

  leaf_litt_C = 0 ; root_litt_C = 0 ; leaf_litt_N = 0 ; root_litt_N = 0
  do i = 1,vegn%n_cohorts
     associate ( cc => vegn%cohorts(i),   &
                 sp => spdata(vegn%cohorts(i)%species) )

     if(is_watch_point())then
        write(*,*)'####### vegn_phenology #######'
        __DEBUG4__(vegn%theta_av_phen, wilt, sp%cnst_crit_phen, sp%fact_crit_phen)
        __DEBUG2__(vegn%psist_av, sp%psi_stress_crit_phen)
        __DEBUG1__(cc%species)
        __DEBUG2__(vegn%tc_av,sp%tc_crit)
     endif
     ! if drought-deciduous or cold-deciduous species
     ! temp=10 degrees C gives good growing season pattern
     ! spp=0 is c3 grass,1 c3 grass,2 deciduous, 3 evergreen
     ! assumption is that no difference between drought and cold deciduous
     cc%status = LEAF_ON; ! set status to indicate no leaf drop

     if(sp%phent == PHEN_DECIDUOUS) then ! deciduous species
        ! actually either fact_crit_phen or cnst_crit_phen is zero, enforced
        ! by logic in the vegn_data.F90
        theta_crit = sp%cnst_crit_phen &
              + wilt*sp%fact_crit_phen
        theta_crit = max(0.0,min(1.0, theta_crit))
        psi_stress_crit = sp%psi_stress_crit_phen
        if (      (psi_stress_crit <= 0. .and. vegn%theta_av_phen < theta_crit) &
             .or. (psi_stress_crit  > 0. .and. vegn%psist_av > psi_stress_crit) &
             .or. (vegn%tc_av < sp%tc_crit) ) then
           cc%status = LEAF_OFF; ! set status to indicate leaf drop
           cc%leaf_age = 0;

           leaf_litter_C = (1.0-l_fract)*cc%bl*cc%nindivs
           leaf_litt_C(:) = leaf_litt_C(:)+[sp%fsc_liv,1-sp%fsc_liv,0.0]*leaf_litter_C
           leaf_litter_N =  cc%leaf_N*cc%nindivs*(1.0-sp%leaf_retranslocation_frac)
           leaf_litt_N(:) = leaf_litt_N(:)+[sp%fsc_liv,1-sp%fsc_liv,0.0]*leaf_litter_N

           root_litter_C = (1.0-l_fract)*cc%br*cc%nindivs
           root_litter_N = cc%root_N*cc%nindivs*(1.0-sp%froot_retranslocation_frac)
           call cohort_root_litter_profile(cc, dz, profile)
           do l = 1, num_l
              root_litt_C(l,:) = root_litt_C(l,:) + profile(l)* &
                   [sp%fsc_froot, 1-sp%fsc_froot, 0.0]*root_litter_C
              root_litt_N(l,:) = root_litt_N(l,:) + profile(l)* &
                   [sp%fsc_froot, 1-sp%fsc_froot, 0.0]*root_litter_N
           enddo

           vegn%litter = vegn%litter + leaf_litter_C + root_litter_C
           vegn%veg_out = vegn%veg_out + leaf_litter_C + root_litter_C

           cc%blv = cc%blv + l_fract*(cc%bl+cc%br);
           cc%bl  = 0.0;
           cc%br  = 0.0;
           cc%lai = 0.0;

           cc%stored_N = cc%stored_N + (cc%leaf_N*sp%leaf_retranslocation_frac+cc%root_N*sp%froot_retranslocation_frac)
           cc%leaf_N = 0.0
           cc%root_N = 0.0

           ! update state
           cc%bliving = cc%blv + cc%br + cc%bl + cc%bsw;
        endif
     endif ! phenology type
     call update_biomass_pools(cc); ! BNS: Just do this every time to avoid things getting messed up
     end associate
  enddo

  ! add litter accumulated over the cohorts
  call add_soil_carbon(soil, vegn, leaf_litter_C=leaf_litt_C, root_litter_C=root_litt_C, &
                                   leaf_litter_N=leaf_litt_N, root_litter_N=root_litt_N  )

end subroutine vegn_phenology_lm3


! =============================================================================
! Added by Weng 2012-02-29
subroutine vegn_phenology_ppa(tile)
  type(land_tile_type), intent(inout) :: tile

  ! ---- local vars
  integer :: i,l
  real    :: leaf_litter_C, leaf_litter_N, leaf_litt_C(N_C_TYPES), leaf_litt_N(N_C_TYPES)
  real, dimension(num_l,n_c_types) :: &
             root_litt_C, root_litt_N ! root litter per soil layer, kgC/m2 and kgN/m2
  real    :: profile(num_l) ! storage for vertical profile of root litter
  real    :: root_litter_C, root_litter_N
  real    :: dead_leaves_C, dead_leaves_N
  real    :: dead_roots_C, dead_roots_N
  real    :: dead_stem_C, dead_stem_N ! isa 20170705
  real    :: dheat ! heat residual due to cohort merging

  ! rates of various tissues decay after leaf drop, per day
  real, parameter :: leaf_fall_rate = 0.075
  real, parameter :: root_mort_rate = 0.0
  real, parameter :: stem_mort_rate = 0.075 ! isa 20170705

  real :: wilt       ! vertically-average soil moisture content at wilting point
  real :: theta_crit ! critical soil moisture for drought-deciduous phenology
  logical :: drought ! drought indicator


  associate(vegn=>tile%vegn, soil=>tile%soil)

  vegn%litter = 0 ;
  leaf_litt_C(:) = 0.0 ; leaf_litt_N(:) = 0.0
  root_litt_C    = 0.0 ; root_litt_N    = 0.0
  do i = 1,vegn%n_cohorts
     associate ( cc => vegn%cohorts(i),   &
                 sp => spdata(vegn%cohorts(i)%species) )
     if (sp%phent==PHEN_EVERGREEN) then
        cc%status = LEAF_ON
        cycle ! do nothing further for evergreens
     endif

     ! if drought-deciduous or cold-deciduous species
     ! temp=10 degrees C gives good growing season pattern
     ! spp=0 is c3 grass,1 c3 grass,2 deciduous, 3 evergreen
     ! assumption is that no difference between drought and cold deciduous

     wilt = soil%w_wilt(1)/soil%pars%vwc_sat
     theta_crit = sp%cnst_crit_phen + wilt*sp%fact_crit_phen
     theta_crit = max(0.0,min(1.0, theta_crit))
     drought    = (sp%psi_stress_crit_phen <= 0 .and. vegn%theta_av_phen < theta_crit) &
             .or. (sp%psi_stress_crit_phen  > 0 .and. vegn%psist_av > sp%psi_stress_crit_phen)

     ! onset of phenology
     select case(cc%status)
     case (LEAF_OFF)
        if (cc%gdd > sp%gdd_crit        .and. &
            vegn%tc_pheno >= sp%tc_crit .and. &
            .not.drought) then
           cc%status = LEAF_ON
           ! isa 201707215 - update target biomass at the satart of the growing season
           !                 for grasses to avoid using previous year targets
           if (sp%lifeform == FORM_GRASS) then
              cc%bl_max = sp%LMA * sp%laimax * cc%crownarea * (1.0-sp%internal_gap_frac)
              cc%br_max = sp%phiRL*cc%bl_max/(sp%LMA*sp%SRA)
           endif
        endif
     case (LEAF_ON)
        ! cold-deciduous: turn leaves off and reset GDD
        if (vegn%tc_pheno < sp%tc_crit) then
           cc%status = LEAF_OFF
           cc%gdd = 0.0
        endif
        ! drought-deciduous: turn leaves off and DO NOT reset GDD
        if (drought) then
           cc%status = LEAF_OFF
           ! Do not reset GDD when leaves drop due to drought, because if we do then
           ! new leaves would not appear until critical GDD is accumulated again.
           ! This would not make sense in many warm droughty regions, e.g. in sub-tropics.
        endif
     end select
     ! slm: there is a potential strange corner case above: when leaves drop due to
     ! autumn drought in a cold region, GDD is not reset, and therefore GDD keeps
     ! accumulating from the previous year. This will result in too-early onset of
     ! phenology the next year, probably instantly after the drought ends.
     !
     ! On the other hand, resetting GDD every time vegn%tc_pheno<sp%tc_crit will not work,
     ! since it will result in GDD never accumulating in temperature range gdd_base_T<T<tc_crit.

     ! leaves falling at the end of a growing season
     if(cc%status==LEAF_OFF .AND. ( cc%bl>0 .OR. ( sp%lifeform==FORM_GRASS .AND. ( cc%bl>0 .OR. cc%bsw>0)))) then
         dead_leaves_C = leaf_fall_rate * max(cc%bl,0.0)
         dead_leaves_N = leaf_fall_rate * min(cc%leaf_N,0.0)
         dead_roots_C  = root_mort_rate * max(cc%br,0.0)
         dead_roots_N  = root_mort_rate * max(cc%root_N,0.0)
         dead_stem_C   = 0.0

         ! isa 20170705
         if (sp%lifeform == FORM_GRASS) then
            dead_stem_C = min(stem_mort_rate * cc%bsw_max, &
                   cc%bsw - sp%rho_wood * sp%alphaBM * ((sp%gammaHT/(sp%alphaHT/sp%seedling_height - 1.0))**(1.0/sp%thetaHT))**2 * sp%seedling_height)
            dead_stem_C = max(dead_stem_C,0.0)
            ! ToDo - it is necessary to implement anonline adjustment of dbh, height and crown area as the plant shrinks
            !        otherwise it can happen that, in a year with a short winter, there is a disadjustment between plant
            !        biomass and its dimensions
            ! just set initial height
            cc%height = sp%seedling_height
            cc%dbh = (sp%gammaHT/(sp%alphaHT/sp%seedling_height - 1.0))**(1.0/sp%thetaHT)
            cc%crownarea = sp%alphaCA * cc%dbh**sp%thetaCA
            ! slm 20170804: reconcile bl_max and crownarea, and drop excess leaf biomass
            !               immediately
            cc%bl_max = sp%LMA * sp%LAImax * cc%crownarea * (1.0-sp%internal_gap_frac)
            dead_leaves_C = max(dead_leaves_C, cc%bl-cc%bl_max)
            if (cc%bl > 0) then
               dead_leaves_N = min(cc%leaf_N * dead_leaves_C/cc%bl, cc%leaf_N)
            else
               dead_leaves_N = cc%leaf_N
            endif
         endif
         ! drop all leaves if LAI is below certain small minimum
         if (cc%lai<min_lai) then
            dead_leaves_C = cc%bl
            dead_leaves_N = cc%leaf_N
         endif
         if (cc%bsw>0) then
            dead_stem_N = min(cc%sapwood_N*dead_stem_C/cc%bsw,cc%sapwood_N)
         else
            dead_stem_N = cc%sapwood_N
         endif
         ! update C and N pools
         cc%nsc      = cc%nsc + l_fract*(dead_leaves_C+dead_roots_C+dead_stem_C) ! isa 20170705
         cc%stored_N = cc%stored_N + (sp%leaf_retranslocation_frac*dead_leaves_N+sp%froot_retranslocation_frac*dead_roots_N)
         cc%bl       = cc%bl - dead_leaves_C ; cc%leaf_N    = cc%leaf_N - dead_leaves_N
         cc%br       = cc%br - dead_roots_C  ; cc%root_N    = cc%root_N - dead_roots_N
         cc%bsw      = cc%bsw - dead_stem_C  ; cc%sapwood_N = cc%sapwood_N - dead_stem_N ! isa 20170705

         cc%lai = leaf_area_from_biomass(cc%bl,cc%species,cc%layer,cc%firstlayer)/ &
                                        (cc%crownarea *(1.0-sp%internal_gap_frac))
         if(cc%bl==0) cc%leaf_age = 0.0
         cc%bliving = cc%blv + cc%br + cc%bl + cc%bsw

         leaf_litter_C = (1-l_fract) * (dead_leaves_C+dead_roots_C+dead_stem_C) * cc%nindivs ! isa20170705, stem and roots become leaf litter
         leaf_litter_N = (1-sp%leaf_retranslocation_frac) * dead_leaves_N * cc%nindivs &
                       + dead_stem_N * cc%nindivs
         leaf_litt_C(:) = leaf_litt_C(:)+[sp%fsc_liv,1-sp%fsc_liv,0.0]*leaf_litter_C
         leaf_litt_N(:) = leaf_litt_N(:)+[sp%fsc_liv,1-sp%fsc_liv,0.0]*leaf_litter_N

         vegn%litter = vegn%litter + leaf_litter_C
         soil%fsc_in(1)  = soil%fsc_in(1) + leaf_litter_C
         vegn%veg_out = vegn%veg_out + leaf_litter_C

         root_litter_C = (1-l_fract) * dead_roots_C * cc%nindivs
         root_litter_N = (1-sp%froot_retranslocation_frac) * dead_roots_N
         call cohort_root_litter_profile(cc, dz, profile)
         do l = 1, num_l
            root_litt_C(l,:) = root_litt_C(l,:) + profile(l)* &
                 [sp%fsc_froot, 1-sp%fsc_froot, 0.0]*root_litter_C
            root_litt_N(l,:) = root_litt_N(l,:) + profile(l)* &
                 [sp%fsc_froot, 1-sp%fsc_froot, 0.0]*root_litter_N
         enddo
     endif
     end associate ! cc, sp
  enddo
  ! add litter accumulated over the cohorts
  call add_soil_carbon(soil, vegn, leaf_litter_C=leaf_litt_C, leaf_litter_N=leaf_litt_N)
  ! phenology can change cohort heights if the grass dies, and therefore change
  ! layers -- therefore we need to relayer, lest cohorts remain in a wrong order
  call vegn_relayer_cohorts_ppa(vegn)
  ! merge similar cohorts, otherwise their number proliferates due to re-layering
  call vegn_mergecohorts_ppa(vegn, dheat)
  tile%e_res_2 = tile%e_res_2 - dheat
  end associate ! vegn, soil
end subroutine vegn_phenology_ppa

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
  integer :: i,k
  real :: delta
  real :: deltafast, deltaslow, deltafast_N, deltaslow_N
  real :: profile(num_l), profile1(num_l), psum ! for depostion profile calculation
  real :: litterC(num_l,N_C_TYPES) ! soil litter C input by layer and type
  real :: litterN(num_l,N_C_TYPES) ! soil litter N input by layer and type
  real, dimension(N_C_TYPES,N_LITTER_POOLS) :: delta_C, delta_N

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
  case (SOILC_CORPSE,SOILC_CORPSE_N)
     ! update fsc input rate so that intermediate fsc pool is never
     ! depleted below zero; on the other hand the pool can be only
     ! depleted, never increased
     vegn%fsc_rate_ag = MAX( 0.0, MIN(vegn%fsc_rate_ag, vegn%fsc_pool_ag/dt_fast_yr));
     deltafast = vegn%fsc_rate_ag*dt_fast_yr;
     vegn%fsc_pool_ag       = vegn%fsc_pool_ag       - deltafast;

     ! update ssc input rate so that intermediate ssc pool is never
     ! depleted below zero; on the other hand the pool can be only
     ! depleted, never increased
     vegn%ssc_rate_ag = MAX(0.0, MIN(vegn%ssc_rate_ag, vegn%ssc_pool_ag/dt_fast_yr));
     deltaslow = vegn%ssc_rate_ag*dt_fast_yr;
     vegn%ssc_pool_ag       = vegn%ssc_pool_ag       - deltaslow;

     vegn%litter_rate_C = MAX(0.0, MIN(vegn%litter_rate_C, vegn%litter_buff_C/dt_fast_yr))
     delta_C = vegn%litter_rate_C*dt_fast_yr

     if(soil_carbon_option == SOILC_CORPSE_N) then
        vegn%litter_rate_N = MAX(0.0, MIN(vegn%litter_rate_N, vegn%litter_buff_N/dt_fast_yr))
     else
        vegn%litter_rate_N = 0.0
     endif
     delta_N = vegn%litter_rate_N*dt_fast_yr

     do i = 1,N_LITTER_POOLS
        call add_litter(soil%litter(i), delta_C(:,i), delta_N(:,i))
     enddo
     vegn%litter_buff_C = vegn%litter_buff_C - delta_C
     vegn%litter_buff_N = vegn%litter_buff_N - delta_N

     ! update fsc input rate so that intermediate fsc pool is never
     ! depleted below zero; on the other hand the pool can be only
     ! depleted, never increased
     vegn%fsc_rate_bg = MAX( 0.0, MIN(vegn%fsc_rate_bg, vegn%fsc_pool_bg/dt_fast_yr));
     deltafast        = vegn%fsc_rate_bg*dt_fast_yr;
     vegn%fsc_pool_bg = vegn%fsc_pool_bg - deltafast;

     ! update ssc input rate so that intermediate ssc pool is never
     ! depleted below zero; on the other hand the pool can be only
     ! depleted, never increased
     vegn%ssc_rate_bg = MAX(0.0, MIN(vegn%ssc_rate_bg, vegn%ssc_pool_bg/dt_fast_yr));
     deltaslow        = vegn%ssc_rate_bg*dt_fast_yr;
     vegn%ssc_pool_bg = vegn%ssc_pool_bg - deltaslow;

     if (soil_carbon_option == SOILC_CORPSE_N) then
        vegn%fsn_rate_bg = MAX( 0.0, MIN(vegn%fsn_rate_bg, vegn%fsn_pool_bg/dt_fast_yr));
        deltafast_N      = vegn%fsn_rate_bg*dt_fast_yr;
        vegn%fsn_pool_bg = vegn%fsn_pool_bg - deltafast_N;

        vegn%ssn_rate_bg = MAX(0.0, MIN(vegn%ssn_rate_bg, vegn%ssn_pool_bg/dt_fast_yr));
        deltaslow_N      = vegn%ssn_rate_bg*dt_fast_yr;
        vegn%ssn_pool_bg = vegn%ssn_pool_bg - deltaslow_N;
     else
        vegn%fsn_rate_bg = 0.0
        deltafast_N      = 0.0
        vegn%fsn_pool_bg = 0.0

        vegn%ssn_rate_bg = 0.0
        deltaslow_N      = 0.0
        vegn%ssn_pool_bg = 0.0
     endif

     ! vertical profile of litter is proportional to the average of liiter profiles
     ! of all cohorts, weighted with biomasses of fine roots. This doesn't seem to
     ! be a very good assumption, since fine roots sometimes die (mass is zero),
     ! but profile should not be zero in this case.
     profile(:) = 0.0
     do i = 1,vegn%n_cohorts
        associate(cc=>vegn%cohorts(i))
        call cohort_root_litter_profile(cc,dz,profile1)
        profile(:) = profile(:) + profile1(:)*cc%br*cc%nindivs
        end associate
     enddo
     psum = sum(profile)
     if (psum>0) then
        profile(:) = profile(:)/psum
     else
        profile(:) = 0.0
        profile(1) = 1.0
     endif
     do k = 1,num_l
        litterC(k,:) = [deltafast,deltaslow,0.0] * profile(k)
        litterN(k,:) = [deltafast_N,deltaslow_N,0.0] * profile(k)
     enddo
     call add_root_litter(soil, vegn, litterC, litterN )
  case default
     call error_mesg('update_soil_pools','The value of soil_carbon_option is invalid. This should never happen. Contact developer.',FATAL)
  end select
end subroutine update_soil_pools


! ============================================================================
function cohort_can_reproduce(cc); logical cohort_can_reproduce
  type(vegn_cohort_type), intent(in) :: cc

  cohort_can_reproduce = (cc%layer == 1 .and. cc%age > spdata(cc%species)%maturalage)
end function cohort_can_reproduce

! ============================================================================
! returns TRUE if cohort can allocate carbon to seeds.
! currently it is the same as cohort_can_reproduce, but it will be different
! in ppa-N
function cohort_makes_seeds(cc, G_WF) result(answer); logical answer
  type(vegn_cohort_type), intent(in) :: cc
  real, intent(in) :: G_WF ! amount of carbon spent on new sapwood growth and seeds

  answer = (cc%layer == 1 .and. cc%age > spdata(cc%species)%maturalage)
  if (soil_carbon_option==SOILC_CORPSE_N.AND.N_limits_live_biomass) then
     associate(sp=>spdata(cc%species))
     answer = answer .AND. &
        .NOT.(cc%nitrogen_stress > sp%max_n_stress_for_seed_production &
              .OR. sp%v_seed*G_WF/sp%seed_c2n>0.1*cc%stored_N )
     end associate
  endif
end function cohort_makes_seeds

! ============================================================================
subroutine vegn_reproduction_ppa(do_seed_transport)
  logical, intent(in) :: do_seed_transport

  type(land_tile_enum_type) :: ce
  type(land_tile_type), pointer :: tile
  real, dimension(lnd%ls:lnd%le,0:nspecies-1) :: &
     ug_dispersed_C,   ug_dispersed_N,   &  ! dispersed seeds, kg per m2 of land
     ug_transported_C, ug_transported_N     ! dispersed seeds, kg per m2 of land
  real, allocatable :: seed_C(:,:), seed_N(:,:)
  integer :: n ! total number of vegetated tiles in the domain
  integer :: k ! vegetated tile index
  integer :: l ! grid cell index in unstructured grid
  integer :: s ! species iterator
  integer :: i ! cohort iterator
  integer :: k1

  real :: btot0, btot1 ! total carbon, for conservation check only
  real :: ntot0, ntot1 ! total nitrogen, for conservation check only

  ! + conservation check part 1
  if (do_check_conservation) then
     btot0 = 0.0; ntot0 = 0.0
     ce = first_elmt(land_tile_map,lnd%ls)
     do while (loop_over_tiles(ce,tile,l))
        btot0 = btot0 + lnd%ug_area(l) * tile%frac * land_tile_carbon(tile)
        ntot0 = ntot0 + lnd%ug_area(l) * tile%frac * land_tile_nitrogen(tile)
     end do
     call mpp_sum(btot0); call mpp_sum(ntot0)
  end if
  ! - conservation check part 1

  ! count total number of vegetation tiles in our domain
  n  = 0
  ce = first_elmt(land_tile_map,lnd%ls)
  do while (loop_over_tiles(ce,tile,l))
     if(associated(tile%vegn)) n = n+1
  end do

  ! calculate amount of seeds (kgC per tile area, by species) for each of the vegetation tiles
  allocate(seed_C(n,0:nspecies-1), seed_N(n,0:nspecies-1))
  seed_C = 0.0; seed_N = 0.0
  ce = first_elmt(land_tile_map, lnd%ls); k = 1
  do while (loop_over_tiles(ce,tile,l))
     if(.not.associated(tile%vegn)) cycle
     do i = 1,tile%vegn%n_cohorts
        associate(cc=>tile%vegn%cohorts(i))
        if (cohort_can_reproduce(cc)) then
           seed_C(k,cc%species) = seed_C(k,cc%species) + cc%nindivs * cc%bseed
           seed_N(k,cc%species) = seed_N(k,cc%species) + cc%nindivs * cc%seed_N
           cc%bseed = 0.0; cc%seed_N = 0.0
        endif
        end associate
     enddo
     k = k+1
  enddo

  ! calculate amount of dispersed seeds, by species, kgC per m2 of land
  ug_dispersed_C = 0.0   ; ug_dispersed_N = 0.0
  ug_transported_C = 0.0 ; ug_transported_N = 0.0
  if (do_seed_transport) then
     ce = first_elmt(land_tile_map, lnd%ls); k = 1
     do while (loop_over_tiles(ce,tile,l)) ! l is the grid cell index in unstructured grid
        if (.not.associated(tile%vegn)) cycle
        do s = 0,nspecies-1
           ug_dispersed_C(l,s)   = ug_dispersed_C(l,s)   + seed_C(k,s)*spdata(s)%frac_seed_dispersed*tile%frac
           ug_dispersed_N(l,s)   = ug_dispersed_N(l,s)   + seed_N(k,s)*spdata(s)%frac_seed_dispersed*tile%frac
           ug_transported_C(l,s) = ug_transported_C(l,s) + seed_C(k,s)*spdata(s)%frac_seed_transported*tile%frac
           ug_transported_N(l,s) = ug_transported_N(l,s) + seed_N(k,s)*spdata(s)%frac_seed_transported*tile%frac
           ! correct the amount of seeds remaining in each tile
           seed_C(k,s) = seed_C(k,s)*(1-spdata(s)%frac_seed_dispersed-spdata(s)%frac_seed_transported)
           seed_N(k,s) = seed_N(k,s)*(1-spdata(s)%frac_seed_dispersed-spdata(s)%frac_seed_transported)
        enddo
        k = k+1
     enddo
     ! diffuse the seeds
     do s = 0,nspecies-1
        call transport_seeds(ug_transported_C(:,s))
        call transport_seeds(ug_transported_N(:,s))
     enddo
  endif

  ! distribute seeds (those that stay in place + dispersed) among tiles
  ce = first_elmt(land_tile_map, lnd%ls); k = 1
  do while (loop_over_tiles(ce,tile,l,k1))
     call set_current_point(l,k1)
     if (.not.associated(tile%vegn)) cycle
     if (is_watch_point()) then
        __DEBUG3__(l,k,k1)
        do s = 0,nspecies-1
           write(*,'(i2.2,x,a32)',advance='NO')s,spdata(s)%name
           call dpri('dispersed_C',ug_dispersed_C(l,s))
           call dpri('transported_C',ug_transported_C(l,s))
           call dpri('seed_C',seed_C(k,s))
           call dpri('dispersed_N',ug_dispersed_N(l,s))
           call dpri('transported_N',ug_transported_N(l,s))
           call dpri('seed_N',seed_N(k,s))
           write(*,*)
        enddo
     endif
     if (tile%vegn%landuse==LU_CROP .and. .not.allow_weeds_on_crops) then
        call add_seedlings_ppa(tile%vegn,tile%soil,(ug_dispersed_C(l,:)+ug_transported_C(l,:))*ug_area_factor(l)+seed_C(k,:), &
                                                   (ug_dispersed_N(l,:)+ug_transported_N(l,:))*ug_area_factor(l)+seed_N(k,:), &
                               germination_factor = 0.0) ! no seeds germinate
        ! This also means that the crops are not allowed to reproduce by themselves.
     else
        call add_seedlings_ppa(tile%vegn,tile%soil,(ug_dispersed_C(l,:)+ug_transported_C(l,:))*ug_area_factor(l)+seed_C(k,:), &
                                                   (ug_dispersed_N(l,:)+ug_transported_N(l,:))*ug_area_factor(l)+seed_N(k,:), &
                               germination_factor = 1.0) ! do not modify natural seed germination
     endif
     k = k+1
  enddo

  ! + conservation check part 2
  if (do_check_conservation) then
     btot1 = 0.0; ntot1 = 0.0
     ce = first_elmt(land_tile_map,lnd%ls)
     do while (loop_over_tiles(ce,tile,l))
        btot1 = btot1 + lnd%ug_area(l) * tile%frac * land_tile_carbon(tile)
        ntot1 = ntot1 + lnd%ug_area(l) * tile%frac * land_tile_nitrogen(tile)
     end do
     call mpp_sum(btot1) ; call mpp_sum(ntot1)
     if (mpp_pe()==mpp_root_pe()) then
        call check_conservation ('vegn_reproduction_ppa','total carbon', &
             btot0/atot, btot1/atot, carbon_cons_tol, severity=FATAL)
        call check_conservation ('vegn_reproduction_ppa','total nitrogen', &
             ntot0/atot, ntot1/atot, nitrogen_cons_tol, severity=FATAL)
     endif
  end if
  ! - conservation check part 2

  deallocate(seed_C, seed_N)
end subroutine vegn_reproduction_ppa

! =======================================================================================
! Given the amount seeds undergoing transport (kg per m2 of land), on unstructured grid,
! updates it to take into account transport among grid cells.
subroutine transport_seeds(ug_bseed)
  real, intent(inout) :: ug_bseed(lnd%ls:lnd%le) ! amount of transported seeds on unstructured
       ! grid, kg per m2 of land

  real :: sg_bseed(lnd%isd:lnd%ied, lnd%jsd:lnd%jed) ! total amount of dispersed seeds per grid cell on structured grid, kg
  real :: tend    (lnd%isd:lnd%ied, lnd%jsd:lnd%jed) ! seed mass tendency due to dispersion, kg
  integer :: i,j,ii,jj

  real :: btot0, btot1 ! total carbon, for conservation check only

  real, parameter :: kernel(-1:1,-1:1) = reshape([ & ! shape of the dispersal function
      0.0,  0.25, 0.0,  &
      0.25, 0.0,  0.25, &
      0.0,  0.25, 0.0   ],[3,3] )

  if(do_check_conservation) then
     btot0 = sum(ug_bseed*lnd%ug_area)
     call mpp_sum(btot0)
  endif
  ! move dispersed seeds to structured grid
  sg_bseed = 0.0
  call mpp_pass_UG_to_SG(lnd%ug_domain,ug_bseed*lnd%ug_area,sg_bseed)
  ! NOTE that bseed is multiplied by area, so it is total per grid cell

  ! update halo
  call mpp_update_domains(sg_bseed,lnd%sg_domain)

  tend = 0.0
  do j = lnd%js, lnd%je
  do i = lnd%is, lnd%ie
     do jj = -1,1
     do ii = -1,1
         tend(i,j) = tend(i,j) + sg_bseed(i-ii,j-jj)*kernel(ii,jj)*sg_soilfrac(i,j)
     enddo
     enddo
  enddo
  enddo
  do j = lnd%js, lnd%je
  do i = lnd%is, lnd%ie
     do jj = -1,1
     do ii = -1,1
         tend(i,j) = tend(i,j) - sg_bseed(i,j)*kernel(ii,jj)*sg_soilfrac(i+ii,j+jj)
     enddo
     enddo
  enddo
  enddo
  ! tendency and updated seed amount are only valid within compute domain, but that's OK
  ! because we only use them there.
  sg_bseed = sg_bseed + tend

  ! move transported field to unstructured grid
  call mpp_pass_SG_to_UG(lnd%ug_domain,sg_bseed,ug_bseed)
  ! renormalize seed amount from total to kg C per unit land area
  ug_bseed(:) = ug_bseed(:)/lnd%ug_area(:)

  if (do_check_conservation) then
     btot1 = sum(ug_bseed*lnd%ug_area)
     call mpp_sum(btot1)
     if (mpp_pe()==mpp_root_pe()) then
        call check_conservation ('transport_seeds','total carbon', &
             btot0/atot, btot1/atot, carbon_cons_tol, severity=FATAL)
     endif
  end if
end subroutine transport_seeds

end module vegn_dynamics_mod
