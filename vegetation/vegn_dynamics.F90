! ============================================================================
! updates carbon pools and rates on the fast time scale
! ============================================================================
module vegn_dynamics_mod

#include "../shared/debug.inc"

use fms_mod, only: error_mesg, FATAL, WARNING
use time_manager_mod, only: time_type

use constants_mod, only : PI,tfreeze
use land_constants_mod, only : days_per_year, seconds_per_year, mol_C
use land_data_mod, only : log_version
use land_debug_mod, only : is_watch_point, check_var_range
use land_tile_diag_mod, only : OP_SUM, OP_AVERAGE, &
     register_tiled_diag_field, send_tile_data, diag_buff_type, &
     register_cohort_diag_field, send_cohort_data, set_default_diag_filter, OP_SUM
use vegn_data_mod, only : spdata, &
     PHEN_DECIDUOUS, LEAF_ON, LEAF_OFF, FORM_WOODY, FORM_GRASS, &
     ALLOM_EW, ALLOM_EW1, ALLOM_HML, &
     agf_bs, l_fract, mcv_min, mcv_lai, do_ppa, tau_seed, &
     understory_lai_factor, wood_fract_min, &
     myc_scav_C_efficiency, myc_mine_C_efficiency, N_fixer_C_efficiency, N_limits_live_biomass, &
     excess_stored_N_leakage_rate, &
     c2n_N_fixer, et_myc, &
     mycorrhizal_turnover_time, N_fixer_turnover_time, spec_data_type
use vegn_tile_mod, only: vegn_tile_type, vegn_tile_carbon, vegn_tile_nitrogen
use soil_tile_mod, only: num_l, dz, soil_tile_type, clw, csw, N_LITTER_POOLS, LEAF, CWOOD, soil_tile_carbon
use vegn_cohort_mod, only : vegn_cohort_type, &
     update_biomass_pools, update_bio_living_fraction, update_species, &
     leaf_area_from_biomass, biomass_of_individual, init_cohort_allometry_ppa, &
     cohort_root_litter_profile, cohort_root_exudate_profile
use vegn_disturbance_mod, only : kill_plants_ppa
use soil_carbon_mod, only: N_C_TYPES, soil_carbon_option, &
    SOILC_CENTURY, SOILC_CENTURY_BY_LAYER, SOILC_CORPSE, SOILC_CORPSE_N, &
    add_litter, debug_pool,soil_NO3_deposition,soil_NH4_deposition,soil_org_N_deposition,deadmic_slow_frac
use soil_util_mod, only: add_soil_carbon, add_root_litter, add_root_exudates
use soil_mod, only: Dsdt, root_N_uptake, myc_scavenger_N_uptake, myc_miner_N_uptake

implicit none
private

! ==== public interfaces =====================================================
public :: vegn_dynamics_init

public :: vegn_carbon_int_lm3   ! fast time-scale integrator of carbon balance
public :: vegn_carbon_int_ppa   ! fast time-scale integrator of carbon balance
public :: vegn_growth           ! slow time-scale redistributor of accumulated carbon
public :: vegn_phenology_lm3
public :: vegn_phenology_ppa
public :: vegn_biogeography

public :: vegn_starvation_ppa   !
public :: vegn_reproduction_ppa ! reproduction for PPA case
public :: kill_small_cohorts_ppa
! ==== end of public interfaces ==============================================

! ==== module constants ======================================================
character(len=*), parameter :: module_name = 'vegn_dynamics_mod'
#include "../shared/version_variable.inc"
character(len=*), parameter :: diag_mod_name = 'vegn'

real, parameter :: GROWTH_RESP=0.333  ! fraction of NPP lost as growth respiration

! ==== module data ===========================================================
real    :: dt_fast_yr ! fast (physical) time step, yr (year is defined as 365 days)

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
    id_exudate


contains

! ============================================================================
subroutine vegn_dynamics_init(id_ug, time, delta_time)
  integer        , intent(in) :: id_ug   !<Unstructured axis id.
  type(time_type), intent(in) :: time       ! initial time for diagnostic fields
  real           , intent(in) :: delta_time ! fast time step, s

  call log_version(version, module_name, &
  __FILE__)

  ! set up global variables
  dt_fast_yr = delta_time/seconds_per_year

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
end subroutine vegn_dynamics_init


subroutine  update_mycorrhizae(cc,sp,dt_fast_yr,&
                           C_allocation_to_N_acq,myc_scav_N_uptake,myc_mine_N_uptake,myc_mine_C_uptake,root_N_uptake,&
                           myc_scav_efficiency,myc_mine_efficiency,&
                           scavenger_myc_C_allocated,miner_myc_C_allocated,N_fixer_C_allocated,&
                           myc_scav_marginal_gain,myc_mine_marginal_gain,N_fix_marginal_gain,rhiz_exud_marginal_gain,&
                           root_exudate_C, root_exudate_N,&
                           scav_N_to_plant,mine_N_to_plant,fix_N_to_plant,&
                           myc_CO2_prod,N_fixation,total_plant_N_uptake,&
                           myc_turnover_C,myc_turnover_N,myc_Nmin)


   type(vegn_cohort_type), intent(inout) :: cc
   type(spec_data_type),   intent(in)    :: sp
   real, intent(in) :: dt_fast_yr
   real, intent(in) :: myc_scav_N_uptake,myc_mine_N_uptake,myc_mine_C_uptake,root_N_uptake  ! N uptake per individual (kgN/m2)
   real, intent(in) :: C_allocation_to_N_acq
   real, intent(in) :: myc_scav_efficiency,myc_mine_efficiency         ! N uptake/myc biomass -- Only used if myc biomass is zero
   real, intent(out):: scavenger_myc_C_allocated,miner_myc_C_allocated,N_fixer_C_allocated ! Allocation per individual (kgC/m2)
   real, intent(out):: myc_scav_marginal_gain,myc_mine_marginal_gain,N_fix_marginal_gain,rhiz_exud_marginal_gain ! kgN/kgC allocated
   real, intent(out):: root_exudate_C, root_exudate_N
   real, intent(out):: scav_N_to_plant,mine_N_to_plant,fix_N_to_plant
   real, intent(out):: myc_CO2_prod, myc_Nmin  ! kgC/m2/individual
   real, intent(out):: N_fixation, total_plant_N_uptake
   real, intent(out):: myc_turnover_C,myc_turnover_N

   real :: myc_scav_exudate_frac,myc_mine_exudate_frac,N_fixer_exudate_frac,rhiz_exud_frac
   real :: scavenger_myc_N_allocated,miner_myc_N_allocated,N_fixer_N_allocated
   real :: fixer_biomass_N_fixation
   real :: d_scav_C_reservoir,d_scav_N_reservoir,d_mine_C_reservoir,d_mine_N_reservoir,d_N_fixer_C_reservoir,d_N_fixer_N_reservoir
   real :: scavenger_myc_growth, miner_myc_growth, N_fixer_growth
   real :: lim_factor
   real :: wood_n2c
   real :: reservoir_C_leakage, maint_resp

   myc_CO2_prod = 0.0
   myc_Nmin = 0.0
   reservoir_C_leakage = 0.0


   if (soil_carbon_option == SOILC_CORPSE_N) then

     ! Update reservoirs with N uptake for this cohort
     cc%scav_myc_N_reservoir = cc%scav_myc_N_reservoir+myc_scav_N_uptake
     cc%mine_myc_C_reservoir = cc%mine_myc_C_reservoir+myc_mine_C_uptake
     cc%mine_myc_N_reservoir = cc%mine_myc_N_reservoir+myc_mine_N_uptake
     N_fixation = cc%N_fixer_biomass_C*sp%N_fixation_rate*dt_fast_yr
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



         call check_var_range(cc%myc_scavenger_biomass_C, 0.0, HUGE(1.0), 'vegn_carbon_int_lm3', 'myc_scavenger_C', FATAL)
         call check_var_range(cc%myc_miner_biomass_C,     0.0, HUGE(1.0), 'vegn_carbon_int_lm3', 'myc_miner_C', FATAL)
         call check_var_range(cc%N_fixer_biomass_C,       0.0, HUGE(1.0), 'vegn_carbon_int_lm3', 'N_fixer_C', FATAL)


         if(is_watch_point()) then
            __DEBUG3__(cc%myc_scavenger_biomass_C,cc%myc_miner_biomass_C,cc%N_fixer_biomass_C)
         endif

      if(cc%myc_scavenger_biomass_C<0) call error_mesg('vegn_carbon_int','Mycorrhizal scavenger biomass < 0',FATAL)
      if(cc%myc_miner_biomass_C<0) then
        __DEBUG3__(cc%myc_miner_biomass_C,myc_mine_C_uptake,miner_myc_C_allocated)
        call error_mesg('vegn_carbon_int','Mycorrhizal miner biomass < 0',FATAL)
      endif
      if(cc%N_fixer_biomass_C<0) call error_mesg('vegn_carbon_int','N fixer biomass < 0',FATAL)

      if(cc%myc_scavenger_biomass_N<0) call error_mesg('vegn_carbon_int','Mycorrhizal scavenger biomass < 0',FATAL)
      if(cc%myc_miner_biomass_N<0) then
        __DEBUG3__(cc%myc_miner_biomass_N,myc_mine_N_uptake,miner_myc_N_allocated)
        call error_mesg('vegn_carbon_int','Mycorrhizal miner biomass N < 0',FATAL)
      endif
      if(cc%N_fixer_biomass_N<0) call error_mesg('vegn_carbon_int','N fixer biomass < 0',FATAL)



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
      reservoir_C_leakage = reservoir_C_leakage + (cc%scav_myc_C_reservoir + cc%mine_myc_C_reservoir + cc%N_fixer_C_reservoir)*dt_fast_yr*365
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
      if(cc%scav_myc_N_reservoir<-1e-10) call error_mesg('vegn_carbon_int','Mycorrhizal scavenger N reservoir < 0',FATAL)
      if(cc%mine_myc_N_reservoir<-1e-10) then
          __DEBUG4__(cc%mine_myc_N_reservoir,cc%mine_myc_C_reservoir,cc%myc_miner_biomass_C,cc%myc_miner_biomass_N)
          __DEBUG4__(cc%species,cc%total_N,cc%leaf_N,cc%stored_N)
          call error_mesg('vegn_carbon_int','Mycorrhizal miner N reservoir < 0',FATAL)
      endif
      if(cc%N_fixer_N_reservoir<-1e-10) call error_mesg('vegn_carbon_int','N fixer N reservoir < 0',FATAL)

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
         rhiz_exud_marginal_gain = (root_N_uptake/dt_fast_yr)/(C_allocation_to_N_acq)!+(myc_mine_marginal_gain+myc_scav_marginal_gain)*0.5
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

      ! Calculate relative fractions
      if (myc_scav_marginal_gain+N_fix_marginal_gain+myc_mine_marginal_gain+rhiz_exud_marginal_gain>0) then
         myc_scav_exudate_frac = myc_scav_marginal_gain/(myc_scav_marginal_gain+myc_mine_marginal_gain+N_fix_marginal_gain+rhiz_exud_marginal_gain)
         myc_mine_exudate_frac = myc_mine_marginal_gain/(myc_scav_marginal_gain+myc_mine_marginal_gain+N_fix_marginal_gain+rhiz_exud_marginal_gain)
         N_fixer_exudate_frac  = N_fix_marginal_gain   /(myc_scav_marginal_gain+myc_mine_marginal_gain+N_fix_marginal_gain+rhiz_exud_marginal_gain)
         rhiz_exud_frac = rhiz_exud_marginal_gain/(myc_scav_marginal_gain+myc_mine_marginal_gain+N_fix_marginal_gain+rhiz_exud_marginal_gain)
      else
         ! Divide evenly if there is no marginal gain.  But this probably only happens if C_allocation_to_N_acq is zero
         myc_scav_exudate_frac = 0.4*0.7
         myc_mine_exudate_frac = 0.3*0.7
         N_fixer_exudate_frac = 0.3*0.7
         rhiz_exud_frac = 0.3
      endif

      N_fixer_C_allocated = C_allocation_to_N_acq*N_fixer_exudate_frac*dt_fast_yr
      miner_myc_C_allocated=C_allocation_to_N_acq*myc_mine_exudate_frac*dt_fast_yr
      scavenger_myc_C_allocated=C_allocation_to_N_acq*myc_scav_exudate_frac*dt_fast_yr

      N_fixer_N_allocated = 0.0
      miner_myc_N_allocated = miner_myc_C_allocated*sp%root_exudate_N_frac
      scavenger_myc_N_allocated = scavenger_myc_C_allocated*sp%root_exudate_N_frac

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


      total_plant_N_uptake = scav_N_to_plant + mine_N_to_plant + fix_N_to_plant + root_N_uptake


   cc%stored_N = cc%stored_N + total_plant_N_uptake

   root_exudate_N = C_allocation_to_N_acq*dt_fast_yr*sp%root_exudate_N_frac - scavenger_myc_N_allocated - N_fixer_N_allocated - miner_myc_N_allocated

 else  ! If nitrogen turned off, everything is zero

      scavenger_myc_C_allocated = 0.0
      miner_myc_C_allocated = 0.0
      N_fixer_C_allocated = 0.0
      scavenger_myc_N_allocated = 0.0
      miner_myc_N_allocated = 0.0
      N_fixer_N_allocated = 0.0
      myc_scav_marginal_gain = 0.0
      myc_mine_marginal_gain = 0.0
      N_fix_marginal_gain = 0.0
      rhiz_exud_marginal_gain = 0.0
      root_exudate_N = 0.0
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

   !root_exudate_N=(current_root_exudation*root_exudate_N_frac*dt_fast_yr - scavenger_myc_N_allocated - N_fixer_N_allocated - miner_myc_N_allocated)/dt_fast_yr
   root_exudate_C = C_allocation_to_N_acq*dt_fast_yr - scavenger_myc_C_allocated - miner_myc_C_allocated - N_fixer_C_allocated + reservoir_C_leakage


   if(root_exudate_N<0) then
     __DEBUG2__(root_exudate_C,root_exudate_N)
     __DEBUG3__(scavenger_myc_N_allocated , N_fixer_N_allocated , miner_myc_N_allocated)
     call error_mesg('update_mycorrhizae','Root exudate N < 0',FATAL)
   endif

   cc%stored_N=cc%stored_N - root_exudate_N - scavenger_myc_N_allocated - N_fixer_N_allocated - miner_myc_N_allocated


   if (is_watch_point()) then
      ! __DEBUG1__(current_root_exudation*dt_fast_yr)
      __DEBUG1__(root_exudate_C)
      __DEBUG1__(scavenger_myc_C_allocated)
      __DEBUG1__(miner_myc_C_allocated)
      __DEBUG1__(myc_mine_C_uptake)
      __DEBUG1__(N_fixer_C_allocated)
      __DEBUG1__(myc_CO2_prod)
   endif


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
  real,dimension(vegn%n_cohorts) :: myc_scav_marginal_gain,myc_mine_marginal_gain, N_fix_marginal_gain, rhiz_exud_marginal_gain
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
        resg(i) = GROWTH_RESP*cc%npp_previous_day
        npp(i)  = npp(i)  - GROWTH_RESP*cc%npp_previous_day
        resp(i) = resp(i) + GROWTH_RESP*cc%npp_previous_day
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
     call update_mycorrhizae(cc,sp,dt_fast_yr,&
                                C_allocation_to_N_acq=C_allocation_to_N_acq,&
                                myc_scav_N_uptake=myc_scav_N_uptake(i),myc_mine_C_uptake=myc_mine_C_uptake(i),myc_mine_N_uptake=myc_mine_N_uptake(i),&
                                root_N_uptake=root_active_N_uptake(i),&
                                myc_scav_efficiency=myc_scav_efficiency,myc_mine_efficiency=myc_mine_efficiency,&
                                scavenger_myc_C_allocated=scavenger_myc_C_allocated(i),miner_myc_C_allocated=miner_myc_C_allocated(i),&
                                N_fixer_C_allocated=N_fixer_C_allocated(i),&
                                myc_scav_marginal_gain=myc_scav_marginal_gain(i),myc_mine_marginal_gain=myc_mine_marginal_gain(i),&
                                N_fix_marginal_gain=N_fix_marginal_gain(i),rhiz_exud_marginal_gain=rhiz_exud_marginal_gain(i),&
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
  call send_cohort_data(id_myc_scav_marginal_gain,diag,c(1:N),myc_scav_marginal_gain(1:N),weight=c(1:N)%nindivs, op=OP_SUM)
  call send_cohort_data(id_myc_mine_marginal_gain,diag,c(1:N),myc_mine_marginal_gain(1:N),weight=c(1:N)%nindivs, op=OP_SUM)
  call send_cohort_data(id_N_fix_marginal_gain,diag,c(1:N),N_fix_marginal_gain(1:N),weight=c(1:N)%nindivs, op=OP_SUM)
  call send_cohort_data(id_rhiz_exudation,diag,c(1:N),root_exudate_C(1:N)/dt_fast_yr,weight=c(1:N)%nindivs, op=OP_SUM)
  call send_cohort_data(id_nitrogen_stress,diag,c(1:N),c(1:N)%nitrogen_stress,weight=c(1:N)%nindivs, op=OP_SUM)
  call send_cohort_data(id_total_plant_N_uptake,diag,c(1:N),total_plant_N_uptake(1:N)/dt_fast_yr,weight=c(1:N)%nindivs, op=OP_SUM)
  call send_cohort_data(id_rhiz_exud_marginal_gain,diag,c(1:N),rhiz_exud_marginal_gain(1:N),weight=c(1:N)%nindivs, op=OP_SUM)

  call send_cohort_data(id_myc_scavenger_N_uptake,diag,c(1:N),myc_scav_N_uptake(1:N)/dt_fast_yr,weight=c(1:N)%nindivs, op=OP_SUM)
  call send_cohort_data(id_myc_miner_N_uptake,diag,c(1:N),myc_mine_N_uptake(1:N)/dt_fast_yr,weight=c(1:N)%nindivs, op=OP_SUM)
  call send_cohort_data(id_symbiotic_N_fixation,diag,c(1:N),N_fixation(1:N)/dt_fast_yr,weight=c(1:N)%nindivs, op=OP_SUM)
  call send_cohort_data(id_active_root_N_uptake,diag,c(1:N), root_active_N_uptake(1:N)/dt_fast_yr,weight=c(1:N)%nindivs, op=OP_SUM)

  call send_cohort_data(id_scav_plant_N_uptake,diag,c(1:N),scav_N_to_plant(1:N)/dt_fast_yr,weight=c(1:N)%nindivs, op=OP_SUM)
  call send_cohort_data(id_mine_plant_N_uptake,diag,c(1:N),mine_N_to_plant(1:N)/dt_fast_yr,weight=c(1:N)%nindivs, op=OP_SUM)
  call send_cohort_data(id_fix_plant_N_uptake,diag,c(1:N),fix_N_to_plant(1:N)/dt_fast_yr,weight=c(1:N)%nindivs, op=OP_SUM)

end subroutine vegn_carbon_int_lm3


! ============================================================================
subroutine vegn_carbon_int_ppa (vegn, soil, tsoil, theta, diag)
  ! TODO: possibly get rid of tsoil, theta, since they can be calculated here
  type(vegn_tile_type), intent(inout) :: vegn
  type(soil_tile_type), intent(inout) :: soil
  real, intent(in) :: tsoil ! average temperature of soil for soil carbon decomposition, deg K
  real, intent(in) :: theta ! average soil wetness, unitless
  type(diag_buff_type), intent(inout) :: diag

  real :: md_wood, md_branch_sw ! in ppa we are losing and replacing branchwood
  real :: deltaBL, deltaBR ! leaf and fine root carbon tendencies
  integer :: i, l, M
  real :: NSC_supply,LR_demand,LR_deficit
  real :: NSCtarget
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
     __DEBUG1__(c%bl)
     __DEBUG1__(c%blv)
     __DEBUG1__(c%br)
     __DEBUG1__(c%bsw)
     __DEBUG1__(c%bwood)
     __DEBUG1__(c%nsc)
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

     cc%carbon_gain  = cc%carbon_gain + gpp(i)*dt_fast_yr - resp(i)*dt_fast_yr

     resp(i) = resp(i) + resg(i) ! add growth respiration and maintenance respiration
     npp(i)  = gpp(i) - resp(i)

     !cc%carbon_gain = cc%carbon_gain-root_exudate_C(i)*dt_fast_yr
     !print *, ' npp ',  npp(i),' root_ex ',  root_exudate_C(i), 'cgain ', cc%carbon_gain

     ! Weng, 2013-01-28
     ! Turnover regardless of STATUS
     deltaBL = cc%bl * sp%alpha_leaf * dt_fast_yr
     deltaBR = cc%br * sp%alpha_root * dt_fast_yr

     cc%bl = cc%bl - deltaBL
     cc%br = cc%br - deltaBR

     ! 20170617: retranslocate part of N back to storage and reduce nitrogen pools; put
     !           retranslocated N into storage
     cc%leaf_N   = cc%leaf_N - deltaBL /sp%leaf_live_c2n*(1.0-sp%leaf_retranslocation_frac)
     cc%root_N   = cc%root_N - deltaBR /sp%froot_live_c2n*(1.0-sp%froot_retranslocation_frac)
     cc%stored_N = cc%stored_N + deltaBL/sp%leaf_live_c2n*sp%leaf_retranslocation_frac   &
                               + deltaBR/sp%froot_live_c2n*sp%froot_retranslocation_frac
     ! compute branch and coarse wood losses for tree types
     md_wood = 0.0
     md_branch_sw = 0.0
     if (spdata(cc%species)%lifeform == FORM_WOODY) then
        md_wood      = sp%branch_wood_frac * Max(cc%bwood,0.0) * sp%alpha_wood * dt_fast_yr
        md_branch_sw = sp%branch_wood_frac * Max(cc%bsw,0.0)   * sp%alpha_wood * dt_fast_yr
     endif
     cc%branch_sw_loss   = cc%branch_sw_loss + md_branch_sw !remember how much was lost over the day
     cc%branch_wood_loss = cc%branch_wood_loss + md_wood

     cc%bsw = cc%bsw - md_branch_sw
     cc%bwood = cc%bwood - md_wood
     ! 20170617: update N pools. For now, assuming that there is no N re-translocation during
     ! either branch loss or turnover of wood and sapwood.
     cc%wood_N    = cc%wood_N    - md_wood/sp%wood_c2n
     cc%sapwood_N = cc%sapwood_N - md_branch_sw/sp%sapwood_c2n

     !reduce nsc by the amount of root exudates during the date
     ! FIXME: when merging with N code take exudates from NSC not carbon gained
     ! FIXME: need different than a fraction of NPP formulation
     !cc%nsc = cc%nsc-root_exudate_C(i)*dt_fast_yr

     ! accumulate liter and soil carbon inputs across all cohorts
     ! 20170617: deposit lost N into litter and soil
     leaf_litt_C(:) = leaf_litt_C(:) + [sp%fsc_liv,  1-sp%fsc_liv,  0.0]*deltaBL*cc%nindivs
     leaf_litt_N(:) = leaf_litt_N(:) + [sp%fsc_liv,  1-sp%fsc_liv,  0.0]*deltaBL/sp%leaf_live_c2n*cc%nindivs*(1.0-sp%leaf_retranslocation_frac)
     wood_litt_C(:) = wood_litt_C(:) + [sp%fsc_wood, 1-sp%fsc_wood, 0.0]*(md_wood+md_branch_sw)*cc%nindivs
!     wood_litt_N(:) = wood_litt_N(:) + cc%nindivs * agf_bs * & -- FIXME: do we need agf_bs in both C and N wood litter?
     wood_litt_N(:) = wood_litt_N(:) + cc%nindivs * &
             [sp%fsc_wood, 1-sp%fsc_wood, 0.0]*(md_wood/sp%wood_c2n+md_branch_sw/sp%sapwood_c2n)
     call cohort_root_litter_profile(cc, dz, profile)
     do l = 1, num_l
        root_litt_C(l,:) = root_litt_C(l,:) + profile(l)*cc%nindivs* &
             [sp%fsc_froot, 1-sp%fsc_froot, 0.0 ]*deltaBR
! 20170617:
        root_litt_N(l,:) = root_litt_N(l,:) + profile(l)*cc%nindivs*( &
             [sp%fsc_froot, 1-sp%fsc_froot, 0.0 ]*deltaBR/sp%froot_live_c2n*(1.0-sp%froot_retranslocation_frac) &
! we are assuming md_wood comes from the branch turnover, so there is nothing
! below ground
!            + [sp%fsc_wood,  1-sp%fsc_wood,  0.0 ]*md_wood*(1-agf_bs)/sp%wood_c2n &
             )
     enddo

     ! resg(i) = 0 ! slm: that doesn't make much sense to me... why?
! 20170617:
     ! Mycorrhizal N uptake

!     C_allocation_to_N_acq = max(npp(i),0.0)*sp%root_exudate_frac
     C_allocation_to_N_acq = cc%nsc/sp%tau_nsc_exudate !This is a rate per year, not per time step
     if (sp%dynamic_root_exudation .AND. soil_carbon_option==SOILC_CORPSE_N) then
        ! 20170617: modify frac for exudate depending on N state
         C_allocation_to_N_acq = C_allocation_to_N_acq*cc%nitrogen_stress
         ! we possibly should impose lower and upper limits on exudate stress factor. Currently,
         ! there is a lower limit on nitrogen_stress, so exudates are never 0. Also nitrogen_stress
         ! cannot rise above 2 in current formulation, but we should be careful if the definition
         ! of stress changes.
         ! N exudates are calculated in update_mycorrhizae.
     endif

     ! this updates N storage in the plant
     ! 20170617: This includes allocation to all N acquisition plus exudates
     call update_mycorrhizae(cc, sp, dt_fast_yr, &
              C_allocation_to_N_acq=C_allocation_to_N_acq, &
              myc_scav_N_uptake=myc_scav_N_uptake(i), myc_mine_C_uptake=myc_mine_C_uptake(i), myc_mine_N_uptake=myc_mine_N_uptake(i), &
              root_N_uptake=root_active_N_uptake(i), &
              myc_scav_efficiency=myc_scav_efficiency, myc_mine_efficiency=myc_mine_efficiency, &
              scavenger_myc_C_allocated=scavenger_myc_C_allocated(i), miner_myc_C_allocated=miner_myc_C_allocated(i), &
              N_fixer_C_allocated=N_fixer_C_allocated(i), &
              myc_scav_marginal_gain=myc_scav_marginal_gain(i), myc_mine_marginal_gain=myc_mine_marginal_gain(i), &
              N_fix_marginal_gain=N_fix_marginal_gain(i), rhiz_exud_marginal_gain=rhiz_exud_marginal_gain(i), &
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
     if(cc%nitrogen_stress<=0.05 .AND. cc%stored_N>0 .AND. soil_carbon_option == SOILC_CORPSE_N) then
       total_N_leakage(:) = total_N_leakage(:) + cc%stored_N*excess_stored_N_leakage_rate*profile(:)*cc%nindivs
       cc%stored_N=cc%stored_N - cc%stored_N*excess_stored_N_leakage_rate*dt_fast_yr
     endif
     cc%nsc = cc%nsc - root_exudate_C(i)

     total_myc_CO2_prod = total_myc_CO2_prod + myc_CO2_prod(i)*cc%nindivs
     total_myc_Nmin(:) = total_myc_Nmin(:) + profile(:)*myc_Nmin*cc%nindivs

     ! increment cohort age
     cc%age = cc%age + dt_fast_yr

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
  vegn%rh = vegn%rh + total_myc_CO2_prod


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
  call send_cohort_data(id_gpp, diag, c(1:M), gpp(1:M), weight=c(1:M)%nindivs, op=OP_SUM)
  call send_cohort_data(id_npp, diag, c(1:M), npp(1:M), weight=c(1:M)%nindivs, op=OP_SUM)
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

!   call send_cohort_data(id_mycorrhizal_scav_allocation,diag,c(1:M),scavenger_myc_C_allocated(1:M)/dt_fast_yr,weight=c(1:M)%nindivs, op=OP_SUM)
!   call send_cohort_data(id_mycorrhizal_mine_allocation,diag,c(1:M),miner_myc_C_allocated(1:M)/dt_fast_yr,weight=c(1:M)%nindivs, op=OP_SUM)
!   call send_cohort_data(id_N_fixer_allocation,diag,c(1:M),N_fixer_C_allocated(1:M)/dt_fast_yr,weight=c(1:M)%nindivs, op=OP_SUM)
!   call send_cohort_data(id_myc_scav_marginal_gain,diag,c(1:M),myc_scav_marginal_gain(1:M),weight=c(1:M)%nindivs, op=OP_SUM)
!   call send_cohort_data(id_myc_mine_marginal_gain,diag,c(1:M),myc_mine_marginal_gain(1:M),weight=c(1:M)%nindivs, op=OP_SUM)
!   call send_cohort_data(id_N_fix_marginal_gain,diag,c(1:M),N_fix_marginal_gain(1:M),weight=c(1:M)%nindivs, op=OP_SUM)
!   call send_cohort_data(id_rhiz_exudation,diag,c(1:M),root_exudate_C(1:M)/dt_fast_yr,weight=c(1:M)%nindivs, op=OP_SUM)
!   call send_cohort_data(id_nitrogen_stress,diag,c(1:M),c(1:M)%nitrogen_stress,weight=c(1:M)%nindivs, op=OP_SUM)
!   call send_cohort_data(id_total_plant_N_uptake,diag,c(1:M),total_plant_N_uptake(1:M)/dt_fast_yr,weight=c(1:M)%nindivs, op=OP_SUM)
!   call send_cohort_data(id_rhiz_exud_marginal_gain,diag,c(1:M),rhiz_exud_marginal_gain(1:M),weight=c(1:M)%nindivs, op=OP_SUM)
!
!   call send_cohort_data(id_myc_scavenger_N_uptake,diag,c(1:M),myc_scav_N_uptake(1:M)/dt_fast_yr,weight=c(1:M)%nindivs, op=OP_SUM)
!   call send_cohort_data(id_myc_miner_N_uptake,diag,c(1:M),myc_mine_N_uptake(1:M)/dt_fast_yr,weight=c(1:M)%nindivs, op=OP_SUM)
!   call send_cohort_data(id_symbiotic_N_fixation,diag,c(1:M),N_fixation(1:M)/dt_fast_yr,weight=c(1:M)%nindivs, op=OP_SUM)
!   call send_cohort_data(id_active_root_N_uptake,diag,c(1:M), root_active_N_uptake(1:M)/dt_fast_yr,weight=c(1:M)%nindivs, op=OP_SUM)
!
!   call send_cohort_data(id_scav_plant_N_uptake,diag,c(1:M),scav_N_to_plant(1:M)/dt_fast_yr,weight=c(1:M)%nindivs, op=OP_SUM)
!   call send_cohort_data(id_mine_plant_N_uptake,diag,c(1:M),mine_N_to_plant(1:M)/dt_fast_yr,weight=c(1:M)%nindivs, op=OP_SUM)
!   call send_cohort_data(id_fix_plant_N_uptake,diag,c(1:M),fix_N_to_plant(1:M)/dt_fast_yr,weight=c(1:M)%nindivs, op=OP_SUM)
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

  if (is_watch_point()) then
     cmass0 = vegn_tile_carbon(vegn)
  endif

 ! write(*,*) 'counting cohorts: ', vegn%n_cohorts
  do i = 1, vegn%n_cohorts
     cc => vegn%cohorts(i)

     if (do_ppa) then

     ! write(*,*)"############## in vegn_growth before #################"
     !__DEBUG1__(cc%carbon_gain)

        cc%nsc=cc%nsc+cc%carbon_gain
        call biomass_allocation_ppa(cc,wood_prod(i),leaf_root_gr(i),sw_seed_gr(i),deltaDBH(i))


     !write(*,*)"############## in vegn_growth #################"

     !__DEBUG1__(cc%carbon_gain)
     !__DEBUG1__(cc%bl)
     !__DEBUG1__(cc%blv)
     !__DEBUG1__(cc%br)
     !__DEBUG1__(cc%bsw)
     !__DEBUG1__(cc%bwood)
     !__DEBUG1__(cc%nsc)

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
       call kill_plants_ppa(cc, vegn, soil, deadtrees, 0.0, &
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
! updates cohort vegetation structure, biomass pools, LAI, SAI, and height
! using accumulated carbon_gain
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
  real :: BSWmax ! max sapwood biomass, kg C/individual
  real :: G_LFR  ! amount of carbon spent on leaf and root growth
  real :: G_WF   ! amount of carbon spent on new sapwood growth and seeds
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
  real :: delta_wood_branch
  real :: br_max_Nstress
  real :: f


  real, parameter :: DBHtp = 0.8 ! m
  logical :: hydraulics_repair = .FALSE.

  !ens
  real :: fsf = 0.2 ! in Weng et al 205, page 2677 units are 0.2 day
  ! new constants for HML equations
!  real ::   alpha_z = 60.97697    ! alphaHT
!  real ::   b_z = 0.8631476       ! thetaHT
!  real ::   gamma_z = 0.6841742   ! gammaHT
!  real ::   alpha_ca = 243.7808   ! alphaCA
!  real ::   b_ca = 1.18162        ! thetaCA
!  real ::   alpha_s = 0.559       ! alphaBM


  associate (sp => spdata(cc%species)) ! F2003

  ! Stress increases as stored N declines relative to total biomass N demand
  ! N stress is calculated based on "potential pools" without N limitation
  ! Elena suggests using 2*(root N + leaf N) as storage target
  ! bliving is already increased after cgain


   NSCtarget = 4.0*cc%bl_max
   N_storage_target = NSCtarget/sp%leaf_live_c2n

  ! 20170617: calculate N stress (function of stored N) here
  !           N traget possibly in relation with NSC target
  if (soil_carbon_option==SOILC_CORPSE_N) then
     cc%nitrogen_stress = (N_storage_target-cc%stored_N)/N_storage_target
     ! if leaves are off, our N stress would be lower, since we re-translocate some of leaf+root N
     ! but plant should be saving enough to rebuild leaf+root next time.
     ! In this case, we should adjust the calculation:
     ! cc%nitrogen_stress = (N_storage_target-cc%bl_max/sp%leaf_live_c2n-cc%br_max/sp%froot_live_c2n-cc%stored_N)/N_storage_target
  else
     cc%nitrogen_stress = 0.0
  endif

  ! TODO: what if carbon_gain is not 0, but leaves are OFF (marginal case? or
  ! typical in lm3?)
  delta_bsw_branch = 0.
  delta_wood_branch  = 0.
  ! set init values for diag output, to be returned when the actual calulations are bypassed:
  wood_prod = 0.0 ; leaf_root_gr = 0.0 ; sw_seed_gr = 0.0 ; deltaDBH = 0.0
  if (cc%status == LEAF_ON) then
     ! update targets based on the nitrogen stress
     ! N_stress_root_factor is something like 0.05 so br_max increases 25% when N stress is 0.5, i.e. stored N is half of target value
     ! N stress or this root factor should not go above some max value, so all biomass doesn't end up in roots in some weird situation
     ! This assumes carbon to build more roots comes from sapwood growth, not leaf growth. But we can make this flexible or have different options
     br_max_Nstress = cc%br_max * (1+cc%nitrogen_stress*sp%N_stress_root_factor)


     ! calculate the carbon spent on growth of leaves and roots
     G_LFR    = max(0.0, min(cc%bl_max+br_max_Nstress-cc%bl-cc%br,  &
                            0.1*cc%nsc/(1+GROWTH_RESP))) ! don't allow more than 0.1/(1+GROWTH_RESP) of nsc per day to spend

     ! and distribute it between roots and leaves
     deltaBL  = min(G_LFR, max(0.0, &
          (G_LFR*cc%bl_max + cc%bl_max*cc%br - br_max_Nstress*cc%bl)/(cc%bl_max + br_max_Nstress) &
          ))
     deltaBR  = G_LFR - deltaBL
     deltaNL = deltaBL*sp%leaf_live_c2n
     deltaNR = deltaBR*sp%froot_live_c2n
     if (deltaNL+deltaNR > 0.1*cc%stored_N.and.N_limits_live_biomass) then
        ! BNS: make this the minimum of 10% NSC depletion OR 10% stored N depletion
        ! Do this every time we do this limitation.
        ! This does not prioritize spending on leaves and roots
        f = 0.1*cc%stored_N/(deltaNL+deltaNR)
        deltaBR = deltaBR*f; deltaBL = deltaBL*f
        deltaNR = deltaNR*f; deltaNL = deltaNL*f
     endif
     cc%stored_N = cc%stored_N - deltaNL - deltaNR
     ! calculate carbon spent on seeds and sapwood growth
     ! Should sapwood be considered as important as leaves and roots?

     !ens first replace lost branch sapwood
     delta_bsw_branch  = 0.0
     if (cc%branch_sw_loss > 0) then
        delta_bsw_branch = max (min(cc%branch_sw_loss, 0.1*cc%nsc/(1+GROWTH_RESP)), 0.0)
        delta_nsw_branch = delta_bsw_branch/sp%sapwood_c2n
        if (delta_nsw_branch>0.1*cc%stored_N.and.N_limits_live_biomass) then
           delta_bsw_branch = delta_bsw_branch*0.1*cc%stored_N/delta_nsw_branch
           delta_nsw_branch = 0.1*cc%stored_N
        endif
        cc%bsw = cc%bsw+delta_bsw_branch
        cc%sapwood_N = cc%sapwood_N+delta_nsw_branch
        cc%branch_sw_loss = cc%branch_sw_loss - delta_bsw_branch
        cc%stored_N = cc%stored_N - delta_nsw_branch
     endif

     ! replace lost branch wood from sapwood pool- do we need respiration here?
     if (cc%branch_wood_loss > 0) then
        delta_wood_branch = max(min(cc%branch_wood_loss, 0.25*cc%bsw), 0.0)
        cc%bsw = cc%bsw - delta_wood_branch
        cc%bwood = cc%bwood + delta_wood_branch
        cc%wood_N = cc%wood_N + delta_wood_branch/sp%wood_c2n
        cc%sapwood_N = cc%wood_N - delta_wood_branch/sp%sapwood_c2n
        cc%stored_N = cc%stored_N - delta_wood_branch/sp%wood_c2n + delta_wood_branch/sp%sapwood_c2n
        cc%branch_wood_loss = cc%branch_wood_loss - delta_wood_branch
     endif
     ! in principle re-building wood and sapwood of branches should be taken into account together
     ! since they are the same process and should be limited by the upper limit of
     ! N storage spending rate together.


     G_WF=0.0
     if (NSCtarget < cc%nsc) then ! ens change this
        G_WF = max (0.0, fsf*(cc%nsc - NSCtarget)/(1+GROWTH_RESP))

     ! make sure that if nsc is not spent because N is limited, sapwood is limited too
     endif
     ! change maturity threashold to a diameter threash-hold
     if (cc%layer == 1 .AND. cc%age > sp%maturalage) then
        ! deltaSeed=      sp%v_seed * (cc%carbon_gain - G_LFR)
        ! deltaBSW = (1.0-sp%v_seed)* (cc%carbon_gain - G_LFR)
        if (cc%nitrogen_stress > sp%max_n_stress_for_seed_production.and.N_limits_live_biomass) then
           deltaSeed = 0
        else
           deltaSeed = sp%v_seed * G_WF
        endif
        deltaBSW = G_WF - deltaSeed
     else
        deltaSeed= 0.0
        deltaBSW = G_WF
     endif
     if (deltaBSW/sp%sapwood_c2n>0.1*cc%stored_N.and.N_limits_live_biomass) then
         deltaBSW = 0.1*cc%stored_N*sp%sapwood_c2n
     endif
     ! update biomass pools due to growth
     cc%bl     = cc%bl    + deltaBL;   cc%leaf_N = cc%leaf_N + deltaNL ! updated in vegn_int_ppa
     cc%br     = cc%br    + deltaBR;   cc%root_N = cc%root_N + deltaNR
     cc%bsw    = cc%bsw   + deltaBSW;  cc%sapwood_N = cc%sapwood_N + deltaBSW/sp%sapwood_c2n
     cc%bseed  = cc%bseed + deltaSeed; cc%seed_N = cc%seed_N + deltaSeed/sp%seed_c2n ! TODO: calculate seed_c2n as derived quantity in vegn_data initialization
     cc%nsc    = cc%nsc - (deltaBL + deltaBR + deltaSeed + deltaBSW + delta_bsw_branch)*(1.+GROWTH_RESP)
     cc%stored_N = cc%stored_N - deltaBSW/sp%sapwood_c2n - deltaSeed/sp%seed_c2n

     wood_prod = deltaBSW*365.0 ! conversion from kgC/day to kgC/year
     leaf_root_gr = G_LFR*365.0 ! conversion from kgC/day to kgC/year
     sw_seed_gr = (G_WF+delta_bsw_branch )*GROWTH_RESP*days_per_year ! conversion from kgC/day to kgC/year

     call check_var_range(cc%nsc,0.0,HUGE(1.0),'vegn_dynamics','cc%nsc', FATAL)
!     call check_var_range(cc%stored_N,0.0,HUGE(1.0),'vegn_dynamics','cc%stored_N', FATAL)

     !ens --compute daily growth to compute respiration, apply it next day, use npp_previous day variable, units kg C/(m2 *year)
     cc%growth_previous_day = cc%growth_previous_day+(max(0., G_LFR+G_WF)+delta_bsw_branch)*GROWTH_RESP ! this is for growth respiration to come from nsc

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
!         __DEBUG4__(cc%bl,cc%br,cc%bsw,cc%bwood)
!         __DEBUG3__(cc%DBH,cc%height,cc%crownarea)
!         __DEBUG4__(deltaBL,deltaBR,deltaBSW,deltaSeed)
!         __DEBUG3__(deltaDBH,deltaHeight,deltaCA)
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
     CSAsw    = sp%alphaCSASW * cc%DBH**sp%thetaCSASW ! sapwood cross-section
     CSAtot   = PI * (cc%DBH/2.0)**2 ! total trunk cross-section
     CSAwd    = max(0.0, CSAtot - CSAsw) ! cross-section of heartwood
     DBHwd    = 2*sqrt(CSAwd/PI) ! DBH of heartwood
     select case(sp%allomt)
     case (ALLOM_EW,ALLOM_EW1)
        BSWmax = sp%alphaBM * sp%rho_wood * (cc%DBH**sp%thetaBM - DBHwd**sp%thetaBM)
     case (ALLOM_HML)
        BSWmax = sp%alphaBM * sp%rho_wood * cc%height * (cc%DBH**2 - DBHwd**2)
     end select

     ! Update Kxa, stem conductance if we are tracking past damage
     ! TODO: make hydraulics_repair a namelist parameter?
     if (.not.hydraulics_repair) then
        deltaCSAsw = CSAsw - (sp%alphaCSASW * (cc%DBH - deltaDBH)**sp%thetaCSASW)
        cc%Kxa = (cc%Kxa*CSAsw + sp%Kxam*deltaCSAsw)/(CSAsw + deltaCSAsw)
     endif

     ! slm 20160523: are we retiring sapwood to wood only if it exceeds max sapwood? why?
     deltaBwood = max(cc%bsw - BSWmax, 0.0)
     cc%bwood   = cc%bwood + deltaBwood
     cc%bsw     = cc%bsw   - deltaBwood
     cc%sapwood_N = cc%sapwood_N - deltaBwood/sp%sapwood_c2n
     cc%wood_N    = cc%wood_N    + deltaBwood/sp%wood_c2n
     cc%stored_N  = cc%stored_N  + deltaBwood/sp%sapwood_c2n - deltaBwood/sp%wood_c2n

     ! update bl_max and br_max daily
     ! slm: why are we updating topyear only when the leaves are displayed? The paper
     !      never mentions this fact (see eq. A6).
     BL_u = sp%LMA * sp%LAImax * cc%crownarea * (1.0-sp%internal_gap_frac) * understory_lai_factor
     BL_c = sp%LMA * sp%LAImax * cc%crownarea * (1.0-sp%internal_gap_frac)
     if (cc%layer > 1 .and. cc%firstlayer == 0) then ! changed back, Weng 2014-01-23
        cc%topyear = 0.0
        cc%bl_max = BL_u
        cc%br_max = 0.8*cc%bl_max/(sp%LMA*sp%SRA)
     else
        if(cc%layer == 1)cc%topyear = cc%topyear + 1.0/365.0  ! daily step
        if(cc%layer>1)cc%firstlayer = 0 ! Just for the first year, those who were
                  ! pushed to understory have the characteristics of canopy trees
        if(sp%lifeform == FORM_GRASS) then
           cc%bl_max = BL_c
        else
           cc%bl_max = BL_u + min(cc%topyear/5.0,1.0)*(BL_c - BL_u)
        endif
        cc%br_max = sp%phiRL*cc%bl_max/(sp%LMA*sp%SRA)
     endif

     if(is_watch_point()) then
        __DEBUG2__(cc%bl_max, cc%br_max)
        __DEBUG2__(cc%carbon_gain, G_LFR)
        __DEBUG5__(deltaBL, deltaBR, deltaBSW, deltaSeed, deltaBwood)
        __DEBUG5__(cc%bl, cc%br, cc%bsw, cc%bseed, cc%bwood)
     endif
  else  ! cc%status == LEAF_OFF
    ! ens, after hsc was updated after carbon gain
    ! cc%nsc = cc%nsc + cc%carbon_gain
    !should some nsc go into sapwood and wood or seed

  endif ! cc%status == LEAF_ON

  ! calculate spending rate of growth respiration, to distribute it uniformly
  ! in time over the next day:
  cc%growth_previous_day_tmp = max(0.0,cc%growth_previous_day)*365.0
  ! factor 365.0 converts the rate of growth respiration release to atmosphere
  ! from kgC/day (frequency of this subroutine calls) to kgC/year, the units we
  ! use for other vegetation fluxes

  ! reset carbon accumulation terms
  cc%carbon_gain = 0
  cc%carbon_loss = 0

  end associate ! F2003

end subroutine biomass_allocation_ppa


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
subroutine vegn_phenology_ppa(vegn, soil)
  type(vegn_tile_type), intent(inout) :: vegn
  type(soil_tile_type), intent(inout) :: soil

  ! ---- local vars
  integer :: i
  real    :: leaf_litter_C, leaf_litter_N, leaf_litt_C(N_C_TYPES), leaf_litt_N(N_C_TYPES)
  real    :: dead_leaves_C, dead_leaves_N
  real    :: dead_roots_C, dead_roots_N

  real, parameter :: leaf_fall_rate = 0.075 ! per day
  real, parameter :: root_mort_rate = 0.0

  vegn%litter = 0 ; leaf_litt_C(:) = 0.0 ; leaf_litt_N(:) = 0.0
  do i = 1,vegn%n_cohorts
     associate ( cc => vegn%cohorts(i),   &
                 sp => spdata(vegn%cohorts(i)%species) )

     ! if drought-deciduous or cold-deciduous species
     ! temp=10 degrees C gives good growing season pattern
     ! spp=0 is c3 grass,1 c3 grass,2 deciduous, 3 evergreen
     ! assumption is that no difference between drought and cold deciduous

!    onset of phenology
     if(cc%status==LEAF_OFF .and. cc%gdd>sp%gdd_crit .and. vegn%tc_pheno>sp%tc_crit ) then
        cc%status = LEAF_ON
     else if ( cc%status == LEAF_ON .and. vegn%tc_pheno < sp%tc_crit ) then
        cc%status = LEAF_OFF
        cc%gdd = 0.0
     endif

!    leaf falling at the end of a growing season
! slm: why leaf fall and root mortality are computed from max bl and br, not
!      not actual bl and br?
     if(cc%status == LEAF_OFF .AND. cc%bl > 0.)then
         dead_leaves_C = min(leaf_fall_rate * cc%bl_max, cc%bl)
         dead_leaves_N = cc%leaf_N * dead_leaves_C/cc%bl
         dead_roots_C = min(root_mort_rate * cc%br_max, cc%br)
         dead_roots_N = cc%root_N * dead_roots_C/cc%br
         ! update C and N pools
         cc%nsc      = cc%nsc      + l_fract * (dead_leaves_C+dead_roots_C)
         cc%stored_N = cc%stored_N + l_fract * (dead_leaves_N+dead_roots_N)
         cc%bl       = cc%bl - dead_leaves_C
         cc%leaf_N   = cc%leaf_N - dead_leaves_N
         cc%br       = cc%br - dead_roots_C
         cc%root_N   = cc%root_N - dead_roots_N

         cc%lai = leaf_area_from_biomass(cc%bl,cc%species,cc%layer,cc%firstlayer)/ &
                                        (cc%crownarea *(1.0-sp%internal_gap_frac))
         if(cc%bl == 0.)cc%leaf_age = 0.0
         cc%bliving = cc%blv + cc%br + cc%bl + cc%bsw

         leaf_litter_C = (1-l_fract) * (dead_leaves_C+dead_roots_C) * cc%nindivs
         leaf_litter_N = (1-l_fract) * (dead_leaves_N+dead_roots_N) * cc%nindivs
         leaf_litt_C(:) = leaf_litt_C(:)+[sp%fsc_liv,1-sp%fsc_liv,0.0]*leaf_litter_C
         leaf_litt_N(:) = leaf_litt_N(:)+[sp%fsc_liv,1-sp%fsc_liv,0.0]*leaf_litter_N
         vegn%litter = vegn%litter + leaf_litter_C
         soil%fsc_in(1)  = soil%fsc_in(1) + leaf_litter_C
         vegn%veg_out = vegn%veg_out + leaf_litter_C
     endif
     end associate
  enddo
  ! add litter accumulated over the cohorts
  call add_soil_carbon(soil, vegn, leaf_litter_C=leaf_litt_C, leaf_litter_N=leaf_litt_N)
end subroutine vegn_phenology_ppa


! ===========================================================================
! leaf falling at LEAF_OFF -- it is unused; why is it here? is there anything
! in new Ensheng's code that uses it?
subroutine vegn_leaf_fall_ppa(vegn,soil)
  type(vegn_tile_type), intent(inout) :: vegn
  type(soil_tile_type), intent(inout) :: soil

  ! ---- local vars
  integer :: i
  real    :: leaf_litt_C(N_C_TYPES), leaf_litt_N(N_C_TYPES)
  real    :: dead_leaves_C, dead_leaves_N

  real, parameter :: leaf_fall_rate = 0.075 ! per day

  vegn%litter = 0 ! slm: why is it set to 0 here?
  leaf_litt_C(:) = 0.0; leaf_litt_N(:) = 0.0
  do i = 1,vegn%n_cohorts
     associate ( cc => vegn%cohorts(i),   &
                 sp => spdata(vegn%cohorts(i)%species) )
     if(cc%status == LEAF_OFF)then
        cc%nsc = cc%nsc + cc%carbon_gain
        dead_leaves_C = MIN(leaf_fall_rate * cc%bl_max, cc%bl)
        dead_leaves_N = cc%leaf_N * dead_leaves_C/cc%bl
        cc%nsc = cc%nsc + l_fract * dead_leaves_C
        cc%stored_N = cc%stored_N + l_fract * dead_leaves_N
        cc%bl  = cc%bl - dead_leaves_C ; cc%leaf_N  = cc%leaf_N - dead_leaves_N
        cc%lai = leaf_area_from_biomass(cc%bl,cc%species,cc%layer,cc%firstlayer)/ &
                                       (cc%crownarea *(1.0-sp%internal_gap_frac))
        if(cc%bl == 0.) cc%leaf_age = 0.0
        cc%bliving = cc%blv + cc%br + cc%bl + cc%bsw

        leaf_litt_C(:) = leaf_litt_C(:) + [sp%fsc_liv,1-sp%fsc_liv,0.0]*(1-l_fract) * dead_leaves_C * cc%nindivs
        leaf_litt_N(:) = leaf_litt_N(:) + [sp%fsc_liv,1-sp%fsc_liv,0.0]*(1-l_fract) * dead_leaves_N * cc%nindivs
     endif
     end associate
  enddo
  vegn%litter  = vegn%litter  + sum(leaf_litt_C)
  vegn%veg_out = vegn%veg_out + sum(leaf_litt_C)
  call add_soil_carbon(soil, vegn, leaf_litter_C=leaf_litt_C, leaf_litter_N=leaf_litt_N)
end subroutine vegn_leaf_fall_ppa


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
! the reproduction of each canopy cohort, yearly time step
! calculate the new cohorts added in this step and states:
! tree density, DBH, woody and living biomass
subroutine vegn_reproduction_ppa (vegn,soil)
  type(vegn_tile_type), intent(inout) :: vegn
  type(soil_tile_type), intent(inout) :: soil

! ---- local vars
  type(vegn_cohort_type), pointer :: ccold(:)   ! pointer to old cohort array
  real :: failed_seeds !, prob_g, prob_e
  integer :: newcohorts ! number of new cohorts to be created
  integer :: i, k ! cohort indices
  real :: litt(N_C_TYPES)

!  write(*,*)'vegn_reproduction_ppa n_cohorts before: ', vegn%n_cohorts

! Check if reproduction happens
  newcohorts = 0
  do i = 1, vegn%n_cohorts
     if (cohort_can_reproduce(vegn%cohorts(i))) newcohorts=newcohorts + 1
  enddo

  if (newcohorts == 0) return ! do nothing if no cohorts are ready for reproduction

!  write(*,*) 'creating new cohorts: ', newcohorts

  ccold => vegn%cohorts ! keep old cohort information
  allocate(vegn%cohorts(vegn%n_cohorts+newcohorts))
  vegn%cohorts(1:vegn%n_cohorts) = ccold(1:vegn%n_cohorts) ! copy old cohort information
  deallocate (ccold)

  litt(:) = 0.0
  ! set up new cohorts
  k = vegn%n_cohorts
  do i = 1,vegn%n_cohorts
    if (.not.cohort_can_reproduce(vegn%cohorts(i))) cycle ! nothing to do

    k = k+1 ! increment new cohort index

    ! update child cohort parameters
    associate (cc     => vegn%cohorts(k), &  ! F2003
               parent => vegn%cohorts(i), &
               sp     => spdata(vegn%cohorts(i)%species) )
    ! copy all parent values into the new cohort then change some of them below
    cc            = parent
    cc%status     = LEAF_OFF
    cc%firstlayer = 0
    cc%age        = 0.0
    cc%topyear    = 0.0
    call init_cohort_allometry_ppa(cc, sp%seedling_height, sp%seedling_nsc_frac)

    ! added germination probability (prob_g) and establishment probability ((prob_e), Weng 2014-01-06
    cc%nindivs = parent%bseed*parent%nindivs * sp%prob_g * sp%prob_e   &
                 /biomass_of_individual(cc)
!    __DEBUG3__(cc%age, cc%layer, cc%nindivs)

    failed_seeds = (1.-sp%prob_g*sp%prob_e) * parent%bseed * parent%nindivs
    vegn%litter = vegn%litter + failed_seeds
    litt(:) = litt(:) + [sp%fsc_liv,1-sp%fsc_liv,0.0]*failed_seeds
    vegn%veg_out = vegn%veg_out + failed_seeds

    parent%bseed = 0.0

    ! we assume that the newborn cohort is dry; since nindivs of the parent
    ! doesn't change we don't need to do anything with its Wl and Ws to
    ! conserve water (since Wl and Ws are per individual)
    cc%Wl = 0 ; cc%Ws = 0
    ! TODO: make sure that energy is conserved in reproduction
    cc%Tv = parent%Tv

    end associate   ! F2003
  enddo
  ! FIXME slm: add nitrogen to seedlings, and respectively to the litter
  call add_soil_carbon(soil, vegn, leaf_litter_C=litt)

  vegn%n_cohorts = k
  if(is_watch_point()) then
     write(*,*)'##### vegn_reproduction_ppa #####'
     __DEBUG2__(newcohorts, vegn%n_cohorts)
     __DEBUG1__(vegn%cohorts%nindivs)
     __DEBUG1__(vegn%cohorts%Tv)
  endif
end subroutine vegn_reproduction_ppa


! ============================================================================
subroutine kill_small_cohorts_ppa(vegn,soil)
  type(vegn_tile_type), intent(inout) :: vegn
  type(soil_tile_type), intent(inout) :: soil

  ! ---- local vars
  real, parameter :: mindensity = 1.0E-6

  type(vegn_cohort_type), pointer :: cc(:) ! array to hold new cohorts
  real, dimension(N_C_TYPES) :: &
     leaf_litt_C, leaf_litt_N, & ! fine surface litter per tile, kgC/m2 and kgN/m2
     wood_litt_C, wood_litt_N    ! coarse surface litter per tile, kgC/m2 and kgN/m2
  real, dimension(num_l, N_C_TYPES) :: &
     root_litt_C, root_litt_N ! root litter per soil layer, kgC/m2 and kgN/m2
  integer :: i,k

!  write(*,*)'kill_small_cohorts_ppa n_cohorts before: ', vegn%n_cohorts

 ! Weng, 2013-09-07
 ! calculate the number of remaining cohorts
  k = 0
  do i = 1, vegn%n_cohorts
     if (vegn%cohorts(i)%nindivs >  mindensity) k=k+1
  enddo

  ! if (k==0) call error_mesg('vegn_mergecohorts_ppa','All cohorts died',WARNING)

  ! exclude cohorts that have zero individuals
  leaf_litt_C = 0 ; wood_litt_C = 0; root_litt_C = 0
  leaf_litt_N = 0 ; wood_litt_N = 0; root_litt_N = 0
  if (k < vegn%n_cohorts) then
     allocate(cc(max(k,1)))
     k=0
     do i = 1,vegn%n_cohorts
        if (vegn%cohorts(i)%nindivs > mindensity) then
           k=k+1
           cc(k) = vegn%cohorts(i)
        else
           call kill_plants_ppa(vegn%cohorts(i), vegn, soil, vegn%cohorts(i)%nindivs, 0.0, &
                                leaf_litt_C, wood_litt_C, root_litt_C, &
                                leaf_litt_N, wood_litt_N, root_litt_N  )
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

  if (is_watch_point()) then
     write(*,*) '##### vegn_mergecohorts_ppa output #####'
     __DEBUG1__(vegn%n_cohorts)
     __DEBUG1__(vegn%cohorts%nindivs)
     __DEBUG1__(vegn%cohorts%Wl)
     __DEBUG1__(vegn%cohorts%Ws)
     __DEBUG1__(vegn%cohorts%mcv_dry)
     __DEBUG1__(vegn%cohorts%Tv)
  endif
!  write(*,*)'kill_small_cohorts_ppa n_cohorts after: ', vegn%n_cohorts

end subroutine kill_small_cohorts_ppa

end module vegn_dynamics_mod
