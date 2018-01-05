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
     do_check_conservation, carbon_cons_tol, set_current_point, land_error_message
use land_tile_mod, only : land_tile_map, land_tile_type, land_tile_enum_type, &
     first_elmt, loop_over_tiles, land_tile_carbon
use land_tile_diag_mod, only : OP_SUM, OP_AVERAGE, cmor_name, &
     register_tiled_diag_field, add_tiled_diag_field_alias, send_tile_data, diag_buff_type, &
     register_cohort_diag_field, send_cohort_data, set_default_diag_filter
use vegn_data_mod, only : spdata, nspecies, do_ppa, &
     PHEN_DECIDUOUS, PHEN_EVERGREEN, LEAF_ON, LEAF_OFF, FORM_WOODY, FORM_GRASS, &
     ALLOM_EW, ALLOM_EW1, ALLOM_HML, LU_CROP, &
     fsc_liv, fsc_wood, fsc_froot, agf_bs, l_fract, understory_lai_factor, min_lai, &
     use_light_saber, nsc_starv_frac
use vegn_tile_mod, only: vegn_tile_type, vegn_tile_carbon, vegn_relayer_cohorts_ppa, &
     vegn_mergecohorts_ppa
use soil_tile_mod, only: num_l, dz, soil_tile_type, LEAF, CWOOD, N_LITTER_POOLS
use vegn_cohort_mod, only : vegn_cohort_type, &
     update_biomass_pools, update_bio_living_fraction, update_species, &
     leaf_area_from_biomass, cohort_root_litter_profile, cohort_root_exudate_profile
use vegn_util_mod, only : kill_plants_ppa, add_seedlings_ppa
use vegn_harvesting_mod, only : allow_weeds_on_crops
use soil_carbon_mod, only: N_C_TYPES, soil_carbon_option, &
    SOILC_CENTURY, SOILC_CENTURY_BY_LAYER, SOILC_CORPSE, &
    add_litter
use soil_mod, only: Dsdt
use soil_util_mod, only: add_soil_carbon, add_root_litter, add_root_exudates

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
integer :: id_resp, id_resl, id_resr, id_ress, id_resg
integer :: id_soilt, id_theta, id_litter, id_age, id_exudate, id_dbh_growth
! CMIP/CMOR diagnostic field IDs
integer :: id_gpp_cmor, id_npp_cmor, id_nep_cmor, id_ra, id_rgrowth

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

  ! set the default sub-sampling filter for CMOR variables
  call set_default_diag_filter('land')
  id_gpp_cmor = register_tiled_diag_field ( cmor_name, 'gpp', (/id_ug/), &
       time, 'Carbon Mass Flux out of Atmosphere due to Gross Primary Production on Land', &
       'kg m-2 s-1', missing_value=-1.0, fill_missing=.TRUE., &
       standard_name='gross_primary_productivity_of_biomass_expressed_as_carbon')
  call add_tiled_diag_field_alias ( id_gpp_cmor, cmor_name, 'gppLut', (/id_ug/), &
       time, 'Gross Primary Productivity on Land Use Tile', 'kg m-2 s-1', missing_value=-1.0, &
       standard_name='gross_primary_productivity_of_biomass_expressed_as_carbon', fill_missing=.FALSE.)
  id_npp_cmor = register_tiled_diag_field ( cmor_name, 'npp', (/id_ug/), &
       time, 'Carbon Mass Flux out of Atmosphere due to Net Primary Production on Land', &
       'kg m-2 s-1', missing_value=-1.0, fill_missing=.TRUE., &
       standard_name='net_primary_productivity_of_biomass_expressed_as_carbon')
  call add_tiled_diag_field_alias ( id_npp_cmor, cmor_name, 'nppLut', (/id_ug/), &
       time, 'Carbon Mass Flux out of Atmosphere due to Net Primary Production on Land', &
       'kg m-2 s-1', missing_value=-1.0, fill_missing=.FALSE., &
       standard_name='net_primary_productivity_of_biomass_expressed_as_carbon')
  id_nep_cmor = register_tiled_diag_field ( cmor_name, 'nep', (/id_ug/), &
       time, 'Net Carbon Mass Flux out of Atmophere due to Net Ecosystem Productivity on Land.', &
       'kg m-2 s-1', missing_value=-1.0, &
       standard_name='surface_net_downward_mass_flux_of_carbon_dioxide_expressed_as_carbon_due_to_all_land_processes_excluding_anthropogenic_land_use_change', &
       fill_missing=.TRUE.)
  id_ra = register_tiled_diag_field ( cmor_name, 'ra', (/id_ug/), &
       time, 'Carbon Mass Flux into Atmosphere due to Autotrophic (Plant) Respiration on Land', &
       'kg m-2 s-1', missing_value=-1.0, &
       standard_name='plant_respiration_carbon_flux', fill_missing=.TRUE.)
  call add_tiled_diag_field_alias (id_ra, cmor_name, 'raLut', (/id_ug/), &
       time, 'Carbon Mass Flux into Atmosphere due to Autotrophic (Plant) Respiration on Land', &
       'kg m-2 s-1', missing_value=-1.0, &
       standard_name='plant_respiration_carbon_flux', fill_missing=.FALSE.)
  id_rgrowth = register_tiled_diag_field ( cmor_name, 'rGrowth', (/id_ug/), &
       time, 'Growth Autotrophic Respiration', 'kg m-2 s-1', missing_value=-1.0, &
       standard_name='surface_upward_carbon_mass_flux_due_to_plant_respiration_for_biomass_growth',&
       fill_missing=.TRUE.)
end subroutine vegn_dynamics_init

! ============================================================================
subroutine vegn_dynamics_end()
   deallocate(sg_soilfrac, ug_area_factor)
end subroutine vegn_dynamics_end


! ============================================================================
subroutine vegn_carbon_int_lm3(vegn, soil, soilt, theta, diag)
  type(vegn_tile_type), intent(inout) :: vegn
  type(soil_tile_type), intent(inout) :: soil
  real, intent(in) :: soilt ! average temperature of soil for soil carbon decomposition, deg K
  real, intent(in) :: theta ! average soil wetness, unitless
  type(diag_buff_type), intent(inout) :: diag

  ! TODO: possibly move soil-related calculations from calling procedure here,
  !       now that we have soil passed as an argument

  type(vegn_cohort_type), pointer :: cc
  type(vegn_cohort_type), pointer :: c(:) ! for debug only
  real :: md_leaf, md_wood, md_froot ! component of maintenance demand
  real :: md ! plant tissue maintenance, kg C/timestep
  real :: root_exudate_C       ! root exudate, kgC/(year indiv)
  real :: total_root_exudate_C(num_l) ! total root exudate per tile, kgC/m2
  real :: leaf_litt(n_c_types) ! fine surface litter per tile, kgC/m2
  real :: wood_litt(n_c_types) ! coarse surface litter per tile, kgC/m2
  real :: root_litt(num_l, n_c_types) ! root litter per soil layer, kgC/m2
  real :: profile(num_l) ! storage for vertical profile of exudates and root litter
  real, dimension(vegn%n_cohorts) :: resp, resl, resr, ress, resg, gpp, npp
  integer :: sp ! shorthand for current cohort specie
  integer :: i, l, N

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
  leaf_litt = 0 ; wood_litt = 0; root_litt = 0 ; total_root_exudate_C = 0
  do i = 1, vegn%n_cohorts
     cc => vegn%cohorts(i)
     sp = cc%species

     call plant_respiration(cc,soilt, resp(i), resl(i), resr(i), ress(i))

     gpp(i) = (cc%An_op - cc%An_cl)*mol_C*cc%leafarea;
     npp(i) = gpp(i) - resp(i)

     if(cc%npp_previous_day > 0) then
        resg(i) = spdata(sp)%GROWTH_RESP*cc%npp_previous_day
        npp(i)  = npp(i)  - spdata(sp)%GROWTH_RESP*cc%npp_previous_day
        resp(i) = resp(i) + spdata(sp)%GROWTH_RESP*cc%npp_previous_day
     else
        resg(i) = 0
     endif
     ! accumulate npp for the current day
     cc%npp_previous_day_tmp = cc%npp_previous_day_tmp + npp(i);

     root_exudate_C = max(npp(i),0.0)*spdata(sp)%root_exudate_frac
     call cohort_root_exudate_profile(cc,dz,profile)
     total_root_exudate_C = total_root_exudate_C + profile*root_exudate_C*cc%nindivs*dt_fast_yr
     cc%carbon_gain = cc%carbon_gain + (npp(i)-root_exudate_C)*dt_fast_yr

     ! check if leaves/roots are present and need to be accounted in maintenance
     if(cc%status == LEAF_ON) then
        md_leaf = cc%Pl * spdata(sp)%alpha_leaf*cc%bliving*dt_fast_yr
        md_froot= cc%Pr * spdata(sp)%alpha_root*cc%bliving*dt_fast_yr
     else
        md_leaf  = 0
        md_froot = 0
     endif

     ! compute branch and coarse wood losses for tree types
     if (spdata(sp)%lifeform==FORM_WOODY) then
        md_wood = 0.6 * cc%bwood * spdata(sp)%alpha_wood*dt_fast_yr
     else
        md_wood = 0
     endif

     md = md_leaf + md_froot + cc%Psw_alphasw * cc%bliving * dt_fast_yr

     cc%bwood_gain = cc%bwood_gain + cc%Psw_alphasw * cc%bliving * dt_fast_yr;
     cc%bwood_gain = cc%bwood_gain - md_wood;
!     if (cc%bwood_gain < 0.0) cc%bwood_gain=0.0; ! potential non-conservation ?
     cc%carbon_gain = cc%carbon_gain - md;
     cc%carbon_loss = cc%carbon_loss + md; ! used in diagnostics only

     ! add maintenance demand from leaf and root pools to fast soil carbon
     leaf_litt(:) = leaf_litt(:) + (/fsc_liv,  1-fsc_liv,  0.0/)*md_leaf*cc%nindivs
     wood_litt(:) = wood_litt(:) + (/fsc_wood, 1-fsc_wood, 0.0/)*md_wood*agf_bs*cc%nindivs
     call cohort_root_litter_profile(cc, dz, profile)
     do l = 1, num_l
        root_litt(l,:) = root_litt(l,:) + profile(l)*cc%nindivs*(/ &
             fsc_froot    *md_froot + fsc_wood    *md_wood*(1-agf_bs), &
             (1-fsc_froot)*md_froot + (1-fsc_wood)*md_wood*(1-agf_bs), &
             0.0/)
     enddo

     vegn%veg_in  = vegn%veg_in  + npp(i)*cc%nindivs*dt_fast_yr;
     vegn%veg_out = vegn%veg_out + (md_leaf+md_froot+md_wood)*cc%nindivs;

     ! update cohort age
     cc%age = cc%age + dt_fast_yr
  enddo

  ! fsc_in and ssc_in updated in add_root_exudates
  call add_root_exudates(soil,total_root_exudate_C)

  ! add litter accumulated over the cohorts
  call add_soil_carbon(soil, leaf_litt, wood_litt, root_litt)

  if(is_watch_point()) then
     write(*,*)'#### vegn_carbon_int_lm3 output ####'
     __DEBUG1__(c%species)
     __DEBUG1__(c%bl)
     __DEBUG1__(c%br)
     __DEBUG1__(c%bsw)
     __DEBUG1__(c%bwood)
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
  endif

  ! update soil carbon
  call Dsdt(vegn, soil, diag, soilt, theta)

  ! NEP is equal to NNP minus soil respiration
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

  ! ---- CMOR diagnostics
  if (id_gpp_cmor>0) call send_tile_data(id_gpp_cmor, sum(gpp(1:N)*c(1:N)%nindivs)/seconds_per_year, diag)
  if (id_npp_cmor>0) call send_tile_data(id_npp_cmor, sum(npp(1:N)*c(1:N)%nindivs)/seconds_per_year, diag)
  if (id_nep_cmor>0) call send_tile_data(id_nep_cmor, vegn%nep/seconds_per_year, diag)
  if (id_ra>0) call send_tile_data(id_ra, sum((resp(1:N)-resg(1:N))*c(1:N)%nindivs)/seconds_per_year, diag)
  if (id_rgrowth>0) call send_tile_data(id_rgrowth, sum(resg(1:N)*c(1:N)%nindivs)/seconds_per_year, diag)

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
  real :: deltaBL, deltaBR ! leaf and fine root carbon tendencies
  integer :: i, l, N
  real :: NSC_supply,LR_demand,LR_deficit
  real :: R_days,fNSC,fLFR,fStem,fSeed
  type(vegn_cohort_type), pointer :: c(:) ! for debug only
  ! accumulators of total input to root litter and soil carbon
  real :: leaf_litt(N_C_TYPES) ! fine surface litter per tile, kgC/m2
  real :: wood_litt(N_C_TYPES) ! coarse surface litter per tile, kgC/m2
  real :: root_litt(num_l, N_C_TYPES) ! root litter per soil layer, kgC/m2
  real :: profile(num_l) ! storage for vertical profile of exudates and root litter
  real, dimension(vegn%n_cohorts) :: resp, resl, resr, ress, resg, gpp, npp
  real, dimension(vegn%n_cohorts) :: root_exudate_C ! root exudate, kgC/(year indiv)
  real :: total_root_exudate_C(num_l) ! total root exudate per tile, kgC/m2


  c=>vegn%cohorts(1:vegn%n_cohorts)
  N = vegn%n_cohorts

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
  leaf_litt = 0 ; wood_litt = 0; root_litt = 0 ; total_root_exudate_C = 0
  ! TODO: add root exudates to the balance. in LM3 it's a fraction of npp;
  !       in PPA perhaps it might be proportional to NSC, with some decay
  !       time?
  do i = 1, vegn%n_cohorts
     associate ( cc => vegn%cohorts(i), sp => spdata(vegn%cohorts(i)%species))

     ! that was in eddy_npp_PPA
     call plant_respiration(cc,tsoil,resp(i),resl(i),resr(i),ress(i))

     gpp(i)  = (cc%An_op - cc%An_cl)*mol_C*cc%leafarea
     ! growth respiration comes from nsc pool not gpp now
     resg(i) = cc%growth_previous_day_tmp ! kg C per individual per year

     cc%growth_previous_day = cc%growth_previous_day - resg(i)*dt_fast_yr

     if (is_watch_point()) then
        write(*,'("vegn_carbon_int_ppa #0, cohort ",i2.2)',advance='NO') i
        __DEBUG2__(cc%nsc,cc%carbon_gain)
     endif
     cc%carbon_gain  = cc%carbon_gain + gpp(i)*dt_fast_yr - resp(i)*dt_fast_yr
     if (is_watch_point()) then
        write(*,'("vegn_carbon_int_ppa #1, cohort ",i2.2)',advance='NO') i
        __DEBUG2__(cc%nsc,cc%carbon_gain)
     endif

     resp(i) = resp(i) + resg(i) ! add growth respiration and maintenance respiration
     npp(i)  = gpp(i) - resp(i)
     root_exudate_C(i) = max(npp(i),0.0)*sp%root_exudate_frac
     call cohort_root_exudate_profile(cc, dz, profile)
     do l = 1,num_l
        total_root_exudate_C(l) = total_root_exudate_C(l) + &
             profile(l)*root_exudate_C(i)*cc%nindivs*dt_fast_yr
     enddo
     npp(i) = npp(i) - root_exudate_C(i)

     cc%carbon_gain = cc%carbon_gain-root_exudate_C(i)*dt_fast_yr
     if (is_watch_point()) then
        write(*,'("vegn_carbon_int_ppa #2, cohort ",i2.2)',advance='NO') i
        __DEBUG2__(cc%nsc,cc%carbon_gain)
     endif
     !print *, ' npp ',  npp(i),' root_ex ',  root_exudate_C(i), 'cgain ', cc%carbon_gain

     ! Weng, 2013-01-28
     ! Turnover regardless of STATUS
     deltaBL = cc%bl * sp%alpha_leaf * dt_fast_yr
     deltaBR = cc%br * sp%alpha_root * dt_fast_yr

     cc%bl = cc%bl - deltaBL
     cc%br = cc%br - deltaBR

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

     !reduce nsc by the amount of root exudates during the date
     ! TODO: when merging with N code take exudates from NSC not carbon gained
     ! TODO: need different than a fraction of NPP formulation
     !cc%nsc = cc%nsc-root_exudate_C(i)*dt_fast_yr

     ! accumulate liter and soil carbon inputs across all cohorts
     leaf_litt(:) = leaf_litt(:) + (/fsc_liv,  1.-fsc_liv,  0.0/)*deltaBL*cc%nindivs
     wood_litt(:) = wood_litt(:) + (/fsc_wood, 1.-fsc_wood, 0.0/)*(md_branch_sw+md_bsw)*cc%nindivs

     call cohort_root_litter_profile(cc, dz, profile)
     do l = 1, num_l
        root_litt(l,:) = root_litt(l,:) + profile(l)*cc%nindivs*(/ &
             fsc_froot    *deltaBR , & ! fast
             (1-fsc_froot)*deltaBR , & ! slow
             0.0/) ! microbes
     enddo


     ! increment cohort age
     cc%age = cc%age + dt_fast_yr
     ! resg(i) = 0 ! slm: that doesn't make much sense to me... why?
     end associate
  enddo

  ! add litter and exudates accumulated over the cohorts
  call add_root_exudates(soil, total_root_exudate_C)
  call add_soil_carbon(soil, leaf_litt, wood_litt, root_litt)
  ! update soil carbon
  call Dsdt(vegn, soil, diag, tsoil, theta)

  ! NEP is equal to GPP minus autotrophic and soil respiration
!  vegn%nep = sum(npp(1:N)*c(1:N)%nindivs) - vegn%rh
  vegn%nep = sum((gpp(1:N)-resp(1:N))*c(1:N)%nindivs) - vegn%rh

  call update_soil_pools(vegn, soil)
  vegn%age = vegn%age + dt_fast_yr;

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
     __DEBUG1__(vegn%rh)
     write(*,*)'#### end of vegn_carbon_int_ppa output ####'
  endif

! ------ diagnostic section
  call send_cohort_data(id_gpp, diag, c(1:N), gpp(1:N), weight=c(1:N)%nindivs, op=OP_SUM)
  call send_cohort_data(id_npp, diag, c(1:N), npp(1:N), weight=c(1:N)%nindivs, op=OP_SUM)
  call send_tile_data(id_nep,vegn%nep,diag)
  call send_tile_data(id_litter,vegn%litter,diag)
  call send_cohort_data(id_resp, diag, c(1:N), resp(1:N), weight=c(1:N)%nindivs, op=OP_SUM)
  call send_cohort_data(id_resl, diag, c(1:N), resl(1:N), weight=c(1:N)%nindivs, op=OP_SUM)
  call send_cohort_data(id_resr, diag, c(1:N), resr(1:N), weight=c(1:N)%nindivs, op=OP_SUM)
  call send_cohort_data(id_ress, diag, c(1:N), ress(1:N), weight=c(1:N)%nindivs, op=OP_SUM)
  call send_cohort_data(id_resg, diag, c(1:N), resg(1:N), weight=c(1:N)%nindivs, op=OP_SUM)
  call send_tile_data(id_soilt,tsoil,diag)
  call send_tile_data(id_theta,theta,diag)
  call send_cohort_data(id_age, diag, c(1:N), c(1:N)%age, weight=c(1:N)%nindivs, op=OP_AVERAGE)
  call send_cohort_data(id_exudate, diag, c(1:N), root_exudate_C(1:N), weight=c(1:N)%nindivs, op=OP_AVERAGE)

  ! ---- CMOR diagnostics
  if (id_gpp_cmor>0) call send_tile_data(id_gpp_cmor, sum(gpp(1:N)*c(1:N)%nindivs)/seconds_per_year, diag)
  if (id_npp_cmor>0) call send_tile_data(id_npp_cmor, sum(npp(1:N)*c(1:N)%nindivs)/seconds_per_year, diag)
  if (id_nep_cmor>0) call send_tile_data(id_nep_cmor, vegn%nep/seconds_per_year, diag)
  if (id_ra>0) call send_tile_data(id_ra, sum((resp(1:N)-resg(1:N))*c(1:N)%nindivs)/seconds_per_year, diag)
  if (id_rgrowth>0) call send_tile_data(id_rgrowth, sum(resg(1:N)*c(1:N)%nindivs)/seconds_per_year, diag)

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

  if (is_watch_point()) then
     cmass1 = vegn_tile_carbon(vegn)
     write(*,*)"############## vegn_growth #################"
  endif
end subroutine vegn_growth

!***********************************************************************
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
  real :: leaf_litt(N_C_TYPES) ! fine surface litter per tile, kgC/m2
  real :: wood_litt(N_C_TYPES) ! coarse surface litter per tile, kgC/m2
  real :: root_litt(num_l, N_C_TYPES) ! root litter per soil layer, kgC/m2
  real :: profile(num_l) ! storage for vertical profile of exudates and root litter

!  write(*,*)'vegn_starvation_ppa n_cohorts before: ', vegn%n_cohorts

  leaf_litt = 0 ; wood_litt = 0; root_litt = 0
  do i = 1, vegn%n_cohorts
     associate ( cc => vegn%cohorts(i)   , &
                 sp => spdata(vegn%cohorts(i)%species)  )  ! F2003

    ! Mortality due to starvation
    if (cc%bsw<0 .or. cc%nsc < nsc_starv_frac*cc%bl_max) then
       deathrate = 1.0

       deadtrees = min(cc%nindivs*deathrate,cc%nindivs) ! individuals / m2
       ! kill starved plants and add dead C from leaf and root pools to soil carbon
       call kill_plants_ppa(cc, vegn, deadtrees, 0.0, leaf_litt, wood_litt, root_litt)

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
  call add_soil_carbon(soil, leaf_litt, wood_litt, root_litt)

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
  real :: G_LFR  ! amount of carbon spent on leaf and root growth
  real :: G_WF   ! amount of carbon spent on new sapwood growth and seeds
  real :: G_WF_max ! max spending rate
  real :: deltaSeed
  real :: deltaBL, deltaBR ! tendencies of leaf and root biomass, kgC/individual
  real :: deltaBSW ! tendency of sapwood biomass, kgC/individual
  real :: deltaBwood ! tendency of wood biomass, kgC/individual
  real :: deltaCA ! tendency of crown area, m2/individual
  real :: deltaHeight ! tendency of vegetation height
  real :: deltaCSAsw ! tendency of sapwood area
  real :: BL_c, BL_u ! canopy and understory target leaf biomasses, kgC/individual
  real :: NSCtarget ! target NSC storage
  real :: delta_bsw_branch
  real :: delta_wood_branch

  real, parameter :: DBHtp = 0.8 ! m

  !ens
  real :: fsf = 0.2 ! in Weng et al 205, page 2677 units are 0.2 day

  real, parameter :: deltaDBH_max = 0.01 ! max growth rate of the trunk, meters per day
  real, parameter :: deltaLAI_max = 1.0  ! max leaf growth rate, m2/(m2 day)
  real  :: nsctmp ! temporal nsc budget to avoid biomass loss, KgC/individual

  associate (sp => spdata(cc%species)) ! F2003

  ! TODO: what if carbon_gain is not 0, but leaves are OFF (marginal case? or
  ! typical in lm3?)
  !delta_wood_branch  = 0.
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
     ! calculate the carbon spent on growth of leaves and roots
     ! don't allow more than 0.1/(1+GROWTH_RESP) of nsc per day to spend
     G_LFR = max(0.0, 0.1*cc%nsc/(1+sp%GROWTH_RESP))
     if (use_light_saber) then
        if (cc%bl > 0 .and. cc%An_newleaf_daily <= 0) G_LFR = 0.0 ! do not grow more leaves if they would not increase An
     else
        ! do not grow more than bl_max, br_max of leaves and roots
        G_LFR = max(0.0, min(cc%bl_max+cc%br_max-cc%bl-cc%br, G_LFR))
     endif

     ! and distribute it between roots and leaves
     deltaBL  = min(G_LFR, max(0.0, &
          (G_LFR*cc%bl_max + cc%bl_max*cc%br - cc%br_max*cc%bl)/(cc%bl_max + cc%br_max) &
          ))
     deltaBR  = G_LFR - deltaBL
     ! calculate carbon spent on seeds and sapwood growth

     !ens first replace lost branch sapwood
     !delta_bsw_branch  = 0.0
     !if (cc%branch_sw_loss > 0) then
     !   delta_bsw_branch = max (min(cc%branch_sw_loss, 0.1*cc%nsc/(1+GROWTH_RESP)), 0.0)
     !   cc%bsw = cc%bsw+delta_bsw_branch
     !   cc%branch_sw_loss = cc%branch_sw_loss - delta_bsw_branch
     !endif

     ! replace lost branch wood from sapwood pool- do we need respiration here?
     !if (cc%branch_wood_loss > 0) then
     !   delta_wood_branch = max(min(cc%branch_wood_loss, 0.25*cc%bsw), 0.0)
     !   cc%bsw = cc%bsw - delta_wood_branch
     !   cc%bwood = cc%bwood + delta_wood_branch
     !   cc%branch_wood_loss = cc%branch_wood_loss - delta_wood_branch
     !endif

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
     cc%brsw = cc%brsw+delta_bsw_branch
     cc%bsw = cc%bsw+delta_bsw_branch

     call check_var_range(cc%bsw, 0.0,HUGE(1.0),'biomass_allocation_ppa #2', 'cc%bsw',FATAL)
     call check_var_range(cc%brsw,0.0,cc%bsw,   'biomass_allocation_ppa #2', 'cc%brsw',FATAL)

     !ens 02/14/17
     !update seed and sapwood biomass pools with the new growth (branches were updated above)
     NSCtarget = sp%NSC2targetbl*cc%bl_max
     G_WF = max (0.0, fsf*(cc%nsc - NSCtarget)/(1+sp%GROWTH_RESP))
     if (cohort_can_reproduce(cc)) then
        deltaSeed = sp%v_seed * G_WF
     else
        deltaSeed = 0.0
     endif
     deltaBSW = G_WF - deltaSeed

     if (sp%lifeform == FORM_GRASS) then ! isa 20170705
        ! 20170724 - new scheme
        nsctmp = cc%nsc
        NSCtarget = sp%NSC2targetbl*cc%bl_max
        G_WF = max(0.0, fsf*(nsctmp - NSCtarget)/(1+sp%GROWTH_RESP))
        ! note that it is only for HML allometry now
        G_WF_max = deltaDBH_max/((sp%gammaHT+cc%DBH**sp%thetaHT)**2/(sp%rho_wood * sp%alphaHT * sp%alphaBM * &
                          (cc%DBH**(1.+sp%thetaHT)*(2.*(sp%gammaHT+cc%DBH**sp%thetaHT)+sp%gammaHT*sp%thetaHT))))
        G_WF = min(G_WF, G_WF_max)
        if (cohort_can_reproduce(cc)) then
           deltaSeed = sp%v_seed * G_WF
        else
           deltaSeed = 0.0
        endif
        deltaBSW = G_WF - deltaSeed

        nsctmp = nsctmp - G_WF * (1+sp%GROWTH_RESP)

        ! Allocation to leaves and fine root growth
!         G_LFR = max(0.0, min( cc%bl_max - cc%bl + max( 0.0, cc%br_max - cc%br ),  &
!                               0.1*nsctmp/(1+sp%GROWTH_RESP)))
        ! don't allow more than 0.1/(1+GROWTH_RESP) of nsc per day to spend
        G_LFR = max(0.0, 0.1*nsctmp/(1+sp%GROWTH_RESP))
        if (use_light_saber) then
           if (cc%bl > 0 .and. cc%An_newleaf_daily <= 0) G_LFR = 0.0 ! do not grow more leaves if they would not increase An
           G_LFR = min(G_LFR, deltaLAI_max*sp%LMA*cc%crownarea)
        else
           G_LFR = max(0.0, min(cc%bl_max-cc%bl+max(0.0, cc%br_max-cc%br), G_LFR))
        endif

        if ( cc%br > cc%br_max ) then      ! roots larger than target
            deltaBL = G_LFR
        else                               ! standard scheme
            deltaBL  = min(G_LFR, max(0.0, &
                          (G_LFR*cc%bl_max + cc%bl_max*cc%br - cc%br_max*cc%bl)/(cc%bl_max + cc%br_max) ))
        end if
        deltaBR  = G_LFR - deltaBL
     endif

     ! update biomass pools due to growth
     cc%bl     = cc%bl    + deltaBL  ! updated in vegn_int_ppa
     cc%br     = cc%br    + deltaBR
     cc%bsw    = cc%bsw   + deltaBSW
     cc%bseed  = cc%bseed + deltaSeed
     ! reduce storage by  growth respiration
     cc%nsc = cc%nsc - (G_LFR + G_WF + delta_bsw_branch)*(1+sp%GROWTH_RESP)

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

     call check_var_range(cc%bsw,  0.0,HUGE(1.0),'biomass_allocation_ppa #4', 'cc%bsw',FATAL)
     call check_var_range(cc%brsw, 0.0,cc%bsw,   'biomass_allocation_ppa #4', 'cc%brsw',FATAL)

     ! update bl_max and br_max daily
     ! slm: why are we updating topyear only when the leaves are displayed? The paper
     !      never mentions this fact (see eq. A6).
     BL_c = sp%LMA * sp%LAImax * cc%crownarea * (1.0-sp%internal_gap_frac)
     BL_u = Bl_c * understory_lai_factor
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
  ! reset accumulated dAn/dLAI derivative
  cc%An_newleaf_daily = 0.0

  end associate ! F2003

end subroutine biomass_allocation_ppa


! ============================================================================
! The combined reduction in decomposition rate as a function of TEMP and MOIST
! Based on CENTURY Parton et al 1993 GBC 7(4):785-809 and Bolker's copy of
! CENTURY code
elemental function A_function(soilt, theta) result(A)
  real :: A                 ! return value, resulting reduction in decomposition rate
  real, intent(in) :: soilt ! effective temperature for soil carbon decomposition
  real, intent(in) :: theta

  real :: soil_temp; ! temperature of the soil, deg C
  real :: Td; ! rate multiplier due to temp
  real :: Wd; ! rate reduction due to moisture

  ! coefficients and terms used in temperature term
  real :: Topt,Tmax,t1,t2,tshl,tshr;

  soil_temp = soilt-273.16;

  ! EFFECT OF TEMPERATURE
  ! from Bolker's century code
  Tmax=45.0;
  if (soil_temp > Tmax) soil_temp = Tmax;
  Topt=35.0;
  tshr=0.2; tshl=2.63;
  t1=(Tmax-soil_temp)/(Tmax-Topt);
  t2=exp((tshr/tshl)*(1.-t1**tshl));
  Td=t1**tshr*t2;

  if (soil_temp > -10) Td=Td+0.05;
  if (Td > 1.) Td=1.;

  ! EFFECT OF MOISTURE
  ! Linn and Doran, 1984, Soil Sci. Amer. J. 48:1267-1272
  ! This differs from the Century Wd
  ! was modified by slm/ens based on the figures from the above paper
  !     (not the reported function)

  if(theta <= 0.3) then
     Wd = 0.2;
  else if(theta <= 0.6) then
     Wd = 0.2+0.8*(theta-0.3)/0.3;
  else
     Wd = exp(2.3*(0.6-theta));
  endif

  A = (Td*Wd); ! the combined (multiplicative) effect of temp and water
               ! on decomposition rates
end function A_function


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
  type(vegn_cohort_type), pointer :: cc
  real :: leaf_litter,root_litter;
  real :: theta_crit; ! critical ratio of average soil water to sat. water
  real :: psi_stress_crit ! critical soil-water-stress index
  real :: wilt ! ratio of wilting to saturated water content
  real :: leaf_litt(n_c_types) ! fine surface litter per tile, kgC/m2
  real :: root_litt(num_l,n_c_types) ! root litter per soil layer, kgC/m2
  real :: profile(num_l) ! storage for vertical profile of root litter
  integer :: i, l

  wilt = soil%w_wilt(1)/soil%pars%vwc_sat
  vegn%litter = 0

  leaf_litt = 0 ; root_litt = 0
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

           leaf_litter = (1.0-l_fract)*cc%bl*cc%nindivs;
           root_litter = (1.0-l_fract)*cc%br*cc%nindivs;
           leaf_litt(:) = leaf_litt(:)+(/fsc_liv*leaf_litter,(1-fsc_liv)*leaf_litter,0.0/)
           call cohort_root_litter_profile(cc, dz, profile)
           do l = 1, num_l
              root_litt(l,:) = root_litt(l,:) + profile(l)*(/ &
                   fsc_froot*root_litter, &
                   (1-fsc_froot)*root_litter, &
                   0.0/)
           enddo

           vegn%litter = vegn%litter + leaf_litter + root_litter
           vegn%veg_out = vegn%veg_out + leaf_litter+root_litter;

           cc%blv = cc%blv + l_fract*(cc%bl+cc%br);
           cc%bl  = 0.0;
           cc%br  = 0.0;
           cc%lai = 0.0;

           ! update state
           cc%bliving = cc%blv + cc%br + cc%bl + cc%bsw;
           call update_bio_living_fraction(cc);
        endif
     endif ! phenology type
     end associate
  enddo

  ! add litter accumulated over the cohorts
  call add_soil_carbon(soil, leaf_litter=leaf_litt, root_litter=root_litt)

end subroutine vegn_phenology_lm3


! =============================================================================
! Added by Weng 2012-02-29
subroutine vegn_phenology_ppa(tile)
  type(land_tile_type), intent(inout) :: tile

  ! ---- local vars
  integer :: i
  real    :: leaf_litter, leaf_litt(N_C_TYPES)
  real    :: leaf_fall, root_mort, stem_mort ! biomass losses from various tissues
  real    :: dheat ! heat residual due to cohort merging

  ! rates of various tissues decay after leaf drop, per day
  real, parameter :: leaf_fall_rate = 0.075
  real, parameter :: root_mort_rate = 0.0
  real, parameter :: stem_mort_rate = 0.075 ! isa 20170705

  ! isa 20170709 - note that in biomass_allocation_ppa deltaDBH is an output to
  ! subroutine vegn_growth, where it is saved ...
  real :: wilt       ! vertically-average soil moisture content at wilting point
  real :: theta_crit ! critical soil moisture for drought-deciduous phenology
  logical :: drought ! drought indicator


  associate(vegn=>tile%vegn, soil=>tile%soil)
  vegn%litter = 0; leaf_litt(:) = 0.0
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

     ! leaf falling at the end of a growing season
     if(cc%status == LEAF_OFF) then
         leaf_fall = leaf_fall_rate * max(cc%bl,0.0)
         root_mort = root_mort_rate * max(cc%br,0.0)
         stem_mort = 0.0 ! isa 20170705
         if (sp%lifeform == FORM_GRASS) then
            ! grasses also lose their stems
            stem_mort = min(stem_mort_rate * cc%bsw, &
                   cc%bsw - sp%rho_wood * sp%alphaBM * ((sp%gammaHT/(sp%alphaHT/sp%seedling_height - 1.0))**(1.0/sp%thetaHT))**2 * sp%seedling_height)
            stem_mort = max(stem_mort,0.0)
            ! ToDo - it is necessary to implement anonline adjustment of dbh, height and crown area as the plant shrinks
            !        otherwise it can happen that, in a year with a short winter, there is a disadjustment between plant
            !        biomass and its dimensions
            ! just set initial height
            cc%height = sp%seedling_height
            cc%dbh = (sp%gammaHT/(sp%alphaHT/sp%seedling_height - 1.0))**(1.0/sp%thetaHT)
            cc%crownarea = sp%alphaCA * cc%dbh**sp%thetaCA
            ! slm 20170804: reconcile bl_max and crownarea, and drop excess leaf biomass immediately
            cc%bl_max = sp%LMA * sp%laimax * cc%crownarea * (1.0-sp%internal_gap_frac)
            leaf_fall = max(leaf_fall, cc%bl-cc%bl_max)
         endif
         ! drop all leaves if LAI is below certain small minimum
         if (cc%lai<min_lai) leaf_fall = cc%bl

         cc%nsc = cc%nsc + l_fract * (leaf_fall+ root_mort + stem_mort) ! isa 20170705
         cc%bl  = cc%bl - leaf_fall
         cc%br  = cc%br - root_mort
         cc%bsw = cc%bsw - stem_mort ! isa 20170705
         cc%lai = leaf_area_from_biomass(cc%bl,cc%species,cc%layer,cc%firstlayer)/ &
                                        (cc%crownarea *(1.0-sp%internal_gap_frac))
         if(cc%bl == 0.0) cc%leaf_age = 0.0
         cc%bliving = cc%blv + cc%br + cc%bl + cc%bsw

         leaf_litter = (1.0-l_fract) * (leaf_fall+root_mort+stem_mort) * cc%nindivs ! isa20170705, stem and roots become leaf litter
         ! call check_var_range(leaf_litter,0.0,HUGE(1.0),'vegn_phenology_ppa','leaf_litter',FATAL)
         leaf_litt(:) = leaf_litt(:)+[fsc_liv,1-fsc_liv,0.0]*leaf_litter

         vegn%litter    = vegn%litter     + leaf_litter
         soil%fsc_in(1) = soil%fsc_in(1)  + leaf_litter
         vegn%veg_out   = vegn%veg_out    + leaf_litter
     endif
     end associate ! cc, sp
  enddo
  ! add litter accumulated over the cohorts
  call add_soil_carbon(soil, leaf_litter=leaf_litt)
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
  integer :: i

  do i = 1, vegn%n_cohorts
     call update_species(vegn%cohorts(i), vegn%t_ann, vegn%t_cold, &
          vegn%p_ann*seconds_per_year, vegn%ncm, vegn%landuse)
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
  real :: deltafast, deltaslow;
  real, dimension(N_C_TYPES,N_LITTER_POOLS) :: delta_C

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
  case (SOILC_CORPSE)
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

     ! update litter rate so that litter pool is never depleted below zero
     vegn%litter_rate_C = MAX(0.0, MIN(vegn%litter_rate_C, vegn%litter_buff_C/dt_fast_yr))
     ! move carbon from litter pool to litter
     delta_C = vegn%litter_rate_C*dt_fast_yr
     call add_litter(soil%leafLitter,       delta_C(:,LEAF))
     call add_litter(soil%coarsewoodLitter, delta_C(:,CWOOD))
     vegn%litter_buff_C = vegn%litter_buff_C - delta_C

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
     !soil%prog(1)%slow_soil_C = soil%prog(1)%slow_soil_C + delta;
     vegn%ssc_pool_bg    = vegn%ssc_pool_bg    - deltaslow;
     ! TODO: figure out how to distribute the additional soil litter vertically.
     ! Strictly speaking we need to keep pools and rates distributed vertically,
     ! For now, add belowground litter using vertical profile of the first cohort.
     ! perhaps we can use average cohort root profile...
     call add_root_litter(soil, vegn%cohorts(1), (/deltafast,deltaslow,0.0/))
  case default
     call error_mesg('update_soil_pools','The value of soil_carbon_option is invalid. This should never happen. Contact developer.',FATAL)
  end select
end subroutine update_soil_pools


! ============================================================================
logical function cohort_can_reproduce(cc)
  type(vegn_cohort_type), intent(in) :: cc

  cohort_can_reproduce = (cc%layer==1.or.spdata(cc%species)%reproduces_in_understory) &
                         .and. cc%age > spdata(cc%species)%maturalage
end function cohort_can_reproduce

! ============================================================================
subroutine vegn_reproduction_ppa(do_seed_transport)
  logical, intent(in) :: do_seed_transport

  type(land_tile_enum_type) :: ce
  type(land_tile_type), pointer :: tile
  real :: ug_dispersed   (lnd%ls:lnd%le,0:nspecies-1)  ! dispersed seeds, kgC per m2 of land
  real :: ug_transported (lnd%ls:lnd%le,0:nspecies-1)  ! dispersed seeds, kgC per m2 of land
  real, allocatable :: bseed(:,:)
  integer :: n ! total number of vegetated tiles in the domain
  integer :: k ! vegetated tile index
  integer :: l ! grid cell index in unstructured grid
  integer :: s ! species iterator
  integer :: i ! cohort iterator
  integer :: k1

  real :: btot0, btot1 ! total carbon, for conservation check only

  ! + conservation check part 1
  if (do_check_conservation) then
     btot0 = 0.0
     ce = first_elmt(land_tile_map,lnd%ls)
     do while (loop_over_tiles(ce,tile,l))
        btot0 = btot0 + land_tile_carbon(tile) * lnd%ug_area(l) * tile%frac
     end do
     call mpp_sum(btot0)
  end if
  ! - conservation check part 1

  ! count total number of vegetation tiles in our domain
  n  = 0
  ce = first_elmt(land_tile_map,lnd%ls)
  do while (loop_over_tiles(ce,tile,l))
     if(associated(tile%vegn)) n = n+1
  end do

  ! calculate amount of seeds (kgC per tile area, by species) for each of the vegetation tiles
  allocate(bseed(n,0:nspecies-1))
  bseed = 0.0
  ce = first_elmt(land_tile_map, lnd%ls); k = 1
  do while (loop_over_tiles(ce,tile,l))
     if(.not.associated(tile%vegn)) cycle
     do i = 1,tile%vegn%n_cohorts
        associate(cc=>tile%vegn%cohorts(i))
        if (cohort_can_reproduce(cc)) then
           bseed(k,cc%species) = bseed(k,cc%species) + cc%bseed * cc%nindivs
           cc%bseed = 0.0
        endif
        end associate
     enddo
     k = k+1
  enddo

  ! calculate amount of dispersed seeds, by species, kgC per m2 of land
  ug_dispersed = 0.0
  ug_transported = 0.0
  if (do_seed_transport) then
     ce = first_elmt(land_tile_map, lnd%ls); k = 1
     do while (loop_over_tiles(ce,tile,l)) ! l is the grid cell index in unstructured grid
        if (.not.associated(tile%vegn)) cycle
        do s = 0,nspecies-1
           ug_dispersed(l,s)   = ug_dispersed(l,s)   + bseed(k,s)*spdata(s)%frac_seed_dispersed*tile%frac
           ug_transported(l,s) = ug_transported(l,s) + bseed(k,s)*spdata(s)%frac_seed_transported*tile%frac
           ! correct the amount of seeds remaining in each tile
           bseed(k,s) = bseed(k,s)*(1-spdata(s)%frac_seed_dispersed-spdata(s)%frac_seed_transported)
        enddo
        k = k+1
     enddo
     ! diffuse the seeds
     do s = 0,nspecies-1
        call transport_seeds(ug_transported(:,s))
     enddo
  endif

  ! distribute seeds (those that stay in place + dispersed) among tiles
  ce = first_elmt(land_tile_map, lnd%ls); k = 1
  do while (loop_over_tiles(ce,tile,l,k1))
     call set_current_point(l,k1)
     if (.not.associated(tile%vegn)) cycle
     if (is_watch_point()) then
        __DEBUG3__(l,k,k1)
        __DEBUG1__(ug_dispersed(l,:))
        __DEBUG1__(ug_transported(l,:))
        __DEBUG1__(bseed(k,:))
     endif
     if (tile%vegn%landuse==LU_CROP .and. .not.allow_weeds_on_crops) then
        call add_seedlings_ppa(tile%vegn,tile%soil,(ug_dispersed(l,:)+ug_transported(l,:))*ug_area_factor(l)+bseed(k,:), &
                               germination_factor = 0.0) ! no seeds germinate
        ! This also means that the crops are not allowed to reproduce by themselves.
     else
        call add_seedlings_ppa(tile%vegn,tile%soil,(ug_dispersed(l,:)+ug_transported(l,:))*ug_area_factor(l)+bseed(k,:), &
                               germination_factor = 1.0) ! do not modify natural seed germination
     endif
     k = k+1
  enddo

  ! + conservation check part 2
  if (do_check_conservation) then
     btot1 = 0.0
     ce = first_elmt(land_tile_map,lnd%ls)
     do while (loop_over_tiles(ce,tile,l))
        btot1 = btot1 + land_tile_carbon(tile) * lnd%ug_area(l) * tile%frac
     end do
     call mpp_sum(btot1)
     if (mpp_pe()==mpp_root_pe()) then
        call check_conservation ('vegn_reproduction_ppa','total carbon', &
             btot0/atot, btot1/atot, carbon_cons_tol, severity=FATAL)
     endif
  end if
  ! TODO: add nitrogen conservation check
  ! - conservation check part 2
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
