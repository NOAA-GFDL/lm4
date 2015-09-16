! ============================================================================
! updates carbon pools and rates on the fast time scale
! ============================================================================
module vegn_dynamics_mod

#include "../shared/debug.inc"

use fms_mod, only: write_version_number, error_mesg, FATAL
use time_manager_mod, only: time_type

use land_constants_mod, only : seconds_per_year, mol_C
use land_tile_diag_mod, only : &
     register_tiled_diag_field, send_tile_data, diag_buff_type
use vegn_data_mod, only : spdata, &
     CMPT_VLEAF, CMPT_SAPWOOD, CMPT_ROOT, CMPT_WOOD, CMPT_LEAF, LEAF_ON, LEAF_OFF, &
     fsc_liv, fsc_wood, fsc_froot, soil_carbon_depth_scale, C2B, agf_bs, &
     l_fract, &
     leaf_fast_c2n,leaf_slow_c2n,froot_fast_c2n,froot_slow_c2n,wood_fast_c2n,wood_slow_c2n, root_exudate_N_frac,& !x2z - ens: lets get rid of c2n?
     ! BNS: C2N ratios should be temporary fix, which we can get rid of once N is integrated into vegetation code
     root_exudate_frac_max, dynamic_root_exudation, c2n_mycorrhizae, mycorrhizal_turnover_time, myc_C_efficiency,&
     N_fixer_turnover_time, N_fixer_C_efficiency, N_fixation_rate, c2n_N_fixer

use vegn_tile_mod, only: vegn_tile_type
use soil_tile_mod, only: soil_tile_type
use vegn_cohort_mod, only : vegn_cohort_type, &
     update_biomass_pools, update_bio_living_fraction, update_species
use soil_carbon_mod, only: soil_carbon_option, &
    SOILC_CENTURY, SOILC_CENTURY_BY_LAYER, SOILC_CORPSE, SOILC_CORPSE_N, &
    add_litter, poolTotals, debug_pool,soil_NO3_deposition,soil_NH4_deposition

use soil_mod, only: add_root_litter, add_root_exudates, Dsdt, myc_scavenger_N_uptake, hypothetical_myc_scavenger_N_uptake

use land_debug_mod, only : is_watch_point

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
character(len=*), private, parameter :: &
   version = '$Id$', &
   tagname = '$Name$' ,&
   module_name = 'vegn'

real, parameter :: GROWTH_RESP=0.333  ! fraction of npp lost as growth respiration

! ==== module data ===========================================================
real    :: dt_fast_yr ! fast (physical) time step, yr (year is defined as 365 days)

! diagnostic field IDs
integer :: id_npp, id_nep, id_gpp, id_resp, id_resl, id_resr, id_resg, &
    id_soilt, id_theta, id_litter, &
    id_mycorrhizal_allocation, id_mycorrhizal_immobilization, id_N_fixer_allocation, id_N_fix_marginal_gain, id_myc_marginal_gain


contains

! ============================================================================
subroutine vegn_dynamics_init(id_lon, id_lat, time, delta_time)
  integer        , intent(in) :: id_lon ! ID of land longitude (X) axis
  integer        , intent(in) :: id_lat ! ID of land latitude (Y) axis
  type(time_type), intent(in) :: time       ! initial time for diagnostic fields
  real           , intent(in) :: delta_time ! fast time step, s

  call write_version_number(version, tagname)

  ! set up global variables
  dt_fast_yr = delta_time/seconds_per_year

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
  id_mycorrhizal_allocation = register_tiled_diag_field ( module_name, 'mycorrhizal_allocation',  &
       (/id_lon,id_lat/), time, 'C allocation to mycorrhizae', 'kg C/(m2 year)', &
       missing_value=-100.0 )
   id_mycorrhizal_immobilization = register_tiled_diag_field ( module_name, 'mycorrhizal_immobilization',  &
        (/id_lon,id_lat/), time, 'N immobilization by mycorrhizae', 'kg N/(m2 year)', &
        missing_value=-100.0 )
    id_N_fixer_allocation = register_tiled_diag_field ( module_name, 'N_fixer_allocation',  &
         (/id_lon,id_lat/), time, 'C allocation to N fixers', 'kg C/(m2 year)', &
         missing_value=-100.0 )
     id_N_fix_marginal_gain = register_tiled_diag_field ( module_name, 'N_fix_marginal_gain',  &
          (/id_lon,id_lat/), time, 'Extra N fixation per unit C allocation', 'kg N/(m2 year)', &
          missing_value=-100.0 )
    id_myc_marginal_gain = register_tiled_diag_field ( module_name, 'myc_marginal_gain',  &
         (/id_lon,id_lat/), time, 'Extra N acquisition per unit C allocation to mycorrhizae', 'kg N/(m2 year)', &
         missing_value=-100.0 )
end subroutine vegn_dynamics_init


! ============================================================================
subroutine vegn_carbon_int(vegn, soil, soilt, theta, diag)
  type(vegn_tile_type), intent(inout) :: vegn
  type(soil_tile_type), intent(inout) :: soil
  real, intent(in) :: soilt ! average temperature of soil for soil carbon decomposition, deg K
  real, intent(in) :: theta ! average soil wetness, unitless
  type(diag_buff_type), intent(inout) :: diag

  ! TODO: possibly move soil-related calculations from calling procedure here,
  !       now that we have soil passed as an argument

  type(vegn_cohort_type), pointer :: cc
  real :: resp, resl, resr, resg ! respiration terms accumulated for all cohorts
  real :: cgain, closs ! carbon gain and loss accumulated for entire tile
  real :: md_alive, md_leaf, md_wood, md_froot ! component of maintenance demand
  real :: gpp ! gross primary productivity per tile
  real :: root_exudate_C, total_root_exudate_C, root_exudate_frac, myc_exudate_frac, mycorrhizal_N_immob
  real :: myc_N_uptake,total_scavenger_myc_C_allocated,scavenger_myc_C_allocated,total_myc_immob
  real :: N_fixation, total_N_fixation, total_N_fixer_C_allocated, N_fixer_C_allocated, N_fixer_exudate_frac
  real :: myc_marginal_gain, N_fix_marginal_gain,rhiz_exud_frac
  real :: N_fixation_2, myc_N_uptake_2
  integer :: sp ! shorthand for current cohort specie
  integer :: i

  if(is_watch_point()) then
     write(*,*)'#### vegn_carbon_int ####'
     __DEBUG2__(soilt,theta)
  endif

total_scavenger_myc_C_allocated = 0.0
total_myc_immob = 0.0
total_N_fixation = 0.0
total_N_fixer_C_allocated = 0.0

  !  update plant carbon
  vegn%npp = 0
  resp = 0 ; resl = 0 ; resr = 0 ; resg = 0 ; gpp = 0
  cgain = 0 ; closs = 0
  total_root_exudate_C = 0
  soil%myc_min_N_uptake = 0
  soil%symbiotic_N_fixation = 0
  do i = 1, vegn%n_cohorts
     cc => vegn%cohorts(i)
     sp = cc%species

     call eddy_npp(cc,soilt);
     ! npp2 is for diagnostics and comparison
     cc%npp2 = cc%miami_npp;  ! treat miami npp as above+below npp

     if(dynamic_root_exudation .AND. soil_carbon_option==SOILC_CORPSE_N) then
       ! Initial allocation scheme: root exudation/mycorrhizal allocation depends on ratio of leaf biomass to max (as determined by N uptake)
       ! Root exudation fraction of NPP limited by some maximum value. Note that bl+br shouldn't exceed max_live_biomass
       root_exudate_frac = root_exudate_frac_max*(cc%bl+cc%br)/cc%max_live_biomass
     else
       root_exudate_frac = spdata(sp)%root_exudate_frac
     endif
     root_exudate_C = max(cc%npp,0.0)*root_exudate_frac
     cc%carbon_gain = cc%carbon_gain + (cc%npp-root_exudate_C)*dt_fast_yr


     ! check if leaves/roots are present and need to be accounted in maintenance
     if(cc%status == LEAF_ON) then
        md_alive = (cc%Pl * spdata(sp)%alpha(CMPT_LEAF) + &
                    cc%Pr * spdata(sp)%alpha(CMPT_ROOT))* &
              cc%bliving*dt_fast_yr
        md_leaf = cc%Pl * spdata(sp)%alpha(CMPT_LEAF)*cc%bliving*dt_fast_yr
        md_froot= cc%Pr * spdata(sp)%alpha(CMPT_ROOT)*cc%bliving*dt_fast_yr
        ! NOTE that mathematically. md_alive = md_leaf + md_froot. Unfortunately,
        ! order of operation matters for the bit-wise reproducibility, so all
        ! three need to be calculated separately
     else
        md_alive = 0
        md_leaf  = 0
	md_froot = 0
     endif

     ! compute branch and coarse wood losses for tree types
     md_wood =0;
     if (sp > 1) then
        md_wood = 0.6 *cc%bwood * spdata(sp)%alpha(CMPT_WOOD)*dt_fast_yr;
     endif

     ! this silliness with the mathematically equivalent formulas is
     ! solely to bit-reproduce old results in both CENTURY and CORPSE
     ! modes: the results are extremely sensitive to the order of operations,
     ! and diverge with time, even in stand-alone land model runs.
     select case (soil_carbon_option)
     case (SOILC_CENTURY, SOILC_CENTURY_BY_LAYER)
        cc%md = md_alive + cc%Psw_alphasw * cc%bliving * dt_fast_yr
    case (SOILC_CORPSE, SOILC_CORPSE_N)
        cc%md = md_leaf + md_froot + cc%Psw_alphasw * cc%bliving * dt_fast_yr
     end select

     cc%bwood_gain = cc%bwood_gain + cc%Psw_alphasw * cc%bliving * dt_fast_yr;
     cc%bwood_gain = cc%bwood_gain - md_wood;
     if (cc%bwood_gain < 0.0) cc%bwood_gain=0.0; ! potential non-conservation ?
     cc%carbon_gain = cc%carbon_gain - cc%md;
     cc%carbon_loss = cc%carbon_loss + cc%md; ! used in diagnostics only

     ! add md from leaf and root pools to fast soil carbon
     select case (soil_carbon_option)
     case (SOILC_CENTURY, SOILC_CENTURY_BY_LAYER)
        soil%fast_soil_C(1) = soil%fast_soil_C(1) +    fsc_liv *md_alive +    fsc_wood *md_wood;
        soil%slow_soil_C(1) = soil%slow_soil_C(1) + (1-fsc_liv)*md_alive + (1-fsc_wood)*md_wood;
     ! for budget tracking
!/*     cp->fsc_in(1)+= data->fsc_liv*md_alive+data->fsc_wood*md_wood; */
!/*     cp->ssc_in(1)+= (1.- data->fsc_liv)*md_alive+(1-data->fsc_wood)*md_wood; */
        soil%fsc_in(1)  = soil%fsc_in(1) + 1*md_alive+0*md_wood;
        soil%ssc_in(1)  = soil%ssc_in(1) + (1.- 1)*md_alive+(1-0)*md_wood;
    case (SOILC_CORPSE_N)
        if (is_watch_point()) then
           call debug_pool(soil%leafLitter,      'leafLitter (before)'      )
           call debug_pool(soil%coarseWoodLitter,'coarseWoodLitter (before)')
        endif
        call add_litter(soil%leafLitter,(/fsc_liv *md_leaf, (1-fsc_liv)*md_leaf ,0.0/),(/fsc_liv *md_leaf/leaf_fast_c2n, (1-fsc_liv)*md_leaf/leaf_slow_c2n ,0.0/))
        call add_litter(soil%coarseWoodLitter,(/fsc_wood *md_wood*agf_bs, (1-fsc_wood)*md_wood*agf_bs, 0.0/),(/fsc_wood *md_wood*agf_bs/wood_fast_c2n, (1-fsc_wood)*md_wood*agf_bs/wood_slow_c2n, 0.0/))
        soil%leaflitter_fsc_in=soil%leaflitter_fsc_in+fsc_liv *md_leaf
        soil%leaflitter_ssc_in=soil%leaflitter_ssc_in+(1-fsc_liv)*md_leaf
        soil%coarsewoodlitter_fsc_in=soil%coarsewoodlitter_fsc_in +    fsc_wood *md_wood*agf_bs
        soil%coarsewoodlitter_ssc_in=soil%coarsewoodlitter_ssc_in + (1-fsc_wood)*md_wood*agf_bs
        soil%leaflitter_fsn_in=soil%leaflitter_fsn_in+fsc_liv *md_leaf/leaf_fast_c2n
        soil%leaflitter_ssn_in=soil%leaflitter_ssn_in+(1-fsc_liv)*md_leaf/leaf_slow_c2n
        soil%coarsewoodlitter_fsn_in=soil%coarsewoodlitter_fsn_in +    fsc_wood *md_wood*agf_bs/wood_fast_c2n
        soil%coarsewoodlitter_ssn_in=soil%coarsewoodlitter_ssn_in + (1-fsc_wood)*md_wood*agf_bs/wood_slow_c2n
	!ssc_in and fsc_in updated in add_root_litter
        call add_root_litter(soil,vegn,(/fsc_froot*md_froot + fsc_wood*md_wood*(1-agf_bs),&
					(1-fsc_froot)*md_froot + (1-fsc_wood)*md_wood*(1-agf_bs),0.0/), &
                    (/fsc_froot*md_froot/froot_fast_c2n + fsc_wood*md_wood*(1-agf_bs)/wood_fast_c2n,&
            					(1-fsc_froot)*md_froot/froot_slow_c2n + (1-fsc_wood)*md_wood*(1-agf_bs)/wood_slow_c2n,0.0/))
	if (is_watch_point()) then
           call debug_pool(soil%leafLitter,      'leafLitter (after)'      )
           call debug_pool(soil%coarseWoodLitter,'coarseWoodLitter (after)')
        endif
    case (SOILC_CORPSE)
        if (is_watch_point()) then
            call debug_pool(soil%leafLitter,      'leafLitter (before)'      )
            call debug_pool(soil%coarseWoodLitter,'coarseWoodLitter (before)')
        endif
        call add_litter(soil%leafLitter,(/fsc_liv *md_leaf, (1-fsc_liv)*md_leaf ,0.0/),(/0.0,0.0 ,0.0/))
        call add_litter(soil%coarseWoodLitter,(/fsc_wood *md_wood*agf_bs, (1-fsc_wood)*md_wood*agf_bs, 0.0/),(/0.0,0.0, 0.0/))
        soil%leaflitter_fsc_in=soil%leaflitter_fsc_in+fsc_liv *md_leaf
        soil%leaflitter_ssc_in=soil%leaflitter_ssc_in+(1-fsc_liv)*md_leaf
        soil%coarsewoodlitter_fsc_in=soil%coarsewoodlitter_fsc_in +    fsc_wood *md_wood*agf_bs
        soil%coarsewoodlitter_ssc_in=soil%coarsewoodlitter_ssc_in + (1-fsc_wood)*md_wood*agf_bs

        !ssc_in and fsc_in updated in add_root_litter
        call add_root_litter(soil,vegn,(/fsc_froot*md_froot + fsc_wood*md_wood*(1-agf_bs),&
        (1-fsc_froot)*md_froot + (1-fsc_wood)*md_wood*(1-agf_bs),0.0/), &
        (/0.0,0.0,0.0/))
        if (is_watch_point()) then
            call debug_pool(soil%leafLitter,      'leafLitter (after)'      )
            call debug_pool(soil%coarseWoodLitter,'coarseWoodLitter (after)')
        endif

     case default
        call error_mesg('vegn_carbon_int','The value of soil_carbon_option is invalid. This should never happen. Contact developer.',FATAL)
     end select


     vegn%veg_in  = vegn%veg_in  + cc%npp*dt_fast_yr;
     vegn%veg_out = vegn%veg_out + md_alive+md_wood;

     if(is_watch_point()) then
        __DEBUG3__(md_wood,md_leaf,md_froot)
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


         ! Calculate return on investment for each N acquisition strategy
         ! Determined by calculating marginal change in N uptake if C allocation to each strategy increased by 10%
         ! Mycorrhizal:
         call hypothetical_myc_scavenger_N_uptake(soil,vegn,cc%myc_scavenger_biomass_C,myc_N_uptake,dt_fast_yr)
         call hypothetical_myc_scavenger_N_uptake(soil,vegn,cc%myc_scavenger_biomass_C+root_exudate_C*0.1*dt_fast_yr*myc_C_efficiency,myc_N_uptake_2,dt_fast_yr)
         myc_marginal_gain = (myc_N_uptake_2-myc_N_uptake)/0.1

         ! N fixer
         N_fixation = cc%N_fixer_biomass_C*N_fixation_rate*dt_fast_yr
         N_fixation_2 = (cc%N_fixer_biomass_C+root_exudate_C*0.1*dt_fast_yr*N_fixer_C_efficiency)*N_fixation_rate*dt_fast_yr
         N_fix_marginal_gain = (N_fixation_2-N_fixation)/0.1

         ! Just exudation
         ! Need to think about a strategy for calculating this.  Leaving at fixed value for now

         ! Calculate relative fractions
         rhiz_exud_frac = 0.3  ! Using fixed value for now

         if (myc_marginal_gain+N_fix_marginal_gain>0) then
           myc_exudate_frac = (1.0-rhiz_exud_frac)*myc_marginal_gain/(myc_marginal_gain+N_fix_marginal_gain)
           N_fixer_exudate_frac = (1.0-rhiz_exud_frac)*N_fix_marginal_gain/(myc_marginal_gain+N_fix_marginal_gain)
         else
           ! Divide evenly if there is no marginal gain.  But this probably only happens if root_exudate_C is zero
           myc_exudate_frac = (1.0-rhiz_exud_frac)*0.5
           N_fixer_exudate_frac = (1.0-rhiz_exud_frac)*0.5
         endif

         ! Do actual mycorrhizal mineral N uptake now
         call myc_scavenger_N_uptake(soil,vegn,cc%myc_scavenger_biomass_C,myc_N_uptake,dt_fast_yr)

        !  myc_exudate_frac=0.5  ! This will eventually be dynamic based on N acquisition cost function
           ! Watch out that mycorrhizae aren't too N limited (for instance if starting from zero biomass)
         scavenger_myc_C_allocated=root_exudate_C*myc_exudate_frac*dt_fast_yr

         mycorrhizal_N_immob = max(0.0,scavenger_myc_C_allocated/c2n_mycorrhizae-scavenger_myc_C_allocated*root_exudate_N_frac)
         mycorrhizal_N_immob = min(mycorrhizal_N_immob,myc_N_uptake)

        !  N_fixer_exudate_frac = 0.2 ! This will eventually be dynamic based on N acquisition cost function
         N_fixer_C_allocated = root_exudate_C*N_fixer_exudate_frac*dt_fast_yr



     else
         scavenger_myc_C_allocated=0.0
         myc_N_uptake=0.0
         mycorrhizal_N_immob = 0.0
         N_fixer_C_allocated = 0.0
     endif


     soil%myc_min_N_uptake=soil%myc_min_N_uptake+myc_N_uptake-mycorrhizal_N_immob

     N_fixation = cc%N_fixer_biomass_C*N_fixation_rate*dt_fast_yr
     soil%symbiotic_N_fixation=soil%symbiotic_N_fixation+N_fixation

     if(cc%myc_scavenger_biomass_C<0) call error_mesg('vegn_carbon_int','Mycorrhizal biomass < 0',FATAL)
     if(cc%N_fixer_biomass_C<0) call error_mesg('vegn_carbon_int','N fixer biomass < 0',FATAL)

     ! First add mycorrhizal and N fixer turnover to soil C pools
     call add_root_litter(soil,vegn,(/0.0,0.0,cc%myc_scavenger_biomass_C/mycorrhizal_turnover_time/)*dt_fast_yr,(/0.0,0.0,cc%myc_scavenger_biomass_N/mycorrhizal_turnover_time/)*dt_fast_yr)
     call add_root_litter(soil,vegn,(/0.0,0.0,cc%N_fixer_biomass_C/N_fixer_turnover_time/)*dt_fast_yr,(/0.0,0.0,cc%N_fixer_biomass_N/N_fixer_turnover_time/)*dt_fast_yr)
     ! Then update biomass of scavenger mycorrhizae and N fixers
     cc%myc_scavenger_biomass_C = cc%myc_scavenger_biomass_C + scavenger_myc_C_allocated*myc_C_efficiency - cc%myc_scavenger_biomass_C/mycorrhizal_turnover_time*dt_fast_yr
     cc%myc_scavenger_biomass_N = cc%myc_scavenger_biomass_N + scavenger_myc_C_allocated*myc_C_efficiency/c2n_mycorrhizae - cc%myc_scavenger_biomass_N/mycorrhizal_turnover_time*dt_fast_yr
     cc%N_fixer_biomass_C = cc%N_fixer_biomass_C + N_fixer_C_allocated*N_fixer_C_efficiency - cc%N_fixer_biomass_C/N_fixer_turnover_time*dt_fast_yr
     cc%N_fixer_biomass_N = cc%N_fixer_biomass_N + N_fixer_C_allocated*N_fixer_C_efficiency/c2n_N_fixer - cc%N_fixer_biomass_N/N_fixer_turnover_time*dt_fast_yr

     total_root_exudate_C = total_root_exudate_C + root_exudate_C*dt_fast_yr - scavenger_myc_C_allocated - N_fixer_C_allocated
     total_scavenger_myc_C_allocated = total_scavenger_myc_C_allocated + scavenger_myc_C_allocated
     total_myc_immob = total_myc_immob + mycorrhizal_N_immob

     total_N_fixation = total_N_fixation + N_fixation
     total_N_fixer_C_allocated = total_N_fixer_C_allocated + N_fixer_C_allocated

   enddo

   ! fsc_in and ssc_in updated in add_root_exudates
   call add_root_exudates(soil,vegn,total_root_exudate_C,total_root_exudate_C*root_exudate_N_frac)


  ! update soil carbon
  call Dsdt(vegn, soil, diag, soilt, theta)

  ! Add respiration/C waste from mycorrhizae and N fixers
  vegn%rh = vegn%rh + (total_scavenger_myc_C_allocated*(1-myc_C_efficiency) + total_N_fixer_C_allocated*(1-N_fixer_C_efficiency))/dt_fast_yr

  ! NEP is equal to NPP minus soil respiration
  vegn%nep = vegn%npp - vegn%rh
  if(is_watch_point()) then
     __DEBUG3__(vegn%npp,vegn%rh,vegn%nep)
  endif

  call update_soil_pools(vegn, soil)
  vegn%age = vegn%age + dt_fast_yr;


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
  call send_tile_data(id_mycorrhizal_allocation,total_scavenger_myc_C_allocated/dt_fast_yr,diag)
  call send_tile_data(id_mycorrhizal_immobilization,total_myc_immob/dt_fast_yr,diag)
  call send_tile_data(id_N_fixer_allocation,total_N_fixer_C_allocated/dt_fast_yr,diag)
  call send_tile_data(id_myc_marginal_gain,myc_marginal_gain,diag)
  call send_tile_data(id_N_fix_marginal_gain,N_fix_marginal_gain,diag)

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

     ! reset carbon acculmulation terms
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
  integer :: i

  wilt = soil%w_wilt(1)/soil%pars%vwc_sat
  vegn%litter = 0

  do i = 1,vegn%n_cohorts
     cc => vegn%cohorts(i)

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
              call add_litter(soil%leafLitter,(/fsc_liv*leaf_litter,(1-fsc_liv)*leaf_litter,0.0/),&
                                                (/fsc_liv*leaf_litter/leaf_fast_c2n,(1-fsc_liv)*leaf_litter/leaf_slow_c2n,0.0/))
              soil%leaflitter_fsc_in=soil%leaflitter_fsc_in+fsc_liv*leaf_litter
              soil%leaflitter_ssc_in=soil%leaflitter_ssc_in+(1-fsc_liv)*leaf_litter
              soil%leaflitter_fsn_in=soil%leaflitter_fsn_in+fsc_liv*leaf_litter/leaf_slow_c2n
              soil%leaflitter_ssn_in=soil%leaflitter_ssn_in+(1-fsc_liv)*leaf_litter/leaf_slow_c2n
              !ssc_in and fsc_in updated in add_root_litter
              call add_root_litter(soil, vegn, (/fsc_froot*root_litter,(1-fsc_froot)*root_litter,0.0/),&
                                (/fsc_froot*root_litter/froot_fast_c2n,(1-fsc_froot)*root_litter/froot_slow_c2n,0.0/))
            case(SOILC_CORPSE)
                call add_litter(soil%leafLitter,(/fsc_liv*leaf_litter,(1-fsc_liv)*leaf_litter,0.0/),&
                            (/0.0,0.0,0.0/))
                soil%leaflitter_fsc_in=soil%leaflitter_fsc_in+fsc_liv*leaf_litter
                soil%leaflitter_ssc_in=soil%leaflitter_ssc_in+(1-fsc_liv)*leaf_litter
                !ssc_in and fsc_in updated in add_root_litter
                call add_root_litter(soil, vegn, (/fsc_froot*root_litter,(1-fsc_froot)*root_litter,0.0/),&
                            (/0.0,0.0,0.0/))
           case default
              call error_mesg('vegn_phenology','The value of soil_carbon_option is invalid. This should never happen. Contact developer.',FATAL)
           end select
           vegn%veg_out = vegn%veg_out + leaf_litter+root_litter;

           cc%blv = cc%blv + l_fract*(cc%bl+cc%br);
           cc%bl  = 0.0;
           cc%br  = 0.0;
           cc%lai = 0.0;

           ! update state
           cc%bliving = cc%blv + cc%br + cc%bl + cc%bsw;
           cc%b = cc%bliving + cc%bwood ;
           call update_bio_living_fraction(cc);
        endif
     endif
  enddo
end subroutine vegn_phenology


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
  real :: deltafast, deltaslow, deltafast_N, deltaslow_N;

  !!!!!!!!!!!!!!!x2z /ADD deposition CHECK
  real ::   TOT_annual_NH4_dep=0.0001
  real ::   TOT_annual_NO3_dep=0.0001
  real ::   TOT_annual_NH4_fert=0.0
  real ::  TOT_annual_NO3_fert=0.0
!  real ::  input_time_fert=130.0

if(TOT_annual_NH4_dep.gt.dfloat(0))then
		call soil_NH4_deposition(TOT_annual_NH4_dep*dt_fast_yr,soil%leafLitter)
	endif

	if(TOT_annual_NO3_dep.gt.dfloat(0))then
		call soil_NO3_deposition(TOT_annual_NO3_dep*dt_fast_yr,soil%leafLitter)
	endif

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
     vegn%leaflitter_buffer_rate_fast_N = MAX(0.0, MIN(vegn%leaflitter_buffer_rate_fast_N, vegn%leaflitter_buffer_fast_N/dt_fast_yr))

     deltafast = vegn%leaflitter_buffer_rate_fast*dt_fast_yr
     deltaslow = vegn%leaflitter_buffer_rate_slow*dt_fast_yr

     if(soil_carbon_option == SOILC_CORPSE_N) then
         vegn%leaflitter_buffer_rate_slow = MAX(0.0, MIN(vegn%leaflitter_buffer_rate_slow, vegn%leaflitter_buffer_slow/dt_fast_yr))
         vegn%leaflitter_buffer_rate_slow_N = MAX(0.0, MIN(vegn%leaflitter_buffer_rate_slow_N, vegn%leaflitter_buffer_slow_N/dt_fast_yr))
         deltafast_N = vegn%leaflitter_buffer_rate_fast_N*dt_fast_yr
         deltaslow_N = vegn%leaflitter_buffer_rate_slow_N*dt_fast_yr
     else
         vegn%leaflitter_buffer_rate_slow =0.0
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
        vegn%coarsewoodlitter_buffer_rate_slow = 0.0
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
