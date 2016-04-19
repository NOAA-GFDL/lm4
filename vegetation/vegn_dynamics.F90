! ============================================================================
! updates carbon pools and rates on the fast time scale
! ============================================================================
module vegn_dynamics_mod

#include "../shared/debug.inc"

use fms_mod, only: error_mesg, FATAL, WARNING
use time_manager_mod, only: time_type

use constants_mod, only : PI,tfreeze
use land_constants_mod, only : seconds_per_year, mol_C
use land_data_mod, only : log_version
use land_debug_mod, only : is_watch_point, check_var_range
use land_tile_diag_mod, only : OP_SUM, OP_MEAN, &
     register_tiled_diag_field, send_tile_data, diag_buff_type, &
     register_cohort_diag_field, send_cohort_data, set_default_diag_filter
use vegn_data_mod, only : spdata, &
     CMPT_NSC, CMPT_SAPWOOD, CMPT_LEAF, CMPT_ROOT, CMPT_VLEAF, CMPT_WOOD, &
     PHEN_DECIDUOUS, LEAF_ON, LEAF_OFF, FORM_WOODY, FORM_GRASS, &
     fsc_liv, fsc_wood, fsc_froot, agf_bs, &
     l_fract, mcv_min, mcv_lai, do_ppa, tau_seed, &
     understory_lai_factor, bseed_distr, wood_fract_min, do_alt_allometry
use vegn_tile_mod, only: vegn_tile_type, vegn_tile_carbon
use soil_tile_mod, only: num_l, dz, soil_tile_type, clw, csw
use vegn_cohort_mod, only : vegn_cohort_type, &
     update_biomass_pools, update_bio_living_fraction, update_species, &
     leaf_area_from_biomass, init_cohort_allometry_ppa, &
     cohort_root_litter_profile, cohort_root_exudate_profile
use vegn_disturbance_mod, only : kill_plants_ppa
use soil_carbon_mod, only: N_C_TYPES, soil_carbon_option, &
    SOILC_CENTURY, SOILC_CENTURY_BY_LAYER, SOILC_CORPSE, &
    add_litter
use soil_mod, only: add_soil_carbon, add_root_litter, add_root_exudates, Dsdt

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
character(len=*), parameter :: module_name = 'vegn_dynamics'
#include "../shared/version_variable.inc"
character(len=*), parameter :: diag_mod_name = 'vegn'

real, parameter :: GROWTH_RESP=0.333  ! fraction of NPP lost as growth respiration

! ==== module data ===========================================================
real    :: dt_fast_yr ! fast (physical) time step, yr (year is defined as 365 days)

! diagnostic field IDs
integer :: id_npp, id_nep, id_gpp
integer :: id_rsoil, id_rsoil_fast, id_rsoil_slow
integer :: id_resp, id_resl, id_resr, id_ress, id_resg, id_asoil
integer :: id_soilt, id_theta, id_litter, id_age


contains

! ============================================================================
subroutine vegn_dynamics_init(id_lon, id_lat, time, delta_time)
  integer        , intent(in) :: id_lon ! ID of land longitude (X) axis 
  integer        , intent(in) :: id_lat ! ID of land latitude (Y) axis
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
       (/id_lon,id_lat/), time, 'gross primary productivity', 'kg C/(m2 year)', &
       missing_value=-100.0)
  id_npp = register_cohort_diag_field ( diag_mod_name, 'npp',  &
       (/id_lon,id_lat/), time, 'net primary productivity', 'kg C/(m2 year)', &
       missing_value=-100.0)
  id_nep = register_tiled_diag_field ( diag_mod_name, 'nep',  &
       (/id_lon,id_lat/), time, 'net ecosystem productivity', 'kg C/(m2 year)', &
       missing_value=-100.0 )
  id_litter = register_tiled_diag_field (diag_mod_name, 'litter', (/id_lon,id_lat/), &
       time, 'litter productivity', 'kg C/(m2 year)', missing_value=-100.0)
  id_resp = register_cohort_diag_field ( diag_mod_name, 'resp', (/id_lon,id_lat/), &
       time, 'respiration', 'kg C/(m2 year)', missing_value=-100.0)
  id_resl = register_cohort_diag_field ( diag_mod_name, 'resl', (/id_lon,id_lat/), &
       time, 'leaf respiration', 'kg C/(m2 year)', missing_value=-100.0)
  id_resr = register_cohort_diag_field ( diag_mod_name, 'resr', (/id_lon,id_lat/), &
       time, 'root respiration', 'kg C/(m2 year)', missing_value=-100.0)
  id_ress = register_cohort_diag_field ( diag_mod_name, 'ress', (/id_lon,id_lat/), &
       time, 'stem respiration', 'kg C/(m2 year)', missing_value=-100.0)
  id_resg = register_cohort_diag_field ( diag_mod_name, 'resg', (/id_lon,id_lat/), &
       time, 'growth respiration', 'kg C/(m2 year)', missing_value=-100.0)
  id_soilt = register_tiled_diag_field ( diag_mod_name, 'tsoil_av',  &
       (/id_lon,id_lat/), time, 'average soil temperature for carbon decomposition', 'degK', &
       missing_value=-100.0 )
  id_theta = register_tiled_diag_field ( diag_mod_name, 'theta',  &
       (/id_lon,id_lat/), time, 'average soil wetness for carbon decomposition', 'm3/m3', &
       missing_value=-100.0 )
  id_age = register_cohort_diag_field ( diag_mod_name, 'age',  &
       (/id_lon,id_lat/), time, 'average cohort age', 'years', &
       missing_value=-100.0)
end subroutine vegn_dynamics_init


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
        resg(i) = GROWTH_RESP*cc%npp_previous_day
        npp(i)  = npp(i)  - GROWTH_RESP*cc%npp_previous_day
        resp(i) = resp(i) + GROWTH_RESP*cc%npp_previous_day
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
        md_leaf = cc%Pl * spdata(sp)%alpha(CMPT_LEAF)*cc%bliving*dt_fast_yr
        md_froot= cc%Pr * spdata(sp)%alpha(CMPT_ROOT)*cc%bliving*dt_fast_yr
     else
        md_leaf  = 0
        md_froot = 0
     endif
     
     ! compute branch and coarse wood losses for tree types
     if (spdata(sp)%lifeform==FORM_WOODY) then
        md_wood = 0.6 * cc%bwood * spdata(sp)%alpha(CMPT_WOOD)*dt_fast_yr
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
  
end subroutine vegn_carbon_int_lm3


! ============================================================================
subroutine vegn_carbon_int_ppa (vegn, soil, tsoil, theta, diag)
  ! TODO: possibly get rid of tsoil, theta, since they can be calculated here 
  type(vegn_tile_type), intent(inout) :: vegn
  type(soil_tile_type), intent(inout) :: soil
  real, intent(in) :: tsoil ! average temperature of soil for soil carbon decomposition, deg K
  real, intent(in) :: theta ! average soil wetness, unitless
  type(diag_buff_type), intent(inout) :: diag

  real :: md_wood;
  real :: deltaBL, deltaBR ! leaf and fine root carbon tendencies
  integer :: i, l, N
  real :: NSC_supply,LR_demand,LR_deficit
  real :: NSCtarget
  real :: R_days,fNSC,fLFR,fStem,fSeed
  type(vegn_cohort_type), pointer :: c(:) ! for debug only
  ! accumulators of total input to root litter and soil carbon
  real :: leaf_litt(N_C_TYPES) ! fine surface litter per tile, kgC/m2
  real :: wood_litt(N_C_TYPES) ! coarse surface litter per tile, kgC/m2
  real :: root_litt(num_l, N_C_TYPES) ! root litter per soil layer, kgC/m2
  real :: profile(num_l) ! storage for vertical profile of exudates and root litter
  real, dimension(vegn%n_cohorts) :: resp, resl, resr, ress, resg, gpp, npp

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
  leaf_litt = 0 ; wood_litt = 0; root_litt = 0 ! ; total_root_exudate_C = 0
  ! TODO: add root exudates to the balance. in LM3 it's a fraction of npp; 
  !       in PPA perhaps it might be proportional to NSC, with some decay
  !       time?
  do i = 1, vegn%n_cohorts
     associate ( cc => vegn%cohorts(i), sp => spdata(vegn%cohorts(i)%species))

     ! that was in eddy_npp_PPA
     call plant_respiration(cc,tsoil,resp(i),resl(i),resr(i),ress(i))
     gpp(i)  = (cc%An_op - cc%An_cl)*mol_C*cc%leafarea

!    A new scheme for plant growth, modified 9/3/2013 based on Steve's
!    suggestions
     R_days = 5.0
     fNSC=0.02
     fLFR=0.2
     fStem= dt_fast_yr/sp%tauNSC  ! 0.05
     fSeed= dt_fast_yr/tau_seed   ! 0.05
     NSCtarget = 4.0*cc%bl_max    ! 1.5*(cc%bl_max + cc%br_max)
     
!    '365.*dt_fast_yr/R_days' is '1/hours', eg. 1/120
!    the fraction filled per hour
     LR_demand = 0.0
     NSC_supply= 0.0

     IF(cc%nsc > 0 .AND. cc%status == LEAF_ON)then
         LR_deficit=max(cc%bl_max+cc%br_max-cc%bl-cc%br, 0.0)
         LR_demand =min(LR_deficit*fLFR*365.*dt_fast_yr, fNSC*cc%nsc)
         NSC_supply=max((cc%nsc-0.5*NSCtarget)*fStem,0.0) ! Weng 2014-01-23 for smoothing deltaDBH
         ! slm: signs are weird: LR_demand>0, NSC_supply<0
     endif

     resg(i) = GROWTH_RESP * (LR_demand + NSC_supply)/dt_fast_yr 
     cc%carbon_gain = cc%carbon_gain + (LR_demand + NSC_supply) 

     resp(i) = resp(i) + resg(i)
     npp(i)  = gpp(i) - resp(i)
     cc%nsc  = cc%nsc + gpp(i)*dt_fast_yr             &
                      - (LR_demand + NSC_supply)  &
                      - resp(i)*dt_fast_yr

     ! Weng, 2013-01-28
     ! Turnover regardless of STATUS
     deltaBL = cc%bl * sp%alpha(CMPT_LEAF) * dt_fast_yr
     deltaBR = cc%br * sp%alpha(CMPT_ROOT) * dt_fast_yr
     cc%bl = cc%bl - deltaBL
     cc%br = cc%br - deltaBR

     cc%carbon_loss = cc%carbon_loss + deltaBL + deltaBR  ! used in diagnostics only

     ! compute branch and coarse wood losses for tree types
     md_wood = 0.0
     if (spdata(cc%species)%lifeform == FORM_WOODY) then 
        ! md_wood = 0.6 *cc%bwood * sp%alpha(CMPT_WOOD)*dt_fast_yr
        ! set turnoverable wood biomass as a linear function of bl_max (max.foliage
        ! biomass) (Wang, Chuankuan 2006)
        md_wood = MIN(cc%bwood,cc%bl_max) * sp%alpha(CMPT_LEAF)*dt_fast_yr
     endif
     ! Why md_wood is set to 0?
     md_wood = 0.0

     ! accumulate liter and soil carbon inputs across all cohorts
     leaf_litt(:) = leaf_litt(:) + (/fsc_liv,  1-fsc_liv,  0.0/)*deltaBL*cc%nindivs
     wood_litt(:) = wood_litt(:) + (/fsc_wood, 1-fsc_wood, 0.0/)*md_wood*agf_bs*cc%nindivs
     call cohort_root_litter_profile(cc, dz, profile)
     do l = 1, num_l
        root_litt(l,:) = root_litt(l,:) + profile(l)*cc%nindivs*(/ &
             fsc_froot    *deltaBR + fsc_wood    *md_wood*(1-agf_bs), & ! fast
             (1-fsc_froot)*deltaBR + (1-fsc_wood)*md_wood*(1-agf_bs), & ! slow
             0.0/) ! microbes
     enddo

     ! increment cohort age
     cc%age = cc%age + dt_fast_yr
     ! resg(i) = 0 ! slm: that doesn't make much sense to me... why?
     end associate
  enddo 

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

  ! add litter accumulated over the cohorts
  call add_soil_carbon(soil, leaf_litt, wood_litt, root_litt)
  ! update soil carbon
  call Dsdt(vegn, soil, diag, tsoil, theta)

  ! NEP is equal to NPP minus soil respiration
  vegn%nep = sum(npp(1:N)*c(1:N)%nindivs) - vegn%rh

  call update_soil_pools(vegn, soil)
  vegn%age = vegn%age + dt_fast_yr;

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
  call send_cohort_data(id_age, diag, c(1:N), c(1:N)%age, weight=c(1:N)%nindivs, op=OP_MEAN)
end subroutine vegn_carbon_int_ppa


! ============================================================================
! updates cohort biomass pools, LAI, SAI, and height using accumulated 
! carbon_gain and bwood_gain
subroutine vegn_growth (vegn)
  type(vegn_tile_type), intent(inout) :: vegn

  ! ---- local vars
  type(vegn_cohort_type), pointer :: cc    ! current cohort
  real :: cmass0,cmass1 ! for debug only
  integer :: i

  if (is_watch_point()) then
     cmass0 = vegn_tile_carbon(vegn)
  endif

  do i = 1, vegn%n_cohorts   
     cc => vegn%cohorts(i)

     if (do_ppa) then
        call biomass_allocation_ppa(cc)
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
     endif

     if (cc%status == LEAF_ON) then
        cc%leaf_age = cc%leaf_age + 1.0
        ! limit the maximum leaf age by the leaf time span (reciprocal of leaf 
        ! turnover rate alpha) for given species. alpha is in 1/year, factor of
        ! 365 converts the result to days.
        if (spdata(cc%species)%alpha(CMPT_LEAF) > 0) &
             cc%leaf_age = min(cc%leaf_age,365.0/spdata(cc%species)%alpha(CMPT_LEAF))
     endif
  end do

  if (is_watch_point()) then
     cmass1 = vegn_tile_carbon(vegn)
     write(*,*)"############## vegn_growth #################"
     __DEBUG3__(cmass1-cmass0, cmass1, cmass0)
  endif
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
  real :: leaf_litt(N_C_TYPES) ! fine surface litter per tile, kgC/m2
  real :: wood_litt(N_C_TYPES) ! coarse surface litter per tile, kgC/m2
  real :: root_litt(num_l, N_C_TYPES) ! root litter per soil layer, kgC/m2
  real :: profile(num_l) ! storage for vertical profile of exudates and root litter

  leaf_litt = 0 ; wood_litt = 0; root_litt = 0
  do i = 1, vegn%n_cohorts
     associate ( cc => vegn%cohorts(i)   , &
                 sp => spdata(vegn%cohorts(i)%species)  )  ! F2003

    ! Mortality due to starvation
    if (cc%bsw<0 .or. cc%nsc < 0.01*cc%bl_max) then
       deathrate = 1.0

       deadtrees = min(cc%nindivs*deathrate,cc%nindivs) ! individuals / m2
       ! kill starved plants and add dead C from leaf and root pools to soil carbon
       call kill_plants_ppa(cc, vegn, soil, deadtrees, 0.0, leaf_litt, wood_litt, root_litt)

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

end subroutine vegn_starvation_ppa



! ==============================================================================
! updates cohort vegetation structure, biomass pools, LAI, SAI, and height  
! using accumulated carbon_gain
subroutine biomass_allocation_ppa(cc)
  type(vegn_cohort_type), intent(inout) :: cc
  
  real :: CSAtot ! total cross section area, m2
  real :: CSAsw  ! Sapwood cross sectional area, m2
  real :: CSAwd  ! Heartwood cross sectional area, m2
  real :: DBHwd  ! diameter of heartwood at breast height, m
  real :: BSWmax ! max sapwood biomass, kg C/individual
  real :: G_LFR  ! amount of carbon spent on leaf and root growth
  real :: deltaSeed
  real :: deltaBL, deltaBR ! tendencies of leaf and root biomass, kgC/individual
  real :: deltaBSW ! tendency of sapwood biomass, kgC/individual
  real :: deltaBwood ! tendency of wood biomass, kgC/individual
  real :: deltaDBH ! tendency of breast height diameter, m
  real :: deltaCA ! tendency of crown area, m2/individual
  real :: deltaHeight ! tendency of vegetation height
  real :: deltaCSAsw ! tendency of sapwood area
  real :: BL_c, BL_u ! canopy and understory target leaf biomasses, kgC/individual

  real, parameter :: DBHtp = 0.8 ! m
  logical :: hydraulics_repair = .FALSE.

  associate (sp => spdata(cc%species)) ! F2003

  ! TODO: what if carbon_gain is not 0, but leaves are OFF (marginal case? or
  ! typical in lm3?)
  if (cc%status == LEAF_ON) then
     ! calculate the carbon spent on growth of leaves and roots
     G_LFR    = max(0.0, min(cc%bl_max+cc%br_max-cc%bl-cc%br,  &
                            (1.-Wood_fract_min)*cc%carbon_gain))
     ! and distribute it between roots and leaves
     deltaBL  = min(G_LFR, max(0.0, &
          (G_LFR*cc%bl_max + cc%bl_max*cc%br - cc%br_max*cc%bl)/(cc%bl_max + cc%br_max) &
          ))
     deltaBR  = G_LFR - deltaBL
     ! calculate carbon spent on seeds and sapwood growth
     if (cc%layer == 1 .AND. cc%age > sp%maturalage) then
         deltaSeed=      sp%v_seed * (cc%carbon_gain - G_LFR)
         deltaBSW = (1.0-sp%v_seed)* (cc%carbon_gain - G_LFR)
     else
         deltaSeed= 0.0 
         deltaBSW = cc%carbon_gain - G_LFR
     endif
     ! update biomass pools due to growth
     cc%bl     = cc%bl    + deltaBL  ! updated in vegn_int_ppa
     cc%br     = cc%br    + deltaBR
     cc%bsw    = cc%bsw   + deltaBSW
     cc%bseed  = cc%bseed + deltaSeed

     ! calculate tendency of breast height diameter given increase of bsw
     deltaDBH     = deltaBSW / (sp%thetaBM * sp%alphaBM * cc%DBH**(sp%thetaBM-1))
     deltaHeight  = sp%thetaHT * sp%alphaHT * cc%DBH**(sp%thetaHT-1) * deltaDBH
!     deltaCA      = sp%thetaCA * sp%alphaCA * cc%DBH**(sp%thetaCA-1) * deltaDBH
     if(do_alt_allometry) then
         deltaCA  = sp%thetaCA * sp%alphaCA *     &
               (1./(exp(15.*(cc%DBH-DBHtp))+1.))  *     &     
               cc%DBH**(sp%thetaCA-1) * deltaDBH
     else
         deltaCA  = sp%thetaCA * sp%alphaCA * cc%DBH**(sp%thetaCA-1) * deltaDBH
     endif

     cc%DBH       = cc%DBH       + deltaDBH
     cc%height    = cc%height    + deltaHeight
     cc%crownarea = cc%crownarea + deltaCA

     ! calculate DBH, BLmax, BRmax, BSWmax using allometric relationships
     ! Weng 2012-01-31 update_bio_living_fraction
     CSAsw    = sp%alphaCSASW * cc%DBH**sp%thetaCSASW
     CSAtot   = PI * (cc%DBH/2.0)**2
     CSAwd    = max(0.0, CSAtot - CSAsw)
     DBHwd    = 2*sqrt(CSAwd/PI)
     BSWmax   = sp%alphaBM * (cc%DBH**sp%thetaBM - DBHwd**sp%thetaBM)
   
     ! Update Kxa, stem conductance if we are tracking past damage
     ! TODO: make hydraulics_repair a namelist parameter?
     if (.not.hydraulics_repair) then
        deltaCSAsw = CSAsw - (sp%alphaCSASW * (cc%DBH - deltaDBH)**sp%thetaCSASW)
        cc%Kxa = (cc%Kxa*CSAsw + sp%Kxam*deltaCSAsw)/(CSAsw + deltaCSAsw)
     endif

     deltaBwood = max(cc%bsw - BSWmax, 0.0)
     cc%bwood   = cc%bwood + deltaBwood
     cc%bsw     = cc%bsw   - deltaBwood

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
     cc%nsc = cc%nsc + cc%carbon_gain
  endif ! cc%status == LEAF_ON

  ! reset carbon accumulation terms
  cc%carbon_gain = 0
  cc%carbon_loss = 0
  
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
  r_vleaf = spdata(sp)%beta(CMPT_VLEAF) * cc%blv*tf;
  if (do_ppa) then
     ! Stem auto-respiration is proportional to cambium area, not sapwood biomass
     Acambium = PI * cc%DBH * cc%height * 1.2
     r_stem   = spdata(sp)%beta(CMPT_SAPWOOD) * Acambium * tf
  else
     r_stem   = spdata(sp)%beta(CMPT_SAPWOOD) * cc%bsw * tf
  endif
  r_root  = spdata(sp)%beta(CMPT_ROOT) * cc%br*tfs;

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
subroutine vegn_phenology_ppa(vegn, soil)
  type(vegn_tile_type), intent(inout) :: vegn
  type(soil_tile_type), intent(inout) :: soil

  ! ---- local vars
  integer :: i
  real    :: leaf_litter, leaf_litt(N_C_TYPES)
  real    :: leaf_fall, leaf_fall_rate ! per day
  real    :: root_mortality, root_mort_rate
  real    :: BL_u,BL_c

  leaf_fall_rate = 0.075
  root_mort_rate = 0.0
  
  vegn%litter = 0; leaf_litt(:) = 0.0
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
     if(cc%status == LEAF_OFF .AND. cc%bl > 0.)then
         leaf_fall = min(leaf_fall_rate * cc%bl_max, cc%bl)
         root_mortality = min( root_mort_rate* cc%br_max, cc%br)
         cc%nsc = cc%nsc + l_fract * (leaf_fall+ root_mortality)
         cc%bl  = cc%bl - leaf_fall
         cc%br  = cc%br - root_mortality
         cc%lai = leaf_area_from_biomass(cc%bl,cc%species,cc%layer,cc%firstlayer)/ &
                                        (cc%crownarea *(1.0-sp%internal_gap_frac))
         if(cc%bl == 0.)cc%leaf_age = 0.0
         cc%bliving = cc%blv + cc%br + cc%bl + cc%bsw

         leaf_litter = (1.-l_fract) * (leaf_fall+root_mortality) * cc%nindivs
         leaf_litt(:) = leaf_litt(:)+(/fsc_liv,1-fsc_liv,0.0/)*leaf_litter
         vegn%litter = vegn%litter + leaf_litter
         soil%fsc_in(1)  = soil%fsc_in(1)  + leaf_litter
         vegn%veg_out = vegn%veg_out + leaf_litter
     endif
     end associate
  enddo
  ! add litter accumulated over the cohorts
  call add_soil_carbon(soil, leaf_litter=leaf_litt)
end subroutine vegn_phenology_ppa


! ===========================================================================
! leaf falling at LEAF_OFF -- it is unused; why is it here? is there anything 
! in new Ensheng's code that uses it?
  subroutine vegn_leaf_fall_ppa(vegn,soil)
  type(vegn_tile_type), intent(inout) :: vegn
  type(soil_tile_type), intent(inout) :: soil

  ! ---- local vars
  integer :: i
  real    :: leaf_litt(N_C_TYPES)
  real    :: leaf_fall, leaf_fall_rate ! per day

  vegn%litter = 0 ! slm: why is it set to 0 here?
  leaf_fall_rate = 0.075 
  leaf_litt(:) = 0.0
  do i = 1,vegn%n_cohorts
     associate ( cc => vegn%cohorts(i),   &
                 sp => spdata(vegn%cohorts(i)%species) )
     if(cc%status == LEAF_OFF)then
        cc%nsc = cc%nsc + cc%carbon_gain
        leaf_fall = MIN(leaf_fall_rate * cc%bl_max, cc%bl)
        cc%nsc = cc%nsc + l_fract * leaf_fall
        cc%bl  = cc%bl - leaf_fall
        cc%lai = leaf_area_from_biomass(cc%bl,cc%species,cc%layer,cc%firstlayer)/ &
                                       (cc%crownarea *(1.0-sp%internal_gap_frac))
        if(cc%bl == 0.)cc%leaf_age = 0.0
        cc%bliving = cc%blv + cc%br + cc%bl + cc%bsw

        leaf_litt(:) = leaf_litt(:) + [fsc_liv,1-fsc_liv,0.0]*(1-l_fract) * leaf_fall * cc%nindivs
     endif
     end associate
  enddo
  vegn%litter  = vegn%litter  + sum(leaf_litt)
  vegn%veg_out = vegn%veg_out + sum(leaf_litt)
  call add_soil_carbon(soil, leaf_litter=leaf_litt)
end subroutine vegn_leaf_fall_ppa


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

     ! NOTE that this code looks very weird from the point of view of carbon
     ! balance: we are depleting the {fsc,ssc}_pool_ag, but adding to litter
     ! from the pools {leaflitter,coarsewoodLitter}_buffer_ag. The latter two
     ! are not used in the calculations of the total carbon.
     vegn%leaflitter_buffer_rate_ag = MAX(0.0, MIN(vegn%leaflitter_buffer_rate_ag, vegn%leaflitter_buffer_ag/dt_fast_yr))
     vegn%coarsewoodlitter_buffer_rate_ag = MAX(0.0, MIN(vegn%coarsewoodlitter_buffer_rate_ag, vegn%coarsewoodlitter_buffer_ag/dt_fast_yr))
     call add_litter(soil%leafLitter,(/vegn%leaflitter_buffer_rate_ag*dt_fast_yr*fsc_liv,vegn%leaflitter_buffer_rate_ag*dt_fast_yr*(1.0-fsc_liv),0.0/))
     call add_litter(soil%coarsewoodLitter,(/vegn%coarsewoodlitter_buffer_rate_ag*dt_fast_yr*fsc_liv,vegn%coarsewoodlitter_buffer_rate_ag*dt_fast_yr*(1.0-fsc_liv),0.0/))
  
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
  logical :: invasion = .FALSE.
  integer :: newcohorts ! number of new cohorts to be created
  integer :: i, k ! cohort indices
  real :: litt(N_C_TYPES)

! Check if reproduction happens
  newcohorts = 0
  do i = 1, vegn%n_cohorts
     if (cohort_can_reproduce(vegn%cohorts(i))) newcohorts=newcohorts + 1
  enddo

  if (newcohorts == 0) return ! do nothing if no cohorts are ready for reproduction

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
    cc         = parent

    cc%status  = LEAF_OFF
    cc%firstlayer = 0
    cc%age     = 0.0
    cc%bl      = sp%seedlingsize * bseed_distr(CMPT_LEAF)
    cc%br      = sp%seedlingsize * bseed_distr(CMPT_ROOT)
    cc%bsw     = sp%seedlingsize * bseed_distr(CMPT_SAPWOOD) ! sp%seedlingsize*0.1
    cc%bwood   = sp%seedlingsize * bseed_distr(CMPT_WOOD)    ! sp%seedlingsize*0.05
    cc%nsc     = sp%seedlingsize * bseed_distr(CMPT_NSC)     ! sp%seedlingsize*0.85
    cc%blv     = sp%seedlingsize * bseed_distr(CMPT_VLEAF)
    cc%bseed   = 0.0
    cc%topyear = 0.0
!   added germination probability (prob_g) and establishment probability ((prob_e), Weng 2014-01-06
    cc%nindivs = parent%bseed*parent%nindivs * sp%prob_g * sp%prob_e   &
                 /(sp%seedlingsize*sum(bseed_distr(:)))

    failed_seeds = (1.-sp%prob_g*sp%prob_e) * parent%bseed * parent%nindivs
    vegn%litter = vegn%litter + failed_seeds
    litt(:) = litt(:) + (/fsc_liv,1-fsc_liv,0.0/)*failed_seeds
    vegn%veg_out = vegn%veg_out + failed_seeds

    parent%bseed = 0.0

    call init_cohort_allometry_ppa(cc)
    cc%carbon_gain = 0.0
    cc%carbon_loss = 0.0
!    call biomass_allocation_ppa(cc)
    cc%bliving     = cc%br + cc%bl + cc%bsw + cc%blv
    cc%DBH_ys      = cc%DBH
    cc%BM_ys       = cc%bsw + cc%bwood
    cc%npp_previous_day     = 0.0
    cc%npp_previous_day_tmp = 0.0
    cc%leaf_age     = 0.0
    
    ! we assume that the newborn cohort is dry; since nindivs of the parent
    ! doesn't change we don't need to do anything with its Wl and Ws to 
    ! conserve water (since Wl and Ws are per individual)
    cc%Wl = 0 ; cc%Ws = 0
    ! TODO: make sure that energy is conserved in reproduction
    cc%Tv = parent%Tv
    
    end associate   ! F2003
  enddo
  call add_soil_carbon(soil, leaf_litter=litt)
  
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
  real :: leaf_litt(N_C_TYPES) ! fine surface litter per tile, kgC/m2
  real :: wood_litt(N_C_TYPES) ! coarse surface litter per tile, kgC/m2
  real :: root_litt(num_l, N_C_TYPES) ! root litter per soil layer, kgC/m2
  real :: profile(num_l) ! storage for vertical profile of exudates and root litter
  integer :: i,k

 ! Weng, 2013-09-07
 ! calculate the number of remaining cohorts
  k = 0
  do i = 1, vegn%n_cohorts
     if (vegn%cohorts(i)%nindivs >  mindensity) k=k+1
  enddo

  ! if (k==0) call error_mesg('vegn_mergecohorts_ppa','All cohorts died',WARNING)

  ! exclude cohorts that have zero individuals
  leaf_litt = 0 ; wood_litt = 0; root_litt = 0
  if (k < vegn%n_cohorts) then
     allocate(cc(max(k,1)))
     k=0
     do i = 1,vegn%n_cohorts
        if (vegn%cohorts(i)%nindivs > mindensity) then
           k=k+1
           cc(k) = vegn%cohorts(i)
        else
           call kill_plants_ppa(vegn%cohorts(i), vegn, soil, vegn%cohorts(i)%nindivs, 0.0, &
                                leaf_litt, wood_litt, root_litt)
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

  if (is_watch_point()) then
     write(*,*) '##### vegn_mergecohorts_ppa output #####'
     __DEBUG1__(vegn%n_cohorts)
     __DEBUG1__(vegn%cohorts%nindivs)
     __DEBUG1__(vegn%cohorts%Wl)
     __DEBUG1__(vegn%cohorts%Ws)
     __DEBUG1__(vegn%cohorts%mcv_dry)
     __DEBUG1__(vegn%cohorts%Tv)
  endif

end subroutine kill_small_cohorts_ppa

end module vegn_dynamics_mod
