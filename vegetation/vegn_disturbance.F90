! ============================================================================
! vegetation disturbances
! ============================================================================
module vegn_disturbance_mod

#include "../shared/debug.inc"

use fms_mod,         only : error_mesg, WARNING
use constants_mod,   only : tfreeze
use land_constants_mod, only : seconds_per_year
use land_debug_mod,  only : is_watch_point
use vegn_data_mod,   only : spdata, fsc_wood, fsc_liv, agf_bs, do_ppa, LEAF_OFF, &
                            DBH_mort, A_mort, B_mort, mortrate_s
use vegn_tile_mod,   only : vegn_tile_type
use soil_tile_mod,   only : soil_tile_type
use vegn_cohort_mod, only : vegn_cohort_type, update_biomass_pools

implicit none
private

! ==== public interfaces =====================================================
public :: vegn_nat_mortality_lm3
public :: vegn_nat_mortality_ppa
public :: vegn_disturbance
public :: update_fuel
! =====end of public interfaces ==============================================

! ==== module constants ======================================================
character(len=*), parameter :: &
     version = '$Id$', &
     tagname = '$Name$', &
     module_name = 'vegn_disturbance_mod'
! TODO: possibly move all definition of cpw,clw,csw in one place
real, parameter :: &
     cpw = 1952.0, & ! specific heat of water vapor at constant pressure
     clw = 4218.0, & ! specific heat of water (liquid)
     csw = 2106.0    ! specific heat of water (ice)

contains ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

subroutine vegn_disturbance(vegn, soil, dt)
  type(vegn_tile_type), intent(inout) :: vegn ! vegetation data
  type(soil_tile_type), intent(inout) :: soil ! soil data
  real, intent(in) :: dt ! time since last disturbance calculations, s
  
  real, parameter :: BMIN = 1e-10; ! should be the same as in growth function
  ! ---- local vars
  type(vegn_cohort_type), pointer :: cc    ! current cohort
  real :: precip;
  real :: delta;
  real :: fraction_lost;
  real :: drought_month;
  real :: deltat
  integer :: i
  integer :: sp ! shorthand for cohort species

  deltat = dt/seconds_per_year ! convert time interval to years

  !  Disturbance Rates
  precip=vegn%p_ann*86400*365; 
  drought_month = vegn%lambda;
  vegn%disturbance_rate(0) = 0.0;
  vegn%disturbance_rate(1) = 0.0;
  
  call calculate_patch_disturbance_rates(vegn)
   
  do i = 1,vegn%n_cohorts   
     cc => vegn%cohorts(i)
     sp = cc%species

     fraction_lost = 1.0-exp(-vegn%disturbance_rate(1)*deltat);	
      
     ! "dead" biomass : wood + sapwood
     delta = (cc%bwood+cc%bsw)*fraction_lost;
      
     soil%slow_soil_C(1) = soil%slow_soil_C(1) + (1.0-spdata(sp)%smoke_fraction)*delta*(1-fsc_wood);
     soil%fast_soil_C(1) = soil%fast_soil_C(1) + (1.0-spdata(sp)%smoke_fraction)*delta*   fsc_wood;
     cc%bwood = cc%bwood * (1-fraction_lost);
     cc%bsw   = cc%bsw   * (1-fraction_lost);
      
     vegn%csmoke_pool = vegn%csmoke_pool + spdata(sp)%smoke_fraction*delta;
      
     ! for budget tracking - temporarily not keeping wood and the rest separately,ens
     !      soil%ssc_in(1)+=delta*(1.0-spdata(sp)%smoke_fraction)*(1-fsc_wood); */
     !      soil%fsc_in(1)+=delta*(1.0-spdata(sp)%smoke_fraction)*fsc_wood; */
      
     soil%ssc_in(1) = soil%ssc_in(1)+(cc%bwood+cc%bsw)*fraction_lost *(1.0-spdata(sp)%smoke_fraction);
     !     soil%fsc_in(1)+=cc%bsw*fraction_lost *(1.0-spdata(sp)%smoke_fraction);
     vegn%veg_out = vegn%veg_out+delta;
      
     !"alive" biomass: leaves, roots, and virtual pool
     delta = (cc%bl+cc%blv+cc%br)*fraction_lost;
     soil%fast_soil_C(1) = soil%fast_soil_C(1) + (1.0-spdata(sp)%smoke_fraction)*delta*    fsc_liv ;
     soil%slow_soil_C(1) = soil%slow_soil_C(1) + (1.0-spdata(sp)%smoke_fraction)*delta*(1- fsc_liv);
      
     cc%bl  = cc%bl  * (1-fraction_lost);
     cc%blv = cc%blv * (1-fraction_lost);
     cc%br  = cc%br  * (1-fraction_lost);
      
     vegn%csmoke_pool = vegn%csmoke_pool + spdata(sp)%smoke_fraction*delta;
      
     ! for budget tracking- temporarily keeping alive separate ens
     ! /*      soil%fsc_in(1)+=delta* fsc_liv; */
     ! /*      soil%ssc_in(1)+=delta* (1-fsc_liv); */
     soil%fsc_in(1) = soil%fsc_in(1)+delta*(1.0-spdata(sp)%smoke_fraction);
     vegn%veg_out = vegn%veg_out+delta;
      
     !"living" biomass:leaves, roots and sapwood
     delta = cc%bliving*fraction_lost;
     cc%bliving = cc%bliving - delta;

     if(cc%bliving < BMIN) then
        ! remove vegetaion competely 	      
        soil%fast_soil_C(1) = soil%fast_soil_C(1) + fsc_liv*cc%bliving+ fsc_wood*cc%bwood;
        soil%slow_soil_C(1) = soil%slow_soil_C(1) + (1.- fsc_liv)*cc%bliving+ (1-fsc_wood)*cc%bwood;
        
        soil%fsc_in(1) = soil%fsc_in(1) + cc%bwood+cc%bliving;
        vegn%veg_out = vegn%veg_out + cc%bwood+cc%bliving;
        
        cc%bliving = 0.;
        cc%bwood   = 0.;
     endif
     call update_biomass_pools(cc)
  enddo

  vegn%csmoke_rate = vegn%csmoke_pool; ! kg C/(m2 yr)
end subroutine vegn_disturbance

! ============================================================================
subroutine calculate_patch_disturbance_rates(vegn)
  type(vegn_tile_type), intent(inout) :: vegn

  ! ---- local vars
  real :: fire_prob;
  real :: fuel;

  fuel = vegn%fuel

#if SIMPLE_FIRE
  ! CALCULATE FIRE DISTURBANCE RATES
  vegn%disturbance_rate(1)=fire(vegn);
#else

  ! lambda is the number of drought months;
  fire_prob = vegn%lambda/(1.+vegn%lambda); 
  ! compute average fuel during fire months
  if (vegn%lambda > 0.00001 ) fuel = fuel/vegn%lambda;
  vegn%disturbance_rate(1) = fuel * fire_prob;

  ! put a threshold for very dry years for warm places
  if (vegn%t_ann > 273.16 .and. vegn%lambda > 3.)  vegn%disturbance_rate(1)=0.33;
#endif  
  
  if(vegn%disturbance_rate(1) > 0.33) vegn%disturbance_rate(1)=0.33;
  
  ! this is only true for the one cohort per patch case
  vegn%disturbance_rate(0) = spdata(vegn%cohorts(1)%species)%treefall_disturbance_rate;
  
  vegn%fuel = fuel;
end subroutine calculate_patch_disturbance_rates


! ============================================================================
subroutine update_fuel(vegn, wilt)
  type(vegn_tile_type), intent(inout) :: vegn
  real, intent(in) :: wilt ! ratio of wilting to saturated water content

  ! ---- local constants
  !  these three used to be in data 
  real, parameter :: fire_height_threashold = 100;
  real, parameter :: fp1 = 1.; ! disturbance rate per kgC/m2 of fuel
  ! ---- local vars
  type(vegn_cohort_type), pointer :: cc    ! current cohort
  real :: theta_crit; ! critical ratio of average soil water to sat. water
  real :: ignition_rate;
  real ::  babove;
  integer :: i

  do i = 1,vegn%n_cohorts   
     cc => vegn%cohorts(i)
     ! calculate theta_crit: actually either fact_crit_fire or cnst_crit_fire 
     ! is zero, enforced by logic in the vegn_data.F90
     theta_crit = spdata(cc%species)%cnst_crit_fire &
           + wilt*spdata(cc%species)%fact_crit_fire
     theta_crit = max(0.0,min(1.0, theta_crit))
     if((cc%height < fire_height_threashold) &
          .and.(vegn%theta_av_fire < theta_crit)  &
          .and.(vegn%tsoil_av > 278.16)) then
        babove = cc%bl + agf_bs * (cc%bsw + cc%bwood + cc%blv);
        ! this is fuel available during the drought months only
        vegn%fuel = vegn%fuel + spdata(cc%species)%fuel_intensity*babove;	
     endif
  enddo

  ! note that the ignition rate calculation based on the value of theta_crit for 
  ! the last cohort -- currently it doesn't matter since we have just one cohort, 
  ! but something needs to be done about that in the future
  ignition_rate = 0.;
  if ( (vegn%theta_av_fire < theta_crit) &
       .and. (vegn%tsoil_av>278.16)) ignition_rate = 1.;
  vegn%lambda = vegn%lambda + ignition_rate;

end subroutine update_fuel


! ============================================================================
subroutine vegn_nat_mortality_lm3(vegn, soil, deltat)
  type(vegn_tile_type), intent(inout) :: vegn  ! vegetation data
  type(soil_tile_type), intent(inout) :: soil  ! soil data
  real, intent(in) :: deltat ! time since last mortality calculations, s
  
  ! ---- local vars
  type(vegn_cohort_type), pointer :: cc    ! current cohort
  real :: delta;
  real :: fraction_lost;
  real :: bdead, balive; ! combined biomass pools
  integer :: i
  
  vegn%disturbance_rate(0)        = 0.0; 
  
  do i = 1,vegn%n_cohorts
     cc => vegn%cohorts(i)
     ! Treat treefall disturbance implicitly, i.e. not creating a new tile.
     ! note that this disturbance rate calculation only works for the one cohort per 
     ! tile case -- in case of multiple cohort disturbance rate perhaps needs to be 
     ! accumulated (or averaged? or something else?) over the cohorts.
     vegn%disturbance_rate(0) = spdata(cc%species)%treefall_disturbance_rate;

     ! calculate combined biomass pools
     balive = cc%bl + cc%blv + cc%br;
     bdead  = cc%bsw + cc%bwood;
     ! ens need a daily PATCH_FREQ here, for now it is set to 48
     fraction_lost = 1.0-exp(-vegn%disturbance_rate(0)*deltat/seconds_per_year);     
      
     ! "dead" biomass : wood + sapwood
     delta = bdead*fraction_lost;

     soil%slow_soil_C(1) = soil%slow_soil_C(1) + (1-fsc_wood)*delta;
     soil%fast_soil_C(1) = soil%fast_soil_C(1) +    fsc_wood *delta;

     cc%bwood = cc%bwood * (1-fraction_lost);
     cc%bsw   = cc%bsw   * (1-fraction_lost);

     ! for budget tracking -temporarily
     ! soil%fsc_in(1)+= cc%bsw*fraction_lost;
     soil%ssc_in(1)  = soil%ssc_in(1)  + (cc%bwood+cc%bsw)*fraction_lost;
     vegn%veg_out = vegn%veg_out + delta;
     
     ! note that fast "living" pools are not included into mortality because their 
     ! turnover is calculated separately

     cc%bliving = cc%bsw + cc%bl + cc%br + cc%blv;
     call update_biomass_pools(cc);
  enddo
     
end subroutine vegn_nat_mortality_lm3


! ============================================================================
! TODO: spread the mortality input to fsc and ssc over a year, to avoid spikes 
! in carbon fluxes
subroutine vegn_nat_mortality_ppa (vegn, soil, deltat)
  type(vegn_tile_type), intent(inout) :: vegn
  type(soil_tile_type), intent(inout) :: soil
  real, intent(in) :: deltat ! time since last mortality calculations, s

  ! ---- local vars
  real :: loss_alive,loss_wood
  real :: deathrate ! mortality rate, 1/year
  real :: deadtrees ! number of trees that died over the time step
  integer :: i, k
!  type(vegn_cohort_type), pointer :: cc(:) ! array to hold new cohorts
  
  
  real, parameter :: min_nindivs = 1e-5 ! 1/m2. If nindivs is less that this number, 
  ! then the entire cohort is killed; 2e-15 is approximately 1 individual per Earth 
  ! surface area
  logical :: do_U_shaped_mortality = .FALSE.
  real :: DBHtp ! for U-shaped mortality

  DBHtp = 1.0  ! for U-shaped mortality

  do i = 1, vegn%n_cohorts   
     associate ( cc => vegn%cohorts(i)   , &
                 sp => spdata(vegn%cohorts(i)%species)  )
     ! mortality rate can be a function of growth rate, age, and environmental
     ! conditions. Here, we only used two constants for canopy layer and under-
     ! story layer (mortrate_d_c and mortrate_d_u)
     if(cc%layer > 1) then
        deathrate = sp%mortrate_d_u * &
                 (1 + A_mort*exp(B_mort*(DBH_mort-cc%dbh)) &
                    /(1.0 + exp(B_mort*(DBH_mort-cc%dbh))) &
                 )
     else
        if(do_U_shaped_mortality)then
            deathrate = sp%mortrate_d_c *                    &
                       (1.0 + 12.*exp(8.*(cc%dbh-DBHtp))     &
                             /(1.0 + exp(8.*(cc%dbh-DBHtp))) &
                       )
        else
            deathrate = sp%mortrate_d_c
        endif
     endif

     deadtrees = cc%nindivs * (1.0-exp(-deathrate*deltat/seconds_per_year)) ! individuals / m2
     cc%nindivs = cc%nindivs-deadtrees
     ! water from dead trees goes to intermediate buffers, to be added to the
     ! precipitation reaching ground
     vegn%drop_wl = vegn%drop_wl + cc%wl*deadtrees
     vegn%drop_ws = vegn%drop_ws + cc%ws*deadtrees
     vegn%drop_hl = vegn%drop_hl + clw*cc%wl*deadtrees*(cc%Tv-tfreeze)
     vegn%drop_hs = vegn%drop_hs + csw*cc%ws*deadtrees*(cc%Tv-tfreeze)

     ! add dead C from leaf and root pools to fast soil carbon
     loss_wood  = deadtrees * cc%bwood
     loss_alive = deadtrees * (cc%bl+cc%br+cc%bsw+cc%blv+cc%bseed+cc%nsc)
     soil%fast_soil_C(1) = soil%fast_soil_C(1) +    fsc_liv *loss_alive +    fsc_wood *loss_wood
     soil%slow_soil_C(1) = soil%slow_soil_C(1) + (1-fsc_liv)*loss_alive + (1-fsc_wood)*loss_wood
          
     ! for budget tracking - temporary
     soil%ssc_in(1) = soil%ssc_in(1) + fsc_liv*loss_alive + fsc_wood *loss_wood
     vegn%veg_out = vegn%veg_out + loss_alive + loss_wood

     end associate
     ! note that in contrast to LM3 mortality calculation, living biomasses are
     ! included in loss; does it mean that we double-counting losses in 
     ! maintenance and here?
     ! TODO: ask Ensheng and Elena about possible double-counting of live biomass losses.
  enddo

end subroutine vegn_nat_mortality_ppa


end module vegn_disturbance_mod
