! ============================================================================
! vegetation disturbances
! ============================================================================
module vegn_disturbance_mod

use land_constants_mod, only : seconds_per_year
use vegn_data_mod,   only : spdata, fsc_wood, fsc_liv, smoke_fraction, spdata, &
     agf_bs, LEAF_OFF
use vegn_tile_mod,   only : vegn_tile_type
use cohort_list_mod, only : vegn_cohort_enum_type, first_cohort, tail_cohort, &
     current_cohort, next_cohort, operator(/=)
use vegn_cohort_mod, only : vegn_cohort_type, height_from_biomass, lai_from_biomass, &
     update_bio_living_fraction, update_biomass_pools

implicit none
private

! ==== public interfaces =====================================================
public :: vegn_nat_mortality
public :: vegn_disturbance
public :: update_fuel
! =====end of public interfaces ==============================================

! ==== module constants ======================================================
character(len=*), parameter :: &
     version = '$Id: vegn_disturbance.F90,v 1.1.2.1 2007/09/16 22:14:24 slm Exp $', &
     tagname = '$Name: omsk_2007_12 $', &
     module_name = 'vegn_disturbance_mod'

contains ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

subroutine vegn_disturbance(vegn, dt)
  type(vegn_tile_type), intent(inout) :: vegn
  real, intent(in) :: dt ! time since last disturbance calculations, s
  

  real, parameter :: BMIN = 1e-10; ! should be the same as in growth function
  ! ---- local vars
  type(vegn_cohort_type), pointer :: cc    ! current cohort
  type(vegn_cohort_enum_type)     :: ci,ce ! cohort enumerators
  real :: area_disturbed_by_treefall;
  real :: precip;
  real :: delta;
  real :: fraction_lost;
  real :: bdead, balive; ! combined biomass pools
  real :: drought_month;
  real :: deltat

  deltat = dt/seconds_per_year ! convert time interval to years

  !  Disturbance Rates
  precip=vegn%p_ann*86400*365; 
  drought_month = vegn%lambda;
  vegn%disturbance_rate(0) = 0.0;
  vegn%disturbance_rate(1) = 0.0;
  
  call calculate_patch_disturbance_rates(vegn)
   
  ! Fire disturbance implicitly, i.e.  not patch creating
  vegn%area_disturbed_by_fire = (1.0-exp(-vegn%disturbance_rate(1)*deltat));
  
  ci = first_cohort(vegn%cohorts); ce = tail_cohort(vegn%cohorts);
  do while(ci/=ce)   
     cc => current_cohort(ci) ; ci = next_cohort(ci)

     fraction_lost = 1.0-exp(-vegn%disturbance_rate(1)*deltat);	
      
     ! "dead" biomass : wood + sapwood
     delta = (cc%bwood+cc%bsw)*fraction_lost;
      
     vegn%slow_soil_C = vegn%slow_soil_C + (1.0-smoke_fraction)*delta*(1-fsc_wood);
     vegn%fast_soil_C = vegn%fast_soil_C + (1.0-smoke_fraction)*delta*   fsc_wood;
     cc%bwood = cc%bwood * (1-fraction_lost);
     cc%bsw   = cc%bsw   * (1-fraction_lost);
      
     vegn%csmoke_pool = vegn%csmoke_pool + smoke_fraction*delta;
      
     ! for budget tracking - temporarily not keeping wood and the rest separately,ens
     !      vegn%ssc_in+=delta*(1.0-smoke_fraction)*(1-fsc_wood); */
     !      vegn%fsc_in+=delta*(1.0-smoke_fraction)*fsc_wood; */
      
     vegn%ssc_in = vegn%ssc_in+(cc%bwood+cc%bsw)*fraction_lost *(1.0-smoke_fraction);
     !     vegn%fsc_in+=cc%bsw*fraction_lost *(1.0-smoke_fraction);
     vegn%veg_out = vegn%veg_out+delta;
      
     !"alive" biomass: leaves, roots, and virtual pool
     delta = (cc%bl+cc%blv+cc%br)*fraction_lost;
     vegn%fast_soil_C = vegn%fast_soil_C + (1.0-smoke_fraction)*delta*    fsc_liv ;
     vegn%slow_soil_C = vegn%slow_soil_C + (1.0-smoke_fraction)*delta*(1- fsc_liv);
      
     cc%bl  = cc%bl  * (1-fraction_lost);
     cc%blv = cc%blv * (1-fraction_lost);
     cc%br  = cc%br  * (1-fraction_lost);
      
     vegn%csmoke_pool = vegn%csmoke_pool + smoke_fraction*delta;
      
     ! for budget tracking- temporarily keeping alive separate ens
     ! /*      vegn%fsc_in+=delta* fsc_liv; */
     ! /*      vegn%ssc_in+=delta* (1-fsc_liv); */
     vegn%fsc_in = vegn%fsc_in+delta*(1.0-smoke_fraction);
     vegn%veg_out = vegn%veg_out+delta;
      
     !"living" biomass:leaves, roots and sapwood
     delta = cc%bliving*fraction_lost;
     cc%bliving = cc%bliving - delta;

     if(cc%bliving < BMIN) then
        ! remove vegetaion competely 	      
        vegn%fast_soil_C = vegn%fast_soil_C + fsc_liv*cc%bliving+ fsc_wood*cc%bwood;
        vegn%slow_soil_C = vegn%slow_soil_C + (1.- fsc_liv)*cc%bliving+ (1-fsc_wood)*cc%bwood;
        
        vegn%fsc_in = vegn%fsc_in + cc%bwood+cc%bliving;
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
  type(vegn_cohort_type), pointer :: cc

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
  cc => current_cohort(first_cohort(vegn%cohorts))
  vegn%disturbance_rate(0) = spdata(cc%species)%treefall_disturbance_rate;
  vegn%total_disturbance_rate = vegn%disturbance_rate(1)+vegn%disturbance_rate(0);
  
  vegn%fuel = fuel;
end subroutine calculate_patch_disturbance_rates


! ============================================================================
function fire(vegn) result(fireterm)
  real :: fireterm; ! return value
  type(vegn_tile_type), intent(inout) :: vegn

  ! ---- local vars
  real :: precip_av ! average precipitation, mm/year

  fireterm = 0
  precip_av = vegn%p_ann * seconds_per_year;
!!$  vegn%ignition_rate = 0.00;             
  vegn%fuel = vegn%total_biomass;

  if(vegn%fuel>0.0) then
     if(precip_av < 400.+40.*(vegn%t_ann-273.16)) then
        fireterm = vegn%fuel*(400. + 40.*(vegn%t_ann-273.16) - precip_av);
     endif
  endif
end function


! ============================================================================
subroutine update_fuel(vegn)
  type(vegn_tile_type), intent(inout) :: vegn

  ! ---- local constants
  !  these three used to be in data 
  real, parameter :: fire_height_threashold = 100;
  real, parameter :: k_fw = 0.1;
  real, parameter :: fp1 = 1.; ! disturbance rate per kgC/m2 of fuel
  ! ---- local vars
  type(vegn_cohort_type), pointer :: cc    ! current cohort
  type(vegn_cohort_enum_type)     :: ci,ce ! cohort enumerators
  real :: ignition_rate;
  real ::  babove;

  ci = first_cohort(vegn%cohorts); ce = tail_cohort(vegn%cohorts);
  do while(ci/=ce)   
     cc => current_cohort(ci) ; ci = next_cohort(ci)
     if((cc%height < fire_height_threashold).and.(vegn%theta_av < k_fw).and.(vegn%tsoil_av>278.16)) then
        babove = cc%bl + agf_bs * (cc%bsw + cc%bwood + cc%blv);
        ! this is fuel available durng the drought months only
        vegn%fuel = vegn%fuel + spdata(cc%species)%fuel_intensity*babove;	
     endif
  enddo

  ignition_rate = 0.;
  if (vegn%theta_av<k_fw .and. vegn%tsoil_av>278.16) ignition_rate = 1.;
  vegn%lambda = vegn%lambda + ignition_rate;

end subroutine update_fuel




! ============================================================================
subroutine vegn_nat_mortality(vegn, deltat)
  type(vegn_tile_type), intent(inout) :: vegn
  real, intent(in) :: deltat ! time since last mortality calculations, s
  
  ! ---- local vars
  type(vegn_cohort_type), pointer :: cc    ! current cohort
  type(vegn_cohort_enum_type)     :: ci,ce ! cohort enumerators
  real :: delta;
  real :: fraction_lost;
  real :: bdead, balive; ! combined biomass pools
  
  vegn%disturbance_rate(0)        = 0.0; 
  vegn%area_disturbed_by_treefall = 0.0;
  
  ci = first_cohort(vegn%cohorts); ce = tail_cohort(vegn%cohorts);
  do while(ci/=ce)   
     cc => current_cohort(ci) ; ci = next_cohort(ci)
     ! Treat treefall disturbance implicitly, i.e. not creating a new tile.
     ! note that this disturbance rate calculation only works for the one cohort per 
     ! tile case -- in case of multiple cohort disturbance rate pehaps needs to be 
     ! accumulated (or averaged? or something else?) over the cohorts.
     vegn%disturbance_rate(0) = spdata(cc%species)%treefall_disturbance_rate;
     vegn%area_disturbed_by_treefall = &
          1.0-exp(-vegn%disturbance_rate(0)*deltat/seconds_per_year);

     ! calculate combined biomass pools
     balive = cc%bl + cc%blv + cc%br;
     bdead  = cc%bsw + cc%bwood;
     ! ens need a daily PATCH_FREQ here, for now it is set to 48
     fraction_lost = 1.0-exp(-vegn%disturbance_rate(0)*deltat/seconds_per_year);     
      
     ! "dead" biomass : wood + sapwood
     delta = bdead*fraction_lost;

     vegn%slow_soil_C = vegn%slow_soil_C + (1-fsc_wood)*delta;
     vegn%fast_soil_C = vegn%fast_soil_C +    fsc_wood *delta;

     cc%bwood = cc%bwood * (1-fraction_lost);
     cc%bsw   = cc%bsw   * (1-fraction_lost);

     ! for budget tracking -temporarily
     ! vegn%fsc_in+= cc%bsw*fraction_lost;
     vegn%ssc_in  = vegn%ssc_in  + (cc%bwood+cc%bsw)*fraction_lost;
     vegn%veg_out = vegn%veg_out + delta;
     
     ! note that fast "living" pools are not included into mortality because their 
     ! turnover is calculated separately

     cc%bliving = cc%bsw + cc%bl + cc%br + cc%blv;
     call update_biomass_pools(cc);
  enddo
     
end subroutine vegn_nat_mortality


end module vegn_disturbance_mod
