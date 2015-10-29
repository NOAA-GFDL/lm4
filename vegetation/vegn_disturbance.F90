! ============================================================================
! vegetation disturbances
! ============================================================================
module vegn_disturbance_mod

#include "../shared/debug.inc"

use fms_mod,         only : error_mesg, WARNING, FATAL
use time_manager_mod,only : time_type, get_date, operator(-)
use constants_mod,   only : tfreeze
use land_constants_mod, only : seconds_per_year
use land_debug_mod,  only : is_watch_point, check_var_range, set_current_point, &
     check_conservation, do_check_conservation, water_cons_tol, carbon_cons_tol
use vegn_data_mod,   only : spdata, fsc_wood, fsc_liv, fsc_froot, agf_bs, &
       do_ppa, LEAF_OFF, DBH_mort, A_mort, B_mort, mortrate_s, nat_mortality_splits_tiles
use vegn_tile_mod,   only : vegn_tile_type, relayer_cohorts
use soil_tile_mod,   only : soil_tile_type, num_l, dz
use land_tile_mod,   only : land_tile_type, land_tile_enum_type, &
     land_tile_list_type, land_tile_list_init, land_tile_list_end, &
     empty, first_elmt, tail_elmt, next_elmt, merge_land_tile_into_list, &
     current_tile, operator(==), operator(/=), remove, insert, new_land_tile, &
     land_tile_heat, land_tile_carbon, get_tile_water
use land_data_mod,   only : lnd, land_time
use soil_mod,        only : add_soil_carbon
use soil_carbon_mod, only : add_litter, soil_carbon_option, &
     SOILC_CENTURY, SOILC_CENTURY_BY_LAYER, SOILC_CORPSE, N_C_TYPES, C_CEL
use vegn_cohort_mod, only : vegn_cohort_type, update_biomass_pools, &
     cohort_root_litter_profile

implicit none
private

! ==== public interfaces =====================================================
public :: vegn_nat_mortality_lm3
public :: vegn_nat_mortality_ppa
public :: vegn_disturbance
public :: update_fuel
public :: kill_plants_ppa
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
  real :: precip;
  real :: delta;
  real :: fraction_lost;
  real :: drought_month;
  real :: deltat
  integer :: i,l
  real :: leaf_litt(N_C_TYPES) ! fine surface litter per tile, kgC/m2
  real :: wood_litt(N_C_TYPES) ! coarse surface litter per tile, kgC/m2
  real :: root_litt(num_l, N_C_TYPES) ! root litter per soil layer, kgC/m2
  real :: profile(num_l) ! storage for vertical profile of exudates and root litter

  deltat = dt/seconds_per_year ! convert time interval to years

  !  Disturbance Rates
  precip=vegn%p_ann*86400*365; 
  drought_month = vegn%lambda;
  vegn%disturbance_rate(0) = 0.0;
  vegn%disturbance_rate(1) = 0.0;
  
  call calculate_patch_disturbance_rates(vegn)
  
  leaf_litt = 0 ; wood_litt = 0; root_litt = 0
  do i = 1,vegn%n_cohorts   
     associate (cc => vegn%cohorts(i), &
                sp => spdata(vegn%cohorts(i)%species)) ! F2003

     fraction_lost = 1.0-exp(-vegn%disturbance_rate(1)*deltat);	
     if (do_ppa) then
        call kill_plants_ppa(cc, vegn, soil, cc%nindivs*fraction_lost, sp%smoke_fraction, &
                             leaf_litt, wood_litt, root_litt)
     else ! original LM3 treatment
        ! "dead" biomass : wood + sapwood
        delta = (cc%bwood+cc%bsw)*fraction_lost*cc%nindivs

        ! accumulate liter and soil carbon inputs across all cohorts
        wood_litt(:) = wood_litt(:) + [fsc_wood, 1-fsc_wood, 0.0]*agf_bs*delta*(1-sp%smoke_fraction)
        call cohort_root_litter_profile(cc, dz, profile)
        do l = 1, num_l
           root_litt(l,:) = root_litt(l,:) + [fsc_wood, 1-fsc_wood, 0.0] * &
                            profile(l)*delta*(1-agf_bs)*(1-sp%smoke_fraction)
        enddo
      
        cc%bwood = cc%bwood * (1-fraction_lost);
        cc%bsw   = cc%bsw   * (1-fraction_lost);
      
        vegn%csmoke_pool = vegn%csmoke_pool + sp%smoke_fraction*delta;
      
        ! for budget tracking - temporarily not keeping wood and the rest separately,ens
        !      soil%ssc_in(1)+=delta*(1.0-sp%smoke_fraction)*(1-fsc_wood); */
        !      soil%fsc_in(1)+=delta*(1.0-sp%smoke_fraction)*fsc_wood; */
      
        vegn%veg_out = vegn%veg_out+delta;
      
        !"alive" biomass: leaves, roots, and virtual pool
        delta = (cc%bl+cc%blv+cc%br)*fraction_lost;
        leaf_litt(:) = leaf_litt(:) + [fsc_liv, 1-fsc_liv, 0.0]*(cc%bl+cc%blv)*fraction_lost*(1-sp%smoke_fraction)
        do l = 1, num_l
           root_litt(l,:) = root_litt(l,:) + [fsc_froot, 1-fsc_froot, 0.0] * &
                            profile(l)*cc%br*fraction_lost*(1-sp%smoke_fraction)
        enddo
      
        cc%bl  = cc%bl  * (1-fraction_lost);
        cc%blv = cc%blv * (1-fraction_lost);
        cc%br  = cc%br  * (1-fraction_lost);
      
        vegn%csmoke_pool = vegn%csmoke_pool + sp%smoke_fraction*delta;
      
        vegn%veg_out = vegn%veg_out+delta;
      
        !"living" biomass:leaves, roots and sapwood
        delta = cc%bliving*fraction_lost;
        cc%bliving = cc%bliving - delta;

        if(cc%bliving < BMIN) then
           ! remove vegetaion competely 	      
           ! accumulate liter and soil carbon inputs across all cohorts
           leaf_litt(:) = leaf_litt(:) + [fsc_liv,  1-fsc_liv,  0.0]*(cc%bl+cc%blv)
           wood_litt(:) = wood_litt(:) + [fsc_wood, 1-fsc_wood, 0.0]*agf_bs*(cc%bwood+cc%bsw)
           do l = 1, num_l
              root_litt(l,:) = root_litt(l,:) + profile(l)*[ &
                   fsc_froot    *cc%br + fsc_wood    *(cc%bwood+cc%bsw)*(1-agf_bs), & ! fast
                   (1-fsc_froot)*cc%br + (1-fsc_wood)*(cc%bwood+cc%bsw)*(1-agf_bs), & ! slow
                   0.0] ! microbes
           enddo

           vegn%veg_out = vegn%veg_out + cc%bwood+cc%bliving;
        
           cc%bliving = 0.;
           cc%bwood   = 0.;
        endif
        call update_biomass_pools(cc)
     endif ! LM3 treatment   
     end associate
  enddo

  call add_soil_carbon(soil, leaf_litt, wood_litt, root_litt)

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
  integer :: i,l
  real :: wood_litt(n_c_types) ! coarse surface litter per tile, kgC/m2
  real :: root_litt(num_l, n_c_types) ! root litter per soil layer, kgC/m2
  real :: profile(num_l) ! storage for vertical profile of exudates and root litter
  
  vegn%disturbance_rate(0) = 0.0; 
  wood_litt = 0; root_litt = 0
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

     wood_litt(:) = wood_litt(:) + [fsc_wood, 1-fsc_wood, 0.0]*delta*agf_bs
     call cohort_root_litter_profile(cc, dz, profile)
     do l = 1, num_l
        root_litt(l,:) = root_litt(l,:) + profile(l)*(1-agf_bs)*delta*[fsc_wood, 1-fsc_wood, 0.0]
     enddo

     cc%bwood = cc%bwood * (1-fraction_lost);
     cc%bsw   = cc%bsw   * (1-fraction_lost);

     vegn%veg_out = vegn%veg_out + delta;
     
     ! note that fast "living" pools are not included into mortality because their 
     ! turnover is calculated separately

     cc%bliving = cc%bsw + cc%bl + cc%br + cc%blv;
     call update_biomass_pools(cc);
  enddo
  ! add litter accumulated over the cohorts
  call add_soil_carbon(soil, wood_litter=wood_litt, root_litter=root_litt)     
end subroutine vegn_nat_mortality_lm3


! ============================================================================
! goes through all tiles in the domain, and applies natural mortality to
! each of them, possibly creating new tiles within affected grid cells.
subroutine vegn_nat_mortality_ppa ( )

  integer :: second, minute, hour, day0, day1, month0, month1, year0, year1
  type(land_tile_type), pointer :: t0,t1
  type(land_tile_enum_type) :: ts, te
  type(land_tile_list_type) :: disturbed_tiles
  integer :: i,j,k
  real, allocatable :: ndead(:)

  ! get components of calendar dates for this and previous time step
  call get_date(land_time,             year0,month0,day0,hour,minute,second)
  call get_date(land_time-lnd%dt_slow, year1,month1,day1,hour,minute,second)

  if (.not.(do_ppa.and.year1 /= year0)) return  ! do nothing 

  call land_tile_list_init(disturbed_tiles)

  do j = lnd%js,lnd%je
  do i = lnd%is,lnd%ie
     if(empty(lnd%tile_map(i,j))) cycle ! skip cells where there is no land
     ! set current point for debugging
     call set_current_point(i,j,1)
     ts = first_elmt(lnd%tile_map(i,j)) ; te=tail_elmt(lnd%tile_map(i,j))
     do while (ts /= te)
        t0=>current_tile(ts); ts=next_elmt(ts)
        if (.not.associated(t0%vegn)) cycle ! do nothing for non-vegetated cycles
        ! calculate death rate for each of the cohorts
        allocate(ndead(t0%vegn%n_cohorts))
        do k = 1,t0%vegn%n_cohorts
           call cohort_nat_mortality_ppa(t0%vegn%cohorts(k), seconds_per_year, ndead(k))
        enddo
        ! given death numbers for each cohort, update tile
        call tile_nat_mortality_ppa(t0,ndead,t1)
        deallocate(ndead)
        if (associated(t1)) call insert(t1,disturbed_tiles)
     enddo
     ! merge disturbed tiles back into the original list
     te = tail_elmt(disturbed_tiles)
     do 
        ts = first_elmt(disturbed_tiles)
        if (ts == te) exit ! reached the end of the list
        t0=>current_tile(ts)
        call remove(ts) ; call merge_land_tile_into_list(t0,lnd%tile_map(i,j))
     enddo
     ! at this point the list of disturbed tiles must be empty 
     if (.not.empty(disturbed_tiles)) call error_mesg('vegn_nat_mortality_ppa', &
        'list of disturbed tiles is not empty at the end of the processing', FATAL)
  enddo
  enddo

  call land_tile_list_end(disturbed_tiles)

end subroutine vegn_nat_mortality_ppa


! ============================================================================
! for a given cohort and time interval, calculates how many individuals die 
! naturally during this interval
subroutine cohort_nat_mortality_ppa(cc, deltat, ndead)
  type(vegn_cohort_type), intent(in)  :: cc     ! cohort
  real,                   intent(in)  :: deltat ! time interval, s
  real,                   intent(out) :: ndead  ! number of individuals to die, indiv/m2

  real :: deathrate ! mortality rate, 1/year
  logical :: do_U_shaped_mortality = .FALSE.
  real, parameter :: DBHtp = 1.0 ! for U-shaped mortality

  associate ( sp => spdata(cc%species) )
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

  ndead = cc%nindivs * (1.0-exp(-deathrate*deltat/seconds_per_year)) ! individuals / m2
  end associate
end subroutine cohort_nat_mortality_ppa


! ============================================================================
subroutine tile_nat_mortality_ppa(t0,ndead,t1)
  type(land_tile_type), intent(inout) :: t0 ! original tile
  real,                 intent(in)    :: ndead(:) ! number of dying individuals, indiv/(m2 of t0)
  type(land_tile_type), pointer       :: t1 ! optionally created disturbed tile

  integer :: n_layers ! number of canopy layers in the tile
  real :: dying_crownwarea ! area of the crowns to die, m2/m2
  real :: f0 ! fraction of original tile in the grid cell
  real, dimension(N_C_TYPES) :: &
     leaf_litt0(N_C_TYPES),  leaf_litt1(N_C_TYPES), & ! accumulated leaf litter, kg C/m2
     wood_litt0(N_C_TYPES),  wood_litt1(N_C_TYPES)    ! accumulated wood litter, kg C/m2
  real :: root_litt0(num_l, N_C_TYPES) ! accumulated root litter per soil layer, kgC/m2
  real :: root_litt1(num_l, N_C_TYPES) ! accumulated root litter per soil layer, kgC/m2
  integer :: i
  ! variables for conservation checks:
  real :: lmass0, fmass0, heat0, cmass0
  real :: lmass1, fmass1, heat1, cmass1
  real :: lmass2, fmass2
  character(64) :: tag

  t1=>NULL()
  if(.not.associated(t0%vegn)) &
      call error_mesg('tile_nat_mortality_ppa','attempt to kill plants in non-vegetated tile', FATAL)
  if(size(ndead)/=t0%vegn%n_cohorts) &
      call error_mesg('tile_nat_mortality_ppa','size of input argument does not match n_cohorts', FATAL)

  f0 = t0%frac ! save original tile fraction for future use

  if (do_check_conservation) then
     ! + conservation check, part 1: calculate the pre-transition totals
     call get_tile_water(t0,lmass0,fmass0)
     heat0  = land_tile_heat  (t0)
     cmass0 = land_tile_carbon(t0)
     ! - end of conservation check, part 1
  endif

  ! count layers and determine if any vegetation in the canopy layer dies
  n_layers         = 0
  dying_crownwarea = 0
  do i = 1,t0%vegn%n_cohorts
     n_layers = max(n_layers,t0%vegn%cohorts(i)%layer)
     if (t0%vegn%cohorts(i)%layer==1) then
        dying_crownwarea = dying_crownwarea + ndead(i)*t0%vegn%cohorts(i)%crownarea
     endif
  enddo
  
  ! zero out litter terms, for accumulation across cohorts
  leaf_litt0 = 0.0; wood_litt0 = 0.0; root_litt0 = 0.0
  leaf_litt1 = 0.0; wood_litt1 = 0.0; root_litt1 = 0.0

  if (n_layers>1.and.dying_crownwarea>0.and.nat_mortality_splits_tiles) then
     ! write(*,*) 'splitting tiles'
     ! split the disturbed fraction of vegetation as a new tile. The new tile
     ! will contain all crown area trees that are dying, while the
     ! original tile has all the crown trees that survive.
     ! Note that in both ne and original tiles nindivs of crown trees are different
     ! than in the original tile before the 
     t1 => new_land_tile(t0)
     t1%frac = t0%frac*dying_crownwarea
     t0%frac = t0%frac - t1%frac ! and update it to conserve area
     
     do i = 1,t0%vegn%n_cohorts
        associate(cc0=>t0%vegn%cohorts(i), cc1=>t1%vegn%cohorts(i))
        if (cc0%layer==1) then
           cc1%nindivs = ndead(i)*f0/t1%frac
           cc0%nindivs = (cc0%nindivs-ndead(i))*f0/t0%frac
           ! kill all canopy plants in disturbed cohort
           call kill_plants_ppa(cc1,t1%vegn,t1%soil,cc1%nindivs,0.0,leaf_litt1,wood_litt1,root_litt1)
        else
           ! kill understory plants in both tiles
           call kill_plants_ppa(cc1,t1%vegn,t1%soil,ndead(i),0.0,leaf_litt1,wood_litt1,root_litt1)
           call kill_plants_ppa(cc0,t0%vegn,t0%soil,ndead(i),0.0,leaf_litt0,wood_litt0,root_litt0)
        endif
        end associate
     enddo
  else
     ! just kill the trees in t0
     do i = 1,t0%vegn%n_cohorts
        ! write(*,*) 'not splitting tiles',i,t0%vegn%cohorts(i)%nindivs,ndead(i),t0%vegn%cohorts(i)%carbon_gain,t0%vegn%cohorts(i)%bwood_gain
        call kill_plants_ppa(t0%vegn%cohorts(i),t0%vegn,t0%soil,ndead(i),0.0,&
                             leaf_litt0,wood_litt0,root_litt0)
     enddo
  endif

  call add_soil_carbon(t0%soil, leaf_litt0, wood_litt0, root_litt0)
  if (associated(t1)) &
     call add_soil_carbon(t1%soil, leaf_litt1, wood_litt1, root_litt1)

  if (do_check_conservation) then
     ! + conservation check, part 2: calculate totals in final state, and compare 
     ! with previous totals
     tag = 'tile_nat_mortality_ppa'
     call get_tile_water(t0,lmass1,fmass1); lmass1 = lmass1*t0%frac; fmass1 = fmass1*t0%frac
     heat1  = land_tile_heat  (t0)*t0%frac
     cmass1 = land_tile_carbon(t0)*t0%frac
     if (associated(t1)) then
        call get_tile_water(t1,lmass2,fmass2); 
        lmass1 = lmass1+lmass2*t1%frac; fmass1 = fmass1+fmass2*t1%frac
        heat1  = heat1+land_tile_heat(t1)*t1%frac
        cmass1 = cmass1+land_tile_carbon(t1)*t1%frac
     endif
     call check_conservation (tag,'liquid water', lmass0, lmass1/f0, water_cons_tol)
     call check_conservation (tag,'frozen water', fmass0, fmass1/f0, water_cons_tol)
     call check_conservation (tag,'carbon'      , cmass0, cmass1/f0, carbon_cons_tol)
!     call check_conservation (tag,'heat content', heat0 , heat1/f0 , 1e-16)
     ! - end of conservation check, part 2
  endif

  call relayer_cohorts(t0%vegn)
  if (associated(t1)) call relayer_cohorts(t1%vegn)

end subroutine tile_nat_mortality_ppa


! ==============================================================================
! given a cohort, number of individuals to kill, and fraction of smoke, kills
! specified fraction of individuals, putting the biomass from the rest into
! soil and smoke pools. Updates carbon pools and water pools (if there is water
! on killed trees) in soil and vegn tiles.

! NOTE that in contrast to LM3 mortality calculation, living biomass is
! included in losses; does it mean that we double-counting losses in 
! maintenance and here?

! TODO: ask Ensheng and Elena about possible double-counting of live biomass losses.

! we are burning all of biomass (like we did in LM3), including the below-ground 
! part. Perhaps we should leave the below-ground part alone, and add it to soil 
! carbon pools?

! TODO: ask Elena about burning of below-ground biomass

subroutine kill_plants_ppa(cc, vegn, soil, ndead, fsmoke, leaf_litt, wood_litt, root_litt)
  type(vegn_cohort_type), intent(inout) :: cc
  type(vegn_tile_type),   intent(inout) :: vegn
  type(soil_tile_type),   intent(inout) :: soil
  real,                   intent(in)    :: ndead ! number of individuals to kill, indiv./m2 
  real,                   intent(in)    :: fsmoke ! fraction of biomass lost to fire, unitless 
  real, intent(inout) :: leaf_litt(N_C_TYPES) ! accumulated leaf litter, kg C/m2
  real, intent(inout) :: wood_litt(N_C_TYPES) ! accumulated wood litter, kg C/m2
  real, intent(inout) :: root_litt(num_l, N_C_TYPES) ! accumulated root litter per soil layer, kgC/m2

  ! ---- local vars
  real :: lost_wood, lost_alive, burned_wood, burned_alive
  real :: profile(num_l) ! storage for vertical profile of exudates and root litter
  integer :: l
  
  call check_var_range(ndead,  0.0, cc%nindivs, 'kill_plants_ppa', 'ndead',  FATAL)
  call check_var_range(fsmoke, 0.0, 1.0,        'kill_plants_ppa', 'fsmoke', FATAL)
  
  ! reduce the number of individuals
  cc%nindivs = cc%nindivs-ndead
  
  ! water from dead trees goes to intermediate buffers, to be added to the
  ! precipitation reaching ground on the next physical time step
  vegn%drop_wl = vegn%drop_wl + cc%wl*ndead
  vegn%drop_ws = vegn%drop_ws + cc%ws*ndead
  vegn%drop_hl = vegn%drop_hl + clw*cc%wl*ndead*(cc%Tv-tfreeze)
  vegn%drop_hs = vegn%drop_hs + csw*cc%ws*ndead*(cc%Tv-tfreeze)

  ! calculate total carbon losses, kgC/m2
  lost_wood  = ndead * (cc%bwood+cc%bsw+cc%bwood_gain)
  lost_alive = ndead * (cc%bl+cc%br+cc%blv+cc%bseed+cc%nsc+cc%carbon_gain)
  ! loss to fire
  burned_wood  = fsmoke*lost_wood
  burned_alive = fsmoke*lost_alive
  ! loss to the soil pools
  lost_wood  = lost_wood  - burned_wood
  lost_alive = lost_alive - burned_alive

  ! add fire carbon losses to smoke pool
  vegn%csmoke_pool = vegn%csmoke_pool + burned_wood + burned_alive
    
  ! add remaining lost C to soil carbon pools
  leaf_litt(:) = leaf_litt(:) + [fsc_liv,  1-fsc_liv,  0.0]*(cc%bl+cc%bseed+cc%carbon_gain)*(1-fsmoke)*ndead
  wood_litt(:) = wood_litt(:) + [fsc_wood, 1-fsc_wood, 0.0]*(cc%bwood+cc%bsw+cc%bwood_gain)*(1-fsmoke)*agf_bs*ndead
  wood_litt(C_CEL) = wood_litt(C_CEL)+cc%nsc*(1-fsmoke)*agf_bs*ndead
  call cohort_root_litter_profile(cc, dz, profile)
  do l = 1, num_l
     root_litt(l,:) = root_litt(l,:) + profile(l)*ndead*(1-fsmoke)*(/ &
          fsc_froot    *cc%br + fsc_wood    *(cc%bsw+cc%bwood+cc%bwood_gain)*(1-agf_bs) + cc%nsc*(1-agf_bs), &
          (1-fsc_froot)*cc%br + (1-fsc_wood)*(cc%bsw+cc%bwood+cc%bwood_gain)*(1-agf_bs), &
          0.0/)
  enddo
       
  ! for budget tracking - temporary
  vegn%veg_out = vegn%veg_out + lost_alive + lost_wood + burned_alive + burned_wood
end subroutine kill_plants_ppa

end module vegn_disturbance_mod
