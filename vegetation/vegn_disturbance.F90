!***********************************************************************
!*                   GNU Lesser General Public License
!*
!* This file is part of the GFDL Land Model 4 (LM4).
!*
!* LM4 is free software: you can redistribute it and/or modify it under
!* the terms of the GNU Lesser General Public License as published by
!* the Free Software Foundation, either version 3 of the License, or (at
!* your option) any later version.
!*
!* LM4 is distributed in the hope that it will be useful, but WITHOUT
!* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
!* FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
!* for more details.
!*
!* You should have received a copy of the GNU Lesser General Public
!* License along with LM4.  If not, see <http://www.gnu.org/licenses/>.
!***********************************************************************
! ============================================================================
! vegetation disturbances
! ============================================================================
module vegn_disturbance_mod

#include "../shared/debug.inc"

use fms_mod,         only : error_mesg, FATAL
use time_manager_mod,only : get_date, operator(-)
use constants_mod,   only : Tfreeze
use land_constants_mod, only : seconds_per_year
use land_debug_mod,  only : is_watch_point, is_watch_cell, set_current_point, &
     check_conservation, do_check_conservation, water_cons_tol, carbon_cons_tol, &
     heat_cons_tol, nitrogen_cons_tol, check_var_range, land_error_message
use vegn_data_mod,   only : do_ppa, nat_mortality_splits_tiles, spdata, agf_bs, &
     FORM_GRASS, FORM_WOODY, LEAF_OFF, DBH_mort, A_mort, B_mort, cold_mort, treeline_mort, &
     treeline_base_T, treeline_thresh_T, treeline_season_length, &
     COLD_INTOLERANT, WARM_INTOLERANT
use land_tile_diag_mod, only : set_default_diag_filter, register_tiled_diag_field, send_tile_data
use vegn_tile_mod,   only : vegn_tile_type, vegn_relayer_cohorts_ppa, vegn_tile_bwood, &
     vegn_mergecohorts_ppa
use snow_tile_mod,   only : snow_active
use soil_tile_mod,   only : soil_tile_type, num_l, dz
use soil_util_mod,   only : add_soil_carbon
use land_tile_mod,   only : land_tile_map, land_tile_type, land_tile_enum_type, &
     land_tile_list_type, land_tile_list_init, land_tile_list_end, &
     empty, first_elmt, tail_elmt, merge_land_tile_into_list, loop_over_tiles, &
     current_tile, operator(==), operator(/=), remove, insert, new_land_tile, &
     land_tile_heat, land_tile_carbon, land_tile_nitrogen, get_tile_water, nitems
use land_data_mod,   only : lnd, log_version
use soil_carbon_mod, only : add_litter, soil_carbon_option, &
     SOILC_CENTURY, SOILC_CENTURY_BY_LAYER, SOILC_CORPSE, N_C_TYPES, C_FAST
use vegn_cohort_mod, only : vegn_cohort_type, update_biomass_pools, &
     cohort_root_litter_profile, cohort_root_exudate_profile
use vegn_util_mod, only : kill_plants_ppa

implicit none
private

! ==== public interfaces =====================================================
public :: vegn_disturbance_init
public :: vegn_nat_mortality_lm3
public :: vegn_nat_mortality_ppa
public :: vegn_disturbance
public :: update_fuel
! =====end of public interfaces ==============================================

! ==== module constants ======================================================
character(len=*), parameter :: module_name = 'vegn_disturbance_mod'
#include "../shared/version_variable.inc"

! IDs of diagnostic variables
integer :: id_treeline_T, id_treeline_season

contains ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

! =============================================================================
subroutine vegn_disturbance_init(id_ug)
  integer, intent(in) :: id_ug   !<Unstructured axis id.

  call log_version(version, module_name, &
  __FILE__)

  ! set the default sub-sampling filter for the fields below
  call set_default_diag_filter('soil')

  id_treeline_T = register_tiled_diag_field ('vegn', 'treeline_T',(/id_ug/), &
       lnd%time, 'average temperature of growing season for treeline calculations', 'degK', &
       missing_value=-999.0)
  id_treeline_season = register_tiled_diag_field ('vegn', 'treeline_season',(/id_ug/), &
       lnd%time, 'length of growing season for treeline calculations', 'degK', &
       missing_value=-999.0)
end subroutine vegn_disturbance_init

subroutine vegn_disturbance(vegn, soil, dt)
  type(vegn_tile_type), intent(inout) :: vegn ! vegetation data
  type(soil_tile_type), intent(inout) :: soil ! soil data
  real, intent(in) :: dt ! time since last disturbance calculations, s

  real, parameter :: BMIN = 1e-10; ! should be the same as in growth function
  ! ---- local vars
  real :: precip;
  real :: delta_C, delta_N
  real :: fraction_lost;
  real :: drought_month;
  real :: deltat
  integer :: i,l
  real :: leaf_litt_C(N_C_TYPES),leaf_litt_N(N_C_TYPES) ! fine surface litter per tile, kgC/m2
  real :: wood_litt_C(N_C_TYPES),wood_litt_N(N_C_TYPES) ! coarse surface litter per tile, kgC/m2
  real :: root_litt_C(num_l, N_C_TYPES),root_litt_N(num_l, N_C_TYPES) ! root litter per soil layer, kgC/m2
  real :: profile(num_l) ! storage for vertical profile of exudates and root litter

  deltat = dt/seconds_per_year ! convert time interval to years

  !  Disturbance Rates
  precip=vegn%p_ann*86400*365;
  drought_month = vegn%lambda;
  vegn%disturbance_rate(0) = 0.0;
  vegn%disturbance_rate(1) = 0.0;

  call calculate_patch_disturbance_rates(vegn)

  leaf_litt_C = 0 ; wood_litt_C = 0; root_litt_C = 0
  leaf_litt_N = 0 ; wood_litt_N = 0; root_litt_N = 0
  do i = 1,vegn%n_cohorts
     associate (cc => vegn%cohorts(i), &
                sp => spdata(vegn%cohorts(i)%species)) ! F2003

     fraction_lost = 1.0-exp(-vegn%disturbance_rate(1)*deltat);
     if (do_ppa) then
        call kill_plants_ppa(cc, vegn, cc%nindivs*fraction_lost, sp%smoke_fraction, &
                             leaf_litt_C, wood_litt_C, root_litt_C, &
                             leaf_litt_N, wood_litt_N, root_litt_N)
     else ! original LM3 treatment
        ! "dead" biomass : wood + sapwood
        delta_C = fraction_lost*cc%nindivs*(cc%bwood+cc%bsw)
        delta_N = fraction_lost*cc%nindivs*(cc%wood_N+cc%sapwood_N)

        ! accumulate liter and soil carbon inputs across all cohorts
        wood_litt_C(:) = wood_litt_C(:) + agf_bs*(1-sp%smoke_fraction)* &
                        delta_C*[sp%fsc_wood, 1-sp%fsc_wood, 0.0]
        wood_litt_N(:) = wood_litt_N(:) + agf_bs*(1-sp%smoke_fraction)* &
                        delta_N*[sp%fsc_wood, 1-sp%fsc_wood, 0.0]
        call cohort_root_litter_profile(cc, dz, profile)
        do l = 1, num_l
           root_litt_C(l,:) = root_litt_C(l,:) + profile(l)*delta_C*(1-agf_bs)*(1-sp%smoke_fraction) * &
                       [sp%fsc_wood, 1-sp%fsc_wood, 0.0]
           root_litt_N(l,:) = root_litt_N(l,:) + profile(l)*delta_N*(1-agf_bs)*(1-sp%smoke_fraction) * &
                       [sp%fsc_wood, 1-sp%fsc_wood, 0.0]
        enddo

        cc%bwood     = cc%bwood     * (1-fraction_lost)
        cc%bsw       = cc%bsw       * (1-fraction_lost)
        cc%wood_N    = cc%wood_N    * (1-fraction_lost)
        cc%sapwood_N = cc%sapwood_N * (1-fraction_lost)

        vegn%Nsmoke_pool = vegn%Nsmoke_pool + sp%smoke_fraction*delta_N
        vegn%csmoke_pool = vegn%csmoke_pool + sp%smoke_fraction*delta_C
        vegn%veg_out = vegn%veg_out+delta_C

        !"alive" biomass: leaves, roots, and virtual pool
        delta_C = (cc%bl+cc%blv+cc%br)*fraction_lost*cc%nindivs;
        delta_N = (cc%leaf_N+cc%stored_N+cc%root_N)*fraction_lost*cc%nindivs
        leaf_litt_C(:) = leaf_litt_C(:) + [sp%fsc_liv, 1-sp%fsc_liv, 0.0]*(cc%bl+cc%blv)*fraction_lost*(1-sp%smoke_fraction)*cc%nindivs
        leaf_litt_N(:) = leaf_litt_N(:) + [sp%fsc_liv, 1-sp%fsc_liv, 0.0]*cc%leaf_N*fraction_lost*(1-sp%smoke_fraction)*cc%nindivs
        do l = 1, num_l
           root_litt_C(l,:) = root_litt_C(l,:) + [sp%fsc_froot, 1-sp%fsc_froot, 0.0] * &
                            profile(l)*cc%br*fraction_lost*(1-sp%smoke_fraction)*cc%nindivs
           root_litt_N(l,:) = root_litt_N(l,:) + [sp%fsc_froot, 1-sp%fsc_froot, 0.0] * &
                            profile(l)*(cc%root_N+cc%stored_N)*fraction_lost*(1-sp%smoke_fraction)*cc%nindivs
        enddo

        cc%bl     = cc%bl      * (1-fraction_lost)
        cc%blv    = cc%blv     * (1-fraction_lost)
        cc%br     = cc%br      * (1-fraction_lost)
        cc%leaf_N = cc%leaf_N  * (1-fraction_lost)
        cc%root_N = cc%root_N  * (1-fraction_lost)
        cc%stored_N = cc%stored_N * (1-fraction_lost)

        vegn%Nsmoke_pool = vegn%Nsmoke_pool + sp%smoke_fraction*delta_N
        vegn%csmoke_pool = vegn%csmoke_pool + sp%smoke_fraction*delta_C
        vegn%veg_out = vegn%veg_out+delta_C

        !"living" biomass:leaves, roots and sapwood
        delta_C = cc%bliving*fraction_lost;
        cc%bliving = cc%bliving - delta_C

        if(cc%bliving < BMIN) then
           ! remove vegetaion competely
           ! accumulate liter and soil carbon inputs across all cohorts
           leaf_litt_C(:) = leaf_litt_C(:) + [sp%fsc_liv,  1-sp%fsc_liv,  0.0]*(cc%bl+cc%blv)
           leaf_litt_N(:) = leaf_litt_N(:) + [sp%fsc_liv,  1-sp%fsc_liv,  0.0]*cc%leaf_N
           wood_litt_C(:) = wood_litt_C(:) + [sp%fsc_wood, 1-sp%fsc_wood, 0.0]*agf_bs*(cc%bwood+cc%bsw)
           wood_litt_N(:) = wood_litt_N(:) + [sp%fsc_wood, 1-sp%fsc_wood, 0.0]*agf_bs*(cc%wood_N+cc%sapwood_N)
           do l = 1, num_l
              root_litt_C(l,:) = root_litt_C(l,:) + profile(l)*[ &
                   sp%fsc_froot    *cc%br + sp%fsc_wood    *(cc%bwood+cc%bsw)*(1-agf_bs), & ! fast
                   (1-sp%fsc_froot)*cc%br + (1-sp%fsc_wood)*(cc%bwood+cc%bsw)*(1-agf_bs), & ! slow
                   0.0] ! microbes
              root_litt_N(l,:) = root_litt_N(l,:) + profile(l)*[ &
                   cc%root_N*   sp%fsc_froot  + (cc%wood_N+cc%sapwood_N)*   sp%fsc_wood *(1-agf_bs), & ! fast
                   cc%root_N*(1-sp%fsc_froot) + (cc%wood_N+cc%sapwood_N)*(1-sp%fsc_wood)*(1-agf_bs), & ! slow
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

  call add_soil_carbon(soil, vegn, leaf_litt_C, wood_litt_C, root_litt_C, &
                                   leaf_litt_N, wood_litt_N, root_litt_N  )

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
  real :: delta_C, delta_N
  real :: fraction_lost
  real :: bdead, balive ! combined biomass pools
  integer :: i,l
  real :: wood_litt_C(n_c_types) ! coarse surface litter per tile, kgC/m2
  real :: wood_litt_N(n_c_types) ! coarse surface litter nitrogen per tile, kgN/m2
  real :: leaf_litt_C(n_c_types) ! surface leaf litter per tile, kgC/m2
  real :: leaf_litt_N(n_c_types) ! surface leaf litter nitrogen per tile, kgN/m2
  real :: root_litt_C(num_l, n_c_types) ! root litter per soil layer, kgC/m2
  real :: root_litt_N(num_l, n_c_types) ! root litter nitrogen per soil layer, kgN/m2
  real :: profile(num_l) ! storage for vertical profile of exudates and root litter

  vegn%disturbance_rate(0) = 0.0;
  wood_litt_C = 0; root_litt_C = 0; leaf_litt_C = 0
  wood_litt_N = 0; root_litt_N = 0; leaf_litt_N = 0
  do i = 1,vegn%n_cohorts
     associate(cc => vegn%cohorts(i), sp => spdata(vegn%cohorts(i)%species))
     ! Treat treefall disturbance implicitly, i.e. not creating a new tile.
     ! note that this disturbance rate calculation only works for the one cohort per
     ! tile case -- in case of multiple cohort disturbance rate perhaps needs to be
     ! accumulated (or averaged? or something else?) over the cohorts.
     vegn%disturbance_rate(0) = spdata(cc%species)%treefall_disturbance_rate;

     ! calculate combined biomass pools
     balive = cc%bl + cc%blv + cc%br;
     bdead  = cc%bsw + cc%bwood;
     ! ens need a daily PATCH_FREQ here, for now it is set to 48
     fraction_lost = 1.0-exp(-vegn%disturbance_rate(0)*deltat/seconds_per_year)

     ! "dead" biomass : wood + sapwood
     delta_C = bdead*fraction_lost
     delta_N = (cc%wood_N+cc%sapwood_N)*fraction_lost

     wood_litt_C(:) = wood_litt_C(:) + [sp%fsc_wood, 1-sp%fsc_wood, 0.0]*cc%bwood*fraction_lost*agf_bs + &
        [sp%fsc_liv, 1-sp%fsc_liv, 0.0]*cc%bsw*fraction_lost*agf_bs
     wood_litt_N(:) = wood_litt_N(:) + [sp%fsc_wood, 1-sp%fsc_wood, 0.0]*cc%wood_N*fraction_lost*agf_bs + &
        [sp%fsc_liv, 1-sp%fsc_liv, 0.0]*cc%sapwood_N*fraction_lost*agf_bs
     if (sp%mortality_kills_balive) then
        leaf_litt_C(:) = leaf_litt_C(:) + [sp%fsc_liv, 1-sp%fsc_liv, 0.0]*fraction_lost*(cc%bl+cc%blv)
        leaf_litt_N(:) = leaf_litt_N(:) + [sp%fsc_liv, 1-sp%fsc_liv, 0.0]*fraction_lost*(cc%leaf_N+cc%stored_N)
     endif
     call cohort_root_litter_profile(cc, dz, profile)
     do l = 1, num_l
        root_litt_C(l,:) = root_litt_C(l,:) + profile(l)*(1-agf_bs)*fraction_lost* &
            (cc%bwood*[sp%fsc_wood, 1-sp%fsc_wood, 0.0]+cc%bsw*[sp%fsc_liv, 1-sp%fsc_liv, 0.0])
        root_litt_N(l,:) = root_litt_N(l,:) + profile(l)*(1-agf_bs)*fraction_lost* &
            (cc%wood_N*[sp%fsc_wood, 1-sp%fsc_wood, 0.0]+cc%sapwood_N*[sp%fsc_liv, 1-sp%fsc_liv, 0.0])
        if (sp%mortality_kills_balive) then
           root_litt_C(l,:) = root_litt_C(l,:) + profile(l)*[sp%fsc_froot, 1-sp%fsc_froot, 0.0]*fraction_lost*cc%br
           root_litt_N(l,:) = root_litt_N(l,:) + profile(l)*[sp%fsc_froot, 1-sp%fsc_froot, 0.0]*fraction_lost*cc%root_N
        endif
     enddo

     cc%bwood     = cc%bwood     * (1-fraction_lost)
     cc%bsw       = cc%bsw       * (1-fraction_lost)
     cc%wood_N    = cc%wood_N    * (1-fraction_lost)
     cc%sapwood_N = cc%sapwood_N * (1-fraction_lost)
     if (sp%mortality_kills_balive) then
        delta_C = delta_C+(cc%br+cc%bl+cc%blv)*(1-fraction_lost) ! For correct vegn%veg_out bookkeeping
        delta_N = delta_N+(cc%leaf_N+cc%stored_N+cc%root_N)*(1-fraction_lost)
        cc%br       = cc%br       * (1-fraction_lost)
        cc%bl       = cc%bl       * (1-fraction_lost)
        cc%blv      = cc%blv      * (1-fraction_lost)
        cc%leaf_N   = cc%leaf_N   * (1-fraction_lost)
        cc%root_N   = cc%root_N   * (1-fraction_lost)
        cc%stored_N = cc%stored_N * (1-fraction_lost)
     endif

     vegn%veg_out = vegn%veg_out + delta_C

     ! note that fast "living" pools are not included into mortality because their
     ! turnover is calculated separately

     cc%bliving = cc%bsw + cc%bl + cc%br + cc%blv;
     call update_biomass_pools(cc);
     end associate
  enddo
  ! add litter accumulated over the cohorts
  call add_soil_carbon(soil, vegn, wood_litter_C=wood_litt_C, leaf_litter_C=leaf_litt_C, root_litter_C=root_litt_C, &
                                   wood_litter_N=wood_litt_N, leaf_litter_N=leaf_litt_N, root_litter_N=root_litt_N  )
end subroutine vegn_nat_mortality_lm3


! ============================================================================
! goes through all tiles in the domain, and applies natural mortality to
! each of them, possibly creating new tiles within affected grid cells.
subroutine vegn_nat_mortality_ppa ( )

  integer :: second, minute, hour, day0, day1, month0, month1, year0, year1
  type(land_tile_type), pointer :: t0,t1
  type(land_tile_enum_type) :: ts, te
  type(land_tile_list_type) :: disturbed_tiles
  integer :: k,l
  real, allocatable :: ndead(:)
  ! variable used for conservation check:
  type(land_tile_type), pointer :: ptr
  real :: lmass0, fmass0, cmass0, heat0 ! pre-transition values
  real :: lmass1, fmass1, cmass1, heat1 ! post-transition values
  real :: lm, fm
  real :: treeline_T ! average temperature for tree line calculations
  character(*),parameter :: tag = 'vegn_nat_mortality_ppa'

  ! get components of calendar dates for this and previous time step
  call get_date(lnd%time-lnd%dt_slow, year1,month1,day1,hour,minute,second)
  call get_date(lnd%time,             year0,month0,day0,hour,minute,second)

  ! update accumulated values for tree line calculations
  if (day0/=day1) then
     ts = first_elmt(land_tile_map)
     do while (loop_over_tiles(ts,t0))
        if (.not.associated(t0%vegn)) cycle ! do nothing for non-vegetated tiles
        if (t0%vegn%tc_daily > treeline_base_T.and..not.snow_active(t0%snow)) then
           ! accumulate average T over growing season, for treeline/mortality calculations
           t0%vegn%treeline_T_accum = t0%vegn%treeline_T_accum + t0%vegn%tc_daily
           t0%vegn%treeline_N_accum = t0%vegn%treeline_N_accum + 1
        endif
        call send_tile_data(id_treeline_T, t0%vegn%treeline_T_accum/max(1.0,t0%vegn%treeline_N_accum), t0%diag)
        call send_tile_data(id_treeline_season, t0%vegn%treeline_N_accum, t0%diag)
     enddo
  endif

  if (.not.(do_ppa.and.year1 /= year0)) return  ! do nothing

  call land_tile_list_init(disturbed_tiles)

  do l = lnd%ls,lnd%le
     if(empty(land_tile_map(l))) cycle ! skip cells where there is no land
     ! set current point for debugging
     call set_current_point(l,1)

     if (do_check_conservation) then
        ! conservation check code, part 1: calculate the pre-transition grid
        ! cell totals
        lmass0 = 0 ; fmass0 = 0 ; cmass0 = 0 ; heat0 = 0
        ts = first_elmt(land_tile_map(l))
        do while (loop_over_tiles(ts,ptr))
           call get_tile_water(ptr,lm,fm)
           lmass0 = lmass0 + lm*ptr%frac ; fmass0 = fmass0 + fm*ptr%frac
           heat0  = heat0  + land_tile_heat  (ptr)*ptr%frac
           cmass0 = cmass0 + land_tile_carbon(ptr)*ptr%frac
        enddo
     endif

     ts = first_elmt(land_tile_map(l))
     do while (loop_over_tiles(ts,t0))
        if (.not.associated(t0%vegn)) cycle ! do nothing for non-vegetated cycles
        ! calculate average T in growing season, for tree line calculations
        if (t0%vegn%treeline_N_accum>0) then
            treeline_T = t0%vegn%treeline_T_accum/t0%vegn%treeline_N_accum
        else
            treeline_T = Tfreeze
        endif
        ! calculate death rate for each of the cohorts
        allocate(ndead(t0%vegn%n_cohorts))
        do k = 1,t0%vegn%n_cohorts
           call cohort_nat_mortality_ppa(t0%vegn%cohorts(k), seconds_per_year, t0%vegn%t_cold, &
                treeline_T, t0%vegn%treeline_N_accum, ndead(k))
        enddo
        if (is_watch_point()) then
           write(*,*)'##### vegn_nat_mortalty_ppa #####'
           call dpri('frac',t0%frac)
           call dpri('coldest_month_T',t0%vegn%t_cold)
           call dpri('treeline_T',treeline_T)
           call dpri('treeline_season',t0%vegn%treeline_N_accum)
           write(*,*)
           do k = 1,t0%vegn%n_cohorts
              associate(cc=>t0%vegn%cohorts(k),sp=>spdata(t0%vegn%cohorts(k)%species))
              write(*,'(i2.2," L ",i2.2,2x,a12)',advance='NO') k, cc%layer, spdata(cc%species)%name
              write(*,'(99(2x,a,g11.4))') &
                  'nindivs',cc%nindivs, &
                  'ndead',ndead(k), &
                  'DBH',cc%DBH, &
                  'H',cc%height, &
                  'Carea',cc%crownarea, &
                  'frac',cc%layerfrac, &
                  'LAI',cc%lai
              end associate
           enddo
        endif
        ! reset treeline_T_accum and treeline_N_accum for the next year
        t0%vegn%treeline_T_accum = 0.0; t0%vegn%treeline_N_accum = 0
        ! given death numbers for each cohort, update tile
        call tile_nat_mortality_ppa(t0,ndead,t1)
        deallocate(ndead)
        if (associated(t1)) call insert(t1,disturbed_tiles)
     enddo
     if (is_watch_cell()) then
        write(*,*) '#### vegn_nat_mortality_ppa ####'
        write(*,*) 'N tiles before merge = ', nitems(land_tile_map(l))
        write(*,*) 'N of disturbed tiles = ', nitems(disturbed_tiles)
     endif
     ! merge disturbed tiles back into the original list
     te = tail_elmt(disturbed_tiles)
     do
        ts = first_elmt(disturbed_tiles)
        if (ts == te) exit ! reached the end of the list
        t0=>current_tile(ts)
        call remove(ts) ; call merge_land_tile_into_list(t0,land_tile_map(l))
     enddo
     if (is_watch_cell()) then
        write(*,*) 'N tiles after merge = ', nitems(land_tile_map(l))
     endif
     ! at this point the list of disturbed tiles must be empty
     if (.not.empty(disturbed_tiles)) call error_mesg('vegn_nat_mortality_ppa', &
        'list of disturbed tiles is not empty at the end of the processing', FATAL)

     if (do_check_conservation) then
        ! conservation check part 2: calculate grid cell totals in final state, and
        ! compare them with pre-transition totals
        lmass1 = 0 ; fmass1 = 0 ; cmass1 = 0 ; heat1 = 0
        ts = first_elmt(land_tile_map(l))
        do while (loop_over_tiles(ts,ptr))
           call get_tile_water(ptr,lm,fm)
           lmass1 = lmass1 + lm*ptr%frac ; fmass1 = fmass1 + fm*ptr%frac
           heat1  = heat1  + land_tile_heat  (ptr)*ptr%frac
           cmass1 = cmass1 + land_tile_carbon(ptr)*ptr%frac
        enddo
        call check_conservation (tag,'liquid water', lmass0, lmass1, water_cons_tol)
        call check_conservation (tag,'frozen water', fmass0, fmass1, water_cons_tol)
        call check_conservation (tag,'carbon'      , cmass0, cmass1, carbon_cons_tol)
        call check_conservation (tag,'heat content', heat0 , heat1 , heat_cons_tol)
     endif
  enddo

  call land_tile_list_end(disturbed_tiles)

end subroutine vegn_nat_mortality_ppa


! ============================================================================
! for a given cohort and time interval, calculates how many individuals die
! naturally during this interval
subroutine cohort_nat_mortality_ppa(cc, deltat, coldest_month_T, treeline_T, treeline_season, ndead)
  type(vegn_cohort_type), intent(in)  :: cc     ! cohort
  real,                   intent(in)  :: deltat ! time interval, s
  real,                   intent(in)  :: coldest_month_T ! temperature of coldest month, degK
  real,                   intent(in)  :: treeline_T ! temperature for treeline calculations, degK
  real,                   intent(in)  :: treeline_season ! growing season length for treeline calculations, days
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

  ! cold mortality
  select case (sp%T_tolerance_type)
  case (COLD_INTOLERANT)
     if (coldest_month_T < sp%Tmin_mort) &
           deathrate = max(deathrate, cold_mort) ! at least 1-exp(-cold_mort) trees die per year
  case (WARM_INTOLERANT)
     if (coldest_month_T > sp%Tmin_mort) &
           deathrate = max(deathrate, cold_mort) ! at least 1-exp(-cold_mort) trees die per year
  end select
  ! tree line
  if (sp%lifeform == FORM_WOODY .and. &
        (treeline_T < treeline_thresh_T.or.treeline_season < treeline_season_length)) &
        deathrate = max(deathrate, treeline_mort) ! at least 1-exp(-treeline_mort) trees die per year

  ndead = cc%nindivs * (1.0-exp(-deathrate*deltat/seconds_per_year)) ! individuals / m2
  end associate
end subroutine cohort_nat_mortality_ppa


! ============================================================================
subroutine tile_nat_mortality_ppa(t0,ndead,t1)
  type(land_tile_type), intent(inout) :: t0 ! original tile
  real,                 intent(in)    :: ndead(:) ! number of dying individuals, indiv/(m2 of t0)
  type(land_tile_type), pointer       :: t1 ! optionally created disturbed tile

  integer :: n_layers ! number of canopy layers in the tile
  real :: dying_crownwarea ! area of the plant crowns to die, m2/m2
  real :: f0 ! fraction of original tile in the grid cell
  real, dimension(N_C_TYPES) :: &
     leaf_litt0_C,  leaf_litt1_C, & ! accumulated leaf litter, kg C/m2
     leaf_litt0_N,  leaf_litt1_N, & ! accumulated leaf litter nitrogen, kg N/m2
     wood_litt0_C,  wood_litt1_C, & ! accumulated wood litter, kg C/m2
     wood_litt0_N,  wood_litt1_N    ! accumulated wood litter nitrogen, kg N/m2
  real, dimension(num_l, N_C_TYPES) :: &
     root_litt0_C, & ! accumulated root litter per soil layer, kg C/m2
     root_litt0_N, & ! accumulated root litter nitrogen per soil layer, kg N/m2
     root_litt1_C, & ! accumulated root litter per soil layer, kg C/m2
     root_litt1_N    ! accumulated root litter nitrogen per soil layer, kg N/m2
  integer :: i
  ! variables for conservation checks:
  real :: lmass0, fmass0, heat0, cmass0, nmass0
  real :: lmass1, fmass1, heat1, cmass1, nmass1
  real :: lmass2, fmass2
  character(*),parameter :: tag = 'tile_nat_mortality_ppa'
  real :: dheat ! heat residual due to cohort merge

  t1=>NULL()
  if(.not.associated(t0%vegn)) &
      call error_mesg('tile_nat_mortality_ppa','attempt to kill plants in non-vegetated tile', FATAL)
  if(size(ndead)/=t0%vegn%n_cohorts) &
      call error_mesg('tile_nat_mortality_ppa','size of input argument does not match n_cohorts', FATAL)

  if (is_watch_point()) then
     write(*,*) '#### tile_mortality_ppa input ####'
     do i = 1,t0%vegn%n_cohorts
        write(*,'(i2.2)', advance='NO') i
        call dpri('Nindivs=',t0%vegn%cohorts(i)%nindivs)
        call dpri('Ndead=',ndead(i))
        call dpri('layer=',t0%vegn%cohorts(i)%layer)
        write(*,*)
     enddo
     write(*,*) '#### end of tile_mortality_ppa input ####'
  endif

  f0 = t0%frac ! save original tile fraction for future use
  ! write (*,*) 'in natural mortality'
  if (do_check_conservation) then
     ! + conservation check, part 1: calculate the pre-transition totals
     call get_tile_water(t0,lmass0,fmass0)
     lmass0 = lmass0*t0%frac; fmass0 = fmass0*t0%frac
     heat0  = land_tile_heat  (t0)*t0%frac
     cmass0 = land_tile_carbon(t0)*t0%frac
     nmass0 = land_tile_nitrogen(t0)*t0%frac
     ! - end of conservation check, part 1
  endif

  ! count layers and determine if any trees in the canopy layer die. Dying grass does not
  ! split tiles, so it is not counted in the crownarea that dies. See NOTE below
  n_layers         = 0
  dying_crownwarea = 0
  do i = 1,t0%vegn%n_cohorts
     associate ( cc => t0%vegn%cohorts(i),   &
                 sp => spdata(t0%vegn%cohorts(i)%species) )
     n_layers = max(n_layers,cc%layer)
     if (cc%layer==1 .and. sp%lifeform/=FORM_GRASS) then
        dying_crownwarea = dying_crownwarea + ndead(i)*cc%crownarea
     endif
     end associate
  enddo

  ! zero out litter terms, for accumulation across cohorts
  leaf_litt0_C = 0.0; wood_litt0_C = 0.0; root_litt0_C = 0.0
  leaf_litt0_N = 0.0; wood_litt0_N = 0.0; root_litt0_N = 0.0
  leaf_litt1_C = 0.0; wood_litt1_C = 0.0; root_litt1_C = 0.0
  leaf_litt1_N = 0.0; wood_litt1_N = 0.0; root_litt1_N = 0.0

  ! NOTE: Perhaps we should not create new tiles at all if trees do not form
  ! closed canopy -- that is, if not all plants in layer 1 are trees. This is
  ! tricky though, in case some grasses are very tall, or if due to numerics a tiny
  ! fraction of grasses happen to be in the canopy layer. For now we use different
  ! approach: we exclude grasses from dying_crownwarea calculations, and in the
  ! tile splitting treat them together with understory plants.
  ! NUMERICS: in rare cases of very low numbers, the product of t0%frac and dying_crownwarea
  ! is zero even if operands are not
  if (n_layers>1.and.dying_crownwarea*t0%frac>0.and.nat_mortality_splits_tiles) then
     ! split the disturbed fraction of vegetation as a new tile. The new tile
     ! will contain all crown area trees that are dying, while the
     ! original tile has all the crown trees that survive.
     ! Note that in both ne and original tiles nindivs of crown trees are different
     ! than in the original tile before the
     t1 => new_land_tile(t0)
     t1%frac = t0%frac*dying_crownwarea
     t0%frac = t0%frac - t1%frac ! and update it to conserve area
     if (is_watch_point()) then
        write(*,*) 'Splitting a disturbed tile'
        __DEBUG3__(f0,t0%frac,t1%frac)
     endif
     ! reset time elapsed since last disturbance in the new tile; do not change age_since_landuse
     t1%vegn%age_since_disturbance = 0.0
     ! change the density of individuals in the cohorts that we split
     do i = 1,t0%vegn%n_cohorts
        associate(cc0=>t0%vegn%cohorts(i), cc1=>t1%vegn%cohorts(i),sp0=>spdata(t0%vegn%cohorts(i)%species))
        if (cc0%layer==1 .and. sp0%lifeform/=FORM_GRASS) then
           cc1%nindivs = ndead(i)*f0/t1%frac
           cc0%nindivs = (cc0%nindivs-ndead(i))*f0/t0%frac
        endif
        end associate
     enddo
     if (is_watch_point()) then
        write(*,*) '#### tile_mortality_ppa: before killing plants ####'
        do i = 1,t0%vegn%n_cohorts
           write(*,'(i2.2)', advance='NO') i
           call dpri('Nindivs0=',t0%vegn%cohorts(i)%nindivs)
           call dpri('Nindivs1=',t1%vegn%cohorts(i)%nindivs)
           call dpri('layer=',t1%vegn%cohorts(i)%layer)
           write(*,*)
        enddo
     endif
     ! kill the individulas affected by natural mortality
     do i = 1,t0%vegn%n_cohorts
        associate(cc0=>t0%vegn%cohorts(i), cc1=>t1%vegn%cohorts(i),sp0=>spdata(t0%vegn%cohorts(i)%species))
        if (cc0%layer==1 .and. sp0%lifeform/=FORM_GRASS) then
           ! kill all canopy plants in disturbed cohort, but do not touch undisturbed one
           call kill_plants_ppa(cc1, t1%vegn, cc1%nindivs, 0.0, &
                       leaf_litt1_C, wood_litt1_C, root_litt1_C, &
                       leaf_litt1_N, wood_litt1_N, root_litt1_N  )
        else
           ! kill understory plants in both tiles
           call kill_plants_ppa(cc1, t1%vegn, ndead(i), 0.0, &
                       leaf_litt1_C, wood_litt1_C, root_litt1_C, &
                       leaf_litt1_N, wood_litt1_N, root_litt1_N  )
           call kill_plants_ppa(cc0, t0%vegn, ndead(i), 0.0, &
                       leaf_litt0_C, wood_litt0_C, root_litt0_C, &
                       leaf_litt0_N, wood_litt0_N, root_litt0_N  )
        endif
        end associate
     enddo
     ! write(*,*)'bwood0=',vegn_tile_bwood(t0%vegn),'bwood1=',vegn_tile_bwood(t1%vegn)
  else
     ! just kill the trees in t0
     if (is_watch_point()) then
        write(*,*) 'NOT splitting a disturbed tile'
     endif
     do i = 1,t0%vegn%n_cohorts
        call kill_plants_ppa(t0%vegn%cohorts(i), t0%vegn, ndead(i), 0.0,&
                       leaf_litt0_C, wood_litt0_C, root_litt0_C, &
                       leaf_litt0_N, wood_litt0_N, root_litt0_N  )
     enddo
  endif

  call add_soil_carbon(t0%soil, t0%vegn, leaf_litt0_C, wood_litt0_C, root_litt0_C, &
                                         leaf_litt0_N, wood_litt0_N, root_litt0_N  )
  if (associated(t1)) &
     call add_soil_carbon(t1%soil, t1%vegn, leaf_litt1_C, wood_litt1_C, root_litt1_C, &
                                            leaf_litt1_N, wood_litt1_N, root_litt1_N  )

  if (is_watch_point()) then
     write(*,*) '#### tile_mortality_ppa output (before relayering cohorts) ####'
     do i = 1,t0%vegn%n_cohorts
        write(*,'(i2.2)', advance='NO') i
        call dpri('Nindivs0=',t0%vegn%cohorts(i)%nindivs)
        call dpri('layer=',t0%vegn%cohorts(i)%layer)
        if(associated(t1)) then
           call dpri('Nindivs1=',t1%vegn%cohorts(i)%nindivs)
           call dpri('layer=',t1%vegn%cohorts(i)%layer)
        endif
        write(*,*)
     enddo
     write(*,*) '#### end of tile_mortality_ppa output ####'
  endif

  if (do_check_conservation) then
     ! + conservation check, part 2: calculate totals in final state, and compare
     ! with previous totals
     call get_tile_water(t0,lmass1,fmass1)
     lmass1 = lmass1*t0%frac; fmass1 = fmass1*t0%frac
     heat1  = land_tile_heat  (t0)*t0%frac
     cmass1 = land_tile_carbon(t0)*t0%frac
     nmass1 = land_tile_nitrogen(t0)*t0%frac
     if (associated(t1)) then
        call get_tile_water(t1,lmass2,fmass2);
        lmass1 = lmass1 + lmass2*t1%frac; fmass1 = fmass1 + fmass2*t1%frac
        heat1  = heat1  + land_tile_heat  (t1)*t1%frac
        cmass1 = cmass1 + land_tile_carbon(t1)*t1%frac
        nmass1 = nmass1 + land_tile_nitrogen(t1)*t1%frac
     endif
     call check_conservation (tag,'liquid water', lmass0, lmass1, water_cons_tol)
     call check_conservation (tag,'frozen water', fmass0, fmass1, water_cons_tol)
     call check_conservation (tag,'carbon'      , cmass0, cmass1, carbon_cons_tol)
     call check_conservation (tag,'nitrogen'    , nmass0, nmass1, nitrogen_cons_tol)
     call check_conservation (tag,'heat content', heat0 , heat1 , heat_cons_tol)
     ! - end of conservation check, part 2
  endif

  call vegn_relayer_cohorts_ppa(t0%vegn)
  call vegn_mergecohorts_ppa(t0%vegn, dheat)
  t0%e_res_2 = t0%e_res_2 - dheat

  if (associated(t1)) then
     call vegn_relayer_cohorts_ppa(t1%vegn)
     call vegn_mergecohorts_ppa(t1%vegn, dheat)
     t1%e_res_2 = t1%e_res_2 - dheat
  endif
end subroutine tile_nat_mortality_ppa

end module vegn_disturbance_mod
