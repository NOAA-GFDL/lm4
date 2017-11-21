module vegn_util_mod

#include "../shared/debug.inc"

use constants_mod,   only : tfreeze
use fms_mod, only : string, WARNING, FATAL

use land_debug_mod, only : is_watch_point, check_var_range, land_error_message, carbon_cons_tol
use soil_carbon_mod, only : N_C_TYPES, C_CEL
use soil_tile_mod, only : soil_tile_type, num_l, dz
use soil_util_mod, only : add_soil_carbon
use vegn_data_mod, only : LEAF_OFF, spdata, nspecies, agf_bs, fsc_liv, fsc_wood, fsc_froot
use vegn_tile_mod, only : vegn_tile_type
use vegn_cohort_mod, only : vegn_cohort_type, biomass_of_individual, &
      cohort_root_litter_profile, cohort_root_exudate_profile, init_cohort_hydraulics, &
      init_cohort_allometry_ppa

implicit none
private

public :: kill_plants_ppa
public :: kill_small_cohorts_ppa
public :: add_seedlings_ppa

! ==== module constants ======================================================
character(len=*), parameter :: module_name = 'vegn_util'
#include "../shared/version_variable.inc"

! TODO: possibly move all definition of cpw,clw,csw in one place
real, parameter :: &
     clw = 4218.0, & ! specific heat of water (liquid)
     csw = 2106.0    ! specific heat of water (ice)

contains

! ==============================================================================
! given a cohort, number of individuals to kill, and fraction of smoke, kills
! specified number of individuals in cohort, putting their biomass into litter
! and smoke pools.
!
! In vegn, this subroutine updates csmoke and intermediate heat and water pools (drop_*)
! that are taken into account in energy/mass solver on the next time step.
! *_litt arrays are incremented by this subroutine, so that consecutive calls
! calculate cumulative amount of litter.

! NOTE that in contrast to LM3 mortality calculation, living biomass is
! included in losses; does it mean that we double-counting losses in
! maintenance and here?

! TODO: ask Ensheng and Elena about possible double-counting of live biomass losses.

! we are burning all of biomass (like we did in LM3), including the below-ground
! part. Perhaps we should leave the below-ground part alone, and add it to soil
! carbon pools?

! TODO: ask Elena about burning of below-ground biomass

subroutine kill_plants_ppa(cc, vegn, ndead, fsmoke, leaf_litt, wood_litt, root_litt)
  type(vegn_cohort_type), intent(inout) :: cc
  type(vegn_tile_type),   intent(inout) :: vegn
  real,                   intent(in)    :: ndead ! number of individuals to kill, indiv./m2
  real,                   intent(in)    :: fsmoke ! fraction of biomass lost to fire, unitless
  real, intent(inout) :: leaf_litt(N_C_TYPES) ! accumulated leaf litter, kg C/m2
  real, intent(inout) :: wood_litt(N_C_TYPES) ! accumulated wood litter, kg C/m2
  real, intent(inout) :: root_litt(num_l, N_C_TYPES) ! accumulated root litter per soil layer, kgC/m2

  ! ---- local vars
  real :: lost_wood, lost_alive, burned_wood, burned_alive
  real :: profile(num_l) ! storage for vertical profile of exudates and root litter
  integer :: l
  real :: bp, bn ! positive and negative amounts of litter, used in adjusting litter pools if one is negative 

  call check_var_range(ndead,  0.0, cc%nindivs, 'kill_plants_ppa', 'ndead',  FATAL)
  call check_var_range(fsmoke, 0.0, 1.0,        'kill_plants_ppa', 'fsmoke', FATAL)

  ! water from dead trees goes to intermediate buffers, to be added to the
  ! precipitation reaching ground on the next physical time step
  vegn%drop_wl = vegn%drop_wl + cc%wl*ndead
  vegn%drop_ws = vegn%drop_ws + cc%ws*ndead
  vegn%drop_hl = vegn%drop_hl + clw*cc%wl*ndead*(cc%Tv-tfreeze)
  vegn%drop_hs = vegn%drop_hs + csw*cc%ws*ndead*(cc%Tv-tfreeze)

  ! calculate total carbon losses, kgC/m2
  lost_wood  = ndead * (cc%bwood+cc%bsw+cc%bwood_gain)
  lost_alive = ndead * (cc%bl+cc%br+cc%blv+cc%bseed+cc%nsc+cc%carbon_gain+cc%growth_previous_day)
  ! loss to fire
  burned_wood  = fsmoke*lost_wood
  burned_alive = fsmoke*lost_alive
  ! loss to the soil pools
  lost_wood  = lost_wood  - burned_wood
  lost_alive = lost_alive - burned_alive

  ! add fire carbon losses to smoke pool
  vegn%csmoke_pool = vegn%csmoke_pool + burned_wood + burned_alive

  ! add remaining lost C to soil carbon pools
  leaf_litt(:) = leaf_litt(:) + [fsc_liv,  1-fsc_liv,  0.0]*(cc%bl+cc%bseed+cc%carbon_gain+cc%growth_previous_day)*(1-fsmoke)*ndead
  wood_litt(:) = wood_litt(:) + [fsc_wood, 1-fsc_wood, 0.0]*(cc%bwood+cc%bsw+cc%bwood_gain)*(1-fsmoke)*agf_bs*ndead
  wood_litt(C_CEL) = wood_litt(C_CEL)+cc%nsc*(1-fsmoke)*agf_bs*ndead
  call cohort_root_litter_profile(cc, dz, profile)
  do l = 1, num_l
     root_litt(l,:) = root_litt(l,:) + profile(l)*ndead*(1-fsmoke)*(/ &
          fsc_froot    *cc%br + fsc_wood    *(cc%bsw+cc%bwood+cc%bwood_gain)*(1-agf_bs) + cc%nsc*(1-agf_bs), &
          (1-fsc_froot)*cc%br + (1-fsc_wood)*(cc%bsw+cc%bwood+cc%bwood_gain)*(1-agf_bs), &
          0.0/)
  enddo

  ! leaf_litt can be below zero if biomasses are very small and carbon_gain is negative:
  ! try to borrow carbon from wood litter.
  do l = 1,N_C_TYPES
     if (leaf_litt(l)<0) then
        wood_litt(l) = wood_litt(l) + leaf_litt(l)
        leaf_litt(l) = 0.0
     endif
  enddo
  call check_var_range(wood_litt, 0.0, HUGE(1.0), 'kill_plants_ppa', 'wood_litt',  WARNING)
  if (any(wood_litt<0.0)) then
     ! if some wood litter components are negative, try to borrow carbon
     ! from positive components so that the total carbon is conserved
     bp = 0.0; bn=0.0
     do l = 1, N_C_TYPES
        if (wood_litt(l)>0) bp = bp+wood_litt(l)
        if (wood_litt(l)<0) bn = bn+abs(wood_litt(l))
     enddo
     if (bp<bn) call land_error_message(&
        'kill_plants_ppa: total wood litter amount is negative ('//string(sum(wood_litt))//')', FATAL)
     do l = 1, N_C_TYPES
        if (wood_litt(l)>0) wood_litt(l) = wood_litt(l)+(bp-bn)/bp
        if (wood_litt(l)<0) wood_litt(l) = 0.0
     enddo
  endif

  ! reduce the number of individuals in cohort
  cc%nindivs = cc%nindivs-ndead

  ! for budget tracking - temporary
  vegn%veg_out = vegn%veg_out + lost_alive + lost_wood + burned_alive + burned_wood
end subroutine kill_plants_ppa

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

!  write(*,*)'kill_small_cohorts_ppa n_cohorts before: ', vegn%n_cohorts

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
           call kill_plants_ppa(vegn%cohorts(i), vegn, vegn%cohorts(i)%nindivs, 0.0, &
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
     write(*,*) '##### kill_small_cohorts_ppa output #####'
     __DEBUG1__(vegn%n_cohorts)
     __DEBUG1__(vegn%cohorts%nindivs)
     __DEBUG1__(vegn%cohorts%Wl)
     __DEBUG1__(vegn%cohorts%Ws)
     __DEBUG1__(vegn%cohorts%mcv_dry)
     __DEBUG1__(vegn%cohorts%Tv)
  endif
!  write(*,*)'kill_small_cohorts_ppa n_cohorts after: ', vegn%n_cohorts

end subroutine kill_small_cohorts_ppa

! ============================================================================
! Given seed biomass for each species (kgC per m2 of tile), add a seedling cohort
! for each species for which bseed is greater than zero.
subroutine add_seedlings_ppa(vegn, soil, bseed)
  type(vegn_tile_type), intent(inout) :: vegn
  type(soil_tile_type), intent(inout) :: soil
  real, intent(in) :: bseed(0:nspecies-1)

  type(vegn_cohort_type), pointer :: ccold(:)   ! pointer to old cohort array
  real    :: failed_seeds
  real    :: litt(N_C_TYPES)
  integer :: newcohorts ! number of new cohorts to be created
  integer :: k ! seedling cohort index
  integer :: i ! species index
  real    :: Tv ! temperature assigned to seedlings

  if(is_watch_point()) then
     write(*,*)'##### add_seedlings_ppa input #####'
     __DEBUG1__(bseed)
  endif
  call check_var_range(bseed,-carbon_cons_tol,HUGE(1.0),'add_seedlings_ppa','bseed', FATAL)

  newcohorts = count(bseed>0)
  if (newcohorts == 0) return ! do nothing if no cohorts are ready for reproduction

  if (vegn%n_cohorts+newcohorts>size(vegn%cohorts)) then
     ! increase the size of cohorts array
     ccold => vegn%cohorts
     allocate(vegn%cohorts(vegn%n_cohorts+newcohorts))
     vegn%cohorts(1:vegn%n_cohorts) = ccold(1:vegn%n_cohorts)
     deallocate (ccold)
  endif

  ! use the temperature of the last (shortest) existing cohort for seedlings
  Tv = vegn%cohorts(vegn%n_cohorts)%Tv

  litt(:) = 0.0
  ! set up new cohorts
  k = vegn%n_cohorts
  do i = 0,nspecies-1
    if (bseed(i)<=0) cycle ! no seeds for this species, nothing to do

    k = k+1 ! increment new cohort index
    ! set seedling cohort parameters
    associate (cc => vegn%cohorts(k), sp=>spdata(i))
    cc%species    = i
    cc%status     = LEAF_OFF
    cc%firstlayer = 0
    cc%age        = 0.0
    cc%topyear    = 0.0
    call init_cohort_allometry_ppa(cc, sp%seedling_height, sp%seedling_nsc_frac)
    call init_cohort_hydraulics(cc, soil%pars%psi_sat_ref)

    ! added germination probability (prob_g) and establishment probability ((prob_e), Weng 2014-01-06
    cc%nindivs = bseed(i) * sp%prob_g * sp%prob_e/biomass_of_individual(cc)
!    __DEBUG3__(cc%age, cc%layer, cc%nindivs)

    failed_seeds = (1.0-sp%prob_g*sp%prob_e) * bseed(i)
    vegn%litter = vegn%litter + failed_seeds
    litt(:) = litt(:) + (/fsc_liv,1-fsc_liv,0.0/)*failed_seeds
    vegn%veg_out = vegn%veg_out + failed_seeds

    ! we assume that the newborn cohort is dry; since nindivs of the parent
    ! doesn't change we don't need to do anything with its Wl and Ws to
    ! conserve water (since Wl and Ws are per individual)
    cc%Wl = 0 ; cc%Ws = 0
    ! TODO: make sure that energy is conserved in reproduction
    cc%Tv = Tv

    end associate   ! F2003
  enddo
  call add_soil_carbon(soil, leaf_litter=litt)

  vegn%n_cohorts = k
  if(is_watch_point()) then
     write(*,*)'##### add_seedlings_ppa output #####'
     __DEBUG2__(newcohorts, vegn%n_cohorts)
     do k = vegn%n_cohorts-newcohorts+1, vegn%n_cohorts
        write(*,'(a,i2.2)',advance='NO') 'cohort=', k
        call dpri(' nindivs=', vegn%cohorts(k)%nindivs)
        call dpri(' species=', vegn%cohorts(k)%species)
        call dpri(' Tv=',      vegn%cohorts(k)%Tv)
        write(*,*)
     enddo
  endif
end subroutine add_seedlings_ppa

end module vegn_util_mod
