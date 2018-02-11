module vegn_util_mod

#include "../shared/debug.inc"

use constants_mod,   only : tfreeze
use fms_mod, only : string, WARNING, FATAL

use land_debug_mod, only : is_watch_point, check_var_range, land_error_message, carbon_cons_tol
use soil_carbon_mod, only : N_C_TYPES, C_FAST, deadmic_slow_frac
use soil_tile_mod, only : soil_tile_type, num_l, dz
use soil_util_mod, only : add_soil_carbon
use vegn_data_mod, only : LEAF_OFF, spdata, nspecies, agf_bs, N_limits_live_biomass, &
      min_cohort_nindivs
use vegn_tile_mod, only : vegn_tile_type
use vegn_cohort_mod, only : vegn_cohort_type, plant_C, &
      cohort_root_litter_profile, cohort_root_exudate_profile, init_cohort_hydraulics, &
      init_cohort_allometry_ppa, cohort_can_reproduce

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

subroutine kill_plants_ppa(cc, vegn, ndead, fsmoke, &
                           leaf_litt_C, wood_litt_C, root_litt_C, &
                           leaf_litt_N, wood_litt_N, root_litt_N  )
  type(vegn_cohort_type), intent(inout) :: cc
  type(vegn_tile_type),   intent(inout) :: vegn
  real,                   intent(in)    :: ndead ! number of individuals to kill, indiv./m2
  real,                   intent(in)    :: fsmoke ! fraction of biomass lost to fire, unitless
  real, intent(inout) :: leaf_litt_C(N_C_TYPES), leaf_litt_N(N_C_TYPES) ! accumulated leaf litter, kg C/m2 and kg N/m2
  real, intent(inout) :: wood_litt_C(N_C_TYPES), wood_litt_N(N_C_TYPES) ! accumulated wood litter, kg C/m2 and kg N/m2
  real, intent(inout) :: root_litt_C(num_l, N_C_TYPES),root_litt_N(num_l, N_C_TYPES) ! accumulated root litter per soil layer, kgC/m2 and kg N/m2

  ! ---- local vars
  real :: lost_wood, lost_alive, burned_wood, burned_alive, burned_N
  real :: profile(num_l) ! storage for vertical profile of exudates and root litter
  integer :: l

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
  ! loss to fire: note that we are burning roots, which we probably should not (slm).
  burned_wood  = fsmoke*lost_wood
  burned_alive = fsmoke*lost_alive
  burned_N     = fsmoke*ndead*(cc%leaf_N + cc%seed_N + cc%wood_N + cc%sapwood_N + cc%stored_N + cc%root_N)

  ! loss to the soil pools
  lost_wood  = lost_wood  - burned_wood
  lost_alive = lost_alive - burned_alive

  ! add fire carbon losses to smoke pool
  vegn%csmoke_pool = vegn%csmoke_pool + burned_wood + burned_alive
  vegn%nsmoke_pool = vegn%nsmoke_pool + burned_N

  ! add remaining lost C to soil carbon pools
  associate(sp=>spdata(cc%species))
  leaf_litt_C(:) = leaf_litt_C(:) + [sp%fsc_liv,  1-sp%fsc_liv,  0.0]*(cc%bl+cc%carbon_gain+cc%growth_previous_day)*(1-fsmoke)*ndead
  wood_litt_C(:) = wood_litt_C(:) + [sp%fsc_wood, 1-sp%fsc_wood, 0.0]*(cc%bwood+cc%bsw+cc%bwood_gain)*(1-fsmoke)*agf_bs*ndead
  wood_litt_C(C_FAST) = wood_litt_C(C_FAST)+cc%nsc*(1-fsmoke)*agf_bs*ndead

  leaf_litt_N(:) = leaf_litt_N(:) + [sp%fsc_liv,  1-sp%fsc_liv,  0.0]*cc%leaf_N*(1-fsmoke)*ndead

  wood_litt_N(:) = wood_litt_N(:) + [sp%fsc_wood, 1-sp%fsc_wood, 0.0]*(cc%wood_N+cc%sapwood_N)*(1-fsmoke)*agf_bs*ndead
  wood_litt_N(C_FAST) = wood_litt_N(C_FAST)+cc%stored_N*(1-fsmoke)*agf_bs*ndead

  if ((.not.sp%mortality_kills_seeds).and.cohort_can_reproduce(cc)) then
     ! save seeds in a temporary tile-level buffer
     vegn%drop_seed_C(cc%species) = vegn%drop_seed_C(cc%species) + cc%bseed*(1-fsmoke)*ndead
     vegn%drop_seed_N(cc%species) = vegn%drop_seed_N(cc%species) + cc%seed_N*(1-fsmoke)*ndead
  else
     ! seeds are killed: add them to leaf litter
     leaf_litt_C(:) = leaf_litt_C(:) + [sp%fsc_liv,  1-sp%fsc_liv,  0.0]*cc%bseed*(1-fsmoke)*ndead
     leaf_litt_N(:) = leaf_litt_N(:) + [sp%fsc_liv,  1-sp%fsc_liv,  0.0]*cc%seed_N*(1-fsmoke)*ndead
  endif

  ! accumulate root litter
  call cohort_root_litter_profile(cc, dz, profile)
  do l = 1, num_l
     root_litt_C(l,:) = root_litt_C(l,:) + profile(l)*ndead*(1-fsmoke)*[ &
          sp%fsc_froot    *cc%br + sp%fsc_wood    *(cc%bsw+cc%bwood+cc%bwood_gain)*(1-agf_bs) + cc%nsc*(1-agf_bs), &
          (1-sp%fsc_froot)*cc%br + (1-sp%fsc_wood)*(cc%bsw+cc%bwood+cc%bwood_gain)*(1-agf_bs), &
          0.0]
     root_litt_N(l,:) = root_litt_N(l,:) + profile(l)*ndead*(1-fsmoke)*[ &
             sp%fsc_froot *cc%root_N +    sp%fsc_wood *(cc%sapwood_N+cc%wood_N)*(1-agf_bs) + cc%stored_N*(1-agf_bs), &
          (1-sp%fsc_froot)*cc%root_N + (1-sp%fsc_wood)*(cc%sapwood_N+cc%wood_N)*(1-agf_bs), &
          0.0]
  enddo

  ! C and N content of mycorrhizae and associated buffer reservoirs
  call cohort_root_exudate_profile(cc,dz,profile)
  do l = 1, num_l
     root_litt_C(l,:) = root_litt_C(l,:) + profile(l)*ndead*[0.0, deadmic_slow_frac, (1-deadmic_slow_frac)]* &
          (cc%myc_scavenger_biomass_C+cc%myc_miner_biomass_C+cc%N_fixer_biomass_C+cc%scav_myc_C_reservoir+cc%mine_myc_C_reservoir+cc%N_fixer_C_reservoir)
     root_litt_N(l,:) = root_litt_N(l,:) + profile(l)*ndead*[0.0, deadmic_slow_frac, (1-deadmic_slow_frac)]* &
          (cc%myc_scavenger_biomass_N+cc%myc_miner_biomass_N+cc%N_fixer_biomass_N+cc%scav_myc_N_reservoir+cc%mine_myc_N_reservoir+cc%N_fixer_N_reservoir)
  enddo

  ! reduce the number of individuals in cohort
  cc%nindivs = cc%nindivs-ndead
  end associate
  ! for budget tracking - temporary
  vegn%veg_out = vegn%veg_out + lost_alive + lost_wood + burned_alive + burned_wood
end subroutine kill_plants_ppa

! ============================================================================
subroutine kill_small_cohorts_ppa(vegn,soil)
  type(vegn_tile_type), intent(inout) :: vegn
  type(soil_tile_type), intent(inout) :: soil

  ! ---- local vars
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
     if (vegn%cohorts(i)%nindivs >  min_cohort_nindivs) k=k+1
  enddo

  ! do nothing if there are no cohorts below minimum density
  if (k==vegn%n_cohorts) return

  ! if (k==0) call land_error_message('vegn_mergecohorts_ppa: All cohorts died',WARNING)

  ! exclude cohorts that have zero individuals
  leaf_litt_C = 0 ; wood_litt_C = 0; root_litt_C = 0
  leaf_litt_N = 0 ; wood_litt_N = 0; root_litt_N = 0
  if (k < vegn%n_cohorts) then
     allocate(cc(max(k,1)))
     k=0
     do i = 1,vegn%n_cohorts
        if (vegn%cohorts(i)%nindivs > min_cohort_nindivs) then
           k=k+1
           cc(k) = vegn%cohorts(i)
        else
           call kill_plants_ppa(vegn%cohorts(i), vegn, vegn%cohorts(i)%nindivs, 0.0, &
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
! for each species for which seed_C is greater than zero.
subroutine add_seedlings_ppa(vegn, soil, seed_C, seed_N, germination_factor)
  type(vegn_tile_type), intent(inout) :: vegn
  type(soil_tile_type), intent(inout) :: soil
  real, intent(in) :: seed_C(0:nspecies-1), seed_N(0:nspecies-1)
  real, intent(in), optional :: germination_factor ! additional multiplier for
      ! seed germination, use 0.0 to kill weed seeds on cropland

  type(vegn_cohort_type), pointer :: ccold(:)   ! pointer to old cohort array
  integer :: newcohorts ! number of new cohorts to be created
  real :: failed_seed_C,failed_seed_N !, prob_g, prob_e
  real :: litt_C(N_C_TYPES), litt_N(N_C_types)
  integer :: k ! seedling cohort index
  integer :: i ! species index
  integer :: l ! seedling layer index
  integer :: nlayers ! total number of layers in the canopy
  real, allocatable :: Tv(:)     ! temperature of vegatation in each layer
  real, allocatable :: height(:) ! height of tallest vegetation in each layer
  real    :: germ_f

  germ_f = 1.0
  if (present(germination_factor)) germ_f = 1.0

  if(is_watch_point()) then
     write(*,*)'##### add_seedlings_ppa input #####'
     __DEBUG1__(seed_C)
  endif
  call check_var_range(seed_C,-carbon_cons_tol,HUGE(1.0),'add_seedlings_ppa','seed_C', FATAL)

  newcohorts = count(seed_C>0)
  if (newcohorts == 0) return ! do nothing if no cohorts are ready for reproduction

  if (vegn%n_cohorts+newcohorts>size(vegn%cohorts)) then
     ! increase the size of cohorts array
     ccold => vegn%cohorts
     allocate(vegn%cohorts(vegn%n_cohorts+newcohorts))
     vegn%cohorts(1:vegn%n_cohorts) = ccold(1:vegn%n_cohorts)
     deallocate (ccold)
  endif

  ! store the height and temperature of tallest vegetation within each canopy layer
  nlayers = maxval(vegn%cohorts(1:vegn%n_cohorts)%layer)
  allocate(height(nlayers),Tv(nlayers))
  height = 0.0
  do k = 1,vegn%n_cohorts
     if (vegn%cohorts(k)%height>height(vegn%cohorts(k)%layer)) then
        height(vegn%cohorts(k)%layer) = vegn%cohorts(k)%height
        Tv    (vegn%cohorts(k)%layer) = vegn%cohorts(k)%Tv
     endif
  enddo
  height(1) = HUGE(1.0) ! guard value so that seedlings are always shorter than the first layer

  litt_C(:) = 0.0
  litt_N(:) = 0.0
  ! set up new cohorts
  k = vegn%n_cohorts
  do i = 0,nspecies-1
    if (seed_C(i)<=0) cycle ! no seeds for this species, nothing to do

    k = k+1 ! increment new cohort index
    ! set seedling cohort parameters
    associate (cc => vegn%cohorts(k), sp=>spdata(i))
    cc%species    = i
    cc%status     = LEAF_OFF
    cc%firstlayer = 0
    cc%age        = 0.0
    cc%topyear    = 0.0
    ! BNS: Start cohort with zero stored N and then calculate it after we balance the carbon
    call init_cohort_allometry_ppa(cc, sp%seedling_height, sp%seedling_nsc_frac, 0.0)
    call init_cohort_hydraulics(cc, soil%pars%psi_sat_ref)

    ! added germination probability (prob_g) and establishment probability ((prob_e), Weng 2014-01-06
    cc%nindivs = seed_C(i) * sp%prob_g * germ_f * sp%prob_e/plant_C(cc)
!    __DEBUG3__(cc%age, cc%layer, cc%nindivs)

    ! Nitrogen needs to be adjusted at this point so it's conserved, since seedling N isn't necessarily consistent with initial C vals
    ! cc%total_N is set in init_cohort_allometry_ppa so it should be correct
    if(cc%nindivs>0) &
        cc%stored_N = seed_N(i)*sp%prob_g*sp%prob_e/cc%nindivs - cc%total_N
    ! If nindivs is zero, seedling should be killed by kill_small_cohorts. Otherwise there could be balance problems
    if(cc%stored_N<0 .AND. N_limits_live_biomass) then
       __DEBUG3__(seed_N(i)/cc%nindivs,cc%total_N,cc%stored_N)
       call land_error_message('add_seedlings_ppa: Not enough N in seeds to make seedling',FATAL)
    endif
    cc%total_N = cc%total_N + cc%stored_N

    ! print *,'Reproduction:'
    ! __DEBUG4__(parent%nindivs,parent%seed_C,parent%seed_N,cc%nsc)
    ! __DEBUG4__(cc%nindivs,cc%stored_N,cc%total_N,cc%total_N-cc%stored_N)

    failed_seed_C = (1-sp%prob_g*germ_f*sp%prob_e) * seed_C(i)
    failed_seed_N = (1-sp%prob_g*germ_f*sp%prob_e) * seed_N(i)

    vegn%litter = vegn%litter + failed_seed_C
    litt_C(:) = litt_C(:) + [sp%fsc_liv,1-sp%fsc_liv,0.0]*failed_seed_C
    litt_N(:) = litt_N(:) + [sp%fsc_liv,1-sp%fsc_liv,0.0]*failed_seed_N
    vegn%veg_out = vegn%veg_out + failed_seed_C

    ! Find shortest layer that is still taller than the seedling and assign the seedling
    ! initial layer.
    do l=nlayers,1, -1
       if (cc%height<=height(l)) then
          cc%layer = l
          cc%Tv = Tv(l) ! TODO: make sure that energy is conserved in reproduction
          exit ! from loop
       endif
    enddo

    ! we assume that the newborn cohort is dry; since nindivs of the parent
    ! doesn't change we don't need to do anything with its Wl and Ws to
    ! conserve water (since Wl and Ws are per individual)
    cc%Wl = 0 ; cc%Ws = 0

    end associate   ! F2003
  enddo

  call add_soil_carbon(soil, vegn, leaf_litter_C=litt_C, leaf_litter_N=litt_N)

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

  deallocate(Tv,height)
end subroutine add_seedlings_ppa

end module vegn_util_mod