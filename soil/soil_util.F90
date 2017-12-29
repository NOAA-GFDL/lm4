module soil_util_mod

use constants_mod, only: PI
use fms_mod, only: error_mesg, FATAL

use land_data_mod, only: log_version
use soil_carbon_mod, only: N_C_TYPES, SOILC_CENTURY, SOILC_CENTURY_BY_LAYER, &
     SOILC_CORPSE, C_FAST, C_SLOW, soil_carbon_option, add_litter, add_carbon_to_cohorts
use soil_tile_mod, only: soil_tile_type, dz, num_l, LEAF, CWOOD, N_LITTER_POOLS
use vegn_cohort_mod, only : vegn_cohort_type, cohort_root_litter_profile
use vegn_tile_mod, only: vegn_tile_type
implicit none
private

public :: soil_util_init

public :: add_root_litter
public :: add_root_exudates
public :: add_soil_carbon

interface add_root_litter
   module procedure add_root_litter_0
   module procedure add_root_litter_1
end interface add_root_litter


! ==== module constants ======================================================
character(len=*), parameter :: module_name = 'soil_util'
#include "../shared/version_variable.inc"

contains

! ============================================================================
subroutine soil_util_init(r)
  real, intent(in) :: r
  call log_version(version, module_name, &
  __FILE__)
end subroutine soil_util_init

! ============================================================================
! Spread new root C through profile, using vertical root profile from vegn_uptake_profile.
! If there are negative values in input litter, they are not added to soil pools, but
! accumulated in negativeInput argument. This should be taken into account by the caller
! to avoid carbon conservation issues.
subroutine add_root_litter_0(soil, litterC, negativeInput)
  type(soil_tile_type) , intent(inout) :: soil
  real                 , intent(in)    :: litterC(num_l,N_C_TYPES) ! kg C/(m2 of soil)
  real, optional       , intent(inout) :: negativeInput(N_C_TYPES)

  integer :: l

  do l = 1,num_l
     call add_litter(soil%soil_C(l), litterC(l,:), negativeInput)
  enddo
end subroutine add_root_litter_0


! ============================================================================
! if there are negative values in input litter, they are not added to soil pools, but
! accumulated in negativeInput argument. This should be taken into account by the caller
! to avoid carbon conservation issues.
subroutine add_root_litter_1(soil, cohort, newlitterC, negativeInput)
  type(soil_tile_type)   , intent(inout) :: soil
  type(vegn_cohort_type) , intent(in)    :: cohort
  real                   , intent(in)    :: newlitterC(N_C_TYPES) ! kg C/m2 of tile
  real, optional         , intent(inout) :: negativeInput(N_C_TYPES)

  real    :: profile(num_l)
  integer :: n

  call cohort_root_litter_profile (cohort, dz(1:num_l), profile)
  do n=1,num_l
      call add_litter(soil%soil_C(n), newLitterC*profile(n), negativeInput)
  enddo
end subroutine add_root_litter_1

! ============================================================================
! Spread root exudate C through profile, using vertical root profile from vegn_uptake_profile
! Differs from add_root_litter -- C is distributed through existing cohorts, not deposited as new cohort
subroutine add_root_exudates(soil,exudateC)
  type(soil_tile_type), intent(inout)  :: soil
  real                , intent(in)     :: exudateC(num_l) ! kgC/(m2 of soil)

  integer :: l

  select case (soil_carbon_option)
  case (SOILC_CENTURY)
     soil%fast_soil_C(1) = soil%fast_soil_C(1) + sum(exudateC(:))
     soil%fsc_in(1) = soil%fsc_in(1) + sum(exudateC(:)) ! for budget tracking
  case (SOILC_CENTURY_BY_LAYER)
     do l = 1,num_l
        soil%fast_soil_C(l) = soil%fast_soil_C(l) + exudateC(l)
        soil%fsc_in(l) = soil%fsc_in(l) + exudateC(l) ! for budget tracking
     enddo
  case (SOILC_CORPSE)
     do l = 1,num_l
        call add_carbon_to_cohorts(soil%soil_C(l),litterC=(/exudateC(l),0.0,0.0/))
     enddo
  end select
end subroutine add_root_exudates

! ============================================================================
subroutine add_soil_carbon(soil,leaf_litter,wood_litter,root_litter)
  type(soil_tile_type)   , intent(inout) :: soil
  real, intent(in), optional :: leaf_litter(N_C_TYPES)
  real, intent(in), optional :: wood_litter(N_C_TYPES)
  real, intent(in), optional :: root_litter(num_l,N_C_TYPES)

  integer :: l
  real :: fsc, ssc
  real :: leaf_litt(N_C_TYPES)
  real :: wood_litt(N_C_TYPES)
  real :: root_litt(num_l,N_C_TYPES)

  if (present(leaf_litter)) then
     leaf_litt(:) = leaf_litter(:)
  else
     leaf_litt(:) = 0.0
  endif
  if (present(wood_litter)) then
     wood_litt(:) = wood_litter(:)
  else
     wood_litt(:) = 0.0
  endif
  if (present(root_litter)) then
     root_litt(:,:) = root_litter(:,:)
  else
     root_litt(:,:) = 0.0
  endif

  ! CEL=cellulose (fast); LIG=lignin (slow); this function reasonably assumes
  ! that there are no microbes in litter

  select case (soil_carbon_option)
  case (SOILC_CENTURY)
     fsc = leaf_litt(C_FAST) + wood_litt(C_FAST) + sum(root_litt(:,C_FAST))
     ssc = leaf_litt(C_SLOW) + wood_litt(C_SLOW) + sum(root_litt(:,C_SLOW))
     soil%fast_soil_C(1) = soil%fast_soil_C(1) + fsc
     soil%slow_soil_C(1) = soil%slow_soil_C(1) + ssc
     ! for budget tracking
     soil%fsc_in(1) = soil%fsc_in(1) + fsc
     soil%ssc_in(1) = soil%ssc_in(1) + ssc
  case (SOILC_CENTURY_BY_LAYER)
     fsc = leaf_litt(C_FAST) + wood_litt(C_FAST)
     ssc = leaf_litt(C_SLOW) + wood_litt(C_SLOW)
     soil%fast_soil_C(1) = soil%fast_soil_C(1) + fsc
     soil%slow_soil_C(1) = soil%slow_soil_C(1) + ssc
     ! for budget tracking
     soil%fsc_in(1) = soil%fsc_in(1) + fsc
     soil%ssc_in(1) = soil%ssc_in(1) + ssc
     do l = 1,num_l
        soil%fast_soil_C(l) = soil%fast_soil_C(l) + root_litt(l,C_FAST)
        soil%slow_soil_C(l) = soil%slow_soil_C(l) + root_litt(l,C_SLOW)
        ! for budget tracking
        soil%fsc_in(l) = soil%fsc_in(l) + root_litt(l,C_FAST)
        soil%ssc_in(l) = soil%ssc_in(l) + root_litt(l,C_SLOW)
     enddo
  case (SOILC_CORPSE)
     call borrow_to_negatives(wood_litt,soil%neg_litt_C) ! borrow from wood litter first
     call borrow_to_negatives(leaf_litt,soil%neg_litt_C) ! and from leaf litter second
     call add_litter(soil%leafLitter,       leaf_litt, soil%neg_litt_C)
     call add_litter(soil%coarseWoodLitter, wood_litt, soil%neg_litt_C)
     ! we do not borrow negatives from soil carbon input, because it might do something
     ! bad to slow processes, even if the amounts are tiny.
     do l = 1,num_l
        call add_litter(soil%soil_C(l), root_litt(l,:), soil%neg_litt_C)
     enddo
  end select

contains

  ! given litter and amount of negative litter from previous time step, attempts to borrow
  ! positive carbon to reduce the amount of negativs
  subroutine borrow_to_negatives(litt, negatives)
    real, intent(inout) :: litt(:), negatives(:)

    litt      = litt + negatives
    negatives = min(litt,0.0)
    litt      = max(litt,0.0)
  end subroutine borrow_to_negatives

end subroutine add_soil_carbon

end module soil_util_mod
