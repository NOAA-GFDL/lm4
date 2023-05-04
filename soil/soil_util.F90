module soil_util_mod

use constants_mod, only: PI
use fms_mod, only: error_mesg, FATAL

use land_constants_mod, only: N_LITTER_POOLS, LITT_LEAF, LITT_CWOOD, &
     N_C_TYPES, C_FAST, C_SLOW
use land_data_mod, only: log_version
use soilc_CENT_type_mod, only: soilc_CENT_t
use soil_carbon_mod, only: soilc_t, soilc_CORPSE_t, &
     SOILC_CENTURY, SOILC_CENTURY_BY_LAYER, &
     soil_carbon_option, add_litter, add_C_N_to_rhizosphere
use soil_tile_mod, only: soil_tile_type, dz, num_l
use vegn_cohort_mod, only : vegn_cohort_type, cohort_root_exudate_profile
use vegn_data_mod, only: spdata, tau_lflitt_transfer, tau_cwlitt_transfer
use vegn_tile_mod, only: vegn_tile_type
implicit none
private

public :: soil_util_init

public :: add_root_litter
public :: add_root_exudates
public :: add_soil_carbon
public :: rhizosphere_frac

! ==== module constants ======================================================
character(len=*), parameter :: module_name = 'soil_util'
#include "../shared/version_variable.inc"

! ==== module variables ======================================================
real :: r_rhiz = -1.0

contains

subroutine soil_util_init(r)
  real, intent(in) :: r
  call log_version(version, module_name, &
  __FILE__)

  r_rhiz = r
end subroutine soil_util_init

! ============================================================================
! Spread new root C through profile
subroutine add_root_litter(soil, vegn, litterC, litterN, negativeInputC, negativeInputN)
  class(soilc_t)       , intent(inout) :: soil
  type(vegn_tile_type) , intent(in)    :: vegn
  real, intent(in) :: litterC(num_l,N_C_TYPES) ! kg C/(m2 of soil)
  real, intent(in) :: litterN(num_l,N_C_TYPES) ! kg C/(m2 of soil)
  real, intent(inout), optional :: negativeInputC(N_C_TYPES), negativeInputN(N_C_TYPES)

  integer :: k
  real :: rhiz_frac(num_l)  ! fraction of rhizosphere in each layer

  select type (soil)
  class is (soilc_CORPSE_t)
     call rhizosphere_frac(vegn, rhiz_frac)
     do k = 1,num_l
        call add_litter(soil%org_matter(k), litterC(k,:), litterN(k,:), rhiz_frac(k), &
                        negativeInputC, negativeInputN)
     enddo
  class default
     call error_mesg('add_root_litter','called for incorrect soil_carbon_option -- this should never happen', FATAL)
  end select

end subroutine add_root_litter

! ============================================================================
! Spread root exudate C through profile, using vertical root profile from vegn_uptake_profile
! Differs from add_root_litter -- C is distributed through existing cohorts, not deposited as new cohort
subroutine add_root_exudates(soil,exudateC,exudateN,ammonium,nitrate)
    class(soilc_t),       intent(inout)  :: soil
    real,dimension(num_l),intent(in) :: exudateC,exudateN
    real,dimension(num_l),intent(in),optional :: ammonium,nitrate

  integer :: k
  real,dimension(num_l) :: NH4,NO3
  real :: fsc

  NH4(:)=0.0
  NO3(:)=0.0
  if(present(ammonium)) NH4=ammonium
  if(present(nitrate))  NO3=nitrate

  ! ignore exudateN for CENTURY-like soil carbon options
  select type (soil)
  class is (soilC_CENT_t)
     select case(soil_carbon_option)
     case (SOILC_CENTURY)
        fsc = sum(exudateC(:))
        soil%fast_soil_C(1) = soil%fast_soil_C(1) + fsc
        soil%fsc_in(1)      = soil%fsc_in(1)      + fsc ! for soil carbon equilibration
     case (SOILC_CENTURY_BY_LAYER)
        do k = 1, num_l
           soil%fast_soil_C(k) = soil%fast_soil_C(k) + exudateC(k)
           soil%fsc_in(k)      = soil%fsc_in(k)      + exudateC(k) ! for soil carbon equilibration
        enddo
     end select
  class is (soilc_CORPSE_t)
     do k=1,num_l
        call add_C_N_to_rhizosphere(soil%org_matter(k),   &
                                newCarbon=[exudateC(k),0.0,0.0], &
                                newNitrogen=[exudateN(k),0.0,0.0]  )
        soil%org_matter(k)%ammonium = soil%org_matter(k)%ammonium+NH4(k)
        soil%org_matter(k)%nitrate = soil%org_matter(k)%nitrate+NO3(k)
     enddo
  end select
end subroutine add_root_exudates

! ============================================================================
subroutine add_soil_carbon(soil,vegn,leaf_litter_C,wood_litter_C,root_litter_C,&
                                     leaf_litter_N,wood_litter_N,root_litter_N)
  class(soilc_t),       intent(inout) :: soil
  type(vegn_tile_type), intent(inout) :: vegn
  real, intent(in), optional :: leaf_litter_C(N_C_TYPES)
  real, intent(in), optional :: wood_litter_C(N_C_TYPES)
  real, intent(in), optional :: root_litter_C(num_l,N_C_TYPES)
  real, intent(in), optional :: leaf_litter_N(N_C_TYPES)
  real, intent(in), optional :: wood_litter_N(N_C_TYPES)
  real, intent(in), optional :: root_litter_N(num_l,N_C_TYPES)

  integer :: l
  real :: fsc, ssc
  real :: leaf_litt_C(N_C_TYPES), leaf_litt_N(N_C_TYPES)
  real :: wood_litt_C(N_C_TYPES), wood_litt_N(N_C_TYPES)
  real :: root_litt_C(num_l,N_C_TYPES), root_litt_N(num_l,N_C_TYPES)
  real :: rhiz_frac(num_l)

  if (present(leaf_litter_C)) then
     leaf_litt_C(:) = leaf_litter_C(:)
  else
     leaf_litt_C(:) = 0.0
  endif
  if (present(wood_litter_C)) then
     wood_litt_C(:) = wood_litter_C(:)
  else
     wood_litt_C(:) = 0.0
  endif
  if (present(root_litter_C)) then
     root_litt_C(:,:) = root_litter_C(:,:)
  else
     root_litt_C(:,:) = 0.0
  endif
  if (present(leaf_litter_N)) then
     leaf_litt_N(:) = leaf_litter_N(:)
  else
     leaf_litt_N(:) = 0.0
  endif
  if (present(wood_litter_N)) then
     wood_litt_N(:) = wood_litter_N(:)
  else
     wood_litt_N(:) = 0.0
  endif
  if (present(root_litter_N)) then
     root_litt_N(:,:) = root_litter_N(:,:)
  else
     root_litt_N(:,:) = 0.0
  endif

  ! CEL=cellulose (fast); LIG=lignin (slow); this function reasonably assumes
  ! that there are no microbes in litter

  select type (soil)
  class is (soilc_CENT_t)
     select case (soil_carbon_option)
     case (SOILC_CENTURY)
        if (tau_cwlitt_transfer>0.or.tau_lflitt_transfer>0) then
           ! put litterfall in litter pools
           soil%litter_century_C(:,LITT_LEAF)  = soil%litter_century_C(:,LITT_LEAF)  + leaf_litt_C(:)
           soil%litter_century_C(:,LITT_CWOOD) = soil%litter_century_C(:,LITT_CWOOD) + wood_litt_C(:)
           fsc = sum(root_litt_C(:,C_FAST))
           ssc = sum(root_litt_C(:,C_SLOW))
        else
           ! add litterfall to soil carbon directly. This is mostly to preserve bitwise
           ! reproducibility with older code versions
           fsc = leaf_litt_C(C_FAST) + wood_litt_C(C_FAST) + sum(root_litt_C(:,C_FAST))
           ssc = leaf_litt_C(C_SLOW) + wood_litt_C(C_SLOW) + sum(root_litt_C(:,C_SLOW))
        endif
        soil%fast_soil_C(1) = soil%fast_soil_C(1) + fsc
        soil%slow_soil_C(1) = soil%slow_soil_C(1) + ssc
        ! for budget tracking
        soil%fsc_in(1) = soil%fsc_in(1) + fsc
        soil%ssc_in(1) = soil%ssc_in(1) + ssc
     case (SOILC_CENTURY_BY_LAYER)
        if (tau_cwlitt_transfer>0.or.tau_lflitt_transfer>0) then
           ! put litterfall in litter pools
           soil%litter_century_C(:,LITT_LEAF)  = soil%litter_century_C(:,LITT_LEAF)  + leaf_litt_C(:)
           soil%litter_century_C(:,LITT_CWOOD) = soil%litter_century_C(:,LITT_CWOOD) + wood_litt_C(:)
           fsc = 0.0; ssc = 0.0
        else
           ! add litterfall to soil carbon directly. This is mostly to preserve bitwise
           ! reproducibility with older code versions
           fsc = leaf_litt_C(C_FAST) + wood_litt_C(C_FAST)
           ssc = leaf_litt_C(C_SLOW) + wood_litt_C(C_SLOW)
        endif
        soil%fast_soil_C(1) = soil%fast_soil_C(1) + fsc
        soil%slow_soil_C(1) = soil%slow_soil_C(1) + ssc
        ! for budget tracking
        soil%fsc_in(1) = soil%fsc_in(1) + fsc
        soil%ssc_in(1) = soil%ssc_in(1) + ssc
        do l = 1,num_l
           soil%fast_soil_C(l) = soil%fast_soil_C(l) + root_litt_C(l,C_FAST)
           soil%slow_soil_C(l) = soil%slow_soil_C(l) + root_litt_C(l,C_SLOW)
           ! for budget tracking
           soil%fsc_in(l) = soil%fsc_in(l) + root_litt_C(l,C_FAST)
           soil%ssc_in(l) = soil%ssc_in(l) + root_litt_C(l,C_SLOW)
        enddo
     end select
  class is (soilc_CORPSE_t)
     call borrow_to_negatives(wood_litt_C,soil%neg_litt_C) ! borrow from wood litter first
     call borrow_to_negatives(wood_litt_N,soil%neg_litt_N) ! borrow from wood litter first
     call borrow_to_negatives(leaf_litt_C,soil%neg_litt_C) ! and from leaf litter second
     call borrow_to_negatives(leaf_litt_N,soil%neg_litt_N) ! and from leaf litter second
     call add_litter(soil%litter_corpse(LITT_LEAF),  leaf_litt_C, leaf_litt_N, negativeInputC=soil%neg_litt_C, negativeInputN=soil%neg_litt_N)
     call add_litter(soil%litter_corpse(LITT_CWOOD), wood_litt_C, wood_litt_N, negativeInputC=soil%neg_litt_C, negativeInputN=soil%neg_litt_N)
     call rhizosphere_frac(vegn, rhiz_frac)
     do l = 1,num_l
        call add_litter(soil%org_matter(l), root_litt_C(l,:), root_litt_N(l,:), rhiz_frac(l), &
                        negativeInputC=soil%neg_litt_C, negativeInputN=soil%neg_litt_N)
     enddo
  class default
     call error_mesg('add_soil_carbon','unrecognized soil carbon option -- this should never happen', FATAL)
  end select

  ! accumulate litterfall diagnostics: it is sent to diag and then reset at every time step
  vegn%litterfall_C(:,LITT_LEAF)  = vegn%litterfall_C(:,LITT_LEAF)  + leaf_litt_C(:)
  vegn%litterfall_C(:,LITT_CWOOD) = vegn%litterfall_C(:,LITT_CWOOD) + wood_litt_C(:)

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

! ============================================================================
subroutine rhizosphere_frac(vegn, rhiz_frac)
  type(vegn_tile_type) , intent(in)    :: vegn
  real, intent(out) :: rhiz_frac(:)  ! volumentric fraction of rhizosphere

  real :: rhiz_vol(num_l)  ! volume of rhizosphere in each layer, m3/m2
  integer :: i

  ! first calculate the volume of rhizosphere
  rhiz_vol(:) = 0.0
  do i = 1,vegn%n_cohorts
     associate(cc=>vegn%cohorts(i),sp=>spdata(vegn%cohorts(i)%species))
     rhiz_vol(:) = rhiz_vol(:) + &
         PI*((r_rhiz+sp%root_r)**2-sp%root_r**2)*cc%root_length(1:num_l)*cc%nindivs
     end associate
  enddo
  ! rhiz_frac(1:num_l) = min(1.0,rhiz_vol(:)/dz(1:num_l))
  ! If root_length is m/m3, then we should not divide by dz here
  rhiz_frac(1:num_l) = min(1.0,rhiz_vol(:))
end subroutine rhizosphere_frac

end module soil_util_mod
