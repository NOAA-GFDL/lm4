module soilc_CENT_type_mod

use fms_mod, only: input_nml_file, check_nml_error, file_exist, close_file, &
            stdlog, mpp_pe, mpp_root_pe, error_mesg, FATAL, NOTE
use time_manager_mod, only: time_type_to_real

use land_constants_mod, only : N_C_TYPES, N_LITTER_POOLS, seconds_per_year
use land_data_mod, only : log_version, lnd
use land_debug_mod, only: land_error_message

use soilc_type_mod, only : soilc_t
use soil_tile_mod, only: soil_tile_type, num_l, soil_theta

implicit none; private

! ---- public items
public :: soilc_CENT_t
public :: new_soilc_CENT
public :: read_soilc_CENT_namelist

! ---- interfces
interface new_soilc_CENT
   module procedure soilc_CENT_ctor
   module procedure soilc_CENT_copy
end interface

! ---- constants
character(len=*), parameter :: module_name = 'soilc_CENT_type_mod'
#include "../../shared/version_variable.inc"

! ----  types
!> @brief soil carbon data container for simplified CENTURY-like soil carbon model
type, extends(soilc_t) :: soilc_CENT_t
  real, dimension(N_C_TYPES, N_LITTER_POOLS) :: litter_century_C !< surface litter (kgC/m2)
  real, allocatable :: &
      fast_soil_C(:), & !< fast soil carbon pool, (kg C/m2), per layer
      slow_soil_C(:), & !< slow soil carbon pool, (kg C/m2), per layer
  ! values for the diagnostic of carbon budget and soil carbon acceleration
      asoil_in(:),    & !< input decomposition rate
      fsc_in(:),      & !< input of fast soil C
      ssc_in(:)         !< input of slow soil C
contains
  procedure :: merge   => merge_CENT   ! merge this soil carbon with another
  procedure :: total_C => total_C_CENT ! returns total C [kgC/m2]
  procedure :: total_N => total_N_CENT ! returns total N [kgN/m2]
  procedure :: rav_C   => rav_C_CENT   ! returns amounts of fas, slow, and (dead) microbial C [kgC/m2]

  procedure :: get_DOC => get_DOC_CENT
  procedure :: get_DON => get_DOC_CENT
  procedure :: get_DIN => get_DIN_CENT
end type soilc_CENT_t

! ---- module data

! namelist
logical :: bulk = .TRUE. !< if True, bulk soil treatment carbon is used, otherwise
                         !! it is treated separately for each layer
real, protected :: K1 = 10.0, K2 = 0.05 !< soil carbon decomposition parameters
real, protected :: tau_lflitt_transfer = 0.0 !< e-folding time scale of leaf litter transfer to soil pools in CENTURY mode, yr; 0 means instant transfer
real, protected :: tau_cwlitt_transfer = 0.0 !< e-folding time scale of coarse wood litter transfer to soil pools in CENTURY mode, yr; 0 means instant transfer

namelist /soil_carbon_CENT_nml/ bulk, K1, K2, tau_lflitt_transfer, tau_cwlitt_transfer

real :: delta_time ! fast (physical) time step, s
real :: dt_fast_yr ! fast (physical) time step, yr (year is defined as 365 days)

contains

!> read namelist
subroutine read_soilc_CENT_namelist()
  integer :: unit         ! unit for namelist i/o
  integer :: io           ! i/o status for the namelist
  integer :: ierr         ! error code, returned by i/o routines

  call log_version(version, module_name, &
  __FILE__)
  read (input_nml_file, nml=soil_carbon_CENT_nml, iostat=io)
  ierr = check_nml_error(io, 'soil_carbon_CENT_nml')
  if (mpp_pe() == mpp_root_pe()) then
     unit=stdlog()
     write(unit, nml=soil_carbon_CENT_nml)
  endif

  delta_time = time_type_to_real(lnd%dt_fast)
  dt_fast_yr = delta_time/seconds_per_year

end subroutine

! constructors
!> @brief Create new (empty) soil carbon representation
!! @return Pointer to new soil carbon data structure
function soilc_CENT_ctor(soil) result(ptr)
  class(soilc_CENT_t), pointer :: ptr
  type(soil_tile_type), intent(in) :: soil !< soil tile data

  allocate(ptr)
  allocate( &
      ptr%fast_soil_C(num_l), &
      ptr%slow_soil_C(num_l), &
      ptr%asoil_in   (num_l), &
      ptr%fsc_in     (num_l), &
      ptr%ssc_in     (num_l)  )

  ptr%litter_century_C(:,:)  = 0.0
  ptr%fast_soil_C(:)         = 0.0
  ptr%slow_soil_C(:)         = 0.0
  ptr%asoil_in(:)            = 0.0
  ptr%fsc_in(:)              = 0.0
  ptr%ssc_in(:)              = 0.0
end function

!> @brief Create a copy of existing soil carbon representation
!! @return Pointer to new soil carbon data structure
function soilc_CENT_copy(soilc) result(ptr)
  type(soilc_CENT_t), pointer :: ptr
  type(soilc_CENT_t), intent(in) :: soilc !< soil carbon data to copy

  allocate(ptr)
  ptr = soilc
end function

!> @brief merge s1 into current soil carbon type s2, with given weights
subroutine merge_CENT(s2,w2,s1,w1)
  class(soilc_CENT_t), intent(inout) :: s2    !< current soil carbon state
  class(soilc_t)     , intent(in)    :: s1    !< soil carbon state to be merged into current
  real               , intent(in)    :: w2,w1 !< merging weights

  real    :: x1, x2 ! normalized relative weights

  ! calculate normalized weights
  x1 = w1/(w1+w2)
  x2 = 1.0 - x1

  select type(s1)
  type is (soilc_CENT_t)
     ! merge soil carbon
     s2%fast_soil_C(:) = s1%fast_soil_C(:)*x1 + s2%fast_soil_C(:)*x2
     s2%slow_soil_C(:) = s1%slow_soil_C(:)*x1 + s2%slow_soil_C(:)*x2
     s2%litter_century_C(:,:) = s1%litter_century_C(:,:)*x1 + s2%litter_century_C(:,:)*x2

     s2%asoil_in(:)    = s1%asoil_in(:)*x1 + s2%asoil_in(:)*x2
     s2%fsc_in(:)      = s1%fsc_in(:)*x1 + s2%fsc_in(:)*x2
     s2%ssc_in(:)      = s1%ssc_in(:)*x1 + s2%ssc_in(:)*x2
  class default
     call land_error_message('merge_CENT: attempt to merge incompatible soil carbon types', FATAL)
  end select
end subroutine


!> @brief Given soil carbon state, return total soil C
!! @return total soil carbon, kgC/m2
real function total_C_CENT(soilc) result(tot_C)
  class(soilc_CENT_t), intent(in)  :: soilc !< soil carbon data structure
  tot_C = sum(soilc%fast_soil_C(:))+sum(soilc%slow_soil_C(:)) &
        + sum(soilc%litter_century_C(:,:))
end function

!> @brief Given soil carbon state, return total soil nitrogen
!! @return total soil nitrogen, kgN/m2
real function total_N_CENT(soilc) result(tot_N)
  class(soilc_CENT_t), intent(in)  :: soilc ! soil carbon data structure
  tot_N = 0.0
end function

!> @brief Given soil carbon state, return carbon amount of litter relevant for surface
!! resistance calculations in legacy treatment of soil surface resistance
subroutine rav_C_CENT(soilc,fast_C,slow_C,dmic_C)
  class(soilc_CENT_t), intent(in)  :: soilc !< soil carbon data structure
  real, intent(out) :: &
     fast_C,    & !< fast litter carbon, [kgC/m2]
     slow_C,    & !< slow litter carbon, [kgC/m2]
     dmic_C       !< mass of dead microbes in litter, [kgC/m2]
! why is that soil C and not litter? originally there were no litter in LM3-like soil
! carbon model, but now that there is, should we switch to using litter? At least when
! tau_*litt_transfer are not zero? or depending on some new namelist parameter?
  fast_C = soilc%fast_soil_C(1)
  slow_C = soilc%slow_soil_C(1)
  dmic_C = 0.0
end subroutine

subroutine get_DOC_CENT(soilC, values)
  class(soilc_CENT_t), intent(in)  :: soilc ! soil carbon data structure
  real,                intent(out) :: values(:,:) ! (N_C_TYPES, num_l) ! [kg C/m^2] dissolved organic carbon

  values(:,:) = 0.0
end subroutine

subroutine get_DIN_CENT(soilC, nitrate, ammonium)
  class(soilc_CENT_t), intent(in)  :: soilc ! soil carbon data structure
  real,                intent(out) :: nitrate(:),ammonium(:) ! [kg N/m^2] dissolved nitrate and ammonium

  nitrate(:)  = 0.0
  ammonium(:) = 0.0
end subroutine

end module
