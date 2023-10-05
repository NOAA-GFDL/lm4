module soil_BGC_SIMPLE_type_mod

use fms_mod, only: input_nml_file, check_nml_error, file_exist, close_file, &
        stdlog, mpp_pe, mpp_root_pe, error_mesg, FATAL, NOTE
use time_manager_mod, only: time_type_to_real

use land_constants_mod, only : N_LITTER_POOLS, seconds_per_year, &
        N_C_TYPES, C_FAST, C_SLOW, C_MIC, LITT_LEAF, LITT_CWOOD

use land_data_mod, only : log_version, lnd
use land_debug_mod, only: land_error_message

use tile_diag_buff_mod, only : diag_buff_type
use tile_diag_base_mod, only : set_default_diag_filter, &
        register_tiled_diag_field, send_tile_data, add_tiled_diag_field_alias, CMOR_NAME

use soil_BGC_type_mod, only : soil_BGC_t, deplete_pool
use soil_BGC_util_mod, only : register_soilc_diag_fields, register_litter_diag_fields, &
        register_litter_soilc_diag_fields
use soil_tile_mod, only: soil_tile_type, num_l, soil_theta, dz
use vegn_tile_mod, only: vegn_tile_type

implicit none; private

! ---- public items
public :: soil_BGC_SIMPLE_t
public :: new_soilc_SIMPLE
public :: read_soil_BGC_SIMPLE_namelist, soil_BGC_diag_init_SIMPLE
public :: save_equilibration_data

! ---- interfces
interface new_soilc_SIMPLE
   module procedure soilc_SIMPLE_ctor
   module procedure soilc_SIMPLE_copy
end interface

! ---- constants
character(len=*), parameter :: module_name = 'soil_BCG_SIMPLE_type_mod'
#include "../../shared/version_variable.inc"

! ----  types
!> @brief soil carbon data container for simple soil carbon model
type, extends(soil_BGC_t) :: soil_BGC_SIMPLE_t
  real, dimension(N_C_TYPES, N_LITTER_POOLS) :: litter_SIMPLE_C !< surface litter (kgC/m2)
  real, allocatable :: &
      fast_soil_C(:), & !< fast soil carbon pool, (kg C/m2), per layer
      slow_soil_C(:), & !< slow soil carbon pool, (kg C/m2), per layer
  ! values for the diagnostic of carbon budget and soil carbon acceleration
      asoil_in(:),    & !< input decomposition rate
      fsc_in(:),      & !< input of fast soil C
      ssc_in(:)         !< input of slow soil C
contains
  procedure :: merge   => merge_SIMPLE   ! merge this soil carbon with another
  procedure :: total_C => total_C_SIMPLE ! returns total C [kgC/m2]
  procedure :: total_N => total_N_SIMPLE ! returns total N [kgN/m2]
  procedure :: rav_C   => rav_C_SIMPLE   ! returns amounts of fast, slow, and (dead) microbial C [kgC/m2]
                                       ! for litter evaporation resistance calculations
  procedure :: get_DOC => get_zero_2D
  procedure :: get_DON => get_zero_2D
  procedure :: get_nit => get_zero_1D
  procedure :: get_amm => get_zero_1D
  procedure :: get_littC => get_littC_SIMPLE

  procedure :: add_soil_carbon   => add_soil_carbon_SIMPLE
  procedure :: add_root_exudates => add_root_exudates_SIMPLE
  procedure :: burn_litter_frac  => burn_litter_frac_SIMPLE
  procedure :: tracer_leaching   => tracer_leaching_SIMPLE

  procedure :: deposit_N              => deposit_N_SIMPLE
  procedure :: active_root_N_uptake   => active_root_N_uptake_SIMPLE
  procedure :: myc_scavenger_N_uptake => myc_scavenger_N_uptake_SIMPLE
  procedure :: myc_miner_N_uptake     => myc_miner_N_uptake_SIMPLE

  procedure :: spend_intermediate_pools => spend_intermediate_pools_SIMPLE
  procedure :: dsdt              => dsdt_SIMPLE
  procedure :: step3             => step3_SIMPLE
  procedure :: redistribute_peat_carbon => redistribute_peat_carbon_SIMPLE
end type soil_BGC_SIMPLE_t

! ---- module data

! namelist
logical :: bulk = .TRUE. !< if True, bulk soil treatment carbon is used, otherwise
                         !! it is treated separately for each layer
real, protected :: K1 = 10.0, K2 = 0.05 !< soil carbon decomposition parameters
real, protected :: tau_lflitt_transfer = 0.0 !< e-folding time scale of leaf litter transfer to soil pools in SIMPLE mode, yr; 0 means instant transfer
real, protected :: tau_cwlitt_transfer = 0.0 !< e-folding time scale of coarse wood litter transfer to soil pools in SIMPLE mode, yr; 0 means instant transfer

logical, protected :: save_equilibration_data = .FALSE. !< if TRUE, information for
                         !! soil BGC equilibration acceleration is saved to disk

namelist /soil_BGC_SIMPLE_nml/ bulk, K1, K2, tau_lflitt_transfer, tau_cwlitt_transfer, &
     save_equilibration_data

real :: delta_time ! fast (physical) time step, s
real :: dt_fast_yr ! fast (physical) time step, yr (year is defined as 365 days)

! diag fields ID
integer :: id_fsc, id_ssc, id_total_soil_C, id_asoil, id_rsoil
integer, dimension(N_C_TYPES) :: id_soil_C, id_rsoil_C
integer, dimension(N_LITTER_POOLS) :: id_litter_total_C
integer, dimension(N_LITTER_POOLS,N_C_TYPES) :: id_litter_C

! diag IDs of CMOR variables
integer :: id_csoil, id_nSoil, id_rh, &
    id_csoilfast, id_csoilmedium, id_csoilslow, id_cSoilLevels, &
    id_cLitter, id_cLitterCwd, id_cLitterLeaf

contains

!> read namelist
subroutine read_soil_BGC_SIMPLE_namelist()
  integer :: unit         ! unit for namelist i/o
  integer :: io           ! i/o status for the namelist
  integer :: ierr         ! error code, returned by i/o routines

  call log_version(version, module_name, &
  __FILE__)
  read (input_nml_file, nml=soil_BGC_SIMPLE_nml, iostat=io)
  ierr = check_nml_error(io, 'soil_BGC_SIMPLE_nml')
  if (mpp_pe() == mpp_root_pe()) then
     unit=stdlog()
     write(unit, nml=soil_BGC_SIMPLE_nml)
  endif

  delta_time = time_type_to_real(lnd%dt_fast)
  dt_fast_yr = delta_time/seconds_per_year
end subroutine

!> @brief Register diagnostic fields
subroutine soil_BGC_diag_init_SIMPLE(id_ug, id_zfull)
  integer,intent(in)  :: id_ug    !< Unstructured axis id
  integer,intent(in)  :: id_zfull !< Vertical (depth) axis id

  character(*), parameter :: diag_mod_name = 'soil'
  integer :: axes(2)

  ! define array of axis indices
  axes = [ id_ug,id_zfull ]

  ! set the default sub-sampling filter for the fields below
  call set_default_diag_filter('soil')

  id_fsc = register_tiled_diag_field ( diag_mod_name, 'fsc', axes(1:1),  &
       lnd%time, 'total fast soil carbon, including soil and litter pools', 'kg C/m2', missing_value=-100.0 )
  id_ssc = register_tiled_diag_field ( diag_mod_name, 'ssc', axes(1:1),  &
       lnd%time, 'total slow soil carbon, including soil and litter pools', 'kg C/m2', missing_value=-100.0 )
  id_soil_C(:) = register_soilc_diag_fields(diag_mod_name, '<ctype>_soil_C', &
       axes, lnd%time, '<ctype> soil carbon', 'kg C/m3', missing_value=-100.0 )
  id_total_soil_C = register_tiled_diag_field ( diag_mod_name, 'tot_soil_C', axes(1:1),  &
       lnd%time, 'total carbon, including soil and litter pools', 'kg C/m2', missing_value=-100.0 )
  id_litter_C(:,:) = register_litter_soilc_diag_fields ( diag_mod_name, '<ltype>litt_<ctype>_C', &
       axes(1:1), lnd%time, '<ctype> <ltype> litter carbon', 'kg C/m2', missing_value=-100.0 )
  id_litter_total_C(:) = register_litter_diag_fields ( diag_mod_name, '<ltype>litt_C', &
       axes(1:1),  lnd%time, '<ltype> litter total carbon', 'kg C/m2', missing_value=-100.0 )

  id_asoil = register_tiled_diag_field ( diag_mod_name, 'asoil', &
       axes(1:1), lnd%time, 'aerobic activity modifier', &
       missing_value=-100.0 )
  id_rsoil = register_tiled_diag_field ( diag_mod_name, 'rsoil',  &
       axes(1:1), lnd%time, 'soil respiration', 'kg C/(m2 year)', &
       missing_value=-100.0 )

  id_rsoil_C(:) = register_soilc_diag_fields( diag_mod_name, 'rsoil_<ctype>', &
       axes, lnd%time, '<ctype> soil carbon respiration', 'kg C/(m3 year)', &
       missing_value=-100.0 )


  ! set the default sub-sampling filter for the fields below
  call set_default_diag_filter('land')

  id_rh = register_tiled_diag_field ( CMOR_NAME, 'rh', [ id_ug ], &
       lnd%time, 'Heterotrophic Respiration', 'kg m-2 s-1', missing_value=-1.0, &
       standard_name='surface_upward_mass_flux_of_carbon_dioxide_expressed_as_carbon_due_to_heterotrophic_respiration', &
       fill_missing=.TRUE.)
  id_csoil = register_tiled_diag_field ( CMOR_NAME, 'cSoil', axes(1:1),  &
       lnd%time, 'Carbon in Soil Pool', 'kg m-2', missing_value=-100.0, &
       standard_name='soil_mass_content_of_carbon', fill_missing=.TRUE.)
  call add_tiled_diag_field_alias ( id_csoil, CMOR_NAME, 'cSoilLut', axes(1:1),  &
       lnd%time, 'Carbon  In Soil Pool On Land Use Tiles', 'kg m-2', missing_value=-100.0, &
       standard_name='soil_mass_content_of_carbon', fill_missing=.FALSE.)
  id_csoilfast = register_tiled_diag_field ( CMOR_NAME, 'cSoilFast', axes(1:1),  &
       lnd%time, 'Carbon Mass in Fast Soil Pool', 'kg m-2', missing_value=-100.0, &
       standard_name='fast_soil_pool_mass_content_of_carbon', fill_missing=.TRUE.)
  id_csoilmedium = register_tiled_diag_field ( CMOR_NAME, 'cSoilMedium', axes(1:1),  &
       lnd%time, 'Carbon Mass in Medium Soil Pool', 'kg m-2', missing_value=-100.0, &
       standard_name='medium_soil_pool_mass_content_of_carbon', fill_missing=.TRUE.)
  id_csoilslow = register_tiled_diag_field ( CMOR_NAME, 'cSoilSlow', axes(1:1),  &
       lnd%time, 'Carbon Mass in Slow Soil Pool', 'kg m-2', missing_value=-100.0, &
       standard_name='slow_soil_pool_mass_content_of_carbon', fill_missing=.TRUE.)
  id_cSoilLevels = register_tiled_diag_field ( CMOR_NAME, 'cSoilLevels', axes,  lnd%time, &
       'Carbon mass in each model soil level (summed over all soil carbon pools in that level)', &
       'kg m-2', missing_value=-100.0, &
       standard_name='soil_mass_content_of_carbon', fill_missing=.TRUE.)
  id_cLitter = register_tiled_diag_field ( CMOR_NAME, 'cLitter', axes(1:1), &
       lnd%time, 'Carbon Mass in Litter Pool', 'kg m-2', &
       missing_value=-100.0, standard_name='litter_mass_content_of_carbon', &
       fill_missing=.TRUE.)
  call add_tiled_diag_field_alias ( id_cLitter, CMOR_NAME, 'cLitterLut', axes(1:1),  &
       lnd%time, 'carbon in above and belowground litter pools on land use tiles', &
       'kg m-2', missing_value=-100.0, &
       standard_name='litter_mass_content_of_carbon', fill_missing=.FALSE.)
  id_cLitterCwd = register_tiled_diag_field ( CMOR_NAME, 'cLitterCwd', axes(1:1), &
       lnd%time, 'Carbon Mass in Coarse Woody Debris', 'kg m-2', &
       missing_value=-100.0, standard_name='wood_debris_mass_content_of_carbon', &
       fill_missing=.TRUE.)
  id_cLitterLeaf = register_tiled_diag_field ( CMOR_NAME, 'cLitterLeaf', axes(1:1), &
       lnd%time, 'Carbon Mass in Leaf Debris', 'kg m-2', &
       missing_value=-100.0, standard_name='leaf_debris_mass_content_of_carbon', &
       fill_missing=.TRUE.)
  call add_tiled_diag_field_alias ( id_rh, CMOR_NAME, 'rhLut', [ id_ug ],  &
       lnd%time, 'Soil Heterotrophic Respiration On Land Use Tile', 'kg m-2 s-1', &
       standard_name='surface_upward_mass_flux_of_carbon_dioxide_expressed_as_carbon_due_to_heterotrophic_respiration', &
       fill_missing=.FALSE., missing_value=-100.0)

end subroutine

! constructors
!> @brief Create new (empty) soil carbon representation
!! @return Pointer to new soil carbon data structure
function soilc_SIMPLE_ctor(soil) result(ptr)
  class(soil_BGC_SIMPLE_t), pointer :: ptr
  type(soil_tile_type), intent(in) :: soil !< soil tile data

  allocate(ptr)
  allocate( &
      ptr%fast_soil_C(num_l), &
      ptr%slow_soil_C(num_l), &
      ptr%asoil_in   (num_l), &
      ptr%fsc_in     (num_l), &
      ptr%ssc_in     (num_l)  )

  ptr%litter_SIMPLE_C(:,:)   = 0.0
  ptr%fast_soil_C(:)         = 0.0
  ptr%slow_soil_C(:)         = 0.0
  ptr%asoil_in(:)            = 0.0
  ptr%fsc_in(:)              = 0.0
  ptr%ssc_in(:)              = 0.0
end function

!> @brief Create a copy of existing soil carbon representation
!! @return Pointer to new soil carbon data structure
function soilc_SIMPLE_copy(soilc) result(ptr)
  type(soil_BGC_SIMPLE_t), pointer :: ptr
  type(soil_BGC_SIMPLE_t), intent(in) :: soilc !< soil carbon data to copy

  allocate(ptr)
  ptr = soilc
end function

!> @brief merge s1 into current soil carbon type s2, with given weights
subroutine merge_SIMPLE(s2,w2,s1,w1)
  class(soil_BGC_SIMPLE_t), intent(inout) :: s2    !< current soil carbon state
  class(soil_BGC_t)  , intent(in)    :: s1    !< soil carbon state to be merged into current
  real               , intent(in)    :: w2,w1 !< merging weights

  real    :: x1, x2 ! normalized relative weights

  ! calculate normalized weights
  x1 = w1/(w1+w2)
  x2 = 1.0 - x1

  select type(s1)
  type is (soil_BGC_SIMPLE_t)
     ! merge soil carbon
     s2%fast_soil_C(:) = s1%fast_soil_C(:)*x1 + s2%fast_soil_C(:)*x2
     s2%slow_soil_C(:) = s1%slow_soil_C(:)*x1 + s2%slow_soil_C(:)*x2
     s2%litter_SIMPLE_C(:,:) = s1%litter_SIMPLE_C(:,:)*x1 + s2%litter_SIMPLE_C(:,:)*x2

     s2%asoil_in(:)    = s1%asoil_in(:)*x1 + s2%asoil_in(:)*x2
     s2%fsc_in(:)      = s1%fsc_in(:)*x1 + s2%fsc_in(:)*x2
     s2%ssc_in(:)      = s1%ssc_in(:)*x1 + s2%ssc_in(:)*x2
  class default
     call land_error_message('merge_SIMPLE: attempt to merge incompatible soil carbon types', FATAL)
  end select
end subroutine


!> @brief Given soil carbon state, return total soil C
!! @return total soil carbon, kgC/m2
real function total_C_SIMPLE(soilc) result(tot_C)
  class(soil_BGC_SIMPLE_t), intent(in)  :: soilc !< soil carbon data structure
  tot_C = sum(soilc%fast_soil_C(:))+sum(soilc%slow_soil_C(:)) &
        + sum(soilc%litter_SIMPLE_C(:,:))
end function

!> @brief Given soil carbon state, return total soil nitrogen
!! @return total soil nitrogen, kgN/m2
real function total_N_SIMPLE(soilc) result(tot_N)
  class(soil_BGC_SIMPLE_t), intent(in)  :: soilc ! soil carbon data structure
  tot_N = 0.0
end function

!> @brief Given soil carbon state, return carbon amount of litter relevant for surface
!! resistance calculations in legacy treatment of soil surface resistance
subroutine rav_C_SIMPLE(soilc,fast_C,slow_C,dmic_C)
  class(soil_BGC_SIMPLE_t), intent(in)  :: soilc !< soil carbon data structure
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

!> @brief Given soil carbon state, returns total litter C per litter pool
subroutine get_littC_SIMPLE(soilc, values)
  class(soil_BGC_SIMPLE_t), intent(in)  :: soilc     !< soil carbon data structure
  real,                intent(out) :: values(:) !< total C in litter, by litter pool, kgC/m2

  integer :: i
  do i = 1, N_LITTER_POOLS
     values(i) = sum(soilc%litter_SIMPLE_C(:,i))
  enddo
end subroutine

!> @brief return zeros in 1D array
!! This subroutine is used to retrieve substances and values that are not present in
!! soil carbon models, e.g. DOC when the dissolved carbon is not implemented
subroutine get_zero_2D(soilC, values)
  class(soil_BGC_SIMPLE_t), intent(in)  :: soilc ! soil carbon data structure
  real,                intent(out) :: values(:,:) ! (N_C_TYPES, num_l) ! [kg C/m^2] dissolved organic carbon

  values(:,:) = 0.0
end subroutine

!> @brief return zeros in 1D array
!! This function is used to retrieve substances and values that are not present in
!! soil carbon models, e.g. nitrate or ammonium when the nitrogen dynamics is not
!! implemented
subroutine get_zero_1D(soilC, values)
  class(soil_BGC_SIMPLE_t), intent(in)  :: soilc !< soil carbon data structure (unused)
  real,                intent(out) :: values(:) !< returned values
  values(:) = 0.0
end subroutine

subroutine add_soil_carbon_SIMPLE(soilc, vegn, &
        leaf_litter_C, wood_litter_C, root_litter_C, &
        leaf_litter_N, wood_litter_N, root_litter_N  )
  class(soil_BGC_SIMPLE_t),   intent(inout) :: soilc
  type(vegn_tile_type), intent(inout) :: vegn
  real, intent(in), optional :: leaf_litter_C(:)   ! (N_C_TYPES)
  real, intent(in), optional :: wood_litter_C(:)   ! (N_C_TYPES)
  real, intent(in), optional :: root_litter_C(:,:) ! (num_l,N_C_TYPES)
  real, intent(in), optional :: leaf_litter_N(:)   ! (N_C_TYPES)
  real, intent(in), optional :: wood_litter_N(:)   ! (N_C_TYPES)
  real, intent(in), optional :: root_litter_N(:,:) ! (num_l,N_C_TYPES)

  ! TODO: check array sizes

  integer :: k
  real :: fsc, ssc
  real :: leaf_litt_C(N_C_TYPES)
  real :: wood_litt_C(N_C_TYPES)
  real :: root_litt_C(size(soilc%fast_soil_C),N_C_TYPES)

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

  ! CEL=cellulose (fast); LIG=lignin (slow); this function reasonably assumes
  ! that there are no microbes in litter

  if (bulk) then
     if (tau_cwlitt_transfer>0.or.tau_lflitt_transfer>0) then
        ! put litterfall in litter pools
        soilc%litter_SIMPLE_C(:,LITT_LEAF)  = soilc%litter_SIMPLE_C(:,LITT_LEAF)  + leaf_litt_C(:)
        soilc%litter_SIMPLE_C(:,LITT_CWOOD) = soilc%litter_SIMPLE_C(:,LITT_CWOOD) + wood_litt_C(:)
        fsc = sum(root_litt_C(:,C_FAST))
        ssc = sum(root_litt_C(:,C_SLOW))
     else
        ! add litterfall to soil carbon directly. This is mostly to preserve bitwise
        ! reproducibility with older code versions
        fsc = leaf_litt_C(C_FAST) + wood_litt_C(C_FAST) + sum(root_litt_C(:,C_FAST))
        ssc = leaf_litt_C(C_SLOW) + wood_litt_C(C_SLOW) + sum(root_litt_C(:,C_SLOW))
     endif
     soilc%fast_soil_C(1) = soilc%fast_soil_C(1) + fsc
     soilc%slow_soil_C(1) = soilc%slow_soil_C(1) + ssc
     ! for budget tracking
     soilc%fsc_in(1) = soilc%fsc_in(1) + fsc
     soilc%ssc_in(1) = soilc%ssc_in(1) + ssc
  else ! by-layer soil carbon model
     if (tau_cwlitt_transfer>0.or.tau_lflitt_transfer>0) then
        ! put litterfall in litter pools
        soilc%litter_SIMPLE_C(:,LITT_LEAF)  = soilc%litter_SIMPLE_C(:,LITT_LEAF)  + leaf_litt_C(:)
        soilc%litter_SIMPLE_C(:,LITT_CWOOD) = soilc%litter_SIMPLE_C(:,LITT_CWOOD) + wood_litt_C(:)
        fsc = 0.0; ssc = 0.0
     else
        ! add litterfall to soil carbon directly. This is mostly to preserve bitwise
        ! reproducibility with older code versions
        fsc = leaf_litt_C(C_FAST) + wood_litt_C(C_FAST)
        ssc = leaf_litt_C(C_SLOW) + wood_litt_C(C_SLOW)
     endif
     soilc%fast_soil_C(1) = soilc%fast_soil_C(1) + fsc
     soilc%slow_soil_C(1) = soilc%slow_soil_C(1) + ssc
     ! for budget tracking
     soilc%fsc_in(1) = soilc%fsc_in(1) + fsc
     soilc%ssc_in(1) = soilc%ssc_in(1) + ssc
     do k = 1,size(soilc%fast_soil_C)
        soilc%fast_soil_C(k) = soilc%fast_soil_C(k) + root_litt_C(k,C_FAST)
        soilc%slow_soil_C(k) = soilc%slow_soil_C(k) + root_litt_C(k,C_SLOW)
        ! for budget tracking
        soilc%fsc_in(k) = soilc%fsc_in(k) + root_litt_C(k,C_FAST)
        soilc%ssc_in(k) = soilc%ssc_in(k) + root_litt_C(k,C_SLOW)
     enddo
  endif

  ! accumulate litterfall diagnostics: it is sent to diag and then reset at every time step
  vegn%litterfall_C(:,LITT_LEAF)  = vegn%litterfall_C(:,LITT_LEAF)  + leaf_litt_C(:)
  vegn%litterfall_C(:,LITT_CWOOD) = vegn%litterfall_C(:,LITT_CWOOD) + wood_litt_C(:)

end subroutine add_soil_carbon_SIMPLE

!> @brief Add root exudates to vertical profile
subroutine add_root_exudates_SIMPLE(soilc, exudateC, exudateN, ammonium, nitrate)
  class(soil_BGC_SIMPLE_t), intent(inout) :: soilc !< soil carbon data
  real, intent(in)           :: exudateC(:) !< (num_l) amount of C in exudate, kgC/m2 per layer
  real, intent(in), optional :: exudateN(:) !< (num_l) amount of N in exudate, kgN/m2 per layer
  real, intent(in), optional :: ammonium(:) !< (num_l) amount of ammonium in exudate, kgN/m2(?) per layer
  real, intent(in), optional :: nitrate (:) !< (num_l) amount of  nitrate in exudate, kgN/m2(?) per layer

  ! NOTE: nitrogen-related exudates are ignored in current implementation of LM3-like soil carbon

  integer :: k  ! iterator across layers
  real    :: fsc

  if (bulk) then
     fsc = sum(exudateC(:))
     soilc%fast_soil_C(1) = soilc%fast_soil_C(1) + fsc
     soilc%fsc_in(1)      = soilc%fsc_in(1)      + fsc ! for soil carbon equilibration
  else
     do k = 1, size(soilc%fast_soil_C(:))
        soilc%fast_soil_C(k) = soilc%fast_soil_C(k) + exudateC(k)
        soilc%fsc_in(k)      = soilc%fsc_in(k)      + exudateC(k) ! for soil carbon equilibration
     enddo
  endif
end subroutine

!> @brief Burn given fraction of litter and return amounts of burned C and N
subroutine burn_litter_frac_SIMPLE(soilc, frac, burned_C, burned_N)
  class(soil_BGC_SIMPLE_t), intent(inout) :: soilc !< soil carbon data structure
  real, intent(in)  :: frac(:)            !< fraction of litter to burn [0,1]
  real, intent(out) :: burned_C, burned_N !< amounts of burned carbon and nitrogen

  integer :: i

  burned_C = 0.0; burned_N = 0.0
  do i = 1,N_LITTER_POOLS
     burned_C = burned_C + frac(i)*sum(soilc%litter_SIMPLE_C(:,i))
     soilc%litter_SIMPLE_C(:,i) = (1-frac(i))*soilc%litter_SIMPLE_C(:,i)
  enddo
end subroutine

subroutine tracer_leaching_SIMPLE(soilC, diag, &
     wl, flow, div, &
     div_hlsp_DOC, div_hlsp_DON, div_hlsp_NO3, div_hlsp_NH4, &
     ! output
     total_DOC_div, total_DON_div, total_NO3_div, total_NH4_div )
  class(soil_BGC_SIMPLE_t),  intent(inout) :: soilC
  type(diag_buff_type), intent(inout) :: diag

  real, intent(in) :: flow(:), div(:), wl(:) ! flow (into layer) and wl in units of mm, downward is >0  !!!xz check the unit of dz (should be m in this subroutine), flow (shoul be mm)
  real, intent(in) :: div_hlsp_DOC(:,:) ! (N_C_TYPES, num_l) [kg C/m^2/s] net divergence loss from tile calculated in hlsp_hydrology
  real, intent(in) :: div_hlsp_DON(:,:) ! (N_C_TYPES, num_l) [kg N/m^2/s] net divergence
  real, intent(in) :: div_hlsp_NO3(:),div_hlsp_NH4(:) ! (num_l) [kg N/m^2/s] net divergence loss from tile calculated in hlsp_hydrology

  real, intent(out) :: total_DOC_div, total_DON_div, total_NO3_div, total_NH4_div

  total_DOC_div=0.0; total_DON_div=0.0; total_NO3_div=0.0; total_NH4_div=0.0
end subroutine

! ============================================================================
! Deposition of nitrogen
subroutine deposit_N_SIMPLE(soilc, NH4, NO3, N_org)
  class(soil_BGC_SIMPLE_t), intent(inout) :: soilc
  real, intent(in) :: NH4, NO3, N_org ! amounts of NH4, NO3, and organic nitrogen to deposit, kg N/m2
  ! do nothing now: nitrogen deposition is ignored
end subroutine

! ============================================================================
! Nitrogen uptake from the rhizosphere by roots (active transport across root-soil interface)
subroutine active_root_N_uptake_SIMPLE(soilc, vegn, N_uptake, dt, update_pools)
  class(soil_BGC_SIMPLE_t), intent(inout) :: soilc
  type(vegn_tile_type), intent(in)    :: vegn
  real,    intent(out) :: N_uptake(:) ! Nitrogen uptake, kg N per individual
  real,    intent(in)  :: dt ! in years
  logical, intent(in)  :: update_pools

  N_uptake       = 0.0
end subroutine

! ============================================================================
! Uptake of mineral N by mycorrhizal "scavengers" -- Should correspond to Arbuscular mycorrhizae
subroutine myc_scavenger_N_uptake_SIMPLE(soilc,vegn,N_uptake_cohorts,myc_efficiency,dt,update_pools)
  class(soil_BGC_SIMPLE_t),  intent(inout) :: soilc
  type(vegn_tile_type), intent(in) :: vegn
  real,intent(out) :: N_uptake_cohorts(:) ! Units: kgN/m2 per individual
  real, intent(in) :: dt  ! dt in years
  logical, intent(in) :: update_pools
  real, intent(out) :: myc_efficiency ! units: kgN/kg myc biomass C. Should give N uptake efficiency even when myc biomass is zero

  N_uptake_cohorts = 0.0
  myc_efficiency   = 0.0
end subroutine

! ============================================================================
! Uptake of mineral N by mycorrhizal "miners" -- Should correspond to Ecto mycorrhizae
subroutine myc_miner_N_uptake_SIMPLE(soilc,soil,vegn,N_uptake_cohorts,C_uptake_cohorts,total_CO2prod,myc_efficiency,dt,update_pools)
  class(soil_BGC_SIMPLE_t), intent(inout) :: soilc
  type(soil_tile_type), intent(in) :: soil
  type(vegn_tile_type), intent(in) :: vegn
  real,    intent(out) :: N_uptake_cohorts(:), C_uptake_cohorts(:)  ! Units kg/m2 of per individual
  real,    intent(out) :: total_CO2prod ! Units of kgC/m2 (not per individual)
  real,    intent(in)  :: dt  ! dt in years
  logical, intent(in)  :: update_pools
  real,    intent(out) :: myc_efficiency  ! units: kgN/kg myc biomass C. Should give N uptake efficiency even when myc biomass is zero
end subroutine

! ============================================================================
! transfer fraction of intermediate pools (used for smoothing out contributions
! of spiky processes, e.g. harvesting) defined in vegetation tile data structure
! to soil BGC pools.
subroutine spend_intermediate_pools_SIMPLE(soilc, vegn)
  class(soil_BGC_SIMPLE_t), intent(inout) :: soilc
  type(vegn_tile_type) , intent(inout) :: vegn

  integer :: i,k
  ! move carbon from intermediate spike-process buffers to litter
  do i = 1,N_C_TYPES
     do k = 1, N_LITTER_POOLS
        call deplete_pool(vegn%litter_buff_C(i,k), vegn%litter_rate_C(i,k), soilc%litter_SIMPLE_C(i,k),vegn%litterfall_C(i,k))
     enddo
  enddo

  call deplete_pool(vegn%fsc_pool_bg, vegn%fsc_rate_bg, soilc%fast_soil_C(1), soilc%fsc_in(1))
  call deplete_pool(vegn%ssc_pool_bg, vegn%ssc_rate_bg, soilc%slow_soil_C(1), soilc%ssc_in(1))

  ! transfer litter to soilc pools, with constant time scales
  call deplete_pool1(soilc%litter_SIMPLE_C(C_FAST, LITT_LEAF),  tau_lflitt_transfer, soilc%fast_soil_C(1), soilc%fsc_in(1))
  call deplete_pool1(soilc%litter_SIMPLE_C(C_MIC,  LITT_LEAF),  tau_lflitt_transfer, soilc%fast_soil_C(1), soilc%fsc_in(1))
  call deplete_pool1(soilc%litter_SIMPLE_C(C_SLOW, LITT_LEAF),  tau_lflitt_transfer, soilc%slow_soil_C(1), soilc%ssc_in(1))

  call deplete_pool1(soilc%litter_SIMPLE_C(C_FAST, LITT_CWOOD), tau_cwlitt_transfer, soilc%fast_soil_C(1), soilc%fsc_in(1))
  call deplete_pool1(soilc%litter_SIMPLE_C(C_MIC,  LITT_CWOOD), tau_cwlitt_transfer, soilc%fast_soil_C(1), soilc%fsc_in(1))
  call deplete_pool1(soilc%litter_SIMPLE_C(C_SLOW, LITT_CWOOD), tau_cwlitt_transfer, soilc%slow_soil_C(1), soilc%ssc_in(1))
end subroutine

! ============================================================================
!> @brief Move substance from one pool to another
!!
!! Given an intermediate pool of matter (e.g. C or N), and its e-folding time
!! scale, move the amount of mass corresponding to one fats time step from the
!! pool to the destination.
!!
!! @note
!! In contrast to "deplete_pool" subroutime it accepts the time scale tau,
!! instead of spending rate
subroutine deplete_pool1(pool, tau, dest, accum)
   real, intent(inout) :: pool !< C or N intermediate pool, kg
   real, intent(in)    :: tau  !< C or N e-folding time scale, years
   real, intent(inout) :: dest !< C or N destination pool, kg
   real, intent(inout), optional :: accum !< accumulator for soil carbon equilibration, e.g. fs_in or ssc_in

   real :: rate ! rate of depletion, kgC/m2/year

   if (tau > 0) then
      rate = pool/tau
   else
      rate = pool/dt_fast_yr
   endif
   call deplete_pool(pool, rate, dest, accum)
end subroutine deplete_pool1

subroutine dsdt_SIMPLE(soilc, soil, vegn, diag, soilt, theta)
  class(soil_BGC_SIMPLE_t)  , intent(inout) :: soilc
  type(vegn_tile_type), intent(inout) :: vegn
  type(soil_tile_type), intent(inout) :: soil
  type(diag_buff_type), intent(inout) :: diag
  real                , intent(in)    :: soilt ! average soil temperature, deg K
  real                , intent(in)    :: theta ! average soil moisture

  real :: fast_C_loss(size(soilc%fast_soil_C))
  real :: slow_C_loss(size(soilc%slow_soil_C))
  real :: A          (size(soilc%slow_soil_C)) ! decomp rate reduction due to moisture and temperature

  if (bulk) then
      A(:) = A_function(soilt, theta)
  else
      A(:) = A_function(soil%T, soil_theta(soil))
  endif

  fast_C_loss = soilc%fast_soil_C(:)*A*K1*dt_fast_yr;
  slow_C_loss = soilc%slow_soil_C(:)*A*K2*dt_fast_yr;

  soilc%fast_soil_C = soilc%fast_soil_C - fast_C_loss;
  soilc%slow_soil_C = soilc%slow_soil_C - slow_C_loss;

  ! for budget check
  vegn%fsc_out = vegn%fsc_out + sum(fast_C_loss(:));
  vegn%ssc_out = vegn%ssc_out + sum(slow_C_loss(:));

  ! loss of C to atmosphere and leaching
  vegn%rh = sum(fast_C_loss(:)+slow_C_loss(:))/dt_fast_yr;

  ! accumulate decomposition rate reduction for the soil carbon restart output
  soilc%asoil_in(:) = soilc%asoil_in(:) + A(:)

  ! ---- diagnostic section
  call send_tile_data(id_rsoil_C(C_FAST), fast_C_loss(:)/(dz(1:num_l)*dt_fast_yr), diag)
  call send_tile_data(id_rsoil_C(C_SLOW), slow_C_loss(:)/(dz(1:num_l)*dt_fast_yr), diag)
  call send_tile_data(id_rsoil, vegn%rh, diag)

  ! TODO: arithmetic averaging of A does not seem correct; we need to invent something better,
  !       e.g. weight it with the carbon loss, or something like that
  if (id_asoil>0) call send_tile_data(id_asoil, sum(A(:))/size(A(:)), diag)
  call send_tile_data(id_rh, vegn%rh/seconds_per_year, diag)
end subroutine dsdt_SIMPLE

! ============================================================================
subroutine step3_SIMPLE(soilc, diag)
  class(soil_BGC_SIMPLE_t),   intent(inout) :: soilc
  type(diag_buff_type), intent(inout) :: diag

  integer :: i, k

  associate (soil=>soilc) ! to avoid renaming
  call send_tile_data(id_fsc, sum(soil%fast_soil_C(:))+sum(soil%litter_SIMPLE_C(C_FAST,:)), diag)
  call send_tile_data(id_ssc, sum(soil%slow_soil_C(:))+sum(soil%litter_SIMPLE_C(C_SLOW,:)), diag)
  call send_tile_data(id_soil_C(C_FAST), soil%fast_soil_C(:)/dz(1:num_l), diag)
  call send_tile_data(id_soil_C(C_SLOW), soil%slow_soil_C(:)/dz(1:num_l), diag)
  call send_tile_data(id_total_soil_C, sum(soil%fast_soil_C(:))+sum(soil%slow_soil_C(:))+sum(soil%litter_SIMPLE_C(:,:)), diag)
  do k = 1, N_LITTER_POOLS
     if (id_litter_total_C(k)>0) call send_tile_data(id_litter_total_C(k), sum(soil%litter_SIMPLE_C(:,k)), diag)
     do i = 1, N_C_TYPES
        call send_tile_data(id_litter_C(k,i), soil%litter_SIMPLE_C(i,k), diag)
     enddo
  enddo

  ! --- CMOR vars
  if (id_csoilfast>0)   call send_tile_data(id_csoilfast,   sum(soil%fast_soil_C(:)), diag)
  if (id_csoilmedium>0) call send_tile_data(id_csoilmedium, sum(soil%slow_soil_C(:)), diag)
  call send_tile_data(id_csoilslow, 0.0, diag)
  if (id_csoil>0)       call send_tile_data(id_csoil, sum(soil%fast_soil_C(:))+sum(soil%slow_soil_C(:)), diag)
  if (id_cSoilLevels>0) call send_tile_data(id_cSoilLevels, soil%fast_soil_C(:)+soil%slow_soil_C(:), diag)
  if (id_cLitter>0)     call send_tile_data(id_cLitter, sum(soil%litter_SIMPLE_C(:,:)), diag)
  if (id_cLitterCwd>0)  call send_tile_data(id_cLitterCwd, sum(soil%litter_SIMPLE_C(:,LITT_CWOOD)), diag)
  if (id_cLitterLeaf>0) call send_tile_data(id_cLitterLeaf, sum(soil%litter_SIMPLE_C(:,LITT_LEAF)), diag)
  ! --- end of CMOR vars
  end associate

end subroutine

! ============================================================================
! The combined reduction in decomposition rate as a funciton of TEMP and MOIST
! Based on CENTURY Parton et al 1993 GBC 7(4):785-809 and Bolker's copy of
! CENTURY code
elemental function A_function(soilt, theta) result(A)
  real :: A                 ! return value, resulting reduction in decomposition rate
  real, intent(in) :: soilt ! effective temperature for soil carbon decomposition
  real, intent(in) :: theta

  real :: soil_temp; ! temperature of the soil, deg C
  real :: Td; ! rate multiplier due to temp
  real :: Wd; ! rate reduction due to mositure

  ! coefficeints and terms used in temperaturex term
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
subroutine redistribute_peat_carbon_SIMPLE(soilc)
  class(soil_BGC_SIMPLE_t), intent(inout) :: soilc

  call error_mesg('redistribute_peat_carbon_SIMPLE','not implemented; should it be?', FATAL)
end subroutine


end module
