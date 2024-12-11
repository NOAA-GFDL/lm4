module gl_snow_tile_mod
#include <fms_platform.h>
#include "../../shared/debug.inc"

use mpp_mod, only: input_nml_file

use fms_mod, only : FATAL
use constants_mod, only : tfreeze, hlf

use land_constants_mod, only : NBANDS
use land_tile_selectors_mod, only : tile_selector_type
use land_data_mod, only : log_version, lnd
use land_debug_mod, only : is_watch_point, is_watch_cell, land_error_message

use tile_diag_buff_mod, only : diag_buff_type
use tile_diag_base_mod, only : set_default_diag_filter, &
        register_tiled_diag_field, send_tile_data

use snowpack_mod, only : snow_layer_type, snowpack_t
use snow_tile_mod, only: snow_tile_type, NTRACERS, z0_momentum, &
     k_over_B, cpw, clw, csw, snow_data_area, snow_lw_properties
use snow_evolution_mod, only : gl_sweep_tiny_snow, assign_substrate_sw_to_surface, &
     albedo_option, ALBEDO_SNICAR, use_internal_sources, thresh_snow_depth_swheat, gl_compute_snow_albedo, &
     gl_snow_step_2_ev => gl_snow_step_2, delta_time, do_mgimplicit, gl_sweep_huge_snow


implicit none
private

! ==== public interfaces =====================================================
public :: gl_snow_tile_type
public :: new_gl_snow_tile
public :: gl_snow_diag_init


! ==== module constants ======================================================
character(*), parameter :: module_name = 'gl_snow_tile_mod'
#include "../../shared/version_variable.inc"

! ---- module interfaces
interface new_gl_snow_tile
   module procedure gl_snow_tile_ctor
   module procedure gl_snow_tile_copy
end interface

! ==== types =================================================================
type, extends(snow_tile_type) :: gl_snow_tile_type
    type(snowpack_t) :: sp ! structure with data for glass snow model
contains
    procedure :: merge_snow_tiles => gl_merge_snow_tiles_wrapper

    procedure :: n_layers => gl_snow_nlayers
    procedure :: snow_is_selected => gl_snow_is_selected
    procedure :: snow_roughness => gl_snow_roughness
    procedure :: radiative_properties => gl_snow_rad_prop
    procedure :: snow_active => gl_snow_active
    procedure :: snow_tile_heat => gl_snow_tile_heat
    procedure :: sfc_temp => gl_snow_get_sfc_temp


    procedure :: get_Ti => gl_snow_get_Ti
    procedure :: get_wli => gl_snow_get_wli
    procedure :: get_wsi => gl_snow_get_wsi

    procedure :: set_Ti =>  gl_snow_set_Ti
    procedure :: set_wli => gl_snow_set_wli
    procedure :: set_wsi => gl_snow_set_wsi

    procedure :: ice => gl_snow_get_total_ice
    procedure :: liq => gl_snow_get_total_liq
    procedure :: get_depth_area => gl_snow_get_depth_area
    procedure :: lai_im => gl_snow_lai_im
    procedure :: lai_em => gl_snow_lai_em

    procedure :: sweep => gl_sweep_snow
    procedure :: partition_sw => gl_partition_sw

    procedure :: step1 => gl_snow_step_1
    procedure :: step2 => gl_snow_step_2
    procedure :: send_diag => gl_snow_send_diag
end type gl_snow_tile_type

! diagnostic field IDs
integer :: id_snow_avrg_optd, id_snow_avrg_sph, id_snow_avrg_dendr, id_snow_density, &
    id_snow_avrg_age, id_snow_nearsurf_bceq_tot, &
    id_snow_nearsurf_bceq_im, id_snow_nearsurf_bceq_em, id_snow_avrg_bceq_tot, &
    id_snow_avrg_bc_tot, id_snow_avrg_md_tot, id_snow_avrg_om_tot, &
    id_snow_avrg_bceq_im, id_snow_avrg_bceq_em, id_snow_nearsurf_optd, &
    id_snow_nearsurf_sph, id_snow_nearsurf_dendr, id_snow_nearsurf_age, &
    id_snow_nearsurf_density, id_snow_liq, id_snow_ice, &
    id_snow_topwater, id_snow_topsnowdeficit, id_snow_topwheat, &
    id_snow_topsnowheatdeficit

contains ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

! Register diagnostic fields
subroutine gl_snow_diag_init(id_ug)
  integer,intent(in)  :: id_ug    !< Unstructured axis id

  character(*), parameter :: diag_mod_name = 'land' ! name of the component used for diagnostic fields


  call log_version(version, module_name, &
  __FILE__)

  ! set the default sub-sampling filter for the fields below
  call set_default_diag_filter('land')

  ! ------------------------------------ EZSNOW new snowpack model added fields ----------
  ! // TODO fix missing values, and add fix for non-extensive variables [e.g., snow grain properties]
  id_snow_avrg_optd = register_tiled_diag_field ( diag_mod_name, 'snow_avrg_optd', (/id_ug/), lnd%time, &
     'Snowpack average optical diameter', 'm', missing_value=-9999.0) !
  id_snow_avrg_sph = register_tiled_diag_field ( diag_mod_name, 'snow_avrg_sph', (/id_ug/), lnd%time, &
     'Snowpack average sphericity', 'dimless', missing_value=-9999.0) !
  id_snow_avrg_dendr = register_tiled_diag_field ( diag_mod_name, 'snow_avrg_dendr', (/id_ug/), lnd%time, &
     'Snowpack average dendricity', 'dimless', missing_value=-9999.0) !
  id_snow_density = register_tiled_diag_field ( diag_mod_name, 'snow_density', (/id_ug/), lnd%time, &
     'Snowpack density', 'kg/m3', missing_value=-9999.0)
  id_snow_avrg_age = register_tiled_diag_field ( diag_mod_name, 'snow_avrg_age', (/id_ug/), lnd%time, &
     'Snowpack average age', 'days', missing_value=-9999.0) !
!   id_snow_avrg_T = register_tiled_diag_field ( diag_mod_name, 'snow_avrg_T', (/id_ug/), lnd%time, &
!      'Snowpack average temperature', 'degK', missing_value=-9999.0) !

  ! TODO: add axis = 3 impurities
  ! id_snow_lai_im = register_tiled_diag_field ( diag_mod_name, 'snow_lai_im', (/id_ug/), lnd%time, &
     ! 'Snowpack content of internally mixed light-absorbing impurities', 'ppm', missing_value=-9999.0)
  ! id_snow_lai_em = register_tiled_diag_field ( diag_mod_name, 'snow_lai_em', (/id_ug/), lnd%time, &
     ! 'Snowpack content of externally mixed light-absorbing impurities', 'ppm', missing_value=-9999.0)

  id_snow_nearsurf_bceq_tot = register_tiled_diag_field ( diag_mod_name, 'snow_nearsurf_bceq_tot', (/id_ug/), lnd%time, &
     'Snowpack total (im + em) near-surface conc. of light-absorbing impurities', 'ppm', missing_value=-9999.0)
  id_snow_nearsurf_bceq_im = register_tiled_diag_field ( diag_mod_name, 'snow_nearsurf_bceq_im', (/id_ug/), lnd%time, &
     'Snowpack near-surface conc. of internally mixed light-absorbing impurities', 'ppm', missing_value=-9999.0)
  id_snow_nearsurf_bceq_em = register_tiled_diag_field ( diag_mod_name, 'snow_nearsurf_bceq_em', (/id_ug/), lnd%time, &
     'Snowpack near-surface conc. of externally mixed light-absorbing impurities', 'ppm', missing_value=-9999.0)
  id_snow_avrg_bceq_tot = register_tiled_diag_field ( diag_mod_name, 'snow_avrg_bceq_tot', (/id_ug/), lnd%time, &
     'Snowpack total (im + em) average conc. of light-absorbing impurities', 'ppm', missing_value=-9999.0)
  id_snow_avrg_bc_tot = register_tiled_diag_field ( diag_mod_name, 'snow_avrg_bc_tot', (/id_ug/), lnd%time, &
     'Snowpack total (im + em) average conc. of black carbon', 'ppm', missing_value=-9999.0)
  id_snow_avrg_md_tot = register_tiled_diag_field ( diag_mod_name, 'snow_avrg_md_tot', (/id_ug/), lnd%time, &
     'Snowpack total (im + em) average conc. of mineral dust', 'ppm', missing_value=-9999.0)
  id_snow_avrg_om_tot = register_tiled_diag_field ( diag_mod_name, 'snow_avrg_om_tot', (/id_ug/), lnd%time, &
     'Snowpack total (im + em) average conc. of organic carbon', 'ppm', missing_value=-9999.0)
  id_snow_avrg_bceq_im = register_tiled_diag_field ( diag_mod_name, 'snow_avrg_bceq_im', (/id_ug/), lnd%time, &
     'Snowpack average conc. of internally mixed light-absorbing impurities', 'ppm', missing_value=-9999.0)
  id_snow_avrg_bceq_em = register_tiled_diag_field ( diag_mod_name, 'snow_avrg_bceq_em', (/id_ug/), lnd%time, &
     'Snowpack average conc. of externally mixed light-absorbing impurities', 'ppm', missing_value=-9999.0)

  id_snow_nearsurf_optd = register_tiled_diag_field ( diag_mod_name, 'snow_nearsurf_optd', (/id_ug/), lnd%time, &
     'Snowpack near-surface optical diameter', 'm', missing_value=-9999.0)
  id_snow_nearsurf_sph = register_tiled_diag_field ( diag_mod_name, 'snow_nearsurf_sph', (/id_ug/), lnd%time, &
     'Snowpack near-surface grain sphericity', 'dimless', missing_value=-9999.0)
  id_snow_nearsurf_dendr = register_tiled_diag_field ( diag_mod_name, 'snow_nearsurf_dendr', (/id_ug/), lnd%time, &
     'Snowpack near-surface grain dendricity', 'dimless', missing_value=-9999.0) !
  id_snow_nearsurf_age = register_tiled_diag_field ( diag_mod_name, 'snow_nearsurf_age', (/id_ug/), lnd%time, &
     'Snowpack near-surface age', 'days', missing_value=-9999.0) !
  id_snow_nearsurf_density = register_tiled_diag_field ( diag_mod_name, 'snow_nearsurf_density', (/id_ug/), lnd%time, &
     'Snowpack near-surface density', 'kg/m3', missing_value=-9999.0)

  ! id_snow_area_frac = register_tiled_diag_field ( diag_mod_name, 'snow_area_frac', (/id_ug/), lnd%time, &
     ! 'Frcational snow-covered area', 'dimless', missing_value=-9999.0)
!   id_snow_depth = register_tiled_diag_field ( diag_mod_name, 'snow_depth', (/id_ug/), lnd%time, &
!      'Snow depth', 'm', missing_value=-9999.0)
  id_snow_liq = register_tiled_diag_field ( diag_mod_name, 'snow_liq', (/id_ug/), lnd%time, &
     'Snowpack total liquid', 'kg/m2', missing_value=-9999.0)
  id_snow_ice = register_tiled_diag_field ( diag_mod_name, 'snow_ice', (/id_ug/), lnd%time, &
     'Snowpack total ice', 'kg/m2', missing_value=-9999.0)

  id_snow_topwater = register_tiled_diag_field ( diag_mod_name, 'snow_topwater', (/id_ug/), lnd%time, &
     'Snowpack topwater', 'kg/m2', missing_value=-1.0e+20)
  id_snow_topsnowdeficit = register_tiled_diag_field ( diag_mod_name, 'snow_topsnowdeficit', (/id_ug/), lnd%time, &
     'Snowpack topsnowdeficit', 'kg/m2', missing_value=-1.0e+20)
  id_snow_topwheat = register_tiled_diag_field ( diag_mod_name, 'snow_topwheat', (/id_ug/), lnd%time, &
     'Snowpack topwheat', 'J/m2', missing_value=-1.0e+20)
  id_snow_topsnowheatdeficit = register_tiled_diag_field ( diag_mod_name, 'snow_topsnowheatdeficit', (/id_ug/), lnd%time, &
     'Snowpack topsnowheatdeficit', 'J/m2', missing_value=-1.0e+20)

end subroutine gl_snow_diag_init


! ============================================================================
function gl_snow_tile_ctor(tag) result(ptr)
  type(gl_snow_tile_type), pointer :: ptr ! return value
  integer, optional, intent(in) :: tag ! kind of tile

  allocate(ptr)
  ptr%tag = 0 ; if(present(tag)) ptr%tag = tag

end function gl_snow_tile_ctor

! ============================================================================
function gl_snow_tile_copy(snow) result(ptr)
  type(gl_snow_tile_type), pointer :: ptr ! return value
  type(gl_snow_tile_type), intent(in) :: snow ! tile to copy

  allocate(ptr)
  ! copy all non-pointer members
  ptr = snow
  ! no need to allocate storage for allocatable components of the type, because
  ! F2003 takes care of that, and also takes care of copying data
end function

! ============================================================================
subroutine gl_delete_snow_tile(snow)
  type(gl_snow_tile_type), pointer :: snow

  ! no need to deallocate components of tile, because F2003 takes care of
  ! allocatable components deallocation when tile is deallocated
  ! deallocate(snow)
  nullify(snow) ! slm: not sure why this is done. soes it not prevent tile from cleaning up?
end subroutine gl_delete_snow_tile

! =============================================================================
function gl_snow_tiles_can_be_merged(snow1,snow2) result(response)
  logical :: response
  type(gl_snow_tile_type), intent(in) :: snow1,snow2
  response = .TRUE.
end function gl_snow_tiles_can_be_merged


! subroutine cm1_merge_snow_tiles(snow1, w1, snow2, w2)
subroutine gl_merge_snow_tiles_wrapper(snow2, w2, snow1, w1)
  class(snow_tile_type), intent(in)    :: snow1
  class(gl_snow_tile_type), intent(inout) :: snow2
  real                , intent(in)    :: w1, w2 ! relative weights
  select type (snow1)
  type is (gl_snow_tile_type)
      call gl_merge_snow_tiles(snow2, w2, snow1, w1)
  class default
  call land_error_message('gl_merge_snow_tiles_wrapper in gl_snow_tile_mod: type is incorrect!', FATAL)
  end select
end subroutine gl_merge_snow_tiles_wrapper


! =============================================================================
subroutine gl_merge_snow_tiles(snow2ez, w2, snow1ez, w1)
  type(gl_snow_tile_type), intent(in)    :: snow1ez
  type(gl_snow_tile_type), intent(inout) :: snow2ez
  type(snowpack_t)  :: snow1
  type(snowpack_t)  :: snow2
  real                , intent(in)    :: w1, w2 ! relative weights


  real wat1, wat2, wat3, dwat, heat1, heat2, heat3, dheat, lai1, lai2, lai3, dlai


  ! ---- local vars
  real    :: x1, x2 ! normalized weights
  ! real    :: HEAT1, HEAT2
  integer :: il, it
  integer :: i1, i2, count
  integer new_nlayers, minc_nlayers
  ! type(snow_layer_type), ALLOCATABLE :: snowt(:)
  type(snowpack_t) :: snow3

  real snow_mass
  real total_lai_loss, total_heat_loss, total_water_loss
  ! below this value simplify the mergining procedure for numerical stability - case of very little snow
  real, parameter :: min_snow_mass_merging=1E-8
  real total_heat, total_mass
  logical, allocatable :: layers2remove(:)
  integer final_nlayers
  ! remove/merge single snow layers below this soild ice content
  real,parameter::thresh_ws=0.001
  real min_ws
  integer current_nlayers
  type(snow_layer_type), allocatable :: snowF(:)
  logical merged
  real, dimension(3) :: temp_lost_im, temp_lost_em
  real, dimension(3) :: final_lost_im, final_lost_em
  real temp_lost_wl, temp_lost_ws, temp_lost_heat, temp_density, temp_lost_dz
  real temp_lost_age_w, temp_lost_sph_w, temp_lost_optd_w, temp_lost_dendr_w
  real true_heat_3, fict_heat_3, orig_heat_1, orig_heat_2
  integer ill
  logical use_first_tile


  snow1 = snow1ez%sp
  snow2 = snow2ez%sp

  lai1 = sum(snow1ez%sp%lai_im() + snow1ez%sp%lai_em())
  lai2 = sum(snow2ez%sp%lai_im() + snow2ez%sp%lai_em())
  wat1 = snow1ez%sp%SWE()
  wat2 = snow2ez%sp%SWE()
  heat1 = snow1ez%sp%heat()
  heat2 = snow2ez%sp%heat()

  ! calculate normalized weights
  x1 = w1/(w1+w2)
  x2 = 1.0 - x1

  ! new_nlayers = snow1%nlayers + snow2%nlayers
  new_nlayers = max(snow1%nlayers, snow2%nlayers)
  minc_nlayers = min(snow1%nlayers, snow2%nlayers)

  snow3%nlayers = new_nlayers
  snow3%topwater = x1 * snow1%topwater + x2 * snow2%topwater
  snow3%topwheat = x1 * snow1%topwheat + x2 * snow2%topwheat
  snow3%topsnowdeficit = x1 * snow1%topsnowdeficit + x2 * snow2%topsnowdeficit
  snow3%topsnowheatdeficit = x1 * snow1%topsnowheatdeficit + x2 * snow2%topsnowheatdeficit

  snow3%beta_rad = x1 * snow1%beta_rad + x2 * snow2%beta_rad

    allocate(snow3%snow(new_nlayers)) ! ok also in case of zero size? move later?
    ! allocate(snowt(new_nlayers)) ! ok also in case of zero size? move later?
    allocate(layers2remove(new_nlayers))
    do il = 1, new_nlayers
      layers2remove(il) = .FALSE.
    enddo


  ! if allocated, deallocate s%e, s%f, s%swheat
  ! albedo not needed anymore? will be computed when needed for new tile
  ! nearsurf properties: not needed, will be computed when needed  for new tile
  if (new_nlayers == 0) then
      ! all done here, just deallocate snow array?
  ! case of only one snow with layers, the other, or both
  else if ((snow1%nlayers>0).and.(snow2%nlayers==0)) then
    snow3%snow = snow1%snow
    do il = 1, new_nlayers ! rescale quantities from tile 1 to new tile area
      snow3%snow(il)%dz = snow3%snow(il)%dz * x1
      snow3%snow(il)%ws = snow3%snow(il)%ws * x1
      snow3%snow(il)%wl = snow3%snow(il)%wl * x1
      do it = 1, NTRACERS ! rescale quantities from tile 1 to new tile area
        snow3%snow(il)%wc_em(it) = snow3%snow(il)%wc_em(it) * x1
        snow3%snow(il)%wc_im(it) = snow3%snow(il)%wc_im(it) * x1
      enddo
    enddo


  else if ((snow1%nlayers==0).and.(snow2%nlayers>0)) then
    snow3%snow = snow2%snow
    do il = 1, new_nlayers ! rescale quantities from tile 1 to new tile area
      snow3%snow(il)%dz = snow3%snow(il)%dz * x2
      snow3%snow(il)%ws = snow3%snow(il)%ws * x2
      snow3%snow(il)%wl = snow3%snow(il)%wl * x2
      do it = 1, NTRACERS ! rescale quantities from tile 1 to new tile area
        snow3%snow(il)%wc_em(it) = snow3%snow(il)%wc_em(it) * x2
        snow3%snow(il)%wc_im(it) = snow3%snow(il)%wc_im(it) * x2
      enddo
    enddo

  else ! both snowpacks have non-zero number of layers

    if(is_watch_cell()) then
      write(*,*) '#### gl_snow_tile - merging tiles glass checkpoint 1 ####'
      write(*,*) "case of both snowpacks with non-zero snow layers"
      write(*,*) "fractions: w1, w2, x1, x2 = ", w1, w2, x1, x2
      ! call snow1%print()
      ! write(*,*) "Second snow object:"
      ! call snow2%print()
      ! write(*,*) "Resulting snow object:"
    endif
    ! stack layers starting from the top
    ! there is probably a better way to do this but this seems effective
    do il=1,new_nlayers
      ! both tiles haved this layer
      if (il .le. minc_nlayers) then
        snow3%snow(il)%dz = snow1%snow(il)%dz * x1 + snow2%snow(il)%dz * x2
        snow3%snow(il)%ws = snow1%snow(il)%ws * x1 + snow2%snow(il)%ws * x2
        snow3%snow(il)%wl = snow1%snow(il)%wl * x1 + snow2%snow(il)%wl * x2
        do it = 1, NTRACERS ! rescale quantities from tile 1 to new tile area
          snow3%snow(il)%wc_em(it) = snow1%snow(il)%wc_em(it) * x1 + snow2%snow(il)%wc_em(it) * x2
          snow3%snow(il)%wc_im(it) = snow1%snow(il)%wc_im(it) * x1 + snow2%snow(il)%wc_im(it) * x2
        enddo
        ! intensive quantities to be weight - averaged based on solid mass:
        snow3%snow(il)%age = ( snow1%snow(il)%age * snow1%snow(il)%ws * x1 + snow2%snow(il)%age * snow2%snow(il)%ws * x2) / (snow1%snow(il)%ws * x1 + snow2%snow(il)%ws * x2)
        snow3%snow(il)%sph = ( snow1%snow(il)%sph * snow1%snow(il)%ws * x1 + snow2%snow(il)%sph * snow2%snow(il)%ws * x2) / (snow1%snow(il)%ws * x1 + snow2%snow(il)%ws * x2)
        snow3%snow(il)%optd = ( snow1%snow(il)%optd * snow1%snow(il)%ws * x1 + snow2%snow(il)%optd * snow2%snow(il)%ws * x2) / (snow1%snow(il)%ws * x1 + snow2%snow(il)%ws * x2)
        snow3%snow(il)%dendr = ( snow1%snow(il)%dendr * snow1%snow(il)%ws * x1 + snow2%snow(il)%dendr * snow2%snow(il)%ws * x2) / (snow1%snow(il)%ws * x1 + snow2%snow(il)%ws * x2)
        ! heat balance of the two layers:
        ! do not change phases here -
        ! snow3%snow(il)%T = TFREEZE + ( (snow1%snow(il)%T - TFREEZE) * snow1%snow(il)%hCap() * x1 + (snow2%snow(il)%T-TFREEZE) * snow2%snow(il)%hCap() * x2) / (snow1%snow(il)%hCap() * x1 + snow2%snow(il)%hCap() * x2)
        !
        ! snow3%snow(il)%T = TFREEZE + ( (snow1%snow(il)%T - TFREEZE) * snow1%snow(il)%hCap() * x1 + (snow2%snow(il)%T-TFREEZE) * snow2%snow(il)%hCap() * x2) / &
        !            (snow1%snow(il)%hCap() * x1 + snow2%snow(il)%hCap() * x2 + snow1%snow(il)%wl * x1 * HLF + snow2%snow(il)%wl * x2 * HLF )
        !
                  !  snow3%snow(il)%T = TFREEZE + ( (snow1%snow(il)%T - TFREEZE) * snow1%snow(il)%hCap() * x1 + snow1%snow(il)%wl * x1 * HLF + (snow2%snow(il)%T-TFREEZE) * snow2%snow(il)%hCap() * x2 + snow2%snow(il)%wl * x2 * HLF ) / &
                  !  (snow1%snow(il)%hCap() * x1 + snow2%snow(il)%hCap() * x2 + snow1%snow(il)%wl * x1 * HLF + snow2%snow(il)%wl * x2 * HLF )

                   snow3%snow(il)%T = TFREEZE + ( (snow1%snow(il)%T - TFREEZE) * snow1%snow(il)%hCap() * x1 + (snow2%snow(il)%T-TFREEZE) * snow2%snow(il)%hCap() * x2 ) / &
                   (snow1%snow(il)%hCap() * x1 + snow2%snow(il)%hCap() * x2  )
        ! end case of shared layer index
      else
        ! only one of two tiles has this layer
        if (il .le. snow1%nlayers) then
          snow3%snow(il)    = snow1%snow(il)
          snow3%snow(il)%T    = snow1%snow(il)%T
          snow3%snow(il)%age    = snow1%snow(il)%age
          snow3%snow(il)%sph    = snow1%snow(il)%sph
          snow3%snow(il)%optd    = snow1%snow(il)%optd
          snow3%snow(il)%dendr    = snow1%snow(il)%dendr
          snow3%snow(il)%dz = snow3%snow(il)%dz * x1
          snow3%snow(il)%ws = snow3%snow(il)%ws * x1
          snow3%snow(il)%wl = snow3%snow(il)%wl * x1
          do it = 1, NTRACERS ! rescale quantities from tile 1 to new tile area
            snow3%snow(il)%wc_em(it) = snow3%snow(il)%wc_em(it) * x1
            snow3%snow(il)%wc_im(it) = snow3%snow(il)%wc_im(it) * x1
          enddo

        else if (il .le. snow2%nlayers) then
          snow3%snow(il)    = snow2%snow(il)
          snow3%snow(il)%T    = snow2%snow(il)%T
          snow3%snow(il)%age    = snow2%snow(il)%age
          snow3%snow(il)%sph    = snow2%snow(il)%sph
          snow3%snow(il)%optd    = snow2%snow(il)%optd
          snow3%snow(il)%dendr    = snow2%snow(il)%dendr
          snow3%snow(il)%dz = snow3%snow(il)%dz * x2
          snow3%snow(il)%ws = snow3%snow(il)%ws * x2
          snow3%snow(il)%wl = snow3%snow(il)%wl * x2
          do it = 1, NTRACERS ! rescale quantities from tile 1 to new tile area
            snow3%snow(il)%wc_em(it) = snow3%snow(il)%wc_em(it) * x2
            snow3%snow(il)%wc_im(it) = snow3%snow(il)%wc_im(it) * x2
          enddo
        else
          call land_error_message("Error in gl_merge_snow_tiles in gl_snow_tile_mod :: something is wrong with the layering", FATAL)
        endif
      endif
    enddo
    !
  endif

  ! ---------
  ! No need to do anything else if there is a single snow layer
  ! if multiple layers, proceed to merge those that are too thin

  if(is_watch_cell()) then
    write(*,*) '#### gl_snow_tile - merging tiles glass checkpoint 2 ####'
    write(*,*) "2 - print resulting snowpack before removing the layers that are too thin"
    call snow3%print()
  endif

  final_nlayers = snow3%nlayers
  ! if ((snow3%nlayers>1).and.(snow3%snow(1)%ws > 1E-6)) then
  if (snow3%nlayers>1) then

    ! init variables to store snow from layers to be removed
    temp_lost_em = 0.0
    temp_lost_im = 0.0
    temp_lost_wl = 0.0
    temp_lost_ws = 0.0
    temp_lost_heat = 0.0
    temp_lost_age_w = 0.0
    temp_lost_sph_w = 0.0
    temp_lost_optd_w = 0.0
    temp_lost_dendr_w = 0.0
    temp_lost_dz = 0.0

    do il = 1, new_nlayers
      if (snow3%snow(il)%ws < thresh_ws) then
        layers2remove(il)=.TRUE.
        final_nlayers = final_nlayers - 1
        ! note: this can give rise to unpysical ice / liquid water equilibrium T
        ! The equilibrium balance between phases will be done at next snow step
        ! not here (it would potentially change the amount of liquid and ice in the snowpack)
        temp_lost_dz = temp_lost_dz + snow3%snow(il)%dz
        temp_lost_ws = temp_lost_ws + snow3%snow(il)%ws
        temp_lost_wl = temp_lost_wl + snow3%snow(il)%wl
        temp_lost_em = temp_lost_em + snow3%snow(il)%wc_em
        temp_lost_im = temp_lost_im + snow3%snow(il)%wc_im
        ! temp_lost_heat = temp_lost_heat + snow3%snow(il)%heat()
        temp_lost_heat = temp_lost_heat + snow3%snow(il)%hCap() * (snow3%snow(il)%T-TFREEZE)
        temp_lost_age_w = temp_lost_age_w + snow3%snow(il)%ws * snow3%snow(il)%age ! ws - weighted average
        temp_lost_sph_w = temp_lost_sph_w + snow3%snow(il)%ws * snow3%snow(il)%sph ! ws - weighted average
        temp_lost_optd_w = temp_lost_optd_w + snow3%snow(il)%ws * snow3%snow(il)%optd ! ws - weighted average
        temp_lost_dendr_w = temp_lost_dendr_w + snow3%snow(il)%ws * snow3%snow(il)%dendr ! ws - weighted average
      endif
    enddo

    ! now do actually remove the layers
    if (final_nlayers > 0) then



      allocate(snowF(final_nlayers))
      count = 1
      do il=1,new_nlayers
        if (.not.layers2remove(il)) then
          snowF(count) = snow3%snow(il)
          count = count + 1
        endif
      enddo
      snow3%snow(1:final_nlayers) = snowF
      snow3%nlayers = final_nlayers

      ! now add all temporary removed things to the [remaining] bottom layer
      ! there will always be at least one layer still here
      snow3%snow(final_nlayers)%age = (snow3%snow(final_nlayers)%age * snow3%snow(final_nlayers)%ws + temp_lost_age_w) / (snow3%snow(final_nlayers)%ws + temp_lost_ws)
      snow3%snow(final_nlayers)%sph = (snow3%snow(final_nlayers)%sph * snow3%snow(final_nlayers)%ws + temp_lost_sph_w) / (snow3%snow(final_nlayers)%ws + temp_lost_ws)
      snow3%snow(final_nlayers)%optd = (snow3%snow(final_nlayers)%optd * snow3%snow(final_nlayers)%ws + temp_lost_optd_w) / (snow3%snow(final_nlayers)%ws + temp_lost_ws)
      snow3%snow(final_nlayers)%dendr = (snow3%snow(final_nlayers)%dendr * snow3%snow(final_nlayers)%ws + temp_lost_dendr_w) / (snow3%snow(final_nlayers)%ws + temp_lost_ws)
      ! snow3%snow(final_nlayers)%T = TFREEZE + (snow3%snow(final_nlayers)%heat() + temp_lost_heat - HLF*temp_lost_wl - HLF*snow3%snow(final_nlayers)%wl  )/( snow3%snow(final_nlayers)%hCap() + CSW*temp_lost_ws + CLW*temp_lost_wl)
      snow3%snow(final_nlayers)%T = TFREEZE + (snow3%snow(final_nlayers)%hcap()*(snow3%snow(final_nlayers)%T-TFREEZE) + temp_lost_heat   )/( snow3%snow(final_nlayers)%hCap() + CSW*temp_lost_ws + CLW*temp_lost_wl)
      snow3%snow(final_nlayers)%ws = snow3%snow(final_nlayers)%ws + temp_lost_ws
      snow3%snow(final_nlayers)%wl = snow3%snow(final_nlayers)%wl + temp_lost_wl
      snow3%snow(final_nlayers)%wc_em = snow3%snow(final_nlayers)%wc_em + temp_lost_em
      snow3%snow(final_nlayers)%wc_im = snow3%snow(final_nlayers)%wc_im + temp_lost_im
      snow3%snow(final_nlayers)%dz = snow3%snow(final_nlayers)%dz + temp_lost_dz
      ! snow3%snow(final_nlayers)%T = TFREEZE + (snow3%snow(final_nlayers)%heat() + temp_lost_heat - HLF*temp_lost_wl - HLF*snow3%snow(final_nlayers)%wl  )/()

    else

      if(is_watch_cell()) then
        write(*,*) '#### gl_snow_tile - merging tiles glass checkpoint 3 ####'
        write(*,*) "3 - case of all layers being too thin and being removed"
        write(*,*) "3 - fractions: w1, w2, x1, x2 = ", w1, w2, x1, x2
      endif
      ! case in which all layers are being removed ( ws < threshold for all)
      ! final_nlayers=1
      allocate(snowF(1))
      snow3%snow(1) = snowF(1)
      ! snow3%snow(1)%age = (temp_lost_age_w) / (temp_lost_ws)
      ! snow3%snow(1)%sph = (temp_lost_sph_w) / (temp_lost_ws)
      ! snow3%snow(1)%optd = (temp_lost_optd_w) / (temp_lost_ws)
      ! snow3%snow(1)%dendr = ( temp_lost_dendr_w) / ( temp_lost_ws)
      ! snow3%snow(1)%T = TFREEZE + ( temp_lost_heat   )/( CSW*temp_lost_ws + CLW*temp_lost_wl)
      snow3%snow(1)%ws = temp_lost_ws
      snow3%snow(1)%wl = temp_lost_wl
      snow3%snow(1)%wc_em = temp_lost_em
      snow3%snow(1)%wc_im = temp_lost_im
      snow3%snow(1)%dz = temp_lost_dz
      snow3%nlayers = 1
      ! now to avoid numerical instabilities, if there is very little snow
      ! assign to resulting snow the properties of first [largest] tile
      ! and add the remaining heat to the topwater heat (topwheat)
      if (temp_lost_ws > min_snow_mass_merging) then
        snow3%snow(1)%age = (temp_lost_age_w) / (temp_lost_ws)
        snow3%snow(1)%sph = (temp_lost_sph_w) / (temp_lost_ws)
        snow3%snow(1)%optd = (temp_lost_optd_w) / (temp_lost_ws)
        snow3%snow(1)%dendr = ( temp_lost_dendr_w) / ( temp_lost_ws)
        snow3%snow(1)%T = TFREEZE + ( temp_lost_heat   )/( CSW*temp_lost_ws + CLW*temp_lost_wl)
      else

        ! determine from which tile properties are taken
        if (snow2%nlayers==0) then
          use_first_tile = .true.
        else if (snow1%nlayers==0) then
          use_first_tile = .false.
        else
          ! case of both non-zero number of layers
          if (x1 > x2) then
            use_first_tile = .true.
          else
            use_first_tile = .false.
          endif
        endif
        if (use_first_tile) then ! get the properties from first tile, surface layer
          snow3%snow(1)%age   = snow1%snow(1)%age
          snow3%snow(1)%sph   = snow1%snow(1)%sph
          snow3%snow(1)%optd  = snow1%snow(1)%optd
          snow3%snow(1)%dendr = snow1%snow(1)%dendr
          snow3%snow(1)%T     = snow1%snow(1)%T
          ! snow3%snow(1)%T     = TFREEZE
        else ! (main_tile==2) ! get the properties from the second tile, surface layer
          snow3%snow(1)%age   = snow2%snow(1)%age
          snow3%snow(1)%sph   = snow2%snow(1)%sph
          snow3%snow(1)%optd  = snow2%snow(1)%optd
          snow3%snow(1)%dendr = snow2%snow(1)%dendr
          snow3%snow(1)%T     = snow2%snow(1)%T
          ! snow3%snow(1)%T     = TFREEZE
        endif
        orig_heat_1 = 0.0
        orig_heat_2 = 0.0
        if (snow1%nlayers>0) then
          do ill = 1, snow1%nlayers
            orig_heat_1 = orig_heat_1 + (snow1%snow(ill)%T - TFREEZE)* (CSW * snow1%snow(ill)%ws + CLW * snow1%snow(ill)%wl )
          enddo
        endif
        if (snow2%nlayers>0) then
          do ill = 1, snow2%nlayers
            orig_heat_2 = orig_heat_2 + (snow2%snow(ill)%T - TFREEZE)* (CSW * snow2%snow(ill)%ws + CLW * snow2%snow(ill)%wl )
          enddo
        endif
        true_heat_3 = x1 * orig_heat_1 + x2 * orig_heat_2 ! no change of phase, only Ts, wl is conserved
        fict_heat_3 = (snow3%snow(1)%T - TFREEZE) * (CSW * snow3%snow(1)%ws + CLW * snow3%snow(1)%wl )
        ! fict_heat_3 = 0.0 ! T = TFREEZE
        ! add / subtract excess heat from topwater heat (topwheat)
        snow3%topwheat = snow3%topwheat + true_heat_3 - fict_heat_3

      endif
    endif

  endif ! end case of nlayers > 1 [merge small layers]

  ! if a single remains and does not contain enough snow,
  ! add heat and water to top water [does not conserve ice and liquid, but only total water]
  ! and remove layer



  ! ! deal now with the remaining case of the top layer potentially being too thin
    ! if ((snow3%nlayers>1).and.(snow3%snow(1)%ws<thresh_ws)) then
    !   snow3%snow(2)%T = TFREEZE  + ((snow3%snow(2)%T - TFREEZE) * snow3%snow(2)%hCap() +  (snow3%snow(1)%T - TFREEZE) * snow3%snow(1)%hCap() ) / (snow3%snow(2)%hCap() + snow3%snow(1)%hCap())
    !   snow3%snow(2)%age = (snow3%snow(2)%age * snow3%snow(2)%ws + snow3%snow(1)%age * snow3%snow(1)%ws )/(snow3%snow(2)%ws + snow3%snow(1)%ws)
    !   snow3%snow(2)%sph = (snow3%snow(2)%sph * snow3%snow(2)%ws + snow3%snow(1)%sph * snow3%snow(1)%ws )/(snow3%snow(2)%ws + snow3%snow(1)%ws)
    !   snow3%snow(2)%optd = (snow3%snow(2)%optd * snow3%snow(2)%ws + snow3%snow(1)%optd * snow3%snow(1)%ws )/(snow3%snow(2)%ws + snow3%snow(1)%ws)
    !   snow3%snow(2)%dendr = (snow3%snow(2)%dendr * snow3%snow(2)%ws + snow3%snow(1)%dendr * snow3%snow(1)%ws )/(snow3%snow(2)%ws + snow3%snow(1)%ws)
    !   snow3%snow(2)%dz = snow3%snow(2)%dz + snow3%snow(1)%dz
    !   snow3%snow(2)%ws = snow3%snow(2)%ws + snow3%snow(1)%ws
    !   snow3%snow(2)%wl = snow3%snow(2)%wl + snow3%snow(1)%wl
    !   snow3%snow(2)%wc_im = snow3%snow(2)%wc_im + snow3%snow(1)%wc_im
    !   snow3%snow(2)%wc_em = snow3%snow(2)%wc_em + snow3%snow(1)%wc_em
    !   snow3%snow(1:snow3%nlayers-1) = snow3%snow(2:snow3%nlayers)
    !   snow3%nlayers = snow3%nlayers - 1
    ! endif









  ! IF TOO LITTLE SNOW AFTER MERGING, SWEEP ENTIRE SNOW TO TOPWATER
  snow_mass  = snow3%ice() - snow3%topsnowdeficit ! SOLID MASS OF SNOW LAYERS ONLY, WITHOUT DEFICIT IF ANY
  total_lai_loss = 0.0
  total_water_loss = 0.0
  total_heat_loss = 0.0
  ! total_lai_loss = sum(final_lost_em + final_lost_im)

  ! ! Instead of here, do relayering and sweep tiny snow at the beginning of the next update_land_model_fast_0d call

  if ((snow3%nlayers>0).and.( snow_mass < min_snow_mass_merging)) then
    ! delete entire snowpack, and store water mass and heat as topwater and topwheat
    ! total_mass = snow3%SWE() - snow3%topsnowdeficit ! do not modify this term
    ! total_heat = snow3%heat() - snow3%topsnowheatdeficit ! do not modify this term
    total_mass = snow3%topwater
    total_heat = snow3%topwheat
    total_lai_loss = sum(snow3%lai_im() + snow3%lai_em())
    do il=1,snow3%nlayers
      total_heat = total_heat + snow3%snow(il)%heat()
      total_mass = total_mass + snow3%snow(il)%ws + snow3%snow(il)%wl
    enddo
    ! total_heat_loss =
    DEALLOCATE(snow3%snow)
    snow3%nlayers = 0
    ! snow3%topsnowdeficit = 0.0 ! do not modify this term
    ! snow3%topsnowheatdeficit = 0.0 ! do not modify this term
    snow3%topwater = total_mass
    snow3%topwheat = total_heat
  ! else
  ! ! INSTEAD, IF ENOUGH SNOW: PERFORM SNOW RELAYERING TO LIMIT NUMBER OF LAYERS
  !   if (snow3%nlayers > 0) then
  !     ! write(*,*) "numbers of snow layers before relayering = ", s%nlayers
  !     ! if (do_merge) call snow3%attempt_merge_layers()
  !     call snow3%attempt_merge_layers()
  !     ! write(*,*) "numbers of snow layers during relayering (before split, after merge) = ", s%nlayers
  !     ! if (do_split) call snow3%attempt_split_layers()
  !     call snow3%attempt_split_layers()
  !     ! write(*,*) "numbers of snow layers after relayering = ", s%nlayers
  !   endif
  endif

  if(is_watch_cell()) then

    write(*,*) '#### gl_snow_tile - merging tiles glass checkpoint 4 ####'
    write(*,*) "after having removed the layers that are too thin ..."
    write(*,*) "fractions: w1, w2, x1, x2 = ", w1, w2, x1, x2
    write(*,*) "SWE1, SWE2, SWE3 = ", snow1%SWE(), snow2%SWE(), snow3%SWE()
    write(*,*) "SWE3 - SWE2*x2 - SWE1*x1 = ", - snow1%SWE() * x1 - snow2%SWE() * x2 + snow3%SWE()
    write(*,*) "heat1, heat2, heat3 = ", snow1%heat(), snow2%heat(), snow3%heat()
    write(*,*) "heat3 - heat2*x2 - heat1*x1 = ", - snow1%heat() * x1 - snow2%heat() * x2 + snow3%heat()
    write(*,*) "snow1%nlayers, snow2%nlayers, minc_nlayers, new_nlayers = ", snow1%nlayers, snow2%nlayers, minc_nlayers, new_nlayers
    write(*,*) "[number of layers in snow3 including thin layers to be removed] new_nlayers = ",new_nlayers
    write(*,*) "[number of layers in snow3 after removing thin layers] final_nlayers = ", final_nlayers
    write(*,*) "temp_lost_ws = ", temp_lost_ws
    write(*,*) "temp_lost_wl = ", temp_lost_wl
    write(*,*) "temp_lost_heat = ", temp_lost_heat
    write(*,*) "First snow object:"
    call snow1%print()
    write(*,*) "Second snow object:"
    call snow2%print()
    write(*,*) "Final snow object:"
    call snow3%print()
  endif

  ! return tile with merged quantities
  snow2ez%sp = snow3

  ! heat, water, em and im conservation checks
  wat3 = snow3%SWE()
  heat3 = snow3%heat()
  lai3 = sum(snow3%lai_im() + snow3%lai_em())
  ! check conservation
  dwat = abs(wat3 - (wat1*x1 + wat2*x2))
  dheat = abs(heat3 - (heat1*x1 + heat2*x2))
  dlai = abs(lai3 + total_lai_loss - (lai1*x1 + lai2*x2))
  if(dlai>1E-6) then
    write(*,*) "snow1%nlayers, snow2%nlayers", snow1%nlayers, snow2%nlayers
    write(*,*) "Initial numbers of layers: snow1%nl, snow2%nl, tot_lai_loss, lai3 = ", snow1%nlayers, snow2%nlayers, total_lai_loss, lai3
    write(*,*) "fractions: w1, w2, x1, x2 = ", w1, w2, x1, x2
    call snow1%print()
    call snow2%print()
    call snow3%print()
    call land_error_message("error in gl_merge_snow_tiles in gl_snow_tile_mod : total LAIs non conserved merging snow tiles!", FATAL)
  endif
  if(dwat>1E-6) then
    write(*,*) "snow1%nlayers, snow2%nlayers", snow1%nlayers, snow2%nlayers
    write(*,*) "Initial numbers of layers: snow1%nl, snow2%nl, tot_lai_loss = ", snow1%nlayers, snow2%nlayers, total_lai_loss
    write(*,*) "fractions: w1, w2, x1, x2 = ", w1, w2, x1, x2
    call snow1%print()
    call snow2%print()
    call snow3%print()
    call land_error_message("error in gl_merge_snow_tiles in gl_snow_tile_mod : total water non conserved merging snow tiles!", FATAL)
  endif
  if(dheat>1E-6) then
    write(*,*) "fractions: w1, w2, x1, x2 = ", w1, w2, x1, x2
    write(*,*) "snow1%nlayers, snow2%nlayers", snow1%nlayers, snow2%nlayers
    write(*,*) "Initial numbers of layers: snow1%nl, snow2%nl, tot_lai_loss = ", snow1%nlayers, snow2%nlayers, total_lai_loss
    write(*,*) "heat1, heat2, heat3 = ", snow1%heat(), snow2%heat(), snow3%heat()
    write(*,*) "heat1, heat2, heat3 = ", snow1%heat() * x1, snow2%heat() * x2, snow3%heat()
    write(*,*) "heat3 - heat2*x2 - heat1*x1 = ", - snow1%heat() * x1 - snow2%heat() * x2 + snow3%heat()
    write(*,*) "snow top heat deficit:", snow1%topsnowheatdeficit*x1, snow2%topsnowheatdeficit*x2, snow3%topsnowheatdeficit
    write(*,*) "snow top heat deficit diff:", - snow1%topsnowheatdeficit*x1 - snow2%topsnowheatdeficit*x2 + snow3%topsnowheatdeficit
    write(*,*) "snow topwater heat:", snow1%topwheat*x1, snow2%topwheat*x2, snow3%topwheat
    write(*,*) "snow topwater heat diff:", - snow1%topwheat*x1 - snow2%topwheat*x2 + snow3%topwheat
    call snow1%print()
    call snow2%print()
    call snow3%print()
    call land_error_message("error in gl_merge_snow_tiles in gl_snow_tile_mod : total heat non conserved merging snow tiles!", FATAL)
  endif
  if(is_watch_point()) then
    write(*,*)'#### gl_merge_snow_tiles: tiles before and after (snow1, snow2, snow3): ####'
    call snow1%print()
    call snow2%print()
    call snow3%print()
endif
end subroutine gl_merge_snow_tiles

! =============================================================================
! returns true if tile fits the specified selector
function gl_snow_is_selected(snow, sel)
  logical gl_snow_is_selected
  type(tile_selector_type),  intent(in) :: sel
  class(gl_snow_tile_type),      intent(in) :: snow

  gl_snow_is_selected = .TRUE.
end function gl_snow_is_selected

! ============================================================================
subroutine gl_snow_roughness(snow, snow_z0s, snow_z0m)
  class(gl_snow_tile_type), intent(in) :: snow ! not used
  real, intent(out):: snow_z0s, snow_z0m
  snow_z0m =  z0_momentum
  snow_z0s =  z0_momentum * exp(-k_over_B)
end subroutine gl_snow_roughness

! returns snow radiative properties: short-wave refletances (by spectral band),
! long-wave reflecatanc, emissivity
subroutine gl_snow_rad_prop (snow, cosz, subs_refl_dif, p_atm, on_glacier, &
                             snow_refl_dir, snow_refl_dif, snow_refl_lw, snow_emis)
  class(gl_snow_tile_type), intent(inout) :: snow
  real, intent(in) :: cosz
  real, intent(in) :: subs_refl_dif(:)
  real, intent(in) :: p_atm
  logical, intent(in) :: on_glacier
  real, intent(out) :: snow_refl_dir(:), snow_refl_dif(:)
  real, intent(out) :: snow_refl_lw, snow_emis

  real :: snow_top_temp

  if (snow%snow_active()) then
      snow_top_temp = snow%sfc_temp()
  else
      snow_top_temp = TFREEZE ! NOT used in this case
  endif
  call snow_lw_properties ( snow_top_temp, snow_refl_lw, snow_emis )
  call gl_compute_snow_albedo ( &
             snow%sp, snow_top_temp, cosz, on_glacier, p_atm, subs_refl_dif, & ! input
             snow_refl_dir, snow_refl_dif )
end subroutine gl_snow_rad_prop

real function gl_snow_tile_heat (snow) result(heat)
  class(gl_snow_tile_type),  intent(in)    :: snow
  integer :: il
  heat = 0
  heat = heat + snow%sp%topsnowheatdeficit ! add deficit
  heat = heat - snow%sp%topsnowdeficit * HLF ! convention LM4p2 heat wrt liquid at T=TF
  heat = heat + snow%sp%topwheat ! all liquid water heat
  heat = heat - snow%sp%topwater * HLF ! subtract latent heat from top liquid ! //FIXME
  if (snow%sp%nlayers > 0) then
    do il = 1,snow%sp%nlayers
      heat = heat - snow%sp%snow(il)%ws * HLF &
          + ( CLW * snow%sp%snow(il)%wl + CSW * snow%sp%snow(il)%ws )  &
                                        * (snow%sp%snow(il)%T-TFREEZE)
    enddo
  endif
end function gl_snow_tile_heat

! returns true if snow plays a role
function gl_snow_active(snow) ; logical gl_snow_active
  class(gl_snow_tile_type), intent(in)  :: snow
  gl_snow_active = (snow%sp%nlayers > 0)
end function gl_snow_active

real function gl_snow_get_sfc_temp(snow)
  class(gl_snow_tile_type), intent(in) :: snow

  if (snow%sp%nlayers > 0) then
    gl_snow_get_sfc_temp = snow%sp%snow(1)%T
    ! call snow%sp%nearsurf_properties()
    ! snow_T = snow%sp%nearsurf_T
  else
    call land_error_message("gl_snow_get_sfc_temp in gl_snow_tile.F90:: sfc temperature requested, but no snow on the ground!", FATAL)
  endif
end function

integer function gl_snow_nlayers(snow)
  class(gl_snow_tile_type), intent(in) :: snow
  gl_snow_nlayers = snow%sp%nlayers
end function


real function gl_snow_get_Ti(snow, i) result(res)
  class(gl_snow_tile_type), intent(in) :: snow
  integer, intent(in) :: i
  res = snow%sp%snow(i)%T
end function

real function gl_snow_get_wsi(snow, i) result(res)
  class(gl_snow_tile_type), intent(in) :: snow
  integer, intent(in) :: i
  res = snow%sp%snow(i)%ws
end function

real function gl_snow_get_wli(snow, i) result(res)
  class(gl_snow_tile_type), intent(in) :: snow
  integer, intent(in) :: i
  res = snow%sp%snow(i)%wl
end function


subroutine gl_snow_set_Ti(snow, i, v)
  class(gl_snow_tile_type), intent(inout) :: snow
  integer, intent(in) :: i
  real, intent(in) :: v
  snow%sp%snow(i)%T = v
end subroutine

subroutine gl_snow_set_wsi(snow, i, v)
  class(gl_snow_tile_type), intent(inout) :: snow
  integer, intent(in) :: i
  real, intent(in) :: v
  snow%sp%snow(i)%ws = v
end subroutine

subroutine gl_snow_set_wli(snow, i, v)
  class(gl_snow_tile_type), intent(inout) :: snow
  integer, intent(in) :: i
  real, intent(in) :: v
  snow%sp%snow(i)%wl = v
end subroutine


real function gl_snow_get_total_ice(snow) result(ice)
  class(gl_snow_tile_type), intent(in) :: snow
  ! ice = sum(snow%ws(:))
  ice = snow%sp%ice()
end function

real function gl_snow_get_total_liq(snow) result(liq)
  class(gl_snow_tile_type), intent(in) :: snow
  ! liq = sum(snow%wl(:))
  liq = snow%sp%liq()
end function

subroutine gl_snow_get_depth_area(snow, snow_depth, snow_area)
  class(gl_snow_tile_type), intent(in) :: snow
  real, intent(out) :: snow_depth, snow_area
  snow_depth = snow%sp%depth()
  call snow_data_area ( snow_depth, snow_area )
end subroutine

subroutine gl_snow_lai_em(snow, tracers)
  class(gl_snow_tile_type), intent(in) :: snow
  real, intent(out)                    :: tracers(:)

  tracers(1:NTRACERS) = snow%sp%lai_em()
end subroutine

subroutine gl_snow_lai_im(snow, tracers)
  class(gl_snow_tile_type), intent(in) :: snow
  real, intent(out)                    :: tracers(:)

  tracers(1:NTRACERS) = snow%sp%lai_im()
end subroutine

subroutine gl_sweep_snow(snow, lrunf, frunf, hlrunf, hfrunf, lost_wc_em, lost_wc_im)
  class(gl_snow_tile_type), intent(inout) :: snow
  real, intent(out) :: lrunf, frunf   ! liquid and solid runoff generated by sweeper
  real, intent(out) :: hlrunf, hfrunf ! heat carried by liquid ans solid runoff
  real, intent(out) :: lost_wc_em(:), lost_wc_im(:) ! tracer losses

  real :: lswept, fswept, hlswept, hfswept
  real :: swept_wc_em(NTRACERS), swept_wc_im(NTRACERS)

  call gl_sweep_tiny_snow(snow%sp, lrunf,  frunf,  hlrunf,  hfrunf,  lost_wc_em,  lost_wc_im)  ! sweeping tiny snow
  call gl_sweep_huge_snow(snow%sp, lswept, fswept, hlswept, hfswept, swept_wc_em, swept_wc_im) ! sweeping huge snow
  lrunf  = lrunf + lswept
  frunf  = frunf + fswept
  hlrunf = hlrunf + hlswept
  hfrunf = hfrunf + hfswept
  lost_wc_em = lost_wc_em + swept_wc_em
  lost_wc_im = lost_wc_im + swept_wc_im
  hlrunf = hlrunf - lrunf * HLF      ! switch to lm4p2 energy conv

  ! DO A RELAYERING HERE - BEST IF TILE MERGING OCCURRED AND THIN LAYERS HAVE BEEN PRODUCED
  call snow%sp%attempt_merge_layers() ! // TODO
end subroutine

subroutine gl_partition_sw( &
   snow, fswg, fswg_dir, fswg_dif, & ! input
   fswg_substrate, fswg_surface) ! output
   !
   ! Given the shortwave radiation absorbed by snow + substrate (fswg) [W/m2]
   ! as well its direct and diffuse components (fswg_dir, fswg_dif)
   ! partition it between surface of snow (where it was absorbed entirely in old cm snow model)
   ! and, if requested, absorption within the snowpack
   ! andabsoirption in the underlying substrate (lake/soil/glacier)
   !
   class(gl_snow_tile_type), intent(inout) :: snow !< state of snowpack
   real, intent(in)  :: fswg ! total sw absorbed by snow + substrate [W/m2]
   real, intent(in)  :: fswg_dir(:), fswg_dif(:) ! total sw absorbed by snow + substrate (dir only, dif only, by spectral band) [W/m2]
   real, intent(out) :: fswg_substrate ! sw radiation passed to substrate [W/m2]
   real, intent(out) :: fswg_surface   ! sw radiation to be absorbed at surface [W/m2]

   integer il
   real, dimension(NBANDS) :: sum_sw_frac_dir, sum_sw_frac_dif

   ! SNICAR computed flux absorbed in each snow layer for unit of incident flux
   ! snow%sp%sw_frac_dir(il, 1) = snow%sp%sw_frac_dir(il, 1) * fswg_dir(1)
   ! snow%sp%sw_frac_dir(il, 2) = snow%sp%sw_frac_dir(il, 2) * fswg_dir(2)
   ! snow%sp%sw_frac_dif(il, 1) = snow%sp%sw_frac_dif(il, 1) * fswg_dif(1)
   ! snow%sp%sw_frac_dif(il, 2) = snow%sp%sw_frac_dif(il, 2) * fswg_dif(2)

   if (ALLOCATED(snow%sp%swheat)) DEALLOCATE(snow%sp%swheat)

   ! slm: does the code below assume that there are always will be 2 bands?
   if (albedo_option==ALBEDO_SNICAR) then
      if ((use_internal_sources) .and. ((snow%sp%depth() > thresh_snow_depth_swheat) &
                                 .and. (snow%sp%nlayers > 0))) then
         ALLOCATE(snow%sp%swheat(snow%sp%nlayers))
         sum_sw_frac_dir = 0.0 ! init total fractions of sw down absorbed by snowpack
         sum_sw_frac_dif = 0.0 ! init total fractions of sw down absorbed by snowpack
         do il=1,snow%sp%nlayers
            snow%sp%swheat(il) = &
                  fswg_dir(1) * snow%sp%sw_frac_dir(il, 1) + &
                  fswg_dif(1) * snow%sp%sw_frac_dif(il, 1) + &
                  fswg_dir(2) * snow%sp%sw_frac_dir(il, 2) + &
                  fswg_dif(2) * snow%sp%sw_frac_dif(il, 2)
            sum_sw_frac_dir(1)  = sum_sw_frac_dir(1) + snow%sp%sw_frac_dir(il, 1)
            sum_sw_frac_dif(1)  = sum_sw_frac_dif(1) + snow%sp%sw_frac_dif(il, 1)
            sum_sw_frac_dir(2)  = sum_sw_frac_dir(2) + snow%sp%sw_frac_dir(il, 2)
            sum_sw_frac_dif(2)  = sum_sw_frac_dif(2) + snow%sp%sw_frac_dif(il, 2)
         enddo
         if ((sum_sw_frac_dir(1)>1.0+1E-7).or. (sum_sw_frac_dif(1)>1.0+1E-7) .or. &
            (sum_sw_frac_dir(2) >1.0+1E-7).or. (sum_sw_frac_dif(2) >1.0+1E-7)  ) then
            write(*,*) "sum of sw_frac_dir(1):", sum_sw_frac_dir(1)
            write(*,*) "sum of sw_frac_dir(2):", sum_sw_frac_dir(2)
            write(*,*) "sum of sw_frac_dif(1):", sum_sw_frac_dif(1)
            write(*,*) "sum of sw_frac_dif(2):", sum_sw_frac_dif(2)
            call land_error_message("Error in sw sources from SNICAR: a total is larger than 1!", severity=FATAL)
         endif
         if ((sum_sw_frac_dir(1)<0.0-1E-7).or. (sum_sw_frac_dif(1)<0.0-1E-7) .or. &
            (sum_sw_frac_dir(2) <0.0-1E-7).or. (sum_sw_frac_dif(2) <0.0-1E-7)  ) then
            write(*,*) "sum of sw_frac_dir(1):", sum_sw_frac_dir(1)
            write(*,*) "sum of sw_frac_dir(2):", sum_sw_frac_dir(2)
            write(*,*) "sum of sw_frac_dif(1):", sum_sw_frac_dif(1)
            write(*,*) "sum of sw_frac_dif(2):", sum_sw_frac_dif(2)
            call land_error_message("Error in sw sources from SNICAR: a total is below 0!", severity=FATAL)
         endif
         fswg_surface = 0.0
         fswg_substrate = &
               fswg_dir(1) * (1.0 - sum_sw_frac_dir(1)) + &
               fswg_dif(1) * (1.0 - sum_sw_frac_dif(1)) + &
               fswg_dir(2) * (1.0 - sum_sw_frac_dir(2)) + &
               fswg_dif(2) * (1.0 - sum_sw_frac_dif(2))
      else ! albedo = snicar, but do not use internal sw sources
         if (snow%sp%nlayers>0) then
            ALLOCATE(snow%sp%swheat(snow%sp%nlayers))
            snow%sp%swheat = 0.0 ! do not change fswg in this case
         else
            ALLOCATE(snow%sp%swheat(1))
            snow%sp%swheat = 0.0 ! do not change fswg in this case
         endif
         fswg_surface = fswg
         fswg_substrate = 0.0
      endif
   else ! albedo model not SNICAR
      if ((use_internal_sources) .and. ((snow%sp%depth() > thresh_snow_depth_swheat) &
                                 .and. (snow%sp%nlayers > 0))) then
         call snow%sp%sw_sources(fswg_dir, fswg_dif, fswg_substrate)
         fswg_surface = 0.0
      else
         ALLOCATE(snow%sp%swheat(snow%sp%nlayers))
         snow%sp%swheat = 0.0
         fswg_surface = fswg
         fswg_substrate = 0.0
      endif
   endif ! end albedo choice for GL snow model option

   if (assign_substrate_sw_to_surface) then
      fswg_surface=fswg_surface + fswg_substrate
      fswg_substrate = 0.0
   endif

   if (is_watch_point()) then
      write(*,*) "##### gl_partition_sw checkpoint 1: #####"
      __DEBUG1__(fswg)
      __DEBUG1__(fswg_surface)
      __DEBUG1__(fswg_substrate)
      __DEBUG1__(snow%sp%swheat)
   endif
end subroutine gl_partition_sw


subroutine gl_snow_step_1( snow, p_surf, grnd_T, snow_G_Z, snow_G_TZ, &
       snow_active, snow_T, snow_rh, snow_liq, snow_ice, &
       snow_subl, snow_area, snow_G0, snow_DGDT, snow_E_max )
  class(gl_snow_tile_type), intent(inout) :: snow
  real,    intent(in) :: p_surf
  real,    intent(in) :: grnd_T
  real,    intent(in) :: snow_G_Z
  real,    intent(in) :: snow_G_TZ
  logical, intent(out):: snow_active
  real,    intent(out):: &
       snow_T, snow_rh, snow_liq, snow_ice, &
       snow_subl, snow_area, snow_G0, snow_DGDT, &
       snow_E_max

  real, parameter :: atmos_T = 273.15 ! slm: per Enrico's comments, it is not needed, remove later

  call snow%sp%step1a(  &              ! input
     snow_active, snow_T, snow_rh, snow_liq, snow_ice, &   ! output
     snow_subl, snow_area, snow_E_max, delta_time, do_mgimplicit, grnd_T)

  if(is_watch_point()) then
     write(*,*) "##### Check after glass snow step 1a #####"
     __DEBUG4__( snow_active, snow_T, snow_liq, snow_ice)
  endif

  ! NOTE: moved it here now that albedo pre-calculation is done before
  call snow%sp%step1b( snow_G_Z, snow_G_TZ,   snow_G0, snow_DGDT,  atmos_T,  p_surf,    delta_time ) ! for all ez models, regardless of albedo
end subroutine

subroutine gl_snow_step_2 ( snow, snow_subl,                     &
     vegn_lprec, vegn_fprec, vegn_hlprec, vegn_hfprec, &
     DTg,  Mg_imp,  evapg,  fswg,  flwg,  sensg,  &
     use_tfreeze_in_grnd_latent, &
     ! output
     subs_DT, &
     subs_M_imp, subs_evap, subs_fsw, subs_flw, subs_sens,  &
     snow_fsw, snow_flw, snow_sens, &
     snow_levap, snow_fevap, snow_melt, &
     snow_lprec, snow_hlprec, snow_lrunf, snow_frunf, &
     snow_hlrunf, snow_hfrunf, snow_Tbot, snow_Cbot, snow_C, &
     snow_avrg_T , &
     ! additional input/output added by Enrico for standalone model only
     !    snow_rho, snow_age, snow_sph, snow_optd, & ! average snow properties
     ! heat1, verbose, hfevap, dt, wind_atm, t_atm, &
     dt, wind_atm, t_atm, p_surf, &
     wetdep, drydep, grnd_T_preprec, &
     ! for conservation checks only :
     begw_check, begh_check, &
     G0, DGDTg, snow_G_Z, snow_G_TZ, &
     mass_lai_em_1, mass_lai_im_1, &
     lost_wc_em_st, lost_wc_im_st, &
     lost_wc_em, lost_wc_im)
  class(gl_snow_tile_type), intent(inout) :: snow
  real, intent(in) :: &
     snow_subl, vegn_lprec, vegn_fprec, vegn_hlprec, vegn_hfprec
  real, intent(in) :: &
     DTg, Mg_imp, evapg, fswg, flwg, sensg
  logical, intent(in) :: use_tfreeze_in_grnd_latent
  real, intent(out) :: &
         subs_DT, subs_M_imp, subs_evap, subs_fsw, subs_flw, subs_sens, &
         snow_fsw, snow_flw, snow_sens, &
         snow_levap, snow_fevap, snow_melt, &
         snow_lprec, snow_hlprec, snow_lrunf, snow_frunf, &
         snow_hlrunf, snow_hfrunf, snow_Tbot, snow_Cbot, snow_C, snow_avrg_T
   real, intent(out) :: grnd_T_preprec
   !  additional in-out variables
   !  real, intent(out) :: heat1
   !  real, intent(out) :: hfevap ! heat released by subl [kg m^-2 s^-1]
   !  logical, intent(in) :: verbose
   !  real, intent(out) :: delta_heat_DTg
   real, intent(in) :: dt ! delta time step
   real, intent(in) :: wind_atm, t_atm, p_surf
   real, intent(in) :: wetdep(:) ! wet deposition rate of tracers from atmosphere [ppm]
   real, intent(in) :: drydep(:) ! dry deposition rate of tracers from atmosphere [mg/m2/s]
   real, intent(out):: lost_wc_em(:), lost_wc_im(:)
   real, intent(in) :: mass_lai_em_1(:), mass_lai_im_1(:) ! mass of LAIs at beginning of step, for mass cons checks
   real, intent(in) :: lost_wc_em_st(:), lost_wc_im_st(:)
   real, intent(in) :: begw_check, begh_check
   real, intent(in) :: G0, DGDTg, snow_G_Z, snow_G_TZ

   call gl_snow_step_2_ev (snow%sp, snow_subl,                     &
                     vegn_lprec, vegn_fprec, vegn_hlprec, vegn_hfprec, &
                     DTg,  Mg_imp,  evapg,  fswg,  flwg,  sensg,  &
                     use_tfreeze_in_grnd_latent, &
                     ! output
                     subs_DT, &
                     subs_M_imp, subs_evap, subs_fsw, subs_flw, subs_sens,  &
                     snow_fsw, snow_flw, snow_sens, &
                     snow_levap, snow_fevap, snow_melt, &
                     snow_lprec, snow_hlprec, snow_lrunf, snow_frunf, &
                     snow_hlrunf, snow_hfrunf, snow_Tbot, snow_Cbot, snow_C, &
                     snow_avrg_T , &
                     ! additional input/output added by Enrico for standalone model only
                     !    snow_rho, snow_age, snow_sph, snow_optd, & ! average snow properties
                     ! heat1, verbose, hfevap, dt, wind_atm, t_atm, &
                     dt, wind_atm, t_atm, p_surf, &
                     wetdep, drydep, grnd_T_preprec, &
                     ! for conservation checks only :
                     begw_check, begh_check, &
                     G0, DGDTg, snow_G_Z, snow_G_TZ, &
                     mass_lai_em_1, mass_lai_im_1, &
                     lost_wc_em_st, lost_wc_im_st, &
                     lost_wc_em, lost_wc_im)
end subroutine

! send the diagnostics
subroutine gl_snow_send_diag(snow, diag)
  class(gl_snow_tile_type),  intent(inout) :: snow !< snow data structure
  type(diag_buff_type), intent(inout) :: diag !< diagnostic buffer

  real :: snow_area

  snow_area = snow%sp%area()
  call snow%sp%nearsurf_properties() ! slm: this updates snow state somehow

   ! if(snow%nlayers > 0) then
  call send_tile_data(id_snow_avrg_optd, snow_area * snow%sp%avrg_optd(), diag)
  call send_tile_data(id_snow_avrg_sph, snow_area * snow%sp%avrg_sph(), diag)
  call send_tile_data(id_snow_avrg_age, snow_area * snow%sp%avrg_age(), diag)
  call send_tile_data(id_snow_avrg_dendr, snow_area * snow%sp%avrg_dendr(), diag)
  call send_tile_data(id_snow_density, snow_area * snow%sp%density(), diag)
!   call send_tile_data(id_snow_avrg_T, snow_area * snow_avrg_T, diag)
  call send_tile_data(id_snow_avrg_bceq_tot, snow_area * snow%sp%avrg_bceq_tot(), diag)
  call send_tile_data(id_snow_avrg_bc_tot, snow_area * snow%sp%avrg_bc_tot(), diag)
  call send_tile_data(id_snow_avrg_md_tot, snow_area * snow%sp%avrg_md_tot(), diag)
  call send_tile_data(id_snow_avrg_om_tot, snow_area * snow%sp%avrg_om_tot(), diag)
  call send_tile_data(id_snow_avrg_bceq_im, snow_area * snow%sp%avrg_bceq_im(), diag)
  call send_tile_data(id_snow_avrg_bceq_em, snow_area * snow%sp%avrg_bceq_em(), diag)
  call send_tile_data(id_snow_nearsurf_bceq_tot, snow_area * snow%sp%nearsurf_bceq_tot, diag)
  call send_tile_data(id_snow_nearsurf_bceq_im, snow_area * snow%sp%nearsurf_bceq_im, diag)
  call send_tile_data(id_snow_nearsurf_bceq_em, snow_area * snow%sp%nearsurf_bceq_em, diag)
  call send_tile_data(id_snow_nearsurf_optd, snow_area * snow%sp%nearsurf_optd, diag)
  call send_tile_data(id_snow_nearsurf_sph, snow_area * snow%sp%nearsurf_sph, diag)
  call send_tile_data(id_snow_nearsurf_density, snow_area * snow%sp%nearsurf_rho, diag)
  call send_tile_data(id_snow_nearsurf_age, snow_area * snow%sp%nearsurf_age, diag)
  call send_tile_data(id_snow_nearsurf_dendr, snow_area * snow%sp%nearsurf_dendr, diag)
  ! endif
  ! snow-related quantities defined also when snow depth = 0 (= no snow layers)
  ! do the follwing vars in update_land_bc_fast, as done in old model version
  ! call send_tile_data(id_snow_area_frac,snow_area_frac,snow%area())
  ! call send_tile_data(id_snow_depth, snow%sp%depth(), diag)
  call send_tile_data(id_snow_liq, snow%sp%liq(), diag)
  call send_tile_data(id_snow_ice, snow%sp%ice(), diag)
  call send_tile_data(id_snow_topwater, snow%sp%topwater, diag)
  call send_tile_data(id_snow_topwheat, snow%sp%topwheat, diag)
  call send_tile_data(id_snow_topsnowdeficit, snow%sp%topsnowdeficit, diag)
  call send_tile_data(id_snow_topsnowheatdeficit, snow%sp%topsnowheatdeficit, diag)

end subroutine gl_snow_send_diag


end module gl_snow_tile_mod
