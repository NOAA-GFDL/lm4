module gl_snow_tile_mod
#include <fms_platform.h>

#ifdef INTERNAL_FILE_NML
use mpp_mod, only: input_nml_file
#else
use fms_mod, only: open_namelist_file
#endif

use fms_mod, only : file_exist, check_nml_error, close_file, stdlog, error_mesg, FATAL, NOTE
use constants_mod,only: tfreeze, hlf
use land_constants_mod, only : NBANDS
use land_tile_selectors_mod, only : tile_selector_type
use land_data_mod, only : log_version

use snow_constants_mod, only: NTRACERS

use snowpack_mod, only : snow_layer_type, snowpack_t, merge_layers
use parent_snow_tile_mod, only : snow_tile_type, mc_fict, z0_momentum, k_over_B, num_l, &
                                  cpw, clw, csw, dz

use land_debug_mod, only : is_watch_point, is_watch_cell, land_error_message

implicit none
private

! ==== public interfaces =====================================================
public :: gl_snow_tile_type
public :: gl_snow_tile_ctor
public :: gl_snow_tile_copy_ctor
public :: gl_delete_snow_tile
public :: gl_snow_tiles_can_be_merged
public :: gl_merge_snow_tiles
public :: gl_snow_is_selected
public :: gl_get_snow_tile_tag
public :: gl_snow_tile_stock_pe
public :: gl_snow_tile_heat ! use that defined in snowpack
public :: gl_snow_active
public :: gl_snow_roughness
public :: gl_snow_get_sfc_temp


! ==== module constants ======================================================
character(len=*), parameter :: module_name = 'gl_snow_tile_mod'
#include "../../shared/version_variable.inc"


! ==== types =================================================================



type, extends(snow_tile_type) :: gl_snow_tile_type

   contains
    procedure :: merge_snow_tiles => gl_merge_snow_tiles_wrapper
    procedure :: get_snow_tile_tag => gl_get_snow_tile_tag

    procedure :: snow_is_selected => gl_snow_is_selected
    procedure :: snow_roughness => gl_snow_roughness
    procedure :: stock_pe => gl_snow_tile_stock_pe
    procedure :: snow_active => gl_snow_active
    procedure :: snow_tile_heat => gl_snow_tile_heat
    procedure :: snow_get_sfc_temp => gl_snow_get_sfc_temp


    procedure :: get_Ti => gl_snow_get_Ti
    procedure :: get_wli => gl_snow_get_wli
    procedure :: get_wsi => gl_snow_get_wsi

    procedure :: set_Ti =>  gl_snow_set_Ti
    procedure :: set_wli => gl_snow_set_wli
    procedure :: set_wsi => gl_snow_set_wsi

    procedure :: ice => gl_snow_get_total_ice
    procedure :: liq => gl_snow_get_total_liq



end type gl_snow_tile_type


contains ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-



! ============================================================================
function gl_snow_tile_ctor(tag) result(ptr)
  type(gl_snow_tile_type), pointer :: ptr ! return value
  integer, optional, intent(in) :: tag ! kind of tile

  allocate(ptr)
  ptr%tag = 0 ; if(present(tag)) ptr%tag = tag

  ! allocate also fields for old snow just in case we need to read them 
  ! from old restart - deallocate after reading restart
  ! allocate(ptr%ws(num_l))
  ! allocate(ptr%wl(num_l))
  ! allocate(ptr%T(num_l))
  ! ptr%nlayers = num_l

  ! allocate(ptr)
  allocate(ptr%ws(num_l))
  allocate(ptr%wl(num_l))
  allocate(ptr%T(num_l))
  allocate(ptr%e(num_l))
  allocate(ptr%f(num_l))
  ptr%nlayers = num_l

end function gl_snow_tile_ctor

! ============================================================================
function gl_snow_tile_copy_ctor(snow) result(ptr)
  type(gl_snow_tile_type), pointer :: ptr ! return value
  type(gl_snow_tile_type), intent(in) :: snow ! tile to copy

  allocate(ptr)
  ! copy all non-pointer members
  ptr = snow
  ! no need to allocate storage for allocatable components of the type, because
  ! F2003 takes care of that, and also takes care of copying data
end function gl_snow_tile_copy_ctor

! ============================================================================
subroutine gl_delete_snow_tile(snow)
  type(gl_snow_tile_type), pointer :: snow

  ! no need to deallocate components of tile, because F2003 takes care of
  ! allocatable components deallocation when tile is deallocated
  ! deallocate(snow)
  nullify(snow)
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
  integer ill, use_first_tile


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
! retruns tag of the tile
function gl_get_snow_tile_tag(snow) result(tag)
  integer :: tag
  class(gl_snow_tile_type), intent(in) :: snow
  tag = snow%tag
end function gl_get_snow_tile_tag

! ============================================================================
subroutine gl_snow_roughness(snow, snow_z0s, snow_z0m)
  class(gl_snow_tile_type), intent(in) :: snow ! not used
  real, intent(out):: snow_z0s, snow_z0m
  snow_z0m =  z0_momentum
  snow_z0s =  z0_momentum * exp(-k_over_B)
end subroutine gl_snow_roughness

! ============================================================================
subroutine gl_snow_tile_stock_pe (snow, twd_liq, twd_sol  )
  class(gl_snow_tile_type),  intent(in)    :: snow
  real,                  intent(out)   :: twd_liq, twd_sol
  twd_liq = snow%sp%liq()
  twd_sol = snow%sp%ice()
end subroutine gl_snow_tile_stock_pe


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
  ! snow_active = ( sum(snow%ws(1:num_l)) > 0 )
  gl_snow_active = (snow%sp%nlayers > 0)
end function gl_snow_active

subroutine gl_snow_get_sfc_temp(snow, snow_T)
  class(gl_snow_tile_type), intent(in) :: snow
  real, intent(out) :: snow_T
  if (snow%sp%nlayers > 0) then
    snow_T = snow%sp%snow(1)%T 
    ! call snow%sp%nearsurf_properties() 
    ! snow_T = snow%sp%nearsurf_T
  else
    call land_error_message("gl_snow_get_sfc_temp in gl_snow_tile.F90:: sfc temperature requested, but no snow on the ground!", FATAL)
  endif
end subroutine


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

end module gl_snow_tile_mod
