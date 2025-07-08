module land_tile_mod

use fms_mod, only : mpp_pe, mpp_root_pe, check_nml_error, error_mesg, stdlog, &
                  & FATAL
use mpp_mod, only: input_nml_file
use land_constants_mod, only : NBANDS
use glac_tile_mod, only : &
     glac_tile_type, new_glac_tile, delete_glac_tile, glac_is_selected, &
     glac_tiles_can_be_merged, merge_glac_tiles, get_glac_tile_tag, &
     glac_tile_stock_pe, glac_tile_heat
use lake_tile_mod, only : &
     lake_tile_type, new_lake_tile, delete_lake_tile, lake_is_selected, &
     lake_tiles_can_be_merged, merge_lake_tiles, get_lake_tile_tag, &
     lake_tile_stock_pe, lake_tile_heat
use soil_tile_mod, only : &
     soil_tile_type, new_soil_tile, delete_soil_tile, soil_is_selected, &
     soil_tiles_can_be_merged, merge_soil_tiles, get_soil_tile_tag, &
     soil_tile_stock_pe, soil_tile_heat
use hillslope_tile_mod, only : hlsp_is_selected
use cana_tile_mod, only : &
     cana_tile_type, new_cana_tile, delete_cana_tile, cana_is_selected, &
     cana_tiles_can_be_merged, merge_cana_tiles, get_cana_tile_tag, &
     cana_tile_stock_pe, cana_tile_carbon, cana_tile_heat
use vegn_tile_mod, only : &
     vegn_tile_type, new_vegn_tile, delete_vegn_tile, vegn_is_selected, &
     vegn_tiles_can_be_merged, vegn_tile_lu_match, merge_vegn_tiles, vegn_tile_tag, &
     vegn_tile_stock_pe, vegn_tile_carbon, vegn_tile_heat, vegn_tile_nitrogen, &
     vegn_tile_bwood
use vegn_util_mod, only : kill_small_cohorts_ppa
use vegn_data_mod, only : landuse_name
! ##### EZSNOW - new snow model #####
use snow_base_mod, only : new_snow_tile, delete_snow_tile, snow_tiles_can_be_merged
use snow_tile_mod, only : snow_tile_type
! ##### end new snow model ######
use soil_BGC_type_mod, only : soil_BGC_t
use soil_BGC_base_mod, only : new_soilc, delete_soilc

use land_tile_selectors_mod, only : tile_selector_type, &
     SEL_SOIL, SEL_VEGN, SEL_LAKE, SEL_GLAC, SEL_SNOW, SEL_CANA, SEL_HLSP
use tile_diag_buff_mod, only : &
     diag_buff_type, init_diag_buff
use land_data_mod, only : lnd, log_version
use land_debug_mod, only : &
     is_watch_cell, check_conservation, &
     water_cons_tol, carbon_cons_tol, nitrogen_cons_tol, heat_cons_tol

implicit none
private
! ==== public interfaces =====================================================
public :: land_tile_type
public :: land_tile_list_type
public :: land_tile_enum_type
public :: diag_buff_type

! operations with tile map
public :: init_tile_map, free_tile_map
public :: max_n_tiles

! operations with tile
public :: new_land_tile, delete_land_tile
public :: merge_land_tiles, merge_land_tile_into_list
public :: remerge_tile_list ! reduces number of tiles by merging all that can be merged

public :: get_tile_water ! returns liquid and frozen water masses
public :: land_tile_carbon ! returns total carbon in the tile
public :: land_tile_nitrogen ! returns total nitrogen in the tile
public :: land_tile_heat ! returns tile heat content
public :: land_tile_grnd_T ! returns temperature of the ground surface

! operations with tile lists and tile list enumerators
public :: land_tile_list_init, land_tile_list_end
public :: first_elmt, tail_elmt
public :: elmt_at_index ! given list and k, returns list[k]
public :: operator(==), operator(/=) ! comparison of two enumerators
public :: next_elmt, prev_elmt ! enumerator advance operations
public :: loop_over_tiles ! provides simple way to iterate over a list of tiles
public :: current_tile ! returns pointer to the tile at a position
public :: insert  ! inserts a tile at a given position, or appends it to a list
public :: erase   ! erases tile at current position
public :: remove  ! removes tile at current position, but does not delete it
public :: get_elmt_indices ! returns i,j,k of current element

public :: empty   ! returns true if the list of tiles is empty
public :: nitems  ! count of items in list

public :: tile_is_selected

! abstract interfaces for accessor functions
public :: tile_test_func, fptr_i0, fptr_i0i, fptr_r0, fptr_r0i, fptr_i0ij, fptr_r0ij, fptr_r0ijk

public :: land_tile_map ! array of tile lists
! ==== end of public interfaces ==============================================

interface new_land_tile
   module procedure land_tile_ctor
   module procedure land_tile_copy_ctor
end interface

interface first_elmt
   module procedure land_tile_list_begin_0d
   module procedure land_tile_list_begin_1d
end interface
interface tail_elmt
   module procedure land_tile_list_end_0d
   module procedure land_tile_list_end_1d
end interface

interface operator(==)
   module procedure enums_are_equal
end interface
interface operator(/=)
   module procedure enums_are_not_equal
end interface

interface insert
   module procedure insert_at_position, append_to_list
end interface
interface remove
   module procedure remove_at_position, remove_all_from_list
end interface
interface erase
   module procedure erase_at_position, erase_all_from_list
end interface
interface nitems
   module procedure n_items_in_list
end interface

! ==== module constants ======================================================
character(len=*), parameter :: module_name = 'land_tile_mod'
#include "shared/version_variable.inc"

! ==== data types ============================================================
! land_tile_type describes the structure of the land model tile; basically
! it is a container for tile-specific data, plus some information common to
! all of them: fraction of tile area, etc.
type :: land_tile_type
   integer :: tag = 0   ! defines type of the tile

   real    :: frac      ! fractional tile area, dimensionless
   type(glac_tile_type),  pointer :: glac  => NULL() ! glacier model data
   type(lake_tile_type),  pointer :: lake  => NULL() ! lake model data
   type(soil_tile_type),  pointer :: soil  => NULL() ! soil model data
   class(snow_tile_type), pointer :: snow  => NULL() ! snow data EZSNOW
   type(cana_tile_type),  pointer :: cana  => NULL() ! canopy air data
   type(vegn_tile_type),  pointer :: vegn  => NULL() ! vegetation model data
   class(soil_BGC_t),     pointer :: soilc => NULL() ! soil carbon data

   type(diag_buff_type) :: diag ! diagnostic data storage

   ! data that are carried from update_land_bc_fast to update_land_fast:
   real :: Sg_dir(NBANDS), Sg_dif(NBANDS) ! fractions of downward direct and
       ! diffuse short-wave radiation absorbed by ground and snow
   ! fractions of downward direct and diffuse radiation absorbed by the
   ! vegetation; dimensions are (NCOHORTS,NBANDS).
   real, allocatable :: Sv_dir(:,:), Sv_dif(:,:)
   ! fractions of downward direct and diffuse radiation on top of each cohort
   ! dimensions are (NCOHORTS,NBANDS).
   real, allocatable :: Sdn_dir(:,:), Sdn_dif(:,:)
   real :: land_refl_dir(NBANDS), land_refl_dif(NBANDS)

   real :: land_d, land_z0m, land_z0s, land_RSL, grnd_z0m, grnd_z0s
   real :: bstar = 0.0 ! turbulent buoyancy scale, m/s2. It is a member of tile structure
   ! and stored in the restarts only because update_land_bc_fast [where it is used to
   ! calculate roughness sublayer depth] is called on initialization before atmos is able
   ! to pass stability data down, so to reproduce across restarts land has to retrieve
   ! previous value from the exiting IC.
   real :: surf_refl_lw ! long-wave reflectivity of the ground surface (possibly snow-covered)
   ! black-background long-wave radiative properties of the vegetation cohorts
   real, allocatable :: vegn_refl_lw(:)  ! reflectance
   real, allocatable :: vegn_tran_lw(:)  ! transmittance

   real :: lwup     = 200.0  ! upward long-wave flux from the entire land (W/m2), the result of
           ! the implicit time step -- used in update_land_bc_fast to return to the flux exchange.
   real :: e_res_1  = 0.0 ! energy residual in canopy air EB equation
   real :: e_res_2  = 0.0 ! energy residual in canopy EB equation
end type land_tile_type

! tile_list_type provides a container for the tiles
type :: land_tile_list_type
   private
   type(land_tile_list_node_type), pointer :: head => NULL()
end type land_tile_list_type

! land_tile_enum_type provides a enumerator of tiles -- a data structure
! that allows to walk through all tiles in a container (or a 2D array of
! containers) without bothering with details of container implementation
type :: land_tile_enum_type
   private
   type(land_tile_list_type), pointer :: &
        tiles(:) => NULL()  ! pointer to array of tiles to walk -- may be disassociated
   integer :: i=0,j=0 ! i,j indices in the above array
   integer :: l=0     ! index in the unstructured domain
   integer :: k=0     ! number of the current tile in its container
   integer :: lo=0    ! offsets of indices (to keep track of non-1 ubounds of tiles array)
   type(land_tile_list_node_type), pointer :: node => NULL() ! pointer to the current container node
end type land_tile_enum_type

! private type -- used internally to implement tile lists
type :: land_tile_list_node_type
   type(land_tile_list_node_type), pointer :: prev => NULL()
   type(land_tile_list_node_type), pointer :: next => NULL()
   type(land_tile_type), pointer :: data => NULL()
end type land_tile_list_node_type

! ==== abstract interfaces ===================================================
abstract interface
  ! the following interface describes the "detector function", which is passed
  ! through the argument list and must return true for any tile to be written
  ! to the specific restart, false otherwise
  logical function tile_test_func(tile)
     import land_tile_type
     type(land_tile_type), pointer :: tile
  end function tile_test_func
  ! the following interfaces describe various accessor subroutines, used to access
  ! data im massive operations on tiles, such as i/o or (sometimes) diagnostics

  ! given land tile, returns pointer to some scalar real data
  ! within this tile, or an unassociated pointer if there is no data
  subroutine fptr_r0(tile, ptr)
     import land_tile_type
     type(land_tile_type), pointer :: tile ! input
     real                , pointer :: ptr  ! returned pointer to the data
  end subroutine fptr_r0
  ! given land tile and an index, returns pointer to some scalar real data
  ! within this tile, or an unassociated pointer if there is no data
  subroutine fptr_r0i(tile, i, ptr)
     import land_tile_type
     type(land_tile_type), pointer :: tile ! input
     integer             , intent(in) :: i ! index in the array
     real                , pointer :: ptr  ! returned pointer to the data
  end subroutine fptr_r0i
  ! given land tile and an 2 indices, returns pointer to some scalar real data
  ! within this tile, or an unassociated pointer if there is no data
  subroutine fptr_r0ij(tile, i,j, ptr)
     import land_tile_type
     type(land_tile_type), pointer :: tile ! input
     integer             , intent(in) :: i,j ! indices in the array
     real                , pointer :: ptr  ! returned pointer to the data
  end subroutine fptr_r0ij
  subroutine fptr_i0ij(tile, i,j, ptr)
     import land_tile_type
     type(land_tile_type), pointer :: tile ! input
     integer             , intent(in) :: i,j ! indices in the array
     integer             , pointer :: ptr  ! returned pointer to the data
  end subroutine fptr_i0ij
  ! given land tile and an 3 indices, returns pointer to some scalar real data
  ! within this tile, or an unassociated pointer if there is no data
  subroutine fptr_r0ijk(tile, i,j,k, ptr)
     import land_tile_type
     type(land_tile_type), pointer :: tile ! input
     integer             , intent(in) :: i,j,k ! indices in the array
     real                , pointer :: ptr  ! returned pointer to the data
  end subroutine fptr_r0ijk

  ! given land tile, returns pointer to some scalar integer data
  ! within this tile, or an unassociated pointer if there is no data
  subroutine fptr_i0(tile, ptr)
     import land_tile_type
     type(land_tile_type), pointer :: tile ! input
     integer             , pointer :: ptr  ! returned pointer to the data
  end subroutine fptr_i0
  ! given land tile and an index, returns pointer to some scalar integer data
  ! within this tile, or an unassociated pointer if there is no data
  subroutine fptr_i0i(tile, i, ptr)
     import land_tile_type
     type(land_tile_type), pointer :: tile ! input
     integer             , intent(in) :: i ! index in the array
     integer             , pointer :: ptr  ! returned pointer to the data
  end subroutine fptr_i0i
  ! NOTE: import statements are needed because in FORTRAN interface blocks
  ! do not have access to their environment by host association, so without
  ! "import" they do not know the definition of land_tile_type, and compilation
  ! fails

  ! given two vegetation tiles, returns TRUE is they are allowed to merge, FALSE
  ! otherwise
  logical function vegn_tiles_merge_check(vegn1, vegn2)
     import vegn_tile_type
     type(vegn_tile_type), intent(in) :: vegn1, vegn2
  end function vegn_tiles_merge_check
end interface

! ==== module data ===========================================================
type(land_tile_list_type), allocatable :: land_tile_map(:) ! map of tiles

real    :: min_tile_frac = 0.0 ! minimum fraction of tile land area that is not
   ! aggressively merged during re-merging of the tiles in remerge_tile_list
! chatacter(32) :: vegn_merge_policy = 'agressive' ! or 'strict-LU-match'
   ! 'agressive' means that tiny tiles can be merged into any large natural, secondary, or
   ! rangeland tile, regardless of their own land use
   ! 'strict-LU-match' means that tiny tiles can be merged only
namelist /tile_merge_nml/ min_tile_frac!, vegn_merge_policy

contains

! #### land_tile_type and operations #########################################

! ============================================================================
! initialize land tile map
subroutine init_tile_map()
  integer :: l
  integer :: unit, ierr, io

  call log_version(version, module_name, &
  __FILE__)

  read (input_nml_file, nml=tile_merge_nml, iostat=io)
  ierr = check_nml_error(io, 'tile_merge_nml')

  if (mpp_pe() == mpp_root_pe()) then
     unit = stdlog()
     write (unit, nml=tile_merge_nml)
  endif

  allocate(land_tile_map(lnd%ls:lnd%le))
  do l = lnd%ls,lnd%le
     call land_tile_list_init(land_tile_map(l))
  enddo
end subroutine init_tile_map

! ============================================================================
! deallocate land tile map
subroutine free_tile_map()
  integer :: l

  do l = lnd%ls,lnd%le
     call land_tile_list_end(land_tile_map(l))
  enddo
end subroutine free_tile_map

! ============================================================================
! get max number of tiles in the domain
function max_n_tiles() result(n)
  integer :: n
  integer :: l

  n=1
  do l=lnd%ls,lnd%le
     n=max(n, nitems(land_tile_map(l)))
  enddo
end function max_n_tiles

! ============================================================================
! tile constructor: given a list of sub-model tile tags, creates a land tile
! calls sub-tile constructors from individual component models
function land_tile_ctor(frac,glac,lake,soil,vegn,tag,htag_j,htag_k) result(tile)
  real   , optional, intent(in) :: frac ! fractional area of tile
  integer, optional, intent(in) :: &
               glac,lake,soil,vegn ! kinds of respective tiles
  integer, optional, intent(in) :: tag  ! general tile tag
  integer, optional, intent(in) :: htag_j  ! optional hillslope position tag
  integer, optional, intent(in) :: htag_k  ! optional hillslope parent tag
  type(land_tile_type), pointer :: tile ! return value

  ! ---- local vars
  integer :: glac_, lake_, soil_, vegn_

  ! initialize internal variables
  glac_ = -1 ; if(present(glac)) glac_ = glac
  lake_ = -1 ; if(present(lake)) lake_ = lake
  soil_ = -1 ; if(present(soil)) soil_ = soil
  vegn_ = -1 ; if(present(vegn)) vegn_ = vegn

  allocate(tile)
  ! fill common fields
  tile%frac = 0.0 ; if(present(frac)) tile%frac = frac
  tile%tag  = 0   ; if(present(tag))  tile%tag  = tag

  ! create sub-model tiles
  tile%cana => new_cana_tile()
  if(glac_>=0) tile%glac => new_glac_tile(glac_)
  if(lake_>=0) tile%lake => new_lake_tile(lake_)
  tile%snow => new_snow_tile()
  if(soil_>=0) then
    if (present(htag_j) .and. present(htag_k)) then
        tile%soil => new_soil_tile(soil_, htag_j, htag_k)
    else
        tile%soil => new_soil_tile(soil_, 0, 0) ! Hillslope model is inactive or
        ! these indices will be set in hlsp_init.
    end if
    tile%soilc => new_soilc(tile%soil)
  end if
  if(vegn_>=0) tile%vegn => new_vegn_tile(vegn_)

  ! create a buffer for diagnostic output
  call init_diag_buff(tile%diag)

end function land_tile_ctor


! ============================================================================
function land_tile_copy_ctor(t) result(tile)
  type(land_tile_type), intent(in) :: t    ! tile to copy
  type(land_tile_type), pointer :: tile ! return value

  allocate(tile)
  tile = t ! copy all non-pointer members
  if (associated(t%glac))  tile%glac=>new_glac_tile(t%glac)
  if (associated(t%lake))  tile%lake=>new_lake_tile(t%lake)
  if (associated(t%soil))  tile%soil=>new_soil_tile(t%soil)
  if (associated(t%snow))  tile%snow=>new_snow_tile(t%snow)
  if (associated(t%cana))  tile%cana=>new_cana_tile(t%cana)
  if (associated(t%vegn))  tile%vegn=>new_vegn_tile(t%vegn)
  if (associated(t%soilc)) tile%soilc=>new_soilc(t%soilc)
end function land_tile_copy_ctor


! ============================================================================
! tile destructor -- releases memory occupied by the tile;
! calls sub-model tile destructors to free the memory of the components
subroutine delete_land_tile(tile)
  type(land_tile_type), pointer :: tile ! tile to delete

  if (.not.associated(tile)) return

  if (associated(tile%glac))  call delete_glac_tile(tile%glac)
  if (associated(tile%lake))  call delete_lake_tile(tile%lake)
  if (associated(tile%soil))  call delete_soil_tile(tile%soil)
  if (associated(tile%snow))  call delete_snow_tile(tile%snow)
  if (associated(tile%cana))  call delete_cana_tile(tile%cana)
  if (associated(tile%vegn))  call delete_vegn_tile(tile%vegn)
  if (associated(tile%soilc)) call delete_soilc(tile%soilc)

  ! release the tile memory
  deallocate(tile)
end subroutine delete_land_tile


! ============================================================================
! returns totals water and ice masses associated with tile
subroutine get_tile_water(tile, lmass, fmass)
  type(land_tile_type), intent(in) :: tile
  real, intent(out) :: lmass, fmass ! liquid and solid water masses, kg/m2

  ! ---- local vars
  real :: lm, fm

  lmass = 0; fmass = 0
  if (associated(tile%cana)) then
     call cana_tile_stock_pe(tile%cana, lm, fm)
     lmass = lmass+lm ; fmass = fmass + fm
  endif
  if (associated(tile%glac)) then
     call glac_tile_stock_pe(tile%glac, lm, fm)
     lmass = lmass+lm ; fmass = fmass + fm
  endif
  if (associated(tile%lake)) then
     call lake_tile_stock_pe(tile%lake, lm, fm)
     lmass = lmass+lm ; fmass = fmass + fm
  endif
  if (associated(tile%soil)) then
     call soil_tile_stock_pe(tile%soil, lm, fm)
     lmass = lmass+lm ; fmass = fmass + fm
  endif
  if (associated(tile%snow)) then
     lm = tile%snow%liq() ; fm = tile%snow%ice()
     lmass = lmass+lm ; fmass = fmass + fm
  endif
  if (associated(tile%vegn)) then
     call vegn_tile_stock_pe(tile%vegn, lm, fm)
     lmass = lmass+lm ; fmass = fmass + fm
  endif

end subroutine get_tile_water


! ============================================================================
! returns total tile carbon, kg C/m2
function land_tile_carbon(tile) result(carbon) ; real carbon
  type(land_tile_type), intent(in) :: tile

  carbon = 0
  if (associated(tile%cana)) &
     carbon = carbon + cana_tile_carbon(tile%cana)
  if (associated(tile%vegn)) &
     carbon = carbon + vegn_tile_carbon(tile%vegn)
  if (associated(tile%soilc)) &
     carbon = carbon + tile%soilc%total_C()
end function land_tile_carbon


! ============================================================================
! returns total tile nitrogen, kg N/m2
function land_tile_nitrogen(tile) result(nitrogen) ; real nitrogen
  type(land_tile_type), intent(in) :: tile

  nitrogen = 0
  ! Not implemented for canopy yet
  ! if (associated(tile%cana)) &
  !    nitrogen = nitrogen + cana_tile_nitrogen(tile%cana)
  if (associated(tile%vegn)) &
     nitrogen = nitrogen + vegn_tile_nitrogen(tile%vegn)
  if (associated(tile%soilc)) &
     nitrogen = nitrogen + tile%soilc%total_N()
end function land_tile_nitrogen


! ============================================================================
! returns total heat content of the tile
function land_tile_heat(tile) result(heat) ; real heat
  type(land_tile_type), intent(in) :: tile

  heat = tile%e_res_1 + tile%e_res_2
  if (associated(tile%cana)) &
       heat = heat+cana_tile_heat(tile%cana)
  if (associated(tile%glac)) &
       heat = heat+glac_tile_heat(tile%glac)
  if (associated(tile%lake)) &
       heat = heat+lake_tile_heat(tile%lake)
  if (associated(tile%soil)) &
       heat = heat+soil_tile_heat(tile%soil)
  if (associated(tile%snow)) &
       heat = heat+tile%snow%snow_tile_heat() ! EZSNOW
  if (associated(tile%vegn)) &
       heat = heat+vegn_tile_heat(tile%vegn)
end function land_tile_heat

! ============================================================================
! returns ground surface temperature
function land_tile_grnd_T(tile) result(T) ; real T
  type(land_tile_type), intent(in) :: tile

  if (tile%snow%snow_active()) then ! always associated
     T = tile%snow%sfc_temp()  ! EZSNOW
  else if (associated(tile%soil)) then
     T = tile%soil%T(1)
  else if (associated(tile%glac)) then
     T = tile%glac%T(1)
  else if (associated(tile%lake)) then
     T = tile%lake%T(1)
  endif
end function land_tile_grnd_T

! ============================================================================
! returns true if tile1 can be merged into tile2
function land_tiles_can_be_merged(tile1,tile2,vegn_merge_check) result (answer)
   logical :: answer ! returned value
   type(land_tile_type), intent(in) :: tile1, tile2
   procedure(vegn_tiles_merge_check), optional :: vegn_merge_check

   ! make sure that the two tiles have the same components. For
   ! uniformity every component is checked, even though snow and
   ! cana are always present in current design
   answer = (associated(tile1%glac).eqv.associated(tile2%glac)).and. &
            (associated(tile1%lake).eqv.associated(tile2%lake)).and. &
            (associated(tile1%soil).eqv.associated(tile2%soil)).and. &
            (associated(tile1%snow).eqv.associated(tile2%snow)).and. &
            (associated(tile1%cana).eqv.associated(tile2%cana)).and. &
            (associated(tile1%vegn).eqv.associated(tile2%vegn))

   if (answer.and.associated(tile1%glac)) &
      answer = answer.and.glac_tiles_can_be_merged(tile1%glac,tile2%glac)
   if (answer.and.associated(tile1%lake)) &
      answer = answer.and.lake_tiles_can_be_merged(tile1%lake,tile2%lake)
   if (answer.and.associated(tile1%soil)) &
      answer = answer.and.soil_tiles_can_be_merged(tile1%soil,tile2%soil)
   if (answer.and.associated(tile1%cana)) &
      answer = answer.and.cana_tiles_can_be_merged(tile1%cana,tile2%cana)
   if (answer.and.associated(tile1%snow)) &
      answer = answer.and.snow_tiles_can_be_merged(tile1%snow,tile2%snow)
   if (answer.and.associated(tile1%vegn)) then
      if (present(vegn_merge_check)) then
          answer = answer.and.vegn_merge_check(tile1%vegn,tile2%vegn)
      else
          ! use traditional check: match between land use types and biomass bins
          answer = answer.and.vegn_tiles_can_be_merged(tile1%vegn,tile2%vegn)
      endif
   endif
end function land_tiles_can_be_merged

! ============================================================================
! merges the two tiles, putting merged state into the second tile. The first
! tile is unchanged
subroutine merge_land_tiles(tile1,tile2)
  type(land_tile_type), intent(in)    :: tile1
  type(land_tile_type), intent(inout) :: tile2

  ! ---- local vars
  real :: x1,x2
  real :: dheat

  if(associated(tile1%glac)) &
       call merge_glac_tiles(tile1%glac, tile1%frac, tile2%glac, tile2%frac)
  if(associated(tile1%lake)) &
       call merge_lake_tiles(tile1%lake, tile1%frac, tile2%lake, tile2%frac)
  if(associated(tile1%soil)) &
       call merge_soil_tiles(tile1%soil, tile1%frac, tile2%soil, tile2%frac)
  if(associated(tile1%soilc)) &
       call tile2%soilc%merge(tile2%frac, tile1%soilc, tile1%frac)

  if(associated(tile1%cana)) &
       call merge_cana_tiles(tile1%cana, tile1%frac, tile2%cana, tile2%frac)
  if(associated(tile1%snow)) &
       call tile2%snow%merge_snow_tiles(tile2%frac, tile1%snow, tile1%frac) ! EZSNOW

  dheat = 0.0
  if (associated(tile1%vegn)) then
     call merge_vegn_tiles(tile1%vegn, tile1%frac, tile2%vegn, tile2%frac, dheat)
     call kill_small_cohorts_ppa(tile2%vegn,tile2%soilc)
  endif

  ! calculate normalized weights
  x1 = tile1%frac/(tile1%frac+tile2%frac)
  x2 = 1.0 - x1

#define __MERGE__(field) tile2%field = x1*tile1%field + x2*tile2%field
  __MERGE__(lwup)
  __MERGE__(e_res_1)
  __MERGE__(e_res_2)
#undef __MERGE__

  tile2%e_res_2 = tile2%e_res_2 - dheat
  tile2%frac = tile1%frac + tile2%frac
end subroutine merge_land_tiles

! ==============================================================================
! given a pointer to a tile and a tile list, insert the tile into the list so that
! if tile can be merged with any one already present, it is merged; otherwise
! the tile is added to the list
subroutine merge_land_tile_into_list(tile, list)
  type(land_tile_type), pointer :: tile
  type(land_tile_list_type), intent(inout) :: list

  ! ---- local vars
  type(land_tile_type), pointer :: ptr
  type(land_tile_enum_type) :: ct

  ! try to find a tile that we can merge to
  ct = first_elmt(list)
  do while(loop_over_tiles(ct,ptr))
     if (land_tiles_can_be_merged(tile,ptr).and.ptr%frac>0.0) then
        call merge_land_tiles(tile,ptr)
        call delete_land_tile(tile)
        return ! break out of the subroutine
     endif
  enddo
  ! we reach here only if no suitable files was found in the list
  ! if no suitable tile was found, just insert given tile into the list.
  call insert(tile,list)
end subroutine merge_land_tile_into_list

! ============================================================================
! given tile list, tries to re-merge as many tiles as possible, to reduce
! computational burden.
subroutine remerge_tile_list(list)
  type(land_tile_list_type), intent(inout) :: list

  type(land_tile_type), pointer :: tile, tile1, tile2, dst
  type(land_tile_enum_type) :: ce, co
  type(land_tile_list_type) :: tmp, tmp1 ! temporary list to hold large and small tiles, respectively
  integer :: i
  real :: d, dmin ! "distance" between vegetation tiles in biomass
  real, parameter :: eps = 0.001 ! small number to make sure "distance" is not unreasonable
        ! when bwood is close to 0

  ! for conservation checks:
  real :: lmass0,fmass0,cmass0,nmass0,heat0
  real :: lmass1,fmass1,cmass1,nmass1,heat1
  real :: lmass,fmass,cmass,nmass,heat

  ! + conservation check part 1
  lmass0=0.0 ; fmass0=0.0 ; cmass0=0.0 ; nmass0=0.0 ; heat0=0.0
  ce=first_elmt(list)
  do while (loop_over_tiles(ce,tile))
     ! tile values, per unit area
     call get_tile_water(tile,lmass,fmass)
     cmass  = land_tile_carbon(tile)
     nmass  = land_tile_nitrogen(tile)
     heat   = land_tile_heat(tile)
     ! accumulate grid cell values
     lmass0 = lmass0 + lmass*tile%frac
     fmass0 = fmass0 + fmass*tile%frac
     cmass0 = cmass0 + cmass*tile%frac
     nmass0 = nmass0 + nmass*tile%frac
     heat0  =  heat0 +  heat*tile%frac
  enddo
  ! - conservation check part 1

  if (is_watch_cell()) then
     write (*,*)'##### remerge_tile_list input #####'
     ce = first_elmt(list); i = 1
     do while(loop_over_tiles(ce, tile))
        write(*,'(i3, 2x)',advance='no') i
        call print_land_tile_info(tile)
        i = i+1
     enddo
  endif

  call land_tile_list_init(tmp)
  call land_tile_list_init(tmp1)
  ! move all tiles into two temporary list: very small soil tiles stored in tmp1,
  ! while tiles with non-negligible land area fraction are merged into tmp
  do while (.not.empty(list))
     ce=first_elmt(list)
     tile=>current_tile(ce)
     call remove(ce)
     if (associated(tile%vegn).and.tile%frac < min_tile_frac) then
        call append_to_list(tile,tmp1)
     else
        ! this merges individual tiles, according to general criteria
        call merge_land_tile_into_list(tile,tmp)
     endif
  enddo

  ! merge all small soil/vegn tiles into larger tiles, using relaxed merge criteria
  do while (.not.empty(tmp1))
     ce=first_elmt(tmp1)
     tile1=>current_tile(ce)
     call remove(ce)
     ! select the best larger tile that tile1 can be merged into
     co = first_elmt(tmp); dmin = HUGE(1.0); dst=>NULL()
     do while (loop_over_tiles(co, tile2))
        ! the check below returns true if the tiles are compatible in all non-vegetation
        ! respects (e.g. soil type, etc.) and their land use types are the same, regardless
        ! of the vegetation state. This loop selects the tiles that are closest in bwood.
        if (land_tiles_can_be_merged(tile1,tile2,vegn_merge_check=vegn_tile_lu_match)) then
            ! this hard-coded rule can be replaced with a more sophisticated function,
            ! if desired
            d = (abs(vegn_tile_bwood(tile1%vegn))+eps)/ &
                (abs(vegn_tile_bwood(tile2%vegn))+eps)
            if (d<1) d = 1.0/d
            if (d<dmin) then
               dst=>tile2; dmin = d
            endif
        endif
     enddo
     if (associated(dst)) then
        call merge_land_tiles(tile1,dst)
        call delete_land_tile(tile1)
     else
        ! we get here only if there are no matching tiles, e.g. tiny cropland tile,
        ! which is the only cropland tile in the grid cell and therefore cannot be
        ! merged with anything.
        call append_to_list(tile1,tmp)
     endif
  enddo

  ! move all tiles from temporary list to the tile map
  do while (.not.empty(tmp))
     ce=first_elmt(tmp)
     tile=>current_tile(ce)
     call remove(ce)
     call append_to_list(tile,list)
  enddo
  call land_tile_list_end(tmp)
  call land_tile_list_end(tmp1)

  if (is_watch_cell()) then
     write (*,*)'##### remerge_tile_list output #####'
     ce = first_elmt(list); i = 1
     do while(loop_over_tiles(ce, tile))
        write(*,'(i3, 2x)',advance='no') i
        call print_land_tile_info(tile)
        i = i+1
     enddo
  endif

  ! + conservation check part 2
  lmass1=0.0 ; fmass1=0.0 ; cmass1=0.0 ; nmass1=0.0 ; heat1=0.0
  ce=first_elmt(list)
  do while (loop_over_tiles(ce,tile))
     ! tile values, per unit area
     call get_tile_water(tile,lmass,fmass)
     cmass  = land_tile_carbon(tile)
     nmass  = land_tile_nitrogen(tile)
     heat   = land_tile_heat(tile)
     ! accumulate grid cell values
     lmass1 = lmass1 + lmass*tile%frac
     fmass1 = fmass1 + fmass*tile%frac
     cmass1 = cmass1 + cmass*tile%frac
     nmass1 = nmass1 + nmass*tile%frac
     heat1  =  heat1 +  heat*tile%frac
  enddo
  call check_conservation ('remerge_tile_list', 'liquid water', lmass0, lmass1, water_cons_tol)
  call check_conservation ('remerge_tile_list', 'frozen water', fmass0, fmass1, water_cons_tol)
  call check_conservation ('remerge_tile_list', 'carbon'      , cmass0, cmass1, carbon_cons_tol)
  call check_conservation ('remerge_tile_list', 'nitrogen'    , nmass0, nmass1, nitrogen_cons_tol)
  call check_conservation ('remerge_tile_list', 'heat'        , heat0,  heat1,  heat_cons_tol)
  ! - conservation check part 2
end subroutine remerge_tile_list


! #### tile container ########################################################

! ============================================================================
! tile list constructor: initializes essential innards of tile collection
! for future use. In current implementation, it is safe to call this function
! on a tile list more then once
subroutine land_tile_list_init(list)
  type(land_tile_list_type), intent(inout) :: list

  if (.not.associated(list%head)) then
     allocate(list%head)
     list%head%prev=>list%head
     list%head%next=>list%head
  endif
end subroutine land_tile_list_init

! ============================================================================
! tile list destructor: destroys the list of tiles. NOTE that it also destroys
! all the tiles that are still in the list.
subroutine land_tile_list_end(list)
  type(land_tile_list_type), intent(inout) :: list

  if(associated(list%head)) then
     call erase(list)
     deallocate(list%head)
  endif
end subroutine land_tile_list_end

! ============================================================================
subroutine check_tile_list_inited(list)
  type(land_tile_list_type), intent(in) :: list

  if (.not.associated(list%head)) &
     call error_mesg('land_tile_mod','tile container was not initialized before use', FATAL)

end subroutine check_tile_list_inited


! ============================================================================
! returns true is the list is empty
function empty(list)
  logical empty
  type(land_tile_list_type), intent(in) :: list

  empty = .not.associated(list%head)
  if (.not.empty) &
       empty = associated(list%head%next,list%head)

end function empty

! ============================================================================
! returns the number of items currently stored in the list
function n_items_in_list(list) result (n)
  type(land_tile_list_type), intent(in) :: list
  integer :: n

  type(land_tile_list_node_type), pointer :: node

  n=0;
  if(.not.associated(list%head)) return

  node => list%head%next
  do while ( .not.(associated(node,list%head)) )
     n = n+1
     node => node%next
  enddo
end function n_items_in_list

! ============================================================================
function elmt_at_index(list,k) result(ptr)
  type(land_tile_list_type), intent(in) :: list ! list of tiles
  integer,                   intent(in) :: k    ! index
  type(land_tile_type), pointer :: ptr ! return value

  type(land_tile_enum_type) :: ct
  integer :: i
  ct = first_elmt(list); i = 1
  do while (loop_over_tiles(ct, ptr))
     if (i==k) exit ! from loop
     i = i+1
  enddo
end function elmt_at_index

! ============================================================================
subroutine append_to_list(tile,list)
  type(land_tile_type),            pointer :: tile
  type(land_tile_list_type), intent(inout) :: list

  call insert_at_position(tile,tail_elmt(list))
end subroutine append_to_list


! ============================================================================
subroutine remove_all_from_list(list)
  type(land_tile_list_type), intent(inout) :: list

  type(land_tile_enum_type) :: ce
  ce=first_elmt(list)
  do while(ce/=tail_elmt(list))
     call remove_at_position(ce)
  enddo
end subroutine remove_all_from_list


! ============================================================================
subroutine erase_all_from_list(list)
  type(land_tile_list_type), intent(inout) :: list

  type(land_tile_enum_type) :: ce
  ce=first_elmt(list)
  do while(ce/=tail_elmt(list))
     call erase_at_position(ce)
  enddo
end subroutine erase_all_from_list



! #### tile container enumerator #############################################

! ============================================================================
! returns enumerator pointing to the first element of the container
function land_tile_list_begin_0d(list) result(ce)
  type(land_tile_enum_type) :: ce  ! return value
  type(land_tile_list_type), intent(in) :: list

  call check_tile_list_inited(list)
  ce%node=>list%head%next
  ce%i = 1 ; ce%j = 1 ; ce%k = 1 ; ce%l = 1
end function land_tile_list_begin_0d


! ============================================================================
! returns enumerator pointing to the first element of the 2D array of
! containers
function land_tile_list_begin_1d(tiles, ls) result(ce)
  type(land_tile_enum_type) :: ce  ! return value
  type(land_tile_list_type), intent(in), target :: tiles(:)
  integer, intent(in), optional :: ls ! origin of the array

  integer :: l

  ! list up pointer to the array of containers
  ce%tiles=>tiles

  ! initialize offsets of indices
  ce%lo = 0
  if(present(ls)) ce%lo = ls-1

  ! initialize current position in the array of containers -- find
  ! first non-empty container and list the pointer to the current
  ! container node
  ce%k = 1
  do l = 1,size(tiles(:))
     call check_tile_list_inited(tiles(l))
     ce%node => tiles(l)%head%next
     ce%l = l
     ce%i = lnd%i_index(l+lnd%ls-1)
     ce%j = lnd%j_index(l+lnd%ls-1)
     if(associated(ce%node%data)) return
  enddo
end function land_tile_list_begin_1d


! ============================================================================
! returns enumerator pointing to the end of container: actually the next element
! behind the last element of the container
function land_tile_list_end_0d(list) result (ce)
  type(land_tile_enum_type) :: ce ! return value
  type(land_tile_list_type), intent(in) :: list

  call check_tile_list_inited(list)
  ce%node=>list%head
  ce%i = 1 ; ce%j = 1 ; ce%l=1 ; ce%k = nitems(list)+1
end function land_tile_list_end_0d


! ============================================================================
! returns enumerator pointing to the end of 2D array of containers: actually
! the next element behind the last element of the last container
function land_tile_list_end_1d(tiles) result (ce)
  type(land_tile_enum_type) :: ce ! return value
  type(land_tile_list_type), intent(in), target :: tiles(:)

  ! list up pointer to the array of containers
  ce%tiles=>tiles

  ! initialize offsets of indices
  ce%lo = 0

  ! initialize current position in the array of containers
  ce%l = ubound(tiles,1)
  ce%i = lnd%i_index(ce%l+lnd%ls-1)
  ce%j = lnd%j_index(ce%l+lnd%ls-1)
  ce%k = nitems(tiles(ce%l))+1

  ! list the pointer to the current tile
  call check_tile_list_inited(tiles(ce%l))
  ce%node=>tiles(ce%l)%head

end function land_tile_list_end_1d


! ============================================================================
! returns enumerator pointing to the next element of the container.
function next_elmt(pos0) result(ce)
  type(land_tile_enum_type) :: ce ! return value
  type(land_tile_enum_type), intent(in) :: pos0

  integer :: le

  ce = pos0
  ce%node => ce%node%next ; ce%k = ce%k+1
  if(associated(ce%tiles)) then
     le = ubound(ce%tiles,1)
     do while(.not.associated(ce%node%data))
        ce%k = 1; ! reset tile index
        if(ce%l<le)then
           ce%l = ce%l+1
           ce%i = lnd%i_index(ce%l+lnd%ls-1)
           ce%j = lnd%j_index(ce%l+lnd%ls-1)
        else
           return
        endif
        call check_tile_list_inited(ce%tiles(ce%l))
        ce%node => ce%tiles(ce%l)%head%next
     enddo
  endif
end function next_elmt


! ============================================================================
! returns enumerator pointing to the previous element of the container.
function prev_elmt(pos0) result(ce)
  type(land_tile_enum_type) :: ce ! return value
  type(land_tile_enum_type), intent(in) :: pos0

  integer :: ls

  ce = pos0
  ce%node => ce%node%prev ; ce%k = ce%k - 1
  if(associated(ce%tiles)) then
     ls = lbound(ce%tiles,1)
     do while(.not.associated(ce%node%data))
        ce%k = 1; ! reset tile index
        if(ce%l>ls)then
           ce%l = ce%l - 1
           ce%i = lnd%i_index(ce%l+lnd%ls-1)
           ce%j = lnd%j_index(ce%l+lnd%ls-1)
        else
           return
        endif
        call check_tile_list_inited(ce%tiles(ce%l))
        ce%node => ce%tiles(ce%l)%head%prev
        ce%k    =  nitems(ce%tiles(ce%l))
     enddo
  endif

end function prev_elmt

! ============================================================================
! returns TRUE if both enums refer to the same list node (and, therefore, tile)
! or if both do not refer to anything.
function enums_are_equal(pos1,pos2) result(ret)
  logical :: ret ! return value
  type(land_tile_enum_type), intent(in) :: pos1,pos2

  if(associated(pos1%node)) then
     ret = associated(pos1%node,pos2%node)
  else
     ret = .not.associated(pos2%node)
  endif
end function enums_are_equal

! ============================================================================
! returns TRUE if two enumerators are not equal
function enums_are_not_equal(pos1,pos2) result(ret)
  logical :: ret ! return value
  type(land_tile_enum_type), intent(in) :: pos1,pos2

  ret=.not.enums_are_equal(pos1,pos2)
end function enums_are_not_equal

! ============================================================================
! returns pointer to the tile currently addressed by the enumerator
function current_tile(ce) result(ptr)
  type(land_tile_type), pointer :: ptr ! return value
  type(land_tile_enum_type), intent(in) :: ce

  ptr => ce%node%data
end function current_tile

! ============================================================================
! returns indices corresponding to the enumerator; for enumerator associated
! with a single tile list (not with 2D array of lists) returned i and j are
! equal to 1
subroutine get_elmt_indices(ce,i,j,k,l)
  type(land_tile_enum_type), intent(in) :: ce
  integer, intent(out), optional :: i, j, l, k

  if (present(i)) i = lnd%i_index(ce%l+lnd%ls-1)
  if (present(j)) j = lnd%j_index(ce%l+lnd%ls-1)
  if (present(k)) k = ce%k
  if (present(l)) l = ce%l+ce%lo

end subroutine get_elmt_indices

! ============================================================================
! given an enumerator, sets tile pointer to the current tile and its indices, and
! attempts to advance enumerator to the next tile. If enumerator was already at
! the end of the tile list, returns FALSE; in this case pointer "tile" and
! indices i,j,k are not defined.
function loop_over_tiles(ce, tile, l, k, i, j) result(R); logical R
  type(land_tile_enum_type), intent(inout) :: ce
  type(land_tile_type)     , pointer, optional :: tile
  integer, intent(out), optional :: i,j,l,k ! indices of the tile

  type(land_tile_type), pointer :: tile_

  tile_=>current_tile(ce)
  if (present(tile)) tile=>tile_
  call get_elmt_indices(ce,i=i,j=j,l=l,k=k)
  ! advance enumerator to the next element
  ce = next_elmt(ce)
  R  = associated(tile_)
end function loop_over_tiles

! ============================================================================
! inserts tile at the position indicated by enumerator: in fact right in front
! of it.
subroutine insert_at_position(tile,ce)
  type(land_tile_type),         pointer :: tile
  type(land_tile_enum_type), intent(in) :: ce

  ! local vars
  type(land_tile_list_node_type), pointer :: node,n,p

  allocate(node)
  node%data=>tile

  n=>ce%node  ; p=>n%prev

  node%next=>n ; node%prev=>p
  n%prev=>node ; p%next=>node

end subroutine insert_at_position

! ============================================================================
subroutine remove_at_position(enum)
  type(land_tile_enum_type), intent(inout) :: enum

  type(land_tile_list_node_type),pointer :: n,p
  type(land_tile_enum_type) :: next

  if(.not.associated(enum%node)) &
     call error_mesg('remove_at_position','attempt to remove tail element of a list', FATAL)

  next = next_elmt(enum)

  n => enum%node%next
  p => enum%node%prev

  n%prev=>p ; p%next=>n
  deallocate(enum%node)

  enum=next
  if(enum%k>1) enum%k = enum%k-1

end subroutine remove_at_position

! ============================================================================
subroutine erase_at_position(ce)
  type(land_tile_enum_type), intent(inout) :: ce

  type(land_tile_type), pointer :: tile

  tile=>current_tile(ce)
  call remove_at_position(ce)
  call delete_land_tile(tile)

end subroutine erase_at_position


! ============================================================================
function tile_is_selected(tile, sel)
! returns true if the tile fits specified selector
  logical :: tile_is_selected
  type(land_tile_type)    , intent(in) :: tile
  type(tile_selector_type), intent(in) :: sel

  tile_is_selected = .FALSE.
  select case(sel%tag)
  case(SEL_SOIL)
     if(associated(tile%soil)) &
          tile_is_selected = soil_is_selected(tile%soil,sel)
  case(SEL_VEGN)
     if(associated(tile%vegn)) &
          tile_is_selected = vegn_is_selected(tile%vegn,sel)
  case(SEL_LAKE)
     if(associated(tile%lake)) &
          tile_is_selected = lake_is_selected(tile%lake,sel)
  case(SEL_GLAC)
     if(associated(tile%glac)) &
          tile_is_selected = glac_is_selected(tile%glac,sel)
  case(SEL_SNOW)
     if(associated(tile%snow)) &
          tile_is_selected = tile%snow%snow_is_selected(sel) ! EZSNOW
  case(SEL_CANA)
     if(associated(tile%cana)) &
          tile_is_selected = cana_is_selected(tile%cana,sel)
  case(SEL_HLSP)
     if(associated(tile%soil)) &
          tile_is_selected = hlsp_is_selected(tile%soil,sel)
  case default
     tile_is_selected=.true.
  end select

end function tile_is_selected


! ============================================================================
subroutine print_land_tile_info(tile)
  type(land_tile_type), intent(in) :: tile

  write(*,'("(tag =",i3,", frac =",g23.16)',advance='no') tile%tag, tile%frac
  if(associated(tile%lake)) write(*,'(a,i3)',advance='no')', lake =',tile%lake%tag
  if(associated(tile%soil)) write(*,'(a,i3)',advance='no')', soil =',tile%soil%tag
  if(associated(tile%glac)) write(*,'(a,i3)',advance='no')', glac =',tile%glac%tag
!  if(associated(tile%snow)) write(*,'(a,i3)',advance='no')', snow =',tile%snow%tag
!  if(associated(tile%cana)) write(*,'(a)',advance='no')', cana'
  if(associated(tile%vegn)) then
       write(*,'(a)',advance='no')', vegn LU = '//landuse_name(tile%vegn%landuse)
       write(*,'(a,g23.16)',advance='no') ', bwood = ',vegn_tile_bwood(tile%vegn)
  endif
  write(*,'(")")')

end subroutine print_land_tile_info

end module land_tile_mod
