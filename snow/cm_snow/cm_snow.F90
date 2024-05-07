! ============================================================================
! snow model module
! ============================================================================
module cm_snow_mod

#include "../../shared/debug.inc"

use mpp_mod, only: input_nml_file

use fms_mod, only : error_mesg, check_nml_error, &
     stdlog, mpp_pe, mpp_root_pe, FATAL, NOTE
use constants_mod,      only: tfreeze, hlv, hlf, PI
use land_tile_mod, only : land_tile_map, land_tile_type, land_tile_enum_type, &
     first_elmt, loop_over_tiles
use land_data_mod, only : lnd, log_version
use land_tile_io_mod, only: land_restart_type, &
     init_land_restart, open_land_restart, save_land_restart, free_land_restart, &
     add_restart_axis, add_tile_data, get_tile_data
use land_debug_mod, only : is_watch_point

use cm_snow_tile_mod, only : cm_snow_tile_type, read_snow_cm_namelist, &
     ! namelist variables:
     snow_density, retro_heat_capacity, albedo_to_use, init_temp, &
     init_pack_wl, init_pack_ws
use snow_tile_mod, only : read_snow_data_namelist, &
     max_lev, use_brdf
use snowpack_mod, only : clw, csw, read_snowpack_namelist


implicit none
private

! ==== public interfaces =====================================================
public :: cm_read_snow_namelist
public :: cm_snow_init
public :: cm_snow_end
public :: cm_save_snow_restart
! =====end of public interfaces ==============================================


! ==== module constants ======================================================
character(len=*), parameter :: module_name = 'cm_snow_mod'
#include "../../shared/version_variable.inc"

! ==== module variables ======================================================

logical         :: module_is_initialized =.FALSE.
integer         :: num_l    ! # of snow layers
! next three 'z' variables are all normalized by total snow pack depth
real            :: dz (max_lev) ! relative thicknesses of layers
real            :: z  (max_lev) ! relative depths of layer bounds
real            :: zz (max_lev) ! relative depths of layer centers
real            :: mc_fict

! ==== end of module variables ===============================================

contains


! ============================================================================
subroutine cm_read_snow_namelist()
  ! ---- local vars
  integer :: unit         ! unit for namelist i/o
  integer :: io           ! i/o status for the namelist
  integer :: ierr         ! error code, returned by i/o routines
  integer :: l            ! layer iterator

  call read_snow_data_namelist(num_l,dz,mc_fict)
  call read_snowpack_namelist()  ! need to read some variables from snowpack module
  call read_snow_cm_namelist()

  call log_version(version, module_name, &
  __FILE__)

  ! set up vertical discretization
  zz(1) = 0
  do l = 1, num_l
     zz(l+1) = zz(l) + dz(l)
     z(l)    = 0.5*(zz(l+1) + zz(l))
  enddo

end subroutine cm_read_snow_namelist

! ============================================================================
! initialize snow model
subroutine cm_snow_init(id_ug)
  integer,intent(in) :: id_ug    !< Unstructured axis id. Currently unused, but
            !! can be used in the future to register model-specific diagnostics

  ! ---- local vars ----------------------------------------------------------
  integer :: k
  type(land_tile_enum_type)     :: ce    ! tile list enumerator
  type(land_tile_type), pointer :: tile  ! pointer to current tile
!   character(*), parameter :: restart_file_name='INPUT/snow.res.nc' ! OLD VERSION
  character(*), parameter :: restart_file_name='INPUT/snow.nc' ! EZSNOW-2022SC
  type(land_restart_type) :: restart
  logical :: restart_exists

  module_is_initialized = .TRUE.

  ! -------- initialize snow state --------
  call open_land_restart(restart,restart_file_name,restart_exists)
  if (restart_exists) then
     call error_mesg('snow_init', 'reading NetCDF restart "'//trim(restart_file_name)//'"', NOTE)
     call get_tile_data(restart, 'temp', 'zfull', snow_temp_ptr)
     call get_tile_data(restart, 'wl'  , 'zfull', snow_wl_ptr)
     call get_tile_data(restart, 'ws'  , 'zfull', snow_ws_ptr)
  else
     call error_mesg('snow_init', 'cold-starting snow', NOTE)
     ce = first_elmt(land_tile_map)
     do while(loop_over_tiles(ce, tile))
        if (.not.associated(tile%snow)) cycle
        do k = 1,num_l
           call tile%snow%set_wli(k,init_pack_wl * dz(k) ) ! EZSNOW
           call tile%snow%set_wsi(k,init_pack_ws * dz(k) )
           call tile%snow%set_ti(k,init_temp)
        enddo
     enddo
  endif
  call free_land_restart(restart)

  if (trim(albedo_to_use)=='') then
     use_brdf = .false.
  elseif (trim(albedo_to_use)=='brdf-params') then
     use_brdf = .true.
  else
     call error_mesg('snow_init',&
          'option albedo_to_use="'//&
          trim(albedo_to_use)//'" is invalid, use "" or "brdf-params"',&
          FATAL)
  endif

end subroutine cm_snow_init


! ============================================================================
subroutine cm_snow_end ()

  module_is_initialized =.FALSE.

end subroutine cm_snow_end


! ============================================================================
subroutine cm_save_snow_restart (tile_dim_length, timestamp)
  integer, intent(in) :: tile_dim_length ! length of tile dim. in the output file
  character(*), intent(in) :: timestamp ! timestamp to add to the file name

  ! ---- local vars
  character(267) :: filename
  type(land_restart_type) :: restart ! restart file i/o object

  call error_mesg('snow_end','writing NetCDF restart',NOTE)
! Note that filename is updated for tile & rank numbers during file creation
  ! filename = trim(timestamp)//'snow.res.nc' ! OLD VERSION
  filename = 'RESTART/'//trim(timestamp)//'snow.nc' ! EZSNOW-2022SC
  call init_land_restart(restart, filename, snow_tile_exists, tile_dim_length)
  ! call add_restart_axis(restart,'zfull',zz(1:num_l), 'Z',longname='depth of level centers',sense=-1) ! OLD VERSION
  call add_restart_axis(restart,'zfull',zz(1:num_l),.false., 'Z',longname='depth of level centers',sense=-1) ! EZSNOW-2022SC

  call add_tile_data(restart,'temp','zfull', snow_temp_ptr, 'snow temperature','degrees_K')
  call add_tile_data(restart,'wl'  ,'zfull', snow_wl_ptr,   'snow liquid water content','kg/m2')
  call add_tile_data(restart,'ws'  ,'zfull', snow_ws_ptr,   'snow solid water content','kg/m2')

  call save_land_restart(restart)
  call free_land_restart(restart)

end subroutine cm_save_snow_restart


! ============================================================================
! tile existence detector: returns a logical value indicating wether component
! model tile exists or not
logical function snow_tile_exists(tile)
   type(land_tile_type), pointer :: tile
   snow_tile_exists = associated(tile%snow)
end function snow_tile_exists

! ============================================================================
! accessor functions: given a pointer to a land tile, they return pointer
! to the desired member of the land tile, of NULL if this member does not
! exist.
subroutine snow_temp_ptr(tile, i, ptr)
   type(land_tile_type), pointer :: tile ! input
   integer             , intent(in) :: i ! index in the array
   real                , pointer :: ptr  ! returned pointer to the data
   ptr=>NULL()
   if(associated(tile)) then
      if(associated(tile%snow)) ptr => tile%snow%T(i)
   endif
end subroutine snow_temp_ptr

subroutine snow_wl_ptr(tile, i, ptr)
   type(land_tile_type), pointer :: tile ! input
   integer             , intent(in) :: i ! index in the array
   real                , pointer :: ptr  ! returned pointer to the data
   ptr=>NULL()
   if(associated(tile)) then
      if(associated(tile%snow)) ptr => tile%snow%wl(i)
   endif
end subroutine snow_wl_ptr

subroutine snow_ws_ptr(tile, i, ptr)
   type(land_tile_type), pointer :: tile ! input
   integer             , intent(in) :: i ! index in the array
   real                , pointer :: ptr  ! returned pointer to the data
   ptr=>NULL()
   if(associated(tile)) then
      if(associated(tile%snow)) ptr => tile%snow%ws(i)
   endif
end subroutine snow_ws_ptr


end module cm_snow_mod



