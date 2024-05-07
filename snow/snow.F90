! ============================================================================
! snow model module
! ============================================================================
module snow_mod

#include "../shared/debug.inc"

use fms_mod, only : error_mesg, FATAL, NOTE, lowercase

use land_constants_mod, only : NBANDS
use land_data_mod, only : log_version
use land_debug_mod, only : is_watch_point, land_error_message

use snow_tile_mod, only : snow_tile_type
use cm_snow_tile_mod, only: cm_snow_tile_type
use gl_snow_tile_mod, only: gl_snow_tile_type
use cm_snow_mod, only: cm_read_snow_namelist, cm_snow_init, cm_snow_end, &
    cm_save_snow_restart
use gl_snow_mod, only: gl_read_snow_namelist, gl_snow_init, gl_snow_end, &
    gl_save_snow_restart
use snow_evolution_mod, only: gl_compute_snow_albedo, albedo_to_use
use snow_base_mod, only : read_snow_model_namelist, snow_option, SNOW_CM, SNOW_GL

implicit none
private

! ==== public interfaces =====================================================
public :: read_snow_namelist
public :: snow_init
public :: snow_end
public :: save_snow_restart

! re-export snow model selector
public :: snow_option, SNOW_CM, SNOW_GL

! ==== module constants ======================================================
character(len=*), parameter :: module_name = 'snow_mod'
#include "../shared/version_variable.inc"

contains


subroutine read_snow_namelist()
  call log_version(version, module_name, &
  __FILE__)

  call read_snow_model_namelist()

  select case(snow_option)
  case(SNOW_CM)
     call cm_read_snow_namelist()
  case(SNOW_GL)
     call gl_read_snow_namelist()
  case default
     call land_error_message('read_snow_namelist: The value of snow_option is invalid. This should never happen. See developer', FATAL)
  end select
end subroutine read_snow_namelist


! initialize snow model
subroutine snow_init(id_ug)
  integer,intent(in) :: id_ug !< ID of "unstructured" diagnostic axis: it encodes
                              !! latitude and longtude

  select case(snow_option)
  case(SNOW_CM)
     call cm_snow_init(id_ug)
  case(SNOW_GL)
     call gl_snow_init(id_ug)
  case default
     call land_error_message('snow_init: The value of snow_option is invalid. This should never happen. See developer', FATAL)
  end select
end subroutine snow_init


! initialize snow model
subroutine snow_end()
  select case(snow_option)
  case(SNOW_CM)
     call cm_snow_end()
  case(SNOW_GL)
     call gl_snow_end()
  case default
     call land_error_message('snow_end: The value of snow_option is invalid. This should never happen. See developer', FATAL)
  end select
end subroutine snow_end


! save snow model restart file
subroutine save_snow_restart(tile_dim_length, timestamp)
  integer, intent(in) :: tile_dim_length ! length of tile dim. in the output file
  character(*), intent(in) :: timestamp ! timestamp to add to the file name

  select case(snow_option)
  case(SNOW_CM)
     call cm_save_snow_restart(tile_dim_length, timestamp)
  case(SNOW_GL)
     call gl_save_snow_restart(tile_dim_length, timestamp)
  case default
     call land_error_message('snow_end: The value of snow_option is invalid. This should never happen. See developer', FATAL)
  end select
end subroutine save_snow_restart


end module snow_mod
