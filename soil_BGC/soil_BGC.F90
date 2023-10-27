module soil_BGC_mod

use land_data_mod, only : log_version
use soil_BGC_base_mod,   only : soil_BGC_option, SOIL_BGC_SIMPLE, SOIL_BGC_CORPSE, soil_BGC_GIMICS
use soil_BGC_SIMPLE_mod, only : soil_BGC_init_SIMPLE, soil_BGC_save_restart_SIMPLE
use soil_BGC_CORPSE_mod, only : soil_BGC_init_CORPSE, soil_BGC_save_restart_CORPSE
use soil_BGC_GIMICS_mod, only : soil_BGC_init_GIMICS, soil_BGC_save_restart_GIMICS

implicit none; private

public :: soil_BGC_init
public :: save_soil_BGC_restart

! ==== module constants ======================================================
character(len=*), parameter :: module_name = 'soil_BGC_mod'
#include "../shared/version_variable.inc"

contains ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

! ============================================================================
subroutine soil_BGC_init( id_ug, id_zfull )
  integer,intent(in) :: id_ug    !< Unstructured axis id
  integer,intent(in) :: id_zfull !< Vertical (depth) axis id

  call log_version(version, module_name, __FILE__)

  select case (soil_BGC_option)
  case(SOIL_BGC_SIMPLE)
    call soil_BGC_init_SIMPLE( id_ug, id_zfull )
  case(SOIL_BGC_CORPSE)
    call soil_BGC_init_CORPSE( id_ug, id_zfull )
  case(SOIL_BGC_GIMICS)
    call soil_BGC_init_GIMICS( id_ug, id_zfull )
  end select
end subroutine

! ============================================================================
subroutine save_soil_BGC_restart(tile_dim_length, timestamp)
  integer, intent(in) :: tile_dim_length ! length of tile dim. in the output file
  character(*), intent(in) :: timestamp ! timestamp to add to the file name

  select case (soil_BGC_option)
  case(SOIL_BGC_SIMPLE)
    call soil_BGC_save_restart_SIMPLE(tile_dim_length, timestamp)
  case(SOIL_BGC_CORPSE)
    call soil_BGC_save_restart_CORPSE(tile_dim_length, timestamp)
  case(SOIL_BGC_GIMICS)
    call soil_BGC_save_restart_GIMICS(tile_dim_length, timestamp)
  end select
end subroutine

end module