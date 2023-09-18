module soilc_restart_mod

use soilc_mod,        only : soil_carbon_option, SOIL_BGC_SIMPLE, SOIL_BGC_CORPSE
use soil_BGC_SIMPLE_mod, only : soil_BGC_init_SIMPLE, soil_BGC_save_restart_SIMPLE
use soilc_CORPSE_mod, only : soilc_init_CORPSE, save_soilc_CORPSE_restart

implicit none; private

! ==== public interfaces =====================================================
public :: soilc_init
public :: save_soilc_restart

contains ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

! ============================================================================
subroutine soilc_init( id_ug, id_zfull )
  integer,intent(in) :: id_ug    !< Unstructured axis id
  integer,intent(in) :: id_zfull !< Vertical (depth) axis id

  select case (soil_carbon_option)
  case(SOIL_BGC_SIMPLE)
    call soil_BGC_init_SIMPLE( id_ug, id_zfull )
  case(SOIL_BGC_CORPSE)
    call soilc_init_CORPSE( id_ug, id_zfull )
  end select
end subroutine

! ============================================================================
subroutine save_soilc_restart(tile_dim_length, timestamp)
  integer, intent(in) :: tile_dim_length ! length of tile dim. in the output file
  character(*), intent(in) :: timestamp ! timestamp to add to the file name

  select case (soil_carbon_option)
  case(SOIL_BGC_SIMPLE)
    call soil_BGC_save_restart_SIMPLE(tile_dim_length, timestamp)
  case(SOIL_BGC_CORPSE)
    call save_soilc_CORPSE_restart(tile_dim_length, timestamp)
  end select
end subroutine

end module