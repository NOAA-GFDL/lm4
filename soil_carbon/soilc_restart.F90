module soilc_restart_mod

use soil_carbon_mod,  only : soil_carbon_option, &
    SOILC_CENTURY, SOILC_CENTURY_BY_LAYER, &
    SOILC_CORPSE, SOILC_CORPSE_N
use soilc_CENT_mod,   only : read_soilc_CENT_restart,   save_soilc_CENT_restart
use soilc_CORPSE_mod, only : read_soilc_CORPSE_restart, save_soilc_CORPSE_restart

implicit none; private

! ==== public interfaces =====================================================
public :: read_soilc_restart
public :: save_soilc_restart

contains ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

! ============================================================================
subroutine read_soilc_restart()
  select case (soil_carbon_option)
  case(SOILC_CENTURY, SOILC_CENTURY_BY_LAYER)
    call read_soilc_CENT_restart()
  case(SOILC_CORPSE, SOILC_CORPSE_N)
    call read_soilc_CORPSE_restart()
  end select
end subroutine

! ============================================================================
subroutine save_soilc_restart(tile_dim_length, timestamp)
  integer, intent(in) :: tile_dim_length ! length of tile dim. in the output file
  character(*), intent(in) :: timestamp ! timestamp to add to the file name

  select case (soil_carbon_option)
  case(SOILC_CENTURY, SOILC_CENTURY_BY_LAYER)
    call save_soilc_CENT_restart(tile_dim_length, timestamp)
  case(SOILC_CORPSE, SOILC_CORPSE_N)
    call save_soilc_CORPSE_restart(tile_dim_length, timestamp)
  end select
end subroutine

end module