module predefined_tiles_mod

use fms_mod, only : input_nml_file, check_nml_error, stdlog, error_mesg, FATAL, &
      mpp_pe, mpp_root_pe
use land_data_mod, only : log_version

implicit none
private

! ==== public interfaces =====================================================
public :: read_predefined_tiles_namelist
! ==== end of public interfaces ==============================================

! ==== module constants ======================================================
character(len=*), parameter :: module_name = 'predefined_tiles_mod'
#include "../shared/version_variable.inc"

!---- namelist ---------------------------------------------------------------
logical, public, protected :: use_predefined_tiles = .FALSE. ! If true, the tiles for
              ! each grid cell and their parameters are read from predefined tiles
              ! database. Otherwise, all other  namelist parameters defined here are
              ! irrelevant.
logical, public, protected :: use_predefined_biomass = .FALSE. ! if true, initialize
              ! plant biomass from predefined tiles database
logical, public, protected :: use_predefined_landuse = .TRUE. ! if true, land use
              ! information for each tile comes from the
logical, public, protected :: downscale_surface_meteorology = .FALSE. ! If true, the
              ! downscaling weights in the predefined tiles database will be used to
              ! downscale prec and sw (could be eventually be used for others as well...)

namelist /predefined_tiles_nml/ use_predefined_tiles, use_predefined_biomass, &
  use_predefined_landuse, downscale_surface_meteorology

contains ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

subroutine read_predefined_tiles_namelist()
  ! ---- local vars
  integer :: unit         ! unit for namelist i/o
  integer :: io           ! i/o status for the namelist
  integer :: ierr         ! error code, returned by i/o routines

  call log_version(version, module_name, &
  __FILE__)

  read (input_nml_file, nml=predefined_tiles_nml, iostat=io)
  ierr = check_nml_error(io, 'predefined_tiles_nml')

  if (mpp_pe() == mpp_root_pe()) then
     unit=stdlog()
     write(unit, nml=predefined_tiles_nml)
  endif
end subroutine read_predefined_tiles_namelist

end module predefined_tiles_mod
