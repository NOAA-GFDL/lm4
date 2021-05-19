module transitions_input_mod
#include <fms_platform.h>

#ifdef INTERNAL_FILE_NML
use mpp_mod, only: input_nml_file
#else
use fms_mod, only: open_namelist_file
#endif
! TODO: clean up unused imports
use fms_mod, only : string, error_mesg, FATAL, WARNING, NOTE, &
     mpp_pe, lowercase, file_exist, close_file, read_data, &
     check_nml_error, stdlog, mpp_root_pe, fms_error_handler


use land_data_mod, only : lnd, log_version, horiz_interp_ug

implicit none
private

! ==== public subroutines ===================================================
public :: read_transitions_namelist
! ==== end of public subroutines ============================================

! ==== module constants =====================================================
character(len=*), parameter :: module_name = 'transitions_input_mod'
#include "../shared/version_variable.inc"

! ---- namelist variables ---------------------------------------------------
logical, public, protected :: do_landuse_change = .FALSE. ! if true, then the landuse changes with time
character(len=1024),public, protected :: input_file  = '' ! input data set of transition dates
character(len=1024),public, protected :: state_file  = '' ! input data set of LU states (for initial transition only)
character(len=1024),public, protected :: static_file = '' ! static data file, for input land fraction
character(len=16),public, protected  :: data_type  = 'luh1' ! or 'luh2'
! distribute_transitions sets how the land use transitions are distributed among
! tiles within grid cells. 'lm3' is traditional (transitions applied to every
! tile in equal measure, except secondary-to-secondary); 'min-tiles' applies
! transitions to tiles in the order of priority, thereby minimizing the number
! of resulting tiles
logical, public, protected :: rangeland_is_pasture = .FALSE. ! if true, rangeland is combined with pastures.
! This only applies to luh2 transitions, since there is no rangeland in luh1 anyway.
character(len=16), public, protected  :: distribute_transitions  = 'lm3' ! or 'min-n-tiles'
! sets how to handle transition overshoot: that is, the situation when transition
! is larger than available area of the given land use type.
character(len=16), public, protected :: overshoot_handling = 'report' ! or 'stop', or 'ignore'
real, public, protected :: overshoot_tolerance = 1e-4 ! tolerance interval for overshoots
! specifies how to handle non-conservation
character(len=16), public, protected :: conservation_handling = 'stop' ! or 'report', or 'ignore'

! for irrigation
character(len=1024), public, protected :: irrigation_file = '' ! irrigation data file
logical, public, protected :: do_irrigation = .FALSE.

character(len=1024), public, protected :: input_file_lake  = '' ! input data set of lake transition dates
character(len=1024), public, protected :: state_file_lake  = '' ! input data set of LU states (for initial transition only)
character(len=1024), public, protected :: depth_file_rsv   = '' ! reservoir construction depth
logical, protected, public :: do_lake_change = .FALSE.

namelist/landuse_nml/do_landuse_change, do_irrigation, data_type, &
     input_file, state_file, static_file, irrigation_file,&
     rangeland_is_pasture, distribute_transitions, &
     overshoot_handling, overshoot_tolerance, &
     conservation_handling, &
     input_file_lake, state_file_lake, depth_file_rsv, do_lake_change

! ==== module data ==========================================================
logical :: module_is_initialized = .FALSE.

contains ! ###################################################################

subroutine read_transitions_namelist ()
  integer :: unit, ierr, io

  if(module_is_initialized) return
  module_is_initialized = .TRUE.
  call log_version(version, module_name, &
  __FILE__)

#ifdef INTERNAL_FILE_NML
  read (input_nml_file, nml=landuse_nml, iostat=io)
  ierr = check_nml_error(io, 'landuse_nml')
#else
  if (file_exist('input.nml')) then
     unit = open_namelist_file ( )
     ierr = 1;
     do while (ierr /= 0)
        read (unit, nml=landuse_nml, iostat=io, end=10)
        ierr = check_nml_error (io, 'landuse_nml')
     enddo
10   continue
     call close_file (unit)
  endif
#endif

  if (mpp_pe() == mpp_root_pe()) then
     unit=stdlog()
     write(unit, nml=landuse_nml)
  endif

  ! TODO: perhaps move namelist options parsing here? -- make respective namelist
  !       members private, but options public?
end subroutine read_transitions_namelist

end module transitions_input_mod