!***********************************************************************
!*                   GNU Lesser General Public License
!*
!* This file is part of the GFDL Land Model 4 (LM4).
!*
!* LM4 is free software: you can redistribute it and/or modify it under
!* the terms of the GNU Lesser General Public License as published by
!* the Free Software Foundation, either version 3 of the License, or (at
!* your option) any later version.
!*
!* LM4 is distributed in the hope that it will be useful, but WITHOUT
!* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
!* FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
!* for more details.
!*
!* You should have received a copy of the GNU Lesser General Public
!* License along with LM4.  If not, see <http://www.gnu.org/licenses/>.
!***********************************************************************
module land_tracers_mod

use fms_mod           , only : mpp_npes, string, error_mesg, FATAL, NOTE, stdout
use field_manager_mod , only : MODEL_LAND
use tracer_manager_mod, only : register_tracers, get_number_tracers, &
   get_tracer_index, get_tracer_names, NO_TRACER
use land_data_mod     , only : log_version
implicit none
private

! ==== public interfaces =====================================================
public :: land_tracers_init, land_tracers_end

integer, protected, public :: ntcana ! number of prognostic land tracers in canopy air
integer, protected, public :: isphum, ico2 ! indices of specific humidity and CO2
! ==== end of public interfaces ==============================================

! ---- module constants ------------------------------------------------------
character(len=*), parameter :: module_name = 'land_tracers_mod'
#include "../shared/version_variable.inc"

! ---- private module variables ----------------------------------------------
logical :: module_is_initialized = .FALSE.

contains ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

! ============================================================================
subroutine land_tracers_init()
  integer :: ntracers, ndiag, i, tr
  character(32) :: name

  ! write the version and tag name to the logfile
  call log_version(version, module_name, &
  __FILE__)

  call register_tracers ( MODEL_LAND, ntracers, ntcana, ndiag )

  do i = 1,ntcana
     call get_tracer_names(MODEL_LAND,i,name)
     call error_mesg('land_tracers_init','land tracer ('//trim(string(i)) &
                     //')= "'//trim(name)//'"', NOTE)
  enddo

  ! check that required tracers are present
  isphum = get_tracer_index ( MODEL_LAND, 'sphum' )
  if (isphum==NO_TRACER) then
     call error_mesg('land_model_init','required land tracer "sphum" not found',FATAL)
  endif
  ico2 = get_tracer_index ( MODEL_LAND, 'co2' )
  if (ico2==NO_TRACER) then
     call error_mesg('land_model_init','required land tracer "co2" not found',FATAL)
  endif
  call error_mesg('land_tracers_init','isphum='//string(isphum),NOTE)
  call error_mesg('land_tracers_init','ico2='//string(ico2),NOTE)

  module_is_initialized = .TRUE.
end subroutine land_tracers_init

! ============================================================================
subroutine land_tracers_end()
   module_is_initialized = .FALSE.
end subroutine land_tracers_end

end module land_tracers_mod
