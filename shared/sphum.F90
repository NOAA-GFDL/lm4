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
module sphum_mod

use constants_mod,      only: rdgas, rvgas
use sat_vapor_pres_mod, only: escomp
use fms_mod,            only: WARNING
use land_debug_mod,     only: check_temp_range

implicit none
private

public :: qscomp

! ==== module constants ======================================================
real, parameter :: d622 = rdgas/rvgas
real, parameter :: d378 = 1.0-d622
real, parameter :: del_temp = 0.1 ! temperature increment for q_sat derivative calc.

contains ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

subroutine qscomp(T, p, qsat, DqsatDT )
  real, intent(in) :: T    ! temperature
  real, intent(in) :: p    ! pressure
  real, intent(out):: qsat ! saturated specific humidity
  real, intent(out), optional :: DqsatDT ! deriv of specific humidity w.r.t. T

  real :: esat ! sat. water vapor pressure

  call check_temp_range(T,'qscomp','temperature')

  ! calculate saturated specific humidity
  call escomp(T,esat)
  qsat = d622*esat /(p-d378*esat )

  ! if requested, calculate the derivative of qsat w.r.t. temperature
  if (present(DqsatDT)) then
     call escomp(T+del_temp,esat)
     DqsatDT = (d622*esat/(p-d378*esat)-qsat)/del_temp
  endif
end subroutine qscomp

end module sphum_mod
