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
module land_constants_mod

use constants_mod, only : rdgas, rvgas, wtmair, dens_h2o, grav

implicit none
private

! ==== public interfaces =====================================================
integer, public, parameter :: &
     NBANDS   = 2, & ! number of spectral bands for short-wave radiation calculations
     BAND_VIS = 1, & ! visible radiation (wavelenght range?)
     BAND_NIR = 2    ! near infra-red radiation (wavelenght range?)

real, public, parameter :: d622 = rdgas/rvgas
real, public, parameter :: d378 = 1.0-d622
real, public, parameter :: d608 = d378/d622

real, public, parameter :: Rugas = 8.314472 ! universal gas constant, J K-1 mol-1
real, public, parameter :: kBoltz= 1.3807e-23 ! Boltzmann's constant, J K-1 Rugas/Avogadro number

real, public, parameter :: days_per_year = 365.0
real, public, parameter :: seconds_per_year = 86400.0*days_per_year
real, public, parameter :: mol_C = 12.0e-3 ! molar mass of carbon, kg
real, public, parameter :: mol_air = wtmair/1000.0 ! molar mass of air, kg
real, public, parameter :: mol_CO2 = 44.00995e-3 ! molar mass of CO2,kg
real, public, parameter :: mol_h2o = 18.0e-3 ! molar mass of water, kg

real, public, parameter :: MPa_per_m = dens_h2o*grav*1.0e-6 ! pressure of one meter of water, Mega Pascal

end module
