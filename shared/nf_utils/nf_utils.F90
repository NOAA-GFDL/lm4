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
module nf_utils_mod

use nfu_mod ! netcdf utilities
use nfc_mod ! netcdf utilities for compressed files

implicit none
private

! ==== public interfaces =====================================================
! export stuff from nfu_mod
public :: nfu_inq_dim, nfu_inq_var, nfu_inq_att
public :: nfu_def_dim, nfu_def_var
public :: nfu_put_att
public :: nfu_get_dim, nfu_get_dim_bounds
public :: nfu_put_var, nfu_put_rec
public :: nfu_get_var, nfu_get_rec
public :: nfu_get_valid_range, nfu_is_valid, nfu_validtype, nfu_validtype2ascii
! export stuff from nfc_mod
public :: nfu_inq_compressed_dim, nfu_inq_compressed_var
public :: nfu_get_compressed_var
public :: nfu_put_compressed_var
public :: nfu_get_compressed_rec
! ==== end of public interfaces ==============================================

end module nf_utils_mod
