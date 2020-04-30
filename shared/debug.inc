!***********************************************************************
!*                   GNU Lesser General Public License
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* FMS is free software: you can redistribute it and/or modify it under
!* the terms of the GNU Lesser General Public License as published by
!* the Free Software Foundation, either version 3 of the License, or (at
!* your option) any later version.
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
!* FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
!* for more details.
!*
!* You should have received a copy of the GNU Lesser General Public
!* License along with FMS.  If not, see <http://www.gnu.org/licenses/>.
!***********************************************************************
! -*-f90-*-
! $Id$

use land_debug_mod, only : dpri

#define __DEBUG1__(x) call dpri(#x,x);write(*,*)
#define __DEBUG2__(x1,x2) call dpri(#x1,x1);call dpri(#x2,x2);write(*,*)
#define __DEBUG3__(x1,x2,x3) call dpri(#x1,x1);call dpri(#x2,x2);call dpri(#x3,x3);write(*,*)
#define __DEBUG4__(x1,x2,x3,x4) call dpri(#x1,x1);call dpri(#x2,x2);call dpri(#x3,x3);call dpri(#x4,x4);write(*,*)
#define __DEBUG5__(x1,x2,x3,x4,x5) call dpri(#x1,x1);call dpri(#x2,x2);call dpri(#x3,x3);call dpri(#x4,x4);call dpri(#x5,x5);write(*,*)

#define __DEBUG___(x) call dpri(#x,x)
