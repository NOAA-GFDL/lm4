module river_type_mod
!-----------------------------------------------------------------------
!                   GNU General Public License                        
!                                                                      
! This program is free software; you can redistribute it and/or modify it and  
! are expected to follow the terms of the GNU General Public License  
! as published by the Free Software Foundation; either version 2 of   
! the License, or (at your option) any later version.                 
!                                                                      
! MOM is distributed in the hope that it will be useful, but WITHOUT    
! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY  
! or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public    
! License for more details.                                           
!                                                                      
! For the full text of the GNU General Public License,                
! write to: Free Software Foundation, Inc.,                           
!           675 Mass Ave, Cambridge, MA 02139, USA.                   
! or see:   http://www.gnu.org/licenses/gpl.html                      
!-----------------------------------------------------------------------
  ! <CONTACT EMAIL="klf@gfdl.noaa.gov"> Kirsten Findell </CONTACT> 
  ! <CONTACT EMAIL="z1l@gfdl.noaa.gov"> Zhi Liang </CONTACT> 

  use time_manager_mod, only : time_type

  implicit none
  private

  !--- version information ---------------------------------------------
  character(len=128) :: version = '$Id: river_type.F90,v 15.0 2007/08/14 03:59:37 fms Exp $'
  character(len=128) :: tagname = '$Name: omsk $'

  !--- public interface ------------------------------------------------
  public :: river_type, Leo_Mad_trios

  !--- public data type ------------------------------------------------

  type river_type
     integer, dimension(:),     pointer :: is    => NULL(), ie    => NULL()
     integer, dimension(:),     pointer :: js    => NULL(), je    => NULL()
     real, dimension(:),        pointer :: lon   => NULL(), lat   => NULL()
     real, dimension(:),        pointer :: lonb  => NULL(), latb  => NULL()
     real, dimension(:,:),      pointer :: celllength    => NULL()
     real, dimension(:,:),      pointer :: cellarea      => NULL()
     integer, dimension(:,:),   pointer :: basinid       => NULL() 
     integer, dimension(:,:),   pointer :: fromcell      => NULL()
     integer, dimension(:,:),   pointer :: travel        => NULL()
     logical, dimension(:,:),   pointer :: pemask        => NULL()
     logical, dimension(:,:),   pointer :: gmask         => NULL()
     integer, dimension(:,:,:), pointer :: fromcell_coef => NULL()
     real, dimension(:,:),      pointer :: runoff        => NULL()
     real, dimension(:,:),      pointer :: runoff_s      => NULL()
     real, dimension(:,:),      pointer :: runoff_h      => NULL()
     real, dimension(:,:,:),    pointer :: runoff_c      => NULL()
     real, dimension(:,:),      pointer :: storage       => NULL()     
     real, dimension(:,:),      pointer :: storage_s     => NULL()     
     real, dimension(:,:),      pointer :: storage_h     => NULL()     
     real, dimension(:,:,:),    pointer :: storage_c     => NULL()     
     real, dimension(:,:),      pointer :: outflow       => NULL()
     real, dimension(:,:),      pointer :: outflow_s     => NULL()
     real, dimension(:,:),      pointer :: outflow_h     => NULL()
     real, dimension(:,:,:),    pointer :: outflow_c     => NULL()
     real, dimension(:,:),      pointer :: outflowmean   => NULL()
     real, dimension(:,:),      pointer :: So            => NULL()
     real, dimension(:,:),      pointer :: width         => NULL()
     real, dimension(:,:),      pointer :: depth         => NULL()
     real, dimension(:,:),      pointer :: vel           => NULL()
     real, dimension(:,:),      pointer :: snow_frac     => NULL()
     real, dimension(:,:),      pointer :: wtr_temp      => NULL()
     real, dimension(:,:,:),    pointer :: conc          => NULL()
     real, dimension(:,:),      pointer :: air_temp      => NULL()
     type (time_type)                   :: Time
     integer                            :: dt_fast, dt_slow
     integer                            :: nlon, nlat
  end type river_type

type Leo_Mad_trios
   real :: on_V, on_d, on_w
end type Leo_Mad_trios


end module river_type_mod
