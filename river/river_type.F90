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
  character(len=128) :: version = ''
  character(len=128) :: tagname = ''

  !--- public interface ------------------------------------------------
  public :: river_type, Leo_Mad_trios

  !--- public data type ------------------------------------------------

  type river_type
     integer, dimension(:),     pointer :: is    => NULL(), ie    => NULL()
     integer, dimension(:),     pointer :: js    => NULL(), je    => NULL()
     real, dimension(:),        pointer :: lon   => NULL(), lat   => NULL()
     real, dimension(:),        pointer :: lonb  => NULL(), latb  => NULL()
     real, dimension(:,:),      pointer :: celllength    => NULL()
     real, dimension(:,:),      pointer :: landfrac      => NULL()
     real, dimension(:,:),      pointer :: cellarea      => NULL()
     integer, dimension(:,:),   pointer :: basinid       => NULL() 
     integer, dimension(:,:),   pointer :: tocell        => NULL()
     integer, dimension(:,:),   pointer :: travel        => NULL()
     logical, dimension(:,:),   pointer :: pemask        => NULL()
     logical, dimension(:,:),   pointer :: gmask         => NULL()
     real, dimension(:,:),      pointer :: storage       => NULL()     
     real, dimension(:,:,:),    pointer :: storage_c     => NULL()     
     real, dimension(:,:),      pointer :: inflow        => NULL()
     real, dimension(:,:,:),    pointer :: inflow_c      => NULL()
     real, dimension(:,:),      pointer :: infloc        => NULL()
     real, dimension(:,:,:),    pointer :: infloc_c      => NULL()
     real, dimension(:,:),      pointer :: outflow       => NULL()
     real, dimension(:,:,:),    pointer :: outflow_c     => NULL()
     real, dimension(:,:),      pointer :: disw2o        => NULL()
!     real, dimension(:,:),      pointer :: diss2o        => NULL()
     real, dimension(:,:),      pointer :: disw2l        => NULL()
!     real, dimension(:,:),      pointer :: diss2l        => NULL()
     real, dimension(:,:,:),    pointer :: disc2o        => NULL()
     real, dimension(:,:,:),    pointer :: disc2l        => NULL()
     real, dimension(:,:,:),    pointer :: removal_c     => NULL()
     real, dimension(:,:),      pointer :: outflowmean   => NULL()
     real, dimension(:,:),      pointer :: o_coef        => NULL()
     real, dimension(:,:),      pointer :: o_exp         => NULL()
     real, dimension(:,:),      pointer :: d_coef        => NULL()
     real, dimension(:,:),      pointer :: d_exp         => NULL()
     real, dimension(:,:,:),    pointer :: source_conc   => NULL()
     real, dimension(:,:,:),    pointer :: source_flux   => NULL()
     real, dimension(:,:),      pointer :: So            => NULL()
     real, dimension(:,:),      pointer :: depth         => NULL()
!     real, dimension(:,:),      pointer :: width         => NULL()
!     real, dimension(:,:),      pointer :: vel           => NULL()
     real, dimension(:),        pointer :: t_ref         => NULL()
     real, dimension(:),        pointer :: vf_ref        => NULL()
     real, dimension(:),        pointer :: q10           => NULL()
     real, dimension(:),        pointer :: kinv          => NULL()
     type (time_type)                   :: Time
     integer                            :: dt_fast, dt_slow
     integer                            :: nlon, nlat, num_species, num_c
     integer                            :: num_phys
     logical                            :: do_age
  end type river_type

type Leo_Mad_trios
   real :: on_V, on_d, on_w
end type Leo_Mad_trios


end module river_type_mod
