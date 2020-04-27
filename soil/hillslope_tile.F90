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
! ============================================================================
! hillslope_tile_module
! ============================================================================
module hillslope_tile_mod

use land_tile_selectors_mod, only : &
     tile_selector_type, SEL_HLSP, register_tile_selector
use soil_tile_mod, only : soil_tile_type

implicit none
private

! ==== public interfaces =====================================================
public :: register_hlsp_selectors
public :: hlsp_is_selected
! =====end of public interfaces ==============================================

contains

! ============================================================================
! Define hillslope selector types to be used in diagnostic output.
subroutine register_hlsp_selectors(max_num_topo_hlsps, num_vertclusters, diagnostics_by_cluster)

   integer, intent(in)  :: max_num_topo_hlsps, num_vertclusters
   logical, intent(in)  :: diagnostics_by_cluster
   character(len=16)    :: selector_name
   integer              :: j  ! cluster level
   character(len=16)    :: jchar ! cluster level as char

   if (max_num_topo_hlsps > 1) then
      call register_tile_selector('domhslope', long_name='tiles from the dominant hillslope in each'//&
          ' gridcell', tag = SEL_HLSP, idata1 = 0 )
   end if

   if (num_vertclusters > 1) then
      call register_tile_selector('lowland', long_name='tiles closest to stream in each hillslope',&
          tag = SEL_HLSP, idata1 = 1 )
      call register_tile_selector('upland', long_name='tiles farthest from stream in each hillslope',&
          tag = SEL_HLSP, idata1 = num_vertclusters )
   end if

   if (diagnostics_by_cluster) then ! create selector for each vertical cluster level in hillslope
      do j=1,num_vertclusters
         write(jchar,*) j
         selector_name= 'hlspcluster' // trim(adjustl(jchar))
         call register_tile_selector(selector_name, long_name='tiles in hillslope cluster '//trim(adjustl(jchar)),&
                tag = SEL_HLSP, idata1 = j )

      end do
   end if

end subroutine register_hlsp_selectors


! =============================================================================
! returns true if tile fits the specified selector
function hlsp_is_selected(soil, sel)
  logical hlsp_is_selected
  type(tile_selector_type),  intent(in) :: sel
  type(soil_tile_type),      intent(in) :: soil

  hlsp_is_selected = (sel%idata1==0 .and. soil%hidx_k==1) .or. (sel%idata1==soil%hidx_j)
end function


end module hillslope_tile_mod
