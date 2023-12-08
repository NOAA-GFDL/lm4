module soil_BGC_GIMICS_mod

use fms_mod, only: error_mesg, NOTE

use land_constants_mod, only : N_C_TYPES, N_LITTER_POOLS, &
     c_shortname, c_longname, l_shortname, l_longname
use land_data_mod, only : log_version, lnd
use land_tile_mod, only : land_tile_type, land_tile_map, land_tile_enum_type, &
     first_elmt, loop_over_tiles
use land_tile_io_mod, only : land_restart_type, &
     init_land_restart, open_land_restart, save_land_restart, free_land_restart, &
     add_tile_data, get_tile_data, add_restart_axis, field_exists
use soil_BGC_GIMICS_type_mod, only : soil_BGC_GIMICS_t, init_GIMICS_state, &
     soil_BGC_diag_init_GIMICS, save_equilibration_data
use soil_tile_mod, only: num_l, zfull

implicit none; private

public :: soil_BGC_save_restart_GIMICS
public :: soil_BGC_init_GIMICS

! ---- module constants
character(len=*), parameter :: module_name = 'soil_BGC_GIMICS_mod'
#include "../../shared/version_variable.inc"

character(len=*), parameter :: filename_base = 'soil_BGC_GIMICS'

! to simplify the restart i/o, in restart i/o we treat two different parts of the soil
! (rhizosphere and bulk) as a two-element array,  wth indices and names defined below
integer, parameter :: N_S_PARTS    = 2 ! number of separate soil partitions (i.e. rhizosphere, bulk)
integer, parameter :: S_PART_RHIZ = 1, S_PART_BULK = 2 ! indices of soi partitions
character(16) :: s_part_name(N_S_PARTS) = [ 'rhiz', 'bulk' ] ! corresponding names of soil partitions
contains ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

! ============================================================================
subroutine soil_BGC_init_GIMICS( id_ug, id_zfull )
  integer,intent(in)  :: id_ug    !< Unstructured axis id
  integer,intent(in)  :: id_zfull !< Vertical (depth) axis id

  character(267)          :: filename ! restart file name
  type(land_restart_type) :: restart  ! restart file i/o object
  logical                 :: restart_exists
  integer :: i,k,l
  type(land_tile_enum_type)     :: ce    ! current tile list element
  type(land_tile_type), pointer :: tile  ! pointer to current tile

  call log_version(version, module_name, __FILE__)

  ! initialize diagnostics
  call soil_BGC_diag_init_GIMICS( id_ug, id_zfull )

  filename = 'INPUT/'//trim(filename_base)//'.nc'
  call open_land_restart(restart,filename,restart_exists)
  if (restart_exists) then
     call error_mesg('soil_BGC_init_GIMICS', 'reading NetCDF restart "'//trim(filename)//'"', NOTE)

     ! surface litter data
     do k = 1,N_LITTER_POOLS
        call get_tile_data(restart,trim('sfc_'//l_shortname(k))//'_dz',litt_dz_ptr,k)
        call get_tile_data(restart,trim('sfc_'//l_shortname(k))//'_metabolicLitterC',litt_metabolicLitterC_ptr,k)
        call get_tile_data(restart,trim('sfc_'//l_shortname(k))//'_structuralLitterC',litt_structuralLitterC_ptr,k)
        call get_tile_data(restart,trim('sfc_'//l_shortname(k))//'_chemResistantC',litt_chemResistantC_ptr,k)
        call get_tile_data(restart,trim('sfc_'//l_shortname(k))//'_availableC',litt_availableC_ptr,k)
        call get_tile_data(restart,trim('sfc_'//l_shortname(k))//'_microbesR',litt_microbesR_ptr,k)
        call get_tile_data(restart,trim('sfc_'//l_shortname(k))//'_microbesK',litt_microbesK_ptr,k)
     enddo

     do k = 1, N_C_TYPES
        if(field_exists(restart, 'negative_litter_C_'//trim(c_shortname(k)))) then
           call get_tile_data(restart,'negative_litter_C_'//trim(c_shortname(k)),litt_negative_C_ptr,k)
        endif
     enddo

     ! soil data
     call get_tile_data(restart,'fRhiz','zfull',soil_fRhiz_ptr)
     do k = 1, N_S_PARTS
        call get_tile_data(restart, trim(s_part_name(k))//'_metabolicLitterC','zfull', soil_metabolicLitterC_ptr, k)
        call get_tile_data(restart, trim(s_part_name(k))//'_structuralLitterC','zfull', soil_structuralLitterC_ptr, k)
        call get_tile_data(restart, trim(s_part_name(k))//'_protectedC','zfull', soil_protectedC_ptr, k)
        call get_tile_data(restart, trim(s_part_name(k))//'_chemResistantC','zfull', soil_chemResistantC_ptr, k)
        call get_tile_data(restart, trim(s_part_name(k))//'_availableC','zfull', soil_availableC_ptr, k)
        call get_tile_data(restart, trim(s_part_name(k))//'_microbesR','zfull', soil_microbesR_ptr, k)
        call get_tile_data(restart, trim(s_part_name(k))//'_microbesK','zfull', soil_microbesK_ptr, k)
     enddo

     call free_land_restart(restart)
  else
     call error_mesg('soil_BGC_init_GIMICS', 'cold-starting soil_BGC_GIMICS', NOTE)
     ! Go through all tiles and initialize GIMICS soil carbon state.
     ce = first_elmt(land_tile_map, ls=lnd%ls)
     do while(loop_over_tiles(ce,tile,l))
        if (.not.associated(tile%soilc)) cycle
        select type(s=>tile%soilc)
        class is (soil_BGC_GIMICS_t)
           call init_GIMICS_state(s)
        end select
     enddo
  endif

  filename = 'INPUT/'//trim(filename_base)//'_eq.nc'
  call open_land_restart(restart,filename,restart_exists)
  if (restart_exists) then
     call error_mesg('soil_BGC_init_GIMICS', 'reading NetCDF restart "'//trim(filename)//'"', NOTE)
!      call get_tile_data(restart,'asoil_in','zfull',soil_asoil_in_ptr)
!      call get_tile_data(restart,'fsc_in','zfull',soil_fsc_in_ptr)
!      call get_tile_data(restart,'ssc_in','zfull',soil_ssc_in_ptr)
     call free_land_restart(restart)
  endif

end subroutine

! ============================================================================
subroutine soil_BGC_save_restart_GIMICS(tile_dim_length, timestamp)
  integer, intent(in) :: tile_dim_length ! length of tile dim. in the output file
  character(*), intent(in) :: timestamp ! timestamp to add to the file name

  character(267) :: filename
  type(land_restart_type) :: restart ! restart file i/o object
  integer :: k

  filename = 'RESTART/'//trim(timestamp)//trim(filename_base)//'.nc'
  call init_land_restart(restart, filename, soilc_tile_exists, tile_dim_length)
  call add_restart_axis(restart,'zfull',zfull(1:num_l),.false.,"Z",'m','full level',sense=-1)

  ! surface litter data
  do k = 1,N_LITTER_POOLS
     call add_tile_data(restart,trim('sfc_'//l_shortname(k))//'_dz',&
        litt_dz_ptr,k,'Thickness of '//trim(l_longname(k))//' surface litter','m')

     call add_tile_data(restart,trim('sfc_'//l_shortname(k))//'_metabolicLitterC',&
        litt_metabolicLitterC_ptr,k,'Metabolic carbon density in '//trim(l_longname(k))//' surface litter','kg/m3')
     call add_tile_data(restart,trim('sfc_'//l_shortname(k))//'_structuralLitterC',&
        litt_structuralLitterC_ptr,k,'Structural carbon density in '//trim(l_longname(k))//' surface litter','kg/m3')
     call add_tile_data(restart,trim('sfc_'//l_shortname(k))//'_chemResistantC',&
        litt_chemResistantC_ptr,k,'Chemically resistant carbon density in '//trim(l_longname(k))//' surface litter','kg/m3')
     call add_tile_data(restart,trim('sfc_'//l_shortname(k))//'_availableC',&
        litt_availableC_ptr,k,'Available carbon density in '//trim(l_longname(k))//' surface litter','kg/m3')
     call add_tile_data(restart,trim('sfc_'//l_shortname(k))//'_microbesR',&
        litt_microbesR_ptr,k,'R microbes carbon density in '//trim(l_longname(k))//' surface litter','kg/m3')
     call add_tile_data(restart,trim('sfc_'//l_shortname(k))//'_microbesK',&
        litt_microbesK_ptr,k,'K microbes carbon density in '//trim(l_longname(k))//' surface litter','kg/m3')
  enddo

  ! soil data
  call add_tile_data(restart,'fRhiz', 'zfull', soil_fRhiz_ptr ,'fraction of rhizosphere in soil', 'm3/m3')
  do k = 1, N_S_PARTS
     call add_tile_data(restart, trim(s_part_name(k))//'_metabolicLitterC', 'zfull', &
         soil_metabolicLitterC_ptr, k, 'Metabolic carbon density in '//trim(s_part_name(k)), 'kg/m3')
     call add_tile_data(restart, trim(s_part_name(k))//'_structuralLitterC', 'zfull', &
         soil_structuralLitterC_ptr, k, 'Structural carbon density in '//trim(s_part_name(k)), 'kg/m3')
     call add_tile_data(restart, trim(s_part_name(k))//'_protectedC', 'zfull', &
         soil_protectedC_ptr, k, 'Protected carbon density in '//trim(s_part_name(k)), 'kg/m3')
     call add_tile_data(restart, trim(s_part_name(k))//'_chemResistantC', 'zfull', &
         soil_chemResistantC_ptr, k, 'Chemically resistant carbon density in '//trim(s_part_name(k)), 'kg/m3')
     call add_tile_data(restart, trim(s_part_name(k))//'_availableC', 'zfull', &
         soil_availableC_ptr, k, 'Available carbon density in '//trim(s_part_name(k)), 'kg/m3')
     call add_tile_data(restart, trim(s_part_name(k))//'_microbesR', 'zfull', &
         soil_microbesR_ptr, k, 'R microbes carbon density in '//trim(s_part_name(k)), 'kg/m3')
     call add_tile_data(restart, trim(s_part_name(k))//'_microbesK', 'zfull', &
         soil_microbesK_ptr, k, 'K microbes carbon density in '//trim(s_part_name(k)), 'kg/m3')
  enddo

  do k = 1, N_C_TYPES
     call add_tile_data(restart,'negative_litter_C_'//trim(c_shortname(k)),litt_negative_C_ptr,k,'accumulated negative '//trim(c_longname(k))//' C litter input','kg/m2')
  enddo

  call save_land_restart(restart)
  call free_land_restart(restart)

  if (save_equilibration_data) then
     filename = 'RESTART/'//trim(timestamp)//trim(filename_base)//'_eq.nc'
     call init_land_restart(restart, filename, soilc_tile_exists, tile_dim_length)
     call add_restart_axis(restart,'zfull',zfull(1:num_l),.false.,"Z",'m','full level',sense=-1)

!      call add_tile_data(restart,'asoil_in','zfull',soil_asoil_in_ptr,'aerobic activity modifier', 'unitless')
!      call add_tile_data(restart,'fsc_in','zfull',soil_fsc_in_ptr,'fast soil carbon input', 'kg C/m2')
!      call add_tile_data(restart,'ssc_in','zfull',soil_ssc_in_ptr,'slow soil carbon input', 'kg C/m2')
!
!      call save_land_restart(restart)
!      call free_land_restart(restart)
  endif
end subroutine

! ============================================================================
! tile existence detector: returns a logical value indicating whether component
! model tile exists or not
logical function soilc_tile_exists(tile)
  type(land_tile_type), pointer :: tile
  soilc_tile_exists = associated(tile%soilc)
end function

! ============================================================================
! accessor functions for GIMICS model soil BGC data

! surface litter properties

subroutine litt_dz_ptr(t,i,p)
  type(land_tile_type), pointer :: t
  integer, intent(in) :: i
  real, pointer :: p
  p=>NULL()
  if(.not.associated(t))       return
  if(.not.associated(t%soilc)) return
  select type(s=>t%soilc)
  class is (soil_BGC_GIMICS_t)
     p=>s%litt(i)%dz
  end select
end subroutine

subroutine litt_metabolicLitterC_ptr(t,i,p)
  type(land_tile_type), pointer :: t
  integer, intent(in) :: i
  real, pointer :: p
  p=>NULL()
  if(.not.associated(t))       return
  if(.not.associated(t%soilc)) return
  select type(s=>t%soilc)
  class is (soil_BGC_GIMICS_t)
     p=>s%litt(i)%metabolicLitterC
  end select
end subroutine

subroutine litt_structuralLitterC_ptr(t,i,p)
  type(land_tile_type), pointer :: t
  integer, intent(in) :: i
  real, pointer :: p
  p=>NULL()
  if(.not.associated(t))       return
  if(.not.associated(t%soilc)) return
  select type(s=>t%soilc)
  class is (soil_BGC_GIMICS_t)
     p=>s%litt(i)%structuralLitterC
  end select
end subroutine

subroutine litt_chemResistantC_ptr(t,i,p)
  type(land_tile_type), pointer :: t
  integer, intent(in) :: i
  real, pointer :: p
  p=>NULL()
  if(.not.associated(t))       return
  if(.not.associated(t%soilc)) return
  select type(s=>t%soilc)
  class is (soil_BGC_GIMICS_t)
     p=>s%litt(i)%chemResistantC
  end select
end subroutine

subroutine litt_availableC_ptr(t,i,p)
  type(land_tile_type), pointer :: t
  integer, intent(in) :: i
  real, pointer :: p
  p=>NULL()
  if(.not.associated(t))       return
  if(.not.associated(t%soilc)) return
  select type(s=>t%soilc)
  class is (soil_BGC_GIMICS_t)
     p=>s%litt(i)%availableC
  end select
end subroutine

subroutine litt_microbesR_ptr(t,i,p)
  type(land_tile_type), pointer :: t
  integer, intent(in) :: i
  real, pointer :: p
  p=>NULL()
  if(.not.associated(t))       return
  if(.not.associated(t%soilc)) return
  select type(s=>t%soilc)
  class is (soil_BGC_GIMICS_t)
     p=>s%litt(i)%microbesR
  end select
end subroutine

subroutine litt_microbesK_ptr(t,i,p)
  type(land_tile_type), pointer :: t
  integer, intent(in) :: i
  real, pointer :: p
  p=>NULL()
  if(.not.associated(t))       return
  if(.not.associated(t%soilc)) return
  select type(s=>t%soilc)
  class is (soil_BGC_GIMICS_t)
     p=>s%litt(i)%microbesK
  end select
end subroutine

subroutine litt_negative_C_ptr(t,i,p)
  type(land_tile_type),pointer::t; integer,intent(in)::i; real,pointer::p; p=>NULL()
  if(.not.associated(t))       return
  if(.not.associated(t%soilc)) return
  select type(s=>t%soilc);
  class is (soil_BGC_GIMICS_t)
    p=>s%neg_litt_C(i)
  end select
end subroutine

! ---- soil properties

subroutine soil_fRhiz_ptr(t,i,p)
  type(land_tile_type), pointer :: t
  integer, intent(in) :: i
  real, pointer :: p
  p=>NULL()
  if(.not.associated(t))       return
  if(.not.associated(t%soilc)) return
  select type(s=>t%soilc)
  class is (soil_BGC_GIMICS_t)
     p=>s%fRhiz(i)
  end select
end subroutine

subroutine soil_metabolicLitterC_ptr(t,i,j,p)
  type(land_tile_type), pointer :: t
  integer, intent(in) :: i
  integer, intent(in) :: j ! rhizosphere or bulk
  real, pointer :: p
  p=>NULL()
  if(.not.associated(t))       return
  if(.not.associated(t%soilc)) return
  select type(s=>t%soilc)
  class is (soil_BGC_GIMICS_t)
    select case (j)
      case (S_PART_RHIZ); p=>s%rhiz(i)%metabolicLitterC
      case (S_PART_BULK); p=>s%bulk(i)%metabolicLitterC
    end select
  end select
end subroutine

subroutine soil_structuralLitterC_ptr(t,i,j,p)
  type(land_tile_type), pointer :: t
  integer, intent(in) :: i
  integer, intent(in) :: j ! rhizosphere or bulk
  real, pointer :: p
  p=>NULL()
  if(.not.associated(t))       return
  if(.not.associated(t%soilc)) return
  select type(s=>t%soilc)
  class is (soil_BGC_GIMICS_t)
    select case (j)
      case (S_PART_RHIZ); p=>s%rhiz(i)%structuralLitterC
      case (S_PART_BULK); p=>s%bulk(i)%structuralLitterC
    end select
  end select
end subroutine

subroutine soil_protectedC_ptr(t,i,j,p)
  type(land_tile_type), pointer :: t
  integer, intent(in) :: i
  integer, intent(in) :: j ! rhizosphere or bulk
  real, pointer :: p
  p=>NULL()
  if(.not.associated(t))       return
  if(.not.associated(t%soilc)) return
  select type(s=>t%soilc)
  class is (soil_BGC_GIMICS_t)
    select case (j)
      case (S_PART_RHIZ); p=>s%rhiz(i)%protectedC
      case (S_PART_BULK); p=>s%bulk(i)%protectedC
    end select
  end select
end subroutine

subroutine soil_chemResistantC_ptr(t,i,j,p)
  type(land_tile_type), pointer :: t
  integer, intent(in) :: i
  integer, intent(in) :: j ! rhizosphere or bulk
  real, pointer :: p
  p=>NULL()
  if(.not.associated(t))       return
  if(.not.associated(t%soilc)) return
  select type(s=>t%soilc)
  class is (soil_BGC_GIMICS_t)
    select case (j)
      case (S_PART_RHIZ); p=>s%rhiz(i)%chemResistantC
      case (S_PART_BULK); p=>s%bulk(i)%chemResistantC
    end select
  end select
end subroutine

subroutine soil_availableC_ptr(t,i,j,p)
  type(land_tile_type), pointer :: t
  integer, intent(in) :: i
  integer, intent(in) :: j ! rhizosphere or bulk
  real, pointer :: p
  p=>NULL()
  if(.not.associated(t))       return
  if(.not.associated(t%soilc)) return
  select type(s=>t%soilc)
  class is (soil_BGC_GIMICS_t)
    select case (j)
      case (S_PART_RHIZ); p=>s%rhiz(i)%availableC
      case (S_PART_BULK); p=>s%bulk(i)%availableC
    end select
  end select
end subroutine

subroutine soil_microbesR_ptr(t,i,j,p)
  type(land_tile_type), pointer :: t
  integer, intent(in) :: i
  integer, intent(in) :: j ! rhizosphere or bulk
  real, pointer :: p
  p=>NULL()
  if(.not.associated(t))       return
  if(.not.associated(t%soilc)) return
  select type(s=>t%soilc)
  class is (soil_BGC_GIMICS_t)
    select case (j)
      case (S_PART_RHIZ); p=>s%rhiz(i)%microbesR
      case (S_PART_BULK); p=>s%bulk(i)%microbesR
    end select
  end select
end subroutine

subroutine soil_microbesK_ptr(t,i,j,p)
  type(land_tile_type), pointer :: t
  integer, intent(in) :: i
  integer, intent(in) :: j ! rhizosphere or bulk
  real, pointer :: p
  p=>NULL()
  if(.not.associated(t))       return
  if(.not.associated(t%soilc)) return
  select type(s=>t%soilc)
  class is (soil_BGC_GIMICS_t)
    select case (j)
      case (S_PART_RHIZ); p=>s%rhiz(i)%microbesK
      case (S_PART_BULK); p=>s%bulk(i)%microbesK
    end select
  end select
end subroutine

end module
