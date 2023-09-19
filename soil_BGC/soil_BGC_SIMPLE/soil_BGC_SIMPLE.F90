module soil_BGC_SIMPLE_mod

use fms_mod, only: error_mesg, NOTE

use land_constants_mod, only : N_C_TYPES, N_LITTER_POOLS, &
     c_shortname, c_longname, l_shortname, l_longname
use land_tile_mod, only : land_tile_type
use land_tile_io_mod, only : land_restart_type, &
     init_land_restart, open_land_restart, save_land_restart, free_land_restart, &
     add_tile_data, get_tile_data, add_restart_axis
use soil_BGC_simple_type_mod, only : soil_BGC_SIMPLE_t, soil_BGC_diag_init_SIMPLE, &
     save_equilibration_data
use soil_tile_mod, only: num_l, zfull

implicit none; private

public :: soil_BGC_save_restart_SIMPLE
public :: soil_BGC_init_SIMPLE

contains ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

! ============================================================================
subroutine soil_BGC_save_restart_SIMPLE(tile_dim_length, timestamp)
  integer, intent(in) :: tile_dim_length ! length of tile dim. in the output file
  character(*), intent(in) :: timestamp ! timestamp to add to the file name

  character(267) :: filename
  type(land_restart_type) :: restart ! restart file i/o object
  integer :: i,k

  filename = 'RESTART/'//trim(timestamp)//'soilc_CENT.nc'
  call init_land_restart(restart, filename, soilc_tile_exists, tile_dim_length)
  call add_restart_axis(restart,'zfull',zfull(1:num_l),.false.,"Z",'m','full level',sense=-1)

  call add_tile_data(restart,'fsc', 'zfull', soil_fast_soil_C_ptr ,'fast soil carbon', 'kg C/m2')
  call add_tile_data(restart,'ssc', 'zfull', soil_slow_soil_C_ptr ,'slow soil carbon', 'kg C/m2')
  do i = 1, N_C_TYPES
     do k = 1,N_LITTER_POOLS
        call add_tile_data(restart,trim(l_shortname(k))//'_litt_'//trim(c_shortname(i))//'_C',&
           litter_C_ptr,i,k,trim(l_longname(k))//' litter '//trim(c_longname(i))//' C','kg/m2')
     enddo
  enddo
  call save_land_restart(restart)
  call free_land_restart(restart)

  if (save_equilibration_data) then
     filename = 'RESTART/'//trim(timestamp)//'soilc_CENT_eq.nc'
     call init_land_restart(restart, filename, soilc_tile_exists, tile_dim_length)
     call add_restart_axis(restart,'zfull',zfull(1:num_l),.false.,"Z",'m','full level',sense=-1)

     call add_tile_data(restart,'asoil_in','zfull',soil_asoil_in_ptr,'aerobic activity modifier', 'unitless')
     call add_tile_data(restart,'fsc_in','zfull',soil_fsc_in_ptr,'fast soil carbon input', 'kg C/m2')
     call add_tile_data(restart,'ssc_in','zfull',soil_ssc_in_ptr,'slow soil carbon input', 'kg C/m2')

     call save_land_restart(restart)
     call free_land_restart(restart)
  endif
end subroutine

! ============================================================================
subroutine soil_BGC_init_SIMPLE( id_ug, id_zfull )
  integer,intent(in)  :: id_ug    !< Unstructured axis id
  integer,intent(in)  :: id_zfull !< Vertical (depth) axis id

  character(*), parameter :: restart_file_name = 'INPUT/soilc_CENT.nc'
  type(land_restart_type) :: restart ! restart file i/o object
  logical                 :: restart_exists
  integer :: i,k

  ! initialize diagnostics
  call soil_BGC_diag_init_SIMPLE( id_ug, id_zfull )

  call open_land_restart(restart,restart_file_name,restart_exists)
  if (restart_exists) then
     call error_mesg('soil_BGC_init_SIMPLE', 'reading NetCDF restart "'//trim(restart_file_name)//'"', NOTE)
     call get_tile_data(restart,'fsc','zfull',soil_fast_soil_C_ptr)
     call get_tile_data(restart,'ssc','zfull',soil_slow_soil_C_ptr)
     ! name is deliberately different from similar CORPSE fields so that we can start
     ! with CORPSE restarts with zero litter
     do i = 1, N_C_TYPES
        do k = 1, N_LITTER_POOLS
           call get_tile_data(restart,trim(l_shortname(k))//'_litt_'//trim(c_shortname(i))//'_C',litter_C_ptr,i,k)
        enddo
     enddo
  else
     call error_mesg('soil_BGC_init_SIMPLE', 'cold-starting soilc_CENT', NOTE)
  endif

  call open_land_restart(restart,'INPUT/soilc_CENT_eq.res.nc',restart_exists)
  if (restart_exists) then
     call error_mesg('soil_BGC_init_SIMPLE', 'reading NetCDF restart "soilc_CENT_eq.nc"', NOTE)
     call get_tile_data(restart,'asoil_in','zfull',soil_asoil_in_ptr)
     call get_tile_data(restart,'fsc_in','zfull',soil_fsc_in_ptr)
     call get_tile_data(restart,'ssc_in','zfull',soil_ssc_in_ptr)
     call free_land_restart(restart)
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
! accessor functions for SIMPLE model soil BGC data
subroutine soil_fast_soil_C_ptr(t,i,p)
  type(land_tile_type), pointer :: t
  integer, intent(in) :: i
  real, pointer :: p
  p=>NULL()
  if(.not.associated(t))       return
  if(.not.associated(t%soilc)) return
  select type(s=>t%soilc)
  class is (soil_BGC_SIMPLE_t)
    p=>s%fast_soil_C(i)
  end select
end subroutine

subroutine soil_slow_soil_C_ptr(t,i,p)
  type(land_tile_type), pointer :: t
  integer, intent(in) :: i
  real, pointer :: p
  p=>NULL()
  if(.not.associated(t))       return
  if(.not.associated(t%soilc)) return
  select type(s=>t%soilc)
  class is (soil_BGC_SIMPLE_t)
    p=>s%slow_soil_C(i)
  end select
end subroutine

subroutine litter_C_ptr(t,i,k,p)
  type(land_tile_type),pointer::t
  integer,intent(in)::i,k
  real,pointer::p
  p=>NULL()
  if(.not.associated(t))       return
  if(.not.associated(t%soilc)) return
  select type(s=>t%soilc)
  class is (soil_BGC_SIMPLE_t)
      p=>s%litter_SIMPLE_C(i,k)
  end select
end subroutine

subroutine soil_asoil_in_ptr(t,i,p)
  type(land_tile_type),pointer::t
  integer,intent(in)::i
  real,pointer::p
  p=>NULL()
  if(.not.associated(t))       return
  if(.not.associated(t%soilc)) return
  select type(s=>t%soilc)
  class is (soil_BGC_SIMPLE_t)
      p=>s%asoil_in(i)
  end select
end subroutine

subroutine soil_fsc_in_ptr(t,i,p)
  type(land_tile_type),pointer::t
  integer,intent(in)::i
  real,pointer::p
  p=>NULL()
  if(.not.associated(t))       return
  if(.not.associated(t%soilc)) return
  select type(s=>t%soilc)
  class is (soil_BGC_SIMPLE_t)
      p=>s%fsc_in(i)
  end select
end subroutine

subroutine soil_ssc_in_ptr(t,i,p)
  type(land_tile_type),pointer::t
  integer,intent(in)::i
  real,pointer::p
  p=>NULL()
  if(.not.associated(t))       return
  if(.not.associated(t%soilc)) return
  select type(s=>t%soilc)
  class is (soil_BGC_SIMPLE_t)
      p=>s%ssc_in(i)
  end select
end subroutine

end module
