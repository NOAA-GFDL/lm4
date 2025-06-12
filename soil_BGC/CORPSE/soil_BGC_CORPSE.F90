module soil_BGC_CORPSE_mod

use fms_mod, only: error_mesg, NOTE

use land_constants_mod, only : N_C_TYPES, N_LITTER_POOLS, &
     c_shortname, c_longname, l_shortname, l_longname
use land_data_mod, only : log_version
use land_tile_mod, only: land_tile_type, land_tile_enum_type, first_elmt, land_tile_map, loop_over_tiles
use land_tile_io_mod, only: land_restart_type, &
     init_land_restart, open_land_restart, save_land_restart, free_land_restart, &
     add_tile_data, add_int_tile_data, get_tile_data, get_int_tile_data, &
     add_restart_axis, field_exists
use soil_tile_mod, only: num_l, zfull
use soil_BGC_CORPSE_type_mod, only: soil_BGC_CORPSE_t, adjust_pool_ncohorts, do_nitrogen, &
     soilMaxCohorts, soil_BGC_diag_init_CORPSE, save_equilibration_data

implicit none; private

public :: soil_BGC_init_CORPSE
public :: soil_BGC_save_restart_CORPSE

! ==== module constants ======================================================
character(len=*), parameter :: module_name = 'soil_BGC_CORPSE_mod'
#include "../../shared/version_variable.inc"
character(len=*), parameter :: filename_base = 'soil_BGC_CORPSE'

contains ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

! ============================================================================
subroutine soil_BGC_init_CORPSE( id_ug, id_zfull )
  integer,intent(in) :: id_ug    !< Unstructured axis id
  integer,intent(in) :: id_zfull !< Vertical (depth) axis id

  character(267)          :: filename ! restart file name
  type(land_restart_type) :: restart ! restart file i/o object
  logical                 :: restart_exists
  type(land_tile_enum_type)     :: ce   ! tile list enumerator
  type(land_tile_type), pointer :: tile   ! pointer to current tile
  integer :: i,k

  call log_version(version, module_name, __FILE__)

  ! initialize diagnostics
  call soil_BGC_diag_init_CORPSE( id_ug, id_zfull )

  filename = 'INPUT/'//trim(filename_base)//'.nc'
  call open_land_restart(restart,filename,restart_exists)
  if (restart_exists) then
     call error_mesg('soil_BGC_init_CORPSE', 'reading NetCDF restart "'//trim(filename)//'"', NOTE)

     ce = first_elmt(land_tile_map)
     do while(loop_over_tiles(ce,tile))
         if (.not.associated(tile%soil)) cycle
         select type (sc=>tile%soilc)
         class is (soil_BGC_CORPSE_t)
            do i = 1,N_LITTER_POOLS
               call adjust_pool_ncohorts(sc%litter_corpse(i))
            enddo
            do i = 1,num_l
               call adjust_pool_ncohorts(sc%org_matter(i))
            enddo
         end select
     end do
     do i = 1, N_C_TYPES
        call get_tile_data(restart,trim(c_shortname(i))//'_soil_C', 'zfull','soilCCohort', sc_soil_C_ptr,i)
        call get_tile_data(restart,trim(c_shortname(i))//'ProtectedC', 'zfull','soilCCohort', sc_protected_C_ptr,i)
        call get_tile_data(restart,'soil_DOC_'//trim(c_shortname(i)), 'zfull', sc_DOC_ptr,i)

        do k = 1, N_LITTER_POOLS
           call get_tile_data(restart,trim(l_shortname(k))//'_litter_'//trim(c_shortname(i))//'_C','litterCCohort',sc_litter_litterC_ptr,i,k)
           call get_tile_data(restart,trim(l_shortname(k))//'_litter_'//trim(c_shortname(i))//'ProtectedC','litterCCohort',sc_litter_protectedC_ptr,i,k)
           call get_tile_data(restart,trim(l_shortname(k))//'_litter_DOC_'//trim(c_shortname(i)),sc_litter_dissolved_carbon_ptr,i,k)
        enddo
     enddo
     call get_tile_data(restart,'liveMic', 'zfull','soilCCohort',sc_livingMicrobeC_ptr)
     call get_tile_data(restart,'CO2', 'zfull','soilCCohort',sc_CO2_ptr)

     do i = 1,N_LITTER_POOLS
        call get_tile_data(restart, trim(l_shortname(i))//'_litter_liveMic_C', 'litterCCohort', sc_litter_livingMicrobeC_ptr, i)
        call get_tile_data(restart, trim(l_shortname(i))//'_litter_CO2',       'litterCCohort', sc_litter_CO2_ptr, i)
     enddo

     if(field_exists(restart, 'gross_nitrogen_flux_into_tile')) then
        call get_tile_data(restart,'gross_nitrogen_flux_into_tile', soil_gross_nitrogen_flux_into_tile_ptr)
        call get_tile_data(restart,'gross_nitrogen_flux_out_of_tile', soil_gross_nitrogen_flux_out_of_tile_ptr)
     endif

     if(field_exists(restart, 'is_peat')) then
        call get_int_tile_data(restart, 'is_peat','zfull', soil_is_peat_ptr)
     endif
     if (field_exists(restart,'fast_soil_N')) then
        do i = 1, N_C_TYPES
           call get_tile_data(restart,trim(c_shortname(i))//'_soil_N', 'zfull','soilCCohort', sc_soil_N_ptr,i)
           call get_tile_data(restart,trim(c_shortname(i))//'ProtectedN', 'zfull','soilCCohort', sc_protected_N_ptr,i)
           call get_tile_data(restart,'soil_DON_'//trim(c_shortname(i)), 'zfull', sc_DON_ptr,i)

           do k = 1, N_LITTER_POOLS
              call get_tile_data(restart,trim(l_shortname(k))//'_litter_'//trim(c_shortname(i))//'_N','litterCCohort',sc_litter_litterN_ptr,i,k)
              call get_tile_data(restart,trim(l_shortname(k))//'_litter_'//trim(c_shortname(i))//'ProtectedN','litterCCohort',sc_litter_protectedN_ptr,i,k)
              call get_tile_data(restart,trim(l_shortname(k))//'_litter_DON_'//trim(c_shortname(i)),sc_litter_dissolved_nitrogen_ptr,i,k)
           enddo
        enddo
        call get_tile_data(restart,'liveMicN', 'zfull','soilCCohort', sc_livingMicrobeN_ptr)
        call get_tile_data(restart,'soil_NO3', 'zfull', sc_nitrate_ptr)
        call get_tile_data(restart,'soil_NH4', 'zfull', sc_ammonium_ptr)
        call get_tile_data(restart,'soil_nitrif', 'zfull', sc_nitrif_ptr)
        call get_tile_data(restart,'soil_denitrif', 'zfull', sc_denitrif_ptr)

        ! Leaving out cohort-level immobilization and mineralization fields for now -- BNS

        do k = 1,N_LITTER_POOLS
           call get_tile_data(restart, trim(l_shortname(k))//'_litter_liveMic_N', 'litterCCohort', sc_litter_livingMicrobeN_ptr,k)
           call get_tile_data(restart, trim(l_shortname(k))//'_litter_NO3', sc_litter_nitrate_ptr,k)
           call get_tile_data(restart, trim(l_shortname(k))//'_litter_NH4', sc_litter_ammonium_ptr,k)
           call get_tile_data(restart, trim(l_shortname(k))//'_litter_nitrif', sc_litter_nitrif_ptr,k)
           call get_tile_data(restart, trim(l_shortname(k))//'_litter_denitrif', sc_litter_denitrif_ptr,k)
        enddo
     endif
     do i = 1, N_C_TYPES
        if(field_exists(restart, 'negative_litter_C_'//trim(c_shortname(i)))) then
           call get_tile_data(restart,'negative_litter_C_'//trim(c_shortname(i)),sc_negative_litter_C_ptr,i)
        endif
        if(field_exists(restart, 'negative_litter_N_'//trim(c_shortname(i)))) then
           call get_tile_data(restart,'negative_litter_N_'//trim(c_shortname(i)),sc_negative_litter_N_ptr,i)
        endif
     enddo
     call free_land_restart(restart)
  else
     call error_mesg('soil_init', 'cold-starting soilc_CORPSE', NOTE)
  endif

  filename = 'INPUT/'//trim(filename_base)//'_eq.nc'
  call open_land_restart(restart,filename,restart_exists)
  if (restart_exists) then
     call error_mesg('soil_BGC_init_CORPSE', 'reading NetCDF restart "'//trim(filename)//'"', NOTE)
     do i = 1,N_C_TYPES
        ! C inputs
        call get_tile_data(restart,trim(c_shortname(i))//'_soil_C_in','zfull',sc_C_in_ptr, i)
        call get_tile_data(restart,trim(c_shortname(i))//'_soil_C_turnover','zfull',sc_C_turnover_ptr, i)
        call get_tile_data(restart,trim(c_shortname(i))//'ProtectedC_in','zfull',sc_protected_C_in_ptr, i)
        call get_tile_data(restart,trim(c_shortname(i))//'ProtectedC_turnover','zfull',sc_protected_C_turnover_ptr, i)
        ! N inputs
        call get_tile_data(restart,trim(c_shortname(i))//'_soil_N_in','zfull',sc_N_in_ptr, i)
        call get_tile_data(restart,trim(c_shortname(i))//'_soil_N_turnover','zfull',sc_N_turnover_ptr, i)
        call get_tile_data(restart,trim(c_shortname(i))//'ProtectedN_in','zfull',sc_protected_N_in_ptr, i)
        call get_tile_data(restart,trim(c_shortname(i))//'ProtectedN_turnover','zfull',sc_protected_N_turnover_ptr, i)
        ! C and N litter inputs
        do k = 1,N_LITTER_POOLS
           call get_tile_data(restart,trim(l_shortname(k))//'_litter_'//trim(c_shortname(i))//'_C_in',sc_litter_C_in_ptr,i,k)
           call get_tile_data(restart,trim(l_shortname(k))//'_litter_'//trim(c_shortname(i))//'_N_in',sc_litter_N_in_ptr,i,k)
           call get_tile_data(restart,trim(l_shortname(k))//'_litter_'//trim(c_shortname(i))//'_C_turnover', sc_litter_C_turnover_ptr, i,k)
           call get_tile_data(restart,trim(l_shortname(k))//'_litter_'//trim(c_shortname(i))//'_N_turnover', sc_litter_N_turnover_ptr, i,k)
        enddo
     enddo
     call free_land_restart(restart)
  endif
end subroutine

! ============================================================================
subroutine soil_BGC_save_restart_CORPSE(tile_dim_length, timestamp)
  integer, intent(in) :: tile_dim_length ! length of tile dim. in the output file
  character(*), intent(in) :: timestamp ! timestamp to add to the file name

  character(267) :: filename
  type(land_restart_type) :: restart ! restart file i/o object
  type(land_tile_enum_type)     :: ce   ! tile list enumerator
  type(land_tile_type), pointer :: tile   ! pointer to current tile
  integer :: i,k

  filename = 'RESTART/'//trim(timestamp)//trim(filename_base)//'.nc'
  call init_land_restart(restart, filename, soilc_tile_exists, tile_dim_length)
  call add_restart_axis(restart,'zfull',zfull(1:num_l),.false.,"Z",'m','full level',sense=-1)
  call add_restart_axis(restart,'soilCCohort',(/(float(i),i=1,soilMaxCohorts)/), .false.)
  call add_restart_axis(restart,'litterCCohort',(/1.0/),.false.)

  ! make sure all arrays of carbon cohorts are of the same length
  ce = first_elmt(land_tile_map)
  do while (loop_over_tiles(ce,tile))
      if (.not.associated(tile%soilc)) cycle
      select type (sc=>tile%soilc)
      class is (soil_BGC_CORPSE_t)
         do i = 1,N_LITTER_POOLS
            call adjust_pool_ncohorts(sc%litter_corpse(i))
         enddo
         do i = 1,num_l
            call adjust_pool_ncohorts(sc%org_matter(i))
         enddo
      end select
  end do
  do i = 1, N_C_TYPES
     call add_tile_data(restart,trim(c_shortname(i))//'_soil_C','zfull','soilCCohort',sc_soil_C_ptr,i,trim(c_longname(i))//' soil carbon','kg/m2')
     call add_tile_data(restart,trim(c_shortname(i))//'ProtectedC','zfull','soilCCohort',sc_protected_C_ptr,i,'Protected '//trim(c_longname(i))//' carbon','kg/m2')
     call add_tile_data(restart,'soil_DOC_'//trim(c_shortname(i)),'zfull',sc_DOC_ptr,i,'Dissolved '//trim(c_longname(i))//' carbon','kg/m2')
     do k = 1,N_LITTER_POOLS
        call add_tile_data(restart,trim(l_shortname(k))//'_litter_'//trim(c_shortname(i))//'_C','litterCCohort',sc_litter_litterC_ptr,i,k,trim(l_longname(k))//' litter '//trim(c_longname(i))//' C','kg/m2')
        call add_tile_data(restart,trim(l_shortname(k))//'_litter_'//trim(c_shortname(i))//'ProtectedC','litterCCohort',sc_litter_protectedC_ptr,i,k,trim(l_longname(k))//' litter '//trim(c_longname(i))//' protected C','kg/m2')
        call add_tile_data(restart,trim(l_shortname(k))//'_litter_DOC_'//trim(c_shortname(i)),sc_litter_dissolved_carbon_ptr,i,k,'Dissolved '//trim(l_longname(k))//' litter '//trim(c_longname(i))//' carbon','kg/m2')
     enddo
  enddo
  call add_tile_data(restart,'liveMic' ,'zfull','soilCCohort',sc_livingMicrobeC_ptr,'Living microbial carbon','kg/m2')
  call add_tile_data(restart,'CO2', 'zfull','soilCCohort',sc_CO2_ptr,'Cohort CO2 generated','kg/m2')

  do k = 1,N_LITTER_POOLS
     call add_tile_data(restart,trim(l_shortname(k))//'_litter_liveMic_C','litterCCohort',sc_litter_livingMicrobeC_ptr,k,trim(l_longname(k))//' litter live microbe C','kg/m2')
     call add_tile_data(restart,trim(l_shortname(k))//'_litter_CO2','litterCCohort',sc_litter_CO2_ptr,k,trim(l_longname(k))//' litter CO2 generated','kg/m2')
  enddo

  do i = 1, N_C_TYPES
     call add_tile_data(restart,'negative_litter_C_'//trim(c_shortname(i)),sc_negative_litter_C_ptr,i,'accumulated negative '//trim(c_longname(i))//' C litter input','kg/m2')
     call add_tile_data(restart,'negative_litter_N_'//trim(c_shortname(i)),sc_negative_litter_N_ptr,i,'accumulated negative '//trim(c_longname(i))//' N litter input','kg/m2')
  enddo

  call add_int_tile_data(restart,'is_peat','zfull',soil_is_peat_ptr,'Is layer peat?','Boolean')

  if (do_nitrogen) then
     do i = 1, N_C_TYPES
        call add_tile_data(restart,trim(c_shortname(i))//'_soil_N', 'zfull','soilCCohort', sc_soil_N_ptr,i,trim(c_longname(i))//' soil nitrogen','kg/m2')
        call add_tile_data(restart,trim(c_shortname(i))//'ProtectedN', 'zfull','soilCCohort', sc_protected_N_ptr,i,'Protected '//trim(c_longname(i))//' soil nitrogen','kg/m2')
        call add_tile_data(restart,'soil_DON_'//trim(c_shortname(i)), 'zfull', sc_DON_ptr,i,'Dissolved '//trim(c_longname(i))//' nitrogen','kg/m2')

        do k = 1,N_LITTER_POOLS
           call add_tile_data(restart,trim(l_shortname(k))//'_litter_'//trim(c_shortname(i))//'_N','litterCCohort',sc_litter_litterN_ptr,i,k,trim(l_longname(k))//' litter '//trim(c_longname(i))//' N','kg/m2')
           call add_tile_data(restart,trim(l_shortname(k))//'_litter_'//trim(c_shortname(i))//'ProtectedN','litterCCohort',sc_litter_protectedN_ptr,i,k,trim(l_longname(k))//' litter '//trim(c_longname(i))//' protected N','kg/m2')
           call add_tile_data(restart,trim(l_shortname(k))//'_litter_DON_'//trim(c_shortname(i)),sc_litter_dissolved_nitrogen_ptr,i,k,'Dissolved '//trim(l_longname(k))//' litter '//trim(c_longname(i))//' nitrogen','kg/m2')
        enddo
     enddo
     call add_tile_data(restart,'liveMicN', 'zfull','soilCCohort', sc_livingMicrobeN_ptr,'Living microbial nitrogen','kg/m2')
     call add_tile_data(restart,'soil_NO3', 'zfull', sc_nitrate_ptr,'Soil nitrate content','kg/m2')
     call add_tile_data(restart,'soil_NH4', 'zfull', sc_ammonium_ptr,'Soil ammonium content','kg/m2')
     call add_tile_data(restart,'soil_nitrif', 'zfull', sc_nitrif_ptr,'Soil cumulative nitrification','kg/m2')
     call add_tile_data(restart,'soil_denitrif', 'zfull', sc_denitrif_ptr,'Soil cumulative denitrification','kg/m2')

     do k = 1,N_LITTER_POOLS
        call add_tile_data(restart,trim(l_shortname(k))//'_litter_liveMic_N', 'litterCCohort', sc_litter_livingMicrobeN_ptr, k, trim(l_longname(k))//' litter live microbe N','kg/m2')
        call add_tile_data(restart,trim(l_shortname(k))//'_litter_NO3', sc_litter_nitrate_ptr, k, trim(l_longname(k))//' litter nitrate content','kg/m2')
        call add_tile_data(restart,trim(l_shortname(k))//'_litter_NH4', sc_litter_ammonium_ptr, k, trim(l_longname(k))//' litter ammonium content','kg/m2')
        call add_tile_data(restart,trim(l_shortname(k))//'_litter_nitrif', sc_litter_nitrif_ptr, k, trim(l_longname(k))//' litter cumulative nitrification','kg/m2')
        call add_tile_data(restart,trim(l_shortname(k))//'_litter_denitrif', sc_litter_denitrif_ptr, k, trim(l_longname(k))//' litter cumulative denitrification','kg/m2')
     enddo

     call add_tile_data(restart,'gross_nitrogen_flux_into_tile',soil_gross_nitrogen_flux_into_tile_ptr,'Cumulative nitrogen flux into tile','kg/m2')
     call add_tile_data(restart,'gross_nitrogen_flux_out_of_tile',soil_gross_nitrogen_flux_out_of_tile_ptr,'Cumulative nitrogen flux out of tile','kg/m2')

  endif

  call save_land_restart(restart)
  call free_land_restart(restart)

  if (save_equilibration_data) then
     filename = 'RESTART/'//trim(timestamp)//trim(filename_base)//'_eq.nc'
     call init_land_restart(restart, filename, soilc_tile_exists, tile_dim_length)
     call add_restart_axis(restart,'zfull',zfull(1:num_l),.false.,"Z",'m','full level',sense=-1)

     do i = 1,N_C_TYPES
        ! C inputs
        call add_tile_data(restart,trim(c_shortname(i))//'_soil_C_in','zfull',&
                 sc_C_in_ptr, i, trim(c_longname(i))//' soil carbon input', 'kg C/m2')
        call add_tile_data(restart,trim(c_shortname(i))//'_soil_C_turnover','zfull', &
                 sc_C_turnover_ptr, i, trim(c_longname(i))//' soil carbon turnover','year-1')
        call add_tile_data(restart,trim(c_shortname(i))//'ProtectedC_in','zfull',&
                 sc_protected_C_in_ptr, i,'protected '//trim(c_longname(i))//' soil carbon input', 'kg C/m2')
        call add_tile_data(restart,trim(c_shortname(i))//'ProtectedC_turnover','zfull', &
                 sc_protected_C_turnover_ptr, i, trim(c_longname(i))//' protected soil carbon turnover', 'year-1')
        ! N_inputs
        call add_tile_data(restart,trim(c_shortname(i))//'_soil_N_in','zfull',&
                 sc_N_in_ptr, i, trim(c_longname(i))//' soil nitrogen input', 'kg C/m2')
        call add_tile_data(restart,trim(c_shortname(i))//'_soil_N_turnover','zfull', &
                 sc_N_turnover_ptr, i, trim(c_longname(i))//' soil carbon turnover','year-1')
        call add_tile_data(restart,trim(c_shortname(i))//'ProtectedN_in','zfull',&
                 sc_protected_N_in_ptr, i,'protected '//trim(c_longname(i))//' soil nitrogen input', 'kg C/m2')
        call add_tile_data(restart,trim(c_shortname(i))//'ProtectedN_turnover','zfull', &
                 sc_protected_N_turnover_ptr, i, trim(c_longname(i))//' protected soil carbon turnover', 'year-1')
        ! C and N litter inputs
        do k = 1,N_LITTER_POOLS
           call add_tile_data(restart,trim(l_shortname(k))//'_litter_'//trim(c_shortname(i))//'_C_in', &
                 sc_litter_C_in_ptr, i, k, trim(c_longname(i))//' '//trim(l_longname(k))//' litter carbon input','kg C/m2')
           call add_tile_data(restart,trim(l_shortname(k))//'_litter_'//trim(c_shortname(i))//'_C_turnover', &
                 sc_litter_C_turnover_ptr, i, k, trim(c_longname(i))//' '//trim(l_longname(k))//' litter carbon turnover', 'year-1')
           call add_tile_data(restart,trim(l_shortname(k))//'_litter_'//trim(c_shortname(i))//'_N_in', &
                 sc_litter_N_in_ptr, i, k, trim(c_longname(i))//' '//trim(l_longname(k))//' litter nitrogen input','kg C/m2')
           call add_tile_data(restart,trim(l_shortname(k))//'_litter_'//trim(c_shortname(i))//'_N_turnover', &
                 sc_litter_C_turnover_ptr, i, k, trim(c_longname(i))//' '//trim(l_longname(k))//' litter nitrogen turnover', 'year-1')
        enddo
     enddo
     call save_land_restart(restart)
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
! accessor functions for CORPSE soil BGC data
subroutine soil_gross_nitrogen_flux_into_tile_ptr(t,p)
    type(land_tile_type),pointer::t
    real,pointer::p
    p=>NULL()
    if(associated(t))then
        if(associated(t%soil))p=>t%soil%gross_nitrogen_flux_into_tile
    endif
end subroutine

subroutine soil_gross_nitrogen_flux_out_of_tile_ptr(t,p)
    type(land_tile_type),pointer::t
    real,pointer::p
    p=>NULL()
    if(associated(t))then
        if(associated(t%soil))p=>t%soil%gross_nitrogen_flux_out_of_tile
    endif
end subroutine

subroutine soil_is_peat_ptr(t,i,p)
  type(land_tile_type),pointer::t; integer,intent(in)::i; integer,pointer::p
  p=>NULL()
  if(.not.associated(t))       return
  if(.not.associated(t%soilc)) return
  select type(s=>t%soilc); class is (soil_BGC_CORPSE_t)
     p=>s%is_peat(i)
  end select
end subroutine

subroutine sc_soil_C_ptr(t,i,j,k,p)
  type(land_tile_type),pointer::t; integer,intent(in)::i,j,k;real,pointer::p
  p=>NULL()
  if(.not.associated(t))       return
  if(.not.associated(t%soilc)) return
  select type(s=>t%soilc); class is (soil_BGC_CORPSE_t)
     p=>s%org_matter(i)%litterCohorts(j)%litterC(k)
  end select
end subroutine

subroutine sc_soil_N_ptr(t,i,j,k,p)
  type(land_tile_type),pointer::t; integer,intent(in)::i,j,k;real,pointer::p
  p=>NULL()
  if(.not.associated(t))       return
  if(.not.associated(t%soilc)) return
  select type(s=>t%soilc); class is (soil_BGC_CORPSE_t)
     p=>s%org_matter(i)%litterCohorts(j)%litterN(k)
  end select
end subroutine

subroutine sc_nitrate_ptr(t,i,p)
  type(land_tile_type),pointer::t; integer,intent(in)::i;real,pointer::p
  p=>NULL()
  if(.not.associated(t))       return
  if(.not.associated(t%soilc)) return
  select type(s=>t%soilc); class is (soil_BGC_CORPSE_t)
     p=>s%org_matter(i)%nitrate
  end select
end subroutine

subroutine sc_ammonium_ptr(t,i,p)
  type(land_tile_type),pointer::t; integer,intent(in)::i;real,pointer::p
  p=>NULL()
  if(.not.associated(t))       return
  if(.not.associated(t%soilc)) return
  select type(s=>t%soilc); class is (soil_BGC_CORPSE_t)
     p=>s%org_matter(i)%ammonium
  end select
end subroutine

subroutine sc_livingMicrobeC_ptr(t,i,j,p)
  type(land_tile_type),pointer::t; real,pointer::p; integer,intent(in)::i,j
  p=>NULL()
  if(.not.associated(t))       return
  if(.not.associated(t%soilc)) return
  select type(s=>t%soilc); class is (soil_BGC_CORPSE_t)
     p=>s%org_matter(i)%litterCohorts(j)%livingMicrobeC
  end select
end subroutine

subroutine sc_livingMicrobeN_ptr(t,i,j,p)
  type(land_tile_type),pointer::t; real,pointer::p; integer,intent(in)::i,j
  p=>NULL()
  if(.not.associated(t))       return
  if(.not.associated(t%soilc)) return
  select type(s=>t%soilc); class is (soil_BGC_CORPSE_t)
     p=>s%org_matter(i)%litterCohorts(j)%livingMicrobeN
  end select
end subroutine

subroutine sc_CO2_ptr(t,i,j,p)
  type(land_tile_type),pointer::t; real,pointer::p; integer,intent(in)::i,j
  p=>NULL()
  if(.not.associated(t))       return
  if(.not.associated(t%soilc)) return
  select type(s=>t%soilc); class is (soil_BGC_CORPSE_t)
     p=>s%org_matter(i)%litterCohorts(j)%CO2
  end select
end subroutine

subroutine sc_litter_nitrate_ptr(t,i,p)
  type(land_tile_type),pointer::t; integer,intent(in)::i; real,pointer::p
  p=>NULL()
  if(.not.associated(t))       return
  if(.not.associated(t%soilc)) return
  select type(s=>t%soilc); class is (soil_BGC_CORPSE_t)
     p=>s%litter_corpse(i)%nitrate
  end select
end subroutine

subroutine sc_litter_ammonium_ptr(t,i,p)
  type(land_tile_type),pointer::t; integer,intent(in)::i; real,pointer::p
  p=>NULL()
  if(.not.associated(t))       return
  if(.not.associated(t%soilc)) return
  select type(s=>t%soilc); class is (soil_BGC_CORPSE_t)
     p=>s%litter_corpse(i)%ammonium
  end select
end subroutine

subroutine sc_litter_CO2_ptr(t,i,j,p)
  type(land_tile_type),pointer::t; integer,intent(in)::i,j; real,pointer::p
  p=>NULL()
  if(.not.associated(t))       return
  if(.not.associated(t%soilc)) return
  select type(s=>t%soilc); class is (soil_BGC_CORPSE_t)
     p=>s%litter_corpse(j)%litterCohorts(i)%CO2
  end select
end subroutine

subroutine sc_nitrif_ptr(t,i,p)
  type(land_tile_type),pointer::t; integer,intent(in)::i;real,pointer::p
  p=>NULL()
  if(.not.associated(t))       return
  if(.not.associated(t%soilc)) return
  select type(s=>t%soilc); class is (soil_BGC_CORPSE_t)
     p=>s%org_matter(i)%nitrif
  end select
end subroutine

subroutine sc_denitrif_ptr(t,i,p)
  type(land_tile_type),pointer::t; integer,intent(in)::i;real,pointer::p
  p=>NULL()
  if(.not.associated(t))       return
  if(.not.associated(t%soilc)) return
  select type(s=>t%soilc); class is (soil_BGC_CORPSE_t)
     p=>s%org_matter(i)%denitrif
  end select
end subroutine

subroutine sc_litter_nitrif_ptr(t,i,p)
  type(land_tile_type),pointer::t; integer,intent(in)::i; real,pointer::p
  p=>NULL()
  if(.not.associated(t))       return
  if(.not.associated(t%soilc)) return
  select type(s=>t%soilc); class is (soil_BGC_CORPSE_t)
     p=>s%litter_corpse(i)%nitrif
  end select
end subroutine

subroutine sc_litter_denitrif_ptr(t,i,p)
  type(land_tile_type),pointer::t; integer,intent(in)::i; real,pointer::p
  p=>NULL()
  if(.not.associated(t))       return
  if(.not.associated(t%soilc)) return
  select type(s=>t%soilc); class is (soil_BGC_CORPSE_t)
     p=>s%litter_corpse(i)%denitrif
  end select
end subroutine


subroutine sc_C_in_ptr(t,i,j,p)
  type(land_tile_type),pointer::t; integer,intent(in)::i,j;real,pointer::p
  p=>NULL()
  if(.not.associated(t))       return
  if(.not.associated(t%soilc)) return
  select type(s=>t%soilc); class is (soil_BGC_CORPSE_t)
     p=>s%org_matter(i)%C_in(j)
  end select
end subroutine

subroutine sc_N_in_ptr(t,i,j,p)
  type(land_tile_type),pointer::t; integer,intent(in)::i,j;real,pointer::p
  p=>NULL()
  if(.not.associated(t))       return
  if(.not.associated(t%soilc)) return
  select type(s=>t%soilc); class is (soil_BGC_CORPSE_t)
     p=>s%org_matter(i)%N_in(j)
  end select
end subroutine

subroutine sc_litter_C_in_ptr(t,i,j,p)
  type(land_tile_type),pointer::t; integer,intent(in)::i,j; real,pointer::p
  p=>NULL()
  if(.not.associated(t))       return
  if(.not.associated(t%soilc)) return
  select type(s=>t%soilc); class is (soil_BGC_CORPSE_t)
     p=>s%litter_corpse(j)%C_in(i)
  end select
end subroutine

subroutine sc_litter_N_in_ptr(t,i,j,p)
  type(land_tile_type),pointer::t; integer,intent(in)::i,j; real,pointer::p
  p=>NULL()
  if(.not.associated(t))       return
  if(.not.associated(t%soilc)) return
  select type(s=>t%soilc); class is (soil_BGC_CORPSE_t)
     p=>s%litter_corpse(j)%N_in(i)
  end select
end subroutine

subroutine sc_litter_dissolved_carbon_ptr(t,i,j,p)
  type(land_tile_type),pointer::t; integer,intent(in)::i,j; real,pointer::p
  p=>NULL()
  if(.not.associated(t))       return
  if(.not.associated(t%soilc)) return
  select type(s=>t%soilc); class is (soil_BGC_CORPSE_t)
     p=>s%litter_corpse(j)%dissolved_carbon(i)
  end select
end subroutine

subroutine sc_litter_dissolved_nitrogen_ptr(t,i,j,p)
  type(land_tile_type),pointer::t; integer,intent(in)::i,j; real,pointer::p
  p=>NULL()
  if(.not.associated(t))       return
  if(.not.associated(t%soilc)) return
  select type(s=>t%soilc); class is (soil_BGC_CORPSE_t)
     p=>s%litter_corpse(j)%dissolved_nitrogen(i)
  end select
end subroutine

subroutine sc_litter_livingMicrobeC_ptr(t,i,j,p)
  type(land_tile_type),pointer::t; integer,intent(in)::i,j; real,pointer::p
  p=>NULL()
  if(.not.associated(t))       return
  if(.not.associated(t%soilc)) return
  select type(s=>t%soilc); class is (soil_BGC_CORPSE_t)
     p=>s%litter_corpse(j)%litterCohorts(i)%livingMicrobeC
  end select
end subroutine

subroutine sc_litter_livingMicrobeN_ptr(t,i,j,p)
  type(land_tile_type),pointer::t; integer,intent(in)::i,j; real,pointer::p
  p=>NULL()
  if(.not.associated(t))       return
  if(.not.associated(t%soilc)) return
  select type(s=>t%soilc); class is (soil_BGC_CORPSE_t)
     p=>s%litter_corpse(j)%litterCohorts(i)%livingMicrobeN
  end select
end subroutine

subroutine sc_protected_C_in_ptr(t,i,j,p)
  type(land_tile_type),pointer::t; integer,intent(in)::i,j;real,pointer::p
  p=>NULL()
  if(.not.associated(t))       return
  if(.not.associated(t%soilc)) return
  select type(s=>t%soilc); class is (soil_BGC_CORPSE_t)
     p=>s%org_matter(i)%protected_C_in(j)
  end select
end subroutine

subroutine sc_protected_N_in_ptr(t,i,j,p)
  type(land_tile_type),pointer::t; integer,intent(in)::i,j;real,pointer::p
  p=>NULL()
  if(.not.associated(t))       return
  if(.not.associated(t%soilc)) return
  select type(s=>t%soilc); class is (soil_BGC_CORPSE_t)
     p=>s%org_matter(i)%protected_N_in(j)
  end select
end subroutine

subroutine sc_C_turnover_ptr(t,i,j,p)
  type(land_tile_type),pointer::t; integer,intent(in)::i,j;real,pointer::p
  p=>NULL()
  if(.not.associated(t))       return
  if(.not.associated(t%soilc)) return
  select type(s=>t%soilc); class is (soil_BGC_CORPSE_t)
     p=>s%org_matter(i)%C_turnover(j)
  end select
end subroutine

subroutine sc_N_turnover_ptr(t,i,j,p)
  type(land_tile_type),pointer::t; integer,intent(in)::i,j;real,pointer::p
  p=>NULL()
  if(.not.associated(t))       return
  if(.not.associated(t%soilc)) return
  select type(s=>t%soilc); class is (soil_BGC_CORPSE_t)
     p=>s%org_matter(i)%N_turnover(j)
  end select
end subroutine

subroutine sc_protected_C_turnover_ptr(t,i,j,p)
  type(land_tile_type),pointer::t; integer,intent(in)::i,j;real,pointer::p
  p=>NULL()
  if(.not.associated(t))       return
  if(.not.associated(t%soilc)) return
  select type(s=>t%soilc); class is (soil_BGC_CORPSE_t)
     p=>s%org_matter(i)%protected_C_turnover(j)
  end select
end subroutine

subroutine sc_protected_N_turnover_ptr(t,i,j,p)
  type(land_tile_type),pointer::t; integer,intent(in)::i,j;real,pointer::p
  p=>NULL()
  if(.not.associated(t))       return
  if(.not.associated(t%soilc)) return
  select type(s=>t%soilc); class is (soil_BGC_CORPSE_t)
     p=>s%org_matter(i)%protected_N_turnover(j)
  end select
end subroutine

subroutine sc_litter_C_turnover_ptr(t,i,j,p)
  type(land_tile_type),pointer::t; integer,intent(in)::i,j; real,pointer::p
  p=>NULL()
  if(.not.associated(t))       return
  if(.not.associated(t%soilc)) return
  select type(s=>t%soilc); class is (soil_BGC_CORPSE_t)
     p=>s%litter_corpse(j)%C_turnover(i)
  end select
end subroutine

subroutine sc_litter_N_turnover_ptr(t,i,j,p)
  type(land_tile_type),pointer::t; integer,intent(in)::i,j; real,pointer::p
  p=>NULL()
  if(.not.associated(t))       return
  if(.not.associated(t%soilc)) return
  select type(s=>t%soilc); class is (soil_BGC_CORPSE_t)
     p=>s%litter_corpse(j)%N_turnover(i)
  end select
end subroutine

subroutine sc_protected_C_ptr(t,i,j,k,p)
  type(land_tile_type),pointer::t; integer,intent(in)::i,j,k;real,pointer::p
  p=>NULL()
  if(.not.associated(t))       return
  if(.not.associated(t%soilc)) return
  select type(s=>t%soilc); class is (soil_BGC_CORPSE_t)
     p=>s%org_matter(i)%litterCohorts(j)%protectedC(k)
  end select
end subroutine

subroutine sc_protected_N_ptr(t,i,j,k,p)
  type(land_tile_type),pointer::t; integer,intent(in)::i,j,k;real,pointer::p
  p=>NULL()
  if(.not.associated(t))       return
  if(.not.associated(t%soilc)) return
  select type(s=>t%soilc); class is (soil_BGC_CORPSE_t)
     p=>s%org_matter(i)%litterCohorts(j)%protectedN(k)
  end select
end subroutine

subroutine sc_litter_litterC_ptr(t,i,j,k,p)
  type(land_tile_type),pointer::t
  integer,intent(in)::i,j,k
  real,pointer::p
  p=>NULL()
  if(.not.associated(t))       return
  if(.not.associated(t%soilc)) return
  select type(s=>t%soilc); class is (soil_BGC_CORPSE_t)
     p=>s%litter_corpse(k)%litterCohorts(i)%litterC(j)
  end select
end subroutine

subroutine sc_litter_litterN_ptr(t,i,j,k,p)
  type(land_tile_type),pointer::t
  integer,intent(in)::i,j,k
  real,pointer::p
  p=>NULL()
  if(.not.associated(t))       return
  if(.not.associated(t%soilc)) return
  select type(s=>t%soilc); class is (soil_BGC_CORPSE_t)
     p=>s%litter_corpse(k)%litterCohorts(i)%litterN(j)
  end select
end subroutine

subroutine sc_litter_protectedC_ptr(t,i,j,k,p)
  type(land_tile_type),pointer::t
  integer,intent(in)::i,j,k
  real,pointer::p
  p=>NULL()
  if(.not.associated(t))       return
  if(.not.associated(t%soilc)) return
  select type(s=>t%soilc); class is (soil_BGC_CORPSE_t)
     p=>s%litter_corpse(k)%litterCohorts(i)%protectedC(j)
  end select
end subroutine

subroutine sc_litter_protectedN_ptr(t,i,j,k,p)
  type(land_tile_type),pointer::t
  integer,intent(in)::i,j,k
  real,pointer::p
  p=>NULL()
  if(.not.associated(t))       return
  if(.not.associated(t%soilc)) return
  select type(s=>t%soilc); class is (soil_BGC_CORPSE_t)
     p=>s%litter_corpse(k)%litterCohorts(i)%protectedN(j)
  end select
end subroutine

subroutine sc_DOC_ptr(t,i,j,p)
  type(land_tile_type),pointer::t; integer,intent(in)::i,j;real,pointer::p
  p=>NULL()
  if(.not.associated(t))       return
  if(.not.associated(t%soilc)) return
  select type(s=>t%soilc); class is (soil_BGC_CORPSE_t)
     p=>s%org_matter(i)%dissolved_carbon(j)
  end select
end subroutine

subroutine sc_DON_ptr(t,i,j,p)
  type(land_tile_type),pointer::t; integer,intent(in)::i,j;real,pointer::p
  p=>NULL()
  if(.not.associated(t))       return
  if(.not.associated(t%soilc)) return
  select type(s=>t%soilc); class is (soil_BGC_CORPSE_t)
     p=>s%org_matter(i)%dissolved_nitrogen(j)
  end select
end subroutine

subroutine sc_negative_litter_C_ptr(t,i,p)
  type(land_tile_type),pointer::t; integer,intent(in)::i; real,pointer::p; p=>NULL()
  if(.not.associated(t))       return
  if(.not.associated(t%soilc)) return
  select type(s=>t%soilc); class is (soil_BGC_CORPSE_t)
    p=>s%neg_litt_C(i)
  end select
end subroutine

subroutine sc_negative_litter_N_ptr(t,i,p)
  type(land_tile_type),pointer::t; integer,intent(in)::i; real,pointer::p; p=>NULL()
  if(.not.associated(t))       return
  if(.not.associated(t%soilc)) return
  select type(s=>t%soilc); class is (soil_BGC_CORPSE_t)
    p=>s%neg_litt_N(i)
  end select
end subroutine

end module
