module land_transitions_mod
#include <fms_platform.h>

#include "../shared/debug.inc"

use netcdf, only: nf90_max_name
use constants_mod, only : PI

use mpp_mod, only: input_nml_file
use fms_mod, only : string, error_mesg, FATAL, WARNING, NOTE, &
     mpp_pe, lowercase, &
     check_nml_error, stdlog, mpp_root_pe, fms_error_handler
use fms2_io_mod, only: FmsNetcdfFile_t, file_exists
use time_manager_mod, only : time_type, set_date, get_date, set_time, &
     operator(+), operator(-), operator(>), operator(<), operator(<=), operator(/), &
     operator(//), operator(==), days_in_year, get_time
use horiz_interp_mod, only : horiz_interp_init
use time_interp_mod, only : time_interp
use diag_manager_mod, only : register_diag_field, send_data, diag_field_add_attribute

use vegn_data_mod, only : &
     N_LU_TYPES, LU_PAST, LU_CROP, LU_NTRL, LU_SCND, LU_RANGE, LU_URBN, &
     landuse_name, landuse_longname

use cana_tile_mod, only : cana_tile_heat
use vegn_tile_mod, only : vegn_tile_heat, vegn_tile_type, vegn_tile_bwood, crop_type
use soil_tile_mod, only : soil_tile_heat

use land_tile_mod, only : land_tile_map, &
     land_tile_type, land_tile_list_type, land_tile_enum_type, new_land_tile, &
     first_elmt, tail_elmt, loop_over_tiles, operator(==), current_tile, &
     land_tile_list_init, land_tile_list_end, nitems, elmt_at_index, &
     erase, remove, insert, merge_land_tile_into_list, delete_land_tile, &
     get_tile_water, land_tile_carbon, land_tile_heat
use land_tile_diag_mod, only : cmor_name

use land_data_mod, only : lnd, log_version, horiz_interp_ug
use vegn_harvesting_mod, only : vegn_cut_forest, clear_all_on_conversion_to_crop

use land_debug_mod, only : set_current_point, is_watch_cell, &
     get_current_point, check_var_range, log_date, string_from_time
use land_numerics_mod, only : rank_descending

use transition_io_mod, only : transition_io_init, infile_T, varset_T

use vegn_debug_crop_mod, only: debug_crop

implicit none
private

! ==== public interface =====================================================
public :: land_transitions_init
public :: land_transitions_end
public :: save_land_transitions_restart

public :: land_transitions
! ==== end of public interface ==============================================

! ==== module constants =====================================================
character(len=*), parameter :: module_name = 'land_transitions_mod'
character(len=*), parameter :: diag_mod_name = 'landuse'
#include "../shared/version_variable.inc"

! selectors for overshoot handling options, for efficiency
integer, parameter :: &
     OPT_IGNORE = 0, &
     OPT_STOP   = 1, &
     OPT_REPORT = 2
integer, parameter :: &
     DISTR_LM3 = 0, &
     DISTR_MIN = 1
! order of transitions (resulting land use types, hight to low priority) for the
! min-n-tiles transition distribution option. ALL land use types MUST be present in this
! array, otherwise some transitions may be missed -- except perhaps LU_NTRL, since we
! assume there are no transitions to LU_NTRL
integer, parameter :: tran_order(N_LU_TYPES) = (/LU_URBN, LU_CROP, LU_PAST, LU_RANGE, LU_SCND, LU_NTRL/)

! TODO: describe differences between data sets

! ==== data types ===========================================================
! a description of single transition
type :: tran_type
   integer :: donor    = 0  ! kind of donor tile
   integer :: acceptor = 0  ! kind of acceptor tile
   real    :: frac     = 0  ! area of transition
end type tran_type

! ==== module data ==========================================================
logical :: module_is_initialized = .FALSE.

integer :: nlon_in, nlat_in

type(infile_T), target :: ftran, fstate
type(varset_T) :: input_tran  (N_LU_TYPES,N_LU_TYPES) ! input transition rate fields
type(varset_T) :: input_state (N_LU_TYPES,N_LU_TYPES) ! input state field (for initial transition only)

integer :: diag_ids  (N_LU_TYPES,N_LU_TYPES)
real, allocatable :: norm_in  (:,:) ! normalizing factor to convert input data to
        ! units of [fractions of vegetated area per year]
type(time_type) :: time0 ! time of previous transition calculations

type(crop_type) :: crop
logical :: crop_exists

integer :: tran_distr_opt = -1 ! selector for transition distribution option, for efficiency
integer :: overshoot_opt = -1 ! selector for overshoot handling options, for efficiency
integer :: conservation_opt = -1 ! selector for non-conservation handling options, for efficiency

! translation of luh2 names and LM3 landuse types
character(5) :: luh2name(13)
integer      :: luh2type(13)
integer :: idata
data (luh2name(idata), luh2type(idata), idata = 1, 13) / &
   'primf', LU_NTRL, &
   'primn', LU_NTRL, &
   'secdf', LU_SCND, &
   'secdn', LU_SCND, &
   'pltns', LU_SCND, &
   'urban', LU_CROP, &
   'c3ann', LU_CROP, &
   'c4ann', LU_CROP, &
   'c3per', LU_CROP, &
   'c4per', LU_CROP, &
   'c3nfx', LU_CROP, &
   'pastr', LU_PAST, &
   'range', LU_RANGE /

! variables for LUMIP diagnostics
integer, parameter :: N_LUMIP_TYPES = 4, &
   LUMIP_PSL = 1, LUMIP_PST = 2, LUMIP_CRP = 3, LUMIP_URB = 4
character(4), parameter :: lumip_name(N_LUMIP_TYPES) = ['psl ','pst ','crop','urbn']
integer :: &
   id_frac_in (N_LUMIP_TYPES) = -1, &
   id_frac_out(N_LUMIP_TYPES) = -1
! translation table: model land use types -> LUMIP types: for each of the model
! LU types it lists the corresponding LUMIP type.
integer, parameter :: lu2lumip(N_LU_TYPES) = [LUMIP_PST, LUMIP_CRP, LUMIP_PSL, LUMIP_PSL, LUMIP_URB, LUMIP_PST]
logical :: close_state_file = .false.

! ---- namelist variables ---------------------------------------------------
logical, protected, public :: do_landuse_change = .FALSE. ! if true, then the landuse changes with time
character(len=1024) :: input_file  = '' ! input data set of transition dates
character(len=1024) :: state_file  = '' ! input data set of LU states (for initial transition only)
character(len=1024) :: static_file = '' ! static data file, for input land fraction
character(len=16)  :: data_type  = 'luh1' ! or 'luh2'
! distribute_transitions sets how the land use transitions are distributed among
! tiles within grid cells. 'lm3' is traditional (transitions applied to every
! tile in equal measure, except secondary-to-secondary); 'min-tiles' applies
! transitions to tiles in the order of priority, thereby minimizing the number
! of resulting tiles
logical :: rangeland_is_pasture = .FALSE. ! if true, rangeland is combined with pastures.
! This only applies to luh2 transitions, since there is no rangeland in luh1 anyway.
character(len=16)  :: distribute_transitions  = 'lm3' ! or 'min-n-tiles'
! sets how to handle transition overshoot: that is, the situation when transition
! is larger than available area of the given land use type.
character(len=16) :: overshoot_handling = 'report' ! or 'stop', or 'ignore'
real :: overshoot_tolerance = 1e-4 ! tolerance interval for overshoots
! specifies how to handle non-conservation
character(len=16) :: conservation_handling = 'stop' ! or 'report', or 'ignore'

namelist/landuse_nml/do_landuse_change, input_file, state_file, static_file, data_type, &
     rangeland_is_pasture, distribute_transitions, &
     overshoot_handling, overshoot_tolerance, &
     conservation_handling


contains ! ###################################################################

! ============================================================================
subroutine land_transitions_init(id_ug, id_cellarea)
  integer, intent(in) :: id_ug !<Unstructured axis id.
  integer, intent(in) :: id_cellarea !<id of cell area diagnostic fields

  ! ---- local vars
  integer        :: unit, ierr, io
  integer        :: year,month,day,hour,min,sec
  integer        :: k1,k2,k3, n1,n2
  character(12)  :: fieldname

  type(land_tile_type), pointer :: tile
  type(land_tile_enum_type) :: ce
  logical :: exists
  character(len=nf90_max_name) :: name
  type(FmsNetcdfFile_t) :: fileobj_static
  integer :: ndims

  if(module_is_initialized) return
  module_is_initialized = .TRUE.
  call log_version(version, module_name, &
  __FILE__)

  call horiz_interp_init()
  call transition_io_init()

  read (input_nml_file, nml=landuse_nml, iostat=io)
  ierr = check_nml_error(io, 'landuse_nml')
  if (mpp_pe() == mpp_root_pe()) then
     unit=stdlog()
     write(unit, nml=landuse_nml)
  endif

  ! read restart file, if any
  if (file_exists('INPUT/landuse.res')) then
     call error_mesg('land_transitions_init','reading restart "INPUT/landuse.res"',&
          NOTE)
     open(newunit=unit, file='INPUT/landuse.res', action="read")
     read(unit,*) year,month,day,hour,min,sec
     time0 = set_date(year,month,day,hour,min,sec)
     close(unit)
  else
     call error_mesg('land_transitions_init','cold-starting land transitions',&
          NOTE)
     time0 = set_date(0001,01,01);
  endif

  ! parse the transition distribution option
  select case(trim(lowercase(distribute_transitions)))
  case ('lm3')
     tran_distr_opt = DISTR_LM3
  case ('min-n-tiles')
     tran_distr_opt = DISTR_MIN
  case default
     call error_mesg('land_transitions_init','distribute_transitions value "'//&
          trim(distribute_transitions)//'" is incorrect, use "lm3" or "min-n-tiles"',&
          FATAL)
  end select

  ! parse the overshoot handling option
  if (trim(overshoot_handling)=='stop') then
     overshoot_opt = OPT_STOP
  else if (trim(overshoot_handling)=='ignore') then
     overshoot_opt = OPT_IGNORE
  else if (trim(overshoot_handling)=='report') then
     overshoot_opt = OPT_REPORT
  else
     call error_mesg('land_transitions_init','overshoot_handling value "'//&
          trim(overshoot_handling)//'" is illegal, use "stop", "report", or "ignore"',&
          FATAL)
  endif

  ! parse the non-conservation handling option
  if (trim(conservation_handling)=='stop') then
     conservation_opt = OPT_STOP
  else if (trim(conservation_handling)=='ignore') then
     conservation_opt = OPT_IGNORE
  else if (trim(conservation_handling)=='report') then
     conservation_opt = OPT_REPORT
  else
     call error_mesg('land_transitions_init','conservation_handling value "'//&
          trim(conservation_handling)//'" is illegal, use "stop", "report", or "ignore"',&
          FATAL)
  endif

  ! initialize diagnostics
  diag_ids(:,:) = 0

  do k1 = 1,size(diag_ids,1)
  do k2 = 1,size(diag_ids,2)
     ! skip unnamed tiles
     if(landuse_name(k1)=='')cycle
     if(landuse_name(k2)=='')cycle
     ! construct a name of input field and register the field
     fieldname = trim(landuse_name(k1))//'2'//trim(landuse_name(k2))
     diag_ids(k1,k2) = register_diag_field(diag_mod_name,fieldname,(/id_ug/), lnd%time, &
          'rate of transition from '//trim(landuse_longname(k1))//' to '//trim(landuse_longname(k2)),&
          units='1/year', missing_value=-1.0)
  enddo
  enddo
  ! register CMIP/LUMIP transition fields
  do k1 = 1,N_LUMIP_TYPES
     id_frac_in(k1) = register_diag_field(cmor_name, &
         'fracInLut_'//trim(lumip_name(k1)), (/id_ug/), lnd%time, &
         'Gross Fraction That Was Transferred into This Tile From Other Land Use Tiles', &
         units='%', standard_name='area_fraction', area = id_cellarea)
     call diag_field_add_attribute(id_frac_in(k1),'ocean_fillvalue',0.0)
     id_frac_out(k1) = register_diag_field(cmor_name, &
         'fracOutLut_'//trim(lumip_name(k1)), (/id_ug/), lnd%time, &
         'Gross Fraction of Land Use Tile That Was Transferred into Other Land Use Tiles', &
         units='%', standard_name='area_fraction', area = id_cellarea)
     call diag_field_add_attribute(id_frac_out(k1),'ocean_fillvalue',0.0)
  enddo

  ! change rangeland to pasture.
  if (rangeland_is_pasture) then
     ! change the type of transitions
     do k1 = 1,size(luh2type)
        if (luh2type(k1)==LU_RANGE) luh2type(k1)=LU_PAST
     enddo
     ! change land use type in existing tiles
     ce = first_elmt(land_tile_map, ls=lnd%ls )
     do while(loop_over_tiles(ce,tile))
        if (.not.associated(tile%vegn)) cycle
        if (tile%vegn%landuse == LU_RANGE) tile%vegn%landuse = LU_PAST
     enddo
  endif

  if (.not.do_landuse_change) return ! do nothing more if no land use requested

  ! stop if landuse.res looks inconsistent
  if (time0>lnd%time) then
     call error_mesg('land_transitions_init',&
          'current model time ('//trim(string_from_time(lnd%time))// &
          ') must be after the time of last land use transition application ('// &
          trim(string_from_time(time0))//')',&
          FATAL)
  endif

  ! check that we are starting from potential vegetation
  if (time0==set_date(0001,01,01)) then
     ce = first_elmt(land_tile_map, ls=lnd%ls )
     do while(loop_over_tiles(ce,tile))
        if (.not.associated(tile%vegn)) cycle
        if (tile%vegn%landuse /= LU_NTRL) then
            call error_mesg('land_transitions_init', &
                'starting land use transitions, but the land use tiles already exist in the initial conditions', &
                FATAL)
        endif
     enddo
  endif

  if (trim(input_file)=='') call error_mesg('land_transitions_init', &
       'do_landuse_change is requested, but landuse transition file is not specified', &
       FATAL)

  ! initialize data structure representing input file and horizontal interpolator
  call ftran%init(input_file,static_file,data_type)

  ! initialize arrays of input fields
  select case (trim(lowercase(data_type)))
  case('luh1')
     do k1 = 1,size(input_tran,1)
     do k2 = 1,size(input_tran,2)
        ! construct a name of input field and register the field
        fieldname = trim(landuse_name(k1))//'2'//trim(landuse_name(k2))
        if(trim(fieldname)=='2') cycle ! skip unspecified tiles
        input_tran(k1,k2)%name=fieldname
        call input_tran(k1,k2)%addvar(ftran,fieldname)
     enddo
     enddo

  case('luh2')
     ! LUH2 data set has more land use types and transitions than LM3,
     ! therefore several transitions need to be aggregated on input to get
     ! the transitions among LM3 land use types
     do n1 = 1,size(luh2type)
     do n2 = 1,size(luh2type)
        k1 = luh2type(n1)
        k2 = luh2type(n2)
        input_tran(k1,k2)%name=trim(landuse_name(k1))//'2'//trim(landuse_name(k2))
        if (k1==k2.and.k1/=LU_SCND) cycle ! skip transitions to the same LM3 LU type, except scnd2scnd
        call input_tran(k1,k2)%addvar(ftran,trim(luh2name(n1))//'_to_'//trim(luh2name(n2)))
     enddo
     enddo
     ! add transitions that are not part "state1_to_state2" variable set
     call input_tran(LU_NTRL,LU_SCND)%addvar(ftran,'primf_harv')
     call input_tran(LU_NTRL,LU_SCND)%addvar(ftran,'primn_harv')
     call input_tran(LU_SCND,LU_SCND)%addvar(ftran,'secmf_harv')
     call input_tran(LU_SCND,LU_SCND)%addvar(ftran,'secyf_harv')
     call input_tran(LU_SCND,LU_SCND)%addvar(ftran,'secnf_harv')
     call input_tran(LU_SCND,LU_SCND)%addvar(ftran,'pltns_harv')

     if (time0==set_date(0001,01,01)) then
        call error_mesg('land_transitions_init','setting up initial land use transitions', NOTE)
        ! initialize state input for initial transition from all-natural state.
        if (trim(state_file)=='') call error_mesg('land_transitions_init',&
            'starting land use transitions, but land use state file is not specified',FATAL)

        ! open state file
        call fstate%init(state_file,static_file,data_type)

        ! initialize state variable array
        do n2 = 1,size(luh2type)
           k2 = luh2type(n2)
           if (k2==LU_NTRL) cycle
           input_state(LU_NTRL,k2)%name='initial '//trim(landuse_name(LU_NTRL))//'2'//trim(landuse_name(k2))
           call input_state(LU_NTRL,k2)%addvar(fstate,luh2name(n2))
        enddo
     endif
  case default
     call error_mesg('land_transitions_init','unknown data_type "'&
                    //trim(data_type)//'", use "luh1" or "luh2"', FATAL)
  end select
  if (mpp_pe()==mpp_root_pe()) then
     write(*,*)'land_transitions_init: summary of land use transitions'
     do k1 = 1,size(input_tran,1)
     do k2 = 1,size(input_tran,2)
        if(input_tran(k1,k2)%name/='') &
             write(*,'(a)') input_tran(k1,k2)%descr()
     enddo
     enddo
  endif

end subroutine land_transitions_init

! ============================================================================
subroutine land_transitions_end()
  module_is_initialized=.FALSE.
  ! close files and deallocate associated memory
  call ftran%destroy()
  call fstate%destroy()
end subroutine land_transitions_end

! ============================================================================
subroutine save_land_transitions_restart(timestamp)
  character(*), intent(in) :: timestamp ! timestamp to add to the file name

  integer :: unit,year,month,day,hour,min,sec

  if (mpp_pe() == mpp_root_pe()) then
     open(newunit=unit, file='RESTART/'//trim(timestamp)//'landuse.res', action="write")
     call get_date(time0, year,month,day,hour,min,sec)
     write(unit,'(6i6,8x,a)') year,month,day,hour,min,sec, &
          'Time of previous landuse transition calculation'
     close(unit)
  endif

end subroutine save_land_transitions_restart


! =============================================================================
subroutine land_transitions (time)
  type(time_type), intent(in) :: time

  ! ---- local vars.
  integer :: i,k1,k2,i1,i2,l
  real    :: frac(lnd%ls:lnd%le)
  type(tran_type), pointer :: transitions(:,:)
  integer :: second, minute, hour, day0, day1, month0, month1, year0, year1
  real    :: w
  real    :: diag(lnd%ls:lnd%le)
  logical :: used
  type(land_tile_type), pointer :: tile
  type(land_tile_enum_type) :: ce

  if (.not.do_landuse_change) &
       return ! do nothing if landuse change not requested
  ! NB: in this case file/interp/data are not initialized, so it is
  ! not even possible to use the code below

  call get_date(time,             year0,month0,day0,hour,minute,second)
  call get_date(time-lnd%dt_slow, year1,month1,day1,hour,minute,second)
  if(year0 == year1) &
!!$  if(day0 == day1) &
       return ! do nothing during a year

  if (mpp_pe()==mpp_root_pe()) &
       call log_date('land_transitions: applying land use transitions on ', time)

  ! get transition rates for current time: read map of transitions, and accumulate
  ! as many time steps in array of transitions as necessary. Note that "transitions"
  ! array gets reallocated inside add_to_transitions as necessary, it has only as many
  ! layers as the max number of transitions occurring at a point at the time.
  transitions => NULL()
  do k1 = 1,N_LU_TYPES
  do k2 = 1,N_LU_TYPES
     ! get transition rate for this specific transition
     frac(:) = 0.0
     if (time0==set_date(0001,01,01).and.fstate%initialized) then
        ! read initial transition from state file
        call time_interp(time, fstate%time_in, w, i1,i2)
        call input_state(k1,k2)%get_data(i1,frac)
     else
        if (input_tran(k1,k2)%nvars>0) then
           call integral_transition(time0,time,input_tran(k1,k2),frac)
        endif
     endif
     call add_to_transitions(frac,time0,time,k1,k2,transitions)
  enddo
  enddo

  ! save the "in" and "out" diagnostics for the transitions
  do k1 = 1, N_LUMIP_TYPES
     if (id_frac_out(k1) > 0) then
        diag(:) = 0.0
        do k2 = 1, size(transitions,2)
        do i = lnd%ls,lnd%le
           if (transitions(i,k2)%donor>0) then
              if (lu2lumip(transitions(i,k2)%donor) == k1) &
                    diag(i) = diag(i) + transitions(i,k2)%frac
           endif
        enddo
        enddo
        used=send_data(id_frac_out(k1), diag*lnd%ug_landfrac*100.0, time)
     endif
  enddo
  do k1 = 1, N_LUMIP_TYPES
     if (id_frac_in(k1) > 0) then
        diag(:) = 0.0
        do k2 = 1, size(transitions,2)
        do i = lnd%ls,lnd%le
           if (transitions(i,k2)%acceptor>0) then
              if (lu2lumip(transitions(i,k2)%acceptor) == k1) &
                    diag(i) = diag(i) + transitions(i,k2)%frac
           endif
        enddo
        enddo
        used=send_data(id_frac_in(k1), diag*lnd%ug_landfrac*100.0, time)
     endif
  enddo

  ! perform the transitions
  do l = lnd%ls,lnd%le
     ! set current point for debugging
     call set_current_point(l,1)
     ce = first_elmt(land_tile_map(l))
     crop_exists = .FALSE.
     do while (loop_over_tiles(ce,tile))
        if (associated(tile%vegn)) then
           if (tile%vegn%landuse == LU_CROP) then
              crop = tile%vegn%Crop
              crop_exists = .TRUE.
!             call debug_crop(tile%vegn,'preexisting crop tile before transitions')
           endif
        endif
     enddo
     ! transition land area between different tile types
     call land_transitions_0d(land_tile_map(l), &
          transitions(l,:)%donor, &
          transitions(l,:)%acceptor,&
          transitions(l,:)%frac )
     ce = first_elmt(land_tile_map(l))
     do while (loop_over_tiles(ce,tile))
        if (associated(tile%vegn)) then
           if (tile%vegn%landuse == LU_CROP) then
              if(crop_exists) then
                 tile%vegn%Crop = crop
!                call debug_crop(tile%vegn,'preexisting crop tile after transitions')
              else
!                call debug_crop(tile%vegn,'new crop tile after transitions')
              endif
           endif
        endif
     enddo
  enddo

  ! deallocate array of transitions
  if (associated(transitions)) deallocate(transitions)

  ! store current time for future reference
  time0=time

end subroutine land_transitions


! =============================================================================
! performs tile transitions in a given grid cell
subroutine land_transitions_0d(d_list,d_kinds,a_kinds,area)
  type(land_tile_list_type), intent(inout) :: d_list ! list of tiles
  integer, intent(in) :: d_kinds(:) ! array of donor tile kinds
  integer, intent(in) :: a_kinds(:) ! array of acceptor tile kinds
  real   , intent(in) :: area(:)    ! array of areas changing from donor tiles to acceptor tiles

  ! ---- local vars
  integer :: i, k
  type(land_tile_type), pointer :: ptr
  type(land_tile_list_type) :: a_list
  type(land_tile_enum_type) :: ts, te
  real :: atot ! total fraction of tiles that can be involved in transitions
  real :: htot ! total fraction heat, for debugging only
  ! variable used for conservation check:
  real :: lmass0, fmass0, cmass0, heat0, &
       soil_heat0, vegn_heat0, cana_heat0, snow_heat0 ! pre-transition values
  real :: lmass1, fmass1, cmass1, heat1, &
       soil_heat1, vegn_heat1, cana_heat1, snow_heat1 ! post-transition values
  real :: lm, fm ! buffers for transition calculations

  ! conservation check code, part 1: calculate the pre-transition grid
  ! cell totals
  lmass0 = 0 ; fmass0 = 0 ; cmass0 = 0 ; heat0 = 0
  soil_heat0 = 0 ;  vegn_heat0 = 0 ; cana_heat0 = 0 ; snow_heat0 = 0
  ts = first_elmt(d_list)
  do while (loop_over_tiles(ts, ptr))
     call get_tile_water(ptr,lm,fm)
     lmass0 = lmass0 + lm*ptr%frac ; fmass0 = fmass0 + fm*ptr%frac

     heat0  = heat0  + land_tile_heat  (ptr)*ptr%frac
     cmass0 = cmass0 + land_tile_carbon(ptr)*ptr%frac

     if(associated(ptr%soil)) soil_heat0 = soil_heat0 + soil_tile_heat(ptr%soil)*ptr%frac
     if(associated(ptr%vegn)) vegn_heat0 = vegn_heat0 + vegn_tile_heat(ptr%vegn)*ptr%frac
     if(associated(ptr%cana)) cana_heat0 = cana_heat0 + cana_tile_heat(ptr%cana)*ptr%frac
     if(associated(ptr%snow)) snow_heat0 = snow_heat0 + ptr%snow%snow_tile_heat()*ptr%frac ! EZSNOW
  enddo

  ! calculate the area that can participate in land transitions
  atot = 0 ; ts = first_elmt(d_list)
  do while (loop_over_tiles(ts,ptr))
     if (associated(ptr%vegn)) atot = atot + ptr%frac
  enddo

  if (is_watch_cell()) then
     write(*,*)'### land_transitions_0d: input parameters ###'
     do i = 1, size(d_kinds)
        write(*,'(i2.2,2x)', advance='no') i
        call dpri('from LU',landuse_name(d_kinds(i)))
        call dpri('to LU',  landuse_name(a_kinds(i)))
        call dpri('frac',   landuse_name(area(i)))
!         __DEBUG4__(i,d_kinds(i),a_kinds(i),area(i))
        write(*,*)
     enddo

     write(*,*)'### land_transitions_0d: land fractions before transitions (initial state) ###'
     ts = first_elmt(d_list); htot = 0.0; k = 1
     do while (loop_over_tiles(ts,ptr))
        if (associated(ptr%vegn)) then
            write(*,'(i2.2,2x)', advance='no') k; k = k+1
            call dpri('LU',landuse_name(ptr%vegn%landuse))
            call dpri('frac',ptr%frac)
            call dpri('heat',vegn_tile_heat(ptr%vegn))
            call dpri('heat*frac',vegn_tile_heat(ptr%vegn)*ptr%frac)
            write(*,*)
            htot = htot+vegn_tile_heat(ptr%vegn)*ptr%frac
        endif
     enddo
     call dpri('total area=',atot)
     call dpri('total heat=',htot)
     write(*,*)
  endif

  ! split each donor tile and gather the parts that undergo a
  ! transition into a separate list. Note that the kind of the landuse is
  ! changed during this transition, including forest harvesting if necessary.
  ! This has to occur at some time before the tiles are merged, and it seems
  ! to be the most convenient place as both original and final landuse kind
  ! is known for each part.
  call land_tile_list_init(a_list)
  select case (tran_distr_opt)
  case (DISTR_LM3)
     do i = 1,size(d_kinds)
        call split_changing_tile_parts(d_list,d_kinds(i),a_kinds(i),area(i)*atot,a_list)
        ! the factor atot normalizes the transitions to the total area in the grid cell
        ! available for the land use, that is, the area of land excluding lakes and glaciers
     enddo
  case (DISTR_MIN)
     ! d_kinds and a_kinds are the arrays of initial and final LU types for each of
     ! the transitions. The arrays are of equal size. For each initial and final
     ! LU types src and dst, there is only one element src->dst in these arrays.
     !
     ! We go in order (U,C,P,R,S) through the final LU types, and apply all transitions
     ! that convert land to this type. Since initial type for each of these transitions
     ! are different, there should not be dependence on the order of operations.
     !
     ! An alternative algorithm would be to arrange d_kinds, a_kinds, and area in
     ! the above order (x->U, x->C, x->P, x->S for any x), and go through the
     ! arranged array.
     do k = 1,size(tran_order)
        do i = 1,size(a_kinds)
           if (a_kinds(i)==tran_order(k)) then
              call split_changing_tile_parts_by_priority( &
                         d_list,d_kinds(i),a_kinds(i),area(i)*atot,a_list)
           endif
        enddo
     enddo
  end select
  if (is_watch_cell()) then
     write(*,*)'### land_transitions_0d: land fractions after splitting changing parts ###'
     atot = 0 ; ts = first_elmt(d_list); k = 1
     do while (loop_over_tiles(ts,ptr))
        if (.not.associated(ptr%vegn)) cycle
        write(*,'(i2.2,2x)', advance='no') k; k = k+1
        call dpri('donor LU',landuse_name(ptr%vegn%landuse))
        call dpri('frac',ptr%frac)
        write(*,*)
!         write(*,'(2(a,g23.16,2x))')'   donor: LU = '//landuse_name(ptr%vegn%landuse),' frac=',ptr%frac
        atot = atot + ptr%frac
     enddo
     ts = first_elmt(a_list); k = 1
     do while (loop_over_tiles(ts, ptr))
        if (.not.associated(ptr%vegn)) cycle
        write(*,'(i2.2,2x)', advance='no') k; k = k+1
        call dpri('acceptor LU',landuse_name(ptr%vegn%landuse))
        call dpri('frac',ptr%frac)
        write(*,*)
!         write(*,'(2(a,g23.16,2x))')'acceptor: LU = '//landuse_name(ptr%vegn%landuse),' frac=',ptr%frac
        atot = atot + ptr%frac
     enddo
     call dpri('total area=',atot)
     write(*,*)
  endif

  ! move all tiles from the donor list to the acceptor list -- this will ensure
  ! that all the tiles that can be merged at this time will be
  te = tail_elmt(d_list)
  do
     ts=first_elmt(d_list)
     if(ts==te) exit ! reached the end of the list
     ptr=>current_tile(ts)
     if(ptr%frac <= 0.0) then
        call erase(ts) ! if area of the tile is zero, free it
     else
        ! otherwise, move it to a_list
        call remove(ts)
        call insert(ptr,a_list)
     endif
  enddo
  ! d_list is empty at this point

  ! merge all generated tiles into the source (donor) list
  te = tail_elmt(a_list)
  do
     ts=first_elmt(a_list)
     if(ts==te) exit ! break out of loop
     ptr=>current_tile(ts)
     call remove(ts)
     if (ptr%frac > 0.0) then
         call merge_land_tile_into_list(ptr,d_list)
     else
         call delete_land_tile(ptr)
     endif
  enddo
  ! a_list is empty at this point
  call land_tile_list_end(a_list)

  if (is_watch_cell()) then
     write(*,*)'### land_transitions_0d: land fractions final state ###'
     ts = first_elmt(d_list); htot = 0.0; k=0
     do while (loop_over_tiles(ts,ptr))
        if (associated(ptr%vegn)) then
            write(*,'(i2.2,2x)', advance='no') k; k = k+1
            call dpri('LU',landuse_name(ptr%vegn%landuse))
            call dpri('frac',ptr%frac)
            call dpri('heat',vegn_tile_heat(ptr%vegn))
            call dpri('heat*frac',vegn_tile_heat(ptr%vegn)*ptr%frac)
            write(*,*)
            htot = htot+vegn_tile_heat(ptr%vegn)*ptr%frac
        endif
     enddo
     call dpri('total area=',atot)
     call dpri('total heat=',htot)
     write(*,*)
  endif

  ! conservation check part 2: calculate grid cell totals in final state, and
  ! compare them with pre-transition totals
  lmass1 = 0 ; fmass1 = 0 ; cmass1 = 0 ; heat1 = 0
  soil_heat1 = 0 ;  vegn_heat1 = 0 ; cana_heat1 = 0 ; snow_heat1 = 0
  ts = first_elmt(d_list)
  do while (loop_over_tiles(ts,ptr))
     call get_tile_water(ptr,lm,fm)
     lmass1 = lmass1 + lm*ptr%frac ; fmass1 = fmass1 + fm*ptr%frac

     heat1  = heat1  + land_tile_heat  (ptr)*ptr%frac
     cmass1 = cmass1 + land_tile_carbon(ptr)*ptr%frac

     if(associated(ptr%soil)) soil_heat1 = soil_heat1 + soil_tile_heat(ptr%soil)*ptr%frac
     if(associated(ptr%vegn)) vegn_heat1 = vegn_heat1 + vegn_tile_heat(ptr%vegn)*ptr%frac
     if(associated(ptr%cana)) cana_heat1 = cana_heat1 + cana_tile_heat(ptr%cana)*ptr%frac
     if(associated(ptr%snow)) snow_heat1 = snow_heat1 + ptr%snow%snow_tile_heat()*ptr%frac ! EZSNOW
  enddo
    ! EZSNOW //FIXME: I have temporarily removed checks as snow merging tiles can chance ice and water, but not their total [not currently used]
  call check_conservation ('liquid + frozen water', lmass0+fmass0, lmass1+fmass1, 1e-6) ! EZSNOW
!   call check_conservation ('liquid water', lmass0, lmass1, 1e-6)
!   call check_conservation ('frozen water', fmass0, fmass1, 1e-6)
  call check_conservation ('carbon'      , cmass0, cmass1, 1e-6)
  call check_conservation ('canopy air heat content', cana_heat0 , cana_heat1 , 1e-6)
! heat content of vegetation may not conserve because of the cohort merging issues
!  call check_conservation ('vegetation heat content', vegn_heat0 , vegn_heat1 , 1e-6)
  call check_conservation ('snow heat content',       snow_heat0 , snow_heat1 , 1e-6)
  call check_conservation ('soil heat content',       soil_heat0 , soil_heat1 , 1e-4)
  call check_conservation ('heat content', heat0 , heat1 , 1e-4)

end subroutine land_transitions_0d

! =============================================================================
! check that the requested area of transitions is not larger than available area
! in tiles
subroutine check_area_overshoot(area, d_kind, a_kind, dfrac)
  real,    intent(in) :: area   ! total area of donor tiles
  integer, intent(in) :: d_kind ! LU type of donor tiles
  integer, intent(in) :: a_kind ! LU type of acceptor tiles
  real,    intent(in) :: dfrac  ! fraction of land area that changes LU type

  integer :: severity ! severity of overshoot errors
  integer :: i,j,k,face ! coordinates of current point, for overshoot diagnostics

  ! check for overshoot situation: that is, a case where the transition area is
  ! larger than the available area
  if(overshoot_opt /= OPT_IGNORE.and.dfrac>area+overshoot_tolerance) then
     severity = WARNING
     if (overshoot_opt==OPT_STOP) severity = FATAL
     call get_current_point(i,j,k,face)
     call error_mesg('landuse',&
          'transition at ('//trim(string(i))//','//trim(string(j))//&
          ',face='//trim(string(face))//&
          ') from "'//trim(landuse_name(d_kind))// &
          '" to "'  //trim(landuse_name(a_kind))//&
          '" ('//trim(string(dfrac))//') is larger than area of "'&
          //trim(landuse_name(d_kind))//'" ('//trim(string(area))//')', &
          severity)
  endif
end subroutine check_area_overshoot

! =============================================================================
! splits changing parts of donor tiles into a separate tile list, performing
! land use changes in the process
subroutine split_changing_tile_parts_by_priority(d_list,d_kind,a_kind,dfrac,a_list)
  type(land_tile_list_type), intent(in) :: d_list ! list of donor tiles
  integer, intent(in) :: d_kind ! LU type of donor tiles
  integer, intent(in) :: a_kind ! LU type of acceptor tiles
  real,    intent(in) :: dfrac  ! fraction of land area that changes LU type
  type(land_tile_list_type), intent(inout) :: a_list ! list of acceptors

  ! ---- local vars
  type(land_tile_enum_type) :: ct
  type(land_tile_type), pointer :: tile, temp
  real :: area, darea, tfrac
  real,    allocatable :: priority(:) ! priority of the land use transition fro each tile
  integer, allocatable :: idx(:)      ! array of tile indices in the descending priority order
  integer :: k
  integer :: ntiles ! number of tiles in d_list

  ! calculate total area of the tiles that should be transitioned to another kind
  ct = first_elmt(d_list); area = 0.0
  do while (loop_over_tiles(ct, tile))
     if (.not.associated(tile%vegn)) cycle
     if (tile%vegn%landuse == d_kind) area = area + tile%frac
  enddo

  call check_area_overshoot(area,d_kind,a_kind,dfrac)

  ! calculate transition priorities
  ntiles = nitems(d_list)
  allocate(priority(ntiles), idx(ntiles))
  priority(:) = -HUGE(1.0)
  k = 0; ct = first_elmt(d_list)
  do while (loop_over_tiles(ct,tile))
     k = k+1
     if(.not.associated(tile%vegn))  cycle ! skip non-vegetated tiles
     if(tile%vegn%landuse /= d_kind) cycle ! skip tiles that do not match donor LU type
     priority(k) = landuse_priority(tile, a_kind)
  enddo

  ! sort landuse transition priorities in descending order
  call rank_descending(priority, idx)

  ! transition cannot be more than current total area of specified kind
  tfrac = min(dfrac,area)
  do k = 1, ntiles
     if (tfrac==0) exit ! from loop, no more area to transition
     tile=>elmt_at_index(d_list, idx(k))
     if (.not.associated(tile%vegn)) cycle ! landuse cannot be applied to non-vegetated tiles
     if(tile%vegn%landuse /= d_kind) cycle ! skip tiles that do not match donor LU type
     darea = min(tile%frac, tfrac)
     if (darea>0) then
        ! make a copy of current tile
        temp => new_land_tile(tile)
        temp%frac = darea
        tile%frac = tile%frac-darea
        ! convert land use type of the tile: cut the forest, if necessary
        if( temp%vegn%landuse==LU_NTRL.or.  &
            temp%vegn%landuse==LU_SCND.or.  &
            temp%vegn%landuse==LU_RANGE.or. &
           ((temp%vegn%landuse/=LU_CROP.and.a_kind==LU_CROP).and.clear_all_on_conversion_to_crop) &
          ) then
           call vegn_cut_forest(temp, a_kind)
        endif
        ! change landuse type of the tile
        temp%vegn%landuse = a_kind
        ! reset time elapsed since last disturbance and time elapsed since last land use
        ! event in the new tile
        temp%vegn%age_since_disturbance = 0.0
        temp%vegn%age_since_landuse     = 0.0
        ! add the new tile to the resulting list
        call insert(temp, a_list) ! insert tile into output list
!       call debug_crop(temp%vegn,'transition from "'//landuse_name(tile%vegn%landuse)//'"')
        ! calculate remaining area of transition
        tfrac = tfrac-darea
     endif
  enddo

end subroutine split_changing_tile_parts_by_priority

! ============================================================================
! returns priority of the land use tile: tiles with highest number will be
! consumed first by the land use transition
function landuse_priority(tile, dst) result(P); real P
  type(land_tile_type), intent(in) :: tile
  integer, intent(in) :: dst ! land use types we are transitioning to

  integer :: src ! land use type of the tile

  P = -HUGE(1.0) ! very low priority
  if (.not.associated(tile%vegn)) return
  src = tile%vegn%landuse

  if ((src==LU_SCND.or.src==LU_NTRL).and.dst==LU_SCND) then
     ! for wood harvesting (NTRL->SCND or SCND->SCND), first
     ! consume tiles with highest wood biomass
     P = vegn_tile_bwood(tile%vegn)
  else if (dst==LU_SCND) then
     ! for abandonment, we first consume top-of-the-hill tiles
     ! hidx_j is the index of the hillslope tile; the higher the index the
     ! higher the tile in the hillslope
     P = tile%soil%hidx_j
  else if (src==LU_CROP.and.dst==LU_PAST) then
     ! for CROP->PAST conversion, start from the top of the hill
     P = tile%soil%hidx_j
  else
     ! for everything else, start from the bottom
     P = -tile%soil%hidx_j
  endif
end function landuse_priority

! =============================================================================
! splits changing parts of donor tiles into a separate tile list, performing
! land use changes in the process
subroutine split_changing_tile_parts(d_list,d_kind,a_kind,dfrac,a_list)
  type(land_tile_list_type), intent(in) :: d_list ! list of donor tiles
  integer, intent(in) :: d_kind ! LU type of donor tiles
  integer, intent(in) :: a_kind ! LU type of acceptor tiles
  real,    intent(in) :: dfrac  ! fraction of land area that changes LU type
  type(land_tile_list_type), intent(inout) :: a_list ! list of acceptors

  ! ---- local vars
  type(land_tile_enum_type) :: ct
  type(land_tile_type), pointer :: tile, temp
  real :: area, darea, area0, area1
  real :: x0,x1,x2 ! values of transition intensity
  real, parameter :: eps = 1e-6 ! area calculation precision
  real, parameter :: factor = 1.6 ! multiplier for solution bracketing
  integer :: iter

  ! calculate total area of the tiles that should be transitioned to another kind
  ct = first_elmt(d_list); area = 0.0
  do while (loop_over_tiles(ct, tile))
     if (.not.associated(tile%vegn)) cycle
     if (tile%vegn%landuse == d_kind)  &
          area = area + tile%frac
  enddo

  call check_area_overshoot(area,d_kind,a_kind,dfrac)

  ! if area of the tiles of requested kind is zero we cannot transition
  ! anything, so just return
  if (area==0) return

  ! transition cannot be more than current total area of specified kind
  darea = min(dfrac, area)

  ! solve equation to get transition intensity
  ! (1) bracket transition intensity interval so that requested area is within it
  x0=0.0; area0 = total_transition_area(d_list, d_kind, a_kind, x0)
  x1=1.0; area1 = total_transition_area(d_list, d_kind, a_kind, x1)
  iter = 0
  do
     if ((area0<=darea).and.(area1>=darea)) exit
     if (area0>darea) then
        x0 = x0-(x1-x0)*factor
        area0 = total_transition_area(d_list, d_kind, a_kind, x0)
     else
        x1 = x1+(x1-x0)*factor
        area1 = total_transition_area(d_list, d_kind, a_kind, x1)
     endif
     iter = iter+1
     if (iter>50) then
        call error_mesg('veg_tile_transitions',&
             'cannot braket transition intensity interval after 50 iterations',&
             FATAL)
     endif
  enddo

  ! find solution for transition intensity by binary search
  do iter = 1,50
     x2 = (x0+x1)/2
     area = total_transition_area(d_list, d_kind, a_kind, x2)
     if (abs(x1-x2)<eps) exit
     if (area>darea) then
        x1=x2
     else
        x0=x2
     endif
  enddo

  ! do tile transitions to destination list
  ct = first_elmt(d_list)
  do while (loop_over_tiles(ct, tile))
     if(.not.associated(tile%vegn))  cycle ! skip all non-vegetation tiles
     if(tile%vegn%landuse /= d_kind) cycle ! skip all tiles that do not match "donor" LU kind
     darea = vegn_tran_priority(tile%vegn, a_kind, x2)
     if(tile%frac*darea > 0) then
        ! make a copy of current tile
        temp => new_land_tile(tile)
        temp%frac = tile%frac*darea
        tile%frac = tile%frac*(1.0-darea)
        ! convert land use type of the tile: cut the forest, if necessary
        if( temp%vegn%landuse==LU_NTRL.or.  &
            temp%vegn%landuse==LU_SCND.or.  &
            temp%vegn%landuse==LU_RANGE.or. &
           ((temp%vegn%landuse/=LU_CROP.and.a_kind==LU_CROP).and.clear_all_on_conversion_to_crop) &
          ) then
           call vegn_cut_forest(temp, a_kind)
        endif
        ! change landuse type of the tile
        temp%vegn%landuse = a_kind
        ! reset time elapsed since last disturbance and time elapsed since last land use
        ! event in the new tile
        temp%vegn%age_since_disturbance = 0.0
        temp%vegn%age_since_landuse     = 0.0

!       call debug_crop(temp%vegn,'transition from "'//landuse_name(tile%vegn%landuse)//'"')
        ! add the new tile to the resulting list
        call insert(temp, a_list) ! insert tile into output list
     endif
  enddo

end subroutine split_changing_tile_parts


! ============================================================================
! calculates total area (fraction of grid cell area) participating in
! vegetation transition from src_kind to dst_kind for given transition
! intensity tau
function total_transition_area(list,src_kind,dst_kind,tau) result (total_area)
  real :: total_area
  type(land_tile_list_type), intent(in) :: list ! list of tiles
  integer , intent(in) :: src_kind, dst_kind ! source and destination kinds
  real    , intent(in) :: tau                ! transition intensity

  ! ---- local vars
  type(land_tile_enum_type) :: ct
  type(land_tile_type), pointer :: tile

  total_area = 0
  ct = first_elmt(list)
  do while (loop_over_tiles(ct, tile))
     if (.not.associated(tile%vegn)) cycle ! skip non-vegetated tiles
     if(tile%vegn%landuse == src_kind) &
          total_area = total_area + tile%frac*vegn_tran_priority(tile%vegn,dst_kind,tau)
  enddo

end function total_transition_area


! ============================================================================
! given a vegetation patch, destination kind of transition, and "transition
! intensity" value, this function returns a fraction of tile that will parti-
! cipate in transition.
!
! this function must be contiguous, monotonic, its value must be within
! interval [0,1]
!
! this function is used to determine what part of each tile is to be converted
! to another land use kind; the equation is solved to get "transition intensity"
! tau for which total area is equal to requested. Tau is, therefore, a dummy
! parameter, and only relative values of the priority functions for tiles
! participating in transition have any meaning. For most transitions the priority
! function is just equal to tau: therefore there is no preference, and all tiles
! contribute equally to converted area. For secondary vegetation harvesting,
! however, priority also depends on wood biomass, and therefore tiles
! with high wood biomass are harvested first.
function vegn_tran_priority(vegn, dst_kind, tau) result(P); real :: P
  type(vegn_tile_type), intent(in) :: vegn
  integer             , intent(in) :: dst_kind
  real                , intent(in) :: tau

  real :: vegn_bwood

  if (vegn%landuse==LU_SCND.and.dst_kind==LU_SCND) then ! secondary biomass harvesting
     vegn_bwood = vegn_tile_bwood(vegn)
     P = max(min(tau+vegn_bwood,1.0),0.0)
  else
     P = max(min(tau,1.0),0.0)
  endif
end function vegn_tran_priority


! ============================================================================

subroutine add_to_transitions(frac, time0,time1,k1,k2,tran)
  real, intent(in) :: frac(lnd%ls:lnd%le)
  type(time_type), intent(in) :: time0       ! time of previous calculation of
    ! transitions (the integral transitions will be calculated between time0
    ! and time)
  type(time_type), intent(in) :: time1       ! current time
  integer, intent(in) :: k1,k2               ! kinds of tiles
  type(tran_type), pointer :: tran(:,:)    ! transition info

  ! ---- local vars
  integer :: k,sec,days,l
  type(tran_type), pointer :: ptr(:,:) => NULL()
  real    :: part_of_year
  logical :: used

  ! allocate array of transitions, if necessary
  if (.not.associated(tran)) allocate(tran(lnd%ls:lnd%le,1))

  do l = lnd%ls, lnd%le
     if(frac(l) == 0) cycle ! skip points where transition rate is zero
     ! find the first empty transition element for the current indices
     k = 1
     do while ( k <= size(tran,2) )
        if(tran(l,k)%donor == 0) exit
        k = k+1
     enddo

     if (k>size(tran,2)) then
        ! if there is no room, make the array of transitions larger
        allocate(ptr(lnd%ls:lnd%le,size(tran,2)*2))
        ptr(:,1:size(tran,2)) = tran
        deallocate(tran)
        tran => ptr
        nullify(ptr)
     end if

     ! store the transition element
     tran(l,k) = tran_type(k1,k2,frac(l))
  enddo

  ! send transition data to diagnostics
  if(diag_ids(k1,k2)>0) then
     call get_time(time1-time0, sec,days)
     part_of_year = (days+sec/86400.0)/days_in_year(time0)
     used = send_data(diag_ids(k1,k2), frac/part_of_year, time1)
  endif

end subroutine add_to_transitions


! ==============================================================================
! given boundaries of time interval [t1,t2], calculates total transition (time
! integral of transition rates) over the specified interval
subroutine integral_transition(t1, t2, tran, frac, err_msg)
  type(time_type), intent(in)  :: t1,t2 ! time boundaries
  type(varset_T),  intent(in)  :: tran ! id of the field
  real           , intent(out) :: frac(:)
  character(len=*),intent(out), optional :: err_msg

  ! ---- local vars
  integer :: n ! size of time axis
  type(time_type) :: ts,te
  integer         :: i1,i2
  real :: w  ! time interpolation weight
  real :: dt ! current time interval, in years
  real :: sum(size(frac(:)))
  integer :: l
  character(len=256) :: msg

  msg = ''
  ! adjust the integration limits, in case they are out of range
  associate(time_in => tran%file%time_in)
  n = size(time_in)
  ts = t1;
  if (ts<time_in(1)) ts = time_in(1)
  if (ts>time_in(n)) ts = time_in(n)
  te = t2
  if (te<time_in(1)) te = time_in(1)
  if (te>time_in(n)) te = time_in(n)

  call time_interp(ts, time_in, w, i1,i2, err_msg=msg)
  if(msg /= '') then
    if(fms_error_handler('integral_transition','Message from time_interp: '//trim(msg),err_msg)) return
  endif
  call tran%get_data(i1,frac)

  dt = (time_in(i2)-time_in(i1))//set_time(0,days_in_year((time_in(i2)+time_in(i1))/2))
  sum = -frac*w*dt
  do while(time_in(i2)<=te)
     call tran%get_data(i1,frac)
     dt = (time_in(i2)-time_in(i1))//set_time(0,days_in_year((time_in(i2)+time_in(i1))/2))
     sum = sum+frac*dt
     i2 = i2+1
     i1 = i2-1
     if(i2>size(time_in)) exit ! from loop
  enddo

  call time_interp(te,time_in,w,i1,i2, err_msg=msg)
  if(msg /= '') then
    if(fms_error_handler('integral_transition','Message from time_interp: '//trim(msg),err_msg)) return
  endif
  call tran%get_data(i1,frac)
  dt = (time_in(i2)-time_in(i1))//set_time(0,days_in_year((time_in(i2)+time_in(i1))/2))
  frac = sum+frac*w*dt
  end associate
  ! check the transition rate validity
  do l = 1,size(frac(:))
     call set_current_point(l+lnd%ls-1,1)
     call check_var_range(frac(l),0.0,HUGE(1.0),'integral_transition',tran%name, FATAL)
  enddo
end subroutine integral_transition


! ==============================================================================
! checks conservation and aborts with fatal error if tolerance is exceeded
subroutine check_conservation(name, d1, d2, tolerance)
  character(*), intent(in) :: name ! name of the component
  real, intent(in) :: d1,d2 ! values to check
  real, intent(in) :: tolerance ! tolerance of the test

  integer :: curr_i, curr_j, face
  integer :: severity ! severity of the generated message
  character(256) :: message

  if (conservation_opt == OPT_IGNORE) return ! do nothing

  severity = WARNING
  if (conservation_opt==OPT_STOP) severity = FATAL

  if (abs(d1-d2)>tolerance) then
     call get_current_point(i=curr_i,j=curr_j,face=face)
     write(message,'(a,3(x,a,i4), 3(x,a,g23.16))')&
          'conservation of '//trim(name)//' is violated', &
          'at i=',curr_i,'j=',curr_j,'face=',face, &
          'value before=', d1, 'after=', d2, 'diff=',d2-d1
     call error_mesg('land_transitions',message,severity)
  endif
end subroutine check_conservation

end module
