module land_transitions_mod
#include <fms_platform.h>

#include "../shared/debug.inc"

use constants_mod, only : PI

use mpp_io_mod, only : mpp_open, mpp_close, MPP_ASCII, MPP_RDONLY

use fms_mod, only : string, error_mesg, FATAL, WARNING, NOTE, &
     mpp_pe, lowercase, file_exist, close_file, read_data, &
     check_nml_error, stdlog, mpp_root_pe, fms_error_handler
use fms_io_mod, only : get_file_name

use time_manager_mod, only : time_type, set_date, get_date, set_time, &
     operator(+), operator(-), operator(>), operator(<), operator(<=), operator(/), &
     operator(//), operator(==), days_in_year, print_date, increment_date, get_time, &
     valid_calendar_types, get_calendar_type
use get_cal_time_mod, only : get_cal_time
use horiz_interp_mod, only : horiz_interp_type, horiz_interp_init, &
     horiz_interp_new, horiz_interp_del
use time_interp_mod, only : time_interp
use diag_manager_mod, only : register_diag_field, send_data, diag_field_add_attribute

use nfu_mod, only : nfu_validtype, nfu_inq_var, nfu_get_dim_bounds, nfu_get_rec, &
     nfu_get_dim, nfu_get_var, nfu_get_valid_range, nfu_is_valid

use vegn_data_mod, only : &
     N_LU_TYPES, M_LU_TYPES, LU_PAST, LU_RAINF, LU_IRRIG, LU_NTRL, LU_SCND, LU_RANGE, LU_URBN, &
     landuse_name, landuse_longname

use cana_tile_mod, only : cana_tile_heat, cana_tile_stock_pe, cana_tile_carbon, new_cana_tile, merge_cana_tiles, cana_tile_type, delete_cana_tile
use snow_tile_mod, only : snow_tile_heat, snow_tile_stock_pe
use vegn_tile_mod, only : vegn_tile_heat, vegn_tile_stock_pe, vegn_tile_carbon, vegn_tile_type, vegn_tile_bwood
use soil_tile_mod, only : soil_tile_heat, soil_tile_stock_pe, soil_tile_carbon
use lake_tile_mod, only : lake_tile_heat, lake_tile_type
use glac_tile_mod, only : glac_tile_heat

use land_tile_mod, only : land_tile_map, &
     land_tile_type, land_tile_list_type, land_tile_enum_type, new_land_tile, delete_land_tile, &
     first_elmt, tail_elmt, loop_over_tiles, operator(==), current_tile, &
     land_tile_list_init, land_tile_list_end, nitems, elmt_at_index, &
     erase, remove, insert, merge_land_tile_into_list, &
     get_tile_water, land_tile_carbon, land_tile_heat, &
     empty
use land_tile_io_mod, only : print_netcdf_error
use land_tile_diag_mod, only : cmor_name

use land_data_mod, only : lnd, log_version, horiz_interp_ug
use vegn_harvesting_mod, only : vegn_cut_forest

use land_debug_mod, only : set_current_point, is_watch_cell, &
     get_current_point, check_var_range, log_date, land_error_message
use land_numerics_mod, only : rank_descending
use lake_mod, only : prohibit_shallow_lake, is_rsv_restart, use_reservoir
use transitions_input_mod
use transition_io_mod, only : transition_io_init, infile_T, varset_T, &
     new_infile_LUH1, new_infile_LUH2, new_infile_CS

implicit none
private

! ==== public interface =====================================================
public :: land_transitions_init
public :: land_transitions_end
public :: save_land_transitions_restart

public :: lake_transitions_init
public :: lake_transitions_end
public :: save_lake_transitions_restart

public :: land_transitions
public :: lake_transitions
! TODO: check the role of re-exports for conflicts with the rest of the transitions
public :: do_landuse_change ! re-exporting from transitions_input_mod
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
integer, parameter :: tran_order(M_LU_TYPES) = [ LU_URBN, LU_RAINF, LU_IRRIG, LU_PAST, LU_RANGE, LU_SCND, LU_NTRL ]

! TODO: describe differences between data sets

! ==== module data ==========================================================
logical :: module_is_initialized = .FALSE.

class(infile_T), pointer :: ftran=>NULL(), fstate=>NULL(), firrig=>NULL()
type(varset_T) :: input_tran  (N_LU_TYPES,N_LU_TYPES) ! input transition rate fields
type(varset_T) :: input_state (N_LU_TYPES,N_LU_TYPES) ! input state field (for initial transition only)

integer :: diag_ids  (N_LU_TYPES,N_LU_TYPES)
type(time_type) :: time0 ! time of previous transition calculations

integer :: tran_distr_opt = -1 ! selector for transition distribution option, for efficiency
integer :: overshoot_opt = -1 ! selector for overshoot handling options, for efficiency
integer :: conservation_opt = -1 ! selector for non-conservation handling options, for efficiency

! translation of luh2 names and LM3 landuse types
character(5) :: luh2name(12)
integer      :: luh2type(12)
integer :: idata
data (luh2name(idata), luh2type(idata), idata = 1, 12) / &
   'primf', LU_NTRL, &
   'primn', LU_NTRL, &
   'secdf', LU_SCND, &
   'secdn', LU_SCND, &
   'urban', LU_RAINF, &
   'c3ann', LU_RAINF, &
   'c4ann', LU_RAINF, &
   'c3per', LU_RAINF, &
   'c4per', LU_RAINF, &
   'c3nfx', LU_RAINF, &
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
! do we need to add irrigation here?
integer, parameter :: lu2lumip(N_LU_TYPES) = [LUMIP_PST, LUMIP_CRP, LUMIP_PSL, LUMIP_PSL, LUMIP_URB, LUMIP_PST]

! variables for irrigation
type(varset_T) :: input_irrig!, input_flood(1) ! input irrigation area
type(varset_T) :: input_crop  ! input crop area
! variables for reservoir
logical :: module_is_initialized_lake = .FALSE.
type(time_type) :: timel0 ! time of previous lake transition calculations
class(infile_T), pointer :: infile_lake_tran  => NULL()
class(infile_T), pointer :: infile_lake_state => NULL()
class(infile_T), pointer :: infile_depth_rsv  => NULL()
type(varset_T) :: input_tran_lake  (2,2) ! input transition rate fields
type(varset_T) :: input_state_lake (2,2) ! input state field (for initial transition only)
type(varset_T) :: input_depth_rsv
character(len=5), parameter  :: &
     landuse_name_lake (2) = (/ 'lake','soil'/)
real :: rsv_depth_min = 2.
character(len=1024) :: input_lake_file  = '' ! input data set of lake transition dates
character(len=1024) :: state_lake_file  = '' ! input data set of LU states (for initial transition only)
character(len=1024) :: depth_rsv_file   = '' ! reservoir construction depth

contains ! ###################################################################

! ============================================================================
subroutine land_transitions_init(id_ug, id_cellarea)
  integer, intent(in) :: id_ug !<Unstructured axis id.
  integer, intent(in) :: id_cellarea !<id of cell area diagnostic fields

  ! ---- local vars
  integer        :: unit
  integer        :: year,month,day,hour,min,sec
  integer        :: k1,k2,k3, n1,n2
  character(12)  :: fieldname

  type(land_tile_type), pointer :: tile
  type(land_tile_enum_type) :: ce

  if(module_is_initialized) return
  module_is_initialized = .TRUE.
  call log_version(version, module_name, &
  __FILE__)

  call horiz_interp_init()

  ! read restart file, if any
  if (file_exist('INPUT/landuse.res')) then
     call error_mesg('land_transitions_init','reading restart "INPUT/landuse.res"',&
          NOTE)
     call mpp_open(unit,'INPUT/landuse.res', action=MPP_RDONLY, form=MPP_ASCII)
     read(unit,*) year,month,day,hour,min,sec
     time0 = set_date(year,month,day,hour,min,sec)
     call mpp_close(unit)
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
          trim(overshoot_handling)//'" is incorrect, use "stop", "report", or "ignore"',&
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
          trim(conservation_handling)//'" is incorrect, use "stop", "report", or "ignore"',&
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
          'current time is earlier than the time of last land use transition application',&
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

  ! initialize arrays of input fields
  select case (trim(lowercase(data_type)))
  case('luh1')
     ftran => new_infile_LUH1(input_file)
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
     ftran => new_infile_LUH2(input_file,static_file)
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

     if (time0==set_date(0001,01,01)) then
        call error_mesg('land_transitions_init','setting up initial land use transitions', NOTE)
        ! initialize state input for initial transition from all-natural state.
        if (trim(state_file)=='') call error_mesg('land_transitions_init',&
            'starting land use transitions, but land use state file is not specified',FATAL)

        ! open state file
        fstate=>new_infile_LUH2(state_file,static_file)

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

  if (do_irrigation) then
     if (trim(irrigation_file)=='') call error_mesg('land_transitions_init', &
          'irrigation transitions are turned on, but irrigation input file is not specified', &
          FATAL)
     if (trim(state_file)=='') call error_mesg('land_transitions_init',&
         'irrigation transitions are turned on, but land use state file is not specified',FATAL)

     ! initialize data structure representing input file and horizontal interpolator
     firrig=>new_infile_LUH2(irrigation_file,static_file)
     ! open state file, if necessary
     if (.not. associated(fstate)) &
        fstate=>new_infile_LUH2(state_file,static_file)

     ! create input variable set for irrigated fraction. Note that currently we sum up
     ! irrigation areas for all crops and use the total.
     input_irrig%name='irrigation fraction'
     input_crop %name='cropland fraction'
     do n2 = 1,size(luh2type)
        if ((luh2type(n2)==LU_RAINF).or.(luh2type(n2)==LU_IRRIG)) then
           call input_irrig%addvar(firrig,trim(luh2name(n2))//'_irrig')
           call input_crop %addvar(fstate,trim(luh2name(n2)))
        endif
     enddo
     if (mpp_pe()==mpp_root_pe()) then
        write(*,*)'land_transitions_init: summary of irrigation-related input'
        write(*,'(a)') input_irrig%descr()
        write(*,'(a)') input_crop %descr()
     endif
  endif

end subroutine land_transitions_init

!===========================================================================
subroutine lake_transitions_init(id_ug)
  integer, intent(in) :: id_ug !<Unstructured axis id.

  ! ---- local vars
  integer        :: unit
  integer        :: year,month,day,hour,mi,sec
  integer        :: k1,k2,k3, n1,n2, l
  real           :: frac(lnd%ls:lnd%le)

  type(land_tile_enum_type)     :: ce    ! land tile enumerator
  type(land_tile_type), pointer :: tile  ! pointer to current tile
  real :: frac2land, whole_lake_area, Afrac_rsv_bak
  logical :: read_dist, io_domain_exist, found_file

  if(module_is_initialized_lake) return
  module_is_initialized_lake = .TRUE.

  if (file_exist('INPUT/laketran.res')) then
     call error_mesg('lake_transitions_init','reading restart "INPUT/laketran.res"',&
          NOTE)
     call mpp_open(unit,'INPUT/laketran.res', action=MPP_RDONLY, form=MPP_ASCII)
     read(unit,*) year,month,day,hour,mi,sec
     timel0 = set_date(year,month,day,hour,mi,sec)
     call mpp_close(unit)
  else
     call error_mesg('lake_transitions_init','cold-starting lake transitions',&
          NOTE)
     timel0 = set_date(0001,01,01)
  endif

  found_file = get_file_name(state_file_lake, state_lake_file, read_dist, io_domain_exist, domain=lnd%sg_domain)
  if (file_exist(state_lake_file)) then
      infile_lake_state => new_infile_CS(input_file_lake)
      ! initialize state variable array
      ! 1 means lake, 2 means soil
      !input_state_lake(2,1)%name='initial '//trim(landuse_name_lake(2))//'2'//trim(landuse_name_lake(1))

      input_state_lake(2,1)%name='lake'
      call input_state_lake(2,1)%addvar(infile_lake_state,'lake')
  endif

  found_file = get_file_name(depth_file_rsv, depth_rsv_file, read_dist, io_domain_exist, domain=lnd%sg_domain)
  if(file_exist(depth_rsv_file))then
      ! open reservoir depth file
      infile_depth_rsv => new_infile_CS(depth_file_rsv)

      input_depth_rsv%name = 'rsv_depth'
      call input_depth_rsv%addvar(infile_depth_rsv,'rsv_depth')
  else
      if(use_reservoir) call error_mesg('lake_transitions_init','use_reservoir mod requires depth_rsv_file', FATAL)
  endif

  if(do_lake_change) then
    found_file = get_file_name(input_file_lake, input_lake_file, read_dist, io_domain_exist, domain=lnd%sg_domain)
    if (trim(input_lake_file)=='') call error_mesg('lake_transitions_init', &
       'do_lake_change is requested, but lake transition file is not specified', &
        FATAL)
    infile_lake_tran => new_infile_CS(input_file_lake)

    !do k1 = 1,2
    !do k2 = 1,2
      input_tran_lake(2,1)%name='soil_to_lake'
      !if (k1==k2) cycle
      call input_tran_lake(2,1)%addvar(infile_lake_tran,'soil_to_lake')
    !enddo
    !enddo
    if(.not.file_exist(depth_rsv_file)) &
      call error_mesg('lake_transitions_init','depth_rsv_file must exist with do_lake_change turned on', FATAL)
    if(timel0==set_date(0001,01,01).and.(.not.file_exist(state_lake_file))) &
      call error_mesg('lake_transitions_init','state_lake_file must exist when do_lake_change start', FATAL)
  endif

! initialize reservoir (rsv_depth, Afrac_rsv, Vfrac_rsv)
  if(.not.timel0==set_date(0001,01,01).and.is_rsv_restart)then
    call error_mesg('lake_transitions_init','reservoir warm start from previous laketran run', NOTE)
    call check_rsv_depth()
    return
  else if(.not.timel0==set_date(0001,01,01).and..not.is_rsv_restart)then !we have laketran.res but restart files are not complete
    call error_mesg('lake_transitions_init','restart files are not complete', FATAL)
  endif

  !timel0==set_date(0001,01,01)
  call error_mesg('lake_transitions_init', &
      'reservoir cold start, rsv_depth and Afrac_rsv are from input, though Vfrac_rsv may be from restart', NOTE)
  if (do_lake_change) then
    call error_mesg('lake_transitions_init','do_lake_change mod, abandon all restarts for rsv if there are any', NOTE)
    call rsv_set_zero()
    return
  endif

  if (use_reservoir) &
    call error_mesg('lake_transitions_init','static reservoir mod', NOTE)
  !use_reservoir%.not.do_lake_change  or  .not.use_reservoir&.not.do_lake_change
  !if Afrac_rsv and rsv_depth files exist, we read data anyway, and then determine Vfrac_rsv by restart or Afrac_rsv
  if (file_exist(state_lake_file).and.file_exist(depth_rsv_file)) then
    n1 = size(infile_depth_rsv%time_in)
    call read_rsv_depth(n1)
    frac(:) = 0.0
    n1 = size(infile_lake_state%time_in)
    call input_state_lake(2,1)%get_data(n1,frac)
    do l=lnd%ls, lnd%le
      ce = first_elmt(land_tile_map(l))
      do while(loop_over_tiles(ce,tile))
        if (.not.associated(tile%lake)) cycle
        Afrac_rsv_bak = tile%lake%Afrac_rsv
        if(lnd%ug_area(l)>0.)then
          frac2land = frac(l)*(lnd%ug_cellarea(l)/lnd%ug_area(l))
        else
          frac2land = 0.
        endif
        frac2land = max(0.,min(frac2land, tile%frac))
        tile%lake%Afrac_rsv = frac2land/tile%frac
        if(tile%lake%rsv_depth <= 0.) tile%lake%Afrac_rsv = 0.
        ! we use Vfrac_rsv from restart if there is any, otherwise set it to Afrac_rsv
        if(.not.is_rsv_restart.or.&
           (is_rsv_restart.and.Afrac_rsv_bak==0..and.tile%lake%Afrac_rsv>0.)) &
          tile%lake%Vfrac_rsv = tile%lake%Afrac_rsv
      enddo
    enddo
    if(use_reservoir) call adjust_whole_lake_area()
  else !if files are incomplete, we keep restart if there is any, otherwise set rsv to zero
    if(use_reservoir) call error_mesg('lake_transitions_init', &
                           'reservoir cold start requires both state_lake_file and depth_rsv_file', FATAL)
    if(.not.is_rsv_restart) call rsv_set_zero()
  endif
end subroutine lake_transitions_init


!===========================================================================
subroutine rsv_set_zero()
  type(land_tile_enum_type)     :: ce    ! land tile enumerator
  type(land_tile_type), pointer :: tile  ! pointer to current tile
  integer :: l

  do l=lnd%ls, lnd%le
    ce = first_elmt(land_tile_map(l))
    do while(loop_over_tiles(ce,tile))
      if (.not.associated(tile%lake)) cycle
      tile%lake%rsv_depth = 0.
      tile%lake%Afrac_rsv = 0.
      tile%lake%Vfrac_rsv = 0.
    enddo
  enddo
end subroutine rsv_set_zero

!===========================================================================
subroutine read_rsv_depth(i)
  integer, intent(in) :: i

  type(land_tile_enum_type)     :: ce    ! land tile enumerator
  type(land_tile_type), pointer :: tile  ! pointer to current tile
  integer :: l
  real, dimension(lnd%ls:lnd%le) :: rsv_depth

  rsv_depth(:) = 0.0
  call input_depth_rsv%get_data(i,rsv_depth)

  do l=lnd%ls, lnd%le
    ce = first_elmt(land_tile_map(l))
    do while(loop_over_tiles(ce,tile))
      if (.not.associated(tile%lake)) cycle
      tile%lake%rsv_depth = rsv_depth(l)
      if(tile%lake%rsv_depth > 0.) tile%lake%rsv_depth = max(rsv_depth_min,tile%lake%rsv_depth)
    enddo
  enddo
end subroutine read_rsv_depth

!===========================================================================
subroutine adjust_whole_lake_area()

  type(land_tile_enum_type)     :: ce    ! land tile enumerator
  type(land_tile_type), pointer :: tile  ! pointer to current tile
  integer :: l
  character(len=512) :: mesg ! for error message

  do l=lnd%ls, lnd%le
    ce = first_elmt(land_tile_map(l))
    do while(loop_over_tiles(ce,tile))
      if (.not.associated(tile%lake)) cycle
      tile%lake%pars%whole_area = tile%lake%pars%whole_area - tile%lake%Afrac_rsv*tile%frac*lnd%ug_area(l)
      if(tile%lake%pars%whole_area<=0.)then
        write(mesg,*)'tile%lake%pars%whole_area=',tile%lake%pars%whole_area,' tile%lake%Afrac_rsv=',tile%lake%Afrac_rsv
        call error_mesg('adjust_whole_lake_area',mesg, NOTE)
        if(tile%lake%Afrac_rsv<1.)then
          call error_mesg('adjust_whole_lake_area','lake_whole_area<0..and.Afrac_rsv<1.',FATAL)
        endif
        tile%lake%pars%whole_area = 0.
      endif
    enddo
  enddo

end subroutine adjust_whole_lake_area
!===========================================================================
subroutine check_rsv_depth()
  type(land_tile_enum_type)     :: ce    ! land tile enumerator
  type(land_tile_type), pointer :: tile  ! pointer to current tile
  integer :: l

    if (.not.use_reservoir) return
    do l=lnd%ls, lnd%le
      ce = first_elmt(land_tile_map(l))
      do while(loop_over_tiles(ce,tile))
        if (.not.associated(tile%lake)) cycle
        if(tile%lake%rsv_depth < 0.) &
          call error_mesg('check_rsv_depth','previous run must also use reservoir when this run use reservoir',FATAL)
      enddo
    enddo

end subroutine check_rsv_depth

! ============================================================================
subroutine land_transitions_end()
  module_is_initialized=.FALSE.
  ! close files and deallocate associated memory
  if (associated(ftran))  deallocate(ftran)
  if (associated(fstate)) deallocate(fstate)
  if (associated(firrig)) deallocate(firrig)
end subroutine land_transitions_end

! ============================================================================
subroutine lake_transitions_end()

  module_is_initialized_lake=.FALSE.
  if (associated(infile_lake_tran)) deallocate(infile_lake_tran)
  if (associated(infile_lake_state)) deallocate(infile_lake_state)
  if (associated(infile_depth_rsv)) deallocate(infile_depth_rsv)

end subroutine lake_transitions_end

! ============================================================================
subroutine save_land_transitions_restart(timestamp)
  character(*), intent(in) :: timestamp ! timestamp to add to the file name

  integer :: unit,year,month,day,hour,min,sec

  call mpp_open( unit, 'RESTART/'//trim(timestamp)//'landuse.res', nohdrs=.TRUE. )
  if (mpp_pe() == mpp_root_pe()) then
     call get_date(time0, year,month,day,hour,min,sec)
     write(unit,'(6i6,8x,a)') year,month,day,hour,min,sec, &
          'Time of previous landuse transition calculation'
  endif
  call mpp_close(unit)

end subroutine save_land_transitions_restart

! ============================================================================
subroutine save_lake_transitions_restart(timestamp)
  character(*), intent(in) :: timestamp ! timestamp to add to the file name

  integer :: unit,year,month,day,hour,min,sec

  call mpp_open( unit, 'RESTART/'//trim(timestamp)//'laketran.res', nohdrs=.TRUE. )
  if (mpp_pe() == mpp_root_pe()) then
     call get_date(timel0, year,month,day,hour,min,sec)
     write(unit,'(6i6,8x,a)') year,month,day,hour,min,sec, &
          'Time of previous lake transition calculation'
  endif
  call mpp_close(unit)

end subroutine save_lake_transitions_restart

! =============================================================================
subroutine land_transitions (time)
  type(time_type), intent(in) :: time

  ! ---- local vars
  integer :: k,k1,k2,l,n
  integer :: second, minute, hour, day0, day1, month0, month1, year0, year1
  integer :: i1, i2; real :: w ! indices and weight from time interpolation

  real, allocatable :: tran(:,:,:) ! (ncells, M_LU_TYPES, M_LU_TYPES) array of transitions among various land use types
  ! variables for optional irrigation transitions
  real    :: irr_area(lnd%ls:lnd%le), crop_area(lnd%ls:lnd%le)
  real    :: irr_frac(lnd%ls:lnd%le)
  real    :: area0 (M_LU_TYPES) ! fraction of each land use type before transitions
  real    :: atot ! total fraction of tiles that can be involved in transitions
  real    :: tran0 (N_LU_TYPES, N_LU_TYPES) ! array of transitions
  ! input arguments for land_transitions_0d
  integer :: src(M_LU_TYPES*M_LU_TYPES), dst(M_LU_TYPES*M_LU_TYPES) ! source and destination LU types
  real    :: frac(M_LU_TYPES*M_LU_TYPES) ! fraction of area undergoing transition
  ! variables for diagnostics
  real    :: diag(lnd%ls:lnd%le), part_of_year
  integer :: sec, days
  logical :: used

  type(land_tile_enum_type) :: ce
  type(land_tile_type), pointer :: tile

  if (.not.do_landuse_change) &
       return ! do nothing if landuse change not requested
  ! NB: in this case file/interp/data are not initialized, so it is
  ! not even possible to use the code below

  call get_date(time,             year0,month0,day0,hour,minute,second)
  call get_date(time-lnd%dt_slow, year1,month1,day1,hour,minute,second)
  if(year0 == year1) &
       return ! do nothing during a year

  if (mpp_pe()==mpp_root_pe()) &
       call log_date('land_transitions: applying land use transitions on ', time)

  ! allocate data for transitions
  allocate(tran(lnd%ls:lnd%le,M_LU_TYPES,M_LU_TYPES))

  ! get transition rates for current time: read map of transitions, and accumulate
  ! as many time steps in array of transitions as necessary. The units of transition
  ! rates are [fraction of vegetated area per year]
  tran(:,:,:) = 0.0
  ! calculate time interval, for diagnostics
  call get_time(time-time0, sec,days)
  part_of_year = (days+sec/86400.0)/days_in_year(time0)
  do k1 = 1,N_LU_TYPES
  do k2 = 1,N_LU_TYPES
     ! get transition rate for this specific transition
     if (time0==set_date(0001,01,01).and.associated(fstate)) then
        ! read initial transition from state file
        call time_interp(time, fstate%time_in, w, i1,i2)
        call input_state(k1,k2)%get_data(i1,tran(:,k1,k2))
     else
        if (input_tran(k1,k2)%nvars>0) then
           call input_tran(k1,k2)%integrate(time0,time,tran(:,k1,k2))
        endif
     endif
     if(diag_ids(k1,k2)>0) then
        used = send_data(diag_ids(k1,k2), frac/part_of_year, time)
     endif
  enddo
  enddo

  ! save the "in" and "out" diagnostics for the transitions
  do k1 = 1, N_LUMIP_TYPES
     if (id_frac_out(k1) <= 0) cycle
     diag(:) = 0.0
     do k2 = 1, M_LU_TYPES
        if(lu2lumip(k2) == k1) then
            diag(:) = diag(:) + sum(tran(:,k1,:),2)
        endif
     enddo
     used=send_data(id_frac_out(k1), diag*lnd%ug_landfrac*100.0, time)
  enddo
  do k1 = 1, N_LUMIP_TYPES
     if (id_frac_out(k1) <= 0) cycle
     diag(:) = 0.0
     do k2 = 1, M_LU_TYPES
        if(lu2lumip(k2) == k1) then
            diag(:) = diag(:) + sum(tran(:,:,k2),2)
        endif
     enddo
     used=send_data(id_frac_in(k1), diag*lnd%ug_landfrac*100.0, time)
  enddo

  if (do_irrigation) then
     ! calculate irrigated fraction of crops. Using irrigated fraction
     ! of crops instead of irrigated area allows to use irrigation data
     ! with different land use data sets that may have a different total crop
     ! area, not necessarily consistent with input irrigation areas.

     ! interpolate irrigation and crop areas from irrigation and crop data
     irr_area(:)  = 0.0
     crop_area(:) = 0.0
     call input_irrig % interpolate(time, irr_area)
     call input_crop  % interpolate(time, crop_area)
     where (crop_area > 0)
        irr_frac = irr_area/crop_area
     elsewhere
        irr_frac = 0.0
     end where
     irr_frac = min(1.0,max(0.0, irr_frac))

     do l = lnd%ls,lnd%le
        call set_current_point(l,1) ! for debug
        ! calculate areas
        atot     = 0.0  ! total area of all vegetated tiles
        area0(:) = 0.0  ! area of each of the land use types
        ce = first_elmt(land_tile_map(l))
        do while(loop_over_tiles(ce,tile))
           if (.not.associated(tile%vegn)) cycle ! skip non-vegetated tiles
           n = tile%vegn%landuse
           atot     = atot     + tile%frac
           area0(n) = area0(n) + tile%frac
        enddo

        if ((area0(LU_IRRIG).ne.0).or.(irr_frac(l).ne.0)) then
           tran0(:,:) = tran(l,1:N_LU_TYPES,1:N_LU_TYPES)
           call add_irrigation_transitions(area0(:), tran0, irr_frac(l), atot, &
                   tran(l,:,:), verbose=is_watch_cell())
        endif
     enddo
  endif ! irrigation

  ! perform the transitions
  do l = lnd%ls,lnd%le
      if(empty(land_tile_map(l))) cycle ! skip cells where there is no land
      ! set current point for debugging
      call set_current_point(l,1)
      ! assemble arrays of LU types involved in transition, and transition rates
      k = 0
      do k1 = 1,M_LU_TYPES
      do k2 = 1,M_LU_TYPES
         if (tran(l,k1,k2)<=0) cycle
         k = k+1; src(k) = k1; dst(k) = k2; frac(k) = tran(l,k1,k2)
      enddo
      enddo

      ! transition land area between different tile types
      call land_transitions_0d(land_tile_map(l), src(1:k), dst(1:k), frac(1:k))
  enddo

  ! deallocate transition array
  deallocate(tran)

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
     if(associated(ptr%snow)) snow_heat0 = snow_heat0 + snow_tile_heat(ptr%snow)*ptr%frac
  enddo

  ! calculate the area that can participate in land transitions
  atot = 0 ; ts = first_elmt(d_list)
  do while (loop_over_tiles(ts,ptr))
     if (associated(ptr%vegn)) atot = atot + ptr%frac
  enddo

  if (is_watch_cell()) then
     write(*,*)'### land_transitions_0d: input parameters ###'
     do i = 1, size(d_kinds)
        __DEBUG4__(i,d_kinds(i),a_kinds(i),area(i))
     enddo

     write(*,*)'### land_transitions_0d: land fractions before transitions (initial state) ###'
     ts = first_elmt(d_list); htot = 0.0; k = 1
     do while (loop_over_tiles(ts,ptr))
        if (associated(ptr%vegn)) then
            write(*,'(i2.2,2x)', advance='no') k; k = k+1
            call dpri('landuse',ptr%vegn%landuse)
            call dpri('area',ptr%frac)
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
     atot = 0 ; ts = first_elmt(d_list)
     do while (loop_over_tiles(ts,ptr))
        if (.not.associated(ptr%vegn)) cycle
        write(*,'(2(a,g23.16,2x))')'   donor: landuse=',ptr%vegn%landuse,' area=',ptr%frac
        atot = atot + ptr%frac
     enddo
     ts = first_elmt(a_list)
     do while (loop_over_tiles(ts, ptr))
        if (.not.associated(ptr%vegn)) cycle
        write(*,'(2(a,g23.16,2x))')'acceptor: landuse=',ptr%vegn%landuse,' area=',ptr%frac
        atot = atot + ptr%frac
     enddo
     write(*,'(a,g23.16)')'total area=',atot
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
     call merge_land_tile_into_list(ptr,d_list)
  enddo
  ! a_list is empty at this point
  call land_tile_list_end(a_list)

  if (is_watch_cell()) then
     write(*,*)'### land_transitions_0d: land fractions final state ###'
     ts = first_elmt(d_list); htot = 0.0; k=0
     do while (loop_over_tiles(ts,ptr))
        if (associated(ptr%vegn)) then
            write(*,'(i2.2,2x)', advance='no') k; k = k+1
            call dpri('landuse',ptr%vegn%landuse)
            call dpri('area',ptr%frac)
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
     if(associated(ptr%snow)) snow_heat1 = snow_heat1 + snow_tile_heat(ptr%snow)*ptr%frac
  enddo
  call check_conservation ('liquid water', lmass0, lmass1, 1e-6)
  call check_conservation ('frozen water', fmass0, fmass1, 1e-6)
  call check_conservation ('carbon'      , cmass0, cmass1, 1e-6)
  call check_conservation ('canopy air heat content', cana_heat0 , cana_heat1 , 1e-6)
! heat content of vegetation may not conserve because of the cohort merging issues
!  call check_conservation ('vegetation heat content', vegn_heat0 , vegn_heat1 , 1e-6)
  call check_conservation ('snow heat content',       snow_heat0 , snow_heat1 , 1e-6)
  call check_conservation ('soil heat content',       soil_heat0 , soil_heat1 , 1e-4)
  call check_conservation ('heat content', heat0 , heat1 , 1e-4)

end subroutine land_transitions_0d

! =============================================================================
subroutine lake_transitions (time)
  type(time_type), intent(in) :: time

  ! ---- local vars.
  integer :: k1,k2,i1,i2,l
  real, dimension(lnd%ls:lnd%le) :: frac
  integer :: second, minute, hour, day0, day1, month0, month1, year0, year1
  real    :: w
  type(land_tile_enum_type) :: ce
  type(land_tile_type), pointer :: tile
  real, dimension(lnd%ls:lnd%le) :: &
       atots, & ! total soil area in grid cell
       atotl    ! total lake area in grid cell
  logical :: is_laketran = .True.
  type(time_type) :: time_adj

  if(.not.do_lake_change) return

  call get_date(time,             year0,month0,day0,hour,minute,second)
  call get_date(time-lnd%dt_slow, year1,month1,day1,hour,minute,second)
  if(year0 == year1) &
       return ! do nothing during a year

  if (mpp_pe()==mpp_root_pe()) &
       call log_date('lake_transitions: applying lake area transitions on ', time)

  !update rsv depth
  if(use_reservoir)then
  ! adjust the integration limits, in case they are out of range
    time_adj = time_adjust(time, infile_depth_rsv%time_in)
    call time_interp(time_adj, infile_depth_rsv%time_in, w, i1,i2)
    call read_rsv_depth(i1)
  endif

  atots = 0.0 ; atotl = 0.0 !; rsv_depth = 0.
  do l = lnd%ls,lnd%le
     ce = first_elmt(land_tile_map(l))
     do while (loop_over_tiles(ce,tile))
        if (associated(tile%soil)) atots(l) = atots(l) + tile%frac !soil frac
        if (associated(tile%lake))then
          atotl(l) = atotl(l) + tile%frac !lake frac
          if(.not.use_reservoir) tile%lake%rsv_depth = -1.
        endif
     enddo
  enddo

  !do k1 = 2,2 ! 1 is lake, and 2 is soil
  !do k2 = 1,1
     ! get transition rate for this specific transition
     frac(:) = 0.0
     if (timel0==set_date(0001,01,01).and.infile_lake_state%ncid>0) then
        ! read initial transition from state file
        call input_state_lake(2,1)%interpolate(time,frac,interp='before')
     else
        if (input_tran_lake(2,1)%nvars>0) then
           call input_tran_lake(2,1)%integrate(timel0,time,frac)
        endif
     endif

     do l=lnd%ls,lnd%le
       if(lnd%ug_area(l)>0.)then
         frac(l) = frac(l)*(lnd%ug_cellarea(l)/lnd%ug_area(l)) !frac2land
       else
         frac(l) = 0.0
       endif
       !if(k1==2.and.k2==1)then !from soil to lake
         frac(l) = min(atots(l), frac(l))
         if(atotl(l)<=0.0) frac(l)=0.0 !Currently we do not build reservoir if there is no original lake in the gridcell.
       !else !from lake to soil
         !frac(l) = min(atotl(l), frac(l))
         !if(atots(l)==0.) frac(l)=0.
       !endif
     enddo
  !enddo
  !enddo

  do l = lnd%ls,lnd%le
     if(empty(land_tile_map(l))) cycle ! skip cells where there is no land  !irr code
     ! set current point for debugging
     call set_current_point(l,1)
     ! transition land area between different tile types
     call lake_transitions_0d(land_tile_map(l), &
          2, &
          1,&
          frac(l))
  enddo

  ! store current time for future reference
  timel0=time

end subroutine lake_transitions
! =============================================================================
! performs tile transitions in a given grid cell
subroutine lake_transitions_0d(d_list,d_kinds,a_kinds,area)
  type(land_tile_list_type), intent(inout) :: d_list ! list of tiles
  integer, intent(in) :: d_kinds ! array of donor tile kinds
  integer, intent(in) :: a_kinds ! array of acceptor tile kinds
  real   , intent(in) :: area    ! array of areas changing from donor tiles to acceptor tiles

  ! ---- local vars
  integer :: i, k
  type(land_tile_type), pointer :: ptr
  type(land_tile_list_type) :: a_list
  type(land_tile_enum_type) :: ts, te
  real :: atots, atotl, atot, frac0
  real :: htot ! total fraction heat, for debugging only
  ! variable used for conservation check:
  real :: lmass0, fmass0, cmass0, heat0, & !lmass fmass:kg/m2, cmass:kgC/m2, heat:J/m2, ! pre-transition values
       soil_heat0, vegn_heat0, cana_heat0, snow_heat0, lake_heat0, glac_heat0, e_res_heat0, &
       soil_lmass0, vegn_lmass0, cana_lmass0, snow_lmass0, &
       soil_fmass0, vegn_fmass0, cana_fmass0, snow_fmass0, &
       soil_cmass0, vegn_cmass0, cana_cmass0
  real :: cana_lmass0_soil, cana_fmass0_soil, cana_heat0_soil, cana_cmass0_soil, &
          snow_lmass0_soil, snow_fmass0_soil, snow_heat0_soil
  real :: cana_lmass0_lake, cana_fmass0_lake, cana_heat0_lake, cana_cmass0_lake, &
          snow_lmass0_lake, snow_fmass0_lake, snow_heat0_lake
  real :: lmass1, fmass1, cmass1, heat1, &
       soil_heat1, vegn_heat1, cana_heat1, snow_heat1, lake_heat1, glac_heat1, e_res_heat1 ! post-transition values
  real :: lwup0, lwup1
  real :: lm, fm, diff ! buffers for transition calculations
  real :: thres = 1.e-8
  type(cana_tile_type), pointer :: cana_soil
  real :: lwup_soil, e_res_1_soil, e_res_2_soil
  real :: frac_ac
  real :: fict_heat_dif = 0.


  if(area<=0.) return  !do nothing if no area change
  ! conservation check code, part 1: calculate the pre-transition grid
  ! cell totals
  lmass0 = 0 ; fmass0 = 0 ; cmass0 = 0 ; heat0 = 0
  soil_heat0  = 0 ;  vegn_heat0  = 0 ; cana_heat0  = 0 ; snow_heat0  = 0
  lake_heat0  = 0 ;  glac_heat0  = 0 ; e_res_heat0 = 0
  soil_lmass0 = 0 ;  vegn_lmass0 = 0 ; cana_lmass0 = 0 ; snow_lmass0 = 0
  soil_fmass0 = 0 ;  vegn_fmass0 = 0 ; cana_fmass0 = 0 ; snow_fmass0 = 0
  soil_cmass0 = 0 ;  vegn_cmass0 = 0 ; cana_cmass0 = 0

  cana_lmass0_soil=0; cana_fmass0_soil=0; cana_heat0_soil=0; cana_cmass0_soil=0
  snow_lmass0_soil=0; snow_fmass0_soil=0; snow_heat0_soil=0
  cana_lmass0_lake=0; cana_fmass0_lake=0; cana_heat0_lake=0; cana_cmass0_lake=0
  snow_lmass0_lake=0; snow_fmass0_lake=0; snow_heat0_lake=0

  lwup0 = 0

  ts = first_elmt(d_list)
  do while (loop_over_tiles(ts, ptr))
     call get_tile_water(ptr,lm,fm)
     lmass0 = lmass0 + lm*ptr%frac ; fmass0 = fmass0 + fm*ptr%frac
     if(associated(ptr%lake))then
       lmass0 = lmass0 + ptr%lake%sub_lmass*ptr%frac
       fmass0 = fmass0 + ptr%lake%sub_fmass*ptr%frac
     endif

     if (associated(ptr%soil))then
       call soil_tile_stock_pe(ptr%soil, lm, fm)
       soil_lmass0 = soil_lmass0 + lm*ptr%frac ; soil_fmass0 = soil_fmass0 + fm*ptr%frac
     endif
     if (associated(ptr%vegn))then
       call vegn_tile_stock_pe(ptr%vegn, lm, fm)
       vegn_lmass0 = vegn_lmass0 + lm*ptr%frac ; vegn_fmass0 = vegn_fmass0 + fm*ptr%frac
     endif
     if (associated(ptr%snow))then
       call snow_tile_stock_pe(ptr%snow, lm, fm)
       snow_lmass0 = snow_lmass0 + lm*ptr%frac ; snow_fmass0 = snow_fmass0 + fm*ptr%frac
     endif
     if (associated(ptr%cana))then
       call cana_tile_stock_pe(ptr%cana, lm, fm)
       cana_lmass0 = cana_lmass0 + lm*ptr%frac ; cana_fmass0 = cana_fmass0 + fm*ptr%frac
     endif

     if (associated(ptr%snow).and.associated(ptr%soil))then
       call snow_tile_stock_pe(ptr%snow, lm, fm)
       snow_lmass0_soil = snow_lmass0_soil + lm*ptr%frac ; snow_fmass0_soil = snow_fmass0_soil + fm*ptr%frac
     endif
     if (associated(ptr%cana).and.associated(ptr%soil))then
       call cana_tile_stock_pe(ptr%cana, lm, fm)
       cana_lmass0_soil = cana_lmass0_soil + lm*ptr%frac ; cana_fmass0_soil = cana_fmass0_soil + fm*ptr%frac
     endif

     if (associated(ptr%snow).and.associated(ptr%lake))then
       call snow_tile_stock_pe(ptr%snow, lm, fm)
       snow_lmass0_lake = snow_lmass0_lake + lm*ptr%frac ; snow_fmass0_lake = snow_fmass0_lake + fm*ptr%frac
     endif
     if (associated(ptr%cana).and.associated(ptr%lake))then
       call cana_tile_stock_pe(ptr%cana, lm, fm)
       cana_lmass0_lake = cana_lmass0_lake + lm*ptr%frac ; cana_fmass0_lake = cana_fmass0_lake + fm*ptr%frac
     endif


     heat0  = heat0  + land_tile_heat  (ptr)*ptr%frac
     if(associated(ptr%lake)) heat0 = heat0 + ptr%lake%sub_heat*ptr%frac

     if(associated(ptr%soil)) soil_heat0 = soil_heat0 + soil_tile_heat(ptr%soil)*ptr%frac
     if(associated(ptr%vegn)) vegn_heat0 = vegn_heat0 + vegn_tile_heat(ptr%vegn)*ptr%frac
     if(associated(ptr%cana)) cana_heat0 = cana_heat0 + cana_tile_heat(ptr%cana)*ptr%frac
     if(associated(ptr%snow)) snow_heat0 = snow_heat0 + snow_tile_heat(ptr%snow)*ptr%frac
     if(associated(ptr%lake)) lake_heat0 = lake_heat0 + lake_tile_heat(ptr%lake)*ptr%frac + ptr%lake%sub_heat*ptr%frac
     if(associated(ptr%glac)) glac_heat0 = glac_heat0 + glac_tile_heat(ptr%glac)*ptr%frac
     e_res_heat0 = e_res_heat0 + (ptr%e_res_1+ptr%e_res_2)*ptr%frac

     if(associated(ptr%cana).and.associated(ptr%soil)) cana_heat0_soil = cana_heat0_soil + cana_tile_heat(ptr%cana)*ptr%frac
     if(associated(ptr%snow).and.associated(ptr%soil)) snow_heat0_soil = snow_heat0_soil + snow_tile_heat(ptr%snow)*ptr%frac

     if(associated(ptr%cana).and.associated(ptr%lake)) cana_heat0_lake = cana_heat0_lake + cana_tile_heat(ptr%cana)*ptr%frac
     if(associated(ptr%snow).and.associated(ptr%lake)) snow_heat0_lake = snow_heat0_lake + snow_tile_heat(ptr%snow)*ptr%frac

     cmass0 = cmass0 + land_tile_carbon(ptr)*ptr%frac
     if(associated(ptr%lake)) cmass0 = cmass0 + ptr%lake%sub_cmass*ptr%frac

     if(associated(ptr%soil)) soil_cmass0 = soil_cmass0 + soil_tile_carbon(ptr%soil)*ptr%frac
     if(associated(ptr%vegn)) vegn_cmass0 = vegn_cmass0 + vegn_tile_carbon(ptr%vegn)*ptr%frac
     if(associated(ptr%cana)) cana_cmass0 = cana_cmass0 + cana_tile_carbon(ptr%cana)*ptr%frac

     if(associated(ptr%cana).and.associated(ptr%soil)) cana_cmass0_soil = cana_cmass0_soil + cana_tile_carbon(ptr%cana)*ptr%frac
     if(associated(ptr%cana).and.associated(ptr%lake)) cana_cmass0_lake = cana_cmass0_lake + cana_tile_carbon(ptr%cana)*ptr%frac

     lwup0 = lwup0 + ptr%lwup*ptr%frac
  enddo

  atots = 0.; atotl = 0.
  ts = first_elmt(d_list)
  do while (loop_over_tiles(ts,ptr))
    if (associated(ptr%soil)) atots = atots + ptr%frac
    if (associated(ptr%lake)) atotl = atotl + ptr%frac
  enddo

  if(d_kinds==2.and.a_kinds==1)then  !from soil to lake

    !before lake transition, we get the cana and other trivials on soil tile first.
    cana_soil => NULL()
    lwup_soil=0.; e_res_1_soil=0.; e_res_2_soil=0.
    frac_ac = 0.
    ts = first_elmt(d_list)
    do while (loop_over_tiles(ts,ptr))
      if(associated(ptr%soil))then !We assume that all soil tiles have cana here
        if(.not.associated(ptr%cana)) &
          call error_mesg('lake_transitions_0d','soil tile must have cana', FATAL)
        if(.not.associated(cana_soil))then
          cana_soil => new_cana_tile(ptr%cana)
          frac_ac = ptr%frac
        else
          call merge_cana_tiles(ptr%cana, ptr%frac, cana_soil, frac_ac)
          frac_ac = frac_ac + ptr%frac
        endif
        e_res_1_soil = e_res_1_soil + ptr%e_res_1*ptr%frac
        e_res_2_soil = e_res_2_soil + ptr%e_res_2*ptr%frac
        lwup_soil = lwup_soil + ptr%lwup*ptr%frac
      endif
    enddo

    !conduct the transition
    ts = first_elmt(d_list)
    do while (loop_over_tiles(ts,ptr))

      if(associated(ptr%lake))then
        if(use_reservoir.and.ptr%lake%rsv_depth<=0.) &
          call error_mesg('lake_transitions_0d','rsv_depth cannnot be less than 0 when using reservoir', FATAL)
        frac0 = ptr%frac
        ptr%frac = ptr%frac + area
        if(ptr%frac>1.+thres) &
          call error_mesg('lake_transitions_0d','lake fraction could not be greater than 1.', FATAL)
        if(ptr%lake%Afrac_rsv<1.) ptr%lake%Afrac_rsv = (ptr%lake%Afrac_rsv*frac0 + area)/ptr%frac

        ptr%lake%wl = ptr%lake%wl*frac0/ptr%frac
        ptr%lake%ws = ptr%lake%ws*frac0/ptr%frac
        ptr%lake%dz = ptr%lake%dz*frac0/ptr%frac
        ptr%lake%sub_lmass = ptr%lake%sub_lmass*frac0/ptr%frac
        ptr%lake%sub_fmass = ptr%lake%sub_fmass*frac0/ptr%frac
        ptr%lake%sub_heat  = ptr%lake%sub_heat *frac0/ptr%frac
        ptr%lake%sub_cmass = ptr%lake%sub_cmass*frac0/ptr%frac

        ptr%e_res_1 = ptr%e_res_1*frac0/ptr%frac
        ptr%e_res_2 = ptr%e_res_2*frac0/ptr%frac
        ptr%lwup = (ptr%lwup*frac0+lwup_soil*(area/atots))/ptr%frac

        if(associated(ptr%snow))then
          ptr%snow%wl = ptr%snow%wl*frac0/ptr%frac
          ptr%snow%ws = ptr%snow%ws*frac0/ptr%frac
          fict_heat_dif = snow_tile_heat(ptr%snow)*ptr%frac - snow_heat0_lake
          ptr%lake%sub_heat  = ptr%lake%sub_heat - fict_heat_dif/ptr%frac
        endif

        if(associated(ptr%cana).and.associated(cana_soil))then
          call merge_cana_tiles(cana_soil, area, ptr%cana, frac0)
          call delete_cana_tile(cana_soil)
        endif

        ptr%lake%sub_lmass = ptr%lake%sub_lmass + (soil_lmass0+vegn_lmass0+snow_lmass0_soil)*(area/atots)/ptr%frac !kg/m2
        ptr%lake%sub_fmass = ptr%lake%sub_fmass + (soil_fmass0+vegn_fmass0+snow_fmass0_soil)*(area/atots)/ptr%frac !kg/m2
        ptr%lake%sub_cmass = ptr%lake%sub_cmass + (soil_cmass0+vegn_cmass0)*(area/atots)/ptr%frac !kg/m2

        ptr%lake%sub_heat  = ptr%lake%sub_heat  + (soil_heat0+vegn_heat0+snow_heat0_soil )*(area/atots)/ptr%frac !J/m2
        ptr%lake%sub_heat  = ptr%lake%sub_heat  + (e_res_1_soil+e_res_2_soil)*(area/atots)/ptr%frac !J/m2

        call prohibit_shallow_lake(ptr%lake)
      endif

      if(associated(ptr%soil))then
        ptr%frac = ptr%frac*(1.-area/atots)
      endif

    enddo

    call land_tile_list_init(a_list)
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
      call merge_land_tile_into_list(ptr,d_list)
    enddo
    ! a_list is empty at this point
    call land_tile_list_end(a_list)


  else
    call error_mesg('lake_transitions_0d','transitions from lake to soil are not supported at this moment', FATAL)
  endif

  atot = 0.
  ts = first_elmt(d_list)
  do while (loop_over_tiles(ts,ptr))
     atot = atot + ptr%frac
  enddo
  if(abs(atot-1.)>thres) &
    call error_mesg('lake_transitions_0d','sum of tile frac must be zero', FATAL)

  ! conservation check part 2: calculate grid cell totals in final state, and
  ! compare them with pre-transition totals
  lmass1 = 0 ; fmass1 = 0 ; cmass1 = 0 ; heat1 = 0
  soil_heat1 = 0 ;  vegn_heat1 = 0 ; cana_heat1 = 0 ; snow_heat1 = 0
  lake_heat1 = 0 ;  glac_heat1 = 0 ; e_res_heat1 = 0
  lwup1 = 0
  ts = first_elmt(d_list)
  do while (loop_over_tiles(ts,ptr))
     call get_tile_water(ptr,lm,fm)
     lmass1 = lmass1 + lm*ptr%frac ; fmass1 = fmass1 + fm*ptr%frac
     if(associated(ptr%lake))then
       lmass1 = lmass1 + ptr%lake%sub_lmass*ptr%frac
       fmass1 = fmass1 + ptr%lake%sub_fmass*ptr%frac
     endif

     heat1  = heat1  + land_tile_heat  (ptr)*ptr%frac
     if(associated(ptr%lake)) heat1 = heat1 + ptr%lake%sub_heat*ptr%frac

     cmass1 = cmass1 + land_tile_carbon(ptr)*ptr%frac
     if(associated(ptr%lake)) cmass1 = cmass1 + ptr%lake%sub_cmass*ptr%frac

     if(associated(ptr%soil)) soil_heat1 = soil_heat1 + soil_tile_heat(ptr%soil)*ptr%frac
     if(associated(ptr%vegn)) vegn_heat1 = vegn_heat1 + vegn_tile_heat(ptr%vegn)*ptr%frac
     if(associated(ptr%cana)) cana_heat1 = cana_heat1 + cana_tile_heat(ptr%cana)*ptr%frac
     if(associated(ptr%snow)) snow_heat1 = snow_heat1 + snow_tile_heat(ptr%snow)*ptr%frac
     if(associated(ptr%lake)) lake_heat1 = lake_heat1 + lake_tile_heat(ptr%lake)*ptr%frac + ptr%lake%sub_heat*ptr%frac
     if(associated(ptr%glac)) glac_heat1 = glac_heat1 + glac_tile_heat(ptr%glac)*ptr%frac
     e_res_heat1 = e_res_heat1 + (ptr%e_res_1+ptr%e_res_2)*ptr%frac

     lwup1 = lwup1 + ptr%lwup*ptr%frac
  enddo
  call check_conservation ('liquid water', lmass0, lmass1, 1e-6)
  call check_conservation ('frozen water', fmass0, fmass1, 1e-6)
  call check_conservation ('carbon'      , cmass0, cmass1, 1e-6)
  call check_conservation ('canopy air heat content', cana_heat0, cana_heat1, 1e-4)
! heat content of vegetation may not conserve because of the cohort merging issues
  call check_conservation ('vegetation heat content', vegn_heat0*(1.-area/atots), vegn_heat1 , 1e-4)
  call check_conservation ('snow heat content',       snow_heat0-snow_heat0_soil*area/atots+fict_heat_dif, snow_heat1, 1e-4)
  call check_conservation ('soil heat content',       soil_heat0*(1.-area/atots), soil_heat1 , 1e-4)
  call check_conservation ('lake heat content',       lake_heat0+(soil_heat0+vegn_heat0+snow_heat0_soil+e_res_1_soil+e_res_2_soil)*area/atots-fict_heat_dif, lake_heat1, 1e-4)
  call check_conservation ('glac heat content',       glac_heat0, glac_heat1, 1e-4)
  call check_conservation ('e_res heat content',      e_res_heat0-(e_res_1_soil+e_res_2_soil)*area/atots, e_res_heat1, 1e-4)
  call check_conservation ('heat content', heat0 , heat1 , 1e-4)
  call check_conservation ('upward longwave', lwup0, lwup1, 1e-4)

end subroutine lake_transitions_0d
! =============================================================================
function time_adjust(time, time_list)
  type(time_type) :: time_adjust
  type(time_type), intent(in) :: time, time_list(:)

  integer :: n

  n = size(time_list)
  time_adjust = time
  if (time_adjust<time_list(1)) time_adjust = time_list(1)
  if (time_adjust>time_list(n)) time_adjust = time_list(n)

end function time_adjust
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
     call land_error_message(&
          'land use transition from "'//trim(landuse_name(d_kind))// &
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
        if(temp%vegn%landuse==LU_NTRL.or.temp%vegn%landuse==LU_SCND.or.temp%vegn%landuse==LU_RANGE) &
                call vegn_cut_forest(temp, a_kind)
        ! change landuse type of the tile
        temp%vegn%landuse = a_kind
        ! reset time elapsed since last disturbance and time elapsed since last land use
        ! event in the new tile
        temp%vegn%age_since_disturbance = 0.0
        temp%vegn%age_since_landuse     = 0.0
        ! add the new tile to the resulting list
        call insert(temp, a_list) ! insert tile into output list
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
  else if ((src==LU_RAINF.or.src==LU_IRRIG).and.dst==LU_PAST) then
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
        ! convert land use type of the tile:
        ! cut the forest, if necessary
        if(temp%vegn%landuse==LU_NTRL.or.temp%vegn%landuse==LU_SCND.or.temp%vegn%landuse==LU_RANGE) &
             call vegn_cut_forest(temp, a_kind)
        ! change landuse type of the tile
        temp%vegn%landuse = a_kind
        ! reset time elapsed since last disturbance and time elapsed since last land use
        ! event in the new tile
        temp%vegn%age_since_disturbance = 0.0
        temp%vegn%age_since_landuse     = 0.0
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

!=================================================================
subroutine add_irrigation_transitions(area0,tranI,fi1,atot,tran1,verbose)
  real, intent(in)  :: area0(M_LU_TYPES) ! area of each of the land use types [frac of land area]
  real, intent(in)  :: tranI(N_LU_TYPES, N_LU_TYPES)   ! initial transition matrix [frac of vegetated area per year]
  real, intent(in)  :: fi1 ! fraction of irrigated area after transition
  real, intent(in)  :: atot ! total area of vegetated tiles (i.e. the area that can be involved in transitions)
  real, intent(out) :: tran1(M_LU_TYPES, M_LU_TYPES) ! resulting transition matrix [frac of vegetated area per year]
  logical, intent(in), optional :: verbose

  ! local constants
  real, parameter :: tol = 1e-14 ! minimum value of non-zero transitions:
       ! for 1x1 degree grid, area is roughly 1e10 m2, so the area involved in transitions
       ! below tol would be below 1 cm2, which is probably safe to ignore
  integer, parameter :: IR=1, II=2, IZ=3 ! indices of aggregated LU types (rain-fed, irrigated, other)
  character(1), parameter :: tname(3) = ['r','i','z']

  type map1_t
     integer :: i, j;
     character(16) :: name
  end type map1_t

  ! local vars
  real    :: tran0(N_LU_TYPES, N_LU_TYPES) ! input transitions [frac of soil area per year]
  integer :: n  ! number of variables
  integer :: m1 ! Number of <= inequalities
  integer :: m3 ! Number of == equalities
  integer :: eq ! equation number
  real    :: area0c, area0z ! total area of crops and everything else before transitions
  real    :: area1i, area1r ! area of irrigated and rain-fed crops after transitions
  real    :: area1c ! total area of crops after transition
  real    :: c2z, z2c ! total transitions from and to crops, respectively
  integer :: map2(3,3) ! mapping from aggregated transition indices to variable index
  type(map1_t) :: map1(6) ! mapping var number to transitions
  real, allocatable :: c(:) ! coefficients of cost function
  real, allocatable :: A_ub(:,:), b_ub(:) ! input matrix and RHS for inequality (<=) conditions
  real, allocatable :: A_eq(:,:), b_eq(:) ! input matrix and RHS for equality conditions
  real, allocatable :: x(:) ! solution of linear optimization problem
  real    :: area00(N_LU_TYPES) ! fraction of each land use type before transition,
                                ! with irrigated and rain-fed crop added together
  integer :: i,j,ierr
  logical :: verbose_
  real    :: fi0 ! fraction of irrigated area before transition, for diagnostics only
  real    :: s   ! accumulator value for various calculations


  verbose_ = .FALSE.
  if (present(verbose)) verbose_ = verbose

  ! check that the area involved in transitions is greater then zero,
  ! and if it is not, return copy of input transitions
  if (atot<=0.0) then
     tran1(:,:) = 0.0
     do i = 1,N_LU_TYPES
     do j = 1,N_LU_TYPES
        tran1(i,j) = tranI(i,j)
     enddo
     enddo
     return
  endif

  ! calculate the fraction each land use type with rain-fed and irrigated areas combined
  area00(:) = area0(1:N_LU_TYPES)
  area00(LU_RAINF) = area00(LU_RAINF)+area0(LU_IRRIG)

  tran0 = tranI*atot ! convert transitions to [frac of land area per year]: it is easier to
                     ! work in this units since tile area units are [fractions of land area].

  ! filter out negatives in input transition matrix
  do i = 1,N_LU_TYPES
  do j = 1,N_LU_TYPES
     if (tran0(i,j) < tol) tran0(i,j)= 0.0
  enddo
  enddo

  ! Ensure that the total transitions from every land use type do not exceed area of that
  ! type of land use (no overshoots). This is necessary because linear programming
  ! optimization algorithms are quite sensitive to the overshoots: they lead to empty
  ! feasible solution region, and consequent failure of the optimization.
  do i = 1,N_LU_TYPES
     s = sum(tran0(i,:))
     if (s>area00(i)) then
        tran0(i,:) = tran0(i,:)*area00(i)/s
     endif
  enddo

  if (verbose_) then
     write(*,*)'INPUT DATA'
     write(*,*)'initial land use fractions'
     do i = 1,M_LU_TYPES
        write(*,'(a," : ",g12.4)') landuse_name(i),area0(i)
     enddo
     write(*,*)
     write(*,*)'atot:', atot

     write(*,*)
     write(*,*)'initial transition matrix:'
     call print_transitions(tran0)

     write(*,*)
     do i = 1,N_LU_TYPES
        write(*,'(99(a,g10.3))') 'sum of transitions from '//landuse_name(i)//':', sum(tran0(i,:)), &
                   ' to '//landuse_name(i)//':', sum(tran0(:,i))
     enddo

     if (area00(LU_RAINF) > 0) then
        fi0 = area0(LU_IRRIG)/area00(LU_RAINF)
     else
        fi0 = 0.0
     endif
     write(*,*)
     write(*,'(x,a,99(2x,a,g12.4:))') 'irrigated cropland fraction','before :', &
               fi0,'after :',fi1

  endif

  ! calculate combined transitions (assuming crop->crop transitions are zero)
  c2z = 0.0; z2c = 0.0
  do i = 1,N_LU_TYPES
     if (i == LU_RAINF) continue
     z2c = z2c + tran0(i,LU_RAINF)
     c2z = c2z + tran0(LU_RAINF,i)
  enddo

  ! calculate areas after transition
  area0c = area0(LU_RAINF) + area0(LU_IRRIG)
  area0z = sum(area0) - area0c

  area1c = area0c + z2c - c2z
  area1i = area1c*fi1
  area1r = area1c - area1i

  if (verbose_) then
     write(*,*)
     write(*,*)'total cropland area after transitions     :',area1c
     write(*,*)'irrigated cropland area after transitions :',area1i
     write(*,*)'rain-fed cropland area after transitions  :',area1r
  endif

  ! initialize mapping of aggregated transitions to variables. This is a constant array;
  ! the set-up can be moved to module initialization code
  n = 0
  do i = 1,3
  do j = 1,3
     if (i.ne.j) then
        n = n+1
        map2(i,j)  = n
        map1(n)%i = i
        map1(n)%j = j
        map1(n)%name = tname(i)//'2'//tname(j)
!         write(*,*) n, i, j, map1(n)%name
     else
        map2(i,j) = 0
     endif
  enddo
  enddo

  ! n  = 6 ! number of variables (z2i,z2r,i2z,r2z,i2r,r2i)
  m1 = 2 ! number of <= constraints: available area constraints
  m3 = 4 ! number of == constraints

  if (verbose_) then
     write(*,*)
     write(*,'(99(a,I2))') 'Number of variables      (n) :',n
     write(*,'(99(a,I2))') 'Number of <= constraints (m1):',m1
     write(*,'(99(a,I2))') 'Number of == constraints (m3):',m3
  endif

  ! The cost function
  allocate (c(n)) ; c(:) = 0.0
!  c(map2(IZ,II)) = 1;  c(map2(II,IZ)) = 1
  c(map2(IZ,II)) = 1.01;  c(map2(II,IZ)) = 1.01
  c(map2(IR,II)) = 1;     c(map2(II,IR)) = 1
  call print_equation(c,label='cost function:')

  ! set up <= conditions
  ! available area constraints
  allocate (A_ub(m1,n),b_ub(m1)) ; A_ub(:,:) = 0.0 ; b_ub(:) = 0.0
  A_ub(1,map2(II,IZ)) = 1 ; A_ub(1,map2(II,IR)) = 1 ; b_ub(1) = area0(LU_IRRIG)
  A_ub(2,map2(IR,IZ)) = 1 ; A_ub(2,map2(IR,II)) = 1 ; b_ub(2) = area0(LU_RAINF)
  call print_equation(A_ub(1,:), '<=', b_ub(1))
  call print_equation(A_ub(2,:), '<=', b_ub(2))

  ! set up == conditions
  allocate (A_eq(m3,n),b_eq(m3)) ; A_eq(:,:) = 0.0 ; b_eq(:) = 0.0
  eq = 1 ! sum of transitons must be equal to total irrigated area change
  b_eq(eq) = area1i - area0(LU_IRRIG)
  A_eq(eq,map2(IZ,II)) = 1 ; A_eq(eq,map2(II,IZ)) = -1
  A_eq(eq,map2(IR,II)) = 1 ; A_eq(eq,map2(II,IR)) = -1
  call print_equation(A_eq(eq,:), '==', b_eq(eq))

  eq = 2 ! sum of transitons must be equal to total rain-fed area change
  b_eq(eq) = area1r - area0(LU_RAINF)
  A_eq(eq,map2(IZ,IR)) = 1 ; A_eq(eq,map2(IR,IZ)) = -1
  A_eq(eq,map2(II,IR)) = 1 ; A_eq(eq,map2(IR,II)) = -1
  call print_equation(A_eq(eq,:), '==', b_eq(eq))

  ! equations (e6): sum of transitions to irrigated and rain-fed cropland is equal
  ! to the total transition to cropland
  eq=3
  A_eq(eq,map2(IZ,II)) = 1 ;  A_eq(eq,map2(IZ,IR)) = 1 ; b_eq(eq) = z2c
  call print_equation(A_eq(eq,:), '==', b_eq(eq))
  eq=4
  A_eq(eq,map2(II,IZ)) = 1 ;  A_eq(eq,map2(IR,IZ)) = 1 ; b_eq(eq) = c2z
  call print_equation(A_eq(eq,:), '==', b_eq(eq))

  if (verbose_) then
     write(*,*)' Input Tables for linprog:'
     do i = 1,n
        write(*,'(5x,a3)',advance='NO') map1(i)%name
     enddo
     write(*,*)
     do i = 1, size(A_ub,1)
        do j = 1,n
!            if (a(i,j).ne.0) then
               write(*,'(f8.4)',advance='NO')A_ub(i,j)
!            else
!                write(*,'(8x)',advance='NO')
!            endif
        enddo
        write(*,'(a,f8.4)') ' <=', b_ub(i)
     enddo
     do i = 1, size(A_eq,1)
        do j = 1,n
!            if (a(i,j).ne.0) then
               write(*,'(f8.4)',advance='NO')A_eq(i,j)
!            else
!                write(*,'(8x)',advance='NO')
!            endif
        enddo
        write(*,'(a,f8.4)') ' ==', b_eq(i)
     enddo
  endif

  allocate (x(n)) ; x(:) = 0.0
  call linprog(c,A_ub,b_ub,A_eq,b_eq,x,ierr)
  if (ierr == 0) then
     if (verbose_) write(*,*) 'simplx finished successfully'
  else
     call land_error_message('add_irrigation_transitions: simplx failed with code '//string(-ierr),FATAL)
     ! if (verbose_) write(*,*) 'simplx finished un-successfully, with code',icase
  endif

  if (verbose_) then
     write(*,*)
!      write(*,*) ' Maximum of objective function = ', A(1,1)
     do i=1,n
        write(*,'("  x",i2.2," : ", a," = ",g12.5)') i, trim(map1(i)%name), x(i)
     enddo
  endif

  ! unpack the solution into final transition matrix
  tran1(:,:) = 0.0
  tran1(1:N_LU_TYPES,1:N_LU_TYPES) = tran0(:,:)
  if (c2z>0) then
     do i = 1,N_LU_TYPES
        tran1(LU_IRRIG,i) = x(map2(II,IZ))*tran0(LU_RAINF,i)/c2z
        tran1(LU_RAINF,i) = x(map2(IR,IZ))*tran0(LU_RAINF,i)/c2z
     enddo
  endif
  if (z2c>0) then
     do i = 1,N_LU_TYPES
        tran1(i,LU_IRRIG) = x(map2(IZ,II))*tran0(i,LU_RAINF)/z2c
        tran1(i,LU_RAINF) = x(map2(IZ,IR))*tran0(i,LU_RAINF)/z2c
     enddo
  endif
  tran1(LU_IRRIG,LU_RAINF) = x(map2(II,IR))
  tran1(LU_RAINF,LU_IRRIG) = x(map2(IR,II))

  if (verbose_) then
    write(*,*)
    write(*,*)'final transition matrix:'
    call print_transitions(tran1)

    write(*,*)
    write(*,'(a/5x,2a23)')'change in land use fractions'
    do i = 1,M_LU_TYPES
       write(*,'(a5,g12.4)',advance='NO') landuse_name(i),area0(i)
       s = area0(i)
       do j = 1,M_LU_TYPES
          s = s + tran1(j,i) - tran1(i,j)
       enddo
       write(*,'("->",g12.4)') s
    enddo

    do i = 1,M_LU_TYPES
       s = sum(tran1(i,:))
       if (s > area0(i)) then
          write (*,'("sum of transitions from ",a," (",g10.3,") exceeds initial area (",g10.3,") by", g23.16)') &
              landuse_name(i),s,area0(i), s-area0(i)
       endif
    enddo
  endif

  ! convert transition units back to [fraction of vegetated area per year]
  tran1 = tran1/atot

  deallocate (A_ub, b_ub, A_eq, b_eq, x)

contains

  subroutine print_equation(a,op,b,label)
     real,         intent(in)           :: a(:)
     character(*), intent(in), optional :: op
     real,         intent(in), optional :: b
     character(*), intent(in), optional :: label

     integer :: i
     character :: sign

     if (.NOT.verbose_) return

     if (present(label)) write(*,'(a)',advance='NO') label
     do i = 1,size(a)
        if (a(i)==0) cycle
        sign = '+'
        if (a(i) < 0) sign = '-'
        write(*,'(x,a,x)',advance='NO') sign
        if (abs(a(i)).ne.1.0)  write(*,'(f8.4,x)',advance='NO') abs(a(i))
        write(*,'(a3)',advance='NO')map1(i)%name
     enddo
     if (present(op)) write(*,'(x,a)',advance='NO') op
     if (present(b))  write(*,'(f8.4)',advance='NO') b
     write(*,*)
  end subroutine print_equation

  subroutine print_transitions(tran)
     real, intent(in) :: tran(:,:)
     integer :: i, j

     do i = 1, size(tran,1)
     do j = 1, size(tran,2)
         if (tran(i,j) <= 0) cycle
         write(*,'(a4,"->",a4,g12.3)') landuse_name(i), landuse_name(j), tran(i,j)
     enddo
     enddo
  end subroutine print_transitions

end subroutine add_irrigation_transitions

! optimize linear programming problem
subroutine linprog(c,A_ub,b_ub,A_eq,b_eq,x,ierr)
  real, intent(in)  :: c(:)      ! coefficients of linear functions to minimize
  real, intent(in)  :: A_ub(:,:) ! the inequality constrain matrix
  real, intent(in)  :: b_ub(:)   ! the inequality constrain vector: A_ub*x <= b_ub
  real, intent(in)  :: A_eq(:,:) ! the equality constrain matrix
  real, intent(in)  :: b_eq(:)   ! the equality constrain vector: A_eq*x == b_eq

  real, intent(out) :: x(:)      ! solution vector
  integer, intent(out) :: ierr

  integer :: n ! number of variables
  integer :: m ! total number of constraints
  integer :: m1 ! number of <= constraints
  integer :: m2 ! number of >= constraints
  integer :: m3 ! number of == constraints
  real, allocatable :: tableau(:,:) ! "tableau" for the simplex module
  integer, allocatable :: izrov(:), iposv(:)
  integer :: i,j,k,eq

  n = size(c)
  ! sanity checks
  if (size(A_ub,2).ne.n)          call error_mesg('linprog','size of A_ub is inconsistent with number of variables',FATAL)
  if (size(A_ub,1).ne.size(b_ub)) call error_mesg('linprog','sizes of A_ub and b_ub are inconsistent',FATAL)
  if (size(A_eq,2).ne.n)          call error_mesg('linprog','size of A_eq is inconsistent with number of variables',FATAL)
  if (size(A_eq,1).ne.size(b_eq)) call error_mesg('linprog','sizes of A_eq and b_eq are inconsistent',FATAL)
  if (size(c).ne.size(x))         call error_mesg('linprog','sizes of c and x are inconsistent',FATAL)

  m = size(A_ub,1) + size(A_eq,1)

  m1 = count(b_ub(:)>=0) ! count <= conditions
  m2 = size(b_ub) - m1   ! and >= conditions
  m3 = size(b_eq)        ! and the number of == conditions

  ! size of the tableau is m+2 because one line is reserved for cost function, and another
  ! for internal needs of simplex method implementation
  allocate(tableau(m+2,n+1))
  tableau(:,:) = 0.0
  ! function to optimize: we are changing the sign because simplx maximizes, not minimizes cost
  tableau(1,2:) = -c(:)
  eq = 2
  do k = 1,size(b_ub)
     if (b_ub(k)>=0) then
        tableau(eq,1)  =  b_ub(k)
        tableau(eq,2:) = -A_ub(k,:)
        eq = eq+1
     endif
  enddo
  do k = 1,size(b_ub)
     if (b_ub(k)<0) then
        tableau(eq, 1)  = -b_ub(k)
        tableau(eq, 2:) =  A_ub(k,:)
        eq = eq+1
     endif
  enddo
  do k = 1,size(b_eq)
     if (b_eq(k)>=0) then
        tableau(eq, 1)  =  b_eq(k)
        tableau(eq, 2:) = -A_eq(k,:)
     else
        tableau(eq, 1)  = -b_eq(k)
        tableau(eq, 2:) =  A_eq(k,:)
     endif
     eq = eq+1
  enddo

!   write(*,*)'m1,m2,m3:', m1, m2, m3
!   write(*,*)'Tableau:'
!   do k = 1,size(tableau,1)
!      write(*,'(999(f8.4))') tableau(k,:)
!   enddo

  allocate(izrov(n),iposv(m))
  call simplx(tableau,m,n,m1,m2,m3,ierr,izrov,iposv)
  do i=1, n
     do j=1, m
        if (IPOSV(J).eq.I) then
           x(i) = tableau(j+1,1)
           goto 3
        endif
     enddo
3 enddo

  deallocate(izrov,iposv,tableau)
end subroutine linprog

!****************************************************************
!*           LINEAR PROGRAMMING: THE SIMPLEX METHOD
!* ------------------------------------------------------------
!* SAMPLE RUN:
!* Maximize z = x1 + x2 + 3x3 -0.5x4 with conditions:
!*          x1  + 2x3 <= 740
!*          2x2 - 7x4 <= 0
!*          x2  - x3 + 2x4 >= 0.5
!*          x1 + x2 + x3 +x4 = 9
!*          and all variables xN  >=0.
!*
!* Number of variables in E.F.: 4
!* Number of <= inequalities..: 2
!* Number of >= inequalities..: 1
!* Number of = equalities.....: 1
!*
!*  Input Table:
!*    0.00    1.00    1.00    3.00   -0.50
!*  740.00   -1.00    0.00   -2.00    0.00
!*    0.00    0.00   -2.00    0.00    7.00
!*    0.50    0.00   -1.00    1.00   -2.00
!*    9.00   -1.00   -1.00   -1.00   -1.00
!*
!*  Maximum of E.F. =    17.02500
!*  X 1 =    0.000000
!*  X 2 =    3.325000
!*  X 3 =    4.725000
!*  X 4 =    0.950000
!*
!* ------------------------------------------------------------
!* Reference: "Numerical Recipes By W.H. Press, B. P. Flannery,
!*             S.A. Teukolsky and W.T. Vetterling, Cambridge
!*             University Press, 1986"
!****************************************************************
!----------------------------------------------------------------------------------------
! USES simp1,simp2,simp3
!Simplex method for linear programming. Input parameters a, m, n, mp, np, m1, m2, and m3,
!and output parameters a, icase, izrov, and iposv are described above.
!number of variables expected; EPS is the absolute precision, which should be adjusted to
!the scale of your variables.
subroutine simplx(a,m,n,m1,m2,m3,icase,izrov,iposv)
  real,    intent(inout) :: a(:,:)
  integer, intent(in) :: n  ! number of variables
  integer, intent(in) :: m  ! total number of constraints (must be equal to m1+m2+m3)
  integer, intent(in) :: m1 ! number of <= constraints
  integer, intent(in) :: m2 ! number of >= constraints
  integer, intent(in) :: m3 ! number of == constraints
  integer, intent(out) :: icase
  integer, intent(out) :: izrov(:)
  integer, intent(out) :: iposv(:)

  real, parameter :: EPS = 1e-6

  integer :: m12, nl2, ir
  INTEGER :: i,ip,is,k,kh,kp,nl1,l1(n),l2(m),l3(m)
  REAL bmax,q1
  character(len=512) :: mesg ! for error message


  ! TODO: check size of input tableau a
  if (m.ne.m1+m2+m3) call error_mesg('simplx','Bad input constraint counts in simplx.', FATAL)
  if (size(izrov)<n) call error_mesg('simplx','ERROR: size of izrov is less then number of variables (N)', FATAL)
  if (size(iposv)<m) call error_mesg('simplx','ERROR: size of ipozv is less then number of constraints (M)', FATAL)

  nl1=n
  do k=1,n
    l1(k)=k    !Initialize index list of columns admissible for exchange.
    izrov(k)=k !Initially make all variables right-hand.
  end do
  nl2=m
  do i=1,m
    if(a(i+1,1).lt.0.) then
        call check_var_range(a(i+1,1),0.0,HUGE(1.0),'Bad input tableau in simplx, Constants bi must be nonnegative.','a(i+1,1)', FATAL)
        write(mesg,*) 'Bad input tableau in simplx, Constants bi must be nonnegative.', a(i+1,1)
        call error_mesg('simplx',mesg, FATAL)!pause ' Bad input tableau in simplx, Constants bi must be nonnegative.'
    endif
    l2(i)=i
    iposv(i)=n+i
  !-------------------------------------------------------------------------------------------------
  !Initial left-hand variables. m1 type constraints are represented by having their slackv ariable
  !initially left-hand, with no artificial variable. m2 type constraints have their slack
  !variable initially left-hand, with a minus sign, and their artificial variable handled implicitly
  !during their first exchange. m3 type constraints have their artificial variable initially
  !left-hand.
  !-------------------------------------------------------------------------------------------------
  end do
  do i=1,m2
    l3(i)=1
  end do
  ir=0
  if(m2+m3.eq.0) goto 30 !The origin is a feasible starting solution. Go to phase two.
  ir=1
  do k=1,n+1             !Compute the auxiliary objective function.
    q1=0.
    do i=m1+1,m
      q1=q1+a(i+1,k)
    end do
    a(m+2,k)=-q1
  end do
  10 call simp1(a,m+1,l1,nl1,0,kp,bmax) !Find max. coeff. of auxiliary objective fn
  if(bmax.le.EPS.and.a(m+2,1).lt.-EPS)then
    !write(*,*) 'bmax', bmax, 'a', a(m+2,1)
    icase=-1        !Auxiliary objective function is still negative and cannot be improved,
    return          !hence no feasible solution exists.
  else if(bmax.le.EPS.and.a(m+2,1).le.EPS)then
  !Auxiliary objective function is zero and cannot be improved; we have a feasible starting vector.
  !Clean out the artificial variables corresponding to any remaining equality constraints by
  !goto 1 and then move on to phase two by goto 30.
    m12=m1+m2+1
    if (m12.le.m) then
      do ip=m12,m
        if(iposv(ip).eq.ip+n)then !Found an artificial variable for an equalityconstraint.
          call simp1(a,ip,l1,nl1,1,kp,bmax)
          if(bmax.gt.EPS) goto 1  !Exchange with column corresponding to maximum
        end if                    !pivot element in row.
      end do
    end if
    ir=0
    m12=m12-1
    if (m1+1.gt.m12) goto 30
    do i=m1+1,m1+m2               !Change sign of row for any m2 constraints
                                  !still present from the initial basis.
      if(l3(i-m1).eq.1)then
        do k=1,n+1
          a(i+1,k)=-a(i+1,k)
        end do
      end if
    end do
    goto 30                        !Go to phase two.
  end if
  call simp2(a,m,n,l2,nl2,ip,kp,q1) !Locate a pivot element (phase one).
  if(ip.eq.0)then                  !Maximum of auxiliary objective function is
                                   !unbounded, so no feasible solution exists.
    icase=-2
    return
  end if
  1 call simp3(a,m+1,n,ip,kp)
  !Exchange a left- and a right-hand variable (phase one), then update lists.
  if(iposv(ip).ge.n+m1+m2+1)then   !Exchanged out an artificial variable for an
                                   !equality constraint. Make sure it stays
                                   !out by removing it from the l1 list.
    do k=1,nl1
      if(l1(k).eq.kp) goto 2
    end do
  2 nl1=nl1-1
    do is=k,nl1
      l1(is)=l1(is+1)
    end do
  else
    if(iposv(ip).lt.n+m1+1) goto 20
    kh=iposv(ip)-m1-n
    if(l3(kh).eq.0) goto 20      !Exchanged out an m2 type constraint.
    l3(kh)=0                     !If it is the first time, correct the pivot column
                                 !or the minus sign and the implicit
                                 !artificial variable.
  end if
  a(m+2,kp+1)=a(m+2,kp+1)+1.
  do i=1,m+2
    a(i,kp+1)=-a(i,kp+1)
  end do
  20 is=izrov(kp)                !Update lists of left- and right-hand variables.
  izrov(kp)=iposv(ip)
  iposv(ip)=is
  if (ir.ne.0) goto 10           !if still in phase one, go back to 10.
  !End of phase one code for finding an initial feasible solution. Now, in phase two, optimize it.
  30 call simp1(a,0,l1,nl1,0,kp,bmax) !Test the z-row for doneness.
  if(bmax.le.EPS)then            !Done. Solution found. Return with the good news.
    icase=0
    return
  end if
  call simp2(a,m,n,l2,nl2,ip,kp,q1)  !Locate a pivot element (phase two).
  if(ip.eq.0)then                !Objective function is unbounded. Report and return.
    icase=1
    return
  end if
  call simp3(a,m,n,ip,kp)  !Exchange a left- and a right-hand variable (phase two),
  goto 20                        !update lists of left- and right-hand variables and
                                 !return for another iteration.
end subroutine simplx

!The preceding routine makes use of the following utility subroutines:

! ==============================================================================
! Determines the maximum of those elements whose index is contained in the
! supplied list ll, either with or without taking the absolute value, as flagged
! by iabf.
SUBROUTINE simp1(a,mm,ll,nll,iabf,kp,bmax)
  real,    intent(in)  :: a(:,:) ! tableau
  integer, intent(in)  :: mm,ll(:),nll,iabf
  integer, intent(out) :: kp
  real,    intent(out) :: bmax

  integer :: k
  real    :: test

  kp=ll(1)
  bmax=a(mm+1,kp+1)
  if(nll.lt.2) return
  do k=2,nll
    if(iabf.eq.0)then
      test=a(mm+1,ll(k)+1)-bmax
    else
      test=abs(a(mm+1,ll(k)+1))-abs(bmax)
    endif
    if(test.gt.0.)then
      bmax=a(mm+1,ll(k)+1)
      kp=ll(k)
    endif
  end do
end subroutine simp1

! ==============================================================================
! Locate a pivot element, taking degeneracy into account.
subroutine simp2(a,m,n,l2,nl2,ip,kp,q1)
  real,    intent(in)  :: a(:,:)
  integer, intent(in)  :: m,n,l2(:),nl2, kp
  integer, intent(out) :: ip
  real,    intent(out) :: q1

  real, parameter :: EPS=1.e-6
  integer :: i,k,ii
  real    :: q,q0,qp
  ip=0
  if(nl2.lt.1) return
  do i=1,nl2
    if(a(i+1,kp+1).lt.-EPS) goto 2
  end do
  return  ! No possible pivots. Return with message.
2 q1=-a(l2(i)+1,1)/a(l2(i)+1,kp+1)
  ip=l2(i)
  if(i+1.gt.nl2) return
  do i=i+1, nl2
    ii=l2(i)
    if(a(ii+1,kp+1).lt.-EPS)then
      q=-a(ii+1,1)/a(ii+1,kp+1)
      if(q.lt.q1)then
        ip=ii
        q1=q
        else if (q.eq.q1) then !We have a degeneracy.
        do k=1,n
          qp=-a(ip+1,k+1)/a(ip+1,kp+1)
          q0=-a(ii+1,k+1)/a(ii+1,kp+1)
          if(q0.ne.qp)goto 6
        end do
6       if(q0.lt.qp) ip=ii
      end if
    end if
  end do
end subroutine simp2

! ==============================================================================
! Matrix operations to exchange a left-hand and right-hand variable (see text).
subroutine simp3(a,i1,k1,ip,kp)
  real, intent(inout) :: a(:,:)
  integer, intent(in) :: i1,k1,ip,kp

  integer :: ii,kk
  real    :: piv

  piv=1./a(ip+1,kp+1)
  if (i1.ge.0) then
    do ii=1,i1+1
      if(ii-1.ne.ip)then
        a(ii,kp+1)=a(ii,kp+1)*piv
        do kk=1,k1+1
          if(kk-1.ne.kp)then
            a(ii,kk)=a(ii,kk)-a(ip+1,kk)*a(ii,kp+1)
          end if
        end do
      end if
    end do
  end if
  do kk=1,k1+1
    if(kk-1.ne.kp) a(ip+1,kk)=-a(ip+1,kk)*piv
  end do
  a(ip+1,kp+1)=piv
end subroutine

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
