module nitrogen_sources_mod

#include "shared/debug.inc"

use constants_mod, only : PI
use time_manager_mod, only : time_type, get_date, operator(/=), operator(-), &
     operator(<), valid_calendar_types, get_calendar_type, time_type_to_real
use fms_mod, only : &
     file_exist, open_namelist_file, &
     check_nml_error, close_file, stdlog, mpp_pe, mpp_root_pe, error_mesg, &
     FATAL, NOTE, string
use mpp_domains_mod, only : mpp_pass_SG_to_UG
use time_interp_external_mod, only : init_external_field, time_interp_external, &
     time_interp_external_init
use horiz_interp_mod, only : horiz_interp_type, horiz_interp_init, &
     horiz_interp_new, horiz_interp_del, horiz_interp
use time_interp_mod, only : time_interp
use get_cal_time_mod, only : get_cal_time

use nfu_mod, only : nfu_validtype, nfu_inq_var, nfu_get_dim_bounds, nfu_get_rec, &
     nfu_get_dim, nfu_get_valid_range, nfu_is_valid

use land_constants_mod, only : seconds_per_year
use land_data_mod, only: land_state_type, lnd, log_version
use land_io_mod, only : read_field
use land_tile_io_mod, only : print_netcdf_error
use land_tile_diag_mod, only : cmor_name, &
     register_tiled_diag_field, send_tile_data, diag_buff_type, set_default_diag_filter, &
     add_tiled_diag_field_alias
use land_debug_mod, only : is_watch_point
use vegn_data_mod, only : LU_PAST, LU_CROP
use soil_carbon_mod, only : soil_carbon_option, SOILC_CORPSE_N

implicit none
private

! ==== public interfaces =====================================================
public :: nitrogen_sources_init
public :: nitrogen_sources_end
public :: update_nitrogen_sources
public :: nitrogen_sources
public :: do_nitrogen_deposition
! ==== end of public interfaces ===============================================

! ==== module constants =======================================================
character(len=*), parameter :: module_name = 'nitrogen_sources'
#include "shared/version_variable.inc"
character(len=*), private, parameter :: &
   diag_mod_name = 'vegn'

integer, parameter :: & ! types of deposition calculations
   NDEP_PRESCRIBED = 1, & ! single prescribed value
   NDEP_INTERP2    = 2, & ! interpolated between natural and anthropogenic
   NDEP_TIMESERIES = 3    ! regular time series

! ==== NetCDF declarations ===================================================
include 'netcdf.inc'
#define __NF_ASRT__(x) call print_netcdf_error((x),module_name,__LINE__)

! ==== module variables ======================================================

! ---- namelist --------------------------------------------------------------
real    :: wetdepfrac = 0.0  ! fraction of deposition associated with rainfall
real    :: ndeprate   = 0.0  ! prescribed deposition rate
character(len=32) :: n_deposition_type = 'interpolate-two-maps' ! or 'prescribed',
                             ! or 'time-series'
integer :: year_nat   = 1860 ! year before which natural deposition rate is used
integer :: year_ant   = 1990 ! year after each anthropogenic deposition rate is used

real    :: ndep_nit_frac   = 1.0 ! fraction of nitrate in deposition
real    :: ndep_amm_frac   = 0.0 ! fraction of ammonium in deposition

real    :: fert_nit_frac   = 1.0 ! fraction of nitrate in fertilizer
real    :: fert_amm_frac   = 0.0 ! fraction of ammonium in fertilizer

real    :: manure_nit_frac = 1.0 ! fraction of nitrate in manure
real    :: manure_amm_frac = 0.0 ! fraction of ammonium in manure

real    :: ndep_amm_volat  = 0.0 ! fraction of ammonium volatilization for deposition
real    :: fert_amm_volat  = 0.0 ! fraction of ammonium volatilization for fertilizer application
real    :: manure_amm_volat= 0.0 ! fraction of ammonium volatilization for manure application

logical :: do_nitrogen_deposition = .FALSE. ! Whether to do N dep at all

character(256) :: &
   ndep_nit_file = 'INPUT/ndep_nit.nc', & ! input file name for nitrate deposition
   ndep_amm_file = 'INPUT/ndep_amm.nc', &              ! input file name for ammonium deposition
   ndep_org_file = 'INPUT/ndep_org.nc'                 ! input file name for DON deposition
character(32) :: &
   ndep_nit_field = 'ndep_nit', &         ! input field name for nitrate deposition
   ndep_amm_field = 'ndep_amm', &             ! input field name for ammonium deposition
   ndep_org_field = 'ndep_org'                ! input field name for ammonium deposition
character(256) :: &
   fml_nit_crop_file   = '', &
   fml_amm_crop_file   = '', &
   fml_org_crop_file   = '', &
   fml_nit_past_file   = '', &
   fml_amm_past_file   = '', &
   fml_org_past_file   = '', &
   fert_crop_file   = '', &           ! input file for cropland fertilization
   fert_past_file   = '', &           ! input file for pasture fertilization
   manure_crop_file = '', &           ! input file for cropland manure application
   manure_past_file = ''              ! input file for pasture manure application
namelist/nitrogen_deposition_nml/n_deposition_type, &
   ndep_nit_file, ndep_nit_field, &
   ndep_amm_file, ndep_amm_field, &
   ndep_org_file, ndep_org_field, &
   wetdepfrac, year_nat,year_ant,ndeprate, &
   ndep_nit_frac,   ndep_amm_frac,   &
   fml_nit_crop_file, fml_amm_crop_file, fml_org_crop_file, &
   fml_nit_past_file, fml_amm_past_file, fml_org_past_file, &
   fert_crop_file,    fert_past_file,    fert_nit_frac,   fert_amm_frac,   &
   manure_crop_file,  manure_past_file,  manure_nit_frac, manure_amm_frac, &
   ndep_amm_volat, fert_amm_volat, manure_amm_volat, do_nitrogen_deposition

! ---- end of namelist -------------------------------------------------------

real :: ndep_org_frac, fert_org_frac, manure_org_frac ! fractions of organic
        ! nitrogen in deposition, fertilizer, and manure, respectively
integer :: n_deposition_option = 0 !selector for deposition options
real, allocatable ::   & ! buffers for input data
   data_SG(:,:),       & ! structured grid input buffer
   ndep_nat(:),        & ! natural deposition rate
   ndep_ant(:),        & ! anthropogenic deposition rate
   ndep_nit_buffer(:), & ! nitrate deposition
   ndep_amm_buffer(:), & ! ammonium deposition
   ndep_org_buffer(:), & ! DON source
   nfml_nit_crop(:), nfml_amm_crop(:), nfml_org_crop(:),  &
   nfml_nit_past(:), nfml_amm_past(:), nfml_org_past(:),  &
   nfert_crop(:),    nfert_past(:), &  ! fertilization for crops and pastures
   nmanure_crop(:),  nmanure_past(:)   ! manure nitrogen for crops and pastures

integer :: &
   id_ndep_nit_in = -1, & ! id of external field for time interpolation
   id_ndep_amm_in = -1, & ! id of external field for time interpolation
   id_ndep_org_in = -1    ! id of external field for time interpolation

type :: fert_data_type ! type describing input fertilization data
   integer :: ncid = -1  ! netcdf ID of the input data file
   integer :: nlon_in = -1, nlat_in = -1 ! input grid size
   type(time_type), pointer :: time_in(:) => NULL() ! input time axis
   real, pointer :: lonb_in(:,:) => NULL() ! input longitude boundaries
   real, pointer :: latb_in(:,:) => NULL() ! input latitude boundaries

   integer :: curr_rec = -1 ! number of the current record
   real, pointer :: nfert(:) => NULL() ! amount of fertilizer
   real, pointer :: tfert(:) => NULL() ! time of application
   real, pointer :: mask (:) => NULL() ! mask of valid data
   type(nfu_validtype) :: nfert_v ! valid range for amount of fertilizer
   type(nfu_validtype) :: tfert_v ! valid range for time of fertilizer
end type

type(fert_data_type) :: &
   fml_nit_crop, fml_amm_crop, fml_org_crop, &
   fml_nit_past, fml_amm_past, fml_org_past, &
   fert_crop, manure_crop, & ! crop fertilization
   fert_past, manure_past    ! pasture fertilization

integer :: & ! diag field IDs
   id_ndep_nit, id_ndep_amm, id_ndep_org, id_ndep, &
   id_nfert_nit, id_nfert_amm, id_nfert_org, id_nfert, &
   id_nmanure_nit, id_nmanure_amm, id_nmanure_org, id_nmanure, &
   id_nfml_nit, id_nfml_amm, id_nfml_org, id_nfml, &
   id_amm_volat

contains ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

! ============================================================================
subroutine nitrogen_sources_init(time, id_ug)
  type(time_type), intent(in) :: time  ! current time
  integer,         intent(in) :: id_ug ! ids of diagnostic axes

  ! ---- local vars
  integer :: unit, ierr, io

  call log_version(version, module_name, &
  __FILE__)
  if (file_exist('input.nml')) then
     unit = open_namelist_file()
     ierr = 1;
     do while (ierr /= 0)
        read (unit, nml=nitrogen_deposition_nml, iostat=io, end=10)
        ierr = check_nml_error (io, 'nitrogen_deposition_nml')
     enddo
10   continue
     call close_file (unit)
  endif
  if (mpp_pe() == mpp_root_pe()) then
     unit = stdlog()
     write (unit, nml=nitrogen_deposition_nml)
     call close_file (unit)
  endif

  if (trim(n_deposition_type)=='interpolate-two-maps') then
    n_deposition_option = NDEP_INTERP2
  else if (trim(n_deposition_type)=='prescribed') then
    n_deposition_option = NDEP_PRESCRIBED
  else if (trim(n_deposition_type)=='time-series') then
    n_deposition_option = NDEP_TIMESERIES
  else
    call error_mesg('nitrogen_sources_init','n_deposition_type = "'//trim(n_deposition_type)&
          //'" is invalid, use "prescribed" or "interpolate-two-maps"',&
          FATAL)
  endif

  ! do nothing further if no nitrogen in the model
  if (do_nitrogen_deposition) then
     if(soil_carbon_option .NE. SOILC_CORPSE_N) call error_mesg('nitrogen_sources_init',&
         'WARNING: do_nitrogen_deposition is TRUE but soil_carbon_option is not SOILC_CORPSE_N. N dep will likely have no effect.',NOTE)

     ! calculate fractions of organic nitrogen in deposition, fertilizer, and  manure
     ndep_org_frac   = 1.0 - ndep_nit_frac   - ndep_amm_frac
     fert_org_frac   = 1.0 - fert_nit_frac   - fert_amm_frac
     manure_org_frac = 1.0 - manure_nit_frac - manure_amm_frac

     ! allocate and initialize buffer for data
     allocate(data_SG(lnd%is:lnd%ie,lnd%js:lnd%je))
     allocate(ndep_nit_buffer(lnd%ls:lnd%le))
     allocate(ndep_amm_buffer(lnd%ls:lnd%le))
     allocate(ndep_org_buffer(lnd%ls:lnd%le))
     allocate(nfml_nit_crop  (lnd%ls:lnd%le))
     allocate(nfml_amm_crop  (lnd%ls:lnd%le))
     allocate(nfml_org_crop  (lnd%ls:lnd%le))
     allocate(nfml_nit_past  (lnd%ls:lnd%le))
     allocate(nfml_amm_past  (lnd%ls:lnd%le))
     allocate(nfml_org_past  (lnd%ls:lnd%le))
     allocate(nfert_crop     (lnd%ls:lnd%le))
     allocate(nfert_past     (lnd%ls:lnd%le))
     allocate(nmanure_crop   (lnd%ls:lnd%le))
     allocate(nmanure_past   (lnd%ls:lnd%le))

     ! read the data, if necessary
     select case (n_deposition_option)
     case (NDEP_INTERP2)
       allocate( ndep_nat(lnd%ls:lnd%le), &
                 ndep_ant(lnd%ls:lnd%le)  )
!                  print *,'lon',lnd%lon(lnd%is:lnd%ie,lnd%js:lnd%je)
!                  print *,'lat',lnd%lat(lnd%is:lnd%ie,lnd%js:lnd%je)
       call read_field( 'INPUT/nbiodata_nat.nc','ndep', ndep_nat, interp='nearest')
       call read_field( 'INPUT/nbiodata_ant.nc','ndep', ndep_ant, interp='nearest')
     case (NDEP_TIMESERIES)
       ! initialize external fields
       call time_interp_external_init()
       if (trim(ndep_nit_file) /= '') then
          ! initilaize external field
          id_ndep_nit_in = init_external_field(trim(ndep_nit_file),trim(ndep_nit_field),&
                                         domain=lnd%sg_domain,use_comp_domain=.TRUE.,verbose=.FALSE.)
       endif
       if (trim(ndep_amm_file) /= '') then
          ! initilaize external field
          id_ndep_amm_in = init_external_field(trim(ndep_amm_file),trim(ndep_amm_field),&
                                            domain=lnd%sg_domain,use_comp_domain=.TRUE.)
       endif
       if (trim(ndep_org_file) /= '') then
          ! initilaize external field
          id_ndep_org_in = init_external_field(trim(ndep_org_file),trim(ndep_org_field),&
                                            domain=lnd%sg_domain,use_comp_domain=.TRUE.)
       endif
     case default
       ! do nothing
     end select

     call fert_data_init(  fml_nit_crop_file,   fml_nit_crop)
     call fert_data_init(  fml_amm_crop_file,   fml_amm_crop)
     call fert_data_init(  fml_org_crop_file,   fml_org_crop)
     call fert_data_init(  fml_nit_past_file,   fml_nit_past)
     call fert_data_init(  fml_amm_past_file,   fml_amm_past)
     call fert_data_init(  fml_org_past_file,   fml_org_past)
     call fert_data_init(  fert_crop_file,   fert_crop)
     call fert_data_init(manure_crop_file, manure_crop)
     call fert_data_init(  fert_past_file,   fert_past)
     call fert_data_init(manure_past_file, manure_past)

  else ! If do_nitrogen_deposition is not .TRUE.
     call error_mesg('nitrogen_sources_init',&
       'do_nitrogen_deposition=.FALSE.: NOT initializing the nitrogen deposition sources', &
       NOTE)
  endif

  ! ---- initialize the diagnostics --------------------------------------------
  call set_default_diag_filter('soil')

  id_ndep_nit = register_tiled_diag_field (diag_mod_name, 'Ndep_nit', (/id_ug/), &
       time, 'nitrate deposition', 'kg N/(m2 year)', missing_value=-999.0)
  id_ndep_amm = register_tiled_diag_field (diag_mod_name, 'Ndep_amm', (/id_ug/), &
       time, 'ammonium deposition', 'kg N/(m2 year)', missing_value=-999.0)
  id_ndep_org = register_tiled_diag_field (diag_mod_name, 'Ndep_org', (/id_ug/), &
       time, 'organic nitrogen deposition', 'kg N/(m2 year)', missing_value=-999.0)
  id_ndep    = register_tiled_diag_field (diag_mod_name, 'Ndep', (/id_ug/), &
       time, 'total nitrogen deposition', 'kg N/(m2 year)', missing_value=-999.0)

  id_nfert_nit = register_tiled_diag_field (diag_mod_name, 'Nfert_nit', (/id_ug/), &
       time, 'nitrate fertilization', 'kg N/(m2 year)', missing_value=-999.0)
  id_nfert_amm = register_tiled_diag_field (diag_mod_name, 'Nfert_amm', (/id_ug/), &
       time, 'ammonium fertilization', 'kg N/(m2 year)', missing_value=-999.0)
  id_nfert_org = register_tiled_diag_field (diag_mod_name, 'Nfert_org', (/id_ug/), &
       time, 'organic nitrogen fertilization', 'kg N/(m2 year)', missing_value=-999.0)
  id_nfert     = register_tiled_diag_field (diag_mod_name, 'Nfert', (/id_ug/), &
       time, 'total nitrogen fertilization', 'kg N/(m2 year)', missing_value=-999.0)

  id_nmanure_nit = register_tiled_diag_field (diag_mod_name, 'Nmanure_nit', (/id_ug/), &
       time, 'manure nitrate', 'kg N/(m2 year)', missing_value=-999.0)
  id_nmanure_amm = register_tiled_diag_field (diag_mod_name, 'Nmanure_amm', (/id_ug/), &
       time, 'manure ammonium', 'kg N/(m2 year)', missing_value=-999.0)
  id_nmanure_org = register_tiled_diag_field (diag_mod_name, 'Nmanure_org', (/id_ug/), &
       time, 'manure organic nitrogen', 'kg N/(m2 year)', missing_value=-999.0)
  id_nmanure    = register_tiled_diag_field (diag_mod_name, 'Nmanure', (/id_ug/), &
       time, 'total manure nitrogen deposition', 'kg N/(m2 year)', missing_value=-999.0)

  id_nfml_nit = register_tiled_diag_field (diag_mod_name, 'Nfml_nit', (/id_ug/), &
       time, 'nitrate fertilization', 'kg N/(m2 year)', missing_value=-999.0)
  id_nfml_amm = register_tiled_diag_field (diag_mod_name, 'Nfml_amm', (/id_ug/), &
       time, 'ammonium fertilization', 'kg N/(m2 year)', missing_value=-999.0)
  id_nfml_org = register_tiled_diag_field (diag_mod_name, 'Nfml_org', (/id_ug/), &
       time, 'organic nitrogen fertilization', 'kg N/(m2 year)', missing_value=-999.0)
  id_nfml     = register_tiled_diag_field (diag_mod_name, 'Nfml', (/id_ug/), &
       time, 'total nitrogen fertilization', 'kg N/(m2 year)', missing_value=-999.0)

  id_amm_volat = register_tiled_diag_field (diag_mod_name, 'amm_volat', (/id_ug/), &
       time, 'rate of ammonium input volatilization', 'kg N/(m2 year)', missing_value=-999.0)

  ! CMOR/CMIP variables
  call set_default_diag_filter('land')

  call add_tiled_diag_field_alias(id_nfert, cmor_name, 'fNfert', [id_ug], time, &
      'total N added for cropland fertilisation (artificial and manure)', 'kg m-2 s-1', missing_value=-999.0, &
      standard_name='tendency_of_soil_mass_content_of_nitrogen_compounds_expressed_as_nitrogen_due_to_fertilization', &
      fill_missing=.TRUE.)
  call add_tiled_diag_field_alias(id_ndep, cmor_name, 'fNdep', [id_ug], time, &
      'Dry and Wet Deposition of Reactive Nitrogen onto Land', 'kg m-2 s-1', missing_value=-999.0, &
      standard_name='tendency_of_atmosphere_mass_content_of_nitrogen_compounds_expressed_as_nitrogen_due_to_deposition', &
      fill_missing=.TRUE.)

end subroutine nitrogen_sources_init


! ============================================================================
#define __DEALLOC__(x) if (allocated(x)) deallocate(x)
subroutine nitrogen_sources_end()
  ! clean up deposition fields
  __DEALLOC__(ndep_nat)
  __DEALLOC__(ndep_ant)
  __DEALLOC__(ndep_nit_buffer); __DEALLOC__(ndep_amm_buffer); __DEALLOC__(ndep_org_buffer)
  __DEALLOC__(nfml_nit_crop); __DEALLOC__(nfml_amm_crop); __DEALLOC__(nfml_org_crop);
  __DEALLOC__(nfert_crop); __DEALLOC__(nmanure_crop)
  __DEALLOC__(nfml_nit_past); __DEALLOC__(nfml_amm_past); __DEALLOC__(nfml_org_past);
  __DEALLOC__(nfert_past); __DEALLOC__(nmanure_past)
  call fert_data_end(  fml_nit_crop)
  call fert_data_end(  fml_amm_crop)
  call fert_data_end(  fml_org_crop)
  call fert_data_end(  fml_nit_past)
  call fert_data_end(  fml_amm_past)
  call fert_data_end(  fml_org_past)
  call fert_data_end(  fert_crop)
  call fert_data_end(manure_crop)
  call fert_data_end(  fert_past)
  call fert_data_end(manure_past)
end subroutine
#undef __DEALLOC__

! =============================================================================
subroutine update_nitrogen_sources (t0,t1)
  type(time_type)     , intent(in)    :: t0,t1 ! time interval

  ! ---- local vars
  integer :: year,month,day,hour,minute,second
  real :: w

  if (.not.do_nitrogen_deposition) return ! do nothing in this case

  ! update nitrogen deposition
  select case(n_deposition_option)
  case (NDEP_PRESCRIBED)
     ndep_nit_buffer(:) = ndeprate*ndep_nit_frac
     ndep_amm_buffer(:) = ndeprate*ndep_amm_frac
     ndep_org_buffer(:) = ndeprate*ndep_org_frac
  case (NDEP_INTERP2)
     call get_date(t0, year,month,day,hour,minute,second)
     if (year_ant > year_nat) then
       w = real(year - year_nat)/real(year_ant-year_nat)
       w = min(w,1.0); w = max(w,0.0)
     else ! They are the same year, so use the nat one -BNS
       w = 0.0
     endif
     ndep_nit_buffer(:) = (1-w)*ndep_nat(:)+w*ndep_ant(:)
     ndep_amm_buffer(:) = ndep_nit_buffer(:)*ndep_amm_frac
     ndep_org_buffer(:) = ndep_nit_buffer(:)*ndep_org_frac
     ndep_nit_buffer(:) = ndep_nit_buffer(:)*ndep_nit_frac
  case (NDEP_TIMESERIES)
     if (id_ndep_nit_in>0) then
        call time_interp_external(id_ndep_nit_in,t0,data_SG)
        call mpp_pass_SG_to_UG(lnd%ug_domain, data_SG, ndep_nit_buffer)
     endif
     if (id_ndep_amm_in>0) then
        call time_interp_external(id_ndep_amm_in,t0,data_SG)
        call mpp_pass_SG_to_UG(lnd%ug_domain, data_SG, ndep_amm_buffer)
     endif
     if (id_ndep_org_in>0) then
        call time_interp_external(id_ndep_org_in,t0,data_SG)
        call mpp_pass_SG_to_UG(lnd%ug_domain, data_SG, ndep_org_buffer)
     endif
  case default
    call error_mesg('nitrogen_deposition','n_deposition_option '//&
          string(n_deposition_option)//' is invalid"', FATAL)
  end select

  call fert_data_get( fml_nit_crop, t0,t1, nfml_nit_crop )
  call fert_data_get( fml_amm_crop, t0,t1, nfml_amm_crop )
  call fert_data_get( fml_org_crop, t0,t1, nfml_org_crop )
  call fert_data_get( fml_nit_past, t0,t1, nfml_nit_past )
  call fert_data_get( fml_amm_past, t0,t1, nfml_amm_past )
  call fert_data_get( fml_org_past, t0,t1, nfml_org_past )
  call fert_data_get(    fert_crop, t0,t1, nfert_crop )
  call fert_data_get(  manure_crop, t0,t1, nmanure_crop )
  call fert_data_get(    fert_past, t0,t1, nfert_past )
  call fert_data_get(  manure_past, t0,t1, nmanure_past )
end subroutine update_nitrogen_sources


! ============================================================================
subroutine nitrogen_sources(time, l, p_ann, precip, lu, input_nit, input_amm, input_org, diag)
  type(time_type), intent(in)  :: time ! current time, for time interpolation
  integer, intent(in)  :: l       ! index of current point in unstructured grid
  real   , intent(in)  :: p_ann   ! annual-mean precipitation, kg/(m2 s)
  real   , intent(in)  :: precip  ! current precipitation, kg/(m2 s)
  integer, intent(in)  :: lu      ! land use type
  real   , intent(out) :: input_nit ! total nitrate input, kg N/(m2 yr)
  real   , intent(out) :: input_amm ! total ammonium input, kg N/(m2 yr)
  real   , intent(out) :: input_org ! total organic nitrogen input, kg N/(m2 yr)
  type(diag_buff_type), intent(inout) :: diag

  ! ---- local vars
  real :: ndep_nit, ndep_amm, ndep_org ! nitrogen deposition rates, kg N/(m2 yr)
  real :: fert_nit, fert_amm, fert_org ! nitrogen fertilization rates, kg N/(m2 yr)
  real :: fml_nit, fml_amm, fml_org ! nitrogen fertilization rates, kg N/(m2 yr)
  real :: manure_nit, manure_amm, manure_org ! nitrogen manure application rates, kg N/(m2 yr)
  real :: p, w

  if (.not.do_nitrogen_deposition) then
     input_nit=0.0 ; input_amm=0.0 ; input_org = 0.0
     return
  endif

  ! scale the deposition by precipitation
  if (p_ann>0) then
     p = precip/p_ann ; w = wetdepfrac
  else
     p = 0.0          ; w = 0.0
  endif
  ndep_nit = ndep_nit_buffer(l)*(1-w) + w*ndep_nit_buffer(l)*p
  ndep_amm = ndep_amm_buffer(l)*(1-w) + w*ndep_amm_buffer(l)*p
  ndep_org = ndep_org_buffer(l)*(1-w) + w*ndep_org_buffer(l)*p

  ! add fertilization
  select case(lu)
  case (LU_CROP)
     fert_nit   = nfert_crop(l)   * fert_nit_frac
     fert_amm   = nfert_crop(l)   * fert_amm_frac
     fert_org   = nfert_crop(l)   * fert_org_frac
     manure_nit = nmanure_crop(l) * manure_nit_frac
     manure_amm = nmanure_crop(l) * manure_amm_frac
     manure_org = nmanure_crop(l) * manure_org_frac
     fml_nit    = nfml_nit_crop(l)
     fml_amm    = nfml_amm_crop(l)
     fml_org    = nfml_org_crop(l)
  case (LU_PAST)
     fert_nit   = nfert_past(l)   * fert_nit_frac
     fert_amm   = nfert_past(l)   * fert_amm_frac
     fert_org   = nfert_past(l)   * fert_org_frac
     manure_nit = nmanure_past(l) * manure_nit_frac
     manure_amm = nmanure_past(l) * manure_amm_frac
     manure_org = nmanure_past(l) * manure_org_frac
     fml_nit    = nfml_nit_past(l)
     fml_amm    = nfml_amm_past(l)
     fml_org    = nfml_org_past(l)
  case default
       fert_nit = 0.0 ;   fert_amm = 0.0 ;   fert_org = 0.0
     manure_nit = 0.0 ; manure_amm = 0.0 ; manure_org = 0.0
        fml_nit = 0.0 ;    fml_amm = 0.0 ;    fml_org = 0.0
  end select

  input_nit = ndep_nit + fert_nit + manure_nit + fml_nit
  input_amm = ndep_amm*(1-ndep_amm_volat) + fert_amm*(1-fert_amm_volat) + manure_amm*(1-manure_amm_volat) + fml_amm
  input_org = ndep_org + fert_org + manure_org + fml_org

  if (is_watch_point()) then
     call dpri('i=',lnd%i_index(l)); call dpri('j=',lnd%j_index(l))
     __DEBUG3__(ndep_nit_buffer(l),ndep_amm_buffer(l),ndep_org_buffer(l))
     __DEBUG2__(w,p)
     write(*,*) '#### nitrogen_sources output'
     __DEBUG5__(ndep_nit, fert_nit, manure_nit, input_nit, fml_nit)
     __DEBUG5__(ndep_amm, fert_amm, manure_amm, input_amm, fml_amm)
     __DEBUG5__(ndep_org, fert_org, manure_org, input_org, fml_org)
  endif


  ! ---- diagnostic section ----------------------------------------------------
  call send_tile_data(id_ndep_nit, ndep_nit, diag)
  call send_tile_data(id_ndep_amm, ndep_amm, diag)
  call send_tile_data(id_ndep_org, ndep_org, diag)
  call send_tile_data(id_ndep,     ndep_nit+ndep_amm+ndep_org, diag)

  call send_tile_data(id_nfert_nit, fert_nit, diag)
  call send_tile_data(id_nfert_amm, fert_amm, diag)
  call send_tile_data(id_nfert_org, fert_org, diag)
  call send_tile_data(id_nfert,     fert_nit+fert_amm+fert_org, diag)

  call send_tile_data(id_nfml_nit, fml_nit, diag)
  call send_tile_data(id_nfml_amm, fml_amm, diag)
  call send_tile_data(id_nfml_org, fml_org, diag)
  call send_tile_data(id_nfml,     fml_nit+fml_amm+fml_org, diag)

  call send_tile_data(id_nmanure_nit, manure_nit, diag)
  call send_tile_data(id_nmanure_amm, manure_amm, diag)
  call send_tile_data(id_nmanure_org, manure_org, diag)
  call send_tile_data(id_nmanure,     manure_nit+manure_amm+manure_org, diag)

  call send_tile_data(id_amm_volat,   &
    ndep_amm*ndep_amm_volat + fert_amm*fert_amm_volat + manure_amm*manure_amm_volat, &
    diag)
end subroutine nitrogen_sources


! ============================================================================
subroutine fert_data_init(filename, fert)
  character(*), intent(in) :: filename
  type(fert_data_type), intent(inout) :: fert

  ! ---- local vars
  integer :: dimids(NF_MAX_VAR_DIMS) ! IDs of dimensions
  integer :: dimlens(NF_MAX_VAR_DIMS) ! sizes of dimensions
  integer :: timedim ! id of the record (time) dimension
  integer :: timevar ! id of the time variable
  character(len=NF_MAX_NAME) :: timename  ! name of the time variable
  character(len=256)         :: timeunits ! units of time in the file
  character(len=24) :: calendar ! model calendar
  real, allocatable :: time(:)  ! time line from fertilization file
  integer :: i
  integer :: ierr ! error code
  integer :: nrec ! number of records in the fertilization data

  if (trim(filename)=='') return ! do nothing if the file name is empty string

  ierr=nf_open(filename,NF_NOWRITE,fert%ncid)

  if(ierr/=NF_NOERR) call error_mesg('nitrogen_sources_init', &
       'do_fertilization is requested, but fertilization "'// &
       trim(filename)//'" could not be opened because '//nf_strerror(ierr), &
       FATAL)

  __NF_ASRT__(nfu_inq_var(fert%ncid,'nfert',dimids=dimids,dimlens=dimlens,nrec=nrec))
  allocate(time(nrec), fert%time_in(nrec))

  ! read fertilization time line
  ! get the time axis
  __NF_ASRT__(nf_inq_unlimdim(fert%ncid, timedim))
  __NF_ASRT__(nfu_get_dim(fert%ncid, timedim, time))
  ! get units of time
  __NF_ASRT__(nf_inq_dimname(fert%ncid, timedim, timename))
  __NF_ASRT__(nf_inq_varid(fert%ncid,timename, timevar))
  timeunits = ' '
  __NF_ASRT__(nf_get_att_text(fert%ncid,timevar,'units',timeunits))
  ! get model calendar
  calendar=valid_calendar_types(get_calendar_type())

  ! loop through the time axis and store time_type values in time_in
  do i = 1,size(time)
     fert%time_in(i) = get_cal_time(time(i),timeunits,calendar)
  end do

  ! get the boundaries of the horizontal axes
  allocate(fert%lonb_in(dimlens(1)+1,1), &
           fert%latb_in(1,dimlens(2)+1)  )
  __NF_ASRT__(nfu_get_dim_bounds(fert%ncid, dimids(1), fert%lonb_in(:,1)))
  __NF_ASRT__(nfu_get_dim_bounds(fert%ncid, dimids(2), fert%latb_in(1,:)))
  ! convert them to radian
  fert%lonb_in = fert%lonb_in*PI/180.0
  fert%latb_in = fert%latb_in*PI/180.0
  ! set up the input buffer sizes
  fert%nlon_in = dimlens(1)
  fert%nlat_in = dimlens(2)
  ! set up the valid ranges
  __NF_ASRT__(nfu_get_valid_range(fert%ncid,'nfert',fert%nfert_v))
  ierr = nfu_get_valid_range(fert%ncid,'tfert',fert%tfert_v) ! ignore errors here

  ! allocate buffers for data
  allocate(fert%nfert(lnd%ls:lnd%le))
  allocate(fert%tfert(lnd%ls:lnd%le))
  ! and mask
  allocate(fert%mask (lnd%ls:lnd%le))
  ! set the current record indicator to non-existing record number
  fert%curr_rec = -1

  ! deallocate temporary variables
  deallocate(time)

end subroutine


! ============================================================================
#define __DEALLOC__(x) if (associated(x)) deallocate(x)
subroutine fert_data_end(fert)
  type(fert_data_type), intent(inout) :: fert

  integer :: ierr ! error code

  __DEALLOC__(fert%time_in)
  __DEALLOC__(fert%lonb_in)
  __DEALLOC__(fert%latb_in)
  __DEALLOC__(fert%nfert)
  __DEALLOC__(fert%tfert)
  __DEALLOC__(fert%mask)
  fert%curr_rec=-1

  ierr = nf_close(fert%ncid) ! ignore errors here
  fert%ncid = -1
end subroutine
#undef __DEALLOC__

! =============================================================================
! given interval [t0,t1), returns for each grid point fertilization that is to
! be applied during this interval, in kg N/(m2 year)
subroutine fert_data_get(fert, t0,t1,fert_rate)
  type(fert_data_type), intent(inout) :: fert  ! fertilization data
  type(time_type)     , intent(in)    :: t0,t1 ! time interval
  real                , intent(out)   :: fert_rate(:)

  ! ---- local vars
  integer :: i1,i2 ! record indices
  integer :: ierr ! error code
  integer :: i,j
  real    :: w
  real    :: t0r, t1r ! real values of t0,t1, days since the beginning of
                      ! current time interval in the input data
  real, allocatable :: tfert_in(:,:) ! input buffer for tfert
  real, allocatable :: nfert_in(:,:) ! input buffer for nfert
  real, allocatable :: mask_in(:,:)  ! mask for interpolation
  type(horiz_interp_type) :: interp

  ! set initial value
  fert_rate(:) = 0

  ! do nothing further if the input was not initialized (or initialized with
  ! empty filename)
  if (fert%ncid<0) return

  ! do nothing if we are before the start of the data: assume there is no
  ! fertilization
  if(t0<fert%time_in(1)) return

  ! get the time input data time interval for the current date
  call time_interp(t0, fert%time_in, w, i1, i2)
  ! make sure we have the proper record in the memory
  if (fert%curr_rec/=i1) then
     ! we need to read/interpolate  the record.
     allocate(tfert_in(fert%nlon_in,fert%nlat_in))
     allocate(nfert_in(fert%nlon_in,fert%nlat_in))
     allocate(mask_in (fert%nlon_in,fert%nlat_in))
     ! We have to construct interpolator every time, because the mask may change
     ! for every record

     ! read fertilization time
     ierr = nfu_get_rec(fert%ncid,"tfert",i1,tfert_in)
     if(ierr==NF_NOERR) then
        where (nfu_is_valid(tfert_in,fert%tfert_v))
          mask_in = 1
        elsewhere
          mask_in = 0
        end where
     else
        tfert_in = 0
        mask_in  = 1
     endif

     ! read fertilization amount
     __NF_ASRT__(nfu_get_rec(fert%ncid,"nfert",i1,nfert_in))
     ! mask is 1 only where both tfert and nfert are valid
     where (.not.nfu_is_valid(nfert_in,fert%nfert_v)) &
           mask_in = mask_in*0

     ! construct interpolator
     call horiz_interp_new(interp, fert%lonb_in, fert%latb_in, &
          lnd%sg_lonb(lnd%is:lnd%ie+1,lnd%js:lnd%je+1), &
          lnd%sg_latb(lnd%is:lnd%ie+1,lnd%js:lnd%je+1), &
          interp_method='conservative',&
          mask_in=mask_in, is_latlon_in=.TRUE. )
     ! note that mask_out does not seem to be set if both input and output grids
     ! are lat-lon (horiz_interp_new_2d, line 546). So it's pretty useless in
     ! this case.

     ! interpolate input fields
     call horiz_interp(interp,nfert_in,data_SG)
     call mpp_pass_SG_to_UG(lnd%ug_domain, data_SG, fert%nfert)
     call horiz_interp(interp,tfert_in,data_SG)
     call mpp_pass_SG_to_UG(lnd%ug_domain, data_SG, fert%tfert)
     call horiz_interp(interp,mask_in,data_SG)
     call mpp_pass_SG_to_UG(lnd%ug_domain, data_SG, fert%mask)
     ! convert time to seconds (from days)
     fert%tfert = fert%tfert*86400
     where(fert%mask==0)
        fert%tfert = -9999.0
        fert%nfert = 0.0
     end where

     ! update the current record number
     fert%curr_rec = i1

     ! clean up memory
     deallocate(tfert_in,nfert_in,mask_in)
     call horiz_interp_del(interp)
  endif

  ! convert time intervals to real numbers: note that in FMS time is positive
  ! by definition [why, oh why???] so that t0-t1 == t1-t0, and an extra step is
  ! necessary to get the right sign of the real time intervals...
  t0r = time_type_to_real(t0-fert%time_in(i1))
  if (t0<fert%time_in(i1)) t0r = -t0r
  t1r = time_type_to_real(t1-fert%time_in(i1))
  if (t1<fert%time_in(i1)) t1r = -t1r

  ! set up the fertilization values
  where(t0r<=fert%tfert.and.fert%tfert<t1r) &
      fert_rate = fert%nfert
  ! convert amount (kg N/m2) to rate (kg N /(m2 year)), same as deposition
  fert_rate = fert_rate/(t1r-t0r)*seconds_per_year
end subroutine

end module nitrogen_sources_mod
