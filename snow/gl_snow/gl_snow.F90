
! snow model module
! ============================================================================
module gl_snow_mod

#include "../../shared/debug.inc"

#ifdef INTERNAL_FILE_NML
use mpp_mod, only: input_nml_file
#else
use fms_mod, only: open_namelist_file
#endif

! EZSNOW - for 
! use mpp_io_mod, only : fieldtype, mpp_get_info, mpp_get_fields, mpp_get_axis_data, &
!      mpp_read, validtype, mpp_get_atts, MPP_RDONLY, MPP_NETCDF, MPP_MULTI, MPP_SINGLE, &
!      axistype, mpp_open, mpp_close, mpp_is_valid, mpp_get_file_name, mpp_get_field_index

use fms_mod, only : error_mesg, file_exist, check_nml_error, &
     stdlog, close_file, mpp_pe, mpp_root_pe, FATAL, WARNING, NOTE
use time_manager_mod,   only: time_type_to_real
use constants_mod,      only: tfreeze, hlv, hlf, PI

use land_constants_mod, only : NBANDS
! use snow_tile_mod, only : &
!      snow_tile_type, read_snow_data_namelist, &
!      snow_data_thermodynamics, snow_data_area, &
!      snow_active, &
!      snow_data_hydraulics, max_lev, cpw, clw, csw, use_brdf

! use parent_snow_tile_mod, only : &
!      snow_tile_type, read_snow_data_namelist, &
!      snow_data_thermodynamics, snow_data_area, &
!      snow_active, &
!      snow_data_hydraulics, max_lev, cpw, clw, csw, use_brdf

! use parent_snow_tile_mod, only : snow_tile_type, max_lev, cpw, clw, csw, use_brdf

use parent_snow_tile_mod, only : &
     ! snow_tile_type, &
     read_snow_data_namelist, &
     snow_data_thermodynamics, snow_data_area, &
   !   snow_active, &
     snow_data_hydraulics, max_lev, cpw, clw, csw, use_brdf

use gl_snow_tile_mod, only: gl_snow_tile_type


use snicar_mod, only: read_snicar_optics_data, read_snow_snicar_namelist

use land_tile_mod,    only : land_tile_map, land_tile_type, land_tile_list_type, &
     land_tile_enum_type, first_elmt, tail_elmt, next_elmt, &
     current_tile, operator(/=), nitems, loop_over_tiles

! use land_tile_mod, only : land_tile_map, land_tile_type, land_tile_enum_type, &
!      first_elmt, loop_over_tiles
use land_data_mod, only : lnd, log_version
use land_tile_io_mod, only: land_restart_type, &
     init_land_restart, open_land_restart, save_land_restart, free_land_restart, &
     add_restart_axis, add_tile_data, get_tile_data, get_tile_by_idx, &
     add_int_tile_data, get_int_tile_data, field_exists
use land_debug_mod, only : is_watch_point, check_var_range

! use nf_utils_mod,     only : nfu_inq_dim, nfu_get_var, nfu_put_var, &
!      nfu_get_rec, nfu_put_rec, nfu_def_dim, nfu_def_var, nfu_put_att, &
!      nfu_inq_var

use snowpack_mod, only: snowpack_init_lm4p2, read_snowpack_namelist


use snow_accessors_mod ! use everything

! use cohort_io_mod, only :  read_create_cohorts, create_cohort_dimension, &
   !   add_cohort_data, add_int_cohort_data, get_cohort_data, get_int_cohort_data
     use snowlayers_io_mod, only :  read_create_snowlayers, create_snowlayer_dimension, &
     add_snowlayer_data, add_int_snowlayer_data, get_snowlayer_data, get_int_snowlayer_data
   !   use snow_io_mod, only :  read_create_cohorts, create_cohort_dimension, &
   !   add_cohort_data, add_int_cohort_data, get_cohort_data, get_int_cohort_data

! use mpp_mod,          only : mpp_pe, mpp_max, mpp_send, mpp_recv, mpp_sync, &
!                              COMM_TAG_1, COMM_TAG_2, COMM_TAG_3, COMM_TAG_4, &
!                              mpp_sync_self, stdout

use snow_evolution_mod, only : gl_snow_step_2, gl_sweep_tiny_snow, &
         read_F06_data, gl_compute_snow_albedo, read_snow_evolution_namelist, use_internal_sources
use snowpack_mod, only : snow_layer_type, MAX_OPT_LAYERS

use fms_io_mod,       only : restart_file_type, get_instance_filename

! use fms_io_mod, only: fms_io_unstructured_register_restart_axis
! use fms_io_mod, only: fms_io_unstructured_register_restart_field
! use fms_io_mod, only: HIDX
! use fms_io_mod, only: fms_io_unstructured_read

implicit none
private

! ==== public interfaces =====================================================
public :: gl_read_snow_namelist
public :: gl_snow_init
public :: gl_snow_end
public :: gl_save_snow_restart
! public :: snow_get_sfc_temp
public :: gl_snow_get_depth_area
public :: gl_sweep_tiny_snow
! public :: snow_step_1
public :: gl_snow_step_2
! public :: compute_snow_albedo
! public :: use_internal_sources
! =====end of public interfaces ==============================================


! ==== module constants ======================================================
character(len=*), parameter :: module_name = 'gl_snow_mod'
#include "../../shared/version_variable.inc"

! name of the "compressed" dimension (and dimension variable) in the output
! netcdf files -- that is, the dimensions written out using compression by
! gathering, as described in CF conventions.
character(len=*),   parameter :: snowlayers_index_name   = 'snow_layer_index'

! ==== module variables ======================================================


abstract interface
  ! given land snow layer, returns pointer to some scalar real data
  ! within this snow layer, or an unassociated pointer if there is no data
  subroutine cptr_r0(tile, ptr)
     import snow_layer_type
     type(snow_layer_type), pointer :: tile ! input
     real                , pointer :: ptr  ! returned pointer to the data
  end subroutine cptr_r0
  ! given land snow layer, returns pointer to some scalar real data
  ! within this snow layer, or an unassociated pointer if there is no data
  subroutine cptr_i0(tile, ptr)
     import snow_layer_type
     type(snow_layer_type), pointer :: tile ! input
     integer               , pointer :: ptr  ! returned pointer to the data
  end subroutine cptr_i0
end interface

! ==== NetCDF declarations ===================================================
include 'netcdf.inc'
#define __NF_ASRT__(x) call print_netcdf_error((x),module_name,__LINE__)


!---- namelist ---------------------------------------------------------------
! logical :: retro_heat_capacity  = .false.
! logical :: lm2  = .false.
! logical :: steal = .false.
! character(len=16):: albedo_to_use = ''  ! or 'brdf-params'
! real :: max_snow       = 1000.
! real :: wet_max        = 0.0  ! TEMP, move to snow_data
! real :: snow_density   = 300. ! TEMP, move to snow_data and generalize
! real :: init_temp = 260.   ! cold-start snow T
! real :: init_pack_ws   =   0.
! real :: init_pack_wl   =   0.
! real :: min_snow_mass = 0.
! logical :: prevent_tiny_snow = .FALSE. ! if true, tiny snow is removed at the
   ! beginning of fast time step to avoid numerical issues. There is no harm
   ! in doing that, but it changes answers, so for compatibility with older code
   ! turn it off.

! namelist /gl_snow_nml/ retro_heat_capacity, lm2, steal, albedo_to_use, &
!                     max_snow, wet_max, snow_density, &
!                     init_temp, init_pack_ws, init_pack_wl, &
!                     min_snow_mass, prevent_tiny_snow
!---- end of namelist --------------------------------------------------------

logical         :: module_is_initialized =.FALSE.
real            :: delta_time
integer         :: num_l    ! # of snow layers
! next three 'z' variables are all normalized by total snow pack depth
real            :: dz (max_lev) ! relative thicknesses of layers
real            :: z  (max_lev) ! relative depths of layer bounds
real            :: zz (max_lev) ! relative depths of layer centers
real            :: heat_capacity_retro = 1.6e6
real            :: mc_fict

! ==== end of module variables ===============================================

contains

! ============================================================================
subroutine gl_read_snow_namelist()
  ! ---- local vars
  integer :: unit         ! unit for namelist i/o
  integer :: io           ! i/o status for the namelist
  integer :: ierr         ! error code, returned by i/o routines
  integer :: l            ! layer iterator

  call read_snow_data_namelist(num_l,dz,mc_fict)

!   call log_version(version, module_name, &
!   __FILE__)
! #ifdef INTERNAL_FILE_NML
!   read (input_nml_file, nml=gl_snow_nml, iostat=io)
!   ierr = check_nml_error(io, 'gl_snow_nml')
! #else
!   if (file_exist('input.nml')) then
!      unit = open_namelist_file()
!      ierr = 1;
!      do while (ierr /= 0)
!         read (unit, nml=gl_snow_nml, iostat=io, end=10)
!         ierr = check_nml_error (io, 'gl_snow_nml')
!      enddo
! 10   continue
!      call close_file (unit)
!   endif
! #endif
!   if (mpp_pe() == mpp_root_pe()) then
!      unit=stdlog()
!      write(unit, nml=gl_snow_nml)
!   endif

  ! -------- set up vertical discretization --------
  ! EZSNOW: add additional parameters here
!   EZSNOW : not needed anymore
!   zz(1) = 0
!   do l = 1, num_l
!      zz(l+1) = zz(l) + dz(l)
!      z(l)    = 0.5*(zz(l+1) + zz(l))
!   enddo

end subroutine gl_read_snow_namelist


! ============================================================================
! initialize snow model
subroutine gl_snow_init()

  ! ---- local vars ----------------------------------------------------------
  integer :: k
  type(land_tile_enum_type)     :: ce    ! tile list enumerator
  type(land_tile_type), pointer :: tile  ! pointer to current tile
  character(*), parameter :: restart_file_name='INPUT/snow.res.nc'
  type(land_restart_type) :: restart
  logical :: restart_exists
  integer ib, ik, ic, counter
  
  logical read_old_snow_restart
  real old_init_snow_density

  old_init_snow_density = 250.0 ! kg/m3, density for old model snow to input ! // TODO read from nml?
  module_is_initialized = .TRUE.
  delta_time = time_type_to_real(lnd%dt_fast)


  ! initialize snowpack -> namelist moved to snowpack_init
!   write(*,*) "read snowpack nml..."
  call read_snowpack_namelist()
!   write(*,*) "init snowpack module..."
  call snowpack_init_lm4p2()
  !read nml for snowpack mod

  ! read nml for snow evolution mod
!   write(*,*) "init snow evolution nml..."
  call read_snow_evolution_namelist()
  ! read data needed for Flanner and Zender 2006 snow aging
!   call read_snow_evolution_namelist()
!   write(*,*) "read Flanner data..."
   call read_snow_snicar_namelist()

   call read_snicar_optics_data()

  call read_F06_data()


  ! -------- initialize snow state --------
  call open_land_restart(restart,restart_file_name,restart_exists)
  if (restart_exists) then ! // TODO EZSNOW read restart with snowpack

   call error_mesg('gl_snow_init', 'reading NetCDF restart "'//trim(restart_file_name)//'"', NOTE)

   read_old_snow_restart = field_exists(restart,'temp')

   if (read_old_snow_restart) then
      ! READ OLD SNOW RESTART FILE
      write(*,*) "start reading old CM restart"
      call get_tile_data(restart, 'temp', 'zfull', cm_snow_temp_ptr)
      call get_tile_data(restart, 'wl'  , 'zfull', cm_snow_wl_ptr)
      call get_tile_data(restart, 'ws'  , 'zfull', cm_snow_ws_ptr)
      write(*,*) "done reading old CM restart"




      ! NOW PASS SNOW VARS TO NEW SNOW STRUCTURE
      ce = first_elmt(land_tile_map)
      do while(loop_over_tiles(ce, tile))
         ! if (.not.associated(tile%snow)) cycle
         if (.not.associated(tile%snow)) cycle
            ! tile%snow%nlayers = 0
            ! call tile%snow%sp%start(1)
               if (minval(tile%snow%wl) < -1E-7) then
                  write(*,*) "minval(tile%snow%wl) = ", minval(tile%snow%wl)
                  call error_mesg('gl_snow_init', 'Found wl < 0, routine should not be used in this case!', FATAL)
               endif

               tile%snow%sp%topwater = 0.0
               tile%snow%sp%topwheat = 0.0
               tile%snow%sp%snow_refl_dir = (/0.9,   0.9/)
               tile%snow%sp%snow_refl_dir = (/0.6,   0.6/)
               tile%snow%sp%beta_rad =      (/100.0, 40.0/)
            if (sum(tile%snow%ws)<0.0) then
               tile%snow%sp%nlayers = 0
               tile%snow%sp%topsnowdeficit = sum(tile%snow%ws)
               tile%snow%sp%topsnowheatdeficit =CSW * sum(tile%snow%ws * (tile%snow%T-TFREEZE))
            else
            ! tile%snow%sp%nlayers = 1
            tile%snow%sp%nlayers = 0 
            do ik=1,size(tile%snow%ws)
               if (tile%snow%ws(ik)>0.0) then
                  tile%snow%sp%nlayers = tile%snow%sp%nlayers + 1
               endif
            enddo

            if (tile%snow%sp%nlayers  > 0) then
               if (allocated(tile%snow%sp%snow)) deallocate(tile%snow%sp%snow) ! EZSNOW
               allocate(tile%snow%sp%snow(tile%snow%sp%nlayers))
               tile%snow%sp%topsnowdeficit = 0.0
               tile%snow%sp%topsnowheatdeficit = 0.0

               ! do ik = 1, tile%snow%sp%nlayers
               counter = 1
               do ik = 1, size(tile%snow%ws)

                  if (tile%snow%ws(ik) > 0.0) then
                 tile%snow%sp%snow(counter)%dz = tile%snow%ws(ik)/300.0 ! div by rho // use 250, TODO read from snow data
                 tile%snow%sp%snow(counter)%ws = tile%snow%ws(ik)
                 tile%snow%sp%snow(counter)%wl = tile%snow%wl(ik)
                 tile%snow%sp%snow(counter)%T = tile%snow%T(ik)
               !   tile%snow%sp%snow(ik)%dz = sum(tile%snow%ws)/300.0 ! div by rho
               !   tile%snow%sp%snow(ik)%ws = sum(tile%snow%ws)
               !   tile%snow%sp%snow(ik)%wl = sum(tile%snow%wl)
               !   tile%snow%sp%snow(ik)%dz = 0.1/300.0
               !   tile%snow%sp%snow(ik)%ws = 0.1
               !   tile%snow%sp%snow(ik)%wl = 0.0
               !   tile%snow%sp%snow(counter)%T = 273.15
                 tile%snow%sp%snow(counter)%dendr = 0.5
                 tile%snow%sp%snow(counter)%optd = 1E-4
                 tile%snow%sp%snow(counter)%sph = 0.5
                 tile%snow%sp%snow(counter)%age = 0.0
                 do ic = 1, 3
                   tile%snow%sp%snow(counter)%wc_em(ic) = 0.0
                   tile%snow%sp%snow(counter)%wc_im(ic) = 0.0
                 enddo
                  counter = counter + 1
               else
                  tile%snow%sp%topsnowdeficit = tile%snow%sp%topsnowdeficit + tile%snow%ws(ik)
                  tile%snow%sp%topsnowheatdeficit = tile%snow%sp%topsnowheatdeficit + CSW * tile%snow%ws(ik) * (tile%snow%T(ik) - TFREEZE)
               endif

               enddo
               endif
            endif
            ! call tile%snow%sp%start(0)
         ! else
            ! tile%snow%nlayers = 0
            ! tile%snow%sp%nlayers = 0
         ! call tile%snow%sp%start(size(tile%snow%ws), dz = tile%snow%ws/old_init_snow_density, &
         !       ws = tile%snow%ws, wl = tile%snow%wl, T = tile%snow%T) 
         deallocate(tile%snow%ws)
         deallocate(tile%snow%wl)
         deallocate(tile%snow%T)
         ! endif
      enddo
   else
      ! READ NEW SNOW RESTART FILE
      write(*,*) "start reading new GLASS restart"
      call read_create_snowlayers(restart)
      call get_snowlayer_data(restart, 'T',    snowlayer_T_ptr)
      call get_snowlayer_data(restart, 'wl',   snowlayer_wl_ptr)
      call get_snowlayer_data(restart, 'ws',   snowlayer_ws_ptr)
      call get_snowlayer_data(restart, 'dz',   snowlayer_dz_ptr)
      call get_snowlayer_data(restart, 'dendr',   snowlayer_dendr_ptr)
      call get_snowlayer_data(restart, 'optd',   snowlayer_optd_ptr)
      call get_snowlayer_data(restart, 'sph',   snowlayer_sph_ptr)
      call get_snowlayer_data(restart, 'age',   snowlayer_age_ptr)

      call get_snowlayer_data(restart,'wc_im_bc',snowlayer_wc_im_bc_ptr)
      call get_snowlayer_data(restart,'wc_im_md',snowlayer_wc_im_md_ptr)
      call get_snowlayer_data(restart,'wc_im_om',snowlayer_wc_im_om_ptr)
      call get_snowlayer_data(restart,'wc_em_bc',snowlayer_wc_em_bc_ptr)
      call get_snowlayer_data(restart,'wc_em_md',snowlayer_wc_em_md_ptr)
      call get_snowlayer_data(restart,'wc_em_om',snowlayer_wc_em_om_ptr)

      call get_int_tile_data(restart,'nlayers', snowtile_nlayers_ptr)
      call get_tile_data(restart,'topwater',snowtile_topwater_ptr)
      call get_tile_data(restart,'topwheat',snowtile_topwheat_ptr)
      call get_tile_data(restart,'topsnowdeficit',snowtile_topsnowdeficit_ptr)
      call get_tile_data(restart,'topsnowheatdeficit',snowtile_topsnowheatdeficit_ptr)

      call get_tile_data(restart, 'beta_rad', 'bands', beta_rad_ptr)

      call get_tile_data(restart,'nearsurf_bceq_tot', snowtile_nearsurf_bceq_tot_ptr)
      call get_tile_data(restart,'nearsurf_bceq_em', snowtile_nearsurf_bceq_em_ptr)
      call get_tile_data(restart,'nearsurf_bceq_im', snowtile_nearsurf_bceq_im_ptr)
      call get_tile_data(restart,'nearsurf_optd', snowtile_nearsurf_optd_ptr)
      call get_tile_data(restart,'nearsurf_sph', snowtile_nearsurf_sph_ptr)
      call get_tile_data(restart,'nearsurf_rho', snowtile_nearsurf_rho_ptr)
      call get_tile_data(restart,'nearsurf_age', snowtile_nearsurf_age_ptr)

      ! // TODO: Added other vars
      ! call get_tile_data(restart,'nearsurf_T', snowtile_nearsurf_T_ptr)
      ! call get_tile_data(restart,'nearsurf_dendr', snowtile_nearsurf_dendr_ptr)
      ! call get_tile_data(restart,'preprec', snowtile_nearsurf_T_ptr)
   endif

  else
     call error_mesg('gl_snow_init', 'cold-starting snow', NOTE)
     ce = first_elmt(land_tile_map)
     do while(loop_over_tiles(ce, tile))
         if (.not.associated(tile%snow)) cycle
            ! if (read_old_snow_restart) then
            !     if (.not.associated(tile%snow)) then
            !        write(*,*) "old CM restart - snow not associated"
            !        call tile%snow%sp%start(0) ! EZSNOW: pass old to new snowpack variables ...
            !     else
            !        write(*,*) "old CM restart - snow is associated"
            !        write(*,*) "old CM restart - snow is associated, num_l = ", num_l
            !        write(*,*) "old CM restart - snow is associated, len(ws) = ", size(tile%snow%ws)
            !        write(*,*) "old CM restart - snow is associated, len(wl) = ", size(tile%snow%wl)
            !        write(*,*) "old CM restart - snow is associated, len(T) = ", size(tile%snow%T)
            !        call tile%snow%sp%start(num_l, dz = tile%snow%ws/old_init_snow_density, &
            !                                       ws = tile%snow%ws, wl = tile%snow%wl, T = tile%snow%T) 
            !       !  write(*,*) "old CM restart - snow is associated - done reading values"
            !       !  write(*,*) "after reading: nlayers = ", tile%snow%sp%nlayers
            !     endif
            ! else
            tile%snow%sp%nlayers = 0 
            tile%snow%sp%topsnowdeficit = 0.0
            tile%snow%sp%topsnowheatdeficit = 0.0
            tile%snow%sp%topwater = 0.0
            tile%snow%sp%topwheat = 0.0
            tile%snow%sp%snow_refl_dir = (/0.9,   0.9/)
            tile%snow%sp%snow_refl_dir = (/0.6,   0.6/)
            tile%snow%sp%beta_rad =      (/100.0, 40.0/)
            ! call tile%snow%sp%start(0) ! EZSNOW: init empty snowpack ...
         ! endif
         ! Now remove the fields from old cm snow object - were needed only for reading old restart
         ! deallocate(tile%snow%ws)
         ! deallocate(tile%snow%wl)
         ! deallocate(tile%snow%T)
     enddo
  endif
  call free_land_restart(restart)

!   if (trim(albedo_to_use)=='') then
!      use_brdf = .false.
!   elseif (trim(albedo_to_use)=='brdf-params') then
!      use_brdf = .true.
!   else
!      call error_mesg('snow_init',&
!           'option albedo_to_use="'//&
!           trim(albedo_to_use)//'" is invalid, use "" or "brdf-params"',&
!           FATAL)
!   endif





end subroutine gl_snow_init


! ============================================================================
subroutine gl_snow_end ()

  module_is_initialized =.FALSE.

end subroutine gl_snow_end


! ============================================================================
subroutine gl_save_snow_restart(tile_dim_length,timestamp)
  integer, intent(in) :: tile_dim_length ! length of tile dim. in the output file
  character(*), intent(in) :: timestamp ! timestamp to add to the file name

  ! ---- local vars
  integer ::  i, j
  type(land_tile_enum_type) :: ce
  type(land_tile_type), pointer :: tile
!   integer :: n_accum, nmn_acm

  character(267) :: filename
  type(land_restart_type) :: restart1 ! restart file i/o object
!   type(land_restart_type) :: restart1, restart2 ! restart file i/o object
!   character:: spnames(fm_field_name_len, nspecies) ! names of the species

  call error_mesg('snow_end','writing NetCDF restart',NOTE)

! Note that filename is updated for tile & rank numbers during file creation
  filename = trim(timestamp)//'snow.res.nc'
  call init_land_restart(restart1, filename, snow_tile_exists, tile_dim_length)

  ! create output file, including internal structure necessary for tile output
!   filename = trim(timestamp)//'vegn1.res.nc'
!   call init_land_restart(restart1, filename, vegn_tile_exists, tile_dim_length)

  ! create compressed dimension for vegetation cohorts -- must be called even
  ! if restart has not been created, because it calls mpp_max and that should
  ! be called on all PEs to work
!   call create_cohort_dimension(restart1)
  call create_snowlayer_dimension(restart1)

   Write(*,*) "writing snow restart ..."
   Write(*,*) "writing snow restart : status of snow object:"
   ! call tile%snow%print()
   ! Write(*,*) "snow restart = ", restart1
   Write(*,*) "snow restart ... cidx = ", restart1%cidx

  call add_snowlayer_data(restart1,'T',snowlayer_T_ptr,'layer temperature','degrees_K')
  call add_snowlayer_data(restart1,'wl',snowlayer_wl_ptr,'layer liquid water content','kg/m2')
  call add_snowlayer_data(restart1,'ws',snowlayer_ws_ptr,'layer ice water content','kg/m2')
  call add_snowlayer_data(restart1,'dz',snowlayer_dz_ptr,'layer thickness','m')
  call add_snowlayer_data(restart1,'age',snowlayer_age_ptr,'layer age','days')
  call add_snowlayer_data(restart1,'optd',snowlayer_optd_ptr,'layer optical diameter','m')
  call add_snowlayer_data(restart1,'dendr',snowlayer_dendr_ptr,'layer snow dendricity','dimless')
  call add_snowlayer_data(restart1,'sph',snowlayer_sph_ptr,'layer snow sphericity','dimless')
  call add_snowlayer_data(restart1,'wc_im_bc',snowlayer_wc_im_bc_ptr,'LAI content', 'mg/m2')
  call add_snowlayer_data(restart1,'wc_im_md',snowlayer_wc_im_md_ptr,'LAI content', 'mg/m2')
  call add_snowlayer_data(restart1,'wc_im_om',snowlayer_wc_im_om_ptr,'LAI content', 'mg/m2')
  call add_snowlayer_data(restart1,'wc_em_bc',snowlayer_wc_em_bc_ptr,'LAI content', 'mg/m2')
  call add_snowlayer_data(restart1,'wc_em_md',snowlayer_wc_em_md_ptr,'LAI content', 'mg/m2')
  call add_snowlayer_data(restart1,'wc_em_om',snowlayer_wc_em_om_ptr,'LAI content', 'mg/m2')


   !   call add_snowlayer_data(restart1, 'e',        snowtile_e_ptr, 'e vector', 'units')
   !   call add_snowlayer_data(restart1, 'f',        snowtile_f_ptr, 'f vector', 'units')
   !   call add_snowlayer_data(restart1, 'swheat',   snowtile_swheat_ptr, 'swheat vector', 'units')
   !   call add_tile_data(restart1, 'e',        snowtile_e_ptr, 'e vector', 'units')
   !   call add_tile_data(restart1, 'f',        snowtile_f_ptr, 'f vector', 'units')
   !   call add_tile_data(restart1, 'swheat',   snowtile_swheat_ptr, 'swheat vector', 'units')


   ! add additional scalar quantities for the snowpack - snow tile element
   ! call add_scalar_data(restart1,'nlayers', snowtile_nlayers_ptr,'number of snow layers', 'number')
   ! call add_scalar_data(restart1,'topwater',snowtile_topwater_ptr,'snowpack topwater', 'kg/m2')
   ! call add_scalar_data(restart1,'topwheat',snowtile_topwheat_ptr,'snowpack topwheat', 'J/m2')
   ! call add_scalar_data(restart1,'topsnowdeficit',snowtile_topsnowdeficit_ptr,'snowpack top snow deficit', 'kg/m2')
   ! call add_scalar_data(restart1,'topsnowheatdeficit',snowtile_topsnowheatdeficit_ptr,'snowpack top snow heat deficit', 'J/m2')

      call add_restart_axis(restart1,'bands',(/ 1.0, 2.0 /),'NB',longname='shortwave bands',sense=-1)
      ! call add_tile_data(restart,'fpdir3d','bands', fpdir3d_ptr, 'direct flux correction','dimless')

      call add_tile_data(restart1,'beta_rad', 'bands', beta_rad_ptr, 'snow optical thickness', 'm^-1')
      ! call add_tile_data(restart1,'beta_rad', snowtile_beta_rad_ptr, 2, 'snow optical thickness NIR', 'm^-1')

   ! call add_tile_data(restart1,'betarad_VIS', snowtile_betarad_VIS_ptr,'snow optical thickness VIS', 'm^-1')
   ! call add_tile_data(restart1,'betarad_NIR', snowtile_betarad_NIR_ptr,'snow optical thickness NIR', 'm^-1')

   call add_int_tile_data(restart1,'nlayers', snowtile_nlayers_ptr,'number of snow layers', 'number')
   call add_tile_data(restart1,'topwater',snowtile_topwater_ptr,'snowpack topwater', 'kg/m2')
   call add_tile_data(restart1,'topwheat',snowtile_topwheat_ptr,'snowpack topwheat', 'J/m2')
   call add_tile_data(restart1,'topsnowdeficit',snowtile_topsnowdeficit_ptr,'snowpack top snow deficit', 'kg/m2')
   call add_tile_data(restart1,'topsnowheatdeficit',snowtile_topsnowheatdeficit_ptr,'snowpack top snow heat deficit', 'J/m2')

   call add_tile_data(restart1,'nearsurf_bceq_tot', snowtile_nearsurf_bceq_tot_ptr, 'nearsurf_bceq_tot','ppm')
   call add_tile_data(restart1,'nearsurf_bceq_em', snowtile_nearsurf_bceq_em_ptr, 'nearsurf_bceq_em','ppm')
   call add_tile_data(restart1,'nearsurf_bceq_im', snowtile_nearsurf_bceq_im_ptr, 'nearsurf_bceq_im','ppm')
   call add_tile_data(restart1,'nearsurf_optd', snowtile_nearsurf_optd_ptr, 'nearsurf_optd','m')
   call add_tile_data(restart1,'nearsurf_sph', snowtile_nearsurf_sph_ptr, 'nearsurf_sph','dimless')
   call add_tile_data(restart1,'nearsurf_rho', snowtile_nearsurf_rho_ptr, 'nearsurf_rho','kg/m3')
   call add_tile_data(restart1,'nearsurf_age', snowtile_nearsurf_age_ptr, 'nearsurf_age','days')


   ! call add_tile_data(restart1,'nearsurf_T', snowtile_nearsurf_T_ptr, 'nearsurf_T','K')
   ! call add_tile_data(restart1,'nearsurf_dendr', snowtile_nearsurf_dendr_ptr, 'nearsurf_dendr','dimless')

   call save_land_restart(restart1)
   call free_land_restart(restart1)

end subroutine gl_save_snow_restart

! ============================================================================
subroutine gl_snow_get_depth_area(snow, snow_depth, snow_area)
  type(gl_snow_tile_type), intent(in) :: snow
  real, intent(out) :: snow_depth, snow_area

  integer :: l

snow_depth = snow%sp%depth()

  call snow_data_area (snow_depth, snow_area )
end subroutine




! ============================================================================
subroutine gl_get_snow_integrals(snow, snow_LMASS, snow_FMASS, snow_HEAT)
  type(gl_snow_tile_type), intent(in) :: snow
  real, intent(out) :: snow_LMASS, snow_FMASS, snow_HEAT
!   integer :: l
snow_LMASS = snow%sp%liq() 
snow_FMASS = snow%sp%ice()
snow_HEAT = snow%sp%heat()
!   snow_LMASS = 0; snow_FMASS = 0; snow_HEAT = 0
!   do l = 1, num_l;
!     snow_LMASS = snow_LMASS + snow%wl(l)
!     snow_FMASS = snow_FMASS + snow%ws(l)
!     snow_HEAT = snow_HEAT + &
!       (mc_fict*dz(l) + clw*snow%wl(l) + csw*snow%ws(l))  &
!                                             * (snow%T(l)-tfreeze)
!   enddo
end subroutine gl_get_snow_integrals

! ============================================================================
subroutine gl_print_snow_integrals(snow)
  type(gl_snow_tile_type), intent(in) :: snow

  real    :: snow_LMASS, snow_FMASS, snow_HEAT
  call gl_get_snow_integrals(snow, snow_LMASS, snow_FMASS, snow_HEAT)
  __DEBUG3__(snow_LMASS, snow_FMASS, snow_HEAT)
end subroutine gl_print_snow_integrals

! ============================================================================
! tile existence detector: returns a logical value indicating wether component
! model tile exists or not
logical function snow_tile_exists(tile)
   type(land_tile_type), pointer :: tile
   snow_tile_exists = associated(tile%snow)
end function snow_tile_exists

! ============================================================================
! accessor functions: given a pointer to a land tile, they return pointer
! to the desired member of the land tile, of NULL if this member does not
! exist.
subroutine snow_temp_ptr(tile, i, ptr)
   type(land_tile_type), pointer :: tile ! input
   integer             , intent(in) :: i ! index in the array
   real                , pointer :: ptr  ! returned pointer to the data
   ptr=>NULL()
   if(associated(tile)) then
      ! if(associated(tile%snow)) ptr => tile%snow%T(i)
      if(associated(tile%snow)) ptr => tile%snow%sp%snow(i)%T ! EZSNOW
   endif
end subroutine snow_temp_ptr

subroutine snow_wl_ptr(tile, i, ptr)
   type(land_tile_type), pointer :: tile ! input
   integer             , intent(in) :: i ! index in the array
   real                , pointer :: ptr  ! returned pointer to the data
   ptr=>NULL()
   if(associated(tile)) then
      ! if(associated(tile%snow)) ptr => tile%snow%wl(i)
      if(associated(tile%snow)) ptr => tile%snow%sp%snow(i)%wl ! EZSNOW
   endif
end subroutine snow_wl_ptr

subroutine snow_ws_ptr(tile, i, ptr)
   type(land_tile_type), pointer :: tile ! input
   integer             , intent(in) :: i ! index in the array
   real                , pointer :: ptr  ! returned pointer to the data
   ptr=>NULL()
   if(associated(tile)) then
      ! if(associated(tile%snow)) ptr => tile%snow%ws(i)
      if(associated(tile%snow)) ptr => tile%snow%sp%snow(i)%ws ! EZSNOW
   endif
end subroutine snow_ws_ptr

! EZSNOW: pointers to additional snow layer fields

subroutine snow_dz_ptr(tile, i, ptr)
   type(land_tile_type), pointer :: tile ! input
   integer             , intent(in) :: i ! index in the array
   real                , pointer :: ptr  ! returned pointer to the data
   ptr=>NULL()
   if(associated(tile)) then
      if(associated(tile%snow)) ptr => tile%snow%sp%snow(i)%dz ! EZSNOW
   endif
end subroutine snow_dz_ptr

subroutine snow_optd_ptr(tile, i, ptr)
   type(land_tile_type), pointer :: tile ! input
   integer             , intent(in) :: i ! index in the array
   real                , pointer :: ptr  ! returned pointer to the data
   ptr=>NULL()
   if(associated(tile)) then
      if(associated(tile%snow)) ptr => tile%snow%sp%snow(i)%optd ! EZSNOW
   endif
end subroutine snow_optd_ptr

subroutine snow_dendr_ptr(tile, i, ptr)
   type(land_tile_type), pointer :: tile ! input
   integer             , intent(in) :: i ! index in the array
   real                , pointer :: ptr  ! returned pointer to the data
   ptr=>NULL()
   if(associated(tile)) then
      if(associated(tile%snow)) ptr => tile%snow%sp%snow(i)%dendr ! EZSNOW
   endif
end subroutine snow_dendr_ptr

subroutine snow_age_ptr(tile, i, ptr)
   type(land_tile_type), pointer :: tile ! input
   integer             , intent(in) :: i ! index in the array
   real                , pointer :: ptr  ! returned pointer to the data
   ptr=>NULL()
   if(associated(tile)) then
      if(associated(tile%snow)) ptr => tile%snow%sp%snow(i)%age ! EZSNOW
   endif
end subroutine snow_age_ptr

subroutine snow_sph_ptr(tile, i, ptr)
   type(land_tile_type), pointer :: tile ! input
   integer             , intent(in) :: i ! index in the array
   real                , pointer :: ptr  ! returned pointer to the data
   ptr=>NULL()
   if(associated(tile)) then
      if(associated(tile%snow)) ptr => tile%snow%sp%snow(i)%sph ! EZSNOW
   endif
end subroutine snow_sph_ptr

subroutine beta_rad_ptr(tile, i, ptr)
   type(land_tile_type), pointer :: tile ! input
   integer             , intent(in) :: i ! index in the array
   real                , pointer :: ptr  ! returned pointer to the data
   ptr=>NULL()
   if(associated(tile)) then
      if(associated(tile%snow)) ptr => tile%snow%sp%beta_rad(i) ! EZSNOW
   endif
end subroutine beta_rad_ptr


! ============================================================================
! old snow accessor functions: given a pointer to a land tile, they return pointer
! to the desired member of the land tile, of NULL if this member does not
! exist.
subroutine cm_snow_temp_ptr(tile, i, ptr)
   type(land_tile_type), pointer :: tile ! input
   integer             , intent(in) :: i ! index in the array
   real                , pointer :: ptr  ! returned pointer to the data
   ptr=>NULL()
   if(associated(tile)) then
      if(associated(tile%snow)) ptr => tile%snow%T(i)
   endif
end subroutine cm_snow_temp_ptr

subroutine cm_snow_wl_ptr(tile, i, ptr)
   type(land_tile_type), pointer :: tile ! input
   integer             , intent(in) :: i ! index in the array
   real                , pointer :: ptr  ! returned pointer to the data
   ptr=>NULL()
   if(associated(tile)) then
      if(associated(tile%snow)) ptr => tile%snow%wl(i)
   endif
end subroutine cm_snow_wl_ptr

subroutine cm_snow_ws_ptr(tile, i, ptr)
   type(land_tile_type), pointer :: tile ! input
   integer             , intent(in) :: i ! index in the array
   real                , pointer :: ptr  ! returned pointer to the data
   ptr=>NULL()
   if(associated(tile)) then
      if(associated(tile%snow)) ptr => tile%snow%ws(i)
   endif
end subroutine cm_snow_ws_ptr

#define F90_TYPE real
#define NF_TYPE NF_DOUBLE
#define NF_FILL_VALUE NF_FILL_DOUBLE
#define READ_0D_FPTR read_snowlayer_data_r0d_fptr
#define WRITE_0D_FPTR write_snowlayer_data_r0d_fptr
#define WRITE_0D write_snowlayer_data_r0d
! #include "vegn_snowlayer_io.inc"

#define F90_TYPE integer
#define NF_TYPE NF_INT
#define NF_FILL_VALUE NF_FILL_INT
#define READ_0D_FPTR read_snowlayer_data_i0d_fptr
#define WRITE_0D_FPTR write_snowlayer_data_i0d_fptr
#define WRITE_0D write_snowlayer_data_i0d
! #include "vegn_snowlayer_io.inc"

end module gl_snow_mod



