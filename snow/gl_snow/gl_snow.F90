
module gl_snow_mod

#include "../../shared/debug.inc"

use fms_mod, only : error_mesg, FATAL, NOTE
use time_manager_mod,   only: time_type_to_real
use constants_mod,      only: tfreeze
use land_tile_mod,    only : land_tile_map, land_tile_type, land_tile_list_type, &
     land_tile_enum_type, first_elmt, loop_over_tiles
use land_data_mod, only : lnd, log_version
use land_tile_io_mod, only: land_restart_type, &
     init_land_restart, open_land_restart, save_land_restart, free_land_restart, &
     add_restart_axis, add_tile_data, get_tile_data, get_tile_by_idx, &
     add_int_tile_data, get_int_tile_data, field_exists

use snow_tile_mod, only : &
     read_snow_data_namelist, snow_data_area, max_lev
use gl_snow_tile_mod, only: gl_snow_tile_type
use snowpack_mod, only: snowpack_init_lm4p2, read_snowpack_namelist, &
     snow_layer_type, csw
use snowlayers_io_mod, only :  read_create_snowlayers, create_snowlayer_dimension, &
     add_snowlayer_data, add_int_snowlayer_data, get_snowlayer_data, get_int_snowlayer_data
use snow_evolution_mod, only : &
     read_F06_data, read_snow_evolution_namelist
use snicar_mod, only: read_snicar_optics_data, read_snow_snicar_namelist


implicit none
private

! ==== public interfaces =====================================================
public :: gl_read_snow_namelist
public :: gl_snow_init
public :: gl_snow_end
public :: gl_save_snow_restart
! =====end of public interfaces ==============================================


! ==== module constants ======================================================
character(len=*), parameter :: module_name = 'gl_snow_mod'
#include "../../shared/version_variable.inc"

! name of the "compressed" dimension (and dimension variable) in the output
! netcdf files -- that is, the dimensions written out using compression by
! gathering, as described in CF conventions.
character(len=*),   parameter :: snowlayers_index_name   = 'snow_layer_index'

! ---- module variables
logical         :: module_is_initialized =.FALSE.
! ---- end of module variables

contains

! ============================================================================
subroutine gl_read_snow_namelist()
  ! local variables only to satisfy interface of read_snow_data_namelist
  integer :: num_l    ! # of snow layers
  real    :: dz (max_lev) ! relative thicknesses of layers
  real    :: mc_fict

  call read_snow_data_namelist(num_l,dz,mc_fict)
end subroutine gl_read_snow_namelist


! ============================================================================
! initialize snow model
subroutine gl_snow_init()

  ! ---- local vars ----------------------------------------------------------
  type(land_tile_enum_type)     :: ce    ! tile list enumerator
  type(land_tile_type), pointer :: tile  ! pointer to current tile
  ! character(*), parameter :: restart_file_name='INPUT/snow.res.nc' ! OLD-VERSION
  character(*), parameter :: restart_file_name='INPUT/snow.nc' ! EZSNOW-2022SC
  type(land_restart_type) :: restart
  logical :: restart_exists
  integer ik, ic, counter

  logical read_old_snow_restart
  real old_init_snow_density

  old_init_snow_density = 250.0 ! kg/m3, density for old model snow to input ! // TODO read from nml?
  module_is_initialized = .TRUE.

  ! initialize snowpack -> namelist moved to snowpack_init
  call read_snowpack_namelist()
  call snowpack_init_lm4p2()

  call read_snow_evolution_namelist() ! read snow evolution namelist parameters
  call read_snow_snicar_namelist() ! read SNICAR namelist parameters
!   call read_snicar_optics_data() ! read SNICAR snow optical data ! already done in read nml
  call read_F06_data() ! Read flanner 2006 parameter table


  ! -------- initialize snow state --------
  call open_land_restart(restart,restart_file_name,restart_exists)
  if (restart_exists) then

   call error_mesg('gl_snow_init', 'reading NetCDF restart "'//trim(restart_file_name)//'"', NOTE)

   read_old_snow_restart = field_exists(restart,'temp')

   if (read_old_snow_restart) then
      write(*,*) "Snow GLASS :: reading old snow CM model restart"
      call get_tile_data(restart, 'temp', 'zfull', cm_snow_temp_ptr)
      call get_tile_data(restart, 'wl'  , 'zfull', cm_snow_wl_ptr)
      call get_tile_data(restart, 'ws'  , 'zfull', cm_snow_ws_ptr)




      ! now pass snow variables to new snow structure
      ce = first_elmt(land_tile_map)
      do while(loop_over_tiles(ce, tile))
         if (.not.associated(tile%snow)) cycle
            if (minval(tile%snow%wl) < -1E-7) then
               write(*,*) "minval(tile%snow%wl) = ", minval(tile%snow%wl)
               call error_mesg('gl_snow_init', 'Found wl < 0, routine should not be used in this case!', FATAL)
            endif
            tile%snow%sp%topwater = 0.0
            tile%snow%sp%topwheat = 0.0
            tile%snow%sp%snow_refl_dir = (/0.9,   0.9/) ! //TODO clean up
            tile%snow%sp%snow_refl_dir = (/0.6,   0.6/)
            tile%snow%sp%beta_rad =      (/100.0, 40.0/)
            if (sum(tile%snow%ws)<0.0) then
               tile%snow%sp%nlayers = 0
               tile%snow%sp%topsnowdeficit = sum(tile%snow%ws)
               tile%snow%sp%topsnowheatdeficit =CSW * sum(tile%snow%ws * (tile%snow%T-TFREEZE))
            else
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

               counter = 1
               do ik = 1, size(tile%snow%ws)
                  if (tile%snow%ws(ik) > 0.0) then
                     tile%snow%sp%snow(counter)%dz = tile%snow%ws(ik)/300.0 ! div by rho // use 250, //TODO read from snow data
                     tile%snow%sp%snow(counter)%ws = tile%snow%ws(ik)
                     tile%snow%sp%snow(counter)%wl = tile%snow%wl(ik)
                     tile%snow%sp%snow(counter)%T = tile%snow%T(ik)
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
         deallocate(tile%snow%ws)
         deallocate(tile%snow%wl)
         deallocate(tile%snow%T)
      enddo
   else
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
            tile%snow%sp%nlayers = 0
            tile%snow%sp%topsnowdeficit = 0.0
            tile%snow%sp%topsnowheatdeficit = 0.0
            tile%snow%sp%topwater = 0.0
            tile%snow%sp%topwheat = 0.0
            tile%snow%sp%snow_refl_dir = (/0.9,   0.9/)
            tile%snow%sp%snow_refl_dir = (/0.6,   0.6/)
            tile%snow%sp%beta_rad =      (/100.0, 40.0/)
     enddo
  endif
  call free_land_restart(restart)


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

  character(267) :: filename
  type(land_restart_type) :: restart1 ! restart file i/o object

  call error_mesg('snow_end','writing NetCDF restart',NOTE)

! Note that filename is updated for tile & rank numbers during file creation
  ! filename = trim(timestamp)//'snow.res.nc' ! OLD VERSION
  filename = 'RESTART/'//trim(timestamp)//'snow.nc' ! EZSNOW-2022SC
  call init_land_restart(restart1, filename, snow_tile_exists, tile_dim_length)

  ! create compressed dimension for snow layers -- must be called even
  ! if restart has not been created, because it calls mpp_max and that should
  ! be called on all PEs to work
  call create_snowlayer_dimension(restart1)

  Write(*,*) "writing snow restart : status of snow object:" !//TODO clean up
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

  ! call add_restart_axis(restart1,'bands',(/ 1.0, 2.0 /),'NB',longname='shortwave bands',sense=-1) ! OLD VERSION
  call add_restart_axis(restart1,'bands',(/ 1.0, 2.0 /),.false., 'NB',longname='shortwave bands',sense=-1) ! EZSNOW-2022SC

  call add_tile_data(restart1,'beta_rad', 'bands', beta_rad_ptr, 'snow optical thickness', 'm^-1')
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
      if(associated(tile%snow)) ptr => tile%snow%sp%snow(i)%T
   endif
end subroutine snow_temp_ptr

subroutine snow_wl_ptr(tile, i, ptr)
   type(land_tile_type), pointer :: tile ! input
   integer             , intent(in) :: i ! index in the array
   real                , pointer :: ptr  ! returned pointer to the data
   ptr=>NULL()
   if(associated(tile)) then
      if(associated(tile%snow)) ptr => tile%snow%sp%snow(i)%wl
   endif
end subroutine snow_wl_ptr

subroutine snow_ws_ptr(tile, i, ptr)
   type(land_tile_type), pointer :: tile ! input
   integer             , intent(in) :: i ! index in the array
   real                , pointer :: ptr  ! returned pointer to the data
   ptr=>NULL()
   if(associated(tile)) then
      if(associated(tile%snow)) ptr => tile%snow%sp%snow(i)%ws
   endif
end subroutine snow_ws_ptr

subroutine snow_dz_ptr(tile, i, ptr)
   type(land_tile_type), pointer :: tile ! input
   integer             , intent(in) :: i ! index in the array
   real                , pointer :: ptr  ! returned pointer to the data
   ptr=>NULL()
   if(associated(tile)) then
      if(associated(tile%snow)) ptr => tile%snow%sp%snow(i)%dz
   endif
end subroutine snow_dz_ptr

subroutine snow_optd_ptr(tile, i, ptr)
   type(land_tile_type), pointer :: tile ! input
   integer             , intent(in) :: i ! index in the array
   real                , pointer :: ptr  ! returned pointer to the data
   ptr=>NULL()
   if(associated(tile)) then
      if(associated(tile%snow)) ptr => tile%snow%sp%snow(i)%optd
   endif
end subroutine snow_optd_ptr

subroutine snow_dendr_ptr(tile, i, ptr)
   type(land_tile_type), pointer :: tile ! input
   integer             , intent(in) :: i ! index in the array
   real                , pointer :: ptr  ! returned pointer to the data
   ptr=>NULL()
   if(associated(tile)) then
      if(associated(tile%snow)) ptr => tile%snow%sp%snow(i)%dendr
   endif
end subroutine snow_dendr_ptr

subroutine snow_age_ptr(tile, i, ptr)
   type(land_tile_type), pointer :: tile ! input
   integer             , intent(in) :: i ! index in the array
   real                , pointer :: ptr  ! returned pointer to the data
   ptr=>NULL()
   if(associated(tile)) then
      if(associated(tile%snow)) ptr => tile%snow%sp%snow(i)%age
   endif
end subroutine snow_age_ptr

subroutine snow_sph_ptr(tile, i, ptr)
   type(land_tile_type), pointer :: tile ! input
   integer             , intent(in) :: i ! index in the array
   real                , pointer :: ptr  ! returned pointer to the data
   ptr=>NULL()
   if(associated(tile)) then
      if(associated(tile%snow)) ptr => tile%snow%sp%snow(i)%sph
   endif
end subroutine snow_sph_ptr

subroutine beta_rad_ptr(tile, i, ptr)
   type(land_tile_type), pointer :: tile ! input
   integer             , intent(in) :: i ! index in the array
   real                , pointer :: ptr  ! returned pointer to the data
   ptr=>NULL()
   if(associated(tile)) then
      if(associated(tile%snow)) ptr => tile%snow%sp%beta_rad(i)
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

! ============================================================================
! snowlayer accessor functions: given a pointer to a snowlayer, return a pointer to a
! specific member of the snowlayer structure

#define DEFINE_SNOWPACK_ACCESSOR_0D(xtype,x) subroutine snowtile_ ## x ## _ptr(t,p);\
type(land_tile_type),pointer::t;xtype,pointer::p;p=>NULL();if(associated(t))then;if(associated(t%snow))p=>t%snow%sp%x;endif;end subroutine

#define DEFINE_SNOWLAYER_ACCESSOR(xtype,x) subroutine snowlayer_ ## x ## _ptr(c,p);\
type(snow_layer_type),pointer::c;xtype,pointer::p;p=>NULL();if(associated(c))p=>c%x;end subroutine

#define DEFINE_SNOWLAYER_ACCESSOR_BC(xtype,x) subroutine snowlayer_ ## x ## _bc_ptr(c,p);\
type(snow_layer_type),pointer::c;xtype,pointer::p;p=>NULL();if(associated(c))p=>c%x(1);end subroutine

#define DEFINE_SNOWLAYER_ACCESSOR_MD(xtype,x) subroutine snowlayer_ ## x ## _md_ptr(c,p);\
type(snow_layer_type),pointer::c;xtype,pointer::p;p=>NULL();if(associated(c))p=>c%x(2);end subroutine

#define DEFINE_SNOWLAYER_ACCESSOR_OM(xtype,x) subroutine snowlayer_ ## x ## _om_ptr(c,p);\
type(snow_layer_type),pointer::c;xtype,pointer::p;p=>NULL();if(associated(c))p=>c%x(3);end subroutine



DEFINE_SNOWPACK_ACCESSOR_0D(integer, nlayers)
DEFINE_SNOWPACK_ACCESSOR_0D(real, topwater)
DEFINE_SNOWPACK_ACCESSOR_0D(real, topwheat)
DEFINE_SNOWPACK_ACCESSOR_0D(real, topsnowdeficit)
DEFINE_SNOWPACK_ACCESSOR_0D(real, topsnowheatdeficit)

DEFINE_SNOWPACK_ACCESSOR_0D(real, nearsurf_bceq_tot)
DEFINE_SNOWPACK_ACCESSOR_0D(real, nearsurf_bceq_em)
DEFINE_SNOWPACK_ACCESSOR_0D(real, nearsurf_bceq_im)
DEFINE_SNOWPACK_ACCESSOR_0D(real, nearsurf_dendr)
DEFINE_SNOWPACK_ACCESSOR_0D(real, nearsurf_optd)
DEFINE_SNOWPACK_ACCESSOR_0D(real, nearsurf_sph)
DEFINE_SNOWPACK_ACCESSOR_0D(real, nearsurf_rho)
DEFINE_SNOWPACK_ACCESSOR_0D(real, nearsurf_age)
DEFINE_SNOWPACK_ACCESSOR_0D(real, nearsurf_T)


DEFINE_SNOWLAYER_ACCESSOR(real,T)
DEFINE_SNOWLAYER_ACCESSOR(real,wl)
DEFINE_SNOWLAYER_ACCESSOR(real,ws)
DEFINE_SNOWLAYER_ACCESSOR(real,dz)
DEFINE_SNOWLAYER_ACCESSOR(real,age)
DEFINE_SNOWLAYER_ACCESSOR(real,dendr)
DEFINE_SNOWLAYER_ACCESSOR(real,optd)
DEFINE_SNOWLAYER_ACCESSOR(real,sph)



DEFINE_SNOWLAYER_ACCESSOR_BC(real,wc_im)
DEFINE_SNOWLAYER_ACCESSOR_MD(real,wc_im)
DEFINE_SNOWLAYER_ACCESSOR_OM(real,wc_im)
DEFINE_SNOWLAYER_ACCESSOR_BC(real,wc_em)
DEFINE_SNOWLAYER_ACCESSOR_MD(real,wc_em)
DEFINE_SNOWLAYER_ACCESSOR_OM(real,wc_em)

end module gl_snow_mod



