#include <fms_platform.h>

module snow_tile_mod

use fms_mod, only : &
     write_version_number, file_exist, open_namelist_file, check_nml_error, &
     close_file, stdlog
use constants_mod,only: tfreeze
use land_constants_mod, only : &
     NBANDS
use land_tile_selectors_mod, only : &
     tile_selector_type

implicit none
private

! ==== public interfaces =====================================================
public :: snow_prog_type
public :: snow_tile_type

public :: new_snow_tile, delete_snow_tile
public :: snow_tiles_can_be_merged, merge_snow_tiles
public :: snow_is_selected
public :: get_snow_tile_tag

public :: read_snow_data_namelist

public :: snow_data_thermodynamics
public :: snow_data_hydraulics
public :: snow_data_area
public :: snow_data_radiation
public :: snow_data_diffusion

public :: max_lev
public :: cpw, clw, csw
! ==== end of public interfaces ==============================================
interface new_snow_tile
   module procedure snow_tile_ctor
   module procedure snow_tile_copy_ctor
end interface


! ==== module constants ======================================================
character(len=*), parameter :: &
     module_name = 'snow_tile_mod' ,&
     version     = '$Id: snow_tile.F90,v 15.0.2.2 2007/12/05 19:41:35 slm Exp $' ,&
     tagname     = '$Name: omsk_2008_03 $'
integer, parameter :: max_lev = 10
real   , parameter :: t_range = 10.0 ! degK

! ==== types =================================================================
type :: snow_prog_type
  real wl
  real ws
  real T
end type snow_prog_type


type :: snow_tile_type
   integer :: tag ! kind of the tile
   type(snow_prog_type), pointer :: prog(:)
   real,                 pointer :: e(:), f(:)
end type snow_tile_type

! ==== module data ===========================================================

!---- namelist ---------------------------------------------------------------
logical :: use_mcm_masking       = .false.   ! MCM snow mask fn
real    :: w_sat                 = 670.
real    :: psi_sat               = -0.06
real    :: k_sat                 = 0.02
real    :: chb                   = 3.5
real    :: thermal_cond_ref      = 0.3
real    :: depth_crit            = 0.0167
real    :: z0_momentum           = 0.001
real    :: refl_snow_max_dir(NBANDS) = (/ 0.8,  0.8  /) ! reset to 0.6 for MCM
real    :: refl_snow_max_dif(NBANDS) = (/ 0.8,  0.8  /) ! reset to 0.6 for MCM
real    :: refl_snow_min_dir(NBANDS) = (/ 0.65, 0.65 /) ! reset to 0.45 for MCM
real    :: refl_snow_min_dif(NBANDS) = (/ 0.65, 0.65 /) ! reset to 0.45 for MCM
real    :: emis_snow_max         = 0.95      ! reset to 1 for MCM
real    :: emis_snow_min         = 0.90      ! reset to 1 for MCM
real    :: k_over_B              = 2         ! reset to 0 for MCM
integer :: num_l                 = 3         ! number of snow levels
real    :: dz(max_lev)           = (/0.1,0.8,0.1,0.,0.,0.,0.,0.,0.,0./)
                                              ! rel. thickness of model layers,
                                              ! from top down
real    :: cpw = 1952.  ! specific heat of water vapor at constant pressure
real    :: clw = 4218.  ! specific heat of water (liquid)
real    :: csw = 2106.  ! specific heat of water (ice)
namelist /snow_data_nml/use_mcm_masking,    w_sat,                 &
                    psi_sat,                k_sat,                 &
                    chb,                                           &
                    thermal_cond_ref,       depth_crit,            &
                    z0_momentum,                                   &
                    refl_snow_max_dir,    refl_snow_min_dir,   &
                    refl_snow_max_dif,    refl_snow_min_dif,   &
                    emis_snow_max,          emis_snow_min,         &
                    k_over_B,             &
                    num_l,                   dz, cpw, clw, csw
!---- end of namelist --------------------------------------------------------

contains ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

! ============================================================================
subroutine read_snow_data_namelist(snow_num_l, snow_dz)
  integer, intent(out) :: snow_num_l
  real,    intent(out) :: snow_dz(:)
  ! ---- local vars
  integer :: unit         ! unit for namelist i/o
  integer :: io           ! i/o status for the namelist
  integer :: ierr         ! error code, returned by i/o routines

  call write_version_number(version, tagname)
  if (file_exist('input.nml')) then
     unit = open_namelist_file()
     ierr = 1;  
     do while (ierr /= 0)
        read (unit, nml=snow_data_nml, iostat=io, end=10)
        ierr = check_nml_error (io, 'snow_data_nml')
     enddo
10   continue
     call close_file (unit)
  endif
  write (stdlog(), nml=snow_data_nml)

  ! initialize global module data here

  ! set up output arguments
  snow_num_l = num_l
  snow_dz    = dz

end subroutine 

! ============================================================================
function snow_tile_ctor(tag) result(ptr)
  type(snow_tile_type), pointer :: ptr ! return value
  integer, optional, intent(in) :: tag ! kind of tile

  allocate(ptr)
  ptr%tag = 0 ; if(present(tag)) ptr%tag = tag
  ! allocate storage for tile data
  allocate(ptr%prog(num_l))
  allocate(ptr%e(num_l))
  allocate(ptr%f(num_l))

end function snow_tile_ctor

! ============================================================================
function snow_tile_copy_ctor(snow) result(ptr)
  type(snow_tile_type), pointer :: ptr ! return value
  type(snow_tile_type), intent(in) :: snow ! tile to copy

  allocate(ptr)
  ! copy all non-pointer members
  ptr = snow
  ! allocate storage for tile data
  allocate(ptr%prog(num_l))
  allocate(ptr%e(num_l))
  allocate(ptr%f(num_l))
  ! copy all pointer members
  ptr%prog(:) = snow%prog(:)
  ptr%e(:) = snow%e(:)
  ptr%f(:) = snow%f(:)
end function snow_tile_copy_ctor

! ============================================================================
subroutine delete_snow_tile(snow)
  type(snow_tile_type), pointer :: snow

  deallocate(snow%prog)
  deallocate(snow%e)
  deallocate(snow%f)
  deallocate(snow)
end subroutine delete_snow_tile

! =============================================================================
function snow_tiles_can_be_merged(snow1,snow2) result(response)
  logical :: response
  type(snow_tile_type), intent(in) :: snow1,snow2

  response = .TRUE.
end function

! =============================================================================
subroutine merge_snow_tiles(snow1, w1, snow2, w2)
  type(snow_tile_type), intent(in)    :: snow1
  type(snow_tile_type), intent(inout) :: snow2
  real                , intent(in)    :: w1, w2 ! relative weights
  
  ! ---- local vars
  real    :: x1, x2 ! normalized weights
  real    :: heat
  integer :: i
  
  ! calculate normalized weights
  x1 = w1/(w1+w2)
  x2 = 1-x1
  
  do i = 1, num_l
    heat = &
      (clw*snow1%prog(i)%wl+csw*snow1%prog(i)%ws)*(snow1%prog(i)%T-tfreeze)*x1 + &
      (clw*snow2%prog(i)%wl+csw*snow2%prog(i)%ws)*(snow2%prog(i)%T-tfreeze)*x2
    snow2%prog(i)%wl = snow1%prog(i)%wl*x1 + snow2%prog(i)%wl*x2
    snow2%prog(i)%ws = snow1%prog(i)%ws*x1 + snow2%prog(i)%ws*x2
    if (snow2%prog(i)%wl/=0.or.snow2%prog(i)%ws/=0) then
       snow2%prog(i)%T  = heat/(clw*snow2%prog(i)%wl+csw*snow2%prog(i)%ws)+tfreeze
    else
       snow2%prog(i)%T  = snow1%prog(i)%T*x1 + snow2%prog(i)%T*x2
    endif
  enddo
end subroutine

! =============================================================================
! returns true if tile fits the specified selector
function snow_is_selected(snow, sel)
  logical snow_is_selected
  type(tile_selector_type),  intent(in) :: sel
  type(snow_tile_type),      intent(in) :: snow

  snow_is_selected = .TRUE.
end function

! ============================================================================
! retruns tag of the tile
function get_snow_tile_tag(snow) result(tag)
  integer :: tag
  type(snow_tile_type), intent(in) :: snow
  
  tag = snow%tag
end function

! ============================================================================
! compute snow thermodynmamic properties.
subroutine snow_data_thermodynamics ( snow_rh, thermal_cond)
  real, intent(out) :: snow_rh
  real, intent(out) :: thermal_cond(:)

  ! snow surface assumed to have air at saturation
  snow_rh = 1

  ! these will eventually be functions of water contents and T.
  thermal_cond  = thermal_cond_ref

end subroutine 


! ============================================================================
! compute snow hydraulic properties (assumed dependent only on wl)
subroutine snow_data_hydraulics (wl, ws, psi, hyd_cond )
  real, intent(in),  dimension(:) :: wl, ws
  real, intent(out), dimension(:) :: psi, hyd_cond

  ! ---- local vars 
  integer :: l
  
  do l = 1, num_l
    psi     (l) = psi_sat *(w_sat/(wl(l)+ws(l)))**chb
    hyd_cond(l) = k_sat*(wl(l)/w_sat)**(3+2*chb)
  enddo

end subroutine snow_data_hydraulics


! ============================================================================
! compute snow area
subroutine snow_data_area ( snow_depth, snow_area )
    real, intent(in)  :: snow_depth
    real, intent(out) :: snow_area

  snow_area = 0.
  if (use_mcm_masking) then
     snow_area = min(1., 0.5*sqrt(max(0.,snow_depth)/depth_crit))
  else
     snow_area = max(0.,snow_depth) / (max(0.,snow_depth) + depth_crit)
  endif

end subroutine


! ============================================================================
subroutine snow_data_radiation(snow_T, snow_refl_dir, snow_refl_dif, snow_emis)
  real, intent(in) :: snow_T
  real, intent(out):: snow_refl_dir(NBANDS), snow_refl_dif(NBANDS), snow_emis

  ! ---- local vars
  real :: blend

  blend = max(0.,min(1.,1.-(tfreeze-snow_T)/t_range))
  snow_refl_dir = refl_snow_max_dir + blend*(refl_snow_min_dir-refl_snow_max_dir)
  snow_refl_dif = refl_snow_max_dif + blend*(refl_snow_min_dif-refl_snow_max_dif)
  snow_emis     = emis_snow_max + blend*(emis_snow_min-emis_snow_max  )
end subroutine


! ============================================================================
subroutine snow_data_diffusion(snow_z0s, snow_z0m)
  real, intent(out):: snow_z0s, snow_z0m

  snow_z0m =  z0_momentum
  snow_z0s =  z0_momentum * exp(-k_over_B)
end subroutine

end module snow_tile_mod
