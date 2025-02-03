module snow_opt_layers_mod

#include <fms_platform.h>
#include "../../shared/debug.inc"

use mpp_mod, only: input_nml_file
use fms_mod, only: error_mesg, check_nml_error, stdlog, mpp_pe, mpp_root_pe, lowercase, &
       string, FATAL, WARNING, NOTE

use land_data_mod,  only : log_version
use land_debug_mod, only : is_watch_point, land_error_message


implicit none
private

public :: dzopt_t
public :: snow_opt_layers_init

! ==== module constants ======================================================
character(len=*), parameter :: module_name = 'snowpack_mod' ! lm4p2
#include "../../shared/version_variable.inc"

! ---- types
type :: dzopt_t
    integer :: n ! number of elements of z, l
    real, ALLOCATABLE :: z(:) !< boundaries of the layers, m
    real, ALLOCATABLE :: l(:) !< layer number corresponding to layer boundaries
contains
    procedure :: init  => dzopt_init  !< set up optimal layer calculations, for given snow depth
    procedure :: clear => dzopt_clear !< free all allocated memory and reset to initial state
    procedure :: dz    => dzopt_dz    !< calculate optimal layer thickness at given depth, m
    procedure :: layer => dzopt_layer !< (fractional) layer number at the given depth, unitless
    procedure :: depth => dzopt_depth !< depth for given (fractional) layer number, m
    procedure :: penalty => dzopt_penalty !< calculate distance from the optimal vertical discretization
    procedure :: print => dzopt_print
end type dzopt_t


! ---- namelist
real :: opt_layer_top = 0.01 !< thickness of the top optimal layer, m
real :: opt_layer_bot = 0.03 !< thickness of the optimal layer at the bottom of the snowpack, m
real :: opt_layer_max = 1.0  !< maximum optimum layer thickness, m
real :: opt_layer_R   = 1.5  !< ratio of thicknesses of two adjacent layers in the middle of the snowpack, unitless
namelist /snow_opt_layers_nml/ &
        opt_layer_top, opt_layer_bot, opt_layer_max, opt_layer_R
! ---- end of namelist

contains

subroutine snow_opt_layers_init()
  integer :: io, ierr
  character(256) :: msg

  call log_version(version, module_name, &
  __FILE__)

  read (input_nml_file, snow_opt_layers_nml, iostat=io)
  ierr = check_nml_error(io, 'snow_opt_layers_nml :: '//trim(msg))
  if (mpp_pe() == mpp_root_pe()) then
     write(stdlog(), nml=snow_opt_layers_nml)
  endif

  ! check optimal layer parameters for sanity
  if (.not. opt_layer_top>0.0) call error_mesg('snow_opt_layers_init', &
       'opt_layer_top ='//string(opt_layer_top)//' in snow_opt_layers_nml in invalid: must be > 0.0', FATAL)
  if (.not. opt_layer_bot>0.0) call error_mesg('snow_opt_layers_init', &
       'opt_layer_bot ='//string(opt_layer_bot)//' in snow_opt_layers_nml in invalid: must be > 0.0', FATAL)
  if (.not. opt_layer_max>0.0) call error_mesg('snow_opt_layers_init', &
       'opt_layer_max ='//string(opt_layer_max)//' in snow_opt_layers_nml in invalid: must be > 0.0', FATAL)
  if (.not. opt_layer_R>=1.0) call error_mesg('snow_opt_layers_init', &
       'opt_layer_R ='//string(opt_layer_R)//' in snow_opt_layers_nml in invalid: must be >= 1.0'// &
       ' for the thickness of layers to increase with depth', FATAL)

end subroutine

!> initialize optimal layer thickness calculations
! Initialization sets up data arrays for calculation of z(layer) and its inverse layer(z),
! where "layer" is a real number, not an integer. these functions are used to calculate
! optimal thickness of layers for any given depth within the snowpack.
!
! TODO: the thickness of the lowest layer can be increased with total snowpack depth,
! since we do not expect it to matter when the snow is very deep (and therefore
! likely to be in equilibrium with underlying substrate)
subroutine dzopt_init(dzopt, depth)
  class(dzopt_t), intent(inout) :: dzopt
  real,           intent(in)    :: depth !< snow depth

  real    :: z, dz, d1, l1
  integer :: k, n

  ! Determine the number of layers needed for a given depth, so that the entire
  ! snowpack is covered by the array of optimal layers:
  z = 0.0; dz = opt_layer_top; n = 1
  do
     z  = z + dz
     n  = n + 1
     if (z > depth) exit ! from loop
     dz = min(dz*opt_layer_R, opt_layer_max)
  enddo
  ! z > depth, and n is at least 2
  dzopt%n = n+1 ! reserve space for one more layers in case a thin layer at the bottom
                ! needs to be inserted
  if (allocated(dzopt%z).or.allocated(dzopt%l)) &
       call land_error_message('dzopt_init :: arrays "z" and/or "l" are already allocated', FATAL)
  allocate(dzopt%z(dzopt%n)) ! boundaries of the layers, m
  allocate(dzopt%l(dzopt%n)) ! layer number corresponding to layer boundaries

  ! initialize depths of optimal layer tops and layer numbers. NOTE that calculation
  ! of layer boundaries must be exactly the same as in the above code where the number
  ! of layers is calculated
  dzopt%l(1) = 0.0; dzopt%z(1) = 0.0; dz = opt_layer_top;
  do k = 2,dzopt%n
     dzopt%l(k) = k-1
     dzopt%z(k) = dzopt%z(k-1) + dz
     dz         = min(dz*opt_layer_R, opt_layer_max)
  enddo

  d1 = max(depth-opt_layer_bot,0.0) ! depth to the near-soil layer
  k  = bisect(dzopt%z(:), d1)       ! k is between 1 and size(dzopt%z(:))-1
  dz = dzopt%z(k+1) - dzopt%z(k) ! thickness of optimal layer at the depth d1

  ! If necessary, add a thin layer at the bottom to better resolve gradients
  ! at the soil-snow interface
  if (dz>opt_layer_bot) then
      l1 = dzopt%layer(d1);
      dzopt%z(k+1) = d1     ; dzopt%l(k+1) = l1
      dzopt%z(k+2) = depth  ; dzopt%l(k+2) = l1+1
      dzopt%n = k+2
      endif

  ! Eliminate duplicate values that can arise when the boundary of the
  ! added bottom layer happens to be precisely equal to the boundary
  ! of the existing one.
  z = dzopt%z(1); n = 1
  do k = 2,dzopt%n
     if (dzopt%z(k) == z) then
        ! skip this array item
      else
        n = n+1
        dzopt%z(n) = dzopt%z(k)
        dzopt%l(n) = dzopt%l(k)
        z = dzopt%z(n)
      endif
  enddo
  dzopt%n = n
end subroutine dzopt_init

!>\brief  Given an 1D array xx of interval bounds in ascending order, and a value x,
!!  finds a position of point x in array xx. Returns i, such that x is
!!  between xx(i) and xx(i+1); for x below the lowest bound in xx returns 1;
!!  for x above the upper bound returns size(xx)-1
real function bisect(xx, x)
  real, intent(in)              :: xx(:)     !< array of boundaries
  real, intent(in)              :: x         !< point to locate

  ! ---- local vars
  integer :: low, high, mid
  integer :: n              ! size of the input array

  n = size(xx)

  ! find the coordinates
  if (x >= xx(1).and.x<=xx(n)) then
     low = 1; high = n
     do while (high-low > 1)
        mid = (low+high)/2
        if (xx(mid) <= x) then
           low = mid
        else
           high = mid
        endif
     enddo
     bisect = low
  else if (x < xx(1)) then
     bisect = 1
  else if (x > xx(n)) then
     bisect = n-1
  endif
end function bisect

!\> given depth within snowpack, calculate corresponding (fractional) optimal layer number
real function dzopt_layer(dzopt, d) result (layer)
  class(dzopt_t), intent(in) :: dzopt
  real,           intent(in) :: d

  integer :: i
  ! TODO: check that dzopt is initialized
  associate(z=>dzopt%z, l=>dzopt%l, n=>dzopt%n)
     i = bisect(dzopt%z(1:n),d)
     layer = l(i) + (d-z(i)) * (l(i+1)-l(i))/(z(i+1)-z(i))
  end associate
end function dzopt_layer

!> given layer (perhaps fractional), calculate corresponding depth
real function dzopt_depth(dzopt, layer) result(depth)
  class(dzopt_t), intent(in) :: dzopt
  real,           intent(in) :: layer !< layer number, perhaps fractional

  integer :: i
  ! TODO: check that dzopt is initialized
  associate(z=>dzopt%z, l=>dzopt%l, n=>dzopt%n)
     i = bisect(dzopt%l(1:n),layer)
     depth = z(i) + (layer-l(i)) * (z(i+1)-z(i))/(l(i+1)-l(i))
  end associate
end function dzopt_depth

!> calculate optimal thickness of layer with the top at z
real function dzopt_dz(dzopt,z) result(dz)
  class(dzopt_t), intent(in) :: dzopt
  real,           intent(in) :: z

  dz = dzopt%depth(dzopt%layer(z)+1.0) - z
end function dzopt_dz

!> \brief calculate penalty value for given layer distribution
real function dzopt_penalty(dzopt, z) result(p)
  class(dzopt_t), intent(in) :: dzopt
  real,           intent(in) :: z(:) !< boundaries of the snow layers

  integer :: k
  real :: dz_opt, dz

  p = 0
  do k = 1, size(z)-1
     dz_opt = dzopt%dz(z(k))
     dz     = z(k+1)-z(k)
     p = p + (dz-dz_opt)**2
  enddo
end function dzopt_penalty

!> \brief print optimal vertical discretization
subroutine dzopt_print(dzopt)
  class(dzopt_t), intent(in) :: dzopt

  integer :: k
  write(*,'(a2,99(",",a9,:))') "k", "z","l","dz_opt"
  do k = 1, dzopt%n
     write(*,'(i2.2,99(",",f9.4,:))') k, dzopt%z(k), dzopt%l(k), dzopt%dz(dzopt%z(k))
  enddo
end subroutine dzopt_print

subroutine dzopt_clear(dzopt)
  class(dzopt_t), intent(inout) :: dzopt

  if (.not.allocated(dzopt%z)) &
        call land_error_message('dzopt%clear :: attempt to deallocate array z that was not allocated', FATAL)
  if (.not.allocated(dzopt%l)) &
        call land_error_message('dzopt%clear :: attempt to deallocate array l that was not allocated', FATAL)
  deallocate(dzopt%z, dzopt%l)
  dzopt%n = 0
end subroutine dzopt_clear

end module snow_opt_layers_mod