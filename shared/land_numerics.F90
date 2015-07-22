! ============================================================================
! module numerics: a collection of useful general-purpose routines
! ============================================================================
#define __ERROR__(message) \
call my_error(mod_name,message,FATAL,__FILE__,__LINE__)
#define __NOTE__(message) \
call my_error(mod_name,message,NOTE,__FILE__,__LINE__)
#define __ASSERT__(x, message) \
if(.NOT.(x))call my_error(mod_name,message,FATAL,__FILE__,__LINE__)

module numerics_mod

use fms_mod, only: error_mesg, FATAL, NOTE, write_version_number

implicit none
private

! ==== public interfaces =====================================================
public :: bisect    ! finds a position of point in array of bounds
public :: lin_int   ! linear interpolation
public :: ludcmp, lubksb ! LU decomposition and back substitution
public :: tridiag   ! tridiagonal system solver
public :: numerics_init
! ==== end of public interfaces ==============================================


!     Linear interpolation.
interface lin_int
   module procedure lin_int0
   module procedure lin_int1
   module procedure lin_int2
   module procedure lin_int1m
   module procedure lin_int2m
end interface

logical :: module_is_initialized =.FALSE.
! module constants
character(len=*), parameter :: &
     mod_name = 'numerics_mod', &
     version  = '$Id$', &
     tagname  = '$Name$'

contains


! ============================================================================
! Initializes the numerics module.
subroutine numerics_init()

  module_is_initialized =.TRUE. 
  call write_version_number(version,tagname)

end subroutine numerics_init


! ============================================================================
!  Finds a position of point in array of bounds. Returns i, such that x is
!  between xx(i) and xx(i+1).
!  Usage:
!     value=bisect( xx, x1, periodic )
function bisect(xx, x1, periodic)
  real, intent(in)              :: xx(:)     ! array of boundaries
  real, intent(in)              :: x1        ! point to locate
  logical, intent(in), optional :: periodic  ! if present and true, the data
                                             ! domain is assumed to be periodic
  ! ---- result type ---------------------------------------------------
  integer bisect

  ! ---- local vars ----------------------------------------------------
  real    :: x              ! duplicate of input value
  integer :: low, high, mid
  integer :: n              ! size of the input array
  logical :: ascending      ! if true, the coordinates are in ascending order

  n = size(xx)
  x = x1

  ! bring the point inside bounds of the period
  if (present(periodic).AND.periodic) then
     __ASSERT__(xx(n)-xx(1)/=0,"periodic bisect: period equal to zero")
     x = modulo(x-min(xx(1),xx(n)),abs(xx(n)-xx(1)))+min(xx(1),xx(n))
  endif

  ! find the coordinates
  if (x >= xx(1).and.x<=xx(n)) then
     low = 1; high = n
     ascending = xx(n) > xx(1)
     do while (high-low > 1)
        mid = (low+high)/2
        if (ascending.eqv.xx(mid) <= x) then
           low = mid
        else
           high = mid
        endif
     enddo
     bisect = low
  else
     bisect = -1
  endif

end function bisect


! ==============================================================================
!     Linearly interpolates 1-D data.
subroutine lin_int0(data, xx, x, res)
  real, intent(in) :: data(:)    ! data to interpolate
  real, intent(in) :: xx(:)      ! coord. corresponding to data
  real, intent(in) :: x          ! coord to interpolate to
  real, intent(inout) :: res     ! result of interpolation

  ! ---- local vars ----------------------------------------------------------
  integer :: i1, i2
  real    :: f1, f2

  ! find where is our time point and calculate weights
  i1 = bisect(xx,x)
  i2 = i1+1
  __ASSERT__(i1>0.AND.i1<size(xx),"Coordinate is out of range")

  f2 = (x-xx(i1))/(xx(i2)-xx(i1))
  f1 = 1.0-f2

  ! update the result
  res = data(i1)*f1+data(i2)*f2

end subroutine lin_int0


! ==============================================================================
!     Linearly interpolates 1-D data.
subroutine lin_int1(data, xx, x, res)

  real, intent(in) :: data(:,:)    ! data to interpolate
  real, intent(in) :: xx(:)        ! coord. corresponding to data
  real, intent(in) :: x            ! coord to interpolate to
  real, intent(inout) :: res(:)    ! result of interpolation

  ! ---- local vars ----------------------------------------------------------
  integer :: i1, i2
  real    :: f1, f2

  ! find where is our time point and calculate weights
  i1 = bisect(xx,x)
  i2 = i1+1

  __ASSERT__(i1>0.AND.i1<size(xx),"Coordinate is out of range")

  f2 = (x-xx(i1))/(xx(i2)-xx(i1))
  f1 = 1.0-f2

  ! update the result
  res = data(:,i1)*f1+data(:,i2)*f2

end subroutine lin_int1


! ==============================================================================
!     Interpolates prescribed over time.
subroutine lin_int2(data, tt, t, res)

  real, intent(in) :: data(:,:,:)  ! data to interpolate
  real, intent(in) :: tt(:)        ! time moments corresponding to data points
  real, intent(in) :: t            ! time to interpolate to
  real, intent(inout) :: res(:,:)  ! result

  ! ---- local vars ----------------------------------------------------------
  integer :: i1, i2
  real    :: f1, f2

  ! find where is our time point and calculate weights
  i1 = bisect(tt,t); i2 = i1+1
  __ASSERT__(i1>0.AND.i1<size(tt),"Coordinate is out of range")

  f2 = (t-tt(i1))/(tt(i2)-tt(i1))
  f1 = 1-f2

  ! update the result
  res = data(:,:,i1)*f1+data(:,:,i2)*f2

end subroutine lin_int2


!     Linearly interpolates 1-D data.
subroutine lin_int1m(data, xx, x, res, mask)
  real, intent(in) :: data(:,:)    ! data to interpolate
  real, intent(in) :: xx(:)        ! coord. corresponding to data
  real, intent(in) :: x            ! coord to interpolate to
  real, intent(inout) :: res(:)    ! result of interpolation
  logical, intent(in) :: mask(:)   ! valid data mask

  ! ---- local vars ----------------------------------------------------------
  integer :: i1, i2
  real    :: f1, f2

  ! find where is our time point and calculate weights
  i1 = bisect(xx,x)
  i2 = i1+1
  __ASSERT__(i1>0.AND.i1<size(xx),"Coordinate is out of range")

  f2 = (x-xx(i1))/(xx(i2)-xx(i1))
  f1 = 1.0-f2

  ! finally, update the result
  where (mask) 
     res = data(:,i1)*f1+data(:,i2)*f2
  endwhere

end subroutine lin_int1m


! ==============================================================================
!     Interpolates prescribed over time.
subroutine lin_int2m(data, tt, t, res, mask)
  real, intent(in) :: data(:,:,:)  ! data to interpolate
  real, intent(in) :: tt(:)        ! time moments corresponding to data points
  real, intent(in) :: t            ! time to interpolate to
  real, intent(inout) :: res(:,:)  ! result
  logical, intent(in) :: mask(:,:) ! interpolation mask

  ! ---- local vars ----------------------------------------------------------
  integer :: i1, i2
  real    :: f1, f2

  ! find where is our time point and calculate weights
  i1 = bisect(tt,t); i2 = i1+1
  __ASSERT__(i1>0.AND.i1<size(tt),"Coordinate is out of range")

  f2 = (t-tt(i1))/(tt(i2)-tt(i1))
  f1 = 1-f2

  ! update the result
  where (mask) 
     res = data(:,:,i1)*f1+data(:,:,i2)*f2
  endwhere
end subroutine lin_int2m


! ==============================================================================
! given a matrix a(n,n) replaces it by LU decomposition of a rowwise permutation
! of itself indx(n) is an output that records the permuatation. This routine is
! used in combination with lubksb to solve linear equations or invert a matrix
! example:
!    call ludcmp(a,indx)
!    call lubksb(a,indx,b1)
!    call lubksb(a,indx,b2)
subroutine ludcmp(a,indx,status)
  real,    intent(inout) :: a(:,:) ! matrix that gets replaced by its LU decomposition
  integer, intent(out)   :: indx(:) ! index of row permutations effecte by partial pivoting
  integer, intent(out), optional :: status

  integer, parameter :: TINY = 1.0e-20
  integer :: n ! size of the matrix 
  integer :: i,j,k,imax
  real    :: aamax,dum,sum
  real    :: vv(size(a,1)) ! implicit scaling for each row

  n = size(a,1)
  if(present(status))status = 0

  ! find largest element in each row and calculate scaling 
  do i = 1,n
     aamax = 0.0
     do j = 1,n
        if(abs(a(i,j)) > aamax)aamax = abs(a(i,j))
     enddo
     if(.not.(aamax /= 0.0)) then 
        if(present(status))then
           status = -1; aamax = TINY
        else
           call error_mesg('ludcmp','Matrix is singular', FATAL)
        endif
     endif
     vv(i) = 1.0/aamax
  enddo

  ! loop over the columns of Crout's method
  do j=1,n
     do i = 1,j-1
        sum = a(i,j)
        do k = 1,i-1
           sum = sum-a(i,k)*a(k,j)
        enddo
        a(i,j) = sum
     enddo
     aamax = 0.0 ! initialize the search for the largest pivot element
     do i = j,n
        sum = a(i,j)
        do k = 1,j-1
           sum = sum-a(i,k)*a(k,j)
        enddo
        a(i,j) = sum
        dum = vv(i)*abs(sum) ! figure of merit for the pivot
        if (dum >= aamax) then ! is it better than the best so far?
           imax = i
           aamax = dum
        endif
     enddo
     if (j /= imax) then ! do we need to interchange rows?
        ! Yes, do so
        do k=1,n
           dum=a(imax,k)
           a(imax,k) = a(j,k)
           a(j,k) = dum
        enddo
        vv(imax) = vv(j)
     endif
     indx(j) = imax
     ! if the pivot element is zero, then the matrix is singular (at least to the
     ! precision of the algorithm). For some applications on singular matrices, it 
     ! is desirable to substitute TINY for zero
     if(a(j,j)==0.0) a(j,j) = TINY

     if (j/=n)then
        ! Finally, divide by the pivot element
        dum = 1.0/a(j,j)
        do i = j+1,n
           a(i,j) = a(i,j)*dum
        enddo
     endif
  enddo ! loop over the columns
end subroutine ludcmp


! ==============================================================================
! given a LU decomposition of matrix a(n,n), permuation vector indx, and right-
! hand side b, solves the set of linear equations A*X = B 
subroutine lubksb(a,indx,b)
  real,    intent(in)    :: a(:,:)  ! LU-decomposed matrix
  integer, intent(in)    :: indx(:) ! permutation vector, as returned by the ludcmp
  real,    intent(inout) :: b(:)    ! right-hand side vector, returns solution

  integer :: i,ii,j,ll,n
  real    :: sum

  n = size(a,1)
  ii = 0
  do i = 1,n
     ll = indx(i)
     sum = b(ll)
     b(ll) = b(i)
     if (ii/=0) then
        do j = ii,i-1
           sum = sum - a(i,j)*b(j)
        enddo
     else if (sum /= 0.0) then
        ii = i
     endif
     b(i) = sum
  enddo
  do i = n, 1, -1
     sum = b(i)
     do j = i+1,n
        sum = sum-a(i,j)*b(j)
     enddo
     b(i) = sum/a(i,i)
  enddo
end subroutine lubksb


! ============================================================================
! given values of the triadiagonal matrix coefficients, computes a solution
subroutine tridiag(a,b,c,r,u)
  real, intent(in)  :: a(:),b(:),c(:),r(:)
  real, intent(out) :: u(:)

  integer :: j
  real :: bet, gam(size(a))
  
  ! check that the sizes are the same
  if(size(a)/=size(b).or.size(a)/=size(c).or.size(a)/=size(r)) &
       call error_mesg('tridiag','sizes of input arrays are not equal',FATAL)
  if(size(u)<size(a)) &
       call error_mesg('tridiag','size of the result is insufficient',FATAL)
  ! check that a(1)==0 and c(N)==0
  if(a(1)/=0.or.c(size(a))/=0) &
       call error_mesg('tridiag','a(1) and c(N) must be equal to 0',FATAL)
  ! decomposition and forward substitution
  bet = b(1)
  u(1) = r(1)/bet
  do j = 2,size(a)
     gam(j) = c(j-1)/bet
     bet = b(j)-a(j)*gam(j)
     if(bet==0) &
          call error_mesg('tridiag','system is ill-defined',FATAL)
     u(j) = (r(j)-a(j)*u(j-1))/bet
  enddo
  ! backward substitution
  do j = size(a)-1,1,-1
     u(j) = u(j)-gam(j+1)*u(j+1)
  enddo
end subroutine tridiag


! ==============================================================================
! Reports error, including file name and line.
subroutine my_error(mod_name, message, mode, file, line)

  character(len=*), intent(in) :: mod_name
  character(len=*), intent(in) :: message
  integer,          intent(in) :: mode
  character(len=*), intent(in) :: file
  integer,          intent(in) :: line

  ! ---- local vars ----------------------------------------------------------
  character(len=512) :: mesg
  write(mesg,'("File ",a," Line ",i4.4," :: ",a)')&
       file, line, trim(message)
  call error_mesg(mod_name, mesg, mode)
end subroutine


end module numerics_mod
