subroutine phi ( il, x, phii, phiix, xleft, xrite )

!*****************************************************************************80
!
!! PHI evaluates a linear basis function and its derivative.
!
!  Discussion:
!
!    The functions are evaluated at a point X in an interval.  In any
!    interval, there are just two basis functions.  The first
!    basis function is a line which is 1 at the left endpoint
!    and 0 at the right.  The second basis function is 0 at
!    the left endpoint and 1 at the right.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 November 2006
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IL, the local index of the basis function.
!
!    Input, real ( kind = 8 ) X, the evaluation point.
!
!    Output, real ( kind = 8 ) PHII, PHIIX, the value of the basis function
!    and its derivative.
!
!    Input, real ( kind = 8 ) XLEFT, XRITE, the endpoints of the interval.
!
  implicit none

  integer ( kind = 4 ) il
  real ( kind = 8 ) phii
  real ( kind = 8 ) phiix
  real ( kind = 8 ) x
  real ( kind = 8 ) xleft
  real ( kind = 8 ) xrite

  if ( xleft <= x .and. x <= xrite ) then

    if ( il == 1 ) then
      phii = ( xrite - x ) / ( xrite - xleft )
      phiix = -1.0D+00 / ( xrite - xleft )
    else
      phii = ( x - xleft ) / ( xrite - xleft )
      phiix = 1.0D+00 / ( xrite - xleft )
    end if
!
!  If X is outside of the interval, then the basis function
!  is always zero.
!
  else
    phii = 0.0D+00
    phiix = 0.0D+00
  end if

  return
end