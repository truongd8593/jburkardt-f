function pp ( x )

!*****************************************************************************80
!
!! PP evaluates the function P in the differential equation.
!
!  Discussion:
!
!    The function P(X) occurs in the differential equation:
!
!      -d/dx ( P(x) * dU(x)/dx ) + Q(x) * U(x)  =  F(x)
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
!    Input, real ( kind = 8 ) X, the evaluation point.
!
!    Output, real ( kind = 8 ) PP, the value of P(X).
!
  implicit none

  integer ( kind = 4 ) problem
  real ( kind = 8 ) pp
  real ( kind = 8 ) x
!
!  Find out which problem we're working on.
!
  call get_problem ( problem )

  if ( problem == 1 ) then
    pp = 1.0D+00
  else if ( problem == 2 ) then
    pp = 1.0D+00
  else if ( problem == 3 ) then
    pp = 1.0D+00
  else if ( problem == 4 ) then
    pp = 1.0D+00
  else if ( problem == 5 ) then
    pp = 1.0D+00
  else if ( problem == 6 ) then
    pp = 1.0D+00
  end if

  return
end