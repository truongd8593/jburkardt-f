function u_exact ( x )

!*****************************************************************************80
!
!! U_EXACT returns the value of the exact solution at any point X.
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
!    Output, real ( kind = 8 ) U_EXACT, the value of the exact solution at X.
!
  implicit none

  real ( kind = 8 ) alpha
  real ( kind = 8 ) beta
  integer ( kind = 4 ) problem
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) u_exact
  real ( kind = 8 ) x
!
!  Find out which problem we're working on.
!
  call get_problem ( problem )

  if ( problem == 1 ) then
    u_exact = x
  else if ( problem == 2 ) then
    u_exact = x**2
  else if ( problem == 3 ) then
    u_exact = sin ( pi * x / 2.0D+00 )
  else if ( problem == 4 ) then
    u_exact = cos ( pi * x / 2.0D+00 )
  else if ( problem == 5 ) then
    call get_beta ( beta )
    u_exact = ( x**( beta + 2.0D+00 ) ) &
      / ( ( beta + 2.0D+00 ) * ( beta + 1.0D+00 ) )
  else if ( problem == 6 ) then
    call get_alpha ( alpha )
    u_exact = atan ( ( x - 0.5D+00 ) / alpha )
  else
    u_exact = 0.0D+00
  end if

  return
end