function ff ( x )

!*****************************************************************************80
!
!! FF evaluates the function F in the differential equation.
!
!  Discussion:
!
!    This is the function F(X) that appears on the right hand
!    side of the equation:
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
!    Output, real ( kind = 8 ) FF, the value of F(X).
!
  implicit none

  real ( kind = 8 ) alpha
  real ( kind = 8 ) beta
  real ( kind = 8 ) ff
  integer ( kind = 4 ) problem
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) x
!
!  Find out which problem we're working on.
!
  call get_problem ( problem )

  if ( problem == 1 ) then

    ff = 0.0D+00

  else if ( problem == 2 ) then

    ff = -2.0D+00 * x

  else if ( problem == 3 ) then

    ff = 0.25D+00 * pi**2 * sin ( 0.5D+00 * pi * x )

  else if ( problem == 4 ) then

    ff = 0.25D+00 * pi**2 * cos ( 0.5D+00 * pi * x )

  else if ( problem == 5 ) then

    call get_beta ( beta )

    ff = - ( x**beta ) + ( x**( beta + 2.0D+00 ) ) &
      / ( ( beta + 2.0D+00 ) * ( beta + 1.0D+00 ) )

  else if ( problem == 6 ) then

    call get_alpha ( alpha )
    ff = 2.0D+00 * alpha * ( x - 0.5D+00 ) &
      / ( alpha**2 + ( x - 0.5D+00 )**2 )**2

  end if

  return
end