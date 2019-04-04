subroutine get_alpha ( alpha )

!*****************************************************************************80
!
!! GET_ALPHA returns the value of ALPHA, for use by problem 6.
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
!    Output, real ( kind = 8 ) ALPHA, the value of ALPHA.
!
  implicit none

  real ( kind = 8 ) alpha

  alpha = 0.01D+00

  return
end