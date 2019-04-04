subroutine get_beta ( beta )

!*****************************************************************************80
!
!! GET_BETA returns the value of BETA, for use by problem 5.
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
!    Output, real ( kind = 8 ) BETA, the value of BETA.
!
  implicit none

  real ( kind = 8 ) beta

  beta = -0.9D+00

  return
end