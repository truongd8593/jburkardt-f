subroutine get_problem ( problem )

!*****************************************************************************80
!
!! GETPRB returns the value of the current problem number.
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
!    Output, integer ( kind = 4 ) PROBLEM, the index of the problem.
!
  implicit none

  integer ( kind = 4 ) problem

  problem = 6

  return
end