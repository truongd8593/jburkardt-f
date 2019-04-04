subroutine solve ( adiag, aleft, arite, f, nu )

!*****************************************************************************80
!
!! SOLVE solves a tridiagonal matrix system of the form A*x = b.
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
!    Input/output, real ( kind = 8 ) ADIAG(NU), ALEFT(NU), ARITE(NU).
!    On input, ADIAG, ALEFT, and ARITE contain the diagonal,
!    left and right entries of the equations.
!    On output, ADIAG and ARITE have been changed in order
!    to compute the solution.
!    Note that for the first equation, there is no ALEFT
!    coefficient, and for the last, there is no ARITE.
!    So there is no need to store a value in ALEFT(1), nor
!    in ARITE(NU).
!
!    Input/output, real ( kind = 8 ) F(NU).
!    On input, F contains the right hand side of the linear
!    system to be solve
!    On output, F contains the solution of the linear system.
!
!    Input, integer ( kind = 4 ) NU, the number of equations to be solved.
!
  implicit none

  integer ( kind = 4 ) nu

  real ( kind = 8 ) adiag(nu)
  real ( kind = 8 ) aleft(nu)
  real ( kind = 8 ) arite(*)
  real ( kind = 8 ) f(nu)
  integer ( kind = 4 ) i
!
!  Handle the special case of a single equation.
!
  if ( nu == 1 ) then
    f(1) = f(1) / adiag(1)
!
!  The general case, when NU is greater than 1.
!
  else
    arite(1) = arite(1) / adiag(1)
    do i = 2, nu-1
      adiag(i) = adiag(i) - aleft(i) * arite(i-1)
      arite(i) = arite(i) / adiag(i)
    end do
    adiag(nu) = adiag(nu) - aleft(nu) * arite(nu-1)

    f(1) = f(1) / adiag(1)
    do i = 2, nu
      f(i) = ( f(i) - aleft(i) * f(i-1) ) / adiag(i)
    end do

    do i = nu-1, 1, -1
      f(i) = f(i) - arite(i) * f(i+1)
    end do

  end if

  return
end