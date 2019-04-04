subroutine prsys ( adiag, aleft, arite, f, nu )

!*****************************************************************************80
!
!! PRSYS prints out the tridiagonal linear system.
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
!    Input, real ( kind = 8 ) ADIAG(NU).
!    ADIAG(I) is the "diagonal" coefficient of the I-th
!    equation in the linear system.  That is, ADIAG(I) is
!    the coefficient of the I-th unknown in the I-th equation.
!
!    Input, real ( kind = 8 ) ALEFT(NU).
!    ALEFT(I) is the "left hand" coefficient of the I-th
!    equation in the linear system.  That is, ALEFT(I) is the
!    coefficient of the (I-1)-th unknown in the I-th equation.
!    There is no value in ALEFT(1), since the first equation
!    does not refer to a "0-th" unknown.
!
!    Input, real ( kind = 8 ) ARITE(NU).
!    ARITE(I) is the "right hand" coefficient of the I-th
!    equation in the linear system.  ARITE(I) is the coefficient
!    of the (I+1)-th unknown in the I-th equation.  There is
!    no value in ARITE(NU) because the NU-th equation does not
!    refer to an "NU+1"-th unknown.
!
!    Input, real ( kind = 8 ) F(NU).
!    ASSEMBLE stores into F the right hand side of the linear
!    equations.
!    SOLVE replaces those values of F by the solution of the
!    linear equations.
!
!    integer ( kind = 4 ) NU.
!    NU is the number of unknowns in the linear system.
!    Depending on the value of IBC, there will be N-1,
!    N, or N+1 unknown values, which are the coefficients
!    of basis functions.
!
  implicit none

  integer ( kind = 4 ) nu

  real ( kind = 8 ) adiag(nu)
  real ( kind = 8 ) aleft(nu)
  real ( kind = 8 ) arite(nu-1)
  real ( kind = 8 ) f(nu)
  integer ( kind = 4 ) i

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Printout of tridiagonal linear system:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Equation  A-Left  A-Diag  A-Rite  RHS'
  write ( *, '(a)' ) ' '

  do i = 1, nu
    if ( i == 1 ) then
      write(*,'(i3,14x,3g14.6)')i,adiag(i),arite(i),f(i)
    else if ( i < nu ) then
      write(*,'(i3,4g14.6)')i,aleft(i),adiag(i),arite(i),f(i)
    else
      write(*,'(i3,2g14.6,14x,g14.6)')i,aleft(i),adiag(i),f(i)
    end if
  end do

  return
end