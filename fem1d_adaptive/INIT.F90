subroutine init ( ibc, n, tol, ul, ur, xl, xn, xr )

!*****************************************************************************80
!
!! INIT initializes some parameters that define the problem.
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
!    Output, integer ( kind = 4 ) IBC.
!    IBC declares what the boundary conditions are.
!    1, at the left endpoint, U has the value UL,
!       at the right endpoint, U' has the value UR.
!    2, at the left endpoint, U' has the value UL,
!       at the right endpoint, U has the value UR.
!    3, at the left endpoint, U has the value UL,
!       and at the right endpoint, U has the value UR.
!    4, at the left endpoint, U' has the value UL,
!       at the right endpoint U' has the value UR.
!
!    Input, integer ( kind = 4 ) N
!    The number of subintervals into which the interval
!    [XL,XR] is broken.
!
!    Output, real ( kind = 8 ) TOL.
!    A tolerance that is used to determine whether the estimated
!    error in an interval is so large that it should be subdivided
!    and the problem solved again.
!
!    Output, real ( kind = 8 ) UL.
!    If IBC is 1 or 3, UL is the value that U is required
!    to have at X = XL.
!    If IBC is 2 or 4, UL is the value that U' is required
!    to have at X = XL.
!
!    Output, real ( kind = 8 ) UR.
!    If IBC is 2 or 3, UR is the value that U is required
!    to have at X = XR.
!    If IBC is 1 or 4, UR is the value that U' is required
!    to have at X = XR.
!
!    Output, real ( kind = 8 ) XL.
!    XL is the left endpoint of the interval over which the
!    differential equation is being solved.
!
!    Output, real ( kind = 8 ) XN(0:N).
!    XN(I) is the location of the I-th node.  XN(0) is XL,
!    and XN(N) is XR.
!
!    Output, real ( kind = 8 ) XR.
!    XR is the right endpoint of the interval over which the
!    differential equation is being solved.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) alpha
  real ( kind = 8 ) beta
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ibc
  integer ( kind = 4 ) problem
  real ( kind = 8 ) tol
  real ( kind = 8 ) u_exact
  real ( kind = 8 ) ul
  real ( kind = 8 ) ur
  real ( kind = 8 ) xl
  real ( kind = 8 ) xn(0:n)
  real ( kind = 8 ) xr

  tol = 0.01D+00
!
!  Find out which problem we're working on.
!
  call get_problem ( problem )
!
!  Set the boundary conditions for the problem, and
!  print out its title.
!
  if ( problem == 1 ) then

    ibc = 3
    ul = 0.0D+00
    ur = 1.0D+00
    xl = 0.0D+00
    xr = 1.0D+00
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Exact solution is U = X'

  else if ( problem == 2 ) then

    ibc = 3
    ul = 0.0D+00
    ur = 1.0D+00
    xl = 0.0D+00
    xr = 1.0D+00
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Exact solution is U = X*X'

  else if ( problem == 3 ) then

    ibc = 3
    ul = 0.0D+00
    ur = 1.0D+00
    xl = 0.0D+00
    xr = 1.0D+00
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Exact solution is U = SIN(PI*X/2)'

  else if ( problem == 4 ) then

    ibc = 3
    ul = 1.0D+00
    ur = 0.0D+00
    xl = 0.0D+00
    xr = 1.0D+00
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Exact solution is U = COS(PI*X/2)'

  else if ( problem == 5 ) then

    ibc = 3
    call get_beta ( beta )
    ul = 0.0D+00
    ur = 1.0D+00 / ( ( beta + 2.0D+00 ) * ( beta + 1.0D+00 ) )
    xl = 0.0D+00
    xr = 1.0D+00
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Rheinboldt problem'

  else if ( problem == 6 ) then

    ibc = 3
    call get_alpha ( alpha )
    xl = 0.0D+00
    xr = 1.0D+00
    ul = u_exact ( xl )
    ur = u_exact ( xr )
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Arctangent problem'

  end if
!
!  The nodes are defined here, and not in the geometry routine.
!  This is because each new iteration chooses the location
!  of the new nodes in a special way.
!
  do i = 0, n
    xn(i) = ( real ( n - i, kind = 8 ) * xl   &
            + real (     i, kind = 8 ) * xr ) &
            / real ( n,     kind = 8 )
  end do

  write ( *, '(a)' ) 'The equation is to be solved for '
  write ( *, '(a,g14.6)' ) 'X greater than ', xl
  write ( *, '(a,g14.6)' ) ' and less than ', xr
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'The boundary conditions are:'
  write ( *, '(a)' ) ' '

  if ( ibc == 1 .or. ibc == 3 ) then
    write ( *, '(a,g14.6)' ) '  At X = XL, U=', ul
  else
    write ( *, '(a,g14.6)' ) '  At X = XL, U''=', ul
  end if

  if ( ibc == 2 .or. ibc == 3 ) then
    write ( *, '(a,g14.6)' ) '  At X = XR, U=', ur
  else
    write ( *, '(a,g14.6)' ) '  At X = XR, U''=', ur
  end if

  return
end