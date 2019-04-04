subroutine solvex ( f, h, ibc, kount, n, nl, nmax, nu, ul, ur, xn )

!*****************************************************************************80
!
!! SOLVEX discretizes and solves a differential equation given the nodes.
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
!    Output, real ( kind = 8 ) F(NU).
!    ASSEMBLE stores into F the right hand side of the linear
!    equations.
!    SOLVE replaces those values of F by the solution of the
!    linear equations.
!
!    Output, real ( kind = 8 ) H(N)
!    H(I) is the length of subinterval I.  This code uses
!    equal spacing for all the subintervals.
!
!    Input, integer ( kind = 4 ) IBC.
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
!    Input, integer ( kind = 4 ) KOUNT, the number of adaptive steps that 
!    have been taken.
!
!    Input, integer ( kind = 4 ) N
!    The number of subintervals into which the interval
!    [XL,XR] is broken.
!
!    Input, integer ( kind = 4 ) NL.
!    The number of basis functions used in a single
!    subinterval.  (NL-1) is the degree of the polynomials
!    used.  For this code, NL is fixed at 2, meaning that
!    piecewise linear functions are used as the basis.
!
!    Input, integer ( kind = 4 ) NMAX, the maximum number of unknowns that 
!    can be handled.
!
!    Output, integer ( kind = 4 ) NU.
!    NU is the number of unknowns in the linear system.
!    Depending on the value of IBC, there will be N-1,
!    N, or N+1 unknown values, which are the coefficients
!    of basis functions.
!
!    Input, real ( kind = 8 ) UL.
!    If IBC is 1 or 3, UL is the value that U is required
!    to have at X = XL.
!    If IBC is 2 or 4, UL is the value that U' is required
!    to have at X = XL.
!
!    Input, real ( kind = 8 ) UR.
!    If IBC is 2 or 3, UR is the value that U is required
!    to have at X = XR.
!    If IBC is 1 or 4, UR is the value that U' is required
!    to have at X = XR.
!
!    Input, real ( kind = 8 ) XN(0:N).
!    XN(I) is the location of the I-th node.  XN(0) is XL,
!    and XN(N) is XR.
!
!  Local parameters:
!
!    Local, real ( kind = 8 ) ADIAG(NU).
!    ADIAG(I) is the "diagonal" coefficient of the I-th
!    equation in the linear system.  That is, ADIAG(I) is
!    the coefficient of the I-th unknown in the I-th equation.
!
!    Local, real ( kind = 8 ) ALEFT(NU).
!    ALEFT(I) is the "left hand" coefficient of the I-th
!    equation in the linear system.  That is, ALEFT(I) is the
!    coefficient of the (I-1)-th unknown in the I-th equation.
!    There is no value in ALEFT(1), since the first equation
!    does not refer to a "0-th" unknown.
!
!    Local, real ( kind = 8 ) ARITE(NU).
!    ARITE(I) is the "right hand" coefficient of the I-th
!    equation in the linear system.  ARITE(I) is the coefficient
!    of the (I+1)-th unknown in the I-th equation.  There is
!    no value in ARITE(NU) because the NU-th equation does not
!    refer to an "NU+1"-th unknown.
!
!    Local, integer ( kind = 4 ) INDX(0:N).
!    For a node I, INDX(I) is the index of the unknown
!    associated with node I.
!    If INDX(I) is equal to -1, then no unknown is associated
!    with the node, because a boundary condition fixing the
!    value of U has been applied at the node instead.
!    Unknowns are numbered beginning with 1.
!    If IBC is 2 or 4, then there is an unknown value of U
!    at node 0, which will be unknown number 1.  Otherwise,
!    unknown number 1 will be associated with node 1.
!    If IBC is 1 or 4, then there is an unknown value of U
!    at node N, which will be unknown N or N+1,
!    depending on whether there was an unknown at node 0.
!
!    Local, integer ( kind = 4 ) NODE(NL,N).
!    For each subinterval I:
!    NODE(1,I) is the number of the left node, and
!    NODE(2,I) is the number of the right node.
!
!    Local, integer ( kind = 4 ) NQUAD
!    The number of quadrature points used in a subinterval.
!
!    Local, real ( kind = 8 ) WQUAD(NQUAD).
!    WQUAD(I) is the weight associated with the I-th point
!    of an NQUAD point Gaussian quadrature rule.
!
!    Local, real ( kind = 8 ) XQUAD(NQUAD,NMAX), the I-th quadrature point
!    in interval J.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nl
  integer ( kind = 4 ) nmax
  integer ( kind = 4 ), parameter :: nquad = 2

  real ( kind = 8 ) adiag(nmax)
  real ( kind = 8 ) aleft(nmax)
  real ( kind = 8 ) arite(nmax)
  real ( kind = 8 ) f(nmax)
  real ( kind = 8 ) h(n)
  integer ( kind = 4 ) ibc
  integer ( kind = 4 ) indx(0:n)
  integer ( kind = 4 ) kount
  integer ( kind = 4 ) node(nl,nmax)
  integer ( kind = 4 ) nu
  real ( kind = 8 ) ul
  real ( kind = 8 ) ur
  real ( kind = 8 ) wquad(nquad)
  real ( kind = 8 ) xn(0:n)
  real ( kind = 8 ) xquad(nquad,nmax)
!
!  Given a set of N nodes (where N increases on each iteration),
!  compute the other geometric information.
!
  call geometry ( h, ibc, indx, n, nl, nmax, node, nquad, nu, wquad, xn, xquad )
!
!  Assemble the linear system.
!
  call assemble ( adiag, aleft, arite, f, h, n, indx, node, nu, nl, &
    nquad, nmax, ul, ur, wquad, xn, xquad )
!
!  Print out the linear system, just once.
!
  if ( kount == 1 ) then
    call prsys ( adiag, aleft, arite, f, nu )
  end if
!
!  Solve the linear system.
!
  call solve ( adiag, aleft, arite, f, nu )
!
!  Print out the solution.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Basic solution'

  !!! call output ( f, ibc, indx, n, nu, ul, ur, xn )
  call write_solution_tecplot ( f, ibc, indx, n, nu, ul, ur, xn )

  return
end