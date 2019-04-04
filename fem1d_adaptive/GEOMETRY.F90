subroutine geometry ( h, ibc, indx, n, nl, nmax, node, nquad, nu, wquad, xn, &
  xquad )

!*****************************************************************************80
!
!! GEOMETRY sets up some of the geometric information for the problem.
!
!  Discussion:
!
!    Note, however, that the location of the nodes
!    is done outside of this routine, and, in fact, before this
!    routine is called.
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
!    Output, integer ( kind = 4 ) INDX(0:N).
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
!    Output, integer ( kind = 4 ) NODE(NL,N).
!    For each subinterval I:
!    NODE(1,I) is the number of the left node, and
!    NODE(2,I) is the number of the right node.
!
!    Input, integer ( kind = 4 ) NQUAD
!    The number of quadrature points used in a subinterval.
!
!    Output, integer ( kind = 4 ) NU.
!    NU is the number of unknowns in the linear system.
!    Depending on the value of IBC, there will be N-1,
!    N, or N+1 unknown values, which are the coefficients
!    of basis functions.
!
!    Output, real ( kind = 8 ) WQUAD(NQUAD).
!    WQUAD(I) is the weight associated with the I-th point
!    of an NQUAD point Gaussian quadrature rule.
!
!    Input, real ( kind = 8 ) XN(0:N).
!    XN(I) is the location of the I-th node.  XN(0) is XL,
!    and XN(N) is XR.
!
!    Output, real ( kind = 8 ) XQUAD(NQUAD,NMAX), the I-th quadrature point
!    in interval J.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nl
  integer ( kind = 4 ) nmax
  integer ( kind = 4 ) nquad

  real ( kind = 8 ) alfa
  real ( kind = 8 ) h(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ibc
  integer ( kind = 4 ) igl
  integer ( kind = 4 ) igr
  integer ( kind = 4 ) indx(0:n)
  integer ( kind = 4 ) node(nl,nmax)
  integer ( kind = 4 ) nu
  real ( kind = 8 ) wquad(nquad)
  real ( kind = 8 ) xl
  real ( kind = 8 ) xn(0:n)
  real ( kind = 8 ) xquad(nquad,nmax)
  real ( kind = 8 ) xr
!
!  Store in NODE the fact that interval I has node I-1
!  as its left endpoint, and node I as its right endpoint.
!
  do i = 1, n
    node(1,i) = i-1
    node(2,i) = i
  end do
!
!  For every node that is associated with an unknown, we
!  record the number of the unknown in INDX.
!
  nu = 0
  do i = 0, n

    if ( i == 0 .and. ( ibc == 1 .or. ibc == 3 ) ) then
      indx(i) = -1
    else if ( i == n .and. ( ibc == 2 .or. ibc == 3 ) ) then
      indx(i) = -1
    else
      nu = nu+1
      indx(i) = nu
    end if

  end do
!
!  We compute the width of each interval.
!
  do i = 1, n
    igl = node(1,i)
    igr = node(2,i)
    h(i) = xn(igr) - xn(igl)
  end do
!
!  We compute the location of the quadrature points in each
!  interval.
!
  do i = 1, n

    xl = xn(node(1,i))
    xr = xn(node(2,i))

    if ( nquad == 1 ) then
      xquad(1,i) = 0.5D+00 * ( xl + xr )
    else if ( nquad == 2 ) then
      alfa = -0.577350D+00
      xquad(1,i) = ( ( 1.0D+00 - alfa ) * xl   &
                   + ( 1.0D+00 + alfa ) * xr ) &
                   /   2.0D+00
      alfa = +0.577350D+00
      xquad(2,i) = ( ( 1.0D+00 - alfa ) * xl   &
                   + ( 1.0D+00 + alfa ) * xr ) &
                   /   2.0D+00
    else if ( nquad == 3 ) then
      alfa = -0.774597D+00
      xquad(1,i) = ( ( 1.0D+00 - alfa ) * xl   &
                   + ( 1.0D+00 + alfa ) * xr ) &
                   /   2.0D+00
      xquad(2,i) = 0.5D+00 * ( xl + xr )
      alfa = +0.774597D+00
      xquad(3,i) = ( ( 1.0D+00 - alfa ) * xl   &
                   + ( 1.0D+00 + alfa ) * xr ) &
                   /   2.0D+00
    end if

  end do
!
!  Store the weights for the quadrature rule.
!
  if ( nquad == 1 ) then
    wquad(1) = 1.0D+00
  else if ( nquad == 2 ) then
    wquad(1) = 0.5D+00
    wquad(2) = 0.5D+00
  else if ( nquad == 3 ) then
    wquad(1) = 4.0D+00 / 9.0D+00
    wquad(2) = 5.0D+00 / 18.0D+00
    wquad(3) = 4.0D+00 / 9.0D+00
  end if

  return
end