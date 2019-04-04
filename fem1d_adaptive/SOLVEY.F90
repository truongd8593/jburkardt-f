subroutine solvey ( eta, f, h, n, nu, ul, ur, xn )

!*****************************************************************************80
!
!! SOLVEY computes error estimators for a finite element solution.
!
!  Discussion:
!
!    SOLVEY accepts information about the solution of a finite element
!    problem on a grid of nodes with coordinates XN.  It then starts
!    at node 0, and for each node, computes two "error estimators",
!    one for the left, and one for the right interval associated with the
!    node.  These estimators are found by solving a finite element problem
!    over the two intervals, using the known values of the original
!    solution as boundary data, and using a mesh that is "slightly"
!    refined over the original one.
!
!    Note that the computations at the 0-th and N-th nodes only involve
!    a single interval.
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
!    Output, real ( kind = 8 ) ETA(N).
!    ETA(I) is the error estimate for interval I.  It is computed
!    as the sum of two quantities, one associated with the left
!    and one with the right node of the interval.
!
!    Input, real ( kind = 8 ) F(NU).
!    ASSEMBLE stores into F the right hand side of the linear
!    equations.
!    SOLVE replaces those values of F by the solution of the
!    linear equations.
!
!    Input, real ( kind = 8 ) H(N)
!    H(I) is the length of subinterval I.  This code uses
!    equal spacing for all the subintervals.
!
!    Input, integer ( kind = 4 ) N
!    The number of subintervals into which the interval
!    [XL,XR] is broken.
!
!    Input, integer ( kind = 4 ) NU.
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
  implicit none

  integer ( kind = 4 ), parameter :: nl = 2
  integer ( kind = 4 ), parameter :: ny = 2
  integer ( kind = 4 ), parameter :: nquad = 2

  integer ( kind = 4 ), parameter :: nmay = 2 * ny

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nu

  real ( kind = 8 ) adiag(nmay)
  real ( kind = 8 ) aleft(nmay)
  real ( kind = 8 ) arite(nmay)
  real ( kind = 8 ) eta(n)
  real ( kind = 8 ) f(nu)
  real ( kind = 8 ) fy(nmay)
  real ( kind = 8 ) h(n)
  real ( kind = 8 ) hy(nmay)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ibcy
  integer ( kind = 4 ) indy(0:nmay)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  integer ( kind = 4 ) jmid
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) nodey(nl,nmay)
  integer ( kind = 4 ) nuy
  real ( kind = 8 ) pp
  real ( kind = 8 ) qq
  real ( kind = 8 ) total
  real ( kind = 8 ) ul
  real ( kind = 8 ) uleft
  real ( kind = 8 ) ulval
  real ( kind = 8 ) uly
  real ( kind = 8 ) uprime
  real ( kind = 8 ) ur
  real ( kind = 8 ) urite
  real ( kind = 8 ) urval
  real ( kind = 8 ) ury
  real ( kind = 8 ) uval
  real ( kind = 8 ) vlval
  real ( kind = 8 ) vprime
  real ( kind = 8 ) vrval
  real ( kind = 8 ) vval
  real ( kind = 8 ) wquad(nquad)
  real ( kind = 8 ) xn(0:n)
  real ( kind = 8 ) xquady(nquad,nmay)
  real ( kind = 8 ) y
  real ( kind = 8 ) yl
  real ( kind = 8 ) ym
  real ( kind = 8 ) yn(0:nmay)
  real ( kind = 8 ) yr
!
!  Initialize the error estimators to zero.
!
  eta(1:n) = 0.0D+00
!
!  Set the boundary conditions for each subproblem to be
!  known values of U at the left and right.
!
!
!  For each node, subdivide its left and right hand intervals
!  into NY subintervals.
!
!  Set up and solve the differential equation again on this
!  smaller region.
!
!  The 0-th and N-th nodes are special cases.
!
  ibcy = 3

  do j = 0, n

    if ( j == 0 ) then
      m = ny
      jlo = j
      jmid = j + 1
      jhi = j + 1
    else if ( j == n ) then
      m = ny
      jlo = j - 1
      jmid = j
      jhi = j
    else
      m = 2 * ny
      jlo = j - 1
      jmid = j
      jhi = j + 1
    end if
!
!  Set the location of the nodes in the subintervals.
!
    yl = xn(jlo)
    ym = xn(jmid)
    yr = xn(jhi)

    do i = 0, ny
      yn(i) = ( real ( ny - i, kind = 8 ) * yl   &
              + real (      i, kind = 8 ) * ym ) &
              / real ( ny,     kind = 8 )
    end do

    do i = ny + 1, m
      yn(i) = ( real ( m - i,      kind = 8 ) * ym   &
              + real (     i - ny, kind = 8 ) * yr ) &
              / real ( m -     ny, kind = 8 )
    end do
!
!  Set up the geometry of the sub-problem.
!
    call geometry ( hy, ibcy, indy, m, nl, nmay, nodey, nquad, nuy, &
      wquad, yn, xquady )
!
!  Set the boundary values for the sub-problem.
!
    if ( j <= 1 ) then
      uly = ul
    else
      uly = f(j-1)
    end if

    if ( n - 1 <= j ) then
      ury = ur
    else
      ury = f(j+1)
    end if
!
!  Assemble the matrix for the sub-problem.
!
    call assemble ( adiag, aleft, arite, fy, hy, m, indy, nodey, nuy, nl, &
      nquad, nmay, uly, ury, wquad, yn, xquady )
!
!  Solve the system.
!
    call solve ( adiag, aleft, arite, fy, nuy )
!
!  Compute the weighted sum of the squares of the differences
!  of the original computed slope and the refined computed slopes.
!
!  Calculation for left interval.
!
    if ( 1 <= j ) then

      if ( j <= 1 ) then
        uleft = ul
        urite = f(1)
      else if ( j == n ) then
        uleft = f(j-1)
        urite = ur
      else
        uleft = f(j-1)
        urite = f(j)
      end if

      uprime = ( urite - uleft ) / h(j)

      total = 0.0D+00
      do i = 1, ny

        yl = yn(i-1)
        yr = yn(i)

        if ( i == 1 ) then
          vlval = uly
          vrval = fy(i)
        else if ( i == m ) then
          vlval = fy(i-1)
          vrval = ury
        else
          vlval = fy(i-1)
          vrval = fy(i)
        end if

        vprime = ( vrval - vlval ) / hy(i)

        ulval = ( real ( ny - i + 1, kind = 8 ) * uleft   &
                + real (      i - 1, kind = 8 ) * urite ) &
                / real ( ny,         kind = 8 )

        urval = ( real ( ny - i, kind = 8 ) * uleft   &
                + real (      i, kind = 8 ) * urite ) &
                / real ( ny,     kind = 8 )
!
!  Compute the integral of
!
!    p(x)*(u'(x)-v'(x))**2 + q(x)*(u(x)-v(x))**2
!
        do k = 1, nquad

          y  =  xquady(k,i)

          uval = ( ( yl - y      ) * urval   &
                 + (      y - yr ) * ulval ) &
                 / ( yl     - yr )

          vval = ( ( yl - y      ) * vrval   &
                 + (      y - yr ) * vlval ) &
                 / ( yl     - yr )

          total = total + 0.5D+00 * wquad(k) * hy(i) * &
            ( pp ( y ) * ( uprime - vprime )**2 &
            + qq ( y ) * ( uval - vval )**2 )

        end do

      end do

      eta(j) = eta(j) + 0.5D+00 * sqrt ( total )

    end if
!
!  Calculation for right interval.
!
    if ( j <= n - 1 ) then

      if ( j == 0 ) then
        uleft = ul
        urite = f(j+1)
      else if ( n - 1 <= j ) then
        uleft = f(j)
        urite = ur
      else
        uleft = f(j)
        urite = f(j+1)
      end if

      uprime = ( urite - uleft ) / h(j+1)

      total = 0.0D+00
      do i = m+1-ny, m

        yl = yn(i-1)
        yr = yn(i)

        if ( i == 1 ) then
          vlval = uly
          vrval = fy(i)
        else if ( i == m ) then
          vlval = fy(i-1)
          vrval = ury
        else
          vlval = fy(i-1)
          vrval = fy(i)
        end if

        vprime = ( vrval - vlval ) / hy(i)

        ulval = ( real (      m - i + 1, kind = 8 ) * uleft   &
                + real ( ny - m + i - 1, kind = 8 ) * urite ) &
                / real ( ny,             kind = 8 )

        urval = ( real (      m - i, kind = 8 ) * uleft   &
                + real ( ny - m + i, kind = 8 ) * urite ) &
                / real ( ny,         kind = 8 )
!
!  Compute the integral of
!
!    p(x)*(u'(x)-v'(x))**2 + q(x)*(u(x)-v(x))**2
!
        do k = 1, nquad

          y  =  xquady(k,i)

          uval = ( ( yl - y      ) * urval   &
                 + (      y - yr ) * ulval ) &
                 / ( yl     - yr )

          vval = ( ( yl - y      ) * vrval   &
                 + (      y - yr ) * vlval ) &
                 / ( yl     - yr )

          total = total + 0.5D+00 * wquad(k) * hy(i) * &
            ( pp ( y ) * ( uprime - vprime )**2 &
            + qq ( y ) * ( uval - vval )**2 )

        end do

      end do

      eta(j+1) = eta(j+1) + 0.5D+00 * sqrt ( total )

    end if

  end do
!
!  Print out the error estimators.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ETA'
  write ( *, '(a)' ) ' '
  do j = 1, n
    write ( *, '(g14.6)' ) eta(j)
  end do

  return
end