subroutine subdiv ( eta, kount, n, nmax, tol, xn, status )

!*****************************************************************************80
!
!! SUBDIV decides which intervals should be subdivided.
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
!    Input, real ( kind = 8 ) ETA(N).
!    ETA(I) is the error estimate for interval I.  It is computed
!    as the sum of two quantities, one associated with the left
!    and one with the right node of the interval.
!
!    Input, integer ( kind = 4 ) KOUNT, the number of adaptive steps that 
!    have been taken.
!
!    Input/output, integer ( kind = 4 ) N
!    The number of subintervals into which the interval
!    [XL,XR] is broken.
!
!    Input, integer ( kind = 4 ) NMAX, the maximum number of unknowns that 
!    can be handled.
!
!    Input, real ( kind = 8 ) TOL.
!    A tolerance that is used to determine whether the estimated
!    error in an interval is so large that it should be subdivided
!    and the problem solved again.
!
!    Input/output, real ( kind = 8 ) XN(0:N).
!    XN(I) is the location of the I-th node.  XN(0) is XL,
!    and XN(N) is XR.
!
!    Output, integer ( kind = 4 ) STATUS, reports status of subdivision.
!    0, a new subdivision was carried out.
!    1, no more subdivisions are needed.
!    -1, no more subdivisions can be carried out.
!
!  Local Parameters:
!
!    Local, integer ( kind = 4 ) JADD(N).
!    JADD(I) is 1 if the error estimates show that interval I
!    should be subdivided.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nmax

  real ( kind = 8 ) ave
  real ( kind = 8 ) eta(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jadd(n)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kount
  integer ( kind = 4 ) status
  real ( kind = 8 ) temp
  real ( kind = 8 ) tol
  real ( kind = 8 ) xn(0:nmax)
  real ( kind = 8 ) xt(0:nmax)

  status = 0
!
!  Add up the ETA's, and get their average.
!
  ave = sum ( eta(1:n) ) / real ( n, kind = 8 )
!
!  Look for intervals whose ETA value is relatively large,
!  and note in JADD that these intervals should be subdivided.
!
  k = 0
  temp = max ( 1.2D+00 * ave + 0.00001D+00, tol**2 / real ( n, kind = 8 ) )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) 'Tolerance = ', temp
  write ( *, '(a)' ) ' '

  do j = 1, n

    if ( temp < eta(j) ) then
      k = k + 1
      jadd(j) = 1
      write ( *, '(a,i8)' ) 'Subdivide interval ', j
    else
      jadd(j) = 0
    end if

  end do
!
!  If no subdivisions needed, we're done.
!
  if ( k <= 0 ) then
    write ( *, '(a,i8)' ) 'Success on step ', kount
    status = 1
    return
  end if
!
!  See if we're about to go over our limit.
!
  if ( nmax < n + k ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'The iterations did not reach their goal.'
    write ( *, '(a,i8)' ) 'The next value of N is ', n + k
    write ( *, '(a,i8)' ) 'which exceeds NMAX = ', nmax
    status = -1
    return
  end if
!
!  Insert new nodes where needed.
!
  k = 0
  xt(0) = xn(0)
  do j = 1, n

    if ( 0 < jadd(j) ) then
      xt(j+k) = 0.5D+00 * ( xn(j) + xn(j-1) )
      k = k + 1
    end if

    xt(j+k) = xn(j)

  end do
!
!  Update the value of N, and copy the new nodes into XN.
!
  n = n + k

  xn(0:n) = xt(0:n)

  return
end