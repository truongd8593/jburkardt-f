!  FEM1D_ADAPTIVE.f90 
!
!  FUNCTIONS:
!	FEM1D_ADAPTIVE      - Entry point of console application.
!
!	Example of displaying 'Hello World' at execution time.
!

!****************************************************************************
!
!  PROGRAM: FEM1D_ADAPTIVE
!
!  PURPOSE:  Entry point for FEM1D_ADAPTIVE console application.
!
!****************************************************************************

	program FEM1D_ADAPTIVE

!*****************************************************************************80
!
!  main program for FEM1D_ADAPTIVE.
!
!  Discussion:
!
!    FEM1D_ADAPTIVE solves a 1D problem using an adaptive finite element method.
!
!    The equation to be treated is:
!
!      -d/dx ( P(x) * dU(x)/dx ) + Q(x) * U(x)  =  F(x)
!
!    by the finite-element method using piecewise linear basis
!    functions.
!
!    An adaptive method is used to try to reduce the maximum
!    error by refining the mesh in certain places.
!
!    Here U is an unknown scalar function of X defined on the
!    interval [XL,XR], and P, Q and F are given functions of X.
!
!    The values of U at XL and XR are also specified.
!
!    The interval [XL,XR] is "meshed" with N+1 points,
!
!      XN(0) = XL, XN(1) = XL+H, XN(2) = XL+2*H, ..., XN(N) = XR.
!
!    This creates N subintervals, with interval I having endpoints
!    XN(I-1) and XN(I).
!
!
!    The algorithm tries to guarantee a certain amount
!    of accuracy by examining the current solution, estimating the error
!    in each subinterval, and, if necessary, subdividing one or more
!    subintervals and repeating the calculation.
!
!    We can think of the adaptive part of the algorithm as a refined
!    problem.  The program re-solves the problem on the pair of
!    intervals J and J+1, which extend from node J-1 to node J+1.
!    The values of U that were just computed at nodes J-1 and J+1
!    will be used as the boundary values for this refined problem.
!    The intervals J and J+1 will each be evenly divided into NY
!    smaller subintervals.  This boundary value problem is solved,
!    and the derivatives of the original and refined solutions are
!    then compared to get an estimate of the error.
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
!    real ( kind = 8 ) ADIAG(NU).
!    ADIAG(I) is the "diagonal" coefficient of the I-th
!    equation in the linear system.  That is, ADIAG(I) is
!    the coefficient of the I-th unknown in the I-th equation.
!
!    real ( kind = 8 ) ALEFT(NU).
!    ALEFT(I) is the "left hand" coefficient of the I-th
!    equation in the linear system.  That is, ALEFT(I) is the
!    coefficient of the (I-1)-th unknown in the I-th equation.
!    There is no value in ALEFT(1), since the first equation
!    does not refer to a "0-th" unknown.
!
!    real ( kind = 8 ) ARITE(NU).
!    ARITE(I) is the "right hand" coefficient of the I-th
!    equation in the linear system.  ARITE(I) is the coefficient
!    of the (I+1)-th unknown in the I-th equation.  There is
!    no value in ARITE(NU) because the NU-th equation does not
!    refer to an "NU+1"-th unknown.
!
!    real ( kind = 8 ) ETA(N).
!    ETA(I) is the error estimate for interval I.  It is computed
!    as the sum of two quantities, one associated with the left
!    and one with the right node of the interval.
!
!    real ( kind = 8 ) F(NU).
!    ASSEMBLE stores into F the right hand side of the linear
!    equations.
!    SOLVE replaces those values of F by the solution of the
!    linear equations.
!
!    real ( kind = 8 ) FY(M).
!    FY is the right hand side of the linear system of the refined
!    problem.
!
!    real ( kind = 8 ) H(N)
!    H(I) is the length of subinterval I.  This code uses
!    equal spacing for all the subintervals.
!
!    real ( kind = 8 ) HY(M).
!    HY(I) is the length of subinterval I in the refined problem.
!
!    integer ( kind = 4 ) IBC.
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
!    integer ( kind = 4 ) IBCY.
!    IBCY declares the boundary conditions for the refined problem
!    which should always be that the value of U is specified at
!    both the left and right endpoints.  This corresponds to a
!    value of IBCY = 3.
!
!    integer ( kind = 4 ) INDX(0:N).
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
!    integer ( kind = 4 ) INDY(0:M).
!    INDY(I) records the index of the unknown associated with
!    node I for the refined problem.
!
!    integer ( kind = 4 ) JADD(N).
!    JADD(I) is 1 if the error estimates show that interval I
!    should be subdivided.
!
!    integer ( kind = 4 ) KOUNT, the number of adaptive steps that have 
!    been taken.
!
!    integer ( kind = 4 ) M.
!    M is the number of subintervals used in the refined problem.
!    M is equal to NY for computations centered at node 0 or node N,
!    and otherwise, M is equal to 2*NY.
!
!    integer ( kind = 4 ) N
!    The number of subintervals into which the interval
!    [XL,XR] is broken.
!
!    integer ( kind = 4 ) NL.
!    The number of basis functions used in a single
!    subinterval.  (NL-1) is the degree of the polynomials
!    used.  For this code, NL is fixed at 2, meaning that
!    piecewise linear functions are used as the basis.
!
!    integer ( kind = 4 ) NMAX, the maximum number of unknowns that can be
!    handled.
!
!    integer ( kind = 4 ) NODE(NL,N).
!    For each subinterval I:
!    NODE(1,I) is the number of the left node, and
!    NODE(2,I) is the number of the right node.
!
!    integer ( kind = 4 ) NODEY(NL,M).
!    NODEY performs the same function for the refined problem that
!    NODE performs for the full problem, recording the node numbers
!    associated with a particular subinterval.
!
!    integer ( kind = 4 ) NQUAD
!    The number of quadrature points used in a subinterval.
!
!    integer ( kind = 4 ) NU.
!    NU is the number of unknowns in the linear system.
!    Depending on the value of IBC, there will be N-1,
!    N, or N+1 unknown values, which are the coefficients
!    of basis functions.
!
!    integer ( kind = 4 ) NUY.
!    The number of unknowns in the refined problem.
!
!    integer ( kind = 4 ) NY.
!    NY is the number of subintervals into which a given interval
!    will be subdivided, before solving the refined probelm.
!
!    integer ( kind = 4 ) PROBLEM, chooses the problem to be solved.
!    The user must choose this value by setting it in routine GETPRB.
!    * 1, u = x, p = 1, q = 0, f = 0, ibc = 3, ul = 0, ur = 1.
!    The program should find the solution exactly, and the
!    adaptive code should find that there is no reason to
!    subdivide any interval.
!    * 2, u = x*x, p = 1, q = 0, f = -2, ibc = 3, ul = 0, ur = 1.
!    This problem should find the solution exactly, and
!    the adaptive code should again find there is nothing
!    to do.
!    *3, u = sin(pi*x/2), p = 1, q = 0, ibc = 3, f = 0.25*pi*pi*sin(pi*x/2),
!    ul = 0, ur = 1.
!    *4, u = cos(pi*x/2), p = 1, q = 0, ibc = 3, f = 0.25*pi*pi*cos(pi*x/2),
!    ul = 1, ur = 0.
!    *5: u = x**(beta+2)/((beta+2)*(beta+1)), p = 1, q = 1, ibc = 3,
!    f = -x**beta + (x**(beta+2))/((beta+2)*(beta+1)),
!    ul = 0, ur = 1/((beta+2)*(beta+1))
!    (beta must be greater than -2, and not equal to -1)
!    *6: u = atan((x-0.5)/alpha), p = 1, q = 0, ibc = 3,
!    f =  2*alpha*(x-0.5) / (alpha**2 + (x-0.5)**2) **2,
!    ul = u(0), ur = u(1)
!
!    integer ( kind = 4 ) STATUS, reports status of subdivision.
!    0, a new subdivision was carried out.
!    1, no more subdivisions are needed.
!    -1, no more subdivisions can be carried out.
!
!    real ( kind = 8 ) TOL.
!    A tolerance that is used to determine whether the estimated
!    error in an interval is so large that it should be subdivided
!    and the problem solved again.
!
!    real ( kind = 8 ) UL.
!    If IBC is 1 or 3, UL is the value that U is required
!    to have at X = XL.
!    If IBC is 2 or 4, UL is the value that U' is required
!    to have at X = XL.
!
!    real ( kind = 8 ) UR.
!    If IBC is 2 or 3, UR is the value that U is required
!    to have at X = XR.
!    If IBC is 1 or 4, UR is the value that U' is required
!    to have at X = XR.
!
!    real ( kind = 8 ) WQUAD(NQUAD).
!    WQUAD(I) is the weight associated with the I-th point
!    of an NQUAD point Gaussian quadrature rule.
!
!    real ( kind = 8 ) XL.
!    XL is the left endpoint of the interval over which the
!    differential equation is being solved.
!
!    real ( kind = 8 ) XN(0:N).
!    XN(I) is the location of the I-th node.  XN(0) is XL,
!    and XN(N) is XR.
!
!    real ( kind = 8 ) XQUAD(NQUAD,NMAX), the I-th quadrature point
!    in interval J.
!
!    real ( kind = 8 ) XQUADY(NQUAD,NMAY ), the I-th quadrature point
!    in subinterval J of the refined problem.
!
!    real ( kind = 8 ) XR.
!    XR is the right endpoint of the interval over which the
!    differential equation is being solved.
!
!    Workspace, double precision XT(0:NMAX), used to compute a new
!    set of nodes.
!
!    real ( kind = 8 ) YN(0:M).
!    YN(I) is the location of the I-th node in the refined
!    problem.
!
  implicit none

  integer ( kind = 4 ), parameter :: nl = 2
  integer ( kind = 4 ), parameter :: nmax = 30
  integer ( kind = 4 ), parameter :: nquad = 2

  real ( kind = 8 ) adiag(nmax)
  real ( kind = 8 ) aleft(nmax)
  real ( kind = 8 ) alpha
  real ( kind = 8 ) arite(nmax)
  real ( kind = 8 ) beta
  real ( kind = 8 ) eta(nmax)
  real ( kind = 8 ) f(nmax)
  real ( kind = 8 ) h(nmax)
  integer ( kind = 4 ) ibc
  integer ( kind = 4 ) indx(0:nmax)
  integer ( kind = 4 ) jadd(nmax)
  integer ( kind = 4 ) kount
  integer ( kind = 4 ) n
  integer ( kind = 4 ) node(nl,nmax)
  integer ( kind = 4 ) nu
  integer ( kind = 4 ) problem
  integer ( kind = 4 ) status
  real ( kind = 8 ) tol
  real ( kind = 8 ) ul
  real ( kind = 8 ) ur
  real ( kind = 8 ) wquad(nquad)
  real ( kind = 8 ) xl
  real ( kind = 8 ) xn(0:nmax)
  real ( kind = 8 ) xquad(nquad,nmax)
  real ( kind = 8 ) xr
  real ( kind = 8 ) xt(0:nmax)

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FEM1D_ADAPTIVE'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Solve the two-point boundary value problem:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  -d/dx ( P(x) * dU(x)/dx ) + Q(x) * U(x)  =  F(x)'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'on the interval [0,1], specifying the value'
  write ( *, '(a)' ) 'of U at each endpoint.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) 'The number of basis functions per element is ', nl
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) 'The number of quadrature points per element is ', nquad

  call get_problem ( problem )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Problem index = ', problem
  write ( *, '(a)' ) ' '

  if ( problem == 1 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  "Linear" problem:'
    write ( *, '(a)' ) '  (No refinement needed)'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  U(X) =  X'
    write ( *, '(a)' ) '  P(X) =  1.0'
    write ( *, '(a)' ) '  Q(X) =  0.0'
    write ( *, '(a)' ) '  F(X) =  0.0'
    write ( *, '(a)' ) '  IBC  =  3'
    write ( *, '(a)' ) '  UL   =  0.0'
    write ( *, '(a)' ) '  UR   =  1.0'

  else if ( problem == 2 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  "Quadratic" problem:'
    write ( *, '(a)' ) '  (No refinement needed)'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  U(X) =  X*X'
    write ( *, '(a)' ) '  P(X) =  1.0'
    write ( *, '(a)' ) '  Q(X) =  0.0'
    write ( *, '(a)' ) '  F(X) = -2.0'
    write ( *, '(a)' ) '  IBC  =  3'
    write ( *, '(a)' ) '  UL   =  0.0'
    write ( *, '(a)' ) '  UR   =  1.0'

  else if ( problem == 3 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  "SINE" problem:'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  U(X) =  SIN(PI*X/2)'
    write ( *, '(a)' ) '  P(X) =  1.0'
    write ( *, '(a)' ) '  Q(X) =  0.0'
    write ( *, '(a)' ) '  F(X) =  PI*PI*SIN(PI*X/2)/4'
    write ( *, '(a)' ) '  IBC  =  3'
    write ( *, '(a)' ) '  UL   =  0.0'
    write ( *, '(a)' ) '  UR   =  1.0'

  else if ( problem == 4 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  "COSINE" problem:'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  U(X) =  COS(PI*X/2)'
    write ( *, '(a)' ) '  P(X) =  1.0'
    write ( *, '(a)' ) '  Q(X) =  0.0'
    write ( *, '(a)' ) '  F(X) =  PI*PI*COS(PI*X/2)/4'
    write ( *, '(a)' ) '  IBC  =  3'
    write ( *, '(a)' ) '  UL   =  0.0'
    write ( *, '(a)' ) '  UR   =  1.0'

  else if ( problem == 5 ) then

    call get_beta ( beta )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  "RHEINBOLDT" problem:'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  U(X) =  X**(B+2)/((B+2)*(B+1))'
    write ( *, '(a)' ) '  P(X) =  1.0'
    write ( *, '(a)' ) '  Q(X) =  1.0'
    write ( *, '(a)' ) '  F(X) =  -X**B+(X**B+2))/((B+2)*(B+1))'
    write ( *, '(a)' ) '  IBC  =  3'
    write ( *, '(a)' ) '  UL   =  0.0'
    write ( *, '(a)' ) '  UR   =  1/((B+2)*(B+1))'
    write ( *, '(a,g14.6)' ) '  B    = ', beta

  else if ( problem == 6 ) then

    call get_alpha ( alpha )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  "ARCTAN" problem:'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  U(X) =  ATAN((X-0.5)/A)'
    write ( *, '(a)' ) '  P(X) =  1.0'
    write ( *, '(a)' ) '  Q(X) =  0.0'
    write ( *, '(a)' ) '  F(X) =  2*A*(X-0.5)/(A**2+(X-0.5)**2)**2'
    write ( *, '(a)' ) '  IBC  =  3'
    write ( *, '(a)' ) '  UL   =  ATAN(-0.5/A)'
    write ( *, '(a)' ) '  UR   =  ATAN( 0.5/A)'
    write ( *, '(a,g14.6)' ) '  A    = ', alpha

  end if
!
!  Start out with just 4 subintervals.
!
  n = 4
!
!  Initialize values that define the problem.
!
  call init ( ibc, n, tol, ul, ur, xl, xn, xr )
!
!  Start the iteration counter off at 0.
!
  kount = 0
!
!  Begin the next iteration.
!
  do

    kount = kount + 1

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8,a)' ) 'Begin new iteration with ', n, ' nodes.'
    write ( *, '(a)' ) ' '
!
!  Solve the regular problem.
!
    call solvex ( f, h, ibc, kount, n, nl, nmax, nu, ul, ur, xn )
!
!  Solve N subproblems to get the error estimators.
!
    call solvey ( eta, f, h, n, nu, ul, ur, xn )
!
!  Examine the error estimators, and see how many intervals should
!  be subdivided.
!
    call subdiv ( eta, kount, n, nmax, tol, xn, status )

    if ( status /= 0 ) then
      exit
    end if
!
!  Solve the problem again, with the new nodes.
!
  end do
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FEM1D_ADAPTIVE:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end program FEM1D_ADAPTIVE

