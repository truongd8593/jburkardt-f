subroutine write_solution_tecplot ( f, ibc, indx, n, nu, ul, ur, xn )

!*****************************************************************************80
!
!! WRITE_SOLUTION_TECPLOT exports the computed solution in Tecplot .dat file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 March 2019
!
!  Author:
!
!    FORTRAN90 version by Truong Dang - modified version of the file OUTPUT.F90
!
!  Parameters:
!
!    Input, real ( kind = 8 ) F(NU).
!    ASSEMBLE stores into F the right hand side of the linear
!    equations.
!    SOLVE replaces those values of F by the solution of the
!    linear equations.
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
!    Input, integer ( kind = 4 ) INDX(0:N).
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

  integer ( kind = 4 ), intent(in) :: n
  integer ( kind = 4 ), intent(in) :: nu

  integer ( kind = 4 ), intent(in) :: ibc
  integer ( kind = 4 ), intent(in) :: indx(0:n)

  real ( kind = 8 ), intent(in) :: ul
  real ( kind = 8 ), intent(in) :: ur
  real ( kind = 8 ), intent(in) :: xn(0:n)
  real ( kind = 8 ), intent(in) :: f(nu)

  integer ( kind = 4 ) :: i
  
  real ( kind = 8 )  :: u
  real ( kind = 8 )  :: uex
  real ( kind = 8 )  :: u_exact
  real ( kind = 8 )  :: error

  character (4)      :: chr_nb_node
  

  write ( chr_nb_node, '(i4.4)' ) n

  open ( unit = 20, file = 'solutionTecplot_'//trim(chr_nb_node)//'.dat', status = 'new' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Node    X(I)        U(X(I))        U exact       Error'
  write ( *, '(a)' ) ' '

  write ( 20, '(a)' ) 'VARIABLES=X,U_APPROXIMATE,U_EXACT'

  do i = 0, n

    if ( i == 0 ) then

      if ( ibc == 1 .or. ibc == 3 ) then
        u = ul
      else
        u = f(indx(i))
      end if

    else if ( i == n ) then

      if ( ibc == 2 .or. ibc == 3 ) then
        u = ur
      else
        u = f(indx(i))
      end if

    else

      u = f(indx(i))

    end if

    uex = u_exact ( xn(i) )
    error = u - uex

    write(*,'(i4,4g14.6)') i, xn(i), u, uex, error
	write( 20, '(3f14.6)') xn(i), u, uex

  end do

  close ( unit = 20 )

  return
end