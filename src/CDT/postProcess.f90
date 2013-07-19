function areap ( x, y, nb, nodes )

!*****************************************************************************80
!
!! AREAP computes the signed area of a polygonal curve.
!
!  Discussion:
!
!    Given a sequence of NB points in the plane, this function
!    computes the signed area bounded by the closed polygonal
!    curve which passes through the points in the
!    specified order.  Each simple closed curve is positively
!    oriented (bounds positive area) if and only if the points
!    are specified in counterclockwise order.  The last point
!    of the curve is taken to be the first point specified, and
!    this point should therefore not be specified twice.
!
!    The area of a triangulation may be computed by calling
!    AREAP with values of NB and NODES determined by subroutine
!    BNODES.
!
!  Modified:
!
!    16 June 2007
!
!  Author:
!
!    Robert Renka,
!    Department of Computer Science,
!    University of North Texas,
!    renka@cs.unt.edu
!
!  Reference:
!
!    Robert Renka,
!    Algorithm 751: TRIPACK, 
!    A Constrained Two-Dimensional Delaunay Triangulation Package,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 1, 1996.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X(*), Y(*), the Cartesian coordinates of a set 
!    of points.
!
!    Input, integer ( kind = 4 ) NB, the number of points in the curve.
!
!    Input, integer ( kind = 4 ) NODES(NB), the indices of the points that
!    make up the closed curve.
!
!    Output, real ( kind = 8 ) AREAP, the signed area bounded by the curve.
!
  implicit none

  integer ( kind = 4 ) nb

  real ( kind = 8 ) a
  real ( kind = 8 ) areap
  integer ( kind = 4 ) i
  integer ( kind = 4 ) nd1
  integer ( kind = 4 ) nd2
  integer ( kind = 4 ) nodes(nb)
  real ( kind = 8 ) x(*)
  real ( kind = 8 ) y(*)

  a = 0.0D+00

  if ( nb < 3 ) then
    areap = 0.0D+00
    return
  end if

  nd2 = nodes(nb)
!
!  Loop on line segments NODES(I-1) -> NODES(I), where
!  NODES(0) = NODES(NB), adding twice the signed trapezoid
!  areas (integrals of the linear interpolants) to A.
!
  do i = 1, nb
    nd1 = nd2
    nd2 = nodes(i)
    a = a + ( x(nd2) - x(nd1) ) * ( y(nd1) + y(nd2) )
  end do
!
!  A contains twice the negative signed area of the region.
!
  areap = -a / 2.0D+00

  return

end function areap




subroutine bnodes ( n, list, lptr, lend, nodes, nb, na, nt )

!*****************************************************************************80
!
!! BNODES returns a list of the boundary nodes.
!
!  Discussion:
!
!    Given a triangulation of N points in the plane, this
!    subroutine returns an array containing the indexes, in
!    counterclockwise order, of the nodes on the boundary of
!    the convex hull of the set of points.
!
!  Modified:
!
!    16 June 2007
!
!  Author:
!
!    Robert Renka,
!    Department of Computer Science,
!    University of North Texas,
!    renka@cs.unt.edu
!
!  Reference:
!
!    Robert Renka,
!    Algorithm 751: TRIPACK, 
!    A Constrained Two-Dimensional Delaunay Triangulation Package,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 1, 1996.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of nodes in the triangulation.
!    3 <= N.
!
!    Input, integer ( kind = 4 ) LIST(*), LPTR(*), LEND(N), the data structure
!    defining the triangulation.  Refer to subroutine TRMESH.
!
!    Output, integer ( kind = 4 ) NODES(NB), ordered sequence of boundary node
!    indexes in the range 1 to N.
!
!    Output, integer ( kind = 4 ) NB, the number of boundary nodes.
!
!    Output, integer ( kind = 4 ) NA, NT, the number of arcs and triangles, 
!    respectively, in the triangulation.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) k
  integer ( kind = 4 ) lend(n)
  integer ( kind = 4 ) list(*)
  integer ( kind = 4 ) lp
  integer ( kind = 4 ) lptr(*)
  integer ( kind = 4 ) n0
  integer ( kind = 4 ) na
  integer ( kind = 4 ) nb
  integer ( kind = 4 ) nodes(*)
  integer ( kind = 4 ) nst
  integer ( kind = 4 ) nt
!
!  Set NST to the first boundary node encountered.
!
  nst = 1

  do

    lp = lend(nst)

    if ( list(lp) < 0 ) then
      exit
    end if

    nst = nst + 1

  end do
!
!  Initialization.
!
  nodes(1) = nst
  k = 1
  n0 = nst
!
!  Traverse the boundary in counterclockwise order.
!
  do

    lp = lend(n0)
    lp = lptr(lp)
    n0 = list(lp)

    if ( n0 == nst ) then
      exit
    end if

    k = k + 1
    nodes(k) = n0

  end do
!
!  Termination.
!
  nb = k
  nt = 2 * n - nb - 2
  na = nt + n - 1

  return

end subroutine bnodes




subroutine circum ( x1, y1, x2, y2, x3, y3, ratio, xc, yc, cr, sa, ar )

!*****************************************************************************80
!
!! CIRCUM determines the circumcenter (and more) of a triangle.
!
!  Discussion:
!
!    Given three vertices defining a triangle, this routine
!    returns the circumcenter, circumradius, signed
!    triangle area, and, optionally, the aspect ratio of the
!    triangle.
!
!  Modified:
!
!    16 June 2007
!
!  Author:
!
!    Robert Renka,
!    Department of Computer Science,
!    University of North Texas,
!    renka@cs.unt.edu
!
!  Reference:
!
!    Robert Renka,
!    Algorithm 751: TRIPACK, 
!    A Constrained Two-Dimensional Delaunay Triangulation Package,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 1, 1996.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X1, Y1, X2, Y2, X3, Y3, the coordinates of
!    the vertices.
!
!    Input, logical RATIO, is TRUE if and only if the aspect ratio is 
!    to be computed.
!
!    Output, real ( kind = 8 ) XC, YC, coordinates of the circumcenter (center
!    of the circle defined by the three points) unless SA = 0, in which XC 
!    and YC are not altered.
!
!    Output, real ( kind = 8 ) CR, the circumradius (radius of the circle
!    defined by the three points) unless SA = 0 (infinite radius), in which
!    case CR is not altered.
!
!    Output, real ( kind = 8 ) SA, the signed triangle area with positive value
!    if and only if the vertices are specified in counterclockwise order:  
!    (X3,Y3) is strictly to the left of the directed line from (X1,Y1)
!    toward (X2,Y2).
!
!    Output, real ( kind = 8 ) AR, the aspect ratio r/CR, where r is the 
!    radius of the inscribed circle, unless RATIO = FALSE, in which case AR
!    is not altered.  AR is in the range 0 to 0.5, with value 0 iff SA = 0 and
!    value 0.5 iff the vertices define an equilateral triangle.
!
  implicit none

  real ( kind = 8 ) ar
  real ( kind = 8 ) cr
  real ( kind = 8 ) ds(3)
  real ( kind = 8 ) fx
  real ( kind = 8 ) fy
  logical ratio
  real ( kind = 8 ) sa
  real ( kind = 8 ) u(3)
  real ( kind = 8 ) v(3)
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) x3
  real ( kind = 8 ) xc
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2
  real ( kind = 8 ) y3
  real ( kind = 8 ) yc
!
!  Set U(K) and V(K) to the x and y components, respectively,
!  of the directed edge opposite vertex K.
!
  u(1) = x3 - x2
  u(2) = x1 - x3
  u(3) = x2 - x1
  v(1) = y3 - y2
  v(2) = y1 - y3
  v(3) = y2 - y1
!
!  Set SA to the signed triangle area.
!
  sa = ( u(1) * v(2) - u(2) * v(1) ) / 2.0D+00

  if ( sa == 0.0D+00 ) then
    if ( ratio ) then
      ar = 0.0D+00
    end if
    return
  end if
!
!  Set DS(K) to the squared distance from the origin to vertex K.
!
  ds(1) = x1 * x1 + y1 * y1
  ds(2) = x2 * x2 + y2 * y2
  ds(3) = x3 * x3 + y3 * y3
!
!  Compute factors of XC and YC.
!
  fx = - dot_product ( ds(1:3), v(1:3) )
  fy =   dot_product ( ds(1:3), u(1:3) )

  xc = fx / ( 4.0D+00 * sa )
  yc = fy / ( 4.0D+00 * sa )
  cr = sqrt ( ( xc - x1 )**2 + ( yc - y1 )**2 )

  if ( .not. ratio ) then
    return
  end if
!
!  Compute the squared edge lengths and aspect ratio.
!
  ds(1:3) = u(1:3)**2 + v(1:3)**2

  ar = 2.0D+00 * abs ( sa ) / &
       ( ( sqrt ( ds(1) ) + sqrt ( ds(2) ) + sqrt ( ds(3) ) ) * cr )

  return

end subroutine circum
