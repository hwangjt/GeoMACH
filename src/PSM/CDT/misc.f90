subroutine insert ( k, lp, list, lptr, lnew )

!*****************************************************************************80
!
!! INSERT inserts K as a neighbor of N1.
!
!  Discussion:
!
!    This subroutine inserts K as a neighbor of N1 following
!    N2, where LP is the LIST pointer of N2 as a neighbor of
!    N1.  Note that, if N2 is the last neighbor of N1, K will
!    become the first neighbor (even if N1 is a boundary node).
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
!    Input, integer ( kind = 4 ) K, the index of the node to be inserted.
!
!    Input, integer ( kind = 4 ) LP, the LIST pointer of N2 as a neighbor of N1.
!
!    Input/output, integer ( kind = 4 ) LIST(*), LPTR(*), LNEW, the data 
!    structure defining the triangulation.  Refer to subroutine TRMESH.  On
!    output, the data structure has been updated to include node K.
!
  implicit none

  integer ( kind = 4 ) k
  integer ( kind = 4 ) list(*)
  integer ( kind = 4 ) lnew
  integer ( kind = 4 ) lp
  integer ( kind = 4 ) lptr(*)
  integer ( kind = 4 ) lsav

  lsav = lptr(lp)
  lptr(lp) = lnew
  list(lnew) = k
  lptr(lnew) = lsav
  lnew = lnew + 1

  return

end subroutine insert




function jrand ( n, ix, iy, iz )

!*****************************************************************************80
!
!! JRAND returns a uniformly distributed random integer between 1 and N.
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
!    Brian Wichmann, David Hill, 
!    An Efficient and Portable Pseudo-random Number Generator,
!    Applied Statistics, 
!    Volume 31, Number 2, 1982, pages 188-190.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the maximum value to be returned.
!
!    Input/output, integer ( kind = 4 ) IX, IY, IZ, seeds initialized to values
!    in the range 1 to 30,000 before the first call to JRAND, and not altered 
!    by the user between subsequent calls (unless a sequence of random 
!    numbers is to be repeated by reinitializing the seeds).
!
!    Output, integer ( kind = 4 ) JRAND, random integer in the range 1 to N.
!
!  Local parameters:
!
!    U = Pseudo-random number uniformly distributed in the interval (0,1).
!    X = Pseudo-random number in the range 0 to 3 whose fractional part is U.
!
  implicit none

  integer ( kind = 4 ) ix
  integer ( kind = 4 ) iy
  integer ( kind = 4 ) iz
  integer ( kind = 4 ) jrand
  integer ( kind = 4 ) n
  real ( kind = 8 ) u
  real ( kind = 8 ) x

  ix = mod ( 171 * ix, 30269 )
  iy = mod ( 172 * iy, 30307 )
  iz = mod ( 170 * iz, 30323 )

  x = ( real ( ix, kind = 8 ) / 30269.0D+00 ) &
    + ( real ( iy, kind = 8 ) / 30307.0D+00 ) &
    + ( real ( iz, kind = 8 ) / 30323.0D+00 )
 
  u = x - int ( x )
  jrand = real ( n, kind = 8 ) * u + 1.0D+00

  return

end function jrand




function left ( x1, y1, x2, y2, x0, y0 )

!*****************************************************************************80
!
!! LEFT determines whether a node is to the left of a line.
!
!  Discussion:
!
!    This function determines whether node N0 is to the left
!    or to the right of the line through N1-N2 as viewed by an
!    observer at N1 facing N2.
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
!    Input, real ( kind = 8 ) X1, Y1, coordinates of N1.
!
!    Input, real ( kind = 8 ) X2, Y2, coordinates of N2.
!
!    Input, real ( kind = 8 ) X0, Y0, coordinates of N0.
!
!    Output, logical LEFT, is .TRUE. if and only if (X0,Y0) is on or 
!    to the left of the directed line N1->N2.
!
!  Local parameters:
!
!    DX1,DY1 = X,Y components of the vector N1->N2
!    DX2,DY2 = X,Y components of the vector N1->N0
!
  implicit none

  real ( kind = 8 ) dx1
  real ( kind = 8 ) dx2
  real ( kind = 8 ) dy1
  real ( kind = 8 ) dy2
  logical left
  real ( kind = 8 ) x0
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) y0
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2

  dx1 = x2 - x1
  dy1 = y2 - y1
  dx2 = x0 - x1
  dy2 = y0 - y1
!
!  If the sign of the vector cross product of N1->N2 and
!  N1->N0 is positive, then sin(A) > 0, where A is the
!  angle between the vectors, and thus A is in the range
!  (0,180) degrees.
!
  left = ( dx1 * dy2 >= dx2 * dy1 )

  return

end function left




function lstptr ( lpl, nb, list, lptr )

!*****************************************************************************80
!
!! LSTPTR returns the index of NB in the adjacency list for N0.
!
!  Discussion:
!
!    This function returns the index (LIST pointer) of NB in
!    the adjacency list for N0, where LPL = LEND(N0).
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
!    Input, integer ( kind = 4 ) LPL = LEND(N0).
!
!    Input, integer ( kind = 4 ) NB, the index of the node whose pointer is to 
!    be returned.  NB must be connected to N0.
!
!    Input, integer ( kind = 4 ) LIST(*), LPTR(*), the data structure defining 
!    the triangulation.  Refer to subroutine TRMESH.
!
!    Output, integer ( kind = 4 ) LSTPTR, pointer such that LIST(LSTPTR) = NB or
!    LIST(LSTPTR) = -NB, unless NB is not a neighbor of N0, in which 
!    case LSTPTR = LPL.
!
  implicit none

  integer ( kind = 4 ) list(*)
  integer ( kind = 4 ) lp
  integer ( kind = 4 ) lpl
  integer ( kind = 4 ) lptr(*)
  integer ( kind = 4 ) lstptr
  integer ( kind = 4 ) nb
  integer ( kind = 4 ) nd

  lp = lptr(lpl)

  do

    nd = list(lp)

    if ( nd == nb ) then
      exit
    end if

    lp = lptr(lp)

    if ( lp == lpl ) then
      exit
    end if

  end do

  lstptr = lp

  return

end function lstptr




function nbcnt ( lpl, lptr )

!*****************************************************************************80
!
!! NBCNT returns the number of neighbors of a node.
!
!  Discussion:
!
!    This function returns the number of neighbors of a node
!    in a triangulation created by TRMESH or TRMSHR.
!
!  Modified:
!
!    25 November 2002
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
!    Input, integer ( kind = 4 ) LPL, the LIST pointer to the last neighbor 
!    of N0.  LPL = LEND(N0).
!
!    Input, integer ( kind = 4 ) LPTR(*), pointers associated with LIST.
!
!    Output, integer ( kind = 4 ) NBCNT, the  number of neighbors of N0.
!
  implicit none

  integer ( kind = 4 ) k
  integer ( kind = 4 ) lp
  integer ( kind = 4 ) lpl
  integer ( kind = 4 ) lptr(*)
  integer ( kind = 4 ) nbcnt

  lp = lpl
  k = 1

  do

    lp = lptr(lp)
    if ( lp == lpl ) then
      exit
    end if

    k = k + 1

  end do

  nbcnt = k

  return

end function nbcnt




subroutine optim ( x, y, na, list, lptr, lend, nit, iwk, ier )

!*****************************************************************************80
!
!! OPTIM optimizes the quadrilateral portion of a triangulation.
!
!  Discussion:
!
!    Given a set of NA triangulation arcs, this subroutine
!    optimizes the portion of the triangulation consisting of
!    the quadrilaterals (pairs of adjacent triangles) which
!    have the arcs as diagonals by applying the circumcircle
!    test and appropriate swaps to the arcs.
!
!    An iteration consists of applying the swap test and
!    swaps to all NA arcs in the order in which they are
!    stored.  The iteration is repeated until no swap occurs
!    or NIT iterations have been performed.  The bound on the
!    number of iterations may be necessary to prevent an
!    infinite loop caused by cycling (reversing the effect of a
!    previous swap) due to floating point inaccuracy when four
!    or more nodes are nearly cocircular.
!
!  Modified:
!
!    20 October 2005
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
!    Input, real ( kind = 8 ) X(*), Y(*), the nodal coordinates.
!
!    Input, integer ( kind = 4 ) NA, the number of arcs in the set.  0 <= NA.
!
!    Input/output, integer ( kind = 4 ) LIST(*), LPTR(*), LEND(N), data 
!    structure defining the triangulation.  Refer to subroutine TRMESH.
!    On output, updated to reflect the swaps.
!
!    Input/output, integer ( kind = 4 ) NIT.  On input, the maximum number of 
!    iterations to be performed.  A reasonable value is 3*NA.  1 <= NIT.  On 
!    output, the number of iterations performed.
!
!    Input/output, integer ( kind = 4 ) IWK(2,NA), containing the nodal indexes 
!    of the arc endpoints (pairs of endpoints are stored in columns).  On 
!    output, the information has been updated to reflect the swaps.
!
!    Output, integer ( kind = 4 ) IER = Error indicator:
!    0, if no errors were encountered.
!    1, if a swap occurred on the last of MAXIT iterations, where MAXIT is 
!      the value of NIT on input.  The new set of arcs in not necessarily
!      optimal in this case.
!    2, if NA < 0 or NIT < 1 on input.
!    3, if IWK(2,I) is not a neighbor of IWK(1,I) for some I in the range 1
!      to NA.  A swap may have occurred in this case.
!    4, if a zero pointer was returned by subroutine SWAP.
!
!  Local parameters:
!
!    I =       Column index for IWK
!    IO1,IO2 = Nodal indexes of the endpoints of an arc in IWK
!    ITER =    Iteration count
!    LP =      LIST pointer
!    LP21 =    Parameter returned by SWAP (not used)
!    LPL =     Pointer to the last neighbor of IO1
!    LPP =     Pointer to the node preceding IO2 as a neighbor of IO1
!    MAXIT =   Input value of NIT
!    N1,N2 =   Nodes opposite IO1->IO2 and IO2->IO1, respectively
!    NNA =     Local copy of NA
!    SWP =     Flag set to TRUE iff a swap occurs in the optimization loop
!
  implicit none

  integer ( kind = 4 ) na

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ier
  integer ( kind = 4 ) io1
  integer ( kind = 4 ) io2
  integer ( kind = 4 ) iter
  integer ( kind = 4 ) iwk(2,na)
  integer ( kind = 4 ) lend(*)
  integer ( kind = 4 ) list(*)
  integer ( kind = 4 ) lp
  integer ( kind = 4 ) lp21
  integer ( kind = 4 ) lpl
  integer ( kind = 4 ) lpp
  integer ( kind = 4 ) lptr(*)
  integer ( kind = 4 ) maxit
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) nit
  integer ( kind = 4 ) nna
  logical swp
  logical swptst
  real ( kind = 8 ) x(*)
  real ( kind = 8 ) y(*)

  nna = na
  maxit = nit

  if ( nna < 0 .or. maxit < 1 ) then
    nit = 0
    ier = 2
    return
  end if
!
!  Initialize iteration count ITER and test for NA = 0.
!
  iter = 0

  if ( nna == 0 ) then
    nit = 0
    ier = 0
    return
  end if
!
!  Top of loop.
!  SWP = TRUE iff a swap occurred in the current iteration.
!
  do

    if ( iter == maxit ) then
      nit = maxit
      ier = 1
      return
    end if

    iter = iter + 1
    swp = .false.
!
!  Inner loop on arcs IO1-IO2.
!
    do i = 1, nna

      io1 = iwk(1,i)
      io2 = iwk(2,i)
!
!  Set N1 and N2 to the nodes opposite IO1->IO2 and
!  IO2->IO1, respectively.  Determine the following:
!
!  LPL = pointer to the last neighbor of IO1,
!  LP = pointer to IO2 as a neighbor of IO1, and
!  LPP = pointer to the node N2 preceding IO2.
!
      lpl = lend(io1)
      lpp = lpl
      lp = lptr(lpp)

      do

        if ( list(lp) == io2 ) then
          go to 3
        end if

        lpp = lp
        lp = lptr(lpp)
        if ( lp == lpl ) then
          exit
        end if

      end do
!
!  IO2 should be the last neighbor of IO1.  Test for no
!  arc and bypass the swap test if IO1 is a boundary node.
!
      if ( abs ( list(lp) ) /= io2 ) then
        nit = iter
        ier = 3
        return
      end if

      if ( list(lp) < 0 ) then
        go to 4
      end if
!
!  Store N1 and N2, or bypass the swap test if IO1 is a
!  boundary node and IO2 is its first neighbor.
!
3     continue

      n2 = list(lpp)

      if ( n2 < 0 ) then
        go to 4
      end if

      lp = lptr(lp)
      n1 = abs ( list(lp) )
!
!  Test IO1-IO2 for a swap, and update IWK if necessary.
!
      if ( .not. swptst ( n1, n2, io1, io2, x, y ) ) then
        go to 4
      end if

      call swap ( n1, n2, io1, io2, list, lptr, lend, lp21 )

      if ( lp21 == 0 ) then
        nit = iter
        ier = 4
        return
      end if

      swp = .true.
      iwk(1,i) = n1
      iwk(2,i) = n2

4     continue

    end do

    if ( .not. swp ) then
      exit
    end if

  end do
!
!  Successful termination.
!
5 continue

  nit = iter
  ier = 0

  return

end subroutine optim




subroutine swap ( in1, in2, io1, io2, list, lptr, lend, lp21 )

!*****************************************************************************80
!
!! SWAP adjusts a triangulation by swapping a diagonal arc.
!
!  Discussion:
!
!    Given a triangulation of a set of points on the unit
!    sphere, this subroutine replaces a diagonal arc in a
!    strictly convex quadrilateral (defined by a pair of adja-
!    cent triangles) with the other diagonal.  Equivalently, a
!    pair of adjacent triangles is replaced by another pair
!    having the same union.
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
!    Input, integer ( kind = 4 ) IN1, IN2, IO1, IO2, the nodal indexes of the 
!    vertices of the quadrilateral.  IO1-IO2 is replaced by IN1-IN2.  
!    (IO1,IO2,IN1) and (IO2,IO1,IN2) must be triangles on input.
!
!    Input/output, integer ( kind = 4 ) LIST(*), LPTR(*), LEND(N), the data 
!    structure defining the triangulation.  Refer to subroutine TRMESH.  On 
!    output, updated with the swap; triangles (IO1,IO2,IN1) and (IO2,IO1,IN2) 
!    are replaced by (IN1,IN2,IO2) and (IN2,IN1,IO1) unless LP21 = 0.
!
!    Output, integer ( kind = 4 ) LP21, the index of IN1 as a neighbor of IN2 
!    after the swap is performed unless IN1 and IN2 are adjacent on input, in 
!    which case LP21 = 0.
!
!  Local parameters:
!
!    LP, LPH, LPSAV = LIST pointers
!
  implicit none

  integer ( kind = 4 ) in1
  integer ( kind = 4 ) in2
  integer ( kind = 4 ) io1
  integer ( kind = 4 ) io2
  integer ( kind = 4 ) lend(*)
  integer ( kind = 4 ) list(*)
  integer ( kind = 4 ) lp
  integer ( kind = 4 ) lp21
  integer ( kind = 4 ) lph
  integer ( kind = 4 ) lpsav
  integer ( kind = 4 ) lptr(*)
  integer ( kind = 4 ) lstptr
!
!  Test for IN1 and IN2 adjacent.
!
  lp = lstptr(lend(in1),in2,list,lptr)

  if ( abs ( list(lp) ) == in2 ) then
    lp21 = 0
    return
  end if
!
!  Delete IO2 as a neighbor of IO1.
!
  lp = lstptr(lend(io1),in2,list,lptr)
  lph = lptr(lp)
  lptr(lp) = lptr(lph)
!
!  If IO2 is the last neighbor of IO1, make IN2 the last neighbor.
!
  if ( lend(io1) == lph ) then
    lend(io1) = lp
  end if
!
!  Insert IN2 as a neighbor of IN1 following IO1
!  using the hole created above.
!
  lp = lstptr(lend(in1),io1,list,lptr)
  lpsav = lptr(lp)
  lptr(lp) = lph
  list(lph) = in2
  lptr(lph) = lpsav
!
!  Delete IO1 as a neighbor of IO2.
!
  lp = lstptr(lend(io2),in1,list,lptr)
  lph = lptr(lp)
  lptr(lp) = lptr(lph)
!
!  If IO1 is the last neighbor of IO2, make IN1 the last neighbor.
!
  if ( lend(io2) == lph ) then
    lend(io2) = lp
  end if
!
!  Insert IN1 as a neighbor of IN2 following IO2.
!
  lp = lstptr(lend(in2),io2,list,lptr)
  lpsav = lptr(lp)
  lptr(lp) = lph
  list(lph) = in1
  lptr(lph) = lpsav
  lp21 = lph

  return

end subroutine swap




function swptst ( in1, in2, io1, io2, x, y )

!*****************************************************************************80
!
!! SWPTST applies the circumcircle test to a quadrilateral.
!
!  Discussion:
!
!    This function applies the circumcircle test to a quadri-
!    lateral defined by a pair of adjacent triangles.  The
!    diagonal arc (shared triangle side) should be swapped for
!    the other diagonl if and only if the fourth vertex is
!    strictly interior to the circumcircle of one of the
!    triangles (the decision is independent of the choice of
!    triangle).  Equivalently, the diagonal is chosen to maxi-
!    mize the smallest of the six interior angles over the two
!    pairs of possible triangles (the decision is for no swap
!    if the quadrilateral is not strictly convex).
!
!    When the four vertices are nearly cocircular (the
!    neutral case), the preferred decision is no swap -- in
!    order to avoid unnecessary swaps and, more important, to
!    avoid cycling in subroutine OPTIM which is called by
!    DELNOD and EDGE.  Thus, a tolerance SWTOL (stored in
!    SWPCOM by TRMESH or TRMSHR) is used to define 'nearness'
!    to the neutral case.
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
!    Input, integer ( kind = 4 ) IN1, IN2, IO1, IO2, the nodal indexes of the 
!    vertices of the quadrilateral.  IO1-IO2 is the triangulation arc (shared 
!    triangle side) to be replaced by IN1-IN2 if the decision is to swap.  The
!    triples (IO1,IO2,IN1) and (IO2,IO1,IN2) must define triangles (be
!    in counterclockwise order) on input.
!
!    Input, real ( kind = 8 ) X(*), Y(*), the nodal coordinates.
!
!    Output, logical SWPTST, .TRUE. if and only if the arc connecting
!    IO1 and IO2 is to be replaced.
!
!  Local parameters:
!
!    DX11,DY11 = X,Y components of the vector IN1->IO1
!    DX12,DY12 = X,Y components of the vector IN1->IO2
!    DX22,DY22 = X,Y components of the vector IN2->IO2
!    DX21,DY21 = X,Y components of the vector IN2->IO1
!    SIN1 =      Cross product of the vectors IN1->IO1 and
!                IN1->IO2 -- proportional to sin(T1), where
!                T1 is the angle at IN1 formed by the vectors
!    COS1 =      Inner product of the vectors IN1->IO1 and
!                IN1->IO2 -- proportional to cos(T1)
!    SIN2 =      Cross product of the vectors IN2->IO2 and
!                IN2->IO1 -- proportional to sin(T2), where
!                T2 is the angle at IN2 formed by the vectors
!    COS2 =      Inner product of the vectors IN2->IO2 and
!                IN2->IO1 -- proportional to cos(T2)
!    SIN12 =     SIN1*COS2 + COS1*SIN2 -- proportional to sin(T1+T2)
!
  implicit none

  real ( kind = 8 ) cos1
  real ( kind = 8 ) cos2
  real ( kind = 8 ) dx11
  real ( kind = 8 ) dx12
  real ( kind = 8 ) dx21
  real ( kind = 8 ) dx22
  real ( kind = 8 ) dy11
  real ( kind = 8 ) dy12
  real ( kind = 8 ) dy21
  real ( kind = 8 ) dy22
  integer ( kind = 4 ) in1
  integer ( kind = 4 ) in2
  integer ( kind = 4 ) io1
  integer ( kind = 4 ) io2
  real ( kind = 8 ) sin1
  real ( kind = 8 ) sin12
  real ( kind = 8 ) sin2
  logical swptst
  real ( kind = 8 ) swtol
  real ( kind = 8 ) x(*)
  real ( kind = 8 ) y(*)
!
!  Tolerance stored by TRMESH or TRMSHR.
!
  common /swpcom/ swtol
!
!  Compute the vectors containing the angles T1 and T2.
!
  dx11 = x(io1) - x(in1)
  dx12 = x(io2) - x(in1)
  dx22 = x(io2) - x(in2)
  dx21 = x(io1) - x(in2)

  dy11 = y(io1) - y(in1)
  dy12 = y(io2) - y(in1)
  dy22 = y(io2) - y(in2)
  dy21 = y(io1) - y(in2)
!
!  Compute inner products.
!
  cos1 = dx11 * dx12 + dy11 * dy12
  cos2 = dx22 * dx21 + dy22 * dy21
!
!  The diagonals should be swapped iff 180 < (T1+T2)
!  degrees.  The following two tests ensure numerical
!  stability:  the decision must be FALSE when both
!  angles are close to 0, and TRUE when both angles
!  are close to 180 degrees.
!
  if ( 0.0D+00 <= cos1 .and. 0.0D+00 <= cos2 ) then
    swptst = .false.
    return
  end if

  if ( cos1 < 0.0D+00 .and. cos2 < 0.0D+00 ) then
    swptst = .true.
    return
  end if
!
!  Compute vector cross products (Z-components).
!
  sin1 = dx11 * dy12 - dx12 * dy11
  sin2 = dx22 * dy21 - dx21 * dy22
  sin12 = sin1 * cos2 + cos1 * sin2

  if ( -swtol <= sin12 ) then
    swptst = .false.
  else
    swptst = .true.
  end if

  return
end function swptst




subroutine timestamp ( )

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    31 May 2001   9:45:54.872 AM
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 May 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none

  character ( len = 8 ) ampm
  integer ( kind = 4 ) d
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  integer ( kind = 4 ) values(8)
  integer ( kind = 4 ) y

  call date_and_time ( values = values )

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  if ( h < 12 ) then
    ampm = 'AM'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
      ampm = 'Noon'
    else
      ampm = 'PM'
    end if
  else
    h = h - 12
    if ( h < 12 ) then
      ampm = 'PM'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'Midnight'
      else
        ampm = 'AM'
      end if
    end if
  end if

  write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return

end subroutine timestamp




subroutine trfind ( nst, px, py, n, x, y, list, lptr, lend, i1, i2, i3 )

!*****************************************************************************80
!
!! TRFIND locates a point relative to a triangulation.
!
!  Discussion:
!
!    This subroutine locates a point P relative to a triangu-
!    lation created by subroutine TRMESH or TRMSHR.  If P is
!    contained in a triangle, the three vertex indexes are
!    returned.  Otherwise, the indexes of the rightmost and
!    leftmost visible boundary nodes are returned.
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
!    Input, integer ( kind = 4 ) NST, the index of a node at which TRFIND begins
!    the search.  Search time depends on the proximity of this node to P.
!
!    Input, real ( kind = 8 ) PX, PY, the coordinates of the point P to be
!    located.
!
!    Input, integer ( kind = 4 ) N, the number of nodes in the triangulation.  
!    3 <= N.
!
!    Input, real ( kind = 8 ) X(N), Y(N), the coordinates of the nodes in
!    the triangulation.
!
!    Input, integer ( kind = 4 ) LIST(*), LPTR(*), LEND(N), the data structure 
!    defining the triangulation.  Refer to subroutine TRMESH.
!
!    Output, integer ( kind = 4 ) I1, I2, I3, nodal indexes, in counterclockwise
!    order, of the vertices of a triangle containing P if P is contained in a 
!    triangle.  If P is not in the convex hull of the nodes, I1 indexes 
!    the rightmost visible boundary node, I2 indexes the leftmost visible
!    boundary node, and I3 = 0.  Rightmost and leftmost are defined from 
!    the perspective of P, and a pair of points are visible from each 
!    other if and only if the line segment joining them intersects no 
!    triangulation arc.  If P and all of the nodes lie on a common line, 
!    then I1 = I2 = I3 = 0 on output.
!
!  Local parameters:
!
!    B1,B2 =    Unnormalized barycentric coordinates of P with respect 
!               to (N1,N2,N3)
!    IX,IY,IZ = Integer seeds for JRAND
!    LP =       LIST pointer
!    N0,N1,N2 = Nodes in counterclockwise order defining a
!               cone (with vertex N0) containing P
!    N1S,N2S =  Saved values of N1 and N2
!    N3,N4 =    Nodes opposite N1->N2 and N2->N1, respectively
!    NB =       Index of a boundary node -- first neighbor of
!               NF or last neighbor of NL in the boundary traversal loops
!    NF,NL =    First and last neighbors of N0, or first
!               (rightmost) and last (leftmost) nodes
!               visible from P when P is exterior to the triangulation
!    NP,NPP =   Indexes of boundary nodes used in the boundary traversal loops
!    XA,XB,XC = Dummy arguments for FRWRD
!    YA,YB,YC = Dummy arguments for FRWRD
!    XP,YP =    Local variables containing the components of P
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) b1
  real ( kind = 8 ) b2
  logical frwrd
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i3
  integer ( kind = 4 ), save :: ix = 1
  integer ( kind = 4 ), save :: iy = 2
  integer ( kind = 4 ), save :: iz = 3
  integer ( kind = 4 ) jrand
  logical left
  integer ( kind = 4 ) lend(n)
  integer ( kind = 4 ) list(*)
  integer ( kind = 4 ) lp
  integer ( kind = 4 ) lptr(*)
  integer ( kind = 4 ) lstptr
  integer ( kind = 4 ) n0
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n1s
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) n2s
  integer ( kind = 4 ) n3
  integer ( kind = 4 ) n4
  integer ( kind = 4 ) nb
  integer ( kind = 4 ) nf
  integer ( kind = 4 ) nl
  integer ( kind = 4 ) np
  integer ( kind = 4 ) npp
  integer ( kind = 4 ) nst
  real ( kind = 8 ) px
  real ( kind = 8 ) py
  real ( kind = 8 ) store
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xa
  real ( kind = 8 ) xb
  real ( kind = 8 ) xc
  real ( kind = 8 ) xp
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) ya
  real ( kind = 8 ) yb
  real ( kind = 8 ) yc
  real ( kind = 8 ) yp
!
!  Statement function:
!
!  FRWRD = TRUE iff C is forward of A->B iff <A->B,A->C> >= 0.
!
  frwrd(xa,ya,xb,yb,xc,yc) = (xb-xa)*(xc-xa) + (yb-ya)*(yc-ya) >= 0.0D+00
!
!  Initialize variables.
!
  xp = px
  yp = py
  n0 = nst

  if ( n0 < 1  .or.  n < n0 ) then
    n0 = jrand ( n, ix, iy, iz )
  end if
!
!  Set NF and NL to the first and last neighbors of N0, and
!  initialize N1 = NF.
!
1 continue

  lp = lend(n0)
  nl = list(lp)
  lp = lptr(lp)
  nf = list(lp)
  n1 = nf
!
!  Find a pair of adjacent neighbors N1,N2 of N0 that define
!  a wedge containing P:  P LEFT N0->N1 and P RIGHT N0->N2.
!
  if ( 0 < nl ) then
    go to 2
  end if
!
!   N0 is a boundary node.  Test for P exterior.
!
  nl = -nl

  if ( .not. left ( x(n0), y(n0), x(nf), y(nf), xp, yp ) ) then
    nl = n0
    go to 9
  end if

  if ( .not. left(x(nl),y(nl),x(n0),y(n0),xp,yp) ) then
    nb = nf
    nf = n0
    np = nl
    npp = n0
    go to 11
  end if

  go to 3
!
!  N0 is an interior node.  Find N1.
!
2 continue

    do

      if ( left(x(n0),y(n0),x(n1),y(n1),xp,yp) ) then
        exit
      end if

      lp = lptr(lp)
      n1 = list(lp)

      if ( n1 == nl ) then
        go to 6
      end if

    end do
!
!  P is to the left of edge N0->N1.  Initialize N2 to the
!  next neighbor of N0.
!
3 continue

    lp = lptr(lp)
    n2 = abs ( list(lp) )

    if ( .not. left(x(n0),y(n0),x(n2),y(n2),xp,yp) ) then
      go to 7
    end if

    n1 = n2
    if ( n1 /= nl ) then
      go to 3
    end if

  if ( .not. left(x(n0),y(n0),x(nf),y(nf),xp,yp) ) then
    go to 6
  end if

  if (xp == x(n0) .and. yp == y(n0)) then
    go to 5
  end if
!
!  P is left of or on edges N0->NB for all neighbors NB of N0.
!  All points are collinear iff P is left of NB->N0 for
!  all neighbors NB of N0.  Search the neighbors of N0.
!  NOTE: N1 = NL and LP points to NL.
!
4   continue

    if ( .not. left(x(n1),y(n1),x(n0),y(n0),xp,yp) ) then
      go to 5
    end if

    lp = lptr(lp)
    n1 = abs ( list(lp) )

    if ( n1 == nl ) then
      i1 = 0
      i2 = 0
      i3 = 0
      return
    end if

    go to 4
!
!  P is to the right of N1->N0, or P=N0.  Set N0 to N1 and start over.
!
5 continue

  n0 = n1
  go to 1
!
!  P is between edges N0->N1 and N0->NF.
!
6 continue

  n2 = nf
!
!  P is contained in the wedge defined by line segments
!  N0->N1 and N0->N2, where N1 is adjacent to N2.  Set
!  N3 to the node opposite N1->N2, and save N1 and N2 to
!  test for cycling.
!
7 continue

  n3 = n0
  n1s = n1
  n2s = n2
!
!  Top of edge hopping loop.  Test for termination.
!
8 continue

  if ( left ( x(n1), y(n1), x(n2), y(n2), xp, yp ) ) then
!
!  P LEFT N1->N2 and hence P is in (N1,N2,N3) unless an
!  error resulted from floating point inaccuracy and
!  collinearity.  Compute the unnormalized barycentric
!  coordinates of P with respect to (N1,N2,N3).
!
    b1 = (x(n3)-x(n2))*(yp-y(n2)) - (xp-x(n2))*(y(n3)-y(n2))
    b2 = (x(n1)-x(n3))*(yp-y(n3)) - (xp-x(n3))*(y(n1)-y(n3))

    if ( store ( b1 + 1.0D+00 ) >= 1.0D+00  .and. &
         store ( b2 + 1.0D+00 ) >= 1.0D+00 ) then
      go to 16
    end if
!
!  Restart with N0 randomly selected.
!
    n0 = jrand ( n, ix, iy, iz )
    go to 1

  end if
!
!  Set N4 to the neighbor of N2 which follows N1 (node
!  opposite N2->N1) unless N1->N2 is a boundary edge.
!
  lp = lstptr(lend(n2),n1,list,lptr)

  if ( list(lp) < 0 ) then
    nf = n2
    nl = n1
    go to 9
  end if

  lp = lptr(lp)
  n4 = abs ( list(lp) )
!
!  Select the new edge N1->N2 which intersects the line
!  segment N0-P, and set N3 to the node opposite N1->N2.
!
  if ( left(x(n0),y(n0),x(n4),y(n4),xp,yp) ) then
    n3 = n1
    n1 = n4
    n2s = n2
    if (n1 /= n1s  .and.  n1 /= n0) go to 8
  else
    n3 = n2
    n2 = n4
    n1s = n1
    if ( n2 /= n2s  .and.  n2 /= n0 ) then
      go to 8
    end if
  end if
!
!  The starting node N0 or edge N1-N2 was encountered
!  again, implying a cycle (infinite loop).  Restart
!  with N0 randomly selected.
!
  n0 = jrand ( n, ix, iy, iz )
  go to 1
!
!  Boundary traversal loops.  NL->NF is a boundary edge and
!  P RIGHT NL->NF.  Save NL and NF.

9 continue

  np = nl
  npp = nf
!
!  Find the first (rightmost) visible boundary node NF.  NB
!  is set to the first neighbor of NF, and NP is the last neighbor.
!
10 continue

  lp = lend(nf)
  lp = lptr(lp)
  nb = list(lp)

  if ( .not. left(x(nf),y(nf),x(nb),y(nb),xp,yp) ) then
    go to 12
  end if
!
!  P LEFT NF->NB and thus NB is not visible unless an error
!  resulted from floating point inaccuracy and collinear-
!  ity of the 4 points NP, NF, NB, and P.
!
11 continue

  if ( frwrd(x(nf),y(nf),x(np),y(np),xp,yp)  .or. &
       frwrd(x(nf),y(nf),x(np),y(np),x(nb),y(nb)) ) then
    i1 = nf
    go to 13
  end if
!
!  Bottom of loop.
!
12 continue

  np = nf
  nf = nb
  go to 10
!
!  Find the last (leftmost) visible boundary node NL.  NB
!  is set to the last neighbor of NL, and NPP is the first
!  neighbor.
!
13 continue

  lp = lend(nl)
  nb = -list(lp)

  if ( .not. left(x(nb),y(nb),x(nl),y(nl),xp,yp) ) then
    go to 14
  end if
!
!  P LEFT NB->NL and thus NB is not visible unless an error
!  resulted from floating point inaccuracy and collinear-
!  ity of the 4 points P, NB, NL, and NPP.
!
  if ( frwrd(x(nl),y(nl),x(npp),y(npp),xp,yp)  .or. &
       frwrd(x(nl),y(nl),x(npp),y(npp),x(nb),y(nb)) ) then
    go to 15
  end if
!
!  Bottom of loop.
!
14 continue

  npp = nl
  nl = nb
  go to 13
!
!  NL is the leftmost visible boundary node.
!
15 continue

  i2 = nl
  i3 = 0
  return
!
!  P is in the triangle (N1,N2,N3).
!
16 continue
 
  i1 = n1
  i2 = n2
  i3 = n3

  return
end subroutine trfind
