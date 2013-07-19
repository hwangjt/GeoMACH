subroutine addnod ( k, xk, yk, ist, ncc, lcc, n, x, y, list, lptr, lend, lnew, &
  ier )

!*****************************************************************************80
!
!! ADDNOD adds a node to a triangulation.
!
!  Discussion:
!
!    Given a triangulation of N nodes in the plane created by
!    subroutine TRMESH or TRMSHR, this subroutine updates the
!    data structure with the addition of a new node in position
!    K.  If node K is inserted into X and Y (K <= N) rather
!    than appended (K = N+1), then a corresponding insertion
!    must be performed in any additional arrays associated
!    with the nodes.  For example, an array of data values Z
!    must be shifted down to open up position K for the new
!    value:  set Z(I+1) to Z(I) for I = N,N-1,...,K.  For
!    optimal efficiency, new nodes should be appended whenever
!    possible.  Insertion is necessary, however, to add a non-
!    constraint node when constraints are present (refer to
!    subroutine ADDCST).
!
!    Note that a constraint node cannot be added by this
!    routine.  In order to insert a constraint node, it is
!    necessary to add the node with no constraints present
!    (call this routine with NCC = 0), update LCC by increment-
!    ing the appropriate entries, and then create (or restore)
!    the constraints by a call to ADDCST.
!
!    The algorithm consists of the following steps:  node K
!    is located relative to the triangulation (TRFIND), its
!    index is added to the data structure (INTADD or BDYADD),
!    and a sequence of swaps (SWPTST and SWAP) are applied to
!    the arcs opposite K so that all arcs incident on node K
!    and opposite node K (excluding constraint arcs) are local-
!    ly optimal (satisfy the circumcircle test).  Thus, if a
!    (constrained) Delaunay triangulation is input, a (con-
!    strained) Delaunay triangulation will result.  All indexes
!    are incremented as necessary for an insertion.
!
!  Modified:
!
!    29 March 2002
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
!    Input, integer ( kind = 4 ) K, the nodal index (index for X, Y, and LEND)
!    of the new node to be added.  1 <= K <= LCC(1).  (K <= N+1 if NCC=0).
!
!    Input, real ( kind = 8 ) XK, YK, the coordinates of the new node (to be
!    stored in X(K) and Y(K)).  The node must not lie in a constraint region.
!
!    Input, integer ( kind = 4 ) IST, the index of a node at which TRFIND 
!    begins the search.  Search time depends on the proximity
!    of this node to node K.  1 <= IST <= N.
!
!    Input, integer ( kind = 4 ) NCC, the number of constraint curves.  
!    0 <= NCC.
!
!    Input/output, integer ( kind = 4 ) LCC(*), list of constraint curve 
!    starting indexes (or dummy array of length 1 if NCC = 0).  Refer to 
!    subroutine ADDCST.  On output, starting indexes incremented by 1 to 
!    reflect the insertion of K unless NCC = 0 or (IER /= 0 and IER /= -4).
!
!    Input/output, integer ( kind = 4 ) N, the number of nodes in the 
!    triangulation.  3 <= N.  Note that N will be incremented following the 
!    addition of node K.
!
!    Input, real ( kind = 8 ) X(N+1), real Y(N+1), containing the coordinates 
!    of the nodes in the first N positions with non-constraint nodes
!    in the first LCC(1)-1 locations if 0 < NCC.  On output, updated with
!    the insertion of XK and YK in the K-th positions (node I+1 was node 
!    I before the insertion for I = K to N if K <= N)
!    unless IER /= 0 and IER /= -4.
!
!    Input/output, integer ( kind = 4 ) LIST(*), LPTR(*), LEND(N), LNEW, the 
!    data structure associated with the triangulation of nodes 1 to N.  The 
!    arrays must have sufficient length for N+1 nodes.  Refer to TRMESH.
!    On output, updated with the addition of node K unless
!    IER /= 0 and IER /= -4.
!
!    Output, integer ( kind = 4 ) IER = Error indicator:
!     0 if no errors were encountered.
!    -1 if K, IST, NCC, N, or an LCC entry is outside its valid range on input.
!    -2 if all nodes (including K) are collinear.
!     L if nodes L and K coincide for some L.
!    -3 if K lies in a constraint region.
!    -4 if an error flag is returned by SWAP implying that the triangulation
!      (geometry) was bad on input.
!
  implicit none

  logical crtri
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i3
  integer ( kind = 4 ) ibk
  integer ( kind = 4 ) ier
  integer ( kind = 4 ) in1
  integer ( kind = 4 ) indxcc
  integer ( kind = 4 ) io1
  integer ( kind = 4 ) io2
  integer ( kind = 4 ) ist
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kk
  integer ( kind = 4 ) l
  integer ( kind = 4 ) lcc(*)
  integer ( kind = 4 ) lccip1
  integer ( kind = 4 ) lend(*)
  integer ( kind = 4 ) list(*)
  integer ( kind = 4 ) lnew
  integer ( kind = 4 ) lp
  integer ( kind = 4 ) lpf
  integer ( kind = 4 ) lpo1
  integer ( kind = 4 ) lptr(*)
  integer ( kind = 4 ) lstptr
  integer ( kind = 4 ) n
  integer ( kind = 4 ) ncc
  integer ( kind = 4 ) nm1
  logical swptst
  real ( kind = 8 ) x(*)
  real ( kind = 8 ) xk
  real ( kind = 8 ) y(*)
  real ( kind = 8 ) yk

  kk = k
!
!  Test for an invalid input parameter.
!
  if ( kk < 1  .or.  ist < 1  .or.  n < ist &
      .or.  ncc < 0  .or.  n < 3 ) then
    ier = -1
    return
  end if

  lccip1 = n + 1

  do i = ncc, 1, -1
    if ( lccip1-lcc(i) < 3 ) then
      ier = -1
      return
    end if
    lccip1 = lcc(i)
  end do

  if ( lccip1 < kk ) then
    ier = -1
    return
  end if
!
!  Find a triangle (I1,I2,I3) containing K or the rightmost
!  (I1) and leftmost (I2) visible boundary nodes as viewed from node K.
!
  call trfind ( ist, xk, yk, n, x, y, list, lptr, lend, i1, i2, i3 )
!
!  Test for collinear nodes, duplicate nodes, and K lying in
!  a constraint region.
!
  if ( i1 == 0 ) then
    ier = -2
    return
  end if

  if ( i3 /= 0 ) then

    l = i1
    if ( xk == x(l)  .and.  yk == y(l) ) then
      ier = l
      return
    end if

    l = i2
    if ( xk == x(l)  .and.  yk == y(l) ) then
      ier = l
      return
    end if

    l = i3
    if ( xk == x(l)  .and.  yk == y(l) ) then
      ier = l
      return
    end if

    if ( 0 < ncc .and.  crtri(ncc,lcc,i1,i2,i3) ) then
      ier = -3
      return
    end if

  else
!
!  K is outside the convex hull of the nodes and lies in a
!  constraint region iff an exterior constraint curve is present.
!
    if ( 0 < ncc .and. indxcc(ncc,lcc,n,list,lend) /= 0 ) then
      ier = -3
      return
    end if

  end if
!
!  No errors encountered.
!
  ier = 0
  nm1 = n
  n = n + 1

  if (kk < n) then
!
!  Open a slot for K in X, Y, and LEND, and increment all
!  nodal indexes which are greater than or equal to K.
!
!  Note that LIST, LPTR, and LNEW are not yet updated with
!  either the neighbors of K or the edges terminating on K.
!
    do ibk = nm1, kk, -1
      x(ibk+1) = x(ibk)
      y(ibk+1) = y(ibk)
      lend(ibk+1) = lend(ibk)
    end do

    do i = 1, ncc
      lcc(i) = lcc(i) + 1
    end do

    l = lnew - 1

    do i = 1, l

      if ( kk <= list(i) ) then
        list(i) = list(i) + 1
      end if

      if ( list(i) <= -kk ) then
        list(i) = list(i) - 1
      end if

    end do

    if ( kk <= i1 ) then
      i1 = i1 + 1
    end if

    if ( kk <= i2 ) then
      i2 = i2 + 1
    end if

    if ( kk <= i3 ) then
      i3 = i3 + 1
    end if

  end if
!
!  Insert K into X and Y, and update LIST, LPTR, LEND, and
!  LNEW with the arcs containing node K.
!
  x(kk) = xk
  y(kk) = yk

  if ( i3 == 0 ) then
    call bdyadd ( kk, i1, i2, list, lptr, lend, lnew )
  else
    call intadd ( kk, i1, i2, i3, list, lptr, lend, lnew )
  end if
!
!  Initialize variables for optimization of the triangulation.
!
  lp = lend(kk)
  lpf = lptr(lp)
  io2 = list(lpf)
  lpo1 = lptr(lpf)
  io1 = abs ( list(lpo1) )
!
!  Begin loop:  find the node opposite K.
!
  do

    lp = lstptr ( lend(io1), io2, list, lptr )

    if ( 0 <= list(lp) ) then

      lp = lptr(lp)
      in1 = abs ( list(lp) )
!
!  Swap test:  if a swap occurs, two new arcs are
!  opposite K and must be tested.
!
      if ( .not. crtri ( ncc, lcc, io1, io2, in1 ) ) then

        if ( swptst(in1,kk,io1,io2,x,y) ) then

          call swap ( in1, kk, io1, io2, list, lptr, lend, lpo1 )

          if ( lpo1 == 0 ) then
            ier = -4
            exit
          end if
  
          io1 = in1

          cycle

        end if

      end if

    end if
!
!  No swap occurred.  Test for termination and reset IO2 and IO1.
!
    if ( lpo1 == lpf .or. list(lpo1) < 0 ) then
      exit
    end if

    io2 = io1
    lpo1 = lptr(lpo1)
    io1 = abs ( list(lpo1) )

  end do

  return

end subroutine addnod




subroutine bdyadd ( kk, i1, i2, list, lptr, lend, lnew )

!*****************************************************************************80
!
!! BDYADD adds a boundary node to a triangulation.
!
!  Discussion:
!
!    This subroutine adds a boundary node to a triangulation
!    of a set of points in the plane.  The data structure is
!    updated with the insertion of node KK, but no optimization
!    is performed.
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
!    Input, integer ( kind = 4 ) KK, the index of a node to be connected to the
!    sequence of all visible boundary nodes.  1 <= KK and
!    KK must not be equal to I1 or I2.
!
!    Input, integer ( kind = 4 ) I1, the first (rightmost as viewed from KK) 
!    boundary node in the triangulation which is visible from
!    node KK (the line segment KK-I1 intersects no arcs.
!
!    Input, integer ( kind = 4 ) I2, the last (leftmost) boundary node which is
!    visible from node KK.  I1 and I2 may be determined by subroutine TRFIND.
!
!    Input/output, integer ( kind = 4 ) LIST(*), LPTR(*), LEND(N), LNEW.  The 
!    triangulation data structure created by TRMESH or TRMSHR.
!    On input, nodes I1 and I2 must be included in the triangulation.
!    On output, the data structure has been updated with the addition 
!    of node KK.  Node KK is connected to I1, I2, and all boundary 
!    nodes in between.
!
  implicit none

  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kk
  integer ( kind = 4 ) lend(*)
  integer ( kind = 4 ) list(*)
  integer ( kind = 4 ) lnew
  integer ( kind = 4 ) lp
  integer ( kind = 4 ) lptr(*)
  integer ( kind = 4 ) lsav
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) next
  integer ( kind = 4 ) nsav

  k = kk
  n1 = i1
  n2 = i2
!
!  Add K as the last neighbor of N1.
!
  lp = lend(n1)
  lsav = lptr(lp)
  lptr(lp) = lnew
  list(lnew) = -k
  lptr(lnew) = lsav
  lend(n1) = lnew
  lnew = lnew + 1
  next = -list(lp)
  list(lp) = next
  nsav = next
!
!  Loop on the remaining boundary nodes between N1 and N2,
!  adding K as the first neighbor.
!
  do

    lp = lend(next)

    call insert ( k, lp, list, lptr, lnew )

    if ( next == n2 ) then
      exit
    end if

    next = -list(lp)
    list(lp) = next

  end do
!
!  Add the boundary nodes between N1 and N2 as neighbors
!  of node K.
!
  lsav = lnew
  list(lnew) = n1
  lptr(lnew) = lnew + 1
  lnew = lnew + 1
  next = nsav

  do

    if ( next == n2 ) then
      exit
    end if

    list(lnew) = next
    lptr(lnew) = lnew + 1
    lnew = lnew + 1
    lp = lend(next)
    next = list(lp)

  end do

  list(lnew) = -n2
  lptr(lnew) = lsav
  lend(k) = lnew
  lnew = lnew + 1

  return

end subroutine bdyadd




function crtri ( ncc, lcc, i1, i2, i3 )

!*****************************************************************************80
!
!! CRTRI determines if a triangle lies in a constraint region.
!
!  Discussion:
!
!    This function returns TRUE if and only if triangle (I1,
!    I2,I3) lies in a constraint region.
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
!    Input, integer ( kind = 4 ) NCC, LCC(*), onstraint data structure.  
!    Refer to ADDCST.
!
!    Input, integer ( kind = 4 ) I1, I2, I3, Nodal indexes of the counterclockwise-
!    ordered vertices of a triangle.
!
!    Output, logical CRTRI, is TRUE if and only if (I1,I2,I3) is a 
!    constraint region triangle.
!
  implicit none

  logical crtri
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i3
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) imin
  integer ( kind = 4 ) lcc(*)
  integer ( kind = 4 ) ncc

  imax = max ( i1, i2, i3 )
!
!  Find the index I of the constraint containing IMAX.
!
  i = ncc + 1

  do

    i = i - 1

    if ( i <= 0 ) then
      crtri = .false.
      return
    end if

    if ( lcc(i) <= imax ) then
      exit
    end if

  end do

  imin = min ( i1, i2, i3 )
!
!  P lies in a constraint region iff I1, I2, and I3 are nodes
!  of the same constraint (LCC(I) <= IMIN), and (IMIN,IMAX)
!  is (I1,I3), (I2,I1), or (I3,I2).
!
  crtri = lcc(i) <= imin .and. ( &
    ( imin == i1 .and. imax == i3 ) .or.  &
    ( imin == i2 .and. imax == i1 ) .or.  &
    ( imin == i3 .and. imax == i2 ) )

  return

end function crtri




function indxcc ( ncc, lcc, n, list, lend )

!*****************************************************************************80
!
!! INDXCC returns the index of an exterior constraint curve.
!
!  Discussion:
!
!    Given a constrained Delaunay triangulation, this 
!    function returns the index, if any, of an exterior constraint
!    curve (an unbounded constraint region).  An exterior 
!    constraint curve is assumed to be present if and only if the
!    clockwise-ordered sequence of boundary nodes is a 
!    subsequence of a constraint node sequence.  The triangulation
!    adjacencies corresponding to constraint edges may or may
!    not have been forced by a call to ADDCST, and the 
!    constraint region may or may not be valid (contain no nodes).
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
!    Input, integer ( kind = 4 ) NCC, the number of constraints.  NCC >= 0.
!
!    Input, integer ( kind = 4 ) LCC(*), list of constraint curve starting 
!    indexes (or dummy array of length 1 if NCC = 0).  Refer to subroutine 
!    ADDCST.
!
!    Input, integer ( kind = 4 ) N, the number of nodes in the triangulation.  
!    N >= 3.
!
!    Input, integer ( kind = 4 ) LIST(*), LEND(N), the data structure defining
!    the triangulation.  Refer to subroutine TRMESH.
!
!    Output, integer ( kind = 4 ) INDXCC, index of the exterior constraint 
!    curve, if present, or 0 otherwise.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ifrst
  integer ( kind = 4 ) ilast
  integer ( kind = 4 ) indxcc
  integer ( kind = 4 ) lcc(*)
  integer ( kind = 4 ) lend(n)
  integer ( kind = 4 ) list(*)
  integer ( kind = 4 ) lp
  integer ( kind = 4 ) n0
  integer ( kind = 4 ) ncc
  integer ( kind = 4 ) nst
  integer ( kind = 4 ) nxt

  indxcc = 0

  if ( ncc < 1 ) then
    return
  end if
!
!  Set N0 to the boundary node with smallest index.
!
  n0 = 0

  do

    n0 = n0 + 1
    lp = lend(n0)

    if ( list(lp) <= 0 ) then
      exit
    end if

  end do
!
!  Search in reverse order for the constraint I, if any, that
!  contains N0.  IFRST and ILAST index the first and last
!  nodes in constraint I.
!
  i = ncc
  ilast = n

  do

    ifrst = lcc(i)

    if ( ifrst <= n0 ) then
      exit
    end if

    if ( i == 1 ) then
      return
    end if

    i = i - 1
    ilast = ifrst - 1

  end do
!
!  N0 is in constraint I which indexes an exterior constraint
!  curve iff the clockwise-ordered sequence of boundary
!  node indexes beginning with N0 is increasing and bounded
!  above by ILAST.
!
  nst = n0

  do

    nxt = -list(lp)

    if ( nxt == nst ) then
      exit
    end if

    if ( nxt <= n0  .or. ilast < nxt ) then
      return
    end if

    n0 = nxt
    lp = lend(n0)

  end do
!
!  Constraint I contains the boundary node sequence as a subset.
!
  indxcc = i

  return

end function indxcc




subroutine intadd ( kk, i1, i2, i3, list, lptr, lend, lnew )

!*****************************************************************************80
!
!! INTADD adds an interior point to a triangulation.
!
!  Discussion:
!
!    This subroutine adds an interior node to a triangulation
!    of a set of points in the plane.  The data structure is
!    updated with the insertion of node KK into the triangle
!    whose vertices are I1, I2, and I3.  No optimization of the
!    triangulation is performed.
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
!    Input, integer ( kind = 4 ) KK, the index of the node to be inserted.  
!    1 <= KK and KK must not be equal to I1, I2, or I3.
!
!    Input, integer ( kind = 4 ) I1, I2, I3, indexes of the 
!    counterclockwise-ordered sequence of vertices of a triangle which 
!    contains node KK.
!
!    Input/output, integer ( kind = 4 ) LIST(*), LPTR(*), LEND(N), LNEW, the
!    data structure defining the triangulation.  Refer to subroutine TRMESH. 
!    Triangle (I1,I2,I3) must be included in the triangulation.
!    On output, updated with the addition of node KK.  KK
!    will be connected to nodes I1, I2, and I3.
!
  implicit none

  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i3
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kk
  integer ( kind = 4 ) lend(*)
  integer ( kind = 4 ) list(*)
  integer ( kind = 4 ) lnew
  integer ( kind = 4 ) lp
  integer ( kind = 4 ) lptr(*)
  integer ( kind = 4 ) lstptr
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) n3

  k = kk
!
!  Initialization.
!
  n1 = i1
  n2 = i2
  n3 = i3
!
!  Add K as a neighbor of I1, I2, and I3.
!
  lp = lstptr(lend(n1),n2,list,lptr)
  call insert (k,lp,list,lptr,lnew)
  lp = lstptr(lend(n2),n3,list,lptr)
  call insert (k,lp,list,lptr,lnew)
  lp = lstptr(lend(n3),n1,list,lptr)
  call insert (k,lp,list,lptr,lnew)
!
!  Add I1, I2, and I3 as neighbors of K.
!
  list(lnew) = n1
  list(lnew+1) = n2
  list(lnew+2) = n3
  lptr(lnew) = lnew + 1
  lptr(lnew+1) = lnew + 2
  lptr(lnew+2) = lnew
  lend(k) = lnew + 2
  lnew = lnew + 3

  return

end subroutine intadd
