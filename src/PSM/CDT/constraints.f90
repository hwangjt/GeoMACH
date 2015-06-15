subroutine edge ( in1, in2, x, y, lwk, iwk, list, lptr, lend, ier )

!*****************************************************************************80
!
!! EDGE swaps arcs to force two nodes to be adjacent.
!
!  Discussion:
!
!    Given a triangulation of N nodes and a pair of nodal
!    indexes IN1 and IN2, this routine swaps arcs as necessary
!    to force IN1 and IN2 to be adjacent.  Only arcs which
!    intersect IN1-IN2 are swapped out.  If a Delaunay triangu-
!    lation is input, the resulting triangulation is as close
!    as possible to a Delaunay triangulation in the sense that
!    all arcs other than IN1-IN2 are locally optimal.
!
!    A sequence of calls to EDGE may be used to force the
!    presence of a set of edges defining the boundary of a non-
!    convex and/or multiply connected region (refer to Subrou-
!    tine ADDCST), or to introduce barriers into the triangula-
!    tion.  Note that Subroutine GETNP will not necessarily
!    return closest nodes if the triangulation has been con-
!    strained by a call to EDGE.  However, this is appropriate
!    in some applications, such as triangle-based interpolation
!    on a nonconvex domain.
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
!    Input, integer ( kind = 4 ) IN1, IN2, indexes (of X and Y) in the range
!    1 to N defining a pair of nodes to be connected by an arc.
!
!    Input, real ( kind = 8 ) X(N), Y(N), the coordinates of the nodes.
!
!    Input/output, integer ( kind = 4 ) LWK.  On input, the number of columns 
!    reserved for IWK.  This must be at least NI, the number of arcs which 
!    intersect IN1-IN2.  (NI is bounded by N-3.)  On output, the number of arcs
!    which intersect IN1-IN2 (but not more than the input value of LWK) unless
!    IER = 1 or IER = 3.  LWK = 0 if and only if IN1 and IN2 were adjacent 
!    (or LWK=0) on input.
!
!    Output, integer ( kind = 4 ) IWK(2*LWK), the indexes of the endpoints of 
!    the new arcs other than IN1-IN2 unless IER > 0 or LWK = 0.  New arcs to 
!    the left of IN2-IN1 are stored in the first K-1 columns (left portion of 
!    IWK), column K contains zeros, and new arcs to the right of IN2-IN1
!    occupy columns K+1,...,LWK.  (K can be determined by searching IWK 
!    for the zeros.)
!
!    Input, integer ( kind = 4 ) LIST(*), LPTR(*), LEND(N), the data structure
!    defining the triangulation.  Refer to subroutine TRMESH.  On output, 
!    updated if necessary to reflect the presence of an arc connecting IN1 and 
!    IN2 unless IER /= 0.  The data structure has been altered if IER = 4.
!
!    Output, integer ( kind = 4 ) IER = Error indicator:
!    0, if no errors were encountered.
!    1, if IN1 < 1, IN2 .LT. 1, IN1 = IN2, or LWK < 0 on input.
!    2, if more space is required in IWK.
!    3, if IN1 and IN2 could not be connected due to either an invalid 
!      data structure or collinear nodes (and floating point error).
!    4, if an error flag was returned by OPTIM.
!
!  Local parameters:
!
!    DX,DY =   Components of arc N1-N2.
!
!    I =       DO-loop index and column index for IWK
!    IERR =    Error flag returned by Subroutine OPTIM
!    IWC =     IWK index between IWF and IWL -- NL->NR is
!              stored in IWK(1,IWC)->IWK(2,IWC)
!    IWCP1 =   IWC + 1
!    IWEND =   Input or output value of LWK
!    IWF =     IWK (column) index of the first (leftmost) arc
!              which intersects IN1->IN2
!    IWL =     IWK (column) index of the last (rightmost) are
!              which intersects IN1->IN2
!    LFT =     Flag used to determine if a swap results in the
!              new arc intersecting IN1-IN2 -- LFT = 0 iff
!              N0 = IN1, LFT = -1 implies N0 LEFT IN1->IN2,
!              and LFT = 1 implies N0 LEFT IN2->IN1
!    LP21 =    Unused parameter returned by SWAP
!    LP =      List pointer (index) for LIST and LPTR
!    LPL =     Pointer to the last neighbor of IN1 or NL
!    N0 =      Neighbor of N1 or node opposite NR->NL
!    N1,N2 =   Local copies of IN1 and IN2
!    N1FRST =  First neighbor of IN1
!    N1LST =   (Signed) last neighbor of IN1
!    NEXT =    Node opposite NL->NR
!    NIT =     Flag or number of iterations employed by OPTIM
!    NL,NR =   Endpoints of an arc which intersects IN1-IN2
!              with NL LEFT IN1->IN2
!    X0,Y0 =   Coordinates of N0
!    X1,Y1 =   Coordinates of IN1
!    X2,Y2 =   Coordinates of IN2
!
  implicit none

  real ( kind = 8 ) dx
  real ( kind = 8 ) dy
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ier
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) in1
  integer ( kind = 4 ) in2
  integer ( kind = 4 ) iwc
  integer ( kind = 4 ) iwcp1
  integer ( kind = 4 ) iwend
  integer ( kind = 4 ) iwf
  integer ( kind = 4 ) iwk(2,*)
  integer ( kind = 4 ) iwl
  logical left
  integer ( kind = 4 ) lend(*)
  integer ( kind = 4 ) lft
  integer ( kind = 4 ) list(*)
  integer ( kind = 4 ) lp
  integer ( kind = 4 ) lp21
  integer ( kind = 4 ) lpl
  integer ( kind = 4 ) lptr(*)
  integer ( kind = 4 ) lwk
  integer ( kind = 4 ) n0
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n1frst
  integer ( kind = 4 ) n1lst
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) next
  integer ( kind = 4 ) nit
  integer ( kind = 4 ) nl
  integer ( kind = 4 ) nr
  real ( kind = 8 ) x(*)
  real ( kind = 8 ) x0
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) y(*)
  real ( kind = 8 ) y0
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2
!
!  Store IN1, IN2, and LWK in local variables and test for errors.
!
  n1 = in1
  n2 = in2
  iwend = lwk

  if ( n1 < 1  .or.  n2 < 1  .or.  n1 == n2  .or. iwend < 0 ) then
    ier = 1
    return
  end if
!
!  Test for N2 as a neighbor of N1.  LPL points to the last neighbor of N1.
!
  lpl = lend(n1)
  n0 = abs ( list(lpl) )
  lp = lpl

  do
!
!  IN1 and IN2 were adjacent on input.
!
    if ( n0 == n2 ) then
      ier = 0
      return
    end if

    lp = lptr(lp)
    n0 = list(lp)

    if ( lp == lpl ) then
      exit
    end if

  end do
!
!  Initialize parameters.
!
  iwl = 0
  nit = 0
!
!  Store the coordinates of N1 and N2.
!
2 continue

  x1 = x(n1)
  y1 = y(n1)
  x2 = x(n2)
  y2 = y(n2)
!
!  Set NR and NL to adjacent neighbors of N1 such that
!  NR LEFT N2->N1 and NL LEFT N1->N2,
!  (NR Forward N1->N2 or NL Forward N1->N2), and
!  (NR Forward N2->N1 or NL Forward N2->N1).
!
!  Initialization:  Set N1FRST and N1LST to the first and
!  (signed) last neighbors of N1, respectively, and
!  initialize NL to N1FRST.
!
  lpl = lend(n1)
  n1lst = list(lpl)
  lp = lptr(lpl)
  n1frst = list(lp)
  nl = n1frst

  if ( n1lst < 0 ) then
    go to 4
  end if
!
!  N1 is an interior node.  Set NL to the first candidate
!  for NR (NL LEFT N2->N1).
!
3   continue

    if ( left(x2,y2,x1,y1,x(nl),y(nl)) ) then
      go to 4
    end if

    lp = lptr(lp)
    nl = list(lp)

    if ( nl /= n1frst ) then
      go to 3
    end if
!
!  All neighbors of N1 are strictly left of N1->N2.
!
  go to 5
!
!  NL = LIST(LP) LEFT N2->N1.  Set NR to NL and NL to the
!  following neighbor of N1.
!
4   continue

    nr = nl
    lp = lptr(lp)
    nl = abs ( list(lp) ) 
!
!  NL LEFT N1->N2 and NR LEFT N2->N1.  The Forward tests
!  are employed to avoid an error associated with
!  collinear nodes.
!
    if ( left(x1,y1,x2,y2,x(nl),y(nl)) ) then

      dx = x2-x1
      dy = y2-y1
      if ((dx*(x(nl)-x1)+dy*(y(nl)-y1) >= 0.  .or. &
           dx*(x(nr)-x1)+dy*(y(nr)-y1) >= 0.)  .and. &
          (dx*(x(nl)-x2)+dy*(y(nl)-y2) <= 0.  .or. &
           dx*(x(nr)-x2)+dy*(y(nr)-y2) <= 0.)) go to 6
!
!  NL-NR does not intersect N1-N2.  However, there is
!  another candidate for the first arc if NL lies on
!  the line N1-N2.
!
      if ( .not. left(x2,y2,x1,y1,x(nl),y(nl)) ) go to 5
    end if
!
!  Bottom of loop.
!
    if ( nl /= n1frst ) then
      go to 4
    end if
!
!  Either the triangulation is invalid or N1-N2 lies on the
!  convex hull boundary and an edge NR->NL (opposite N1 and
!  intersecting N1-N2) was not found due to floating point
!  error.  Try interchanging N1 and N2 -- NIT > 0 iff this
!  has already been done.
!
5 continue

  if ( 0 < nit ) then
    go to 33
  end if

  nit = 1
  n1 = n2
  n2 = in1
  go to 2
!
!  Store the ordered sequence of intersecting edges NL->NR in
!  IWK(1,IWL)->IWK(2,IWL).
!
6 continue

  iwl = iwl + 1

  if (iwl > iwend) then
     print *, iwl, iwend
    ier = 2
    return
  end if

  iwk(1,iwl) = nl
  iwk(2,iwl) = nr
!
!  Set NEXT to the neighbor of NL which follows NR.
!
  lpl = lend(nl)
  lp = lptr(lpl)
!
!  Find NR as a neighbor of NL.  The search begins with
!  the first neighbor.
!
7   continue

    if (list(lp) == nr) go to 8
    lp = lptr(lp)
    if (lp /= lpl) go to 7
!
!  NR must be the last neighbor, and NL->NR cannot be a
!  boundary edge.
!
  if (list(lp) /= nr) then
    go to 33
  end if
!
!  Set NEXT to the neighbor following NR, and test for
!  termination of the store loop.
!
8 continue

  lp = lptr(lp)
  next = abs ( list(lp) )

  if (next == n2) then
    go to 9
  end if
!
!  Set NL or NR to NEXT.
!
  if ( left(x1,y1,x2,y2,x(next),y(next)) ) then
    nl = next
  else
    nr = next
  end if

  go to 6
!
!  IWL is the number of arcs which intersect N1-N2.
!  Store LWK.
!
9 continue

  lwk = iwl
  iwend = iwl
!
!  Initialize for edge swapping loop -- all possible swaps
!  are applied (even if the new arc again intersects
!  N1-N2), arcs to the left of N1->N2 are stored in the
!  left portion of IWK, and arcs to the right are stored in
!  the right portion.  IWF and IWL index the first and last
!  intersecting arcs.
!
  iwf = 1
!
!  Top of loop -- set N0 to N1 and NL->NR to the first edge.
!  IWC points to the arc currently being processed.  LFT
!  <= 0 iff N0 LEFT N1->N2.
!
10 continue

  lft = 0
  n0 = n1
  x0 = x1
  y0 = y1
  nl = iwk(1,iwf)
  nr = iwk(2,iwf)
  iwc = iwf
!
!  Set NEXT to the node opposite NL->NR unless IWC is the last arc.
!
11 continue

  if (iwc == iwl) go to 21
  iwcp1 = iwc + 1
  next = iwk(1,iwcp1)
  if (next /= nl) go to 16
  next = iwk(2,iwcp1)
!
!  NEXT RIGHT N1->N2 and IWC < IWL.  Test for a possible swap.
!
  if ( .not. left(x0,y0,x(nr),y(nr),x(next),y(next)) ) then
    go to 14
  end if

  if (lft >= 0) then
    go to 12
  end if

  if ( .not. left(x(nl),y(nl),x0,y0,x(next),y(next)) ) then
    go to 14
  end if
!
!  Replace NL->NR with N0->NEXT.
!
  call swap (next,n0,nl,nr, list,lptr,lend, lp21)
  iwk(1,iwc) = n0
  iwk(2,iwc) = next
  go to 15
!
!  Swap NL-NR for N0-NEXT, shift columns IWC+1,...,IWL to
!  the left, and store N0-NEXT in the right portion of IWK.
!
12 continue

  call swap (next,n0,nl,nr, list,lptr,lend, lp21)

  do i = iwcp1,iwl
    iwk(1,i-1) = iwk(1,i)
    iwk(2,i-1) = iwk(2,i)
  end do

  iwk(1,iwl) = n0
  iwk(2,iwl) = next
  iwl = iwl - 1
  nr = next
  go to 11
!
!  A swap is not possible.  Set N0 to NR.
!
14 continue

  n0 = nr
  x0 = x(n0)
  y0 = y(n0)
  lft = 1
!
!  Advance to the next arc.
!
15 continue

  nr = next
  iwc = iwc + 1
  go to 11
!
!  NEXT LEFT N1->N2, NEXT /= N2, and IWC < IWL.
!  Test for a possible swap.
!
16 continue

  if ( .not. left(x(nl),y(nl),x0,y0,x(next),y(next)) ) then
    go to 19
  end if

  if (lft <= 0) then
    go to 17
  end if

  if ( .not. left(x0,y0,x(nr),y(nr),x(next),y(next)) ) then
    go to 19
  end if
!
!  Replace NL->NR with NEXT->N0.
!
  call swap (next,n0,nl,nr, list,lptr,lend, lp21)
  iwk(1,iwc) = next
  iwk(2,iwc) = n0
  go to 20
!
!  Swap NL-NR for N0-NEXT, shift columns IWF,...,IWC-1 to
!  the right, and store N0-NEXT in the left portion of IWK.
!
17 continue

  call swap (next,n0,nl,nr, list,lptr,lend, lp21)

  do i = iwc-1,iwf,-1
    iwk(1,i+1) = iwk(1,i)
    iwk(2,i+1) = iwk(2,i)
  end do

  iwk(1,iwf) = n0
  iwk(2,iwf) = next
  iwf = iwf + 1
  go to 20
!
!  A swap is not possible.  Set N0 to NL.
!
19 continue

  n0 = nl
  x0 = x(n0)
  y0 = y(n0)
  lft = -1
!
!  Advance to the next arc.
!
20 continue

  nl = next
  iwc = iwc + 1
  go to 11
!
!  N2 is opposite NL->NR (IWC = IWL).
!
21 continue

  if (n0 == n1) go to 24
  if (lft < 0) go to 22
!
!  N0 RIGHT N1->N2.  Test for a possible swap.
!
  if ( .not. left(x0,y0,x(nr),y(nr),x2,y2) ) go to 10
!
!  Swap NL-NR for N0-N2 and store N0-N2 in the right
!  portion of IWK.
!
  call swap (n2,n0,nl,nr, list,lptr,lend, lp21)
  iwk(1,iwl) = n0
  iwk(2,iwl) = n2
  iwl = iwl - 1
  go to 10
!
!  N0 LEFT N1->N2.  Test for a possible swap.
!
22 continue

  if ( .not. left(x(nl),y(nl),x0,y0,x2,y2) ) then
    go to 10
  end if
!
!  Swap NL-NR for N0-N2, shift columns IWF,...,IWL-1 to the
!  right, and store N0-N2 in the left portion of IWK.
!
  call swap ( n2, n0, nl, nr, list, lptr, lend, lp21 )
  i = iwl

23 continue

  iwk(1,i) = iwk(1,i-1)
  iwk(2,i) = iwk(2,i-1)
  i = i - 1

  if (i > iwf) then
    go to 23
  end if

  iwk(1,iwf) = n0
  iwk(2,iwf) = n2
  iwf = iwf + 1
  go to 10
!
!  IWF = IWC = IWL.  Swap out the last arc for N1-N2 and
!  store zeros in IWK.
!
24 continue

  call swap (n2,n1,nl,nr, list,lptr,lend, lp21)
  iwk(1,iwc) = 0
  iwk(2,iwc) = 0
!
!  Optimization procedure.
!
  if ( 1 < iwc ) then
!
!  Optimize the set of new arcs to the left of IN1->IN2.
!
    nit = 3*(iwc-1)
    call optim ( x, y, iwc-1, list, lptr, lend, nit, iwk, ierr )
    if ( ierr /= 0 ) then
      go to 34
    end if
  end if

  if ( iwc < iwend ) then
!
!  Optimize the set of new arcs to the right of IN1->IN2.
!
    nit = 3*(iwend-iwc)
    call optim ( x, y, iwend-iwc, list, lptr, lend, nit, &
      iwk(1,iwc+1), ierr )

    if (ierr /= 0) then
      go to 34
    end if

  end if
!
!  Successful termination.
!
  ier = 0
  return
!
!  Invalid triangulation data structure or collinear nodes
!  on convex hull boundary.
!
   33 ier = 3
  write (*,130) in1, in2
  130 format (//5x,'*** error in edge:  invalid triangula', &
          'tion or null triangles on boundary'/ &
          9x,'in1 =',i4,', in2=',i4/)
  return
!
! Error flag returned by OPTIM.
!
   34 ier = 4
  write (*,140) nit, ierr
  140 format (//5x,'*** error in optim:  nit = ',i4, &
          ', ier = ',i1,' ***'/)
  return

end subroutine edge




subroutine addcst ( ncc, lcc, n, x, y, lwk, iwk, list, lptr, lend, ier )

!*****************************************************************************80
!
!! ADDCST adds constraint curves to a Delaunay triangulation.
!
!  Discussion:
!
!    This subroutine provides for creation of a constrained
!    Delaunay triangulation which, in some sense, covers an
!    arbitrary connected region R rather than the convex hull
!    of the nodes.  This is achieved simply by forcing the
!    presence of certain adjacencies (triangulation arcs) 
!    corresponding to constraint curves.  The union of triangles
!    coincides with the convex hull of the nodes, but triangles
!    in R can be distinguished from those outside of R.  The
!    only modification required to generalize the definition of
!    the Delaunay triangulation is replacement of property 5
!    (refer to TRMESH) by the following:
!
!    5')  If a node is contained in the interior of the 
!         circumcircle of a triangle, then every interior point
!         of the triangle is separated from the node by a
!         constraint arc.
!
!    In order to be explicit, we make the following definitions.  
!    A constraint region is the open interior of a
!    simple closed positively oriented polygonal curve defined
!    by an ordered sequence of three or more distinct nodes
!    (constraint nodes) P(1),P(2),...,P(K), such that P(I) is
!    adjacent to P(I+1) for I = 1,...,K with P(K+1) = P(1).
!    Thus, the constraint region is on the left (and may have
!    nonfinite area) as the sequence of constraint nodes is
!    traversed in the specified order.  The constraint regions
!    must not contain nodes and must not overlap.  The region
!    R is the convex hull of the nodes with constraint regions
!    excluded.
!
!    Note that the terms boundary node and boundary arc are
!    reserved for nodes and arcs on the boundary of the convex
!    hull of the nodes.
!
!    The algorithm is as follows:  given a triangulation
!    which includes one or more sets of constraint nodes, the
!    corresponding adjacencies (constraint arcs) are forced to
!    be present (Subroutine EDGE).  Any additional new arcs
!    required are chosen to be locally optimal (satisfy the
!    modified circumcircle property).
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
!    Input, integer ( kind = 4 ) NCC, the number of constraint curves 
!    (constraint regions).  0 <= NCC.
!
!    Input, integer ( kind = 4 ) LCC(NCC) (or dummy array of length 1 if 
!    NCC = 0) containing the index (for X, Y, and LEND) of the first node of 
!    constraint I in LCC(I) for I = 1 to NCC.  Thus, constraint I
!    contains K = LCC(I+1) - LCC(I) nodes, K >= 3, stored in (X,Y) 
!    locations LCC(I), ..., LCC(I+1)-1, where LCC(NCC+1) = N+1.
!
!    Input, integer ( kind = 4 ) N, the number of nodes in the triangulation, 
!    including constraint nodes.  3 <= N.
!
!    Input, real ( kind = 8 ) X(N), Y(N), the coordinates of the nodes with
!    non-constraint nodes in the first LCC(1)-1 locations, followed by NCC
!    sequences of constraint nodes.  Only one of these sequences may be
!    specified in clockwise order to represent an exterior constraint curve (a 
!    constraint region with nonfinite area).
!
!    Input/output, integer ( kind = 4 ) LWK.  On input, the length of IWK.  
!    This must be at least 2*NI where NI is the maximum number of arcs which 
!    intersect a constraint arc to be added.  NI is bounded by N-3.  On output,
!    the required length of IWK unless IER = 1 or IER = 3.  In the case of 
!    IER = 1, LWK is not altered from its input value.
!
!    Output, integer ( kind = 4 ) IWK(LWK), the endpoint indexes of the new 
!    arcs which were swapped in by the last call to subroutine EDGE.
!
!    Input/output, integer ( kind = 4 ) LIST(*), LPTR(*), LEND(N), the data 
!    structure defining the triangulation.  Refer to subroutine TRMESH.  On 
!    output, the structure has all constraint arcs present unless IER /= 0.  
!    These arrays are not altered if IER = 1.
!
!    Output, integer ( kind = 4 ) IER = Error indicator:
!    0 if no errors were encountered.
!    1 if NCC, N, or an LCC entry is outside its valid range, or 
!      LWK < 0 on input.
!    2 if more space is required in IWK.
!    3 if the triangulation data structure is invalid, or failure (in 
!      EDGE or OPTIM) was caused by collinear nodes on the convex hull 
!      boundary.  An error message is written to logical unit 6 in
!      this case.
!    4 if intersecting constraint arcs were encountered.
!    5 if a constraint region contains a node.
!
  implicit none

  integer ( kind = 4 ) lwk
  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ier
  integer ( kind = 4 ) ifrst
  integer ( kind = 4 ) ilast
  integer ( kind = 4 ) iwk(lwk)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kbak
  integer ( kind = 4 ) kfor
  integer ( kind = 4 ) kn
  integer ( kind = 4 ) lcc(*)
  integer ( kind = 4 ) lccip1
  integer ( kind = 4 ) lend(n)
  integer ( kind = 4 ) list(*)
  integer ( kind = 4 ) lp
  integer ( kind = 4 ) lpb
  integer ( kind = 4 ) lpf
  integer ( kind = 4 ) lpl
  integer ( kind = 4 ) lptr(*)
  integer ( kind = 4 ) lw
  integer ( kind = 4 ) lwd2
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) ncc
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)

  lwd2 = lwk / 2
!
!  Test for errors in input parameters.
!
  ier = 0

  if ( ncc < 0 .or. lwk < 0 ) then
    ier = 1
    return
  end if

  if ( ncc == 0 ) then

    if ( n < 3 ) then
      ier = 1
    else
      ier = 0
      lwk = 0
    end if

    return

  else

    lccip1 = n + 1

    do i = ncc, 1, -1

      if ( lccip1 - lcc(i) < 3 ) then
        ier = 1
        return
      end if

      lccip1 = lcc(i)

    end do

    if ( lccip1 < 1 ) then
      ier = 1
      return
    end if

  end if
!
!  Force the presence of constraint arcs.  The outer loop is
!  on constraints in reverse order.  IFRST and ILAST are
!  the first and last nodes of constraint I.
!
  lwk = 0
  ifrst = n + 1

  do i = ncc, 1, -1

    ilast = ifrst - 1
    ifrst = lcc(i)
!
!  Inner loop on constraint arcs N1-N2 in constraint I.
!
    n1 = ilast

    do n2 = ifrst, ilast

      lw = lwd2

      call edge ( n1, n2, x, y, lw, iwk, list, lptr, lend, ier )

      lwk = max ( lwk, 2 * lw )

      if ( ier == 4 ) then
        ier = 3
      end if

      if ( ier /= 0 ) then
        return
      end if

      n1 = n2
    end do

  end do
!
!  Test for errors.  The outer loop is on constraint I with
!  first and last nodes IFRST and ILAST, and the inner loop
!  is on constraint nodes K with (KBAK,K,KFOR) a subsequence 
!  of constraint I.
!
  ifrst = n + 1

  do i = ncc, 1, -1

    ilast = ifrst - 1
    ifrst = lcc(i)
    kbak = ilast

    do k = ifrst,ilast

      kfor = k + 1

      if ( k == ilast ) then
        kfor = ifrst
      end if
!
!  Find the LIST pointers LPF and LPB of KFOR and KBAK as neighbors of K.
!
      lpf = 0
      lpb = 0
      lpl = lend(k)
      lp = lpl

      do

        lp = lptr(lp)
        kn = abs ( list(lp) )

        if ( kn == kfor ) then
          lpf = lp
        end if

        if ( kn == kbak ) then
          lpb = lp
        end if

        if ( lp == lpl ) then
          exit
        end if

      end do
!
!  A pair of intersecting constraint arcs was encountered
!  if and only if a constraint arc is missing (introduction 
!  of the second caused the first to be swapped out).
!
      if ( lpf == 0 .or. lpb == 0 ) then
        ier = 4
        return
      end if
!
!  Loop on neighbors KN of node K which follow KFOR and
!  precede KBAK.  The constraint region contains no nodes
!  if and only if all such nodes KN are in constraint I.
!
      lp = lpf

      do

        lp = lptr(lp)

        if ( lp == lpb ) then
          exit
        end if

        kn = abs ( list(lp) )

        if ( kn < ifrst .or. ilast < kn ) then
          ier = 5
          return
        end if

      end do

      kbak = k

    end do

  end do

  ier = 0

  return

end subroutine addcst
