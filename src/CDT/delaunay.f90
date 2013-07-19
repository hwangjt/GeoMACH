subroutine trmesh ( n, x, y, list, lptr, lend, lnew, near, next, dist, ier )

!*****************************************************************************80
!
!! TRMESH triangulates a set of points in the plane.
!
!  Discussion:
!
!    This subroutine creates a Delaunay triangulation of a
!    set of N arbitrarily distributed points in the plane
!    referred to as nodes.  The Delaunay triangulation is defined
!    as a set of triangles with the following five properties:
!
!    1)  The triangle vertices are nodes.
!    2)  No triangle contains a node other than its vertices.
!    3)  The interiors of the triangles are pairwise disjoint.
!    4)  The union of triangles is the convex hull of the set
!        of nodes (the smallest convex set which contains
!        the nodes).
!    5)  The interior of the circumcircle of each triangle
!        contains no node.
!
!    The first four properties define a triangulation, and the
!    last property results in a triangulation which is as close
!    as possible to equiangular in a certain sense and which is
!    uniquely defined unless four or more nodes lie on a common
!    circle.  This property makes the triangulation well-suited
!    for solving closest point problems and for triangle-based
!    interpolation.
!
!    The triangulation can be generalized to a constrained
!    Delaunay triangulation by a call to Subroutine ADDCST.
!    This allows for user-specified boundaries defining a non-
!    convex and/or multiply connected region.
!
!    The algorithm for constructing the triangulation has
!    expected time complexity O(N*log(N)) for most nodal dis-
!    tributions.  Also, since the algorithm proceeds by adding
!    nodes incrementally, the triangulation may be updated with
!    the addition (or deletion) of a node very efficiently.
!    The adjacency information representing the triangulation
!    is stored as a linked list requiring approximately 13N
!    storage locations.
!
!
!    The following is a list of the software package modules
!    which a user may wish to call directly:
!
!    ADDCST - Generalizes the Delaunay triangulation to allow
!             for user-specified constraints.
!
!    ADDNOD - Updates the triangulation by appending or
!             inserting a new node.
!
!    AREAP  - Computes the area bounded by a closed polygonal
!             curve such as the boundary of the triangula-
!             tion or of a constraint region.
!
!    BNODES - Returns an array containing the indexes of the
!             boundary nodes in counterclockwise order.
!             Counts of boundary nodes, triangles, and arcs
!             are also returned.
!
!    CIRCUM - Computes the area, circumcenter, circumradius,
!             and, optionally, the aspect ratio of a trian-
!             gle defined by user-specified vertices.
!
!    DELARC - Deletes a boundary arc from the triangulation.
!
!    DELNOD - Updates the triangulation with the deletion of a
!             node.
!
!    EDGE   - Forces a pair of nodes to be connected by an arc
!             in the triangulation.
!
!    GETNP  - Determines the ordered sequence of L closest
!             nodes to a given node, along with the associ-
!             ated distances.  The distance between nodes is
!             taken to be the length of the shortest connec-
!             ting path which intersects no constraint
!             region.
!
!    INTSEC - Determines whether or not an arbitrary pair of
!             line segments share a common point.
!
!    JRAND  - Generates a uniformly distributed pseudo-random
!             integer ( kind = 4 ).
!
!    LEFT   - Locates a point relative to a line.
!
!    NEARND - Returns the index of the nearest node to an
!             arbitrary point, along with its squared
!             distance.
!
!    STORE  - Forces a value to be stored in main memory so
!             that the precision of floating point numbers
!             in memory locations rather than registers is
!             computed.
!
!    TRLIST - Converts the triangulation data structure to a
!             triangle list more suitable for use in a fin-
!             ite element code.
!
!    TRLPRT - Prints the triangle list created by TRLIST.
!
!    TRMESH - Creates a Delaunay triangulation of a set of nodes.
!
!    TRMSHR - Creates a Delaunay triangulation (more effici-
!             ently than TRMESH) of a set of nodes lying at
!             the vertices of a (possibly skewed) rectangu-
!             lar grid.
!
!    TRPLOT - Creates a level-2 Encapsulated Postscript (EPS)
!             file containing a triangulation plot.
!
!    TRPRNT - Prints the triangulation data structure and,
!             optionally, the nodal coordinates.
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
!    N >= 3.
!
!    Input, real ( kind = 8 ) X(N), Y(N), the coordinates of the nodes. 
!    (X(K),Y(K)) is referred to as node K, and K is referred to as a nodal
!    index.  The first three nodes must not be collinear.
!
!    Output, integer ( kind = 4 ) LIST(6*N-12), nodal indexes which, along with LPTR,
!    LEND, and LNEW, define the triangulation as a set of N adjacency 
!    lists; counterclockwise-ordered sequences of neighboring nodes such
!    that the first and last neighbors of a boundary node are boundary 
!    nodes (the first neighbor of an interior node is arbitrary).  In
!    order to distinguish between interior and boundary nodes, the last 
!    neighbor of each boundary node is represented by the negative
!    of its index.
!
!    Output, integer ( kind = 4 ) LPTR(6*N-12), pointers (LIST indexes) in one-to-one
!    correspondence with the elements of LIST.  LIST(LPTR(I)) indexes the 
!    node which follows LIST(I) in cyclical counterclockwise order
!    (the first neighbor follows the last neighbor).
!
!    Output, integer ( kind = 4 ) LEND(N), pointers to adjacency lists.  LEND(K)
!    points to the last neighbor of node K for K = 1,...,N.  Thus, 
!    LIST(LEND(K)) < 0 if and only if K is a boundary node.
!
!    Output, integer ( kind = 4 ) LNEW, pointer to the first empty location in LIST
!    and LPTR (list length plus one).  LIST, LPTR, LEND, and LNEW are 
!    not altered if IER < 0, and are incomplete if IER > 0.
!
!    Workspace NEAR(N), NEXT(N), DIST(N).  The space is used to efficiently
!    determine the nearest triangulation node to each unprocessed node for 
!    use by ADDNOD.
!
!    Output, integer ( kind = 4 ) IER = Error indicator:
!     0 if no errors were encountered.
!    -1 if N < 3 on input.
!    -2 if the first three nodes are collinear.
!    -4 if an error flag was returned by a call to SWAP in ADDNOD.  This is 
!      an internal error and should be reported to the programmer.
!     L if nodes L and M coincide for some M > L.  The linked list represents
!      a triangulation of nodes 1 to M-1 in this case.
!
!  Local parameters:
!
!    D =        Squared distance from node K to node I
!    D1,D2,D3 = Squared distances from node K to nodes 1, 2,
!               and 3, respectively
!    EPS =      Half the machine precision
!    I,J =      Nodal indexes
!    I0 =       Index of the node preceding I in a sequence of
!               unprocessed nodes:  I = NEXT(I0)
!    K =        Index of node to be added and DO-loop index: K > 3
!    KM1 =      K-1
!    LCC(1) =   Dummy array
!    LP =       LIST index (pointer) of a neighbor of K
!    LPL =      Pointer to the last neighbor of K
!    NCC =      Number of constraint curves
!    NEXTI =    NEXT(I)
!    NN =       Local copy of N
!    SWTOL =    Tolerance for function SWPTST
!
  implicit none
  
  integer ( kind = 4 ) n
  
  real ( kind = 8 ) d
  real ( kind = 8 ) d1
  real ( kind = 8 ) d2
  real ( kind = 8 ) d3
  real ( kind = 8 ) dist(n)
  real ( kind = 8 ) eps
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i0
  integer ( kind = 4 ) ier
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) km1
  integer ( kind = 4 ) lcc(1)
  logical left
  integer ( kind = 4 ) lend(n)
  integer ( kind = 4 ) list(*)
  integer ( kind = 4 ) lnew
  integer ( kind = 4 ) lp
  integer ( kind = 4 ) lpl
  integer ( kind = 4 ) lptr(*)
  integer ( kind = 4 ) ncc
  integer ( kind = 4 ) near(n)
  integer ( kind = 4 ) next(n)
  integer ( kind = 4 ) nexti
  integer ( kind = 4 ) nn
  real ( kind = 8 ) swtol
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)
  
  common /swpcom/ swtol
  
  nn = n
  
  if ( nn < 3 ) then
     ier = -1
     return
  end if
!
!  Compute a tolerance for function SWPTST:  SWTOL = 10*
!  (machine precision)
!
  eps = epsilon ( eps )

  swtol = eps * 20.0D+00
!
!  Store the first triangle in the linked list.
!
  if ( .not. left ( x(1), y(1), x(2), y(2), x(3), y(3) ) ) then
!
!  The initial triangle is (3,2,1) = (2,1,3) = (1,3,2).
!
     list(1) = 3
     lptr(1) = 2
     list(2) = -2
     lptr(2) = 1
     lend(1) = 2
     
     list(3) = 1
     lptr(3) = 4
     list(4) = -3
     lptr(4) = 3
     lend(2) = 4
     
     list(5) = 2
     lptr(5) = 6
     list(6) = -1
     lptr(6) = 5
     lend(3) = 6
     
  else if ( .not. left(x(2),y(2),x(1),y(1),x(3),y(3)) ) then
!
!  The initial triangle is (1,2,3).
!
     list(1) = 2
     lptr(1) = 2
     list(2) = -3
     lptr(2) = 1
     lend(1) = 2
     
     list(3) = 3
     lptr(3) = 4
     list(4) = -1
     lptr(4) = 3
     lend(2) = 4
     
     list(5) = 1
     lptr(5) = 6
     list(6) = -2
     lptr(6) = 5
     lend(3) = 6

  else
!
!  The first three nodes are collinear.
!
     ier = -2
     return
  end if
!
!  Initialize LNEW and test for N = 3.
!
  lnew = 7
  if (nn == 3) then
     ier = 0
     return
  end if
!
!  A nearest-node data structure (NEAR, NEXT, and DIST) is
!  used to obtain an expected-time (N*log(N)) incremental
!  algorithm by enabling constant search time for locating
!  each new node in the triangulation.
!
!  For each unprocessed node K, NEAR(K) is the index of the
!  triangulation node closest to K (used as the starting
!  point for the search in Subroutine TRFIND) and DIST(K)
!  is an increasing function of the distance between nodes
!  K and NEAR(K).
!
!  Since it is necessary to efficiently find the subset of
!  unprocessed nodes associated with each triangulation
!  node J (those that have J as their NEAR entries), the
!  subsets are stored in NEAR and NEXT as follows:  for
!  each node J in the triangulation, I = NEAR(J) is the
!  first unprocessed node in J's set (with I = 0 if the
!  set is empty), L = NEXT(I) (if I > 0) is the second,
!  NEXT(L) (if L > 0) is the third, etc.  The nodes in each
!  set are initially ordered by increasing indexes (which
!  maximizes efficiency) but that ordering is not main-
!  tained as the data structure is updated.
!
!  Initialize the data structure for the single triangle.
!
  near(1) = 0
  near(2) = 0
  near(3) = 0
  
  do k = nn, 4, -1

     d1 = ( x(k) - x(1) )**2 + ( y(k) - y(1) )**2
     d2 = ( x(k) - x(2) )**2 + ( y(k) - y(2) )**2
     d3 = ( x(k) - x(3) )**2 + ( y(k) - y(3) )**2
     
     if ( d1 <= d2  .and.  d1 <= d3 ) then
        near(k) = 1
        dist(k) = d1
        next(k) = near(1)
        near(1) = k
     else if (d2 <= d1  .and.  d2 <= d3) then
        near(k) = 2
        dist(k) = d2
        next(k) = near(2)
        near(2) = k
     else
        near(k) = 3
        dist(k) = d3
        next(k) = near(3)
        near(3) = k
     end if

  end do
!
!  Add the remaining nodes.  Parameters for ADDNOD are as follows:
!
!  K = Index of the node to be added.
!  NEAR(K) = Index of the starting node for the search in TRFIND.
!  NCC = Number of constraint curves.
!  LCC = Dummy array (since NCC = 0).
!  KM1 = Number of nodes in the triangulation.
!
  ncc = 0

  do k = 4, nn

     km1 = k-1
     
     call addnod ( k, x(k), y(k), near(k), ncc, lcc, km1, x, y, &
          list, lptr, lend, lnew, ier )
     
     if ( ier /= 0 ) then
        return
     end if
!
!  Remove K from the set of unprocessed nodes associated with NEAR(K).
!
     i = near(k)
     
     if (near(i) == k) then
        
        near(i) = next(k)
        
     else
        
        i = near(i)
        
        do
           
           i0 = i
           i = next(i0)
           if (i == k) then
              exit
           end if
           
        end do
        
        next(i0) = next(k)

     end if

     near(k) = 0
!
!  Loop on neighbors J of node K.
!
     lpl = lend(k)
     lp = lpl

4    continue
     
     lp = lptr(lp)
     j = abs ( list(lp) )
!
!  Loop on elements I in the sequence of unprocessed nodes
!  associated with J:  K is a candidate for replacing J
!  as the nearest triangulation node to I.  The next value
!  of I in the sequence, NEXT(I), must be saved before I
!  is moved because it is altered by adding I to K's set.
!
     i = near(j)

5    continue

     if ( i == 0 ) go to 6
     nexti = next(i)
!
!  Test for the distance from I to K less than the distance
!  from I to J.
!
     d = (x(k)-x(i))**2 + (y(k)-y(i))**2
!
!  Replace J by K as the nearest triangulation node to I:
!  update NEAR(I) and DIST(I), and remove I from J's set
!  of unprocessed nodes and add it to K's set.
!
     if ( d < dist(i) ) then
        near(i) = k
        dist(i) = d
        if (i == near(j)) then
           near(j) = nexti
        else
           next(i0) = nexti
        end if
        next(i) = near(k)
        near(k) = i
     else
        i0 = i
     end if
!
!  Bottom of loop on I.
!
     i = nexti
     go to 5
!
!  Bottom of loop on neighbors J.
!
6    continue

     if ( lp /= lpl ) then
        go to 4
     end if
     
  end do
  
  return

end subroutine trmesh




subroutine trmshr ( n, nx, x, y, nit, list, lptr, lend, lnew, ier )

  !*****************************************************************************80
  !
  !! TRMSHR triangulates logically rectangular data.
  !
  !  Discussion:
  !
  !    This subroutine creates a Delaunay triangulation of a
  !    set of N nodes in the plane, where the nodes are the vert-
  !    ices of an NX by NY skewed rectangular grid with the
  !    natural ordering.  Thus, N = NX*NY, and the nodes are
  !    ordered from left to right beginning at the top row so
  !    that adjacent nodes have indexes which differ by 1 in the
  !    x-direction and by NX in the y-direction.  A skewed rec-
  !    tangular grid is defined as one in which each grid cell is
  !    a strictly convex quadrilateral (and is thus the convex
  !    hull of its four vertices).  Equivalently, any transfor-
  !    mation from a rectangle to a grid cell which is bilinear
  !    in both components has an invertible Jacobian.
  !
  !    If the nodes are not distributed and ordered as defined
  !    above, Subroutine TRMESH must be called in place of this
  !    routine.  Refer to Subroutine ADDCST for the treatment of
  !    constraints.
  !
  !    The first phase of the algorithm consists of construc-
  !    ting a triangulation by choosing a diagonal arc in each
  !    grid cell.  If NIT = 0, all diagonals connect lower left
  !    to upper right corners and no error checking or additional
  !    computation is performed.  Otherwise, each diagonal arc is
  !    chosen to be locally optimal, and boundary arcs are added
  !    where necessary in order to cover the convex hull of the
  !    nodes.  (This is the first iteration.)  If NIT > 1 and no
  !    error was detected, the triangulation is then optimized by
  !    a sequence of up to NIT-1 iterations in which interior
  !    arcs of the triangulation are tested and swapped if appro-
  !    priate.  The algorithm terminates when an iteration
  !    results in no swaps and/or when the allowable number of
  !    iterations has been performed.  NIT = 0 is sufficient to
  !    produce a Delaunay triangulation if the original grid is
  !    actually rectangular, and NIT = 1 is sufficient if it is
  !    close to rectangular.  Note, however, that the ordering
  !    and distribution of nodes is not checked for validity in
  !    the case NIT = 0, and the triangulation will not be valid
  !    unless the rectangular grid covers the convex hull of the
  !    nodes.
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
  !    Input, integer ( kind = 4 ) N, the number of nodes in the grid.  N = NX*NY for some
  !    NY >= 2.
  !
  !    Input, integer ( kind = 4 ) NX, the number of grid points in the x-direction.
  !    NX >= 2.
  !
  !    Input, real ( kind = 8 ) X(N), Y(N), the coordinates of the nodes with 
  !    the ordering and distribution defined in the header comments above.
  !    (X(K),Y(K)) is referred to as node K.
  !
  !    Input/output, integer ( kind = 4 ) NIT.  On input, the maximum number of iterations
  !    to be employed.  Refer to the header comments above.  On output,
  !    the actual number of iterations.
  !
  !    Output, integer ( kind = 4 ) LIST(6*N-12), LPTR(6*N-12), LEND(N), LNEW, 
  !    data structure defining the triangulation.  Refer to subroutine TRMESH.
  !
  !    Output, integer ( kind = 4 ) IER = Error indicator:
  !     0 if no errors were encountered.
  !     K if the grid element with upper left corner at node K is not a strictly
  !      convex quadrilateral.  The algorithm is terminated when the first such
  !      occurrence is detected.  Note that this test is not performed if
  !      NIT = 0 on input.
  !    -1 if N, NX, or NIT is outside its valid range on input.
  !    -2 if NIT > 1 on input, and the optimization loop failed to converge
  !      within the allowable number of iterations.  The triangulation is
  !      valid but not optimal in this case.
  !
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) eps
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ier
  integer ( kind = 4 ) iter
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kp1
  logical left
  integer ( kind = 4 ) lend(n)
  integer ( kind = 4 ) list(*)
  integer ( kind = 4 ) lnew
  integer ( kind = 4 ) lp
  integer ( kind = 4 ) lpf
  integer ( kind = 4 ) lpk
  integer ( kind = 4 ) lpl
  integer ( kind = 4 ) lpp
  integer ( kind = 4 ) lptr(*)
  integer ( kind = 4 ) lstptr
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) m3
  integer ( kind = 4 ) m4
  integer ( kind = 4 ) maxit
  integer ( kind = 4 ) n0
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) n3
  integer ( kind = 4 ) n4
  integer ( kind = 4 ) nbcnt
  integer ( kind = 4 ) ni
  integer ( kind = 4 ) nit
  integer ( kind = 4 ) nj
  integer ( kind = 4 ) nm1
  integer ( kind = 4 ) nn
  integer ( kind = 4 ) nnb
  integer ( kind = 4 ) nx
  logical swptst
  real ( kind = 8 ) swtol
  logical tst
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)

  common /swpcom/ swtol
  !
  !  Store local variables and test for errors in input parameters.
  !
  ni = nx
  nj = n / ni
  nn = ni * nj
  maxit = nit
  nit = 0

  if ( n /= nn .or. nj < 2 .or. ni < 2 .or. maxit < 0 ) then
     ier = -1
     return
  end if

  ier = 0
  !
  !  Compute a tolerance for function SWPTST.
  !
  eps = epsilon ( eps )
  swtol = eps * 20.0D+00
  !
  !  Loop on grid points (I,J) corresponding to nodes K =
  !  (J-1)*NI + I.  TST = TRUE iff diagonals are to be
  !  chosen by the swap test.  M1, M2, M3, and M4 are the
  !  slopes (-1, 0, or 1) of the diagonals in quadrants 1
  !  to 4 (counterclockwise beginning with the upper right)
  !  for a coordinate system with origin at node K.
  !
  tst = maxit > 0
  m1 = 0
  m4 = 0
  lp = 0
  kp1 = 1

  do 6 j = 1,nj

     do 5 i = 1,ni

        m2 = m1
        m3 = m4
        k = kp1
        kp1 = k + 1
        lpf = lp + 1

        if ( j == nj .and. i /= ni ) go to 2

        if ( i /= 1 ) then

           if ( j /= 1 ) then
              !
              !  K is not in the top row, leftmost column, or bottom row
              !  (unless K is the lower right corner).  Take the first
              !  neighbor to be the node above K.
              !
              lp = lp + 1
              list(lp) = k - ni
              lptr(lp) = lp + 1

              if ( m2 <= 0 ) then
                 lp = lp + 1
                 list(lp) = k - 1 - ni
                 lptr(lp) = lp + 1
              end if

           end if
           !
           !  K is not in the leftmost column.  The next (or first)
           !  neighbor is to the left of K.
           !
           lp = lp + 1
           list(lp) = k - 1
           lptr(lp) = lp + 1
           if (j == nj) go to 3

           if (m3 >= 0) then
              lp = lp + 1
              list(lp) = k - 1 + ni
              lptr(lp) = lp + 1
           end if

        end if
        !
        !  K is not in the bottom row.  The next (or first) neighbor is below K.
        !
        lp = lp + 1
        list(lp) = k + ni
        lptr(lp) = lp + 1
        !
        !  Test for a negative diagonal in quadrant 4 unless K is
        !  in the rightmost column.  The quadrilateral associated
        !  with the quadrant is tested for strict convexity un-
        !  less NIT = 0 on input.
        !
        if ( i == ni ) then
           go to 3
        end if

        m4 = 1
        if ( .not. tst ) go to 2

        if ( left(x(kp1),y(kp1),x(k+ni),y(k+ni),x(k),y(k)) .or. &
             left(x(k),y(k),x(kp1+ni),y(kp1+ni),x(k+ni),y(k+ni)) .or. &
             left(x(k+ni),y(k+ni),x(kp1),y(kp1), x(kp1+ni),y(kp1+ni)) .or. &
             left(x(kp1+ni),y(kp1+ni),x(k),y(k), x(kp1),y(kp1)) ) then
           ier = k
           return
        end if

        if ( swptst ( kp1, k+ni, k, kp1+ni, x, y ) ) go to 2

        m4 = -1
        lp = lp + 1
        list(lp) = kp1 + ni
        lptr(lp) = lp + 1
        !
        !  The next (or first) neighbor is to the right of K.
        !
2       continue

        lp = lp + 1
        list(lp) = kp1
        lptr(lp) = lp + 1
        !
        !  Test for a positive diagonal in quadrant 1 (the neighbor
        !  of K-NI which follows K is not K+1) unless K is in the
        !  top row.
        !
        if (j == 1) go to 3

        if (tst) then
           m1 = -1
           lpk = lstptr(lend(k-ni),k,list,lptr)
           lpk = lptr(lpk)

           if ( list(lpk) /= kp1 ) then
              m1 = 1
              lp = lp + 1
              list(lp) = kp1 - ni
              lptr(lp) = lp + 1
           end if

        end if
        !
        !  If K is in the leftmost column (and not the top row) or
        !  in the bottom row (and not the rightmost column), then
        !  the next neighbor is the node above K.
        !
        if ( i /= 1 .and. j /= nj ) go to 4

        lp = lp + 1
        list(lp) = k - ni
        lptr(lp) = lp + 1
        if ( i == 1 ) go to 3
        !
        !  K is on the bottom row (and not the leftmost or rightmost column).
        !
        if ( m2 <= 0 ) then
           lp = lp + 1
           list(lp) = k - 1 - ni
           lptr(lp) = lp + 1
        end if

        lp = lp + 1
        list(lp) = k - 1
        lptr(lp) = lp + 1
        !
        !  K is a boundary node.
        !
3       continue

        list(lp) = -list(lp)
        !
        !  Bottom of loop.  Store LEND and correct LPTR(LP).
        !  LPF and LP point to the first and last neighbors of K.
        !
4       continue

        lend(k) = lp
        lptr(lp) = lpf

5       continue
     !end do !Added 7/4/2013 [JH]
6    continue
  !end do !Added 7/4/2013 [JH]
  !
  !  Store LNEW, and terminate the algorithm if NIT = 0 on input.
  !
  lnew = lp + 1

  if ( maxit == 0 ) then
     return
  end if
  !
  !  Add boundary arcs where necessary in order to cover the
  !  convex hull of the nodes.  N1, N2, and N3 are consecu-
  !  tive boundary nodes in counterclockwise order, and N0
  !  is the starting point for each loop around the boundary.
  !
  n0 = 1
  n1 = n0
  n2 = ni + 1
  !
  !  TST is set to TRUE if an arc is added.  The boundary
  !  loop is repeated until a traversal results in no
  !  added arcs.
  !
7 continue

  tst = .false.
  !
  !  Top of boundary loop.  Set N3 to the first neighbor of
  !  N2, and test for N3 LEFT N1 -> N2.
  !
8 continue

  lpl = lend(n2)

  lp = lptr(lpl)
  n3 = list(lp)

  if ( left(x(n1),y(n1),x(n2),y(n2),x(n3),y(n3)) ) then
     n1 = n2
  end if

  if (n1 /= n2) then
     !
     !  Add the boundary arc N1-N3.  If N0 = N2, the starting
     !  point is changed to N3, since N2 will be removed from
     !  the boundary.  N3 is inserted as the first neighbor of
     !  N1, N2 is changed to an interior node, and N1 is
     !  inserted as the last neighbor of N3.
     !
     tst = .true.
     if (n2 == n0) n0 = n3
     lp = lend(n1)
     call insert (n3,lp, list,lptr,lnew )
     list(lpl) = -list(lpl)
     lp = lend(n3)
     list(lp) = n2
     call insert (-n1,lp, list,lptr,lnew )
     lend(n3) = lnew - 1
  end if
  !
  !  Bottom of loops.  Test for termination.
  !
  n2 = n3
  if (n1 /= n0) go to 8
  if (tst) go to 7
  !
  !  Terminate the algorithm if NIT = 1 on input.
  !
  nit = 1

  if ( maxit == 1 ) then
     return
  end if
  !
  !  Optimize the triangulation by applying the swap test and
  !  appropriate swaps to the interior arcs.  The loop is
  !  repeated until no swaps are performed or MAXIT itera-
  !  tions have been applied.  ITER is the current iteration,
  !  and TST is set to TRUE if a swap occurs.
  !
  iter = 1
  nm1 = nn - 1

9 continue

  iter = iter + 1
  tst = .false.
  !
  !  Loop on interior arcs N1-N2, where N2 > N1 and
  !  (N1,N2,N3) and (N2,N1,N4) are adjacent triangles.
  !
  !  Top of loop on nodes N1.
  !
  do 11 n1 = 1,nm1

     lpl = lend(n1)
     n4 = list(lpl)
     lpf = lptr(lpl)
     n2 = list(lpf)
     lp = lptr(lpf)
     n3 = list(lp)
     nnb = nbcnt(lpl,lptr)
     !
     !  Top of loop on neighbors N2 of N1.  NNB is the number of
     !  neighbors of N1.
     !
     do i = 1,nnb
        !
        !  Bypass the swap test if N1 is a boundary node and N2 is
        !  the first neighbor (N4 < 0), N2 < N1, or N1-N2 is a
        !  diagonal arc (already locally optimal) when ITER = 2.
        !
        if ( n4 > 0  .and.  n2 > n1  .and. &
             ( iter /= 2  .or.  abs ( n1+ni-n2 ) /= 1 ) ) then

           if (swptst(n3,n4,n1,n2,x,y) ) then
              !
              !  Swap diagonal N1-N2 for N3-N4, set TST to TRUE, and set
              !  N2 to N4 (the neighbor preceding N3).
              !
              call swap (n3,n4,n1,n2, list,lptr,lend, lpp)
              if (lpp /= 0) then
                 tst = .true.
                 n2 = n4
              end if
           end if
        end if
        !
        !  Bottom of neighbor loop.
        !
        if (list(lpl) == -n3) then
           go to 11
        end if

        n4 = n2
        n2 = n3
        lp = lstptr(lpl,n2,list,lptr)
        lp = lptr(lp)
        n3 = abs ( list(lp) )

     end do

11   continue
  !end do !Added 7/4/2013 [JH]
  !
  !  Test for termination.
  !
  if (tst  .and.  iter < maxit) go to 9

  nit = iter

  if ( tst ) then
     ier = -2
  end if

  return

end subroutine trmshr




function store ( x )

  !*****************************************************************************80
  !
  !! STORE forces its argument to be stored.
  !
  !  Discussion:
  !
  !    This function forces its argument X to be stored in a
  !    memory location, thus providing a means of determining
  !    floating point number characteristics (such as the machine
  !    precision) when it is necessary to avoid computation in
  !    high precision registers.
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
  !    Input, real ( kind = 8 ) X, the value to be stored.
  !
  !    Output, real ( kind = 8 ) STORE, the value of X after it has been stored
  !    and possibly truncated or rounded to the single precision word length.
  !
 implicit none

 real ( kind = 8 ) store
 real ( kind = 8 ) x
 real ( kind = 8 ) y

 common /stcom/ y

 y = x
 store = y

 return

end function store
