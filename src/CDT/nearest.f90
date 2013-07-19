subroutine getnp ( ncc, lcc, n, x, y, list, lptr, lend, l, npts, ds, ier )

!*****************************************************************************80
!
!! GETNP sets the next nearest node to a given node.
!
!  Discussion:
!
!    Given a triangulation of N nodes and an array NPTS con-
!    taining the indexes of L-1 nodes ordered by distance from
!    NPTS(1), this subroutine sets NPTS(L) to the index of the
!    next node in the sequence -- the node, other than NPTS(1),
!    ...,NPTS(L-1), which is closest to NPTS(1).  Thus, the
!    ordered sequence of K closest nodes to N1 (including N1)
!    may be determined by K-1 calls to GETNP with NPTS(1) = N1
!    and L = 2,3,...,K for K >= 2.  Note that NPTS must 
!    include constraint nodes as well as non-constraint nodes.
!    Thus, a sequence of K1 closest non-constraint nodes to N1
!    must be obtained as a subset of the closest K2 nodes to N1
!    for some K2 >= K1.
!
!    The terms closest and distance have special definitions
!    when constraint nodes are present in the triangulation.
!    Nodes N1 and N2 are said to be visible from each other if
!    and only if the line segment N1-N2 intersects no constraint
!    arc (except possibly itself) and is not an interi-
!    or constraint arc (arc whose interior lies in a constraint
!    region).  A path from N1 to N2 is an ordered sequence of
!    nodes, with N1 first and N2 last, such that adjacent path
!    elements are visible from each other.  The path length is
!    the sum of the Euclidean distances between adjacent path
!    nodes.  Finally, the distance from N1 to N2 is defined to
!    be the length of the shortest path from N1 to N2.
!
!    The algorithm uses the property of a Delaunay triangulation
!    that the K-th closest node to N1 is a neighbor of one
!    of the K-1 closest nodes to N1.  With the definition of
!    distance used here, this property holds when constraints
!    are present as long as non-constraint arcs are locally
!    optimal.
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
!    Input, integer ( kind = 4 ) LCC(*), a list of constraint curve starting 
!    indexes (or dummy array of length 1 if NCC = 0).  Refer to subroutine 
!    ADDCST.
!
!    Input, integer ( kind = 4 ) N, the number of nodes in the triangulation.  
!    N >= 3.
!
!    Input, real ( kind = 8 ) X(N), Y(N), the coordinates of the nodes with 
!    non-constraint nodes in the first LCC(1)-1 locations if NCC > 0.
!
!    Input, integer ( kind = 4 ) LIST(*), LPTR(*), LEND(N), the triangulation 
!    data structure.  Refer to subroutine TRMESH.
!
!    Input, integer ( kind = 4 ) L, the number of nodes in the sequence 
!    on output.  2 <= L <= N.
!
!    Input/output, integer ( kind = 4 ) NPTS(L), on input, the indexes of the 
!    L-1 closest nodes to NPTS(1) in the first L-1 locations.  On output, 
!    updated with the index of the L-th closest node to NPTS(1) in position L 
!    unless IER /= 0.
!
!    Input/output, real ( kind = 8 ) DS(L), the distance (defined above) 
!    between NPTS(1) and NPTS(I) in the I-th position for I = 1,...,L-1.  
!    Thus, DS(1) = 0.  On output, updated with the distance between NPTS(1) 
!    and NPTS(L) in position L unless IER /= 0.
!
!    Output, integer ( kind = 4 ) IER = Error indicator:
!     0 if no errors were encountered.
!    -1 if NCC, N, L, or an LCC entry is outside its valid range on input.
!     K if NPTS(K) is not a valid index in the range 1 to N.
!
  implicit none

  integer ( kind = 4 ) l
  integer ( kind = 4 ) n

  real ( kind = 8 ) dc
  real ( kind = 8 ) dl
  real ( kind = 8 ) ds(l)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ier
  integer ( kind = 4 ) ifrst
  integer ( kind = 4 ) ilast
  logical intsec
  logical isw
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) km1
  integer ( kind = 4 ) lcc(*)
  integer ( kind = 4 ) lcc1
  integer ( kind = 4 ) lend(n)
  logical lft1
  logical lft2
  logical lft12
  integer ( kind = 4 ) list(*)
  integer ( kind = 4 ) lm1
  integer ( kind = 4 ) lp
  integer ( kind = 4 ) lpcl
  integer ( kind = 4 ) lpk
  integer ( kind = 4 ) lpkl
  integer ( kind = 4 ) lptr(*)
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) nc
  integer ( kind = 4 ) ncc
  logical ncf
  integer ( kind = 4 ) nf1
  integer ( kind = 4 ) nf2
  integer ( kind = 4 ) nj
  logical njf
  integer ( kind = 4 ) nk
  integer ( kind = 4 ) nkbak
  integer ( kind = 4 ) nkfor
  integer ( kind = 4 ) nl
  integer ( kind = 4 ) nn
  integer ( kind = 4 ) npts(l)
  logical skip
  logical sksav
  logical vis
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) x1
  real ( kind = 8 ) xc
  real ( kind = 8 ) xj
  real ( kind = 8 ) xk
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) y1
  real ( kind = 8 ) yc
  real ( kind = 8 ) yj
  real ( kind = 8 ) yk
!
!  Store parameters in local variables and test for errors.
!  LCC1 indexes the first constraint node.
!
  ier = -1
  nn = n
  lcc1 = nn + 1
  lm1 = l - 1

  if ( ncc < 0 ) then
     return
  end if
  
  if ( lm1 < 1  .or.  lm1 >= nn) then
     return
  end if
  
  if (ncc == 0) then
     if (nn < 3) then
        return
     end if
  else
     do i = ncc,1,-1
        if (lcc1 - lcc(i) < 3) return
        lcc1 = lcc(i)
     end do
     
     if (lcc1 < 1) then
        return
     end if
     
  end if
!
!  Test for an invalid index in NPTS.
!
  do k = 1,lm1
     nk = npts(k)
     if ( nk < 1  .or.  nk > nn ) then
        ier = k
        return
     end if
  end do
!
!  Store N1 = NPTS(1) and mark the elements of NPTS.
!
  n1 = npts(1)
  x1 = x(n1)
  y1 = y(n1)
  
  do k = 1,lm1
     nk = npts(k)
     lend(nk) = -lend(nk)
  end do
!
!  Candidates NC for NL = NPTS(L) are the unmarked visible
!  neighbors of nodes NK in NPTS.  ISW is an initialization
!  switch set to .TRUE. when NL and its distance DL from N1
!  have been initialized with the first candidate encountered.
!
  isw = .false.
  dl = 0.0D+00
!
!  Loop on marked nodes NK = NPTS(K).  LPKL indexes the last
!  neighbor of NK in LIST.
!
  do k = 1,lm1
     
     km1 = k - 1
     nk = npts(k)
     xk = x(nk)
     yk = y(nk)
     lpkl = -lend(nk)
     nkfor = 0
     nkbak = 0
     vis = .true.
!
!  NK is a constraint node.  Set NKFOR and NKBAK to the
!  constraint nodes which follow and precede NK.  IFRST
!  and ILAST are set to the first and last nodes in the
!  constraint containing NK.
!
     if (nk >= lcc1) then
        
        ifrst = nn + 1
        
        do i = ncc,1,-1
           ilast = ifrst - 1
           ifrst = lcc(i)
           if ( nk >= ifrst ) then
              exit
           end if
        end do
        
        if (nk < ilast) then
           nkfor = nk + 1
        else
           nkfor = ifrst
        end if

        if (nk > ifrst) then
           nkbak = nk - 1
        else
           nkbak = ilast
        end if
!
!  Initialize VIS to TRUE iff NKFOR precedes NKBAK in the
!  adjacency list for NK -- the first neighbor is visible and is not NKBAK.
!
        lpk = lpkl

        do

           lpk = lptr(lpk)
           nc = abs ( list(lpk) )
           
           if ( nc == nkfor .or. nc == nkbak ) then
              exit
           end if
           
        end do
        
        vis = nc == nkfor
        
     end if
!
!  Loop on neighbors NC of NK, bypassing marked and nonvisible neighbors.
!
     lpk = lpkl

7    continue

     lpk = lptr(lpk)
     nc = abs ( list(lpk) )
     
     if ( nc == nkbak ) then
        vis = .true.
     end if
!
!  VIS = .FALSE. iff NK-NC is an interior constraint arc
!  (NK is a constraint node and NC lies strictly between
!  NKFOR and NKBAK).
!
     if ( .not. vis ) go to 15

     if ( nc == nkfor ) then
        vis = .false.
     end if

     if ( lend(nc) < 0 ) go to 15
!
!  Initialize distance DC between N1 and NC to Euclidean distance.
!
     xc = x(nc)
     yc = y(nc)
     dc = sqrt((xc-x1)*(xc-x1) + (yc-y1)*(yc-y1))
     if (isw  .and.  dc >= dl) go to 15

     if (k == 1) then
        go to 14
     end if
!
!  K >= 2.  Store the pointer LPCL to the last neighbor of NC.
!
     lpcl = lend(nc)
!
!  Set DC to the length of the shortest path from N1 to NC
!  which has not previously been encountered and which is
!  a viable candidate for the shortest path from N1 to NL.
!  This is Euclidean distance iff NC is visible from N1.
!  Since the shortest path from N1 to NL contains only ele-
!  ments of NPTS which are constraint nodes (in addition to
!  N1 and NL), only these need be considered for the path
!  from N1 to NC.  Thus, for distance function D(A,B) and
!  J = 1,...,K, DC = min(D(N1,NJ) + D(NJ,NC)) over con-
!  straint nodes NJ = NPTS(J) which are visible from NC.
!
     do j = 1,km1
        
        nj = npts(j)

        if ( 1 < j .and.  nj < lcc1 ) then
           go to 13
        end if
!
!  If NC is a visible neighbor of NJ, a path from N1 to NC
!  containing NJ has already been considered.  Thus, NJ may
!  be bypassed if it is adjacent to NC.
!
        lp = lpcl

8       continue

        lp = lptr(lp)

        if ( nj == abs ( list(lp) ) ) then
           go to 12
        end if
        
        if (lp /= lpcl) then
           go to 8
        end if
!
!  NJ is a constraint node (unless J=1) not adjacent to NC,
!  and is visible from NC iff NJ-NC is not intersected by
!  a constraint arc.  Loop on constraints I in reverse
!  order.
!
        xj = x(nj)
        yj = y(nj)
        ifrst = nn+1
        
        do 11 i = ncc,1,-1
           ilast = ifrst - 1
           ifrst = lcc(i)
           nf1 = ilast
           ncf = nf1 == nc
           njf = nf1 == nj
           skip = ncf  .or.  njf
!
!  Loop on boundary constraint arcs NF1-NF2 which contain
!  neither NC nor NJ.  NCF and NJF are TRUE iff NC (or NJ)
!  has been encountered in the constraint, and SKIP =
!  .TRUE. iff NF1 = NC or NF1 = NJ.
!
           do nf2 = ifrst,ilast
              
              if (nf2 == nc) ncf = .true.
              if (nf2 == nj) njf = .true.
              sksav = skip
              skip = nf2 == nc  .or.  nf2 == nj
!
!  The last constraint arc in the constraint need not be
!  tested if none of the arcs have been skipped.
!
              if ( sksav  .or.  skip  .or. &
                   (nf2 == ilast  .and. &
                   .not. ncf  .and.  .not. njf) ) then
                 go to 9
              end if
              
              if ( intsec(x(nf1),y(nf1),x(nf2),y(nf2), &
                   xc,yc,xj,yj) ) then
                 go to 12
              end if
              
9             continue

              nf1 = nf2

           end do
!
!  NC and NJ are constraint nodes in the same constraint.
!  NC-NJ is intersected by an interior constraint arc iff
!  1)  NC LEFT NF2->NF1 and (NJ LEFT NF1->NC and NJ LEFT NC->NF2) or
!  2)  NC .NOT. LEFT NF2->NF1 and (NJ LEFT NF1->NC or NJ LEFT NC->NF2),
!  where NF1, NC, NF2 are consecutive constraint nodes.
!
           if (.not. ncf  .or.  .not. njf) go to 11

           if (nc /= ifrst) then
              nf1 = nc - 1
           else
              nf1 = ilast
           end if
           
           if (nc /= ilast) then
              nf2 = nc + 1
           else
              nf2 = ifrst
           end if
           
           lft1 = (xc-x(nf1))*(yj-y(nf1)) >= (xj-x(nf1))*(yc-y(nf1))
           lft2 = (x(nf2)-xc)*(yj-yc) >= (xj-xc)*(y(nf2)-yc)
           lft12 = (x(nf1)-x(nf2))*(yc-y(nf2)) >= (xc-x(nf2))*(y(nf1)-y(nf2))
           
           if ( (lft1  .and.  lft2)  .or.  (.not. lft12 &
                .and.  (lft1  .or.  lft2)) ) go to 12
           
11         continue

        !end do !ADDED 7/4/2013 [JTH]
!
!  NJ is visible from NC.  Exit the loop with DC = Euclidean
!  distance if J = 1.
!
        if (j == 1) then
           go to 14
        end if
        
        dc = min(dc,ds(j) + sqrt((xc-xj)*(xc-xj) + &
             (yc-yj)*(yc-yj)))
        go to 13
!
!  NJ is not visible from NC or is adjacent to NC.  Initialize DC 
!  with D(N1,NK) + D(NK,NC) if J = 1.
!
12      continue
        
        if (j == 1) then
           dc = ds(k) + sqrt((xc-xk)*(xc-xk) &
                + (yc-yk)*(yc-yk))
        end if
        
13      continue

     end do
!
!  Compare DC with DL.
!
     if ( isw  .and.  dc >= dl) then
        go to 15
     end if
!
!  The first (or a closer) candidate for NL has been encountered.
!
14   continue
     
     nl = nc
     dl = dc
     isw = .true.
     
15   continue
     
     if (lpk /= lpkl) then
        go to 7
     end if
     
  end do
!
!  Unmark the elements of NPTS and store NL and DL.
!
  do k = 1,lm1
     nk = npts(k)
     lend(nk) = -lend(nk)
  end do
  
  npts(l) = nl
  ds(l) = dl
  ier = 0
  
  return

end subroutine getnp




function intsec ( x1, y1, x2, y2, x3, y3, x4, y4 )

!*****************************************************************************80
!
!! INTSEC determines if two line segments intersect.
!
!  Discussion:
!
!    Given a pair of line segments P1-P2 and P3-P4, this
!    function returns the value .TRUE. if and only if P1-P2
!    shares one or more points with P3-P4.  The line segments
!    include their endpoints, and the four points need not be
!    distinct.  Thus, either line segment may consist of a
!    single point, and the segments may meet in a V (which is
!    treated as an intersection).  Note that an incorrect
!    decision may result from floating point error if the four
!    endpoints are nearly collinear.
!
!  Modified:
!
!    19 May 2005
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
!    Input, real ( kind = 8 ) X1, Y1 = Coordinates of P1.
!
!    Input, real ( kind = 8 ) X2, Y2 = Coordinates of P2.
!
!    Input, real ( kind = 8 ) X3, Y3 = Coordinates of P3.
!
!    Input, real ( kind = 8 ) X4, Y4 = Coordinates of P4.
!
!    Output, logical INTSEC, the logical value defined above.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) d
  real ( kind = 8 ) dx12
  real ( kind = 8 ) dx31
  real ( kind = 8 ) dx34
  real ( kind = 8 ) dy12
  real ( kind = 8 ) dy31
  real ( kind = 8 ) dy34
  logical intsec
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) x3
  real ( kind = 8 ) x4
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2
  real ( kind = 8 ) y3
  real ( kind = 8 ) y4
!
!  Test for overlap between the smallest rectangles that
!  contain the line segments and have sides parallel to
!  the axes.
!
  if (( x1 < x3 .and.  x1 < x4 .and.  x2 < x3 .and.  x2 < x4)  .or. &
      ( x3 < x1 .and.  x1 > x4 .and.  x2 > x3 .and.  x2 > x4)  .or. &
      ( y1 < y3 .and.  y1 < y4 .and.  y2 < y3 .and.  y2 < y4)  .or. &
      ( y3 < y1 .and.  y4 < y1 .and.  y2 > y3 .and.  y2 > y4)) then
    intsec = .false.
    return
  end if
!
!  Compute A = P4-P3 X P1-P3, B = P2-P1 X P1-P3, and
!  D = P2-P1 X P4-P3 (Z components).
!
  dx12 = x2 - x1
  dy12 = y2 - y1
  dx34 = x4 - x3
  dy34 = y4 - y3
  dx31 = x1 - x3
  dy31 = y1 - y3
  a = dx34 * dy31 - dx31 * dy34
  b = dx12 * dy31 - dx31 * dy12
  d = dx12 * dy34 - dx34 * dy12
!
!  D /= 0 and the point of intersection of the lines defined by the line
!  segments is P = P1 + (A/D)*(P2-P1) = P3 + (B/D)*(P4-P3).
!
  if ( d /= 0.0D+00 ) then

    intsec = &
      0.0D+00 <= a / d .and. a / d <= 1.0D+00 .and. &
      0.0D+00 <= b / d .and. b / d <= 1.0D+00
!
!  D == 0 and thus either the line segments are parallel,
!  or one (or both) of them is a single point.
!
  else

    intsec = ( a == 0.0D+00 .and. b == 0.0D+00 )

  end if

  return

end function intsec




function nearnd ( xp, yp, ist, n, x, y, list, lptr, lend, dsq )

!*****************************************************************************80
!
!! NEARND finds the nearest triangulation node to a point.
!
!  Discussion:
!
!    Given a point P in the plane and a Delaunay triangulation created by
!    subroutine TRMESH or TRMSHR, this function returns the index of the 
!    nearest triangulation node to P.
!
!    The algorithm consists of implicitly adding P to the
!    triangulation, finding the nearest neighbor to P, and
!    implicitly deleting P from the triangulation.  Thus, it
!    is based on the fact that, if P is a node in a Delaunay
!    triangulation, the nearest node to P is a neighbor of P.
!
!    Note that the number of candidates for NEARND
!    (neighbors of P) is limited to LMAX defined in
!    the PARAMETER statement below.
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
!    Input, real ( kind = 8 ) XP, YP, the coordinates of the point P.
!
!    Input, integer ( kind = 4 ) IST, the index of a node at which TRFIND begins 
!    the search.  Search time depends on the proximity
!    of this node to P.
!
!    Input, integer ( kind = 4 ) N, the number of nodes in the triangulation.
!
!    Input, real ( kind = 8 ) X(N), Y(N), the coordinates of the nodes.
!
!    Input, integer ( kind = 4 ) LIST(*), LPTR(*), LEND(N), a data structure 
!    defining the triangulation.  Refer to TRMESH.
!
!    Output, real ( kind = 8 ) DSQ, the square of the distance between P and
!    node NEARND.
!
!    Output, integer ( kind = 4 ) NEARND, the index of the nearest node to P, 
!    or 0 if N < 3 or the triangulation data structure is invalid.
!
  implicit none

  integer ( kind = 4 ), parameter :: lmax = 25
  integer ( kind = 4 ) n

  real ( kind = 8 ) cos1
  real ( kind = 8 ) cos2
  real ( kind = 8 ) ds1
  real ( kind = 8 ) dsq
  real ( kind = 8 ) dsr
  real ( kind = 8 ) dx11
  real ( kind = 8 ) dx12
  real ( kind = 8 ) dx21
  real ( kind = 8 ) dx22
  real ( kind = 8 ) dy11
  real ( kind = 8 ) dy12
  real ( kind = 8 ) dy21
  real ( kind = 8 ) dy22
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i3
  integer ( kind = 4 ) ist
  integer ( kind = 4 ) l
  integer ( kind = 4 ) lend(n)
  integer ( kind = 4 ) list(*)
  integer ( kind = 4 ) listp(lmax)
  integer ( kind = 4 ) lp
  integer ( kind = 4 ) lp1
  integer ( kind = 4 ) lp2
  integer ( kind = 4 ) lpl
  integer ( kind = 4 ) lptr(*)
  integer ( kind = 4 ) lptrp(lmax)
  integer ( kind = 4 ) lstptr
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) n3
  integer ( kind = 4 ) nearnd
  integer ( kind = 4 ) nr
  integer ( kind = 4 ) nst
  real ( kind = 8 ) sin1
  real ( kind = 8 ) sin2
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xp
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) yp
!
!  Store local parameters and test for N invalid.
!
  nearnd = 0
  dsq = -1.0D+00

  if ( n < 3 ) then
    return
  end if

  nst = ist

  if ( nst < 1 .or. n < nst ) then
    nst = 1
  end if
!
!  Find a triangle (I1,I2,I3) containing P, or the rightmost
!  (I1) and leftmost (I2) visible boundary nodes as viewed from P.
!
  call trfind ( nst, xp, yp, n, x, y, list, lptr, lend, i1, i2, i3 )
!
!  Test for collinear nodes.
!
  if ( i1 == 0 ) then
    return
  end if
!
!  Store the linked list of 'neighbors' of P in LISTP and
!  LPTRP.  I1 is the first neighbor, and 0 is stored as
!  the last neighbor if P is not contained in a triangle.
!  L is the length of LISTP and LPTRP, and is limited to LMAX.
!
  if ( i3 /= 0 ) then

    listp(1) = i1
    lptrp(1) = 2
    listp(2) = i2
    lptrp(2) = 3
    listp(3) = i3
    lptrp(3) = 1
    l = 3

  else

    n1 = i1
    l = 1
    lp1 = 2
    listp(l) = n1
    lptrp(l) = lp1
!
!  Loop on the ordered sequence of visible boundary nodes
!  N1 from I1 to I2.
!
    do

      lpl = lend(n1)
      n1 = -list(lpl)
      l = lp1
      lp1 = l+1
      listp(l) = n1
      lptrp(l) = lp1

      if ( n1 == i2 .or. lmax <= lp1 ) then
        exit
      end if

    end do

    l = lp1
    listp(l) = 0
    lptrp(l) = 1

  end if
!
!  Initialize variables for a loop on arcs N1-N2 opposite P
!  in which new 'neighbors' are 'swapped' in.  N1 follows
!  N2 as a neighbor of P, and LP1 and LP2 are the LISTP
!  indexes of N1 and N2.
!
  lp2 = 1
  n2 = i1
  lp1 = lptrp(1)
  n1 = listp(lp1)
!
!  Begin loop:  find the node N3 opposite N1->N2.
!
  do

    lp = lstptr ( lend(n1), n2, list, lptr )

    if ( list(lp) < 0 ) then
      go to 4
    end if

    lp = lptr(lp)
    n3 = abs ( list(lp) )
!
!  Swap test:  Exit the loop if L = LMAX.
!
    if ( lmax <= l ) then
      exit
    end if

    dx11 = x(n1) - x(n3)
    dx12 = x(n2) - x(n3)
    dx22 = x(n2) - xp
    dx21 = x(n1) - xp

    dy11 = y(n1) - y(n3)
    dy12 = y(n2) - y(n3)
    dy22 = y(n2) - yp
    dy21 = y(n1) - yp

    cos1 = dx11 * dx12 + dy11 * dy12
    cos2 = dx22 * dx21 + dy22 * dy21

    if ( 0.0D+00 <= cos1 .and. cos2 >= 0.0D+00 ) then
      go to 4
    end if

    if ( cos1 < 0.0D+00 .and. cos2 < 0.0D+00 ) then
      go to 3
    end if

    sin1 = dx11 * dy12 - dx12 * dy11
    sin2 = dx22 * dy21 - dx21 * dy22

    if ( sin1 * cos2 + cos1 * sin2 >= 0.0D+00 ) then
      go to 4
    end if
!
!  Swap:  Insert N3 following N2 in the adjacency list for P.
!  The two new arcs opposite P must be tested.
!
3   continue

    l = l+1
    lptrp(lp2) = l
    listp(l) = n3
    lptrp(l) = lp1
    lp1 = l
    n1 = n3
    cycle
!
!  No swap:  Advance to the next arc and test for termination
!  on N1 = I1 (LP1 = 1) or N1 followed by 0.
!
4   continue

    if ( lp1 == 1 ) then
      exit
    end if

    lp2 = lp1
    n2 = n1
    lp1 = lptrp(lp1)
    n1 = listp(lp1)

    if ( n1 == 0 ) then
      exit
    end if

  end do
!
!  Set NR and DSR to the index of the nearest node to P and
!  its squared distance from P, respectively.
!
  nr = i1
  dsr = ( x(nr) - xp )**2 + ( y(nr) - yp )**2

  do lp = 2, l

    n1 = listp(lp)

    if ( n1 == 0 ) then
      cycle
    end if

    ds1 = ( x(n1) - xp )**2 + ( y(n1) - yp )**2

    if ( ds1 < dsr ) then
      nr = n1
      dsr = ds1
    end if

  end do

  dsq = dsr
  nearnd = nr

  return

end function nearnd
