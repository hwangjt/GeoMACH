subroutine trlist ( ncc, lcc, n, list, lptr, lend, nrow, nt, ltri, lct, ier )

!*****************************************************************************80
!
!! TRLIST converts a triangulation to triangle list form.
!
!  Discussion:
!
!    This subroutine converts a triangulation data structure
!    from the linked list created by subroutine TRMESH or
!    TRMSHR to a triangle list.
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
!    Input, integer ( kind = 4 ) LIST(*), LPTR(*), LEND(N), linked list data structure 
!    defining the triangulation.  Refer to subroutine TRMESH.
!
!    Input, integer ( kind = 4 ) NROW, the number of rows (entries per triangle) 
!    reserved for the triangle list LTRI.  The value must be 6 if only 
!    the vertex indexes and neighboring triangle indexes are to be
!    stored, or 9 if arc indexes are also to be assigned and stored.  
!    Refer to LTRI.
!
!    Input, integer ( kind = 4 ) LTRI(NROW*NT), where NT is at most 2N-5.  (A 
!    sufficient length is 12 * N if NROW=6 or 18*N if NROW=9.)
!
!    Output, integer ( kind = 4 ) NT, the number of triangles in the triangulation unless
!    IER /= 0, in which case NT = 0.  NT = 2N - NB- 2, where NB is the number 
!    of boundary nodes.
!
!    Output, integer ( kind = 4 ) LTRI(NROW,NT), whose J-th column contains the vertex nodal
!    indexes (first three rows), neighboring triangle indexes (second three
!    rows), and, if NROW = 9, arc indexes (last three rows) associated with
!    triangle J for J = 1,...,NT.  The vertices are ordered counterclockwise
!    with the first vertex taken to be the one with smallest index.  Thus,
!    LTRI(2,J) and LTRI(3,J) are larger than LTRI(1,J) and index adjacent
!    neighbors of node LTRI(1,J).  For I = 1,2,3, LTRI(I+3,J) and LTRI(I+6,J)
!    index the triangle and arc, respectively, which are opposite (not shared
!    by) node LTRI(I,J), with LTRI(I+3,J) = 0 if LTRI(I+6,J) indexes a boundary
!    arc.  Vertex indexes range from 1 to N, triangle indexes from 0 to NT,
!    and, if included, arc indexes from 1 to NA = NT+N-1.  The triangles are 
!    ordered on first (smallest) vertex indexes, except that the sets of
!    constraint triangle (triangles contained in the closure of a constraint
!    region) follow the non-constraint triangles.
!
!    Output, integer ( kind = 4 ) LCT(NCC), containing the triangle index of the first
!    triangle of constraint J in LCT(J).  Thus, the number of non-constraint
!    triangles is LCT(1)-1, and constraint J contains LCT(J+1)-LCT(J) 
!    triangles, where LCT(NCC+1) = NT+1.
!
!    Output, integer ( kind = 4 ) IER = Error indicator.
!    0, if no errors were encountered.
!    1, if NCC, N, NROW, or an LCC entry is outside its valid range on input.
!    2, if the triangulation data structure (LIST,LPTR,LEND) is invalid.  
!
!  Local Parameters:
!
!    ARCS = TRUE iff arc indexes are to be stored.
!    KA,KT = Numbers of currently stored arcs and triangles.
!    N1ST = Starting index for the loop on nodes (N1ST = 1 on
!           pass 1, and N1ST = LCC1 on pass 2).
!    NM2 = Upper bound on candidates for N1.
!    PASS2 = TRUE iff constraint triangles are to be stored.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nrow

  logical arcs
  logical cstri
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i3
  integer ( kind = 4 ) ier
  integer ( kind = 4 ) isv
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jlast
  integer ( kind = 4 ) ka
  integer ( kind = 4 ) kn
  integer ( kind = 4 ) kt
  integer ( kind = 4 ) l
  integer ( kind = 4 ) lcc(*)
  integer ( kind = 4 ) lcc1
  integer ( kind = 4 ) lct(*)
  integer ( kind = 4 ) lend(n)
  integer ( kind = 4 ) list(*)
  integer ( kind = 4 ) lp
  integer ( kind = 4 ) lp2
  integer ( kind = 4 ) lpl
  integer ( kind = 4 ) lpln1
  integer ( kind = 4 ) lptr(*)
  integer ( kind = 4 ) ltri(nrow,*)
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n1st
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) n3
  integer ( kind = 4 ) ncc
  integer ( kind = 4 ) nm2
  integer ( kind = 4 ) nn
  integer ( kind = 4 ) nt
  logical pass2
!
!  Test for invalid input parameters and store the index
!  LCC1 of the first constraint node (if any).
!
  nn = n

  if ( ncc < 0 .or. ( nrow /= 6  .and. nrow /= 9 ) ) then
    nt = 0
    ier = 1
    return
  end if

  lcc1 = nn+1

  if (ncc == 0) then

    if ( nn < 3 ) then
      nt = 0
      ier = 1
      return
    end if

  else

    do i = ncc, 1, -1
      if ( lcc1 - lcc(i) < 3 ) then
        nt = 0
        ier = 1
        return
      end if
      lcc1 = lcc(i)
    end do

    if ( lcc1 < 1 ) then
      nt = 0
      ier = 1
      return
    end if

  end if
!
!  Initialize parameters for loop on triangles KT = (N1,N2,
!  N3), where N1 < N2 and N1 < N3.  This requires two
!  passes through the nodes with all non-constraint
!  triangles stored on the first pass, and the constraint
!  triangles stored on the second.
!
  arcs = nrow == 9
  ka = 0
  kt = 0
  n1st = 1
  nm2 = nn - 2
  pass2 = .false.
!
!  Loop on nodes N1:  
!  J = constraint containing N1,
!  JLAST = last node in constraint J.
!
2 continue

  j = 0
  jlast = lcc1 - 1

  do n1 = n1st, nm2

    if ( jlast < n1 ) then
!
!  N1 is the first node in constraint J+1.  Update J and
!  JLAST, and store the first constraint triangle index
!  if in pass 2.
!
      j = j + 1

      if ( j < ncc ) then
        jlast = lcc(j+1) - 1
      else
        jlast = nn
      end if

      if ( pass2 ) then
        lct(j) = kt + 1
      end if

    end if
!
!  Loop on pairs of adjacent neighbors (N2,N3).  LPLN1 points
!  to the last neighbor of N1, and LP2 points to N2.
!
    lpln1 = lend(n1)
    lp2 = lpln1

    3 continue

      lp2 = lptr(lp2)
      n2 = list(lp2)
      lp = lptr(lp2)
      n3 = abs ( list(lp) )

      if ( n2 < n1 .or. n3 < n1 ) then
        go to 10
      end if
!
!  (N1,N2,N3) is a constraint triangle iff the three nodes
!  are in the same constraint and N2 < N3.  Bypass con-
!  straint triangles on pass 1 and non-constraint triangles
!  on pass 2.
!
      cstri = n1 >= lcc1  .and.  n2 < n3  .and. n3 <= jlast

      if ( ( cstri  .and.  .not. pass2 )  .or. &
          ( .not. cstri  .and.  pass2 ) ) then
        go to 10
      end if
!
!  Add a new triangle KT = (N1,N2,N3).
!
      kt = kt + 1
      ltri(1,kt) = n1
      ltri(2,kt) = n2
      ltri(3,kt) = n3
!
!  Loop on triangle sides (I1,I2) with neighboring triangles
!  KN = (I1,I2,I3).
!
      do i = 1,3

        if ( i == 1 ) then
          i1 = n3
          i2 = n2
        else if ( i == 2 ) then
          i1 = n1
          i2 = n3
        else
          i1 = n2
          i2 = n1
        end if
!
!  Set I3 to the neighbor of I1 which follows I2 unless
!  I2->I1 is a boundary arc.
!
        lpl = lend(i1)
        lp = lptr(lpl)

4       continue

          if (list(lp) == i2) then
            go to 5
          end if

          lp = lptr(lp)

          if ( lp /= lpl ) then
            go to 4
          end if
!
!  I2 is the last neighbor of I1 unless the data structure
!  is invalid.  Bypass the search for a neighboring
!  triangle if I2->I1 is a boundary arc.
!
        if ( abs ( list(lp) ) /= i2 ) then
          go to 13
        end if

        kn = 0

        if (list(lp) < 0) then
          go to 8
        end if
!
!  I2->I1 is not a boundary arc, and LP points to I2 as
!  a neighbor of I1.
!
5   continue

        lp = lptr(lp)
        i3 = abs ( list(lp) )
!
!  Find L such that LTRI(L,KN) = I3 (not used if KN > KT),
!  and permute the vertex indexes of KN so that I1 is
!  smallest.
!
        if ( i1 < i2  .and.  i1 < i3 ) then
          l = 3
        else if (i2 < i3) then
          l = 2
          isv = i1
          i1 = i2
          i2 = i3
          i3 = isv
        else
          l = 1
          isv = i1
          i1 = i3
          i3 = i2
          i2 = isv
        end if
!
!  Test for KN > KT (triangle index not yet assigned).
!
        if ( i1 > n1  .and.  .not. pass2 ) then
          go to 9
        end if
! 
!  Find KN, if it exists, by searching the triangle list in
!  reverse order.
!
        do kn = kt-1,1,-1
          if ( ltri(1,kn) == i1  .and.  ltri(2,kn) == &
              i2 .and. ltri(3,kn) == i3 ) then
            go to 7
          end if
        end do

        go to 9
!
!  Store KT as a neighbor of KN.
!
7       continue

        ltri(l+3,kn) = kt
!
!  Store KN as a neighbor of KT, and add a new arc KA.
!
8       continue

        ltri(i+3,kt) = kn

        if (arcs) then
          ka = ka + 1
          ltri(i+6,kt) = ka
          if ( kn /= 0 ) then
            ltri(l+6,kn) = ka
          end if
        end if

9       continue

    end do
! 
!  Bottom of loop on triangles.
!
10  continue

    if ( lp2 /= lpln1 ) then
      go to 3
    end if

  end do
!
!  Bottom of loop on nodes.
!
  if ( .not. pass2 .and. 0 < ncc ) then
    pass2 = .true.
    n1st = lcc1
    go to 2
  end if
!
!  No errors encountered.
!
  nt = kt
  ier = 0
  return
!
!  Invalid triangulation data structure:  I1 is a neighbor of
!  I2, but I2 is not a neighbor of I1.
!
   13 continue

  nt = 0
  ier = 2

  return

end subroutine trlist




subroutine trlprt ( ncc, lct, n, x, y, nrow, nt, ltri, prntx )

!*****************************************************************************80
!
!! TRLPRT prints the triangles in a triangulation.
!
!  Discussion:
!
!    Given a triangulation of a set of points in the plane,
!    this subroutine prints the triangle list created by
!    subroutine TRLIST and, optionally, the nodal coordinates
!    on logical unit LOUT.  The numbers of boundary nodes,
!    triangles, and arcs, and the constraint region triangle
!    indexes, if any, are also printed.
!
!    All parameters other than PRNTX should be
!    unaltered from their values on output from TRLIST.
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
!    Input, integer ( kind = 4 ) NCC, the number of constraints.
!
!    Input, integer ( kind = 4 ) LCT(NCC), the list of constraint triangle 
!    starting indexes (or dummy array of length 1 if NCC = 0).
!
!    Input, integer ( kind = 4 ) N, the number of nodes in the triangulation.
!    3 <= N <= 9999.
!
!    Input, real ( kind = 8 ) X(N), Y(N), the coordinates of the nodes in the 
!    triangulation; not used unless PRNTX = TRUE.
!
!    Input, integer ( kind = 4 ) NROW, the number of rows (entries per triangle) 
!    reserved for the triangle list LTRI.  The value must be 6 if only 
!    the vertex indexes and neighboring triangle indexes are stored, 
!    or 9 if arc indexes are also stored.
!
!    Input, integer ( kind = 4 ) NT, the number of triangles in the 
!    triangulation.  1 <= NT <= 9999.
!
!    Input, integer ( kind = 4 ) LTRI(NROW,NT), array whose J-th column contains
!    the vertex nodal indexes (first three rows), neighboring triangle 
!    indexes (second three rows), and, if NROW = 9, arc indexes (last
!    three rows) associated with triangle J for J = 1,...,NT.
!
!    Input, logical PRNTX, is TRUE if and only if X and Y are to be printed.
!
!  Local parameters:
!
!    I = DO-loop, nodal index, and row index for LTRI
!    K = DO-loop and triangle index
!    LUN = Logical unit number for output
!    NA = Number of triangulation arcs
!    NB = Number of boundary nodes
!    NL = Number of lines printed on the current page
!    NLMAX = Maximum number of print lines per page
!    NMAX = Maximum value of N and NT (4-digit format)
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nrow
  integer ( kind = 4 ) nt

  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  integer ( kind = 4 ) lct(*)
  integer ( kind = 4 ) ltri(nrow,nt)
  integer ( kind = 4 ) na
  integer ( kind = 4 ) nb
  integer ( kind = 4 ) nl
  integer ( kind = 4 ), parameter :: nlmax = 60
  integer ( kind = 4 ), parameter :: nmax = 9999
  integer ( kind = 4 ) ncc
  logical prntx
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)
!
!  Print a heading and test for invalid input.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TRIPACK (TRLIST) Output:'

  nl = 1
!
!  Print an error message and bypass the loops.
!
  if ( n < 3  .or.  n > nmax  .or. &
      (nrow /= 6  .and.  nrow /= 9)  .or. &
      nt < 1  .or.  nt > nmax) then
    write (*,110) n, nrow, nt
    go to 3
  end if
!
!  Print X and Y.
!
  if (prntx) then

    write (*,101)
    nl = 6

    do i = 1,n
      if (nl >= nlmax) then
        write (*,106)
        nl = 0
      end if
      write (*,102) i, x(i), y(i)
      nl = nl + 1
    end do

  end if
!
!  Print the triangulation LTRI.
!
  if ( nlmax/2 < nl ) then
    write (*,106)
    nl = 0
  end if

  if ( nrow == 6 ) then
    write (*,103)
  else
    write (*,104)
  end if

  nl = nl + 5

  do k = 1, nt
    if ( nlmax <= nl ) then
      write (*,106)
      nl = 0
    end if
    write (*,105) k, (ltri(i,k), i = 1,nrow)
    nl = nl + 1
  end do
!
!  Print NB, NA, and NT (boundary nodes, arcs, and triangles).
!
  nb = 2 * n - nt - 2
  na = nt + n - 1

  if ( nlmax-6 < nl ) then
    write (*,106)
  end if

  write (*,107) nb, na, nt
!
!  Print NCC and LCT.
!
3 continue

  write ( *, '(a,i3)' ) '  Number of constraint curves, NCC = ', ncc

  if ( 0 < ncc ) then
    write (*,109) lct(1:ncc)
  end if

  return
!
!  Print formats:
!
  101 format (//16x,'node',7x,'x(node)',10x,'y(node)'//)
  102 format (16x,i4,2e17.6)
  103 format (//1x,'triangle',8x,'vertices',12x,'neighbors'/ &
          4x,'kt',7x,'n1',5x,'n2',5x,'n3',4x,'kt1',4x, &
          'kt2',4x,'kt3'/)
  104 format (//1x,'triangle',8x,'vertices',12x,'neighbors', &
          14x,'arcs'/ &
          4x,'kt',7x,'n1',5x,'n2',5x,'n3',4x,'kt1',4x, &
          'kt2',4x,'kt3',4x,'ka1',4x,'ka2',4x,'ka3'/)
  105 format (2x,i4,2x,6(3x,i4),3(2x,i5))
  106 format (///)
  107 format (/1x,'nb = ',i4,' boundary nodes',5x, &
          'na = ',i5,' arcs',5x,'nt = ',i5, &
          ' triangles')
  109 format (1x,9x,14i5)
  110 format (//1x,10x,'*** invalid parameter:  n =',i5, &
          ', nrow =',i5,', nt =',i5,' ***')

end subroutine trlprt




subroutine trmtst ( n, x, y, list, lptr, lend, lnew, tol, armax, ier )

!*****************************************************************************80
!
!! TRMTST tests a data structure representing a Delaunay triangulation.
!
!  Discussion:
!
!    This subroutine tests the validity of the data structure
!    representing a Delaunay triangulation created by subrou-
!    tine TRMESH.  The following properties are tested:
!
!    1)  Each interior node has at least three neighbors, and
!        each boundary node has at least two neighbors.
!
!    2)  abs ( LIST(LP) ) is a valid nodal index in the range
!        1 to N and LIST(LP) > 0 unless LP = LEND(K) for some
!        nodal index K.
!
!    3)  Each pointer LEND(K) for K = 1 to N and LPTR(LP) for
!        LP = 1 to LNEW-1 is a valid LIST index in the range
!        1 to LNEW-1.
!
!    4)  N .GE. NB .GE. 3, NT = 2*N-NB-2, and NA = 3*N-NB-3 =
!        (LNEW-1)/2, where NB, NT, and NA are the numbers of
!        boundary nodes, triangles, and arcs, respectively.
!
!    5)  Each circumcircle defined by the vertices of a tri-
!        angle contains no nodes in its interior.  This prop-
!        erty distinguishes a Delaunay triangulation from an
!        arbitrary triangulation of the nodes.
!
!    Note that no test is made for the property that a triangulation 
!    covers the convex hull of the nodes.
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
!    Input, integer ( kind = 4 ) N, the number of nodes.  N .GE. 3.
!
!    Input, real ( kind = 8 ) X(N), Y(N), the nodal coordinates.
!
!    Input, integer ( kind = 4 ) LIST(*), LPTR(*), LEND(N), the data structure 
!    containing the triangulation.  Refer to subroutine TRMESH.
!
!    Input, real ( kind = 8 ) TOL, nonnegative tolerance to allow for
!    floating-point errors in the circumcircle test.  An error situation
!    is defined as 
!      (R**2 - D**2) / R**2 > TOL, 
!    where R is the radius of a circumcircle and D is the distance from the
!    circumcenter to the nearest node.  A reasonable value for TOL is 
!    10*EPS, where EPS is the machine precision.  The test is effectively
!    bypassed by making TOL large.  If TOL < 0, the tolerance is taken 
!    to be 0.
!
!    Output, real ( kind = 8 ) ARMAX, maximum aspect ratio (radius of inscribed
!    circle divided by circumradius) of a triangle in the triangulation 
!    unless 0 < IER.
!
!    Output, integer ( kind = 4 ) IER = Error indicator:
!    -1 if one or more null triangles (area = 0) are present but no (other)
!      errors were encountered.  A null triangle is an error only if it 
!      occurs in the the interior.
!     0 if no errors or null triangles were encountered.
!     1 if a node has too few neighbors.
!     2 if a LIST entry is outside its valid range.
!     3 if a LPTR or LEND entry is outside its valid range.
!     4 if the triangulation parameters (N, NB, NT, NA, and LNEW) are
!      inconsistent (or N < 3 or LNEW is invalid).
!     5 if a triangle contains a node interior to its circumcircle.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) ar
  real ( kind = 8 ) armax
  real ( kind = 8 ) cr
  real ( kind = 8 ) cx
  real ( kind = 8 ) cy
  integer ( kind = 4 ) ier
  integer ( kind = 4 ) k
  integer ( kind = 4 ) lend(n)
  integer ( kind = 4 ) list(*)
  integer ( kind = 4 ) lnew
  integer ( kind = 4 ) lp
  integer ( kind = 4 ) lpl
  integer ( kind = 4 ) lpmax
  integer ( kind = 4 ) lpn
  integer ( kind = 4 ) lptr(*)
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) n3
  integer ( kind = 4 ) na
  integer ( kind = 4 ) nb
  integer ( kind = 4 ) nfail
  integer ( kind = 4 ) nn
  integer ( kind = 4 ) nnb
  integer ( kind = 4 ) nt
  integer ( kind = 4 ) null
  logical ratio
  real ( kind = 8 ) rs
  real ( kind = 8 ) rtol
  real ( kind = 8 ) sa
  real ( kind = 8 ) tol
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)
!
!  Store local variables, test for errors in input, and
!  initialize counts.
!
  nn = n
  lpmax = lnew - 1
  rtol = tol
  rtol = max ( rtol, 0.0D+00 )
  ratio = .true.
  armax = 0.0D+00

  if ( nn < 3 ) then
    go to 14
  end if

  nb = 0
  nt = 0
  null = 0
  nfail = 0
!
!  Loop on triangles (N1,N2,N3) such that N2 and N3 index
!  adjacent neighbors of N1 and are both larger than N1
!  (each triangle is associated with its smallest index).
!  NNB is the neighbor count for N1.
!
  do n1 = 1,nn

    nnb = 0
    lpl = lend(n1)

    if ( lpl < 1  .or.  lpmax < lpl ) then
      lp = lpl
      go to 13
    end if

    lp = lpl
!
!  Loop on neighbors of N1.
!
1   continue

      lp = lptr(lp)
      nnb = nnb + 1

      if (lp < 1  .or.  lp > lpmax) then
        go to 13
      end if

      n2 = list(lp)

      if (n2 < 0) then
        if (lp /= lpl) then
          go to 12
        end if
        if (n2 == 0  .or.  -n2 > nn) go to 12
        nb = nb + 1
        go to 4
      end if

      if (n2 < 1  .or.  n2 > nn) then
        go to 12
      end if

      lpn = lptr(lp)
      n3 = abs ( list(lpn) )
      if (n2 < n1  .or.  n3 < n1) go to 4
      nt = nt + 1
!
!  Compute the coordinates of the circumcenter of (N1,N2,N3).
!
      call circum (x(n1),y(n1),x(n2),y(n2),x(n3),y(n3), &
                   ratio, cx,cy,cr,sa,ar)
      if (sa == 0.) then
        null = null + 1
        go to 4
      end if

      armax = max(armax,ar)
!
!  Test for nodes within the circumcircle.
!
      rs = cr*cr*(1.-rtol)

      do k = 1,nn
        if ( k == n1  .or.  k == n2  .or. &
            k == n3 ) go to 2
        if ((cx-x(k))**2 + (cy-y(k))**2 < rs) go to 3
2       continue
      end do

      go to 4
!
!  Node K is interior to the circumcircle of (N1,N2,N3).
!
3     continue

      nfail = nfail + 1
!
!  Bottom of loop on neighbors.
!
4     continue

      if (lp /= lpl) go to 1
    if (nnb < 2  .or.  (nnb == 2  .and. &
        list(lpl) > 0)) go to 11

  end do
!
!  Test parameters for consistency and check for NFAIL = 0.
!
  na = lpmax/2
  if (nb < 3  .or.  nt /= 2 * nn - nb - 2  .or. &
      na /= 3*nn-nb-3) go to 14
  if (nfail /= 0) go to 15
!
!  No errors were encountered.
!
  ier = 0
  if (null == 0) return
  ier = -1
  write (*,100) null
  100 format (//5x,'*** trmtst -- ',i5,' null triangles ', &
          'are present'/19x,'(null triangles ', &
          'on the boundary are unavoidable) ***'//)
  return
!
!  Node N1 has fewer than three neighbors.
!
11 continue

  ier = 1
  write (*,110) n1, nnb
  110 format (//5x,'*** trmtst -- node ',i5, &
          ' has only ',i5,' neighbors ***'/)
  return
!
!  N2 = LIST(LP) is outside its valid range.
!
12 continue

  ier = 2
  write (*,120) n2, lp, n1
  120 format (//5x,'*** trmtst -- list(lp) =',i5, &
          ', for lp =',i5,','/19x, &
  'is not a valid neighbor of ',i5,' ***'/)
  return
!
!  LIST pointer LP is outside its valid range.
!
13 continue

  ier = 3
  write (*,130) lp, lnew, n1
  130 format (//5x,'*** trmtst -- lp =',i5,' is not in the', &
          ' range 1 to lnew-1 for lnew = ',i5/ &
          19x,'lp points to a neighbor of ',i5, &
          ' ***'/)
  return
!
!  Inconsistent triangulation parameters encountered.
!
14 continue

  ier = 4
  write (*,140) n, lnew, nb, nt, na
  140 format (//5x,'*** trmtst -- inconsistent parameters', &
          ' ***'/19x,'n = ',i5,' nodes',12x,'lnew =',i5/ &
          19x,'nb = ',i5,' boundary nodes'/ &
          19x,'nt = ',i5,' triangles'/ &
          19x,'na = ',i5,' arcs'/)
  return
!
!  Circumcircle test failure.
!
15 continue

  ier = 5
  write (*,150) nfail
  150 format (//5x,'*** trmtst -- ',i5,' circumcircles ', &
          'contain nodes in their interiors ***'/)
  return

end subroutine trmtst




subroutine trplot ( lun, pltsiz, wx1, wx2, wy1, wy2, ncc, lcc, &
  n, x, y, list, lptr, lend, title, numbr, ier )

!*****************************************************************************80
!
!! TRPLOT plots a triangulation in an EPS file.
!
!  Discussion:
!
!    This subroutine creates a level-2 Encapsulated Postscript (EPS) file 
!    containing a triangulation plot.
!
!    Various plotting options can be controlled by altering
!    the data statement below.
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
!    Input, integer ( kind = 4 ) LUN, the logical unit number in the range 0 to 99.
!    The unit should be opened with an appropriate
!    file name before the call to this routine.
!
!    Input, real ( kind = 8 ) PLTSIZ, the plot size in inches.  The window is
!    mapped, with aspect ratio preserved, to a rectangular viewport with maximum
!    side-length equal to .88*PLTSIZ (leaving room for labels outside 
!    the viewport).  The viewport is centered on the 8.5 by 11 inch page, 
!    and its boundary is drawn.  1.0 <= PLTSIZ <= 8.5.
!
!    Input, real ( kind = 8 ) WX1, WX2, WY1, WY2, parameters defining a
!    rectangular window against which the triangulation is clipped.  (Only the
!    portion of the triangulation that lies in the window is drawn.)
!    (WX1,WY1) and (WX2,WY2) are the lower left and upper right 
!    corners, respectively.  WX1 < WX2 and WY1 < WY2.
!
!    Input, integer ( kind = 4 ) NCC, the number of constraint curves.  Refer to 
!    subroutine ADDCST.  NCC >= 0.
!
!    Input, integer ( kind = 4 ) LCC(NCC) (or dummy parameter if NCC = 0) containing the
!    index of the first node of constraint I in LCC(I).  For I = 1 to
!    NCC, LCC(I+1)-LCC(I) >= 3, where LCC(NCC+1) = N+1.
!
!    Input, integer ( kind = 4 ) N, the number of nodes in the triangulation.  
!    3 <= N.
!
!    Input, real ( kind = 8 ) X(N), Y(N), the coordinates of the nodes with
!    non-constraint nodes in the first LCC(1)-1 locations.
!
!    Input, integer ( kind = 4 ) LIST(*), LPTR(*), LEND(N), data structure 
!    defining the triangulation.  Refer to subroutine TRMESH.
!
!    Input, character ( len = * ) TITLE, a string to be centered above the 
!    plot.  The string must be enclosed in parentheses; i.e., the first and 
!    last characters must be '(' and ')', respectively, but these are not
!    displayed.  TITLE may have at most 80 characters including the parentheses.
!
!    Input, logical NUMBR, option indicator:  If NUMBR = TRUE, the
!    nodal indexes are plotted next to the nodes.
!
!    Output, integer ( kind = 4 ) IER = Error indicator:
!    0, if no errors were encountered.
!    1, if LUN, PLTSIZ, NCC, or N is outside its valid range.  LCC is not 
!      tested for validity.
!    2, if WX1 >= WX2 or WY1 >= WY2.
!    3, if an error was encountered in writing to unit LUN.
!
!  Local parameters:
!
!    ANNOT =     Logical variable with value TRUE iff the plot
!                is to be annotated with the values of WX1,
!                WX2, WY1, and WY2
!    CNSTR       Logical variable used to flag constraint arcs:
!                TRUE iff N0-N1 lies in a constraint region
!    DASHL =     Length (in points, at 72 points per inch) of
!                dashes and spaces in a dashed line pattern
!                used for drawing constraint arcs
!    DX =        Window width WX2-WX1
!    DY =        Window height WY2-WY1
!    FSIZN =     Font size in points for labeling nodes with
!                their indexes if NUMBR = TRUE
!    FSIZT =     Font size in points for the title (and
!                annotation if ANNOT = TRUE)
!    I =         Constraint index (1 to NCC)
!    IFRST =     Index of the first node in constraint I
!    IH =        Height of the viewport in points
!    ILAST =     Index of the last node in constraint I
!    IPX1,IPY1 = X and y coordinates (in points) of the lower
!                left corner of the bounding box or viewport
!    IPX2,IPY2 = X and y coordinates (in points) of the upper
!                right corner of the bounding box or viewport
!    IW =        Width of the viewport in points
!    LP =        LIST index (pointer)
!    LPL =       Pointer to the last neighbor of N0
!    N0 =        Nodal index and DO-loop index
!    N0BAK =     Predecessor of N0 in a constraint curve
!                (sequence of adjacent constraint nodes)
!    N0FOR =     Successor to N0 in a constraint curve
!    N1 =        Index of a neighbor of N0
!    NLS =       Index of the last non-constraint node
!    PASS1 =     Logical variable used to flag the first pass
!                through the constraint nodes
!    R =         Aspect ratio DX/DY
!    SFX,SFY =   Scale factors for mapping world coordinates
!                (window coordinates in [WX1,WX2] X [WY1,WY2])
!                to viewport coordinates in [IPX1,IPX2] X [IPY1,IPY2]
!    T =         Temporary variable
!    TX,TY =     Translation vector for mapping world coordi-
!                nates to viewport coordinates
!    X0,Y0 =     X(N0),Y(N0) or label location
!
  implicit none

  integer ( kind = 4 ) n

  logical, parameter :: annot = .true.
  logical cnstr
  real ( kind = 8 ), parameter :: dashl = 4.0D+00
  real ( kind = 8 ) dx
  real ( kind = 8 ) dy
  real ( kind = 8 ), parameter :: fsizn = 10.0D+00
  real ( kind = 8 ), parameter :: fsizt = 16.0D+00
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ier
  integer ( kind = 4 ) ifrst
  integer ( kind = 4 ) ih
  integer ( kind = 4 ) ilast
  integer ( kind = 4 ) ipx1
  integer ( kind = 4 ) ipx2
  integer ( kind = 4 ) ipy1
  integer ( kind = 4 ) ipy2
  integer ( kind = 4 ) iw
  integer ( kind = 4 ) lcc(*)
  integer ( kind = 4 ) lend(n)
  integer ( kind = 4 ) list(*)
  integer ( kind = 4 ) lp
  integer ( kind = 4 ) lpl
  integer ( kind = 4 ) lptr(*)
  integer ( kind = 4 ) lun
  integer ( kind = 4 ) n0
  integer ( kind = 4 ) n0bak
  integer ( kind = 4 ) n0for
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) ncc
  integer ( kind = 4 ) nls
  logical numbr
  logical pass1
  real ( kind = 8 ) pltsiz
  real ( kind = 8 ) r
  real ( kind = 8 ) sfx
  real ( kind = 8 ) sfy
  real ( kind = 8 ) t
  character ( len = * ) title
  real ( kind = 8 ) tx
  real ( kind = 8 ) ty
  real ( kind = 8 ) wx1
  real ( kind = 8 ) wx2
  real ( kind = 8 ) wy1
  real ( kind = 8 ) wy2
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) x0
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) y0
!
!  Test for error 1.
!
  if (lun < 0  .or.  lun > 99  .or. &
      pltsiz < 1.0D+00  .or.  pltsiz > 8.5  .or. &
      ncc < 0  .or.  n < 3) then
    ier = 1
    return
  end if
!
!  Set NLS to the last non-constraint node.
!
  if ( ncc > 0 ) then
    nls = lcc(1)-1
  else
    nls = n
  end if
!
!  Compute the aspect ratio of the window.
!
  dx = wx2 - wx1
  dy = wy2 - wy1

  if ( dx <= 0.0D+00 .or. dy <= 0.0D+00 ) then
    ier = 2
    return
  end if

  r = dx / dy
!
!  Compute the lower left (IPX1,IPY1) and upper right
!  (IPX2,IPY2) corner coordinates of the bounding box.
!  The coordinates, specified in default user space units
!  (points, at 72 points/inch with origin at the lower
!  left corner of the page), are chosen to preserve the
!  aspect ratio R, and to center the plot on the 8.5 by 11
!  inch page.  The center of the page is (306,396), and
!  T = PLTSIZ/2 in points.
!
  t = 36.0D+00 * pltsiz

  if ( 1.0D+00 <= r ) then
    ipx1 = 306 - nint(t)
    ipx2 = 306 + nint(t)
    ipy1 = 396 - nint(t/r)
    ipy2 = 396 + nint(t/r)
  else
    ipx1 = 306 - nint(t*r)
    ipx2 = 306 + nint(t*r)
    ipy1 = 396 - nint(t)
    ipy2 = 396 + nint(t)
  end if
!
!  Output header comments.
!
  write (lun,100,err=13) ipx1, ipy1, ipx2, ipy2
  100 format ('%!ps-adobe-3.0 epsf-3.0'/ &
          '%%boundingbox:',4i4/ &
          '%%title:  triangulation'/ &
          '%%creator:  tripack'/ &
          '%%endcomments')
!
!  Set (IPX1,IPY1) and (IPX2,IPY2) to the corner coordinates
!  of a viewport obtained by shrinking the bounding box by
!  12% in each dimension.
!
  iw = nint ( 0.88D+00 * real ( ipx2 - ipx1, kind = 8 ) )
  ih = nint ( 0.88D+00 * real ( ipy2 - ipy1, kind = 8 ) )
  ipx1 = 306 - iw/2
  ipx2 = 306 + iw/2
  ipy1 = 396 - ih/2
  ipy2 = 396 + ih/2
!
!  Set the line thickness to 2 points, and draw the viewport boundary.
!
  t = 2.0D+00
  write (lun,110,err=13) t
  write (lun,120,err=13) ipx1, ipy1
  write (lun,130,err=13) ipx1, ipy2
  write (lun,130,err=13) ipx2, ipy2
  write (lun,130,err=13) ipx2, ipy1
  write (lun,140,err=13)
  write (lun,150,err=13)
  110 format (f12.6,' setlinewidth')
  120 format (2i4,' moveto')
  130 format (2i4,' lineto')
  140 format ('closepath')
  150 format ('stroke')
!
!  Set up a mapping from the window to the viewport.
!
  sfx = real ( iw, kind = 8 ) / dx
  sfy = real ( ih, kind = 8 ) / dy
  tx = ipx1 - sfx*wx1
  ty = ipy1 - sfy*wy1
  write (lun,160,err=13) tx, ty, sfx, sfy
  160 format (2f12.6,' translate'/ &
          2f12.6,' scale')
!
!  The line thickness (believe it or not) must be
!  changed to reflect the new scaling which is applied to
!  all subsequent output.  Set it to 1.0 point.
!
  t = 2.0D+00 / (sfx+sfy)
  write (lun,110,err=13) t
!
!  Save the current graphics state, and set the clip path to
!  the boundary of the window.
!
  write (lun,170,err=13)
  write (lun,180,err=13) wx1, wy1
  write (lun,190,err=13) wx2, wy1
  write (lun,190,err=13) wx2, wy2
  write (lun,190,err=13) wx1, wy2
  write (lun,200,err=13)
  170 format ('gsave')
  180 format (2f12.6,' moveto')
  190 format (2f12.6,' lineto')
  200 format ('closepath clip newpath')
!
!  Draw the edges N0->N1, where N1 > N0, beginning with a
!  loop on non-constraint nodes N0.  LPL points to the
!  last neighbor of N0.
!
  do n0 = 1,nls

    x0 = x(n0)
    y0 = y(n0)
    lpl = lend(n0)
    lp = lpl
!
!  Loop on neighbors N1 of N0.
!
    do

      lp = lptr(lp)
      n1 = abs ( list(lp) )

      if ( n0 < n1 ) then
        write (lun,210,err=13) x0, y0, x(n1), y(n1)
  210       format (2f12.6,' moveto',2f12.6,' lineto')
      end if

      if ( lp == lpl) then
        exit
      end if

    end do

  end do
!
!  Loop through the constraint nodes twice.  The non-constraint arcs 
!  incident on constraint nodes are drawn (with solid lines) on the first 
!  pass, and the constraint arcs (both boundary and interior, if any)
!  are drawn (with dashed lines) on the second pass.
!
  pass1 = .true.
!
!  Loop on constraint nodes N0 with (N0BAK,N0,N0FOR) a subsequence of 
!  constraint I.  The outer loop is on constraints I with first and last 
!  nodes IFRST and ILAST.
!
4 continue

  ifrst = n+1

  do i = ncc, 1, -1

    ilast = ifrst - 1
    ifrst = lcc(i)
    n0bak = ilast

    do n0 = ifrst, ilast

      n0for = n0 + 1
      if (n0 == ilast) n0for = ifrst
      lpl = lend(n0)
      x0 = x(n0)
      y0 = y(n0)
      lp = lpl
!
!  Loop on neighbors N1 of N0.  CNSTR = TRUE iff N0-N1 is a
!  constraint arc.
!
!  Initialize CNSTR to TRUE iff the first neighbor of N0
!  strictly follows N0FOR and precedes or coincides with
!  N0BAK (in counterclockwise order).
!
      do

        lp = lptr(lp)
        n1 = abs ( list(lp) )
        if ( n1 == n0for .or. n1 == n0bak ) then
          exit
        end if

      end do

      cnstr = n1 == n0bak
      lp = lpl
!
!  Loop on neighbors N1 of N0.  Update CNSTR and test for N1 > N0.
!
6 continue

        lp = lptr(lp)
        n1 = abs ( list(lp) )

        if (n1 == n0for) then
          cnstr = .true.
        end if
!
!  Draw the edge iff (PASS1=TRUE and CNSTR=FALSE) or
!  (PASS1=FALSE and CNSTR=TRUE); i.e., CNSTR and PASS1
!  have opposite values.
!
        if ( n0 < n1 ) then

          if ( cnstr .neqv. pass1 ) then
            write (lun,210,err=13) x0, y0, x(n1), y(n1)
          end if

        end if

        if ( n1 == n0bak ) then
          cnstr = .false.
        end if
!
!  Bottom of loops.
!
        if ( lp /= lpl ) then
          go to 6
        end if

      n0bak = n0

    end do

  end do

  if (pass1) then
!
!  End of first pass:  paint the path and change to dashed
!  lines for subsequent drawing.  Since the scale factors
!  are applied to everything, the dash length must be
!  specified in world coordinates.
!
    pass1 = .false.
    write (lun,150,err=13)
    t = dashl * 2.0D+00 / ( sfx + sfy )
    write (lun,220,err=13) t
  220   format ('[',f12.6,'] 0 setdash')
    go to 4

  end if
!
!  Paint the path and restore the saved graphics state (with
!  no clip path).
!
  write (lun,150,err=13)
  write (lun,230,err=13)
  230 format ('grestore')

  if (numbr) then
!
!  Nodes in the window are to be labeled with their indexes.
!  Convert FSIZN from points to world coordinates, and
!  output the commands to select a font and scale it.
!
    t = fsizn * 2.0D+00 / ( sfx + sfy )
    write (lun,240,err=13) t
  240   format ('/Helvetica findfont'/ &
            f12.6,' scalefont setfont')
!
!  Loop on nodes N0 with coordinates (X0,Y0).
!
    do n0 = 1, n

      x0 = x(n0)
      y0 = y(n0)
!
!  Move to (X0,Y0), and draw the label N0.  The first character will 
!  have its lower left corner about one
!  character width to the right of the nodal position.
!
      if ( x0 >= wx1  .and.  x0 <= wx2  .and. &
           y0 >= wy1  .and.  y0 <= wy2 ) then
        write (lun,180,err=13) x0, y0
        write (lun,250,err=13) n0
  250   format ('(',i3,') show')
      end if

    end do

  end if
!
!  Convert FSIZT from points to world coordinates, and output
!  the commands to select a font and scale it.
!
  t = fsizt * 2.0D+00 / ( sfx + sfy )
  write (lun,240,err=13) t
!
!  Display TITLE centered above the plot:
!
  y0 = wy2 + 3.0D+00 * t
  write (lun,260,err=13) title, ( wx1 + wx2 ) / 2.0D+00, y0
  260 format (a80/'  stringwidth pop 2 div neg ',f12.6, &
          ' add ',f12.6,' moveto')
  write (lun,270,err=13) title
  270 format (a80/'  show')
  if (annot) then
!
!  Display the window extrema below the plot.
!
    x0 = wx1
    y0 = wy1 - 100.0D+00 / ( sfx + sfy )
    write (lun,180,err=13) x0, y0
    write (lun,280,err=13) wx1, wx2
    y0 = y0 - 2.0D+00 * t
    write (lun,290,err=13) x0, y0, wy1, wy2
  280   format ('(window:   wx1 = ',e9.3,',   wx2 = ',e9.3, &
            ') show')
  290   format ('(window:  ) stringwidth pop ',f12.6,' add', &
            f12.6,' moveto'/ &
            '( wy1 = ',e9.3,',   wy2 = ',e9.3,') show')
  end if
!
!  Paint the path and output the showpage command and
!  end-of-file indicator.
!
  write (lun,300,err=13)
  300 format ('stroke'/ &
          'showpage'/ &
          '%%eof')
!
!  HP's interpreters require a one-byte End-of-PostScript-Job
!  indicator (to eliminate a timeout error message): ASCII 4.
!
  write (lun,310,err=13) char(4)
  310 format (a1)
!
!  No error encountered.
!
  ier = 0
  return
!
!  Error writing to unit LUN.
!
   13 ier = 3
  return

end subroutine trplot




subroutine trprnt ( ncc, lcc, n, x, y, list, lptr, lend, prntx )

!*****************************************************************************80
!
!! TRPRNT prints information about a planar triangulation.
!
!  Discussion:
!
!    Given a triangulation of a set of points in the plane,
!    this subroutine prints the adjacency lists and, optionally, 
!    the nodal coordinates.  The list of neighbors of a boundary 
!    node is followed by index 0.  The numbers of boundary nodes, 
!    triangles, and arcs, and the constraint curve starting indexes, 
!    if any, are also printed.
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
!    Input, integer ( kind = 4 ) NCC, the number of constraints.
!
!    Input, integer ( kind = 4 ) LCC(*), list of constraint curve starting 
!    indexes (or dummy array of length 1 if NCC = 0).
!
!    Input, integer ( kind = 4 ) N, the number of nodes in the triangulation.
!    3 <= N <= 9999.
!
!    Input, real ( kind = 8 ) X(N), Y(N), the coordinates of the nodes in the 
!    triangulation; not used unless PRNTX = TRUE.
!
!    Input, integer ( kind = 4 ) LIST(*), LPTR(*), LEND(N), data structure 
!    defining the triangulation.  Refer to subroutine TRMESH.
!
!    Input, logical PRNTX, TRUE if and only if X and Y are to be printed.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) k
  integer ( kind = 4 ) lcc(*)
  integer ( kind = 4 ) lend(n)
  integer ( kind = 4 ) list(*)
  integer ( kind = 4 ) lp
  integer ( kind = 4 ) lpl
  integer ( kind = 4 ) lptr(*)
  integer ( kind = 4 ) na
  integer ( kind = 4 ) nabor(100)
  integer ( kind = 4 ) nb
  integer ( kind = 4 ) ncc
  integer ( kind = 4 ) nd
  integer ( kind = 4 ) nl
  integer ( kind = 4 ), parameter :: nlmax = 60
  integer ( kind = 4 ), parameter :: nmax = 9999
  integer ( kind = 4 ) nn
  integer ( kind = 4 ) node
  integer ( kind = 4 ) nt
  logical prntx
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)

  nn = n
!
!  Print a heading and test the range of N.
!
  write ( *,100) nn
  if ( nn < 3 .or. nmax < nn ) then
!
!  N is outside its valid range.
!
    write (*,110)
    go to 5
  end if
!
!  Initialize NL (the number of lines printed on the current
!  page) and NB (the number of boundary nodes encountered).
!
  nl = 6
  nb = 0
  if (.not. prntx) then
!
!  Print LIST only.  K is the number of neighbors of NODE
!  which are stored in NABOR.
!
    write (*,101)

    do node = 1, nn

      lpl = lend(node)
      lp = lpl
      k = 0

      do

        k = k + 1
        lp = lptr(lp)
        nd = list(lp)
        nabor(k) = nd

        if ( lp == lpl ) then
          exit
        end if

      end do

      if ( nd <= 0 ) then
!
!  NODE is a boundary node.  Correct the sign of the last
!  neighbor, add 0 to the end of the list, and increment NB.
!
        nabor(k) = -nd
        k = k + 1
        nabor(k) = 0
        nb = nb + 1
      end if
!
!  Increment NL and print the list of neighbors.
!
      inc = (k-1) / 14 + 2
      nl = nl + inc

      if ( nlmax < nl ) then
        write (*,106)
        nl = inc
      end if

      write (*,103) node, nabor(1:k)
      if ( k /= 14 ) then
        write (*,105)
      end if

    end do

  else
!
!  Print X, Y, and LIST.
!
    write (*,102)

    do node = 1,nn

      lpl = lend(node)
      lp = lpl
      k = 0

      do

        k = k + 1
        lp = lptr(lp)
        nd = list(lp)
        nabor(k) = nd

        if ( lp == lpl ) then
          exit
        end if

      end do

      if ( nd <= 0 ) then
!
!  NODE is a boundary node.
!
        nabor(k) = -nd
        k = k + 1
        nabor(k) = 0
        nb = nb + 1
      end if
!
!  Increment NL and print X, Y, and NABOR.
!
      inc = (k-1) / 8 + 2
      nl = nl + inc

      if ( nlmax < nl ) then
        write (*,106)
        nl = inc
      end if

      write (*,104) node, x(node), y(node), &
                      (nabor(i), i = 1,k)
      if (k /= 8) write (*,105)

    end do

  end if
!
!  Print NB, NA, and NT (boundary nodes, arcs, and triangles).
!
  nt = 2 * nn - nb - 2
  na = nt + nn - 1

  if ( nlmax - 6 < nl ) then
    write (*,106)
  end if

  write (*,107) nb, na, nt
!
!  Print NCC and LCC.
!
5 continue

  write ( *, '(a)' ) ' '
  write ( *, 108 ) ncc
  write ( *, 109 ) lcc(1:ncc)

  return

  100 format (///,26x,'adjacency sets,    n = ',i5//)
  101 format (1x,'node',32x,'neighbors of node'//)
  102 format (1x,'node',5x,'x(node)',8x,'y(node)', &
          20x,'neighbors of node'//)
  103 format (1x,i4,5x,14i5/(1x,9x,14i5))
  104 format (1x,i4,2e15.6,5x,8i5/(1x,39x,8i5))
  105 format (1x)
  106 format (///)
  107 format (/1x,'nb = ',i4,' boundary nodes',5x, &
          'na = ',i5,' arcs',5x,'nt = ',i5, &
          ' triangles')
  108 format ('ncc =',i3,' constraint curves')
  109 format (1x,9x,14i5)
  110 format (1x,10x,'*** N is outside its valid range ***')

end subroutine trprnt
