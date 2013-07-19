subroutine delarc ( n, io1, io2, list, lptr, lend, lnew, ier )

!*****************************************************************************80
!
!! DELARC deletes a boundary arc from a triangulation.
!
!  Discussion:
!
!    This subroutine deletes a boundary arc from a triangula-
!    tion.  It may be used to remove a null triangle from the
!    convex hull boundary.  Note, however, that if the union of
!    triangles is rendered nonconvex, Subroutines DELNOD, EDGE,
!    and TRFIND may fail.  Thus, Subroutines ADDCST, ADDNOD,
!    DELNOD, EDGE, and NEARND should not be called following
!    an arc deletion.
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
!    4 <= N.
!
!    Input, integer ( kind = 4 ) IO1, IO2, the indexes (in the range 1 to N) of
!    a pair of adjacent boundary nodes defining the arc to be removed.
!
!    Input/output, integer ( kind = 4 ) LIST(*), LPTR(*), LEND(N), LNEW, the
!    triangulation data structure created by TRMESH or TRMSHR.
!    On output, updated with the removal of arc IO1-IO2 unless 0 < IER.
!
!    Output, integer ( kind = 4 ) IER = Error indicator:
!    0, if no errors were encountered.
!    1, if N, IO1, or IO2 is outside its valid range, or IO1 = IO2.
!    2, if IO1-IO2 is not a boundary arc.
!    3, if the node opposite IO1-IO2 is already a boundary node, and 
!      thus IO1 or IO2 has only two neighbors or a deletion would result 
!      in two triangulations sharing a single node.
!    4, if one of the nodes is a neighbor of the other, but not vice 
!      versa, implying an invalid triangulation data structure.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) ier
  integer ( kind = 4 ) io1
  integer ( kind = 4 ) io2
  integer ( kind = 4 ) lend(n)
  integer ( kind = 4 ) list(*)
  integer ( kind = 4 ) lnew
  integer ( kind = 4 ) lp
  integer ( kind = 4 ) lph
  integer ( kind = 4 ) lpl
  integer ( kind = 4 ) lptr(*)
  integer ( kind = 4 ) lstptr
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) n3

  n1 = io1
  n2 = io2
!
!  Test for errors, and set N1->N2 to the directed boundary
!  edge associated with IO1-IO2:  (N1,N2,N3) is a triangle
!  for some N3.
!
  if ( n < 4  .or.  n1 < 1  .or.  n < n1 .or. &
      n2 < 1  .or.  n < n2 .or.  n1 == n2 ) then
    ier = 1
    return
  end if

  lpl = lend(n2)

  if ( -list(lpl) /= n1 ) then

    n1 = n2
    n2 = io1
    lpl = lend(n2)

    if ( -list(lpl) /= n1 ) then
      ier = 2
      return
    end if

  end if
!
!  Set N3 to the node opposite N1->N2 (the second neighbor
!  of N1), and test for error 3 (N3 already a boundary node).
!
  lpl = lend(n1)
  lp = lptr(lpl)
  lp = lptr(lp)
  n3 = abs ( list(lp) )
  lpl = lend(n3)

  if ( list(lpl) <= 0 ) then
    ier = 3
    return
  end if
!
!  Delete N2 as a neighbor of N1, making N3 the first
!  neighbor, and test for error 4 (N2 not a neighbor
!  of N1).  Note that previously computed pointers may
!  no longer be valid following the call to DELNB.
!
  call delnb ( n1, n2, n, list, lptr, lend, lnew, lph )

  if ( lph < 0 ) then
    ier = 4
    return
  end if
!
!  Delete N1 as a neighbor of N2, making N3 the new last neighbor.
!
  call delnb ( n2, n1, n, list, lptr, lend, lnew, lph )
!
!  Make N3 a boundary node with first neighbor N2 and last neighbor N1.
!
  lp = lstptr ( lend(n3), n1, list, lptr )
  lend(n3) = lp
  list(lp) = -n1
!
!  No errors encountered.
!
  ier = 0

  return

end subroutine delarc




subroutine delnb ( n0, nb, n, list, lptr, lend, lnew, lph )

!*****************************************************************************80
!
!! DELNB deletes a neighbor from an adjacency list.
!
!  Discussion:
!
!    This subroutine deletes a neighbor NB from the adjacency
!    list of node N0 (but N0 is not deleted from the adjacency
!    list of NB) and, if NB is a boundary node, makes N0 a
!    boundary node.  For pointer (LIST index) LPH to NB as a
!    neighbor of N0, the empty LIST,LPTR location LPH is filled
!    in with the values at LNEW-1, pointer LNEW-1 (in LPTR and
!    possibly in LEND) is changed to LPH, and LNEW is decremented
!    This requires a search of LEND and LPTR entailing an
!    expected operation count of O(N).
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
!    Input, integer ( kind = 4 ) N0, NB, indexes, in the range 1 to N, of a 
!    pair of nodes such that NB is a neighbor of N0.
!    (N0 need not be a neighbor of NB.)
!
!    Input, integer ( kind = 4 ) N, the number of nodes in the triangulation.  
!    3 <= N.
!
!    Input/output, integer ( kind = 4 ) LIST(*), LPTR(*), LEND(N), LNEW, the 
!    data structure defining the triangulation.  On output, updated with
!    the removal of NB from the adjacency list of N0 unless LPH < 0.
!
!    Output, integer ( kind = 4 ) LPH, list pointer to the hole (NB as a 
!    neighbor of N0) filled in by the values at LNEW-1 or error indicator:
!    >0, if no errors were encountered.
!    -1, if N0, NB, or N is outside its valid range.
!    -2, if NB is not a neighbor of N0.
!
!  Local parameters:
!
!    I =   DO-loop index
!    LNW = LNEW-1 (output value of LNEW)
!    LP =  LIST pointer of the last neighbor of NB
!    LPB = Pointer to NB as a neighbor of N0
!    LPL = Pointer to the last neighbor of N0
!    LPP = Pointer to the neighbor of N0 that precedes NB
!    NN =  Local copy of N
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) lend(n)
  integer ( kind = 4 ) list(*)
  integer ( kind = 4 ) lnew
  integer ( kind = 4 ) lnw
  integer ( kind = 4 ) lp
  integer ( kind = 4 ) lpb
  integer ( kind = 4 ) lph
  integer ( kind = 4 ) lpl
  integer ( kind = 4 ) lpp
  integer ( kind = 4 ) lptr(*)
  integer ( kind = 4 ) n0
  integer ( kind = 4 ) nb
  integer ( kind = 4 ) nn

  nn = n
!
!  Test for error 1.
!
  if ( n0 < 1  .or. nn < n0 .or.  nb < 1  .or. &
      nn < nb .or. nn < 3 ) then
    lph = -1
    return
  end if
!
!  Find pointers to neighbors of N0:
!
!  LPL points to the last neighbor,
!  LPP points to the neighbor NP preceding NB, and
!  LPB points to NB.
!
  lpl = lend(n0)
  lpp = lpl
  lpb = lptr(lpp)

  do

    if ( list(lpb) == nb ) then
      go to 2
    end if

    lpp = lpb
    lpb = lptr(lpp)

    if ( lpb == lpl ) then
      exit
    end if

  end do
!
!  Test for error 2 (NB not found).
!
  if ( abs ( list(lpb) ) /= nb ) then
    lph = -2
    return
  end if
!
!  NB is the last neighbor of N0.  Make NP the new last
!  neighbor and, if NB is a boundary node, then make N0
!  a boundary node.
!
  lend(n0) = lpp
  lp = lend(nb)
  if ( list(lp) < 0 ) then
    list(lpp) = -list(lpp)
  end if
  go to 3
!
!  NB is not the last neighbor of N0.  If NB is a boundary
!  node and N0 is not, then make N0 a boundary node with
!  last neighbor NP.
!
2 continue

  lp = lend(nb)

  if ( list(lp) < 0  .and.  0 < list(lpl) ) then
    lend(n0) = lpp
    list(lpp) = -list(lpp)
  end if
!
!  Update LPTR so that the neighbor following NB now fol-
!  lows NP, and fill in the hole at location LPB.
!
3 continue

  lptr(lpp) = lptr(lpb)
  lnw = lnew - 1
  list(lpb) = list(lnw)
  lptr(lpb) = lptr(lnw)

  do i = nn, 1, -1
    if (lend(i) == lnw) then
      lend(i) = lpb
      exit
    end if
  end do

  do i = 1, lnw-1
    if (lptr(i) == lnw) then
      lptr(i) = lpb
    end if
  end do
!
!  No errors encountered.
!
  lnew = lnw
  lph = lpb

  return

end subroutine delnb




subroutine delnod ( k, ncc, lcc, n, x, y, list, lptr, lend, lnew, lwk, iwk, &
  ier )

!*****************************************************************************80
!
!! DELNOD deletes a node from a triangulation.
!
!  Discussion:
!
!    This subroutine deletes node K (along with all arcs
!    incident on node K) from a triangulation of N nodes in the
!    plane, and inserts arcs as necessary to produce a triangu-
!    lation of the remaining N-1 nodes.  If a Delaunay triangu-
!    lation is input, a Delaunay triangulation will result, and
!    thus, DELNOD reverses the effect of a call to subroutine
!    ADDNOD.
!
!    Note that a constraint node cannot be deleted by this
!    routine.  In order to delete a constraint node, it is
!    necessary to call this routine with NCC = 0, decrement the
!    appropriate LCC entries (LCC(I) such that K < LCC(I) ), and
!    then create (or restore) the constraints by a call to sub-
!    routine ADDCST.
!
!    Note that the deletion may result in all remaining nodes
!    being collinear.  This situation is not flagged.
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
!    Input, integer ( kind = 4 ) K, the index (for X and Y) of the node to be
!    deleted.  1 <= K < LCC(1).  (K <= N if NCC=0).
!
!    Input, integer ( kind = 4 ) NCC, the number of constraint curves.  
!    0 <= NCC.
!
!    Input/output, integer ( kind = 4 ) LCC(*), list of constraint curve 
!    starting indexes (or dummy array of length 1 if NCC = 0).  Refer to 
!    subroutine ADDCST.  On output, decremented by 1 to reflect the deletion 
!    of K unless NCC = 0 or 1 <= IER <= 4.
!
!    Input/output, integer ( kind = 4 ) N, the number of nodes in the 
!    triangulation.  4 <= N.  Note that N will be decremented following the 
!    deletion.
!
!    Input/output, real ( kind = 8 ) X(N), Y(N), the coordinates of the nodes
!    with non-constraint nodes in the first LCC(1)-1 locations if 0 < NCC.
!    On output, updated arrays of length N-1 containing nodal coordinates 
!    (with elements K+1,...,N shifted a position and thus overwriting 
!    element K) unless 1 <= IER <= 4.  (N here denotes the input value.)
!
!    Input/output, integer ( kind = 4 ) LIST(*), LPTR(*), LEND(N), LNEW, the 
!    data structure defining the triangulation.  Refer to subroutine TRMESH.
!    On output, updated to reflect the deletion unless IER /= 0.  Note
!    that the data structure may have been altered if 3 <= IER.
!
!    Input/output, integer ( kind = 4 ) LWK.  On input, the number of columns 
!    reserved for IWK.  LWK must be at least NNB-3, where NNB is the number of
!    neighbors of node K, including an extra pseudo-node if K is a 
!    boundary node.  On output, the number of IWK columns required unless 
!    IER = 1 or IER = 3.
!
!    Output, integer ( kind = 4 ) IWK(2,LWK), indexes of the endpoints of the 
!    new arcs added unless LWK = 0 or 1 <= IER <= 4.  (Arcs are associated with
!    columns, or pairs of adjacent elements if IWK is declared as a 
!    singly-subscripted array.)
!
!    Output, integer ( kind = 4 ) IER = Error indicator:
!    0, if no errors were encountered.
!    1, if K, NCC, N, or an LCC entry is outside its valid range or LWK < 0 
!      on input.
!    2, if more space is required in IWK.  Refer to LWK.
!    3, if the triangulation data structure is invalid on input.
!    4, if K is an interior node with 4 or more neighbors, and the number of
!      neighbors could not be reduced to 3 by swaps.  This could be caused by
!      floating point errors with collinear nodes or by an invalid data
!      structure.
!    5, if an error flag was returned by OPTIM.  An error message is written
!      to the standard output unit in this event.
!
  implicit none

  logical bdry
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ier
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) iwk(2,*)
  integer ( kind = 4 ) iwl
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) lcc(*)
  integer ( kind = 4 ) lccip1
  logical left
  integer ( kind = 4 ) lend(*)
  integer ( kind = 4 ) list(*)
  integer ( kind = 4 ) lnew
  integer ( kind = 4 ) lnw
  integer ( kind = 4 ) lp
  integer ( kind = 4 ) lp21
  integer ( kind = 4 ) lpf
  integer ( kind = 4 ) lph
  integer ( kind = 4 ) lpl
  integer ( kind = 4 ) lpl2
  integer ( kind = 4 ) lpn
  integer ( kind = 4 ) lptr(*)
  integer ( kind = 4 ) lstptr
  integer ( kind = 4 ) lwk
  integer ( kind = 4 ) lwkl
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) nbcnt
  integer ( kind = 4 ) ncc
  integer ( kind = 4 ) nfrst
  integer ( kind = 4 ) nit
  integer ( kind = 4 ) nl
  integer ( kind = 4 ) nn
  integer ( kind = 4 ) nnb
  integer ( kind = 4 ) nr
  real ( kind = 8 ) x(*)
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) xl
  real ( kind = 8 ) xr
  real ( kind = 8 ) y(*)
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2
  real ( kind = 8 ) yl
  real ( kind = 8 ) yr
!
!  Set N1 to K and NNB to the number of neighbors of N1 (plus
!  one if N1 is a boundary node), and test for errors.  LPF
!  and LPL are LIST indexes of the first and last neighbors
!  of N1, IWL is the number of IWK columns containing arcs,
!  and BDRY is TRUE iff N1 is a boundary node.
!
  ier = 0
  n1 = k
  nn = n

  if ( ncc < 0  .or.  n1 < 1  .or.  nn < 4  .or. lwk < 0 ) then
    ier = 1
    return
  end if

  lccip1 = nn + 1

  do i = ncc, 1, -1

    if ( lccip1 - lcc(i) < 3 ) then
      ier = 1
      return
    end if

    lccip1 = lcc(i)

  end do

  if ( lccip1 <= n1 ) then
    ier = 1
    return
  end if

  lpl = lend(n1)
  lpf = lptr(lpl)
  nnb = nbcnt ( lpl, lptr )
  bdry = list(lpl) < 0

  if ( bdry ) then
    nnb = nnb + 1
  end if

  if ( nnb < 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DELNOD - Fatal error!'
    write ( *, '(a)' ) '  NNB < 3.'
    write ( *, '(a,i6)' ) '  NNB = ', nnb
    ier = 3
    return
  end if

  lwkl = lwk
  lwk = nnb - 3

  if ( lwkl < lwk ) then
    ier = 2
    return
  end if

  iwl = 0

  if ( nnb == 3 ) then
    go to 5
  end if
!
!  Initialize for loop on arcs N1-N2 for neighbors N2 of N1,
!  beginning with the second neighbor.  NR and NL are the
!  neighbors preceding and following N2, respectively, and
!  LP indexes NL.  The loop is exited when all possible
!  swaps have been applied to arcs incident on N1.  If N1
!  is interior, the number of neighbors will be reduced
!  to 3.
!
  x1 = x(n1)
  y1 = y(n1)
  nfrst = list(lpf)
  nr = nfrst
  xr = x(nr)
  yr = y(nr)
  lp = lptr(lpf)
  n2 = list(lp)
  x2 = x(n2)
  y2 = y(n2)
  lp = lptr(lp)
!
!  Top of loop:  set NL to the neighbor following N2.
!
2 continue

  nl = abs ( list(lp) )

  if ( nl == nfrst .and. bdry ) then
    go to 5
  end if

  xl = x(nl)
  yl = y(nl)
!
!  Test for a convex quadrilateral.  To avoid an incorrect
!  test caused by collinearity, use the fact that if N1
!  is a boundary node, then N1 LEFT NR->NL and if N2 is
!  a boundary node, then N2 LEFT NL->NR.
!
  lpl2 = lend(n2)

  if ( (bdry  .or.  left(xr,yr,xl,yl,x1,y1))  .and. &
       (list(lpl2) < 0  .or. &
        left(xl,yl,xr,yr,x2,y2)) ) then
    go to 3
  end if
!
!  Nonconvex quadrilateral -- no swap is possible.
!
  nr = n2
  xr = x2
  yr = y2
  go to 4
!
!  The quadrilateral defined by adjacent triangles
!  (N1,N2,NL) and (N2,N1,NR) is convex.  Swap in
!  NL-NR and store it in IWK.  Indexes larger than N1
!  must be decremented since N1 will be deleted from
!  X and Y.
!
3 continue

  call swap ( nl, nr, n1, n2, list, lptr, lend, lp21 )
  iwl = iwl + 1

  if ( nl <= n1 ) then
    iwk(1,iwl) = nl
  else
    iwk(1,iwl) = nl - 1
  end if

  if ( nr <= n1 ) then
    iwk(2,iwl) = nr
  else
    iwk(2,iwl) = nr - 1
  end if
!
!  Recompute the LIST indexes LPL,LP and decrement NNB.
!
  lpl = lend(n1)
  nnb = nnb - 1

  if ( nnb == 3 ) then
    go to 5
  end if

  lp = lstptr ( lpl, nl, list, lptr )

  if ( nr == nfrst ) then
    go to 4
  end if
!
!  NR is not the first neighbor of N1.
!  Back up and test N1-NR for a swap again:  Set N2 to
!  NR and NR to the previous neighbor of N1 -- the
!  neighbor of NR which follows N1.  LP21 points to NL
!  as a neighbor of NR.
!
  n2 = nr
  x2 = xr
  y2 = yr
  lp21 = lptr(lp21)
  lp21 = lptr(lp21)
  nr = abs ( list(lp21) )
  xr = x(nr)
  yr = y(nr)
  go to 2
!
!  Bottom of loop -- test for invalid termination.
!
4 continue

  if ( n2 == nfrst ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DELNOD - Fatal error!'
    write ( *, '(a)' ) '  N2 = NFRST.'
    ier = 4
    return
  end if

  n2 = nl
  x2 = xl
  y2 = yl
  lp = lptr(lp)
  go to 2
!
!  Delete N1 from the adjacency list of N2 for all neighbors
!  N2 of N1.  LPL points to the last neighbor of N1.
!  LNEW is stored in local variable LNW.
!
5 continue

  lp = lpl
  lnw = lnew
!
!  Loop on neighbors N2 of N1, beginning with the first.
!
6 continue

    lp = lptr(lp)
    n2 = abs ( list(lp) )

    call delnb ( n2, n1, n, list, lptr, lend, lnw, lph )

    if ( lph < 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DELNOD - Fatal error!'
      write ( *, '(a)' ) '  LPH < 0.'
      ier = 3
      return
    end if
!
!  LP and LPL may require alteration.
!
    if ( lpl == lnw ) then
      lpl = lph
    end if

    if ( lp == lnw ) then
      lp = lph
    end if

    if ( lp /= lpl ) then
      go to 6
    end if
!
!  Delete N1 from X, Y, and LEND, and remove its adjacency
!  list from LIST and LPTR.  LIST entries (nodal indexes)
!  which are larger than N1 must be decremented.
!
  nn = nn - 1

  if ( nn < n1 ) then
    go to 9
  end if

  do i = n1,nn
    x(i) = x(i+1)
    y(i) = y(i+1)
    lend(i) = lend(i+1)
  end do

  do i = 1,lnw-1

    if (list(i) > n1) then
      list(i) = list(i) - 1
    end if

    if (list(i) < -n1) then
      list(i) = list(i) + 1
    end if

  end do
!
!  For LPN = first to last neighbors of N1, delete the
!  preceding neighbor (indexed by LP).
!
!  Each empty LIST,LPTR location LP is filled in with the
!  values at LNW-1, and LNW is decremented.  All pointers
!  (including those in LPTR and LEND) with value LNW-1
!  must be changed to LP.
!
!  LPL points to the last neighbor of N1.
!
9 continue

  if ( bdry ) then
    nnb = nnb - 1
  end if

  lpn = lpl

  do j = 1, nnb

    lnw = lnw - 1
    lp = lpn
    lpn = lptr(lp)
    list(lp) = list(lnw)
    lptr(lp) = lptr(lnw)

    if ( lptr(lpn) == lnw ) then
      lptr(lpn) = lp
    end if

    if ( lpn == lnw ) then
      lpn = lp
    end if

    do i = nn,1,-1
      if (lend(i) == lnw) then
        lend(i) = lp
        exit
      end if
    end do

    do i = lnw-1, 1, -1
      if ( lptr(i) == lnw ) then
        lptr(i) = lp
      end if
    end do

  end do
!
!  Decrement LCC entries.
!
  lcc(1:ncc) = lcc(1:ncc) - 1
!
!  Update N and LNEW, and optimize the patch of triangles
!  containing K (on input) by applying swaps to the arcs in IWK.
!
  n = nn
  lnew = lnw

  if ( 0 < iwl ) then
    nit = 4 * iwl
    call optim ( x, y, iwl, list, lptr, lend, nit, iwk, ierr )
    if ( ierr /= 0 ) then
      ier = 5
      write (*,100) nit, ierr
  100 format (//5x,'*** error in optim:  nit = ',i4, &
          ', ier = ',i1,' ***'/)
      return
    end if
  end if
!
!  Successful termination.
!
  ier = 0

  return

end subroutine delnod
