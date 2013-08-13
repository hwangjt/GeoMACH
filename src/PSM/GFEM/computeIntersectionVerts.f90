subroutine computeIntersectionVerts(nvert0, nedge, ngroup, nint, nsplit, nvert, &
     verts0, edges, edge_group, groupIntPtr, groupInts, groupSplitPtr, groupSplits, verts)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nvert0, nedge, ngroup, nint, nsplit, nvert, verts0, edges, edge_group, groupIntPtr, groupInts, groupSplitPtr, groupSplits
  !f2py intent(out) verts
  !f2py depend(nvert0) verts0
  !f2py depend(nedge) edges, edge_group
  !f2py depend(ngroup) groupIntPtr, groupSplitPtr
  !f2py depend(nint) groupInts
  !f2py depend(nsplit) groupSplits
  !f2py depend(nvert) verts

  !Input
  integer, intent(in) ::  nvert0, nedge, ngroup, nint, nsplit, nvert
  double precision, intent(in) ::  verts0(nvert0,2)
  integer, intent(in) ::  edges(nedge,2), edge_group(nedge)
  integer, intent(in) ::  groupIntPtr(ngroup,2), groupSplitPtr(ngroup,2)
  double precision, intent(in) ::  groupInts(nint)
  double precision, intent(in) ::  groupSplits(nsplit)

  !Output
  double precision, intent(out) ::  verts(nvert,2)

  !Working
  double precision t, A(2), B(2)
  integer ivert, iedge, igroup, int, isplit

  verts(1:nvert0,:) = verts0
  
  ivert = nvert0 + 1
  do iedge=1,nedge
     igroup = abs(edge_group(iedge))
     if (edge_group(iedge) .gt. 0) then
        A = verts0(edges(iedge,1),:)
        B = verts0(edges(iedge,2),:)
     else
        A = verts0(edges(iedge,2),:)
        B = verts0(edges(iedge,1),:)
     end if
     do int=groupIntPtr(igroup,1),groupIntPtr(igroup,2)
        t = groupInts(int)
        verts(ivert,:) = (1-t)*A + t*B
        ivert = ivert + 1
     end do
     do isplit=groupSplitPtr(igroup,1),groupSplitPtr(igroup,2)
        t = groupSplits(isplit)
        verts(ivert,:) = (1-t)*A + t*B
        ivert = ivert + 1
     end do
  end do

end subroutine computeIntersectionVerts




subroutine countIntersectionVerts(nedge, ngroup, edge_group, &
     groupIntPtr, groupSplitPtr, nvert)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nedge, ngroup, edge_group, groupIntPtr, groupSplitPtr
  !f2py intent(out) nvert
  !f2py depend(nedge) edge_group
  !f2py depend(ngroup) groupIntPtr, groupSplitPtr

  !Input
  integer, intent(in) ::  nedge, ngroup, edge_group(nedge)
  integer, intent(in) ::  groupIntPtr(ngroup,2), groupSplitPtr(ngroup,2)

  !Output
  integer, intent(out) ::  nvert

  !Working
  integer iedge, igroup

  nvert = 0
  do iedge=1,nedge
     igroup = abs(edge_group(iedge))
     nvert = nvert + groupIntPtr(igroup,2) - groupIntPtr(igroup,1) + 1
     nvert = nvert + groupSplitPtr(igroup,2) - groupSplitPtr(igroup,1) + 1
  end do

end subroutine countIntersectionVerts
