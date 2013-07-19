subroutine computeIntersectionVerts(nvert0, nedge, ngroup, nint, nvert, &
     verts0, edges, edge_group, groupIntPtr, groupInts, verts)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nvert0, nedge, ngroup, nint, nvert, verts0, edges, edge_group, groupIntPtr, groupInts
  !f2py intent(out) verts
  !f2py depend(nvert0) verts0
  !f2py depend(nedge) edges, edge_group
  !f2py depend(ngroup) groupIntPtr
  !f2py depend(nint) groupInts
  !f2py depend(nvert) verts

  !Input
  integer, intent(in) ::  nvert0, nedge, ngroup, nint, nvert
  double precision, intent(in) ::  verts0(nvert0,2)
  integer, intent(in) ::  edges(nedge,2), edge_group(nedge), groupIntPtr(ngroup,2)
  double precision, intent(in) ::  groupInts(nint)

  !Output
  double precision, intent(out) ::  verts(nvert,2)

  !Working
  double precision t, A(2), B(2)
  integer ivert, iedge, igroup, int

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
  end do

end subroutine computeIntersectionVerts




subroutine countIntersectionVerts(nedge, ngroup, edge_group, groupIntPtr, nvert)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nedge, ngroup, edge_group, groupIntPtr
  !f2py intent(out) nvert
  !f2py depend(nedge) edge_group
  !f2py depend(ngroup) groupIntPtr

  !Input
  integer, intent(in) ::  nedge, ngroup, edge_group(nedge), groupIntPtr(ngroup,2)

  !Output
  integer, intent(out) ::  nvert

  !Working
  integer iedge, igroup

  nvert = 0
  do iedge=1,nedge
     igroup = abs(edge_group(iedge))
     nvert = nvert + groupIntPtr(igroup,2) - groupIntPtr(igroup,1) + 1
  end do

end subroutine countIntersectionVerts
