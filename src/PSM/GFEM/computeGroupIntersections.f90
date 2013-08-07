subroutine countGroupIntersections(nvert, nedge, ngroup, verts, edges, &
     edge_group, groupIntCount0, groupIntCount)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nvert, nedge, ngroup, verts, edges, edge_group, groupIntCount0
  !f2py intent(out) groupIntCount
  !f2py depend(nvert) verts
  !f2py depend(nedge) edges, edge_group
  !f2py depend(ngroup) groupIntCount0, groupIntCount

  !Input
  integer, intent(in) ::  nvert, nedge, ngroup
  double precision, intent(in) ::  verts(nvert,2)
  integer, intent(in) ::  edges(nedge,2), edge_group(nedge), groupIntCount0(ngroup)

  !Output
  integer, intent(out) ::  groupIntCount(ngroup)

  !Working
  integer iedge, ivert
  logical validSplit

  groupIntCount = groupIntCount0

  do iedge=1,nedge
     do ivert=1,nvert
        if (validSplit(verts(edges(iedge,1),:), verts(edges(iedge,2),:), verts(ivert,:))) then
           groupIntCount(abs(edge_group(iedge))) = groupIntCount(abs(edge_group(iedge))) + 1
        end if
     end do
  end do

end subroutine countGroupIntersections




subroutine computeGroupIntPtr(ngroup, groupIntCount, groupIntPtr)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) ngroup, groupIntCount
  !f2py intent(out) groupIntPtr
  !f2py depend(ngroup) groupIntCount, groupIntPtr

  !Input
  integer, intent(in) ::  ngroup, groupIntCount(ngroup)

  !Output
  integer, intent(out) ::  groupIntPtr(ngroup,2)

  !Working
  integer igroup

  groupIntPtr(1,1) = 1
  groupIntPtr(1,2) = groupIntCount(1)
  do igroup=2,ngroup
     groupIntPtr(igroup,1) = groupIntPtr(igroup-1,2) + 1
     groupIntPtr(igroup,2) = groupIntPtr(igroup-1,2) + groupIntCount(igroup)
  end do

end subroutine computeGroupIntPtr





subroutine computeGroupIntersections(nvert, nedge, ngroup, nint, verts, edges, &
     edge_group, groupIntPtr, groupInts0, groupInts)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nvert, nedge, ngroup, nint, verts, edges, edge_group, groupIntPtr, groupInts0
  !f2py intent(out) groupInts
  !f2py depend(nvert) verts
  !f2py depend(nedge) edges, edge_group
  !f2py depend(ngroup) groupIntPtr
  !f2py depend(nint) groupInts0, groupInts

  !Input
  integer, intent(in) ::  nvert, nedge, ngroup, nint
  double precision, intent(in) ::  verts(nvert,2)
  integer, intent(in) ::  edges(nedge,2), edge_group(nedge), groupIntPtr(ngroup,2)
  double precision, intent(in) ::  groupInts0(nint)

  !Output
  double precision, intent(out) ::  groupInts(nint)

  !Working
  integer iedge, ivert, igroup, index, counter(ngroup)
  logical validSplit
  double precision split

  groupInts = groupInts0

  counter(:) = 0
  do iedge=1,nedge
     do ivert=1,nvert
        if (validSplit(verts(edges(iedge,1),:), verts(edges(iedge,2),:), verts(ivert,:))) then
           igroup = abs(edge_group(iedge))
           index = groupIntPtr(igroup,1) + counter(igroup)
           if (edge_group(iedge) .gt. 0) then
              groupInts(index) = split(verts(edges(iedge,1),:), verts(edges(iedge,2),:), verts(ivert,:))
           else
              groupInts(index) = 1 - split(verts(edges(iedge,1),:), verts(edges(iedge,2),:), verts(ivert,:))
           end if
           counter(igroup) = counter(igroup) + 1
        end if
     end do
  end do

end subroutine computeGroupIntersections




function validSplit(A, B, C)

  implicit none

  !Input
  double precision, intent(in) ::  A(2), B(2), C(2)

  !Output
  logical validSplit

  !Working
  double precision v1(2), v2(2), det

  v1 = C - A
  v2 = C - B
  det = abs(v1(1)*v2(2) - v1(2)*v2(1))

  if ((det .lt. 1e-10) .and. (dot_product(v1,v2) .lt. 0)) then
     validSplit = .True.
  else
     validSplit = .False.
  end if

end function validSplit




function split(A, B, C)

  implicit none

  !Input
  double precision, intent(in) ::  A(2), B(2), C(2)

  !Output
  double precision split

  !Working
  double precision d, d2

  d = sqrt(dot_product(C - A, C - A))
  d2 = sqrt(dot_product(B - A, B - A))

  split = d/d2

end function split
