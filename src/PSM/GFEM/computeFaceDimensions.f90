subroutine addGroupLengths(ni, nj, nsurf, nedge, ngroup, surfs, &
     surf_edge, edge_group, surfEdgeLengths, groupLengths0, groupCount0, groupLengths, groupCount)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) ni, nj, nsurf, nedge, ngroup, surfs, surf_edge, edge_group, surfEdgeLengths, groupLengths0, groupCount0
  !f2py intent(out) groupLengths, groupCount
  !f2py depend(ni,nj) surfs
  !f2py depend(nsurf) surf_edge, surfEdgeLengths
  !f2py depend(nedge) edge_group
  !f2py depend(ngroup) groupLengths, groupCount, groupLengths0, groupCount0

  !Input
  integer, intent(in) ::  ni, nj, nsurf, nedge, ngroup, surfs(ni,nj)
  integer, intent(in) ::  surf_edge(nsurf,2,2), edge_group(nedge)
  double precision, intent(in) ::  surfEdgeLengths(nsurf,2,2)
  double precision, intent(in) ::  groupLengths0(ngroup)
  integer, intent(in) ::  groupCount0(ngroup)

  !Output
  double precision, intent(out) ::  groupLengths(ngroup)
  integer, intent(out) ::  groupCount(ngroup)

  !Working
  integer s, i, j, k, ugroup, vgroup


  groupLengths = groupLengths0
  groupCount = groupCount0
  do i=1,ni
     do j=1,nj
        s = surfs(i,j)
        if (s .gt. 0) then
           ugroup = edge_group(abs(surf_edge(s,1,1)))
           vgroup = edge_group(abs(surf_edge(s,2,1)))
           do k=1,2
              groupLengths(ugroup) = groupLengths(ugroup) + surfEdgeLengths(s,1,k)
              groupLengths(vgroup) = groupLengths(vgroup) + surfEdgeLengths(s,2,k)
           end do
           groupCount(ugroup) = groupCount(ugroup) + 2
           groupCount(vgroup) = groupCount(vgroup) + 2
        end if
     end do
  end do

end subroutine addGroupLengths




subroutine computeFaceDimensions(ni, nj, nsurf, nedge, ngroup, &
     surfs, surf_edge, edge_group, groupLengths, idims, jdims)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) ni, nj, nsurf, nedge, ngroup, surfs, surf_edge, edge_group, groupLengths
  !f2py intent(out) idims, jdims
  !f2py depend(ni,nj) surfs
  !f2py depend(nsurf) surf_edge
  !f2py depend(nedge) edge_group
  !f2py depend(ngroup) groupLengths
  !f2py depend(ni) idims
  !f2py depend(nj) jdims

  !Input
  integer, intent(in) ::  ni, nj, nsurf, nedge, ngroup
  integer, intent(in) ::  surfs(ni,nj), surf_edge(nsurf,2,2), edge_group(nedge)
  double precision, intent(in) ::  groupLengths(ngroup)

  !Output
  double precision, intent(out) ::  idims(ni+1), jdims(nj+1)

  !Working
  integer s, i, j, ugroup, vgroup

  idims(1) = 0.
  do i=1,ni
     s = surfs(i,1)
     if (s .gt. 0) then
        ugroup = edge_group(abs(surf_edge(s,1,1)))
        idims(i+1) = idims(i) + groupLengths(ugroup)
     else
        idims(i+1) = idims(i)
     end if
  end do

  jdims(1) = 0.
  do j=1,nj
     s = surfs(1,j)
     if (s .gt. 0) then
        vgroup = edge_group(abs(surf_edge(s,2,1)))
        jdims(j+1) = jdims(j) + groupLengths(vgroup)
     else
        jdims(j+1) = jdims(j)
     end if
  end do

  idims(:) = idims(:)/idims(ni+1)
  jdims(:) = jdims(:)/jdims(nj+1)

  idims(1) = 0.
  jdims(1) = 0.
  idims(ni+1) = 1.
  jdims(nj+1) = 1.

end subroutine computeFaceDimensions




subroutine computeFaceDimensions0(ni, nj, nsurf, &
     surfs, surfEdgeLengths, idims, jdims)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) ni, nj, nsurf, surfs, surfEdgeLengths
  !f2py intent(out) idims, jdims
  !f2py depend(ni,nj) surfs
  !f2py depend(nsurf) surfEdgeLengths
  !f2py depend(ni) idims
  !f2py depend(nj) jdims

  !Input
  integer, intent(in) ::  ni, nj, nsurf, surfs(ni,nj)
  double precision, intent(in) ::  surfEdgeLengths(nsurf,2,2)

  !Output
  double precision, intent(out) ::  idims(ni+1), jdims(nj+1)

  !Working
  double precision idims1(ni+1), idims2(ni+1)
  double precision jdims1(nj+1), jdims2(nj+1)
  integer i, j

  idims1(1) = 0.
  idims2(1) = 0.
  do i=1,ni
     idims1(i+1) = idims1(i) + surfEdgeLengths(surfs(i,1),1,1)
     idims2(i+1) = idims2(i) + surfEdgeLengths(surfs(i,nj),1,2)
  end do
  idims1(:) = idims1(:)/idims1(ni+1)
  idims2(:) = idims2(:)/idims2(ni+1)

  jdims1(1) = 0.
  jdims2(1) = 0.
  do j=1,nj
     jdims1(j+1) = jdims1(j) + surfEdgeLengths(surfs(1,j),2,1)
     jdims2(j+1) = jdims2(j) + surfEdgeLengths(surfs(ni,j),2,2)
  end do
  jdims1(:) = jdims1(:)/jdims1(nj+1)
  jdims2(:) = jdims2(:)/jdims2(nj+1)

  idims = 0.5*(idims1 + idims2)
  jdims = 0.5*(jdims1 + jdims2)

  idims(1) = 0.
  jdims(1) = 0.
  idims(ni+1) = 1.
  jdims(nj+1) = 1.

end subroutine computeFaceDimensions0
