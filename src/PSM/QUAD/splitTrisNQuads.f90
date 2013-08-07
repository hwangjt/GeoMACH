subroutine splitTris(nvert0, nedge0, ntri, nvert, nedge, &
     verts0, edges0, edgeCon0, triangles, verts, edges, edgeCon)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nvert0, nedge0, ntri, nvert, nedge, verts0, edges0, edgeCon0, triangles
  !f2py intent(out) verts, edges, edgeCon
  !f2py depend(nvert0) verts0
  !f2py depend(nedge0) edges0, edgeCon0
  !f2py depend(ntri) triangles
  !f2py depend(nvert) verts
  !f2py depend(nedge) edges, edgeCon

  !Input
  integer, intent(in) ::  nvert0, nedge0, ntri, nvert, nedge
  double precision, intent(in) ::  verts0(nvert0,2)
  integer, intent(in) ::  edges0(nedge0,2), triangles(ntri,3)
  logical, intent(in) ::  edgeCon0(nedge0)

  !Output
  double precision, intent(out) ::  verts(nvert,2)
  integer, intent(out) ::  edges(nedge,2)
  logical, intent(out) ::  edgeCon(nedge)

  !Working
  integer ivert, iedge, itri
  integer i1, i2, i3

  edgeCon(:) = .False.

  verts(1:nvert0,:) = verts0(:,:)
  edges(1:nedge0,:) = edges0(:,:)
  edgeCon(1:nedge0) = edgeCon0

  ivert = nvert0
  iedge = nedge0
  do itri=1,ntri
     i1 = triangles(itri,1)
     i2 = triangles(itri,2)
     i3 = triangles(itri,3)
     verts(ivert+4,:) = (verts(i1,:) + verts(i2,:) + verts(i3,:))/3.0
     verts(ivert+1,:) = (verts(i2,:) + verts(i3,:))/2.0
     verts(ivert+2,:) = (verts(i3,:) + verts(i1,:))/2.0
     verts(ivert+3,:) = (verts(i1,:) + verts(i2,:))/2.0
     edges(iedge+1,1) = ivert+1
     edges(iedge+2,1) = ivert+2
     edges(iedge+3,1) = ivert+3
     edges(iedge+1:iedge+3,2) = ivert+4
     ivert = ivert + 4
     iedge = iedge + 3
  end do
  if (ivert .ne. nvert) then
     print *, 'Error in splitTris', ivert,nvert
     call exit(1)
  end if
  if (iedge .ne. nedge) then
     print *, 'Error in splitTris', iedge,nedge
     call exit(1)
  end if

end subroutine splitTris




subroutine splitQuads(nvert0, nedge0, nquad, nvert, nedge, &
     verts0, edges0, edgeCon0, quads, verts, edges, edgeCon)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nvert0, nedge0, nquad, nvert, nedge, verts0, edges0, edgeCon0, quads
  !f2py intent(out) verts, edges, edgeCon
  !f2py depend(nvert0) verts0
  !f2py depend(nedge0) edges0, edgeCon0
  !f2py depend(nquad) quads
  !f2py depend(nvert) verts
  !f2py depend(nedge) edges, edgeCon

  !Input
  integer, intent(in) ::  nvert0, nedge0, nquad, nvert, nedge
  double precision, intent(in) ::  verts0(nvert0,2)
  integer, intent(in) ::  edges0(nedge0,2), quads(nquad,4)
  logical, intent(in) ::  edgeCon0(nedge0)

  !Output
  double precision, intent(out) ::  verts(nvert,2)
  integer, intent(out) ::  edges(nedge,2)
  logical, intent(out) ::  edgeCon(nedge)

  !Working
  integer ivert, iedge, iquad
  integer i1, i2, i3, i4

  edgeCon(:) = .False.

  verts(1:nvert0,:) = verts0(:,:)
  edges(1:nedge0,:) = edges0(:,:)
  edgeCon(1:nedge0) = edgeCon0

  ivert = nvert0
  iedge = nedge0
  do iquad=1,nquad
     i1 = quads(iquad,1)
     i2 = quads(iquad,2)
     i3 = quads(iquad,3)
     i4 = quads(iquad,4)
     verts(ivert+5,:) = (verts(i1,:) + verts(i2,:) + verts(i3,:) + verts(i4,:))/4.0
     verts(ivert+1,:) = (verts(i1,:) + verts(i2,:))/2.0
     verts(ivert+2,:) = (verts(i2,:) + verts(i3,:))/2.0
     verts(ivert+3,:) = (verts(i3,:) + verts(i4,:))/2.0
     verts(ivert+4,:) = (verts(i4,:) + verts(i1,:))/2.0
     edges(iedge+1,1) = ivert+1
     edges(iedge+2,1) = ivert+2
     edges(iedge+3,1) = ivert+3
     edges(iedge+4,1) = ivert+4
     edges(iedge+1:iedge+4,2) = ivert+5
     ivert = ivert + 5
     iedge = iedge + 4
  end do
  if (ivert .ne. nvert) then
     print *, 'Error in splitQuads', ivert,nvert
     call exit(1)
  end if
  if (iedge .ne. nedge) then
     print *, 'Error in splitQuads', iedge,nedge
     call exit(1)
  end if

end subroutine splitQuads




subroutine splitTrisNQuads(nvert0, nedge0, ntri, nquad, nvert, nedge, &
     verts0, edges0, edgeCon0, triangles, quads, verts, edges, edgeCon)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nvert0, nedge0, ntri, nquad, nvert, nedge, verts0, edges0, edgeCon0, triangles, quads
  !f2py intent(out) verts, edges, edgeCon
  !f2py depend(nvert0) verts0
  !f2py depend(nedge0) edges0, edgeCon0
  !f2py depend(ntri) triangles
  !f2py depend(nquad) quads
  !f2py depend(nvert) verts
  !f2py depend(nedge) edges, edgeCon

  !Input
  integer, intent(in) ::  nvert0, nedge0, ntri, nquad, nvert, nedge
  double precision, intent(in) ::  verts0(nvert0,2)
  integer, intent(in) ::  edges0(nedge0,2), triangles(ntri,3), quads(nquad,4)
  logical, intent(in) ::  edgeCon0(nedge0)

  !Output
  double precision, intent(out) ::  verts(nvert,2)
  integer, intent(out) ::  edges(nedge,2)
  logical, intent(out) ::  edgeCon(nedge)

  !Working
  integer ivert, iedge, itri, iquad
  integer i1, i2, i3, i4

  edgeCon(:) = .False.

  verts(1:nvert0,:) = verts0(:,:)
  edges(1:nedge0,:) = edges0(:,:)
  edgeCon(1:nedge0) = edgeCon0

  ivert = nvert0
  iedge = nedge0
  do itri=1,ntri
     i1 = triangles(itri,1)
     i2 = triangles(itri,2)
     i3 = triangles(itri,3)
     verts(ivert+4,:) = (verts(i1,:) + verts(i2,:) + verts(i3,:))/3.0
     verts(ivert+1,:) = (verts(i2,:) + verts(i3,:))/2.0
     verts(ivert+2,:) = (verts(i3,:) + verts(i1,:))/2.0
     verts(ivert+3,:) = (verts(i1,:) + verts(i2,:))/2.0
     edges(iedge+1,1) = ivert+1
     edges(iedge+2,1) = ivert+2
     edges(iedge+3,1) = ivert+3
     edges(iedge+1:iedge+3,2) = ivert+4
     ivert = ivert + 4
     iedge = iedge + 3
  end do
  do iquad=1,nquad
     i1 = quads(iquad,1)
     i2 = quads(iquad,2)
     i3 = quads(iquad,3)
     i4 = quads(iquad,4)
     verts(ivert+5,:) = (verts(i1,:) + verts(i2,:) + verts(i3,:) + verts(i4,:))/4.0
     verts(ivert+1,:) = (verts(i1,:) + verts(i2,:))/2.0
     verts(ivert+2,:) = (verts(i2,:) + verts(i3,:))/2.0
     verts(ivert+3,:) = (verts(i3,:) + verts(i4,:))/2.0
     verts(ivert+4,:) = (verts(i4,:) + verts(i1,:))/2.0
     edges(iedge+1,1) = ivert+1
     edges(iedge+2,1) = ivert+2
     edges(iedge+3,1) = ivert+3
     edges(iedge+4,1) = ivert+4
     edges(iedge+1:iedge+4,2) = ivert+5
     ivert = ivert + 5
     iedge = iedge + 4
  end do
  if (ivert .ne. nvert) then
     print *, 'Error in splitTrisNQuads', ivert,nvert
     call exit(1)
  end if
  if (iedge .ne. nedge) then
     print *, 'Error in splitTrisNQuads', iedge,nedge
     call exit(1)
  end if

end subroutine splitTrisNQuads
