subroutine computeTri2Quad(nvert0, nedge0, ntri, nvert, nedge, &
     verts0, edges0, triangles, verts, edges)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nvert0, nedge0, ntri, nvert, nedge, verts0, edges0, triangles
  !f2py intent(out) verts, edges
  !f2py depend(nvert0) verts0
  !f2py depend(nedge0) edges0
  !f2py depend(ntri) triangles
  !f2py depend(nvert) verts
  !f2py depend(nedge) edges

  !Input
  integer, intent(in) ::  nvert0, nedge0, ntri, nvert, nedge
  double precision, intent(in) ::  verts0(nvert0,2)
  integer, intent(in) ::  edges0(nedge0,2), triangles(ntri,3)

  !Output
  double precision, intent(out) ::  verts(nvert,2)
  integer, intent(out) ::  edges(nedge,2)

  !Working
  integer itri, ivert, iedge
  integer i1, i2, i3

  verts(1:nvert0,:) = verts0(:,:)
  edges(1:nedge0,:) = edges0(:,:)

  ivert = nvert0 + 1
  iedge = nedge0 + 0
  do itri=1,ntri
     i1 = triangles(itri,1)
     i2 = triangles(itri,2)
     i3 = triangles(itri,3)
     verts(ivert,:) = (verts(i1,:) + verts(i2,:) + verts(i3,:))/3.0
     verts(ivert+1,:) = (verts(i2,:) + verts(i3,:))/2.0
     verts(ivert+2,:) = (verts(i3,:) + verts(i1,:))/2.0
     verts(ivert+3,:) = (verts(i1,:) + verts(i2,:))/2.0
     edges(iedge+1,1) = ivert+1
     edges(iedge+2,1) = ivert+2
     edges(iedge+3,1) = ivert+3
     edges(iedge+1:iedge+3,2) = ivert
     ivert = ivert + 4
     iedge = iedge + 3
  end do

end subroutine computeTri2Quad
