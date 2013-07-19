subroutine computeSurfaceEdges(nvert, nedge0, nedge, &
     xmin, xmax, ymin, ymax, verts, edges0, edges)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nvert, nedge0, nedge, xmin, xmax, ymin, ymax, verts, edges0
  !f2py intent(out) edges
  !f2py depend(nvert) verts
  !f2py depend(nedge0) edges0
  !f2py depend(nedge) edges

  !Input
  integer, intent(in) ::  nvert, nedge0, nedge
  double precision, intent(in) ::  xmin, xmax, ymin, ymax
  double precision, intent(in) ::  verts(nvert,2)
  integer, intent(in) ::  edges0(nedge0,2)

  !Output
  double precision, intent(out) ::  edges(nedge,2,2)

  !Working
  logical xTrue, yTrue
  integer iedge0, iedge, j
  double precision x, y

  iedge = 1
  do iedge0=1,nedge0
     x = 0.5*verts(edges0(iedge0,1),1) + 0.5*verts(edges0(iedge0,2),1)
     y = 0.5*verts(edges0(iedge0,1),2) + 0.5*verts(edges0(iedge0,2),2)
     xTrue = (xmin - 1e-10 .le. x) .and. (x .le. xmax + 1e-10)
     yTrue = (ymin - 1e-10 .le. y) .and. (y .le. ymax + 1e-10)
     if (xTrue .and. yTrue) then
        edges(iedge,1,:) = verts(edges0(iedge0,1),:)
        edges(iedge,2,:) = verts(edges0(iedge0,2),:)
        iedge = iedge + 1
     end if
  end do

  do iedge=1,nedge
     do j=1,2
        edges(iedge,j,1) = (edges(iedge,j,1)-xmin)/(xmax-xmin)
        edges(iedge,j,2) = (edges(iedge,j,2)-ymin)/(ymax-ymin)
     end do
  end do

end subroutine computeSurfaceEdges




subroutine countSurfaceEdges(nvert, nedge0, &
     xmin, xmax, ymin, ymax, verts, edges, nedge)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nvert, nedge0, xmin, xmax, ymin, ymax, verts, edges
  !f2py intent(out) nedge
  !f2py depend(nvert) verts
  !f2py depend(nedge0) edges

  !Input
  integer, intent(in) ::  nvert, nedge0
  double precision, intent(in) ::  xmin, xmax, ymin, ymax
  double precision, intent(in) ::  verts(nvert,2)
  integer, intent(in) ::  edges(nedge0,2)

  !Output
  integer, intent(out) ::  nedge

  !Working
  logical xTrue, yTrue
  integer iedge
  double precision x, y

  nedge = 0
  do iedge=1,nedge0
     x = 0.5*verts(edges(iedge,1),1) + 0.5*verts(edges(iedge,2),1)
     y = 0.5*verts(edges(iedge,1),2) + 0.5*verts(edges(iedge,2),2)
     xTrue = (xmin - 1e-10 .le. x) .and. (x .le. xmax + 1e-10)
     yTrue = (ymin - 1e-10 .le. y) .and. (y .le. ymax + 1e-10)
     if (xTrue .and. yTrue) then
        nedge = nedge + 1
     end if
  end do

end subroutine countSurfaceEdges
