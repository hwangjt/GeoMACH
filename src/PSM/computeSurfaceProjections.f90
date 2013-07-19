subroutine computeSurfaceProjections(nvert, verts, P, Q)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nvert, verts
  !f2py intent(out) P, Q
  !f2py depend(nvert) verts, P, Q

  !Input
  integer, intent(in) ::  nvert
  double precision, intent(in) ::  verts(nvert,2)

  !Output
  double precision, intent(out) ::  P(nvert,3), Q(nvert,3)

  !Working
  integer i
  double precision xmin, xmax, ymin, ymax
  
  xmin = minval(verts(:,1))
  xmax = maxval(verts(:,1))
  ymin = minval(verts(:,2))
  ymax = maxval(verts(:,2))

  P(:,:) = 0.
  Q(:,:) = 0.
  Q(:,3) = 1.
  do i=1,nvert
     P(i,1) = (verts(i,1)-xmin)/(xmax-xmin)
     P(i,2) = (verts(i,2)-ymin)/(ymax-ymin)
  end do

end subroutine computeSurfaceProjections
