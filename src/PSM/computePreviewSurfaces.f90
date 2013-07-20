subroutine computePreviewSurfaces(nnode, nsurf, quads, s, u, v)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nnode, nsurf
  !f2py intent(out) quads, s, u, v
  !f2py depend(nsurf) quads
  !f2py depend(nnode) s, u, v

  !Input
  integer, intent(in) ::  nnode, nsurf

  !Output
  integer, intent(out) ::  quads(nsurf,4), s(nnode)
  double precision, intent(out) ::  u(nnode), v(nnode)

  !Working
  integer isurf, i1, i2

  quads(:,:) = 0
  s(:) = 0
  u(:) = 0.
  v(:) = 0.
  i1 = 1
  i2 = 4
  do isurf=1,nsurf
     quads(isurf,:) = (/ i1, i1+1, i1+3, i1+2 /)
     s(i1:i2) = isurf - 1
     u(i1:i2) = (/ dble(0), dble(1), dble(0), dble(1) /)
     v(i1:i2) = (/ dble(0), dble(0), dble(1), dble(1) /)
     i1 = i1 + 4
     i2 = i2 + 4
  end do

end subroutine computePreviewSurfaces
