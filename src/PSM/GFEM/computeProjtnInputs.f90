subroutine computeProjtnInputs(nvert, ni, nj, verts, face_surf, P0, surfs, Q)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nvert, ni, nj, verts, face_surf
  !f2py intent(out) P0, surfs, Q
  !f2py depend(nvert) verts, P0, surfs, Q
  !f2py depend(ni,nj) face_surf

  !Input
  integer, intent(in) ::  nvert, ni, nj
  double precision, intent(in) ::  verts(nvert,2)
  integer, intent(in) ::  face_surf(ni,nj)

  !Output
  double precision, intent(out) ::  P0(nvert,3)
  integer, intent(out) ::  surfs(nvert)
  double precision, intent(out) ::  Q(nvert,3)

  !Working
  integer ivert, i, j

  P0(:,1:2) = verts
  P0(:,3) = 0.0

  Q(:,1:2) = 0.0
  Q(:,3) = 1.0

  do ivert=1,nvert
     i = ceiling(ni*verts(ivert,1))
     j = ceiling(nj*verts(ivert,2))
     if (i .eq. 0) then
        i = 1
     end if
     if (j .eq. 0) then
        j = 1
     end if
     surfs(ivert) = face_surf(i,j)
  end do

end subroutine computeProjtnInputs
