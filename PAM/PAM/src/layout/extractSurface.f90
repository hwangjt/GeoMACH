subroutine extractFlattened(n, nP, nvert, npoly, verts, poly_vert, P)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) n, nP, nvert, npoly, verts, poly_vert
  !f2py intent(out) P
  !f2py depend(nvert) verts
  !f2py depend(npoly) poly_vert
  !f2py depend(nP) P

  !Input
  integer, intent(in) ::  n, nP, nvert, npoly
  double precision, intent(in) ::  verts(nvert,2)
  integer, intent(in) ::  poly_vert(npoly,5)

  !Output
  double precision, intent(out) ::  P(nP,3)

  !Working
  integer iP, u, v, poly
  double precision C1(2), C2(2), C3(2), C4(2)

  P(:,:) = 0.0

  iP = 0
  do poly=1,npoly
     C1 = verts(poly_vert(poly,1),:)
     C2 = verts(poly_vert(poly,2),:)
     C3 = verts(poly_vert(poly,3),:)
     C4 = verts(poly_vert(poly,4),:)
     do u=1,n
        do v=1,n
           iP = iP + 1
           P(iP,1:2) = (n-u)*(n-v)*C1 + (n-u)*(v-1)*C2 &
                + (u-1)*(v-1)*C3 + (u-1)*(n-v)*C4 
           P(iP,:) = P(iP,:)/(n-1)**2
        end do
     end do
  end do

end subroutine extractFlattened



subroutine extractSurface(poly, n, nvert, npoly, verts, poly_vert, P)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) poly, n, nvert, npoly, verts, poly_vert
  !f2py intent(out) P
  !f2py depend(nvert) verts
  !f2py depend(npoly) poly_vert
  !f2py depend(n) P

  !Input
  integer, intent(in) ::  poly, n, nvert, npoly
  double precision, intent(in) ::  verts(nvert,2)
  integer, intent(in) ::  poly_vert(npoly,5)

  !Output
  double precision, intent(out) ::  P(n,n,2)

  !Working
  integer u, v
  double precision C1(2), C2(2), C3(2), C4(2)

  C1 = verts(poly_vert(poly,1),:)
  C2 = verts(poly_vert(poly,2),:)
  C3 = verts(poly_vert(poly,3),:)
  C4 = verts(poly_vert(poly,4),:)
  do u=1,n
     do v=1,n
        P(u,v,:) = (n-u)*(n-v)*C1 + (n-u)*(v-1)*C2 &
             + (u-1)*(v-1)*C3 + (u-1)*(n-v)*C4 
        P(u,v,:) = P(u,v,:)/(n-1)**2
     end do
  end do

end subroutine extractSurface
  
  
  
