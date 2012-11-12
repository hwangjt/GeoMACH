subroutine createSurfaces(n, du, dv, d, u1, u2, v1, v2, P)

  ! Computes an array, P, representing a rectangular plane
  ! n: number of points per edge
  ! du, dv: {1,2,3} maps to {x,y,z}; negative sign means reverse order
  ! d: translate the remaining coordinate by d
  ! u1,u2,v1,v2: linspace(u1,u2,n) and linspace(v1,v2,n)

  !Fortran-python interface directives
  !f2py intent(in) n, du, dv, d, u1, u2, v1, v2
  !f2py intent(out) P
  !f2py depend(n) P

  !Input
  integer, intent(in) ::  n, du, dv
  double precision, intent(in) ::  d, u1, u2, v1, v2

  !Output
  double precision, intent(out) ::  P(n,n,3)

  !Working
  integer u, v
  double precision uu,vv

  uu = (u2 - u1)/(n-1)
  vv = (v2 - v1)/(n-1)

  P(:,:,:) = d
  do u=1,n
     do v=1,n
        P(u,v,abs(du)) = u1 + uu*(u-1)
        P(u,v,abs(dv)) = v1 + vv*(v-1)
     end do
  end do

end subroutine createSurfaces



subroutine createInterface(nu, nv, edge1, edge2, P)

  !Fortran-python interface directives
  !f2py intent(in) nu, nv, edge1, edge2
  !f2py intent(out) P
  !f2py depend(nv) edge1, edge2
  !f2py depend(nu,nv) P

  !Input
  integer, intent(in) ::  nu, nv
  double precision, intent(in) ::  edge1(nv,3), edge2(nv,3)

  !Output
  double precision, intent(out) ::  P(nu,nv,3)

  !Working
  integer u, v
  double precision e1(3), ee(3)

  do v=1,nv
     e1 = edge1(v,:)
     ee = (edge2(v,:) - edge1(v,:))/(nu-1)
     do u=1,nu
        P(u,v,:) = e1 + ee*(u-1)
     end do
  end do

end subroutine createInterface
