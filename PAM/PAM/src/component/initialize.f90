subroutine createSurfaces(nu, nv, du0, dv0, d, u1, u2, v1, v2, P)

  !Fortran-python interface directives
  !f2py intent(in) nu, nv, du0, dv0, d, u1, u2, v1, v2
  !f2py intent(out) P
  !f2py depend(nu,nv) P

  !Input
  integer, intent(in) ::  nu, nv, du0, dv0
  double precision, intent(in) ::  d, u1, u2, v1, v2

  !Output
  double precision, intent(out) ::  P(nu,nv,3)

  !Working
  integer u, v, du, dv
  double precision uu,vv

  du = abs(du0)
  dv = abs(dv0)
  uu = (u2 - u1)/(nu-1)
  vv = (v2 - v1)/(nv-1)

  P(:,:,:) = d
  do u=1,nu
     do v=1,nv
        P(u,v,du) = u1 + uu*(u-1)
        P(u,v,dv) = v1 + vv*(v-1)
     end do
  end do

  if (du0 .lt. 0) then
     P(:,:,du) = 1 - P(:,:,du)
  end if
  if (dv0 .lt. 0) then
     P(:,:,dv) = 1 - P(:,:,dv)
  end if

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
