subroutine createSurfaces(n, du0, dv0, x0, x1, z0, z1, d, u1, u2, v1, v2, P)

  ! Computes an array, P, representing a rectangular plane
  ! n: number of points per edge
  ! du0, dv0: {1,2,3} maps to {x,y,z}; negative sign means reverse order
  ! x0,x1,z0,z1: set y to zero if x or z is 0 or 1, respectively
  ! d: translate the remaining coordinate by d
  ! u1,u2,v1,v2: linspace(u1,u2,n) and linspace(v1,v2,n)

  !Fortran-python interface directives
  !f2py intent(in) n, du0, dv0, x0, x1, z0, z1, d, u1, u2, v1, v2
  !f2py intent(out) P
  !f2py depend(n) P

  !Input
  integer, intent(in) ::  n, du0, dv0, x0, x1, z0, z1
  double precision, intent(in) ::  d, u1, u2, v1, v2

  !Output
  double precision, intent(out) ::  P(n,n,3)

  !Working
  integer u, v, du, dv
  double precision uu,vv
  logical val

  du = abs(du0)
  dv = abs(dv0)
  uu = (u2 - u1)/(n-1)
  vv = (v2 - v1)/(n-1)

  P(:,:,:) = d
  do u=1,n
     do v=1,n
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

  do u=1,n
     do v=1,n
        val = .False.
        val = val .or. ((x0 .eq. 1) .and. (P(u,v,1) .eq. 0.0))
        val = val .or. ((x1 .eq. 1) .and. (P(u,v,1) .eq. 1.0))
        val = val .or. ((z0 .eq. 1) .and. (P(u,v,3) .eq. 0.0))
        val = val .or. ((z1 .eq. 1) .and. (P(u,v,3) .eq. 1.0))
        if (val) then
           P(u,v,2) = 0.0
        end if
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
