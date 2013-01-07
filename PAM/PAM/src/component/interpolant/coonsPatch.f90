subroutine coonsPatch(nu, nv, Pu0, Pu1, P0v, P1v, P, dPdw)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nu, nv, Pu0, Pu1, P0v, P1v
  !f2py intent(out) P
  !f2py depend(nu) Pu0, Pu1
  !f2py depend(nv) P0v, P1v
  !f2py depend(nu,nv) P, dPdw

  !Input
  integer, intent(in) ::  nu, nv
  double precision, intent(in) ::  Pu0(nu,3), Pu1(nu,3), P0v(nv,3), P1v(nv,3)

  !Output
  double precision, intent(out) ::  P(nu,nv,3), dPdw(nu,nv,3)

  !Working
  double precision P00(3), P01(3), P10(3), P11(3)
  double precision dPdu(3), dPdv(3), norm
  double precision denu, denv
  double precision u, v
  integer i, j

  P00 = P0v(1,:)
  P10 = P1v(1,:)
  P01 = P0v(nv,:)
  P11 = P1v(nv,:)

  denu = 1.0/(nu-1)
  denv = 1.0/(nv-1)

  do i=1,nu
     u = (i-1)*denu
     do j=1,nv
        v = (j-1)*denv
        P(i,j,:) = (1-u)*P0v(j,:) + u*P1v(j,:) + (1-v)*Pu0(i,:) + v*Pu1(i,:)
        P(i,j,:) = P(i,j,:) - (1-u)*(1-v)*P00 - u*(1-v)*P10 - (1-u)*v*P01 - u*v*P11
        dPdu = -P0v(j,:) + P1v(j,:) + (1-v)*P00 - (1-v)*P10 + v*P01 - v*P11
        dPdv = -Pu0(i,:) + Pu1(i,:) + (1-u)*P00 + u*P10 - (1-u)*P01 - u*P11
        call cross_product(dPdu, dPdv, dPdw(i,j,:))
        norm = dot_product(dPdw(i,j,:),dPdw(i,j,:))**0.5
        if (norm .gt. 1e-10) then
           dPdw(i,j,:) = dPdw(i,j,:)/norm
        end if
     end do
  end do

end subroutine coonsPatch



subroutine cross_product(u, v, w)

  implicit none

  !Input
  double precision, intent(in) ::  u(3), v(3)

  !Output
  double precision, intent(out) ::  w(3)

  !Working
  double precision cross_ij

  w(1) = cross_ij(2, 3, u, v)
  w(2) = cross_ij(3, 1, u, v)
  w(3) = cross_ij(1, 2, u, v)

end subroutine cross_product



function cross_ij(i, j, u, v)
  
  implicit none

  integer, intent(in) ::  i, j
  double precision, intent(in) ::  u(3), v(3)
  double precision ::  cross_ij

  cross_ij = u(i)*v(j) - u(j)*v(i)

end function cross_ij
