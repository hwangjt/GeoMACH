subroutine computeParameter(mu, mv, nu, nv, P, Tu, Tv, Q)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) mu, mv, nu, nv, P, Tu, Tv
  !f2py intent(out) Q
  !f2py depend(mu,mv) P
  !f2py depend(mu) Tu
  !f2py depend(mv) Tv
  !f2py depend(nu,nv) Q

  !Input
  integer, intent(in) ::  mu, mv, nu, nv
  double precision, intent(in) ::  P(mu,mv,5), Tu(mu), Tv(mv)

  !Output
  double precision, intent(out) ::  Q(nu,nv)

  !Working
  integer i, j, i0, j0
  double precision u, v, u0, v0
  double precision Pu0, Pu1, P0v, P1v
  double precision evalCubic

  Q(:,:) = 0.0

  do i0=1,nu
     u0 = 1.0*(i0-1)/(nu-1)
     call locateParameter(mu, u0, Tu, i, u)
     do j0=1,nv
        v0 = 1.0*(j0-1)/(nv-1)
        call locateParameter(mv, v0, Tv, j, v)
        if ((i .gt. 0) .and. (j .gt. 0)) then
           Pu0 = evalCubic(1, u, P(i,j,:), P(i+1,j,:))
           Pu1 = evalCubic(1, u, P(i,j+1,:), P(i+1,j+1,:))
           P0v = evalCubic(2, v, P(i,j,:), P(i,j+1,:))
           P1v = evalCubic(2, v, P(i+1,j,:), P(i+1,j+1,:))
           Q(i0,j0) = (1-u)*P0v + u*P1v + (1-v)*Pu0 + v*Pu1 &
                - (1-u)*(1-v)*P(i,j,1) - u*(1-v)*P(i+1,j,1) &
                - (1-u)*v*P(i,j+1,1) - u*v*P(i+1,j+1,1)
        else if ((i .eq. 0) .and. (j .eq. 0)) then
           Q(i0,j0) = P(1,1,1)
        else if ((i .eq. 0) .and. (j .gt. 0)) then
           i = 1
           Q(i0,j0) = evalCubic(2, v, P(i,j,:), P(i,j+1,:))
        else if ((i .gt. 0) .and. (j .eq. 0)) then
           j = 1
           Q(i0,j0) = evalCubic(1, u, P(i,j,:), P(i+1,j,:))
        end if
     end do
  end do

end subroutine computeParameter



function evalCubic(dim, t, V0, V1)
  
  implicit none

  integer, intent(in) ::  dim
  double precision, intent(in) ::  t, V0(5), V1(5)
  double precision ::  evalCubic
  double precision P0, P1, D0, D1
  logical B0, B1
  double precision A, B, C, D

  P0 = V0(1)
  D0 = V0(1+dim)
  B0 = V0(3+dim) .gt. 0.5

  P1 = V1(1)
  D1 = V1(1+dim)
  B1 = V1(3+dim) .gt. 0.5

  if (B0 .and. B1) then
     A = P0
     B = D0/3.0 + P0
     C = -D1/3.0 + P1
     D = P1
  else if (B0) then
     A = P0
     B = D0/3.0 + P0
     C = D0/3.0 + 2.0/3.0*P0 + P1/3.0
     D = P1
  else if (B1) then
     A = P0
     B = -D1/3.0 + 2.0/3.0*P1 + P0/3.0
     C = -D1/3.0 + P1
     D = P1
  else
     A = P0
     B = 2.0/3.0*P0 + P1/3.0
     C = 2.0/3.0*P1 + P0/3.0
     D = P1
  end if
  
  evalCubic = (1-t)**3*A + 3*t*(1-t)**2*B + 3*t**2*(1-t)*C + t**3*D

end function evalCubic



subroutine locateParameter(n, u0, T, i, u)

  implicit none

  !Input
  integer, intent(in) ::  n
  double precision, intent(in) ::  u0, T(n)

  !Output
  integer, intent(out) ::  i
  double precision, intent(out) ::  u

  !Working
  integer k

  i = -1
  u = 0.0
  do k=1,n-1
     if ((T(k) .le. u0) .and. (u0 .le. T(k+1))) then
        i = k
        u = (u0-T(k))/(T(k+1)-T(k))
     end if
  end do
  if (n .eq. 1) then
     i = 0
  end if

end subroutine locateParameter
