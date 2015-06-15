subroutine computeAngles(ax1, ax2, t, nor, rot, drot_dt)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) ax1, ax2, t, nor
  !f2py intent(out) rot, drot_dt

  !Input
  integer, intent(in) ::  ax1, ax2
  double precision, intent(in) ::  t(3), nor(3)

  !Output
  double precision, intent(out) ::  rot(3), drot_dt(3,3)

  !Working
  integer ax3
  double precision vnorm, v(2), w(2), dv_dt(2,3), dw_dt(2,3)
  double precision p, q, dp_dv(2), dq_dw(2)

  ax3 = 6 - ax1 - ax2
  v = (/t(ax1),t(ax2)/)
  vnorm = dot_product(v,v)**0.5
  if (vnorm .lt. 1e-12) then
     vnorm = 1.0
  end if
  w = (/vnorm,-t(ax3)/)
  dv_dt(:,:) = 0.0
  dw_dt(:,:) = 0.0
  dv_dt(1,ax1) = 1.0
  dv_dt(2,ax2) = 1.0
  dw_dt(1,ax1) = t(ax1)/vnorm
  dw_dt(1,ax2) = t(ax2)/vnorm
  dw_dt(2,ax3) = -1.0
  if ((ax2-ax1 .eq. -1) .or. (ax2-ax1 .eq. 2)) then
     v(2) = -v(2)
     w(2) = -w(2)
     dv_dt(2,:) = -dv_dt(2,:)
     dw_dt(2,:) = -dw_dt(2,:)
  end if
  call arctan2pi(v, p, dp_dv)
  call arctan2pi(w, q, dq_dw)
  rot(1) = p*nor(1)
  rot(2) = q*nor(2)
  rot(3) = 0.0
  drot_dt(1,:) = matmul(dp_dv, dv_dt)*nor(1)
  drot_dt(2,:) = matmul(dq_dw, dw_dt)*nor(2)
  drot_dt(3,:) = 0.0

end subroutine computeAngles



subroutine arctan2pi(P, t, dt_dP)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) P
  !f2py intent(out) t, dt_dP

  !Input
  double precision, intent(in) ::  P(2)

  !Output
  double precision, intent(out) ::  t, dt_dP(2)

  !Working
  double precision x, y, pi

  pi = 2*acos(0.0)

  x = P(1)
  y = P(2)

  if (dot_product(P,P)**0.5 .lt. 1e-10) then
     t = 0
  else if (x .eq. 0) then
     if (y .gt. 0) then
        t = pi/2.0
     else if (y .lt. 0) then
        t = 3*pi/2.0
     end if
  else if (y .eq. 0) then
     if (x .gt. 0) then
        t = 0
     else if (x .lt. 0) then
        t = pi
     end if
  else if (x .lt. 0) then
     t = atan(y/x) + pi   
  else if (y .lt. 0) then
     t = atan(y/x) + 2*pi    
  else if (y .gt. 0) then
     t = atan(y/x)
  else
     t = 0
  end if

  if (dot_product(P,P)**0.5 .lt. 1e-10) then
     dt_dP(:) = 0.0
  else
     dt_dP(1) = -y/(x**2+y**2)
     dt_dP(2) = x/(x**2+y**2)
  end if

end subroutine arctan2pi
