subroutine interpolateFrames(i, nu0, nv0, iu, iv, Q, dQdw)

  implicit none

  !Input
  integer, intent(in) ::  i, nu0, nv0, iu(4), iv(4)

  !Output
  double precision, intent(inout) ::  Q(nu0,nv0,3)
  double precision, intent(out) ::  dQdw(nu0,nv0,3)

  !Working
  integer j, nu, nv
  
  dQdw(:,:,:) = 0.0
  do j=1,3
     nv = iv(j+1)-iv(j)
     nu = iu(i+1)-iu(i)
     if ((nu.gt.0) .and. (nv.gt.0)) then
        call coonsPatch(nu+1, nv+1, &
             Q(iu(i):iu(i+1),iv(j),:), &
             Q(iu(i):iu(i+1),iv(j+1),:), &
             Q(iu(i),iv(j):iv(j+1),:), &
             Q(iu(i+1),iv(j):iv(j+1),:), &
             Q(iu(i):iu(i+1),iv(j):iv(j+1),:), &
             dQdw(iu(i):iu(i+1),iv(j):iv(j+1),:))
     end if
  end do  

end subroutine interpolateFrames



subroutine bezierCurve(n, A, S, B, T, P)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) n, A, S, B, T
  !f2py intent(out) P
  !f2py depend(n) P

  !Input
  integer, intent(in) ::  n
  double precision, intent(in) ::  A(3), S(3), B(3), T(3)

  !Output
  double precision, intent(out) ::  P(n,3)

  !Working
  double precision den, u
  double precision P0(3), P1(3), P2(3), P3(3)
  logical Sgiven, Tgiven
  integer i

  den = 1.0/(n-1)

  if (dot_product(S,S)**0.5 .lt. 1e-10) then
     Sgiven = .false.
  else
     Sgiven = .true.
  end if

  if (dot_product(T,T)**0.5 .lt. 1e-10) then
     Tgiven = .false.
  else
     Tgiven = .true.
  end if

  if ((.not. Sgiven) .and. (.not. Tgiven)) then
     P0 = A
     P1 = B
     P2 = 0
     P3 = 0
  elseif ((Sgiven) .and. (.not. Tgiven)) then
     P0 = A
     P1 = S/2.0 + A
     P2 = B
     P3 = 0
  elseif ((.not. Sgiven) .and. (Tgiven)) then
     P0 = A
     P1 = T/2.0 + B
     P2 = B
     P3 = 0
  else
     P0 = A
     P1 = S/3.0 + A
     P2 = T/3.0 + B
     P3 = B
  end if

  do i=1,n
     u = (i-1)*den
     if ((.not. Sgiven) .and. (.not. Tgiven)) then
        P(i,:) = (1-u)*P0 + u*P1
     elseif (Sgiven .and. Tgiven) then
        P(i,:) = (1-u)**3*P0 + 3*u*(1-u)**2*P1 + 3*u**2*(1-u)*P2 + u**3*P3
     else
        P(i,:) = (1-u)**2*P0 + 2*u*(1-u)*P1 + u**2*P2
     end if
  end do

end subroutine bezierCurve
