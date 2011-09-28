subroutine evalcurve(k, m, n, B, i0, C, P)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) k,m,n,B,i0,C
  !f2py intent(out) P
  !f2py depend(n,k) B
  !f2py depend(n) i0
  !f2py depend(m) C
  !f2py depend(n) P

  !Input
  integer, intent(in) ::  k, m, n
  double precision, intent(in) ::  B(n,k)
  integer, intent(in) ::  i0(n)
  double precision, intent(in) ::  C(m)

  !Output
  double precision, intent(out) ::  P(n)

  !Working
  integer u, i 

  P(:) = 0
  do i=1,k
     do u=1,n
        P(u) = P(u) + B(u,i)*C(i0(u)+i)
     end do
  end do

end subroutine evalcurve


subroutine curvejacob(k, m, n, B, i0, dPdC)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) k,m,n,B,i0
  !f2py intent(out) dPdC
  !f2py depend(n,k) B
  !f2py depend(n) i0
  !f2py depend(n,m) dPdC

  !Input
  integer, intent(in) ::  k, m, n
  double precision, intent(in) ::  B(n,k)
  integer, intent(in) ::  i0(n)

  !Output
  double precision, intent(out) ::  dPdC(n,m)

  !Working
  integer u,i

  dPdC(:,:) = 0
  do i=1,k
     do u=1,n
        dPdC(u,i0(u)+i) = B(u,i)
     end do
  end do

end subroutine curvejacob


subroutine curvefit(k, m, n, B, i0, P0, dPdC, P)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) k,m,n,B,i0,P0
  !f2py intent(out) dPdC,P
  !f2py depend(n,k) B
  !f2py depend(n) i0
  !f2py depend(n) P0
  !f2py depend(n,m) dPdC
  !f2py depend(n) P

  !Input
  integer, intent(in) ::  k, m, n
  double precision, intent(in) ::  B(n,k)
  integer, intent(in) ::  i0(n)
  double precision, intent(in) ::  P0(n,3)

  !Output
  double precision, intent(out) ::  dPdC(n-2,m-2)
  double precision, intent(out) ::  P(n-2,3)

  !Working
  integer u,i
  double precision dPdC0(n,m)

  dPdC0(:,:) = 0.0
  do i=1,k
     do u=1,n
        dPdC0(u,i0(u)+i) = B(u,i)
     end do
  end do
  do i=1,m-2
     do u=1,n-2
        dPdC(u,i) = dPdC0(u+1,i+1)
     end do
  end do
  do i=1,3
     do u=1,n-2
        P(u,i) = P0(u+1,i) - P0(1,i)*dPdC0(u+1,1) - P0(n,i)*dPdC0(u+1,m)
     end do
  end do

end subroutine curvefit
