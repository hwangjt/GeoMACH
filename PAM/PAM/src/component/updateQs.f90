subroutine updateQs(nQ, ni, nj, Nf, Qf, Q)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nQ, ni, nj, Nf, Qf
  !f2py intent(in,out) Q
  !f2py depend(ni,nj) Nf
  !f2py depend(ni,nj) Qf
  !f2py depend(nQ) Q

  !Input
  integer, intent(in) ::  nQ, ni, nj
  integer, intent(in) ::  Nf(ni, nj, 5)
  complex*16, intent(in) ::  Qf(ni, nj, 3)

  !Output
  double precision, intent(inout) ::  Q(nQ,3)

  !Working
  integer i, j, k

  do j=1,nj
     do i=1,ni
        if (Nf(i,j,1) .ne. -1) then
           do k=1,3
              Q(Nf(i,j,1)+1,k) = realpart(Qf(i,j,k))
           end do
        end if
     end do
  end do

end subroutine updateQs



subroutine computeRtnMtx(rot, T)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) rot
  !f2py intent(out) T

  !Input
  double precision, intent(in) ::  rot(3)

  !Output
  double precision, intent(out) ::  T(3,3)

  !Working
  double precision p, q, r, T0(3,3), pi
  integer k

  pi = 2*acos(0.0)

  p = rot(1)*pi/180.0
  q = rot(2)*pi/180.0
  r = rot(3)*pi/180.0

  T(:,:) = 0.0
  do k=1,3
     T(k,k) = 1.0
  end do

  T0(:,:) = 0.0
  T0(1,1) = 1.0
  T0(2,2) = cos(p)
  T0(2,3) = sin(p)
  T0(3,2) = -sin(p)
  T0(3,3) = cos(p)
  T = matmul(T0, T)

  T0(:,:) = 0.0
  T0(2,2) = 1.0
  T0(1,1) = cos(p)
  T0(1,3) = sin(p)
  T0(3,1) = -sin(p)
  T0(3,3) = cos(p)
  T = matmul(T0, T)

  T0(:,:) = 0.0
  T0(3,3) = 1.0
  T0(1,1) = cos(p)
  T0(1,2) = sin(p)
  T0(2,1) = -sin(p)
  T0(2,2) = cos(p)
  T = matmul(T0, T)

end subroutine computeRtnMtx
