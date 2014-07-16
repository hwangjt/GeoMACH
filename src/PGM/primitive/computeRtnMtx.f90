subroutine computeRtnMtx(ax1, ax2, ax3, rot, T, dT_drot)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) ax1, ax2, ax3, rot
  !f2py intent(out) T, dT_drot

  !Input
  integer, intent(in) ::  ax1, ax2, ax3
  double precision, intent(in) ::  rot(3)

  !Output
  double precision, intent(out) ::  T(3,3), dT_drot(3,3,3)

  !Working
  double precision T1(3,3), T2(3,3), T3(3,3)
  double precision dT1(3,3), dT2(3,3), dT3(3,3)

  call computeRtn(ax3, rot(1), T1, dT1)
  call computeRtn(ax2, rot(2), T2, dT2)
  call computeRtn(ax1, rot(3), T3, dT3)

  T = matmul(matmul(T1,T2),T3)
  dT_drot(:,:,1) = matmul(matmul(dT1,T2),T3)
  dT_drot(:,:,2) = matmul(matmul(T1,dT2),T3)
  dT_drot(:,:,3) = matmul(matmul(T1,T2),dT3)

end subroutine computeRtnMtx



subroutine computeRtn(k, p, T, dT)

  implicit none

  !Input
  integer, intent(in) ::  k
  double precision, intent(in) ::  p

  !Output
  double precision, intent(out) ::  T(3,3), dT(3,3)

  !Working
  integer i, j

  i = k + 1
  j = k + 2
  if (i .gt. 3) then
     i = i - 3
  end if
  if (j .gt. 3) then
     j = j - 3
  end if

  T(:,:) = 0.0
  T(k,k) = 1.0
  T(i,i) = cos(p)
  T(i,j) = -sin(p)
  T(j,i) = sin(p)
  T(j,j) = cos(p)

  dT(:,:) = 0.0
  dT(i,i) = -sin(p)
  dT(i,j) = -cos(p)
  dT(j,i) = cos(p)
  dT(j,j) = -sin(p)

end subroutine computeRtn
