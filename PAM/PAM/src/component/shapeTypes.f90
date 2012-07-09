subroutine getShapeRect(n, n0, n1, u0, v0)

  implicit none

  !Input
  integer, intent(in) ::  n, n0, n1

  !Output
  double precision, intent(out) ::  u0(n), v0(n)

  !Working
  integer i, j, index
  double precision uu, vv

  uu = 1.0/(n0-1)
  vv = 1.0/(n1-1)

  index = 0
  do j=1,n1
     do i=1,n0
        index = index + 1
        u0(index) = (i-1)*uu
        v0(index) = (j-1)*vv
     end do
  end do

end subroutine getShapeRect



subroutine getShapeHole(n, n0, n1, n2, u0, v0)

  implicit none

  !Input
  integer, intent(in) ::  n, n0, n1, n2

  !Output
  double precision, intent(out) ::  u0(n), v0(n)

  !Working
  integer i

  u0(:) = 0.0
  v0(:) = 0.0

end subroutine getShapeHole
