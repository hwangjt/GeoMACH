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
  integer index, i,j
  double precision x0(n0), y0(n0), x1(n1), y1(n1)
  double precision rz, ry, ta, tb, t1, t2

  rz = 0.25
  ry = 0.25

  u0(:) = 0.0
  v0(:) = 0.0

  ta = 0.0
  tb = 1.0
  
  index = 0

  t1 = -0.25
  t2 = 0.25
  call computeRoundedSection(n1, rz, ry, ta, tb, ta, tb, t1, t2, x1, y1)
  x1(:) = x1(:) + 0.5
  y1(:) = y1(:) + 0.5
  do j=1,n2
     do i=1,n1
        index = index + 1
        u0(index) = x1(i) + (1.0 - x1(i))*(j-1)/(n2-1)
        v0(index) = y1(i) + (1.0*(i-1)/(n1-1) - y1(i))*(j-1)/(n2-1)
     end do
  end do

  t1 = 0.25
  t2 = 0.75
  call computeRoundedSection(n0, rz, ry, ta, tb, ta, tb, t1, t2, x0, y0)
  x0(:) = x0(:) + 0.5
  y0(:) = y0(:) + 0.5
  do j=1,n2
     do i=1,n0
        index = index + 1
        u0(index) = x0(i) + (1.0*(n0-i)/(n0-1) - x0(i))*(j-1)/(n2-1)
        v0(index) = y0(i) + (1.0 - y0(i))*(j-1)/(n2-1)
     end do
  end do

  t1 = 0.75
  t2 = 1.25
  call computeRoundedSection(n1, rz, ry, ta, tb, ta, tb, t1, t2, x1, y1)
  x1(:) = x1(:) + 0.5
  y1(:) = y1(:) + 0.5
  do j=1,n2
     do i=1,n1
        index = index + 1
        u0(index) = x1(i) + (0.0 - x1(i))*(j-1)/(n2-1)
        v0(index) = y1(i) + (1.0*(n1-i)/(n1-1) - y1(i))*(j-1)/(n2-1)
     end do
  end do

  t1 = 1.25
  t2 = 1.75
  call computeRoundedSection(n0, rz, ry, ta, tb, ta, tb, t1, t2, x0, y0)
  x0(:) = x0(:) + 0.5
  y0(:) = y0(:) + 0.5
  do j=1,n2
     do i=1,n0
        index = index + 1
        u0(index) = x0(i) + (1.0*(i-1)/(n0-1) - x0(i))*(j-1)/(n2-1)
        v0(index) = y0(i) + (0.0 - y0(i))*(j-1)/(n2-1)
     end do
  end do

end subroutine getShapeHole
