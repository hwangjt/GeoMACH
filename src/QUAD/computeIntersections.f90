subroutine computeIntersections(nvert0, nedge, nvert, verts0, edges, verts)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nvert0, nedge, nvert, verts0, edges
  !f2py intent(out) verts
  !f2py depend(nvert0) verts0
  !f2py depend(nedge) edges
  !f2py depend(nvert) verts

  !Input
  integer, intent(in) ::  nvert0, nedge, nvert
  double precision, intent(in) ::  verts0(nvert0,2)
  integer, intent(in) ::  edges(nedge,2)

  !Output
  double precision, intent(out) ::  verts(nvert,2)

  !Working
  integer i, i1, i2
  logical valid
  double precision v(2)

  verts(1:nvert0,:) = verts0(:,:)

  i = nvert0 + 1
  do i1=1,nedge
     do i2=i1+1,nedge
        call intersect(verts0(edges(i1,1),:), verts0(edges(i1,2),:), &
             verts0(edges(i2,1),:), verts0(edges(i2,2),:), valid, v)
        if (valid) then
           verts(i,:) = v
           i = i + 1
        end if
     end do
  end do

end subroutine computeIntersections




subroutine countIntersections(nvert, nedge, verts, edges, nint)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nvert, nedge, verts, edges
  !f2py intent(out) nint
  !f2py depend(nvert) verts
  !f2py depend(nedge) edges

  !Input
  integer, intent(in) ::  nvert, nedge
  double precision, intent(in) ::  verts(nvert,2)
  integer, intent(in) ::  edges(nedge,2)

  !Output
  integer, intent(out) ::  nint

  !Working
  integer i1, i2
  logical valid
  double precision v(2)

  nint = 0
  do i1=1,nedge
     do i2=i1+1,nedge
        call intersect(verts(edges(i1,1),:), verts(edges(i1,2),:), &
             verts(edges(i2,1),:), verts(edges(i2,2),:), valid, v)
        if (valid) then
           nint = nint + 1
        end if
     end do
  end do

end subroutine countIntersections




subroutine intersect(a1, b1, a2, b2, valid, v)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) a1, b1, a2, b2
  !f2py intent(out) valid, v

  !Input
  double precision, intent(in) ::  a1(2), b1(2), a2(2), b2(2)

  !Output
  logical valid
  double precision v(2)

  !Working
  double precision d1(2), d2(2)
  double precision det, mtx(2,2), inv(2,2), rhs(2)

  !a1 + v1 (b1-a1) = a2 + v2 (b2-a2)

  d1 = b1 - a1
  d2 = b2 - a2
  rhs = a2 - a1
  mtx(:,1) = d1
  mtx(:,2) = -d2
  det = mtx(1,1)*mtx(2,2) - mtx(1,2)*mtx(2,1)

  valid = .False.
  v(:) = 0.0
  if (abs(det) .gt. 1e-12) then
     inv(1,1) = mtx(2,2)/det
     inv(2,2) = mtx(1,1)/det
     inv(1,2) = -mtx(1,2)/det
     inv(2,1) = -mtx(2,1)/det
     v(:) = matmul(inv,rhs)
     if ((v(1) .ge. -1e-12) .and. (v(1) .le. 1+1e-12) .and. &
          (v(2) .ge. -1e-12) .and. (v(2) .le. 1+1e-12)) then
        valid = .True.
        v(:) = a1 + v(1)*d1
     end if
  end if

end subroutine intersect
