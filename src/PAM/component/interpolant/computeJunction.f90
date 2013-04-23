subroutine computeJunction(nu, nv, nu1, nu2, nu3, nv1, nv2, nv3, &
     f0, m0, W, E, N, S, fQ, shape0, Q)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nu, nv, nu1, nu2, nu3, nv1, nv2, nv3, f0, m0, W, E, N, S, fQ, shape0
  !f2py intent(out) Q
  !f2py depend(nu2) W, E
  !f2py depend(nv2) N, S
  !f2py depend(nu,nv) fQ, shape0, Q

  !Input
  integer, intent(in) ::  nu, nv, nu1, nu2, nu3, nv1, nv2, nv3
  double precision, intent(in) ::  f0, m0
  double precision, intent(in) ::  W(nu2,2,3), E(nu2,2,3)
  double precision, intent(in) ::  N(nv2,2,3), S(nv2,2,3)
  double precision, intent(in) ::  fQ(nu,nv,3), shape0(nu,nv)

  !Output
  double precision, intent(out) ::  Q(nu,nv,3)

  !Working
  integer i, j, k, iu(4), iv(4)
  double precision verts(4,4,2,3), dQdw(nu,nv,3)

  iu(1) = 2
  iu(2) = (nu1-1) + 1
  iu(3) = (nu1-1) + (nu2-1) + 1
  iu(4) = (nu1-1) + (nu2-1) + (nu3-1) + 0

  iv(1) = 2
  iv(2) = (nv1-1) + 1
  iv(3) = (nv1-1) + (nv2-1) + 1
  iv(4) = (nv1-1) + (nv2-1) + (nv3-1) + 0

  Q = fQ
  Q(iu(2),iv(2):iv(3),:) = N(:,1,:)
  Q(iu(3),iv(2):iv(3),:) = S(:,1,:)
  if (nu2 .gt. 1) then
     Q(iu(2):iu(3),iv(2),:) = W(:,1,:)
     Q(iu(2):iu(3),iv(3),:) = E(:,1,:)
  end if
  dQdw(:,:,:) = 0.0

  verts(:,:,:,:) = 0.0
  do i=1,4
     do j=1,4
        verts(i,j,1,:) = Q(iu(i),iv(j),:)
     end do
  end do
  j = 1
  do i=2,3
     verts(i,j,2,:) = f0*(Q(iu(i),iv(j),:) - Q(iu(i),iv(j)-1,:))
  end do
  j = 4
  do i=2,3
     verts(i,j,2,:) = f0*(Q(iu(i),iv(j),:) - Q(iu(i),iv(j)+1,:))
  end do
  i = 1
  do j=2,3
     verts(i,j,2,:) = f0*(Q(iu(i),iv(j),:) - Q(iu(i)-1,iv(j),:))
  end do
  i = 4
  do j=2,3
     verts(i,j,2,:) = f0*(Q(iu(i),iv(j),:) - Q(iu(i)+1,iv(j),:))
  end do
  verts(2,2,2,:) = m0*(N(1,1,:) - N(1,2,:))
  verts(2,3,2,:) = m0*(N(nv2,1,:) - N(nv2,2,:))
  verts(3,2,2,:) = m0*(S(1,1,:) - S(1,2,:))
  verts(3,3,2,:) = m0*(S(nv2,1,:) - S(nv2,2,:))

  i = 1
  do j=2,3
     call bezierCurve(nu1-1, verts(i,j,1,:), verts(i,j,2,:), &
          verts(i+1,j,1,:), verts(i+1,j,2,:), &
          Q(iu(i):iu(i+1),iv(j),:))
  end do
  i = 3
  do j=2,3
     call bezierCurve(nu3-1, verts(i,j,1,:), verts(i,j,2,:), &
          verts(i+1,j,1,:), verts(i+1,j,2,:), &
          Q(iu(i):iu(i+1),iv(j),:))
  end do
  j = 1
  do i=2,3
     call bezierCurve(nv1-1, verts(i,j,1,:), verts(i,j,2,:), &
          verts(i,j+1,1,:), verts(i,j+1,2,:), &
          Q(iu(i),iv(j):iv(j+1),:))
  end do
  j = 3
  do i=2,3
     call bezierCurve(nv3-1, verts(i,j,1,:), verts(i,j,2,:), &
          verts(i,j+1,1,:), verts(i,j+1,2,:), &
          Q(iu(i),iv(j):iv(j+1),:))
  end do

  call interpolateFrames(1, nu, nv, iu, iv, Q, dQdw)
  call interpolateFrames(2, nu, nv, iu, iv, Q, dQdw)
  call interpolateFrames(3, nu, nv, iu, iv, Q, dQdw)
  if (nu2 .eq. 1) then
     Q(iu(2),iv(2):iv(3),:) = N(:,1,:)
     call interpolateFrames(1, nu, nv, iu, iv, Q, dQdw)
  end if

  do k=1,3
     Q(:,:,k) = Q(:,:,k) + shape0(:,:)*dQdw(:,:,k)
  end do

end subroutine computeJunction
