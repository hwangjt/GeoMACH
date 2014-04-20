subroutine computeCoons(nD, nu, nv, nu1, nu2, nu3, nv1, nv2, nv3, &
     inds, Da, Di, Dj)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nD, nu, nv, nu1, nu2, nu3, nv1, nv2, nv3, inds
  !f2py intent(out) Da, Di, Dj
  !f2py depend(nu, nv) inds
  !f2py depend(nD) Da, Di, Dj

  !Input
  integer, intent(in) ::  nD, nu, nv, nu1, nu2, nu3, nv1, nv2, nv3
  integer, intent(in) ::  inds(nu,nv)

  !Output
  double precision, intent(out) ::  Da(nD)
  integer, intent(out) ::  Di(nD), Dj(nD)

  !Working
  integer i, j, iu(4), iv(4), iD, ii, jj
  double precision denu, denv, u, v

  iu(1) = 1
  iu(2) = (nu1-1) + 1
  iu(3) = (nu1-1) + (nu2-1) + 1
  iu(4) = (nu1-1) + (nu2-1) + (nu3-1) + 1

  iv(1) = 1
  iv(2) = (nv1-1) + 1
  iv(3) = (nv1-1) + (nv2-1) + 1
  iv(4) = (nv1-1) + (nv2-1) + (nv3-1) + 1
  
  iD = 0

  do jj=1,3
     denv = 1.0 / (iv(jj+1)-iv(jj))
     do ii=1,3
        denu = 1.0 / (iu(ii+1)-iu(ii))
        if ((nu2 .ne. 1) .or. (ii .ne. 2)) then
           do j=iv(jj)+1,iv(jj+1)-1
              v = (j - iv(jj)) * denv
              do i=iu(ii)+1,iu(ii+1)-1
                 u = (i - iu(ii)) * denu

                 Di(iD+1:iD+8) = inds(i, j)
                 Da(iD+1) = -(1-u) * (1-v)
                 Dj(iD+1) = inds(iu(ii), iv(jj))
                 Da(iD+2) = -u * (1-v)
                 Dj(iD+2) = inds(iu(ii+1), iv(jj))
                 Da(iD+3) = -(1-u) * v
                 Dj(iD+3) = inds(iu(ii), iv(jj+1))
                 Da(iD+4) = -u * v
                 Dj(iD+4) = inds(iu(ii+1), iv(jj+1))
                 Da(iD+5) = 1-u
                 Dj(iD+5) = inds(iu(ii), j)
                 Da(iD+6) = u
                 Dj(iD+6) = inds(iu(ii+1), j)
                 Da(iD+7) = 1-v
                 Dj(iD+7) = inds(i, iv(jj))
                 Da(iD+8) = v
                 Dj(iD+8) = inds(i, iv(jj+1))
                 iD = iD + 8
              end do
           end do
        end if
     end do     
  end do

  if (iD .ne. nD) then
     print *, 'Error in computeCoons', iD, nD
  end if

end subroutine computeCoons



subroutine computeWireframe(nD, nu, nv, nu1, nu2, nu3, nv1, nv2, nv3, &
     f0, m0, W, E, N, S, fInds, inds, Da, Di, Dj)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nD, nu, nv, nu1, nu2, nu3, nv1, nv2, nv3, f0, m0, W, E, N, S, fInds, inds
  !f2py intent(out) Da, Di, Dj
  !f2py depend(nu2) W, E
  !f2py depend(nv2) N, S
  !f2py depend(nu,nv) fInds, inds
  !f2py depend(nD) Da, Di, Dj

  !Input
  integer, intent(in) ::  nD, nu, nv, nu1, nu2, nu3, nv1, nv2, nv3
  double precision, intent(in) ::  f0, m0
  integer, intent(in) ::  W(nu2,2), E(nu2,2)
  integer, intent(in) ::  N(nv2,2), S(nv2,2)
  integer, intent(in) ::  fInds(nu,nv), inds(nu,nv)

  !Output
  double precision, intent(out) ::  Da(nD)
  integer, intent(out) ::  Di(nD), Dj(nD)

  !Working
  integer i, j, iu(4), iv(4), iD
  double precision den, C(4)

  iu(1) = 1
  iu(2) = (nu1-1) + 1
  iu(3) = (nu1-1) + (nu2-1) + 1
  iu(4) = (nu1-1) + (nu2-1) + (nu3-1) + 1

  iv(1) = 1
  iv(2) = (nv1-1) + 1
  iv(3) = (nv1-1) + (nv2-1) + 1
  iv(4) = (nv1-1) + (nv2-1) + (nv3-1) + 1
  
  iD = 0

  ! Border: N, S
  do i=1,nu,nu-1
     do j=1,nv
        iD = iD + 1
        Da(iD) = 1.0
        Di(iD) = inds(i, j)
        Dj(iD) = fInds(i, j)
     end do
  end do

  ! Border: W, E
  do j=1,nv,nv-1
     do i=2,nu-1
        iD = iD + 1
        Da(iD) = 1.0
        Di(iD) = inds(i, j)
        Dj(iD) = fInds(i, j)
     end do
  end do

  ! Upper middle 2 corners
  Da(iD+1:iD+2) = 1.0
  Di(iD+1) = inds(iu(2), iv(2))
  Dj(iD+1) = N(1,1)
  Di(iD+2) = inds(iu(2), iv(3))
  Dj(iD+2) = N(nv2,1)
  iD = iD + 2

  ! Lower middle 2 corners
  if (nu2 .ne. 1) then
     Da(iD+1:iD+2) = 1.0
     Di(iD+1) = inds(iu(3), iv(2))
     Dj(iD+1) = S(1,1)
     Di(iD+2) = inds(iu(3), iv(3))
     Dj(iD+2) = S(nv2,1)
     iD = iD + 2
  end if

  ! North
  i = iu(2)
  do j=iv(2)+1,iv(3)-1
     iD = iD + 1
     Da(iD) = 1.0
     Di(iD) = inds(i, j)
     Dj(iD) = N(j,1)
  end do

  ! South
  i = iu(3)
  do j=iv(2)+1,iv(3)-1
     iD = iD + 1
     Da(iD) = 1.0
     Di(iD) = inds(i, j)
     Dj(iD) = S(j,1)
  end do

  if (nu2 .ne. 1) then  
     ! West
     j = iv(2)
     do i=iu(2)+1,iu(3)-1
        iD = iD + 1
        Da(iD) = 1.0
        Di(iD) = inds(i, j)
        Dj(iD) = W(i,1)
     end do

     ! East
     j = iv(3)
     do i=iu(2)+1,iu(3)-1
        iD = iD + 1
        Da(iD) = 1.0
        Di(iD) = inds(i, j)
        Dj(iD) = E(i,1)
     end do
  end if

  j = iv(2)
  den = 1.0 / (iu(2)-iu(1))
  do i=iu(1)+1,iu(2)-1
     call sparseBezier(den * (i-iu(1)), f0, m0, C)
     Da(iD+1:iD+4) = C(:)
     Di(iD+1:iD+4) = inds(i, j)
     Dj(iD+1) = fInds(iu(1), j)
     Dj(iD+2) = fInds(iu(1)+1, j)
     Dj(iD+3) = N(1,1)
     Dj(iD+4) = N(1,2)
     iD = iD + 4
  end do

  j = iv(3)
  den = 1.0 / (iu(2)-iu(1))
  do i=iu(1)+1,iu(2)-1
     call sparseBezier(den * (i-iu(1)), f0, m0, C)
     Da(iD+1:iD+4) = C(:)
     Di(iD+1:iD+4) = inds(i, j)
     Dj(iD+1) = fInds(iu(1), j)
     Dj(iD+2) = fInds(iu(1)+1, j)
     Dj(iD+3) = N(nv2,1)
     Dj(iD+4) = N(nv2,2)
     iD = iD + 4
  end do

  j = iv(2)
  den = 1.0 / (iu(4)-iu(3))
  do i=iu(3)+1,iu(4)-1
     call sparseBezier(den * (iu(4)-i), f0, m0, C)
     Da(iD+1:iD+4) = C(:)
     Di(iD+1:iD+4) = inds(i, j)
     Dj(iD+1) = fInds(iu(4), j)
     Dj(iD+2) = fInds(iu(4)-1, j)
     Dj(iD+3) = S(1,1)
     Dj(iD+4) = S(1,2)
     iD = iD + 4
  end do

  j = iv(3)
  den = 1.0 / (iu(4)-iu(3))
  do i=iu(3)+1,iu(4)-1
     call sparseBezier(den * (iu(4)-i), f0, m0, C)
     Da(iD+1:iD+4) = C(:)
     Di(iD+1:iD+4) = inds(i, j)
     Dj(iD+1) = fInds(iu(4), j)
     Dj(iD+2) = fInds(iu(4)-1, j)
     Dj(iD+3) = S(nv2,1)
     Dj(iD+4) = S(nv2,2)
     iD = iD + 4
  end do

  i = iu(2)
  den = 1.0 / (iv(2)-iv(1))
  do j=iv(1)+1,iv(2)-1
     call sparseBezier(den * (j-iv(1)), f0, m0, C)
     Da(iD+1:iD+4) = C(:)
     Di(iD+1:iD+4) = inds(i, j)
     Dj(iD+1) = fInds(i, iv(1))
     Dj(iD+2) = fInds(i, iv(1)+1)
     Dj(iD+3) = N(1,1)
     Dj(iD+4) = N(1,2)
     iD = iD + 4
  end do

  i = iu(2)
  den = 1.0 / (iv(4)-iv(3))
  do j=iv(3)+1,iv(4)-1
     call sparseBezier(den * (iv(4)-j), f0, m0, C)
     Da(iD+1:iD+4) = C(:)
     Di(iD+1:iD+4) = inds(i, j)
     Dj(iD+1) = fInds(i, iv(4))
     Dj(iD+2) = fInds(i, iv(4)-1)
     Dj(iD+3) = N(nv2,1)
     Dj(iD+4) = N(nv2,2)
     iD = iD + 4
  end do

  if (nu2 .ne. 1) then
     i = iu(3)
     den = 1.0 / (iv(2)-iv(1))
     do j=iv(1)+1,iv(2)-1
        call sparseBezier(den * (j-iv(1)), f0, m0, C)
        Da(iD+1:iD+4) = C(:)
        Di(iD+1:iD+4) = inds(i, j)
        Dj(iD+1) = fInds(i, iv(1))
        Dj(iD+2) = fInds(i, iv(1)+1)
        Dj(iD+3) = S(1,1)
        Dj(iD+4) = S(1,2)
        iD = iD + 4
     end do
     
     i = iu(3)
     den = 1.0 / (iv(4)-iv(3))
     do j=iv(3)+1,iv(4)-1
        call sparseBezier(den * (iv(4)-j), f0, m0, C)
        Da(iD+1:iD+4) = C(:)
        Di(iD+1:iD+4) = inds(i, j)
        Dj(iD+1) = fInds(i, iv(4))
        Dj(iD+2) = fInds(i, iv(4)-1)
        Dj(iD+3) = S(nv2,1)
        Dj(iD+4) = S(nv2,2)
        iD = iD + 4
     end do
  end if

  if (iD .ne. nD) then
     print *, 'Error in computeWireframe', iD, nD
  end if

end subroutine computeWireframe



subroutine sparseBezier(u, s0, s1, C)

  implicit none

  !Input
  double precision, intent(in) ::  u, s0, s1

  !Output
  double precision, intent(out) ::  C(4)

  !Working
  double precision u1, u2, u3, t1, t2, t3

  u1 = u
  u2 = u*u1
  u3 = u*u2
  t1 = (1-u)
  t2 = (1-u)*t1
  t3 = (1-u)*t2

  C(1) = t3 + (3-s0) * u1*t2
  C(2) = s0 * u1*t2
  C(3) = u3 + (s1+3) * u2*t1
  C(4) = -s1 * u2*t1

end subroutine sparseBezier



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
