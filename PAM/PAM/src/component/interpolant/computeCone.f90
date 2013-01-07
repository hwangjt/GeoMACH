subroutine computeCone(nu, nv, nu1, nu2, nu3, nv1, nv2, nv3, &
     scale0, f0, m0, W, E, N, S, shape0, Q)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nu, nv, nu1, nu2, nu3, nv1, nv2, nv3, scale0, f0, m0, W, E, N, S, shape0
  !f2py intent(out) Q
  !f2py depend(nu) W, E
  !f2py depend(nv) N, S
  !f2py depend(nu,nv) shape0, Q

  !Input
  integer, intent(in) ::  nu, nv, nu1, nu2, nu3, nv1, nv2, nv3
  double precision, intent(in) ::  scale0, f0, m0
  double precision, intent(in) ::  W(nu,2,3), E(nu,2,3) 
  double precision, intent(in) ::  N(nv,2,3), S(nv,2,3) 
  double precision, intent(in) ::  shape0(nu,nv)

  !Output
  double precision, intent(out) ::  Q(nu,nv,3)

  !Working
  integer k, iu(4), iv(4)
  double precision C(3), C1(3), C2(3), H(3), V(3)
  double precision tW(3), tE(3), tN(3), tS(3)
  double precision dQdw(nu,nv,3)

  iu(1) = 1
  iu(2) = nu1 + 1
  iu(3) = nu1 + nu2 + 1
  iu(4) = nu1 + nu2 + nu3 + 1

  iv(1) = 1
  iv(2) = nv1 + 1
  iv(3) = nv1 + nv2 + 1
  iv(4) = nv1 + nv2 + nv3 + 1

  C1 = 0.5*(W(iu(2),1,:) + E(iu(2),1,:))
  C2 = 0.5*(W(iu(2),2,:) + E(iu(2),2,:))
  C = C1 + (C1-C2)*scale0
  H = 0.5*m0*(E(iu(2),1,:) - W(iu(2),1,:))
  if (nu3 .gt. 0) then
     V = 0.5*m0*(S(iv(2),1,:) - N(iv(2),1,:))
  else
     V = m0*(C1 - N(iv(2),1,:))
  end if

  tW = scale0*f0*(W(iu(2),1,:) - W(iu(2),2,:))
  tE = scale0*f0*(E(iu(2),2,:) - E(iu(2),1,:))
  tN = scale0*f0*(N(iv(2),1,:) - N(iv(2),2,:))
  tS = scale0*f0*(S(iv(2),2,:) - S(iv(2),1,:))

  Q(:,:,:) = 0.0
  Q(:,iv(1),:) = W(:,1,:)
  Q(:,iv(4),:) = E(:,1,:)
  Q(iu(1),:,:) = N(:,1,:)
  if (nu3 .gt. 0) then
     Q(iu(4),:,:) = S(:,1,:)
  end if

  call bezierCurve(nv1+1, W(iu(2),1,:), tW, C, -H, Q(iu(2),iv(1):iv(2),:))
  call bezierCurve(nv3+1, C, H, E(iu(2),1,:), -tE, Q(iu(2),iv(3):iv(4),:))
  call bezierCurve(nu1+1, N(iv(2),1,:), tN, C, -V, Q(iu(1):iu(2),iv(2),:))
  if (nu3 .gt. 0) then
     call bezierCurve(nu3+1, C, V, S(iv(2),1,:), -tS, Q(iu(3):iu(4),iv(2),:))
  end if

  Q(iu(3),:,:) = Q(iu(2),:,:)
  Q(:,iv(3),:) = Q(:,iv(2),:)

  call interpolateFrames(nu, nv, iu, iv, Q, dQdw)

  do k=1,3
     Q(:,:,k) = Q(:,:,k) + shape0(:,:)*dQdw(:,:,k)
  end do

end subroutine computeCone
