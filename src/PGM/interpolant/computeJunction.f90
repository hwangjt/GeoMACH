subroutine computeJunctionCoons(nD, nu, nv, nu1, nu2, nu3, nv1, nv2, nv3, &
     inds, Da, Di, Dj)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nD, nu, nv, nu1, nu2, nu3, nv1, nv2, nv3, inds
  !f2py intent(out) Da, Di, Dj
  !f2py depend(nu, nv) inds
  !f2py depend(nD) Da, Di, Dj

  !Input
  integer, intent(in) ::  nD, nu, nv, nu1, nu2, nu3, nv1, nv2, nv3
  integer, intent(in) ::  inds(nu,nv,3)

  !Output
  double precision, intent(out) ::  Da(nD)
  integer, intent(out) ::  Di(nD), Dj(nD)

  !Working
  integer i, j, k, iu(4), iv(4), iD, ii, jj
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
                 do k=1,3
                    call computeCoons(u, v, inds(i, j, k), &
                         inds(iu(ii), j, k), inds(iu(ii+1), j, k), &
                         inds(i, iv(jj), k), inds(i, iv(jj+1), k), &
                         inds(iu(ii), iv(jj), k), inds(iu(ii+1), iv(jj), k), &
                         inds(iu(ii), iv(jj+1), k), inds(iu(ii+1), iv(jj+1), k), &
                         Da(iD+1:iD+8), Di(iD+1:iD+8), Dj(iD+1:iD+8))
                    iD = iD + 8
                 end do
              end do
           end do
        end if
     end do     
  end do

  if (iD .ne. nD) then
     print *, 'Error in computeJunctionCoons', iD, nD
  end if

end subroutine computeJunctionCoons



subroutine computeJunctionWireframe(nD, nu, nv, nu1, nu2, nu3, nv1, nv2, nv3, &
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
  integer, intent(in) ::  W(nu2,2,3), E(nu2,2,3)
  integer, intent(in) ::  N(nv2,2,3), S(nv2,2,3)
  integer, intent(in) ::  fInds(nu,nv,3), inds(nu,nv,3)

  !Output
  double precision, intent(out) ::  Da(nD)
  integer, intent(out) ::  Di(nD), Dj(nD)

  !Working
  integer i, j, k, iu(4), iv(4), iD
  double precision den, C(4)

  Da(:) = 0.0
  Di(:) = 0
  Dj(:) = 0

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
        do k=1,3
           iD = iD + 1
           Da(iD) = 1.0
           Di(iD) = inds(i, j, k)
           Dj(iD) = fInds(i, j, k)
        end do
     end do
  end do

  ! Border: W, E
  do j=1,nv,nv-1
     do i=2,nu-1
        do k=1,3
           iD = iD + 1
           Da(iD) = 1.0
           Di(iD) = inds(i, j, k)
           Dj(iD) = fInds(i, j, k)
        end do
     end do
  end do

  ! Upper middle 2 corners
  do k=1,3
     Da(iD+1:iD+2) = 1.0
     Di(iD+1) = inds(iu(2), iv(2), k)
     Dj(iD+1) = N(1,1,k)
     Di(iD+2) = inds(iu(2), iv(3), k)
     Dj(iD+2) = N(nv2,1,k)
     iD = iD + 2
  end do

  ! Lower middle 2 corners
  if (nu2 .ne. 1) then
     do k=1,3
        Da(iD+1:iD+2) = 1.0
        Di(iD+1) = inds(iu(3), iv(2), k)
        Dj(iD+1) = S(1,1,k)
        Di(iD+2) = inds(iu(3), iv(3), k)
        Dj(iD+2) = S(nv2,1,k)
        iD = iD + 2
     end do
  end if

  ! North
  i = iu(2)
  do j=iv(2)+1,iv(3)-1
     do k=1,3
        iD = iD + 1
        Da(iD) = 1.0
        Di(iD) = inds(i, j, k)
        Dj(iD) = N(j-iv(2)+1,1,k)
     end do
  end do

  ! South
  i = iu(3)
  do j=iv(2)+1,iv(3)-1
     do k=1,3
        iD = iD + 1
        Da(iD) = 1.0
        Di(iD) = inds(i, j, k)
        Dj(iD) = S(j-iv(2)+1,1,k)
     end do
  end do

  if (nu2 .ne. 1) then  
     ! West
     j = iv(2)
     do i=iu(2)+1,iu(3)-1
        do k=1,3
           iD = iD + 1
           Da(iD) = 1.0
           Di(iD) = inds(i, j, k)
           Dj(iD) = W(i-iu(2)+1,1,k)
        end do
     end do

     ! East
     j = iv(3)
     do i=iu(2)+1,iu(3)-1
        do k=1,3
           iD = iD + 1
           Da(iD) = 1.0
           Di(iD) = inds(i, j, k)
           Dj(iD) = E(i-iu(2)+1,1,k)
        end do
     end do
  end if

  j = iv(2)
  den = 1.0 / (iu(2)-iu(1))
  do i=iu(1)+1,iu(2)-1
     do k=1,3
        call sparseBezier(den * (i-iu(1)), f0, m0, C)
        Da(iD+1:iD+4) = C(:)
        Di(iD+1:iD+4) = inds(i, j, k)
        Dj(iD+1) = fInds(iu(1), j, k)
        Dj(iD+2) = fInds(iu(1)+1, j, k)
        Dj(iD+3) = N(1,1,k)
        Dj(iD+4) = N(1,2,k)
        iD = iD + 4
     end do
  end do

  j = iv(3)
  den = 1.0 / (iu(2)-iu(1))
  do i=iu(1)+1,iu(2)-1
     do k=1,3
        call sparseBezier(den * (i-iu(1)), f0, m0, C)
        Da(iD+1:iD+4) = C(:)
        Di(iD+1:iD+4) = inds(i, j, k)
        Dj(iD+1) = fInds(iu(1), j, k)
        Dj(iD+2) = fInds(iu(1)+1, j, k)
        Dj(iD+3) = N(nv2,1,k)
        Dj(iD+4) = N(nv2,2,k)
        iD = iD + 4
     end do
  end do

  j = iv(2)
  den = 1.0 / (iu(4)-iu(3))
  do i=iu(3)+1,iu(4)-1
     do k=1,3
        call sparseBezier(den * (iu(4)-i), f0, m0, C)
        Da(iD+1:iD+4) = C(:)
        Di(iD+1:iD+4) = inds(i, j, k)
        Dj(iD+1) = fInds(iu(4), j, k)
        Dj(iD+2) = fInds(iu(4)-1, j, k)
        Dj(iD+3) = S(1,1,k)
        Dj(iD+4) = S(1,2,k)
        iD = iD + 4
     end do
  end do

  j = iv(3)
  den = 1.0 / (iu(4)-iu(3))
  do i=iu(3)+1,iu(4)-1
     do k=1,3
        call sparseBezier(den * (iu(4)-i), f0, m0, C)
        Da(iD+1:iD+4) = C(:)
        Di(iD+1:iD+4) = inds(i, j, k)
        Dj(iD+1) = fInds(iu(4), j, k)
        Dj(iD+2) = fInds(iu(4)-1, j, k)
        Dj(iD+3) = S(nv2,1,k)
        Dj(iD+4) = S(nv2,2,k)
        iD = iD + 4
     end do
  end do

  i = iu(2)
  den = 1.0 / (iv(2)-iv(1))
  do j=iv(1)+1,iv(2)-1
     do k=1,3
        call sparseBezier(den * (j-iv(1)), f0, m0, C)
        Da(iD+1:iD+4) = C(:)
        Di(iD+1:iD+4) = inds(i, j, k)
        Dj(iD+1) = fInds(i, iv(1), k)
        Dj(iD+2) = fInds(i, iv(1)+1, k)
        Dj(iD+3) = N(1,1,k)
        Dj(iD+4) = N(1,2,k)
        iD = iD + 4
     end do
  end do

  i = iu(2)
  den = 1.0 / (iv(4)-iv(3))
  do j=iv(3)+1,iv(4)-1
     do k=1,3
        call sparseBezier(den * (iv(4)-j), f0, m0, C)
        Da(iD+1:iD+4) = C(:)
        Di(iD+1:iD+4) = inds(i, j, k)
        Dj(iD+1) = fInds(i, iv(4), k)
        Dj(iD+2) = fInds(i, iv(4)-1, k)
        Dj(iD+3) = N(nv2,1,k)
        Dj(iD+4) = N(nv2,2,k)
        iD = iD + 4
     end do
  end do

  if (nu2 .ne. 1) then
     i = iu(3)
     den = 1.0 / (iv(2)-iv(1))
     do j=iv(1)+1,iv(2)-1
        do k=1,3
           call sparseBezier(den * (j-iv(1)), f0, m0, C)
           Da(iD+1:iD+4) = C(:)
           Di(iD+1:iD+4) = inds(i, j, k)
           Dj(iD+1) = fInds(i, iv(1), k)
           Dj(iD+2) = fInds(i, iv(1)+1, k)
           Dj(iD+3) = S(1,1,k)
           Dj(iD+4) = S(1,2,k)
           iD = iD + 4
        end do
     end do
     
     i = iu(3)
     den = 1.0 / (iv(4)-iv(3))
     do j=iv(3)+1,iv(4)-1
        do k=1,3
           call sparseBezier(den * (iv(4)-j), f0, m0, C)
           Da(iD+1:iD+4) = C(:)
           Di(iD+1:iD+4) = inds(i, j, k)
           Dj(iD+1) = fInds(i, iv(4), k)
           Dj(iD+2) = fInds(i, iv(4)-1, k)
           Dj(iD+3) = S(nv2,1,k)
           Dj(iD+4) = S(nv2,2,k)
           iD = iD + 4
        end do
     end do
  end if

  if (iD .ne. nD) then
     print *, 'Error in computeJunctionWireframe', iD, nD
  end if

end subroutine computeJunctionWireframe
