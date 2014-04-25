subroutine computeConeCoons(nD, nu, nv, inds, Da, Di, Dj)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nD, nu, nv, inds
  !f2py intent(out) Da, Di, Dj
  !f2py depend(nu,nv) inds
  !f2py depend(nD) Da, Di, Dj

  !Input
  integer, intent(in) ::  nD, nu, nv
  integer, intent(in) ::  inds(nu,nv)

  !Output
  double precision, intent(out) ::  Da(nD)
  integer, intent(out) ::  Di(nD), Dj(nD)

  !Working
  integer i, j, iu(4), iv(4), iD, ii, jj
  double precision denu, denv, u, v

  iu(1) = 1
  iu(2) = int((nu-1)/2) + 1
  iu(3) = nu

  iv(1) = 1
  iv(2) = int((nv-1)/2) + 1
  iv(3) = nv
  
  iD = 0

  do jj=1,2
     denv = 1.0 / (iv(jj+1)-iv(jj))
     do ii=1,2
        denu = 1.0 / (iu(ii+1)-iu(ii))
        do j=iv(jj)+1,iv(jj+1)-1
           v = (j - iv(jj)) * denv
           do i=iu(ii)+1,iu(ii+1)-1
              u = (i - iu(ii)) * denu
              call computeCoons(u, v, inds(i, j), &
                   inds(iu(ii), j), inds(iu(ii+1), j), &
                   inds(i, iv(jj)), inds(i, iv(jj+1)), &
                   inds(iu(ii), iv(jj)), inds(iu(ii+1), iv(jj)), &
                   inds(iu(ii), iv(jj+1)), inds(iu(ii+1), iv(jj+1)), &
                   Da(iD+1:iD+8), Di(iD+1:iD+8), Dj(iD+1:iD+8))
              iD = iD + 8
           end do
        end do
     end do
  end do
    
  if (iD .ne. nD) then
     print *, 'Error in computeConeCoons', iD, nD
  end if    

end subroutine computeConeCoons



subroutine computeConeWireframe(nD, nu, nv, scale0, f0, m0, W, E, N, S, inds, Da, Di, Dj)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nD, nu, nv, scale0, f0, m0, W, E, N, S, inds
  !f2py intent(out) Da, Di, Dj
  !f2py depend(nu) W, E
  !f2py depend(nv) N, S
  !f2py depend(nu,nv), inds
  !f2py depend(nD) Da, Di, Dj

  !Input
  integer, intent(in) ::  nD, nu, nv
  double precision, intent(in) ::  scale0, f0, m0
  integer, intent(in) ::  W(nu,2), E(nu,2)
  integer, intent(in) ::  N(nv,2), S(nv,2)
  integer, intent(in) ::  inds(nu,nv)

  !Output
  double precision, intent(out) ::  Da(nD)
  integer, intent(out) ::  Di(nD), Dj(nD)

  !Working
  integer i, j, iD, nuC, nvC
  double precision den, C(4), u, v

  nuC = int((nu-1)/2) + 1
  nvC = int((nv-1)/2) + 1

  Da(:) = 0.0
  Di(:) = 0
  Dj(:) = 0
  
  iD = 0

  ! Border: N
  i = 1
  do j=1,nv
     iD = iD + 1
     Da(iD) = 1.0
     Di(iD) = inds(i, j)
     Dj(iD) = N(j, 1)
  end do

  ! Border: S
  i = nu
  do j=1,nv
     iD = iD + 1
     Da(iD) = 1.0
     Di(iD) = inds(i, j)
     Dj(iD) = S(j, 1)
  end do

  ! Border: W
  j = 1
  do i=2,nu-1
     iD = iD + 1
     Da(iD) = 1.0
     Di(iD) = inds(i, j)
     Dj(iD) = W(i, 1)
  end do

  ! Border: E
  j = nv
  do i=2,nu-1
     iD = iD + 1
     Da(iD) = 1.0
     Di(iD) = inds(i, j)
     Dj(iD) = E(i, 1)
  end do

  call sparseBezier(dble(0.5), -f0, f0, C)

  i = nuC
  j = nvC
  Da(iD+1:iD+4) = 0.5 * C(:)
  Da(iD+5:iD+8) = 0.5 * C(:)
  Di(iD+1:iD+8) = inds(i, j)
  Dj(iD+1) = N(j,1)
  Dj(iD+2) = N(j,2)
  Dj(iD+3) = S(j,1)
  Dj(iD+4) = S(j,2)
  Dj(iD+5) = W(i,1)
  Dj(iD+6) = W(i,2)
  Dj(iD+7) = E(i,1)
  Dj(iD+8) = E(i,2)
  iD = iD + 8

  i = nuC
  den = 1.0 / (nvC-1)
  do j=2,nvC-1
     v = den * (j-1)
     call sparseBezier(v, -f0, f0, C)
     Da(iD+1:iD+4) = (1-v) * C(:)
     Di(iD+1:iD+4) = inds(i, j)
     Dj(iD+1) = W(i,1)
     Dj(iD+2) = W(i,2)
     Dj(iD+3) = E(i,1)
     Dj(iD+4) = E(i,2)
     iD = iD + 4
     call sparseBezier(dble(0.5), -f0, f0, C)
     Da(iD+1:iD+4) = v * C(:)
     Di(iD+1:iD+4) = inds(i, j)
     Dj(iD+1) = N(j,1)
     Dj(iD+2) = N(j,2)
     Dj(iD+3) = S(j,1)
     Dj(iD+4) = S(j,2)
     iD = iD + 4
  end do
  do j=nvC+1,nv-1
     v = den * (j-nvC)
     call sparseBezier(v, -f0, f0, C)
     Da(iD+1:iD+4) = v * C(:)
     Di(iD+1:iD+4) = inds(i, j)
     Dj(iD+1) = W(i,1)
     Dj(iD+2) = W(i,2)
     Dj(iD+3) = E(i,1)
     Dj(iD+4) = E(i,2)
     iD = iD + 4
     call sparseBezier(dble(0.5), -f0, f0, C)
     Da(iD+1:iD+4) = (1-v) * C(:)
     Di(iD+1:iD+4) = inds(i, j)
     Dj(iD+1) = N(j,1)
     Dj(iD+2) = N(j,2)
     Dj(iD+3) = S(j,1)
     Dj(iD+4) = S(j,2)
     iD = iD + 4
  end do

  j = nvC
  den = 1.0 / (nuC-1)
  do i=2,nuC-1
     u = den * (i-1)
     call sparseBezier(u, -f0, f0, C)
     Da(iD+1:iD+4) = (1-u) * C(:)
     Di(iD+1:iD+4) = inds(i, j)
     Dj(iD+1) = N(j,1)
     Dj(iD+2) = N(j,2)
     Dj(iD+3) = S(j,1)
     Dj(iD+4) = S(j,2)
     iD = iD + 4
     call sparseBezier(dble(0.5), -f0, f0, C)
     Da(iD+1:iD+4) = u * C(:)
     Di(iD+1:iD+4) = inds(i, j)
     Dj(iD+1) = W(i,1)
     Dj(iD+2) = W(i,2)
     Dj(iD+3) = E(i,1)
     Dj(iD+4) = E(i,2)
     iD = iD + 4
  end do
  do i=nuC+1,nu-1
     u = den * (i-nuC)
     call sparseBezier(u, -f0, f0, C)
     Da(iD+1:iD+4) = u * C(:)
     Di(iD+1:iD+4) = inds(i, j)
     Dj(iD+1) = N(j,1)
     Dj(iD+2) = N(j,2)
     Dj(iD+3) = S(j,1)
     Dj(iD+4) = S(j,2)
     iD = iD + 4
     call sparseBezier(dble(0.5), -f0, f0, C)
     Da(iD+1:iD+4) = (1-u) * C(:)
     Di(iD+1:iD+4) = inds(i, j)
     Dj(iD+1) = W(i,1)
     Dj(iD+2) = W(i,2)
     Dj(iD+3) = E(i,1)
     Dj(iD+4) = E(i,2)
     iD = iD + 4
  end do

  if (iD .ne. nD) then
     print *, 'Error in computeConeWireframe', iD, nD
  end if

end subroutine computeConeWireframe
