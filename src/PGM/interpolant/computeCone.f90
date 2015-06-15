subroutine computeConeCoons(nD, nu, nv, inds, Da, Di, Dj)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nD, nu, nv, inds
  !f2py intent(out) Da, Di, Dj
  !f2py depend(nu,nv) inds
  !f2py depend(nD) Da, Di, Dj

  !Input
  integer, intent(in) ::  nD, nu, nv
  integer, intent(in) ::  inds(nu,nv,3)

  !Output
  double precision, intent(out) ::  Da(nD)
  integer, intent(out) ::  Di(nD), Dj(nD)

  !Working
  integer i, j, k, iu(4), iv(4), iD, ii, jj
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
     end do
  end do
    
  if (iD .ne. nD) then
     print *, 'Error in computeConeCoons', iD, nD
  end if    

end subroutine computeConeCoons



subroutine computeConeWireframe(nD, nu, nv, a, f0, W, E, N, S, inds, Da, Di, Dj)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nD, nu, nv, a, f0, W, E, N, S, inds
  !f2py intent(out) Da, Di, Dj
  !f2py depend(nu) W, E
  !f2py depend(nv) N, S
  !f2py depend(nu,nv), inds
  !f2py depend(nD) Da, Di, Dj

  !Input
  integer, intent(in) ::  nD, nu, nv
  double precision, intent(in) ::  a, f0
  integer, intent(in) ::  W(nu,2,3), E(nu,2,3)
  integer, intent(in) ::  N(nv,2,3), S(nv,2,3)
  integer, intent(in) ::  inds(nu,nv,3)

  !Output
  double precision, intent(out) ::  Da(nD)
  integer, intent(out) ::  Di(nD), Dj(nD)

  !Working
  integer i, j, k, iD, nuC, nvC
  double precision den, C(4), u, v, t, pi

  pi = 2*acos(0.0)

  nuC = int((nu-1)/2) + 1
  nvC = int((nv-1)/2) + 1

  Da(:) = 0.0
  Di(:) = 0
  Dj(:) = 0
  
  iD = 0

  ! Border: N
  i = 1
  do j=1,nv
     do k=1,3
        iD = iD + 1
        Da(iD) = 1.0
        Di(iD) = inds(i, j, k)
        Dj(iD) = N(j, 1, k)
     end do
  end do

  ! Border: S
  i = nu
  do j=1,nv
     do k=1,3
        iD = iD + 1
        Da(iD) = 1.0
        Di(iD) = inds(i, j, k)
        Dj(iD) = S(j, 1, k)
     end do
  end do

  ! Border: W
  j = 1
  do i=2,nu-1
     do k=1,3
        iD = iD + 1
        Da(iD) = 1.0
        Di(iD) = inds(i, j, k)
        Dj(iD) = W(i, 1, k)
     end do
  end do

  ! Border: E
  j = nv
  do i=2,nu-1
     do k=1,3
        iD = iD + 1
        Da(iD) = 1.0
        Di(iD) = inds(i, j, k)
        Dj(iD) = E(i, 1, k)
     end do
  end do

  ! Center point
  call sparseBezier(dble(0.5), -f0, f0, C)
  i = nuC
  j = nvC
  do k=1,3
     Da(iD+1:iD+4) = 0.5 * C(:)
     Da(iD+5:iD+8) = 0.5 * C(:)
     Di(iD+1:iD+8) = inds(i, j, k)
     Dj(iD+1) = N(j,1,k)
     Dj(iD+2) = N(j,2,k)
     Dj(iD+3) = S(j,1,k)
     Dj(iD+4) = S(j,2,k)
     Dj(iD+5) = W(i,1,k)
     Dj(iD+6) = W(i,2,k)
     Dj(iD+7) = E(i,1,k)
     Dj(iD+8) = E(i,2,k)
     iD = iD + 8
  end do

  i = nuC
  den = 1.0 / (nv-1)
  ! Bezier: left
  do j=2,nvC-1
     v = den * (j-1)
     t = 0.5 * tanh(a*(v-0.5)) / tanh(a*0.5) + 0.5
     do k=1,3
        call sparseBezier(t, -f0, f0, C)
        Da(iD+1:iD+4) = (1-v) * C(:)
        Di(iD+1:iD+4) = inds(i, j, k)
        Dj(iD+1) = W(i,1,k)
        Dj(iD+2) = W(i,2,k)
        Dj(iD+3) = E(i,1,k)
        Dj(iD+4) = E(i,2,k)
        iD = iD + 4
        call sparseBezier(dble(0.5), -f0, f0, C)
        Da(iD+1:iD+4) = v * C(:)
        Di(iD+1:iD+4) = inds(i, j, k)
        Dj(iD+1) = N(j,1,k)
        Dj(iD+2) = N(j,2,k)
        Dj(iD+3) = S(j,1,k)
        Dj(iD+4) = S(j,2,k)
        iD = iD + 4
     end do
  end do
  ! Bezier: right
  do j=nvC+1,nv-1
     v = den * (j-1)
     t = 0.5 * tanh(a*(v-0.5)) / tanh(a*0.5) + 0.5
     do k=1,3
        call sparseBezier(t, -f0, f0, C)
        Da(iD+1:iD+4) = v * C(:)
        Di(iD+1:iD+4) = inds(i, j, k)
        Dj(iD+1) = W(i,1,k)
        Dj(iD+2) = W(i,2,k)
        Dj(iD+3) = E(i,1,k)
        Dj(iD+4) = E(i,2,k)
        iD = iD + 4
        call sparseBezier(dble(0.5), -f0, f0, C)
        Da(iD+1:iD+4) = (1-v) * C(:)
        Di(iD+1:iD+4) = inds(i, j, k)
        Dj(iD+1) = N(j,1,k)
        Dj(iD+2) = N(j,2,k)
        Dj(iD+3) = S(j,1,k)
        Dj(iD+4) = S(j,2,k)
        iD = iD + 4
     end do
  end do

  j = nvC
  den = 1.0 / (nu-1)
  ! Bezier: top
  do i=2,nuC-1
     u = den * (i-1)
     t = 0.5 * tanh(a*(u-0.5)) / tanh(a*0.5) + 0.5
     do k=1,3
        call sparseBezier(t, -f0, f0, C)
        Da(iD+1:iD+4) = (1-u) * C(:)
        Di(iD+1:iD+4) = inds(i, j, k)
        Dj(iD+1) = N(j,1,k)
        Dj(iD+2) = N(j,2,k)
        Dj(iD+3) = S(j,1,k)
        Dj(iD+4) = S(j,2,k)
        iD = iD + 4
        call sparseBezier(dble(0.5), -f0, f0, C)
        Da(iD+1:iD+4) = u * C(:)
        Di(iD+1:iD+4) = inds(i, j, k)
        Dj(iD+1) = W(i,1,k)
        Dj(iD+2) = W(i,2,k)
        Dj(iD+3) = E(i,1,k)
        Dj(iD+4) = E(i,2,k)
        iD = iD + 4
     end do
  end do
  ! Bezier: bottom
  do i=nuC+1,nu-1
     u = den * (i-1)
     t = 0.5 * tanh(a*(u-0.5)) / tanh(a*0.5) + 0.5
     do k=1,3
        call sparseBezier(t, -f0, f0, C)
        Da(iD+1:iD+4) = u * C(:)
        Di(iD+1:iD+4) = inds(i, j, k)
        Dj(iD+1) = N(j,1,k)
        Dj(iD+2) = N(j,2,k)
        Dj(iD+3) = S(j,1,k)
        Dj(iD+4) = S(j,2,k)
        iD = iD + 4
        call sparseBezier(dble(0.5), -f0, f0, C)
        Da(iD+1:iD+4) = (1-u) * C(:)
        Di(iD+1:iD+4) = inds(i, j, k)
        Dj(iD+1) = W(i,1,k)
        Dj(iD+2) = W(i,2,k)
        Dj(iD+3) = E(i,1,k)
        Dj(iD+4) = E(i,2,k)
        iD = iD + 4
     end do
  end do

  if (iD .ne. nD) then
     print *, 'Error in computeConeWireframe', iD, nD
  end if

end subroutine computeConeWireframe
