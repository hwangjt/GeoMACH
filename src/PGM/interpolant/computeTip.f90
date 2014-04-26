subroutine computeTip(nD, nu, nv, f0, m0, N, S, inds, Da, Di, Dj)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nD, nu, nv, f0, m0, N, S, inds
  !f2py intent(out) Da, Di, Dj
  !f2py depend(nv) N, S
  !f2py depend(nu,nv) inds
  !f2py depend(nD) Da, Di, Dj

  !Input
  integer, intent(in) ::  nD, nu, nv
  double precision, intent(in) ::  f0, m0
  integer, intent(in) ::  N(nv,2), S(nv,2)
  integer, intent(in) ::  inds(nu,nv)

  !Output
  double precision, intent(out) ::  Da(nD)
  integer, intent(out) ::  Di(nD), Dj(nD)

  !Working
  integer i, j, iD
  double precision den, C(4), u

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

  den = 1.0 / (nu-1)
  do i=2,nu-1
     u = den * (i-1)
     call sparseBezier(u, -f0, f0, C)
     do j=1,nv
        Da(iD+1:iD+4) = C(:)
        Di(iD+1:iD+4) = inds(i, j)
        Dj(iD+1) = N(j,1)
        Dj(iD+2) = N(j,2)
        Dj(iD+3) = S(j,1)
        Dj(iD+4) = S(j,2)
        iD = iD + 4
     end do
  end do

  if (iD .ne. nD) then
     print *, 'Error in computeTip', iD, nD
  end if


end subroutine computeTip
