subroutine computeProps(nD, mu, mv, nu, nv, B, Tu, Tv, &
     param_inds, prop_inds, Da, Di, Dj)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nD, mu, mv, nu, nv, B, Tu, Tv, param_inds, prop_inds
  !f2py intent(out) Da, Di, Dj
  !f2py depend(mu,mv) B
  !f2py depend(mu) Tu
  !f2py depend(mv) Tv
  !f2py depend(mu,mv) param_inds
  !f2py depend(nu,nv) prop_inds
  !f2py depend(nD) Da, Di, Dj

  !Input
  integer, intent(in) ::  nD, mu, mv, nu, nv
  logical, intent(in) ::  B(mu,mv,2)
  double precision, intent(in) ::  Tu(mu), Tv(mv)
  integer, intent(in) ::  param_inds(mu,mv,3), prop_inds(nu,nv)

  !Output
  double precision, intent(out) ::  Da(nD)
  integer, intent(out) ::  Di(nD), Dj(nD)

  !Working
  integer i, j, i0, j0, iD
  double precision u, v, u0, v0, denu, denv, C(4)

  iD = 0

  Da(:) = 0.0
  Di(:) = 0
  Dj(:) = 0

  denu = 1.0 / (nu-1)
  denv = 1.0 / (nv-1)
  do i0=1,nu
     u0 = 1.0*(i0-1) / (nu-1)
     call locateParameter(mu, u0, Tu, i, u)
     do j0=1,nv
        v0 = 1.0*(j0-1) / (nv-1)
        call locateParameter(mv, v0, Tv, j, v)
        Da(iD+1:iD+12) = 0.0
        Di(iD+1:iD+12) = prop_inds(i0,j0)
        if ((i .gt. 0) .and. (j .gt. 0) .and. &
             (i .lt. mu) .and. (j .lt. mv)) then
           Dj(iD+1:iD+3) = param_inds(i,j,:)
           Dj(iD+4:iD+6) = param_inds(i+1,j,:)
           Dj(iD+7:iD+9) = param_inds(i,j+1,:)
           Dj(iD+10:iD+12) = param_inds(i+1,j+1,:)

           call sparseCubic(1, u, B(i,j,:), B(i+1,j,:), C)
           Da(iD+1) = Da(iD+1) + C(1) * (1-v)
           Da(iD+2) = Da(iD+2) + C(2) * (1-v)
           Da(iD+4) = Da(iD+4) + C(3) * (1-v)
           Da(iD+5) = Da(iD+5) + C(4) * (1-v)
           call sparseCubic(1, u, B(i,j+1,:), B(i+1,j+1,:), C)
           Da(iD+7) = Da(iD+7) + C(1) * v
           Da(iD+8) = Da(iD+8) + C(2) * v
           Da(iD+10) = Da(iD+10) + C(3) * v
           Da(iD+11) = Da(iD+11) + C(4) * v
           call sparseCubic(2, v, B(i,j,:), B(i,j+1,:), C)
           Da(iD+1) = Da(iD+1) + C(1) * (1-u)
           Da(iD+3) = Da(iD+3) + C(2) * (1-u)
           Da(iD+7) = Da(iD+7) + C(3) * (1-u)
           Da(iD+9) = Da(iD+9) + C(4) * (1-u)
           call sparseCubic(2, v, B(i+1,j,:), B(i+1,j+1,:), C)
           Da(iD+4) = Da(iD+4) + C(1) * u
           Da(iD+6) = Da(iD+6) + C(2) * u
           Da(iD+10) = Da(iD+10) + C(3) * u
           Da(iD+12) = Da(iD+12) + C(4) * u

           Da(iD+1) = Da(iD+1) - (1-u)*(1-v)
           Da(iD+4) = Da(iD+4) - u*(1-v)
           Da(iD+7) = Da(iD+7) - (1-u)*v
           Da(iD+10) = Da(iD+10) - u*v
        else if ((i .eq. 0) .and. (j .eq. 0) .or. &
             ((mu .eq. 1) .and. (mv .eq. 1))) then
           Da(iD+1) = 1.0
           Dj(iD+1) = param_inds(1,1,1)
        else if ((i .eq. 0) .and. (j .gt. 0) .or. (mu .eq. 1)) then
           i = 1
           Dj(iD+1:iD+3) = param_inds(i,j,:)
           Dj(iD+7:iD+9) = param_inds(i,j+1,:)
           call sparseCubic(2, v, B(i,j,:), B(i,j+1,:), C)
           Da(iD+1) = C(1)
           Da(iD+3) = C(2)
           Da(iD+7) = C(3)
           Da(iD+9) = C(4)
        else if ((i .gt. 0) .and. (j .eq. 0) .or. (mv .eq. 1)) then
           j = 1
           Dj(iD+1:iD+3) = param_inds(i,j,:)
           Dj(iD+4:iD+6) = param_inds(i+1,j,:)
           call sparseCubic(1, u, B(i,j,:), B(i+1,j,:), C)
           Da(iD+1) = C(1)
           Da(iD+2) = C(2)
           Da(iD+4) = C(3)
           Da(iD+5) = C(4)
        !else
        !   print *, 'Error in computeProps'
        end if
        iD = iD + 12
     end do
  end do

  if (iD .ne. nD) then
     print *, 'Error in computeProps', iD, nD
  end if

end subroutine computeProps



subroutine sparseCubic(dim, t, W0, W1, C)

  implicit none

  integer, intent(in) ::  dim
  double precision, intent(in) ::  t
  logical, intent(in) ::  W0(2), W1(2)
  double precision, intent(out) ::  C(4)
  double precision tA, tB, tC, tD

  tA = (1-t)**3
  tB = 3*t*(1-t)**2
  tC = 3*t**2*(1-t)
  tD = t**3

  if (W0(dim) .and. W1(dim)) then
     C(1) = tA + tB
     C(2) = tB/3.0
     C(3) = tC + tD
     C(4) = -tC/3.0
  else if (W0(dim)) then
     C(1) = tA + tB + 2.0/3.0*tC
     C(2) = tB/3.0 + tC/3.0
     C(3) = tC/3.0 + tD
     C(4) = 0.0
  else if (W1(dim)) then
     C(1) = tA + tB/3.0
     C(2) = 0.0
     C(3) = 2.0/3.0*tB + tC + tD
     C(4) = -tB/3.0 - tC/3.0
  else
     C(1) = tA + 2.0/3.0*tB + tC/3.0
     C(2) = 0.0
     C(3) = tB/3.0 + 2.0/3.0*tC + tD
     C(4) = 0.0
  end if  

end subroutine sparseCubic



subroutine locateParameter(n, u0, T, i, u)

  implicit none

  !Input
  integer, intent(in) ::  n
  double precision, intent(in) ::  u0, T(n)

  !Output
  integer, intent(out) ::  i
  double precision, intent(out) ::  u

  !Working
  integer k

  i = -1
  u = 0.0
  do k=1,n-1
     if ((T(k) .le. u0) .and. (u0 .le. T(k+1))) then
        i = k
        u = (u0-T(k))/(T(k+1)-T(k))
     end if
  end do
  if (n .eq. 1) then
     i = 0
  end if

end subroutine locateParameter
