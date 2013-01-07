subroutine updateQs(nQ, ni, nj, nvar, Nf, Qf, Q)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nQ, ni, nj, nvar, Nf, Qf
  !f2py intent(in,out) Q
  !f2py depend(ni,nj) Nf
  !f2py depend(ni,nj) Qf
  !f2py depend(nQ,nvar) Q

  !Input
  integer, intent(in) ::  nQ, ni, nj, nvar
  integer, intent(in) ::  Nf(ni, nj, 5)
  double precision, intent(in) ::  Qf(ni, nj, 3)

  !Output
  double precision, intent(inout) ::  Q(nQ, nvar)

  !Working
  integer i, j, k

  do j=1,nj
     do i=1,ni
        if (Nf(i,j,1) .ne. -1) then
           do k=1,3
              Q(Nf(i,j,1)+1,k) = Qf(i,j,k)
           end do
        end if
     end do
  end do

end subroutine updateQs



subroutine inflateVector(ni, nj, nV, V, Q)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) ni, nj, nV, V
  !f2py intent(out) Q
  !f2py depend(nV) V
  !f2py depend(ni,nj) Q

  !Input
  integer, intent(in) ::  ni, nj, nV
  double precision, intent(in) ::  V(nV)

  !Output
  double precision, intent(out) ::  Q(ni,nj,3)

  !Working
  integer i, j, k

  do k=1,3
     do j=1,nj
        do i=1,ni
           Q(i,j,k) = V(ni*nj*(k-1) + ni*(j-1) + i)
        end do
     end do
  end do

end subroutine inflateVector



subroutine bilinearInterp(n, ni, nj, i, j, C, P)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) n, ni, nj, i, j, C
  !f2py intent(out) P
  !f2py depend(n) P

  !Input
  integer, intent(in) ::  n, ni, nj, i, j
  double precision, intent(in) ::  C(2,2,3)

  !Output
  double precision, intent(out) ::  P(n,n,3)

  !Working
  integer u, v, uu, vv, nu, nv

  nu = ni*(n-1)
  nv = nj*(n-1)
  do v=1,n
     vv = (j-1)*(n-1) + v-1
     do u=1,n
        uu = (i-1)*(n-1) + u-1
        P(u,v,:) = (nu-uu)*(nv-vv)*C(1,1,:) + &
             (nu-uu)*vv*C(1,2,:) + &
             uu*(nv-vv)*C(2,1,:) + &
             uu*vv*C(2,2,:)
     end do
  end do
  P(:,:,:) = P(:,:,:)/(nu*nv)

end subroutine bilinearInterp
