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
