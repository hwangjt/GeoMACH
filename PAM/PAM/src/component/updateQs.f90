subroutine updateQs(nQ, ni, nj, Nf, Qf, Q)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nQ, ni, nj, Nf, Qf
  !f2py intent(in,out) Q
  !f2py depend(ni,nj) Nf
  !f2py depend(ni,nj) Qf
  !f2py depend(nQ) Q

  !Input
  integer, intent(in) ::  nQ, ni, nj
  integer, intent(in) ::  Nf(ni, nj, 5)
  double precision, intent(in) ::  Qf(ni, nj, 3)

  !Output
  double precision, intent(inout) ::  Q(nQ,3)

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
