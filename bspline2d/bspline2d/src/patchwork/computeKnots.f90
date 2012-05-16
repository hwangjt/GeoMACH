subroutine getD(ngroup, numD, group_k, group_m, group_d)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) ngroup, numD, group_k, group_m
  !f2py intent(out) group_d
  !f2py depend(ngroup) group_k
  !f2py depend(ngroup) group_m
  !f2py depend(numD) group_d

  !Input
  integer, intent(in) ::  ngroup, numD
  integer, intent(in) ::  group_k(ngroup), group_m(ngroup)

  !Output
  double precision, intent(out) ::  group_d(numD)

  !Working
  integer group, d1, d2
  integer k,m
  double precision, allocatable, dimension(:) ::  buffer
  
  d1= 1
  d2 = 1
  do group=1,ngroup
     k = group_k(group)
     m = group_m(group)
     d2 = d1 + k + m - 1
     allocate(buffer(k+m))
     call knotopen(k,k+m,buffer)
     group_d(d1:d2) = buffer(:)
     deallocate(buffer)
     d1 = d1 + k + m
  end do

end subroutine getD
