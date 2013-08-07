subroutine computeFaceDimensions(ni, nj, nsurf, &
     surfs, surfEdgeLengths, idims, jdims)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) ni, nj, nsurf, surfs, surfEdgeLengths
  !f2py intent(out) idims, jdims
  !f2py depend(ni,nj) surfs
  !f2py depend(nsurf) surfEdgeLengths
  !f2py depend(ni) idims
  !f2py depend(nj) jdims

  !Input
  integer, intent(in) ::  ni, nj, nsurf, surfs(ni,nj)
  double precision, intent(in) ::  surfEdgeLengths(nsurf,2,2)

  !Output
  double precision, intent(out) ::  idims(ni+1), jdims(nj+1)

  !Working
  double precision idims1(ni+1), idims2(ni+1)
  double precision jdims1(nj+1), jdims2(nj+1)
  integer i, j

  idims1(1) = 0.
  idims2(1) = 0.
  do i=1,ni
     idims1(i+1) = idims1(i) + surfEdgeLengths(surfs(i,1),1,1)
     idims2(i+1) = idims2(i) + surfEdgeLengths(surfs(i,nj),1,2)
  end do
  idims1(:) = idims1(:)/idims1(ni+1)
  idims2(:) = idims2(:)/idims2(ni+1)

  jdims1(1) = 0.
  jdims2(1) = 0.
  do j=1,nj
     jdims1(j+1) = jdims1(j) + surfEdgeLengths(surfs(1,j),2,1)
     jdims2(j+1) = jdims2(j) + surfEdgeLengths(surfs(ni,j),2,2)
  end do
  jdims1(:) = jdims1(:)/jdims1(nj+1)
  jdims2(:) = jdims2(:)/jdims2(nj+1)

  idims = 0.5*(idims1 + idims2)
  jdims = 0.5*(jdims1 + jdims2)

  idims(1) = 0.
  jdims(1) = 0.
  idims(ni+1) = 1.
  jdims(nj+1) = 1.

end subroutine computeFaceDimensions
