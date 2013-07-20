subroutine computeEdgeLengths(nvar, nquad, nodes, quads, lengths)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nvar, nquad, nodes, quads
  !f2py intent(out) lengths
  !f2py depend(nquad,nvar) nodes
  !f2py depend(nquad) quads, lengths

  !Input
  integer, intent(in) ::  nvar, nquad
  double precision, intent(in) ::  nodes(4*nquad,nvar)
  integer, intent(in) ::  quads(nquad,4)

  !Output
  double precision, intent(out) ::  lengths(nquad,2,2)

  !Working
  integer iquad
  double precision getEdgeLength

  do iquad=1,nquad
     lengths(iquad,1,1) = getEdgeLength(quads(iquad,2),quads(iquad,1),nvar,nquad,nodes)
     lengths(iquad,1,2) = getEdgeLength(quads(iquad,4),quads(iquad,3),nvar,nquad,nodes)
     lengths(iquad,2,1) = getEdgeLength(quads(iquad,4),quads(iquad,1),nvar,nquad,nodes)
     lengths(iquad,2,2) = getEdgeLength(quads(iquad,3),quads(iquad,2),nvar,nquad,nodes)
  end do

end subroutine computeEdgeLengths




function getEdgeLength(i1, i2, nvar, nquad, nodes)

  implicit none

  !Input
  integer, intent(in) ::  i1, i2, nvar, nquad
  double precision, intent(in) ::  nodes(4*nquad,nvar)

  !Output
  double precision getEdgeLength

  !Working
  double precision V(nvar)

  V = nodes(i2,:) - nodes(i1,:)
  getEdgeLength = dot_product(V,V)**0.5

end function getEdgeLength
