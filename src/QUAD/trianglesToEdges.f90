subroutine trianglesToEdges(ntri0, ntri, nedge, triangles, edges)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) ntri0, ntri, nedge, triangles
  !f2py intent(out) edges
  !f2py depend(ntri0) triangles
  !f2py depend(nedge) edges

  !Input
  integer, intent(in) ::  ntri0, ntri, nedge, triangles(ntri0,3)

  !Output
  integer, intent(out) ::  edges(nedge,2)

  !Working
  integer itri

  do itri=1,ntri
     edges(3*itri-2,1) = triangles(itri,1)
     edges(3*itri-2,2) = triangles(itri,2)
     edges(3*itri-1,1) = triangles(itri,2)
     edges(3*itri-1,2) = triangles(itri,3)
     edges(3*itri-0,1) = triangles(itri,3)
     edges(3*itri-0,2) = triangles(itri,1)
  end do

end subroutine trianglesToEdges
