subroutine computeAdjMap(nvert, nedge, nadj, edges, adjPtr, adjMap)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nvert, nedge, nadj, edges
  !f2py intent(out) adjPtr, adjMap
  !f2py depend(nedge) edges
  !f2py depend(nvert) adjPtr
  !f2py depend(nadj) adjMap

  !Input
  integer, intent(in) ::  nvert, nedge, nadj, edges(nedge,2)

  !Output
  integer, intent(out) ::  adjPtr(nvert,2), adjMap(nadj)

  !Working
  integer i, j, index, counter(nvert)

  adjPtr(:,:) = 0
  do i=1,nedge
     do j=1,2
        adjPtr(edges(i,j),2) = adjPtr(edges(i,j),2) + 1
     end do
  end do

  adjPtr(1,1) = 1
  do i=2,nvert
     adjPtr(i,1) = adjPtr(i-1,2) + 1
     adjPtr(i,2) = adjPtr(i-1,2) + adjPtr(i,2)
  end do

  counter(:) = 0
  do i=1,nedge
     do j=1,2
        index = adjPtr(edges(i,j),1) + counter(edges(i,j))
        adjMap(index) = edges(i,3-j)
        counter(edges(i,j)) = counter(edges(i,j)) + 1
     end do
  end do  

end subroutine computeAdjMap
