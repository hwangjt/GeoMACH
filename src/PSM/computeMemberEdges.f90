subroutine computeMemberEdges(imem, nmem, mem_group, edges, edge_group)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) imem, nmem, mem_group
  !f2py intent(out) edges, edge_group
  !f2py depend(nmem) mem_group

  !Input
  integer, intent(in) ::  imem, nmem, mem_group(nmem,2,2)

  !Output
  double precision, intent(out) ::  edges(4,2,2)
  integer, intent(out) ::  edge_group(4)

  !Working
  double precision verts(2,2,2)
  integer i, j

  do i=1,2
     do j=1,2
        verts(i,j,1) = dble(i-1)
        verts(i,j,2) = dble(j-1)
     end do
  end do

  edges(1,1,:) = verts(1,1,:)
  edges(1,2,:) = verts(2,1,:)
  edges(2,1,:) = verts(1,2,:)
  edges(2,2,:) = verts(2,2,:)
  edges(3,1,:) = verts(1,1,:)
  edges(3,2,:) = verts(1,2,:)
  edges(4,1,:) = verts(2,1,:)
  edges(4,2,:) = verts(2,2,:)

  edge_group(1) = mem_group(imem,1,1)
  edge_group(2) = mem_group(imem,1,2)
  edge_group(3) = mem_group(imem,2,1)
  edge_group(4) = mem_group(imem,2,2)

end subroutine computeMemberEdges
