subroutine identifySymmNodes(nnode, nodes, symm)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nnode, nodes
  !f2py intent(out) symm
  !f2py depend(nnode) nodes

  !Input
  integer, intent(in) ::  nnode
  double precision, intent(in) ::  nodes(nnode,3)

  !Output
  logical, intent(out) ::  symm(nnode)

  !Working
  integer inode

  do inode=1,nnode
     symm(inode) = abs(nodes(inode,3)) .lt. 1e-3
  end do

end subroutine identifySymmNodes
