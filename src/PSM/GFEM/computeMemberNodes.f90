subroutine computeMemberNodes(imem, nmem, nnode, &
     membersInt, membersFlt, nodes, nodesInt, nodesFlt)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) imem, nmem, nnode, membersInt, membersFlt, nodes
  !f2py intent(out) nodesInt, nodesFlt
  !f2py depend(nmem) membersInt, membersFlt
  !f2py depend(nnode) nodes, nodesInt, nodesFlt

  !Input
  integer, intent(in) ::  imem, nmem, nnode, membersInt(nmem,4,2)
  double precision, intent(in) ::  membersFlt(nmem,4,2,2,6), nodes(nnode,2)

  !Output
  integer, intent(out) ::  nodesInt(nnode,4,3)
  double precision, intent(out) ::  nodesFlt(nnode,4,3)

  !Working
  integer inode, isrc
  double precision u, v, maxu, maxv, w(2,2), a(2,2), b(2,2), bilinear

  maxu = maxval(nodes(:,1))
  maxv = maxval(nodes(:,2))

  do inode=1,nnode
     u = nodes(inode,1)/maxu
     v = nodes(inode,2)/maxv
     nodesInt(inode,:,1:2) = membersInt(imem,:,:)
     nodesInt(inode,:,3) = 0
     do isrc=1,4
        w = membersFlt(imem,isrc,:,:,1)
        a = membersFlt(imem,isrc,:,:,2)
        b = membersFlt(imem,isrc,:,:,3)
        nodesFlt(inode,isrc,1) = bilinear(u,v,membersFlt(imem,isrc,:,:,1))
        nodesFlt(inode,isrc,2) = bilinear(u,v,membersFlt(imem,isrc,:,:,2))
        nodesFlt(inode,isrc,3) = bilinear(u,v,membersFlt(imem,isrc,:,:,3))
     end do
  end do

end subroutine computeMemberNodes



function bilinear(u,v,A)

  implicit none

  !Input
  double precision u, v, A(2,2)

  !Output
  double precision bilinear

  bilinear = (1-u)*(1-v)*A(1,1) + u*(1-v)*A(2,1) &
       + (1-u)*v*A(1,2) + u*v*A(2,2)

end function bilinear
