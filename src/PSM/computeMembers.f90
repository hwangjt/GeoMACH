subroutine computeMemberProj(isurf, isrc, nnode, npts, &
     nodesInt, nodesFlt, inds, P, Q)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) isurf, isrc, nnode, pts, nodesInt, nodesFlt
  !f2py intent(out) inds, P, Q
  !f2py depend(nnode) nodesInt, nodesFlt
  !f2py depend(npts) inds, P, Q

  !Input
  integer, intent(in) ::  isurf, isrc, nnode, npts, nodesInt(nnode,4,3)
  double precision, intent(in) ::  nodesFlt(nnode,4,3)

  !Output
  integer, intent(out) ::  inds(npts)
  double precision, intent(out) ::  P(npts,3), Q(npts,3)

  !Working
  integer inode, ipts

  P(:,:) = 0.
  Q(:,:) = 0.
  Q(:,3) = 1.

  ipts = 1
  do inode=1,nnode
     if (nodesInt(inode,isrc,3) .eq. isurf) then
        inds(ipts) = inode
        P(ipts,1) = nodesFlt(inode,isrc,2)
        P(ipts,2) = nodesFlt(inode,isrc,3)
        ipts = ipts + 1
     end if
  end do

end subroutine computeMemberProj




subroutine countMembers(isurf, isrc, nnode, nodesInt, npts)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) isurf, isrc, nnode, nodesInt
  !f2py intent(out) npts
  !f2py depend(nnode) nodesInt

  !Input
  integer, intent(in) ::  isurf, isrc, nnode, nodesInt(nnode,4,3)

  !Output
  integer, intent(out) ::  npts

  !Working
  integer inode

  npts = 0
  do inode=1,nnode
     if (nodesInt(inode,isrc,3) .eq. isurf) then
        npts = npts + 1
     end if
  end do

end subroutine countMembers




subroutine computeMemberLocalCoords(comp, face, ni, nj, nnode, idims, jdims, surfs, nodesInt, nodesFlt)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) comp, face, ni, nj, nnode, idims, jdims, surfs
  !f2py intent(inout) nodesInt, nodesFlt
  !f2py depend(ni) idims
  !f2py depend(nj) jdims
  !f2py depend(ni,nj) surfs
  !f2py depend(nnode) nodesInt, nodesFlt

  !Input
  integer, intent(in) ::  comp, face, ni, nj, nnode
  double precision, intent(in) ::  idims(ni+1), jdims(nj+1)
  integer, intent(in) ::  surfs(ni,nj)

  !Output
  integer, intent(inout) ::  nodesInt(nnode,4,3)
  double precision, intent(inout) ::  nodesFlt(nnode,4,3)

  !Working
  double precision a, b, u, v
  integer inode, isrc, i, j, ifound, jfound

  do inode=1,nnode
     do isrc=1,4
        if ((nodesInt(inode,isrc,1) .eq. comp) .and. &
             (nodesInt(inode,isrc,2) .eq. face)) then
           a = nodesFlt(inode,isrc,2)
           b = nodesFlt(inode,isrc,3)
           iloop: do i=1,ni
              if (a .le. idims(i+1)) then
                 ifound = i
                 u = (a-idims(i))/(idims(i+1)-idims(i))
                 exit iloop
              end if
           end do iloop
           jloop: do j=1,nj
              if (b .le. jdims(j+1)) then
                 jfound = j
                 v = (b-jdims(j))/(jdims(j+1)-jdims(j))
                 exit jloop
              end if
           end do jloop
           nodesInt(inode,isrc,3) = surfs(ifound,jfound)
           nodesFlt(inode,isrc,2) = u
           nodesFlt(inode,isrc,3) = v
        end if
     end do
  end do

end subroutine computeMemberLocalCoords
