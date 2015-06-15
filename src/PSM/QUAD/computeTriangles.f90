subroutine computeTriangles(nvert, nadj, ntri, adjPtr, adjMap, triangles)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nvert, nadj, ntri, adjPtr, adjMap
  !f2py intent(out) triangles
  !f2py depend(nvert) adjPtr
  !f2py depend(nadj) adjMap
  !f2py depend(ntri) triangles

  !Input
  integer, intent(in) ::  nvert, nadj, ntri, adjPtr(nvert,2), adjMap(nadj)

  !Output
  integer, intent(out) ::  triangles(ntri,3)

  !Working
  integer itri, ivert0, ivert1, ivert2, ivert3, i1, i2, i3

  itri = 0
  do ivert0=1,nvert
     do i1=adjPtr(ivert0,1),adjPtr(ivert0,2)
        ivert1 = adjMap(i1)
        do i2=adjPtr(ivert1,1),adjPtr(ivert1,2)
           ivert2 = adjMap(i2)
           do i3=adjPtr(ivert2,1),adjPtr(ivert2,2)
              ivert3 = adjMap(i3)
              if (ivert0 .eq. ivert3) then
                 itri = itri + 1
                 triangles(itri,1) = ivert1
                 triangles(itri,2) = ivert2
                 triangles(itri,3) = ivert3
              end if
           end do
        end do
     end do
  end do
  if (itri .ne. ntri) then
     print *, 'Error in computeTriangles', itri, ntri
     call exit(1)
  end if

end subroutine computeTriangles




subroutine countTriangles(nvert, nadj, adjPtr, adjMap, ntri)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nvert, nadj, adjPtr, adjMap
  !f2py intent(out) ntri
  !f2py depend(nvert) adjPtr
  !f2py depend(nadj) adjMap
  
  !Input
  integer, intent(in) ::  nvert, nadj, adjPtr(nvert,2), adjMap(nadj)

  !Output
  integer, intent(out) ::  ntri

  !Working
  integer ivert0, ivert1, ivert2, ivert3, i1, i2, i3

  ntri = 0
  do ivert0=1,nvert
     do i1=adjPtr(ivert0,1),adjPtr(ivert0,2)
        ivert1 = adjMap(i1)
        do i2=adjPtr(ivert1,1),adjPtr(ivert1,2)
           ivert2 = adjMap(i2)
           do i3=adjPtr(ivert2,1),adjPtr(ivert2,2)
              ivert3 = adjMap(i3)
              if (ivert0 .eq. ivert3) then
                 ntri = ntri + 1
              end if
           end do
        end do
     end do
  end do

end subroutine countTriangles




subroutine rotateTriangles(nvert, ntri, verts, triangles)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nvert, ntri, verts
  !f2py intent(inout) triangles
  !f2py depend(nvert) verts
  !f2py depend(ntri) triangles

  !Input
  integer, intent(in) ::  nvert, ntri
  double precision, intent(in) ::  verts(nvert,2)

  !Output
  integer, intent(inout) ::  triangles(ntri,3)

  !Working
  double precision v1(2), v2(2)
  integer itri

  do itri=1,ntri
     v1 = verts(triangles(itri,2),:) - verts(triangles(itri,1),:)
     v2 = verts(triangles(itri,3),:) - verts(triangles(itri,2),:)
     if ((v1(1)*v2(2) - v1(2)*v2(1)) .lt. 0) then
        triangles(itri,:) = triangles(itri,3:1:-1)
     end if
  end do

end subroutine rotateTriangles




subroutine computeTri2Edge(nedge, ntri, edges, triangles, edge2tri, tri2edge)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nedge, ntri, edges, triangles
  !f2py intent(out) edge2tri, tri2edge
  !f2py depend(nedge) edges, edge2tri
  !f2py depend(ntri) triangles, tri2edge

  !Input
  integer, intent(in) ::  nedge, ntri, edges(nedge,2), triangles(ntri,3)

  !Output
  integer, intent(out) ::  edge2tri(nedge,2), tri2edge(ntri,3)

  !Working
  integer itri, iedge, edgeT(3,2), sameEdge, k, res

  tri2edge(:,:) = 0
  edge2tri(:,:) = 0
  do itri=1,ntri
     edgeT(1,:) = (/ triangles(itri,1) , triangles(itri,2) /)
     edgeT(2,:) = (/ triangles(itri,2) , triangles(itri,3) /)
     edgeT(3,:) = (/ triangles(itri,3) , triangles(itri,1) /)
     do iedge=1,nedge
        do k=1,3
           res = sameEdge(edgeT(k,:), edges(iedge,:))
           if (res .eq. 1) then
              tri2edge(itri,k) = iedge
              edge2tri(iedge,1) = itri
           else if (res .eq. -1) then
              tri2edge(itri,k) = -iedge
              edge2tri(iedge,2) = itri
           end if
        end do
     end do
  end do

end subroutine computeTri2Edge
