subroutine mergeTriangles(nedge0, ntri0, nedge, ntri, nquad, &
     edgeCon0, edges0, triangles0, edge2tri, removeEdge, removeTri, &
     edgeCon, edges, triangles, quads)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nedge0, ntri0, nedge, ntri, nquad, edgeCon0, edges0, triangles0, edge2tri, removeEdge, removeTri
  !f2py intent(out) edgeCon, edges, triangles, quads
  !f2py depend(nedge0) edgeCon0, edges0, edge2tri, removeEdge
  !f2py depend(ntri0) triangles0, removeTri
  !f2py depend(nedge) edgeCon, edges
  !f2py depend(ntri) triangles
  !f2py depend(nquad) quads

  !Input
  integer, intent(in) ::  nedge0, ntri0, nedge, ntri, nquad
  logical, intent(in) ::  edgeCon0(nedge0)
  integer, intent(in) ::  edges0(nedge0,2), triangles0(ntri0,3)
  integer, intent(in) ::  edge2tri(nedge0,2)
  logical, intent(in) ::  removeEdge(nedge0), removeTri(ntri0)

  !Output
  logical, intent(out) ::  edgeCon(nedge)
  integer, intent(out) ::  edges(nedge,2), triangles(ntri,3), quads(nquad,4)

  !Working
  integer iedge0, itri0
  integer iedge, itri, iquad
  integer otherVert

  iedge = 0
  do iedge0=1,nedge0
     if (.not. removeEdge(iedge0)) then
        iedge = iedge + 1
        edges(iedge,:) = edges0(iedge0,:)
        edgeCon(iedge) = edgeCon0(iedge0)
     end if
  end do
  if (iedge .ne. nedge) then
     print *, 'Error in mergeTriangles: iedge'
     call exit(1)
  end if

  triangles(1,:) = 0
  itri = 1
  do itri0=1,ntri0
     if (.not. removeTri(itri0)) then
        itri = itri + 1
        triangles(itri,:) = triangles0(itri0,:)
     end if
  end do
  if (itri .ne. ntri) then
     print *, 'Error in mergeTriangles: itri'
     call exit(1)
  end if

  quads(1,:) = 0
  iquad = 1
  do iedge0=1,nedge0
     if (removeEdge(iedge0)) then
        iquad = iquad + 1
        quads(iquad,1) = edges0(iedge0,1)
        quads(iquad,2) = otherVert(edges0(iedge0,:),triangles0(edge2tri(iedge0,1),:))
        quads(iquad,3) = edges0(iedge0,2)
        quads(iquad,4) = otherVert(edges0(iedge0,:),triangles0(edge2tri(iedge0,2),:))
     end if
  end do

end subroutine mergeTriangles




function otherVert(edge, tri)

  implicit none
  integer, intent(in) ::  edge(2), tri(3)
  integer otherVert, sameEdge

  if (sameEdge(edge,tri(1:2)) .ne. 0) then
     otherVert = tri(3)
  else if (sameEdge(edge,tri(2:3)) .ne. 0) then
     otherVert = tri(1)
  else if (sameEdge(edge,tri(1:3:2)) .ne. 0) then
     otherVert = tri(2)
  else
     print *, 'Error in otherVert', edge, tri
     call exit(1)
  end if

end function otherVert




subroutine computeQuadDominant(nvert, nedge, ntri, &
     verts, edges, edgeCon, triangles, tri2edge, edge2tri, &
     nrem, removeEdge, removeTri)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nvert, nedge, ntri, verts, edges, edgeCon, triangles, tri2edge, edge2tri
  !f2py intent(out) nrem, removeEdge, removeTri
  !f2py depend(nvert) verts
  !f2py depend(nedge) edges, edgeCon, edge2tri, removeEdge
  !f2py depend(ntri) triangles, tri2edge, removeTri

  !Input
  integer, intent(in) ::  nvert, nedge, ntri
  double precision, intent(in) ::  verts(nvert,2)
  integer, intent(in) ::  edges(nedge,2)
  logical, intent(in) ::  edgeCon(nedge)
  integer, intent(in) ::  triangles(ntri,3), tri2edge(ntri,3), edge2tri(nedge,2)

  !Output
  integer, intent(out) ::  nrem
  logical, intent(out) ::  removeEdge(nedge), removeTri(ntri)

  !Working
  double precision anglesT(ntri,3), angles1(3), angles2(3), angles(4), maxa(nedge)
  double precision L1, L2, L3, getDist, pi
  integer itri, k, iedge, itri1, itri2

  pi = 2*acos(0.0)

  do itri=1,ntri
     L1 = getDist(verts(triangles(itri,1),:), verts(triangles(itri,2),:))
     L2 = getDist(verts(triangles(itri,2),:), verts(triangles(itri,3),:))
     L3 = getDist(verts(triangles(itri,3),:), verts(triangles(itri,1),:))
     anglesT(itri,1) = acos((L1**2 + L3**2 - L2**2)/2/L1/L3)
     anglesT(itri,2) = acos((L2**2 + L1**2 - L3**2)/2/L2/L1)
     anglesT(itri,3) = acos((L3**2 + L2**2 - L1**2)/2/L3/L2)
  end do

  do iedge=1,nedge
     if (edgeCon(iedge)) then
        maxa(iedge) = 2*pi
     else
        itri = edge2tri(iedge,1)
        if (edges(iedge,1) .eq. triangles(itri,1)) then
           angles1 = (/ anglesT(itri,1), anglesT(itri,2), anglesT(itri,3) /)
        else if (edges(iedge,1) .eq. triangles(itri,2)) then
           angles1 = (/ anglesT(itri,2), anglesT(itri,3), anglesT(itri,1) /)
        else if (edges(iedge,1) .eq. triangles(itri,3)) then
           angles1 = (/ anglesT(itri,3), anglesT(itri,1), anglesT(itri,2) /)
        else
           print *, 'Error in computeQuadDominant1: edge2tri wrong', edges(iedge,1), triangles(itri,:)
           call exit(1)
        end if
        itri = edge2tri(iedge,2)
        if (edges(iedge,1) .eq. triangles(itri,1)) then
           angles2 = (/ anglesT(itri,1), anglesT(itri,2), anglesT(itri,3) /)
        else if (edges(iedge,1) .eq. triangles(itri,2)) then
           angles2 = (/ anglesT(itri,2), anglesT(itri,3), anglesT(itri,1) /)
        else if (edges(iedge,1) .eq. triangles(itri,3)) then
           angles2 = (/ anglesT(itri,3), anglesT(itri,1), anglesT(itri,2) /)
        else
           print *, 'Error in computeQuadDominant2: edge2tri wrong', edges(iedge,1), triangles(itri,:)
           call exit(1)
        end if
        angles(1) = angles1(1) + angles2(1)
        angles(2) = angles1(2) + angles2(3)
        angles(3) = angles1(3)
        angles(4) = angles2(2)
        maxa(iedge) = maxval(angles)
     end if
  end do

  nrem = 0
  removeEdge(:) = .False.
  removeTri(:) = .False.
  main: do k=1,1000
     iedge = minloc(maxa,1)
     if (maxa(iedge) .lt. 3*pi/4.0) then
        removeEdge(iedge) = .True.
        if (edge2tri(iedge,1) .ne. 0) then
           removeTri(edge2tri(iedge,1)) = .True.
        end if
        if (edge2tri(iedge,2) .ne. 0) then
           removeTri(edge2tri(iedge,2)) = .True.
        end if
        nrem = nrem + 1
        itri1 = edge2tri(iedge,1)
        itri2 = edge2tri(iedge,2)
        do iedge=1,3
           maxa(abs(tri2edge(itri1,iedge))) = 2*pi
           maxa(abs(tri2edge(itri2,iedge))) = 2*pi
        end do
     else
        exit main
     end if
  end do main

end subroutine computeQuadDominant




function getDist(A,B)

  implicit none
  double precision, intent(in) ::  A(2), B(2)
  double precision getDist
  
  getDist = sqrt(dot_product(B-A,B-A))

end function getDist




subroutine removeEdges(nedge0, nedge, remove, edges0, edgeCon0, edges, edgeCon)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nedge0, nedge, remove, edges0, edgeCon0
  !f2py intent(out) edges, edgeCon
  !f2py depend(nedge0) remove, edges0, edgeCon0
  !f2py depend(nedge) edges, edgeCon

  !Input
  integer, intent(in) ::  nedge0, nedge
  logical, intent(in) ::  remove(nedge0)
  integer, intent(in) ::  edges0(nedge0,2)
  logical, intent(in) ::  edgeCon0(nedge0)

  !Outpout
  integer, intent(out) ::  edges(nedge,2)
  logical, intent(out) ::  edgeCon(nedge)

  !Working
  integer iedge, iedge0

  iedge = 0
  do iedge0=1,nedge0
     if (.not. remove(iedge0)) then
        iedge = iedge + 1
        edges(iedge,:) = edges0(iedge0,:)
        edgeCon(iedge) = edgeCon0(iedge0)
     end if
  end do
  if (iedge .ne. nedge) then
     print *, 'Error in removeEdges', iedge, nedge
     call exit(1)
  end if

end subroutine removeEdges
