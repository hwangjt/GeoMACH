subroutine countQuadSplits(nedge, ngroup, npoly, poly_edge, &
     edge_group, group_split, nsplit)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nedge, ngroup, npoly, poly_edge, edge_group, group_split
  !f2py intent(out) nsplit
  !f2py depend(npoly) poly_edge
  !f2py depend(nedge) edge_group
  !f2py depend(ngroup) group_split

  !Input
  integer, intent(in) ::  nedge, ngroup, npoly
  integer, intent(in) ::  poly_edge(npoly,5), edge_group(nedge)
  logical, intent(in) ::  group_split(ngroup)

  !Output
  integer, intent(out) ::  nsplit

  !Working
  integer p

  nsplit = 0
  do p=1,npoly
     if (poly_edge(p,4) .ne. 0) then
        if (group_split(edge_group(poly_edge(p,1)))) then
           nsplit = nsplit + 1
        end if
        if (group_split(edge_group(poly_edge(p,2)))) then
           nsplit = nsplit + 1
        end if
     end if
  end do

end subroutine countQuadSplits



subroutine computeTriSplits(nedge, ngroup, npoly, edge_group, poly_edge, group_split)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nedge, ngroup, npoly, edge_group, poly_edge
  !f2py intent(out) group_split
  !f2py depend(nedge) edge_group
  !f2py depend(npoly) poly_edge
  !f2py depend(ngroup) group_split

  !Input
  integer, intent(in) ::  nedge, ngroup, npoly
  integer, intent(in) :: edge_group(nedge), poly_edge(npoly,5)

  !Output
  logical, intent(out) ::  group_split(ngroup)

  !Working
  integer p, k

  group_split(:) = .False.
  do p=1,npoly
     if (poly_edge(p,4) .eq. 0) then
        do k=1,3
           group_split(edge_group(poly_edge(p,k))) = .True.
        end do
     end if
  end do

end subroutine computeTriSplits



subroutine addPolySplits(nvert, nedge, nvert0, nedge0, ngroup, npoly, &
     verts0, edges0, edge_group, poly_vert, poly_edge, group_split, verts, edges)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nvert, nedge, nvert0, nedge0, ngroup, npoly, verts0, edges0, edge_group, poly_vert, poly_edge, group_split
  !f2py intent(out) verts, edges
  !f2py depend(nvert0) verts0
  !f2py depend(nedge0) edges0
  !f2py depend(nedge0) edge_group
  !f2py depend(npoly) poly_vert, poly_edge
  !f2py depend(ngroup) group_split
  !f2py depend(nvert) verts
  !f2py depend(nedge) edges

  !Input
  integer, intent(in) ::  nvert, nedge, nvert0, nedge0, ngroup, npoly
  double precision, intent(in) ::  verts0(nvert0,2), edges0(nedge0,5)
  integer, intent(in) :: edge_group(nedge0), poly_vert(npoly,5), poly_edge(npoly,5)
  logical, intent(in) ::  group_split(ngroup)

  !Output
  double precision, intent(out) ::  verts(nvert,2), edges(nedge,5)

  !Working
  integer v, e, p, k
  double precision vert(4,2), avg(2), mid(4,2)

  do v=1,nvert0
     verts(v,:) = verts0(v,:)
  end do
  v = nvert0

  do e=1,nedge0
     edges(e,:) = edges0(e,:)
  end do
  e = nedge0

  do p=1,npoly
     if (poly_edge(p,4) .eq. 0) then
        do k=1,3
           vert(k,:) = verts(poly_vert(p,k),:)
        end do
        avg = 1.0/3.0*vert(1,:) + 1.0/3.0*vert(2,:) + 1.0/3.0*vert(3,:)
        mid(1,:) = 0.5*vert(1,:) + 0.5*vert(2,:)
        mid(2,:) = 0.5*vert(2,:) + 0.5*vert(3,:)
        mid(3,:) = 0.5*vert(3,:) + 0.5*vert(1,:)
        verts(v+1,:) = avg(:)
        do k=1,3
           verts(v+1+k,:) = mid(k,:)
        end do
        do k=1,3
           edges(e+k,1) = v+1
           edges(e+k,2) = v+1+k
        end do
        v = v + 4
        e = e + 3
     end if
  end do

  do p=1,npoly
     if (poly_edge(p,4) .ne. 0) then
        do k=1,4
           vert(k,:) = verts(poly_vert(p,k),:)
        end do
        mid(1,:) = 0.5*vert(1,:) + 0.5*vert(2,:)
        mid(2,:) = 0.5*vert(2,:) + 0.5*vert(3,:)
        mid(3,:) = 0.5*vert(3,:) + 0.5*vert(4,:)
        mid(4,:) = 0.5*vert(4,:) + 0.5*vert(1,:)
        if (group_split(edge_group(poly_edge(p,1)))) then
           verts(v+1,:) = mid(2,:)
           verts(v+2,:) = mid(4,:)
           edges(e+1,1) = v+1
           edges(e+1,2) = v+2
           v = v + 2
           e = e + 1
        end if
        if (group_split(edge_group(poly_edge(p,2)))) then
           verts(v+1,:) = mid(1,:)
           verts(v+2,:) = mid(3,:)
           edges(e+1,1) = v+1
           edges(e+1,2) = v+2
           v = v + 2
           e = e + 1
        end if
     end if
  end do

end subroutine addPolySplits



subroutine splitPentagons(nedge, nvert, nedge0, npoly, Lx, Ly, &
     verts, edges0, poly_vert, edges)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nedge, nvert, nedge0, npoly, Lx, Ly, verts, edges0, poly_vert
  !f2py intent(out) edges
  !f2py depend(nvert) verts
  !f2py depend(nedge0) edges0
  !f2py depend(npoly) poly_vert
  !f2py depend(nedge) edges

  !Input
  integer, intent(in) ::  nedge, nvert, nedge0, npoly
  double precision, intent(in) ::  Lx, Ly
  double precision, intent(in) ::  verts(nvert,2), edges0(nedge0,5)
  integer, intent(in) ::  poly_vert(npoly,5)

  !Output
  double precision, intent(out) ::  edges(nedge,5)

  !Working
  integer p, e, i, index, vi(5)
  double precision v(6,2), t(6), a(5), L2(5)
  double precision pi

  pi = 2*acos(0.0)

  do e=1,nedge0
     edges(e,:) = edges0(e,:)
  end do
  e = nedge0 + 1

  v(:,:) = 0
  t(:) = 0
  do p=1,npoly
     if (poly_vert(p,5) .ne. 0) then
        vi(:) = poly_vert(p,:)
        do i=1,5
           v(i,:) = verts(vi(i),:)
        end do
        v(6,:) = v(1,:)
        do i=1,5
           L2(i) = dot_product(v(i+1,:)-v(i,:),v(i+1,:)-v(i,:))
           call arc_tan(v(i+1,:)-v(i,:), Lx, Ly, t(i+1))
           t(i+1) = t(i+1)/pi
        end do
        t(1) = t(6)
        do i=1,5
           a(i) = t(i+1) - t(i)
           if (a(i) .lt. 0) then
              a(i) = a(i) + 2.0
           end if
        end do
        index = maxloc(a,1)
        vi = cshift(vi,index-1)
        a = cshift(a,index-1)
        L2 = cshift(L2,index-1)
        edges(e,1) = vi(1)
        if (L2(1) .gt. L2(5)) then
           edges(e,2) = vi(3)
        else
           edges(e,2) = vi(4)
        end if
        e = e + 1
     end if
  end do

end subroutine splitPentagons



subroutine computeGroups(nedge, npoly, poly_edge, edge_group)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nedge, npoly, poly_edge
  !f2py intent(out) edge_group
  !f2py depend(npoly) poly_edge
  !f2py depend(nedge) edge_group

  !Input
  integer, intent(in) ::  nedge, npoly
  integer, intent(in) ::  poly_edge(npoly,5)

  !Output
  integer, intent(out) ::  edge_group(nedge)

  !Working
  integer e, p, k
  integer g(4)
  integer e1, e2
  integer edge_group0(nedge)
  integer index

  do e=1,nedge
     edge_group0(e) = e
  end do

  do p=1,npoly
     if ((poly_edge(p,4).ne.0) .and. (poly_edge(p,5).eq.0)) then
        do k=1,4
           g(k) = edge_group0(poly_edge(p,k))
        end do
        call setAll(nedge, g(1), g(3), edge_group0)
        call setAll(nedge, g(2), g(4), edge_group0)
     end if
  end do

  edge_group(:) = 0
  
  index = 0
  do e1=1,nedge
     if (edge_group0(e1) .ne. 0) then
        index = index + 1
        edge_group(e1) = index
        do e2=e1+1,nedge
           if (edge_group0(e1) .eq. edge_group0(e2)) then
              edge_group(e2) = index
              edge_group0(e2) = 0
           end if
        end do
        edge_group0(e1) = 0
     end if
  end do

end subroutine computeGroups



subroutine setAll(nedge, g1, g2, edge_group)

  implicit none

  !Input
  integer, intent(in) ::  nedge, g1, g2
  
  !Input/Output
  integer, intent(inout) ::  edge_group(nedge)

  !Working
  integer e
  
  do e=1,nedge
     if (edge_group(e) .eq. g1) then
        edge_group(e) = g2
     end if
  end do

end subroutine setAll
