subroutine countDuplicateEdges(nedge, edges, count)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nedge, edges
  !f2py intent(out) count
  !f2py depend(nedge) edges

  !Input
  integer, intent(in) ::  nedge
  double precision, intent(in) ::  edges(nedge,5)

  !Output
  integer, intent(out) ::  count

  !Working
  integer e1, e2
  integer edgeID(nedge), index
  logical found1, found2

  edgeID(:) = 0

  index = 0
  do e1=1,nedge
     if (edgeID(e1) .eq. 0) then
        index = index + 1
        edgeID(e1) = index
        do e2=e1+1,nedge
           if (edgeID(e2) .eq. 0) then
              found1 = (edges(e1,1).eq.edges(e2,1)).and.(edges(e1,2).eq.edges(e2,2))
              found2 = (edges(e1,1).eq.edges(e2,2)).and.(edges(e1,2).eq.edges(e2,1))
              if (found1 .or. found2) then
                 edgeID(e2) = index                 
              end if
           end if
        end do
     end if
  end do
  count = nedge - index

end subroutine countDuplicateEdges



subroutine deleteDuplicateEdges(nedge, nedge0, edges0, edges)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nedge, nedge0, edges0
  !f2py intent(out) edges
  !f2py depend(nedge0) edges0
  !f2py depend(nedge) edges

  !Input
  integer, intent(in) ::  nedge, nedge0
  double precision, intent(in) ::  edges0(nedge0,5)

  !Output
  double precision, intent(out) ::  edges(nedge,5)

  !Working
  integer e1, e2
  integer edgeID(nedge0), index
  logical found1, found2

  edgeID(:) = 0

  index = 0
  do e1=1,nedge0
     if (edgeID(e1) .eq. 0) then
        index = index + 1
        edgeID(e1) = index
        edges(index,:) = edges0(e1,:)
        do e2=e1+1,nedge0
           if (edgeID(e2) .eq. 0) then
              found1 = (edges0(e1,1).eq.edges0(e2,1)).and.(edges0(e1,2).eq.edges0(e2,2))
              found2 = (edges0(e1,1).eq.edges0(e2,2)).and.(edges0(e1,2).eq.edges0(e2,1))
              if (found1 .or. found2) then
                 edgeID(e2) = index
                 if (edges(index,3) .gt. edges0(e2,3)) then
                    edges(index,:) = edges0(e2,:)
                 end if
              end if
           end if
        end do
     end if
  end do

end subroutine deleteDuplicateEdges



subroutine countDuplicateVerts(nvert, verts, count)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nvert, verts
  !f2py intent(out) count
  !f2py depend(nvert) verts

  !Input
  integer, intent(in) ::  nvert
  double precision, intent(in) ::  verts(nvert,2)

  !Output
  integer, intent(out) ::  count

  !Working
  integer v1, v2
  integer vertID(nvert), index
  double precision dP(2), norm

  vertID(:) = 0
  
  index = 0
  do v1=1,nvert
     if (vertID(v1) .eq. 0) then
        index = index + 1
        vertID(v1) = index
        do v2=v1+1,nvert
           if (vertID(v2) .eq. 0) then
              dP = verts(v1,:) - verts(v2,:)
              norm = dot_product(dP,dP)**0.5
              if (norm .lt. 1e-14) then
                 vertID(v2) = index
              end if
           end if
        end do
     end if
  end do
  count = nvert - index

end subroutine countDuplicateVerts







subroutine deleteDuplicateVerts(nvert, nvert0, nedge0, verts0, edges0, verts, edges)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nvert, nvert0, nedge0, verts0, edges0
  !f2py intent(out) verts, edges
  !f2py depend(nvert0) verts0
  !f2py depend(nedge0) edges0
  !f2py depend(nvert) verts
  !f2py depend(nedge0) edges

  !Input
  integer, intent(in) ::  nvert, nvert0, nedge0
  double precision, intent(in) ::  verts0(nvert0,2), edges0(nedge0,5)

  !Output
  double precision, intent(out) ::  verts(nvert,2), edges(nedge0,5)

  !Working
  integer v1, v2, e, d
  integer vertID(nvert0), index
  double precision dP(2), norm

  vertID(:) = 0

  index = 0
  do v1=1,nvert0
     if (vertID(v1) .eq. 0) then
        index = index + 1
        vertID(v1) = index
        verts(index,:) = verts0(v1,:)
        do v2=v1+1,nvert0
           if (vertID(v2) .eq. 0) then
              dP = verts0(v1,:) - verts0(v2,:)
              norm = dot_product(dP,dP)**0.5
              if (norm .lt. 1e-14) then
                 vertID(v2) = index
              end if
           end if
        end do
     end if
  end do

  do e=1,nedge0
     do d=1,2
        edges(e,d) = vertID(int(edges0(e,d)))
     end do
  end do

end subroutine deleteDuplicateVerts
