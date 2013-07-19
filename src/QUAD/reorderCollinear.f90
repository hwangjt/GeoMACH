subroutine reorderCollinear(nvert, nedge, verts0, edges0, verts, edges)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nvert, nedge, verts0, edges0
  !f2py intent(out) verts, edges
  !f2py depend(nvert) verts0, verts
  !f2py depend(nedge) edges0, edges

  !Input
  integer, intent(in) ::  nvert, nedge
  double precision, intent(in) ::  verts0(nvert,2)
  integer, intent(in) ::  edges0(nedge,2)

  !Output
  double precision, intent(out) ::  verts(nvert,2)
  integer, intent(out) ::  edges(nedge,2)

  !Working
  logical collinear
  integer ivert, iedge, index, j
  double precision temp(2)

  verts(:,:) = verts0(:,:)
  edges(:,:) = edges0(:,:)
  
  if (collinear(3, nvert, verts)) then
     index = 0
     findVert: do ivert=3,nvert
        if (.not. collinear(ivert, nvert, verts)) then
           index = ivert
           exit findVert
        end if
     end do findVert
     if (index .eq. 0) then
        print *, 'Error: all points collinear'
     else
        temp = verts(3,:)
        verts(3,:) = verts(index,:)
        verts(index,:) = temp
        checkEdges: do iedge=1,nedge
           do j=1,2
              if (edges(iedge,j) .eq. 3) then
                 edges(iedge,j) = index
              else if (edges(iedge,j) .eq. index) then
                 edges(iedge,j) = 3
              end if
           end do
        end do checkEdges
     end if
  end if

end subroutine reorderCollinear




function collinear(i, nvert, verts)

  implicit none

  !Input
  integer, intent(in) ::  i, nvert
  double precision, intent(in) ::  verts(nvert,2)

  !Output
  logical collinear

  !Working
  double precision v1(2), v2(2), det

  v1 = verts(2,:) - verts(1,:)
  v2 = verts(i,:) - verts(2,:)
  det = abs(v1(1)*v2(2) - v1(2)*v2(1))

  if (det .lt. 1e-10) then
     collinear = .True.
  else
     collinear = .False.
  end if

end function collinear
