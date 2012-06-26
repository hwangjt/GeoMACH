subroutine getEdges(nedge, edges0, edges)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nedge, edges0
  !f2py intent(out) edges
  !f2py depend(nedge) edges0, edges

  !Input
  integer, intent(in) ::  nedge
  double precision, intent(in) ::  edges0(nedge,4)

  !Output
  double precision, intent(out) ::  edges(nedge,5)

  !Working
  integer v, e, e1, e2, d1, d2
  double precision dP(2), norm

  edges(:,:) = 0

  v = 0
  do e1=1,nedge
     do d1=1,2
        if (edges(e1,d1) .eq. 0) then
           v = v + 1
           edges(e1,d1) = v
           do e2=e1+1,nedge
              do d2=1,2
                 if (edges(e2,d2) .eq. 0) then
                    dP = edges0(e1,2*d1-1:2*d1) - edges0(e2,2*d2-1:2*d2)
                    norm = (dP(1)**2 + dP(2)**2)**0.5
                    if (norm .lt. 1e-14) then
                       edges(e2,d2) = v
                    end if
                 end if
              end do
           end do
        end if
     end do
  end do
  do e=1,nedge
     edges(e,3) = e
     edges(e,4) = 0.0
     edges(e,5) = 1.0
  end do

end subroutine getEdges



subroutine getVerts(nvert, nedge, edges0, edges, verts)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nvert, nedge, edges0, edges
  !f2py intent(out) verts
  !f2py depend(nedge) edges0, edges
  !f2py depend(nvert) verts

  !Input
  integer, intent(in) ::  nvert, nedge
  double precision, intent(in) ::  edges0(nedge,4), edges(nedge,5)

  !Output
  double precision, intent(out) ::  verts(nvert,2)

  !Working
  integer e, d, v

  do e=1,nedge
     do d=1,2
        v = int(edges(e,d))
        verts(v,:) = edges0(e,2*d-1:2*d)
     end do
  end do

end subroutine getVerts
