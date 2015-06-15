subroutine computeQuad2Edge(nedge, nquad, edges, quads, quad2edge)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nedge, nquad, edges, quads
  !f2py intent(out) quad2edge
  !f2py depend(nedge) edges
  !f2py depend(nquad) quads, quad2edge

  !Input
  integer, intent(in) ::  nedge, nquad, edges(nedge,2), quads(nquad,4)

  !Output
  integer, intent(out) ::  quad2edge(nquad,2,2)

  !Working
  integer iquad, iedge, sameEdge, s(2,2)

  quad2edge(:,:,:) = 0
  do iquad=1,nquad
     do iedge=1,nedge
        s(1,1) = sameEdge(quads(iquad,1:2),edges(iedge,:))
        s(1,2) = -sameEdge(quads(iquad,3:4),edges(iedge,:))
        s(2,1) = sameEdge(quads(iquad,1:4:3),edges(iedge,:))
        s(2,2) = sameEdge(quads(iquad,2:3),edges(iedge,:))
        if (s(1,1) .ne. 0) then
           quad2edge(iquad,1,1) = iedge * s(1,1)
        else if (s(1,2) .ne. 0) then
           quad2edge(iquad,1,2) = iedge * s(1,2)
        else if (s(2,1) .ne. 0) then
           quad2edge(iquad,2,1) = iedge * s(2,1)
        else if (s(2,2) .ne. 0) then
           quad2edge(iquad,2,2) = iedge * s(2,2)
        end if
     end do
  end do

end subroutine computeQuad2Edge
