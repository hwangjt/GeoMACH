subroutine getJacobianF(nJf, nJ, nC, nsurf, nedge, ngroup, nvert, surf_vert, surf_edge, edge_group, group_m, vert_count, edge_count, surf_index_C, edge_index_C, vert_symm, edge_symm, Ja, Ji, Jj, Jfa, Jfi, Jfj)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nJf, nJ, nC, nsurf, nedge, ngroup, nvert, surf_vert, surf_edge, edge_group, group_m, vert_count, edge_count, surf_index_C, edge_index_C, Ja, Ji, Jj
  !f2py intent(out) Jfa, Jfi, Jfj
  !f2py depend(nsurf) surf_vert, surf_edge
  !f2py depend(nedge) edge_group
  !f2py depend(ngroup) group_m
  !f2py depend(nvert) vert_count
  !f2py depend(nedge) edge_count
  !f2py depend(nsurf) surf_index_C
  !f2py depend(nedge) edge_index_C
  !f2py depend(nvert) vert_symm
  !f2py depend(nedge) edge_symm
  !f2py depend(nJ) Ja, Ji, Jj
  !f2py depend(nJf) Jfa, Jfi, Jfj

  !Input
  integer, intent(in) ::  nJf, nJ, nC, nsurf, nedge, ngroup, nvert
  integer, intent(in) ::  surf_vert(nsurf,2,2), surf_edge(nsurf,2,2), edge_group(nedge), group_m(ngroup), vert_count(nvert), edge_count(nedge), surf_index_C(nsurf,2), edge_index_C(nedge,2)
  logical, intent(in) ::  vert_symm(nvert), edge_symm(nedge)
  double precision, intent(in) ::  Ja(nJ)
  integer, intent(in) ::  Ji(nJ), Jj(nJ)

  !Output
  double precision, intent(out) ::  Jfa(nJf)
  integer, intent(out) ::  Jfi(nJf), Jfj(nJf)

  !Working
  integer surf, iJ, iJf
  integer mu, mv, ugroup, vgroup
  integer c_count(nC)
  integer, allocatable, dimension(:,:) ::  mappingC
  logical interior(nJ)

  call getCcount(nC, nedge, ngroup, nvert, edge_group, group_m, vert_count, edge_count, edge_index_C, c_count)
  interior(:) = .True.
  iJf = 1
  do surf=1,nsurf
     ugroup = edge_group(abs(surf_edge(surf,1,1)))
     vgroup = edge_group(abs(surf_edge(surf,2,1)))
     mu = group_m(ugroup)
     mv = group_m(vgroup)
     allocate(mappingC(mu,mv))
     call getMapping(surf, mu, mv, nsurf, nedge, ngroup, nvert, surf_vert, surf_edge, edge_group, group_m, surf_index_C, edge_index_C, mappingC)
     do iJ=1,nJ
        if (Jj(iJ).eq.surf_vert(surf,1,1)) then
           Jfa(iJf) = Ja(iJ)/c_count(Jj(iJ))
           if (vert_symm(surf_vert(surf,1,1))) then
              Jfa(iJf) = 0
           end if
           Jfi(iJf) = Ji(iJ)
           Jfj(iJf) = mappingC(2,2)
           interior(iJ) = .False.
           iJf = iJf + 1
        else if (Jj(iJ).eq.surf_vert(surf,1,2)) then
           Jfa(iJf) = Ja(iJ)/c_count(Jj(iJ))
           if (vert_symm(surf_vert(surf,1,2))) then
              Jfa(iJf) = 0
           end if
           Jfi(iJf) = Ji(iJ)
           Jfj(iJf) = mappingC(2,mv-1)
           interior(iJ) = .False.
           iJf = iJf + 1
        else if (Jj(iJ).eq.surf_vert(surf,2,1)) then
           Jfa(iJf) = Ja(iJ)/c_count(Jj(iJ))
           if (vert_symm(surf_vert(surf,2,1))) then
              Jfa(iJf) = 0
           end if
           Jfi(iJf) = Ji(iJ)
           Jfj(iJf) = mappingC(mu-1,2)
           interior(iJ) = .False.
           iJf = iJf + 1
        else if (Jj(iJ).eq.surf_vert(surf,2,2)) then
           Jfa(iJf) = Ja(iJ)/c_count(Jj(iJ))
           if (vert_symm(surf_vert(surf,2,2))) then
              Jfa(iJf) = 0
           end if
           Jfi(iJf) = Ji(iJ)
           Jfj(iJf) = mappingC(mu-1,mv-1)
           interior(iJ) = .False.
           iJf = iJf + 1
        else if ((nvert+edge_index_C(abs(surf_edge(surf,1,1)),1).lt.Jj(iJ)).and.(Jj(iJ).le.nvert+edge_index_C(abs(surf_edge(surf,1,1)),2))) then
           Jfa(iJf) = Ja(iJ)/c_count(Jj(iJ))
           if (edge_symm(abs(surf_edge(surf,1,1)))) then
              Jfa(iJf) = 0
           end if
           Jfi(iJf) = Ji(iJ)
           if (surf_edge(surf,1,1).gt.0) then
              Jfj(iJf) = mappingC(Jj(iJ)-nvert-edge_index_C(abs(surf_edge(surf,1,1)),1)+1,2)
           else
              Jfj(iJf) = mappingC(nvert+edge_index_C(abs(surf_edge(surf,1,1)),2)-Jj(iJ)+2,2)
           end if
           interior(iJ) = .False.
           iJf = iJf + 1 
        else if ((nvert+edge_index_C(abs(surf_edge(surf,1,2)),1).lt.Jj(iJ)).and.(Jj(iJ).le.nvert+edge_index_C(abs(surf_edge(surf,1,2)),2))) then
           Jfa(iJf) = Ja(iJ)/c_count(Jj(iJ))
           if (edge_symm(abs(surf_edge(surf,1,2)))) then
              Jfa(iJf) = 0
           end if
           Jfi(iJf) = Ji(iJ)
           if (surf_edge(surf,1,2).gt.0) then
              Jfj(iJf) = mappingC(Jj(iJ)-nvert-edge_index_C(abs(surf_edge(surf,1,2)),1)+1,mv-1)
           else
              Jfj(iJf) = mappingC(nvert+edge_index_C(abs(surf_edge(surf,1,2)),2)-Jj(iJ)+2,mv-1)
           end if
           interior(iJ) = .False.
           iJf = iJf + 1   
        else if ((nvert+edge_index_C(abs(surf_edge(surf,2,1)),1).lt.Jj(iJ)).and.(Jj(iJ).le.nvert+edge_index_C(abs(surf_edge(surf,2,1)),2))) then
           Jfa(iJf) = Ja(iJ)/c_count(Jj(iJ))
           if (edge_symm(abs(surf_edge(surf,2,1)))) then
              Jfa(iJf) = 0
           end if
           Jfi(iJf) = Ji(iJ)
           if (surf_edge(surf,2,1).gt.0) then
              Jfj(iJf) = mappingC(2,Jj(iJ)-nvert-edge_index_C(abs(surf_edge(surf,2,1)),1)+1)
           else
              Jfj(iJf) = mappingC(2,nvert+edge_index_C(abs(surf_edge(surf,2,1)),2)-Jj(iJ)+2)
           end if
           interior(iJ) = .False.
           iJf = iJf + 1     
        else if ((nvert+edge_index_C(abs(surf_edge(surf,2,2)),1).lt.Jj(iJ)).and.(Jj(iJ).le.nvert+edge_index_C(abs(surf_edge(surf,2,2)),2))) then
           Jfa(iJf) = Ja(iJ)/c_count(Jj(iJ))
           if (edge_symm(abs(surf_edge(surf,2,2)))) then
              Jfa(iJf) = 0
           end if
           Jfi(iJf) = Ji(iJ)
           if (surf_edge(surf,2,2).gt.0) then
              Jfj(iJf) = mappingC(mu-1,Jj(iJ)-nvert-edge_index_C(abs(surf_edge(surf,2,2)),1)+1)
           else
              Jfj(iJf) = mappingC(mu-1,nvert+edge_index_C(abs(surf_edge(surf,2,2)),2)-Jj(iJ)+2)
           end if
           interior(iJ) = .False.
           iJf = iJf + 1                
        end if
     end do
     deallocate(mappingC)
  end do
  do iJ=1,nJ
     if (interior(iJ)) then
        Jfa(iJf) = Ja(iJ)
        Jfi(iJf) = Ji(iJ)
        Jfj(iJf) = Jj(iJ)
        iJf = iJf + 1
     end if
  end do
  do iJ=1,nJf
     Jfi(iJ) = Jfi(iJ) - 1
     Jfj(iJ) = Jfj(iJ) - edge_index_C(nedge,2) - nvert - 1
  end do

end subroutine getJacobianF



subroutine getJfnnz(nJ, nC, nedge, ngroup, nvert, edge_group, group_m, vert_count, edge_count, edge_index_C, Jj, nJf)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nJ, nC, nedge, ngroup, nvert, edge_group, group_m, vert_count, edge_count, edge_index_C, Jj
  !f2py intent(out) nJf
  !f2py depend(nedge) edge_group
  !f2py depend(ngroup) group_m
  !f2py depend(nvert) vert_count
  !f2py depend(nedge) edge_count
  !f2py depend(nedge) edge_index_C
  !f2py depend(nJ) Jj

  !Input
  integer, intent(in) ::  nJ, nC, nedge, ngroup, nvert
  integer, intent(in) ::  edge_group(nedge), group_m(ngroup), vert_count(nvert), edge_count(nedge), edge_index_C(nedge,2)
  integer, intent(in) ::  Jj(nJ)

  !Output
  integer, intent(out) :: nJf

  !Working
  integer iJ
  integer c_count(nC)

  call getCcount(nC, nedge, ngroup, nvert, edge_group, group_m, vert_count, edge_count, edge_index_C, c_count)
  nJf = 0
  do iJ=1,nJ
     nJf = nJf + c_count(Jj(iJ)+1)
  end do

end subroutine getJfnnz



subroutine getCcount(nC, nedge, ngroup, nvert, edge_group, group_m, vert_count, edge_count, edge_index_C, c_count)

  implicit none

  !Input
  integer, intent(in) ::  nC, nedge, ngroup, nvert
  integer, intent(in) ::  edge_group(nedge), group_m(ngroup), vert_count(nvert), edge_count(nedge), edge_index_C(nedge,2)

  !Output
  integer, intent(out) ::  c_count(nC)

  !Working
  integer vert, edge, u, iC

  c_count(:) = 1
  iC = 1
  do vert=1,nvert
     c_count(iC) = vert_count(vert)
     iC = iC + 1
  end do
  do edge=1,nedge
     do u=1,group_m(edge_group(edge))-2
        c_count(iC) = edge_count(edge)
        iC = iC + 1        
     end do
  end do

end subroutine getCcount
