subroutine getMnnz(nsurf, nedge, ngroup, nvert, surf_edge, edge_group, group_m,&
           surf_index, edge_index, vert_index, edge_count, surf_c1, edge_c1, nM)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nsurf, nedge, ngroup, nvert, surf_edge, edge_group, group_m, surf_index, edge_index, vert_index, edge_count, surf_c1, surf_c2
  !f2py intent(out) nM
  !f2py depend(nsurf) surf_edge
  !f2py depend(nedge) edge_group
  !f2py depend(ngroup) group_m
  !f2py depend(nsurf) surf_index
  !f2py depend(nedge) edge_index
  !f2py depend(nvert) vert_index
  !f2py depend(nedge) edge_count
  !f2py depend(nsurf) surf_c1
  !f2py depend(nedge) edge_c1
  
  !Input
  integer, intent(in) ::  nsurf, nedge, ngroup, nvert
  integer, intent(in) ::  surf_edge(nsurf,2,2), edge_group(nedge), &
                          group_m(ngroup), surf_index(nsurf,2), & 
                          edge_index(nedge,2), vert_index(nvert), & 
                          edge_count(nedge)
  logical, intent(in) ::  surf_c1(nsurf,3,3), edge_c1(nedge,2)

  !Output
  integer, intent(out) ::  nM

  !Working
  integer vert, edge, surf
  integer i,j

  nM = 0
  do vert=1,nvert
     if (vert_index(vert) .ne. 0) then
        nM = nM + 1
     end if
  end do
  do edge=1,nedge
     if (edge_index(edge,2) .ne. 0) then
        nM = nM + group_m(edge_group(edge)) - 2
     end if
  end do
  do surf=1,nsurf
     do i=1,2
        if (surf_c1(surf,2,2*i-1)) then
           nM = nM + group_m(edge_group(abs(surf_edge(surf,1,i)))) - 2
        end if
        if (surf_c1(surf,2*i-1,2)) then
           nM = nM + group_m(edge_group(abs(surf_edge(surf,2,i)))) - 2
        end if
     end do
     do i=1,3,2
        do j=1,3,2
           if (surf_c1(surf,i,j)) then
              nM = nM + 1
           end if
        end do
     end do
  end do
  do edge=1,nedge
     do i=1,2
        if (edge_c1(edge,i)) then
           nM = nM + edge_count(edge)
        end if
     end do
  end do
  nM = nM + surf_index(nsurf,2)

end subroutine getMnnz



subroutine getDOFmapping(nM, nsurf, nedge, ngroup, nvert, surf_vert, surf_edge,&
           edge_group, group_m, surf_index_C, edge_index_C, edge_index_Q, &
           vert_index_Q, edge_count, surf_c1, edge_c1, Ma, Mi, Mj)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nM, nsurf, nedge, ngroup, nvert, surf_vert, surf_edge, edge_group, group_m, surf_index_C, edge_index_C, edge_index_Q, vert_index_Q, edge_count, surf_c1, edge_c1
  !f2py intent(out) Ma, Mi, Mj
  !f2py depend(nsurf) surf_vert
  !f2py depend(nsurf) surf_edge
  !f2py depend(nedge) edge_group
  !f2py depend(ngroup) group_m
  !f2py depend(nsurf) surf_index_C
  !f2py depend(nedge) edge_index_C
  !f2py depend(nedge) edge_index_Q
  !f2py depend(nvert) vert_index_Q
  !f2py depend(nedge) edge_count
  !f2py depend(nsurf) surf_c1
  !f2py depend(nedge) edge_c1
  !f2py depend(nM) Ma, Mi, Mj
  
  !Input
  integer, intent(in) ::  nM, nsurf, nedge, ngroup, nvert
  integer, intent(in) ::  surf_vert(nsurf,2,2), surf_edge(nsurf,2,2), &
           edge_group(nedge), group_m(ngroup), surf_index_C(nsurf,2), &
           edge_index_C(nedge,2), edge_index_Q(nedge,2), vert_index_Q(nvert), &
           edge_count(nedge)
  logical, intent(in) ::  surf_c1(nsurf,3,3), edge_c1(nedge,2)

  !Output
  double precision, intent(out) ::  Ma(nM)
  integer, intent(out) ::  Mi(nM), Mj(nM)

  !Working
  integer surf, edge, vert
  integer iM, u, v, i, j
  integer ugroup, vgroup, mu, mv
  double precision edge_c1count(nedge), vert_c1count(nvert)

  Ma(:) = 0.0
  Mi(:) = 0
  Mj(:) = 0

  edge_c1count(:) = 0
  vert_c1count(:) = 0
  do surf=1,nsurf
     do i=1,2
        if (surf_c1(surf,2,2*i-1)) then
           edge_c1count(abs(surf_edge(surf,1,i))) = &
           edge_c1count(abs(surf_edge(surf,1,i))) + 1
        end if
        if (surf_c1(surf,2*i-1,2)) then
           edge_c1count(abs(surf_edge(surf,2,i))) = &
           edge_c1count(abs(surf_edge(surf,2,i))) + 1
        end if
     end do
     do i=1,2
        do j=1,2
           if (surf_c1(surf,2*i-1,2*j-1)) then
              vert_c1count(surf_vert(surf,i,j)) = & 
              vert_c1count(surf_vert(surf,i,j)) + 1
           end if
        end do
     end do
     do i=1,2
        do j=1,2
           if ((edge_c1(abs(surf_edge(surf,1,i)),j)) .and. &
              (surf_edge(surf,1,i) .gt. 0)) then
              vert_c1count(surf_vert(surf,j,i)) = &
              vert_c1count(surf_vert(surf,j,i)) + &
              1.0/edge_count(abs(surf_edge(surf,1,i)))
           end if
           if ((edge_c1(abs(surf_edge(surf,1,i)),j)) .and. & 
              (surf_edge(surf,1,i) .lt. 0)) then
              vert_c1count(surf_vert(surf,3-j,i)) = &
              vert_c1count(surf_vert(surf,3-j,i)) + &
              1.0/edge_count(abs(surf_edge(surf,1,i)))
           end if
           if ((edge_c1(abs(surf_edge(surf,2,i)),j)) .and. &
              (surf_edge(surf,2,i) .gt. 0)) then
              vert_c1count(surf_vert(surf,i,j)) = &
              vert_c1count(surf_vert(surf,i,j)) + &
              1.0/edge_count(abs(surf_edge(surf,2,i)))
           end if
           if ((edge_c1(abs(surf_edge(surf,2,i)),j)) .and. &
              (surf_edge(surf,2,i) .lt. 0)) then
              vert_c1count(surf_vert(surf,i,3-j)) = &
              vert_c1count(surf_vert(surf,i,3-j)) + &
              1.0/edge_count(abs(surf_edge(surf,2,i)))
           end if
        end do
     end do
  end do

  iM = 1
  do vert=1,nvert
     if (vert_c1count(vert) .eq. 0) then
        Ma(iM) = 1.0
        Mi(iM) = vert
        Mj(iM) = vert_index_Q(vert)
        iM = iM + 1
     end if
  end do
  do edge=1,nedge
     if (edge_c1count(edge) .eq. 0) then
        do u = 1,group_m(edge_group(edge)) - 2
           Ma(iM) = 1.0
           Mi(iM) = nvert + edge_index_C(edge,1) + u
           Mj(iM) = maxval(vert_index_Q) + edge_index_Q(edge,1) + u
           iM = iM + 1
        end do
     end if
  end do
  do surf=1,nsurf
     ugroup = edge_group(abs(surf_edge(surf,1,1)))
     vgroup = edge_group(abs(surf_edge(surf,2,1)))
     mu = group_m(ugroup)
     mv = group_m(vgroup)
     do i=1,2
        if (surf_c1(surf,2,2*i-1)) then
           edge = abs(surf_edge(surf,1,i))
           do u = 1,mu-2
              Ma(iM) = 1.0/edge_c1count(edge)
              Mi(iM) = nvert + edge_index_C(edge,1) + u
              if (surf_edge(surf,1,i) .gt. 0) then
                 Mj(iM) = maxval(vert_index_Q) + maxval(edge_index_Q(:,2)) + &
                 surf_index_C(surf,1) + u + (i-1)*(mv-3)*(mu-2)
              else
                 Mj(iM) = maxval(vert_index_Q) + maxval(edge_index_Q(:,2)) + &
                 surf_index_C(surf,1) + mu - 1 - u + (i-1)*(mv-3)*(mu-2)
              end if
              iM = iM + 1 
           end do
        end if
        if (surf_c1(surf,2*i-1,2)) then
           edge = abs(surf_edge(surf,2,i))
           do v = 1,mv - 2
              Ma(iM) = 1.0/edge_c1count(edge)
              Mi(iM) = nvert + edge_index_C(edge,1) + v
              if (surf_edge(surf,2,i) .gt. 0) then
                 Mj(iM) = maxval(vert_index_Q) + maxval(edge_index_Q(:,2)) + &
                 surf_index_C(surf,1) + 1 + (v-1)*(mu-2) + (i-1)*(mu-3)
              else
                 Mj(iM) = maxval(vert_index_Q) + maxval(edge_index_Q(:,2)) + &
                 surf_index_C(surf,1) + 1 + (mv-2-v)*(mu-2) + (i-1)*(mu-3)
              end if
              iM = iM + 1 
           end do
        end if
     end do
     do i=1,2
        do j=1,2
           if (surf_c1(surf,2*i-1,2*j-1)) then
              vert = surf_vert(surf,i,j)
              Ma(iM) = 1.0/vert_c1count(vert)
              Mi(iM) = vert
              Mj(iM) = maxval(vert_index_Q) + maxval(edge_index_Q(:,2)) + &
              surf_index_C(surf,1) + 1 + (i-1)*(mu-3) + (j-1)*(mv-3)*(mu-2)
              iM = iM + 1
           end if
        end do
     end do
     do i=1,2
        do j=1,2
           edge = surf_edge(surf,1,i)
           if ((edge_c1(abs(edge),j)) .and. (edge .gt. 0)) then
              vert = surf_vert(surf,j,i)
              Ma(iM) = 1.0/vert_c1count(vert)/edge_count(abs(edge))
              Mi(iM) = vert
              Mj(iM) = maxval(vert_index_Q) + edge_index_Q(abs(edge),1) + 1 + &
              (j-1)*(group_m(edge_group(abs(edge))) - 3)
              iM = iM + 1
           end if
           if ((edge_c1(abs(edge),j)) .and. (edge .lt. 0)) then
              vert = surf_vert(surf,3-j,i)
              Ma(iM) = 1.0/vert_c1count(vert)/edge_count(abs(edge))
              Mi(iM) = vert
              Mj(iM) = maxval(vert_index_Q) + edge_index_Q(abs(edge),1) + 1 + &
              (j-1)*(group_m(edge_group(abs(edge))) - 3)
              iM = iM + 1
           end if
           edge = surf_edge(surf,2,i)
           if ((edge_c1(abs(edge),j)) .and. (edge .gt. 0)) then
              vert = surf_vert(surf,i,j)
              Ma(iM) = 1.0/vert_c1count(vert)/edge_count(abs(edge))
              Mi(iM) = vert
              Mj(iM) = maxval(vert_index_Q) + edge_index_Q(abs(edge),1) + 1 + &
              (j-1)*(group_m(edge_group(abs(edge))) - 3)
              iM = iM + 1
           end if
           if ((edge_c1(abs(edge),j)) .and. (edge .lt. 0)) then
              vert = surf_vert(surf,i,3-j)
              Ma(iM) = 1.0/vert_c1count(vert)/edge_count(abs(edge))
              Mi(iM) = vert
              Mj(iM) = maxval(vert_index_Q) + edge_index_Q(abs(edge),1) + 1 + &
              (j-1)*(group_m(edge_group(abs(edge))) - 3)
              iM = iM + 1
           end if
        end do
     end do
  end do
  u = nvert + edge_index_C(nedge,2) + 1
  v = maxval(vert_index_Q) + maxval(edge_index_Q(:,2)) + 1
  do i = 1,nM-iM+1
     Ma(iM) = 1
     Mi(iM) = u
     Mj(iM) = v
     u = u + 1
     v = v + 1
     iM = iM + 1
  end do

  do iM=1,nM
     Mi(iM) = Mi(iM) - 1
     Mj(iM) = Mj(iM) - 1
  end do

end subroutine getDOFmapping
