subroutine initializeConnectivities(nsurf, vtol, etol, P, nvert, nedge, surf_vert, & 
           surf_edge)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nsurf, vtol, etol, P
  !f2py intent(out) nvert, nedge, surf_vert, surf_edge
  !f2py depend(nsurf) P
  !f2py depend(nsurf) surf_vert
  !f2py depend(nsurf) surf_edge

  !Input
  integer, intent(in) ::  nsurf
  double precision, intent(in) ::  vtol, etol
  double precision, intent(in) ::  P(nsurf,3,3,3)

  !Output
  integer, intent(out) ::  nvert, nedge
  integer, intent(out) ::  surf_vert(nsurf,2,2)
  integer, intent(out) ::  surf_edge(nsurf,2,2)

  !Working
  integer surf_ptrs(nsurf,3,3)
  integer surf,surf1,surf2
  integer u1,v1,u2,v2
  integer du1,dv1,du2,dv2
  integer vert, edge
  double precision dP(3), norm

  surf_ptrs(:,:,:) = 0

  vert = 0
  do surf1=1,nsurf
     do v1=1,3,2
        do u1=1,3,2
           if (surf_ptrs(surf1,u1,v1) .eq. 0) then
              vert = vert + 1
              surf_ptrs(surf1,u1,v1) = vert
              do surf2=surf1+1,nsurf
                 do v2=1,3,2
                    do u2=1,3,2
                       if (surf_ptrs(surf2,u2,v2) .eq. 0) then
                          dP(:) = P(surf1,u1,v1,:) - P(surf2,u2,v2,:)
                          norm = (dP(1)**2+dP(2)**2+dP(3)**2)**0.5
                          if (norm .lt. vtol) then
                             surf_ptrs(surf2,u2,v2) = vert
                          end if
                       end if
                    end do
                 end do
              end do
           end if
        end do
     end do
  end do

  edge = 0
  do surf1=1,nsurf
     do v1=2,2
        do u1=1,3,2
           du1 = 0
           dv1 = 1
           if (surf_ptrs(surf1,u1,v1) .eq. 0) then
              edge = edge + 1
              surf_ptrs(surf1,u1,v1) = edge
              do surf2=surf1+1,nsurf
                 do v2=2,2
                    do u2=1,3,2
                       du2 = 0
                       dv2 = 1
                       call checkEdge(nsurf,edge,surf1,surf2,u1,v1,u2,v2,du1, &
                            dv1,du2,dv2,etol,P,surf_ptrs)
                    end do
                 end do
                 do v2=1,3,2
                    do u2=2,2
                       du2 = 1
                       dv2 = 0
                       call checkEdge(nsurf,edge,surf1,surf2,u1,v1,u2,v2,du1, &
                            dv1,du2,dv2,etol,P,surf_ptrs)
                    end do
                 end do
              end do
           end if
        end do
     end do
     do v1=1,3,2
        do u1=2,2
           du1 = 1
           dv1 = 0
           if (surf_ptrs(surf1,u1,v1) .eq. 0) then
              edge = edge + 1
              surf_ptrs(surf1,u1,v1) = edge
              do surf2=surf1+1,nsurf
                 do v2=2,2
                    do u2=1,3,2
                       du2 = 0
                       dv2 = 1
                       call checkEdge(nsurf,edge,surf1,surf2,u1,v1,u2,v2,du1, & 
                            dv1,du2,dv2,etol,P,surf_ptrs)
                    end do
                 end do
                 do v2=1,3,2
                    do u2=2,2
                       du2 = 1
                       dv2 = 0
                       call checkEdge(nsurf,edge,surf1,surf2,u1,v1,u2,v2,du1, & 
                            dv1,du2,dv2,etol,P,surf_ptrs)
                    end do
                 end do
              end do
           end if
        end do
     end do
  end do

  do surf=1,nsurf
     surf_vert(surf,1,1) = surf_ptrs(surf,1,1)
     surf_vert(surf,1,2) = surf_ptrs(surf,1,3)
     surf_vert(surf,2,1) = surf_ptrs(surf,3,1)
     surf_vert(surf,2,2) = surf_ptrs(surf,3,3)
     surf_edge(surf,1,1) = surf_ptrs(surf,2,1)
     surf_edge(surf,1,2) = surf_ptrs(surf,2,3)
     surf_edge(surf,2,1) = surf_ptrs(surf,1,2)
     surf_edge(surf,2,2) = surf_ptrs(surf,3,2)
  end do

  nvert = vert
  nedge = edge

end subroutine initializeConnectivities



subroutine checkEdge(nsurf, edge, surf1, surf2, u1, v1, u2, v2, du1, dv1, du2, & 
           dv2, etol, P, surf_ptrs)

  implicit none

  !Input
  integer, intent(in) ::  nsurf, edge, surf1, surf2
  integer, intent(in) ::  u1, v1, u2, v2, du1, dv1, du2, dv2
  double precision, intent (in) ::  etol, P(nsurf,3,3,3)

  !Output
  integer, intent (inout) ::  surf_ptrs(nsurf,3,3)

  !Working
  double precision dP(3), norm

  if (surf_ptrs(surf2,u2,v2) .eq. 0) then
     dP(:) = P(surf1,u1,v1,:) - P(surf2,u2,v2,:)
     norm = (dP(1)**2+dP(2)**2+dP(3)**2)**0.5
     if ((norm .lt. etol) .and. (surf_ptrs(surf1,u1-du1,v1-dv1) .eq.  & 
        surf_ptrs(surf2,u2-du2,v2-dv2)) .and. (surf_ptrs(surf1,u1+du1,v1+dv1) & 
        .eq. surf_ptrs(surf2,u2+du2,v2+dv2))) then
        surf_ptrs(surf2,u2,v2) = edge
        
     else if ((norm .lt. etol) .and. (surf_ptrs(surf1,u1-du1,v1-dv1) .eq. &
        surf_ptrs(surf2,u2+du2,v2+dv2)) .and. (surf_ptrs(surf1,u1+du1,v1+dv1) & 
        .eq. surf_ptrs(surf2,u2-du2,v2-dv2))) then
        surf_ptrs(surf2,u2,v2) = -edge
     end if
  end if

end subroutine checkEdge
