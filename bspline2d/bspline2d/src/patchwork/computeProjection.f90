subroutine computeProjection(surf,nu,nv,nD,nT,nC,nP0u,nP0v,nP,nsurf,nedge,ngroup,nvert,surf_vert,surf_edge,edge_group,group_k,group_m,group_n,group_d,surf_index_P,edge_index_P,surf_index_C,edge_index_C,T,C,P0,P,minP,minu,minv)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) surf,nu,nv,nD,nT,nC,nP0u,nP0v,nP,nsurf,nedge,ngroup,nvert,surf_vert,surf_edge,edge_group,group_k,group_m,group_n,group_d,surf_index_P,edge_index_P,surf_index_C,edge_index_C,T,C,P0,P
  !f2py intent(out) minP,minu,minv
  !f2py depend(nsurf) surf_vert
  !f2py depend(nsurf) surf_edge
  !f2py depend(nedge) edge_group
  !f2py depend(ngroup) group_k, group_m, group_n
  !f2py depend(nD) group_d
  !f2py depend(nsurf) surf_index_P
  !f2py depend(nedge) edge_index_P
  !f2py depend(nsurf) surf_index_C
  !f2py depend(nedge) edge_index_C
  !f2py depend(nT) T
  !f2py depend(nC) C
  !f2py depend(nP0u,nP0v) P0
  !f2py depend(nP) P
  !f2py depend(nP0u,nP0v) minP
  !f2py depend(nP0u,nP0v) minu
  !f2py depend(nP0u,nP0v) minv

  !Input
  integer, intent(in) ::  surf,nu,nv,nD,nT,nC,nP0u,nP0v,nP,nsurf,nedge,ngroup,nvert
  integer, intent(in) ::  surf_vert(nsurf,2,2), surf_edge(nsurf,2,2), edge_group(nedge)
  integer, intent(in) ::  group_k(ngroup), group_m(ngroup), group_n(ngroup)
  double precision, intent(in) ::  group_d(nD)
  integer, intent(in) ::  surf_index_P(nsurf,2), edge_index_P(nedge,2), surf_index_C(nsurf,2), edge_index_C(nedge,2)
  double precision, intent(in) ::  T(nT), C(nC,3), P0(nP0u,nP0v,3), P(nP,3)
  
  !Output
  double precision, intent(out) ::  minP(nP0u,nP0v,3), minu(nP0u,nP0v), minv(nP0u,nP0v)

  !Working
  integer iP0u,iP0v, i
  integer counter
  integer u, v
  double precision bufferT(nu,nv,2), bufferP(nu,nv,3)
  double precision mind,d
  integer u0, v0
  double precision x(2), dx(2), g(2), H(2,2), W(2,2), norm, det
  double precision Pc(3), Pu(3), Pv(3), Puu(3), Puv(3), Pvv(3)

  call extractSurfaceT(surf, nu, nv, nT, nsurf, nedge, surf_edge, surf_index_P, edge_index_P, T, bufferT)
  call getSurfaceP(surf, nP, nu, nv, nsurf, nedge, ngroup, nvert, surf_vert, surf_edge, edge_group, group_n, surf_index_P, edge_index_P, P, bufferP)
  do iP0u=1,nP0u
     do iP0v=1,nP0v
        mind = 1e10
        do u=1,nu
           do v=1,nv
              call getNorm(bufferP(u,v,:)-P0(iP0u,iP0v,:), d)
              if (d .lt. mind) then
                 mind = d
                 u0 = u
                 v0 = v
              end if
           end do
        end do
        x(1) = bufferT(u0,v0,1)
        x(2) = bufferT(u0,v0,2)
        do counter=0,100
           call computePt(surf,0,0,nD,nC,nsurf,nedge,ngroup,nvert,x(1),x(2),surf_vert,surf_edge,edge_group,group_k,group_m,group_n,group_d,surf_index_P,edge_index_P,surf_index_C,edge_index_C,C,Pc)
           call computePt(surf,1,0,nD,nC,nsurf,nedge,ngroup,nvert,x(1),x(2),surf_vert,surf_edge,edge_group,group_k,group_m,group_n,group_d,surf_index_P,edge_index_P,surf_index_C,edge_index_C,C,Pu)
           call computePt(surf,0,1,nD,nC,nsurf,nedge,ngroup,nvert,x(1),x(2),surf_vert,surf_edge,edge_group,group_k,group_m,group_n,group_d,surf_index_P,edge_index_P,surf_index_C,edge_index_C,C,Pv)
           call computePt(surf,2,0,nD,nC,nsurf,nedge,ngroup,nvert,x(1),x(2),surf_vert,surf_edge,edge_group,group_k,group_m,group_n,group_d,surf_index_P,edge_index_P,surf_index_C,edge_index_C,C,Puu)
           call computePt(surf,1,1,nD,nC,nsurf,nedge,ngroup,nvert,x(1),x(2),surf_vert,surf_edge,edge_group,group_k,group_m,group_n,group_d,surf_index_P,edge_index_P,surf_index_C,edge_index_C,C,Puv)
           call computePt(surf,0,2,nD,nC,nsurf,nedge,ngroup,nvert,x(1),x(2),surf_vert,surf_edge,edge_group,group_k,group_m,group_n,group_d,surf_index_P,edge_index_P,surf_index_C,edge_index_C,C,Pvv)
           g(1) = 2*dot_product(Pc-P0(iP0u,iP0v,:),Pu)
           g(2) = 2*dot_product(Pc-P0(iP0u,iP0v,:),Pv)
           H(1,1) = 2*dot_product(Pu,Pu) + 2*dot_product(Pc-P0(iP0u,iP0v,:),Puu)
           H(1,2) = 2*dot_product(Pu,Pv) + 2*dot_product(Pc-P0(iP0u,iP0v,:),Puv)
           H(2,1) = 2*dot_product(Pu,Pv) + 2*dot_product(Pc-P0(iP0u,iP0v,:),Puv)
           H(2,2) = 2*dot_product(Pv,Pv) + 2*dot_product(Pc-P0(iP0u,iP0v,:),Pvv)
           do i=1,2
              if (((x(i).eq.0).and.(g(i).gt.0)).or.((x(i).eq.1).and.(g(i).lt.0))) then
                 g(i) = 0.0
                 H(1,2) = 0.0
                 H(2,1) = 0.0
                 H(i,i) = 1.0
              end if
           end do
           det = H(1,1)*H(2,2) - H(1,2)*H(2,1)
           W(1,1) = H(2,2)/det
           W(1,2) = -H(1,2)/det
           W(2,1) = -H(2,1)/det
           W(2,2) = H(1,1)/det
           norm = (g(1)**2 + g(2)**2)**0.5
           dx(1) = -dot_product(W(1,:),g)
           dx(2) = -dot_product(W(2,:),g)
           do i=1,2
              if (x(i)+dx(i).lt.0) then
                 dx(i) = -x(i)
              else if (x(i)+dx(i).gt.1) then
                 dx(i) = 1-x(i)
              end if
           end do
           !print *, counter,norm
           if ((norm.lt.1e-13).or.((dx(1)**2 + dx(2)**2)**0.5.lt.1e-13)) then
              exit
           end if
           x = x + dx
        end do
        call computePt(surf,0,0,nD,nC,nsurf,nedge,ngroup,nvert,x(1),x(2),surf_vert,surf_edge,edge_group,group_k,group_m,group_n,group_d,surf_index_P,edge_index_P,surf_index_C,edge_index_C,C,Pc)
        minP(iP0u,iP0v,:) = Pc(:)
        minu(iP0u,iP0v) = x(1)
        minv(iP0u,iP0v) = x(2)
        call getNorm(bufferP(u0,v0,:)-P0(iP0u,iP0v,:), d)
        if (d .lt. ((Pc(1)-P0(iP0u,iP0v,1))**2 + (Pc(2)-P0(iP0u,iP0v,2))**2 + (Pc(3)-P0(iP0u,iP0v,3))**2)**0.5) then
           if ((norm.gt.1e-13) .and. ((dx(1)**2 + dx(2)**2)**0.5.gt.1e-13)) then
              minP(iP0u,iP0v,:) = bufferP(u0,v0,:)
              minu(iP0u,iP0v) = bufferT(u0,v0,1)
              minv(iP0u,iP0v) = bufferT(u0,v0,2)
              print *, 'Newton-based projection failed'
           end if
        end if
     end do
  end do

end subroutine computeProjection



subroutine getOuter(a,b,C)

  implicit none

  double precision, intent(in) ::  a(2),b(2)
  double precision, intent(out) ::  C(2,2)
  integer i,j

  do i=1,2
     do j=1,2
        C(i,j) = a(i)*b(j)
     end do
  end do

end subroutine getOuter



subroutine getNorm(P,norm)

  implicit none

  double precision, intent(in) ::  P(3)
  double precision, intent(out) ::  norm

  norm = (P(1)**2 + P(2)**2 + P(3)**2)**0.5

end subroutine getNorm



subroutine computePt(surf,uder,vder,nD,nC,nsurf,nedge,ngroup,nvert,u,v,surf_vert,surf_edge,edge_group,group_k,group_m,group_n,group_d,surf_index_P,edge_index_P,surf_index_C,edge_index_C,C,P)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) surf,uder,vder,nD,nC,nsurf,nedge,ngroup,nvert,u,v,surf_vert,surf_edge,edge_group,group_k,group_m,group_n,group_d,surf_index_P,edge_index_P,surf_index_C,edge_index_C,C
  !f2py intent(out) P
  !f2py depend(nsurf) surf_vert
  !f2py depend(nsurf) surf_edge
  !f2py depend(nedge) edge_group
  !f2py depend(ngroup) group_k
  !f2py depend(ngroup) group_m
  !f2py depend(ngroup) group_n
  !f2py depend(nD) group_d
  !f2py depend(nsurf) surf_index_P
  !f2py depend(nedge) edge_index_P
  !f2py depend(nsurf) surf_index_C
  !f2py depend(nedge) edge_index_C
  !f2py depend(nC) C

  !Input
  integer, intent(in) ::  surf,uder,vder,nD,nC,nsurf,nedge,ngroup,nvert
  double precision, intent(in) ::  u,v
  integer, intent(in) ::  surf_vert(nsurf,2,2), surf_edge(nsurf,2,2), edge_group(nedge),group_k(ngroup),group_m(ngroup),group_n(ngroup)
  double precision, intent(in) ::  group_d(nD)
  integer, intent(in) ::  surf_index_P(nsurf,2), edge_index_P(nedge,2), surf_index_C(nsurf,2), edge_index_C(nedge,2)
  double precision, intent(in) ::  C(nC,3)

  !Output
  double precision, intent(out) ::  P(3)

  !Working
  integer ugroup, vgroup, ku, kv
  integer iB, nB, k
  double precision t(2)
  double precision, allocatable, dimension(:) ::  Ba
  integer, allocatable, dimension(:) ::  Bi, Bj

  ugroup = edge_group(abs(surf_edge(surf,1,1)))
  vgroup = edge_group(abs(surf_edge(surf,2,1)))
  ku = group_k(ugroup)
  kv = group_k(vgroup)
  nB = ku*kv
  t(1) = u
  t(2) = v
  allocate(Ba(nB))
  allocate(Bi(nB))
  allocate(Bj(nB))
  call computeDerivativeSingle(surf, uder, vder, nB, nD, nsurf, nedge, ngroup, nvert, t, surf_vert, surf_edge, edge_group, group_k, group_m, group_n, group_d, surf_index_P, edge_index_P, surf_index_C, edge_index_C, Ba, Bi, Bj)
  P(:) = 0.0
  do iB=1,nB
     do k=1,3
        P(k) = P(k) + Ba(iB)*C(Bj(iB)+1,k)
     end do
  end do
  deallocate(Ba)
  deallocate(Bi)
  deallocate(Bj)

end subroutine computePt



