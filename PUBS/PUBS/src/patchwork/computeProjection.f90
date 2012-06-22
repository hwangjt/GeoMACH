subroutine computeProjection(nP0, ns, nD, nT, nC, nP, nsurf, nedge, ngroup, &
     nvert, surfs, surf_vert, surf_edge, edge_group, group_k, group_m, group_n,&
     group_d, surf_index_P, edge_index_P, surf_index_C, edge_index_C, &
     knot_index, T, C, P, P0, mins, minu, minv)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nP0, ns, nD, nT, nC, nP, nsurf, nedge, ngroup, nvert, surfs, surf_vert, surf_edge, edge_group, group_k, group_m, group_n, group_d, surf_index_P, edge_index_P, surf_index_C, edge_index_C, knot_index, T, C, P, P0
  !f2py intent(out) mins, minu, minv
  !f2py depend(ns) surfs
  !f2py depend(nsurf) surf_vert
  !f2py depend(nsurf) surf_edge
  !f2py depend(nedge) edge_group
  !f2py depend(ngroup) group_k,group_m,group_n
  !f2py depend(nD) group_d
  !f2py depend(nsurf) surf_index_P
  !f2py depend(nedge) edge_index_P
  !f2py depend(nsurf) surf_index_C
  !f2py depend(nedge) edge_index_C
  !f2py depend(ngroup) knot_index
  !f2py depend(nT) T
  !f2py depend(nC) C
  !f2py depend(nP) P
  !f2py depend(nP0) P0
  !f2py depend(nP0) mins
  !f2py depend(nP0) minu
  !f2py depend(nP0) minv

  !Input
  integer, intent(in) ::  nP0, ns, nD, nT, nC, nP
  integer, intent(in) ::  nsurf, nedge, ngroup, nvert
  integer, intent(in) ::  surfs(ns), surf_vert(nsurf,2,2),&
       surf_edge(nsurf,2,2), edge_group(nedge)
  integer, intent(in) ::  group_k(ngroup),group_m(ngroup),group_n(ngroup)
  double precision, intent(in) ::  group_d(nD)
  integer, intent(in) ::  surf_index_P(nsurf,2), edge_index_P(nedge,2), &
                          surf_index_C(nsurf,2), edge_index_C(nedge,2), &
                          knot_index(ngroup,2)
  double precision, intent(in) ::  T(nT), C(nC,3), P(nP,3), P0(nP0,3)
  
  !Output
  integer, intent(out) ::  mins(nP0)
  double precision, intent(out) ::  minu(nP0), minv(nP0)

  !Working
  integer surf, ugroup, vgroup, ku, kv, mu, mv, nu, nv, nB
  double precision, allocatable, dimension(:,:,:) ::  bufferT, bufferP
  integer k, i, s, u, v, u0, v0
  double precision x(2), dx(2), g(2), H(2,2), W(2,2), norm, det
  double precision Pc(3), Pu(3), Pv(3), Puu(3), Puv(3), Pvv(3), f(3)
  integer counter
  double precision minP(nP0), d, mind

  minP(:) = 1e10
  mins(:) = 1
  minu(:) = 1
  minv(:) = 1
  do s=1,ns
     surf = surfs(s)
     ugroup = edge_group(abs(surf_edge(surf,1,1)))
     vgroup = edge_group(abs(surf_edge(surf,2,1)))
     ku = group_k(ugroup)
     kv = group_k(vgroup)
     mu = group_m(ugroup)
     mv = group_m(vgroup)
     nu = group_n(ugroup)
     nv = group_n(vgroup)
     nB = ku*kv
     allocate(bufferT(nu,nv,2))
     allocate(bufferP(nu,nv,3))     
     call extractSurfaceT(surf, nu, nv, nT, nsurf, nedge, surf_edge, &
          surf_index_P, edge_index_P, T, bufferT)
     call getSurfaceP(surf, nP, nu, nv, nsurf, nedge, nvert, surf_vert, &
          surf_edge, surf_index_P, edge_index_P, P, bufferP)
     do k=1,nP0
        mind = 1e10
        u0 = 1
        v0 = 1
        do u=1,nu
           do v=1,nv
              call getNorm(bufferP(u,v,:)-P0(k,:), d)
              if (d .lt. mind) then
                 mind = d
                 u0 = u
                 v0 = v
              end if
           end do
        end do
        if (mind .lt. minP(k)) then
           minP(k) = mind
           mins(k) = surf
           minu(k) = bufferT(u0,v0,1)
           minv(k) = bufferT(u0,v0,2)
        end if
        x(1) = bufferT(u0,v0,1)
        x(2) = bufferT(u0,v0,2)
        do counter=0,100
           call computePt(surf,0,0,ku,kv,mu,mv,nB,nD,nC,nsurf,nedge,ngroup,nvert,x(1),x(2),& 
                surf_vert,surf_edge,edge_group,group_d,& 
                surf_index_C,edge_index_C,knot_index,C,Pc)
           call computePt(surf,1,0,ku,kv,mu,mv,nB,nD,nC,nsurf,nedge,ngroup,nvert,x(1),x(2),& 
                surf_vert,surf_edge,edge_group,group_d,& 
                surf_index_C,edge_index_C,knot_index,C,Pu)
           call computePt(surf,0,1,ku,kv,mu,mv,nB,nD,nC,nsurf,nedge,ngroup,nvert,x(1),x(2),& 
                surf_vert,surf_edge,edge_group,group_d,& 
                surf_index_C,edge_index_C,knot_index,C,Pv)
           call computePt(surf,2,0,ku,kv,mu,mv,nB,nD,nC,nsurf,nedge,ngroup,nvert,x(1),x(2),& 
                surf_vert,surf_edge,edge_group,group_d,& 
                surf_index_C,edge_index_C,knot_index,C,Puu)
           call computePt(surf,1,1,ku,kv,mu,mv,nB,nD,nC,nsurf,nedge,ngroup,nvert,x(1),x(2),& 
                surf_vert,surf_edge,edge_group,group_d,& 
                surf_index_C,edge_index_C,knot_index,C,Puv)
           call computePt(surf,0,2,ku,kv,mu,mv,nB,nD,nC,nsurf,nedge,ngroup,nvert,x(1),x(2),& 
                surf_vert,surf_edge,edge_group,group_d,& 
                surf_index_C,edge_index_C,knot_index,C,Pvv)
           f = Pc - P0(k,:)
           g(1) = 2*dot_product(f,Pu)
           g(2) = 2*dot_product(f,Pv)
           H(1,1) = 2*dot_product(Pu,Pu) + 2*dot_product(f,Puu)
           H(1,2) = 2*dot_product(Pu,Pv) + 2*dot_product(f,Puv)
           H(2,1) = 2*dot_product(Pu,Pv) + 2*dot_product(f,Puv)
           H(2,2) = 2*dot_product(Pv,Pv) + 2*dot_product(f,Pvv)
           do i=1,2
              if (((x(i).eq.0).and.(g(i).gt.0)).or.((x(i).eq.1).and. & 
                 (g(i).lt.0))) then
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
           print *, counter,norm,g
           if ((norm.lt.1e-13).or.((dx(1)**2 + dx(2)**2)**0.5.lt.1e-13)) then
              exit
           end if
           x = x + dx
        end do
        call computePt(surf,0,0,ku,kv,mu,mv,nB,nD,nC,nsurf,nedge,ngroup,nvert,x(1),x(2),& 
             surf_vert,surf_edge,edge_group,group_d,& 
             surf_index_C,edge_index_C,knot_index,C,Pc)
        call getNorm(Pc(:)-P0(k,:), d)
        if (d .lt. minP(k)) then
           minP(k) = d
           mins(k) = surf
           minu(k) = x(1)
           minv(k) = x(2)
        end if
     end do
     deallocate(bufferT)
     deallocate(bufferP)
  end do

end subroutine computeProjection



subroutine computePjtnAlongQ(nP0, ns, nD, nT, nC, nP, nsurf, nedge, ngroup, &
     nvert, surfs, surf_vert, surf_edge, edge_group, group_k, group_m, group_n,&
     group_d, surf_index_P, edge_index_P, surf_index_C, edge_index_C, &
     knot_index, T, C, P, P0, Q, mins, minu, minv)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nP0, ns, nD, nT, nC, nP, nsurf, nedge, ngroup, nvert, surfs, surf_vert, surf_edge, edge_group, group_k, group_m, group_n, group_d, surf_index_P, edge_index_P, surf_index_C, edge_index_C, knot_index, T, C, P, P0, Q
  !f2py intent(out) mins, minu, minv
  !f2py depend(ns) surfs
  !f2py depend(nsurf) surf_vert
  !f2py depend(nsurf) surf_edge
  !f2py depend(nedge) edge_group
  !f2py depend(ngroup) group_k,group_m,group_n
  !f2py depend(nD) group_d
  !f2py depend(nsurf) surf_index_P
  !f2py depend(nedge) edge_index_P
  !f2py depend(nsurf) surf_index_C
  !f2py depend(nedge) edge_index_C
  !f2py depend(ngroup) knot_index
  !f2py depend(nT) T
  !f2py depend(nC) C
  !f2py depend(nP) P
  !f2py depend(nP0) P0
  !f2py depend(nP0) Q
  !f2py depend(nP0) mins
  !f2py depend(nP0) minu
  !f2py depend(nP0) minv

  !Input
  integer, intent(in) ::  nP0, ns, nD, nT, nC, nP
  integer, intent(in) ::  nsurf, nedge, ngroup, nvert
  integer, intent(in) ::  surfs(ns), surf_vert(nsurf,2,2),&
       surf_edge(nsurf,2,2), edge_group(nedge)
  integer, intent(in) ::  group_k(ngroup),group_m(ngroup),group_n(ngroup)
  double precision, intent(in) ::  group_d(nD)
  integer, intent(in) ::  surf_index_P(nsurf,2), edge_index_P(nedge,2), &
                          surf_index_C(nsurf,2), edge_index_C(nedge,2), &
                          knot_index(ngroup,2)
  double precision, intent(in) ::  T(nT), C(nC,3), P(nP,3), P0(nP0,3), Q(nP0,3)
  
  !Output
  integer, intent(out) ::  mins(nP0)
  double precision, intent(out) ::  minu(nP0), minv(nP0)

  !Working
  integer surf, ugroup, vgroup, ku, kv, mu, mv, nu, nv, nB
  double precision, allocatable, dimension(:,:,:) ::  bufferT, bufferP
  integer k, i, s, u, v, u0, v0
  double precision x(3), dx(3), g(3), H(3,3), W(3,3), norm, det
  double precision Pc(3), Pu(3), Pv(3), Puu(3), Puv(3), Pvv(3), f(3)
  integer counter
  double precision minP(nP0), d, mind

  minP(:) = 1e10
  mins(:) = 1
  minu(:) = 1
  minv(:) = 1
  do s=1,ns
     surf = surfs(s)
     ugroup = edge_group(abs(surf_edge(surf,1,1)))
     vgroup = edge_group(abs(surf_edge(surf,2,1)))
     ku = group_k(ugroup)
     kv = group_k(vgroup)
     mu = group_m(ugroup)
     mv = group_m(vgroup)
     nu = group_n(ugroup)
     nv = group_n(vgroup)
     nB = ku*kv
     allocate(bufferT(nu,nv,2))
     allocate(bufferP(nu,nv,3))     
     call extractSurfaceT(surf, nu, nv, nT, nsurf, nedge, surf_edge, &
          surf_index_P, edge_index_P, T, bufferT)
     call getSurfaceP(surf, nP, nu, nv, nsurf, nedge, nvert, surf_vert, &
          surf_edge, surf_index_P, edge_index_P, P, bufferP)
     do k=1,nP0
        mind = 1e10
        u0 = 1
        v0 = 1
        do u=1,nu
           do v=1,nv
              d = abs(dot_product(bufferP(u,v,:)-P0(k,:),Q(k,:)))
              d = d/(Q(k,1)**2 + Q(k,2)**2 + Q(k,3)**2)**0.5
              if (d .lt. mind) then
                 mind = d
                 u0 = u
                 v0 = v
              end if
           end do
        end do
        if (mind .lt. minP(k)) then
           minP(k) = mind
           mins(k) = surf
           minu(k) = bufferT(u0,v0,1)
           minv(k) = bufferT(u0,v0,2)
        end if
        x(1) = bufferT(u0,v0,1)
        x(2) = bufferT(u0,v0,2)
        x(3) = 0
        do counter=0,100
           call computePt(surf,0,0,ku,kv,mu,mv,nB,nD,nC,nsurf,nedge,ngroup,nvert,x(1),x(2),& 
                surf_vert,surf_edge,edge_group,group_d,& 
                surf_index_C,edge_index_C,knot_index,C,Pc)
           call computePt(surf,1,0,ku,kv,mu,mv,nB,nD,nC,nsurf,nedge,ngroup,nvert,x(1),x(2),& 
                surf_vert,surf_edge,edge_group,group_d,& 
                surf_index_C,edge_index_C,knot_index,C,Pu)
           call computePt(surf,0,1,ku,kv,mu,mv,nB,nD,nC,nsurf,nedge,ngroup,nvert,x(1),x(2),& 
                surf_vert,surf_edge,edge_group,group_d,& 
                surf_index_C,edge_index_C,knot_index,C,Pv)
           call computePt(surf,2,0,ku,kv,mu,mv,nB,nD,nC,nsurf,nedge,ngroup,nvert,x(1),x(2),& 
                surf_vert,surf_edge,edge_group,group_d,& 
                surf_index_C,edge_index_C,knot_index,C,Puu)
           call computePt(surf,1,1,ku,kv,mu,mv,nB,nD,nC,nsurf,nedge,ngroup,nvert,x(1),x(2),& 
                surf_vert,surf_edge,edge_group,group_d,& 
                surf_index_C,edge_index_C,knot_index,C,Puv)
           call computePt(surf,0,2,ku,kv,mu,mv,nB,nD,nC,nsurf,nedge,ngroup,nvert,x(1),x(2),& 
                surf_vert,surf_edge,edge_group,group_d,& 
                surf_index_C,edge_index_C,knot_index,C,Pvv)
           f = Pc - P0(k,:) + x(3)*Q(k,:)
           g(1) = 2*dot_product(f,Pu)
           g(2) = 2*dot_product(f,Pv)
           g(3) = 2*dot_product(f,Q(k,:))
           H(1,1) = 2*dot_product(Pu,Pu) + 2*dot_product(f,Puu)
           H(1,2) = 2*dot_product(Pu,Pv) + 2*dot_product(f,Puv)
           H(2,2) = 2*dot_product(Pv,Pv) + 2*dot_product(f,Pvv)
           H(1,3) = 2*dot_product(Pu,Q(k,:))
           H(2,3) = 2*dot_product(Pv,Q(k,:))
           H(3,3) = 2*dot_product(Q(k,:),Q(k,:))
           H(2,1) = H(1,2)
           H(3,1) = H(1,3)
           H(3,2) = H(2,3)
           do i=1,2
              !if (((x(i).lt.1e-13).and.(g(i).gt.0)).or.((x(i).gt.1.0-1e-13).and. & 
              if (((x(i).eq.0).and.(g(i).gt.0)).or.((x(i).eq.1).and. & 
                 (g(i).lt.0))) then
                 g(i) = 0.0
                 H(i,:) = 0.0
                 H(:,i) = 0.0
                 H(i,i) = 1.0
              end if
           end do
           W(1,1) = H(2,2)*H(3,3) - H(2,3)*H(3,2)
           W(1,2) = H(2,3)*H(3,1) - H(2,1)*H(3,3)
           W(1,3) = H(2,1)*H(3,2) - H(2,2)*H(3,1)
           W(2,2) = H(1,1)*H(3,3) - H(1,3)*H(3,1)
           W(2,3) = H(1,2)*H(3,1) - H(1,1)*H(3,2)
           W(3,3) = H(1,1)*H(2,2) - H(1,2)*H(2,1)
           W(2,1) = W(1,2)
           W(3,1) = W(1,3)
           W(3,2) = W(2,3)
           det = H(1,1)*W(1,1) + H(1,2)*W(1,2) + H(1,3)*W(1,3)
           W(:,:) = W(:,:)/det
           dx(1) = -dot_product(W(1,:),g)
           dx(2) = -dot_product(W(2,:),g)
           dx(3) = -dot_product(W(3,:),g)
           do i=1,2
              if (x(i)+dx(i).lt.0) then
                 dx(i) = -x(i)
              else if (x(i)+dx(i).gt.1) then
                 dx(i) = 1-x(i)
              end if
           end do
           call getNorm(g, norm)
           print *, counter,norm,g,dx,x
           if ((norm.lt.1e-13).or.((dx(1)**2 + dx(2)**2 + dx(3)**2)**0.5.lt.1e-13)) then
           !if ((norm.lt.1e-13).or.((dx(1)**2 + dx(2)**2)**0.5.lt.1e-13)) then
              exit
           end if
           x = x + dx
        end do
        call computePt(surf,0,0,ku,kv,mu,mv,nB,nD,nC,nsurf,nedge,ngroup,nvert,x(1),x(2),& 
             surf_vert,surf_edge,edge_group,group_d,& 
             surf_index_C,edge_index_C,knot_index,C,Pc)
        call getNorm(Pc(:)-P0(k,:)+x(3)*Q(k,:), d)
        if (d .lt. minP(k)) then
           minP(k) = d
           mins(k) = surf
           minu(k) = x(1)
           minv(k) = x(2)
        end if
     end do
     deallocate(bufferT)
     deallocate(bufferP)
  end do

end subroutine computePjtnAlongQ



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

  !norm = (P(1)**2 + P(2)**2 + P(3)**2)**0.5
  norm = dot_product(P,P)**0.5

end subroutine getNorm



subroutine computePt(surf,uder,vder,ku,kv,mu,mv,nB,nD,nC,nsurf,nedge,ngroup,nvert,u,v,& 
           surf_vert,surf_edge,edge_group,group_d,& 
           surf_index_C,edge_index_C,knot_index,C,P)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) surf,uder,vder,ku,kv,mu,mv,nB,nD,nC,nsurf,nedge,ngroup,nvert,u,v,surf_vert,surf_edge,edge_group,group_d,surf_index_C,edge_index_C,knot_index,C
  !f2py intent(out) P
  !f2py depend(nsurf) surf_vert
  !f2py depend(nsurf) surf_edge
  !f2py depend(nedge) edge_group
  !f2py depend(nD) group_d
  !f2py depend(nsurf) surf_index_C
  !f2py depend(nedge) edge_index_C
  !f2py depend(ngroup) knot_index
  !f2py depend(nC) C

  !Input
  integer, intent(in) ::  surf,uder,vder,ku,kv,mu,mv,nB,nD,nC,nsurf,nedge,ngroup,nvert
  double precision, intent(in) ::  u,v
  integer, intent(in) ::  surf_vert(nsurf,2,2), surf_edge(nsurf,2,2), &
                          edge_group(nedge)
  double precision, intent(in) ::  group_d(nD)
  integer, intent(in) ::  surf_index_C(nsurf,2), edge_index_C(nedge,2), &
                          knot_index(ngroup,2)
  double precision, intent(in) ::  C(nC,3)

  !Output
  double precision, intent(out) ::  P(3)

  !Working
  integer iB, k
  double precision t(2)
  double precision Ba(nB)
  integer Bi(nB), Bj(nB)

  t(1) = u
  t(2) = v
  call computeDerivativeSingle(surf, uder, vder, ku, kv, mu, mv, ku+mu, &
       kv+mv, nB, nD, nsurf, nedge, ngroup, nvert, t, surf_vert, &
       surf_edge, edge_group, group_d, surf_index_C, &
       edge_index_C, knot_index, Ba, Bi, Bj)
  P(:) = 0.0
  do iB=1,nB
     do k=1,3
        P(k) = P(k) + Ba(iB)*C(Bj(iB)+1,k)
     end do
  end do

end subroutine computePt



