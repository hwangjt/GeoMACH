subroutine evaluateProjection(nP0, ns, nD, nT, nC, nP, nsurf, nedge, ngroup, &
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
  integer surf, ugroup, vgroup, ku, kv, mu, mv, nu, nv
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
     allocate(bufferT(nu,nv,2))
     allocate(bufferP(nu,nv,3))     
     call getSurfaceT(surf, nu, nv, nT, nsurf, nedge, surf_edge, &
          surf_index_P, edge_index_P, T, bufferT)
     call getSurfaceP(surf, nP, nu, nv, 3, nsurf, nedge, nvert, surf_vert, &
          surf_edge, surf_index_P, edge_index_P, P, bufferP)
     do k=1,nP0
        mind = 1e10
        u0 = 1
        v0 = 1
        do u=1,nu!,ceiling(nu/100.0)
           do v=1,nv!,ceiling(nv/100.0)
              Pc = bufferP(u,v,:)
              d = abs(dot_product(Pc-P0(k,:),Pc-P0(k,:)))
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
        do counter=0,40
           call evaluatePoint(surf,0,0,ku,kv,mu,mv,3,nD,nC,nsurf,nedge,ngroup,nvert,x(1),x(2),& 
                surf_vert,surf_edge,edge_group,group_d,& 
                surf_index_C,edge_index_C,knot_index,C,Pc)
           call evaluatePoint(surf,1,0,ku,kv,mu,mv,3,nD,nC,nsurf,nedge,ngroup,nvert,x(1),x(2),& 
                surf_vert,surf_edge,edge_group,group_d,& 
                surf_index_C,edge_index_C,knot_index,C,Pu)
           call evaluatePoint(surf,0,1,ku,kv,mu,mv,3,nD,nC,nsurf,nedge,ngroup,nvert,x(1),x(2),& 
                surf_vert,surf_edge,edge_group,group_d,&  
                surf_index_C,edge_index_C,knot_index,C,Pv)
           call evaluatePoint(surf,2,0,ku,kv,mu,mv,3,nD,nC,nsurf,nedge,ngroup,nvert,x(1),x(2),& 
                surf_vert,surf_edge,edge_group,group_d,& 
                surf_index_C,edge_index_C,knot_index,C,Puu)
           call evaluatePoint(surf,1,1,ku,kv,mu,mv,3,nD,nC,nsurf,nedge,ngroup,nvert,x(1),x(2),& 
                surf_vert,surf_edge,edge_group,group_d,& 
                surf_index_C,edge_index_C,knot_index,C,Puv)
           call evaluatePoint(surf,0,2,ku,kv,mu,mv,3,nD,nC,nsurf,nedge,ngroup,nvert,x(1),x(2),& 
                surf_vert,surf_edge,edge_group,group_d,& 
                surf_index_C,edge_index_C,knot_index,C,Pvv)
           f = Pc - P0(k,:)
           g(1) = 2*dot_product(f,Pu)
           g(2) = 2*dot_product(f,Pv)
           H(1,1) = 2*dot_product(Pu,Pu) + 2*dot_product(f,Puu)
           H(1,2) = 2*dot_product(Pu,Pv) + 2*dot_product(f,Puv)
           H(2,2) = 2*dot_product(Pv,Pv) + 2*dot_product(f,Pvv)
           H(2,1) = H(1,2)
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
           W(2,2) = H(1,1)/det
           W(2,1) = W(1,2)
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
        call evaluatePoint(surf,0,0,ku,kv,mu,mv,3,nD,nC,nsurf,nedge,ngroup,nvert,x(1),x(2),& 
             surf_vert,surf_edge,edge_group,group_d,& 
             surf_index_C,edge_index_C,knot_index,C,Pc)
        d = abs(dot_product(Pc-P0(k,:),Pc-P0(k,:)))
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

end subroutine evaluateProjection



subroutine evaluatePjtnAlongQ(nP0, ns, nD, nT, nC, nP, nsurf, nedge, ngroup, &
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
  integer surf, ugroup, vgroup, ku, kv, mu, mv, nu, nv
  double precision, allocatable, dimension(:,:,:) ::  bufferT, bufferP
  integer k, i, s, u, v, u0, v0
  double precision x(2), dx(2), g(2), H(2,2), W(2,2), norm, det
  double precision Pc(3), Pu(3), Pv(3), Puu(3), Puv(3), Pvv(3)
  double precision f(3), fu(3), fv(3), fuu(3), fuv(3), fvv(3)
  double precision R(3)
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
     allocate(bufferT(nu,nv,2))
     allocate(bufferP(nu,nv,3)) 
     call getSurfaceT(surf, nu, nv, nT, nsurf, nedge, surf_edge, &
          surf_index_P, edge_index_P, T, bufferT)
     call getSurfaceP(surf, nP, nu, nv, 3, nsurf, nedge, nvert, surf_vert, &
          surf_edge, surf_index_P, edge_index_P, P, bufferP)
     do k=1,nP0
        R(:) = Q(k,:)
        mind = 1e10
        u0 = 1
        v0 = 1
        do u=1,nu,ceiling(nu/100.0)
           do v=1,nv,ceiling(nv/100.0)
              Pc = bufferP(u,v,:)
              f = Pc - (P0(k,:) + R*dot_product(Pc-P0(k,:),R)/dot_product(R,R))
              d = abs(dot_product(f,f))
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
        do counter=0,40
           call evaluatePoint(surf,0,0,ku,kv,mu,mv,3,nD,nC,nsurf,nedge,ngroup,nvert,x(1),x(2),& 
                surf_vert,surf_edge,edge_group,group_d,& 
                surf_index_C,edge_index_C,knot_index,C,Pc)
           call evaluatePoint(surf,1,0,ku,kv,mu,mv,3,nD,nC,nsurf,nedge,ngroup,nvert,x(1),x(2),& 
                surf_vert,surf_edge,edge_group,group_d,& 
                surf_index_C,edge_index_C,knot_index,C,Pu)
           call evaluatePoint(surf,0,1,ku,kv,mu,mv,3,nD,nC,nsurf,nedge,ngroup,nvert,x(1),x(2),& 
                surf_vert,surf_edge,edge_group,group_d,& 
                surf_index_C,edge_index_C,knot_index,C,Pv)
           call evaluatePoint(surf,2,0,ku,kv,mu,mv,3,nD,nC,nsurf,nedge,ngroup,nvert,x(1),x(2),& 
                surf_vert,surf_edge,edge_group,group_d,& 
                surf_index_C,edge_index_C,knot_index,C,Puu)
           call evaluatePoint(surf,1,1,ku,kv,mu,mv,3,nD,nC,nsurf,nedge,ngroup,nvert,x(1),x(2),& 
                surf_vert,surf_edge,edge_group,group_d,& 
                surf_index_C,edge_index_C,knot_index,C,Puv)
           call evaluatePoint(surf,0,2,ku,kv,mu,mv,3,nD,nC,nsurf,nedge,ngroup,nvert,x(1),x(2),& 
                surf_vert,surf_edge,edge_group,group_d,& 
                surf_index_C,edge_index_C,knot_index,C,Pvv)
           f = Pc - (P0(k,:) + R*dot_product(Pc-P0(k,:),R)/dot_product(R,R))
           fu = Pu - R*dot_product(Pu,R)/dot_product(R,R)
           fv = Pv - R*dot_product(Pv,R)/dot_product(R,R)
           fuu = Puu - R*dot_product(Puu,R)/dot_product(R,R)
           fuv = Puv - R*dot_product(Puv,R)/dot_product(R,R)
           fvv = Pvv - R*dot_product(Pvv,R)/dot_product(R,R)
           g(1) = 2*dot_product(f,fu)
           g(2) = 2*dot_product(f,fv)
           H(1,1) = 2*dot_product(fu,fu) + 2*dot_product(f,fuu)
           H(1,2) = 2*dot_product(fu,fv) + 2*dot_product(f,fuv)
           H(2,2) = 2*dot_product(fv,fv) + 2*dot_product(f,fvv)
           H(2,1) = H(1,2)
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
           W(2,2) = H(1,1)/det
           W(2,1) = W(1,2)
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
        call evaluatePoint(surf,0,0,ku,kv,mu,mv,3,nD,nC,nsurf,nedge,ngroup,nvert,x(1),x(2),& 
             surf_vert,surf_edge,edge_group,group_d,& 
             surf_index_C,edge_index_C,knot_index,C,Pc)
        f = Pc - (P0(k,:) + R*dot_product(Pc-P0(k,:),R)/dot_product(R,R))
        d = abs(dot_product(f,f))
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

end subroutine evaluatePjtnAlongQ



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
