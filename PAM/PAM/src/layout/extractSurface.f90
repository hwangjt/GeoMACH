subroutine extractEdges(nP, n, nu1, nu2, nv1, nv2, nvert, u1, u2, v1, v2, verts, P)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nP, n, nu1, nu2, nv1, nv2, nvert, u1, u2, v1, v2, verts
  !f2py intent(out) P
  !f2py depend(nvert) verts
  !f2py depend(nP) P

  !Input
  integer, intent(in) ::  nP, n, nu1, nu2, nv1, nv2, nvert
  double precision, intent(in) ::  u1, u2, v1, v2
  double precision, intent(in) ::  verts(nvert,2)

  !Output
  double precision, intent(out) ::  P(nP,3)

  !Working
  integer iP, i1, i2, j1, j2
  double precision Cu1(nu1+1,2), Cu2(nu2+1,2), Cv1(nv1+1,2), Cv2(nv2+1,2)
  integer i, u, v
  logical val
  double precision den(2)

  i1 = 0
  i2 = 0
  j1 = 0
  j2 = 0
  do v=1,nvert
     call isAlongEdge(u1, u1, v1, v2, verts(v,:), val)
     if (val) then
        i1 = i1 + 1
        Cu1(i1,:) = verts(v,:)
     end if 
    
     call isAlongEdge(u2, u2, v1, v2, verts(v,:), val)
     if (val) then
        i2 = i2 + 1
        Cu2(i2,:) = verts(v,:)
     end if 
    
     call isAlongEdge(u1, u2, v1, v1, verts(v,:), val)
     if (val) then
        j1 = j1 + 1
        Cv1(j1,:) = verts(v,:)
     end if 
    
     call isAlongEdge(u1, u2, v2, v2, verts(v,:), val)
     if (val) then
        j2 = j2 + 1
        Cv2(j2,:) = verts(v,:)
     end if
  end do

  call sortVerts(nu1+1, Cu1)
  call sortVerts(nu2+1, Cu2)
  call sortVerts(nv1+1, Cv1)
  call sortVerts(nv2+1, Cv2)
  
  P(:,:) = 0.0

  iP = 0
  do i=1,nu1
     den = (Cu1(i+1,:) - Cu1(i,:))/(n-1)
     do u=1,n
        iP = iP + 1
        P(iP,1:2) = Cu1(i,:) + (u-1)*den
     end do
  end do
  do i=1,nu2
     den = (Cu2(i+1,:) - Cu2(i,:))/(n-1)
     do u=1,n
        iP = iP + 1
        P(iP,1:2) = Cu2(i,:) + (u-1)*den
     end do
  end do
  do i=1,nv1
     den = (Cv1(i+1,:) - Cv1(i,:))/(n-1)
     do u=1,n
        iP = iP + 1
        P(iP,1:2) = Cv1(i,:) + (u-1)*den
     end do
  end do
  do i=1,nv2
     den = (Cv2(i+1,:) - Cv2(i,:))/(n-1)
     do u=1,n
        iP = iP + 1
        P(iP,1:2) = Cv2(i,:) + (u-1)*den
     end do
  end do

end subroutine extractEdges



subroutine sortVerts(nC, C)

  implicit none

  !Input
  integer, intent(in) ::  nC

  !Output
  double precision, intent(inout) ::  C(nC,2)

  !Working
  double precision C0(nC,2)
  integer i, k, index

  C0(:,:) = C(:,:)

  if (abs(C0(1,1) - C0(2,1)) .gt. 1e-14) then
     k = 1
  else
     k = 2
  end if

  do i=1,nC
     index = maxloc(C0(:,k),1)
     C(nC+1-i,:) = C0(index,:)
     C0(index,:) = -1.0
  end do

end subroutine sortVerts



subroutine countJunctionEdges(nJQ, nvert, nquad, u1, u2, v1, v2, JQ, verts, poly_vert, nu1, nu2, nv1, nv2)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nJQ, nvert, nquad, u1, u2, v1, v2, JQ, verts, poly_vert
  !f2py intent(out) nu1, nu2, nv1, nv2
  !f2py depend(nJQ) JQ
  !f2py depend(nvert) verts
  !f2py depend(nquad) poly_vert

  !Input
  integer, intent(in) ::  nJQ, nvert, nquad
  double precision, intent(in) ::  u1, u2, v1, v2
  integer, intent(in) ::  JQ(nJQ)
  double precision, intent(in) ::  verts(nvert,2)
  integer, intent(in) ::  poly_vert(nquad,5)

  !Output
  integer, intent(out) ::  nu1, nu2, nv1, nv2

  !Working
  integer j, q, k
  double precision C(5,2)
  logical val1, val2

  nu1 = 0
  nu2 = 0
  nv1 = 0
  nv2 = 0
  do j=1,nJQ
     q = JQ(j)
     if (q .eq. -1) then
        continue
     end if

     do k=1,4
        C(k,:) = verts(poly_vert(q,k),:)
     end do
     C(5,:) = C(1,:)

     do k=1,4
        call isAlongEdge(u1, u1, v1, v2, C(k,:), val1)
        call isAlongEdge(u1, u1, v1, v2, C(k+1,:), val2)
        if (val1 .and. val2) then
           nu1 = nu1 + 1
        end if

        call isAlongEdge(u2, u2, v1, v2, C(k,:), val1)
        call isAlongEdge(u2, u2, v1, v2, C(k+1,:), val2)
        if (val1 .and. val2) then
           nu2 = nu2 + 1
        end if

        call isAlongEdge(u1, u2, v1, v1, C(k,:), val1)
        call isAlongEdge(u1, u2, v1, v1, C(k+1,:), val2)
        if (val1 .and. val2) then
           nv1 = nv1 + 1
        end if

        call isAlongEdge(u1, u2, v2, v2, C(k,:), val1)
        call isAlongEdge(u1, u2, v2, v2, C(k+1,:), val2)
        if (val1 .and. val2) then
           nv2 = nv2 + 1
        end if
     end do
  end do

end subroutine countJunctionEdges



subroutine isAlongEdge(u1, u2, v1, v2, C, val)

  implicit none

  !Input
  double precision, intent(in) ::  u1, u2, v1, v2, C(2)

  !Output
  logical, intent(out) ::  val

  !Working
  double precision tol

  tol = 1e-14
  val = ((u1 - tol) .lt. C(1)) .and. (C(1) .lt. (u2 + tol)) &
       .and. ((v1 - tol) .lt. C(2)) .and. (C(2) .lt. (v2 + tol))

end subroutine isAlongEdge



subroutine extractFlattened(n, nJQ, nP, nvert, npoly, JQ, verts, poly_vert, P)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) n, nJQ, nP, nvert, npoly, JQ, verts, poly_vert
  !f2py intent(out) P
  !f2py depend(nJQ) JQ
  !f2py depend(nvert) verts
  !f2py depend(npoly) poly_vert
  !f2py depend(nP) P

  !Input
  integer, intent(in) ::  n, nJQ, nP, nvert, npoly
  integer, intent(in) ::  JQ(nJQ)
  double precision, intent(in) ::  verts(nvert,2)
  integer, intent(in) ::  poly_vert(npoly,5)

  !Output
  double precision, intent(out) ::  P(nP,3)

  !Working
  integer iP, u, v, poly, q
  double precision C1(2), C2(2), C3(2), C4(2)
  logical found

  P(:,:) = 0.0

  iP = 0
  do poly=1,npoly
     found = .False.
     do q=1,nJQ
        if (poly .eq. JQ(q)) then
           found = .True.
        end if
     end do
     if (.not. found) then
        C1 = verts(poly_vert(poly,1),:)
        C2 = verts(poly_vert(poly,2),:)
        C3 = verts(poly_vert(poly,3),:)
        C4 = verts(poly_vert(poly,4),:)
        do v=1,n
           do u=1,n
              iP = iP + 1
              P(iP,1:2) = (n-u)*(n-v)*C1 + (n-u)*(v-1)*C2 &
                   + (u-1)*(v-1)*C3 + (u-1)*(n-v)*C4 
              P(iP,:) = P(iP,:)/(n-1)**2
           end do
        end do
     end if
  end do

end subroutine extractFlattened



subroutine extractSurface(poly, n, nvert, npoly, verts, poly_vert, P)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) poly, n, nvert, npoly, verts, poly_vert
  !f2py intent(out) P
  !f2py depend(nvert) verts
  !f2py depend(npoly) poly_vert
  !f2py depend(n) P

  !Input
  integer, intent(in) ::  poly, n, nvert, npoly
  double precision, intent(in) ::  verts(nvert,2)
  integer, intent(in) ::  poly_vert(npoly,5)

  !Output
  double precision, intent(out) ::  P(n,n,2)

  !Working
  integer u, v
  double precision C1(2), C2(2), C3(2), C4(2)

  C1 = verts(poly_vert(poly,1),:)
  C2 = verts(poly_vert(poly,2),:)
  C3 = verts(poly_vert(poly,3),:)
  C4 = verts(poly_vert(poly,4),:)
  do u=1,n
     do v=1,n
        P(u,v,:) = (n-u)*(n-v)*C1 + (n-u)*(v-1)*C2 &
             + (u-1)*(v-1)*C3 + (u-1)*(n-v)*C4 
        P(u,v,:) = P(u,v,:)/(n-1)**2
     end do
  end do

end subroutine extractSurface
  


subroutine getQuadIndices(nJQ, npoly, JQ, quad_index)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nJQ, npoly, JQ
  !f2py intent(out) quad_index
  !f2py depend(nJQ) JQ
  !f2py depend(npoly) quad_index

  !Input
  integer, intent(in) ::  nJQ, npoly
  integer, intent(in) ::  JQ(nJQ)

  !Output
  integer, intent(out) ::  quad_index(npoly)

  !Working
  integer poly, index, q
  logical found

  quad_index(:) = 0

  index = 0
  do poly=1,npoly
     found = .False.
     do q=1,nJQ
        if (poly .eq. JQ(q)) then
           found = .True.
        end if
     end do
     if (.not. found) then
        index = index + 1
        quad_index(poly) = index
     end if
  end do

end subroutine getQuadIndices
  
  
