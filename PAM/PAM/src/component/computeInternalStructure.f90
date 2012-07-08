subroutine assembleAmtx(nA, n1, nP, nJ1, nJ2, nJQ1, nJQ2, nvert, nquad, &
     Jn1, Jn2, Ju1, Ju2, JQ1, JQ2, quad_indices, verts, poly_vert, P, Aa, Ai, Aj)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nA, n1, nP, nJ1, nJ2, nJQ1, nJQ2, nvert, nquad, Jn1, Jn2, Ju1, Ju2, JQ1, JQ2, quad_indices, verts, poly_vert, P
  !f2py intent(out) Aa, Ai, Aj
  !f2py depend(nJ1) Jn1
  !f2py depend(nJ2) Jn2
  !f2py depend(nJ1) Ju1
  !f2py depend(nJ2) Ju2
  !f2py depend(nJQ1) JQ1
  !f2py depend(nJQ2) JQ2
  !f2py depend(nquad) quad_indices
  !f2py depend(nvert) verts
  !f2py depend(nquad) poly_vert
  !f2py depend(nP) P
  !f2py depend(nA) Aa, Ai, Aj

  !Input
  integer, intent(in) ::  nA, n1, nP, nJ1, nJ2, nJQ1, nJQ2, nvert, nquad
  integer, intent(in) ::  Jn1(nJ1,4), Jn2(nJ2,4)
  double precision, intent(in) ::  Ju1(nJ1,4), Ju2(nJ2,4)
  integer, intent(in) ::  JQ1(nJQ1), JQ2(nJQ2)
  integer, intent(in) ::  quad_indices(nquad,2)
  double precision, intent(in) ::  verts(nvert,2)
  integer, intent(in) ::  poly_vert(nquad,5)
  double precision, intent(in) ::  P(nP,3)

  !Output
  double precision, intent(out) ::  Aa(nA)
  integer, intent(out) ::  Ai(nA), Aj(nA)

  !Working
  integer iA, s, iP
  integer quad, jctn
  double precision uu, vv, Q(4,2), Z(2)

  Aa(:) = 0.0
  Ai(:) = 1
  Aj(:) = 1

  iA = 0
  do iP=1,nP
     call findQuaduv(nvert, nquad, P(iP,1), P(iP,2), verts, poly_vert, quad)
     do s=1,2
        if (quad_indices(quad,s) .eq. 0) then
           if (s .eq. 1) then
              call findJunctionuv(quad, nvert, nquad, nJ1, &
                   verts, poly_vert, Ju1, jctn, Q)
           else
              call findJunctionuv(quad, nvert, nquad, nJ2, &
                   verts, poly_vert, Ju2, jctn, Q)
           end if
           call invBilinearMap(P(iP,1), P(iP,2), Q, uu, vv)
           call appendJunctionA(s, quad, nA, n1, nquad, nJ1, &
                quad_indices, Jn1, uu, vv, P(iP,3), iP, iA, Aa, Ai, Aj)
        else
           Q(1,:) = verts(poly_vert(quad,1),:)
           Q(2,:) = verts(poly_vert(quad,2),:)
           Q(3,:) = verts(poly_vert(quad,3),:)
           Q(4,:) = verts(poly_vert(quad,4),:)
           call invBilinearMap(P(iP,1), P(iP,2), Q, uu, vv)
           call appendQuadA(s, quad, nA, n1, nquad, nJ1, &
                quad_indices, Jn1, uu, vv, P(iP,3), iP, iA, Aa, Ai, Aj)
        end if
     end do
  end do

  do iA=1,nA
     Ai(iA) = Ai(iA) - 1
     Aj(iA) = Aj(iA) - 1
  end do

end subroutine assembleAmtx



subroutine countAnnz(nP, nvert, nquad, verts, poly_vert, quad_indices, P, nA)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nP, nvert, nquad, verts, poly_vert, quad_indices, P
  !f2py intent(out) nA
  !f2py depend(nvert) verts
  !f2py depend(nquad) poly_vert, quad_indices
  !f2py depend(nP) P

  !Input
  integer, intent(in) ::  nP, nvert, nquad
  double precision, intent(in) ::  verts(nvert,2)
  integer, intent(in) ::  poly_vert(nquad,5), quad_indices(nquad,2)
  double precision, intent(in) ::  P(nP,3)

  !Output
  integer, intent(out) ::  nA

  !Working
  integer iP, quad, s

  nA = 0
  do iP=1,nP
     call findQuaduv(nvert, nquad, P(iP,1), P(iP,2), verts, poly_vert, quad)
     do s=1,2
        if (quad_indices(quad,s) .eq. 0) then
           nA = nA + 12
        else
           nA = nA + 4
        end if
     end do
  end do

end subroutine countAnnz



subroutine computeInternalNodes(nP, nS, n1, nM, members, P, S)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nP, nS, n1, nM, members
  !f2py intent(out) P, S
  !f2py depend(nM) members
  !f2py depend(nP) P
  !f2py depend(nS) S

  !Input
  integer, intent(in) ::  nP, nS, n1, nM
  double precision, intent(in) ::  members(nM,33)

  !Output
  double precision, intent(out) ::  P(nP,3)
  integer, intent(out) ::  S(nS,2)

  !Working
  integer iM, mem, div, iP, iS, i
  integer shape, nmem, ndiv, n2, n
  double precision mm, dd
  double precision A(3), B(3), C(3), D(3)
  double precision A0(3), B0(3), C0(3), D0(3)
  double precision A1(3), B1(3), C1(3), D1(3)
  double precision A2(3), B2(3), C2(3), D2(3)
  double precision, allocatable, dimension(:) ::  u0, v0, u, v, w
  
  P(:,:) = 0.0

  iP = 0
  iS = 0
  do iM=1,nM
     shape = int(members(iM,2))
     nmem = int(members(iM,3))
     ndiv = int(members(iM,4))
     n2 = int(members(iM,33))
     A1 = members(iM,5:7)
     B1 = members(iM,8:10)
     C1 = members(iM,11:13)
     D1 = members(iM,14:16)
     A2 = members(iM,17:19)
     B2 = members(iM,20:22)
     C2 = members(iM,23:25)
     D2 = members(iM,26:28)
     if (nmem .eq. 1) then
        mm = 0
     else
        mm = 1.0/(nmem-1)
     end if
     dd = 1.0/ndiv
     do mem=1,nmem
        call weightedAvg((mem-1)*mm, A1, A2, A0)
        call weightedAvg((mem-1)*mm, B1, B2, B0)
        call weightedAvg((mem-1)*mm, C1, C2, C0)
        call weightedAvg((mem-1)*mm, D1, D2, D0)
        do div=1,ndiv
           call weightedAvg((div-1)*dd, A0, D0, A)
           call weightedAvg((div-1)*dd, B0, C0, B)
           call weightedAvg(div*dd, B0, C0, C)
           call weightedAvg(div*dd, A0, D0, D)
           if (shape .eq. 1) then
              n = n1*n2
              iS = iS + 1
              S(iS,:) = (/ n1 , n2 /)
           else
              n = 0
           end if
           if (n .ne. 0) then
              allocate(u0(n))
              allocate(v0(n))
              allocate(u(n))
              allocate(v(n))
              allocate(w(n))
              if (shape .eq. 1) then   
                 call getShapeRect(n, n1, n2, u0, v0)
              else
                 print *, 'Shape not found'
              end if
              call projectuvw(n, A, B, C, D, u0, v0, u, v, w)
              do i=1,n
                 iP = iP + 1
                 P(iP,1) = u(i)
                 P(iP,2) = v(i)
                 P(iP,3) = w(i)
              end do
              deallocate(u0)
              deallocate(v0)
              deallocate(u)
              deallocate(v)
              deallocate(w)
           end if
        end do
     end do
  end do

end subroutine computeInternalNodes



subroutine countInternalNodes(n1, nM, members, nP, nS)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) n1, nM, members
  !f2py intent(out) nP, nS
  !f2py depend(nM) members

  !Input
  integer, intent(in) ::  n1, nM
  double precision, intent(in) ::  members(nM,33)

  !Output
  integer, intent(out) ::  nP, nS

  !Working
  integer iM, mem, div
  integer shape, nmem, ndiv, n2
  
  nP = 0
  nS = 0
  do iM=1,nM
     shape = int(members(iM,2))
     nmem = int(members(iM,3))
     ndiv = int(members(iM,4))
     n2 = int(members(iM,33))
     do mem=1,nmem
        do div=1,ndiv
           if (shape .eq. 1) then
              nP = nP + n1*n2
              nS = nS + 1
           end if
        end do
     end do
  end do

end subroutine countInternalNodes



subroutine appendJunctionA(s, quad, nA, n1, nquad, nJ1, &
     quad_indices, Jn1, uu, vv, ww, iP, iA, Aa, Ai, Aj)

  implicit none

  !Input
  integer, intent(in) ::  s, quad, nA, n1, nquad, nJ1
  integer, intent(in) ::  quad_indices(nquad,2)
  integer, intent(in) ::  Jn1(nJ1,4)
  double precision, intent(in) ::  uu, vv, ww
  integer, intent(in) ::  iP

  !Output
  integer, intent(inout) ::  iA
  double precision, intent(inout) ::  Aa(nA)
  integer, intent(inout) ::  Ai(nA), Aj(nA)

  !Working
  double precision u1, u2, v1, v2, w
  integer i1, i2, j1, j2
  integer index

  u1 = floor(uu*(n1-1))/(n1-1)
  u2 = ceiling(uu*(n1-1))/(n1-1)
  v1 = floor(vv*(n1-1))/(n1-1)
  v2 = ceiling(vv*(n1-1))/(n1-1)

  i1 = floor(uu*(n1-1)) + 1
  i2 = ceiling(uu*(n1-1)) + 1
  j1 = floor(vv*(n1-1)) + 1
  j2 = ceiling(vv*(n1-1)) + 1

  if (s .eq. 1) then
     w = ww
  else
     w = 1-ww
  end if

  iA = iA + 12
  !call getIndexQuad(j, quad, i1, j1, &
  !     n1, nquad, nJ1, quad_indices, Jn1, index)
  !Aa(iA) = w*(u2 - uu)*(v2 - vv)*(n1-1)**2
  !Ai(iA) = iA
  !Aj(iA) = index

end subroutine appendJunctionA



subroutine appendQuadA(s, quad, nA, n1, nquad, nJ1, &
     quad_indices, Jn1, uu, vv, ww, iP, iA, Aa, Ai, Aj)

  implicit none

  !Input
  integer, intent(in) ::  s, quad, nA, n1, nquad, nJ1
  integer, intent(in) ::  quad_indices(nquad,2), Jn1(nJ1,4)
  double precision, intent(in) ::  uu, vv, ww
  integer, intent(in) ::  iP

  !Output
  integer, intent(inout) ::  iA
  double precision, intent(inout) ::  Aa(nA)
  integer, intent(inout) ::  Ai(nA), Aj(nA)

  !Working
  double precision u1, u2, v1, v2, w, den
  integer ii, jj, i1, i2, j1, j2
  integer index

  den = 1.0/(n1-1)

  if (uu .ge. (1 - 1e-14)) then
     ii = n1 - 1
  else
     ii = floor(uu*(n1-1)) + 1
  end if

  if (vv .ge. (1 - 1e-14)) then
     jj = n1 - 1
  else
     jj = floor(vv*(n1-1)) + 1
  end if
  
  i1 = ii
  i2 = ii + 1
  j1 = jj
  j2 = jj + 1

  u1 = (i1-1)*den
  u2 = (i2-1)*den
  v1 = (j1-1)*den
  v2 = (j2-1)*den

  if (s .eq. 1) then
     w = ww
  else
     w = 1-ww
  end if

  iA = iA + 1
  call getIndexQuad(s, quad, i1, j1, &
       n1, nquad, nJ1, quad_indices, Jn1, index)
  Aa(iA) = w*(u2 - uu)*(v2 - vv)*(n1-1)**2
  Ai(iA) = iP
  Aj(iA) = index

  iA = iA + 1
  call getIndexQuad(s, quad, i2, j1, &
       n1, nquad, nJ1, quad_indices, Jn1, index)
  Aa(iA) = w*(uu - u1)*(v2 - vv)*(n1-1)**2
  Ai(iA) = iP
  Aj(iA) = index

  iA = iA + 1
  call getIndexQuad(s, quad, i1, j2, &
       n1, nquad, nJ1, quad_indices, Jn1, index)
  Aa(iA) = w*(u2 - uu)*(vv - v1)*(n1-1)**2
  Ai(iA) = iP
  Aj(iA) = index

  iA = iA + 1
  call getIndexQuad(s, quad, i2, j2, &
       n1, nquad, nJ1, quad_indices, Jn1, index)
  Aa(iA) = w*(uu - u1)*(vv - v1)*(n1-1)**2
  Ai(iA) = iP
  Aj(iA) = index

end subroutine appendQuadA



subroutine getIndexQuad(s, quad, i, j, &
     n1, nquad, nJ1, quad_indices, Jn1, index)

  implicit none

  !Input
  integer, intent(in) ::  s, quad, i, j
  integer, intent(in) ::  n1, nquad, nJ1
  integer, intent(in) ::  quad_indices(nquad,2)
  integer, intent(in) ::  Jn1(nJ1,4)

  !Output
  integer, intent(out) ::  index

  if (s .eq. 1) then
     index = n1**2*(quad_indices(quad,1)-1) + i + (j-1)*n1
  else     
     index = n1**2*(maxval(quad_indices(:,1))) &
          + n1**2*(quad_indices(quad,2)-1) + i + (j-1)*n1
     if (Jn1(1,1) .ne. -1) then
        index = index + sum(Jn1(:,:))
     end if
  end if

end subroutine getIndexQuad



subroutine getIndexJunction(s, jctn, edge, i, &
     n1, nquad, nJ1, nJ2, quad_indices, Jn1, Jn2, index)

  implicit none

  !Input
  integer, intent(in) ::  s, jctn, edge, i
  integer, intent(in) ::  n1, nquad, nJ1, nJ2
  integer, intent(in) ::  quad_indices(nquad,2)
  integer, intent(in) ::  Jn1(nJ1,4), Jn2(nJ2,4)

  !Output
  integer, intent(out) ::  index

  if (s .eq. 1) then
     index = n1**2*(maxval(quad_indices(:,1))) + sum(Jn1(1:jctn-1,:)) &
          + sum(Jn1(jctn,1:edge-1)) + i
  else     
     index = n1**2*(maxval(quad_indices(:,1)) + maxval(quad_indices(:,2))) &
          + sum(Jn1(:,:)) + sum(Jn2(1:jctn-1,:)) + sum(Jn2(jctn,1:edge-1)) + i
  end if

end subroutine getIndexJunction



subroutine findJunctionuv(quad, nvert, nquad, nJ, &
     verts, poly_vert, Ju, jctn, Q)

  implicit none

  !Input
  integer, intent(in) ::  quad, nvert, nquad, nJ
  double precision, intent(in) ::  verts(nvert,2)
  integer, intent(in) ::  poly_vert(nquad,5)
  double precision, intent(in) ::  Ju(nJ,4)

  !Output
  integer, intent(out) ::  jctn
  double precision, intent(out) ::  Q(4,2)

  !Working
  integer i
  double precision P(2)

  P(:) = 0.0
  do i=1,4
     P(:) = P(:) + 0.25*verts(poly_vert(quad,i),:)
  end do

  jctn = 0
  do i=1,nJ
     if ((Ju(i,1) .le. P(1)) .and. (P(1) .le. Ju(i,2)) .and. &
          (Ju(i,3) .le. P(2)) .and. (P(2) .le. Ju(i,4))) then
        jctn = i
        exit
     end if
  end do
  if (jctn .eq. 0) then
     print *, 'Error: junction not found'
  end if

  Q(1,:) = (/ Ju(jctn,1), Ju(jctn,3) /)
  Q(2,:) = (/ Ju(jctn,1), Ju(jctn,4) /)
  Q(3,:) = (/ Ju(jctn,2), Ju(jctn,4) /)
  Q(4,:) = (/ Ju(jctn,2), Ju(jctn,3) /)

end subroutine findJunctionuv



subroutine invBilinearMap(x, y, Q, u, v)

  implicit none

  !Input
  double precision, intent(in) ::  x, y, Q(4,2)

  !Output
  double precision, intent(out) ::  u, v

  !Working
  double precision P0(2), P1(2), P2(2), P3(2), P4(2)
  double precision A, B, C, B1, B2
  double precision u1, u2
  double precision P14(2), P23(2)

  P0(:) = (/ x , y /)
  P1(:) = Q(1,:)
  P2(:) = Q(2,:)
  P3(:) = Q(3,:)
  P4(:) = Q(4,:)

  call crossproduct(P1 - P4, P3 - P2, A)
  call crossproduct(P0 - P1, P3 - P2, B1)
  call crossproduct(P0 - P2, P1 - P4, B2)
  call crossproduct(P0 - P1, P2 - P1, C)
  B = B1 + B2

  if (abs(A) .lt. 1e-14) then
     u = -C/B
  else
     u1 = (-B - (B**2 - 4*A*C)**0.5)/2.0/A
     u2 = (-B + (B**2 - 4*A*C)**0.5)/2.0/A
     if ((0 .le. u1) .and. (u1 .le. 1)) then
        u = u1
     else
        u = u2   
     end if
  end if

  P14 = P1 + (P4-P1)*u
  P23 = P2 + (P3-P2)*u
  if (abs(P23(1) - P14(1)) .gt. abs(P23(2) - P14(2))) then
     v = (P0(1) - P14(1))/(P23(1) - P14(1))
  else
     v = (P0(2) - P14(2))/(P23(2) - P14(2))
  end if

end subroutine invBilinearMap



subroutine findQuaduv(nvert, nquad, u, v, verts, poly_vert, quad)

  implicit none

  !Input
  integer, intent(in) ::  nvert, nquad
  double precision, intent(in) ::  u, v
  double precision, intent(in) ::  verts(nvert,2)
  integer, intent(in) ::  poly_vert(nquad,5)

  !Output
  integer, intent(out) ::  quad

  !Working
  integer q, k
  double precision C(5,2), E(2), D(2), P0(2), cross
  logical contained

  P0(1) = u
  P0(2) = v

  quad = 0
  do q=1,nquad
     do k=1,4
        C(k,:) = verts(poly_vert(q,k),:)
     end do
     C(5,:) = C(1,:)
     
     contained = .True.
     do k=1,4
        E = C(k+1,:) - C(k,:)
        D = P0 - C(k,:)
        call crossproduct(D,E,cross)
        if (cross .lt. 0) then
           contained = .False.
        end if
     end do

     if (contained) then
        quad = q
        exit
     end if
  end do
  if (quad .eq. 0) then
     print *, 'Error: quad not found'
  end if

end subroutine findQuaduv



subroutine crossproduct(D,E,cross)

  implicit none

  !Input
  double precision, intent(in) ::  D(2), E(2)

  !Output
  double precision, intent(out) ::  cross

  cross = D(1)*E(2) - D(2)*E(1)

end subroutine crossproduct



subroutine projectuvw(n, A, B, C, D, u0, v0, u, v, w)

  implicit none

  !Input
  integer, intent(in) ::  n
  double precision, intent(in) ::  A(3), B(3), C(n), D(n)
  double precision, intent(in) ::  u0(n), v0(n)

  !Output
  double precision, intent(out) :: u(n), v(n), w(n)

  !Working
  double precision P(3)
  integer k

  u(:) = 0.0
  v(:) = 0.0
  w(:) = 0.0
  do k=1,n
     P(:) = 0.0
     P = P + A*(1 - u0(k))*(1 - v0(k))
     P = P + B*(1 - u0(k))*(v0(k) - 0)
     P = P + C*(u0(k) - 0)*(v0(k) - 0)
     P = P + D*(u0(k) - 0)*(1 - v0(k))
     u(k) = P(1)
     v(k) = P(2)
     w(k) = P(3)
     if (u(k) .lt. 0) then
        u(k) = 0
     else if (u(k) .gt. 1) then
        u(k) = 1
     end if
     if (v(k) .lt. 0) then
        v(k) = 0
     else if (v(k) .gt. 1) then
        v(k) = 1
     end if
     if (w(k) .lt. 0) then
        w(k) = 0
     else if (w(k) .gt. 1) then
        w(k) = 1
     end if
  end do

end subroutine projectuvw



subroutine weightedAvg(r, v1, v2, v)

  implicit none

  !Input
  double precision, intent(in) ::  r, v1(3), v2(3)

  !Output
  double precision, intent(out) ::  v(3)

  v = v1*(1-r) + v2*r

end subroutine weightedAvg
