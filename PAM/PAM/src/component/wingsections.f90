subroutine eye(n, A)

  implicit none

  !Input
  integer, intent(in) ::  n
 
  !Output
  double precision, intent(out) ::  A(n,n)

  !Working
  integer i

  A(:,:) = 0.0
  do i=1,n
     A(i,i) = 1.0
  end do

end subroutine eye



subroutine outer(n, x, y, A)

  implicit none

  !Input
  integer, intent(in) ::  n
  double precision, intent(in) ::  x(n), y(n)

  !Output
  double precision, intent(out) ::  A(n,n)

  !Working
  integer i, j

  do i=1,n
     do j=1,n
        A(i,j) = x(i)*y(j)
     end do
  end do

end subroutine outer



subroutine computeRotations(ax1, ax2, nj, nD, pos, nor, rot, Da, Di, Dj)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) ax1, ax2, nj, nD, pos, nor
  !f2py intent(out) rot, Da, Di, Dj
  !f2py depend(nj) pos, nor, rot
  !f2py depend(nD) Da, Di, Dj

  !Input
  integer, intent(in) ::  ax1, ax2, nj, nD
  double precision, intent(in) ::  pos(nj,3), nor(nj,3)

  !Output
  double precision, intent(out) ::  rot(nj,3)
  double precision, intent(out) ::  Da(nD)
  integer, intent(out) ::  Di(nD), Dj(nD)

  !Working
  integer j, k, l, iD
  double precision z, one, pi, t(3), ta(3), tb(3), ra, rb
  double precision dt_dta(3,3), dt_dtb(3,3)
  double precision dt_dpos_jp1(3,3), dt_dpos_j(3,3), dt_dpos_jm1(3,3)
  double precision I(3,3), A(3,3), B(3,3), drotj_dt(3,3)

  pi = 2*acos(0.0)
  call eye(3,I)
  iD = 1
  z = 0.0
  one = 1.0
  do j=1,nj
     dt_dpos_jm1(:,:) = 0.0
     dt_dpos_j(:,:) = 0.0
     dt_dpos_jp1(:,:) = 0.0
     if (j .eq. 1) then
        t = pos(j+1,:) - pos(j,:)
        dt_dpos_jp1(:,:) = I
        dt_dpos_j(:,:) = -I
     else if (j .eq. nj) then
        t = pos(j,:) - pos(j-1,:)
        dt_dpos_j(:,:) = I
        dt_dpos_jm1(:,:) = -I
     else
        ta = pos(j,:) - pos(j-1,:)
        tb = pos(j+1,:) - pos(j,:)
        ra = dot_product(ta,ta)**0.5
        rb = dot_product(tb,tb)**0.5
        if (ra .lt. 1e-12) then
           ra = 1.0
        end if
        if (rb .lt. 1e-12) then
           rb = 1.0
        end if
        t = ta/ra + tb/rb
        call outer(3,ta,ta,A)
        call outer(3,tb,tb,B)
        dt_dta = (ra**2*I - A)/ra**3
        dt_dtb = (rb**2*I - B)/rb**3
        dt_dpos_jm1(:,:) = -dt_dta
        dt_dpos_j(:,:) = dt_dta - dt_dtb
        dt_dpos_jp1(:,:) = dt_dtb
     end if
     call computeAngles(ax1, ax2, t, nor(j,:), rot(j,:), drotj_dt)
     do k=1,3
        do l=1,3
           if (j .ne. 1) then
              Da(iD) = dot_product(drotj_dt(k,:),dt_dpos_jm1(:,l))
              Di(iD) = (k-1)*nj + j
              Dj(iD) = (l-1)*nj + j-1
              iD = iD + 1
           end if
           Da(iD) = dot_product(drotj_dt(k,:),dt_dpos_j(:,l))
           Di(iD) = (k-1)*nj + j
           Dj(iD) = (l-1)*nj + j
           iD = iD + 1
           if (j .ne. nj) then
              Da(iD) = dot_product(drotj_dt(k,:),dt_dpos_jp1(:,l))
              Di(iD) = (k-1)*nj + j
              Dj(iD) = (l-1)*nj + j+1
              iD = iD + 1
           end if
        end do
     end do
  end do

  Di(:) = Di(:) - 1
  Dj(:) = Dj(:) - 1

end subroutine computeRotations



subroutine computeAngles(ax1, ax2, t, nor, rot, drot_dt)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) ax1, ax2, t, nor
  !f2py intent(out) rot, drot_dt

  !Input
  integer, intent(in) ::  ax1, ax2
  double precision, intent(in) ::  t(3), nor(3)

  !Output
  double precision, intent(out) ::  rot(3), drot_dt(3,3)

  !Working
  integer ax3
  double precision vnorm, v(2), w(2), dv_dt(2,3), dw_dt(2,3)
  double precision p, q, dp_dv(2), dq_dw(2)

  ax3 = 6 - ax1 - ax2
  v = (/t(ax1),t(ax2)/)
  vnorm = dot_product(v,v)**0.5
  if (vnorm .lt. 1e-12) then
     vnorm = 1.0
  end if
  w = (/vnorm,-t(ax3)/)
  dv_dt(:,:) = 0.0
  dw_dt(:,:) = 0.0
  dv_dt(1,ax1) = 1.0
  dv_dt(2,ax2) = 1.0
  dw_dt(1,ax1) = t(ax1)/vnorm
  dw_dt(1,ax2) = t(ax2)/vnorm
  dw_dt(2,ax3) = -1.0
  if ((ax2-ax1 .eq. -1) .or. (ax2-ax1 .eq. 2)) then
     v(2) = -v(2)
     w(2) = -w(2)
     dv_dt(2,:) = -dv_dt(2,:)
     dw_dt(2,:) = -dw_dt(2,:)
  end if
  call arctan2pi(v, p, dp_dv)
  call arctan2pi(w, q, dq_dw)
  rot(1) = p*nor(1)
  rot(2) = q*nor(2)
  rot(3) = 0.0
  drot_dt(1,:) = matmul(dp_dv, dv_dt)*nor(1)
  drot_dt(2,:) = matmul(dq_dw, dw_dt)*nor(2)
  drot_dt(3,:) = 0.0

end subroutine computeAngles



subroutine computeSections(ax1, ax2, f, ni, nj, nD, ishape, r, offset, chord, &
     pos, rot, shape0, Q, Da, Di, Dj)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) ax1, ax2, f, ni, nj, nD, ishape, r, offset, chord, pos, rot, shape0
  !f2py intent(out) Q, Da, Di, Dj
  !f2py depend(nj) chord, pos, rot
  !f2py depend(ni,nj) shape0, Q
  !f2py depend(nD) Da, Di, Dj

  !Input
  integer, intent(in) ::  ax1, ax2, f, ni, nj, nD, ishape
  double precision, intent(in) ::  r(3), offset(3)
  double precision, intent(in) ::  chord(nj), pos(nj,3), rot(nj,3)
  double precision, intent(in) ::  shape0(ni,nj,3)

  !Output
  double precision, intent(out) ::  Q(ni,nj,3), Da(nD)
  integer, intent(out) ::  Di(nD), Dj(nD)

  !Working
  integer i, j, k, l, iD, index, ax3
  double precision T(3,3), dT_drot(3,3,3)

  ax3 = 6 - ax1 - ax2
  iD = 1
  do j=1,nj
     call computeRtnMtx(ax1, ax2, ax3, rot(j,:), T, dT_drot)
     do i=1,ni
        Q(i,j,:) = (matmul(T,shape0(i,j,:)-r) + r)*chord(j) + pos(j,:) + offset
        do k=1,3
           index = ni*nj*(k-1) + ni*(j-1) + i
           Da(iD) = dot_product(T(k,:),shape0(i,j,:)-r) + r(k)
           Di(iD) = index
           Dj(iD) = j
           iD = iD + 1
           do l=1,3
              Da(iD) = dot_product(dT_drot(k,:,l),shape0(i,j,:)-r)*chord(j)
              Di(iD) = index
              Dj(iD) = nj + nj*(l-1) + j
              iD = iD + 1
              Da(iD) = T(k,l)*chord(j)
              Di(iD) = index
              Dj(iD) = 4*nj + ishape + ni*nj*(l-1) + ni*(j-1) + i
              iD = iD + 1
           end do
        end do
     end do
  end do

  Di(:) = Di(:) - 1
  Dj(:) = Dj(:) - 1

end subroutine computeSections



subroutine computeRtnMtx(ax1, ax2, ax3, rot, T, dT_drot)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) ax1, ax2, ax3, rot
  !f2py intent(out) T, dT_drot

  !Input
  integer, intent(in) ::  ax1, ax2, ax3
  double precision, intent(in) ::  rot(3)

  !Output
  double precision, intent(out) ::  T(3,3), dT_drot(3,3,3)

  !Working
  double precision T1(3,3), T2(3,3), T3(3,3)
  double precision dT1(3,3), dT2(3,3), dT3(3,3)

  call computeRtn(ax3, rot(1), T1, dT1)
  call computeRtn(ax2, rot(2), T2, dT2)
  call computeRtn(ax1, rot(3), T3, dT3)

  T = matmul(matmul(T1,T2),T3)
  dT_drot(:,:,1) = matmul(matmul(dT1,T2),T3)
  dT_drot(:,:,2) = matmul(matmul(T1,dT2),T3)
  dT_drot(:,:,3) = matmul(matmul(T1,T2),dT3)

end subroutine computeRtnMtx



subroutine computeRtn(k, p, T, dT)

  implicit none

  !Input
  integer, intent(in) ::  k
  double precision, intent(in) ::  p

  !Output
  double precision, intent(out) ::  T(3,3), dT(3,3)

  !Working
  integer i, j

  i = k + 1
  j = k + 2
  if (i .gt. 3) then
     i = i - 3
  end if
  if (j .gt. 3) then
     j = j - 3
  end if

  T(:,:) = 0.0
  T(k,k) = 1.0
  T(i,i) = cos(p)
  T(i,j) = -sin(p)
  T(j,i) = sin(p)
  T(j,j) = cos(p)

  dT(:,:) = 0.0
  dT(i,i) = -sin(p)
  dT(i,j) = -cos(p)
  dT(j,i) = cos(p)
  dT(j,j) = -sin(p)

end subroutine computeRtn



subroutine arctan2pi(P, t, dt_dP)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) P
  !f2py intent(out) t, dt_dP

  !Input
  double precision, intent(in) ::  P(2)

  !Output
  double precision, intent(out) ::  t, dt_dP(2)

  !Working
  double precision x, y, pi

  pi = 2*acos(0.0)

  x = P(1)
  y = P(2)

  if (dot_product(P,P)**0.5 .lt. 1e-10) then
     t = 0
  else if (x .eq. 0) then
     if (y .gt. 0) then
        t = pi/2.0
     else if (y .lt. 0) then
        t = 3*pi/2.0
     end if
  else if (y .eq. 0) then
     if (x .gt. 0) then
        t = 0
     else if (x .lt. 0) then
        t = pi
     end if
  else if (x .lt. 0) then
     t = atan(y/x) + pi   
  else if (y .lt. 0) then
     t = atan(y/x) + 2*pi    
  else if (y .gt. 0) then
     t = atan(y/x)
  else
     t = 0
  end if

  if (dot_product(P,P)**0.5 .lt. 1e-10) then
     dt_dP(:) = 0.0
  else
     dt_dP(1) = -y/(x**2+y**2)
     dt_dP(2) = x/(x**2+y**2)
  end if

end subroutine arctan2pi
