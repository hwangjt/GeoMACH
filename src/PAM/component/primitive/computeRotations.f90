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
