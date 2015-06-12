subroutine computeRotations(ax1, ax2, nj, nor, pos, rot, drot_dpos)

  implicit none

  !Fortran-python interface directive
  !f2py intent(in) ax1, ax2, nj, nor, pos
  !f2py intent(out) rot, drot_dpos
  !f2py depend(nj) nor, pos, rot, drot_dpos

  !Input
  integer, intent(in) ::  ax1, ax2, nj
  double precision, intent(in) ::  nor(nj,3), pos(nj,3)

  !Output
  double precision, intent(out) ::  rot(nj,3), drot_dpos(nj,3,3,3)

  !Working
  integer j, k, l
  double precision z, one, pi, t(3), ta(3), tb(3), ra, rb
  double precision dt_dta(3,3), dt_dtb(3,3)
  double precision dt_dpos_jp1(3,3), dt_dpos_j(3,3), dt_dpos_jm1(3,3)
  double precision I(3,3), A(3,3), B(3,3), drotj_dt(3,3)

  I(:,:) = 0.0
  do k=1,3
     I(k,k) = 1.0
  end do

  drot_dpos(:,:,:,:) = 0.0
  do j=1,nj
     if (j .eq. 1) then
        t = pos(j+1,:) - pos(j,:)
        dt_dpos_jp1(:,:) = I
        dt_dpos_j(:,:) = -I
        dt_dpos_jm1(:,:) = 0.0
     else if (j .eq. nj) then
        t = pos(j,:) - pos(j-1,:)
        dt_dpos_jp1(:,:) = 0.0
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
              drot_dpos(j,k,l,1) = dot_product(drotj_dt(k,:),dt_dpos_jm1(:,l))
           end if
           drot_dpos(j,k,l,2) = dot_product(drotj_dt(k,:),dt_dpos_j(:,l))
           if (j .ne. nj) then
              drot_dpos(j,k,l,3) = dot_product(drotj_dt(k,:),dt_dpos_jp1(:,l))
           end if
        end do
     end do
  end do

end subroutine computeRotations



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
