subroutine computeRotations2(ax1, ax2, nj, nor, pos, rot, rot_nor, drot_dpos)

  implicit none

  !Fortran-python interface directive
  !f2py intent(in) ax1, ax2, nj, nor, pos, rot
  !f2py intent(out) rot_nor, drot_dpos
  !f2py depend(nj) nor, pos, rot, rot_nor, drot_dpos

  !Input
  integer, intent(in) ::  ax1, ax2, nj
  double precision, intent(in) ::  nor(nj,3), pos(nj,3), rot(nj,3)

  !Output
  double precision, intent(out) ::  rot_nor(nj,3), drot_dpos(nj,3,3,3)

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
     call computeAngles(ax1, ax2, t, nor(j,:), rot_nor(j,:), drotj_dt)
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

end subroutine computeRotations2



subroutine computeSections2(ax1, ax2, ni, nj, nD, &
     ogn, nor, pos, rot, scl, shX, shY, shZ, &
     ogn_inds, nor_inds, pos_inds, rot_inds, scl_inds, &
     shX_inds, shY_inds, shZ_inds, cp_inds, &
     cp_array, Da, Di, Dj)

  implicit none

  !Fortran-python interface directive
  !f2py intent(in) ax1, ax2, ni, nj, nD, ogn, nor, pos, rot, scl, shX, shY, shZ, ogn_inds, nor_inds, pos_inds, rot_inds, scl_inds, shX_inds, shY_inds, shZ_inds, cp_inds
  !f2py intent(out) cp_array, Da, Di, Dj
  !f2py depend(nj) ogn, nor, pos, rot, scl
  !f2py depend(nj) ogn_inds, nor_inds, pos_inds, rot_inds, scl_inds
  !f2py depend(ni,nj) shX, shY, shZ, shX_inds, shY_inds, shZ_inds, cp_inds, cp_array
  !f2py depend(nD) Da, Di, Dj

  !Input
  integer, intent(in) ::  ax1, ax2, ni, nj, nD
  double precision, intent(in) ::  ogn(nj,3), nor(nj,3)
  double precision, intent(in) ::  pos(nj,3), rot(nj,3), scl(nj,3)
  double precision, intent(in) ::  shX(ni,nj), shY(ni,nj), shZ(ni,nj)
  integer, intent(in) ::  ogn_inds(nj,3), nor_inds(nj,3)
  integer, intent(in) ::  pos_inds(nj,3), rot_inds(nj,3), scl_inds(nj,3)
  integer, intent(in) ::  shX_inds(ni,nj), shY_inds(ni,nj), shZ_inds(ni,nj), cp_inds(ni,nj,3)

  !Output
  double precision, intent(out) ::  cp_array(ni,nj,3), Da(nD)
  integer, intent(out) ::  Di(nD), Dj(nD)

  !Working
  integer i, j, k, l, iD, ax3
  integer shp_inds(ni,nj,3)
  double precision T(3,3), dT_drot(3,3,3), shp(ni,nj,3), pi
  double precision rot_nor(nj,3), rot_tot(nj,3)
  double precision drot_dpos(nj,3,3,3)

  pi = 2 * acos(0.0)
  shp(:,:,1) = shX(:,:)
  shp(:,:,2) = shY(:,:)
  shp(:,:,3) = shZ(:,:)
  shp_inds(:,:,1) = shX_inds(:,:)
  shp_inds(:,:,2) = shY_inds(:,:)
  shp_inds(:,:,3) = shZ_inds(:,:)

  ax3 = 6 - ax1 - ax2

  iD = 0

  ! X_ij,k = T_j,kl * ((shp_ij,l - ogn_j,l) * scl_j,l) + pos_j,k
  do j=1,nj
     call computeRotations2(ax1, ax2, nj, nor, pos, rot, &
          rot_nor, drot_dpos)
     rot_tot(:,:) = rot(:,:) * pi / 180.0 + rot_nor(:,:)
     call computeRtnMtx(ax1, ax2, ax3, rot_tot(j,:), T, dT_drot)
     do i=1,ni
        cp_array(i,j,:) = matmul(T,(shp(i,j,:)-ogn(j,:))*scl(j,:)) + pos(j,:)
        do k=1,3
           do l=1,3
              Di(iD+1:iD+4) = cp_inds(i,j,k)
              Da(iD+1) = T(k,l)*(shp(i,j,l)-ogn(j,l))
              Dj(iD+1) = scl_inds(j,l)
              Da(iD+2) = dot_product(dT_drot(k,:,l),(shp(i,j,:)-ogn(j,:))*scl(j,:)) * pi / 180.0
              Dj(iD+2) = rot_inds(j,l)
              Da(iD+3) = T(k,l)*scl(j,l)
              Dj(iD+3) = shp_inds(i,j,l)
              Da(iD+4) = 1.0
              if (j .ne. 1) then
                 Da(iD+4) = Da(iD+4) + dot_product(dT_drot(k,:,1),(shp(i,j,:)-ogn(j,:))*scl(j,:)) * drot_dpos(j-1,1,l,3)
              end if
              Da(iD+4) = Da(iD+4) + dot_product(dT_drot(k,:,1),(shp(i,j,:)-ogn(j,:))*scl(j,:)) * drot_dpos(j,1,l,2)
              if (j .ne. nj) then
                 Da(iD+4) = Da(iD+4) + dot_product(dT_drot(k,:,1),(shp(i,j,:)-ogn(j,:))*scl(j,:)) * drot_dpos(j+1,1,l,1)
              end if
              Dj(iD+4) = pos_inds(j,l)
              iD = iD + 4
           end do
        end do
     end do
  end do

  if (iD .ne. nD) then
     print *, 'Error in computeSections', iD, nD
  end if

end subroutine computeSections2



subroutine computeSections(ax1, ax2, ni, nj, nD, ishape, origin, scale0, &
     pos, rot, shape0, Q, Da, Di, Dj)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) ax1, ax2, ni, nj, nD, ishape, origin, scale0, pos, rot, shape0
  !f2py intent(out) Q, Da, Di, Dj
  !f2py depend(nj) origin, scale0, pos, rot
  !f2py depend(ni,nj) shape0, Q
  !f2py depend(nD) Da, Di, Dj

  !Input
  integer, intent(in) ::  ax1, ax2, ni, nj, nD, ishape
  double precision, intent(in) ::  origin(nj,3), scale0(nj,3)
  double precision, intent(in) ::  pos(nj,3), rot(nj,3)
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
        Q(i,j,:) = matmul(T,(shape0(i,j,:)-origin(j,:))*scale0(j,:)) + pos(j,:)
        do k=1,3
           index = ni*nj*(k-1) + ni*(j-1) + i
           do l=1,3
              Da(iD) = T(k,l)*(shape0(i,j,l)-origin(j,l))
              Di(iD) = index
              Dj(iD) = nj*(l-1) + j
              iD = iD + 1
              Da(iD) = dot_product(dT_drot(k,:,l),(shape0(i,j,:)-origin(j,:))*scale0(j,:))
              Di(iD) = index
              Dj(iD) = 6*nj + nj*(l-1) + j
              iD = iD + 1
              Da(iD) = T(k,l)*scale0(j,l)
              Di(iD) = index
              Dj(iD) = 9*nj + ishape + ni*nj*(l-1) + ni*(j-1) + i
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
