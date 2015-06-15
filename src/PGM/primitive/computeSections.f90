subroutine computeSections(ax1, ax2, ni, nj, nD, &
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

  Da(:) = 0.0
  Di(:) = 0
  Dj(:) = 0

  ! X_ij,k = T_j,kl * ((shp_ij,l - ogn_j,l) * scl_j,l) + pos_j,k
  ! (1) scl   (2) rot   (3) shp   (4) pos
  do j=1,nj
     call computeRotations(ax1, ax2, nj, nor, pos, &
          rot_nor, drot_dpos)
     rot_tot(:,:) = rot(:,:) * pi / 180.0 + rot_nor(:,:)
     call computeRtnMtx(ax1, ax2, ax3, rot_tot(j,:), T, dT_drot)
     do i=1,ni
        cp_array(i,j,:) = matmul(T,(shp(i,j,:)-ogn(j,:))*scl(j,:)) + pos(j,:)
        do k=1,3
           do l=1,3
              Di(iD+1:iD+6) = cp_inds(i,j,k)

              Da(iD+1) = T(k,l)*(shp(i,j,l)-ogn(j,l))
              Dj(iD+1) = scl_inds(j,l)
              Da(iD+2) = dot_product(dT_drot(k,:,l),(shp(i,j,:)-ogn(j,:))*scl(j,:)) * pi / 180.0
              Dj(iD+2) = rot_inds(j,l)
              Da(iD+3) = T(k,l)*scl(j,l)
              Dj(iD+3) = shp_inds(i,j,l)

              Dj(iD+4) = pos_inds(j,l)
              if (k .eq. l) then
                 Da(iD+4) = Da(iD+4) + 1.0
              end if
              Da(iD+4) = Da(iD+4) + dot_product(dT_drot(k,:,1),(shp(i,j,:)-ogn(j,:))*scl(j,:)) &
                   * drot_dpos(j,1,l,2)
              Da(iD+4) = Da(iD+4) + dot_product(dT_drot(k,:,2),(shp(i,j,:)-ogn(j,:))*scl(j,:)) &
                   * drot_dpos(j,2,l,2)
              Da(iD+4) = Da(iD+4) + dot_product(dT_drot(k,:,3),(shp(i,j,:)-ogn(j,:))*scl(j,:)) &
                   * drot_dpos(j,3,l,2)

              if (j .ne. 1) then
                 Dj(iD+5) = pos_inds(j-1,l)
                 Da(iD+5) = Da(iD+5) + dot_product(dT_drot(k,:,1),(shp(i,j,:)-ogn(j,:))*scl(j,:)) &
                      * drot_dpos(j,1,l,1)
                 Da(iD+5) = Da(iD+5) + dot_product(dT_drot(k,:,2),(shp(i,j,:)-ogn(j,:))*scl(j,:)) &
                      * drot_dpos(j,2,l,1)
                 Da(iD+5) = Da(iD+5) + dot_product(dT_drot(k,:,3),(shp(i,j,:)-ogn(j,:))*scl(j,:)) &
                      * drot_dpos(j,3,l,1)
              end if

              if (j .ne. nj) then
                 Dj(iD+6) = pos_inds(j+1,l)
                 Da(iD+6) = Da(iD+6) + dot_product(dT_drot(k,:,1),(shp(i,j,:)-ogn(j,:))*scl(j,:)) &
                      * drot_dpos(j,1,l,3)
                 Da(iD+6) = Da(iD+6) + dot_product(dT_drot(k,:,2),(shp(i,j,:)-ogn(j,:))*scl(j,:)) &
                      * drot_dpos(j,2,l,3)
                 Da(iD+6) = Da(iD+6) + dot_product(dT_drot(k,:,3),(shp(i,j,:)-ogn(j,:))*scl(j,:)) &
                      * drot_dpos(j,3,l,3)
              end if
              iD = iD + 6
           end do
        end do
     end do
  end do

  if (iD .ne. nD) then
     print *, 'Error in computeSections', iD, nD
  end if

end subroutine computeSections
