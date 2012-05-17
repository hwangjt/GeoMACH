subroutine nsurfaces(filename, n)

  implicit none
  include 'cgnslib_f.h'

  !Fortran-python interface directives
  !f2py intent(in) filename
  !f2py intent(out) n
  
  !Input
  character*32, intent(in) ::  filename

  !Output
  integer, intent(out) ::  n

  !Working
  integer errorstatus
  integer cg
  integer cg_nbases, cg_base
  integer cg_nzones, cg_zone
  integer cg_nbocos, cg_boco
  character*32 cg_base_name
  integer cg_base_celldim, cg_base_physdim
  character*32 cg_zone_name
  integer cg_zone_size(3)
  character*32 cg_boco_name
  integer cg_boco_type, cg_boco_ptset_type, cg_boco_npnts
  integer cg_boco_NormalIndex(3)
  integer cg_boco_NormalListFlag, cg_boco_NormalDataType, cg_boco_ndataset

  call cg_open_f(filename, CG_MODE_READ, cg, errorstatus)
  if (errorstatus .eq. CG_ERROR) call cg_error_exit_f

  call cg_nbases_f(cg, cg_nbases, errorstatus)
  if (errorstatus .eq. CG_ERROR) call cg_error_exit_f
  if (cg_nbases .lt. 1) then
     print *, 'Error: No bases found in CGNS file'
     stop
  end if

  cg_base = 1
  call cg_base_read_f(cg, cg_base, cg_base_name, cg_base_celldim, & 
       cg_base_physdim, errorstatus)
  if (errorstatus .eq. CG_ERROR) call cg_error_exit_f
  if (cg_base_celldim .ne. 3 .or. cg_base_physdim .ne. 3) then
     print *, 'The cell and physical dimensions must be 3 and 3, respectively'
     !stop
  end if

  ! Goto Base Node
  call cg_goto_f(cg, cg_base, errorstatus, 'end')
  if (errorstatus .eq. CG_ERROR) call cg_error_exit_f

  call cg_nzones_f(cg, cg_base, cg_nzones, errorstatus)
  if (errorstatus .eq. CG_ERROR) call cg_error_exit_f

  n = 0
  do cg_zone=1, cg_nzones
     call cg_zone_read_f(cg, cg_base, cg_zone, cg_zone_name, cg_zone_size, &
          errorstatus)
     if (errorstatus .eq. CG_ERROR) call cg_error_exit_f
     call cg_nbocos_f(cg, cg_base, cg_zone, cg_nbocos, errorstatus)
     do cg_boco=1, cg_nbocos
        call cg_boco_info_f(cg, cg_base, cg_zone, cg_boco, cg_boco_name, & 
             cg_boco_type, cg_boco_ptset_type, cg_boco_npnts, & 
             cg_boco_NormalIndex, cg_boco_NormalListFlag, &
             cg_boco_NormalDataType, cg_boco_ndataset, errorstatus)
        if (cg_boco_type .eq. 21 .or. cg_boco_type .eq. 23) then
           n = n + 1
        end if
        print *, cg_zone, cg_boco_type
     end do
  end do

  call cg_close_f(cg, errorstatus)

end subroutine nsurfaces



subroutine surfacesizes(filename, n, z, sizes, normal, location)

  implicit none
  include 'cgnslib_f.h'

  !Fortran-python interface directives
  !f2py intent(in) filename, n
  !f2py intent(out) z, sizes, normal, location
  !f2py depend(n) z
  !f2py depend(n) sizes
  !f2py depend(n) normal
  !f2py depend(n) location

  !Input
  character*32, intent(in) ::  filename
  integer, intent(in) ::  n

  !Output
  integer, intent(out) ::  z(n)
  integer, intent(out) ::  sizes(n,2)
  integer, intent(out) ::  normal(n)
  logical, intent(out) ::  location(n)

  !Working
  integer errorstatus
  integer cg
  integer cg_nbases, cg_base
  integer cg_nzones, cg_zone
  integer cg_nbocos, cg_boco
  character*32 cg_base_name
  integer cg_base_celldim, cg_base_physdim
  character*32 cg_zone_name
  integer cg_zone_size(3)
  character*32 cg_boco_name
  integer cg_boco_type, cg_boco_ptset_type, cg_boco_npnts
  integer cg_boco_NormalIndex(3)
  integer cg_boco_NormalListFlag, cg_boco_NormalDataType, cg_boco_ndataset
  integer counter
  

  call cg_open_f(filename, CG_MODE_READ, cg, errorstatus)
  if (errorstatus .eq. CG_ERROR) call cg_error_exit_f

  call cg_nbases_f(cg, cg_nbases, errorstatus)
  if (errorstatus .eq. CG_ERROR) call cg_error_exit_f
  if (cg_nbases .lt. 1) then
     print *, 'Error: No bases found in CGNS file'
     stop
  end if

  cg_base = 1
  call cg_base_read_f(cg, cg_base, cg_base_name, cg_base_celldim, &
       cg_base_physdim, errorstatus)
  if (errorstatus .eq. CG_ERROR) call cg_error_exit_f
  if (cg_base_celldim .ne. 3 .or. cg_base_physdim .ne. 3) then
     print *, 'The cell and physical dimensions must be 3 and 3, respectively'
     !stop
  end if

  ! Goto Base Node
  call cg_goto_f(cg, cg_base, errorstatus, 'end')
  if (errorstatus .eq. CG_ERROR) call cg_error_exit_f

  call cg_nzones_f(cg, cg_base, cg_nzones, errorstatus)
  if (errorstatus .eq. CG_ERROR) call cg_error_exit_f

  z(:) = 0
  sizes(:,:) = 0
  normal(:) = 0
  location(:) = .True.
  counter = 1
  do cg_zone=1, cg_nzones
     call cg_zone_read_f(cg, cg_base, cg_zone, cg_zone_name, cg_zone_size, &
          errorstatus)
     if (errorstatus .eq. CG_ERROR) call cg_error_exit_f
     call cg_nbocos_f(cg, cg_base, cg_zone, cg_nbocos, errorstatus)
     do cg_boco=1, cg_nbocos
        call cg_boco_info_f(cg, cg_base, cg_zone, cg_boco, cg_boco_name, &
             cg_boco_type, cg_boco_ptset_type, cg_boco_npnts, &
             cg_boco_NormalIndex, cg_boco_NormalListFlag, &
             cg_boco_NormalDataType, cg_boco_ndataset, errorstatus)
        if (cg_boco_type .eq. 21 .or. cg_boco_type .eq. 23) then
           if (cg_boco_name(1:1)=='I') then
              sizes(counter,1) = cg_zone_size(2)
              sizes(counter,2) = cg_zone_size(3)
              normal(counter) = 1
           else if (cg_boco_name(1:1)=='J') then
              sizes(counter,1) = cg_zone_size(3)
              sizes(counter,2) = cg_zone_size(1)
              normal(counter) = 2
           else if (cg_boco_name(1:1)=='K') then
              sizes(counter,1) = cg_zone_size(1)
              sizes(counter,2) = cg_zone_size(2)
              normal(counter) = 3
           else
              print *, 'Error in boco_name - normal'
           end if
           if (cg_boco_name(2:3)=='hi') then
              location(counter) = .True.
           else if (cg_boco_name(2:3)=='lo') then
              location(counter) = .False.
           else
              print *, 'Error in boco_name - location'
           end if
           z(counter) = cg_zone
           counter = counter + 1
        end if
     end do
  end do

  call cg_close_f(cg, errorstatus)

end subroutine surfacesizes 



subroutine getsurface(filename, z, normal, location, ni, nj, P)

  implicit none
  include 'cgnslib_f.h'

  !Fortran-python interface directives
  !f2py intent(in) filename, z, normal, location, ni, nj
  !f2py intent(out) P
  !f2py depend(ni,nj) P
  
  !Input
  character*32, intent(in) ::  filename
  integer, intent(in) ::  z
  integer, intent(in) ::  normal
  logical, intent(in) ::  location
  integer, intent(in) ::  ni
  integer, intent(in) ::  nj

  !Output
  double precision, intent(out) ::  P(ni,nj,3)

  !Working
  integer errorstatus
  integer cg
  integer cg_nbases, cg_base
  integer cg_nzones, cg_zone
  integer cg_nbocos, cg_boco
  character*32 cg_base_name
  integer cg_base_celldim, cg_base_physdim
  character*32 cg_zone_name
  integer cg_zone_size(3)
  character*32 cg_boco_name
  integer cg_boco_type, cg_boco_ptset_type, cg_boco_npnts
  integer cg_boco_NormalIndex(3)
  integer cg_boco_NormalListFlag, cg_boco_NormalDataType, cg_boco_ndataset
  double precision, allocatable,dimension(:,:,:,:) :: coorX,coorY,coorZ
  integer start(3)
  integer end(3)
  integer i,j

  call cg_open_f(filename, CG_MODE_READ, cg, errorstatus)
  if (errorstatus .eq. CG_ERROR) call cg_error_exit_f

  call cg_nbases_f(cg, cg_nbases, errorstatus)
  if (errorstatus .eq. CG_ERROR) call cg_error_exit_f
  if (cg_nbases .lt. 1) then
     print *, 'Error: No bases found in CGNS file'
     stop
  end if

  cg_base = 1
  call cg_base_read_f(cg, cg_base, cg_base_name, cg_base_celldim, &
       cg_base_physdim, errorstatus)
  if (errorstatus .eq. CG_ERROR) call cg_error_exit_f
  if (cg_base_celldim .ne. 3 .or. cg_base_physdim .ne. 3) then
     print *, 'The cell and physical dimensions must be 3 and 3, respectively'
     !stop
  end if

  ! Goto Base Node
  call cg_goto_f(cg, cg_base, errorstatus, 'end')
  if (errorstatus .eq. CG_ERROR) call cg_error_exit_f

  call cg_nzones_f(cg, cg_base, cg_nzones, errorstatus)
  if (errorstatus .eq. CG_ERROR) call cg_error_exit_f

  call cg_zone_read_f(cg, cg_base, z, cg_zone_name, cg_zone_size, errorstatus)
  if (errorstatus .eq. CG_ERROR) call cg_error_exit_f
  allocate(coorX(cg_zone_size(1),cg_zone_size(2),cg_zone_size(3),3))
  allocate(coorY(cg_zone_size(1),cg_zone_size(2),cg_zone_size(3),3))
  allocate(coorZ(cg_zone_size(1),cg_zone_size(2),cg_zone_size(3),3))
  start(:) = 1
  call cg_coord_read_f(cg, cg_base, z,'CoordinateX',RealDouble,start,&
       cg_zone_size,coorX,errorstatus)
  call cg_coord_read_f(cg, cg_base, z,'CoordinateY',RealDouble,start,&
       cg_zone_size,coorY,errorstatus)
  call cg_coord_read_f(cg, cg_base, z,'CoordinateZ',RealDouble,start,&
       cg_zone_size,coorZ,errorstatus)

  if (location .eqv. .True.) then
     if (normal==1) then
        do j=1,nj
           do i=1,ni
              P(i,j,1) = coorX(cg_zone_size(1),i,j,1)
              P(i,j,2) = coorY(cg_zone_size(1),i,j,1)
              P(i,j,3) = coorZ(cg_zone_size(1),i,j,1)
           end do
        end do
     else if (normal==2) then
        do j=1,nj
           do i=1,ni
              P(i,j,1) = coorX(j,cg_zone_size(2),i,1)
              P(i,j,2) = coorY(j,cg_zone_size(2),i,1)
              P(i,j,3) = coorZ(j,cg_zone_size(2),i,1)
           end do
        end do
     else if (normal==3) then
        do j=1,nj
           do i=1,ni
              P(i,j,1) = coorX(i,j,cg_zone_size(3),1)
              P(i,j,2) = coorY(i,j,cg_zone_size(3),1)
              P(i,j,3) = coorZ(i,j,cg_zone_size(3),1)
           end do
        end do
     end if
  else if (location .eqv. .False.) then
     if (normal==1) then
        do j=1,nj
           do i=1,ni
              P(i,j,1) = coorX(1,i,j,1)
              P(i,j,2) = coorY(1,i,j,1)
              P(i,j,3) = coorZ(1,i,j,1)
           end do
        end do
     else if (normal==2) then
        do j=1,nj
           do i=1,ni
              P(i,j,1) = coorX(j,1,i,1)
              P(i,j,2) = coorY(j,1,i,1)
              P(i,j,3) = coorZ(j,1,i,1)
           end do
        end do
     else if (normal==3) then
        do j=1,nj
           do i=1,ni
              P(i,j,1) = coorX(i,j,1,1)
              P(i,j,2) = coorY(i,j,1,1)
              P(i,j,3) = coorZ(i,j,1,1)
           end do
        end do
     end if
  end if

  deallocate(coorX,coorY,coorZ)

  call cg_close_f(cg, errorstatus)

end subroutine getsurface
