subroutine nsurfaces2(filename, n)

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
  character*32 cg_base_name
  integer cg_base_celldim, cg_base_physdim
  character*32 cg_zone_name
  integer cg_zone_size(3)

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
  if (cg_base_celldim .ne. 2 .or. cg_base_physdim .ne. 3) then
     print *, 'The cell and physical dimensions must be 2 and 3, respectively'
     stop
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
     if ((cg_zone_name(1:6) == 'NSWall').or.(cg_zone_name(1:9) == 'EulerWall'))&
        then
        n = n + 1
     end if
  end do

  call cg_close_f(cg, errorstatus)

end subroutine nsurfaces2



subroutine surfacesizes2(filename, n, z, sizes)

  implicit none
  include 'cgnslib_f.h'

  !Fortran-python interface directives
  !f2py intent(in) filename, n
  !f2py intent(out) z, sizes
  !f2py depend(n) z
  !f2py depend(n) sizes

  !Input
  character*32, intent(in) ::  filename
  integer, intent(in) ::  n

  !Output
  integer, intent(out) ::  z(n)
  integer, intent(out) ::  sizes(n,2)

  !Working
  integer errorstatus
  integer cg
  integer cg_nbases, cg_base
  integer cg_nzones, cg_zone
  character*32 cg_base_name
  integer cg_base_celldim, cg_base_physdim
  character*32 cg_zone_name
  integer cg_zone_size(3)
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
  if (cg_base_celldim .ne. 2 .or. cg_base_physdim .ne. 3) then
     print *, 'The cell and physical dimensions must be 2 and 3, respectively'
     stop
  end if

  ! Goto Base Node
  call cg_goto_f(cg, cg_base, errorstatus, 'end')
  if (errorstatus .eq. CG_ERROR) call cg_error_exit_f

  call cg_nzones_f(cg, cg_base, cg_nzones, errorstatus)
  if (errorstatus .eq. CG_ERROR) call cg_error_exit_f

  z(:) = 0
  sizes(:,:) = 0
  counter = 1
  do cg_zone=1, cg_nzones
     call cg_zone_read_f(cg, cg_base, cg_zone, cg_zone_name, & 
          cg_zone_size, errorstatus)
     if (errorstatus .eq. CG_ERROR) call cg_error_exit_f
     if ((cg_zone_name(1:6) == 'NSWall').or.(cg_zone_name(1:9) == 'EulerWall'))&
        then
        z(counter) = cg_zone
        sizes(counter,1) = cg_zone_size(1)
        sizes(counter,2) = cg_zone_size(2)
        counter = counter + 1
     end if
  end do

  call cg_close_f(cg, errorstatus)

end subroutine surfacesizes2 



subroutine getsurface2(filename, z, ni, nj, P)

  implicit none
  include 'cgnslib_f.h'

  !Fortran-python interface directives
  !f2py intent(in) filename, z, ni, nj
  !f2py intent(out) P
  !f2py depend(ni,nj) P
  
  !Input
  character*32, intent(in) ::  filename
  integer, intent(in) ::  z
  integer, intent(in) ::  ni
  integer, intent(in) ::  nj

  !Output
  double precision, intent(out) ::  P(ni,nj,3)

  !Working
  integer errorstatus
  integer cg
  integer cg_nbases, cg_base
  integer cg_nzones, cg_zone
  character*32 cg_base_name
  integer cg_base_celldim, cg_base_physdim
  character*32 cg_zone_name
  integer cg_zone_size(3)
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
  if (cg_base_celldim .ne. 2 .or. cg_base_physdim .ne. 3) then
     print *, 'The cell and physical dimensions must be 2 and 3, respectively'
     stop
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

  P(:,:,1) = coorX(:,:,1,1)
  P(:,:,2) = coorY(:,:,1,1)
  P(:,:,3) = coorZ(:,:,1,1)

  deallocate(coorX,coorY,coorZ)

  call cg_close_f(cg, errorstatus)

end subroutine getsurface2
