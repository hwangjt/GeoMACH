subroutine splitEdges(nvert, nedge0, nedge, &
     verts, edges0, edgeCon0, edges, edgeCon)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nvert, nedge0, nedge, verts, edges0, edgeCon0
  !f2py intent(out) edges, edgeCon
  !f2py depend(nvert) verts
  !f2py depend(nedge0) edges0, edgeCon0
  !f2py depend(nedge) edges, edgeCon

  !Input
  integer, intent(in) ::  nvert, nedge0, nedge
  double precision, intent(in) ::  verts(nvert,2)
  integer, intent(in) ::  edges0(nedge0,2)
  logical, intent(in) ::  edgeCon0(nedge0)

  !Output
  integer, intent(out) ::  edges(nedge,2)
  logical, intent(out) ::  edgeCon(nedge)

  !Working
  integer iedge0, iedge, i, i1, i2, isplit, nsplit
  logical validSplit
  double precision split
  double precision, allocatable, dimension(:) ::  t
  integer, allocatable, dimension(:) ::  ti, lo, hi

  iedge = 1
  do iedge0=1,nedge0
     i1 = edges0(iedge0,1)
     i2 = edges0(iedge0,2)
     nsplit = 0
     do i=1,nvert
        if ((i.ne.i1) .and. (i.ne.i2)) then
           if (validSplit(verts(i1,:),verts(i2,:),verts(i,:))) then
              nsplit = nsplit + 1
           end if
        end if
     end do
     if (nsplit .eq. 0) then
        edges(iedge,:) = edges0(iedge0,:)
        edgeCon(iedge) = edgeCon0(iedge0)
        iedge = iedge + 1
     else
        allocate(t(nsplit))
        allocate(ti(nsplit))
        allocate(lo(nsplit+1))
        allocate(hi(nsplit+1))
        isplit = 1
        do i=1,nvert
           if ((i.ne.i1) .and. (i.ne.i2)) then
              if (validSplit(verts(i1,:),verts(i2,:),verts(i,:))) then
                 t(isplit) = split(verts(i1,:),verts(i2,:),verts(i,:))
                 ti(isplit) = i
                 isplit = isplit + 1                 
              end if
           end if
        end do
        call sort(nsplit, t, ti)
        lo(1) = edges0(iedge0,1)
        lo(2:) = ti
        hi(1:nsplit) = ti
        hi(nsplit+1) = edges0(iedge0,2)
        do i=1,nsplit+1
           edges(iedge,1) = lo(i)
           edges(iedge,2) = hi(i)
           edgeCon(iedge) = edgeCon0(iedge0)
           iedge = iedge + 1
        end do
        deallocate(t)
        deallocate(ti)
        deallocate(lo)
        deallocate(hi)
     end if
  end do

end subroutine splitEdges




subroutine sort(n, t, ti)

  implicit none

  !Input
  integer, intent(in) ::  n

  !Output
  double precision, intent(inout) ::  t(n)
  integer, intent(inout) ::  ti(n)

  !Working
  integer i, ti0(n)

  ti0 = ti
  do i=1,n
     ti(i) = ti0(minloc(t,1))
     t(minloc(t,1)) = 2.0
  end do

end subroutine sort




subroutine countSplits(nvert, nedge, verts, edges, nsplit)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nvert, nedge, verts, edges
  !f2py intent(out) nsplit
  !f2py depend(nvert) verts
  !f2py depend(nedge) edges

  !Input
  integer, intent(in) ::  nvert, nedge
  double precision, intent(in) ::  verts(nvert,2)
  integer, intent(in) ::  edges(nedge,2)

  !Output
  integer, intent(out) ::  nsplit

  !Working
  integer iedge, i, i1, i2
  logical validSplit

  nsplit = 0
  do iedge=1,nedge
     i1 = edges(iedge,1)
     i2 = edges(iedge,2)
     do i=1,nvert
        if ((i.ne.i1) .and. (i.ne.i2)) then
           if (validSplit(verts(i1,:),verts(i2,:),verts(i,:))) then
              nsplit = nsplit + 1
           end if
        end if
     end do
  end do

end subroutine countSplits




function validSplit(A, B, C)

  implicit none

  !Input
  double precision, intent(in) ::  A(2), B(2), C(2)

  !Output
  logical validSplit

  !Working
  double precision v1(2), v2(2), det

  v1 = C - A
  v2 = C - B
  det = abs(v1(1)*v2(2) - v1(2)*v2(1))

  if ((det .lt. 1e-10) .and. (dot_product(v1,v2) .lt. 0)) then
     validSplit = .True.
  else
     validSplit = .False.
  end if

end function validSplit




function split(A, B, C)

  implicit none

  !Input
  double precision, intent(in) ::  A(2), B(2), C(2)

  !Output
  double precision split

  !Working
  double precision d, d2

  d = sqrt(dot_product(C - A, C - A))
  d2 = sqrt(dot_product(B - A, B - A))

  split = d/d2

end function split
