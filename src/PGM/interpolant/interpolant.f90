subroutine computeCoons(u, v, induv, ind0v, ind1v, indu0, indu1, &
     ind00, ind10, ind01, ind11, Da, Di, Dj)

  implicit none

  !Input
  double precision, intent(in) ::  u, v
  integer, intent(in) ::  induv
  integer, intent(in) ::  ind0v, ind1v, indu0, indu1
  integer, intent(in) ::  ind00, ind10, ind01, ind11

  !Output
  double precision, intent(out) ::  Da(8)
  integer, intent(out) ::  Di(8), Dj(8)

  Di(1:8) = induv
  Da(1) = -(1-u) * (1-v)
  Dj(1) = ind00
  Da(2) = -u * (1-v)
  Dj(2) = ind10
  Da(3) = -(1-u) * v
  Dj(3) = ind01
  Da(4) = -u * v
  Dj(4) = ind11
  Da(5) = 1-u
  Dj(5) = ind0v
  Da(6) = u
  Dj(6) = ind1v
  Da(7) = 1-v
  Dj(7) = indu0
  Da(8) = v
  Dj(8) = indu1

end subroutine computeCoons



subroutine sparseBezier(u, s0, s1, C)

  implicit none

  !Input
  double precision, intent(in) ::  u, s0, s1

  !Output
  double precision, intent(out) ::  C(4)

  !Working
  double precision u1, u2, u3, t1, t2, t3

  u1 = u
  u2 = u*u1
  u3 = u*u2
  t1 = (1-u)
  t2 = (1-u)*t1
  t3 = (1-u)*t2

  ! A -> B         C -> D
  ! *--------------*
  C(1) = t3 + (3-s0) * u1*t2
  C(2) = s0 * u1*t2
  C(3) = u3 + (s1+3) * u2*t1
  C(4) = -s1 * u2*t1

end subroutine sparseBezier
