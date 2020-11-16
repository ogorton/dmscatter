!-------------------------------------------------------------------------------
subroutine chgm ( a, b, x, hg )

!*****************************************************************************80
!
!! CHGM computes the confluent hypergeometric function M(a,b,x).
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    27 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, parameters.
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) HG, the value of M(a,b,x).
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) a0
  real ( kind = 8 ) a1
  !real ( kind = 8 ) aa
  real ( kind = 8 ) b
  real ( kind = 8 ) hg
  real ( kind = 8 ) hg1
  real ( kind = 8 ) hg2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) la
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nl
  real ( kind = 8 ) pi
  real ( kind = 8 ) r
  real ( kind = 8 ) r1
  real ( kind = 8 ) r2
  real ( kind = 8 ) rg
  real ( kind = 8 ) sum1
  real ( kind = 8 ) sum2
  real ( kind = 8 ) ta
  real ( kind = 8 ) tb
  real ( kind = 8 ) tba
  real ( kind = 8 ) x
  real ( kind = 8 ) x0
  real ( kind = 8 ) xg
  real ( kind = 8 ) y0
  real ( kind = 8 ) y1

  pi = 3.141592653589793D+00
  a0 = a
  a1 = a
  x0 = x
  hg = 0.0D+00

  if ( b == 0.0D+00 .or. b == - abs ( int ( b ) ) ) then
    hg = 1.0D+300
  else if ( a == 0.0D+00 .or. x == 0.0D+00 ) then
    hg = 1.0D+00
  else if ( a == -1.0D+00 ) then
    hg = 1.0D+00 - x / b
  else if ( a == b ) then
    hg = exp ( x )
  else if ( a - b == 1.0D+00 ) then
    hg = ( 1.0D+00 + x / b ) * exp ( x )
  else if ( a == 1.0D+00 .and. b == 2.0D+00 ) then
    hg = ( exp ( x ) - 1.0D+00 ) / x
  else if ( a == int ( a ) .and. a < 0.0D+00 ) then
    m = int ( - a )
    r = 1.0D+00
    hg = 1.0D+00
    do k = 1, m
      r = r * ( a + k - 1.0D+00 ) / k / ( b + k - 1.0D+00 ) * x
      hg = hg + r
    end do
  end if

  if ( hg /= 0.0D+00 ) then
    return
  end if

  if ( x < 0.0D+00 ) then
    a = b - a
    a0 = a
    x = abs ( x )
  end if

  if ( a < 2.0D+00 ) then
    nl = 0
  end if

  if ( 2.0D+00 <= a ) then
    nl = 1
    la = int ( a )
    a = a - la - 1.0D+00
  end if

  do n = 0, nl

    if ( 2.0D+00 <= a0 ) then
      a = a + 1.0D+00
    end if

    if ( x <= 30.0D+00 + abs ( b ) .or. a < 0.0D+00 ) then

      hg = 1.0D+00
      rg = 1.0D+00
      do j = 1, 500
        rg = rg * ( a + j - 1.0D+00 ) &
          / ( j * ( b + j - 1.0D+00 ) ) * x
        hg = hg + rg
        if ( abs ( rg / hg ) < 1.0D-15 ) then
          exit
        end if
      end do

    else

      call gamma ( a, ta )
      call gamma ( b, tb )
      xg = b - a
      call gamma ( xg, tba )
      sum1 = 1.0D+00
      sum2 = 1.0D+00
      r1 = 1.0D+00
      r2 = 1.0D+00
      do i = 1, 8
        r1 = - r1 * ( a + i - 1.0D+00 ) * ( a - b + i ) / ( x * i )
        r2 = - r2 * ( b - a + i - 1.0D+00 ) * ( a - i ) / ( x * i )
        sum1 = sum1 + r1
        sum2 = sum2 + r2
      end do
      hg1 = tb / tba * x ** ( - a ) * cos ( pi * a ) * sum1
      hg2 = tb / ta * exp ( x ) * x ** ( a - b ) * sum2
      hg = hg1 + hg2

    end if

    if ( n == 0 ) then
      y0 = hg
    else if ( n == 1 ) then
      y1 = hg
    end if

  end do

  if ( 2.0D+00 <= a0 ) then
    do i = 1, la - 1
      hg = ( ( 2.0D+00 * a - b + x ) * y1 + ( b - a ) * y0 ) / a
      y0 = y1
      y1 = hg
      a = a + 1.0D+00
    end do
  end if

  if ( x0 < 0.0D+00 ) then
    hg = hg * exp ( x0 )
  end if

  a = a1
  x = x0

  return
end




!-------------------------------------------------------------------------------
subroutine gamma ( x, ga )

!*****************************************************************************80
!
!! GAMMA evaluates the Gamma function.
!
!  Licensing:
!
!    The original FORTRAN77 version of this routine is copyrighted by 
!    Shanjie Zhang and Jianming Jin.  However, they give permission to 
!    incorporate this routine into a user program that the copyright 
!    is acknowledged.
!
!  Modified:
!
!    08 September 2007
!
!  Author:
!
!    Original FORTRAN77 version by Shanjie Zhang, Jianming Jin.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!    X must not be 0, or any negative integer.
!
!    Output, real ( kind = 8 ) GA, the value of the Gamma function.
!
  implicit none

  real ( kind = 8 ), dimension ( 26 ) :: g = (/ &
    1.0D+00, &
    0.5772156649015329D+00, &
   -0.6558780715202538D+00, &
   -0.420026350340952D-01, &
    0.1665386113822915D+00, &
   -0.421977345555443D-01, &
   -0.96219715278770D-02, &
    0.72189432466630D-02, &
   -0.11651675918591D-02, &
   -0.2152416741149D-03, &
    0.1280502823882D-03, &
   -0.201348547807D-04, &
   -0.12504934821D-05, &
    0.11330272320D-05, &
   -0.2056338417D-06, &
    0.61160950D-08, &
    0.50020075D-08, &
   -0.11812746D-08, &
    0.1043427D-09, &
    0.77823D-11, &
   -0.36968D-11, &
    0.51D-12, &
   -0.206D-13, &
   -0.54D-14, &
    0.14D-14, &
    0.1D-15 /)
  real ( kind = 8 ) ga
  real ( kind = 8 ) gr
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m1
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r
  real ( kind = 8 ) x
  real ( kind = 8 ) z

  if ( x == aint ( x ) ) then

    if ( 0.0D+00 < x ) then
      ga = 1.0D+00
      m1 = int ( x ) - 1
      do k = 2, m1
        ga = ga * k
      end do
    else
      ga = 1.0D+300
    end if

  else

    if ( 1.0D+00 < abs ( x ) ) then
      z = abs ( x )
      m = int ( z )
      r = 1.0D+00
      do k = 1, m
        r = r * ( z - real ( k, kind = 8 ) )
      end do
      z = z - real ( m, kind = 8 )
    else
      z = x
    end if

    gr = g(26)
    do k = 25, 1, -1
      gr = gr * z + g(k)
    end do

    ga = 1.0D+00 / ( gr * z )

    if ( 1.0D+00 < abs ( x ) ) then
      ga = ga * r
      if ( x < 0.0D+00 ) then
        ga = - pi / ( x* ga * sin ( pi * x ) )
      end if
    end if

  end if

  return
end
