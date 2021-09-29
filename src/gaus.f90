module gaus
contains

function g8(func, x, h, tid)

    implicit none
    real ( kind = 8 ), external :: func
    real ( kind = 8 ) h
    real ( kind = 8 ) x
    integer tid

    real (kind = 8 ) :: g8

    real ( kind = 8 ) :: w1 = 3.62683783378361983D-01
    real ( kind = 8 ) :: w2 = 3.13706645877887287D-01
    real ( kind = 8 ) :: w3 = 2.22381034453374471D-01
    real ( kind = 8 ) :: w4 = 1.01228536290376259D-01
    real ( kind = 8 ) :: x1 = 1.83434642495649805D-01
    real ( kind = 8 ) :: x2 = 5.25532409916328986D-01
    real ( kind = 8 ) :: x3 = 7.96666477413626740D-01
    real ( kind = 8 ) :: x4 = 9.60289856497536232D-01

    g8 = h * ( ( & 
                w1 * ( func ( x - x1 * h, tid ) + func ( x + x1 * h, tid ) )   &   
              + w2 * ( func ( x - x2 * h, tid ) + func ( x + x2 * h, tid ) ) ) & 
            + ( w3 * ( func ( x - x3 * h, tid ) + func ( x + x3 * h, tid ) )   &   
              + w4 * ( func ( x - x4 * h, tid ) + func ( x + x4 * h, tid ) ) ) )

end function g8          

subroutine gaus8_threadsafe ( func, a, b, err, result, ierr, tid )

!*****************************************************************************80
!
!! GAUS8 estimates the integral of a function.
!   Modified by O. Gorton c. May 2021 for use with openMP.
!
!  Discussion:
!
!    GAUS8 integrates real functions of one variable over finite
!    intervals using an adaptive 8-point Legendre-Gauss
!    algorithm.
!
!    GAUS8 is intended primarily for high accuracy integration or
!    integration of smooth functions.
!
!  Modified:
!
!    30 October 2000
!
!  Author:
!
!    Ron Jones,
!    Sandia National Laboratory,
!    Los Alamos, New Mexico
!
!  Reference:
!
!    Philip Davis, Philip Rabinowitz,
!    Methods of Numerical Integration,
!    Second Edition,
!    Dover, 2007,
!    ISBN: 0486453391,
!    LC: QA299.3.D28.
!
!  Parameters:
!
!    Input, real ( kind = 8 ), external FUNC, name of external function to
!    be integrated.  This name must be in an external statement in the
!    calling program.  FUNC must be a function of one real argument.  The value
!    of the argument to FUNC is the variable of integration
!    which ranges from A to B.
!
!    Input, real ( kind = 8 ) A, the lower limit of integration.
!
!    Input, real ( kind = 8 ) B, the upper limit of integration.
!
!    Input/output, real ( kind = 8 ) ERR.
!    On input, ERR is a requested pseudorelative error
!    tolerance.  Normally pick a value of ABS ( ERR ) so that
!    STOL < ABS ( ERR ) <= 1.0D-3 where STOL is the single precision
!    unit roundoff.
!    RESULT will normally have no more error than
!    ABS ( ERR ) times the integral of the absolute value of
!    FUN(X).  Usually, smaller values for ERR yield more
!    accuracy and require more function evaluations.
!    A negative value for ERR causes an estimate of the
!    absolute error in RESULT to be returned in ERR.  Note that
!    ERR must be a variable (not a constant) in this case.
!    Note also that the user must reset the value of ERR
!    before making any more calls that use the variable ERR.
!    On output, ERR will be an estimate of the absolute error
!    in RESULT if the input value of ERR was negative.  ERR is
!    unchanged if the input value of ERR was non-negative.
!    The estimated error is solely for information to the user
!    and should not be used as a correction to the computed integral.
!
!    Output, real ( kind = 8 ) RESULT, the computed value of the integral.
!
!    Output, integer ( kind = 4 ) IERR, a status code.
!    Normal Codes:
!     1 RESULT most likely meets requested error tolerance, or A = B.
!    -1 A and B are too nearly equal to allow normal
!        integration.  RESULT is set to zero.
!     Abnormal Code:
!     2 RESULT probably does not meet requested error tolerance.
!
!    input, real arg : extra argument for func
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) aa(30)
  real ( kind = 8 ) ae
  real ( kind = 8 ) anib
  real ( kind = 8 ) area
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) ce
  real ( kind = 8 ) ee
  real ( kind = 8 ) ef
  real ( kind = 8 ) eps
  real ( kind = 8 ) err
  real ( kind = 8 ) est
  real ( kind = 8 ), external :: func
  real ( kind = 8 ) gl
  real ( kind = 8 ) glr
  real ( kind = 8 ) gr(30)
  real ( kind = 8 ) hh(30)
  integer ( kind = 4 ), save :: icall = 0
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) k
  integer ( kind = 4 ), save :: kml = 6
  integer ( kind = 4 ), save :: kmx = 5000
  integer ( kind = 4 ) l
  integer ( kind = 4 ) lmn
  integer ( kind = 4 ) lmx
  integer ( kind = 4 ) lr(30)
  integer ( kind = 4 ) mxl
  integer ( kind = 4 ) nbits
  integer ( kind = 4 ) nib
  integer ( kind = 4 ), save :: nlmn = 1
  integer ( kind = 4 ) nlmx
  real ( kind = 8 ) result
  real ( kind = 8 ) tol
  real ( kind = 8 ) vl(30)
  real ( kind = 8 ) vr

  integer tid

  if ( a == b ) then
    err = 0.0D+00
    result = 0.0D+00
    return
  end if
 
!  if ( icall /= 0 ) then
!    write ( *, '(a)' ) ' '
!    write ( *, '(a)' ) 'GAUS8 - Fatal error!'
!    write ( *, '(a)' ) '  GAUS8 was called recursively.'
!    stop 1
!  end if 

  icall = 1
!
!  DIGITS ( X ) = number of base 2 digits in representation of X.
!
  k = digits ( result )

  anib = log10 ( 2.0D+00 ) * real ( k, kind = 8 ) / 0.30102000D+00
  nbits = int ( anib )
  nlmx = min ( 30, ( nbits * 5 ) / 8 )
  result = 0.0D+00
  ierr = 1
  ce = 0.0D+00
  result = 0.0D+00
  lmx = nlmx
  lmn = nlmn
 
  if ( b /= 0.0D+00 ) then

    if ( sign ( 1.0D+00, b ) * a <= 0.0D+00 ) then
      go to 10
    end if

    c = abs ( 1.0D+00 - a / b )
    if ( 0.1D+00 < c ) then
      go to 10
    end if

    if ( c <= 0.0D+00 ) then
      icall = 0
      if ( err < 0.0D+00 ) then
        err = ce
      end if
      return
    end if

    anib = 0.5D+00 - log ( c ) / 0.69314718D+00
    nib = int ( anib )
    lmx = min ( nlmx, nbits-nib-7 )

    if ( lmx < 1 ) then
      ierr = -1
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'GAUS8 - Warning!'
      write ( *, '(a)' ) '  A and B are too close to carry out integration.'
      icall = 0
      if ( err < 0.0D+00 ) then
        err = ce
      end if
      return
    end if

    lmn = min ( lmn, lmx )

  end if
 
10    continue
 
  tol = max ( abs ( err ), 2.0D+00**(5-nbits)) / 2.0D+00
  if ( err == 0.0D+00 ) then
    tol = sqrt ( epsilon ( 1.0D+00 ) )
  end if

  eps = tol
  hh(1) = ( b - a ) / 4.0D+00
  aa(1) = a
  lr(1) = 1
  l = 1
  est = g8 (func, aa(l) + 2.0D+00 * hh(l), 2.0D+00 * hh(l), tid )
  k = 8
  area = abs ( est )
  ef = 0.5D+00
  mxl = 0
!
!  Compute refined estimates, estimate the error, etc.
!
20 continue
 
  gl = g8 (func, aa(l) + hh(l), hh(l), tid )
  gr(l) = g8 (func, aa(l) + 3.0D+00 * hh(l), hh(l), tid )
  k = k + 16
  area = area + ( abs ( gl ) + abs ( gr(l) ) - abs ( est ) )
 
  glr = gl + gr(l)
  ee = abs ( est - glr ) * ef
  ae = max ( eps * area, tol * abs ( glr ) )

  if ( ee - ae <= 0.0D+00 ) then
    go to 40
  else
    go to 50
  end if
 
30 continue
 
  mxl = 1
 
40 continue
 
  ce = ce + ( est - glr )
 
  if ( lr(l) <= 0 ) then
    go to 60
  else
    go to 80
  end if
!
!  Consider the left half of this level
!
50 continue

  if ( kmx < k ) then
    lmx = kml
  end if

  if ( lmx <= l ) then
    go to 30
  end if

  l = l + 1
  eps = eps * 0.5D+00
  ef = ef / sqrt ( 2.0D+00 )
  hh(l) = hh(l-1) * 0.5D+00
  lr(l) = -1
  aa(l) = aa(l-1)
  est = gl
  go to 20
!
!  Proceed to right half at this level
!
60 continue

  vl(l) = glr

70 continue

  est = gr(l-1)
  lr(l) = 1
  aa(l) = aa(l) + 4.0D+00 * hh(l)
  go to 20
!
!  Return one level
!
80 continue

  vr = glr

90 continue

  if ( l <= 1 ) then
    go to 120
  end if

  l = l - 1
  eps = eps * 2.0D+00
  ef = ef * sqrt ( 2.0D+00 )
 
  if ( lr(l) <= 0 ) then
    vl(l) = vl(l+1) + vr
    go to 70
  else
    vr = vl(l+1) + vr
    go to 90
  end if
!
!  Exit
!
120   continue
 
  result = vr

  if ( mxl == 0 .or. abs ( ce ) <= 2.0D+00 * tol * area ) then 
    icall = 0 
    if ( err < 0.0D+00 ) then
      err = ce
    end if
    return
  end if

  ierr = 2
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'GAUS8 - Warning!'
  write ( *, '(a)' ) '  RESULT is probably insufficiently accurate.'
  icall = 0

  if ( err < 0.0D+00 ) then
    err = ce
  end if

  return
end subroutine gaus8_threadsafe

end module gaus
