module gausquad

contains

function gaussquad(func,n,a,b,tid) result(z)
  !
  ! This function computes the indefinite integral of the function func
  ! from a to b using an n-point Gauss-Legendre quadrature algorithm.
  !
  ! This function is thread-safe. 
  !
  ! The function func must be an external function which takes two arguments:
  !          func(x, tid)
  ! where x is a ready kind-8 and tid is an integer. Typically tid is used
  ! for accessing data unique to a thread, available in shared memory.
  !
  ! This function is modified from its original form found on Rosettacode.org:
  !     https://rosettacode.org
  !     https://rosettacode.org/wiki/Numerical_integration/Gauss-Legendre_Quadrature
  ! which is freely available under GNU free documentation license 1.2.
  !
  implicit none

  real (kind=8), external :: func
  integer, parameter :: p=8
  integer                 :: n
  real(kind=p), parameter :: pi = 4*atan(1._p)
  real(kind=p)            :: r(2, n), x, f, df, dx
  integer                 :: i,  iter, k
  real (kind = p) :: pc0(n+1), pc1(n+1), pctmp(n+1)
  real(kind=p)            :: z, a, b, zsum, vv
  integer tid

  pc0(1) = 1._p
  pc1(1:2) = [1._p, 0._p]
  do k = 2, n
    pctmp(1:k+1) = ((2*k-1)*[pc1(1:k),0._p]-(k-1)*[0._p, 0._p,pc0(1:k-1)])/k
    pc0 = pc1; pc1 = pctmp
  end do
  do i = 1, n
    x = cos(pi*(i-0.25_p)/(n+0.5_p))
    do iter = 1, 10
      f = pc1(1); df = 0._p
      do k = 2, size(pc1)
        df = f + x*df
        f  = pc1(k) + x * f
      end do
      dx =  f / df
      x = x - dx
      if (abs(dx)<10*epsilon(dx)) exit
    end do
    r(1,i) = x
    r(2,i) = 2/((1-x**2)*df**2)
  end do

  zsum = 0._p
  do i = 1, n
      vv = (a+b)/2 + r(1,i)*(b-a)/2
      zsum = zsum + r(2,i) * func(vv, tid)
  end do

  z = (b-a)/2 * zsum

end function

end module gausquad     
