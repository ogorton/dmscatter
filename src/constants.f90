module constants
    
    implicit none
    real(kind=8) :: pi = datan(1d0)*4d0!3.14159265358979323846264338327950288419716939937510
    real(kind=8) :: GeV = 1d0
    real(kind=8) :: kev = 10d0**(-6)
    real(kind=8) :: femtometer = 5.0677d0 ! /GeV
    real(kind=8) :: centimeter = 5.0677d0*10d0**13d0
    real(kind=8) :: kilometerpersecond = 3.33564d0 * 10d0**(-6d0)
    real(kind=8) :: kilogramday = 7.3634d0*10d0**55d0
    real(kind=8) :: mn = 0.938272d0
    real(kind=8) :: mv = 246.2d0

    real(kind=8) :: vearth  = 232.0d0!km/s
    REAL(kind=8) :: vescape = 12 * 232.d0! km/s
    REAL(kind=8) :: vscale  = 220.0d0! km/s 

    real(kind=8) :: ntscale = 1d0

end module constants
