module norm
    contains
!-------------------------------------------------------------------------------
function Jnorm(j)
    implicit none
    REAL(kind=8), INTENT(IN)  :: j
    REAL(kind=8) :: Jnorm
    Jnorm = 2.0 * j + 1.0
end function Jnorm


!-------------------------------------------------------------------------------
function Qnorm(j)
    implicit none
    REAL(kind=8), INTENT(IN) :: j
    Real(kind=8)  :: Qnorm
    if (j<0) STOP "qnorm error: negative j" ! O.C.G: avoid NAN
    Qnorm = sqrt(2.0 * j + 1.0)
end function Qnorm
end module norm
