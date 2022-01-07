module norm
contains

    function Jnorm(j)
        implicit none
        integer, INTENT(IN)  :: j
        real(kind=8) :: Jnorm

        Jnorm = 2d0 * j + 1d0

    end function Jnorm


    function j2norm(two_j)
        implicit none
        integer, intent(in) :: two_j
        integer :: j2norm
    
        j2norm = two_j + 1
    
    end function j2norm
    

    function Qnorm(j)
        implicit none
        integer, INTENT(IN) :: j
        real(kind=8)  :: Qnorm

        if (j<0) STOP "qnorm error: negative j"

        Qnorm = sqrt(2d0 * j + 1d0)

    end function Qnorm


    function q2norm(two_j)
        implicit none
        integer, intent(in) :: two_j
        real(kind=8) :: q2norm
        if (two_j<0) STOP "qnorm error: negative j"

        q2norm = sqrt(two_j + 1d0)
    
    end function q2norm    

end module norm
