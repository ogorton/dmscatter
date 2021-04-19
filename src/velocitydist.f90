module Mmaxbolt
    contains
function maxbolt(v,v0)
    use kinds
    use constants
    implicit none

    real(doublep), intent(in) :: v 
    real(doublep), intent(in) :: v0
    real(doublep) :: maxbolt

    maxbolt = exp(-(v/v0)**2.0d0) / ((pi)**1.5d0 * v0**3.0d0)

end function maxbolt

end module

