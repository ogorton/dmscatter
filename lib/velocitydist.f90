function maxwell_boltzmann(v,v0)
    use kinds
    use constants
    implicit none

    real(doublep), intent(in) :: v 
    real(doublep), intent(in) :: v0
    real(doublep) :: maxwell_boltzmann

    maxwell_boltzmann = exp(-(v/v0)**2.0) / ((pi)**1.5d0 * v0**3.0)

end function maxwell_boltzmann



