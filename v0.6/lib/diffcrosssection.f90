function diffCrossSection(v, q, jchi, y, Mtiso)
    ! Computes the differential cross section per recoil energy ds/dEr
    use kinds
    use masses
    use constants
    implicit none
    interface
         function transition_probability(q,v,jchi,y,Mtiso)
            use kinds
            implicit none
            REAL(doublep), INTENT(IN) :: q
            REAL(doublep), INTENT(IN) :: v
            REAL(doublep), INTENT(IN) :: jchi
            REAL(doublep), INTENT(IN) :: y
            integer, intent(in) :: Mtiso
            real(doublep) :: transition_probability
        end function
    end interface
    real(doublep) :: diffCrossSection
    real(doublep) :: v ! velocity of DM particle in lab frame
    real(doublep) :: q
    real(doublep) :: y
    real(doublep) :: jchi
    integer :: Mtiso

    !q = sqrt(2.0 * Miso * (mN/GeV) * ERkeV) * 1E-3 * GeV
    !y = ((q * bfm)/2.0)**2.0

!    diffCrossSection = (2.0*Miso/(4.0*pi*v*v)) * transition_probability(q,v,jchi,y,Mtiso)
    diffCrossSection = (Miso*mN/(2.0*pi*v*v)) * transition_probability(q,v,jchi,y,Mtiso)

end function diffCrossSection
