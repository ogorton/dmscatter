function diffCrossSection(v, q, wimp, nucl, eft)
    ! Computes the differential cross section per recoil energy ds/dEr
    use kinds
    use constants
    use parameters
    implicit none
    interface
         function transition_probability(q,v,wimp,nucl,eft)
            use kinds
            use parameters
            implicit none
            REAL(doublep), INTENT(IN) :: q
            REAL(doublep), INTENT(IN) :: v
            type(particle) :: wimp
            type(nucleus) :: nucl
            real(kind=8), allocatable, intent(in) :: eft(:,:)
            real(doublep) :: transition_probability
        end function
    end interface
    type(particle) :: wimp
    type(nucleus) :: nucl
    real(kind=8), allocatable, intent(in) :: eft(:,:)

    real(doublep) :: diffCrossSection
    real(doublep) :: v ! velocity of DM particle in lab frame
    real(doublep) :: q

    diffCrossSection = (nucl%mass*mN/(2d0*pi*v*v)) * transition_probability(q,v,wimp,nucl,eft)

end function diffCrossSection
