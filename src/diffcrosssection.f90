module mdiffcrosssection
    contains
function diffCrossSection(v, q, wimp, nucl, eft)
    ! Computes the differential cross section per recoil energy ds/dEr
    use kinds
    use constants
    use parameters
    use Mtransition_probability
    implicit none
    type(particle) :: wimp
    type(nucleus) :: nucl
    type(eftheory) :: eft

    real(doublep) :: diffCrossSection
    real(doublep) :: v ! velocity of DM particle in lab frame
    real(doublep) :: q
    
    diffCrossSection = (nucl%mass*mN/(2d0*pi*v*v)) * transition_probability(q,v,wimp,nucl,eft)

end function diffCrossSection
end module
