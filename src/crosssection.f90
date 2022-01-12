module crosssection
contains
    function diffCrossSection(v, q, wimp, nucl, eft)
        ! Computes the differential cross section per recoil energy ds/dEr
        use kinds
        use constants
        use types
        use transition
        implicit none
        type(particle) :: wimp
        type(nucleus) :: nucl
        type(eftheory) :: eft
    
        real(dp) :: diffCrossSection
        real(dp) :: v ! velocity of DM particle in lab frame
        real(dp) :: q
        
        diffCrossSection = (nucl%mass*mN/(2d0*pi*v*v)) * transition_probability(q,v,wimp,nucl,eft)
    
    end function diffCrossSection
end module
