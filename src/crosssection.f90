module crosssection
    use main
    use constants
    use transition
    contains
    function diffCrossSection(v, q)
        ! Computes the differential cross section per recoil energy ds/dEr
        implicit none
        real(kind=8) :: diffCrossSection
        real(kind=8) :: v ! velocity of DM particle in lab frame
        real(kind=8) :: q
        
        diffCrossSection = (nuc_target%mass*mN/(2d0*pi*v*v)) * transition_probability(v, q)
    
    end function diffCrossSection
end module
