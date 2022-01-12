module distributions
    contains
    function maxbolt(v,v0)
        use kinds
        use constants
        implicit none
    
        real(dp), intent(in) :: v 
        real(dp), intent(in) :: v0
        real(dp) :: maxbolt
    
        maxbolt = exp(-(v/v0)**2.0d0) / ((pi)**1.5d0 * v0**3.0d0)
    
    end function maxbolt

end module

