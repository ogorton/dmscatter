module distributions
    contains

    function maxbolt(v,v0)
        ! Maxwell-Boltzmann (MB) probability distribution 
        use kinds
        use constants
        implicit none
    
        real(dp), intent(in) :: v 
        real(dp), intent(in) :: v0
        real(dp) :: maxbolt
    
        maxbolt = exp(-(v/v0)**2.0d0) / ((pi)**1.5d0 * v0**3.0d0)
    
    end function maxbolt

    function shm(v, v0, vesc)
        ! Standard halo model probability distribution, i.e.
        ! a MB distribution with a sharp cutoff to recognize the
        ! finite escape speed of the galaxy.
        use kinds
        use constants
        implicit none

        real(dp), intent(in) :: v, v0, vesc
        real(dp) :: shm
        real(dp) :: Nresc, z, y

        if (v > vesc) then
            shm = 0d0
            return
        end if

        z = vesc / v0
        y = v / v0

        Nresc = erf(z) - 2d0*z * exp(-z*z)/sqrt(pi)
        shm = exp(-y*y) / (Nresc * (pi)**1.5d0 * v0**3d0)

    end function shm

    function sshm(v, v0, vesc)
        ! Smooth standard halo model. Basically a SHM with an extra constant 
        ! term meanth to smooth the distribution near vesc.
        use kinds
        use constants
        implicit none

        real(dp), intent(in) :: v, v0, vesc
        real(dp) :: sshm
        real(dp) :: Nresc, z, y

        if (v > vesc) then
            sshm = 0d0
            return
        end if        

        z = vesc / v0
        y = v / v0

        Nresc = erf(z) - 2*z*(1d0 + z*z/1.5d0)*exp(-z*z)/sqrt(pi)
        sshm = (exp(-y*y) - exp(-z*z)) / (Nresc * (pi)**1.5d0 * v0**3d0)

    end function sshm

end module

