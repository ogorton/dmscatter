module distributions
    use constants

    contains

    function maxbolt(v,v0)
        ! Maxwell-Boltzmann (MB) probability distribution 
        
        implicit none
    
        real(kind=8), intent(in) :: v 
        real(kind=8), intent(in) :: v0
        real(kind=8) :: maxbolt
    
        maxbolt = exp(-(v/v0)**2.0d0) / ((pi)**1.5d0 * v0**3.0d0)
    
    end function maxbolt

    function shm(v, v0, vesc)
        ! Standard halo model probability distribution, i.e.
        ! a MB distribution with a sharp cutoff to recognize the
        ! finite escape speed of the galaxy.
        
        implicit none

        real(kind=8), intent(in) :: v, v0, vesc
        real(kind=8) :: shm
        real(kind=8) :: Nresc, z, y

        if (v > vesc) then
            shm = 0d0
            return
        end if

        z = vesc / v0
        y = v / v0

        Nresc = erf(z) - 2d0*z * exp(-z*z)/sqrt(pi)
        shm = exp(-y*y) / (Nresc * (pi)**1.5d0 * v0**3d0)

    end function shm

    function sshm(v, ve, v0, vesc)
        ! Smooth standard halo model. Basically a SHM with an extra constant 
        ! term meanth to smooth the distribution near vesc.
        
        implicit none

        real(kind=8), intent(in) :: v, ve, v0, vesc
        real(kind=8) :: sshm
        real(kind=8) :: Nesc, z

        z = vesc / v0

        Nesc = (erf(z)-2d0*z*(1d0+z*z/1.5d0)*exp(-z*z)/sqrt(pi))*pi**1.5d0*v0**3d0
        sshm = (Imbcutoff(v, ve, v0, vesc) - Isccutoff(v, ve, v0, vesc))/Nesc

    end function sshm

    function Imbcutoff(v, ve, v0, vesc)
        ! Integrand for radial part of Maxwell Boltzmann distribution with cutoff
        
        implicit none

        real(kind=8), intent(in) :: v, ve, v0, vesc
        real(kind=8) :: Imbcutoff

        if (v<vesc-ve) then
            Imbcutoff = g(v-ve,v0)-g(v+ve,v0)
        else ! v>vesc-ve
            Imbcutoff = g(v-ve,v0)-g(vesc,v0)
        end if
        Imbcutoff = Imbcutoff * pi*(v0**2d0/ve)

    end function Imbcutoff

    function Isccutoff(v, ve, v0, vesc)
        ! Integrand for radial part of smooth component with cutoff
        
        implicit none

        real(kind=8), intent(in) :: v, ve, v0, vesc
        real(kind=8) :: Isccutoff

        if (v<vesc-ve) then
            Isccutoff = 2d0 * v
        else ! v>vesc-ve
            Isccutoff = (vesc**2d0 - (v-ve)**2d0)/(2d0*ve)
        end if
        Isccutoff = Isccutoff * 2d0*pi*g(vesc,v0)

    end function Isccutoff

    function g(v, v0)
        ! Gaussian term
        
        implicit none

        real(kind=8), intent(in) :: v
        real(kind=8), intent(in) :: v0
        real(kind=8) :: g

        g = exp(-(v/v0)**2d0) 

    end function

end module

