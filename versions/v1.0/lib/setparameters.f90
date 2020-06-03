subroutine setparameters

    use constants
    use masses
    use dmparticles
    use velocities
    use stateinfo
    use quadrature 
    use momenta

    implicit none

    ! Constants
    GeV = 1.0
    femtometer = 5.0677/GeV


    ! Masses
    mV = 246.2
    mN = 0.938272 
    Mtiso = num_p-num_n
    Miso = num_p+num_n
    muT = mchi * Miso * mN / (mchi+Miso*mN)

    ! Velocities
    ve = 232.0 !km/s earths velocity in the galactic rest frame.
    v0 = 220.0 !km/s velocity distribution scaling
    vesc = 12.0 * v0 !km/s ! infinity

    bfm = (41.467/(45.*(num_n+num_p)**(-1./3) - 25.*(num_n+num_p)**(-2./3)))**0.5 * femtometer 
    y = (q*bfm/2.0)**2.0

    ! dm particles
    Nt = 1.0
    jchi=0.5
    mchi=50.
    rhochi = 1.0

    ! Quadrature
    quadrature_type = 1
    lattice_points = 1000
    vdist_min = q/(2.0*muT)
    vdist_max = 12*v0


    print*,'b[dimless]=',bfm/femtometer
    print*,'b[fm]=',bfm
    print*,'y=',y
    print*,'mN',mN
    print*,'jchi',jchi
    print*,'mchi',mchi
    print*,'Jiso, Tiso=',Jiso,Tiso
    print*,'Mtiso=',Mtiso
    print*,'Miso=',Miso
    print*,'muT=',muT
    print*,'q=',q
    print*,'vdist_min = ',vdist_min,'vdist_max = ',vdist_max
    print*,'nt',nt
    print*,'rhochi',rhochi
    print*,'mchi',mchi
    print*,'v0=',v0
    print*,'ve=',ve
    print*,'Integral lattice size = ',lattice_points

end subroutine setparameters
