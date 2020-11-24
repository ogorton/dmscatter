subroutine setparameters(nuc_target)

    use constants
    use quadrature 
    use momenta
    use spspace
    use parameters

    implicit none

    type(nucleus) :: nuc_target

    print*,'Setting default parameter values.'

    ! Constants
    GeV = 1.0
    femtometer = 5.0677/GeV

    ! Masses
    mV = 246.2
    mN = 0.938272 

    !Mtiso = num_p-num_n
    nuc_target%Mt = nuc_target%Z - nuc_target%N
    !mtarget = num_p+num_n
    nuc_target%mass = nuc_target%Z + nuc_target%N
    nuc_target%nt = 1.0
    !muT = mchi * mtarget * mN / (mchi+mtarget*mN)

    ! Velocities
    vdist_t%vearth = 232.0 !km/s earths velocity in the galactic rest frame.
    vdist_t%vscale = 220.0 !km/s velocity distribution scaling
    vdist_t%vescape = 12.0 * vdist_t%vscale !km/s ! infinity

    bfm = (41.467/(45d0*(nuc_target%mass)**(-1./3) - 25.*(nuc_target%mass)**(-2./3)))**0.5 * femtometer 

    ! Quadrature
    quadrature_type = 1
    lattice_points = 1000
    vdist_max = 12*vdist_t%vscale

end subroutine setparameters

subroutine printparameters(wimp,nuc_target,eft,detector_t)

    use constants
    use quadrature
    use momenta
    use spspace
    use parameters

    implicit none

    type(nucleus) :: nuc_target
    type(eftheory) :: eft
    type(particle) :: wimp
    type(detector) :: detector_t    

    print*,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    print*,'Parameters used in this calculation:'
    print*,'femtometer =',femtometer
    print*,'GeV = ',gev
    print*,'kev = ',kev
    print*,'b[dimless]=',bfm/femtometer
    print*,'b[fm]=',bfm
    print*,'mN',mN
    print*,'jchi',wimp%j
    print*,'mchi',wimp%mass
    print*,'J_target=',nuc_target%groundstate%Jx2
    print*,'T_target=',nuc_target%groundstate%Tx2
    print*,'Mt_targe=',nuc_target%Mt
    print*,'Target Z,N',nuc_target%Z, nuc_target%N
    print*,'mtarget=',nuc_target%mass
    print*,'muT=',wimp%mass * nuc_target%mass * mN/(wimp%mass+nuc_target%mass*mn)!mchi * mtarget * mN / (mchi+mtarget*mN)
    print*,'vdist_max = ',vdist_max
    print*,'Nt',nuc_target%nt
    print*,'rhochi',wimp%localdensity
    print*,'v0=',vdist_t%vscale
    print*,'ve=',vdist_t%vearth
    print*,'Integral lattice size = ',lattice_points
    print*,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'

end subroutine printparameters
