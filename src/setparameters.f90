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
    GeV = 1.0d0
    femtometer = 5.0677d0/GeV

    ! Masses
    mV = 246.2d0
    mN = 0.938272d0 

    !Mtiso = num_p-num_n
    nuc_target%Mt = nuc_target%N - nuc_target%Z
    !mtarget = num_p+num_n
    nuc_target%mass = nuc_target%Z + nuc_target%N
    nuc_target%nt = 1.0d0 / ( nuc_target%mass * mN)
    !muT = mchi * mtarget * mN / (mchi+mtarget*mN)

    ! Velocities
    vdist_t%vearth = 232.0d0 !km/s earths velocity in the galactic rest frame.
    vdist_t%vscale = 220.0d0 !km/s velocity distribution scaling
    vdist_t%vescape = 12.0d0 * vdist_t%vscale !km/s ! infinity

    bfm = (41.467d0/(45d0*(nuc_target%mass)**(-1.d0/3) &
            - 25.d0*(nuc_target%mass)**(-2.d0/3)))**0.5d0 * femtometer 

    ! Quadrature
    quadrature_type = 1
    lattice_points = 1000

end subroutine setparameters

subroutine printparameters(wimp,nuc_target,eft)

    use constants
    use quadrature
    use momenta
    use spspace
    use parameters

    implicit none
    integer i
    type(nucleus) :: nuc_target
    type(eftheory) :: eft
    type(particle) :: wimp

    print*,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    print*,'Parameters used in this calculation:'
    print*,'femtometer =',femtometer
    print*,'GeV = ',gev
    print*,'kev = ',kev
    print*,'b (fm)',bfm/femtometer
    print*,'b (1/GeV)',bfm
    print*,'Nucleon mass (GeV)',mN
    print*,'WIMP spin',wimp%j
    print*,'WIMP mass (GeV)',wimp%mass
    print*,'Target J',nuc_target%groundstate%Jx2
    print*,'Target T',nuc_target%groundstate%Tx2
    print*,'Target Mt',nuc_target%Mt
    print*,'Target Z,N',nuc_target%Z, nuc_target%N
    print*,'Target atomic mass',nuc_target%mass
    print*,'System reduced mass',wimp%mass * nuc_target%mass * mN/(wimp%mass+nuc_target%mass*mn)!mchi * mtarget * mN / (mchi+mtarget*mN)
    print*,'Target mass density (1/GeV)',nuc_target%nt
    print*,'Local WIMP density (GeV/cm^3)',wimp%localdensity
    print*,'v0 (km/s)',vdist_t%vscale
    print*,'ve (km/s)',vdist_t%vearth
    print*,'v escape (km/s)',vdist_t%vescape
    print*,'Integral lattice size = ',lattice_points
    print*,'EFT coupling coefficients:'
    print*,'i    p    n    s    v'
    write(6,"(A,T12,A,T24,A,T36,A,T48,A)")'i','p','n','s','v'
    do i = 1, num_response_coef
        if ((eft%xpnc(0)%c(i).ne.0) &
            .or. (eft%xpnc(1)%c(i).ne.0) &
            .or. (eft%isoc(0)%c(i).ne.0) &
            .or. (eft%isoc(1)%c(i).ne.0)) then
                write(6,"(I2,x2e11.4,e11.4,e11.4,e11.4)")i,&
                        eft%xpnc(0)%c(i),eft%xpnc(1)%c(i),eft%isoc(0)%c(i),eft%isoc(1)%c(i)
        end if
    end do
    print*,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'

end subroutine printparameters
