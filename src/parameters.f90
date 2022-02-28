module parameters

    use main
    use constants
    use quadrature
    use orbitals, only: bfm
    use settings
    use keywords

    contains
    subroutine setparameters
    
        implicit none
    
        print *,' '
        print *,'Enter the target proton number '
        read(5,*) nuc_target%Z
        print *,' '
        print *,'Enter the target neutron number '
        read(5,*) nuc_target%N
    
        nuc_target%A = nuc_target%Z + nuc_target%N
    
        print*,'Setting default parameter values.'
    
        nuc_target%Mt = + nuc_target%Z - nuc_target%N ! note non-standard definition
        nuc_target%mass = nuc_target%Z + nuc_target%N
        nuc_target%nt = 1.0d0 / ( nuc_target%mass * mN)
    
        bfm = (41.467d0/(45d0*(nuc_target%mass)**(-1.d0/3) &
                - 25.d0*(nuc_target%mass)**(-2.d0/3)))**0.5d0 * femtometer 
    
    end subroutine setparameters
    
    subroutine printparameters
    
        implicit none
        integer i
    
        print*,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
        print*,'Control file settings used:'
        do i = 1, numkeywords
            print '(A, "    ", F20.10)', keywordpairs(i)%key, keywordpairs(i)%val
        end do
    
        print*,''
        print*,'Parameters used in this calculation:'
    
        print'("kilometerpersecond",T30,ES20.5)',kilometerpersecond
        print'("kilogamday",T30,ES20.5)',kilogramday
        print'("femtometer =",T30,F10.6)',femtometer
        print'("GeV = ",T30,F10.6)',gev
        print'("kev = ",T30,F10.6)',kev
        print'("b (fm)",T30,F10.6)',bfm/femtometer
        print'("b (1/GeV)",T30,F10.6)',bfm
        print'("Nucleon mass (GeV)",T30,F10.6)',mN
        print'("WIMP spin",T30,F10.6)',wimp%j
        print'("WIMP mass (GeV)",T30,F10.6)',wimp%mass
        print'("Target J",T30,I2)',nuc_target%groundstate%Jx2
        print'("Target T",T30,I2)',nuc_target%groundstate%Tx2
        print'("Target Mt",T30,I2)',nuc_target%Mt
        print'("Target Z,N",T30,I2,"  ",I2)',nuc_target%Z, nuc_target%N
        print'("Target atomic mass",T30,F10.6)',nuc_target%mass
        print'("System reduced mass",T30,F10.6)',wimp%mass * nuc_target%mass * mN/(wimp%mass+nuc_target%mass*mn)!mchi * mtarget * mN / (mchi+mtarget*mN)
        print'("Target mass density (1/GeV)",T30,F10.6)',nuc_target%nt
        print'("Target experimental mass (kg*days)",T40,F20.6)',ntscale
        print'("Local WIMP density (GeV/cm^3)",T30,F10.6)',wimp%localdensity
        print'("v0 (km/s)",T30,F10.6)',vscale
        print'("ve (km/s)",T30,F10.6)',vearth
        print'("v escape (km/s)",T30,F15.6)',vescape
        print'("quadtype",T30,I2)',quadrature_type
        select case(quadrature_type)
        case(1)
        print'("quadrelerr",T30,F10.6)',quadrature_relerr
        case(2)
        print'("gaussorder",T30,I3)',gaussorder
        case default
         continue
        end select
        print'("useenergyfile",T30,L2)',useenergyfile   
    
        print*,''
        print*,'EFT coupling coefficients:'
        write(6,"(A,T12,A,T24,A,T36,A,T48,A)")'i','p','n','s','v'
        do i = 1, size(eft%xpnc(0)%c)
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

end module parameters
