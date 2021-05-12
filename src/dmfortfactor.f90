program dmfortfactor

    implicit none
    integer :: computeoption

    print*,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    print*,'Welcome to our FORTAN 90 implementation of the WIMP-nucleon scattering code.'
    print*,'Based on the Mathematica script described in arXiv:1308.6288 (2003).'
    print*,'  VERSION 1.4 UPDATED: 2021.01.15 @ SDSU'
    print*,'  Dev. contact: cjohnson@sdsu.edu'
    print*,'                 ogorton@sdsu.edu'
    print*,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    print*,''
    print*,'Select an option:'
    print*,'[1] Event rate per unit recoil energy (spectra)'
    print*,'[2] Scattering probability'
    print*,'[3] Differential cross section per recoil energy'
    print*,'[4] (*) Total cross section'
    print*,'[5] Total scattering events per detector (does not produce spectra data)'
    print*,'[6] Nuclear response function test'
    print*," "
    print*,"(*) Not available in this release."
    read*,computeoption
    if (computeoption>6) stop "Invalid option."

    call setupandrun(computeoption)

    contains

        subroutine setupandrun(computeoption)
            use main
            use mod_spectra, only: eventrate_spectra, velocity_curve
            use wignerfunctions, only: tableJmin, tableJmax, sj2itable
            use nucresponse, only: test_nucresponse
            implicit none
            integer, intent(in) :: computeoption
            integer(kind=8) :: clock_rate, tstart, tstop
        
            call system_clock(count_rate = clock_rate)
            print*,' '
            print*,'Enter the target proton number '
            read(5,*) nuc_target%Z
            print*,' '
            print*,'Enter the target neutron number '
            read(5,*) nuc_target%N
        
            nuc_target%A = nuc_target%Z + nuc_target%N
        
            call setparameters(nuc_target)
            call setupcoef(eft)
            call opencontrolfile(2)
            call readcontrolfile(2, eft, wimp)
            call convertisospinform(eft)
            call normalizecoeffs(eft, wimp)
            call setup_nuclearinputs(nuc_target)
            call printparameters(wimp,nuc_target,eft)
            call sj2itable
        
            call system_clock(tstart)
            select case(computeoption)
              case(5)
                call totaleventrate(nuc_target)
              case(1)
                call eventrate_spectra(wimp, nuc_target, eft)
              case(2)
                  call velocity_curve(wimp, nuc_target, eft,2)
              case(3)
                  call velocity_curve(wimp, nuc_target, eft,3)
              case(6)
                  call test_nucresponse(nuc_target)                 
              case default
                print*,"Invalid compute option."
            end select
            call system_clock(tstop)
            print*,'Compute time (s)',real(tstop-tstart)/real(clock_rate)
            print*,"6-J Table Lookup min/max requested:",tableJmin, tableJmax
        
        end subroutine setupandrun

end program dmfortfactor
