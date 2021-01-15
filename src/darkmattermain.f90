program darkmattermain

    use parameters
    use kinds
    use espectra

    implicit none

    type(nucleus) :: nuc_target
    type(eftheory) :: eft
    type(particle) :: wimp

    real(kind=8) :: v
    integer :: computeoption
    integer(kind=8) :: clock_rate, tstart, tstop

    call system_clock(count_rate = clock_rate)

    print*,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    print*,'Welcome to our FORTAN 90 implementation of the WIMP-nucleon scattering code.'
    print*,'Based on the Mathematica script described in arXiv:1308.6288 (2003).'
    print*,'  VERSION 1.1 UPDATED: 2020.05.23 @ SDSU'
    print*,'  Dev. contact: cjohnson@sdsu.edu'
    print*,'                 ogorton@sdsu.edu'
    print*,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    print*,''
    print*,'Select an option:'
    print*,'[1] Event rate per detector per target nucleus per unit recoil energy spectra'
    print*,'[2] Scattering probability'
    print*,'[3] Differential cross section per recoil energy'
    print*,'[4] (Future feature) Total cross section'
    print*,'[5] (Future feature) Total scattering rate per detector'
    read*,computeoption

    if (computeoption==2.or.computeoption==3) then
        print*,"Enter darkmatter velocity:"
        read*,v
    endif

    print*,' '
    print*,'Enter the target proton number '
    read(5,*) nuc_target%Z
    print*,' '
    print*,'Enter the target neutron number '
    read(5,*) nuc_target%N

    call setparameters(nuc_target)
    call setupcoef(eft)

    call opencontrolfile(2)
    call readcontrolfile(2, eft, wimp)
    call convertisospinform(eft)

    call normalizecoeffs(eft, wimp)

    call setup_nuclearinputs(nuc_target)

    call printparameters(wimp,nuc_target,eft)

    if (computeoption == 1) then
        call system_clock(tstart)
        call eventrate_spectra(wimp, nuc_target, eft)
        call system_clock(tstop)
        print*,'Compute time (s)',real(tstop-tstart)/real(clock_rate)

    else if (computeoption == 2) then
        stop 'Not implemented'
        !output = transition_probability(q,v,jchi,y)
        !print*,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
        !print*,'Scattering probability = ',output
        !print*,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    else if (computeoption == 3) then
        stop 'Not implemented'    
        !output = diffCrossSection(v, q, jchi, y, Mtiso)
        !print*,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
        !print*,'Differential cross section = ',output
        !print*,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    endif

end program  
