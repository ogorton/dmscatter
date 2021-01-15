module main

    use parameters
    implicit none

    type(nucleus) :: nuc_target
    type(eftheory) :: eft
    type(particle) :: wimp

    real(kind=8) :: v    

contains

subroutine computefunction(computeoption)
    use parameters
    use kinds
    use espectra
    implicit none
    integer, intent(in) :: computeoption
    integer(kind=8) :: clock_rate, tstart, tstop    

    call system_clock(count_rate = clock_rate)
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

    select case(computeoption)
      case(5)
        call system_clock(tstart)
        call totaleventrate
        call system_clock(tstop)

      case(1)
        call system_clock(tstart)
        call eventrate_spectra(wimp, nuc_target, eft)
        call system_clock(tstop)
      case(2)
        stop 'Not implemented'
        !output = transition_probability(q,v,jchi,y)
        !print*,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
        !print*,'Scattering probability = ',output
        !print*,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
      case(3)
        stop 'Not implemented'
        !output = diffCrossSection(v, q, jchi, y, Mtiso)
        !print*,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
        !print*,'Differential cross section = ',output
        !print*,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
      case default
        print*,"Invalid compute option."
    end select    

    print*,'Compute time (s)',real(tstop-tstart)/real(clock_rate)

end subroutine computefunction
end module main
