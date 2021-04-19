subroutine computefunction(computeoption)
    use main
    use parameters
    use kinds
    use espectra
    use sj2iref
    use spspace, only: bfm
    use constants, only: kev, mN, kilometerpersecond

    use Mtransition_probability
    use Mdiffcrosssection
    implicit none
    integer, intent(in) :: computeoption
    integer(kind=8) :: clock_rate, tstart, tstop    
    real(kind=8) :: e, q, v, output

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
      case default
        print*,"Invalid compute option."
    end select    
    call system_clock(tstop)
    print*,'Compute time (s)',real(tstop-tstart)/real(clock_rate)
    print*,"6-J Table Lookup min/max requested:",tableJmin, tableJmax

end subroutine computefunction
