program dmfortfactor

    use main
    use mod_spectra, only: eventrate_spectra, velocity_curve
    use wigner
    use nucresponse, only: test_nucresponse

    implicit none
    integer(kind=8) :: clock_rate, tstart, tstop    
    character(len=2) :: option
    logical :: usingopenMP
    usingopenMP = .false.
    !$ usingopenMP = .true.

    print *,'================================================================================='
    print *,'Welcome to our modern Fortran implementation of the WIMP-nucleon scattering code.'
    print *,'Based on the Mathematica script described in arXiv:1308.6288 (2003).'
    print *,'  VERSION 1.9 UPDATED: 2021.10.21 @ SDSU'
    print *,'  Developer contacts:'
    print *,'            cjohnson@sdsu.edu'
    print *,'             ogorton@sdsu.edu'
    print *,' openMP:',usingopenMP
    print *,'================================================================================='
    print *,''
    print *,'Select an option:'
    print *,'[er] Event rate per unit recoil energy (spectra)'
    print *,'[cs] Differential cross section per recoil energy'
    print *,'[tp] Scattering probability'
    print *,'[te] Total scattering events per detector (does not produce spectra data)'
    print *,'[wd] Nuclear response function'
    print *,'[ws] Nuclear response spectra'
    print *,''

    read *,option

    call system_clock(count_rate = clock_rate)
    print *,' '
    print *,'Enter the target proton number '
    read(5,*) nuc_target%Z
    print *,' '
    print *,'Enter the target neutron number '
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
    call threej_table_init()
    call sixj_table_init()
    
    call system_clock(tstart)
    select case(option)
      case('te')
        call totaleventrate(nuc_target)
      case('er')
        call eventrate_spectra(wimp, nuc_target, eft)
      case('cs')
        call velocity_curve(wimp, nuc_target, eft,2)
      case('tp')
        call velocity_curve(wimp, nuc_target, eft,3)
      case('wd')
        call test_nucresponse(nuc_target)  
      case('ws')
        call nucresponse_spectra(nuc_target)        
      case default
        print*,"Invalid compute option."
    end select
    call system_clock(tstop)

    print*,'Compute time (s)',real(tstop-tstart)/real(clock_rate)
    print*,"6-J Table Lookup min/max requested:", tablemin_used, tablemax_used
        
end program dmfortfactor
