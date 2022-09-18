program dmscatter

    use spectra, only: eventrate_spectra, velocity_curve, nucresponse_spectra
    use wigner
    use nucresponse, only: test_nucresponse
    use keywords, only: controlfile
    use totaleventrate
    use parameters

    implicit none
    integer(kind=8) :: clock_rate, tstart, tstop    
    character(len=2) :: option
    logical :: usingopenMP

    call system_clock(count_rate = clock_rate)
    usingopenMP = .false.
    !$ usingopenMP = .true.

    print '(10(A,/))',&
        '=================================================================================',&
        'Welcome to dmscatter, a modern Fortran WIMP-nucleon scattering code.',&
        'Based on the Mathematica script described in arXiv:1308.6288 (2003).',&
        '  VERSION 1.0 RELEASED 2022 from SDSU',&
        '  MIT License',&
        '  Copyright (c) 2022, dmscatter development team',&
        '  Developer contacts:',&
        '            cjohnson@sdsu.edu',&
        '             ogorton@sdsu.edu',&
        '================================================================================='
    print '(8(A,/))',&
        'Select an option:',&
        '[er] Event rate per unit recoil energy (spectra)',&
        '[cs] Differential cross section per recoil energy',&
        '[tp] Scattering probability',&
        '[te] Total scattering events per detector (does not produce spectra data)',&
        '[wd] Nuclear response function',&
        '[ws] Nuclear response spectra',&
        ''

    read *,option

    call setparameters
    call controlfile
    call nucinputs

    call printparameters

    call threej_table_init
    call sixj_table_init
    
    call system_clock(tstart)
    select case(option)
      case('te')
        call totaler
      case('er')
        call eventrate_spectra
      case('cs')
        call velocity_curve(2)
      case('tp')
        call velocity_curve(3)
      case('wd')
        call test_nucresponse
      case('ws')
        call nucresponse_spectra     
      case default
        print*,"Invalid compute option."
    end select
    call system_clock(tstop)

    print*,'Compute time (s)',real(tstop-tstart)/real(clock_rate)
    print*,"6-J Table Lookup min/max requested:", tablemin_used, tablemax_used
        
end program dmscatter
