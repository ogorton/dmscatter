module settings
    
    logical :: usemomentum = .false.
    logical :: useenergyfile = .false.
    logical :: fillcore = .true.
    logical :: printdens = .false.
    logical :: pnresponse = .false.
    integer :: nlocalstates !,maxlocalstates
    integer :: nsporb   ! # of single-particle orbits

end module settings
