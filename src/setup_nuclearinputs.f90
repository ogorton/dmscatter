subroutine setup_nuclearinputs(nuc_target)

    use parameters
    use spspace
    use constants
    implicit none

    type(nucleus) :: nuc_target
    integer :: resfile = 100
    character :: yn

    call GetSPS
    call openresults(resfile)
    call setupdensities(nuc_target)
    call readheaderv2(nuc_target,resfile)
    call readalldensities(nuc_target,resfile)

    print*,'Fill core? [y/n]'
    read*,yn
    if (yn=='y') then
        call coredensity(nuc_target)
    end if

    bfm = (41.467/(45.*(nuc_target%mass)**(-1./3) &
                - 25.*(nuc_target%mass)**(-2./3)))**0.5 * femtometer

end subroutine
