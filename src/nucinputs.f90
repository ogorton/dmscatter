subroutine setup_nuclearinputs(nuc_target)

    use parameters
    use spspace
    use constants
    implicit none

    type(nucleus) :: nuc_target
    integer :: resfile = 100

    call GetSPS
    call openresults(resfile)
    call setupdensities(nuc_target)
    call readheaderv2(nuc_target,resfile)
    call readalldensities(nuc_target,resfile)

    if (fillcore) then
        call coredensity(nuc_target)
    end if
 
    !call printdensities(nuc_target)

end subroutine
