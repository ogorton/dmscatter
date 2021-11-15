subroutine setup_nuclearinputs(nuc_target)

    use parameters
    use orbitals
    use constants
    use densities
    implicit none

    type(nucleus) :: nuc_target
    integer :: resfile = 100
    integer :: Atarget

    Atarget = nuc_target%A

    call openresults(resfile)
    call getorbitals(resfile,Atarget)
    call setupdensities(nuc_target)
    call readheaderv2(nuc_target,resfile)
    call readalldensities(nuc_target,resfile)

    if (fillcore) then
        call coredensity(nuc_target)
    end if
 
    if (printdens) call printdensities(nuc_target)

end subroutine
