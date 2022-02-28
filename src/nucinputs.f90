subroutine nucinputs

    use main
    use settings, only: fillcore, printdens
    use orbitals
    use densities
    implicit none

    integer :: resfile = 100
    integer :: Atarget

    Atarget = nuc_target%A

    call openresults(resfile)
    call getorbitals(resfile,Atarget)
    call setupdensities
    call readheaderv2(resfile)
    call readalldensities(resfile)

    if (fillcore) then
        call coredensity
    end if
 
    if (printdens) call printdensities

end subroutine
