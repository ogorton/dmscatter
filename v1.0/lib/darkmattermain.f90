!  ctrlchar = 'c'  count up # of states
!  ctrlchar = 'f'  fill up info on states
!  ctrlchar = 'p'  fill up parent reference states
!  ctrlchar = 'd'  fill up daughter reference states

!================================================

program darkmattermain
    use response
    use spspace
    use stateinfo
    use masses
    use kinds
    use constants
    use dmparticles
    use velocities
    use momenta
    implicit none

    interface
        subroutine openresults(resfile)
            use stateinfo
            implicit none
            integer resfile
            character(22):: filename
        end subroutine openresults
        function transition_probability(q,v,jchi,y)
            implicit none
            REAL(kind=8), INTENT(IN) :: q
            REAL(kind=8), INTENT(IN) :: v
            REAL(kind=8), INTENT(IN) :: jchi
            REAL(kind=8), INTENT(IN) :: y
            real(kind=8) :: transition_probability
        end function
        function eventrate(Nt, rhochi, ve, v0, q, jchi, y)
            use kinds
            implicit none
            real(doublep) :: EventRate
            real(doublep), intent(in) :: Nt ! Number of target nuclei
            real(doublep), intent(in) :: rhochi ! local dark matter density
            real(doublep), intent(in) :: ve ! Earth's velocity in the galactic rest frame
            real(doublep), intent(in) :: v0 ! rms velocity of the visible matter distribution 
            real(doublep), intent(in) :: q
            real(doublep), intent(in) :: y
            real(doublep), intent(in) :: jchi
        end function
    end interface

    integer, parameter :: resfile = 33

    real(kind=8) :: output
    character :: yn

    call setparameters
    call setupcoef
    call GetSPS

    nsporb = norb(1)

    call openresults(resfile)
    call setupdensities
    call readheaderv2(resfile)
    call readalldensities(resfile)

    print*,'Fill core? [y/n]'
    read*,yn
    if (yn=='y') then
        call coredensity
    end if
    call printdensities

    print*,' '
    print*,' Enter the neutron number '
    read(5,*)num_n
    print*,' '
    print*,' Enter the proton number '
    read(5,*)num_p

    print*,'Enter q, the three-momentum transfer of the scattering reaction:'
    read*,q

    call setparameters

    output= eventrate(Nt, rhochi, ve, v0, q, jchi, y)
    print*,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    print*,'Event rate = ',output
    print*,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'

end program  
