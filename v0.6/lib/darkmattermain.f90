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
    integer :: ap, an

    ! TEST BLOCK <<<<<a
    real(doublep) :: Nt ! Number of target nuclei
    real(doublep) :: rhochi ! local dark matter density
    real(doublep) :: ve ! Earth's velocity in the galactic rest frame
    real(doublep) :: v0 ! rms velocity of the visible matter distribution<
    real(kind=8) :: output
    real(kind=8) :: q,v,jchi,y,yy
    REAL(kind=8) :: nucResponse
    real(kind=8) :: femtometer, GeV, diffcrosssection
    character :: yn
    integer :: i

    GeV = 1.0
    femtometer = 5.0677/GeV
    q=1.*GeV
    v=1.
    jchi=0.5
    mchi=50.

    allocate(cvec(1)%c(17))
    allocate(cvec(0)%c(17))
    cvec(0)%c = 0.0
    cvec(0)%c = 0.0

    call opencoeffmatrix(2)
    call readcoeffmatrix(2)
    call normalizecoeffs

    ! >>>>>>> END TEST BLOCK
    
    call GetSPS

    print*,' '
    print*,' Reading one-body density matrix from Bigstick .res file '
    print*,' '

    nsporb = norb(1)

    call openresults(resfile)
    call setupdensities
    call readheaderv2(resfile)
    call readalldensities(resfile)
    print*,'Fill core? [y/n]'
    read*,yn
    if (yn=='y') then
        print*,'Filling core shell density matrix elements.'
        call coredensity
    end if
    call printdensities

    print*,' '
    print*,' Enter the neutron number '
    read(5,*)an
    print*,' '
    print*,' Enter the proton number '
    read(5,*)ap

    bfm = (41.467/(45.*(an+ap)**(-1./3) - 25.*(ap+an)**(-2./3)))**0.5 * femtometer
    y = (q*bfm/2.0)**2.0

    print*,'b[dimless]=',bfm/femtometer
    print*,'b[fm]=',bfm
    print*,'y=',y
    print*,'mN',mN
    print*,'jchi',jchi
    print*,'mchi',mchi

    Mtiso = (ap-an)
    Miso = ap+an
    muT = mchi * Miso * mN / (mchi+Miso*mN)

    print*,'Jiso, Tiso=',Jiso,Tiso

    print*,'Mtiso=',Mtiso
    print*,'Miso=',Miso
    print*,'muT=',muT
    print*,'v=',v
    print*,'q=',q

    print*,'denom',(4.0*mN*mchi)**2.0
    print*,'cvec',cvec(0)%c(:)
    output = nucResponse(0,0,1,y)
    print*,'Nuclear response =',output

    output = transition_probability(q,v,jchi,y)
    print*,"Transition probability =",output

    output = diffCrossSection(v, q, jchi, y, Mtiso)
    print*,'Differential cross section =',output

    Nt = 1.
    rhochi = 1.
    ve = 232.
    v0 = 220.
    output= eventrate(Nt, rhochi, ve, v0, q, jchi, y)

end program  
