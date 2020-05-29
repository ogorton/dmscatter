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
    use quadrature
    implicit none

    integer :: resfile = 100
    real(kind=8) :: output, eventrate
    character :: yn

    print*,'Enter q, the three-momentum transfer of the scattering reaction:'
    read*,q
    print*,' '
    print*,' Enter the neutron number '
    read(5,*)num_n
    print*,' '
    print*,' Enter the proton number '
    read(5,*)num_p

    call setparameters

    call setupcoef

    call opencontrolfile(2)
    call readcontrolfile(2)

    call normalizecoeffs
    call convertisospinform

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

    print*,'b[dimless]=',bfm/femtometer
    print*,'b[fm]=',bfm
    print*,'y=',y
    print*,'mN',mN
    print*,'jchi',jchi
    print*,'mchi',mchi
    print*,'Jiso=',Jiso
    print*,'Tiso=',Tiso
    print*,'Mtiso=',Mtiso
    print*,'ap,an',num_p,num_n
    print*,'Miso=',Miso
    print*,'muT=',muT
    print*,'q=',q
    print*,'vdist_min = ',vdist_min
    print*,'vdist_max = ',vdist_max
    print*,'nt',nt
    print*,'rhochi',rhochi
    print*,'mchi',mchi
    print*,'v0=',v0
    print*,'ve=',ve
    print*,'Integral lattice size = ',lattice_points

    output= eventrate(Nt, rhochi, ve, v0, q, jchi, y)
    print*,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    print*,'Event rate = ',output
    print*,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'

end program  
