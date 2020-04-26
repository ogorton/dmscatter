function EventRate(Nt, rhochi, ve, v0, q, jchi, y, Mtiso)
    ! Computes the differential cross section per recoil energy ds/dEr
    use kinds
    use masses
        ! Uses: 
        !    mchi : Dark matter mass
        !    Miso : nuclear isotope mass
    use constants ! pi
    
    implicit none
    interface
!         function diffCrossSection(ERkev, v, q, jchi, y, Mtiso)
!            use kinds
!            implicit none
!            real(doublep) :: diffCrossSection
!            real(doublep) :: ERkev
!            real(doublep) :: v
!            real(doublep) :: q
!            real(doublep) :: y
!            real(doublep) :: jchi
!            integer :: Mtiso
!        end function
        function maxwell_boltzmann(v,v0)
            use kinds
            real(doublep), intent(in) :: v
            real(doublep), intent(in) :: v0
            real(doublep) :: maxwell_boltzmann            
        end function
    end interface
    real(doublep) :: EventRate
    real(doublep) :: Nt ! Number of target nuclei
    real(doublep) :: rhochi ! local dark matter density
    real(doublep) :: ve ! Earth's velocity in the galactic rest frame
    real(doublep) :: v0 ! rms velocity of the visible matter distribution 
    real(doublep) :: q
    real(doublep) :: y
    real(doublep) :: jchi
    integer :: Mtiso

    real(doublep) :: ERkeV ! Recoil energy of the nucleus
   
    real(doublep) :: v  ! DM velocity variable
    real(doublep) :: dv ! DM differential velocity / lattive spacing
    real(doublep),allocatable :: v_lattice(:) ! velocity domain
    real(doublep), allocatable :: EventRate_integrand(:)
    integer :: i, Nv
    real(doublep) :: tmp, diffcrosssection

    dv = v0 / 100
    Nv = (12*v0 )/dv

    print*,'Integral lattice size = ',Nv

    allocate(v_lattice(Nv))
    allocate(EventRate_integrand(Nv))

    print*,'vmin = ',-Nv*dv/2,'vmax = ',Nv*dv/2

    ERkev = 1.0
    v = 1.
    q = 1.
    jchi = 0.5
    y = 1.
    mtiso = 0.
    print*,'aa0'
    tmp = diffCrossSection(ERkev,v,q,jchi,y,Mtiso)
    print*,'aa1'

    do i = 1, Nv

        print*,i
        v = (i-Nv/2-1) * dv
        print*,'g'
        v_lattice(i) = v
        print*,'k',ERkev,v,q,jchi,y,mtiso
        tmp = diffCrossSection(ERkev, v, q, jchi, y, Mtiso)
        print*,'h'
        
        EventRate_integrand(i) = diffCrossSection(ERkev, v, q, jchi, y, Mtiso)&
                    * v * maxwell_boltzmann(v,v0)
        print*,'g'
    end do
    print*,'here'
    
    ERkev = q*q/(2.0*Miso)

    call simpson(Nv,Nv,EventRate_integrand,dv,EventRate)

    EventRate = Nt * rhochi * EventRate
  
    print*,'Event Rate = ',EventRate

end function EventRate
