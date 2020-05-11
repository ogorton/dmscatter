function EventRate(Nt, rhochi, ve, v0, q, jchi, y)
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
    real(doublep) :: tmp, diffcrosssection, vmin

    dv = 12*v0 / 1000
    Nv = (12*v0 )/dv
    vmin = q/(2.0*muT)

    print*,'Integral lattice size = ',Nv

    allocate(v_lattice(Nv))
    allocate(EventRate_integrand(Nv))
    
    print*,'vmin = ',vmin,'vmax = ',Nv*dv+vmin
    print*,'dv=',dv
    print*,'nt',nt
    print*,'rhochi',rhochi
    print*,'mchi',mchi

    do i = 1, Nv

        v = vmin + (i-1) * dv + 0.00001
        v_lattice(i) = v
        
        EventRate_integrand(i) = diffCrossSection(v, q, jchi, y, Mtiso)&
                    * v * maxwell_boltzmann(v,v0) *v*v
    end do
    print*,'here',v
    
    call boole(Nv,EventRate_integrand,dv,EventRate)

    EventRate = Nt * rhochi * EventRate * 4*pi/mchi
  
    print*,'Event Rate = ',EventRate

end function EventRate
