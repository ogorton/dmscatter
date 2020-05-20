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
         function transition_probability(q,v,jchi,y,Mtiso)
            use kinds
            implicit none
            REAL(doublep), INTENT(IN) :: q
            REAL(doublep), INTENT(IN) :: v
            REAL(doublep), INTENT(IN) :: jchi
            REAL(doublep), INTENT(IN) :: y
            integer, intent(in) :: Mtiso
            real(doublep) :: transition_probability
        end function
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
    real(doublep) :: tmp, diffcrosssection, vmin, vmax

    Nv = 1000
    vmin = q/(2.0*muT)

    vmax = 12*v0
    dv = (vmax-vmin)/Nv

    print*,'Integral lattice size = ',Nv

    allocate(EventRate_integrand(Nv))
    
    print*,'vmin = ',vmin,'vmax = ',Nv*dv+vmin
    print*,'dv=',dv
    print*,'nt',nt
    print*,'rhochi',rhochi
    print*,'mchi',mchi
    print*,'v0=',v0
    print*,'ve=',ve

    do i = 1, Nv
 
        v = vmin + (i-1) * dv
 
        write(22,*)v,maxwell_boltzmann(v+ve,v0)
        
        EventRate_integrand(i) = diffCrossSection(v, q, jchi, y, Mtiso)&
                    * v * v * ( + maxwell_boltzmann(v-ve,v0) &
                            - maxwell_boltzmann(v+ve,v0) ) &
                    * (pi*v0**2/(ve))
        write(23,*),v,EventRate_integrand(i)/Mchi
    end do
    
    call simpson(Nv,Nv,EventRate_integrand,dv,EventRate)

    EventRate = Nt * (rhochi/Mchi) * EventRate
  
end function EventRate
