function EventRate(Nt, rhochi, ve, v0, q, jchi, y)
    ! Computes the differential cross section per recoil energy ds/dEr
    use quadrature
    use kinds
    use masses
        ! Uses: 
        !    mchi : Dark matter mass
        !    Miso : nuclear isotope mass
    use constants ! pi
    
    implicit none
    interface
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

    real(doublep) :: v  ! DM velocity variable
    real(doublep) :: dv ! DM differential velocity / lattive spacing
    real(doublep),allocatable :: v_lattice(:) ! velocity domain
    real(doublep), allocatable :: EventRate_integrand(:)
    integer :: i
    real(doublep) :: tmp, diffcrosssection


    dv = (vdist_max-vdist_min)/lattice_points

    allocate(EventRate_integrand(lattice_points))
    
    do i = 1, lattice_points
 
        v = vdist_min + (i-1) * dv
 
        write(22,*)v,maxwell_boltzmann(v+ve,v0)
        
        EventRate_integrand(i) = diffCrossSection(v, q, jchi, y, Mtiso)&
                    * v * v * ( + maxwell_boltzmann(v-ve,v0) &
                            - maxwell_boltzmann(v+ve,v0) ) 

        write(23,*),v,EventRate_integrand(i)/Mchi

    end do

    if (quadrature_type == 1) then    
        call simpson(lattice_points,lattice_points,EventRate_integrand,dv,EventRate)
    else
        EventRate = -1
   end if

    EventRate = Nt * (rhochi/Mchi) * EventRate *  (pi*v0**2/(ve))
  
end function EventRate
