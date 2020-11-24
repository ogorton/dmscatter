function EventRate(q, wimp, nuc_target, eft, detector_t)
    ! Computes the differential cross section per recoil energy ds/dEr
    use quadrature
    use kinds
    use constants
    use parameters
    
    implicit none
    interface
        function diffCrossSection(v, q, wimp, nucl, eft)
            use kinds
            use constants
            use parameters
            implicit none
            type(particle) :: wimp
            type(nucleus) :: nucl
            type(eftheory) :: eft
            real(doublep) :: diffCrossSection
            real(doublep) :: v ! velocity of DM particle in lab frame
            real(doublep) :: q
        end function diffcrossSection    
        function maxwell_boltzmann(v,v0)
            use kinds
            real(doublep), intent(in) :: v
            real(doublep), intent(in) :: v0
            real(doublep) :: maxwell_boltzmann            
        end function
    end interface
    real(doublep) :: EventRate
    real(doublep) :: q
    type(particle) :: wimp
    type(nucleus) :: nuc_target
    type(eftheory) :: eft
    type(detector) :: detector_t

    real(doublep) :: mchi, muT
    real(doublep) :: Nt
    real(doublep) :: rhochi

    real(doublep) :: v  ! DM velocity variable
    real(doublep) :: dv ! DM differential velocity / lattive spacing
    real(doublep), allocatable :: EventRate_integrand(:)
    integer :: i
    real(doublep) :: ve, v0

    ve = vdist_t%vearth
    v0 = vdist_t%vscale

    muT = wimp%mass * nuc_target%mass * mN / (wimp%mass + nuc_target%mass * mN)
    vdist_min = q/(2d0*muT)

!    print*,'Integrating dv from',vdist_min,'to',vdist_max
    dv = (vdist_max-vdist_min)/lattice_points

    allocate(EventRate_integrand(lattice_points))
   
!$OMP parallel do private(v) shared(wimp, nuc_target, eft)
    do i = 1, lattice_points
 
        v = vdist_min + (i-1) * dv
        EventRate_integrand(i) = diffCrossSection(v, q, wimp, nuc_target, eft) &
                    * v * v * ( maxwell_boltzmann(v-ve,v0) &
                            - maxwell_boltzmann(v+ve,v0) ) 
    end do
!$OMP end parallel do
!$OMP barrier

    if (quadrature_type == 1) then    
        call boole(lattice_points,EventRate_integrand,dv,EventRate)
    else
        EventRate = -1
    end if

    Nt = detector_t%Nt
    rhochi = wimp%localdensity
    Mchi = wimp%mass

    EventRate = Nt * (rhochi/Mchi) * EventRate *  (pi*v0**2/(ve))
  
end function EventRate
