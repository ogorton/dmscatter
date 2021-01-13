function EventRate(q, wimp, nuc_target, eft)
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
            real(doublep), allocatable :: eft(:,:)
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

    real(doublep) :: mchi, muT
    real(doublep) :: Nt, units
    real(doublep) :: rhochi
    real(doublep), allocatable :: eftsmall(:,:)

    real(doublep) :: v  ! DM velocity variable
    real(doublep) :: dv ! DM differential velocity / lattive spacing
    real(doublep), allocatable :: EventRate_integrand(:), xtab(:)
    integer :: i
    real(doublep) :: ve, v0, vesc, intscale, error

    allocate(eftsmall(0:1,num_response_coef))
    eftsmall(0,:) = eft%isoc(0)%c
    eftsmall(1,:) = eft%isoc(1)%c

    muT = wimp%mass * nuc_target%mass * mN / (wimp%mass + nuc_target%mass * mN)

    ve = vdist_t%vearth * kilometerpersecond
    v0 = vdist_t%vscale * kilometerpersecond
    vesc = vdist_t%vescape * kilometerpersecond
    vdist_min = q/(2d0*muT)
    dv = (vesc-vdist_min)/lattice_points 

    allocate(EventRate_integrand(lattice_points))
    allocate(xtab(lattice_points))
   
!$OMP parallel do private(v) shared(wimp, nuc_target, eft)
    do i = 1, lattice_points
        v = (vdist_min + (i-1) * dv)
        xtab(i) = v
        EventRate_integrand(i) = diffCrossSection(v, q, wimp, nuc_target, eftsmall) &
                    * v * v * ( maxwell_boltzmann(v-ve,v0) &
                            - maxwell_boltzmann(v+ve,v0) ) 
    end do
!$OMP end parallel do
!$OMP barrier

    if (quadrature_type == 1) then    
        call cubint ( lattice_points, xtab, EventRate_integrand, 1, &
                lattice_points, EventRate, error )
    else
        EventRate = -1
    end if

    Nt = nuc_target%Nt
    rhochi = wimp%localdensity / centimeter**3d0
    Mchi = wimp%mass
    EventRate = kilogramday * Nt * (rhochi/Mchi) * EventRate *  (pi*v0**2/ve)
  
end function EventRate
