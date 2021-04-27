module eventrate
    implicit none
    real (kind=8) :: globalq
    contains
function deventrate(q, wimp, nuc_target, eft)
    use quadrature
    use kinds
    use constants
    use parameters
   
    use crosssection
    use Mmaxbolt
    implicit none
    real(doublep) :: dEventRate
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
    integer, allocatable :: indx(:)
    integer :: i, j, itmp, ind
    real(doublep) :: ve, v0, vesc, vmin, error

    real(doublep) :: abserror, relerror

    logical :: adaptive = .true.

    globalq = q
    muT = wimp%mass * nuc_target%mass * mN / (wimp%mass + nuc_target%mass * mN)

    ve = vearth * kilometerpersecond
    v0 = vscale * kilometerpersecond
    vesc = vescape * kilometerpersecond
    vmin = q/(2d0*muT)

    if (vmin > vesc) then
        dEventRate = 0d0
        return
    end if
    dv = (vesc-vmin)/lattice_points 

    allocate(EventRate_integrand(lattice_points))
    allocate(xtab(lattice_points))

    abserror=dspectra(vesc) ! This line prevents a segfault... idk why
    relerror = 1e-12
    abserror = 1e-12


    if (adaptive) then
        call gaus8 ( dspectra, vmin, vesc, abserror, deventrate, ind )
    else
        !$OMP parallel do private(v) shared(wimp, nuc_target, eftsmall) &
        !$OMP schedule(dynamic,10)
            do i = 1, lattice_points
                v = (vmin + (i-1) * dv)
                xtab(i) = v
                EventRate_integrand(i) = diffCrossSection(v, q, wimp, nuc_target, eft) &
                            * v * v * ( maxbolt(v-ve,v0) - maxbolt(v+ve,v0) ) 
            end do
        !$OMP end parallel do
        !$OMP barrier
!            do i = 1, lattice_points
!                write(1099,*)xtab(i)/kilometerpersecond,EventRate_integrand(i)
!            end do
        
            if (quadrature_type == 1) then    
                if (all(EventRate_integrand == 0)) then
                    dEventrate = 0.0
                else
                    call cubint ( lattice_points, xtab, EventRate_integrand, 1, &
                        lattice_points, dEventRate, error )
                end if
            else
                dEventRate = -1
            end if
    end if
    Nt = nuc_target%Nt
    rhochi = wimp%localdensity / centimeter**3d0
    Mchi = wimp%mass
    dEventRate = ntscale * kilogramday * Nt * (rhochi/Mchi) * dEventRate *  (pi*v0**2/ve)
  
end function deventrate

function dspectra(vv)
    ! The purpose of this function is to create a callable func for the
    ! adaptive integral routine. Other functions in this program use
    ! (the good practice of) explicit data dependency. An exception
    ! was made for this fuction for compatibility with the adaptive
    ! integration routine used to integrate the event rate spectra.

    use main ! This is the only function allowed to use main.
    use kinds
    use constants
    use crosssection
    use mmaxbolt
    implicit none
    real(doublep) :: vv
    real(doublep) :: dspectra, ve, v0
    ve = vearth * kilometerpersecond
    v0 = vscale * kilometerpersecond

    dspectra = diffCrossSection(vv, globalq, wimp, nuc_target, eft) &
        * vv * vv * ( maxbolt(vv-ve,v0) - maxbolt(vv+ve,v0) )

end function

end module eventrate
