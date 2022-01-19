module eventrate

    implicit none
    real(kind=8) :: qglobal(128)
    real(kind=8) :: ve, v0, vesc
    contains

    function deventrate(q, wimp, nuc_target, eft)
        use OMP_LIB
        use quadrature
        use kinds
        use constants
        use types

        implicit none
        real(dp) :: dEventRate
        real(dp) :: q
        type(particle) :: wimp
        type(nucleus) :: nuc_target
        type(eftheory) :: eft
    
        real(dp) :: mchi, muT
        real(dp) :: Nt
        real(dp) :: rhochi
    
        integer :: ind, tid
        real(dp) :: vmin, error, vmax
        real(dp) :: relerror
    
        tid = 1
        ! Don't delete the following line; it's an openMP command, not a comment.
        !$  tid = omp_get_thread_num() + 1
        qglobal(tid) = q 
   
        Mchi = wimp%mass 
        muT = Mchi * nuc_target%mass * mN / (Mchi + nuc_target%mass * mN)
        ve = vearth * kilometerpersecond
        v0 = vscale * kilometerpersecond
        vesc = vescape * kilometerpersecond
        Nt = nuc_target%Nt
        rhochi = wimp%localdensity / centimeter**3d0
        vmin = q/(2d0*muT)
        relerror = quadrature_relerr
        vmax = (vesc + ve)
    
        select case(quadrature_type)
          case(1)
            call gaus8_threadsafe(spectraintegrand1d, vmin, vmax, &
                relerror, deventrate, ind, tid )
          case(2)
            deventrate = gaussquad(spectraintegrand1d, gaussorder, vmin, vmax, tid)
          case default
            stop "Invalid quadrature option."
        end select
    
        dEventRate = ntscale * kilogramday * Nt * (rhochi/Mchi) &
            * dEventRate *  (pi*v0**2/ve)
      
    end function deventrate
    
    function spectraintegrand1d(vv, tid)
        ! The purpose of this function is to create a callable func for the
        ! adaptive integral routine. Other functions in this program use
        ! (the good practice of) explicit data dependency. An exception
        ! was made for this fuction for compatibility with the adaptive
        ! integration routine used to integrate the event rate spectra.
    
        use main ! This is the only function allowed to use main.
        use kinds
        use crosssection
        use distributions
    
        implicit none
        real(dp) :: vv, qq
        real(dp) :: spectraintegrand1d
        integer tid
    
        qq = qglobal(tid)
    
        spectraintegrand1d = diffCrossSection(vv, qq, wimp, nuc_target, eft) &
            * (sshm(vv-ve,v0,vesc) - sshm(vv+ve,v0,vesc)) * vv**2
    
    end function

end module eventrate
