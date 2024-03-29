module eventrate
    use OMP_LIB
    use quadrature
    use constants
    use crosssection
    use distributions
    use main

    implicit none
    real(kind=8) :: qglobal(128)
    real(kind=8) :: ve, v0, vesc
    contains

    function deventrate(q)

        implicit none
        real(kind=8) :: dEventRate
        real(kind=8) :: q
    
        real(kind=8) :: mchi, muT, mnuc
        real(kind=8) :: Nt
        real(kind=8) :: rhochi
    
        integer :: ind, tid
        real(kind=8) :: vmin, vmax
        real(kind=8) :: relerror
    
        tid = 1
        ! Don't delete the following line; it's an openMP command, not a comment.
        !$  tid = omp_get_thread_num() + 1
        qglobal(tid) = q 
   
        Mchi = wimp%mass 
        Mnuc = nuc_target%mass * mN
        muT = Mchi * Mnuc / (Mchi + Mnuc)
        ve = vearth * kilometerpersecond
        v0 = vscale * kilometerpersecond
        vesc = vescape * kilometerpersecond
        Nt = nuc_target%Nt
        rhochi = wimp%localdensity / centimeter**3d0
        vmin = q/(2d0*muT)
        relerror = quadrature_relerr
        vmax = vesc + ve
        !if (vmin > vmax) then
        !    deventrate = 0d0
        !    return
        !end if
    
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
            * dEventRate
      
    end function deventrate
    
    function spectraintegrand1d(vv, tid)
        ! The purpose of this function is to create a callable func for the
        ! adaptive integral routine. Other functions in this program use
        ! (the good practice of) explicit data dependency. An exception
        ! was made for this fuction for compatibility with the adaptive
        ! integration routine used to integrate the event rate spectra.
    
        implicit none
        real(kind=8) :: vv, qq
        real(kind=8) :: spectraintegrand1d
        integer tid
    
        qq = qglobal(tid)
    
        spectraintegrand1d = diffCrossSection(vv, qq) &
            * sshm(vv,ve,v0,vesc) * vv**2
    
    end function

end module eventrate
