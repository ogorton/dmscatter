subroutine totaleventrate
    use kinds
    use parameters
    implicit none
    real(doublep) :: ER_start, ER_stop, error, relerror, result, abserror
    integer :: ind
    real(doublep) :: fspectra, tmp

    print*,"Using adaptive numerical integration to determine &
         total integrated event rate."
    call get_energy_limits(ER_start, ER_stop, relerror)

    print*,"Computing integral..."
    abserror=fspectra(ER_stop) ! This line prevents a segfault... idk why

    call chinsp ( fspectra, ER_start, ER_stop, relerror, error, result )

    print*,"Error estimate:"
    print*,error

    print*,"Total integrated events:"
    print*,result

end subroutine totaleventrate

function fspectra(q)
    ! The purpose of this function is to create a callable func for the
    ! adaptive integral routine. Other functions in this program use
    ! (the good practice of) explicit data dependency. An exception
    ! was made for this fuction for compatibility with the adaptive
    ! integration routine used to integrate the event rate spectra.

    use main ! This is the only function allowed to use main.
    use kinds
    use constants
    implicit none
    real(doublep), intent(in) :: q
    real(doublep) :: fspectra
    interface
        function dEventRate(q, wimp, nuc_target, eft)
            use kinds
            use parameters
            implicit none
            real(doublep) :: q
            type(particle) :: wimp
            type(nucleus) :: nuc_target
            type(eftheory) :: eft
            real(doublep) :: deventrate
        end function deventrate
    end interface        

    fspectra = q * dEventRate(q, wimp, nuc_target, eft) /(mn*nuc_target%mass)

end function 

subroutine get_energy_limits(ER_start, ER_stop, error)
    use kinds
    use momenta
    implicit none
    real(doublep) :: ER_start, ER_stop, error
    character(len=40) :: filename

    ! Get recoil energy grid from user or from file
        if (usemomentum) then
            print*,'What are the limits of integration for transfer momenta?'
            print*,'Enter starting momentum, stopping momentum, relative error:'
        else
            print*,'What are the limits of integration for recoil energies in kev?'
            print*,'Enter starting energy, stoping energy, relative error:'
        end if

        read*,ER_start, ER_stop, error

        ER_start = ER_start + 1E-12

        if (usemomentum) then
            print*,"q min  (gev/c)",ER_start
            print*,"q max  (gev/c)",ER_stop
        else
            print*,"E min  (kev)",ER_start
            print*,"E max  (kev)",ER_stop
        end if
        print*,"Rel. error tolerance",error

end subroutine
