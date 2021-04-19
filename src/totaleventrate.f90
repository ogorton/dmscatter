subroutine totaleventrate(nuc_target)
    use kinds
    use parameters
    implicit none
    type(nucleus) :: nuc_target
    real(doublep) :: q_start, q_stop, error, relerror, result, abserror
    real(doublep) :: fspectra, mtarget

    mtarget = nuc_target%mass

    print*,"Using adaptive numerical integration to determine &
         & total integrated event rate."
    call get_q_limits(mtarget, q_start, q_stop, relerror)

    print*,"Computing integral..."
    abserror=fspectra(q_stop) ! This line prevents a segfault... idk why

    call chinsp ( fspectra, q_start, q_stop, relerror, error, result )

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
    use integral
    implicit none
    real(doublep), intent(in) :: q
    real(doublep) :: fspectra

    fspectra = q * dEventRate(q, wimp, nuc_target, eft) /(mn*nuc_target%mass)

end function 

subroutine get_q_limits(mtarget,q_start, q_stop, error)
    use kinds
    use momenta
    use constants
    implicit none
    real(doublep) :: mtarget, E_start, E_stop, error
    real(doublep) :: q_start, q_stop

    ! Get recoil energy grid from user or from file
    if (usemomentum) then
        print*,'What are the limits of integration for transfer momenta?'
        print*,'Enter starting momentum, stopping momentum, relative error:'
        read*,q_start, q_stop, error        
    else
        print*,'What are the limits of integration for recoil energies in kev?'
        print*,'Enter starting energy, stoping energy, relative error:'
        read*,E_start, E_stop, error
    end if

    E_start = E_start + 1e-12
    q_start = q_start + 1e-12

    if (.not.usemomentum) then
        q_start = sqrt(2d0*mtarget*mN*E_start*kev)
        q_stop =sqrt(2d0*mtarget*mN*E_stop*kev)
    end if

    print*,"q min  (gev/c)",q_start
    print*,"q max  (gev/c)",q_stop
    print*,"E min  (kev)",e_start
    print*,"E max  (kev)",e_stop
    print*,"Rel. error tolerance",error

end subroutine
