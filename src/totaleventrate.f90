module totaleventrate
    use quadrature
    use constants
    use eventrate
    use settings
    use main

    contains
    subroutine totaler()
        
        implicit none
        real(kind=8) :: q_start, q_stop, error, relerror, result, abserror
        real(kind=8) :: mtarget
        integer :: ind
    
        mtarget = nuc_target%mass
    
        print*,"Using adaptive numerical integration to determine &
             & total integrated event rate."
        call get_q_limits(mtarget, q_start, q_stop, relerror)
    
        print*,"Computing integral..."
        abserror=fspectra(q_stop) ! This line prevents a segfault... idk why
        relerror = quadrature_relerr
    
        call gaus8_threadsafe(fspectra,  q_start, q_stop, relerror, result, ind, 1)
    
        print*,"Error estimate:"
        print*,error
    
        print*,"Total integrated events:"
        print*,result
    
    end subroutine totaler
    
    function fspectra(q)
        ! The purpose of this function is to create a callable func for the
        ! adaptive integral routine. Other functions in this program use
        ! (the good practice of) explicit data dependency. An exception
        ! was made for this fuction for compatibility with the adaptive
        ! integration routine used to integrate the event rate spectra.
    
        implicit none
        real(kind=8), intent(in) :: q
        real(kind=8) :: fspectra
    
        fspectra = q * dEventRate(q) /(mn*nuc_target%mass)
    
    end function 
    
    subroutine get_q_limits(mtarget,q_start, q_stop, error)
        
        implicit none
        real(kind=8) :: mtarget, E_start, E_stop, error
        real(kind=8) :: q_start, q_stop
    
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
end module totaleventrate
