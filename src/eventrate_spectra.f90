module espectra

    use kinds
    implicit none

    real(doublep) :: ER_start
    real(doublep) :: ER_stop
    real(doublep) :: ER_step

contains        

subroutine eventrate_spectra
    use masses
    use constants
    use momenta
    implicit none
    interface
        function eventrate(q)
            use kinds
            implicit none
            real(doublep) :: q
            real(doublep) :: eventrate
        end function eventrate
    end interface
    integer :: calc_num, num_calc
    real(doublep) :: recoil_energy, momentum_transfer
    real(doublep), allocatable :: event_rate_spectra(:)

    if (usemomentum) then
        print*,'What is the range of transfer momenta?'
        print*,'Enter starting momentum, stopping momentum, setp size:'
    else
        print*,'What is the range of recoil energies in kev?'
        print*,'Enter starting energy, stoping energy, step size:'
    end if

    read*,ER_start, ER_stop, ER_step
    num_calc = int((ER_stop - ER_start) / ER_step) + 1
    print*,'Number of event rates to compute:',num_calc
    allocate(event_rate_spectra(num_calc))

!!!$OMP parallel do private(recoil_energy, momentum_transfer)
    do calc_num = 1, num_calc
        if (usemomentum) then
            momentum_transfer = ER_start + (calc_num - 1) * ER_step
            recoil_energy = momentum_transfer**2d0 / (2d0*mtarget)
        else
            recoil_energy = ER_start + (calc_num - 1) * ER_step
            momentum_transfer = sqrt(2d0 * mtarget * recoil_energy * kev)
        end if
        event_rate_spectra(calc_num) = EventRate(momentum_transfer)
        print*,'E_recoil = ',recoil_energy, 'q=',momentum_transfer,'ER=',event_rate_spectra(calc_num)
    end do
!!!$OMP end parallel do
!!!$OMP barrier

    open(unit=157, file='eventrate_spectra.dat')
    write(157,*)'#Recoil energy (kev)    Event rate (Events/second/)'
    do calc_num = 1, num_calc
        recoil_energy = ER_start + (calc_num - 1) * ER_step
        write(157,*)recoil_energy,event_rate_spectra(calc_num)
    end do

    close(157)

    print*,"Event rate spectra written to evenrate_spectra.dat"

end subroutine eventrate_spectra

end module espectra

