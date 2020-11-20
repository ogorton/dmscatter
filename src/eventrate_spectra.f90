module espectra

    use kinds
    implicit none

    real(doublep) :: ER_start
    real(doublep) :: ER_stop
    real(doublep) :: ER_step
    real(doublep), allocatable :: energy_grid(:)
    integer :: energy_grid_size

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
    character(len=40) :: filename

    ! Get recoil energy grid from user or from file
    if (useenergyfile) then
        print*,'Recoil energies will be read from file. Enter filename:'
        read*,filename
        call read_energy_grid(filename)
    else
        if (usemomentum) then
            print*,'What is the range of transfer momenta?'
            print*,'Enter starting momentum, stopping momentum, setp size:'
        else
            print*,'What is the range of recoil energies in kev?'
            print*,'Enter starting energy, stoping energy, step size:'
        end if

        read*,ER_start, ER_stop, ER_step
        energy_grid_size = int((ER_stop - ER_start) / ER_step) + 1
        allocate(energy_grid(energy_grid_size))
        do calc_num = 1, energy_grid_size
            energy_grid(calc_num) = ER_start + (calc_num - 1) * ER_step
        end do
            
    end if

    ! Setup calculation and compute
    print*,'Number of event rates to compute:',energy_grid_size
    allocate(event_rate_spectra(energy_grid_size))
    
    print*,'E_recoil       q_transfer        ER'
    do calc_num = 1, energy_grid_size
        if (usemomentum) then
            momentum_transfer = ER_start + (calc_num - 1) * ER_step
            recoil_energy = momentum_transfer**2d0 / (2d0*mtarget)
        else
            recoil_energy = energy_grid(calc_num)
            momentum_transfer = sqrt(2d0 * mtarget * recoil_energy * kev)
        end if
        event_rate_spectra(calc_num) = EventRate(momentum_transfer)
        print*,recoil_energy,momentum_transfer,event_rate_spectra(calc_num)
    end do

    open(unit=157, file='eventrate_spectra.dat')
    write(157,*)'# Recoil energy (kev)    Event rate (Events/second/)'
    do calc_num = 1, num_calc
        if (usemomentum) then
            momentum_transfer = ER_start + (calc_num - 1) * ER_step
            recoil_energy = momentum_transfer**2d0 / (2d0*mtarget)
        else
            recoil_energy = ER_start + (calc_num - 1) * ER_step
            momentum_transfer = sqrt(2d0 * mtarget * recoil_energy * kev)
        end if        
        write(157,*)recoil_energy,event_rate_spectra(calc_num)
    end do

    close(157)

    print*,"Event rate spectra written to evenrate_spectra.dat"

end subroutine eventrate_spectra

subroutine read_energy_grid(filename)

    use kinds
    implicit none
    character(len=40) :: filename
    integer :: i, io

    open(unit=159,file=trim(filename))

    do 
        read(159,*,iostat=io)
        if (io/=0) exit
        energy_grid_size = energy_grid_size + 1
    end do

    allocate(energy_grid(energy_grid_size))

    rewind(159)

    do i = 1, energy_grid_size
        read(159,*) energy_grid(i)
    end do

    close(159)

end subroutine read_energy_grid

end module espectra

