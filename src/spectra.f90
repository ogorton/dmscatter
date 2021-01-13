module espectra

    use kinds
    implicit none

    real(doublep) :: ER_start
    real(doublep) :: ER_stop
    real(doublep) :: ER_step
    real(doublep), allocatable :: energy_grid(:)
    integer :: energy_grid_size

contains        

subroutine eventrate_spectra(wimp, nuc_target, eft)
    use constants
    use momenta
    use parameters

    implicit none

    interface
        function EventRate(q, wimp, nuc_target, eft)
            use kinds
            use parameters
            implicit none
            real(doublep) :: q
            type(particle) :: wimp
            type(nucleus) :: nuc_target
            type(eftheory) :: eft
            real(doublep) :: eventrate
        end function eventrate
    end interface

    type(particle) :: wimp
    type(nucleus) :: nuc_target
    type(eftheory) :: eft
    
    integer :: calc_num, num_calc
    real(doublep) :: recoil_energy, momentum_transfer
    real(doublep), allocatable :: event_rate_spectra(:)
    character(len=40) :: filename

    real(doublep) :: mchi, jchi
    real(doublep) :: mtarget, jtarget
    real(doublep) :: Nt
    real(doublep) :: rhochi 

    ! Get parameters
    mchi = wimp%mass
    jchi = wimp%j
    rhochi = wimp%localdensity
    mtarget = nuc_target%mass
    jtarget = nuc_target%groundstate%jx2
    nt = nuc_target%nt

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
    
    write(6,"(A,T30,A,T56,A)")" E-recoil (kev)","q-transfer (gev/c)","Eventrate (events/gev)"
    do calc_num = 1, energy_grid_size
        if (usemomentum) then
            momentum_transfer = ER_start + (calc_num - 1) * ER_step
            recoil_energy = momentum_transfer**2d0 / (2d0*mtarget*mN*kev)
            energy_grid(calc_num) = recoil_energy
        else
            recoil_energy = energy_grid(calc_num)
            momentum_transfer = sqrt(2d0 * mtarget*mN * recoil_energy * kev)
        end if
        event_rate_spectra(calc_num) = EventRate(momentum_transfer,&
                wimp, nuc_target, eft)
        print*,recoil_energy,momentum_transfer,event_rate_spectra(calc_num)
    end do

    ! Write results to file
    open(unit=157, file='eventrate_spectra.dat')
    write(157,"(A,T30,A)")"# Recoil energy (kev)","Event rate (events/gev)"
    do calc_num = 1, energy_grid_size
        recoil_energy = energy_grid(calc_num)
        write(157,*)recoil_energy,event_rate_spectra(calc_num)
    end do

    close(157)

    print*,"Event rate spectra written to eventrate_spectra.dat"

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

