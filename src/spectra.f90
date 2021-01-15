module espectra

    use kinds
    implicit none

    real(doublep) :: ER_start
    real(doublep) :: ER_stop
    real(doublep) :: ER_step
    real(doublep), allocatable :: energy_grid(:), momentum_grid(:)
    integer :: energy_grid_size

contains        

function spectra(momenta, wimp, nuc_target, eft)
    use parameters
    implicit none
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
    real(doublep), dimension(:) :: momenta
    type(particle) :: wimp
    type(nucleus) :: nuc_target
    type(eftheory) :: eft
    real(doublep), dimension(size(momenta)) :: spectra

    integer :: i, N
    real(doublep) :: q

    N = size(momenta)
    do i = 1, N
        q = momenta(i)
        spectra(i) = dEventRate(q, wimp, nuc_target, eft)
    end do

end function spectra

subroutine eventrate_spectra(wimp, nuc_target, eft)
    use constants
    use momenta
    use parameters

    implicit none

    type(particle) :: wimp
    type(nucleus) :: nuc_target
    type(eftheory) :: eft
    
    integer :: calc_num
    real(doublep) :: recoil_energy, momentum_transfer
    real(doublep), allocatable :: event_rate_spectra(:)

    real(doublep) :: mtarget, totaleventrate, error

    print*,"Computing differential event rate spectra"

    ! Get parameters
    mtarget = nuc_target%mass

    ! Get domain
    call get_energy_grid
    allocate(momentum_grid(energy_grid_size))
    do calc_num = 1, energy_grid_size
        if (usemomentum) then
            momentum_transfer = ER_start + (calc_num - 1) * ER_step
            recoil_energy = momentum_transfer**2d0 / (2d0*mtarget*mN*kev)
            momentum_grid(calc_num) = momentum_transfer
            energy_grid(calc_num) = recoil_energy
        else
            recoil_energy = energy_grid(calc_num)
            momentum_grid(calc_num) = sqrt(2d0*mtarget*mN*recoil_energy*kev)
        end if
    end do    

    ! Setup calculation and compute
    print*,'Number of event rates to compute:',energy_grid_size
    allocate(event_rate_spectra(energy_grid_size))

    event_rate_spectra = spectra(momentum_grid, wimp, nuc_target, eft)

    call cubint ( energy_grid_size, momentum_grid, momentum_grid*Event_rate_spectra/(mn*mtarget), 1, &
                energy_grid_size, totaleventrate, error )

    write(6,"(A,T30,A,T56,A)")" E-recoil (kev)","q-transfer (gev/c)","Eventrate (events/gev)"
    do calc_num = 1, energy_grid_size
        momentum_transfer = momentum_grid(calc_num)
        recoil_energy = energy_grid(calc_num)
        print*,recoil_energy,momentum_transfer,event_rate_spectra(calc_num)
    end do

    print*,'Total integrated eventrate (events)',totaleventrate,'pm',error

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

subroutine get_energy_grid
    use momenta
    implicit none
    integer :: calc_num
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

        if (usemomentum) then
            print*,"q min  (gec/c)",ER_start
            print*,"q max  (gev/c)",ER_stop
            print*,"q step (gev/c)",ER_step
        else
            print*,"E min  (kev)",ER_start
            print*,"E max  (kev)",ER_stop
            print*,"E step (kev)",ER_step
        end if        

        energy_grid_size = int((ER_stop - ER_start) / ER_step) + 1

        allocate(energy_grid(energy_grid_size))
        do calc_num = 1, energy_grid_size
            energy_grid(calc_num) = ER_start + (calc_num - 1) * ER_step
        end do
    end if    
end subroutine

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

