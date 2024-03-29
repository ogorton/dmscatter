module spectra

    use main
    use crosssection
    use transition
    use constants, only: kev, mN, kilometerpersecond
    use eventrate, only: deventrate
    use quadrature, only: boole
    use settings
    use orbitals, only: bfm
    use nucresponse, only: nucFormFactor, nucFormFactor_transform
    use densities, only: pndens
        
    implicit none

    real(kind=8) :: x_start
    real(kind=8) :: x_stop
    real(kind=8) :: x_step
    real(kind=8), allocatable :: x_grid(:), energy_grid(:), momentum_grid(:)
    integer :: energy_grid_size

    contains

    function velocitycurve(vlist, q, option)
        implicit none
        real(kind=8) :: q, v
        real(kind=8), dimension(:) :: vlist
        integer :: option
    
        real(kind=8), dimension(size(vlist)) :: velocitycurve
    
        integer :: i
    
        select case(option)
        case(3)
            do i = 1, size(vlist)
                v = vlist(i)
                velocitycurve(i) = transition_probability(v, q)
            end do
        case(2)
            do i = 1, size(vlist)
                v = vlist(i)
                velocitycurve(i) = diffCrossSection(v, q)
            end do        
        case default
            stop "Not a velocity curve option."
        end select
    end function velocitycurve
    
    
    subroutine velocity_curve(option)
        implicit none
    
        real(kind=8) :: Er, Qr, vstart, vstop, vstep
        integer :: sizevlist, i
        real(kind=8), allocatable :: vlist(:), cslist(:)
        integer :: option, iunit
    
        print*,"Enter recoil E (keV):"
        read*,Er
        print '("Er = ",F8.4," (keV)")',Er
        Qr = sqrt(2d0*nuc_target%mass*mN*er*kev)
        print '("Qr = ",F8.4,"(GeV/c)")',qr
        print '("Vmin = ",F8.4,"(km/s)")',1/kilometerpersecond&
            *qr/(2.0* wimp%mass * nuc_target%mass * mN / (wimp%mass + nuc_target%mass * mN))
    
        print*,"Enter start v (km/s):"
        read*,vstart
        print '("v start = ",F8.4," (km/s)")',vstart
        vstart = vstart * kilometerpersecond
        print '("v start = ",ES12.4," (c)")',vstart
    
        print*,"Enter stop v (km/s):"
        read*,vstop
        print '("v stop = ",F12.4," (km/s)")',vstop
        vstop = vstop * kilometerpersecond
        print '("v stop = ",ES12.4," (c)")',vstop
    
        print*,"Enter step v (km/s):"
        read*,vstep
        print '("v step = ",F8.4," (km/s)")',vstep
        vstep = vstep * kilometerpersecond
        print '("v step = ",ES12.4," (c)")',vstep
    
        sizevlist = int(abs(vstop - vstart)/vstep)
        allocate(vlist(sizevlist))
        allocate(cslist(sizevlist))
        do i = 1, sizevlist
            vlist(i) = vstart + (i-1) * vstep
        end do
    
    
        cslist = velocitycurve(vlist, qr, option)
        select case(option)
        case(3)
            if (outfile=='default') outfile = "transition_probability.dat"
            open(newunit=iunit, file=trim(outfile))
            do i = 1, sizevlist
                write(iunit,*) vlist(i)/kilometerpersecond, cslist(i)
            end do
            close(iunit)
        case(2)
            if (outfile=='default') outfile = "crosssection.dat"
            open(newunit=iunit, file=trim(outfile))
            do i = 1, sizevlist
                write(iunit,*) vlist(i)/kilometerpersecond, cslist(i)
            end do
            close(iunit)
        case default
            STOP "Not a velocity curve option."
        end select        
        
    end subroutine
    
    subroutine eventrate_spectra
    
        implicit none
    
        integer :: calc_num
        integer :: iunit
        integer :: i, N
        real(kind=8) :: q, dq

        real(kind=8) :: recoil_energy, momentum_transfer
        real(kind=8), allocatable :: event_rate_spectra(:)
        real(kind=8) :: mtarget, totaleventrate
    
        print*,"Computing differential event rate spectra"
    
        ! Get parameters
        mtarget = nuc_target%mass
    
        ! Get domain
        call get_energy_grid(mtarget)
    
        ! Setup calculation and compute
        print*,'Number of event rates to compute:',energy_grid_size
        allocate(event_rate_spectra(energy_grid_size))
    
        N = size(momentum_grid)
        !$OMP parallel do private(q) schedule(dynamic, 1)
        do i = 1, N
            q = momentum_grid(i)
            Event_rate_spectra(i) = dEventRate(q)
        end do
        !$OMP end parallel do
        !$OMP barrier        
    
        if (.not.useenergyfile) then
            dq = abs(momentum_grid(2) - momentum_grid(1))
            call boole(energy_grid_size, momentum_grid*Event_rate_spectra/(mn*mtarget),&
                dq, totaleventrate)
            print*,'Total integrated eventrate (events)',totaleventrate
        end if 
    
        write(6,"(A,T30,A,T56,A)")" E-recoil (kev)","q-transfer (gev/c)","Eventrate (events/gev)"
        do calc_num = 1, energy_grid_size
            momentum_transfer = momentum_grid(calc_num)
            recoil_energy = energy_grid(calc_num)
            print*,recoil_energy,momentum_transfer,event_rate_spectra(calc_num)
        end do
    
        ! Write results to file
        if (outfile=='default') outfile = "eventrate_spectra.dat"
        open(newunit=iunit, file=trim(outfile))
        write(iunit,"(A,T30,A)")"# Recoil energy (kev)","Event rate (events/gev)"
        do calc_num = 1, energy_grid_size
            recoil_energy = energy_grid(calc_num)
            write(iunit,*)recoil_energy,event_rate_spectra(calc_num)
        end do
    
        close(iunit)
        print*,"Event rate spectra written to eventrate_spectra.dat"
    
    end subroutine eventrate_spectra
    
    subroutine get_energy_grid(mtarget)
        implicit none
        integer :: calc_num
        character(len=100) :: filename
        real(kind=8) :: momentum_transfer, mtarget, recoil_energy
    
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
    
            read*,x_start, x_stop, x_step
    
            if (usemomentum) then
                print*,"q min  (gec/c)",x_start
                print*,"q max  (gev/c)",x_stop
                print*,"q step (gev/c)",x_step
            else
                print*,"E min  (kev)",x_start
                print*,"E max  (kev)",x_stop
                print*,"E step (kev)",x_step
            end if        
    
            energy_grid_size = int((x_stop - x_start) / x_step) + 1
    
            allocate(energy_grid(energy_grid_size))
            do calc_num = 1, energy_grid_size
                energy_grid(calc_num) = x_start + (calc_num - 1) * x_step
            end do
        end if
    
        allocate(momentum_grid(energy_grid_size))
        do calc_num = 1, energy_grid_size
            if (usemomentum) then
                momentum_transfer = x_start + (calc_num - 1) * x_step
                recoil_energy = momentum_transfer**2d0 / (2d0*mtarget*mN*kev)
                momentum_grid(calc_num) = momentum_transfer
                energy_grid(calc_num) = recoil_energy
            else
                recoil_energy = energy_grid(calc_num)
                momentum_grid(calc_num) = sqrt(2d0*mtarget*mN*recoil_energy*kev)
            end if
        end do    
    end subroutine
    
    subroutine read_energy_grid(filename)
    
        
        implicit none
        character(len=100) :: filename
        integer :: i, io, iunit
    
        open(newunit=iunit,file=trim(filename))
        do 
            read(iunit,*,iostat=io)
            if (io/=0) exit
            energy_grid_size = energy_grid_size + 1
        end do
    
        allocate(energy_grid(energy_grid_size))
        rewind(iunit)
    
        do i = 1, energy_grid_size
            read(iunit,*) energy_grid(i)
        end do
        close(iunit)
    end subroutine read_energy_grid
    
    !-----------------------------------------------------------------------------80
    subroutine nucresponse_spectra

        implicit none
        integer :: iunit
        integer :: ioperator, tau, tau_prime
        integer :: iqq
    
        real(kind=8) :: yy, qq, xx
        real(kind=8), allocatable :: Wlist(:,:,:,:)
        real(kind=8) :: mtarget
        real(kind=8) :: tolerance, factor
        logical :: dotransform
    
        tolerance = epsilon(qq)
    
        mtarget = nuc_target%mass
        call get_energy_grid(mtarget)
    
        
        allocate(Wlist(energy_grid_size, 8, 0:1, 0:1))

        dotransform = (pndens .neqv. pnresponse)
        if (pndens .and. dotransform) then
          factor = 0.25d0
        else
          factor = 1d0
        end if
    
        if (dotransform) then
          !$OMP parallel do private(qq, yy, tau, tau_prime, ioperator) schedule(dynamic, 10)
          do iqq = 1, energy_grid_size
            qq = momentum_grid(iqq)
            yy = (qq * bfm / 2d0) ** 2d0        
            do tau = 0, 1
              do tau_prime = 0, 1
                do ioperator = 1, 8
                  Wlist(iqq, ioperator, tau, tau_prime) = factor * nucFormFactor_transform(&
                      tau, tau_prime, ioperator, yy, nuc_target%densitymats%rho,&
                      nuc_target%groundstate%Tx2, nuc_target%Mt )
                end do
              end do
            end do
          end do
        else
          !$OMP parallel do private(qq, yy, tau, tau_prime, ioperator) schedule(dynamic, 10)
          do iqq = 1, energy_grid_size
            qq = momentum_grid(iqq)
            yy = (qq * bfm / 2d0) ** 2d0
            do tau = 0, 1
              do tau_prime = 0, 1
                do ioperator = 1, 8
                  Wlist(iqq, ioperator, tau, tau_prime) = factor * nucFormFactor(&
                      tau, tau_prime, ioperator, yy, nuc_target%densitymats%rho,&
                      nuc_target%groundstate%Tx2, nuc_target%Mt )
                end do
              end do
            end do
          end do
        end if      
    
        if (outfile=='default') outfile = 'nucresponse_spectra.dat'
        open(newunit=iunit, file=trim(outfile))
        write(iunit,'(a,i3)')'# Z = ',nuc_target%Z
        write(iunit,'(a,i3)')'# N = ',nuc_target%N
        write(iunit,'(a)')'# Nuclear response functions W(operator,tau, tau_prime, q)'
        write(iunit,'(a)')'# operator = 1, ..., 8:'
        write(iunit,'(a)')"# M, \Phi'', \tilde{\Phi}', \Delta, \Sigma', \Sigma'', \Phi''M, \Delta \Sigma'"  
        write(iunit,'(a)')'# tau = 0, 1 isospin coupling'
        if (usemomentum) then
            write(iunit,'(a)')'# x = momentum transfer (gev/c)'
        else
            write(iunit,'(a)')'# x = recoil energy (kev)'
        end if
        write(iunit,'(a)')  '# Format:'
        write(iunit,'(a)')  '# x, { W(i, 0, 0, x) W(i, 1, 0, x) W(i, 0, 1, x) W(i, 1, 1, x); i = 1, ..., 8}'
    
        do iqq = 1, energy_grid_size
          if (usemomentum) then
              xx = momentum_grid(iqq)
          else
              xx = energy_grid(iqq)
          end if
            write(iunit,'(ES15.5E3,32(ES15.5E3))') xx,  Wlist(iqq, 1, 0:1, 0:1),&
                                                        Wlist(iqq, 2, 0:1, 0:1),&
                                                        Wlist(iqq, 3, 0:1, 0:1),&
                                                        Wlist(iqq, 4, 0:1, 0:1),&
                                                        Wlist(iqq, 5, 0:1, 0:1),&
                                                        Wlist(iqq, 6, 0:1, 0:1),&
                                                        Wlist(iqq, 7, 0:1, 0:1),&
                                                        Wlist(iqq, 8, 0:1, 0:1)
          print*,Wlist(iqq, 1, 0, 0),Wlist(iqq, 1, 1, 0),Wlist(iqq, 1, 0, 1),Wlist(iqq, 1, 1, 1)
        end do
        close(iunit)
    
        deallocate(Wlist)
    end subroutine nucresponse_spectra

end module spectra
