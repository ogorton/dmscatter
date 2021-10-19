!-----------------------------------------------------------------------------80
subroutine nucresponse_spectra(nuc_target)
    use momenta
    use mod_spectra
    use orbitals, only: bfm
    use parameters, only: nucleus, pndens
    use nucresponse, only: nucFormFactor, nucFormFactor_transform
    implicit none
    integer :: n_xvalues
    integer :: iunit
    integer :: ioperator, tau, tau_prime, iy
    integer :: iqq
    type(nucleus) :: nuc_target   

    real(doublep) :: yy, qq
    real(doublep), allocatable :: xlist(:), Wlist(:,:,:,:)
    real(doublep) :: mtarget
    real(doublep) :: tolerance

    tolerance = epsilon(qq)

    mtarget = nuc_target%mass
    call get_energy_grid(mtarget)

    
    allocate(Wlist(energy_grid_size, 8, 0:1, 0:1))

    if (pndens) then
      !$OMP parallel do private(qq, yy, tau, tau_prime, ioperator) schedule(dynamic, 10)
      do iqq = 1, energy_grid_size
        qq = momentum_grid(iqq)
        yy = (qq * bfm / 2d0) ** 2d0        
        do tau = 0, 1
          do tau_prime = 0, 1
            do ioperator = 1, 8
              Wlist(iqq, ioperator, tau, tau_prime) = 0.25d0 * nucFormFactor_transform(&
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
              Wlist(iqq, ioperator, tau, tau_prime) = 0.25d0 * nucFormFactor(&
                  tau, tau_prime, ioperator, yy, nuc_target%densitymats%rho,&
                  nuc_target%groundstate%Tx2, nuc_target%Mt )
            end do
          end do
        end do
      end do
    end if      

    open(newunit=iunit, file='nucresponse_spectra.dat')
    write(iunit,'(a,i3)')'# Z = ',nuc_target%Z
    write(iunit,'(a,i3)')'# N = ',nuc_target%N
    write(iunit,'(a)')'# Nuclear response functions W(operator,tau, tau_prime, q)'
    write(iunit,'(a)')'# operator = 1, ..., 8:'
    write(iunit,'(a)')"# M, \Phi'', \tilde{\Phi}', \Delta, \Sigma', \Sigma'', \Phi''M, \Delta \Sigma'"  
    write(iunit,'(a)')'# tau = 0, 1 isospin coupling'
    write(iunit,'(a)')'# q = momentum transfer'
    write(iunit,'(a)')  '# Format:'
    write(iunit,'(a)')  '# q, { W(i, 0, 0, q) W(i, 1, 0, q) W(i, 0, 1, q) W(i, 1, 1, q); i = 1, ..., 8}'

    do iqq = 1, energy_grid_size
      qq = momentum_grid(iqq)
      write(iunit,'(ES15.5E3,32(ES15.5E3))') qq, Wlist(iqq, 1:8, 0:1, 0:1)
    end do
    close(iunit)

    deallocate(Wlist)

end subroutine nucresponse_spectra
