function transition_probability(q,v,wimp,nucl,eft)
    use kinds
    use constants
    use parameters
    use spspace, only: bfm
    implicit none
    interface 
        function nucResponse(tau1,tau2,ioption,y,densmat,Tiso,Mtiso)
            use kinds
            use parameters
            integer :: tau1, tau2
            integer :: ioption
            real(doublep) :: y
            integer :: Mtiso, Tiso
            real(doublep), allocatable :: densmat(:,:,:,:)
            REAL(doublep) :: nucResponse
        end function
        function dmresponsefun(eft, term, tau1, tau2, q, v, jchi, muT)
            use kinds
            use parameters
            real(doublep), allocatable, intent(in) :: eft(:,:)
            integer :: term
            integer :: tau1, tau2
            real(doublep) :: q, v, jchi, muT
            REAL(doublep) :: dmresponsefun
        end function
    end interface
    REAL(doublep) :: q
    REAL(doublep) :: v
    type(particle) :: wimp
    type(nucleus) :: nucl
    real(kind=8), allocatable, intent(in) :: eft(:,:)
    !
    REAL(doublep) :: y
    REAL(doublep) :: mchi, jchi, mtarget, jtarget
    REAL(doublep) :: muT
    REAL(doublep) :: transition_probability

    integer :: tau1, tau2, term
    integer :: Mtiso, Tiso

    bfm = (41.467/(45.*(nucl%mass)**(-1./3) &
                - 25.*(nucl%mass)**(-2./3)))**0.5 * femtometer
    y = (q*bfm/2d0)**2d0

    mchi = wimp%mass
    jchi = wimp%j
    mtarget = nucl%mass
    Tiso = nucl%groundstate%Tx2
    Mtiso = nucl%Mt
    jtarget = nucl%groundstate%jx2

    muT = mchi * mtarget * mN / (mchi+mtarget*mN)

    transition_probability = 0.0

    do tau1 = 0, 1
        do tau2 = 0, 1
            do term = 1, 8
                transition_probability = transition_probability &
                    + dmresponsefun(eft, term, tau1, tau2, q, v, jchi, muT)&
                    * nucResponse(tau1,tau2,term,y,nucl%densitymats%rho,Tiso,Mtiso)
            end do ! term
        end do ! tau2
    end do ! tau1

    transition_probability = transition_probability &
        * (4.0*pi/(jtarget+1)) / ((4.0*mN*mchi)**2.0)
    
end function transition_probability
