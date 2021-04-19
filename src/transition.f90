module mtransition_probability
    contains
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
    type(eftheory) :: eft
    real(kind=8), allocatable :: eftsmall(:,:)
    !
    REAL(doublep) :: y
    REAL(doublep) :: muT
    REAL(doublep) :: transition_probability, tmpprod

    integer :: tau1, tau2, term

    allocate(eftsmall(0:1,num_response_coef))
    eftsmall(0,:) = eft%xpnc(0)%c
    eftsmall(1,:) = eft%xpnc(1)%c    

    y = (q*bfm/2d0)**2d0
    muT = wimp%mass * mN * nucl%mass / (wimp%mass + mN * nucl%mass)

    transition_probability = 0.0
    if (v <  q/(2d0*muT)) return

    do tau1 = 0, 1
        do tau2 = 0, 1
            do term = 1, 8
                tmpprod = dmresponsefun(eftsmall, term, tau1, tau2, q, v, wimp%j, muT)
                if (tmpprod.ne.0) tmpprod = tmpprod * &
                    nucResponse(tau1, tau2, term, y, nucl%densitymats%rho, nucl%groundstate%Tx2, nucl%Mt)
                transition_probability = transition_probability + tmpprod
            end do ! term
        end do ! tau2
    end do ! tau1

    transition_probability = transition_probability &
        * (4.0*pi/(nucl%groundstate%jx2+1)) / ((4.0*mN*wimp%mass)**2.0)
    
end function transition_probability
end module
