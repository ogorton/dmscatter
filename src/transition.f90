module transition
contains
function transition_probability(q,v,wimp,nucl,eft)
    use kinds
    use constants
    use orbitals, only: bfm
    use dmresponse, only: dmresponsefun
    use nucresponse
    use types
    implicit none
    REAL(dp) :: q
    REAL(dp) :: v
    type(particle) :: wimp
    type(nucleus) :: nucl
    type(eftheory) :: eft
    real(kind=8), allocatable :: eftsmall(:,:)
    !
    REAL(dp) :: y
    REAL(dp) :: muT
    REAL(dp) :: transition_probability, tmpprod

    integer :: tau1, tau2, term

    allocate(eftsmall(0:1,size(eft%xpnc(0)%c)))
    eftsmall(0,:) = eft%xpnc(0)%c
    eftsmall(1,:) = eft%xpnc(1)%c    

    y = (q*bfm/2d0)**2d0
    muT = wimp%mass * mN * nucl%mass / (wimp%mass + mN * nucl%mass)

    transition_probability = 0.0
    if (v <  q/(2d0*muT)) then
        transition_probability = 0.0
    else
        do tau1 = 0, 1
            do tau2 = 0, 1
                do term = 1, 8
                    tmpprod = dmresponsefun(eftsmall, term, tau1, tau2, &
                        q, v, wimp%j, muT)

                    if (tmpprod.eq.0) cycle

                    tmpprod = tmpprod * nucFormFactor(tau1, tau2, term, y, &
                        nucl%densitymats%rho, nucl%groundstate%Tx2, nucl%Mt)

                    transition_probability = transition_probability + tmpprod
                end do ! term
            end do ! tau2
        end do ! tau1

        transition_probability = transition_probability &
            * (4d0*pi/(nucl%groundstate%jx2+1d0)) / ((4d0*mN*wimp%mass)**2)

    end if
end function transition_probability
end module
