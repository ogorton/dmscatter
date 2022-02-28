module transition
    use constants
    use orbitals, only: bfm
    use dmresponse, only: dmresponsefun
    use nucresponse
    use main

    contains
    function transition_probability(v, q)
        
        implicit none
        REAL(kind=8) :: q
        REAL(kind=8) :: v
        real(kind=8), allocatable :: eftsmall(:,:)
        !
        REAL(kind=8) :: y
        REAL(kind=8) :: muT
        REAL(kind=8) :: transition_probability, tmpprod
    
        integer :: tau1, tau2, term
    
        allocate(eftsmall(0:1,size(eft%xpnc(0)%c)))
        eftsmall(0,:) = eft%xpnc(0)%c
        eftsmall(1,:) = eft%xpnc(1)%c    
    
        y = (q*bfm/2d0)**2d0
        muT = wimp%mass * mN * nuc_target%mass / (wimp%mass + mN * nuc_target%mass)
    
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
                            nuc_target%densitymats%rho, nuc_target%groundstate%Tx2, nuc_target%Mt)
    
                        transition_probability = transition_probability + tmpprod
                    end do ! term
                end do ! tau2
            end do ! tau1
    
            transition_probability = transition_probability &
                * (4d0*pi/(nuc_target%groundstate%jx2+1d0)) / ((4d0*mN*wimp%mass)**2)
    
        end if
    end function transition_probability
end module
