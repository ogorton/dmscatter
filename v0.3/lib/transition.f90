function transition_probability(q,v,jchi,y,Mtiso)
    use dmresponse
    
    implicit none
    interface 
        function nucResponse(tau1,tau2,ioption,y,Mtiso)
            integer :: tau1, tau2
            integer :: ioption
            real(kind=8) :: y
            integer :: Mtiso
            REAL(kind=8) :: nucResponse
        end function
    end interface
    REAL(kind=8), INTENT(IN) :: q
    REAL(kind=8), INTENT(IN) :: v
    REAL(kind=8), INTENT(IN) :: jchi
    REAL(kind=8), INTENT(IN) :: y
    integer, intent(in) :: Mtiso
    REAL(kind=8) :: transition_probability

    integer :: tau1, tau2

    transition_probability = 0.0

    do tau1 = 0, 1
        do tau2 = 0, 1
            transition_probability = transition_probability &
                + dmrMJ(tau1, tau2, q, v, jchi) * nucResponse(tau1,tau2,1,y,Mtiso) &
                + dmrPhiPPJ(tau1, tau2, q, v, jchi)* nucResponse(tau1,tau2,2,y,Mtiso) &
                + dmrPhiTPJ(tau1, tau2, q, v, jchi) * nucResponse(tau1,tau2,3,y,Mtiso) &
                + dmrDeltaJ(tau1, tau2, q, v, jchi) * nucResponse(tau1,tau2,4,y,Mtiso) &
                + dmrSigmaPJ(tau1, tau2, q, v, jchi) * nucResponse(tau1,tau2,5,y,Mtiso) &
                + dmrSigmaPPJ(tau1, tau2, q, v, jchi) * nucResponse(tau1,tau2,6,y,Mtiso) &
                + dmrPhiPPJMJ(tau1, tau2, q, v, jchi) * nucResponse(tau1,tau2,7,y,Mtiso) &
                + dmrSigmaPJDeltaJ(tau1, tau2, q, v, jchi) * nucResponse(tau1,tau2,8,y,Mtiso)
        end do ! tau2
    end do ! tau1
    
end function transition_probability
