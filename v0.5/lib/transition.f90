function transition_probability(q,v,jchi,y,Mtiso)

    use dmresponse
    use stateinfo 
    implicit none
    interface 
        function nucResponse(tau1,tau2,ioption,y,Mtiso)
            integer :: tau1, tau2
            integer :: ioption
            real(kind=8) :: y
            integer :: Mtiso
            REAL(kind=8) :: nucResponse
        end function
        function dmresponsecoef(ifunc, tau1, tau2, q, v, jchi)
            integer :: ifunc
            integer :: tau1, tau2
            real(kind=8) :: q, v, jchi
            REAL(kind=8) :: dmresponsecoef
        end function
    end interface
    REAL(kind=8), INTENT(IN) :: q
    REAL(kind=8), INTENT(IN) :: v
    REAL(kind=8), INTENT(IN) :: jchi
    REAL(kind=8), INTENT(IN) :: y
    integer, intent(in) :: Mtiso
    REAL(kind=8) :: transition_probability
    real(kind=8) :: pi = 3.14159265358979, tmp

    integer :: tau1, tau2, ifunc

    transition_probability = 0.0

    do tau1 = 0, 1
        do tau2 = 0, 1
            do ifunc = 1, 8
                tmp = dmresponsecoef(ifunc, tau1, tau2, q, v, jchi)
                print*,'tmp',tmp
                transition_probability = transition_probability &
                    + dmresponsecoef(ifunc, tau1, tau2, q, v, jchi) &
                    * nucResponse(tau1,tau2,ifunc,y,Mtiso) 
            end do ! ifunc
        end do ! tau2
    end do ! tau1
    transition_probability = transition_probability*4.0*pi/(Jiso+1.0)
    
end function transition_probability
