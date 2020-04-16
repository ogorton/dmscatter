function transition_probability(q,v,jchi,y,Mtiso)
    use kinds
    use dmresponse
    use stateinfo 
    implicit none
    interface 
        function nucResponse(tau1,tau2,ioption,y,Mtiso)
            use kinds
            integer :: tau1, tau2
            integer :: ioption
            real(doublep) :: y
            integer :: Mtiso
            REAL(doublep) :: nucResponse
        end function
        function dmresponsecoef(ifunc, tau1, tau2, q, v, jchi)
            use kinds
            integer :: ifunc
            integer :: tau1, tau2
            real(doublep) :: q, v, jchi
            REAL(doublep) :: dmresponsecoef
        end function
    end interface
    REAL(doublep), INTENT(IN) :: q
    REAL(doublep), INTENT(IN) :: v
    REAL(doublep), INTENT(IN) :: jchi
    REAL(doublep), INTENT(IN) :: y
    integer, intent(in) :: Mtiso
    REAL(doublep) :: transition_probability
    real(doublep) :: pi = 3.14159265358979, tmp

    integer :: tau1, tau2, ifunc

    transition_probability = 0.0

    do tau1 = 0, 1
        do tau2 = 0, 1
            do ifunc = 1, 8
                tmp = dmresponsecoef(ifunc, tau1, tau2, q, v, jchi)
!                print*,'tmp',tmp
                transition_probability = transition_probability &
                    + dmresponsecoef(ifunc, tau1, tau2, q, v, jchi) &
                    * nucResponse(tau1,tau2,ifunc,y,Mtiso) 
            end do ! ifunc
        end do ! tau2
    end do ! tau1
    transition_probability = transition_probability*4.0*pi/(Jiso+1.0)
    
end function transition_probability
