module phi
    use norm
    use bessel
    use wigner
    use sigma, only: SigmaJ
    contains

    function PhiPPJ(y,np,lp,jp,n,l,j,bigJ) result(operatorME)
    
        implicit none
        integer, intent(in) :: np,lp,jp,n,l,j,bigJ
        real(kind=8), intent(in) :: y
    
        real(kind=8) :: operatorME
    
        operatorMe = Qnorm(bigJ+1)*SQRT(bigJ+1d0) &
            * (PhiPPsummand1(y,np,lp,jp,n,l,j,bigJ) &
            + PhiPPsummand2(y,np,lp,jp,n,l,j,bigJ))

        if (bigJ .ne. 0) then
            operatorMe = operatorMe + Qnorm(bigJ-1)*SQRT(dble(bigJ)) &
                 & *(PhiPPsummand3(y,np,lp,jp,n,l,j,bigJ)    &
                 & +PhiPPsummand4(y,np,lp,jp,n,l,j,bigJ))
        end if
        operatorME = operatorME * PhiPPoverall(lp,jp,j)    
    
    end function PhiPPJ
   

    function PhiTPJ(y,np,lp,jp,n,l,j,bigJ) result(operatorME)

        implicit none
        integer, intent(in) :: np,lp,jp,n,l,j,bigJ
        real(kind=8), intent(in) :: y
        real(kind=8) :: operatorME
    
        operatorME = -Qnorm(bigJ+1)*SQRT(dble(bigJ)) &
            * (PhiPPsummand1(y,np,lp,jp,n,l,j,bigJ) &
            + PhiPPsummand2(y,np,lp,jp,n,l,j,bigJ))
        if (bigJ .ne. 0) then
            operatorMe = operatorMe + Qnorm(bigJ-1)*SQRT(bigJ+1d0) &
                 & *(PhiPPsummand3(y,np,lp,jp,n,l,j,bigJ)    &
                 & +PhiPPsummand4(y,np,lp,jp,n,l,j,bigJ))
        end if
        operatorME = operatorME * PhiPPoverall(lp,jp,j)
        operatorME = operatorME + SigmaJ(y,np,lp,jp,n,l,j,bigJ,bigJ)/2d0    
    
    end function PhiTPJ    
    

    function PhiPPsummand1(y,np,lp,jp,n,l,j,bigJ)
    
        implicit none
        integer, intent(in) :: np,lp,jp,n,l,j,bigJ
        real(kind=8), intent(in) :: y
        integer :: bigL
        real(kind=8) :: PhiPPsummand1
    
        PhiPPsummand1 = 0.D0
    
        do bigL = bigJ,bigJ+1
            PhiPPsummand1= PhiPPsummand1+(-1d0)**(bigJ-bigL+1d0)*Jnorm(bigL) &
                    * sixj_lookup(2*(bigJ+1),2,2*bigL,2,2*bigJ,2) &
                    * sixj_lookup(2*(bigJ+1),2,2*bigL,2*l,2*lp,2*(l+1)) &
                    * ninej(2*lp,2*l,2*bigL,1,1,2,jp,j,2*bigJ)
        end do
    
        PhiPPsummand1= PhiPPsummand1*Qnorm(l + 1)*SQRT(l + 1d0) &
             * threej_lookup(2*lp,2*(bigJ+1),2*(l+1),0,0,0) &
             * BesselElementminus(y,np,lp,n,l,bigJ+1)
    
    end function PhiPPsummand1
    
    
    function PhiPPsummand2(y,np,lp,jp,n,l,j,bigJ)
    
        implicit none
        integer, intent(in) :: np,lp,jp,n,l,j
        integer, intent(in) :: bigJ
        integer :: bigL
        real(kind=8), intent(in) :: y
        real(kind=8) :: PhiPPsummand2
    
        PhiPPsummand2 = 0.0
    
        if (l .ne. 0) then
    
            do bigL = bigJ,bigJ+1
                PhiPPsummand2= PhiPPsummand2+(-1.0)**(bigJ-bigL)*Jnorm(bigL) &
                    * sixj_lookup(2*(bigJ+1),2,2*bigL,2,2*bigJ,2) &
                    * sixj_lookup(2*(bigJ+1),2,2*bigL,2*l,2*lp,2*(l-1)) &
                    * ninej(2*lp,2*l,2*bigL,1,1,2,jp,j,2*bigJ)
            end do
    
            PhiPPsummand2= PhiPPsummand2 * SQRT(dble(l)) &
                * threej_lookup(2*lp,2*(bigJ+1),2*(l-1),0,0,0) &
                * BesselElementplus(y,np,lp,n,l,bigJ+1)
    
            if (l.ne.0) PhiPPsummand2= PhiPPsummand2*Qnorm(l - 1)
    
        end if
    
    end function PhiPPsummand2
    
    
    function PhiPPsummand3(y,np,lp,jp,n,l,j,bigJ)
        implicit none
        integer, intent(in) :: np,lp,jp,n,l,j
        integer, intent(in) :: bigJ
        integer :: bigL
        real(kind=8), intent(in) :: y
        real(kind=8) :: PhiPPsummand3
    
        PhiPPsummand3 = 0.0
    
        if (bigJ .ne. 0) then
            do bigL = bigJ-1,bigJ
                PhiPPsummand3= PhiPPsummand3+(-1.0)**(bigJ-bigL+1.0)*Jnorm(bigL)&
                     * sixj_lookup(2*(bigJ-1),2,2*bigL,2,2*bigJ,2) &
                     * sixj_lookup(2*(bigJ-1),2,2*bigL,2*l,2*lp,2*(l-1)) &
                     * ninej(2*lp,2*l,2*bigL,1,1,2,jp,j,2*bigJ)
            end do
            
            PhiPPsummand3= PhiPPsummand3 * SQRT(dble(l))
    
            if (l.ne.0) PhiPPsummand3 = PhiPPsummand3 * Qnorm(l - 1) &
                    * threej_lookup(2*lp,2*(bigJ-1),2*(l-1),0,0,0) &
                    * BesselElementplus(y,np,lp,n,l,bigJ-1)
    
        end if
    
    end function PhiPPsummand3
    
    
    function PhiPPsummand4(y,np,lp,jp,n,l,j,bigJ)
    
        implicit none
        integer, intent(in) :: np,lp,jp,n,l,j,bigJ
        real(kind=8), intent(in) :: y
        integer :: bigL
        real(kind=8) :: PhiPPsummand4
    
        PhiPPsummand4 = 0.0
    
        if (bigJ .ne. 0 ) then
            do bigL = bigJ-1,bigJ
                PhiPPsummand4= PhiPPsummand4+(-1d0)**(bigJ-bigL)*Jnorm(bigL)&
                    * sixj_lookup(2*(bigJ-1),2,2*bigL,2,2*bigJ,2) &
                    * sixj_lookup(2*(bigJ-1),2,2*bigL,2*l,2*lp,2*(l+1)) &
                    * ninej(2*lp,2*l,2*bigL,1,1,2,jp,j,2*bigJ)
            end do
    
            PhiPPsummand4 = PhiPPsummand4*Qnorm(l + 1)*SQRT(l+1d0)&
                * threej_lookup(2*lp,2*(bigJ-1),2*(l+1),0,0,0) &
                * BesselElementminus(y,np,lp,n,l,bigJ-1)
    
        end if
    
    end function PhiPPsummand4
    
    
    function PhiPPoverall(lp,jp,j)
    
        implicit none
        integer, intent(in) :: lp,jp,j
        real(kind=8) :: Pi = 3.1415926535897932
        real(kind=8) :: PhiPPoverall
    
        PhiPPoverall = (-1)**(lp+1) * 6d0 * SQRT(J2norm(j)*J2norm(jp) &
            * Jnorm(lp)/(4d0*Pi))
    
    end function PhiPPoverall


end module phi
