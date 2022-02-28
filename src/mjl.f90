module mjl

    ! Note: j values are integer number of half-quanta, i.e. j = 2 x j.
    ! bigJ values are not multiplied by 2, nor are l or n values.

    use norm
    use wigner
    use bessel
    use constants, only: pi
    contains

    function MJ(y,np,lp,jp,n,l,j,bigJ)
        implicit none
        integer, intent(in) :: np,lp,jp,n,l,j,bigJ
        real(kind=8), intent(in) :: y
        real(kind=8) :: Pi = 3.1415926535897932
        real(kind=8) :: MJ
    
        MJ = (-1d0)**(0.5d0+j/2d0+dble(bigJ))&
            * SQRT(J2norm(j)*J2norm(jp) &
            * Jnorm(l)*Jnorm(lp)*Jnorm(bigJ)/(4*Pi)) &
              * threej_lookup(2*lp,2*bigJ,2*l,0,0,0) &
              * sixj_lookup(2*lp,jp,1,j,2*l,2*bigJ)  &
              * BesselElement(y,np,lp,n,l,bigJ)
    
    end function MJ
    
    
    function MJLDivQoverall(lp,jp,l,j,bigJ,bigL)
        implicit none
        integer, intent(in) :: lp,jp,l,j,bigJ,bigL
        real(kind=8) :: MJLDivQoverall
    
        MJLDivQoverall = (-1.0)**(bigL+(j+1)/2)* Qnorm(lp) * q2norm(jp) &
            * q2norm(j)*Qnorm(bigJ)*Qnorm(bigL) &
            * sixj_lookup(2*lp,jp,1,j,2*l,2*bigJ)/SQRT(4d0*Pi)
    
    end function MJLDivQoverall
    
    
    function MJLDivQsummand1(y,np,lp,n,l,bigJ,bigL)
        implicit none
        integer, intent(in) :: np,lp,n,l,bigJ,bigL
        real(kind=8), intent(in) :: y
        real(kind=8) :: MJLDivQsummand1
    
        MJLDivQsummand1 = -SQRT(l+1d0)*Qnorm(l+1) &
            * sixj_lookup(2*bigL,2,2*bigJ,2*l,2*lp,2*(l+1)) &
            * threej_lookup(2*lp,2*bigL,2*(l+1),0,0,0) &
            * BesselElementminus(y,np,lp,n,l,bigL)
    
    end function MJLDivQsummand1
    
    
    function MJLDivQsummand2(y,np,lp,n,l,bigJ,bigL)
        implicit none
        integer, intent(in) :: np,lp,n,l,bigJ,bigL
        real(kind=8), intent(in) :: y
        real(kind=8) :: MJLDivQsummand2
    
        if (l.eq.0) then ! O.C.G: avoid Qnorm(-1)
            MJLDivQsummand2 = 0d0
        else
            MJLDivQsummand2 = SQRT(dble(l))*Qnorm(l-1) &
                * sixj_lookup(2*bigL,2,2*bigJ,2*l,2*lp,2*(l-1))  &
                * threej_lookup(2*lp,2*bigL,2*(l-1),0,0,0) &
                * BesselElementplus(y,np,lp,n,l,bigL)
        endif
    
    end function MJLDivQsummand2
    
    
    function MJLDivQ(y,np,lp,jp,n,l,j,bigJ,bigL)
        implicit none
        integer, intent(in) :: np,lp,jp,n,l,j,bigJ,bigL
        real(kind=8), intent(in) :: y
        real(kind=8) :: MJLDivQ
    
        MJLDivQ = MJLDivQoverall(lp,jp,l,j,bigJ,bigL) &            
            * (MJLDivQsummand1(y,np,lp,n,l,bigJ,bigL) &
            + MJLDivQsummand2(y,np,lp,n,l,bigJ,bigL))
    
    end function MJLDivQ

end module mjl
