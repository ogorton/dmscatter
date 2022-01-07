module sigma
    use norm
    use bessel
    use wigner
    contains

    function SigmaJ(y,np,lp,jp,n,l,j,bigJ,bigL)
        implicit none
        integer, intent(in) :: np,lp,jp,n,l,j
        integer, intent(in) :: bigJ,bigL
        real(kind=8), intent(in) :: y
        real(kind=8) :: Pi = 3.14159265358979323
        real(kind=8) :: SigmaJ
    
        SigmaJ = (-1d0)**lp * SQRT(J2norm(j)*j2norm(jp)*Jnorm(l)*Jnorm(lp)  &
              * Jnorm(bigJ) * Jnorm(bigL)/(4d0*Pi))     &
              * SQRT(6d0) * threej_lookup(2*lp,2*bigL,2*l,0,0,0) &
              * BesselElement(y, np, lp, n, l, bigL) &
              * ninej(2*lp, 2*l, 2*bigL, 1, 1, 2, jp, j, 2*bigJ) 
    
    end function SigmaJ
    
    
    function SigmaPJ(y,np,lp,jp,n,l,j,bigJ)
      integer, intent(in) :: np,lp,jp,n,l,j
      integer, intent(in) :: bigJ
      real(kind=8), intent(in) :: y
      real(kind=8) :: SigmaPJ
    
      SigmaPJ = -SQRT(dble(bigJ))/Qnorm(bigJ) &
          * SigmaJ(y,np,lp,jp,n,l,j,bigJ,bigJ+1) &
          + SQRT(bigJ+1d0)/Qnorm(bigJ) & 
          * SigmaJ(y,np,lp,jp,n,l,j,bigJ,bigJ-1)
    
    end function SigmaPJ
    
    
    function SigmaPPJ(y,np,lp,jp,n,l,j,bigJ)
        integer, intent(in) :: np,lp,jp,n,l,j
        integer, intent(in) :: bigJ
        real(kind=8), intent(in) :: y
        real(kind=8) :: SigmaPPJ
    
        SigmaPPJ = SQRT(bigJ+1d0)/Qnorm(bigJ) &
            * SigmaJ(y,np,lp,jp,n,l,j,bigJ,bigJ+1) &
            + SQRT(dble(bigJ))/Qnorm(bigJ) &
            * SigmaJ(y,np,lp,jp,n,l,j,bigJ,bigJ-1)
    
    end function SigmaPPJ

end module sigma
