module bessel
    use wigner
    contains
    !-------------------------------------------------------------------------------
    function BesselF1(y,np,lp,n,l,bigL)
        implicit none
    
        INTEGER, INTENT(IN) :: np,lp,n,l,bigL
        REAL(kind=8),INTENT(IN) :: y
        REAL(kind=8) :: BesselF1
    
        BesselF1 = 2.0**(dble(bigL)) / exp(logdoublefac(2*bigL+1)) &
            * y**(bigL/2.d0) &
            * exp(-y) * sqrt(EXP(logfac(np-1)+logfac(n-1)))
    
    end function BesselF1
    
    
    !-------------------------------------------------------------------------------
    function BesselF1_plusminus(y,np,lp,n,l,bigL)
        implicit none
    
        INTEGER, INTENT(IN) :: np,lp,n,l,bigL
        REAL(kind=8),INTENT(IN) :: y
        REAL(kind=8) :: BesselF1_plusminus
    
        BesselF1_plusminus = 2**(dble(bigL)-1.0) / exp(logdoublefac(2*bigL+1)) &
            * y**((dble(bigL)-1.0)/2.d0) &
            * exp(-y) * sqrt (EXP(logfac(np-1)+logfac(n-1)))
    
    end function BesselF1_plusminus
    
    
    !-------------------------------------------------------------------------------
    function BesselF2(np,lp,n,l)
      implicit none
      INTEGER, INTENT(IN) :: np,lp,n,l
      REAL(kind=8) :: BesselF2

      BesselF2 = sqrt(gamma(np+lp+0.5d0) * gamma(n+l+0.5d0))

    end function BesselF2
    
    
    !-------------------------------------------------------------------------------
    function BesselF3(y,np,lp,n,l,bigL)
        implicit none
    
        INTEGER, INTENT(IN) :: np,lp,n,l, bigL
        REAL(kind=8),INTENT(IN) :: y
        INTEGER :: m,mp!ip, i, m, mp
        REAL(kind=8) :: BesselF3
    
        BesselF3 = 0.0d0
    
        do m = 0, n-1
            do mp = 0, np-1
                BesselF3 = BesselF3 + Summand1(m,mp,n,np) &
                    * Gammas_regular(m,mp,l,lp,bigL) * CHGF_regular(y,m,mp,l,lp,bigL)
            end do
        end do
    
    end function BesselF3
    
    
    !-------------------------------------------------------------------------------
    function BesselF3_minus(y,np,lp,n,l,bigL)
        implicit none
    
        INTEGER, INTENT(IN) :: np,lp,n,l, bigL
        REAL(kind=8),INTENT(IN) :: y
        INTEGER :: mp, m
        REAL(kind=8) :: BesselF3_minus
    
        BesselF3_minus = 0.d0
    
        do m = 0, n-1
            do mp = 0, np-1
                BesselF3_minus = BesselF3_minus + Summand1(m,mp,n,np) &
                    * Gammas_plusminus(m,mp,l,lp,bigL) * CHGF_minus(y,m,mp,l,lp,bigL)
            end do
        end do
    
    end function BesselF3_minus
    
    
    !-------------------------------------------------------------------------------
    function BesselF3_plus(y,np,lp,n,l,bigL)
        implicit none
    
        INTEGER, INTENT(IN) :: np,lp,n,l, bigL
        REAL(kind=8),INTENT(IN) :: y
        INTEGER :: mp, m
        REAL(kind=8) :: BesselF3_plus
    
        BesselF3_plus = 0.d0
    
        do m = 0, n-1
            do mp = 0, np-1
                BesselF3_plus = BesselF3_plus + Summand1(m,mp,n,np) &
                    * Gammas_plusminus(m,mp,l,lp,bigL) * CHGF_plus(y,m,mp,l,lp,bigL)
            end do
        end do
    
    end function BesselF3_plus
    
    
    !-------------------------------------------------------------------------------
    function BesselElement(y,np,lp,n,l,bigL)
    
        implicit none
    
        INTEGER, INTENT(IN) :: np,lp,n,l, bigL
        REAL(kind=8),INTENT(IN) :: y
        REAL(kind=8) :: BesselElement
    
        BesselElement = BesselF1(y,np,lp,n,l,bigL) &
            * BesselF2(np,lp,n,l) * BesselF3(y,np,lp,n,l,bigL)
    
    end function BesselElement
    
    
    !-------------------------------------------------------------------------------
    function BesselElementminus(y,np,lp,n,l,bigL)
    
        implicit none
    
        INTEGER, INTENT(IN) :: np,lp,n,l, bigL
        REAL(kind=8),INTENT(IN) :: y
        REAL(kind=8) :: BesselElementminus
    
        BesselElementminus = BesselF1_plusminus(y,np,lp,n,l,bigL) &
            * BesselF2(np,lp,n,l) * BesselF3_minus(y,np,lp,n,l,bigL)
    
    end function BesselElementminus
    
    
    !-------------------------------------------------------------------------------
    function BesselElementplus(y,np,lp,n,l,bigL)
    
        implicit none
    
        INTEGER, INTENT(IN) :: np,lp,n,l, bigL
        REAL(kind=8),INTENT(IN) :: y
        REAL(kind=8) :: BesselElementplus
    
        BesselElementplus = BesselF1_plusminus(y,np,lp,n,l,bigL) &
            * BesselF2(np,lp,n,l) * BesselF3_plus(y,np,lp,n,l,bigL)
    
    end function BesselElementplus
    
    !-------------------------------------------------------------------------------
    function Summand1(m,mp,n,np)
        implicit none
        INTEGER, INTENT(IN) :: m, mp, n, np
        REAL(kind=8) :: Summand1
    
        Summand1 = (-1)**(m+mp) / exp(logfac(m)+logfac(mp) &
            + logfac(n-1-m)+logfac(np-1-mp))
    
    end function Summand1
    
    
    !-------------------------------------------------------------------------------
    function Gammas_regular(m,mp,l,lp,bigL)
        implicit none
        INTEGER, INTENT(IN) :: m, mp, l, lp,bigL
        REAL(kind=8) :: Gammas_regular

        Gammas_regular = GAMMA((l+lp+bigL+2*m+2*mp+3.d0)/2.d0) &
            / GAMMA(l+m+3.d0/2.d0) / GAMMA(lp+mp+3.d0/2.d0)

    end function Gammas_regular
    
    
    !-------------------------------------------------------------------------------
    function Gammas_plusminus(m,mp,l,lp,bigL)
        implicit none
        INTEGER, INTENT(IN) :: m, mp, l, lp, bigL
        REAL(kind=8) :: Gammas_plusminus

        Gammas_plusminus =  GAMMA((l+lp+bigL+2*m+2*mp+2.d0)/2.d0) &
            / GAMMA(l+m+3.d0/2.d0) / GAMMA(lp+mp+3.d0/2.d0)

    end function Gammas_plusminus
    
    
    !-------------------------------------------------------------------------------
    function CHGF_regular(y,m,mp,l,lp,bigL)
        implicit none
        INTEGER, INTENT(IN) :: m, mp, l, lp,bigL
        REAL(kind=8), INTENT(IN) :: y
        REAL(kind=8) :: a, b
        REAL(kind=8) :: CHGF_regular
    
        a = dble(bigL-lp-l-2*mp-2*m)/2.d0
        b = dble(bigL)+3.d0/2.d0

        !  CHGF_regular = chg(a, b, y) 
        call chgm(a,b,y,CHGF_regular)

    end function CHGF_regular
    
    
    !-------------------------------------------------------------------------------
    function CHGF_minus(y,m,mp,l,lp,bigL)
        implicit none
        INTEGER, INTENT(IN) :: m, mp, l, lp, bigL
        REAL(kind=8), INTENT(IN) :: y
        REAL(kind=8) :: a, b, ap,hg1,hg2
        REAL(kind=8) :: CHGF_minus
    
        a = (bigL-lp-l-2*mp-2*m-1.d0)/2.d0
        ap= (bigL-lp-l-2*mp-2*m+1.d0)/2.d0
        b = bigL+3.d0/2.d0
    
        call chgm(a,b,y,hg1)
        call chgm(ap,b,y,hg2)
    
        ! CHGF_minus = -(l+lp+bigL+2*m+2*mp+2.d0)/2.d0*chg(a, b, y)+2*m*chg(ap,b,y)
        CHGF_minus = -(l+lp+bigL+2*m+2*mp+2.d0)/2.d0*hg1+2*m*hg2
    
    end function CHGF_minus
    
    
    !-------------------------------------------------------------------------------
    function CHGF_plus(y,m,mp,l,lp,bigL)
        implicit none
    
        INTEGER, INTENT(IN) :: m, mp, l, lp, bigL
        REAL(kind=8), INTENT(IN) :: y
        REAL(kind=8) :: a, b, ap,hg1,hg2
        REAL(kind=8) :: CHGF_plus
    
        a = (bigL-lp-l-2*mp-2*m-1.d0)/2.d0
        ap= (bigL-lp-l-2*mp-2*m+1.d0)/2.d0
        b = bigL+3.d0/2.d0
    
        call chgm(a,b,y,hg1)
        call chgm(ap,b,y,hg2)
    
        ! CHGF_plus = -(l+lp+bigL+2*m+2*mp+2.d0)/2.d0*chg(a, b, y)+(2*l+2*m+1.d0)*chg(ap,b,y)
        CHGF_plus = -(l+lp+bigL+2*m+2*mp+2.d0)/2.d0*hg1+(2*l+2*m+1.d0)*hg2
    
    end function CHGF_plus

end module bessel
