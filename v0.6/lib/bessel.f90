!-------------------------------------------------------------------------------
function BesselF1(y,np,lp,n,l,bigL)
    implicit none
    INTERFACE

        function FacLOG(N)
            implicit none
            INTEGER, INTENT(IN) :: N
            REAL(kind = 8) ::  FacLOG
        end function FacLOG

        function DBLEFacLOG(N)
            implicit none
            INTEGER, INTENT(IN) :: N
            REAL(kind = 8) :: DBLEFacLOG
        end function DBLEFacLOG

    end INTERFACE

    INTEGER, INTENT(IN) :: np,lp,n,l,bigL
    REAL(kind=8),INTENT(IN) :: y
    REAL(kind=8) :: BesselF1

    BesselF1 = 2.0**(dble(bigL)) / exp(DBLEFacLOG(2*bigL+1)) * y**(bigL/2.d0) &
        & * exp(-y) * sqrt (EXP(FacLOG(np-1)+FacLOG(n-1)))

end function BesselF1


!-------------------------------------------------------------------------------
function BesselF1A(y,np,lp,n,l,bigL)
    implicit none
    INTERFACE

         function FacLOG(N)
             implicit none
             INTEGER, INTENT(IN) :: N
             REAL(kind=8) ::  FacLOG
         end function FacLOG

         function DBLEFacLOG(N)
             implicit none
             INTEGER, INTENT(IN) :: N
             REAL(kind=8) :: DBLEFacLOG
         end function DBLEFacLOG

    end INTERFACE

    INTEGER, INTENT(IN) :: np,lp,n,l,bigL
    REAL(kind=8),INTENT(IN) :: y
    REAL(kind=8) :: BesselF1A

    BesselF1A = 2**(dble(bigL)-1.0) / exp(DBLEFacLOG(2*bigL+1)) * y**((dble(bigL)-1.0)/2.d0) &
        & * exp(-y) * sqrt (EXP(FacLOG(np-1)+FacLOG(n-1)))

end function BesselF1A


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
    INTERFACE

        function Summand1(m,mp,n,np)
            implicit none
            INTEGER, INTENT(IN) :: m, mp, n, np
            Real(kind=8) ::  Summand1
        end function Summand1

        function Summand2(m,mp,l,lp,bigL)
            implicit none
            INTEGER, INTENT(IN) :: m, mp, l, lp, bigL
            Real(kind=8) :: Summand2
        end function Summand2

        function Summand3(y,m,mp,l,lp,bigL)
            implicit none
            INTEGER, INTENT(IN) :: m, mp, l, lp, bigL
            REAL(kind=8), INTENT(IN) :: y
            Real(kind=8) ::  Summand3
        end function Summand3

    end INTERFACE

    INTEGER, INTENT(IN) :: np,lp,n,l, bigL
    REAL(kind=8),INTENT(IN) :: y
    INTEGER :: m,mp!ip, i, m, mp
    REAL(kind=8) :: BesselF3

    BesselF3 = 0.0d0

    do m = 0, n-1
        do mp = 0, np-1
            BesselF3 = BesselF3 + Summand1(m,mp,n,np) * Summand2(m,mp,l,lp,bigL) * Summand3(y,m,mp,l,lp,bigL)
        end do
    end do

end function BesselF3


!-------------------------------------------------------------------------------
function BesselF3A(y,np,lp,n,l,bigL)
    implicit none
    INTERFACE

        function Summand1(m,mp,n,np)
            implicit none
            INTEGER, INTENT(IN) :: m, mp, n, np
            Real(kind=8) ::  Summand1
        end function Summand1

        function Summand2A(m,mp,l,lp,bigL)
            implicit none
            INTEGER, INTENT(IN) :: m, mp, l, lp, bigL
            Real(kind=8) ::  Summand2A
        end function Summand2A

        function Summand3A(y,m,mp,l,lp,bigL)
            implicit none
            INTEGER, INTENT(IN) :: m, mp, l, lp, bigL
            real(kind=8), INTENT(IN) :: y
            Real(kind=8) ::  Summand3A
        end function Summand3A

    end INTERFACE

    INTEGER, INTENT(IN) :: np,lp,n,l, bigL
    REAL(kind=8),INTENT(IN) :: y
    INTEGER :: mp, m
    REAL(kind=8) :: BesselF3A

    BesselF3A = 0.d0

    do m = 0, n-1
        do mp = 0, np-1
            BesselF3A = BesselF3A + Summand1(m,mp,n,np) * Summand2A(m,mp,l,lp,bigL) * Summand3A(y,m,mp,l,lp,bigL)
        end do
    end do

end function BesselF3A


!-------------------------------------------------------------------------------
function BesselF4A(y,np,lp,n,l,bigL)
    implicit none
    INTERFACE

        function Summand1(m,mp,n,np)
            implicit none
            INTEGER, INTENT(IN) :: m, mp, n, np
            Real(kind=8) ::  Summand1
        end function Summand1

        function Summand2A(m,mp,l,lp,bigL)
            implicit none
            INTEGER, INTENT(IN) :: m, mp, l, lp, bigL
            Real(kind=8) :: Summand2A
        end function Summand2A

        function Summand4A(y,m,mp,l,lp,bigL)
            implicit none
            INTEGER, INTENT(IN) :: m, mp, l, lp, bigL
            REAL(kind=8), INTENT(IN) :: y
            Real(kind=8) :: Summand4A
        end function Summand4A

    end INTERFACE

    INTEGER, INTENT(IN) :: np,lp,n,l, bigL
    REAL(kind=8),INTENT(IN) :: y
    INTEGER :: mp, m
    REAL(kind=8) :: BesselF4A

    BesselF4A = 0.d0

    do m = 0, n-1
        do mp = 0, np-1
            BesselF4A = BesselF4A + Summand1(m,mp,n,np) * Summand2A(m,mp,l,lp,bigL) * Summand4A(y,m,mp,l,lp,bigL)
        end do
    end do

end function BesselF4A


!-------------------------------------------------------------------------------
function BesselElement(y,np,lp,n,l,bigL)

    implicit none
    INTERFACE

        function BesselF1(y,np,lp,n,l,bigL)
            implicit none
            INTEGER, INTENT(IN) :: n, np, l, lp, bigL
            REAL(kind=8), INTENT(IN) :: y
            REAL(kind=8) :: BesselF1
        end function BesselF1

        function BesselF2(np,lp,n,l)
            implicit none
            INTEGER, INTENT(IN) :: np, lp, n, l
            Real(kind=8) ::  BesselF2
        end function BesselF2

        function BesselF3(y,np,lp,n,l,bigL)
            implicit none
            INTEGER, INTENT(IN) :: n, np, l, lp, bigL
            REAL(kind=8), INTENT(IN) :: y
            REAL(kind=8) :: BesselF3
        end function BesselF3

    end INTERFACE

    INTEGER, INTENT(IN) :: np,lp,n,l, bigL
    REAL(kind=8),INTENT(IN) :: y
    REAL(kind=8) :: BesselElement

    BesselElement = BesselF1(y,np,lp,n,l,bigL) * BesselF2(np,lp,n,l) * BesselF3(y,np,lp,n,l,bigL)

end function BesselElement


!-------------------------------------------------------------------------------
function BesselElementminus(y,np,lp,n,l,bigL)

    implicit none
    INTERFACE

        function BesselF1A(y,np,lp,n,l,bigL)
            implicit none
            INTEGER, INTENT(IN) :: n, np, l, lp, bigL
            REAL(kind=8), INTENT(IN) :: y
            REAL(kind=8) :: BesselF1A
        end function BesselF1A

        function BesselF2(np,lp,n,l)
            implicit none
            INTEGER, INTENT(IN) :: np, lp, n, l
            Real(kind=8) ::  BesselF2
        end function BesselF2

        function BesselF3A(y,np,lp,n,l,bigL)
            implicit none
            INTEGER, INTENT(IN) :: n, np, l, lp, bigL
            REAL(kind=8), INTENT(IN) :: y
            REAL(kind=8) :: BesselF3A
        end function BesselF3A

    end INTERFACE

    INTEGER, INTENT(IN) :: np,lp,n,l, bigL
    REAL(kind=8),INTENT(IN) :: y
    REAL(kind=8) :: BesselElementminus

    BesselElementminus = BesselF1A(y,np,lp,n,l,bigL) * BesselF2(np,lp,n,l) * BesselF3A(y,np,lp,n,l,bigL)

end function BesselElementminus


!-------------------------------------------------------------------------------
function BesselElementplus(y,np,lp,n,l,bigL)

    implicit none
    INTERFACE

        function BesselF1A(y,np,lp,n,l,bigL)
            implicit none
            INTEGER, INTENT(IN) :: n, np, l, lp, bigL
            REAL(kind=8), INTENT(IN) :: y
            REAL(kind=8) :: BesselF1A
        end function BesselF1A

        function BesselF2(np,lp,n,l)
            implicit none
            INTEGER, INTENT(IN) :: np, lp, n, l
            Real(kind=8) ::  BesselF2
        end function BesselF2

        function BesselF4A(y,np,lp,n,l,bigL)
            implicit none
            INTEGER, INTENT(IN) :: n, np, l, lp, bigL
            REAL(kind=8), INTENT(IN) :: y
            REAL(kind=8) :: BesselF4A
        end function BesselF4A

    end INTERFACE

    INTEGER, INTENT(IN) :: np,lp,n,l, bigL
    REAL(kind=8),INTENT(IN) :: y
    REAL(kind=8) :: BesselElementplus

    BesselElementplus = BesselF1A(y,np,lp,n,l,bigL) * BesselF2(np,lp,n,l) * BesselF4A(y,np,lp,n,l,bigL)

end function BesselElementplus

