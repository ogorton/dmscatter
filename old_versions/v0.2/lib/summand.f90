!-------------------------------------------------------------------------------
function Summand1(m,mp,n,np)
    implicit none
    INTEGER, INTENT(IN) :: m, mp, n, np
    REAL(kind=8) :: Summand1
    INTERFACE
        function FacLOG(N)
            implicit none
            INTEGER, INTENT(IN) :: N
            REAL(kind=8) ::  FacLOG
        end function FacLOG
    end INTERFACE

    Summand1 = (-1)**(m+mp) / exp(FacLOG(m)+FacLOG(mp)+FacLOG(n-1-m)+FacLOG(np-1-mp))

end function Summand1


!-------------------------------------------------------------------------------
function Summand2(m,mp,l,lp,bigL)
    implicit none
    INTEGER, INTENT(IN) :: m, mp, l, lp,bigL
    REAL(kind=8) :: Summand2
    Summand2 = GAMMA((l+lp+bigL+2*m+2*mp+3.d0)/2.d0) / GAMMA(l+m+3.d0/2.d0) / GAMMA(lp+mp+3.d0/2.d0)
end function Summand2


!-------------------------------------------------------------------------------
function Summand2A(m,mp,l,lp,bigL)
    implicit none
    INTEGER, INTENT(IN) :: m, mp, l, lp, bigL
    REAL(kind=8) :: Summand2A
    Summand2A =  GAMMA((l+lp+bigL+2*m+2*mp+2.d0)/2.d0) / GAMMA(l+m+3.d0/2.d0) / GAMMA(lp+mp+3.d0/2.d0)
end function Summand2A


!-------------------------------------------------------------------------------
function Summand3(y,m,mp,l,lp,bigL)
    implicit none
    INTEGER, INTENT(IN) :: m, mp, l, lp,bigL
    REAL(kind=8), INTENT(IN) :: y
    REAL(kind=8) :: a, b
    REAL(kind=8) :: Summand3,hg

    a = dble(bigL-lp-l-2*mp-2*m)/2.d0
    b = dble(bigL)+3.d0/2.d0
    !  Summand3 = chg(a, b, y) 
    call chgm(a,b,y,Summand3)
    call chgm(0.0d0,1.5d0,4.0d0,hg)
    !  if (abs(a) .le. 1.0d-7) Summand3 = -chg(a,b,y)
    !  print*, 'here ', chg(0.0d0,1.5d0,4.d0), chg(-1.d0,1.5d0,4.d0), hg
end function Summand3


!-------------------------------------------------------------------------------
function Summand3A(y,m,mp,l,lp,bigL)
    implicit none
    INTEGER, INTENT(IN) :: m, mp, l, lp, bigL
    REAL(kind=8), INTENT(IN) :: y
    REAL(kind=8) :: a, b, ap,hg1,hg2
    REAL(kind=8) :: Summand3A

    a = (bigL-lp-l-2*mp-2*m-1.d0)/2.d0
    ap= (bigL-lp-l-2*mp-2*m+1.d0)/2.d0
    b = bigL+3.d0/2.d0

    call chgm(a,b,y,hg1)
    call chgm(ap,b,y,hg2)

    ! Summand3A = -(l+lp+bigL+2*m+2*mp+2.d0)/2.d0*chg(a, b, y)+2*m*chg(ap,b,y)
    Summand3A = -(l+lp+bigL+2*m+2*mp+2.d0)/2.d0*hg1+2*m*hg2

end function Summand3A


!-------------------------------------------------------------------------------
function Summand4A(y,m,mp,l,lp,bigL)
    implicit none
!    INTERFACE
!  
!       function chg(a,b,x)
!         
!         implicit none
!         REAL(kind=8) :: a, b
!         REAL(kind=8) :: x
!         REAL(kind=8) ::  chg
!       end function chg
!  
!    end INTERFACE

    INTEGER, INTENT(IN) :: m, mp, l, lp, bigL
    REAL(kind=8), INTENT(IN) :: y
    REAL(kind=8) :: a, b, ap,hg1,hg2
    REAL(kind=8) :: Summand4A

    a = (bigL-lp-l-2*mp-2*m-1.d0)/2.d0
    ap= (bigL-lp-l-2*mp-2*m+1.d0)/2.d0
    b = bigL+3.d0/2.d0

    call chgm(a,b,y,hg1)
    call chgm(ap,b,y,hg2)

    ! Summand4A = -(l+lp+bigL+2*m+2*mp+2.d0)/2.d0*chg(a, b, y)+(2*l+2*m+1.d0)*chg(ap,b,y)
    Summand4A = -(l+lp+bigL+2*m+2*mp+2.d0)/2.d0*hg1+(2*l+2*m+1.d0)*hg2

end function Summand4A


