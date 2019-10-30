!!-----------------------------------------------------------------------------------
!! ResponseNuclear(y,i,tau1,tau2) outputs the nuclear response functions, specified 
!! 0<=i<=8 (which correspond to the eight nuclear response functions and interferences 
!! in Eqn. 39) and the isospins tau1 and tau2.
!!-----------------------------------------------------------------------------------
!Program ResponseNuclear !(y, i, tau1, tau2)
!
!  IMPLICIT NONE
!
!  INTEGER :: i
!  INTEGER :: np,lp,jp,n,l,j,bigJ
!  REAL (kind=8) :: y
!  REAL (kind=8) :: operatorME
!
!  operatorME = 0.d0
!
!  i = 1
!  np= 1
!  lp= 1
!  jp= 1
!  n = 1
!  l = 1
!  j = 1
!  bigJ = 1
!  y = 1.0
!
!  Call OperME(i,y,np,lp,jp,n,l,j,bigJ,operatorME)
!
!
!end 


!-----------------------------------------------------------------------------------
! OperatorsME(j, operatorME) outputs the reduced matrix element of jth operator 
! (which correspond to the six operators and interferences in Eqn. 39).
!-----------------------------------------------------------------------------------
subroutine OperME(i,y,np,lp,jp,n,l,j,bigJ,operatorME)

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: i
  INTEGER, INTENT(IN) :: np,lp,jp,n,l,j,bigJ
  REAL (kind=8), INTENT(IN)  :: y
  REAL (kind=8), INTENT(OUT) :: operatorME
  REAL (kind=8) :: MJ,SigmaJ,PhiPPoverall,Qnorm,PhiPPsummand1,PhiPPsummand2,PhiPPsummand3,PhiPPsummand4,&
                  &MJLDivQ,SigmaPJ,SigmaPPJ 


  If (i .LT. 1 .or. i .GT. 7) then 
     print *, 'The choice of operator j should be from 1 to 6'
     return

! i = 1, for operator MJ

  else if (i .eq. 1) then

  operatorME = MJ(y,np,lp,jp,n,l,j,bigJ)

! i = 2, for operator SigmaJ

  else if (i .eq. 2) then

  operatorME = SigmaJ(y,np,lp,jp,n,l,j,bigJ,bigJ)

! i = 3, for operator PhiPPJ

else if (i .eq. 3) then

  operatorMe = Qnorm(dble(bigJ)+1.0)*SQRT(dble(bigJ)+1.0)*(PhiPPsummand1(y,np,lp,jp,n,l,j,bigJ) &
             & +PhiPPsummand2(y,np,lp,jp,n,l,j,bigJ))
  if (bigJ .ne. 0) then
  operatorMe = operatorMe + Qnorm(dble(bigJ)-1)*SQRT(dble(bigJ)) &
             & *(PhiPPsummand3(y,np,lp,jp,n,l,j,bigJ)    &
             & +PhiPPsummand4(y,np,lp,jp,n,l,j,bigJ))
  end if

  operatorME = operatorME * PhiPPoverall(lp,jp,l,j)

! i = 4, for operator PhiTPJ

else if (i .eq. 4) then

  operatorME = -Qnorm(dble(bigJ)+1.0)*SQRT(dble(bigJ))*(PhiPPsummand1(y,np,lp,jp,n,l,j,bigJ) &
             & +PhiPPsummand2(y,np,lp,jp,n,l,j,bigJ))
  if (bigJ .ne. 0) then
  operatorMe = operatorMe + Qnorm(dble(bigJ)-1.0)*SQRT(dble(bigJ)+1.0) &
             & *(PhiPPsummand3(y,np,lp,jp,n,l,j,bigJ)    &
             & +PhiPPsummand4(y,np,lp,jp,n,l,j,bigJ))
  end if

  operatorME = operatorME * PhiPPoverall(lp,jp,l,j)

  operatorME = operatorME + SigmaJ(y,np,lp,jp,n,l,j,bigJ,bigJ)/2.0

! i = 5, for operator DeltaJ

else if (i .eq. 5) then

  operatorME = MJLDivQ(y,np,lp,jp,n,l,j,bigJ,bigJ)

! i = 6, for operator SigmaPJ

else if (i .eq. 6) then
 
  operatorME = SigmaPJ(y,np,lp,jp,n,l,j,bigJ) 

! i = 7, for operator SigmaPPJ

else if (i .eq. 7) then

  operatorME = SigmaPPJ(y,np,lp,jp,n,l,j,bigJ)

else
  

end if

end 


Function Jnorm(j)
  IMPLICIT NONE

  REAL(kind=8), INTENT(IN)  :: j
  REAL(kind=8) :: Jnorm

  Jnorm = 2.0 * j + 1.0

End Function Jnorm

Function Qnorm(j)

  IMPLICIT NONE

  REAL(kind=8), INTENT(IN) :: j
  Real(kind=8)  :: Qnorm

  Qnorm = sqrt(2.0 * j + 1.0)

End Function Qnorm

Function BesselF1(y,np,lp,n,l,bigL)

  IMPLICIT NONE

  INTERFACE

     FUNCTION FacLOG(N)

       IMPLICIT NONE
       INTEGER, INTENT(IN) :: N
       REAL(kind = 8) ::  FacLOG
     END FUNCTION FacLOG

     FUNCTION DBLEFacLOG(N)
  
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: N
       REAL(kind = 8) :: DBLEFacLOG
     END FUNCTION DBLEFacLOG

  END INTERFACE

  INTEGER, INTENT(IN) :: np,lp,n,l,bigL
  REAL(kind=8),INTENT(IN) :: y
  REAL(kind=8) :: BesselF1

  BesselF1 = 2.0**(dble(bigL)) / exp(DBLEFacLOG(2*bigL+1)) * y**(bigL/2.d0) &
           & * exp(-y) * sqrt (EXP(FacLOG(np-1)+FacLOG(n-1)))

END FUNCTION BesselF1

Function BesselF1A(y,np,lp,n,l,bigL)

  IMPLICIT NONE

  INTERFACE

     FUNCTION FacLOG(N)
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: N
       REAL(kind=8) ::  FacLOG
     END FUNCTION FacLOG

     FUNCTION DBLEFacLOG(N)
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: N
       REAL(kind=8) :: DBLEFacLOG
     END FUNCTION DBLEFacLOG

  END INTERFACE

  INTEGER, INTENT(IN) :: np,lp,n,l,bigL
  REAL(kind=8),INTENT(IN) :: y
  REAL(kind=8) :: BesselF1A

  BesselF1A = 2**(dble(bigL)-1.0) / exp(DBLEFacLOG(2*bigL+1)) * y**((dble(bigL)-1.0)/2.d0) &
           & * exp(-y) * sqrt (EXP(FacLOG(np-1)+FacLOG(n-1)))

END FUNCTION BesselF1A


Function BesselF2(np,lp,n,l)

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: np,lp,n,l
  REAL(kind=8) :: BesselF2

  BesselF2 = sqrt(gamma(np+lp+0.5d0) * gamma(n+l+0.5d0))

END FUNCTION BesselF2

Function Summand1(m,mp,n,np)

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: m, mp, n, np
  REAL(kind=8) :: Summand1

  INTERFACE

     FUNCTION FacLOG(N)
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: N
       REAL(kind=8) ::  FacLOG
     END FUNCTION FacLOG

  End INTERFACE

  Summand1 = (-1)**(m+mp) / exp(FacLOG(m)+FacLOG(mp)+FacLOG(n-1-m)+FacLOG(np-1-mp))

End Function Summand1

Function Summand2(m,mp,l,lp,bigL)

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: m, mp, l, lp,bigL

  REAL(kind=8) :: Summand2

  Summand2 = GAMMA((l+lp+bigL+2*m+2*mp+3.d0)/2.d0) / GAMMA(l+m+3.d0/2.d0) / GAMMA(lp+mp+3.d0/2.d0)

End Function Summand2

Function Summand2A(m,mp,l,lp,bigL)

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: m, mp, l, lp, bigL

  REAL(kind=8) :: Summand2A

  Summand2A =  GAMMA((l+lp+bigL+2*m+2*mp+2.d0)/2.d0) / GAMMA(l+m+3.d0/2.d0) / GAMMA(lp+mp+3.d0/2.d0)

END FUNCTION Summand2A

Function Summand3(y,m,mp,l,lp,bigL)

  IMPLICIT NONE

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
End Function Summand3

Function Summand3A(y,m,mp,l,lp,bigL)

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: m, mp, l, lp, bigL
  REAL(kind=8), INTENT(IN) :: y
  REAL(kind=8) :: a, b, ap,hg1,hg2
  REAL(kind=8) :: Summand3A

  a = (bigL-lp-l-2*mp-2*m-1.d0)/2.d0
  ap= (bigL-lp-l-2*mp-2*m+1.d0)/2.d0
  b = bigL+3.d0/2.d0

  call chgm(a,b,y,hg1)
  call chgm(ap,b,y,hg2)

!  Summand3A = -(l+lp+bigL+2*m+2*mp+2.d0)/2.d0*chg(a, b, y)+2*m*chg(ap,b,y)
  Summand3A = -(l+lp+bigL+2*m+2*mp+2.d0)/2.d0*hg1+2*m*hg2

End Function Summand3A

Function Summand4A(y,m,mp,l,lp,bigL)

  IMPLICIT NONE

!  INTERFACE
!
!     FUNCTION chg(a,b,x)
!       
!       IMPLICIT NONE
!       REAL(kind=8) :: a, b
!       REAL(kind=8) :: x
!       REAL(kind=8) ::  chg
!     END FUNCTION chg
!
!  End INTERFACE

  INTEGER, INTENT(IN) :: m, mp, l, lp, bigL
  REAL(kind=8), INTENT(IN) :: y
  REAL(kind=8) :: a, b, ap,hg1,hg2
  REAL(kind=8) :: Summand4A

  a = (bigL-lp-l-2*mp-2*m-1.d0)/2.d0
  ap= (bigL-lp-l-2*mp-2*m+1.d0)/2.d0
  b = bigL+3.d0/2.d0

  call chgm(a,b,y,hg1)
  call chgm(ap,b,y,hg2)

!  Summand4A = -(l+lp+bigL+2*m+2*mp+2.d0)/2.d0*chg(a, b, y)+(2*l+2*m+1.d0)*chg(ap,b,y)
  Summand4A = -(l+lp+bigL+2*m+2*mp+2.d0)/2.d0*hg1+(2*l+2*m+1.d0)*hg2

End Function Summand4A


FUNCTION BesselF3(y,np,lp,n,l,bigL)

  IMPLICIT NONE

  INTERFACE

     FUNCTION Summand1(m,mp,n,np)

       IMPLICIT NONE
       INTEGER, INTENT(IN) :: m, mp, n, np
       Real(kind=8) ::  Summand1
     END FUNCTION Summand1

     FUNCTION Summand2(m,mp,l,lp,bigL)

       IMPLICIT NONE
       INTEGER, INTENT(IN) :: m, mp, l, lp, bigL
       Real(kind=8) :: Summand2
     END FUNCTION Summand2

     FUNCTION Summand3(y,m,mp,l,lp,bigL)

       IMPLICIT NONE
       INTEGER, INTENT(IN) :: m, mp, l, lp, bigL
       REAL(kind=8), INTENT(IN) :: y
       Real(kind=8) ::  Summand3
     END FUNCTION Summand3

  End INTERFACE

  INTEGER, INTENT(IN) :: np,lp,n,l, bigL
  REAL(kind=8),INTENT(IN) :: y
  INTEGER :: ip, i, m, mp
  REAL(kind=8) :: BesselF3

  BesselF3 = 0.0d0

  Do m = 0, n-1
    Do mp = 0, np-1
      BesselF3 = BesselF3 + Summand1(m,mp,n,np) * Summand2(m,mp,l,lp,bigL) * Summand3(y,m,mp,l,lp,bigL)
    end do
  end do

end FUNCTION BesselF3

FUNCTION BesselF3A(y,np,lp,n,l,bigL)

  IMPLICIT NONE

  INTERFACE

     FUNCTION Summand1(m,mp,n,np)

       IMPLICIT NONE
       INTEGER, INTENT(IN) :: m, mp, n, np
       Real(kind=8) ::  Summand1
     END FUNCTION Summand1

     FUNCTION Summand2A(m,mp,l,lp,bigL)

       IMPLICIT NONE
       INTEGER, INTENT(IN) :: m, mp, l, lp, bigL
       Real(kind=8) ::  Summand2A
     END FUNCTION Summand2A

     FUNCTION Summand3A(y,m,mp,l,lp,bigL)

       IMPLICIT NONE
       INTEGER, INTENT(IN) :: m, mp, l, lp, bigL
       real(kind=8), INTENT(IN) :: y
       Real(kind=8) ::  Summand3A
     END FUNCTION Summand3A

  End INTERFACE

  INTEGER, INTENT(IN) :: np,lp,n,l, bigL
  REAL(kind=8),INTENT(IN) :: y
  INTEGER :: mp, m
  REAL(kind=8) :: BesselF3A

  BesselF3A = 0.d0

  Do m = 0, n-1
    Do mp = 0, np-1
      BesselF3A = BesselF3A + Summand1(m,mp,n,np) * Summand2A(m,mp,l,lp,bigL) * Summand3A(y,m,mp,l,lp,bigL)
    end do
  end do

end FUNCTION BesselF3A

FUNCTION BesselF4A(y,np,lp,n,l,bigL)

  IMPLICIT NONE

  INTERFACE

     FUNCTION Summand1(m,mp,n,np)

       IMPLICIT NONE
       INTEGER, INTENT(IN) :: m, mp, n, np
       Real(kind=8) ::  Summand1
     END FUNCTION Summand1

     FUNCTION Summand2A(m,mp,l,lp,bigL)

       IMPLICIT NONE
       INTEGER, INTENT(IN) :: m, mp, l, lp, bigL
       Real(kind=8) :: Summand2A
     END FUNCTION Summand2A

     FUNCTION Summand4A(y,m,mp,l,lp,bigL)

       IMPLICIT NONE
       INTEGER, INTENT(IN) :: m, mp, l, lp, bigL
       REAL(kind=8), INTENT(IN) :: y
       Real(kind=8) :: Summand4A
     END FUNCTION Summand4A

  End INTERFACE

  INTEGER, INTENT(IN) :: np,lp,n,l, bigL
  REAL(kind=8),INTENT(IN) :: y
  INTEGER :: mp, m
  REAL(kind=8) :: BesselF4A

  BesselF4A = 0.d0

  Do m = 0, n-1
    Do mp = 0, np-1
      BesselF4A = BesselF4A + Summand1(m,mp,n,np) * Summand2A(m,mp,l,lp,bigL) * Summand4A(y,m,mp,l,lp,bigL)
    end do
  end do

end FUNCTION BesselF4A

Function BesselElement(y,np,lp,n,l,bigL)

  IMPLICIT NONE

  INTERFACE

     FUNCTION BesselF1(y,np,lp,n,l,bigL)

       IMPLICIT NONE
       INTEGER, INTENT(IN) :: n, np, l, lp, bigL
       REAL(kind=8), INTENT(IN) :: y
       REAL(kind=8) :: BesselF1
     END FUNCTION BesselF1

     FUNCTION BesselF2(np,lp,n,l)

       IMPLICIT NONE
       INTEGER, INTENT(IN) :: np, lp, n, l
       Real(kind=8) ::  BesselF2
     END FUNCTION BesselF2

     FUNCTION BesselF3(y,np,lp,n,l,bigL)

       IMPLICIT NONE
       INTEGER, INTENT(IN) :: n, np, l, lp, bigL
       REAL(kind=8), INTENT(IN) :: y
       REAL(kind=8) :: BesselF3
     END FUNCTION BesselF3

  End INTERFACE

  INTEGER, INTENT(IN) :: np,lp,n,l, bigL
  REAL(kind=8),INTENT(IN) :: y
  REAL(kind=8) :: BesselElement

  BesselElement = BesselF1(y,np,lp,n,l,bigL) * BesselF2(np,lp,n,l) * BesselF3(y,np,lp,n,l,bigL)

END FUNCTION BesselElement

Function BesselElementminus(y,np,lp,n,l,bigL)

  IMPLICIT NONE

  INTERFACE

     FUNCTION BesselF1A(y,np,lp,n,l,bigL)

       IMPLICIT NONE
       INTEGER, INTENT(IN) :: n, np, l, lp, bigL
       REAL(kind=8), INTENT(IN) :: y
       REAL(kind=8) :: BesselF1A
     END FUNCTION BesselF1A

     FUNCTION BesselF2(np,lp,n,l)

       IMPLICIT NONE
       INTEGER, INTENT(IN) :: np, lp, n, l
       Real(kind=8) ::  BesselF2
     END FUNCTION BesselF2
  
     FUNCTION BesselF3A(y,np,lp,n,l,bigL)
  
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: n, np, l, lp, bigL
       REAL(kind=8), INTENT(IN) :: y
       REAL(kind=8) :: BesselF3A
     END FUNCTION BesselF3A

  End INTERFACE

  INTEGER, INTENT(IN) :: np,lp,n,l, bigL
  REAL(kind=8),INTENT(IN) :: y
  REAL(kind=8) :: BesselElementminus

  BesselElementminus = BesselF1A(y,np,lp,n,l,bigL) * BesselF2(np,lp,n,l) * BesselF3A(y,np,lp,n,l,bigL)

END FUNCTION BesselElementminus

Function BesselElementplus(y,np,lp,n,l,bigL)

  IMPLICIT NONE

  INTERFACE

     FUNCTION BesselF1A(y,np,lp,n,l,bigL)

       IMPLICIT NONE
       INTEGER, INTENT(IN) :: n, np, l, lp, bigL
       REAL(kind=8), INTENT(IN) :: y
       REAL(kind=8) :: BesselF1A
     END FUNCTION BesselF1A

     FUNCTION BesselF2(np,lp,n,l)

       IMPLICIT NONE
       INTEGER, INTENT(IN) :: np, lp, n, l
       Real(kind=8) ::  BesselF2
     END FUNCTION BesselF2
  
     FUNCTION BesselF4A(y,np,lp,n,l,bigL)
  
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: n, np, l, lp, bigL
       REAL(kind=8), INTENT(IN) :: y
       REAL(kind=8) :: BesselF4A
     END FUNCTION BesselF4A

  End INTERFACE

  INTEGER, INTENT(IN) :: np,lp,n,l, bigL
  REAL(kind=8),INTENT(IN) :: y
  REAL(kind=8) :: BesselElementplus

  BesselElementplus = BesselF1A(y,np,lp,n,l,bigL) * BesselF2(np,lp,n,l) * BesselF4A(y,np,lp,n,l,bigL)

END FUNCTION BesselElementplus

Function MJ(y,np,lp,jp,n,l,j,bigJ)

  IMPLICIT NONE

  INTERFACE
    Function BesselElement(y,np,lp,n,l,bigL)

       implicit none

       INTEGER, INTENT(IN) :: n, np, l, lp, bigL
       REAL(kind=8), INTENT(IN) :: y
       REAL(kind=8) :: BesselElement
    end function BesselElement

    Function Jnorm(j)
      IMPLICIT NONE
    
      REAL(kind=8), INTENT(IN)  :: j
      REAL(kind=8) :: Jnorm
    end Function Jnorm

  End INTERFACE

  INTEGER, INTENT(IN) :: np,lp,jp,n,l,j,bigJ
  REAL(kind=8), INTENT(IN) :: y
  REAL(kind=8) :: xjp,xj
  REAL(kind=8) :: Pi = 3.1415926535897932
  REAL(kind=8) :: Wigner_3j,Wigner_6j
  REAL(kind=8) :: DBLEFacLOG
  REAL(kind=8) :: MJ


  xjp = dble(jp)/2.0
  xj  = dble(j )/2.0

  MJ = (-1)**(0.5d0+xj+dble(bigJ))*SQRT(Jnorm(xj)*Jnorm(xjp)*Jnorm(dble(l))*Jnorm(dble(lp))*Jnorm(dble(bigJ))/(4*Pi)) &
             & * Wigner_3j(2*lp,2*bigJ,2*l,0,0,0) * Wigner_6j(2*lp,jp,1,j,2*l,2*bigJ)   &
             & * BesselElement(y,np,lp,n,l,bigJ)

End function MJ

Function SigmaJ(y,np,lp,jp,n,l,j,bigJ,bigL)

  implicit none

  INTERFACE

    Function Jnorm(j)
      IMPLICIT NONE

      REAL(kind=8), INTENT(IN)  :: j
      REAL(kind=8) :: Jnorm
    end Function Jnorm

    Function BesselElement(y,np,lp,n,l,bigL)

       implicit none

       INTEGER, INTENT(IN) :: n, np, l, lp, bigL
       REAL(kind=8), INTENT(IN) :: y
       REAL(kind=8) :: BesselElement
    end function BesselElement

  End INTERFACE

  INTEGER, INTENT(IN) :: np,lp,jp,n,l,j
  INTEGER, INTENT(IN) :: bigJ,bigL
  REAL(kind=8), INTENT(IN) :: y
  REAL(kind=8) :: xjp,xj
  REAL(kind=8) :: Pi = 3.14159265358979323
  REAL(kind=8) :: Wigner_3j, Wigner_9j
  REAL(kind=8) :: SigmaJ

  xjp = dble(jp)/2.0
  xj  = dble(j )/2.0

  SigmaJ = (-1.0)**lp * SQRT(Jnorm(xj)*Jnorm(xjp)*Jnorm(dble(l))*Jnorm(dble(lp))  &
             & *Jnorm(dble(bigJ))*Jnorm(dble(bigL))/(4*Pi))                       &
             & * SQRT(6.0) * Wigner_3j(2*lp,2*bigL,2*l,0,0,0) * Wigner_9j(2*lp,2*l,2*bigL,1,1,2,jp,j,2*bigJ) &
             & * BesselElement(y,np,lp,n,l,bigL)


end function SigmaJ

Function PhiPPsummand1(y,np,lp,jp,n,l,j,bigJ)
  implicit none

  INTERFACE 

    Function Jnorm(j)
      IMPLICIT NONE

      REAL(kind=8), INTENT(IN)  :: j
      REAL(kind=8) :: Jnorm
    end Function Jnorm

    Function Qnorm(j)
      IMPLICIT NONE

      REAL(kind=8), INTENT(IN)  :: j
      REAL(kind=8) :: Qnorm
    end Function Qnorm

    Function BesselElementminus(y,np,lp,n,l,bigL)

       implicit none

       INTEGER, INTENT(IN) :: n, np, l, lp, bigL
       REAL(kind=8), INTENT(IN) :: y
       REAL(kind=8) :: BesselElementminus
    end function BesselElementminus

  End INTERFACE

  INTEGER, INTENT(IN) :: np,lp,jp,n,l,j,bigJ
  REAL(kind=8), INTENT(IN) :: y
  INTEGER :: bigL
  REAL(kind=8) :: xjp,xj
  REAL(kind=8) :: Wigner_3j,Wigner_6j,Wigner_9j
  REAL(kind=8) :: PhiPPsummand1

  xjp = dble(jp)/2.0
  xj  = dble(j )/2.0

  PhiPPsummand1 = 0.D0

  Do bigL = bigJ,bigJ+1
   PhiPPsummand1= PhiPPsummand1+(-1.0)**(bigJ-bigL+1.0)*Jnorm(dble(bigL))*Wigner_6j(2*(bigJ+1),2,2*bigL,2,2*bigJ,2) &
                & *Wigner_6j(2*(bigJ+1),2,2*bigL,2*l,2*lp,2*(l+1))                                         &
                & *Wigner_9j(2*lp,2*l,2*bigL,1,1,2,jp,j,2*bigJ)
  end do
   PhiPPsummand1= PhiPPsummand1*Qnorm(dble(l)+1.0)*SQRT(dble(l)+1.0)*Wigner_3j(2*lp,2*(bigJ+1),2*(l+1),0,0,0) &
                &*BesselElementminus(y,np,lp,n,l,bigJ+1)

end function PhiPPsummand1

Function PhiPPsummand2(y,np,lp,jp,n,l,j,bigJ)
  implicit none

  INTERFACE

    Function Jnorm(j)
      IMPLICIT NONE

      REAL(kind=8), INTENT(IN)  :: j
      REAL(kind=8) :: Jnorm
    end Function Jnorm

    Function Qnorm(j)
      IMPLICIT NONE

      REAL(kind=8), INTENT(IN)  :: j
      REAL(kind=8) :: Qnorm
    end Function Qnorm

    Function BesselElementplus(y,np,lp,n,l,bigL)

       implicit none

       INTEGER, INTENT(IN) :: n, np, l, lp, bigL
       REAL(kind=8),    INTENT(IN) :: y
       REAL(kind=8) :: BesselElementplus
    end function BesselElementplus

  End INTERFACE

  INTEGER, INTENT(IN) :: np,lp,jp,n,l,j
  INTEGER, INTENT(IN) :: bigJ
  INTEGER :: bigL
  REAL(kind=8), INTENT(IN) :: y
  REAL(kind=8) :: xjp,xj
  REAL(kind=8) :: Wigner_3j,Wigner_6j,Wigner_9j
  REAL(kind=8) :: PhiPPsummand2

  PhiPPsummand2 = 0.0

  If (l .ne. 0) then
  
  Do bigL = bigJ,bigJ+1
   PhiPPsummand2= PhiPPsummand2+(-1.0)**(bigJ-bigL)*Jnorm(dble(bigL))*Wigner_6j(2*(bigJ+1),2,2*bigL,2,2*bigJ,2) &
                & *Wigner_6j(2*(bigJ+1),2,2*bigL,2*l,2*lp,2*(l-1))                                         &
                & *Wigner_9j(2*lp,2*l,2*bigL,1,1,2,jp,j,2*bigJ)
  end do
   PhiPPsummand2= PhiPPsummand2*Qnorm(dble(l)-1.0)*SQRT(dble(l))*Wigner_3j(2*lp,2*(bigJ+1),2*(l-1),0,0,0) &
                &*BesselElementplus(y,np,lp,n,l,bigJ+1)
  end if

end function PhiPPsummand2

Function PhiPPsummand3(y,np,lp,jp,n,l,j,bigJ)
  implicit none

  INTERFACE

    Function Jnorm(j)
      IMPLICIT NONE

      REAL(kind=8), INTENT(IN)  :: j
      REAL(kind=8) :: Jnorm
    end Function Jnorm

    Function Qnorm(j)
      IMPLICIT NONE

      REAL(kind=8), INTENT(IN)  :: j
      REAL(kind=8) :: Qnorm
    end Function Qnorm 

    Function BesselElementplus(y,np,lp,n,l,bigL)

       implicit none

       INTEGER, INTENT(IN) :: n, np, l, lp, bigL
       REAL(kind=8), INTENT(IN) :: y
       REAL(kind=8) :: BesselElementplus
    end function BesselElementplus

  End INTERFACE

  INTEGER, INTENT(IN) :: np,lp,jp,n,l,j
  INTEGER, INTENT(IN) :: bigJ
  INTEGER :: bigL
  REAL(kind=8), INTENT(IN) :: y
  REAL(kind=8) :: Wigner_3j,Wigner_6j,Wigner_9j
  REAL(kind=8) :: PhiPPsummand3

  PhiPPsummand3 = 0.0

  If (bigJ .ne. 0) then

  Do bigL = bigJ-1,bigJ
   PhiPPsummand3= PhiPPsummand3+(-1.0)**(bigJ-bigL+1.0)*Jnorm(dble(bigL))*Wigner_6j(2*(bigJ-1),2,2*bigL,2,2*bigJ,2) &
                & *Wigner_6j(2*(bigJ-1),2,2*bigL,2*l,2*lp,2*(l-1))                                         &
                & *Wigner_9j(2*lp,2*l,2*bigL,1,1,2,jp,j,2*bigJ)
  end do
   PhiPPsummand3= PhiPPsummand3*Qnorm(dble(l)-1.0)*SQRT(dble(l))*Wigner_3j(2*lp,2*(bigJ-1),2*(l-1),0,0,0) &
                &*BesselElementplus(y,np,lp,n,l,bigJ-1)

  end if

end function PhiPPsummand3

Function PhiPPsummand4(y,np,lp,jp,n,l,j,bigJ)
  implicit none

  INTERFACE

    Function Jnorm(j)
      IMPLICIT NONE

      REAL(kind=8), INTENT(IN)  :: j
      REAL(kind=8) :: Jnorm
    end Function Jnorm

   Function Qnorm(j)
      IMPLICIT NONE

      REAL(kind=8), INTENT(IN)  :: j
      REAL(kind=8) :: Qnorm
    end Function Qnorm

    Function BesselElementminus(y,np,lp,n,l,bigL)

       implicit none

       INTEGER, INTENT(IN) :: n, np, l, lp, bigL
       REAL(kind=8),    INTENT(IN) :: y
       REAL(kind=8) :: BesselElementminus
    end function BesselElementminus

  End INTERFACE

  INTEGER, INTENT(IN) :: np,lp,jp,n,l,j,bigJ
  REAL(kind=8), INTENT(IN) :: y
  INTEGER :: bigL
  REAL(kind=8) :: Wigner_3j,Wigner_6j,Wigner_9j
  REAL(kind=8) :: PhiPPsummand4

  PhiPPsummand4 = 0.0

  if (bigJ .ne. 0 .and. l .ne. 0) then
  Do bigL = bigJ-1,bigJ
   PhiPPsummand4= PhiPPsummand4+(-1.0)**(bigJ-bigL)*Jnorm(dble(bigL))*Wigner_6j(2*(bigJ-1),2,2*bigL,2,2*bigJ,2) &
                & *Wigner_6j(2*(bigJ-1),2,2*bigL,2*l,2*lp,2*(l+1))                                         &
                & *Wigner_9j(2*lp,2*l,2*bigL,1,1,2,jp,j,2*bigJ)
  end do
   PhiPPsummand4= PhiPPsummand4*Qnorm(dble(l)+1.0)*SQRT(dble(l)+1.0)*Wigner_3j(2*lp,2*(bigJ-1),2*(l+1),0,0,0) &
                &*BesselElementminus(y,np,lp,n,l,bigJ-1)
  end if

end function PhiPPsummand4

Function PhiPPoverall(lp,jp,l,j)

  implicit none

  INTERFACE

  Function Jnorm(j)
  IMPLICIT NONE

  REAL(kind=8), INTENT(IN)  :: j
  REAL(kind=8) :: Jnorm
  END FUNCTION Jnorm

  END INTERFACE
  
  INTEGER, INTENT(IN) :: lp,jp,l,j
  REAL(kind=8) :: Pi = 3.1415926535897932
  REAL(kind=8) :: xjp,xj
  REAL(kind=8) :: PhiPPoverall

  xjp = dble(jp)/2.0
  xj  = dble(j )/2.0
  PhiPPoverall = (-1)**(lp+1) * 6.0 * SQRT(Jnorm(xj)*Jnorm(xjp)*Jnorm(dble(lp))/(4*Pi))

end function PhiPPoverall

Function MJLDivQoverall(lp,jp,l,j,bigJ,bigL)

  implicit none

  INTERFACE

  Function Qnorm(j)
  IMPLICIT NONE

  REAL(kind=8), INTENT(IN)  :: j
  REAL(kind=8) :: Qnorm
  END FUNCTION Qnorm

  END INTERFACE

  INTEGER, INTENT(IN) :: lp,jp,l,j,bigJ,bigL
  REAL(kind=8) :: xjp,xj
  REAL(kind=8) :: Pi = 3.1415926535897932
  REAL(kind=8) :: Wigner_6j
  REAL(kind=8) :: MJLDivQoverall

  xjp = dble(jp)/2.0
  xj  = dble(j )/2.0

  MJLDivQoverall = (-1)**(bigL+j+0.5)* Qnorm(dble(lp))*Qnorm(xjp)*Qnorm(xj)*Qnorm(dble(bigJ))*Qnorm(dble(bigL)) &
                 & *Wigner_6j(2*lp,jp,1,j,2*l,2*bigJ)/SQRT(4*Pi)

end function MJLDivQoverall

Function MJLDivQsummand1(y,np,lp,jp,n,l,j,bigJ,bigL)

  implicit none

  INTERFACE

    Function Qnorm(j)
      IMPLICIT NONE

      REAL(kind=8), INTENT(IN)  :: j
      REAL(kind=8) :: Qnorm
    end Function Qnorm

    Function BesselElementminus(y,np,lp,n,l,bigL)

       implicit none

       INTEGER, INTENT(IN) :: n, np, l, lp
       integer, INTENT(IN) :: bigL
       REAL(kind=8),    INTENT(IN) :: y
       REAL(kind=8) :: BesselElementminus
    end function BesselElementminus

  END INTERFACE

  INTEGER, INTENT(IN) :: np,lp,jp,n,l,j,bigJ,bigL
  REAL(kind=8), INTENT(IN) :: y 
  REAL(kind=8) :: Wigner_3j,Wigner_6j
  REAL(kind=8) :: MJLDivQsummand1

  MJLDivQsummand1 = -SQRT(dble(l)+1.0)*Qnorm(dble(l)+1.0)*Wigner_6j(2*bigL,2,2*bigJ,2*l,2*lp,2*(l+1))  &
                  & *Wigner_3j(2*lp,2*bigL,2*(l+1),0,0,0)*BesselElementminus(y,np,lp,n,l,bigL) 

END Function MJLDivQsummand1

Function MJLDivQsummand2(y,np,lp,jp,n,l,j,bigJ,bigL)

  implicit none

  INTERFACE

    Function Qnorm(j)
      IMPLICIT NONE

      REAL(kind=8), INTENT(IN)  :: j
      REAL(kind=8) :: Qnorm
    end Function Qnorm

    Function BesselElementplus(y,np,lp,n,l,bigL)

       implicit none

       INTEGER, INTENT(IN) :: n, np, l, lp
       integer, INTENT(IN) :: bigL
       REAL(kind=8),    INTENT(IN) :: y
       REAL(kind=8) :: BesselElementplus
    end function BesselElementplus

  END INTERFACE

  INTEGER, INTENT(IN) :: np,lp,jp,n,l,j,bigJ,bigL
  REAL(kind=8), INTENT(IN) :: y
  REAL(kind=8) :: Wigner_3j,Wigner_6j
  REAL(kind=8) :: MJLDivQsummand2

  MJLDivQsummand2 = SQRT(dble(l))*Qnorm(dble(l)-1.0)*Wigner_6j(2*bigL,2,2*bigJ,2*l,2*lp,2*(l-1))  &
                  & *Wigner_3j(2*lp,2*bigL,2*(l-1),0,0,0)*BesselElementplus(y,np,lp,n,l,bigL)

END Function MJLDivQsummand2

Function MJLDivQ(y,np,lp,jp,n,l,j,bigJ,bigL)

  implicit none

   INTERFACE

    Function MJLDivQoverall(lp,jp,l,j,bigJ,bigL)
       implicit none

       INTEGER, INTENT(IN) :: lp,jp,l,j,bigJ,bigL
       REAL(kind=8) :: Pi = 3.1415926535897932
       REAL(kind=8) :: MJLDivQoverall

    end Function MJLDivQoverall

    Function MJLDivQsummand1(y,np,lp,jp,n,l,j,bigJ,bigL)

       implicit none

       INTEGER, INTENT(IN) :: np,lp,jp,n,l,j,bigJ,bigL
       REAL(kind=8), INTENT(IN) :: y
       REAL(kind=8) :: MJLDivQsummand1

    end function MJLDivQsummand1

    Function MJLDivQsummand2(y,np,lp,jp,n,l,j,bigJ,bigL)

       implicit none

       INTEGER, INTENT(IN) :: np,lp,jp,n,l,j,bigJ,bigL
       REAL(kind=8), INTENT(IN) :: y
       REAL(kind=8) :: MJLDivQsummand2

    end function MJLDivQsummand2

  END INTERFACE

       INTEGER, INTENT(IN) :: np,lp,jp,n,l,j,bigJ,bigL
       REAL(kind=8), INTENT(IN) :: y
       REAL(kind=8) :: MJLDivQ

       MJLDivQ = MJLDivQoverall(lp,jp,l,j,bigJ,bigL) * (MJLDivQsummand1(y,np,lp,jp,n,l,j,bigJ,bigL) &
               &+MJLDivQsummand2(y,np,lp,jp,n,l,j,bigJ,bigL))

end Function MJLDivQ

Function SigmaPJ(y,np,lp,jp,n,l,j,bigJ)

  INTERFACE

    Function Qnorm(j)
      IMPLICIT NONE

      REAL(kind=8), INTENT(IN)  :: j
      REAL(kind=8) :: Qnorm
    end Function Qnorm

    Function SigmaJ(y,np,lp,jp,n,l,j,bigJ,bigL)

      implicit none
      INTEGER, INTENT(IN) :: np,lp,jp,n,l,j
      INTEGER, INTENT(IN) :: bigJ,bigL
      REAL(kind=8), INTENT(IN) :: y
      REAL(kind=8) :: SigmaJ
   END Function SigmaJ

  END INTERFACE

      INTEGER, INTENT(IN) :: np,lp,jp,n,l,j
      INTEGER, INTENT(IN) :: bigJ
      REAL(kind=8), INTENT(IN) :: y
      REAL(kind=8) :: SigmaPJ

  SigmaPJ = -SQRT(dble(bigJ))/Qnorm(dble(bigJ))*SigmaJ(y,np,lp,jp,n,l,j,bigJ,bigJ+1)   &
          & +SQRT(dble(bigJ)+1.0)/Qnorm(dble(bigJ))*SigmaJ(y,np,lp,jp,n,l,j,bigJ,bigJ-1)

end Function SigmaPJ

Function SigmaPPJ(y,np,lp,jp,n,l,j,bigJ)

  INTERFACE

    Function Qnorm(j)
      IMPLICIT NONE

      REAL(kind=8), INTENT(IN)  :: j
      REAL(kind=8) :: Qnorm
    end Function Qnorm

    Function SigmaJ(y,np,lp,jp,n,l,j,bigJ,bigL)

      implicit none
      INTEGER, INTENT(IN) :: np,lp,jp,n,l,j
      INTEGER, INTENT(IN) :: bigJ,bigL
      REAL(kind=8), INTENT(IN) :: y
      REAL(kind=8) :: SigmaJ
   END Function SigmaJ

  END INTERFACE

      INTEGER, INTENT(IN) :: np,lp,jp,n,l,j
      INTEGER, INTENT(IN) :: bigJ
      REAL(kind=8), INTENT(IN) :: y
      REAL(kind=8) :: SigmaPPJ

  SigmaPPJ = SQRT(dble(bigJ)+1.0)/Qnorm(dble(bigJ))*SigmaJ(y,np,lp,jp,n,l,j,bigJ,bigJ+1)   &
          & +SQRT(dble(bigJ))/Qnorm(dble(bigJ))*SigmaJ(y,np,lp,jp,n,l,j,bigJ,bigJ-1)

end Function SigmaPPJ


!FUNCTION Triangle(a,b,c)
!  USE PRECISION
!  IMPLICIT NONE
!
!  ! \Delta(a/2, b/2, c/2)
!  ! \Delta(a/2, b/2, c/2) = sqrt( (a/2+b/2-c/2)! (a/2-b/2+c/2)! (-a/2+b/2+c/2)! / (a/2+b/2+c/2+1)! )
!  ! if a/2,b/2,c/2 do not satisfy the triangle inequality, it returns zero.
!
!  INTERFACE
!     FUNCTION FacLOG(N)
!       USE PRECISION
!       IMPLICIT NONE
!       INTEGER, INTENT(IN) :: N
!       REAL(DP) ::  FacLOG
!     END FUNCTION FacLOG
!  END INTERFACE
!
!  INTEGER, INTENT(IN) :: a,b,c
!
!  REAL(DP) :: Triangle, temp
!
!  REAL(DP) :: F1, F2, F3, F4
!
!  Triangle = 0.0D0
!  IF (a+b-c .LT. 0 .OR. a-b+c .LT. 0 .OR. -a+b+c .LT. 0 .OR. a+b+c+1 .LT. 0) RETURN
!
!  IF( MOD( a+b-c,2) .NE. 0) RETURN
!  IF( MOD( a-b+c,2) .NE. 0) RETURN
!  IF( MOD(-a+b+c,2) .NE. 0) RETURN
!  IF( MOD( a+b+c,2) .NE. 0) RETURN
!
!  F1 = FacLOG(( a + b - c    )/2)
!  F2 = FacLOG(( a - b + c    )/2)
!  F3 = FacLOG((-a + b + c    )/2)
!  F4 = FacLOG(( a + b + c)/2 + 1)
!
!  temp = F1 + F2 + F3 - F4
!  Triangle = EXP(0.50D0 * temp)
!
!!  Triangle = F1 * F2 * F3 / F4
!
!  RETURN
!
!END FUNCTION Triangle

subroutine chgm ( a, b, x, hg )

!*****************************************************************************80
!
!! CHGM computes the confluent hypergeometric function M(a,b,x).
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    27 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, parameters.
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) HG, the value of M(a,b,x).
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) a0
  real ( kind = 8 ) a1
  real ( kind = 8 ) aa
  real ( kind = 8 ) b
  real ( kind = 8 ) hg
  real ( kind = 8 ) hg1
  real ( kind = 8 ) hg2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) la
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nl
  real ( kind = 8 ) pi
  real ( kind = 8 ) r
  real ( kind = 8 ) r1
  real ( kind = 8 ) r2
  real ( kind = 8 ) rg
  real ( kind = 8 ) sum1
  real ( kind = 8 ) sum2
  real ( kind = 8 ) ta
  real ( kind = 8 ) tb
  real ( kind = 8 ) tba
  real ( kind = 8 ) x
  real ( kind = 8 ) x0
  real ( kind = 8 ) xg
  real ( kind = 8 ) y0
  real ( kind = 8 ) y1

  pi = 3.141592653589793D+00
  a0 = a
  a1 = a
  x0 = x
  hg = 0.0D+00

  if ( b == 0.0D+00 .or. b == - abs ( int ( b ) ) ) then
    hg = 1.0D+300
  else if ( a == 0.0D+00 .or. x == 0.0D+00 ) then
    hg = 1.0D+00
  else if ( a == -1.0D+00 ) then
    hg = 1.0D+00 - x / b
  else if ( a == b ) then
    hg = exp ( x )
  else if ( a - b == 1.0D+00 ) then
    hg = ( 1.0D+00 + x / b ) * exp ( x )
  else if ( a == 1.0D+00 .and. b == 2.0D+00 ) then
    hg = ( exp ( x ) - 1.0D+00 ) / x
  else if ( a == int ( a ) .and. a < 0.0D+00 ) then
    m = int ( - a )
    r = 1.0D+00
    hg = 1.0D+00
    do k = 1, m
      r = r * ( a + k - 1.0D+00 ) / k / ( b + k - 1.0D+00 ) * x
      hg = hg + r
    end do
  end if

  if ( hg /= 0.0D+00 ) then
    return
  end if

  if ( x < 0.0D+00 ) then
    a = b - a
    a0 = a
    x = abs ( x )
  end if

  if ( a < 2.0D+00 ) then
    nl = 0
  end if

  if ( 2.0D+00 <= a ) then
    nl = 1
    la = int ( a )
    a = a - la - 1.0D+00
  end if

  do n = 0, nl

    if ( 2.0D+00 <= a0 ) then
      a = a + 1.0D+00
    end if

    if ( x <= 30.0D+00 + abs ( b ) .or. a < 0.0D+00 ) then

      hg = 1.0D+00
      rg = 1.0D+00
      do j = 1, 500
        rg = rg * ( a + j - 1.0D+00 ) &
          / ( j * ( b + j - 1.0D+00 ) ) * x
        hg = hg + rg
        if ( abs ( rg / hg ) < 1.0D-15 ) then
          exit
        end if
      end do

    else

      call gamma ( a, ta )
      call gamma ( b, tb )
      xg = b - a
      call gamma ( xg, tba )
      sum1 = 1.0D+00
      sum2 = 1.0D+00
      r1 = 1.0D+00
      r2 = 1.0D+00
      do i = 1, 8
        r1 = - r1 * ( a + i - 1.0D+00 ) * ( a - b + i ) / ( x * i )
        r2 = - r2 * ( b - a + i - 1.0D+00 ) * ( a - i ) / ( x * i )
        sum1 = sum1 + r1
        sum2 = sum2 + r2
      end do
      hg1 = tb / tba * x ** ( - a ) * cos ( pi * a ) * sum1
      hg2 = tb / ta * exp ( x ) * x ** ( a - b ) * sum2
      hg = hg1 + hg2

    end if

    if ( n == 0 ) then
      y0 = hg
    else if ( n == 1 ) then
      y1 = hg
    end if

  end do

  if ( 2.0D+00 <= a0 ) then
    do i = 1, la - 1
      hg = ( ( 2.0D+00 * a - b + x ) * y1 + ( b - a ) * y0 ) / a
      y0 = y1
      y1 = hg
      a = a + 1.0D+00
    end do
  end if

  if ( x0 < 0.0D+00 ) then
    hg = hg * exp ( x0 )
  end if

  a = a1
  x = x0

  return
end

subroutine gamma ( x, ga )

!*****************************************************************************80
!
!! GAMMA evaluates the Gamma function.
!
!  Licensing:
!
!    The original FORTRAN77 version of this routine is copyrighted by 
!    Shanjie Zhang and Jianming Jin.  However, they give permission to 
!    incorporate this routine into a user program that the copyright 
!    is acknowledged.
!
!  Modified:
!
!    08 September 2007
!
!  Author:
!
!    Original FORTRAN77 version by Shanjie Zhang, Jianming Jin.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!    X must not be 0, or any negative integer.
!
!    Output, real ( kind = 8 ) GA, the value of the Gamma function.
!
  implicit none

  real ( kind = 8 ), dimension ( 26 ) :: g = (/ &
    1.0D+00, &
    0.5772156649015329D+00, &
   -0.6558780715202538D+00, &
   -0.420026350340952D-01, &
    0.1665386113822915D+00, &
   -0.421977345555443D-01, &
   -0.96219715278770D-02, &
    0.72189432466630D-02, &
   -0.11651675918591D-02, &
   -0.2152416741149D-03, &
    0.1280502823882D-03, & 
   -0.201348547807D-04, &
   -0.12504934821D-05, &
    0.11330272320D-05, &
   -0.2056338417D-06, & 
    0.61160950D-08, &
    0.50020075D-08, &
   -0.11812746D-08, &
    0.1043427D-09, & 
    0.77823D-11, &
   -0.36968D-11, &
    0.51D-12, &
   -0.206D-13, &
   -0.54D-14, &
    0.14D-14, &
    0.1D-15 /)
  real ( kind = 8 ) ga
  real ( kind = 8 ) gr
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m1
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r
  real ( kind = 8 ) x
  real ( kind = 8 ) z

  if ( x == aint ( x ) ) then

    if ( 0.0D+00 < x ) then
      ga = 1.0D+00
      m1 = int ( x ) - 1
      do k = 2, m1
        ga = ga * k
      end do
    else
      ga = 1.0D+300
    end if

  else

    if ( 1.0D+00 < abs ( x ) ) then
      z = abs ( x )
      m = int ( z )
      r = 1.0D+00
      do k = 1, m
        r = r * ( z - real ( k, kind = 8 ) )
      end do
      z = z - real ( m, kind = 8 )
    else
      z = x
    end if

    gr = g(26)
    do k = 25, 1, -1
      gr = gr * z + g(k)
    end do

    ga = 1.0D+00 / ( gr * z )

    if ( 1.0D+00 < abs ( x ) ) then
      ga = ga * r
      if ( x < 0.0D+00 ) then
        ga = - pi / ( x* ga * sin ( pi * x ) )
      end if
    end if

  end if

  return
end
