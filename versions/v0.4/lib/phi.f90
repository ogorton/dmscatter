!-------------------------------------------------------------------------------
function PhiPPsummand1(y,np,lp,jp,n,l,j,bigJ)

    implicit none
    INTERFACE

        function Jnorm(j)
            implicit none
            REAL(kind=8), INTENT(IN)  :: j
            REAL(kind=8) :: Jnorm
        end function Jnorm

        function Qnorm(j)
            implicit none
            REAL(kind=8), INTENT(IN)  :: j
            REAL(kind=8) :: Qnorm
        end function Qnorm

        function BesselElementminus(y,np,lp,n,l,bigL)
            implicit none
            INTEGER, INTENT(IN) :: n, np, l, lp, bigL
            REAL(kind=8), INTENT(IN) :: y
            REAL(kind=8) :: BesselElementminus
        end function BesselElementminus

    end INTERFACE

    INTEGER, INTENT(IN) :: np,lp,jp,n,l,j,bigJ
    REAL(kind=8), INTENT(IN) :: y
    INTEGER :: bigL
    REAL(kind=8) :: xjp,xj
    REAL(kind=8) :: Wigner_3j,Wigner_6j,Wigner_9j
    REAL(kind=8) :: PhiPPsummand1

    xjp = dble(jp)/2.0
    xj  = dble(j )/2.0

    PhiPPsummand1 = 0.D0

    do bigL = bigJ,bigJ+1
        PhiPPsummand1= PhiPPsummand1+(-1.0)**(bigJ-bigL+1.0)*Jnorm(dble(bigL))*Wigner_6j(2*(bigJ+1),2,2*bigL,2,2*bigJ,2) &
              & *Wigner_6j(2*(bigJ+1),2,2*bigL,2*l,2*lp,2*(l+1))                                         &
              & *Wigner_9j(2*lp,2*l,2*bigL,1,1,2,jp,j,2*bigJ)
    end do
        PhiPPsummand1= PhiPPsummand1*Qnorm(dble(l)+1.0)*SQRT(dble(l)+1.0)*Wigner_3j(2*lp,2*(bigJ+1),2*(l+1),0,0,0) &
            &*BesselElementminus(y,np,lp,n,l,bigJ+1)

end function PhiPPsummand1


!-------------------------------------------------------------------------------
function PhiPPsummand2(y,np,lp,jp,n,l,j,bigJ)

    implicit none
    INTERFACE

        function Jnorm(j)
            implicit none
            REAL(kind=8), INTENT(IN)  :: j
            REAL(kind=8) :: Jnorm
        end function Jnorm

        function Qnorm(j)
            implicit none
            REAL(kind=8), INTENT(IN)  :: j
            REAL(kind=8) :: Qnorm
        end function Qnorm

        function BesselElementplus(y,np,lp,n,l,bigL)
            implicit none
            INTEGER, INTENT(IN) :: n, np, l, lp, bigL
            REAL(kind=8),    INTENT(IN) :: y
            REAL(kind=8) :: BesselElementplus
        end function BesselElementplus

    end INTERFACE

    INTEGER, INTENT(IN) :: np,lp,jp,n,l,j
    INTEGER, INTENT(IN) :: bigJ
    INTEGER :: bigL
    REAL(kind=8), INTENT(IN) :: y
    !REAL(kind=8) :: xjp,xj
    REAL(kind=8) :: Wigner_3j,Wigner_6j,Wigner_9j
    REAL(kind=8) :: PhiPPsummand2

    PhiPPsummand2 = 0.0

    if (l .ne. 0) then

        do bigL = bigJ,bigJ+1
            PhiPPsummand2= PhiPPsummand2+(-1.0)**(bigJ-bigL)*Jnorm(dble(bigL))*Wigner_6j(2*(bigJ+1),2,2*bigL,2,2*bigJ,2) &
                & *Wigner_6j(2*(bigJ+1),2,2*bigL,2*l,2*lp,2*(l-1))                                         &
                & *Wigner_9j(2*lp,2*l,2*bigL,1,1,2,jp,j,2*bigJ)
        end do

        PhiPPsummand2= PhiPPsummand2*Qnorm(dble(l)-1.0)*SQRT(dble(l))*Wigner_3j(2*lp,2*(bigJ+1),2*(l-1),0,0,0) &
            &*BesselElementplus(y,np,lp,n,l,bigJ+1)

    end if

end function PhiPPsummand2


!-------------------------------------------------------------------------------
function PhiPPsummand3(y,np,lp,jp,n,l,j,bigJ)
    implicit none
    INTERFACE

        function Jnorm(j)
            implicit none
            REAL(kind=8), INTENT(IN)  :: j
            REAL(kind=8) :: Jnorm
        end function Jnorm

        function Qnorm(j)
            implicit none
            REAL(kind=8), INTENT(IN)  :: j
            REAL(kind=8) :: Qnorm
        end function Qnorm

        function BesselElementplus(y,np,lp,n,l,bigL)
            implicit none
            INTEGER, INTENT(IN) :: n, np, l, lp, bigL
            REAL(kind=8), INTENT(IN) :: y
            REAL(kind=8) :: BesselElementplus
        end function BesselElementplus

    end INTERFACE

    INTEGER, INTENT(IN) :: np,lp,jp,n,l,j
    INTEGER, INTENT(IN) :: bigJ
    INTEGER :: bigL
    REAL(kind=8), INTENT(IN) :: y
    REAL(kind=8) :: Wigner_3j,Wigner_6j,Wigner_9j
    REAL(kind=8) :: PhiPPsummand3

    PhiPPsummand3 = 0.0

    if (bigJ .ne. 0) then
        do bigL = bigJ-1,bigJ
            PhiPPsummand3= PhiPPsummand3+(-1.0)**(bigJ-bigL+1.0)*Jnorm(dble(bigL))*Wigner_6j(2*(bigJ-1),2,2*bigL,2,2*bigJ,2) &
                & *Wigner_6j(2*(bigJ-1),2,2*bigL,2*l,2*lp,2*(l-1))                                         &
                & *Wigner_9j(2*lp,2*l,2*bigL,1,1,2,jp,j,2*bigJ)
        end do

        PhiPPsummand3= PhiPPsummand3*Qnorm(dble(l)-1.0)*SQRT(dble(l))*Wigner_3j(2*lp,2*(bigJ-1),2*(l-1),0,0,0) &
            &*BesselElementplus(y,np,lp,n,l,bigJ-1)

    end if

end function PhiPPsummand3


!-------------------------------------------------------------------------------
function PhiPPsummand4(y,np,lp,jp,n,l,j,bigJ)

    implicit none
    INTERFACE

        function Jnorm(j)
            implicit none
            REAL(kind=8), INTENT(IN)  :: j
            REAL(kind=8) :: Jnorm
        end function Jnorm

        function Qnorm(j)
            implicit none
            REAL(kind=8), INTENT(IN)  :: j
            REAL(kind=8) :: Qnorm
        end function Qnorm

        function BesselElementminus(y,np,lp,n,l,bigL)
            implicit none
            INTEGER, INTENT(IN) :: n, np, l, lp, bigL
            REAL(kind=8),    INTENT(IN) :: y
            REAL(kind=8) :: BesselElementminus
      end function BesselElementminus

    end INTERFACE

    INTEGER, INTENT(IN) :: np,lp,jp,n,l,j,bigJ
    REAL(kind=8), INTENT(IN) :: y
    INTEGER :: bigL
    REAL(kind=8) :: Wigner_3j,Wigner_6j,Wigner_9j
    REAL(kind=8) :: PhiPPsummand4

    PhiPPsummand4 = 0.0

    if (bigJ .ne. 0 .and. l .ne. 0) then
        do bigL = bigJ-1,bigJ
            PhiPPsummand4= PhiPPsummand4+(-1.0)**(bigJ-bigL)*Jnorm(dble(bigL))*Wigner_6j(2*(bigJ-1),2,2*bigL,2,2*bigJ,2) &
                  & *Wigner_6j(2*(bigJ-1),2,2*bigL,2*l,2*lp,2*(l+1))                                         &
                  & *Wigner_9j(2*lp,2*l,2*bigL,1,1,2,jp,j,2*bigJ)
        end do

        PhiPPsummand4= PhiPPsummand4*Qnorm(dble(l)+1.0)*SQRT(dble(l)+1.0)*Wigner_3j(2*lp,2*(bigJ-1),2*(l+1),0,0,0) &
            &*BesselElementminus(y,np,lp,n,l,bigJ-1)

    end if

end function PhiPPsummand4


!-------------------------------------------------------------------------------
function PhiPPoverall(lp,jp,l,j)

    implicit none
    INTERFACE

        function Jnorm(j)
            implicit none
            REAL(kind=8), INTENT(IN)  :: j
            REAL(kind=8) :: Jnorm
        end function Jnorm

    end INTERFACE

    INTEGER, INTENT(IN) :: lp,jp,l,j
    REAL(kind=8) :: Pi = 3.1415926535897932
    REAL(kind=8) :: xjp,xj
    REAL(kind=8) :: PhiPPoverall

    xjp = dble(jp)/2.0
    xj  = dble(j )/2.0
    PhiPPoverall = (-1)**(lp+1) * 6.0 * SQRT(Jnorm(xj)*Jnorm(xjp)*Jnorm(dble(lp))/(4*Pi))

end function PhiPPoverall



