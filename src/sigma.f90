!-------------------------------------------------------------------------------
function SigmaJ(y,np,lp,jp,n,l,j,bigJ,bigL)
    use sj2iref
    implicit none
    INTERFACE

        function Jnorm(j)
            implicit none
            REAL(kind=8), INTENT(IN)  :: j
            REAL(kind=8) :: Jnorm
        end function Jnorm

        function BesselElement(y,np,lp,n,l,bigL)
            implicit none
            INTEGER, INTENT(IN) :: n, np, l, lp, bigL
            REAL(kind=8), INTENT(IN) :: y
            REAL(kind=8) :: BesselElement
        end function BesselElement

    end INTERFACE

    INTEGER, INTENT(IN) :: np,lp,jp,n,l,j
    INTEGER, INTENT(IN) :: bigJ,bigL
    REAL(kind=8), INTENT(IN) :: y
    REAL(kind=8) :: xjp,xj
    REAL(kind=8) :: Pi = 3.14159265358979323
    REAL(kind=8) :: Wigner_9j
    REAL(kind=8) :: SigmaJ

    xjp = dble(jp)/2.0
    xj  = dble(j )/2.0

    SigmaJ = (-1.0)**lp * SQRT(Jnorm(xj)*Jnorm(xjp)*Jnorm(dble(l))*Jnorm(dble(lp))  &
         & *Jnorm(dble(bigJ))*Jnorm(dble(bigL))/(4*Pi))     &
         & * SQRT(6.0) * tj2i_lookup(2*lp,2*bigL,2*l,0,0,0) &
         & * Wigner_9j(2*lp,2*l,2*bigL,1,1,2,jp,j,2*bigJ) &
         & * BesselElement(y,np,lp,n,l,bigL)

end function SigmaJ


!-------------------------------------------------------------------------------
function SigmaPJ(y,np,lp,jp,n,l,j,bigJ)

    INTERFACE

        function Qnorm(j)
            implicit none
            REAL(kind=8), INTENT(IN)  :: j
            REAL(kind=8) :: Qnorm
        end function Qnorm

        function SigmaJ(y,np,lp,jp,n,l,j,bigJ,bigL)
            implicit none
            INTEGER, INTENT(IN) :: np,lp,jp,n,l,j
            INTEGER, INTENT(IN) :: bigJ,bigL
            REAL(kind=8), INTENT(IN) :: y
            REAL(kind=8) :: SigmaJ
        end function SigmaJ

  end INTERFACE

  INTEGER, INTENT(IN) :: np,lp,jp,n,l,j
  INTEGER, INTENT(IN) :: bigJ
  REAL(kind=8), INTENT(IN) :: y
  REAL(kind=8) :: SigmaPJ

  SigmaPJ = -SQRT(dble(bigJ))/Qnorm(dble(bigJ))*SigmaJ(y,np,lp,jp,n,l,j,bigJ,bigJ+1)   &
      & +SQRT(dble(bigJ)+1.0)/Qnorm(dble(bigJ))*SigmaJ(y,np,lp,jp,n,l,j,bigJ,bigJ-1)

end function SigmaPJ


!-------------------------------------------------------------------------------
function SigmaPPJ(y,np,lp,jp,n,l,j,bigJ)

    INTERFACE

        function Qnorm(j)
            implicit none
            REAL(kind=8), INTENT(IN)  :: j
            REAL(kind=8) :: Qnorm
        end function Qnorm

        function SigmaJ(y,np,lp,jp,n,l,j,bigJ,bigL)
            implicit none
            INTEGER, INTENT(IN) :: np,lp,jp,n,l,j
            INTEGER, INTENT(IN) :: bigJ,bigL
            REAL(kind=8), INTENT(IN) :: y
            REAL(kind=8) :: SigmaJ
       end function SigmaJ

    end INTERFACE

    INTEGER, INTENT(IN) :: np,lp,jp,n,l,j
    INTEGER, INTENT(IN) :: bigJ
    REAL(kind=8), INTENT(IN) :: y
    REAL(kind=8) :: SigmaPPJ

    SigmaPPJ = SQRT(dble(bigJ)+1.0)/Qnorm(dble(bigJ))*SigmaJ(y,np,lp,jp,n,l,j,bigJ,bigJ+1)   &
        & +SQRT(dble(bigJ))/Qnorm(dble(bigJ))*SigmaJ(y,np,lp,jp,n,l,j,bigJ,bigJ-1)

end function SigmaPPJ



