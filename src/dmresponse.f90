! File contains dark matter response (dmr) functions roughly corresponding to
! equation (40) in the paper.
!

module dmresponse

    use constants
    use parameters
    implicit none

contains

    ! Darkmatter particle response function:
    !   $R_M^{\tau\tau'}(\vec{v}_T^{\perp 2}, \frac{\vec{q}^2}{2m_N^2})$
    !   Eqn. 38 
    !
    !   This function comes from "dmformfactor-V6.m: FFfinal[DMmatrix_,J_,T)]".
    !   This is the coef. for $W_M^{\tau\tau'}(y)$.

    ! dmrMJ(eft, tau1, tau2, q, v, jchi, muT)
    ! dmrPhiPPJ(eft, tau1, tau2, q, v, jchi, muT)
    ! dmrPhiPPJMJ(eft, tau1, tau2, q, v, jchi, muT)
    ! dmrPhiTPJ(eft, tau1, tau2, q, v, jchi, muT)
    ! dmrSigmaPPJ(eft, tau1, tau2, q, v, jchi, muT)
    ! dmrSigmaPJ(eft, tau1, tau2, q, v, jchi, muT)
    ! dmrDeltaJ(eft, tau1, tau2, q, v, jchi, muT)
    ! dmrSigmaPJDeltaJ(eft, tau1, tau2, q, v, jchi, muT)

    ! From dmformfactor-V6.m:
    ! (*These response coefficients differ in their definition from those in the paper 
    ! because they have c's instead of a's, resulting in an overall factor of 
    ! 1/(4mN mchi)^2.  This will show up as an additional factor that we multiply
    ! FF by before printing out StrucFunction, TransitionProbability, and 
    ! ResponseNucl *)

!-------------------------------------------------------------------------------
function dmrMJ(eft, tau1, tau2, q, v, jchi, muT)
    ! O. Gorton 2020.01.
    ! import modules from modules.f90

    implicit none

    ! outputs
    REAL(kind=8) :: dmrMJ

    ! inputs
    type(eftheory) :: eft
    integer, INTENT(IN) :: tau1, tau2
    REAL(kind=8), INTENT(IN) :: q
    REAL(kind=8), INTENT(IN) :: v
    REAL(kind=8), INTENT(IN) :: jchi, muT

    ! functions called
    !REAL(kind=8) :: Cl

    dmrMJ = 0.25*Cl(jchi) &
        * ( &
            (eft%isoc(tau1)%c(5)*eft%isoc(tau2)%c(5)*q*q &
                + eft%isoc(tau1)%c(8)*eft%isoc(tau2)%c(8)) &
            * (v*v - q*q/(4*muT*muT)) &
            + eft%isoc(tau1)%c(11)*eft%isoc(tau2)%c(11)*q*q &
        ) &
        + (eft%isoc(tau1)%c(1) + eft%isoc(tau1)%c(2) * (v*v - q*q/(4*muT*muT))) &
        * ( &
            eft%isoc(tau2)%c(1) + eft%isoc(tau2)%c(2) * (v*v - q*q/(4*muT*muT)) &
        )

end function dmrMJ

!-------------------------------------------------------------------------------
function dmrPhiPPJ(eft, tau1, tau2, q, v, jchi, muT)
    ! O. Gorton 2020.01.

    implicit none

    ! outputs
    REAL(kind=8) :: dmrPhiPPJ

    ! inputs
    type(eftheory) :: eft
    REAL(kind=8), INTENT(IN) :: q
    REAL(kind=8), INTENT(IN) :: v
    REAL(kind=8), INTENT(IN) :: jchi, muT
    integer, INTENT(IN) :: tau1, tau2

    ! functions called
    !REAL(kind=8) :: Cl

    dmrPhiPPJ = q*q/(16*mN*mN) * Cl(jchi) * (eft%isoc(tau1)%c(12) &
            - eft%isoc(tau1)%c(15)*q*q) &
            * (eft%isoc(tau2)%c(12) - eft%isoc(tau2)%c(15)*q*q) &
            + q**4.0/(4*mN*mN) * eft%isoc(tau1)%c(3)*eft%isoc(tau2)%c(3)

end function dmrPhiPPJ

!-------------------------------------------------------------------------------
function dmrPhiPPJMJ(eft, tau1, tau2, q, v, jchi, muT)
    ! O. Gorton 2020.01.

    implicit none

    ! outputs
    REAL(kind=8) :: dmrPhiPPJMJ

    ! inputs
    type(eftheory) :: eft
    REAL(kind=8), INTENT(IN) :: q
    REAL(kind=8), INTENT(IN) :: v
    REAL(kind=8), INTENT(IN) :: jchi, muT
    integer, INTENT(IN) :: tau1, tau2

    ! functions called
    !REAL(kind=8) :: Cl

    dmrPhiPPJMJ = q*q/(4*mN) * Cl(jchi) * eft%isoc(tau1)%c(11) * ( &
                eft%isoc(tau2)%c(12) - eft%isoc(tau2)%c(15)*q*q &
            ) + q*q/(mN)*eft%isoc(tau2)%c(3) * ( &
                eft%isoc(tau1)%c(1) + eft%isoc(tau1)%c(2) * ( &
                    v*v - q*q/(4*muT*muT) &
                ) &
        )

end function dmrPhiPPJMJ

!-------------------------------------------------------------------------------
function dmrPhiTPJ(eft, tau1, tau2, q, v, jchi, muT)
    ! O. Gorton 2020.01.

    implicit none

    ! outputs
    REAL(kind=8) :: dmrPhiTPJ

    ! inputs
    type(eftheory) :: eft
    REAL(kind=8), INTENT(IN) :: q
    REAL(kind=8), INTENT(IN) :: v
    REAL(kind=8), INTENT(IN) :: jchi, muT
    integer, INTENT(IN) :: tau1, tau2

    ! functions called
    !REAL(kind=8) :: Cl

    dmrPhiTPJ = q*q/(16*mN*mN)*Cl(jchi) * ( &
            eft%isoc(tau1)%c(12)*eft%isoc(tau2)%c(12)*q*q &
            + eft%isoc(tau1)%c(12)*eft%isoc(tau2)%c(12) &
        )

end function dmrPhiTPJ

!-------------------------------------------------------------------------------
function dmrSigmaPPJ(eft, tau1, tau2, q, v, jchi, muT)
    ! O. Gorton 2020.01.

    implicit none

    ! outputs
    REAL(kind=8) :: dmrSigmaPPJ

    ! inputs
    type(eftheory) :: eft
    REAL(kind=8), INTENT(IN) :: q
    REAL(kind=8), INTENT(IN) :: v
    REAL(kind=8), INTENT(IN) :: jchi, muT
    integer, INTENT(IN) :: tau1, tau2

    ! functions called
    !REAL(kind=8) :: Cl

    dmrSigmaPPJ = (1./16) * Cl(jchi) * ( &
            eft%isoc(tau1)%c(6)*eft%isoc(tau2)%c(6)*q**4 + ( &
                eft%isoc(tau1)%c(13)*eft%isoc(tau2)%c(13)*q*q &
                + eft%isoc(tau1)%c(12)*eft%isoc(tau2)%c(12) &
            ) * (v*v - q*q/(4*muT*muT)) &
            + 2*eft%isoc(tau1)%c(4)*eft%isoc(tau2)%c(6)*q*q &
            + eft%isoc(tau1)%c(4)*eft%isoc(tau2)%c(4) &
        ) + (1./4)*eft%isoc(tau1)%c(10)*eft%isoc(tau2)%c(10)*q*q

end function dmrSigmaPPJ

!-------------------------------------------------------------------------------
function dmrSigmaPJ(eft, tau1, tau2, q, v, jchi, muT)
    ! O. Gorton 2020.01.

    implicit none

    ! outputs
    REAL(kind=8) :: dmrSigmaPJ

    ! inputs
    type(eftheory) :: eft
    REAL(kind=8), INTENT(IN) :: q
    REAL(kind=8), INTENT(IN) :: v
    REAL(kind=8), INTENT(IN) :: jchi, muT
    integer, INTENT(IN) :: tau1, tau2

    ! functions called
    !REAL(kind=8) :: Cl

    dmrSigmaPJ = (1./32) * Cl(jchi) * ( &
            2*eft%isoc(tau1)%c(9)*eft%isoc(tau2)%c(9)*q*q + ( &
                eft%isoc(tau1)%c(15)*eft%isoc(tau2)%c(15)*q**4 &
                + eft%isoc(tau1)%c(14)*eft%isoc(tau2)%c(14)*q*q &
                - 2*eft%isoc(tau1)%c(12)*eft%isoc(tau2)%c(15)*q*q &
                + eft%isoc(tau1)%c(12)*eft%isoc(tau2)%c(12) &
            ) * (v*v - q*q/(4*muT*muT)) &
            + 2*eft%isoc(tau1)%c(4)*eft%isoc(tau2)%c(4) &
        ) + (1./8) * ( &
            eft%isoc(tau1)%c(3)*eft%isoc(tau2)%c(3)*q*q &
            + eft%isoc(tau1)%c(7)*eft%isoc(tau2)%c(7) &
        ) * (v*v - q*q/(4*muT*muT))

end function dmrSigmaPJ

!-------------------------------------------------------------------------------
function dmrDeltaJ(eft, tau1, tau2, q, v, jchi, muT)
    ! O. Gorton 2020.01.
    implicit none

    ! outputs
    REAL(kind=8) :: dmrDeltaJ

    ! inputs
    type(eftheory) :: eft
    REAL(kind=8), INTENT(IN) :: q
    REAL(kind=8), INTENT(IN) :: v
    REAL(kind=8), INTENT(IN) :: jchi, muT
    integer, INTENT(IN) :: tau1, tau2

    ! functions called
    !REAL(kind=8) :: Cl

    dmrDeltaJ = q*q/(4*mN*mN) * Cl(jchi) * ( &
            eft%isoc(tau1)%c(5)*eft%isoc(tau2)%c(5)*q*q &
            + eft%isoc(tau1)%c(8)*eft%isoc(tau2)%c(8) &
        ) + 2*q*q/(mN*mN) * eft%isoc(tau1)%c(2)*eft%isoc(tau2)%c(2) &
            * (v*v - q*q/(4*muT*muT))

end function dmrDeltaJ

!-------------------------------------------------------------------------------
function dmrSigmaPJDeltaJ(eft, tau1, tau2, q, v, jchi, muT)
    ! O. Gorton 2020.01.

    implicit none

    ! outputs
    REAL(kind=8) :: dmrSigmaPJDeltaJ

    ! inputs
    type(eftheory) :: eft
    REAL(kind=8), INTENT(IN) :: q
    REAL(kind=8), INTENT(IN) :: v
    REAL(kind=8), INTENT(IN) :: jchi, muT
    integer, INTENT(IN) :: tau1, tau2

    ! functions called
    !REAL(kind=8) :: Cl

    dmrSigmaPJDeltaJ = q*q/(4*mN*mN) * Cl(jchi) * ( &
            eft%isoc(tau1)%c(4)*eft%isoc(tau2)%c(5) - eft%isoc(tau2)%c(8)*eft%isoc(tau1)%c(9) &
        ) - q*q/(mN)*eft%isoc(tau2)%c(2)*eft%isoc(tau1)%c(3) * (v*v - q*q/(4*muT*muT))

end function dmrSigmaPJDeltaJ

!-------------------------------------------------------------------------------
function Cl(j)
    ! O. Gorton 2020.01.
    ! Abbreviation for 4j(j+1)/3
    implicit none
    REAL(kind=8), INTENT(IN) :: j
    REAL(kind=8) Cl

    Cl = 4.0 * j * (j + 1) / 3.0

end function Cl

end module dmresponse
