! File contains dark matter response (dmr) functions roughly corresponding to
! equation (40) in the paper.
!

module dmresponse

    use response

    implicit none

contains

    ! Darkmatter particle response function:
    !   $R_M^{\tau\tau'}(\vec{v}_T^{\perp 2}, \frac{\vec{q}^2}{2m_N^2})$
    !   Eqn. 38 
    !
    !   This function comes from "dmformfactor-V6.m: FFfinal[DMmatrix_,J_,T)]".
    !   This is the coef. for $W_M^{\tau\tau'}(y)$.

    ! dmrMJ(tau1, tau2, q, v, jchi)
    ! dmrPhiPPJ(tau1, tau2, q, v, jchi)
    ! dmrPhiPPJMJ(tau1, tau2, q, v, jchi)
    ! dmrPhiTPJ(tau1, tau2, q, v, jchi)
    ! dmrSigmaPPJ(tau1, tau2, q, v, jchi)
    ! dmrSigmaPJ(tau1, tau2, q, v, jchi)
    ! dmrDeltaJ(tau1, tau2, q, v, jchi)
    ! dmrSigmaPJDeltaJ(tau1, tau2, q, v, jchi)

    ! From dmformfactor-V6.m:
    ! (*These response coefficients differ in their definition from those in the paper 
    ! because they have c's instead of a's, resulting in an overall factor of 
    ! 1/(4mN mchi)^2.  This will show up as an additional factor that we multiply
    ! FF by before printing out StrucFunction, TransitionProbability, and 
    ! ResponseNucl *)

!-------------------------------------------------------------------------------
function dmrMJ(tau1, tau2, q, v, jchi)
    ! O. Gorton 2020.01.
    ! import modules from modules.f90
    use masses

    implicit none

    ! outputs
    REAL(kind=8) :: dmrMJ

    ! inputs
    integer, INTENT(IN) :: tau1, tau2
    REAL(kind=8), INTENT(IN) :: q
    REAL(kind=8), INTENT(IN) :: v
    REAL(kind=8), INTENT(IN) :: jchi

    ! functions called
    !REAL(kind=8) :: Cl

    dmrMJ = 0.25*Cl(jchi) * ( &
            (cvec(tau1)%c(5)*cvec(tau2)%c(5)*q*q + cvec(tau1)%c(8)*cvec(tau2)%c(8)) &
            * (v*v - q*q/(4*muT*muT)) &
            + cvec(tau1)%c(11)*cvec(tau2)%c(11)*q*q &
        ) + (cvec(tau1)%c(1) + cvec(tau1)%c(2) * (v*v - q*q/(4*muT*muT))) * ( &
            cvec(tau2)%c(1) + cvec(tau2)%c(2) * (v*v - q*q/(4*muT*muT)) &
        )

end function dmrMJ

!-------------------------------------------------------------------------------
function dmrPhiPPJ(tau1, tau2, q, v, jchi)
    ! O. Gorton 2020.01.

    ! import modules from modules.f90
    use masses

    implicit none

    ! outputs
    REAL(kind=8) :: dmrPhiPPJ

    ! inputs
    REAL(kind=8), INTENT(IN) :: q
    REAL(kind=8), INTENT(IN) :: v
    REAL(kind=8), INTENT(IN) :: jchi
    integer, INTENT(IN) :: tau1, tau2

    ! functions called
    !REAL(kind=8) :: Cl

    dmrPhiPPJ = q*q/(16*mN*mN) * Cl(jchi) * (cvec(tau1)%c(12) - cvec(tau1)%c(15)*q*q) &
            * (cvec(tau2)%c(12) - cvec(tau2)%c(15)*q*q) + q**4.0/(4*mN*mN) * cvec(tau1)%c(3)*cvec(tau2)%c(3)

end function dmrPhiPPJ

!-------------------------------------------------------------------------------
function dmrPhiPPJMJ(tau1, tau2, q, v, jchi)
    ! O. Gorton 2020.01.

    ! import modules from modules.f90
    use masses

    implicit none

    ! outputs
    REAL(kind=8) :: dmrPhiPPJMJ

    ! inputs
    REAL(kind=8), INTENT(IN) :: q
    REAL(kind=8), INTENT(IN) :: v
    REAL(kind=8), INTENT(IN) :: jchi
    integer, INTENT(IN) :: tau1, tau2

    ! functions called
    !REAL(kind=8) :: Cl

    dmrPhiPPJMJ = q*q/(4*mN) * Cl(jchi) * cvec(tau1)%c(11) * ( &
                cvec(tau2)%c(12) - cvec(tau2)%c(15)*q*q &
            ) + q*q/(mN)*cvec(tau2)%c(3) * ( &
                cvec(tau1)%c(1) + cvec(tau1)%c(2) * ( &
                    v*v - q*q/(4*muT*muT) &
                ) &
        )

end function dmrPhiPPJMJ

!-------------------------------------------------------------------------------
function dmrPhiTPJ(tau1, tau2, q, v, jchi)
    ! O. Gorton 2020.01.

    ! import modules from modules.f90
    use masses

    implicit none

    ! outputs
    REAL(kind=8) :: dmrPhiTPJ

    ! inputs
    REAL(kind=8), INTENT(IN) :: q
    REAL(kind=8), INTENT(IN) :: v
    REAL(kind=8), INTENT(IN) :: jchi
    integer, INTENT(IN) :: tau1, tau2

    ! functions called
    !REAL(kind=8) :: Cl

    dmrPhiTPJ = q*q/(16*mN*mN)*Cl(jchi) * ( &
            cvec(tau1)%c(12)*cvec(tau2)%c(12)*q*q + cvec(tau1)%c(12)*cvec(tau2)%c(12) &
        )

end function dmrPhiTPJ

!-------------------------------------------------------------------------------
function dmrSigmaPPJ(tau1, tau2, q, v, jchi)
    ! O. Gorton 2020.01.

    ! import modules from modules.f90
    use masses

    implicit none

    ! outputs
    REAL(kind=8) :: dmrSigmaPPJ

    ! inputs
    REAL(kind=8), INTENT(IN) :: q
    REAL(kind=8), INTENT(IN) :: v
    REAL(kind=8), INTENT(IN) :: jchi
    integer, INTENT(IN) :: tau1, tau2

    ! functions called
    !REAL(kind=8) :: Cl

    dmrSigmaPPJ = (1./16) * Cl(jchi) * ( &
            cvec(tau1)%c(6)*cvec(tau2)%c(6)*q**4 + ( &
                cvec(tau1)%c(13)*cvec(tau2)%c(13)*q*q + cvec(tau1)%c(12)*cvec(tau2)%c(12) &
            ) * (v*v - q*q/(4*muT*muT)) &
            + 2*cvec(tau1)%c(4)*cvec(tau2)%c(6)*q*q + cvec(tau1)%c(4)*cvec(tau2)%c(4) &
        ) + (1./4)*cvec(tau1)%c(10)*cvec(tau2)%c(10)*q*q

end function dmrSigmaPPJ

!-------------------------------------------------------------------------------
function dmrSigmaPJ(tau1, tau2, q, v, jchi)
    ! O. Gorton 2020.01.

    ! import modules from modules.f90
    use masses

    implicit none

    ! outputs
    REAL(kind=8) :: dmrSigmaPJ

    ! inputs
    REAL(kind=8), INTENT(IN) :: q
    REAL(kind=8), INTENT(IN) :: v
    REAL(kind=8), INTENT(IN) :: jchi
    integer, INTENT(IN) :: tau1, tau2

    ! functions called
    !REAL(kind=8) :: Cl

    dmrSigmaPJ = (1./32) * Cl(jchi) * ( &
            2*cvec(tau1)%c(9)*cvec(tau2)%c(9)*q*q + ( &
                cvec(tau1)%c(15)*cvec(tau2)%c(15)*q**4 &
                + cvec(tau1)%c(14)*cvec(tau2)%c(14)*q*q &
                - 2*cvec(tau1)%c(12)*cvec(tau2)%c(15)*q*q &
                + cvec(tau1)%c(12)*cvec(tau2)%c(12) &
            ) * (v*v - q*q/(4*muT*muT)) + 2*cvec(tau1)%c(4)*cvec(tau2)%c(4) &
        ) + (1./8) * ( &
            cvec(tau1)%c(3)*cvec(tau2)%c(3)*q*q &
            + cvec(tau1)%c(7)*cvec(tau2)%c(7) &
        ) * (v*v - q*q/(4*muT*muT))

end function dmrSigmaPJ

!-------------------------------------------------------------------------------
function dmrDeltaJ(tau1, tau2, q, v, jchi)
    ! O. Gorton 2020.01.

    ! import modules from modules.f90
    use masses

    implicit none

    ! outputs
    REAL(kind=8) :: dmrDeltaJ

    ! inputs
    REAL(kind=8), INTENT(IN) :: q
    REAL(kind=8), INTENT(IN) :: v
    REAL(kind=8), INTENT(IN) :: jchi
    integer, INTENT(IN) :: tau1, tau2

    ! functions called
    !REAL(kind=8) :: Cl

    dmrDeltaJ = q*q/(4*mN*mN) * Cl(jchi) * ( &
            cvec(tau1)%c(5)*cvec(tau2)%c(5)*q*q &
            + cvec(tau1)%c(8)*cvec(tau2)%c(8) &
        ) + 2*q*q/(mN*mN) * cvec(tau1)%c(2)*cvec(tau2)%c(2) &
            * (v*v - q*q/(4*muT*muT))

end function dmrDeltaJ

!-------------------------------------------------------------------------------
function dmrSigmaPJDeltaJ(tau1, tau2, q, v, jchi)
    ! O. Gorton 2020.01.

    ! import modules from modules.f90
    use masses

    implicit none

    ! outputs
    REAL(kind=8) :: dmrSigmaPJDeltaJ

    ! inputs
    REAL(kind=8), INTENT(IN) :: q
    REAL(kind=8), INTENT(IN) :: v
    REAL(kind=8), INTENT(IN) :: jchi
    integer, INTENT(IN) :: tau1, tau2

    ! functions called
    !REAL(kind=8) :: Cl

    dmrSigmaPJDeltaJ = q*q/(4*mN*mN) * Cl(jchi) * ( &
            cvec(tau1)%c(4)*cvec(tau2)%c(5) - cvec(tau2)%c(8)*cvec(tau1)%c(9) &
        ) - q*q/(mN)*cvec(tau2)%c(2)*cvec(tau1)%c(3) * (v*v - q*q/(4*muT*muT))

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
