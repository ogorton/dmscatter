! File contains dark matter response (dmr) functions roughly corresponding to
! equation (40) in the paper.
!

module dmresponse
    implicit none

contains

    ! From dmformfactor-V6.m:
    ! (*These response coefficients differ in their definition from those in the paper 
    ! because they have c's instead of a's, resulting in an overall factor of 
    ! 1/(4mN mchi)^2.  This will show up as an additional factor that we multiply
    ! FF by before printing out StrucFunction, TransitionProbability, and 
    ! ResponseNucl *)

!-------------------------------------------------------------------------------
function dmrMJ(q, v, jchi)
    ! O. Gorton 2020.01.
    ! Darkmatter particle response function:
    !	$R_M^{\tau\tau'}(\vec{v}_T^{\perp 2}, \frac{\vec{q}^2}{2m_N^2})$
    !	Eqn. 38 
    !
    !	This function comes from "dmformfactor-V6.m: FFfinal[DMmatrix_,J_,T)]".
    !	This is the coef. for $W_M^{\tau\tau'}(y)$.

    ! import modules from modules.f90
    use masses

    implicit none

    ! outputs
    REAL(kind=8) :: dmrMJ

    ! inputs
    REAL(kind=8) :: c5, c5p, c8, c8p, c11, c11p, c1, c1p, c2, c2p
    REAL(kind=8), INTENT(IN) :: q
    REAL(kind=8), INTENT(IN) :: v
    REAL(kind=8), INTENT(IN) :: jchi

    ! functions called
    !REAL(kind=8) :: Cl

    dmrMJ = 0.25*Cl(jchi) * ( &
            (c5*c5p*q*q + c8*c8p) &
            * (v*v - q*q/(4*muT*muT)) &
            + c11*c11p*q*q &
        ) + (c1 + c2 * (v*v - q*q/(4*muT*muT))) * ( &
            c1p + c2p * (v*v - q*q/(4*muT*muT)) &
        )

end function dmrMJ

!-------------------------------------------------------------------------------
function dmrPhiPPJ(q, v, jchi)
    ! O. Gorton 2020.01.
    ! Darkmatter particle response function:
    !   $R_{\Phi''}^{\tau\tau'}(\vec{v}_T^{\perp 2}, \frac{\vec{q}^2}{2m_N^2})$
    !   Eqn. 38 
    !
    !   This function comes from "dmformfactor-V6.m: FFfinal[DMmatrix_,J_,T)]".
    !   This is the coef. for $W_{\Phi''}^{\tau\tau'}(y)$.

    ! import modules from modules.f90
    use masses

    implicit none

    ! outputs
    REAL(kind=8) :: dmrPhiPPJ

    ! inputs
    REAL(kind=8) :: c12, c12p, c15, c15p, c3, c3p
    REAL(kind=8), INTENT(IN) :: q
    REAL(kind=8), INTENT(IN) :: v
    REAL(kind=8), INTENT(IN) :: jchi

    ! functions called
    !REAL(kind=8) :: Cl

    dmrPhiPPJ = q*q/(16*mN*mN) * Cl(jchi) * (c12 - c15*q*q) &
            * (c12p - c15p*q*q) + q**4.0/(4*mN*mN) * c3*c3p

end function dmrPhiPPJ

!-------------------------------------------------------------------------------
function dmrPhiPPJMJ(q, v, jchi)
    ! O. Gorton 2020.01.
    ! Darkmatter particle response function:
    !   $R_{\Phi''M}^{\tau\tau'}(\vec{v}_T^{\perp 2}, \frac{\vec{q}^2}{2m_N^2})$
    !   Eqn. 38 
    !
    !   This function comes from "dmformfactor-V6.m: FFfinal[DMmatrix_,J_,T)]".
    !   This is the coef. for $W_{\Phi''M}^{\tau\tau'}(y)$.

    ! import modules from modules.f90
    use masses

    implicit none

    ! outputs
    REAL(kind=8) :: dmrPhiPPJMJ

    ! inputs
    REAL(kind=8) :: c1, c2, c12, c12p, c15, c15p, c3, c3p, c11, c11p
    REAL(kind=8), INTENT(IN) :: q
    REAL(kind=8), INTENT(IN) :: v
    REAL(kind=8), INTENT(IN) :: jchi

    ! functions called
    !REAL(kind=8) :: Cl

    dmrPhiPPJMJ = q*q/(4*mN) * Cl(jchi) * c11 * ( &
                c12p - c15p*q*q &
            ) + q*q/(mN)*c3p * ( &
                c1 + c2 * ( &
                    v*v - q*q/(4*muT*muT) &
                ) &
            )

end function dmrPhiPPJMJ

!-------------------------------------------------------------------------------
function dmrPhiTPJ(q, v, jchi)
    ! O. Gorton 2020.01.
    ! Darkmatter particle response function:
    !   $R_{\tilde{\Phi}'M}^{\tau\tau'}(\vec{v}_T^{\perp 2}, \frac{\vec{q}^2}{2m_N^2})$
    !   Eqn. 38 
    !   
    !   This function comes from "dmformfactor-V6.m: FFfinal[DMmatrix_,J_,T)]".
    !   This is the coef. for $W_{\tilde{\Phi}'M}^{\tau\tau'}(y)$.

    ! import modules from modules.f90
    use masses

    implicit none

    ! outputs
    REAL(kind=8) :: dmrPhiTPJ

    ! inputs
    REAL(kind=8) :: c1, c2, c12, c12p, c15, c15p, c3, c3p, c11, c11p
    REAL(kind=8), INTENT(IN) :: q
    REAL(kind=8), INTENT(IN) :: v
    REAL(kind=8), INTENT(IN) :: jchi

    ! functions called
    !REAL(kind=8) :: Cl

    dmrPhiTPJ = q*q/(16*mN*mN)*Cl(jchi) * ( &
            c12*c12p*q*q + c12*c12p &
        )

end function dmrPhiTPJ

!-------------------------------------------------------------------------------
function dmrSigmaPPJ(q, v, jchi)
    ! O. Gorton 2020.01.
    ! Darkmatter particle response function:
    !   $R_{\tilde{\Phi}'M}^{\tau\tau'}(\vec{v}_T^{\perp 2}, \frac{\vec{q}^2}{2m_N^2})$
    !   Eqn. 38 
    !   
    !   This function comes from "dmformfactor-V6.m: FFfinal[DMmatrix_,J_,T)]".
    !   This is the coef. for $W_{\tilde{\Phi}'M}^{\tau\tau'}(y)$.

    ! import modules from modules.f90
    use masses

    implicit none

    ! outputs
    REAL(kind=8) :: dmrSigmaPPJ

    ! inputs
    REAL(kind=8) :: c6, c6p, c13, c13p, c12, c12p, c4, c4p, c10, c10p
    REAL(kind=8), INTENT(IN) :: q
    REAL(kind=8), INTENT(IN) :: v
    REAL(kind=8), INTENT(IN) :: jchi

    ! functions called
    !REAL(kind=8) :: Cl

    dmrSigmaPPJ = (1./16) * Cl(jchi) * ( &
            c6*c6p*q**4 + ( &
                c13*c13p*q*q + c12*c12p &
            ) * (v*v - q*q/(4*muT*muT)) &
            + 2*c4*c6p*q*q + c4*c4p &
        ) + (1./4)*c10*c10p*q*q

end function dmrSigmaPPJ

!-------------------------------------------------------------------------------
function dmrSigmaPJ(q, v, jchi)
    ! O. Gorton 2020.01.
    ! Darkmatter particle response function:
    !   $R_{\tilde{\Phi}'M}^{\tau\tau'}(\vec{v}_T^{\perp 2}, \frac{\vec{q}^2}{2m_N^2})$
    !   Eqn. 38 
    !   
    !   This function comes from "dmformfactor-V6.m: FFfinal[DMmatrix_,J_,T)]".
    !   This is the coef. for $W_{\tilde{\Phi}'M}^{\tau\tau'}(y)$.

    ! import modules from modules.f90
    use masses

    implicit none

    ! outputs
    REAL(kind=8) :: dmrSigmaPJ

    ! inputs
    REAL(kind=8) :: c9, c9p, c15, c15p, c14, c14p, c12, c12p, c4, c4p, c3, c3p
    REAL(kind=8) :: c7, c7p
    REAL(kind=8), INTENT(IN) :: q
    REAL(kind=8), INTENT(IN) :: v
    REAL(kind=8), INTENT(IN) :: jchi

    ! functions called
    !REAL(kind=8) :: Cl

    dmrSigmaPJ = (1./32) * Cl(jchi) * ( &
            2*c9*c9p*q*q + ( &
                c15*c15p*q**4 + c14*c14p*q*q - 2*c12*c15p*q*q + c12*c12p &
            ) * (v*v - q*q/(4*muT*muT)) + 2*c4*c4p &
        ) + (1./8) * ( &
            c3*c3p*q*q + c7*c7p &
        ) * (v*v - q*q/(4*muT*muT))

end function dmrSigmaPJ

!-------------------------------------------------------------------------------
function dmrDeltaJ(q, v, jchi)
    ! O. Gorton 2020.01.
    ! Darkmatter particle response function:
    !   $R_{\tilde{\Phi}'M}^{\tau\tau'}(\vec{v}_T^{\perp 2}, \frac{\vec{q}^2}{2m_N^2})$
    !   Eqn. 38 
    !   
    !   This function comes from "dmformfactor-V6.m: FFfinal[DMmatrix_,J_,T)]".
    !   This is the coef. for $W_{\tilde{\Phi}'M}^{\tau\tau'}(y)$.

    ! import modules from modules.f90
    use masses

    implicit none

    ! outputs
    REAL(kind=8) :: dmrDeltaJ

    ! inputs
    REAL(kind=8) :: c5, c5p, c8, c8p, c2, c2p
    REAL(kind=8), INTENT(IN) :: q
    REAL(kind=8), INTENT(IN) :: v
    REAL(kind=8), INTENT(IN) :: jchi

    ! functions called
    !REAL(kind=8) :: Cl

    dmrDeltaJ = q*q/(4*mN*mN) * Cl(jchi) * ( &
            c5*c5p*q*q + c8*c8p &
        ) + 2*q*q/(mN*mN) * c2*c2p * (v*v - q*q/(4*muT*muT))

end function dmrDeltaJ

!-------------------------------------------------------------------------------
function dmrSigmaPJDeltaJ(q, v, jchi)
    ! O. Gorton 2020.01.
    ! Darkmatter particle response function:
    !   $R_{\tilde{\Phi}'M}^{\tau\tau'}(\vec{v}_T^{\perp 2}, \frac{\vec{q}^2}{2m_N^2})$
    !   Eqn. 38 
    !   
    !   This function comes from "dmformfactor-V6.m: FFfinal[DMmatrix_,J_,T)]".
    !   This is the coef. for $W_{\tilde{\Phi}'M}^{\tau\tau'}(y)$.

    ! import modules from modules.f90
    use masses

    implicit none

    ! outputs
    REAL(kind=8) :: dmrSigmaPJDeltaJ

    ! inputs
    REAL(kind=8) :: c4, c5p, c8p, c9, c2p, c3
    REAL(kind=8), INTENT(IN) :: q
    REAL(kind=8), INTENT(IN) :: v
    REAL(kind=8), INTENT(IN) :: jchi

    ! functions called
    !REAL(kind=8) :: Cl

    dmrSigmaPJDeltaJ = q*q/(4*mN*mN) * Cl(jchi) * ( &
            c4*c5p - c8p*c9 &
        ) - q*q/(mN)*c2p*c3 * (v*v - q*q/(4*muT*muT))

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
