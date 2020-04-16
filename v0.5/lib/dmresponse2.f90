function dmresponsecoef(ifunc, tau1, tau2, q, v, jchi)
    use dmresponse
    implicit none
    ! output
    real(kind=8) :: dmresponsecoef
    ! inputs 
    integer, intent(in) :: ifunc ! integer indicating response function
    ! ifunc	nuclear operator combination which coef belongs to
    ! 1 	Mj Mj
    ! 2		PhiPPJ PhiPPJ
    ! 3		PhiTPJ PhiTPJ
    ! 4		DeltaJ DeltaJ
    ! 5		SigmaPJ SigmaPJ
    ! 6		SigmaPPJ SigmaPPJ
    ! 7		PhiPPJ MJ
    ! 8		DeltaJ SigmaPJ
    ! (passed to the real response coefficient functions)
    integer, INTENT(IN) :: tau1, tau2
    REAL(kind=8), INTENT(IN) :: q
    REAL(kind=8), INTENT(IN) :: v
    REAL(kind=8), INTENT(IN) :: jchi
    ! functions in dmresponse.f90
!    real(kind=8) :: dmrMJ, dmrPhiPPJ, dmrPhiTPJ, dmrDeltaJ
!    real(kind=8) :: dmrSigmaPJ, dmrSigmaPPJ, dmrPhiPPJMJ, dmrSigmaPJDeltaJ

    dmresponsecoef = -1

    if (ifunc == 1) dmresponsecoef = dmrMJ(tau1, tau2, q, v, jchi)
    if (ifunc == 2) dmresponsecoef = dmrPhiPPJ(tau1, tau2, q, v, jchi)
    if (ifunc == 3) dmresponsecoef = dmrPhiTPJ(tau1, tau2, q, v, jchi)
    if (ifunc == 4) dmresponsecoef = dmrDeltaJ(tau1, tau2, q, v, jchi)
    if (ifunc == 5) dmresponsecoef = dmrSigmaPJ(tau1, tau2, q, v, jchi)
    if (ifunc == 6) dmresponsecoef = dmrSigmaPPJ(tau1, tau2, q, v, jchi)
    if (ifunc == 7) dmresponsecoef = dmrPhiPPJMJ(tau1, tau2, q, v, jchi)
    if (ifunc == 8) dmresponsecoef = dmrSigmaPJDeltaJ(tau1, tau2, q, v, jchi)

    if (dmresponsecoef == -1) then
        print*,"No dark matter response coefficient was set!"
        STOP -1
    end if

end function dmresponsecoef
