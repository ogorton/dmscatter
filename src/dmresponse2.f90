function dmresponsecoef(eft, ifunc, tau1, tau2, q, v, jchi, muT)
    use parameters
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
    real(kind=8), allocatable, intent(in) :: eft(:,:)
    integer, INTENT(IN) :: tau1, tau2
    REAL(kind=8), INTENT(IN) :: q
    REAL(kind=8), INTENT(IN) :: v
    REAL(kind=8), INTENT(IN) :: jchi
    REAL(kind=8), INTENT(IN) :: muT

    dmresponsecoef = -1

    if (ifunc == 1) then
         dmresponsecoef = dmrMJ(eft, tau1, tau2, q, v, jchi, muT)
    else if (ifunc == 2) then
         dmresponsecoef = dmrPhiPPJ(eft, tau1, tau2, q, v, jchi, muT)
    else if (ifunc == 3) then
         dmresponsecoef = dmrPhiTPJ(eft, tau1, tau2, q, v, jchi, muT)
    else if (ifunc == 4) then
         dmresponsecoef = dmrDeltaJ(eft, tau1, tau2, q, v, jchi, muT)
    else if (ifunc == 5) then
         dmresponsecoef = dmrSigmaPJ(eft, tau1, tau2, q, v, jchi, muT)
    else if (ifunc == 6) then
         dmresponsecoef = dmrSigmaPPJ(eft, tau1, tau2, q, v, jchi, muT)
    else if (ifunc == 7) then
         dmresponsecoef = dmrPhiPPJMJ(eft, tau1, tau2, q, v, jchi, muT)
    else if (ifunc == 8) then
        dmresponsecoef = dmrSigmaPJDeltaJ(eft, tau1, tau2, q, v, jchi, muT)
    end if

    if (dmresponsecoef == -1) then
        print*,"No dark matter response coefficient was set!"
        STOP -1
    end if

end function dmresponsecoef
