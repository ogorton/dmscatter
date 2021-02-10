function dmresponse_pn(eft, ifunc, tau1, tau2, q, v, jchi, muT) 
    use parameters
    use dmresponse
    implicit none
    ! output
    real(kind=8) :: dmresponse_pn
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

    select case (ifunc)
      case(1)
         dmresponse_pn = dmrMJ(eft, tau1, tau2, q, v, jchi, muT)
      case(2)
         dmresponse_pn = dmrPhiPPJ(eft, tau1, tau2, q, v, jchi, muT)
      case(3)
         dmresponse_pn = dmrPhiTPJ(eft, tau1, tau2, q, v, jchi, muT)
      case(4)
         dmresponse_pn = dmrDeltaJ(eft, tau1, tau2, q, v, jchi, muT)
      case(5)
         dmresponse_pn = dmrSigmaPJ(eft, tau1, tau2, q, v, jchi, muT)
      case(6)
         dmresponse_pn = dmrSigmaPPJ(eft, tau1, tau2, q, v, jchi, muT)
      case(7)
         dmresponse_pn = dmrPhiPPJMJ(eft, tau1, tau2, q, v, jchi, muT)
      case(8)
        dmresponse_pn = dmrSigmaPJDeltaJ(eft, tau1, tau2, q, v, jchi, muT)
      case default
        print*,"No dark matter response coefficient was set!"
        STOP -1
    end select        

end function dmresponse_pn
