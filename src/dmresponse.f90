module dmresponse
    contains

function dmresponsefun(eft, term, tau1, tau2, q, v, jchi, muT)
    use kinds
    use parameters
    implicit none
    real(doublep), allocatable, intent(in) :: eft(:,:)
    integer :: term
    integer :: tau1, tau2
    real(doublep) :: q, v, jchi, muT
    REAL(doublep) :: dmresponsefun

    if (pndens) then
       dmresponsefun = dmresponse_pn(eft, term, tau1, tau2, q, v, jchi, muT)
    else
       dmresponsefun = dmresponse_iso(eft, term, tau1, tau2, q, v, jchi, muT)
    end if
    
end function

!----------------------------------------------------------------------------80
function dmresponse_iso(eft, ifunc, tau1, tau2, q, v, jchi, muT)
    ! Computes the dm response function in the proton-neutron coupling scheme
    ! Does so by transforming the isospin response function.

    implicit none
    real(kind=8) :: dmresponse_iso
    ! inputs 
    integer, intent(in) :: ifunc ! integer indicating response function
    ! ifunc     nuclear operator combination which coef belongs to
    ! 1         Mj Mj
    ! 2         PhiPPJ PhiPPJ
    ! 3         PhiTPJ PhiTPJ
    ! 4         DeltaJ DeltaJ
    ! 5         SigmaPJ SigmaPJ
    ! 6         SigmaPPJ SigmaPPJ
    ! 7         PhiPPJ MJ
    ! 8         DeltaJ SigmaPJ
    ! (passed to the real response coefficient functions)
    real(kind=8), allocatable, intent(in) :: eft(:,:)
    integer, INTENT(IN) :: tau1, tau2
    REAL(kind=8), INTENT(IN) :: q
    REAL(kind=8), INTENT(IN) :: v
    REAL(kind=8), INTENT(IN) :: jchi
    REAL(kind=8), INTENT(IN) :: muT    

    character(len=2) :: coupleto

    dmresponse_iso = -1000

    ! Determine isospin coupling by recombining the proton and neutron
    ! couplings. Recall that for p/n coupling:
    !   0 = proton
    !   1 = neutron
    write(coupleto,'(I1,I1)') tau1, tau2
    select case(coupleto)
      case('00')    
        dmresponse_iso = &
                   + dmresponse_pn(eft, ifunc, 0, 0, q, v, jchi, muT) &
                   + dmresponse_pn(eft, ifunc, 0, 1, q, v, jchi, muT) &
                   + dmresponse_pn(eft, ifunc, 1, 0, q, v, jchi, muT) &
                   + dmresponse_pn(eft, ifunc, 1, 1, q, v, jchi, muT)
      case('01')
        dmresponse_iso = &
                   + dmresponse_pn(eft, ifunc, 0, 0, q, v, jchi, muT) &
                   - dmresponse_pn(eft, ifunc, 0, 1, q, v, jchi, muT) &
                   + dmresponse_pn(eft, ifunc, 1, 0, q, v, jchi, muT) &
                   - dmresponse_pn(eft, ifunc, 1, 1, q, v, jchi, muT)
      case('10')
        dmresponse_iso = &
                   + dmresponse_pn(eft, ifunc, 0, 0, q, v, jchi, muT) &
                   + dmresponse_pn(eft, ifunc, 0, 1, q, v, jchi, muT) &
                   - dmresponse_pn(eft, ifunc, 1, 0, q, v, jchi, muT) &
                   - dmresponse_pn(eft, ifunc, 1, 1, q, v, jchi, muT)
      case('11')
        dmresponse_iso = &
                   + dmresponse_pn(eft, ifunc, 0, 0, q, v, jchi, muT) &
                   - dmresponse_pn(eft, ifunc, 0, 1, q, v, jchi, muT) &
                   - dmresponse_pn(eft, ifunc, 1, 0, q, v, jchi, muT) &
                   + dmresponse_pn(eft, ifunc, 1, 1, q, v, jchi, muT)
      case default
        stop "Invalid coupleto var."
    end select
    dmresponse_iso = 0.25d0*dmresponse_iso
    return

end function dmresponse_iso    

!----------------------------------------------------------------------------80
function dmresponse_pn(eft, ifunc, tau1, tau2, q, v, jchi, muT) 
    use parameters
    use dmresponselib
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

end module dmresponse
