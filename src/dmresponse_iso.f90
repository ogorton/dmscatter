function dmresponse_iso(eft, ifunc, tau1, tau2, q, v, jchi, muT)
    ! Computes the dm response function in the proton-neutron coupling scheme
    ! Does so by transforming the isospin response function.

    implicit none
    interface
        function dmresponse_pn(eft, term, tau1, tau2, q, v, jchi, muT)
            use kinds
            use parameters
            real(doublep), allocatable, intent(in) :: eft(:,:)
            integer :: term
            integer :: tau1, tau2
            real(doublep) :: q, v, jchi, muT
            REAL(doublep) :: dmresponse_pn
        end function
    end interface
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
    dmresponse_iso = -1000

    ! Determine isospin coupling by recombining the proton and neutron
    ! couplings. Recall that for p/n coupling:
    !   0 = proton
    !   1 = neutron

    if (tau1==0 .and. tau2==0) then
        dmresponse_iso = &
                     dmresponse_pn(eft, ifunc, 0, 0, q, v, jchi, muT) &
                   + dmresponse_pn(eft, ifunc, 0, 1, q, v, jchi, muT) &
                   + dmresponse_pn(eft, ifunc, 1, 0, q, v, jchi, muT) &
                   + dmresponse_pn(eft, ifunc, 1, 1, q, v, jchi, muT)
    end if
    if (tau1==0 .and. tau2==1) then
        dmresponse_iso = &
                     dmresponse_pn(eft, ifunc, 0, 0, q, v, jchi, muT) &
                   - dmresponse_pn(eft, ifunc, 0, 1, q, v, jchi, muT) &
                   + dmresponse_pn(eft, ifunc, 1, 0, q, v, jchi, muT) &
                   - dmresponse_pn(eft, ifunc, 1, 1, q, v, jchi, muT)
    end if
    if (tau1==1 .and. tau2==0) then
        dmresponse_iso = &
                     dmresponse_pn(eft, ifunc, 0, 0, q, v, jchi, muT) &
                   + dmresponse_pn(eft, ifunc, 0, 1, q, v, jchi, muT) &
                   - dmresponse_pn(eft, ifunc, 1, 0, q, v, jchi, muT) &
                   - dmresponse_pn(eft, ifunc, 1, 1, q, v, jchi, muT)
    end if
    if (tau1==1 .and. tau2==1) then
        dmresponse_iso = &
                     dmresponse_pn(eft, ifunc, 0, 0, q, v, jchi, muT) &
                   - dmresponse_pn(eft, ifunc, 0, 1, q, v, jchi, muT) &
                   - dmresponse_pn(eft, ifunc, 1, 0, q, v, jchi, muT) &
                   + dmresponse_pn(eft, ifunc, 1, 1, q, v, jchi, muT)
    end if
    dmresponse_iso = 0.25d0*dmresponse_iso
    return

end function dmresponse_iso    
