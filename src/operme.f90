
!-------------------------------------------------------------------------------
! OperatorsME(j, operatorME) outputs the reduced matrix element of jth operator 
! (which correspond to the six operators and interferences in Eqn. 39).
!-----------------------------------------------------------------------------------

subroutine OperME(i,y,np,lp,jp,n,l,j,bigJ,operatorME)

    implicit none

    INTEGER, INTENT(IN) :: i
    INTEGER, INTENT(IN) :: np,lp,jp,n,l,j,bigJ
    REAL (kind=8), INTENT(IN)  :: y
    REAL (kind=8), INTENT(OUT) :: operatorME
    REAL (kind=8) :: MJ,SigmaJ,PhiPPoverall,Qnorm,PhiPPsummand1,PhiPPsummand2,PhiPPsummand3,PhiPPsummand4,&
                    &MJLDivQ,SigmaPJ,SigmaPPJ


    if (i .LT. 1 .or. i .GT. 7) then
        print *, 'The choice of operator j should be from 1 to 7'
        return
    ! i = 1, for operator MJ
    else if (i .eq. 1) then
        operatorME = MJ(y,np,lp,jp,n,l,j,bigJ)
    ! i = 2, for operator SigmaJ
    else if (i .eq. 2) then
        operatorME = SigmaJ(y,np,lp,jp,n,l,j,bigJ,bigJ)
    ! i = 3, for operator PhiPPJ
    else if (i .eq. 3) then
        operatorMe = Qnorm(dble(bigJ)+1.0)*SQRT(dble(bigJ)+1.0)*(PhiPPsummand1(y,np,lp,jp,n,l,j,bigJ) &
            & +PhiPPsummand2(y,np,lp,jp,n,l,j,bigJ))
        if (bigJ .ne. 0) then
            operatorMe = operatorMe + Qnorm(dble(bigJ)-1)*SQRT(dble(bigJ)) &
                 & *(PhiPPsummand3(y,np,lp,jp,n,l,j,bigJ)    &
                 & +PhiPPsummand4(y,np,lp,jp,n,l,j,bigJ))
        end if
    
        operatorME = operatorME * PhiPPoverall(lp,jp,l,j)
    
    ! i = 4, for operator PhiTPJ
    else if (i .eq. 4) then

        operatorME = -Qnorm(dble(bigJ)+1.0)*SQRT(dble(bigJ))*(PhiPPsummand1(y,np,lp,jp,n,l,j,bigJ) &
            & +PhiPPsummand2(y,np,lp,jp,n,l,j,bigJ))
        if (bigJ .ne. 0) then
            operatorMe = operatorMe + Qnorm(dble(bigJ)-1.0)*SQRT(dble(bigJ)+1.0) &
                 & *(PhiPPsummand3(y,np,lp,jp,n,l,j,bigJ)    &
                 & +PhiPPsummand4(y,np,lp,jp,n,l,j,bigJ))
        end if

        operatorME = operatorME * PhiPPoverall(lp,jp,l,j)
        operatorME = operatorME + SigmaJ(y,np,lp,jp,n,l,j,bigJ,bigJ)/2.0

    ! i = 5, for operator DeltaJ
    else if (i .eq. 5) then
        operatorME = MJLDivQ(y,np,lp,jp,n,l,j,bigJ,bigJ)

    ! i = 6, for operator SigmaPJ
    else if (i .eq. 6) then
        operatorME = SigmaPJ(y,np,lp,jp,n,l,j,bigJ)

    ! i = 7, for operator SigmaPPJ
    else if (i .eq. 7) then
        operatorME = SigmaPPJ(y,np,lp,jp,n,l,j,bigJ)

    else

    end if

end subroutine OperME
