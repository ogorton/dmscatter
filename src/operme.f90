
!-------------------------------------------------------------------------------
! OperatorsME(j, operatorME) outputs the reduced matrix element of jth operator 
! (which correspond to the six operators and interferences in Eqn. 39).
!-----------------------------------------------------------------------------------

subroutine OperME(i,y,np,lp,jp,n,l,j,bigJ,operatorME)

    use mjl
    use norm
    use phi
    use sigma

    implicit none

    INTEGER, INTENT(IN) :: i
    INTEGER, INTENT(IN) :: np,lp,jp,n,l,j,bigJ
    REAL (kind=8), INTENT(IN)  :: y
    REAL (kind=8), INTENT(OUT) :: operatorME

    select case(i)
    case(1)
        operatorME = MJ(y,np,lp,jp,n,l,j,bigJ)
    case(2)
        operatorME = SigmaJ(y,np,lp,jp,n,l,j,bigJ,bigJ)
    case(3)
        operatorMe = Qnorm(dble(bigJ)+1.0)*SQRT(dble(bigJ)+1.0)*(PhiPPsummand1(y,np,lp,jp,n,l,j,bigJ) &
            & +PhiPPsummand2(y,np,lp,jp,n,l,j,bigJ))
        if (bigJ .ne. 0) then
            operatorMe = operatorMe + Qnorm(dble(bigJ)-1)*SQRT(dble(bigJ)) &
                 & *(PhiPPsummand3(y,np,lp,jp,n,l,j,bigJ)    &
                 & +PhiPPsummand4(y,np,lp,jp,n,l,j,bigJ))
        end if
        operatorME = operatorME * PhiPPoverall(lp,jp,l,j)
    ! i = 4, for operator PhiTPJ
    case(4)
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
    case(5)
        operatorME = MJLDivQ(y,np,lp,jp,n,l,j,bigJ,bigJ)
    case(6)
        operatorME = SigmaPJ(y,np,lp,jp,n,l,j,bigJ)
    case(7)
        operatorME = SigmaPPJ(y,np,lp,jp,n,l,j,bigJ)
    case default
        stop 'The choice of operator j should be from 1 to 7'
    end select
    
end subroutine OperME
