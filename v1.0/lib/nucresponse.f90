
!  ctrlchar = 'c'  count up # of states
!  ctrlchar = 'f'  fill up info on states
!  ctrlchar = 'p'  fill up parent reference states
!  ctrlchar = 'd'  fill up daughter reference states

!================================================

function nucResponse(tau1,tau2,ioption,y)
    use kinds
    use spspace
    use stateinfo
    implicit none

    integer :: tau1, tau2
    integer :: ioption
    real(doublep) :: y

    integer :: j,a,b!,ap,an
    integer :: jmin, jmax
    
    integer :: op1, op2

    REAL(doublep) :: Wigner_3j
    REAL(doublep) :: spOME1,spOME2
    REAL(doublep), dimension (0:6) :: DRME1, DRME2
    REAL(doublep)  :: nucResponse

    if (ioption .eq. 1) then   
        ! Mj Mj
        op1 = 1; op2 = 1
        jmin = 0; jmax = 6
    else if (ioption .eq. 2) then
        ! PhiPPJ PhiPPJ
        op1 = 3; op2 = 3
        jmin = 0; jmax = 6
    else if (ioption .eq. 3) then
        ! PhiTPJ PhiTPJ
        op1 = 4; op2 = 4
        jmin = 2; jmax = 6
    else if (ioption .eq. 4) then
        ! DeltaJ DeltaJ
        op1 = 5; op2 = 5
        jmin = 1; jmax = 5
    else if (ioption .eq. 5) then
        ! SigmaPJ SigmaPJ
        op1 = 6; op2 = 6
        jmin = 1; jmax = 5
    else if (ioption .eq. 6) then
        ! SigmaPPJ SigmaPPJ
        op1 = 7; op2 = 7
        jmin = 1; jmax = 5
    else if (ioption .eq. 7) then
        ! PhiPPJ MJ
        op1 = 3; op2 = 1
        jmin = 0; jmax = 6
    else if (ioption .eq. 8) then
        ! DeltaJ SigmaPJ
        op1 = 5; op2 = 6
        jmin = 1; jmax = 5
    end if

    !print*,'operators: ',op1, op2
!    print*,'debug tiso mtiso',tiso,mtiso
 
    DRME1(0:6) = 0.d0
    DRME2(0:6) = 0.d0

    nucResponse = 0.d0

    Do j = jmin,jmax,2
        Do a = 1, ntotal(1)
            Do b = 1, ntotal(1)
                If (abs(densitymats%rho(j,tau1,a,b)) .ge. 1.0e-9) then
                If (abs(densitymats%rho(j,tau2,a,b)) .ge. 1.0e-9) then

                    ! Operator 1 with tau2 <j| op1,tau1 |j>
                    call OperME(op1,y,nodal(a),lorb(a),jorb(a),nodal(b),lorb(b),jorb(b),j,spOME1)
                    DRME1(j) = DRME1(j) + densitymats%rho(j,tau1,a,b) * spOME1

                    ! Operator 2 with tau2 <j| op2,tau2 |j>
                    call OperME(op2,y,nodal(a),lorb(a),jorb(a),nodal(b),lorb(b),jorb(b),j,spOME2)
                    DRME2(j) = DRME2(j) + densitymats%rho(j,tau2,a,b) * spOME2

                end if
                end if
            end do
        end do
        nucResponse = nucResponse + 1 &
!                         * sqrt(2.d0) *sqrt(2*dble(tau1)+1.0) &
!                         * sqrt(2.d0) *sqrt(2*dble(tau2)+1.0) &
!                         * (-1.0)**((Tiso - Mtiso)/2)*Wigner_3j(Tiso,2*tau1,Tiso,-Mtiso,0,Mtiso) &
!                         * (-1.0)**((Tiso - Mtiso)/2)*Wigner_3j(Tiso,2*tau2,Tiso,-Mtiso,0,Mtiso) &
                         * DRME1(j) * DRME2(j) 
    end do

    nucResponse = nucResponse * 2.0*sqrt( (2*dble(tau1)+1)*(2*dble(tau2)+1) ) &
                  * Wigner_3j(Tiso,2*tau1,Tiso,-Mtiso,0,Mtiso) &
                  * Wigner_3j(Tiso,2*tau2,Tiso,-Mtiso,0,Mtiso)

end function nucResponse
