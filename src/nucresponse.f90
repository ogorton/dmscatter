!================================================

function nucResponse(tau1,tau2,term,y,densmat,Tiso,Mtiso)

    use kinds
    use spspace
    use parameters
    use sj2iref

    implicit none

    integer :: tau1, tau2
    integer :: term
    real(doublep) :: y
    real(doublep), allocatable, intent(in) :: densmat(:,:,:,:)

    integer :: j,a,b!,ap,an
    integer :: jmin, jmax
    
    integer :: op1, op2

    integer :: Mtiso, Tiso

    REAL(doublep) :: spOME1,spOME2
    REAL(doublep) :: DRME1, DRME2
    REAL(doublep)  :: nucResponse

    jmin = -1
    jmax = -1

    if (term .eq. 1) then   
        ! Mj Mj
        op1 = 1; op2 = 1
        jmin = 0; jmax = 6
    else if (term .eq. 2) then
        ! PhiPPJ PhiPPJ
        op1 = 3; op2 = 3
        jmin = 0; jmax = 6
    else if (term .eq. 3) then
        ! PhiTPJ PhiTPJ
        op1 = 4; op2 = 4
        jmin = 2; jmax = 6
    else if (term .eq. 4) then
        ! DeltaJ DeltaJ
        op1 = 5; op2 = 5
        jmin = 1; jmax = 5
    else if (term .eq. 5) then
        ! SigmaPJ SigmaPJ
        op1 = 6; op2 = 6
        jmin = 1; jmax = 5
    else if (term .eq. 6) then
        ! SigmaPPJ SigmaPPJ
        op1 = 7; op2 = 7
        jmin = 1; jmax = 5
    else if (term .eq. 7) then
        ! PhiPPJ MJ
        op1 = 3; op2 = 1
        jmin = 0; jmax = 6
    else if (term .eq. 8) then
        ! DeltaJ SigmaPJ
        op1 = 5; op2 = 6
        jmin = 1; jmax = 5
    end if

    nucResponse = 0.d0

    do j = jmin,jmax,2
       DRME1 = 0.d0
       DRME2 = 0.d0
      do a = 1, ntotal(1)
        do b = 1, ntotal(1)
          if (abs(densmat(j,tau1,a,b)) .ge. 1.0e-5 &
                .or. abs(densmat(j,tau2,a,b)) .ge. 1.0e-5) then

            ! Operator 1 with tau2 <j| op1,tau1 |j>
            call OperME(op1,y,nodal(a),lorb(a),jorb(a),nodal(b),lorb(b),&
                    jorb(b),j,spOME1)
            DRME1 = DRME1 + densmat(j,tau1,a,b) * spOME1 

            ! Operator 2 with tau2 <j| op2,tau2 |j>
            call OperME(op2,y,nodal(a),lorb(a),jorb(a),nodal(b),lorb(b),&
                    jorb(b),j,spOME2)
            DRME2 = DRME2 + densmat(j,tau2,a,b) * spOME2 

          end if
        end do
      end do
      nucResponse = nucResponse + DRME1 * DRME2

    end do

    if (.not.pndens) then
        nucResponse = nucResponse * sqrt(2*dble(tau1)+1.0) * sqrt(2*dble(tau2)+1.0) &
            * tj2i_lookup(Tiso,2*tau1,Tiso,-Mtiso,0,Mtiso) &
            * tj2i_lookup(Tiso,2*tau2,Tiso,-Mtiso,0,Mtiso) &
            * 2.0
    end if

end function nucResponse
