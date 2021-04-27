!================================================
module nucresponse
    contains
function nucResponseFun(tau1,tau2,term,y,densmat,Tiso,Mtiso)

    use kinds
    use orbitals
    use parameters
    use wignerfunctions

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
    REAL(doublep)  :: nucResponsefun

    jmin = -1
    jmax = -1

    select case(term)
      case(1)
        ! Mj Mj
        op1 = 1; op2 = 1
        jmin = 0; jmax = 6            
      case(2)
        ! PhiPPJ PhiPPJ
        op1 = 3; op2 = 3
        jmin = 0; jmax = 6
      case(3)
        ! PhiTPJ PhiTPJ
        op1 = 4; op2 = 4
        jmin = 2; jmax = 6
      case(4)
        ! DeltaJ DeltaJ
        op1 = 5; op2 = 5
        jmin = 1; jmax = 5
      case(5)
        ! SigmaPJ SigmaPJ
        op1 = 6; op2 = 6
        jmin = 1; jmax = 5
      case(6)
        ! SigmaPPJ SigmaPPJ
        op1 = 7; op2 = 7
        jmin = 1; jmax = 5
      case(7)
        ! PhiPPJ MJ
        op1 = 3; op2 = 1
        jmin = 0; jmax = 6
      case(8) 
        ! DeltaJ SigmaPJ
        op1 = 5; op2 = 6
        jmin = 1; jmax = 5
      case default
        stop "Invalid term."
    end select        

    nucResponsefun = 0.d0

    do j = jmin,jmax,2
      DRME1 = 0.d0
      DRME2 = 0.d0
      do a = 1, ntotal(1)
        do b = 1, ntotal(1)
          !if (abs(densmat(j,tau1,a,b)) < 1.0e-5 &
          !      .and. abs(densmat(j,tau2,a,b)) < 1.0e-5) cycle

          ! Operator 1 with tau2 <j| op1,tau1 |j>
          if (abs(densmat(j,tau1,a,b)) > 1.0e-5 ) then
              call OperME(op1,y,nodal(a),lorb(a),jorb(a),nodal(b),lorb(b),jorb(b),j,spOME1)
              DRME1 = DRME1 + densmat(j,tau1,a,b) * spOME1 
          end if

          ! Operator 2 with tau2 <j| op2,tau2 |j>
          if (abs(densmat(j,tau2,a,b)) > 1.0e-5 ) then
              call OperME(op2,y,nodal(a),lorb(a),jorb(a),nodal(b),lorb(b),jorb(b),j,spOME2)
              DRME2 = DRME2 + densmat(j,tau2,a,b) * spOME2 
          end if

        end do
      end do
      nucResponsefun = nucResponsefun + DRME1 * DRME2

    end do

    if (.not.pndens) then
        nucResponsefun = nucResponsefun * sqrt(2*dble(tau1)+1.0) * sqrt(2*dble(tau2)+1.0) &
            * tj2i_lookup(Tiso,2*tau1,Tiso,-Mtiso,0,Mtiso) &
            * tj2i_lookup(Tiso,2*tau2,Tiso,-Mtiso,0,Mtiso) &
            * 2.0
    end if

end function nucResponsefun
subroutine test_nucresponse(nuc_target)

    use parameters
    implicit none

    type(nucleus) :: nuc_target
    integer :: ioption, op1, op2, Mtiso
    integer :: tt1,tt2
    REAL(kind=8) :: y !, spOME, spOME1,spOME2
    REAL(kind=8) :: Response
    logical :: success

    success = .false.
    print*,' '
    do while(.not.success)

       print*," Enter options of operator (1-8) "
       ! J=0,2,...
       print*," (1) W_M "
       print*," (2) W_{\PhiPP}"
       ! J=2,4,...
       print*," (3) W_{\tilde{\Phi}'} "
       ! J=1,3,...
       print*," (4) W_{\Delta} "
       print*," (5) W_{\SigmaP} "
       print*," (6) W_{\SigmaPP} "
       ! Two diff operators
       print*," (7) W_{\PhiPP M} "
       print*," (8) W_{\Delta \SigmaP} "

       read(5,*)ioption

       if (ioption .gt. 8 .or. ioption .lt. 1) goto 2
       success = .true.
       goto 1

2      continue
       print*,' options should be 1 to 8 '
1      continue

    end do

    print*,' '
    print*,' Enter the value of y '
    read(5,*)y

    Response = 0.d0

    Write (*,*) "Non-zero Response for option",ioption
    Do tt1 = 0,1
        Do tt2 = 0,1
            Response = nucResponseFun(tt1,tt2,ioption,y,&
                nuc_target%densitymats%rho, nuc_target%groundstate%Tx2, &
                nuc_target%Mt)
            print*, tt1, tt2, Response
        end do
    end do

    print*,'test_nucresponse COMPLETE'

end subroutine test_nucresponse

end module nucresponse
