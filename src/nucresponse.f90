!================================================
module nucresponse
    contains
function nucFormFactor(tau1,tau2,term,y,densmat,Tiso,Mtiso)

    ! Computes the nuclear form factor terms for a given 
    ! pn-coupling (iso-coupling) using pn-formalism (iso-formalism)
    ! density matrices.

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
    REAL(doublep)  :: nucFormFactor, isofactor

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

    nucFormFactor = 0.d0

    do j = jmin,jmax,2
      DRME1 = 0.d0
      DRME2 = 0.d0
      do a = 1, ntotal(1)
        do b = 1, ntotal(1)

          ! Operator 1 with tau2 <j| op1,tau1 |j>
          if (abs(densmat(j,tau1,a,b)) > 1.0e-5 ) then
              call OperME(op1,y,nodal(a),lorb(a),jorb(a), &
                                nodal(b),lorb(b),jorb(b), &
                                j,spOME1)
              DRME1 = DRME1 + densmat(j,tau1,a,b) * spOME1 
          end if

          ! Operator 2 with tau2 <j| op2,tau2 |j>
          if (abs(densmat(j,tau2,a,b)) > 1.0e-5 ) then
              call OperME(op2,y,nodal(a),lorb(a),jorb(a),& 
                                nodal(b),lorb(b),jorb(b),&
                                j,spOME2)
              DRME2 = DRME2 + densmat(j,tau2,a,b) * spOME2 
          end if

        end do
      end do
      nucFormFactor = nucFormFactor + DRME1 * DRME2

    end do

    if (.not.pndens) then
        isofactor =  2*sqrt(2*dble(tau1)+1.0) * sqrt(2*dble(tau2)+1.0) &
            * (-1.0)**((Tiso - Mtiso)/2) * tj2i_lookup(Tiso,2*tau1,Tiso,Mtiso,0,-Mtiso) &
            * (-1.0)**((Tiso - Mtiso)/2) * tj2i_lookup(Tiso,2*tau2,Tiso,Mtiso,0,-Mtiso) 
        nucFormFactor = nucFormFactor * isofactor
    end if

end function nucFormFactor
!-----------------------------------------------------------------------------80
function nucFormFactor_transform(tau1,tau2,term,y,densmat,Tiso,Mtiso) result(FF)

    ! Provides the nuclear form factors for a given isospin coupling, computed 
    ! from pn-formalism density matrices, and vice versa. 
    ! For pn -> iso, an additional factor of (1/4) is required.
    ! For iso -> pn, no additional factor is required.

    use kinds
    use orbitals
    use parameters
    use wignerfunctions
    use constants, only: pi

    implicit none

    integer :: tau1, tau2
    integer :: term
    real(doublep) :: y
    real(doublep), allocatable, intent(in) :: densmat(:,:,:,:)

    integer :: Mtiso, Tiso, Jx2iso
    real(doublep) :: FF
    character(len=2) :: coupleto

    write(coupleto,'(I1,I1)')tau1, tau2
    select case(coupleto)
        case('00')
            FF = ( &
                + nucFormFactor(0,0,term,y,densmat,Tiso,Mtiso) &
                + nucFormFactor(0,1,term,y,densmat,Tiso,Mtiso) &
                + nucFormFactor(1,0,term,y,densmat,Tiso,Mtiso) &
                + nucFormFactor(1,1,term,y,densmat,Tiso,Mtiso) )
        case('01')
            FF = ( &
                + nucFormFactor(0,0,term,y,densmat,Tiso,Mtiso) &
                - nucFormFactor(0,1,term,y,densmat,Tiso,Mtiso) &
                + nucFormFactor(1,0,term,y,densmat,Tiso,Mtiso) &
                - nucFormFactor(1,1,term,y,densmat,Tiso,Mtiso) )
        case('10')
            FF = ( &
                + nucFormFactor(0,0,term,y,densmat,Tiso,Mtiso) &
                + nucFormFactor(0,1,term,y,densmat,Tiso,Mtiso) &
                - nucFormFactor(1,0,term,y,densmat,Tiso,Mtiso) &
                - nucFormFactor(1,1,term,y,densmat,Tiso,Mtiso) )
        case('11')
            FF = ( &
                + nucFormFactor(0,0,term,y,densmat,Tiso,Mtiso) &
                - nucFormFactor(0,1,term,y,densmat,Tiso,Mtiso) &
                - nucFormFactor(1,0,term,y,densmat,Tiso,Mtiso) &
                + nucFormFactor(1,1,term,y,densmat,Tiso,Mtiso) )            
        case default
            stop "Invalid coupleto var."
    end select 

end function nucFormFactor_transform
!-----------------------------------------------------------------------------80
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

       print*," Enter your choice of response form factor (1-8) "
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

    Write (*,*) "Response for option",ioption
    if (pndens) then
        print*,"From pn-density:"
        print*,"x,   x',   response"
        Do tt1 = 0,1
          Do tt2 = 0,1
            Response = nucFormFactor(tt1, tt2, ioption, y, &
                nuc_target%densitymats%rho,&
                nuc_target%groundstate%Tx2, &
                nuc_target%Mt )
            print*, tt1, tt2, Response
          end do
        end do
        print*,"tau,   tau',   response"
        do tt1 = 0, 1
          Do tt2 = 0,1
            Response = 0.25d0 * nucFormFactor_transform(tt1, tt2, ioption, y, &
                nuc_target%densitymats%rho,&
                nuc_target%groundstate%Tx2, &
                nuc_target%Mt)
            print*, tt1, tt2, Response
          end do
        end do
    else
        print*,"From iso-density:"
        print*,"x,   x',   response"
        do tt1 = 0, 1
          Do tt2 = 0,1
            Response = 0.25d0 * nucFormFactor_transform(tt1, tt2, ioption, y, &
                nuc_target%densitymats%rho,&
                nuc_target%groundstate%Tx2, &
                nuc_target%Mt)
            print*, tt1, tt2, Response
          end do
        end do

        print*,"tau,   tau',   response"
        Do tt1 = 0,1
          Do tt2 = 0,1
            Response = 0.25d0 * nucFormFactor(tt1, tt2, ioption, y, &
                nuc_target%densitymats%rho,&
                nuc_target%groundstate%Tx2, &
                nuc_target%Mt )
            print*, tt1, tt2, Response
          end do
        end do        
    end if

    print*,'test_nucresponse COMPLETE'

end subroutine test_nucresponse

end module nucresponse
