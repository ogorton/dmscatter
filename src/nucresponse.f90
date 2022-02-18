module nucresponse
    contains

    function OperME(opid,y,a,b,bigJ) result(operatorME)
    
        use mjl
        use norm
        use phi
        use sigma
        use orbitals, only: nodal, lorb, jorb
    
        implicit none
    
        character(len=8), intent(in) :: opid
        integer, intent(in) :: a, b
        integer :: np,lp,jp,n,l,j,bigJ
        REAL (kind=8), INTENT(IN)  :: y
        REAL (kind=8) :: operatorME
    
        np = nodal(a)
        lp = lorb(a)
        jp = jorb(a)
        n = nodal(b)
        l = lorb(b)
        j = jorb(b)
    
        select case(opid)
        case('mj')
            operatorME = MJ(y,np,lp,jp,n,l,j,bigJ)
        case('sigmaj')
            operatorME = SigmaJ(y,np,lp,jp,n,l,j,bigJ,bigJ)
        case('phippj')
            operatorME = PhiPPJ(y,np,lp,jp,n,l,j,bigJ)
        case('phitpj')
            operatorME = PhiTPJ(y,np,lp,jp,n,l,j,bigJ)
        case('deltaj')
            operatorME = MJLDivQ(y,np,lp,jp,n,l,j,bigJ,bigJ)
        case('sigmapj')
            operatorME = SigmaPJ(y,np,lp,jp,n,l,j,bigJ)
        case('sigmappj')
            operatorME = SigmaPPJ(y,np,lp,jp,n,l,j,bigJ)
        case default
            stop 'The choice of operator is invalid.'
        end select
    
    end function OperME
    

    function nucFormFactor(tau1, tau2, term, y, densmat, Tiso, Mtiso)
    
        ! Computes the nuclear form factor terms for a given 
        ! pn-coupling (iso-coupling) using pn-formalism (iso-formalism)
        ! density matrices.
    
        
        use orbitals
        use wigner, only: threej_lookup, vector_couple
        use densities, only: maxJt, pndens
    
        implicit none
    
        integer :: tau1, tau2
        integer :: term
        real(kind=8) :: y
        real(kind=8), allocatable, intent(in) :: densmat(:,:,:,:)
        real(kind=8) :: dme1, dme2
    
        integer :: j,a,b!,ap,an
        integer :: jmin, jmax
        
        character(len=8) :: op1, op2
        logical :: same
    
        integer :: Mtiso, Tiso
    
        real(kind=8) :: dmepsilon
        real(kind=8) :: spome
        real(kind=8) :: wfn_me1, wfn_me2
        real(kind=8) :: nucFormFactor, isofactor
    
        jmin = -1
        jmax = -1
        jmax = maxJt
        dmepsilon = 1e-5
    
        same = .true.
        select case(term)
          case(1)
            op1 = 'mj'
            op2 = op1
            jmin = 0
          case(2)
            op1 = 'phippj'
            op2 = op1
            jmin = 0
          case(3)
            op1 = 'phitpj'
            op2 = op1
            jmin = 2
          case(4)
            op1 = 'deltaj'
            op2 = op1
            jmin = 1
          case(5)
            op1 = 'sigmapj'
            op2 = op1
            jmin = 1
          case(6)
            op1 = 'sigmappj'
            op2 = op1
            jmin = 1
          case(7)
            op1 = 'mj'
            op2 = 'phippj'
            jmin = 0
            same = .false.
          case(8) 
            op1 = 'sigmapj'
            op2 = 'deltaj'
            jmin = 1
            same=.false.
          case default
            stop "Invalid term."
        end select        
    
        nucFormFactor = 0.d0
    
        do j = jmin, jmax, 2 ! j of the one-body (transition) operator
            wfn_me1 = 0.d0
            wfn_me2 = 0.d0
            do a = 1, ntotal(1)
                do b = 1, ntotal(1)
                    dme1 = densmat(j,tau1,a,b)
                    dme2 = densmat(j,tau2,a,b)
    
                    if (same) then
                        if (abs(dme1) >= dmepsilon .or. abs(dme2) >= dmepsilon) then
                            spome = OperME(op1, y, a, b, j)
                            wfn_me1 = wfn_me1 + dme1 * spome
                            wfn_me2 = wfn_me2 + dme2 * spome
                        end if
                    else
                        ! Operator 1 with tau2 <j| op1,tau1 |j>
                        if (abs(dme1) >= dmepsilon ) then
                            wfn_me1 = wfn_me1 + dme1 * OperME(op1, y, a, b, j)
                        end if         
    
                        ! Operator 2 with tau2 <j| op2,tau2 |j>
                        if (abs(dme2) >= dmepsilon ) then
                            wfn_me2 = wfn_me2 + dme2 * OperME(op2, y, a, b, j)
                        end if
                    end if
                end do
            end do
            nucFormFactor = nucFormFactor + wfn_me1 * wfn_me2
    
        end do
    
        if (.not.pndens) then

            ! Note on this factor: arguments of the 3j functions
            ! depend on the definition of Mt. If Mt = (Z - N)/2, use:
            !       T  tau   T
            !      -Mt 0     Mt
            ! Else if Mt = (N - Z)/2, use:
            !       T   tau   T
            !       Mt  0    -Mt

            isofactor =  2d0 * sqrt((2d0*dble(tau1)+1d0) * (2d0*dble(tau2)+1d0)) &
                * threej_lookup(Tiso, 2*tau1, Tiso, -Mtiso, 0, Mtiso) &
                * threej_lookup(Tiso, 2*tau2, Tiso, -Mtiso, 0, Mtiso) 
            nucFormFactor = nucFormFactor * isofactor

        end if
    
    end function nucFormFactor
    

    function nucFormFactor_transform(tau1, tau2, term, y, densmat, Tiso, Mtiso) result(FF)
    
        ! Provides the nuclear form factors for a given isospin coupling, computed 
        ! from pn-formalism density matrices, and vice versa. 
        ! For pn -> iso, an additional factor of (1/4) is required.
        ! For iso -> pn, no additional factor is required.
    
        
        use orbitals
    
        implicit none
    
        integer :: tau1, tau2
        integer :: term
        real(kind=8) :: y
        real(kind=8), allocatable, intent(in) :: densmat(:,:,:,:)
    
        integer :: Mtiso, Tiso
        real(kind=8) :: FF
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
    

    subroutine test_nucresponse(nuc_target)
    
        use orbitals, only: bfm
        use types, only: nucleus
        use densities, only: pndens
        implicit none
    
        type(nucleus) :: nuc_target
        integer :: ioption
        integer :: tt1,tt2
        REAL(kind=8) :: y, q !, spOME, spOME1,spOME2
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
        print*,' Enter nuclear recoil q in GeV:'
        read(5,*)q
    
        y = (q * bfm / 2d0) ** 2d0
    
        if (.true.) then
            print*,'q',q
            print*,'bfm',bfm
            print*,'y',y
        end if
    
        Response = 0.d0
    
        write (*,*) "Response for option",ioption
        if (pndens) then
            print*,""
            print*,"From pn-density:"
            print*,""
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
            print*,""         
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
