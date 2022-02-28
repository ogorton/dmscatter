module coefficients

    use main
    use constants

    contains

    subroutine setupcoef

        implicit none
        integer :: num_response_coef
    
        num_response_coef = 15
    
        ! Explicit proton-neutron coefficients
        allocate(eft%xpnc(0)%c(num_response_coef))
        allocate(eft%xpnc(1)%c(num_response_coef))
    
        ! Isospin formalism coefficients
        allocate(eft%isoc(0)%c(num_response_coef))
        allocate(eft%isoc(1)%c(num_response_coef))
    
        eft%xpnc(0)%c = 0.0
        eft%xpnc(1)%c = 0.0
        eft%isoc(0)%c = 0.0
        eft%isoc(1)%c = 0.0
    
    end subroutine setupcoef
    
    subroutine setpncoeffsnonrel(coupling, op, coeffdimless)
        implicit none
        integer, intent(in) :: op
        real(kind=8), intent(in) :: coeffdimless
        character(1), intent(in) :: coupling
    
        if (op<1 .or. op>15) then
            print*,"Warning, EFT operator must be 1 <= i <= 15."
            stop "Invalid EFT coefficient."
        end if
    
        select case(coupling)
            case('p')
                eft%xpnc(0)%c(op) = coeffdimless
            case('n')
                eft%xpnc(1)%c(op) = coeffdimless
            case('s')
                eft%isoc(0)%c(op) = coeffdimless
            case('v')
                eft%isoc(1)%c(op) = coeffdimless
            case default
                print*,'Invalid coupling',coupling
        end select
    
    end subroutine
    
    
    subroutine convertisospinform
        integer :: i
        integer :: num_response_coef
        
        num_response_coef = size(eft%xpnc(0)%c)
    
        print*,'Converting EFT coefficients isospin <=> proton-neutron formalism.'
    
        do i = 1, num_response_coef
            if (eft%xpnc(0)%c(i) .ne. 0 .or. eft%xpnc(1)%c(i).ne.0) then
                eft%isoc(0)%c(i) = (eft%xpnc(0)%c(i) + eft%xpnc(1)%c(i))
                eft%isoc(1)%c(i) = (eft%xpnc(0)%c(i) - eft%xpnc(1)%c(i))
                print*,'s',i,eft%isoc(0)%c(i)
                print*,'v',i,eft%isoc(1)%c(i)
            end if
        end do
        do i = 1, num_response_coef
            if (eft%isoc(0)%c(i) .ne. 0 .or. eft%isoc(1)%c(i).ne.0) then
                eft%xpnc(0)%c(i) = (eft%isoc(0)%c(i) + eft%isoc(1)%c(i))/2
                eft%xpnc(1)%c(i) = (eft%isoc(0)%c(i) - eft%isoc(1)%c(i))/2
                print*,'p',i,eft%xpnc(0)%c(i)
                print*,'n',i,eft%xpnc(1)%c(i)
            end if
        end do    
    
    end subroutine
    
    subroutine normalizecoeffs()
        
        implicit none
        real(kind=8) :: mchi
        integer :: i, j, op, mNdenomOps(7), mNmNdenomOps(2)
    
        print*,'Normalizing EFT coefficients.'
    
        mchi = wimp%mass
    
        mNdenomOps = [3, 5, 9, 10, 11, 13, 14]
        mNmNdenomOps = [6, 15]
    
        do i=0,1
            eft%xpnc(i)%c(:) = eft%xpnc(i)%c(:) * (4.0*mN*mchi)/(mV*mV)
            eft%isoc(i)%c(:) = eft%isoc(i)%c(:) * (4.0*mN*mchi)/(mV*mV)
    
            do j = 1, 2
                op = mNmNdenomOps(j)
                eft%xpnc(i)%c(op) = eft%xpnc(i)%c(op)/(mN*mN)
                eft%isoc(i)%c(op) = eft%isoc(i)%c(op)/(mN*mN)
            enddo
            do j = 1, 7
                op = mNdenomOps(j)
                eft%xpnc(i)%c(op) = eft%xpnc(i)%c(op)/(mN)
                eft%isoc(i)%c(op) = eft%isoc(i)%c(op)/(mN)
            enddo
        enddo
    
    end subroutine normalizecoeffs

end module coefficients
