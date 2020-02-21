subroutine test_nucresponse(Mtiso)

    implicit none

    interface
       function nucresponse(tau1, tau2, ioption, y, Mtiso)
           implicit none
           REAL(kind=8) :: nucResponse
           integer :: ioption
           integer :: tau1, tau2
           real(kind=8) :: y
           integer :: Mtiso
       end function nucresponse
    end interface

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
            Response = 0.d0
            Response = nucResponse(tt1, tt2, ioption, y, Mtiso)
            print*, tt1, tt2, Response
        end do
    end do

    print*,'test_nucresponse COMPLETE'

end subroutine test_nucresponse
