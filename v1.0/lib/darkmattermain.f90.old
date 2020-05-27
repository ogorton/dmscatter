
!  ctrlchar = 'c'  count up # of states
!  ctrlchar = 'f'  fill up info on states
!  ctrlchar = 'p'  fill up parent reference states
!  ctrlchar = 'd'  fill up daughter reference states

!================================================

program darkmattermain
    use spspace
    use stateinfo
    implicit none

    interface
       subroutine openresults(resfile)
           use stateinfo
           implicit none
           integer resfile
           character(22):: filename
       end subroutine openresults
    end interface

    integer, parameter :: resfile = 33
    integer :: i,j,a,b,ap,an
    character(23) :: tmpline
    
    integer :: ioption, op1, op2, Mtiso
    integer :: t,tt1,tt2
    REAL(kind=8) :: Wigner_3j
    REAL(kind=8) :: Pi = 3.1415926535897932
    REAL(kind=8) :: y, spOME, spOME1,spOME2
    REAL(kind=8), dimension (0:1) :: IsoME,IsoME1,IsoME2
    REAL(kind=8), dimension (0:6,0:1) :: DRME, SRME, DRME1, SRME1, DRME2, SRME2 
    REAL(kind=8), dimension (0:1,0:1) :: Response
    logical :: success

    call GetSPS

    print*,' '
    print*,' Reading one-body density matrix from Bigstick .res file '
    print*,' '

    nsporb = norb(1)

    call openresults(resfile)
    call setupdensities
    call readheaderv2(resfile)
    call readalldensities(resfile)

    print*,' '
    print*,' Enter the neutron number '
    read(5,*)an

    print*,' '
    print*,' Enter the proton number '
    read(5,*)ap

    Mtiso = ap-an

    success = .false.
    print*,' '
    do while(.not.success)

       print*," Enter options of operator (1-8) "
       print*," (1) W_M "
       print*," (2) W_{\PhiPP}"
       print*," (3) W_{\tilde{\Phi}'} "
       print*," (4) W_{\Delta} "
       print*," (5) W_{\SigmaP} "
       print*," (6) W_{\SigmaPP} "
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
 
    if (ioption .eq. 1) then   
        op1 = 1; op2 = 1
    else if (ioption .eq. 2) then
        op1 = 3; op2 = 3
    else if (ioption .eq. 3) then
        op1 = 4; op2 = 4
    else if (ioption .eq. 4) then
        op1 = 5; op2 = 5
    else if (ioption .eq. 5) then
        op1 = 6; op2 = 6
    else if (ioption .eq. 6) then
        op1 = 7; op2 = 7
    else if (ioption .eq. 7) then
        op1 = 3; op2 = 1
    else if (ioption .eq. 8) then
        op1 = 5; op2 = 6
 
    end if
 
    print*,' '
    print*,' Enter the value of y '
    read(5,*)y

    IsoME(0:1)=0.d0
    DRME(0:6,0:1) = 0.d0
    SRME(0:6,0:1) = 0.d0

    IsoME1(0:1)=0.d0
    DRME1(0:6,0:1) = 0.d0
    SRME1(0:6,0:1) = 0.d0

    IsoME2(0:1)=0.d0
    DRME2(0:6,0:1) = 0.d0
    SRME2(0:6,0:1) = 0.d0

    Response(0:1,0:1) = 0.d0

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (ioption .le. 3) then
        Do j = 0,6,2
            Do t = 0, 1
                Do a = 1, ntotal(1)
                    Do b = 1, ntotal(1)
                        If (abs(densitymats%rho(j,t,a,b)) .ge. 1.0e-9) then
                            call OperME(op1,y,nodal(a),lorb(a),jorb(a),nodal(b),lorb(b),jorb(b),j,spOME)
                            IsoME(t)= spOME * sqrt(2.d0) *sqrt(2*dble(t)+1.0)
                            DRME(j,t) = DRME(j,t) + densitymats%rho(j,t,a,b)*IsoME(t)
                            SRME(j,t) = (-1.0)**((Tiso - Mtiso)/2)*Wigner_3j(Tiso,2*t,Tiso,-Mtiso,0,Mtiso)*DRME(j,t)    
                        end if
                    end do !b
                end do !a
            end do !t
        end do !j

        Do j = 0, 6, 2
            Do tt1 = 0,1
                Do tt2 = 0,1
                    Response(tt1,tt2) = Response(tt1,tt2)+ SRME(j,tt1)*SRME(j,tt2)
                end do !tt2
            end do !tt1
        end do !j
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    else if (ioption .ge. 4 .and. ioption .le. 6) then
        Do j = 1,5,2
            Do t = 0, 1
                Do a = 1, ntotal(1)
                    Do b = 1, ntotal(1)
                        If (abs(densitymats%rho(j,t,a,b)) .ge. 1.0e-9) then
                            call OperME(op1,y,nodal(a),lorb(a),jorb(a),nodal(b),lorb(b),jorb(b),j,spOME)
                            IsoME(t)= spOME * sqrt(2.d0) *sqrt(2*dble(t)+1.0)
                            DRME(j,t) = DRME(jt,t) + densitymats%rho(j,t,a,b)*IsoME(t)
                            SRME(j,t) = (-1.0)**((Tiso - Mtiso)/2)*Wigner_3j(Tiso,2*t,Tiso,-Mtiso,0,Mtiso)*DRME(j,t)
                        end if
                    end do !b
                end do !a
            end do !t
        end do !j

        Do j = 1, 5, 2
            Do tt1 = 0,1
                Do tt2 = 0,1
                    Response(tt1,tt2) = Response(tt1,tt2)+ SRME(j,tt1)*SRME(j,tt2)
              end do !tt2
            end do !tt1
        end do !j
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    else if (ioption .eq. 7) then
        Do j = 0,6,2
            Do t = 0, 1
                Do a = 1, ntotal(1)
                    Do b = 1, ntotal(1)
                        If (abs(densitymats%rho(j,t,a,b)) .ge. 1.0e-9) then
                            call OperME(op1,y,nodal(a),lorb(a),jorb(a),nodal(b),lorb(b),jorb(b),j,spOME1)
                            call OperME(op2,y,nodal(a),lorb(a),jorb(a),nodal(b),lorb(b),jorb(b),j,spOME2)

                            IsoME1(t)= spOME1 * sqrt(2.d0) *sqrt(2*dble(t)+1.0)
                            DRME1(j,t) = DRME1(j,t) + densitymats%rho(j,t,a,b)*IsoME1(t)
                            SRME1(j,t) = (-1.0)**((Tiso - Mtiso)/2)*Wigner_3j(Tiso,2*t,Tiso,-Mtiso,0,Mtiso)*DRME1(j,t)

                            IsoME2(t)= spOME2 * sqrt(2.d0) *sqrt(2*dble(t)+1.0)
                            DRME2(j,t) = DRME2(j,t) + densitymats%rho(j,t,a,b)*IsoME2(t)
                            SRME2(j,t) = (-1.0)**((Tiso - Mtiso)/2)*Wigner_3j(Tiso,2*t,Tiso,-Mtiso,0,Mtiso)*DRME2(j,t)

                        end if
                    end do
                end do
            end do
        end do
        Do j = 0, 6, 2
            Do tt1 = 0,1
                Do tt2 = 0,1
                    Response(tt1,tt2) = Response(tt1,tt2)+ SRME1(j,tt1)*SRME2(j,tt2)
                end do
            end do
        end do
    end if
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Write (*,*) "Non-zero Response for option",ioption 
    Do tt1 = 0,1
        Do tt2 = 0,1
            if (abs(Response(tt1,tt2)) .gt. 1.0d-15) then
                Write(*,"(i4,i4,f25.17)") tt1, tt2, Response(tt1,tt2)
            end if
        end do
    end do
end program  
