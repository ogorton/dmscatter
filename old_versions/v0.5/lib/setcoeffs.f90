subroutine setcoeffsnonrel(op, coeffdimless, nucleon)
    use response
    use masses
    implicit none
    integer, intent(in) :: op
    real(kind=8), intent(in) :: coeffdimless
    integer, intent(in) :: nucleon 
    ! nucleon = 0 : proton, 1 : neutron

    cvec(nucleon)%c(op) = coeffdimless 

end subroutine

subroutine opencoeffmatrix(resfile)
    implicit none
    integer, intent(in) :: resfile
    character(22) :: filename
    integer :: ilast
    logical :: success

    success = .false.
    print*,' '
    do while(.not.success)

        print*,' Enter name of coefficient matrix file (.mat) '

        read(5,'(a)')filename
        ilast = index(filename,' ')-1
        open(unit=resfile,file=filename(1:ilast)//'.mat',status='old',err=2)
        success = .true.
        return
2       continue
        print*,filename(1:ilast),'.res does not exist '

    end do    

end subroutine opencoeffmatrix

subroutine readcoeffmatrix(resfile)
    implicit none
    integer, intent(in) :: resfile
    character(20) :: line
    integer :: op
    real(kind=8) :: coef
    integer :: nucleon


    logical :: EOF
    EOF = .false.

    do while (.not. EOF)

        ! Read past comments
        read(resfile,'(a20)',end=111) line

        if (line(1:1).eq.'#' .or. line(1:1).eq.'!') then 
            print*,line
            cycle
        end if

        backspace(resfile)
    
        ! Read in coefficient matrix
        do while(.not.EOF)

            read(resfile,*,end=111) nucleon, op, coef
            call setcoeffsnonrel(op, coef, nucleon)

        end do

        return
111     continue

        EOF = .true.
        return

    end do

end subroutine readcoeffmatrix

!!  If[Op==1,  coeff=(4mN*mchiFORMAL/mV^2)coeffdimless;];
!!  If[Op==2,  coeff=(4mN*mchiFORMAL/mV^2)coeffdimless;];
!!  If[Op==3,  coeff=(4mN*mchiFORMAL/mV^2)coeffdimless/mN;];
!!  If[Op==4,  coeff=(4mN*mchiFORMAL/mV^2)coeffdimless;];
!!  If[Op==5,  coeff=(4mN*mchiFORMAL/mV^2)coeffdimless/mN;];
!!  If[Op==6,  coeff=(4mN*mchiFORMAL/mV^2)coeffdimless/mN^2;];
!!  If[Op==7,  coeff=(4mN*mchiFORMAL/mV^2)coeffdimless;];
!!  If[Op==8,  coeff=(4mN*mchiFORMAL/mV^2)coeffdimless;];
!!  If[Op==9,  coeff=(4mN*mchiFORMAL/mV^2)coeffdimless/mN;];
!!  If[Op==10, coeff=(4mN*mchiFORMAL/mV^2)coeffdimless/mN;];
!!  If[Op==11, coeff=(4mN*mchiFORMAL/mV^2)coeffdimless/mN;];
!!  If[Op==12, coeff=(4mN*mchiFORMAL/mV^2)coeffdimless;];
!!  If[Op==13, coeff=(4mN*mchiFORMAL/mV^2)coeffdimless/mN;];
!!  If[Op==14, coeff=(4mN*mchiFORMAL/mV^2)coeffdimless/mN;];
!!  If[Op==15, coeff=(4mN*mchiFORMAL/mV^2)coeffdimless/mN^2;];
!!  
subroutine normalizecoeffs
    use masses
    use response
    implicit none
    integer :: i, j, op, mNdenomOps(7), mNmNdenomOps(2)

    mNdenomOps = [3, 5, 9, 10, 11, 13, 14]
    mNmNdenomOps = [6, 15]

    do i=0,1
        cvec(i)%c = cvec(i)%c * (4.0*mN*mchi)/(mV*mV)

        do j = 1, 2
            op = mNmNdenomOps(j)
            cvec(i)%c(op) = cvec(i)%c(op)/(mN*mN)
        enddo
        do j = 1, 7
            op = mNdenomOps(j)
            cvec(i)%c(op) = cvec(i)%c(op)/(mN)
        enddo
    enddo

end subroutine normalizecoeffs



