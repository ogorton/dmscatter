subroutine setcoeffsnonrel(op, coeffdimless, nucleon)
    use response
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
            print*,nucleon, op, coef
            call setcoeffsnonrel(op, coef, nucleon)

        end do

        return
111     continue

        EOF = .true.
        return

    end do

end subroutine readcoeffmatrix
