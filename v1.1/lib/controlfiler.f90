
subroutine opencontrolfile(resfile)
    implicit none
    integer, intent(in) :: resfile
    character(22) :: filename
    integer :: ilast
    logical :: success

    success = .false.
    print*,' '
    do while(.not.success)

        print*,'Enter name of control file (.control):'

        read(5,'(a)')filename
        ilast = index(filename,' ')-1
        open(unit=resfile,file=filename(1:ilast)//'.control',status='old',err=2)
        success = .true.
        return
2       continue
        print*,filename(1:ilast),'.control does not exist '

    end do    

end subroutine opencontrolfile

subroutine readcontrolfile(resfile)
    use keywords
    implicit none
    integer, intent(in) :: resfile
    character(20) :: line
    character(20) :: keyword
    integer :: op
    real(kind=8) :: coef,keyvalue
    integer :: nucleon
    character(20) :: coefkey

    integer :: i

    logical :: EOF

    print*,"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
    print*,"Reading control file."

    coefkey = "coefnonrel"
    EOF = .false.

    print*,'Possible keywords:'
    do i =1,5
        print*,keyword_array(i)
    enddo
    print*,''

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

            read(resfile,*,end=111)keyword, keyvalue
            if (keyword == coefkey) then
                ! Coefficients have more arguments.
                ! Go back to read them.
                backspace(resfile)
                read(resfile,*,end=111)keyword, nucleon, op, coef
                call setpncoeffsnonrel(op, coef, nucleon)
                print*,'Set non-relativistic coefficient: op',op,'p/n',nucleon,'c',coef
            else
                call setkeyword(keyword,keyvalue)
            endif

        end do

        return
111     continue

        EOF = .true.
        return

    end do
  
    print*,'End of control file.'
    print*,''

end subroutine readcontrolfile
subroutine setkeyword(keyword, keyvalue)

    use kinds
    use keywords
    use masses
    use velocities
    use constants
    use dmparticles
    use quadrature

    implicit none

    character (len=20) :: keyword
    real(doublep) :: keyvalue

    if (trim(keyword) == trim(keyword_array(2))) then
        ve = keyvalue
        print*,trim(keyword),": Set velocity of earth in galactic frame set to",ve
    else if (trim(keyword) == trim(keyword_array(3))) then
        rhochi = keyvalue
        print*,trim(keyword),": Set local dark matter density to",rhochi
    else if (trim(keyword) == trim(keyword_array(4))) then
        print*,keyword,'not yet implemented.'
    else if (trim(keyword) == trim(keyword_array(5))) then
        lattice_points = int( keyvalue )
        print*,trim(keyword),": Set number of integral lattice points to",lattice_points
    else
        print*,'Invalid keyword "',trim(keyword),'". Ignoring.'
    endif


end subroutine setkeyword





