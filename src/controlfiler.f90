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

subroutine readcontrolfile(resfile, eft, wimp, detector_t)
    use keywords
    use parameters
    implicit none
    type(particle) :: wimp
    type(detector) :: detector_t
    type(eftheory) :: eft
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
    do i =1,14
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
                call setpncoeffsnonrel(eft, op, coef, nucleon)
                print*,'Set non-relativistic coefficient: op',op,'p/n',nucleon,'c',coef
            else
                call setkeyword(keyword,keyvalue, wimp, detector_t)
            endif

        end do

        print*,'End of control file.'
        print*,''

        return
111     continue

        EOF = .true.
        print*,'End of control file.'
        print*,''
        return

    end do

end subroutine readcontrolfile
subroutine setkeyword(keyword, keyvalue, wimp, detector_t)

    use kinds
    use keywords
    use constants
    use momenta
    use quadrature
    use spspace
    use parameters

    implicit none

    type(particle) :: wimp
    type(detector) :: detector_t

    character (len=20) :: keyword
    real(doublep) :: keyvalue

    select case (keyword)

    case('vearth')
        vdist_t%vearth = keyvalue
        print*,trim(keyword),": Set velocity of earth in galactic frame set to",vdist_t%vearth

    case('dmdens')
        wimp%localdensity = keyvalue
        print*,trim(keyword),": Set local dark matter density to",keyvalue

    case('quadtype')
        print*,keyword,': not implemented.'

    case('intpoints')
        lattice_points = int( keyvalue )
        print*,trim(keyword),": Set number of integral lattice points to",lattice_points

    case('gev')
        gev = keyvalue
        print*,keyword,': Set GeV units to',gev
        femtometer = 5.0677/GeV
        print*,'femtometer updated'

    case('femtometer')
        femtometer = keyvalue
        print*,keyword,': Set femtometer units to',femtometer

    case('wimpmass')
        wimp%mass = keyvalue
        print*,keyword,': Set dark matter particle mass to',keyvalue

    case('vescape')
        vdist_t%vescape = keyvalue
        print*,keyword,': Set escape velocity to', vdist_t%vescape
        vdist_max = vdist_t%vescape

    case('ntarget')
        detector_t%nt = keyvalue
        print*,keyword,': Set number of target nuclei to', keyvalue

    case('weakmscale')
        mV = keyvalue
        print*,keyword,': Set standard-model weak interaction mass scale to',mv

    case('maxwellv0')
        vdist_t%vscale = keyvalue
        print*,keyword,': Set velocity distribution scaling to',vdist_t%vscale

    case('mnucleon')
        mn = keyvalue
        print*,keyword,': Set nucleon mass to',mn

    case('dmspin')
        wimp%j = keyvalue
        print*,keyword,': Set dark matter particle spin to',keyvalue

    case('usemomentum')
        if (keyvalue==1) then
            usemomentum = .true.
            print*,'Set usemomentum to true; using momentum transfer instead& 
                    & of recoil energy'
        end if

    case('useenergyfile')
        if (keyvalue==1) then
            useenergyfile = .true.
            print*,'Set useenergyfile to true; using recoil energies in file&
                    & instead of linear grid specification.'
        end if

    case default
        print*,'Invalid keyword "',trim(keyword),'". Ignoring.'

    end select

end subroutine setkeyword





