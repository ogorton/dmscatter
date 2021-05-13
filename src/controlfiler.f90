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
        stop "Missing control file."

    end do    

end subroutine opencontrolfile

subroutine readcontrolfile(resfile, eft, wimp)
    use keywords
    use parameters
    implicit none
    type(particle) :: wimp
    type(eftheory) :: eft
    integer, intent(in) :: resfile
    character(100) :: line
    character(20) :: keyword
    integer :: op
    real(kind=8) :: coef,keyvalue
    character(1) :: coupling
    character(20) :: coefkey

    integer :: i

    logical :: EOF

    print*,"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
    print*,"Reading control file."

    coefkey = "coefnonrel"
    EOF = .false.

    do while (.not. EOF)

        ! Read past comments
        read(resfile,'(a100)',end=111) line

        if (line(1:1).eq.'#' .or. line(1:1).eq.'!') then 
            print*,trim(line)
            cycle
        end if

        backspace(resfile)
    
        ! Read in coefficient matrix
        do while(.not.EOF)
            read(resfile,'(a100)',end=111) line
            backspace(resfile)
            read(resfile,*,err=110)keyword, keyvalue
            if (keyword == coefkey) then
                ! Coefficients have more arguments.
                ! Go back to read them.
                backspace(resfile)
                read(resfile,*,end=111)keyword, op, coupling, coef
                call setpncoeffsnonrel(eft, coupling, op, coef)
            else
                call addKeywordpair(keyword, keyvalue)
                call setkeyword(keyword,keyvalue, wimp)
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

110     continue
        print*,"Problem reading in control command. This command will be ignored:"
        print*,line

    end do

end subroutine readcontrolfile
subroutine setkeyword(keyword, keyvalue, wimp)

    use kinds
    use keywords
    use constants
    use momenta
    use quadrature
    use parameters
    use orbitals, only: bfm

    implicit none

    type(particle) :: wimp

    character (len=20) :: keyword
    real(doublep) :: keyvalue

    select case (keyword)

    case('vearth')
        vearth = keyvalue
        print*,trim(keyword),": Set velocity of earth in galactic frame set to",vearth

    case('dmdens')
        wimp%localdensity = keyvalue
        print*,trim(keyword),": Set local dark matter density to",keyvalue

    case('quadtype')
        quadrature_type = keyvalue
        print*,trim(keyword),': not implemented.'

    case('quadrelerr')
        quadrature_relerr = dble( keyvalue )
        print*,trim(keyword),": Set adaptive quadrature routine to seek relative error:",keyvalue

    case('gev')
        gev = keyvalue
        print*,trim(keyword),': Set GeV units to',gev
        femtometer = 5.0677/GeV
        print*,'femtometer updated'

    case('femtometer')
        femtometer = keyvalue
        print*,trim(keyword),': Set femtometer units to',femtometer

    case('wimpmass')
        wimp%mass = keyvalue
        print*,trim(keyword),': Set dark matter particle mass to',keyvalue

    case('vescape')
        vescape = keyvalue
        print*,trim(keyword),': Set escape velocity to', vescape

    case('weakmscale')
        mV = keyvalue
        print*,trim(keyword),': Set standard-model weak interaction mass scale to',mv

    case('maxwellv0')
        vscale = keyvalue
        print*,trim(keyword),': Set velocity distribution scaling to',vscale

    case('mnucleon')
        mn = keyvalue
        print*,trim(keyword),': Set nucleon mass to',mn

    case('dmspin')
        wimp%j = keyvalue
        print*,trim(keyword),': Set dark matter particle spin to',keyvalue

    case('usemomentum')
        if (keyvalue==1) then
            usemomentum = .true.
            print*,trim(keyword),': Set to use momentum transfer instead& 
                    & of recoil energy'
        else
            usemomentum = .false.
        end if

    case('useenergyfile')
        if (keyvalue==1) then
            useenergyfile = .true.
            print*,trim(keyword),': Set to use recoil energies in file&
                    & instead of linear grid specification.'
        end if

    case('fillnuclearcore')
        if (keyvalue==1) then
            fillcore = .true.
            print*,trim(keyword),': Set to fill density matrix core orbitals.'
        else
            fillcore = .false.
            print*,trim(keyword),': Set to leave denisty matrix core orbitals unfilled.'
        end if
    case('ntscale')
        ntscale = keyvalue
        print*,trim(keyword),': Set target nuclei mass density scaling to',keyvalue

    case('hofrequency')
        bfm = 6.43d0/sqrt(keyvalue) * femtometer
        print*,trim(keyword),': Set harmonic oscillator frequency to',keyvalue,'and b(1/GeV) to',bfm

    case('hoparameter')
        bfm = keyvalue * femtometer
        print*,trim(keyword),': Set harmonic oscillator parameter to',keyvalue,'and b(1/GeV) to',bfm

    case('printdensities')
        if (keyvalue==1) then
            printdens = .true.
            print*,trim(keyword),': Set to print density matrix.'
        else
            printdens = .false.
            print*,trim(keyword),': Set to not print density matrix.'
        end if        

    case('gaussorder')
        gaussorder = int(keyvalue)
        print*,trim(keyword),': Set order of Gauss-Legendre quadrature to',keyvalue

    case default
        print*,'Invalid keyword "',trim(keyword),'". Ignoring.'

    end select

end subroutine setkeyword





