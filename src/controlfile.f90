module keywords

    implicit none
    type keywordpair
        character(len=20) :: key
        real(kind=8) :: val
    end type keywordpair
    integer, parameter :: maxkeywords = 100
    integer :: numkeywords = 0
    type(keywordpair) :: keywordpairs(maxkeywords)
    integer :: controlfileid

    contains

    subroutine controlfile(eft, wimp)

        use types, only: eftheory, particle
        implicit none
        type(eftheory) :: eft
        type(particle) :: wimp

        call setupcoef(eft)
        call opencontrolfile()
        call readcontrolfile(eft, wimp)
        call convertisospinform(eft)
        call normalizecoeffs(eft, wimp)        

    end subroutine controlfile

    subroutine addKeywordpair(key, val)
        implicit none
        character(len=20) :: key
        real(kind=8) :: val
        numkeywords = numkeywords + 1
        keywordpairs(numkeywords)%key = key
        keywordpairs(numkeywords)%val = val
    end subroutine addKeywordpair
        
    subroutine opencontrolfile()
        implicit none
        character(22) :: filename
        integer :: ilast
        logical :: success
    
        success = .false.
        print*,' '
        do while(.not.success)
    
            print*,'Enter name of control file (.control):'
    
            read(5,'(a)')filename
            ilast = index(filename,' ')-1
            open(newunit=controlfileid,file=filename(1:ilast)//'.control',status='old',err=2)
            success = .true.
            return
    2       continue
            print*,filename(1:ilast),'.control does not exist '
    
        end do    
    
    end subroutine opencontrolfile
    
    subroutine readcontrolfile(eft, wimp)
        
        use types, only: particle, eftheory
        implicit none
        type(particle) :: wimp
        type(eftheory) :: eft
        character(100) :: line
        character(20) :: keyword
        integer :: op
        real(kind=8) :: coef,keyvalue
        character(1) :: coupling
        character(20) :: coefkey
    
        print*,"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
        print*,"Reading control file."
    
        coefkey = "coefnonrel"
    
        do while(.true.)
    
            read(controlfileid,'(a100)',end=111) line
    
            if (line(1:1).eq.'#' .or. line(1:1).eq.'!') then
                print*,'Comment: ',trim(line)
                cycle
            end if            
    
            backspace(controlfileid)
            read(controlfileid,*,err=110)keyword, keyvalue
            if (keyword == coefkey) then
                ! Coefficients have more arguments.
                ! Go back to read them.
                backspace(controlfileid)
                read(controlfileid,*,end=111)keyword, op, coupling, coef
                call setpncoeffsnonrel(eft, coupling, op, coef)
            else
                call addKeywordpair(keyword, keyvalue)
                call setkeyword(keyword,keyvalue, wimp)
            endif
    
        end do
    
        print*,'End of control file.'
        print*,''
    
        return
    111 continue
    
        print*,'End of control file.'
        print*,''
        return
    
    110 continue
        print*,"Problem reading in control command. This command will be ignored:"
        print*,line
    
    end subroutine readcontrolfile
    
    !=============================================
    subroutine setkeyword(keyword, keyvalue, wimp)
    
        use kinds
        
        use constants
        use settings 
        use quadrature
        use types, only: particle
        use wigner, only: tablemin2j, tablemax2j
        use orbitals, only: bfm
    
        implicit none
    
        type(particle) :: wimp
    
        character (len=20) :: keyword
        real(dp) :: keyvalue
    
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
            print*,trim(keyword),': Set experimental mass scaling to',keyvalue
    
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
    
        case('sj2tablemin')
            tablemin2j = int(keyvalue)
            print*,trim(keyword),': Set minimum Wigner 6-J table value to Jx2 =',keyvalue
    
        case('sj2tablemax')
            tablemax2j = int(keyvalue)
            print*,trim(keyword),': Set mmaximum Wigner 6-J table value to Jx2 =',keyvalue
    
        case default
            print*,''
            print*,'WARNING: Invalid keyword "',trim(keyword),'". Ignoring.'
            print*,''
    
            stop "Invalid keyword in control file."
    
        end select
    
    end subroutine setkeyword

end module keywords