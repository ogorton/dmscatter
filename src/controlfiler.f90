
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
                call setpncoeffsnonrel(op, coef, nucleon)
                print*,'Set non-relativistic coefficient: op',op,'p/n',nucleon,'c',coef
            else
                call setkeyword(keyword,keyvalue)
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
subroutine setkeyword(keyword, keyvalue)

    use kinds
    use keywords
    use masses
    use velocities
    use constants
    use targetinfo
    use momenta
    use dmparticles
    use quadrature
    use spspace

    implicit none

    character (len=20) :: keyword
    real(doublep) :: keyvalue

    select case (keyword)

    case('vearth')
        ve = keyvalue
        print*,trim(keyword),": Set velocity of earth in galactic frame set to",ve

    case('dmdens')
        rhochi = keyvalue
        print*,trim(keyword),": Set local dark matter density to",rhochi

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
        bfm = (41.467/(45.*(num_n+num_p)**(-1./3) - 25.*(num_n+num_p)**(-2./3)))**0.5 * femtometer
        print*,'bfm updated.'

    case('femtometer')
        femtometer = keyvalue
        print*,keyword,': Set femtometer units to',femtometer
        bfm = (41.467/(45.*(num_n+num_p)**(-1./3) - 25.*(num_n+num_p)**(-2./3)))**0.5 * femtometer
        print*,'bfm updated.'

    case('wimpmass')
        mchi = keyvalue
        print*,keyword,': Set dark matter particle mass to',mchi
        muT = mchi * mtarget * mN / (mchi+mtarget*mN)
        print*,'muT reduced mass updated.'
        vdist_min = q/(2.0*muT)
        print*,'vmin updated.'

    case('vescape')
        vesc = keyvalue
        print*,keyword,': Set escape velocity to', vesc
        vdist_max = vesc

    case('ntarget')
        Nt = keyvalue
        print*,keyword,': Set number of target nuclei to', Nt

    case('weakmscale')
        mV = keyvalue
        print*,keyword,': Set standard-model weak interaction mass scale to',mv

    case('maxwellv0')
        v0 = keyvalue
        print*,keyword,': Set velocity distribution scaling to',v0

    case('mnucleon')
        mn = keyvalue
        print*,keyword,': Set nucleon mass to',mn
        muT = mchi * mtarget * mN / (mchi+mtarget*mN)
        print*,'muT reduced mass udpated.'
        vdist_min = q/(2.0*muT)
        print*,'vmin updated'

    case('dmspin')
        jchi = keyvalue
        print*,keyword,': Set dark matter particle spin to',jchi

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





