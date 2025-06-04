module orbitals

    use main
    implicit none
    ! Harmonic oscillator parameter
    real(kind=8) :: bfm = 0.0
    !------------------SINGLE-PARTICLE STATES---------------------------   

    integer norb(2)                ! # of single-particle j-orbits
                                   ! 1 = p,  2 = n
    integer nsps(2)                ! # of single-particle m-states
                                   ! 1 = p, 2 = n
    integer ncore(2)               ! # of single-particle j-orbits in the core
                                   ! 1 = p,  2 = n
    integer ntotal(2)
    integer,allocatable  :: jorb(:)     ! 2 x J of orbit (note: p,n different)
    integer,allocatable  :: lorb(:)
    integer,allocatable  :: nrorb(:)    ! radial quantum number N of j-orbit
    integer,allocatable  :: torb(:)     ! 2 x Tz of orbit ( p/n = +/-1 )
    integer,allocatable  :: nprincipal(:)
    integer,allocatable  :: nodal(:)
    logical spinless

    type orbit
        integer :: nradial
        integer :: l
        integer :: jx2
        integer :: nodal
    end type orbit    

    type(orbit) :: cores(11)
    ! nodal = nradial + 1:
    data cores(1)%nradial/0/, cores(1)%l/0/, cores(1)%jx2/1/, cores(1)%nodal/1/
    data cores(2)%nradial/0/, cores(2)%l/1/, cores(2)%jx2/3/, cores(2)%nodal/1/
    data cores(3)%nradial/0/, cores(3)%l/1/, cores(3)%jx2/1/, cores(3)%nodal/1/
    data cores(4)%nradial/1/, cores(4)%l/0/, cores(4)%jx2/1/, cores(4)%nodal/2/
    data cores(5)%nradial/0/, cores(5)%l/2/, cores(5)%jx2/5/, cores(5)%nodal/1/
    data cores(6)%nradial/0/, cores(6)%l/2/, cores(6)%jx2/3/, cores(6)%nodal/1/
    data cores(7)%nradial/0/, cores(7)%l/3/, cores(7)%jx2/7/, cores(7)%nodal/1/
    data cores(8)%nradial/1/, cores(8)%l/1/, cores(8)%jx2/3/, cores(8)%nodal/2/
    data cores(9)%nradial/0/, cores(9)%l/3/, cores(9)%jx2/5/, cores(9)%nodal/1/
    data cores(10)%nradial/1/, cores(10)%l/1/, cores(10)%jx2/1/, cores(10)%nodal/2/
    data cores(11)%nradial/0/, cores(11)%l/4/, cores(11)%jx2/9/, cores(11)%nodal/1/
    contains

    subroutine getorbitals(resfile, A)
        implicit none
        integer :: resfile
        integer :: A ! Z + N of target
    
        integer :: Acore
        integer :: i
        logical :: done
        character(len=100) :: tmpline
        integer :: maxorbits
        integer, allocatable :: orbiti(:), qnN(:), qnL(:), qn2J(:)
    
        integer :: Zval, Nval
        integer :: maxvalence, ncoreorb
    
        maxorbits = 1000
        allocate(orbiti(maxorbits))
        allocate(qnN(maxorbits))
        allocate(qnL(maxorbits))
        allocate(qn2J(maxorbits))
        orbiti = -1
    
        rewind(resfile)
        read(resfile,*)
        read(resfile,*)Zval, Nval
        print*,'Valence Z, N:'
        print*,Zval, Nval
    
        rewind(resfile)
        done = .false.
        print*,'Looking for orbitals.'
        do while (.not. done)
            read(resfile,'(a)')tmpline
            if (tmpline(1:5)=='ORBIT') then
                print*,'Found single particle state quantum numbers.'
                do i = 1, maxorbits
                                         ! orbit, n radial, l, 2J
                    read(resfile,*,err=1)orbiti(i), qnN(i), qnL(i), qn2J(i)
                end do
            end if
            cycle
    1       done = .true.
        end do
        print*,'Valence orbitals:'
        print*,'ORBIT      N      L      2 x J'
        norb(1) = 0
        maxvalence = 0
        do i = 1, maxorbits
            if (orbiti(i)<0) exit
            norb(1) = norb(1) + 1
            maxvalence = maxvalence + (qn2J(i) + 1)
            print*,orbiti(i), qnN(i), qnL(i), qn2J(i)
        end do 
        norb(2) = norb(1)
        print*,'Number of valence orbitals:'
        print*,norb(1)
        print*,'Max. valence occupancy:'
        print*,maxvalence
    
        print*,"Corrected valence particles Z, N:"
        if (zval < 0) zval = maxvalence + zval
        if (nval < 0) nval = maxvalence + nval
        print*,Zval, Nval
    
        nuc_target%Zval = Zval
        nuc_target%Nval = Nval

        call infercore(A, zval, nval, Acore, ncoreorb)
    
        print*,'Number of core orbitals:'
        print*,ncoreorb
    
        ntotal(1) = ncoreorb + norb(1)
        ntotal(2) = ncoreorb + norb(2)
    
        allocate(nrorb(ncoreorb + norb(1)))
        allocate(jorb(ncoreorb + norb(1)))
        allocate(lorb(ncoreorb + norb(1)))
        allocate(nodal(ncoreorb + norb(1)))
    
        ! Assign valence orbital quantum numbers
        nrorb(1:norb(1)) = qnN(1:norb(1))
        lorb(1:norb(1)) = qnL(1:norb(1))
        jorb(1:norb(1)) = qn2J(1:norb(1))
        nodal(1:norb(1)) = qnN(1:norb(1)) + 1
    
        ! Assign core orbital quantum numbers
        nrorb(norb(1)+1:ntotal(1)) = cores(1:ncoreorb)%nradial
        lorb(norb(1)+1:ntotal(1)) = cores(1:ncoreorb)%l
        jorb(norb(1)+1:ntotal(1)) = cores(1:ncoreorb)%jx2
        nodal(norb(1)+1:ntotal(1)) = cores(1:ncoreorb)%nodal
    
        print*,'Core orbitals:'
        print*,'ORBIT      N      L      2 x J'    
        do i = norb(1) + 1, ntotal(1)
            print*,i,nrorb(i),lorb(i),jorb(i)
        end do
    
    end subroutine getorbitals
    
    subroutine infercore(A, zval, nval, Acore, ncoreorb)
    
        integer :: A, zval, nval
        integer :: Acore, ncoreorb
    
        Acore = A - zval - nval
        print*,'Inert core inferred from valence particles:'
        select case(Acore)
            case(0)
                print*,'No core.'
                ncoreorb = 0
            case(4)
                print*,'He-4 core.'
                ncoreorb = 1
            case(16)
                print*,'O-16 core.'
                ncoreorb = 3
            case(40)
                print*,'Ca-40 core.'
                ncoreorb = 6
            case(56)
                print*,'Ni-56 core.'
                ncoreorb = 7
            case(100)
                print*,'Sn-100 core.'
                ncoreorb = 11
            case default
                print*,'Non-standard core: A-core=',Acore
                ncoreorb = -1
        end select
    
    end subroutine infercore
end module orbitals
