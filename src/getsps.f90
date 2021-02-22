subroutine getorbitals(resfile, A)
    use spspace
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

    maxorbits = 100
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

    call infercore(A, zval, nval, maxvalence, Acore, ncoreorb)

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

subroutine infercore(A, zval, nval, maxvalence, Acore, ncoreorb)

    integer :: A, zval, nval, maxvalence
    integer :: Acore, ncoreorb

    if (zval < 0) zval = maxvalence + zval
    if (nval < 0) nval = maxvalence + nval

    Acore = A - zval - nval
    print*,'Inert core inferred from valence particles:'
    select case(Acore)
        case(0)
            print*,'No core.'
            ncoreorb = 0
        case(16)
            print*,'O-16 core.'
            ncoreorb = 3
        case(40)
            print*,'Ca-40 core.'
            ncoreorb = 6
        case(44)
            print*,'Ni-56 core.'
            ncoreorb = 7
        case(100)
            print*,'Sn-100 core.'
            ncoreorb = 11
        case default
            print*,'Non-standard core: A=',Acore
            ncoreorb = -1
    end select

end subroutine infercore

