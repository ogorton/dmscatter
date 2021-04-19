
module sj2iref
    implicit none

    real(kind=8), allocatable :: sj2i_table(:,:,:,:,:,:)
    real(kind=8), allocatable :: tj2i_table(:,:,:,:,:,:)
    integer :: maxsj2i=12, minsj2i=-2
    integer :: tableJmin, tableJmax
    integer(kind=8) :: sj2i_dim
contains
subroutine sj2itable

! Reads in lookup table for SJ2I symbols
! that were created using the PNISM utility
! sj2i.f90. File must follow the corresponding
! formatting.
    implicit none
    interface 
        FUNCTION Wigner_6j(a,b,c,d,e,f)
            implicit none
            INTEGER, INTENT(IN) :: a,b,c,d,e,f
            REAL(kind=8) :: Wigner_6j            
        end function 
        FUNCTION Wigner_3j(a,b,c,d,e,f)
            implicit none
            INTEGER, INTENT(IN) :: a,b,c,d,e,f
            REAL(kind=8) :: Wigner_3j
        end function        
    end interface
    logical :: timeit
    integer (kind=8) :: ti, tf, clock_rate
    real :: rn
    integer :: i,j,k,l,m,n, a,b

    tableJmin = maxsj2i
    tableJmax = 0
    timeit = .true.
    call system_clock(count_rate = clock_rate)

    if (timeit) call system_clock(count = ti)
    write(6,'(a)')repeat("-",80)
    print*,'Generating Wigner six-J look-up table'
    write(6,'(a)')repeat("-",80)
    print*,'Table max 2J:',maxsj2i

    a = minsj2i
    b = maxsj2i

    allocate(sj2i_table(a:b,a:b,a:b,a:b,a:b,a:b))
    sj2i_table=0.0

    allocate(tj2i_table(a:b,a:b,a:b,a:b,a:b,a:b))
    tj2i_table=0.0    

!$OMP PARALLEL DO PRIVATE(j,k,l,m,n)
    do n = minsj2i, maxsj2i
        print*,'Jx2=',n
        do m = minsj2i, maxsj2i
            do l = minsj2i, maxsj2i
                do k = minsj2i, maxsj2i
                    do j = minsj2i, maxsj2i
                        do i = minsj2i, maxsj2i
                           sj2i_table(i,j,k,l,m,n)=Wigner_6j(I,J,K,L,M,N)
                           tj2i_table(i,j,k,l,m,n)=Wigner_3j(I,J,K,L,M,N)
                        enddo
                    enddo
                enddo
            enddo
        enddo
    enddo
!$OMP end parallel do

    print*,'Tables have been saved to memory. Mem. used (MB):',real(sizeof(rn)*2*(maxsj2i+1)**6,4)*10.**(-6.)
    if (timeit) then
        call system_clock(count = tf)
        print*,'Time:',real((tf-ti))/real(clock_rate)
    endif

end subroutine sj2itable

function sj2i_lookup(I,J,K,L,M,N) result(sj2i)
    implicit none
    real(kind=8) sj2i
    interface
        FUNCTION Wigner_6j(a,b,c,d,e,f)
            implicit none
            INTEGER, INTENT(IN) :: a,b,c,d,e,f
            REAL(kind=8) :: Wigner_6j
        end function
    end interface
    integer, intent(in) :: I,J,K,L,M,N
    integer :: maxJ, minJ

    maxJ = max(I,J,K,L,M,N)
    minJ = min(I,J,K,L,M,N)
    tableJmin = min(tableJmin, minJ)
    tableJmax = max(tableJmax, maxJ)    

    if (maxJ > maxsj2i .or. minJ<0) then
        sj2i = real(Wigner_6j(I,J,K,L,M,N))
    else
        sj2i = sj2i_table(i,j,k,l,m,n)
    endif

    return

end function

function tj2i_lookup(I,J,K,L,M,N) result(tj2i)
    implicit none
    real(kind=8) tj2i
    interface
        FUNCTION Wigner_3j(a,b,c,d,e,f)
            implicit none
            INTEGER, INTENT(IN) :: a,b,c,d,e,f
            REAL(kind=8) :: Wigner_3j
        end function
    end interface
    integer, intent(in) :: I,J,K,L,M,N
    integer :: maxJ

    maxJ = max(I,J,K,L,M,N)

    if (maxJ > maxsj2i) then
        tj2i = real(Wigner_3j(I,J,K,L,M,N))
    else
        tj2i = tj2i_table(i,j,k,l,m,n)
    endif

    return

end function

end module sj2iref


