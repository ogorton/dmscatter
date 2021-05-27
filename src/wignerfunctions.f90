
module wignerfunctions
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
    print '("min Jx2:  ",I3)',minsj2i
    print '("max Jx2:  ",I3)',maxsj2i
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
!    interface
!        FUNCTION Wigner_6j(a,b,c,d,e,f)
!            implicit none
!            INTEGER, INTENT(IN) :: a,b,c,d,e,f
!            REAL(kind=8) :: Wigner_6j
!        end function
!    end interface
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
    !interface
    !    FUNCTION Wigner_3j(a,b,c,d,e,f)
    !        implicit none
    !        INTEGER, INTENT(IN) :: a,b,c,d,e,f
    !        REAL(kind=8) :: Wigner_3j
    !    end function
    !end interface
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



FUNCTION CG(j1, m1, j2, m2, J, M)

  IMPLICIT NONE

!  INTERFACE
!     FUNCTION FacLOG(N)
!
!       IMPLICIT NONE
!       INTEGER, INTENT(IN) :: N  
!       REAL(kind=8) ::  FacLOG
!     END FUNCTION FacLOG
!  END INTERFACE
  
  INTEGER, INTENT(IN) :: j1, j2, m1, m2, J, M
  REAL(kind=8) :: CG

  INTEGER :: H_Min, H_Max, H
  REAL(kind=8) :: CG_1, CG_2

  REAL(kind=8) :: fact1, fact2, fact3, fact4, fact5, fact6, fact7, fact8
  REAL(kind=8) :: fact9, fact10

  REAL(kind=8) :: fact11
  
  IF (  &
       & (m1 + m2   .NE. M )        .OR. &
       & (ABS(m1)   .GT. j1)        .OR. &
       & (ABS(m2)   .GT. j2)        .OR. &
       & (ABS(M)    .GT. J )        .OR. &  
       & (MOD(j1+m1,2) .NE. 0) .OR. &
       & (MOD(j2+m2,2) .NE. 0) .OR. &  
       & (MOD(J + M,2 ).NE. 0) .OR. &
       & (j1 + J  - j2 .LT. 0)       .OR. &
       & (J  - j1 + j2 .LT. 0)       .OR. &
       & (j1 + j2 - J  .LT. 0)       .OR. &
       & (j1 .LT. 0) .OR. (j2 .LT. 0) .OR. (J .LT. 0))THEN
     CG = 0 
     RETURN
  END IF

  H_Max = MIN(J-j1+j2, J+m1+m2, J+j2+m1) / 2
  H_Min = MAX(0, m1-j1, -j1+j2+m1+m2) / 2

  IF (H_Min .GT. H_Max) STOP "Factorial: Hmin Hmax error"

  fact1 = FacLOG((j1 + J - j2) / 2) 
  fact2 = FacLOG((J - j1 + j2) / 2) 
  fact3 = FacLOG((j1 + j2 - J) / 2) 
  fact4 = FacLOG((J + m1 + m2) / 2) 
  fact5 = FacLOG((J - m1 - m2) / 2) 
  fact6 = FacLOG((J + j1 + j2) / 2 + 1) 
  fact7 = FacLOG((j1 - m1)     / 2) 
  fact8 = FacLOG((j1 + m1)     / 2) 
  fact9 = FacLOG((j2 - m2)     / 2) 
  fact10= FacLOG((j2 + m2)     / 2) 

 !  IF (   fact1 .LE. 0.0 .OR. &
 !      & fact2 .LE. 0.0 .OR. &
 !      & fact3 .LE. 0.0 .OR. &
 !      & fact4 .LE. 0.0 .OR. &
 !      & fact5 .LE. 0.0 .OR. &
 !      & fact6 .LE. 0.0 .OR. &
 !      & fact7 .LE. 0.0 .OR. &
 !      & fact8 .LE. 0.0 .OR. &
 !      & fact9 .LE. 0.0 .OR. &
 !      & fact10 .LE. 0.0 ) THEN
 !    WRITE(0,*) j1, j2, J, m1, m2
 !    WRITE(0,*) fact1, fact2, fact3, fact4, fact5
 !    WRITE(0,*) fact6, fact7, fact8, fact9, fact10
 !    STOP
  !END IF

 CG_1 = fact1 + fact2 + fact3 + fact4 + fact5 &
       & - &
       &(fact6 + fact7 + fact8 + fact9 + fact10)

 CG_2 = 0.0D0

  DO H = H_Min, H_Max

     fact1 = FacLOG((J + j2 + m1) / 2 - H) 
     fact2 = FacLOG((j1 - m1) / 2 + H) 
     fact3 = FacLOG((J - j1 + j2) / 2 - H) 
     fact4 = FacLOG((J + m1 + m2) / 2 - H) 
     fact5 = FacLOG(H) 
     fact6 = FacLOG((j1 - j2 - m1 - m2) / 2 + H) 

     fact11 = fact1 + fact2 - (fact3 + fact4 + fact5 + fact6)

     CG_2 = CG_2 + &
          & (-1)**(MOD(H+(j2+m2)/2,2)) * SQRT(DBLE(J+1)) * EXP(fact11)

  END DO

  CG = exp(0.50D0 * CG_1) * CG_2

  RETURN

END FUNCTION CG

FUNCTION Wigner_3j (j1, j2, J, m1, m2, M)

  IMPLICIT NONE
  
  ! Wigner 3j-Symbol 
  
  ! j1/2 j2/2 J/2 !
  ! m1/2 m2/2 M/2 !

!  INTERFACE     
!     FUNCTION CG(j1, m1, j2, m2, J, M)
!
!       IMPLICIT NONE
!       INTEGER, INTENT(IN) :: j1, j2, m1, m2, J, M
!       REAL(kind=8) :: CG
!     END FUNCTION CG
!  END INTERFACE
!
  INTEGER, INTENT(IN) :: j1, j2, J, m1, m2, M

  REAL (kind=8) :: Wigner_3j

  Wigner_3j = (-1)**((j1-j2-M)/2) / SQRT(DBLE(J)+1.0D0) * CG(j1,m1,j2,m2,J,-M)

  RETURN

END FUNCTION Wigner_3j


FUNCTION FacLOG(N)

  IMPLICIT NONE

  ! return LOG(N!), not N!

  INTEGER, INTENT(IN) :: N
  
  REAL(kind=8) ::  FacLOG

  INTEGER :: i

  FacLOG = 1.0D0

  IF ( N .EQ. 0) THEN
     FacLOG = 0.0D0
     RETURN
  END IF
 
  IF ( N .LT. 0) THEN
     PRINT *, N
     STOP "FacLOG: Impossible to calculate N! N is negative"
  END IF
  
  DO i = 1, N
     FacLOG = FacLOG * DBLE(i)
  END DO

  FacLOG = LOG(FacLOG)
  
  RETURN

END FUNCTION FacLOG

FUNCTION DBLEFacLOG(N)

  IMPLICIT NONE
  ! returns LOG(N!!), not N!!
  INTEGER, INTENT(IN) :: N
  REAL(kind=8) ::  DBLEFacLOG

  INTEGER :: i

  DBLEFacLOG = 1.0D0

  IF ( N .EQ. 0) THEN
     DBLEFacLOG = 0.0D0
     RETURN
  END IF
 
  IF ( N .LT. 0) THEN
     PRINT *, N
     STOP "FacLOG: Impossible to calculate N! N is negative"
  END IF

  IF (MOD(N, 2) .EQ. 0 ) THEN ! N = even
     DO i = 2, N, 2
        DBLEFacLOG = DBLEFacLOG * DBLE(i)
     END DO
  ELSE ! N = odd
     DO i = 1, N, 2
        DBLEFacLOG = DBLEFacLOG * DBLE(i)
     END DO
  END IF

  DBLEFacLOG = LOG(DBLEFacLOG)
  
  RETURN

END FUNCTION DBLEFacLOG

FUNCTION Triangle(a,b,c)

  IMPLICIT NONE

  ! \Delta(a/2, b/2, c/2)
  ! \Delta(a/2, b/2, c/2) = sqrt( (a/2+b/2-c/2)! (a/2-b/2+c/2)! (-a/2+b/2+c/2)! / (a/2+b/2+c/2+1)! )
  ! if a/2,b/2,c/2 do not satisfy the triangle inequality, it returns zero.

!  INTERFACE
!     FUNCTION FacLOG(N)
!
!       IMPLICIT NONE
!       INTEGER, INTENT(IN) :: N
!       REAL(kind=8) ::  FacLOG
!     END FUNCTION FacLOG
!  END INTERFACE

  INTEGER, INTENT(IN) :: a,b,c

  REAL(kind=8) :: Triangle, temp

  REAL(kind=8) :: F1, F2, F3, F4

  Triangle = 0.0D0
  IF (a+b-c .LT. 0 .OR. a-b+c .LT. 0 .OR. -a+b+c .LT. 0 .OR. a+b+c+1 .LT. 0) RETURN 

  IF( MOD( a+b-c,2) .NE. 0) RETURN
  IF( MOD( a-b+c,2) .NE. 0) RETURN
  IF( MOD(-a+b+c,2) .NE. 0) RETURN
  IF( MOD( a+b+c,2) .NE. 0) RETURN

  F1 = FacLOG(( a + b - c    )/2)
  F2 = FacLOG(( a - b + c    )/2)
  F3 = FacLOG((-a + b + c    )/2)
  F4 = FacLOG(( a + b + c)/2 + 1)

  temp = F1 + F2 + F3 - F4
  Triangle = EXP(0.50D0 * temp)

!  Triangle = F1 * F2 * F3 / F4

  RETURN

END FUNCTION Triangle

FUNCTION Wigner_6j(a,b,c,d,e,f)

  IMPLICIT NONE

  ! Wigner 6j Symbol

  ! a b c !
  ! d e f !
 
!  INTERFACE
!
!     FUNCTION FacLOG(N)
!
!       IMPLICIT NONE
!       INTEGER, INTENT(IN) :: N
!       REAL(kind=8) ::  FacLOG
!     END FUNCTION FacLOG
!
!     FUNCTION Triangle(a,b,c)
!
!       IMPLICIT NONE
!       INTEGER, INTENT(IN) :: a,b,c
!       REAL(kind=8) :: Triangle
!     END FUNCTION Triangle
!
!  END INTERFACE

 INTEGER, INTENT(IN) :: a,b,c,d,e,f

  REAL(kind=8) :: Wigner_6j

  REAL(kind=8),DIMENSION(1:4) :: Delta  
  REAL(kind=8),DIMENSION(1:7) :: Fact
  REAL(kind=8) :: W6J1, W6J2, temp

  INTEGER :: t, t_min, t_max

  Wigner_6j = 0.0D0
  
  IF ( a < 0) RETURN
  IF ( b < 0) RETURN
  IF ( c < 0) RETURN
  IF ( d < 0) RETURN
  IF ( e < 0) RETURN
  IF ( f < 0) RETURN
  IF ( a < ABS(b-c) .OR. a > b+c ) RETURN
  IF ( a < ABS(e-f) .OR. a > e+f ) RETURN
  IF ( d < ABS(b-f) .OR. d > b+f ) RETURN
  IF ( d < ABS(e-c) .OR. d > e+c ) RETURN


  Delta(1) = DBLE(Triangle(a,b,c))
  Delta(2) = DBLE(Triangle(a,e,f))
  Delta(3) = DBLE(Triangle(d,b,f))
  Delta(4) = DBLE(Triangle(d,e,c))
  
  W6J1 = Delta(1)*Delta(2)*Delta(3)*Delta(4)

!  IF ( ABS(W6J1) .LT. 1.0D-10) THEN
!     Wigner_6j = 0.0D0
!     RETURN
!  END IF

  t_max = MIN((a+b+d+e)/2, (b+c+e+f)/2, (c+a+f+d)/2)

  t_min = MAX((a+b+c)/2, (a+e+f)/2, (d+b+f)/2, (d+e+c)/2)

!  IF ( t_min .GT. t_max) THEN
!     Wigner_6j = 0.0D0
!     RETURN
!  END IF

  W6J2 = 0.0D0

  DO t = t_min, t_max

     fact(1) = FacLOG( t + (- a - b - c)/2 )
     fact(2) = FacLOG( t + (- a - e - f)/2 )
     fact(3) = FacLOG( t + (- d - b - f)/2 )
     fact(4) = FacLOG( t + (- d - e - c)/2 )
     fact(5) = FacLOG(-t + (  a + b + d + e )/2 )
     fact(6) = FacLOG(-t + (  b + c + e + f )/2 )
     fact(7) = FacLOG(-t + (  c + a + f + d )/2 )

     temp = FacLOG(t+1) - ( fact(1) + fact(2) + fact(3) + fact(4) + fact(5) + fact(6) + fact(7))

     W6J2 = W6J2 + (-1)**(t) * EXP(temp)

!     W6J2 = W6J2 + (-1)**(t) * DBLE(FacLOG(t+1)) / &
!          & (fact(1)*fact(2)*fact(3)*fact(4)*fact(5)*fact(6)*fact(7))

  END DO

  Wigner_6j = W6J1 * W6J2

  RETURN

END FUNCTION Wigner_6j

FUNCTION Wigner_9j(j1,j2,j3,j4,j5,j6,j7,j8,j9)

  IMPLICIT NONE

! j1/2 j2/2 j3/2 !
! j4/2 j5/2 j6/2 !
! j7/2 j8/2 j9/2 !

!  INTERFACE
!
!  FUNCTION Wigner_6j(a,b,c,d,e,f)
!
!       IMPLICIT NONE
!       INTEGER, INTENT(IN) :: a,b,c,d,e,f
!       REAL(kind=8) :: Wigner_6j
!     END FUNCTION Wigner_6j
!
!  END INTERFACE

  INTEGER, INTENT(IN) :: j1,j2,j3,j4,j5,j6,j7,j8,j9
  REAL(kind=8) :: Wigner_9j

  INTEGER :: x_min, x_max, x

  Wigner_9j = 0.0D0

  IF(j1 .LT. 0) RETURN
  IF(j2 .LT. 0) RETURN
  IF(j3 .LT. 0) RETURN
  IF(j4 .LT. 0) RETURN
  IF(j5 .LT. 0) RETURN
  IF(j6 .LT. 0) RETURN
  IF(j7 .LT. 0) RETURN
  IF(j8 .LT. 0) RETURN
  IF(j9 .LT. 0) RETURN

  x_min = MAX( ABS(j1-j9), ABS(j4-j8), ABS(j2-j6))
  x_max = MIN(     j1+j9,      j4+j8,      j2+j6)

  DO x = x_min, x_max

     Wigner_9j = Wigner_9j + (-1)**(x) * (x + 1.0D0) * &
          & sj2i_lookup(j1,j4,j7,j8,j9, x) *&
          & sj2i_lookup(j2,j5,j8,j4, x,j6) *&
          & sj2i_lookup(j3,j6,j9, x,j1,j2)
          !& Wigner_6j(j1,j4,j7,j8,j9, x)  * &
          !& Wigner_6j(j2,j5,j8,j4, x,j6) * &
          !& Wigner_6j(j3,j6,j9, x,j1,j2)

  END DO

  RETURN

END FUNCTION Wigner_9j
end module wignerfunctions
