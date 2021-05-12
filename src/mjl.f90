module mjl
    use norm
    use wignerfunctions
    use bessel
    contains

!-------------------------------------------------------------------------------
function MJ(y,np,lp,jp,n,l,j,bigJ)
    implicit none

    INTEGER, INTENT(IN) :: np,lp,jp,n,l,j,bigJ
    REAL(kind=8), INTENT(IN) :: y
    REAL(kind=8) :: xjp,xj
    REAL(kind=8) :: Pi = 3.1415926535897932
    !REAL(kind=8) :: DBLEFacLOG
    REAL(kind=8) :: MJ

    xjp = dble(jp)/2.0
    xj  = dble(j )/2.0

    MJ = (-1)**(0.5d0+xj+dble(bigJ))&
        *SQRT(Jnorm(xj)*Jnorm(xjp) &
        *Jnorm(dble(l))*Jnorm(dble(lp))*Jnorm(dble(bigJ))/(4*Pi)) &
          * tj2i_lookup(2*lp,2*bigJ,2*l,0,0,0) &
          * sj2i_lookup(2*lp,jp,1,j,2*l,2*bigJ)   &
          * BesselElement(y,np,lp,n,l,bigJ)

end function MJ


!-------------------------------------------------------------------------------
function MJLDivQoverall(lp,jp,l,j,bigJ,bigL)
    implicit none

    INTEGER, INTENT(IN) :: lp,jp,l,j,bigJ,bigL
    REAL(kind=8) :: xjp,xj
    REAL(kind=8) :: Pi = 3.1415926535897932
    REAL(kind=8) :: MJLDivQoverall

    xjp = dble(jp)/2.0
    xj  = dble(j )/2.0

    ! Error trap
    ! (-1.0)**(bigL+j/2+1/2)
    MJLDivQoverall = (-1.0)**(bigL+j/2+1/2)* Qnorm(dble(lp))*Qnorm(xjp) &
        *Qnorm(xj)*Qnorm(dble(bigJ))*Qnorm(dble(bigL)) &
        *sj2i_lookup(2*lp,jp,1,j,2*l,2*bigJ)/SQRT(4*Pi)

end function MJLDivQoverall


!-------------------------------------------------------------------------------
function MJLDivQsummand1(y,np,lp,jp,n,l,j,bigJ,bigL)
    implicit none

    INTEGER, INTENT(IN) :: np,lp,jp,n,l,j,bigJ,bigL
    REAL(kind=8), INTENT(IN) :: y
    REAL(kind=8) :: MJLDivQsummand1

    MJLDivQsummand1 = -SQRT(dble(l)+1.0)*Qnorm(dble(l)+1.0) &
        * sj2i_lookup(2*bigL,2,2*bigJ,2*l,2*lp,2*(l+1))  &
        * tj2i_lookup(2*lp,2*bigL,2*(l+1),0,0,0) &
        * BesselElementminus(y,np,lp,n,l,bigL)

end function MJLDivQsummand1


!-------------------------------------------------------------------------------
function MJLDivQsummand2(y,np,lp,jp,n,l,j,bigJ,bigL)

    ! l must be .ge. 1 for Qnorm function.
    use wignerfunctions
    implicit none

    INTEGER, INTENT(IN) :: np,lp,jp,n,l,j,bigJ,bigL
    REAL(kind=8), INTENT(IN) :: y
    REAL(kind=8) :: MJLDivQsummand2

    if (l.eq.0) then ! O.C.G: avoid Qnorm(-1)
        MJLDivQsummand2 = 0.0
    else

        MJLDivQsummand2 = SQRT(dble(l))*Qnorm(dble(l)-1.0) &
            *sj2i_lookup(2*bigL,2,2*bigJ,2*l,2*lp,2*(l-1))  &
            *tj2i_lookup(2*lp,2*bigL,2*(l-1),0,0,0) &
            *BesselElementplus(y,np,lp,n,l,bigL)

    endif

end function MJLDivQsummand2


!-------------------------------------------------------------------------------
function MJLDivQ(y,np,lp,jp,n,l,j,bigJ,bigL)
    use wignerfunctions
    implicit none

     INTEGER, INTENT(IN) :: np,lp,jp,n,l,j,bigJ,bigL
     REAL(kind=8), INTENT(IN) :: y
     REAL(kind=8) :: MJLDivQ, tmp

     MJLDivQ = MJLDivQoverall(lp,jp,l,j,bigJ,bigL) &
         * (MJLDivQsummand1(y,np,lp,jp,n,l,j,bigJ,bigL) &
         +MJLDivQsummand2(y,np,lp,jp,n,l,j,bigJ,bigL))

end function MJLDivQ


end module mjl
