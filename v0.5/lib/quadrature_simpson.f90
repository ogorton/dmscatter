subroutine simpson(Np,N,array,h,ans)
!   subroutine to compute integral using simpson's rule
!   CWJ - SDSU - August 2006
!
!  INPUT:
!    Np: declared dimension of array
!    N : used dimension of array (note, must be even!) 
!    array(0:Np): array of points of function to be integrated
!    h :  dx 
!
!  OUTPUT:
!    ans : integral
!
      implicit none

!..........INPUT...............
      integer Np, N        ! dimensions
      real(kind=8) array(Np)     ! this is legal, since being fed from outside
      real(kind=8) h               ! = dx

!.........OUTPUT................
      real(kind=8) ans

!.........INTERMEDIATE..........
      integer i

!............. ERROR TRAP...........
!..... check that N is even.........
      if( (n/2)*2 .ne. n)then
        print*,' Dimension of array is not even ',N
        stop
      endif
!................. END ERROR TRAP

      ans = 0.0            ! initialize output
      do i = 2,N-1,2 
        ans = ans+array(i-1)+ 4*array(i)+array(i+1)
      enddo
      ans = ans*h/3.

      return

end subroutine
