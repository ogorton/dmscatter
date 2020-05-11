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

!.............................................Boole's rule
      subroutine boole(N,arr,dx,ans)
!INPUT: 
!       N: DIMENSION OF FUNCTION ARRAY
!     arr: ARRAY OF FUNCTION VALUES TO BE INTEGRATED
!      dx: STEP SIZE OF FUNCTION VALUES
!OUTPUT:
!     ans: DEFINITE INTEGRAL OF FUNCTION ARRAY FROM BOOLES RULE
!CALLS: 
!     NONE   
! 
      implicit none            
      integer N, i, m
      real(kind=8) arr(N),ans
      integer x0
      real(kind=8) dboole, boole_out,dx  , intfx 
                         
      m = N/4
      boole_out = 0 
           
      do i = 1, m-1, 1 
        x0 = 4*i                       
        dboole = ( 7.*arr(x0+4) &
                   + 32.*arr(x0+3) &
                   + 12.*arr(x0+2) &
                   + 32.*arr(x0+1) &
                   + 7.*arr(x0) )        
        boole_out = boole_out + dboole
      enddo  
                 
      ans=dx*(2./45)*boole_out
         
      end subroutine  



