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
      real(kind=8) dboole, boole_out,dx
                         
      m = N/4
      boole_out = 0 
           
      do i = 1, m-1, 1 
        x0 = 4*i                       
        dboole = ( 7d0*arr(x0+4) &
                   + 32d0*arr(x0+3) &
                   + 12d0*arr(x0+2) &
                   + 32d0*arr(x0+1) &
                   + 7d0*arr(x0) )        
        boole_out = boole_out + dboole
      enddo  
                 
      ans=dx*(2d0/45d0)*boole_out
         
      end subroutine  
