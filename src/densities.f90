!===================================================================
subroutine openresults(resfile)
   use parameters
   implicit none
   integer resfile

   character(22):: filename
   integer ilast

   logical success

   success = .false.
   print*,' '
   do while(.not.success)

       print*,' Enter name of one-body density file (.dres) '

       read(5,'(a)')filename
       ilast = index(filename,' ')-1
       open(unit=resfile,file=filename(1:ilast)//'.dres',status='old',err=2)
       success = .true.
       return
2      continue
       print*,filename(1:ilast),'.res does not exist '

   end do

   return
end subroutine openresults


!===================================================================
subroutine setupdensities(nuc_target)

    use spspace
    use parameters
!    use op_info
    implicit none

    type(nucleus) :: nuc_target


    ! densities(J,iso,a,b)
    allocate(nuc_target%densitymats%rho( 0:10,0:1,1:ntotal(1),1:ntotal(1)) )
    print*,'Density matrix allocated, size:',size(nuc_target%densitymats%rho)
    nuc_target%densitymats%rho(:,:,:,:) = 0.0

end subroutine setupdensities


!===================================================================
subroutine coredensity(nuc_target)

  use spspace
  use parameters

  implicit none

  type(nucleus) :: nuc_target
  integer :: i
  integer :: Jiso, Tiso

  Jiso = nuc_target%groundstate%Jx2
  Tiso = nuc_target%groundstate%Tx2

  print*,'Filling core orbitals.'

  if (pndens) then
      do i = norb(1)+1, ntotal(1)
          nuc_target%densitymats%rho(0,0,i,i) = sqrt((jorb(i)+1.0) * (Jiso+1.0))
          nuc_target%densitymats%rho(0,1,i,i) = sqrt((jorb(i)+1.0) * (Jiso+1.0))
      end do
  else
      do i = norb(1)+1, ntotal(1) ! core comes after valence. ntotal = ncore+nval.
         nuc_target%densitymats%rho(0,0,i,i) = sqrt(2.0*(jorb(i)+1.0)*(Jiso+1.0)*(Tiso+1.0))
         nuc_target%densitymats%rho(0,1,i,i) = 0.0
      end do
  end if

end subroutine coredensity

!===================================================================

subroutine printdensities(nuc_target)
    use parameters
    use spspace
    implicit none
    type(nucleus) :: nuc_target
    integer J,a, b
    print*,'Printing density matrix.'
    print*,'# spo =', ntotal(1)
    do J=0,10
        print*,'J=',J
        do a=1,ntotal(1)
            do b=1,ntotal(1)
                if (nuc_target%densitymats%rho(J,0,a,b)+nuc_target%densitymats%rho(J,1,a,b) .eq. 0) cycle
                print*,a,b,nuc_target%densitymats%rho(J,0,a,b),nuc_target%densitymats%rho(J,1,a,b)
            enddo
        enddo
    enddo

end subroutine printdensities

!===================================================================


!================================================
!
!  function to force conversion of unconverged xJ to integer J
!  that is, odd 2 x J for odd A, and even for even A
!
  function closest2J(evenA,xj)

  implicit none
  integer closest2J
  real(kind=8) xj
  logical evenA

  if(evenA)then
     closest2J = 2*nint(xj)
     if(closest2J < 0)closest2J = 0
  else
     closest2J = 2*nint(xj-0.5)+1
     if(closest2J < 1)closest2J = 1
  end if

  return
  end function closest2J
!================================================
