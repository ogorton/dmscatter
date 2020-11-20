!===================================================================
subroutine openresults(resfile)
   use stateinfo
   implicit none
   integer resfile

   character(22):: filename
   integer ilast

   logical success

   success = .false.
   print*,' '
   do while(.not.success)

       print*,' Enter name of one-body density file (.res) '

       read(5,'(a)')filename
       ilast = index(filename,' ')-1
       open(unit=resfile,file=filename(1:ilast)//'.res',status='old',err=2)
       success = .true.
       return
2      continue
       print*,filename(1:ilast),'.res does not exist '

   end do

   return
end subroutine openresults


!===================================================================
subroutine setupdensities

    use spspace
    use stateinfo
!    use op_info
    implicit none

    print*, ntotal(1)

!                 densitymats%good = .true.
    ! densities(J,iso,a,b)
    allocate(densitymats%rho( 0:10,0:1,1:ntotal(1),1:ntotal(1)) )
    densitymats%rho(:,:,:,:) = 0.0
    allocate(densitymats%rhop(0:10,1:ntotal(1),1:ntotal(1)))
    allocate(densitymats%rhon(0:10,1:ntotal(1),1:ntotal(1)))
    densitymats%rhop(:,:,:) = 0.0
    densitymats%rhon(:,:,:) = 0.0

end subroutine setupdensities


!===================================================================
subroutine coredensity

  use spspace
  use stateinfo

  integer :: i

  print*,'Filling core orbitals.'

  do i = norb(1)+1, ntotal(1) ! core comes after valence. ntotal = ncore+nval.
     densitymats%rho(0,0,i,i) = sqrt(2.0*(jorb(i)+1.0)*(Jiso+1.0)*(Tiso+1.0))
     densitymats%rho(0,1,i,i) = sqrt(2.0*(jorb(i)+1.0)*(Jiso+1.0)*(Tiso+1.0))
  end do

end subroutine coredensity

!===================================================================

subroutine printdensities
    use stateinfo
    use spspace
    implicit none
    integer J,a, b
    print*,'Printing density matrix.'
    print*,'# spo =', ntotal(1)
    do J=0,10
        print*,'J=',J
        do a=1,ntotal(1)
            do b=1,ntotal(1)
                if (densitymats%rho(J,0,a,b)+densitymats%rho(J,1,a,b) .eq. 0) cycle
                print*,a,b,densitymats%rho(J,0,a,b),densitymats%rho(J,1,a,b)
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
  real xj
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
