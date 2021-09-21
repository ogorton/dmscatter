!===================================================================
subroutine openresults(resfile)
   use parameters
   implicit none
   integer resfile

   character(100):: filename
   integer ilast

   logical success

   success = .false.
   print*,' '
   do while(.not.success)

       print*,' Enter name of one-body density file (.dres) '

       read(5,'(a)')filename
       ilast = index(filename,' ')-1
       print*,filename(1:ilast)//'.dres'
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

    use orbitals
    use parameters
    implicit none

    type(nucleus) :: nuc_target

    ! densities(J,iso,a,b)
    allocate(nuc_target%densitymats%rho( 0:10,0:1,1:ntotal(1),1:ntotal(1)) )
    print*,'Density matrix allocated, size:',size(nuc_target%densitymats%rho)
    nuc_target%densitymats%rho(:,:,:,:) = 0.0

end subroutine setupdensities


!===================================================================
subroutine coredensity(nuc_target)

  use orbitals
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
          print*,"Orbital",i,sqrt((jorb(i)+1.0) * (Jiso+1.0))
          nuc_target%densitymats%rho(0,0,i,i) = sqrt((jorb(i)+1.0) * (Jiso+1.0))
          nuc_target%densitymats%rho(0,1,i,i) = sqrt((jorb(i)+1.0) * (Jiso+1.0))
      end do
  else
      do i = norb(1)+1, ntotal(1) ! core comes after valence. ntotal = ncore+nval.
         print*,"Orbital",i
         nuc_target%densitymats%rho(0,0,i,i) = sqrt(2.0*(jorb(i)+1.0)*(Jiso+1.0)*(Tiso+1.0))
         nuc_target%densitymats%rho(0,1,i,i) = 0.0
      end do
  end if

end subroutine coredensity

!===================================================================

subroutine printdensities(nuc_target)
    use parameters
    use orbitals
    implicit none
    type(nucleus) :: nuc_target
    integer J,a, b
    print*,'Printing density matrix.'
    print*,'# spo =', ntotal(1)
    do J=0,10
        print*,'J=',J
        do a=1,ntotal(1)
            do b=1,ntotal(1)
                if (nuc_target%densitymats%rho(J,0,a,b)&
                    +nuc_target%densitymats%rho(J,1,a,b) .eq. 0) cycle
                print*,a,b,nuc_target%densitymats%rho(J,0,a,b),&
                    nuc_target%densitymats%rho(J,1,a,b)
            enddo
        enddo
    enddo

end subroutine printdensities


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
subroutine readheaderv2(nuc_target, resfile)
   use parameters
   implicit none
   type(nucleus) :: nuc_target
   integer resfile
   character(23) :: tmpline
   integer i,j,n
   real(kind=8) e,ex,xj,xt
   integer np,zp
   integer closest2J   ! function to convert poorly converged values
      rewind(resfile)
      read(resfile,'(a)')tmpline
      read(resfile,*)zp,np
      if( mod(np+zp,2) == 1)then
         evenA = .false.
      else
         evenA = .true.
      end if
      rewind(resfile)
   do i = 1,20
      read(resfile,'(a)')tmpline
      if(tmpline(3:7)=='State')then

            nlocalstates = 0
            do j = 1,50000

               read(resfile,*,err=3)n

               if(n==1)then
                 backspace(resfile)
                 read(resfile,*,err=3)n,e,ex,xj,xt

               end if

              if(n==j)then
                   nlocalstates = nlocalstates +1
               else
                   exit
               end if
            end do

3           continue

           nuc_target%groundstate%Jx2 = closest2J(evenA,xj)
           nuc_target%groundstate%Tx2 = closest2J(evenA,xt)
         backspace(resfile)
         return
      end if

   end do
   print*,' Did not find header '
   stop

end subroutine readheaderv2


!-------------------------------------------------------------------------------
subroutine read2state(resfile,locchar,n,found,finished)
   implicit none
   integer resfile
   character(1) :: locchar
   integer n
   logical found,finished
   character(4) :: myloc
   !integer m
   found = .false.
   select case (locchar)
     case ('i')
     do while(.not.found)
        read(resfile,1,end=111)myloc
1 format(a4)
        if(myloc(2:4)=='Ini')then
        backspace(resfile)
        read(resfile,11)myloc,n
!       print*, 'int', n
11 format(a4,13x,i4)

           found = .true.
           exit

        endif
      end do

      case('f')

        read(resfile,11)myloc,n
        if(myloc(2:4)=='Fin')found=.true.

!        print*, 'fin',n
   end select
   finished = .false.
   return
111 continue
   finished = .true.
   return

end subroutine read2state


subroutine read2Jtrans(resfile,found)
   use parameters
!   use op_info,only:pndens
   implicit none
   integer resfile
   logical found
   character(3) :: tmpchar
   !integer j

   read(resfile,'(a3)',end=111)tmpchar
   if(tmpchar(2:3) == ' ' .or. tmpchar(2:3)=='In' .or. tmpchar(1:2)=='++')then
      found = .false.
      return
   endif
   if(tmpchar(2:3) == 'Jt')then
     backspace(resfile)
     read(resfile,'(5x,i4)')jt

   end if
   found = .true.
!..... ADDED IN VERSION 9 MAR 2017.... CHECK IF ISOSPIN OR PN FORMALISM...
   backspace(resfile)
   read(resfile,'(11x,a3)')tmpchar

   if(tmpchar=='pro')then
           pndens=.true.
           print*,'Density matrix is in explicit proton-neutron format.'
   else
           print*,"Density matrix is in isospin format."
           print*,"This format is not currently supported."
!           STOP "Incompatible density matrix format."
           pndens=.false.
   end if

   return
111 continue
   found=.false.
   return
end subroutine read2Jtrans


subroutine readdensity(nuc_target, resfile,success)
   use parameters
   use orbitals
   implicit none
   type(nucleus) :: nuc_target
   integer resfile
   integer a,b,i
   real(kind=8) ops,opv
   logical :: success

   success=.false.

    do i = 1,norb(1)*norb(1)!nsporb*nsporb
        read(resfile,*,err=1,end=1)a,b,ops,opv

             if(ops /= 0.0)then
               nuc_target%densitymats%rho(jt,0,a,b)= ops
               success=.true.
             end if

             if(opv /= 0.0)then
               nuc_target%densitymats%rho(jt,1,a,b) = opv
               success=.true.
             end if

        !end if
   end do

   return

1  continue
   backspace(resfile)
   return
end subroutine readdensity


subroutine readalldensities(nuc_target,resfile)
!   use op_info
!   use spspace
   use parameters
   implicit none
   type(nucleus) :: nuc_target
   integer resfile
   integer istate,fstate
   logical foundi,foundf,foundjt
   logical endoffile,endoflist
   logical finished,success
   logical nodensities  ! flag to make sure densities are not empty

   endoffile = .false.
   nodensities=.true.

   do while(.not.endoffile)
      call read2state(resfile,'i',istate,foundi,finished)
      if(finished)exit
      if(.not.foundi)then
           endoffile = .true.
           exit
      end if
      call read2state(resfile,'f',fstate,foundf,finished)

      if(.not.foundf)then
           endoffile = .true.
           exit
      end if
      if(finished)exit

      if(istate .ne. 1 .or. fstate .ne. 1) cycle
      endoflist = .false.
      do while(.not.endoflist)

          call read2Jtrans(resfile,foundjt)
          if(.not.foundjt)then
              endoflist = .true.
              backspace(resfile)
              exit
          end if

          call readdensity(nuc_target, resfile,success)
          if(success)nodensities=.false.
!          end if
      end do ! endoflist
      exit

   end do  ! endoffile

   if(nodensities)then
          print*,' Wait! That density file held no densities ! '
   end if

   return
end subroutine readalldensities
