subroutine readheaderv2(nuc_target, resfile)
   use parameters
!   use op_info,only:pndens,pnops
   implicit none
   type(nucleus) :: nuc_target
   integer resfile
!   character(1) :: ctrlchar
   character(23) :: tmpline
   integer i,j,n
   real(kind=8) e,ex,xj,xt
!   real(kind=8) etol
!   real(kind=8), pointer :: ee(:)
!   integer, pointer :: jx2(:),tx2(:)
!   integer nmax
   integer np,zp
   integer closest2J   ! function to convert poorly converged values
!   logical askfortshift

!............ check whether EVEN or ODD..............
!   if(ctrlchar=='p')then
      read(resfile,'(a)')tmpline
      read(resfile,'(a)')tmpline
      read(resfile,*)zp,np
!      print*,zp,np
!          if(.not.readinZN)then
!                  Z0 = zp
!                  N0 = np
!          else
!                  if(Z0/=zp .or. N0/=np)then
!                          print*,' Mismatch in valence particles '
!                          print*,' Expect Z,N = ',Z0,N0
!                          print*,' Found ',zp,np
!                          stop
!                  end if
!          end if
      if( mod(np+zp,2) == 1)then
         evenA = .false.
      else
         evenA = .true.
      end if
      rewind(resfile)
!   end if
   do i = 1,20
      read(resfile,'(a)')tmpline
      if(tmpline(3:7)=='State')then

!         select case (ctrlchar)
!         case ('c')

            nlocalstates = 0
            do j = 1,50000

               read(resfile,*,err=3)n

               if(n==1)then
                 backspace(resfile)
                 read(resfile,*,err=3)n,e,ex,xj,xt
         !      print*, n,e,ex,xj,xt

               end if

              if(n==j)then
                   nlocalstates = nlocalstates +1
!                   print*, n,j
               else
!                   print*, n,j
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
   else
           pndens=.false.
   end if

   return
111 continue
   found=.false.
   return
end subroutine read2Jtrans


subroutine readdensity(nuc_target, resfile,success)
   use parameters
   use spspace
   implicit none
   type(nucleus) :: nuc_target
   integer resfile
   integer a,b,i
   real(kind=8) ops,opv
   logical :: success

   success=.false.

    do i = 1,norb(1)*norb(1)!nsporb*nsporb
        read(resfile,*,err=1,end=1)a,b,ops,opv
        print*,a,b,ops,opv

        if(pndens)then
            if(ops/=0.0)then
              nuc_target%densitymats%rhop(jt,a,b)= ops
              success=.true.
            end if
            if(opv/=0.0)then
              nuc_target%densitymats%rhon(jt,a,b)= opv
              success=.true.
            end if
        else

             if(ops /= 0.0)then
               nuc_target%densitymats%rho(jt,0,a,b)= ops
               success=.true.
             end if

             if(opv /= 0.0)then
               nuc_target%densitymats%rho(jt,1,a,b) = opv
               success=.true.
             end if

        end if
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
!      print *, 'foundi= ', foundi, 'finishedi = ',finished
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
!CFJIAO
      exit

   end do  ! endoffile

   if(nodensities)then
          print*,' Wait! That density file held no densities ! '
   end if

   return
end subroutine readalldensities


