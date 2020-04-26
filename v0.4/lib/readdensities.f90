subroutine readheaderv2(resfile)
   use stateinfo
!   use op_info,only:pndens,pnops
   implicit none
   integer resfile
!   character(1) :: ctrlchar
   character(23) :: tmpline
   integer i,j,n
   real e,ex,xj,xt
!   real etol
!   real, pointer :: ee(:)
!   integer, pointer :: jx2(:),tx2(:)
!   integer nmax
   integer np,zp
   integer closest2J   ! function to convert poorly converged values
!   logical askfortshift

!   askfortshift=.false.
!
!   print*,' reading file ',ctrlchar
!   tshift = 0.0
!   select case (ctrlchar)
!     case ('f')
!        ee => energy
!        Jx2 => Jx2state
!        Tx2 => Tx2state
!        nmax = nlocalstates
!                tshift = 0.0
!                if(askfortshift)then
!        print*,' Enter any shift for isospin ( x T(T+1) for all states )'
!        print*,' (Typical value = 0, no shift )'
!        read*,tshift
!            end if
!     case ('p')
!        ee => energy_parent
!        Jx2 => Jx2state_parent
!        Tx2 => Tx2state_parent
!        nmax = nparent_ref
!     case ('d')
!        ee => energy_daughter
!        Jx2 => Jx2state_daughter
!        Tx2 => Tx2state_daughter
!        nmax = ndaughter_ref
!
!     case default
!
!   end select
!
!   etol = 1.0e-3
!
!   rewind(resfile)
!............ check whether EVEN or ODD..............
!   if(ctrlchar=='p')then
      read(resfile,'(a)')tmpline
      read(resfile,'(a)')tmpline
      read(resfile,*)zp,np
      print*,zp,np
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
               print*, n,e,ex,xj,xt

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

           Jiso = closest2J(evenA,xj)
           Tiso = closest2J(evenA,xt)
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
       print*, 'int', n
11 format(a4,13x,i4)

           found = .true.
           exit

        endif
      end do

      case('f')

        read(resfile,11)myloc,n
        if(myloc(2:4)=='Fin')found=.true.

        print*, 'fin',n
   end select
   finished = .false.
   return
111 continue
   finished = .true.
   return

end subroutine read2state


subroutine read2Jtrans(resfile,found)
   use stateinfo
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


subroutine readdensity(resfile,success)
   use stateinfo
!   use op_info
   implicit none
   integer resfile
   !integer istate,fstate
   integer a,b,i
   real ops,opv
   !real fact0t,fact1t  ! isospin factors
   logical :: success
   !real cleb !       ! function from LIBRA.f

   success=.false.

!   fact0t = cleb(Tx2state_parent(istate),Mzi,0,0,Tx2state_daughter(fstate),Mzf)*sqrt(2.)/sqrt(Tx2state_daughter(fstate)+1.)
!   fact1t = cleb(Tx2state_parent(istate),Mzi,2,0,Tx2state_daughter(fstate),Mzf)*sqrt(6.)/sqrt(Tx2state_daughter(fstate)+1.)

    do i = 1,nsporb*nsporb
        read(resfile,*,err=1,end=1)a,b,ops,opv

!          if(j==jt .and. istate >0 .and. fstate > 0)then
!             if(densitymats(istate,fstate)%good)then
        if(pndens)then
            if(ops/=0.0)then
              densitymats%rhop(jt,a,b)= ops
              success=.true.
            end if
            if(opv/=0.0)then
              densitymats%rhon(jt,a,b)= opv
              success=.true.
            end if
        else

             if(ops /= 0.0)then
               densitymats%rho(jt,0,a,b)= ops
               success=.true.
             end if

             if(opv /= 0.0)then
               densitymats%rho(jt,1,a,b) = opv
               success=.true.
             end if

!.......... CONVERT DENSITIES FROM ISOSPIN TO PROTON-NEUTRON......................
!                       if(success)then
!                                                densitymats(istate,fstate)%rhop(a,b) = 0.5*(fact0t*ops+fact1t*opv)
!                                                densitymats(istate,fstate)%rhon(a,b) = 0.5*(fact0t*ops-fact1t*opv)
!                                          end if                                
!                                  else  ! MUST CONVERT
!                                          if(ops/=0.0 .or. opv /=0.0)then
!                                                  success=.true.
!                                                densitymats(istate,fstate)%rhop(a,b) = 0.5*(fact0t*ops+fact1t*opv)
!                                                densitymats(istate,fstate)%rhon(a,b) = 0.5*(fact0t*ops-fact1t*opv)                                            
!  
!                                          end if
!                                        
!                                  end if        
!                                
        end if
!             end if
!          end if
   end do

   return

1  continue
   backspace(resfile)
   return
end subroutine readdensity


subroutine readalldensities(resfile)
!   use op_info
!   use spspace
   use stateinfo
   implicit none
   integer resfile
   integer istate,fstate
   !integer a,b,i,j
   !real ops,opv
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
          print *, endoflist, jt
!          if(map2parent(istate) > 0 .and. map2daughter(fstate) > 0)then
!                         print*,' states ',istate,map2parent(istate),fstate,map2daughter(fstate)

          call readdensity(resfile,success)
          if(success)nodensities=.false.
!          end if
      end do ! endoflist
!CFJIAO
      exit

   end do  ! endoffile

   if(nodensities)then
          print*,' Wait! That density file held no densities ! '
   end if

!  fill the one-body density matrix of the core

   call coredensity

   return
end subroutine readalldensities





