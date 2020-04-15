module stateinfo
   implicit none
   logical :: evenA

   integer nlocalstates !,maxlocalstates

   integer nsporb   ! # of single-particle orbits

   integer jt, tt   ! spin, isospin of transition operator
   integer Jiso, Tiso   ! ang, isospin of ground state
   logical ::  pndens 
   type onebdenmat
        logical good
        real, allocatable :: rho(:,:,:,:)
                real, allocatable :: rhop(:,:,:),rhon(:,:,:)  ! added in version 9

   end type onebdenmat
   type (onebdenmat) :: densitymats

end module stateinfo

!===================================================================

      module spspace
         implicit none
!------------------SINGLE-PARTICLE STATES---------------------------   

         integer norb(2)                ! # of single-particle j-orbits
                                        ! 1 = p,  2 = n
         integer nsps(2)                ! # of single-particle m-states
                                        ! 1 = p, 2 = n

         integer ncore(2)               ! # of single-particle j-orbits in the core
                                        ! 1 = p,  2 = n
         integer ntotal(2)

         integer,allocatable  :: jorb(:)        ! 2 x J of orbit (note: p,n different)
         integer,allocatable  :: lorb(:)
         integer,allocatable  :: nrorb(:)       ! radial quantum number N of j-orbit
         integer,allocatable  :: torb(:)        ! 2 x Tz of orbit ( p/n = +/-1 )
         integer,allocatable  :: nprincipal(:)
         integer,allocatable  :: nodal(:)
                                
        logical spinless
      end module spspace

!===================================================================

subroutine setupdensities

   use spspace
   use stateinfo
!   use op_info
   implicit none

   print*, ntotal(1)

!                 densitymats%good = .true.
                 allocate(densitymats%rho( 0:10,0:1,1:ntotal(1),1:ntotal(1)) )
                 densitymats%rho(:,:,:,:) = 0.0
                 allocate(densitymats%rhop(0:10,1:ntotal(1),1:ntotal(1)))
                 allocate(densitymats%rhon(0:10,1:ntotal(1),1:ntotal(1)))
                 densitymats%rhop(:,:,:) = 0.0
                 densitymats%rhon(:,:,:) = 0.0

end subroutine setupdensities


!  ctrlchar = 'c'  count up # of states
!  ctrlchar = 'f'  fill up info on states
!  ctrlchar = 'p'  fill up parent reference states
!  ctrlchar = 'd'  fill up daughter reference states

subroutine readheaderv2(resfile)
   use stateinfo
!   use op_info,only:pndens,pnops
   implicit none
   integer resfile
!   character(1) :: ctrlchar
   character(23) :: tmpline
   integer i,j,n,k
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

subroutine read2state(resfile,locchar,n,found,finished)
   implicit none
   integer resfile
   character(1) :: locchar
   integer n
   logical found,finished
   character(4) :: myloc
   integer m
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
   integer j

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
   integer istate,fstate
   integer a,b,i
   real ops,opv
   real fact0t,fact1t  ! isospin factors
   logical :: success
   real cleb !       ! function from LIBRA.f

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
   integer a,b,i,j
   real ops,opv
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


Program darkmattermain
    use spspace
    use stateinfo
    implicit none

    interface
       subroutine openresults(resfile)
       use stateinfo
       implicit none
       integer resfile
       character(22):: filename

       end subroutine openresults

    end interface

    integer, parameter :: resfile = 33
    integer :: i,j,a,b,ap,an
    character(23) :: tmpline
    
    integer :: ioption, op1, op2, Mtiso
    integer :: t,tt1,tt2
    REAL(kind=8) :: Wigner_3j
    REAL(kind=8) :: Pi = 3.1415926535897932
    REAL(kind=8) :: y, spOME, spOME1,spOME2
    REAL(kind=8), dimension (0:1) :: IsoME,IsoME1,IsoME2
    REAL(kind=8), dimension (0:6,0:1) :: DRME, SRME, DRME1, SRME1, DRME2, SRME2 
    REAL(kind=8), dimension (0:1,0:1) :: Response
    logical :: success

    call GetSPS

    print*,' '
    print*,' Reading one-body density matrix from Bigstick .res file '
    print*,' '

    nsporb = norb(1)

    call openresults(resfile)
    call setupdensities

    call readheaderv2(resfile)

    call readalldensities(resfile)

   print*,' '

       print*,' Enter the neutron number '

       read(5,*)an

   print*,' '

       print*,' Enter the proton number '

       read(5,*)ap

   Mtiso = ap-an

   success = .false.
   print*,' '
   do while(.not.success)

       print*," Enter options of operator (1-8) "
       print*," (1) W_M "
       print*," (2) W_{\PhiPP}"
       print*," (3) W_{\tilde{\Phi}'} "
       print*," (4) W_{\Delta} "
       print*," (5) W_{\SigmaP} "
       print*," (6) W_{\SigmaPP} "
       print*," (7) W_{\PhiPP M} "
       print*," (8) W_{\Delta \SigmaP} "

       read(5,*)ioption
       
       if (ioption .gt. 8 .or. ioption .lt. 1) goto 2
       success = .true.
       goto 1     

2      continue
       print*,' options should be 1 to 8 '
1      continue

   end do

   if (ioption .eq. 1) then   
       op1 = 1; op2 = 1
   else if (ioption .eq. 2) then
       op1 = 3; op2 = 3
   else if (ioption .eq. 3) then
       op1 = 4; op2 = 4
   else if (ioption .eq. 4) then
       op1 = 5; op2 = 5
   else if (ioption .eq. 5) then
       op1 = 6; op2 = 6
   else if (ioption .eq. 6) then
       op1 = 7; op2 = 7
   else if (ioption .eq. 7) then
       op1 = 3; op2 = 1
   else if (ioption .eq. 8) then
       op1 = 5; op2 = 6

   end if

   print*,' '

       print*,' Enter the value of y '

       read(5,*)y

    IsoME(0:1)=0.d0
    DRME(0:6,0:1) = 0.d0
    SRME(0:6,0:1) = 0.d0

    IsoME1(0:1)=0.d0
    DRME1(0:6,0:1) = 0.d0
    SRME1(0:6,0:1) = 0.d0

    IsoME2(0:1)=0.d0
    DRME2(0:6,0:1) = 0.d0
    SRME2(0:6,0:1) = 0.d0

    Response(0:1,0:1) = 0.d0

   if (ioption .le. 3) then
      Do j = 0,6,2
       Do t = 0, 1
        Do a = 1, ntotal(1)
         Do b = 1, ntotal(1)
      If (abs(densitymats%rho(j,t,a,b)) .ge. 1.0e-9) then
!           write(*,"(i4,i4,i4,i4,i4,i4,f9.5)") &
!           & jt,tt,nprincipal(a),jorb(a),nprincipal(b),jorb(b), densitymats%rho(jt,tt,a,b)

           call OperME(op1,y,nodal(a),lorb(a),jorb(a),nodal(b),lorb(b),jorb(b),j,spOME)
!           print*, nodal(a),lorb(a),jorb(a),nodal(b),lorb(b),jorb(b), spOME

           IsoME(t)= spOME * sqrt(2.d0) *sqrt(2*dble(t)+1.0)
           DRME(j,t) = DRME(j,t) + densitymats%rho(j,t,a,b)*IsoME(t)
           SRME(j,t) = (-1.0)**((Tiso - Mtiso)/2)*Wigner_3j(Tiso,2*t,Tiso,-Mtiso,0,Mtiso)*DRME(j,t)    
      end if
         end do 
       end do
      end do
     end do

     Do j = 0, 6, 2
       Do tt1 = 0,1
         Do tt2 = 0,1
          Response(tt1,tt2) = Response(tt1,tt2)+ SRME(j,tt1)*SRME(j,tt2)
         end do
       end do
     end do

   else if (ioption .ge. 4 .and. ioption .le. 6) then
      Do j = 1,5,2
       Do t = 0, 1
        Do a = 1, ntotal(1)
         Do b = 1, ntotal(1)
      If (abs(densitymats%rho(j,t,a,b)) .ge. 1.0e-9) then
!           write(*,"(i4,i4,i4,i4,i4,i4,f9.5)") &
!           & jt,tt,nprincipal(a),jorb(a),nprincipal(b),jorb(b), densitymats%rho(jt,tt,a,b)
         
           call OperME(op1,y,nodal(a),lorb(a),jorb(a),nodal(b),lorb(b),jorb(b),j,spOME)
!           print*, nodal(a),lorb(a),jorb(a),nodal(b),lorb(b),jorb(b), spOME
         
           IsoME(t)= spOME * sqrt(2.d0) *sqrt(2*dble(t)+1.0)
           DRME(j,t) = DRME(jt,t) + densitymats%rho(j,t,a,b)*IsoME(t)
           SRME(j,t) = (-1.0)**((Tiso - Mtiso)/2)*Wigner_3j(Tiso,2*t,Tiso,-Mtiso,0,Mtiso)*DRME(j,t)
      end if
         end do
       end do
      end do
     end do

     Do j = 1, 5, 2
       Do tt1 = 0,1
         Do tt2 = 0,1
          Response(tt1,tt2) = Response(tt1,tt2)+ SRME(j,tt1)*SRME(j,tt2)
         end do
       end do
     end do

   else if (ioption .eq. 7) then

      Do j = 0,6,2
       Do t = 0, 1
        Do a = 1, ntotal(1)
         Do b = 1, ntotal(1)
      If (abs(densitymats%rho(j,t,a,b)) .ge. 1.0e-9) then
!           write(*,"(i4,i4,i4,i4,i4,i4,f9.5)") &
!           & jt,tt,nprincipal(a),jorb(a),nprincipal(b),jorb(b), densitymats%rho(jt,tt,a,b)

           call OperME(op1,y,nodal(a),lorb(a),jorb(a),nodal(b),lorb(b),jorb(b),j,spOME1)
           call OperME(op2,y,nodal(a),lorb(a),jorb(a),nodal(b),lorb(b),jorb(b),j,spOME2)

           IsoME1(t)= spOME1 * sqrt(2.d0) *sqrt(2*dble(t)+1.0)
           DRME1(j,t) = DRME1(j,t) + densitymats%rho(j,t,a,b)*IsoME1(t)
           SRME1(j,t) = (-1.0)**((Tiso - Mtiso)/2)*Wigner_3j(Tiso,2*t,Tiso,-Mtiso,0,Mtiso)*DRME1(j,t)

           IsoME2(t)= spOME2 * sqrt(2.d0) *sqrt(2*dble(t)+1.0)
           DRME2(j,t) = DRME2(j,t) + densitymats%rho(j,t,a,b)*IsoME2(t)
           SRME2(j,t) = (-1.0)**((Tiso - Mtiso)/2)*Wigner_3j(Tiso,2*t,Tiso,-Mtiso,0,Mtiso)*DRME2(j,t)

      end if
         end do
       end do
      end do
     end do

     Do j = 0, 6, 2
       Do tt1 = 0,1
         Do tt2 = 0,1
          Response(tt1,tt2) = Response(tt1,tt2)+ SRME1(j,tt1)*SRME2(j,tt2)
         end do
       end do
     end do

   end if

   Write (*,*) "Non-zero Response for option",ioption 
   Do tt1 = 0,1
    Do tt2 = 0,1
      if (abs(Response(tt1,tt2)) .gt. 1.0d-15) then
       Write(*,"(i4,i4,f25.17)") tt1, tt2, Response(tt1,tt2)
      end if
    end do
   end do
end  

subroutine coredensity

  use spspace
  use stateinfo

  integer :: i

  do i = norb(1)+1, ntotal(1)
     densitymats%rho(0,0,i,i) = sqrt(2.0*(jorb(i)+1.0)*(Jiso+1.0)*(Tiso+1.0))
  end do

end 

      subroutine GetSPS
!
!  reads in .sps file information on j-orbits 
!
! REVISION 8 April 2004: reads in REDSTICK-compatible .sps files
! Revision 13 May 2004: use allocatable arrays
!
!
!  OUTPUT:     norb = # of j-orbits
!              jorb(i) = 2*j of ith j-orbit 
!              nrorb(i)=  radial quantum number of ith j-orbit
!              torb(i) = 2*Tz of j-orbit
!              nsps = # of m-states
!              orb_qn = quantum numbers of j-orbits, stored differently
!              spsqn  = quantum numbers of m-states
!
!  CALLED BY: main

      use spspace
!      use chf_dim
      implicit none

!..............INTERMEDIATES........................

      integer,allocatable   :: ind(:)
      real,allocatable      ::  xn(:),xl(:),xj(:)
      integer,allocatable   ::  ilabel(:,:)
      integer i,j,m
      integer ns

!...............FINDING CLEBSCH GORDON ARRAY.....................

      real clebr
      integer ia,ja,k
      real xji,xjj,xmi,xmj,xka
      real xxn,xxl,xxj
      integer iphase

!..............FILE HANDLING..........................

      character spsfil*15               ! name of .sps file 
      integer ilast
      character isoread*3
      logical isoflag
      real yy                           ! dummy

!================================================================      


!      if(allocated(spsqn))deallocate(spsqn)
!      if(allocated(orb_qn))deallocate(orb_qn)
      if(allocated(jorb))deallocate(jorb)
      if(allocated(torb))deallocate(torb)
      if(allocated(nrorb))deallocate(nrorb)
!      if(allocated(jcore))deallocate(jcore)
!      if(allocated(tcore))deallocate(tcore)
!      if(allocated(nrcore))deallocate(nrcore)
!      if(allocated(spsqn))deallocate(spsqn)
!      if(allocated(kmax))deallocate(kmax)
!      if(allocated(kmin))deallocate(kmin)
!      if(allocated(indx))deallocate(indx)
!      if(allocated(clb))deallocate(clb)

1     continue
      write(6,*)' Enter shell-model space file name (.sps)'
      read(5,*)spsfil
      ilast=index(spsfil,' ')
      if(ilast.ne.0)then
         ilast=ilast-1
      else
         ilast=15
      endif
!----------  Now get orbit information for the calculation
!----------  njl is the number of single-particle states

      open(unit=1,file=spsfil(1:ilast)//'.sps',status='old',err=3)
      goto 4
3     continue
      write(6,*)' That file does not exist '
      goto 1

4     continue
      write(6,'('' Shell-Model space file name'',1x,a15)')spsfil

!.............................READS IN REDSTICK-COMPATIBLE FORMAT........

      read(1,'(a3)')isoread
      if(isoread.eq.'iso' .or. isoread.eq.'ISO')then
        isoflag=.true.

      elseif(isoread.eq.'pn' .or. isoread.eq.'PN')then
        isoflag=.false.
      else
        write(6,*)' You do not have a redstick compatible .sps file '
        write(6,*)' Nonetheless, I will try to read with old format '
        rewind(1)
        goto 1011
      endif


!....................ISOSPIN FORMALISM..........................

      if(isoflag)then

        read(1,*)norb(1)
        allocate(xn(norb(1)),xl(norb(1)),xj(norb(1)))
        do i = 1,norb(1)
        read(1,*,err=44)xn(i),xl(i),xj(i),yy
        enddo
        read(1,*)ncore(1)
        ntotal(1) = norb(1)+ncore(1)
         deallocate(xn,xj,xl)

        rewind(1)
    
        read(1,'(a3)')isoread
        read(1,*)norb(1)
        allocate(xn(norb(1)),xl(norb(1)),xj(norb(1)))
        allocate(nrorb(2*ntotal(1)),lorb(2*ntotal(1)),    &
                & jorb(2*ntotal(1)),torb(2*ntotal(1)),    &
                & nprincipal(2*ntotal(1)), nodal(2*ntotal(1)))
        do i = 1,norb(1)

        read(1,*,err=44)xn(i),xl(i),xj(i),yy

        nrorb(i)=int(xn(i))
        jorb(i)=int(2*xj(i))

        torb(i)=1
        lorb(i) = int(xl(i))
!Add by CFJ
        nprincipal(i) = 2*nrorb(i)+lorb(i) 
        nodal(i)      = nrorb(i)+1 
!................ERROR TRAPS ...............

        spinless = .false.
        if(xl(i) == xj(i))spinless = .true.

!..............END ERROR TRAPS.............
        
         enddo

        deallocate(xn,xj,xl)
!..............START TO READ CORE..........

        read(1,*)ncore(1)         
        allocate(xn(ncore(1)),xl(ncore(1)),xj(ncore(1)))
!        allocate(nrcore(2*ncore(1)),lcore(2*ncore(1)),        &
!                & jcore(2*ncore(1)),tcore(2*ncore(1)))
        do i = 1,ncore(1)
           j = norb(1)+i
        read(1,*,err=44)xn(i),xl(i),xj(i),yy

        nrorb(j)=int(xn(i))
        jorb(j) =int(2*xj(i))

        torb(j) =1
        lorb(j) = int(xl(i))
!Add by CFJ
        nprincipal(j) = 2*nrorb(j)+lorb(j)
        nodal(j)      = nrorb(j)+1

        enddo

         goto 35
   44    continue
         write(6,*)' some error in reading file ',i
         stop

      endif


!.......................PN FORMALISM.............................


        if(.not.isoflag)then
                write(6,*)' This formalism not yet implemented '
                goto 1
        endif

      return
!.............................THIS SECTION IS FOR OLD FORMAT.........

 1011 continue

      write(6,*)' This version assumes equal proton and neutron orbits '
      norb(1) = 0

      i=0
      do
        read(1,*,err=3033,end=34)i,xxn,xxl,xxj
        norb(1)=norb(1)+1
!................ERROR TRAPS ...............
        spinless = .false.
        if(xl(i) == xl(j))spinless = .true.


!..............END ERROR TRAPS.............
        
      enddo

 34   continue


      allocate(ind(norb(1)),xn(norb(1)),xl(norb(1)),xj(norb(1)))
      allocate(nrorb(2*norb(1)),lorb(2*norb(1)),jorb(2*norb(1)),  &
     &                                          torb(2*norb(1)))

      rewind(1)


      do i=1,norb(1)
        read(1,*,end=35,err=3033)ind(i),xn(i),xl(i),xj(i)
        write(6,*)i,ind(i),xn(i),xl(i),xj(i)
        nrorb(ind(i))=int(xn(i))
        jorb(ind(i))=int(2*xj(i))

        torb(ind(i))=1
        
        lorb(ind(i)) = int(xl(i))


!..............END ERROR TRAPS.............
        
      enddo

      deallocate(xn,xj,xl,ind)

   35 continue
      close(unit=1)

      norb(2) = norb(1)
      ncore(2)= ncore(1)
      ntotal(2)=ntotal(1)

      do i =1,ntotal(1)
        write(6,*)i,nrorb(i),jorb(i),lorb(i)

!-----------FILL NEUTRON SPACES --------------

        nrorb(i+ntotal(1))=nrorb(i)
        jorb(i+ntotal(1)) = jorb(i)
        lorb(i+ntotal(1)) = lorb(i)
        torb(i+ntotal(1)) = -1
        
      enddo

!!-------- FILL IN WEO ARRAYS
!
!      allocate(orb_qn(norb(1)+norb(2),3))
!
!      do i =1,norb(1)+norb(2)
!        orb_qn(i,1) = float(nrorb(i))
!        orb_qn(i,2) = float(lorb(i))
!        orb_qn(i,3) = float(jorb(i))/2.
!      enddo
!
!      ns = 0
!      do i=1,norb(1)
!         do m=-jorb(i),jorb(i),2
!           ns=ns+1
!         enddo
!      enddo
!
!      allocate(spsqn(2,ns,6))
!
!      j2max=0
!      ns=0
!      do i =1,norb(1)
!        j2max=max(j2max,jorb(i))
!        do m = -jorb(i),jorb(i),2
!            ns=ns+1
!            spsqn(1,ns,1) = i
!            spsqn(1,ns,2) = nrorb(i)
!            spsqn(1,ns,3) = lorb(i)
!            spsqn(1,ns,4) = jorb(i)
!            spsqn(1,ns,5) = m
!            spsqn(1,ns,6) = 1
!        enddo
!      enddo
!
!      ns=0
!      do i=norb(1)+1,norb(1)+norb(2)
!        do m=-jorb(i),jorb(i),2
!           ns=ns+1
!            spsqn(2,ns,1) = i
!            spsqn(2,ns,2) = nrorb(i)
!            spsqn(2,ns,3) = lorb(i)
!            spsqn(2,ns,4) = jorb(i)
!            spsqn(2,ns,5) = m
!            spsqn(2,ns,6) = -1
!        enddo
!      enddo
!
!      nsps(1) = ns
!      nsps(2) = ns
!
!
!      nindx=0
!weo---  Set up all possible two particle combinations along with spin
!weo---  Also find min and max angular momentum that can be coupled
!weo---  Do all like particles first

! CWJ -- NOTE NOTE I assume proton neutron have identical orbits --

!      allocate(kmax(norb(1),norb(1)),kmin(norb(1),norb(1)))
!
!      do i=1,norb(1)
!        do j=1,norb(1)
!            if(torb(i).eq.torb(j))then
!                  nindx=nindx+1
!            end if
!        end do
!      end do
!
!      allocate(indx(2,nindx))
!
!      nindx=0
!      allocate(ilabel(norb(1),norb(1)))
!
!      do i=1,norb(1)
!        do j=1,norb(1)
!            if(torb(i).eq.torb(j))then
!                  nindx=nindx+1
!                  indx(1,nindx)=i
!                  indx(2,nindx)=j
!                  ilabel(i,j) = nindx
!
!                  kmax(i,j)=int(xj(i)+xj(j))
!                  kmin(i,j)=int(abs(xj(i)-xj(j)))
!            end if
!        end do
!      end do
!Cweo---- Now different types  ---- Note xtz(i)=0. in isospin formalism
!C      do i=1,norb(1)+norb(2)
!C        do j=1,norb(1)+norb(2)
!C           if(torb(i).ne.torb(j))then
!C                  nindx=nindx+1
!C                  indx(1,nindx)=i
!C                  indx(2,nindx)=j
!C                 ilabel(i,j) = nindx
!C                  kmax(i,j)=int(xj(i)+xj(j))
!C                  kmin(i,j)=int(abs(xj(i)-xj(j)))
!C             end if
!C        end do
!C      end do
!C      write(6,*)'  nindx ',nindx
!C
!C---- compute phased Clebsch-Gordon array
!      allocate(clb(0:j2max,1:nsps(1),1:nsps(1)))
!      do i = 1,nsps(1)
!         ia = spsqn(1,i,1)
!         xji = float(spsqn(1,i,4))/2.
!         xmi = float(spsqn(1,i,5))/2.
!         do j = 1,nsps(1)
!            ja = spsqn(1,j,1)
!            xjj=float(spsqn(1,j,4))/2.
!            xmj=float(spsqn(1,j,5))/2.
!            iphase =  (spsqn(1,j,4)-spsqn(1,j,5))/2
!            do k = kmin(ia,ja),kmax(ia,ja)
!               xk = float(k)
!               clb(k,i,j)=(-1)**iphase*
!     &               clebr(xji,xmi,xjj,-xmj,xK,xmi-xmj)
!            enddo
!
!         enddo
!      enddo

      return
 3033 continue
      write(6,*)' I cannot understand this file. ', &
     & 'Please choose another '
      goto 1

      end

