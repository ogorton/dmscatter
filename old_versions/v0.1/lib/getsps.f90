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
      read(5,'(a)')spsfil
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

end subroutine GetSPS

