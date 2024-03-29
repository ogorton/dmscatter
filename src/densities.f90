module densities
    use main
    use orbitals
    implicit none
    integer :: jt, tt
    integer :: maxJt=-1
    integer :: minJt=999
    logical :: pndens
    logical :: evenA
    
    contains
    !===================================================================
    subroutine openresults(resfile)
    
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
           print*,filename(1:ilast),'.dres does not exist '
    
       end do
    
       return
    end subroutine openresults
    
    
    !===================================================================
    subroutine setupdensities
    
        implicit none
    
        ! densities(J,iso,a,b)
        allocate(nuc_target%densitymats%rho( 0:10,0:1,1:ntotal(1),1:ntotal(1)) )
        print*,'Density matrix allocated, size:',size(nuc_target%densitymats%rho)
        nuc_target%densitymats%rho(:,:,:,:) = 0.0
    
    end subroutine setupdensities


    !===================================================================
    subroutine testValenceCount
      use norm
      implicit none
      integer :: a, jax2
      integer :: J0x2 ! ground state angular momentum
      integer :: Jt ! transition operator angular momentum
      real(kind=8) :: Zval, Nval
      real(kind=8) :: numberp, numbern

      Jt = 0
      J0x2 = nuc_target%groundstate%Jx2

      Zval = 0
      Nval = 0

      print '("Validating valence particle count against declared values")'
      if (pndens) then
        do a = 1, norb(1)
          jax2 = jorb(a)
          numberp = nuc_target%densitymats%rho(Jt,0,a,a) * q2norm(jax2) / q2norm(J0x2)
          Zval = Zval + numberp

          numbern = nuc_target%densitymats%rho(Jt,1,a,a) * q2norm(jax2) / q2norm(J0x2)
          Nval = Nval + numbern
          print '("a ",i3," jx2 ",i3," protons ",f10.4," neutrons ",f10.4)',a,jax2,numberp,numbern
        end do
      else
        print '("Test valence particle count not implemented for isospin densities.")'
      end if

      print '(" Protons: ",f10.4," Expected: ",i3)',Zval,nuc_target%Zval
      print '(" Neutrons: ",f10.4," Expected: ",i3)',Nval,nuc_target%Nval

      Zval = 0
      Nval = 0

      print '("Validating total particle count against declared values")'
      if (pndens) then
        do a = 1, ntotal(1)
          jax2 = jorb(a)
          numberp = nuc_target%densitymats%rho(Jt,0,a,a) * q2norm(jax2) / q2norm(J0x2)
          Zval = Zval + numberp

          numbern = nuc_target%densitymats%rho(Jt,1,a,a) * q2norm(jax2) / q2norm(J0x2)
          Nval = Nval + numbern
        end do
      else
        print '("Test valence particle count not implemented for isospin densities.")'
      end if

      print '(" Protons: ",f10.4," Expected: ",i3)',Zval,nuc_target%Z
      print '(" Neutrons: ",f10.4," Expected: ",i3)',Nval,nuc_target%N

      if (pndens.and.abs(Zval - nuc_target%Z)>1e-3) STOP "ERROR: Density matrix does not satisfy proton number definition"
      if (pndens.and.abs(Nval - nuc_target%N)>1e-3) STOP "ERROR: Density matrix does not satisfy neutron number definition"      

    end subroutine testValenceCount
    
    
    !===================================================================
    subroutine coredensity
    
      implicit none
    
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
    
    
    subroutine printdensities
        implicit none
        integer J,a, b
        print*,'Printing density matrix.'
        print*,"PN-format:",pndens
        print*,"min Jt:",minJt
        print*,"max Jt:",maxJt
        print*,'# spo =', ntotal(1)
        do J=minJt,maxJt
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
    
    
    subroutine readheaderv2(resfile)
       implicit none
       integer resfile
       character(23) :: tmpline
       integer i,j,n
       integer :: nlocalstates
       real(kind=8) e,ex,xj,xt
       integer np,zp

       rewind(resfile)
       read(resfile,'(a)')tmpline
       read(resfile,*)zp,np
       if( mod(np+zp,2) /= 0)then
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
                     nuc_target%groundstate%Jx2 = closest2J(evenA,xj)
                     nuc_target%groundstate%Tx2 = closest2J(evenA,xt)
                     print*,"evenA",evenA
                     print*,"Ground state Jx2",nuc_target%groundstate%Jx2
                     print*,"Grounod state Tx2",nuc_target%groundstate%Tx2
    
                   end if
    
                  if(n==j)then
                       nlocalstates = nlocalstates +1
                   else
                       exit
                   end if
                end do
    
    3           continue
    
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
       !integer m
       found = .false.
       select case (locchar)
         case ('i')
           do while(.not.found)
              read(resfile,'(a4)',end=111)myloc
              if(myloc(2:4)=='Ini')then
                  backspace(resfile)
                  read(resfile,'(a4,13x,i4)')myloc,n
                  found = .true.
                  exit
              endif
           end do
         case('f')
           read(resfile,'(a4,13x,i4)')myloc,n
           if(myloc(2:4)=='Fin')found=.true.
       end select
       finished = .false.
       return
    111 continue
       finished = .true.
       return
    end subroutine read2state
    
    
    subroutine read2Jtrans(resfile,found)
    
        implicit none
        integer resfile
        logical found
        character(3) :: tmpchar
     
        read(resfile,'(a3)',end=111)tmpchar
        if(tmpchar(2:3) == ' ' .or. tmpchar(2:3)=='In' .or. tmpchar(1:2)=='++')then
           found = .false.
           return
        endif
        if(tmpchar(2:3) == 'Jt')then
          backspace(resfile)
          read(resfile,'(5x,i4)')jt
     
          !update min/max jt
          minJt = min(minJt, jt)
          maxJt = max(maxJt, jt)
     
        end if
        found = .true.
     !..... ADDED IN VERSION 9 MAR 2017.... CHECK IF ISOSPIN OR PN FORMALISM...
        backspace(resfile)
        read(resfile,'(11x,a3)')tmpchar
     
        if(tmpchar=='pro')then
                pndens=.true.
                print*,'Density matrix is in explicit proton-neutron format.'
        else
                pndens=.false.
        end if
     
        return
     111 continue
        found=.false.
        return
    
    end subroutine read2Jtrans
    
    
    subroutine readdensity(resfile, success)
       implicit none
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
       end do
    
       return
    
    1  continue
       backspace(resfile)
       return
    end subroutine readdensity
    
    
    subroutine readalldensities(resfile)
       implicit none
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
    
              call readdensity(resfile, success)
              if(success)nodensities=.false.
          end do ! endoflist
          exit
    
       end do  ! endoffile
    
       if(nodensities)then
           print*,' Wait! That density file held no densities ! '
       end if
    
       return
    end subroutine readalldensities

end module densities
