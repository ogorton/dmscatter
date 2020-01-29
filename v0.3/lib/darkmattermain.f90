
!  ctrlchar = 'c'  count up # of states
!  ctrlchar = 'f'  fill up info on states
!  ctrlchar = 'p'  fill up parent reference states
!  ctrlchar = 'd'  fill up daughter reference states

!================================================

program darkmattermain
    use response
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
        function transition_probability(q,v,jchi,y,Mtiso)
            implicit none
            REAL(kind=8), INTENT(IN) :: q
            REAL(kind=8), INTENT(IN) :: v
            REAL(kind=8), INTENT(IN) :: jchi
            REAL(kind=8), INTENT(IN) :: y
            integer, intent(in) :: Mtiso
            real(kind=8) :: transition_probability
        end function
    end interface

    integer, parameter :: resfile = 33
    integer :: ap, an, Mtiso

    ! TEST BLOCK <<<<<<
    real(kind=8) :: transprob
    real(kind=8) :: q,v,jchi,y
    q=1.
    v=1.
    jchi=1.
    y=1.
    allocate(cvec(1)%c(17))
    allocate(cvec(0)%c(17))
    cvec(1)%c=1.0
    cvec(0)%c=1.0
    ! >>>>>>> END TEST BLOCK
    
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

    call test_nucresponse(Mtiso) 

    transprob = transition_probability(q,v,jchi,y,Mtiso)
   
    print*,transprob

end program  

