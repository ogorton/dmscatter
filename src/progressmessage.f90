subroutine progressmessage(percent)
    ! Prints the current percent completion of a task.
    ! Use case: place this in a slow loop which otherwise doesn't print 
    ! progress information. This subroutine will continuously print the percent
    ! progress to standard output, but it will only use one line. This is 
    ! accomplished by printing backspace characters.
    ! The desired effect is only reached if no other print statements are used
    ! between calls to progressmessage.
    use iso_c_binding, only: c_backspace
    implicit none
    real :: percent
    character(50) :: message
    integer :: nback
    write(message, '(f5.1, a)') percent,'% complete.'
    nback = len(trim(message))
    write(6, '(a,a)', advance='no') repeat(c_backspace,nback), trim(message)
    if (percent==100.0) print*,""
    
end subroutine progressmessage
