module folium

    ! File object class
    ! Lets you open files given a file name character string. 
    ! Automatically finds an available unit number, or you can specify
    ! your own.
    !
    ! Example usage:
    !     type (foli) :: myfile
    !     real :: x, y
    !     myfile = foli("my_file_name.ext")
    !     call sopen(myfile)
    !     read(myfile%iunit,*) x, y
    !     call sclose(myfile)
    implicit none
    private 
    public :: foli, sopen, sclose, sopennew, sopenold
    logical :: lopen = .false.
    logical :: lexists = .false.
    logical :: lused = .false.

    type foli
        character(100) :: cfname
        integer :: iunit = 1
    end type

    interface foli
        procedure :: userfile
    end interface foli

contains

    subroutine sopenold(this)
        type(foli) :: this
        call sopen(this, new=.false.)
    end subroutine

    subroutine sopennew(this)
        type(foli) :: this
        call sopen(this, new=.true.)
    end subroutine    

    subroutine sopen(this,new)
        implicit none
        type(foli) :: this
        logical :: new

        print '("Attempting to open ",A)',trim(this%cfname)
        if (new) then
            lexists = .true.
        else
            inquire(file = trim(this%cfname), exist = lexists)
        end if
        inquire(unit = this%iunit, opened=lused)        

        do while (lused)
            print '("Warning: the requested unit ",I4," is already in use.")',&
                this%iunit
            print '("Incrementing unit number and trying again")'
            this%iunit = this%iunit + 1
            inquire(unit = this%iunit, opened=lused)
        end do

        if (lexists) then
            open(unit = this%iunit, file = trim(this%cfname))
            print '(A," has been opened.")', trim(this%cfname)
            lopen = .true.
        else
            print '("File does not exist.")'
        end if
    end subroutine sopen

    subroutine sclose(this)
        implicit none      
        type(foli), intent(in) :: this

        print '("Attempting to close ",A)',trim(this%cfname)
        inquire(unit = this%iunit, opened = lopen)

        if (lopen) then
            close(this%iunit)
            print '(A," has been closed.")',&
                trim(this%cfname)
        else
            print '(A," (unit",I4,") is already closed.")',&
                trim(this%cfname), this%iunit
        end if
    end subroutine sclose

    type(foli) function userfile()
        implicit none
        print '("Enter filename")'
        read*,userfile%cfname
    end function userfile


end module folium
