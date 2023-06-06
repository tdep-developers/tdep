submodule(type_forceconstant_firstorder) type_forceconstant_firstorder_io
use konstanter, only: lo_forceconstant_1st_HartreeBohr_to_eVA
use gottochblandat, only: open_file
implicit none
contains

!> write the forceconstant to file.
module subroutine writetofile(fc, p, fn)
    !> second order force constant
    class(lo_forceconstant_firstorder), intent(in) :: fc
    !> crystal structure
    type(lo_crystalstructure), intent(in) :: p
    !> filename
    character(len=*), intent(in) :: fn

    integer :: ii, uu

    ! Dump it
    uu = open_file('out', trim(fn))

    write (uu, "('# first order force constants for',I7,' atoms')") p%na

    do ii = 1, p%na
        associate (phi => fc%atom(ii)%m*lo_forceconstant_1st_HartreeBohr_to_eVA)
            write (uu, "(1X,3(1X,F20.15))") phi
        end associate
    end do

    close (uu)
end subroutine

end submodule
