#include "precompilerdefinitions"
module type_forceconstant_firstorder
use konstanter, only: flyt, lo_huge, lo_hugeint
use type_crystalstructure, only: lo_crystalstructure
implicit none

private
public :: lo_forceconstant_firstorder

type lo_fc1_atom
    !> index in the unit cell to the atom
    integer :: i1 = -lo_hugeint
    !> absolute vectors positioning the atoms
    real(flyt), dimension(3) :: v1 = lo_huge
    !> lattice vectors positioning the unit cell
    real(flyt), dimension(3) :: lv1 = lo_huge
    !> the force constant matrix
    real(flyt), dimension(3) :: m = lo_huge
end type

type lo_forceconstant_firstorder
    !> number of atoms
    integer :: na = -lo_hugeint
    !> info about each atom
    type(lo_fc1_atom), allocatable, dimension(:) :: atom

contains
    !> write to file
    procedure :: writetofile
end type

! Interfaces to type_forceconstant_firstorder_io
interface
    module subroutine writetofile(fc, p, fn)
        class(lo_forceconstant_firstorder), intent(in) :: fc
        type(lo_crystalstructure), intent(in) :: p
        character(len=*), intent(in) :: fn
    end subroutine
end interface

contains

end module
