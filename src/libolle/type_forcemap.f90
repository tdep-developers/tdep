#include "precompilerdefinitions"
module type_forcemap
use konstanter, only: r8, i8, lo_iou, lo_huge, lo_hugeint, lo_tol, lo_sqtol, lo_status, lo_exitcode_blaslapack, &
                      lo_exitcode_param, lo_exitcode_baddim, lo_exitcode_symmetry
use gottochblandat, only: tochar, walltime, lo_sqnorm, lo_mean, lo_chop, lo_flattentensor, &
                          lo_unflatten_2tensor, lo_frobnorm, lo_invert3x3matrix, lo_enforce_linear_constraints, &
                          lo_progressbar_init, lo_progressbar, lo_return_unique, qsort, open_file, lo_determ
use type_crystalstructure, only: lo_crystalstructure
use type_symmetryoperation, only: lo_operate_on_vector, lo_spacegroup_operation
use type_blas_lapack_wrappers, only: lo_gemm, lo_gemv, lo_dgesvd
use type_forceconstant_firstorder, only: lo_forceconstant_firstorder
use type_forceconstant_secondorder, only: lo_forceconstant_secondorder
use type_forceconstant_thirdorder, only: lo_forceconstant_thirdorder
use type_forceconstant_fourthorder, only: lo_forceconstant_fourthorder
use lo_dielectric_interaction, only: lo_dielectric_tensor
use type_jij_secondorder, only: lo_jij_secondorder
use hdf5_wrappers, only: lo_h5_store_data, lo_h5_read_data, lo_h5_store_attribute, lo_h5_read_attribute, HID_T, &
                         H5F_ACC_TRUNC_F, H5F_ACC_RDONLY_F, h5open_f, h5fopen_f, h5fcreate_f, h5fclose_f, h5close_f, &
                         h5gopen_f, h5gcreate_f, h5gclose_f
use mpi_wrappers, only: lo_mpi_helper, lo_stop_gracefully
use lo_memtracker, only: lo_mem_helper
use lo_symmetry_of_interactions, only: lo_interaction_tensors, lo_tensor_shell, &
                                       lo_tensor_singletop, lo_tensor_pairop, lo_tensor_tripletop, lo_tensor_quartetop

implicit none

private
public :: lo_forcemap
! generate linear constraints
public :: lo_secondorder_rot_herm_huang
! generate coefficient matrices
public :: lo_coeffmatrix_singlet
public :: lo_coeffmatrix_pair
public :: lo_coeffmatrix_triplet
public :: lo_coeffmatrix_quartet
!public :: lo_coeffmatrix_magnetic_pair
!public :: lo_coeffmatrix_magnetic_longitudinal
!public :: lo_coeffmatrix_magnetic_crossterm_energy
!public :: lo_coeffmatrix_magnetic_crossterm_forces
public :: lo_coeffmatrix_unitcell_Z_singlet
public :: lo_coeffmatrix_supercell_Z_singlet
public :: lo_coeffmatrix_Z_pair
public :: lo_coeffmatrix_Z_triplet
public :: lo_coeffmatrix_eps_singlet
public :: lo_coeffmatrix_eps_pair

! What is the forcemap? To solve for the forceconstants in one go, we need the
! coefficient matrices for each tuplet, and the information on how to build a system
! of equations from those. Also it should contain the information on how to construct
! the linear constraints for rotational and translational invariance.
!     It also contains the information how to translate a set of irreducible
! forceconstants to usable forceconstants, as well as the capabilities to read the
! class from file and enforce the linear constraints for an arbitrary set of linear
! constraints. Schematically, this is how it is arranged:
!
!
!                              nfc (number of irreducible ifcs)
!                              relfc (which are the irreducible)
!                         relfcind (global index of the irreducible)
!                                      irreducible_shell
!                                              ▲
!                                              │
!                                              │      ┌─▶M (3x3 coefficient matrix)
!                                        ┌─────┴────┐ │
!                                    ┌──▶│ singlet  │─┘
!                                    │   └─────┬────┘    M (9x9 coefficient matrix)
!                         ┌───────┐  │   ┌─────┴────┐    i2 (index to atom 2 in the pair)
!                         │  uc   │──┼──▶│   pair   │───▶r (vector to atom 2)
!                         └───────┘  │   └─────┬────┘    lv (lattice vector to the unitcell of atom 2)
!                             ▲      │   ┌─────┴────┐
!                             │      └──▶│ triplet  │─┐
! ┌───────────┐   generate    │          └──────────┘ │  M (27x27 coefficient matrix)
! │lo_forcemap│───────────────┤                       │  i2 (index to atom 2 in the triplet)
! └───────────┘read_from_file │          ┌──────────┐ │  i3 (index to atom 3 in the triplet)
!       │                     │      ┌──▶│ singlet  │ └─▶v2 (vector to atom 2)
!       │                     ▼      │   └──────────┘    v3 (vector to atom 3)
!       │                 ┌───────┐  │   ┌──────────┐    lv2 (lattice vector to the unitcell of atom 2)
!       │                 │  ss   │──┼──▶│   pair   │    lv3 (lattice vector to the unitcell of atom 3)
!       │                 └───────┘  │   └──────────┘
!       │                            │   ┌──────────┐
!       │                            └──▶│ triplet  │
!       │                                └──────────┘
!       │
!       │                        eq_rot1 (first order rotational constraints)
!       │                        eq_rot2 (second order rotational constraints)
!       └──create_constraints───▶eq_rot3 (third order rotational constraints)
!                                eq_asr3 (third order translational constraints)
!

!> a general forcemap tuplet
type, abstract :: lo_forcemap_atom_tuplet
    !> which irreducible shell
    integer :: irreducible_shell = -lo_hugeint
    !> which operation transforms the irreducible here?
    integer :: operation_from_shell = -lo_hugeint
end type

!> forcemap singlet
type, extends(lo_forcemap_atom_tuplet) :: lo_forcemap_atom_singlet
end type

!> a forcemap pair
type, extends(lo_forcemap_atom_tuplet) :: lo_forcemap_atom_pair
    !> index to atom 1
    integer :: i1 = -lo_hugeint
    !> index to atom 2
    integer :: i2 = -lo_hugeint
    !> pair vector
    real(r8), dimension(3) :: r = lo_huge
    !> lattice vector to atom 2
    real(r8), dimension(3) :: lv = lo_huge
    !> lattice vectors to atom 2 in fractional coordinates
    real(r8), dimension(3) :: flv = lo_huge
    !> is it a self-term?
    logical :: selfterm = .false.
end type

!> a forcemap triplet
type, extends(lo_forcemap_atom_tuplet) :: lo_forcemap_atom_triplet
    !> index to atom number 2 and 3
    integer :: i1 = -lo_hugeint, i2 = -lo_hugeint, i3 = -lo_hugeint
    !> vector to atom number 2 and 3
    real(r8), dimension(3) :: v2 = lo_huge, v3 = lo_huge
    !> lattice vector to the atom 2 and 3
    real(r8), dimension(3) :: lv2 = lo_huge, lv3 = lo_huge
    !> fractional lattice vector to the atom 2 and 3
    real(r8), dimension(3) :: flv2 = lo_huge, flv3 = lo_huge
    !> is it a self-term?
    logical :: selfterm = .false.
end type

!> a forcemap quartet
type, extends(lo_forcemap_atom_tuplet) :: lo_forcemap_atom_quartet
    !> index to atom number 1,2,3,4
    integer :: i1 = -lo_hugeint, i2 = -lo_hugeint, i3 = -lo_hugeint, i4 = -lo_hugeint
    !> vector to atom number 2,3 and 4
    real(r8), dimension(3) :: v2 = lo_huge, v3 = lo_huge, v4 = lo_huge
    !> lattice vector to the atom 2,3 and 4
    real(r8), dimension(3) :: lv2 = lo_huge, lv3 = lo_huge, lv4 = lo_huge
    !> fractional lattice vector to the atom 2,3 and 4
    real(r8), dimension(3) :: flv2 = lo_huge, flv3 = lo_huge, flv4 = lo_huge
    !> is it a self-term?
    logical :: selfterm = .false.
end type

!> generic shell with interactions
type lo_forcemap_shell
    !> how many independent things do we have for this shell
    integer :: nx = -lo_hugeint
    !> in a global sense, what are the indices of these independent numbers?
    integer, dimension(:), allocatable :: ind_global
    !> in a local sense, which columns do thetaind represent?
    integer, dimension(:), allocatable :: ind_local
    !> the coefficient matrix
    real(r8), dimension(:, :), allocatable :: coeff
end type
type :: lo_forcemap_magpairshell
    !> how many independent things do we have for this shell
    integer :: ntheta_jij = -lo_hugeint, ntheta_tij = -lo_hugeint, ntheta_qij = -lo_hugeint
    !> in a global sense, what are the indices of these independent numbers?
    integer, dimension(:), allocatable :: thetaind_jij, thetaind_tij, thetaind_qij
    !> in a local sense, which columns do thetaind represent?
    integer, dimension(:), allocatable :: thetalocind_jij, thetalocind_tij, thetalocind_qij
    !> the coefficient matrix
    real(r8), dimension(:, :), allocatable :: coeff_jij, coeff_tij, coeff_qij
end type

!> encapsulate the constraints
type lo_forcemap_constraints
    !> total number of constraints
    integer :: nconstr_tot = -lo_hugeint
    !> number of "clean" constraints that don't mix orders
    integer :: neq1 = -lo_hugeint
    integer :: neq2 = -lo_hugeint
    integer :: neq3 = -lo_hugeint
    integer :: neq4 = -lo_hugeint
    !> C in Cx=D
    real(r8), dimension(:, :), allocatable :: eq1, eq2, eq3, eq4
    !> D in Cx=D
    real(r8), dimension(:), allocatable :: d1, d2, d3, d4

    !> number of Z triplet constraints
    integer :: neqz3 = -lo_hugeint
    real(r8), dimension(:, :), allocatable :: eqz3
    real(r8), dimension(:), allocatable :: dz
end type

!> pre-defined groups of tuplets to speed things up later
type lo_tupletgroup
    !> number of tuplets in this group
    integer :: n = -lo_hugeint
    !> indices to the tuplets of this group
    integer, dimension(:), allocatable :: ind
    !> does this group contain the self-term?
    logical :: contains_selfterm = .false.
end type

!> unitcell things
type lo_forcemap_unitcell
    ! force constant tuplets:
    integer :: n_fc_singlet = -lo_hugeint
    integer :: n_fc_pair = -lo_hugeint
    integer :: n_fc_triplet = -lo_hugeint
    integer :: n_fc_quartet = -lo_hugeint
    type(lo_forcemap_atom_singlet), dimension(:), allocatable :: fc_singlet
    type(lo_forcemap_atom_pair), dimension(:), allocatable :: fc_pair
    type(lo_forcemap_atom_triplet), dimension(:), allocatable :: fc_triplet
    type(lo_forcemap_atom_quartet), dimension(:), allocatable :: fc_quartet
    integer :: nx_fc_singlet = -lo_hugeint
    integer :: nx_fc_pair = -lo_hugeint
    integer :: nx_fc_triplet = -lo_hugeint
    integer :: nx_fc_quartet = -lo_hugeint
    real(r8), dimension(:), allocatable :: x_fc_singlet
    real(r8), dimension(:), allocatable :: x_fc_pair
    real(r8), dimension(:), allocatable :: x_fc_triplet
    real(r8), dimension(:), allocatable :: x_fc_quartet
    type(lo_tupletgroup), dimension(:), allocatable :: fc_triplet_group
    type(lo_tupletgroup), dimension(:), allocatable :: fc_quartet_group

    ! dielectric things
    integer :: n_eps_singlet = -lo_hugeint
    integer :: n_eps_pair = -lo_hugeint
    type(lo_forcemap_atom_singlet), dimension(:), allocatable :: eps_singlet
    type(lo_forcemap_atom_pair), dimension(:), allocatable :: eps_pair
    integer :: nx_eps_global = -lo_hugeint
    integer :: nx_eps_singlet = -lo_hugeint
    integer :: nx_eps_pair = -lo_hugeint
    real(r8), dimension(:), allocatable :: x_eps_global
    real(r8), dimension(:), allocatable :: x_eps_global_deriv
    real(r8), dimension(:), allocatable :: x_eps_singlet
    real(r8), dimension(:), allocatable :: x_eps_pair

    integer :: n_Z_singlet = -lo_hugeint
    integer :: n_Z_pair = -lo_hugeint
    integer :: n_Z_triplet = -lo_hugeint
    type(lo_forcemap_atom_singlet), dimension(:), allocatable :: Z_singlet
    type(lo_forcemap_atom_pair), dimension(:), allocatable :: Z_pair
    type(lo_forcemap_atom_triplet), dimension(:), allocatable :: Z_triplet
    integer :: nx_Z_singlet = -lo_hugeint
    integer :: nx_Z_pair = -lo_hugeint
    integer :: nx_Z_triplet = -lo_hugeint
    real(r8), dimension(:), allocatable :: x_Z_singlet
    real(r8), dimension(:), allocatable :: x_Z_pair
    real(r8), dimension(:), allocatable :: x_Z_triplet
    type(lo_tupletgroup), dimension(:), allocatable :: Z_triplet_group

    ! Magnetic things will go here
end type
!> supercell things
type lo_forcemap_supercell
    !> force constant tuplets:
    integer :: n_fc_singlet = -lo_hugeint
    integer :: n_fc_pair = -lo_hugeint
    integer :: n_fc_triplet = -lo_hugeint
    integer :: n_fc_quartet = -lo_hugeint
    integer, dimension(:, :), allocatable :: ind_fc_singlet
    integer, dimension(:, :), allocatable :: ind_fc_pair
    integer, dimension(:, :), allocatable :: ind_fc_triplet
    integer, dimension(:, :), allocatable :: ind_fc_quartet
    !> dielectric tuplers
    integer :: n_Z_singlet = -lo_hugeint
    integer :: n_Z_pair = -lo_hugeint
    integer :: n_Z_triplet = -lo_hugeint
    integer, dimension(:, :), allocatable :: ind_Z_singlet
    integer, dimension(:, :), allocatable :: ind_Z_pair
    integer, dimension(:, :), allocatable :: ind_Z_triplet
    integer :: n_eps_singlet = -lo_hugeint
    integer :: n_eps_pair = -lo_hugeint
    integer, dimension(:, :), allocatable :: ind_eps_singlet
    integer, dimension(:, :), allocatable :: ind_eps_pair
    ! magnetic stuff goes here
end type

!> map forces and stuff every which way
type lo_forcemap
    !> number of unitcell atoms
    integer :: n_atom_uc = -lo_hugeint
    !> number of supercell atoms
    integer :: n_atom_ss = -lo_hugeint
    !> unitcell things
    type(lo_forcemap_unitcell) :: xuc
    !> supercell things
    type(lo_forcemap_supercell) :: xss
    !> developer mode
    logical :: devmode

    !> forceconstant things
    logical :: have_fc_singlet = .false.
    logical :: have_fc_pair = .false.
    logical :: have_fc_triplet = .false.
    logical :: have_fc_quartet = .false.
    integer :: n_fc_singlet_shell = -lo_hugeint
    integer :: n_fc_pair_shell = -lo_hugeint
    integer :: n_fc_triplet_shell = -lo_hugeint
    integer :: n_fc_quartet_shell = -lo_hugeint
    type(lo_tensor_shell), dimension(:), allocatable :: fc_singlet_shell
    type(lo_tensor_shell), dimension(:), allocatable :: fc_pair_shell
    type(lo_tensor_shell), dimension(:), allocatable :: fc_triplet_shell
    type(lo_tensor_shell), dimension(:), allocatable :: fc_quartet_shell

    !> dielectric things
    integer :: polar = -lo_hugeint
    logical :: have_Z_pair = .false.
    logical :: have_Z_triplet = .false.
    logical :: have_eps_singlet = .false.
    logical :: have_eps_pair = .false.
    !> what type of polar correction
    integer :: polarcorrectiontype = -lo_hugeint
    integer :: n_Z_singlet_shell = -lo_hugeint
    integer :: n_Z_pair_shell = -lo_hugeint
    integer :: n_Z_triplet_shell = -lo_hugeint
    integer :: n_eps_singlet_shell = -lo_hugeint
    integer :: n_eps_pair_shell = -lo_hugeint
    type(lo_tensor_shell), dimension(:), allocatable :: Z_singlet_shell
    type(lo_tensor_shell), dimension(:), allocatable :: Z_pair_shell
    type(lo_tensor_shell), dimension(:), allocatable :: Z_triplet_shell
    type(lo_tensor_shell)                            :: eps_global_shell
    type(lo_tensor_shell), dimension(:), allocatable :: eps_singlet_shell
    type(lo_tensor_shell), dimension(:), allocatable :: eps_pair_shell

    !> magnetic things
    ! logical :: magnetic_pair_interactions=.false.
    ! logical :: magnetic_singlet_interactions=.false.
    ! integer :: nmagpairshells=-lo_hugeint
    ! integer :: nmagsingletshells=-lo_hugeint
    ! type(lo_forcemap_magpairshell), dimension(:), allocatable :: magpairshell
    ! type(lo_forcemap_shell),        dimension(:), allocatable :: magsingletshell

    !> Some operations
    type(lo_tensor_singletop), dimension(:), allocatable :: op_singlet
    type(lo_tensor_pairop),    dimension(:), allocatable :: op_pair
    type(lo_tensor_tripletop), dimension(:), allocatable :: op_triplet
    type(lo_tensor_quartetop), dimension(:), allocatable :: op_quartet

    !> Linear constraints to enforce proper invariances
    type(lo_forcemap_constraints) :: constraints
contains
    procedure :: generate
    procedure :: get_firstorder_forceconstant
    procedure :: get_secondorder_forceconstant
    procedure :: get_thirdorder_forceconstant
    procedure :: get_fourthorder_forceconstant
    procedure :: get_dielectric_tensors
    !procedure :: get_secondorder_jij
    procedure :: enforce_constraints
    procedure :: write_to_hdf5
    procedure :: read_from_hdf5
    procedure :: forceconstant_constraints
end type

! Interfaces to type_forcemap_generate
interface
    module subroutine generate(map, uc, ss, polarcorrectiontype, st, mw, mem, verbosity, devmode)
        class(lo_forcemap), intent(out) :: map
        type(lo_crystalstructure), intent(in) :: uc
        type(lo_crystalstructure), intent(in) :: ss
        integer, intent(in) :: polarcorrectiontype
        type(lo_interaction_tensors), intent(inout) :: st
        type(lo_mpi_helper), intent(inout) :: mw
        type(lo_mem_helper), intent(inout) :: mem
        integer, intent(in) :: verbosity
        logical, intent(in), optional :: devmode
    end subroutine
end interface
! Interfaces to type_forcemap_constraints
interface
    module subroutine forceconstant_constraints(map, uc, rotational, huanginvariances, hermitian, hermitian_rhs, huang_rhs, rotational_rhs, verbosity)
        class(lo_forcemap), intent(inout) :: map
        type(lo_crystalstructure), intent(in) :: uc
        logical, intent(in) :: rotational
        logical, intent(in) :: huanginvariances
        logical, intent(in) :: hermitian
        real(r8), dimension(:), intent(in) :: hermitian_rhs
        real(r8), dimension(:), intent(in) :: huang_rhs
        real(r8), dimension(:), intent(in) :: rotational_rhs
        integer, intent(in) :: verbosity
    end subroutine
    module subroutine enforce_constraints(map)
        class(lo_forcemap), intent(inout) :: map
    end subroutine
    module subroutine lo_secondorder_rot_herm_huang(map, uc, eq2, vD, neq2, rotational, huang, hermitian, hermitian_rhs, huang_rhs, rotational_rhs)
        class(lo_forcemap), intent(in) :: map
        type(lo_crystalstructure), intent(in) :: uc
        real(r8), dimension(:, :), allocatable, intent(out) :: eq2
        real(r8), dimension(:), allocatable, intent(out) :: vD
        integer, intent(out) :: neq2
        logical, intent(in) :: rotational, huang, hermitian
        real(r8), dimension(:), intent(in) :: hermitian_rhs
        real(r8), dimension(:), intent(in) :: huang_rhs
        real(r8), dimension(:), intent(in) :: rotational_rhs
    end subroutine
end interface
! Interfaces to the coefficient matrices
interface
    module subroutine lo_coeffmatrix_singlet(CM, map)
        real(r8), dimension(:, :), intent(out) :: CM
        type(lo_forcemap), intent(in) :: map
    end subroutine
    module subroutine lo_coeffmatrix_pair(UM, CM, map)
        real(r8), dimension(:, :), intent(in) :: UM
        real(r8), dimension(:, :), intent(out) :: CM
        type(lo_forcemap), intent(in) :: map
    end subroutine
    module subroutine lo_coeffmatrix_triplet(UM, CM, map)
        real(r8), dimension(:, :), intent(in) :: UM
        real(r8), dimension(:, :), intent(out) :: CM
        type(lo_forcemap), intent(in) :: map
    end subroutine
    module subroutine lo_coeffmatrix_quartet(UM, CM, map, relntheta)
        real(r8), dimension(:, :), intent(in) :: UM
        real(r8), dimension(:, :), intent(out) :: CM
        type(lo_forcemap), intent(in) :: map
        logical, dimension(81), intent(in) :: relntheta
    end subroutine
    ! module subroutine lo_coeffmatrix_magnetic_longitudinal(UM,CM,map)
    !     real(r8), dimension(:,:), intent(in) :: UM
    !     real(r8), dimension(:,:), intent(out) :: CM
    !     type(lo_forcemap), intent(in) :: map
    ! end subroutine
    ! module subroutine lo_coeffmatrix_magnetic_crossterm_energy(UM,MM,CM,map)
    !     real(r8), dimension(:,:), intent(in) :: UM
    !     real(r8), dimension(:,:), intent(in) :: MM
    !     real(r8), dimension(:), intent(out) :: CM
    !     type(lo_forcemap), intent(in) :: map
    ! end subroutine
    ! module subroutine lo_coeffmatrix_magnetic_crossterm_forces(UM,CM,map)
    !     real(r8), dimension(:,:), intent(in) :: UM
    !     real(r8), dimension(:,:), intent(out) :: CM
    !     type(lo_forcemap), intent(in) :: map
    ! end subroutine
    ! module subroutine lo_coeffmatrix_magnetic_pair(UM,CM,map)
    !     real(r8), dimension(:,:), intent(in) :: UM
    !     real(r8), dimension(:), intent(out) :: CM
    !     type(lo_forcemap), intent(in) :: map
    ! end subroutine
    module subroutine lo_coeffmatrix_unitcell_Z_singlet(map, CM)
        class(lo_forcemap), intent(in) :: map
        real(r8), dimension(:, :), intent(out) :: CM
    end subroutine
    module subroutine lo_coeffmatrix_supercell_Z_singlet(map, ss, CM)
        class(lo_forcemap), intent(in) :: map
        type(lo_crystalstructure), intent(in) :: ss
        real(r8), dimension(:, :), intent(out) :: CM
    end subroutine
    module subroutine lo_coeffmatrix_Z_pair(map, UM, CM)
        type(lo_forcemap), intent(in) :: map
        real(r8), dimension(:, :), intent(in) :: UM
        real(r8), dimension(:, :), intent(out) :: CM
    end subroutine
    module subroutine lo_coeffmatrix_Z_triplet(map, UM, CM)
        type(lo_forcemap), intent(in) :: map
        real(r8), dimension(:, :), intent(in) :: UM
        real(r8), dimension(:, :), intent(out) :: CM
    end subroutine
    module subroutine lo_coeffmatrix_eps_singlet(map, UM, CM)
        class(lo_forcemap), intent(in) :: map
        real(r8), dimension(:, :), intent(in) :: UM
        real(r8), dimension(:, :), intent(out) :: CM
    end subroutine
    module subroutine lo_coeffmatrix_eps_pair(map, UM, CM)
        type(lo_forcemap), intent(in) :: map
        real(r8), dimension(:, :), intent(in) :: UM
        real(r8), dimension(:, :), intent(out) :: CM
    end subroutine
end interface
! Interfaces to type_forcemap_io
interface
    module subroutine write_to_hdf5(map, filename, verbosity)
        class(lo_forcemap), intent(in) :: map
        character(len=*), intent(in) :: filename
        integer, intent(in) :: verbosity
    end subroutine
    module subroutine read_from_hdf5(map, uc, filename, verbosity)
        class(lo_forcemap), intent(out) :: map
        type(lo_crystalstructure), intent(in) :: uc
        character(len=*), intent(in) :: filename
        integer, intent(in) :: verbosity
    end subroutine
end interface
! Interfaces to type_forcemap_returntensors
interface
    module subroutine get_firstorder_forceconstant(map, p, fc)
        class(lo_forcemap), intent(in) :: map
        type(lo_crystalstructure), intent(in) :: p
        type(lo_forceconstant_firstorder), intent(out) :: fc
    end subroutine
    module subroutine get_secondorder_forceconstant(map, p, fc, mem, verbosity)
        class(lo_forcemap), intent(in) :: map
        type(lo_crystalstructure), intent(in) :: p
        type(lo_forceconstant_secondorder), intent(out) :: fc
        type(lo_mem_helper), intent(inout) :: mem
        integer, intent(in) :: verbosity
    end subroutine
    module subroutine get_thirdorder_forceconstant(map, p, fc)
        class(lo_forcemap), intent(in) :: map
        type(lo_crystalstructure), intent(in) :: p
        type(lo_forceconstant_thirdorder), intent(out) :: fc
    end subroutine
    module subroutine get_fourthorder_forceconstant(map, p, fc)
        class(lo_forcemap), intent(in) :: map
        type(lo_crystalstructure), intent(in) :: p
        type(lo_forceconstant_fourthorder), intent(out) :: fc
    end subroutine
    module subroutine get_dielectric_tensors(map, p, di)
        class(lo_forcemap), intent(in) :: map
        type(lo_crystalstructure), intent(in) :: p
        type(lo_dielectric_tensor), intent(out) :: di
    end subroutine
end interface

! A ton of interfaces, needed to make compilation bearable
interface
    module pure function coeff81x79(cn, u2, u3, u4) result(m)
        real(r8), dimension(81, 79), intent(in) :: cn
        real(r8), dimension(3), intent(in) :: u2, u3, u4
        real(r8), dimension(79, 3) :: m
        real(r8) :: u2x, u2y, u2z
        real(r8) :: u3x, u3y, u3z
        real(r8) :: u4x, u4y, u4z
    end function
    module pure function coeff81x80(cn, u2, u3, u4) result(m)
        real(r8), dimension(81, 80), intent(in) :: cn
        real(r8), dimension(3), intent(in) :: u2, u3, u4
        real(r8), dimension(80, 3) :: m
        real(r8) :: u2x, u2y, u2z
        real(r8) :: u3x, u3y, u3z
        real(r8) :: u4x, u4y, u4z
    end function
    module pure function coeff81x81(cn, u2, u3, u4) result(m)
        real(r8), dimension(81, 81), intent(in) :: cn
        real(r8), dimension(3), intent(in) :: u2, u3, u4
        real(r8), dimension(81, 3) :: m
        real(r8) :: u2x, u2y, u2z
        real(r8) :: u3x, u3y, u3z
        real(r8) :: u4x, u4y, u4z
    end function
    module pure function coeff81x76(cn, u2, u3, u4) result(m)
        real(r8), dimension(81, 76), intent(in) :: cn
        real(r8), dimension(3), intent(in) :: u2, u3, u4
        real(r8), dimension(76, 3) :: m
        real(r8) :: u2x, u2y, u2z
        real(r8) :: u3x, u3y, u3z
        real(r8) :: u4x, u4y, u4z
    end function
    module pure function coeff81x77(cn, u2, u3, u4) result(m)
        real(r8), dimension(81, 77), intent(in) :: cn
        real(r8), dimension(3), intent(in) :: u2, u3, u4
        real(r8), dimension(77, 3) :: m
        real(r8) :: u2x, u2y, u2z
        real(r8) :: u3x, u3y, u3z
        real(r8) :: u4x, u4y, u4z
    end function
    module pure function coeff81x78(cn, u2, u3, u4) result(m)
        real(r8), dimension(81, 78), intent(in) :: cn
        real(r8), dimension(3), intent(in) :: u2, u3, u4
        real(r8), dimension(78, 3) :: m
        real(r8) :: u2x, u2y, u2z
        real(r8) :: u3x, u3y, u3z
        real(r8) :: u4x, u4y, u4z
    end function
    module pure function coeff81x72(cn, u2, u3, u4) result(m)
        real(r8), dimension(81, 72), intent(in) :: cn
        real(r8), dimension(3), intent(in) :: u2, u3, u4
        real(r8), dimension(72, 3) :: m
        real(r8) :: u2x, u2y, u2z
        real(r8) :: u3x, u3y, u3z
        real(r8) :: u4x, u4y, u4z
    end function
    module pure function coeff81x73(cn, u2, u3, u4) result(m)
        real(r8), dimension(81, 73), intent(in) :: cn
        real(r8), dimension(3), intent(in) :: u2, u3, u4
        real(r8), dimension(73, 3) :: m
        real(r8) :: u2x, u2y, u2z
        real(r8) :: u3x, u3y, u3z
        real(r8) :: u4x, u4y, u4z
    end function
    module pure function coeff81x74(cn, u2, u3, u4) result(m)
        real(r8), dimension(81, 74), intent(in) :: cn
        real(r8), dimension(3), intent(in) :: u2, u3, u4
        real(r8), dimension(74, 3) :: m
        real(r8) :: u2x, u2y, u2z
        real(r8) :: u3x, u3y, u3z
        real(r8) :: u4x, u4y, u4z
    end function
    module pure function coeff81x75(cn, u2, u3, u4) result(m)
        real(r8), dimension(81, 75), intent(in) :: cn
        real(r8), dimension(3), intent(in) :: u2, u3, u4
        real(r8), dimension(75, 3) :: m
        real(r8) :: u2x, u2y, u2z
        real(r8) :: u3x, u3y, u3z
        real(r8) :: u4x, u4y, u4z
    end function
    module pure function coeff81x63(cn, u2, u3, u4) result(m)
        real(r8), dimension(81, 63), intent(in) :: cn
        real(r8), dimension(3), intent(in) :: u2, u3, u4
        real(r8), dimension(63, 3) :: m
        real(r8) :: u2x, u2y, u2z
        real(r8) :: u3x, u3y, u3z
        real(r8) :: u4x, u4y, u4z
    end function
    module pure function coeff81x64(cn, u2, u3, u4) result(m)
        real(r8), dimension(81, 64), intent(in) :: cn
        real(r8), dimension(3), intent(in) :: u2, u3, u4
        real(r8), dimension(64, 3) :: m
        real(r8) :: u2x, u2y, u2z
        real(r8) :: u3x, u3y, u3z
        real(r8) :: u4x, u4y, u4z
    end function
    module pure function coeff81x65(cn, u2, u3, u4) result(m)
        real(r8), dimension(81, 65), intent(in) :: cn
        real(r8), dimension(3), intent(in) :: u2, u3, u4
        real(r8), dimension(65, 3) :: m
        real(r8) :: u2x, u2y, u2z
        real(r8) :: u3x, u3y, u3z
        real(r8) :: u4x, u4y, u4z
    end function
    module pure function coeff81x66(cn, u2, u3, u4) result(m)
        real(r8), dimension(81, 66), intent(in) :: cn
        real(r8), dimension(3), intent(in) :: u2, u3, u4
        real(r8), dimension(66, 3) :: m
        real(r8) :: u2x, u2y, u2z
        real(r8) :: u3x, u3y, u3z
        real(r8) :: u4x, u4y, u4z
    end function
    module pure function coeff81x67(cn, u2, u3, u4) result(m)
        real(r8), dimension(81, 67), intent(in) :: cn
        real(r8), dimension(3), intent(in) :: u2, u3, u4
        real(r8), dimension(67, 3) :: m
        real(r8) :: u2x, u2y, u2z
        real(r8) :: u3x, u3y, u3z
        real(r8) :: u4x, u4y, u4z
    end function
    module pure function coeff81x68(cn, u2, u3, u4) result(m)
        real(r8), dimension(81, 68), intent(in) :: cn
        real(r8), dimension(3), intent(in) :: u2, u3, u4
        real(r8), dimension(68, 3) :: m
        real(r8) :: u2x, u2y, u2z
        real(r8) :: u3x, u3y, u3z
        real(r8) :: u4x, u4y, u4z
    end function
    module pure function coeff81x69(cn, u2, u3, u4) result(m)
        real(r8), dimension(81, 69), intent(in) :: cn
        real(r8), dimension(3), intent(in) :: u2, u3, u4
        real(r8), dimension(69, 3) :: m
        real(r8) :: u2x, u2y, u2z
        real(r8) :: u3x, u3y, u3z
        real(r8) :: u4x, u4y, u4z
    end function
    module pure function coeff81x70(cn, u2, u3, u4) result(m)
        real(r8), dimension(81, 70), intent(in) :: cn
        real(r8), dimension(3), intent(in) :: u2, u3, u4
        real(r8), dimension(70, 3) :: m
        real(r8) :: u2x, u2y, u2z
        real(r8) :: u3x, u3y, u3z
        real(r8) :: u4x, u4y, u4z
    end function
    module pure function coeff81x71(cn, u2, u3, u4) result(m)
        real(r8), dimension(81, 71), intent(in) :: cn
        real(r8), dimension(3), intent(in) :: u2, u3, u4
        real(r8), dimension(71, 3) :: m
        real(r8) :: u2x, u2y, u2z
        real(r8) :: u3x, u3y, u3z
        real(r8) :: u4x, u4y, u4z
    end function
    module pure function coeff81x52(cn, u2, u3, u4) result(m)
        real(r8), dimension(81, 52), intent(in) :: cn
        real(r8), dimension(3), intent(in) :: u2, u3, u4
        real(r8), dimension(52, 3) :: m
        real(r8) :: u2x, u2y, u2z
        real(r8) :: u3x, u3y, u3z
        real(r8) :: u4x, u4y, u4z
    end function
    module pure function coeff81x53(cn, u2, u3, u4) result(m)
        real(r8), dimension(81, 53), intent(in) :: cn
        real(r8), dimension(3), intent(in) :: u2, u3, u4
        real(r8), dimension(53, 3) :: m
        real(r8) :: u2x, u2y, u2z
        real(r8) :: u3x, u3y, u3z
        real(r8) :: u4x, u4y, u4z
    end function
    module pure function coeff81x54(cn, u2, u3, u4) result(m)
        real(r8), dimension(81, 54), intent(in) :: cn
        real(r8), dimension(3), intent(in) :: u2, u3, u4
        real(r8), dimension(54, 3) :: m
        real(r8) :: u2x, u2y, u2z
        real(r8) :: u3x, u3y, u3z
        real(r8) :: u4x, u4y, u4z
    end function
    module pure function coeff81x55(cn, u2, u3, u4) result(m)
        real(r8), dimension(81, 55), intent(in) :: cn
        real(r8), dimension(3), intent(in) :: u2, u3, u4
        real(r8), dimension(55, 3) :: m
        real(r8) :: u2x, u2y, u2z
        real(r8) :: u3x, u3y, u3z
        real(r8) :: u4x, u4y, u4z
    end function
    module pure function coeff81x56(cn, u2, u3, u4) result(m)
        real(r8), dimension(81, 56), intent(in) :: cn
        real(r8), dimension(3), intent(in) :: u2, u3, u4
        real(r8), dimension(56, 3) :: m
        real(r8) :: u2x, u2y, u2z
        real(r8) :: u3x, u3y, u3z
        real(r8) :: u4x, u4y, u4z
    end function
    module pure function coeff81x57(cn, u2, u3, u4) result(m)
        real(r8), dimension(81, 57), intent(in) :: cn
        real(r8), dimension(3), intent(in) :: u2, u3, u4
        real(r8), dimension(57, 3) :: m
        real(r8) :: u2x, u2y, u2z
        real(r8) :: u3x, u3y, u3z
        real(r8) :: u4x, u4y, u4z
    end function
    module pure function coeff81x58(cn, u2, u3, u4) result(m)
        real(r8), dimension(81, 58), intent(in) :: cn
        real(r8), dimension(3), intent(in) :: u2, u3, u4
        real(r8), dimension(58, 3) :: m
        real(r8) :: u2x, u2y, u2z
        real(r8) :: u3x, u3y, u3z
        real(r8) :: u4x, u4y, u4z
    end function
    module pure function coeff81x59(cn, u2, u3, u4) result(m)
        real(r8), dimension(81, 59), intent(in) :: cn
        real(r8), dimension(3), intent(in) :: u2, u3, u4
        real(r8), dimension(59, 3) :: m
        real(r8) :: u2x, u2y, u2z
        real(r8) :: u3x, u3y, u3z
        real(r8) :: u4x, u4y, u4z
    end function
    module pure function coeff81x60(cn, u2, u3, u4) result(m)
        real(r8), dimension(81, 60), intent(in) :: cn
        real(r8), dimension(3), intent(in) :: u2, u3, u4
        real(r8), dimension(60, 3) :: m
        real(r8) :: u2x, u2y, u2z
        real(r8) :: u3x, u3y, u3z
        real(r8) :: u4x, u4y, u4z
    end function
    module pure function coeff81x61(cn, u2, u3, u4) result(m)
        real(r8), dimension(81, 61), intent(in) :: cn
        real(r8), dimension(3), intent(in) :: u2, u3, u4
        real(r8), dimension(61, 3) :: m
        real(r8) :: u2x, u2y, u2z
        real(r8) :: u3x, u3y, u3z
        real(r8) :: u4x, u4y, u4z
    end function
    module pure function coeff81x62(cn, u2, u3, u4) result(m)
        real(r8), dimension(81, 62), intent(in) :: cn
        real(r8), dimension(3), intent(in) :: u2, u3, u4
        real(r8), dimension(62, 3) :: m
        real(r8) :: u2x, u2y, u2z
        real(r8) :: u3x, u3y, u3z
        real(r8) :: u4x, u4y, u4z
    end function
    module pure function coeff81x38(cn, u2, u3, u4) result(m)
        real(r8), dimension(81, 38), intent(in) :: cn
        real(r8), dimension(3), intent(in) :: u2, u3, u4
        real(r8), dimension(38, 3) :: m
        real(r8) :: u2x, u2y, u2z
        real(r8) :: u3x, u3y, u3z
        real(r8) :: u4x, u4y, u4z
    end function
    module pure function coeff81x39(cn, u2, u3, u4) result(m)
        real(r8), dimension(81, 39), intent(in) :: cn
        real(r8), dimension(3), intent(in) :: u2, u3, u4
        real(r8), dimension(39, 3) :: m
        real(r8) :: u2x, u2y, u2z
        real(r8) :: u3x, u3y, u3z
        real(r8) :: u4x, u4y, u4z
    end function
    module pure function coeff81x40(cn, u2, u3, u4) result(m)
        real(r8), dimension(81, 40), intent(in) :: cn
        real(r8), dimension(3), intent(in) :: u2, u3, u4
        real(r8), dimension(40, 3) :: m
        real(r8) :: u2x, u2y, u2z
        real(r8) :: u3x, u3y, u3z
        real(r8) :: u4x, u4y, u4z
    end function
    module pure function coeff81x41(cn, u2, u3, u4) result(m)
        real(r8), dimension(81, 41), intent(in) :: cn
        real(r8), dimension(3), intent(in) :: u2, u3, u4
        real(r8), dimension(41, 3) :: m
        real(r8) :: u2x, u2y, u2z
        real(r8) :: u3x, u3y, u3z
        real(r8) :: u4x, u4y, u4z
    end function
    module pure function coeff81x42(cn, u2, u3, u4) result(m)
        real(r8), dimension(81, 42), intent(in) :: cn
        real(r8), dimension(3), intent(in) :: u2, u3, u4
        real(r8), dimension(42, 3) :: m
        real(r8) :: u2x, u2y, u2z
        real(r8) :: u3x, u3y, u3z
        real(r8) :: u4x, u4y, u4z
    end function
    module pure function coeff81x43(cn, u2, u3, u4) result(m)
        real(r8), dimension(81, 43), intent(in) :: cn
        real(r8), dimension(3), intent(in) :: u2, u3, u4
        real(r8), dimension(43, 3) :: m
        real(r8) :: u2x, u2y, u2z
        real(r8) :: u3x, u3y, u3z
        real(r8) :: u4x, u4y, u4z
    end function
    module pure function coeff81x44(cn, u2, u3, u4) result(m)
        real(r8), dimension(81, 44), intent(in) :: cn
        real(r8), dimension(3), intent(in) :: u2, u3, u4
        real(r8), dimension(44, 3) :: m
        real(r8) :: u2x, u2y, u2z
        real(r8) :: u3x, u3y, u3z
        real(r8) :: u4x, u4y, u4z
    end function
    module pure function coeff81x45(cn, u2, u3, u4) result(m)
        real(r8), dimension(81, 45), intent(in) :: cn
        real(r8), dimension(3), intent(in) :: u2, u3, u4
        real(r8), dimension(45, 3) :: m
        real(r8) :: u2x, u2y, u2z
        real(r8) :: u3x, u3y, u3z
        real(r8) :: u4x, u4y, u4z
    end function
    module pure function coeff81x46(cn, u2, u3, u4) result(m)
        real(r8), dimension(81, 46), intent(in) :: cn
        real(r8), dimension(3), intent(in) :: u2, u3, u4
        real(r8), dimension(46, 3) :: m
        real(r8) :: u2x, u2y, u2z
        real(r8) :: u3x, u3y, u3z
        real(r8) :: u4x, u4y, u4z
    end function
    module pure function coeff81x47(cn, u2, u3, u4) result(m)
        real(r8), dimension(81, 47), intent(in) :: cn
        real(r8), dimension(3), intent(in) :: u2, u3, u4
        real(r8), dimension(47, 3) :: m
        real(r8) :: u2x, u2y, u2z
        real(r8) :: u3x, u3y, u3z
        real(r8) :: u4x, u4y, u4z
    end function
    module pure function coeff81x48(cn, u2, u3, u4) result(m)
        real(r8), dimension(81, 48), intent(in) :: cn
        real(r8), dimension(3), intent(in) :: u2, u3, u4
        real(r8), dimension(48, 3) :: m
        real(r8) :: u2x, u2y, u2z
        real(r8) :: u3x, u3y, u3z
        real(r8) :: u4x, u4y, u4z
    end function
    module pure function coeff81x49(cn, u2, u3, u4) result(m)
        real(r8), dimension(81, 49), intent(in) :: cn
        real(r8), dimension(3), intent(in) :: u2, u3, u4
        real(r8), dimension(49, 3) :: m
        real(r8) :: u2x, u2y, u2z
        real(r8) :: u3x, u3y, u3z
        real(r8) :: u4x, u4y, u4z
    end function
    module pure function coeff81x50(cn, u2, u3, u4) result(m)
        real(r8), dimension(81, 50), intent(in) :: cn
        real(r8), dimension(3), intent(in) :: u2, u3, u4
        real(r8), dimension(50, 3) :: m
        real(r8) :: u2x, u2y, u2z
        real(r8) :: u3x, u3y, u3z
        real(r8) :: u4x, u4y, u4z
    end function
    module pure function coeff81x51(cn, u2, u3, u4) result(m)
        real(r8), dimension(81, 51), intent(in) :: cn
        real(r8), dimension(3), intent(in) :: u2, u3, u4
        real(r8), dimension(51, 3) :: m
        real(r8) :: u2x, u2y, u2z
        real(r8) :: u3x, u3y, u3z
        real(r8) :: u4x, u4y, u4z
    end function
    module pure function coeff81x23(cn, u2, u3, u4) result(m)
        real(r8), dimension(81, 23), intent(in) :: cn
        real(r8), dimension(3), intent(in) :: u2, u3, u4
        real(r8), dimension(23, 3) :: m
        real(r8) :: u2x, u2y, u2z
        real(r8) :: u3x, u3y, u3z
        real(r8) :: u4x, u4y, u4z
    end function
    module pure function coeff81x24(cn, u2, u3, u4) result(m)
        real(r8), dimension(81, 24), intent(in) :: cn
        real(r8), dimension(3), intent(in) :: u2, u3, u4
        real(r8), dimension(24, 3) :: m
        real(r8) :: u2x, u2y, u2z
        real(r8) :: u3x, u3y, u3z
        real(r8) :: u4x, u4y, u4z
    end function
    module pure function coeff81x25(cn, u2, u3, u4) result(m)
        real(r8), dimension(81, 25), intent(in) :: cn
        real(r8), dimension(3), intent(in) :: u2, u3, u4
        real(r8), dimension(25, 3) :: m
        real(r8) :: u2x, u2y, u2z
        real(r8) :: u3x, u3y, u3z
        real(r8) :: u4x, u4y, u4z
    end function
    module pure function coeff81x26(cn, u2, u3, u4) result(m)
        real(r8), dimension(81, 26), intent(in) :: cn
        real(r8), dimension(3), intent(in) :: u2, u3, u4
        real(r8), dimension(26, 3) :: m
        real(r8) :: u2x, u2y, u2z
        real(r8) :: u3x, u3y, u3z
        real(r8) :: u4x, u4y, u4z
    end function
    module pure function coeff81x27(cn, u2, u3, u4) result(m)
        real(r8), dimension(81, 27), intent(in) :: cn
        real(r8), dimension(3), intent(in) :: u2, u3, u4
        real(r8), dimension(27, 3) :: m
        real(r8) :: u2x, u2y, u2z
        real(r8) :: u3x, u3y, u3z
        real(r8) :: u4x, u4y, u4z
    end function
    module pure function coeff81x28(cn, u2, u3, u4) result(m)
        real(r8), dimension(81, 28), intent(in) :: cn
        real(r8), dimension(3), intent(in) :: u2, u3, u4
        real(r8), dimension(28, 3) :: m
        real(r8) :: u2x, u2y, u2z
        real(r8) :: u3x, u3y, u3z
        real(r8) :: u4x, u4y, u4z
    end function
    module pure function coeff81x29(cn, u2, u3, u4) result(m)
        real(r8), dimension(81, 29), intent(in) :: cn
        real(r8), dimension(3), intent(in) :: u2, u3, u4
        real(r8), dimension(29, 3) :: m
        real(r8) :: u2x, u2y, u2z
        real(r8) :: u3x, u3y, u3z
        real(r8) :: u4x, u4y, u4z
    end function
    module pure function coeff81x30(cn, u2, u3, u4) result(m)
        real(r8), dimension(81, 30), intent(in) :: cn
        real(r8), dimension(3), intent(in) :: u2, u3, u4
        real(r8), dimension(30, 3) :: m
        real(r8) :: u2x, u2y, u2z
        real(r8) :: u3x, u3y, u3z
        real(r8) :: u4x, u4y, u4z
    end function
    module pure function coeff81x31(cn, u2, u3, u4) result(m)
        real(r8), dimension(81, 31), intent(in) :: cn
        real(r8), dimension(3), intent(in) :: u2, u3, u4
        real(r8), dimension(31, 3) :: m
        real(r8) :: u2x, u2y, u2z
        real(r8) :: u3x, u3y, u3z
        real(r8) :: u4x, u4y, u4z
    end function
    module pure function coeff81x32(cn, u2, u3, u4) result(m)
        real(r8), dimension(81, 32), intent(in) :: cn
        real(r8), dimension(3), intent(in) :: u2, u3, u4
        real(r8), dimension(32, 3) :: m
        real(r8) :: u2x, u2y, u2z
        real(r8) :: u3x, u3y, u3z
        real(r8) :: u4x, u4y, u4z
    end function
    module pure function coeff81x33(cn, u2, u3, u4) result(m)
        real(r8), dimension(81, 33), intent(in) :: cn
        real(r8), dimension(3), intent(in) :: u2, u3, u4
        real(r8), dimension(33, 3) :: m
        real(r8) :: u2x, u2y, u2z
        real(r8) :: u3x, u3y, u3z
        real(r8) :: u4x, u4y, u4z
    end function
    module pure function coeff81x34(cn, u2, u3, u4) result(m)
        real(r8), dimension(81, 34), intent(in) :: cn
        real(r8), dimension(3), intent(in) :: u2, u3, u4
        real(r8), dimension(34, 3) :: m
        real(r8) :: u2x, u2y, u2z
        real(r8) :: u3x, u3y, u3z
        real(r8) :: u4x, u4y, u4z
    end function
    module pure function coeff81x35(cn, u2, u3, u4) result(m)
        real(r8), dimension(81, 35), intent(in) :: cn
        real(r8), dimension(3), intent(in) :: u2, u3, u4
        real(r8), dimension(35, 3) :: m
        real(r8) :: u2x, u2y, u2z
        real(r8) :: u3x, u3y, u3z
        real(r8) :: u4x, u4y, u4z
    end function
    module pure function coeff81x36(cn, u2, u3, u4) result(m)
        real(r8), dimension(81, 36), intent(in) :: cn
        real(r8), dimension(3), intent(in) :: u2, u3, u4
        real(r8), dimension(36, 3) :: m
        real(r8) :: u2x, u2y, u2z
        real(r8) :: u3x, u3y, u3z
        real(r8) :: u4x, u4y, u4z
    end function
    module pure function coeff81x37(cn, u2, u3, u4) result(m)
        real(r8), dimension(81, 37), intent(in) :: cn
        real(r8), dimension(3), intent(in) :: u2, u3, u4
        real(r8), dimension(37, 3) :: m
        real(r8) :: u2x, u2y, u2z
        real(r8) :: u3x, u3y, u3z
        real(r8) :: u4x, u4y, u4z
    end function
    module pure function coeff81x1(cn, u2, u3, u4) result(m)
        real(r8), dimension(81, 1), intent(in) :: cn
        real(r8), dimension(3), intent(in) :: u2, u3, u4
        real(r8), dimension(1, 3) :: m
        real(r8) :: u2x, u2y, u2z
        real(r8) :: u3x, u3y, u3z
        real(r8) :: u4x, u4y, u4z
    end function
    module pure function coeff81x2(cn, u2, u3, u4) result(m)
        real(r8), dimension(81, 2), intent(in) :: cn
        real(r8), dimension(3), intent(in) :: u2, u3, u4
        real(r8), dimension(2, 3) :: m
        real(r8) :: u2x, u2y, u2z
        real(r8) :: u3x, u3y, u3z
        real(r8) :: u4x, u4y, u4z
    end function
    module pure function coeff81x3(cn, u2, u3, u4) result(m)
        real(r8), dimension(81, 3), intent(in) :: cn
        real(r8), dimension(3), intent(in) :: u2, u3, u4
        real(r8), dimension(3, 3) :: m
        real(r8) :: u2x, u2y, u2z
        real(r8) :: u3x, u3y, u3z
        real(r8) :: u4x, u4y, u4z
    end function
    module pure function coeff81x4(cn, u2, u3, u4) result(m)
        real(r8), dimension(81, 4), intent(in) :: cn
        real(r8), dimension(3), intent(in) :: u2, u3, u4
        real(r8), dimension(4, 3) :: m
        real(r8) :: u2x, u2y, u2z
        real(r8) :: u3x, u3y, u3z
        real(r8) :: u4x, u4y, u4z
    end function
    module pure function coeff81x5(cn, u2, u3, u4) result(m)
        real(r8), dimension(81, 5), intent(in) :: cn
        real(r8), dimension(3), intent(in) :: u2, u3, u4
        real(r8), dimension(5, 3) :: m
        real(r8) :: u2x, u2y, u2z
        real(r8) :: u3x, u3y, u3z
        real(r8) :: u4x, u4y, u4z
    end function
    module pure function coeff81x6(cn, u2, u3, u4) result(m)
        real(r8), dimension(81, 6), intent(in) :: cn
        real(r8), dimension(3), intent(in) :: u2, u3, u4
        real(r8), dimension(6, 3) :: m
        real(r8) :: u2x, u2y, u2z
        real(r8) :: u3x, u3y, u3z
        real(r8) :: u4x, u4y, u4z
    end function
    module pure function coeff81x7(cn, u2, u3, u4) result(m)
        real(r8), dimension(81, 7), intent(in) :: cn
        real(r8), dimension(3), intent(in) :: u2, u3, u4
        real(r8), dimension(7, 3) :: m
        real(r8) :: u2x, u2y, u2z
        real(r8) :: u3x, u3y, u3z
        real(r8) :: u4x, u4y, u4z
    end function
    module pure function coeff81x8(cn, u2, u3, u4) result(m)
        real(r8), dimension(81, 8), intent(in) :: cn
        real(r8), dimension(3), intent(in) :: u2, u3, u4
        real(r8), dimension(8, 3) :: m
        real(r8) :: u2x, u2y, u2z
        real(r8) :: u3x, u3y, u3z
        real(r8) :: u4x, u4y, u4z
    end function
    module pure function coeff81x9(cn, u2, u3, u4) result(m)
        real(r8), dimension(81, 9), intent(in) :: cn
        real(r8), dimension(3), intent(in) :: u2, u3, u4
        real(r8), dimension(9, 3) :: m
        real(r8) :: u2x, u2y, u2z
        real(r8) :: u3x, u3y, u3z
        real(r8) :: u4x, u4y, u4z
    end function
    module pure function coeff81x10(cn, u2, u3, u4) result(m)
        real(r8), dimension(81, 10), intent(in) :: cn
        real(r8), dimension(3), intent(in) :: u2, u3, u4
        real(r8), dimension(10, 3) :: m
        real(r8) :: u2x, u2y, u2z
        real(r8) :: u3x, u3y, u3z
        real(r8) :: u4x, u4y, u4z
    end function
    module pure function coeff81x11(cn, u2, u3, u4) result(m)
        real(r8), dimension(81, 11), intent(in) :: cn
        real(r8), dimension(3), intent(in) :: u2, u3, u4
        real(r8), dimension(11, 3) :: m
        real(r8) :: u2x, u2y, u2z
        real(r8) :: u3x, u3y, u3z
        real(r8) :: u4x, u4y, u4z
    end function
    module pure function coeff81x12(cn, u2, u3, u4) result(m)
        real(r8), dimension(81, 12), intent(in) :: cn
        real(r8), dimension(3), intent(in) :: u2, u3, u4
        real(r8), dimension(12, 3) :: m
        real(r8) :: u2x, u2y, u2z
        real(r8) :: u3x, u3y, u3z
        real(r8) :: u4x, u4y, u4z
    end function
    module pure function coeff81x13(cn, u2, u3, u4) result(m)
        real(r8), dimension(81, 13), intent(in) :: cn
        real(r8), dimension(3), intent(in) :: u2, u3, u4
        real(r8), dimension(13, 3) :: m
        real(r8) :: u2x, u2y, u2z
        real(r8) :: u3x, u3y, u3z
        real(r8) :: u4x, u4y, u4z
    end function
    module pure function coeff81x14(cn, u2, u3, u4) result(m)
        real(r8), dimension(81, 14), intent(in) :: cn
        real(r8), dimension(3), intent(in) :: u2, u3, u4
        real(r8), dimension(14, 3) :: m
        real(r8) :: u2x, u2y, u2z
        real(r8) :: u3x, u3y, u3z
        real(r8) :: u4x, u4y, u4z
    end function
    module pure function coeff81x15(cn, u2, u3, u4) result(m)
        real(r8), dimension(81, 15), intent(in) :: cn
        real(r8), dimension(3), intent(in) :: u2, u3, u4
        real(r8), dimension(15, 3) :: m
        real(r8) :: u2x, u2y, u2z
        real(r8) :: u3x, u3y, u3z
        real(r8) :: u4x, u4y, u4z
    end function
    module pure function coeff81x16(cn, u2, u3, u4) result(m)
        real(r8), dimension(81, 16), intent(in) :: cn
        real(r8), dimension(3), intent(in) :: u2, u3, u4
        real(r8), dimension(16, 3) :: m
        real(r8) :: u2x, u2y, u2z
        real(r8) :: u3x, u3y, u3z
        real(r8) :: u4x, u4y, u4z
    end function
    module pure function coeff81x17(cn, u2, u3, u4) result(m)
        real(r8), dimension(81, 17), intent(in) :: cn
        real(r8), dimension(3), intent(in) :: u2, u3, u4
        real(r8), dimension(17, 3) :: m
        real(r8) :: u2x, u2y, u2z
        real(r8) :: u3x, u3y, u3z
        real(r8) :: u4x, u4y, u4z
    end function
    module pure function coeff81x18(cn, u2, u3, u4) result(m)
        real(r8), dimension(81, 18), intent(in) :: cn
        real(r8), dimension(3), intent(in) :: u2, u3, u4
        real(r8), dimension(18, 3) :: m
        real(r8) :: u2x, u2y, u2z
        real(r8) :: u3x, u3y, u3z
        real(r8) :: u4x, u4y, u4z
    end function
    module pure function coeff81x19(cn, u2, u3, u4) result(m)
        real(r8), dimension(81, 19), intent(in) :: cn
        real(r8), dimension(3), intent(in) :: u2, u3, u4
        real(r8), dimension(19, 3) :: m
        real(r8) :: u2x, u2y, u2z
        real(r8) :: u3x, u3y, u3z
        real(r8) :: u4x, u4y, u4z
    end function
    module pure function coeff81x20(cn, u2, u3, u4) result(m)
        real(r8), dimension(81, 20), intent(in) :: cn
        real(r8), dimension(3), intent(in) :: u2, u3, u4
        real(r8), dimension(20, 3) :: m
        real(r8) :: u2x, u2y, u2z
        real(r8) :: u3x, u3y, u3z
        real(r8) :: u4x, u4y, u4z
    end function
    module pure function coeff81x21(cn, u2, u3, u4) result(m)
        real(r8), dimension(81, 21), intent(in) :: cn
        real(r8), dimension(3), intent(in) :: u2, u3, u4
        real(r8), dimension(21, 3) :: m
        real(r8) :: u2x, u2y, u2z
        real(r8) :: u3x, u3y, u3z
        real(r8) :: u4x, u4y, u4z
    end function
    module pure function coeff81x22(cn, u2, u3, u4) result(m)
        real(r8), dimension(81, 22), intent(in) :: cn
        real(r8), dimension(3), intent(in) :: u2, u3, u4
        real(r8), dimension(22, 3) :: m
        real(r8) :: u2x, u2y, u2z
        real(r8) :: u3x, u3y, u3z
        real(r8) :: u4x, u4y, u4z
    end function
end interface

contains

!> measure size in memory
function size_in_mem(map) result(mem)
    !> forcemap
    class(lo_forcemap), intent(inout) :: map
    !> size in memory, in bytes
    integer(i8) :: mem

end function

!> destroy
subroutine destroy(map)
    !> forcemap
    class(lo_forcemap), intent(inout) :: map

end subroutine

end module
