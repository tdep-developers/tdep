module lo_symmetry_of_interactions
!!
!! Determine all the different types of tensorial interactions between ions. Not the
!! numerical values, but the analytical structure of those interactions -- a precursor
!! to solving it numerically.
!!
use konstanter, only: r8,lo_iou,lo_tol,lo_sqtol,lo_huge,lo_hugeint,lo_exitcode_symmetry,lo_exitcode_param,lo_bohr_to_A
use gottochblandat, only: walltime,tochar
use mpi_wrappers, only: lo_mpi_helper,lo_stop_gracefully
use lo_memtracker, only: lo_mem_helper
use type_crystalstructure, only: lo_crystalstructure
use type_distancetable, only: lo_distancetable
use type_voronoi_distancetable, only: lo_voronoi_distancetable

implicit none
private
public :: lo_interaction_tensors
public :: lo_tensor_shell
public :: lo_tensor_singletop
public :: lo_tensor_pairop
public :: lo_tensor_tripletop
public :: lo_tensor_quartetop

!> a generic tuplet
type, abstract :: lo_tuplet
    !> which unique tuplet is this
    integer :: unique=-lo_hugeint
    !> which operation brings the unique to this pair
    integer :: operation=-lo_hugeint
end type
    !> generic singlet
    type, extends(lo_tuplet) :: lo_tuplet_singlet
    end type
        !> singlet in the unitcell
        type, extends(lo_tuplet_singlet) :: lo_tuplet_ucsinglet
            !> unit cell index
            integer :: i1=-lo_hugeint
            !> index to representative singlet in the supercell
            integer :: ind_ss_singlet=-lo_hugeint
        end type
        !> singlet in the supercell
        type, extends(lo_tuplet_singlet) :: lo_tuplet_sssinglet
            !> unit cell index
            integer :: i1=-lo_hugeint
            !> supercell index
            integer :: j1=-lo_hugeint
        end type
        !> prototype singlet
        type, extends(lo_tuplet_ucsinglet) :: lo_tuplet_protsinglet
            !> how many atoms does this transform to
            integer :: n_unfold=-lo_hugeint
            !> which tuplet in does it unfold to
            integer, dimension(:), allocatable :: unfold_index
            !> which operation unfolds it
            integer, dimension(:), allocatable :: unfold_operation
            !> which operations leave this singlet invariant?
            integer, dimension(:), allocatable :: invariant_operation
        end type

    !> generic pair
    type, extends(lo_tuplet) :: lo_tuplet_pair
        !> index to atoms in the unitcell
        integer :: i1=-lo_hugeint,i2=-lo_hugeint
        !> index to unique atoms
        integer :: ui1=-lo_hugeint,ui2=-lo_hugeint
        !> vector to second atom
        real(r8), dimension(3) :: v=lo_huge
        !> operation that takes the unique shell to this pair
        !integer :: operation_from_shell_to_pair=-lo_hugeint
    end type
        !> pair in the unitcell
        type, extends(lo_tuplet_pair) :: lo_tuplet_ucpair
            !> index to representative pair in supercell
            integer :: ind_ss_pair=-lo_hugeint
        end type
        !> pair in the supercell
        type, extends(lo_tuplet_pair) :: lo_tuplet_sspair
            !> index to atoms in the supercell
            integer :: j1=-lo_hugeint,j2=-lo_hugeint
        end type
        !> prototype pair
        type, extends(lo_tuplet_ucpair) :: lo_tuplet_protpair
            !> how many atoms does this transform to
            integer :: n_unfold=-lo_hugeint
            !> which pairs does it unfold to
            integer, dimension(:), allocatable :: unfold_index
            !> which operation unfolds it
            integer, dimension(:), allocatable :: unfold_operation
            !> which operations leave this pair invariant?
            integer, dimension(:), allocatable :: invariant_operation
        end type

    !> generic triplet
    type, extends(lo_tuplet) :: lo_tuplet_triplet
        !> index to atoms
        integer :: i1=-lo_hugeint,i2=-lo_hugeint,i3=-lo_hugeint
        !> index to unique atoms
        integer :: ui1=-lo_hugeint,ui2=-lo_hugeint,ui3=-lo_hugeint
        !> vector to other atoms
        real(r8), dimension(3) :: v1=lo_huge,v2=lo_huge,v3=lo_huge
        !> operation that takes the unique shell to this pair
        integer :: operation_from_shell_to_triplet=-lo_hugeint
    end type
        !> pair in the unitcell
        type, extends(lo_tuplet_triplet) :: lo_tuplet_uctriplet
            !> indices to the supercell
            integer :: ind_ss_triplet=-lo_hugeint
        end type
        !> pair in the supercell
        type, extends(lo_tuplet_triplet) :: lo_tuplet_sstriplet
            !> indices in supercell
            integer :: j1=-lo_hugeint,j2=-lo_hugeint,j3=-lo_hugeint
        end type
        !> prototype triplet
        type, extends(lo_tuplet_uctriplet) :: lo_tuplet_prottriplet
            !> how many atoms does this transform to
            integer :: n_unfold=-lo_hugeint
            !> which triplets does it unfold to
            integer, dimension(:), allocatable :: unfold_index
            !> which operation unfolds it
            integer, dimension(:), allocatable :: unfold_operation
            !> which operations leave this triplet invariant?
            integer, dimension(:), allocatable :: invariant_operation
        end type

    !> generic quartet
    type, extends(lo_tuplet) :: lo_tuplet_quartet
        !> index to atoms
        integer :: i1=-lo_hugeint,i2=-lo_hugeint,i3=-lo_hugeint,i4=-lo_hugeint
        !> index to unique atoms
        integer :: ui1=-lo_hugeint,ui2=-lo_hugeint,ui3=-lo_hugeint,ui4=-lo_hugeint
        !> vector to other atoms
        real(r8), dimension(3) :: v1=lo_huge,v2=lo_huge,v3=lo_huge,v4=lo_huge
        !> operation that takes the unique shell to this pair
        integer :: operation_from_shell_to_quartet=-lo_hugeint
    end type
        !> pair in the unitcell
        type, extends(lo_tuplet_quartet) :: lo_tuplet_ucquartet
            !> indices to the supercell
            integer :: ind_ss_quartet=-lo_hugeint
        end type
        !> pair in the supercell
        type, extends(lo_tuplet_quartet) :: lo_tuplet_ssquartet
            !> indices in supercell
            integer :: j1=-lo_hugeint,j2=-lo_hugeint,j3=-lo_hugeint,j4=-lo_hugeint
        end type
        !> prototype triplet
        type, extends(lo_tuplet_ucquartet) :: lo_tuplet_protquartet
            !> how many atoms does this transform to
            integer :: n_unfold=-lo_hugeint
            !> which quartets does it unfold to
            integer, dimension(:), allocatable :: unfold_index
            !> which operation unfolds it
            integer, dimension(:), allocatable :: unfold_operation
            !> which operations leave this quartet invariant?
            integer, dimension(:), allocatable :: invariant_operation
        end type

! So, I should explain what I mean by expanded
! operations. Consider a normal symmetry operation
! on a rank 2 tensor, in matrix form it is
!
! ┌───────┐    ┌───────┐┌───────┐┌───────┐
! │       │    │       ││       ││       │
! │ Trot  │ =  │   R   ││   T   ││  R^T  │
! │       │    │       ││       ││       │
! └───────┘    └───────┘└───────┘└───────┘
!
! What I need is the same operation but as a single
! matrix multiplication. That can be done by
! representing the tensors as vectors:
!
!   ┌─┐   ┌──────────────┐┌─┐
!   │ │   │              ││ │
!   │ │   │              ││ │
!   │T│   │              ││ │
!   │r│ = │      R       ││T│
!   │o│   │              ││ │
!   │t│   │              ││ │
!   │ │   │              ││ │
!   │ │   │              ││ │
!   └─┘   └──────────────┘└─┘
!
! Which is just as general, and easier to work with.
! The types defined below are the representation in terms of large
! matrices.

!> container for operations on flattened tensors
type, abstract :: lo_tensor_tupletop
    !> which transposition, in the original notation
    integer :: permind=-lo_hugeint
    !> which operation, in original notation
    integer :: opind=-lo_hugeint
    !> determinant of the operation, for pseudovectors
    integer :: detop=-lo_hugeint
    !> forward operation in Cartesian coordinates
    real(r8), dimension(3,3) :: m3=-lo_huge
    !> inverse operation in Cartesian coordinates
    real(r8), dimension(3,3) :: im3=-lo_huge
    !> permutation array for the pure operation, not considering permutation
    integer, dimension(:), allocatable :: fmap
end type
    !> singlet operations
    type, extends(lo_tensor_tupletop) :: lo_tensor_singletop
    end type
    !> pair operations
    type, extends(lo_tensor_tupletop) :: lo_tensor_pairop
        !> Tr*Op matrix
        real(r8), dimension(9,9) :: sotr=lo_huge
        !> array that permutes indices
        integer, dimension(2) :: perm=-lo_hugeint
    end type
    !> triplet operations
    type, extends(lo_tensor_tupletop) :: lo_tensor_tripletop
        !> Tr*Op matrix
        real(r8), dimension(27,27) :: sotr=lo_huge
        !> array that permutes indices
        integer, dimension(3) :: perm=-lo_hugeint
    end type
    !> quartet operations
    type, extends(lo_tensor_tupletop) :: lo_tensor_quartetop
        !> Tr*Op matrix
        real(r8), dimension(81,81) :: sotr=lo_huge
        !> array that permutes indices
        integer, dimension(4) :: perm=-lo_hugeint
    end type

! Let me quickly explain what I mean by a shell.
! Any tuplet, like a pair or triplet, can be
! defined as some prototype + a transform:
!
!  ┌──┐   ┌───────┐┌──┐
!  │  │   │       ││  │
!  │a1│ = │  R1   ││a0│
!  │  │   │       ││  │
!  └──┘   └───────┘└──┘
!  ┌──┐   ┌───────┐┌──┐
!  │  │   │       ││  │
!  │a2│ = │  R2   ││a0│
!  │  │   │       ││  │
!  └──┘   └───────┘└──┘
!
! All the tuplets that share the same prototype a0,
! are defined as being in the same shell. So I only
! have to figure out the symmetry properties of the
! prototype of each shell, the rest can be
! constructed by transformations. That is what the
! indices 'operation_from_blabla' index in the tuplet
! types above, and the shells below define the
! prototypes.

!> generic shell
type lo_tensor_shell
    !> number of irreducible components
    integer :: nx
    !> coefficient matrix
    real(r8), dimension(:,:), allocatable :: coeff
    !> global index of the irreducible components
    integer, dimension(:), allocatable :: ind_global
    !> local index to undo the sparseness
    integer, dimension(:), allocatable :: ind_local
    !> how many tuplets does this shell unfold to
    integer :: n_unfold=-lo_hugeint
    !> indices to tuplets to unfold to
    integer, dimension(:), allocatable :: unfold_index
    !> operation that unfolds
    integer, dimension(:), allocatable :: unfold_operation
end type

!> The complete symmetry information. Details which pairs are equal to which pairs, and so on.
type lo_interaction_tensors
    !> symmetry operations in a variety of representations
    integer :: n_singlet_operation=-lo_hugeint
    integer :: n_pair_operation=-lo_hugeint
    integer :: n_triplet_operation=-lo_hugeint
    integer :: n_quartet_operation=-lo_hugeint
    type(lo_tensor_singletop), dimension(:), allocatable :: singletop
    type(lo_tensor_pairop),    dimension(:), allocatable :: pairop
    type(lo_tensor_tripletop), dimension(:), allocatable :: tripletop
    type(lo_tensor_quartetop), dimension(:), allocatable :: quartetop

    !> tuplets
    integer :: n_uc_pair=-lo_hugeint
    integer :: n_ss_pair=-lo_hugeint
    integer :: n_uc_triplet=-lo_hugeint
    integer :: n_ss_triplet=-lo_hugeint
    integer :: n_uc_quartet=-lo_hugeint
    integer :: n_ss_quartet=-lo_hugeint
    type(lo_tuplet_ucsinglet), dimension(:), allocatable :: uc_singlet
    type(lo_tuplet_sssinglet), dimension(:), allocatable :: ss_singlet
    type(lo_tuplet_ucpair),    dimension(:), allocatable :: uc_pair
    type(lo_tuplet_sspair),    dimension(:), allocatable :: ss_pair
    type(lo_tuplet_uctriplet), dimension(:), allocatable :: uc_triplet
    type(lo_tuplet_sstriplet), dimension(:), allocatable :: ss_triplet
    type(lo_tuplet_ucquartet), dimension(:), allocatable :: uc_quartet
    type(lo_tuplet_ssquartet), dimension(:), allocatable :: ss_quartet

    !> what to consider, in terms of general tuplets
    logical :: consider_pair   =.false.
    logical :: consider_triplet=.false.
    logical :: consider_quartet=.false.

    !> forceconstant shells
    logical :: have_fc_singlet    =.false.
    logical :: have_fc_pair       =.false.
    logical :: have_fc_triplet    =.false.
    logical :: have_fc_quartet    =.false.
    integer :: n_fc_singlet_shell =-lo_hugeint
    integer :: n_fc_pair_shell    =-lo_hugeint
    integer :: n_fc_triplet_shell =-lo_hugeint
    integer :: n_fc_quartet_shell =-lo_hugeint
    integer :: nx_fc_singlet      =-lo_hugeint
    integer :: nx_fc_pair         =-lo_hugeint
    integer :: nx_fc_triplet      =-lo_hugeint
    integer :: nx_fc_quartet      =-lo_hugeint
    type(lo_tensor_shell), dimension(:), allocatable :: fc_singlet_shell
    type(lo_tensor_shell), dimension(:), allocatable :: fc_pair_shell
    type(lo_tensor_shell), dimension(:), allocatable :: fc_triplet_shell
    type(lo_tensor_shell), dimension(:), allocatable :: fc_quartet_shell
    !> dielectric shells
    logical :: have_eps_global     =.false.
    logical :: have_eps_singlet    =.false.
    logical :: have_eps_pair       =.false.
    logical :: have_Z_singlet      =.false.
    logical :: have_Z_pair         =.false.
    logical :: have_Z_triplet      =.false.
    integer :: nx_eps_global       =-lo_hugeint
    integer :: nx_eps_singlet      =-lo_hugeint
    integer :: nx_eps_pair         =-lo_hugeint
    integer :: n_eps_singlet_shell =-lo_hugeint
    integer :: n_eps_pair_shell    =-lo_hugeint
    integer :: nx_Z_singlet        =-lo_hugeint
    integer :: nx_Z_pair           =-lo_hugeint
    integer :: nx_Z_triplet        =-lo_hugeint
    integer :: n_Z_singlet_shell   =-lo_hugeint
    integer :: n_Z_pair_shell      =-lo_hugeint
    integer :: n_Z_triplet_shell   =-lo_hugeint
    type(lo_tensor_shell)                            :: eps_global_shell
    type(lo_tensor_shell), dimension(:), allocatable :: eps_singlet_shell
    type(lo_tensor_shell), dimension(:), allocatable :: eps_pair_shell
    type(lo_tensor_shell), dimension(:), allocatable :: Z_singlet_shell
    type(lo_tensor_shell), dimension(:), allocatable :: Z_pair_shell
    type(lo_tensor_shell), dimension(:), allocatable :: Z_triplet_shell
    !> magnetic shells
    !logical :: have_mag_singlet=.false.
    !logical :: have_mag_pair   =.false.
    !logical :: have_mag_triplet=.false.
    !logical :: have_mag_quartet=.false.

    contains
        !> create the whole thing
        procedure :: generate
end type

!> Helper type with intermediate quantities, only exists temporarily
type lo_symtabhelper
    !> Distance table for the unit cell
    class(lo_distancetable), allocatable :: dtuc
    !> Distance table for the supercell
    class(lo_distancetable), allocatable :: dtss
    !> list of unique tuplets
    type(lo_tuplet_protsinglet), dimension(:), allocatable :: singlet
    type(lo_tuplet_protpair),    dimension(:), allocatable :: pair
    type(lo_tuplet_prottriplet), dimension(:), allocatable :: triplet
    type(lo_tuplet_protquartet), dimension(:), allocatable :: quartet

    !> count tuplets per atom
    integer, dimension(:), allocatable :: uc_ctr_pair
    integer, dimension(:), allocatable :: uc_ctr_triplet
    integer, dimension(:), allocatable :: uc_ctr_quartet
    integer, dimension(:), allocatable :: ss_ctr_pair
    integer, dimension(:), allocatable :: ss_ctr_triplet
    integer, dimension(:), allocatable :: ss_ctr_quartet

    !> offset in tuplet indices depending on atom1 indices
    integer, dimension(:), allocatable :: uc_offset_pair
    integer, dimension(:), allocatable :: uc_offset_triplet
    integer, dimension(:), allocatable :: uc_offset_quartet
    integer, dimension(:), allocatable :: ss_offset_pair
    integer, dimension(:), allocatable :: ss_offset_triplet
    integer, dimension(:), allocatable :: ss_offset_quartet
end type

interface
    module subroutine check_cutoff(ss,cutoff,tol,dt,wraparound,mw)
        type(lo_crystalstructure), intent(in) :: ss
        real(r8), intent(inout) :: cutoff
        real(r8), intent(in) :: tol
        type(lo_distancetable), intent(in) :: dt
        logical, intent(in) :: wraparound
        type(lo_mpi_helper), intent(inout) :: mw
    end subroutine
    module subroutine align_distance_tables(dtuc,dtss,ss,tol,mw,mem)
        class(lo_distancetable), intent(in) :: dtuc
        class(lo_distancetable), intent(inout) :: dtss
        type(lo_crystalstructure), intent(in) :: ss
        real(r8), intent(in) :: tol
        type(lo_mpi_helper), intent(inout) :: mw
        type(lo_mem_helper), intent(inout) :: mem
    end subroutine
    module subroutine expandoperations(sl,uc,spacegroup,transposition,mem)
        type(lo_interaction_tensors), intent(inout) :: sl
        type(lo_crystalstructure), intent(in) :: uc
        logical, intent(in) :: spacegroup
        logical, intent(in) :: transposition
        type(lo_mem_helper), intent(inout) :: mem
    end subroutine
    module pure function expandoperation_pair(o) result(bigo)
        real(r8), dimension(3,3), intent(in) :: o
        real(r8), dimension(9,9) :: bigo
    end function
    module pure function expandoperation_triplet(o) result(bigo)
        real(r8), dimension(3,3), intent(in) :: o
        real(r8), dimension(27,27) :: bigo
    end function
    module pure function expandoperation_quartet(o) result(bigo)
        real(r8), dimension(3,3), intent(in) :: o
        real(r8), dimension(81,81) :: bigo
    end function
end interface

interface
    module subroutine construct_tuplets(sl,uc,ss,sh,rc3,rc4,dc3,verbosity)
        class(lo_interaction_tensors), intent(inout) :: sl
        type(lo_crystalstructure), intent(in) :: uc,ss
        type(lo_symtabhelper), intent(inout) :: sh
        real(r8), intent(in) :: rc3,rc4,dc3
        integer, intent(in) :: verbosity
    end subroutine
    module subroutine tuplets_to_shells(sl,uc,ss,sh,mw,mem,verbosity)
        class(lo_interaction_tensors), intent(inout) :: sl
        type(lo_crystalstructure), intent(in) :: uc,ss
        type(lo_symtabhelper), intent(inout) :: sh
        type(lo_mpi_helper), intent(inout) :: mw
        type(lo_mem_helper), intent(inout) :: mem
        integer, intent(in) :: verbosity
    end subroutine
    module function compare_tuplets_after_one_is_rotated(sl,sh,ss,t1,t2,op,tol) result(match)
        type(lo_interaction_tensors), intent(in) :: sl
        type(lo_symtabhelper), intent(in) :: sh
        type(lo_crystalstructure), intent(in) :: ss
        class(lo_tuplet), intent(in) :: t1,t2
        class(lo_tensor_tupletop), intent(in) :: op
        real(r8), intent(in) :: tol
        logical :: match
    end function
end interface

interface
    module subroutine nullspace_fc_singlet(sl,sh,uc,mw,mem)
        type(lo_interaction_tensors), intent(inout) :: sl
        type(lo_symtabhelper), intent(in) :: sh
        type(lo_crystalstructure), intent(in) :: uc
        type(lo_mpi_helper), intent(inout) :: mw
        type(lo_mem_helper), intent(inout) :: mem
    end subroutine
    module subroutine nullspace_fc_pair(sl,sh,ss,mw,mem)
        type(lo_interaction_tensors), intent(inout) :: sl
        type(lo_symtabhelper), intent(in) :: sh
        type(lo_crystalstructure), intent(in) :: ss
        type(lo_mpi_helper), intent(inout) :: mw
        type(lo_mem_helper), intent(inout) :: mem
    end subroutine
    module subroutine nullspace_fc_triplet(sl,sh,ss,mw,mem)
        type(lo_interaction_tensors), intent(inout) :: sl
        type(lo_symtabhelper), intent(in) :: sh
        type(lo_crystalstructure), intent(in) :: ss
        type(lo_mpi_helper), intent(inout) :: mw
        type(lo_mem_helper), intent(inout) :: mem
    end subroutine
    module subroutine nullspace_fc_quartet(sl,sh,ss,mw,mem)
        type(lo_interaction_tensors), intent(inout) :: sl
        type(lo_symtabhelper), intent(in) :: sh
        type(lo_crystalstructure), intent(in) :: ss
        type(lo_mpi_helper), intent(inout) :: mw
        type(lo_mem_helper), intent(inout) :: mem
    end subroutine
    module subroutine nullspace_eps_global(sl,mw,mem)
        type(lo_interaction_tensors), intent(inout) :: sl
        type(lo_mpi_helper), intent(inout) :: mw
        type(lo_mem_helper), intent(inout) :: mem
    end subroutine
    module subroutine nullspace_eps_singlet(sl,sh,uc,mw,mem)
        type(lo_interaction_tensors), intent(inout) :: sl
        type(lo_symtabhelper), intent(in) :: sh
        type(lo_crystalstructure), intent(in) :: uc
        type(lo_mpi_helper), intent(inout) :: mw
        type(lo_mem_helper), intent(inout) :: mem
    end subroutine
    module subroutine nullspace_eps_pair(sl,sh,uc,ss,cutoff,mw,mem)
        type(lo_interaction_tensors), intent(inout) :: sl
        type(lo_symtabhelper), intent(in) :: sh
        type(lo_crystalstructure), intent(in) :: uc
        type(lo_crystalstructure), intent(in) :: ss
        real(r8), intent(in) :: cutoff
        type(lo_mpi_helper), intent(inout) :: mw
        type(lo_mem_helper), intent(inout) :: mem
    end subroutine
    module subroutine nullspace_Z_singlet(sl,sh,uc,mw,mem)
        type(lo_interaction_tensors), intent(inout) :: sl
        type(lo_symtabhelper), intent(in) :: sh
        type(lo_crystalstructure), intent(in) :: uc
        type(lo_mpi_helper), intent(inout) :: mw
        type(lo_mem_helper), intent(inout) :: mem
    end subroutine
    module subroutine nullspace_Z_pair(sl,sh,uc,ss,cutoff,mw,mem)
        type(lo_interaction_tensors), intent(inout) :: sl
        type(lo_symtabhelper), intent(in) :: sh
        type(lo_crystalstructure), intent(in) :: uc
        type(lo_crystalstructure), intent(in) :: ss
        real(r8), intent(in) :: cutoff
        type(lo_mpi_helper), intent(inout) :: mw
        type(lo_mem_helper), intent(inout) :: mem
    end subroutine
    module subroutine nullspace_Z_triplet(sl,sh,uc,ss,cutoff,mw,mem)
        type(lo_interaction_tensors), intent(inout) :: sl
        type(lo_symtabhelper), intent(in) :: sh
        type(lo_crystalstructure), intent(in) :: uc
        type(lo_crystalstructure), intent(in) :: ss
        real(r8), intent(in) :: cutoff
        type(lo_mpi_helper), intent(inout) :: mw
        type(lo_mem_helper), intent(inout) :: mem
    end subroutine
    ! module subroutine nullspace_magnetic_pair(sl,sh,uc,cutoff,mw,mem)
    !     type(lo_interaction_tensors), intent(inout) :: sl
    !     type(lo_symtabhelper), intent(in) :: sh
    !     type(lo_crystalstructure), intent(in) :: uc
    !     real(r8), intent(in) :: cutoff
    !     type(lo_mpi_helper), intent(inout) :: mw
    !     type(lo_mem_helper), intent(inout) :: mem
    ! end subroutine
end interface

contains

!> create the symmetry list
subroutine generate(sl,uc,ss,cutoff2,cutoff3,cutoff4,polar,mw,mem,verbosity,&
                    transposition,spacegroup,&
                    firstorder,wzdim,nj2,nj3,nj4,tol,&
                    wraparound,magcutoff2,magsinglet,&
                    dielcutoff2,dielcutoff3)
    !> full symmetry table thing
    class(lo_interaction_tensors), intent(out) :: sl
    !> unitcell crystal structure
    type(lo_crystalstructure), intent(inout) :: uc
    !> supercell crystal structure
    type(lo_crystalstructure), intent(inout) :: ss
    !> cutoff for second order
    real(r8), intent(in) :: cutoff2
    !> cutoff for third order
    real(r8), intent(in) :: cutoff3
    !> cutoff for fourth order
    real(r8), intent(in) :: cutoff4
    !> is it a polar material
    logical, intent(in) :: polar
    !> mpi helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem
    !> how much to talk
    integer, intent(in) :: verbosity
    !> generate transposition symmetries
    logical, intent(in), optional :: transposition
    !> generate spacegroup symmetries
    logical, intent(in), optional :: spacegroup
    !> get firstorder forceconstants
    logical, intent(in), optional :: firstorder
    !> WZ supercell dimensions
    integer, dimension(3), optional, intent(in) :: wzdim
    !> graph cutoff
    integer, optional, intent(in) :: nj2,nj3,nj4
    !> tolerance
    real(r8), intent(in), optional :: tol
    !> allow for really long interaction thingies
    logical, intent(in), optional :: wraparound
    !> pair magnetic interactions
    real(r8), intent(in), optional :: magcutoff2
    !> magnetic on-site term
    logical, intent(in), optional :: magsinglet
    !> dielectric expansion cutoff for pairs
    real(r8), intent(in), optional :: dielcutoff2
    !> dielectric expansion cutoff for triplets
    real(r8), intent(in), optional :: dielcutoff3

    type(lo_symtabhelper) :: sh
    real(r8) :: rc2,rc3,rc4,dc2,dc3,mc2
    real(r8) :: t0,t1,timer

    ! Start timers
    call mem%tick()
    timer=walltime()
    t0=timer
    t1=timer

    ! Try to kill conflicting options
    init: block
        if ( verbosity .gt. 0 ) then
            write(lo_iou,*) ''
            write(lo_iou,*) 'DETERMINING SYMMETRY OF INTERACTION TENSORS'
        endif

        ! Make sure we know how the supercell and unitcell are related
        call ss%classify('supercell',uc)
    end block init

    ! Generate the distance tables, easier said than done.
    distancetable: block
        type(lo_distancetable) :: dt
        real(r8) :: f0,f1
        integer, dimension(:), allocatable :: nbrctr
        integer, dimension(3) :: wzd
        integer :: itr,i,j,l
        logical :: wrap,wzcutoff,groupsym,transsym

        ! Pick up the cutoffs from the input
        rc2=cutoff2
        rc3=cutoff3
        rc4=cutoff4
        if ( present(magcutoff2) ) then
            mc2=magcutoff2
        else
            mc2=-1.0_r8
        endif
        if ( present(dielcutoff2) ) then
            dc2=dielcutoff2
        else
            dc2=-1.0_r8
        endif
        if ( present(dielcutoff3) ) then
            dc3=dielcutoff3
        else
            dc3=-1.0_r8
        endif

        ! Perhaps I specify number of neighbours instead of a cutoff.
        if ( present(nj2) .or. present(nj3) .or. present(nj4) ) then
            ! get a distance table that is as long as possible
            call dt%generate(uc%r,uc%latticevectors,cutoff=ss%maxcutoff(),mw=mw,verbosity=0)
            call mem%allocate(nbrctr,uc%na,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
            nbrctr=0
            ! increase the cutoff until I have the desired number of neighbours within it
            if ( present(nj2) ) then; if ( nj2 .gt. 0 ) then
                f0=ss%mincutoff()   ! shortest possible cutoff
                f1=f0*0.01_r8     ! step to increase the cutoff with
                do itr=1,1000000
                    f0=f0+f1
                    call dt%count_neighbours_within_cutoff(f0,nbrctr)
                    l=minval(nbrctr) ! smallest number of neighbours
                    if ( l .ge. nj2 ) exit
                    if ( f0 .gt. ss%maxcutoff() ) exit
                enddo
                rc2=min(f0,ss%maxcutoff())
            endif; endif
            if ( present(nj3) ) then; if ( nj3 .gt. 0 ) then;
                f0=ss%mincutoff()   ! shortest possible cutoff
                f1=f0*0.01_r8     ! step to increase the cutoff with
                do itr=1,1000000
                    f0=f0+f1
                    call dt%count_neighbours_within_cutoff(f0,nbrctr)
                    l=minval(nbrctr) ! smallest number of neighbours
                    if ( l .ge. nj3 ) exit
                    if ( f0 .gt. ss%maxcutoff() ) exit
                enddo
                rc3=min(f0,ss%maxcutoff())
            endif; endif
            if ( present(nj4) ) then; if ( nj4 .gt. 0 ) then
                f0=ss%mincutoff()   ! shortest possible cutoff
                f1=f0*0.01_r8     ! step to increase the cutoff with
                do itr=1,1000000
                    f0=f0+f1
                    call dt%count_neighbours_within_cutoff(f0,nbrctr)
                    l=minval(nbrctr) ! smallest number of neighbours
                    if ( l .ge. nj4 ) exit
                    if ( f0 .gt. ss%maxcutoff() ) exit
                enddo
                rc4=min(f0,ss%maxcutoff())
            endif; endif
            call mem%deallocate(nbrctr,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        endif

        ! Build a too large distance table to check the cutoffs
        if ( present(wraparound) ) then
            wrap=wraparound
        else
            wrap=.false.
        endif
        if ( wrap ) then
            f0=maxval([rc2,rc3,rc3,mc2,dc2,dc3])+0.5_r8
        else
            f0=ss%maxcutoff()+0.5_r8
        endif
        call dt%generate(uc%r,uc%latticevectors,cutoff=f0,mw=mw,verbosity=0)
        ! Now actually check the cutoffs
        call check_cutoff(ss,rc2,2*lo_tol,dt,wrap,mw)
        call check_cutoff(ss,rc3,2*lo_tol,dt,wrap,mw)
        call check_cutoff(ss,rc4,2*lo_tol,dt,wrap,mw)
        call check_cutoff(ss,mc2,2*lo_tol,dt,wrap,mw)
        call check_cutoff(ss,dc2,2*lo_tol,dt,wrap,mw)
        call check_cutoff(ss,dc3,2*lo_tol,dt,wrap,mw)

        ! Use Wigner-Seitz cutoff?
        if ( present(wzdim) ) then
            if ( wzdim(1) .gt. 0 ) then
                wzd=wzdim
                wzcutoff=.true.
            elseif ( wzdim(1) .eq. -2 ) then
                wzcutoff=.true.
                wzd=-2
            else
                wzd=-1
                wzcutoff=.false.
            endif
        else
            wzd=-1
            wzcutoff=.false.
        endif

        ! generate distance tables for the unit cell and supercell
        if ( wzcutoff ) then
            ! Use wigner-seitz cutoff for the second order
            allocate(lo_voronoi_distancetable::sh%dtuc)
            allocate(lo_voronoi_distancetable::sh%dtss)
            select type(a=>sh%dtuc); type is(lo_voronoi_distancetable)
                if ( wzd(1) .eq. -2 ) then
                    ! Secret option don't bother thinking about it
                    call a%generate_wzcutoff(uc%r,uc%latticevectors,[1,1,1],ss%latticevectors)
                else
                    call a%generate_wzcutoff(uc%r,uc%latticevectors,wzd,uc%latticevectors,vorona=uc%na)
                endif
                ! some sanity checks for this. It will get weird if the WZ cutoff is smaller than
                ! any of the other cutoffs, and makes no sense.
                if ( a%cutoff .lt. rc3 ) call lo_stop_gracefully(['WZ cutoff smaller than third order'],lo_exitcode_symmetry,__FILE__,__LINE__,mw%comm)
                if ( a%cutoff .lt. rc4 ) call lo_stop_gracefully(['WZ cutoff smaller than fourth order'],lo_exitcode_symmetry,__FILE__,__LINE__,mw%comm)
                if ( a%cutoff .lt. mc2 ) call lo_stop_gracefully(['WZ cutoff smaller than magnetic pair'],lo_exitcode_symmetry,__FILE__,__LINE__,mw%comm)
                if ( a%cutoff .lt. dc2 ) call lo_stop_gracefully(['WZ cutoff smaller than dielectric pair'],lo_exitcode_symmetry,__FILE__,__LINE__,mw%comm)
                if ( a%cutoff .lt. dc3 ) call lo_stop_gracefully(['WZ cutoff smaller than dielectric triplet'],lo_exitcode_symmetry,__FILE__,__LINE__,mw%comm)
            end select
            select type(a=>sh%dtss); type is(lo_voronoi_distancetable)
                if ( wzd(1) .eq. -2 ) then
                    call a%generate_wzcutoff(ss%r,ss%latticevectors,[1,1,1],ss%latticevectors)
                else
                    call a%generate_wzcutoff(ss%r,ss%latticevectors,wzd,uc%latticevectors,vorona=uc%na)
                endif
            end select
        else
            ! make perfectly normal distance tables
            allocate(lo_distancetable::sh%dtuc)
            allocate(lo_distancetable::sh%dtss)
            f0=max(rc2,rc3,rc4,mc2,dc2,dc3)
            select type(a=>sh%dtuc); type is(lo_distancetable)
                call a%generate(uc%r,uc%latticevectors,f0,verbosity=0)
            end select
            select type(a=>sh%dtss); type is(lo_distancetable)
                call a%generate(ss%r,ss%latticevectors,f0,verbosity=0)
            end select
        endif

        ! What will make life a lot easier later on is if I align the distance tables.
        call align_distance_tables(sh%dtuc,sh%dtss,ss,lo_tol,mw,mem)
        if ( verbosity .gt. 0 ) then
            t1=walltime()
            write(lo_iou,*) '... constructed distance tables (',tochar(t1-t0),'s)'
            t0=t1
        endif

        ! While I have the information handy, construct the empty tuplets
        call construct_tuplets(sl,uc,ss,sh,rc3,rc4,dc3,verbosity)
        if ( verbosity .gt. 0 ) then
            t1=walltime()
            write(lo_iou,*) '... built prototype tuplets (',tochar(t1-t0),'s)'
            t0=t1
        endif

        ! Now I can kind of know what I should consider
        sl%consider_pair=.false.
        if ( rc2 .gt. 0.0_r8 ) sl%consider_pair=.true.
        if ( mc2 .gt. 0.0_r8 ) sl%consider_pair=.true.
        if ( dc2 .gt. 0.0_r8 ) sl%consider_pair=.true.
        if ( wzcutoff ) sl%consider_pair=.true.
        sl%consider_triplet=.false.
        if ( rc3 .gt. 0.0_r8 ) sl%consider_triplet=.true.
        if ( dc3 .gt. 0.0_r8 ) sl%consider_triplet=.true.
        sl%consider_quartet=.false.
        if ( rc4 .gt. 0.0_r8 ) sl%consider_quartet=.true.

        ! Build all the relevant symmetry operations
        if ( present(transposition) ) then
            transsym=transposition
        else
            transsym=.true.
        endif
        if ( present(spacegroup) ) then
            groupsym=spacegroup
        else
            groupsym=.true.
        endif
        call expandoperations(sl,uc,groupsym,transsym,mem)
        if ( verbosity .gt. 0 ) then
            t1=walltime()
            write(lo_iou,*) '... expanded symmetry operations (',tochar(t1-t0),'s)'
            t0=t1
        endif

        ! And get the unique tuplets
        call tuplets_to_shells(sl,uc,ss,sh,mw,mem,verbosity)

        ! This is a sensible time to say what we got:
        if ( verbosity .gt. 0 ) then

            l=20
            if ( allocated(sh%singlet) ) then
                j=size(sh%singlet)
                write(lo_iou,*) 'N singlet shells:',j
                do i=1,j
                    write(lo_iou,"(1X,A15,I3,1X,I3)") 'shell:',i,sh%singlet(i)%n_unfold
                    if ( i .ge. l ) then
                        if ( j .gt. l ) write(lo_iou,*) '    ... and ',tochar(j-l),' more shells'
                        exit
                    endif
                enddo
            endif
            if ( allocated(sh%pair) ) then
                j=size(sh%pair)
                write(lo_iou,*) 'N pair shells:',j
                do i=1,j
                    write(lo_iou,*) '   shell:',i,sh%pair(i)%n_unfold,norm2(sh%pair(i)%v)*lo_bohr_to_A
                    if ( i .ge. l ) then
                        if ( j .gt. l ) write(lo_iou,*) '    ... and ',tochar(j-l),' more shells'
                        exit
                    endif
                enddo
            endif
            if ( allocated(sh%triplet) ) then
                j=size(sh%triplet)
                write(lo_iou,*) 'N triplet shells:',j
                do i=1,j
                    write(lo_iou,*) '   shell:',i,sh%triplet(i)%n_unfold
                    if ( i .ge. l ) then
                        if ( j .gt. l ) write(lo_iou,*) '    ... and ',tochar(j-l),' more shells'
                        exit
                    endif
                enddo
            endif
            if ( allocated(sh%quartet) ) then
                j=size(sh%quartet)
                write(lo_iou,*) 'N quartet shells:',j
                do i=1,j
                    write(lo_iou,*) '   shell:',i,sh%quartet(i)%n_unfold
                    if ( i .ge. l ) then
                        if ( j .gt. l ) write(lo_iou,*) '    ... and ',tochar(j-l),' more shells'
                        exit
                    endif
                enddo
            endif
        endif
    end block distancetable

    ! Now start counting components for a variety of different kind of tuplets
    nspace: block
        ! Sort out singlet things, always do this since it's really fast
        call nullspace_fc_singlet(sl,sh,uc,mw,mem)
        if ( verbosity .gt. 0 ) write(*,*) '... determined singlet nullspace'

        ! Sort out forceconstant pairs
        call nullspace_fc_pair(sl,sh,ss,mw,mem)
        if ( verbosity .gt. 0 ) write(*,*) '... determined pair nullspace'
        ! Switch off if not desired.
        if ( present(firstorder) ) then
            if ( .not.firstorder ) then
                sl%have_fc_singlet=.false.
                sl%nx_fc_singlet=0
                sl%n_fc_singlet_shell=0
            endif
        endif
        ! Sort out triplets, maybe?
        if ( rc3 .gt. 0.0_r8 ) then
            call nullspace_fc_triplet(sl,sh,ss,mw,mem)
            if ( verbosity .gt. 0 ) write(*,*) '... determined triplet nullspace'
        else
            sl%have_fc_triplet=.false.
            sl%nx_fc_triplet=0
            sl%n_fc_triplet_shell=0
        endif

        ! And perhaps quartets
        if ( rc4 .gt. 0.0_r8 ) then
            call nullspace_fc_quartet(sl,sh,ss,mw,mem)
            if ( verbosity .gt. 0 ) write(*,*) '... determined quartet nullspace'
        else
            sl%have_fc_quartet=.false.
            sl%nx_fc_quartet=0
            sl%n_fc_quartet_shell=0
        endif

        ! Do the polar things? To not get lost in the heuristics I do it simple:
        ! always get the nullspace for the dielectric constant and Born charges.
        ! if not polar, switch them off.
        call nullspace_eps_global(sl,mw,mem)
        if ( verbosity .gt. 0 ) write(*,*) '... determined dielectric tensor nullspace'
        call nullspace_Z_singlet(sl,sh,uc,mw,mem)
        if ( verbosity .gt. 0 ) write(*,*) '... determined Born charge nullspace'

        if ( polar ) then
            if ( .not.sl%have_Z_singlet ) then
                ! Switch off the polar stuff in case all Z disappeared. Maybe say that?
                if ( mw%talk ) then
                    write(lo_iou,*) ''
                    write(lo_iou,*) 'WARNING: you specified polar corrections, but due to symmetry'
                    write(lo_iou,*) 'all born charges are forced to be zero. Make sure that this is'
                    write(lo_iou,*) 'what is expected. Will continue, but switching off polar things.'
                    write(lo_iou,*) ''
                endif
                sl%have_Z_singlet=.false.
                sl%have_eps_global=.false.
                sl%nx_Z_singlet=0
                sl%nx_eps_global=0
            endif
        elseif ( dc2 .lt. 0.0_r8 ) then
            ! Switch off the polar stuff
            sl%have_Z_singlet=.false.
            sl%have_eps_global=.false.
            sl%nx_Z_singlet=0
            sl%nx_eps_global=0
        endif
        ! Cost nothing to get the first order response of eps so I do it all the time. Or not.
        ! Actually not, this behaves strange when there are 100+ atoms in the unit cell. Some
        ! memory thing I think.
        ! FK: this breaks Raman detection in any material, check for no. of atoms instead
        ! if ( max(dc2,dc3) .gt. 0.0_r8 ) then
        if ( uc%na < 33 ) then
            call nullspace_eps_singlet(sl,sh,uc,mw,mem)
        else
            if ( mw%talk ) then
                write(lo_iou,*) 'FIXME: Large unitcell with ', uc%na, ' atoms used'
                write(lo_iou,*) 'FIXME: Raman activity will NOT be reported'
            end if
            sl%have_eps_singlet=.false.
            sl%nx_eps_singlet=0
            sl%n_eps_singlet_shell=0
        endif

        ! And the second order response, a bit fancier
        if ( dc3 .gt. 0.0_r8 ) then
            call nullspace_eps_pair(sl,sh,uc,ss,dc3,mw,mem)
        else
            sl%have_eps_pair=.false.
            sl%nx_eps_pair=0
            sl%n_eps_pair_shell=0
        endif

        if ( dc2 .gt. 0.0_r8 ) then
            call nullspace_Z_pair(sl,sh,uc,ss,dc2,mw,mem)
        else
            sl%have_Z_pair=.false.
            sl%nx_Z_pair=0
            sl%n_Z_pair_shell=0
        endif

        if ( dc3 .gt. 0.0_r8 ) then
            call nullspace_Z_triplet(sl,sh,uc,ss,dc3,mw,mem)
        else
            sl%have_Z_triplet=.false.
            sl%nx_Z_triplet=0
            sl%n_Z_triplet_shell=0
        endif

        !! Magnetic exchange interactions?
        !if ( mc2 .gt. 0.0_r8 ) then
        !    call nullspace_magnetic_pair(sl,sh,uc,mc2,mw,mem)
        !else
        !    sl%have_mag_pair=.false.
        !endif
    end block nspace

    if ( verbosity .gt. 0 ) then
    infodump: block
        write(lo_iou,*) 'Interactions:'
        write(lo_iou,*) '       firstorder forceconstant:',sl%n_fc_singlet_shell,sl%nx_fc_singlet
        write(lo_iou,*) '      secondorder forceconstant:',sl%n_fc_pair_shell,sl%nx_fc_pair
        write(lo_iou,*) '       thirdorder forceconstant:',sl%n_fc_triplet_shell,sl%nx_fc_triplet
        write(lo_iou,*) '      fourthorder forceconstant:',sl%n_fc_quartet_shell,sl%nx_fc_quartet
        write(lo_iou,*) '            dielectric constant:',0,sl%nx_eps_global
        write(lo_iou,*) '             dielectric singlet:',sl%n_eps_singlet_shell,sl%nx_eps_singlet
        write(lo_iou,*) '                dielectric pair:',sl%n_eps_pair_shell,sl%nx_eps_pair
        write(lo_iou,*) '  born effective charge singlet:',sl%n_Z_singlet_shell,sl%nx_Z_singlet
        write(lo_iou,*) '     born effective charge pair:',sl%n_Z_pair_shell,sl%nx_Z_pair
        write(lo_iou,*) '  born effective charge triplet:',sl%n_Z_triplet_shell,sl%nx_Z_triplet
    end block infodump
    endif

    call mem%tock(__FILE__,__LINE__)
end subroutine

end module
