#include "precompilerdefinitions"
module type_symmetrylist
    !! Type that contains all tuplets in a cell, and how they are related to each other.
    use konstanter, only: flyt,lo_tol,lo_sqtol,lo_huge,lo_hugeint,lo_status,lo_exitcode_symmetry,lo_exitcode_param,lo_bohr_to_A
    use gottochblandat, only: lo_chop,tochar,lo_sqnorm,walltime,lo_progressbar,lo_outerproduct,lo_fetch_tolerance,&
                       lo_flattentensor,lo_frobnorm,lo_real_gram_schmidt,lo_return_tensor_transpositions,&
                       lo_return_unique,lo_identitymatrix,lo_progressbar_init,lo_progressbar,qsort,lo_stop_gracefully,&
                       lo_nullspace_coefficient_matrix,lo_determ,lo_transpositionmatrix,lo_make_coeffmatrix_tidy
    use type_crystalstructure, only: lo_crystalstructure
    use type_symmetryoperation, only: lo_symset,lo_spacegroup_operation,lo_expandoperation_pair
    use type_distancetable, only: lo_distancetable,lo_distancetable_particle
    use type_voronoi_distancetable, only: lo_voronoi_distancetable
    use type_blas_lapack_wrappers, only: lo_dgesdd,lo_gemm

    implicit none
    private
    public :: lo_symlist
    public :: lo_symlist_singletop
    public :: lo_symlist_pairop
    public :: lo_symlist_tripletop
    public :: lo_symlist_quartetop

    !> a generic tuplet
    type, abstract :: lo_symlist_tuplet
        !> which unique shell is it part of
        integer :: unique_shell=-lo_hugeint
    end type
        !> generic singlet
        type, extends(lo_symlist_tuplet) :: lo_symlist_singlet
            !> operation that takes the unique shell to this pair
            integer :: operation_from_shell_to_atom=-lo_hugeint
        end type
            !> singlet in the unitcell
            type, extends(lo_symlist_singlet) :: lo_symlist_ucsinglet
                !> indices to the supercell
                integer :: ssa=-lo_hugeint
            end type
            !> singlet in the supercell
            type, extends(lo_symlist_singlet) :: lo_symlist_sssinglet
                !> index to atoms in the supercell
                integer :: j1=-lo_hugeint
            end type

        !> generic pair
        type, extends(lo_symlist_tuplet) :: lo_symlist_pair
            !> index to atoms in the unitcell
            integer :: i1=-lo_hugeint,i2=-lo_hugeint
            !> unique classifier for the atoms
            integer :: ui1=-lo_hugeint,ui2=-lo_hugeint
            !> vector to second atom
            real(flyt), dimension(3) :: v=lo_huge
            !> operation that takes the unique shell to this pair
            integer :: operation_from_shell_to_pair=-lo_hugeint
        end type
            !> pair in the unitcell
            type, extends(lo_symlist_pair) :: lo_symlist_ucpair
                !> indices to the supercell
                integer :: ssa=-lo_hugeint,ssi=-lo_hugeint
            end type
            !> pair in the supercell
            type, extends(lo_symlist_pair) :: lo_symlist_sspair
                !> index to atoms in the supercell
                integer :: j1=-lo_hugeint,j2=-lo_hugeint
            end type

        !> generic triplet
        type, extends(lo_symlist_tuplet) :: lo_symlist_triplet
            !> index to atoms
            integer :: i1=-lo_hugeint,i2=-lo_hugeint,i3=-lo_hugeint
            !> unique classifier for the atoms
            integer :: ui1=-lo_hugeint,ui2=-lo_hugeint,ui3=-lo_hugeint
            !> vector to other atoms
            real(flyt), dimension(3) :: v1=lo_huge,v2=lo_huge,v3=lo_huge
            !> operation that takes the unique shell to this pair
            integer :: operation_from_shell_to_triplet=-lo_hugeint
        end type
            !> pair in the unitcell
            type, extends(lo_symlist_triplet) :: lo_symlist_uctriplet
                !> indices to the supercell
                integer :: ssa=-lo_hugeint,ssi=-lo_hugeint
            end type
            !> pair in the supercell
            type, extends(lo_symlist_triplet) :: lo_symlist_sstriplet
                !> indices in supercell
                integer :: j1=-lo_hugeint,j2=-lo_hugeint,j3=-lo_hugeint
            end type

        !> generic quartet
        type, extends(lo_symlist_tuplet) :: lo_symlist_quartet
            !> index to atoms
            integer :: i1=-lo_hugeint,i2=-lo_hugeint,i3=-lo_hugeint,i4=-lo_hugeint
            !> unique classifier for the atoms
            integer :: ui1=-lo_hugeint,ui2=-lo_hugeint,ui3=-lo_hugeint,ui4=-lo_hugeint
            !> vector to other atoms
            real(flyt), dimension(3) :: v1=lo_huge,v2=lo_huge,v3=lo_huge,v4=lo_huge
            !> operation that takes the unique shell to this pair
            integer :: operation_from_shell_to_quartet=-lo_hugeint
        end type
            !> pair in the unitcell
            type, extends(lo_symlist_quartet) :: lo_symlist_ucquartet
                !> indices to the supercell
                integer :: ssa=-lo_hugeint,ssi=-lo_hugeint
            end type
            !> pair in the supercell
            type, extends(lo_symlist_quartet) :: lo_symlist_ssquartet
                !> indices in supercell
                integer :: j1=-lo_hugeint,j2=-lo_hugeint,j3=-lo_hugeint,j4=-lo_hugeint
            end type

        !> generic magnetic pair
        type, extends(lo_symlist_tuplet) :: lo_symlist_magpair
            !> index to atoms in the unitcell
            integer :: i1=-lo_hugeint,i2=-lo_hugeint
            !> unique classifier for the atoms
            integer :: ui1=-lo_hugeint,ui2=-lo_hugeint
            !> vector to second atom
            real(flyt), dimension(3) :: v=lo_huge
            !> operation that takes the unique shell to this pair
            integer :: operation_from_shell_to_jij=-lo_hugeint
            integer :: operation_from_shell_to_bird=-lo_hugeint
            integer :: operation_from_shell_to_dude=-lo_hugeint
            integer :: operation_from_shell_to_qij=-lo_hugeint
        end type
            !> magnetic pair in the unitcell
            type, extends(lo_symlist_magpair) :: lo_symlist_ucmagpair
                !> indices to the supercell
                integer :: ssa=-lo_hugeint,ssi=-lo_hugeint
            end type
            !> magnetic pair in the supercell
            type, extends(lo_symlist_magpair) :: lo_symlist_ssmagpair
                !> index to atoms in the supercell
                integer :: j1=-lo_hugeint,j2=-lo_hugeint
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
    type, abstract :: lo_symlist_tupletop
        !> which transposition, in the original notation
        integer :: permind=-lo_hugeint
        !> which operation, in original notation
        integer :: opind=-lo_hugeint
    end type
    !> singlet operations
    type, extends(lo_symlist_tupletop) :: lo_symlist_singletop
        !> 3x3 matrix to rotate stuff
        real(flyt), dimension(3,3) :: m3=lo_huge,im3=lo_huge
    end type
    !> pair operations
    type, extends(lo_symlist_tupletop) :: lo_symlist_pairop
        !> Tr*Op matrix
        real(flyt), dimension(9,9) :: sotr=lo_huge
        !> Tr*Op*det(Op) matrix, for spin-lattice terms
        real(flyt), dimension(9,9) :: msotr=lo_huge
        !> 3x3 matrix to rotate stuff
        real(flyt), dimension(3,3) :: m3=lo_huge,im3=lo_huge
        !> array that permutes indices
        integer, dimension(2) :: perm=-lo_hugeint
    end type
    !> triplet operations
    type, extends(lo_symlist_tupletop) :: lo_symlist_tripletop
        !> Tr*Op matrix
        real(flyt), dimension(27,27) :: sotr=lo_huge
        !> 3x3 matrix to rotate stuff
        real(flyt), dimension(3,3) :: m3=lo_huge,im3=lo_huge
        !> array that permutes indices
        integer, dimension(3) :: perm=-lo_hugeint
    end type
    !> quartet operations
    type, extends(lo_symlist_tupletop) :: lo_symlist_quartetop
        !> Tr*Op matrix
        real(flyt), dimension(81,81) :: sotr=lo_huge
        !> 3x3 matrix to rotate stuff
        real(flyt), dimension(3,3) :: m3=lo_huge,im3=lo_huge
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
    type, abstract :: lo_symlist_shell
        !> number of things that belong to this shell
        integer :: n
    end type
    !> shell for generic forceconstants
    type, extends(lo_symlist_shell) :: lo_symlist_latticeshell
        !> number of force constants
        integer :: nfc=-lo_hugeint
        !> which are the relevant forceconstants
        integer, dimension(:), allocatable :: relfc,relfcind
        !> coefficient matrix
        real(flyt), dimension(:,:), allocatable :: coeff_M
    end type
            !> singlet coordination shell
            type, extends(lo_symlist_latticeshell) :: lo_symlist_singletshell
                !> the prototype atom, index in the unit cell
                integer :: protatom=-lo_hugeint
                !> coefficient matrix
                real(flyt), dimension(3,3) :: coeffM=lo_huge
            end type
            !> pair coordination shell
            type, extends(lo_symlist_latticeshell) :: lo_symlist_pairshell
                !> the prototype pair
                type(lo_symlist_ucpair) :: protpair
                !> coefficient matrix for this shell
                real(flyt), dimension(9,9) :: coeffM=lo_huge,invarM=0.0_flyt
                !> indices to the atoms of the pairs in this shell
                integer, dimension(:), allocatable :: atomind
                !> indices to the pairs for those atoms for all pairs in this shell
                integer, dimension(:), allocatable :: pairind
                !> how many operations leave this shell invariant
                integer :: noperations=-lo_hugeint
                !> which operations leave it invariant
                integer, dimension(:), allocatable :: operations
            end type
            !> triplet coordination shell
            type, extends(lo_symlist_latticeshell) :: lo_symlist_tripletshell
                !> the prototype triplet
                type(lo_symlist_uctriplet) :: prottriplet
                !> coefficient matrix
                real(flyt), dimension(27,27) :: coeffM=lo_huge
            end type
            !> quartet coordination shell
            type, extends(lo_symlist_latticeshell) :: lo_symlist_quartetshell
                !> the prototype quartet
                type(lo_symlist_ucquartet) :: protquartet
                !> coefficient matrix
                real(flyt), dimension(81,81) :: coeffM=lo_huge
            end type

    !> coordination shell for Born charges
    type, extends(lo_symlist_shell) :: lo_symlist_Zshell
        !> coefficent matrix
        real(flyt), dimension(:,:), allocatable :: coeff_M
    end type
    type, extends(lo_symlist_shell) :: lo_symlist_epsshell
        !> coefficient matrix
        real(flyt), dimension(:,:), allocatable :: coeff_M
    end type
    !> magnetic pair coordination shell
    type, extends(lo_symlist_shell) :: lo_symlist_magpairshell
        !> number of irreducible components
        integer :: nfc_jij=-lo_hugeint,nfc_tij=-lo_hugeint,nfc_qij=-lo_hugeint
        !> the prototype pair
        type(lo_symlist_ucmagpair) :: protpair
        !> coefficient matrices
        real(flyt), dimension(:,:), allocatable :: coeffM_jij,coeffM_tij,coeffM_qij
        !> which are the relevant forceconstants
        integer, dimension(:), allocatable :: relfc_jij,relfcind_jij
        integer, dimension(:), allocatable :: relfc_tij,relfcind_tij
        integer, dimension(:), allocatable :: relfc_qij,relfcind_qij

        !> indices to the atoms of the pairs in this shell
        integer, dimension(:), allocatable :: atomind
        !> indices to the pairs for those atoms for all pairs in this shell
        integer, dimension(:), allocatable :: pairind
    end type
    type, extends(lo_symlist_shell) :: lo_symlist_magsingletshell
        !> the prototype atom
        integer :: protatom=-lo_hugeint
    end type
    type, extends(lo_symlist_shell) :: lo_symlist_epssingletshell
        !> the prototype atom
        integer :: protatom=-lo_hugeint
        !> coefficient matrix
        real(flyt), dimension(:,:), allocatable :: coeff_M
    end type
    type, extends(lo_symlist_shell) :: lo_symlist_epspairshell
        !> prototype pair
        type(lo_symlist_ucpair) :: protpair
        !> coefficient matrix
        real(flyt), dimension(:,:), allocatable :: coeff_M
        !> number of independent
        integer :: nx=-lo_hugeint
        !> global offset in x-indices
        integer :: offset=-lo_hugeint
    end type
    type, extends(lo_symlist_shell) :: lo_symlist_Zpairshell
        !> prototype pair
        type(lo_symlist_ucpair) :: protpair
        !> coefficient matrix
        real(flyt), dimension(:,:), allocatable :: coeff_M
    end type
    type, extends(lo_symlist_shell) :: lo_symlist_Ztripletshell
        !> the prototype triplet
        type(lo_symlist_uctriplet) :: prottriplet
        !> coefficient matrix
        real(flyt), dimension(:,:), allocatable :: coeff_M
    end type


    ! For each atom in the system, I define all the tuplets it belongs to.
    type lo_symlist_ucatom
        integer :: npair=-lo_hugeint
        integer :: ntriplet=-lo_hugeint
        integer :: nquartet=-lo_hugeint
        integer :: nmagpair=-lo_hugeint
        integer :: nepspair=-lo_hugeint
        integer :: nZpair=-lo_hugeint
        integer :: nZtriplet=-lo_hugeint
        type(lo_symlist_ucsinglet)                            :: singlet
        type(lo_symlist_ucpair),    dimension(:), allocatable :: pair
        type(lo_symlist_uctriplet), dimension(:), allocatable :: triplet
        type(lo_symlist_ucquartet), dimension(:), allocatable :: quartet
        type(lo_symlist_ucmagpair), dimension(:), allocatable :: magpair
        type(lo_symlist_ucsinglet)                            :: magsinglet
        type(lo_symlist_ucsinglet)                            :: epssinglet
        type(lo_symlist_ucpair),    dimension(:), allocatable :: epspair
        type(lo_symlist_ucpair),    dimension(:), allocatable :: Zpair
        type(lo_symlist_uctriplet), dimension(:), allocatable :: Ztriplet
    end type
    type lo_symlist_ssatom
        integer :: npair=-lo_hugeint
        integer :: ntriplet=-lo_hugeint
        integer :: nquartet=-lo_hugeint
        integer :: nmagpair=-lo_hugeint
        integer :: nepspair=-lo_hugeint
        integer :: nZpair=-lo_hugeint
        integer :: nZtriplet=-lo_hugeint
        type(lo_symlist_sssinglet)                            :: singlet
        type(lo_symlist_sspair),    dimension(:), allocatable :: pair
        type(lo_symlist_sstriplet), dimension(:), allocatable :: triplet
        type(lo_symlist_ssquartet), dimension(:), allocatable :: quartet
        type(lo_symlist_ssmagpair), dimension(:), allocatable :: magpair
        type(lo_symlist_sssinglet)                            :: magsinglet
        type(lo_symlist_sssinglet)                            :: epssinglet
        type(lo_symlist_sspair),    dimension(:), allocatable :: epspair
        type(lo_symlist_sspair),    dimension(:), allocatable :: Zpair
        type(lo_symlist_sstriplet), dimension(:), allocatable :: Ztriplet
    end type

    !> The complete symmetry information. Details which pairs are equal to which pairs, and so on.
    type lo_symlist
        !> How much to talk
        integer :: verbosity=-lo_hugeint
        !> cutoff radii
        real(flyt) :: rc2=-lo_huge,rc3=-lo_huge,rc4=-lo_huge
        !> cutoff in terms of jumps
        integer :: nj2=-lo_hugeint,nj3=-lo_hugeint,nj4=-lo_hugeint
        !> how many neighbours are defined to be in the first coordination shell?
        !integer, dimension(:), allocatable :: coordination
        !> which orders to consider
        logical :: firstorder=.false.,secondorder=.false.,thirdorder=.false.,fourthorder=.false.
        !> dimensions of WZ cell
        logical :: wzcutoff=.false.
        integer, dimension(3) :: wzdim=-lo_hugeint
        !> when testing for symmetry, consider other things as well
        logical :: testcoll=.false.,testnoncoll=.false.,transposition=.false.,spacegroup=.false.
        !> tolerances
        real(flyt) :: sqtol_cart=lo_huge,sqtol_frac=lo_huge,sqtol_mag=lo_huge
        real(flyt) :: tol_cart=lo_huge,tol_frac=lo_huge,tol_mag=lo_huge
        !> determine nullspace
        logical :: nspace=.false.
        !> consider a polar material
        logical :: polar=.false.
        !> consider interactions in a magnetic material
        logical :: magnetic_pair_interactions=.false.,magnetic_singlet_interactions=.false.
        !> consider dielectric expansion
        logical :: dielectric_pair=.false.,dielectric_triplet=.false.
        !> cutoff radii for magnetic interactions
        real(flyt) :: mc2=-lo_huge
        !> cutoff radii for dielectric interactions
        real(flyt) :: dc2=-lo_hugeint,dc3=-lo_hugeint
        !> list of symmetry operations
        type(lo_symset) :: sym
        integer :: npairop=-lo_hugeint,ntripletop=-lo_hugeint,nquartetop=-lo_hugeint
        type(lo_symlist_pairop), dimension(:), allocatable :: pairop
        type(lo_symlist_tripletop), dimension(:), allocatable :: tripletop
        type(lo_symlist_quartetop), dimension(:), allocatable :: quartetop

        !> information about all the atoms
        integer :: nuc=-lo_hugeint,nss=-lo_hugeint
        type(lo_symlist_ucatom), dimension(:), allocatable :: uc
        type(lo_symlist_ssatom), dimension(:), allocatable :: ss

        !> unique coordination shells of each kind
        integer :: nsingletshells=-lo_hugeint
        integer :: npairshells=-lo_hugeint
        integer :: ntripletshells=-lo_hugeint
        integer :: nquartetshells=-lo_hugeint
        integer :: nmagpairshells=-lo_hugeint
        integer :: nmagsingletshells=-lo_hugeint
        integer :: nzshellpair=-lo_hugeint
        integer :: nzshelltriplet=-lo_hugeint
        integer :: nepsshellsinglet=-lo_hugeint
        integer :: nepsshellpair=-lo_hugeint

        type(lo_symlist_singletshell),    dimension(:), allocatable :: singletshell
        type(lo_symlist_pairshell),       dimension(:), allocatable :: pairshell
        type(lo_symlist_tripletshell),    dimension(:), allocatable :: tripletshell
        type(lo_symlist_quartetshell),    dimension(:), allocatable :: quartetshell
        type(lo_symlist_Zshell)                                     :: Zshell
        type(lo_symlist_epsshell)                                   :: epsshell
        type(lo_symlist_magpairshell),    dimension(:), allocatable :: magpairshell
        type(lo_symlist_magsingletshell), dimension(:), allocatable :: magsingletshell
        type(lo_symlist_epssingletshell)                            :: epssingletshell
        type(lo_symlist_epspairshell),    dimension(:), allocatable :: epspairshell
        type(lo_symlist_Zpairshell),      dimension(:), allocatable :: Zpairshell
        type(lo_symlist_Ztripletshell),   dimension(:), allocatable :: Ztripletshell
        !type(lo_symlist_epssingletshell),    dimension(:), allocatable :: epssingletshell
        integer :: nsingletifc=-lo_hugeint
        integer :: npairifc=-lo_hugeint
        integer :: ntripletifc=-lo_hugeint
        integer :: nquartetifc=-lo_hugeint
        integer :: nepstheta=-lo_hugeint
        integer :: nepstheta_singlet=-lo_hugeint
        integer :: nepstheta_pair=-lo_hugeint
        integer :: nZtheta=-lo_hugeint
        integer :: nZtheta_pair=-lo_hugeint
        integer :: nZtheta_triplet=-lo_hugeint
        integer :: nmagsinglettheta=-lo_hugeint
        integer :: nmagpairtheta_jij=-lo_hugeint
        integer :: nmagpairtheta_tij=-lo_hugeint
        integer :: nmagpairtheta_qij=-lo_hugeint
        contains
            !> create the whole thing
            procedure :: generate
#ifdef AGRESSIVE_SANITY
            procedure :: size_in_mem
#endif
    end type

    !> Helper type with intermediate quantities, only exists temporarily
    type lo_symtabhelper
        ! Distance table for the unit cell
        class(lo_distancetable), allocatable :: dtuc
        ! Distance table for the supercell
        class(lo_distancetable), allocatable :: dtss
        ! operation connection list between unitcell atoms
        logical, dimension(:,:,:), allocatable :: clist

        logical, dimension(:), allocatable :: is_unitcell_atom_prototype
        integer, dimension(:), allocatable :: unique_unitcell_atom
        integer, dimension(:), allocatable :: unique_unitcell_operation
        integer, dimension(:), allocatable :: unique_supercell_atom
        integer, dimension(:), allocatable :: unique_supercell_operation

        ! Lists of all the unique tuplets, intermediate values
        type(lo_symlist_ucpair), dimension(:), allocatable :: allpairs
        type(lo_symlist_uctriplet), dimension(:), allocatable :: alltriplets
        type(lo_symlist_ucquartet), dimension(:), allocatable :: allquartets
        type(lo_symlist_ucmagpair), dimension(:), allocatable :: allmagpairs
        type(lo_symlist_ucpair), dimension(:), allocatable :: allepspairs
        type(lo_symlist_ucpair), dimension(:), allocatable :: allZpairs
        type(lo_symlist_uctriplet), dimension(:), allocatable :: allZtriplets
        ! Hashes for all the tuplets, for fast comparison
        real(flyt), dimension(:), allocatable :: pair_allhash,pair_uchash,pair_sshash
        real(flyt), dimension(:), allocatable :: triplet_allhash,triplet_uchash,triplet_sshash
        real(flyt), dimension(:), allocatable :: quartet_allhash,quartet_uchash,quartet_sshash
        real(flyt), dimension(:), allocatable :: magpair_allhash,magpair_uchash,magpair_sshash
        real(flyt), dimension(:), allocatable :: epspair_allhash,epspair_uchash,epspair_sshash
        real(flyt), dimension(:), allocatable :: Zpair_allhash,Zpair_uchash,Zpair_sshash
        real(flyt), dimension(:), allocatable :: Ztriplet_allhash,Ztriplet_uchash,Ztriplet_sshash
    end type

    ! interfaces to the helper submodule
    interface
        module subroutine expandoperations(sl)
            type(lo_symlist), intent(inout) :: sl
        end subroutine
        module subroutine check_cutoffs(uc,ss,rc2,rc3,rc4,wraparound)
            type(lo_crystalstructure), intent(in) :: uc,ss
            real(flyt), intent(inout) :: rc2,rc3,rc4
            logical, intent(in) :: wraparound
        end subroutine
        module function hashtuplet(t) result(hash)
            class(lo_symlist_tuplet), intent(in) :: t
            real(flyt) :: hash
        end function
        module pure function expandoperation_triplet(o) result(bigo)
            real(flyt), dimension(3,3), intent(in) :: o
            real(flyt), dimension(27,27) :: bigo
        end function
    end interface
    ! interfaces to the unique submodule
    interface
        module subroutine findunique(sl,sh,uc,ss)
            type(lo_symlist), intent(inout) :: sl
            type(lo_symtabhelper), intent(inout) :: sh
            type(lo_crystalstructure), intent(in) :: uc
            type(lo_crystalstructure), intent(in) :: ss
        end subroutine
        module subroutine setuptuplets(sl,sh,uc,ss)
            type(lo_symlist), intent(inout) :: sl
            type(lo_symtabhelper), intent(inout) :: sh
            type(lo_crystalstructure), intent(in) :: uc,ss
        end subroutine
        module subroutine get_unique_tuplets(sl,sh,uc,ss)
            type(lo_symlist), intent(inout) :: sl
            type(lo_symtabhelper), intent(inout) :: sh
            type(lo_crystalstructure), intent(in) :: uc
            type(lo_crystalstructure), intent(in) :: ss
        end subroutine
        module subroutine validops_pairshells(sl,sh,ss)
            type(lo_symlist), intent(inout) :: sl
            type(lo_symtabhelper), intent(inout) :: sh
            type(lo_crystalstructure), intent(in) :: ss
        end subroutine
        module function compare_tuplets_after_one_is_rotated(sl,sh,ss,t1,t2,op,testperm) result(match)
            type(lo_symlist), intent(in) :: sl
            type(lo_symtabhelper), intent(in) :: sh
            type(lo_crystalstructure), intent(in) :: ss
            class(lo_symlist_tuplet), intent(in) :: t1,t2
            class(lo_symlist_tupletop), intent(in) :: op
            logical, intent(in), optional :: testperm
            logical :: match
        end function
    end interface
    ! interfaces to the coefficientmatrix submodule
    interface
        module subroutine nullspace_singletshells(sl,sh)
            type(lo_symlist), intent(inout) :: sl
            type(lo_symtabhelper), intent(inout) :: sh
        end subroutine
        module subroutine nullspace_pairshells(sl,sh,ss)
            type(lo_symlist), intent(inout) :: sl
            type(lo_symtabhelper), intent(inout) :: sh
            type(lo_crystalstructure), intent(in) :: ss
        end subroutine
        module subroutine nullspace_tripletshells(sl,sh,ss)
            type(lo_symlist), intent(inout) :: sl
            type(lo_symtabhelper), intent(inout) :: sh
            type(lo_crystalstructure), intent(in) :: ss
        end subroutine
        module subroutine nullspace_quartetshells(sl,sh,ss)
            type(lo_symlist), intent(inout) :: sl
            type(lo_symtabhelper), intent(inout) :: sh
            type(lo_crystalstructure), intent(in) :: ss
        end subroutine
        module subroutine nullspace_magpairshells(sl,sh,ss)
            type(lo_symlist), intent(inout) :: sl
            type(lo_symtabhelper), intent(inout) :: sh
            type(lo_crystalstructure), intent(in) :: ss
        end subroutine
        module subroutine nullspace_Zshells(sl,uc)
            type(lo_symlist), intent(inout) :: sl
            type(lo_crystalstructure), intent(in) :: uc
        end subroutine
        module subroutine nullspace_epsshell(sl,uc)
            type(lo_symlist), intent(inout) :: sl
            type(lo_crystalstructure), intent(in) :: uc
        end subroutine
        module subroutine nullspace_epssinglet(sl,uc)
            type(lo_symlist), intent(inout) :: sl
            type(lo_crystalstructure), intent(in) :: uc
        end subroutine
    end interface

contains

!> create the symmetry list
subroutine generate(sl,uc,ss,cutoff2,cutoff3,cutoff4,polar,transposition,spacegroup,&
                    verbosity,firstorder,wzdim,nj2,nj3,nj4,tol,&
                    wraparound,magcutoff2,magsinglet,&
                    dielcutoff2,dielcutoff3)
    !> the full symmetry table thing
    class(lo_symlist), intent(out) :: sl
    !> unitcell crystal structure
    type(lo_crystalstructure), intent(inout) :: uc
    !> supercell crystal structure
    type(lo_crystalstructure), intent(inout) :: ss
    !> cutoff for second order
    real(flyt), intent(in) :: cutoff2
    !> cutoff for third order
    real(flyt), intent(in) :: cutoff3
    !> cutoff for fourth order
    real(flyt), intent(in) :: cutoff4
    !> is it a polar material
    logical, intent(in) :: polar
    !> generate transposition symmetries
    logical, intent(in), optional :: transposition
    !> generate spacegroup symmetries
    logical, intent(in), optional :: spacegroup
    !> get firstorder forceconstants
    logical, intent(in), optional :: firstorder
    !> how much to talk
    integer, intent(in), optional :: verbosity
    !> the WZ supercell dimensions
    integer, dimension(3), optional, intent(in) :: wzdim
    !> graph cutoff
    integer, optional, intent(in) :: nj2,nj3,nj4
    !> tolerance
    real(flyt), intent(in), optional :: tol
    !> allow for really long interaction thingies
    logical, intent(in), optional :: wraparound
    !> pair magnetic interactions
    real(flyt), intent(in), optional :: magcutoff2
    !> magnetic on-site term
    logical, intent(in), optional :: magsinglet
    !> dielectric expansion cutoff for pairs
    real(flyt), intent(in), optional :: dielcutoff2
    !> dielectric expansion cutoff for triplets
    real(flyt), intent(in), optional :: dielcutoff3

    type(lo_symtabhelper) :: sh
    real(flyt) :: t0
    real(flyt) :: timer_init,timer_sanity1,timer_unique

    ! Fix all the general settings first
    setoptions: block
        type(lo_distancetable) :: dt
        real(flyt), dimension(:), allocatable :: dum1,dum2
        real(flyt) :: f0,f1
        integer, dimension(:), allocatable :: nbrctr
        integer :: i,j,l,itr

        if ( present(verbosity) ) then
            sl%verbosity=verbosity
        else
            sl%verbosity=0
        endif
        if ( sl%verbosity .gt. 0 ) then
            t0=walltime()
            timer_init=walltime()
            write(*,*) ''
            write(*,*) 'Building the tuplet symmetry table'
        else
            t0=0.0_flyt
            timer_init=0.0_flyt
        endif
        ! Set the cutoffs and check them for numerical stability
        sl%rc2=cutoff2
        sl%rc3=cutoff3
        sl%rc4=cutoff4
        ! Can not do polar things if there is only one atom per cell
        if ( uc%na .gt. 1 ) then
            sl%polar=polar
        else
            sl%polar=.false.
        endif

        ! Perhaps I specify number of neighbours instead of a cutoff.
        if ( present(nj2) .or. present(nj3) .or. present(nj4) ) then
            ! get a distance table that is as long as possible
            call dt%generate(uc%r,uc%latticevectors,cutoff=ss%maxcutoff(),verbosity=0)
            lo_allocate(nbrctr(uc%na))
            ! increase the cutoff until I have the desired number of neighbours within it
            if ( present(nj2) ) then; if ( nj2 .gt. 0 ) then
                f0=ss%mincutoff()   ! shortest possible cutoff
                f1=f0*0.01_flyt     ! step to increase the cutoff with
                do itr=1,1000000
                    f0=f0+f1
                    call dt%count_neighbours_within_cutoff(f0,nbrctr)
                    l=minval(nbrctr) ! smallest number of neighbours
                    if ( l .ge. nj2 ) exit
                    if ( f0 .gt. ss%maxcutoff() ) exit
                enddo
                sl%rc2=min(f0,ss%maxcutoff())
            endif; endif
            if ( present(nj3) ) then; if ( nj3 .gt. 0 ) then;
                f0=ss%mincutoff()   ! shortest possible cutoff
                f1=f0*0.01_flyt     ! step to increase the cutoff with
                do itr=1,1000000
                    f0=f0+f1
                    call dt%count_neighbours_within_cutoff(f0,nbrctr)
                    l=minval(nbrctr) ! smallest number of neighbours
                    if ( l .ge. nj3 ) exit
                    if ( f0 .gt. ss%maxcutoff() ) exit
                enddo
                sl%rc3=min(f0,ss%maxcutoff())
            endif; endif
            if ( present(nj4) ) then; if ( nj4 .gt. 0 ) then
                f0=ss%mincutoff()   ! shortest possible cutoff
                f1=f0*0.01_flyt     ! step to increase the cutoff with
                do itr=1,1000000
                    f0=f0+f1
                    call dt%count_neighbours_within_cutoff(f0,nbrctr)
                    l=minval(nbrctr) ! smallest number of neighbours
                    if ( l .ge. nj4 ) exit
                    if ( f0 .gt. ss%maxcutoff() ) exit
                enddo
                sl%rc4=min(f0,ss%maxcutoff())
            endif; endif
            lo_deallocate(nbrctr)
        endif

        ! I might allow for wraparaound thingies? Anyway, now that the cutoffs are
        ! settled store what type of tuplets are relevant.
        if ( present(wraparound) ) then
            call check_cutoffs(uc,ss,sl%rc2,sl%rc3,sl%rc4,wraparound=wraparound)
        else
            call check_cutoffs(uc,ss,sl%rc2,sl%rc3,sl%rc4,wraparound=.false.)
        endif
        if ( present(firstorder) ) then
            sl%firstorder=firstorder
        else
            sl%firstorder=.false.
        endif
        if ( sl%rc2 .gt. 0.0_flyt ) sl%secondorder=.true.
        if ( sl%rc3 .gt. 0.0_flyt ) sl%thirdorder=.true.
        if ( sl%rc4 .gt. 0.0_flyt ) sl%fourthorder=.true.
        if ( present(wzdim) ) then
            if ( wzdim(1) .gt. 0 ) then
                sl%wzdim=wzdim
                sl%wzcutoff=.true.
            elseif ( wzdim(1) .eq. -2 ) then
                sl%wzcutoff=.true.
                sl%wzdim=-2
            else
                sl%wzdim=-1
                sl%wzcutoff=.false.
            endif
        else
            sl%wzdim=-1
            sl%wzcutoff=.false.
        endif

        ! Magnetic thingies
        if ( present(magcutoff2) ) then
            if ( magcutoff2 .gt. 0.0_flyt ) then
                sl%magnetic_pair_interactions=.true.
                sl%mc2=magcutoff2
                ! Check the cutoff so that it makes sense
                f0=-1.0_flyt
                f1=-1.0_flyt
                call check_cutoffs(uc,ss,sl%mc2,f0,f1,wraparound=.false.)
            else
                sl%mc2=-1.0_flyt
                sl%magnetic_pair_interactions=.false.
            endif
        else
            sl%mc2=-1.0_flyt
            sl%magnetic_pair_interactions=.false.
        endif
        if ( present(magsinglet) ) then
            sl%magnetic_singlet_interactions=magsinglet
        else
            sl%magnetic_singlet_interactions=.false.
        endif

        ! Dielectric things
        if ( present(dielcutoff2) ) then
            if ( dielcutoff2 .gt. 0.0_flyt ) then
                sl%dielectric_pair=.true.
                sl%dc2=dielcutoff2
            else
                sl%dielectric_pair=.false.
                sl%dc2=-1.0_flyt
            endif
        else
            sl%dielectric_pair=.false.
            sl%dc2=-1.0_flyt
        endif
        if ( present(dielcutoff3) ) then
            if ( dielcutoff2 .gt. 0.0_flyt ) then
                sl%dielectric_triplet=.true.
                sl%dc3=dielcutoff3
            else
                sl%dielectric_triplet=.false.
                sl%dc3=-1.0_flyt
            endif
        else
            sl%dielectric_triplet=.false.
            sl%dc3=-1.0_flyt
        endif
        f0=-1.0_flyt
        if ( sl%dielectric_pair ) then
            call check_cutoffs(uc,ss,sl%dc2,sl%dc3,f0,wraparound=.false.)
        endif

        ! Ok, decided on order, cutoffs and stuff. Now decide what to test
        ! for, including magnetic symmetry stuff
        if ( uc%info%collmag ) then
            sl%testcoll=.true.
            sl%testnoncoll=.false.
        elseif ( uc%info%noncollmag ) then
            sl%testcoll=.false.
            sl%testnoncoll=.true.
        else
            sl%testcoll=.false.
            sl%testnoncoll=.false.
        endif
        if ( present(transposition) ) then
            sl%transposition=transposition
        else
            sl%transposition=.true.
        endif
        if ( present(spacegroup) ) then
            sl%spacegroup=spacegroup
        else
            sl%spacegroup=.true.
        endif
        ! Set the tolerances
        if ( present(tol) ) then
            f0=tol
        else
            f0=lo_tol
        endif
        call lo_fetch_tolerance(f0,uc%latticevectors,realspace_cart_tol=sl%tol_cart,&
                                realspace_fractional_tol=sl%tol_frac)
        sl%tol_mag=f0
        sl%sqtol_cart=sl%tol_cart**2
        sl%sqtol_frac=sl%tol_frac**2
        sl%sqtol_mag=sl%tol_mag**2

        ! Make a little space for arrays, and initialize lots of stuff to zero
        sl%nuc=uc%na
        sl%nss=ss%na
        allocate(sl%uc(sl%nuc))
        allocate(sl%ss(sl%nss))
        do i=1,sl%nuc
            sl%uc(i)%npair=0
            sl%uc(i)%ntriplet=0
            sl%uc(i)%nquartet=0
        enddo
        do i=1,sl%nss
            sl%ss(i)%npair=0
            sl%ss(i)%ntriplet=0
            sl%ss(i)%nquartet=0
        enddo
        sl%npairop=0
        sl%ntripletop=0
        sl%nquartetop=0

        sl%nsingletshells=0
        sl%npairshells=0
        sl%ntripletshells=0
        sl%nquartetshells=0
        sl%nmagpairshells=0
        sl%nmagsingletshells=0
        sl%nzshellpair=0
        sl%nzshelltriplet=0
        sl%nepsshellsinglet=0
        sl%nepsshellpair=0

        sl%nsingletifc=0
        sl%npairifc=0
        sl%ntripletifc=0
        sl%nquartetifc=0
        sl%nepstheta=0
        sl%nepstheta_singlet=0
        sl%nepstheta_pair=0
        sl%nZtheta=0
        sl%nZtheta_pair=0
        sl%nZtheta_triplet=0
        sl%nmagsinglettheta=0
        sl%nmagpairtheta_jij=0
        sl%nmagpairtheta_tij=0
        sl%nmagpairtheta_qij=0

        ! Print what we did
        if ( sl%verbosity .gt. 0 ) then
            timer_init=walltime()-timer_init
            write(*,*) '... decided on the settings:'
            write(*,*) '                 consider first order: ',sl%firstorder
            write(*,*) '                  second order cutoff: ',tochar(sl%rc2*lo_bohr_to_A),'A'
            write(*,*) '                   third order cutoff: ',tochar(sl%rc3*lo_bohr_to_A),'A'
            write(*,*) '                  fourth order cutoff: ',tochar(sl%rc4*lo_bohr_to_A),'A'
            write(*,*) '                       Polar material: ',sl%polar
            write(*,*) '              use Wigner-Seitz cutoff: ',sl%wzcutoff
            write(*,*) '              Wigner-Seitz dimensions: ',tochar(sl%wzdim)
            write(*,*) '        check non-collinear magnetism: ',sl%testnoncoll
            write(*,*) '                 check transpositions: ',sl%transposition
            write(*,*) '                     check spacegroup: ',sl%spacegroup
            write(*,*) '      calculate magnetic interactions: ',sl%magnetic_pair_interactions
            write(*,*) '     cutoff for magnetic interactions: ',tochar(sl%mc2*lo_bohr_to_A),'A'
            write(*,*) '    calculate dielectric interactions: ',sl%dielectric_pair
            write(*,*) '   cutoff for dielectric interactions: ',tochar(sl%dc2*lo_bohr_to_A),'A'
            write(*,*) '   cutoff for dielectric interactions: ',tochar(sl%dc3*lo_bohr_to_A),'A'
            call dt%generate(uc%r,uc%latticevectors,cutoff=maxval([sl%rc2,sl%rc3,sl%rc4,sl%mc2,sl%dc2,sl%dc3]),verbosity=0)
            ! fetch a list of all distances
            l=0
            do i=1,dt%np
            do j=2,dt%particle(i)%n
                l=l+1
            enddo
            enddo
            lo_allocate(dum1(l))
            dum1=0.0_flyt
            l=0
            do i=1,dt%np
            do j=2,dt%particle(i)%n
                l=l+1
                dum1(l)=dt%particle(i)%d(j)
            enddo
            enddo
            ! get the unique, within a rather rough tolerance
            f0=ss%mincutoff()*0.01_flyt
            call lo_return_unique(dum1,dum2,sl%tol_cart)
            call qsort(dum2)
            dum2=dum2+sl%tol_cart
            lo_allocate(nbrctr(uc%na))
            write(*,*) '    List of the number of neighbours included, per atom, for a specific cutoff:'
            do i=1,size(dum2)
                call dt%count_neighbours_within_cutoff(dum2(i),nbrctr)
                if ( uc%na .lt. 20 ) then
                    write(*,*) '      cutoff:',tochar(dum2(i)*lo_bohr_to_A,frmt='(F10.5)'),' min # ',tochar(minval(nbrctr),-4),'  per atom:',tochar(nbrctr)
                else
                    write(*,*) '      cutoff:',tochar(dum2(i)*lo_bohr_to_A,frmt='(F10.5)'),' min # ',tochar(minval(nbrctr),-4),'  per atom:',tochar(nbrctr(1:20))
                endif
            enddo
        endif
    end block setoptions

    ! Do some sanity checks right away, easier to find discrepancies now than later
    sanitycheck1: block
        integer :: i
        real(flyt), dimension(3) :: v1
        real(flyt), dimension(3,3) :: m0,m1
        real(flyt) :: f0,f1

        if ( sl%verbosity .gt. 1 ) then
            timer_sanity1=walltime()
            write(*,*) '... checking wether the unit and supercell really match'
        endif

        ! First, match the unit and supercell, if possible. This should
        ! also halt stupid input choices.
        uc%info%verbosity=0
        call ss%classify('supercell',uc)
        m0=matmul(uc%latticevectors,ss%info%supercellmatrix)
        m1=m0-ss%latticevectors
        f0=lo_frobnorm(m1)
        if ( f0 .gt. sl%tol_cart/100 ) then
            write(*,*) 'ERROR: too large mismatch between unit and supercell'
            write(*,*) '  found these supercell lattice vectors:'
            do i=1,3
                write(*,*) ss%latticevectors(:,i)
            enddo
            write(*,*) '   when I expected these (fourth column is the mismatch):'
            do i=1,3
                write(*,*) m0(:,i),norm2(m1(:,i))
            enddo
            write(*,*) '   there is no valid reason to ever have a mismatch here, so you'
            write(*,*) '   probably did something wrong.'
            stop
        endif
        if ( sl%verbosity .gt. 1 ) then
            write(*,*) '... found this supercell matrix:'
            do i=1,3
                write(*,"(8X,3(1X,I4),8X,A,F12.8)") ss%info%supercellmatrix(:,i),'mismatch:',norm2(m1(:,i))
            enddo
        endif

        ! However, it is a little too sloppy when it comes to matching.
        ! Perhaps it's a decent idea to check for sure: I will rebuild the supercell
        ! from the unitcell, and check the mismatch
        f0=0.0_flyt
        f1=0.0_flyt
        do i=1,ss%na
            v1=(ss%info%cellindex(:,i)-1)*1.0_flyt+uc%r(:,ss%info%index_in_unitcell(i))
            v1=uc%fractional_to_cartesian(v1)
            v1=ss%cartesian_to_fractional(v1)-ss%r(:,i)
            v1=ss%displacement_fractional_to_cartesian(v1)
            f0=max(f0,norm2(v1))  ! store largest mismatch
            f1=f1+norm2(v1)/ss%na ! store average mismatch
        enddo

        if ( f0 .gt. sl%tol_cart/100 ) then
            write(*,*) 'ERROR: too large mismatch between unit and supercell'
            write(*,*) 'the largest error was ',f0,'A, this is too large.'
            write(*,*) 'with this large mismatch, I can not guarantee that the symmetry'
            write(*,*) 'analysis is correct. Also, there is no valid reason to have any mismatch'
            write(*,*) 'between the cells whatsoever, so you probably did something wrong.'
            stop
        endif

        if ( sl%verbosity .gt. 1 ) then
            timer_sanity1=walltime()-timer_sanity1
            write(*,*) '       largest mismatch unit/supercell: ',f0*lo_bohr_to_A,'A'
            write(*,*) '       average mismatch unit/supercell: ',f1*lo_bohr_to_A,'A'
            write(*,*) '... unit and supercell seems to match, good! (mellantid ',tochar(timer_sanity1),'s)'
        endif
    end block sanitycheck1

    ! Now do the actual symmetry reduction thing
    finduniquetuplets: block
        real(flyt) :: tt0

        if ( sl%verbosity .gt. 1 ) then
            timer_unique=walltime()
            write(*,*) '... finding irreducible representation'
        endif
        ! start with a list of symmetry operations.
        call sl%sym%generate(uc%latticevectors,.false.,uc%r,uc%species,possible=.true.)
        ! if so desired, switch off symmetry
        if ( sl%spacegroup .eqv. .false. ) sl%sym%n=1
        if ( sl%verbosity .gt. 1 ) write(*,*) '... got symmetryoperations'
        ! identify the unique atoms
        tt0=walltime()
        call findunique(sl,sh,uc,ss)
        if ( sl%verbosity .gt. 1 ) write(*,*) '... identified unique atoms (',tochar(walltime()-tt0),'s)'
        ! create the massive operations that are symop+transposition
        tt0=walltime()
        call expandoperations(sl)
        if ( sl%verbosity .gt. 1 ) write(*,*) '... generated compound operations (',tochar(walltime()-tt0),'s)'
        ! build all the empty tuplets
        call setuptuplets(sl,sh,uc,ss)
        tt0=walltime()
        if ( sl%verbosity .gt. 1 ) write(*,*) '... generated empty tuplets (',tochar(walltime()-tt0),'s)'
        ! find the unique tuplets
        call get_unique_tuplets(sl,sh,uc,ss)

        if ( sl%firstorder  )                call nullspace_singletshells(sl,sh)
        if ( sl%secondorder )                call nullspace_pairshells(sl,sh,ss)
        if ( sl%thirdorder  )                call nullspace_tripletshells(sl,sh,ss)
        if ( sl%fourthorder )                call nullspace_quartetshells(sl,sh,ss)
        if ( sl%polar )                      call nullspace_Zshells(sl,uc)
        if ( sl%polar )                      call nullspace_epsshell(sl,uc)
        if ( sl%magnetic_pair_interactions ) call nullspace_magpairshells(sl,sh,ss)
        if ( sl%dielectric_pair )            call nullspace_epssinglet(sl,uc)
        ! does not cost anything to store the invariant operations, at least for the pairs, so I do that by default
        ! instead of adding another input parameter, I think.
        if ( sl%secondorder ) call validops_pairshells(sl,sh,ss)

        if ( sl%verbosity .gt. 0 ) then
            timer_unique=walltime()-timer_unique
            write(*,*) '... found a unique representation'
            if ( sl%firstorder )  write(*,*) '             number of singlet shells: ',tochar(sl%nsingletshells),' (N independent=',tochar(sl%nsingletifc),')'
            if ( sl%secondorder ) write(*,*) '                number of pair shells: ',tochar(sl%npairshells),' (N independent=',tochar(sl%npairifc),')'
            if ( sl%thirdorder )  write(*,*) '             number of triplet shells: ',tochar(sl%ntripletshells),' (N independent=',tochar(sl%ntripletifc),')'
            if ( sl%fourthorder ) write(*,*) '             number of quartet shells: ',tochar(sl%nquartetshells),' (N independent=',tochar(sl%nquartetifc),')'
            if ( sl%polar )       write(*,*) '        number of unique Born charges: ',tochar(sl%nZtheta)
            if ( sl%polar )       write(*,*) '      number of components of epsilon: ',tochar(sl%nepstheta)
            if ( sl%magnetic_pair_interactions ) write(*,*) '       number of magnetic pair shells: ',tochar(sl%nmagpairshells),' (N independent=',tochar(sl%nmagpairtheta_jij),'+',tochar(sl%nmagpairtheta_tij),')'
            if ( sl%dielectric_pair ) write(*,*) '     number of components of epsilon": ',tochar(sl%nepstheta_singlet)
        endif
    end block finduniquetuplets

    if ( sl%verbosity .gt. 1 ) then
        write(*,*) 'Done building symmetry table in ',tochar(walltime()-t0),'s'
    endif
end subroutine

#ifdef AGRESSIVE_SANITY
    !> measure approximate size in memory, in bytes
    pure function size_in_mem(sl) result(mem)
        !> symmetry table
        class(lo_symlist), intent(in) :: sl
        !> size in bytes
        integer :: mem

        integer :: i,j,l
        mem=0
        l=0
        l=l+storage_size(sl)/8
        l=l+sl%sym%size_in_mem()
        mem=mem+l
        ! operations
        l=0
        if ( allocated(sl%pairop) ) l=l+sl%npairop*storage_size(sl%pairop)/8
        if ( allocated(sl%tripletop) ) l=l+sl%ntripletop*storage_size(sl%tripletop)/8
        if ( allocated(sl%quartetop) ) l=l+sl%nquartetop*storage_size(sl%quartetop)/8
        ! shells
        mem=mem+l
        l=0
        if ( allocated(sl%singletshell) ) then
            l=l+sl%nsingletshells*storage_size(sl%singletshell)/8
            do i=1,sl%nsingletshells
                j=sl%singletshell(i)%nfc
                if ( j .gt. 0 ) l=l+j*storage_size(sl%singletshell(i)%relfc)/8
                if ( j .gt. 0 ) l=l+j*storage_size(sl%singletshell(i)%relfcind)/8
                if ( allocated(sl%Zshell%coeff_M) ) l=l+storage_size(sl%Zshell%coeff_M)*size(sl%Zshell%coeff_M)
            enddo
        endif
        mem=mem+l
        l=0
        if ( allocated(sl%pairshell) ) then
            l=l+sl%npairshells*storage_size(sl%pairshell)/8
            do i=1,sl%npairshells
                j=sl%pairshell(i)%nfc
                if ( j .gt. 0 ) l=l+j*storage_size(sl%pairshell(i)%relfc)/8
                if ( j .gt. 0 ) l=l+j*storage_size(sl%pairshell(i)%relfcind)/8
                l=l+storage_size(sl%pairshell(i)%protpair)/8
            enddo
        endif
        mem=mem+l
        l=0
        if ( allocated(sl%tripletshell) ) then
            l=l+sl%ntripletshells*storage_size(sl%tripletshell)/8
            do i=1,sl%ntripletshells
                j=sl%tripletshell(i)%nfc
                if ( j .gt. 0 ) l=l+j*storage_size(sl%tripletshell(i)%relfc)/8
                if ( j .gt. 0 ) l=l+j*storage_size(sl%tripletshell(i)%relfcind)/8
                l=l+storage_size(sl%tripletshell(i)%prottriplet)/8
            enddo
        endif
        mem=mem+l
        l=0
        if ( allocated(sl%quartetshell) ) then
            l=l+sl%nquartetshells*storage_size(sl%quartetshell)/8
            do i=1,sl%nquartetshells
                j=sl%quartetshell(i)%nfc
                if ( j .gt. 0 ) l=l+j*storage_size(sl%quartetshell(i)%relfc)/8
                if ( j .gt. 0 ) l=l+j*storage_size(sl%quartetshell(i)%relfcind)/8
                l=l+storage_size(sl%quartetshell(i)%protquartet)/8
            enddo
        endif
        mem=mem+l

        l=0
        l=l+storage_size(sl%Zshell)
        if ( allocated(sl%Zshell%coeff_M) ) l=l+storage_size(sl%Zshell%coeff_M)*size(sl%Zshell%coeff_M)
        mem=mem+l/8

        l=0
        l=l+storage_size(sl%epsshell)
        if ( allocated(sl%epsshell%coeff_M) ) l=l+storage_size(sl%epsshell%coeff_M)*size(sl%epsshell%coeff_M)
        mem=mem+l/8

        l=0
        if ( allocated(sl%magpairshell) ) then
            l=l+sl%nmagpairshells*storage_size(sl%magpairshell)/8
            do i=1,sl%nmagpairshells
                j=sl%magpairshell(i)%nfc_jij
                if ( j .gt. 0 ) l=l+j*storage_size(sl%magpairshell(i)%relfc_jij)/8
                if ( j .gt. 0 ) l=l+j*storage_size(sl%magpairshell(i)%relfcind_jij)/8
                j=sl%magpairshell(i)%nfc_tij
                if ( j .gt. 0 ) l=l+j*storage_size(sl%magpairshell(i)%relfc_tij)/8
                if ( j .gt. 0 ) l=l+j*storage_size(sl%magpairshell(i)%relfcind_tij)/8
            enddo
        endif
        mem=mem+l
        ! atoms
        l=0
        if ( allocated(sl%uc) ) then
            l=l+sl%nuc*storage_size(sl%uc)/8
            do i=1,sl%nuc
                l=l+storage_size(sl%uc(i)%singlet)/8
                if ( sl%uc(i)%npair .gt. 0 ) l=l+sl%uc(i)%npair*storage_size(sl%uc(i)%pair)/8
                if ( sl%uc(i)%ntriplet .gt. 0 ) l=l+sl%uc(i)%ntriplet*storage_size(sl%uc(i)%triplet)/8
                if ( sl%uc(i)%nquartet .gt. 0 ) l=l+sl%uc(i)%nquartet*storage_size(sl%uc(i)%quartet)/8
                if ( sl%uc(i)%nmagpair .gt. 0 ) l=l+sl%uc(i)%nmagpair*storage_size(sl%uc(i)%magpair)/8
            enddo
        endif
        if ( allocated(sl%ss) ) then
            l=l+sl%nss*storage_size(sl%ss)/8
            do i=1,sl%nss
                l=l+storage_size(sl%ss(i)%singlet)/8
                if ( sl%ss(i)%npair .gt. 0 ) l=l+sl%ss(i)%npair*storage_size(sl%ss(i)%pair)/8
                if ( sl%ss(i)%ntriplet .gt. 0 ) l=l+sl%ss(i)%ntriplet*storage_size(sl%ss(i)%triplet)/8
                if ( sl%ss(i)%nquartet .gt. 0 ) l=l+sl%ss(i)%nquartet*storage_size(sl%ss(i)%quartet)/8
                if ( sl%ss(i)%nmagpair .gt. 0 ) l=l+sl%ss(i)%nmagpair*storage_size(sl%ss(i)%magpair)/8
            enddo
        endif
        mem=mem+l
    end function
#endif

end module
