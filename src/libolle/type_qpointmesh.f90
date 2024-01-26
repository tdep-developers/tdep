module type_qpointmesh
!!
!! Extensive module that handles meshes and paths in reciprocal space, generated in a variety of ways.
!! Uncharacteristically, this module exposes quite a few things. But there was really no other way, and
!! it looks far more complicated than it actually is. If you don't care too much, you can forget all about
!! how the mesh is generated and just use it.
!!
use konstanter, only: r8,i8,lo_iou,lo_huge,lo_hugeint,lo_tol,lo_sqtol,lo_pi,lo_twopi,lo_exitcode_param,&
    lo_exitcode_symmetry,lo_exitcode_io,lo_degenvector,lo_freqtol,lo_bohr_to_A,lo_radiantol
use mpi_wrappers, only: lo_mpi_helper,lo_stop_gracefully
use lo_memtracker, only: lo_mem_helper
use gottochblandat, only: walltime,tochar,lo_sqnorm,lo_clean_fractional_coordinates,&
    lo_chop,lo_index_in_periodic_array,lo_unsigned_tetrahedron_volume,lo_gauss,&
    open_file,lo_determ,lo_trueNtimes,lo_progressbar,lo_progressbar_init
use geometryfunctions, only: lo_linesegment,lo_plane,lo_inscribed_sphere_in_box,lo_angle_between_vectors
use type_crystalstructure, only: lo_crystalstructure
use type_symmetryoperation, only: lo_operate_on_vector
use hdf5

implicit none
private

! derived types
public :: lo_qpoint
public :: lo_qpoint_irrwedge
public :: lo_qpoint_mesh
public :: lo_qpoint_grid
public :: lo_monkhorst_pack_mesh
public :: lo_fft_mesh
public :: lo_wedge_mesh
public :: lo_commensurate_mesh
public :: lo_bandstructure

! some routines
public :: lo_read_qmesh_from_file
public :: lo_generate_qmesh
public :: lo_integration_weights_for_one_tetrahedron
public :: lo_get_small_group_of_qpoint
public :: lo_get_wedge_mesh
public :: lo_LV_tetrahedron_weights
public :: lo_LV_complex_tetrahedron_weights
public :: lo_LV_tetrahedron_heaviside
public :: lo_LV_tetrahedron_fermi

!> A tetrahedron in a q-grid mesh
type lo_qptetrahedron
    !> indices to irreducible points
    integer, dimension(4) :: irreducible_index=-lo_hugeint
    !> indices to full points
    integer, dimension(4) :: full_index=-lo_hugeint
    !> integration weight of this tetrahedron. Sums to 1.
    real(r8) :: integration_weight=-lo_huge
    !> center of mass of tetrahedron
    real(r8), dimension(3) :: center_of_mass=-lo_huge
    !> lattice vector shifts of tetrahedron corners
    integer, dimension(3,4) :: lattice_vector_shift=-lo_hugeint
end type

!> A general q-point
type lo_qpoint
    !> Cartesian coordinate, folded into the first BZ
    real(r8), dimension(3) :: r=lo_huge
    !> how many operations are there in the small group of q?
    integer :: n_invariant_operation=-lo_hugeint
    !> Which are these operations?
    integer, dimension(:), allocatable :: invariant_operation
    contains
        !> measure the size in memory
        procedure :: size_in_mem=>qpoint_size_in_mem
end type

!> A q-point on a grid in the BZ
type, extends(lo_qpoint) :: lo_qpoint_fullgrid
    !> which irreducible k-point is it equivalent to
    integer :: irreducible_index=-lo_hugeint
    !> which operation transforms the irreducible to this point
    integer :: operation_from_irreducible=-lo_hugeint
    !> What is the integration weight. Normalized wuch that sum(weight)=1
    real(r8) :: integration_weight=-lo_huge
    !> Radius of sphere that gives the same volume as the volume element the q-point represents, not taking weight into account
    real(r8) :: radius=-lo_huge
end type

!> A q-point in the irreducible wedge
type, extends(lo_qpoint) :: lo_qpoint_irrwedge
    !> The relative weight. Normalized so that sum(weight)=1
    real(r8) :: integration_weight=lo_huge
    !> which index in the full grid does it correspond to?
    integer :: full_index=-lo_hugeint
    !> how many points in the full mesh can this k-point transform to
    integer :: n_full_point=-lo_hugeint
    !> which full points can it transform to. Negative index means add inversion!!
    integer, dimension(:), allocatable :: index_full_point
    !> which symmetry operation does the transformation to the index specified above
    integer, dimension(:), allocatable :: operation_full_point
    !> Radius of sphere that gives the same volume as the volume element the q-point represents, not taking weight into account
    real(r8) :: radius=-lo_huge
end type

!> A q-point along a path in the Brillouin zone
type, extends(lo_qpoint) :: lo_qpoint_bandstructure
    !> absolute position, as generated
    real(r8), dimension(3) :: r_abs=-lo_huge
    !> which of the path segments does it sit on?
    integer :: path=-lo_hugeint
end type

!> A general q-point mesh in the Brillouin zone.
type :: lo_qpoint_mesh
    !> number of irreducible k-points on this rank
    integer :: n_irr_point=-lo_hugeint
    !> number of full k-points on this rank
    integer :: n_full_point=-lo_hugeint
    !> irreducible q-points
    type(lo_qpoint_irrwedge), allocatable, dimension(:) :: ip
    !> all the points in the BZ
    type(lo_qpoint_fullgrid), allocatable, dimension(:) :: ap
    !> number of irreducible tetrahedron
    integer :: n_irr_tet=-lo_hugeint
    !> number of full tetrahedron
    integer :: n_full_tet=-lo_hugeint
    !> irreducible tetrahedrons
    type(lo_qptetrahedron), allocatable, dimension(:) :: it
    !> full set of tetrahedrons
    type(lo_qptetrahedron), allocatable, dimension(:) :: at

    ! Some private things

    !> scaled reciprocal basis for adaptive gaussian smearing
    real(r8), dimension(3,3), private :: scaledrecbasis
    !> does the property on the mesh have time reversal symmetry
    logical :: timereversal=.false.
    !> Which operations are ok to test with? Only relevant in weird cases
    logical, dimension(:), allocatable :: operationok
    contains
        !> Get the adaptive gaussian smearing parameter
        procedure :: smearingparameter
        !> Updated way of determining smearing parameter
        procedure, nopass :: adaptive_sigma
        !> Write it to file
        procedure :: write_to_file
        !> Write human readable information to an hdf5 handle
        procedure :: write_metadata_to_hdf5
        !> update the integration weights
        procedure :: update_integration_weight=>voronoi_integration_weights
        !> Measure size in memory, in bytes
        procedure, nopass :: size_in_mem=>qmesh_size_in_mem
        !> Destroy
        procedure, nopass :: destroy
end type

!> A q-point mesh that is a regular grid, but not completely clear which kind of regular grid.
type, abstract, extends(lo_qpoint_mesh) :: lo_qpoint_grid
    !> dimensions of grid
    integer, dimension(3) :: griddensity=-lo_hugeint
    !> conversion from linear index to three indices
    integer, dimension(:,:), allocatable :: ind2gridind
    !> conversion from three indices to a linear
    integer, dimension(:,:,:), allocatable :: gridind2ind
    contains
        procedure(coordinate_from_index), deferred :: coordinate_from_index
        procedure(index_from_coordinate), deferred :: index_from_coordinate
        procedure(is_point_on_grid), deferred :: is_point_on_grid
end type

! Abstract interfaces
interface
    pure function coordinate_from_index(qp,i) result(r)
        import :: lo_qpoint_grid
        import :: r8
        class(lo_qpoint_grid), intent(in) :: qp
        integer, dimension(3), intent(in) :: i
        real(r8), dimension(3) :: r
    end function
    pure function index_from_coordinate(qp,r) result(i)
        import :: lo_qpoint_grid
        import :: r8
        class(lo_qpoint_grid), intent(in) :: qp
        real(r8), dimension(3), intent(in) :: r
        integer, dimension(3) :: i
    end function
    pure function is_point_on_grid(qp,r) result(ongrid)
        import :: lo_qpoint_grid
        import :: r8
        class(lo_qpoint_grid), intent(in) :: qp
        real(r8), dimension(3), intent(in) :: r
        logical :: ongrid
    end function
end interface

!> A monkhorst pack mesh
type, extends(lo_qpoint_grid) :: lo_monkhorst_pack_mesh
    contains
        procedure :: coordinate_from_index => mp_coordinate_from_index
        procedure :: index_from_coordinate => mp_index_from_coordinate
        procedure :: is_point_on_grid => mp_is_point_on_grid
end type
!> q-grid commensurate with the Fourier transform of a cell
type, extends(lo_qpoint_grid) :: lo_fft_mesh
    contains
        procedure :: coordinate_from_index => fft_coordinate_from_index
        procedure :: index_from_coordinate => fft_index_from_coordinate
        procedure :: is_point_on_grid => fft_is_point_on_grid
end type

!> A q-point mesh generated from the irreducible wedge
type, extends(lo_qpoint_mesh) :: lo_wedge_mesh
    contains
        !> chop into tetrahedrons
        procedure :: tesselate_wedge_mesh
end type

!> A q-point mesh built from a commensurate supercell
type, extends(lo_qpoint_mesh) :: lo_commensurate_mesh
end type

!> A bandstructure in the Brillouin zone
type lo_bandstructure
    !> number of points in total
    integer :: n_point=-lo_hugeint
    !> number of path segments
    integer :: n_path=-lo_hugeint
    !> number of points along each path
    integer :: n_point_per_path=-lo_hugeint
    !> stride (calculate values every stride points, then interpolate)
    integer :: stride=-lo_hugeint
    !> q-points on the path
    type(lo_qpoint_bandstructure), dimension(:), allocatable :: q

    !> if using a stride, what q-points should be calculated exactly?
    integer, dimension(:,:), allocatable :: stride_q_ind

    ! Mostly auxiliary information below

    !> line segments
    type(lo_linesegment), dimension(:), allocatable :: segment
    !> symbols for the q-points on at the start of the paths
    character(len=100), dimension(:), allocatable :: symb_q_start
    !> symbols for the q-points on at the end of the paths
    character(len=100), dimension(:), allocatable :: symb_q_end
    !> q-axis
    real(r8), dimension(:), allocatable :: q_axis
    !> q-axis tick marks
    real(r8), dimension(:), allocatable :: q_axis_ticks
    !> q-axis tick labels
    character(len=100), dimension(:), allocatable :: q_axis_tick_labels
    contains
        !> read a custom path from a file
        procedure :: read_path_from_file
        !> get the standard path
        procedure :: standardpath
end type

! Interfaces to bandstructure things
interface
    module subroutine interpolate_kpoint_path(bs,p,mw,verbosity)
        class(lo_bandstructure), intent(inout) :: bs
        type(lo_crystalstructure), intent(in) :: p
        type(lo_mpi_helper), intent(inout) :: mw
        integer, intent(in) :: verbosity
    end subroutine
    module subroutine read_path_from_file(bs,p,mw,verbosity)
        class(lo_bandstructure), intent(inout) :: bs
        type(lo_crystalstructure), intent(inout) :: p
        type(lo_mpi_helper), intent(inout) :: mw
        integer, intent(in) :: verbosity
    end subroutine
    module subroutine standardpath(bs,p,mw,verbosity)
        class(lo_bandstructure), intent(inout) :: bs
        type(lo_crystalstructure), intent(in) :: p
        type(lo_mpi_helper), intent(inout) :: mw
        integer, intent(in) :: verbosity
    end subroutine
end interface

! Interfaces to integration weights
interface
    module pure function lo_LV_complex_tetrahedron_weights(ein,z,tol) result(weights)
        real(r8), dimension(4), intent(in) :: ein
        complex(r8), intent(in) :: z
        real(r8), intent(in) :: tol
        complex(r8), dimension(4) :: weights
    end function
    module pure function lo_LV_tetrahedron_weights(ein,z,tol,sigma) result(weights)
        real(r8), intent(in), dimension(4) :: ein
        real(r8), intent(in) :: z
        real(r8), intent(in) :: tol
        real(r8), intent(in) :: sigma
        real(r8), dimension(4) :: weights
    end function
    module pure function lo_LV_tetrahedron_heaviside(ein,z,tol) result(weight)
        real(r8), intent(in), dimension(4) :: ein
        real(r8), intent(in) :: z
        real(r8), intent(in) :: tol
        real(r8) :: weight
    end function
    module pure function lo_LV_tetrahedron_fermi(ein,z,temperature,tol) result(weight)
        real(r8), intent(in), dimension(4) :: ein
        real(r8), intent(in) :: z
        real(r8), intent(in) :: temperature
        real(r8), intent(in) :: tol
        real(r8) :: weight
    end function
    module pure function smearingparameter(qp,gradient,sigma,adaptiveparameter) result(w)
        class(lo_qpoint_mesh), intent(in) :: qp
        real(r8), dimension(3), intent(in) :: gradient
        real(r8), intent(in) :: sigma
        real(r8), intent(in) :: adaptiveparameter
        real(r8) :: w
    end function
    module pure function adaptive_sigma(radius,gradient,default_sigma,scale) result(sigma)
        real(r8), intent(in) :: radius
        real(r8), dimension(3), intent(in) :: gradient
        real(r8), intent(in) :: default_sigma
        real(r8), intent(in) :: scale
        real(r8) :: sigma
    end function
    module subroutine lo_integration_weights_for_one_tetrahedron(tet,c,eps,sigma,tol,weights,blochlcorrections,hweights)
        type(lo_qptetrahedron), intent(in) :: tet
        real(r8), dimension(4), intent(in) :: c
        real(r8), intent(in) :: eps
        real(r8), intent(in) :: sigma
        real(r8), intent(in) :: tol
        real(r8), intent(out), dimension(4) :: weights
        logical, intent(in), optional :: blochlcorrections
        real(r8), intent(out), dimension(4), optional :: hweights
    end subroutine
    module subroutine voronoi_integration_weights(qp,p,mw,mem,integration_error,verbosity)
        class(lo_qpoint_mesh), intent(inout) :: qp
        type(lo_crystalstructure), intent(in) :: p
        type(lo_mpi_helper), intent(inout) :: mw
        type(lo_mem_helper), intent(inout) :: mem
        real(r8), intent(out) :: integration_error
        integer, intent(in) :: verbosity
    end subroutine
end interface

! interfaces to type_qpointmesh_io
interface
    module subroutine write_to_file(qp,p,filename,mem,verbosity)
        class(lo_qpoint_mesh), intent(in) :: qp
        type(lo_crystalstructure), intent(in) :: p
        character(len=*), intent(in) :: filename
        type(lo_mem_helper), intent(inout) :: mem
        integer, intent(in) :: verbosity
    end subroutine
    module subroutine lo_read_qmesh_from_file(qp,p,filename,mem,verbosity)
        class(lo_qpoint_mesh), intent(out), allocatable :: qp
        type(lo_crystalstructure), intent(inout) :: p
        character(len=*), intent(in) :: filename
        type(lo_mem_helper), intent(inout) :: mem
        integer, intent(in) :: verbosity
    end subroutine
    module subroutine write_metadata_to_hdf5(qp,filename,input_id)
        class(lo_qpoint_mesh), intent(in) :: qp
        character(len=*), intent(in), optional :: filename
        integer(HID_T), intent(in), optional :: input_id
    end subroutine
end interface

! interfaces to lo_qpointmesh_commensurate
interface
    module subroutine lo_generate_commensurate_qmesh(qp,uc,supercellmatrix,mw,mem,verbosity)
        class(lo_commensurate_mesh), intent(out) :: qp
        type(lo_crystalstructure), intent(in) :: uc
        real(r8), dimension(3,3), intent(in) :: supercellmatrix
        type(lo_mpi_helper), intent(inout) :: mw
        type(lo_mem_helper), intent(inout) :: mem
        integer, intent(in) :: verbosity
    end subroutine
    module subroutine get_supercell_matrix(uc,griddensity,supercellmatrix)
        type(lo_crystalstructure), intent(in) :: uc
        integer, dimension(3), intent(in) :: griddensity
        real(r8), dimension(3,3), intent(out) :: supercellmatrix
    end subroutine
end interface

! Interfaces to type_qpointmesh_gridgeneration
interface
    module subroutine build_grid(qp,p,mw,mem,verbosity,nosymmetry)
        class(lo_qpoint_grid), intent(inout) :: qp
        type(lo_crystalstructure), intent(in) :: p
        type(lo_mpi_helper), intent(inout) :: mw
        type(lo_mem_helper), intent(inout) :: mem
        integer, intent(in) :: verbosity
        logical, intent(in) :: nosymmetry
    end subroutine
    module subroutine tesselate_qgrid(qp,p,mw,mem,verbosity)
        class(lo_qpoint_grid), intent(inout) :: qp
        type(lo_crystalstructure), intent(in) :: p
        type(lo_mpi_helper), intent(inout) :: mw
        type(lo_mem_helper), intent(inout) :: mem
        integer, intent(in) :: verbosity
    end subroutine
    module pure function mp_coordinate_from_index(qp,i) result(r)
        class(lo_monkhorst_pack_mesh), intent(in) :: qp
        integer, dimension(3), intent(in) :: i
        real(r8), dimension(3) :: r
    end function
    module pure function mp_index_from_coordinate(qp,r) result(i)
        class(lo_monkhorst_pack_mesh), intent(in) :: qp
        real(r8), dimension(3), intent(in) :: r
        integer, dimension(3) :: i
    end function
    module pure function mp_is_point_on_grid(qp,r) result(ongrid)
        class(lo_monkhorst_pack_mesh), intent(in) :: qp
        real(r8), dimension(3), intent(in) :: r
        logical :: ongrid
    end function
    module pure function fft_is_point_on_grid(qp,r) result(ongrid)
        class(lo_fft_mesh), intent(in) :: qp
        real(r8), dimension(3), intent(in) :: r
        logical :: ongrid
    end function
    module pure function fft_index_from_coordinate(qp,r) result(i)
        class(lo_fft_mesh), intent(in) :: qp
        real(r8), dimension(3), intent(in) :: r
        integer, dimension(3) :: i
    end function
    module pure function fft_coordinate_from_index(qp,i) result(r)
        class(lo_fft_mesh), intent(in) :: qp
        integer, dimension(3), intent(in) :: i
        real(r8), dimension(3) :: r
    end function
end interface

! Interfaces to type_qpointmesh_wedgegeneration
interface
    module subroutine lo_get_wedge_mesh(qp,uc,griddensity,mw,mem,verbosity)
        type(lo_wedge_mesh), intent(inout) :: qp
        type(lo_crystalstructure), intent(in) :: uc
        integer, dimension(3), intent(in) :: griddensity
        type(lo_mpi_helper), intent(inout) :: mw
        type(lo_mem_helper), intent(inout) :: mem
        integer, intent(in) :: verbosity
    end subroutine
    module subroutine tesselate_wedge_mesh(qp,uc,wedgevertices,prunetol,splittol,mw,mem,verbosity)
        class(lo_wedge_mesh), intent(inout) :: qp
        type(lo_crystalstructure), intent(in) :: uc
        real(r8), dimension(:,:), allocatable, intent(inout) :: wedgevertices
        real(r8), intent(in) :: prunetol
        real(r8), intent(in) :: splittol
        type(lo_mpi_helper), intent(inout) :: mw
        type(lo_mem_helper), intent(inout) :: mem
        integer, intent(in) :: verbosity
    end subroutine
end interface

contains

!> Generate q-point meshes.
subroutine lo_generate_qmesh(qp,uc,griddensity,meshtype,timereversal,headrankonly,mw,mem,verbosity,nosym)
    !> q-point mesh
    class(lo_qpoint_mesh), intent(out), allocatable :: qp
    !> crystal structure
    type(lo_crystalstructure), intent(inout) :: uc
    !> grid density, if applicable
    integer, dimension(3) :: griddensity
    !> what kind of grid
    character(len=*), intent(in) :: meshtype
    !> should the mesh reflect time-reveral symmetry?
    logical, intent(in) :: timereversal
    !> create the mesh on the head rank only (rank 0 on the provided communicator)?
    logical, intent(in) :: headrankonly
    !> MPI helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> Memory tracker
    type(lo_mem_helper), intent(inout) :: mem
    !> How much information should be printed to stdout?
    integer, intent(in) :: verbosity
    !> Skip all symmetry stuff
    logical, intent(in), optional :: nosym

    real(r8) :: timer
    logical :: nosymmetry

    ! First of all, start the timer
    timer=walltime()

    ! Insert a tick for memory tracking
    call mem%tick()

    ! Set some initial things
    init: block
        if ( verbosity .gt. 0 ) then
            write(lo_iou,*) ' '
            write(lo_iou,*) 'Generating Brillouin zone mesh'
        endif

        ! Check that we have enough information to create the mesh.
        if ( uc%info%havespacegroup .eqv. .false. ) then
            call lo_stop_gracefully(['The space group needs to be determined before generating meshes.'],lo_exitcode_param,__FILE__,__LINE__,mw%comm)
        else
            if ( uc%sym%timereversal .neqv. timereversal ) then
                call lo_stop_gracefully(['Conflicting instructions regarding time-reversal symmetry!'],lo_exitcode_param,__FILE__,__LINE__,mw%comm)
            endif
        endif

        ! Allocate the kind of mesh on each rank. When we are finalizing later on
        ! it will be deallocated except for where it is supposed to be. We create
        ! the mesh on a single rank, always, to save as much memory as possible.
        select case(trim(meshtype))
            case('monkhorst')
                ! Requires more symmetry than the minimal
                if ( uc%info%havebz .eqv. .false. ) then
                    call lo_stop_gracefully(['Need the Brillouin zone to generate q-mesh'],lo_exitcode_param,__FILE__,__LINE__,mw%comm)
                endif
                allocate(lo_monkhorst_pack_mesh::qp)
            case('fft')
                ! Requires more symmetry than the minimal
                if ( uc%info%havebz .eqv. .false. ) then
                    call lo_stop_gracefully(['Need the Brillouin zone to generate q-mesh'],lo_exitcode_param,__FILE__,__LINE__,mw%comm)
                endif
                allocate(lo_fft_mesh::qp)
            case('wedge')
                ! This requires slightly more information about the symmetry.
                if ( uc%info%havewedge .eqv. .false. ) then
                    call lo_stop_gracefully(['The irreducible wedge needs to be determined before generating meshes.'],lo_exitcode_param,__FILE__,__LINE__,mw%comm)
                endif
                allocate(lo_wedge_mesh::qp)
            case('commensurate')
                allocate(lo_commensurate_mesh::qp)
            case default
                call lo_stop_gracefully(['Unknown mesh type'],lo_exitcode_param,__FILE__,__LINE__)
        end select

        ! And make a note that we have decided on the time reversal thing.
        qp%timereversal=timereversal

        ! Here is a very ugly thing, but it is mostly for debugging: temporarily switch
        ! off any symmetries. Think hard before you use it.
        if ( present(nosym) ) then
            nosymmetry=nosym
        else
            nosymmetry=.false.
        endif
    end block init

    ! Now build the actual mesh. It is done quite differently depending
    ! on what kind of mesh we want.
    buildmesh: block
        type(lo_mpi_helper) :: ml
        real(r8), dimension(3,3) :: supercellmatrix
        real(r8) :: interr

        if ( headrankonly ) then
            ! Split the communicator into tiny pieces
            call mw%split(ml,mw%r+1,__FILE__,__LINE__)
        endif

        !@TODO Here I can insert my thing: In case I only want to
        ! get stuff on the head rank, split the communicator and
        ! use only a single communicator on a single rank and pass
        ! that through to the routines.

        if ( headrankonly ) then
            if ( mw%r .eq. 0 ) then
                ! Unusual case, note non-standard communicator
                select type(qp)
                type is (lo_monkhorst_pack_mesh)
                    ! Get a normal MP mesh
                    qp%griddensity=griddensity
                    call build_grid(qp,uc,ml,mem,verbosity,nosymmetry)
                    ! Tesselate the grid into tetrahedrons
                    call tesselate_qgrid(qp,uc,ml,mem,verbosity)
                type is (lo_fft_mesh)
                    ! Get a mesh that is commensurate with a discrete Fourier transform
                    qp%griddensity=griddensity
                    call build_grid(qp,uc,ml,mem,verbosity,nosymmetry)
                    ! And tesselate it
                    call tesselate_qgrid(qp,uc,ml,mem,verbosity)
                type is (lo_wedge_mesh)
                    ! Generate the fancy CGAL wedge-based mesh.
                    call lo_get_wedge_mesh(qp,uc,griddensity,ml,mem,verbosity)
                    ! Update weights
                    call voronoi_integration_weights(qp,uc,ml,mem,interr,verbosity)
                type is(lo_commensurate_mesh)
                    ! In case we have a grid-density as input, we have to guesstimate the best
                    ! supercell from this. Or the supercell could be input, also an option.
                    call get_supercell_matrix(uc,griddensity,supercellmatrix)
                    call lo_generate_commensurate_qmesh(qp,uc,supercellmatrix,mw,mem,verbosity)
                    call lo_stop_gracefully(['COMMENSURATE MESHES NOT QUITE DONE.'],lo_exitcode_param,__FILE__,__LINE__,mw%comm)
                class default
                    call lo_stop_gracefully(['unknown mesh type'],lo_exitcode_param,__FILE__,__LINE__,mw%comm)
                end select
            endif
        else
            ! Normal case, normal communicator.
            select type(qp)
            type is (lo_monkhorst_pack_mesh)
                ! Get a normal MP mesh
                qp%griddensity=griddensity
                call build_grid(qp,uc,mw,mem,verbosity,nosymmetry)
                ! Tesselate the grid into tetrahedrons
                call tesselate_qgrid(qp,uc,mw,mem,verbosity)
            type is (lo_fft_mesh)
                ! Get a mesh that is commensurate with a discrete Fourier transform
                qp%griddensity=griddensity
                call build_grid(qp,uc,mw,mem,verbosity,nosymmetry)
                ! And tesselate it
                call tesselate_qgrid(qp,uc,mw,mem,verbosity)
            type is (lo_wedge_mesh)
                ! Generate the fancy CGAL wedge-based mesh.
                call lo_get_wedge_mesh(qp,uc,griddensity,mw,mem,verbosity)
                ! Update weights
                call voronoi_integration_weights(qp,uc,mw,mem,interr,verbosity)
            type is(lo_commensurate_mesh)
                ! In case we have a grid-density as input, we have to guesstimate the best
                ! supercell from this. Or the supercell could be input, also an option.
                call get_supercell_matrix(uc,griddensity,supercellmatrix)
                call lo_generate_commensurate_qmesh(qp,uc,supercellmatrix,mw,mem,verbosity)
                call lo_stop_gracefully(['COMMENSURATE MESHES NOT QUITE DONE.'],lo_exitcode_param,__FILE__,__LINE__,mw%comm)
            class default
                call lo_stop_gracefully(['unknown mesh type'],lo_exitcode_param,__FILE__,__LINE__,mw%comm)
            end select
        endif

        ! Store the scaled basis for adaptive gaussian smearing
        ! Should be done in a better place, perhaps. Also should
        ! be done better in general.
        select type(qp)
        type is (lo_monkhorst_pack_mesh)
            qp%scaledrecbasis=uc%reciprocal_latticevectors*lo_twopi
            qp%scaledrecbasis(:,1)=qp%scaledrecbasis(:,1)/griddensity(1)
            qp%scaledrecbasis(:,2)=qp%scaledrecbasis(:,2)/griddensity(2)
            qp%scaledrecbasis(:,3)=qp%scaledrecbasis(:,3)/griddensity(3)
        type is (lo_fft_mesh)
            qp%scaledrecbasis=uc%reciprocal_latticevectors*lo_twopi
            qp%scaledrecbasis(:,1)=qp%scaledrecbasis(:,1)/griddensity(1)
            qp%scaledrecbasis(:,2)=qp%scaledrecbasis(:,2)/griddensity(2)
            qp%scaledrecbasis(:,3)=qp%scaledrecbasis(:,3)/griddensity(3)
        type is (lo_wedge_mesh)
            ! Do nothing
        type is(lo_commensurate_mesh)
            ! Do nothing
        class default
            call lo_stop_gracefully(['unknown mesh type'],lo_exitcode_param,__FILE__,__LINE__,mw%comm)
        end select

        ! In case we only wanted the mesh on the head rank, sort some thing out.
        if ( headrankonly ) then
            ! Destroy all meshes not on head rank. Makes sense.
            if ( mw%r .ne. 0 ) then
                if ( allocated(qp) ) deallocate(qp)
            endif
            ! Also destroy the temporary communicators I've created, they are no longer needed.
            call mw%barrier(__FILE__,__LINE__)
            call ml%free(__FILE__,__LINE__)
        endif
    end block buildmesh

    ! Check that memory is cleared propely? There should be nothing trailing behind at this point.
    call mem%tock(__FILE__,__LINE__,mw%comm)

    ! And we are all done
    if ( verbosity .gt. 0 ) write(lo_iou,*) 'Done building mesh (',tochar(walltime()-timer),'s)'
end subroutine

!> Generate the small point group of q
subroutine lo_get_small_group_of_qpoint(qpoint,p)
    !> q-point
    class(lo_qpoint), intent(inout) :: qpoint
    !> crystal structure
    type(lo_crystalstructure), intent(in) :: p

    real(r8), dimension(3) :: v0,v1
    integer, dimension(p%sym%n*2) :: ind
    integer :: i,ctr

    ind=0
    ctr=0
    v0=qpoint%r
    ! First check perfectly normal operations
    do i=1,p%sym%n
        v1=v0-lo_operate_on_vector(p%sym%op(i),v0,reciprocal=.true.)
        v1=matmul(p%inv_reciprocal_latticevectors,v1)
        v1=v1-anint(v1)
        if ( lo_sqnorm(v1) .lt. lo_sqtol ) then
            ctr=ctr+1
            ind(ctr)=i
        endif
    enddo

    ! ! Then see if we can do the time-reversal thing?
    ! do i=1,p%sym%n
    !     v1=v0+lo_operate_on_vector(p%sym%op(i),v0,reciprocal=.true.)
    !     v1=matmul(p%inv_reciprocal_latticevectors,v1)
    !     v1=v1-anint(v1)
    !     if ( lo_sqnorm(v1) .lt. lo_sqtol ) then
    !         ctr=ctr+1
    !         ind(ctr)=-i
    !     endif
    ! enddo

    ! Store list of operations
    qpoint%n_invariant_operation=ctr
    if ( allocated(qpoint%invariant_operation) ) deallocate(qpoint%invariant_operation)
    allocate(qpoint%invariant_operation(ctr))
    qpoint%invariant_operation=ind(1:ctr)
end subroutine

!> destroy and free all memory
subroutine destroy(qp)
    !> dispersions
    class(lo_qpoint_mesh), allocatable, intent(inout) :: qp

    integer :: i

    if ( allocated(qp) ) then
        if ( allocated(qp%operationok) ) deallocate(qp%operationok)
        if ( allocated(qp%ip) ) then
            do i=1,size(qp%ip)
                if ( allocated(qp%ip(i)%invariant_operation) ) deallocate(qp%ip(i)%invariant_operation)
                if ( allocated(qp%ip(i)%index_full_point) ) deallocate(qp%ip(i)%index_full_point)
                if ( allocated(qp%ip(i)%operation_full_point) ) deallocate(qp%ip(i)%operation_full_point)
            enddo
            deallocate(qp%ip)
        endif
        if ( allocated(qp%ap) ) then
            do i=1,size(qp%ap)
                if ( allocated(qp%ap(i)%invariant_operation) ) deallocate(qp%ap(i)%invariant_operation)
            enddo
            deallocate(qp%ap)
        endif
        if ( allocated(qp%it) ) deallocate(qp%it)
        if ( allocated(qp%at) ) deallocate(qp%at)
        deallocate(qp)
    endif
end subroutine

!> measure size in memory, in bytes
function qmesh_size_in_mem(qp) result(mem)
    !> dispersions
    class(lo_qpoint_mesh), allocatable, intent(in) :: qp
    !> memory in bytes
    integer :: mem,i

    ! No memory if not allocated yet
    if ( allocated(qp) .eqv. .false. ) then
        mem=0
        return
    endif

    mem=0
    select type(qp)
    type is(lo_wedge_mesh)
        mem=mem+storage_size(qp)
    type is(lo_fft_mesh)
        mem=mem+storage_size(qp)
    type is(lo_monkhorst_pack_mesh)
        mem=mem+storage_size(qp)
    type is(lo_commensurate_mesh)
        mem=mem+storage_size(qp)
    end select

    if ( allocated(qp%operationok) ) mem=mem+storage_size(qp%operationok)*size(qp%operationok)
    if ( allocated(qp%it) ) then
        mem=mem+storage_size(qp%it)*size(qp%it)
    endif
    if ( allocated(qp%at) ) then
        mem=mem+storage_size(qp%at)*size(qp%at)
    endif
    mem=mem/8

    if ( allocated(qp%ip) ) then
        do i=1,size(qp%ip)
            mem=mem+qp%ip(i)%size_in_mem()
        enddo
    endif
    if ( allocated(qp%ap) ) then
        do i=1,size(qp%ap)
            mem=mem+qp%ap(i)%size_in_mem()
        enddo
    endif
end function

!> measure size in memory, in bytes
function qpoint_size_in_mem(q) result(mem)
    !> q-point
    class(lo_qpoint), intent(in) :: q
    !> memory in bytes
    integer :: mem

    mem=0
    select type(q)
    type is(lo_qpoint)
        mem=mem+storage_size(q)
        if ( allocated(q%invariant_operation) ) mem=mem+storage_size(q%invariant_operation)*size(q%invariant_operation)
    type is(lo_qpoint_fullgrid)
        mem=mem+storage_size(q)
        if ( allocated(q%invariant_operation) ) mem=mem+storage_size(q%invariant_operation)*size(q%invariant_operation)
    type is(lo_qpoint_irrwedge)
        mem=mem+storage_size(q)
        if ( allocated(q%invariant_operation) ) mem=mem+storage_size(q%invariant_operation)*size(q%invariant_operation)
        if ( allocated(q%index_full_point) ) mem=mem+storage_size(q%index_full_point)*size(q%index_full_point)
        if ( allocated(q%operation_full_point) ) mem=mem+storage_size(q%operation_full_point)*size(q%operation_full_point)
    type is(lo_qpoint_bandstructure)
        mem=mem+storage_size(q)
        if ( allocated(q%invariant_operation) ) mem=mem+storage_size(q%invariant_operation)*size(q%invariant_operation)
    end select
    mem=mem/8
end function

! !> Sanity test a q-point mesh, to see if something went wrong somewhere
! subroutine sanitytest(qp,p)
!     !> the mesh
!     class(lo_qpoint_mesh), intent(in) :: qp
!     !> crystal structure
!     type(lo_crystalstructure), intent(in) :: p
!
!     integer :: i,j,k,l
!     real(r8), dimension(3) :: v0,v1,v2
!     real(r8) :: f0
!     logical :: err
!
!     write(*,*) '... starting sanity test of q-point mesh'
!
!     ! Test the symmetry operations
!     err=.false.
!     do i=1,qp%n_full_point
!         j=qp%ap(i)%irreducible_index
!         k=qp%ap(i)%operation_from_irreducible
!         ! The idea is that
!         ! operation*irreducible = qpoint i
!         ! I should really test that, both in Cartesian and fractional coordinates. Start with Cartesian.
!         if ( k .gt. 0 ) then
!             v0=qp%ap(i)%r
!             v1=lo_operate_on_vector(p%sym%op(k),qp%ip(j)%r,reciprocal=.true.,inverse=.false.)
!             v2=matmul(p%inv_reciprocal_latticevectors,v0-v1)
!             v2=lo_clean_fractional_coordinates(v2+0.5_r8)-0.5_r8
!             if ( lo_sqnorm(v2) .gt. lo_sqrectol ) then
!                 write(*,*) 'ERROR, bad symmetry operation in q-point mesh, Cartesian coordinates'
!                 err=.true.
!             endif
!             ! Then fractional
!             v0=matmul(p%inv_reciprocal_latticevectors,qp%ap(i)%r)
!             v1=matmul(p%inv_reciprocal_latticevectors,qp%ip(j)%r)
!             v1=lo_operate_on_vector(p%sym%op(k),v1,reciprocal=.true.,fractional=.true.,inverse=.false.)
!             v2=lo_clean_fractional_coordinates(v0-v1+0.5_r8)-0.5_r8
!             if ( lo_sqnorm(v2) .gt. lo_sqrectol ) then
!                 write(*,*) 'ERROR, bad symmetry operation in q-point mesh, fractional coordinates'
!                 err=.true.
!             endif
!         else
!             v0=qp%ap(i)%r
!             v1=lo_operate_on_vector(p%sym%op(-k),qp%ip(j)%r,reciprocal=.true.,inverse=.false.)
!             v2=matmul(p%inv_reciprocal_latticevectors,v0+v1)
!             v2=lo_clean_fractional_coordinates(v2+0.5_r8)-0.5_r8
!             if ( lo_sqnorm(v2) .gt. lo_sqrectol ) then
!                 write(*,*) 'ERROR, bad symmetry operation in q-point mesh, Cartesian coordinates'
!                 err=.true.
!             endif
!             v0=matmul(p%inv_reciprocal_latticevectors,qp%ap(i)%r)
!             v1=matmul(p%inv_reciprocal_latticevectors,qp%ip(j)%r)
!             v1=lo_operate_on_vector(p%sym%op(-k),v1,reciprocal=.true.,fractional=.true.,inverse=.false.)
!             v2=lo_clean_fractional_coordinates(v0+v1+0.5_r8)-0.5_r8
!             if ( lo_sqnorm(v2) .gt. lo_sqrectol ) then
!                 write(*,*) 'ERROR, bad symmetry operation in q-point mesh, fractional coordinates'
!                 err=.true.
!             endif
!         endif
!     enddo
!
!     ! Test the outward operations as well.
!
!     ! Test weights
!     f0=0.0_r8
!     do i=1,qp%n_irr_point
!         f0=f0+qp%ip(i)%integration_weight
!     enddo
!     if ( abs(f0-1.0_r8) .gt. lo_sqtol ) then
!         write(*,*) 'ERROR, irreducible weights add up to',f0
!         err=.true.
!     endif
!     f0=0.0_r8
!     do i=1,qp%n_full_point
!         f0=f0+qp%ap(i)%integration_weight
!     enddo
!     if ( abs(f0-1.0_r8) .gt. lo_sqtol ) then
!         write(*,*) 'ERROR, total weights add up to',f0
!         err=.true.
!     endif
!
!     ! Make sure no index is out of bounds
!     do i=1,qp%n_irr_point
!         j=qp%ip(i)%full_index
!         if ( j .lt. 1 .or. j .gt. qp%n_full_point ) then
!             write(*,*) 'ERROR, bad full grid index',j
!             err=.true.
!         endif
!     enddo
!     do i=1,qp%n_full_point
!         j=qp%ap(i)%irreducible_index
!         k=qp%ap(i)%operation_from_irreducible
!         if ( k .lt. 0 ) then
!         if ( qp%timereversal .eqv. .false. ) then
!             write(*,*) 'ERROR, negative operation but no timereversal'
!             err=.true.
!         endif
!         endif
!         if ( abs(k) .lt. 1 .or. abs(k) .gt. p%sym%n ) then
!             write(*,*) 'ERROR, index to operation that does not exist'
!             err=.true.
!         endif
!         if ( j .lt. 1 .or. j .gt. qp%n_irr_point ) then
!             write(*,*) 'ERROR, bad irreducible index',j
!             err=.true.
!         endif
!     enddo
!
!     do i=1,qp%n_irr_tet
!     do j=1,4
!         k=qp%it(i)%irreducible_index(j)
!         l=qp%it(i)%full_index(j)
!         if ( k .lt. 0 .or. k .gt. qp%n_irr_point ) then
!             write(*,*) 'ERROR, bad irreducible index in tetrahedron'
!             err=.true.
!         endif
!         if ( l .lt. 0 .or. l .gt. qp%n_full_point ) then
!             write(*,*) 'ERROR, bad gridindex in tetrahedron'
!             err=.true.
!         endif
!     enddo
!     enddo
!     do i=1,qp%n_full_tet
!     do j=1,4
!         k=qp%at(i)%irreducible_index(j)
!         l=qp%at(i)%full_index(j)
!         if ( k .lt. 0 .or. k .gt. qp%n_irr_point ) then
!             write(*,*) 'ERROR, bad irreducible index in tetrahedron'
!             err=.true.
!         endif
!         if ( l .lt. 0 .or. l .gt. qp%n_full_point ) then
!             write(*,*) 'ERROR, bad gridindex in tetrahedron'
!             err=.true.
!         endif
!     enddo
!     enddo
!
!
!     if ( err ) then
!         call lo_stop_gracefully(['Trouble with q-point mesh'],lo_exitcode_symmetry,__FILE__,__LINE__)
!     endif
! end subroutine

end module
