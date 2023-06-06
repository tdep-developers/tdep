#include "precompilerdefinitions"
module type_crystalstructure
!! Some text about the crystalstructure module
use konstanter, only: flyt,i8,lo_huge,lo_tiny,lo_pi,lo_hugeint,lo_tol,lo_sqtol,lo_radiantol,lo_status,lo_bohr_to_A,&
                      lo_A_to_bohr,lo_velocity_Afs_to_au,lo_amu_to_emu,lo_velocity_au_to_Afs,lo_exitcode_param,&
                      lo_exitcode_symmetry,lo_exitcode_io
use gottochblandat, only: open_file,tochar,qsort,walltime,lo_clean_fractional_coordinates,lo_chop,lo_determ,lo_choplarge,&
                   lo_frobnorm,lo_sqnorm,lo_reciprocal_basis,lo_mean,lo_unsigned_tetrahedron_volume,&
                   lo_index_in_periodic_array,lo_cross,lo_get_axis_angles,lo_invert3x3matrix,lo_stop_gracefully,&
                   lo_permutations,lo_return_unique
use geometryfunctions, only: lo_plane,lo_inscribed_sphere_in_box,lo_bounding_sphere_of_box
use type_distancetable, only: lo_distancetable
use type_voronoi, only: lo_voronoi_diagram,lo_voronoi_diagram_cell
use type_symmetryoperation, only: lo_symset, lo_operate_on_vector
use hdf5

implicit none
private
public :: lo_crystalstructure

!> define an atom with partial occupancy of different species, as in a random alloy
type lo_alloyatom
    !> how many components
    integer :: n=-lo_hugeint
    !> what are the components
    character(len=40), dimension(:), allocatable :: atomic_symbol
    !> concentration
    real(flyt), dimension(:), allocatable :: concentration
    !> atomic number of component
    integer, dimension(:), allocatable :: atomic_number
    !> mass of the components
    real(flyt), dimension(:), allocatable :: mass
end type

!> some shorthand to keep track of how disordered magnetic states are to be created
type lo_maginfo
    !> keep track of which atoms have magnetic moments
    logical, dimension(:), allocatable :: atom_has_moment
    !> collinear moment on each atom
    real(flyt), dimension(:), allocatable :: collinear_moment
    !> noncollinear moment on each atom
    real(flyt), dimension(:,:), allocatable :: noncollinear_moment
    !> species list is a little bit different in case of magnetism
    integer, dimension(:), allocatable :: magspecies
    !> number of (magnetic) species
    integer :: nmagspecies=-lo_hugeint
    !> counter for the magnetic species
    integer, dimension(:), allocatable :: magspeciescounter
end type

!> Distribution of isotopes for an atom.
type lo_isotope_distribution
    !> how many isotopes
    integer :: n=-lo_hugeint
    !> concentration of isotopes
    real(flyt), dimension(:), allocatable :: conc
    !> mass of isotope
    real(flyt), dimension(:), allocatable :: mass
    !> average mass
    real(flyt) :: mean_mass=lo_huge
    !> disorderparameter
    real(flyt) :: disorderparameter=lo_huge
    contains
        !> fetch the natural distribution
        procedure :: naturaldistribution
        !> calculate the mass disorder parameter
        procedure :: mass_disorder_parameter
end type

!> Brilluoin zone, stored as a polyhedron.
type, extends(lo_voronoi_diagram_cell) :: lo_brillouin_zone
    !> how many high symmetry points
    integer :: nhighsymmetrypoints=-lo_hugeint
    !> high-symmetry points, all of them
    real(flyt), dimension(:,:), allocatable :: highsymmetrypoints
    !> labels for the high-symmetry points
    character(len=10), dimension(:), allocatable :: label
    contains
        !> reciprocal lattice vector that shifts a point to the first bz
        procedure :: gshift
        !> calculates the distance to the zone edge.
        procedure :: distance_to_zone_edge
end type

!> A face in the wedge in the Brilluoin zone.
type lo_brillouin_zone_wedge_face
    !> how many points on this face?
    integer :: n=-lo_hugeint
    !> which nodes make up this face
    integer, dimension(:), allocatable :: ind
    !> what plane defines this face
    type(lo_plane) :: plane
end type

!> Irreducible wedge in the Brilluoin zone, polyhedron in reciprocal space.
type :: lo_brillouin_zone_irreducible_wedge
    !> how many irreducible points
    integer :: nnodes=-lo_hugeint
    !> how many faces in the polyhedron
    integer :: nfaces=-lo_hugeint
    !> faces
    type(lo_brillouin_zone_wedge_face), dimension(:), allocatable :: face
    !> coordinates to the high symmetry points
    real(flyt), dimension(:,:), allocatable :: r
    !> the labels of these points
    character(len=10), dimension(:), allocatable :: label
end type

!> Subcontainer with classification, none of which is particularly important.
type lo_crystalstructure_classification
    !> Just a title for this crystalstructure
    character(len=1000) :: title='empty header'
    ! Heuristics

    !> Is this an alloy?
    logical :: alloy=.false.
    !> have I specified collinear magnetic moments?
    logical :: collmag=.false.
    !> have I specified noncollinear magnetic moments?
    logical :: noncollmag=.false.
    !> have I calculated the symmetry operations?
    logical :: havespacegroup=.false.
    !> have I calculated the Brillouin zone
    logical :: havebz=.false.
    !> have I calculated the irreducible wedge
    logical :: havewedge=.false.
    !> have I figured out which Bravais lattice it is?
    logical :: havebravais=.false.
    !> Are all high symmetry points labelled ok?
    logical :: pointslabelled=.false.
    !> Have I deced on time-reversal symmetry?
    logical :: decidedtimereversal=.false.
    !> Have I calculated the character table?
    logical :: havecharactertable=.false.

    ! Classification things

    !> unique representation of the primitive lattice
    real(flyt), dimension(3,3) :: unique_primitive_basis=lo_huge
    !> unique representation of the conventional lattice
    real(flyt), dimension(3,3) :: unique_conventional_basis=lo_huge
    !> Permutation matrix from the basis to the unique primitive basis
    real(flyt), dimension(3,3) :: permutation_to_unique=lo_huge
    !> Transformation from the basis to the unique primitive basis
    real(flyt), dimension(3,3) :: transformation_to_unique=lo_huge
    !> the "a" lattice parameter, for pretty output.
    real(flyt) :: unitcell_lattice_parameter=lo_huge
    !> What Bravais lattice is it?
    character(len=10) :: bravaislattice='nothing'
    character(len=50) :: bravaislatticelongname='nothing'

    ! Supercell stuff

    !> Is this a supercell?
    logical :: supercell=.false.
    !> What are the dimensions?
    integer, dimension(3,3) :: supercellmatrix=-lo_hugeint
    !> For Fourier transforms and things like that it's good to know the index to the unitcell in the supercell
    integer, dimension(:,:), allocatable :: cellindex
    !> Might also be nice to have an index in the unit cell that it can be related to
    integer, dimension(:), allocatable :: index_in_unitcell
    !> How much should it talk?
    integer :: verbosity=-lo_hugeint
end type

!> Crystal structure class. Contains where the atoms are, how many and so on, basically everything about the crystal structure. When initialized by reading a structure from file, a whole bunch of stuff is calculated and gathered from tables, such as atomic masses, atomic numbers, space group, Brillouin zone and so on.
type lo_crystalstructure
    !> Isotope distributions
    type(lo_isotope_distribution), dimension(:), allocatable :: isotope
    !> Info about a species defined as a random alloy
    type(lo_alloyatom), dimension(:), allocatable :: alloyspecies
    !> Brillouin zone
    type(lo_brillouin_zone) :: bz
    !> irreducible wedge of the Brillouin zone
    type(lo_brillouin_zone_irreducible_wedge) :: irrw
    !> The space group for this structure
    type(lo_symset) :: sym
    !> Some extra classification information, if that is useful
    type(lo_crystalstructure_classification) :: info
    !> Information about the magnetic state
    type(lo_maginfo) :: mag
    !> number of atoms in the cell
    integer :: na=-lo_hugeint
    !> basis vectors
    real(flyt), dimension(3,3) :: latticevectors=lo_huge
    !> inverse of basis vectors
    real(flyt), dimension(3,3) :: inv_latticevectors=lo_huge
    !> reciprocal lattice vectors
    real(flyt), dimension(3,3) :: reciprocal_latticevectors=lo_huge
    !> inverse reciprocal lattice vectors
    real(flyt), dimension(3,3) :: inv_reciprocal_latticevectors=lo_huge
    !> Volume of the cell
    real(flyt) :: volume=lo_huge
    !> How many different elements
    integer :: nelements=-lo_hugeint
    !> Counter for each element type
    integer, dimension(:), allocatable :: element_counter
    !> What are the atomic symbols of these elements, e.g. Li,Fe,O
    character(len=8), dimension(:), allocatable :: atomic_symbol
    !> Atomic number for each atom
    integer, dimension(:), allocatable :: atomic_number
    !> What species is atom i?
    integer, allocatable, dimension(:) :: species
    !> Positions, fractional coordinates
    real(flyt), allocatable, dimension(:,:) :: r
    !> Positions, cartesian coordinates
    real(flyt), allocatable, dimension(:,:) :: rcart
    !> Velocities
    real(flyt), allocatable, dimension(:,:) :: v
    !> Forces
    real(flyt), allocatable, dimension(:,:) :: f
    !> Displacements
    real(flyt), allocatable, dimension(:,:) :: u
    !> Mass of each atom
    real(flyt), allocatable, dimension(:) :: mass
    !> Inverse square root of mass of each atom
    real(flyt), allocatable, dimension(:) :: invsqrtmass
    !> The inelastic neutron cross-section
    real(flyt), allocatable, dimension(:) :: inelastic_neutron_cross_section
    !> Flavor on each atom, that is to be preserved
    integer, allocatable, dimension(:) :: flavor
    contains
        !> create the structure
        procedure :: generate
        !> classify it in different ways
        procedure :: classify
        !> coordinate conversion
        procedure :: fractional_to_cartesian
        procedure :: cartesian_to_fractional
        procedure :: displacement_fractional_to_cartesian
        !> get the kinetic energy from velocities
        procedure :: kinetic_energy
        !> build a supercell
        procedure :: build_supercell
        !> largest pair cutoff in cell
        procedure :: maxcutoff
        !> smallest pair cutoff in cell
        procedure :: mincutoff
        !> nearest neighbour distance
        procedure :: nearest_neighbour_distance
        !> convert a coordinate to a high symmetry point
        procedure :: coordinate_from_high_symmetry_point_label
        !> return the permutation of an SQS
        procedure :: alloy_site_permutation
        !> permute the positions in a structure
        procedure :: permute_positions
        !> return the names of the unique atoms
        procedure :: unique_atom_label

        !> Initialize from file
        procedure :: readfromfile
        !> write to file, various formats
        procedure :: writetofile
        !> read isotope distribution from file
        procedure :: readisotopefromfile
        !> write the Brillouin zone to hdf5
        procedure :: write_bz_to_hdf5
        !> write structure information to hdf5
        procedure :: write_structure_to_hdf5

        !> measure size in memory
        procedure :: size_in_mem=>structure_size_in_mem
end type

! Global max for the max number of components in an alloy
integer, parameter :: max_n_components=6

! Interfaces to things in type_crystalstructure_atomdata
interface
    module function z_to_symbol(z_nucleus) result(symbol)
        integer, intent(in) :: z_nucleus
        character(len=2) :: symbol
    end function
    module function symbol_to_z(symbol) result(z_nucleus)
        character(len=*), intent(in) :: symbol
        integer :: z_nucleus
    end function
    module function neutron_cross_section(atomic_number) result(xs)
        integer, intent(in) :: atomic_number
        real(flyt) :: xs
    end function
end interface
! Interfaces to things in type_crystalstructure_io
interface
    module subroutine readfromfile(p,filename,verbosity)
        class(lo_crystalstructure), intent(out) :: p
        character(len=*), intent(in) :: filename
        integer, intent(in), optional :: verbosity
    end subroutine
    module subroutine writetofile(p,filename,output_format,write_velocities,transformationmatrix)
        class(lo_crystalstructure), intent(in) :: p
        character(len=*), intent(in) :: filename
        integer, intent(in) :: output_format
        logical, intent(in), optional :: write_velocities
        real(flyt), dimension(3,3), intent(out), optional :: transformationmatrix
    end subroutine
    module subroutine readisotopefromfile(p)
        class(lo_crystalstructure), intent(inout) :: p
    end subroutine
    module subroutine write_bz_to_hdf5(p,filename,input_id)
        class(lo_crystalstructure), intent(in) :: p
        character(len=*), intent(in) :: filename
        integer(HID_T), intent(in), optional :: input_id
    end subroutine
    module subroutine write_structure_to_hdf5(p,filename,input_id)
        class(lo_crystalstructure), intent(in) :: p
        character(len=*), intent(in), optional :: filename
        integer(HID_T), intent(in), optional :: input_id
    end subroutine
end interface
! Interfaces to things in type_crystalstructure_symmetry
interface
    module function coordinate_from_high_symmetry_point_label(p,label,previous) result(qpoint)
        class(lo_crystalstructure), intent(in) :: p
        character(len=*), intent(in) :: label
        real(flyt), dimension(3), intent(in), optional :: previous
        real(flyt), dimension(3) :: qpoint
    end function
    module subroutine classify(p,how,uc,tolerance,refine,timereversal)
        class(lo_crystalstructure), intent(inout) :: p
        character(len=*), intent(in) :: how
        type(lo_crystalstructure), intent(in), optional :: uc
        real(flyt), intent(in), optional :: tolerance
        logical, intent(in), optional :: refine
        logical, intent(in), optional :: timereversal
    end subroutine
end interface
! Interfaces to utility alloy routines
interface
    module subroutine alloy_site_permutation(ss,sqs,permutation)
        class(lo_crystalstructure), intent(in) :: ss
        type(lo_crystalstructure), intent(in) :: sqs
        integer, dimension(:), intent(out) :: permutation
    end subroutine
    module subroutine permute_positions(p,permutation,forward)
        class(lo_crystalstructure), intent(inout) :: p
        integer, dimension(:), intent(in) :: permutation
        logical, intent(in) :: forward
    end subroutine
end interface

contains

!> Transform lists of stuff to a proper object. How a crystal structure should be initialized.
subroutine generate(p,latticevectors,positions,atomic_numbers,enhet,velocities,header,&
                    collmag,noncollmag,cmatom,collmagmom,noncollmagmom,verbosity,alloy,&
                    alloy_componentcounter,alloy_components,alloy_concentrations)
    !> crystal structure
    class(lo_crystalstructure), intent(out) :: p
    !> basis
    real(flyt), dimension(3,3), intent(in) :: latticevectors
    !> positions of atoms, in fractional coordinates
    real(flyt), dimension(:,:), intent(in) :: positions
    !> the atomic numbers of each atom
    integer, dimension(:), intent(in) :: atomic_numbers
    !> in what unit is the lattice vectors?
    integer, intent(in) :: enhet
    !> specify velocities?
    real(flyt), dimension(:,:), intent(in), optional :: velocities
    !> give it a name?
    character(len=*), intent(in), optional :: header
    !> collinear magnetism?
    logical, intent(in), optional :: collmag
    !> non-collinear magnetism?
    logical, intent(in), optional :: noncollmag
    !> if so, which atoms are to be randomized?
    logical, dimension(:), intent(in), optional :: cmatom
    !> some magnetic moments
    real(flyt), dimension(:), intent(in), optional :: collmagmom
    !> and some non-collinear magnetic moments
    real(flyt), dimension(:,:), intent(in), optional :: noncollmagmom
    !> is it an alloy
    logical, intent(in), optional :: alloy
    !> number of alloy components per site
    integer, dimension(:), intent(in), optional :: alloy_componentcounter
    !> what species are there per site
    integer, dimension(:,:), intent(in), optional :: alloy_components
    !> concentration of these species
    real(flyt), dimension(:,:), intent(in), optional :: alloy_concentrations
    !> talk a lot?
    integer, intent(in), optional :: verbosity

    ! Since it is intent(out), nothing is allocated. Set everything to nothing, to be on the
    ! safe side, or something if it is obvious. Also, have some heuristics to figure out if the
    ! input is stupid.
    if ( size(positions,2) .eq. size(atomic_numbers,1) ) then
        p%na=size(positions,2)
    else
        call lo_stop_gracefully(['You need the same number of positions as atomic numbers when creating a structure.'],&
                                lo_exitcode_param,__FILE__,__LINE__)
    endif

    ! Right now, I have done nothing, and know nothing.
    if ( present(verbosity) ) then
        p%info%verbosity=verbosity
    else
        p%info%verbosity=0
    endif
    if ( present(header) ) then
        p%info%title=trim(adjustl(header))
    else
        p%info%title='cell' ! This has to be something.
    endif
    if ( present(alloy) ) then
        p%info%alloy=alloy
    else
        p%info%alloy=.false.
    endif
    if ( present(collmag) ) then
        p%info%collmag=collmag
    else
        p%info%collmag=.false.
    endif
    if ( present(noncollmag) ) then
        p%info%noncollmag=noncollmag
    else
        p%info%noncollmag=.false.
    endif
    p%info%havespacegroup=.false.
    p%info%havebz=.false.
    p%info%havewedge=.false.
    p%info%havebravais=.false.
    p%info%supercell=.false.
    p%info%pointslabelled=.false.
    p%info%decidedtimereversal=.false.
    p%info%supercellmatrix=-1
    p%info%unique_primitive_basis=-1
    p%info%unique_conventional_basis=-1
    p%info%permutation_to_unique=-1
    p%info%transformation_to_unique=-1
    p%info%bravaislattice='dunno'
    p%info%bravaislatticelongname='dunno'
    p%info%unitcell_lattice_parameter=-1.0_flyt

    ! Get the lattice vectors every which way, with some safety checks, as well as setting the positions
    ! the volume and stuff like that.
    fixlattice: block
        real(flyt) :: posunitfactor,velunitfactor
        integer :: i
        ! Make sure we have nonzero volume
        if ( abs(lo_determ(latticevectors))/lo_frobnorm(latticevectors) .lt. lo_tol ) then
            call lo_stop_gracefully(['The determinant of the basis is really small: '//tochar(lo_determ(latticevectors))],&
                                    lo_exitcode_param,__FILE__,__LINE__)
        endif
        ! What unit are they in?
        select case(enhet)
        case(1)
            ! This is vasp-ish units, positions in A and velocitied in A/fs
            ! convert that to Hartree atomic units
            posunitfactor=lo_A_to_bohr
            velunitfactor=lo_velocity_Afs_to_au
        case(2)
            ! Input is in atomic units. Distances in Bohr, velocities in something weird.
            posunitfactor=1.0_flyt
            velunitfactor=1.0_flyt
        case default
            call lo_stop_gracefully(['Unknown unit for positions and velocities'],&
                                    lo_exitcode_param,__FILE__,__LINE__)
        end select
        ! Store them
        p%latticevectors=latticevectors*posunitfactor
        p%inv_latticevectors=lo_invert3x3matrix( p%latticevectors )
        p%reciprocal_latticevectors=lo_reciprocal_basis(p%latticevectors)
        p%inv_reciprocal_latticevectors=lo_invert3x3matrix( p%reciprocal_latticevectors )
        ! Clean a little, for cosmetic reasons
        p%latticevectors=lo_chop(p%latticevectors,lo_sqtol)
        p%inv_latticevectors=lo_chop(p%inv_latticevectors,lo_sqtol)
        p%reciprocal_latticevectors=lo_chop(p%reciprocal_latticevectors,lo_sqtol)
        p%inv_reciprocal_latticevectors=lo_chop(p%inv_reciprocal_latticevectors,lo_sqtol)
        ! Get the volume
        p%volume=abs(lo_determ(p%latticevectors))
        ! store the positions. Do the cleaning thing just to be sure.
        lo_allocate(p%r(3,p%na))
        lo_allocate(p%rcart(3,p%na))
        lo_allocate(p%u(3,p%na))
        lo_allocate(p%v(3,p%na))
        lo_allocate(p%f(3,p%na))
        p%r=0.0_flyt
        p%rcart=0.0_flyt
        p%u=0.0_flyt
        p%v=0.0_flyt
        p%f=0.0_flyt
        do i=1,p%na
            p%r(:,i)=lo_chop( lo_clean_fractional_coordinates(positions(:,i)) , lo_sqtol )
            p%rcart(:,i)=lo_chop( p%fractional_to_cartesian(p%r(:,i)) , lo_sqtol )
        enddo
        ! also, velocities and stuff
        if ( present(velocities) ) then
            p%v=velocities*velunitfactor
        endif
    end block fixlattice

    if ( p%info%alloy ) then
    ! If I have an alloy this gets rather tedious.
    fixalloyatoms: block
        type(lo_isotope_distribution), dimension(:), allocatable :: dumiso
        integer :: i,j,k,l,nc

        ! Make some space for the things that need to be specified
        lo_allocate(p%species(p%na))
        lo_allocate(p%isotope(p%na))
        lo_allocate(p%mass(p%na))
        lo_allocate(p%invsqrtmass(p%na))
        lo_allocate(p%atomic_number(p%na))
        lo_allocate(p%inelastic_neutron_cross_section(p%na))

        ! First step is to figure out how many distinct species I have
        p%species=0
        p%nelements=0
        do i=1,p%na
            if ( p%species(i) .ne. 0 ) cycle
            ! new species!
            p%nelements=p%nelements+1
            p%species(i)=p%nelements
            do j=i+1,p%na
                if ( alloy_componentcounter(i) .eq. alloy_componentcounter(j) ) then
                if ( sum(abs(alloy_components(:,i)-alloy_components(:,j))) .eq. 0 ) then
                if ( sum(abs(alloy_concentrations(:,i)-alloy_concentrations(:,j))) .lt. lo_tol ) then
                    p%species(j)=p%species(i)
                endif
                endif
                endif
            enddo
        enddo
        ! Get a symbol for each species
        lo_allocate(p%atomic_symbol(p%nelements))
        do i=1,p%nelements
            p%atomic_symbol(i)='dunno'
        enddo

        ! Not sure if atomic number make sense here
        do i=1,p%na
            p%atomic_number(i)=-1
        enddo

        ! Count the number of components per species
        lo_allocate(p%element_counter(p%nelements))
        p%element_counter=0
        do i=1,p%na
            p%element_counter(p%species(i))=p%element_counter(p%species(i))+1
        enddo
        ! Get the masses. This is far too careful to be sane, but whatever.
        do i=1,p%na
            nc=alloy_componentcounter(i)
            lo_allocate(dumiso(nc))
            ! Get the natural isotope distribution for all components of the alloy
            do j=1,nc
                call dumiso(j)%naturaldistribution( trim(z_to_symbol(alloy_components(j,i))) )
            enddo
            ! Count total number of isotopes
            l=0
            do j=1,nc
                l=l+dumiso(j)%n
            enddo
            ! Make space
            p%isotope(i)%n=l
            lo_allocate(p%isotope(i)%conc( p%isotope(i)%n ))
            lo_allocate(p%isotope(i)%mass( p%isotope(i)%n ))
            ! Store concentrations and masses
            l=0
            do j=1,nc
            do k=1,dumiso(j)%n
                l=l+1
                p%isotope(i)%conc(l)=dumiso(j)%conc(k)*alloy_concentrations(j,i)
                p%isotope(i)%mass(l)=dumiso(j)%mass(k)
            enddo
            enddo
            ! ensure it adds up to 1
            p%isotope(i)%conc=p%isotope(i)%conc/sum(p%isotope(i)%conc)
            p%isotope(i)%mean_mass=sum(p%isotope(i)%conc*p%isotope(i)%mass)
            p%isotope(i)%disorderparameter=p%isotope(i)%mass_disorder_parameter()
            lo_deallocate(dumiso)
            ! Store mean masses
            p%mass(i)=p%isotope(i)%mean_mass
            p%invsqrtmass(i)=1.0_flyt/sqrt(p%isotope(i)%mean_mass)
            p%inelastic_neutron_cross_section(i)=0.0_flyt !neutron_cross_section(p%atomic_number(i))
        enddo

        ! Perhaps store all the alloy information somehow.
        lo_allocate(p%alloyspecies(p%nelements))
        lo_allocate(dumiso(1))
        do i=1,p%nelements
            ! locate an atom of this species
            k=-1
            do j=1,p%na
                if ( p%species(j) .eq. i ) then
                    k=j
                    exit
                endif
            enddo
            p%alloyspecies(i)%n=alloy_componentcounter(k)
            lo_allocate(p%alloyspecies(i)%mass( p%alloyspecies(i)%n ))
            lo_allocate(p%alloyspecies(i)%concentration( p%alloyspecies(i)%n ))
            lo_allocate(p%alloyspecies(i)%atomic_symbol( p%alloyspecies(i)%n ))
            lo_allocate(p%alloyspecies(i)%atomic_number( p%alloyspecies(i)%n ))
            do j=1,p%alloyspecies(i)%n
                p%alloyspecies(i)%atomic_number(j)=alloy_components(j,k)
                p%alloyspecies(i)%atomic_symbol(j)=trim(z_to_symbol( alloy_components(j,k) ))
                p%alloyspecies(i)%concentration(j)=alloy_concentrations(j,k)
                call dumiso(1)%naturaldistribution( trim(z_to_symbol(alloy_components(j,k))) )
                p%alloyspecies(i)%mass(j)=dumiso(1)%mean_mass
            enddo
        enddo
        lo_deallocate(dumiso)

        ! There might be non-alloy sublattices, fix the species and atomic number of those
        do i=1,p%nelements
            if ( p%alloyspecies(i)%n .ne. 1 ) cycle
            j=p%alloyspecies(i)%atomic_number(1)
            do k=1,p%na
                if ( p%species(k) .eq. i ) p%atomic_number(k)=j
            enddo
            p%atomic_symbol(i)=trim(z_to_symbol( j ))
        enddo

    end block fixalloyatoms
    else
    ! Figure out how classify the species in a neat way. Once I fix alloys, I have to do something else.
    ! This gets the symbols, species in the right way, masses and so on.
    fixatoms: block
        integer, dimension(:), allocatable :: dum
        integer :: i,j
        character(len=1000) :: trams
        ! start by getting the union of the atomic numbers
        call lo_return_unique(atomic_numbers,dum)
        ! number of elements in the system
        p%nelements=size(dum,1)
        lo_allocate(p%isotope(p%na))
        lo_allocate(p%atomic_symbol(p%nelements))
        lo_allocate(p%element_counter(p%nelements))
        lo_allocate(p%species(p%na))
        lo_allocate(p%mass(p%na))
        lo_allocate(p%invsqrtmass(p%na))
        lo_allocate(p%atomic_number(p%na))
        lo_allocate(p%inelastic_neutron_cross_section(p%na))
        ! Set the species
        p%element_counter=0
        do i=1,p%na
            do j=1,p%nelements
                if ( dum(j) .eq. atomic_numbers(i) ) then
                    p%species(i)=j
                    p%element_counter(j)=p%element_counter(j)+1
                endif
            enddo
        enddo
        ! Get the symbols
        do i=1,p%nelements
            p%atomic_symbol(i)=trim(adjustl(z_to_symbol( dum(i) )))
        enddo
        ! And the natural isotope distributions, and masses, and atomic numbers. Also neutron cross-section.
        do i=1,p%na
            call p%isotope(i)%naturaldistribution( p%atomic_symbol(p%species(i)) )
            p%mass(i)=p%isotope(i)%mean_mass
            p%invsqrtmass(i)=1.0_flyt/sqrt(p%mass(i))
            p%atomic_number(i)=atomic_numbers(i)
            p%inelastic_neutron_cross_section(i)=neutron_cross_section(p%atomic_number(i))
        enddo
        ! And maybe talk a little
        if ( p%info%verbosity .gt. 0 ) then
            write(*,*) 'Generating structure, lattice vectors:'
            write(*,"(1X,A,3(2X,F18.12))") '    a1:',p%latticevectors(:,1)
            write(*,"(1X,A,3(2X,F18.12))") '    a2:',p%latticevectors(:,2)
            write(*,"(1X,A,3(2X,F18.12))") '    a3:',p%latticevectors(:,3)
            write(*,*) '... reciprocal lattice vectors:'
            write(*,"(1X,A,3(2X,F18.12))") '    b1:',p%reciprocal_latticevectors(:,1)
            write(*,"(1X,A,3(2X,F18.12))") '    b2:',p%reciprocal_latticevectors(:,2)
            write(*,"(1X,A,3(2X,F18.12))") '    b3:',p%reciprocal_latticevectors(:,3)
            trams=''
            do i=1,p%nelements
                trams=trim(trams)//' '//trim(p%atomic_symbol(i))//': '//tochar(p%element_counter(i))
            enddo
            write(*,*) '... with composition '//trim(trams)
        endif
    end block fixatoms
    ! And for now, I am satisfied. I don't classify unless I really have to, I suppose.
    ! maybe add an option to agressively classify later on.
    endif

    ! If I specify magnetic moments, that is useful to take care of.
    if ( p%info%collmag .or. p%info%noncollmag ) then
    fixmag: block
        real(flyt), dimension(:,:), allocatable :: unmom,dum
        integer, dimension(:), allocatable :: dumi1,dumi2
        integer :: i,j,ii
        logical :: coll
        ! Quick check of one or three magnetic moments
        if ( p%info%collmag ) then
            coll=.true.
        else
            coll=.false.
        endif
        ! Small sanity check
        if ( p%info%collmag .and. p%info%noncollmag ) then
            call lo_stop_gracefully(['Please specifify either collinear or noncollinear magnetism, not both.'],&
                                    lo_exitcode_param,__FILE__,__LINE__)
        endif
        ! Count the number of magnetic atoms, and get their moments
        lo_allocate(dum(3,p%na))
        dum=0.0_flyt
        ii=0
        do i=1,p%na
            if ( cmatom(i) ) then
                ii=ii+1
                if ( coll ) then
                    dum(1,ii)=collmagmom(i)
                else
                    dum(:,ii)=noncollmagmom(:,i)
                endif
            endif
        enddo
        ! Get the unique moments
        call lo_return_unique(dum(:,1:ii),unmom,lo_tol)
        ! Build a new kind of species
        lo_allocate(dumi1(p%na))
        do i=1,p%na
            if ( cmatom(i) ) then
                do j=1,size(unmom,2)
                    if ( lo_sqnorm(unmom(:,j)-dum(:,i)) .lt. lo_sqtol ) then
                        dumi1(i)=atomic_numbers(i)*1000+j
                    endif
                enddo
            else
                dumi1(i)=atomic_numbers(i)
            endif
        enddo
        ! The new, unique species
        call lo_return_unique(dumi1,dumi2)
        ! Start storing some things
        p%mag%nmagspecies=size(dumi2,1)
        lo_allocate(p%mag%magspeciescounter( p%mag%nmagspecies ))
        lo_allocate(p%mag%magspecies( p%na ))
        p%mag%magspeciescounter=0
        do i=1,p%na
            do j=1,p%mag%nmagspecies
                if ( dumi1(i) .eq. dumi2(j) ) then
                    p%mag%magspeciescounter(j)=p%mag%magspeciescounter(j)+1
                    p%mag%magspecies(i)=j
                    exit
                endif
            enddo
        enddo
        lo_allocate(p%mag%atom_has_moment(p%na))
        if ( coll ) then
            lo_allocate(p%mag%collinear_moment(p%na))
            p%mag%collinear_moment=collmagmom
        else
            lo_allocate(p%mag%noncollinear_moment(3,p%na))
            p%mag%noncollinear_moment=noncollmagmom
        endif
        p%mag%atom_has_moment=cmatom
        !
        lo_deallocate(unmom)
        lo_deallocate(dum)
        lo_deallocate(dumi1)
        lo_deallocate(dumi2)
    end block fixmag
    endif
end subroutine

!> Build a supercell. The cell will be returned in ss.
subroutine build_supercell(p,ss,dimensions,nondiagdimensions)
    !> unitcell
    class(lo_crystalstructure), intent(inout) :: p
    !> supercell
    type(lo_crystalstructure), intent(out) :: ss
    !> how many times to repeat the unit cell along a1,a2,a3, the basis vectors.
    integer, dimension(3), intent(in), optional :: dimensions
    !> non-diagonal supercell. Fancy.
    integer, dimension(3,3), intent(in), optional :: nondiagdimensions

    real(flyt), dimension(:,:), allocatable :: r,noncollmagmom,alloy_concentrations
    real(flyt), dimension(:), allocatable :: collmagmom
    real(flyt), dimension(3,3) :: basis
    integer, dimension(:,:), allocatable :: alloy_components
    integer, dimension(:), allocatable :: atomic_numbers,alloy_componentcounter
    logical, dimension(:), allocatable :: magatom
    logical :: diagonal
    character(len=1000) :: hdr

    ! Basic heuristics
    if ( present(dimensions) ) then
        diagonal=.true.
    elseif ( present(nondiagdimensions) ) then
        diagonal=.false.
    else
        call lo_stop_gracefully(['You need to specify dimensions when building a supercell somehow.'],lo_exitcode_param,__FILE__,__LINE__)
    endif

    ! I want to know what Bravais lattice it is
    if ( p%info%havebravais .eqv. .false. ) call p%classify('bravais')

    if ( diagonal ) then
        call diagonal_cell(p,dimensions,basis,atomic_numbers,r,magatom,collmagmom,noncollmagmom,&
                          alloy_componentcounter,alloy_components,alloy_concentrations)
        hdr=tochar(dimensions(1))//'x'//tochar(dimensions(2))//'x'//tochar(dimensions(3))//' supercell'
    else
        call nondiagonal_cell(p,nondiagdimensions,basis,atomic_numbers,r,magatom,collmagmom,noncollmagmom,&
                              alloy_componentcounter,alloy_components,alloy_concentrations)
        hdr='non-diagonal supercell'
    endif

    ! Create a nice object from this
    if ( p%info%collmag ) then
        call ss%generate(basis,r,atomic_numbers,enhet=2,header=trim(hdr),verbosity=p%info%verbosity,collmag=.true.,&
                         cmatom=magatom,collmagmom=collmagmom)
    elseif ( p%info%noncollmag ) then
        call ss%generate(basis,r,atomic_numbers,enhet=2,header=trim(hdr),verbosity=p%info%verbosity,&
                         noncollmag=.true.,cmatom=magatom,noncollmagmom=noncollmagmom)
    elseif ( p%info%alloy ) then
        call ss%generate(basis,r,atomic_numbers,enhet=2,header=trim(hdr),verbosity=p%info%verbosity,&
                        alloy=.true.,alloy_componentcounter=alloy_componentcounter,&
                        alloy_components=alloy_components,alloy_concentrations=alloy_concentrations)
    else
        call ss%generate(basis,r,atomic_numbers,enhet=2,header=trim(hdr),verbosity=p%info%verbosity)
    endif

    ! Store the latticeparameter, perhaps
    ss%info%unitcell_lattice_parameter=p%info%unitcell_lattice_parameter
    contains

    !> build normal diagonal supercell
    subroutine diagonal_cell(p,dimensions,basis,atomic_numbers,r,magatom,collmagmom,noncollmagmom,alloy_componentcounter,alloy_components,alloy_concentrations)
        !> unit cell
        type(lo_crystalstructure), intent(in) :: p
        !> supercell dimensions
        integer, dimension(3), intent(in) :: dimensions
        !> supercell positions
        real(flyt), dimension(:,:), allocatable, intent(out) :: r
        !> supercell species
        integer, dimension(:), allocatable, intent(out) :: atomic_numbers
        !> supercell lattice vectors
        real(flyt), dimension(3,3), intent(out) :: basis
        !> which atoms are magnetic
        logical, dimension(:), allocatable, intent(out) :: magatom
        !> what their collinear magnetic moments
        real(flyt), dimension(:), allocatable, intent(out) :: collmagmom
        !> or the non-collinear magnetic moments
        real(flyt), dimension(:,:), allocatable, intent(out) :: noncollmagmom
        !> alloy component counter
        integer, dimension(:), allocatable, intent(out) :: alloy_componentcounter
        !> what are the components
        integer, dimension(:,:), allocatable, intent(out) :: alloy_components
        !> concentration of said components
        real(flyt), dimension(:,:), allocatable, intent(out) :: alloy_concentrations

        integer :: i,j,k,l,a1,na,ii,jj

        ! build the cell
        basis=p%latticevectors
        basis(:,1)=basis(:,1)*dimensions(1)
        basis(:,2)=basis(:,2)*dimensions(2)
        basis(:,3)=basis(:,3)*dimensions(3)
        na=product(dimensions)*p%na
        lo_allocate(r(3,na))
        lo_allocate(atomic_numbers(na))
        lo_allocate(magatom(na))
        lo_allocate(collmagmom(na))
        lo_allocate(noncollmagmom(3,na))
        lo_allocate(alloy_componentcounter(na))
        lo_allocate(alloy_components(max_n_components,na))
        lo_allocate(alloy_concentrations(max_n_components,na))
        r=0.0_flyt
        atomic_numbers=0
        magatom=.false.
        collmagmom=0.0_flyt
        noncollmagmom=0.0_flyt
        alloy_componentcounter=0
        alloy_components=0
        alloy_concentrations=0.0_flyt
        ! loop and build everything
        l=0
        do a1=1,p%na
            do i=1,dimensions(1)
            do j=1,dimensions(2)
            do k=1,dimensions(3)
                l=l+1
                r(:,l)=[i-1,j-1,k-1]*1.0_flyt+p%r(:,a1)
                do ii=1,3
                    r(ii,l)=r(ii,l)/(1.0_flyt*dimensions(ii))
                enddo
                atomic_numbers(l)=p%atomic_number(a1)
                if ( p%info%collmag ) then
                    magatom(l)=p%mag%atom_has_moment(a1)
                    collmagmom(l)=p%mag%collinear_moment(a1)
                elseif ( p%info%noncollmag ) then
                    magatom(l)=p%mag%atom_has_moment(a1)
                    noncollmagmom(:,l)=p%mag%noncollinear_moment(:,a1)
                elseif ( p%info%alloy ) then
                    ii=p%species(a1)
                    jj=p%alloyspecies(ii)%n
                    alloy_componentcounter(l)=p%alloyspecies(ii)%n
                    alloy_concentrations(1:jj,l)=p%alloyspecies(ii)%concentration(1:jj)
                    alloy_components(1:jj,l)=p%alloyspecies(ii)%atomic_number(1:jj)
                else
                    magatom(l)=.false.
                    collmagmom(l)=0.0_flyt
                    noncollmagmom(:,l)=0.0_flyt
                endif
            enddo
            enddo
            enddo
        enddo
        ! Clean the positions a little:
        r=lo_chop(lo_clean_fractional_coordinates(r),lo_sqtol)
        basis=lo_chop(basis,lo_sqtol)
    end subroutine

    !> build fancy non-diagonal supercell
    subroutine nondiagonal_cell(p,dimensions,basis,atomic_numbers,r,magatom,collmagmom,noncollmagmom,alloy_componentcounter,alloy_components,alloy_concentrations)
        !> unit cell
        type(lo_crystalstructure), intent(in) :: p
        !> supercell dimensions
        integer, dimension(3,3), intent(in) :: dimensions
        !> supercell lattice vectors
        real(flyt), dimension(3,3), intent(out) :: basis
        !> supercell species
        integer, dimension(:), allocatable, intent(out) :: atomic_numbers
        !> supercell positions
        real(flyt), dimension(:,:), allocatable, intent(out) :: r
        !> which atoms are magnetic
        logical, dimension(:), allocatable, intent(out) :: magatom
        !> what are their magnetic moments
        real(flyt), dimension(:), allocatable, intent(out) :: collmagmom
        !> or the non-collinear magnetic moments
        real(flyt), dimension(:,:), allocatable, intent(out) :: noncollmagmom
        !> alloy component counter
        integer, dimension(:), allocatable, intent(out) :: alloy_componentcounter
        !> what are the components
        integer, dimension(:,:), allocatable, intent(out) :: alloy_components
        !> concentration of said components
        real(flyt), dimension(:,:), allocatable, intent(out) :: alloy_concentrations

        real(flyt), dimension(:,:), allocatable :: dumr1,dumr2
        real(flyt), dimension(3,3) :: invbasis,m0,m1
        real(flyt), dimension(3) :: v0
        real(flyt) :: f0
        integer, dimension(:), allocatable :: ind
        integer :: i,j,k,l,a1,na,nrep,ctr,ii,jj

        ! Get the new basis
        i=lo_determ(dimensions)
        if ( i .eq. 0 ) then
            call lo_stop_gracefully(['Determinant of supercell matrix is zero'],lo_exitcode_param,__FILE__,__LINE__)
        elseif ( i .lt. 0 ) then
            write(*,*) 'NOTE: I flipped sign of the supercell matrix to not invert the supercell'
            basis=-lo_chop(matmul(p%latticevectors,dimensions),lo_sqtol)
        else
            basis=lo_chop(matmul(p%latticevectors,dimensions),lo_sqtol)
        endif
        ! Possibly, this could be cleaned in a smart way, to get good precision for
        ! the more obvious cases, such as cubic/tetragonal/hexagonal
        invbasis=lo_invert3x3matrix( basis )

        ! Figure out how many times the unitcell should be repeated
        f0=lo_bounding_sphere_of_box(basis)*2
        nrep=1
        do
            m0=p%latticevectors*(2*nrep+1)
            if ( lo_inscribed_sphere_in_box(m0) .gt. f0 ) then
                ! got it, was enough
                exit
            else
                ! check larger
                nrep=nrep+1
            endif
        enddo
        nrep=nrep+1 ! just to be on the safe side.
        ! Expected number of atoms
        na=abs(int(anint(lo_determ(dimensions*1.0_flyt)))*p%na)
        m1=lo_chop( matmul(invbasis,p%latticevectors) ,lo_sqtol)

        ! Count atoms
        ctr=0
        do a1=1,p%na
            do i=-nrep,nrep
            do j=-nrep,nrep
            do k=-nrep,nrep
                v0=matmul(m1,[i,j,k]+p%r(:,a1)) ! fractional coordinates in new system
                if ( minval(v0) .gt. -lo_sqtol .and. maxval(v0) .lt. 1.0_flyt+lo_sqtol ) ctr=ctr+1
            enddo
            enddo
            enddo
        enddo

        ! Store these
        lo_allocate(dumr1(6,ctr))
        ctr=0
        do a1=1,p%na
            do i=-nrep,nrep
            do j=-nrep,nrep
            do k=-nrep,nrep
                v0=matmul(m1,[i,j,k]+p%r(:,a1)) ! fractional coordinates in new system
                if ( minval(v0) .gt. -lo_sqtol .and. maxval(v0) .lt. 1.0_flyt+lo_sqtol ) then
                    ctr=ctr+1
                    dumr1(1:3,ctr)=lo_chop(lo_clean_fractional_coordinates(v0),lo_sqtol)
                    dumr1(4,ctr)=p%atomic_number(a1)
                    dumr1(5,ctr)=p%species(a1)
                    if ( p%info%collmag .or. p%info%noncollmag ) then
                        dumr1(6,ctr)=p%mag%magspecies(a1)
                    else
                        dumr1(6,ctr)=0
                    endif
                endif
            enddo
            enddo
            enddo
        enddo

        ! get the unique
        call lo_return_unique(dumr1,dumr2)

        ! Simple sanity check
        if ( size(dumr2,2) .ne. na ) then
            call lo_stop_gracefully(['Failed building supercell. Probably a pathological supercell matrix'],lo_exitcode_symmetry,__FILE__,__LINE__)
        endif

        ! A bit more serious sanity check: this is how many of each
        ! atom I should find.
        l=abs(lo_determ(dimensions))
        do i=1,p%nelements
            k=0
            do j=1,size(dumr2,2)
                if ( nint(dumr2(5,j)) .eq. i ) k=k+1
            enddo
            if ( k .ne. l*p%element_counter(i) ) then
                call lo_stop_gracefully(['Got the correct number of atoms, but the composition is off. Not good.'],lo_exitcode_symmetry,__FILE__,__LINE__)
            endif
        enddo

        ! Now I suppose things are ok. Prepare thing for output:
        basis=lo_choplarge(basis/p%info%unitcell_lattice_parameter,lo_sqtol)*p%info%unitcell_lattice_parameter
        lo_allocate(ind(na))
        lo_allocate(r(3,na))
        lo_allocate(atomic_numbers(na))
        lo_allocate(magatom(na))
        lo_allocate(collmagmom(na))
        lo_allocate(noncollmagmom(3,na))
        lo_allocate(alloy_componentcounter(na))
        lo_allocate(alloy_components(max_n_components,na))
        lo_allocate(alloy_concentrations(max_n_components,na))
        ind=0
        r=0.0_flyt
        atomic_numbers=0
        magatom=.false.
        collmagmom=0.0_flyt
        noncollmagmom=0.0_flyt
        alloy_componentcounter=0
        alloy_components=0
        alloy_concentrations=0.0_flyt
        call qsort(dumr2(5,:),ind)
        do i=1,na
            j=ind(i)
            r(:,i)=dumr2(1:3,j)
            atomic_numbers(i)=nint(dumr2(4,j))
            if ( p%info%collmag ) then
                do k=1,p%na
                    if ( nint(dumr2(6,j)) .eq. p%mag%magspecies(k) ) then
                        magatom(i)=p%mag%atom_has_moment(k)
                        collmagmom(i)=p%mag%collinear_moment(k)
                    endif
                enddo
            elseif ( p%info%noncollmag ) then
                do k=1,p%na
                    if ( nint(dumr2(6,j)) .eq. p%mag%magspecies(k) ) then
                        magatom(i)=p%mag%atom_has_moment(k)
                        noncollmagmom(:,i)=p%mag%noncollinear_moment(:,k)
                    endif
                enddo
            elseif ( p%info%alloy ) then
                ii=nint(dumr2(5,i)) ! note i here, not j
                jj=p%alloyspecies(ii)%n
                alloy_componentcounter(i)=p%alloyspecies(ii)%n
                alloy_concentrations(1:jj,i)=p%alloyspecies(ii)%concentration(1:jj)
                alloy_components(1:jj,i)=p%alloyspecies(ii)%atomic_number(1:jj)
            endif
        enddo
        ! and things are done!
        lo_deallocate(ind)
        lo_deallocate(dumr1)
        lo_deallocate(dumr2)
    end subroutine
end subroutine

!> Get the kinetic energy of the atoms in the cell. Will return it in Hartree/cell, not per atom.
real(flyt) function kinetic_energy(p)
    !> the crystal structure
    class(lo_crystalstructure), intent(in) :: p

    integer :: i
    real(flyt) :: mass,normv,ek

    ek=0.0_flyt
    do i=1,p%na
        mass=p%mass(i)
        normv=lo_sqnorm(p%v(:,i))
        ek=ek+mass*normv*0.5_flyt
    enddo
    kinetic_energy=ek
end function

!> The shortest distance in a cell, should be the shortest neighbour distance.
#ifdef AGRESSIVE_SANITY
function mincutoff(p) result(r)
#else
pure function mincutoff(p) result(r)
#endif
    !> the crystal structure
    class(lo_crystalstructure), intent(in) :: p
    !> the cutoff
    real(flyt) :: r

    integer :: i,j
    real(flyt), dimension(3) :: v
    real(flyt) :: f0,f1

    f0=-lo_huge
    do i=1,p%na
        f1=lo_huge
        do j=1,p%na
            v=p%r(:,j)-p%r(:,i)
            v=p%displacement_fractional_to_cartesian(v)
            if ( norm2(v) .gt. lo_tol ) then
                f1=min(f1,norm2(v))
            endif
        enddo
        f0=max(f0,f1)
    enddo
    ! Catch pathological case of feeding a unitcell to this routine
    if ( f0 .lt. lo_tol .or. p%na .eq. 1 ) then
        f0=lo_bounding_sphere_of_box(p%latticevectors)
    endif
    ! And a little bit more for good measure.
    r=f0*1.001_flyt
end function

!> calculate the nearest neighbour distance. Not fast.
#ifdef AGRESSIVE_SANITY
function nearest_neighbour_distance(p) result(r)
#else
pure function nearest_neighbour_distance(p) result(r)
#endif
    !> structure
    class(lo_crystalstructure), intent(in) :: p
    !> nearest neighbour distance
    real(flyt) :: r

    real(flyt), dimension(3,3) :: m0
    real(flyt), dimension(3) :: v0
    real(flyt) :: f0
    integer :: nrep,a1,a2,i,j,k

    f0=lo_bounding_sphere_of_box(p%latticevectors)
    do nrep=1,100
        m0=p%latticevectors*(2*nrep+1)
        if ( lo_inscribed_sphere_in_box(m0) .gt. f0 ) exit
    enddo

    r=lo_huge
    do a1=1,p%na
    do a2=1,p%na
        do i=-nrep,nrep
        do j=-nrep,nrep
        do k=-nrep,nrep
            v0=[i,j,k]+p%r(:,a2)-p%r(:,a1)
            f0=norm2(matmul(p%latticevectors,v0))
            if ( f0 .gt. lo_tol ) r=min(r,f0)
        enddo
        enddo
        enddo
    enddo
    enddo
end function

!> The longest distance in a cell. Will return the radius of the largest sphere that can fit safely inside the cell
#ifdef AGRESSIVE_SANITY
function maxcutoff(p) result(r)
#else
pure function maxcutoff(p) result(r)
#endif
    !> the crystal structure
    class(lo_crystalstructure), intent(in) :: p
    !> the cutoff
    real(flyt) :: r
    ! just get the radius of the inscribed sphere, with some tolerance
    r=lo_inscribed_sphere_in_box(p%latticevectors)-10*lo_tol
end function

!> Converts a vector from cartesian to fractional coordinates. Not fast at all, don't use for something speed-sensitive.
pure function cartesian_to_fractional(p,v,reciprocal,pbc) result(r)
    !> the crystal structure
    class(lo_crystalstructure), intent(in) :: p
    !> the vector in cartesian coordinates
    real(flyt), dimension(3), intent(in) :: v
    !> in reciprocal space?
    logical, intent(in), optional :: reciprocal
    !> should the pbc checks be done
    logical, intent(in), optional :: pbc
    !> the vector in fractional coordinates
    real(flyt), dimension(3) :: r
    !
    logical :: check,reclat

    ! Choice wether to care about periodic boundary conditions
    if ( present(pbc) ) then
        check=pbc
    else
        check=.true.
    endif
    ! Real or reciprocal space?
    if ( present(reciprocal) ) then
        reclat=reciprocal
    else
        reclat=.false.
    endif

    if ( reclat ) then
        r=matmul(p%inv_reciprocal_latticevectors,v)
    else
        r=matmul(p%inv_latticevectors,v)
    endif
    if ( check ) r=lo_clean_fractional_coordinates(r)
end function

!> Convert a displacement vector from fractional to cartesian coordinates. This means that a displacement in fractional coordinates larger than 0.5 will be shifted, as will those smaller than -0.5
pure function displacement_fractional_to_cartesian(p,v,reciprocal) result(r)
    !> the crystal structure
    class(lo_crystalstructure), intent(in) :: p
    !> displacement vector to be converted
    real(flyt), dimension(3), intent(in) :: v
    !> in reciprocal space?
    logical, intent(in), optional :: reciprocal
    !> cartesian displacement vector
    real(flyt), dimension(3) :: r

    ! local
    real(flyt), dimension(3) :: w
    logical :: reclat

    if ( present(reciprocal) ) then
        reclat=reciprocal
    else
        reclat=.false.
    endif

    w=lo_clean_fractional_coordinates(v+0.5_flyt)-0.5_flyt
    if ( reclat ) then
        r=matmul(p%reciprocal_latticevectors,w)
    else
        r=matmul(p%latticevectors,w)
    endif
end function

!> Converts a vector from fractional to cartesian coordinates.
pure function fractional_to_cartesian(p,v,reciprocal) result(r)
    !> the crystal structure
    class(lo_crystalstructure), intent(in) :: p
    !> vector in fractional coordinates
    real(flyt), dimension(3), intent(in) :: v
    !> in reciprocal space?
    logical, intent(in), optional :: reciprocal
    !> vector in cartesian coordinates
    real(flyt), dimension(3) :: r

    logical :: reclat
    if ( present(reciprocal) ) then
        reclat=reciprocal
    else
        reclat=.false.
    endif
    if ( reclat ) then
        r=matmul(p%reciprocal_latticevectors,v)
    else
        r=matmul(p%latticevectors,v)
    endif
end function

!> Gives the reciprocal lattice vector that moves q into the first BZ. Uses the information from the Brilloin zone polyhedron and moves around through edges in a pretty neat way. It is pretty fast.
pure function gshift(bz,point) result(g)
    !> brillouin zone
    class(lo_brillouin_zone), intent(in) :: bz
    !> point
    real(flyt), dimension(3), intent(in) :: point
    !> shift
    real(flyt), dimension(3) :: g

    integer :: i,facectr
    real(flyt), dimension(3) :: sh,pt,nv

    ! First a very cheap check, if the point is inside the inscribed sphere of the BZ, no need
    ! to do expensive checks, the shift is zero.
    if ( lo_sqnorm(point) .lt. bz%rmin**2 ) then
        g=0.0_flyt
        return
    endif

    pt=point
    sh=0.0_flyt
    shiftloop: do
        facectr=0
        faceloop: do i=1,bz%nfaces
            ! Check the signed distance between the point and all the planes constituting
            ! the Brillouin zone.
            if ( bz%face(i)%plane%distance_to_point(pt) .gt. lo_tol ) then
                ! it means the point is outside this face, get the vector pointing to the
                ! Gamma-point on the other side of this face
                nv=bz%face(i)%neighbourvector
                ! add this G-vector to the shift
                sh=sh+nv
                ! move the point with this vector
                pt=pt-nv
                exit faceloop
            else
                ! the point is on the correct side of this face
                facectr=facectr+1
            endif
        enddo faceloop
        ! We might be back in the first BZ if we are on the negative side of all faces.
        if ( facectr .eq. bz%nfaces ) exit shiftloop
    enddo shiftloop
    ! And we are done!
    g=sh
end function

!> Calculates distance to the Brillouin zone edge from a given wave vector. Will return a value between 0 and 1, where 1 is the zone center and 0 is exactly on the edge.
pure function distance_to_zone_edge(bz,q) result(d)
    !> the Brillouin zone
    class(lo_brillouin_zone), intent(in) :: bz
    !> the q-point
    real(flyt), dimension(3), intent(in) :: q
    !> the distance
    real(flyt) :: d
    !
    real(flyt), dimension(3) :: v0
    real(flyt) :: f0,f1,f2,qnorm
    integer :: i

    ! first check if we are at gamma, then we don't have to care at all.
    qnorm=norm2(q)
    if ( qnorm .lt. lo_tol ) then
        ! this is gamma
        d=1.0_flyt
        return
    endif

    ! unit vector for the line from Gamma to q.
    v0=q/qnorm
    ! Check the intersection between this line and all possible planes
    f0=lo_huge*0.5_flyt
    d=lo_huge
    do i=1,bz%nfaces
        f1=dot_product(v0,bz%face(i)%plane%normal)
        if ( abs(f1) .gt. lo_tol ) then
            ! not perpendicular, intersection at finite distance
            f2=bz%face(i)%plane%p/f1
        else
            ! intersection infinitely far away
            f2=lo_huge
        endif
        !
        if ( f2 .lt. f0 .and. f2 .gt. 0.0_flyt ) then
            f0=f2
            d=min(d,f2)
        endif
    enddo
    ! Straight ratio, 0 at zone center, 1 at edge
    d=lo_chop(1.0_flyt-qnorm/d,lo_sqtol)
end function

!> return a string containing the names of unique atoms separated by spaces
subroutine unique_atom_label(p,labelstring)
    !> structure
    class(lo_crystalstructure), intent(in) :: p
    !> labels
    character(len=2000), intent(out) :: labelstring

    character(len=100), dimension(:), allocatable :: unique_atom_names
    integer :: a1,a2,s1,s2,ctr

    ! Get names for the unique atoms?
    allocate(unique_atom_names(p%sym%n_irreducible_atom))
    unique_atom_names='dunno'
    a1l: do a1=1,p%sym%n_irreducible_atom
        ! If already decided, don't bother.
        if ( trim(unique_atom_names(a1)) .ne. 'dunno' ) cycle
        s1=p%species(p%sym%irr_to_all(a1))
        ctr=0
        do a2=1,p%sym%n_irreducible_atom
            s2=p%species(p%sym%irr_to_all(a2))
            if ( s1 .eq. s2 ) ctr=ctr+1
        enddo

        if ( ctr .eq. 1 ) then
            ! If only one, just give it the normal name
            unique_atom_names(a1)=trim(adjustl(p%atomic_symbol(s1)))
            ! and move on to the next atom
            cycle a1l
        endif

        ! More than one unique atom of this kind, decorate the names
        ctr=0
        do a2=1,p%sym%n_irreducible_atom
            s2=p%species(p%sym%irr_to_all(a2))
            if ( s1 .eq. s2 ) then
                ctr=ctr+1
                unique_atom_names(a2)=trim(adjustl(p%atomic_symbol(s2)))//"_"//tochar(ctr)
            endif
        enddo
    enddo a1l

    ! Stuff the names into a string
    labelstring=""
    do a1=1,p%sym%n_irreducible_atom
        labelstring=trim(adjustl(labelstring))//" "//trim(adjustl(unique_atom_names(a1)))
    enddo
    deallocate(unique_atom_names)
end subroutine

!> measure size in memory, in bytes
function structure_size_in_mem(p) result(mem)
    !> dispersions
    class(lo_crystalstructure), intent(in) :: p
    !> memory in bytes
    integer(i8) :: mem

    integer :: i

    mem=0
    mem=mem+storage_size(p)

    if ( allocated(p%flavor                         ) ) mem=mem+storage_size(p%flavor                         )*size(p%flavor                         )
    if ( allocated(p%inelastic_neutron_cross_section) ) mem=mem+storage_size(p%inelastic_neutron_cross_section)*size(p%inelastic_neutron_cross_section)
    if ( allocated(p%invsqrtmass                    ) ) mem=mem+storage_size(p%invsqrtmass                    )*size(p%invsqrtmass                    )
    if ( allocated(p%mass                           ) ) mem=mem+storage_size(p%mass                           )*size(p%mass                           )
    if ( allocated(p%u                              ) ) mem=mem+storage_size(p%u                              )*size(p%u                              )
    if ( allocated(p%f                              ) ) mem=mem+storage_size(p%f                              )*size(p%f                              )
    if ( allocated(p%v                              ) ) mem=mem+storage_size(p%v                              )*size(p%v                              )
    if ( allocated(p%rcart                          ) ) mem=mem+storage_size(p%rcart                          )*size(p%rcart                          )
    if ( allocated(p%r                              ) ) mem=mem+storage_size(p%r                              )*size(p%r                              )
    if ( allocated(p%species                        ) ) mem=mem+storage_size(p%species                        )*size(p%species                        )
    if ( allocated(p%atomic_number                  ) ) mem=mem+storage_size(p%atomic_number                  )*size(p%atomic_number                  )
    if ( allocated(p%atomic_symbol                  ) ) mem=mem+storage_size(p%atomic_symbol                  )*size(p%atomic_symbol                  )
    if ( allocated(p%element_counter                ) ) mem=mem+storage_size(p%element_counter                )*size(p%element_counter                )
    if ( allocated(p%info%cellindex                 ) ) mem=mem+storage_size(p%info%cellindex                 )*size(p%info%cellindex                 )
    if ( allocated(p%info%index_in_unitcell         ) ) mem=mem+storage_size(p%info%index_in_unitcell         )*size(p%info%index_in_unitcell         )
    if ( allocated(p%mag%atom_has_moment            ) ) mem=mem+storage_size(p%mag%atom_has_moment            )*size(p%mag%atom_has_moment            )
    if ( allocated(p%mag%collinear_moment           ) ) mem=mem+storage_size(p%mag%collinear_moment           )*size(p%mag%collinear_moment           )
    if ( allocated(p%mag%noncollinear_moment        ) ) mem=mem+storage_size(p%mag%noncollinear_moment        )*size(p%mag%noncollinear_moment        )
    if ( allocated(p%mag%magspecies                 ) ) mem=mem+storage_size(p%mag%magspecies                 )*size(p%mag%magspecies                 )
    if ( allocated(p%mag%magspeciescounter          ) ) mem=mem+storage_size(p%mag%magspeciescounter          )*size(p%mag%magspeciescounter          )

    if ( allocated(p%isotope) ) then
        do i=1,size(p%isotope)
            mem=mem+storage_size(p%isotope(i))
            if ( allocated(p%isotope(i)%conc ) ) mem=mem+storage_size(p%isotope(i)%conc)*size(p%isotope(i)%conc)
            if ( allocated(p%isotope(i)%mass ) ) mem=mem+storage_size(p%isotope(i)%mass)*size(p%isotope(i)%mass)
        enddo
    endif
    if ( allocated(p%alloyspecies) ) then
        do i=1,size(p%alloyspecies)
            mem=mem+storage_size(p%alloyspecies(i))
            if ( allocated(p%alloyspecies(i)%atomic_symbol) ) mem=mem+storage_size(p%alloyspecies(i)%atomic_symbol)*size(p%alloyspecies(i)%atomic_symbol)
            if ( allocated(p%alloyspecies(i)%concentration) ) mem=mem+storage_size(p%alloyspecies(i)%concentration)*size(p%alloyspecies(i)%concentration)
            if ( allocated(p%alloyspecies(i)%atomic_number) ) mem=mem+storage_size(p%alloyspecies(i)%atomic_number)*size(p%alloyspecies(i)%atomic_number)
            if ( allocated(p%alloyspecies(i)%mass         ) ) mem=mem+storage_size(p%alloyspecies(i)%mass         )*size(p%alloyspecies(i)%mass         )
        enddo
    endif
    ! take care of the bz
    if ( allocated(p%bz%highsymmetrypoints) ) mem=mem+storage_size(p%bz%highsymmetrypoints)*size(p%bz%highsymmetrypoints)
    if ( allocated(p%bz%label             ) ) mem=mem+storage_size(p%bz%label             )*size(p%bz%label             )
    if ( allocated(p%bz%r                 ) ) mem=mem+storage_size(p%bz%r                 )*size(p%bz%r                 )
    if ( allocated(p%bz%face) ) then
        do i=1,size(p%bz%face)
            mem=mem+storage_size(p%bz%face(i))
            if ( allocated(p%bz%face(i)%ind) ) mem=mem+storage_size(p%bz%face(i)%ind)*size(p%bz%face(i)%ind)
        enddo
    endif
    if ( allocated(p%bz%node) ) then
        do i=1,size(p%bz%node)
            mem=mem+storage_size(p%bz%node(i))
            if ( allocated(p%bz%node(i)%faceind        ) ) mem=mem+storage_size(p%bz%node(i)%faceind        )*size(p%bz%node(i)%faceind        )
            if ( allocated(p%bz%node(i)%neighbourind   ) ) mem=mem+storage_size(p%bz%node(i)%neighbourind   )*size(p%bz%node(i)%neighbourind   )
            if ( allocated(p%bz%node(i)%neighbourvector) ) mem=mem+storage_size(p%bz%node(i)%neighbourvector)*size(p%bz%node(i)%neighbourvector)
        enddo
    endif
    if ( allocated(p%bz%edge) ) then
        mem=mem+storage_size(p%bz%edge)*size(p%bz%edge)
    endif
    ! irreducible wedge?
    if ( allocated(p%irrw%label) ) mem=mem+storage_size(p%irrw%label)*size(p%irrw%label)
    if ( allocated(p%irrw%r    ) ) mem=mem+storage_size(p%irrw%r    )*size(p%irrw%r    )
    if ( allocated(p%irrw%face  ) ) then
        do i=1,size(p%irrw%face)
            mem=mem+storage_size(p%irrw%face(i))
            if ( allocated(p%irrw%face(i)%ind) ) mem=mem+storage_size(p%irrw%face(i)%ind)*size(p%irrw%face(i)%ind)
        enddo
    endif

    ! get it to bytes
    mem=mem/8
    ! Finally add the spacegroup
    mem=mem+p%sym%size_in_mem()
end function

!> Mass disorder parameter. This function takes a lo_isotope_distribution and calculates the mass disorder parameter.
pure function mass_disorder_parameter(self) result(g)
    !> the isotope distribution
    class(lo_isotope_distribution), intent(in) :: self
    !> the mass disorder parameter
    real(flyt) :: g
    !
    integer :: i
    real(flyt) :: f

    f=0.0_flyt
    do i=1,self%n
        f=f+self%conc(i)*((self%mass(i)-self%mean_mass)**2)
    enddo
    g=f/(self%mean_mass**2)
end function

!> Returns the natural isotope distribution for the specified element in atomic mass units:
subroutine naturaldistribution(self,symbol)
    !> the isotope distribution
    class(lo_isotope_distribution), intent(out) :: self
    !> Atomic symbol, e.g. "Si"
    character(len=*), intent(in) :: symbol

    ! Make sure they are not already allocated
    if ( allocated(self%conc) ) lo_deallocate(self%conc)
    if ( allocated(self%mass) ) lo_deallocate(self%mass)
    !
    select case(trim(symbol))
        case("H")
            self%n=2
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.9998850000_flyt
            self%conc(2)=0.0001150000_flyt
            self%mass(1)=    1.007825032070_flyt
            self%mass(2)=    2.014101777800_flyt
        case("He")
            self%n=2
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.0000013400_flyt
            self%conc(2)=0.9999986600_flyt
            self%mass(1)=    3.016029319100_flyt
            self%mass(2)=    4.002603254150_flyt
        case("Li")
            self%n=2
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.0759000000_flyt
            self%conc(2)=0.9241000000_flyt
            self%mass(1)=    6.015122795000_flyt
            self%mass(2)=    7.016004550000_flyt
        case("Be")
            self%n=1
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=1.0000000000_flyt
            self%mass(1)=    9.012182200000_flyt
        case("B")
            self%n=2
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.1990000000_flyt
            self%conc(2)=0.8010000000_flyt
            self%mass(1)=   10.012937000000_flyt
            self%mass(2)=   11.009305400000_flyt
        case("C")
            self%n=2
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.9893000000_flyt
            self%conc(2)=0.0107000000_flyt
            self%mass(1)=   12.000000000000_flyt
            self%mass(2)=   13.003354837800_flyt
        case("N")
            self%n=2
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.9963600000_flyt
            self%conc(2)=0.0036400000_flyt
            self%mass(1)=   14.003074004800_flyt
            self%mass(2)=   15.000108898200_flyt
        case("O")
            self%n=3
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.9975700000_flyt
            self%conc(2)=0.0003800000_flyt
            self%conc(3)=0.0020500000_flyt
            self%mass(1)=   15.994914619560_flyt
            self%mass(2)=   16.999131700000_flyt
            self%mass(3)=   17.999161000000_flyt
        case("F")
            self%n=1
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=1.0000000000_flyt
            self%mass(1)=   18.998403220000_flyt
        case("Ne")
            self%n=3
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.9048000000_flyt
            self%conc(2)=0.0027000000_flyt
            self%conc(3)=0.0925000000_flyt
            self%mass(1)=   19.992440175400_flyt
            self%mass(2)=   20.993846680000_flyt
            self%mass(3)=   21.991385114000_flyt
        case("Na")
            self%n=1
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=1.0000000000_flyt
            self%mass(1)=   22.989769280900_flyt
        case("Mg")
            self%n=3
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.7899000000_flyt
            self%conc(2)=0.1000000000_flyt
            self%conc(3)=0.1101000000_flyt
            self%mass(1)=   23.985041700000_flyt
            self%mass(2)=   24.985836920000_flyt
            self%mass(3)=   25.982592929000_flyt
        case("Al")
            self%n=1
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=1.0000000000_flyt
            self%mass(1)=   26.981538630000_flyt
        case("Si")
            self%n=3
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.9222300000_flyt
            self%conc(2)=0.0468500000_flyt
            self%conc(3)=0.0309200000_flyt
            self%mass(1)=   27.976926532500_flyt
            self%mass(2)=   28.976494700000_flyt
            self%mass(3)=   29.973770170000_flyt
        case("P")
            self%n=1
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=1.0000000000_flyt
            self%mass(1)=   30.973761630000_flyt
        case("S")
            self%n=4
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.9499000000_flyt
            self%conc(2)=0.0075000000_flyt
            self%conc(3)=0.0425000000_flyt
            self%conc(4)=0.0001000000_flyt
            self%mass(1)=   31.972071000000_flyt
            self%mass(2)=   32.971458760000_flyt
            self%mass(3)=   33.967866900000_flyt
            self%mass(4)=   35.967080760000_flyt
        case("Cl")
            self%n=2
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.7576000000_flyt
            self%conc(2)=0.2424000000_flyt
            self%mass(1)=   34.968852680000_flyt
            self%mass(2)=   36.965902590000_flyt
        case("Ar")
            self%n=3
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.0033650000_flyt
            self%conc(2)=0.0006320000_flyt
            self%conc(3)=0.9960030000_flyt
            self%mass(1)=   35.967545106000_flyt
            self%mass(2)=   37.962732400000_flyt
            self%mass(3)=   39.962383122500_flyt
        case("K")
            self%n=3
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.9325810000_flyt
            self%conc(2)=0.0001170000_flyt
            self%conc(3)=0.0673020000_flyt
            self%mass(1)=   38.963706680000_flyt
            self%mass(2)=   39.963998480000_flyt
            self%mass(3)=   40.961825760000_flyt
        case("Ca")
            self%n=6
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.9694100000_flyt
            self%conc(2)=0.0064700000_flyt
            self%conc(3)=0.0013500000_flyt
            self%conc(4)=0.0208600000_flyt
            self%conc(5)=0.0000400000_flyt
            self%conc(6)=0.0018700000_flyt
            self%mass(1)=   39.962590980000_flyt
            self%mass(2)=   41.958618010000_flyt
            self%mass(3)=   42.958766600000_flyt
            self%mass(4)=   43.955481800000_flyt
            self%mass(5)=   45.953692600000_flyt
            self%mass(6)=   47.952534000000_flyt
        case("Sc")
            self%n=1
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=1.0000000000_flyt
            self%mass(1)=   44.955911900000_flyt
        case("Ti")
            self%n=5
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.0825000000_flyt
            self%conc(2)=0.0744000000_flyt
            self%conc(3)=0.7372000000_flyt
            self%conc(4)=0.0541000000_flyt
            self%conc(5)=0.0518000000_flyt
            self%mass(1)=   45.952631600000_flyt
            self%mass(2)=   46.951763100000_flyt
            self%mass(3)=   47.947946300000_flyt
            self%mass(4)=   48.947870000000_flyt
            self%mass(5)=   49.944791200000_flyt
        case("V")
            self%n=2
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.0025000000_flyt
            self%conc(2)=0.9975000000_flyt
            self%mass(1)=   49.947158500000_flyt
            self%mass(2)=   50.943959500000_flyt
        case("Cr")
            self%n=4
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.0434500000_flyt
            self%conc(2)=0.8378900000_flyt
            self%conc(3)=0.0950100000_flyt
            self%conc(4)=0.0236500000_flyt
            self%mass(1)=   49.946044200000_flyt
            self%mass(2)=   51.940507500000_flyt
            self%mass(3)=   52.940649400000_flyt
            self%mass(4)=   53.938880400000_flyt
        case("Mn")
            self%n=1
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=1.0000000000_flyt
            self%mass(1)=   54.938045100000_flyt
        case("Fe")
            self%n=4
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.0584500000_flyt
            self%conc(2)=0.9175400000_flyt
            self%conc(3)=0.0211900000_flyt
            self%conc(4)=0.0028200000_flyt
            self%mass(1)=   53.939610500000_flyt
            self%mass(2)=   55.934937500000_flyt
            self%mass(3)=   56.935394000000_flyt
            self%mass(4)=   57.933275600000_flyt
        case("Co")
            self%n=1
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=1.0000000000_flyt
            self%mass(1)=   58.933195000000_flyt
        case("Ni")
            self%n=5
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.6807690000_flyt
            self%conc(2)=0.2622310000_flyt
            self%conc(3)=0.0113990000_flyt
            self%conc(4)=0.0363450000_flyt
            self%conc(5)=0.0092560000_flyt
            self%mass(1)=   57.935342900000_flyt
            self%mass(2)=   59.930786400000_flyt
            self%mass(3)=   60.931056000000_flyt
            self%mass(4)=   61.928345100000_flyt
            self%mass(5)=   63.927966000000_flyt
        case("Cu")
            self%n=2
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.6915000000_flyt
            self%conc(2)=0.3085000000_flyt
            self%mass(1)=   62.929597500000_flyt
            self%mass(2)=   64.927789500000_flyt
        case("Zn")
            self%n=5
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.4826800000_flyt
            self%conc(2)=0.2797500000_flyt
            self%conc(3)=0.0410200000_flyt
            self%conc(4)=0.1902400000_flyt
            self%conc(5)=0.0063100000_flyt
            self%mass(1)=   63.929142200000_flyt
            self%mass(2)=   65.926033400000_flyt
            self%mass(3)=   66.927127300000_flyt
            self%mass(4)=   67.924844200000_flyt
            self%mass(5)=   69.925319300000_flyt
        case("Ga")
            self%n=2
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.6010800000_flyt
            self%conc(2)=0.3989200000_flyt
            self%mass(1)=   68.925573600000_flyt
            self%mass(2)=   70.924701300000_flyt
        case("Ge")
            self%n=5
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.2038000000_flyt
            self%conc(2)=0.2731000000_flyt
            self%conc(3)=0.0776000000_flyt
            self%conc(4)=0.3672000000_flyt
            self%conc(5)=0.0783000000_flyt
            self%mass(1)=   69.924247400000_flyt
            self%mass(2)=   71.922075800000_flyt
            self%mass(3)=   72.923458900000_flyt
            self%mass(4)=   73.921177800000_flyt
            self%mass(5)=   75.921402600000_flyt
        case("As")
            self%n=1
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=1.0000000000_flyt
            self%mass(1)=   74.921596500000_flyt
        case("Se")
            self%n=6
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.0089000000_flyt
            self%conc(2)=0.0937000000_flyt
            self%conc(3)=0.0763000000_flyt
            self%conc(4)=0.2377000000_flyt
            self%conc(5)=0.4961000000_flyt
            self%conc(6)=0.0873000000_flyt
            self%mass(1)=   73.922476400000_flyt
            self%mass(2)=   75.919213600000_flyt
            self%mass(3)=   76.919914000000_flyt
            self%mass(4)=   77.917309100000_flyt
            self%mass(5)=   79.916521300000_flyt
            self%mass(6)=   81.916699400000_flyt
        case("Br")
            self%n=2
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.5069000000_flyt
            self%conc(2)=0.4931000000_flyt
            self%mass(1)=   78.918337100000_flyt
            self%mass(2)=   80.916290600000_flyt
        case("Kr")
            self%n=6
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.0035500000_flyt
            self%conc(2)=0.0228600000_flyt
            self%conc(3)=0.1159300000_flyt
            self%conc(4)=0.1150000000_flyt
            self%conc(5)=0.5698700000_flyt
            self%conc(6)=0.1727900000_flyt
            self%mass(1)=   77.920364800000_flyt
            self%mass(2)=   79.916379000000_flyt
            self%mass(3)=   81.913483600000_flyt
            self%mass(4)=   82.914136000000_flyt
            self%mass(5)=   83.911507000000_flyt
            self%mass(6)=   85.910610730000_flyt
        case("Rb")
            self%n=2
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.7217000000_flyt
            self%conc(2)=0.2783000000_flyt
            self%mass(1)=   84.911789738000_flyt
            self%mass(2)=   86.909180527000_flyt
        case("Sr")
            self%n=4
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.0056000000_flyt
            self%conc(2)=0.0986000000_flyt
            self%conc(3)=0.0700000000_flyt
            self%conc(4)=0.8258000000_flyt
            self%mass(1)=   83.913425000000_flyt
            self%mass(2)=   85.909260200000_flyt
            self%mass(3)=   86.908877100000_flyt
            self%mass(4)=   87.905612100000_flyt
        case("Y")
            self%n=1
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=1.0000000000_flyt
            self%mass(1)=   88.905848300000_flyt
        case("Zr")
            self%n=5
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.5145000000_flyt
            self%conc(2)=0.1122000000_flyt
            self%conc(3)=0.1715000000_flyt
            self%conc(4)=0.1738000000_flyt
            self%conc(5)=0.0280000000_flyt
            self%mass(1)=   89.904704400000_flyt
            self%mass(2)=   90.905645800000_flyt
            self%mass(3)=   91.905040800000_flyt
            self%mass(4)=   93.906315200000_flyt
            self%mass(5)=   95.908273400000_flyt
        case("Nb")
            self%n=1
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=1.0000000000_flyt
            self%mass(1)=   92.906378100000_flyt
        case("Mo")
            self%n=7
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.1477000000_flyt
            self%conc(2)=0.0923000000_flyt
            self%conc(3)=0.1590000000_flyt
            self%conc(4)=0.1668000000_flyt
            self%conc(5)=0.0956000000_flyt
            self%conc(6)=0.2419000000_flyt
            self%conc(7)=0.0967000000_flyt
            self%mass(1)=   91.906811000000_flyt
            self%mass(2)=   93.905088300000_flyt
            self%mass(3)=   94.905842100000_flyt
            self%mass(4)=   95.904679500000_flyt
            self%mass(5)=   96.906021500000_flyt
            self%mass(6)=   97.905408200000_flyt
            self%mass(7)=   99.907477000000_flyt
        case("Tc")
            self%n=1
            allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=1.0_flyt
            self%mass(1)=96.906_flyt
        case("Ru")
            self%n=7
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.0554000000_flyt
            self%conc(2)=0.0187000000_flyt
            self%conc(3)=0.1276000000_flyt
            self%conc(4)=0.1260000000_flyt
            self%conc(5)=0.1706000000_flyt
            self%conc(6)=0.3155000000_flyt
            self%conc(7)=0.1862000000_flyt
            self%mass(1)=   95.907598000000_flyt
            self%mass(2)=   97.905287000000_flyt
            self%mass(3)=   98.905939300000_flyt
            self%mass(4)=   99.904219500000_flyt
            self%mass(5)=  100.905582100000_flyt
            self%mass(6)=  101.904349300000_flyt
            self%mass(7)=  103.905433000000_flyt
        case("Rh")
            self%n=1
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=1.0000000000_flyt
            self%mass(1)=  102.905504000000_flyt
        case("Pd")
            self%n=6
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.0102000000_flyt
            self%conc(2)=0.1114000000_flyt
            self%conc(3)=0.2233000000_flyt
            self%conc(4)=0.2733000000_flyt
            self%conc(5)=0.2646000000_flyt
            self%conc(6)=0.1172000000_flyt
            self%mass(1)=  101.905609000000_flyt
            self%mass(2)=  103.904036000000_flyt
            self%mass(3)=  104.905085000000_flyt
            self%mass(4)=  105.903486000000_flyt
            self%mass(5)=  107.903892000000_flyt
            self%mass(6)=  109.905153000000_flyt
        case("Ag")
            self%n=2
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.5183900000_flyt
            self%conc(2)=0.4816100000_flyt
            self%mass(1)=  106.905097000000_flyt
            self%mass(2)=  108.904752000000_flyt
        case("Cd")
            self%n=8
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.0125000000_flyt
            self%conc(2)=0.0089000000_flyt
            self%conc(3)=0.1249000000_flyt
            self%conc(4)=0.1280000000_flyt
            self%conc(5)=0.2413000000_flyt
            self%conc(6)=0.1222000000_flyt
            self%conc(7)=0.2873000000_flyt
            self%conc(8)=0.0749000000_flyt
            self%mass(1)=  105.906459000000_flyt
            self%mass(2)=  107.904184000000_flyt
            self%mass(3)=  109.903002100000_flyt
            self%mass(4)=  110.904178100000_flyt
            self%mass(5)=  111.902757800000_flyt
            self%mass(6)=  112.904401700000_flyt
            self%mass(7)=  113.903358500000_flyt
            self%mass(8)=  115.904756000000_flyt
        case("In")
            self%n=2
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.0429000000_flyt
            self%conc(2)=0.9571000000_flyt
            self%mass(1)=  112.904058000000_flyt
            self%mass(2)=  114.903878000000_flyt
        case("Sn")
            self%n=10
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.0097000000_flyt
            self%conc(2)=0.0066000000_flyt
            self%conc(3)=0.0034000000_flyt
            self%conc(4)=0.1454000000_flyt
            self%conc(5)=0.0768000000_flyt
            self%conc(6)=0.2422000000_flyt
            self%conc(7)=0.0859000000_flyt
            self%conc(8)=0.3258000000_flyt
            self%conc(9)=0.0463000000_flyt
            self%conc(10)=0.0579000000_flyt
            self%mass(1)=  111.904818000000_flyt
            self%mass(2)=  113.902779000000_flyt
            self%mass(3)=  114.903342000000_flyt
            self%mass(4)=  115.901741000000_flyt
            self%mass(5)=  116.902952000000_flyt
            self%mass(6)=  117.901603000000_flyt
            self%mass(7)=  118.903308000000_flyt
            self%mass(8)=  119.902194700000_flyt
            self%mass(9)=  121.903439000000_flyt
            self%mass(10)=  123.905273900000_flyt
        case("Sb")
            self%n=2
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.5721000000_flyt
            self%conc(2)=0.4279000000_flyt
            self%mass(1)=  120.903815700000_flyt
            self%mass(2)=  122.904214000000_flyt
        case("Te")
            self%n=8
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.0009000000_flyt
            self%conc(2)=0.0255000000_flyt
            self%conc(3)=0.0089000000_flyt
            self%conc(4)=0.0474000000_flyt
            self%conc(5)=0.0707000000_flyt
            self%conc(6)=0.1884000000_flyt
            self%conc(7)=0.3174000000_flyt
            self%conc(8)=0.3408000000_flyt
            self%mass(1)=  119.904020000000_flyt
            self%mass(2)=  121.903043900000_flyt
            self%mass(3)=  122.904270000000_flyt
            self%mass(4)=  123.902817900000_flyt
            self%mass(5)=  124.904430700000_flyt
            self%mass(6)=  125.903311700000_flyt
            self%mass(7)=  127.904463100000_flyt
            self%mass(8)=  129.906224400000_flyt
        case("I")
            self%n=1
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=1.0000000000_flyt
            self%mass(1)=  126.904473000000_flyt
        case("Xe")
            self%n=9
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.0009520000_flyt
            self%conc(2)=0.0008900000_flyt
            self%conc(3)=0.0191020000_flyt
            self%conc(4)=0.2640060000_flyt
            self%conc(5)=0.0407100000_flyt
            self%conc(6)=0.2123240000_flyt
            self%conc(7)=0.2690860000_flyt
            self%conc(8)=0.1043570000_flyt
            self%conc(9)=0.0885730000_flyt
            self%mass(1)=  123.905893000000_flyt
            self%mass(2)=  125.904274000000_flyt
            self%mass(3)=  127.903531300000_flyt
            self%mass(4)=  128.904779400000_flyt
            self%mass(5)=  129.903508000000_flyt
            self%mass(6)=  130.905082400000_flyt
            self%mass(7)=  131.904153500000_flyt
            self%mass(8)=  133.905394500000_flyt
            self%mass(9)=  135.907219000000_flyt
        case("Cs")
            self%n=1
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=1.0000000000_flyt
            self%mass(1)=  132.905451933000_flyt
        case("Ba")
            self%n=7
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.0010600000_flyt
            self%conc(2)=0.0010100000_flyt
            self%conc(3)=0.0241700000_flyt
            self%conc(4)=0.0659200000_flyt
            self%conc(5)=0.0785400000_flyt
            self%conc(6)=0.1123200000_flyt
            self%conc(7)=0.7169800000_flyt
            self%mass(1)=  129.906320800000_flyt
            self%mass(2)=  131.905061300000_flyt
            self%mass(3)=  133.904508400000_flyt
            self%mass(4)=  134.905688600000_flyt
            self%mass(5)=  135.904575900000_flyt
            self%mass(6)=  136.905827400000_flyt
            self%mass(7)=  137.905247200000_flyt
        case("La")
            self%n=2
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.0009000000_flyt
            self%conc(2)=0.9991000000_flyt
            self%mass(1)=  137.907112000000_flyt
            self%mass(2)=  138.906353300000_flyt
        case("Ce")
            self%n=4
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.0018500000_flyt
            self%conc(2)=0.0025100000_flyt
            self%conc(3)=0.8845000000_flyt
            self%conc(4)=0.1111400000_flyt
            self%mass(1)=  135.907172000000_flyt
            self%mass(2)=  137.905991000000_flyt
            self%mass(3)=  139.905438700000_flyt
            self%mass(4)=  141.909244000000_flyt
        case("Pr")
            self%n=1
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=1.0000000000_flyt
            self%mass(1)=  140.907652800000_flyt
        case("Nd")
            self%n=7
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.2720000000_flyt
            self%conc(2)=0.1220000000_flyt
            self%conc(3)=0.2380000000_flyt
            self%conc(4)=0.0830000000_flyt
            self%conc(5)=0.1720000000_flyt
            self%conc(6)=0.0570000000_flyt
            self%conc(7)=0.0560000000_flyt
            self%mass(1)=  141.907723300000_flyt
            self%mass(2)=  142.909814300000_flyt
            self%mass(3)=  143.910087300000_flyt
            self%mass(4)=  144.912573600000_flyt
            self%mass(5)=  145.913116900000_flyt
            self%mass(6)=  147.916893000000_flyt
            self%mass(7)=  149.920891000000_flyt
        case("Pm")
            self%n=1
            allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=1.0_flyt
            self%mass(1)=144.91_flyt
        case("Sm")
            self%n=7
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.0307000000_flyt
            self%conc(2)=0.1499000000_flyt
            self%conc(3)=0.1124000000_flyt
            self%conc(4)=0.1382000000_flyt
            self%conc(5)=0.0738000000_flyt
            self%conc(6)=0.2675000000_flyt
            self%conc(7)=0.2275000000_flyt
            self%mass(1)=  143.911999000000_flyt
            self%mass(2)=  146.914897900000_flyt
            self%mass(3)=  147.914822700000_flyt
            self%mass(4)=  148.917184700000_flyt
            self%mass(5)=  149.917275500000_flyt
            self%mass(6)=  151.919732400000_flyt
            self%mass(7)=  153.922209300000_flyt
        case("Eu")
            self%n=2
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.4781000000_flyt
            self%conc(2)=0.5219000000_flyt
            self%mass(1)=  150.919850200000_flyt
            self%mass(2)=  152.921230300000_flyt
        case("Gd")
            self%n=7
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.0020000000_flyt
            self%conc(2)=0.0218000000_flyt
            self%conc(3)=0.1480000000_flyt
            self%conc(4)=0.2047000000_flyt
            self%conc(5)=0.1565000000_flyt
            self%conc(6)=0.2484000000_flyt
            self%conc(7)=0.2186000000_flyt
            self%mass(1)=  151.919791000000_flyt
            self%mass(2)=  153.920865600000_flyt
            self%mass(3)=  154.922622000000_flyt
            self%mass(4)=  155.922122700000_flyt
            self%mass(5)=  156.923960100000_flyt
            self%mass(6)=  157.924103900000_flyt
            self%mass(7)=  159.927054100000_flyt
        case("Tb")
            self%n=1
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=1.0000000000_flyt
            self%mass(1)=  158.925346800000_flyt
        case("Dy")
            self%n=7
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.0005600000_flyt
            self%conc(2)=0.0009500000_flyt
            self%conc(3)=0.0232900000_flyt
            self%conc(4)=0.1888900000_flyt
            self%conc(5)=0.2547500000_flyt
            self%conc(6)=0.2489600000_flyt
            self%conc(7)=0.2826000000_flyt
            self%mass(1)=  155.924283000000_flyt
            self%mass(2)=  157.924409000000_flyt
            self%mass(3)=  159.925197500000_flyt
            self%mass(4)=  160.926933400000_flyt
            self%mass(5)=  161.926798400000_flyt
            self%mass(6)=  162.928731200000_flyt
            self%mass(7)=  163.929174800000_flyt
        case("Ho")
            self%n=1
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=1.0000000000_flyt
            self%mass(1)=  164.930322100000_flyt
        case("Er")
            self%n=6
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.0013900000_flyt
            self%conc(2)=0.0160100000_flyt
            self%conc(3)=0.3350300000_flyt
            self%conc(4)=0.2286900000_flyt
            self%conc(5)=0.2697800000_flyt
            self%conc(6)=0.1491000000_flyt
            self%mass(1)=  161.928778000000_flyt
            self%mass(2)=  163.929200000000_flyt
            self%mass(3)=  165.930293100000_flyt
            self%mass(4)=  166.932048200000_flyt
            self%mass(5)=  167.932370200000_flyt
            self%mass(6)=  169.935464300000_flyt
        case("Tm")
            self%n=1
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=1.0000000000_flyt
            self%mass(1)=  168.934213300000_flyt
        case("Yb")
            self%n=7
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.0013000000_flyt
            self%conc(2)=0.0304000000_flyt
            self%conc(3)=0.1428000000_flyt
            self%conc(4)=0.2183000000_flyt
            self%conc(5)=0.1613000000_flyt
            self%conc(6)=0.3183000000_flyt
            self%conc(7)=0.1276000000_flyt
            self%mass(1)=  167.933897000000_flyt
            self%mass(2)=  169.934761800000_flyt
            self%mass(3)=  170.936325800000_flyt
            self%mass(4)=  171.936381500000_flyt
            self%mass(5)=  172.938210800000_flyt
            self%mass(6)=  173.938862100000_flyt
            self%mass(7)=  175.942571700000_flyt
        case("Lu")
            self%n=2
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.9741000000_flyt
            self%conc(2)=0.0259000000_flyt
            self%mass(1)=  174.940771800000_flyt
            self%mass(2)=  175.942686300000_flyt
        case("Hf")
            self%n=6
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.0016000000_flyt
            self%conc(2)=0.0526000000_flyt
            self%conc(3)=0.1860000000_flyt
            self%conc(4)=0.2728000000_flyt
            self%conc(5)=0.1362000000_flyt
            self%conc(6)=0.3508000000_flyt
            self%mass(1)=  173.940046000000_flyt
            self%mass(2)=  175.941408600000_flyt
            self%mass(3)=  176.943220700000_flyt
            self%mass(4)=  177.943698800000_flyt
            self%mass(5)=  178.945816100000_flyt
            self%mass(6)=  179.946550000000_flyt
        case("Ta")
            self%n=2
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.0001200000_flyt
            self%conc(2)=0.9998800000_flyt
            self%mass(1)=  179.947464800000_flyt
            self%mass(2)=  180.947995800000_flyt
        case("W")
            self%n=5
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.0012000000_flyt
            self%conc(2)=0.2650000000_flyt
            self%conc(3)=0.1431000000_flyt
            self%conc(4)=0.3064000000_flyt
            self%conc(5)=0.2843000000_flyt
            self%mass(1)=  179.946704000000_flyt
            self%mass(2)=  181.948204200000_flyt
            self%mass(3)=  182.950223000000_flyt
            self%mass(4)=  183.950931200000_flyt
            self%mass(5)=  185.954364100000_flyt
        case("Re")
            self%n=2
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.3740000000_flyt
            self%conc(2)=0.6260000000_flyt
            self%mass(1)=  184.952955000000_flyt
            self%mass(2)=  186.955753100000_flyt
        case("Os")
            self%n=7
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.0002000000_flyt
            self%conc(2)=0.0159000000_flyt
            self%conc(3)=0.0196000000_flyt
            self%conc(4)=0.1324000000_flyt
            self%conc(5)=0.1615000000_flyt
            self%conc(6)=0.2626000000_flyt
            self%conc(7)=0.4078000000_flyt
            self%mass(1)=  183.952489100000_flyt
            self%mass(2)=  185.953838200000_flyt
            self%mass(3)=  186.955750500000_flyt
            self%mass(4)=  187.955838200000_flyt
            self%mass(5)=  188.958147500000_flyt
            self%mass(6)=  189.958447000000_flyt
            self%mass(7)=  191.961480700000_flyt
        case("Ir")
            self%n=2
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.3730000000_flyt
            self%conc(2)=0.6270000000_flyt
            self%mass(1)=  190.960594000000_flyt
            self%mass(2)=  192.962926400000_flyt
        case("Pt")
            self%n=6
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.0001400000_flyt
            self%conc(2)=0.0078200000_flyt
            self%conc(3)=0.3296700000_flyt
            self%conc(4)=0.3383200000_flyt
            self%conc(5)=0.2524200000_flyt
            self%conc(6)=0.0716300000_flyt
            self%mass(1)=  189.959932000000_flyt
            self%mass(2)=  191.961038000000_flyt
            self%mass(3)=  193.962680300000_flyt
            self%mass(4)=  194.964791100000_flyt
            self%mass(5)=  195.964951500000_flyt
            self%mass(6)=  197.967893000000_flyt
        case("Au")
            self%n=1
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=1.0000000000_flyt
            self%mass(1)=  196.966568700000_flyt
        case("Hg")
            self%n=7
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.0015000000_flyt
            self%conc(2)=0.0997000000_flyt
            self%conc(3)=0.1687000000_flyt
            self%conc(4)=0.2310000000_flyt
            self%conc(5)=0.1318000000_flyt
            self%conc(6)=0.2986000000_flyt
            self%conc(7)=0.0687000000_flyt
            self%mass(1)=  195.965833000000_flyt
            self%mass(2)=  197.966769000000_flyt
            self%mass(3)=  198.968279900000_flyt
            self%mass(4)=  199.968326000000_flyt
            self%mass(5)=  200.970302300000_flyt
            self%mass(6)=  201.970643000000_flyt
            self%mass(7)=  203.973493900000_flyt
        case("Tl")
            self%n=2
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.2952000000_flyt
            self%conc(2)=0.7048000000_flyt
            self%mass(1)=  202.972344200000_flyt
            self%mass(2)=  204.974427500000_flyt
        case("Pb")
            self%n=4
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.0140000000_flyt
            self%conc(2)=0.2410000000_flyt
            self%conc(3)=0.2210000000_flyt
            self%conc(4)=0.5240000000_flyt
            self%mass(1)=  203.973043600000_flyt
            self%mass(2)=  205.974465300000_flyt
            self%mass(3)=  206.975896900000_flyt
            self%mass(4)=  207.976652100000_flyt
        case("Bi")
            self%n=1
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=1.0000000000_flyt
            self%mass(1)=  208.980398700000_flyt
        case("Po")
            self%n=1
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=1.0000000000_flyt
            self%mass(1)=208.98_flyt
        case("At")
            self%n=1
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=1.0000000000_flyt
            self%mass(1)=209.99_flyt
        case("Rn")
            self%n=1
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=1.0000000000_flyt
            self%mass(1)=222.02_flyt
        case("Fr")
            self%n=1
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=1.0000000000_flyt
            self%mass(1)=223.02_flyt
        case("Ra")
            self%n=1
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=1.0000000000_flyt
            self%mass(1)=226.03_flyt
        case("Ac")
            self%n=1
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=1.0000000000_flyt
            self%mass(1)=227.03_flyt
        case("Th")
            self%n=1
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=1.0000000000_flyt
            self%mass(1)=  232.038055300000_flyt
        case("Pa")
            self%n=1
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=1.0000000000_flyt
            self%mass(1)=  231.035884000000_flyt
        case("U")
            self%n=3
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.0000540000_flyt
            self%conc(2)=0.0072040000_flyt
            self%conc(3)=0.9927420000_flyt
            self%mass(1)=  234.040952100000_flyt
            self%mass(2)=  235.043929900000_flyt
            self%mass(3)=  238.050788200000_flyt
        case("Np")
            self%n=1
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=1.0000000000_flyt
            self%mass(1)=237.05_flyt
        case("Pu")
            self%n=1
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=1.0000000000_flyt
            self%mass(1)=244.06_flyt
        case("Am")
            self%n=1
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=1.0000000000_flyt
            self%mass(1)=243.06_flyt
        case("Cm")
            self%n=1
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=1.0000000000_flyt
            self%mass(1)=247.07_flyt
        case("Bk")
            self%n=1
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=1.0000000000_flyt
            self%mass(1)=247.07_flyt
        case("Cf")
            self%n=1
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=1.0000000000_flyt
            self%mass(1)=251.08_flyt
        case("Es")
            self%n=1
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=1.0000000000_flyt
            self%mass(1)=252.08_flyt
        case("Fm")
            self%n=1
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=1.0000000000_flyt
            self%mass(1)=257.1_flyt
        case("Md")
            self%n=1
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=1.0000000000_flyt
            self%mass(1)=258.1_flyt
        case("No")
            self%n=1
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=1.0000000000_flyt
            self%mass(1)=259.1_flyt
        case("Lr")
            self%n=1
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=1.0000000000_flyt
            self%mass(1)=262.11_flyt
        case default
            call lo_stop_gracefully(['no isotope distributions available for '//trim(symbol)//', are you sure it is a stable element?'],lo_exitcode_param,__FILE__,__LINE__)
    end select

    ! Convert to atomic units
    self%mass=self%mass*lo_amu_to_emu
    ! Get the actual mass
    self%conc=self%conc/sum(self%conc)
    self%mean_mass=sum(self%conc*self%mass)
    self%disorderparameter=self%mass_disorder_parameter()
end subroutine

end module
