#include "precompilerdefinitions"
module options
use konstanter, only: flyt, lo_status, lo_author, lo_version, lo_licence, &
                      lo_A_to_bohr, lo_huge, lo_hugeint, lo_frequency_THz_to_Hartree
use flap, only: command_line_interface
implicit none
private
public :: lo_opts

type lo_opts
    ! cutoffs
    real(flyt) :: cutoff2 = -lo_huge, cutoff3 = -lo_huge, cutoff4 = -lo_huge
    ! consider first order force constants
    logical :: firstorder = .false.
    ! use WZ cutoff for second order?
    integer, dimension(3) :: wzdimensions = -lo_hugeint
    ! use graph jumps for the cutoff?
    integer :: njump2 = -lo_hugeint, njump3 = -lo_hugeint, njump4 = -lo_hugeint
    ! minimal coordination numbers per atom
    integer, dimension(:), allocatable :: coordnumber
    ! how much to talk
    integer :: verbosity = -lo_hugeint
    ! Info about solver
    integer :: solver = -lo_hugeint
    ! stride, use every n timesteps
    integer :: stride = -lo_hugeint
    ! read the forcemap from file
    logical :: readforcemap = .false.
    ! read the irreducible components from file
    logical :: readirreducible = .false.
    ! dump the potential energy differences
    logical :: ediff = .false.
    ! relax the structure using the first order forceconstants
    logical :: relax = .false.
    ! use rotational constraints on the forceconstant
    logical :: rotationalconstraints = .false.
    ! use Huang invariances
    logical :: huanginvariance = .false.
    ! use spacegroup symmetries
    logical :: spacegroup = .false.
    ! use transposition symmetries
    logical :: transposition = .false.
    ! enforce Hermitian character
    logical :: hermitian = .false.
    ! do the second order fit without reference positions
    logical :: noreference = .false.
    ! is this a polar material?
    logical :: polar = .false.
    ! print the forcemap
    logical :: printforcemap = .false.
    ! what kind of polar correction?
    integer :: polarcorrectiontype = -lo_hugeint
    ! cutoff for magnetic pair interactions
    real(flyt) :: magcutoff2 = -lo_huge
    ! on-site magnetic term?
    logical :: magnetic_onsite = .false.
    ! use only N random steps?
    integer :: nrandom = -lo_hugeint
    ! largest allowed mean square displacement
    real(flyt) :: max_displacement = -lo_huge
    ! cutoff for dielectric pair interactions
    real(flyt) :: dielcutoff2 = -lo_huge
    ! cutoff for dielectric triplet interactions
    real(flyt) :: dielcutoff3 = -lo_huge
    ! temperature
    real(flyt) :: temperature = -lo_huge
    ! dump fake dielectric interactions
    logical :: fake_dielectric = .false.
    ! developer mode
    logical :: devmode = .false.
    ! extract realspace pressure from force constants
    logical :: pressure = .false.

contains
    procedure :: parse
end type

contains

!> parse the command line arguments and set defaults
subroutine parse(opts)
    !> the options
    class(lo_opts), intent(out) :: opts
    !> the helper parser
    type(command_line_interface) :: cli
    !
    logical :: dumlog

    call cli%init(progname='extract_forceconstants', &
                  authors=lo_author, &
                  version=lo_version, &
                  license=lo_licence, &
                  help='Usage: ', &
                  description='The main algorithm of the TDEP method. &
                              &Starting with a symmetry analysis, this code finds the irreducible &
                              &representation of interatomic forceconstants and extracts them &
                              &from position and force data.', &
                  examples=["extract_forceconstants -rc2 5.1           ", &
                            "extract_forceconstants -rc2 4.5 -rc3 3.21 "], &
                  epilog=new_line('a')//"...")

    ! real options
    call cli%add(switch='--secondorder_cutoff', switch_ab='-rc2', &
                 help='Cutoff for the second order force constants', &
                 required=.false., act='store', def='5.0', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--thirdorder_cutoff', switch_ab='-rc3', &
                 help='Cutoff for the third order force constants', &
                 required=.false., act='store', def='-1', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--fourthorder_cutoff', switch_ab='-rc4', &
                 help='Cutoff for the fourth order force constants', &
                 required=.false., act='store', def='-1', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--polar', &
                 help='Add dipole-dipole corrections for polar materials.', &
                 required=.false., act='store_true', def='.false.', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--stride', switch_ab='-s', &
                 help='Use every N configuration instead of all. Useful for long MD &
                      &simulations with linearly dependent configurations.', &
                 required=.false., act='store', def='1', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--firstorder', &
                 help='Include the first order force constants. &
                      &These can be used to find the finite temperature equilibrium structure.', &
                 required=.false., act='store_true', def='.false.', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--readforcemap',hidden=.true., &
                 help='Read `infile.forcemap.hdf5` from file instead of calculating &
                      &all symmetry relations. Useful for sets of calculations with the same structure.', &
                 required=.false., act='store_true', def='.false.', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--readirreducible',hidden=.true., &
                 help='Read the irreducible forceconstants from `infile.irrifc_*` &
                      &instead of solving for them. This option requires an `infile.forcemap.hdf5`, as above.', &
                 required=.false., act='store_true', def='.false.', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--potential_energy_differences', switch_ab='-U0',hidden=.true., &
                 help='Calculate the difference in potential energy from the &
                      &simulation and the forceconstants to determine U0.', &
                 help_markdown='As referenced in the thermodynamics section of &
                               &[phonon dispersion relations](phonon_dispersion_relations.html) &
                               &this is the renormalized baseline for the TDEP free energy: &
                               &$$U_0= \left\langle U^{\textrm{BO}}(t)-\frac{1}{2} &
                               &\sum_{ij}\sum_{\alpha\beta} \Phi_{ij}^{\alpha\beta} &
                               &\mathbf{u}^{\alpha}_i(t) \mathbf{u}^{\beta}_j(t) \right\rangle$$ &
                               &This number should be added to the appropriate phonon free energy.', &
                 required=.false., act='store_true', def='.false.', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--printforcemap',hidden=.true., &
                 help='Print `outfile.forcemap.hdf5` for reuse.', &
                 required=.false., act='store_true', def='.false.', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--temperature',&
                 help='Temperature for self-consistent solver.',&
                 required=.false., act='store', def='-1', error=lo_status)
    if (lo_status .ne. 0) stop

    ! hidden useless options
    call cli%add(switch='--manpage', hidden=.true., help='', &
                 required=.false., act='store_true', def='.false.', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--verbose', hidden=.true., help='', &
                 required=.false., act='store_true', def='.false.', error=lo_status)
    if (lo_status .ne. 0) stop
    ! hidden and TOP SECRET!
    call cli%add(switch='--wigner_seitz_cutoff', switch_ab='-wz', hidden=.true., &
                 help='Use the Voronoi cell of the supercell, a super Wigner-Seitz cell, &
                      &for the second order cutoff. The three numbers are the dimensions of the supercell.', &
                 nargs='3', required=.false., act='store', def='-1 -1 -1', exclude='--secondorder_njump', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--solver', hidden=.true., &
                 required=.false., act='store', def='1', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--relax', hidden=.true., help='', &
                 required=.false., act='store_true', def='.false.', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--notranspose', hidden=.true., help='', &
                 required=.false., act='store_true', def='.false.', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--nospacegroup', hidden=.true., help='', &
                 required=.false., act='store_true', def='.false.', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--norotational', hidden=.false., help='Turn off imposing rotational invariance. Needed for 2D systems.', &
                 required=.false., act='store_true', def='.false.', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--nohuang', hidden=.false., help='Turn off imposing Huang invariances. Useful for 2D systems.', &
                 required=.false., act='store_true', def='.false.', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--nohermitian', hidden=.false., help='', &
                 required=.false., act='store_true', def='.false.', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--noreference', hidden=.true., help='', &
                 required=.false., act='store_true', def='.false.', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--polarcorrectiontype', switch_ab='-pc', hidden=.true., &
                 help='What kind of polar correction to use.', &
                 required=.false., act='store', def='3', choices='1,2,3', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--secondorder_njump', switch_ab='-nj2', &
                 help='Second order neighbour jumps.', &
                 required=.false., act='store', def='-1', exclude='-wzdim', error=lo_status, hidden=.true.)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--thirdorder_njump', switch_ab='-nj3', &
                 help='Third order neighbour jumps', &
                 required=.false., act='store', def='-1', exclude='--thirdorder_cutoff', error=lo_status, hidden=.true.)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--fourthorder_njump', switch_ab='-nj4', &
                 help='Fourth order neighbour jumps', &
                 required=.false., act='store', def='-1', exclude='--fourthorder_cutoff', error=lo_status, hidden=.true.)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--coordination_number', switch_ab='-cn', &
                 help='In conjuntion with jump-based cutoffs, manually define the &
                      &coordination per atom. One integer per atom in the unit cell.', &
                 help_markdown='This will count the N shortest distances from each &
                               &atom as belonging to the first coordination shell. &
                               &If the desired coordination is impossible, e.g. &
                               &coordination 5 in fcc, the code will stop.', &
                 required=.false., nargs='+', act='store', def='-1', error=lo_status, hidden=.true.)

    call cli%add(switch='--magnetic_pair_cutoff', switch_ab='-mc2', hidden=.true., &
                 help='Cutoff for the pair magnetic interactions', &
                 required=.false., act='store', def='-1.0', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--magnetic_onsite', hidden=.true., &
                 help='Add on on-site magnetic term for longitudinal fluctuations.', &
                 required=.false., act='store_true', def='.false.', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--nrand', hidden=.true., &
                 help='Use a random subset of N configuration instead of all.', &
                 required=.false., act='store', def='-1', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--maxdisp', hidden=.true., &
                 help='Largest allowed mean square displacement.', &
                 required=.false., act='store', def='-1', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--dielectric_pair_cutoff', switch_ab='-dc2', hidden=.true., &
                 help='Cutoff for the pair dielectric interactions', &
                 required=.false., act='store', def='-1.0', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--dielectric_triplet_cutoff', switch_ab='-dc3', hidden=.true., &
                 help='Cutoff for the triplet dielectric interactions', &
                 required=.false., act='store', def='-1.0', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--fakediel', hidden=.true., &
                 help='Dump arbitrary fake Raman and IR tensor that still obey the proper symmetry.', &
                 required=.false., act='store_true', def='.false.', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--developermode', switch_ab='-dev', hidden=.true., help='dev. mode', &
                 required=.false., act='store_true', def='.false.', error=lo_status)
    call cli%add(switch='--pressure', &
                 help='Extract TDEP pressure from samples for lattice expansion etc.', &
                 required=.false., act='store_true', def='.false.', error=lo_status)
    if (lo_status .ne. 0) stop

    ! actually parse it
    call cli%parse(error=lo_status)
    if (lo_status .ne. 0) stop

    ! generate manpage?
    call cli%get(switch='--manpage', val=dumlog)
    if (dumlog) then
        call cli%save_man_page(trim(cli%progname)//'.1')
        call cli%save_usage_to_markdown(trim(cli%progname)//'.md')
        write (*, *) 'Wrote manpage for "'//trim(cli%progname)//'"'
        stop
    end if

    ! verbose output?
    opts%verbosity = 0
    call cli%get(switch='--verbose', val=dumlog)
    if (dumlog) opts%verbosity = 1

    ! Parse the rest of stuff
    call cli%get(switch='-rc2', val=opts%cutoff2)
    call cli%get(switch='-rc3', val=opts%cutoff3)
    call cli%get(switch='-rc4', val=opts%cutoff4)
    call cli%get(switch='-mc2', val=opts%magcutoff2)
    call cli%get(switch='--magnetic_onsite', val=opts%magnetic_onsite)
    call cli%get(switch='--firstorder', val=opts%firstorder)
    call cli%get(switch='-wz', val=opts%wzdimensions)
    call cli%get(switch='--solver', val=opts%solver)
    call cli%get(switch='--printforcemap', val=opts%printforcemap)
    call cli%get(switch='--readforcemap', val=opts%readforcemap)
    call cli%get(switch='--readirreducible', val=opts%readirreducible)
    call cli%get(switch='-U0', val=opts%ediff)
    call cli%get(switch='--relax', val=opts%relax)
    call cli%get(switch='--stride', val=opts%stride)
    call cli%get(switch='--norotational', val=dumlog)
    opts%rotationalconstraints = .not. dumlog
    call cli%get(switch='--nohuang', val=dumlog)
    opts%huanginvariance = .not. dumlog
    call cli%get(switch='--notranspose', val=dumlog)
    opts%transposition = .not. dumlog
    call cli%get(switch='--nospacegroup', val=dumlog)
    opts%spacegroup = .not. dumlog
    call cli%get(switch='--nohermitian', val=dumlog)
    opts%hermitian = .not. dumlog
    call cli%get(switch='--noreference', val=opts%noreference)
    call cli%get(switch='-nj2', val=opts%njump2)
    call cli%get(switch='-nj3', val=opts%njump3)
    call cli%get(switch='-nj4', val=opts%njump4)
    call cli%get_varying(switch='-cn', val=opts%coordnumber)
    call cli%get(switch='--polar', val=opts%polar)
    call cli%get(switch='--polarcorrectiontype', val=opts%polarcorrectiontype, error=lo_status)
    call cli%get(switch='--nrand', val=opts%nrandom)
    call cli%get(switch='--maxdisp', val=opts%max_displacement)
    call cli%get(switch='-dc2', val=opts%dielcutoff2)
    call cli%get(switch='-dc3', val=opts%dielcutoff3)
    call cli%get(switch='--temperature', val=opts%temperature)
    call cli%get(switch='--fakediel', val=opts%fake_dielectric)
    call cli%get(switch='--developermode', val=opts%devmode)
    call cli%get(switch='--pressure', val=opts%pressure)

    if (lo_status .ne. 0) stop

    ! Try to resolve conflicting options:
    if (opts%readirreducible) then
        opts%readforcemap = .true.
    end if
    if (opts%njump2 .gt. 0 .and. opts%wzdimensions(1) .gt. 0) then

    end if
    if (opts%relax) opts%firstorder = .true.

    ! Unit conversions of the input to atomic units right away
    opts%cutoff2 = opts%cutoff2*lo_A_to_Bohr
    opts%cutoff3 = opts%cutoff3*lo_A_to_Bohr
    opts%cutoff4 = opts%cutoff4*lo_A_to_Bohr
    opts%magcutoff2 = opts%magcutoff2*lo_A_to_Bohr
    opts%dielcutoff2 = opts%dielcutoff2*lo_A_to_Bohr
    opts%dielcutoff3 = opts%dielcutoff3*lo_A_to_Bohr

end subroutine

end module
