#include "precompilerdefinitions"
module options
use konstanter, only: r8, lo_status, lo_huge, lo_hugeint, lo_author, lo_version, lo_licence, lo_twopi, lo_frequency_Hartree_to_meV, &
                      lo_frequency_Hartree_to_THz, lo_frequency_Hartree_to_icm, lo_A_to_bohr
use flap, only: command_line_interface
implicit none
private
public :: lo_opts

type lo_opts
    character(len=3) :: enhet = '***'
    integer :: nf = -lo_hugeint
    integer, dimension(3) :: qgrid = -lo_hugeint
    integer, dimension(3) :: qgrid_sigma = -lo_hugeint
    real(r8) :: temperature = -lo_huge
    real(r8) :: sigma = -lo_huge
    real(r8) :: maxf = -lo_huge

    logical :: polariton = .false.

    logical :: isotopescattering = .false.
    logical :: thirdorder = .false.
    logical :: fourthorder = .false.
    logical :: slightsmearing = .false.
    integer :: integrationtype = -lo_hugeint
    logical :: readiso = .false.
    integer :: verbosity = -lo_hugeint

    logical :: readqmesh = .false.
    integer :: meshtype = -lo_hugeint

    logical :: timereversal = .false.
    logical :: diagonal = .false.
    logical :: nosym = .false.
    logical :: convolution = .false.

    real(r8) :: unitfactor = -lo_huge
    character(len=4) :: unitname = '****'
contains
    procedure :: parse
end type

contains

subroutine parse(opts)
    !> the options
    class(lo_opts), intent(out) :: opts
    !> the helper parser
    type(command_line_interface) :: cli
    !
    logical :: dumlog

    ! basic info
    call cli%init(progname='lineshape', &
                  authors=lo_author, &
                  version=lo_version, &
                  license=lo_licence, &
                  help='Usage: ', &
                  description='Calculate the frequency-dependent self-energy and phonon spectral function.', &
                  examples=["mpirun lineshape --highsymmetrypoint X --temperature 500      ", &
                            "mpirun lineshape --path -qg 10 10 10 --noisotope              ", &
                            "mpirun lineshape --dos -qg 10 10 10 --dos_qpoint_grid 24 24 24"], &
                  epilog=new_line('a')//"...")

    cli_unit
    cli_temperature
    call cli%add(switch='--n_energies', switch_ab='-ne', &
                 help='Number of energies for the energy-dependent self-energy.', &
                 required=.false., act='store', def='1200', error=lo_status)
    if (lo_status .ne. 0) stop
    cli_qpoint_grid
    call cli%add(switch='--meshtype', help='Type of q-point mesh. 1 Is a Monkhorst-Pack mesh, 2 an FFT mesh and 3 my fancy wedge-based mesh with approximately the same density the grid-based meshes. 4 build the commensurate mesh of an approximately cubic supercell.', &
                 required=.false., act='store', def='2', choices='1,2,3', error=lo_status)
    call cli%add(switch='--integrationtype', switch_ab='-it', &
                 help='Type of integration for phase space integrals. 1 is Gaussian, 2 adaptive Gaussian and 3 Tetrahedron.', &
                 required=.false., act='store', def='2', choices='1,2,3,4,5', error=lo_status)
    if (lo_status .ne. 0) stop
    cli_sigma
    call cli%add(switch='--sigma_qpoint_grid',hidden=.false., &
                 help='Mesh the self-energy is evaluated on.', &
                 nargs='3', required=.false., act='store', def='4 4 4', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--max_energy', &
                 help='Maximum energy where the output is cut off, in multiples of the maximum harmonic frequency.', &
                 required=.false., act='store', def='1.4', error=lo_status)
    if (lo_status .ne. 0) stop
    cli_no_isotope_scattering
    call cli%add(switch='--no_thirdorder_scattering', &
                 help='Switch of three-phonon scattering', &
                 required=.false., act='store_true', def='.false.', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--fourthorder', &
                 help='Consider four-phonon contributions to the real part of the self-energy.', hidden=.true., &
                 required=.false., act='store_true', def='.false.', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--nondiagonal', &
                 help='Consider non-diagonal contributions to the self-energy.', hidden=.true., &
                 required=.false., act='store_true', def='.false.', error=lo_status)
    if (lo_status .ne. 0) stop

    call cli%add(switch='--convolution', &
                 help='Calculate somewhat self-consistent lineshapes and stuff.', hidden=.true., &
                 required=.false., act='store_true', def='.false.', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--polariton', &
                 help='Calculate the self-energy for phonon polaritons.', &
                 required=.false., act='store_true', def='.false.', error=lo_status)

    cli_readiso
    cli_readqmesh
    cli_manpage
    cli_verbose

    ! actually parse it
    call cli%parse(error=lo_status)
    if (lo_status .ne. 0) stop
    ! Should the manpage be generated? In that case, no point to go further than here.
    dumlog = .false.
    call cli%get(switch='--manpage', val=dumlog)
    if (dumlog) then
        call cli%save_man_page(trim(cli%progname)//'.1')
        call cli%save_usage_to_markdown(trim(cli%progname)//'.md')
        write (*, *) 'Wrote manpage for "'//trim(cli%progname)//'"'
        stop
    end if
    opts%verbosity = 0
    call cli%get(switch='--verbose', val=dumlog); if (dumlog) opts%verbosity = 2

    ! get real options
    call cli%get(switch='--unit', val=opts%enhet)
    call cli%get(switch='--qpoint_grid', val=opts%qgrid)
    call cli%get(switch='--sigma_qpoint_grid', val=opts%qgrid_sigma)
    call cli%get(switch='--sigma', val=opts%sigma)
    call cli%get(switch='--integrationtype', val=opts%integrationtype)
    call cli%get(switch='--temperature', val=opts%temperature, error=lo_status)
    call cli%get(switch='--readiso', val=opts%readiso)
    call cli%get(switch='--n_energies', val=opts%nf)
    call cli%get(switch='--max_energy', val=opts%maxf)
    call cli%get(switch='--meshtype', val=opts%meshtype)
    call cli%get(switch='--readqmesh', val=opts%readqmesh)
    call cli%get(switch='--polariton', val=opts%polariton)

    call cli%get(switch='--no_isotope_scattering', val=dumlog)
    opts%isotopescattering = .not. dumlog
    call cli%get(switch='--no_thirdorder_scattering', val=dumlog)
    opts%thirdorder = .not. dumlog
    call cli%get(switch='--fourthorder', val=opts%fourthorder)
    call cli%get(switch='--nondiagonal', val=dumlog)
    opts%diagonal = .not. dumlog
    !call cli%get(switch='--nosym', val=opts%nosym)

    call cli%get(switch='--convolution', val=opts%convolution)

    ! Not really options
    opts%timereversal = .true.
    opts%slightsmearing = .false.

    ! can keep the unit factor here, quite convenient
    select case (opts%enhet)
    case ('thz')
        opts%unitfactor = lo_frequency_Hartree_to_THz
        opts%unitname = 'Thz'
    case ('icm')
        opts%unitfactor = lo_frequency_Hartree_to_icm
        opts%unitname = 'cm-1'
    case ('mev')
        opts%unitfactor = lo_frequency_Hartree_to_meV
        opts%unitname = 'meV'
    case default
        write (*, *) "Unknown unit, try 'thz', 'mev' or 'icm'"
        stop
    end select
end subroutine

end module
