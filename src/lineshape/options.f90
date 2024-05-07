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
    integer, dimension(3) :: qgrid_dos = -lo_hugeint
    real(r8), dimension(3) :: qpoint = -lo_huge
    real(r8) :: temperature = -lo_huge
    real(r8) :: sigma = -lo_huge
    real(r8) :: maxf = -lo_huge

    logical :: isotopescattering = .false.
    logical :: thirdorder = .false.
    logical :: fourthorder = .false.
    logical :: slightsmearing = .false.
    integer :: integrationtype = -lo_hugeint
    logical :: readiso = .false.
    integer :: verbosity = -lo_hugeint

    integer :: nq_on_path = -lo_hugeint
    logical :: readpathfromfile = .false.
    logical :: readqmesh = .false.
    integer :: meshtype = -lo_hugeint
    integer :: stride = -lo_hugeint
    integer :: minsmear = -lo_hugeint

    logical :: oneqpoint = .false.
    logical :: phonondos = .false.
    logical :: qpointpath = .false.
    logical :: timereversal = .false.
    logical :: diagonal = .false.
    logical :: exppath = .false.
    logical :: geninterp = .false.
    logical :: fancyinterp = .false.
    logical :: grid = .false.

    real(r8) :: cutoff = -lo_huge

    logical :: dielectric = .false.
    logical :: nosym = .false.
    real(r8), dimension(3) :: q_in = -lo_huge
    real(r8), dimension(3) :: q_out = -lo_huge
    logical :: dumpgrid = .false.
    logical :: convolution = .false.

    character(len=10) :: highsymmetrypoint = '**********'
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
    call cli%add(switch='--path', &
                 help='Calculate the self-energy and spectral function on a path through the BZ.', &
                 required=.false., act='store_true', def='.false.', error=lo_status)
    if (lo_status .ne. 0) stop
    cli_readpath
    cli_nq_on_path
    call cli%add(switch='--dos', hidden=.true., &
                 help='Calculate the broadened and shifted phonon DOS.', &
                 required=.false., act='store_true', def='.false.', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--dos_qpoint_grid', hidden=.true., &
                 help='Interpolate to a (preferrably) denser q-mesh when calculating the DOS.', &
                 nargs='3', required=.false., act='store', def='-1 -1 -1', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--qpoint', &
                 help='Calculate the self-energy at a single q-point, input in fractional coordinates.', &
                 nargs='3', required=.false., act='store', def='0 0 0', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--highsymmetrypoint', &
                 help='Samy as above, but you can specify the label of a high-symmetry point instead, e.g. "X" or "L".', &
                 required=.false., act='store', def='none', error=lo_status)
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
    call cli%add(switch='--exppath', &
                 help='Add experimental resolution function to a path.', hidden=.true., &
                 required=.false., act='store_true', def='.false.', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--stride', hidden=.true., &
                 help='Add a stride when calculating the lineshape along a path, for nicer-looking pictures.', &
                 required=.false., act='store', def='1', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--minsmear', hidden=.true., &
                 help='Add a baseline to the imaginary self-energy such that the plots look nicer', &
                 required=.false., act='store', def='-1', error=lo_status)
    if (lo_status .ne. 0) stop

    call cli%add(switch='--dielectric', &
                 help='Calculate dielectric things.', hidden=.true., &
                 required=.false., act='store_true', def='.false.', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--nosym', &
                 help='Switch off the use of symmetry.', hidden=.true., &
                 required=.false., act='store_true', def='.false.', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--qdirin', hidden=.false., &
                 help='Incident wavevector (Cartesian coordinates). Determins the behaviour of the non-analytical components at the zone center.', &
                 nargs='3', required=.false., act='store', def='1 0 0', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--qdirout', hidden=.true., &
                 help='Outgoing wavevector (Cartesian coordinates).', &
                 nargs='3', required=.false., act='store', def='1 0 0', error=lo_status)
    if (lo_status .ne. 0) stop

    call cli%add(switch='--dumpgrid', &
                 help='Write harmonic and anharmonic properties on a mesh to file.', hidden=.true., &
                 required=.false., act='store_true', def='.false.', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--convolution', &
                 help='Calculate somewhat self-consistent lineshapes and stuff.', hidden=.true., &
                 required=.false., act='store_true', def='.false.', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--geninterp', hidden=.true., &
                 help='First rule of interpolation is you do not talk about interpolation.', &
                 required=.false., act='store_true', def='.false.', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--cutoff', hidden=.true., &
                 help='Cutoff for forceconstant interpolation.', &
                 required=.false., act='store', def='5.0', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--fancyinterp', hidden=.true., &
                 help='Second rule of interpolation is you do not talk about interpolation.', &
                 required=.false., act='store_true', def='.false.', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--grid', &
                 help='Calculate spectral functions on a grid.', &
                 required=.false., act='store_true', def='.false.', error=lo_status)
    if (lo_status .ne. 0) stop

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
    call cli%get(switch='--nq_on_path', val=opts%nq_on_path)
    call cli%get(switch='--readpath', val=opts%readpathfromfile)
    call cli%get(switch='--qpoint_grid', val=opts%qgrid)
    call cli%get(switch='--dos_qpoint_grid', val=opts%qgrid_dos)
    call cli%get(switch='--sigma', val=opts%sigma)
    call cli%get(switch='--integrationtype', val=opts%integrationtype)
    call cli%get(switch='--dos', val=opts%phonondos)
    call cli%get(switch='--path', val=opts%qpointpath)
    call cli%get(switch='--temperature', val=opts%temperature, error=lo_status)
    call cli%get(switch='--readiso', val=opts%readiso)
    call cli%get(switch='--n_energies', val=opts%nf)
    call cli%get(switch='--qpoint', val=opts%qpoint)
    call cli%get(switch='--highsymmetrypoint', val=opts%highsymmetrypoint)
    call cli%get(switch='--max_energy', val=opts%maxf)
    call cli%get(switch='--meshtype', val=opts%meshtype)
    call cli%get(switch='--readqmesh', val=opts%readqmesh)
    call cli%get(switch='--stride', val=opts%stride)

    call cli%get(switch='--no_isotope_scattering', val=dumlog)
    opts%isotopescattering = .not. dumlog
    call cli%get(switch='--no_thirdorder_scattering', val=dumlog)
    opts%thirdorder = .not. dumlog
    call cli%get(switch='--fourthorder', val=opts%fourthorder)
    call cli%get(switch='--nondiagonal', val=dumlog)
    opts%diagonal = .not. dumlog
    call cli%get(switch='--minsmear', val=opts%minsmear)
    call cli%get(switch='--dielectric', val=opts%dielectric)
    call cli%get(switch='--nosym', val=opts%nosym)

    call cli%get(switch='--qdirin', val=opts%q_in)
    call cli%get(switch='--qdirout', val=opts%q_out)

    call cli%get(switch='--dumpgrid', val=opts%dumpgrid)
    call cli%get(switch='--convolution', val=opts%convolution)
    call cli%get(switch='--cutoff', val=opts%cutoff)
    call cli%get(switch='--fancyinterp', val=opts%fancyinterp)
    call cli%get(switch='--grid', val=opts%grid)
    ! Convert cutoff to atomic units right away.
    opts%cutoff = opts%cutoff*lo_A_to_bohr

    ! figure out what to do, and try to resolve strange conflicts
    opts%oneqpoint = .true.
    opts%phonondos = .false.
    opts%qpointpath = .false.
    opts%exppath = .false.
    opts%grid = .false.
    call cli%get(switch='--path', val=dumlog)
    if (dumlog) then
        opts%oneqpoint = .false.
        opts%phonondos = .false.
        opts%qpointpath = .true.
        opts%exppath = .false.
        opts%grid = .false.
        if (opts%qgrid_dos(1) .ne. -1) then
            write (*, *) 'Makes no sense to specify a dos q-grid and a path simultaneously.'
            stop
        end if
    end if
    call cli%get(switch='--dos', val=dumlog)
    if (dumlog) then
        opts%oneqpoint = .false.
        opts%phonondos = .true.
        opts%qpointpath = .false.
        opts%exppath = .false.
        opts%grid = .false.
        if (opts%stride .ne. 1) then
            write (*, *) 'Makes no sense to specify a stride when not calculating a path.'
            stop
        end if
        if (opts%qgrid_dos(1) .eq. -1) then
            opts%qgrid_dos = opts%qgrid
        end if
    end if
    call cli%get(switch='--geninterp', val=dumlog)
    if (dumlog) then
        opts%oneqpoint = .false.
        opts%phonondos = .false.
        opts%qpointpath = .false.
        opts%grid = .false.
        opts%geninterp = .true.
        if (opts%qgrid_dos(1) .eq. -1) then
            opts%qgrid_dos = opts%qgrid
        end if
    end if
    call cli%get(switch='--grid', val=dumlog)
    if (dumlog) then
        opts%oneqpoint = .false.
        opts%phonondos = .false.
        opts%qpointpath = .false.
        opts%geninterp = .false.
        opts%grid = .true.
        if (opts%qgrid_dos(1) .eq. -1) then
            opts%qgrid_dos = opts%qgrid
        end if
    end if

    if (opts%oneqpoint) then
        if (opts%stride .ne. 1) then
            write (*, *) 'Makes no sense to specify a stride when not calculating a path.'
            stop
        end if
        if (opts%qgrid_dos(1) .ne. -1) then
            write (*, *) 'Makes no sense to specify a dos q-grid and a single point simultaneously'
            stop
        end if
    end if

    ! avoid using --readpath when not calculating a path, see #64
    if (opts%readpathfromfile) then
        if (.not. opts%qpointpath) then
            write (*, *) '*** Makes no sense to specify --readpath when not calculating a path (--path).'
            stop
        end if
        if (opts%dumpgrid) then
            write (*, *) '*** Makes no sense to specify --readpath when calculating on grid (--dumpgrid).'
            stop
        end if
    end if

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
