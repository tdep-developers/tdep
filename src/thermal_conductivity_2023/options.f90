#include "precompilerdefinitions"
module options
use konstanter, only: flyt, lo_status, lo_author, lo_version, lo_licence, lo_m_to_bohr
use flap, only: command_line_interface
implicit none
private
public :: lo_opts

type lo_opts
    integer, dimension(3) :: qgrid  !< the main q-grid
    logical :: readqmesh            !< read q-grid from file
    integer :: trangenpts           !< how many temperatures
    real(flyt) :: trangemin         !< minimum temperature
    real(flyt) :: trangemax         !< max temperature
    logical :: logtempaxis          !< logarithmically spaced temperature points
    real(flyt) :: sigma             !< scaling factor for adaptice gaussian
    real(flyt) :: thres             !< consider Gaussian 0 if x-mu is larger than this number times sigma.
    real(flyt) :: tau_boundary      !< add a constant as boundary scattering
    real(flyt) :: mfp_max           !< add a length as boundary scattering
    logical :: readiso              !< read isotope distribution from file
    integer :: integrationtype      !< gaussian or tetrahedron
    integer :: scfiterations        !< maximum number of self-consistent iterations
    real(flyt) :: scftol            !< tolerance for the SCF cycle

    integer :: correctionlevel      !< how hard to correct
    integer :: mfppts               !< number of points on mfp-plots
    logical :: dumpgrid             !< print everything on a grid
    !logical :: thinfilm             !< Austins thin film thing

    ! Debugging things
    logical :: timereversal
    logical :: qpsymmetry
    logical :: isotopescattering
    !
    integer :: verbosity
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
    real(flyt) :: f0
    real(flyt), dimension(3) :: dumflytv

    ! basic info
    call cli%init(progname='thermal_conductivity_2023', &
                  authors=lo_author, &
                  version=lo_version, &
                  license=lo_licence, &
                  help='Usage: ', &
                  description='Calculates the lattice thermal conductivity from the iterative solution of the &
                              &phonon Boltzmann equation. In addition, cumulative plots and raw data dumps &
                              &of intermediate values are available.', &
                  examples=["mpirun thermal_conductivity_2023 --temperature 300                              ", &
                            "mpirun thermal_conductivity_2023 -qg 15 15 15 --temperature_range 200 600 50    ", &
                            "mpirun thermal_conductivity_2023 --integrationtype 2 -qg 30 30 30 --max_mfp 1E-6"], &
                  epilog=new_line('a')//"...")
    ! real options
    call cli%add(switch='--readiso', &
                 help='Read the isotope distribution from `infile.isotopes`.', &
                 help_markdown='The format is specified [here](../page/files.html#infile.isotopes).', &
                 required=.false., act='store_true', def='.false.', error=lo_status)
    if (lo_status .ne. 0) stop
    cli_qpoint_grid
    call cli%add(switch='--integrationtype', switch_ab='-it', &
                 help='Type of integration for the phonon DOS. 1 is Gaussian, 2 adaptive Gaussian and 3 Tetrahedron.', &
                 required=.false., act='store', def='2', choices='1,2,3', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--sigma', &
                 help='Global scaling factor for adaptive Gaussian smearing.', &
                 required=.false., act='store', def='1.0', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--threshold', &
                 help='Consider a Gaussian distribution to be 0 after this many standard deviations.', &
                 required=.false., act='store', def='4.0', error=lo_status)
    if (lo_status .ne. 0) stop
    cli_readqmesh

    call cli%add(switch='--temperature', &
                 help='Evaluate thermal conductivity at a single temperature.', &
                 required=.false., act='store', def='-1', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--temperature_range', &
                 help='Series of temperatures for thermal conductivity. Specify min, max and the number of points.', &
                 nargs='3', required=.false., act='store', def='100 300 5', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--logtempaxis', &
                 help='Space the temperature points logarithmically instead of linearly.', &
                 required=.false., act='store_true', def='.false.', error=lo_status)
    call cli%add(switch='--max_mfp', &
                 help='Add a limit on the mean free path as an approximation of domain size.', &
                 required=.false., act='store', def='-1', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--dumpgrid', &
                 help='Write files with q-vectors, frequencies, eigenvectors and group velocities for a grid.', &
                 required=.false., act='store_true', def='.false.', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--noisotope', &
                 help='Do not consider isotope scattering.', &
                 required=.false., act='store_true', def='.false.', error=lo_status)
    if (lo_status .ne. 0) stop

    ! hidden
    call cli%add(switch='--tau_boundary', hidden=.true., &
                 help='Add a constant boundary scattering term to the lifetimes.', &
                 required=.false., act='store', def='-1', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--mfppts', hidden=.true., help='', &
                 required=.false., act='store', def='200', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--scfiterations', hidden=.true., help='', &
                 required=.false., act='store', def='200', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--notr', hidden=.true., help='', &
                 required=.false., act='store_true', def='.false.', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--nosym', hidden=.true., help='', &
                 required=.false., act='store_true', def='.false.', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--correctionlevel', hidden=.true., &
                 help='How agressively things are corrected due to broken symmetries.', &
                 required=.false., act='store', def='4', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--scftol', hidden=.true., &
                 help='What tolerance to converge the self-consistent cycle to.', &
                 required=.false., act='store', def='1E-5', error=lo_status)
    if (lo_status .ne. 0) stop
    !call cli%add(switch='--thinfilm',hidden=.true.,&
    !        help='Calculate the suppression of kappa from in a thin film.',&
    !        required=.false.,act='store_true',def='.false.',error=lo_status)
    !        if ( lo_status .ne. 0 ) stop
    cli_manpage
    cli_verbose

    ! actually parse it
    call cli%parse(error=lo_status)
    if (lo_status .ne. 0) stop
    !
    ! Should the manpage be generated? In that case, no point to go further than here.
    !
    dumlog = .false.
    call cli%get(switch='--manpage', val=dumlog)
    if (dumlog) then
        call cli%save_man_page(trim(cli%progname)//'.1')
        call cli%save_usage_to_markdown(trim(cli%progname)//'.md')
        write (*, *) 'Wrote manpage for "'//trim(cli%progname)//'"'
        stop
    end if

    ! store things in the right place

    call cli%get(switch='--temperature', val=f0)
    call cli%get(switch='--temperature_range', val=dumflytv)
    opts%trangemin = dumflytv(1)
    opts%trangemax = dumflytv(2)
    opts%trangenpts = int(anint(dumflytv(3)))
    ! if --temperature is specified, override the range
    if (f0 .gt. 0.0_flyt) then
        opts%trangemin = f0
        opts%trangemax = f0
        opts%trangenpts = 1
    end if
    call cli%get(switch='--qpoint_grid', val=opts%qgrid)
    call cli%get(switch='--sigma', val=opts%sigma)
    call cli%get(switch='--threshold', val=opts%thres)
    call cli%get(switch='--tau_boundary', val=opts%tau_boundary)
    if (opts%tau_boundary .gt. 0.0_flyt) opts%tau_boundary = 1E10_flyt
    call cli%get(switch='--readqmesh', val=opts%readqmesh)
    call cli%get(switch='--integrationtype', val=opts%integrationtype)
    call cli%get(switch='--logtempaxis', val=opts%logtempaxis)
    call cli%get(switch='--readiso', val=opts%readiso)
    call cli%get(switch='--mfppts', val=opts%mfppts)
    call cli%get(switch='--scfiterations', val=opts%scfiterations)
    call cli%get(switch='--max_mfp', val=opts%mfp_max)
    call cli%get(switch='--dumpgrid', val=opts%dumpgrid)
    !call cli%get(switch='--thinfilm',val=opts%thinfilm)
    ! stuff that's not really an option
    call cli%get(switch='--correctionlevel', val=opts%correctionlevel)
    call cli%get(switch='--scftol', val=opts%scftol)
    call cli%get(switch='--notr', val=dumlog)
    opts%timereversal = .not. dumlog
    call cli%get(switch='--nosym', val=dumlog)
    opts%qpsymmetry = .not. dumlog
    call cli%get(switch='--verbose', val=dumlog)
    if (dumlog) then
        opts%verbosity = 2
    else
        opts%verbosity = 0
    end if
    call cli%get(switch='--noisotope', val=dumlog)
    opts%isotopescattering = .not. dumlog

    ! Get things to atomic units
    opts%mfp_max = opts%mfp_max*lo_m_to_Bohr

end subroutine

end module
