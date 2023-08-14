#include "precompilerdefinitions"
module options
use konstanter, only: flyt, lo_author, lo_version, lo_licence, lo_tiny, lo_status, lo_huge, lo_hugeint
use flap, only: command_line_interface
implicit none
private
public :: lo_opts
type lo_opts
    ! common
    integer :: verbosity = -lo_hugeint
    logical :: readiso = .false.
    character(len=3) :: enhet = '   '
    ! for dispersions
    integer :: nq = -lo_hugeint
    logical :: readpathfromfile = .false.
    logical :: gruneisen = .false.
    logical :: fullmesh = .false.
    ! dos stuff
    logical :: dos = .false.
    integer, dimension(3) :: qgrid = -lo_hugeint
    real(flyt) :: sigma = -lo_huge
    integer :: dosintegrationtype = -lo_hugeint
    integer :: dospoints = -lo_hugeint
    logical :: readqmesh = .false.
    integer :: meshtype = -lo_hugeint
    logical :: dumpgrid = .false.
    ! energy stuff
    real(flyt) :: temperature = -lo_huge
    real(flyt) :: trangemin = -lo_huge
    real(flyt) :: trangemax = -lo_huge
    integer :: trangenpts = -lo_hugeint
    ! hidden debugging things
    logical :: timereversal = .false.
    logical :: sortbands = .false.
    logical :: dumpsecret = .false.
    logical :: unfold = .false.
    logical :: inelastic = .false.
    logical :: pdfplot = .false.
    ! hidden unfolding options
    integer :: nenergy = -lo_hugeint
    integer, dimension(3) :: ssdim = -lo_hugeint
    integer :: nsample = -lo_hugeint
    logical :: U0 = .false.
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
    real(flyt), dimension(3) :: dumflytv

    ! basic info
    call cli%init(progname='phonon_dispersion_relations', &
                  authors=lo_author, &
                  version=lo_version, &
                  license=lo_licence, &
                  help='Usage: ', &
                  description='Calculate phonon dispersions and related quantities. &
                              &Per default, only the dispersions along a default path &
                              &will be calculated. Options are available for calculating &
                              &mode gruneisen parameters, phonon density of states projected &
                              &in a variety of ways, thermodynamic quantities and pure data dumps.', &
                  examples=["phonon_dispersion_relations                              ", &
                            "phonon_dispersion_relations --dos --qpoint_grid 24 24 24 "], &
                  epilog=new_line('a')//"...")

! Global things
    cli_unit
    cli_nq_on_path
    cli_readpath

! Dos things
    call cli%add(switch='--dos', &
                 help='Calculate the phonon DOS', &
                 required=.false., act='store_true', def='.false.', error=lo_status)
    if (lo_status .ne. 0) stop

    cli_qpoint_grid
    cli_meshtype
    cli_sigma
    cli_readqmesh

    call cli%add(switch='--integrationtype', switch_ab='-it', &
                 help='Type of integration for the phonon DOS. 1 is Gaussian, 2 adaptive Gaussian and 3 Tetrahedron.', &
                 required=.false., act='store', def='2', choices='1,2,3', error=lo_status)
    if (lo_status .ne. 0) stop

    call cli%add(switch='--dospoints', &
                 help='Number of points on the frequency axis of the phonon dos.', &
                 required=.false., act='store', def='400', error=lo_status)
    if (lo_status .ne. 0) stop

! Free energy things
    call cli%add(switch='--temperature',hidden=.true., &
                 help='Evaluate thermodynamic phonon properties at a single temperature.', &
                 required=.false., act='store', def='-1', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--temperature_range',hidden=.true., &
                 help='Evaluate thermodynamic phonon properties for a series of &
                      &temperatures, specify min, max and the number of points.', &
                 nargs='3', required=.false., act='store', def='-1 -1 -1', error=lo_status)
    if (lo_status .ne. 0) stop

! Third order things
    call cli%add(switch='--gruneisen', &
                 help='Use third order force constants to calculate mode Gruneisen parameters.', &
                 required=.false., act='store_true', def='.false.', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--dumpgrid', &
                 help='Write files with q-vectors, frequencies, eigenvectors and group velocities for a grid.', &
                 required=.false., act='store_true', def='.false.', error=lo_status)
    if (lo_status .ne. 0) stop

    ! Printing options
    call cli%add(switch='--pdf', switch_ab='-p', hidden=.true., &
                 help='Produce gnuplot_pdf output file for printing to pdf.', &
                 required=.false., act='store_true', def='.false.', error=lo_status)
    if (lo_status .ne. 0) stop

! Hidden stuff, not intented for the end-user
    cli_manpage
    cli_verbose

    call cli%add(switch='--alloydos', hidden=.true., help='', &
                 required=.false., act='store_true', def='.false.', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--alloy_ssdim', hidden=.true., help='', &
                 nargs='3', required=.false., act='store', def='-1 -1 -1', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--readiso', hidden=.true., help='', &
                 required=.false., act='store_true', def='.false.', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--notr', hidden=.true., help='', &
                 required=.false., act='store_true', def='.false.', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--sortbands', hidden=.true., help='', &
                 required=.false., act='store_true', def='.false.', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--dumpsecret', hidden=.true., &
                 help='Write file with everything one can deduce from a single point phonon calculation.', &
                 required=.false., act='store_true', def='.false.', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--U0', hidden=.true., &
                 help='Calculate some free energy differences.', &
                 required=.false., act='store_true', def='.false.', error=lo_status)
    if (lo_status .ne. 0) stop

! Unfolding things
    call cli%add(switch='--unfold', hidden=.true., &
                 help='Try to calculate unfolded phonon bandstructures.', &
                 required=.false., act='store_true', def='.false.', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--nunfoldenergy', hidden=.true., &
                 help='Number of points on the frequency axis of the unfolded spectrum', &
                 required=.false., act='store', def='2000', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--nunfoldsample', hidden=.true., &
                 help='Number of points on the frequency axis of the unfolded spectrum', &
                 required=.false., act='store', def='10', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--unfolddimensions', hidden=.true., &
                 help='Dimensions of supercell.', &
                 nargs='3', required=.false., act='store', def='3 3 3', error=lo_status)
    if (lo_status .ne. 0) stop
! Inelastic scattering spectra
    call cli%add(switch='--inelastic', hidden=.true., &
                 help='Calculate the inelastic scattering intensitiy.', &
                 required=.false., act='store_true', def='.false.', error=lo_status)
    if (lo_status .ne. 0) stop

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

! Now for the annoying step, all these things needs to be stuffed into my options type. Otherwise the
! code will get really messy from all the calls to cli, and passing them to functions will be annoying.
    call cli%get(switch='--unit', val=opts%enhet)
    call cli%get(switch='--nq_on_path', val=opts%nq)
    call cli%get(switch='--readpath', val=opts%readpathfromfile)
    call cli%get(switch='--dos', val=opts%dos)
    call cli%get(switch='--qpoint_grid', val=opts%qgrid)
    call cli%get(switch='--meshtype', val=opts%meshtype)
    call cli%get(switch='--sigma', val=opts%sigma)
    call cli%get(switch='--readqmesh', val=opts%readqmesh)
    call cli%get(switch='--integrationtype', val=opts%dosintegrationtype)
    call cli%get(switch='--dospoints', val=opts%dospoints)
    call cli%get(switch='--gruneisen', val=opts%gruneisen)
    call cli%get(switch='--sortbands', val=opts%sortbands)
    call cli%get(switch='--temperature', val=opts%temperature, error=lo_status)
    call cli%get(switch='--temperature_range', val=dumflytv)
    opts%trangemin = dumflytv(1)
    opts%trangemax = dumflytv(2)
    opts%trangenpts = int(anint(dumflytv(3)))
! if --temperature is specified, override the range
    if (opts%temperature .gt. -lo_tiny) then
        opts%trangemin = opts%temperature
        opts%trangemax = opts%temperature
        opts%trangenpts = 1
    end if

    ! plot
    call cli%get(switch='--pdf', val=opts%pdfplot)

    call cli%get(switch='--verbose', val=dumlog)
    if (dumlog) then
        opts%verbosity = 2
    else
        opts%verbosity = 0
    end if
    call cli%get(switch='--readiso', val=opts%readiso)
    call cli%get(switch='--dumpgrid', val=opts%dumpgrid)
    call cli%get(switch='--dumpsecret', val=opts%dumpsecret)

    call cli%get(switch='--unfold', val=opts%unfold)
    call cli%get(switch='--nunfoldenergy', val=opts%nenergy)
    call cli%get(switch='--nunfoldsample', val=opts%nsample)
    call cli%get(switch='--unfolddimensions', val=opts%ssdim)

    call cli%get(switch='--inelastic', val=opts%inelastic)
    call cli%get(switch='--U0', val=opts%U0)

! should the full mesh be calculated?
    opts%fullmesh = .false.
    if (opts%dos) opts%fullmesh = .true.
    if (opts%dumpgrid) opts%fullmesh = .true.
    if (opts%trangenpts .gt. 0) opts%fullmesh = .true.
    if (opts%inelastic) opts%fullmesh = .true.
    if (opts%dumpsecret) then
        opts%fullmesh = .true.
        opts%dos = .true.
    end if

    call cli%get(switch='--notr', val=dumlog)
    opts%timereversal = .not. dumlog

! Stupid switch to make secret dumping faster
    if (opts%dumpsecret) opts%nq = 3

end subroutine

end module
