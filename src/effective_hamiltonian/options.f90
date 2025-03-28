#include "precompilerdefinitions"
module options
use konstanter, only: r8, lo_hugeint, lo_huge, lo_author, lo_version, lo_licence, lo_tiny, lo_status, lo_exitcode_baddim, lo_exitcode_param  
use gottochblandat, only: lo_stop_gracefully
use flap, only: command_line_interface
implicit none
private
public :: lo_opts

type lo_opts
    integer :: verbosity = -lo_hugeint
    logical :: thirdorder = .false.
    logical :: fourthorder = .false.
    integer :: stride = 1
    integer :: nconf = -1
    logical :: quantum = .false.
    real(flyt) :: temperature = -lo_huge
contains
    procedure :: parse
end type

contains

subroutine parse(opts)
    !> the options
    class(lo_opts), intent(out) :: opts
    !> the helper parser
    type(command_line_interface) :: cli

    logical :: dumlog
    real(r8), dimension(3) :: dumr8v
    integer :: errctr

    ! basic info
    call cli%init(progname='effective_hamiltonian', &
                  authors=lo_author, &
                  version=lo_version, &
                  license=lo_licence, &
                  help='Usage: ', &
                  description='Calculates the TDEP effective Hamiltonian from force constants. By default only the potential energy component is calculated.', &
                  examples=["mpirun effective_hamiltonian --thirdorder --fourthorder"], &
                  epilog=new_line('a')//"...")

    ! Specify some options
    call cli%add(switch='--thirdorder', &
                 help='Compute third order anharmonic correction to the free energy', &
                 required=.false., act='store_true', def='.false.', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--fourthorder', &
                 help='Compute fourth order anharmonic correction to the free energy', &
                 required=.false., act='store_true', def='.false.', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--stride', switch_ab='s' &
                 help='Use every N configuration instead of all.', &
                 required=.false., act='store', def='1', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--nconf', switch_ab='-s' &
        help='Automatically generate configurations, as done by the canonical configurations command. If not passed, an infile.positions is required.', &
        required=.false., act='store', def='-1', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--quantum', &
        help='If using --nconf, generate quantum or classical configurations.', &
        required=.false., act='store_true', def='.false.', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--temperature', switch_ab='-t' &
        help='If using --nconf, generate configurations at this temperature.', &
        required=.false., act='store', def='-1', error=lo_status)
    if (lo_status .ne. 0) stop
    cli_manpage
    cli_verbose

    ! actually parse it
    errctr = 0
    call cli%parse(error=lo_status)
    if (lo_status .ne. 0) stop
    ! Should the manpage be generated? In that case, no point to go further than here.
    dumlog = .false.
    call cli%get(switch='--manpage', val=dumlog, error=lo_status); errctr = errctr + lo_status
    if (dumlog) then
        call cli%save_man_page(trim(cli%progname)//'.1')
        call cli%save_usage_to_markdown(trim(cli%progname)//'.md')
        write (*, *) 'Wrote manpage for "'//trim(cli%progname)//'"'
        stop
    end if
    call cli%get(switch='--verbose', val=dumlog, error=lo_status); errctr = errctr + lo_status; if (dumlog) opts%verbosity = 2

    ! get real options
    call cli%get(switch='--thirdorder', val=opts%thirdorder, error=lo_status); errctr = errctr + lo_status
    call cli%get(switch='--fourthorder', val=opts%fourthorder, error=lo_status); errctr = errctr + lo_status
    call cli%get(switch="--stride", val = opts%stride, error = lo_status); errctr = errctr + lo_status
    call cli%get(switch="--nconf", val = opts%nconf, error = lo_status); errctr = errctr + lo_status
    call cli%get(switch="--temperature", val = opts%temperature, error = lo_status); errctr = errctr + lo_status
    call cli%get(switch='--quantum', val = opts%quantum, error=lo_status); errctr = errctr + lo_status

    if (errctr .ne. 0) call lo_stop_gracefully(['Failed parsing the command line options'], lo_exitcode_baddim)

end subroutine

end module
