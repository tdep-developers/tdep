#include "precompilerdefinitions"
module options
use konstanter, only: flyt, lo_status, lo_author, lo_version, lo_licence, lo_hugeint
use flap, only: command_line_interface
implicit none

private
public :: lo_opts

type lo_opts
    !> how much to talk
    integer :: verbosity = -lo_hugeint
    !> use every N timesteps
    integer :: stride = -lo_hugeint
    !> is it a constant pressure simulation?
    logical :: npt = .false.
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
    integer :: i
    logical :: dumlog
    real(flyt), dimension(3) :: dumflytv

    call cli%init(progname='average_structure', &
                  authors=lo_author, &
                  version=lo_version, &
                  license=lo_licence, &
                  help='Usage: ', &
                  description='Utility to fit equations of state to volume-energy data.', &
                  examples=["eosfit"], &
                  epilog=new_line('a')//"...")

    cli_manpage
    cli_verbose
    call cli%add(switch='--stride', switch_ab='-s', &
                 help='Number of dimenions for the fit', &
                 required=.false., act='store', def='1', error=lo_status)
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
    opts%verbosity = 0
    call cli%get(switch='--verbose', val=dumlog)
    if (dumlog) opts%verbosity = 2

    ! Parse the rest
    call cli%get(switch='--stride', val=opts%stride, error=lo_status)
end subroutine

end module
