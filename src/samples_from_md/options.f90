#include "precompilerdefinitions"
module options
use konstanter, only: flyt, lo_author, lo_licence, lo_version, lo_status
use flap, only: command_line_interface
implicit none
private
public :: lo_opts

type lo_opts
    integer :: n
    integer :: maxiter
    real(flyt) :: temp
    integer :: output_format
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

    logical :: dumlog

    call cli%init(progname='samples_from_md', &
                  authors=lo_author, &
                  version=lo_version, &
                  license=lo_licence, &
                  help='Usage: ', &
                  description='Choose representative uncorrelated samples from an MD simulation. The samples are chosen to be approximately evenly spaced, and reproduce the average potential energy, average kinetic energy have the same standard deviation of potential and kinetic energy.', &
                  examples=["samples_from_md -n 100"], &
                  epilog=new_line('a')//"...")

    call cli%add(switch='--nsamples', switch_ab='-n', &
                 help='Number of samples', &
                 required=.false., act='store', def='50', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--output_format', switch_ab='-of', &
                 help='Selects output format. 1 is VASP, 2 is Abinit, 4 is FHI-Aims, 5 is Siesta. Default 1.', &
                 required=.false., act='store', def='1', choices='1,2,4,5', error=lo_status)
    if (lo_status .ne. 0) stop

    cli_manpage
    dumlog = .false.
    call cli%get(switch='--manpage', val=dumlog)
    if (dumlog) then
        call cli%save_man_page(trim(cli%progname)//'.1')
        call cli%save_usage_to_markdown(trim(cli%progname)//'.md')
        write (*, *) 'Wrote manpage for "'//trim(cli%progname)//'"'
        stop
    end if
    ! actually parse it
    call cli%parse(error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%get(switch='--nsamples', val=opts%n)
    call cli%get(switch='--output_format', val=opts%output_format)
    opts%maxiter = 100000
    opts%temp = 1E-4_flyt
end subroutine

end module
