#include "precompilerdefinitions"
module options
use konstanter, only: lo_status, lo_author, lo_version, lo_licence
use flap, only: command_line_interface
implicit none
private
public :: lo_opts

type lo_opts
    logical :: lammps
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

    call cli%init(progname='remap_forceconstant', &
                  authors=lo_author, &
                  version=lo_version, &
                  license=lo_licence, &
                  help='Usage: ', &
                  description='Takes a set of forceconstants and a unit cell and outputs them for a different unit cell. The lattices must be equivalent, and if the new lattice is rotated with respect to old the rotation matrix has to be specified.', &
                  examples=["remap_forceconstant"], &
                  epilog=new_line('a')//"...")

    call cli%add(switch='--lammps', hidden=.true., &
                 help='Generate output compatible with my LAMMPS thing.', &
                 required=.false., act='store_true', def='.false.', error=lo_status)
    if (lo_status .ne. 0) stop

    cli_manpage
    ! actually parse it
    call cli%parse(error=lo_status)
    if (lo_status .ne. 0) stop

    dumlog = .false.
    call cli%get(switch='--manpage', val=dumlog)
    if (dumlog) then
        call cli%save_man_page(trim(cli%progname)//'.1')
        call cli%save_usage_to_markdown(trim(cli%progname)//'.md')
        write (*, *) 'Wrote manpage for "'//trim(cli%progname)//'"'
        stop
    end if

    call cli%get(switch='--lammps', val=opts%lammps)

end subroutine

end module
