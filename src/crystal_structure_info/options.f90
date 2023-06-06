#include "precompilerdefinitions"
module options
use konstanter, only: lo_author, lo_version, lo_licence, lo_status
use flap, only: command_line_interface
implicit none
private
public :: lo_opts

type lo_opts
    !> Print out the symmetry operations
    logical :: printsymmetry
    !> How much to talk
    integer :: verbosity
    !> Consider time reversal symmetry
    logical :: timereversal
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
    !
    call cli%init(progname='crystal_structure_info', &
                  authors=lo_author, &
                  version=lo_version, &
                  license=lo_licence, &
                  help='Usage: ', &
                  description='This code serves as a diagnostic tool to check that symmetry heuristics are working as they should. The code prints which Bravais lattice was identified, which high symmetry points in the BZ are inequivalent, and so on. The Brillouin zone and its irreducible wedges are printed to files as polyhedra for manual inspection, and the symmetry operations of the lattice can be printed.', &
                  examples=["crystal_structure_info                ", &
                            "crystal_structure_info --printsymmetry"], &
                  epilog=new_line('a')//"...")
    !
    call cli%add(switch='--printsymmetry', &
                 help='Also prints the symmetry operations', &
                 required=.false., act='store_true', def='.false.', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--notr', hidden=.true., help='', &
                 required=.false., act='store_true', def='.false.', error=lo_status)
    if (lo_status .ne. 0) stop

    !
    cli_manpage
    cli_verbose
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
    ! verbosity?
    opts%verbosity = 0
    call cli%get(switch='--verbose', val=dumlog)
    if (dumlog) opts%verbosity = 2
    call cli%get(switch='--printsymmetry', val=opts%printsymmetry)
    call cli%get(switch='--notr', val=dumlog)
    opts%timereversal = .not. dumlog
    !
end subroutine

end module
