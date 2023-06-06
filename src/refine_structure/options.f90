#include "precompilerdefinitions"
module options
use konstanter, only: r8, lo_status, lo_author, lo_version, lo_licence, lo_huge
use flap, only: command_line_interface
implicit none

private
public :: lo_opts

type lo_opts
    !> tolerance
    real(r8) :: tolerance_lattice = -lo_huge
    real(r8) :: tolerance_internal = -lo_huge
    !> unitcell filename
    character(len=2000) :: unitcell_filename
    !> supercell filename
    character(len=2000) :: supercell_filename
    !> prototype unitcell filename
    character(len=2000) :: prototype_unitcell
    !> how much to talk
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

    logical :: dumlog

    call cli%init(progname='refine_structure', &
                  authors=lo_author, &
                  version=lo_version, &
                  license=lo_licence, &
                  help='Usage: ', &
                  description='Small utility to ensure that the input structure satisfies all symmetries to high precision.', &
                  examples=["refine_structure"], &
                  epilog=new_line('a')//"...")
    cli_manpage
    cli_verbose
    call cli%add(switch='--tolerance_lattice', switch_ab='-tl', &
                 help='Tolerance for the lattice.', &
                 required=.false., act='store', def='1E-5', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--tolerance_internal', switch_ab='-ti', &
                 help='Tolerance for internal degrees of freedom.', &
                 required=.false., act='store', def='1E-5', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--unitcell', switch_ab='-uc', &
                 help='Filename for the unitcell to refine.', &
                 required=.false., act='store', def='infile.ucposcar', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--supercell', switch_ab='-ss', &
                 help='Filename for supercell to refine.', &
                 required=.false., act='store', def='none', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--prototype', switch_ab='-pf', &
                 help='Prototype unitcell where you know the symmetry is correct.', &
                 required=.false., act='store', def='none', error=lo_status)
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
    call cli%get(switch='--tolerance_internal', val=opts%tolerance_internal, error=lo_status)
    call cli%get(switch='--tolerance_lattice', val=opts%tolerance_lattice, error=lo_status)
    call cli%get(switch='--unitcell', val=opts%unitcell_filename, error=lo_status)
    call cli%get(switch='--supercell', val=opts%supercell_filename, error=lo_status)
    call cli%get(switch='--prototype', val=opts%prototype_unitcell, error=lo_status)

end subroutine

end module
