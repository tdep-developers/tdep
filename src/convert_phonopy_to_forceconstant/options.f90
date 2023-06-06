#include "precompilerdefinitions"
module options
use konstanter, only: lo_hugeint, lo_author, lo_version, lo_licence
use flap, only: command_line_interface

type lo_opts
    !> talk a lot?
    integer :: verbosity = -lo_hugeint
    !> filename for IFCs
    character(len=10000) :: fc_filename
    !> filename for BORN
    character(len=10000) :: born_filename
    !> unitcell filename
    character(len=10000) :: uc_filename
    !> supercell filename
    character(len=10000) :: ss_filename
    !> truncate forceconstant to get better symmetry
    logical :: truncate = .false.
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

    call cli%init(progname='convert_phonopy_to_forceconstant', &
                  authors=lo_author, &
                  version=lo_version, &
                  license=lo_licence, &
                  help='Usage: ', &
                  description='Generates a forceconstant file from phonopy data.', &
                  epilog=new_line('a')//"...")

    cli_verbose
    cli_manpage

    call cli%add(switch='-fc', &
                 help='Specify the filename for the phonopy FORCE_CONSTANTS file.', &
                 required=.false., act='store', def='FORCE_CONSTANTS', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='-fb', &
                 help='Specify the filename for the phonopy BORN file.', &
                 required=.false., act='store', def='BORN', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='-fuc', &
                 help='Specify the filename for unitcell file.', &
                 required=.false., act='store', def='infile.ucposcar', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='-fss', &
                 help='Specify the filename for the supercell file.', &
                 required=.false., act='store', def='infile.ssposcar', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--truncate', &
                 help='Truncate the realspace IFCs.', &
                 required=.false., act='store_true', def='.false.', error=lo_status)
    if (lo_status .ne. 0) stop

    call cli%parse(error=lo_status)
    if (lo_status .ne. 0) stop

    call cli%get(switch='--manpage', val=dumlog)
    if (dumlog) then
        call cli%save_man_page(trim(cli%progname)//'.1')
        call cli%save_usage_to_markdown(trim(cli%progname)//'.md')
        write (*, *) 'Wrote manpage for "'//trim(cli%progname)//'"'
        stop
    end if
    call cli%get(switch='--verbose', val=dumlog)
    if (dumlog) opts%verbosity = 2
    call cli%get(switch='-fc', val=opts%fc_filename)
    call cli%get(switch='-fb', val=opts%born_filename)
    call cli%get(switch='-fuc', val=opts%uc_filename)
    call cli%get(switch='-fss', val=opts%ss_filename)
    call cli%get(switch='--truncate', val=opts%truncate)

end subroutine

end module
