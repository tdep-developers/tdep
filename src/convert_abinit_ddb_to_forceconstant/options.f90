#include "precompilerdefinitions"
module options
use konstanter, only: flyt, lo_hugeint, lo_author, lo_version, lo_licence, lo_status, lo_huge, lo_A_to_bohr
use flap, only: command_line_interface
implicit none

private
public :: lo_opts

type lo_opts
    !> talk a lot
    integer :: verbosity = -lo_hugeint
    !> DDB filename(s)
    character(len=1000), dimension(:), allocatable :: filename
    !> q-mesh for the DFPT
    integer, dimension(3) :: qgrid = -lo_hugeint
    !> truncate the IFCs
    logical :: truncate = .false.
    !> enforce correction of symmetries
    logical :: enforcesym = .false.
    !> manually turn off polar corrections
    logical :: forcenopolar = .false.
    !> force a cutoff
    real(flyt) :: cutoff = -lo_huge
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

    call cli%init(progname='convert_abinit_ddb_to_forceconstant', &
                  authors=lo_author, &
                  version=lo_version, &
                  license=lo_licence, &
                  help='Usage: ', &
                  description='Generates "outfile.ucposcar" and "outfile.forceconstant" from abinit DDB files.', &
                  epilog=new_line('a')//"...")

    cli_verbose
    cli_manpage

    call cli%add(switch='--qpoint_grid', switch_ab='-qg', &
                 help='Phonon q-point mesh used in Abinit.', &
                 nargs='3', required=.false., act='store', def='-1 -1 -1', error=lo_status); 
    call cli%add(switch='--files', switch_ab='-f', &
                 help='DDB filename.', &
                 required=.false., act='store', def='', nargs='+', error=lo_status)
    call cli%add(switch='--truncate', &
                 help='Truncate the realspace IFCs.', &
                 required=.false., act='store_true', def='.false.', error=lo_status)
    call cli%add(switch='--enforcesym', hidden=.true., &
                 help='Force all symmetries to be satisfied.', &
                 required=.false., act='store_true', def='.false.', error=lo_status)
    call cli%add(switch='--nopolar', hidden=.true., &
                 help='Manually turn off polar corrections.', &
                 required=.false., act='store_true', def='.false.', error=lo_status)
    call cli%add(switch='--cutoff', switch_ab='-rc', &
                 help='Force a cutoff on the realspace ifcs.', &
                 required=.false., act='store', def='-1', error=lo_status); 
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
    if (dumlog) then
        opts%verbosity = 2
    else
        opts%verbosity = 0
    end if
    call cli%get_varying(switch='--files', val=opts%filename, error=lo_status)
    call cli%get(switch='--qpoint_grid', val=opts%qgrid, error=lo_status)
    call cli%get(switch='--truncate', val=opts%truncate, error=lo_status)
    call cli%get(switch='--enforcesym', val=opts%enforcesym, error=lo_status)
    call cli%get(switch='--nopolar', val=opts%forcenopolar, error=lo_status)
    call cli%get(switch='--cutoff', val=opts%cutoff, error=lo_status)

    opts%cutoff = opts%cutoff*lo_A_to_bohr

end subroutine

end module
