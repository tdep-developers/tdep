#include "precompilerdefinitions"
module options
use konstanter, only: flyt, lo_status, lo_author, lo_version, lo_licence, lo_A_to_Bohr
use flap, only: command_line_interface
implicit none
private
public :: lo_opts

type lo_opts
    real(flyt) :: cutoff2
    real(flyt) :: cutoff3
    real(flyt) :: cutoff4
    logical :: loto
    logical :: collmagnetism
    logical :: noncollmagnetism
    ! supercell dimensions
    integer, dimension(3) :: ssdim
    integer, dimension(3, 3) :: ndssdim
    ! magnetic stuff
    integer :: magnconf
    integer :: magnbin
    ! sqs stuff
    integer :: nsqs
    integer :: verbosity
    integer :: outputformat
    ! automagic stuff
    integer :: desired_na
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
    integer, dimension(9) :: dumr

    call cli%init(progname='generate_structure', &
                  authors=lo_author, &
                  version=lo_version, &
                  license=lo_licence, &
                  help='Usage: ', &
                  description='Builds supercells, diagonal and non-diagonal. Also has the capability to find the optimal supercells for a given lattice, very handy when you have complicated structures.', &
                  examples=["generate_structure -dim 4 3 5",&
                            "generate_structure -na 200   "], &
                  epilog=new_line('a')//"...")
    !
    call cli%add(switch='--dimensions', switch_ab='-d', &
                 help='Dimensions of supercell.', &
                 nargs='3', required=.false., act='store', def='5 5 5', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--nondiagonal_dimensions', switch_ab='-nd', hidden=.false., &
                 help='Non-diagonal dimensions of supercell.', &
                 nargs='9', required=.false., act='store', def='0 0 0 0 0 0 0 0 0', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--secondorder_cutoff', switch_ab='-rc2', hidden=.true., &
                 help='Cutoff for the second order force constants', &
                 required=.false., act='store', def='5.0', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--output_format', switch_ab='-of', hidden=.false., &
                 help='Output format. 1 is VASP, 2 Abinit, 3 LAMMPS, 4 FHI-Aims and 5 Siesta', &
                 required=.false., act='store', def='1', choices='1,2,3,4,5', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--thirdorder_cutoff', switch_ab='-rc3', hidden=.true., &
                 help='Cutoff for the third order force constants', &
                 required=.false., act='store', def='-1', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--fourthorder_cutoff', switch_ab='-rc4', hidden=.true., &
                 help='Cutoff for the fourth order force constants', &
                 required=.false., act='store', def='-1', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--collmagnetism', hidden=.true., &
                 help='Take collinear magnetism into account for the symmetry analysis. Requires a minimal input file ("infile.collmagnetism"), example:'//new_line('a')//new_line('a')// &
                 '        1.13'//new_line('a')// &
                 '        -1.13'//new_line('a')// &
                 '        2.23'//new_line('a')// &
                 '        -2.23'//new_line('a')//new_line('a')// &
                 '    one line for each atom in "infile.ucposcar" denoting the magnetic moment. The values are only used for comparison with each other, to ensure you can not transform an atom with moment 1 to an atom with moment -1.', &
                 required=.false., act='store_true', def='.false.', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--noncollmagnetism', hidden=.true., &
                 help='Take non-collinear magnetism into account. As above, requires a minimal input file ("infile.noncollmagnetism"), example:  '//new_line('a')//new_line('a')// &
                 '        1.0 -1.0  0.0  '//new_line('a')// &
                 '        0.0  1.0 -1.0  '//new_line('a')// &
                 '        2.3  2.3  0.0  '//new_line('a')// &
                 '        2.3 -2.3  0.0  '//new_line('a')//new_line('a')// &
                 '    one line for each atom in "infile.ucposcar" denoting the magnetic moment. As above, the values are only used for comparisons.', &
                 required=.false., act='store_true', def='.false.', error=lo_status)
    if (lo_status .ne. 0) stop
!    call cli%add(switch='--readinitial',hidden=.true.,&
!            help='Read the initial alloy distribution, and also output the relaxation history towards the completely random configuration.',&
!            required=.false.,act='store_true',def='.false.',error=lo_status)
!            if ( lo_status .ne. 0 ) stop
    call cli%add(switch='--magnetic_bins', switch_ab='-mb', hidden=.true., &
                 help='Number of bins in correlation functions, from order to disorder.', &
                 required=.false., act='store', def='5', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--magnetic_configurations', switch_ab='-mc', hidden=.true., &
                 help='Number of magnetic configurations per bin', &
                 required=.false., act='store', def='10', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--desired_na', switch_ab='-na', &
                 help='Desired number of atoms in supercell. Will try to choose a cell as cubic as possible.', &
                 required=.false., act='store', def='-1', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--nsqs', hidden=.true., &
                 help='Number SQS structures to generate', &
                 required=.false., act='store', def='5', error=lo_status)
    if (lo_status .ne. 0) stop
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
    ! verbose output?
    opts%verbosity = 0
    call cli%get(switch='--verbose', val=dumlog)
    if (dumlog) opts%verbosity = 1

    ! Parse the rest of stuff
    call cli%get(switch='-rc2', val=opts%cutoff2)
    call cli%get(switch='-rc3', val=opts%cutoff3)
    call cli%get(switch='-rc4', val=opts%cutoff4)
    call cli%get(switch='--collmagnetism', val=opts%collmagnetism)
    call cli%get(switch='--noncollmagnetism', val=opts%noncollmagnetism)
    call cli%get(switch='--dimensions', val=opts%ssdim)
    call cli%get(switch='--nondiagonal_dimensions', val=dumr)
    opts%ndssdim = reshape(dumr, [3, 3])
    call cli%get(switch='--output_format', val=opts%outputformat)

    call cli%get(switch='-mb', val=opts%magnbin)
    call cli%get(switch='-mc', val=opts%magnconf)
    call cli%get(switch='-na', val=opts%desired_na)
    call cli%get(switch='--nsqs', val=opts%nsqs)

    ! get input to atomic units right away
    opts%cutoff2 = opts%cutoff2*lo_A_to_Bohr
    opts%cutoff3 = opts%cutoff3*lo_A_to_Bohr
    opts%cutoff4 = opts%cutoff4*lo_A_to_Bohr
end subroutine

end module
