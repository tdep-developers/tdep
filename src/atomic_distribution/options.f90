#include "precompilerdefinitions"
module options
use konstanter, only: flyt, lo_status, lo_author, lo_licence, lo_version, lo_A_to_Bohr, lo_hugeint, lo_huge
use flap, only: command_line_interface
implicit none
private
public :: lo_opts

type lo_opts
    real(flyt) :: cutoff = -lo_huge
    integer :: nbin = -lo_hugeint
    integer :: bintype = -lo_hugeint
    integer :: stride = -lo_hugeint
    logical :: transform = .false.
    integer :: verbosity = -lo_hugeint
    logical :: promise_no_diffusion = .false.
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

    call cli%init(progname='atomic_distribution', &
                  authors=lo_author, &
                  version=lo_version, &
                  license=lo_licence, &
                  help='Usage: ', &
                  description='Calculates properties of the atomic distribution from molecular dynamics, such as mean square displacement, pair distribution function, vector distribution functions and probability densities. Useful for analysing simulations close to instabilities/phase transitions to have some idea where the atoms are.', &
                  examples=["atomic_distribution --cutoff 4.3              ", &
                            "atomic_distribution --cutoff 4.3 --notransform"], &
                  epilog=new_line('a')//"...")

    call cli%add(switch='--cutoff', switch_ab='-r', &
                 help='Consider pairs up to this distance, in A.', &
                 required=.false., act='store', def='5.0', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--nbins', switch_ab='-n', &
                 help='Number of bins in each dimension.', &
                 required=.false., act='store', def='300', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--notransform', hidden=.false., &
                 help='Do no rotate the coordinate systems of the vector distribution. By default, the coordinate system is aligned with the positive x-direction in the direction of the bond.', &
                 required=.false., act='store_true', def='.false.', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--bintype', hidden=.false., &
                 help='Select the binning type for the vector distribution. 1 is straight binning (fastest), 2 is binning with a Gaussian, 3 is binning with a gaussian, but without the subpixel resolution.', &
                 required=.false., act='store', def='1', choices='1,2,3', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--stride', switch_ab='-s', &
                 help='Use every N configuration instead of all.', &
                 required=.false., act='store', def='1', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--nodiffusion', hidden=.false., &
                 help='Treat the system as though it is crystalline and not diffusing, not matter what.', &
                 required=.false., act='store_true', def='.false.', error=lo_status)
    if (lo_status .ne. 0) stop
    cli_manpage
    cli_verbose

    dumlog = .false.
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

    ! actually parse it
    call cli%parse(error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%get(switch='--cutoff', val=opts%cutoff)
    call cli%get(switch='--nbins', val=opts%nbin)
    call cli%get(switch='--notransform', val=dumlog)
    opts%transform = .not. dumlog
    call cli%get(switch='--bintype', val=opts%bintype)
    call cli%get(switch='--stride', val=opts%stride)
    call cli%get(switch='--nodiffusion', val=opts%promise_no_diffusion)

    ! convert to atomic units
    opts%cutoff = opts%cutoff*lo_A_to_Bohr
end subroutine

end module
