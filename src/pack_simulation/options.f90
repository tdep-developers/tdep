#include "precompilerdefinitions"
module options
use konstanter, only: r8, lo_status, lo_author, lo_version, lo_licence, lo_huge, lo_hugeint
use flap, only: command_line_interface
implicit none

private
public :: lo_opts

type lo_opts
    !> how much to talk
    integer :: verbosity = -lo_hugeint
    !> skip a few steps
    integer :: stride = -lo_hugeint
    !> pick a few steps at random
    integer :: nrand = -lo_hugeint
    !> is it a variable cell simulation
    logical :: variable_lattice = .false.
    !> do we have magnetic moments?
    logical :: magnetic_moments = .false.
    !> is it a molecular dynamics simulation
    logical :: molecular_dynamics = .false.
    !> do we have born charges and dielectric constants?
    logical :: dielectric = .false.
    !> specify temperature explicitly
    real(r8) :: temperature_thermostat = -lo_huge
    !> tidy up the simulation
    logical :: make_tidy = .false.

    !> sensible dielectric tensor (diagonal part)
    real(r8) :: eps0 = -lo_huge
    !> sensible standard deviation for dielectric tensor
    real(r8) :: dev0 = -lo_huge
    !> how many sigma to keep?
    integer :: nsigma = -lo_hugeint
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
    call cli%init(progname='pack_simulation', &
                  authors=lo_author, &
                  version=lo_version, &
                  license=lo_licence, &
                  help='Usage: ', &
                  description='Small utility to pack a simulation to hdf5.', &
                  examples=["pack_simulation"], &
                  epilog=new_line('a')//"...")

    cli_manpage
    cli_verbose

    call cli%add(switch='--stride', switch_ab='-s', &
                 help='Pack every N configuration instead of all.', &
                 required=.false., act='store', def='1', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--nrand', &
                 help='Pack N random configurations instead of all.', &
                 required=.false., act='store', def='-1', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--temperature', &
                 help='Override the simulation temperature.', &
                 required=.false., act='store', def='-1', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--variable_cell', switch_ab='-npt', &
                 help='Make sure to store variable cell information.', &
                 required=.false., act='store_true', def='.false.', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--magnetic_moments', switch_ab='-mag', &
                 help='Make sure to store projected magnetic moments.', &
                 required=.false., act='store_true', def='.false.', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--dielectric', &
                 help='Make sure to store dielectric constant and Born charges.', &
                 required=.false., act='store_true', def='.false.', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--molecular_dynamics', switch_ab='-md', &
                 help='Make sure to specify that this is real molecular dynamics.', &
                 required=.false., act='store_true', def='.false.', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--notidy', &
                 help='Per default the simulation is cleaned, drift removed and so on. This switch skips that.', &
                 required=.false., act='store_true', def='.false.', error=lo_status)
    if (lo_status .ne. 0) stop

    call cli%add(switch='--eps0', hidden=.true., &
                 help='Sensible dielectric tensor.', &
                 required=.false., act='store', def='-1', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--dev0', hidden=.true., &
                 help='Sensible deviation of dielectric tensor.', &
                 required=.false., act='store', def='-1', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--nsigma', hidden=.true., &
                 help='Number of deviations to keep.', &
                 required=.false., act='store', def='-1', error=lo_status)
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

    call cli%get(switch='--stride', val=opts%stride)
    call cli%get(switch='--nrand', val=opts%nrand)
    call cli%get(switch='--temperature', val=opts%temperature_thermostat)
    call cli%get(switch='--variable_cell', val=opts%variable_lattice)
    call cli%get(switch='--magnetic_moments', val=opts%magnetic_moments)
    call cli%get(switch='--molecular_dynamics', val=opts%molecular_dynamics)
    call cli%get(switch='--notidy', val=dumlog)
    call cli%get(switch='--dielectric', val=opts%dielectric)
    call cli%get(switch='--eps0', val=opts%eps0)
    call cli%get(switch='--dev0', val=opts%dev0)
    call cli%get(switch='--nsigma', val=opts%nsigma)
    opts%make_tidy = .not. dumlog

end subroutine

end module
