#include "precompilerdefinitions"
module options
use konstanter, only: r8, lo_author, lo_version, lo_licence, lo_status, lo_frequency_THz_to_hartree, lo_huge, lo_hugeint, lo_A_to_bohr
use flap, only: command_line_interface
implicit none
private
public :: lo_opts

type lo_opts
    real(r8) :: temperature = -lo_huge
    real(r8) :: mindist = -lo_huge
    real(r8) :: debye_temperature = -lo_huge
    real(r8) :: maximum_frequency = -lo_huge
    real(r8) :: curie_temperature = -lo_huge
    integer :: output_format = -lo_hugeint
    integer :: nconf = -lo_hugeint
    integer :: verbosity = -lo_hugeint
    logical :: zpm = .false.
    ! For fake magnetic things
    real(r8) :: exchange_J = -lo_huge
    real(r8) :: mean_moment = -lo_huge
    ! Dump a faked simulation for circle tests
    logical :: fakesim = .false.
    ! Generate semi-random configurations for dielectric things
    logical :: semirandom = .false.
    real(r8) :: dielcutoff2 = -lo_huge
    real(r8) :: dielcutoff3 = -lo_huge
    logical :: modes = .false.
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

    call cli%init(progname='canonical_configuration', &
                  authors=lo_author, &
                  version=lo_version, &
                  license=lo_licence, &
                  help='Usage: ', &
                  description='Use forceconstants or a Debye temperature to generate uncorrelated supercell configurations emulating a canonical ensemble. These configurations can be used to either start ensemble runs of molecular dynamics with negligible equilibration time, or be used to directly sample phase space.', &
                  examples=["canonical_configuration -n 10 -t 300                                 ", &
                            "canonical_configuration -n 300 -t 0 --quantum                        ", &
                            "canonical_configuration -n 20 -t 10 --quantum --debye_temperature 400"], &
                  epilog=new_line('a')//"...")

    ! real options
    call cli%add(switch='--temperature', switch_ab='-t', &
                 help='Temperature to emulate', &
                 required=.false., act='store', def='300', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--nconf', switch_ab='-n', &
                 help='Number of configurations to generate', &
                 required=.false., act='store', def='5', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--quantum', &
                 help='Use Bose-Einstein statistics instead of Maxwell-Boltzmann. ', &
                 help_markdown="That is, use \( \sqrt{\frac{\hbar (2n+1) }{2 m \omega}} \) as the mean normal mode amplitudes instead of the classical \( \frac{1}{\omega}\sqrt{\frac{k_BT}{m}} \)", &
                 required=.false., act='store_true', def='.false.', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--output_format', switch_ab='-of', &
                 help='Selects output format. 1 is VASP, 2 is Abinit, 3 is LAMMPS, 4 is FHI-Aims, 5 is Siesta. Default 1.', &
                 required=.false., act='store', def='1', choices='1,2,3,4,5', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--mindist', &
                 help='What is the smallest distance between two atoms allowed, in units of the nearest neighbour distance.', &
                 required=.false., act='store', def='-1', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--debye_temperature', switch_ab='-td', &
                 help='Generate forceconstants that match a Debye temperature, and build displacements according to these.', &
                 help_markdown=' See details below.', &
                 required=.false., act='store', def='-1', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--maximum_frequency', switch_ab='-mf', &
                 help='Generate forceconstants that match a maximum frequency (in THz), and build displacements according to these.', &
                 help_markdown=' See details below.', &
                 required=.false., act='store', def='-1', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--critical_temperature', switch_ab='-tc', hidden=.true., &
                 help='Generate magnetic exchange interactions that match a certain critical temperature.', &
                 help_markdown=' See details below.', &
                 required=.false., act='store', def='-1', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--mean_moment', hidden=.true., &
                 help='Set the mean magnetic moment.', &
                 help_markdown=' See details below.', &
                 required=.false., act='store', def='1', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--Jij', hidden=.true., &
                 help='Strength of magnetic nearest neighbour exchange interactions.', &
                 help_markdown=' See details below.', &
                 required=.false., act='store', def='0', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--fakedataset', &
                 help='Dump a fake simulation to use for debugging.', hidden=.true., &
                 required=.false., act='store_true', def='.false.', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--semirandom', &
                 help='Do not pick configurations completele at random.', hidden=.true., &
                 required=.false., act='store_true', def='.false.', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--dielectric_pair_cutoff', switch_ab='-dc2', hidden=.true., &
                 help='Cutoff for the pair dielectric interactions', &
                 required=.false., act='store', def='-1.0', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--dielectric_triplet_cutoff', switch_ab='-dc3', hidden=.true., &
                 help='Cutoff for the triplet dielectric interactions', &
                 required=.false., act='store', def='-1.0', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--modes', &
                 help='Print displacements for every individual mode.', hidden=.true., &
                 required=.false., act='store_true', def='.false.', error=lo_status)
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
    ! verbosity?
    opts%verbosity = 0
    call cli%get(switch='--verbose', val=dumlog)
    if (dumlog) opts%verbosity = 2

    call cli%get(switch='--temperature', val=opts%temperature)
    call cli%get(switch='--debye_temperature', val=opts%debye_temperature)
    call cli%get(switch='--maximum_frequency', val=opts%maximum_frequency)
    call cli%get(switch='--mindist', val=opts%mindist)
    call cli%get(switch='--nconf', val=opts%nconf)
    call cli%get(switch='--output_format', val=opts%output_format)
    call cli%get(switch='--quantum', val=opts%zpm)
    call cli%get(switch='--mean_moment', val=opts%mean_moment)
    call cli%get(switch='--Jij', val=opts%exchange_J)
    call cli%get(switch='--fakedataset', val=opts%fakesim)
    call cli%get(switch='--semirandom', val=opts%semirandom)
    call cli%get(switch='-dc2', val=opts%dielcutoff2)
    call cli%get(switch='-dc3', val=opts%dielcutoff3)
    call cli%get(switch='--modes', val=opts%modes)

    ! Convert input to atomic units right away
    opts%maximum_frequency = opts%maximum_frequency*lo_frequency_THz_to_hartree
    opts%dielcutoff2 = opts%dielcutoff2*lo_A_to_Bohr
    opts%dielcutoff3 = opts%dielcutoff3*lo_A_to_Bohr

end subroutine

end module
