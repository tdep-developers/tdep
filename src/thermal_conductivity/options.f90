#include "precompilerdefinitions"
module options
use konstanter, only: flyt, lo_status, lo_author, lo_version, lo_licence, lo_m_to_bohr, lo_hugeint
use flap, only: command_line_interface
implicit none
private
public :: lo_opts

type lo_opts
    integer, dimension(3) :: qgrid   !< the main q-grid
    integer, dimension(3) :: qg3ph   !< The grid for the threephonon integration
    integer, dimension(3) :: qg4ph   !< The grid for the fourphonon integration
    logical :: readqmesh             !< read q-grid from file
    real(flyt) :: temperature        !< temperature
    real(flyt) :: sigma              !< scaling factor for adaptive gaussian
    real(flyt) :: tau_boundary       !< add a constant as boundary scattering
    real(flyt) :: mfp_max            !< add a length as boundary scattering
    real(flyt) :: itertol            !< tolerance for the iterative solution
    integer :: itermaxsteps          !< Number of iteration for the Boltzmann equation
    logical :: classical             !< Use a classical formulation
    logical :: readiso               !< read isotope distribution from file
    logical :: thirdorder            !< use fourth order contribution
    logical :: fourthorder           !< use fourth order contribution
    logical :: isotopescattering     !< use isotope scattering
    integer :: integrationtype       !< adaptive or standard gaussian integration
    integer :: seed                  !< seed for the Monte-Carlo grid

    ! Debugging things
    logical :: timereversal
    logical :: qpsymmetry
    !
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
    !
    logical :: dumlog
    real(flyt) :: f0
    real(flyt), dimension(3) :: dumflytv
    integer :: i

    ! basic info
    call cli%init(progname='thermal_conductivity', &
                  authors=lo_author, &
                  version=lo_version, &
                  license=lo_licence, &
                  help='Usage: ', &
                  description='Calculates the lattice thermal conductivity in the&
                              & mode-coupling formalism, including collective and off-diagonal&
                              & contributions up to fourth-order interactions.',&
                  examples=["mpirun thermal_conductivity --temperature 300                        ", &
                            "mpirun thermal_conductivity --fourthorder  -qg 30 30 30 -qg4ph 4 4 4 "], &
                  epilog=new_line('a')//"...")
    ! real options
    call cli%add(switch='--readiso', &
                 help='Read the isotope distribution from `infile.isotopes`.', &
                 help_markdown='The format is specified [here](../page/files.html#infile.isotopes).', &
                 required=.false., act='store_true', def='.false.', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--integrationtype', switch_ab='-it', &
                 help='Type of integration for the phonon DOS. 1 is Gaussian, 2 adaptive Gaussian.', &
                 required=.false., act='store', def='2', choices='1,2,6', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--nothirdorder', &
                 help='Do not consider third order contributions to the scattering.', &
                 required=.false., act='store_true', def='.false.', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--fourthorder', &
                 help='Consider four-phonon contributions to the scattering.', &
                 required=.false., act='store_true', def='.false.', error=lo_status)
    if (lo_status .ne. 0) stop
    cli_qpoint_grid
    call cli%add(switch='--sigma', &
                 help='Global scaling factor for Gaussian/adaptive Gaussian smearing. The default is determined procedurally, and scaled by this number.', &
                 required=.false., act='store', def='1.0', error=lo_status)
    if (lo_status .ne. 0) stop
    cli_readqmesh

    call cli%add(switch='--temperature', &
                 help='Evaluate thermal conductivity at a single temperature.', &
                 required=.false., act='store', def='300', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--max_mfp', &
                 help='Add a limit on the mean free path as an approximation of domain size (in m).', &
                 required=.false., act='store', def='-1', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--iterative_tolerance', &
                 help='Tolerance for the iterative solution.', &
                 required=.false., act='store', def='1e-5', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--noisotope', &
                 help='Do not consider isotope scattering.', &
                 required=.false., act='store_true', def='.false.', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--classical', &
                 help='Use the classical limit for phonon occupation and heat capacity.', &
                 required=.false., act='store_true', def='.false.', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--iterative_maxsteps', &
                 help='Max number of iterations for the iterative solution.', &
                 required=.false., act='store', def='200', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--qpoint_grid3ph', switch_ab='-qg3ph', &
                 help='Dimension of the grid for the threephonon integration.', &
                 nargs='3', required=.false., act='store', def='-1 -1 -1', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--qpoint_grid4ph', switch_ab='-qg4ph', &
                 help='Dimension of the grid for the fourphonon integration.', &
                 nargs='3', required=.false., act='store', def='-1 -1 -1', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--seed', &
                 help='Positive integer to seed the random number generator for the Monte-Carlo grids.', &
                 required=.false., act='store', def='-1', error=lo_status)
    if (lo_status .ne. 0) stop

    ! hidden
    call cli%add(switch='--tau_boundary', hidden=.true., &
                 help='Add a constant boundary scattering term to the lifetimes.', &
                 required=.false., act='store', def='-1', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--notr', hidden=.true., help='', &
                 required=.false., act='store_true', def='.false.', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--nosym', hidden=.true., help='', &
                 required=.false., act='store_true', def='.false.', error=lo_status)
    if (lo_status .ne. 0) stop
    cli_manpage
    cli_verbose

    ! actually parse it
    call cli%parse(error=lo_status)
    if (lo_status .ne. 0) stop
    !
    ! Should the manpage be generated? In that case, no point to go further than here.
    !
    dumlog = .false.
    call cli%get(switch='--manpage', val=dumlog)
    if (dumlog) then
        call cli%save_man_page(trim(cli%progname)//'.1')
        call cli%save_usage_to_markdown(trim(cli%progname)//'.md')
        write (*, *) 'Wrote manpage for "'//trim(cli%progname)//'"'
        stop
    end if

    ! store things in the right place

    call cli%get(switch='--temperature', val=opts%temperature)
    call cli%get(switch='--qpoint_grid', val=opts%qgrid)
    call cli%get(switch='--qpoint_grid3ph', val=opts%qg3ph)
    call cli%get(switch='--qpoint_grid4ph', val=opts%qg4ph)
    call cli%get(switch='--iterative_maxsteps', val=opts%itermaxsteps)
    call cli%get(switch='--sigma', val=opts%sigma)
    call cli%get(switch='--tau_boundary', val=opts%tau_boundary)
    call cli%get(switch='--nothirdorder', val=dumlog)
    opts%thirdorder = .not. dumlog
    call cli%get(switch='--fourthorder', val=opts%fourthorder)
    if (opts%tau_boundary .gt. 0.0_flyt) opts%tau_boundary = 1E10_flyt
    call cli%get(switch='--readqmesh', val=opts%readqmesh)
    call cli%get(switch='--integrationtype', val=opts%integrationtype)
    call cli%get(switch='--readiso', val=opts%readiso)
    call cli%get(switch='--max_mfp', val=opts%mfp_max)
    call cli%get(switch='--iterative_tolerance', val=opts%itertol)
    call cli%get(switch='--classical', val=opts%classical)
    call cli%get(switch='--seed', val=opts%seed)
    ! stuff that's not really an option
    call cli%get(switch='--notr', val=dumlog)
    opts%timereversal = .not. dumlog
    call cli%get(switch='--nosym', val=dumlog)
    opts%qpsymmetry = .not. dumlog
    call cli%get(switch='--verbose', val=dumlog)
    if (dumlog) then
        opts%verbosity = 2
    else
        opts%verbosity = 0
    end if
    call cli%get(switch='--noisotope', val=dumlog)
    opts%isotopescattering = .not. dumlog

    ! Get things to atomic units
    opts%mfp_max = opts%mfp_max*lo_m_to_Bohr

    if (maxval(opts%qg4ph) .gt. 0 .and. .not. opts%fourthorder) then
        write(*, *) 'You have to enable fourthorder to use a fourth order Monte-Carlo grid, stopping calculation.'
        stop
    end if

    ! Set automatic values for Monte-Carlo grids
    if (opts%thirdorder) then
        do i = 1, 3
            if (opts%qg3ph(i) .lt. 0 .or. opts%qg3ph(i) .gt. opts%qgrid(i)) opts%qg3ph(i) = opts%qgrid(i)
        end do
    end if
    if (opts%fourthorder) then
        do i = 1, 3
            if (opts%qg4ph(i) .lt. 0 .or. opts%qg4ph(i) .gt. opts%qgrid(i)) opts%qg4ph(i) = opts%qgrid(i)
        end do
    end if

end subroutine

end module
