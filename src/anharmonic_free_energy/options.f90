#include "precompilerdefinitions"
module options
use konstanter, only: r8,lo_author,lo_version,lo_licence,lo_tiny,lo_status,lo_exitcode_baddim,lo_exitcode_param,&
                      lo_frequency_Hartree_to_THz,lo_frequency_Hartree_to_icm,lo_frequency_Hartree_to_meV
use gottochblandat, only: lo_stop_gracefully
use flap, only: command_line_interface
implicit none
private
public :: lo_opts

type lo_opts
    character(len=3) :: enhet
    integer, dimension(3) :: qgrid
!    integer, dimension(3) :: qgrid_support
    integer :: trangenpts
    real(r8) :: temperature
    real(r8) :: trangemin
    real(r8) :: trangemax
    real(r8) :: sigma
    real(r8) :: reference_energy
    !
    integer :: integrationtype
    logical :: readiso
    integer :: verbosity
    integer :: nq_on_path
    logical :: readpathfromfile
    logical :: readqmesh
    integer :: meshtype
    !
    logical :: fullgrid
    logical :: qpointpath
    !
    real(r8) :: unitfactor
    character(len=4) :: unitname
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
    real(r8), dimension(3) :: dumr8v
    integer :: errctr

    ! basic info
    call cli%init(progname    = 'anharmonic_free_energy',&
              authors     = lo_author,&
              version     = lo_version,&
              license     = lo_licence,&
              help        = 'Usage: ',&
              description = 'End world hunger!',&
              examples    = ["mpirun anharmonic_free_energy"],&
              epilog      = new_line('a')//"...")

    ! Specify some options
    !cli_unit
    cli_temperature
    cli_qpoint_grid
    call cli%add(switch='--integrationtype',switch_ab='-it',&
        help='Type of integration. 1 is Gaussian, 2 adaptive Gaussian and 3 Tetrahedron.',&
        required=.false.,act='store',def='2',choices='1,2,3',error=lo_status)
        if ( lo_status .ne. 0 ) stop
    cli_sigma
    ! call cli%add(switch='--path',&
    !     help='Calculate the anharmonic free energy on a path through the BZ.',&
    !     required=.false.,act='store_true',def='.false.',error=lo_status)
    !     if ( lo_status .ne. 0 ) stop
    !cli_readpath
    !cli_nq_on_path
    ! call cli%add(switch='--support_qpoint_grid',switch_ab='-sqg',&
    !     help='Interpolate to a (preferrably) denser q-mesh when calculating the DOS.',&
    !     nargs='3',required=.false.,act='store',def='-1 -1 -1',error=lo_status)
    !     if ( lo_status .ne. 0 ) stop
    call cli%add(switch='--temperature_range',&
        help='Evaluate thermodynamic phonon properties for a series of temperatures, specify min, max and the number of points.',&
        nargs='3',required=.false.,act='store',def='-1 -1 -1',error=lo_status)
        if ( lo_status .ne. 0 ) stop
    cli_readiso
    cli_manpage
    cli_verbose
    cli_meshtype
    cli_readqmesh
    !call cli%add(switch='--Uref',&
    !    help='Baseline energy to be substracted when calculating U0',&
    !    required=.false.,act='store',def='0',error=lo_status)
    !    if ( lo_status .ne. 0 ) stop

    ! actually parse it
    errctr=0
    call cli%parse(error=lo_status)
    if ( lo_status .ne. 0 ) stop
    ! Should the manpage be generated? In that case, no point to go further than here.
    dumlog=.false.
    call cli%get(switch='--manpage',val=dumlog,error=lo_status); errctr=errctr+lo_status
    if ( dumlog ) then
        call cli%save_man_page(trim(cli%progname)//'.1')
        call cli%save_usage_to_markdown(trim(cli%progname)//'.md')
        write(*,*) 'Wrote manpage for "'//trim(cli%progname)//'"'
        stop
    endif
    call cli%get(switch='--verbose',val=dumlog,error=lo_status); errctr=errctr+lo_status; if ( dumlog ) opts%verbosity=2

    ! get real options
    !call cli%get(switch='--unit',               val=opts%enhet,            error=lo_status); errctr=errctr+lo_status
    !call cli%get(switch='--nq_on_path',         val=opts%nq_on_path,       error=lo_status); errctr=errctr+lo_status
    !call cli%get(switch='--path',               val=opts%qpointpath,       error=lo_status); errctr=errctr+lo_status
    !call cli%get(switch='--readpath',           val=opts%readpathfromfile, error=lo_status); errctr=errctr+lo_status
    call cli%get(switch='--qpoint_grid',        val=opts%qgrid,            error=lo_status); errctr=errctr+lo_status
    !call cli%get(switch='--support_qpoint_grid',val=opts%qgrid_support,    error=lo_status); errctr=errctr+lo_status
    call cli%get(switch='--readqmesh',          val=opts%readqmesh,        error=lo_status); errctr=errctr+lo_status
    call cli%get(switch='--meshtype',           val=opts%meshtype,         error=lo_status); errctr=errctr+lo_status
    call cli%get(switch='--sigma',              val=opts%sigma,            error=lo_status); errctr=errctr+lo_status
    call cli%get(switch='--integrationtype',    val=opts%integrationtype,  error=lo_status); errctr=errctr+lo_status
    call cli%get(switch='--readiso',            val=opts%readiso,          error=lo_status); errctr=errctr+lo_status
    call cli%get(switch='--temperature',        val=opts%temperature,      error=lo_status); errctr=errctr+lo_status
    call cli%get(switch='--temperature_range',  val=dumr8v,                error=lo_status); errctr=errctr+lo_status
    !call cli%get(switch='-U0',                  val=opts%U0,               error=lo_status); errctr=errctr+lo_status
    call cli%get(switch='--Uref',               val=opts%reference_energy, error=lo_status); errctr=errctr+lo_status

    if ( errctr .ne. 0 ) call lo_stop_gracefully(['Failed parsing the command line options'],lo_exitcode_baddim)

    opts%trangemin=dumr8v(1)
    opts%trangemax=dumr8v(2)
    opts%trangenpts=int(anint(dumr8v(3)))
    ! if --temperature is specified, override the range
    if ( opts%temperature .gt. -lo_tiny ) then
        opts%trangemin=opts%temperature
        opts%trangemax=opts%temperature
        opts%trangenpts=1
    endif

    ! figure out what to do, and try to resolve strange conflicts
    opts%fullgrid=.true.
    if ( opts%qpointpath ) opts%fullgrid=.false.

    ! need a support grid if not explicitly specified
    ! if ( opts%fullgrid ) then
    !     if ( opts%qgrid_support(1) .eq. -1 ) then
    !         opts%qgrid_support=opts%qgrid
    !     endif
    ! endif

    ! can keep the unit factor here, quite convenient
    select case(opts%enhet)
    case('thz')
        opts%unitfactor=lo_frequency_Hartree_to_THz
        opts%unitname='Thz'
    case('icm')
        opts%unitfactor=lo_frequency_Hartree_to_icm
        opts%unitname='cm-1'
    case('mev')
        opts%unitfactor=lo_frequency_Hartree_to_mev
        opts%unitname='meV'
    case default
        call lo_stop_gracefully(["Unknown unit, try 'thz', 'mev' or 'icm'"],lo_exitcode_param)
    end select

end subroutine

end module
