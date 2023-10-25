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
    integer, dimension(3) :: qgrid
    integer :: integrationtype
    logical :: readiso
    integer :: verbosity
    logical :: readqmesh
    integer :: meshtype
    logical :: quantum = .false.
    logical :: stochastic = .false.
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
              description = 'Calculates the anharmonic free energy.',&
              examples    = ["mpirun anharmonic_free_energy"],&
              epilog      = new_line('a')//"...")

    ! Specify some options
    cli_qpoint_grid
    call cli%add(switch='--integrationtype',switch_ab='-it',&
        help='Type of integration. 1 is Gaussian, 2 adaptive Gaussian and 3 Tetrahedron.',&
        required=.false.,act='store',def='2',choices='1,2,3',error=lo_status)
        if ( lo_status .ne. 0 ) stop
    call cli%add(switch='--quantum', &
        help='Use Bose-Einstein occupations to compute the free energy', &
        required=.false., act='store_true',def='.false.',error=lo_status)
        if ( lo_status .ne. 0 ) stop
    call cli%add(switch='--stochastic', &
        help='Add second order cumulant contribution to the free energy with a minus sign for self-consistent sampling', &
        required=.false., act='store_true',def='.false.',error=lo_status)
        if ( lo_status .ne. 0 ) stop
    cli_readiso
    cli_manpage
    cli_verbose
    cli_meshtype
    cli_readqmesh

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
    call cli%get(switch='--qpoint_grid',        val=opts%qgrid,            error=lo_status); errctr=errctr+lo_status
    call cli%get(switch='--readqmesh',          val=opts%readqmesh,        error=lo_status); errctr=errctr+lo_status
    call cli%get(switch='--meshtype',           val=opts%meshtype,         error=lo_status); errctr=errctr+lo_status
    call cli%get(switch='--integrationtype',    val=opts%integrationtype,  error=lo_status); errctr=errctr+lo_status
    call cli%get(switch='--readiso',            val=opts%readiso,          error=lo_status); errctr=errctr+lo_status
    call cli%get(switch='--quantum',            val=opts%quantum,          error=lo_status); errctr=errctr+lo_status
    call cli%get(switch='--stochastic',         val=opts%stochastic,       error=lo_status); errctr=errctr+lo_status

    if ( errctr .ne. 0 ) call lo_stop_gracefully(['Failed parsing the command line options'],lo_exitcode_baddim)

end subroutine

end module
