#include "precompilerdefinitions"
module options
!
use konstanter
use flap

type lo_opts
    !> q-grid density
    integer, dimension(3) :: qgrid
    !> what kind of qgrid
    integer :: meshtype
    !> read it from file instead
    logical :: readqpointsfromfile
    contains
        procedure :: parse
end type

contains

subroutine parse(opts)
    !> the options
    class(lo_opts), intent(out) :: opts

    ! the parser
    type(command_line_interface) :: cli

    ! basic stuff, for the help thingy
    call cli%init(progname    = 'dump_dynamical_matrices',&
              authors     = lo_author,&
              version     = lo_version,&
              license     = lo_licence,&
              help        = 'Usage: ',&
              description = 'Dump the dynamical matrices to file, either on a grid or at the points specified by a file.',&
              epilog      = new_line('a')//"...")

    cli_qpoint_grid
    cli_meshtype
    !
    call cli%add(switch='--readqpoints',&
            help='Instead of generating a q-mesh, read it in fractional coordinates from a file called "infile.dynmatqpoints"',&
            required=.false.,act='store_true',def='.false.',error=lo_status)
            if ( lo_status .ne. 0 ) stop
    call cli%get(switch='--qpoint_grid',val=opts%qgrid)
    call cli%get(switch='--meshtype',val=opts%meshtype)    
    call cli%get(switch='--readqpoints',val=opts%readqpointsfromfile)    
end subroutine


end module
