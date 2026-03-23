#include "precompilerdefinitions"
module options
use konstanter, only: lo_status, lo_author, lo_version, lo_licence, lo_hugeint, lo_huge, flyt
use flap, only: command_line_interface
private
public :: lo_opts

type lo_opts
    !> how much to talk
    integer :: verbosity = -lo_hugeint
    !> grid dimensions
    integer, dimension(3) :: qgrid = -lo_hugeint
    !> read path from file
    logical :: readpathfromfile = .false.
    !> How many points should there be along the path?
    integer :: nqpath = -lo_hugeint
    !> what kind of mesh
    integer :: meshtype = -lo_hugeint
    !> refine wedge meshes
    logical :: refinemesh = .false.
    !> force tolerance for wedge mesh
    real(flyt) :: forcetolerance = -lo_huge
    !> size tolerance for wedge mesh
    real(flyt) :: sizetolerance = -lo_huge
    !> determine mesh density automatically
    logical :: autodens = .false.
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

    call cli%init(progname='generate_kpoints', &
                  authors=lo_author, &
                  version=lo_version, &
                  license=lo_licence, &
                  help='Usage: ', &
                  description='Utility to generate k/q-point meshes and paths in the Brillouin zone.', &
                  examples=["generate_kpoints -qg 16 16 16                ", &
                            "generate_kpoints --meshtype 3 -qg 30 30 30   "], &
                  epilog=new_line('a')//"...")

    cli_qpoint_grid
    cli_meshtype
    cli_nq_on_path
    cli_readpath
    cli_manpage
    cli_verbose

    call cli%add(switch='--refine', &
                 help='Refine wedge based mesh to get a more even distribution of tetrahedron angles/volumes/edges.', &
                 required=.false., act='store_true', def='.false.', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--forcetol', &
                 help='When refining the mesh, what is the max force on a vertex.', &
                 required=.false., act='store', def='1E-4', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--sizetol', &
                 help='When refining the mesh, what is the smallest tetrahedron that will be split.', &
                 required=.false., act='store', def='1.5', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--autodensity', &
                 help='Automatically pick the number of grid-points in each direction to get a uniform mesh.', &
                 required=.false., act='store_true', def='.false.', error=lo_status)
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

    ! Parse the rest
    call cli%get(switch='--qpoint_grid', val=opts%qgrid)
    call cli%get(switch='--nq_on_path', val=opts%nqpath)
    call cli%get(switch='--readpath', val=opts%readpathfromfile)
    call cli%get(switch='--meshtype', val=opts%meshtype)
    call cli%get(switch='--refine', val=opts%refinemesh)
    call cli%get(switch='--forcetol', val=opts%forcetolerance)
    call cli%get(switch='--sizetol', val=opts%sizetolerance)
    call cli%get(switch='--autodensity', val=opts%autodens)

end subroutine

end module
