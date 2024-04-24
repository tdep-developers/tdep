#include "precompilerdefinitions"
module options
use konstanter, only: flyt, lo_status, lo_author, lo_licence, lo_version
use flap, only: command_line_interface
implicit none
private
public :: lo_opts

type lo_opts
    integer, dimension(3) :: qgrid
    real(flyt), dimension(3) :: qpoint
    real(flyt), dimension(2) :: thetaphi
    logical :: readiso
    integer :: verbosity
    logical :: thirdorder
    logical :: readqmesh
    logical :: povray
    logical :: timereversal
    integer, dimension(:, :), allocatable :: modespec
    integer :: nmodespec
    integer :: povrayquality
    logical :: intensities
    character(len=10) :: highsymmetrypoint
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
    integer :: i, j, l, n
    integer, dimension(:), allocatable :: dumi
    logical :: dumlog

    ! basic info
    call cli%init(progname='phasespace_surface', &
                  authors=lo_author, &
                  version=lo_version, &
                  license=lo_licence, &
                  help='Usage: ', &
                  description='Calculate the three-phonon scattering surface.', &
                  examples=["phasespace_surface --highsymmetrypoint X  ", &
                            "phasespace_surface --qpoint 0.1 0.2 0.3   "], &
                  epilog=new_line('a')//"...")
    !
    call cli%add(switch='--qpoint', &
                 help='Specify the q-point for which to calculate the scattering surface.', &
                 nargs='3', required=.false., act='store', def='0 0 0', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--highsymmetrypoint', &
                 help='Same as above, but you can specify the label of a high-symmetry point instead, e.g. "X" or "L".', &
                 required=.false., act='store', def='none', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--modespec', &
                 help='Specify which surfaces to look at. 6 integers per surface: band1 band2 band3 plus/minus color opacity. Modes go from 1 to number of bands. plus/minus is 0 for plus events 1 for minus. Color from 1-3, opacity from 0-100.', &
                 nargs='+', required=.false., act='store', def='1 2 3 0 1 100', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--povray', &
                 help='Write POV-Ray script for rendering.', &
                 required=.false., act='store_true', def='.false.', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--intensities', &
                 help='Calculate intensities (matrix elements) as well.', &
                 required=.false., act='store_true', def='.false.', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--pv_theta_phi', &
                 help='Povray specification, theta and phi angles (in degrees) deteriming camera location.', &
                 nargs='2', required=.false., act='store', def='20 30', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--pv_quality', &
                 help='Povray rendering quality, 1-9 where 1 is worst.', &
                 required=.false., act='store', def='9', error=lo_status)
    if (lo_status .ne. 0) stop
    cli_qpoint_grid
    cli_readqmesh
    cli_manpage
    cli_verbose

    ! actually parse it
    call cli%parse(error=lo_status)
    if (lo_status .ne. 0) stop
    ! Should the manpage be generated? In that case, no point to go further than here.
    dumlog = .false.
    call cli%get(switch='--manpage', val=dumlog)
    if (dumlog) then
        call cli%save_man_page(trim(cli%progname)//'.1')
        call cli%save_usage_to_markdown(trim(cli%progname)//'.md')
        write (*, *) 'Wrote manpage for "'//trim(cli%progname)//'"'
        stop
    end if
    call cli%get(switch='--verbose', val=dumlog); if (dumlog) opts%verbosity = 2
    !
    ! get real options
    call cli%get(switch='--qpoint_grid', val=opts%qgrid)
    call cli%get(switch='--qpoint', val=opts%qpoint)
    call cli%get(switch='--highsymmetrypoint', val=opts%highsymmetrypoint)
    call cli%get(switch='--readqmesh', val=opts%readqmesh)
    call cli%get(switch='--povray', val=opts%povray)
    call cli%get(switch='--pv_theta_phi', val=opts%thetaphi)
    call cli%get(switch='--pv_quality', val=opts%povrayquality)
    call cli%get_varying(switch='--modespec', val=dumi)
    call cli%get(switch='--intensities', val=opts%intensities)

    ! Clean up the mode specification thing right away
    n = size(dumi, 1) ! number of arguments
    if (mod(n, 6) .ne. 0) then
        write (*, *) 'Error: the mode specification must be done in multiples of 6 integers.'
        stop
    end if
    opts%nmodespec = n/6
    lo_allocate(opts%modespec(6, opts%nmodespec))
    opts%modespec = 0
    l = 0
    do i = 1, opts%nmodespec
    do j = 1, 6
        l = l + 1
        opts%modespec(j, i) = dumi(l)
    end do
    end do

    opts%timereversal = .true.
    opts%thirdorder = .false.
    !
end subroutine

end module
