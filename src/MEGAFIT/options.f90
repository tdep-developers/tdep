#include "precompilerdefinitions"
module options
use konstanter, only: flyt,lo_status,lo_author,lo_version,lo_licence,lo_huge,lo_hugeint,lo_A_to_bohr
use flap, only: command_line_interface
implicit none

private
public :: lo_opts

type lo_opts
    ! cutoffs
    real(flyt) :: cutoff2=-lo_huge
    real(flyt) :: cutoff3=-lo_huge
    real(flyt) :: cutoff4=-lo_huge
    ! consider first order force constants
    logical :: firstorder=.false.
    ! use WZ cutoff for second order?
    integer, dimension(3) :: wzdimensions=-1
    ! use graph jumps for the cutoff?
    integer :: njump2=-lo_hugeint
    integer :: njump3=-lo_hugeint
    integer :: njump4=-lo_hugeint
    ! how much to talk
    integer :: verbosity=-lo_hugeint
    ! use rotational constraints on the forceconstant
    logical :: rotationalconstraints=.false.
    ! use Huang invariances
    logical :: huanginvariance=.false.
    ! use spacegroup symmetries
    logical :: spacegroup=.false.
    ! use transposition symmetries
    logical :: transposition=.false.
    ! enforce Hermitian character
    logical :: hermitian=.false.
    ! is this a polar material?
    logical :: polar=.false.
    ! what kind of polar correction?
    integer :: polarcorrectiontype=-lo_hugeint
    ! how should I fit the second order
    integer :: pairfittype=-lo_hugeint
    ! cutoff for magnetic pair interactions
    real(flyt) :: magcutoff2=-lo_hugeint
    ! Order of polynomial fits
    integer :: order=-lo_hugeint
    ! Evaluate the free energy on the grid-points
    logical :: evalenergy=.false.
    ! q-mesh
    integer, dimension(3) :: qgrid_harm=-lo_hugeint
    integer, dimension(3) :: qgrid_anharm=-lo_hugeint
    ! Temperature scaling thing
    real(flyt) :: temperature_scale=-lo_huge
    ! Evaluate quasiharmonic energies as well
    logical :: quasiharmonic=.false.
    ! Distance scale parameter
    real(flyt) :: distance_scale=-lo_huge
    ! Dump all the points on the grid
    logical :: dumpfullgrid=.false.
    ! Dump things on the input grid
    logical :: dumpinputgrid=.false.
    ! Calculate a bunch of diagnostics
    logical :: diagnostics=.false.
    ! cutoff for dielectric pair interactions
    real(flyt) :: dielcutoff2=-lo_huge
    ! cutoff for dielectric triplet interactions
    real(flyt) :: dielcutoff3=-lo_huge
    ! weight the fit by size of forces
    logical :: weighted=.false.

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
    integer :: errctr
    logical :: dumlog


    call cli%init(progname    = 'MEGAFIT',&
                  authors     = lo_author,&
                  version     = lo_version,&
                  license     = lo_licence,&
                  help        = 'Usage: ',&
                  description = 'Fit everything to everything all at once.',&
                  examples    = ["MEGAFIT -rc2 5.1           ",&
                                 "MEGAFIT -rc2 4.5 -rc3 3.21 "],&
                  epilog      = new_line('a')//"...")

    ! Real options
    call cli%add(switch='--secondorder_cutoff',switch_ab='-rc2',&
            help='Cutoff for the second order force constants',&
            required=.false.,act='store',def='5.0',error=lo_status)
            if ( lo_status .ne. 0 ) stop
    call cli%add(switch='--thirdorder_cutoff',switch_ab='-rc3',&
            help='Cutoff for the third order force constants',&
            required=.false.,act='store',def='-1',error=lo_status)
            if ( lo_status .ne. 0 ) stop
    call cli%add(switch='--fourthorder_cutoff',switch_ab='-rc4',&
            help='Cutoff for the fourth order force constants',&
            required=.false.,act='store',def='-1',error=lo_status)
            if ( lo_status .ne. 0 ) stop

!    call cli%add(switch='--firstorder',&
!            help='Include the first order force constants. These can be used to find the finite temperature equilibrium structure.',&
!            required=.false.,act='store_true',def='.false.',error=lo_status)
!            if ( lo_status .ne. 0 ) stop
    call cli%add(switch='--magnetic_pair_cutoff',switch_ab='-mc2',&
            help='Cutoff for the pair magnetic interactions',&
            required=.false.,act='store',def='-1.0',error=lo_status)
            if ( lo_status .ne. 0 ) stop

    ! hidden useless options, to print the man page and so on.
    call cli%add(switch='--manpage',hidden=.true.,help='',&
        required=.false.,act='store_true',def='.false.',error=lo_status)
        if ( lo_status .ne. 0 ) stop
    call cli%add(switch='--verbose',hidden=.true.,help='',&
        required=.false.,act='store_true',def='.false.',error=lo_status)
        if ( lo_status .ne. 0 ) stop

    ! hidden and TOP SECRET! (i.e. debugging)
    call cli%add(switch='--wigner_seitz_cutoff',switch_ab='-wz',hidden=.true.,&
        help='Use the Voronoi cell of the supercell, a super Wigner-Seitz cell, for the second order cutoff. The three numbers are the dimensions of the supercell.',&
        nargs='3',required=.false.,act='store',def='-1 -1 -1',exclude='--secondorder_njump',error=lo_status)
        if ( lo_status .ne. 0 ) stop
    call cli%add(switch='--notranspose',hidden=.true.,help='',&
        required=.false.,act='store_true',def='.false.',error=lo_status)
        if ( lo_status .ne. 0 ) stop
    call cli%add(switch='--nospacegroup',hidden=.true.,help='',&
        required=.false.,act='store_true',def='.false.',error=lo_status)
        if ( lo_status .ne. 0 ) stop
    call cli%add(switch='--norotational',hidden=.true.,help='',&
        required=.false.,act='store_true',def='.false.',error=lo_status)
        if ( lo_status .ne. 0 ) stop
    call cli%add(switch='--nohuang',hidden=.true.,help='',&
        required=.false.,act='store_true',def='.false.',error=lo_status)
        if ( lo_status .ne. 0 ) stop
    call cli%add(switch='--nohermitian',hidden=.true.,help='',&
        required=.false.,act='store_true',def='.false.',error=lo_status)
        if ( lo_status .ne. 0 ) stop
    call cli%add(switch='--secondorder_njump',switch_ab='-nj2',&
        help='Second order neighbour jumps.',&
        required=.false.,act='store',def='-1',exclude='-wzdim',error=lo_status,hidden=.true.)
        if ( lo_status .ne. 0 ) stop
    call cli%add(switch='--thirdorder_njump',switch_ab='-nj3',&
        help='Third order neighbour jumps',&
        required=.false.,act='store',def='-1',exclude='--thirdorder_cutoff',error=lo_status,hidden=.true.)
        if ( lo_status .ne. 0 ) stop
    call cli%add(switch='--fourthorder_njump',switch_ab='-nj4',&
        help='Fourth order neighbour jumps',&
        required=.false.,act='store',def='-1',exclude='--fourthorder_cutoff',error=lo_status,hidden=.true.)
        if ( lo_status .ne. 0 ) stop
    call cli%add(switch='--polar',&
        help='Add dipole-dipole corrections for polar materials.',&
        required=.false.,act='store_true',def='.false.',error=lo_status)
        if ( lo_status .ne. 0 ) stop
    call cli%add(switch='--polarcorrectiontype',switch_ab='-pc',&
        help='What kind of polar correction to use.',&
        required=.false.,act='store',def='3',choices='1,2,3',error=lo_status)
        if ( lo_status .ne. 0 ) stop

    ! Specific for MEGAFIT
    call cli%add(switch='--order',switch_ab='-o',&
        help='Order of the polynomials for the grid fitting procedure.',&
        required=.false.,act='store',def='2',error=lo_status)
        if ( lo_status .ne. 0 ) stop
    call cli%add(switch='--evalenergy',&
        help="Evaluate the free energy at the points specified in `infile.evalpoints`",&
        required=.false.,act='store_true',def='.false.',error=lo_status)
        if ( lo_status .ne. 0 ) stop
    call cli%add(switch='--harmonic_qpoint_grid',switch_ab='-qgh',&
        help='Density of q-mesh for harmonic free energy.',nargs='3',required=.false.,&
        act='store',def='26 26 26',error=lo_status);
        if ( lo_status .ne. 0 ) stop
    call cli%add(switch='--anharmonic_qpoint_grid',switch_ab='-qga',&
        help='Density of q-mesh for anharmonic free energy.',&
        nargs='3',required=.false.,act='store',def='10 10 10',error=lo_status)
        if ( lo_status .ne. 0 ) stop
    call cli%add(switch='--quasiharmonic',&
        help='In addition to the full anharmonic free energy, evaluate the quasiharmonic free energy for reference.',&
        required=.false.,act='store_true',def='.false.',error=lo_status)
        if ( lo_status .ne. 0 ) stop
    call cli%add(switch='--temperaturescale',switch_ab='-ts',&
        help='Scale the temperature.',&
        required=.false.,act='store',def='-1',error=lo_status)
        if ( lo_status .ne. 0 ) stop
    call cli%add(switch='--distancescale',switch_ab='-ds',&
        help='Scale distances.',&
        required=.false.,act='store',def='0.1',error=lo_status)
        if ( lo_status .ne. 0 ) stop
    call cli%add(switch='--pairfittype',switch_ab='-pf',&
        help='What method to use then fitting the second order. 1 is a global polynomial, 2 locally adaptive polynomials.',&
        required=.false.,act='store',def='1',choices='1,2',error=lo_status)
        if ( lo_status .ne. 0 ) stop
    call cli%add(switch='--dumpgrid',hidden=.true.,help='',&
        required=.false.,act='store_true',def='.false.',error=lo_status)
        if ( lo_status .ne. 0 ) stop
    call cli%add(switch='--dumpinputgrid',hidden=.true.,help='',&
        required=.false.,act='store_true',def='.false.',error=lo_status)
        if ( lo_status .ne. 0 ) stop
    call cli%add(switch='--diagnostics',hidden=.true.,help='',&
        required=.false.,act='store_true',def='.false.',error=lo_status)
        if ( lo_status .ne. 0 ) stop
    call cli%add(switch='--dielectric_pair_cutoff',switch_ab='-dc2',hidden=.true.,&
        help='Cutoff for the pair dielectric interactions',&
        required=.false.,act='store',def='-1.0',error=lo_status)
        if ( lo_status .ne. 0 ) stop
    call cli%add(switch='--dielectric_triplet_cutoff',switch_ab='-dc3',hidden=.true.,&
        help='Cutoff for the triplet dielectric interactions',&
        required=.false.,act='store',def='-1.0',error=lo_status)
        if ( lo_status .ne. 0 ) stop
    call cli%add(switch='--weighted',hidden=.true.,help='',&
        required=.false.,act='store_true',def='.false.',error=lo_status)
        if ( lo_status .ne. 0 ) stop

    ! actually parse it
    call cli%parse(error=lo_status)
    if ( lo_status .ne. 0 ) stop

    ! generate manpage?
    call cli%get(switch='--manpage',val=dumlog)
    if ( dumlog ) then
        call cli%save_man_page(trim(cli%progname)//'.1')
        call cli%save_usage_to_markdown(trim(cli%progname)//'.md')
        write(*,*) 'Wrote manpage for "'//trim(cli%progname)// '"'
        stop
    endif

    ! verbose output?
    opts%verbosity=0
    call cli%get(switch='--verbose',val=dumlog)
    if ( dumlog ) opts%verbosity=1

    ! Parse the rest of stuff
!    call cli%get(switch='--firstorder',val=opts%firstorder)
    errctr=0
    call cli%get(switch='-rc2',val=opts%cutoff2      ,error=lo_status); errctr=errctr+lo_status
    call cli%get(switch='-rc3',val=opts%cutoff3      ,error=lo_status); errctr=errctr+lo_status
    call cli%get(switch='-rc4',val=opts%cutoff4      ,error=lo_status); errctr=errctr+lo_status
    call cli%get(switch='-wz',val=opts%wzdimensions  ,error=lo_status); errctr=errctr+lo_status
    call cli%get(switch='--norotational',val=dumlog  ,error=lo_status); errctr=errctr+lo_status
    opts%rotationalconstraints=.not.dumlog
    call cli%get(switch='--nohuang',val=dumlog       ,error=lo_status); errctr=errctr+lo_status
    opts%huanginvariance=.not.dumlog
    call cli%get(switch='--notranspose',val=dumlog   ,error=lo_status); errctr=errctr+lo_status
    opts%transposition=.not.dumlog
    call cli%get(switch='--nospacegroup',val=dumlog  ,error=lo_status); errctr=errctr+lo_status
    opts%spacegroup=.not.dumlog
    call cli%get(switch='--nohermitian',val=dumlog   ,error=lo_status); errctr=errctr+lo_status
    opts%hermitian=.not.dumlog
    call cli%get(switch='-nj2',val=opts%njump2       ,error=lo_status); errctr=errctr+lo_status
    call cli%get(switch='-nj3',val=opts%njump3       ,error=lo_status); errctr=errctr+lo_status
    call cli%get(switch='-nj4',val=opts%njump4       ,error=lo_status); errctr=errctr+lo_status
    call cli%get(switch='-mc2',val=opts%magcutoff2   ,error=lo_status); errctr=errctr+lo_status
    call cli%get(switch='--order',val=opts%order     ,error=lo_status); errctr=errctr+lo_status

    call cli%get(switch='--evalenergy',val=opts%evalenergy,                  error=lo_status); errctr=errctr+lo_status
    call cli%get(switch='--temperaturescale',val=opts%temperature_scale     ,error=lo_status); errctr=errctr+lo_status
    call cli%get(switch='--distancescale',val=opts%distance_scale           ,error=lo_status); errctr=errctr+lo_status
    call cli%get(switch='--harmonic_qpoint_grid' ,val=opts%qgrid_harm       ,error=lo_status); errctr=errctr+lo_status
    call cli%get(switch='--anharmonic_qpoint_grid',val=opts%qgrid_anharm    ,error=lo_status); errctr=errctr+lo_status
    call cli%get(switch='--quasiharmonic',val=opts%quasiharmonic            ,error=lo_status); errctr=errctr+lo_status
    call cli%get(switch='--polar',val=opts%polar                            ,error=lo_status); errctr=errctr+lo_status
    call cli%get(switch='--polarcorrectiontype',val=opts%polarcorrectiontype,error=lo_status); errctr=errctr+lo_status
    call cli%get(switch='--pairfittype',val=opts%pairfittype,                error=lo_status); errctr=errctr+lo_status
    call cli%get(switch='--dumpgrid',val=opts%dumpfullgrid,                  error=lo_status); errctr=errctr+lo_status
    call cli%get(switch='--dumpinputgrid',val=opts%dumpinputgrid,            error=lo_status); errctr=errctr+lo_status
    call cli%get(switch='--diagnostics',val=opts%diagnostics,                error=lo_status); errctr=errctr+lo_status
    call cli%get(switch='--weighted',val=opts%weighted,                      error=lo_status); errctr=errctr+lo_status
    call cli%get(switch='-dc2',val=opts%dielcutoff2)
    call cli%get(switch='-dc3',val=opts%dielcutoff3)
    if ( errctr .gt. 0 ) stop

    ! Convert to atomic units right away
    opts%cutoff2=opts%cutoff2*lo_A_to_bohr
    opts%cutoff3=opts%cutoff3*lo_A_to_bohr
    opts%cutoff4=opts%cutoff4*lo_A_to_bohr
    opts%magcutoff2=opts%magcutoff2*lo_A_to_bohr
    opts%dielcutoff2=opts%dielcutoff2*lo_A_to_Bohr
    opts%dielcutoff3=opts%dielcutoff3*lo_A_to_Bohr

end subroutine

end module
