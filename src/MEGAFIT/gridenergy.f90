#include "precompilerdefinitions"
module gridenergy
!! Evaluate forceconstants and associated quantities at arbitrary points
use konstanter, only: flyt,r8,i8,lo_huge,lo_hugeint,lo_pi,lo_twopi,lo_imag,lo_status,lo_exitcode_param,lo_exitcode_symmetry,&
                      lo_tol,lo_sqtol,lo_pressure_HartreeBohr_to_GPa,lo_pressure_GPa_to_HartreeBohr,lo_Hartree_to_eV,&
                      lo_volume_bohr_to_A,lo_volume_A_to_bohr,lo_freqtol
use gottochblandat, only: walltime,tochar,lo_progressbar_init,lo_progressbar,lo_looptimer,lo_sqnorm,lo_mean,&
                          lo_planck,open_file,lo_does_file_exist,lo_flattentensor,lo_linspace,qsort,lo_return_unique
use mpi_wrappers, only: lo_mpi_helper,lo_stop_gracefully,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_IN_PLACE
use lo_memtracker, only: lo_mem_helper
use geometryfunctions, only: lo_linesegment,lo_convex_hull_2d
use type_crystalstructure, only: lo_crystalstructure
use type_forceconstant_firstorder,  only: lo_forceconstant_firstorder
use type_forceconstant_secondorder, only: lo_forceconstant_secondorder
use type_forceconstant_thirdorder,  only: lo_forceconstant_thirdorder
use type_forceconstant_fourthorder, only: lo_forceconstant_fourthorder
use type_jij_secondorder, only: lo_jij_secondorder
use type_forcemap, only: lo_forcemap,lo_secondorder_rot_herm_huang
use type_qpointmesh, only: lo_qpoint_mesh,lo_fft_mesh,lo_generate_qmesh
use type_phonon_dispersions, only: lo_phonon_dispersions
use lo_phonon_bandstructure_on_path, only: lo_phonon_bandstructure
use type_blas_lapack_wrappers, only: lo_dgels
use hdf5_wrappers, only: lo_hdf5_helper,lo_h5_store_data,lo_h5_store_attribute,HID_T,H5F_ACC_TRUNC_F,&
                         h5close_f,h5open_f,h5fclose_f,h5fopen_f,h5fcreate_f,h5gclose_f,h5gopen_f,h5gcreate_f

use type_gridsim, only: lo_gridsim
use type_equation_of_state, only: lo_eos,lo_eos_1d,lo_eos_2d,lo_eos_birch_murnaghan,lo_eos_vinet,lo_eos_2d_birch_murnaghan
use type_polynomial_interpolation, only: lo_grid_interpolation

implicit none

private
public :: lo_gridenergy

!> quasiharmonic volume grid
type lo_gridenergy_vol
    !> how many volumes
    integer :: nv=-lo_hugeint
    integer :: nt=-lo_hugeint
    !> volume axis
    real(flyt), dimension(:), allocatable :: volume
    !> temperature axis
    real(flyt), dimension(:), allocatable :: temperature
    !> energies
    real(flyt), dimension(:,:), allocatable :: U,U0,fph,ah3,ah4
end type

!> volume-temperature grid
type lo_gridenergy_vol_temp
    !> how many volumes
    integer :: nv=-lo_hugeint
    !> how many temperatures
    integer :: nt=-lo_hugeint
    !> volume axis
    real(flyt), dimension(:,:), allocatable :: volume
    !> temperature axis
    real(flyt), dimension(:), allocatable :: temperature
    !> energies
    real(flyt), dimension(:,:), allocatable :: U,U0,fph,ah3,ah4
    real(flyt), dimension(:,:), allocatable :: qh_fph,qh_ah3,qh_ah4
end type

!> volume-temperature-eta grid
type lo_gridenergy_vol_temp_eta
    !> how many volumes
    integer :: nv=-lo_hugeint
    !> how many temperatures
    integer :: nt=-lo_hugeint
    !> how many eta
    integer :: neta=-lo_hugeint
    !> volume axis
    real(flyt), dimension(:), allocatable :: volume
    !> temperature axis
    real(flyt), dimension(:), allocatable :: temperature
    !> eta axis axis
    real(flyt), dimension(:), allocatable :: eta
    !> energies
    real(flyt), dimension(:,:,:), allocatable :: U,U0,fph,ah3,ah4
    real(flyt), dimension(:,:,:), allocatable :: qh_fph,qh_ah3,qh_ah4
end type


!> evaluate the free energy across the mesh
type lo_gridenergy
    !> how many dimensions
    integer :: ndim
    !> gridpoints per dimension
    integer, dimension(:), allocatable :: pts_per_dim
    !> min per dimension
    real(flyt), dimension(:), allocatable :: min_per_dim
    !> min max per dimension
    real(flyt), dimension(:), allocatable :: max_per_dim
    !> evaluate quasiharmonic energyies
    logical :: quasiharmonic=.false.
    !> what kind of grid do we have
    integer :: gridtype
    !> Do it slightly differently depending on the number of dimenisons
    !type(lo_gridenergy_2d) :: grid_2d
    type(lo_gridenergy_vol) :: grid_V
    type(lo_gridenergy_vol_temp) :: grid_VT
    type(lo_gridenergy_vol_temp_eta) :: grid_VTeta
    contains
        !> create everything
        procedure :: generate
end type

! how oftern to report the slow loops
real(flyt), parameter :: timereport=30.0_flyt
integer, parameter :: pm_vgrid=1
integer, parameter :: pm_vtgrid=2
integer, parameter :: pm_vtetagrid=3

contains

#include "gridenergy_aux.f90"
#include "gridenergy_vt_grid.f90"
#include "gridenergy_vteta_grid.f90"
#include "gridenergy_impossible.f90"

!> evaluate the irreducible representation at a certain point
subroutine generate(ge,gs,map,qgrid_harm,qgrid_anharm,quasiharmonic,dumpgrid,mw,mem)
    !> interpolation grid
    class(lo_gridenergy), intent(out) :: ge
    !> simulation grid
    type(lo_gridsim), intent(inout) :: gs
    !> forcemap
    type(lo_forcemap), intent(inout) :: map
    !> q-point mesh dimensions
    integer, dimension(3) :: qgrid_harm
    !> support q-grid mesh dimensions
    integer, dimension(3) :: qgrid_anharm
    !> should I calculate the quasiharmonic as well
    logical, intent(in) :: quasiharmonic
    !> should I dump the input files for the entire grid
    logical, intent(in) :: dumpgrid
    !> MPI helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem

    character(len=3), dimension(gs%ndim) :: pointspacing
    real(flyt), dimension(:,:), allocatable :: gridcoord,gridenergy
    real(flyt), dimension(:), allocatable :: adapt_minvol,adapt_maxvol
    real(flyt) :: pressurestep,timer
    integer, dimension(:,:), allocatable :: gridind
    integer :: verbosity,evalmode

    ! read the input file and figure out some stuff
    init: block
        integer :: i,u
        if ( mw%talk ) then
            write(*,*) ''
            write(*,*) 'INTERPOLATING FREE ENERGY'
            timer=walltime()
            verbosity=1
        else
            timer=walltime()
            verbosity=-1
        endif
        ! grab dimensions and grid stuff from file
        ge%ndim=gs%ndim
        u=open_file('in','infile.evalpoints')
            read(u,*) evalmode
            select case(evalmode)
            case(1)
                ! generic way of definig the mesh
                lo_allocate(ge%min_per_dim(ge%ndim))
                lo_allocate(ge%max_per_dim(ge%ndim))
                lo_allocate(ge%pts_per_dim(ge%ndim))
                ge%min_per_dim=0.0_flyt
                ge%max_per_dim=0.0_flyt
                ge%pts_per_dim=0
                do i=1,ge%ndim
                    read(u,*) ge%pts_per_dim(i),pointspacing(i)
                    read(u,*) ge%min_per_dim(i),ge%max_per_dim(i)
                enddo
                read(u,*) pressurestep
                if ( ge%ndim .eq. 1 .and. gs%info%dim_volume .gt. 0 ) then
                    write(*,*) 'onlyqh'
                endif
                ! Convert to atomic units
                i=gs%info%dim_volume
                if ( i .gt. 0 ) then
                    ge%min_per_dim(i)=ge%min_per_dim(i)*lo_volume_A_to_bohr
                    ge%max_per_dim(i)=ge%max_per_dim(i)*lo_volume_A_to_bohr
                endif
                pressurestep=pressurestep*lo_pressure_GPa_to_HartreeBohr
            case(3)
                ! specialized that only works in 2D now, and with V-T as the grid
                lo_allocate(ge%min_per_dim(ge%ndim))
                lo_allocate(ge%max_per_dim(ge%ndim))
                lo_allocate(ge%pts_per_dim(ge%ndim))
                ge%min_per_dim=0.0_flyt
                ge%max_per_dim=0.0_flyt
                ge%pts_per_dim=0

                i=gs%info%dim_temperature
                read(u,*) ge%pts_per_dim(i),pointspacing(i)
                read(u,*) ge%min_per_dim(i),ge%max_per_dim(i)
                i=gs%info%dim_volume
                read(u,*) ge%pts_per_dim(i),pointspacing(i)
                read(u,*) pressurestep
                ! Instead of reading the volume limits from file, get them per temperature
                ! from the convex hull.
                allocate(adapt_minvol( ge%pts_per_dim(gs%info%dim_temperature) ))
                allocate(adapt_maxvol( ge%pts_per_dim(gs%info%dim_temperature) ))
                adapt_minvol=0.0_flyt
                adapt_minvol=0.0_flyt
                i=gs%info%dim_temperature
                call volume_limits_from_convex_hull(gs%grid_coordinates,gs%info%dim_volume,gs%info%dim_temperature,&
                    ge%min_per_dim(i),ge%max_per_dim(i),ge%pts_per_dim(i),pointspacing(i),adapt_minvol,adapt_maxvol,verbosity)
                ! Convert to atomic units
                pressurestep=pressurestep*lo_pressure_GPa_to_HartreeBohr
            case default
                call lo_stop_gracefully(['NOT DONE'],lo_exitcode_param,__FILE__,__LINE__,mw%comm)
            end select
        close(u)
        if ( mw%talk ) then
            write(*,*) '              ndim:',ge%ndim
        endif
        ! Figure out what kind of grid we have.
        select case(gs%ndim)
        case(pm_vgrid)
            if ( gs%info%dim_volume .gt. 0 ) then
                ! V grid
                ge%gridtype=pm_vgrid !  1
            else
                call lo_stop_gracefully(['1-D NOT DONE'],lo_exitcode_param,__FILE__,__LINE__,mw%comm)
            endif
        case(pm_vtgrid)
            if ( gs%info%dim_temperature .gt. 0 .and. gs%info%dim_volume .gt. 0 ) then
                ! V-T grid
                ge%gridtype=pm_vtgrid !  1
                if ( mw%talk ) then
                    write(*,*) '          gridtype:',ge%gridtype,'(V-T grid)'
                    write(*,*) '   temperature-dim:',gs%info%dim_temperature
                    write(*,*) '        volume-dim:',gs%info%dim_volume
                endif
            else
                call lo_stop_gracefully(['NOT DONE'],lo_exitcode_param,__FILE__,__LINE__,mw%comm)
                stop
            endif
        case(pm_vtetagrid)
            if ( gs%info%dim_temperature .gt. 0 .and. gs%info%dim_volume .gt. 0 .and. gs%info%dim_eta .gt. 0 ) then
                ! V-T-eta grid
                ge%gridtype=pm_vtetagrid !2
                if ( mw%talk ) then
                    write(*,*) '          gridtype:',ge%gridtype,'(V-T-eta grid)'
                    write(*,*) '   temperature-dim:',gs%info%dim_temperature
                    write(*,*) '        volume-dim:',gs%info%dim_volume
                    write(*,*) '           eta-dim:',gs%info%dim_eta
                endif
            else
                call lo_stop_gracefully(['NOT DONE'],lo_exitcode_param,__FILE__,__LINE__,mw%comm)
            endif
        end select
    end block init

    ! Set the flat-ish grid
    setgrid: block
        real(flyt), dimension(:), allocatable :: dumvol
        real(flyt) :: f0,f1
        integer :: npts,i,j,k,l,ii,jj,kk

        npts=product(ge%pts_per_dim)
        lo_allocate(gridcoord(ge%ndim,npts))
        lo_allocate(gridenergy(8,npts))
        lo_allocate(gridind(ge%ndim,npts))
        gridcoord=0.0_flyt
        gridenergy=0.0_flyt
        gridind=0

        ! build the grid
        select case(ge%gridtype)
        case(pm_vgrid) ! V-grid
            ge%grid_V%nv=ge%pts_per_dim( gs%info%dim_volume )
            write(*,*) 'FIXME VOLUME'
            stop
        case(pm_vtgrid) ! V-T grid
            ge%grid_VT%nv=ge%pts_per_dim( gs%info%dim_volume )
            ge%grid_VT%nt=ge%pts_per_dim( gs%info%dim_temperature )
            allocate( ge%grid_VT%U     ( ge%grid_VT%nv, ge%grid_VT%nt ) )
            allocate( ge%grid_VT%U0    ( ge%grid_VT%nv, ge%grid_VT%nt ) )
            allocate( ge%grid_VT%fph   ( ge%grid_VT%nv, ge%grid_VT%nt ) )
            allocate( ge%grid_VT%ah3   ( ge%grid_VT%nv, ge%grid_VT%nt ) )
            allocate( ge%grid_VT%ah4   ( ge%grid_VT%nv, ge%grid_VT%nt ) )
            allocate( ge%grid_VT%qh_fph( ge%grid_VT%nv, ge%grid_VT%nt ) )
            allocate( ge%grid_VT%qh_ah3( ge%grid_VT%nv, ge%grid_VT%nt ) )
            allocate( ge%grid_VT%qh_ah4( ge%grid_VT%nv, ge%grid_VT%nt ) )
            ge%grid_VT%U     =0.0_flyt
            ge%grid_VT%U0    =0.0_flyt
            ge%grid_VT%fph   =0.0_flyt
            ge%grid_VT%ah3   =0.0_flyt
            ge%grid_VT%ah4   =0.0_flyt
            ge%grid_VT%qh_fph=0.0_flyt
            ge%grid_VT%qh_ah3=0.0_flyt
            ge%grid_VT%qh_ah4=0.0_flyt

            lo_allocate(ge%grid_VT%volume( ge%grid_VT%nv,ge%grid_VT%nt ))
            lo_allocate(ge%grid_VT%temperature( ge%grid_VT%nt ))
            ge%grid_VT%volume=0.0_flyt
            ge%grid_VT%temperature=0.0_flyt
            lo_allocate(dumvol( ge%grid_VT%nv ))

            ! fix volume axis
            select case(evalmode)
            case(1)
                ! same volume for all temperatures
                f0=ge%min_per_dim( gs%info%dim_volume )
                f1=ge%max_per_dim( gs%info%dim_volume )
                select case( pointspacing(gs%info%dim_volume) )
                case('lin') ! linearly spaced in volume
                    call lo_linspace(f0,f1,dumvol)
                case('den') ! linearly spaced in density
                    call lo_linspace(1.0_flyt/f1,1.0_flyt/f0,dumvol)
                    dumvol=1.0_flyt/dumvol
                end select
                do i=1,ge%grid_VT%nt
                    ge%grid_VT%volume(:,i)=dumvol
                enddo
            case(3)
                ! different volumes for different temperatures
                do i=1,ge%grid_VT%nt
                    f0=adapt_minvol(i)
                    f1=adapt_maxvol(i)
                    select case( pointspacing(gs%info%dim_volume) )
                    case('lin') ! linearly spaced in volume
                        call lo_linspace(f0,f1,dumvol)
                    case('den') ! linearly spaced in density
                        call lo_linspace(1.0_flyt/f1,1.0_flyt/f0,dumvol)
                        dumvol=1.0_flyt/dumvol
                    end select
                    ge%grid_VT%volume(:,i)=dumvol
                enddo
            case default
            end select

            ! fix temperature axis
            f0=ge%min_per_dim( gs%info%dim_temperature )
            f1=ge%max_per_dim( gs%info%dim_temperature )
            select case( pointspacing(gs%info%dim_temperature) )
            case('lin') ! linearly spaced in temperature
                call lo_linspace(f0,f1,ge%grid_VT%temperature)
            case('log') ! log=spaced
                call logspace(f0,f1,ge%grid_VT%temperature)
            end select
            if ( mw%talk ) then
                write(*,*) '      temperatures: ',tochar(minval(ge%grid_VT%temperature)),' -> ',tochar(maxval(ge%grid_VT%temperature)),' with ',tochar(ge%grid_VT%nt),' points'
                write(*,*) '           volumes: ',tochar(minval(ge%grid_VT%volume*lo_volume_bohr_to_A)),' -> ',tochar(maxval(ge%grid_VT%volume*lo_volume_bohr_to_A)),' with ',tochar(ge%grid_VT%nv),' points'
            endif
            l=0
            ii=gs%info%dim_volume
            jj=gs%info%dim_temperature
            do i=1,ge%grid_VT%nv
            do j=1,ge%grid_VT%nt
                l=l+1
                gridind(:,l)=[i,j]
                gridcoord(ii,l)=ge%grid_VT%volume( i,j )
                gridcoord(jj,l)=ge%grid_VT%temperature( j )
            enddo
            enddo
        case(pm_vtetagrid) ! V-T-eta grid
            ge%grid_VTeta%nv=ge%pts_per_dim( gs%info%dim_volume )
            ge%grid_VTeta%nt=ge%pts_per_dim( gs%info%dim_temperature )
            ge%grid_VTeta%neta=ge%pts_per_dim( gs%info%dim_eta )
            allocate( ge%grid_VTeta%U     ( ge%grid_VTeta%nv, ge%grid_VTeta%nt, ge%grid_VTeta%neta ) )
            allocate( ge%grid_VTeta%U0    ( ge%grid_VTeta%nv, ge%grid_VTeta%nt, ge%grid_VTeta%neta ) )
            allocate( ge%grid_VTeta%fph   ( ge%grid_VTeta%nv, ge%grid_VTeta%nt, ge%grid_VTeta%neta ) )
            allocate( ge%grid_VTeta%ah3   ( ge%grid_VTeta%nv, ge%grid_VTeta%nt, ge%grid_VTeta%neta ) )
            allocate( ge%grid_VTeta%ah4   ( ge%grid_VTeta%nv, ge%grid_VTeta%nt, ge%grid_VTeta%neta ) )
            allocate( ge%grid_VTeta%qh_fph( ge%grid_VTeta%nv, ge%grid_VTeta%nt, ge%grid_VTeta%neta ) )
            allocate( ge%grid_VTeta%qh_ah3( ge%grid_VTeta%nv, ge%grid_VTeta%nt, ge%grid_VTeta%neta ) )
            allocate( ge%grid_VTeta%qh_ah4( ge%grid_VTeta%nv, ge%grid_VTeta%nt, ge%grid_VTeta%neta ) )
            ge%grid_VTeta%U     =0.0_flyt
            ge%grid_VTeta%U0    =0.0_flyt
            ge%grid_VTeta%fph   =0.0_flyt
            ge%grid_VTeta%ah3   =0.0_flyt
            ge%grid_VTeta%ah4   =0.0_flyt
            ge%grid_VTeta%qh_fph=0.0_flyt
            ge%grid_VTeta%qh_ah3=0.0_flyt
            ge%grid_VTeta%qh_ah4=0.0_flyt

            allocate(ge%grid_VTeta%volume( ge%grid_VTeta%nv ))
            allocate(ge%grid_VTeta%temperature( ge%grid_VTeta%nt ))
            allocate(ge%grid_VTeta%eta( ge%grid_VTeta%neta ))
            ! fix volume axis
            f0=ge%min_per_dim( gs%info%dim_volume )
            f1=ge%max_per_dim( gs%info%dim_volume )
            select case( pointspacing(gs%info%dim_volume) )
            case('lin') ! linearly spaced in volume
                call lo_linspace(f0,f1,ge%grid_VTeta%volume)
            case('den') ! linearly spaced in density
                call lo_linspace(1.0_flyt/f1,1.0_flyt/f0,ge%grid_VTeta%volume)
                ge%grid_VTeta%volume=1.0_flyt/ge%grid_VTeta%volume
            end select
            ! fix temperature axis
            f0=ge%min_per_dim( gs%info%dim_temperature )
            f1=ge%max_per_dim( gs%info%dim_temperature )
            select case( pointspacing(gs%info%dim_temperature) )
            case('lin') ! linearly spaced in temperature
                call lo_linspace(f0,f1,ge%grid_VTeta%temperature)
            case('log') ! log=spaced
                call logspace(f0,f1,ge%grid_VTeta%temperature)
            end select
            ! fix eta axis
            f0=ge%min_per_dim( gs%info%dim_eta )
            f1=ge%max_per_dim( gs%info%dim_eta )
            select case( pointspacing(gs%info%dim_eta) )
            case('lin') ! linearly spaced in eta
                call lo_linspace(f0,f1,ge%grid_VTeta%eta)
            case('log') ! log=spaced
                call logspace(f0,f1,ge%grid_VTeta%eta)
            end select
            ! and get the grid coordinates
            l=0
            ii=gs%info%dim_volume
            jj=gs%info%dim_temperature
            kk=gs%info%dim_eta
            do i=1,ge%grid_VTeta%nv
            do j=1,ge%grid_VTeta%nt
            do k=1,ge%grid_VTeta%neta
                l=l+1
                gridind(:,l)=[i,j,k]
                gridcoord(ii,l)=ge%grid_VTeta%volume( i )
                gridcoord(jj,l)=ge%grid_VTeta%temperature( j )
                gridcoord(kk,l)=ge%grid_VTeta%eta( k )
            enddo
            enddo
            enddo

            if ( mw%talk ) then
                write(*,*) '      temperatures: ',tochar(minval(ge%grid_VTeta%temperature)),' -> ',tochar(maxval(ge%grid_VTeta%temperature)),' with ',tochar(ge%grid_VTeta%nt),' points'
                write(*,*) '           volumes: ',tochar(minval(ge%grid_VTeta%volume)),' -> ',tochar(maxval(ge%grid_VTeta%volume)),' with ',tochar(ge%grid_VTeta%nv),' points'
                write(*,*) '               eta: ',tochar(minval(ge%grid_VTeta%eta)),' -> ',tochar(maxval(ge%grid_VTeta%eta)),' with ',tochar(ge%grid_VTeta%neta),' points'
            endif
        case default
            call lo_stop_gracefully(['NOT DONE'],lo_exitcode_param,__FILE__,__LINE__,mw%comm)
        end select
    end block setgrid

    ! if ( dumpgrid .and. mw%talk ) then
    ! dgrid: block
    !     type(lo_crystalstructure) :: p
    !     type(lo_forceconstant_secondorder) :: fc
    !     type(lo_forceconstant_thirdorder) :: fct
    !     type(lo_forceconstant_fourthorder) :: fcf
    !     type(lo_jij_secondorder) :: jij
    !     real(flyt), dimension(:,:), allocatable :: pairconstraints
    !     integer :: i,j,l,nconstr
    !
    !     write(*,*) '... dumping input files across the entire grid'
    !
    !     select case(ge%gridtype)
    !     case(pm_vtgrid)
    !         l=0
    !         do i=1,ge%grid_VT%nv
    !         do j=1,ge%grid_VT%nt
    !             l=l+1
    !             ! Get the structure
    !             call gs%structure%interpolate(gridcoord(:,l),p,gs%info%dim_volume)
    !             call lo_secondorder_rot_herm_huang( map,p,pairconstraints,nconstr,.true.,.true.,.true. )
    !             if ( nconstr .gt. 0 ) then
    !                 call gs%eval(map,gridcoord(:,l),pairconstraints)
    !             else
    !                 call gs%eval(map,gridcoord(:,l))
    !             endif
    !             p%info%title='interp struct V= '//tochar(ge%grid_VT%volume(i,j)*lo_volume_bohr_to_A,ndecimals=11)//' T= '//tochar(ge%grid_VT%temperature(j),ndecimals=10)
    !             call p%writetofile('uc_'//tochar(i)//'_'//tochar(j),1)
    !
    !             ! Now dump all the component guys
    !             if ( gs%info%secondorder ) then
    !                 call map%get_secondorder_forceconstant(p,fc,mem,-1)
    !                 call fc%writetofile(p,'fc2_'//tochar(i)//'_'//tochar(j))
    !             endif
    !             if ( gs%info%thirdorder ) then
    !                 call map%get_thirdorder_forceconstant(p,fct)
    !                 call fct%writetofile(p,'fc3_'//tochar(i)//'_'//tochar(j))
    !             endif
    !             if ( gs%info%fourthorder ) then
    !                 call map%get_fourthorder_forceconstant(p,fcf)
    !                 call fcf%writetofile(p,'fc4_'//tochar(i)//'_'//tochar(j))
    !             endif
    !             ! if ( gs%info%magnetic_pair_interactions ) then
    !             !     call map%get_secondorder_jij(p,jij)
    !             !     call jij%writetofile(p,'jij_'//tochar(i)//'_'//tochar(j))
    !             ! endif
    !         enddo
    !         enddo
    !     case default
    !         write(*,*) 'FIXME DUMPGRID'
    !         stop
    !     end select
    !
    ! end block dgrid
    ! endif

    ! Evaluate stuff
    evalen: block
        type(lo_mpi_helper) :: ml
        real(flyt), dimension(gs%ndim) :: dcrd
        real(flyt) :: t0,timer_anharmonic,temperature,volume,eta
        integer :: i,l,ii,jj,kk,npts

        ! Split communicator to serial things
        call mw%split(ml,mw%r,__FILE__,__LINE__)

        npts=size(gridcoord,2)
        ! First intepolate U0 and static energy
        if ( mw%talk ) then
            t0=walltime()
            call lo_progressbar_init()
        endif
        do i=1,npts
            ! to make it MPI parallel
            if ( mod(i,mw%n) .ne. mw%r ) cycle
            ! interpolate the internal energy
            call gs%energy%interpolate( gridcoord(:,i),gridenergy(2,i) )
            ! add the energy from the provided equation of state
            select type(eos=>gs%eos)
            class is(lo_eos_1d)
                volume=gridcoord(gs%info%dim_volume,i)
                gridenergy(1,i)=eos%energy_from_volume( volume )
            class is(lo_eos_2d)
                volume=gridcoord(gs%info%dim_volume,i)
                eta=gridcoord(gs%info%dim_eta,i)
                gridenergy(1,i)=eos%energy_from_volume_eta( volume,eta )
            class default
                ! with no equation of state, do nothing.
                gridenergy(1,i)=0.0_flyt
            end select
            if ( mw%talk .and. i .lt. npts ) call lo_progressbar(' ... interpolating U0',i,npts,walltime()-t0)
        enddo
        if ( mw%talk ) call lo_progressbar(' ... interpolating U0',npts,npts,walltime()-t0)

        ! Now the phonon free energy
        if ( mw%talk ) then
            t0=walltime()
            call lo_progressbar_init()
        endif
        do i=1,npts
            if ( mod(i,mw%n) .ne. mw%r ) cycle
            temperature=gridcoord( gs%info%dim_temperature,i )
            call phonon_free_energy_for_single_point( gs,map,gridcoord(:,i),qgrid_harm,temperature,ml,mem,gridenergy(3,i) )
            if ( quasiharmonic ) then
                dcrd=gridcoord(:,i)
                dcrd( gs%info%dim_temperature )=0.0_flyt
                call phonon_free_energy_for_single_point( gs,map,dcrd,qgrid_harm,temperature,ml,mem,gridenergy(6,i) )
            endif
            if ( mw%talk .and. i .lt. npts ) call lo_progressbar(' ... phonon free energy',i,npts,walltime()-t0)
        enddo
        if ( mw%talk ) call lo_progressbar(' ... phonon free energy',npts,npts,walltime()-t0)

        ! And destroy the temporary communicators
        call ml%free(__FILE__,__LINE__)

        ! Add it up over ranks, what we have so far
        call mpi_allreduce(MPI_IN_PLACE,gridenergy,8*npts,MPI_DOUBLE_PRECISION,MPI_SUM,mw%comm,mw%error )

        ! Now the anharmonic free energy. This is parallel internally instead.
        if ( map%have_fc_triplet ) then
        if ( map%have_fc_quartet ) then
        if ( qgrid_anharm(1) .gt. 0 ) then
            if ( mw%talk ) then
                timer_anharmonic=walltime()
                t0=walltime()
            endif
            do i=1,npts
                temperature=gridcoord( gs%info%dim_temperature,i )
                call anharmonic_free_energy_for_single_point( gs,map,gridcoord(:,i),qgrid_anharm,temperature,&
                                                              gridenergy(4,i),gridenergy(5,i),mw,mem )
                if ( quasiharmonic ) then
                   dcrd=gridcoord(:,i)
                   dcrd( gs%info%dim_temperature )=0.0_flyt
                   call anharmonic_free_energy_for_single_point( gs,map,dcrd,qgrid_anharm,temperature,&
                                                                 gridenergy(7,i),gridenergy(8,i),mw,mem )
                endif
                if ( mw%talk ) then
                if ( walltime()-t0 .gt. timereport ) then
                    call lo_looptimer('... anharmonic free energy',timer_anharmonic,walltime(),i,npts)
                    t0=walltime()
                endif
                endif
            enddo
        endif
        endif
        endif

        ! Store it organized
        select case(ge%gridtype)
        case(pm_vgrid) ! V-line
            write(*,*) 'FIXME PURE VOLUME GRID'; stop
        case(pm_vtgrid) ! V-T grid
            do l=1,npts
                ii=gridind(1,l)
                jj=gridind(2,l)
                ge%grid_VT%U     ( ii,jj )=gridenergy(1,l)
                ge%grid_VT%U0    ( ii,jj )=gridenergy(2,l)
                ge%grid_VT%fph   ( ii,jj )=gridenergy(3,l)
                ge%grid_VT%ah3   ( ii,jj )=gridenergy(4,l)
                ge%grid_VT%ah4   ( ii,jj )=gridenergy(5,l)
                ge%grid_VT%qh_fph( ii,jj )=gridenergy(6,l)
                ge%grid_VT%qh_ah3( ii,jj )=gridenergy(7,l)
                ge%grid_VT%qh_ah4( ii,jj )=gridenergy(8,l)
            enddo
        case(pm_vtetagrid) ! V-T-eta grid
            do l=1,npts
                ii=gridind(1,l)
                jj=gridind(2,l)
                kk=gridind(3,l)
                ge%grid_VTeta%U     ( ii,jj,kk )=gridenergy(1,l)
                ge%grid_VTeta%U0    ( ii,jj,kk )=gridenergy(2,l)
                ge%grid_VTeta%fph   ( ii,jj,kk )=gridenergy(3,l)
                ge%grid_VTeta%ah3   ( ii,jj,kk )=gridenergy(4,l)
                ge%grid_VTeta%ah4   ( ii,jj,kk )=gridenergy(5,l)
                ge%grid_VTeta%qh_fph( ii,jj,kk )=gridenergy(6,l)
                ge%grid_VTeta%qh_ah3( ii,jj,kk )=gridenergy(7,l)
                ge%grid_VTeta%qh_ah4( ii,jj,kk )=gridenergy(8,l)
            enddo
        end select
    end block evalen

    ! then finalize it
    select case(ge%gridtype)
    case(pm_vgrid)
write(*,*) 'FIXME PURE VOLUME GRID'; stop
    case(pm_vtgrid)
        call VTfinalize( ge%grid_VT,gs,map,'outfile.interpolated_free_energy.hdf5',qgrid_harm,quasiharmonic,pressurestep,dumpgrid,mw,mem )
    case(pm_vtetagrid)
        call VTetafinalize( ge%grid_VTeta,gs,map,'outfile.interpolated_free_energy.hdf5',qgrid_harm,quasiharmonic,pressurestep,mw,mem )
    end select
end subroutine

!> anharmonic free energy at a single point
subroutine anharmonic_free_energy_for_single_point( gs,map,depvar,qgrid,temperature,free_energy_thirdorder,free_energy_fourthorder,mw,mem )
    !> simulation grid
    type(lo_gridsim), intent(inout) :: gs
    !> forcemap
    type(lo_forcemap), intent(inout) :: map
    !> where to evaluate
    real(flyt), dimension(:), intent(in) :: depvar
    !> q-mesh density
    integer, dimension(3), intent(in) :: qgrid
    !> temperature
    real(flyt), intent(in) :: temperature
    !> free energy
    real(flyt), intent(out) :: free_energy_thirdorder,free_energy_fourthorder
    !> MPI helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem

    class(lo_qpoint_mesh), allocatable :: qp
    type(lo_phonon_dispersions) :: dr
    type(lo_crystalstructure) :: p
    type(lo_forceconstant_secondorder) :: fc
    type(lo_forceconstant_thirdorder) :: fct
    type(lo_forceconstant_fourthorder) :: fcf
    real(flyt), dimension(:,:), allocatable :: pairconstraints
    integer :: nconstr

    free_energy_thirdorder  = 0.0_flyt
    free_energy_fourthorder = 0.0_flyt
    ! setup stuffs, first get a structure
    call gs%structure%interpolate(depvar,p,gs%info%dim_volume)
    call p%classify('wedge',timereversal=.true.)
    ! forceconstants, first the constraints
    !call lo_secondorder_rot_herm_huang( map,p,pairconstraints,nconstr,.true.,.true.,.true. )
    nconstr=0
    if ( nconstr .gt. 0 ) then
        call gs%eval(map,depvar,pairconstraints)
    else
        call gs%eval(map,depvar)
    endif
    call map%get_secondorder_forceconstant(p,fc,mem,-1)
    call map%get_thirdorder_forceconstant(p,fct)
    call map%get_fourthorder_forceconstant(p,fcf)
    ! q-points
    call lo_generate_qmesh(qp,p,qgrid,'fft',timereversal=.true.,headrankonly=.false.,mw=mw,mem=mem,verbosity=-1)
    ! dispersions
    call dr%generate(qp,fc,p,mw=mw,mem=mem,verbosity=-1)
    ! check for unstable modes right away
    if ( dr%omega_min .lt. lo_freqtol ) then
        ! kinda dumb, but I have no better idea right now.
        free_energy_thirdorder  = 123456789.0_flyt
        free_energy_fourthorder = 123456789.0_flyt
        return
    endif
    select type(qp); type is(lo_fft_mesh)
        call anharmonic_free_energy(p,fct,fcf,qp,dr,temperature,free_energy_thirdorder,free_energy_fourthorder,mw,mem,verbosity=-1)
    end select
end subroutine

!> phonon free energy at a single point
subroutine phonon_free_energy_for_single_point(gs,map,depvar,qgrid,temperature,mw,mem,fph)
    !> simulation grid
    type(lo_gridsim), intent(inout) :: gs
    !> forcemap
    type(lo_forcemap), intent(inout) :: map
    !> where to evaluate
    real(flyt), dimension(:), intent(in) :: depvar
    !> q-mesh density
    integer, dimension(3), intent(in) :: qgrid
    !> temperature
    real(flyt), intent(in) :: temperature
    !> mpi communicator
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem
    !> free energy
    real(flyt), intent(out) :: fph


    class(lo_qpoint_mesh), allocatable :: qp
    type(lo_phonon_dispersions) :: dr
    type(lo_crystalstructure) :: p
    type(lo_forceconstant_secondorder) :: fc
    real(flyt), dimension(:,:), allocatable :: pairconstraints
    integer :: nconstr

    ! setup stuffs, first get a structure
    call gs%structure%interpolate(depvar,p,gs%info%dim_volume)
    call p%classify('wedge',timereversal=.true.)
    ! constraints
    !call lo_secondorder_rot_herm_huang(map,p,pairconstraints,nconstr,.true.,.true.,.true.)
    nconstr=0
    ! evaluate forceconstants
    if ( nconstr .gt. 0 ) then
        call gs%eval(map,depvar,pairconstraints)
    else
        call gs%eval(map,depvar)
    endif
    ! forceconstant at this point
    call map%get_secondorder_forceconstant(p,fc,mem,-1)
    ! q-points
    call lo_generate_qmesh(qp,p,qgrid,'monkhorst',timereversal=.true.,headrankonly=.false.,mw=mw,mem=mem,verbosity=-1)
    ! phonon dispersions
    call dr%generate(qp,fc,p,mw=mw,mem=mem,verbosity=-1)
    ! and the free energy
    fph=dr%phonon_free_energy(temperature)
end subroutine

end module
