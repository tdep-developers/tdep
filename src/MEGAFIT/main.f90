#include "precompilerdefinitions"
program MEGAFIT
!!{!holding/MEGAFIT/manual.md!}
use konstanter, only: r8,lo_status,lo_exitcode_io,lo_tol,lo_volume_A_to_Bohr,lo_A_to_bohr,lo_bohr_to_A,lo_pressure_HartreeBohr_to_GPa
use gottochblandat, only: tochar,walltime,open_file,lo_determ,lo_mean,lo_chop,&
                   lo_does_file_exist
use mpi_wrappers, only: lo_mpi_helper,lo_stop_gracefully
use lo_memtracker, only: lo_mem_helper
use lo_timetracker, only: lo_timer
use type_crystalstructure, only: lo_crystalstructure
use type_forceconstant_firstorder,  only: lo_forceconstant_firstorder
use type_forceconstant_secondorder, only: lo_forceconstant_secondorder
use type_forceconstant_thirdorder,  only: lo_forceconstant_thirdorder
use type_forceconstant_fourthorder, only: lo_forceconstant_fourthorder
use lo_dielectric_interaction, only: lo_dielectric_tensor
!use type_jij_secondorder, only: lo_jij_secondorder
use lo_symmetry_of_interactions, only: lo_interaction_tensors
use type_forcemap, only: lo_forcemap,lo_secondorder_rot_herm_huang
use hdf5_wrappers, only: lo_h5_read_data,HID_T,H5F_ACC_TRUNC_F,H5f_ACC_RDONLY_F,&
                         h5close_f,h5open_f,h5fclose_f,h5fopen_f,h5fcreate_f,h5gclose_f,h5gopen_f,h5gcreate_f

use options, only: lo_opts
use type_gridsim, only: lo_gridsim
use gridenergy, only: lo_gridenergy
use diagnostics, only: get_diagnostics

implicit none
type(lo_opts) :: opts
type(lo_crystalstructure) :: uc,ss
type(lo_forcemap) :: map
type(lo_gridsim) :: gs
type(lo_mpi_helper) :: mw
type(lo_mem_helper) :: mem
type(lo_timer) :: tmr

real(r8) :: tstart
integer :: evalmode

! Figure out all the symmetry stuff.
init: block
    type(lo_interaction_tensors) :: slt
    real(r8) :: t0

    call mw%init()
    tstart=walltime()
    t0=tstart
    call opts%parse()
    call mem%init()
    call tmr%start()
    call tmr%tick()

    ! First things first: make sure all the input files make sense
    call check_input('infile.simulations','infile.evalpoints',evalmode,mw)

    if ( mw%talk .eqv. .false. ) opts%verbosity=-100
    if ( mw%talk ) write(*,*) '... running on ',tochar(mw%n),' MPI ranks'
    call uc%readfromfile('infile.ucposcar',verbosity=opts%verbosity)
    call uc%classify('wedge',timereversal=.true.)
    call ss%readfromfile('infile.ssposcar',verbosity=opts%verbosity)
    call ss%classify('supercell',uc)
    if ( mw%talk ) write(*,*) '... read structures'

    call tmr%tock('reading input files')

    ! I need a forcemap. For now, I assume it's identical across all simulations.
    ! Get all the symmetries
    call slt%generate(uc,ss,opts%cutoff2,opts%cutoff3,opts%cutoff4,opts%polar,mw,mem,opts%verbosity+1,&
             transposition=opts%transposition,spacegroup=opts%spacegroup,&
             wzdim=opts%wzdimensions,nj2=opts%njump2,nj3=opts%njump3,nj4=opts%njump4,&
             firstorder=.true.,dielcutoff2=opts%dielcutoff2,dielcutoff3=opts%dielcutoff3)
             !dielcutoff2=opts%dielcutoff2,dielcutoff3=opts%dielcutoff3)
    ! And create the map
    call map%generate(uc,ss,polarcorrectiontype=opts%polarcorrectiontype,st=slt,mw=mw,mem=mem,verbosity=opts%verbosity)

    ! Create the constraints
    t0=walltime()
    !call map%forceconstant_constraints(uc,opts%rotationalconstraints,opts%huanginvariance,opts%hermitian,opts%verbosity+10)

    call tmr%tock('determined symmetry')

    if ( mw%talk ) write(*,*) '... got constraints in ',tochar(walltime()-t0),'s'

    if ( mw%talk ) then
        write(*,*) 'Created all symmetry related things (',tochar(walltime()-t0),'s)'
    endif
end block init

! grab all the simulations
grabsims: block
    ! Set all heuristics
    call gs%init('infile.simulations',map,opts%order,opts%temperature_scale,opts%distance_scale,opts%pairfittype,opts%weighted,opts%verbosity+1,mw)
    call tmr%tock('initialized grid')
    ! Read all the simulations from files, and arrange it over MPI ranks.
    call gs%read_simulations_from_file('infile.simulations',mw,mem,opts%verbosity+1)
    call tmr%tock('read simulations')
    ! Interpolate the structure somehow
    call gs%create_structure_interpolation(map,uc)
    call tmr%tock('created structure interpolation')
end block grabsims

! solve for everything
solvestuff: block
    real(r8) :: tt0,t0
    tt0=walltime()

    ! Now I can start solving for things, first deal with the polar stuff.
    select case(map%polar)
    case(1)
        t0=walltime()
        call gs%create_polar_interpolation(map,mw,mem,opts%verbosity+1)
        if ( map%polarcorrectiontype .eq. 3 ) call gs%subtract_polar_forces(map,mw,mem,opts%verbosity+1)
        if ( mw%talk ) write(*,*) 'Fixed the polar things in ',tochar(walltime()-t0),'s'
        call tmr%tock('created polar interpolation interpolation')
    case(2)
        call gs%solve_dielectric(map,ss,mw,mem,opts%verbosity+1)
        call tmr%tock('created dielectric interactions')
    case default
        if ( mw%talk ) write(*,*) '... no polar interactions to worry about.'
    end select

    ! Then the possible magnetic things
    ! if ( gs%info%magnetic_pair_interactions ) then
    !     t0=walltime()
    !     call gs%solve_magnetic_onsite(map,mw,opts%verbosity+1)
    !     call gs%solve_magnetic_pair_gridfit(map,mw,opts%verbosity+1)
    !     call gs%solve_magnetic_pair_polyfit(map,mw,opts%verbosity+1)
    !     call gs%subtract_magnetic_forces(map,mw,opts%verbosity+1)
    !     if ( mw%talk ) write(*,*) 'Fixed magnetic things in ',tochar(walltime()-t0),'s'
    ! endif

    ! Then there are two variants: either polynomials or the adaptive polynomial gridfit thing.
    select case(gs%info%pairfittype)
    case(1)
        ! Solve for secondorder forceconstants.
        if ( map%have_fc_pair ) then
            t0=walltime()
            call gs%solve_secondorder(map,mw,opts%verbosity+1)
            call gs%subtract_secondorder_forces(map,mw,mem,opts%verbosity+1)
            if ( mw%talk ) write(*,*) 'Fixed the pair things in ',tochar(walltime()-t0),'s'
            call tmr%tock('created pair fit')
        endif

        ! Solve for thirdorder forceconstants
        if ( map%have_fc_triplet ) then
            t0=walltime()
            call gs%solve_thirdorder(map,mw,opts%verbosity+1)
            call gs%subtract_thirdorder_forces(map,mw,opts%verbosity+1)
            if ( mw%talk ) write(*,*) 'Fixed the triplet things in ',tochar(walltime()-t0),'s'
            call tmr%tock('created triplet fit')
        endif

        ! Solve for fourthorder forceconstants
        if ( map%have_fc_quartet ) then
            t0=walltime()
            call gs%solve_fourthorder(map,mw,opts%verbosity+1)
            call gs%subtract_fourthorder_forces(map,mw,opts%verbosity+1)
            if ( mw%talk ) write(*,*) 'Fixed the quartet things in ',tochar(walltime()-t0),'s'
            call tmr%tock('created quartet fit')
        endif
    case(2)
        call gs%solve_gridfit(map,mw,mem,opts%verbosity+1)
        call tmr%tock('created fit')
    end select

    ! Interpolate the structure somehow
    call gs%create_energy_interpolation(opts%temperature_scale,mw,opts%verbosity+1)
    call tmr%tock('created energy interpolation')

    ! Get some diagnostics
    if ( opts%diagnostics ) then
        call get_diagnostics(gs,map,mw,mem,opts%verbosity+1)
        call tmr%tock('diagnostics')
    endif

    ! And now the heavy stuff is done!
    if ( mw%talk ) then
        write(*,*) 'Got all relevant polynomials in ',tochar(walltime()-tt0),'s'
    endif
end block solvestuff

! Could be the case that I should dump everything on the input grid.
if ( opts%dumpinputgrid ) then
if ( mw%talk ) then
dumpinput: block
    type(lo_forceconstant_secondorder) :: fc
    type(lo_forceconstant_thirdorder) :: fct
    type(lo_forceconstant_fourthorder) :: fcf
    type(lo_dielectric_tensor) :: di
    type(lo_crystalstructure) :: p
    real(r8), dimension(:,:), allocatable :: pairconstraints
    integer :: i,nconstr

real(r8), dimension(gs%ndim) :: transformed_gridcoord
real(r8) :: f0,f1
integer :: j,k,u

    write(*,*) ''
    write(*,*) 'Dumping data on the input grid.'

    ! get stuffs at points
    do i=1,gs%nsim
        call gs%structure%interpolate(gs%grid_coordinates(:,i),p,gs%info%dim_volume)
        call p%classify('bravais')
        call p%writetofile('ongrid_uc_'//tochar(i),1)
        ! forceconstants, first the constraints
        !call lo_secondorder_rot_herm_huang( map,p,pairconstraints,nconstr,.true.,.true.,.true. )
        nconstr=0
        if ( nconstr .gt. 0 ) then
            call gs%eval(map,gs%grid_coordinates(:,i),pairconstraints)
        else
            call gs%eval(map,gs%grid_coordinates(:,i))
        endif
        if ( map%have_fc_pair ) then
            call map%get_secondorder_forceconstant(p,fc,mem,-1)
            call fc%writetofile(p,'ongrid_fc2_'//tochar(i))
        endif
        if ( map%have_fc_triplet ) then
            call map%get_thirdorder_forceconstant(p,fct)
            call fct%writetofile(p,'ongrid_fc3_'//tochar(i))
        endif
        if ( map%have_fc_quartet ) then
            call map%get_fourthorder_forceconstant(p,fcf)
            call fcf%writetofile(p,'ongrid_fc4_'//tochar(i))
        endif
        if ( map%polar .gt. 1 ) then
            call map%get_dielectric_tensors(uc,di)
            call di%writetofile(p,'ongrid_di_'//tochar(i))
            u=open_file('out','ongrid_loto_'//tochar(i))
                do j=1,3
                    write(u,*) di%eps_inf(:,j)
                enddo
                do k=1,p%na
                do j=1,3
                    write(u,*) di%Z_singlet(:,j,k) !/lo_bohr_to_A
                enddo
                enddo
            close(u)
        endif

        ! Also decent place to dump free energy?

        ! if ( gs%info%magnetic_pair_interactions ) then
        !     call map%get_secondorder_jij(p,jij)
        !     call jij%writetofile(p,'jij_'//tochar(i))
        ! endif
    enddo
end block dumpinput
endif
call tmr%tock('dump data on input grid')
endif

! Now that we have everything it's reasonable to do something with the interpolation.
! Perhaps calculate free energy on a grid?
if ( opts%evalenergy ) then
if ( evalmode .eq. 1 .or. evalmode .eq. 3 ) then
griden: block
    type(lo_gridenergy) :: ge
    call ge%generate(gs,map,opts%qgrid_harm,opts%qgrid_anharm,opts%quasiharmonic,opts%dumpfullgrid,mw,mem)
    call tmr%tock('evaluated fine mesh')
end block griden
endif
endif

! dump data at a list of provided points
if ( evalmode .eq. 2 ) then
if ( mw%talk ) then
dumppts: block
    type(lo_crystalstructure) :: p
    type(lo_forceconstant_secondorder) :: fc
    type(lo_forceconstant_thirdorder) :: fct
    type(lo_forceconstant_fourthorder) :: fcf
    type(lo_dielectric_tensor) :: di
    real(r8), dimension(:,:), allocatable :: gridcoord,pairconstraints
    integer :: u,i,j,npts,nconstr

    ! grab the points
    u=open_file('in','infile.evalpoints')
        read(u,*) j
        read(u,*) npts
        lo_allocate(gridcoord(gs%ndim,npts))
        do i=1,npts
            read(u,*) gridcoord(:,i)
        enddo
        i=gs%info%dim_volume
        if ( i .gt. 0 ) gridcoord(i,:)=gridcoord(i,:)*lo_volume_A_to_Bohr
    close(u)

    write(*,*) ''
    write(*,*) 'Evaluating on a set of '//tochar(npts)//' points'

    ! get stuffs at points
    do i=1,npts
        call gs%structure%interpolate(gridcoord(:,i),p,gs%info%dim_volume)
        call p%classify('bravais')
        call p%writetofile('uc_'//tochar(i),1)
        ! forceconstants, first the constraints
        !call lo_secondorder_rot_herm_huang( map,p,pairconstraints,nconstr,.true.,.true.,.true. )
        nconstr=0
        if ( nconstr .gt. 0 ) then
            call gs%eval(map,gridcoord(:,i),pairconstraints)
        else
            call gs%eval(map,gridcoord(:,i))
        endif
        if ( map%have_fc_pair ) then
            call map%get_secondorder_forceconstant(p,fc,mem,-1)
            call fc%writetofile(p,'fc2_'//tochar(i))
        endif
        if ( map%have_fc_triplet ) then
            call map%get_thirdorder_forceconstant(p,fct)
            call fct%writetofile(p,'fc3_'//tochar(i))
        endif
        if ( map%have_fc_quartet ) then
            call map%get_fourthorder_forceconstant(p,fcf)
            call fcf%writetofile(p,'fc4_'//tochar(i))
        endif
        if ( map%polar .gt. 1 ) then
            call map%get_dielectric_tensors(uc,di)
            call di%writetofile(p,'di_'//tochar(i))
        endif
        ! if ( gs%info%magnetic_pair_interactions ) then
        !     call map%get_secondorder_jij(p,jij)
        !     call jij%writetofile(p,'jij_'//tochar(i))
        ! endif
    enddo
end block dumppts
endif
    call tmr%tock('dump data on grid')
endif

call tmr%stop()
call tmr%dump(mw,'Timings:')
if ( mw%talk ) write(*,*) 'Done in ',tochar(walltime()-tstart),'s'

call mw%destroy()

contains

! I really need the code to die fast if we have bad input.
subroutine check_input(fn_sim,fn_eval,evalmode,mw)
    !> simulation input file
    character(len=*), intent(in) :: fn_sim
    !> evaluation points input file
    character(len=*), intent(in) :: fn_eval
    !> return how things are to be evaluated. Have to fix this later.
    integer, intent(out) :: evalmode
    !> MPI helper, to make sure I kill all instances
    type(lo_mpi_helper), intent(in) :: mw

    ! predefined options that are valid
    character(len=5), dimension(4), parameter :: ok_eos=['Birch','Vinet','2D-Bi','null ']
    character(len=3), dimension(3), parameter :: ok_var=['V  ','T  ','eta']
    character(len=3), dimension(3), parameter :: ok_spacing=['lin','log','den']
    integer, dimension(3), parameter :: ok_evalmode=[1,2,3]

    integer(HID_T) :: file_id
    real(r8), dimension(:,:), allocatable :: gridcoord,latticevectors,positions,evalpts
    real(r8), dimension(:), allocatable :: mincoord,maxcoord
    real(r8), dimension(4) :: eos4p
    real(r8), dimension(9) :: eos9p
    real(r8) :: pressurespacing
    character(len=3), dimension(:), allocatable :: varnames
    character(len=3) :: cspacing
    character(len=5) :: eosname
    character(len=2000) :: simfn
    integer, dimension(:), allocatable :: pts_per_dim,order_per_dim
    integer :: u,ndim,nsim,i,j,l,ctr,npts

    call h5open_f(lo_status)
    if ( lo_status .ne. 0 ) call lo_stop_gracefully(['Could initialize hdf5 library'],lo_exitcode_io,communicator=mw%comm)
    ! start with infile simulations
    ctr=0
    u=open_file('in','infile.simulations')
        ctr=ctr+1; read(u,*,iostat=lo_status) ndim
        if ( lo_status .ne. 0 ) call lo_stop_gracefully(['Could not read number of dimensions from "'//trim(fn_sim)//'"'],lo_exitcode_io,communicator=mw%comm)
        ctr=ctr+1; read(u,*,iostat=lo_status) nsim
        if ( lo_status .ne. 0 ) call lo_stop_gracefully(['Could not read number of simulations from "'//trim(fn_sim)//'"'],lo_exitcode_io,communicator=mw%comm)

        allocate(varnames(ndim))
        ctr=ctr+1; read(u,*,iostat=lo_status) varnames
        if ( lo_status .ne. 0 ) call lo_stop_gracefully(['Could not read variable names from "'//trim(fn_sim)//'"'],lo_exitcode_io,communicator=mw%comm)
        allocate(order_per_dim(ndim))
        ctr=ctr+1; read(u,*,iostat=lo_status) order_per_dim
        if ( lo_status .ne. 0 ) call lo_stop_gracefully(['Could not read polynomial orders from "'//trim(fn_sim)//'"'],lo_exitcode_io,communicator=mw%comm)
        do i=1,ndim
            if ( order_per_dim(i) .gt. 4 ) then
                call lo_stop_gracefully(['Polynomials of order 5 or higher is highly unreliable.'],lo_exitcode_io,communicator=mw%comm)
            elseif ( order_per_dim(i) .lt. 0 ) then
                call lo_stop_gracefully(['Polynomials orders need to be positiv.'],lo_exitcode_io,communicator=mw%comm)
            endif
        enddo

        ! Check with respect to the allowed character names, can not have more than 1 volume for example
        do i=1,size(ok_var)
            l=0
            do j=1,ndim
                if ( trim(lo_lowercase(varnames(j))) .eq. trim(lo_lowercase(ok_var(i))) ) l=l+1
            enddo
            if ( l .gt. 1 ) call lo_stop_gracefully(['You have to specify one label per dimension, and that label needs to be unique'],lo_exitcode_io,communicator=mw%comm)
        enddo
        ! Check that the variables specified are ones that I have thought about
        do i=1,ndim
            l=0
            do j=1,size(ok_var)
                if ( trim(lo_lowercase(varnames(i))) .eq. trim(lo_lowercase(ok_var(j))) ) l=l+1
            enddo
            if ( l .eq. 0 ) call lo_stop_gracefully(['The label "'//trim(varnames(i))//'" is not something I recognize. Consult the manual.'],lo_exitcode_io,communicator=mw%comm)
        enddo

        ! get the name of the equation of state
        ctr=ctr+1; read(u,*,iostat=lo_status) eosname
        if ( lo_status .ne. 0 ) call lo_stop_gracefully(['Could not read label of equation of state from "'//trim(fn_sim)//'"'],lo_exitcode_io,communicator=mw%comm)
        ! check that it's a real EOS
        l=0
        do i=1,size(ok_eos)
            if ( trim(lo_lowercase(ok_eos(i))) .eq. trim(lo_lowercase(eosname)) ) l=i
        enddo
        select case(l)
            case(1)
                ctr=ctr+1; read(u,*,iostat=lo_status) eos4p
                if ( lo_status .ne. 0 ) call lo_stop_gracefully(['Could not read 4 Birch-Murnaghan parameters from "'//trim(fn_sim)//'"'],lo_exitcode_io,communicator=mw%comm)
            case(2)
                ctr=ctr+1; read(u,*,iostat=lo_status) eos4p
                if ( lo_status .ne. 0 ) call lo_stop_gracefully(['Could not read 4 Vinet parameters from "'//trim(fn_sim)//'"'],lo_exitcode_io,communicator=mw%comm)
            case(3)
                ctr=ctr+1; read(u,*,iostat=lo_status) eos9p
                if ( lo_status .ne. 0 ) call lo_stop_gracefully(['Could not read 9 2-D Birch-Murnaghan parameters from "'//trim(fn_sim)//'"'],lo_exitcode_io,communicator=mw%comm)
            case(4)
                ! do nothing
        case default
            call lo_stop_gracefully(['The equation of state "'//trim(eosname)//'" is not recognized. Consult the manual.'],lo_exitcode_io,communicator=mw%comm)
        end select

        ! now go over each simulation, and make sure they exist
        lo_allocate(gridcoord(ndim,nsim))
        gridcoord=0.0_r8
        do i=1,nsim
            ctr=ctr+1; read(u,*,iostat=lo_status) gridcoord(:,i),simfn


            if ( lo_status .ne. 0 ) call lo_stop_gracefully(['Could not read grid coordinates or simulation filename at line '//tochar(ctr)//' from "'//trim(fn_sim)//'"'],lo_exitcode_io,communicator=mw%comm)
            ! convert to atomic units
            do j=1,ndim
                if ( trim(lo_lowercase(varnames(j))) .eq. 'v' ) gridcoord(j,i)=gridcoord(j,i)*lo_volume_A_to_Bohr
            enddo
            if ( lo_does_file_exist(trim(simfn)) .eqv. .false. ) then
                call lo_stop_gracefully(['Simulation file "'//trim(simfn)//'" does not exist.'],lo_exitcode_io,communicator=mw%comm)
            else
                ! if the file does exist, and one of my variables is volume, make sure the volume in the input file
                ! and that of the simulation match
                do j=1,ndim
                    if ( trim(lo_lowercase(varnames(j))) .ne. 'v' ) cycle
                    call h5fopen_f(trim(simfn), H5F_ACC_RDONLY_F, file_id, lo_status)
                    if ( lo_status .ne. 0 ) call lo_stop_gracefully(['Could not open "'//trim(simfn)//'"'],lo_exitcode_io,communicator=mw%comm)

                        call lo_h5_read_data(latticevectors  ,file_id,'unitcell_latticevectors',error=lo_status)
                        if ( lo_status .ne. 0 ) call lo_stop_gracefully(['Could read "unitcell_latticevectors" from "'//trim(simfn)//'"'],lo_exitcode_io,communicator=mw%comm)
                        latticevectors=latticevectors*lo_A_to_bohr
                        call lo_h5_read_data(positions       ,file_id,'unitcell_positions',error=lo_status)
                        if ( lo_status .ne. 0 ) call lo_stop_gracefully(['Could read "unitcell_positions" from "'//trim(simfn)//'"'],lo_exitcode_io,communicator=mw%comm)
                    call h5fclose_f(file_id, lo_status)
                    if ( lo_status .ne. 0 ) call lo_stop_gracefully(['Could not close "'//trim(simfn)//'"'],lo_exitcode_io,communicator=mw%comm)
                    if ( abs(gridcoord(j,i)-abs(lo_determ(latticevectors))/size(positions,2))/gridcoord(j,i) .gt. lo_tol ) then
                        call lo_stop_gracefully(['Volume on line '//tochar(ctr)//' of "'//trim(fn_sim)//'" and volume in "'//trim(simfn)//'" does not match.'],lo_exitcode_io,communicator=mw%comm)
                    endif
                enddo
            endif
        enddo
    close(u)
    call h5close_f(lo_status)
    if ( lo_status .ne. 0 ) call lo_stop_gracefully(['Could not close hdf5 library'],lo_exitcode_io,communicator=mw%comm)

    ! so far, so good. Now check the evaluation grid!
    ctr=0
    u=open_file('in',trim(fn_eval))
        ctr=ctr+1; read(u,*,iostat=lo_status) evalmode
        if ( lo_status .ne. 0 ) call lo_stop_gracefully(['Could not read evaluation mode from "'//trim(fn_eval)//'"'],&
                                     lo_exitcode_io,communicator=mw%comm)
        select case(evalmode)
        case(1) ! a full grid
            lo_allocate(pts_per_dim(ndim))
            lo_allocate(mincoord(ndim))
            lo_allocate(maxcoord(ndim))
            pts_per_dim=0
            mincoord=0.0_r8
            maxcoord=0.0_r8
            do i=1,ndim

                ! get number of points
                ctr=ctr+1; read(u,*,iostat=lo_status) pts_per_dim(i),cspacing
                if ( lo_status .ne. 0 ) call lo_stop_gracefully(&
                    ['Could read number of points and their spacing in dimension '//tochar(i)//' from "'//trim(fn_eval)//'"'],&
                    lo_exitcode_io,communicator=mw%comm)
                if ( pts_per_dim(i) .le. 4 ) call lo_stop_gracefully(&
                    ['Interpolation will be too strange if one of the dimensions has too few points.'],lo_exitcode_io,__FILE__,__LINE__)

                ! is this an acceptable spacing?
                select case(trim(lo_lowercase(varnames(i))))
                case('v')
                    l=0
                    if ( trim(lo_lowercase(cspacing)) .eq. 'lin' ) l=l+1
                    if ( trim(lo_lowercase(cspacing)) .eq. 'den' ) l=l+1
                    if ( l .eq. 0 ) then
                        call lo_stop_gracefully(['In volume allowed spacings are "lin" and "den".'],&
                                                 lo_exitcode_io,communicator=mw%comm)
                    endif
                case('t')
                    l=0
                    if ( trim(lo_lowercase(cspacing)) .eq. 'lin' ) l=l+1
                    if ( trim(lo_lowercase(cspacing)) .eq. 'log' ) l=l+1
                    if ( l .eq. 0 ) then
                        call lo_stop_gracefully(['In temperature allowed spacings are "lin" and "log".'],&
                                                 lo_exitcode_io,communicator=mw%comm)
                    endif
                case('eta')
                    if ( trim(lo_lowercase(cspacing)) .ne. 'lin' ) then
                        call lo_stop_gracefully(['In eta-direction only linear spacing is allowed.'],&
                                                lo_exitcode_io,communicator=mw%comm)
                    endif
                end select
                ! get min and max
                ctr=ctr+1; read(u,*,iostat=lo_status) mincoord(i),maxcoord(i)
                if ( lo_status .ne. 0 ) call lo_stop_gracefully(&
                    ['Could not read min and max in dimension '//tochar(i)//' from "'//trim(fn_eval)//'"'],&
                    lo_exitcode_io,communicator=mw%comm)
                ! Adjust unit to atomic units
                if ( trim(lo_lowercase(varnames(i))) .eq. 'v' ) then
                    mincoord(i)=mincoord(i)*lo_volume_A_to_bohr
                    maxcoord(i)=maxcoord(i)*lo_volume_A_to_bohr
                endif
                ! quick bounds check to make sure that people are not extrapolating.
                if ( minval(gridcoord(i,:)) .gt. mincoord(i)+lo_tol ) then
                    call lo_stop_gracefully(['Lowest value in dimension '//tochar(i)//' ('//trim(varnames(i))//') in "'//trim(fn_eval)//'" outside the grid.'],&
                                            lo_exitcode_io,communicator=mw%comm)
                endif
                if ( maxval(gridcoord(i,:)) .lt. maxcoord(i)-lo_tol ) then
                    call lo_stop_gracefully(['Largest value in dimension '//tochar(i)//' ('//trim(varnames(i))//') in "'//trim(fn_eval)//'" outside the grid.'],&
                         lo_exitcode_io,communicator=mw%comm)
                endif
            enddo
        case(2) ! just a list of points
            ctr=ctr+1; read(u,*,iostat=lo_status) npts
            if ( lo_status .ne. 0 ) call lo_stop_gracefully(['Could not read number of points from "'//trim(fn_eval)//'"'],&
                                     lo_exitcode_io,communicator=mw%comm)
            lo_allocate(evalpts(ndim,npts))
            evalpts=0.0_r8
            do i=1,npts
                read(u,*,iostat=lo_status) evalpts(:,i)

                if ( lo_status .ne. 0 ) then
                    call lo_stop_gracefully(['Could not read coordinates of point '//tochar(i)//' "'//trim(fn_eval)//'"'],&
                                            lo_exitcode_io,communicator=mw%comm)
                endif
                ! Adjust unit to atomic units
                do j=1,ndim
                    if ( trim(lo_lowercase(varnames(j))) .eq. 'v' ) then
                        evalpts(j,i)=evalpts(j,i)*lo_volume_A_to_bohr
                    endif
                enddo
                ! check bounds
                do j=1,ndim
                if ( minval(gridcoord(j,:)) .gt. evalpts(j,i)+lo_tol ) then
                    call lo_stop_gracefully(['Point '//tochar(i)//' is outside the grid in dimension '//tochar(j)],&
                                             lo_exitcode_io,communicator=mw%comm)
                endif
                if ( maxval(gridcoord(j,:)) .lt. evalpts(j,i)-lo_tol ) then
                    call lo_stop_gracefully(['Point '//tochar(i)//' is outside the grid in dimension '//tochar(j)],&
                                             lo_exitcode_io,communicator=mw%comm)
                endif
                enddo
            enddo
        case(3) ! Specialized grid where temperature is specified but volume is determined by the convex hull of the grid data
            lo_allocate(pts_per_dim(ndim))
            lo_allocate(mincoord(ndim))
            lo_allocate(maxcoord(ndim))
            pts_per_dim=0
            mincoord=0.0_r8
            maxcoord=0.0_r8
            ! have to make sure that it's two-dimensional (at least for now)
            if ( ndim .ne. 2 ) call lo_stop_gracefully(['Can only use adaptive mesh for two dimensions'],&
                                    lo_exitcode_io,communicator=mw%comm)
            ! the variables have to be temperature and volume
            l=0
            if ( trim(varnames(1)) .eq. 'V' .or. trim(varnames(2)) .eq. 'V' ) l=l+1
            if ( trim(varnames(1)) .eq. 'T' .or. trim(varnames(2)) .eq. 'T' ) l=l+1
            if ( l .ne. 2 ) call lo_stop_gracefully(['The dimensions have to be "V" and "T".'],&
                                 lo_exitcode_io,communicator=mw%comm)

            i=0
            if ( trim(varnames(1)) .eq. 'T' ) i=1
            if ( trim(varnames(2)) .eq. 'T' ) i=2
            ! grab the temperature and make sure it's an allowed spacing
            ctr=ctr+1; read(u,*,iostat=lo_status) pts_per_dim(i),cspacing
            if ( lo_status .ne. 0 ) call lo_stop_gracefully(&
                ['Could read number of points and their spacing in temperature from "'//trim(fn_eval)//'"'],&
                lo_exitcode_io,communicator=mw%comm)
            if ( pts_per_dim(i) .le. 4 ) call lo_stop_gracefully(&
                ['Interpolation will be too strange if one of the dimensions has too few points.'],lo_exitcode_io,__FILE__,__LINE__)
            l=0
            if ( trim(lo_lowercase(cspacing)) .eq. 'lin' ) l=l+1
            if ( trim(lo_lowercase(cspacing)) .eq. 'log' ) l=l+1
            if ( l .eq. 0 ) then
                call lo_stop_gracefully(['In temperature allowed spacings are "lin" and "log".'],&
                                         lo_exitcode_io,communicator=mw%comm)
            endif
            ! check bounds of temperature?
            ctr=ctr+1; read(u,*,iostat=lo_status) mincoord(i),maxcoord(i)
            if ( lo_status .ne. 0 ) call lo_stop_gracefully(&
                ['Could not read min and max in temperature '//tochar(i)//' from "'//trim(fn_eval)//'"'],&
                lo_exitcode_io,communicator=mw%comm)
            ! quick bounds check to make sure that people are not extrapolating.
            if ( minval(gridcoord(i,:)) .gt. mincoord(i)+lo_tol ) then
                call lo_stop_gracefully(['Lowest value in temperature in "'//trim(fn_eval)//'" outside the grid.'],&
                                        lo_exitcode_io,communicator=mw%comm)
            endif
            if ( maxval(gridcoord(i,:)) .lt. maxcoord(i)-lo_tol ) then
                call lo_stop_gracefully(['Largest value in temperature "'//trim(fn_eval)//'" outside the grid.'],&
                     lo_exitcode_io,communicator=mw%comm)
            endif
            i=0
            if ( trim(varnames(1)) .eq. 'V' ) i=1
            if ( trim(varnames(2)) .eq. 'V' ) i=2
            ! grab the number of volume points and make sure it's an allowed spacing
            ctr=ctr+1; read(u,*,iostat=lo_status) pts_per_dim(i),cspacing
            if ( pts_per_dim(i) .le. 4 ) call lo_stop_gracefully(&
                ['Interpolation will be too strange if one of the dimensions has too few points.'],lo_exitcode_io,__FILE__,__LINE__)
            if ( lo_status .ne. 0 ) call lo_stop_gracefully(&
                ['Could read number of points and their spacing in volume from "'//trim(fn_eval)//'"'],&
                lo_exitcode_io,communicator=mw%comm)
            l=0
            if ( trim(lo_lowercase(cspacing)) .eq. 'lin' ) l=l+1
            if ( trim(lo_lowercase(cspacing)) .eq. 'den' ) l=l+1
            if ( l .eq. 0 ) then
                call lo_stop_gracefully(['In volume allowed spacings are "lin" and "den".'],&
                                         lo_exitcode_io,communicator=mw%comm)
            endif
            ctr=ctr+1; read(u,*,iostat=lo_status) pressurespacing
            if ( lo_status .ne. 0 ) call lo_stop_gracefully(&
                ['Could read the spacing in pressure from "'//trim(fn_eval)//'"'],&
                lo_exitcode_io,communicator=mw%comm)
        case default
        end select
    close(u)
end subroutine

!> convert a string to all lowercase
pure function lo_lowercase(string) result(lowercasestring)
    !> string to convert
    character(len=*), intent(in) :: string
    !> lowercase string
    character(len=:), allocatable :: lowercasestring

    character (len=26), parameter :: L = "abcdefghijklmnopqrstuvwxyz"
    character (len=26), parameter :: U = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    character(len=1) :: a
    integer :: i,j

    allocate(lowercasestring,source=string)
    do i=1,len(lowercasestring)
        a=lowercasestring(i:i)
        j=index(U,a)
        if ( j .gt. 0 ) a=L(j:j)
        lowercasestring(i:i)=a
    enddo
end function

end program
