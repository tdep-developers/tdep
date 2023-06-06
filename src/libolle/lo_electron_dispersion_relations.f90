module lo_electron_dispersion_relations
!!
!! Handle electron energies on a mesh in the Brilloin zone.
!!
use konstanter, only: r8,lo_iou,lo_huge,lo_hugeint,lo_status,lo_tol,lo_sqtol,lo_Hartree_to_eV,lo_eV_to_Hartree,&
    lo_exitcode_param,lo_kb_hartree
use gottochblandat, only: open_file,walltime,tochar,lo_sqnorm,lo_chop,lo_trueNtimes,lo_fermi,&
                   lo_progressbar_init,lo_progressbar
use lo_sorting, only: lo_qsort
use lo_memtracker, only: lo_mem_helper
use lo_brents_method, only: lo_brent_helper
use mpi_wrappers, only: lo_mpi_helper,lo_stop_gracefully
use type_qpointmesh, only: lo_qpoint_mesh,lo_LV_tetrahedron_heaviside,lo_LV_tetrahedron_fermi
use type_effective_hamiltonian, only: lo_effective_hamiltonian
implicit none

private
! types
public :: lo_electron_dispersions_kpoint
public :: lo_electron_dispersions

!> A k-point in dispersion relations
type lo_electron_dispersions_kpoint
    !> KS eigenvalues in eV, per spin channel
    real(r8), dimension(:,:), allocatable :: eigenvalue
    !> KS eigenvectors
    complex(r8), dimension(:,:,:), allocatable :: eigenvector
    !> Electron group velocity
    real(r8), dimension(:,:,:), allocatable :: groupvelocity
    !> Linewidth due to electron-phonon coupling
    real(r8), dimension(:,:), allocatable :: electronphononlinewidth
    !> Electron phonon phasespace
    real(r8), dimension(:,:), allocatable :: electronphononphasespace
    contains
        !> measure size in memory, in bytes
        procedure :: size_in_mem=>kpoint_size_in_mem
end type

!> Electron dispersion relations in the full BZ
type lo_electron_dispersions
    !> how many bands
    integer :: n_band=-lo_hugeint
    !> how many basis functions
    integer :: n_basis=lo_hugeint
    !> how many spin channels
    integer :: n_spin=-lo_hugeint
    !> how many electrons
    real(r8) :: n_electron=-lo_huge
    !> the Fermi energy (or highest occupied state)
    real(r8) :: efermi=-lo_huge
    !> approximate band gap
    real(r8) :: bandgap=-lo_huge
    !> approximate energy where the valence starts
    real(r8) :: evalence=-lo_huge
    !> number of points in the irreducible part
    integer :: n_irr_kpoint=-lo_hugeint
    !> total number of points
    integer :: n_full_kpoint=-lo_hugeint
    !> k-points in the irreducible wedge
    type(lo_electron_dispersions_kpoint), dimension(:), allocatable :: ik
    !> k-points in the full zone. Maybe don't use these.
    ! type(lo_electron_dispersions_kpoint), dimension(:), allocatable :: ak
    !> max energy
    real(r8) :: energy_max=-lo_huge
    !> min energy
    real(r8) :: energy_min=-lo_huge
    !> min energy per band
    real(r8), dimension(:), allocatable :: bandmin
    !> max energy for each band
    real(r8), dimension(:), allocatable :: bandmax
    !> default smearing for each band (iband,ispin)
    real(r8), dimension(:,:), allocatable :: default_smearing
    contains
        !> figure out some basic things
        procedure :: classify
        !> generate via an effective Hamiltonian
        procedure :: generate
        !> measure size in memory, in bytes
        procedure :: size_in_mem=>edr_size_in_mem
        !> destroy
        procedure :: destroy
end type

contains

!> Get electron energies on a mesh
subroutine generate(edr,eh,kp,fermi_smearing_type,sigma,temperature,eigenvectors,groupvelocities,mw,mem,verbosity)
    !> electron dispersions
    class(lo_electron_dispersions), intent(out) :: edr
    !> effective Hamiltonian
    type(lo_effective_hamiltonian), intent(in) :: eh
    !> distributed k-point mesh
    class(lo_qpoint_mesh), intent(in) :: kp
    !> how to smear the occupation numbers
    integer, intent(in) :: fermi_smearing_type
    !> smearing parameter, if applicable
    real(r8), intent(in) :: sigma
    !> temperature, if applicable
    real(r8), intent(in) :: temperature
    !> keep the eigenvectors?
    logical, intent(in) :: eigenvectors
    !> keep group velocities
    logical, intent(in) :: groupvelocities
    !> MPI helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> Memory tracker
    type(lo_mem_helper), intent(inout) :: mem
    !> talk?
    integer, intent(in) :: verbosity

    real(r8) :: timer,t0,t1
    integer :: slvstrat

    ! Start timers
    timer=walltime()
    t0=timer
    t1=timer

    ! Figure out some basic things
    init: block
        integer :: i
        if ( verbosity .gt. 0 ) then
            write(lo_iou,*) ''
            write(lo_iou,*) 'Calculating ks eigenvalues on k-point grid (',tochar(mw%n),' ranks)'
            timer=walltime()
        endif

        edr%n_band=eh%n_band
        edr%n_basis=eh%n_basis
        edr%n_electron=eh%n_electron
        edr%n_spin=eh%n_spin
        edr%efermi=-lo_huge
        edr%evalence=-lo_huge
        edr%n_irr_kpoint=kp%n_irr_point
        edr%n_full_kpoint=-1 ! Disabled the full thing for now, takes too much space. Have to think.
        edr%energy_max=-lo_huge
        edr%energy_min=lo_huge

        ! Figure out how to solve
        slvstrat=0
        if ( eigenvectors ) then
            slvstrat=slvstrat+1
        endif
        if ( groupvelocities ) then
            slvstrat=slvstrat+2
        elseif ( fermi_smearing_type .eq. 5 ) then
            slvstrat=slvstrat+2
        endif

        ! Make some space for the solution
        allocate(edr%ik(edr%n_irr_kpoint))
        do i=1,edr%n_irr_kpoint
            allocate(edr%ik(i)%eigenvalue(edr%n_band,edr%n_spin))
            edr%ik(i)%eigenvalue=0.0_r8

            select case(slvstrat)
            case(0)
                ! nothing extra
            case(1)
                ! only eigenvectors
                allocate(edr%ik(i)%eigenvector(edr%n_basis,edr%n_band,edr%n_spin))
                edr%ik(i)%eigenvector=0.0_r8
            case(2)
                ! only groupvelocities
                allocate(edr%ik(i)%groupvelocity(3,edr%n_band,edr%n_spin))
                edr%ik(i)%groupvelocity=0.0_r8
            case(3)
                ! groupvelocities and eigenvectors
                allocate(edr%ik(i)%eigenvector(edr%n_basis,edr%n_band,edr%n_spin))
                allocate(edr%ik(i)%groupvelocity(3,edr%n_band,edr%n_spin))
                edr%ik(i)%eigenvector=0.0_r8
                edr%ik(i)%groupvelocity=0.0_r8
            end select
        enddo

        ! And space for aux things
        allocate(edr%bandmin(edr%n_band))
        allocate(edr%bandmax(edr%n_band))
        allocate(edr%default_smearing(edr%n_band,edr%n_spin))
        edr%bandmin=0.0_r8
        edr%bandmax=0.0_r8
        edr%default_smearing=0.0_r8
    end block init

    ! Solve the whole thing
    solve: block
        complex(r8), dimension(:,:,:,:), allocatable :: egvbuf
        real(r8), dimension(:,:,:,:), allocatable :: egvelbuf
        real(r8), dimension(:,:,:), allocatable :: evbuf
        real(r8) :: t0
        integer :: i,ispin,ctr
        t0=walltime()

        if ( verbosity .gt. 0 ) call lo_progressbar_init()

        select case(slvstrat)
        case(0)
            call mem%allocate(evbuf,[edr%n_band,edr%n_spin,edr%n_irr_kpoint],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
            evbuf=0.0_r8
        case(1)
            call mem%allocate(evbuf,[edr%n_band,edr%n_spin,edr%n_irr_kpoint],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
            call mem%allocate(egvbuf,[edr%n_basis,edr%n_band,edr%n_spin,edr%n_irr_kpoint],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
            evbuf=0.0_r8
            egvbuf=0.0_r8
        case(2)
            call mem%allocate(evbuf,[edr%n_band,edr%n_spin,edr%n_irr_kpoint],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
            call mem%allocate(egvelbuf,[3,edr%n_band,edr%n_spin,edr%n_irr_kpoint],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
            evbuf=0.0_r8
            egvelbuf=0.0_r8
        case(3)
            call mem%allocate(evbuf,[edr%n_band,edr%n_spin,edr%n_irr_kpoint],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
            call mem%allocate(egvbuf,[edr%n_basis,edr%n_band,edr%n_spin,edr%n_irr_kpoint],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
            call mem%allocate(egvelbuf,[3,edr%n_band,edr%n_spin,edr%n_irr_kpoint],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
            evbuf=0.0_r8
            egvbuf=0.0_r8
            egvelbuf=0.0_r8
        end select

        ctr=0
        do i=1,edr%n_irr_kpoint
            do ispin=1,edr%n_spin
                ctr=ctr+1
                if ( mod(ctr,mw%n) .ne. mw%r ) cycle
                select case(slvstrat)
                case(0)
                    call eh%solve_eigenproblem( kp%ip(i),evbuf(:,:,i),mem )
                case(1)
                    call eh%solve_eigenproblem( kp%ip(i),evbuf(:,:,i),mem, eigenvectors=egvbuf(:,:,:,i) )
                case(2)
                    call eh%solve_eigenproblem( kp%ip(i),evbuf(:,:,i),mem, groupvelocities=egvelbuf(:,:,:,i) )
                case(3)
                    call eh%solve_eigenproblem( kp%ip(i),evbuf(:,:,i),mem, eigenvectors=egvbuf(:,:,:,i),groupvelocities=egvelbuf(:,:,:,i) )
                end select
            enddo
            if ( verbosity .gt. 0 .and. lo_trueNtimes(i,5,edr%n_irr_kpoint) ) then
                call lo_progressbar(' ... solving KS eigenproblem',i,edr%n_irr_kpoint,walltime()-t0)
            endif
        enddo

        ! Sync things!
        select case(slvstrat)
        case(0)
            call mw%allreduce('sum',evbuf)
        case(1)
            call mw%allreduce('sum',evbuf)
            call mw%allreduce('sum',egvbuf)
        case(2)
            call mw%allreduce('sum',evbuf)
            call mw%allreduce('sum',egvelbuf)
        case(3)
            call mw%allreduce('sum',evbuf)
            call mw%allreduce('sum',egvbuf)
            call mw%allreduce('sum',egvelbuf)
        end select

        ! Copy it everywhere
        do i=1,edr%n_irr_kpoint
        do ispin=1,edr%n_spin
            select case(slvstrat)
            case(0)
                edr%ik(i)%eigenvalue=evbuf(:,:,i)
            case(1)
                edr%ik(i)%eigenvalue=evbuf(:,:,i)
                edr%ik(i)%eigenvector=egvbuf(:,:,:,i)
            case(2)
                edr%ik(i)%eigenvalue=evbuf(:,:,i)
                edr%ik(i)%groupvelocity=egvelbuf(:,:,:,i)
            case(3)
                edr%ik(i)%eigenvalue=evbuf(:,:,i)
                edr%ik(i)%eigenvector=egvbuf(:,:,:,i)
                edr%ik(i)%groupvelocity=egvelbuf(:,:,:,i)
            end select
        enddo
        enddo

        ! And some cleanup
        select case(slvstrat)
        case(0)
            call mem%deallocate(evbuf,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        case(1)
            call mem%deallocate(evbuf,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
            call mem%deallocate(egvbuf,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        case(2)
            call mem%deallocate(evbuf,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
            call mem%deallocate(egvelbuf,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        case(3)
            call mem%deallocate(evbuf,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
            call mem%deallocate(egvbuf,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
            call mem%deallocate(egvelbuf,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        end select

        if ( verbosity .gt. 0 ) call lo_progressbar(' ... solving KS eigenproblem',edr%n_irr_kpoint,edr%n_irr_kpoint,walltime()-t0)
    end block solve

    ! And get some info about the alignment, bandgap and stuff like that
    call edr%classify(kp,fermi_smearing_type,sigma,temperature,mw,mem,verbosity)

    ! And we are done!
    if ( verbosity .gt. 0 ) then
        write(lo_iou,*) 'Done calculating eigenvalues on grid (',tochar(walltime()-timer),'s)'
    endif
end subroutine

!> sort out some basic things about the electronic structure
subroutine classify(edr,kp,integrationtype,sigma,temperature,mw,mem,verbosity)
    !> electron dispersions
    class(lo_electron_dispersions), intent(inout) :: edr
    !> k-mesh
    class(lo_qpoint_mesh), intent(in) :: kp
    !> how to integrate
    integer, intent(in) :: integrationtype
    !> smearing parameter, if applicable
    real(r8), intent(in) :: sigma
    !> temperature, if applicable
    real(r8), intent(in) :: temperature
    !> MPI helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem
    !> talk?
    integer, intent(in) :: verbosity

    ! Some rough constants I just made up
    real(r8), parameter :: valenceseparation=4.0_r8*lo_eV_to_Hartree
    real(r8) :: t0,t1

    t0=walltime()
    t1=t0

    ! First get the largest and smallest energies, bandwidths and things like that
    minmax: block
        integer :: i,j,k

        ! Get the largest and smallest energies per band
        edr%bandmin=lo_huge
        edr%bandmax=-lo_huge
        do i=1,edr%n_irr_kpoint
            if ( mod(i,mw%n) .ne. mw%r ) cycle
            do k=1,edr%n_spin
            do j=1,edr%n_band
                edr%bandmin(j)=min(edr%bandmin(j),edr%ik(i)%eigenvalue(j,k))
                edr%bandmax(j)=max(edr%bandmax(j),edr%ik(i)%eigenvalue(j,k))
            enddo
            enddo
        enddo
        call mw%allreduce('min',edr%bandmin)
        call mw%allreduce('max',edr%bandmax)
        edr%energy_min=minval(edr%bandmin)
        edr%energy_max=maxval(edr%bandmax)
        if ( verbosity .gt. 0 ) write(lo_iou,*) '... determined bandwidths'
    end block minmax

    ! Determine default smearing parameter, per band.
    defaultsmearing: block
        real(r8), dimension(:), allocatable :: x
        integer :: ispin,iband,ikp,i
        integer :: solrnk
        character(len=1000) :: opf

        ! Which rank to sort on
        solrnk=mw%n-1
        ! dummy space
        call mem%allocate(x,kp%n_irr_point,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        x=0.0_r8

        edr%default_smearing=0.0_r8
        do ispin=1,edr%n_spin
        do iband=1,edr%n_band
            ! Collect values
            x=0.0_r8
            do ikp=1,kp%n_irr_point
                if ( mod(ikp,mw%n) .ne. mw%r ) cycle
                x(ikp)=edr%ik(ikp)%eigenvalue(iband,ispin)
            enddo
            call mw%allreduce('sum',x)
            if ( mw%r .eq. solrnk ) then
                call lo_qsort(x)
                do i=1,size(x)-1
                    edr%default_smearing(iband,ispin)=max(edr%default_smearing(iband,ispin),(x(i+1)-x(i))*0.5_r8)
                enddo
            endif
        enddo
        enddo
        call mem%deallocate(x,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        ! Sync across ranks
        call mw%bcast(edr%default_smearing,from=solrnk)

        ! And report
        if ( verbosity .gt. 0 ) then
            t1=walltime()
            write(lo_iou,*) '... default smearing (',tochar(t1-t0),'s)'
            opf="(1X,A14,"//tochar(edr%n_spin)//"(1X,F13.6),' (meV)')"
            do iband=1,edr%n_band
                if ( edr%bandmax(iband) .lt. edr%efermi-2*lo_eV_to_Hartree ) cycle
                if ( edr%bandmin(iband) .gt. edr%efermi+2*lo_eV_to_Hartree ) cycle
                write(lo_iou,opf) 'band '//tochar(iband)//': ',edr%default_smearing(iband,:)*lo_Hartree_to_eV*1000
            enddo
            t0=walltime()
        endif
    end block defaultsmearing

    ! Figure out wether there is a bandgap or not
    bngap: block
        real(r8) :: f0,f1,efermi,evalence,stepsize
        integer :: i,j,ctr,iter
        logical :: bgap

        ! First figure out if there is a bandgap. If there is, I am fairly certain that one of
        ! bandwidth maxima must give ebt%xloctly the correct number of electrons!
        f1=lo_huge
        j=-1
        do i=1,edr%n_band
            f0=roughly_count_electrons(edr,kp,edr%bandmax(i),mw)
            if ( abs(f0-edr%n_electron) .lt. f1 ) then
                j=i
                f1=abs(f0-edr%n_electron)
            endif
        enddo
        if ( abs( roughly_count_electrons(edr,kp,edr%bandmax(j),mw) - edr%n_electron ) .lt. lo_tol ) then
            ! might be a bandgap here!
            bgap=.true.
            efermi=edr%bandmax(j)
            ! are there any bandwidths that straddle the fermi level, it does not have a gap.
            do i=1,edr%n_band
                f0=lo_chop(edr%bandmin(i)-efermi,lo_tol)
                f1=lo_chop(edr%bandmax(i)-efermi,lo_tol)
                if ( lo_chop(f0*f1,lo_sqtol) .lt. -lo_sqtol ) then
                    bgap=.false.
                    exit
                endif
            enddo
        else
            ! probably not a bandgap.
            bgap=.false.
        endif

        ! Now get a decent estimate of the Fermi level:
        if ( bgap ) then
            ! should be easy enough to find the gap now:
            f0=lo_huge
            do i=1,edr%n_band
                if ( edr%bandmin(i)+lo_sqtol .lt. efermi ) cycle
                f0=min(f0,edr%bandmin(i)-efermi)
            enddo
            ! store bandgap and fermi energy
            edr%bandgap=f0
            edr%efermi=efermi
            if ( verbosity .gt. 0 ) write(lo_iou,*) '... found a band gap: ',tochar(edr%bandgap*lo_Hartree_to_eV),&
                ' eV, hos: ',tochar(edr%efermi*lo_Hartree_to_eV),' eV'
        else
            ! Get a rough fermi level
            stepsize=(maxval(edr%bandmax)-minval(edr%bandmin))/20.0_r8
            efermi=minval(edr%bandmin)-stepsize*0.1_r8
            ctr=0
            do iter=1,5000
                ! count electrons
                f1=roughly_count_electrons(edr,kp,efermi,mw)
                ! decide on something depending on how many there where
                if ( abs(f1-edr%n_electron) .lt. lo_tol ) then
                    ctr=ctr+1
                    if ( ctr .lt. 50 ) then
                        efermi=efermi-stepsize
                        stepsize=stepsize*0.5_r8
                    else
                        exit
                    endif
                endif
                ! Determine what to do: if too small, make another step. If too small, make a backwards step
                ! and decrease the stepsize.
                if ( f1 .lt. edr%n_electron ) then
                    efermi=efermi+stepsize
                elseif ( f1 .gt. edr%n_electron ) then
                    efermi=efermi-stepsize
                    stepsize=stepsize*0.5_r8
                endif
            enddo
            if ( verbosity .gt. 0 ) write(lo_iou,*) '... got a rough fermi level: ',tochar(efermi*lo_Hartree_to_eV),' eV'
            edr%efermi=efermi
            edr%bandgap=0.0_r8
        endif

        ! Now that I know approximately where the highest occupied state is, try to figure out there the
        ! valence starts. Mostly used for density of states plots later or something. How do I define the valence?
        ! I just start at the highest occupied state, step down to find a spot where there are no states 4eV above
        ! or below. Stupid, but surprisingly robust.
        evalence=edr%efermi
        stepsize=0.25_r8*lo_eV_to_Hartree
        ctr=0
        do iter=1,10000
            ! decrease
            f0=roughly_count_electrons(edr,kp,evalence-valenceseparation,mw)
            f1=roughly_count_electrons(edr,kp,evalence+valenceseparation,mw)
            if ( abs(f0-f1) .gt. lo_sqtol ) then
                evalence=evalence-stepsize
            else
                exit
            endif
        enddo
        edr%evalence=evalence
    end block bngap

    ! Get the Fermi level
    efermi: block
        real(r8), parameter :: eftol=1E-10_r8
        integer, parameter :: maxiter=10000
        type(lo_brent_helper) :: bt
        real(r8), dimension(:,:,:,:), allocatable :: tetenergies
        real(r8) :: x0,x1,y0,y1,x,y
        integer :: iter,iband,ispin,itet,ikp,icrn

        select case(integrationtype)
        case(3:4)
            ! we are going to tetrahedron integrate, prefetch some energies
            call mem%allocate(tetenergies,[4,kp%n_irr_tet,edr%n_spin,edr%n_band],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
            tetenergies=0.0_r8
            do iband=1,edr%n_band
            do ispin=1,edr%n_spin
            do itet=1,kp%n_irr_tet
            do icrn=1,4
                ikp=kp%it(itet)%irreducible_index(icrn)
                tetenergies(icrn,itet,ispin,iband)=edr%ik(ikp)%eigenvalue(iband,ispin)
            enddo
            enddo
            enddo
            enddo
        end select

        ! Find safe bounds for brents method. These are really really safe.
        x0=edr%energy_min-1.0_r8
        x1=edr%energy_max+1.0_r8
        y0=sensibly_count_electrons(edr,kp,x0,integrationtype,sigma,temperature,tetenergies,mw)-edr%n_electron
        y1=sensibly_count_electrons(edr,kp,x1,integrationtype,sigma,temperature,tetenergies,mw)-edr%n_electron
        ! Initialize Brent's method
        call bt%init( x0,x1,y0,y1,eftol )
        ! Iterate through to find the Fermi level
        do iter=1,maxiter
            !if ( mw%talk ) write(*,*) iter,bt%current_x()*lo_Hartree_to_eV,bt%current_y()
            ! check for convergence
            if ( abs(bt%current_y()) .lt. eftol ) then
                edr%efermi=bt%current_x()
                exit
            endif
            ! pick next x-value
            x=bt%next_x()
            ! evaluate function
            y=sensibly_count_electrons(edr,kp,x,integrationtype,sigma,temperature,tetenergies,mw)-edr%n_electron
            ! tell the solver what happened
            call bt%update(x,y)
        enddo

        ! And cleanup, possibly
        select case(integrationtype)
        case(3:4)
            ! we are going to tetrahedron integrate, prefetch some energies
            call mem%deallocate(tetenergies,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        end select

    end block efermi
end subroutine

!> shorthand to count electrons in a rather rough way
function roughly_count_electrons(edr,kp,efermi,mw) result(nelec)
    !> electron dispersions
    type(lo_electron_dispersions), intent(in) :: edr
    !> k-points
    class(lo_qpoint_mesh), intent(in) :: kp
    !> approximate fermi level
    real(r8), intent(in) :: efermi
    !> mpi helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> how many electrons for this fermi level?
    real(r8) :: nelec

    integer :: i,j,k

    nelec=0.0_r8
    do j=1,edr%n_band
        ! make it parallel, not really for speed but to make sure things are synced.
        if ( mod(j,mw%n) .ne. mw%r ) cycle
        ! check if it's a relevant band
        if ( edr%bandmax(j) .lt. efermi+lo_sqtol ) then
            ! this band is completely below
            select case(edr%n_spin)
            case(1)
                nelec=nelec+2.0_r8
            case(2)
                nelec=nelec+1.0_r8
            end select
        elseif ( edr%bandmin(j) .gt. efermi+lo_sqtol ) then
            ! completely above
            cycle
        else
            ! somewhere in between
            do k=1,edr%n_spin
            do i=1,edr%n_irr_kpoint
                if ( edr%ik(i)%eigenvalue(j,k) .lt. efermi+lo_sqtol ) then
                    select case(edr%n_spin)
                    case(1)
                        nelec=nelec+kp%ip(i)%integration_weight*2.0_r8
                    case(2)
                        nelec=nelec+kp%ip(i)%integration_weight
                    end select
                endif
            enddo
            enddo
        endif
    enddo
    call mw%allreduce('sum',nelec)
end function

!> count electrons in a somewhat accurate way.
function sensibly_count_electrons(edr,kp,efermi,it,sigma,temperature,tetenergies,mw) result(nelec)
    !> electron dispersions
    class(lo_electron_dispersions), intent(inout) :: edr
    !> k-mesh
    class(lo_qpoint_mesh), intent(in) :: kp
    !> approximate fermi level
    real(r8), intent(in) :: efermi
    !> integration type
    integer, intent(in) :: it
    !> smearing parameter
    real(r8), intent(in) :: sigma
    !> temperature
    real(r8), intent(in) :: temperature
    !> energies at tetrahedron corners, if needed
    real(r8), dimension(:,:,:,:), intent(in) :: tetenergies
    !> MPI helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> how many electrons for this fermi level?
    real(r8) :: nelec

    select case(it)
    case(1)
    ! The usual Gaussian smearing thing, nothing fancy.
    gaussiansmearing: block
        real(r8) :: sixsigma,isig,f0,f1,f2
        integer :: iband,ikp,ispin,ctr

        sixsigma=6*maxval(edr%default_smearing)
        ctr=0
        ! Count the easy ones
        bandloop1: do iband=1,edr%n_band
            ! check if it's a relevant band
            if ( edr%bandmax(iband) .lt. efermi-sixsigma ) then
                ! this band is completely below
                select case(edr%n_spin)
                case(1)
                    ctr=ctr+2
                case(2)
                    ctr=ctr+1
                end select
            endif
        enddo bandloop1

        nelec=0.0_r8
        spinloop: do ispin=1,edr%n_spin
        bandloop2: do iband=1,edr%n_band
            if ( edr%bandmax(iband) .lt. efermi-sixsigma ) cycle bandloop2
            if ( edr%bandmin(iband) .gt. efermi+sixsigma ) cycle bandloop2
            ! smearing parameter for this band
            isig=edr%default_smearing(iband,ispin)*sigma
            isig=1.0_r8/isig/sqrt(2.0_r8)
            do ikp=1,kp%n_irr_point
                if ( mod(ikp,mw%n) .ne. mw%r ) cycle
                f0=edr%ik(ikp)%eigenvalue(iband,ispin)
                f1=kp%ip(ikp)%integration_weight
                ! Gaussian smearing thing
                f2=(1-erf(isig*(f0-efermi) ))*0.5_r8*f1
                select case(edr%n_spin)
                case(1)
                    nelec=nelec+2*f2
                case(2)
                    nelec=nelec+f2
                end select
            enddo
        enddo bandloop2
        enddo spinloop
        ! Sync
        call mw%allreduce('sum',nelec)
        nelec=nelec+ctr
    end block gaussiansmearing
    case(2)
    ! Fermi-dirac smearing, nothing fancy
    fermidiracsmearing: block
        real(r8) :: thirtysixsigma,f0,f1,f2
        integer :: iband,ikp,ispin,ctr

        thirtysixsigma=36*temperature*lo_kb_Hartree
        ctr=0
        ! Count the easy ones
        bandloop1: do iband=1,edr%n_band
            ! check if it's a relevant band
            if ( edr%bandmax(iband) .lt. efermi-thirtysixsigma ) then
                ! this band is completely below
                select case(edr%n_spin)
                case(1)
                    ctr=ctr+2
                case(2)
                    ctr=ctr+1
                end select
            endif
        enddo bandloop1

        ! Count the less trivial ones
        nelec=0.0_r8
        bandloop2: do iband=1,edr%n_band
            ! check if it's a relevant band
            if ( edr%bandmax(iband) .lt. efermi-thirtysixsigma ) then
                cycle bandloop2
            elseif ( edr%bandmin(iband) .gt. efermi+thirtysixsigma ) then
                ! completely above, do nothing
                cycle bandloop2
            else
                ! somewhere in between
                do ispin=1,edr%n_spin
                do ikp=1,kp%n_irr_point
                    if ( mod(ikp,mw%n) .ne. mw%r ) cycle
                    f0=edr%ik(ikp)%eigenvalue(iband,ispin)
                    f1=kp%ip(ikp)%integration_weight
                    ! Fermi smearing thing
                    f2=lo_fermi(f0,efermi,temperature)*f1
                    select case(edr%n_spin)
                    case(1)
                        nelec=nelec+2*f2
                    case(2)
                        nelec=nelec+f2
                    end select
                enddo
                enddo
            endif
        enddo bandloop2
        ! Sync
        call mw%allreduce('sum',nelec)
        nelec=nelec+ctr
    end block fermidiracsmearing
    case(3)
    ! Tetrahedron integrate instead
    tetrahedronsmearing: block
        real(r8), parameter :: energytol=1E-7_r8*lo_eV_to_Hartree
        real(r8), parameter :: fixsigma=1E-5_r8*lo_eV_to_Hartree ! not used?
        real(r8), dimension(4) :: tete
        real(r8) :: weight
        integer :: iband,itet,ispin
        integer :: ctr


        ! Count the easy ones
        ctr=0
        bandloop1: do iband=1,edr%n_band
            ! check if it's a relevant band
            if ( edr%bandmax(iband) .lt. efermi-energytol ) then
                ! this band is completely below
                select case(edr%n_spin)
                case(1)
                    ctr=ctr+2
                case(2)
                    ctr=ctr+1
                end select
            endif
        enddo bandloop1

        ! Then the tricky ones
        nelec=0.0_r8
        bandloop2: do iband=1,edr%n_band
            ! check if it's a relevant band
            if ( edr%bandmax(iband) .lt. efermi-energytol ) cycle bandloop2
            if ( edr%bandmin(iband) .gt. efermi+energytol ) cycle bandloop2
            ! somewhere in between
            do ispin=1,edr%n_spin
            do itet=1,kp%n_irr_tet
                if ( mod(itet,mw%n) .ne. mw%r ) cycle
                tete=tetenergies(:,itet,ispin,iband)
                weight=lo_LV_tetrahedron_heaviside(tete,efermi,energytol)
                weight=weight*kp%it(itet)%integration_weight
                select case(edr%n_spin)
                case(1)
                    nelec=nelec+2*weight
                case(2)
                    nelec=nelec+weight
                end select
            enddo
            enddo
        enddo bandloop2
        ! Sync
        call mw%allreduce('sum',nelec)
        nelec=nelec+ctr
    end block tetrahedronsmearing
    case(4)
    ! Tetrahedron integrate with fancy smearing
    tetrahedronsmearingfermi: block
        real(r8), parameter :: energytol=1E-13_r8*lo_eV_to_Hartree
        real(r8), dimension(4) :: tete
        real(r8) :: weight,thirtysixsigma
        integer :: iband,itet,ispin
        integer :: ctr

        thirtysixsigma=36*temperature*lo_kb_Hartree

        ! Count the easy ones
        ctr=0
        bandloop1: do iband=1,edr%n_band
            ! check if it's a relevant band
            if ( edr%bandmax(iband) .lt. efermi-thirtysixsigma ) then
                ! this band is completely below
                select case(edr%n_spin)
                case(1)
                    ctr=ctr+2
                case(2)
                    ctr=ctr+1
                end select
            endif
        enddo bandloop1

        ! Then the tricky ones
        nelec=0.0_r8
        bandloop2: do iband=1,edr%n_band
            ! check if it's a relevant band
            if ( edr%bandmax(iband) .lt. efermi-thirtysixsigma ) cycle bandloop2
            if ( edr%bandmin(iband) .gt. efermi+thirtysixsigma ) cycle bandloop2
            ! somewhere in between
            do ispin=1,edr%n_spin
            do itet=1,kp%n_irr_tet
                if ( mod(itet,mw%n) .ne. mw%r ) cycle
                tete=tetenergies(:,itet,ispin,iband)
                weight=lo_LV_tetrahedron_fermi(tete,efermi,temperature,energytol)
                weight=weight*kp%it(itet)%integration_weight
                select case(edr%n_spin)
                case(1)
                    nelec=nelec+2*weight
                case(2)
                    nelec=nelec+weight
                end select
            enddo
            enddo
        enddo bandloop2
        ! Sync
        call mw%allreduce('sum',nelec)
        nelec=nelec+ctr
    end block tetrahedronsmearingfermi
    case(5)
    ! Linear sphere model
    linearsphere: block
        real(r8) :: sixsigma,isig,f0,f1,f2
        integer :: iband,ikp,ispin,ctr

        sixsigma=6*maxval(edr%default_smearing)
        ctr=0
        ! Count the easy ones
        bandloop1: do iband=1,edr%n_band
            ! check if it's a relevant band
            if ( edr%bandmax(iband) .lt. efermi-sixsigma ) then
                ! this band is completely below
                select case(edr%n_spin)
                case(1)
                    ctr=ctr+2
                case(2)
                    ctr=ctr+1
                end select
            endif
        enddo bandloop1

        nelec=0.0_r8
        spinloop: do ispin=1,edr%n_spin
        bandloop2: do iband=1,edr%n_band
            if ( edr%bandmax(iband) .lt. efermi-sixsigma ) cycle bandloop2
            if ( edr%bandmin(iband) .gt. efermi+sixsigma ) cycle bandloop2
            ! smearing parameter for this band
            isig=edr%default_smearing(iband,ispin)*sigma
            isig=1.0_r8/isig/sqrt(2.0_r8)
            do ikp=1,kp%n_irr_point
                if ( mod(ikp,mw%n) .ne. mw%r ) cycle
                f0=edr%ik(ikp)%eigenvalue(iband,ispin)
                f1=kp%ip(ikp)%integration_weight
                ! Gaussian smearing thing
                f2=(1-erf(isig*(f0-efermi) ))*0.5_r8*f1
                select case(edr%n_spin)
                case(1)
                    nelec=nelec+2*f2
                case(2)
                    nelec=nelec+f2
                end select
            enddo
        enddo bandloop2
        enddo spinloop
        ! Sync
        call mw%allreduce('sum',nelec)
        nelec=nelec+ctr
    end block linearsphere
    case default
        call lo_stop_gracefully(['Unknown smearing type'],lo_exitcode_param,__FILE__,__LINE__,mw%comm)
    end select
end function

!> measure size in memory, in bytes
function kpoint_size_in_mem(p) result(mem)
    !> one qpoint from the dispersions
    class(lo_electron_dispersions_kpoint), intent(in) :: p
    !> memory in bytes
    integer :: mem

    mem=0
    mem=mem+storage_size(p)

    if ( allocated( p%eigenvalue               ) ) mem=mem+storage_size(p%eigenvalue              )*size(p%eigenvalue              )
    if ( allocated( p%eigenvector              ) ) mem=mem+storage_size(p%eigenvector             )*size(p%eigenvector             )
    if ( allocated( p%groupvelocity            ) ) mem=mem+storage_size(p%groupvelocity           )*size(p%groupvelocity           )
    if ( allocated( p%electronphononlinewidth  ) ) mem=mem+storage_size(p%electronphononlinewidth )*size(p%electronphononlinewidth )
    if ( allocated( p%electronphononphasespace ) ) mem=mem+storage_size(p%electronphononphasespace)*size(p%electronphononphasespace)
    mem=mem/8
end function

!> measure size in memory, in bytes
function edr_size_in_mem(edr) result(mem)
    !> dispersions
    class(lo_electron_dispersions), intent(in) :: edr
    !> memory in bytes
    integer :: mem

    integer :: i
    mem=0
    ! easy things
    mem=mem+storage_size(edr)
    if ( allocated(edr%bandmin) ) mem=mem+storage_size(edr%bandmin)*size(edr%bandmin)
    if ( allocated(edr%bandmax) ) mem=mem+storage_size(edr%bandmax)*size(edr%bandmax)
    mem=mem/8
    ! more annoying things
    if ( allocated(edr%ik) ) then
        do i=1,size(edr%ik)
             mem=mem+edr%ik(i)%size_in_mem()
        enddo
    endif
end function

!> destroy!
subroutine destroy(edr)
    !> dispersions
    class(lo_electron_dispersions), intent(inout) :: edr

    if ( allocated(edr%bandmin) ) deallocate(edr%bandmin)
    if ( allocated(edr%bandmax) ) deallocate(edr%bandmax)
    if ( allocated(edr%ik     ) ) deallocate(edr%ik     )
    edr%n_band=-lo_hugeint
    edr%n_basis=lo_hugeint
    edr%n_spin=-lo_hugeint
    edr%n_electron=-lo_huge
    edr%efermi=-lo_huge
    edr%bandgap=-lo_huge
    edr%evalence=-lo_huge
    edr%n_irr_kpoint=-lo_hugeint
    edr%n_full_kpoint=-lo_hugeint
    edr%energy_max=-lo_huge
    edr%energy_min=-lo_huge
end subroutine

! !> Read the electron eigenvalues from a VASP EIGENVAL file
! subroutine readfromfile(edr,kp,filename,p,verbosity)
!     !> the electronic dispersion relations
!     class(lo_electron_dispersions), intent(out) :: edr
!     !> the k-point mesh
!     class(lo_qpoint_mesh), intent(in) :: kp
!     !> the filename
!     character(len=*), intent(in) :: filename
!     !> the crystal structure
!     type(lo_crystalstructure), intent(in) :: p
!     !> number of spin channels
!     !> talk?
!     integer, intent(in), optional :: verbosity
!
!     ! grab everything from file
!     initandread: block
!         character(len=1) :: trams
!         real(r8), dimension(3) :: v0
!         real(r8) :: f0
!         integer, dimension(:), allocatable :: dum
!         integer :: i,j,k,l,u
!
!         ! set some standard stuff
!         edr%n_irr_kpoint=kp%nq_irr
!         edr%n_full_kpoint=kp%nq_tot
!         if ( present(verbosity) ) then
!             verbosity=verbosity
!         else
!             verbosity=0
!         endif
!
!         if ( verbosity .gt. 0 ) then
!             write(lo_iou,*) ''
!             write(lo_iou,*) 'Reading electron bandstructure from vasp eigenval file'
!         endif
!         ! start reading the file
!         u=open_file('in',trim(filename))
!             ! skip the header
!             read(u,*) i,j,k,l
!             ! I hope the third number is number of spin channels?
!             edr%n_spin=k
!             read(u,*) trams
!             read(u,*) trams
!             read(u,*) trams
!             read(u,*) trams
!             read(u,*) i,j,k
!             ! sanity check:
!             if ( j .ne. edr%n_irr_kpoint ) then
!                 write(lo_iou,*) 'Found',j,'irreducible points, expected',edr%n_irr_kpoint
!                 stop
!             endif
!             ! store number of bands and electrons
!             edr%n_band=k
!             edr%n_electron=i
!             if ( verbosity .gt. 0 ) write(lo_iou,*) '... parsed header, # bands: ',tochar(edr%n_band),'  # electrons: ',tochar(edr%n_electron)
!             ! Allocate stuff
!             lo_allocate(edr%ik(edr%n_irr_kpoint))
!             do i=1,edr%n_irr_kpoint
!                 lo_allocate(edr%ik(i)%eigenvalue(edr%n_band,edr%n_spin))
!                 edr%ik(i)%eigenvalue=0.0_r8
!             enddo
!             lo_allocate(dum(edr%n_irr_kpoint))
!             dum=0
!             ! Read stuff
!             do i=1,edr%n_irr_kpoint
!                 ! Got the vector in fractional coordinates
!                 read(u,*) v0,f0
!                 ! Sanity check, vasp screws around with the points
!                 v0=v0-matmul(p%inv_reciprocal_latticevectors,kp%ip(i)%w)
!                 v0=p%displacement_fractional_to_cartesian(v0)
!                 if ( lo_sqnorm(v0) .gt. lo_sqtol*100 ) then
!                     write(lo_iou,*) 'incompatible qgrid and eigenval'
!                     stop
!                 endif
!                 ! read the energies
!                 do j=1,edr%n_band
!                     select case(edr%n_spin)
!                     case(1)
!                         read(u,*) k,edr%ik(i)%eigenvalue(j,1)
!                     case(2)
!                         read(u,*) k,edr%ik(i)%eigenvalue(j,1),f0,edr%ik(i)%eigenvalue(j,2)
!                     end select
!                 enddo
!                 ! sort the energies by size, just to be on the safe side
!                 select case(edr%n_spin)
!                 case(1)
!                     call qsort(edr%ik(i)%eigenvalue(:,1))
!                 case(2)
!                     call qsort(edr%ik(i)%eigenvalue(:,1))
!                     call qsort(edr%ik(i)%eigenvalue(:,2))
!                 end select
!             enddo
!         close(u)
!         if ( verbosity .gt. 0 ) write(lo_iou,*) '... read energies for ',tochar(edr%n_irr_kpoint),' irreducible k-points'
!
!         ! fold out the energies to all k-points, and figure out min and max energies
!         edr%energy_min=lo_huge
!         edr%energy_max=-lo_huge
!         lo_allocate(edr%ak(edr%n_full_kpoint))
!         do i=1,edr%n_full_kpoint
!             lo_allocate(edr%ak(i)%eigenvalue(edr%n_band,edr%n_spin))
!             j=kp%ap(i)%irrind
!             edr%ak(i)%eigenvalue=edr%ik(j)%eigenvalue
!             edr%energy_min=min(edr%energy_min,minval(edr%ak(i)%eigenvalue))
!             edr%energy_max=max(edr%energy_max,maxval(edr%ak(i)%eigenvalue))
!         enddo
!         if ( verbosity .gt. 0 ) write(lo_iou,*) '... folded out to ',tochar(edr%n_full_kpoint),' k-points'
!
!         ! It's neat with a bandwith of sorts.
!         lo_allocate(edr%bandmin(edr%n_band))
!         lo_allocate(edr%bandmax(edr%n_band))
!         edr%bandmin=lo_huge
!         edr%bandmax=-lo_huge
!         do i=1,edr%n_irr_kpoint
!         do k=1,edr%n_spin
!         do j=1,edr%n_band
!             edr%bandmin(j)=min(edr%bandmin(j),edr%ik(i)%eigenvalue(j,k))
!             edr%bandmax(j)=max(edr%bandmax(j),edr%ik(i)%eigenvalue(j,k))
!         enddo
!         enddo
!         enddo
!         if ( verbosity .gt. 0 ) write(lo_iou,*) '... determined bandwidths'
!     end block initandread
! end subroutine

end module
