module lo_distributed_phonon_dispersion_relations
! Same as phonon dispersions, only distributed over ranks
! Should be more conservative with respect to memory
use konstanter, only: r8,lo_iou,lo_huge,lo_hugeint,lo_exitcode_param,lo_sqtol,lo_degenvector,lo_freqtol
use gottochblandat, only: lo_sqnorm,lo_progressbar,lo_progressbar_init,walltime,lo_trueNtimes,tochar
use type_qpointmesh, only: lo_qpoint_mesh,lo_qpoint,lo_get_small_group_of_qpoint
use type_crystalstructure, only: lo_crystalstructure
use type_forceconstant_secondorder, only: lo_forceconstant_secondorder
use type_blas_lapack_wrappers, only: lo_gemm,lo_zheev
use type_symmetryoperation, only: lo_eigenvector_transformation_matrix, lo_operate_on_vector
use mpi_wrappers, only: lo_mpi_helper, lo_stop_gracefully
use lo_memtracker, only: lo_mem_helper
use lo_sorting, only: lo_qsort

implicit none
private
public :: lo_generate_distributed_dispersions
public :: lo_distributed_phonon_dispersions
public :: lo_distributed_phonon_dispersions_qpoint

!> A q-point in distributed phonon dispersion
type lo_distributed_phonon_dispersions_qpoint
    !> frequencies
    real(r8), dimension(:), allocatable :: omega
    !> group velocities
    real(r8), dimension(:, :), allocatable :: vel
    !> generalized group velocities
    real(r8), dimension(:, :, :), allocatable :: genvel
    !> mode eigenvectors
    complex(r8), dimension(:, :), allocatable :: egv
    !> eigenvectors / sqrt( 2 omega_{qs} m_i)
    complex(r8), dimension(:, :), allocatable :: nuvec
    !> how degenerate is the mode?
    integer, dimension(:), allocatable :: degeneracy
    !> with what mode is it degenerate?
    integer, dimension(:, :), allocatable :: degenmode
    !> smearing parameter for adaptive integrations
    real(r8), dimension(:), allocatable :: sigma
    contains
        procedure :: generate=>harmonic_things_for_single_point
        procedure :: destroy=>destroy_ompoint
end type

!> irreducible point
type, extends(lo_distributed_phonon_dispersions_qpoint) :: lo_distributed_phonon_dispersions_qpoint_irr
    !> what is the global irreducible index, in the q-mesh
    integer :: global_irreducible_index=-lo_hugeint
    !> how many full points does this point transform to
    integer :: n_full_point=-lo_hugeint
    !> index (local) to full points
    integer, dimension(:), allocatable :: local_full_index
    !> operation to (local) full point
    integer, dimension(:), allocatable :: local_operation_full
end type

!> full point
type, extends(lo_distributed_phonon_dispersions_qpoint) :: lo_distributed_phonon_dispersions_qpoint_full
    !> what is the global index, in the q-mesh
    integer :: global_full_index=-lo_hugeint
    !> what is the local irreducible index
    integer :: local_irr_index=-lo_hugeint
    !> operation that takes the irreducible here
    integer :: operation_from_irreducible=-lo_hugeint
end type

type :: lo_distributed_phonon_dispersions
    !> how many phonon modes
    integer :: n_mode = -lo_hugeint
    !> how many irreducible q-points (on this rank)
    integer :: n_irr_qpoint_local=-lo_hugeint
    !> how many full q-points (on this rank)
    integer :: n_full_qpoint_local=-lo_hugeint
    !> largest frequency
    real(r8) :: omega_max = -lo_huge
    !> smallest nonzero frequency
    real(r8) :: omega_min = -lo_huge
    !> default smearing parameter, per mode
    real(r8), dimension(:), allocatable :: default_smearing
    !> irreducible points
    type(lo_distributed_phonon_dispersions_qpoint_irr), dimension(:), allocatable :: iq
    !> full points
    type(lo_distributed_phonon_dispersions_qpoint_full), dimension(:), allocatable :: aq
    contains
        procedure :: generate=>lo_generate_distributed_dispersions
end type

contains

subroutine lo_generate_distributed_dispersions(dr,qp,p,fc,smearing_prefactor,mw,mem,verbosity)
    !> distributed dispersions
    class(lo_distributed_phonon_dispersions), intent(out) :: dr
    !> q-point mesh
    class(lo_qpoint_mesh), intent(in) :: qp
    !> structure
    type(lo_crystalstructure), intent(in) :: p
    !> second order force constant
    type(lo_forceconstant_secondorder), intent(inout) :: fc
    !> adaptive smearing prefactor
    real(r8), intent(in) :: smearing_prefactor
    !> MPI communicator
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem
    !> how much to talk
    integer, intent(in) :: verbosity

    real(r8) :: t0,t1,timer
    integer, dimension(:,:), allocatable :: assignment_to_rank

    init: block
        t0=walltime()
        timer=t0
        if ( verbosity .gt. 0 ) then
            write(lo_iou,*) ''
            write(lo_iou,*) 'Generating distributed phonon dispersions'
        endif

        ! Set number of modes, at least I know that much
        dr%n_mode=p%na*3
    end block init

    distribute: block
        integer, dimension(:), allocatable :: points_per_rank

        ! This is to have things reasonably load-balanced
        call assign_rank_to_qpoints(qp,mw,assignment_to_rank,points_per_rank)
        dr%n_irr_qpoint_local = points_per_rank(mw%r+1)
        deallocate(points_per_rank)

        if ( dr%n_irr_qpoint_local .gt. 0 ) allocate(dr%iq(dr%n_irr_qpoint_local))
    end block distribute

    ! Calculate only the irreducible ones
    evaluate: block
        integer :: iq,jq

        if (verbosity .gt. 0) call lo_progressbar_init()

        do iq=1,dr%n_irr_qpoint_local
            jq=assignment_to_rank(iq,mw%r+1)
            call dr%iq(iq)%generate(fc,p,mem,qpoint=qp%ip(jq))

            if (verbosity .gt. 0) then
                if (lo_trueNtimes(iq, 127, dr%n_irr_qpoint_local)) call lo_progressbar(' ... irreducible dispersions', iq, dr%n_irr_qpoint_local, walltime() - t0)
            end if
        enddo

        if (verbosity .gt. 0) then
            t1 = walltime()
            call lo_progressbar(' ... irreducible dispersions', dr%n_irr_qpoint_local, dr%n_irr_qpoint_local, t1 - t0)
            t0 = t1
        end if
    end block evaluate

    ! Expand irreducible to full set
    expand: block
        complex(r8), dimension(:,:), allocatable :: rotmat
        integer :: i,j,k,l,iq,jq,ctr,iop

        ! First make a little space for the full
        j=0
        do i=1,dr%n_irr_qpoint_local
            iq=assignment_to_rank(i,mw%r+1)
            j=j+qp%ip(iq)%n_full_point
        enddo
        dr%n_full_qpoint_local=j
        if ( j .gt. 0 ) allocate(dr%aq(dr%n_full_qpoint_local))

        ! Start storing some information
        ctr=0
        do i=1,dr%n_irr_qpoint_local
            ! Store which q-point this is
            iq=assignment_to_rank(i,mw%r+1)
            dr%iq(i)%global_irreducible_index=iq
            dr%iq(i)%n_full_point = qp%ip(iq)%n_full_point
            allocate(dr%iq(i)%local_full_index( dr%iq(i)%n_full_point ))
            allocate(dr%iq(i)%local_operation_full( dr%iq(i)%n_full_point ))
            dr%iq(i)%local_full_index=0
            dr%iq(i)%local_operation_full=0
            do j=1,dr%iq(i)%n_full_point
                ctr=ctr+1
                jq=qp%ip(iq)%index_full_point(j)
                dr%iq(i)%local_operation_full(j) = qp%ip(iq)%operation_full_point( j )
                dr%iq(i)%local_full_index(j) = ctr
                dr%aq(ctr)%global_full_index = jq
                dr%aq(ctr)%local_irr_index = i
                dr%aq(ctr)%operation_from_irreducible = qp%ap(jq)%operation_from_irreducible
            enddo
        enddo

        ! Now I can start transforming
        allocate(rotmat(p%na*3,p%na*3))
        rotmat=0.0_r8

        do i=1,dr%n_full_qpoint_local
            allocate(dr%aq(i)%omega( p%na*3 ))
            allocate(dr%aq(i)%egv( p%na*3,p%na*3 ))
            allocate(dr%aq(i)%vel( 3,p%na*3 ))
            allocate(dr%aq(i)%genvel( 3,p%na*3,p%na*3 ))
            allocate(dr%aq(i)%nuvec( p%na*3,p%na*3 ))
            allocate(dr%aq(i)%sigma( p%na*3 ))
            allocate(dr%aq(i)%degeneracy( p%na*3 ))
            j=dr%aq(i)%local_irr_index
            k=size(dr%iq(j)%degenmode,1)
            allocate(dr%aq(i)%degenmode(k, p%na*3 ))

            dr%aq(i)%omega=0.0_r8
            dr%aq(i)%egv=0.0_r8
            dr%aq(i)%vel=0.0_r8
            dr%aq(i)%genvel=0.0_r8
            dr%aq(i)%nuvec=0.0_r8
            dr%aq(i)%degeneracy=0
            dr%aq(i)%degenmode=0

            ! Start filling in
            iop=dr%aq(i)%operation_from_irreducible
            iq=dr%iq(j)%global_irreducible_index
            jq=dr%aq(i)%global_full_index
            if (iop .gt. 0) then
                ! the eigenvector transformation matrix
                call lo_eigenvector_transformation_matrix(rotmat, p%rcart, qp%ip(iq)%r, p%sym%op(abs(iop)))
                ! rotate the eigenvectors
                call lo_gemm(rotmat, dr%iq(j)%egv, dr%aq(i)%egv)
                do k = 1, dr%n_mode
                    dr%aq(i)%omega(k) = dr%iq(j)%omega(k)
                    ! Make sure the velocities have the proper symmetry
                    dr%aq(i)%vel(:,k) = lo_operate_on_vector(p%sym%op(abs(iop)), dr%iq(j)%vel(:, k), reciprocal=.true.)
                    do l=1,dr%n_mode
                        dr%aq(i)%genvel(:,k,l) = lo_operate_on_vector(p%sym%op(abs(iop)), dr%iq(j)%genvel(:, k, l), reciprocal=.true.)
                    enddo
                end do
            else
                ! Same thing, but now we have to time-reverse it
                call lo_eigenvector_transformation_matrix(rotmat, p%rcart, qp%ip(iq)%r, p%sym%op(abs(iop)))
                ! rotate the eigenvectors
                call lo_gemm(rotmat, dr%iq(j)%egv, dr%aq(i)%egv)
                ! and conjugate them, that's what Maradudin says.
                dr%aq(i)%egv = conjg(dr%aq(i)%egv)
                do k = 1, dr%n_mode
                    dr%aq(i)%omega(k) = dr%iq(j)%omega(k)
                    ! Make sure the velocities have the proper symmetry. Note negative sign here.
                    dr%aq(i)%vel(:,k) = -lo_operate_on_vector(p%sym%op(abs(iop)), dr%iq(j)%vel(:,k), reciprocal=.true.)
                end do
            end if
        enddo

        deallocate(rotmat)
    end block expand

    ! Get the default smearing parameters, per band
    baselinesmearing: block
        real(r8), dimension(:), allocatable :: x
        real(r8) :: f0
        integer :: imode, iqp,i

        ! Find min and max frequency, always useful.
        dr%omega_max = -lo_huge
        dr%omega_min = lo_huge
        do iqp = 1, dr%n_irr_qpoint_local
        do imode = 1, dr%n_mode
            f0 = dr%iq(iqp)%omega(imode)
            if (abs(f0) .gt. lo_freqtol) dr%omega_min = min(dr%omega_min, f0)
            dr%omega_max = max(dr%omega_max, f0)
        end do
        end do
        call mw%allreduce('max',dr%omega_max)
        call mw%allreduce('min',dr%omega_min)

        ! Get the default smearing parameter for Gaussian integrations, per band.
        allocate (dr%default_smearing(dr%n_mode))
        dr%default_smearing = 0.0_r8

        allocate(x(qp%n_irr_point))
        x = 0.0_r8

        do imode=1,dr%n_mode
            x = 0.0_r8
            do iqp = 1, dr%n_irr_qpoint_local
                i=dr%iq(iqp)%global_irreducible_index
                x(i)=dr%iq(iqp)%omega(imode)
            enddo
            call mw%allreduce('max',x)
            call lo_qsort(x)
            f0=0.0_r8
            do iqp = 1, qp%n_irr_point - 1
                f0 = max(f0, abs(x(iqp + 1) - x(iqp)))
            end do
            dr%default_smearing(imode) = f0
        enddo

        ! Now, if I havey very flat bands this can become slightly problematic
        ! and the DOS gets extremely spiky. So I set some sensible minimum
        ! value for the smearing.
        f0 = maxval(dr%default_smearing)/5.0_r8
        dr%default_smearing = max(dr%default_smearing, f0)
    end block baselinesmearing

    ! Pre-calculate some useful things
    precalc: block
        real(r8) :: f0,f1
        integer :: i,j,k,iq,imode,iatom,ialpha,ii

        ! Get the adaptive smearing parameter for integrations
        do i=1,dr%n_irr_qpoint_local
            iq=dr%iq(i)%global_irreducible_index
            do imode=1,dr%n_mode
                if ( dr%iq(i)%omega(imode) .gt. lo_freqtol ) then
                    dr%iq(i)%sigma(imode) = qp%adaptive_sigma( qp%ip(iq)%radius,dr%iq(i)%vel(:,imode),dr%default_smearing(imode),smearing_prefactor)
                else
                    dr%iq(i)%sigma(imode) = -1.0_r8
                endif
            enddo
            ! and fill this out to the rest
            do j=1,dr%iq(i)%n_full_point
                k=dr%iq(i)%local_full_index(j)
                dr%aq(k)%sigma = dr%iq(i)%sigma
            enddo
        enddo

        ! Pre-stuff all prefactors into the eigenvectors, for the irreducible it is already done.
        do i=1,dr%n_full_qpoint_local
            do imode=1,dr%n_mode
                if ( dr%aq(i)%omega(imode) .gt. lo_freqtol ) then
                    f0=1.0_r8/sqrt(2.0_r8*dr%aq(i)%omega(imode))
                else
                    f0=0.0_r8
                endif
                do iatom=1,p%na
                    f1=p%invsqrtmass(iatom)
                    do ialpha=1,3
                        ii=(iatom-1)*3+ialpha
                        dr%aq(i)%nuvec(ii,imode)=dr%aq(i)%egv(ii,imode)*f0*f1
                    enddo
                enddo
            enddo
        enddo
    end block precalc

    deallocate(assignment_to_rank)

    if (verbosity .gt. 0) then
        write (lo_iou, *) 'Calculated the full dispersion relations in ', tochar(walltime() - timer), 's'
    end if
end subroutine

!> first-shot something at distributing q-points across ranks somewhat evenly.
subroutine assign_rank_to_qpoints(qp,mw,assignment_to_rank,points_per_rank)
    !> q-point mesh
    class(lo_qpoint_mesh), intent(in) :: qp
    !> MPI communicator
    type(lo_mpi_helper), intent(inout) :: mw
    !> Assignment to rank
    integer, dimension(:,:), allocatable, intent(out) :: assignment_to_rank
    !> Number of points per rank
    integer, dimension(:), allocatable, intent(out) :: points_per_rank

    ! When distributing the q-point I want to distribute the
    ! irreducible set, and then unfold the irreducible. But
    ! since each irreducible corresponds to a different number
    ! of full q-points it might get very unbalanced in the end.
    ! so I'm going to need some kind of heuristics how to go
    ! about this.

    integer, parameter :: maxiter=500
    integer, dimension(:,:), allocatable :: assignment
    integer, dimension(:), allocatable :: weight
    integer, dimension(:), allocatable :: di,dj
    integer :: i,j,k,iq,ctr,iter,imin,imax,ii,jj,kk,mismatch

    allocate(points_per_rank(mw%n))
    allocate(weight(mw%n))
    points_per_rank=0
    weight=0

    ! Count points per rank
    points_per_rank=0
    do i=1,qp%n_irr_point
        j=mod(i,mw%n)+1
        points_per_rank(j)=points_per_rank(j)+1
    enddo
    j=maxval(points_per_rank)
    allocate(assignment(j,mw%n))
    allocate(assignment_to_rank(j,mw%n))
    allocate(di(j))
    allocate(dj(j))
    assignment=-lo_hugeint
    assignment_to_rank=-lo_hugeint
    di=-lo_hugeint
    dj=-lo_hugeint

    ! Start with a stupid assignment
    points_per_rank=0
    do i=1,qp%n_irr_point
        j=mod(i,mw%n)+1
        points_per_rank(j)=points_per_rank(j)+1
        assignment( points_per_rank(j),j )=i
    enddo

    ! Get the weight per rank, also sort assignment so
    ! that the most expensive q-point is first
    weight=0
    do i=1,mw%n
        di=1
        do j=1,points_per_rank(i)
            iq=assignment(j,i)
            di(j)=-qp%ip(iq)%n_full_point
            weight(i)=weight(i) + qp%ip(iq)%n_full_point
        enddo
        call lo_qsort(di(1:points_per_rank(i)), order=dj(1:points_per_rank(i)))
        di(1:points_per_rank(i))=assignment(1:points_per_rank(i),i)
        do j=1,points_per_rank(i)
            assignment(j,i)=di(dj(j))
        enddo
    enddo

    mismatch=maxval(weight)-minval(weight)
    assignment_to_rank=assignment

    ! First we do a simple worst-to-best flip for a little while
    ctr=0
    iterloop: do iter=1,maxiter

        ! Find the two extreme ranks
        imin=-1
        imax=-1
        ii=lo_hugeint
        jj=-lo_hugeint
        do i=1,mw%n
            if ( weight(i) .eq. 0 ) cycle
            if ( weight(i) .le. ii ) then
                ii=weight(i)
                imin=i
            endif
            if ( weight(i) .ge. jj ) then
                jj=weight(i)
                imax=i
            endif
        enddo

        if ( weight(imin) .eq. weight(imax) ) then
            ! we are perfectly balanced, highly unlikely but whatever
            exit
        endif

        ! Flip the largest to the smallest?
        ii=weight(imax)-weight(imin)
        kk=min(points_per_rank(imax),points_per_rank(imin))
        do k=1,kk
            i=assignment(k,imax)
            j=assignment(points_per_rank(imin)+1-k,imin)

            assignment(k,imax)=j
            assignment(points_per_rank(imin)+1-k,imin)=i
            ! Re-calculate weight
            jj=0
            do i=1,points_per_rank(imax)
                iq=assignment(i,imax)
                jj=jj+qp%ip(iq)%n_full_point
            enddo
            do i=1,points_per_rank(imin)
                iq=assignment(i,imin)
                jj=jj-qp%ip(iq)%n_full_point
            enddo
            if ( jj .le. 0 ) exit
        enddo

        ! Re-sort
        i=imin
        do j=1,points_per_rank(i)
            iq=assignment(j,i)
            di(j)=-qp%ip(iq)%n_full_point
        enddo
        call lo_qsort(di(1:points_per_rank(i)), order=dj(1:points_per_rank(i)))
        di(1:points_per_rank(i))=assignment(1:points_per_rank(i),i)
        do j=1,points_per_rank(i)
            assignment(j,i)=di(dj(j))
        enddo
        i=imax

        do j=1,points_per_rank(i)
            iq=assignment(j,i)
            di(j)=-qp%ip(iq)%n_full_point
        enddo
        call lo_qsort(di(1:points_per_rank(i)), order=dj(1:points_per_rank(i)))
        di(1:points_per_rank(i))=assignment(1:points_per_rank(i),i)
        do j=1,points_per_rank(i)
            assignment(j,i)=di(dj(j))
        enddo
        ! Re-calculate weights
        weight=0
        do i=1,mw%n
            do j=1,points_per_rank(i)
                iq=assignment(j,i)
                weight(i)=weight(i) + qp%ip(iq)%n_full_point
            enddo
        enddo

        i=maxval(weight)-minval(weight)
        if ( i .lt. mismatch ) then
            ! This was a good flip
            assignment_to_rank=assignment
            mismatch=i
        else
            ! bad flip
            ctr=ctr+1
            if ( ctr .ge. 5 ) exit iterloop
        endif
    enddo iterloop

    if ( allocated(assignment) ) deallocate(assignment)
    if ( allocated(weight    ) ) deallocate(weight    )
    if ( allocated(di        ) ) deallocate(di        )
    if ( allocated(dj        ) ) deallocate(dj        )

    ! If needed I can further refine here, but so far this seems
    ! like an ok scheme. Like a monte-carlo flip thingy could be added.
end subroutine

!> Calculate all the harmonic things for a q-point
subroutine harmonic_things_for_single_point(ompoint, fc, p, mem, qpoint, qvec, qdirection)
    !> point in the dispersions
    class(lo_distributed_phonon_dispersions_qpoint), intent(inout) :: ompoint
    !> forceconstant
    type(lo_forceconstant_secondorder), intent(inout) :: fc
    !> structure
    type(lo_crystalstructure), intent(in) :: p
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem
    !> q-point
    class(lo_qpoint), intent(in), optional :: qpoint
    !> q-point as just a vector
    real(r8), dimension(3), intent(in), optional :: qvec
    !> q-direction?
    real(r8), dimension(3), intent(in), optional :: qdirection

    real(r8), dimension(3), parameter :: default_direction = [1.0_r8, 0.0_r8, 0.0_r8]
    type(lo_qpoint) :: dumqpoint
    real(r8), dimension(3) :: qdir
    integer :: nb
    logical :: skipna

    ! So this is slightly too complicated for something this simple, one would think.
    ! But it's really not that bad. Just makes sense to make the thing flexible.
    ! Sort out some simple things first:
    init: block
        ! Really basic sanity tests. Will make more when I think of it.
        if (p%na .ne. fc%na) then
            call lo_stop_gracefully(['Different number of atoms in structure and forceconstant.'],lo_exitcode_param, __FILE__, __LINE__)
        end if
        ! Make sure we start fresh
        call ompoint%destroy()
        ! Get the number of modes
        nb = fc%na*3

        allocate (ompoint%omega(nb))
        allocate (ompoint%egv(nb, nb))
        allocate (ompoint%nuvec(nb, nb))
        allocate (ompoint%vel(3, nb))
        allocate (ompoint%genvel(3, nb, nb))
        allocate (ompoint%sigma(nb))
        ompoint%omega = 0.0_r8
        ompoint%egv = 0.0_r8
        ompoint%nuvec = 0.0_r8
        ompoint%vel = 0.0_r8
        ompoint%genvel = 0.0_r8
        ompoint%sigma = 0.0_r8

        ! Figure out what to do with the q-point
        if (present(qpoint) .and. present(qvec)) then
            call lo_stop_gracefully(['Either provide a q-point, or a q-vector, not both'], &
                                    lo_exitcode_param, __FILE__, __LINE__)
        end if

        ! Is there a direction specified?
        if (present(qdirection)) then
            if (lo_sqnorm(qdirection) .gt. lo_sqtol) then
                qdir = qdirection
                skipna = .false.
            else
                qdir = 0.0_r8
                skipna = .true.
            end if
        else
            qdir = default_direction
            skipna = .false.
        end if
        qdir = qdir/norm2(qdir)

        ! If there is no q-point supplied, construct one from the input vector
        if (present(qvec)) then
            dumqpoint%r = qvec - p%bz%gshift(qvec + lo_degenvector)
            call lo_get_small_group_of_qpoint(dumqpoint, p)
        end if
    end block init

    ! Solve for harmonic things
    slv: block
        complex(r8), dimension(:, :, :), allocatable :: Dq
        complex(r8), dimension(:, :), allocatable :: D,m0,m1,m2,m3
        real(r8), dimension(:), allocatable :: v0,v1
        real(r8) :: f0,f1
        integer, dimension(:, :), allocatable :: di
        logical, dimension(:), allocatable :: modefixed
        integer :: l,b1,b2,i,j,iatom,ialpha,ii,jj

        ! Get the frequencies
        call mem%allocate(Dq, [nb, nb, 3], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(D,  [nb, nb], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(m0, [nb, nb], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(di, [nb, nb], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        Dq = 0.0_r8
        D = 0.0_r8
        m0 = 0.0_r8
        di = 0
        if (present(qpoint)) then
            call fc%dynamicalmatrix( &
                p, qpoint, D, mem, Dq, qdirection=qdir, skipnonanalytical=skipna)
            call fc%frequencies_eigenvectors_groupvelocities( &
                D, ompoint%omega, mem, Dq, eigenvectors=ompoint%egv, groupvelocities=ompoint%vel, qpoint=qpoint)
        else
            call fc%dynamicalmatrix( &
                p, dumqpoint, D, mem, Dq, qdirection=qdir, skipnonanalytical=skipna)
            call fc%frequencies_eigenvectors_groupvelocities( &
                D, ompoint%omega, mem, Dq, eigenvectors=ompoint%egv, groupvelocities=ompoint%vel, qpoint=dumqpoint)
        end if

        ! Sort out degeneracies, first count largest amount of degeneracies
        i=0
        do b1 = 1, nb
            j=0
            do b2 = b1, nb
                if (abs(ompoint%omega(b1) - ompoint%omega(b2)) .lt. lo_freqtol) j=j+1
            enddo
            i=max(i,j)
        enddo

        allocate(ompoint%degeneracy(nb))
        allocate(ompoint%degenmode(i, nb))
        ompoint%degeneracy=0
        ompoint%degenmode=0

        ! Then we store the actual degeneracies
        di = 0
        do b1 = 1, nb
        do b2 = b1, nb
            if (abs(ompoint%omega(b1) - ompoint%omega(b2)) .lt. lo_freqtol) then
                di(b1, b2) = 1
                di(b2, b1) = 1
            end if
        end do
        end do
        do b1 = 1, nb
            ompoint%degeneracy(b1) = sum(di(:, b1))
        end do
        ompoint%degenmode = 0
        do b1 = 1, nb
            l = 0
            do b2 = 1, nb
                if (di(b1, b2) .eq. 1) then
                    l = l + 1
                    ompoint%degenmode(l, b1) = b2
                end if
            end do
        end do

        ! Pre-calculate the scaled eigenvectors
        do b1=1, nb
            if ( ompoint%omega(b1) .gt. lo_freqtol ) then
                f0=1.0_r8/sqrt(2.0_r8*ompoint%omega(b1))
            else
                f0=0.0_r8
            endif
            do iatom=1,p%na
                f1=p%invsqrtmass(iatom)
                do ialpha=1,3
                    ii=(iatom-1)*3+ialpha
                    ompoint%nuvec(ii,b1)=ompoint%egv(ii,b1)*f0*f1
                enddo
            enddo
        enddo

        ! Evaluate off-diagonal group velocities
        allocate(modefixed(nb))
        allocate(v1(nb))
        v1=0.0_r8
        modefixed=.false.
        do ialpha=1,3
            ! Find the properly rotated eigenvectors
            modefixed=.false.
            v1=0.0_r8
            m0=0.0_r8
            do b1=1,nb
                if ( modefixed(b1) ) cycle

                if ( ompoint%omega(b1) .lt. lo_freqtol ) then
                    ! Don't care for acoustic modes
                    m0(:,b1)=0.0_r8
                    modefixed(b1)=.true.
                elseif ( ompoint%degeneracy(b1) .eq. 1 ) then
                    ! Not degenerate, nothing to worry about
                    m0(:,b1)=ompoint%egv(:,b1)
                    modefixed(b1)=.true.
                else
                    allocate(m1(ompoint%degeneracy(b1),ompoint%degeneracy(b1)))
                    allocate(m2(nb,ompoint%degeneracy(b1)))
                    allocate(m3(nb,ompoint%degeneracy(b1)))
                    allocate(v0(ompoint%degeneracy(b1)))
                    m1=0.0_r8
                    m2=0.0_r8
                    m3=0.0_r8
                    v0=0.0_r8

                    do i=1,ompoint%degeneracy(b1)
                    do j=1,ompoint%degeneracy(b1)
                        ii=ompoint%degenmode(i,b1)
                        jj=ompoint%degenmode(j,b1)
                        m1(i,j)=dot_product(ompoint%egv(:,ii),matmul(Dq(:,:,ialpha),ompoint%egv(:,jj)))
                    enddo
                    enddo
                    ! Figure out the mean eigenvalue
                    call lo_zheev(m1,v0,jobz='V')
                    f0=sum(v0)/real(ompoint%degeneracy(b1),r8)
                    ! Rotate eigenvectors
                    do i=1,ompoint%degeneracy(b1)
                        ii=ompoint%degenmode(i,b1)
                        m2(:,i)=ompoint%egv(:,ii)
                    enddo
                    m3=matmul(m2,m1)
                    ! Store correctly rotated eigenvectors
                    do i=1,ompoint%degeneracy(b1)
                        b2=ompoint%degenmode(i,b1)
                        m0(:,b2)=m3(:,i)
                        v1(b2)=f0
                        modefixed(b2)=.true.
                    enddo

                    deallocate(m1)
                    deallocate(m2)
                    deallocate(m3)
                    deallocate(v0)
                endif
            enddo

            ! Sandwich with rotated eigenvectors
            allocate(m2(nb,nb))
            m2=0.0_r8
            do b1=1,nb
            do b2=1,nb
                m2(b1,b2)=dot_product(m0(:,b1),matmul(Dq(:,:,ialpha),m0(:,b2)))
            enddo
            enddo
            ! Make sure the degeneracies hold
            do b1=1,nb
                if ( ompoint%degeneracy(b1) .eq. 1 ) cycle
                do i=1,ompoint%degeneracy(b1)
                    b2=ompoint%degenmode(i,b1)
                    if ( b1 .eq. b2 ) then
                        m2(b1,b1)=v1(b1)
                    else
                        m2(b1,b2)=0.0_r8
                    endif
                enddo
            enddo
            ! Store as generalized group velocity
            do b1=1,nb
            do b2=1,nb
                if ( ompoint%omega(b1) .lt. lo_freqtol ) cycle
                if ( ompoint%omega(b2) .lt. lo_freqtol ) cycle
                ompoint%genvel(ialpha,b1,b2)=real( m2(b1,b2)/(ompoint%omega(b1)+ompoint%omega(b2)), r8)
            enddo
            enddo

            deallocate(m2)
        enddo


        ! Cleanup
        call mem%deallocate(Dq, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(D,  persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(m0, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(di, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        if (present(qvec)) then
            deallocate (dumqpoint%invariant_operation)
        end if
    end block slv
end subroutine

subroutine destroy_ompoint(self)
    class(lo_distributed_phonon_dispersions_qpoint), intent(inout) :: self

    if ( allocated(self%omega      ) ) deallocate(self%omega      )
    if ( allocated(self%vel        ) ) deallocate(self%vel        )
    if ( allocated(self%genvel     ) ) deallocate(self%genvel     )
    if ( allocated(self%egv        ) ) deallocate(self%egv        )
    if ( allocated(self%nuvec      ) ) deallocate(self%nuvec      )
    if ( allocated(self%degeneracy ) ) deallocate(self%degeneracy )
    if ( allocated(self%degenmode  ) ) deallocate(self%degenmode  )
    if ( allocated(self%sigma      ) ) deallocate(self%sigma      )
    select type(self)
    type is(lo_distributed_phonon_dispersions_qpoint_irr)
        self%global_irreducible_index=-lo_hugeint
        self%n_full_point=-lo_hugeint
        if ( allocated(self%local_full_index    ) ) deallocate(self%local_full_index    )
        if ( allocated(self%local_operation_full) ) deallocate(self%local_operation_full)
    type is(lo_distributed_phonon_dispersions_qpoint_full)
        self%global_full_index=-lo_hugeint
        self%local_irr_index=-lo_hugeint
        self%operation_from_irreducible=-lo_hugeint
    end select
end subroutine

subroutine destroy_dispersions(self)
    class(lo_distributed_phonon_dispersions), intent(inout) :: self

    integer :: i

    self%n_mode = -lo_hugeint
    self%n_irr_qpoint_local=-lo_hugeint
    self%n_full_qpoint_local=-lo_hugeint
    self%omega_max = -lo_huge
    self%omega_min = -lo_huge
    if ( allocated(self%iq) ) then
        do i=1,size(self%iq)
            call self%iq(i)%destroy()
        enddo
        deallocate(self%iq)
    endif
    if ( allocated(self%aq) ) then
        do i=1,size(self%aq)
            call self%aq(i)%destroy()
        enddo
        deallocate(self%aq)
    endif
end subroutine

end module