module phononevents
!!
use konstanter, only: r8, i8, lo_status, lo_tol, lo_sqtol, lo_freqtol, lo_huge, lo_hugeint, lo_tiny, &
                      lo_exitcode_symmetry, lo_exitcode_param, lo_frequency_THz_to_Hartree
use gottochblandat, only: lo_gauss, lo_stddev, lo_sqnorm, walltime, lo_identitymatrix, &
                          lo_progressbar_init, lo_progressbar, &
                          lo_nullspace_coefficient_matrix, lo_real_pseudoinverse
use mpi_wrappers, only: lo_mpi_helper, lo_stop_gracefully
use lo_memtracker, only: lo_mem_helper
use type_crystalstructure, only: lo_crystalstructure
use type_forceconstant_secondorder, only: lo_forceconstant_secondorder
use type_forceconstant_thirdorder, only: lo_forceconstant_thirdorder
use type_qpointmesh, only: lo_qpoint_mesh, lo_fft_mesh, lo_LV_tetrahedron_weights, lo_integration_weights_for_one_tetrahedron
use type_phonon_dispersions, only: lo_phonon_dispersions
use type_symmetryoperation, only: lo_operate_on_vector

implicit none
private
public :: lo_threephononevents
public :: lo_find_all_scattering_events

! Container for isotope integrationweights
type lo_iso4
    integer :: gi2 = -lo_hugeint
    real(r8) :: deltafunction, scatterstrength, W
end type
type lo_iso3
    integer :: n = -lo_hugeint
    type(lo_iso4), dimension(:), allocatable :: e
end type
type lo_iso2
    integer :: gi1 = -lo_hugeint
    type(lo_iso3), dimension(:, :), allocatable :: band
end type

! For three-phonon integrationweights
type lo_3phqp4
    ! On the Monkhorst-pack grid, which three points are involved?
    integer :: gi2 = -lo_hugeint
    integer :: gi3 = -lo_hugeint
    ! The integration weight, matrix element and the extra parameter to
    ! keep track of Pierls-Boltzmann equation stuff.
    real(r8) :: deltafunction = -lo_huge
    real(r8) :: psisquare = -lo_huge
    real(r8) :: W = -lo_huge
    ! if I have weird things, store frequencies and eigenvectors
    real(r8) :: omega3 = -lo_huge
    complex(r8), dimension(:), allocatable :: egv3
end type
type lo_3phqp3
    integer :: n = -lo_hugeint
    type(lo_3phqp4), dimension(:), allocatable :: e
end type
type lo_3phqp2
    type(lo_3phqp3), dimension(:, :, :), allocatable :: plus, minus
    integer :: gi1 = -lo_hugeint
end type

!> type to keep track of three-phonon and two-phonon events
type lo_threephononevents
    !> Number of q-points on this rank
    integer :: n_local_qpoint = -lo_hugeint
    ! Plus corresponds to
    ! omega-omega'-omega''=0 or om1=om2+om3 or om1-om3 = om2
    !
    ! Minus corresponds to
    ! omega-omega'+omega''=0 or om1=om3-om2 or om1-om3 = -om2
    type(lo_3phqp2), dimension(:), allocatable :: q
    !> Isotope events
    type(lo_iso2), dimension(:), allocatable :: iq
    !> where to cut Gaussian tails
    real(r8) :: gaussian_threshold = -lo_huge
    !> how to adjust the automatic smearing
    real(r8) :: smearing_prefactor
    !> How to integrate
    integer :: integrationtype = -lo_hugeint
    !> should I include isotope scattering?
    logical :: isotopescattering = .false.
    !> Boundary scattering another way
    real(r8) :: mfpmax = -lo_huge
    !> Correction level. Keep it at 2, that's safe. The others are experimental.
    integer :: correctionlevel = -lo_hugeint
end type
contains

#include "phononevents_gaussian.f90"
#include "phononevents_tetrahedron.f90"

!> Count all scattering events, threephonon and isotope
subroutine lo_find_all_scattering_events(sc, qp, dr, p, mw, mem, sigma, thres, integrationtype, &
                                         correctionlevel, mfpmax, isotopescattering)
    !> scattering events
    type(lo_threephononevents), intent(out) :: sc
    !> q-point mesh
    class(lo_qpoint_mesh), intent(in) :: qp
    !> phonon dispersions
    type(lo_phonon_dispersions), intent(in) :: dr
    !> structure
    type(lo_crystalstructure), intent(in) :: p
    !> mpi helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem
    !> smearing parameter
    real(r8), intent(in) :: sigma
    !> threshold for gaussians
    real(r8), intent(in) :: thres
    !> integration type
    integer, intent(in) :: integrationtype
    !> correction level
    integer, intent(in) :: correctionlevel
    !> max mfp for boundary scattering
    real(r8), intent(in) :: mfpmax
    !> consider isotope scattering
    logical, intent(in) :: isotopescattering

    real(r8) :: t0, t1

    t0 = walltime()
    t1 = t0

    ! Set some basic things
    init: block
        if (mw%talk) then
            write (*, *) ' '
            write (*, *) 'Counting scattering events and calculating integration weights'
        end if

        ! Store some parameters
        sc%gaussian_threshold = thres
        sc%smearing_prefactor = sigma
        sc%integrationtype = integrationtype
        sc%isotopescattering = isotopescattering
        sc%mfpmax = mfpmax
        sc%correctionlevel = correctionlevel

        ! what kind of integration?
        if (mw%talk) then
            select case (sc%integrationtype)
            case (1)
                write (*, *) '... fixed gaussian smearing'
            case (2)
                write (*, *) '... adaptive gaussian smearing, scalingfactor=', sc%smearing_prefactor
            case (3)
                write (*, *) '... tetrahedron integration'
            case default
                call lo_stop_gracefully(['Unknown integration type'], lo_exitcode_param, __FILE__, __LINE__, mw%comm)
            end select

            write (*, *) 'sc%gaussian_threshold', sc%gaussian_threshold
            write (*, *) 'sc%smearing_prefactor', sc%smearing_prefactor
            write (*, *) 'sc%integrationtype   ', sc%integrationtype
            write (*, *) 'sc%isotopescattering ', sc%isotopescattering
            write (*, *) 'sc%mfpmax            ', sc%mfpmax
            write (*, *) 'sc%correctionlevel   ', sc%correctionlevel

        end if
    end block init

    ! Start making a little space
    makespace: block
        integer :: i

        ! Count q-points on this rank.
        sc%n_local_qpoint = 0
        do i = 1, qp%n_irr_point
            if (mod(i, mw%n) .eq. mw%r) sc%n_local_qpoint = sc%n_local_qpoint + 1
        end do

        ! I only allocate the q-point for this rank
        if (sc%n_local_qpoint .gt. 0) then
            allocate (sc%q(sc%n_local_qpoint))
            allocate (sc%iq(sc%n_local_qpoint))
            do i = 1, sc%n_local_qpoint
                allocate (sc%q(i)%plus(dr%n_mode, dr%n_mode, dr%n_mode))
                allocate (sc%q(i)%minus(dr%n_mode, dr%n_mode, dr%n_mode))
                allocate (sc%iq(i)%band(dr%n_mode, dr%n_mode))
            end do
        end if
    end block makespace

    ! Start measuring phase space
    measurephasespace: block
        integer :: i, lqp, gi1

        if (mw%talk) call lo_progressbar_init()
        lqp = 0
        do i = 1, qp%n_irr_point
            ! Make it parallel
            if (mod(i, mw%n) .ne. mw%r) cycle
            ! Increment counter
            lqp = lqp + 1
            ! first index on the grid
            gi1 = qp%ip(i)%full_index

            ! Three-phonon scattering
            select case (sc%integrationtype)
            case (1:2) ! gaussian integration
                call threephonon_gaussian_oneqp(qp, dr, sc%q(lqp), gi1, sc%gaussian_threshold, &
                                                sc%smearing_prefactor, sc%integrationtype, mem)
            case (3) ! tetrahedron
                call threephonon_tetrahedron_oneqp(qp, dr, sc%q(lqp), gi1, mem)
            end select

            ! Isotope scattering
            if (sc%isotopescattering) then
                select case (sc%integrationtype)
                case (1:2) ! fix gaussian
                    call iso_gaussian_oneqp(qp, dr, sc%iq(lqp), gi1, sc%gaussian_threshold, &
                                            sc%smearing_prefactor, sc%integrationtype, mem)
                case (3) ! tetrahedron
                    call iso_tetrahedron_oneqp(qp, dr, sc%iq(lqp), gi1, mem)
                end select
            end if

            if (mw%talk .and. lqp .lt. sc%n_local_qpoint) then
                call lo_progressbar(' ... counting scattering events', lqp, sc%n_local_qpoint, walltime() - t0)
            end if
        end do
        call mw%barrier()
        t1 = walltime()
        if (mw%talk) then
            call lo_progressbar(' ... counting scattering events', sc%n_local_qpoint, sc%n_local_qpoint, t1 - t0)
        end if
        t0 = t1
    end block measurephasespace

    ! Some diagnostics and final things
    diagnostics: block
        integer(i8) :: ct_plus, ct_minus, ct_iso, ct_possible_3ph, ct_possible_iso
        real(r8) :: f0, f1, f2
        integer :: i, lqp, b1, b2, b3

        ct_plus = 0
        ct_minus = 0
        lqp = 0
        do i = 1, qp%n_irr_point
            ! Make it parallel
            if (mod(i, mw%n) .ne. mw%r) cycle
            ! Increment counter
            lqp = lqp + 1
            do b1 = 1, dr%n_mode
            do b2 = 1, dr%n_mode
            do b3 = 1, dr%n_mode
                ct_plus = ct_plus + sc%q(lqp)%plus(b1, b2, b3)%n
                ct_minus = ct_minus + sc%q(lqp)%minus(b1, b2, b3)%n
            end do
            end do
            end do
        end do

        ct_iso = 0
        lqp = 0
        do i = 1, qp%n_irr_point
            ! Make it parallel
            if (mod(i, mw%n) .ne. mw%r) cycle
            ! Increment counter
            lqp = lqp + 1
            do b1 = 1, dr%n_mode
            do b2 = 1, dr%n_mode
                ct_iso = ct_iso + sc%iq(lqp)%band(b1, b2)%n
            end do
            end do
        end do

        ! reduce it
        call mw%allreduce('sum', ct_plus)
        call mw%allreduce('sum', ct_minus)
        call mw%allreduce('sum', ct_iso)

        ! dump some info
        if (mw%talk) then
            ! calculate the possible number of events
            ct_possible_3ph = qp%n_irr_point*qp%n_full_point
            ct_possible_3ph = ct_possible_3ph*2*(dr%n_mode**3)
            ct_possible_iso = qp%n_irr_point*qp%n_full_point*dr%n_mode**2
            write (*, *) ''
            write (*, *) '   number of + events:', ct_plus
            write (*, *) '   number of - events:', ct_minus
            write (*, *) '   number of i events:', ct_iso
            write (*, *) '                % 3ph:', 100*real(ct_plus + ct_minus, r8)/real(ct_possible_3ph, r8)
            write (*, *) '                % iso:', 100*real(ct_iso, r8)/real(ct_possible_iso, r8)
        end if

        ! do i=1,qp%n_irr_point
        !
        ! enddo
        ! if ( mw%r .eq. 0 ) then
        !     totevents=qp%n_irr_point
        !     totevents=totevents*qp%n_full_point
        !     totevents=totevents*(dr%n_mode**3)*2
        !     write(*,"(1X,A,F6.2,A,E13.5,A,E13.5,A)") &
        !     '... found',100.0*(ctr_plus+ctr_minus)/totevents,&
        !     '% of threephonon events to be relevant,',ctr_plus,' +',ctr_minus,' events'
        !
        !
        !
        !     totevents=qp%n_full_point
        !     totevents=totevents*(dr%n_mode**2)
        !     totevents=totevents*qp%n_irr_point
        !
        !     write(*,"(1X,A,F6.2,A,E13.5,A)") &
        !     '... found',100.0*ctr_iso/totevents,'% of isotope events to be relevant, ',ctr_iso,' events'
        ! endif
    end block diagnostics

    ! ! For the more stringent corrections I need some more things
    ! if ( sc%correctionlevel .ge. 4 ) then
    ! symthings: block
    !     real(r8), dimension(3,3) :: I3,m0
    !     real(r8), dimension(3) :: v0,v1,v2
    !     real(r8) :: f0,f1
    !     integer, dimension(3) :: gi
    !     integer :: o
    !     ! Also build the q-point rotation scheme while I'm here.
    !     lo_allocate(sc%qprotp(uc%sym%n,qp%n_full_point))
    !     lo_allocate(sc%qprotn(uc%sym%n,qp%n_full_point))
    !     sc%qprotp=0
    !     sc%qprotn=0
    !     do i=1,qp%n_full_point
    !         v0=qp%ap(i)%r ! q-point
    !         do j=1,uc%sym%n
    !             if ( qp%operationok(j) .eqv. .false. ) cycle
    !             v1=lo_operate_on_vector(uc%sym%op(j),v0,reciprocal=.true.)
    !             v1=uc%cartesian_to_fractional(v1,reciprocal=.true.)
    !             v2=-lo_operate_on_vector(uc%sym%op(j),v0,reciprocal=.true.)
    !             v2=uc%cartesian_to_fractional(v2,reciprocal=.true.)
    !             select type(qp); type is(lo_fft_mesh)
    !                 if ( qp%is_point_on_grid(v1) ) then
    !                     gi=qp%index_from_coordinate(v1)
    !                     sc%qprotp(j,i)=qp%gridind2ind(gi(1),gi(2),gi(3))
    !                 else
    !                     !write(*,*) 'BAD symmetry operation, q-mesh not symmetric or something'
    !                 endif
    !                 if ( qp%is_point_on_grid(v2) ) then
    !                     gi=qp%index_from_coordinate(v2)
    !                     sc%qprotn(j,i)=qp%gridind2ind(gi(1),gi(2),gi(3))
    !                 else
    !                     !write(*,*) 'BAD symmetry operation, q-mesh not symmetric or something'
    !                 endif
    !             end select
    !         enddo
    !     enddo
    !
    !     ! Also neat to keep the nullspace matrices for each q-point
    !     call lo_identitymatrix(I3)
    !     lqp=0
    !     do i=1,qp%n_irr_point
    !         ! Make it parallel
    !         if ( mod(i,mw%n) .ne. mw%r ) cycle
    !         ! Increment counter
    !         lqp=lqp+1
    !         m0=0.0_r8
    !         do o=1,uc%sym%n
    !             if ( qp%operationok(o) .eqv. .false. ) cycle
    !             f0=0.0_r8
    !             f1=0.0_r8
    !             do j=1,dr%n_mode
    !                 v0=dr%iq(i)%vel(:,j)
    !                 v1=lo_operate_on_vector(uc%sym%op(o),v0,reciprocal=.true.)
    !                 f0=f0+lo_sqnorm(v0-v1)
    !                 f1=f1+lo_sqnorm(v0+v1)
    !             enddo
    !             if ( f0 .lt. dr%n_mode*lo_sqtol ) then
    !                 m0=m0+uc%sym%op(o)%m-I3
    !             endif
    !             if ( f1 .lt. dr%n_mode*lo_sqtol ) then
    !                 m0=m0-uc%sym%op(o)%m-I3
    !             endif
    !         enddo
    !         if ( sum(abs(m0)) .gt. lo_tol ) then
    !             call lo_nullspace_coefficient_matrix(m0,sc%q(lqp)%invariant_operator,sc%q(lqp)%nvar)
    !             ! So, how does one use this. Hmm. Call the matrix M.
    !             ! To get the least squares solution we get
    !             ! x = M+ v
    !             ! vs = M x
    !             ! vs = MM+ v
    !             ! seems easy enough?
    !             if ( sc%q(lqp)%nvar .gt. 0 ) then
    !                 allocate(sc%q(lqp)%pinv_invariant_operator(sc%q(lqp)%nvar,3))
    !                 call lo_real_pseudoinverse(sc%q(lqp)%invariant_operator,sc%q(lqp)%pinv_invariant_operator)
    !             endif
    !         else
    !             sc%q(lqp)%nvar=0
    !         endif
    !     enddo
    ! end block symthings
    ! endif
end subroutine

!> returns the index on the grid that gives q3=-q1-q2
pure function fft_third_grid_index(i1, i2, dims) result(i3)
    !> index to q1
    integer, intent(in) :: i1
    !> index to q2
    integer, intent(in) :: i2
    !> dimensions of the grid
    integer, dimension(3), intent(in) :: dims
    !> index to q3
    integer :: i3

    integer, dimension(3) :: gi1, gi2, gi3
    integer :: l, k

    ! Convert triplet to singlet
    gi1 = singlet_to_triplet(i1, dims(2), dims(3))
    gi2 = singlet_to_triplet(i2, dims(2), dims(3))
    do l = 1, 3
        gi3(l) = 3 - gi1(l) - gi2(l)
    end do
    do k = 1, 3
    do l = 1, 3
        if (gi3(l) .lt. 1) gi3(l) = gi3(l) + dims(l)
        if (gi3(l) .gt. dims(l)) gi3(l) = gi3(l) - dims(l)
    end do
    end do
    ! convert it back to a singlet
    i3 = triplet_to_singlet(gi3, dims(2), dims(3))

contains
    !> convert a linear index to a triplet
    pure function singlet_to_triplet(l, ny, nz) result(gi)
        !> linear index
        integer, intent(in) :: l
        !> second dimension
        integer, intent(in) :: ny
        !> third dimension
        integer, intent(in) :: nz
        !> grid-index
        integer, dimension(3) :: gi

        integer :: i, j, k

        k = mod(l, nz)
        if (k .eq. 0) k = nz
        j = mod((l - k)/nz, ny) + 1
        i = (l - k - (j - 1)*nz)/(nz*ny) + 1
        gi = [i, j, k]
    end function
    !> convert a triplet index to a singlet
    pure function triplet_to_singlet(gi, ny, nz) result(l)
        !> grid-index
        integer, dimension(3), intent(in) :: gi
        !> second dimension
        integer, intent(in) :: ny
        !> third dimension
        integer, intent(in) :: nz
        !> linear index
        integer :: l

        l = (gi(1) - 1)*ny*nz + (gi(2) - 1)*nz + gi(3)
    end function
end function

end module
