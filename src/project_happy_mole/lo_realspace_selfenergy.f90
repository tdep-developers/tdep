module lo_realspace_selfenergy
!!
!! Realspace representation of phonon self-energies
!!
use konstanter, only: r8, lo_iou, lo_hugeint, lo_huge, lo_sqtol, lo_exitcode_symmetry, lo_twopi, &
                      lo_bohr_to_A, lo_freqtol, lo_exitcode_param, lo_pi, lo_imag, lo_groupvel_ms_to_Hartreebohr
use gottochblandat, only: tochar, walltime, lo_progressbar_init, lo_progressbar, lo_chop, &
                          lo_linear_least_squares, lo_rsquare, lo_linspace, lo_linear_interpolation, lo_trapezoid_integration
use geometryfunctions, only: lo_inscribed_sphere_in_box
use mpi_wrappers, only: lo_mpi_helper, lo_stop_gracefully
use lo_memtracker, only: lo_mem_helper
use hdf5_wrappers, only: lo_hdf5_helper
use type_crystalstructure, only: lo_crystalstructure
use type_forceconstant_secondorder, only: lo_forceconstant_secondorder
use type_qpointmesh, only: lo_qpoint_mesh, lo_fft_mesh, lo_qpoint
use type_phonon_dispersions, only: lo_phonon_dispersions, lo_phonon_dispersions_qpoint
use lo_symmetry_of_interactions, only: lo_interaction_tensors
use type_symmetryoperation, only: lo_eigenvector_transformation_matrix
use type_blas_lapack_wrappers, only: lo_gemm, lo_dgels
implicit none

private
public :: lo_interpolated_selfenergy

type lo_interpolated_selfenergy_pair
    !> which atoms are involved?
    integer :: a1 = -lo_hugeint
    integer :: a2 = -lo_hugeint
    !> lattice vector?
    real(r8), dimension(3) :: lv = -lo_huge
    !> coefficient matrix
    real(r8), dimension(:, :), allocatable :: coeff
    !> index to irreducible ifc
    integer, dimension(:), allocatable :: xind
    !> how many irreducible?
    integer :: nx = -lo_hugeint
end type

type lo_interpolated_selfenergy
    !> how many pairs?
    integer :: n_pair = -lo_hugeint
    !> how many energies?
    integer :: n_energy = -lo_hugeint
    !> how many irreducible numbers?
    integer :: n_x = -lo_hugeint
    !> how many atoms
    integer :: n_atom = -lo_hugeint
    !> what temperature is it calculated for?
    real(r8) :: temperature = -lo_huge
    !> shift added for correct scaling
    real(r8) :: omega_shift = -lo_huge
    !> actual pairs
    type(lo_interpolated_selfenergy_pair), dimension(:), allocatable :: pair
    !> energy-dependent irreducible representation
    real(r8), dimension(:), allocatable :: energy
    !> energy-dependent irreducible representation
    real(r8), dimension(:, :), allocatable :: x_re
    real(r8), dimension(:, :), allocatable :: x_im
contains
    !> create interpolation from gridded self-energy
    procedure :: generate
    !> evaluate self-energy
    procedure :: diagonal_selfenergy
    !> evaluate self-energy and velocities
    procedure :: diagonal_selfenergy_and_velocity
    !> write to hdf5
    procedure :: write_to_hdf5
    !> read from hdf5
    procedure :: read_from_hdf5
end type

! Constant for now, to ensure sensible interpolation
real(r8), parameter :: omega_normalizing_shift = 1.E-1_r8 !1E-2_r8

contains

#include "lo_realspace_selfenergy_io.f90"

!> evaluate diagonal self-energy at arbitrary q-point. Can be made a lot faster.
subroutine diagonal_selfenergy_and_velocity(ise, p, fc, qpoint, ompoint, omega, egv, sigmaRe, sigmaIm, vel, mem)
    !> self-energy guy
    class(lo_interpolated_selfenergy), intent(in) :: ise
    !> structure
    type(lo_crystalstructure), intent(in) :: p
    !> forceconstant
    type(lo_forceconstant_secondorder), intent(inout) :: fc
    !> q-vector
    class(lo_qpoint), intent(in) :: qpoint
    !> harmonic properties
    type(lo_phonon_dispersions_qpoint), intent(in) :: ompoint
    !> frequencies
    real(r8), dimension(:), intent(in) :: omega
    !> phonon eigenvectors at this q-point
    complex(r8), dimension(:, :), intent(in) :: egv
    !> real part of self-energy (energy,mode)
    real(r8), dimension(:, :), intent(inout) :: sigmaRe
    !> imaginary part of self-energy (energy,mode)
    real(r8), dimension(:, :), intent(inout) :: sigmaIm
    !> energy-dependent group velocities (xyz,energy,mode,mode)
    complex(r8), dimension(:, :, :, :), intent(inout) :: vel
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem

    !fastversion: block
    complex(r8), dimension(:, :, :), allocatable :: gcm0, gcm1, grad_dynmat
    complex(r8), dimension(:, :), allocatable :: coeff, kronegv, cm0, cm1, cm2
    complex(r8), dimension(:, :), allocatable :: buf_egv, buf_egw
    complex(r8), dimension(:, :), allocatable :: buf_x, buf_y
    complex(r8), dimension(3, 3) :: cr0
    complex(r8), dimension(3) :: cv0
    integer :: nb, a1, a2, ix, iy, iz, ia, ib, ic, i, j, ii, jj, ie, iop, k

    ! Some buffers
    nb = ise%n_atom*3
    allocate (gcm0(nb**2, ise%n_x, 3))
    allocate (cm0(nb**2, ise%n_x))
    allocate (cm1(nb**2, ise%n_x))
    allocate (buf_egv(nb, nb))
    allocate (buf_egw(nb, nb))
    allocate (kronegv(nb**2, nb**2))
    allocate (coeff(nb, ise%n_x))

    ! First step is to Kronecker expand the eigenvectors. But first convert
    ! the eigenvectors to nu-vectors.
    kronegv = 0.0_r8
    do i = 1, nb
        if (omega(i) .gt. lo_freqtol) then
            do a1 = 1, ise%n_atom
            do ix = 1, 3
                ib = (a1 - 1)*3 + ix
                buf_egw(ib, i) = egv(ib, i)*p%invsqrtmass(a1)/sqrt(omega(i))
                !buf_egw(ib,i)=egv(ib,i)*p%invsqrtmass(a1)*sqrt( omega(i) )
            end do
            end do
        else
            buf_egw(:, i) = 0.0_r8
        end if
    end do
    ! Then the proper Kronecker expansion, compatible with the flattening in the coefficient matrix
    do i = 1, nb
    do j = 1, nb
        buf_egv = buf_egw*conjg(buf_egw(j, i))
        do ii = 1, nb
        do jj = 1, nb
            ia = (i - 1)*nb + ii
            ib = (j - 1)*nb + jj
            kronegv(ia, ib) = buf_egv(jj, ii)
        end do
        end do
    end do
    end do

    ! Obtain raw coefficient matrix
    call coefficient_matrix_per_q(ise, qpoint%r, cm0, gcm0)
    ! Apply the transformation to the coefficient matrix
    call lo_gemm(kronegv, cm0, cm1)
    ! Select the right lines for the self-energy, only want the diagonal part.
    do i = 1, nb
        j = (i - 1)*nb + i
        coeff(i, :) = cm1(j, :)
    end do
    ! Get the self-energies
    allocate (buf_x(ise%n_x, ise%n_energy))
    allocate (buf_y(ise%n_energy, nb))
    buf_x = ise%x_im
    call lo_gemm(buf_x, coeff, buf_y, transa='T', transb='T')
    sigmaIm = real(buf_y, r8)

    buf_x = ise%x_re
    call lo_gemm(buf_x, coeff, buf_y, transa='T', transb='T')
    sigmaRe = real(buf_y, r8)

    ! A little cleanup
    deallocate (cm0)
    deallocate (cm1)

    ! That sorted out the self-energies, now tackle the group velocities. This can get
    ! a little hairy, but bear with me. Think it will work. I want to full-full thing
    ! so first get the derivative of the dynamical matrix.

    allocate (grad_dynmat(nb, nb, 3))
    allocate (cm0(nb, nb))
    allocate (cm2(nb**2, 3))
    grad_dynmat = 0.0_r8
    cm0 = 0.0_r8
    cm2 = 0.0_r8
    ! Get dynamical matrix?
    call fc%dynamicalmatrix(p, qpoint, cm0, mem, grad_dynmat, qdirection=[1.0_r8, 0.0_r8, 0.0_r8])
    ! Flatten and scale with masses
    do iz = 1, 3
        do a1 = 1, ise%n_atom
        do a2 = 1, ise%n_atom
        do ix = 1, 3
        do iy = 1, 3
            ib = (a1 - 1)*3 + ix
            ia = (a2 - 1)*3 + iy
            ic = flattenind(a1, a2, ix, iy, nb)
            cm2(ic, iz) = grad_dynmat(ia, ib, iz)/(p%invsqrtmass(a1)*p%invsqrtmass(a2))
        end do
        end do
        end do
        end do
    end do

    deallocate (grad_dynmat)

    ! Convert from irreducible to proper
    allocate (gcm1(nb**2, ise%n_energy, 3))
    gcm1 = 0.0_r8
    buf_x = ise%x_re
    do iz = 1, 3
        call lo_gemm(gcm0(:, :, iz), buf_x, gcm1(:, :, iz), beta=cmplx(1.0_r8, 0.0_r8, r8))
    end do
    buf_x = ise%x_im
    do iz = 1, 3
        call lo_gemm(gcm0(:, :, iz), buf_x, gcm1(:, :, iz), beta=cmplx(0.0_r8, 1.0_r8, r8))
    end do
    ! Add the harmonic
    do iz = 1, 3
    do ie = 1, ise%n_energy
        gcm1(:, ie, iz) = gcm1(:, ie, iz) + cm2(:, iz)
    end do
    end do

    ! Average over the group. Can do this quite neatly I think.
    kronegv = 0.0_r8
    do k = 1, qpoint%n_invariant_operation
        iop = qpoint%invariant_operation(k)
        call lo_eigenvector_transformation_matrix(cm0, p%rcart, qpoint%r, p%sym%op(abs(iop)))
        if (iop .lt. 0) then
            !buf_egv=conjg(ompoint%egv)
            call lo_gemm(cm0, ompoint%egv, buf_egw)
            buf_egw = conjg(buf_egw)
        else
            call lo_gemm(cm0, ompoint%egv, buf_egw)
            !buf_egv=ompoint%egv
        end if
        !call lo_gemm(cm0,buf_egv,buf_egw)
        do i = 1, nb
            if (omega(i) .gt. lo_freqtol) then
                do a1 = 1, ise%n_atom
                do ix = 1, 3
                    ib = (a1 - 1)*3 + ix
                    buf_egw(ib, i) = buf_egw(ib, i)*p%invsqrtmass(a1)/sqrt(omega(i)*2.0_r8)
                end do
                end do
            else
                buf_egw(:, i) = 0.0_r8
            end if
        end do
        do i = 1, nb
        do j = 1, nb
            buf_egv = buf_egw*conjg(buf_egw(j, i))
            do ii = 1, nb
            do jj = 1, nb
                ia = (i - 1)*nb + ii
                ib = (j - 1)*nb + jj
                kronegv(ia, ib) = kronegv(ia, ib) + buf_egv(jj, ii)
            end do
            end do
        end do
        end do
    end do
    kronegv = kronegv/real(qpoint%n_invariant_operation, r8)

    deallocate (gcm0)
    ! And finally convert!
    allocate (gcm0(nb**2, ise%n_energy, 3))
    gcm0 = 0.0_r8
    do iz = 1, 3
        call lo_gemm(kronegv, gcm1(:, :, iz), gcm0(:, :, iz))
    end do

    ! Return in a sensible format.
    vel = 0.0_r8
    do ie = 1, ise%n_energy
        do i = 1, nb
        do j = 1, nb
            ii = (i - 1)*nb + j
            jj = (j - 1)*nb + i
            cv0 = gcm0(ii, ie, :) !+conjg(gcm0(jj,ie,:)) )*0.5_r8
            !cv0=matmul(ompoint%invariance_fixer,cv0)
            vel(:, ie, i, j) = cv0 !gcm0(ii,ie,:)
            !if ( i .ne. j ) vel(:,ie,j,i)=conjg(cv0)
        end do
        end do
    end do

    ! Make sure it is invariant as it should be?
    deallocate (cm0)
    allocate (cm0(3, ise%n_energy))

    do j = 1, nb
    do i = 1, nb
        cm0 = 0.0_r8
        do k = 1, qpoint%n_invariant_operation
            iop = qpoint%invariant_operation(k)
            cr0 = p%sym%op(iop)%m
            call lo_gemm(cr0, vel(:, :, i, j), cm0, beta=cmplx(1.0_r8, 0.0_r8, r8))
        end do
        vel(:, :, i, j) = cm0/real(qpoint%n_invariant_operation, r8)
    end do
    end do
    vel = lo_chop(vel, lo_groupvel_ms_to_Hartreebohr*1E-7_r8)

!
! write(*,*) ''
! write(*,*) ompoint%invariance_fixer
!
!         do ie=100,130
!             do i=1,nb
!             do j=1,nb
!                 ii=(i-1)*nb+j
!                 if ( i .eq. j .and. i .eq. 1 ) then
! write(*,*) ie,i,aimag(gcm0(ii,ie,:))*1000
!                 endif
!             enddo
!             enddo
!         enddo

    !
    ! do i=1,3
    !     cm2(:,i)=matmul(kronegv,cm2(:,i))
    ! enddo
    ! do ix=1,3
    ! do i=1,nb
    ! do j=1,nb
    !     ii=(i-1)*nb+j
    !     grad_dynmat(i,j,ix)=cm2(ii,ix)
    ! enddo
    ! enddo
    ! enddo
    !
    ! write(*,*) ''
    ! do i=1,nb
    !     write(*,*) real(grad_dynmat(i,i,:),r8)*1000
    ! enddo
    ! write(*,*) ''
    ! do i=1,nb
    !     write(*,*) aimag(grad_dynmat(i,i,:))
    ! enddo

    !end block fastversion

end subroutine

!> Get the spectral functions per mode, somehow
subroutine generate(ise, p, dr, qvec, wp, maxe, sigmaRe, sigmaIm, cutoff, temperature, mw, mem, verbosity)
    !> list of spectral functions
    class(lo_interpolated_selfenergy), intent(out) :: ise
    !> structure
    type(lo_crystalstructure), intent(inout) :: p
    !> dispersions on grid
    type(lo_phonon_dispersions), intent(in) :: dr
    !> q-vectors
    real(r8), dimension(:, :), intent(in) :: qvec
    !> qpoint mesh
    type(lo_phonon_dispersions_qpoint), dimension(:), intent(in) :: wp
    !> maximum energy
    real(r8), intent(in) :: maxe
    !> real part of self-energy (energy,mode,q-point)
    real(r8), dimension(:, :, :), intent(in) :: sigmaRe
    !> imaginary part of self-energy (energy,mode,q-point)
    real(r8), dimension(:, :, :), intent(in) :: sigmaIm
    !> cutoff for the interpolation thing
    real(r8), intent(in) :: cutoff
    !> make a note about which temperature we are doing
    real(r8), intent(in) :: temperature
    !> MPI communicator
    type(lo_mpi_helper), intent(inout) :: mw
    !> Memory tracker
    type(lo_mem_helper), intent(inout) :: mem
    !> how much to talk
    integer, intent(in) :: verbosity

    integer :: outeriter, nq
    real(r8), dimension(:, :), allocatable :: relerr

    init: block
        if (verbosity .gt. 0) then
            write (lo_iou, *) ''
            write (lo_iou, *) 'BUILDING INTERPOLATED SELFENERGY'
        end if

        ! Number of q-points
        nq = size(wp)

        write (*, *) 'nq', nq

        allocate (relerr(dr%n_mode, nq))
        relerr = 1.0_r8

    end block init

    ! Generate symmetry-related stuff. Can treat this like it's a
    ! second order forceconstant but energy-dependent.
    gensym: block
        type(lo_interaction_tensors) :: slt
        type(lo_crystalstructure) :: ss
        real(r8), dimension(3) :: v0
        integer :: ipair, ishell, iop

        ! Generate symmetry stuff. First check the supercell the q-mesh defines:
        ! select type(qp)
        ! type is(lo_fft_mesh)
        !     do i=1,3
        !         m0(:,i)=p%latticevectors(:,i)*qp%griddensity(i)
        !     enddo
        !     if ( verbosity .gt. 0 ) write(*,*) '... max cutoff:',tochar(lo_inscribed_sphere_in_box(m0)*lo_bohr_to_A),'A'
        ! ! class default
        ! !     call lo_stop_gracefully(['Only use FFT meshes for now.'],lo_exitcode_param,__FILE__,__LINE__,mw%comm)
        ! end select

        ss = p
        call ss%classify('supercell', p)
        call slt%generate(p, ss, cutoff, -1.0_r8, -1.0_r8, .false., mw, mem, verbosity, &
                          transposition=.true., spacegroup=.true., wraparound=.true.)

        ! Set up the skeleton?
        ise%n_pair = slt%n_uc_pair
        allocate (ise%pair(ise%n_pair))
        do ipair = 1, ise%n_pair
            ! Indices and lattice vector?
            ise%pair(ipair)%a1 = slt%uc_pair(ipair)%i1
            ise%pair(ipair)%a2 = slt%uc_pair(ipair)%i2
            v0 = slt%uc_pair(ipair)%v - p%rcart(:, ise%pair(ipair)%a2) + p%rcart(:, ise%pair(ipair)%a1)
            v0 = matmul(p%inv_latticevectors, v0)
            if (sum(abs(v0 - anint(v0))) .gt. lo_sqtol) then
                call lo_stop_gracefully(['Clearly I do not understand vectors'], lo_exitcode_symmetry, __FILE__, __LINE__, mw%comm)
            else
                ise%pair(ipair)%lv = lo_chop(matmul(p%latticevectors, anint(v0)), 1E-12_r8)
            end if

            ! Coefficient things
            ishell = slt%uc_pair(ipair)%unique
            iop = slt%uc_pair(ipair)%operation
            ise%pair(ipair)%nx = slt%fc_pair_shell(ishell)%nx

            if (ise%pair(ipair)%nx .gt. 0) then
                allocate (ise%pair(ipair)%coeff(9, ise%pair(ipair)%nx))
                allocate (ise%pair(ipair)%xind(ise%pair(ipair)%nx))
                ise%pair(ipair)%coeff = matmul(slt%pairop(iop)%sotr, slt%fc_pair_shell(ishell)%coeff)
                ise%pair(ipair)%xind = slt%fc_pair_shell(ishell)%ind_global
            else
                ! Dummy allocation
                allocate (ise%pair(ipair)%coeff(1, 1))
                allocate (ise%pair(ipair)%xind(1))
                ise%pair(ipair)%coeff = -lo_huge
                ise%pair(ipair)%xind = -lo_hugeint
            end if
        end do

        ! Store generic things
        ise%n_energy = size(sigmaRe, 1)
        ise%n_x = slt%nx_fc_pair
        ise%n_atom = p%na
        ise%temperature = temperature
        allocate (ise%energy(ise%n_energy))
        allocate (ise%x_re(ise%n_x, ise%n_energy))
        allocate (ise%x_im(ise%n_x, ise%n_energy))
        call lo_linspace(0.0_r8, maxe, ise%energy)
        ise%x_re = 0.0_r8
        ise%x_im = 0.0_r8

        ! Select a sensible omega shift
        ise%omega_shift = dr%omega_max*omega_normalizing_shift

        if (verbosity .gt. 0) then
            write (lo_iou, *) '... sorted out symmetry'
        end if
    end block gensym

    do outeriter = 1, 30

        ! set up equations for the irreducible and solve.
        buildeq: block
            complex(r8), dimension(:, :), allocatable :: buf_CM0, buf_egv
            complex(r8), dimension(:, :), allocatable :: buf_dm0, buf_dm1, buf_dm2
            real(r8), dimension(:, :), allocatable :: buf_cm1, buf_cm2, buf_b1, buf_b2
            real(r8), dimension(:), allocatable :: buf_wt
            integer :: nb, solrnk
            integer :: iq, i, j, ie
            integer :: a1, a2, ix, iy, ia, ib, ic

            nb = dr%n_mode**2
            solrnk = mw%n - 1

            call mem%allocate(buf_cm0, [nb, ise%n_x], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%allocate(buf_cm1, [nq*2*nb, ise%n_x], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%allocate(buf_cm2, [nq*2*nb, ise%n_x], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%allocate(buf_b1, [nq*2*nb, ise%n_energy], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%allocate(buf_b2, [nq*2*nb, ise%n_energy], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%allocate(buf_dm0, [dr%n_mode, dr%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%allocate(buf_dm1, [dr%n_mode, dr%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%allocate(buf_dm2, [dr%n_mode, dr%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%allocate(buf_egv, [dr%n_mode, dr%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%allocate(buf_wt, nq, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            buf_cm0 = 0.0_r8
            buf_cm1 = 0.0_r8
            buf_cm2 = 0.0_r8
            buf_b1 = 0.0_r8
            buf_b2 = 0.0_r8
            buf_dm0 = 0.0_r8
            buf_dm1 = 0.0_r8
            buf_dm2 = 0.0_r8
            buf_egv = 0.0_r8

            buf_wt = 1.0_r8
            ! Get weights per q-point?
            do iq = 1, nq
                ! f1=0.0_r8
                ! do ia=1,dr%n_mode
                !     if ( dr%iq(iq)%omega(ia) .gt. lo_freqtol ) then
                !         !f0=1.0_r8/lo_trapezoid_integration(ise%energy,sigmaIm(:,ia,iq))
                !         f0=lo_trapezoid_integration(ise%energy,sigmaIm(:,ia,iq))
                !     else
                !         f0=0.0_r8
                !     endif
                !     f1=f1+f0
                ! enddo
                buf_wt(iq) = norm2(relerr(:, iq)) !*qp%ip(iq)%integration_weight
                !write(*,*) iq,buf_wt(iq),f1
            end do
            buf_wt = sqrt(buf_wt)
            buf_wt = buf_wt/sum(buf_wt)

            ! Build huge coefficient matrix
            do iq = 1, nq !qp%n_irr_point
                ! Get coefficient for this q-vector
                call coefficient_matrix_per_q(ise, qvec(:, iq), buf_cm0)
                buf_cm0 = buf_cm0*buf_wt(iq)
                ! Accumulate to large matrix
                do i = 1, nb
                    j = (iq - 1)*2*nb + i
                    buf_cm1(j, :) = real(buf_cm0(i, :), r8)
                    j = (iq - 1)*2*nb + nb + i
                    buf_cm1(j, :) = aimag(buf_cm0(i, :))
                end do
            end do
            buf_cm2 = buf_cm1

            if (verbosity .gt. 0) then
                write (lo_iou, *) '... built coefficient'
            end if

            ! Build large thing to match
            do ie = 1, ise%n_energy
                if (mod(ie, mw%n) .ne. mw%r) cycle
                do iq = 1, nq
                    buf_dm0 = 0.0_r8
                    buf_dm1 = 0.0_r8
                    do i = 1, dr%n_mode
                        buf_dm0(i, i) = sigmaIm(ie, i, iq)
                        buf_dm1(i, i) = sigmaRe(ie, i, iq)
                    end do
                    ! rotate away from normal mode basis to Cartesian. For that I need
                    ! the inverse matrix of nu-vectors
                    buf_egv = wp(iq)%egv
                    do ia = 1, dr%n_mode
                        if (wp(iq)%omega(ia) .gt. lo_freqtol) then
                            do a1 = 1, p%na
                            do ix = 1, 3
                                ib = (a1 - 1)*3 + ix
                                !buf_egv(ib,ia)=buf_egv(ib,ia)*sqrt( dr%iq(iq)%omega(ia)+1E-7_r8 )/p%invsqrtmass(a1)
                                buf_egv(ib, ia) = buf_egv(ib, ia)*sqrt(wp(iq)%omega(ia) + ise%omega_shift)/p%invsqrtmass(a1)
                                !buf_egv(ib,ia)=buf_egv(ib,ia)/sqrt( dr%iq(iq)%omega(ia) )/p%invsqrtmass(a1)
                            end do
                            end do
                        else
                            buf_egv(:, ia) = 0.0_r8
                        end if
                    end do
                    ! Actual transformation
                    call lo_gemm(buf_egv, buf_dm0, buf_dm2)
                    call lo_gemm(buf_dm2, buf_egv, buf_dm0, transb='C')
                    call lo_gemm(buf_egv, buf_dm1, buf_dm2)
                    call lo_gemm(buf_dm2, buf_egv, buf_dm1, transb='C')

                    buf_dm0 = buf_dm0*buf_wt(iq)
                    buf_dm1 = buf_dm1*buf_wt(iq)

                    ! Store in the right place
                    do a1 = 1, p%na
                    do a2 = 1, p%na
                    do ix = 1, 3
                    do iy = 1, 3
                        ib = (a1 - 1)*3 + ix
                        ia = (a2 - 1)*3 + iy

                        ic = flattenind(a1, a2, ix, iy, 3*p%na)

                        j = (iq - 1)*2*nb + ic
                        buf_b1(j, ie) = real(buf_dm0(ia, ib), r8)
                        buf_b2(j, ie) = real(buf_dm1(ia, ib), r8)

                        j = (iq - 1)*2*nb + nb + ic
                        buf_b1(j, ie) = aimag(buf_dm0(ia, ib))
                        buf_b2(j, ie) = aimag(buf_dm1(ia, ib))
                    end do
                    end do
                    end do
                    end do
                end do
            end do
            call mw%allreduce('sum', buf_b1)
            call mw%allreduce('sum', buf_b2)

            if (verbosity .gt. 0) then
                write (lo_iou, *) '... built target'
            end if

            ! Solve for reals
            if (mw%r .eq. solrnk) then
                call lo_dgels(buf_cm1, buf_b1)
                call lo_dgels(buf_cm2, buf_b2)
                ! Store solution in the right place.
                do ie = 1, ise%n_energy
                    ise%x_im(1:ise%n_x, ie) = buf_b1(1:ise%n_x, ie)
                    ise%x_re(1:ise%n_x, ie) = buf_b2(1:ise%n_x, ie)
                end do
            else
                ise%x_im = 0.0_r8
                ise%x_re = 0.0_r8
            end if
            call mw%bcast(ise%x_im, from=solrnk)
            call mw%bcast(ise%x_re, from=solrnk)

            if (verbosity .gt. 0) then
                write (lo_iou, *) '... solved for irreducible'
            end if

            call mem%deallocate(buf_cm0, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%deallocate(buf_cm1, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%deallocate(buf_cm2, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%deallocate(buf_b1, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%deallocate(buf_b2, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%deallocate(buf_dm0, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%deallocate(buf_dm1, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%deallocate(buf_dm2, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%deallocate(buf_egv, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%deallocate(buf_wt, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        end block buildeq

        ! Here I can insert a block that does some sanity checks
        ! and massages things a little.

        ! Calculate an R^2 to see if things turned out alright.
        verify: block
            real(r8), dimension(:, :, :), allocatable :: bRe, bIm
            real(r8) :: f0, f1, f2, f3
            integer :: iq, i

            allocate (bRe(ise%n_energy, dr%n_mode, nq))
            allocate (bIm(ise%n_energy, dr%n_mode, nq))
            bRe = 0.0_r8
            bIm = 0.0_r8
            do iq = 1, nq
                if (mod(iq, mw%n) .ne. mw%r) cycle
                call diagonal_selfenergy(ise, p, qvec(:, iq), wp(iq)%omega, wp(iq)%egv, bRe(:, :, iq), bIm(:, :, iq))
            end do
            call mw%allreduce('sum', bRe)
            call mw%allreduce('sum', bIm)

            ! Calculate the relative error
            f2 = 0.0_r8
            f3 = 0.0_r8
            do iq = 1, nq
            do i = 1, dr%n_mode
                if (wp(iq)%omega(i) .gt. lo_freqtol) then
                    f0 = lo_linear_interpolation(ise%energy, sigmaIm(:, i, iq), wp(iq)%omega(i))
                    f1 = lo_linear_interpolation(ise%energy, bIm(:, i, iq), wp(iq)%omega(i))
! if ( mw%talk ) then
!     write(*,*) iq,i,f0,f1
! endif
                    if (f0 .gt. lo_freqtol) then
                        f2 = f2 + ((f0 - f1)/f0)**2
                        f3 = max(f3, f2)
                        relerr(i, iq) = (relerr(i, iq)*0.6_r8 + 0.4_r8*((f0 - f1)/f0)**2)
                    else
                        relerr(i, iq) = 0.0_r8
                    end if
                else
                    relerr(i, iq) = 0.0_r8
                end if
            end do
            end do

            if (verbosity .gt. 0) then
                write (lo_iou, *) '... evaluated'
                write (lo_iou, *) '... Im:', lo_rsquare(bIm, sigmaIm)
                write (lo_iou, *) '... Re:', lo_rsquare(bRe, sigmaRe)
                write (lo_iou, *) '... err', f2/nq/dr%n_mode, f3
            end if

        end block verify

    end do
end subroutine

!> evaluate diagonal self-energy at arbitrary q-point. Can be made a lot faster.
subroutine diagonal_selfenergy(ise, p, qvec, omega, egv, sigmaRe, sigmaIm)
    !> self-energy guy
    class(lo_interpolated_selfenergy), intent(in) :: ise
    !> structure
    type(lo_crystalstructure), intent(in) :: p
    !> q-vector
    real(r8), dimension(3), intent(in) :: qvec
    !> frequencies
    real(r8), dimension(:), intent(in) :: omega
    !> phonon eigenvectors at this q-point
    complex(r8), dimension(:, :), intent(in) :: egv
    !> real part of self-energy (energy,mode)
    real(r8), dimension(:, :), intent(inout) :: sigmaRe
    !> imaginary part of self-energy (energy,mode)
    real(r8), dimension(:, :), intent(inout) :: sigmaIm

    fastversion: block
        complex(r8), dimension(:, :), allocatable :: coeff, kronegv, cm0, cm1
        complex(r8), dimension(:, :), allocatable :: buf_egv, buf_egw
        complex(r8), dimension(:, :), allocatable :: buf_x, buf_y
        integer :: nb, a1, ix, ia, ib, i, j, ii, jj

        ! Some buffers
        nb = ise%n_atom*3
        allocate (cm0(nb**2, ise%n_x))
        allocate (cm1(nb**2, ise%n_x))
        allocate (buf_egv(nb, nb))
        allocate (buf_egw(nb, nb))
        allocate (kronegv(nb**2, nb**2))
        allocate (coeff(nb, ise%n_x))

        ! First step is to Kronecker expand the eigenvectors. But first convert
        ! the eigenvectors to nu-vectors.
        kronegv = 0.0_r8
        do i = 1, nb
            if (omega(i) .gt. lo_freqtol) then
                do a1 = 1, ise%n_atom
                do ix = 1, 3
                    ib = (a1 - 1)*3 + ix
                    buf_egw(ib, i) = egv(ib, i)*p%invsqrtmass(a1)/sqrt(omega(i) + ise%omega_shift)
                    !buf_egw(ib,i)=egv(ib,i)*p%invsqrtmass(a1)*sqrt( omega(i) )
                end do
                end do
            else
                buf_egw(:, i) = 0.0_r8
            end if
        end do
        ! Then the proper Kronecker expansion, compatible with the flattening in the coefficient matrix
        do i = 1, nb
        do j = 1, nb
            buf_egv = buf_egw*conjg(buf_egw(j, i))
            do ii = 1, nb
            do jj = 1, nb
                ia = (i - 1)*nb + ii
                ib = (j - 1)*nb + jj
                kronegv(ia, ib) = buf_egv(jj, ii)
            end do
            end do
        end do
        end do

        ! Obtain raw coefficient matrix
        call coefficient_matrix_per_q(ise, qvec, cm0)
        ! Apply the transformation to the coefficient matrix
        call lo_gemm(kronegv, cm0, cm1)
        ! Select the right lines
        do i = 1, nb
            j = (i - 1)*nb + i
            coeff(i, :) = cm1(j, :)
        end do
        ! Get the self-energies
        allocate (buf_x(ise%n_x, ise%n_energy))
        allocate (buf_y(ise%n_energy, nb))
        buf_x = ise%x_im
        call lo_gemm(buf_x, coeff, buf_y, transa='T', transb='T')
        sigmaIm = real(buf_y, r8)
        buf_x = ise%x_re
        call lo_gemm(buf_x, coeff, buf_y, transa='T', transb='T')
        sigmaRe = real(buf_y, r8)
    end block fastversion

    ! slowversion: block
    !     complex(r8), dimension(:,:), allocatable :: coeff,bRe,bIm
    !     complex(r8), dimension(:,:), allocatable :: buf_dm0,buf_dm1,buf_dm2
    !     integer :: nb,ie
    !     integer :: a1,a2,ix,iy,ia,ib,ic
    !
    !     ! Some buffers
    !     nb=ise%n_atom*3
    !     allocate(coeff(nb**2,ise%n_x))
    !     allocate(bRe(nb**2,ise%n_energy))
    !     allocate(bIm(nb**2,ise%n_energy))
    !     allocate(buf_dm0(nb,nb))
    !     allocate(buf_dm1(nb,nb))
    !     allocate(buf_dm2(nb,nb))
    !
    !     ! Build the coefficient matrix
    !     call coefficient_matrix_per_q(ise,qvec,coeff)
    !     bRe=matmul(coeff,ise%x_re)
    !     bIm=matmul(coeff,ise%x_im)
    !     ! Apply rotation thingy
    !     sigmaRe=0.0_r8
    !     sigmaIm=0.0_r8
    !     do ie=1,ise%n_energy
    !         ! Rearrange slightly
    !         do a1=1,ise%n_atom
    !         do a2=1,ise%n_atom
    !         do ix=1,3
    !         do iy=1,3
    !             ib=(a1-1)*3+ix
    !             ia=(a2-1)*3+iy
    !             ic=flattenind(a1,a2,ix,iy,nb)
    !             buf_dm0(ia,ib)=bRe(ic,ie)
    !             buf_dm1(ia,ib)=bIm(ic,ie)
    !         enddo
    !         enddo
    !         enddo
    !         enddo
    !
    !         ! Rotate
    !         call lo_gemm(egv,buf_dm0,buf_dm2,transa='C')
    !         call lo_gemm(buf_dm2,egv,buf_dm0)
    !         call lo_gemm(egv,buf_dm1,buf_dm2,transa='C')
    !         call lo_gemm(buf_dm2,egv,buf_dm1)
    !
    !         ! Store
    !         do ia=1,nb
    !             sigmaRe(ie,ia)=real( buf_dm0(ia,ia),r8)
    !             sigmaIm(ie,ia)=real( buf_dm1(ia,ia),r8)
    !         enddo
    !     enddo
    !
    !     deallocate(coeff)
    !     deallocate(bRe)
    !     deallocate(bIm)
    !     deallocate(buf_dm0)
    !     deallocate(buf_dm1)
    !     deallocate(buf_dm2)
    ! end block slowversion
end subroutine

!> get the coefficient matrix for a q-point?
subroutine coefficient_matrix_per_q(ise, qvec, coeff, gradcoeff)
    !> self-energy guy
    type(lo_interpolated_selfenergy), intent(in) :: ise
    !> q-vector
    real(r8), dimension(3), intent(in) :: qvec
    !> coefficient matrix (complex)
    complex(r8), dimension(:, :), intent(out) :: coeff
    !> coefficient matrix for derivative w.r.t q (complex)
    complex(r8), dimension(:, :, :), intent(out), optional :: gradcoeff

    real(r8), dimension(3) :: v0
    real(r8) :: qdotr
    complex(r8) :: expiqr
    integer :: ipair, a1, a2, ix, iy, ia, ib, ic, id, imu, nb, iz

    nb = ise%n_atom*3
    coeff = 0.0_r8

    if (present(gradcoeff)) then
        gradcoeff = 0.0_r8
    end if

    do ipair = 1, ise%n_pair
        if (ise%pair(ipair)%nx .le. 0) cycle

        ! Complex phase
        v0 = ise%pair(ipair)%lv
        qdotr = dot_product(v0, qvec)*lo_twopi
        expiqr = cmplx(cos(qdotr), sin(qdotr), r8)

        ! Store in the right place?
        a1 = ise%pair(ipair)%a1
        a2 = ise%pair(ipair)%a2
        do iy = 1, 3
        do ix = 1, 3
            ia = (iy - 1)*3 + ix
            ib = flattenind(a1, a2, ix, iy, nb)
            id = flattenind(a1, a1, ix, iy, nb)
            do imu = 1, ise%pair(ipair)%nx
                ic = ise%pair(ipair)%xind(imu)
                ! First line is normal, second line translational invariance
                coeff(ib, ic) = coeff(ib, ic) + ise%pair(ipair)%coeff(ia, imu)*expiqr
                coeff(id, ic) = coeff(id, ic) - ise%pair(ipair)%coeff(ia, imu)
            end do
        end do
        end do

        if (present(gradcoeff)) then
            do iz = 1, 3
            do iy = 1, 3
            do ix = 1, 3
                ia = (iy - 1)*3 + ix
                ib = flattenind(a1, a2, ix, iy, nb)
                id = flattenind(a1, a1, ix, iy, nb)
                do imu = 1, ise%pair(ipair)%nx
                    ic = ise%pair(ipair)%xind(imu)
                    gradcoeff(ib, ic, iz) = gradcoeff(ib, ic, iz) + lo_imag*v0(iz)*ise%pair(ipair)%coeff(ia, imu)*expiqr
                end do
            end do
            end do
            end do
        end if
    end do
end subroutine

! Consistent index flattening? Impossibru to get consistent.
function flattenind(a1, a2, ix, iy, nb) result(i)
    integer, intent(in) :: a1, a2, ix, iy, nb
    integer :: i

    integer :: ia, ib

    ia = (a1 - 1)*3 + ix
    ib = (a2 - 1)*3 + iy
    i = (ib - 1)*nb + ia
end function

end module
