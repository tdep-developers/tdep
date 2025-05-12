module scatteringrates
use konstanter, only: r8, lo_freqtol, lo_twopi, lo_sqtol
use gottochblandat, only: tochar, walltime, lo_trueNtimes, lo_progressbar_init, lo_progressbar, lo_sqnorm, lo_clean_fractional_coordinates
use mpi_wrappers, only: lo_mpi_helper, lo_stop_gracefully
use lo_memtracker, only: lo_mem_helper
use lo_timetracker, only: lo_timer
use type_crystalstructure, only: lo_crystalstructure
use type_forceconstant_secondorder, only: lo_forceconstant_secondorder
use type_forceconstant_thirdorder, only: lo_forceconstant_thirdorder
use type_forceconstant_fourthorder, only: lo_forceconstant_fourthorder
use type_qpointmesh, only: lo_qpoint_mesh, lo_qpoint, lo_fft_mesh
use type_phonon_dispersions, only: lo_phonon_dispersions, lo_phonon_dispersions_qpoint
use lineshape_helper, only: index_on_grid
implicit none

private
public :: lo_listofscatteringrates

type lo_listofscatteringrates
    !> frequency at q''
    !real(r8), dimension(:, :), allocatable :: omega3
    !> group velocity at q''
    !real(r8), dimension(:, :, :), allocatable :: vel3
    !> three-phonon matrix elements, not squared (mode1,mode2,mode3,q)
    complex(r8), dimension(:, :, :, :), allocatable :: psi_3ph
    !> isotope matrix elements (mode1,mode2,q)
    real(r8), dimension(:, :, :), allocatable :: psi_iso
    !> make a note if q1 is gamma
    logical :: atgamma = .false.

    !> q-vectors are needed for evaluation of spectral functions
    real(r8), dimension(:, :), allocatable :: qvec3
    !> eigenvectors are needed for spectral functions
    !complex(r8), dimension(:, :, :), allocatable :: egv3
    !> harmonic properties at the third q-vector
    type(lo_phonon_dispersions_qpoint), dimension(:), allocatable :: wp3
    contains
        !> calculate matrix elements
        procedure :: generate
end type

!> build a flattened helper to make things fast
type flathelper
    !> eigenvectors?
    complex(r8), dimension(:, :), allocatable :: ugv1
    complex(r8), dimension(:, :, :), allocatable :: ugv2, ugv3
    !> pre-transformed matrix element
    complex(r8), dimension(:), allocatable :: ptf_phi
    complex(r8), dimension(:), allocatable :: evp1
    complex(r8), dimension(:), allocatable :: evp2
end type

contains

!> Calculate all scattering rates
subroutine generate(sr, qpoint, ompoint, gpoint, qp, dr, uc, fc, fct, isoscatter, threephononscatter, skipsym, eigenvectors, closedgrid, tmr, mw, mem, verbosity)
    !> scattering rates
    class(lo_listofscatteringrates), intent(out) :: sr
    !> q-point in question
    type(lo_qpoint) :: qpoint
    !> harmonic properties at the relevant q-point
    type(lo_phonon_dispersions_qpoint), intent(in) :: ompoint
    !> harmonic properties at Gamma with the right direction
    type(lo_phonon_dispersions_qpoint), intent(in) :: gpoint
    !> qpoint mesh
    class(lo_qpoint_mesh), intent(in) :: qp
    !> phonon dispersions
    type(lo_phonon_dispersions), intent(in) :: dr
    !> crystal structure
    type(lo_crystalstructure), intent(in) :: uc
    !> second order forceconstant
    type(lo_forceconstant_secondorder), intent(inout) :: fc
    !> third order force constant
    type(lo_forceconstant_thirdorder), intent(in) :: fct
    !> calculate everything?
    logical, intent(in) :: isoscatter, threephononscatter, skipsym, eigenvectors, closedgrid
    !> timer
    type(lo_timer), intent(inout) :: tmr
    !> MPI communicator
    type(lo_mpi_helper), intent(inout) :: mw
    !> Memory tracker
    type(lo_mem_helper), intent(inout) :: mem
    !> how much to talk
    integer, intent(in) :: verbosity

    type(flathelper) :: fh
    real(r8) :: timer, t0, t1
    integer :: ind_gamma_full, ind_gamma_irr

    ! Start timer
    timer = walltime()
    t0 = timer
    t1 = timer

    call tmr%tick()

    init: block
        integer :: i

        ! Check if we are at gamma?
        if (norm2(qpoint%r) .lt. lo_sqtol) then
            sr%atgamma = .true.
        else
            sr%atgamma = .false.
        end if
        ! Skip using symmetry for paths, for example.
        if (skipsym) sr%atgamma = .false.

        ! Locate the index for Gamma on the q-grid, always useful
        ! to make sure we are consistent in everything, and that
        ! whenever I reference Gamma I pick a Gamma with the correct
        ! q-direction attached?
        ind_gamma_irr = -1
        ind_gamma_full = -1
        do i = 1, qp%n_irr_point
            if (lo_sqnorm(qp%ip(i)%r) .lt. lo_sqtol) then
                ! found Gamma
                ind_gamma_irr = i
                ind_gamma_full = qp%ip(i)%full_index
                exit
            end if
        end do
    end block init

    ! First get the third frequency at all q-points
    if (threephononscatter) then
    if (.not. sr%atgamma) then
        thirdfreq: block
            type(lo_phonon_dispersions_qpoint) :: op3
            complex(r8), dimension(:, :, :), allocatable :: egv3
            real(r8), dimension(:, :), allocatable :: omega3
            real(r8), dimension(:, :, :), allocatable :: vel3
            real(r8), dimension(3) :: qv1, qv2, qv3, v0
            integer, dimension(:,:,:), allocatable :: degenmode3
            integer, dimension(:,:), allocatable :: degeneracy3
            integer :: q, jq, j

            ! Make some space
            allocate (sr%qvec3(3, qp%n_full_point))
            sr%qvec3=0.0_r8

            call mem%allocate(egv3, [dr%n_mode, dr%n_mode, qp%n_full_point], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%allocate(omega3, [dr%n_mode, qp%n_full_point], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%allocate(vel3, [3,dr%n_mode, qp%n_full_point], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%allocate(degenmode3, [6,dr%n_mode, qp%n_full_point], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%allocate(degeneracy3, [dr%n_mode, qp%n_full_point], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            egv3=0.0_r8
            omega3=0.0_r8
            vel3=0.0_r8
            degenmode3=0
            degeneracy3=0

            ! sr%omega3 = 0.0_r8
            ! sr%vel3 = 0.0_r8
            ! egv3 = 0.0_r8
            ! qvec3 = 0.0_r8

            ! Calculate all the frequencies and stuff for the third q-point
            do q = 1, qp%n_full_point
                ! Make it parallel
                if (mod(q, mw%n) .ne. mw%r) cycle

                if (closedgrid) then
                    ! A closed grid, fetch data from the grid for the third point
                    qv1 = qpoint%r
                    qv2 = qp%ap(q)%r
                    qv3 = -qv1 - qv2
                    qv3 = matmul(uc%inv_reciprocal_latticevectors, qv3)
                    qv3 = lo_clean_fractional_coordinates(qv3)
                    jq = index_on_grid(qp, qv3)
                    omega3(:, q) = dr%aq(jq)%omega
                    vel3(:, :, q) = dr%aq(jq)%vel
                    egv3(:, :, q) = dr%aq(jq)%egv
                    sr%qvec3(:, q) = -qv1 - qv2
                    j=size(dr%aq(jq)%degenmode,1)
                    degenmode3(1:j,:,q) = dr%aq(jq)%degenmode(1:j,:)
                    degeneracy3(:,q) = dr%aq(jq)%degeneracy(:)
                else
                    ! Not a closed grid, just fetch everything.
                    ! Get the q-vectors
                    qv1 = qpoint%r
                    qv2 = qp%ap(q)%r
                    qv3 = -qv1 - qv2
                    ! Get harmonic things
                    call op3%generate(fc, uc, mem, qvec=qv3)
                    ! Store some harmonic things


                    v0 = matmul(uc%inv_reciprocal_latticevectors, qv3)
                    if (sum(abs(v0 - anint(v0))) .lt. lo_sqtol) then
                        ! we are at Gamma. Replace with proper Gamma.
                        sr%qvec3(:, q) = 0.0_r8
                        egv3(:, :, q) = gpoint%egv
                        vel3(:, :, q) = gpoint%vel
                        omega3(:, q)  = gpoint%omega
                        j=size(gpoint%degenmode,1)
                        degenmode3(1:j,:,q) = gpoint%degenmode(1:j,:)
                        degeneracy3(:,q)    = gpoint%degeneracy(:)
                    else
                        ! Not at gamma, continue normally
                        omega3(:, q) = op3%omega
                        vel3(:, :, q) = op3%vel
                        egv3(:, :, q) = op3%egv
                        sr%qvec3(:, q) = qv3
                        j=size(op3%degenmode,1)
                        degenmode3(1:j,:,q) = op3%degenmode(1:j,:)
                        degeneracy3(:,q)    = op3%degeneracy(:)
                    end if
                end if
            end do
            ! sync
            call mw%allreduce('sum', sr%qvec3)
            call mw%allreduce('sum', omega3)
            call mw%allreduce('sum', egv3)
            call mw%allreduce('sum', vel3)
            call mw%allreduce('sum', degeneracy3)
            call mw%allreduce('sum', degenmode3)

            ! Store as proper dispersions q-point object?
            allocate(sr%wp3(qp%n_full_point))
            do q=1,qp%n_full_point
                allocate(sr%wp3(q)%omega(uc%na*3))
                allocate(sr%wp3(q)%egv(uc%na*3,uc%na*3))
                allocate(sr%wp3(q)%vel(3,uc%na*3))
                allocate(sr%wp3(q)%degeneracy(uc%na*3))
                j=maxval(degeneracy3(:,q))
                allocate(sr%wp3(q)%degenmode(j,uc%na*3))
                sr%wp3(q)%omega = omega3(:,q)
                sr%wp3(q)%egv = egv3(:,:,q)
                sr%wp3(q)%vel = vel3(:,:,q)
                sr%wp3(q)%degeneracy = degeneracy3(:,q)
                sr%wp3(q)%degenmode = degenmode3(1:j,:,q)
            enddo

            ! Intermediate cleanup
            call mem%deallocate(omega3, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%deallocate(egv3, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%deallocate(vel3, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%deallocate(degeneracy3, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%deallocate(degenmode3, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)

            if (verbosity .gt. 0) then
                t1 = walltime()
                write (*, *) ''
                write (*, *) "... got q'' frequencies (", tochar(t1 - t0), "s)"
                t0 = t1
            end if
        end block thirdfreq
    end if
    end if

    call tmr%tock("q'' harmonic properties")

    ! Pre-scale things
    if (threephononscatter) then
        pretransform: block
            real(r8) :: f0, f1
            integer :: iq, imode, iatom, ix, ialpha, ctr

            ! Some temporary space for scaled eigenvectors
            if (sr%atgamma) then
                allocate (fh%ugv1(dr%n_mode, dr%n_mode))
                allocate (fh%ugv2(dr%n_mode, dr%n_mode, qp%n_irr_point))
                fh%ugv1 = 0.0_r8
                fh%ugv2 = 0.0_r8
            else
                allocate (fh%ugv1(dr%n_mode, dr%n_mode))
                allocate (fh%ugv2(dr%n_mode, dr%n_mode, qp%n_full_point))
                allocate (fh%ugv3(dr%n_mode, dr%n_mode, qp%n_full_point))
                fh%ugv1 = 0.0_r8
                fh%ugv2 = 0.0_r8
                fh%ugv3 = 0.0_r8
            end if

            ! First fix the one at the main q-point, always the same
            ! just pre-multiply with masses and frequencies
            do imode = 1, dr%n_mode
            do iatom = 1, uc%na
            do ix = 1, 3
                if (norm2(qpoint%r) .lt. lo_sqtol) then
                    ! Get consistent Gamma
                    if (gpoint%omega(imode) .gt. lo_freqtol) then
                        f0 = 1.0_r8/sqrt(gpoint%omega(imode))
                    else
                        f0 = 0.0_r8
                    end if
                    ialpha = (iatom - 1)*3 + ix
                    f1 = uc%invsqrtmass(iatom)
                    fh%ugv1(ialpha, imode) = gpoint%egv(ialpha, imode)*f0*f1
                else
                    ! Just normal copy
                    if (ompoint%omega(imode) .gt. lo_freqtol) then
                        f0 = 1.0_r8/sqrt(ompoint%omega(imode))
                    else
                        f0 = 0.0_r8
                    end if
                    ialpha = (iatom - 1)*3 + ix
                    f1 = uc%invsqrtmass(iatom)
                    fh%ugv1(ialpha, imode) = ompoint%egv(ialpha, imode)*f0*f1
                end if
            end do
            end do
            end do

            ! Pre-multiply eigenvectors with masses and frequencies
            if (sr%atgamma) then
                ! Case that will use symmetry eventually
                ctr = 0
                do iq = 1, qp%n_irr_point
                do imode = 1, dr%n_mode
                    ctr = ctr + 1
                    if (mod(ctr, mw%n) .ne. mw%r) cycle
                    if (iq .eq. ind_gamma_irr) then
                        if (gpoint%omega(imode) .gt. lo_freqtol) then
                            f0 = 1.0_r8/sqrt(gpoint%omega(imode))
                        else
                            f0 = 0.0_r8
                        end if

                        do iatom = 1, uc%na
                            f1 = uc%invsqrtmass(iatom)
                            do ix = 1, 3
                                ialpha = (iatom - 1)*3 + ix
                                fh%ugv2(ialpha, imode, iq) = gpoint%egv(ialpha, imode)*f0*f1
                            end do
                        end do
                    else
                        if (dr%iq(iq)%omega(imode) .gt. lo_freqtol) then
                            f0 = 1.0_r8/sqrt(dr%iq(iq)%omega(imode))
                        else
                            f0 = 0.0_r8
                        end if

                        do iatom = 1, uc%na
                            f1 = uc%invsqrtmass(iatom)
                            do ix = 1, 3
                                ialpha = (iatom - 1)*3 + ix
                                fh%ugv2(ialpha, imode, iq) = dr%iq(iq)%egv(ialpha, imode)*f0*f1
                            end do
                        end do
                    end if
                end do
                end do
                call mw%allreduce('sum', fh%ugv2)
            else
                ! Case that will not use symmetry
                ctr = 0
                do iq = 1, qp%n_full_point
                do imode = 1, dr%n_mode
                    ctr = ctr + 1
                    if (mod(ctr, mw%n) .ne. mw%r) cycle

                    ! Third guy
                    if (sr%wp3(iq)%omega(imode) .gt. lo_freqtol) then
                        f0 = 1.0_r8/sqrt( sr%wp3(iq)%omega(imode) )
                    else
                        f0 = 0.0_r8
                    end if
                    do iatom = 1, uc%na
                        f1 = uc%invsqrtmass(iatom)
                        do ix = 1, 3
                            ialpha = (iatom - 1)*3 + ix
                            fh%ugv3(ialpha, imode, iq) = sr%wp3(iq)%egv(ialpha,imode)*f0*f1
                        end do
                    end do

                    if (iq .eq. ind_gamma_full) then
                        ! pick consistent Gamma
                        if (gpoint%omega(imode) .gt. lo_freqtol) then
                            f0 = 1.0_r8/sqrt(gpoint%omega(imode))
                        else
                            f0 = 0.0_r8
                        end if
                        do iatom = 1, uc%na
                            f1 = uc%invsqrtmass(iatom)
                            do ix = 1, 3
                                ialpha = (iatom - 1)*3 + ix
                                fh%ugv2(ialpha, imode, iq) = gpoint%egv(ialpha, imode)*f0*f1
                            end do
                        end do
                    else
                        ! Not at Gamma, don't have to care
                        if (dr%aq(iq)%omega(imode) .gt. lo_freqtol) then
                            f0 = 1.0_r8/sqrt(dr%aq(iq)%omega(imode))
                        else
                            f0 = 0.0_r8
                        end if
                        do iatom = 1, uc%na
                            f1 = uc%invsqrtmass(iatom)
                            do ix = 1, 3
                                ialpha = (iatom - 1)*3 + ix
                                fh%ugv2(ialpha, imode, iq) = dr%aq(iq)%egv(ialpha, imode)*f0*f1
                            end do
                        end do
                    end if

                end do
                end do
                call mw%allreduce('sum', fh%ugv2)
                call mw%allreduce('sum', fh%ugv3)
            end if

            ! Space for the pre-transformed phi
            allocate (fh%ptf_phi(dr%n_mode**3))
            fh%ptf_phi = 0.0_r8

            ! Space for outer-producted eigenvectors?
            allocate (fh%evp1(dr%n_mode**2))
            allocate (fh%evp2(dr%n_mode**3))
            fh%evp1 = 0.0_r8
            fh%evp2 = 0.0_r8

            if (verbosity .gt. 0) then
                t1 = walltime()
                write (*, *) "... pretransformed and made temporary space (", tochar(t1 - t0), "s)"
                t0 = t1
            end if
        end block pretransform
    end if

    call tmr%tock("eigenvector scaling")


    ! Then the isotope scattering rates
    if (isoscatter) then
        isosc: block
            integer :: q, b1, b2, l

            if (verbosity .gt. 0) call lo_progressbar_init()

            if (sr%atgamma) then
                allocate (sr%psi_iso(dr%n_mode, dr%n_mode, qp%n_irr_point))
                sr%psi_iso = 0.0_r8
                l = 0
                do q = 1, qp%n_irr_point
                    do b1 = 1, dr%n_mode
                    do b2 = 1, dr%n_mode
                        l = l + 1
                        ! this is how I distribute over MPI
                        if (mod(l, mw%n) .ne. mw%r) cycle
                        sr%psi_iso(b1, b2, q) = isotope_scattering_strength(uc, ompoint%egv(:, b1), dr%iq(q)%egv(:, b2), ompoint%omega(b1), dr%iq(q)%omega(b2))
                    end do
                    end do
                    if (verbosity .gt. 0) then
                        if (lo_trueNtimes(q, 127, qp%n_full_point)) call lo_progressbar('... isotope matrixelements', q, qp%n_full_point)
                    end if
                end do
            else
                allocate (sr%psi_iso(dr%n_mode, dr%n_mode, qp%n_full_point))
                sr%psi_iso = 0.0_r8
                l = 0
                do q = 1, qp%n_full_point
                    do b1 = 1, dr%n_mode
                    do b2 = 1, dr%n_mode
                        l = l + 1
                        ! this is how I distribute over MPI
                        if (mod(l, mw%n) .ne. mw%r) cycle
                        sr%psi_iso(b1, b2, q) = isotope_scattering_strength(uc, ompoint%egv(:, b1), dr%aq(q)%egv(:, b2), ompoint%omega(b1), dr%aq(q)%omega(b2))
                    end do
                    end do
                    if (verbosity .gt. 0) then
                        if (lo_trueNtimes(q, 127, qp%n_full_point)) call lo_progressbar('... isotope matrixelements', q, qp%n_full_point)
                    end if
                end do
            end if

            ! sum it up
            call mw%allreduce('sum', sr%psi_iso)

            if (verbosity .gt. 0) call lo_progressbar('... isotope matrixelements', qp%n_full_point, qp%n_full_point, walltime() - t0)
        end block isosc
    end if

    call tmr%tock("isotope matrix elements")

    ! get the threephonon matrix elements
    if (threephononscatter) then
        threephsc2: block
            real(r8) :: t_up, t_tot
            real(r8), dimension(3) :: qv2, qv3
            complex(r8), dimension(:, :, :), allocatable :: psi_3ph_tmp
            integer :: q, b1, b2, b3, ctr, l

            ! Some space
            t0 = walltime()
            if (verbosity .gt. 0) call lo_progressbar_init()

            if (sr%atgamma) then
                ! At Gamma, we can use symmetry and only calculate a couple of matrix elements.
                ! Space for matrix elements
                allocate (sr%psi_3ph(dr%n_mode, dr%n_mode, dr%n_mode, qp%n_irr_point))
                sr%psi_3ph = 0.0_r8

                ! temporary array for scattering rates
                allocate (psi_3ph_tmp(dr%n_mode, dr%n_mode, dr%n_mode))
                psi_3ph_tmp = 0.0_r8

                ctr = 0
                l = 0
                do q = 1, qp%n_irr_point

                    ! fkdev check for l overflow
                    if (l < 0) then
                        write (*, *) 'FIXME MPI: l overflow in lineshape/scatteringrates.f90'
                        stop
                    end if

                    ! fetch q-vectors
                    qv2 = qp%ip(q)%r
                    qv3 = -qv2
                    ! do the half-transform
                    call pretransform_phi(fct, qv2, qv3, fh%ptf_phi)
                    !call fct%pretransform(qv2,qv3,fh%ptf_phi)

                    ! get matrix elements for all modes
                    psi_3ph_tmp = 0.0_r8
                    do b1 = 1, dr%n_mode

                        do b2 = 1, dr%n_mode
                            ! MPI division
                            l = l + 1
                            if (mod(l, mw%n) .ne. mw%r) cycle
                            fh%evp1 = 0.0_r8
                            call zgeru(dr%n_mode, dr%n_mode, (1.0_r8, 0.0_r8), fh%ugv2(:, b2, q), 1, fh%ugv1(:, b1), 1, fh%evp1, dr%n_mode)

                            do b3 = 1, dr%n_mode
                                fh%evp2 = 0.0_r8
                                ! complicated conjugation thingy compared with the unsymmetric below.
                                call zgerc(dr%n_mode, dr%n_mode*dr%n_mode, (1.0_r8, 0.0_r8), fh%ugv2(:, b3, q), 1, fh%evp1, 1, fh%evp2, dr%n_mode)
                                psi_3ph_tmp(b1, b2, b3) = dot_product(fh%evp2, fh%ptf_phi)
                            end do
                        end do

                        if (verbosity .gt. 0) then
                            ctr = ctr + 1
                            if (lo_trueNtimes(ctr, 127, qp%n_irr_point*dr%n_mode)) then
                                call lo_progressbar('... threephonon matrixelements', ctr, qp%n_irr_point*dr%n_mode, walltime() - t0)
                            end if
                        end if
                    end do

                    ! fkdev: MPI barrier here to avoid parallelizing over q
                    ! fixed a bug for certain q-meshes on Dardel
                    call mw%allreduce('sum', psi_3ph_tmp)
                    sr%psi_3ph(:, :, :, q) = psi_3ph_tmp

                    ! fkdev: report on scattering times
                    if ((mw%talk) .and. (verbosity .gt. 0)) then
                        t_up = walltime() - t0
                        t_tot = t_up/q*qp%n_irr_point
                        write (*, '(A,I7,A,I7,A,F10.1,A,F10.1,A)') &
                            ' ... 3ph scattering: ', q, ' of ', qp%n_irr_point, ' qpoints done, time: ', t_up, 's , projected total time: ', t_tot, 's'
                    end if
                end do
                deallocate (psi_3ph_tmp)
            else
                ! Not at Gamma, have to calculate all of them
                ! Space for matrix elements
                allocate (sr%psi_3ph(dr%n_mode, dr%n_mode, dr%n_mode, qp%n_full_point))
                sr%psi_3ph = 0.0_r8

                ctr = 0
                l = 0
                do q = 1, qp%n_full_point
                    ! fetch q-vectors
                    qv2 = qp%ap(q)%r
                    qv3 = sr%qvec3(:, q)
                    ! do the half-transform
                    call pretransform_phi(fct, qv2, qv3, fh%ptf_phi)
                    !call fct%pretransform(qv2,qv3,fh%ptf_phi)
                    ! get matrix elements for all modes
                    do b1 = 1, dr%n_mode
                        do b2 = 1, dr%n_mode
                            ! MPI division
                            l = l + 1
                            if (mod(l, mw%n) .ne. mw%r) cycle
                            fh%evp1 = 0.0_r8
                            call zgeru(dr%n_mode, dr%n_mode, (1.0_r8, 0.0_r8), fh%ugv2(:, b2, q), 1, fh%ugv1(:, b1), 1, fh%evp1, dr%n_mode)
                            do b3 = 1, dr%n_mode
                                fh%evp2 = 0.0_r8
                                call zgeru(dr%n_mode, dr%n_mode*dr%n_mode, (1.0_r8, 0.0_r8), fh%ugv3(:, b3, q), 1, fh%evp1, 1, fh%evp2, dr%n_mode)
                                fh%evp2 = conjg(fh%evp2)
                                sr%psi_3ph(b1, b2, b3, q) = dot_product(fh%evp2, fh%ptf_phi)
                            end do
                        end do

                        if (verbosity .gt. 0) then
                            ctr = ctr + 1
                            if (lo_trueNtimes(ctr, 127, qp%n_full_point*dr%n_mode)) then
                                call lo_progressbar('... threephonon matrixelements', ctr, qp%n_full_point*dr%n_mode, walltime() - t0)
                            end if
                        end if
                    end do
                end do
                ! sync across ranks
                call mw%allreduce('sum', sr%psi_3ph)

            end if

            ! Cleanup
            if ( allocated(fh%ugv1   ) ) deallocate(fh%ugv1   )
            if ( allocated(fh%ugv2   ) ) deallocate(fh%ugv2   )
            if ( allocated(fh%ugv3   ) ) deallocate(fh%ugv3   )
            if ( allocated(fh%ptf_phi) ) deallocate(fh%ptf_phi)
            if ( allocated(fh%evp1   ) ) deallocate(fh%evp1   )
            if ( allocated(fh%evp2   ) ) deallocate(fh%evp2   )

            if (verbosity .gt. 0) then
                t1 = walltime()
                call lo_progressbar('... threephonon matrixelements', qp%n_full_point*dr%n_mode, qp%n_full_point*dr%n_mode, t1 - t0)
                t0 = t1
            end if
        end block threephsc2

        call tmr%tock("three-phonon matrix elements")
    end if

    call mem%tock(__FILE__, __LINE__, mw%comm)

    ! if ( mw%talk ) write(*,*) 'done here for now ',__FILE__,__LINE__
    ! call mw%destroy()
    ! stop
end subroutine

!> isotope scattering strength
pure function isotope_scattering_strength(p, egv1, egv2, om1, om2) result(lambda)
    !> crystal structure
    type(lo_crystalstructure), intent(in) :: p
    !> first frequency
    real(r8), intent(in) :: om1
    !> second frequency
    real(r8), intent(in) :: om2
    !> first eigenvector
    complex(r8), dimension(:), intent(in) :: egv1
    !> second eigenvector
    complex(r8), dimension(:), intent(in) :: egv2
    !> scattering strength
    real(r8) :: lambda

    integer :: i, j
    real(r8) :: f0, f1
    complex(r8), dimension(3) :: cv0, cv1

    f1 = 0.0_r8
    do i = 1, p%na
        cv0 = egv1((i - 1)*3 + 1:(i*3))
        cv1 = egv2((i - 1)*3 + 1:(i*3))
        f0 = 0.0_r8
        do j = 1, 3
            f0 = f0 + abs(conjg(cv0(j))*cv1(j))
        end do
        f0 = f0**2
        f1 = f1 + f0*p%isotope(i)%disorderparameter
    end do
    lambda = f1*om1*om2
end function

!> pre-transform to get half of the matrix element
subroutine pretransform_phi(fct, q2, q3, ptf)
    !> third order forceconstant
    type(lo_forceconstant_thirdorder), intent(in) :: fct
    !> q-vectors
    real(r8), dimension(3), intent(in) :: q2, q3
    !> flattened, pretransformed matrix element
    complex(r8), dimension(:), intent(out) :: ptf

    integer :: i, j, k, l

    complex(r8) :: expiqr
    real(r8), dimension(3) :: rv2, rv3
    real(r8) :: iqr
    integer :: a1, a2, a3, ia, ib, ic, t, nb

    nb = fct%na*3
    ptf = 0.0_r8
    do a1 = 1, fct%na
    do t = 1, fct%atom(a1)%n
        a2 = fct%atom(a1)%triplet(t)%i2
        a3 = fct%atom(a1)%triplet(t)%i3

        rv2 = fct%atom(a1)%triplet(t)%lv2
        rv3 = fct%atom(a1)%triplet(t)%lv3

        iqr = dot_product(q2, rv2) + dot_product(q3, rv3)
        iqr = -iqr*lo_twopi
        expiqr = cmplx(cos(iqr), sin(iqr), r8)
        do i = 1, 3
        do j = 1, 3
        do k = 1, 3
            ia = (a1 - 1)*3 + i
            ib = (a2 - 1)*3 + j
            ic = (a3 - 1)*3 + k
            ! Now for the grand flattening scheme, consistent with the zgeru operations above.
            l = (ia - 1)*nb*nb + (ib - 1)*nb + ic
            ptf(l) = ptf(l) + fct%atom(a1)%triplet(t)%m(i, j, k)*expiqr
        end do
        end do
        end do
    end do
    end do
end subroutine

!> derivative of half-transformed three-phonon matrix element

! !> returns the index on the grid that gives q3=-q1-q2
! pure function fft_third_grid_index(i1,i2,dims) result(i3)
!     !> index to q1
!     integer, intent(in) :: i1
!     !> index to q2
!     integer, intent(in) :: i2
!     !> dimensions of the grid
!     integer, dimension(3), intent(in) :: dims
!     !> index to q3
!     integer :: i3
!
!     integer, dimension(3) :: gi1,gi2,gi3
!     integer :: l,k
!
!     ! Convert triplet to singlet
!     gi1=singlet_to_triplet(i1,dims(2),dims(3))
!     gi2=singlet_to_triplet(i2,dims(2),dims(3))
!     do l=1,3
!         gi3(l)=3-gi1(l)-gi2(l)
!     enddo
!     do k=1,3
!     do l=1,3
!         if ( gi3(l) .lt. 1 ) gi3(l)=gi3(l)+dims(l)
!         if ( gi3(l) .gt. dims(l) ) gi3(l)=gi3(l)-dims(l)
!     enddo
!     enddo
!     ! convert it back to a singlet
!     i3=triplet_to_singlet(gi3,dims(2),dims(3))
!
!     contains
!     !> convert a linear index to a triplet
!     pure function singlet_to_triplet(l,ny,nz) result(gi)
!         !> linear index
!         integer, intent(in) :: l
!         !> second dimension
!         integer, intent(in) :: ny
!         !> third dimension
!         integer, intent(in) :: nz
!         !> grid-index
!         integer, dimension(3) :: gi
!
!         integer :: i,j,k
!
!         k=mod(l,nz)
!         if ( k .eq. 0 ) k=nz
!         j=mod( (l-k)/nz , ny )+1
!         i=(l-k-(j-1)*nz)/(nz*ny)+1
!         gi=[i,j,k]
!     end function
!     !> convert a triplet index to a singlet
!     pure function triplet_to_singlet(gi,ny,nz) result(l)
!         !> grid-index
!         integer, dimension(3), intent(in) :: gi
!         !> second dimension
!         integer, intent(in) :: ny
!         !> third dimension
!         integer, intent(in) :: nz
!         !> linear index
!         integer :: l
!
!         l=(gi(1)-1)*ny*nz+(gi(2)-1)*nz+gi(3)
!     end function
! end function

end module
