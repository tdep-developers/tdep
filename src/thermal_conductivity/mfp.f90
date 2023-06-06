#include "precompilerdefinitions"
module mfp
!!
!! Generates mean free path plots and other useful things once the calculations are finished.
!!
use konstanter, only: r8, lo_huge, lo_hugeint, lo_pi, lo_twopi, lo_tol, lo_sqtol, lo_status, lo_exitcode_param, lo_freqtol, &
                      lo_kappa_SI_to_au, lo_groupvel_ms_to_Hartreebohr, lo_bohr_to_m, lo_kappa_au_to_SI, lo_kb_Hartree, &
                      lo_frequency_hartree_to_icm, lo_frequency_hartree_to_meV, lo_frequency_hartree_to_THz
use gottochblandat, only: lo_planck, lo_gauss, lo_trapezoid_integration, lo_outerproduct, walltime, tochar, lo_chop, &
                          lo_harmonic_oscillator_cv, lo_linear_interpolation, lo_linspace, lo_return_unique, &
                          lo_logspace, qsort, lo_trace
use mpi_wrappers, only: lo_mpi_helper, lo_stop_gracefully, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_IN_PLACE
use lo_memtracker, only: lo_mem_helper
use type_crystalstructure, only: lo_crystalstructure
use type_qpointmesh, only: lo_qpoint_mesh, lo_LV_tetrahedron_weights, lo_integration_weights_for_one_tetrahedron
use type_phonon_dispersions, only: lo_phonon_dispersions
use type_phonon_dos, only: lo_phonon_dos
use type_symmetryoperation, only: lo_operate_on_vector
use hdf5_wrappers, only: lo_hdf5_helper, lo_h5_store_data, lo_h5_store_attribute, HID_T, H5F_ACC_TRUNC_F, &
                         h5close_f, h5open_f, h5fclose_f, h5fopen_f, h5fcreate_f, h5gclose_f, h5gopen_f, h5gcreate_f
implicit none
private
public :: lo_mfp, get_cumulative_plots, write_cumulative_plots

!> the cumulative kappa plot for one temperature
type lo_mfp_temperature
    !> the temperature
    real(r8) :: temperature = -lo_huge
    !> total kappa
    real(r8), dimension(9) :: kappa_total = -lo_huge
    !> mean free path
    real(r8), dimension(:), allocatable :: mfaxis
    real(r8), dimension(:), allocatable :: mf_kappa
    real(r8), dimension(:, :), allocatable :: mf_kappa_band
    real(r8), dimension(:, :), allocatable :: mf_kappa_atom
    !> frequency
    real(r8), dimension(:, :), allocatable :: fq_kappa
    real(r8), dimension(:, :, :), allocatable :: fq_kappa_band
    real(r8), dimension(:, :, :), allocatable :: fq_kappa_atom
    !> boundary scattering
    real(r8), dimension(:), allocatable :: boundary_xaxis
    real(r8), dimension(:, :, :), allocatable :: boundary_kappa
    !> angular momentum matrix
    real(r8), dimension(3, 3) :: angmomalpha
    !> spectral angular momentum matrix
    real(r8), dimension(:, :), allocatable :: fq_angmom
    real(r8), dimension(:, :, :), allocatable :: fq_angmom_band
    real(r8), dimension(:, :, :), allocatable :: fq_angmom_atom
end type

!> thin film kappa
type lo_mfp_thinfilm
    !> temperature
    real(r8) :: temperature = -lo_huge
    !> the film normal
    real(r8), dimension(3) :: filmnormal = -lo_huge
    !> the temperature gradient
    real(r8), dimension(3) :: tempgradient = -lo_huge
    !> the x-axis
    real(r8), dimension(:), allocatable :: thicknessaxis
    !> the total kappa
    real(r8), dimension(:), allocatable :: kc_tot
    !> band decomposed cumulative kappa
    real(r8), dimension(:, :), allocatable :: kc_band
end type

!> place to hold mean free paths
type lo_mfp
    !> how many points on the x-axis
    integer :: np = -lo_hugeint
    !> for how many temperatures is this done?
    integer :: nt = -lo_hugeint
    !> the cumulative kappa plot, one for each temperature
    type(lo_mfp_temperature), dimension(:), allocatable :: temp
    !> the thin film kappa plot
    type(lo_mfp_thinfilm), dimension(:), allocatable :: film
end type

contains

!> get the mfp histogram thing for one temperature
subroutine get_cumulative_plots(mf, qp, dr, pd, uc, np, temperature, sigma, kappa, mw, mem)
    !> the cumulative plot
    type(lo_mfp_temperature), intent(inout) :: mf
    !> the q-grid
    class(lo_qpoint_mesh), intent(in) :: qp
    !> the dispersions and kappa
    type(lo_phonon_dispersions), intent(in) :: dr
    !> phonon density of states
    type(lo_phonon_dos), intent(inout) :: pd
    !> structure
    type(lo_crystalstructure), intent(in) :: uc
    !> the number of points on the x-axis
    integer, intent(in) :: np
    !> the current temperature
    real(r8), intent(in) :: temperature
    !> additional smearing for mean free path plots?
    real(r8), intent(in) :: sigma
    !> total kappa
    real(r8), dimension(3, 3) :: kappa
    !> MPI helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory helper
    type(lo_mem_helper), intent(inout) :: mem

    real(r8), parameter :: kappatol = 1E-5_r8*lo_kappa_SI_to_au
    real(r8), dimension(:, :, :), allocatable :: kappa_chop, kappa_flat
    real(r8) :: totkappa

    ! Chop up kappa into irreducible contributions
    chopkappa: block
        complex(r8), dimension(3) :: cv0
        real(r8), dimension(:, :, :), allocatable :: kf
        real(r8), dimension(uc%na) :: prj
        real(r8) :: f0, f1
        integer :: q, atom, mode, i, ii

        ! Make an extremely flat thing. Start by adding everything together into the irreducible wedge.
        lo_allocate(kappa_flat(qp%n_full_point, dr%n_mode, 9))
        lo_allocate(kf(9, qp%n_irr_point, dr%n_mode))
        kappa_flat = 0.0_r8
        kf = 0.0_r8

        do q = 1, qp%n_full_point
            if (mod(q, mw%n) .ne. mw%r) cycle
            ii = qp%ap(q)%irreducible_index
            do mode = 1, dr%n_mode
                kf(1, ii, mode) = kf(1, ii, mode) + dr%aq(q)%kappa(1, 1, mode)
                kf(2, ii, mode) = kf(2, ii, mode) + dr%aq(q)%kappa(1, 2, mode)
                kf(3, ii, mode) = kf(3, ii, mode) + dr%aq(q)%kappa(1, 3, mode)
                kf(4, ii, mode) = kf(4, ii, mode) + dr%aq(q)%kappa(2, 1, mode)
                kf(5, ii, mode) = kf(5, ii, mode) + dr%aq(q)%kappa(2, 2, mode)
                kf(6, ii, mode) = kf(6, ii, mode) + dr%aq(q)%kappa(2, 3, mode)
                kf(7, ii, mode) = kf(7, ii, mode) + dr%aq(q)%kappa(3, 1, mode)
                kf(8, ii, mode) = kf(8, ii, mode) + dr%aq(q)%kappa(3, 2, mode)
                kf(9, ii, mode) = kf(9, ii, mode) + dr%aq(q)%kappa(3, 3, mode)
            end do
        end do
        call mw%allreduce('sum', kf)
        !call mpi_allreduce(MPI_IN_PLACE,kf,qp%n_irr_point*dr%n_mode*9,MPI_DOUBLE_PRECISION,MPI_SUM,mw%comm,mw%error)
        ! Now rotate this out again
        do q = 1, qp%n_full_point
            ii = qp%ap(q)%irreducible_index
            f0 = 1.0_r8/qp%ip(ii)%integration_weight
            do mode = 1, dr%n_mode
                kappa_flat(q, mode, 1) = kf(1, ii, mode)*f0
                kappa_flat(q, mode, 2) = kf(2, ii, mode)*f0
                kappa_flat(q, mode, 3) = kf(3, ii, mode)*f0
                kappa_flat(q, mode, 4) = kf(4, ii, mode)*f0
                kappa_flat(q, mode, 5) = kf(5, ii, mode)*f0
                kappa_flat(q, mode, 6) = kf(6, ii, mode)*f0
                kappa_flat(q, mode, 7) = kf(7, ii, mode)*f0
                kappa_flat(q, mode, 8) = kf(8, ii, mode)*f0
                kappa_flat(q, mode, 9) = kf(9, ii, mode)*f0
            end do
        end do
        kappa_flat = lo_chop(kappa_flat/qp%n_full_point, kappatol)
        ! And grab the total Kappa
        do i = 1, 9
            mf%kappa_total(i) = sum(kappa_flat(:, :, i))/qp%n_full_point
        end do

        lo_allocate(kappa_chop(qp%n_irr_point, dr%n_mode, uc%na))
        kappa_chop = 0.0_r8
        do q = 1, qp%n_irr_point
            if (mod(q, mw%n) .ne. mw%r) cycle
            do mode = 1, dr%n_mode
                ! get the projection per atom
                prj = 0.0_r8
                do atom = 1, uc%na
                    cv0 = dr%iq(q)%egv((atom - 1)*3 + 1:atom*3, mode)
                    prj(atom) = prj(atom) + abs(dot_product(conjg(cv0), cv0))
                end do
                prj = lo_chop(prj/sum(prj), lo_sqtol)
                ! the scalar thermal conductivity
                f0 = norm2(kf(:, q, mode))/sqrt(3.0_r8)
                do atom = 1, uc%na
                    kappa_chop(q, mode, atom) = f0*prj(atom)
                end do
            end do
        end do
        call mw%allreduce('sum', kappa_chop)
        !call mpi_allreduce(MPI_IN_PLACE,kappa_chop,qp%n_irr_point*dr%n_mode*uc%na,MPI_DOUBLE_PRECISION,MPI_SUM,mw%comm,mw%error)

        ! Make sure it adds up to the total
        totkappa = sum(kappa)/3
        f1 = sum(kappa_chop)
        kappa_chop = kappa_chop*totkappa/f1
        lo_deallocate(kf)

        ! Also, if it's too small, take care of that
        if (totkappa*lo_kappa_au_to_SI .lt. 10*lo_tol) then
            ! Just put zeros everywhere
            mf%temperature = temperature
            lo_allocate(mf%mfaxis(np))
            lo_allocate(mf%mf_kappa(np))
            lo_allocate(mf%mf_kappa_band(np, dr%n_mode))
            lo_allocate(mf%mf_kappa_atom(np, uc%na))
            mf%mfaxis = 0.0_r8
            mf%mf_kappa = 0.0_r8
            mf%mf_kappa_band = 0.0_r8
            mf%mf_kappa_atom = 0.0_r8
            ! put zeros for the frequency part as well
            lo_allocate(mf%fq_kappa(np, 9))
            lo_allocate(mf%fq_kappa_band(np, dr%n_mode, 9))
            lo_allocate(mf%fq_kappa_atom(np, uc%na, 9))
            mf%fq_kappa = 0.0_r8
            mf%fq_kappa_band = 0.0_r8
            mf%fq_kappa_atom = 0.0_r8
            ! and the boundary scattering
            lo_allocate(mf%boundary_xaxis(np))
            lo_allocate(mf%boundary_kappa(3, 3, np))
            mf%boundary_xaxis = 0.0_r8
            mf%boundary_kappa = 0.0_r8
            return
        end if
    end block chopkappa

    ! First fix cumulative vs mean free path
    cmf: block
        real(r8), parameter :: mfptol = 1E-12_r8*lo_groupvel_ms_to_Hartreebohr
        integer, parameter :: refinefactor = 20
        real(r8), dimension(:, :, :), allocatable :: ym, yn
        real(r8), dimension(:, :, :), allocatable :: xm, xn
        real(r8), dimension(:, :), allocatable :: sitebuf
        real(r8), dimension(:), allocatable :: x, y

        real(r8) :: minx, maxx, sig1, sig2, isig, f1, prefactor
        real(r8) :: min_mfp, max_mfp
        real(r8) :: a, b, logsig
        integer, dimension(:), allocatable :: ind
        integer :: npts, ii, jj, mode, atom, ctr, q, i, j, k
        real(r8) :: tt0, tt1, tt2, f0
        !
        tt0 = walltime()
        tt1 = 0.0_r8
        tt2 = 0.0_r8
        npts = np*refinefactor ! tighter grid to interpolate later

        ! range of mean free paths, to determine the smearing parameter
        lo_allocate(x(qp%n_irr_point*dr%n_mode))
        min_mfp = lo_huge
        max_mfp = 0.0_r8
        i = 0
        do q = 1, qp%n_irr_point
            min_mfp = min(min_mfp, minval(dr%iq(q)%scalar_mfp, dr%iq(q)%scalar_mfp > mfptol))
            max_mfp = max(max_mfp, maxval(dr%iq(q)%scalar_mfp))
            do j = 1, dr%n_mode
                i = i + 1
                if (dr%iq(q)%scalar_mfp(j) > mfptol) then
                    x(i) = log(dr%iq(q)%scalar_mfp(j))
                else
                    x(i) = log(mfptol)
                end if
            end do
        end do
        call qsort(x)
        a = x(ceiling(size(x)*0.5_r8))
        ii = min(size(x), ceiling(size(x)*0.98_r8))
        b = x(ii)
        logsig = 0.0_r8
        j = 0
        do i = 1, size(x) - 1
            if (x(i) .gt. a .and. x(i) .lt. b) then
                logsig = max(logsig, x(i + 1) - x(i))
                j = j + 1
            end if
        end do
        logsig = logsig*sigma
        lo_deallocate(x)

        ! temporary buffers
        lo_allocate(ym(npts, dr%n_mode, uc%na))
        lo_allocate(yn(npts, dr%n_mode, uc%na))
        lo_allocate(xm(npts, dr%n_mode, uc%na))
        lo_allocate(xn(npts, dr%n_mode, uc%na))
        lo_allocate(x(qp%n_irr_point))
        lo_allocate(y(qp%n_irr_point))
        lo_allocate(ind(qp%n_irr_point))
        xm = 0.0_r8
        ym = 0.0_r8
        xn = 0.0_r8
        yn = 0.0_r8
        x = 0.0_r8
        y = 0.0_r8
        ind = 0

        ctr = 0
        do atom = 1, uc%na
        do mode = 1, dr%n_mode
            ctr = ctr + 1
            ! Do the fancy MPI-division here!
            if (mod(ctr, mw%n) .ne. mw%r) cycle

            ! Fetch mfp and kappa, convert mfp to nm, and sort the values according to mean free path
            y = 0.0_r8
            do q = 1, qp%n_irr_point
                x(q) = dr%iq(q)%scalar_mfp(mode)
                y(q) = kappa_chop(q, mode, atom)
            end do
            call qsort(x, ind)
            y = y(ind)
            minx = minval(x, x > lo_tol)*0.1_r8
            maxx = maxval(x)*3.0_r8

            ! Calculate the actual mean free path thing, the integrated quantity
            call lo_linspace(minx, maxx, xm(:, mode, atom))

            ! Actually sum things up
            do q = 1, qp%n_irr_point
                if (x(q) .lt. minx) cycle
                if (x(q) .gt. maxx) cycle

                ! Weird log-gaussian thing
                ! exp( -( log(x)-log(mu) )^2/b )
                ! log(x)-log(mu)=-4*b
                a = log(x(q))  !< mu, sort of
                b = logsig
                sig1 = exp(a - 4*b) - x(q) ! for integration ranges
                sig2 = exp(a + 4*b) - x(q) ! for integration ranges
                prefactor = sqrt(b*lo_pi)*exp(a + b*0.25_r8)
                isig = 1.0_r8/b

                ii = max(1, floor(npts*(x(q) + sig1 - minx)/maxx))
                jj = min(npts, ceiling(npts*(x(q) + sig2 - minx)/maxx))
                if (jj .le. ii) cycle
                do k = ii, jj
                    f1 = (1.0_r8 + erf((log(xm(k, mode, atom)) - log(x(q)))*isig))*0.5_r8
                    ym(k, mode, atom) = ym(k, mode, atom) + f1*y(q)
                end do
                do k = jj + 1, npts
                    ym(k, mode, atom) = ym(k, mode, atom) + y(q)
                end do
            end do
        end do
        end do
        ! Sum it over MPI
        call mw%allreduce('sum', xm, xn)
        call mw%allreduce('sum', ym, yn)
        !call mpi_allreduce(xm,xn,npts*dr%n_mode*uc%na,MPI_DOUBLE_PRECISION,MPI_SUM,mw%comm,mw%error)
        !call mpi_allreduce(ym,yn,npts*dr%n_mode*uc%na,MPI_DOUBLE_PRECISION,MPI_SUM,mw%comm,mw%error)

        ! Now add it together properly
        mf%temperature = temperature
        lo_allocate(mf%mfaxis(np))
        minx = minval(xn)
        maxx = maxval(xn)*2
        call logspace(minx, maxx, mf%mfaxis, smallest_separation=minx*10)
        allocate (mf%mf_kappa(np))
        allocate (mf%mf_kappa_band(np, dr%n_mode))
        allocate (mf%mf_kappa_atom(np, uc%na))
        mf%mf_kappa = 0.0_r8
        mf%mf_kappa_band = 0.0_r8
        mf%mf_kappa_atom = 0.0_r8

        ! Sum up things
        do mode = 1, dr%n_mode
        do atom = 1, uc%na
        do i = 1, np
            if (mf%mfaxis(i) .lt. xn(1, mode, atom)*1.01_r8) then
                f0 = 0.0_r8
            elseif (mf%mfaxis(i) .gt. xn(npts, mode, atom)*0.999_r8) then
                f0 = yn(npts, mode, atom)
            else
                f0 = lo_linear_interpolation(xn(:, mode, atom), yn(:, mode, atom), mf%mfaxis(i))
            end if
            mf%mf_kappa(i) = mf%mf_kappa(i) + f0
            mf%mf_kappa_band(i, mode) = mf%mf_kappa_band(i, mode) + f0
            mf%mf_kappa_atom(i, atom) = mf%mf_kappa_atom(i, atom) + f0
        end do
        end do
        end do

        ! Make sure it's monotonically increasing
        do i = 1, np - 1
            if (mf%mf_kappa(i + 1) .lt. mf%mf_kappa(i)) mf%mf_kappa(i + 1) = mf%mf_kappa(i)
            do atom = 1, uc%na
                if (mf%mf_kappa_atom(i + 1, atom) .lt. mf%mf_kappa_atom(i, atom)) then
                    mf%mf_kappa_atom(i + 1, atom) = mf%mf_kappa_atom(i, atom)
                end if
            end do
            do mode = 1, dr%n_mode
                if (mf%mf_kappa_band(i + 1, mode) .lt. mf%mf_kappa_band(i, mode)) then
                    mf%mf_kappa_band(i + 1, mode) = mf%mf_kappa_band(i, mode)
                end if
            end do
        end do

        ! Fix site degeneracy
        lo_allocate(sitebuf(np, uc%na))
        sitebuf = 0.0_r8
        do i = 1, uc%na
            ! enfore site degeneracy
            do j = 1, uc%sym%degeneracy(i)
                k = uc%sym%degenerate_atom(j, i)
                sitebuf(:, i) = sitebuf(:, i) + mf%mf_kappa_atom(:, k)/(1.0_r8*uc%sym%degeneracy(j))
            end do
        end do
        mf%mf_kappa_atom = sitebuf
        lo_deallocate(sitebuf)

        ! Now make sure everything is normalized properly.
        mf%mf_kappa = mf%mf_kappa*totkappa/mf%mf_kappa(np)
        do i = 1, np
            ! fix the band-decomposed
            f0 = sum(mf%mf_kappa_band(i, :))
            f1 = mf%mf_kappa(i)
            if (f1 .gt. kappatol .and. f0 .gt. kappatol) then
                mf%mf_kappa_band(i, :) = mf%mf_kappa_band(i, :)*mf%mf_kappa(i)/f0
            else
                mf%mf_kappa_band(i, :) = 0.0_r8
            end if
            ! fix the atom-decomposed
            f0 = sum(mf%mf_kappa_atom(i, :))
            f1 = mf%mf_kappa(i)
            if (f1 .gt. kappatol .and. f0 .gt. kappatol) then
                mf%mf_kappa_atom(i, :) = mf%mf_kappa_atom(i, :)*mf%mf_kappa(i)/f0
            else
                mf%mf_kappa_atom(i, :) = 0.0_r8
            end if
        end do
        lo_deallocate(ind)
        lo_deallocate(xm)
        lo_deallocate(ym)
        lo_deallocate(xn)
        lo_deallocate(yn)
        lo_deallocate(x)
        lo_deallocate(y)
    end block cmf

    ! Now do the cumulative kappa vs frequency, with a proper tetrahedron integration.
    ! Moved it to it's own place to make this less cluttered.
    call pd%spectral_kappa(uc, qp, dr, mw, mem, spec_kappa=mf%fq_kappa, &
                           spec_kappa_band=mf%fq_kappa_band, spec_kappa_atom=mf%fq_kappa_atom)

    ! Then do the whole angular momentum guy.
    call dr%phonon_angular_momentum_matrix(qp, uc, temperature, mf%angmomalpha, mw)

    call pd%spectral_angular_momentum(uc, qp, dr, temperature, mw, mem, spec_angmom=mf%fq_angmom, &
                                      spec_angmom_band=mf%fq_angmom_band, spec_angmom_atom=mf%fq_angmom_atom)

    ! calculate kappa as a function of boundary scattering, continously.
    kpvssize: block
        real(r8), dimension(:, :, :), allocatable :: kpa, kpb
        real(r8), dimension(:, :), allocatable :: bsc
        real(r8), dimension(:), allocatable :: x
        real(r8), dimension(3, 3) :: k0
        real(r8), dimension(3) :: v0, v1
        real(r8) :: minx, maxx, f0, f1, enhet, omega, cv, tau, omthres, n
        integer :: i, j, k, l, ll, q, mode
        integer :: npts

        ! Figure out the smallest and largest scalar mean free paths
        minx = lo_huge
        maxx = -lo_huge
        do q = 1, qp%n_irr_point
        do mode = 1, dr%n_mode
            f0 = dr%iq(q)%scalar_mfp(mode)
            if (f0 .gt. lo_tol) then
                minx = min(minx, f0)
                maxx = max(maxx, f0)
            end if
        end do
        end do
        minx = minx*1E-1_r8
        maxx = maxx*1E6_r8

        ! Get an axis for this
        npts = np*2
        lo_allocate(x(npts))
        call lo_logspace(minx*1E9_r8, maxx*1E9_r8, x)
        x = x*1E-9_r8
        lo_allocate(bsc(dr%n_mode, qp%n_irr_point))
        lo_allocate(kpa(3, 3, npts))
        lo_allocate(kpb(3, 3, npts))
        bsc = 0.0_r8
        kpa = 0.0_r8
        kpb = 0.0_r8

        ! Calculate kappa, where I for each q add an extra scattering rate.
        omthres = dr%omega_min*0.5_r8
        do ll = 1, npts
            ! should be parallel over mfp points
            if (mod(ll, mw%n) .ne. mw%r) cycle
            ! slightly different depending on wether i
            if (allocated(dr%iq(1)%Fn)) then
                ! In this case we have run self-consistent calculations.
                do i = 1, qp%n_irr_point
                do mode = 1, dr%n_mode
                    omega = dr%iq(i)%omega(mode)
                    n = lo_planck(temperature, omega)
                    f1 = n*(n + 1)*norm2(dr%iq(i)%vel(:, mode))/x(ll)
                    f0 = dr%iq(i)%Qs(mode)
                    if (f0 + f1 .gt. lo_freqtol) then
                        bsc(mode, i) = f0/(f0 + f1)
                    else
                        bsc(mode, i) = 1.0_r8
                    end if
                end do
                end do

                enhet = 1.0_r8/(lo_kb_Hartree*temperature*uc%volume)
                k0 = 0.0_r8
                do i = 1, qp%n_full_point
                    k = qp%ap(i)%operation_from_irreducible
                    do j = 1, dr%n_mode
                        ! Skip gamma for acoustic branches
                        if (dr%aq(i)%omega(j) .lt. omthres) cycle
                        ! Which operation takes this point from the wedge to here
                        l = qp%ap(i)%irreducible_index
                        ! Rotate things to this points. Negative is the time reversal thingy, but does not really matter here.
                        if (k .gt. 0) then
                            v0 = lo_operate_on_vector(uc%sym%op(k), dr%iq(l)%Fn(:, j), reciprocal=.true.)
                            v1 = lo_operate_on_vector(uc%sym%op(k), dr%iq(l)%vel(:, j), reciprocal=.true.)
                        else
                            v0 = -lo_operate_on_vector(uc%sym%op(abs(k)), dr%iq(l)%Fn(:, j), reciprocal=.true.)
                            v1 = -lo_operate_on_vector(uc%sym%op(abs(k)), dr%iq(l)%vel(:, j), reciprocal=.true.)
                        end if
                        ! Get kappa for this q-point
                        omega = dr%iq(l)%omega(j)
                        n = lo_planck(temperature, omega)
                        f0 = omega*(n + 1)*n
                        k0 = k0 + enhet*f0*lo_outerproduct(v0, v1)*bsc(j, l)*qp%ap(i)%integration_weight
                    end do
                end do
            else
                ! Non-selfconsistent, get the additional boundary scattering
                do i = 1, qp%n_irr_point
                do mode = 1, dr%n_mode
                    bsc(mode, i) = 0.5_r8*norm2(dr%iq(i)%vel(:, mode))/x(ll)
                end do
                end do
                ! Calculate kappa with this added scattering
                enhet = 1.0_r8/uc%volume
                k0 = 0.0_r8
                do i = 1, qp%n_full_point
                    k = qp%ap(i)%operation_from_irreducible
                    l = qp%ap(i)%irreducible_index
                    do j = 1, dr%n_mode
                        ! Skip gamma for acoustic branches
                        if (dr%aq(i)%omega(j) .lt. omthres) cycle
                        ! Which operation takes this point from the wedge to here
                        if (k .gt. 0) then
                            v0 = lo_operate_on_vector(uc%sym%op(k), dr%iq(l)%vel(:, j), reciprocal=.true.)
                        else
                            v0 = -lo_operate_on_vector(uc%sym%op(abs(k)), dr%iq(l)%vel(:, j), reciprocal=.true.)
                        end if
                        omega = dr%iq(l)%omega(j)
                        cv = lo_harmonic_oscillator_cv(temperature, omega)
                        tau = 0.5_r8/(dr%iq(l)%linewidth(j) + bsc(j, l))
                        k0 = k0 + lo_outerproduct(v0, v0)*cv*tau*enhet*qp%ap(i)%integration_weight
                    end do
                end do
            end if
            kpa(:, :, ll) = k0
        end do

        ! Add it up over ranks
        call mpi_allreduce(kpa, kpb, npts*9, MPI_DOUBLE_PRECISION, MPI_SUM, mw%comm, mw%error)
        ! nice to store this
        lo_allocate(mf%boundary_xaxis(npts))
        lo_allocate(mf%boundary_kappa(3, 3, npts))
        mf%boundary_xaxis = x
        mf%boundary_kappa = kpb
    end block kpvssize
end subroutine

!> print the cumulative plots
subroutine write_cumulative_plots(mf, pd, uc, enhet, filename, verbosity)
    !> plots
    type(lo_mfp), intent(in) :: mf
    !> phonon dos
    type(lo_phonon_dos), intent(in) :: pd
    !> unitcell
    type(lo_crystalstructure), intent(in) :: uc
    !> frequency unit
    character(len=3) :: enhet
    !> filename
    character(len=*), intent(in) :: filename
    !> talk a lot?
    integer, intent(in) :: verbosity

    type(lo_hdf5_helper) :: h5
    real(r8), dimension(:, :, :), allocatable :: d3
    real(r8), dimension(:, :), allocatable :: d2
    real(r8), dimension(:), allocatable :: d1
    real(r8), dimension(3, 3) :: m0
    real(r8) :: unitfactor
    integer :: i, j, k, t, nt, nf, na, nb
    character(len=1000) :: gname, omstr, dosstr

    ! integer :: n_species
    ! integer :: n_unique
    ! integer, dimension(:), allocatable :: map_to_unique
    ! character(len=20), dimension(:), allocatable :: name_of_unique

    if (verbosity .gt. 0) write (*, *) 'Writing cumulative plots'

    ! The unit thingy for frequencies
    select case (enhet)
    case ('thz')
        unitfactor = lo_frequency_hartree_to_THz
        omstr = 'THz'
        dosstr = 'States/THz'
    case ('mev')
        unitfactor = lo_frequency_hartree_to_meV
        omstr = 'meV'
        dosstr = 'States/meV'
    case ('icm')
        unitfactor = lo_frequency_hartree_to_icm
        omstr = 'cm^-1'
        dosstr = 'States/cm^-1'
    case default
        call lo_stop_gracefully(['Unknown unit'], lo_exitcode_param, __FILE__, __LINE__)
    end select

    ! xy(1)='xx'
    ! xy(2)='xy'
    ! xy(3)='xz'
    ! xy(4)='yx'
    ! xy(5)='yy'
    ! xy(6)='yz'
    ! xy(7)='zx'
    ! xy(8)='zy'
    ! xy(9)='zz'

    ! Number of temperatures
    nt = size(mf%temp)

    ! Get the uniqueness contraction
    !call return_mapping_to_unique(uc,n_unique,map_to_unique,name_of_unique)

    ! Initialize hdf5
    call h5%init(__FILE__, __LINE__)
    call h5%open_file('write', trim(filename))

    ! some metadata, first how many temperature there are
    call h5%store_attribute(nt, h5%file_id, 'number_of_temperatures')
    ! then which species there are
    omstr = ""
    do i = 1, uc%nelements
        omstr = trim(omstr)//' '//trim(uc%atomic_symbol(i))
    end do
    omstr = trim(adjustl(omstr))
    call h5%store_attribute(omstr, h5%file_id, 'atomic_species')

    ! store the phonon DOS
    lo_allocate(d1(pd%n_dos_point))
    d1 = pd%omega*unitfactor
    call h5%store_data(d1, h5%file_id, 'frequencies', enhet=trim(omstr), dimensions='frequency')
    d1 = pd%dos/unitfactor
    call h5%store_data(d1, h5%file_id, 'dos', enhet=trim(dosstr), dimensions='frequency')
    lo_deallocate(d1)

    lo_allocate(d2(pd%n_dos_point, uc%na))
    d2 = pd%pdos_site/unitfactor
    call h5%store_data(d2, h5%file_id, 'dos_per_site', enhet=trim(dosstr), dimensions='site,frequency')
    lo_deallocate(d2)

    lo_allocate(d2(pd%n_dos_point, pd%n_mode))
    d2 = pd%pdos_mode/unitfactor
    call h5%store_data(d2, h5%file_id, 'dos_per_mode', enhet=trim(dosstr), dimensions='mode,frequency')
    lo_deallocate(d2)

    lo_allocate(d2(pd%n_dos_point, uc%nelements))
    d2 = 0.0_r8
    do i = 1, uc%na
        j = uc%species(i)
        d2(:, j) = d2(:, j) + pd%pdos_site(:, i)/unitfactor
    end do
    call h5%store_data(d2, h5%file_id, 'dos_per_species', enhet=trim(dosstr), dimensions='species,frequency')
    lo_deallocate(d2)

    if (verbosity .gt. 0) write (*, *) '... stored phonon DOS'

    ! create a group per temperature and stuff things there.
    do t = 1, nt
        gname = 'temperature_'//tochar(t)
        if (verbosity .gt. 0) write (*, *) '... storing temperature ', tochar(t)
        call h5%open_group('write', trim(gname))
        ! Store the temperature
        call h5%store_attribute(mf%temp(t)%temperature, h5%group_id, 'temperature')

        ! store cumulative kappa vs mean free path
        call h5%store_data(mf%temp(t)%mfaxis*lo_bohr_to_m, h5%group_id, 'mean_free_path_axis', enhet='m', dimensions='mfp')
        call h5%store_data(mf%temp(t)%mf_kappa*lo_kappa_au_to_SI, h5%group_id, &
                           'cumulative_kappa_vs_mean_free_path_total', enhet='W/mK', dimensions='mfp')
        call h5%store_data(mf%temp(t)%mf_kappa_band*lo_kappa_au_to_SI, h5%group_id, &
                           'cumulative_kappa_vs_mean_free_path_per_mode', enhet='W/mK', dimensions='mode,mfp')
        call h5%store_data(mf%temp(t)%mf_kappa_atom*lo_kappa_au_to_SI, h5%group_id, &
                           'cumulative_kappa_vs_mean_free_path_per_atom', enhet='W/mK', dimensions='atom,mfp')
        lo_allocate(d2(size(mf%temp(t)%mfaxis), uc%nelements))
        d2 = 0.0_r8
        do i = 1, uc%na
            j = uc%species(i)
            d2(:, j) = d2(:, j) + mf%temp(t)%mf_kappa_atom(:, i)*lo_kappa_au_to_SI
        end do
        call h5%store_data(d2, h5%group_id, &
                           'cumulative_kappa_vs_mean_free_path_per_species', enhet='W/mK', dimensions='species,mfp')
        lo_deallocate(d2)

        ! store cumulative kappa vs frequency
        nf = size(mf%temp(t)%fq_kappa, 1)
        nb = size(mf%temp(t)%fq_kappa_band, 2)
        na = size(mf%temp(t)%fq_kappa_atom, 2)

        ! first the frequency axis
        call h5%store_data(pd%omega*unitfactor, h5%group_id, 'frequency_axis', enhet=trim(omstr), dimensions='frequency')

        ! total kappa vs frequency
        lo_allocate(d1(nf))
        do i = 1, nf
            m0 = reshape(mf%temp(t)%fq_kappa(i, :), [3, 3])
            d1(i) = lo_trace(m0)/3.0_r8/unitfactor
        end do
        d1 = d1*lo_kappa_au_to_SI
        call h5%store_data(d1, h5%group_id, 'spectral_kappa_vs_frequency_total', enhet='W/mK/THz', dimensions='frequency')
        lo_deallocate(d1)

        ! per direction
        lo_allocate(d3(nf, 3, 3))
        d3 = 0.0_r8
        do i = 1, nf
            d3(i, :, :) = reshape(mf%temp(t)%fq_kappa(i, :), [3, 3])
        end do
        d3 = d3/unitfactor
        d3 = d3*lo_kappa_au_to_SI
        call h5%store_data(d3, h5%group_id, 'spectral_kappa_vs_frequency_per_direction', &
                           enhet='W/mK/THz', dimensions='xyz,xyz,frequency')
        lo_deallocate(d3)

        ! per mode, average over directions
        lo_allocate(d2(nf, nb))
        d2 = 0.0_r8
        do j = 1, nb
        do i = 1, nf
            d2(i, j) = norm2(mf%temp(t)%fq_kappa_band(i, j, :))/unitfactor/sqrt(3.0_r8)
        end do
        end do
        d2 = d2*lo_kappa_au_to_SI
        call h5%store_data(d2, h5%group_id, 'spectral_kappa_vs_frequency_per_mode', &
                           enhet='W/mK/THz', dimensions='mode,frequency')
        lo_deallocate(d2)

        ! per atom, average over directions
        lo_allocate(d2(nf, na))
        d2 = 0.0_r8
        do j = 1, na
        do i = 1, nf
            d2(i, j) = norm2(mf%temp(t)%fq_kappa_atom(i, j, :))/unitfactor/sqrt(3.0_r8)
        end do
        end do
        d2 = d2*lo_kappa_au_to_SI
        call h5%store_data(d2, h5%group_id, 'spectral_kappa_vs_frequency_per_atom', &
                           enhet='W/mK/THz', dimensions='atom,frequency')
        lo_deallocate(d2)

        ! per species, average over directions
        lo_allocate(d2(nf, uc%nelements))
        d2 = 0.0_r8
        do j = 1, na
        do i = 1, nf
            k = uc%species(j)
            d2(i, k) = d2(i, k) + norm2(mf%temp(t)%fq_kappa_atom(i, j, :))/unitfactor/sqrt(3.0_r8)
        end do
        end do
        d2 = d2*lo_kappa_au_to_SI
        call h5%store_data(d2, h5%group_id, 'spectral_kappa_vs_frequency_per_species', &
                           enhet='W/mK/THz', dimensions='species,frequency')
        lo_deallocate(d2)

        ! store the boundary scattering things
        call h5%store_data(mf%temp(t)%boundary_xaxis*lo_bohr_to_m, h5%group_id, &
                           'boundary_scattering_lengths', enhet='m', dimensions='domainsize')
        call h5%store_data(mf%temp(t)%boundary_kappa*lo_kappa_au_to_SI, h5%group_id, &
                           'boundary_scattering_kappa', enhet='W/mK', dimensions='xyz,xyz,domainsize')

        ! Angular momentum matrix, per direction
        lo_allocate(d3(nf, 3, 3))
        d3 = 0.0_r8
        do i = 1, nf
            d3(i, :, :) = reshape(mf%temp(t)%fq_angmom(i, :), [3, 3])
        end do
        d3 = d3/unitfactor
        !d3=d3*lo_kappa_au_to_SI
        call h5%store_data(d3, h5%group_id, 'spectral_angmom_vs_frequency_per_direction', &
                           enhet='dunno', dimensions='xyz,xyz,frequency')
        lo_deallocate(d3)

        ! Store just the tensor
        call h5%store_data(mf%temp(t)%angmomalpha, h5%group_id, 'angular_momentum_tensor', &
                           enhet='dunno', dimensions='xyz,xyz')

        call h5%close_group()
    end do

    call h5%close_file()
    call h5%destroy(__FILE__, __LINE__)
end subroutine

!> logarithmically spaced points, slightly smarter in how the shifts/spacings are decided
subroutine logspace(minv, maxv, x, smallest_separation)
    !> lower bound
    real(r8), intent(in) :: minv
    !> upper bound
    real(r8), intent(in) :: maxv
    !> separation
    real(r8), intent(in), optional :: smallest_separation
    !> logarithmically spaced points
    real(r8), dimension(:), intent(out) :: x

    integer :: n, i
    real(r8) :: dl
    real(r8) :: f0, f1, f2, f3
    real(r8) :: shift, rng, target_dl, inc

    n = size(x, 1)
    if (n .eq. 1) then
        x(1) = (minv + maxv)*0.5_r8
        return
    end if

    rng = maxv - minv
    dl = 1.0_r8/(n - 1.0_r8)
    if (present(smallest_separation)) then
        target_dl = smallest_separation
    else
        target_dl = rng/(n*10)
    end if
    ! decide on a reasonable shift
    shift = 1.0_r8
    inc = 2.0_r8

    shift = lo_sqtol
    ! what separation do I get between the first two points from this?
    do i = 1, 100
        f0 = log(shift)
        f1 = log(rng + shift)
        f2 = (f1 - f0)*dl
        f3 = exp(f0 + f2) - exp(f0)

        if (abs(f3 - target_dl) .lt. 1E-14_r8) then
            exit
        end if

        if (f3 .lt. target_dl) then
            shift = shift*inc
        else
            shift = shift/inc
            inc = (inc + 1)*0.5_r8
        end if
    end do

    f0 = log(shift)
    f1 = log(rng + shift)
    do i = 1, n
        x(i) = (f1 - f0)*dl*(i - 1) + f0
        x(i) = exp(x(i)) - shift + minv
    end do
    x(1) = minv
    x(n) = maxv
end subroutine

! !> Get reduction to unique atoms
! subroutine return_mapping_to_unique(uc,n_unique,map_to_unique,name_of_unique)
!     !> crystal structure
!     type(lo_crystalstructure), intent(in) :: uc
!     !> how many unique atoms
!     integer, intent(out) :: n_unique
!     !> map to unique
!     integer, dimension(:), allocatable, intent(out) :: map_to_unique
!     !> name of unique
!     character(len=20), dimension(:), allocatable, intent(out) :: name_of_unique
!
!     integer, dimension(uc%na) :: di0,di1
!     integer, dimension(:), allocatable :: dj0
!     integer :: a1,a2,o,i,j,k,l
!
!     di0=1
!     do a1=1,uc%na
!         di1(a1)=a1
!     enddo
!     do a1=1,uc%na
!         if ( di0(a1) .eq. 0 ) cycle
!         do o=1,uc%sym%n
!             a2=uc%sym%op(o)%fmap(a1)
!             if ( a2 .gt. a1 ) then
!                 di0(a2)=0
!                 di1(a2)=a1
!             endif
!         enddo
!     enddo
!     call lo_return_unique(di1,dj0)
!     n_unique=size(dj0)
!     ! Now I know how many there are, then get how to map things there
!     lo_allocate(map_to_unique(uc%na))
!     map_to_unique=0
!     atl1: do a1=1,uc%na
!         do o=1,uc%sym%n
!             a2=uc%sym%op(o)%fmap(a1)
!             do i=1,n_unique
!                 if ( a2 .eq. dj0(i) ) then
!                     map_to_unique(a1)=i
!                     cycle atl1
!                 endif
!             enddo
!         enddo
!     enddo atl1
!     ! Finally, get the names?
!     lo_allocate(name_of_unique(n_unique))
!     name_of_unique=''
!     do i=1,n_unique
!         name_of_unique(i)=uc%atomic_symbol(uc%species(dj0(i)))
!         write(*,*) i,name_of_unique(i)
!     enddo
!
!     l=0
!     do i=1,uc%na
!         write(*,*) i,di0(i),di1(i),map_to_unique(i)
!     enddo
!
! stop
!
! end subroutine

end module

! !> calculate the thin film suppression thing that Austin likes.
! subroutine calculate_thinfilm_kappa(mf,qp,dr,np,uc,temperature,filmnormal,tempgradient)
!     !> the cumulative plot
!     type(lo_mfp_thinfilm), intent(inout) :: mf
!     !> q-grid
!     type(lo_fft_mesh), intent(in) :: qp
!     !> dispersions, including kappa
!     type(lo_phonon_dispersions), intent(in) :: dr
!     !> how many points
!     integer, intent(in) :: np
!     !> unit cell
!     type(lo_crystalstructure), intent(in) :: uc
!     !> temperature
!     real(r8), intent(in) :: temperature
!     !> what is the normal of the film?
!     real(r8), dimension(3), intent(in) :: filmnormal,tempgradient
!
!     real(r8), dimension(:,:,:), allocatable :: Lambda
!     real(r8), dimension(:,:), allocatable :: lky,lkz,prekappa
!     real(r8) :: minLambda,maxLambda,d,id,bigC,bigD,bigE,f0
!     real(r8) :: enhet,omega,n,omthres,mfpthres,lwtol
!     integer, dimension(:,:), allocatable :: pm
!     integer :: i,j,k,l
!
!     omthres=dr%omega_min*0.5_r8    ! tolerance for zero frequency
!     mfpthres=lo_tol*1E-9_r8        ! tolerance for zero mfp
!     lwtol=lo_tol*1E-12_r8*lo_twopi ! tolerance for zero linewidth
!
!     ! Calculate all the mean free paths
!     lo_allocate(Lambda(3,dr%n_mode,dr%nq_tot))
!     Lambda=0.0_r8
!     minLambda=lo_huge
!     maxLambda=0.0_r8
!     do i=1,dr%nq_tot
!         k=qp%ap(i)%irreducible_index
!         do j=1,dr%n_mode
!             if ( dr%aq(i)%omega(j) .gt. omthres ) then
!                 f0=dr%iq(k)%linewidth(j)
!                 if ( f0 .gt. lwtol ) then
!                     f0=1.0_r8/f0
!                 else
!                     f0=0.0_r8
!                 endif
!                 Lambda(:,j,i)=dr%aq(i)%vel(:,j)*f0
!                 if ( norm2(Lambda(:,j,i)) .gt. mfpthres ) minLambda=min(minLambda,norm2(Lambda(:,j,i)))
!                 maxLambda=max(maxLambda,norm2(Lambda(:,j,i)))
!             else
!                 Lambda(:,j,i)=0.0_r8
!             endif
!         enddo
!     enddo
!
!     ! Precalculate everything that does not depend on film thickness.
!     lo_allocate(lky(dr%n_mode,dr%nq_tot))
!     lo_allocate(lkz(dr%n_mode,dr%nq_tot))
!     lo_allocate(prekappa(dr%n_mode,dr%nq_tot))
!     lo_allocate(pm(dr%n_mode,dr%nq_tot))
!     pm=0
!     lky=0.0_r8
!     lkz=0.0_r8
!     prekappa=0.0_r8
!
!     ! Unit for kappa, along with some other constant factors.
!     enhet=lo_hbar_J/lo_kb_J
!     enhet=enhet*1E30_r8/uc%volume
!     enhet=enhet*lo_hbar_J
!     enhet=enhet/(temperature*temperature*qp%n_full_point)
!
!     do i=1,dr%nq_tot
!     do j=1,dr%n_mode
!         ! figure out if this is positive or negative or 0
!         lky(j,i)=dot_product(tempgradient,Lambda(:,j,i))
!         lkz(j,i)=dot_product(filmnormal,Lambda(:,j,i))
!         if ( lky(j,i) .lt. -mfpthres ) then
!             pm(j,i)=1
!         elseif ( lky(j,i) .gt. mfpthres ) then
!             pm(j,i)=2
!         else
!             pm(j,i)=3
!         endif
!         ! the common prefactor thing
!         omega=dr%aq(i)%omega(j)
!         n=lo_planck(temperature,omega)
!         f0=omega*omega*(n+1)*n
!         if ( dr%aq(i)%omega(j) .gt. omthres ) then
!             prekappa(j,i)=enhet*f0*lkz(j,i)*dot_product(dr%aq(i)%vel(:,j),filmnormal)
!         else
!             prekappa(j,i)=0.0_r8
!         endif
!     enddo
!     enddo
!
!     ! Get an axis with relevant film thicknesses
!     lo_allocate(mf%thicknessaxis(np))
!     call lo_logspace(minLambda*1E9_r8/10,100*maxLambda*1E9_r8,mf%thicknessaxis)
!     mf%thicknessaxis=mf%thicknessaxis*1E-9_r8
!     lo_allocate(mf%kc_tot(np))
!     lo_allocate(mf%kc_band(dr%n_mode,np))
!     mf%kc_tot=0.0_r8
!     mf%kc_band=0.0_r8
!
!     ! Calculate and store things
!     mf%filmnormal=filmnormal
!     mf%tempgradient=tempgradient
!     mf%temperature=temperature
!     do l=1,np
!         d=mf%thicknessaxis(l)
!         id=1.0_r8/d
!         do j=1,dr%n_mode
!             bigC=0.0_r8
!             bigD=0.0_r8
!             bigE=0.0_r8
!             do i=1,dr%nq_tot
!                 select case(pm(j,i))
!                 case(1)
!                     f0=(lky(j,i)*exp(d/lky(j,i))-lky(j,i)-d)*id
!                     bigD=bigD-prekappa(j,i)*f0
!                 case(2)
!                     f0=( -lky(j,i)*exp(-d/lky(j,i))+lky(j,i)-d )*id
!                     bigC=bigC-prekappa(j,i)*f0
!                 case(3)
!                     bigE=bigE+prekappa(j,i)
!                 end select
!             enddo
!             mf%kc_band(j,l)=bigC+bigD+bigE
!         enddo
!         mf%kc_tot(l)=sum(mf%kc_band(:,l))
!     enddo
! end subroutine
