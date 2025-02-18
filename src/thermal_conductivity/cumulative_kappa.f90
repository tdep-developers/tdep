#include "precompilerdefinitions"
module cumulative_kappa
use konstanter, only: r8, lo_huge, lo_kappa_SI_to_au, lo_kappa_au_to_SI, lo_sqtol, &
                      lo_groupvel_ms_to_Hartreebohr, lo_pi, lo_tol, lo_bohr_to_m, lo_exitcode_param, &
                      lo_frequency_hartree_to_icm, lo_frequency_hartree_to_THz, lo_frequency_hartree_to_meV
use gottochblandat, only: lo_chop, lo_logspace, lo_linear_interpolation, qsort, lo_linspace, open_file
use type_symmetryoperation, only: lo_operate_on_secondorder_tensor
use type_crystalstructure, only: lo_crystalstructure
use type_qpointmesh, only: lo_qpoint_mesh, lo_fft_mesh
use type_phonon_dispersions, only: lo_phonon_dispersions
use type_phonon_dos, only: lo_phonon_dos
use mpi_wrappers, only: lo_mpi_helper, lo_stop_gracefully
use lo_memtracker, only: lo_mem_helper
use hdf5_wrappers, only: lo_hdf5_helper

implicit none
private
public :: lo_cumulative_kappa

!> The cumulative kappa type
type lo_cumulative_kappa
    !> The temperature
    real(r8) :: temperature = -lo_huge
    !> Total kappa
    real(r8), dimension(6) :: kappa_total = -lo_huge
    !> Mean free path
    real(r8), dimension(:), allocatable :: mfaxis
    real(r8), dimension(:, :), allocatable :: mf_kappa
    real(r8), dimension(:, :, :), allocatable :: mf_kappa_band
    real(r8), dimension(:, :, :), allocatable :: mf_kappa_atom
    !> Frequency resolved kappa
    real(r8), dimension(:, :), allocatable :: fq_kappa
    real(r8), dimension(:, :, :), allocatable :: fq_kappa_band
    real(r8), dimension(:, :, :), allocatable :: fq_kappa_atom

contains
    procedure :: get_cumulative_kappa
    procedure :: get_spectral_kappa
    procedure :: write_to_hdf5
end type

contains
!> Get the cumulative kappa
subroutine get_cumulative_kappa(mf, qp, dr, uc, np, temperature, sigma, mw, mem)
    !> The cumulative plot
    class(lo_cumulative_kappa), intent(out) :: mf
    !> The q-grid
    class(lo_qpoint_mesh), intent(in) :: qp
    !> The dispersion
    type(lo_phonon_dispersions), intent(in) :: dr
!   !> Phonon density of states
!   type(lo_phonon_dos), intent(inout) :: pd
    !> The structure
    type(lo_crystalstructure), intent(in) :: uc
    !> the number of points on the x-axis
    integer, intent(in) :: np
    !> the current temperature
    real(r8), intent(in) :: temperature
    !> additional smearing for mean free path plots?
    real(r8), intent(in) :: sigma
!   !> total kappa
!   real(r8), dimension(3, 3) :: kappa
    !> MPI helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory helper
    type(lo_mem_helper), intent(inout) :: mem


    real(r8), parameter :: kappatol = 1E-5_r8*lo_kappa_SI_to_au
    real(r8), dimension(:, :, :, :), allocatable :: kappa_chop, kappa_flat
    real(r8) :: totkappa

    mf%temperature = temperature

    call mem%allocate(kappa_chop, [qp%n_irr_point, dr%n_mode, uc%na, 6], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    kappa_chop = 0.0_r8
    ! Chop up kappa into irreducible contributions
    chopkappa: block
        real(r8), dimension(:, :, :), allocatable :: kf, kappa_flat
        real(r8), dimension(3, 3) :: buf
        complex(r8), dimension(3) :: cv0
        real(r8), dimension(uc%na) :: prj
        real(r8), dimension(3, 3) :: k0
        real(r8), dimension(6) :: f0
        integer :: iq, atom, imode, i, aq, k

        allocate(kf(6, qp%n_irr_point, dr%n_mode))
        kf = 0.0_r8
        do aq=1, qp%n_full_point
            if (mod(aq, mw%n) .ne. mw%r) cycle
            k = qp%ap(aq)%operation_from_irreducible
            iq = qp%ap(aq)%irreducible_index
            do imode=1, dr%n_mode
                buf = lo_operate_on_secondorder_tensor(uc%sym%op(k), dr%iq(iq)%kappa(:, :, imode))
                kf(1, iq, imode) = kf(1, iq, imode) + buf(1, 1)
                kf(2, iq, imode) = kf(2, iq, imode) + buf(2, 2)
                kf(3, iq, imode) = kf(3, iq, imode) + buf(3, 3)
                kf(4, iq, imode) = kf(4, iq, imode) + buf(2, 3)
                kf(5, iq, imode) = kf(5, iq, imode) + buf(1, 3)
                kf(6, iq, imode) = kf(6, iq, imode) + buf(1, 2)
            end do
        end do
        call mw%allreduce('sum', kf)

        do iq=1, qp%n_irr_point
            if (mod(iq, mw%n) .ne. mw%r) cycle
            do imode=1, dr%n_mode
                prj = 0.0_r8
                do atom=1, uc%na
                    cv0 = dr%iq(iq)%egv((atom-1)*3+1:atom*3, imode)
                    prj(atom) = prj(atom) + abs(dot_product(conjg(cv0), cv0))
                end do
                prj = lo_chop(prj/sum(prj), lo_sqtol)
                f0 = kf(:, iq, imode)
                do atom=1, uc%na
                    kappa_chop(iq, imode, atom, :) = f0(:) * prj(atom) / qp%n_full_point
                end do
            end do
        end do
        call mw%allreduce('sum', kappa_chop)
        do i=1, 6
            mf%kappa_total(i) = sum(kappa_chop(:, :, :, i))
        end do
        mf%kappa_total = lo_chop(mf%kappa_total, kappatol)

        deallocate(kf)
    end block chopkappa

    cmf: block
        real(r8), parameter :: mfptol = 1e-12_r8 * lo_groupvel_ms_to_Hartreebohr
        integer, parameter :: refinefactor = 20
        real(r8), dimension(:, :, :, :), allocatable :: ym, yn
        real(r8), dimension(:, :, :), allocatable :: xm, xn
        real(r8), dimension(:, :), allocatable :: sitebuf
        real(r8), dimension(:, :), allocatable :: y, ybuf
        real(r8), dimension(:), allocatable :: x
        real(r8), dimension(6) :: f0
        integer, dimension(:), allocatable :: ind

        real(r8) :: minx, maxx, sig1, sig2, isig, f1, f2, prefactor
        real(r8) :: min_mfp, max_mfp
        real(r8) :: a, b, logsig
        integer :: npts, ii, jj, imode, atom, ctr, iq, i, j, k
        real(r8) :: tt0, tt1, tt2

        npts = np*refinefactor

        ! First we determine the smearing parameter
        call mem%allocate(x, qp%n_irr_point*dr%n_mode, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        min_mfp = lo_huge
        max_mfp = 0.0_r8
        i = 0
        ! Here we get all the mean free paths
        do iq=1, qp%n_irr_point
            min_mfp = min(min_mfp, minval(dr%iq(iq)%scalar_mfp, dr%iq(iq)%scalar_mfp > mfptol))
            max_mfp = max(max_mfp, maxval(dr%iq(iq)%scalar_mfp))
            do imode=1, dr%n_mode
                i = i + 1
                if (dr%iq(iq)%scalar_mfp(imode) .gt. mfptol) then
                    x(i) = log(dr%iq(iq)%scalar_mfp(imode))
                else
                    x(i) = log(mfptol)
                end if
            end do
        end do
        ! And we sort
        call qsort(x)
        a = x(ceiling(size(x)*0.5_r8))
        ii = min(size(x), ceiling(size(x)*0.98_r8))
        b = x(ii)
        logsig = 0.0_r8
        do i=1, size(x) - 1
            if (x(i) .gt. a .and. x(i) .lt. b) then
                logsig = max(logsig, x(i+1) - x(i))
            end if
        end do
        logsig = logsig * sigma
        call mem%deallocate(x, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)

        ! We allocate some temporary buffers
        call mem%allocate(ym, [npts, dr%n_mode, uc%na, 6], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(yn, [npts, dr%n_mode, uc%na, 6], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(xm, [npts, dr%n_mode, uc%na], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(xn, [npts, dr%n_mode, uc%na], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(x, qp%n_irr_point, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(y, [qp%n_irr_point, 6], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(ybuf, [qp%n_irr_point, 6], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(ind, qp%n_irr_point, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        xm = 0.0_r8
        ym = 0.0_r8
        xn = 0.0_r8
        yn = 0.0_r8
        x = 0.0_r8
        y = 0.0_r8
        ybuf = 0.0_r8
        ind = 0

        ctr = 0
        do atom=1, uc%na
        do imode=1, dr%n_mode
            ctr = ctr + 1
            if(mod(ctr, mw%n) .ne. mw%r) cycle

            ! Here we get the kappa, and sort the values according to mean free paths
            y = 0.0_r8
            do iq=1, qp%n_irr_point
                x(iq) = dr%iq(iq)%scalar_mfp(imode)
                ybuf(iq, :) = kappa_chop(iq, imode, atom, :)
            end do
            call qsort(x, ind)
            ! We sort it for each directions
            do ii=1, 6
                y(:, ii) = ybuf(ind, ii)
            end do
            minx = minval(x, x > lo_tol) * 0.1_r8
            maxx = maxval(x) * 1.5_r8

            ! Now we will calculate the integrated mean free path
            call lo_linspace(minx, maxx, xm(:, imode, atom))

            do iq=1, qp%n_irr_point
                if (x(iq) .lt. minx) cycle
                if (x(iq) .gt. maxx) cycle

                a = log(x(iq))
                b = logsig
                sig1 = exp(a - 4 * b) - x(iq)
                sig2 = exp(a + 4 * b) - x(iq)
                prefactor = sqrt(b * lo_pi) * exp(a + b * 0.25_r8)
                isig = 1.0_r8 / b

                ii = max(1, floor(npts * (x(iq) + sig1 - minx) / maxx))
                jj = min(npts, ceiling(npts * (x(iq) + sig2 - minx) / maxx))
                if (jj .le. ii) cycle
                do k=ii, jj
                    f1 = (1.0_r8 + erf((log(xm(k, imode, atom)) - log(x(iq))) * isig)) * 0.5_r8
                    ym(k, imode, atom, :) = ym(k, imode, atom, :) + f1 * y(iq, :)
                end do
                do k=jj+1, npts
                    ym(k, imode, atom, :) = ym(k, imode, atom, :) + y(iq, :)
                end do
            end do
        end do
        end do
!       call mw%allreduce('sum', xm, xn)
!       call mw%allreduce('sum', ym, yn)
        call mw%allreduce('sum', xm)
        call mw%allreduce('sum', ym)

        allocate(mf%mfaxis(np))
        allocate(mf%mf_kappa(np, 6))
        allocate(mf%mf_kappa_band(np, dr%n_mode, 6))
        allocate(mf%mf_kappa_atom(np, uc%na, 6))
        mf%mf_kappa = 0.0_r8
        mf%mf_kappa_band = 0.0_r8
        mf%mf_kappa_atom = 0.0_r8
        minx = minval(xm)
        maxx = maxval(xm) * 2.0_r8
        call lo_logspace(minx, maxx, mf%mfaxis)

        do imode=1, dr%n_mode
        do atom=1, uc%na
        do i=1, np
            if (mf%mfaxis(i) .lt. xm(1, imode, atom)*1.01_r8) then
                f0 = 0.0_r8
            elseif (mf%mfaxis(i) .gt. xm(npts, imode, atom)*0.999_r8) then
                f0(:) = ym(npts, imode, atom, :)
            else
                do ii=1, 6
                    f0(ii) = lo_linear_interpolation(xm(:, imode, atom), ym(:, imode, atom, ii), mf%mfaxis(i))
                end do
            end if
            mf%mf_kappa(i, :) = mf%mf_kappa(i, :) + f0(:)
            mf%mf_kappa_band(i, imode, :) = mf%mf_kappa_band(i, imode, :) + f0(:)
            mf%mf_kappa_atom(i, atom, :) = mf%mf_kappa_atom(i, atom, :) + f0(:)
        end do
        end do
        end do

        ! Make sure it's monotically increasing
        do i=1, np-1
            do ii=1, 6
                if (mf%mf_kappa(i+1, ii) .lt. mf%mf_kappa(i, ii)) then
                    mf%mf_kappa(i+1, ii) = mf%mf_kappa(i, ii)
                end if
                do atom=1, uc%na
                    if (mf%mf_kappa_atom(i+1, atom, ii) .lt. mf%mf_kappa_atom(i, atom, ii)) then
                        mf%mf_kappa_atom(i+1, atom, ii) = mf%mf_kappa_atom(i, atom, ii)
                    end if
                end do
                do imode=1, dr%n_mode
                    if (mf%mf_kappa_band(i+1, imode, ii) .lt. mf%mf_kappa_band(i, imode, ii)) then
                        mf%mf_kappa_band(i+1, imode, ii) = mf%mf_kappa_band(i, imode, ii)
                    end if
                end do
            end do
        end do

        ! Fix site degeneracy
        call mem%allocate(sitebuf, [np, uc%na], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        do ii=1, 6
            sitebuf = 0.0_r8
            do i=1, uc%na
                do j=1, uc%sym%degeneracy(i)
                    k = uc%sym%degenerate_atom(j, i)
                    sitebuf(:, i) = sitebuf(:, i) + mf%mf_kappa_atom(:, k, ii) / real(uc%sym%degeneracy(j), r8)
                end do
            end do
            mf%mf_kappa_atom(:, :, ii) = sitebuf
        end do
        call mem%deallocate(sitebuf, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)

        ! Now make sure everything is normalized properly
        do i=1, np
            do ii=1, 6
                ! Fix the band decomposed
                f1 = sum(mf%mf_kappa_band(i, :, ii))
                f2 = mf%mf_kappa(i, ii)
                if (f1 .gt. kappatol .and. f2 .gt. kappatol) then
                    mf%mf_kappa_band(i, :, ii) = mf%mf_kappa_band(i, :, ii) * mf%mf_kappa(i, ii) / f1
                else
                    mf%mf_kappa_band(i, :, ii) = 0.0_r8
                end if
                ! Fix the atom decomposed
                f1 = sum(mf%mf_kappa_atom(i, :, ii))
                f2 = mf%mf_kappa(i, ii)
                if (f1 .gt. kappatol .and. f2 .gt. kappatol) then
                    mf%mf_kappa_atom(i, :, ii) = mf%mf_kappa_atom(i, :, ii) * mf%mf_kappa(i, ii) / f1
                else
                    mf%mf_kappa_atom(i, :, ii) = 0.0_r8
                end if
            end do
        end do

        call mem%deallocate(ym, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(yn, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(xm, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(xn, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(x, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(y, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(ybuf, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(ind, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    end block cmf

    call mem%deallocate(kappa_chop, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
end subroutine

subroutine get_spectral_kappa(mf, uc, qp, dr, pd, mw, mem)
    !> The cumulative plot
    class(lo_cumulative_kappa), intent(inout) :: mf
    !> The structure
    type(lo_crystalstructure), intent(in) :: uc
    !> The q-grid
    class(lo_qpoint_mesh), intent(in) :: qp
    !> The dispersion
    type(lo_phonon_dispersions), intent(inout) :: dr
    !> Phonon density of states
    type(lo_phonon_dos), intent(inout) :: pd
    !> MPI helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory helper
    type(lo_mem_helper), intent(inout) :: mem

    getfullgridkappa: block
        real(r8), dimension(3, 3) :: buf
        integer :: aq, iq, k, imode

        ! We need to allocate and get the thermal and conductivity on the full grid
        do aq=1, qp%n_full_point
            allocate(dr%aq(aq)%kappa(3, 3, dr%n_mode))
            iq = qp%ap(aq)%irreducible_index
            k = qp%ap(aq)%operation_from_irreducible
            do imode=1, dr%n_mode
                buf = lo_operate_on_secondorder_tensor(uc%sym%op(k), dr%iq(iq)%kappa(:, :, imode))
                dr%aq(aq)%kappa(:, :, imode) = buf
            end do
        end do

    end block getfullgridkappa

    ! Now do the cumulative kappa vs frequency, with a proper tetrahedron integration.
    ! This is called from the dedicated thing in type_phonon_dos
    call pd%spectral_kappa(uc, qp, dr, mw, mem, spec_kappa=mf%fq_kappa, &
                           spec_kappa_band=mf%fq_kappa_band, spec_kappa_atom=mf%fq_kappa_atom)
end subroutine

subroutine write_to_hdf5(mf, pd, uc, enhet, filename, mem)
    !> The cumulative plot
    class(lo_cumulative_kappa), intent(in) :: mf
    !> phonon dos
    type(lo_phonon_dos), intent(in) :: pd
    !> unitcell
    type(lo_crystalstructure), intent(in) :: uc
    !> frequency unit
    character(len=3) :: enhet
    !> filename
    character(len=*), intent(in) :: filename
    !> memory helper
    type(lo_mem_helper), intent(inout) :: mem

    type(lo_hdf5_helper) :: h5
    real(r8), dimension(:), allocatable :: d1
    real(r8), dimension(:, :), allocatable :: d2
    real(r8), dimension(:, :, :), allocatable :: d3
    real(r8) :: unitfactor
    integer :: i
    character(len=1000) :: spstr, omstr, dosstr

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

    call h5%init(__FILE__, __LINE__)
    call h5%open_file('write', trim(filename))

    spstr = ""
    do i=1, uc%nelements
        spstr = trim(spstr)// ' '//trim(uc%atomic_symbol(i))
    end do
    spstr = trim(adjustl(spstr))
    call h5%store_attribute(trim(spstr), h5%file_id, 'atomic_species')
    call h5%store_attribute(trim('temperature'), h5%file_id, 'temperature')

    ! store the phonon DOS
    call mem%allocate(d1, pd%n_dos_point, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    d1 = pd%omega * unitfactor
    call h5%store_data(d1, h5%file_id, 'frequencies', enhet=trim(omstr), dimensions='frequency')
    d1 = pd%dos / unitfactor
    call h5%store_data(d1, h5%file_id, 'dos', enhet=trim(dosstr), dimensions='frequency')
    call mem%deallocate(d1, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)

    ! Write cumulative kappa vs mfp
    call h5%store_data(mf%mfaxis*lo_bohr_to_m, h5%file_id, 'mean_free_path_axis', enhet='m', dimensions='mfp')
    call h5%store_data(mf%mf_kappa*lo_kappa_au_to_SI, h5%file_id, &
                       'cumulative_kappa_vs_mean_free_path', enhet='W/m/K', dimensions='xx,mfp')
    call h5%store_data(mf%mf_kappa_band*lo_kappa_au_to_SI, h5%file_id, &
                       'cumulative_kappa_vs_mean_free_path_per_mode', enhet='W/m/K', dimensions='xx,mode,mfp')
    call h5%store_data(mf%mf_kappa_atom*lo_kappa_au_to_SI, h5%file_id, &
                       'cumulative_kappa_vs_mean_free_path_per_atom', enhet='W/m/K', dimensions='xx,atom,mfp')

    call mem%allocate(d2, [pd%n_dos_point, 6], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    ! Get the spectral kappa in Voigt notation
    d2(:, 1) = mf%fq_kappa(:, 1)
    d2(:, 2) = mf%fq_kappa(:, 5)
    d2(:, 3) = mf%fq_kappa(:, 9)
    d2(:, 4) = mf%fq_kappa(:, 6)
    d2(:, 5) = mf%fq_kappa(:, 3)
    d2(:, 6) = mf%fq_kappa(:, 2)
    ! Write spectral kappa vs frequency
    call h5%store_data(d2*lo_kappa_au_to_SI/unitfactor, h5%file_id, &
                       'spectral_kappa_vs_frequency', enhet='W/m/K', dimensions='xx,frequency')
    call mem%deallocate(d2, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)

    call mem%allocate(d3, [pd%n_dos_point, 3*uc%na, 6], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    d3(:, :, 1) = mf%fq_kappa_band(:, :, 1)
    d3(:, :, 2) = mf%fq_kappa_band(:, :, 5)
    d3(:, :, 3) = mf%fq_kappa_band(:, :, 9)
    d3(:, :, 4) = mf%fq_kappa_band(:, :, 6)
    d3(:, :, 5) = mf%fq_kappa_band(:, :, 3)
    d3(:, :, 6) = mf%fq_kappa_band(:, :, 2)
    call h5%store_data(d3*lo_kappa_au_to_SI/unitfactor, h5%file_id, &
                       'spectral_kappa_vs_frequency_per_mode', enhet='W/m/K', dimensions='xx,frequency,mode')
    call mem%deallocate(d3, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)

    call mem%allocate(d3, [pd%n_dos_point, uc%na, 6], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    d3(:, :, 1) = mf%fq_kappa_atom(:, :, 1)
    d3(:, :, 2) = mf%fq_kappa_atom(:, :, 5)
    d3(:, :, 3) = mf%fq_kappa_atom(:, :, 9)
    d3(:, :, 4) = mf%fq_kappa_atom(:, :, 6)
    d3(:, :, 5) = mf%fq_kappa_atom(:, :, 3)
    d3(:, :, 6) = mf%fq_kappa_atom(:, :, 2)
    call h5%store_data(d3*lo_kappa_au_to_SI/unitfactor, h5%file_id, &
                       'spectral_kappa_vs_frequency_per_mode', enhet='W/m/K', dimensions='atom,frequency,xx')
    call mem%deallocate(d3, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)

    call h5%close_file()
    call h5%destroy(__FILE__, __LINE__)
end subroutine

end module
