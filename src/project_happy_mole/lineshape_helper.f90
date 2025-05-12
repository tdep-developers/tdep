module lineshape_helper
!!
!! Common helper routines.
!!
use konstanter, only: r8, i8, lo_imag, lo_pi, lo_hugeint, lo_huge, lo_iou, lo_exitcode_param, lo_freqtol, lo_tiny, lo_frequency_Hartree_to_THz
use gottochblandat, only: lo_gauss, lo_trapezoid_integration, lo_planck, lo_return_unique, lo_linspace, open_file
use hdf5_wrappers, only: lo_hdf5_helper
use mpi_wrappers, only: lo_mpi_helper, lo_stop_gracefully
use lo_memtracker, only: lo_mem_helper
use type_qpointmesh, only: lo_qpoint_mesh, lo_fft_mesh
use type_phonon_dispersions, only: lo_phonon_dispersions
use lo_brents_method, only: lo_brent_helper
use quadratures_stencils, only: lo_centraldifference, lo_gaussianquadrature
use lo_sorting, only: lo_qsort
use, intrinsic :: iso_c_binding

implicit none
include 'fftw3.f'

private
public :: diagonal_greensfunction
public :: evaluate_spectral_function
public :: taperfn_im
public :: gaussian_smear_spectral_function
public :: gaussian_smear_distribution
public :: lo_spectralfunction_helper
public :: index_on_grid
public :: interpolated_spectral_function
public :: find_spectral_function_max_and_fwhm
public :: integrate_spectral_function
public :: integrate_two_spectral_functions
public :: lo_convolution_helper

!> container for spectral functions
type lo_spectralfunction_helper
    !> number of (irreducible) q-points
    integer :: n_irr_qpoint = -lo_hugeint
    !> number of energy points
    integer :: n_energy = -lo_hugeint
    !> max energy
    real(r8) :: max_energy = -lo_huge
    !> temperature used to generate
    real(r8) :: temperature = -lo_huge
    !> coordinates of q-points
    real(r8), dimension(:, :), allocatable :: qv
    !> spectral function (energy,mode,qvector)
    real(r8), dimension(:, :, :), allocatable :: spectralfunction
    !> real self-energy (energy,mode,qvector)
    real(r8), dimension(:, :, :), allocatable :: sigma_Re
    !> imaginary self-energy (energy,mode,qvector)
    real(r8), dimension(:, :, :), allocatable :: sigma_Im
contains
    !> write to hdf5
    procedure :: write_to_hdf5 => write_spectralfunction_helper_to_hdf5
    !> read from hdf5
    procedure :: read_from_hdf5 => read_spectralfunction_helper_from_hdf5
    !> smear spectral functions
    !procedure :: smear => smear_spectralfunction
    !> size in memory
    procedure :: size_in_mem => isf_size_in_mem
end type

!> container to help with convolution things
type lo_convolution_helper
    !> Number of points is -n:n for all our quantities
    integer :: n = -lo_hugeint
    !> step in frequency axis
    real(r8) :: deltax = -lo_huge
    !> frequency axis
    real(r8), dimension(:), allocatable :: x
    !> n(omega)+1, planck factor + 1
    real(r8), dimension(:), allocatable :: planck_plus_one
    !> inverse of the above?
    real(r8), dimension(:), allocatable :: inv_planck_plus_one
    !> buffer for frequency domain spectral functions
    real(r8), dimension(:, :), allocatable :: sf2
    real(r8), dimension(:, :), allocatable :: sf3
    !> buffer for spectral functions with prefactors multiplied in
    real(r8), dimension(:, :), allocatable :: pref_sf2
    real(r8), dimension(:, :), allocatable :: pref_sf3

    !> FFTW plans
    type(c_ptr) :: plan_fwd, plan_bwd
    !> FFT buffers
    complex(c_double_complex), dimension(:), allocatable :: fft0, fft1
    !> Fourier transformed spectral function times n+1
    complex(c_double_complex), dimension(:, :), allocatable :: fft_pref_sf2
    complex(c_double_complex), dimension(:, :), allocatable :: fft_pref_sf3
    complex(c_double_complex), dimension(:, :, :), allocatable :: fft_smearing
contains
    !> create all the buffers
    procedure :: generate => generate_convolution_helper
    !> buffer the spectral functions
    procedure :: buffer_spectral_functions
    !> convolute them
    procedure :: convolute_to_sfun
    !> destroy properly
    procedure :: destroy => destroy_convolution_helper
end type

contains

subroutine destroy_convolution_helper(ch)
    !> convolution helper
    class(lo_convolution_helper), intent(out) :: ch

    call dfftw_destroy_plan(ch%plan_fwd)
    call dfftw_destroy_plan(ch%plan_bwd)
end subroutine

!> create the convolution helper
subroutine generate_convolution_helper(ch, energy_axis, temperature, n_mode)
    !> convolution helper
    class(lo_convolution_helper), intent(out) :: ch
    !> energy axis
    real(r8), dimension(:), intent(in) :: energy_axis
    !> temperature
    real(r8), intent(in) :: temperature
    !> number of modes
    integer, intent(in) :: n_mode

    integer :: i
    ! Number of points we are interested in?
    ch%n = size(energy_axis) - 1
    ! Some space for buffers
    allocate (ch%x(-ch%n:ch%n))
    allocate (ch%planck_plus_one(-ch%n:ch%n))
    allocate (ch%inv_planck_plus_one(-ch%n:ch%n))
    allocate (ch%sf2(-ch%n:ch%n, n_mode))
    allocate (ch%sf3(-ch%n:ch%n, n_mode))
    allocate (ch%pref_sf2(-ch%n:ch%n, n_mode))
    allocate (ch%pref_sf3(-ch%n:ch%n, n_mode))

    ch%x = 0.0_r8
    ch%planck_plus_one = 0.0_r8
    ch%inv_planck_plus_one = 0.0_r8
    ch%sf2 = 0.0_r8
    ch%sf3 = 0.0_r8
    ch%pref_sf2 = 0.0_r8
    ch%pref_sf3 = 0.0_r8

    ! Buffers for Fourier transforming
    allocate (ch%fft0(-ch%n:ch%n))
    allocate (ch%fft1(-ch%n:ch%n))
    allocate (ch%fft_pref_sf2(-ch%n:ch%n, n_mode))
    allocate (ch%fft_pref_sf3(-ch%n:ch%n, n_mode))
    allocate (ch%fft_smearing(-ch%n:ch%n, n_mode, n_mode))

    ! Get an x-axis, might come handy. x(0)=0
    ch%x(0) = 0.0_r8
    do i = 1, ch%n
        ch%x(i) = energy_axis(i + 1)
        ch%x(-i) = -ch%x(i)
    end do
    ch%deltax = ch%x(1) - ch%x(0)
    ! Pre-calculate the planck factors?
    do i = 1, ch%n
        ch%planck_plus_one(i) = lo_planck(temperature, ch%x(i)) + 1.0_r8
        ch%planck_plus_one(-i) = -lo_planck(temperature, ch%x(i))
    end do
    do i = 1, ch%n
        ch%inv_planck_plus_one(i) = 1.0_r8/ch%planck_plus_one(i)
        ch%inv_planck_plus_one(-i) = 1.0_r8/ch%planck_plus_one(-i)
    end do
    ! Pre-plan the FFTs?
    call dfftw_plan_dft_1d(ch%plan_fwd, size(ch%fft0), ch%fft0, ch%fft1, FFTW_FORWARD, FFTW_ESTIMATE)
    call dfftw_plan_dft_1d(ch%plan_bwd, size(ch%fft0), ch%fft0, ch%fft1, FFTW_BACKWARD, FFTW_ESTIMATE)
end subroutine

!> fetch and massage spectral functions
subroutine buffer_spectral_functions(ch, sf2, sf3)
    !> convolution helper
    class(lo_convolution_helper), intent(inout) :: ch
    !> spectral functions for q' (energy,mode)
    real(r8), dimension(:, :), intent(in) :: sf2
    !> spectral functions for q'' (energy,mode)
    real(r8), dimension(:, :), intent(in) :: sf3

    integer :: imode, i
    ! Sort the spectral functions into the double-arrays

    ch%sf2 = 0.0_r8
    ch%sf3 = 0.0_r8
    ch%pref_sf2 = 0.0_r8
    ch%pref_sf3 = 0.0_r8
    do imode = 1, size(ch%sf2, 2)
        do i = 1, ch%n
            ch%sf2(i, imode) = sf2(i + 1, imode)
            ch%sf3(i, imode) = sf3(i + 1, imode)
        end do
        do i = 1, ch%n
            ch%sf2(-i, imode) = -ch%sf2(i, imode)
            ch%sf3(-i, imode) = -ch%sf3(i, imode)
        end do
        ch%pref_sf2(:, imode) = ch%sf2(:, imode)*ch%planck_plus_one
        ch%pref_sf3(:, imode) = ch%sf3(:, imode)*ch%planck_plus_one
    end do

    ! Probably a neat idea to pre-transform the spectral functions?
    do imode = 1, size(ch%sf2, 2)
        ch%fft0 = ch%pref_sf2(:, imode)
        call dfftw_execute_dft(ch%plan_fwd, ch%fft0, ch%fft_pref_sf2(:, imode))
        ch%fft0 = ch%pref_sf3(:, imode)
        call dfftw_execute_dft(ch%plan_fwd, ch%fft0, ch%fft_pref_sf3(:, imode))
    end do

    ! While I'm at it, generate all the possible smearing parameters
    ! and Fourier transform them. Maybe not yet.
end subroutine

!> get the sfun
subroutine convolute_to_sfun(ch, imode, jmode, sigma, sfun)
    !> convolution helper
    class(lo_convolution_helper), intent(inout) :: ch
    !> which modes
    integer, intent(in) :: imode, jmode
    !> smearing parameter
    real(r8), intent(in) :: sigma
    !> s-function
    real(r8), dimension(:), intent(out) :: sfun

    sfun = 0.0_r8
    ! Multiply
    ch%fft0 = ch%fft_pref_sf2(:, imode)*ch%fft_pref_sf3(:, jmode)
    ch%fft1 = 0.0_r8
    ! Inverse-transform
    call dfftw_execute_dft(ch%plan_bwd, ch%fft0, ch%fft1)

    ! Do the FFT-shift thingy
    sfun(1:ch%n + 1) = real(ch%fft1(-ch%n:0), r8)
    !sfun(1:ch%n+1)=real(ch%fft1(0:ch%n),r8)
    ! Scale to sensible thing
    sfun = sfun*ch%deltax/real(size(ch%x), r8)
    ! Make sure it's positive

    ! Add the planck factor thing
    sfun = sfun*ch%inv_planck_plus_one(0:ch%n)
    !sfun=sfun / ch%planck_plus_one(0:ch%n)

    ! Make sure it's positive
    sfun = max(0.0_r8, sfun)
    sfun(1) = 0
    ! Smear it
    call gaussian_smear_spectral_function(ch%x(0:ch%n), sigma, sfun)
end subroutine

! !> Smear spectral functions before integrations
! subroutine smear_spectralfunction(sf, qp, dr, maxf, smearing_prefactor, mw, mem)
!     !> helper container
!     class(lo_spectralfunction_helper), intent(inout) :: sf
!     !> q-mesh
!     class(lo_qpoint_mesh), intent(in) :: qp
!     !> phonon dispersions
!     type(lo_phonon_dispersions), intent(in) :: dr
!     !> max frequency
!     real(r8), intent(in) :: maxf
!     !> smearing prefactor
!     real(r8), intent(in) :: smearing_prefactor
!     !> MPI helper
!     type(lo_mpi_helper), intent(inout) :: mw
!     !> memory helper
!     type(lo_mem_helper), intent(inout) :: mem

!     real(r8), dimension(:, :, :), allocatable :: sm_buf
!     real(r8), dimension(:), allocatable :: buf0
!     real(r8) :: sigma, f0
!     integer :: iq, imode, l

!     ! Count number of local q-points
!     l = 0
!     do iq = 1, qp%n_irr_point
!         if (mod(iq, mw%n) .ne. mw%r) cycle
!         l = l + 1
!     end do

!     if (l .gt. 0) then
!         ! Space for smeared spectral function
!         call mem%allocate(sm_buf, [sf%n_energy, dr%n_mode, l], persistent=.false., scalable=.true., file=__FILE__, line=__LINE__)
!         call mem%allocate(buf0, sf%n_energy, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
!         sm_buf = 0.0_r8
!         buf0 = 0.0_r8
!         call lo_linspace(0.0_r8, maxf, buf0)
!     else
!         ! Allocate dummies
!         call mem%allocate(sm_buf, [1, 1, 1], persistent=.false., scalable=.true., file=__FILE__, line=__LINE__)
!         call mem%allocate(buf0, 1, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
!     end if

!     ! Smear spectral functions
!     l = 0
!     do iq = 1, qp%n_irr_point
!         if (mod(iq, mw%n) .ne. mw%r) cycle
!         l = l + 1
!         do imode = 1, dr%n_mode
!             if (dr%iq(iq)%omega(imode) .gt. lo_freqtol) then
!                 sigma = qp%adaptive_sigma(qp%ip(iq)%radius, dr%iq(iq)%vel(:, imode), dr%default_smearing(imode), smearing_prefactor)
!                 sm_buf(:, imode, l) = sf%spectralfunction(:, imode, iq)

!                 call gaussian_smear_spectral_function(buf0, sigma, sm_buf(:, imode, l))

!                 f0 = lo_trapezoid_integration(buf0, sm_buf(:, imode, l))
!                 sm_buf(:, imode, l) = sm_buf(:, imode, l)/f0
!             else
!                 sm_buf(:, imode, l) = 0.0_r8
!             end if
!         end do
!     end do

!     ! Distribute
!     sf%spectralfunction = 0.0_r8
!     l = 0
!     do iq = 1, qp%n_irr_point
!         if (mod(iq, mw%n) .ne. mw%r) cycle
!         l = l + 1
!         sf%spectralfunction(:, :, iq) = sm_buf(:, :, l)
!     end do
!     call mw%allreduce('sum', sf%spectralfunction)

!     ! Cleanup
!     call mem%deallocate(sm_buf, persistent=.false., scalable=.true., file=__FILE__, line=__LINE__)
!     call mem%deallocate(buf0, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
! end subroutine

!> evaluate (diagonal) spectral function at arbitrary point
function interpolated_spectral_function(x, yim, yre, omega, scalefactor, xi) result(yi)
    !> x-axis that defines self-energy
    real(r8), dimension(:), intent(in) :: x
    !> imaginary part of self energy
    real(r8), dimension(:), intent(in) :: yim
    !> real part of self energy
    real(r8), dimension(:), intent(in) :: yre
    !> harmonic frequency
    real(r8), intent(in) :: omega
    !> global scaling factor (for normalization)
    real(r8), intent(in) :: scalefactor
    !> point to evaluate at
    real(r8), intent(in) :: xi
    !> function value
    real(r8) :: yi

    ! Find the correct interval for interpolation.
    real(r8) :: xhi, f0, f1, f2, sigmaIm, sigmaRe
    integer :: n, ilo, ihi

    n = size(x)
    xhi = x(n)

    ! This should put us in the right interval.
    ilo = floor(n*xi/xhi) + 1
    ! ilo=max(ilo,1)
    ! ilo=min(ilo,n-1)
    ihi = ilo + 1
    ! evaluate imaginary and real self-energy
    f0 = (x(ihi) - xi)/(x(ihi) - x(ilo))
    sigmaIm = f0*yim(ilo) + (1.0_r8 - f0)*yim(ihi)
    sigmaRe = f0*yre(ilo) + (1.0_r8 - f0)*yre(ihi)

    ! get the spectral function
    f0 = 2*omega*sigmaIm*2*omega/lo_pi
    f1 = (xi**2 - omega**2 - 2*omega*sigmaRe)**2
    f2 = 4*(omega**2)*(sigmaIm**2)
    yi = scalefactor*f0/(f1 + f2)

    ! try again?
    ! f0=(xi-omega-sigmaRe)**2+sigmaIm**2
    ! f1=(xi+omega+sigmaRe)**2+sigmaIm**2
    ! yi=scalefactor*sigmaIm*(1.0_r8/f0-1.0_r8/f1)/lo_pi
end function

!> read helper from file
subroutine read_spectralfunction_helper_from_hdf5(sf, filename, qp, maxf, temperature, n_energy, mw, verbosity)
    !> helper container
    class(lo_spectralfunction_helper), intent(out) :: sf
    !> filename
    character(len=*), intent(in) :: filename
    !> q-mesh
    class(lo_qpoint_mesh), intent(in) :: qp
    !> max energy
    real(r8), intent(in) :: maxf
    !> temperature
    real(r8), intent(in) :: temperature
    !> number of energy values
    integer, intent(in) :: n_energy
    !> MPI helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> talk a lot?
    integer, intent(in) :: verbosity

    type(lo_hdf5_helper) :: h5
    integer :: readrnk, n_mode, iq

    if (verbosity .gt. 0) write (lo_iou, *) '... reading spectral functions from "'//trim(filename)//'"'

    ! Which rank to read on
    readrnk = mw%n - 1
    n_mode = -1

    ! Initialize hdf5
    if (mw%r .eq. readrnk) then
        call h5%init()
        call h5%open_file('read', trim(filename))

        ! Read metadata
        call h5%read_attribute(sf%temperature, h5%file_id, 'temperature')
        call h5%read_attribute(sf%n_irr_qpoint, h5%file_id, 'n_irr_qpoint')
        call h5%read_attribute(sf%n_energy, h5%file_id, 'n_energy')
        call h5%read_attribute(sf%max_energy, h5%file_id, 'max_energy')
        ! Read actual data
        call h5%read_data(sf%qv, h5%file_id, 'qvectors')
        call h5%read_data(sf%spectralfunction, h5%file_id, 'spectralfunction')
        n_mode = size(sf%spectralfunction, 2)

        call h5%close_file()
        call h5%destroy()

        if (verbosity .gt. 0) write (lo_iou, *) '... read spectral functions'
    end if

    ! Distribute metadata
    call mw%bcast(sf%temperature, from=readrnk)
    call mw%bcast(sf%n_irr_qpoint, from=readrnk)
    call mw%bcast(sf%n_energy, from=readrnk)
    call mw%bcast(sf%max_energy, from=readrnk)
    call mw%bcast(n_mode, from=readrnk)

    ! Make space
    if (mw%r .ne. readrnk) then
        allocate (sf%spectralfunction(sf%n_energy, n_mode, sf%n_irr_qpoint))
        allocate (sf%qv(3, sf%n_irr_qpoint))
        sf%spectralfunction = 0.0_r8
        sf%qv = 0.0_r8
    end if

    ! Distribute data
    call mw%bcast(sf%spectralfunction, from=readrnk)
    call mw%bcast(sf%qv, from=readrnk)

    ! Now sanity test to make sure things agree.
    if (sf%n_irr_qpoint .ne. qp%n_irr_point) then
        call lo_stop_gracefully(['Mismatch in number of q-points'], lo_exitcode_param, __FILE__, __LINE__, mw%comm)
    end if
    if (sf%n_energy .ne. n_energy) then
        call lo_stop_gracefully(['Mismatch in number of energy values'], lo_exitcode_param, __FILE__, __LINE__, mw%comm)
    end if
    if (abs(sf%max_energy - maxf)/sf%max_energy .gt. 1E-10_r8) then
        call lo_stop_gracefully(['Mismatch in energy values'], lo_exitcode_param, __FILE__, __LINE__, mw%comm)
    end if
    do iq = 1, sf%n_irr_qpoint
        if (norm2(sf%qv(:, iq) - qp%ip(iq)%r) .gt. 1E-10_r8) then
            call lo_stop_gracefully(['Mismatch in q-points'], lo_exitcode_param, __FILE__, __LINE__, mw%comm)
        end if
    end do
    if (abs(temperature - sf%temperature) .gt. 1E-5_r8) then
        call lo_stop_gracefully(['Mismatch in temperature'], lo_exitcode_param, __FILE__, __LINE__, mw%comm)
    end if

    if (verbosity .gt. 0) write (lo_iou, *) '... distributed spectral functions'
end subroutine

!> dump the helper to file
subroutine write_spectralfunction_helper_to_hdf5(sf, filename)
    !> helper container
    class(lo_spectralfunction_helper), intent(in) :: sf
    !> filename
    character(len=*), intent(in) :: filename

    type(lo_hdf5_helper) :: h5

    ! Initialize hdf5
    call h5%init(__FILE__, __LINE__)
    call h5%open_file('write', trim(filename))

    ! Store metadata
    call h5%store_attribute(sf%temperature, h5%file_id, 'temperature')
    call h5%store_attribute(sf%n_irr_qpoint, h5%file_id, 'n_irr_qpoint')
    call h5%store_attribute(sf%n_energy, h5%file_id, 'n_energy')
    call h5%store_attribute(sf%max_energy, h5%file_id, 'max_energy')
    ! Store actual data
    call h5%store_data(sf%qv, h5%file_id, 'qvectors', enhet='fractional')
    call h5%store_data(sf%spectralfunction, h5%file_id, 'spectralfunction', enhet='atomic')

    ! Close the file
    call h5%close_file()
    call h5%destroy(__FILE__, __LINE__)
end subroutine

!> get a function that neatly tapers the imaginary part of the self-energy to sensible values.
subroutine taperfn_im(x, maxomega, harmomega, y)
    !> energy axis
    real(r8), dimension(:), intent(in) :: x
    !> max harmonic frequency
    real(r8), intent(in) :: maxomega
    !> current harmonic frequency
    real(r8), intent(in) :: harmomega
    !> tapering function
    real(r8), dimension(:), intent(out) :: y

    real(r8), parameter :: lowtaperonset = 0.25_r8
    real(r8), parameter :: taperonset = 2.1_r8
    integer :: nlowtaper

    integer :: n, ii, jj, i, j
    real(r8) :: xmax, f0

    n = size(x) ! number of energy values
    xmax = x(n) ! largest x
    y = 1.0_r8
    ! Taper the low end, with a sine?
    nlowtaper = ceiling(n*lowtaperonset*harmomega/xmax)

    y(1) = 0.0_r8
    ! do i=2,nlowtaper-1
    !     f0=real(i-1,r8)/real(nlowtaper-1,r8)
    !     f0=f0*lo_pi*0.5_r8
    !     y(i)=sin(f0)
    ! enddo

    ! Taper the high end, with a cosine?
    ii = floor(n*maxomega*taperonset/xmax)
    jj = floor(n*(maxomega*taperonset*4 + xmax)/5.0_r8/xmax)
    do i = ii, jj
        j = i - ii
        f0 = real(j, r8)/real(jj - ii, r8)
        f0 = (1.0_r8 + cos(f0*lo_pi))*0.5_r8
        y(i) = f0
    end do
    do i = jj, n
        y(i) = 0.0_r8
    end do
end subroutine

!> helper function to get the greens-function guy
elemental function diagonal_greensfunction(omega, Z, dl, gm) result(g)
    !> harmonic energy
    real(r8), intent(in) :: omega
    !> probing energy
    real(r8), intent(in) :: Z
    !> real part of self-energy
    real(r8), intent(in) :: dl
    !> imaginary part of self-energy
    real(r8), intent(in) :: gm
    !> greenfunction guy (but not really, just the common prefactor)
    complex(r8) :: g

    complex(r8) :: f0, f1

    f0 = 2*omega
    f1 = omega**2 + 2*omega*(dl - lo_imag*gm) - Z**2
    g = f0/f1
end function

!> evaluate the phonon spectral function
subroutine evaluate_spectral_function(bigOmega, sigmaIm, sigmaRe, omega, spectralfunction)
    !> probing frequency, defined from 0 to some high frequency
    real(r8), dimension(:), intent(in) :: bigOmega
    !> imaginary part of self-energy, defined for all values in bigOmega
    real(r8), dimension(:), intent(in) :: sigmaIm
    !> real part of self-energy, defined for all values in bigOmega
    real(r8), dimension(:), intent(in) :: sigmaRe
    !> harmonic frequency
    real(r8) :: omega
    !> intensity
    real(r8), dimension(:), intent(out) :: spectralfunction

    real(r8) :: f0, f1, f2
    integer :: i

    ! calculate the intensity
    do i = 1, size(bigOmega)
        f0 = 2*omega*sigmaIm(i)*2*omega/lo_pi
        f1 = (bigOmega(i)**2 - omega**2 - 2*omega*sigmaRe(i))**2
        f2 = 4*(omega**2)*(sigmaIm(i)**2)
        spectralfunction(i) = f0/(f1 + f2)
    end do
end subroutine

!> smear the spectral function with a Gaussian
subroutine gaussian_smear_distribution(x, sigma, y)
    !> energy axis
    real(r8), dimension(:), intent(in) :: x
    !> smearing
    real(r8), intent(in) :: sigma
    !> spectral function
    real(r8), dimension(:), intent(inout) :: y

    integer, parameter :: nsigma = 6
    real(r8), dimension(:), allocatable :: bufy, bufg
    real(r8) :: f0
    integer :: n, nstep, i, j, k

    ! First some heuristics when it comes to the smearing. I want to walk 6 sigma
    ! when doing the convolution -- how many steps in the frequency axis is that?
    n = size(x)
    f0 = x(2) - x(1) ! size of step
    nstep = ceiling(nsigma*sigma/f0)
    ! Make space for buffers, now defined for the full range in question?
    allocate (bufy(-(nstep - 1):n + nstep))
    allocate (bufg(-nstep:nstep))
    bufy = 0.0_r8
    bufg = 0.0_r8
    ! copy the intensity, unsmeared
    do i = 1, n
        bufy(i) = y(i)
    end do
    ! get the Gaussian thing for the convolution
    do i = -nstep, nstep
        f0 = i*(x(2) - x(1))
        bufg(i) = lo_gauss(f0, 0.0_r8, sigma)
    end do
    bufg = bufg/sum(bufg)

    ! Do the actual convolution. If I did it right, it should be zero at zero.
    y = 0.0_r8
    do i = 1, n
        do j = -nstep, nstep
            k = i + j
            y(i) = y(i) + bufg(j)*bufy(k)
        end do
    end do
end subroutine

!> smear the spectral function with a Gaussian
subroutine gaussian_smear_spectral_function(bigOmega, sigma, spectralfunction)
    !> energy axis
    real(r8), dimension(:), intent(in) :: bigOmega
    !> smearing
    real(r8), intent(in) :: sigma
    !> spectral function
    real(r8), dimension(:), intent(inout) :: spectralfunction

    integer, parameter :: nsigma = 6
    real(r8), dimension(:), allocatable :: bufy, bufg
    real(r8) :: f0
    integer :: n, nstep, i, j, k

    ! First some heuristics when it comes to the smearing. I want to walk 6 sigma
    ! when doing the convolution -- how many steps in the frequency axis is that?
    n = size(bigOmega)
    f0 = bigOmega(2) - bigOmega(1) ! size of step
    nstep = ceiling(nsigma*sigma/f0)
    ! Make space for buffers, now defined for the full range in question?
    allocate (bufy(-nstep:n + nstep - 1))
    allocate (bufg(-nstep:nstep))
    bufy = 0.0_r8
    bufg = 0.0_r8
    ! copy the intensity, unsmeared
    do i = 1, n - 1
        bufy(i) = spectralfunction(i + 1)
    end do
    ! add the extra tails, remember that the function is odd.
    do i = 1, nstep
        bufy(-i) = -bufy(i)
    end do
    ! get the Gaussian thing for the convolution
    do i = -nstep, nstep
        f0 = i*(bigOmega(2) - bigOmega(1))
        bufg(i) = lo_gauss(f0, 0.0_r8, sigma)
    end do
    bufg = bufg/sum(bufg)

    ! Do the actual convolution. If I did it right, it should be zero at zero.
    spectralfunction = 0.0_r8
    do i = 1, n
        do j = -nstep, nstep
            k = i + j - 1
            spectralfunction(i) = spectralfunction(i) + bufg(j)*bufy(k)
        end do
    end do
    spectralfunction(1) = 0.0_r8
end subroutine

!> locate a point on the FFT grid
function index_on_grid(qp, r) result(ind)
    !> q-mesh
    class(lo_qpoint_mesh), intent(in) :: qp
    !> coordinate (fractional)
    real(r8), dimension(3), intent(in) :: r
    !> linear index to full q-point
    integer :: ind

    integer, dimension(3) :: gi

    select type (qp)
    type is (lo_fft_mesh)
        if (qp%is_point_on_grid(r)) then
            gi = qp%index_from_coordinate(r)
            ind = qp%gridind2ind(gi(1), gi(2), gi(3))
        else
            call lo_stop_gracefully(['Could not locate q-point in mesh.'], lo_exitcode_param, __FILE__, __LINE__)
        end if
    class default
        call lo_stop_gracefully(['Can only locate points on FFT meshes'], lo_exitcode_param, __FILE__, __LINE__)
    end select
end function

!> find peak maximum of spectral function
subroutine find_spectral_function_max_and_fwhm(omega, maxomega, energy, sigmaIm, sigmaRe, energy_max, energy_lo, energy_hi)
    !> harmonic frequency
    real(r8), intent(in) :: omega
    !> largest harmonic frequency
    real(r8), intent(in) :: maxomega
    !> energy axis
    real(r8), dimension(:), intent(in) :: energy
    !> imaginary self-energy
    real(r8), dimension(:), intent(in) :: sigmaIm
    !> real self-energy
    real(r8), dimension(:), intent(in) :: sigmaRe
    !> highest peak
    real(r8), intent(out) :: energy_max
    !> fwhm lower
    real(r8), intent(out) :: energy_lo
    !> fwhm higher
    real(r8), intent(out) :: energy_hi

    type(lo_brent_helper) :: bt
    integer, parameter :: maxiter = 200
    integer, parameter :: derivorder = 3
    real(r8), dimension(derivorder*2 + 1, 4) :: sc
    real(r8), dimension(derivorder*2 + 1) :: fval
    real(r8), dimension(:), allocatable :: y0, taperIm
    real(r8) :: x0, xlo, xhi, ylo, yhi, derivdl, peak_max, f0
    integer :: n, i, iter, ind_max, ctr, ind_left, ind_right

    ! First evaluate the spectral function. Make sure to
    ! add a small amount to the imaginary part so that
    ! we are guaranteed to get a resolveable peak.
    n = size(energy)
    allocate (y0(n))
    allocate (taperIm(n))
    y0 = 0.0_r8
    f0 = energy(2) - energy(1)
    y0 = max(sigmaIm, f0)

    call taperfn_im(energy, maxomega, omega, taperIm)
    taperIm = y0*taperIm
    call evaluate_spectral_function(energy, taperIm, sigmaRe, omega, y0)

    ! Rough search to get a starting point
    x0 = -1.0_r8
    f0 = -lo_huge
    ind_max = -1
    do i = 2, n
        if (y0(i) .gt. f0) then
            f0 = y0(i)
            x0 = energy(i)
            ind_max = i
        end if
    end do

    ! Rough search to get fwhm locations to the left and to the right?
    peak_max = y0(ind_max)

    f0 = lo_huge
    ind_left = -1
    do i = 1, ind_max - 1
        if (abs(peak_max*0.5_r8 - y0(i)) .lt. f0) then
            f0 = abs(peak_max*0.5_r8 - y0(i))
            ind_left = i
        end if
    end do

    f0 = lo_huge
    ind_right = -1
    do i = ind_max + 1, n
        if (abs(peak_max*0.5_r8 - y0(i)) .lt. f0) then
            f0 = abs(peak_max*0.5_r8 - y0(i))
            ind_right = i
        end if
    end do

    ! Values that bracket the maximum, to be on the safe side.
    xlo = energy(ind_max - 1)
    xhi = energy(ind_max + 1)
    derivdl = (energy(2) - energy(1))*1E-2_r8

    ! Try a few times to see if we get something that sticks?
    initloop: do iter = 1, 10
        ! Function values (i.e. derivatives)
        call lo_centraldifference(derivorder, xlo, derivdl, sc)
        do i = 1, 2*derivorder + 1
            fval(i) = interpolated_spectral_function(energy, sigmaIm, sigmaRe, omega, 1.0_r8, sc(i, 1))
        end do
        ylo = sum(fval*sc(:, 2))

        call lo_centraldifference(derivorder, xhi, derivdl, sc)
        do i = 1, 2*derivorder + 1
            fval(i) = interpolated_spectral_function(energy, sigmaIm, sigmaRe, omega, 1.0_r8, sc(i, 1))
        end do
        yhi = sum(fval*sc(:, 2))

        ! If we bracket, we are good!
        if (ylo*yhi .lt. 0.0_r8) then
            exit initloop
        end if

        ! Ok, did not bracket. Try something.
        if (yhi .gt. 0.0_r8) then
            ! Move xhi a little closer?
            xhi = 0.5_r8*xhi + x0*0.5_r8
        end if
        if (ylo .lt. 0.0_r8) then
            ! Move xlo a little closer?
            xlo = 0.5_r8*xlo + x0*0.5_r8
        end if

        ! Give up after a few attempts. This happens when the peak is at zero, I think.
        ! This is a somewhat safe choice anyway.
        if (iter .ge. 10) then
            ! I will just use the rough values instead. Close enough.
            energy_max = energy(ind_max)
            if (ind_left .eq. 1) then
                energy_lo = energy_max*0.5_r8
            else
                energy_lo = energy(ind_left)
            end if
            if (ind_right .eq. n) then
                energy_hi = (energy_max + energy(n))*0.5_r8
            else
                energy_hi = energy(ind_right)
            end if
            ! return early!
            return
        end if
    end do initloop

    call bt%init(xlo, xhi, ylo, yhi, (energy(2) - energy(1))*1E-7_r8)
    x0 = -1.0_r8
    energy_max = -1.0_r8
    peak_max = -1.0_r8
    ctr = 0
    do iter = 1, maxiter
        ! Check for convergence
        if (abs(x0 - bt%next_x()) .lt. bt%tolerance) ctr = ctr + 1
        if (ctr .ge. 3) then
            ! found a zero
            energy_max = bt%current_x()
            peak_max = fval(derivorder + 1)
            exit
        end if
        ! pick next x-value
        x0 = bt%next_x()
        ! evaluate function
        call lo_centraldifference(derivorder, x0, derivdl, sc)
        do i = 1, 2*derivorder + 1
            fval(i) = interpolated_spectral_function(energy, sigmaIm, sigmaRe, omega, 1.0_r8, sc(i, 1))
        end do
        f0 = sum(fval*sc(:, 2))
        ! tell the solver what happened
        call bt%update(x0, f0)
    end do

    ! Now do the same thing again, but look for the FWHM to the left.
    energy_lo = -lo_huge
    xhi = energy_max
    yhi = peak_max - peak_max*0.5_r8
    xlo = -1
    ctr = -1
    do i = n, 1, -1
        if (energy(i) .gt. energy_max - (energy(2) - energy(1))) cycle
        if (y0(i) .lt. peak_max*0.45_r8) then
            xlo = energy(i)
            exit
        end if
    end do
    if (xlo .lt. 0.0_r8) xlo = 0.0_r8
    ylo = interpolated_spectral_function(energy, sigmaIm, sigmaRe, omega, 1.0_r8, xlo)
    ylo = ylo - peak_max*0.5_r8

    ! Check if something went wrong in the bracketing. Could happen with weird lineshapes.
    if (ylo*yhi .gt. 0.0_r8) then
        energy_max = energy(ind_max)
        if (ind_left .eq. 1) then
            energy_lo = energy_max*0.5_r8
        elseif (energy(ind_left) .ge. energy_max) then
            energy_lo = energy_max*0.5_r8
        else
            energy_lo = energy(ind_left)
        end if
    end if

    if (energy_lo .lt. 0.0_r8) then
        call bt%init(xlo, xhi, ylo, yhi, (energy(2) - energy(1))*1E-9_r8)
        x0 = -1.0_r8
        ctr = 0
        do iter = 1, maxiter
            ! Check for convergence
            if (abs(x0 - bt%next_x()) .lt. bt%tolerance) ctr = ctr + 1
            if (ctr .ge. 3) then
                ! found a zero
                energy_lo = bt%current_x()
                exit
            end if
            ! pick next x-value
            x0 = bt%next_x()
            ! evaluate function
            f0 = interpolated_spectral_function(energy, sigmaIm, sigmaRe, omega, 1.0_r8, x0)
            f0 = f0 - peak_max*0.5_r8
            ! tell the solver what happened
            call bt%update(x0, f0)
        end do
    end if

    ! And once again, in the other direction.
    energy_hi = -lo_huge
    xlo = energy_max
    ylo = peak_max - peak_max*0.5_r8
    xhi = -1
    do i = 1, n
        if (energy(i) .lt. energy_max + energy(2) - energy(1)) cycle
        if (y0(i) .lt. peak_max*0.45_r8) then
            xhi = energy(i)
            exit
        end if
    end do
    if (xhi .lt. 0.0_r8) xhi = energy(n - 1)
    yhi = interpolated_spectral_function(energy, sigmaIm, sigmaRe, omega, 1.0_r8, xhi)
    yhi = yhi - peak_max*0.5_r8
    if (ylo*yhi .gt. 0.0_r8) then
        ! This did not really work. Give the rough value instead!
        if (ind_right .eq. n) then
            energy_hi = (energy_max + energy(n))*0.5_r8
        else
            energy_hi = energy(ind_right)
        end if
    end if

    ! Only do this if bracketing worked
    if (energy_hi .lt. 0.0_r8) then
        ! Find the exact value
        call bt%init(xlo, xhi, ylo, yhi, (energy(2) - energy(1))*1E-9_r8)
        x0 = -1.0_r8
        ctr = 0
        do iter = 1, maxiter
            ! Check for convergence
            if (abs(x0 - bt%next_x()) .lt. bt%tolerance) ctr = ctr + 1
            if (ctr .ge. 3) then
                ! found a zero
                energy_hi = bt%current_x()
                exit
            end if
            ! pick next x-value
            x0 = bt%next_x()
            ! evaluate function
            f0 = interpolated_spectral_function(energy, sigmaIm, sigmaRe, omega, 1.0_r8, x0)
            f0 = f0 - peak_max*0.5_r8
            ! tell the solver what happened
            call bt%update(x0, f0)
        end do
    end if

    ! Cleanup
    deallocate (y0)
    deallocate (taperIm)
end subroutine

!> Integrate spectral function in a semi-smart way
subroutine integrate_spectral_function(x, omega, sigmaIm, sigmaRe, xmid, xlo, xhi, scalefactor, temperature, tolerance, integrated_spectralfunction, integrated_sf_squared, integrated_sf_nn)
    !> rough x-axis
    real(r8), dimension(:), intent(in) :: x
    !> harmonic frequency
    real(r8), intent(in) :: omega
    !> imaginary part
    real(r8), dimension(:), intent(in) :: sigmaIm
    !> real part
    real(r8), dimension(:), intent(in) :: sigmaRe
    !> peak location
    real(r8), intent(in) :: xmid
    !> fwhm left location
    real(r8), intent(in) :: xlo
    !> fwhm right location
    real(r8), intent(in) :: xhi
    !> scaling factor for the spectral function
    real(r8), intent(in) :: scalefactor
    !> temperature
    real(r8), intent(in) :: temperature
    !> relative tolerance to integrate to?
    real(r8), intent(in) :: tolerance
    !> integral of spectral function
    real(r8), intent(out) :: integrated_spectralfunction
    !> integrated pi*sf^2
    real(r8), intent(out) :: integrated_sf_squared
    !> integrated sf^2*n*(n+1)
    real(r8), intent(out) :: integrated_sf_nn

    integer, parameter :: maxiter = 50
    integer, parameter :: maxintervals = 10000
    ! Order of quadratures to check.
    integer, parameter :: nlo = 8
    integer, parameter :: nhi = 11
    real(r8), dimension(2, nlo) :: gqlo
    real(r8), dimension(2, nhi) :: gqhi
    real(r8), dimension(nlo) :: lo_sf, lo_sq, lo_kp
    real(r8), dimension(nhi) :: hi_sf, hi_sq, hi_kp

    real(r8), dimension(:), allocatable :: current_nodes, next_nodes
    real(r8), dimension(:), allocatable :: current_sf, current_sq, current_kp
    real(r8), dimension(:), allocatable :: next_sf, next_sq, next_kp
    logical, dimension(:), allocatable :: current_ok, next_ok
    logical :: segmentok
    real(r8) :: xmax, f0, f1, x0
    real(r8) :: int_sf, int_sq, int_kp
    real(r8) :: err_sf, err_sq, err_kp
    real(r8) :: vlo_sf, vlo_sq, vlo_kp
    real(r8) :: vhi_sf, vhi_sq, vhi_kp
    integer :: n, n_curr
    integer :: iter, i, j, ctr, ctrok, ctrnew

    n = size(x)
    xmax = x(n)

    ! Maybe try do do this recursively somehow.
    ! Think I will do it like this, I will chop things into a
    ! few intervals. Then do a nlo and nhi point Gaussian quadrature
    ! on these intervals. If within tolerance, fine, if not split
    ! it in half and go again? Can make this faster later, just
    ! proof of concept for now.

    ! Start building the first intervals.
    n_curr = 9
    allocate (current_nodes(n_curr))
    allocate (current_ok(n_curr))
    allocate (current_sf(n_curr))
    allocate (current_sq(n_curr))
    allocate (current_kp(n_curr))
    current_nodes = 0.0_r8
    current_ok = .false.
    current_sf = 0.0_r8
    current_sq = 0.0_r8
    current_kp = 0.0_r8

    ! Select some sensible starting nodes.
    current_nodes(1) = 0.0_r8
    current_nodes(2) = xlo*0.5_r8
    current_nodes(3) = xlo
    current_nodes(4) = xlo + (xmid - xlo)*0.5_r8
    current_nodes(5) = xmid
    current_nodes(6) = xmid + (xhi - xmid)*0.5_r8
    current_nodes(7) = xhi
    current_nodes(8) = xhi + (xmax - xhi)*0.5_r8
    current_nodes(9) = xmax

    ! Start adaptively integrating
    iterloop: do iter = 1, maxiter
        allocate (next_nodes(2*n_curr))
        allocate (next_ok(2*n_curr))
        allocate (next_sf(2*n_curr))
        allocate (next_sq(2*n_curr))
        allocate (next_kp(2*n_curr))
        next_nodes = -1.0_r8
        next_ok = .false.
        next_sf = 0.0_r8
        next_sq = 0.0_r8
        next_kp = 0.0_r8

        err_sf = 0.0_r8
        err_sq = 0.0_r8
        err_kp = 0.0_r8
        int_sf = 0.0_r8
        int_sq = 0.0_r8
        int_kp = 0.0_r8

        ctrok = 0
        ctrnew = 0
        ctr = 0
        do i = 1, n_curr - 1

            ! Check interval.
            if (current_ok(i)) then
                ! This interval is already fine.
                segmentok = .true.
                int_sf = int_sf + current_sf(i)
                int_sq = int_sq + current_sq(i)
                int_kp = int_kp + current_kp(i)
            else
                ! We have to actually check. Evaluate functions and integrate.
                call lo_gaussianquadrature(nlo, current_nodes(i), current_nodes(i + 1), gqlo)
                call lo_gaussianquadrature(nhi, current_nodes(i), current_nodes(i + 1), gqhi)
                do j = 1, nlo
                    x0 = gqlo(1, j)
                    f0 = interpolated_spectral_function(x, sigmaIm, sigmaRe, omega, scalefactor, x0)
                    f1 = lo_planck(temperature, x0)
                    lo_sf(j) = f0
                    lo_sq(j) = f0*f0*lo_pi
                    lo_kp(j) = f0*f0*f1*(f1 + 1.0_r8)
                end do
                do j = 1, nhi
                    x0 = gqhi(1, j)
                    f0 = interpolated_spectral_function(x, sigmaIm, sigmaRe, omega, scalefactor, x0)
                    f1 = lo_planck(temperature, x0)
                    hi_sf(j) = f0
                    hi_sq(j) = f0*f0*lo_pi
                    hi_kp(j) = f0*f0*f1*(f1 + 1.0_r8)*x0*x0
                end do
                vlo_sf = sum(lo_sf*gqlo(2, :))
                vlo_sq = sum(lo_sq*gqlo(2, :))
                vlo_kp = sum(lo_kp*gqlo(2, :))
                vhi_sf = sum(hi_sf*gqhi(2, :))
                vhi_sq = sum(hi_sq*gqhi(2, :))
                vhi_kp = sum(hi_kp*gqhi(2, :))

                ! Check if fine:
                !if ( abs(f0-f1) .lt. tolerance ) then
                if (abs(vlo_sf - vhi_sf) .lt. tolerance) then
                    current_sf(i) = vhi_sf
                    current_sq(i) = vhi_sq
                    current_kp(i) = vhi_kp
                    segmentok = .true.
                else
                    segmentok = .false.
                end if
                ! Accumulate things
                int_sf = int_sf + vhi_sf
                int_sq = int_sq + vhi_sq
                int_kp = int_kp + vhi_kp

                err_sf = err_sf + abs(vhi_sf - vlo_sf)
                err_sq = err_sq + abs(vhi_sq - vlo_sq)
                err_kp = err_kp + abs(vhi_kp - vlo_kp)
            end if

            if (segmentok) then
                ! Just copy
                ctr = ctr + 1
                next_nodes(ctr) = current_nodes(i)
                next_ok(ctr) = .true.
                next_sf(ctr) = current_sf(i)
                next_sq(ctr) = current_sq(i)
                next_kp(ctr) = current_kp(i)
                ctrok = ctrok + 1
            else
                ! Split in half? Seems sensible.
                ctr = ctr + 1
                next_nodes(ctr) = current_nodes(i)
                next_ok(ctr) = .false.
                next_sf(ctr) = 0.0_r8
                next_sq(ctr) = 0.0_r8
                next_kp(ctr) = 0.0_r8
                ctr = ctr + 1
                next_nodes(ctr) = (current_nodes(i) + current_nodes(i + 1))*0.5_r8
                next_ok(ctr) = .false.
                next_sf(ctr) = 0.0_r8
                next_sq(ctr) = 0.0_r8
                next_kp(ctr) = 0.0_r8
                ! make a note I added a segment
                ctrnew = ctrnew + 1
            end if
        end do

        !write(*,*) iter,'int',int_sf,int_kp,err_sf,ctr,n_curr,ctrnew

        ! Check if we are converged
        if (ctr + 1 .eq. n_curr) exit iterloop

        deallocate (current_nodes)
        deallocate (current_ok)
        deallocate (current_sf)
        deallocate (current_sq)
        deallocate (current_kp)
        n_curr = ctr + 1
        allocate (current_nodes(n_curr))
        allocate (current_ok(n_curr))
        allocate (current_sf(n_curr))
        allocate (current_sq(n_curr))
        allocate (current_kp(n_curr))
        current_nodes = 0.0_r8
        current_ok = .false.
        current_sf = 0.0_r8
        current_sq = 0.0_r8
        current_kp = 0.0_r8
        current_nodes(1:ctr) = next_nodes(1:ctr)
        current_nodes(ctr + 1) = xmax
        current_ok(1:ctr) = next_ok(1:ctr)
        current_sf(1:ctr) = next_sf(1:ctr)
        current_sq(1:ctr) = next_sq(1:ctr)
        current_kp(1:ctr) = next_kp(1:ctr)
        deallocate (next_nodes)
        deallocate (next_ok)
        deallocate (next_sf)
        deallocate (next_sq)
        deallocate (next_kp)
    end do iterloop

    integrated_spectralfunction = int_sf
    integrated_sf_squared = int_sq
    integrated_sf_nn = int_kp

    ! A little bit of cleanup?
    !deallocate(current_nodes)
    !deallocate(current_ok)
    !deallocate(current_sf)
    !deallocate(current_sq)
    !deallocate(current_kp)
    !deallocate(next_nodes)
    !deallocate(next_ok)
    !deallocate(next_sf)
    !deallocate(next_sq)
    !deallocate(next_kp)
end subroutine

!> Integrate two spectral functions in a semi-smart way
subroutine integrate_two_spectral_functions(x, omega1, omega2, sigmaIm1, sigmaIm2, sigmaRe1, sigmaRe2, xmid1, xmid2, xlo1, xlo2, xhi1, xhi2, scalefactor1, scalefactor2, temperature, tolerance, integrated_sf_squared, integrated_sf_nn)
    !> rough x-axis
    real(r8), dimension(:), intent(in) :: x
    !> harmonic frequency
    real(r8), intent(in) :: omega1, omega2
    !> imaginary part
    real(r8), dimension(:), intent(in) :: sigmaIm1, sigmaIm2
    !> real part
    real(r8), dimension(:), intent(in) :: sigmaRe1, sigmaRe2
    !> peak location
    real(r8), intent(in) :: xmid1, xmid2
    !> fwhm left location
    real(r8), intent(in) :: xlo1, xlo2
    !> fwhm right location
    real(r8), intent(in) :: xhi1, xhi2
    !> scaling factor for the spectral function
    real(r8), intent(in) :: scalefactor1, scalefactor2
    !> temperature
    real(r8), intent(in) :: temperature
    !> relative tolerance to integrate to?
    real(r8), intent(in) :: tolerance
    !> integrated pi*sf^2
    real(r8), intent(out) :: integrated_sf_squared
    !> integrated sf^2*n*(n+1)
    real(r8), intent(out) :: integrated_sf_nn

    integer, parameter :: maxiter = 40
    ! Order of quadratures to check.
    integer, parameter :: nlo = 8
    integer, parameter :: nhi = 11
    real(r8), dimension(2, nlo) :: gqlo
    real(r8), dimension(2, nhi) :: gqhi
    real(r8), dimension(nlo) :: lo_sq, lo_kp
    real(r8), dimension(nhi) :: hi_sq, hi_kp

    real(r8), dimension(:), allocatable :: current_nodes, next_nodes
    real(r8), dimension(:), allocatable :: current_sq, current_kp
    real(r8), dimension(:), allocatable :: next_sq, next_kp
    logical, dimension(:), allocatable :: current_ok, next_ok
    logical :: segmentok
    real(r8) :: xmax, f0, f1, f2, x0
    real(r8) :: int_sq, int_kp
    real(r8) :: err_sq, err_kp
    real(r8) :: vlo_sq, vlo_kp
    real(r8) :: vhi_sq, vhi_kp
    integer :: n, n_curr
    integer :: iter, i, j, ctr, ctrok, ctrnew

    n = size(x)
    xmax = x(n)

    ! Maybe try do do this recursively somehow.
    ! Think I will do it like this, I will chop things into a
    ! few intervals. Then do a 10 and 15 point Gaussian quadrature
    ! on these intervals. If within tolerance, fine, if not split
    ! it in half and go again? Can make this faster later, just
    ! proof

    ! Select some sensible starting nodes.
    allocate (next_nodes(16))
    next_nodes(1) = 0.0_r8
    next_nodes(2) = xlo1*0.5_r8
    next_nodes(3) = xlo1
    next_nodes(4) = xlo1 + (xmid1 - xlo1)*0.5_r8
    next_nodes(5) = xmid1
    next_nodes(6) = xmid1 + (xhi1 - xmid1)*0.5_r8
    next_nodes(7) = xhi1
    next_nodes(8) = xhi1 + (xmax - xhi1)*0.5_r8
    next_nodes(9) = xmax
    next_nodes(10) = xlo2*0.5_r8
    next_nodes(11) = xlo2
    next_nodes(12) = xlo2 + (xmid2 - xlo2)*0.5_r8
    next_nodes(13) = xmid2
    next_nodes(14) = xmid2 + (xhi2 - xmid2)*0.5_r8
    next_nodes(15) = xhi2
    next_nodes(16) = xhi2 + (xmax - xhi2)*0.5_r8
    call lo_return_unique(next_nodes, current_nodes, xmax/1000.0_r8)
    deallocate (next_nodes)
    call lo_qsort(current_nodes)
    n_curr = size(current_nodes)
    current_nodes(1) = 0.0_r8
    current_nodes(n_curr) = xmax

    allocate (current_ok(n_curr))
    allocate (current_sq(n_curr))
    allocate (current_kp(n_curr))
    current_ok = .false.
    current_sq = 0.0_r8
    current_kp = 0.0_r8

    ! Start adaptively integrating
    iterloop: do iter = 1, maxiter
        allocate (next_nodes(2*n_curr))
        allocate (next_ok(2*n_curr))
        allocate (next_sq(2*n_curr))
        allocate (next_kp(2*n_curr))
        next_nodes = -1.0_r8
        next_ok = .false.
        next_sq = 0.0_r8
        next_kp = 0.0_r8

        err_sq = 0.0_r8
        err_kp = 0.0_r8
        int_sq = 0.0_r8
        int_kp = 0.0_r8

        ctrok = 0
        ctrnew = 0
        ctr = 0
        do i = 1, n_curr - 1

            ! Check interval.
            if (current_ok(i)) then
                ! This interval is already fine.
                segmentok = .true.
                int_sq = int_sq + current_sq(i)
                int_kp = int_kp + current_kp(i)
            else
                ! We have to actually check. Evaluate functions and integrate.
                call lo_gaussianquadrature(nlo, current_nodes(i), current_nodes(i + 1), gqlo)
                call lo_gaussianquadrature(nhi, current_nodes(i), current_nodes(i + 1), gqhi)
                do j = 1, nlo
                    x0 = gqlo(1, j)
                    f0 = interpolated_spectral_function(x, sigmaIm1, sigmaRe1, omega1, scalefactor1, x0)
                    f1 = interpolated_spectral_function(x, sigmaIm2, sigmaRe2, omega2, scalefactor2, x0)
                    f2 = lo_planck(temperature, x0)
                    lo_sq(j) = f0*f1*lo_pi
                    lo_kp(j) = f0*f1*f2*(f2 + 1.0_r8)*x0*x0
                end do
                do j = 1, nhi
                    x0 = gqhi(1, j)
                    f0 = interpolated_spectral_function(x, sigmaIm1, sigmaRe1, omega1, scalefactor1, x0)
                    f1 = interpolated_spectral_function(x, sigmaIm2, sigmaRe2, omega2, scalefactor2, x0)
                    f2 = lo_planck(temperature, x0)
                    hi_sq(j) = f0*f1*lo_pi
                    hi_kp(j) = f0*f1*f2*(f2 + 1.0_r8)*x0*x0
                end do
                vlo_sq = sum(lo_sq*gqlo(2, :))
                vlo_kp = sum(lo_kp*gqlo(2, :))
                vhi_sq = sum(hi_sq*gqhi(2, :))
                vhi_kp = sum(hi_kp*gqhi(2, :))

                ! Check if fine:
                if (abs(vlo_sq - vhi_sq) .lt. tolerance) then
                    current_sq(i) = vhi_sq
                    current_kp(i) = vhi_kp
                    segmentok = .true.
                else
                    segmentok = .false.
                end if
                ! Accumulate things
                int_sq = int_sq + vhi_sq
                int_kp = int_kp + vhi_kp

                err_sq = err_sq + abs(vhi_sq - vlo_sq)
                err_kp = err_kp + abs(vhi_kp - vlo_kp)
            end if

            if (segmentok) then
                ! Just copy
                ctr = ctr + 1
                next_nodes(ctr) = current_nodes(i)
                next_ok(ctr) = .true.
                next_sq(ctr) = current_sq(i)
                next_kp(ctr) = current_kp(i)
                ctrok = ctrok + 1
            else
                ! Split in half? Seems sensible.
                ctr = ctr + 1
                next_nodes(ctr) = current_nodes(i)
                next_ok(ctr) = .false.
                next_sq(ctr) = 0.0_r8
                next_kp(ctr) = 0.0_r8
                ctr = ctr + 1
                next_nodes(ctr) = (current_nodes(i) + current_nodes(i + 1))*0.5_r8
                next_ok(ctr) = .false.
                next_sq(ctr) = 0.0_r8
                next_kp(ctr) = 0.0_r8
                ! make a note I added a segment
                ctrnew = ctrnew + 1
            end if
        end do

        !write(*,*) iter,'int',int_sq,int_kp,err_sq,ctr,n_curr,ctrnew

        ! Check if we are converged
        if (ctr + 1 .eq. n_curr) exit iterloop

        deallocate (current_nodes)
        deallocate (current_ok)
        deallocate (current_sq)
        deallocate (current_kp)
        n_curr = ctr + 1
        allocate (current_nodes(n_curr))
        allocate (current_ok(n_curr))
        allocate (current_sq(n_curr))
        allocate (current_kp(n_curr))
        current_nodes = 0.0_r8
        current_ok = .false.
        current_sq = 0.0_r8
        current_kp = 0.0_r8
        current_nodes(1:ctr) = next_nodes(1:ctr)
        current_nodes(ctr + 1) = xmax
        current_ok(1:ctr) = next_ok(1:ctr)
        current_sq(1:ctr) = next_sq(1:ctr)
        current_kp(1:ctr) = next_kp(1:ctr)
        deallocate (next_nodes)
        deallocate (next_ok)
        deallocate (next_sq)
        deallocate (next_kp)
    end do iterloop

    integrated_sf_squared = int_sq
    integrated_sf_nn = int_kp

    ! A little bit of cleanup?
    if (allocated(current_nodes)) deallocate (current_nodes)
    if (allocated(current_ok)) deallocate (current_ok)
    if (allocated(current_sq)) deallocate (current_sq)
    if (allocated(current_kp)) deallocate (current_kp)
    if (allocated(next_nodes)) deallocate (next_nodes)
    if (allocated(next_ok)) deallocate (next_ok)
    if (allocated(next_sq)) deallocate (next_sq)
    if (allocated(next_kp)) deallocate (next_kp)
end subroutine

!> size in memory, in bytes
function isf_size_in_mem(isf) result(mem)
    !> atom
    class(lo_spectralfunction_helper), intent(in) :: isf
    !> size in memory, bytes
    integer(i8) :: mem

    mem = 0
    mem = mem + storage_size(isf)

    if (allocated(isf%qv)) mem = mem + storage_size(isf%qv)*size(isf%qv)
    if (allocated(isf%spectralfunction)) mem = mem + storage_size(isf%spectralfunction)*size(isf%spectralfunction)
    mem = mem/8
end function

end module
