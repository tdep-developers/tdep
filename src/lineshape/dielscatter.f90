module dielscatter
!!
!! Handles dielectric response, such as dielectric constants, infrared
!! and Raman scattering
!!
use konstanter, only: i8, r8, lo_huge, lo_iou, lo_freqtol, lo_tiny, lo_exitcode_symmetry, lo_imag, lo_twopi, lo_frequency_Hartree_to_THz, &
                      lo_pi, lo_frequency_hartree_to_Hz, lo_sqtol, lo_exitcode_param
use gottochblandat, only: tochar, walltime, lo_chop, lo_trueNtimes, lo_progressbar_init, lo_progressbar, lo_outerproduct, &
                          lo_gauss, lo_planck, lo_trapezoid_integration, lo_unflatten_2tensor, lo_unflatten_4tensor, lo_linear_interpolation, &
                          lo_points_on_sphere, lo_flattentensor
use mpi_wrappers, only: lo_mpi_helper, lo_stop_gracefully
use lo_memtracker, only: lo_mem_helper
use lo_timetracker, only: lo_timer
use type_crystalstructure, only: lo_crystalstructure
use type_forceconstant_secondorder, only: lo_forceconstant_secondorder
use type_forceconstant_thirdorder, only: lo_forceconstant_thirdorder
use lo_dielectric_interaction, only: lo_dielectric_tensor
use type_qpointmesh, only: lo_qpoint_mesh, lo_qpoint
use type_phonon_dispersions, only: lo_phonon_dispersions, lo_phonon_dispersions_qpoint
use type_symmetryoperation, only: lo_expandoperation_pair, lo_expandoperation_quartet
use hdf5_wrappers, only: lo_hdf5_helper
use fftw_wrappers, only: lo_abcd_convolution

use options, only: lo_opts
use phonondamping, only: lo_phonon_selfenergy
use lineshape_helper, only: diagonal_greensfunction, lo_spectralfunction_helper, lo_convolution_helper
implicit none

private

public :: lo_dielectric_response

type lo_dielectric_response
    !> order 0 dipole response (xyz,xyz,energy,mode)
    complex(r8), dimension(:, :, :, :), allocatable :: chi_tensor_zero
    !> order 1 dipole response (xyz,xyz,energy,mode)
    complex(r8), dimension(:, :, :, :), allocatable :: chi_tensor_one, chi_tensor_one_sph
    !> order 2 dipole response (xyz,xyz,energy)
    complex(r8), dimension(:, :, :), allocatable :: chi_tensor_two
    !> order 3 dipole response (xyz,xyz,energy)
    complex(r8), dimension(:, :, :), allocatable :: chi_tensor_three

    !> order 0 raman response (xyz^4,energy,mode)
    real(r8), dimension(:, :, :), allocatable :: I_tensor_zero
    !> order 1 raman response (xyz^4,energy,mode)
    real(r8), dimension(:, :, :), allocatable :: I_tensor_one
    !> order 2 raman response (xyz^4,energy,mode,mode)
    real(r8), dimension(:, :, :, :), allocatable :: I_tensor_two

    !> baseline dielectric tensor
    real(r8), dimension(3, 3) :: eps_inf = -lo_huge
    !> dynamic correction to dielectric tensor
    real(r8), dimension(3, 3) :: eps_inf_corr = -lo_huge
contains
    procedure :: generate
    procedure :: size_in_mem => dir_size_in_mem
end type

type diel_helper
    !> compound matrix element raman I (xyz^4,s1)
    real(r8), dimension(:, :), allocatable :: cmp_rmI
    !> compound matrix element raman II (xyz^4,s1,s2,q)
    real(r8), dimension(:, :, :, :), allocatable :: cmp_rmII

    !> compound matrix element IR I (xyz^2,s1)
    real(r8), dimension(:, :), allocatable :: cmp_irI
    !> compound matrix element IR I (xyz^2,s1), spherically averaged
    real(r8), dimension(:, :), allocatable :: cmp_irI_sph
    !> compound matrix element IR II (xyz^2,s1,s2,q)
    real(r8), dimension(:, :, :, :), allocatable :: cmp_irII
    !> compound matrix element IR III (xyz^2,s1,s2,s3,q)
    real(r8), dimension(:, :, :, :, :), allocatable :: cmp_irIII

    ! pre-transformed eigenvectors
    complex(r8), dimension(:, :), allocatable :: ugv_gamma
    complex(r8), dimension(:, :, :), allocatable :: ugv

    ! shifts and linewidhts
    real(r8), dimension(:, :), allocatable :: shift
    real(r8), dimension(:, :), allocatable :: linewidth
    real(r8), dimension(:, :, :), allocatable :: spectral_function
    real(r8), dimension(:), allocatable :: omega_axis
contains
    !> create the helper
    procedure :: generate => buildhelper
    procedure :: size_in_mem => dh_size_in_mem
end type diel_helper

interface
    module subroutine buildhelper(dh, wp, di, p, qp, dr, fc, fct, temperature, deltaeps, mw, mem, verbosity)
        class(diel_helper), intent(out) :: dh
        type(lo_phonon_dispersions_qpoint), intent(in) :: wp
        type(lo_dielectric_tensor), intent(in) :: di
        type(lo_crystalstructure), intent(in) :: p
        class(lo_qpoint_mesh), intent(in) :: qp
        type(lo_phonon_dispersions), intent(in) :: dr
        type(lo_forceconstant_secondorder), intent(inout) :: fc
        type(lo_forceconstant_thirdorder), intent(in) :: fct
        real(r8), intent(in) :: temperature
        real(r8), dimension(:, :), intent(out) :: deltaeps
        type(lo_mpi_helper), intent(inout) :: mw
        type(lo_mem_helper), intent(inout) :: mem
        integer, intent(in) :: verbosity
    end subroutine
end interface

contains

!> create all the spectra
subroutine generate(dir, wp, di, p, qp, dr, fc, fct, se, isf, opts, tmr, mw, mem, verbosity)
    !> dielectric response
    class(lo_dielectric_response), intent(out) :: dir
    !> harmonic properties at this q-point (which should be Gamma)
    type(lo_phonon_dispersions_qpoint), intent(in) :: wp
    !> interaction tensors
    type(lo_dielectric_tensor), intent(in) :: di
    !> structure
    type(lo_crystalstructure), intent(in) :: p
    !> q-point mesh
    class(lo_qpoint_mesh), intent(in), allocatable :: qp
    !> phonon dispersions
    type(lo_phonon_dispersions), intent(in) :: dr
    !> forceconstant
    type(lo_forceconstant_secondorder), intent(inout) :: fc
    !> thirdorder forceconstant
    type(lo_forceconstant_thirdorder), intent(in) :: fct
    !> phonon self-energy
    type(lo_phonon_selfenergy), intent(inout) :: se
    !> tabulated spectral functions
    type(lo_spectralfunction_helper), intent(in) :: isf
    !> options
    type(lo_opts), intent(in) :: opts
    !> timer (start and stop outside this routine)
    type(lo_timer), intent(inout) :: tmr
    !> mpi helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem
    !> talk a lot?
    integer, intent(in) :: verbosity

    type(diel_helper) :: dh
    real(r8) :: timer, t0, t1

    call mem%tick()
    ! start timers
    timer = walltime()
    t0 = timer
    t1 = timer

    ! set some general things
    init: block

        if (verbosity .gt. 0) then
            write (lo_iou, *) ''
            write (lo_iou, *) 'CALCULATING DIELECTRIC RESPONSE'
        end if

        ! Sanity test for compatible integration types
        select case (opts%integrationtype)
        case (3)
            call lo_stop_gracefully(['No tetrahedron integration for dielectric response'], lo_exitcode_param, __FILE__, __LINE__, mw%comm)
        case (4)
            call lo_stop_gracefully(['Deprecated for now.'], lo_exitcode_param, __FILE__, __LINE__, mw%comm)
        end select

        ! Start some sanity checks. First make sure we are at Gamma? Makes sense I suppose.
        ! I do this by checking that the eigenvectors are purely real, since that will make
        ! life significantly easier.
        if (sum(abs(aimag(wp%egv))) .gt. lo_tiny) then
            call lo_stop_gracefully(['Eigenvectors at Gamma not real, suspicious'], lo_exitcode_symmetry, __FILE__, __LINE__, mw%comm)
        end if

        ! Some space for the response
        allocate (dir%chi_tensor_zero(3, 3, se%n_energy, dr%n_mode))
        allocate (dir%chi_tensor_one(3, 3, se%n_energy, dr%n_mode))
        allocate (dir%chi_tensor_one_sph(3, 3, se%n_energy, dr%n_mode))
        allocate (dir%chi_tensor_two(3, 3, se%n_energy))
        allocate (dir%chi_tensor_three(3, 3, se%n_energy))
        dir%chi_tensor_zero = 0.0_r8
        dir%chi_tensor_one = 0.0_r8
        dir%chi_tensor_one_sph = 0.0_r8
        dir%chi_tensor_two = 0.0_r8
        dir%chi_tensor_three = 0.0_r8

        allocate (dir%I_tensor_zero(81, se%n_energy, dr%n_mode))
        allocate (dir%I_tensor_one(81, se%n_energy, dr%n_mode))
        allocate (dir%I_tensor_two(81, se%n_energy, dr%n_mode, dr%n_mode))
        dir%I_tensor_zero = 0.0_r8
        dir%I_tensor_one = 0.0_r8
        dir%I_tensor_two = 0.0_r8

        if (verbosity .gt. 0) then
            t1 = walltime()
            write (lo_iou, *) '... prepared space (', tochar(t1 - t0), 's)'
            t0 = t1
        end if

        call tmr%tock('init dielectric')
    end block init

    ! Create the helper with matrix elements and things
    call dh%generate(wp, di, p, qp, dr, fc, fct, opts%temperature, dir%eps_inf_corr, mw, mem, verbosity)
    dir%eps_inf = di%eps_inf

    if (verbosity .gt. 0) then
        write (lo_iou, *) ''
        write (lo_iou, *) 'current memory usage:   avg                  max                  min (MiB)'
    end if
    ! Dump current memory usage
    call memdump('structure', int(p%size_in_mem()), mw, lo_iou)
    call memdump('wp', int(wp%size_in_mem()), mw, lo_iou)
    call memdump('di', int(di%size_in_mem()), mw, lo_iou)
    call memdump('qp', int(qp%size_in_mem(qp)), mw, lo_iou)
    call memdump('dr', int(dr%size_in_mem()), mw, lo_iou)
    call memdump('fc', int(fc%size_in_mem()), mw, lo_iou)
    call memdump('fct', int(fct%size_in_mem()), mw, lo_iou)
    call memdump('se', int(se%size_in_mem()), mw, lo_iou)
    call memdump('isf', int(isf%size_in_mem()), mw, lo_iou)
    call memdump('dir', int(dir%size_in_mem()), mw, lo_iou)
    call memdump('dh', int(dh%size_in_mem()), mw, lo_iou)
    if (verbosity .gt. 0) then
        write (lo_iou, *) ''
    end if

    call tmr%tock('dielectric matrixelements')

    ! get the lowest order stuff
    term1: block
        complex(r8) :: c0
        real(r8), dimension(3, 3) :: mI, mIs
        real(r8) :: om, dl, gm, Z, smallsmearing
        integer :: ie, imode, ctr

        ! miniscule smearing to avoid divergencies?
        smallsmearing = (se%energy_axis(2) - se%energy_axis(1))*0.01_r8

        ctr = 0
        do imode = 1, dr%n_mode
            om = wp%omega(imode)
            if (om .lt. lo_freqtol) cycle

            mI = lo_unflatten_2tensor(dh%cmp_irI(:, imode))
            mIs = lo_unflatten_2tensor(dh%cmp_irI_sph(:, imode))
            do ie = 1, se%n_energy
                ctr = ctr + 1
                if (mod(ctr, mw%n) .ne. mw%r) cycle
                ! fetch probing energy and real and imaginary part of the self-energy
                Z = se%energy_axis(ie)
                dl = se%re_3ph(ie, imode) + se%re_4ph(ie, imode)
                gm = se%im_3ph(ie, imode) + se%im_iso(ie, imode)

                ! Zeroth order IR, selfenergy=0
                c0 = diagonal_greensfunction(om, Z, 0.0_r8, smallsmearing)
                dir%chi_tensor_zero(:, :, ie, imode) = c0*mI
                ! First order IR, with selfenergy
                c0 = diagonal_greensfunction(om, Z, dl, gm)
                dir%chi_tensor_one(:, :, ie, imode) = c0*mI
                ! First order IR, spherically averaged
                dir%chi_tensor_one_sph(:, :, ie, imode) = c0*mIs

                ! Zeroth order Raman, selfenergy=0
                c0 = diagonal_greensfunction(om, Z, 0.0_r8, smallsmearing)
                dir%I_tensor_zero(:, ie, imode) = aimag(c0)*dh%cmp_rmI(:, imode)
                ! First order Raman
                c0 = diagonal_greensfunction(om, Z, dl, gm)
                dir%I_tensor_one(:, ie, imode) = aimag(c0)*dh%cmp_rmI(:, imode)
            end do
        end do

        call mw%allreduce('sum', dir%chi_tensor_zero)
        call mw%allreduce('sum', dir%chi_tensor_one)
        call mw%allreduce('sum', dir%chi_tensor_one_sph)

        call mw%allreduce('sum', dir%I_tensor_zero)
        call mw%allreduce('sum', dir%I_tensor_one)

        if (verbosity .gt. 0) then
            t1 = walltime()
            write (lo_iou, *) '... got lowest order chi (', tochar(t1 - t0), 's)'
            t0 = t1
        end if
    end block term1

    call tmr%tock('lowest order response')

    ! Now for the terms that require actual integrals.
    itt: block
        type(lo_convolution_helper) :: ch
        integer, parameter :: nsigma = 8
        complex(r8), dimension(:), allocatable :: z0
        complex(r8) :: eta
        real(r8), dimension(:, :, :), allocatable :: buf_part_IR_III
        real(r8), dimension(:, :), allocatable :: buf_raman_II
        real(r8), dimension(:, :), allocatable :: buf_ir_IIim, buf_ir_IIre
        real(r8), dimension(:, :), allocatable :: buf_ir_IIIim, buf_ir_IIIre
        real(r8), dimension(:, :), allocatable :: grfun
        real(r8), dimension(:), allocatable :: sabfun, x, xs, y0, y1
        !real(r8), dimension(:,:), allocatable :: buf_spectral2
        !real(r8), dimension(:), allocatable :: buf_fpre,buf_om,buf_conv
        !real(r8), dimension(:), allocatable :: buf_odd_sf2a,buf_odd_sf2b
        !real(r8), dimension(:), allocatable :: buf_odd_sf3a,buf_odd_sf3b

        real(r8), dimension(81) :: buf81a
        real(r8), dimension(9) :: buf9a
        real(r8), dimension(3, 3) :: m3a, m3b, m3c
        real(r8) :: om1, om2, n1, n2, s1, s2, pref, val
        real(r8) :: qp1_radius, invf, xp, Z, dl, gm, sigma, sigma2, sigma3
        integer :: iq, ctr, b1, b2, b3
        integer :: ie, ilo, ihi, ia
        integer :: i, j

        ! radius of first q-point, for determining smearing
        qp1_radius = (3.0_r8/p%volume/real(qp%n_full_point, r8)/4.0_r8/lo_pi)**(1.0_r8/3.0_r8)
        ! buffer for the delta-function thingy
        call mem%allocate(sabfun, se%n_energy, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(grfun, [se%n_energy, dr%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        sabfun = 0.0_r8
        grfun = 0.0_r8
        ! to convert to numbers
        invf = se%n_energy/maxval(se%energy_axis)

        ! Buffers to integrate to, will be KK-transformed later
        allocate (buf_Raman_II(81, se%n_energy))
        allocate (buf_IR_IIim(se%n_energy, 9))
        allocate (buf_IR_IIre(se%n_energy, 9))
        allocate (buf_part_IR_III(se%n_energy, 9, dr%n_mode))
        allocate (buf_IR_IIIim(se%n_energy, 9))
        allocate (buf_IR_IIIre(se%n_energy, 9))
        buf_Raman_II = 0.0_r8
        buf_IR_IIim = 0.0_r8
        buf_IR_IIre = 0.0_r8
        buf_part_IR_III = 0.0_r8
        buf_IR_IIIim = 0.0_r8
        buf_IR_IIIre = 0.0_r8

        ! pre-fetch the imaginary part of the Green's function for each mode at Gamma
        ctr = 0
        grfun = 0.0_r8
        do b1 = 1, dr%n_mode
            om1 = wp%omega(b1)
            if (om1 .lt. lo_freqtol) cycle
            do ie = 1, se%n_energy
                ctr = ctr + 1
                if (mod(ctr, mw%n) .ne. mw%r) cycle
                ! fetch probing energy and real and imaginary part of the self-energy
                Z = se%energy_axis(ie)
                dl = se%re_3ph(ie, b1) + se%re_4ph(ie, b1)
                gm = se%im_3ph(ie, b1) + se%im_iso(ie, b1)
                grfun(ie, b1) = aimag(diagonal_greensfunction(om1, Z, dl, gm))
            end do
        end do
        call mw%allreduce('sum', grfun)

        ! If doing this via convolution things we need more buffers
        select case (opts%integrationtype)
        case (4:5)
            ! Set up the helper guy for the convolutions
            call ch%generate(se%energy_axis, opts%temperature, dr%n_mode)
        end select

        call tmr%tock('prep integration')

        if (verbosity .gt. 0) call lo_progressbar_init()
        ctr = 0
        qploop: do iq = 1, qp%n_irr_point

            ! If applicable, get the (pre-smeared) spectral function
            select case (opts%integrationtype)
            case (5)
                !buf_spectral2=isf%spectralfunction(:,:,iq)
                ! Pre-fetch the spectral functions
                call ch%buffer_spectral_functions(isf%spectralfunction(:, :, iq), isf%spectralfunction(:, :, iq))
            end select

            mode1loop: do b1 = 1, dr%n_mode
            mode2loop: do b2 = b1, dr%n_mode

                ! make it parallel?
                ctr = ctr + 1
                if (mod(ctr, mw%n) .ne. mw%r) cycle

                ! reset the weight function
                sabfun = 0.0_r8

                ! Two ways of creating the phase-space function S_ab: either gaussian-type integration
                ! or convolution of spectral functions. First up is the convolution way of doing it.
                select case (opts%integrationtype)
                case (5)
                    ! Smear?
                    sigma2 = qp%adaptive_sigma(qp%ip(iq)%radius, dr%iq(iq)%vel(:, b1), dr%default_smearing(b1), se%smearing_prefactor)
                    sigma3 = qp%adaptive_sigma(qp%ip(iq)%radius, dr%iq(iq)%vel(:, b2), dr%default_smearing(b2), se%smearing_prefactor)
                    sigma = sqrt(sigma2**2 + sigma3**2)

                    call ch%convolute_to_sfun(b1, b2, sigma, sabfun)
                    ! Adjust sign,weight and prefactor
                    sabfun = sabfun*qp%ip(iq)%integration_weight*lo_pi
                case default
                    ! Most probable option is with a standard Gaussian integration scheme.

                    om1 = dr%iq(iq)%omega(b1)
                    om2 = dr%iq(iq)%omega(b2)
                    ! pick some sigma for the Gaussian integration?
                    select case (se%integrationtype)
                    case (1)
                        s1 = dr%default_smearing(b1)
                        s2 = dr%default_smearing(b2)
                    case (2)
                        s1 = qp%adaptive_sigma(qp%ip(iq)%radius, dr%iq(iq)%vel(:, b1), dr%default_smearing(b1), se%smearing_prefactor)
                        s2 = qp%adaptive_sigma(qp%ip(iq)%radius, dr%iq(iq)%vel(:, b2), dr%default_smearing(b2), se%smearing_prefactor)
                    end select
                    sigma = sqrt(s1**2 + s2**2)

                    ! Planck factors.
                    n1 = lo_planck(opts%temperature, om1)
                    n2 = lo_planck(opts%temperature, om2)

                    ! Start adding terms:
                    ! w1 + w2 - Z = 0
                    ! Z = w1 + w2
                    val = om1 + om2
                    ilo = loind(val - nsigma*sigma, invf, se%n_energy)
                    ihi = hiind(val + nsigma*sigma, invf, se%n_energy)
                    pref = (n1 + n2 + 1)*qp%ip(iq)%integration_weight*lo_pi
                    do ie = ilo, ihi
                        sabfun(ie) = sabfun(ie) + pref*lo_gauss(se%energy_axis(ie), val, sigma)
                    end do

                    ! w1 + w2 + Z = 0
                    ! Z = - w1 - w2
                    ! can not happen!
                    ! actually it can contribute when you take the tail of the
                    ! Gaussian into account. Feels a little off. As long as the optical
                    ! mode (om1+om2) is larger than 4*sigma, it can not contribute. I think
                    ! it is better to not include this term. In the limit of a really dense
                    ! grid it should disappear.
                    ! val=-om1-om2
                    ! ilo=loind(val-4*sigma,invf,se%n_energy)
                    ! ihi=loind(val+4*sigma,invf,se%n_energy)
                    ! pref=(n1+n2+1)
                    ! do ie=ilo,ihi
                    !    sabfun(ie)=sabfun(ie)+pref*lo_gauss( se%energy_axis_selfenergy(ie),val,sigma )
                    ! enddo

                    ! Most likely only one of the following two can contribute
                    ! In the limit sigma -> 0 only one can contribute.
                    ! w2 - w1 + Z = 0
                    ! Z = w1-w2
                    val = om1 - om2
                    ilo = loind(val - nsigma*sigma, invf, se%n_energy)
                    ihi = hiind(val + nsigma*sigma, invf, se%n_energy)
                    pref = -(n1 - n2)*qp%ip(iq)%integration_weight*lo_pi ! note negative sign!
                    do ie = ilo, ihi
                        sabfun(ie) = sabfun(ie) + pref*lo_gauss(se%energy_axis(ie), val, sigma)
                    end do

                    ! w2 - w1 - Z = 0
                    ! Z = w2 - w1
                    val = om2 - om1
                    ilo = loind(val - nsigma*sigma, invf, se%n_energy)
                    ihi = hiind(val + nsigma*sigma, invf, se%n_energy)
                    pref = (n1 - n2)*qp%ip(iq)%integration_weight*lo_pi
                    do ie = ilo, ihi
                        sabfun(ie) = sabfun(ie) + pref*lo_gauss(se%energy_axis(ie), val, sigma)
                    end do
                end select

                ! Make sure there is no numerical noise messing things up.
                sabfun = max(sabfun, 0.0_r8)
                ! At this point we have:
                ! sabfun       = the guy with the delta-functions for (b1,b2) = (b2,b1)
                ! grfun1       = imaginary part of Green's function for all modes at Gamma
                if (di%n_eps_pair .gt. 0) then
                    ! Fetch second order Raman matrix elements
                    if (b1 .eq. b2) then
                        buf81a = dh%cmp_rmII(:, b1, b2, iq)
                    else
                        buf81a = dh%cmp_rmII(:, b1, b2, iq) + dh%cmp_rmII(:, b2, b1, iq)
                    end if
                    ! Accumulate to the spectrum
                    !call dger(81, se%n_energy, 1.0_r8, buf81a, 1, sabfun, 1, buf_Raman_II, 81)
                    call dger(81, se%n_energy, 1.0_r8, buf81a, 1, sabfun, 1, dir%I_tensor_two(:, :, b1, b2), 81)

                    ! The weird second order + anharmonic + firstorder, accumulate into a partial thing

                    !allocate(buf_part_IR_III(se%n_energy,9,dr%n_mode))
                    !@TODO add weird guy here
                end if

                if (di%n_Z_pair .gt. 0) then
                    ! Fetch second order IR matrix elements
                    if (b1 .eq. b2) then
                        buf9a = dh%cmp_irII(:, b1, b2, iq)
                    else
                        buf9a = dh%cmp_irII(:, b1, b2, iq) + dh%cmp_irII(:, b2, b1, iq)
                    end if
                    ! Accumulate to the spectrum
                    call dger(se%n_energy, 9, 1.0_r8, sabfun, 1, buf9a, 1, buf_ir_IIim, se%n_energy)

                    ! The weird second order + anharmonic + firstorder
                    if (di%n_Z_singlet .gt. 0) then
                        do b3 = 1, dr%n_mode
                            if (wp%omega(b3) .lt. lo_freqtol) cycle

                            if (b1 .eq. b2) then
                                buf9a = dh%cmp_irIII(:, b1, b2, b3, iq)
                            else
                                buf9a = dh%cmp_irIII(:, b1, b2, b3, iq) + dh%cmp_irIII(:, b2, b1, b3, iq)
                            end if
                            call dger(se%n_energy, 9, 1.0_r8, sabfun, 1, buf9a, 1, buf_part_IR_III(:, :, b3), se%n_energy)
                        end do
                    end if

                end if

                if (verbosity .gt. 0 .and. iq .lt. qp%n_irr_point) then
                    call lo_progressbar(' ... compound integral', ctr, qp%n_irr_point*dr%n_mode**2, walltime() - t0)
                end if
            end do mode2loop
            end do mode1loop
        end do qploop

        ! Sync across ranks
        call mw%allreduce('sum', dir%I_tensor_two)
        call mw%allreduce('sum', buf_IR_IIim)
        call mw%allreduce('sum', buf_part_IR_III)

        ! Add the final touch to the III term
        buf_IR_IIIim = 0.0_r8
        do i = 1, 9
        do b3 = 1, dr%n_mode
            buf_IR_IIIim(:, i) = buf_IR_IIIim(:, i) + buf_part_IR_III(:, i, b3)*grfun(:, b3) !*dh%cmp_rmI(:,imode)
        end do
        end do

        ! Cleanup in the case of convolutions
        select case (opts%integrationtype)
        case (4:5)
            call ch%destroy()
            ! deallocate(buf_fpre)
            ! deallocate(buf_om)
            ! deallocate(buf_odd_sf2a)
            ! deallocate(buf_odd_sf2b)
            ! deallocate(buf_odd_sf3a)
            ! deallocate(buf_odd_sf3b)
            ! deallocate(buf_conv)
            ! call mem%deallocate(buf_spectral2,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        end select

        if (verbosity .gt. 0) then
            t1 = walltime()
            call lo_progressbar(' ... compound integral', qp%n_irr_point*dr%n_mode**2, qp%n_irr_point*dr%n_mode**2, t1 - t0)
            t0 = t1
        end if

        ! Now I calculated the imaginary parts, get the real parts via KK transforms
        call mem%allocate(x, se%n_energy, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(xs, se%n_energy, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(z0, se%n_energy, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(y0, se%n_energy, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(y1, se%n_energy, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)

        x = se%energy_axis                  ! copy of x-axis
        xs = x**2                           ! x^2, precalculated
        eta = lo_imag*(x(2) - x(1))*1E-10_r8  ! small imaginary thing to avoid divergence
        buf_IR_IIre = 0.0_r8
        buf_IR_IIIre = 0.0_r8

        ctr = 0
        do ia = 1, 9
            do ie = 1, se%n_energy
                ctr = ctr + 1
                if (mod(ctr, mw%n) .ne. mw%r) cycle
                xp = x(ie)
                z0 = xp**2 - xs + eta
                y0 = buf_IR_IIim(:, ia)*x
                y0 = real(y0/z0, r8)
                y1 = buf_IR_IIIim(:, ia)*x
                y1 = real(y1/z0, r8)

                buf_IR_IIre(ie, ia) = lo_trapezoid_integration(x, y0)
                buf_IR_IIIre(ie, ia) = lo_trapezoid_integration(x, y1)
            end do
        end do
        call mw%allreduce('sum', buf_IR_IIre)
        call mw%allreduce('sum', buf_IR_IIIre)
        buf_IR_IIre = buf_IR_IIre*2/lo_pi ! 2/pi is part of the KK transform
        buf_IR_IIIre = buf_IR_IIIre*2/lo_pi

        if (verbosity .gt. 0) then
            t1 = walltime()
            write (lo_iou, *) '... got real part via transform (', tochar(t1 - t0), 's)'
            t0 = t1
        end if

        ! Now I have all the different components. Just make sure to store them in the right place
        ! and add eventual prefactors.

        ! Store in the right place, first the IR
        do ie = 1, se%n_energy
            do i = 1, 3
            do j = 1, 3
                m3a(i, j) = buf_IR_IIre(ie, (i - 1)*3 + j)
                m3b(i, j) = buf_IR_IIim(ie, (i - 1)*3 + j)
            end do
            end do
            m3c = m3c + m3b
            dir%chi_tensor_two(:, :, ie) = cmplx(m3a, m3b, r8)
            do i = 1, 3
            do j = 1, 3
                m3a(i, j) = buf_IR_IIIre(ie, (i - 1)*3 + j)
                m3b(i, j) = buf_IR_IIIim(ie, (i - 1)*3 + j)
            end do
            end do
            dir%chi_tensor_three(:, :, ie) = cmplx(m3a, m3b, r8)
        end do

        ! Fill out all the modes in the Raman
        do b1 = 1, se%n_mode
        do b2 = b1 + 1, se%n_mode
            dir%I_tensor_two(:, :, b2, b1) = dir%I_tensor_two(:, :, b1, b2)
        end do
        end do

        if (verbosity .gt. 0) then
            t1 = walltime()
            write (lo_iou, *) '... rearranged tensors (', tochar(t1 - t0), 's)'
            t0 = t1
        end if

        ! Cleanup
        call mem%deallocate(sabfun, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(grfun, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(x, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(xs, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(z0, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(y0, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(y1, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    end block itt

    call tmr%tock('higher order response')

    ! Make sure the tensors are symmetric and clean and nice? Always a neat idea.
    finalize: block
        dir%chi_tensor_zero = dir%chi_tensor_zero/p%volume
        dir%chi_tensor_one = dir%chi_tensor_one/p%volume
        dir%chi_tensor_one_sph = dir%chi_tensor_one_sph/p%volume
        dir%chi_tensor_two = dir%chi_tensor_two*2.0_r8/p%volume

        dir%I_tensor_zero = dir%I_tensor_zero/lo_twopi
        dir%I_tensor_one = dir%I_tensor_one/lo_twopi
        dir%I_tensor_two = dir%I_tensor_two*2.0_r8/lo_pi
    end block finalize

    call mem%tock(__FILE__, __LINE__, mw%comm)
end subroutine

! just shorthand function for the loop indices
pure function loind(value, invf, nmax) result(i)
    real(r8), intent(in) :: value
    real(r8), intent(in) :: invf
    integer, intent(in) :: nmax
    integer :: i

    i = floor(value*invf) + 1
    i = min(i, nmax)
    i = max(i, 1)
end function
pure function hiind(value, invf, nmax) result(i)
    real(r8), intent(in) :: value
    real(r8), intent(in) :: invf
    integer, intent(in) :: nmax
    integer :: i

    i = ceiling(value*invf) + 1
    i = min(i, nmax)
    i = max(i, 1)
end function

!> dump a message regarding memory usage
subroutine memdump(obj, bytes, mw, iou)
    !> what are we measuring the size of
    character(len=*), intent(in) :: obj
    !> size of object, on this rank, in bytes
    integer, intent(in) :: bytes
    !> MPI helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> unit to write to
    integer, intent(in) :: iou

    real(r8), parameter :: toMiB = 1.0_r8/1024_r8**2
    real(r8), dimension(mw%n) :: mb

    mb = 0.0_r8
    mb(mw%r + 1) = bytes*toMiB
    call mw%allreduce('sum', mb)
    if (mw%talk) then
        write (iou, "(1X,A30,':',3(1X,F20.6))") trim(adjustl(obj)), sum(mb)/mw%n, maxval(mb), minval(mb)
    end if
end subroutine memdump

!> size in memory, in bytes
function dir_size_in_mem(dir) result(mem)
    !> atom
    class(lo_dielectric_response), intent(in) :: dir
    !> size in memory, bytes
    integer(i8) :: mem

    mem = 0
    mem = mem + storage_size(dir)
    if (allocated(dir%chi_tensor_zero)) mem = mem + storage_size(dir%chi_tensor_zero)*size(dir%chi_tensor_zero)
    if (allocated(dir%chi_tensor_one)) mem = mem + storage_size(dir%chi_tensor_one)*size(dir%chi_tensor_one)
    if (allocated(dir%chi_tensor_one_sph)) mem = mem + storage_size(dir%chi_tensor_one_sph)*size(dir%chi_tensor_one_sph)
    if (allocated(dir%chi_tensor_two)) mem = mem + storage_size(dir%chi_tensor_two)*size(dir%chi_tensor_two)
    if (allocated(dir%chi_tensor_three)) mem = mem + storage_size(dir%chi_tensor_three)*size(dir%chi_tensor_three)
    if (allocated(dir%I_tensor_zero)) mem = mem + storage_size(dir%I_tensor_zero)*size(dir%I_tensor_zero)
    if (allocated(dir%I_tensor_one)) mem = mem + storage_size(dir%I_tensor_one)*size(dir%I_tensor_one)
    if (allocated(dir%I_tensor_two)) mem = mem + storage_size(dir%I_tensor_two)*size(dir%I_tensor_two)
    mem = mem/8
end function

!> size in memory, in bytes
function dh_size_in_mem(dh) result(mem)
    !> atom
    class(diel_helper), intent(in) :: dh
    !> size in memory, bytes
    integer(i8) :: mem

    mem = 0
    mem = mem + storage_size(dh)
    if (allocated(dh%cmp_rmI)) mem = mem + storage_size(dh%cmp_rmI)*size(dh%cmp_rmI)
    if (allocated(dh%cmp_rmII)) mem = mem + storage_size(dh%cmp_rmII)*size(dh%cmp_rmII)
    if (allocated(dh%cmp_irI)) mem = mem + storage_size(dh%cmp_irI)*size(dh%cmp_irI)
    if (allocated(dh%cmp_irI_sph)) mem = mem + storage_size(dh%cmp_irI_sph)*size(dh%cmp_irI_sph)
    if (allocated(dh%cmp_irII)) mem = mem + storage_size(dh%cmp_irII)*size(dh%cmp_irII)
    if (allocated(dh%cmp_irIII)) mem = mem + storage_size(dh%cmp_irIII)*size(dh%cmp_irIII)
    if (allocated(dh%ugv_gamma)) mem = mem + storage_size(dh%ugv_gamma)*size(dh%ugv_gamma)
    if (allocated(dh%ugv)) mem = mem + storage_size(dh%ugv)*size(dh%ugv)
    if (allocated(dh%shift)) mem = mem + storage_size(dh%shift)*size(dh%shift)
    if (allocated(dh%linewidth)) mem = mem + storage_size(dh%linewidth)*size(dh%linewidth)
    if (allocated(dh%spectral_function)) mem = mem + storage_size(dh%spectral_function)*size(dh%spectral_function)
    if (allocated(dh%omega_axis)) mem = mem + storage_size(dh%omega_axis)*size(dh%omega_axis)
    mem = mem/8
end function

end module
