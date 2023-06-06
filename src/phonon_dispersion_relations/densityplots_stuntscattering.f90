submodule(densityplots) densityplots_stuntscattering
use lo_randomnumbers, only: lo_mersennetwister
implicit none
contains

!> Inelastic scattering spectra, at the harmonic level
module subroutine inelastic_spectra(bs, uc, fc, qp, dr, mw, mem, verbosity)
    !> bandstructur
    type(lo_phonon_bandstructure), intent(inout) :: bs
    !> unitcell
    type(lo_crystalstructure), intent(in) :: uc
    !> forceconstant
    type(lo_forceconstant_secondorder), intent(inout) :: fc
    !> q-mesh
    class(lo_qpoint_mesh), intent(in) :: qp
    !> phonon dispersions on a grid
    type(lo_phonon_dispersions), intent(in) :: dr
    !> MPI helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem
    !> Talk a lot?
    integer, intent(in) :: verbosity

    ! Some parameters I should make input eventually
    real(r8), parameter :: temperature = 300.0_r8       !< Temperature
    real(r8), parameter :: maxomegafactor = 1.3_r8      !< Where to cut the spectrum in height
    real(r8), parameter :: qsigma = 0.03_r8             !< Size of q-resolution
    real(r8), parameter :: esigma = 0.25_r8/lo_frequency_Hartree_to_meV
    integer, parameter :: ne = 2000                     !< Number of points on the energy axis
    integer, parameter :: nqsph = 150                    !< How many additional q-points per point on the path

    real(r8), dimension(:, :, :), allocatable :: thermal_sigma
    real(r8), dimension(:, :), allocatable :: qvcloud
    real(r8) :: timer

    init: block
        type(lo_mersennetwister) :: tw
        real(r8) :: sigma
        integer :: i, j, k

        if (verbosity .gt. 0) then
            write (*, *) ''
            write (*, *) 'Calculating fake scattering spectrum'
        end if
        timer = walltime()

        call tw%init(iseed=mw%r, rseed=walltime())

        ! Generate a Gaussian cloud of q-points, centered around zero. What
        ! we will use for the q-resolution.
        sigma = uc%bz%rmin*qsigma
        allocate (qvcloud(3, nqsph))
        qvcloud = 0.0_r8
        i = 0
        do j = 1, nqsph
        do k = 1, 3
            i = i + 1
            if (mod(i, mw%n) .ne. mw%r) cycle
            qvcloud(k, j) = tw%rnd_gaussian_real(0.0_r8, sigma)
        end do
        end do
        call mw%allreduce('sum', qvcloud)

        if (verbosity .gt. 0) then
            write (*, *) '... got ', tochar(nqsph*bs%n_point), ' q-points (', tochar(walltime() - timer + 1), 's)'
        end if

        ! Then prepare some space for the intensity
        allocate (bs%energy_axis(ne))
        allocate (bs%spectral_function(bs%n_point, ne))
        allocate (bs%spectral_function_with_prefactor(bs%n_point, ne))
        bs%energy_axis = 0.0_r8
        bs%spectral_function = 0.0_r8
        bs%spectral_function_with_prefactor = 0.0_r8
        call lo_linspace(0.0_r8, dr%omega_max*maxomegafactor, bs%energy_axis)

        ! Also get the thermal displacement factor
        allocate (thermal_sigma(3, 3, uc%na))
        thermal_sigma = 0.0_r8
        call dr%thermal_displacement_matrix(qp, uc, temperature, thermal_sigma, mw, mem)

        if (verbosity .gt. 0) then
            write (*, *) '... set up structurefactor calculation'
        end if
    end block init

    ! Now calculate the actual intensity
    intens: block
        type(lo_phonon_dispersions_qpoint) :: ompoint
        real(r8), dimension(:), allocatable :: dy0, dy1, dy2, dz0
        real(r8), dimension(3) :: v0
        real(r8) :: f0, t0
        integer :: iq, jq, imode

        allocate (dy0(ne))
        allocate (dy1(ne))
        allocate (dy2(ne))
        allocate (dz0(ne))

        t0 = walltime()
        if (verbosity .gt. 0) call lo_progressbar_init()

        do iq = 1, bs%n_point
            if (mod(iq, mw%n) .ne. mw%r) cycle

            ! First thing, get the intensities of all the surrounding q-points of this point
            dy1 = 0.0_r8
            dy2 = 0.0_r8
            do jq = 1, nqsph
                ! Get the absolute q-vector
                v0 = bs%q(iq)%r + qvcloud(:, jq)
                call ompoint%generate(fc, uc, qvec=v0, mem=mem)

                ! Accumulate
                do imode = 1, dr%n_mode
                    f0 = thermal_pref(v0, ompoint%egv(:, imode), uc, temperature, ompoint%omega(imode), thermal_sigma)
                    dy0 = lo_gauss(bs%energy_axis, ompoint%omega(imode), esigma)
                    dy1 = dy1 + dy0
                    dy2 = dy2 + dy0*f0
                end do
            end do

            if (verbosity .gt. 0 .and. iq .lt. bs%n_point) then
                call lo_progressbar(' ... inelastic intensity', iq, bs%n_point, walltime() - t0)
            end if

            ! Store
            bs%spectral_function(iq, :) = dy2
            bs%spectral_function_with_prefactor(iq, :) = dy1
        end do

        call mw%allreduce('sum', bs%spectral_function)
        call mw%allreduce('sum', bs%spectral_function_with_prefactor)

        if (verbosity .gt. 0) then
            call lo_progressbar(' ... inelastic intensity', bs%n_point, bs%n_point, walltime() - t0)
        end if

    end block intens
end subroutine

! The magical thermal prefactor
function thermal_pref(bigQ, eigenvector, uc, temperature, omega, sigma) result(pref)
    real(r8), dimension(3), intent(in) :: bigQ
    complex(r8), dimension(:), intent(in) :: eigenvector
    type(lo_crystalstructure), intent(in) :: uc
    real(r8), intent(in) :: temperature
    real(r8), intent(in) :: omega
    real(r8), dimension(:, :, :), intent(in) :: sigma
    real(r8) :: pref

    integer :: i
    complex(r8) :: c0, c1, c2
    real(r8) :: f0, f1, xs, dw

    pref = 1.0_r8
    ! Add the thermal factor
    pref = pref*(2*lo_planck(temperature, omega) + 1.0_r8)

    ! Cross-product guy
    c2 = 0.0_r8
    do i = 1, uc%na
        ! Grab the cross-section?
        xs = uc%inelastic_neutron_cross_section(i)*uc%invsqrtmass(i)
        ! Debye-Waller factor?
        dw = exp(-0.5_r8*dot_product(bigQ, matmul(sigma(:, :, i), bigQ))*lo_twopi**2)

        f1 = dot_product(lo_twopi*bigQ, uc%rcart(:, i))
        c0 = cmplx(cos(f1), sin(f1), r8)
        c1 = dot_product(lo_twopi*bigQ, eigenvector((i - 1)*3 + 1:i*3))
        c2 = c2 + c0*c1*xs
    end do
    pref = pref*abs(conjg(c2)*c2)
    ! And the 1/omega factor I do not like
    f0 = 0.01_r8*lo_frequency_THz_to_Hartree
    if (omega .gt. f0) then
        pref = pref/omega
    else
        pref = pref/f0
    end if

!    ! Calculate the prefactors
!    do q1=1,bs%nptot
!    do i=1,bs%nb
!        f0=0.0_r8
!        do j=1,uc%na
!            f1=dot_product(lo_twopi*bs%q(q1)%v,uc%rcart(:,j))
!            c0=dcmplx( cos(f1), sin(f1) )
!            ii=(j-1)*3+1
!            jj=j*3
!            c1=dot_product(lo_twopi*bs%q(q1)%v,bs%p(q1)%egv(ii:jj,i))
!            c2=c0*c1*uc%inelastic_neutron_cross_section(j)
!            f0=f0+abs(conjg(c2)*c2)
!        enddo
!        bs%p(q1)%thermal_prefactor(i)=f0
!    enddo
!    enddo
end function

end submodule
