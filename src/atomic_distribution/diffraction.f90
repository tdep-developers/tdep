!> Get the symmetry operations for a structure in a neat way.
module diffraction
use konstanter
!use geometryfunctions
use mpi_wrappers, only: lo_mpi_helper
use gottochblandat !, only : lo_identitymatrix
use type_crystalstructure, only: lo_crystalstructure
!use type_blas_lapack_wrappers, only: lo_dgesvd,lo_dgels,lo_dgglse,lo_dgelss,lo_dsyev
use type_mdsim, only: lo_mdsim
use type_distancetable
use hdf5_wrappers, only: lo_hdf5_helper
!use dump_data, only: lo_dump_gnuplot_2d_real
implicit none

private
public :: lo_powderdiffraction

! Standard Cu k-alpha wavelength, in A
real(r8), parameter :: lambda = 1.5405981_r8*lo_A_to_Bohr
!> Angles I care about
real(r8), parameter :: min2theta = 5.0_r8, max2theta = 60.0_r8
!> How many point on my diffraction axis?
integer, parameter :: npts = 16000

type lo_powderdiffraction
    !> number of points on the plot
    integer :: nx
    !> norm of G-vectors
    real(r8), dimension(:), allocatable :: norm_G_axis
    !> twotheta axis
    real(r8), dimension(:), allocatable :: twotheta_axis
    !> ideal intensity
    real(r8), dimension(:), allocatable :: intensity_ideal
    !> thermal intensity
    real(r8), dimension(:), allocatable :: intensity
    !> Site-projected intensity?
contains
    procedure :: generate
    procedure :: write_to_hdf5
end type

contains

subroutine generate(df, uc, ss, sim, mw, verbosity)
    !> powder diffraction pattern
    class(lo_powderdiffraction), intent(out) :: df
    !> unitcell
    type(lo_crystalstructure), intent(in) :: uc
    !> supercell
    type(lo_crystalstructure), intent(in) :: ss
    !> MD data
    type(lo_mdsim), intent(in) :: sim
    !> MPI helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> talk?
    integer, intent(in) :: verbosity

    real(r8), dimension(:, :), allocatable :: Gvecs
    real(r8), dimension(:), allocatable :: Gnorm
    real(r8) :: Gmin, Gmax, sigma
    real(r8) :: timer, t0, t1
    integer :: nG

    timer = walltime()
    t0 = timer
    t1 = timer

    init: block
        type(lo_distancetable) :: dt
        real(r8), dimension(3, 1) :: dumr
        integer :: i, l
        ! Some space for axes
        df%nx = npts
        allocate (df%norm_G_axis(df%nx))
        allocate (df%twotheta_axis(df%nx))
        allocate (df%intensity_ideal(df%nx))
        allocate (df%intensity(df%nx))
        df%norm_G_axis = 0.0_r8
        df%twotheta_axis = 0.0_r8
        df%intensity = 0.0_r8
        df%intensity_ideal = 0.0_r8

        ! Figure out max and min G-vectors
        Gmin = sin(min2theta*lo_pi/360.0_r8)*2/lambda
        Gmax = sin(max2theta*lo_pi/360.0_r8)*2/lambda
        call lo_linspace(Gmin, Gmax, df%norm_G_axis)
        do i = 1, df%nx
            df%twotheta_axis(i) = asin(df%norm_G_axis(i)*lambda*0.5_r8)*360.0_r8/lo_pi
        end do

        ! Get the list of G-vectors?
        dumr = 0.0_r8
        call dt%generate(dumr, ss%reciprocal_latticevectors, Gmax + lo_tol, verbosity, mw)
        nG = 0
        do i = 1, dt%particle(1)%n
            if (dt%particle(1)%d(i) .lt. Gmin) cycle
            if (dt%particle(1)%d(i) .gt. Gmax) cycle
            nG = nG + 1
        end do

        allocate (Gvecs(3, nG))
        allocate (Gnorm(ng))
        Gvecs = 0.0_r8
        Gnorm = 0.0_r8
        l = 0
        do i = 1, dt%particle(1)%n
            if (dt%particle(1)%d(i) .lt. Gmin) cycle
            if (dt%particle(1)%d(i) .gt. Gmax) cycle
            l = l + 1
            Gvecs(:, l) = dt%particle(1)%v(:, i)
            Gvecs(:, l) = matmul(ss%inv_reciprocal_latticevectors, Gvecs(:, l))*lo_twopi
            Gnorm(l) = dt%particle(1)%d(i)
        end do

        ! Decide on some tiny smearing?
        sigma = (df%norm_G_axis(2) - df%norm_G_axis(1))*5

        if (verbosity .gt. 0) then
            write (*, *) ''
            write (*, *) 'POWDER DIFFRACTION PATTERN'
            write (*, *) '      lambda:', lambda*lo_bohr_to_A
            write (*, *) '        Gmin:', Gmin
            write (*, *) '        Gmax:', Gmax
            write (*, *) '    N G-vecs:', nG
        end if
    end block init

    avgpattern: block
        complex(r8), dimension(:), allocatable :: dc0
        complex(r8) :: c0
        real(r8), dimension(:), allocatable :: sfactor0, sfactor1
        real(r8), dimension(:), allocatable :: dr0
        real(r8) :: f0, f1, invf
        integer :: iG, iatm, i, t, ilo, ihi

        t0 = walltime()

        invf = df%nx/(Gmax - Gmin)
        allocate (dr0(ss%na))
        allocate (dc0(ss%na))
        allocate (sfactor0(nG))
        allocate (sfactor1(nG))
        sfactor0 = 0.0_r8
        sfactor1 = 0.0_r8

        if (verbosity .gt. 0) call lo_progressbar_init()

        ! First do the idealized structure factors
        do iG = 1, nG
            if (mod(iG, mw%n) .ne. mw%r) cycle
            ! get structure factor?
            dr0 = matmul(Gvecs(:, ig), ss%r)
            dc0 = cmplx(cos(dr0), sin(dr0), r8)
            c0 = sum(dc0)
            f0 = abs(c0*conjg(c0)) ! intensity to add
            sfactor0(iG) = sfactor0(iG) + f0
        end do
        call mw%allreduce('sum', sfactor0)

        ! Then the averaged structure factors
        do t = 1, sim%nt
            do iG = 1, nG
                if (mod(iG, mw%n) .ne. mw%r) cycle
                ! get structure factor?
                dr0 = matmul(Gvecs(:, ig), sim%r(:, :, t))
                dc0 = cmplx(cos(dr0), sin(dr0), r8)
                c0 = sum(dc0)
                f0 = abs(c0*conjg(c0)) ! intensity to add
                sfactor1(iG) = sfactor1(iG) + f0
            end do
            if (verbosity .gt. 0 .and. t .lt. sim%nt) then
                call lo_progressbar(' ... diffraction pattern', t, sim%nt, walltime() - t0)
            end if
        end do
        call mw%allreduce('sum', sfactor1)

        if (verbosity .gt. 0) then
            t1 = walltime()
            call lo_progressbar(' ... diffraction pattern', sim%nt, sim%nt, t1 - t0)
            t0 = t1
        end if

        ! And add to histogram
        df%intensity = 0.0_r8
        df%intensity_ideal = 0.0_r8
        do ig = 1, nG
            if (mod(iG, mw%n) .ne. mw%r) cycle
            f0 = (Gnorm(iG) - Gmin - 4*sigma)*invf
            f1 = (Gnorm(iG) - Gmin + 4*sigma)*invf
            ilo = floor(f0 + 1)
            ihi = ceiling(f1) + 1
            ilo = max(ilo, 1)
            ihi = min(ihi, df%nx)
            do i = ilo, ihi
                df%intensity_ideal(i) = df%intensity_ideal(i) + sfactor0(ig)*lo_gauss(df%norm_G_axis(i), Gnorm(iG), sigma)
                df%intensity(i) = df%intensity(i) + sfactor1(ig)*lo_gauss(df%norm_G_axis(i), Gnorm(iG), sigma)
            end do
        end do
        call mw%allreduce('sum', df%intensity_ideal)
        call mw%allreduce('sum', df%intensity)

        df%intensity = df%intensity/lo_trapezoid_integration(df%norm_G_axis, df%intensity)
        df%intensity_ideal = df%intensity_ideal/lo_trapezoid_integration(df%norm_G_axis, df%intensity_ideal)

        if (verbosity .gt. 0) then
            t1 = walltime()
            write (*, *) '... binned pattern (', tochar(t1 - t0), 's)'
            t0 = t1
        end if
    end block avgpattern
end subroutine

!> Dump everything to file
subroutine write_to_hdf5(df, filename)
    !> pair distribution function
    class(lo_powderdiffraction), intent(in) :: df
    !> filename
    character(len=*), intent(in) :: filename

    type(lo_hdf5_helper) :: h5

    call h5%init(__FILE__, __LINE__)
    call h5%open_file('write', trim(filename))

    call h5%store_data(df%twotheta_axis, h5%file_id, 'twotheta')
    call h5%store_data(df%intensity, h5%file_id, 'intensity')
    call h5%store_data(df%intensity_ideal, h5%file_id, 'intensity_ideal')

    call h5%close_file()
    call h5%destroy(__FILE__, __LINE__)
end subroutine

end module
