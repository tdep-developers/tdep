module pair_distribution
use konstanter, only: r8, lo_huge, lo_hugeint, lo_pi, lo_sqtol, lo_tol, lo_status, lo_bohr_to_A, lo_time_au_to_fs
use gottochblandat, only: walltime, lo_progressbar_init, lo_progressbar, lo_trapezoid_integration, lo_gauss, tochar, lo_linspace, &
                          lo_return_unique
use type_crystalstructure, only: lo_crystalstructure
use type_mdsim, only: lo_mdsim
use pairmapping, only: lo_pairmapping
use lo_memtracker, only: lo_mem_helper
use mpi_wrappers, only: lo_mpi_helper, lo_stop_gracefully
use hdf5_wrappers, only: lo_hdf5_helper

implicit none
private
public :: lo_pair_distribution

!> pair distribution per shell
type lo_pair_distribution_shell
    !> smallest distance
    real(r8) :: rmin = -lo_huge
    !> largest distance
    real(r8) :: rmax = -lo_huge
    !> shell-resolved histogram
    real(r8), dimension(:), allocatable :: x
    real(r8), dimension(:), allocatable :: y
end type

!> symmetry decomposed pair distribution function
type lo_pair_distribution
    !> number of bins in the histogram
    integer :: nbin = -lo_hugeint

    !> shell-resolved pdf
    type(lo_pair_distribution_shell), dimension(:), allocatable :: sh
    !> smearing parameter for binning
    real(r8) :: sigma = -lo_huge
contains
    !> Bin the timesteps
    procedure :: bin
    !> Write it to file
    procedure :: write_to_hdf5
end type

contains

!> Bin the pair distribution function
subroutine bin(pdf, uc, ss, pm, sim, nbin, nodiffusion, mw, mem, verbosity)
    !> the histogram stuff
    class(lo_pair_distribution), intent(out) :: pdf
    !> unitcell
    type(lo_crystalstructure), intent(in) :: uc
    !> supercell
    type(lo_crystalstructure), intent(in) :: ss
    !> pair mapping
    type(lo_pairmapping), intent(in) :: pm
    !> the simulation data
    type(lo_mdsim), intent(in) :: sim
    !> how many bins in the histogram
    integer, intent(in) :: nbin
    !> always assume no diffusion?
    logical, intent(in) :: nodiffusion
    !> mpi helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem
    !> talk a lot?
    integer, intent(in) :: verbosity

    real(r8), parameter :: buf = 0.1_r8 ! buffer the histograms a little bit
    real(r8) :: timer, t0, t1

    ! start timers
    timer = walltime()
    t0 = timer
    t1 = timer

    ! Set some basics
    init: block
        integer :: ish

        if (verbosity .gt. 0) then
            write (*, *) ''
            write (*, *) 'Calculating pair distribution function'
        end if

        ! How many bins?
        pdf%nbin = nbin

        ! Space for shells
        allocate (pdf%sh(pm%n_shell))
        do ish = 1, pm%n_shell
            pdf%sh(ish)%rmin = lo_huge
            pdf%sh(ish)%rmax = -lo_huge
            allocate (pdf%sh(ish)%x(pdf%nbin))
            allocate (pdf%sh(ish)%y(pdf%nbin))
            pdf%sh(ish)%x = 0.0_r8
            pdf%sh(ish)%y = 0.0_r8
        end do
    end block init

    ! Do a first sweep
    firstpass: block
        real(r8), dimension(:, :), allocatable :: dr0
        real(r8) :: f0, f1
        integer :: ish, ctr, ipair, t, ctrtot

        ctrtot = 0
        do ish = 1, pm%n_shell
            ctrtot = ctrtot + pm%sh(ish)%n_ss_pair
        end do

        call mem%allocate(dr0, [3, sim%nt], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        dr0 = 0.0_r8

        if (verbosity .gt. 0) call lo_progressbar_init()
        ctr = 0
        do ish = 1, pm%n_shell
        do ipair = 1, pm%sh(ish)%n_ss_pair
            ctr = ctr + 1
            if (mod(ctr, mw%n) .ne. mw%r) cycle
            ! fetch the trajectory?
            call extract_trajectory(pm, ish, ipair, sim, dr0)
            do t = 1, sim%nt
                f0 = norm2(dr0(:, t))
                pdf%sh(ish)%rmin = min(pdf%sh(ish)%rmin, f0)
                pdf%sh(ish)%rmax = max(pdf%sh(ish)%rmax, f0)
            end do

            if (verbosity .gt. 0 .and. ctr .lt. ctrtot) then
                call lo_progressbar(' ... measuring span', ctr, ctrtot, walltime() - t0)
            end if
        end do
        end do

        ! sync and build axis
        f0 = 0.0_r8
        ctr = 0
        do ish = 1, pm%n_shell
            call mw%allreduce('min', pdf%sh(ish)%rmin)
            call mw%allreduce('max', pdf%sh(ish)%rmax)
            if (norm2(pm%sh(ish)%r) .gt. lo_tol) then
                ctr = ctr + 1
                f1 = (pdf%sh(ish)%rmax - pdf%sh(ish)%rmin)/real(pdf%nbin, r8)
                f0 = max(f0, f1)
            end if
        end do
        pdf%sigma = f0*2

        do ish = 1, pm%n_shell
            pdf%sh(ish)%rmin = pdf%sh(ish)%rmin - 7*pdf%sigma
            pdf%sh(ish)%rmax = pdf%sh(ish)%rmax + 7*pdf%sigma
            call lo_linspace(pdf%sh(ish)%rmin, pdf%sh(ish)%rmax, pdf%sh(ish)%x)
        end do

        if (verbosity .gt. 0) then
            t1 = walltime()
            call lo_progressbar(' ... measuring span', ctrtot, ctrtot, t1 - t0)
            t0 = t1
        end if

        call mem%deallocate(dr0, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    end block firstpass

    ! Do a second sweep
    secondpass: block
        real(r8), dimension(:, :), allocatable :: dr0
        real(r8) :: f0, f1, f2, f3, invf
        integer :: ish, ctr, ipair, t, ctrtot
        integer :: i, j, k, l, ilo, ihi

        ctrtot = 0
        do ish = 1, pm%n_shell
            ctrtot = ctrtot + pm%sh(ish)%n_ss_pair
        end do

        call mem%allocate(dr0, [3, sim%nt], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        dr0 = 0.0_r8

        if (verbosity .gt. 0) call lo_progressbar_init()
        ctr = 0
        do ish = 1, pm%n_shell
            invf = pdf%nbin/(pdf%sh(ish)%rmax - pdf%sh(ish)%rmin)
            do ipair = 1, pm%sh(ish)%n_ss_pair
                ctr = ctr + 1
                if (mod(ctr, mw%n) .ne. mw%r) cycle
                ! fetch the trajectory?
                call extract_trajectory(pm, ish, ipair, sim, dr0)
                do t = 1, sim%nt
                    f0 = norm2(dr0(:, t))  ! pair distance
                    f1 = (f0 - pdf%sh(ish)%rmin - 4*pdf%sigma)*invf
                    f2 = (f0 - pdf%sh(ish)%rmin + 4*pdf%sigma)*invf
                    ilo = floor(f1) + 1
                    ihi = ceiling(f2) + 1
                    do i = ilo, ihi
                        f3 = lo_gauss(pdf%sh(ish)%x(i), f0, pdf%sigma)
                        pdf%sh(ish)%y(i) = pdf%sh(ish)%y(i) + f3
                    end do
                end do

                if (verbosity .gt. 0 .and. ctr .lt. ctrtot - 1) then
                    call lo_progressbar(' ... binning', ctr, ctrtot, walltime() - t0)
                end if
            end do
        end do

        ! sync and build axis
        do ish = 1, pm%n_shell
            call mw%allreduce('sum', pdf%sh(ish)%y)
        end do

        if (verbosity .gt. 0) then
            t1 = walltime()
            call lo_progressbar(' ... binning', ctrtot, ctrtot, t1 - t0)
            t0 = t1
        end if

        call mem%deallocate(dr0, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    end block secondpass

    ! Normalize so that it integrates to the correct number of atoms in spherical coordinates
    finalize: block

        real(r8) :: f0
        integer :: ish, t

        do ish = 1, pm%n_shell
            f0 = lo_trapezoid_integration(pdf%sh(ish)%x, pdf%sh(ish)%y*pdf%sh(ish)%x*pdf%sh(ish)%x)*4.0_r8*lo_pi
            pdf%sh(ish)%y = pdf%sh(ish)%y*pm%sh(ish)%n_ss_pair/f0
        end do
    end block finalize

end subroutine

!> Dump everything to file
subroutine write_to_hdf5(pdf, pm, uc, sim, filename)
    !> pair distribution function
    class(lo_pair_distribution), intent(in) :: pdf
    !> pair mapping
    type(lo_pairmapping), intent(in) :: pm
    !> unitcell
    type(lo_crystalstructure), intent(in) :: uc
    !> simulation
    type(lo_mdsim), intent(in) :: sim
    !> filename
    character(len=*), intent(in) :: filename

    type(lo_hdf5_helper) :: h5
    real(r8), dimension(:), allocatable :: dr0,x,y0,y1,y2
    real(r8) :: xlo,xhi,sigma,f0
    integer :: ish, i, nshell, ii, ctr, j
    character(len=2000) :: str

    call h5%init(__FILE__, __LINE__)
    call h5%open_file('write', trim(filename))

    ! Count number of relevant shells?
    nshell = 0
    do i = 1, pm%n_shell
        if (norm2(pm%sh(i)%r) .gt. lo_tol) then
            nshell = nshell + 1
        end if
    end do

    ! Store number of shells?
    if ( nshell .lt. 30 ) then
        call h5%store_attribute(nshell, h5%file_id, 'n_shell')

        allocate (dr0(pdf%nbin))
        dr0 = 0.0_r8

        ! Store the symmetry unique shells

        nshell = 0
        do ish = 1, pm%n_shell
            if (norm2(pm%sh(ish)%r) .lt. lo_tol) cycle
            nshell = nshell + 1
            ! store things per shell
            call h5%open_group('write', 'unique_shell_'//tochar(nshell))

            ! store the prototype vector
            call h5%store_data(pm%sh(ish)%r*lo_bohr_to_A, h5%group_id, 'pair_vector', enhet='A')

            ! store histogram
            call h5%store_data(pdf%sh(ish)%x*lo_bohr_to_A, h5%group_id, 'r', enhet='A')
            call h5%store_data(pdf%sh(ish)%y/lo_bohr_to_A, h5%group_id, 'pdf', enhet='1/A')

            ! create ideal peak at the average position
            dr0 = 0.0_r8
            do i = 1, pdf%nbin
                dr0(i) = lo_gauss(pdf%sh(ish)%x(i), norm2(pm%sh(ish)%r), pdf%sigma)
            end do
            dr0 = dr0/lo_bohr_to_A
            call h5%store_data(dr0, h5%group_id, 'pdf_ideal', enhet='1/A')

            call h5%close_group()
        end do

        deallocate (dr0)
    endif

    ! Store as species shells as well.
    allocate(x(pdf%nbin*10))
    allocate(y0(pdf%nbin*10))
    allocate(y1(pdf%nbin*10))
    allocate(y2(pdf%nbin*10))

    str=trim(pm%species_pair_label(1))
    do i=2,pm%n_species_pair
        str=trim(str)//"|"//trim(pm%species_pair_label(i))
    enddo
    call h5%store_attribute(trim(str),h5%file_id,'pair_shell_label')

    do i=1,pm%n_species_pair
        call h5%open_group('write', 'shell_'//trim(pm%species_pair_label(i)))

        ! Count min and max x?
        xhi=-lo_huge
        xlo=lo_huge
        do ish=1,pm%n_shell
            if ( pm%sh(ish)%index_species_pair .ne. i ) cycle
            if ( norm2(pm%sh(ish)%r) .lt. lo_tol ) cycle
            xhi=max(xhi,maxval(pdf%sh(ish)%x))
            xlo=min(xlo,minval(pdf%sh(ish)%x))
        enddo
        call lo_linspace(xlo,xhi,x)
        sigma=x(3)-x(1)

        y0=0.0_r8
        y1=0.0_r8

        ! Accumulate histograrms
        ctr=0
        j=0
        do ish=1,pm%n_shell
            if ( pm%sh(ish)%index_species_pair .ne. i ) cycle
            if ( norm2(pm%sh(ish)%r) .lt. lo_tol ) cycle
            xhi=max(xhi,maxval(pdf%sh(ish)%x))
            xlo=min(xlo,minval(pdf%sh(ish)%x))
            ctr=ctr + pm%sh(ish)%n_ss_pair
            ! Accumulate pdf
            call put_shape_on_new_mesh(pdf%sh(ish)%x,pdf%sh(ish)%y,x,y2)
            y1=y1+y2
            ! Keep track of the ideal peak locations
            j=j+1
            y0(j)=norm2(pm%sh(ish)%r)
        enddo

        call h5%store_data(x*lo_bohr_to_A, h5%group_id, 'r', enhet='A')
        call h5%store_data(y1/lo_bohr_to_A/4/lo_pi, h5%group_id, 'pdf', enhet='1/A')
        call h5%store_data(y0(1:j)*lo_bohr_to_A, h5%group_id, 'ideal_peak_locations', enhet='A')

        call h5%close_group()
    enddo




    call h5%close_file()
    call h5%destroy(__FILE__, __LINE__)
end subroutine

!> get the trajectory for one pair
subroutine extract_trajectory(pm, ish, ipair, sim, trajectory)
    !> pair mapping
    type(lo_pairmapping), intent(in) :: pm
    !> coordination shell
    integer, intent(in) :: ish
    !> supercell pair
    integer, intent(in) :: ipair
    !> simulation data
    type(lo_mdsim), intent(in) :: sim
    !> trajectory
    real(r8), dimension(:, :), intent(out) :: trajectory

    real(r8), dimension(3) :: v0, v1
    integer :: a1, a2, t

    trajectory = 0.0_r8
    a1 = pm%sh(ish)%pairind(1, ipair)
    a2 = pm%sh(ish)%pairind(2, ipair)
    v0 = pm%sh(ish)%pairvec(:, ipair)
    do t = 1, sim%nt
        v1 = v0 + sim%u(:, a2, t) - sim%u(:, a1, t)
        trajectory(:, t) = v1
    end do

end subroutine

subroutine put_shape_on_new_mesh(x,y,xi,yi)
    real(r8), dimension(:), intent(in) :: x,y
    real(r8), dimension(:), intent(in) :: xi
    real(r8), dimension(:), intent(out) :: yi
    !
    integer :: i,j,n,ni
    integer :: ii,jj
    real(r8) :: f0,sigma,n0,n1
    real(r8) :: ei,eimin,eimax,eirange,foursigma
    !
    n=size(x,1)
    ni=size(xi,1)
    yi=0.0_r8
    sigma=x(2)-x(1)
    !
    n0=lo_trapezoid_integration(x,y)
    eimin=minval(x)
    eimax=maxval(x)
    eirange=1.0_r8/(eimax-eimin)
    foursigma=4*sigma
    do i=1,ni
        ! find lower and upper bound
        ei=xi(i)
        ii=floor(n*(ei-foursigma-eimin)*eirange)
        ii=max(ii,1)
        ii=min(ii,n)
        !
        jj=ceiling(n*(ei+foursigma-eimin)*eirange)
        jj=min(jj,n)
        jj=max(jj,1)
        do j=ii,jj
            f0=lo_gauss(xi(i),x(j),sigma)
            yi(i)=yi(i)+f0*y(j)
        enddo
    enddo
    n1=lo_trapezoid_integration(xi,yi)
    yi=yi*n0/n1
end subroutine

end module
