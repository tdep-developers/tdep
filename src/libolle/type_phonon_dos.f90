module type_phonon_dos
!!
!! Get the phonon density of states, as well as the stuff you can integrate from it!
!! In addition, we can also get other quantities frequency-resolved and projected, while
!! we are doing it.
!!
use konstanter, only: r8, i8, lo_huge, lo_hugeint, lo_freqtol, lo_status, lo_twopi, lo_pi, lo_tol, lo_exitcode_param, &
                      lo_gnuplotterminal, lo_frequency_hartree_to_icm, lo_frequency_hartree_to_thz, &
                      lo_frequency_hartree_to_mev, lo_sqtol
use gottochblandat, only: open_file, walltime, lo_gauss, lo_trueNtimes, lo_trapezoid_integration, tochar, lo_linspace, &
                          lo_progressbar_init, lo_progressbar, lo_chop, lo_mean
use mpi_wrappers, only: lo_mpi_helper, lo_stop_gracefully
use lo_memtracker, only: lo_mem_helper
use hdf5_wrappers, only: lo_hdf5_helper
use dump_data, only: lo_dump_palette_to_gnuplot
use type_crystalstructure, only: lo_crystalstructure
use type_phonon_dispersions, only: lo_phonon_dispersions
use type_qpointmesh, only: lo_qpoint_mesh, lo_integration_weights_for_one_tetrahedron, lo_LV_tetrahedron_weights

implicit none
private
public :: lo_phonon_dos

!> Phonon density of states
type lo_phonon_dos
    !> number of dos-points
    integer :: n_dos_point = -lo_hugeint
    !> number of atoms
    integer :: n_atom = -lo_hugeint
    !> number of modes
    integer :: n_mode = -lo_hugeint
    !> the frequency axis
    real(r8), dimension(:), allocatable :: omega
    !> the DOS
    real(r8), dimension(:), allocatable :: dos
    !> site projected density of states
    real(r8), dimension(:, :), allocatable :: pdos_site
    !> mode projected density of states
    real(r8), dimension(:, :), allocatable :: pdos_mode
    !> maximum frequency
    real(r8) :: dosmax = -lo_huge
    !> minimum frequency
    real(r8) :: dosmin = -lo_huge
    !> which integration scheme to use
    integer :: integrationtype = -lo_hugeint
    !> Should I add a prefactor to the default smearing
    real(r8) :: smearing_prefactor = -lo_huge
contains
    !> Calculate the DOS
    procedure :: generate
    !> Calculate spectral thermal conductivity
    procedure :: spectral_kappa
    !> Calculate spectral angular momentum
    procedure :: spectral_angular_momentum
    !> dump it to file
    procedure :: write_to_file
    !> dump it to hdf5
    procedure :: write_to_hdf5
    !> size in memory
    procedure :: size_in_mem => phonon_dos_size_in_mem
    !> destroy the phonon dos
    procedure :: destroy
end type

! interfaces to type_phonon_dos_integrations
interface
    module subroutine lo_get_phonon_dos_tetrahedron(pd, qp, dr, mw, mem, verbosity)
        type(lo_phonon_dos), intent(inout) :: pd
        type(lo_qpoint_mesh), intent(in) :: qp
        type(lo_phonon_dispersions), intent(in) :: dr
        type(lo_mpi_helper), intent(inout) :: mw
        type(lo_mem_helper), intent(inout) :: mem
        integer, intent(in) :: verbosity
    end subroutine
    module subroutine lo_get_phonon_dos_gaussian(pd, qp, dr, mw, verbosity)
        type(lo_phonon_dos), intent(inout) :: pd
        class(lo_qpoint_mesh), intent(in) :: qp
        type(lo_phonon_dispersions), intent(in) :: dr
        type(lo_mpi_helper), intent(inout) :: mw
        integer, intent(in) :: verbosity
    end subroutine
    module subroutine spectral_kappa(pd, uc, qp, dr, mw, mem, spec_kappa, spec_kappa_band, spec_kappa_atom)
        class(lo_phonon_dos), intent(inout) :: pd
        type(lo_crystalstructure), intent(in) :: uc
        type(lo_qpoint_mesh), intent(in) :: qp
        type(lo_phonon_dispersions), intent(in) :: dr
        type(lo_mpi_helper), intent(inout) :: mw
        type(lo_mem_helper), intent(inout) :: mem
        real(r8), dimension(:, :), allocatable, intent(out) :: spec_kappa
        real(r8), dimension(:, :, :), allocatable, intent(out) :: spec_kappa_band
        real(r8), dimension(:, :, :), allocatable, intent(out) :: spec_kappa_atom
    end subroutine
    module subroutine spectral_angular_momentum(pd, uc, qp, dr, temperature, mw, mem, spec_angmom, spec_angmom_band, spec_angmom_atom)
        class(lo_phonon_dos), intent(inout) :: pd
        type(lo_crystalstructure), intent(in) :: uc
        type(lo_qpoint_mesh), intent(in) :: qp
        type(lo_phonon_dispersions), intent(in) :: dr
        real(r8), intent(in) :: temperature
        type(lo_mpi_helper), intent(inout) :: mw
        type(lo_mem_helper), intent(inout) :: mem
        real(r8), dimension(:, :), allocatable, intent(out) :: spec_angmom
        real(r8), dimension(:, :, :), allocatable, intent(out) :: spec_angmom_band
        real(r8), dimension(:, :, :), allocatable, intent(out) :: spec_angmom_atom
    end subroutine
end interface

contains

!> Calculate the phonon density of states
subroutine generate(pd, dr, qp, uc, mw, mem, verbosity, sigma, n_dos_point, integrationtype)
    !> The phonon dos
    class(lo_phonon_dos), intent(inout) :: pd
    !> the dispersion relations
    type(lo_phonon_dispersions), intent(in) :: dr
    !> the q-point grid
    class(lo_qpoint_mesh), intent(inout) :: qp
    !> the crystal structure
    type(lo_crystalstructure), intent(in) :: uc
    !> MPI helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem
    !> how much to talk
    integer, intent(in) :: verbosity
    !> Non-default smearing parameter
    real(r8), intent(in), optional :: sigma
    !> change the default number of points on the x-axis
    integer, intent(in), optional :: n_dos_point
    !> set the integration type
    integer, intent(in), optional :: integrationtype

    real(r8) :: timer, t0, t1

    ! Set some basic things, defaults and such
    init: block
        integer :: i

        if (verbosity .gt. 0) then
            timer = walltime()
            t0 = timer
            t1 = 0.0_r8
            write (*, *) ''
            write (*, *) 'Calculating phonon density of states'
        end if

        ! Sort out the integration type
        if (present(integrationtype)) then
            pd%integrationtype = integrationtype
        else
            pd%integrationtype = 2
        end if

        ! Sort out smearing
        if (present(sigma)) then
            pd%smearing_prefactor = sigma
        else
            pd%smearing_prefactor = 1.0_r8
        end if

        ! Number of points in the density of states.
        if (present(n_dos_point)) then
            pd%n_dos_point = n_dos_point
        else
            pd%n_dos_point = 1000
        end if

        ! Choose the smallest frequency so that imaginary things will show
        pd%dosmin = lo_huge
        do i = 1, dr%n_irr_qpoint
            pd%dosmin = min(pd%dosmin, minval(dr%iq(i)%omega))
        end do
        if (pd%dosmin .gt. -lo_freqtol) pd%dosmin = 0.0_r8
        ! and the largest frequency is quite straightforward
        pd%dosmax = dr%omega_max + 4.01_r8*maxval(dr%default_smearing)

        ! use the min and max to get the frequency axis
        allocate (pd%omega(pd%n_dos_point))
        call lo_linspace(pd%dosmin, pd%dosmax, pd%omega)
        ! and space for the actual density of states
        pd%n_mode = dr%n_mode
        pd%n_atom = uc%na
        allocate (pd%dos(pd%n_dos_point))
        allocate (pd%pdos_site(pd%n_dos_point, pd%n_atom))
        allocate (pd%pdos_mode(pd%n_dos_point, pd%n_mode))
        pd%dos = 0.0_r8
        pd%pdos_site = 0.0_r8
        pd%pdos_mode = 0.0_r8
    end block init

    ! Call the actual integration to get a raw, not so pretty DOS
    integrate: block
        select case (pd%integrationtype)
        case (1)
            if (verbosity .gt. 0) then

                write (*, "(1X,A)") '... fix gaussian integration, sigma = '//tochar(lo_mean(dr%default_smearing)*pd%smearing_prefactor*lo_Frequency_Hartree_to_THz)//' THz, N energy points: '//tochar(pd%n_dos_point)
            end if
            call lo_get_phonon_dos_gaussian(pd, qp, dr, mw, verbosity)
        case (2)
            if (verbosity .gt. 0) then
                write (*, "(1X,A)") '... adaptive gaussian integration, scalingfactor = '//tochar(pd%smearing_prefactor)//', N energy points: '//tochar(pd%n_dos_point)
            end if
            call lo_get_phonon_dos_gaussian(pd, qp, dr, mw, verbosity)
        case (3)
            if (verbosity .gt. 0) then
                write (*, "(1X,A)") '... tetrahedron integration, N energy points: '//tochar(pd%n_dos_point)
            end if
            call lo_get_phonon_dos_tetrahedron(pd, qp, dr, mw, mem, verbosity)
        case default
            write (*, *) 'unknown integration type'
            stop
        end select
    end block integrate

    ! Now fix the DOS so that it looks nice.
    makenice: block
        integer :: i, j, k
        real(r8), dimension(:, :), allocatable :: sitebuf
        real(r8) :: f0, f1, dostol

        ! I might have gotten a contribution at zero frequency due to the smearing. That has to
        ! be removed. I subtract the line that goes from the dos at 0 to the max frequency.
        do j = 1, pd%n_mode
            f1 = pd%pdos_mode(1, j)
            do i = 1, pd%n_dos_point
                f0 = (pd%n_dos_point - i)*1.0_r8/((pd%n_dos_point - 1)*1.0_r8)
                f0 = f0**2
                pd%pdos_mode(i, j) = pd%pdos_mode(i, j) - f0*f1
                if (pd%pdos_mode(i, j) .lt. lo_freqtol) pd%pdos_mode(i, j) = 0.0_r8
            end do
            pd%pdos_mode(1, j) = 0.0_r8
        end do
        do j = 1, pd%n_atom
            f1 = pd%pdos_site(1, j)
            do i = 1, pd%n_dos_point
                f0 = (pd%n_dos_point - i)*1.0_r8/((pd%n_dos_point - 1)*1.0_r8)
                f0 = f0**2
                pd%pdos_site(i, j) = pd%pdos_site(i, j) - f0*f1
                if (pd%pdos_site(i, j) .lt. lo_freqtol) pd%pdos_site(i, j) = 0.0_r8
            end do
            pd%pdos_site(1, j) = 0.0_r8
        end do

        ! Sum up contributions and normalize things
        pd%dos = 0.0_r8
        do j = 1, pd%n_mode
            f0 = 1.0_r8/lo_trapezoid_integration(pd%omega, pd%pdos_mode(:, j))
            pd%dos = pd%dos + pd%pdos_mode(:, j)
            pd%pdos_mode(:, j) = pd%pdos_mode(:, j)*f0
        end do
        pd%dos(1) = 0.0_r8
        ! Get a tolerance where to slice off the phonon dos.
        dostol = lo_mean(pd%dos)*lo_sqtol
        pd%dos = lo_chop(pd%dos, dostol)
        ! Normalize the total
        f0 = lo_trapezoid_integration(pd%omega, pd%dos)
        if (verbosity .gt. 0) write (*, *) '... raw normalization (should be 1):', f0/dr%n_mode
        pd%dos = pd%dos*dr%n_mode/f0
        ! Adjust the mode projected so that they sum up to the total
        do i = 1, pd%n_dos_point
            f0 = sum(pd%pdos_mode(i, :))
            if (f0 .gt. dostol) then
                pd%pdos_mode(i, :) = pd%pdos_mode(i, :)*pd%dos(i)/f0
            else
                pd%pdos_mode(i, :) = 0.0_r8
            end if
        end do

        ! Fix the degeneracy of the site-projected
        call mem%allocate(sitebuf, [pd%n_dos_point, pd%n_atom], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        sitebuf = 0.0_r8
        do i = 1, pd%n_atom
            ! enfore site degeneracy
            do j = 1, uc%sym%degeneracy(i)
                k = uc%sym%degenerate_atom(j, i)
                sitebuf(:, i) = sitebuf(:, i) + pd%pdos_site(:, k)/(1.0_r8*uc%sym%degeneracy(j))
            end do
            f0 = lo_trapezoid_integration(pd%omega, sitebuf(:, i))
            sitebuf(:, i) = sitebuf(:, i)*3.0_r8/f0
        end do
        pd%pdos_site = sitebuf
        call mem%deallocate(sitebuf, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        ! normalize the projections so that they add up to the total
        do i = 1, pd%n_dos_point
            f0 = sum(pd%pdos_site(i, :))
            if (f0 .gt. dostol) then
                pd%pdos_site(i, :) = pd%pdos_site(i, :)*pd%dos(i)/f0
            else
                pd%pdos_site(i, :) = 0.0_r8
            end if
        end do
    end block makenice

    if (verbosity .gt. 0) then
        write (*, *) 'Finished calculating the phonon density of states (', tochar(walltime() - timer), 's)'
    end if
end subroutine

!> Print the phonon dos to hdf5
subroutine write_to_hdf5(pd, p, enhet, filename, mem, hdftag)
    !> phonon dos
    class(lo_phonon_dos), intent(in) :: pd
    !> crystal structure
    type(lo_crystalstructure), intent(in) :: p
    !> unit, "thz", "mev" or "icm"
    character(len=3), intent(in) :: enhet
    !> filename
    character(len=*), intent(in) :: filename
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem
    !> optionally, write to an already open hdf5 file.
    integer(i8), intent(in), optional :: hdftag

    type(lo_hdf5_helper) :: h5
    real(r8), dimension(:, :), allocatable :: buf
    real(r8) :: unitfactor
    character(len=100) :: omstr, dosstr
    character(len=2000) :: atomnames
    integer :: a1, a2, i

    ! The unit thingy
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

    ! create the file, or use a provided tag
    if (present(hdftag)) then
        ! Write to the tag that was provided
        h5%file_id = hdftag
    else
        ! Create a new file.
        call h5%init(__FILE__, __LINE__)
        call h5%open_file('write', trim(filename))
    end if

    ! store all the data
    call h5%store_data(pd%omega*unitfactor, h5%file_id, 'frequencies', enhet=trim(omstr), dimensions='frequency')
    call h5%store_data(pd%dos/unitfactor, h5%file_id, 'dos', enhet=trim(dosstr), dimensions='dos')
    call h5%store_data(pd%pdos_site/unitfactor, h5%file_id, 'dos_per_site', enhet=trim(dosstr), dimensions='site,dos')
    call h5%store_data(pd%pdos_mode/unitfactor, h5%file_id, 'dos_per_mode', enhet=trim(dosstr), dimensions='mode,dos')

    ! And store the atom labels
    call p%unique_atom_label(atomnames)
    call h5%store_attribute(trim(adjustl(atomnames)), h5%file_id, 'unique_atom_labels')

    ! Store the density of states per unique atom
    call mem%allocate(buf, [pd%n_dos_point, p%sym%n_irreducible_atom], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    buf = 0.0_r8
    do a1 = 1, p%sym%n_irreducible_atom
        do i = 1, p%sym%irr_unfold_ctr(a1)
            a2 = p%sym%irr_unfold_index(i, a1)
            buf(:, a1) = buf(:, a1) + pd%pdos_site(:, a2)
        end do
    end do
    buf = buf/unitfactor
    call h5%store_data(buf, h5%file_id, 'dos_per_unique_atom', enhet=trim(dosstr), dimensions='unique atom,dos')
    call mem%deallocate(buf, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)

    ! close the file and hdf5, if relevant.
    if (present(hdftag) .eqv. .false.) then
        call h5%close_file()
        call h5%destroy(__FILE__, __LINE__)
    end if
end subroutine

!> Writes phonon dos to file
subroutine write_to_file(pd, uc, enhet, filename)
    !> The phonon dos
    class(lo_phonon_dos), intent(in) :: pd
    !> The crystal structure
    type(lo_crystalstructure), intent(in) :: uc
    !> The unit, "thz", "mev" or "icm"
    character(len=3), intent(in) :: enhet
    !> optionally a different filename
    character(len=*), intent(in) :: filename
    !
    integer :: i, j, u
    real(r8), dimension(:), allocatable :: y
    real(r8) :: unitfactor
    character(len=100) :: opf, lgd

    ! Figure out the output format
    i = 2 + size(pd%pdos_site, 2) + size(pd%pdos_mode, 2)
    opf = "("//tochar(i)//"(1X,E18.12))"
    allocate (y(i))
    y = 0.0_r8

    ! The unit thingy
    select case (enhet)
    case ('thz')
        unitfactor = lo_frequency_hartree_to_THz
    case ('mev')
        unitfactor = lo_frequency_hartree_to_meV
    case ('icm')
        unitfactor = lo_frequency_hartree_to_icm
    case default
        call lo_stop_gracefully(['Unknown unit'], lo_exitcode_param, __FILE__, __LINE__)
    end select

    ! Print the raw data
    u = open_file('out', trim(filename))
    do i = 1, pd%n_dos_point
        y(1) = pd%omega(i)*unitfactor
        y(2) = pd%dos(i)/unitfactor
        do j = 1, size(pd%pdos_site, 2)
            y(2 + j) = pd%pdos_site(i, j)/unitfactor
        end do
        do j = 1, size(pd%pdos_mode, 2)
            y(2 + size(pd%pdos_site, 2) + j) = pd%pdos_mode(i, j)/unitfactor
        end do
        write (u, opf) y
    end do
    close (u)
    deallocate (y)

    ! Get the gnuplot file
    u = open_file('out', trim(filename)//'.gnuplot')
    ! Choose terminal
    write (u, *) 'set terminal '//lo_gnuplotterminal//' size 800,400 enhanced font "CMU Serif,8"'
    write (u, *) 'set multiplot'
    write (u, *) 'set origin 0.0, 0.0'
    write (u, *) 'set size 0.5, 1.0'
    write (u, *) 'set border lw 0.5'
    write (u, *) 'set ytics scale 0.5'
    write (u, *) 'set xtics scale 0.5'
    write (u, *) 'set mytics 10'
    write (u, *) 'set mxtics 10'
    select case (enhet)
    case ('thz')
        write (u, *) ' set xlabel "Frequency (THz)" '
    case ('mev')
        write (u, *) ' set xlabel "Frequency (meV)" '
    case ('icm')
        write (u, *) ' set xlabel "Frequency (cm^-^1)" '
    end select
    write (u, *) ' set ylabel "DOS" '
    write (u, *) ' set title "Site-projected DOS" '
    write (u, *) 'set key top left'
    write (u, '(A)', advance='no') 'plot'
    write (u, '(A)') ' "'//trim(filename)//'" u 1:2 w line lc rgb "#318712" t "Total", \'
    do i = 1, pd%n_atom
        lgd = 'site '//tochar(i)//', '//uc%atomic_symbol(uc%species(i))
        write (u, '(A)', advance='no') ' "'//trim(filename)//'" u 1:'//tochar(i + 2)//' w line'
        write (u, '(A)', advance='no') ' t "'//trim(lgd)//'"'
        if (i .lt. pd%n_atom) then
            write (u, '(A)') ',\'
        end if
    end do
    write (u, *) ' '
    write (u, *) 'set origin 0.5, 0.0'
    write (u, *) 'set size 0.5, 1.0'
    write (u, *) 'set border lw 0.5'
    write (u, *) 'set ytics scale 0.5'
    write (u, *) 'set xtics scale 0.5'
    write (u, *) 'set mytics 10'
    write (u, *) 'set mxtics 10'
    select case (enhet)
    case ('thz')
        write (u, *) ' set xlabel "Frequency (THz)" '
    case ('mev')
        write (u, *) ' set xlabel "Frequency (meV)" '
    case ('icm')
        write (u, *) ' set xlabel "Frequency (cm^-^1)" '
    end select
    write (u, *) ' set ylabel "DOS" '
    write (u, *) ' set title "Mode-projected DOS" '
    write (u, *) 'set key top left'
    write (u, '(A)', advance='no') 'plot'
    write (u, '(A)') ' "'//trim(filename)//'" u 1:2 w line lc rgb "#318712" t "Total", \'
    do i = 1, size(pd%pdos_mode, 2)
        write (u, '(A)', advance='no') ' "'//trim(filename)//'" u 1:'//tochar(i + 2 + size(pd%pdos_site, 2))//' w line'
        write (u, '(A)', advance='no') ' t ""'
        if (i .lt. size(pd%pdos_mode, 2)) then
            write (u, '(A)') ',\'
        end if
    end do
    write (u, *) ' '
    call lo_dump_palette_to_gnuplot(u)
    close (u)
end subroutine

!> measure size in memory, in bytes
function phonon_dos_size_in_mem(d) result(mem)
    !> phonon dos
    class(lo_phonon_dos), intent(in) :: d
    !> memory in bytes
    integer :: mem

    mem = 0
    mem = mem + storage_size(d)
    if (allocated(d%omega)) mem = mem + storage_size(d%omega)*size(d%omega)
    if (allocated(d%dos)) mem = mem + storage_size(d%dos)*size(d%dos)
    if (allocated(d%pdos_site)) mem = mem + storage_size(d%pdos_site)*size(d%pdos_site)
    if (allocated(d%pdos_mode)) mem = mem + storage_size(d%pdos_mode)*size(d%pdos_mode)
    mem = mem/8
end function

!> destroy everything
subroutine destroy(pd)
    !> phonon dos
    class(lo_phonon_dos), intent(inout) :: pd

    if (allocated(pd%omega)) deallocate (pd%omega)
    if (allocated(pd%dos)) deallocate (pd%dos)
    if (allocated(pd%pdos_site)) deallocate (pd%pdos_site)
    if (allocated(pd%pdos_mode)) deallocate (pd%pdos_mode)
    pd%n_dos_point = -lo_hugeint
    pd%n_atom = -lo_hugeint
    pd%n_mode = -lo_hugeint
    pd%dosmax = -lo_huge
    pd%dosmin = -lo_huge
    pd%integrationtype = -lo_hugeint
    pd%smearing_prefactor = -lo_huge
end subroutine

end module
