module io
!! Dump things to file, neatly
use konstanter, only: r8, lo_iou, lo_status, lo_sqtol, lo_freqtol, lo_frequency_Hartree_to_meV, lo_exitcode_param, &
                      lo_frequency_Hartree_to_THz, lo_frequency_Hartree_to_icm, lo_exitcode_symmetry, lo_twopi, &
                      lo_pi, lo_bohr_to_A
use gottochblandat, only: lo_chop, lo_trace, lo_planck, lo_linspace, walltime, tochar
use geometryfunctions, only: lo_plane
use lo_memtracker, only: lo_mem_helper
use type_crystalstructure, only: lo_crystalstructure
use lo_dielectric_interaction, only: lo_dielectric_tensor
use mpi_wrappers, only: lo_stop_gracefully
use hdf5_wrappers, only: lo_hdf5_helper, lo_h5_store_data, lo_h5_store_attribute, HID_T, H5F_ACC_TRUNC_F, h5close_f, h5open_f, h5fclose_f, h5fopen_f, h5fcreate_f
use lo_phonon_bandstructure_on_path, only: lo_phonon_bandstructure
use phonondamping, only: lo_phonon_selfenergy !,getintensity
use dielscatter, only: lo_dielectric_response
use type_qpointmesh, only: lo_qpoint, lo_qpoint_mesh
use type_phonon_dispersions, only: lo_phonon_dispersions_qpoint
use type_blas_lapack_wrappers, only: lo_gemm

use lineshape_helper, only: evaluate_spectral_function
implicit none

private
public :: write_lineshape_to_hdf5

contains

!> dump the lineshape to hdf5
subroutine write_lineshape_to_hdf5(se, p, qpt, wp, dir, di, qp, enhet, temperature, filename, mem, verbosity, dielectric)
    !> self-energy
    type(lo_phonon_selfenergy), intent(in) :: se
    !> structure
    type(lo_crystalstructure), intent(in) :: p
    !> q-point
    type(lo_qpoint), intent(in) :: qpt
    !> harmonic properties
    type(lo_phonon_dispersions_qpoint), intent(in) :: wp
    !> dielectric response
    type(lo_dielectric_response), intent(in) :: dir
    !> dielectric interactions
    type(lo_dielectric_tensor), intent(in) :: di
    !> q-point mesh
    class(lo_qpoint_mesh), intent(in) :: qp
    !> the unit
    character(len=3), intent(in) :: enhet
    !> temperature
    real(r8), intent(in) :: temperature
    !> filename
    character(len=*), intent(in) :: filename
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem
    !> talk a lot?
    integer, intent(in) :: verbosity
    !> write dielectric things
    logical, intent(in) :: dielectric

    type(lo_hdf5_helper) :: h5
    real(r8), dimension(:, :), allocatable :: proj
    real(r8) :: unitfactor, t0, t1, timer
    character(len=10) :: unitname

    ! Start timer
    timer = walltime()
    t0 = timer
    t1 = timer

    ! Set some basic things
    init: block
        complex(r8), dimension(:, :), allocatable :: omega, dynmat, wA
        complex(r8), dimension(3) :: v0
        integer :: i, j

        ! Get the units right
        select case (enhet)
        case ('thz')
            unitfactor = lo_frequency_Hartree_to_THz
            unitname = 'Thz'
        case ('icm')
            unitfactor = lo_frequency_Hartree_to_icm
            unitname = '1/cm'
        case ('mev')
            unitfactor = lo_frequency_Hartree_to_meV
            unitname = 'meV'
        case default
            call lo_stop_gracefully(["Unknown unit, try 'thz', 'mev' or 'icm'"], lo_exitcode_param, __FILE__, __LINE__)
        end select

        ! Create a new file.
        call h5%init(__FILE__, __LINE__)
        call h5%open_file('write', trim(filename))

        ! Write some basic things that should always be present, first the temperature
        call h5%store_attribute(temperature, h5%file_id, 'temperature')
        ! which q-point we are on
        call h5%store_data(qpt%r, h5%file_id, 'q-point', enhet='1/A')
        ! what q-mesh did I use? Store the metadata in the file, just to be on the safe side.
        call h5%open_group('write', 'qmesh')
        call qp%write_metadata_to_hdf5(input_id=h5%group_id)
        call h5%close_group()
        ! does not hurt to store the structure as well
        call h5%open_group('write', 'structure')
        call p%write_structure_to_hdf5(input_id=h5%group_id)
        call h5%close_group()

        ! store incident wavecetor
        call h5%store_data(se%qdir, h5%file_id, 'incident_wavevector', enhet='Cartesian')

        ! Reconstruct dynamical matrix? Seems like a sensible thing to do.
        call mem%allocate(omega, [se%n_mode, se%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(dynmat, [se%n_mode, se%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(wA, [se%n_mode, se%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        omega = 0.0_r8
        dynmat = 0.0_r8
        wA = 0.0_r8
        do i = 1, se%n_mode
            omega(i, i) = (wp%omega(i)*unitfactor)**2
        end do
        call lo_gemm(wp%egv, omega, wA)
        call lo_gemm(wA, wp%egv, dynmat, transb='C')

        ! Seems sensible. Store the Harmonic things.
        call h5%open_group('write', 'harmonic')

        call h5%store_data(wp%omega*unitfactor, h5%group_id, 'harmonic_frequencies', enhet=trim(unitname))
        call h5%store_data(real(dynmat), h5%group_id, 'dynamical_matrix_re')
        call h5%store_data(aimag(dynmat), h5%group_id, 'dynamical_matrix_im')
        call h5%store_data(real(wp%egv), h5%group_id, 'eigenvectors_re')
        call h5%store_data(aimag(wp%egv), h5%group_id, 'eigenvectors_im')

        call mem%deallocate(omega, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(dynmat, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(wA, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)

        ! Get the array that can project onto sites
        call mem%allocate(proj, [p%na, se%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        proj = 0.0_r8
        do j = 1, se%n_mode
            do i = 1, p%na
                v0 = wp%egv((i - 1)*3 + 1:i*3, j)
                proj(i, j) = real(dot_product(v0, v0), r8)
            end do
            proj(:, j) = proj(:, j)/sum(abs(proj(:, j)))
            proj(:, j) = lo_chop(proj(:, j), 1E-12_r8)
        end do

        if (verbosity .gt. 0) then
            t1 = walltime()
            write (lo_iou, *) '... wrote harmonic properties (', tochar(t1 - t0), 's)'
        end if

        call h5%close_group()
    end block init

    ! Store phonon self-energy
    phononphonon: block
        real(r8), dimension(:, :), allocatable :: dr0, dr1, dr2
        integer :: i, j, k

        call h5%open_group('write', 'anharmonic')

        call h5%store_data(se%energy_axis*unitfactor, h5%group_id, 'frequency', enhet=trim(unitname))
        call h5%store_data(se%im_iso*unitfactor, h5%group_id, 'imaginary_isotope_selfenergy', enhet=trim(unitname))
        call h5%store_data(se%im_3ph*unitfactor, h5%group_id, 'imaginary_threephonon_selfenergy', enhet=trim(unitname))
        call h5%store_data(se%re_3ph*unitfactor, h5%group_id, 'real_threephonon_selfenergy', enhet=trim(unitname))
        call h5%store_data(se%re_4ph*unitfactor, h5%group_id, 'real_fourphonon_selfenergy', enhet=trim(unitname))

        call mem%allocate(dr0, [se%n_energy, se%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(dr1, [se%n_energy, p%na], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(dr2, [se%n_energy, p%sym%n_irreducible_atom], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)

        dr0 = 0.0_r8
        dr1 = 0.0_r8
        dr2 = 0.0_r8
        ! Imaginary part of Green's function
        do i = 1, se%n_mode
            if (wp%omega(i) .gt. lo_freqtol*10) then
                call evaluate_spectral_function(se%energy_axis, se%im_3ph(:, i) + se%im_iso(:, i), se%re_3ph(:, i) + se%re_4ph(:, i), wp%omega(i), dr0(:, i))
            else
                dr0(:, i) = 0.0_r8
            end if
        end do
        ! project onto sites
        do j = 1, se%n_mode
        do i = 1, p%na
            dr1(:, i) = dr1(:, i) + dr0(:, j)*proj(i, j)
        end do
        end do
        ! project onto unique atoms
        do i = 1, p%sym%n_irreducible_atom
        do j = 1, p%sym%irr_unfold_ctr(i)
            k = p%sym%irr_unfold_index(j, i)
            dr2(:, i) = dr2(:, i) + dr1(:, k)
        end do
        end do
        ! Store the different versions of the intensities
        dr0 = dr0/unitfactor
        dr1 = dr1/unitfactor
        dr2 = dr2/unitfactor
        call h5%store_data(dr0, h5%group_id, 'spectralfunction_per_mode', enhet='states/'//trim(unitname), dimensions='mode,frequency')
        call h5%store_data(dr1, h5%group_id, 'spectralfunction_per_site', enhet='states/'//trim(unitname), dimensions='site,frequency')
        call h5%store_data(dr2, h5%group_id, 'spectralfunction_per_unique_atom', enhet='states/'//trim(unitname), dimensions='atom,frequency')

        ! and cleanup
        call mem%deallocate(dr0, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(dr1, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(dr2, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)

        call h5%close_group()

        if (verbosity .gt. 0) then
            t1 = walltime()
            write (lo_iou, *) '... wrote anharmonic properties (', tochar(t1 - t0), 's)'
            t0 = t1
        end if
    end block phononphonon

    ! Store IR things
    if (dielectric) then
        infrared: block
            complex(r8), dimension(3, 3) :: m0
            complex(r8), dimension(:, :, :), allocatable :: dr0
            real(r8), dimension(:, :, :), allocatable :: dr1
            integer :: i, j

            call h5%open_group('write', 'susceptibility')

            ! Start with the baseline eps_infinity and the thermal correction
            call h5%store_data(dir%eps_inf, h5%group_id, 'epsilon_infinity', enhet='dimensionless')
            call h5%store_data(dir%eps_inf_corr, h5%group_id, 'epsilon_infinity_thermal', enhet='dimensionless')

            ! Born charges, makes sense to have them stashed away

            if (di%n_Z_singlet .gt. 0) then
                call h5%store_data(di%Z_singlet, h5%group_id, 'born_effective_charges', enhet='e/A', dimensions='Exyz,Rxyz,atom')
            end if

            call mem%allocate(dr0, [se%n_energy, 3, 3], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%allocate(dr1, [se%n_energy, 3, 3], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            dr0 = 0.0_r8
            dr1 = 0.0_r8
            ! Ok, there will be a ton of terms here. Not sure how to do it right.
            ! note, I have no idea what the units should be.
            call h5%store_data(se%energy_axis*unitfactor, h5%group_id, 'frequency', enhet=trim(unitname))

            ! first the zeroth order response, just for reference to see how bad it is.
            dr0 = 0.0_r8
            do i = 1, se%n_mode
            do j = 1, se%n_energy
                m0 = dir%chi_tensor_zero(:, :, j, i)
                dr0(j, :, :) = dr0(j, :, :) + m0
            end do
            end do

            dr1 = real(dr0, r8)
            call h5%store_data(dr1, h5%group_id, 'chi_0_re', enhet='dunno', dimensions='xyz,xyz,frequency')
            dr1 = aimag(dr0)
            call h5%store_data(dr1, h5%group_id, 'chi_0_im', enhet='dunno', dimensions='xyz,xyz,frequency')

            ! first order with broadening
            dr0 = 0.0_r8
            !dr1=0.0_r8
            do i = 1, se%n_mode
            do j = 1, se%n_energy
                m0 = dir%chi_tensor_one(:, :, j, i)
                dr0(j, :, :) = dr0(j, :, :) + m0
            end do
            end do
            ! this extra step of putting the real and imaginary parts in another buffer was needed
            ! to avoid segfaults for tight grids + gcc combination. Go figure.
            dr1 = real(dr0, r8)
            call h5%store_data(dr1, h5%group_id, 'chi_1_re', enhet='dunno', dimensions='xyz,xyz,frequency')
            dr1 = aimag(dr0)
            call h5%store_data(dr1, h5%group_id, 'chi_1_im', enhet='dunno', dimensions='xyz,xyz,frequency')

            ! first order spherically averaged
            dr0 = 0.0_r8
            do i = 1, se%n_mode
            do j = 1, se%n_energy
                m0 = dir%chi_tensor_one_sph(:, :, j, i)
                dr0(j, :, :) = dr0(j, :, :) + m0
            end do
            end do
            dr1 = real(dr0, r8)
            call h5%store_data(dr1, h5%group_id, 'chi_1_sph_re', enhet='dunno', dimensions='xyz,xyz,frequency')
            dr1 = aimag(dr0)
            call h5%store_data(dr1, h5%group_id, 'chi_1_sph_im', enhet='dunno', dimensions='xyz,xyz,frequency')

            ! second order
            dr0 = 0.0_r8
            do j = 1, se%n_energy
                m0 = dir%chi_tensor_two(:, :, j)
                dr0(j, :, :) = dr0(j, :, :) + m0
            end do
            dr1 = real(dr0, r8)
            call h5%store_data(dr1, h5%group_id, 'chi_2_re', enhet='dunno', dimensions='xyz,xyz,frequency')
            dr1 = aimag(dr0)
            call h5%store_data(dr1, h5%group_id, 'chi_2_im', enhet='dunno', dimensions='xyz,xyz,frequency')

            dr0 = 0.0_r8
            do j = 1, se%n_energy
                m0 = dir%chi_tensor_three(:, :, j)
                dr0(j, :, :) = dr0(j, :, :) + m0
            end do
            dr1 = real(dr0, r8)
            call h5%store_data(dr1, h5%group_id, 'chi_3_re', enhet='dunno', dimensions='xyz,xyz,frequency')
            dr1 = aimag(dr0)
            call h5%store_data(dr1, h5%group_id, 'chi_3_im', enhet='dunno', dimensions='xyz,xyz,frequency')

            ! cleanup
            call mem%deallocate(dr0, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%deallocate(dr1, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)

            call h5%close_group()

            if (verbosity .gt. 0) then
                t1 = walltime()
                write (lo_iou, *) '... wrote susceptibility (', tochar(t1 - t0), 's)'
                t0 = t1
            end if
        end block infrared
    end if

    if (dielectric) then
        raman: block
            integer, parameter :: nangle = 360
            type(lo_plane) :: pl_in
            real(r8), dimension(:, :, :), allocatable :: dw0
            real(r8), dimension(:, :), allocatable :: du0, du1, du2, du3, du4
            real(r8), dimension(:), allocatable :: dv0
            real(r8), dimension(81) :: dv_par, dv_per, dv_integrate
            real(r8), dimension(3) :: ev1, ev2, ev3
            integer :: i, j, k, l, iang, imode, jmode, ii

            ! Create a group to store Raman things
            call h5%open_group('write', 'Raman')

            ! Vector that performs integration over all angles to get the total. Zeros and ones that select the proper elements.
            ! dv_integrate=0.0_r8
            ! dv_integrate(1)=1.0_r8
            ! dv_integrate(11)=1.0_r8
            ! dv_integrate(21)=1.0_r8
            ! dv_integrate(31)=1.0_r8
            ! dv_integrate(41)=1.0_r8
            ! dv_integrate(51)=1.0_r8
            ! dv_integrate(61)=1.0_r8
            ! dv_integrate(71)=1.0_r8
            ! dv_integrate(81)=1.0_r8
            ! On second thought, make this the || integration, perhaps. Subtle difference.
            dv_integrate = 0.0_r8
            dv_integrate(1) = 0.600000000000000_r8
            dv_integrate(5) = 0.200000000000000_r8
            dv_integrate(9) = 0.200000000000000_r8
            dv_integrate(11) = 0.200000000000000_r8
            dv_integrate(13) = 0.200000000000000_r8
            dv_integrate(21) = 0.200000000000000_r8
            dv_integrate(25) = 0.200000000000000_r8
            dv_integrate(29) = 0.200000000000000_r8
            dv_integrate(31) = 0.200000000000000_r8
            dv_integrate(37) = 0.200000000000000_r8
            dv_integrate(41) = 0.600000000000000_r8
            dv_integrate(45) = 0.200000000000000_r8
            dv_integrate(51) = 0.200000000000000_r8
            dv_integrate(53) = 0.200000000000000_r8
            dv_integrate(57) = 0.200000000000000_r8
            dv_integrate(61) = 0.200000000000000_r8
            dv_integrate(69) = 0.200000000000000_r8
            dv_integrate(71) = 0.200000000000000_r8
            dv_integrate(73) = 0.200000000000000_r8
            dv_integrate(77) = 0.200000000000000_r8
            dv_integrate(81) = 0.600000000000000_r8

            ! Some temporary space
            call mem%allocate(dv0, se%n_energy, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%allocate(du0, [se%n_energy, se%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%allocate(dw0, [se%n_energy, se%n_mode, se%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            dv0 = 0.0_r8
            du0 = 0.0_r8
            dw0 = 0.0_r8

            ! Store frequency axis
            dv0 = se%energy_axis*unitfactor
            call h5%store_data(dv0, h5%group_id, 'frequency', enhet=trim(unitname))

            do i = 1, se%n_energy
                dv0(i) = lo_planck(temperature, se%energy_axis(i))
            end do
            call h5%store_data(dv0, h5%group_id, 'thermal_prefactor', enhet='dimensionless')

            ! Store integrated Raman I, per mode
            du0 = 0.0_r8
            do imode = 1, se%n_mode
                call dgemv('T', 81, se%n_energy, 1.0_r8, dir%I_tensor_one(:, :, imode), 81, dv_integrate, 1, 1.0_r8, du0(:, imode), 1)
            end do
            call h5%store_data(du0, h5%group_id, 'intensity_I', enhet='dimensionless')

            ! Store integrated Raman II, per mode,mode
            dw0 = 0.0_r8
            do jmode = 1, se%n_mode
            do imode = 1, se%n_mode
                call dgemv('T', 81, se%n_energy, 1.0_r8, dir%I_tensor_two(:, :, imode, jmode), 81, dv_integrate, 1, 1.0_r8, dw0(:, imode, jmode), 1)
            end do
            end do
            call h5%store_data(dw0, h5%group_id, 'intensity_II', enhet='dimensionless')

            call mem%deallocate(dv0, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%deallocate(du0, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%deallocate(dw0, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)

            ! Store PO Raman
            call mem%allocate(dv0, nangle, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%allocate(du0, [81, se%n_energy], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%allocate(du1, [se%n_energy, nangle], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%allocate(du2, [se%n_energy, nangle], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%allocate(du3, [se%n_energy, nangle], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%allocate(du4, [se%n_energy, nangle], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            dv0 = 0.0_r8
            du0 = 0.0_r8
            du1 = 0.0_r8
            du2 = 0.0_r8
            du3 = 0.0_r8
            du4 = 0.0_r8

            ! Get the PO angle
            do iang = 1, nangle
                dv0(iang) = real(iang - 1, r8)*lo_twopi/360.0_r8
            end do

            ! Sum up second order over all modes
            du0 = 0.0_r8
            do jmode = 1, se%n_mode
            do imode = 1, se%n_mode
                du0 = du0 + dir%I_tensor_two(:, :, imode, jmode)
            end do
            end do

            ! DANGER ZONE. Temporarily make the rank-two-tensor approximation.
            ! dv_par=0.0_r8
            ! do i=1,3
            ! do j=1,3
            ! do k=1,3
            ! do l=1,3
            !     ii=l + (k-1)*3 + (j-1)*9 + (i-1)*27
            !     ! if ( i .eq. k ) then
            !     ! if ( j .eq. l ) then
            !     if ( i .eq. j ) then
            !     if ( k .eq. l ) then
            !         dv_par(ii)=1.0_r8
            !     endif
            !     endif
            ! enddo
            ! enddo
            ! enddo
            ! enddo
            ! do jj=1,se%n_energy
            !     f0=du0(1,jj)
            !     dv_par=0.0_r8
            !     do i=1,3
            !         ii=i + (i-1)*3 + (i-1)*9 + (i-1)*27
            !         dv_par(ii)=f0
            !         !ii=l + (k-1)*3 + (j-1)*9 + (i-1)*27
            !         ii=i + (i-1)*3 + (1-1)*9 + (1-1)*27
            !         dv_par(ii)=f0
            !         ii=i + (i-1)*3 + (2-1)*9 + (2-1)*27
            !         dv_par(ii)=f0
            !         ii=i + (i-1)*3 + (3-1)*9 + (3-1)*27
            !         dv_par(ii)=f0
            !     enddo
            !     du0(:,jj)=du0(:,jj)*dv_par
            ! enddo

            ! Generate a plane perpendicular to the incident light to get two basis vectors to rotate.
            call pl_in%generate(point=[0.0_r8, 0.0_r8, 0.0_r8], normal=se%qdir)

            do iang = 1, nangle
                ! electric field vector, incident
                ev1 = pl_in%v1*cos(dv0(iang)) + pl_in%v2*sin(dv0(iang))
                ! electric field vector, outgoing (same as above = backscattering)
                ev2 = ev1
                ! electric field vector for perpendicular polarization thing, 90 degrees angle to the vector above.
                ev3 = pl_in%v1*cos(dv0(iang) + lo_pi*0.5_r8) + pl_in%v2*sin(dv0(iang) + lo_pi*0.5_r8)
                do i = 1, 3
                do j = 1, 3
                do k = 1, 3
                do l = 1, 3
                    ii = l + (k - 1)*3 + (j - 1)*9 + (i - 1)*27
                    dv_par(ii) = ev1(i)*ev1(k)*ev2(j)*ev2(l)
                    dv_per(ii) = ev1(i)*ev1(k)*ev3(j)*ev3(l)
                end do
                end do
                end do
                end do
                ! do i=1,3
                ! do j=1,3
                ! do k=1,3
                ! do l=1,3
                !     ii=l + (k-1)*3 + (j-1)*9 + (i-1)*27
                !     dv_par(ii)=ev1(i)*ev1(j)*ev2(k)*ev2(l)
                !     dv_per(ii)=ev1(i)*ev1(j)*ev3(k)*ev3(l)
                ! enddo
                ! enddo
                ! enddo
                ! enddo

                ! Add things together with magig gemv and index algebra.
                do imode = 1, se%n_mode
                    call dgemv('T', 81, se%n_energy, 1.0_r8, dir%I_tensor_one(:, :, imode), 81, dv_par, 1, 1.0_r8, du1(:, iang), 1)
                    call dgemv('T', 81, se%n_energy, 1.0_r8, dir%I_tensor_one(:, :, imode), 81, dv_per, 1, 1.0_r8, du2(:, iang), 1)
                end do

                call dgemv('T', 81, se%n_energy, 1.0_r8, du0, 81, dv_par, 1, 1.0_r8, du3(:, iang), 1)
                call dgemv('T', 81, se%n_energy, 1.0_r8, du0, 81, dv_per, 1, 1.0_r8, du4(:, iang), 1)
            end do

            ! Write the PO stuff, first get the angles in degrees
            do iang = 1, nangle
                dv0(iang) = real(iang - 1, r8)
            end do
            call h5%store_data(dv0, h5%group_id, 'PO_angle', enhet='degree')
            call h5%store_data(du1, h5%group_id, 'PO_one_parallel', enhet='dimensionless')
            call h5%store_data(du2, h5%group_id, 'PO_one_perpendicular', enhet='dimensionless')
            call h5%store_data(du3, h5%group_id, 'PO_two_parallel', enhet='dimensionless')
            call h5%store_data(du4, h5%group_id, 'PO_two_perpendicular', enhet='dimensionless')

            call mem%deallocate(dv0, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%deallocate(du0, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%deallocate(du1, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%deallocate(du2, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%deallocate(du3, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%deallocate(du4, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)

            call h5%close_group()

            if (verbosity .gt. 0) then
                t1 = walltime()
                write (lo_iou, *) '... wrote Raman spectra (', tochar(t1 - t0), 's)'
                t0 = t1
            end if
        end block raman
    end if

    ! cleanup
    call mem%deallocate(proj, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    ! Close the file
    call h5%close_file()
    call h5%destroy(__FILE__, __LINE__)
end subroutine

end module
