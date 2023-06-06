#include "precompilerdefinitions"
module velocitydos
use konstanter, only: r8, lo_status, lo_tiny, lo_tol, lo_sqtol, lo_Hartree_to_eV, lo_bohr_to_A, lo_kb_eV, &
                      lo_huge, lo_pressure_HartreeBohr_to_GPa, lo_freqtol
use gottochblandat, only: open_file, lo_planck, lo_gauss, tochar, lo_frobnorm, lo_sqnorm, lo_linspace, &
                          lo_symmetric_eigensystem_3x3matrix, lo_trapezoid_integration, lo_outerproduct, &
                          lo_invert3x3matrix, lo_return_unique, lo_chop
use mpi_wrappers, only: lo_mpi_helper
use lo_memtracker, only: lo_mem_helper
use quadratures_stencils, only: lo_gaussianquadrature

use type_qpointmesh, only: lo_qpoint_mesh
use type_crystalstructure, only: lo_crystalstructure
use type_forceconstant_secondorder, only: lo_forceconstant_secondorder
use type_forceconstant_thirdorder, only: lo_forceconstant_thirdorder
use type_phonon_dispersions, only: lo_phonon_dispersions
use type_phonon_dos, only: lo_phonon_dos
use lo_phonon_bandstructure_on_path, only: lo_phonon_bandstructure
use type_distancetable, only: lo_distancetable
use type_mdsim, only: lo_mdsim
use hdf5_wrappers, only: HID_T, H5F_ACC_TRUNC_F, h5close_f, h5open_f, h5fclose_f, h5fopen_f, h5fcreate_f, &
                         lo_h5_store_data, lo_h5_store_attribute

implicit none
private
!public :: dump_velocity_dos
!public :: lo_kmesh_density
public :: calculate_everything
public :: lo_calculate_U0
contains

!> Calculate everything I can think of from a single-point calculation
subroutine calculate_everything(uc, bs, dr, pd, qp, fc, mw, mem)
    !> structure
    type(lo_crystalstructure), intent(inout) :: uc
    !> bandstructure
    type(lo_phonon_bandstructure), intent(in) :: bs
    !> dispersions
    type(lo_phonon_dispersions), intent(inout) :: dr
    !> phonon dos
    type(lo_phonon_dos), intent(in) :: pd
    !> q-mesh
    class(lo_qpoint_mesh), intent(in) :: qp
    !> forceconstant
    type(lo_forceconstant_secondorder), intent(inout) :: fc
    !> mpi helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem

    ! Fix parameters for all the secret things
    real(r8), dimension(4), parameter :: temperatures = [0.0_r8, 300.0_r8, 600.0_r8, 2000.0_r8]
    character(len=22) :: filename = 'outfile.topsecret.hdf5'

    ! for hdf5 io
    integer(HID_T) :: file_id

    ! Other stuff
    real(r8), dimension(4) :: da0, da1, da2, da3

    integer, parameter :: manyatoms = 50
    real(r8), dimension(:, :, :, :), allocatable :: thermal_sigma
    real(r8), dimension(:), allocatable :: dy
    real(r8), dimension(3, 3) :: m0
    real(r8), dimension(3) :: v0
    real(r8) :: Tdebye, nndist
    integer, dimension(:), allocatable :: di
    integer :: i, j, ti

    i = bs%n_point

    ! Kill early if it is unstable

    do i = 1, dr%n_irr_qpoint
    do j = 1, dr%n_mode
        if (dr%iq(i)%omega(j) .lt. -lo_freqtol) then
            if (mw%talk) write (*, *) 'Unstable, stopping now'
            call mw%destroy()
            stop
        end if
    end do
    end do

    ! Initialize the output file
    if (mw%talk) then
        call h5open_f(lo_status)
        call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, lo_status)
    end if

    ! Start by printing basic things
    if (mw%talk) then
        call lo_h5_store_data(temperatures, file_id, 'temperatures', lo_status, enhet='K')
        call lo_h5_store_data(uc%latticevectors*lo_bohr_to_A, file_id, 'lattice_vectors', lo_status, 'A')
        call lo_h5_store_data(uc%atomic_number, file_id, 'atomic_numbers', lo_status)
        call lo_h5_store_data(uc%r, file_id, 'atom_positions', lo_status, 'fractional')
    end if

    do i = 1, size(temperatures)
        da0(i) = dr%phonon_free_energy(temperatures(i))*lo_Hartree_to_eV
        da1(i) = dr%phonon_entropy(temperatures(i))*lo_Hartree_to_eV
        da2(i) = dr%phonon_cv(temperatures(i))*lo_Hartree_to_eV/lo_kb_eV
        da3(i) = da0(i) + da1(i)*temperatures(i)
    end do
    Tdebye = da0(1)*8.0_r8/9.0_r8/lo_kb_eV

    ! Store some of the thermodynamic stuff
    if ( mw%talk ) then
! h5store gets confused with pure scalars into attributes, and tries to subscript them.
! should fix the interface or add the attribute writing for scalars more
! cleanly. This works: passing a 1 element array
        call lo_h5_store_attribute([Tdebye],file_id,'debye_temperature',lo_status)
        call lo_h5_store_data(da0,file_id,'phonon_free_energy',lo_status,'eV/atom')
        call lo_h5_store_data(da1,file_id,'phonon_entropy',lo_status,'eV/atom/K')
        call lo_h5_store_data(da2,file_id,'phonon_heat_capacity',lo_status,'kB')
        call lo_h5_store_data(da3,file_id,'phonon_internal_energy',lo_status,'eV/atom')
    endif

    ! Get the commensurate modes?
    if (uc%na .ge. manyatoms) then
        call uc%classify('supercell', uc)
        call fc%commensurate_modes(uc, uc, fc, mw)
    end if

    ! Get the mean square displacement matrix and data from that
    nndist = uc%nearest_neighbour_distance()
    da0 = lo_huge
    da1 = -lo_huge
    da2 = 0.0_r8
    lo_allocate(thermal_sigma(3, 3, uc%na, size(temperatures)))
    thermal_sigma = 0.0_r8
    do ti = 1, size(temperatures)
        if (uc%na .ge. manyatoms) then
            call lo_thermal_displacement_matrix_commensurate(uc, uc, fc%eigenvectors, fc%omega, &
                                                             temperatures(ti), sigma=thermal_sigma(:, :, :, ti))
        else
            call dr%thermal_displacement_matrix(qp, uc, temperatures(ti), thermal_sigma(:, :, :, ti), mw, mem)
        end if
        ! Turn this into some scalar quantities?
        do i = 1, uc%na
            call lo_symmetric_eigensystem_3x3matrix(thermal_sigma(:, :, i, ti), v0, m0)
            da0(ti) = min(da0(ti), sqrt(minval(v0)))
            da1(ti) = max(da0(ti), sqrt(maxval(v0)))
            da2(ti) = da2(ti) + sum(sqrt(v0))/3.0_r8/uc%na
        end do
    end do
    ! Store this
    if (mw%talk) then
        call lo_h5_store_data(thermal_sigma*lo_bohr_to_A**2, file_id, 'thermal_displacement_matrices', lo_status, enhet='A^2')
        call lo_h5_store_data(da0*lo_bohr_to_A, file_id, 'min_displacement', lo_status, enhet='A')
        call lo_h5_store_data(da1*lo_bohr_to_A, file_id, 'max_displacement', lo_status, enhet='A')
        call lo_h5_store_data(da2*lo_bohr_to_A, file_id, 'mean_displacement', lo_status, enhet='A')
        call lo_h5_store_data(da0/nndist, file_id, 'min_displacement_scaled', lo_status, enhet='dimensionless')
        call lo_h5_store_data(da1/nndist, file_id, 'max_displacement_scaled', lo_status, enhet='dimensionless')
        call lo_h5_store_data(da2/nndist, file_id, 'mean_displacement_scaled', lo_status, enhet='dimensionless')
    end if

    ! The moments of the dos
    if (mw%talk) then
        lo_allocate(dy(pd%n_dos_point))
        dy = pd%dos/uc%na
        dy = dy*pd%omega
        da0(1) = lo_trapezoid_integration(pd%omega, dy)*lo_Hartree_to_eV
        dy = dy*pd%omega
        da0(2) = lo_trapezoid_integration(pd%omega, dy)*lo_Hartree_to_eV**2
        dy = dy*pd%omega
        da0(3) = lo_trapezoid_integration(pd%omega, dy)*lo_Hartree_to_eV**3
        dy = dy*pd%omega
        da0(4) = lo_trapezoid_integration(pd%omega, dy)*lo_Hartree_to_eV**4
        call lo_h5_store_data(da0, file_id, 'moments_of_phonon_dos', lo_status, enhet='eV')
    end if

    ! And create the magical fit, perhaps?
    if (mw%talk) then
        call fit_to_forceconstant(fc, uc, di, dy)
        call lo_h5_store_data(di, file_id, 'pairs_for_fc_model', lo_status, enhet='dimensionless')
        call lo_h5_store_data(dy, file_id, 'prefactor_for_fc_model', lo_status, enhet='dimensionless')
        call fake_forceconstant(fc, uc, di, dy)
        call fc%writetofile(uc, 'outfile.fakeforceconstant')
    end if

    ! Elastic constants?
    call fc%get_elastic_constants(uc)
    if (mw%talk) then
        call lo_h5_store_data(fc%elastic_constants_voigt*lo_pressure_HartreeBohr_to_GPa, &
                              file_id, 'elastic_constants', lo_status, enhet='GPa')
    end if

    ! Weird angular momentum thing?
    call dr%phonon_angular_momentum_matrix(qp, uc, 300.0_r8, m0, mw)
    if (mw%talk) then
        call lo_h5_store_data(m0, file_id, 'angular_momentum', lo_status, enhet='dunno')
    end if

    ! And close the output file
    if (mw%talk) then
        call h5fclose_f(file_id, lo_status)
        call h5close_f(lo_status)
    end if
end subroutine

subroutine fake_forceconstant(fc, uc, pairs, prefactors)
    !> unitcell
    type(lo_crystalstructure), intent(in) :: uc
    !> forceconstant
    type(lo_forceconstant_secondorder), intent(inout) :: fc
    !> pairs
    integer, dimension(:), intent(in) :: pairs
    !> prefactors
    real(r8), dimension(:), intent(in) :: prefactors

    integer :: a1, a2, i, j, k

    do a1 = 1, fc%na
    do i = 1, fc%atom(a1)%n
        if (lo_sqnorm(fc%atom(a1)%pair(i)%r) .lt. lo_sqtol) cycle
        a2 = fc%atom(a1)%pair(i)%i2
        j = topairid(uc%atomic_number(a1), uc%atomic_number(a2))
        do k = 1, size(pairs)
            if (pairs(k) .eq. j) then
                fc%atom(a1)%pair(i)%m = m_from_r(fc%atom(a1)%pair(i)%r, prefactors(k), n=4)
            end if
        end do
    end do
    end do
    call fc%setsumtozero()
end subroutine

!> Make my weird fit to forceconstant
subroutine fit_to_forceconstant(fc, uc, uniquepairs, prefactors)
    !> input forceconstant
    type(lo_forceconstant_secondorder), intent(inout) :: fc
    !> crystal structure
    type(lo_crystalstructure), intent(in) :: uc
    !> what are the unique pairs I am working with?
    integer, dimension(:), allocatable, intent(out) :: uniquepairs
    !> prefactor for the magic function
    real(r8), dimension(:), allocatable, intent(out) :: prefactors

    integer :: nun

    ! Get some basic stuff
    init: block
        integer, dimension(:), allocatable :: di
        integer :: i, j, l

        ! Grab the unique number of pairs of atomic numbers
        lo_allocate(di(uc%na**2))
        di = 0
        l = 0
        do i = 1, uc%na
        do j = 1, uc%na
            l = l + 1
            di(l) = topairid(uc%atomic_number(i), uc%atomic_number(j))
        end do
        end do
        call lo_return_unique(di, uniquepairs)
        nun = size(uniquepairs)

        ! Make space for the fitted forceconstants
        lo_allocate(prefactors(nun))
        prefactors = 0.0_r8

    end block init

    fit: block
        real(r8), dimension(:, :, :), allocatable :: forceconstant
        real(r8), dimension(:, :), allocatable :: vector
        integer :: un, ctr
        integer :: a1, a2, i, j

        do un = 1, nun
            ! Count all the pairs from os this kind
            ctr = 0
            do a1 = 1, fc%na
            do i = 1, fc%atom(a1)%n
                if (lo_sqnorm(fc%atom(a1)%pair(i)%r) .lt. lo_sqtol) cycle
                a2 = fc%atom(a1)%pair(i)%i2
                j = topairid(uc%atomic_number(a1), uc%atomic_number(a2))
                if (j .ne. uniquepairs(un)) cycle
                ctr = ctr + 1
            end do
            end do
            if (ctr .eq. 0) then
                ! Throw error maybe? Not sure
            end if
            ! Collect information from all these pairs
            lo_allocate(forceconstant(3, 3, ctr))
            lo_allocate(vector(3, ctr))

            forceconstant = 0.0_r8
            vector = 0.0_r8
            ctr = 0
            do a1 = 1, fc%na
            do i = 1, fc%atom(a1)%n
                if (lo_sqnorm(fc%atom(a1)%pair(i)%r) .lt. lo_sqtol) cycle
                a2 = fc%atom(a1)%pair(i)%i2
                j = topairid(uc%atomic_number(a1), uc%atomic_number(a2))
                if (j .ne. uniquepairs(un)) cycle
                ctr = ctr + 1
                vector(:, ctr) = fc%atom(a1)%pair(i)%r
                forceconstant(:, :, ctr) = fc%atom(a1)%pair(i)%m
            end do
            end do
            ! Now fit my magical function to this
            call minimize_forceconstant_function(vector, forceconstant, prefactors(un))

            lo_deallocate(forceconstant)
            lo_deallocate(vector)
        end do

    end block fit

contains

    subroutine minimize_forceconstant_function(vector, forceconstant, prefactor)
        real(r8), dimension(:, :), intent(in) :: vector
        real(r8), dimension(:, :, :), intent(in) :: forceconstant
        real(r8), intent(out) :: prefactor

        integer, parameter :: n = 4 ! best exponent
        real(r8), dimension(:, :, :), allocatable :: dfc
        integer :: i, np

        np = size(vector, 2)
        lo_allocate(dfc(3, 3, np))
        do i = 1, np
            dfc(:, :, i) = m_from_r(vector(:, i), 1.0_r8, n)
        end do
        prefactor = norm2(forceconstant)/norm2(dfc)
    end subroutine
end subroutine

function m_from_r(v, x, n) result(m)
    real(r8), intent(in), dimension(3) :: v
    real(r8), intent(in) :: x
    integer, intent(in) :: n
    real(r8), dimension(3, 3) :: m

    real(r8) :: r, rx, ry, rz, fpp

    r = norm2(v)
    if (r .lt. lo_tol) then
        m = 0.0_r8
    else
        rx = v(1)
        ry = v(2)
        rz = v(3)
        fpp = x/(r**n)
        m(1, 1) = -rx**2*r*fpp
        m(2, 2) = -ry**2*r*fpp
        m(3, 3) = -rz**2*r*fpp
        m(1, 2) = rx*ry*(-r*fpp)
        m(1, 3) = rx*rz*(-r*fpp)
        m(2, 3) = ry*rz*(-r*fpp)
        m(2, 1) = m(1, 2)
        m(3, 1) = m(1, 3)
        m(3, 2) = m(2, 3)
        m = m/(r**2)
        m = lo_chop(m, 1E-12_r8)
    end if
end function

! Convert two atomic numbers to one index.
function topairid(i, j) result(k)
    integer, intent(in) :: i, j
    integer :: k

    integer, parameter :: shift = 1000
    if (i .gt. j) then
        k = shift*j + i
    else
        k = shift*i + j
    end if
end function

! !> Weighted average of group velocity?
! function average_group_velocity(dr,temperature) result(avgvel)
!     !> dispersions
!     type(lo_phonon_dispersions), intent(in) :: dr
!     !> temperature
!     real(r8), intent(in) :: temperature
!     !> averaged group velocity
!     real(r8) :: avgvel
!
!
! end function

! !> Return a q-mesh density that takes the amount of atoms into account
! function lo_kmesh_density(uc,nq) result(density)
!     !> crystalstructure
!     type(lo_crystalstructure), intent(in) :: uc
!     !> linear density
!     integer, intent(in) :: nq
!     !> resulting q-density
!     integer, dimension(3) :: density
!
!     real(r8) :: f0,a,b,c,avg
!
!     ! So, what do I want. I want that product(density)*natoms = constant
!     ! So this is approximately the total number of points.
!     f0=real(nq**3,r8)/uc%na
!     ! And this would be the sensible number per direction
!     f0=f0**(1.0_r8/3.0_r8)
!
!     ! Then, It should be proportional to the length of the reciprocal lattice vectors
!     a=norm2(uc%reciprocal_latticevectors(:,1))
!     b=norm2(uc%reciprocal_latticevectors(:,2))
!     c=norm2(uc%reciprocal_latticevectors(:,3))
!     avg=(a+b+c)/3.0_r8
!     ! This should be somewhat ok, I guess
!     density(1)=ceiling(f0*a/avg)
!     density(2)=ceiling(f0*b/avg)
!     density(3)=ceiling(f0*c/avg)
! end function

! ! Calculate and dump the velocity dos
! subroutine dump_velocity_dos(dr,qp,temperature)
!     ! the phonon dispersions
!     type(lo_phonon_dispersions), intent(in) :: dr
!     !the qpoint grid
!     class(lo_qpoint_mesh), intent(in) :: qp
!     ! the temperature
!     real(r8), intent(in) :: temperature
!
!     integer :: i,j,k,l
!     integer :: n,u
!     real(r8) :: f0,f1,sigma
!     real(r8), dimension(:,:), allocatable :: dum,dos
!     real(r8), dimension(:), allocatable :: xax,dumdum
!
!     ! Sort all the velocity norms
!     lo_allocate(dum(dr%nq_irr,dr%nb))
!     do i=1,dr%nq_irr
!     do j=1,dr%nb
!         dum(i,j)=norm2(dr%iq(i)%vel(:,j))
!     enddo
!     enddo
!     !
!     u=open_file('out','omega_vs_vel')
!     do i=1,dr%nq_irr
!     do j=1,dr%nb
!         write(u,*) dr%iq(i)%omega(j),norm2(dr%iq(i)%vel(:,j))
!     enddo
!     enddo
!     close(u)
!     return
!     ! Figure out a sigma
!     lo_allocate(dumdum(dr%nq_irr*dr%nb))
!     l=0
!     do i=1,dr%nq_irr
!     do j=1,dr%nb
!         l=l+1
!         dumdum(l)=dum(i,j)
!     enddo
!     enddo
!     !
!     call qsort(dumdum)
!     f1=0.0_r8
!     do i=1,dr%nq_irr*dr%nb-1
!         f0=dumdum(i+1)-dumdum(i)
!         f1=max(f1,f0)
!     enddo
!     sigma=f1
!     ! Get an x-axis
!     n=300
!     lo_allocate(xax(n))
!     call lo_linspace(0.0_r8,maxval(dum)+6*sigma,xax)
!     ! Get the actual velocity DOS
!     lo_allocate(dos(dr%nb,n))
!     dos=0.0_r8
!     do i=1,dr%nq_irr
!     do j=1,dr%nb
!         f0=lo_planck(temperature,dr%iq(i)%omega(j))
!         do k=1,n
!             dos(j,k)=dos(j,k)+lo_gauss(dum(i,j),xax(k),sigma)*qp%ip(i)%weight
!         enddo
!     enddo
!     enddo
!     !
!     lo_deallocate(dum)
!     lo_allocate(dum(2,n))
!     dum=0.0_r8
!     do i=1,n
!         dum(1,i)=xax(i)
!         do j=1,dr%nb
!             dum(2,i)=dum(2,i)+dos(j,i)
!         enddo
!     enddo
!     call lo_dump_gnuplot_2d_real(dum,'outfile.velocity_dos')
!     !
! end subroutine

! Rather general routines, should probably be moved somewhere
subroutine lo_thermal_displacement_matrix_commensurate(ss, uc, eigenvectors, omega, temperature, thres, sigma, maxdisp)
    !> supercell
    type(lo_crystalstructure), intent(in) :: ss
    !> unitcell
    type(lo_crystalstructure), intent(in) :: uc
    !> eigenvectors
    real(r8), dimension(:, :), intent(in) :: eigenvectors
    !> frequencies
    real(r8), dimension(:), intent(in) :: omega
    !> temperature
    real(r8), intent(in) :: temperature
    !> threshold to replace imaginary frequencies with
    real(r8), intent(in), optional :: thres
    !> displacement matrix
    real(r8), dimension(:, :, :), intent(out), optional :: sigma
    !> max spherically averaged mean displacement
    real(r8), intent(out), optional :: maxdisp

    type(lo_distancetable) :: dt
    real(r8), dimension(3, 3) :: m0
    real(r8), dimension(3) :: v0
    real(r8) :: imfreqval, f0, om
    integer :: nuc, uca, i, a1

    !@todo add sanity test
    nuc = maxval(ss%info%index_in_unitcell)

    if (present(thres)) then
        imfreqval = thres
    else
        imfreqval = 0.0_r8
    end if

    if (present(sigma)) sigma = 0.0_r8
    if (present(maxdisp)) then
        call dt%generate(uc%r, uc%latticevectors, ss%nearest_neighbour_distance()*3, verbosity=0)
        maxdisp = 0.0_r8
    end if
    do uca = 1, nuc
        m0 = 0.0_r8
        a1 = -1
        do i = 1, ss%na
            if (ss%info%index_in_unitcell(i) .eq. uca) then
                a1 = i
                exit
            end if
        end do
        ! Now sum over all modes
        do i = 1, ss%na*3
            om = omega(i)
            if (om .lt. -lo_tiny) om = imfreqval
            if (abs(om) .lt. lo_tiny) cycle
            f0 = 0.5_r8*(1.0_r8 + 2.0_r8*lo_planck(temperature, om))/om
            v0 = eigenvectors((a1 - 1)*3 + 1:a1*3, i)
            m0 = m0 + lo_outerproduct(v0, v0)*f0
        end do
        ! And the mass
        m0 = m0/ss%mass(a1)
        ! Maybe store
        if (present(sigma)) sigma(:, :, uca) = m0
        !if ( present(maxdisp) ) then
        !    ! Use the largest mean square displacement and calculate
        !    ! a spherically averaged mean displacement from this.
        !    !call lo_symmetric_eigensystem_3x3matrix(m0,v0,m1)
        !    !f0=sqrt(maxval(v0))*gammafactor
        !    !maxdisp=max(maxdisp,f0)

        !    ! Not spherically averaged, look in the relevant directions instead.
        !    m1=lo_invert3x3matrix(m0)
        !    do i=2,dt%particle(uca)%n
        !        v0=dt%particle(uca)%v(:,i)/dt%particle(uca)%d(i)
        !        f0=1.0_r8/sqrt( dot_product(v0,matmul(m1,v0)) )
        !        f0=f0/dt%particle(uca)%d(i)
        !        maxdisp=max(maxdisp,f0)
        !    enddo
        !endif
    end do
end subroutine

! calculate U0
subroutine lo_calculate_U0(sim, fc, uc, ss, U0, mw, verbosity)
    !> simulation data
    type(lo_mdsim), intent(inout) :: sim
    !> forceconstant
    type(lo_forceconstant_secondorder), intent(inout) :: fc
    !> unitcell
    type(lo_crystalstructure), intent(in) :: uc
    !> supercell
    type(lo_crystalstructure), intent(inout) :: ss
    !> potential energy-ish
    real(r8), intent(out) :: U0
    !> mpi helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> talk a lot?
    integer, intent(in) :: verbosity

    type(lo_forceconstant_secondorder) :: fcss
    real(r8), dimension(:, :, :, :), allocatable :: fcp
    real(r8), dimension(3) :: v0
    real(r8) :: e2
    integer :: t, i1, a1, a2

    U0 = 0.0_r8

    ! First remap things
    if (verbosity .gt. 0) then
        write (*, *) ''
        write (*, *) 'PREPARING POTENTIAL ENERGY DIFFERENCES'
    end if

    call ss%classify('supercell', uc)
    call fc%remap(uc, ss, fcss)
    if (verbosity .gt. 0) then
        write (*, *) '... remapped second order'
    end if

    if (fc%polar) then
        allocate (fcp(3, 3, ss%na, ss%na))
        fcp = 0.0_r8
        call fc%supercell_longrange_dynamical_matrix_at_gamma(ss, fcp, 1E-15_r8)
    end if
    if (verbosity .gt. 0) then
        write (*, *) '... built polar forceconstant'
    end if

    ! Excellent, now calculate energies
    U0 = 0.0_r8
    do t = 1, sim%nt

        if (mod(t, mw%n) .ne. mw%r) cycle

        e2 = 0.0_r8
        do a1 = 1, fcss%na
            v0 = 0.0_r8
            do i1 = 1, fcss%atom(a1)%n
                a2 = fcss%atom(a1)%pair(i1)%i2
                v0 = v0 - matmul(fcss%atom(a1)%pair(i1)%m, sim%u(:, a2, t))
            end do
            e2 = e2 - dot_product(v0, sim%u(:, a1, t))*0.5_r8
        end do

        if (fc%polar) then
            do a1 = 1, sim%na
                v0 = 0.0_r8
                do a2 = 1, sim%na
                    v0 = v0 - matmul(fcp(:, :, a1, a2), sim%u(:, a2, t))
                end do
                e2 = e2 - dot_product(sim%u(:, a1, t), v0)*0.5_r8
            end do
        end if
        U0 = U0 + sim%stat%potential_energy(t) !-e2
    end do
    call mw%allreduce('sum', U0)

    U0 = U0/sim%na/sim%nt

end subroutine

end module
