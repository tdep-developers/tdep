module lo_accumulate_thermal_transport
!! Container for things related to thermal transport. Not how it is actually calculated,
!! that is generated elsewhere, but how it is stored and written to file.
use konstanter, only: r8, lo_iou, lo_hugeint, lo_huge, lo_sqtol, lo_freqtol, lo_pi, lo_exitcode_param, &
                      lo_groupvel_Hartreebohr_to_ms, lo_kb_Hartree, lo_kappa_au_to_SI, lo_frequency_Hartree_to_THz, &
                      lo_frequency_Hartree_to_icm, lo_frequency_Hartree_to_meV, lo_time_au_to_s, lo_bohr_to_A
use gottochblandat, only: tochar, walltime, lo_chop, &
                          lo_real_nullspace_coefficient_matrix, lo_transpositionmatrix, lo_real_pseudoinverse, &
                          lo_flattentensor, lo_trapezoid_integration, lo_linspace, lo_outerproduct, lo_planck, &
                          lo_harmonic_oscillator_cv, lo_unflatten_2tensor, lo_linear_interpolation
use mpi_wrappers, only: lo_mpi_helper, lo_stop_gracefully
use lo_memtracker, only: lo_mem_helper
use hdf5_wrappers, only: lo_hdf5_helper
use type_crystalstructure, only: lo_crystalstructure
use type_forceconstant_secondorder, only: lo_forceconstant_secondorder
use type_qpointmesh, only: lo_qpoint_mesh
use type_phonon_dispersions, only: lo_phonon_dispersions
use lo_distributed_phonon_dispersion_relations, only: lo_distributed_phonon_dispersions
use lo_distributed_phonon_dispersion_relations
use type_symmetryoperation, only: lo_expandoperation_pair, lo_eigenvector_transformation_matrix
use type_blas_lapack_wrappers, only: lo_gemm
use lo_spectralfunction_helpers, only: lo_evaluate_spectral_function, lo_integrate_spectral_function, lo_integrate_two_spectral_functions, lo_tapering_function
use quadratures_stencils, only: lo_centraldifference
implicit none

private
public :: lo_thermal_conductivity

!> container for all things thermal conductivity
type lo_thermal_conductivity
    !> number of modes
    integer, private :: n_mode = -lo_hugeint
    !> number of energy points in spectral plots
    integer, private :: n_energy = -lo_hugeint
    !> number of q-points on this rank?
    integer, private :: n_local_qpoint = -lo_hugeint
    !> temperature
    real(r8), private :: temperature = -lo_huge
    !> which global q-points do the local q-points correspond to?
    integer, dimension(:), allocatable :: qind
    !> energy axis in spectral plots
    real(r8), dimension(:), allocatable, private :: omega
    !> spectral kappa (energy,xyz,xyz,mode,mode)
    real(r8), dimension(:, :, :, :, :), allocatable :: spectral_kappa
    !> 'lifetime' per irreducible q-qpoint (mode,mode,qpoint)
    real(r8), dimension(:, :, :), allocatable :: lifetime
    !> 'mfp' per irreducible q-point (mode,mode,qpoint)
    real(r8), dimension(:, :, :), allocatable :: mean_free_path
    !> kappa per irreducible q-point (xyz,xyz,mode,mode,qpoint)
    real(r8), dimension(:, :, :, :, :), allocatable :: kappa

    !> number of probe frequencies
    integer :: n_omegaT=-lo_hugeint
    !> probe frequencies
    real(r8), dimension(:), allocatable :: omegaT
    !> AC thermal conductivity (xyz,xyz,energy)
    real(r8), dimension(:, :, :), allocatable :: AC_kappa
    real(r8), dimension(:, :, :), allocatable :: AC_kappa_od
    contains
        !> set up thermal conductivity container
        procedure :: initialize
        !> accumulate data from a single q-point
        !procedure :: accumulate
        !> report to screen
        !procedure :: report
        !> dump to file
        !procedure :: write_to_hdf5
        !> destroy
        procedure :: destroy=>destroy_tc
end type

contains

!> initialize thermal conductivity container
subroutine initialize(tc, dr, omega, omegaT, temperature, mw)
    !> thermal conductivity
    class(lo_thermal_conductivity), intent(out) :: tc
    !> phonon dispersions
    type(lo_distributed_phonon_dispersions), intent(in) :: dr
    !> frequency axis
    real(r8), dimension(:), intent(in) :: omega
    !> probe frequency axis
    real(r8), dimension(:), intent(in) :: omegaT
    !> temperature
    real(r8), intent(in) :: temperature
    !> MPI helper
    type(lo_mpi_helper), intent(inout) :: mw

    integer :: i, j

    tc%n_mode = dr%n_mode
    tc%n_energy = size(omega)
    tc%n_omegaT = size(omegaT)
    tc%temperature = temperature

    allocate(tc%omega(tc%n_energy))
    tc%omega=omega
    allocate(tc%omegaT(tc%n_omegaT))

    allocate(tc%AC_kappa(3,3,tc%n_omegaT))
    allocate(tc%AC_kappa_od(3,3,tc%n_omegaT))
    tc%AC_kappa=0.0_r8
    tc%AC_kappa_od=0.0_r8
end subroutine

! !> accumulate data from a single q-point on a grid.
! subroutine accumulate(tc, iq, local_iq, qp, dr, fc, p, spectral_smeared, sigmaIm, sigmaRe, scalefactor, xmid, xlo, xhi, mem)
!     !> thermal conductivity
!     class(lo_thermal_conductivity), intent(inout) :: tc
!     !> index to irreducible q-point we are working on
!     integer, intent(in) :: iq
!     !> local q-point index
!     integer, intent(in) :: local_iq
!     !> q-mesh
!     class(lo_qpoint_mesh), intent(in) :: qp
!     !> dispersion relations
!     type(lo_phonon_dispersions), intent(in) :: dr
!     !> second order force constant
!     type(lo_forceconstant_secondorder), intent(inout) :: fc
!     !> structure
!     type(lo_crystalstructure), intent(in) :: p
!     !> smeared normalized spectral function
!     real(r8), dimension(:, :), intent(in) :: spectral_smeared
!     !> imaginary part of the self-energy
!     real(r8), dimension(:, :), intent(in) :: sigmaIm
!     !> real part of the self-energy
!     real(r8), dimension(:, :), intent(in) :: sigmaRe
!     !> scaling factor per mode
!     real(r8), dimension(:), intent(in) :: scalefactor
!     !> peak location
!     real(r8), dimension(:), intent(in) :: xmid
!     !> left fwhm
!     real(r8), dimension(:), intent(in) :: xlo
!     !> right fwhm
!     real(r8), dimension(:), intent(in) :: xhi
!     !> memory tracker
!     type(lo_mem_helper), intent(inout) :: mem

!     real(r8), dimension(:, :, :), allocatable :: buf_vel
!     real(r8), dimension(:, :, :, :), allocatable :: buf_velsq
!     integer :: n_energy

!     ! Make note which q-point we are storing
!     tc%qind(local_iq)=iq

!     ! Space for group velocity buffers
!     call mem%allocate(buf_vel, [3, dr%n_mode, dr%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
!     call mem%allocate(buf_velsq, [3, 3, dr%n_mode, dr%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
!     buf_vel = 0.0_r8
!     buf_velsq = 0.0_r8
!     n_energy = size(spectral_smeared, 1)

!     ! Calculate the off-diagonal group velocity.
!     groupvel: block
!         complex(r8), dimension(:, :, :), allocatable :: buf_grad_dynmat
!         complex(r8), dimension(:, :), allocatable :: kronegv, buf_egv, buf_egw, buf_cm0, buf_cm1, buf_cm2
!         complex(r8), dimension(3) :: cv0
!         real(r8), dimension(3) :: v0, v1
!         integer :: a1, a2, ia, ib, ic, ix, iy, iz, k, iop, i, ii, j, jj

!         ! Some buffers
!         call mem%allocate(buf_grad_dynmat, [dr%n_mode, dr%n_mode, 3], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
!         call mem%allocate(buf_cm0, [dr%n_mode, dr%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
!         call mem%allocate(buf_cm1, [dr%n_mode**2, 3], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
!         call mem%allocate(buf_cm2, [dr%n_mode**2, 3], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
!         call mem%allocate(buf_egv, [dr%n_mode, dr%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
!         call mem%allocate(buf_egw, [dr%n_mode, dr%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
!         call mem%allocate(kronegv, [dr%n_mode**2, dr%n_mode**2], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)

!         buf_grad_dynmat = 0.0_r8
!         buf_cm0 = 0.0_r8
!         buf_cm1 = 0.0_r8
!         buf_cm2 = 0.0_r8
!         buf_egv = 0.0_r8
!         buf_egw = 0.0_r8
!         kronegv = 0.0_r8

!         ! Dynamical matrix and derivatives
!         call fc%dynamicalmatrix(p, qp%ip(iq), buf_cm0, mem, buf_grad_dynmat, qdirection=[1.0_r8, 0.0_r8, 0.0_r8])

!         ! Flatten gradient of dynamical matrix
!         do iz = 1, 3
!             do a1 = 1, p%na !ise%n_atom
!             do a2 = 1, p%na !ise%n_atom
!             do ix = 1, 3
!             do iy = 1, 3
!                 ib = (a1 - 1)*3 + ix
!                 ia = (a2 - 1)*3 + iy
!                 ic = flattenind(a1, a2, ix, iy, dr%n_mode)
!                 buf_cm1(ic, iz) = buf_grad_dynmat(ia, ib, iz)/(p%invsqrtmass(a1)*p%invsqrtmass(a2))
!             end do
!             end do
!             end do
!             end do
!         end do

!         ! Average over all operations
!         kronegv = 0.0_r8
!         do k = 1, qp%ip(iq)%n_invariant_operation
!             iop = qp%ip(iq)%invariant_operation(k)
!             ! Rotate eigenvectors
!             call lo_eigenvector_transformation_matrix(buf_cm0, p%rcart, qp%ip(iq)%r, p%sym%op(abs(iop)))
!             if (iop .lt. 0) then
!                 call lo_gemm(buf_cm0, dr%iq(iq)%egv, buf_egw)
!                 buf_egw = conjg(buf_egw)
!             else
!                 call lo_gemm(buf_cm0, dr%iq(iq)%egv, buf_egw)
!             end if

!             do i = 1, dr%n_mode
!                 if (dr%iq(iq)%omega(i) .gt. lo_freqtol) then
!                     do a1 = 1, p%na
!                     do ix = 1, 3
!                         ib = (a1 - 1)*3 + ix
!                         buf_egw(ib, i) = buf_egw(ib, i)*p%invsqrtmass(a1)/sqrt(dr%iq(iq)%omega(i)*2.0_r8)
!                     end do
!                     end do
!                 else
!                     buf_egw(:, i) = 0.0_r8
!                 end if
!             end do

!             do i = 1, dr%n_mode
!             do j = 1, dr%n_mode
!                 buf_egv = buf_egw*conjg(buf_egw(j, i))
!                 do ii = 1, dr%n_mode
!                 do jj = 1, dr%n_mode
!                     ia = (i - 1)*dr%n_mode + ii
!                     ib = (j - 1)*dr%n_mode + jj
!                     kronegv(ia, ib) = kronegv(ia, ib) + buf_egv(jj, ii)
!                 end do
!                 end do
!             end do
!             end do
!         end do
!         kronegv = kronegv/real(qp%ip(iq)%n_invariant_operation, r8)
!         ! this means sandwich with eigenvectors,  frequencies,
!         ! prefactors and masses are already in there.
!         call lo_gemm(kronegv, buf_cm1, buf_cm2)

!         ! Keep the group velocities?
!         do i = 1, dr%n_mode
!         do j = 1, dr%n_mode
!             ii = (i - 1)*dr%n_mode + j
!             !jj = (j - 1)*dr%n_mode + i
!             cv0 = buf_cm2(ii, :)
!             ! remove tiny numbers.
!             cv0 = lo_chop(cv0, 1E-10/(lo_groupvel_Hartreebohr_to_ms/1000))
!             ! I can take the real part since at the end we sum over
!             ! both modes and the imaginary components disappear.
!             buf_vel(:, i, j) = real(cv0, r8)
!         end do
!         end do

!         ! Be very careful with degeneracies. Will give very wrong contribution
!         ! to thermal transport if not taken care of properly, I'm pretty sure.
!         ! Feels like there could be significant double-counting otherwise.
!         ! This seems ok for now, but have to keep an eye out.

!         ! Rotate out the group velocities to get the symmetry averaged ones
!         buf_velsq = 0.0_r8
!         do i = 1, dr%n_mode
!         do j = 1, dr%n_mode
!             v0 = buf_vel(:, i, j)
!             do k = 1, qp%ip(iq)%n_full_point
!                 iop = qp%ip(iq)%operation_full_point(k)
!                 v1 = matmul(p%sym%op(abs(iop))%m, v0)
!                 buf_velsq(:, :, i, j) = buf_velsq(:, :, i, j) + lo_outerproduct(v1, v1)
!             end do
!         end do
!         end do
!         buf_velsq = buf_velsq/real(qp%ip(iq)%n_full_point, r8)

!         ! Sanity check that at least the diagonal ones are ok?
!         ! do i=1,dr%n_mode
!         !     write(*,*) i,dr%iq(iq)%vel(:,i)-buf_vel(:,i,i)
!         ! enddo

!         ! cleanup
!         call mem%deallocate(buf_grad_dynmat, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
!         call mem%deallocate(buf_cm0, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
!         call mem%deallocate(buf_cm1, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
!         call mem%deallocate(buf_cm2, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
!         call mem%deallocate(buf_egv, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
!         call mem%deallocate(buf_egw, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
!         call mem%deallocate(kronegv, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
!     end block groupvel

!     ! Integrate spectral functions and get thermal transport things.
!     integrate: block
!         real(r8), parameter :: integraltol = 1E-13_8
!         real(r8), dimension(:, :, :, :, :), allocatable :: buf_spectral
!         real(r8), dimension(:, :, :, :), allocatable :: buf_kappa
!         real(r8), dimension(:, :), allocatable :: buf_lifetime
!         real(r8), dimension(:), allocatable :: buf_thermal, buf0
!         real(r8) :: pref, f0, f1, f2
!         integer :: imode, jmode, ie

!         ! Smaller buffers for spectral kappa
!         call mem%allocate(buf0, n_energy, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
!         call mem%allocate(buf_thermal, n_energy, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
!         call mem%allocate(buf_spectral, [n_energy, 3, 3, tc%n_mode, tc%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
!         call mem%allocate(buf_kappa, [3, 3, tc%n_mode, tc%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
!         call mem%allocate(buf_lifetime, [tc%n_mode, tc%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
!         buf0 = 0.0_r8
!         buf_thermal = 0.0_r8
!         buf_spectral = 0.0_r8
!         buf_kappa = 0.0_r8
!         buf_lifetime = 0.0_r8

!         ! Thermal prefactor (for spectral plots, not integrals)
!         buf_thermal = 0.0_r8
!         do ie = 2, n_energy
!             f0 = lo_planck(tc%temperature, tc%omega(ie))
!             buf_thermal(ie) = f0*(f0 + 1.0_r8)
!         end do

!         do imode = 1, dr%n_mode
!         do jmode = 1, dr%n_mode

!             if (dr%iq(iq)%omega(imode) .lt. lo_freqtol) cycle
!             if (dr%iq(iq)%omega(jmode) .lt. lo_freqtol) cycle

!             if (imode .eq. jmode) then
!                 ! Diagonal part. Can use the easy way of doing things.
!                 call lo_integrate_spectral_function(tc%omega, dr%iq(iq)%omega(imode), sigmaIm(:, imode), sigmaRe(:, imode), xmid(imode), xlo(imode), xhi(imode), &
!                                                  scalefactor(imode), tc%temperature, integraltol, f0, f1, f2)
!             else
!                 ! Nondiagonal, trickier integral I suppose.
!                 call lo_integrate_two_spectral_functions(tc%omega, &
!                                                       dr%iq(iq)%omega(imode), dr%iq(iq)%omega(jmode), &
!                                                       sigmaIm(:, imode), sigmaIm(:, jmode), &
!                                                       sigmaRe(:, imode), sigmaRe(:, jmode), &
!                                                       xmid(imode), xmid(jmode), xlo(imode), xlo(jmode), xhi(imode), xhi(jmode), &
!                                                       scalefactor(imode), scalefactor(jmode), tc%temperature, integraltol*1E4_r8, f1, f2)
!             end if

!             ! Constant prefactor
!             pref = 1.0_r8 ! Omega^2 is added in the integrals, it's the correct way of doing things.
!             pref = pref*lo_pi/(lo_kb_hartree*tc%temperature**2*p%volume)

!             buf_lifetime(imode, jmode) = f1
!             buf_kappa(:, :, imode, jmode) = f2*pref*buf_velsq(:, :, imode, jmode)

!             ! Make sure the smeared ones normalize to the correct number.
!             buf0 = spectral_smeared(:, imode)*spectral_smeared(:, jmode)*buf_thermal
!             f0 = lo_trapezoid_integration(tc%omega, buf0)
!             buf0 = buf0*f2/f0
!             ! Spectral kappa
!             do ie = 2, tc%n_energy
!                 buf_spectral(ie, :, :, imode, jmode) = buf_spectral(ie, :, :, imode, jmode) + buf_velsq(:, :, imode, jmode)*buf0(ie)*qp%ip(iq)%integration_weight
!             end do
!         end do
!         end do

!         ! And store things on the corret rank.
!         tc%spectral_kappa = tc%spectral_kappa + buf_spectral
!         tc%lifetime(:, :, local_iq) = buf_lifetime
!         tc%kappa(:, :, :, :, local_iq) = buf_kappa
!         ! do imode = 1, dr%n_mode
!         !     tc%diagonal_imaginary_selfenergy(imode, local_iq) = lo_linear_interpolation(tc%omega, sigmaIm(:, imode), dr%iq(iq)%omega(imode), lo_freqtol*1E-5_r8)
!         ! end do

!         ! And some cleanup
!         call mem%deallocate(buf0, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
!         call mem%deallocate(buf_thermal, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
!         call mem%deallocate(buf_kappa, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
!         call mem%deallocate(buf_lifetime, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
!         call mem%deallocate(buf_spectral, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
!     end block integrate

!     ! Cleanup
!     call mem%deallocate(buf_vel, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
!     call mem%deallocate(buf_velsq, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)

!     contains
!     ! Consistent index flattening? Impossibru to get consistent.
!     function flattenind(a1, a2, ix, iy, nb) result(i)
!         integer, intent(in) :: a1, a2, ix, iy, nb
!         integer :: i

!         integer :: ia, ib

!         ia = (a1 - 1)*3 + ix
!         ib = (a2 - 1)*3 + iy
!         i = (ib - 1)*nb + ia
!     end function
! end subroutine

! !> report to screen how it's going?
! subroutine report(tc, qp, mw)
!     !> thermal conductivity
!     class(lo_thermal_conductivity), intent(in) :: tc
!     !> q-mesh
!     class(lo_qpoint_mesh), intent(in) :: qp
!     !> MPI helper
!     type(lo_mpi_helper), intent(inout) :: mw

!     real(r8), dimension(3, 3) :: m0, m1, m2
!     integer :: i, iq, s1, s2, ix
!     ! Simple integration

!     m0 = 0.0_r8
!     m1 = 0.0_r8
!     m2 = 0.0_r8
!     do i = 1, tc%n_local_qpoint
!         iq = tc%qind(i)
!         do s1 = 1, tc%n_mode
!         do s2 = 1, tc%n_mode
!             if (s1 .eq. s2) then
!                 m0 = m0 + tc%kappa(:, :, s1, s2, i)*qp%ip(iq)%integration_weight
!             else
!                 m1 = m1 + tc%kappa(:, :, s1, s2, i)*qp%ip(iq)%integration_weight
!             end if
!         end do
!         end do
!     end do
!     m2 = m0 + m1

!     call mw%allreduce('sum', m0)
!     call mw%allreduce('sum', m1)
!     call mw%allreduce('sum', m2)

!     m0 = lo_chop(m0*lo_kappa_au_to_SI, 1E-10_r8)
!     m1 = lo_chop(m1*lo_kappa_au_to_SI, 1E-10_r8)
!     m2 = lo_chop(m2*lo_kappa_au_to_SI, 1E-10_r8)

!     if (mw%talk) then
!         write (*, *) 'total:'
!         do ix = 1, 3
!             write (*, "(1X,3(1X,F17.7))") m2(:, ix)
!         end do
!         write (*, *) 'diagonal:'
!         do ix = 1, 3
!             write (*, "(1X,3(1X,F17.7))") m0(:, ix)
!         end do
!         write (*, *) 'off-diagonal:'
!         do ix = 1, 3
!             write (*, "(1X,3(1X,F17.7))") m1(:, ix)
!         end do
!     end if
! end subroutine

! !> Massage a little and write to file
! subroutine write_to_hdf5(tc, qp, dr, p, filename, enhet, mw, mem)
!     !> thermal conductivity
!     class(lo_thermal_conductivity), intent(inout) :: tc
!     !> q-mesh
!     class(lo_qpoint_mesh), intent(in) :: qp
!     !> phonon dispersions
!     type(lo_phonon_dispersions), intent(in) :: dr
!     !> structure
!     type(lo_crystalstructure), intent(in) :: p
!     !> filename
!     character(len=*), intent(in) :: filename
!     !> energy unit
!     character(len=3), intent(in) :: enhet
!     !> mpi helper
!     type(lo_mpi_helper), intent(inout) :: mw
!     !> memory tracker
!     type(lo_mem_helper), intent(inout) :: mem

!     type(lo_hdf5_helper) :: h5
!     real(r8), dimension(9, 9) :: symmetrizer
!     real(r8), dimension(:, :), allocatable :: coeff, icoeff
!     real(r8), dimension(:), allocatable :: x
!     real(r8) :: kappatol, unitfactor
!     integer :: nx
!     character(len=10) :: unitname

!     call mem%tick()

!     ! Prepare some initial stuffs
!     init: block
!         real(r8), dimension(:, :, :), allocatable :: ops
!         real(r8), dimension(3, 3) :: m0
!         real(r8), dimension(9) :: v9
!         integer :: i, iq, s1, s2, iop

!         if (mw%talk) then
!             write (*, *) ''
!             write (*, *) 'Writing thermal conductivity to file'
!         end if

!         ! First calculate the total Kappa.
!         m0 = 0.0_r8
!         do i = 1, tc%n_local_qpoint
!             iq = tc%qind(i)
!             do s1 = 1, tc%n_mode
!             do s2 = 1, tc%n_mode
!                 m0 = m0 + tc%kappa(:, :, s1, s2, i)*qp%ip(iq)%integration_weight
!             end do
!             end do
!         end do
!         call mw%allreduce('sum', m0)

!         ! if (mw%talk) then
!         !     write (*, *) ''
!         !     do i=1,3
!         !         write(*,*) 'm0',m0(i,:)*lo_kappa_au_to_SI
!         !     enddo
!         !     do iq=1,mw%n
!         !         do i=1,3
!         !             write(*,*) 'pr',iq,dr0(i,:,iq)*lo_kappa_au_to_SI
!         !         enddo
!         !     enddo
!         ! end if

!         ! Decide on a tolerance
!         kappatol = maxval(m0)*1E-14_r8
!         ! Use the tolerance to clean things neatly.
!         if (tc%n_local_qpoint .gt. 0) then
!             tc%kappa = lo_chop(tc%kappa, kappatol)
!             tc%spectral_kappa = lo_chop(tc%spectral_kappa, kappatol)
!         end if

!         ! Then get the matrix that can convert this to something
!         ! irreducible. Always easier to work with the irreducible
!         ! nonzero stuff and convert back later.
!         call mem%allocate(ops, [9, 9, p%sym%n + 1], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
!         ops = 0.0_r8
!         do iop = 1, p%sym%n
!             ops(:, :, iop) = lo_expandoperation_pair(p%sym%op(iop)%m)
!         end do
!         ! Add transposition as the last one, make the kappa tensor symmetric. Suppose it has to be?
!         ! Maybe not. I'm not sure. No it has to be, comes from the outer product stuff.
!         call lo_transpositionmatrix(ops(:, :, p%sym%n + 1))
!         call lo_real_nullspace_coefficient_matrix(invariant_operations=ops, coeff=coeff, nvar=nx)
!         allocate (icoeff(nx, 9))
!         ! And get the pseudoinverse.
!         call lo_real_pseudoinverse(coeff, icoeff)
!         ! And the irreducible representation of the total
!         allocate (x(nx))
!         v9 = lo_flattentensor(m0)
!         x = matmul(icoeff, v9)
!         call mem%deallocate(ops, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)

!         ! Get the matrix that least squares symmetrizes.
!         ! x = icoeff * v
!         ! w = coeff * x
!         ! w = coeff * icoeff * v
!         symmetrizer = lo_chop(matmul(coeff, icoeff), lo_sqtol)

!         ! Get the units right
!         select case (enhet)
!         case ('thz')
!             unitfactor = lo_frequency_Hartree_to_THz
!             unitname = 'Thz'
!         case ('icm')
!             unitfactor = lo_frequency_Hartree_to_icm
!             unitname = '1/cm'
!         case ('mev')
!             unitfactor = lo_frequency_Hartree_to_meV
!             unitname = 'meV'
!         case default
!             call lo_stop_gracefully(["Unknown unit, try 'thz', 'mev' or 'icm'"], lo_exitcode_param, __FILE__, __LINE__)
!         end select

!         if (mw%talk) then
!             write (*, *) '... sorted out symmetry'
!         end if

!         ! And open a file for output
!         if (mw%talk) then
!             ! Create a new file.
!             call h5%init(__FILE__, __LINE__)
!             call h5%open_file('write', trim(filename))

!             ! what q-mesh did I use? Store the metadata in the file, just to be on the safe side.
!             call h5%open_group('write', 'qmesh')
!             call qp%write_metadata_to_hdf5(input_id=h5%group_id)
!             call h5%close_group()
!             ! does not hurt to store the structure as well
!             call h5%open_group('write', 'structure')
!             call p%write_structure_to_hdf5(input_id=h5%group_id)
!             call h5%close_group()
!         end if
!     end block init

!     ! Relaxation time approximation, for reference
!     rta: block
!         real(r8), dimension(:, :), allocatable :: buf_lifetime, buf_mfp, buf_linewidth
!         real(r8), dimension(3, 3) :: m0, m1, m2
!         real(r8), dimension(9) :: v9
!         real(r8), dimension(3) :: v0, v1
!         real(r8) :: f0, cv
!         integer :: i, iq, imode, j, iop

!         ! Calculate lifetimes?
!         call mem%allocate(buf_linewidth, [dr%n_mode, qp%n_irr_point], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
!         call mem%allocate(buf_lifetime, [dr%n_mode, qp%n_irr_point], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
!         call mem%allocate(buf_mfp, [dr%n_mode, qp%n_irr_point], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)

!         buf_linewidth = 0.0_r8
!         buf_lifetime = 0.0_r8
!         buf_mfp = 0.0_r8
!         m1 = 0.0_r8
!         m2 = 0.0_r8
!         do i = 1, tc%n_local_qpoint
!             iq = tc%qind(i)
!             do imode = 1, dr%n_mode

!                 if (dr%iq(iq)%omega(imode) .lt. lo_freqtol) cycle

!                 ! Get the lifetime and mean free path
!                 !f0 = 2*tc%diagonal_imaginary_selfenergy(imode, i)
!                 buf_linewidth(imode, iq) = f0
!                 buf_lifetime(imode, iq) = 1.0_r8/f0
!                 buf_mfp(imode, iq) = norm2(dr%iq(iq)%vel(:, imode))*buf_lifetime(imode, iq)

!                 ! Get the outer product velocity thing
!                 m0 = 0.0_r8
!                 v0 = dr%iq(iq)%vel(:, imode)
!                 do j = 1, qp%ip(iq)%n_full_point
!                     iop = qp%ip(iq)%operation_full_point(j)
!                     v1 = matmul(p%sym%op(abs(iop))%m, v0)
!                     m0 = m0 + lo_outerproduct(v1, v1)
!                 end do
!                 m0 = m0/real(qp%ip(iq)%n_full_point, r8)

!                 ! Heat capacity?
!                 cv = lo_harmonic_oscillator_cv(tc%temperature, dr%iq(iq)%omega(imode))

!                 ! Accumulate thermal conductivity
!                 m1 = m1 + m0*cv*buf_lifetime(imode, iq)*qp%ip(iq)%integration_weight/p%volume
!                 m2 = m2 + m0*cv*tc%lifetime(imode, imode, i)*qp%ip(iq)%integration_weight/p%volume
!             end do
!         end do
!         call mw%allreduce('sum', buf_lifetime)
!         call mw%allreduce('sum', buf_linewidth)
!         call mw%allreduce('sum', buf_mfp)
!         call mw%allreduce('sum', m1)
!         call mw%allreduce('sum', m2)

!         ! make sure they are symmetric
!         v9 = lo_flattentensor(m1)
!         v9 = matmul(symmetrizer, v9)
!         m1 = lo_unflatten_2tensor(v9)

!         v9 = lo_flattentensor(m2)
!         v9 = matmul(symmetrizer, v9)
!         m2 = lo_unflatten_2tensor(v9)

!         if (mw%talk) then
!             write (*, *) 'RTA kappa'
!             do i = 1, 3
!                 write (*, *) m1(:, i)*lo_kappa_au_to_SI
!             end do
!             write (*, *) 'almost RTA kappa'
!             do i = 1, 3
!                 write (*, *) m2(:, i)*lo_kappa_au_to_SI
!             end do
!         end if

!         ! Store things to file. Can add more things later, if needed.
!         if (mw%talk) then
!             call h5%open_group('write', 'RTA')

!             ! store linewidths
!             buf_linewidth = buf_linewidth*unitfactor
!             call h5%store_data(buf_linewidth, h5%group_id, 'linewidth', enhet=unitname)

!             ! store frequencies as well.
!             do i = 1, qp%n_irr_point
!                 buf_linewidth(:, i) = dr%iq(i)%omega*unitfactor
!             end do
!             call h5%store_data(buf_linewidth, h5%group_id, 'frequencies', enhet=unitname)

!             ! Convert lifetimes to sensible unit.
!             ! probably SI, that's the most common one I think. Now it's in seconds.
!             buf_lifetime = buf_lifetime*lo_time_au_to_s
!             call h5%store_data(buf_lifetime, h5%group_id, 'lifetime', enhet='s')
!             ! Convert mean free path to SI as well
!             buf_mfp = buf_mfp*lo_bohr_to_A*1E-10_r8
!             call h5%store_data(buf_lifetime, h5%group_id, 'mean_free_path', enhet='s')

!             call h5%close_group()
!         end if

!         call mem%deallocate(buf_lifetime, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
!         call mem%deallocate(buf_linewidth, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
!         call mem%deallocate(buf_mfp, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
!     end block rta

!     ! Get the spectral representation
!     spectral: block
!         real(r8), dimension(:, :, :, :, :), allocatable :: bf2
!         real(r8), dimension(:, :, :), allocatable :: bf1
!         real(r8), dimension(:, :), allocatable :: sp_tot, bf0
!         real(r8), dimension(9) :: v9
!         real(r8) :: f0, f1
!         integer :: imode, jmode, ix, iy, ia, ctr, ie

!         call mem%allocate(sp_tot, [tc%n_energy, nx], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
!         call mem%allocate(bf0, [tc%n_energy, 9], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
!         bf0 = 0.0_r8
!         sp_tot = 0.0_r8

!         ! Symmetrize the total spectral kappa properly.
!         ctr = 0
!         do imode = 1, tc%n_mode
!         do jmode = 1, tc%n_mode

!             ctr = ctr + 1
!             if (mod(ctr, mw%n) .ne. mw%r) cycle

!             do iy = 1, 3
!             do ix = 1, 3
!                 ia = (ix - 1)*3 + iy
!                 bf0(:, ia) = tc%spectral_kappa(:, ix, iy, imode, jmode)
!             end do
!             end do
!             call lo_gemm(bf0, icoeff, sp_tot, transb='T', beta=1.0_r8)
!         end do
!         end do
!         call mw%allreduce('sum', sp_tot)

!         ! Make sure the total integrates to the right thing:
!         do ix = 1, nx
!             f0 = lo_trapezoid_integration(tc%omega, sp_tot(:, ix))
!             if (f0 .gt. kappatol) then
!                 f0 = x(ix)/f0
!             else
!                 f0 = 0.0_r8
!             end if
!             sp_tot(:, ix) = sp_tot(:, ix)*f0
!         end do

!         ! Now make sure the mode-projected are just as neat.
!         call mem%allocate(bf1, [nx, tc%n_mode, tc%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
!         call mem%allocate(bf2, [tc%n_energy, 3, 3, tc%n_mode, tc%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
!         bf1 = 0.0_r8
!         bf2 = 0.0_r8

!         do ie = 1, tc%n_energy

!             if (mod(ie, mw%n) .ne. mw%r) cycle

!             do jmode = 1, tc%n_mode
!             do imode = 1, tc%n_mode
!                 do ix = 1, 3
!                 do iy = 1, 3
!                     ia = (ix - 1)*3 + iy
!                     v9(ia) = (tc%spectral_kappa(ie, ix, iy, imode, jmode) + tc%spectral_kappa(ie, ix, iy, jmode, imode))*0.5_r8
!                     bf1(:, imode, jmode) = matmul(icoeff, v9)
!                 end do
!                 end do
!             end do
!             end do
!             ! Normalize
!             do ix = 1, nx
!                 f0 = sum(bf1(ix, :, :))
!                 if (f0 .gt. kappatol) then
!                     f1 = sp_tot(ie, ix)/f0
!                 else
!                     f1 = 0.0_r8
!                 end if
!                 bf1(ix, :, :) = bf1(ix, :, :)*f1
!             end do
!             ! Put back in the right place
!             do imode = 1, tc%n_mode
!             do jmode = imode, tc%n_mode
!                 v9 = matmul(coeff, bf1(:, imode, jmode))
!                 do ix = 1, 3
!                 do iy = 1, 3
!                     ia = (ix - 1)*3 + iy
!                     bf2(ie, ix, iy, imode, jmode) = v9(ia)
!                 end do
!                 end do
!             end do
!             end do
!         end do
!         call mw%allreduce('sum', bf2)
!         tc%spectral_kappa = bf2

!         ! Write some spectral things to file
!         if (mw%talk) then
!             bf2 = tc%spectral_kappa
!             bf2 = bf2*lo_kappa_au_to_SI/unitfactor
!             ! f0=0.0_r8
!             ! do imode=1,tc%n_mode
!             !     f0=f0+lo_trapezoid_integration(tc%omega,tc%spectral_kappa(:,1,1,imode,imode))*lo_kappa_au_to_SI
!             !     write(*,*) 'INTEGRATED SPECTRAL',imode,lo_trapezoid_integration(tc%omega,tc%spectral_kappa(:,1,1,imode,imode))*lo_kappa_au_to_SI,f0
!             ! enddo
!             call h5%store_data(tc%omega*unitfactor, h5%file_id, 'energy_axis', enhet=unitname)
!             call h5%store_data(bf2, h5%file_id, 'spectral_kappa')
!             write (*, *) '... stored spectral kappa'
!         end if

!         call mem%deallocate(bf0, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
!         call mem%deallocate(bf1, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
!         call mem%deallocate(bf2, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
!         call mem%deallocate(sp_tot, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
!     end block spectral

!     ! Close file
!     if (mw%talk) then
!         call h5%close_file()
!         call h5%destroy(__FILE__, __LINE__)
!     end if

!     call mem%tock(__FILE__, __LINE__, mw%comm)
! end subroutine

subroutine destroy_tc(tc)
    class(lo_thermal_conductivity), intent(inout) :: tc

    tc%n_mode = -lo_hugeint
    tc%n_energy = -lo_hugeint
    tc%n_local_qpoint = -lo_hugeint
    tc%temperature = -lo_huge

    if ( allocated(tc%qind                         ) ) deallocate(tc%qind                         )
    if ( allocated(tc%omega                        ) ) deallocate(tc%omega                        )
    if ( allocated(tc%spectral_kappa               ) ) deallocate(tc%spectral_kappa               )
    if ( allocated(tc%lifetime                     ) ) deallocate(tc%lifetime                     )
    if ( allocated(tc%mean_free_path               ) ) deallocate(tc%mean_free_path               )
    if ( allocated(tc%kappa                        ) ) deallocate(tc%kappa                        )
    !if ( allocated(tc%diagonal_imaginary_selfenergy) ) deallocate(tc%diagonal_imaginary_selfenergy)
end subroutine

end module
