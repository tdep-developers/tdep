module lo_selfenergy_interpolation
!! Allows us to evaluate the phonon self-energy at arbitrary q-points
use konstanter, only: r8, lo_iou, lo_hugeint, lo_huge, lo_sqtol, lo_exitcode_symmetry, lo_twopi, &
                      lo_bohr_to_A, lo_freqtol, lo_exitcode_param, lo_pi, lo_imag, lo_groupvel_ms_to_Hartreebohr
use gottochblandat, only: tochar, walltime, lo_progressbar_init, lo_progressbar, lo_chop, &
                          lo_linear_least_squares, lo_rsquare, lo_linspace, lo_linear_interpolation, lo_trapezoid_integration,&
                          lo_clean_fractional_coordinates, lo_cross
use geometryfunctions, only: lo_inscribed_sphere_in_box, lo_plane
use mpi_wrappers, only: lo_mpi_helper, lo_stop_gracefully
use lo_memtracker, only: lo_mem_helper
use hdf5_wrappers, only: lo_hdf5_helper
use type_crystalstructure, only: lo_crystalstructure
use type_forceconstant_secondorder, only: lo_forceconstant_secondorder
use type_qpointmesh, only: lo_qpoint_mesh, lo_fft_mesh, lo_wedge_mesh, lo_qpoint, lo_read_qmesh_from_file
use type_phonon_dispersions, only: lo_phonon_dispersions, lo_phonon_dispersions_qpoint
use lo_symmetry_of_interactions, only: lo_interaction_tensors
use type_symmetryoperation, only: lo_eigenvector_transformation_matrix
use type_blas_lapack_wrappers, only: lo_gemm, lo_dgels
use lo_phonon_bandstructure_on_path, only: lo_phonon_bandstructure


use lineshape_helper, only: lo_spectralfunction_helper,evaluate_spectral_function
use lo_tetrahedron_interpolation, only: lo_linear_tetrahedron_interpolation
implicit none

private
public :: lo_interpolated_selfenergy_grid

type lo_interpolated_selfenergy_grid
    !> point locator thingy
    type(lo_linear_tetrahedron_interpolation) :: box
    !> q-point mesh used to generate the self-energy
    class(lo_qpoint_mesh), allocatable :: qp
    !> number of energies
    integer :: n_energy=-lo_hugeint
    !> energy axis
    real(r8), dimension(:), allocatable :: omega
    !> real part of self-energy (xyz,xyz,energy,q)
    real(r8), dimension(:,:,:,:), allocatable :: sigma_Re
    !> imaginary part of self-energy (xyz,xyz,energy,q)
    real(r8), dimension(:,:,:,:), allocatable :: sigma_Im
    contains
        procedure :: read_from_hdf5=>read_interpolated_selfenergy_from_hdf5
        procedure :: evaluate=>evaluate_self_energy
        procedure :: destroy=>destroy_interpolated_selfenergy
        procedure :: spectral_function_along_path=>spectral_function_path_interp
end type

interface ! to evaluate
    module subroutine evaluate_self_energy(ise,p,qv,wp,sigma_Re,sigma_Im)
        !> self-energy interpolation
        class(lo_interpolated_selfenergy_grid), intent(in) :: ise
        type(lo_crystalstructure), intent(in) :: p

        !> q-point
        real(r8), dimension(3), intent(in) :: qv
        !> harmonic phonon properties at this q
        type(lo_phonon_dispersions_qpoint), intent(in) :: wp
        !> real part of self-energy
        real(r8), dimension(:,:), intent(out) :: sigma_Re
        !> imaginary part of self-energy
        real(r8), dimension(:,:), intent(out) :: sigma_Im
    end subroutine
end interface

interface ! to path
    module subroutine spectral_function_path_interp(ise, bs, uc, mw, mem)
        class(lo_interpolated_selfenergy_grid), intent(in) :: ise
        type(lo_phonon_bandstructure), intent(inout) :: bs
        type(lo_crystalstructure), intent(inout) :: uc
        type(lo_mpi_helper), intent(inout) :: mw
        type(lo_mem_helper), intent(inout) :: mem
    end subroutine
end interface

contains

subroutine read_interpolated_selfenergy_from_hdf5(ise,p,filename,mw,mem,verbosity)
    !> self-energy
    class(lo_interpolated_selfenergy_grid), intent(out) :: ise
    !> crystal structure
    type(lo_crystalstructure), intent(inout) :: p
    !> filename
    character(len=*), intent(in) :: filename
    !> mpi helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem
    !> talk a lot?
    integer, intent(in) :: verbosity

    ! First we grab the raw data from file.
    readfile: block
        real(r8), dimension(:,:,:), allocatable :: rbuf
        type(lo_hdf5_helper) :: h5
        integer :: iq

        if ( verbosity .gt. 0 ) then
            write(*,*) ''
            write(*,*) 'Reading interpolated self-energy from file'
        endif

        call h5%init(__FILE__,__LINE__)
        call h5%open_file('read',trim(filename))

        ! First we read the q-point mesh from file
        call h5%open_group('read','qmesh')
        call lo_read_qmesh_from_file(ise%qp,p,'null',mem,0,h5%group_id)
        call h5%close_group()

        if ( verbosity .gt. 0 ) then
            write(*,*) '... read q-mesh'
        endif

        ! Then we grab the energy-axis, seems reasonable
        call h5%read_data(ise%omega,h5%file_id,'omega')
        ise%n_energy = size(ise%omega)

        ! Create buffer space
        allocate(ise%sigma_Im(p%na*3,p%na*3,size(ise%omega),ise%qp%n_irr_point))
        allocate(ise%sigma_Re(p%na*3,p%na*3,size(ise%omega),ise%qp%n_irr_point))
        ise%sigma_Im=0.0_r8
        ise%sigma_Re=0.0_r8

        do iq=1,ise%qp%n_irr_point
            call h5%open_group('read','selfenergy_qpoint_'//tochar(iq))

            call h5%read_data(rbuf,h5%group_id,'sigma_Re')
            ise%sigma_Re(:,:,:,iq)=rbuf
            deallocate(rbuf)
            call h5%read_data(rbuf,h5%group_id,'sigma_Im')
            ise%sigma_Im(:,:,:,iq)=rbuf
            deallocate(rbuf)

            call h5%close_group()
        enddo

        if ( verbosity .gt. 0 ) then
            write(*,*) '... read self-energy'
        endif

        call h5%close_file()
        call h5%destroy()

        if ( verbosity .gt. 0 ) then
            write(*,*) 'Done reading self-energy from file'
        endif
    end block readfile

    ! Generate triangulation thingy
    call ise%box%generate(ise%qp,p)
end subroutine

subroutine destroy_interpolated_selfenergy(ise)
    class(lo_interpolated_selfenergy_grid), intent(inout) :: ise

    call ise%box%destroy()
    if ( allocated(ise%qp) ) then
        call ise%qp%destroy(ise%qp)
        deallocate(ise%qp)
    endif
    ise%n_energy=-lo_hugeint
    if ( allocated(ise%omega   ) ) deallocate(ise%omega   )
    if ( allocated(ise%sigma_Re) ) deallocate(ise%sigma_Re)
    if ( allocated(ise%sigma_Im) ) deallocate(ise%sigma_Im)
end subroutine

! ! Consistent index flattening? Impossibru to get consistent.
! function flattenind(a1, a2, ix, iy, nb) result(i)
!     integer, intent(in) :: a1, a2, ix, iy, nb
!     integer :: i

!     integer :: ia, ib

!     ia = (a1 - 1)*3 + ix
!     ib = (a2 - 1)*3 + iy
!     i = (ib - 1)*nb + ia
! end function

end module
