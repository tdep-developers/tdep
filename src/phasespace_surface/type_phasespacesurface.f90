#include "precompilerdefinitions"
module type_phasespacesurface
use konstanter, only: flyt, lo_pi, lo_huge, lo_freqtol, lo_status, lo_sqtol, lo_twopi
use gottochblandat, only: tochar, lo_cross, open_file, lo_signed_tetrahedron_volume, walltime, lo_sqnorm, &
                          lo_progressbar_init, lo_progressbar, qsort
use lo_memtracker, only: lo_mem_helper
use type_crystalstructure, only: lo_crystalstructure
use type_forceconstant_secondorder, only: lo_forceconstant_secondorder
use type_forceconstant_thirdorder, only: lo_forceconstant_thirdorder
use type_phonon_dispersions, only: lo_phonon_dispersions, lo_phonon_dispersions_qpoint
use type_qpointmesh, only: lo_qpoint_mesh, lo_wedge_mesh, lo_qpoint, lo_get_small_group_of_qpoint
use hdf5_wrappers, only: lo_h5_store_data, lo_h5_store_attribute, HID_T, H5F_ACC_TRUNC_F, h5close_f, h5open_f, &
                         h5fclose_f, h5fopen_f, h5fcreate_f, h5gclose_f, h5gopen_f, h5gcreate_f

implicit none
public :: lo_phasespacesurface

type lo_phasespacesurface_patch
    !> number of triangles
    integer :: ntri
    !> number of points
    integer :: npts
    !> the triangles
    integer, dimension(:, :), allocatable :: tri
    !> the points
    real(flyt), dimension(:, :), allocatable :: pts
    !> the gradient of e(k) over each triangle
    real(flyt), dimension(:, :), allocatable :: grad
    !> the matrix element (intensity) at each point
    real(flyt), dimension(:), allocatable :: psisquared
    !> the frequencies of the modes involved at each point (index: point, modelabel)
    real(flyt), dimension(:, :), allocatable :: omega
    !> the eigenvectors of the modes involved at each point (index: point, :, modelabel)
    real(flyt), dimension(:, :, :), allocatable :: egv
    !> the q-vectors of the modes involved at each point (index: point, :, modelabel)
    real(flyt), dimension(:, :, :), allocatable :: q
end type

type lo_phasespacesurface
    !> the energy of the isosurface
    real(flyt) :: energy
    !> how many bands are in this isosurface
    integer :: nb
    !> what are the global indices to these bands
    integer, dimension(:), allocatable :: bandind
    !> the actual surfaces
    type(lo_phasespacesurface_patch), dimension(:, :, :), allocatable :: plus
    type(lo_phasespacesurface_patch), dimension(:, :, :), allocatable :: minus
contains
    !> build the fermisurface
    procedure :: generate
    !> dump it to file
    procedure :: write_to_hdf5
    !> dump it to povray
    procedure :: write_to_povray
end type

!> helper type to get unique points
type lo_points_in_boxes
    integer :: n
    integer, dimension(:), allocatable :: ind
end type

contains

#include "type_phasespacesurface_cleanmesh.f90"
#include "type_phasespacesurface_energycalculations.f90"
#include "type_phasespacesurface_povray.f90"

!> generate three-phonon scattering surfaces
subroutine generate(ps, qp, uc, fc, fct, point, verbosity, modespec, calcintens, mem)
    !> fermisurface
    class(lo_phasespacesurface), intent(out) :: ps
    !> qpoint mesh
    type(lo_wedge_mesh), intent(in) :: qp
    !> unitcell
    type(lo_crystalstructure), intent(in) :: uc
    !> forceconstant
    type(lo_forceconstant_secondorder), intent(inout) :: fc
    !> thirdorder forceconstant
    type(lo_forceconstant_thirdorder), intent(in) :: fct
    !> the q-point in question
    type(lo_qpoint), intent(in) :: point
    !> how much to talk
    integer, intent(in) :: verbosity
    !> perhaps not all modes
    integer, dimension(:, :), intent(in), optional :: modespec
    !> calculate intensities (matrix elements)?
    logical, intent(in), optional :: calcintens
    type(lo_mem_helper), intent(inout) :: mem
    !
    integer :: i, j, k, l, u, ii, srf
    integer :: b1, b2, b3
    real(flyt), dimension(3, 3) :: tripts
    real(flyt), dimension(3) :: v0, v1, v2, v3
    real(flyt) :: f0
    real(flyt), dimension(:, :, :, :), allocatable :: energy_plus, energy_minus
    logical :: calcintensities
    !
    if (present(calcintens)) then
        calcintensities = calcintens
    else
        calcintensities = .false.
    end if
    !
    if (verbosity .gt. 0) write (*, *) ''
    if (verbosity .gt. 0) write (*, *) 'Calculating isosurfaces'
    ps%nb = uc%na*3
    lo_allocate(ps%plus(ps%nb, ps%nb, ps%nb))
    lo_allocate(ps%minus(ps%nb, ps%nb, ps%nb))
    lo_allocate(energy_plus(qp%n_full_point, ps%nb, ps%nb, ps%nb))
    lo_allocate(energy_minus(qp%n_full_point, ps%nb, ps%nb, ps%nb))
    do b1 = 1, ps%nb
    do b2 = 1, ps%nb
    do b3 = 1, ps%nb
        ps%plus(b1, b2, b3)%ntri = 0
        ps%plus(b1, b2, b3)%npts = 0
        ps%minus(b1, b2, b3)%ntri = 0
        ps%minus(b1, b2, b3)%npts = 0
    end do
    end do
    end do
    if (verbosity .gt. 0) write (*, *) '... made some initial space'

    ! Get the energy values
    call energy_values_on_vertices(point, qp, uc, fc, energy_plus, energy_minus, mem=mem)
    if (verbosity .gt. 0) write (*, *) '... got energies at nodes'

    if (present(modespec)) then
        ! Calculate only some surfaces
        do srf = 1, size(modespec, 2)
            write (*, *) 'building surface ', tochar(srf), ' out of ', tochar(size(modespec, 2))
            b1 = modespec(1, srf)
            b2 = modespec(2, srf)
            b3 = modespec(3, srf)
            ii = modespec(4, srf)
            ! get the isosurfaces
            if (ii .eq. 0) then
                call bz_isosurface(qp, energy_plus(:, b1, b2, b3), 0.0_flyt, ps%plus(b1, b2, b3)%pts, &
                                   ps%plus(b1, b2, b3)%tri, ps%plus(b1, b2, b3)%grad, &
                                   uc, fc, point, .true., b1, b2, b3, verbosity=1, &
                                   fct=fct, calcintens=calcintensities, psisq=ps%plus(b1, b2, b3)%psisquared, &
                                   tot_omega=ps%plus(b1, b2, b3)%omega, tot_egv=ps%plus(b1, b2, b3)%egv, &
                                   tot_q=ps%plus(b1, b2, b3)%q)
            else
                call bz_isosurface(qp, energy_minus(:, b1, b2, b3), 0.0_flyt, ps%minus(b1, b2, b3)%pts, &
                                   ps%minus(b1, b2, b3)%tri, ps%minus(b1, b2, b3)%grad, &
                                   uc, fc, point, .false., b1, b2, b3, verbosity=1, &
                                   fct=fct, calcintens=calcintensities, psisq=ps%plus(b1, b2, b3)%psisquared, &
                                   tot_omega=ps%minus(b1, b2, b3)%omega, tot_egv=ps%minus(b1, b2, b3)%egv, &
                                   tot_q=ps%minus(b1, b2, b3)%q)
            end if
            if (allocated(ps%plus(b1, b2, b3)%tri)) then
                ps%plus(b1, b2, b3)%ntri = size(ps%plus(b1, b2, b3)%tri, 2)
                ps%plus(b1, b2, b3)%npts = size(ps%plus(b1, b2, b3)%pts, 2)
            else
                ps%plus(b1, b2, b3)%ntri = 0
                ps%plus(b1, b2, b3)%npts = 0
            end if
            if (allocated(ps%minus(b1, b2, b3)%tri)) then
                ps%minus(b1, b2, b3)%ntri = size(ps%minus(b1, b2, b3)%tri, 2)
                ps%minus(b1, b2, b3)%npts = size(ps%minus(b1, b2, b3)%pts, 2)
            else
                ps%minus(b1, b2, b3)%ntri = 0
                ps%minus(b1, b2, b3)%npts = 0
            end if
        end do
    else
        ! If not specified, calculate all the surfaces
        l = 0
        call lo_progressbar_init()
        do b1 = 1, ps%nb
        do b2 = 1, ps%nb
        do b3 = 1, ps%nb
            ! get the isosurfaces
            call bz_isosurface(qp, energy_plus(:, b1, b2, b3), 0.0_flyt, ps%plus(b1, b2, b3)%pts, &
                               ps%plus(b1, b2, b3)%tri, ps%plus(b1, b2, b3)%grad, &
                               uc, fc, point, .true., b1, b2, b3, &
                               fct=fct, calcintens=calcintensities, psisq=ps%plus(b1, b2, b3)%psisquared, &
                               tot_omega=ps%plus(b1, b2, b3)%omega, tot_egv=ps%plus(b1, b2, b3)%egv, &
                               tot_q=ps%plus(b1, b2, b3)%q)
            call bz_isosurface(qp, energy_minus(:, b1, b2, b3), 0.0_flyt, ps%minus(b1, b2, b3)%pts, &
                               ps%minus(b1, b2, b3)%tri, ps%minus(b1, b2, b3)%grad, &
                               uc, fc, point, .false., b1, b2, b3, &
                               fct=fct, calcintens=calcintensities, psisq=ps%plus(b1, b2, b3)%psisquared, &
                               tot_omega=ps%minus(b1, b2, b3)%omega, tot_egv=ps%minus(b1, b2, b3)%egv, &
                               tot_q=ps%minus(b1, b2, b3)%q)
            if (allocated(ps%plus(b1, b2, b3)%tri)) then
                ps%plus(b1, b2, b3)%ntri = size(ps%plus(b1, b2, b3)%tri, 2)
                ps%plus(b1, b2, b3)%npts = size(ps%plus(b1, b2, b3)%pts, 2)
            else
                ps%plus(b1, b2, b3)%ntri = 0
                ps%plus(b1, b2, b3)%npts = 0
            end if
            if (allocated(ps%minus(b1, b2, b3)%tri)) then
                ps%minus(b1, b2, b3)%ntri = size(ps%minus(b1, b2, b3)%tri, 2)
                ps%minus(b1, b2, b3)%npts = size(ps%minus(b1, b2, b3)%pts, 2)
            else
                ps%minus(b1, b2, b3)%ntri = 0
                ps%minus(b1, b2, b3)%npts = 0
            end if
            l = l + 1
            call lo_progressbar(' ... isosurfaces', l, ps%nb**3)
        end do
        end do
        end do
    end if

    if (verbosity .gt. 0) write (*, *) '... got isosurfaces'

    ! Align the triangles
    do b1 = 1, ps%nb
    do b2 = 1, ps%nb
    do b3 = 1, ps%nb
        v0 = 0.0_flyt
        do i = 1, ps%plus(b1, b2, b3)%npts
            v0 = v0 + ps%plus(b1, b2, b3)%pts(:, i)/ps%plus(b1, b2, b3)%npts
        end do
        do i = 1, ps%plus(b1, b2, b3)%ntri
            tripts = ps%plus(b1, b2, b3)%pts(:, ps%plus(b1, b2, b3)%tri(:, i))
            v1 = tripts(:, 1) - tripts(:, 3)
            v2 = tripts(:, 2) - tripts(:, 3)
            v3 = lo_cross(v1, v2)
            !do j=1,3
            !    v1(j)=sum(tripts(j,:))/3.0_flyt
            !enddo
            if (dot_product(v3, v0) .lt. 0.0_flyt) then
                ps%plus(b1, b2, b3)%tri(:, i) = ps%plus(b1, b2, b3)%tri([3, 2, 1], i)
            end if
        end do
        v0 = 0.0_flyt
        do i = 1, ps%minus(b1, b2, b3)%npts
            v0 = v0 + ps%minus(b1, b2, b3)%pts(:, i)/ps%minus(b1, b2, b3)%npts
        end do
        do i = 1, ps%minus(b1, b2, b3)%ntri
            tripts = ps%minus(b1, b2, b3)%pts(:, ps%minus(b1, b2, b3)%tri(:, i))
            v1 = tripts(:, 1) - tripts(:, 3)
            v2 = tripts(:, 2) - tripts(:, 3)
            v3 = lo_cross(v1, v2)
            !do j=1,3
            !    v1(j)=sum(tripts(j,:))/3.0_flyt
            !enddo
            if (dot_product(v3, v0) .lt. 0.0_flyt) then
                ps%minus(b1, b2, b3)%tri(:, i) = ps%minus(b1, b2, b3)%tri([3, 2, 1], i)
            end if
        end do
    end do
    end do
    end do

    if (verbosity .gt. 0) write (*, *) '... got aligned triangles'
    !
end subroutine

!> dump the surfaces as triangles
subroutine write_to_hdf5(ps, uc, filename)
    !> fermi surface
    class(lo_phasespacesurface), intent(in) :: ps
    !> unitcell
    type(lo_crystalstructure), intent(in) :: uc
    !> filename
    character(len=*), intent(in) :: filename

    integer :: i, j, k, np, nm
    integer :: b1, b2, b3
    integer(HID_T) :: file_id, group_id
    character(len=1000) :: dum
    !
    write (*, *) ''
    write (*, *) 'Writing surfaces to hdf5 file'
    call h5open_f(lo_status); if (lo_status .ne. 0) stop
    call h5fcreate_f(trim(filename), H5F_ACC_TRUNC_F, file_id, lo_status); if (lo_status .ne. 0) stop
    ! store some metadata
    call lo_h5_store_attribute(ps%nb, file_id, 'number_of_bands')
    do b1 = 1, ps%nb
    do b2 = 1, ps%nb
    do b3 = 1, ps%nb
        np = ps%plus(b1, b2, b3)%ntri
        dum = 'number_of_plus_triangles_'//tochar(b1)//'_'//tochar(b2)//'_'//tochar(b3)
        call lo_h5_store_attribute(np, file_id, trim(dum))
        if (np .gt. 0) then
            dum = 'triangles_plus_'//tochar(b1)//'_'//tochar(b2)//'_'//tochar(b3)
            call lo_h5_store_data(ps%plus(b1, b2, b3)%tri, file_id, trim(dum))
            dum = 'points_plus_'//tochar(b1)//'_'//tochar(b2)//'_'//tochar(b3)
            call lo_h5_store_data(ps%plus(b1, b2, b3)%pts, file_id, trim(dum))
            if (any(ps%plus(b1, b2, b3)%psisquared .gt. 0.0_flyt)) then
                dum = 'psisq_plus_'//tochar(b1)//'_'//tochar(b2)//'_'//tochar(b3)
                call lo_h5_store_data(ps%plus(b1, b2, b3)%psisquared, file_id, trim(dum))
                dum = 'omega_plus_'//tochar(b1)//'_'//tochar(b2)//'_'//tochar(b3)
                call lo_h5_store_data(ps%plus(b1, b2, b3)%omega, file_id, trim(dum))
                dum = 'egv_plus_'//tochar(b1)//'_'//tochar(b2)//'_'//tochar(b3)
                call lo_h5_store_data(ps%plus(b1, b2, b3)%egv, file_id, trim(dum))
                dum = 'q_plus_'//tochar(b1)//'_'//tochar(b2)//'_'//tochar(b3)
                call lo_h5_store_data(ps%plus(b1, b2, b3)%q, file_id, trim(dum))
            end if
        end if
        nm = ps%minus(b1, b2, b3)%ntri
        dum = 'number_of_minus_triangles_'//tochar(b1)//'_'//tochar(b2)//'_'//tochar(b3)
        call lo_h5_store_attribute(nm, file_id, trim(dum))
        if (nm .gt. 0) then
            dum = 'triangles_minus_'//tochar(b1)//'_'//tochar(b2)//'_'//tochar(b3)
            call lo_h5_store_data(ps%minus(b1, b2, b3)%tri, file_id, trim(dum))
            dum = 'points_minus_'//tochar(b1)//'_'//tochar(b2)//'_'//tochar(b3)
            call lo_h5_store_data(ps%minus(b1, b2, b3)%pts, file_id, trim(dum))
            if (any(ps%plus(b1, b2, b3)%psisquared .gt. 0.0_flyt)) then
                dum = 'psisq_minus_'//tochar(b1)//'_'//tochar(b2)//'_'//tochar(b3)
                call lo_h5_store_data(ps%minus(b1, b2, b3)%psisquared, file_id, trim(dum))
                dum = 'omega_minus_'//tochar(b1)//'_'//tochar(b2)//'_'//tochar(b3)
                call lo_h5_store_data(ps%minus(b1, b2, b3)%omega, file_id, trim(dum))
                dum = 'egv_minus_'//tochar(b1)//'_'//tochar(b2)//'_'//tochar(b3)
                call lo_h5_store_data(ps%minus(b1, b2, b3)%egv, file_id, trim(dum))
                dum = 'q_minus_'//tochar(b1)//'_'//tochar(b2)//'_'//tochar(b3)
                call lo_h5_store_data(ps%minus(b1, b2, b3)%q, file_id, trim(dum))
            end if
        end if
    end do
    end do
    end do
    ! create a group for the BZ
    call h5gcreate_f(file_id, 'brillouin_zone', group_id, lo_status)
    if (lo_status .ne. 0) then
        write (*, *) '... error writing Brillouin zone'
    end if
    ! some meta
    call lo_h5_store_attribute(uc%bz%nnodes, group_id, 'number_of_nodes')
    call lo_h5_store_attribute(uc%bz%nfaces, group_id, 'number_of_faces')
    ! store the nodes
    call lo_h5_store_data(uc%bz%r, group_id, 'nodes')
    ! and the faces. Not the prettiest way, but it does not really matter.
    do i = 1, uc%bz%nfaces
        dum = 'face_'//tochar(i)
        call lo_h5_store_data(uc%bz%face(i)%ind, group_id, trim(dum))
    end do
    call h5gclose_f(group_id, lo_status)

    call h5fclose_f(file_id, lo_status); if (lo_status .ne. 0) stop
    call h5close_f(lo_status); if (lo_status .ne. 0) stop
end subroutine

end module
