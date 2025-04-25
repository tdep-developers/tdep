module scf_helper
!!
!! Common helper routines.
!!
use konstanter, only: r8, i8, lo_imag, lo_pi, lo_hugeint, lo_huge, lo_iou, lo_exitcode_param, lo_freqtol, lo_tiny, lo_frequency_Hartree_to_THz
use gottochblandat, only: lo_gauss, lo_trapezoid_integration, lo_planck, lo_return_unique, lo_linspace, open_file, lo_linear_interpolation
use hdf5_wrappers, only: lo_hdf5_helper
use mpi_wrappers, only: lo_mpi_helper, lo_stop_gracefully
use lo_memtracker, only: lo_mem_helper
use type_qpointmesh, only: lo_qpoint_mesh, lo_fft_mesh, lo_qpoint_irrwedge
use type_phonon_dispersions, only: lo_phonon_dispersions
use lo_brents_method, only: lo_brent_helper
use quadratures_stencils, only: lo_centraldifference, lo_gaussianquadrature
use lo_sorting, only: lo_qsort
use type_crystalstructure, only: lo_crystalstructure
use type_forceconstant_secondorder, only: lo_forceconstant_secondorder
use type_forcemap, only: lo_forcemap
use lineshape_helper, only: lo_spectralfunction_helper
use ifc_solvers, only: lo_irreducible_forceconstant_from_qmesh_dynmat

implicit none

private
public :: return_new_bare_phonons

contains

subroutine return_new_bare_phonons(qp,dr,sf,map,uc,fc,mw,mem,verbosity)
    !> q-point mesh
    class(lo_qpoint_mesh), intent(in) :: qp
    !> harmonic properties on this mesh
    type(lo_phonon_dispersions), intent(in) :: dr
    !> spectralfunction helper
    type(lo_spectralfunction_helper), intent(in) :: sf
    !> forcemap
    type(lo_forcemap), intent(inout) :: map
    !> structure
    type(lo_crystalstructure), intent(inout) :: uc
    !> current forceconstant
    type(lo_forceconstant_secondorder), intent(inout) :: fc
    !> mpi helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem
    !> talk a lot?
    integer, intent(in) :: verbosity

    complex(r8), dimension(:,:,:), allocatable :: dynmat,lrdynmat,sedynmat

    allocate(dynmat(3*uc%na,3*uc%na,qp%n_irr_point))
    allocate(lrdynmat(3*uc%na,3*uc%na,qp%n_irr_point))
    allocate(sedynmat(3*uc%na,3*uc%na,qp%n_irr_point))
    dynmat=0.0_r8
    lrdynmat=0.0_r8
    sedynmat=0.0_r8

    standarddm: block
        complex(r8), dimension(:,:,:,:), allocatable :: wdm
        integer :: iq,a1,a2,i,j,ii,jj

        allocate(wdm(3,3,uc%na,uc%na))

        do iq=1,qp%n_irr_point
            if ( mod(iq,mw%n) .ne. mw%r ) cycle
            ! Get the normal dynamical matrix, skipping the non-analytical
            call fc%dynamicalmatrix(uc,qp%ip(iq),dynmat(:,:,iq),mem,skipnonanalytical=.true.)
            ! Get only the long-range
            if ( map%polar .ne. 0 ) then
                call fc%longrange_dynamical_matrix(wdm, uc, qp%ip(iq)%r)
                do a1 = 1, uc%na
                do a2 = 1, uc%na
                    do i = 1, 3
                    do j = 1, 3
                        ii = (a1 - 1)*3 + i
                        jj = (a2 - 1)*3 + j
                        lrdynmat(jj, ii, iq) = wdm(i, j, a1, a2)
                    end do
                    end do
                end do
                end do
            endif
        enddo
        call mw%allreduce('sum',dynmat)
        call mw%allreduce('sum',lrdynmat)
    end block standarddm

    ! Needs to be done serially for now, I think that's fine.

    selfenergycorrection: block
        complex(r8), dimension(:,:), allocatable :: cm0,cm1,cm2
        complex(r8), dimension(:,:,:,:), allocatable :: wdm
        real(r8), dimension(:), allocatable :: x0
        real(r8) :: f0
        integer :: iq,imode

        sedynmat=0.0_r8
        if ( mw%r .eq. mw%n-1 ) then
            allocate(cm0(3*uc%na,3*uc%na))
            allocate(cm1(3*uc%na,3*uc%na))
            allocate(cm2(3*uc%na,3*uc%na))

            allocate(x0(sf%n_energy))
            call lo_linspace(0.0_r8,sf%max_energy,x0)

            do iq=1,qp%n_irr_point
                cm0=0.0_r8
                cm1=0.0_r8
                do imode=1,dr%n_mode
                    cm0(:,imode)=dr%iq(iq)%egv(:,imode)*sqrt( dr%iq(iq)%omega(imode) )
                    ! Let's pick the static one for now. Can adjust later.
                    ! Instead of the static,
                    !cm1(imode,imode)=sf%sigma_Re(1,imode,iq)
                    f0=lo_linear_interpolation(x0,sf%sigma_Re(:,imode,iq),dr%iq(iq)%omega(imode))
                    cm1(imode,imode)=f0 !sf%sigma_Re(1,imode,iq)
                enddo
                !cm0=transpose(conjg(cm0))
                cm2=matmul(cm0,cm1)
                sedynmat(:,:,iq)=matmul(cm2,transpose(conjg(cm0)))
                !sedynmat(:,:,iq)=matmul(cm2,transpose(cm0))
            enddo
        endif
        call mw%allreduce('sum',sedynmat)
        dynmat=dynmat + sedynmat
    end block selfenergycorrection

    ! Solve for new IFCs
    newfc: block
        type(lo_qpoint_irrwedge), dimension(:), allocatable :: dummyq
        complex(r8), dimension(:,:,:), allocatable :: wdm,wlr
        integer :: i,nq

        nq=0
        do i=1,qp%n_irr_point
            if ( mod(i,mw%n) .ne. mw%r ) cycle
            nq=nq+1
        enddo

        if ( nq .gt. 0 ) then
            allocate(wdm(3*uc%na,3*uc%na,nq))
            allocate(wlr(3*uc%na,3*uc%na,nq))
            allocate(dummyq(nq))
            wdm=0.0_r8
            wlr=0.0_r8
        endif

        nq=0
        do i=1,qp%n_irr_point
            if ( mod(i,mw%n) .ne. mw%r ) cycle
            nq=nq+1
            wdm(:,:,nq)=dynmat(:,:,i)
            wlr(:,:,nq)=lrdynmat(:,:,i)
            dummyq(nq)=qp%ip(i)
        enddo

        deallocate(dynmat)
        deallocate(lrdynmat)
        deallocate(sedynmat)

        call lo_irreducible_forceconstant_from_qmesh_dynmat(map,uc,dummyq,nq,wdm,wlr,.true.,.false.,mw,verbosity+2)

    end block newfc

end subroutine

end module
