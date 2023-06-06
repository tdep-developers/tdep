module activity
!! Estimates modes that are Raman or IR active
use konstanter, only: r8, lo_iou, lo_huge, lo_hugeint, lo_sqtol
use gottochblandat, only: lo_chop
use mpi_wrappers, only: lo_mpi_helper
use lo_memtracker, only: lo_mem_helper

use type_crystalstructure, only: lo_crystalstructure
use lo_symmetry_of_interactions, only: lo_interaction_tensors
use lo_phonon_bandstructure_on_path, only: lo_phonon_bandstructure

implicit none
private
public :: estimate_activity

contains

!> estimate Raman and IR activity, lowest order.
subroutine estimate_activity(p, bs, mw, mem, verbosity)
    !> structure
    type(lo_crystalstructure), intent(inout) :: p
    !> phonon dispersions along path
    type(lo_phonon_bandstructure), intent(inout) :: bs
    !> mpi helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem
    !> talk a lot?
    integer, intent(in) :: verbosity

    type(lo_interaction_tensors) :: sl
    real(r8), dimension(:, :, :), allocatable :: coeffRamanMode
    real(r8), dimension(:, :, :), allocatable :: coeffIRMode
    real(r8), dimension(:, :), allocatable :: coeffRaman
    real(r8), dimension(:, :), allocatable :: coeffIR

    ! Prepare some sensible tensors to play with.
    init: block
        real(r8), dimension(:, :), allocatable :: buf
        integer :: i, j, k, ia, ix, im, ishell, iop, nx
        !integer :: l
        ! Symmetry analysis to determine wether it could be Raman or IR active
        call sl%generate(p, p, p%mincutoff(), -1.0_r8, -1.0_r8, .true., mw, mem, verbosity)

        if (mw%talk) then
            write (lo_iou, *) ''
            write (lo_iou, *) 'Determining activity'
            if (sl%nx_eps_singlet .gt. 0) then
                write (lo_iou, *) '... Raman active'
            else
                write (lo_iou, *) '... not Raman active'
            end if
            if (sl%nx_Z_singlet .gt. 0) then
                write (lo_iou, *) '... IR active'
            else
                write (lo_iou, *) '... not IR active'
            end if
        end if

        ! Start with Raman activity. There is a magic realspace tensor that should be
        ! contracted with the eigenvectors: construct that one. But I will do it in the
        ! irreducible space, I think, so that we get a big fancy coefficient matrix
        ! to multiply with the eigenvectors? Yes/no/maybe.
        if (sl%nx_eps_singlet .gt. 0) then
            call mem%allocate(coeffRaman, [3*3*3*p%na, sl%nx_eps_singlet], &
                              persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%allocate(coeffRamanMode, [3*3, sl%nx_eps_singlet, p%na*3], &
                              persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%allocate(buf, [27, sl%nx_eps_singlet], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            coeffRaman = 0.0_r8
            coeffRamanMode = 0.0_r8
            buf = 0.0_r8

            ! Unfold the Raman coefficients to get a large coefficient matrix
            do ishell = 1, sl%n_eps_singlet_shell
            do k = 1, sl%eps_singlet_shell(ishell)%n_unfold

                ia = sl%eps_singlet_shell(ishell)%unfold_index(k)
                iop = sl%eps_singlet_shell(ishell)%unfold_operation(k)
                nx = sl%eps_singlet_shell(ishell)%nx
                if (nx .eq. 0) cycle
                ! Rotate coefficient matrix
                buf(:, 1:nx) = matmul(sl%tripletop(iop)%sotr, sl%eps_singlet_shell(ishell)%coeff)

                ! Store in the right place
                do i = 1, sl%eps_singlet_shell(ishell)%nx
                do j = 1, 27
                    ix = sl%eps_singlet_shell(ishell)%ind_global(i)
                    im = (ia - 1)*27 + j
                    coeffRaman(im, ix) = buf(j, i)
                end do
                end do
            end do
            end do

            if (mw%talk) write (lo_iou, *) '... transformed Raman coefficient'

            call mem%deallocate(buf, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        else
            call mem%allocate(coeffRaman, [1, 1], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%allocate(coeffRamanMode, [1, 1, 1], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            coeffRaman = 0.0_r8
            coeffRamanMode = 0.0_r8
        end if

        if (sl%nx_Z_singlet .gt. 0) then
            call mem%allocate(coeffIR, [3*3*p%na, sl%nx_Z_singlet], &
                              persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%allocate(coeffIRMode, [3, sl%nx_Z_singlet, 3*p%na], &
                              persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%allocate(buf, [9, sl%nx_Z_singlet], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            coeffIR = 0.0_r8
            coeffIRMode = 0.0_r8
            buf = 0.0_r8

            ! Unfold the Raman coefficients to get a large coefficient matrix
            do ishell = 1, sl%n_Z_singlet_shell
            do k = 1, sl%Z_singlet_shell(ishell)%n_unfold

                ia = sl%Z_singlet_shell(ishell)%unfold_index(k)
                iop = sl%Z_singlet_shell(ishell)%unfold_operation(k)
                nx = sl%Z_singlet_shell(ishell)%nx
                if (nx .eq. 0) cycle
                ! Rotate coefficient matrix
                buf(:, 1:nx) = matmul(sl%pairop(iop)%sotr, sl%Z_singlet_shell(ishell)%coeff)

                ! Store in the right place
                do i = 1, sl%Z_singlet_shell(ishell)%nx
                do j = 1, 9
                    ix = sl%Z_singlet_shell(ishell)%ind_global(i)
                    im = (ia - 1)*9 + j
                    coeffIR(im, ix) = buf(j, i)
                end do
                end do
            end do
            end do

            if (mw%talk) write (lo_iou, *) '... transformed IR coefficient'

            call mem%deallocate(buf, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        else
        end if
        ! if ( mw%talk ) then
        !     l=0
        !     do i=1,3
        !     do j=1,3
        !     do k=1,3
        !         l=l+1
        !         ia=(i-1)*9+(j-1)*3+k
        !         write(*,*) i,j,k,l,ia
        !     enddo
        !     enddo
        !     enddo
        ! endif
    end block init

    ! Sort it out, per q and mode, wether things are Raman or IR active.
    checkmodes: block
        integer :: imode, ia, ialpha, ibeta, igamma, imu, inu, ixi
        integer :: iq, jq

        ! First count the number of Gamma-points. Without any Gamma, no activity.
        bs%n_active_qpoint = 0
        do iq = 1, bs%n_point
            if (norm2(bs%q(iq)%r) .lt. lo_sqtol) bs%n_active_qpoint = bs%n_active_qpoint + 1
        end do

        ! If there are no active q-points, we are done here.
        if (bs%n_active_qpoint .eq. 0) then
            if (mw%talk) write (lo_iou, *) '... no zone-center q-points on path'
            call mem%deallocate(coeffRaman, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%deallocate(coeffRamanMode, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%deallocate(coeffIR, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%deallocate(coeffIRMode, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            return
        end if

        ! Store some data
        allocate (bs%active_qpoint(bs%n_active_qpoint))
        jq = 0
        do iq = 1, bs%n_point
            if (norm2(bs%q(iq)%r) .gt. lo_sqtol) cycle
            jq = jq + 1
            bs%active_qpoint(jq)%ind_path = iq
            bs%active_qpoint(jq)%nx_raman = sl%nx_eps_singlet
            bs%active_qpoint(jq)%nx_IR = sl%nx_Z_singlet
            if (bs%active_qpoint(jq)%nx_raman .gt. 0) then
                allocate (bs%active_qpoint(jq)%coeff_Raman(27*p%na, bs%active_qpoint(jq)%nx_raman))
                allocate (bs%active_qpoint(jq)%coeff_Raman_mode(9, bs%active_qpoint(jq)%nx_raman, 3*p%na))
                bs%active_qpoint(jq)%coeff_Raman = coeffRaman
                bs%active_qpoint(jq)%coeff_Raman_mode = 0.0_r8

                ! Contract coefficient matrix per mode.
                do imode = 1, p%na*3
                    do ia = 1, p%na
                    do ialpha = 1, 3
                    do ibeta = 1, 3
                    do igamma = 1, 3

                        imu = (ia - 1)*27 + (ialpha - 1)*9 + (ibeta - 1)*3 + igamma

                        inu = (ia - 1)*3 + igamma

                        ixi = (ialpha - 1)*3 + ibeta

                        bs%active_qpoint(jq)%coeff_Raman_mode(ixi, :, imode) = &
                            bs%active_qpoint(jq)%coeff_Raman_mode(ixi, :, imode) + &
                            bs%active_qpoint(jq)%coeff_Raman(imu, :)*real(bs%p(iq)%egv(inu, imode), r8)
                    end do
                    end do
                    end do
                    end do
                end do
                bs%active_qpoint(jq)%coeff_Raman_mode = lo_chop(bs%active_qpoint(jq)%coeff_Raman_mode, 1E-11_r8)
            end if

            if (bs%active_qpoint(jq)%nx_IR .gt. 0) then
                allocate (bs%active_qpoint(jq)%coeff_IR(9*p%na, bs%active_qpoint(jq)%nx_IR))
                allocate (bs%active_qpoint(jq)%coeff_IR_mode(3, bs%active_qpoint(jq)%nx_IR, 3*p%na))
                bs%active_qpoint(jq)%coeff_IR = coeffIR
                bs%active_qpoint(jq)%coeff_IR_mode = 0.0_r8

                ! Contract coefficient matrix per mode.
                do imode = 1, p%na*3
                    do ia = 1, p%na
                    do ialpha = 1, 3
                    do ibeta = 1, 3
                        ! This might very well be backwards
                        imu = (ia - 1)*9 + (ibeta - 1)*3 + ialpha
                        inu = (ia - 1)*3 + ibeta
                        ixi = (ibeta - 1)*3 + ialpha

                        bs%active_qpoint(jq)%coeff_IR_mode(ialpha, :, imode) = &
                            bs%active_qpoint(jq)%coeff_IR_mode(ialpha, :, imode) + &
                            bs%active_qpoint(jq)%coeff_IR(ixi, :)*real(bs%p(iq)%egv(inu, imode), r8)
                    end do
                    end do
                    end do
                end do
                bs%active_qpoint(jq)%coeff_IR_mode = lo_chop(bs%active_qpoint(jq)%coeff_IR_mode, 1E-11_r8)
            end if
        end do

    end block checkmodes

end subroutine

end module
