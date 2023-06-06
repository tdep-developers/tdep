submodule(type_phonon_dispersions) type_phonon_dispersions_generation
use konstanter, only: lo_iou
use lo_sorting, only: lo_qsort
implicit none
contains

!> Calculate dispersion relations on a q-point mesh
module subroutine generate(dr, qp, fc, p, mw, mem, verbosity)
    !> dispersion relations
    class(lo_phonon_dispersions), intent(out) :: dr
    !> qpoint mesh
    class(lo_qpoint_mesh), intent(inout) :: qp
    !> second order force constant
    type(lo_forceconstant_secondorder), intent(inout) :: fc
    !> crystal structure
    type(lo_crystalstructure), intent(inout) :: p
    !> MPI communicator
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem
    !> how much to talk
    integer, intent(in) :: verbosity

    real(r8) :: timer, t0, t1

    ! start timers
    timer = walltime()
    t0 = timer
    t1 = timer

    ! insert memory checkpoint
    call mem%tick()

    ! start by figuring out what to do
    initandheuristics: block
        integer :: q

        ! Check that all the info is consistent
        if (p%info%havespacegroup .eqv. .false.) then
            call lo_stop_gracefully(['Need spacegroup for generating phonon dispersions.'], &
                                    lo_exitcode_param, __FILE__, __LINE__, mw%comm)
        end if

        ! Say what we are going to do, perhaps?
        if (verbosity .gt. 0) then
            write (lo_iou, *) ''
            write (lo_iou, *) 'Calculating full dispersion relations over ', tochar(mw%n), ' MPI ranks'
        end if

        ! Some general things that are good to store.
        dr%n_mode = p%na*3
        dr%n_irr_qpoint = qp%n_irr_point
        dr%n_full_qpoint = qp%n_full_point

        ! Make some space
        allocate (dr%iq(qp%n_irr_point))
        do q = 1, qp%n_irr_point
            allocate (dr%iq(q)%omega(dr%n_mode))
            allocate (dr%iq(q)%vel(3, dr%n_mode))
            allocate (dr%iq(q)%egv(dr%n_mode, dr%n_mode))
            allocate (dr%iq(q)%degeneracy(dr%n_mode))
            allocate (dr%iq(q)%degenmode(dr%n_mode, dr%n_mode))
            dr%iq(q)%omega = 0.0_r8
            dr%iq(q)%vel = 0.0_r8
            dr%iq(q)%egv = 0.0_r8
            dr%iq(q)%degeneracy = 0
            dr%iq(q)%degenmode = 0
        end do
        allocate (dr%aq(qp%n_full_point))
        do q = 1, qp%n_full_point
            allocate (dr%aq(q)%omega(dr%n_mode))
            allocate (dr%aq(q)%vel(3, dr%n_mode))
            allocate (dr%aq(q)%egv(dr%n_mode, dr%n_mode))
            allocate (dr%aq(q)%degeneracy(dr%n_mode))
            allocate (dr%aq(q)%degenmode(dr%n_mode, dr%n_mode))
            dr%aq(q)%omega = 0.0_r8
            dr%aq(q)%vel = 0.0_r8
            dr%aq(q)%egv = 0.0_r8
            dr%aq(q)%degeneracy = 0
            dr%aq(q)%degenmode = 0
        end do
    end block initandheuristics

    ! Now calculate the dispersions in the irreducible wedge
    irrwedge: block
        integer :: q, qc, ctr

        ! Counter for progressbar
        ctr = 0
        do q = 1, qp%n_irr_point
            if (mod(q, mw%n) .eq. mw%r) ctr = ctr + 1
        end do

        if (verbosity .gt. 0) call lo_progressbar_init()
        qc = 0
        do q = 1, qp%n_irr_point
            ! Make it parallel
            if (mod(q, mw%n) .ne. mw%r) cycle
            ! Get dispersions
            call dr%iq(q)%generate(fc, p, mem, qp%ip(q))
            if (verbosity .gt. 0) then
                qc = qc + 1
                if (lo_trueNtimes(qc, 127, ctr)) call lo_progressbar(' ... irreducible dispersions', qc, ctr, walltime() - t0)
            end if
        end do
        ! And distribute it
        call dr%allgather_irreducible(mw, mem)
        ! maybe say that we are done?
        if (verbosity .gt. 0) then
            t1 = walltime()
            call lo_progressbar(' ... irreducible dispersions', ctr, ctr, t1 - t0)
            t0 = t1
        end if
    end block irrwedge

    ! Now get the rest of the dispersions. That is as simple as unfolding it.
    fullzone: block
        complex(r8), dimension(:, :), allocatable :: rotmat
        integer :: j, q, qc, ctr, iqp, iop

        ! space for rotation matrix
        call mem%allocate(rotmat, [dr%n_mode, dr%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        rotmat = 0.0_r8

        ! Counter for progressbar
        ctr = 0
        do q = 1, qp%n_full_point
            if (mod(q, mw%n) .eq. mw%r) ctr = ctr + 1
        end do

        if (verbosity .gt. 0) call lo_progressbar_init()
        qc = 0
        do q = 1, qp%n_full_point
            ! Make it parallel
            if (mod(q, mw%n) .ne. mw%r) cycle
            ! for each band, get the infomation from the irreducible wedge and transform it here
            ! which irreducible q-point is this
            iqp = qp%ap(q)%irreducible_index
            ! what's the operation that takes the point to the wedge
            ! positive kk means it's a straight operation, negative means it's under reversal
            iop = qp%ap(q)%operation_from_irreducible
            if (iop .gt. 0) then
                ! the eigenvector transformation matrix
                call lo_eigenvector_transformation_matrix(rotmat, p%rcart, qp%ip(iqp)%r, p%sym%op(abs(iop)))
                ! rotate the eigenvectors
                call lo_gemm(rotmat, dr%iq(iqp)%egv, dr%aq(q)%egv)
                do j = 1, dr%n_mode
                    dr%aq(q)%omega(j) = dr%iq(iqp)%omega(j)
                    ! Make sure the velocities have the proper symmetry
                    dr%aq(q)%vel(:, j) = lo_operate_on_vector(p%sym%op(abs(iop)), dr%iq(iqp)%vel(:, j), reciprocal=.true.)
                end do
            else
                ! Same thing, but now we have to time-reverse it
                call lo_eigenvector_transformation_matrix(rotmat, p%rcart, qp%ip(iqp)%r, p%sym%op(abs(iop)))
                ! rotate the eigenvectors
                call lo_gemm(rotmat, dr%iq(iqp)%egv, dr%aq(q)%egv)
                ! and conjugate them, that's what Maradudin says.
                dr%aq(q)%egv = conjg(dr%aq(q)%egv)
                do j = 1, dr%n_mode
                    dr%aq(q)%omega(j) = dr%iq(iqp)%omega(j)
                    ! Make sure the velocities have the proper symmetry. Note negative sign here.
                    dr%aq(q)%vel(:, j) = -lo_operate_on_vector(p%sym%op(abs(iop)), dr%iq(iqp)%vel(:, j), reciprocal=.true.)
                end do
            end if
            ! The degeneracy I can just copy
            dr%aq(q)%degeneracy = dr%iq(iqp)%degeneracy
            dr%aq(q)%degenmode = dr%iq(iqp)%degenmode

            if (verbosity .gt. 0) then
                qc = qc + 1
                if (lo_trueNtimes(qc, 127, ctr)) call lo_progressbar(' ... full dispersions', qc, ctr, walltime() - t0)
            end if
        end do

        ! and communicate it
        call dr%allgather_fullmesh(mw, mem)
        ! Cleanup
        call mem%deallocate(rotmat, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)

        if (verbosity .gt. 0) then
            t1 = walltime()
            call lo_progressbar(' ... full dispersions', qp%n_full_point, qp%n_full_point, t1 - t0)
            t0 = t1
        end if
    end block fullzone

    ! Get the default smearing parameters, per band
    baselinesmearing: block
        real(r8), dimension(:), allocatable :: x
        real(r8) :: f0
        integer :: imode, iqp

        ! Find min and max frequency, always useful.
        dr%omega_max = -lo_huge
        dr%omega_min = lo_huge
        do iqp = 1, dr%n_irr_qpoint
        do imode = 1, dr%n_mode
            f0 = dr%iq(iqp)%omega(imode)
            if (abs(f0) .gt. lo_freqtol) dr%omega_min = min(dr%omega_min, f0)
            dr%omega_max = max(dr%omega_max, f0)
        end do
        end do

        ! Get the default smearing parameter for Gaussian integrations, per band.
        allocate (dr%default_smearing(dr%n_mode))
        dr%default_smearing = 0.0_r8

        call mem%allocate(x, dr%n_irr_qpoint, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        x = 0.0_r8
        do imode = 1, dr%n_mode
            if (mod(imode, mw%n) .ne. mw%r) cycle
            ! get the energies for this band
            do iqp = 1, dr%n_irr_qpoint
                x(iqp) = dr%iq(iqp)%omega(imode)
            end do
            ! sort it?
            call lo_qsort(x)
            ! find the largest separation
            f0 = 0.0_r8
            do iqp = 1, dr%n_irr_qpoint - 1
                f0 = max(f0, abs(x(iqp + 1) - x(iqp)))
            end do
            dr%default_smearing(imode) = f0
        end do
        call mw%allreduce('max', dr%default_smearing)
        call mem%deallocate(x, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)

        ! Now, if I havey very flat bands this can become slightly problematic
        ! and the DOS gets extremely spiky. So I set some sensible minimum
        ! value for the smearing.
        f0 = maxval(dr%default_smearing)/5.0_r8
        dr%default_smearing = max(dr%default_smearing, f0)
    end block baselinesmearing

    ! Check memory
    call mem%tock(__FILE__, __LINE__, mw%comm)

    if (verbosity .gt. 0) then
        write (lo_iou, *) 'Calculated the full dispersion relations in ', tochar(walltime() - timer), 's'
    end if
end subroutine

!> Distribute a phonon dispersion object over MPI. Only distributes the data in the irreducible wedge.
module subroutine allgather_irreducible(dr, mw, mem)
    !> phonon dispersions
    class(lo_phonon_dispersions), intent(inout) :: dr
    !> mpi helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem

    character, dimension(:), allocatable :: sbuf, rbuf
    integer :: i, j, ctr, ctr_tot, pos

    ! Start counting bytes per rank that should be sent.
    ctr = 0
    do i = 1, dr%n_irr_qpoint
        if (mod(i, mw%n) .ne. mw%r) cycle
        ! first want I pack the q-index
        ctr = ctr + storage_size(i)/8
        ! then the rest of the q-point
        ctr = ctr + dr%iq(i)%size_packed()
    end do

    ! Get size per rank and offset
    call mw%size_and_offset(ctr, ctr_tot)

    ! Space for communication buffers
    if (ctr .gt. 0) then
        call mem%allocate(sbuf, ctr, persistent=.false., scalable=.true., file=__FILE__, line=__LINE__)
    else
        call mem%allocate(sbuf, 1, persistent=.false., scalable=.true., file=__FILE__, line=__LINE__)
    end if
    call mem%allocate(rbuf, ctr_tot, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)

    ! Pack things into send buffer
    pos = 0
    do i = 1, dr%n_irr_qpoint
        if (mod(i, mw%n) .ne. mw%r) cycle
        ! pack the index
        call mw%pack(i, sbuf, pos)
        ! pack the point
        call dr%iq(i)%pack_to_buf(sbuf, pos, mw)
    end do
    ! Now I can commence operation actual operation
    call mw%allgatherv(sbuf, rbuf, mw%ctr_size_per_rank, mw%ctr_offset_per_rank, __FILE__, __LINE__)
    ! Clear the send buffer
    call mem%deallocate(sbuf, persistent=.false., scalable=.true., file=__FILE__, line=__LINE__)
    ! Start unpacking
    pos = 0
    do j = 1, dr%n_irr_qpoint
        ! fetch which q-index is stored here
        call mw%unpack(i, rbuf, pos)
        ! unpack the point
        call dr%iq(i)%unpack_from_buf(rbuf, pos, mw)
    end do
    ! Clear the recieve buffer
    call mem%deallocate(rbuf, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
end subroutine

!> Distribute a phonon dispersion object over MPI. Distributes data on the full grid.
module subroutine allgather_fullmesh(dr, mw, mem)
    !> phonon dispersions
    class(lo_phonon_dispersions), intent(inout) :: dr
    !> mpi helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem

    character, dimension(:), allocatable :: sbuf, rbuf
    integer :: i, j, ctr, ctr_tot, pos

    ! Start counting bytes per rank that should be sent.
    ctr = 0
    do i = 1, dr%n_full_qpoint
        if (mod(i, mw%n) .ne. mw%r) cycle
        ! first want I pack the q-index
        ctr = ctr + storage_size(i)/8
        ! then the rest of the q-point
        ctr = ctr + dr%aq(i)%size_packed()
    end do

    ! Get size per rank and offset
    call mw%size_and_offset(ctr, ctr_tot)

    ! Space for communication buffers
    if (ctr .gt. 0) then
        call mem%allocate(sbuf, ctr, persistent=.false., scalable=.true., file=__FILE__, line=__LINE__)
    else
        call mem%allocate(sbuf, 1, persistent=.false., scalable=.true., file=__FILE__, line=__LINE__)
    end if
    call mem%allocate(rbuf, ctr_tot, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)

    ! Pack things into send buffer
    pos = 0
    do i = 1, dr%n_full_qpoint
        if (mod(i, mw%n) .ne. mw%r) cycle
        ! pack the index
        call mw%pack(i, sbuf, pos)
        ! pack the point
        call dr%aq(i)%pack_to_buf(sbuf, pos, mw)
    end do

    ! Now I can commence operation actual operation
    call mw%allgatherv(sbuf, rbuf, mw%ctr_size_per_rank, mw%ctr_offset_per_rank, __FILE__, __LINE__)

    ! Clear the send buffer
    call mem%deallocate(sbuf, persistent=.false., scalable=.true., file=__FILE__, line=__LINE__)
    ! Start unpacking
    pos = 0
    do j = 1, dr%n_full_qpoint
        ! fetch which q-index is stored here
        call mw%unpack(i, rbuf, pos)
        ! unpack the point
        call dr%aq(i)%unpack_from_buf(rbuf, pos, mw)
    end do
    ! Clear the recieve buffer
    call mem%deallocate(rbuf, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
end subroutine

! !> Calculate mode gruneisen parameters on a q-point mesh
! module subroutine gruneisen(dr,qp,fct,p,verbosity,mpi_communicator)
!     !> the dispersion relations
!     class(lo_phonon_dispersions), intent(inout) :: dr
!     !> the qpoint mesh
!     class(lo_qpoint_mesh), intent(inout) :: qp
!     !> the third order force constant
!     type(lo_forceconstant_thirdorder), intent(in) :: fct
!     !> the crystal structure
!     type(lo_crystalstructure), intent(in) :: p
!     !> how much to talk
!     integer, intent(in), optional :: verbosity
!     !> do the calculation distributed
!     integer, intent(in), optional :: mpi_communicator
!     !
!     write(*,*) 'FIXME GRUNEISEN MPI'
!     stop
!
!     integer :: q,rank,nrank,lqp
!     real(r8) :: t0
!     logical :: usempi
!
!     t0=walltime()
!     if ( present(verbosity) ) then
!        verbosity=verbosity
!     else
!        verbosity=0
!     endif
!
!     ! If we have an MPI communicator attached, check that the q-points are distributed accordingly.
!     if ( present(mpi_communicator) ) then
!         ! Are the q-points distributed across a communicator?
!         if ( qp%mpi%initialized ) then
!             ! Yep, distributed. The same communicator that is sent here?
!             if ( qp%mpi%communicator .ne. mpi_communicator ) then
!                 ! This is weird. Not sure what to do. Probably redistribute?
!                 call qp%divide_across_mpi(mpi_communicator)
!             endif
!         else
!             ! Nope, not initialized. Do that!
!             call qp%divide_across_mpi(mpi_communicator)
!         endif
!         usempi=.true.
!         call mpi_comm_rank(mpi_communicator,rank,lo_status)
!         call mpi_comm_size(mpi_communicator,nrank,lo_status)
!         ! remove verbosity except on head rank
!         if ( rank .ne. 0 ) verbosity=0
!         ! if only one rank, don't use MPI at all.
!         if ( nrank .eq. 1 ) usempi=.false.
!     else
!         ! Don't use MPI at all.
!         usempi=.false.
!         rank=0
!         nrank=1
!     endif
!
!     if ( verbosity .gt. 0 ) then
!         if ( usempi ) then
!             write(*,*) '... calculating gruneisen parameters over ',tochar(nrank),' MPI ranks'
!         else
!             write(*,*) '... calculating gruneisen parameters'
!         endif
!     endif
!
!     ! Make some space
!     do q=1,qp%n_irr_point
!         lo_allocate(dr%iq(q)%gruneisen(dr%n_mode))
!         dr%iq(q)%gruneisen=0.0_r8
!     enddo
!     do q=1,qp%n_full_point
!         lo_allocate(dr%aq(q)%gruneisen(dr%n_mode))
!         dr%aq(q)%gruneisen=0.0_r8
!     enddo
!
!     ! And the actual calculation
!     if ( usempi ) then
!         do lqp=1,qp%mpi%nqirr
!             q=qp%mpi%qirr(lqp)
!             call fct%mode_gruneisen_parameter(p,qp%ip(q),dr%iq(q)%omega,dr%iq(q)%egv,dr%iq(q)%gruneisen)
!         enddo
!         call dr%allgather_irreducible(qp,mpi_communicator,fieldlist=[4])
!     else
!         do q=1,qp%n_irr_point
!             call fct%mode_gruneisen_parameter(p,qp%ip(q),dr%iq(q)%omega,dr%iq(q)%egv,dr%iq(q)%gruneisen)
!         enddo
!     endif
!
!     ! and spread them to the full thing
!     do q=1,qp%n_full_point
!         dr%aq(q)%gruneisen=dr%iq( qp%ap(q)%irrind )%gruneisen
!     enddo
!     if ( verbosity .gt. 0 ) then
!         write(*,*) '... got gruneisen parameters in irreducible wedge in '//tochar(walltime()-t0,6)//' s'
!         t0=walltime()
!     endif
! end subroutine

end submodule
