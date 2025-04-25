
!> write thingy to hdf5
subroutine write_to_hdf5(ise, filename)
    !> self-energy guy
    class(lo_interpolated_selfenergy), intent(in) :: ise
    !> filename
    character(len=*), intent(in) :: filename

    type(lo_hdf5_helper) :: h5
    real(r8), dimension(:, :), allocatable :: dr0, dr1
    integer, dimension(:, :), allocatable :: di
    integer :: ipair, nx, i, l

    call h5%init(__FILE__, __LINE__)
    call h5%open_file('write', trim(filename))

    ! Store simple things
    call h5%store_attribute(ise%n_atom, h5%file_id, 'n_atom')
    call h5%store_attribute(ise%n_x, h5%file_id, 'n_x')
    call h5%store_attribute(ise%n_pair, h5%file_id, 'n_pair')
    call h5%store_attribute(ise%n_energy, h5%file_id, 'n_energy')
    call h5%store_attribute(ise%temperature, h5%file_id, 'temperature')
    call h5%store_attribute(ise%omega_shift, h5%file_id, 'omega_shift')

    call h5%store_data(ise%energy, h5%file_id, 'energy')
    call h5%store_data(ise%x_re, h5%file_id, 'x_re')
    call h5%store_data(ise%x_im, h5%file_id, 'x_im')

    allocate (di(3, ise%n_pair))
    nx = -1
    do ipair = 1, ise%n_pair
        di(:, ipair) = [ise%pair(ipair)%a1, ise%pair(ipair)%a2, ise%pair(ipair)%nx]
        nx = max(nx, ise%pair(ipair)%nx)
    end do
    call h5%store_data(di, h5%file_id, 'pair_ind')
    deallocate (di)

    allocate (di(nx, ise%n_pair))
    allocate (dr0(nx, 9*ise%n_pair))
    allocate (dr1(3, ise%n_pair))
    di = 0
    dr0 = 0.0_r8
    dr1 = 0.0_r8

    l = 0
    do ipair = 1, ise%n_pair
        if (ise%pair(ipair)%nx .eq. 0) cycle
        di(1:ise%pair(ipair)%nx, ipair) = ise%pair(ipair)%xind
        dr1(:, ipair) = ise%pair(ipair)%lv
        do i = 1, 9
            l = l + 1
            dr0(1:ise%pair(ipair)%nx, l) = ise%pair(ipair)%coeff(i, :)
        end do
    end do
    call h5%store_data(di, h5%file_id, 'x_ind')
    call h5%store_data(dr0, h5%file_id, 'coeff')
    call h5%store_data(dr1, h5%file_id, 'latvec')

    deallocate (di)
    deallocate (dr0)
    deallocate (dr1)

    call h5%close_file()
    call h5%destroy(__FILE__, __LINE__)
end subroutine

!> write thingy to hdf5
subroutine read_from_hdf5(ise, filename, mw)
    !> self-energy guy
    class(lo_interpolated_selfenergy), intent(out) :: ise
    !> filename
    character(len=*), intent(in) :: filename
    !> mpi helper
    type(lo_mpi_helper), intent(inout) :: mw

    type(lo_hdf5_helper) :: h5
    real(r8), dimension(:, :), allocatable :: dr0, dr1
    integer, dimension(:, :), allocatable :: di, dj
    integer :: readrnk, nx, ipair, i, l

    readrnk = mw%n - 1

    ! Open on head rank
    if (mw%r .eq. readrnk) then
        call h5%init(__FILE__, __LINE__)
        call h5%open_file('read', trim(filename))

        ! Read some basics
        call h5%read_attribute(ise%n_atom, h5%file_id, 'n_atom')
        call h5%read_attribute(ise%n_x, h5%file_id, 'n_x')
        call h5%read_attribute(ise%n_pair, h5%file_id, 'n_pair')
        call h5%read_attribute(ise%n_energy, h5%file_id, 'n_energy')
        call h5%read_attribute(ise%temperature, h5%file_id, 'temperature')
        call h5%read_attribute(ise%omega_shift, h5%file_id, 'omega_shift')
        call h5%read_data(ise%energy, h5%file_id, 'energy')
        call h5%read_data(ise%x_re, h5%file_id, 'x_re')
        call h5%read_data(ise%x_im, h5%file_id, 'x_im')

        call h5%read_data(di, h5%file_id, 'pair_ind')
        call h5%read_data(dj, h5%file_id, 'x_ind')
        call h5%read_data(dr0, h5%file_id, 'coeff')
        call h5%read_data(dr1, h5%file_id, 'latvec')
        nx = size(dr0, 1)

        call h5%close_file()
        call h5%destroy(__FILE__, __LINE__)
    end if

    ! Spread metadata around a bit
    call mw%bcast(ise%n_atom, from=readrnk)
    call mw%bcast(ise%n_x, from=readrnk)
    call mw%bcast(ise%n_pair, from=readrnk)
    call mw%bcast(ise%n_energy, from=readrnk)
    call mw%bcast(ise%temperature, from=readrnk)
    call mw%bcast(nx, from=readrnk)

    if (mw%r .ne. readrnk) then
        allocate (di(3, ise%n_pair))
        allocate (dj(nx, ise%n_pair))
        allocate (dr0(nx, 9*ise%n_pair))
        allocate (dr1(3, ise%n_pair))

        allocate (ise%energy(ise%n_energy))
        allocate (ise%x_re(ise%n_x, ise%n_energy))
        allocate (ise%x_im(ise%n_x, ise%n_energy))
    end if

    call mw%bcast(di, from=readrnk)
    call mw%bcast(dj, from=readrnk)
    call mw%bcast(dr0, from=readrnk)
    call mw%bcast(dr1, from=readrnk)
    call mw%bcast(ise%energy, from=readrnk)
    call mw%bcast(ise%x_re, from=readrnk)
    call mw%bcast(ise%x_im, from=readrnk)

    ! Make space
    allocate (ise%pair(ise%n_pair))

    l = 0
    do ipair = 1, ise%n_pair
        ise%pair(ipair)%a1 = di(1, ipair)
        ise%pair(ipair)%a2 = di(2, ipair)
        ise%pair(ipair)%nx = di(3, ipair)
        ise%pair(ipair)%lv = dr1(:, ipair)
        if (ise%pair(ipair)%nx .eq. 0) cycle
        allocate (ise%pair(ipair)%xind(ise%pair(ipair)%nx))
        allocate (ise%pair(ipair)%coeff(9, ise%pair(ipair)%nx))
        ise%pair(ipair)%xind = 0.0_r8
        ise%pair(ipair)%coeff = 0.0_r8
        ise%pair(ipair)%xind = dj(1:ise%pair(ipair)%nx, ipair)
        do i = 1, 9
            l = l + 1
            ise%pair(ipair)%coeff(i, :) = dr0(1:ise%pair(ipair)%nx, l)
        end do
    end do

    deallocate (di)
    deallocate (dj)
    deallocate (dr0)
    deallocate (dr1)
end subroutine
