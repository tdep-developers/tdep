module lo_fftgrid_helper
use konstanter, only: r8
use lo_randomnumbers, only: lo_mersennetwister

implicit none
private
public :: lo_montecarlo_grid
public :: singlet_to_triplet
public :: triplet_to_singlet
public :: fft_third_grid_index
public :: fft_fourth_grid_index

type lo_montecarlo_grid
    !> The size of the grid
    integer :: npoints
    !> The weight of each point on the Monte-Carlo grid
    real(r8) :: weight
    !> The dimensions of the Monte-Carlo grid
    integer, dimension(3) :: mc_dims
    !> The dimensions of the full grid
    integer, dimension(3) :: full_dims
    !> The ratio between the full and mc grids
    real(r8), dimension(3) :: ratio

contains
    !> Initialize the grid
    procedure :: initialize => initialize_montecarlo_grid
    !> Generate a point from the Monte-Carlo to the full grid
    procedure :: mc_point_to_full
    !> Generate an array with the grid
    procedure :: generate_grid
end type

contains

subroutine initialize_montecarlo_grid(mcg, full_dims, mc_dims)
    !> The Monte-Carlo grid
    class(lo_montecarlo_grid), intent(out) :: mcg
    !> The dimensions of the full grid
    integer, dimension(3), intent(in) :: full_dims
    !> The dimensions of the Monte-Carlo grid
    integer, dimension(3), intent(in) :: mc_dims

    !> Some integer for the do loop
    integer :: i

    mcg%full_dims = full_dims
    do i=1, 3
        mcg%mc_dims(i) = min(mc_dims(i), mcg%full_dims(i))
        mcg%ratio(i) = real(mcg%full_dims(i), r8) / real(mcg%mc_dims(i), r8)
    end do
    mcg%npoints = mcg%mc_dims(1) * mcg%mc_dims(2) * mcg%mc_dims(3)
    mcg%weight = 1.0_r8 / real(mcg%npoints, r8)
end subroutine

function mc_point_to_full(mcg, imc, rng) result(ifull)
    !> The Monte-Carlo grid
    class(lo_montecarlo_grid), intent(in) :: mcg
    !> The index of the point on the Monte-Carlo grid
    integer, intent(in) :: imc
    !> The random number generator
    type(lo_mersennetwister), intent(inout) :: rng
    !> The index of the point on the full grid
    integer :: ifull

    !> The triplets of point on the Monte-Carlo and full grids
    integer, dimension(3) :: gi_mc, gi_full

    gi_mc = singlet_to_triplet(imc, mcg%mc_dims(2), mcg%mc_dims(3))
    ! This way of generating the number makes it work even if the ratio is not an integer
    gi_full(1) = ceiling((real(gi_mc(1), r8) - 1.0_r8) * mcg%ratio(1) + rng%rnd_real() * mcg%ratio(1))
    gi_full(2) = ceiling((real(gi_mc(2), r8) - 1.0_r8) * mcg%ratio(2) + rng%rnd_real() * mcg%ratio(2))
    gi_full(3) = ceiling((real(gi_mc(3), r8) - 1.0_r8) * mcg%ratio(3) + rng%rnd_real() * mcg%ratio(3))
    ifull = triplet_to_singlet(gi_full, mcg%full_dims(2), mcg%full_dims(3))
end function

subroutine generate_grid(mcg, qgrid, rng)
    !> The Monte-Carlo grid
    class(lo_montecarlo_grid), intent(in) :: mcg
    !> The random number generator
    type(lo_mersennetwister), intent(inout) :: rng
    !> The grid to be generated
    integer, dimension(:), intent(out) :: qgrid

    !> Some integers for the do loop
    integer :: qi, qprev, qtest

    ! To improve convergence, we avoid repeating points in the integration grid
    qgrid(1) = mcg%mc_point_to_full(1, rng)
    qprev = qgrid(1)
    qtest = qgrid(1)
    do qi=2, mcg%npoints
        do while(qtest .eq. qprev)
            qtest = mcg%mc_point_to_full(qi, rng)
        end do
        qgrid(qi) = qtest
        qprev = qtest
    end do
end subroutine

!> convert a linear index to a triplet
pure function singlet_to_triplet(l, ny, nz) result(gi)
    !> linear index
    integer, intent(in) :: l
    !> second dimension
    integer, intent(in) :: ny
    !> third dimension
    integer, intent(in) :: nz
    !> grid-index
    integer, dimension(3) :: gi

    integer :: i, j, k

    k = mod(l, nz)
    if (k .eq. 0) k = nz
    j = mod((l - k)/nz, ny) + 1
    i = (l - k - (j - 1)*nz)/(nz*ny) + 1
    gi = [i, j, k]
end function

!> convert a triplet index to a singlet
pure function triplet_to_singlet(gi, ny, nz) result(l)
    !> grid-index
    integer, dimension(3), intent(in) :: gi
    !> second dimension
    integer, intent(in) :: ny
    !> third dimension
    integer, intent(in) :: nz
    !> linear index
    integer :: l

    l = (gi(1) - 1)*ny*nz + (gi(2) - 1)*nz + gi(3)
end function

!> returns the index on the grid that gives q3=-q1-q2
pure function fft_third_grid_index(i1, i2, dims) result(i3)
    !> index to q1
    integer, intent(in) :: i1
    !> index to q2
    integer, intent(in) :: i2
    !> dimensions of the grid
    integer, dimension(3), intent(in) :: dims
    !> index to q3
    integer :: i3

    integer, dimension(3) :: gi1, gi2, gi3
    integer :: l, k

    ! Convert triplet to singlet
    gi1 = singlet_to_triplet(i1, dims(2), dims(3))
    gi2 = singlet_to_triplet(i2, dims(2), dims(3))
    do l = 1, 3
        gi3(l) = 3 - gi1(l) - gi2(l)
    end do
    do k = 1, 3
    do l = 1, 3
        if (gi3(l) .lt. 1) gi3(l) = gi3(l) + dims(l)
        if (gi3(l) .gt. dims(l)) gi3(l) = gi3(l) - dims(l)
    end do
    end do
    ! convert it back to a singlet
    i3 = triplet_to_singlet(gi3, dims(2), dims(3))
end function

!> returns the index on the grid that gives q4=-q3-q2-q1
pure function fft_fourth_grid_index(i1, i2, i3, dims) result(i4)
    !> index to q1
    integer, intent(in) :: i1
    !> index to q2
    integer, intent(in) :: i2
    !> index to q3
    integer, intent(in) :: i3
    !> dimensions of the grid
    integer, dimension(3), intent(in) :: dims
    !> index to q4
    integer :: i4

    integer, dimension(3) :: gi1, gi2, gi3, gi4
    integer :: l, k
    ! Convert triplet to singlet
    gi1 = singlet_to_triplet(i1, dims(2), dims(3))
    gi2 = singlet_to_triplet(i2, dims(2), dims(3))
    gi3 = singlet_to_triplet(i3, dims(2), dims(3))
    do l = 1, 3
         gi4(l) = 4 - gi1(l) - gi2(l) - gi3(l)
   end do
   do k = 1, 3
   do l = 1, 3
       if (gi4(l) .lt. 1) gi4(l) = gi4(l) + dims(l)
       if (gi4(l) .gt. dims(l)) gi4(l) = gi4(l) - dims(l)
   end do
   end do

    ! convert it back to a singlet
    i4 = triplet_to_singlet(gi4, dims(2), dims(3))
end function
end module
