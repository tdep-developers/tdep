#include "precompilerdefinitions"
program samples_from_md
!!{!src/samples_from_md/manual.md!}
use konstanter, only: flyt, lo_huge, lo_tol, lo_hugeint
use gottochblandat, only: open_file, lo_mean, lo_stddev, tochar, lo_progressbar_init, lo_progressbar, walltime, lo_stop_gracefully
use type_crystalstructure, only: lo_crystalstructure
use lo_randomnumbers, only: lo_mersennetwister
use type_mdsim, only: lo_mdsim
use options, only: lo_opts

implicit none
type(lo_opts) :: opts
type(lo_mdsim) :: sim
type(lo_crystalstructure) :: p
type(lo_mersennetwister) :: tw

integer, dimension(:), allocatable  :: sample, indices, testsample
integer :: i, j, k, ii, jj, mindist
real(flyt) :: f0, f1, f2, t0, t1
real(flyt) :: mean_ep, mean_ek, sigma_ep, sigma_ek
real(flyt) :: mp0, mp1, mk0, mk1, sp0, sp1, sk0, sk1
character(len=80) :: fname

! Get options
call opts%parse()
! Seed random numbers
call tw%init(iseed=0, rseed=walltime())
! Get simulation
call sim%read_from_hdf5('infile.sim.hdf5', verbosity=1)

! Stop when velocities are not present
if (.not. allocated(sim%v)) then
    call lo_stop_gracefully(['No velocities present in the simulation. Is this an MD?'], 9)
end if

! Check if the number of samples is not too large
if (opts%n .gt. sim%nt) then
    call lo_stop_gracefully(['Number of samples is larger than the number of timesteps in the simulation.'], 9)
end if

! List of all the samples
allocate (indices(sim%nt))
do i = 1, sim%nt
    indices(i) = i
end do

! Get a first sample
allocate (sample(opts%n))
allocate (testsample(opts%n))
sample = 0
testsample = 0
sample(1) = tw%rnd_int(sim%nt)

! Get the minimum distance between two samples
mindist = sim%nt/opts%n/5
write (*, *) 'minimum distance', mindist

! Choose the rest
do i = 2, opts%n
    j = select_ok_random_sample(indices, sample(1:i - 1), mindist, tw)
    sample(i) = j
end do

! Values for the entire simulation
mean_ep = lo_mean(sim%stat%potential_energy)
mean_ek = lo_mean(sim%stat%kinetic_energy)
sigma_ep = lo_stddev(sim%stat%potential_energy)
sigma_ek = lo_stddev(sim%stat%kinetic_energy)

! Starting values for monte-carlo run
mp0 = lo_mean(sim%stat%potential_energy(sample))
mk0 = lo_mean(sim%stat%kinetic_energy(sample))
sp0 = lo_stddev(sim%stat%potential_energy(sample))
sk0 = lo_stddev(sim%stat%kinetic_energy(sample))
mp0 = abs(mp0 - mean_ep)
mk0 = abs(mk0 - mean_ek)
sp0 = abs(sp0 - sigma_ep)
sk0 = abs(sk0 - sigma_ek)

t0 = mp0 + mk0 + sp0 + sk0

write (*, *) 'Starting Monte-Carlo, distances:'
write (*, "(9X,4(1X,A18))") '|<Ep>-<Ep0>|', '|<Ek>-<Ek0>|', '|S(Ep)-S(Ep0)|', '|S(Ek)-S(Ek0)|'
write (*, "(1X,I8,4(1X,F18.6))") 0, mp0, mk0, sp0, sk0

do i = 1, opts%maxiter
    ! try a new configuration
    testsample = sample
    ! CHANGE PLACES
    k = tw%rnd_int(opts%n)
    j = select_ok_random_sample(indices, sample, mindist, tw)
    testsample(k) = j

    ! Measure distances
    mp1 = lo_mean(sim%stat%potential_energy(testsample))
    mk1 = lo_mean(sim%stat%kinetic_energy(testsample))
    sp1 = lo_stddev(sim%stat%potential_energy(testsample))
    sk1 = lo_stddev(sim%stat%kinetic_energy(testsample))
    mp1 = abs(mp1 - mean_ep)
    mk1 = abs(mk1 - mean_ek)
    sp1 = abs(sp1 - sigma_ep)
    sk1 = abs(sk1 - sigma_ek)

    t1 = mp1 + mk1 + sp1 + sk1
    ! MC check
    f2 = tw%rnd_real()
    if (exp(-(t1 - t0)/opts%temp) .gt. f2) then
        sample = testsample
        mp0 = mp1
        mk0 = mk1
        sp0 = sp1
        sk0 = sk1
        t0 = t1
        write (*, "(1X,I8,4(1X,F18.6))") i, mp0, mk0, sp0, sk0
        if (t0 .lt. lo_tol) exit
    end if
end do

write (*, *) ' '
write (*, *) ' Results of Monte-Carlo run '
write (*, "('                 ',4(1X,A18))") '<Ep>', '<Ek>', 'Sigma(Ep)', 'Sigma(Ek)'
write (*, "('Full simulation:',4(1X,F18.6))") &
    lo_mean(sim%stat%potential_energy), lo_mean(sim%stat%kinetic_energy), lo_stddev(sim%stat%potential_energy), lo_stddev(sim%stat%kinetic_energy)
write (*, "('        Samples:',4(1X,F18.6))") &
    lo_mean(sim%stat%potential_energy(sample)), lo_mean(sim%stat%kinetic_energy(sample)), lo_stddev(sim%stat%potential_energy(sample)), lo_stddev(sim%stat%kinetic_energy(sample))

! Find the minimum, maximum and mean distance between two samples
ii = 0
jj = lo_hugeint
f0 = 0.0_flyt
do i = 1, opts%n
    f1 = lo_huge
    do j = 1, opts%n
        if (i .ne. j) then
            jj = min(abs((sample(i) - sample(j))), jj)
            f1 = min(abs((sample(i) - sample(j))*1.0_flyt), f1)
        end if
    end do
    ii = max(floor(f1), ii)
    f0 = f0 + f1
end do
f0 = f0/opts%n
write (*, *) '  Max distance:', ii
write (*, *) '  Min distance:', jj
write (*, *) ' Mean distance:', int(f0)

! write the samples
p = sim%crystalstructure
call tw%shuffle_int_array(sample)
do i = 1, opts%n
    p%r = sim%r(:, :, sample(i))
    p%v = sim%v(:, :, sample(i))
    p%info%title = 'sample '//tochar(sample(i), 7)
    fname = 'sample'//tochar(i, 5)
    select case (opts%output_format)
    case (1) ! VASP
        call p%writetofile(fname, opts%output_format, .true.)
    case (2) ! Abinit
        call p%writetofile(fname, opts%output_format, .true.)
    case (3) ! LAMMPS
        call lo_stop_gracefully(['Native LAMMPS IO was removed, please use external converters.'], 8)
    case (4) ! AIMS
        call p%writetofile(fname, opts%output_format, .true.)
    case default
        write (*, *) 'Output format not fixe yet'
        stop
    end select
end do

contains

!> choose a sample, but not too stupidly
integer function select_ok_random_sample(indices, forbidden, mindist, tw)
    integer, dimension(:), intent(in) :: indices, forbidden
    integer, intent(in) :: mindist
    type(lo_mersennetwister), intent(inout) :: tw
    !
    integer :: i, ii, j, l, nt, nf

    nf = size(forbidden, 1)
    nt = size(indices, 1)
    do ii = 1, 100
        i = tw%rnd_int(nt)
        l = 0
        do j = 1, nf
            if (abs(indices(i) - indices(forbidden(j))) .le. mindist) l = l + 1
        end do
        if (l .eq. 0) exit
    end do
    select_ok_random_sample = i
end function

end program
