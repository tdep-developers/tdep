program pack_simulation
!!{!src/pack_simulation/manual.md!}
use konstanter, only: r8, lo_exitcode_io
use gottochblandat, only: open_file, lo_does_file_exist, tochar, lo_stop_gracefully, lo_mean, lo_stddev
use mpi_wrappers, only: lo_mpi_helper
use lo_memtracker, only: lo_mem_helper
use type_mdsim, only: lo_mdsim
use type_crystalstructure, only: lo_crystalstructure
use lo_symmetry_of_interactions, only: lo_interaction_tensors
use type_forcemap, only: lo_forcemap
use type_forceconstant_secondorder, only: lo_forceconstant_secondorder
use ifc_solvers, only: lo_solve_for_borncharges

use options, only: lo_opts

implicit none
type(lo_crystalstructure) :: uc, ss
type(lo_forceconstant_secondorder) :: fc
type(lo_interaction_tensors) :: slt
type(lo_forcemap) :: map
type(lo_opts) :: opts
type(lo_mdsim) :: sim, trimsim
type(lo_mpi_helper) :: mw
type(lo_mem_helper) :: mem

real(r8), dimension(3, 3) :: eps
real(r8), dimension(:, :, :), allocatable :: Z
integer :: u, i, j

! Get the command line arguments
call opts%parse()
call mw%init()
call mem%init()

! Sanity check to make sure no genius runs this in parallel
if (mw%n .gt. 1) then
    if (mw%talk) write (*, *) "Don't run this in parallel."
    call mw%destroy()
    stop
end if

write (*, *) 'PACKING SIMULATION TO HDF5'

! Read the simulation from file, first the perfect structures
call uc%readfromfile('infile.ucposcar')
call uc%classify('wedge', timereversal=.true.)
call ss%readfromfile('infile.ssposcar')
call ss%classify('supercell', uc)

write (*, *) '... read structures'

! Then some sanity checks
if (lo_does_file_exist('infile.projected_moments')) then
    if (opts%magnetic_moments) then
        write (*, *) '... found "infile.projected_moments"'
    else
        write (*, *) 'WARNING/OBSERVATION: found "infile.projected_moments", but you have not decided to pack it'
    end if
else
    if (opts%magnetic_moments) then
        call lo_stop_gracefully(['Could not find "infile.projected_moments"'], lo_exitcode_io)
    end if
end if
if (lo_does_file_exist('infile.lattice')) then
    if (opts%variable_lattice) then
        write (*, *) '... found "infile.lattice"'
    else
        write (*, *) 'WARNING/OBSERVATION: found "infile.lattice", but you have not decided to pack it'
    end if
else
    if (opts%variable_lattice) then
        call lo_stop_gracefully(['Could not find "infile.lattice"'], lo_exitcode_io)
    end if
end if
if (lo_does_file_exist('infile.born_charges')) then
    if (opts%dielectric) then
        write (*, *) '... found "infile.born_charges"'
    else
        write (*, *) 'WARNING/OBSERVATION: found "infile.born_charges", but you have not decided to pack it'
    end if
else
    if (opts%dielectric) then
        call lo_stop_gracefully(['Could not find "infile.born_charges"'], lo_exitcode_io)
    end if
end if
if (lo_does_file_exist('infile.dielectric_tensor')) then
    if (opts%dielectric) then
        write (*, *) '... found "infile.dielectric_tensor"'
    else
        write (*, *) 'WARNING/OBSERVATION: found "infile.dielectric_tensor", but you have not decided to pack it'
    end if
else
    if (opts%dielectric) then
        call lo_stop_gracefully(['Could not find "infile.dielectric_tensor"'], lo_exitcode_io)
    end if
end if

! Check if we are dealing with an alloy
if (uc%info%alloy) then
    fixalloy: block
        type(lo_mdsim) :: asim
        type(lo_crystalstructure) :: sqsid, sqsrl
        real(r8), dimension(:, :), allocatable :: buf_rref
        integer, dimension(:), allocatable :: alloy_permutation, buf_atomicnumbers
        integer :: i

        write (*, *) '... this is a random alloy, looking for more files'

        ! Hmmm. What do we need? The perfect SQS, and the relaxed SQS
        if (lo_does_file_exist('infile.sqs_ideal')) then
            write (*, *) '... found "infile.sqs_ideal"'
        else
            call lo_stop_gracefully(['Could not find "infile.sqs_ideal"'], lo_exitcode_io)
        end if

        if (lo_does_file_exist('infile.sqs_relaxed')) then
            write (*, *) '... found "infile.sqs_relaxed"'
        else
            call lo_stop_gracefully(['Could not find "infile.sqs_relaxed"'], lo_exitcode_io)
        end if

        ! We will definitely need the perfect and relaxed SQS
        call sqsid%readfromfile('infile.sqs_ideal', verbosity=opts%verbosity)
        call sqsrl%readfromfile('infile.sqs_relaxed', verbosity=opts%verbosity)

        ! Sort out the magic permutation thingy
        allocate (alloy_permutation(ss%na))
        alloy_permutation = 0
        call ss%alloy_site_permutation(sqsid, alloy_permutation)

        ! Read the simulation into a temporary simulation, since we have to massage it a little
        ! to sort out the alloy stuff.
        call asim%read_from_file(opts%verbosity + 2, magnetic=opts%magnetic_moments, &
                                 variable_lattice=opts%variable_lattice, &
                                 dynamics=opts%molecular_dynamics, &
                                 stride=opts%stride, nrand=opts%nrand, readalloy=.true.)

        write (*, *) '... starting to rearrange simulation'

        ! Generate an empty simulation to pack it into
        call sim%init_empty(uc, ss, asim%nt, asim%temperature_thermostat, &
                            timestep=asim%timestep, &
                            magnetic_moments=asim%have_magnetic_moments, &
                            variable_lattice=asim%variable_lattice)

        ! Permute the simulation we just read so that it aligns with infile.ssposcar in some sense.
        ! do i=1,ss%na
        !     write(*,*) ss%r(:,i)-sqsid%r(:,alloy_permutation(i))
        ! enddo

        ! Permute atomic numbers and relaxed positions
        allocate (buf_atomicnumbers(ss%na))
        allocate (buf_rref(3, ss%na))
        buf_rref = sqsrl%r(:, alloy_permutation)
        buf_atomicnumbers = sqsrl%atomic_number(alloy_permutation)

        ! Permute positions and forces
        do i = 1, asim%nt
            asim%r(:, :, i) = asim%r(:, alloy_permutation, i)
            asim%f(:, :, i) = asim%f(:, alloy_permutation, i)
        end do

        ! Store in clean simulation
        do i = 1, asim%nt
            call sim%add_timestep( &
                asim%r(:, :, i), &
                asim%f(:, :, i), &
                asim%stat%potential_energy(i), &
                asim%stat%kinetic_energy(i), &
                asim%stat%temperature(i), &
                asim%stat%stress(:, :, i), &
                atomic_numbers=buf_atomicnumbers, &
                reference_positions=buf_rref)
        end do
    end block fixalloy
else
    ! Not an alloy, just read it from file
    call sim%read_from_file(opts%verbosity + 2, &
                            magnetic=opts%magnetic_moments, &
                            variable_lattice=opts%variable_lattice, &
                            dynamics=opts%molecular_dynamics, &
                            dielectric=opts%dielectric, &
                            stride=opts%stride, &
                            nrand=opts%nrand)
end if

! Buffers for the Born charges and dielectric tensor
allocate (Z(3, 3, uc%na))
Z = 0.0_r8
eps = 0.0_r8

! Possibly clean the simulation a little
if (opts%make_tidy) then

    ! Remove drift in forces and positions.
    call sim%remove_force_and_center_of_mass_drift()

    if (lo_does_file_exist('infile.lotosplitting')) then
        u = open_file('in', 'infile.lotosplitting')
        do i = 1, 3
            read (u, *) eps(:, i)
        end do
        do j = 1, uc%na
        do i = 1, 3
            read (u, *) Z(:, i, j)
        end do
        end do
        close (u)

        ! Get all the symmetries
        call slt%generate(uc, ss, cutoff2=ss%mincutoff(), cutoff3=-1.0_r8, cutoff4=-1.0_r8, polar=.true., mw=mw, mem=mem, verbosity=-1)
        ! And create the map
        call map%generate(uc, ss, polarcorrectiontype=3, st=slt, mw=mw, mem=mem, verbosity=-1)
        if (map%xuc%nx_Z_singlet .gt. 0) then
            call lo_solve_for_borncharges(map, p=uc, Z=Z, eps=eps, verbosity=opts%verbosity, mw=mw, mem=mem)
            call map%get_secondorder_forceconstant(uc, fc, verbosity=opts%verbosity, mem=mem)
            eps = fc%loto%eps
            Z = fc%loto%born_effective_charges
        end if

        write (*, *) '... got symmetrized dielectric constant and Born effective charges'
    end if
end if

call sim%write_to_hdf5(uc, ss, 'outfile.sim.hdf5', opts%verbosity + 1, eps, Z)

! remove outliers that give really odd dielectric constant. This is a VASP
! artifact, sometimes VASP fails silently and gives dielectric constants that
! are completely bonkers and will screw everything up.
if (opts%nsigma .gt. 0) then
    removeoutliers: block
        real(r8) :: f0
        integer :: i, j, ctr

        trimsim = sim
        trimsim%nt = 0

        ctr = 0
        tsloop: do i = 1, sim%nt
            do j = 1, 3
                f0 = abs(sim%eps(j, j, i) - opts%eps0)
                if (f0 .gt. opts%nsigma*opts%dev0) cycle tsloop
            end do

            ctr = ctr + 1
            trimsim%r(:, :, ctr) = sim%r(:, :, i)
            trimsim%u(:, :, ctr) = sim%u(:, :, i)
            trimsim%f(:, :, ctr) = sim%f(:, :, i)
            trimsim%eps(:, :, ctr) = sim%eps(:, :, i)
            trimsim%Z(:, :, :, ctr) = sim%Z(:, :, :, i)
            trimsim%stat%potential_energy(ctr) = sim%stat%potential_energy(i)
            trimsim%stat%total_energy(ctr) = sim%stat%total_energy(i)
            trimsim%stat%kinetic_energy(ctr) = sim%stat%kinetic_energy(i)
            trimsim%stat%temperature(ctr) = sim%stat%temperature(i)
            trimsim%stat%pressure(ctr) = sim%stat%pressure(i)
            trimsim%stat%stress(:, :, ctr) = sim%stat%stress(:, :, i)
            write (*, *) i, ctr
        end do tsloop
        trimsim%nt = ctr
        call trimsim%write_to_hdf5(uc, ss, 'outfile.sim.hdf5', opts%verbosity + 1, eps, Z)
    end block removeoutliers
end if

write (*, *) ''
write (*, *) 'All done!'

end program
