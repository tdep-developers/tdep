#include "precompilerdefinitions"
module type_mdsim
!! Information about an MD simulation
use konstanter, only: r8, lo_pi, lo_huge, lo_hugeint, lo_sqtol, lo_status, lo_force_eVA_to_HartreeBohr, lo_eV_to_Hartree, &
                      lo_pressure_GPa_to_HartreeBohr, lo_A_to_bohr, lo_exitcode_param, lo_bohr_to_A, lo_Hartree_to_eV, &
                      lo_force_Hartreebohr_to_eVA, lo_pressure_HartreeBohr_to_GPa, lo_exitcode_io, lo_kb_Hartree, &
                      lo_time_fs_to_au, lo_time_au_to_fs, lo_velocity_au_to_Afs
use gottochblandat, only: open_file, lo_clean_fractional_coordinates, lo_mean, tochar, lo_trueNtimes, walltime, lo_sqnorm, &
                          lo_progressbar_init, lo_progressbar, lo_find_rotation_that_makes_strain_diagonal, lo_trace, &
                          lo_mass_from_Z
use mpi_wrappers, only: lo_mpi_helper, lo_stop_gracefully
use hdf5_wrappers, only: lo_hdf5_helper, lo_h5_store_attribute, lo_h5_store_data, lo_h5_read_attribute, &
                         lo_h5_read_data, lo_h5_does_dataset_exist
use lo_randomnumbers, only: lo_mersennetwister
use type_crystalstructure, only: lo_crystalstructure
use type_blas_lapack_wrappers, only: lo_gemm

implicit none
private
public :: lo_mdsim

!> All the energies in an simulation
type lo_mdsim_stat
    !> instantaneous temperature
    real(r8), dimension(:), allocatable :: temperature
    !> total energy
    real(r8), dimension(:), allocatable :: total_energy
    !> potential energy
    real(r8), dimension(:), allocatable :: potential_energy
    !> kinetic energy
    real(r8), dimension(:), allocatable :: kinetic_energy
    !> dipole-dipole energy
    real(r8), dimension(:), allocatable :: polar_potential_energy
    !> harmonic potential energy
    real(r8), dimension(:), allocatable :: secondorder_potential_energy
    !> third order potential energy
    real(r8), dimension(:), allocatable :: thirdorder_potential_energy
    !> fourth order potential energy
    real(r8), dimension(:), allocatable :: fourthorder_potential_energy
    !> magnetic potential energy
    real(r8), dimension(:), allocatable :: magnetic_potential_energy
    !> hydrostatic pressure
    real(r8), dimension(:), allocatable :: pressure
    !> stress tensor
    real(r8), dimension(:, :, :), allocatable :: stress
end type

!> information from the reference starting position, such as Born charges, unitcell and so on
type lo_mdsim_extra
    !> latticevectors
    real(r8), dimension(:, :), allocatable :: unitcell_latticevectors
    real(r8), dimension(:, :), allocatable :: supercell_latticevectors
    !> reference positions
    real(r8), dimension(:, :), allocatable :: unitcell_positions
    real(r8), dimension(:, :), allocatable :: supercell_positions
    !> atomic numbers
    integer, dimension(:), allocatable :: unitcell_atomic_numbers
    integer, dimension(:), allocatable :: supercell_atomic_numbers
    !> dielectric constant
    real(r8), dimension(:, :), allocatable :: dielectric_tensor
    real(r8), dimension(:, :, :), allocatable :: born_effective_charges
    !> Alloy things, if needed
    integer, dimension(:), allocatable ::      unitcell_componentcounter
    integer, dimension(:), allocatable ::      supercell_componentcounter
    integer, dimension(:, :), allocatable ::    unitcell_components
    integer, dimension(:, :), allocatable ::    supercell_components
    real(r8), dimension(:, :), allocatable :: unitcell_concentrations
    real(r8), dimension(:, :), allocatable :: supercell_concentrations
end type

!> A set of force-displacement configurations, typically an MD simulation
type lo_mdsim
    !> energies, stress temperatures
    type(lo_mdsim_stat) :: stat
    !> general information
    type(lo_mdsim_extra) :: extra
    !> reference lattice for simulation
    type(lo_crystalstructure) :: crystalstructure
    !> thermostat temperature
    real(r8) :: temperature_thermostat = -lo_huge
    !> is the shape fixed or not?
    logical :: variable_lattice = .false.
    !> do we have magnetic moments?
    logical :: have_magnetic_moments = .false.
    !> do we have born charges and dielectric tensors?
    logical :: have_dielectric = .false.
    !> Is it stochastic samples or molecular dynamics?
    logical :: stochastic = .false.
    !> Is this a simulation of a random alloy?
    logical :: alloy = .false.
    !> number of timesteps
    integer :: nt = -lo_hugeint
    !> number of atoms
    integer :: na = -lo_hugeint
    !> timestep in fs
    real(r8) :: timestep = -lo_huge
    !> positions in fractional coordinates
    real(r8), allocatable, dimension(:, :, :) :: r
    !> displacements from reference in cartesian coordinates
    real(r8), allocatable, dimension(:, :, :) :: u
    !> velocities in cartesian coordinates
    real(r8), allocatable, dimension(:, :, :) :: v
    !> forces in cartesian coordinates
    real(r8), allocatable, dimension(:, :, :) :: f
    !> positions in fractional coordinates with periodic boundary conditions removed.
    real(r8), allocatable, dimension(:, :, :) :: r_npbc
    !> projected magnetic moments on atoms
    real(r8), allocatable, dimension(:, :, :) :: m
    !> lattice vectors per timestep
    real(r8), allocatable, dimension(:, :, :) :: lattice
    !> atomic numbers, per timestep (or species? Not sure. This is unambiguous at least.)
    integer, allocatable, dimension(:, :) :: atomic_numbers
    !> reference position, per timestep
    real(r8), allocatable, dimension(:, :, :) :: r_ref
    !> born effective charge, per atom and timestep
    real(r8), allocatable, dimension(:, :, :, :) :: Z
    !> dielectric tensor, per timestep
    real(r8), allocatable, dimension(:, :, :) :: eps
contains
    !> read from file
    procedure :: read_from_file
    !> read from hdf5
    procedure :: read_from_hdf5
    !> get non-pbc positions
    procedure :: get_nonpbc_positions
    !> get velocities
    procedure :: get_velocities
    !> displacements
    procedure :: get_displacements
    !> remove drift
    procedure :: remove_force_and_center_of_mass_drift
    !> get the rms displacement
    procedure :: rms_displacement
    !> group timesteps
    procedure :: group_timesteps
    !> allocate and create empty simulation
    procedure :: init_empty
    !> add a timestep to a simulation
    procedure :: add_timestep
    !> write to hdf5
    procedure :: write_to_hdf5
    !> write to plain text
    procedure :: write_to_plaintext
end type

contains

!> Add a timestep to a simulation
subroutine add_timestep(sim, positions, forces, potential_energy, kinetic_energy, temperature, stresstensor, latticevectors, magnetic_moments, velocities, atomic_numbers, reference_positions)
    !> md simulation
    class(lo_mdsim), intent(inout) :: sim
    !> positions, in fractional coordinates
    real(r8), dimension(:, :), intent(in) :: positions
    !> forces, in Cartesian coordinates
    real(r8), dimension(:, :), intent(in) :: forces
    !> potential energy (per cell)
    real(r8), intent(in) :: potential_energy
    !> kinetic energy (per cell)
    real(r8), intent(in) :: kinetic_energy
    !> temperature
    real(r8), intent(in) :: temperature
    !> stress tensor
    real(r8), dimension(3, 3), intent(in) :: stresstensor
    !> lattice vectors
    real(r8), dimension(3, 3), intent(in), optional :: latticevectors
    !> magnetic moment per atom
    real(r8), dimension(:, :), intent(in), optional :: magnetic_moments
    !> velocities
    real(r8), dimension(:, :), intent(in), optional :: velocities
    !> atomic numbers?
    integer, dimension(:), intent(in), optional :: atomic_numbers
    !> reference positions
    real(r8), dimension(:, :), intent(in), optional :: reference_positions

    integer :: i, t, tmax

    ! Current step
    t = sim%nt
    ! Current max number of timesteps
    tmax = size(sim%r, 3)
    ! Sanity tests
    if (t + 1 .gt. tmax) then
        call lo_stop_gracefully(['Not enough space to store timestep'], lo_exitcode_param, __FILE__, __LINE__)
        ! In the future I should grow the arrays instead.
    end if
    if (sim%have_magnetic_moments) then
        if (.not. present(magnetic_moments)) then
            call lo_stop_gracefully(['Have to provide magnetic moments'], lo_exitcode_param, __FILE__, __LINE__)
        end if
    end if
    if (sim%variable_lattice) then
        if (.not. present(latticevectors)) then
            call lo_stop_gracefully(['Have to provide latticevectors'], lo_exitcode_param, __FILE__, __LINE__)
        end if
    end if
    if (sim%stochastic .eqv. .false.) then
        if (.not. present(velocities)) then
            call lo_stop_gracefully(['Have to provide velocities'], lo_exitcode_param, __FILE__, __LINE__)
        end if
    end if
    if (sim%alloy) then
        if (.not. present(atomic_numbers)) then
            call lo_stop_gracefully(['Have to provide atomic numbers'], lo_exitcode_param, __FILE__, __LINE__)
        end if
        if (.not. present(reference_positions)) then
            call lo_stop_gracefully(['Have to provide reference positions'], lo_exitcode_param, __FILE__, __LINE__)
        end if
    end if

    ! Increment the counter for number of steps
    t = t + 1

    ! Start storing things
    sim%r(:, :, t) = positions
    sim%f(:, :, t) = forces
    if (present(latticevectors) .and. sim%variable_lattice) then
        sim%lattice(:, :, t) = latticevectors
    end if
    if (present(magnetic_moments) .and. sim%have_magnetic_moments) then
        sim%m(:, :, t) = magnetic_moments
    end if
    if (present(velocities)) then
    if (sim%stochastic .eqv. .false.) then
        sim%v(:, :, t) = velocities
    end if
    end if
    if (sim%alloy) then
        if (present(atomic_numbers)) then
            sim%atomic_numbers(:, t) = atomic_numbers
        end if
        if (present(reference_positions)) then
            sim%r_ref(:, :, t) = reference_positions
        end if
    end if

    ! Fix displacements
    do i = 1, sim%na
        if (sim%alloy) then
            sim%u(:, i, t) = sim%r(:, i, t) - sim%r_ref(:, i, t)
        else
            sim%u(:, i, t) = sim%r(:, i, t) - sim%crystalstructure%r(:, i)
        end if
        sim%u(:, i, t) = sim%crystalstructure%displacement_fractional_to_cartesian(sim%u(:, i, t))
    end do

    ! Store energies and stuff
    sim%stat%stress(:, :, t) = stresstensor
    sim%stat%potential_energy(t) = potential_energy
    sim%stat%kinetic_energy(t) = kinetic_energy
    sim%stat%total_energy(t) = potential_energy + kinetic_energy
    sim%stat%temperature(t) = temperature

    ! Make a note that we have added another step
    sim%nt = sim%nt + 1
end subroutine

!> create empty sim container, to be filled incrementally
subroutine init_empty(sim, uc, ss, nstep, temperature, timestep, magnetic_moments, variable_lattice)
    !> md simulation
    class(lo_mdsim), intent(out) :: sim
    !> unitcell
    type(lo_crystalstructure), intent(in) :: uc
    !> supercell
    type(lo_crystalstructure), intent(in) :: ss
    !> how many steps
    integer, intent(in) :: nstep
    !> specify temperature
    real(r8), intent(in) :: temperature
    !> specify timestep
    real(r8), intent(in) :: timestep
    !> make space for magnetic moments?
    logical, intent(in) :: magnetic_moments
    !> make space for lattice vectors?
    logical, intent(in) :: variable_lattice

    ! Set metadata
    init: block
        ! Some sanity tests
        if (ss%info%supercell .eqv. .false.) then
            call lo_stop_gracefully(['Need a supercell/unitcell pair to generate empty simulation'], lo_exitcode_param, __FILE__, __LINE__)
        end if
        if (nstep .lt. 1) then
            call lo_stop_gracefully(['Need at least one timestep to generate empty simulation'], lo_exitcode_param, __FILE__, __LINE__)
        end if

        ! Set some basic things
        sim%crystalstructure = ss
        sim%na = ss%na
        sim%nt = 0
        if (timestep .gt. 0.0_r8) then
            sim%timestep = timestep
            sim%stochastic = .false.
        else
            sim%timestep = -1.0_r8
            sim%stochastic = .true.
        end if
        sim%temperature_thermostat = temperature
        sim%have_magnetic_moments = magnetic_moments
        sim%variable_lattice = variable_lattice
        ! attempt at dealing with alloy things
        if (uc%info%alloy) then
            sim%alloy = .true.
        else
            sim%alloy = .false.
        end if
    end block init

    ! Info on the perfect structure, will come handy later
    structure: block
        integer :: i, j, k, l

        lo_allocate(sim%extra%unitcell_positions(3, uc%na))
        lo_allocate(sim%extra%supercell_positions(3, ss%na))
        lo_allocate(sim%extra%unitcell_atomic_numbers(uc%na))
        lo_allocate(sim%extra%supercell_atomic_numbers(ss%na))
        lo_allocate(sim%extra%unitcell_latticevectors(3, 3))
        lo_allocate(sim%extra%supercell_latticevectors(3, 3))
        lo_allocate(sim%extra%dielectric_tensor(3, 3))
        lo_allocate(sim%extra%born_effective_charges(3, 3, uc%na))
        sim%extra%unitcell_positions = uc%r
        sim%extra%supercell_positions = ss%r
        sim%extra%unitcell_atomic_numbers = uc%atomic_number
        sim%extra%supercell_atomic_numbers = ss%atomic_number
        sim%extra%unitcell_latticevectors = uc%latticevectors
        sim%extra%supercell_latticevectors = ss%latticevectors
        sim%extra%dielectric_tensor = 0.0_r8
        sim%extra%born_effective_charges = 0.0_r8
        ! Also, if alloy store a lot of extra things
        if (sim%alloy) then
            ! Figure out the max number of components
            l = 0
            do i = 1, uc%na
                l = max(l, uc%alloyspecies(uc%species(i))%n)
            end do
            ! Space for alloy specification
            lo_allocate(sim%extra%unitcell_componentcounter(uc%na))
            lo_allocate(sim%extra%supercell_componentcounter(ss%na))
            lo_allocate(sim%extra%unitcell_components(l, uc%na))
            lo_allocate(sim%extra%supercell_components(l, ss%na))
            lo_allocate(sim%extra%unitcell_concentrations(l, uc%na))
            lo_allocate(sim%extra%supercell_concentrations(l, ss%na))
            sim%extra%unitcell_componentcounter = -1
            sim%extra%supercell_componentcounter = -1
            sim%extra%unitcell_components = -1
            sim%extra%supercell_components = -1
            sim%extra%unitcell_concentrations = 0.0_r8
            sim%extra%supercell_concentrations = 0.0_r8
            ! Store alloy specification
            do i = 1, uc%na
                j = uc%species(i)
                sim%extra%unitcell_componentcounter(i) = uc%alloyspecies(j)%n
                do k = 1, sim%extra%unitcell_componentcounter(i)
                    sim%extra%unitcell_components(k, i) = uc%alloyspecies(j)%atomic_number(k)
                    sim%extra%unitcell_concentrations(k, i) = uc%alloyspecies(j)%concentration(k)
                end do
            end do
            do i = 1, ss%na
                j = ss%species(i)
                sim%extra%supercell_componentcounter(i) = ss%alloyspecies(j)%n
                do k = 1, sim%extra%supercell_componentcounter(i)
                    sim%extra%supercell_components(k, i) = ss%alloyspecies(j)%atomic_number(k)
                    sim%extra%supercell_concentrations(k, i) = ss%alloyspecies(j)%concentration(k)
                end do
            end do
        else
            ! Allocate dummy empty arrays
            lo_allocate(sim%extra%unitcell_componentcounter(1))
            lo_allocate(sim%extra%supercell_componentcounter(1))
            lo_allocate(sim%extra%unitcell_components(1, 1))
            lo_allocate(sim%extra%supercell_components(1, 1))
            lo_allocate(sim%extra%unitcell_concentrations(1, 1))
            lo_allocate(sim%extra%supercell_concentrations(1, 1))
            sim%extra%unitcell_componentcounter = -lo_hugeint
            sim%extra%supercell_componentcounter = -lo_hugeint
            sim%extra%unitcell_components = -lo_hugeint
            sim%extra%supercell_components = -lo_hugeint
            sim%extra%unitcell_concentrations = -lo_huge
            sim%extra%supercell_concentrations = -lo_huge
        end if
    end block structure

    trajectories: block
        ! Space for trajectories, some are always there
        lo_allocate(sim%r(3, sim%na, nstep))
        lo_allocate(sim%u(3, sim%na, nstep))
        lo_allocate(sim%f(3, sim%na, nstep))
        sim%r = 0.0_r8
        sim%u = 0.0_r8
        sim%f = 0.0_r8

        if (sim%stochastic .eqv. .false.) then
            lo_allocate(sim%v(3, sim%na, nstep))
            lo_allocate(sim%r_npbc(3, sim%na, nstep))
            sim%v = 0.0_r8
            sim%r_npbc = 0.0_r8
        end if
        if (sim%have_magnetic_moments) then
            lo_allocate(sim%m(3, sim%na, nstep))
            sim%m = 0.0_r8
        end if
        if (sim%variable_lattice) then
            lo_allocate(sim%lattice(3, 3, nstep))
            sim%lattice = 0.0_r8
        end if
        if (sim%alloy) then
            lo_allocate(sim%r_ref(3, sim%na, nstep))
            lo_allocate(sim%atomic_numbers(sim%na, nstep))
            sim%r_ref = 0.0_r8
            sim%atomic_numbers = 0
        end if
    end block trajectories

    ! And some space for energies
    energies: block
        lo_allocate(sim%stat%potential_energy(nstep))
        lo_allocate(sim%stat%kinetic_energy(nstep))
        lo_allocate(sim%stat%total_energy(nstep))
        lo_allocate(sim%stat%temperature(nstep))
        lo_allocate(sim%stat%pressure(nstep))
        lo_allocate(sim%stat%magnetic_potential_energy(nstep))
        lo_allocate(sim%stat%polar_potential_energy(nstep))
        lo_allocate(sim%stat%secondorder_potential_energy(nstep))
        lo_allocate(sim%stat%thirdorder_potential_energy(nstep))
        lo_allocate(sim%stat%fourthorder_potential_energy(nstep))
        lo_allocate(sim%stat%stress(3, 3, nstep))
        sim%stat%potential_energy = 0.0_r8
        sim%stat%kinetic_energy = 0.0_r8
        sim%stat%total_energy = 0.0_r8
        sim%stat%temperature = 0.0_r8
        sim%stat%pressure = 0.0_r8
        sim%stat%magnetic_potential_energy = 0.0_r8
        sim%stat%polar_potential_energy = 0.0_r8
        sim%stat%secondorder_potential_energy = 0.0_r8
        sim%stat%thirdorder_potential_energy = 0.0_r8
        sim%stat%fourthorder_potential_energy = 0.0_r8
        sim%stat%stress = 0.0_r8
    end block energies
end subroutine

!> put time steps in groups of constant position
subroutine group_timesteps(sim, groupcounter, groupind, ngroup, nconf)
    !> md simulation
    class(lo_mdsim), intent(in) :: sim
    !> number of steps in each group
    integer, dimension(:), allocatable, intent(out) :: groupcounter
    !> indices to the steps in each group
    integer, dimension(:, :), allocatable, intent(out) :: groupind
    !> number of groups
    integer, intent(out) :: ngroup
    !> total number of differences
    integer, intent(out) :: nconf

    integer, parameter :: grthres = 2 ! how many structures counts as a group
    real(r8) :: f0
    integer :: i, j, ii
    integer, dimension(:, :), allocatable :: ddi
    integer, dimension(:), allocatable :: ctr, di

    lo_allocate(ctr(sim%nt))
    lo_allocate(di(sim%nt))
    lo_allocate(ddi(sim%nt, sim%nt))
    ctr = 0
    di = 1
    ddi = 0
    ! Figure out how many unique positions there are?
    do i = 1, sim%nt
        if (di(i) .eq. 0) cycle
        do j = i, sim%nt
            ! really simple check
            f0 = abs(sim%r(1, 1, j) - sim%r(1, 1, i))
            if (f0 .gt. lo_sqtol) cycle ! not the same
            ! slightly better check
            f0 = sum(abs(sim%r(:, :, j) - sim%r(:, :, i)))/sim%na/3
            if (f0 .gt. lo_sqtol) cycle ! not the same
            ! if I made it here, it's equal
            di(j) = 0
            ctr(i) = ctr(i) + 1
            ddi(ctr(i), i) = j
        end do
    end do

    ! Build the actual groups
    ngroup = 0
    do i = 1, sim%nt
        if (ctr(i) .ge. grthres) ngroup = ngroup + 1
    end do
    lo_allocate(groupcounter(ngroup))
    lo_allocate(groupind(maxval(ctr), ngroup))
    groupcounter = 0
    groupind = 0
    j = 0
    do i = 1, sim%nt
        if (ctr(i) .ge. grthres) then
            j = j + 1
            groupcounter(j) = ctr(i)
            groupind(1:groupcounter(j), j) = ddi(1:ctr(i), i)
        end if
    end do

    ! Count number of differences
    nconf = 0
    do ii = 1, ngroup
        do i = 1, groupcounter(ii)
        do j = i + 1, groupcounter(ii)
            nconf = nconf + 1
        end do
        end do
    end do
end subroutine

!> Remove periodic boundary conditions
subroutine get_nonpbc_positions(sim)
    !> the simulation
    class(lo_mdsim), intent(inout) :: sim

    real(r8), dimension(3) :: v0
    integer :: i, j

    if (sim%stochastic) then
        call lo_stop_gracefully(['Trying to get non-pbc positions for stochastic simulation. Makes no sense.'], lo_exitcode_param, __FILE__, __LINE__)
    end if

    ! Old version
    sim%r_npbc = 0.0_r8
    ! Fix the first timestep
    do j = 1, sim%na
        if (sim%alloy) then
            v0 = sim%r(:, j, 1) - sim%r_ref(:, j, 1)
        else
            v0 = sim%r(:, j, 1) - sim%crystalstructure%r(:, j)
        end if
        ! Wrap it around
        v0 = lo_clean_fractional_coordinates(v0 + 0.5_r8) - 0.5_r8
        ! Store the unwrapped position
        sim%r_npbc(:, j, 1) = sim%crystalstructure%r(:, j) + v0
    end do

    ! Now align all the other timesteps with respect to this
    do i = 2, sim%nt
    do j = 1, sim%na
        ! Shift from previous timestep
        v0 = sim%r(:, j, i) - sim%r_npbc(:, j, i - 1)
        ! Wrap it around
        v0 = lo_clean_fractional_coordinates(v0 + 0.5_r8) - 0.5_r8
        ! Store the unwrapped position
        sim%r_npbc(:, j, i) = sim%r_npbc(:, j, i - 1) + v0
    end do
    end do
end subroutine

!> Calculate the velocities via backwards velocity verlet
subroutine get_velocities(sim)
    !> the simulation
    class(lo_mdsim), intent(inout) :: sim

    real(r8), dimension(3, sim%na) :: dr
    real(r8) :: its
    integer :: i

    ! Can only use this with MD
    if (sim%stochastic) then
        call lo_stop_gracefully(['Trying to get velocities for stochastic simulation. Makes no sense.'], lo_exitcode_param, __FILE__, __LINE__)
    end if
    if (sim%variable_lattice) then
        call lo_stop_gracefully(['Trying to get velocities for npt simulation. Have to think.'], lo_exitcode_param, __FILE__, __LINE__)
    end if
    sim%v = 0.0_r8
    ! If only one timestep don't bother
    if (sim%nt .eq. 1) return

    its = 1.0_r8/(sim%timestep*2.0_r8)
    do i = 2, sim%nt - 1
        dr = sim%r_npbc(:, :, i + 1) - sim%r_npbc(:, :, i - 1)
        dr = dr*its
        call lo_gemm(sim%crystalstructure%latticevectors, dr, sim%v(:, :, i))
    end do
    ! endpoints
    dr = (sim%r_npbc(:, :, 2) - sim%r_npbc(:, :, 1))*2.0_r8*its
    call lo_gemm(sim%crystalstructure%latticevectors, dr, sim%v(:, :, 1))
    dr = (sim%r_npbc(:, :, sim%nt - 1) - sim%r_npbc(:, :, sim%nt))*2.0_r8*its
    call lo_gemm(sim%crystalstructure%latticevectors, dr, sim%v(:, :, sim%nt))

    ! ! get velocities in fractional coordinates
    ! do i=2,sim%nt-1
    !     do j=1,sim%na
    !         do k=1,3
    !             sim%v(k,j,i)=sim%r_npbc(k,j,i+1)-sim%r_npbc(k,j,i-1)
    !         enddo
    !     enddo
    ! enddo
    ! ! endpoints
    ! do j=1,sim%na
    !     do k=1,3
    !         sim%v(k,j,1)=(sim%r_npbc(k,j,2)-sim%r_npbc(k,j,1))*2.0_r8
    !         sim%v(k,j,sim%nt)=(sim%r_npbc(k,j,sim%nt-1)-sim%r_npbc(k,j,sim%nt))*2.0_r8
    !     enddo
    ! enddo
    ! ! convert to cartesian
    ! do i=1,sim%nt
    !     do j=1,sim%na
    !         sim%v(:,j,i)=sim%crystalstructure%fractional_to_cartesian(sim%v(:,j,i))/(sim%timestep*2.0_r8)
    !     enddo
    ! enddo
end subroutine

!> Calculate displacements from equilibrium
subroutine get_displacements(sim)
    !> the simulation
    class(lo_mdsim), intent(inout) :: sim

    real(r8), dimension(3, sim%na) :: dr
    integer :: i !,j

    ! do i=1,sim%nt
    !     do j=1,sim%na
    !         sim%u(:,j,i)=sim%r(:,j,i)-sim%crystalstructure%r(:,j)
    !         sim%u(:,j,i)=sim%crystalstructure%displacement_fractional_to_cartesian(sim%u(:,j,i))
    !     enddo
    ! enddo

    if (sim%stochastic) then
        if (sim%alloy) then
            do i = 1, sim%nt
                dr = lo_clean_fractional_coordinates(sim%r(:, :, i) - sim%r_ref(:, :, i) + 0.5_r8) - 0.5_r8
                call lo_gemm(sim%crystalstructure%latticevectors, dr, sim%u(:, :, i))
            end do
        else
            do i = 1, sim%nt
                dr = lo_clean_fractional_coordinates(sim%r(:, :, i) - sim%crystalstructure%r + 0.5_r8) - 0.5_r8
                call lo_gemm(sim%crystalstructure%latticevectors, dr, sim%u(:, :, i))
            end do
        end if
    else
        do i = 1, sim%nt
            dr = sim%r_npbc(:, :, i) - sim%crystalstructure%r
            call lo_gemm(sim%crystalstructure%latticevectors, dr, sim%u(:, :, i))
        end do
    end if
end subroutine

!> Adjust so that forces sum up to 0, and that the center of mass does not move.
subroutine remove_force_and_center_of_mass_drift(sim)
    !> simulation
    class(lo_mdsim), intent(inout) :: sim

    real(r8), dimension(3) :: ftot, ctot
    real(r8) :: invmass, invna
    integer :: i, j

    invmass = 1.0_r8/sum(sim%crystalstructure%mass)
    invna = 1.0_r8/real(sim%na, r8)

    ! Done slightly differently depending on wether it is stochastic or not.
    if (sim%stochastic) then

        if (sim%alloy) then
            do i = 1, sim%nt
                invmass = 1.0_r8/sum(lo_mass_from_Z(sim%atomic_numbers(:, i)))

                ftot = 0.0_r8
                ctot = 0.0_r8
                do j = 1, sim%na
                    ftot = ftot + sim%f(:, j, i)
                    ctot = ctot + (lo_clean_fractional_coordinates(sim%r(:, j, i) - sim%r_ref(:, j, i) + 0.5_r8) - 0.5_r8)*lo_mass_from_Z(sim%atomic_numbers(j, i))
                end do
                ftot = ftot*invna
                ctot = ctot*invmass
                ! Remove the drift
                do j = 1, sim%na
                    sim%f(:, j, i) = sim%f(:, j, i) - ftot
                    sim%r(:, j, i) = sim%r(:, j, i) - ctot
                    sim%u(:, j, i) = sim%r(:, j, i) - sim%r_ref(:, j, i)
                    sim%u(:, j, i) = sim%crystalstructure%displacement_fractional_to_cartesian(sim%u(:, j, i))
                end do
            end do
        else
            do i = 1, sim%nt
                ! Calculate the drift
                ftot = 0.0_r8
                ctot = 0.0_r8
                do j = 1, sim%na
                    ftot = ftot + sim%f(:, j, i)
                    ctot = ctot + (lo_clean_fractional_coordinates(sim%r(:, j, i) - sim%crystalstructure%r(:, j) + 0.5_r8) - 0.5_r8)*sim%crystalstructure%mass(j)
                end do
                ftot = ftot*invna
                ctot = ctot*invmass
                ! Remove the drift
                do j = 1, sim%na
                    sim%f(:, j, i) = sim%f(:, j, i) - ftot
                    sim%r(:, j, i) = sim%r(:, j, i) - ctot
                    sim%u(:, j, i) = sim%r(:, j, i) - sim%crystalstructure%r(:, j)
                    sim%u(:, j, i) = sim%crystalstructure%displacement_fractional_to_cartesian(sim%u(:, j, i))
                end do
            end do
        end if
        ! Put the coordinates back in the 0-1 box they belong.
        sim%r = lo_clean_fractional_coordinates(sim%r)
    else

        if (sim%alloy) then
            do i = 1, sim%nt
                invmass = 1.0_r8/sum(lo_mass_from_Z(sim%atomic_numbers(:, i)))
                ! Calculate the drift
                ftot = 0.0_r8
                ctot = 0.0_r8
                do j = 1, sim%na
                    ftot = ftot + sim%f(:, j, i)
                    ctot = ctot + (sim%r_npbc(:, j, i) - sim%r_ref(:, j, i))*lo_mass_from_Z(sim%atomic_numbers(j, i))
                    ctot = ctot + (sim%r_npbc(:, j, i) - sim%crystalstructure%r(:, j))*sim%crystalstructure%mass(j)
                end do
                ftot = ftot*invna
                ctot = ctot*invmass
                ! Remove the drift
                do j = 1, sim%na
                    sim%f(:, j, i) = sim%f(:, j, i) - ftot
                    sim%r(:, j, i) = sim%r(:, j, i) - ctot
                    sim%r_npbc(:, j, i) = sim%r_npbc(:, j, i) - ctot
                    sim%u(:, j, i) = sim%r(:, j, i) - sim%r_ref(:, j, i)
                    sim%u(:, j, i) = sim%crystalstructure%displacement_fractional_to_cartesian(sim%u(:, j, i))
                end do
            end do
            ! Put the coordinates back in the 0-1 box they belong.
            sim%r = lo_clean_fractional_coordinates(sim%r)
        else
            do i = 1, sim%nt
                ! Calculate the drift
                ftot = 0.0_r8
                ctot = 0.0_r8
                do j = 1, sim%na
                    ftot = ftot + sim%f(:, j, i)
                    ctot = ctot + (sim%r_npbc(:, j, i) - sim%crystalstructure%r(:, j))*sim%crystalstructure%mass(j)
                end do
                ftot = ftot*invna
                ctot = ctot*invmass
                ! Remove the drift
                do j = 1, sim%na
                    sim%f(:, j, i) = sim%f(:, j, i) - ftot
                    sim%r(:, j, i) = sim%r(:, j, i) - ctot
                    sim%r_npbc(:, j, i) = sim%r_npbc(:, j, i) - ctot
                    sim%u(:, j, i) = sim%r(:, j, i) - sim%crystalstructure%r(:, j)
                    sim%u(:, j, i) = sim%crystalstructure%displacement_fractional_to_cartesian(sim%u(:, j, i))
                end do
            end do
            ! Put the coordinates back in the 0-1 box they belong.
            sim%r = lo_clean_fractional_coordinates(sim%r)
        end if

    end if
end subroutine

!> Mean square displacement per unitcell atom?
subroutine rms_displacement(sim, uc, ss)
    !> simulation
    class(lo_mdsim), intent(in) :: sim
    !> unitcell
    type(lo_crystalstructure), intent(in) :: uc
    !> supercell
    type(lo_crystalstructure), intent(in) :: ss

    real(r8), dimension(:), allocatable :: msd
    integer, dimension(:), allocatable :: ctr
    integer :: i, j, t

    if (ss%info%supercell .eqv. .false.) then
        call lo_stop_gracefully(['Unitcell-supercell mapping needs be be resolved for rms displacement'], &
                                lo_exitcode_param, __FILE__, __LINE__)
    end if

    lo_allocate(msd(uc%na))
    lo_allocate(ctr(uc%na))
    msd = 0.0_r8
    ctr = 0.0_r8
    do t = 1, sim%nt
        do i = 1, sim%na
            j = ss%info%index_in_unitcell(i)
            msd(j) = msd(j) + lo_sqnorm(sim%u(:, i, t))
            ctr(j) = ctr(j) + 1
        end do
    end do
    msd = msd/real(ctr, r8)
end subroutine

!> Reads a simulation from file
subroutine read_from_file(sim, verbosity, stride, magnetic, dielectric, variable_lattice, dynamics, nrand, mw, readalloy)
    !> the md simulation
    class(lo_mdsim), intent(out) :: sim
    !> how much to talk
    integer, intent(in) :: verbosity
    !> read it in strides
    integer, intent(in), optional :: stride
    !> read magnetic moments
    logical, intent(in), optional :: magnetic
    !> read dielectric things
    logical, intent(in), optional :: dielectric
    !> is it a constant pressure simulation
    logical, intent(in), optional :: variable_lattice
    !> is it molecular dynamics or stochastic?
    logical, intent(in), optional :: dynamics
    !> read only N configutations, at random
    integer, intent(in), optional :: nrand
    !> mpi helper
    type(lo_mpi_helper), intent(inout), optional :: mw
    !> are we reading raw data for an alloy, to be processed later
    logical, intent(in), optional :: readalloy

    type(lo_mersennetwister) :: tw
    real(r8), dimension(:, :, :), allocatable :: idealpos, refpos
    real(r8) :: timer
    integer, dimension(:, :), allocatable :: atomic_numbers
    integer :: step, nstep_tot, nrnd, readrank
    logical, dimension(:), allocatable :: readstep
    logical :: readonthisrank, mpiparallel

    timer = walltime()

    ! Check that all files are there, allocate space and so on
    init: block
        integer, dimension(:), allocatable :: ind
        integer :: i, u

        ! talk?
        if (verbosity .gt. 0) then
            write (*, *) ''
            write (*, *) 'READING SIMULATION FROM FILE'
        end if

        ! Which rank will do all the reading?
        if (present(mw)) then
            ! We are reading on one rank, then broadcasting
            readrank = 0
            if (mw%r .eq. readrank) then
                readonthisrank = .true.
            else
                readonthisrank = .false.
            end if
            mpiparallel = .true.
        else
            ! Reading not in parallel.
            readonthisrank = .true.
            mpiparallel = .false.
        end if

        ! read in strides?
        if (present(stride)) then
            step = stride
        else
            step = 1
        end if

        ! read only a few at random
        if (present(nrand)) then
            nrnd = nrand
            call tw%init(iseed=0, rseed=walltime())
        else
            nrnd = -1
        end if

        ! read latticevectors as a function of time?
        if (present(variable_lattice)) then
            sim%variable_lattice = variable_lattice
        else
            sim%variable_lattice = .false.
        end if
        ! get the ideal supercell
        call sim%crystalstructure%readfromfile('infile.ssposcar')

        ! Is it an annoying alloy? If so, we need to have made an active choice, I think.
        if (sim%crystalstructure%info%alloy) then
            if (.not. present(readalloy)) then
                call lo_stop_gracefully(['Simulations of alloys need to be processed with "pack_simulation", see the manual.'], lo_exitcode_io, communicator=mw%comm)
            else
                if (readalloy .eqv. .false.) then
                    call lo_stop_gracefully(['Conflicting instructions wether it is an alloy or not.'], lo_exitcode_io, communicator=mw%comm)
                end if
            end if
        end if
        ! There could be other weird things happening.
        if (present(readalloy)) then
        if (readalloy) then
            if (sim%crystalstructure%info%alloy .eqv. .false.) then
                call lo_stop_gracefully(['Conflicting instructions wether it is an alloy or not.'], lo_exitcode_io, communicator=mw%comm)
            end if
        end if
        end if

        ! Read annoying magnetic moments
        if (present(magnetic)) then
            sim%have_magnetic_moments = magnetic
        else
            sim%have_magnetic_moments = .false.
        end if

        ! Read many born charges and dielectric constants
        if (present(dielectric)) then
            sim%have_dielectric = dielectric
        else
            sim%have_dielectric = .false.
        end if

        ! read the meta file to get some more reasonable things
        if (readonthisrank) then
            u = open_file('in', 'infile.meta')
            read (u, *) sim%na
            read (u, *) nstep_tot
            read (u, *) sim%timestep
            read (u, *) sim%temperature_thermostat
            close (u)
        end if

        ! Tell all the other ranks
        if (mpiparallel) then
            call mw%bcast(sim%na, readrank, __FILE__, __LINE__)
            call mw%bcast(nstep_tot, readrank, __FILE__, __LINE__)
            call mw%bcast(sim%timestep, readrank, __FILE__, __LINE__)
            call mw%bcast(sim%temperature_thermostat, readrank, __FILE__, __LINE__)
        end if

        ! Probably specify if it's dynamics or not manually.
        if (present(dynamics)) then
            sim%stochastic = .not. dynamics
        else
            ! If not explicitly specified, I assume it is stochastic. Think that makes sense.
            sim%stochastic = .true.
        end if

        ! Sanity check that the timestep makes sense.
        if (.not. sim%stochastic) then
        if (sim%timestep .lt. 0.0_r8) then
            call lo_stop_gracefully(['Strange timestep in "infile.meta" for molecular dynamics'], lo_exitcode_io, communicator=mw%comm)
        end if
        end if

        ! sanity checks, perhaps add some more.
        if (sim%na .ne. sim%crystalstructure%na) then
            call lo_stop_gracefully(['Not the same number of atoms in "infile.meta" and "infile.ssposcar"'], lo_exitcode_io, communicator=mw%comm)
        end if
        if (nstep_tot .lt. 1) then
            call lo_stop_gracefully(['Need at least one timestep.'], lo_exitcode_io, communicator=mw%comm)
        end if
        if (nrnd .gt. 0 .and. step .gt. 1) then
            call lo_stop_gracefully(['Makes no sense to choose random subset and stride at the same time.'], lo_exitcode_param, communicator=mw%comm)
        end if
        if (nrnd .gt. 0) then
        if (sim%stochastic .eqv. .false.) then
            call lo_stop_gracefully(['Makes no sense to choose random subset and molecular dynamics at the same time.'], lo_exitcode_param, communicator=mw%comm)
        end if
        end if

        ! make sure the subset things make sense
        if (nrnd .gt. 0) nrnd = min(nrnd, nstep_tot)
        if (step .gt. nstep_tot) step = nstep_tot

        ! Decide which timesteps to read. Most probably all of them, but you never know.
        lo_allocate(readstep(nstep_tot))
        readstep = .false.
        if (nrnd .gt. 0) then
            allocate (ind(nstep_tot))
            do i = 1, nstep_tot
                ind(i) = i
            end do
            ! shuffle timesteps
            call tw%shuffle_int_array(ind)
            if (mpiparallel) then
                call mw%bcast(ind, readrank, __FILE__, __LINE__)
            end if
            ! decide which ones to read
            do i = 1, nrnd
                readstep(ind(i)) = .true.
            end do
        else
            do i = 1, nstep_tot
                if (mod(i, step) .eq. 0) readstep(i) = .true.
            end do
        end if

        ! Now I can readjust and figure out how many steps to actually read, and what the timestep will be.
        sim%nt = 0
        do i = 1, nstep_tot
            if (readstep(i)) sim%nt = sim%nt + 1
        end do
        sim%timestep = sim%timestep*step

        ! Space for trajectories, some are always there
        lo_allocate(sim%r(3, sim%na, sim%nt))
        lo_allocate(sim%u(3, sim%na, sim%nt))
        lo_allocate(sim%f(3, sim%na, sim%nt))
        sim%r = 0.0_r8
        sim%u = 0.0_r8
        sim%f = 0.0_r8
        ! and some things are sometimes needed
        if (sim%stochastic .eqv. .false.) then
            lo_allocate(sim%v(3, sim%na, sim%nt))
            lo_allocate(sim%r_npbc(3, sim%na, sim%nt))
            sim%v = 0.0_r8
            sim%r_npbc = 0.0_r8
        end if
        if (sim%have_magnetic_moments) then
            lo_allocate(sim%m(3, sim%na, sim%nt))
            sim%m = 0.0_r8
        end if
        if (sim%have_dielectric) then
            allocate (sim%Z(3, 3, sim%na, sim%nt))
            allocate (sim%eps(3, 3, sim%nt))
            sim%Z = 0.0_r8
            sim%eps = 0.0_r8
        end if
        if (sim%variable_lattice) then
            allocate (sim%lattice(3, 3, sim%nt))
            sim%lattice = 0.0_r8
        end if
        if (present(readalloy)) then
        if (readalloy) then
            lo_allocate(idealpos(3, sim%na, sim%nt))
            lo_allocate(refpos(3, sim%na, sim%nt))
            lo_allocate(atomic_numbers(sim%na, sim%nt))
            idealpos = 0.0_r8
            refpos = 0.0_r8
            atomic_numbers = 0
        end if
        end if

        ! And some space for energies
        allocate (sim%stat%potential_energy(sim%nt))
        allocate (sim%stat%kinetic_energy(sim%nt))
        allocate (sim%stat%total_energy(sim%nt))
        allocate (sim%stat%temperature(sim%nt))
        allocate (sim%stat%pressure(sim%nt))
        allocate (sim%stat%magnetic_potential_energy(sim%nt))
        allocate (sim%stat%polar_potential_energy(sim%nt))
        allocate (sim%stat%secondorder_potential_energy(sim%nt))
        allocate (sim%stat%thirdorder_potential_energy(sim%nt))
        allocate (sim%stat%fourthorder_potential_energy(sim%nt))
        allocate (sim%stat%stress(3, 3, sim%nt))
        sim%stat%potential_energy = 0.0_r8
        sim%stat%kinetic_energy = 0.0_r8
        sim%stat%total_energy = 0.0_r8
        sim%stat%temperature = 0.0_r8
        sim%stat%pressure = 0.0_r8
        sim%stat%magnetic_potential_energy = 0.0_r8
        sim%stat%polar_potential_energy = 0.0_r8
        sim%stat%secondorder_potential_energy = 0.0_r8
        sim%stat%thirdorder_potential_energy = 0.0_r8
        sim%stat%fourthorder_potential_energy = 0.0_r8
        sim%stat%stress = 0.0_r8

        if (verbosity .gt. 0) write (*, *) '... made space and decided on heuristics'
    end block init

    ! Actually read things
    readnormal: block
        integer :: i, j, k, l, u, nmaxread, nroffset
        real(r8), dimension(3, 3) :: m0, m1, I3
        real(r8), dimension(6) :: sigma
        real(r8) :: f0, t0

        I3 = 0.0_r8
        do i = 1, 3
            I3(i, i) = 1.0_r8
        end do

        ! Figure out how many lines I have to read, in total
        nmaxread = nstep_tot*3
        if (sim%have_magnetic_moments) nmaxread = nmaxread + nstep_tot
        if (sim%have_dielectric) nmaxread = nmaxread + nstep_tot*2
        if (sim%variable_lattice) nmaxread = nmaxread + nstep_tot
        ! Only read on one rank.
        if (readonthisrank) then
            t0 = walltime()
            if (verbosity .gt. 0) call lo_progressbar_init()
            nroffset = 0
            u = open_file('in', 'infile.positions')
            i = 0
            do l = 1, nstep_tot
                if (readstep(l)) then
                    i = i + 1
                    do j = 1, sim%na
                        read (u, *) sim%r(:, j, i)
                    end do
                else
                    do j = 1, sim%na
                        read (u, *)
                    end do
                end if
                if (verbosity .gt. 0) then
                if (lo_trueNtimes(l, 25, nstep_tot)) then
                    call lo_progressbar('... reading simulation', nroffset + l, nmaxread, walltime() - t0)
                end if
                end if
            end do
            close (u)
            ! read forces
            nroffset = nroffset + nstep_tot
            u = open_file('in', 'infile.forces')
            i = 0
            do l = 1, nstep_tot
                if (readstep(l)) then
                    i = i + 1
                    do j = 1, sim%na
                        read (u, *) sim%f(:, j, i)
                    end do
                else
                    do j = 1, sim%na
                        read (u, *)
                    end do
                end if
                if (verbosity .gt. 0) then
                if (lo_trueNtimes(l, 25, nstep_tot)) then
                    call lo_progressbar('... reading simulation', nroffset + l, nmaxread, walltime() - t0)
                end if
                end if
            end do
            close (u)
            ! read stat
            nroffset = nroffset + nstep_tot
            u = open_file('in', 'infile.stat')
            i = 0
            do l = 1, nstep_tot
                if (readstep(l)) then
                    i = i + 1
                    read (u, *) j, f0, sim%stat%total_energy(i), sim%stat%potential_energy(i), &
                        sim%stat%kinetic_energy(i), sim%stat%temperature(i), sim%stat%pressure(i), sigma
                    sim%stat%stress(1, 1, i) = sigma(1)
                    sim%stat%stress(2, 2, i) = sigma(2)
                    sim%stat%stress(3, 3, i) = sigma(3)
                    sim%stat%stress(1, 2, i) = sigma(4)
                    sim%stat%stress(2, 1, i) = sigma(4)
                    sim%stat%stress(2, 3, i) = sigma(5)
                    sim%stat%stress(3, 2, i) = sigma(5)
                    sim%stat%stress(1, 3, i) = sigma(6)
                    sim%stat%stress(3, 1, i) = sigma(6)
                else
                    read (u, *)
                end if
                if (verbosity .gt. 0) then
                if (lo_trueNtimes(l, 25, nstep_tot)) then
                    call lo_progressbar('... reading simulation', nroffset + l, nmaxread, walltime() - t0)
                end if
                end if
            end do
            close (u)
            if (sim%have_magnetic_moments) then
                ! read magnetic moments
                nroffset = nroffset + nstep_tot
                u = open_file('in', 'infile.projected_moments')
                i = 0
                do l = 1, nstep_tot
                    if (readstep(l)) then
                        i = i + 1
                        do j = 1, sim%na
                            read (u, *) sim%m(:, j, i)
                        end do
                    else
                        do j = 1, sim%na
                            read (u, *)
                        end do
                    end if
                    if (verbosity .gt. 0) then
                    if (lo_trueNtimes(l, 25, nstep_tot)) then
                        call lo_progressbar('... reading simulation', nroffset + l, nmaxread, walltime() - t0)
                    end if
                    end if
                end do
                close (u)
            end if
            if (sim%variable_lattice) then
                ! read lattice vectors
                nroffset = nroffset + nstep_tot
                u = open_file('in', 'infile.lattice')
                i = 0
                do l = 1, nstep_tot
                    if (readstep(l)) then
                        i = i + 1
                        do j = 1, 3
                            read (u, *) sim%lattice(:, j, i)
                        end do
                        sim%lattice(:, :, i) = sim%lattice(:, :, i)*lo_A_to_bohr
                    else
                        do j = 1, 3
                            read (u, *)
                        end do
                    end if
                    if (verbosity .gt. 0) then
                    if (lo_trueNtimes(l, 25, nstep_tot)) then
                        call lo_progressbar('... reading simulation', nroffset + l, nmaxread, walltime() - t0)
                    end if
                    end if
                end do
                close (u)
            end if
            if (sim%have_dielectric) then
                ! read born charges and dielectric constants
                nroffset = nroffset + nstep_tot
                u = open_file('in', 'infile.born_charges')
                i = 0
                do l = 1, nstep_tot
                    if (readstep(l)) then
                        i = i + 1
                        do j = 1, sim%na
                        do k = 1, 3
                            read (u, *) sim%Z(:, k, j, i)
                            !read(u,*) sim%Z(k,:,j,i)
                        end do
                        end do
                    else
                        do j = 1, sim%na
                        do k = 1, 3
                            read (u, *)
                        end do
                        end do
                    end if
                    if (verbosity .gt. 0) then
                    if (lo_trueNtimes(l, 25, nstep_tot)) then
                        call lo_progressbar('... reading simulation', nroffset + l, nmaxread, walltime() - t0)
                    end if
                    end if
                end do
                close (u)
                nroffset = nroffset + nstep_tot
                u = open_file('in', 'infile.dielectric_tensor')
                i = 0
                do l = 1, nstep_tot
                    if (readstep(l)) then
                        i = i + 1
                        do j = 1, 3
                            read (u, *) sim%eps(:, j, i)
                            !read(u,*) m0(:,j)
                        end do
                        ! Convert dielectric tensor to derivative wrt electric field
                        !m1=(I3-m0)*sim%crystalstructure%volume/4/lo_pi
                        !sim%eps(:,:,i)=m1
                    else
                        do j = 1, 3
                            read (u, *)
                        end do
                    end if
                    if (verbosity .gt. 0) then
                    if (lo_trueNtimes(l, 25, nstep_tot)) then
                        call lo_progressbar('... reading simulation', nroffset + l, nmaxread, walltime() - t0)
                    end if
                    end if
                end do
                close (u)
            end if
        end if ! end if readrank
        ! Now all positions and things are read. Send data to all ranks, if needed
        if (mpiparallel) then
            call mw%bcast(sim%r, readrank, __FILE__, __LINE__)
            call mw%bcast(sim%f, readrank, __FILE__, __LINE__)
            if (sim%variable_lattice) call mw%bcast(sim%lattice, readrank, __FILE__, __LINE__)
            if (sim%have_magnetic_moments) call mw%bcast(sim%m, readrank, __FILE__, __LINE__)
            if (sim%have_dielectric) call mw%bcast(sim%Z, readrank, __FILE__, __LINE__)
            if (sim%have_dielectric) call mw%bcast(sim%eps, readrank, __FILE__, __LINE__)
            call mw%bcast(sim%stat%temperature, readrank, __FILE__, __LINE__)
            call mw%bcast(sim%stat%total_energy, readrank, __FILE__, __LINE__)
            call mw%bcast(sim%stat%potential_energy, readrank, __FILE__, __LINE__)
            call mw%bcast(sim%stat%kinetic_energy, readrank, __FILE__, __LINE__)
            call mw%bcast(sim%stat%pressure, readrank, __FILE__, __LINE__)
            call mw%bcast(sim%stat%stress, readrank, __FILE__, __LINE__)
        end if
        ! And say we are done reading
        if (verbosity .gt. 0) call lo_progressbar('... reading simulation', nmaxread, nmaxread, walltime() - t0)

        ! Convert units right away. First ensure things are inside the 0-1 box
        sim%r = lo_clean_fractional_coordinates(sim%r)
        sim%timestep = sim%timestep*lo_time_fs_to_au
        sim%f = sim%f*lo_force_eVA_to_HartreeBohr
        sim%stat%total_energy = sim%stat%total_energy*lo_eV_to_Hartree
        sim%stat%potential_energy = sim%stat%potential_energy*lo_eV_to_Hartree
        sim%stat%kinetic_energy = sim%stat%kinetic_energy*lo_eV_to_Hartree
        sim%stat%pressure = sim%stat%pressure*lo_pressure_GPa_to_HartreeBohr
        sim%stat%stress = sim%stat%stress*lo_pressure_GPa_to_HartreeBohr
    end block readnormal

    fixnormal: block
        if (sim%stochastic .eqv. .false.) then
            ! get non-pbc positions
            call sim%get_nonpbc_positions()
            if (verbosity .gt. 0) write (*, *) '... unwrapped positions'
            ! get velocities
            call sim%get_velocities()
            if (verbosity .gt. 0) write (*, *) '... got velocities'
        end if
        ! get displacements
        call sim%get_displacements()
        if (verbosity .gt. 0) write (*, *) '... got displacements'
    end block fixnormal

    ! Calculate some things
    finalize: block
        real(r8), dimension(3, 3) :: m0
        real(r8) :: avg_pressure, avg_temperature

        integer :: i, j
        character(len=1000) :: opf
        ! get averages of stat things
        avg_pressure = lo_mean(sim%stat%pressure)*lo_pressure_HartreeBohr_to_GPa
        avg_temperature = lo_mean(sim%stat%temperature)
        do i = 1, 3
        do j = 1, 3
            m0(j, i) = lo_mean(sim%stat%stress(j, i, :))
        end do
        end do
        m0 = m0*lo_pressure_HartreeBohr_to_GPa
        ! all done, perhaps dump some info
        if (verbosity .gt. 0) then
            write (*, *) '... short summary of simulation:'
            write (*, *) '                      number of atoms: ', tochar(sim%na)
            write (*, *) '             number of configurations: ', tochar(sim%nt)
            write (*, *) '                        timestep (fs): ', tochar(sim%timestep*lo_time_au_to_fs)
            write (*, *) '              average temperature (K): ', tochar(avg_temperature)
            write (*, *) '                thermostat set to (K): ', tochar(sim%temperature_thermostat)
            write (*, *) '               average pressure (GPa): ', tochar(avg_pressure)
            opf = "(1X,'  average stresstensor xx xy xz (GPa): ',3(F12.4,' ') )"
            write (*, opf) m0(:, 1)
            opf = "(1X,'  average stresstensor yx yy yz (GPa): ',3(F12.4,' ') )"
            write (*, opf) m0(:, 2)
            opf = "(1X,'  average stresstensor zx zy zz (GPa): ',3(F12.4,' ') )"
            write (*, opf) m0(:, 3)
            write (*, *) '... finished reading simulation (', tochar(walltime() - timer), 's)'
        end if
    end block finalize
end subroutine

!> write a simulation to hdf5
subroutine write_to_hdf5(sim, uc, ss, filename, verbosity, eps, Z)
    !> md simulation
    class(lo_mdsim), intent(in) :: sim
    !> unitcell
    type(lo_crystalstructure), intent(in) :: uc
    !> supercell
    type(lo_crystalstructure), intent(inout) :: ss
    !> filename
    character(len=*), intent(in) :: filename
    !> Talk a lot?
    integer, intent(in) :: verbosity
    !> dielectric tensor
    real(r8), dimension(3, 3), intent(in), optional :: eps
    !> Born effective charges
    real(r8), dimension(:, :, :), intent(in), optional :: Z

    real(r8), dimension(3, 3) :: buf_eps
    real(r8), dimension(:, :, :), allocatable :: buf_Z
    real(r8) :: timer

    init: block
        if (verbosity .gt. 0) then
            timer = walltime()
            write (*, *) ''
            write (*, *) 'Writing simulation to "'//trim(filename)//'"'
        end if

        ! Check that the supercell really is a supercell
        call ss%classify('supercell', uc)

        ! Buffers for Born charges and dielectric tensor
        if (present(eps)) then
            buf_eps = eps
        else
            buf_eps = 0.0_r8
        end if
        lo_allocate(buf_Z(3, 3, uc%na))
        if (present(Z)) then
            buf_Z = Z
        else
            buf_Z = 0.0_r8
        end if
    end block init

    writefile: block
        type(lo_hdf5_helper) :: h5
        real(r8), dimension(:, :, :, :), allocatable :: dw
        real(r8), dimension(:, :, :), allocatable :: dr
        real(r8), dimension(:), allocatable :: ds
        ! Initialize hdf5 properly
        call h5%init(__FILE__, __LINE__)
        call h5%open_file('write', trim(filename))

        ! Store some metadata. Not sure if effective.
        call lo_h5_store_attribute(sim%na, h5%file_id, 'number_of_atoms')
        call lo_h5_store_attribute(sim%nt, h5%file_id, 'number_of_timesteps')
        call lo_h5_store_attribute(sim%timestep*lo_time_au_to_fs, h5%file_id, 'timestep_in_fs')
        call lo_h5_store_attribute(sim%temperature_thermostat, h5%file_id, 'temperature_thermostat')
        call lo_h5_store_attribute(sim%stochastic, h5%file_id, 'is_stampling_stochastic')
        call lo_h5_store_attribute(sim%alloy, h5%file_id, 'is_simulation_alloy')

        ! Write positions
        lo_allocate(dr(3, sim%na, sim%nt))
        dr = sim%r(:, :, 1:sim%nt)
        call lo_h5_store_data(dr, h5%file_id, 'positions', enhet='fractional')
        if (verbosity .gt. 0) write (*, *) '... wrote positions'

        ! Write displacements
        dr = sim%u(:, :, 1:sim%nt)*lo_bohr_to_A
        call lo_h5_store_data(dr, h5%file_id, 'displacements', enhet='A')
        lo_deallocate(dr)
        if (verbosity .gt. 0) write (*, *) '... wrote displacements'

        ! Write velocities
        if (allocated(sim%v)) then
            lo_allocate(dr(3, sim%na, sim%nt))
            dr = sim%v(:, :, 1:sim%nt)*lo_velocity_au_to_Afs
            call lo_h5_store_data(dr, h5%file_id, 'velocities', enhet='A/fs')
            lo_deallocate(dr)
            if (verbosity .gt. 0) write (*, *) '... wrote velocities'
        end if

        ! Write forces
        lo_allocate(dr(3, sim%na, sim%nt))
        dr = sim%f(:, :, 1:sim%nt)*lo_force_Hartreebohr_to_eVA
        call lo_h5_store_data(dr, h5%file_id, 'forces', enhet='eV/A')
        lo_deallocate(dr)
        if (verbosity .gt. 0) write (*, *) '... wrote forces'

        ! Write magnetic moment per atom, if applicable
        if (sim%have_magnetic_moments) then
            lo_allocate(dr(3, sim%na, sim%nt))
            dr = sim%m(:, :, 1:sim%nt)
            call lo_h5_store_data(dr, h5%file_id, 'magnetic_moments', enhet='mu_B')
            lo_deallocate(dr)
            if (verbosity .gt. 0) write (*, *) '... wrote magnetic moments'
        end if

        ! Write the time-dependent lattice, if applicable
        if (sim%variable_lattice) then
            lo_allocate(dr(3, 3, sim%nt))
            dr = sim%lattice*lo_bohr_to_A
            call lo_h5_store_data(dr, h5%file_id, 'lattice', enhet='A')
            lo_deallocate(dr)
            if (verbosity .gt. 0) write (*, *) '... wrote lattice'
        end if

        ! Write time-dependent Born charges and dielectric constant
        if (sim%have_dielectric) then
            lo_allocate(dr(3, 3, sim%nt))
            dr = sim%eps(:, :, 1:sim%nt)
            call lo_h5_store_data(dr, h5%file_id, 'eps', enhet='dimensionless')
            lo_deallocate(dr)

            lo_allocate(dw(3, 3, sim%na, sim%nt))
            dw = sim%Z(:, :, :, 1:sim%nt)
            call lo_h5_store_data(dw, h5%file_id, 'Z', enhet='charge/bohr')
            lo_deallocate(dw)

            if (verbosity .gt. 0) write (*, *) '... wrote dielectric constant and Born charge'
        end if

        ! Write some alloy things?
        if (sim%alloy) then
            call lo_h5_store_data(sim%r_ref, h5%file_id, 'reference_positions', enhet='fractional')
            call lo_h5_store_data(sim%atomic_numbers, h5%file_id, 'atomic_numbers', enhet='Z')
            ! Then the alloy specification
            call lo_h5_store_data(sim%extra%unitcell_componentcounter, h5%file_id, 'unitcell_componentcounter')
            call lo_h5_store_data(sim%extra%supercell_componentcounter, h5%file_id, 'supercell_componentcounter')
            call lo_h5_store_data(sim%extra%unitcell_components, h5%file_id, 'unitcell_components')
            call lo_h5_store_data(sim%extra%supercell_components, h5%file_id, 'supercell_components')
            call lo_h5_store_data(sim%extra%unitcell_concentrations, h5%file_id, 'unitcell_concentrations')
            call lo_h5_store_data(sim%extra%supercell_concentrations, h5%file_id, 'supercell_concentrations')
        end if

        ! Write a lot of energies
        lo_allocate(ds(sim%nt))
        ds = sim%stat%total_energy(1:sim%nt)*lo_Hartree_to_eV
        call lo_h5_store_data(ds, h5%file_id, 'total_energy', enhet='eV')
        ds = sim%stat%potential_energy(1:sim%nt)*lo_Hartree_to_eV
        call lo_h5_store_data(ds, h5%file_id, 'potential_energy', enhet='eV')
        ds = sim%stat%kinetic_energy(1:sim%nt)*lo_Hartree_to_eV
        call lo_h5_store_data(ds, h5%file_id, 'kinetic_energy', enhet='eV')
        ds = sim%stat%temperature(1:sim%nt)
        call lo_h5_store_data(ds, h5%file_id, 'temperature', enhet='K')
        ds = sim%stat%pressure(1:sim%nt)*lo_pressure_HartreeBohr_to_GPa
        call lo_h5_store_data(ds, h5%file_id, 'pressure', enhet='GPa')
        lo_deallocate(ds)
        lo_allocate(dr(3, 3, sim%nt))
        dr = sim%stat%stress(:, :, 1:sim%nt)*lo_pressure_HartreeBohr_to_GPa
        call lo_h5_store_data(dr, h5%file_id, 'stress', enhet='GPa')
        lo_deallocate(dr)

        ! Maybe some auxiliary stuff
        call lo_h5_store_data(buf_Z, h5%file_id, 'born_effective_charges', enhet='e/A')
        call lo_h5_store_data(buf_eps, h5%file_id, 'dielectric_tensor', enhet='dimensionless')
        call lo_h5_store_data(uc%latticevectors*lo_bohr_to_A, h5%file_id, 'unitcell_latticevectors', enhet='A')
        call lo_h5_store_data(ss%latticevectors*lo_bohr_to_A, h5%file_id, 'supercell_latticevectors', enhet='A')
        call lo_h5_store_data(uc%r, h5%file_id, 'unitcell_positions', enhet='dimensionless')
        call lo_h5_store_data(ss%r, h5%file_id, 'supercell_positions', enhet='dimensionless')
        call lo_h5_store_data(uc%atomic_number, h5%file_id, 'unitcell_atomic_numbers', enhet='e')
        call lo_h5_store_data(ss%atomic_number, h5%file_id, 'supercell_atomic_numbers', enhet='e')
        if (verbosity .gt. 0) write (*, *) '... wrote energies and metadata'
        ! ! And close
        call h5%close_file()
        call h5%destroy(__FILE__, __LINE__)
    end block writefile

    if (verbosity .gt. 0) write (*, *) 'wrote simulation (', tochar(walltime() - timer), 's)'
end subroutine

!> dump simulation to plain text @TODO add lattice, magmom, and other things.
subroutine write_to_plaintext(sim, filepref)
    !> md simulation
    class(lo_mdsim), intent(in) :: sim
    !> prefix for filenames
    character(len=*), intent(in) :: filepref

    real(r8), dimension(6) :: sigma
    real(r8) :: f0
    integer :: u, i, j

    u = open_file('out', trim(filepref)//'.meta')
    write (u, *) sim%na
    write (u, *) sim%nt
    write (u, *) sim%timestep*lo_time_au_to_fs
    write (u, *) sim%temperature_thermostat
    close (u)

    u = open_file('out', trim(filepref)//'.stat')
    do i = 1, sim%nt
        sigma(1) = sim%stat%stress(1, 1, i)
        sigma(2) = sim%stat%stress(2, 2, i)
        sigma(3) = sim%stat%stress(3, 3, i)
        sigma(4) = sim%stat%stress(2, 1, i)
        sigma(5) = sim%stat%stress(2, 3, i)
        sigma(6) = sim%stat%stress(3, 1, i)
        sigma = sigma*lo_pressure_HartreeBohr_to_GPa
        f0 = lo_trace(sim%stat%stress(:, :, i))*lo_pressure_HartreeBohr_to_GPa/3.0_r8

        write (u, "(I6,' ',12(2X,EN16.7))") &
            i, &
            i*sim%timestep*lo_time_au_to_fs, &
            sim%stat%total_energy(i)*lo_Hartree_to_eV, &
            sim%stat%potential_energy(i)*lo_Hartree_to_eV, &
            sim%stat%kinetic_energy(i)*lo_Hartree_to_eV, &
            sim%stat%temperature(i), &
            f0, &
            sigma
    end do
    close (u)

    u = open_file('out', trim(filepref)//'.positions')
    do i = 1, sim%nt
        do j = 1, sim%na
            write (u, *) sim%r(:, j, i)
        end do
    end do
    close (u)

    u = open_file('out', trim(filepref)//'.forces')
    do i = 1, sim%nt
    do j = 1, sim%na
        write (u, *) sim%f(:, j, i)*lo_force_HartreeBohr_to_eVA
    end do
    end do
    close (u)

    ! if ( sim%alloy ) then
    !     u=open_file('out',trim(filepref)//'.alloy_atomic_numbers')
    !         do i=1,sim%nt
    !         do j=1,sim%na
    !             write(u,*) sim%atomic_numbers(j,i)
    !         enddo
    !         enddo
    !     close(u)
    !
    !     u=open_file('out',trim(filepref)//'.alloy_reference_positions')
    !         do i=1,sim%nt
    !         do j=1,sim%na
    !             write(u,*) sim%r_ref(:,j,i)
    !         enddo
    !         enddo
    !     close(u)
    ! endif
end subroutine

!> read simulation from the hdf5 format
subroutine read_from_hdf5(sim, filename, verbosity, stride, nrand, mw)
    !> md simulation
    class(lo_mdsim), intent(out) :: sim
    !> filename
    character(len=*), intent(in) :: filename
    !> how much to talk
    integer, intent(in) :: verbosity
    !> read it in strides
    integer, intent(in), optional :: stride
    !> read N random timesteps
    integer, intent(in), optional :: nrand
    !> MPI communicator
    type(lo_mpi_helper), intent(inout), optional :: mw

    type(lo_mersennetwister) :: tw
    real(r8) :: timer
    integer :: step, nrnd, readrank
    logical :: readonthisrank, mpiparallel

    ! decide on some things
    init: block
        if (verbosity .gt. 0) then
            timer = walltime()
            write (*, *) ''
            write (*, *) 'READING SIMULATION FROM "'//trim(filename)//'"'
        end if

        ! Which rank will do all the reading?
        if (present(mw)) then
            ! We are reading on one rank, then broadcasting
            readrank = 0
            if (mw%r .eq. readrank) then
                readonthisrank = .true.
            else
                readonthisrank = .false.
            end if
            mpiparallel = .true.
        else
            ! Reading not in parallel.
            readrank = 0
            readonthisrank = .true.
            mpiparallel = .false.
        end if
        ! read in strides?
        if (present(stride)) then
            step = stride
        else
            step = 1
        end if
        ! read only a few at random
        if (present(nrand)) then
            nrnd = nrand
            call tw%init(iseed=0, rseed=walltime())
        else
            nrnd = -1
        end if

        ! Just some early sanity checks
        if (mpiparallel) then
            call mw%check_and_sync(readrank, readrank, 2, 'readrank', __FILE__, __LINE__)
            call mw%check_and_sync(step, readrank, 2, 'step', __FILE__, __LINE__)
            call mw%check_and_sync(nrnd, readrank, 2, 'nrnd', __FILE__, __LINE__)
        end if
    end block init

    ! start reading
    readstuff: block
        type(lo_hdf5_helper) :: h5
        real(r8), dimension(:, :, :, :), allocatable :: dw
        real(r8), dimension(:, :, :), allocatable :: dr
        real(r8), dimension(:), allocatable :: ds
        real(r8) :: t0
        integer, dimension(:, :), allocatable :: dj
        integer, dimension(:), allocatable :: di
        integer :: i, l, nstep_tot
        logical, dimension(:), allocatable :: readstep

        ! Initialize hdf5 properly
        call h5%init(__FILE__, __LINE__)
        call h5%open_file('read', trim(filename))

        ! Start reading some metadata, can do this on all ranks since it's very little
        ! data and I don't have to think that much.
        call lo_h5_read_attribute(nstep_tot, h5%file_id, 'number_of_timesteps')
        call lo_h5_read_attribute(sim%na, h5%file_id, 'number_of_atoms')
        call lo_h5_read_attribute(sim%timestep, h5%file_id, 'timestep_in_fs')
        call lo_h5_read_attribute(sim%temperature_thermostat, h5%file_id, 'temperature_thermostat')
        call lo_h5_read_attribute(sim%stochastic, h5%file_id, 'is_stampling_stochastic')
        call lo_h5_read_attribute(sim%alloy, h5%file_id, 'is_simulation_alloy')
        sim%have_magnetic_moments = lo_h5_does_dataset_exist(h5%file_id, 'magnetic_moments')
        sim%variable_lattice = lo_h5_does_dataset_exist(h5%file_id, 'lattice')
        sim%have_dielectric = lo_h5_does_dataset_exist(h5%file_id, 'Z')
        call lo_h5_read_data(sim%extra%born_effective_charges, h5%file_id, 'born_effective_charges')
        call lo_h5_read_data(sim%extra%dielectric_tensor, h5%file_id, 'dielectric_tensor')
        call lo_h5_read_data(sim%extra%unitcell_latticevectors, h5%file_id, 'unitcell_latticevectors')
        call lo_h5_read_data(sim%extra%supercell_latticevectors, h5%file_id, 'supercell_latticevectors')
        call lo_h5_read_data(sim%extra%unitcell_positions, h5%file_id, 'unitcell_positions')
        call lo_h5_read_data(sim%extra%supercell_positions, h5%file_id, 'supercell_positions')
        call lo_h5_read_data(sim%extra%unitcell_atomic_numbers, h5%file_id, 'unitcell_atomic_numbers')
        call lo_h5_read_data(sim%extra%supercell_atomic_numbers, h5%file_id, 'supercell_atomic_numbers')
        sim%extra%unitcell_latticevectors = sim%extra%unitcell_latticevectors*lo_A_to_bohr
        sim%extra%supercell_latticevectors = sim%extra%supercell_latticevectors*lo_A_to_bohr
        if (sim%alloy) then
            call lo_h5_read_data(sim%extra%unitcell_componentcounter, h5%file_id, 'unitcell_componentcounter')
            call lo_h5_read_data(sim%extra%supercell_componentcounter, h5%file_id, 'supercell_componentcounter')
            call lo_h5_read_data(sim%extra%unitcell_components, h5%file_id, 'unitcell_components')
            call lo_h5_read_data(sim%extra%supercell_components, h5%file_id, 'supercell_components')
            call lo_h5_read_data(sim%extra%unitcell_concentrations, h5%file_id, 'unitcell_concentrations')
            call lo_h5_read_data(sim%extra%supercell_concentrations, h5%file_id, 'supercell_concentrations')
        end if

        if (mpiparallel) then
            call mw%check_and_sync(nstep_tot, readrank, 2, 'nstep_tot', __FILE__, __LINE__)
            call mw%check_and_sync(sim%na, readrank, 2, 'sim%na', __FILE__, __LINE__)
            call mw%check_and_sync(sim%timestep, readrank, 2, 'sim%timestep', __FILE__, __LINE__)
        end if
        ! Now I can figure out the actual number of timesteps I'm going to read
        lo_allocate(readstep(nstep_tot))
        readstep = .false.
        if (nrnd .gt. 0) then
            allocate (di(nstep_tot))
            do i = 1, nstep_tot
                di(i) = i
            end do
            ! shuffle timesteps
            call tw%shuffle_int_array(di)
            if (mpiparallel) then
                call mw%bcast(di, readrank, __FILE__, __LINE__)
            end if
            ! decide which ones to read
            do i = 1, nrnd
                readstep(di(i)) = .true.
            end do
            deallocate (di)
        else
            do i = 1, nstep_tot
                if (mod(i, step) .eq. 0) readstep(i) = .true.
            end do
        end if
        sim%nt = 0
        do i = 1, nstep_tot
            if (readstep(i)) sim%nt = sim%nt + 1
        end do
        sim%timestep = sim%timestep*step

        ! Get an index array that denote which steps I want
        ! I should figure out how to read from hdf5 with strides, but it makes my head hurt
        lo_allocate(di(sim%nt))
        di = 0
        l = 0
        do i = 1, nstep_tot
            if (readstep(i)) then
                l = l + 1
                di(l) = i
            end if
        end do

        ! Space for trajectories, some are always there
        lo_allocate(sim%r(3, sim%na, sim%nt))
        lo_allocate(sim%u(3, sim%na, sim%nt))
        lo_allocate(sim%f(3, sim%na, sim%nt))
        sim%r = 0.0_r8
        sim%u = 0.0_r8
        sim%f = 0.0_r8
        ! and some things are sometimes needed
        if (sim%stochastic .eqv. .false.) then
            lo_allocate(sim%v(3, sim%na, sim%nt))
            lo_allocate(sim%r_npbc(3, sim%na, sim%nt))
            sim%v = 0.0_r8
            sim%r_npbc = 0.0_r8
        end if
        if (sim%have_magnetic_moments) then
            lo_allocate(sim%m(3, sim%na, sim%nt))
            sim%m = 0.0_r8
        end if
        if (sim%variable_lattice) then
            lo_allocate(sim%lattice(3, 3, sim%nt))
            sim%lattice = 0.0_r8
        end if
        if (sim%alloy) then
            lo_allocate(sim%r_ref(3, sim%na, sim%nt))
            lo_allocate(sim%atomic_numbers(sim%na, sim%nt))
            sim%r_ref = 0.0_r8
            sim%atomic_numbers = 0
        end if
        if (sim%have_dielectric) then
            allocate (sim%eps(3, 3, sim%nt))
            allocate (sim%Z(3, 3, sim%na, sim%nt))
            sim%eps = 0.0_r8
            sim%Z = 0.0_r8
        end if

        ! And some space for energies
        lo_allocate(sim%stat%potential_energy(sim%nt))
        lo_allocate(sim%stat%kinetic_energy(sim%nt))
        lo_allocate(sim%stat%total_energy(sim%nt))
        lo_allocate(sim%stat%temperature(sim%nt))
        lo_allocate(sim%stat%pressure(sim%nt))
        lo_allocate(sim%stat%magnetic_potential_energy(sim%nt))
        lo_allocate(sim%stat%polar_potential_energy(sim%nt))
        lo_allocate(sim%stat%secondorder_potential_energy(sim%nt))
        lo_allocate(sim%stat%thirdorder_potential_energy(sim%nt))
        lo_allocate(sim%stat%fourthorder_potential_energy(sim%nt))
        lo_allocate(sim%stat%stress(3, 3, sim%nt))
        sim%stat%potential_energy = 0.0_r8
        sim%stat%kinetic_energy = 0.0_r8
        sim%stat%total_energy = 0.0_r8
        sim%stat%temperature = 0.0_r8
        sim%stat%pressure = 0.0_r8
        sim%stat%magnetic_potential_energy = 0.0_r8
        sim%stat%polar_potential_energy = 0.0_r8
        sim%stat%secondorder_potential_energy = 0.0_r8
        sim%stat%thirdorder_potential_energy = 0.0_r8
        sim%stat%fourthorder_potential_energy = 0.0_r8
        sim%stat%stress = 0.0_r8

        if (verbosity .gt. 0) write (*, *) '... read metadata'

        ! Now actually read
        if (readonthisrank) then
            t0 = walltime()
            call lo_h5_read_data(dr, h5%file_id, 'positions')
            do i = 1, sim%nt
                sim%r(:, :, i) = dr(:, :, di(i))
            end do
            if (verbosity .gt. 0) write (*, *) '... read positions (', tochar(walltime() - t0), ')'

            t0 = walltime()
            call lo_h5_read_data(dr, h5%file_id, 'forces')
            do i = 1, sim%nt
                sim%f(:, :, i) = dr(:, :, di(i))
            end do
            if (verbosity .gt. 0) write (*, *) '... read forces (', tochar(walltime() - t0), ')'

            if (sim%have_magnetic_moments) then
                t0 = walltime()
                call lo_h5_read_data(dr, h5%file_id, 'magnetic_moments')
                do i = 1, sim%nt
                    sim%m(:, :, i) = dr(:, :, di(i))
                end do
                if (verbosity .gt. 0) write (*, *) '... read magnetic moments (', tochar(walltime() - t0), ')'
            end if

            if (sim%variable_lattice) then
                t0 = walltime()
                call lo_h5_read_data(dr, h5%file_id, 'lattice')
                do i = 1, sim%nt
                    sim%lattice(:, :, i) = dr(:, :, di(i))
                end do
                if (verbosity .gt. 0) write (*, *) '... read lattice (', tochar(walltime() - t0), ')'
            end if

            if (sim%alloy) then
                t0 = walltime()
                call lo_h5_read_data(dr, h5%file_id, 'reference_positions')
                do i = 1, sim%nt
                    sim%r_ref(:, :, i) = dr(:, :, di(i))
                end do
                if (verbosity .gt. 0) write (*, *) '... read reference positions (', tochar(walltime() - t0), ')'

                t0 = walltime()
                call lo_h5_read_data(dj, h5%file_id, 'atomic_numbers')
                do i = 1, sim%nt
                    sim%atomic_numbers(:, i) = dj(:, di(i))
                end do
                if (verbosity .gt. 0) write (*, *) '... read atomic number (', tochar(walltime() - t0), ')'
            end if

            if (sim%have_dielectric) then
                t0 = walltime()
                call lo_h5_read_data(dr, h5%file_id, 'eps')
                do i = 1, sim%nt
                    sim%eps(:, :, i) = dr(:, :, di(i))
                end do
                if (verbosity .gt. 0) write (*, *) '... read dielectric tensor (', tochar(walltime() - t0), ')'
                t0 = walltime()
                call lo_h5_read_data(dw, h5%file_id, 'Z')
                do i = 1, sim%nt
                    sim%Z(:, :, :, i) = dw(:, :, :, di(i))
                end do
                if (verbosity .gt. 0) write (*, *) '... read Born charges (', tochar(walltime() - t0), ')'
            end if

            ! And read energies and things like that
            t0 = walltime()
            call lo_h5_read_data(ds, h5%file_id, 'total_energy')
            sim%stat%total_energy = ds(di)
            call lo_h5_read_data(ds, h5%file_id, 'potential_energy')
            sim%stat%potential_energy = ds(di)
            call lo_h5_read_data(ds, h5%file_id, 'kinetic_energy')
            sim%stat%kinetic_energy = ds(di)
            call lo_h5_read_data(ds, h5%file_id, 'temperature')
            sim%stat%temperature = ds(di)
            call lo_h5_read_data(ds, h5%file_id, 'pressure')
            sim%stat%pressure = ds(di)
            call lo_h5_read_data(dr, h5%file_id, 'stress')
            sim%stat%stress = dr(:, :, di)
            if (verbosity .gt. 0) write (*, *) '... read energies (', tochar(walltime() - t0), ')'

            ! Convert to atomic units
            sim%f = sim%f*lo_force_eVA_to_HartreeBohr
            sim%timestep = sim%timestep*lo_time_fs_to_au
            if (sim%variable_lattice) sim%lattice = sim%lattice*lo_A_to_bohr
            sim%stat%total_energy = sim%stat%total_energy*lo_eV_to_Hartree
            sim%stat%potential_energy = sim%stat%potential_energy*lo_eV_to_Hartree
            sim%stat%kinetic_energy = sim%stat%kinetic_energy*lo_eV_to_Hartree
            sim%stat%pressure = sim%stat%pressure*lo_pressure_GPa_to_HartreeBohr
            sim%stat%stress = sim%stat%stress*lo_pressure_GPa_to_HartreeBohr
        end if
        ! If it's parallel I have to spread it across rank
        if (mpiparallel) then
            t0 = walltime()
            call mw%bcast(sim%timestep, readrank, __FILE__, __LINE__)
            call mw%bcast(sim%r, readrank, __FILE__, __LINE__)
            call mw%bcast(sim%f, readrank, __FILE__, __LINE__)
            if (sim%have_magnetic_moments) call mw%bcast(sim%m, readrank, __FILE__, __LINE__)
            if (sim%variable_lattice) call mw%bcast(sim%lattice, readrank, __FILE__, __LINE__)
            call mw%bcast(sim%stat%total_energy, readrank, __FILE__, __LINE__)
            call mw%bcast(sim%stat%potential_energy, readrank, __FILE__, __LINE__)
            call mw%bcast(sim%stat%kinetic_energy, readrank, __FILE__, __LINE__)
            call mw%bcast(sim%stat%temperature, readrank, __FILE__, __LINE__)
            call mw%bcast(sim%stat%pressure, readrank, __FILE__, __LINE__)
            call mw%bcast(sim%stat%stress, readrank, __FILE__, __LINE__)
            if (sim%alloy) call mw%bcast(sim%r_ref, readrank, __FILE__, __LINE__)
            if (sim%alloy) call mw%bcast(sim%atomic_numbers, readrank, __FILE__, __LINE__)
            if (sim%have_dielectric) call mw%bcast(sim%eps, readrank, __FILE__, __LINE__)
            if (sim%have_dielectric) call mw%bcast(sim%Z, readrank, __FILE__, __LINE__)
        end if

        call h5%close_file()
        call h5%destroy(__FILE__, __LINE__)
    end block readstuff

    finalize: block
        real(r8), dimension(3, 3) :: m0
        real(r8) :: avg_pressure, avg_temperature
        real(r8) :: avg_displacement, max_displacement, f0
        integer :: i, j
        character(len=1000) :: opf

        ! Get the reference crystal structure
        if (sim%alloy) then
            call sim%crystalstructure%generate(sim%extra%supercell_latticevectors, sim%extra%supercell_positions, sim%extra%supercell_atomic_numbers, 2, &
                                               alloy=.true., &
                                               alloy_componentcounter=sim%extra%supercell_componentcounter, &
                                               alloy_components=sim%extra%supercell_components, &
                                               alloy_concentrations=sim%extra%supercell_concentrations)
        else
            call sim%crystalstructure%generate(sim%extra%supercell_latticevectors, sim%extra%supercell_positions, sim%extra%supercell_atomic_numbers, 2)
        end if

        if (sim%stochastic) then
            call sim%get_displacements()
        else
            call sim%get_nonpbc_positions()
            call sim%get_velocities()
            call sim%get_displacements()
        end if

        ! get averages of stat things
        avg_pressure = lo_mean(sim%stat%pressure)*lo_pressure_HartreeBohr_to_GPa
        avg_temperature = lo_mean(sim%stat%temperature)
        do i = 1, 3
        do j = 1, 3
            m0(j, i) = lo_mean(sim%stat%stress(j, i, :))
        end do
        end do
        m0 = m0*lo_pressure_HartreeBohr_to_GPa
        avg_displacement = 0.0_r8
        max_displacement = 0.0_r8
        do i = 1, sim%nt
        do j = 1, sim%na
            f0 = norm2(sim%u(:, j, i))
            avg_displacement = avg_displacement + f0
            max_displacement = max(max_displacement, f0)
        end do
        end do
        avg_displacement = avg_displacement/sim%na/sim%nt*lo_bohr_to_A
        max_displacement = max_displacement*lo_bohr_to_A
        ! all done, perhaps dump some info
        if (verbosity .gt. 0) then
            write (*, *) '... short summary of simulation:'
            write (*, *) '                      number of atoms: ', tochar(sim%na)
            write (*, *) '             number of configurations: ', tochar(sim%nt)
            write (*, *) '              average temperature (K): ', tochar(avg_temperature)
            write (*, *) '                thermostat set to (K): ', tochar(sim%temperature_thermostat)
            write (*, *) '          displacement (mean,max) (A): ', tochar([avg_displacement, max_displacement])
            write (*, *) '               average pressure (GPa): ', tochar(avg_pressure)
            opf = "(1X,'  average stresstensor xx xy xz (GPa): ',3(F12.4,' ') )"
            write (*, opf) m0(:, 1)
            opf = "(1X,'  average stresstensor yx yy yz (GPa): ',3(F12.4,' ') )"
            write (*, opf) m0(:, 2)
            opf = "(1X,'  average stresstensor zx zy zz (GPa): ',3(F12.4,' ') )"
            write (*, opf) m0(:, 3)
            write (*, *) '... finished reading simulation (', tochar(walltime() - timer), 's)'
        end if

        ! Sanity check that we are synced over MPI
        if (mpiparallel) then
            i = size(sim%v)
            call mw%check_and_sync(i, readrank, 2, 'size(sim%v)', __FILE__, __LINE__)
            call mw%check_and_sync(sim%timestep, readrank, 2, 'sim%timestep', __FILE__, __LINE__)
            call mw%check_and_sync(avg_pressure, readrank, 2, 'avg_pressure', __FILE__, __LINE__)
            call mw%check_and_sync(avg_temperature, readrank, 2, 'avg_temperature', __FILE__, __LINE__)
            call mw%check_and_sync(sim%nt, readrank, 2, 'sim%nt', __FILE__, __LINE__)
            call mw%check_and_sync(sim%na, readrank, 2, 'sim%na', __FILE__, __LINE__)
        end if
    end block finalize

    if (verbosity .gt. 0) write (*, *) '... read and processed simulation (', tochar(walltime() - timer), ')'
end subroutine

end module
