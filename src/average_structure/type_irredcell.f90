#include "precompilerdefinitions"
!> Get the symmetry operations for a structure in a neat way.
module type_irredcell
! My stuff
use konstanter, only: flyt, lo_tol, lo_sqtol, lo_pi, lo_huge, lo_hugeint, lo_exitcode_symmetry, lo_status, &
                      lo_pressure_HartreeBohr_to_GPa, lo_force_HartreeBohr_to_eVA, &
                      lo_pressure_GPa_to_HartreeBohr, lo_pressure_eVA_to_HartreeBohr, &
                      lo_pressure_HartreeBohr_to_eVA, lo_Hartree_to_eV, lo_eV_to_Hartree, &
                      lo_bohr_to_A, lo_time_au_to_fs
use dump_data
use geometryfunctions
use gottochblandat !, only : lo_identitymatrix
use type_crystalstructure, only: lo_crystalstructure
use type_blas_lapack_wrappers, only: lo_dgesvd, lo_dgels, lo_dgglse, lo_dgelss, lo_dsyev
use type_mdsim, only: lo_mdsim
use type_symmetryoperation
implicit none

private
public :: lo_irredcell
!public :: irreducible_trajectory
public :: avg_structure
public :: avg_structure_nvt
public :: remove_rotations

!> symmetry-reduced representation of the unit cell
type lo_irredcell
    !> number of degrees of freedom of the lattice
    integer :: nx_lattice = -lo_hugeint
    !> number of degrees of freedom for the forces
    integer :: nx_internal = -lo_hugeint
    !> coefficient matrix for the lattice
    real(flyt), dimension(:, :), allocatable :: coeffM_lattice
    !> coefficient matrix for the forces
    real(flyt), dimension(:, :), allocatable :: coeffM_internal
    !> reference lattice
    real(flyt), dimension(3, 3) :: reference_lattice = -lo_huge
    !> reference positions
    real(flyt), dimension(:, :), allocatable :: reference_positions
    !> irrep for the lattice
    real(flyt), dimension(:), allocatable :: x_lattice
    !> irrep for the positions
    real(flyt), dimension(:), allocatable :: x_internal
contains
    !> create everything
    procedure :: generate
    !> transform stress
    !procedure :: strain_to_irreducible

    !> Store a configuration in the trajectory
    !procedure :: add_step
    !> Guess new structure
    !procedure :: predict_new_structure
    !> symmetrize stresses and forces
    !procedure :: symmetrize_stress_and_forces
end type
contains

! More stuff
#include "type_irredcell_aux.f90"

!> Get the average of the internal structure
subroutine avg_structure_nvt(uc, ss, ic, sim, nuc, nss, ur)
    !> unitcell
    type(lo_crystalstructure), intent(in) :: uc
    !> supercell
    type(lo_crystalstructure), intent(in) :: ss
    !> irreducible cell
    type(lo_irredcell), intent(inout) :: ic
    !> simulation
    type(lo_mdsim), intent(inout) :: sim
    !> averaged unitcell
    type(lo_crystalstructure), intent(out) :: nuc
    !> averaged supercell
    type(lo_crystalstructure), intent(out) :: nss
    !> Reference unitcell
    type(lo_crystalstructure), intent(in), optional :: ur

    real(flyt), dimension(3, uc%na) :: average_positions
    real(flyt), dimension(3, 3) :: average_lattice

    ! Align the simulation with the reference positions
    init: block
        real(flyt), dimension(:, :), allocatable :: dr
        real(flyt), dimension(3*uc%na) :: flatr
        real(flyt), dimension(3, 3) :: m0
        real(flyt), dimension(3) :: v0, v1
        real(flyt), dimension(ic%nx_internal) :: xavg
        real(flyt) :: t0
        integer :: i
        write (*, *) ''
        write (*, *) 'Averaging structure'

        lo_allocate(dr(ic%nx_internal + 1, sim%nt))
        dr = 0.0_flyt

        t0 = walltime()
        ! Align the simulation with the reference positions? Good idea? Not sure
        xavg = 0.0_flyt
        do i = 1, sim%nt
            call get_internal_irrep_from_supercell(ic, uc, ss, sim%r(:, :, i))
            xavg = xavg + ic%x_internal
            dr(1, i) = (i - 1)*sim%timestep*lo_time_au_to_fs
            dr(2:1 + ic%nx_internal, i) = ic%x_internal
        end do
        xavg = xavg/sim%nt

        flatr = matmul(ic%coeffM_internal, xavg)
        do i = 1, uc%na
            average_positions(:, i) = ic%reference_positions(:, i) + flatr((i - 1)*3 + 1:i*3)
        end do
        average_lattice = uc%latticevectors

        call lo_dump_gnuplot_2d_real(dr, 'outfile.internal')

        write (*, *) 'got pos', walltime() - t0
    end block init

    ! And return clean proper structures
    genclean: block
        real(flyt), dimension(:, :), allocatable :: dr
        real(flyt), dimension(3, 3) :: sslv, isslv
        real(flyt), dimension(3) :: v0
        integer, dimension(3) :: ssd
        integer :: i, j, l

        call nuc%generate(average_lattice, average_positions, uc%atomic_number, enhet=2)
        call nuc%classify('bravais')

        ! Then the supercell. Have to be careful so that it is identical to the generating
        ! supercell, so I will build it a little differently.
        sslv = matmul(nuc%latticevectors, ss%info%supercellmatrix)
        isslv = lo_invert3x3matrix(sslv)

        lo_allocate(dr(3, ss%na))
        do i = 1, ss%na
            v0 = matmul(nuc%latticevectors, ss%info%cellindex(:, i) - 1) + nuc%rcart(:, ss%info%index_in_unitcell(i))
            dr(:, i) = matmul(isslv, v0)
        end do
        call nss%generate(sslv, dr, ss%atomic_number, enhet=2)
        call nss%classify('supercell', nuc)
    end block genclean

end subroutine

! !> Calculate the average structure from a simulation, in a neat symmetrized way.
! subroutine avg_internal_structure(uc,ss,ic,sim,nuc,nss)
!     !> unitcell
!     type(lo_crystalstructure), intent(in) :: uc
!     !> supercell
!     type(lo_crystalstructure), intent(in) :: ss
!     !> irreducible cell
!     type(lo_irredcell), intent(inout) :: ic
!     !> simulation
!     type(lo_mdsim), intent(inout) :: sim
!     !> averaged unitcell
!     type(lo_crystalstructure), intent(out) :: nuc
!     !> averaged supercell
!     type(lo_crystalstructure), intent(out) :: nss
!
!     real(flyt), dimension(3,uc%na) :: average_positions
!     real(flyt), dimension(3,3) :: average_lattice
!
!     ! Get the average lattice
!     if ( allocated(sim%lattice) ) then
!     latt: block
!         real(flyt), dimension(:,:), allocatable :: coeffM
!         real(flyt), dimension(:), allocatable :: etaM
!         real(flyt), dimension(3,3) :: I3,trf1,etaguess,m0,lattice
!         real(flyt) :: f0,f1
!         integer :: t,i
!         ! Some transformations
!         I3=0.0_flyt
!         do i=1,3
!             I3(i,i)=1.0_flyt
!         enddo
!         trf1=lo_invert3x3matrix( real(ss%info%supercellmatrix,flyt) )
!         ! First guess for eta?
!         lo_allocate(coeffM(sim%nt*9,ic%nx_lattice))
!         lo_allocate(etaM(sim%nt*9))
!         coeffM=0.0_flyt
!         etaM=0.0_flyt
!
!         f0=0.0_flyt
!         do t=1,sim%nt
!             ! Current lattice
!             lattice=matmul(sim%lattice(:,:,t),trf1)
!             f0=f0+lo_determ(lattice)
!             ! Convert to a strain and store for solution
!             m0=matmul(lattice,uc%inv_latticevectors)-I3
!             coeffm( (t-1)*9+1:t*9,:)=ic%coeffM_lattice
!             etam( (t-1)*9+1:t*9)=lo_flattentensor(m0)
!         enddo
!         f0=f0/sim%nt
!         ! Solve
!         call lo_linear_least_squares(coeffm,etam,ic%x_lattice)
!         ! Multiply out to get the average
!         m0=lo_unflatten_2tensor(matmul(ic%coeffM_lattice,ic%x_lattice))+I3
!         average_lattice=matmul(m0,ic%reference_lattice)
!         f1=(f0/lo_determ(average_lattice))**(1.0_flyt/3.0_flyt)
!         average_lattice=average_lattice*f1
!
! write(*,*) 'latvec',lo_determ(average_lattice),f1,f0
! do i=1,3
!     write(*,*) average_lattice(:,i)*lo_bohr_to_A
! enddo
!     end block latt
!     endif
!
!     ! And the positions. Maybe a little trickier, maybe not.
!     posi: block
!         integer, parameter :: maxnavg=10
!         integer, parameter :: stride=2
!         real(flyt), dimension(:,:), allocatable :: ucref,ssref,ucpos,sspos,r0,r1
!         real(flyt), dimension(:,:), allocatable :: temp_avgpos_uc,temp_refpos_uc,refpos
!
!
!         real(flyt), dimension(:,:), allocatable :: irr_internal,dr,du
!         real(flyt), dimension(:,:), allocatable :: wA
!         real(flyt), dimension(:), allocatable :: WB
!         real(flyt), dimension(3,3) :: trf1,trf2,lattice,invlattice
!         real(flyt), dimension(3) :: v0,v1,v2,refshift
!         real(flyt) :: f0,f1,f2,f3
!         integer, dimension(:), allocatable :: dctr
!         integer :: nf,nt,tt,t,i,j,k,ii,jj
!         integer :: iter
!
!         ! Space for strange averages and things.
!         lo_allocate(sspos(3,ss%na))
!         lo_allocate(ssref(3,ss%na))
!         lo_allocate(ucref(3,uc%na))
!         lo_allocate(ucpos(3,uc%na))
!         lo_allocate(r0(3,uc%na))
!         lo_allocate(r1(3,uc%na))
!         ucref=uc%r
!         ssref=ss%r
!
!         ! Get unwrapped positions
!         do t=1,sim%nt
!         do i=1,sim%na
!             v0=lo_clean_fractional_coordinates(sim%r(:,i,t)-ssref(:,i)+0.5_flyt)-0.5_flyt
!             sim%r_npbc(:,i,t)=ssref(:,i)+v0
!         enddo
!         enddo
!
!         ! Raw average positions for the supercell
!         sspos=0.0_flyt
!         do t=1,sim%nt
!             do i=1,sim%na
!                 sspos(:,i)=sspos(:,i)+sim%r_npbc(:,i,t)
!             enddo
!         enddo
!         sspos=sspos/sim%nt
!
!         ! Downfolded average positions
!         f0=uc%volume/ss%volume
!         trf1=matmul(uc%inv_latticevectors,ss%latticevectors)
!         ucpos=0.0_flyt
!         do i=1,ss%na
!             j=ss%info%index_in_unitcell(i)
!             v0=ss%info%cellindex(:,i)-1
!             v0=matmul(trf1,sspos(:,i))-v0
!             ucpos(:,j)=ucpos(:,j)+v0*f0
!         enddo
!
!         ! Symmetrize these?
!         r0=ucpos
!         r1=0.0_flyt
!         do iter=1,2
!             do ii=1,uc%sym%n
!                 do i=1,uc%na
!                     j=uc%sym%op(ii)%fmap(i)
!                     r1(:,j)=lo_operate_on_vector(uc%sym%op(ii),r0(:,i),fractional=.true.)
!                 enddo
!                 r1=lo_clean_fractional_coordinates(r1-r0+0.5_flyt,1E-13_flyt)-0.5_flyt
!                 r0=lo_clean_fractional_coordinates(r0+0.5_flyt*r1,1E-13_flyt)
!             enddo
!             r0=lo_chop(r0,1E-13_flyt)
!         enddo
!
!         average_positions=r0
! !        do i=1,uc%na
! !            average_positions(:,i)=r0(:,i)-r0(:,1)
! !        enddo
!     end block posi
!
!     ! And return clean proper structures
!     genclean: block
!         real(flyt), dimension(:,:), allocatable :: dr
!         real(flyt), dimension(3,3) :: sslv,isslv
!         real(flyt), dimension(3) :: v0
!         integer, dimension(3) :: ssd
!         integer :: i,j,l
!
!         call nuc%generate(average_lattice,average_positions,uc%atomic_number,enhet=2)
!         call nuc%classify('bravais')
!
!         ! Then the supercell. Have to be careful so that it is identical to the generating
!         ! supercell, so I will build it a little differently.
!         sslv=matmul(nuc%latticevectors,ss%info%supercellmatrix)
!         isslv=lo_invert3x3matrix(sslv)
!
!         lo_allocate(dr(3,ss%na))
!         do i=1,ss%na
!             v0=matmul(nuc%latticevectors,ss%info%cellindex(:,i)-1)+nuc%rcart(:,ss%info%index_in_unitcell(i))
!             dr(:,i)=matmul(isslv,v0)
!         enddo
!         call nss%generate(sslv,dr,ss%atomic_number,enhet=2)
!         call nss%classify('supercell',nuc)
!         !
!     end block genclean
!
! end subroutine

!> Remove arbitrary precession of the lattice
subroutine remove_rotations(uc, ss, ic, sim)
    !> unitcell
    type(lo_crystalstructure), intent(in) :: uc
    !> supercell
    type(lo_crystalstructure), intent(in) :: ss
    !> irreducible cell
    type(lo_irredcell), intent(inout) :: ic
    !> simulation
    type(lo_mdsim), intent(inout) :: sim

    real(flyt), dimension(3, 3) :: rotation, etaguess, lattice, strain
    real(flyt), dimension(3, 3) :: m0, invlatt, trf1, trf2, I3
    real(flyt), dimension(3) :: v0, v1, com, v2
    real(flyt) :: invmass, invna
    integer :: i, t

    ! Some transformations
    I3 = 0.0_flyt
    do i = 1, 3
        I3(i, i) = 1.0_flyt
    end do
    trf1 = lo_invert3x3matrix(real(ss%info%supercellmatrix, flyt))
    trf2 = real(ss%info%supercellmatrix, flyt)
    ! First guess for eta?
    m0 = matmul(sim%lattice(:, :, 1), trf1)
    m0 = matmul(m0, invlatt)
    etaguess = m0

    ! First, temporarily convert everything to fractional coordinates.
    do t = 1, sim%nt
        lattice = matmul(sim%lattice(:, :, t), trf1)
        invlatt = lo_invert3x3matrix(lattice)
        ! Forces and velocities
        do i = 1, sim%na
            sim%f(:, i, t) = matmul(invlatt, sim%f(:, i, t))
            sim%v(:, i, t) = matmul(invlatt, sim%v(:, i, t))
        end do
        ! Set things to zero that have to be fixed later
        sim%u(:, :, t) = 0.0_flyt
        ! And the stress tensor. Hmm. I think it works like this.
        ! Actually not sure.
        m0 = sim%stat%stress(:, :, t)
        m0 = matmul(matmul(invlatt, m0), transpose(invlatt))
        sim%stat%stress(:, :, t) = m0
    end do

    ! First guess for eta?
    m0 = matmul(sim%lattice(:, :, 1), trf1)
    m0 = matmul(m0, invlatt)
    etaguess = m0
    ! Now find the rotaiton for the lattice, and rotate it
    do t = 1, sim%nt
        ! Get the transform for the lattice
        lattice = matmul(sim%lattice(:, :, t), trf1)
        call lo_find_rotation_that_makes_strain_diagonal(ic%reference_lattice, lattice, rotation, strain, etaguess)
        etaguess = strain
        ! Then transform it and store
        m0 = matmul(rotation, lattice)
        m0 = matmul(m0, trf2)
        sim%lattice(:, :, t) = m0
    end do

    ! Then convert everything back to Cartesian coordinates
    do t = 1, sim%nt
        lattice = matmul(sim%lattice(:, :, t), trf1)
        ! And the stress tensor. Hmm. I think it works like this.
        m0 = sim%stat%stress(:, :, t)
        m0 = matmul(matmul(lattice, m0), transpose(lattice))
        sim%stat%stress(:, :, t) = m0
        ! Forces, velocities and displacements. Not sure what makes more
        ! sense, to use some average or the instantaneous. Will use the
        ! instantaneous, I concluded for no good reason.
        do i = 1, sim%na
            sim%f(:, i, t) = matmul(lattice, sim%f(:, i, t))
            sim%v(:, i, t) = matmul(lattice, sim%v(:, i, t))
            v0 = lo_clean_fractional_coordinates(sim%r(:, i, t) - ss%r(:, i) + 0.5_flyt) - 0.5_flyt
            sim%r_npbc(:, i, t) = sim%r(:, i, t) + v0
            sim%u(:, i, t) = matmul(lattice, v0)
        end do
    end do

    invmass = 1.0_flyt/sum(ss%mass)
    invna = 1.0_flyt/real(ss%na, flyt)

    ! In the new coordinates, also remove any drift.
    do t = 1, sim%nt
        v0 = 0.0_flyt
        v1 = 0.0_flyt
        do i = 1, sim%na
            v2 = lo_clean_fractional_coordinates(sim%r(:, i, t) - ss%r(:, i) + 0.5_flyt) - 0.5_flyt
            v0 = v0 + v2 !*ss%mass(i)
            v1 = v1 + sim%f(:, i, t)
        end do

        v0 = v0*invna !mass
        v1 = v1*invna
        do i = 1, sim%na
            sim%r(:, i, t) = sim%r(:, i, t) - v0
            sim%f(:, i, t) = sim%f(:, i, t) - v1
            v2 = lo_clean_fractional_coordinates(sim%r(:, i, t) - ss%r(:, i) + 0.5_flyt) - 0.5_flyt
            sim%r_npbc(:, i, t) = sim%r(:, i, t) + v2
            sim%u(:, i, t) = matmul(sim%lattice(:, :, t), v2)
        end do
    end do
end subroutine

!> Calculate the average structure from a simulation, in a neat symmetrized way.
subroutine avg_structure(uc, ss, ic, sim, nuc, nss)
    !> unitcell
    type(lo_crystalstructure), intent(in) :: uc
    !> supercell
    type(lo_crystalstructure), intent(in) :: ss
    !> irreducible cell
    type(lo_irredcell), intent(inout) :: ic
    !> simulation
    type(lo_mdsim), intent(inout) :: sim
    !> averaged unitcell
    type(lo_crystalstructure), intent(out) :: nuc
    !> averaged supercell
    type(lo_crystalstructure), intent(out) :: nss

    real(flyt), dimension(3, uc%na) :: average_positions
    real(flyt), dimension(3, 3) :: average_lattice

    ! Get the average lattice
    if (allocated(sim%lattice)) then
        latt: block
            real(flyt), dimension(:, :), allocatable :: coeffM
            real(flyt), dimension(:), allocatable :: etaM
            real(flyt), dimension(3, 3) :: I3, trf1, etaguess, m0, lattice
            real(flyt) :: f0, f1
            integer :: t, i
            ! Some transformations
            I3 = 0.0_flyt
            do i = 1, 3
                I3(i, i) = 1.0_flyt
            end do
            trf1 = lo_invert3x3matrix(real(ss%info%supercellmatrix, flyt))
            ! First guess for eta?
            lo_allocate(coeffM(sim%nt*9, ic%nx_lattice))
            lo_allocate(etaM(sim%nt*9))
            coeffM = 0.0_flyt
            etaM = 0.0_flyt

            f0 = 0.0_flyt
            do t = 1, sim%nt
                ! Current lattice
                lattice = matmul(sim%lattice(:, :, t), trf1)
                f0 = f0 + lo_determ(lattice)
                ! Convert to a strain and store for solution
                m0 = matmul(lattice, uc%inv_latticevectors) - I3
                coeffm((t - 1)*9 + 1:t*9, :) = ic%coeffM_lattice
                etam((t - 1)*9 + 1:t*9) = lo_flattentensor(m0)
            end do
            f0 = f0/sim%nt
            ! Solve
            call lo_linear_least_squares(coeffm, etam, ic%x_lattice)
            ! Multiply out to get the average
            m0 = lo_unflatten_2tensor(matmul(ic%coeffM_lattice, ic%x_lattice)) + I3
            average_lattice = matmul(m0, ic%reference_lattice)
            f1 = (f0/lo_determ(average_lattice))**(1.0_flyt/3.0_flyt)
            average_lattice = average_lattice*f1

            write (*, *) 'latvec', lo_determ(average_lattice), f1, f0
            do i = 1, 3
                write (*, *) average_lattice(:, i)*lo_bohr_to_A
            end do
        end block latt
    else
        average_lattice = uc%latticevectors
    end if

    ! And the positions. Maybe a little trickier, maybe not.
    posi: block
        integer, parameter :: maxnavg = 10
        integer, parameter :: stride = 2
        real(flyt), dimension(:, :), allocatable :: ucref, ssref, ucpos, sspos, r0, r1
        real(flyt), dimension(:, :), allocatable :: temp_avgpos_uc, temp_refpos_uc, refpos

        real(flyt), dimension(:, :), allocatable :: irr_internal, dr, du
        real(flyt), dimension(:, :), allocatable :: wA
        real(flyt), dimension(:), allocatable :: WB
        real(flyt), dimension(3, 3) :: trf1, trf2, lattice, invlattice
        real(flyt), dimension(3) :: v0, v1, v2, refshift
        real(flyt) :: f0, f1, f2, f3
        integer, dimension(:), allocatable :: dctr
        integer :: nf, nt, tt, t, i, j, k, ii, jj
        integer :: iter

        ! Space for strange averages and things.
        lo_allocate(sspos(3, ss%na))
        lo_allocate(ssref(3, ss%na))
        lo_allocate(ucref(3, uc%na))
        lo_allocate(ucpos(3, uc%na))
        lo_allocate(r0(3, uc%na))
        lo_allocate(r1(3, uc%na))
        ucref = uc%r
        ssref = ss%r

        ! Get unwrapped positions
        do t = 1, sim%nt
        do i = 1, sim%na
            v0 = lo_clean_fractional_coordinates(sim%r(:, i, t) - ssref(:, i) + 0.5_flyt) - 0.5_flyt
            sim%r_npbc(:, i, t) = ssref(:, i) + v0
        end do
        end do
        write (*, *) sum(abs(sim%r_npbc))

        ! Raw average positions for the supercell
        sspos = 0.0_flyt
        do t = 1, sim%nt
            do i = 1, sim%na
                sspos(:, i) = sspos(:, i) + sim%r_npbc(:, i, t)
            end do
        end do
        sspos = sspos/sim%nt

        ! Downfolded average positions
        f0 = uc%volume/ss%volume
        trf1 = matmul(uc%inv_latticevectors, ss%latticevectors)
        ucpos = 0.0_flyt
        do i = 1, ss%na
            j = ss%info%index_in_unitcell(i)
            v0 = ss%info%cellindex(:, i) - 1
            v0 = matmul(trf1, sspos(:, i)) - v0
            ucpos(:, j) = ucpos(:, j) + v0*f0
        end do

        ! Symmetrize these?
        r0 = ucpos
        do iter = 1, 2
            do ii = 1, uc%sym%n
                r1 = 0.0_flyt
                do i = 1, uc%na
                    j = uc%sym%op(ii)%fmap(i)
                    r1(:, j) = lo_operate_on_vector(uc%sym%op(ii), r0(:, i), fractional=.true.)
                end do
                r1 = lo_clean_fractional_coordinates(r1 - r0 + 0.5_flyt, 1E-13_flyt) - 0.5_flyt
                r0 = lo_clean_fractional_coordinates(r0 + 0.5_flyt*r1, 1E-13_flyt)
            end do
            r0 = lo_chop(r0, 1E-13_flyt)
        end do

        average_positions = r0
!        do i=1,uc%na
!            average_positions(:,i)=r0(:,i)-r0(:,1)
!        enddo
    end block posi

    ! And return clean proper structures
    genclean: block
        real(flyt), dimension(:, :), allocatable :: dr
        real(flyt), dimension(3, 3) :: sslv, isslv
        real(flyt), dimension(3) :: v0
        integer, dimension(3) :: ssd
        integer :: i, j, l

        call nuc%generate(average_lattice, average_positions, uc%atomic_number, enhet=2)
        call nuc%classify('bravais')

        ! Then the supercell. Have to be careful so that it is identical to the generating
        ! supercell, so I will build it a little differently.
        sslv = matmul(nuc%latticevectors, ss%info%supercellmatrix)
        isslv = lo_invert3x3matrix(sslv)

        lo_allocate(dr(3, ss%na))
        do i = 1, ss%na
            v0 = matmul(nuc%latticevectors, ss%info%cellindex(:, i) - 1) + nuc%rcart(:, ss%info%index_in_unitcell(i))
            dr(:, i) = matmul(isslv, v0)
        end do
        call nss%generate(sslv, dr, ss%atomic_number, enhet=2)
        call nss%classify('supercell', nuc)
        !
    end block genclean

end subroutine

!!> Get the irreducible trajectory
!subroutine irreducible_trajectory(uc,ss,ic,sim)
!    !> unitcell
!    type(lo_crystalstructure), intent(in) :: uc
!    !> supercell
!    type(lo_crystalstructure), intent(in) :: ss
!    !> irreducible cell
!    type(lo_irredcell), intent(inout) :: ic
!    !> simulation
!    type(lo_mdsim), intent(inout) :: sim
!
!    real(flyt), dimension(3,3) :: avg_lattice
!
!
!
!
!    ! First the lattice
!    latt: block
!        real(flyt), dimension(:,:), allocatable :: irr_lattice
!        real(flyt), dimension(9) :: fv0
!        real(flyt), dimension(3,3) :: m0,trf1,I3
!        integer :: i
!real(flyt) :: f0,f1,f2
!
!        I3=0.0_flyt
!        do i=1,3
!            I3(i,i)=1.0_flyt
!        enddo
!
!        ! Transformation matrix to unitcell
!        trf1=lo_invert3x3matrix( real(ss%info%supercellmatrix,flyt) )
!        lo_allocate(irr_lattice(ic%nx_lattice,sim%nt))
!        irr_lattice=0.0_flyt
!
!        f0=0.0_flyt
!        do i=1,sim%nt
!            ! Get the unitcell strain. Always useful. First convert to unitcell lattice vectors
!            m0=matmul(trf1,sim%lattice(:,:,i))
!            call ic%strain_to_irreducible(m0)
!            irr_lattice(:,i)=ic%x_lattice
!            f0=f0+lo_determ(sim%lattice(:,:,i))/sim%na/sim%nt
!        enddo
!        write(*,*) 'avg vol/A',f0
!        do i=1,ic%nx_lattice
!            ic%x_lattice(i)=lo_mean(irr_lattice(i,:))
!        enddo
!        fv0=matmul(ic%coeffM_lattice,ic%x_lattice)
!        avg_lattice=matmul(ic%reference_lattice,I3+lo_unflatten_2tensor(fv0))
!
!f0=lo_determ(avg_lattice)/uc%na
!write(*,*) 'avg vol/A',f0
!!
!!fv0=matmul(ic%coeffM_lattice,[1.0_flyt,0.5_flyt])
!!avg_lattice=lo_unflatten_2tensor(fv0)
!!
!m0=lo_unflatten_2tensor(fv0)
!
!do i=1,3
!    write(*,*)  m0(:,i) !avg_lattice(:,i)*lo_bohr_to_A
!enddo
!do i=1,3
!    write(*,*) ic%reference_lattice(:,i)
!enddo
!do i=1,3
!    write(*,*) avg_lattice(:,i)/norm2(avg_lattice(:,i))
!enddo
!stop
!
!!
!
!    end block latt
!
!    ! Think about center of mass, drift and such.
!    centerofmass: block
!        real(flyt), dimension(3) :: v0,v1,v2
!        integer :: i,j,k,l
!
!        ! First align the centers of mass with each other.
!    !    do i=1,sim%nt
!    !        v0=0.0_flyt
!    !        do j=1,sim%na
!    !            v1=sim%r_npbc(:,j,i)-ss%r(:,j)
!    !            v0=v0+v1*ss%mass(j)
!    !        enddo
!    !        v0=v0/sum(ss%mass)
!    !        do j=1,sim%na
!    !            sim%r(:,j,i)=sim%r(:,j,i)-v0
!    !            sim%r_npbc(:,j,i)=sim%r_npbc(:,j,i)-v0
!    !        enddo
!    !    enddo
!    !    ! Then make sure there is no rigid shift of everything.
!    !    v0=0.0_flyt
!    !    do i=1,sim%nt
!    !    do j=1,sim%na
!    !        v1=lo_clean_fractional_coordinates(sim%r_npbc(:,j,i)-ss%r(:,j)+0.5_flyt)-0.5_flyt
!    !        v0=v0+v1
!    !    enddo
!    !    enddo
!    !    v0=v0/sim%nt/sim%na
!    !    do i=1,sim%nt
!    !        do j=1,sim%na
!    !            sim%r(:,j,i)=sim%r(:,j,i)-v0
!    !            sim%r_npbc(:,j,i)=sim%r_npbc(:,j,i)-v0
!    !        enddo
!    !    enddo
!    end block centerofmass
!
!    ! Then the internal positions. Hmm.
!    posi: block
!        integer, parameter :: maxnavg=10
!        real(flyt), dimension(:,:), allocatable :: irr_internal
!        real(flyt), dimension(:,:), allocatable :: big_coeff,big_ref,wA,WB
!        real(flyt), dimension(:), allocatable :: big_pos
!        real(flyt), dimension(3,3) :: trf1
!        real(flyt), dimension(3) :: v0,v1
!        real(flyt), dimension(:), allocatable :: invctr
!        integer :: i,j,k,l,ii
!        integer :: navg,nf
!
!        ! How many bins for the averaging?
!        navg=min(sim%nt,maxnavg)
!        nf=sim%na*3
!
!        ! First get the reference positions
!        lo_allocate(big_ref(3,ss%na))
!        lo_allocate(big_pos(3*ss%na))
!        lo_allocate(big_coeff(3*ss%na,ic%nx_internal))
!        lo_allocate(irr_internal(ic%nx_internal,sim%nt))
!        irr_internal=0.0_flyt
!        big_ref=0.0_flyt
!
!        ! Get coefficient matrix and reference positions
!        big_coeff=0.0_flyt
!        do i=1,ss%na
!            j=ss%info%index_in_unitcell(i)
!            v0=ic%reference_positions(:,j)
!            v0=v0+(ss%info%cellindex(:,i)-1)
!            big_ref( :,i )=v0
!            big_coeff( (i-1)*3+1:i*3,: )=ic%coeffM_internal( (j-1)*3+1:j*3,:)
!        enddo
!        ! Expand the big matrix
!        lo_allocate(wA(nf*navg,ic%nx_internal))
!        lo_allocate(wB(nf*navg,1))
!        wA=0.0_flyt
!        wB=0.0_flyt
!        do i=1,navg
!            wA( (i-1)*nf+1:i*nf,: )=big_coeff
!        enddo
!
!        ! Get the averaging factor
!        lo_allocate(invctr(navg))
!        invctr=0.0_flyt
!        do i=1,sim%nt
!            j=mod(i,navg)+1
!            invctr(j)=invctr(j)+1.0_flyt
!        enddo
!        invctr=1.0_flyt/invctr
!
!        trf1=real(ss%info%supercellmatrix,flyt)
!        sim%r_npbc(:,:,1)=ss%r
!        do i=1,sim%nt
!
!            do j=1,sim%na
!                v0=matmul(trf1,sim%r_npbc(:,j,i))-big_ref(:,j)
!                big_pos( (j-1)*3+1:j*3 )=v0
!            enddo
!            !big_pos=lo_clean_fractional_coordinates(big_pos+0.5_flyt)-0.5_flyt
!            irr_internal(:,i)=matmul(transpose(big_coeff),big_pos)
!write(*,*) i,irr_internal(:,i)
!        enddo
!
!        do i=1,sim%na
!        do j=1,3
!            big_ref(j,i)=lo_mean(sim%r_npbc(j,i,:))
!        enddo
!        enddo
!        big_ref=lo_clean_fractional_coordinates(big_ref)
!write(*,*) 'Average:'
!        do i=1,sim%na
!            write(*,*) i,lo_clean_fractional_coordinates(big_ref(:,i)-ss%r(:,i)+0.5_flyt)-0.5_flyt
!        enddo
!
!    end block posi
!
!end subroutine

!> create the irreducible representation of the lattice
subroutine generate(ic, p, verbosity, reference_positions)
    !> irreducible representation
    class(lo_irredcell), intent(out) :: ic
    !> normal crystal structure
    type(lo_crystalstructure), intent(inout) :: p
    !> talk?
    integer, intent(in) :: verbosity
    !> Not the normal reference positions?
    real(flyt), dimension(:, :), intent(in), optional :: reference_positions

    ! decide on settings and such
    init: block
        integer :: i
        if (verbosity .gt. 0) then
            write (*, *) ''
            write (*, *) 'BUILDING IRREDUCIBLE REPRESENTATION OF THE LATTICE'
        end if

        ! We need all the symmetry things
        if (p%info%havespacegroup .eqv. .false.) then
            call p%classify('wedge', timereversal=.true.)
        end if

        ! Store the reference positions
        lo_allocate(ic%reference_positions(3, p%na))
        if (present(reference_positions)) then
            ic%reference_positions = reference_positions
        else
            ic%reference_positions = p%r
        end if
        ic%reference_lattice = p%latticevectors

    end block init

    ! Start with the lattice vectors, or more precisely the symmetry of the strain tensor.
    lattice: block
        real(flyt), dimension(9, 9) :: invarM, I9, T, rotM
        real(flyt), dimension(3, 3) :: m0
        real(flyt), dimension(9) :: dv
        integer :: i, j

        ! Transposition operator
        T(:, 1) = [1, 0, 0, 0, 0, 0, 0, 0, 0]
        T(:, 2) = [0, 0, 0, 1, 0, 0, 0, 0, 0]
        T(:, 3) = [0, 0, 0, 0, 0, 0, 1, 0, 0]
        T(:, 4) = [0, 1, 0, 0, 0, 0, 0, 0, 0]
        T(:, 5) = [0, 0, 0, 0, 1, 0, 0, 0, 0]
        T(:, 6) = [0, 0, 0, 0, 0, 0, 0, 1, 0]
        T(:, 7) = [0, 0, 1, 0, 0, 0, 0, 0, 0]
        T(:, 8) = [0, 0, 0, 0, 0, 1, 0, 0, 0]
        T(:, 9) = [0, 0, 0, 0, 0, 0, 0, 0, 1]
        ! identity matrix
        call lo_identitymatrix(I9)
        ! now add all the operations of the lattice
        invarM = 0.0_flyt
        do i = 1, p%sym%n
            rotM = lo_expandoperation_pair(p%sym%op(i)%m)
            invarM = invarM + rotM - I9
        end do
        ! and transposition symmetry
        invarM = invarM + T - I9
        ! and project out the nullspace
        call lo_nullspace_coefficient_matrix(invarM, ic%coeffM_lattice, ic%nx_lattice, tolerance=lo_sqtol)

        ! It's weird if it's not at least 1, or at most 6, since the number of zero singular values
        ! should be the number of degrees of freedom.
        if (ic%nx_lattice .lt. 1 .or. ic%nx_lattice .gt. 6) then
            call lo_stop_gracefully(['Strain tensor must have at least 1 and at most 6 degrees of freedom'], &
                                    lo_exitcode_symmetry, __FILE__, __LINE__)
        end if
        ! and a little space for the solution
        lo_allocate(ic%x_lattice(ic%nx_lattice))
        ic%x_lattice = 0.0_flyt

        if (verbosity .gt. 0) then
            write (*, *) '... found ', tochar(ic%nx_lattice), ' degrees of freedom for the lattice:'
            do i = 1, ic%nx_lattice
                write (*, *) 'x'//tochar(i)//':'
                dv = 0.0_flyt
                dv(i) = 1.0_flyt
                m0 = lo_unflatten_2tensor(matmul(ic%coeffM_lattice, dv(1:ic%nx_lattice)))
                do j = 1, 3
                    write (*, "(3(2X,F16.9))") m0(:, j)
                end do
            end do
        end if

    end block lattice

    ! Now do the same thing, but for the internal positions. And in this case it's the symmetry of the
    ! forces we are interested in.
    bigmatrix: block
        real(flyt), dimension(:, :), allocatable :: invarM, rotM, IM, dr
        real(flyt), dimension(:), allocatable :: dv
        integer :: i, j, o, n, u
        !
        n = p%na*3
        ! Some space for stuff
        lo_allocate(rotM(n, n))
        lo_allocate(IM(n, n))
        lo_allocate(invarM(n, n))
        call lo_identitymatrix(IM)

        ! Build a matrix representation how all the internal positions transform under an
        ! operation. I use large block matrices for this.
        invarM = 0.0_flyt
        do o = 1, p%sym%n
            rotM = 0.0_flyt
            do i = 1, p%na
                ! the semimagic %fmap here keeps track of which atom is transform to which atom
                ! under each operation. I can mimic this switcharoo by putting the small 3x3
                ! rotation matrix in the right place in the big matrix. That way I don't have to
                ! deal with translations explicitly, and life becomes much easier.
                j = p%sym%op(o)%fmap(i)
                !rotM((j-1)*3+1:j*3,(i-1)*3+1:3*i)=transpose(p%sym%op(o)%fm)
                ! Mathematica confused me again, this is how it's supposed to be. Hmm.
                rotM((i - 1)*3 + 1:i*3, (j - 1)*3 + 1:3*j) = transpose(p%sym%op(o)%fm)
            end do
            invarM = invarM + rotM - IM
        end do

        ! Then add the zero-sum constraint? You should convince yourself that this really is the
        ! matrix representation of adding all forces together.
        rotM = 0.0_flyt
        do i = 1, p%na
            do j = 1, 3
                rotM((i - 1)*3 + j, j) = 1.0_flyt
            end do
        end do
        invarM = invarM + rotM

        ! In an ideal world I should probably add rotational invariance as well, but I don't think
        ! it matters, and I also don't know how to do it. Also worth noting is that the symmetry is
        ! done in fractional coordinates, and all shifts from starting positions are also in fractional
        ! coordinates. Otherwise it gets real wonky when the volume changes.

        ! Anyway, now we can project out the nullspace of this matrix
        call lo_nullspace_coefficient_matrix(invarM, ic%coeffM_internal, ic%nx_internal, tolerance=lo_sqtol)
        ! and some space for a solution
        if (ic%nx_internal .gt. 0) then
            lo_allocate(ic%x_internal(ic%nx_internal))
            ic%x_internal = 0.0_flyt
        end if

        if (verbosity .gt. 0) then
            write (*, *) '... found ', tochar(ic%nx_internal), ' internal degrees of freedom'
            lo_allocate(dv(p%na*3))
            lo_allocate(dr(3, p%na))
            do i = 1, ic%nx_internal
                write (*, *) 'x'//tochar(i)//': (fractional, Cartesian)'
                dv = 0.0_flyt
                dv(i) = 1.0_flyt
                dv = matmul(ic%coeffM_internal, dv(1:ic%nx_internal))
                do j = 1, p%na
                    dr(:, j) = dv((j - 1)*3 + 1:j*3)
                    write (*, "(6(2X,F15.9))") dr(:, j), matmul(p%latticevectors, dr(:, j))
                end do
            end do
            write (*, *) ''
            write (*, *) 'Reference positions:'
            do j = 1, p%na
                dr(:, j) = dv((j - 1)*3 + 1:j*3)
                write (*, "(3(2X,F15.9))") ic%reference_positions(:, j)
            end do
        end if
    end block bigmatrix

end subroutine

! !> Turn the raw forces and stress tensor from AIMS into something clean and symmetric. Not used for anything at the moment. But could come handy on a rainy day.
! subroutine symmetrize_stress_and_forces(ic,p,orig_stress,orig_forces,clean_stress,clean_forces)
!     !> irreducible representation
!     class(lo_irredcell), intent(inout) :: ic
!     !> crystal structure
!     type(lo_crystalstructure), intent(in) :: p
!     !> noisy stress tensor
!     real(flyt), dimension(3,3), intent(in) :: orig_stress
!     !> noisy forces
!     real(flyt), dimension(:,:), intent(in) :: orig_forces
!     !> nice stress tensor
!     real(flyt), dimension(3,3), intent(out) :: clean_stress
!     !> nice forces
!     real(flyt), dimension(:,:), intent(out) :: clean_forces
!
!     ! least square force the stress into the irreducible form
!     fixlattice: block
!         real(flyt), dimension(3,3) :: m0
!         real(flyt), dimension(9,ic%nx_lattice) :: wM
!         real(flyt), dimension(9,1) :: wV
!         real(flyt) :: f0
!         integer :: i
!
!         ! Solve for the irreducible components
!         wM=ic%coeffM_lattice
!         wV(:,1)=lo_flatten_2tensor(orig_stress)
!         call lo_dgels(wM,wV,info=lo_status)
!         if ( lo_status .ne. 0 ) then
!             write(*,*) 'Lapack failed for the lattice'
!             stop
!         endif
!         ic%x_lattice=wV(1:ic%nx_lattice,1)
!         ! Store the symetrized stress
!         clean_stress=lo_chop( lo_unflatten_2tensor(matmul(ic%coeffM_lattice,ic%x_lattice)), 1E-12_flyt)
!     end block fixlattice
!
!     ! Same thing for the forces
!     if ( ic%nx_internal .gt. 0 ) then
!     fixforces: block
!         real(flyt), dimension(:,:), allocatable :: mA,mC,mB
!         real(flyt), dimension(:), allocatable :: vB,vD
!         real(flyt), dimension(3) :: v0
!         integer :: i,j,l,n
!
!         n=size(ic%coeffM_internal,1)
!         if ( size(orig_forces,1)*size(orig_forces,2) .ne. n ) then
!             write(*,*) 'Wrong number of forces'
!             stop
!         endif
!
!         ! normal least squares
!         lo_allocate(mA(n,ic%nx_internal))
!         lo_allocate(mB(n,1))
!         l=0
!         do i=1,size(orig_forces,2)
!         do j=1,size(orig_forces,1)
!             l=l+1
!             mB(l,1)=orig_forces(j,i)
!         enddo
!         enddo
!         call lo_dgels(mA,mB,info=lo_status)
!         if ( lo_status .ne. 0 ) then
!             write(*,*) 'Lapack failed for forces'
!             stop
!         endif
!         ic%x_internal=mB(1:ic%nx_internal,1)
!
!         ! Store the symmetrized forces
!         if ( .not.allocated(vB) ) allocate(vB(n))
!         vB=matmul(ic%coeffM_internal,ic%x_internal)
!         do i=1,p%na
!             v0=vb( (i-1)*3+1:i*3 )
!             v0=matmul( p%latticevectors, v0 )
!             clean_forces(:,i)=lo_chop( v0, 1E-12_flyt )
!         enddo
!     end block fixforces
!     else
!         ! In case there are no internal degrees of freedom.
!         clean_forces=0.0_flyt
!     endif
! end subroutine

end module
