program extract_forceconstants
!!{!src/extract_forceconstants/manual.md!}!
use konstanter, only: r8, i8, lo_tol, lo_bohr_to_A, lo_pressure_HartreeBohr_to_GPa, lo_Hartree_to_eV, lo_eV_to_Hartree, &
                      lo_exitcode_param, lo_force_hartreebohr_to_eVa, lo_iou, lo_forceconstant_1st_HartreeBohr_to_eVA
use gottochblandat, only: tochar, walltime, open_file, lo_chop, lo_does_file_exist, lo_trueNtimes, lo_mean, lo_stddev, &
                          lo_outerproduct, lo_trace, lo_frobnorm, lo_linear_least_squares
use mpi_wrappers, only: lo_mpi_helper, lo_stop_gracefully
use type_crystalstructure, only: lo_crystalstructure
use type_forceconstant_firstorder, only: lo_forceconstant_firstorder
use type_forceconstant_secondorder, only: lo_forceconstant_secondorder
use type_forceconstant_thirdorder, only: lo_forceconstant_thirdorder
use type_forceconstant_fourthorder, only: lo_forceconstant_fourthorder
use lo_dielectric_interaction, only: lo_dielectric_tensor
use lo_symmetry_of_interactions, only: lo_interaction_tensors
use type_forcemap, only: lo_forcemap, lo_secondorder_rot_herm_huang
use type_mdsim, only: lo_mdsim
use type_qpointmesh, only: lo_qpoint
use lo_memtracker, only: lo_mem_helper

use options, only: lo_opts
use ifc_solvers, only: lo_solve_for_irreducible_ifc, lo_solve_for_irreducible_ifc_fastugly

implicit none
type(lo_opts) :: opts
type(lo_crystalstructure) :: uc, ss
type(lo_mdsim) :: sim
type(lo_forcemap) :: map
type(lo_forceconstant_firstorder) :: fc1
type(lo_forceconstant_secondorder) :: fc2
type(lo_forceconstant_thirdorder) :: fc3
type(lo_forceconstant_fourthorder) :: fc4
type(lo_dielectric_tensor) :: di
!type(lo_jij_secondorder) :: jij
type(lo_mpi_helper) :: mw
type(lo_mem_helper) :: mem
real(r8) :: t0

! To start, I need at least the unitcell
init: block
    call mw%init()
    call mem%init()
    t0 = walltime()
    call opts%parse()
    if (mw%talk) write (lo_iou, *) 'READ STRUCTURE AND SETUP CUTOFFS'
    if (mw%talk .eqv. .false.) opts%verbosity = -100
    if (mw%talk) write (lo_iou, *) '... reading unitcell'
    call uc%readfromfile('infile.ucposcar', verbosity=opts%verbosity)
    call uc%classify('wedge', timereversal=.true.)
end block init

! We need a forcemap to continue.
getfrcmap: block
    type(lo_interaction_tensors) :: slt

    ! Either read it from file or create it
    if (opts%readforcemap) then
        ! I might need the supercell anyway
        if (opts%readirreducible .eqv. .false.) then
            call ss%readfromfile('infile.ssposcar')
            call ss%classify('supercell', uc)
        end if
        call map%read_from_hdf5(uc, 'infile.forcemap.hdf5', opts%verbosity)
        ! Update constraints that depend on structure
        !call lo_secondorder_rot_herm_huang(map, uc, map%constraints%eq2, map%constraints%neq2, &
        !                                   opts%rotationalconstraints, opts%huanginvariance, opts%hermitian)
    else
        ! Now we have to calculate the whole thing.
        call uc%classify('wedge', timereversal=.true.)
        call ss%readfromfile('infile.ssposcar')
        if (mw%talk) write (*, '(1X,A,1X,F12.5)') '... min cutoff: ', ss%mincutoff()*lo_bohr_to_A
        if (mw%talk) write (*, '(1X,A,1X,F12.5)') '... max cutoff: ', ss%maxcutoff()*lo_bohr_to_A
        if (mw%talk) write (*, '(1X,A,1X,F12.5)') '--> rc2 cutoff: ', opts%cutoff2*lo_bohr_to_A
        if (mw%talk) write (*, '(1X,A,1X,F12.5)') '--> rc3 cutoff: ', opts%cutoff3*lo_bohr_to_A
        if (mw%talk) write (*, '(1X,A,1X,F12.5)') '--> rc4 cutoff: ', opts%cutoff4*lo_bohr_to_A
        ! Die early if there is no infile.lotosplitting
        if (opts%polar) then
        if (lo_does_file_exist('infile.lotosplitting') .eqv. .false.) then
            call lo_stop_gracefully(['You need to provide "infile.lotosplitting" for polar corrections'], &
                                    lo_exitcode_param, __FILE__, __LINE__, mw%comm)
        end if
        end if
        ! Get all the symmetries
        call slt%generate(uc, ss, opts%cutoff2, opts%cutoff3, opts%cutoff4, opts%polar, mw, mem, opts%verbosity + 1, &
                          transposition=opts%transposition, spacegroup=opts%spacegroup, &
                          wzdim=opts%wzdimensions, nj2=opts%njump2, nj3=opts%njump3, nj4=opts%njump4, &
                          firstorder=opts%firstorder, magcutoff2=opts%magcutoff2, magsinglet=opts%magnetic_onsite, &
                          dielcutoff2=opts%dielcutoff2, dielcutoff3=opts%dielcutoff3)
        ! And create the map
        call map%generate(uc, ss, polarcorrectiontype=opts%polarcorrectiontype, &
                          st=slt, mw=mw, mem=mem, verbosity=opts%verbosity, devmode=opts%devmode)

        if (opts%printforcemap .and. mw%talk) then
            call map%write_to_hdf5('outfile.forcemap.hdf5', opts%verbosity)
        end if
    end if

end block getfrcmap

idiot: block
    integer :: i
    if (opts%fake_dielectric .and. map%polar .gt. 0) then
        ! This mean just dump some fake things to have something to work with, for debugging.
        if (map%xuc%nx_Z_singlet .gt. 0) map%xuc%x_Z_singlet = 1.0_r8
        if (map%xuc%nx_Z_pair .gt. 0) map%xuc%x_Z_pair = 1.0_r8
        if (map%xuc%nx_Z_triplet .gt. 0) map%xuc%x_Z_triplet = 1.0_r8
        if (map%xuc%nx_eps_global .gt. 0) map%xuc%x_eps_global = 1.0_r8
        if (map%xuc%nx_eps_singlet .gt. 0) map%xuc%x_eps_singlet = 1.0_r8
        if (map%xuc%nx_eps_pair .gt. 0) map%xuc%x_eps_pair = 1.0_r8

        if (map%xuc%nx_eps_singlet .gt. 0) then
            map%xuc%x_eps_singlet = 1.0_r8
            do i = 1, map%xuc%nx_eps_singlet
                map%xuc%x_eps_singlet(i) = i !1.0_r8
            end do
        end if
        call map%get_dielectric_tensors(uc, di)
        if (mw%talk) call di%writetofile(uc, 'outfile.dielectric_interaction')
    end if
end block idiot

! With the forcemap, we can get the irreducible representation.
if (opts%readirreducible) then
    irrfromfile: block
        integer :: i, u

        if (map%have_fc_singlet) then
            u = open_file('in', 'infile.irrifc_firstorder')
            do i = 1, map%xuc%nx_fc_singlet
                read (u, *) map%xuc%x_fc_singlet(i)
            end do
            close (u)
            call map%get_firstorder_forceconstant(uc, fc1)
        end if
        if (map%have_fc_pair) then
            u = open_file('in', 'infile.irrifc_secondorder')
            do i = 1, map%xuc%nx_fc_pair
                read (u, *) map%xuc%x_fc_pair(i)
            end do
            close (u)
        end if
        if (map%have_fc_triplet) then
            u = open_file('in', 'infile.irrifc_thirdorder')
            do i = 1, map%xuc%nx_fc_triplet
                read (u, *) map%xuc%x_fc_triplet(i)
            end do
            close (u)
        end if
        if (map%have_fc_quartet) then
            u = open_file('in', 'infile.irrifc_fourthorder')
            do i = 1, map%xuc%nx_fc_quartet
                read (u, *) map%xuc%x_fc_quartet(i)
            end do
            close (u)
        end if
        if (map%polar .eq. 1) then
            u = open_file('in', 'infile.irrtheta_loto')
            do i = 1, map%xuc%nx_eps_global
                read (u, *) map%xuc%x_eps_global(i)
            end do
            do i = 1, map%xuc%nx_Z_singlet
                read (u, *) map%xuc%x_Z_singlet(i)
            end do
            close (u)
        end if
        ! if ( map%magnetic_pair_interactions ) then
        !     write(*,*) 'FIXME READ IRREDUCIBLE MAGNETIC PAIRS'
        !     stop
        ! endif
        if (mw%talk) write (lo_iou, *) '... enforce the constraints'
        call map%enforce_constraints()
    end block irrfromfile
else
    getirr: block
        logical :: readmag, readdiel
        integer :: i, u

        ! Now we have to calculate the whole thing, first we need positions and forces.
        if (lo_does_file_exist('infile.sim.hdf5')) then
            call sim%read_from_hdf5('infile.sim.hdf5', verbosity=opts%verbosity + 2, stride=opts%stride)
        else
            ! if ( map%magnetic_pair_interactions ) then
            !     readmag=.true.
            ! else
            !     readmag=.false.
            ! endif
            readmag = .false.
            if (map%polar .gt. 1) then
                if (opts%dielcutoff2 .lt. 0.0_r8) then
                    readdiel = .false.
                else
                    readdiel = .true.
                end if
            else
                readdiel = .false.
            end if
            ! Heuristics are tricky. Hmmm.
            call sim%read_from_file(verbosity=opts%verbosity + 2, stride=opts%stride, &
                                    magnetic=readmag, dielectric=readdiel, nrand=opts%nrandom, mw=mw)
        end if

        ! Report how overdetermined the system is
        overdetermined: block
            real(r8) :: n_equations
            integer(i8) :: nx_fc_total
            character(*), parameter :: fmt1 = '(5X,A,I12)'
            character(*), parameter :: fmt = '(5X,A,I12,3X,A,1X,F12.1)'

            nx_fc_total = 0_i8
            n_equations = 3.0_r8*sim%na*sim%nt

            if (mw%talk) then
                write (lo_iou, *) ''
                write (lo_iou, *) 'REPORT GRADE OF OVERDETERMINATION'
                write (lo_iou, fmt1) 'total number of equations: ', int(n_equations)
                if (map%have_fc_singlet) then
                    nx_fc_total = nx_fc_total + map%xuc%nx_fc_singlet
                    write (lo_iou, fmt) ' number of FC up 1. order: ', &
                        nx_fc_total, '--> overdetermination ratio:', n_equations/nx_fc_total

                end if
                if (map%have_fc_pair) then
                    nx_fc_total = nx_fc_total + map%xuc%nx_fc_pair
                    write (lo_iou, fmt) ' number of FC up 2. order: ', &
                        nx_fc_total, '--> overdetermination ratio:', n_equations/nx_fc_total
                end if
                if (map%have_fc_triplet) then
                    nx_fc_total = nx_fc_total + map%xuc%nx_fc_triplet
                    write (lo_iou, fmt) ' number of FC up 3. order: ', &
                        nx_fc_total, '--> overdetermination ratio:', n_equations/nx_fc_total

                end if
                if (map%have_fc_quartet) then
                    nx_fc_total = nx_fc_total + map%xuc%nx_fc_quartet
                    write (lo_iou, fmt) ' number of FC up 4. order: ', &
                        nx_fc_total, '--> overdetermination ratio:', n_equations/nx_fc_total
                end if
                ! write (lo_iou, *) ''
            end if
        end block overdetermined

        ! Actual solution
        call lo_solve_for_irreducible_ifc(map, sim, uc, ss, opts%solver, mw, mem, opts%verbosity + 2, &
                                          lotofilename='infile.lotosplitting', &
                                          !dev1=dev1,dev2=dev2,dev3=dev3,dev4=dev4,&
                                          temperature=opts%temperature, &
                                          max_displacement=opts%max_displacement, &
                                          noreference=opts%noreference)

        !call lo_solve_for_irreducible_ifc_fastugly(map,sim,mw,mem,opts%verbosity+1)

        ! Solve for magnetic Jij
        ! if ( map%magnetic_pair_interactions ) then
        !     write(*,*) 'FIXME MAGNETIC'
        !     stop
        !     !call map%solve_for_jij(sim,uc,ss,verbosity=2)
        ! endif

        ! And store the irreducible representation
        if (map%have_fc_singlet .and. mw%talk) then
            u = open_file('out', 'outfile.irrifc_firstorder')
            do i = 1, map%xuc%nx_fc_singlet
                write (u, *) map%xuc%x_fc_singlet(i) !,dev1(i)
            end do
            close (u)
        end if
        if (map%have_fc_pair .and. mw%talk) then
            u = open_file('out', 'outfile.irrifc_secondorder')
            do i = 1, map%xuc%nx_fc_pair
                write (u, *) map%xuc%x_fc_pair(i) !,dev2(i)
            end do
            close (u)
        end if
        if (map%have_fc_triplet .and. mw%talk) then
            u = open_file('out', 'outfile.irrifc_thirdorder')
            do i = 1, map%xuc%nx_fc_triplet
                write (u, *) map%xuc%x_fc_triplet(i) !,dev3(i)
            end do
            close (u)
        end if
        if (map%have_fc_quartet .and. mw%talk) then
            u = open_file('out', 'outfile.irrifc_fourthorder')
            do i = 1, map%xuc%nx_fc_quartet
                write (u, *) map%xuc%x_fc_quartet(i) !,dev4(i)
            end do
            close (u)
        end if
        if (map%polar .gt. 0 .and. mw%talk) then
            u = open_file('out', 'outfile.irrtheta_loto')
            do i = 1, map%xuc%nx_eps_global
                write (u, *) map%xuc%x_eps_global(i)
            end do
            do i = 1, map%xuc%nx_Z_singlet
                write (u, *) map%xuc%x_Z_singlet(i)
            end do
            close (u)
        end if
        ! if ( map%magnetic_pair_interactions .and. mw%talk ) then
        !     u=open_file('out','outfile.irrtheta_magsecondorder')
        !         do i=1,map%ntheta_magpair_jij
        !             write(u,*) map%theta_magpair_jij(i)
        !         enddo
        !         do i=1,map%ntheta_magpair_tij
        !             write(u,*) map%theta_magpair_tij(i)
        !         enddo
        !     close(u)
        ! endif
    end block getirr
end if

! Now we can get the actual forceconstants from the forcemap and the irreducible representation
getfc: block
    integer :: i
    real(r8), dimension(3) :: force_norm

    if (map%have_fc_singlet) then
        call map%get_firstorder_forceconstant(uc, fc1)
        if (mw%talk) then
            force_norm = 0.0_r8
            write (*, *) 'residual forces per atom (eV/A):'
            do i = 1, uc%na
                write (*, *) i, fc1%atom(i)%m*lo_forceconstant_1st_HartreeBohr_to_eVA
                force_norm = force_norm + (fc1%atom(i)%m*lo_forceconstant_1st_HartreeBohr_to_eVA)**2
            end do
            write (*, '(90("-"))')
            write (*, *) '      RMSE:', sqrt(force_norm/uc%na)
            write (*, '(A,E12.3,A)') ' RMSE total:', sqrt(sum(force_norm)/uc%na), ' (eV/AA)'
            write (*, *)
            write (*, *) '.. write to outfile.forceconstant_firstorder'
            call fc1%writetofile(uc, 'outfile.forceconstant_firstorder')
            write (*, *)
        end if
    end if
    if (map%have_fc_pair) then
        call map%get_secondorder_forceconstant(uc, fc2, mem, opts%verbosity)
        call fc2%get_elastic_constants(uc, mw, opts%verbosity)
        if (mw%talk) then
            call fc2%writetofile(uc, 'outfile.forceconstant')
            write (*, *) ''
            write (*, *) 'ELASTIC CONSTANTS (short range) (GPa):'
            do i = 1, 6
                write (*, "(6(3X,F15.5))") lo_chop(fc2%elastic_constants_voigt(:, i)*lo_pressure_HartreeBohr_to_GPa, lo_tol)
            end do
            write (*, *) 'ELASTIC CONSTANTS (long range) (GPa):'
            do i = 1, 6
                write (*, "(6(3X,F15.5))") lo_chop(fc2%elastic_constants_voigt_longrange(:, i)*lo_pressure_HartreeBohr_to_GPa, lo_tol)
            end do
            write (*, *) 'ELASTIC CONSTANTS (total) (GPa):'
            do i = 1, 6
                write (*, "(6(3X,F15.5))") lo_chop((fc2%elastic_constants_voigt(:, i) + fc2%elastic_constants_voigt_longrange(:, i))*lo_pressure_HartreeBohr_to_GPa, lo_tol)
            end do
        end if

        !@FIXME add a proper option once it works
        call fc2%write_to_qe(uc,mw,mem)
    end if
    if (map%have_fc_triplet) then
        call map%get_thirdorder_forceconstant(uc, fc3)
        if (mw%talk) call fc3%writetofile(uc, 'outfile.forceconstant_thirdorder')
    end if
    if (map%have_fc_quartet) then
        call map%get_fourthorder_forceconstant(uc, fc4)
        if (mw%talk) call fc4%writetofile(uc, 'outfile.forceconstant_fourthorder')
    end if
    if (map%polar .gt. 1) then
        call map%get_dielectric_tensors(uc, di)
        if (mw%talk) call di%writetofile(uc, 'outfile.dielectric_interaction')
    end if
    ! if ( map%magnetic_pair_interactions .and. mw%talk ) then
    !     call map%get_secondorder_jij(uc,jij)
    !     call jij%writetofile(uc,'outfile.jij')
    ! endif
end block getfc

! And calculate U0
getU0: block
    type(lo_forceconstant_secondorder) :: fc2_ss
    type(lo_forceconstant_thirdorder) :: fc3_ss
    type(lo_forceconstant_fourthorder) :: fc4_ss
    real(r8), dimension(:, :, :, :), allocatable :: polarfc
    real(r8), dimension(:, :, :), allocatable :: f0, fp, f2, f3, f4
    real(r8), dimension(:), allocatable :: e0, ep, e2, e3, e4, ebuf
    real(r8), dimension(3, 3, 3, 3) :: m4
    real(r8), dimension(3, 3, 3) :: m3
    real(r8), dimension(3, 3) :: m2
    real(r8), dimension(3) :: v0, u2, u3, u4
    real(r8) :: energy, baseline, tomev, toev
    real(r8) :: mf0, mfp, mf2, mf3, mf4
    real(r8) :: sf0, sfp, sf2, sf3, sf4
    real(r8) :: dfp, df2, df3, df4
    integer :: t, i, a1, a2, a3, a4, i1, i2, i3, i4, u, ctr, ctrtot

    ! Make sure we have the supercell and MD data
    if (ss%na .lt. 0) then
        call ss%readfromfile('infile.ssposcar')
        call ss%classify('supercell', uc)
    end if
    if (sim%na .lt. 0) then
        call sim%read_from_file(verbosity=2, stride=opts%stride)
        call sim%remove_force_and_center_of_mass_drift()
    end if

    ! Remap the forceconstants to the supercell
    if (map%have_fc_pair) call fc2%remap(uc, ss, fc2_ss)
    if (map%have_fc_triplet) call fc3%remap(uc, ss, fc3_ss)
    if (map%have_fc_quartet) call fc4%remap(uc, ss, fc4_ss)
    if (map%polar .gt. 0) then
        allocate (polarfc(3, 3, ss%na, ss%na))
        polarfc = 0.0_r8
        call fc2%supercell_longrange_dynamical_matrix_at_gamma(ss, polarfc, 1E-15_r8)
    end if

    ! Calculate all the different energies and forces for report and diagnostics.
    allocate (f0(3, ss%na, sim%nt))
    allocate (fp(3, ss%na, sim%nt))
    allocate (f2(3, ss%na, sim%nt))
    allocate (f3(3, ss%na, sim%nt))
    allocate (f4(3, ss%na, sim%nt))
    allocate (e0(sim%nt))
    allocate (ep(sim%nt))
    allocate (e2(sim%nt))
    allocate (e3(sim%nt))
    allocate (e4(sim%nt))
    allocate (ebuf(sim%nt))
    f0 = 0.0_r8
    fp = 0.0_r8
    f2 = 0.0_r8
    f3 = 0.0_r8
    f4 = 0.0_r8
    e0 = 0.0_r8
    ep = 0.0_r8
    e2 = 0.0_r8
    e3 = 0.0_r8
    e4 = 0.0_r8
    ebuf = 0.0_r8

    ! Some unit conversion
    tomev = lo_Hartree_to_eV*1000/ss%na
    toev = lo_Hartree_to_eV/ss%na

    if ((mw%talk) .and. (opts%verbosity > 0)) then
        write (*, *) ''
        write (*, *) 'CALCULATING POTENTIAL ENERGIES (meV/atom)'
        write (*, '(A)') '   conf              Epot                  Epolar                &
        &Epair                 Etriplet              Equartet'
    end if

    ctr = 0
    ctrtot = 0
    do t = 1, sim%nt
        if (mod(t, mw%n) .ne. mw%r) cycle
        ctrtot = ctrtot + 1
    end do

    ! Calculate energies and stuff
    do t = 1, sim%nt
        ! make it parallel to not confuse anyone
        if (mod(t, mw%n) .ne. mw%r) cycle
        ! Copy of DFT force and energy
        e0(t) = sim%stat%potential_energy(t)
        f0(:, :, t) = sim%f(:, :, t)
        ! then the pair term
        if (map%have_fc_pair) then
            energy = 0.0_r8
            do a1 = 1, ss%na
                v0 = 0.0_r8
                do i1 = 1, fc2_ss%atom(a1)%n
                    a2 = fc2_ss%atom(a1)%pair(i1)%i2
                    m2 = fc2_ss%atom(a1)%pair(i1)%m
                    v0 = v0 - matmul(m2, sim%u(:, a2, t))
                end do
                energy = energy - dot_product(sim%u(:, a1, t), v0)*0.5_r8
                f2(:, a1, t) = v0
            end do
            e2(t) = energy
        end if
        ! Possible polar term?
        if (fc2%polar) then
            energy = 0.0_r8
            do a1 = 1, ss%na
                v0 = 0.0_r8
                do a2 = 1, ss%na
                    v0 = v0 - matmul(polarfc(:, :, a1, a2), sim%u(:, a2, t))
                end do
                energy = energy - dot_product(sim%u(:, a1, t), v0)*0.5_r8
                fp(:, a1, t) = v0
            end do
            ep(t) = energy
        end if
        ! triplet term
        if (map%have_fc_triplet) then
            energy = 0.0_r8
            do a1 = 1, fc3_ss%na
                v0 = 0.0_r8
                do i = 1, fc3_ss%atom(a1)%n
                    m3 = fc3_ss%atom(a1)%triplet(i)%m
                    a2 = fc3_ss%atom(a1)%triplet(i)%i2
                    a3 = fc3_ss%atom(a1)%triplet(i)%i3
                    u2 = sim%u(:, a2, t)
                    u3 = sim%u(:, a3, t)
                    do i1 = 1, 3
                    do i2 = 1, 3
                    do i3 = 1, 3
                        v0(i1) = v0(i1) - m3(i1, i2, i3)*u2(i2)*u3(i3)
                    end do
                    end do
                    end do
                end do
                v0 = v0*0.5_r8
                f3(:, a1, t) = v0
                energy = energy - dot_product(v0, sim%u(:, a1, t))/3.0_r8
            end do
            e3(t) = energy
        end if
        ! quartet term
        if (map%have_fc_quartet) then
            energy = 0.0_r8
            do a1 = 1, fc4_ss%na
                v0 = 0.0_r8
                do i = 1, fc4_ss%atom(a1)%n
                    m4 = fc4_ss%atom(a1)%quartet(i)%m
                    a2 = fc4_ss%atom(a1)%quartet(i)%i2
                    a3 = fc4_ss%atom(a1)%quartet(i)%i3
                    a4 = fc4_ss%atom(a1)%quartet(i)%i4
                    u2 = sim%u(:, a2, t)
                    u3 = sim%u(:, a3, t)
                    u4 = sim%u(:, a4, t)
                    do i1 = 1, 3
                    do i2 = 1, 3
                    do i3 = 1, 3
                    do i4 = 1, 3
                        v0(i1) = v0(i1) - m4(i1, i2, i3, i4)*u2(i2)*u3(i3)*u4(i4)
                    end do
                    end do
                    end do
                    end do
                end do
                v0 = v0/6.0_r8
                f4(:, a1, t) = v0
                energy = energy - dot_product(v0, sim%u(:, a1, t))/4.0_r8
            end do
            e4(t) = energy
        end if
        if ((mw%talk) .and. (opts%verbosity > 0)) then
            ctr = ctr + 1
            if (lo_trueNtimes(ctr, 20, ctrtot)) then
                write (*, "(1X,I5,F30.12,4(2X,F20.12))") t, e0(t)*tomev, ep(t)*tomev, e2(t)*tomev, e3(t)*tomev, e4(t)*tomev
            end if
        end if
    end do

    ! sync across ranks
    call mw%allreduce('sum', f0)
    call mw%allreduce('sum', fp)
    call mw%allreduce('sum', f2)
    call mw%allreduce('sum', f3)
    call mw%allreduce('sum', f4)
    call mw%allreduce('sum', e0)
    call mw%allreduce('sum', ep)
    call mw%allreduce('sum', e2)
    call mw%allreduce('sum', e3)
    call mw%allreduce('sum', e4)

    ! dump to outfile.enegies
    if ((mw%talk) .and. (opts%verbosity > 0)) then
        ! file for the energies
        u = open_file('out', 'outfile.energies')
        write (u, '(A,A)') '# Unit:      ', 'eV/atom'
        write (u, '(A,A)') '# no. atoms: ', tochar(ss%na)
        write (u, "(A)") '#  conf    Epot                  Epolar                &
            &Epair                 Etriplet              Equartet'

        do t = 1, sim%nt
            ! Dump it to file
            write (u, "(1X,I5,5(2X,E20.12))") t, e0(t)*toev, ep(t)*toev, e2(t)*toev, e3(t)*toev, e4(t)*toev
        end do

        ! close outfile.energies
        write (*, '(A)') ' ... energies writen to `outfile.energies`'
        close (u)
    end if

    ! Subtract a baseline to get sensible numbers to work with
    ebuf = e0 - e2 - e3 - e4 - ep
    baseline = lo_mean(ebuf)
    e0 = e0 - baseline

    ! Now I can report a little?
    if (mw%talk) then

        write (*, *) ''
        write (*, "(1X,A)") 'ENERGIES (meV/atom):'
        write (*, "(21X,A,10X,A,10X,A)") 'rms', 'std', 'std(res)'
        write (*, '(1X,A15,2(1X,F12.6),5X,A)') 'input:', &
            sqrt(lo_mean(e0**2))*tomev, lo_stddev(e0)*tomev, '-'
        if (map%polar .gt. 0) write (*, '(1X,A15,3(1X,F12.6))') 'polar', &
            sqrt(lo_mean(ep**2))*tomev, lo_stddev(ep)*tomev, lo_stddev(e0 - ep)*tomev
        if (map%have_fc_pair) write (*, '(1X,A15,3(1X,F12.6))') 'second order:', &
            sqrt(lo_mean(e2**2))*tomev, lo_stddev(e2)*tomev, lo_stddev(e0 - ep - e2)*tomev
        if (map%have_fc_triplet) write (*, '(1X,A15,3(1X,F12.6))') 'third order:', &
            sqrt(lo_mean(e3**2))*tomev, lo_stddev(e3)*tomev, lo_stddev(e0 - ep - e2 - e3)*tomev
        if (map%have_fc_quartet) write (*, '(1X,A15,3(1X,F12.6))') 'fourth order:', &
            sqrt(lo_mean(e4**2))*tomev, lo_stddev(e4)*tomev, lo_stddev(e0 - ep - e2 - e3 - e4)*tomev

        ! some things on the forces?
        mf0 = sqrt(sum(f0**2)/sim%na/sim%nt)*lo_force_hartreebohr_to_eVa
        mfp = sqrt(sum(fp**2)/sim%na/sim%nt)*lo_force_hartreebohr_to_eVa
        mf2 = sqrt(sum(f2**2)/sim%na/sim%nt)*lo_force_hartreebohr_to_eVa
        mf3 = sqrt(sum(f3**2)/sim%na/sim%nt)*lo_force_hartreebohr_to_eVa
        mf4 = sqrt(sum(f4**2)/sim%na/sim%nt)*lo_force_hartreebohr_to_eVa

        ! flokno: stddev as well
        sf0 = lo_stddev(f0)*lo_force_hartreebohr_to_eVa
        sfp = lo_stddev(f0 - fp)*lo_force_hartreebohr_to_eVa
        sf2 = lo_stddev(f0 - fp - f2)*lo_force_hartreebohr_to_eVa
        sf3 = lo_stddev(f0 - fp - f2 - f3)*lo_force_hartreebohr_to_eVa
        sf4 = lo_stddev(f0 - fp - f2 - f3 - f4)*lo_force_hartreebohr_to_eVa

        dfp = sqrt(sum((f0 - fp)**2)/sim%na/sim%nt)*lo_force_hartreebohr_to_eVa
        df2 = sqrt(sum((f0 - fp - f2)**2)/sim%na/sim%nt)*lo_force_hartreebohr_to_eVa
        df3 = sqrt(sum((f0 - fp - f2 - f3)**2)/sim%na/sim%nt)*lo_force_hartreebohr_to_eVa
        df4 = sqrt(sum((f0 - fp - f2 - f3 - f4)**2)/sim%na/sim%nt)*lo_force_hartreebohr_to_eVa

        write (*, *) ''
        write (*, "(1X,A)") 'FORCES (eV/A):'
        write (*, "(21X,A,10X,A,5X,A,5X,A,5X,A)") 'rms', 'rms(res)', 'std(res)', 'R^2(res)', 'normalized std(res)'
        write (*, '(1X,A15,3(1X,F12.6),5X,A,12X,A)') 'input:', mf0, mf0, sf0, '-', '-'
        if (map%polar .gt. 0) write (*, '(1X,A15,5(1X,F12.6))') 'polar:', mfp, dfp, sfp, 1 - (sfp/sf0)**2, sfp/sf0
        if (map%have_fc_pair) write (*, '(1X,A15,5(1X,F12.6),13X,A)') 'second order:', mf2, df2, sf2, 1 - (sf2/sf0)**2, sf2/sf0, '<-- anharmonicity measure'
        if (map%have_fc_triplet) write (*, '(1X,A15,5(1X,F12.6))') 'third order:', mf3, df3, sf3, 1 - (sf3/sf0)**2, sf3/sf0
        if (map%have_fc_quartet) write (*, '(1X,A15,5(1X,F12.6))') 'fourth order:', mf4, df4, sf4, 1 - (sf4/sf0)**2, sf4/sf0

        write (*, *) ''
        write (*, *) 'BASELINE ENERGY (eV/atom):'
        write (*, '(1X,A15,1X,F25.6)') '            U0:', baseline*toev

        ! Dump it to file
        u = open_file('out', 'outfile.U0')
        write (u, '(A,A)') '# Unit:      ', 'eV/atom'
        write (u, '(A,A)') '# no. atoms: ', tochar(ss%na)
        write (u, "(A)") '#            mean(Epot)                     mean(Epot - Epolar - E2) &
            &      mean(Epot - Epolar - E2 - E3)  mean(Epot - Epolar - E2 - E3 - E4)'

        write (u, "(4(1X,E30.12))") &
            (baseline + lo_mean(e0))*toev, &
            (baseline + lo_mean(e0 - ep - e2))*toev, &
            (baseline + lo_mean(e0 - ep - e2 - e3))*toev, &
            (baseline + lo_mean(e0 - ep - e2 - e3 - e4))*toev
        close (u)
    end if

end block getU0

! predict new reference structure
if (map%have_fc_singlet) then
    update_reference: block
        !   [ ] worry about long-range later
        !   [ ] maybe there is a more clever way as usual
        integer :: ii, jj
        real(r8), dimension(:), allocatable :: f, s, im
        complex(r8), dimension(:, :), allocatable :: D
        type(lo_qpoint) :: gammapoint

        if (mw%talk) then
            write (*, *) ''
            write (*, *) &
                'Since we have the first order, predict new reference positions.'
        end if

        associate (n => uc%na)
            allocate (f(3*n))
            allocate (s(3*n))
            allocate (im(3*n))
            allocate (D(3*n, 3*n))
        end associate

        ! init
        f = 0.0_r8
        s = 0.0_r8
        im = 0.0_r8
        D = 0.0_r8
        gammapoint%r = 0.0_8

        ! make sure supercell knows about unitcell
        call ss%classify('supercell', uc)

        ! get dynamical matrix
        call fc2%dynamicalmatrix(uc, gammapoint, D, mem)

        ! matrix at gamma should be real, sanity check
        associate (norm => lo_frobnorm(aimag(D)))
            if (norm .gt. lo_tol) then
                call lo_stop_gracefully( &
                    ['Dynamical matrix not real by: ', tochar(norm)], &
                    lo_exitcode_param, __FILE__, __LINE__, mw%comm)
            end if
        end associate

        ! arrange mass-weighted forces and mass weights
        do ii = 1, uc%na
            f(3*(ii - 1) + 1:3*ii) = fc1%atom(ii)%m*uc%invsqrtmass(ii)
            im(3*(ii - 1) + 1:3*ii) = uc%invsqrtmass(ii)
        end do

        ! predict step as D.T @ s = f
        call lo_linear_least_squares(transpose(real(D)), f, s)

        ! undo mass-weighting
        s = s*im

        ! update unitcell positions
        do ii = 1, uc%na
            associate (disp => s(3*(ii - 1) + 1:3*ii))
                ! in fractional coords
                associate (d => uc%cartesian_to_fractional(disp, pbc=.false.), &
                           r => uc%r(:, ii))
                    r = r - d
                end associate
            end associate
        end do

        ! update supercell positions
        do jj = 1, ss%na
            ii = ss%info%index_in_unitcell(jj)
            associate (disp => s(3*(ii - 1) + 1:3*ii))
                ! add displacement in fractional coords
                associate (d => ss%cartesian_to_fractional(disp, pbc=.false.), &
                           r => ss%r(:, jj))
                    r = r - d
                end associate
            end associate
        end do

        ! dump files
        call uc%writetofile('outfile.new_ucposcar', 1)
        call ss%writetofile('outfile.new_ssposcar', 1)

    end block update_reference
end if

if (mw%talk) then
    write (*, *) ''
    write (*, *) 'Done in ', tochar(walltime() - t0), 's'
end if
call mw%destroy()

end program
