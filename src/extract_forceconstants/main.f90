program extract_forceconstants
!!{!src/extract_forceconstants/manual.md!}!
use konstanter, only: r8, i8, lo_iou, lo_tol, lo_bohr_to_A, lo_pressure_HartreeBohr_to_GPa, &
                      lo_Hartree_to_eV, lo_exitcode_param, &
                      lo_forceconstant_1st_HartreeBohr_to_eVA, lo_force_hartreebohr_to_eVa, &
                      lo_kb_Hartree
use gottochblandat, only: tochar, walltime, open_file, lo_chop, lo_does_file_exist, &
                          lo_trueNtimes, lo_mean, lo_stddev, lo_outerproduct, lo_trace, &
                          lo_frobnorm, lo_linear_least_squares, lo_invert_real_matrix
use mpi_wrappers, only: lo_mpi_helper, lo_stop_gracefully
use lo_memtracker, only: lo_mem_helper
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
use symmetrize, only: symmetrize_stress, symmetrize_forces
use helper, only: write_stress_3x3, write_stress_voigt, write_file_stress_voigt_time, &
                  write_file_stress_voigt
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
    if (mw%talk .eqv. .false.) opts%verbosity = -100
    if (mw%talk) write (lo_iou, *) '... reading unitcell'
    call uc%readfromfile('infile.ucposcar', verbosity=opts%verbosity)
    call uc%classify('wedge', timereversal=.true.)
end block init

! We need a forcemap to continue.
getfrcmap: block
    type(lo_interaction_tensors) :: slt
    real(r8) :: t0

    ! Either read it from file or create it
    if (opts%readforcemap) then
        ! I might need the supercell anyway
        if (opts%readirreducible .eqv. .false.) then
            call ss%readfromfile('infile.ssposcar')
            call ss%classify('supercell', uc)
        end if
        call map%read_from_hdf5(uc, 'infile.forcemap.hdf5', opts%verbosity)
        ! Update constraints that depend on structure
        call lo_secondorder_rot_herm_huang(map, uc, map%constraints%eq2, map%constraints%neq2, &
                                           opts%rotationalconstraints, opts%huanginvariance, opts%hermitian)
    else
        ! Now we have to calculate the whole thing.
        call uc%classify('wedge', timereversal=.true.)
        call ss%readfromfile('infile.ssposcar')
        if (mw%talk) write (*, *) '... min cutoff: ', tochar(ss%mincutoff()*lo_bohr_to_A)
        if (mw%talk) write (*, *) '... max cutoff: ', tochar(ss%maxcutoff()*lo_bohr_to_A)
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

        ! Create the constraints
        t0 = walltime()
        call map%forceconstant_constraints(uc, opts%rotationalconstraints, &
                                           opts%huanginvariance, opts%hermitian, opts%verbosity + 10)
        if (mw%talk) write (*, *) '... got constraints in ', tochar(walltime() - t0), 's'
        ! Maybe print forcemap to file.
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
    integer :: i, uu, ii
    real(r8) :: compressibility, bulk_modulus
    real(r8), dimension(3) :: force_norm
    real(r8), dimension(6, 6) :: c_ij, s_ij  ! elastic constants, compliances

    c_ij = 0.0_r8
    s_ij = 0.0_r8
    compressibility = 0.0_r8
    bulk_modulus = 0.0_r8

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
        if (mw%talk) then
            call fc2%writetofile(uc, 'outfile.forceconstant')

            ! elastic properties
            call fc2%get_elastic_constants(uc)
            c_ij = fc2%elastic_constants_voigt*lo_pressure_HartreeBohr_to_GPa
            s_ij = fc2%elastic_compliances_voigt/lo_pressure_HartreeBohr_to_GPa

            write (*, *) ''
            write (*, *) 'ELASTIC CONSTANTS (GPa):'

            do i = 1, 6
                write (*, "(6(3X,F15.5))") lo_chop(c_ij(:, i), lo_tol)
            end do

            ! compressibility and bulk modulus (Nye p. 146)
            compressibility = s_ij(1, 1) + s_ij(2, 2) + s_ij(3, 3)
            compressibility = compressibility + 2.0_r8*(s_ij(1, 2) + s_ij(2, 3) + s_ij(3, 1))
            bulk_modulus = 1.0_r8/compressibility
            write (*, *) '--> Isothermal compressibility (1/GPa):'
            write (*, "(3X,F15.5)") compressibility
            write (*, *) '--> Isothermal bulk modulus (GPa):'
            write (*, "(3X,F15.5)") bulk_modulus

            ! fkdev: write to file
            write (*, *) '.. write compliance tensor to outfile.compliance_tensor'
            uu = open_file('out', 'outfile.compliance_tensor')

            write (uu, *) "# Compliance tensor in Voigt notation (1/GPa)"
            do ii = 1, 6
                write (uu, "(1X,6(1X,F20.15))") s_ij(:, ii)
            end do
            close (uu)

        end if
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
        write (*, '(A)') '   conf        Epot                  Epolar                &
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
                write (*, "(1X,I5,5(2X,F20.12))") t, e0(t)*tomev, ep(t)*tomev, e2(t)*tomev, e3(t)*tomev, e4(t)*tomev
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
        write (u, "(4(1X,E19.12))") &
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
        integer :: ii, jj, ip, a1
        real(r8), dimension(:), allocatable :: f, s, im
        real(r8), dimension(3, sim%na) :: f0, f0_std  ! , f2, f3
        real(r8), dimension(3, uc%na) :: f0_prim, f0_prim_sym  ! , f2_prim, f2_prim_sym, f3_prim, f3_prim_sym
        real(r8), dimension(3, uc%na) :: f0_prim_std, f0_prim_err, f0_prim_err_sym
        complex(r8), dimension(:, :), allocatable :: D
        type(lo_qpoint) :: gammapoint
        type(lo_crystalstructure) :: uc_new, ss_new

        if (mw%talk) then
            write (*, *) ''
            write (*, *) &
                '--> since ie have the first order, predict new reference positions', &
                ' ind write to outfile.??poscar.new_refpos'
        end if

        ! get fresh supercells which we can change
        call uc_new%readfromfile('infile.ucposcar')
        call ss_new%readfromfile('infile.ssposcar')
        call uc_new%classify('wedge', timereversal=.true.)
        call ss_new%classify('supercell', uc_new)

        associate (n => uc_new%na)
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

        ! get dynamical matrix
        call fc2%dynamicalmatrix(uc_new, gammapoint, D, mem)

        ! matrix at gamma should be real, sanity check
        associate (norm => lo_frobnorm(aimag(D)))
            if (norm .gt. lo_tol) then
                call lo_stop_gracefully( &
                    ['Dynamical matrix not real by: ', tochar(norm)], &
                    lo_exitcode_param, __FILE__, __LINE__, mw%comm)
            end if
        end associate

        ! arrange mass-weighted forces and mass weights
        do ii = 1, uc_new%na
            f(3*(ii - 1) + 1:3*ii) = fc1%atom(ii)%m*uc_new%invsqrtmass(ii)
            im(3*(ii - 1) + 1:3*ii) = uc_new%invsqrtmass(ii)
        end do

        ! predict step as D.T @ s = f
        call lo_linear_least_squares(transpose(real(D)), f, s)

        ! undo mass-weighting
        s(:) = s*im

        ! update unitcell positions
        do ii = 1, uc_new%na
            associate (disp => s(3*(ii - 1) + 1:3*ii))
                ! in fractional coords
                associate (d => uc_new%cartesian_to_fractional(disp, pbc=.false.), &
                           r => uc_new%r(:, ii))
                    r = r - d
                end associate
            end associate
        end do

        ! update supercell positions
        do jj = 1, ss_new%na
            ii = ss_new%info%index_in_unitcell(jj)
            associate (disp => s(3*(ii - 1) + 1:3*ii))
                ! add displacement in fractional coords
                associate (d => ss_new%cartesian_to_fractional(disp, pbc=.false.), &
                           r => ss_new%r(:, jj))
                    r = r - d
                end associate
            end associate

        end do
        ! dump files
        call uc_new%writetofile('outfile.ucposcar.new_refpos', 1)
        call ss_new%writetofile('outfile.ssposcar.new_refpos', 1)

        ! some force statistic
        ! I used this to compare to the first-order FC
        f0 = 0.0_r8
        f0_std = 0.0_r8
        ! f2 = 0.0_r8
        ! f3 = 0.0_r8
        f0_prim = 0.0_r8
        f0_prim_std = 0.0_r8
        f0_prim_err = 0.0_r8
        f0_prim_err_sym = 0.0_r8
        ! f2_prim = 0.0_r8
        ! f3_prim = 0.0_r8
        do a1 = 1, ss%na
            ip = ss%info%index_in_unitcell(a1)
            do ii = 1, 3
                f0(ii, a1) = lo_mean(sim%f(ii, a1, :))
                f0_std(ii, a1) = lo_stddev(sim%f(ii, a1, :))
                ! f2(i, a1) = lo_mean(f2t(i, a1, :))
                ! f3(i, a1) = lo_mean(f3t(i, a1, :))
            end do
            f0_prim(:, ip) = f0_prim(:, ip) + f0(:, a1)
            f0_prim_std(:, ip) = f0_prim(:, ip) + f0(:, a1)
            ! f2_prim(:, ip) = f2_prim(:, ip) + f2(:, a1)
            ! f3_prim(:, ip) = f3_prim(:, ip) + f3(:, a1)
        end do

        f0_prim = f0_prim/ss%na*uc%na
        f0_prim_err = f0_prim_std/ss%na*uc%na/sqrt(1.0*sim%nt)
        ! f2_prim = f2_prim/ss%na*uc%na
        ! f3_prim = f3_prim/ss%na*uc%na

        ! write (lo_iou, *) 'Mean force:'
        ! write (*, *) sum(f0_prim, dim=2)/uc%na
        ! write (*, *) sum(f2_prim, dim=2)/uc%na
        ! write (*, *) sum(f3_prim, dim=2)/uc%na

        ! ! symmetrize the averaged forces
        call symmetrize_forces(f0_prim, uc, f0_prim_sym, verbose=.false.)
        call symmetrize_forces(f0_prim_err, uc, f0_prim_err_sym, verbose=.false.)
        ! call symmetrize_forces(f2_prim, uc, f2_prim_sym)
        ! call symmetrize_forces(f3_prim, uc, f3_prim_sym)

        if (mw%talk) then
            write (lo_iou, *)
            write (lo_iou, *) 'Averaged forces in primitive cell (symmetrized, should be == firstorder FC) with 95% confidence interval:'
            do ip = 1, uc%na
                write (*, "(1X,I4,3(ES11.3),'  +/- ',3(ES11.3))") ip, &
                    f0_prim_sym(:, ip)*lo_force_hartreebohr_to_eVa, 1.96*f0_prim_err(:, ip)*lo_force_hartreebohr_to_eVa
            end do
        end if

        ! if (map%have_fc_pair) then
        !     write (lo_iou, *) 'Averaged harmonic forces in primitive cell (symmetrized):'
        !     do ip = 1, uc%na
        !         write (*, *) ip, f2_prim_sym(:, ip)*lo_force_hartreebohr_to_eVa
        !     end do
        ! end if

        ! if (map%have_fc_triplet) then
        !     write (lo_iou, *) 'Averaged third-order forces in primitive cell (symmetrized):'
        !     do ip = 1, uc%na
        !         write (*, *) ip, f3_prim_sym(:, ip)*lo_force_hartreebohr_to_eVa
        !     end do
        ! end if

    end block update_reference
end if

! pressure business
if (opts%pressure) then
    ! get realspace stress from the force constants:
    ! two contributions:
    ! - [x] harmonic contribution \Phi_ij u_i u_j
    ! - [x] quasi-harmonic contribution \Phi'_ij u_i u_j where \Phi' comes from 3. order FC

    pressure_tdep: block
        ! this case is similar to harmonic energy

        type(lo_forceconstant_secondorder) :: fc2_ss
        type(lo_forceconstant_thirdorder) :: fc3_ss
        type(lo_forceconstant_fourthorder) :: fc4_ss
        real(r8), dimension(:, :, :, :), allocatable :: polarfc
        real(r8), dimension(3, 3, sim%nt) :: s0t, s2t, spt, s3t, s4t, srt
        real(r8), dimension(3, sim%na, sim%nt) :: f0t, fpt, f2t, f3t, f4t
        real(r8), dimension(sim%nt) :: e0t, ept, e2t, e3t, e4t, ert, p0t, p2t, ppt, p3t, p4t, prt
        real(r8), dimension(3, 3, 3) :: m3
        real(r8), dimension(3, 3, 3, 3) :: m4
        real(r8), dimension(3, 3) :: m2, s0, s2, sp, s3, s4, sr
        real(r8), dimension(3, 3) :: s0_std, s2_std, s3_std, s4_std, sr_std
        real(r8), dimension(3, 3) :: s0_err, s2_err, s3_err, s4_err, sr_err
        real(r8), dimension(3, 3) :: s0_sym, s2_sym, s3_sym, sr_sym  ! s4_sym
        real(r8), dimension(3, 3) :: s0_std_sym, s2_std_sym, s3_std_sym, sr_std_sym  ! s4_std_sym
        real(r8), dimension(3, 3) :: s0_err_sym, s2_err_sym, s3_err_sym, sr_err_sym  ! s4_err_sym
        real(r8), dimension(3) :: v0, u1, u2, u3, u4, rv1, rv2, rv3, rv4
        !> averaged force
        ! real(r8), dimension(3, sim%na) :: f0, f2, f3, f4, f0_std
        ! real(r8), dimension(3, uc%na) :: f0_prim, f0_prim_sym, f2_prim, f2_prim_sym, f3_prim, f3_prim_sym
        ! real(r8), dimension(3, uc%na) :: f0_prim_std, f0_prim_err
        ! real(r8) :: dfp, df2, df3, df4
        real(r8) :: n_samples, energy, baseline, tomev, toev, volume
        real(r8) :: p0, p2, pp, p3, p4, pr
        real(r8) :: p0_std, p2_std, pp_std, p3_std, p4_std, pr_std
        real(r8) :: p0_err, p2_err, pp_err, p3_err, p4_err, pr_err
        integer :: i, j, t, a1, a2, a3, a4, ii, i1, i2, i3, i4, i5, u, ctr, ctrtot

        ! Make sure we have the supercell and MD data
        if (ss%na .lt. 0) then
            call ss%readfromfile('infile.ssposcar')
            call ss%classify('supercell', uc)
        end if
        if (sim%na .lt. 0) then
            call sim%read_from_file(verbosity=2, stride=opts%stride)
            call sim%remove_force_and_center_of_mass_drift()
        end if

        ! system volume and number of samples
        n_samples = 1.0_r8*sim%nt
        volume = ss%volume

        ! Remap the forceconstants to the supercell
        if (map%have_fc_pair) call fc2%remap(uc, ss, fc2_ss)
        ! fkdev: thirdorder will become sth like secondorder, later problem
        if (map%have_fc_triplet) call fc3%remap(uc, ss, fc3_ss)
        if (map%have_fc_quartet) call fc4%remap(uc, ss, fc4_ss)
        if (map%polar .gt. 0) then
            allocate (polarfc(3, 3, ss%na, ss%na))
            polarfc = 0.0_r8
            call fc2%supercell_longrange_dynamical_matrix_at_gamma(ss, polarfc, 1E-15_r8)
        end if

        ! init
        p0t = 0.; p2t = 0.; ppt = 0.; p3t = 0.; p4t = 0.; prt = 0.; 
        p0_std = 0.; p2_std = 0.; pp_std = 0.; p3_std = 0.; p0_err = 0.; p2_err = 0.; pp_err = 0.; p3_err = 0.; 
        s0t = 0.; spt = 0.; s2t = 0.; s3t = 0.; s4t = 0.; sr = 0.; srt = 0.; 
        s0_std = 0.; s2_std = 0.; s3_std = 0.; sr_std = 0.; 
        s0_err = 0.; s2_err = 0.; s3_err = 0.; sr_err = 0.; 
        s0_sym = 0.; s2_sym = 0.; s3_sym = 0.; sr_sym = 0.; 
        s0_std_sym = 0.; s2_std_sym = 0.; s3_std_sym = 0.; 
        s0_err_sym = 0.; s2_err_sym = 0.; s3_err_sym = 0.; 
        f0t = 0.; fpt = 0.; f2t = 0.; f3t = 0.; f4t = 0.; 
        e0t = 0.; ept = 0.; e2t = 0.; e3t = 0.; e4t = 0.; ert = 0.; 
        p0 = 0.; p2 = 0.; pp = 0.; p3 = 0.; p4 = 0.; p0 = 0.; 
        ! Unit conversion
        tomev = lo_Hartree_to_eV*1000/ss%na
        toev = lo_Hartree_to_eV/ss%na

        ctr = 0
        ctrtot = 0
        do t = 1, sim%nt
            if (mod(t, mw%n) .ne. mw%r) cycle
            ctrtot = ctrtot + 1
        end do

        ! Calculate forces and stresses
        do t = 1, sim%nt
            ! make it parallel to not confuse anyone
            if (mod(t, mw%n) .ne. mw%r) cycle
            ! Copy of DFT energy, force, stress
            e0t(t) = sim%stat%potential_energy(t)
            f0t(:, :, t) = sim%f(:, :, t)
            s0t(:, :, t) = sim%stat%stress(:, :, t)
            ! then the harmonic terms
            if (map%have_fc_pair) then
                energy = 0.0_r8
                s2 = 0.0_r8
                do a1 = 1, ss%na
                    v0 = 0.0_r8
                    u1 = sim%u(:, a1, t)
                    do ii = 1, fc2_ss%atom(a1)%n
                        a2 = fc2_ss%atom(a1)%pair(ii)%i2
                        m2 = fc2_ss%atom(a1)%pair(ii)%m
                        rv2 = fc2_ss%atom(a1)%pair(ii)%r
                        v0 = v0 - matmul(m2, sim%u(:, a2, t))
                        ! fkdev: used fixed reference position
                        do i1 = 1, 3
                        do i2 = 1, 3
                        do i3 = 1, 3
                            associate ( &
                                s => s2(i1, i2), &
                                prefactor => 1.0_r8/volume, &
                                phi_prime => m2(i3, i1)*u1(i3)*rv2(i2) &
                                )
                                s = s + prefactor*phi_prime
                            end associate
                        end do
                        end do
                        end do

                    end do
                    associate (f => v0)
                        energy = energy - 0.5_r8*dot_product(u1, f)
                        f2t(:, a1, t) = f
                        s2 = s2 - lo_outerproduct(u1, f)/volume
                    end associate
                end do
                e2t(t) = energy
                s2t(:, :, t) = s2
            end if
            ! Possible polar term?
            if (fc2%polar) then
                energy = 0.0_r8
                sp = 0.0_r8
                do a1 = 1, ss%na
                    v0 = 0.0_r8
                    do a2 = 1, ss%na
                        v0 = v0 - matmul(polarfc(1:3, 1:3, a1, a2), sim%u(1:3, a2, t))
                    end do
                    associate (u => sim%u(:, a1, t), f => v0)
                        energy = energy - dot_product(u, v0)*0.5_r8
                        fpt(:, a1, t) = f
                        sp = sp - 0.5_r8/volume*(lo_outerproduct(u, f) + lo_outerproduct(f, u))
                    end associate

                end do
                ept(t) = energy
                spt(:, :, t) = sp
            end if
            ! triplet term gives quasi-harmonic contribution
            if (map%have_fc_triplet) then
                energy = 0.0_r8
                s3 = 0.0_r8
                do a1 = 1, fc3_ss%na
                    v0 = 0.0_r8
                    do i = 1, fc3_ss%atom(a1)%n
                        m3 = fc3_ss%atom(a1)%triplet(i)%m
                        a2 = fc3_ss%atom(a1)%triplet(i)%i2
                        a3 = fc3_ss%atom(a1)%triplet(i)%i3
                        u1 = sim%u(:, a1, t)
                        u2 = sim%u(:, a2, t)
                        u3 = sim%u(:, a3, t)

                        ! positions
                        rv1 = fc3_ss%atom(a1)%triplet(i)%rv1
                        rv2 = fc3_ss%atom(a1)%triplet(i)%rv2
                        rv3 = fc3_ss%atom(a1)%triplet(i)%rv3

                        do ii = 1, 3
                        do i2 = 1, 3
                        do i3 = 1, 3
                            v0(ii) = v0(ii) - m3(ii, i2, i3)*u2(i2)*u3(i3)
                            do i4 = 1, 3
                                associate ( &
                                    s => s3(ii, i2), &
                                    prefactor => 0.5_r8/volume, &
                                    ! fkdev: use r_i = r0_i + u_i instead of r0_i
                                    phi_prime1 => m3(i3, i4, i2)*u1(i3)*u2(i4)*rv3(ii) &
                                    ! phi_prime2 => m3(i3, i2, i4)*u1(i3)*u3(i4)*rv2(i1) &
                                    )
                                    ! s = s + prefactor*(phi_prime1 + phi_prime2 )/2.0_r8
                                    s = s + prefactor*phi_prime1
                                end associate
                            end do
                        end do
                        end do
                        end do
                    end do
                    v0 = v0*0.5_r8
                    f3t(:, a1, t) = v0

                    ! fkdev: add energy term to stress similar to 2nd order
                    ! comment: doesn't seem to reduce variance/fluctuations
                    associate (u => sim%u(:, a1, t), f => v0)
                        s3 = s3 - lo_outerproduct(u, f)/volume
                    end associate
                    energy = energy - dot_product(v0, sim%u(:, a1, t))/3.0_r8
                    ! fkdev: I think this cannot work b.c. p.b.c but I tried
                    ! sft(:, :, t) = sft(:, :, t) - lo_outerproduct(ss%rcart(:, a1), v0)/volume/2.0_r8
                end do
                e3t(t) = energy
                s3t(:, :, t) = s3
            end if
            ! quartet term
            if (map%have_fc_quartet) then
                energy = 0.0_r8
                s4 = 0.0_r8
                do a1 = 1, fc4_ss%na
                    v0 = 0.0_r8
                    do i = 1, fc4_ss%atom(a1)%n
                        m4 = fc4_ss%atom(a1)%quartet(i)%m
                        a2 = fc4_ss%atom(a1)%quartet(i)%i2
                        a3 = fc4_ss%atom(a1)%quartet(i)%i3
                        a4 = fc4_ss%atom(a1)%quartet(i)%i4
                        u1 = sim%u(:, a1, t)
                        u2 = sim%u(:, a2, t)
                        u3 = sim%u(:, a3, t)
                        u4 = sim%u(:, a4, t)

                        ! positions
                        rv3 = fc4_ss%atom(a1)%quartet(i)%rv3
                        rv4 = fc4_ss%atom(a1)%quartet(i)%rv4

                        do ii = 1, 3
                        do i2 = 1, 3
                        do i3 = 1, 3
                        do i4 = 1, 3
                            v0(ii) = v0(ii) - m4(ii, i2, i3, i4)*u2(i2)*u3(i3)*u4(i4)
                            do i5 = 1, 3
                                associate ( &
                                    s => s4(ii, i2), &
                                    prefactor => 1.0_r8/volume/6.0_r8, &
                                    phi_prime => u1(i3)*u2(i4)*u3(i5)*m4(i3, i4, i5, ii)*rv4(i2) &
                                    )
                                    s = s + prefactor*phi_prime
                                end associate

                            end do
                        end do
                        end do
                        end do
                        end do
                    end do
                    v0 = v0/6.0_r8
                    f4t(:, a1, t) = v0

                    ! fkdev: add energy term to stress similar to 2nd order
                    ! associate (u => sim%u(:, a1, t), f => v0)
                    !     s4 = s4 - lo_outerproduct(u, f)/volume
                    ! end associate

                    energy = energy - dot_product(v0, sim%u(:, a1, t))/4.0_r8
                end do
                e4t(t) = energy
                s4t(:, :, t) = s4
            end if
        end do

        ! sync across ranks
        call mw%allreduce('sum', f0t)
        call mw%allreduce('sum', fpt)
        call mw%allreduce('sum', f2t)
        call mw%allreduce('sum', f3t)
        call mw%allreduce('sum', f4t)
        call mw%allreduce('sum', e0t)
        call mw%allreduce('sum', ept)
        call mw%allreduce('sum', e2t)
        call mw%allreduce('sum', e3t)
        call mw%allreduce('sum', e4t)
        call mw%allreduce('sum', s0t)
        call mw%allreduce('sum', s2t)
        call mw%allreduce('sum', spt)
        call mw%allreduce('sum', s3t)
        call mw%allreduce('sum', s4t)

        ! Subtract a baseline to get sensible numbers to work with
        ert = e0t - e2t - e3t - e4t - ept
        ! stress residual
        srt = s0t - s2t - spt - s3t - s4t
        baseline = lo_mean(ert)
        e0t = e0t - baseline

        ! do everything else in GPa for digestive purposes
        s0t = s0t*lo_pressure_HartreeBohr_to_GPa
        spt = spt*lo_pressure_HartreeBohr_to_GPa
        s2t = s2t*lo_pressure_HartreeBohr_to_GPa
        s3t = s3t*lo_pressure_HartreeBohr_to_GPa
        s4t = s4t*lo_pressure_HartreeBohr_to_GPa
        srt = srt*lo_pressure_HartreeBohr_to_GPa

        do i = 1, 3
            do j = 1, 3
                s0(i, j) = lo_mean(s0t(i, j, :))
                s2(i, j) = lo_mean(s2t(i, j, :))
                s3(i, j) = lo_mean(s3t(i, j, :))
                s4(i, j) = lo_mean(s4t(i, j, :))
                sr(i, j) = lo_mean(srt(i, j, :))
                s0_std(i, j) = lo_stddev(s0t(i, j, :)) ! /sqrt(n_samples)
                s2_std(i, j) = lo_stddev(s2t(i, j, :)) ! /sqrt(n_samples)
                s3_std(i, j) = lo_stddev(s3t(i, j, :)) ! /sqrt(n_samples)
                s4_std(i, j) = lo_stddev(s4t(i, j, :))
                sr_std(i, j) = lo_stddev(srt(i, j, :)) ! /sqrt(n_samples)
            end do
            p0t = p0t - s0t(i, i, :)/3.0_r8
            p2t = p2t - s2t(i, i, :)/3.0_r8
            p3t = p3t - s3t(i, i, :)/3.0_r8
            p4t = p4t - s4t(i, i, :)/3.0_r8
            ppt = ppt - spt(i, i, :)/3.0_r8
            prt = prt - srt(i, i, :)/3.0_r8
        end do

        ! error
        s0_err = s0_std/sqrt(n_samples)
        s2_err = s2_std/sqrt(n_samples)
        s3_err = s3_std/sqrt(n_samples)
        s4_err = s4_std/sqrt(n_samples)
        sr_err = sr_std/sqrt(n_samples)

        ! pressure
        p0 = lo_mean(p0t)
        p2 = lo_mean(p2t)
        p3 = lo_mean(p3t)
        p4 = lo_mean(p4t)
        pp = lo_mean(ppt)
        pr = lo_mean(prt)
        p0_std = lo_stddev(p0t)
        p2_std = lo_stddev(p2t)
        pp_std = lo_stddev(ppt)
        p3_std = lo_stddev(p3t)
        p4_std = lo_stddev(p4t)
        pr_std = lo_stddev(prt)
        p0_err = p0_std/sqrt(n_samples)
        p2_err = p2_std/sqrt(n_samples)
        pp_err = pp_std/sqrt(n_samples)
        p3_err = p3_std/sqrt(n_samples)
        p4_err = p4_std/sqrt(n_samples)
        pr_err = pr_std/sqrt(n_samples)

        ! symmetrize things
        call symmetrize_stress(s0, uc, s0_sym)
        call symmetrize_stress(s2, uc, s2_sym)
        call symmetrize_stress(s3, uc, s3_sym)
        call symmetrize_stress(sr, uc, sr_sym)
        call symmetrize_stress(s0_std, uc, s0_std_sym)
        call symmetrize_stress(s2_std, uc, s2_std_sym)
        call symmetrize_stress(s3_std, uc, s3_std_sym)
        call symmetrize_stress(sr_std, uc, sr_std_sym)
        call symmetrize_stress(s0_err, uc, s0_err_sym)
        call symmetrize_stress(s2_err, uc, s2_err_sym)
        call symmetrize_stress(s3_err, uc, s3_err_sym)
        call symmetrize_stress(sr_err, uc, sr_err_sym)

        ! Now I can report a little?
        if (mw%talk) then
            write (*, *) ''
            write (*, "(1X,A,11X,A,10X,A,7X,A)") 'ENERGIES:', 'rms', 'stddev', 'stddev(residual) (meV/atom)'
            write (*, '(1X,A15,2(1X,F12.6),5X,A)') '       input:', sqrt(lo_mean(e0t**2))*tomev, lo_stddev(e0t)*tomev, '-'
            write (*, '(1X,A15,3(1X,F12.6))') '       polar:', sqrt(lo_mean(ept**2))*tomev, lo_stddev(ept)*tomev, lo_stddev(e0t - ept - e2t)*tomev
            write (*, '(1X,A15,3(1X,F12.6))') 'second order:', sqrt(lo_mean(e2t**2))*tomev, lo_stddev(e2t)*tomev, lo_stddev(e0t - e2t)*tomev
            write (*, '(1X,A15,3(1X,F12.6))') 'third  order:', sqrt(lo_mean(e3t**2))*tomev, lo_stddev(e3t)*tomev, lo_stddev(e0t - ept - e2t - e3t)*tomev
            write (*, '(1X,A15,3(1X,F12.6))') 'fourth order:', sqrt(lo_mean(e4t**2))*tomev, lo_stddev(e4t)*tomev, lo_stddev(e0t - ept - e2t - e3t - e4t)*tomev

            write (*, "(1X,A,11X,A,10X,A,5X,A)") 'PRESSURE: ', 'mean', 'stddev', 'stdev (residual) (GPa)'
            write (*, '(1X,A15,2(1X,F12.6),5X,A)') '    input:', p0, p0_std, '-'
            write (*, '(1X,A15,3(1X,F12.6),5X,A)') '    polar:', pp, pp_std, lo_stddev(p0t - p2t - ppt)
            write (*, '(1X,A15,3(1X,F12.6),5X,A)') '2nd order:', p2, p2_std, lo_stddev(p0t - p2t)
            write (*, '(1X,A15,3(1X,F12.6),5X,A)') '3rd order:', p3, p3_std, lo_stddev(p0t - p2t - ppt - p3t)
            write (*, '(1X,A15,3(1X,F12.6),5X,A)') '4th order:', p4, p4_std, lo_stddev(p0t - p2t - ppt - p3t - p4t)
            write (*, '(1X,A15,2(1X,F12.6),5X,A)') ' residual:', pr, pr_std, '-'

            write (*, *)
            associate (prefactor => 2.0_r8/sim%na/lo_kb_Hartree/3.0_r8)
                write (*, "(1X,A,F12.6,A,F12.6,A)") 'Harmonic temperature: ', &
                    prefactor*lo_mean(e2t + ept), ' +/- ', prefactor*lo_stddev(e2t + ept)/sqrt(n_samples), ' K'
            end associate

            write (*, *)
            associate (prefactor => 2.0_r8/3.0_r8/volume*lo_pressure_HartreeBohr_to_GPa)
                write (*, "(1X,A)") 'Harmonic pressure from energy p = 2/3 E_pot / V: '
                write (*, "(1X,F12.6,A,F12.6,A)") &
                    prefactor*lo_mean(e2t + ept), ' +/- ', prefactor*lo_stddev(e2t + ept)/sqrt(n_samples), ' GPa'
            end associate
            write (*, *)

            ! call write_stress_3x3(s0, ds0, 'DFT stress')
            ! call write_stress_3x3(s0_sym, ds0_sym, 'DFT stress (symmetrized')
            call write_stress_voigt(s0, s0_std, s0_err, 'DFT stress (GPa)')
            call write_stress_voigt(s0_sym, s0_std_sym, s0_err_sym, 'DFT stress (symmetrized)')
            ! call write_stress_voigt(s2, s2_std, s2_err, 'Harmonic FC stress (2nd order)')
            call write_stress_voigt(s2_sym, s2_std_sym, s2_err_sym, 'Harmonic FC stress (2nd order) (symmetrized)')
            if (map%have_fc_triplet) then
                ! call write_stress_voigt(s3, s3_std, s3_err, 'Anharmonic FC stress (3rd order)')
                call write_stress_voigt(s3_sym, s3_std_sym, s3_err_sym, 'Anharmonic FC stress (3rd order) (symetrized)')
                ! call write_stress_voigt(sr, sr_std, sr_err, 'Residual stress')
                call write_stress_voigt(sr_sym, sr_std_sym, sr_err_sym, 'Residual stress (symmetrized)')
            end if

            ! dump pressure to file
            write (lo_iou, *) '... dump pressure to file'
            u = open_file('out', 'outfile.pressure.csv')
            write (u, '(a)') '# unit: GPa'
            write (u, '(a)') 'timestep,pressure_dft,pressure_fc2,pressure_fc3,pressure_res,pressure_fc4'
            do i = 1, sim%nt
                write (u, '(I5,5(",",E15.8))') i, p0t(i), p2t(i), p3t(i), prt(i), p4t(i)
            end do
            close (u)

            ! dump stress to file for each time step
            write (lo_iou, *) '... dump DFT stress to file'
            call write_file_stress_voigt_time(s0t, 'outfile.stress_dft.csv', '# DFT stress, unit: GPa')
            write (lo_iou, *) '... dump FC2 stress to file'
            call write_file_stress_voigt_time(s2t, 'outfile.stress_fc2.csv', '# FC2 stress, unit: GPa')
            if (map%have_fc_triplet) then
                ! call write_stress_voigt(s3, s3_std, s3_err, 'Anharmonic FC stress (3rd order)')
                write (lo_iou, *) '... dump FC3 stress to file'
                call write_file_stress_voigt_time(s3t, 'outfile.stress_fc3.csv', '# FC3 stress, unit: GPa')
                write (lo_iou, *) '... dump res. stress to file'
                call write_file_stress_voigt_time(srt, 'outfile.stress_res.csv', '# Residual stress (DFT-QHA), unit: GPa')
            end if
            if (map%have_fc_quartet) then
                write (lo_iou, *) '... dump FC4 stress to file'
                call write_file_stress_voigt_time(s4t, 'outfile.stress_fc4.csv', '# FC4 stress, unit: GPa')

            end if

            ! dump symmetrized stress for relaxation incl. std. dev. and std. err.
            ! DFT stress
            write (lo_iou, *) '... dump symmetrized DFT stress to file'
            call write_file_stress_voigt(s0_sym, s0_std_sym, s0_err_sym, 'outfile.stress_dft', '# DFT stress (symmetrized), unit: GPa')
            ! 2nd order FC stress
            write (lo_iou, *) '... dump symmetrized FC2 stress to file'
            call write_file_stress_voigt(s2_sym, s2_std_sym, s2_err_sym, 'outfile.stress_fc2', '# FC2 stress (symmetrized), unit: GPa')
            if (map%have_fc_triplet) then
                ! 3rd order FC stress
                write (lo_iou, *) '... dump symmetrized FC3 stress to file'
                call write_file_stress_voigt(s3_sym, s3_std_sym, s3_err_sym, 'outfile.stress_fc3', '# FC3 stress (symmetrized), unit: GPa')
                ! residual stress (DFT - FC2 - FC3)
                write (lo_iou, *) '... dump symmetrized residual stress to file'
                call write_file_stress_voigt(sr_sym, sr_std_sym, sr_err_sym, 'outfile.stress_res', '# Residual stress (DFT-FC) (symmetrized), unit: GPa')
            end if

        end if

    end block pressure_tdep

end if

if (mw%talk) then
    write (*, *) ''
    write (*, *) 'Done in ', tochar(walltime() - t0), 's'
end if
call mw%destroy()

end program
