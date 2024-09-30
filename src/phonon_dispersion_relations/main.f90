program phonon_dispersion_relations
!!{!src/phonon_dispersion_relations/manual.md!}
use konstanter, only: r8, lo_iou, lo_tol, lo_sqtol, lo_hartree_to_eV, &
                      lo_exitcode_param, lo_frequency_Hartree_to_THz
use gottochblandat, only: walltime, tochar, lo_chop, lo_mean, lo_kmesh_density, &
                          open_file, lo_linspace, lo_frobnorm
use mpi_wrappers, only: lo_mpi_helper, lo_stop_gracefully
use lo_memtracker, only: lo_mem_helper
use dump_data, only: lo_dump_gnuplot_2d_real
use options, only: lo_opts
use type_crystalstructure, only: lo_crystalstructure
use type_forceconstant_secondorder, only: lo_forceconstant_secondorder
use type_forceconstant_thirdorder, only: lo_forceconstant_thirdorder
use type_qpointmesh, only: lo_qpoint_mesh, lo_generate_qmesh, lo_read_qmesh_from_file
use type_phonon_dispersions, only: lo_phonon_dispersions
use type_phonon_dos, only: lo_phonon_dos
use lo_phonon_bandstructure_on_path, only: lo_phonon_bandstructure
use type_mdsim, only: lo_mdsim
! local
use velocitydos, only: calculate_everything, lo_calculate_U0
use densityplots, only: inelastic_spectra
use activity, only: estimate_activity

implicit none
! for normal phonon things
type(lo_opts) :: opts
type(lo_mpi_helper) :: mw
type(lo_mem_helper) :: mem
type(lo_forceconstant_secondorder) :: fc
type(lo_forceconstant_thirdorder) :: fct
type(lo_crystalstructure) :: uc
type(lo_phonon_bandstructure) :: bs
class(lo_qpoint_mesh), allocatable :: qp
type(lo_phonon_dispersions) :: dr
type(lo_phonon_dos) :: pd
real(r8) :: timer

! Read things from file
init: block
    ! Start timer
    timer = walltime()
    ! init MPI
    call mw%init()
    ! Get the command line arguments
    call opts%parse()
    if (mw%talk .eqv. .false.) opts%verbosity = -10
    ! Start memory tracker
    call mem%init()

    ! Read the crystal structure and make sure I have all the symmetry information
    call uc%readfromfile('infile.ucposcar', verbosity=opts%verbosity)
    call uc%classify('wedge', timereversal=opts%timereversal)
    ! Maybe non-default isotope distribution
    if (opts%readiso) then
        if (mw%talk) write (lo_iou, *) '... reading isotope distribution from file'
        call uc%readisotopefromfile()
    end if
    ! Read the force constant
    call fc%readfromfile(uc, 'infile.forceconstant', mem, opts%verbosity)
    ! Perhaps some Gruneisen parameters, in that case we need third order force constants.
    if (opts%gruneisen) then
        call fct%readfromfile(uc, 'infile.forceconstant_thirdorder')
    end if
end block init

! Dispersions along a path?
path: block
    real(r8), dimension(3) :: v1
    real(r8) :: t0
    integer :: i, verb
    character(len=1000) :: opf

    ! Get the dispersions on the path
    t0 = walltime()
    if (uc%na .ge. 4 .and. mw%talk) then
        verb = opts%verbosity + 2
    else
        verb = opts%verbosity
    end if

    if (opts%gruneisen) then
        call bs%generate(uc, fc, timereversal=opts%timereversal, mw=mw, mem=mem, verbosity=verb, &
                         fct=fct, npts=opts%nq, readpathfromfile=opts%readpathfromfile)
    else
        call bs%generate(uc, fc, timereversal=opts%timereversal, mw=mw, mem=mem, verbosity=verb, &
                         npts=opts%nq, readpathfromfile=opts%readpathfromfile)
    end if

    ! Write a small summary
    if (mw%talk) then
        write (lo_iou, *) ''
        write (lo_iou, *) ' Settings: '
        write (lo_iou, *) '      mode Gruneisen parameters: ', opts%gruneisen
        write (lo_iou, *) '                number of bands: ', tochar(bs%n_mode)
        write (lo_iou, *) '                number of paths: ', tochar(bs%n_path)
        write (lo_iou, *) '  number of points on each path: ', tochar(bs%n_point_per_path)
        write (lo_iou, *) ''
        do i = 1, bs%n_path
            write (lo_iou, *) '        path '//tochar(i)//': from '//trim(bs%symb_q_start(i))//' to '//trim(bs%symb_q_end(i))
            opf = "(1X,'  starting point: cartesian',3(F9.6,' '),'fractional',3(F9.6,' '))"
            v1 = lo_chop(uc%cartesian_to_fractional(bs%segment(i)%r1, reciprocal=.true., pbc=.false.), lo_sqtol)
            write (lo_iou, opf) bs%segment(i)%r1, v1
            opf = "(1X,'    ending point: cartesian',3(F9.6,' '),'fractional',3(F9.6,' '))"
            v1 = lo_chop(uc%cartesian_to_fractional(bs%segment(i)%r2, reciprocal=.true., pbc=.false.), lo_sqtol)
            write (lo_iou, opf) bs%segment(i)%r2, v1
        end do
        write (lo_iou, *) ' '
        write (lo_iou, "(' ...            got the phonon bandstructure in ',F12.9,'s' )") walltime() - t0

        ! Dump things to files for plotting
        if (opts%sortbands) call bs%sort_bands(mem, opts%verbosity)

        call bs%write_dispersive_property(opts%enhet, 'frequency', 'outfile.dispersion_relations', opts%pdfplot)
        call bs%write_dispersive_property(opts%enhet, 'groupvelocity', 'outfile.group_velocities', opts%pdfplot)
        if (opts%gruneisen) then
            call bs%write_dispersive_property(opts%enhet, 'gruneisen', 'outfile.mode_gruneisen_parameters', .false.)
        end if
    end if
end block path

! Determine Raman/IR activity and report
activity_raman_ir: block
    integer :: iq, jq, imode, u, active_raman, active_ir, n_raman, n_ir
    real(r8) :: norm
    character(len=*), parameter:: file = 'outfile.mode_activity.csv'

    ! Now estimate Raman and IR activity, always neat I suppose.
    call estimate_activity(uc, bs, mw, mem, opts%verbosity)

    u = open_file('out', file)
    write (u, '(a,I3,a)') '# Raman activity for ', bs%n_active_qpoint, ' active q-points'
    write (u, '(a)') '# iq: q-point index on path'
    write (u, '(a)') '# imode: mode index'
    write (u, '(a)') '# frequency: mode frequency in THz'
    write (u, '(a)') '# active_raman: Raman active yes/no (1/0)'
    write (u, '(a)') '# active_ir: IR active yes/no (1/0)'
    write (u, '(a)') 'iq,imode,frequency,active_raman,active_ir'

    do jq = 1, bs%n_active_qpoint
        n_ir = 0
        n_raman = 0
        iq = bs%active_qpoint(jq)%ind_path
        do imode = 1, size(bs%p(iq)%omega(:))

            ! Raman
            active_raman = 0
            if (allocated(bs%active_qpoint(jq)%coeff_Raman_mode)) then
                norm = lo_frobnorm(bs%active_qpoint(jq)%coeff_Raman_mode(:, :, imode))
                if (norm > lo_tol) active_raman = 1
            end if

            ! IR
            active_ir = 0
            if (allocated(bs%active_qpoint(jq)%coeff_IR_mode)) then
                norm = lo_frobnorm(bs%active_qpoint(jq)%coeff_IR_mode(:, :, imode))
                if (norm > lo_tol) active_ir = 1
            end if

            ! dump
            write (u, '(2(I5,","),E15.8,2(",",I5))') &
                iq, imode, bs%p(iq)%omega(imode)*lo_frequency_Hartree_to_THz, active_raman, active_ir

            ! count
            n_ir = n_ir + active_ir
            n_raman = n_raman + active_raman
        end do
    end do
    close (u)

    ! report
    if (mw%talk) then
        write (lo_iou, *) '... found ', bs%n_active_qpoint, ' active q-points'
        write (lo_iou, *) '... found ', n_raman, ' Raman active modes'
        write (lo_iou, *) '... found ', n_ir, ' IR active modes'
        write (lo_iou, *) '... activity written to ', file
    end if

end block activity_raman_ir

! Get the full mesh and the density of states?
dos: block
    type(lo_mdsim) :: sim
    real(r8), dimension(:), allocatable :: temperatures
    real(r8) :: f0, f1, f2, f3, t0, U0
    integer :: i, u, verb
    ! The dispersions along the path are done. Perhaps we want the full dispersion relations for something?
    if (opts%fullmesh) then

        ! Adjust the number of q-points in a neat way
        if (opts%dumpsecret) then
            ! This should be a mesh equivalent to a 34x34x34 mesh for fcc.
            opts%meshtype = 2
            opts%qgrid = lo_kmesh_density(uc%reciprocal_latticevectors, uc%na, 18)
        end if

        if (mw%talk) then
            write (lo_iou, *) ' '
            write (lo_iou, *) 'Calculating phonon dispersion across the Brillouin zone'
        end if

        ! Get a q-point mesh
        if (opts%readqmesh) then
            call lo_read_qmesh_from_file(qp, uc, 'infile.qgrid.hdf5', mem, verbosity=opts%verbosity)
        else
            select case (opts%meshtype)
            case (1)
                call lo_generate_qmesh(qp, uc, opts%qgrid, 'monkhorst', timereversal=opts%timereversal, &
                                       headrankonly=.false., mw=mw, mem=mem, verbosity=opts%verbosity)
            case (2)
                call lo_generate_qmesh(qp, uc, opts%qgrid, 'fft', timereversal=opts%timereversal, &
                                       headrankonly=.false., mw=mw, mem=mem, verbosity=opts%verbosity)
            case (3)
                call lo_generate_qmesh(qp, uc, opts%qgrid, 'wedge', timereversal=opts%timereversal, &
                                       headrankonly=.false., mw=mw, mem=mem, verbosity=opts%verbosity)
            case default
                call lo_stop_gracefully(['Unknown mesh type'], lo_exitcode_param, __FILE__, __LINE__, mw%comm)
            end select
        end if

        ! Get the full dispersion relations
        t0 = walltime()
        if (uc%na .ge. 4 .and. mw%talk) then
            verb = opts%verbosity + 2
        else
            verb = opts%verbosity
        end if
        call dr%generate(qp, fc, uc, mw=mw, mem=mem, verbosity=verb)
    end if

    ! Perhaps we want the phonon density of states and related things
    if (opts%dos) then
        t0 = walltime()
        if (mw%talk .eqv. .false.) then
            verb = 0
        else
            verb = max(1, opts%verbosity)
        end if
        call pd%generate(dr, qp, uc, mw, mem, verbosity=verb, sigma=opts%sigma, &
                         n_dos_point=opts%dospoints, integrationtype=opts%dosintegrationtype)
        ! Dump it to file
        if (mw%talk) then
            call pd%write_to_file(uc, opts%enhet, 'outfile.phonon_dos')
            call pd%write_to_hdf5(uc, opts%enhet, 'outfile.phonon_dos.hdf5', mem)
        end if

        if (opts%U0) then
            ! Read simulation data
            call sim%read_from_file(verbosity=opts%verbosity + 1, stride=1, &
                                    magnetic=.false., dielectric=.false., nrand=-1, mw=mw)
            call lo_calculate_U0(sim, fc, uc, sim%crystalstructure, U0, mw, opts%verbosity + 1)
        end if

        ! Get the density-plot thing
        !if ( mw%talk ) call omega_vs_norm_q_density(qp,dr,uc)
    end if

    ! A positive number of temperatures means I should dump energy. I do this not in parallel.
    if (opts%trangenpts .gt. 0 .and. mw%talk) then
        if (.not. opts%fullmesh) then
            write (*, *) 'I screwed up the heuristics'
            stop
        end if
        write (*, *) '... calculating free energy'

        allocate (temperatures(opts%trangenpts))
        call lo_linspace(opts%trangemin, opts%trangemax, temperatures)
        ! dump some free energies and entropies
        write (*, *) ''
        write (*, *) '      T(K)     F(eV/atom)         S(eV/K/atom)       Cv(eV/K/atom)'
        u = open_file('out', 'outfile.free_energy')
        do i = 1, opts%trangenpts
            f0 = temperatures(i)
            f1 = dr%phonon_free_energy(f0)*lo_hartree_to_eV
            f2 = dr%phonon_entropy(f0)*lo_hartree_to_eV
            f3 = dr%phonon_cv(f0)*lo_hartree_to_eV
            U0 = U0*lo_hartree_to_eV
            if (abs(f1 - 123456789.0_r8) .lt. lo_tol) then
                ! return NaN if I have imaginary frequencies
                write (u, *) f0, ' NaN  NaN  NaN'
                write (*, *) f0, ' NaN  NaN  NaN'
            else
                if (opts%U0) then
                    write (u, "(1X,F12.5,5(1X,E18.11))") temperatures(i), f1, f2, f3, U0 !f1+U0
                    write (*, "(1X,F12.5,5(1X,E18.11))") temperatures(i), f1, f2, f3, f1 + U0
                else
                    write (u, "(1X,F12.5,4(1X,E18.11))") temperatures(i), f1, f2, f3
                    write (*, "(1X,F12.5,4(1X,E18.11))") temperatures(i), f1, f2, f3
                end if
            end if
        end do
        close (u)
    end if
end block dos

! Maybe I want the fancy unfolded bandstructure
if (opts%unfold) then
    unfold: block
        !call unfold_bandstructure(uc,fc,bs,mw,opts,opts%verbosity+1)
        !call unfold_brutal(uc,fc,bs,mw,opts,opts%verbosity+1)
    end block unfold
end if

! Maybe get the inelastic scattering spectra
if (opts%inelastic) then
    call inelastic_spectra(bs, uc, fc, qp, dr, mw, mem, opts%verbosity + 1)
    if (mw%r .eq. 0) then
        call bs%write_spectral_function_to_hdf5(opts%enhet, 'outfile.phonon_spectral_function.hdf5')
    end if
end if

! Print some file and wrap it up
wrapup: block
    character(len=1000) :: buf
    ! Write everything on the grid to file, for some inexplicable reason.
    if (opts%dumpgrid) then
        ! small sanity check
        if (.not. opts%fullmesh) then
            write (*, *) 'I screwed up the heuristics'
            stop
        end if
        if (opts%gruneisen) then
            !call dr%gruneisen(qp,fct,uc,opts%verbosity)
        end if
        ! actually write it
        if (mw%talk) then
            buf = 'outfile.grid_dispersions_irreducible.hdf5'
            call dr%write_irreducible_to_hdf5(qp, uc, trim(buf), mem)
            write (*, *) "... wrote data for the irreducible grid to '", trim(buf), "'"

            buf = 'outfile.grid_dispersions.hdf5'
            call dr%write_to_hdf5(qp, uc, trim(buf), mem)
            write (*, *) "... wrote data for the full grid to        '", trim(buf), "'"
        end if
    end if

    ! Dump everything to hdf5
    if (mw%talk) then
        if (opts%sortbands) call bs%sort_bands(mem, opts%verbosity)
        call bs%write_to_hdf5(uc, opts%enhet, 'outfile.dispersion_relations.hdf5', mem)
    end if

    ! And finally dump secret things, if needed
    if (opts%dumpsecret) then
        call calculate_everything(uc, bs, dr, pd, qp, fc, mw, mem)
    end if

    ! if ( mw%talk ) write(lo_iou,*) 'Memory usage:'
    ! call memdump('q-mesh'              ,qp%size_in_mem(qp),mw,lo_iou)
    ! call memdump('phonon bandstructure',bs%size_in_mem(),mw,lo_iou)

    if (mw%talk) then
        write (*, *) ' '
        write (*, *) 'All done in ', tochar(walltime() - timer), 's'
    end if
    call mw%destroy()
end block wrapup

contains

! !> dump a message regarding memory usage
! subroutine memdump(obj,bytes,mw,iou)
!     !> what are we measuring the size of
!     character(len=*), intent(in) :: obj
!     !> size of object, on this rank, in bytes
!     integer, intent(in) :: bytes
!     !> MPI helper
!     type(lo_mpi_helper), intent(inout) :: mw
!     !> unit to write to
!     integer, intent(in) :: iou
!
!     real(r8), parameter :: toMiB=1.0_r8/1024_r8**2
!     real(r8), dimension(mw%n) :: mb
!
!     mb=0.0_r8
!     mb(mw%r+1)=bytes*toMiB
!     call mw%allreduce('sum',mb)
!     if ( mw%talk ) then
!         write(iou,"(1X,A30,':',3(1X,F20.6))") trim(adjustl(obj)),sum(mb)/mw%n,maxval(mb),minval(mb)
!     endif
! end subroutine memdump

end program
