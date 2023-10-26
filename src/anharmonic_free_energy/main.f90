program anharmonic_free_energy
!!{!src/anharmonic_free_energy/manual.md!}
use konstanter, only: r8,lo_Hartree_to_eV,lo_kb_Hartree,lo_pressure_HartreeBohr_to_GPa
use gottochblandat, only: open_file,walltime,lo_linspace,lo_progressbar_init,lo_progressbar,tochar,&
    lo_does_file_exist,lo_mean,lo_stddev
use mpi_wrappers, only: lo_mpi_helper
use lo_memtracker, only: lo_mem_helper
use lo_timetracker, only: lo_timer
use type_crystalstructure, only: lo_crystalstructure
use type_forceconstant_secondorder, only: lo_forceconstant_secondorder
use type_forceconstant_thirdorder, only: lo_forceconstant_thirdorder
use type_forceconstant_fourthorder, only: lo_forceconstant_fourthorder
use type_mdsim, only: lo_mdsim
use type_qpointmesh, only: lo_qpoint_mesh,lo_generate_qmesh,lo_read_qmesh_from_file,lo_fft_mesh
use type_phonon_dispersions, only: lo_phonon_dispersions
use lo_phonon_bandstructure_on_path, only: lo_phonon_bandstructure
use type_phonon_dos, only: lo_phonon_dos

use options, only: lo_opts
use energy, only: perturbative_anharmonic_free_energy
use epot, only: lo_energy_differences

implicit none

type(lo_opts) :: opts
type(lo_forceconstant_secondorder) :: fc
type(lo_forceconstant_thirdorder) :: fct
type(lo_forceconstant_fourthorder) :: fcf
type(lo_phonon_dispersions) :: dr
type(lo_crystalstructure) :: uc
class(lo_qpoint_mesh), allocatable :: qp
type(lo_mpi_helper) :: mw
type(lo_mem_helper) :: mem
type(lo_timer) :: tmr
type(lo_mdsim) :: sim

real(r8), dimension(3,5) :: cumulant
real(r8) :: timer_init,timer_total
real(r8) :: U0,U1,energy_unit_factor
logical :: havehighorder

! Init MPI, timers and options
call mw%init()
timer_total=walltime()
timer_init=walltime()
call opts%parse()
call tmr%start()

! set up all possible meshes and get all dispersions
init: block
    integer :: t,u

    ! only be verbose on the first rank
    if ( .not. mw%talk ) opts%verbosity=-10

    ! Read structure
    call uc%readfromfile('infile.ucposcar',verbosity=opts%verbosity)
    call uc%classify('wedge',timereversal=.true.)
    if ( mw%talk ) write(*,*) '... read unitcell'
    if ( opts%readiso ) then
        if ( mw%talk ) write(*,*) '... reading isotope distribution from file'
        call uc%readisotopefromfile()
    endif

    ! Read forceconstants
    call fc%readfromfile(uc,'infile.forceconstant',mem,verbosity=-1)
    if ( mw%talk ) write(*,*) '... read second order forceconstant'

    if (opts%thirdorder) call fct%readfromfile(uc,'infile.forceconstant_thirdorder')
    if (opts%fourthorder) call fcf%readfromfile(uc,'infile.forceconstant_fourthorder')
    havehighorder=.false.
    if (opts%thirdorder .or. opts%fourthorder) havehighorder=.true.

    call tmr%tock('reading input')

    ! Get a q-mesh for the integrations. Always an FFT mesh.
    if ( mw%talk ) write(*,*) '... generating q-mesh'
    call lo_generate_qmesh(qp,uc,opts%qgrid,'fft',timereversal=.true.,headrankonly=.false.,mw=mw,mem=mem,verbosity=opts%verbosity)

    ! Dispersions for everyone!
    call dr%generate(qp,fc,uc,mw=mw,mem=mem,verbosity=opts%verbosity)
    if ( mw%talk ) write(*,*) '... got the full dispersion relations'
    call tmr%tock('harmonic properties')

    ! Check for imaginary modes right away
    if ( dr%omega_min .lt. 0.0_r8 ) then
        ! Dump the free energies
        if ( mw%talk ) then
             write(*,*) 'Found negative eigenvalues. Stopping prematurely since no free energy can be defined.'
        endif
        call mw%destroy()
        stop
    endif
    timer_init=walltime()-timer_init

end block init

! We start with the potential energy terms since those are much faster
epotthings: block
    type(lo_energy_differences) :: pot
    type(lo_crystalstructure) :: ss

    real(r8), dimension(:,:), allocatable :: f2,f3,f4,fp
    real(r8), dimension(:,:), allocatable :: ediff

    real(r8) :: e2,e3,e4,ep,inverse_kbt
    integer :: i

    call ss%readfromfile('infile.ssposcar')
    call ss%classify('supercell',uc)
    call pot%setup(uc,ss,fc,fct,fcf,mw,opts%verbosity+1)

    ! Fetch simulation data from file
    if (lo_does_file_exist('infile.sim.hdf5')) then
        call sim%read_from_hdf5('infile.sim.hdf5', verbosity=opts%verbosity + 2, stride=-1)
    else
        call sim%read_from_file(verbosity=opts%verbosity + 2, stride=1, magnetic=.false., dielectric=.false., nrand=-1, mw=mw)
    end if

    allocate(f2(3,ss%na))
    allocate(f3(3,ss%na))
    allocate(f4(3,ss%na))
    allocate(fp(3,ss%na))

    ! Calculate the baseline energy
    allocate(ediff(sim%nt,5))
    ediff=0.0_r8

    do i=1,sim%nt
        if ( mod(i,mw%n) .ne. mw%r ) cycle
        call pot%energies_and_forces(sim%u(:,:,i),e2,e3,e4,ep,f2,f3,f4,fp)
        ediff(i,1)=sim%stat%potential_energy(i)
        ediff(i,2)=sim%stat%potential_energy(i)-e2
        ediff(i,3)=sim%stat%potential_energy(i)-e2-ep
        ediff(i,4)=sim%stat%potential_energy(i)-e2-ep-e3
        ediff(i,5)=sim%stat%potential_energy(i)-e2-ep-e3-e4
    enddo
    call mw%allreduce('sum',ediff)

    if ( sim%temperature_thermostat .gt. 1E-5_r8 ) then
        inverse_kbt=1.0_r8/lo_kb_Hartree/sim%temperature_thermostat
    else
        inverse_kbt=0.0_r8
    endif

    ! Compute the first and second order cumulants
    do i=1,5
        cumulant(1,i)=lo_mean(ediff(:,i))
        cumulant(2,i)=lo_mean((ediff(:,i)-cumulant(1,i))**2)
        cumulant(2,i)=cumulant(2,i)*inverse_kbt*0.5_r8
        cumulant(3,i)=lo_mean((ediff(:,i)-cumulant(1,i))**3)
        cumulant(3,i)=cumulant(3,i)*inverse_kbt**2/6.0_r8
    enddo

    ! And normalize it to be per atom
    cumulant=cumulant/real(ss%na,r8)


    if ( mw%talk ) then
        write(*,*) 'Temperature (K) (from infile.meta): ',sim%temperature_thermostat
        write(*,*) 'Potential energy:'
        write(*,"(1X,A,E21.14,1X,A,F21.14)") '                  input: ',cumulant(1,1)*lo_Hartree_to_eV,' upper bound:',cumulant(2,1)*lo_Hartree_to_eV
        write(*,"(1X,A,E21.14,1X,A,F21.14)") '                  E-fc2: ',cumulant(1,2)*lo_Hartree_to_eV,' upper bound:',cumulant(2,2)*lo_Hartree_to_eV
        write(*,"(1X,A,E21.14,1X,A,F21.14)") '            E-fc2-polar: ',cumulant(1,3)*lo_Hartree_to_eV,' upper bound:',cumulant(2,3)*lo_Hartree_to_eV
        if (opts%thirdorder) then
            write(*,"(1X,A,E21.14,1X,A,F21.14)") '        E-fc2-polar-fc3: ',cumulant(1,4)*lo_Hartree_to_eV,' upper bound:',cumulant(2,4)*lo_Hartree_to_eV
        end if
        if (opts%fourthorder) then
            write(*,"(1X,A,E21.14,1X,A,F21.14)") '    E-fc2-polar-fc3-fc4: ',cumulant(1,5)*lo_Hartree_to_eV,' upper bound:',cumulant(2,5)*lo_Hartree_to_eV
        end if
    endif
end block epotthings

! Calculate the actual free energy
getenergy: block
    real(r8) :: f_ph,ah3,ah4,fe2_1,fe2_2,fe3_1,fe3_2,fe4_1,fe4_2,pref
    integer :: u

    ! Some heuristics to figure out what the temperature is.
    if (opts%quantum) then
        f_ph=dr%phonon_free_energy(sim%temperature_thermostat)
    else
        f_ph=dr%phonon_free_energy_classical(sim%temperature_thermostat)
    end if

    if ( havehighorder ) then
        select type(qp); type is(lo_fft_mesh)
            call perturbative_anharmonic_free_energy(uc,fct,fcf,qp,dr,sim%temperature_thermostat,ah3,ah4,&
                                                     opts%fourthorder,opts%quantum,mw,mem,opts%verbosity+1)
        end select
    else
        ah3=0.0_r8
        ah4=0.0_r8
    endif

    if (opts%stochastic) then
        pref = -1.0_r8
    else
        pref = 1.0_r8
    end if
    if ( mw%talk ) then
        ! Write on a file
        u = open_file('out', 'outfile.anharmonic_free_energy')
        fe2_1 = (cumulant(1,3) + f_ph)*lo_Hartree_to_eV
        fe2_2 = (cumulant(1,3) + f_ph + pref * cumulant(2,3))*lo_Hartree_to_eV
        write(u, "(1X, A17, F8.2, 1X, A17)") '# Free energy at ', sim%temperature_thermostat, 'K, unit : eV/atom'
        write(u, *) '# Lowest order (1st order cumulant, 2nd order cumulant)'
        write(u, "(1X, 2(F12.5,' '))") fe2_1, fe2_2
        write(*,*) ''
        write(*,*) 'Lowest order approximation to the free energy: (eV/atom)'
        write(*,*) 'Calculated as <U - U_second - U_polar> + F_phonon'
        write(*,*) 'F (eV/atom) = ', fe2_1
        write(*,*) 'F_phonon (eV/atom) =', f_ph*lo_Hartree_to_eV
        write(*,*) 'Second order cumulant (meV/atom):',cumulant(2,3)*lo_Hartree_to_eV*1000
        write(*,*) 'Third order cumulant (meV/atom):',cumulant(3,3)*lo_Hartree_to_eV*1000
        if (opts%thirdorder .or. opts%fourthorder) then
            fe3_1 = (cumulant(1,4) + f_ph + ah3)*lo_Hartree_to_eV
            fe3_2 = (cumulant(1,4) + f_ph + ah3 + pref * cumulant(2,4))*lo_Hartree_to_eV
            write(u, *) '# Third order anharmonic corrections (1st order cumulant, 2nd order cumulant)'
            write(u, "(1X, 2(F12.5,' '))") fe3_1, fe3_2
            write(*,*) ''
            write(*,*) 'Free energy with third order anharmonic corrections: (eV/atom)'
            write(*,*) 'Calculated as <U - U_second - U_polar - U_third> + F_phonon + F_3'
            write(*,*) 'F (eV/atom) = ', fe3_1
            write(*,*) 'F_3 (eV/atom) = ', ah3*lo_Hartree_to_eV
            write(*,*) 'Second order cumulant (meV/atom):',cumulant(2,4)*lo_Hartree_to_eV*1000
            write(*,*) 'Third order cumulant (meV/atom):',cumulant(3,4)*lo_Hartree_to_eV*1000
        end if
        if (opts%fourthorder) then
            fe4_1 = (cumulant(1,5) + f_ph + ah3 + ah4)*lo_Hartree_to_eV
            fe4_2 = (cumulant(1,5) + f_ph + ah3 + ah4 + pref * cumulant(2,5))*lo_Hartree_to_eV
            write(u, *) '# Fourth order anharmonic corrections (1st order cumulant, 2nd order cumulant)'
            write(u, "(1X, 2(F12.5,' '))") fe4_1, fe4_2
            write(*,*) ''
            write(*,*) 'Free energy with fourth order anharmonic corrections: (eV/atom)'
            write(*,*) 'Calculated as <U - U_second - U_polar - U_third - U_fourth> + F_phonon + F_3 + F_4'
            write(*,*) 'F (eV/atom) = ',fe4_1
            write(*,*) 'F_3 (eV/atom) = ', ah3*lo_Hartree_to_eV
            write(*,*) 'F_4 (eV/atom) = ', ah4*lo_Hartree_to_eV
            write(*,*) 'Second order cumulant (meV/atom):',cumulant(2,5)*lo_Hartree_to_eV*1000
            write(*,*) 'Third order cumulant (meV/atom):',cumulant(3,5)*lo_Hartree_to_eV*1000
        end if
    endif

end block getenergy

! Kill MPI
call mw%destroy()

end program
