program anharmonic_free_energy
!!{!src/anharmonic_free_energy/manual.md!}
use konstanter, only: r8,lo_Hartree_to_eV
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

real(r8), dimension(:), allocatable :: temperatures
real(r8) :: timer_init,timer_total
real(r8) :: U0,U1
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

    havehighorder=.true.
    if ( .not.lo_does_file_exist('infile.forceconstant_thirdorder') ) havehighorder=.false.
    if ( .not.lo_does_file_exist('infile.forceconstant_fourthorder') ) havehighorder=.false.

    if ( havehighorder ) then
        call fct%readfromfile(uc,'infile.forceconstant_thirdorder')
        if ( mw%talk ) write(*,*) '... read third order forceconstant'
        call fcf%readfromfile(uc,'infile.forceconstant_fourthorder')
        if ( mw%talk ) write(*,*) '... read fourth order forceconstant'
        call tmr%tock('reading input')
    endif

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
            allocate(temperatures(opts%trangenpts))
            call lo_linspace(opts%trangemin,opts%trangemax,temperatures)
            u=open_file('out','outfile.anharmonic_free_energy')
                do t=1,opts%trangenpts
                    write(u,"(4(1X,E18.12))") temperatures(t),1E9_r8,1E9_r8,1E9_r8
                    write(*,"(4(1X,E18.12))") temperatures(t),1E9_r8,1E9_r8,1E9_r8
                enddo
            close(u)
            write(*,*) 'Found negative eigenvalues. Stopping prematurely.'
        endif
        call mw%destroy()
    endif
    timer_init=walltime()-timer_init
end block init

! We start with the potential energy terms since those are much faster
epotthings: block
    type(lo_energy_differences) :: pot
    type(lo_crystalstructure) :: ss
    type(lo_mdsim) :: sim
    real(r8), dimension(:,:), allocatable :: f2,f3,f4,fp
    real(r8), dimension(:,:), allocatable :: ediff
    real(r8) :: e2,e3,e4,ep
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
        ediff(i,1)=sim%stat%total_energy(i)
        ediff(i,2)=sim%stat%total_energy(i)-e2
        ediff(i,3)=sim%stat%total_energy(i)-e2-ep
        ediff(i,4)=sim%stat%total_energy(i)-e2-ep-e3
        ediff(i,5)=sim%stat%total_energy(i)-e2-ep-e3-e4
    enddo
    call mw%allreduce('sum',ediff)
    ! Make sure it is energy per unit cell
    ediff=ediff*real(uc%na,r8)/real(ss%na,r8)

    if ( mw%talk ) then
        write(*,*) 'Potential energy:'
        write(*,"(1X,A,E21.14,1X,A,F21.14)") '                  input: ',lo_mean(ediff(:,1)),' stddev:',lo_stddev(ediff(:,1))
        write(*,"(1X,A,E21.14,1X,A,F21.14)") '                  E-fc2: ',lo_mean(ediff(:,2)),' stddev:',lo_stddev(ediff(:,2))
        write(*,"(1X,A,E21.14,1X,A,F21.14)") '            E-fc2-polar: ',lo_mean(ediff(:,3)),' stddev:',lo_stddev(ediff(:,3))
        write(*,"(1X,A,E21.14,1X,A,F21.14)") '        E-fc2-polar-fc3: ',lo_mean(ediff(:,4)),' stddev:',lo_stddev(ediff(:,4))
        write(*,"(1X,A,E21.14,1X,A,F21.14)") '    E-fc2-polar-fc3-fc4: ',lo_mean(ediff(:,5)),' stddev:',lo_stddev(ediff(:,5))
    endif
    ! Make a note of the baseline energy
    U0=lo_mean( ediff(:,5) )
    U1=lo_mean( ediff(:,3) )
end block epotthings

! Calculate the actual free energy
getenergy: block
    real(r8) :: f_ph,ah3,ah4

    ! First we get the phonon free energy
    f_ph=dr%phonon_free_energy(opts%temperature)

    if ( havehighorder ) then
        select type(qp); type is(lo_fft_mesh)
            call perturbative_anharmonic_free_energy(uc,fct,fcf,qp,dr,opts%temperature,ah3,ah4,mw,mem,opts%verbosity+1)
        end select
    else
        ah3=0.0_r8
        ah4=0.0_r8
    endif

    if ( mw%talk ) then
        write(*,*) ''
        write(*,*) 'Lowest order approximation to the free energy: (eV/atom)'
        write(*,*) 'Calculated as <U - U_second - U_polar> + F_phonon'
        write(*,*) 'F (eV/atom) = ',(U1 + f_ph)*lo_Hartree_to_eV

        if ( havehighorder ) then
            write(*,*) ''
            write(*,*) 'Free energy with anharmonic corrections: (eV/atom)'
            write(*,*) 'Calculated as <U - U_second - U_polar - U_third - U_fourth> + F_phonon + F_3 + F_4'
            write(*,*) 'F (eV/atom) = ',(U0+f_ph+ah3+ah4)*lo_Hartree_to_eV
        endif
        write(*,*) ''
    endif

end block getenergy

! Kill MPI

call mw%destroy()

end program
