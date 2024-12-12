#include "precompilerdefinitions"
program average_structure
!!{!holding/average_structure/manual.md!}
use konstanter !, only: flyt,lo_status,lo_volume_A_to_Bohr,lo_volume_bohr_to_A,lo_eV_to_Hartree,lo_Hartree_to_eV,lo_pressure_HartreeBohr_to_GPa
use gottochblandat
use options, only: lo_opts
use type_crystalstructure, only: lo_crystalstructure
use type_irredcell !, only: lo_irredcell
use type_mdsim, only: lo_mdsim
use type_symmetryoperation, only: lo_symset, lo_operate_on_vector
implicit none
! for normal phonon things
type(lo_opts) :: opts
type(lo_crystalstructure) :: uc, ss, ur, nuc, nss
type(lo_irredcell) :: ic
type(lo_mdsim) :: sim
real(flyt), dimension(:), allocatable :: volume, energy, eta, pressure, bulkmodulus
real(flyt), dimension(:), allocatable :: energy_model
integer :: u, i, ne
character(len=1000) :: opf

init: block
    integer :: u
    integer :: i, j, k
    real(flyt), dimension(6) :: sigma
    character(len=1000) :: opf
    real(flyt), dimension(3, 3) :: m0
    ! Get the command line arguments
    call opts%parse()

    ! Get the unitcell that determines the symmetry
    call uc%readfromfile('infile.ucposcar')
    ! Get the unitcell that defines reference positions
    call ss%readfromfile('infile.ssposcar')
    call ss%classify('supercell', uc)

    ! Reference positions to determine displacement from
    call ur%readfromfile('infile.refposcar')

    ! Read the MD simulation
    call sim%read_from_file(verbosity=2, stride=opts%stride)
!    lo_deallocate(sim%lattice)
    ! Remove center of mass drift
!    call sim%remove_force_and_center_of_mass_drift()

    call ic%generate(uc, verbosity=2, reference_positions=ur%r)
    call avg_structure_nvt(uc, ss, ic, sim, nuc, nss, ur)

!stop
!    call sim%remove_force_and_center_of_mass_drift()
!    ! Remove arbitrary rotations
!    nuc=uc
!    nss=ss
!    do i=1,1 !20
!        m0=nuc%latticevectors
!        ! Get the symmetric representation
!        if ( i .eq. 1 ) then
!            call ic%generate(nuc,verbosity=2)
!        else
!            call ic%generate(nuc,verbosity=0)
!        endif
!        !
!!        call remove_rotations(nuc,nss,ic,sim)
!        ! Get the averaged structure
!        uc=nuc
!        ss=nss
!        call avg_structure(uc,ss,ic,sim,nuc,nss)
!        !
!        write(*,*) i,sum(abs(m0-nuc%latticevectors))
!        if ( sum(abs(m0-nuc%latticevectors)) .lt. 1E-10_flyt ) exit
!    enddo

    ! Dump the new, clean averaged structure
    call nuc%writetofile('outfile.ucposcar', 1)
    call nss%writetofile('outfile.ssposcar', 1)

    write (*, *) '... wrote averaged structure'

    stop

!    u=open_file('out','outfile.meta')
!        write(u,*) nss%na
!        write(u,*) sim%nt
!        write(u,*) sim%ts
!        write(u,*) sim%t_thermostat
!    close(u)
!    write(*,*) '... wrote new stat'
!
!    u=open_file('out','outfile.positions')
!        do i=1,sim%nt
!        do j=1,sim%na
!            write(u,*) sim%r(:,j,i)
!        enddo
!        enddo
!    close(u)
!    write(*,*) '... wrote new positions'
!
!    u=open_file('out','outfile.forces')
!        do i=1,sim%nt
!        do j=1,sim%na
!            write(u,*) sim%f(:,j,i)*lo_force_HartreeBohr_to_eVA
!        enddo
!        enddo
!    close(u)
!    write(*,*) '... wrote new forces'
!
!    u=open_file('out','outfile.stat')
!        do i=1,sim%nt
!            sigma(1)=sim%stat%stress(1,1,i)
!            sigma(2)=sim%stat%stress(2,2,i)
!            sigma(3)=sim%stat%stress(3,3,i)
!            sigma(4)=sim%stat%stress(1,2,i)
!            sigma(5)=sim%stat%stress(2,3,i)
!            sigma(6)=sim%stat%stress(1,3,i)
!            sigma=sigma*lo_pressure_HartreeBohr_to_Gpa
!            write(u,"(I6,' ',13(2X,EN16.7))") i,i*sim%ts,sim%stat%et(i)*lo_Hartree_to_eV,sim%stat%ep(i)*lo_Hartree_to_eV,&
!                                              sim%stat%ek(i)*lo_Hartree_to_eV,sim%stat%t(i),sim%stat%p(i)*lo_pressure_HartreeBohr_to_Gpa,sigma
!        enddo
!    close(u)
!
!    write(*,*) '... wrote new stat'
!
!    if ( allocated(sim%lattice) ) then
!        u=open_file('out','outfile.lattice')
!            do i=1,sim%nt
!            do j=1,3
!                write(u,*) sim%lattice(:,j,i)*lo_Bohr_to_A
!            enddo
!            enddo
!        close(u)
!        write(*,*) '... wrote new lattice'
!    endif

end block init

! check: block
!     type(lo_symset) :: sym
!     real(flyt), dimension(:,:), allocatable :: r0,r1,r2
!     real(flyt), dimension(:), allocatable :: lv0,lv1,x0,x1
!     real(flyt), dimension(3) :: v0,v1
!     real(flyt) :: f0,f1,f2
!     integer, dimension(:), allocatable :: di
!     integer :: i,j,k,l,o,ii,jj
!     integer :: ix
!
! !    call sym%generate(basis=uc%latticevectors,timereversal=.true.,verbosity=1)
!
!
!     lo_allocate(r0(3,uc%na))
!     lo_allocate(di(uc%na))
!     r0=lo_huge
!     di=1
!     do i=1,uc%na
!         if ( di(i) .eq. 0 ) cycle
!         ! round the coordinates of this atom
!         do j=1,3
!             v0(j)=closest_fraction(uc%r(j,i))
!         enddo
! write(*,*) i,v0
!         do o=1,uc%sym%n
!             v1=lo_operate_on_vector(uc%sym%op(o),v0,fractional=.true.)
!             r0(:,uc%sym%op(o)%fmap(i))=v1
!             di(uc%sym%op(o)%fmap(i))=0
! write(*,*) i,o,v1
!         enddo
!     enddo
!
! write(*,*)
! do i=1,uc%na
!     write(*,*) r0(:,i)
! enddo
!
!
!
! stop
!     lo_allocate(r0(3,uc%na))
!     lo_allocate(r1(3,uc%na))
!     lo_allocate(r2(3,uc%na))
!     lo_allocate(lv0(3*uc%na))
!     lo_allocate(x0(ic%nx_internal))
!     lo_allocate(x1(ic%nx_internal))
!     ! Which degree of freedom?
!     ix=1
!     ! Get the displacement pattern for this degree of freedom
!     do i=1,uc%na
!         lv0( (i-1)*3+1:i*3 )=uc%r(:,i)
!     enddo
!     x0=lo_clean_fractional_coordinates(matmul(transpose(ic%coeffM_internal),lv0) )
!     x0=matmul(transpose(ic%coeffM_internal),lv0)
!     do i=1,ic%nx_internal
!         x1(i)=closest_fraction(x0(i))
! write(*,*) i,x0(i),x1(i)
!     enddo
!
!     lv0=lo_clean_fractional_coordinates(matmul(ic%coeffM_internal,x1))
!     lv0=matmul(ic%coeffM_internal,x0)
!
!     do i=1,uc%na
!         r1(:,i)=lv0( (i-1)*3+1:i*3 )
!         write(*,*) i
!         write(*,*) uc%r(:,i)
!         write(*,*) r1(:,i)
!         write(*,*) uc%r(:,i)+r1(:,i)
!         write(*,*) uc%r(:,i)-r1(:,i)
!     enddo
!
! !write(*,*) shape(ic%coeffM_
!
!
! write(*,*) matmul(transpose(ic%coeffM_internal),lv0)
!
!     ! See if there is a nice one that fits!
!     f0=lo_huge
!     do ii=0,6
!     do jj=1,6
!
!     enddo
!     enddo
!
! end block check

contains

function closest_fraction(x) result(y)
    real(flyt) :: x, y

    integer, parameter :: nmax = 4
    real(flyt) :: f0, f1
    integer :: i, j
    f0 = lo_huge

    do i = 0, nmax
    do j = 1, nmax
        f1 = real(i, flyt)/real(j, flyt)
        if (abs(x - f1) .lt. f0) then
            f0 = abs(x - f1)
            y = f1
        end if
    end do
    end do
!    y=x
end function

end program

