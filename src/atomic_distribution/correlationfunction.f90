#include "precompilerdefinitions"
module correlationfunction
use konstanter, only: flyt
!use dump_data
!use gottochblandat
!use type_crystalstructure
!use type_mdsim
!use pairmapping
!use hdf5_wrappers
!use hdf5
!use fftw_wrappers

implicit none

private
public :: lo_correlationfunction

!> realspace autocorrelation functions
type lo_correlationfunction
    !> number of unique atoms
    integer :: na
    !> number of timesteps
    integer :: nt
    !> time-axis
    real(flyt), dimension(:), allocatable :: x
    !> total autocorrelation function
    real(flyt), dimension(:), allocatable :: y
    !> correlation function per atom
    real(flyt), dimension(:, :), allocatable :: yp
contains
!        !> Calculate the mean square displacement
!        procedure :: generate
!        !> Write it to file
!!        procedure :: write_to_hdf5
!        procedure :: write_to_plaintext
end type

contains

! !> calculate the autocorrelation functions
! subroutine generate(cf,uc,ss,pm,sim)
!     !> autocorrelation functions
!     class(lo_correlationfunction), intent(out) :: cf
!     !> unitcell
!     type(lo_crystalstructure), intent(in) :: uc
!     !> supercell
!     type(lo_crystalstructure), intent(in) :: ss
!     !> pair mapping
!     type(lo_pairmapping), intent(in) :: pm
!     !> the simulation data
!     type(lo_mdsim), intent(in) :: sim
!
!     real(flyt), dimension(:), allocatable :: dy
!     real(flyt) :: t0
!     integer, dimension(:), allocatable :: dum
!     integer :: a1,una,i,j,k,l,t,dt
!
!     ! some constants and space
!     cf%na=pm%na
!     cf%nt=sim%nt
!     lo_allocate(cf%x(cf%nt))
!     lo_allocate(cf%y(cf%nt))
!     lo_allocate(cf%yp(cf%nt,cf%na))
!     cf%x=0.0_flyt
!     cf%y=0.0_flyt
!     cf%yp=0.0_flyt
!     call lo_linspace(0.0_flyt,(cf%nt-1)*sim%ts,cf%x)
!     lo_allocate(dum(cf%na))
!     dum=0
!     do i=1,pm%nass
!         j=pm%ssatom(i)%unique_atom
!         dum(j)=dum(j)+1
!     enddo
!
!     ! start transforming
!     t0=walltime()
!     call lo_progressbar_init()
!     lo_allocate(dy(cf%nt))
!     do i=1,sim%na
!         ! where to store it
!         a1=pm%ssatom(i)%unique_atom
!         do j=1,3
!             dy=sim%v(j,i,:)
!             call lo_acf(dy)
!             cf%yp(:,a1)=cf%yp(:,a1)+dy/(dum(a1)*3.0_flyt)
!         enddo
!         call lo_progressbar(' ... correlation functions',i,sim%na,walltime()-t0)
!     enddo
!
! end subroutine
!
! !> Dump everything to file
! subroutine write_to_hdf5(cf)
!     !> mean square displacement
!     class(lo_correlationfunction), intent(in) :: cf
!
! !    real(flyt), dimension(:), allocatable :: dumy
! !    integer(HID_T) :: file_id,dset_id
! !    integer :: i,j,k,l
! !    integer :: a1,sh,nt
! !    character(len=1000) :: filename,dname
! !
! !    filename='outfile.mean_square_displacement.hdf5'
! !
! !    call h5open_f(lo_status)
! !    call h5fcreate_f(trim(filename), H5F_ACC_TRUNC_F, file_id, lo_status)
! !
! !    call lo_h5_store_attribute(msd%na,file_id,'number_unique_atoms',lo_status)
! !
! !    ! store the time axis
! !    call lo_h5_store_data(msd%x1,file_id,'time_axis')
! !        call h5dopen_f(file_id, 'time_axis', dset_id, lo_status)
! !        call lo_h5_store_attribute('femtosecond',dset_id,'unit',lo_status)
! !        call h5dclose_f(dset_id,lo_status)
! !    ! store the total MSD
! !    lo_allocate(dumy(msd%nt))
! !    dumy=0.0_flyt
! !    do a1=1,msd%na
! !        dumy=dumy+msd%y1(:,a1)
! !    enddo
! !    call lo_h5_store_data(dumy,file_id,'mean_square_displacement')
! !        call h5dopen_f(file_id,'mean_square_displacement', dset_id, lo_status)
! !        call lo_h5_store_attribute('Angstrom^2',dset_id,'unit',lo_status)
! !        call h5dclose_f(dset_id,lo_status)
! !    ! store the projected MSD per atom
! !    if ( msd%na .gt. 1 ) then
! !        do a1=1,msd%na
! !            dumy=msd%y1(:,1)
! !            dname='mean_square_displacement_atom_'//trim(int2char(a1))
! !            call lo_h5_store_data(dumy,file_id,trim(dname))
! !                call h5dopen_f(file_id,trim(dname), dset_id, lo_status)
! !                call lo_h5_store_attribute('Angstrom^2',dset_id,'unit',lo_status)
! !                call h5dclose_f(dset_id,lo_status)
! !        enddo
! !    endif
! !    !
! !    call h5fclose_f(file_id, lo_status)
! !    call h5close_f(lo_status)
!
! end subroutine
!
! !> Dump it to a plaintext file
! subroutine write_to_plaintext(cf,pm,uc)
!     !> mean square displacement
!     class(lo_correlationfunction), intent(in) :: cf
!     !> pair mapping
!     type(lo_pairmapping), intent(in) :: pm
!     !> unitcell
!     type(lo_crystalstructure), intent(in) :: uc
!
!     integer :: a1,i,j
!     character(len=1000), dimension(:), allocatable :: legend
!     real(flyt), dimension(:,:), allocatable :: dum
!
!     lo_allocate(dum(cf%na+1,cf%nt))
!     dum(1,:)=cf%x
!     do a1=1,cf%na
!         dum(a1+1,:)=cf%yp(:,a1)
!     enddo
!     call lo_dump_gnuplot_2d_real(dum,'outfile.velocity_autocorrelation_function',xlabel='Time (fs)',ylabel='something')
!     lo_deallocate(dum)
!
! !    ! Only one kind of atom
! !    if ( msd%na .eq. 1 ) then
! !        ! Only one kind of atom
! !        lo_allocate(dum(2,msd%nt))
! !        dum(1,:)=msd%x1
! !        dum(2,:)=msd%y1(:,1)
! !        call lo_dump_gnuplot_2d_real(dum,'outfile.mean_square_displacement',xlabel='Simulation time (fs)',ylabel='Mean square displacement (Å^2)')
! !        lo_deallocate(dum)
! !    else
! !        ! More than one kind of atom
! !        lo_allocate(dum(2+msd%na,msd%nt))
! !        lo_allocate(legend(msd%na+1))
! !        legend(1)='Total'
! !        dum=0.0_flyt
! !        dum(1,:)=msd%x1
! !        do a1=1,msd%na
! !            dum(2,:)=dum(2,:)+msd%y1(:,a1)
! !            dum(2+a1,:)=msd%y1(:,a1)
! !            ! get the annoying legend reasonable
! !            j=pm%atom(a1)%unitcell_indices(1)
! !            legend(a1+1)=trim(uc%atomic_symbol(uc%species(j)))//':'
! !            do i=1,pm%atom(a1)%n_uc_atoms
! !                j=pm%atom(a1)%unitcell_indices(i)
! !                legend(a1+1)=trim(legend(a1+1))//' atom '//trim(int2char(j))
! !            enddo
! !        enddo
! !        call lo_dump_gnuplot_2d_real(dum,'outfile.mean_square_displacement',xlabel='Simulation time (fs)',ylabel='Mean square displacement (Å^2)',legend=legend)
! !        lo_deallocate(dum)
! !    endif
!
! end subroutine

end module
