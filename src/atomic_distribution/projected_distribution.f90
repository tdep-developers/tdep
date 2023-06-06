#include "precompilerdefinitions"
module timedistance_correlation
use konstanter, only: flyt
!use gottochblandat
use type_crystalstructure, only: lo_crystalstructure
use type_mdsim, only: lo_mdsim
!use pairmapping
!use pair_distribution
use hdf5_wrappers, only: lo_h5_store_data, HID_T, H5F_ACC_TRUNC_F, h5close_f, h5open_f, h5fclose_f, h5fopen_f, h5fcreate_f
!use fftw_wrappers

implicit none

private
public :: lo_timedistance_correlation

!> pair distribution function for one shell
type lo_tdc_atom_shell
!    !> multiplicity of this specific shell
!    real(flyt) :: multiplicity
!    !> shortest distance
!    real(flyt) :: rmin
!    !> longest distance
!    real(flyt) :: rmax
!    !> bin centers, r-axis
!    real(flyt), dimension(:), allocatable :: x
!    !> time axis
!    real(flyt), dimension(:), allocatable :: y
    !> histogram
    real(flyt), dimension(:, :), allocatable :: h
end type

!> pair distribution function for one symmetry-inequvialent pair
type lo_tdc_atom
    !> number of shells
    integer :: nshell
    !> shell
    type(lo_tdc_atom_shell), dimension(:), allocatable :: shell
end type

!> symmetry decomposed pair distribution function
type lo_timedistance_correlation
    !> number of unique atoms
    integer :: na
    !> number of timesteps
    integer :: nt
    !> number of bins
    integer :: nbin
    !> one distribution per unique atom
    type(lo_tdc_atom), dimension(:), allocatable :: atom
    !> the distance axis
    real(flyt), dimension(:), allocatable :: x
    !> time axis
    real(flyt), dimension(:), allocatable :: y
    !> total histogram
    real(flyt), dimension(:, :), allocatable :: h
contains
!        !> Bin the timesteps
!        procedure :: generate
!        !> Write it to file
!        procedure :: write_to_hdf5
end type

contains

! !> Dump everything to file
! subroutine write_to_hdf5(tdc)
!     class(lo_timedistance_correlation), intent(in) :: tdc
!
!     real(flyt), dimension(:,:), allocatable :: dumx,dumy
!     integer(HID_T) :: file_id,dset_id
!     integer :: i,j,k,l
!     integer :: a1,sh,nb
!     character(len=1000) :: filename,dname
!
!     filename='outfile.time_distance_correlation.hdf5'
!
!     call h5open_f(lo_status)
!     call h5fcreate_f(trim(filename), H5F_ACC_TRUNC_F, file_id, lo_status)
!
!     ! store the total
!     call lo_h5_store_data(tdc%x,file_id,'r_axis')
!     call lo_h5_store_data(tdc%y,file_id,'time_axis')
!     call lo_h5_store_data(tdc%h,file_id,'total_radial_correlation_function')
! !        call h5dopen_f(file_id, 'r_axis', dset_id, lo_status)
! !        call lo_h5_store_attribute('Angstrom',dset_id,'unit',lo_status)
! !        call h5dclose_f(dset_id,lo_status)
! !    call lo_h5_store_data(pdf%y,file_id,'radial_pair_distribution_function')
! !
! !    ! store the projected
! !    call lo_h5_store_attribute(pdf%na,file_id,'number_unique_atoms',lo_status)
! !    do a1=1,pdf%na
! !
! !        ! pack the data
! !        if ( pdf%diffusion ) then
! !            nb=size(pdf%atom(a1)%shell(1)%x,1)
! !            lo_allocate(dumx(pdf%atom(a1)%nshell,nb))
! !            lo_allocate(dumy(pdf%atom(a1)%nshell,nb))
! !            do sh=1,pdf%atom(a1)%nshell
! !                dumx(sh,:)=pdf%atom(a1)%shell(sh)%x
! !                dumy(sh,:)=pdf%atom(a1)%shell(sh)%y
! !            enddo
! !        else
! !            nb=size(pdf%atom(a1)%shell(1)%x,1)
! !            lo_allocate(dumx(pdf%atom(a1)%nshell-1,nb))
! !            lo_allocate(dumy(pdf%atom(a1)%nshell-1,nb))
! !            do sh=2,pdf%atom(a1)%nshell
! !                dumx(sh-1,:)=pdf%atom(a1)%shell(sh)%x
! !                dumy(sh-1,:)=pdf%atom(a1)%shell(sh)%y
! !            enddo
! !        endif
! !
! !        dname='projected_r_axis_'//trim(int2char(a1))
! !        call lo_h5_store_data(dumx,file_id,trim(dname))
! !        dname='projected_pair_distribution_function_atom_'//trim(int2char(a1))
! !        call lo_h5_store_data(dumy,file_id,trim(dname))
! !
! !        lo_deallocate(dumx)
! !        lo_deallocate(dumy)
! !    enddo
! !
!     call h5fclose_f(file_id, lo_status)
!     call h5close_f(lo_status)
!
! end subroutine
!
! !> Bin the pair distribution function
! subroutine generate(tdc,pdf,uc,ss,pm,sim,nbin)
!     !> the histogram stuff
!     class(lo_timedistance_correlation), intent(out) :: tdc
!     !> the total pair correlation function
!     class(lo_pair_distribution), intent(in) :: pdf
!     !> unitcell
!     type(lo_crystalstructure), intent(in) :: uc
!     !> supercell
!     type(lo_crystalstructure), intent(in) :: ss
!     !> pair mapping
!     type(lo_pairmapping), intent(in) :: pm
!     !> the simulation data
!     type(lo_mdsim), intent(in) :: sim
!     !> how many bins in the histogram
!     integer, intent(in) :: nbin
!     !
!     real(flyt), dimension(:,:), allocatable :: trajectory
!     real(flyt), dimension(:,:), allocatable :: dum
!     real(flyt), dimension(:), allocatable :: timebinweight
!     real(flyt), dimension(:), allocatable :: dumtr
!     real(flyt), dimension(3) :: r0,v0,v1
!     real(flyt) :: rmin,rmax,f0,f1,f2,f3,t0
!     real(flyt) :: ridiff,tidiff
!     real(flyt) :: binsize,bdx1,bdx2
!     integer :: i,j,k,l,ii,jj,kk,dt
!     integer :: a1,a2,sh,i1,una,unsh,t
!     integer :: ctr,npair,stride,timer,timertot
!     character(len=1000) :: fn
!
!     integer :: ntimebin
!
! !    complex(flyt), dimension(:), allocatable :: x,y,z
! !    integer :: n
! !    n=40
! !    allocate(x(n))
! !    allocate(y(n))
! !    allocate(z(n))
! !    do i=1,n
! !        x(i)=sin(i*1.0_flyt)
! !    enddo
! !    y=x
! !    call lo_fft(y)
! !    z=y
! !    call lo_ifft(z)
! !    !
! !    z=x
! !    call lo_acf(z)
! !
! !    do i=1,n
! !        write(*,*) i,x(i),z(i)
! !    enddo
!
! !stop
!
! write(*,*) '... hello hello, binning and stuff'
!
!     ! Make some space. When I did the normal pair distribution I already
!     ! figured out nice limits for binning.
!     tdc%nt=sim%nt  !500 !floor(sim%nt*0.5_flyt)
!     !tdc%nt=sim%nt !floor(sim%nt*0.5_flyt)
!     tdc%na=pdf%na
!     tdc%nbin=pdf%nbin
!     lo_allocate(tdc%x(tdc%nbin))
!     lo_allocate(tdc%y(tdc%nt))
!     lo_allocate(tdc%atom(tdc%na))
!     tdc%x=pdf%x
!     call lo_linspace(0.0_flyt,(tdc%nt)*sim%ts,tdc%y)
!     do a1=1,tdc%na
!         tdc%atom(a1)%nshell=pdf%atom(a1)%nshell
!         lo_allocate(tdc%atom(a1)%shell( tdc%atom(a1)%nshell ))
!         do sh=1,tdc%atom(a1)%nshell
!             lo_allocate(tdc%atom(a1)%shell(sh)%h( tdc%nbin, tdc%nt ))
!             tdc%atom(a1)%shell(sh)%h=0.0_flyt
!         enddo
!     enddo
!
!     npair=0
!     do a1=1,pm%nass
!         npair=npair+pm%ssatom(a1)%npair
!     enddo
!
!     lo_allocate(timebinweight(tdc%nt))
!     timebinweight=0.0_flyt
!     tidiff=(1.0_flyt*tdc%nt)/(sim%nt*1.0_flyt)
!     do i=1,sim%nt
!         jj=floor( (i-1.0_flyt)*tidiff )+1
!         timebinweight(jj)=timebinweight(jj)+1.0_flyt
!     enddo
!     do i=1,tdc%nt
!         timebinweight(i)=tdc%nt*1.0_flyt/timebinweight(i)
!     enddo
!
! write(*,*) '... made space for correlation function'
!
!     ! Start with the brutal binning thing. This can take a while.
!     t0=walltime()
!     timer=0
!     call lo_progressbar_init()
!     lo_allocate(trajectory(3,sim%nt))
!     lo_allocate(dumtr(sim%nt))
!     do a1=1,pm%nass
!     do i1=1,pm%ssatom(a1)%npair
!         ! where should it be binned?
!         una=pm%ssatom(a1)%unique_atom
!         unsh=pm%ssatom(a1)%pair(i1)%shell
!         ! skip selfterm
!         if ( norm2(pm%atom(una)%shell(unsh)%prototype_vector) .lt. lo_tol ) cycle
!         ! get the relevant trajectory and aux stuff
!         call extract_pair_trajectory(pm,a1,i1,sim,ss,trajectory)
!         ! parameters binning
!         lo_allocate(dum(tdc%nbin,tdc%nt))
!         rmin=0.0_flyt
!         rmax=maxval(tdc%x)
!         ridiff=tdc%nbin/(rmax-rmin)
!         tidiff=(1.0_flyt*tdc%nt)/(sim%nt*1.0_flyt)
!         binsize=tdc%x(2)-tdc%x(1)
!         ! the equilibrium distance thing
!         a2=pm%ssatom(a1)%pair(i1)%i2
!         r0=pm%ssatom(a1)%pair(i1)%lv
!         r0=r0+ss%rcart(:,a2)-ss%rcart(:,a1)
!         ! actual binning
!         dum=0.0_flyt
!         dumtr=0.0_flyt
!         do dt=1,sim%nt ! tdc%nt !sim%nt
!             !
!             v0=r0+sim%u(:,a2,dt)-sim%u(:,a1,dt)
!             dumtr(dt)=norm2(v0)
!             !
! !            ! where to bin in time
! !            jj=dt !floor( (dt-1.0_flyt)*tidiff )+1
! !            ! where to bin in space
! !            v0=r0+sim%u(:,a2,dt)-sim%u(:,a1,dt)
! !            f0=norm2(v0)
! !            ii=floor( (f0-rmin)*ridiff )+1
! !            if ( ii .gt. 1 .and. ii .le. tdc%nbin ) then
! !                dum(ii,jj)=dum(ii,jj)+1.0_flyt*timebinweight(jj)
! !            endif
!
! !            v0=r0+sim%u(:,a2,dt)-sim%u(:,a1,dt)
! !            v1=r0+sim%u(:,a2,dt+1)-sim%u(:,a1,dt+1)
! !            f0=norm2(v0)
! !            f1=norm2(v0)
! !            !
! !            ii=floor( (f0-rmin)*ridiff )+1
! !            jj=floor( (f1-rmin)*ridiff )+1
! !            if ( ii .eq. jj ) then
! !                ! add it in the same box
! !                dum(ii,dt)=dum(ii,dt)+1.0_flyt
! !            elseif ( ii .lt. jj ) then
! !                !
! !                bdx1=binsize-mod(f0,binsize)
! !                bdx2=mod(f1,binsize)
! !                f2=max(jj-ii-1,0)*binsize+bdx1+bdx2
! !                f2=1.0_flyt/f2
! !                !
! !                dum(ii,dt)=dum(ii,dt)+bdx1*f2
! !                dum(jj,dt)=dum(jj,dt)+bdx2*f2
! !                do i=ii+1,jj-1
! !                    dum(i,dt)=dum(i,dt)+binsize*f2
! !                enddo
! !            else
! !                bdx1=mod(f0,binsize)
! !                bdx2=binsize-mod(f1,binsize)
! !                f2=max(ii-jj-1,0)*binsize+bdx1+bdx2
! !                f2=1.0_flyt/f2
! !                !
! !                dum(ii,dt)=dum(ii,dt)+bdx1*f2
! !                dum(jj,dt)=dum(jj,dt)+bdx2*f2
! !                do i=jj+1,ii-1
! !                    dum(i,dt)=dum(i,dt)+binsize*f2
! !                enddo
! !            endif
!
!             !
!             !ctr=0
!             !do t=1,tdc%nt,1 !sim%nt-dt,5 ! add stride here perhaps
!             !!do t=1,tdc%nt-dt,5 !sim%nt-dt,5 ! add stride here perhaps
!             !    !v0=trajectory(:,t)-trajectory(:,t+dt-1)
!             !    v0=r0+sim%u(:,a2,t)-sim%u(:,a1,t+dt-1)
!             !    f0=norm2(v0)
!             !    ii=floor( (f0-rmin)*ridiff )+1
!             !    if ( ii .gt. 1 .and. ii .le. tdc%nbin ) then
!             !        ctr=ctr+1
!             !        dum(ii,dt)=dum(ii,dt)+1.0_flyt
!             !    endif
!             !enddo
!             !if ( ctr .gt. 0 ) then
!             !    dum(:,dt)=dum(:,dt)/ctr
!             !else
!             !    dum(:,dt)=0.0_flyt
!             !endif
!         enddo
!         ! Fourier transform in-place
!         call lo_fft(dumtr)
!
!         ! add it to the total
!         tdc%atom(una)%shell(unsh)%h=tdc%atom(una)%shell(unsh)%h+dum
!         lo_deallocate(dum)
!         !
!         timer=timer+1
!         call lo_progressbar(' ... timespace correlation',timer,npair,walltime()-t0)
!         !
!     enddo
!     enddo
!     call lo_progressbar(' ... timespace correlation',npair,npair,walltime()-t0)
!
!     ! no idea how this should be normalized properly.
!     do a1=1,pm%na
!     do sh=1,pm%atom(a1)%nshell
!         if ( norm2(pm%atom(a1)%shell(sh)%prototype_vector) .lt. lo_tol ) cycle
!         ! normalize per time
!         do i=1,tdc%nt
!             tdc%atom(a1)%shell(sh)%h(:,i)=tdc%atom(a1)%shell(sh)%h(:,i)/sum(tdc%atom(a1)%shell(sh)%h(:,i))
!         enddo
!         ! and by the distance thing
!         f0=f0*maxval(tdc%x)/tdc%nbin
!         f0=f0*npair/ss%volume
!         f0=f0*4*lo_pi
!         do i=1,tdc%nbin
!             f1=f0*(tdc%x(i)**2)
!             if ( f1 .gt. lo_sqtol ) then
!                 tdc%atom(a1)%shell(sh)%h(i,:)=tdc%atom(a1)%shell(sh)%h(i,:)/f1
!             else
!                 tdc%atom(a1)%shell(sh)%h(i,:)=0.0_flyt
!             endif
!         enddo
!     enddo
!     enddo
!
!     timertot=0
!     do a1=1,tdc%na
!         timertot=timertot+tdc%atom(a1)%nshell
!     enddo
!
!     t0=walltime()
!     call lo_progressbar_init()
!     lo_allocate(tdc%h(pdf%nbin,tdc%nt))
!     lo_allocate(dum(pdf%nbin,tdc%nt))
!     tdc%x=pdf%x
!     call lo_linspace(0.0_flyt,(tdc%nt-1)*sim%ts,tdc%y)
!     tdc%h=0.0_flyt
!     timer=0
!     do a1=1,tdc%na
!     do sh=1,tdc%atom(a1)%nshell
!         ! skip selfterm
!         if ( norm2(pm%atom(a1)%shell(sh)%prototype_vector) .lt. lo_tol ) cycle
!         tdc%h=tdc%h+tdc%atom(a1)%shell(sh)%h
!         !
!         timer=timer+1
!         call lo_progressbar(' ... rearranging',timer,timertot,walltime()-t0)
!     enddo
!     enddo
!     call lo_progressbar(' ... rearranging',timertot,timertot,walltime()-t0)
!
!     contains
!
! !    !> take a histogram from one set of x-y axes to another
! !    subroutine put_shape_on_new_mesh(x,y,xi,yi)
! !        real(flyt), dimension(:), intent(in) :: x,y
! !        real(flyt), dimension(:), intent(in) :: xi
! !        real(flyt), dimension(:), intent(out) :: yi
! !        !
! !        integer :: i,j,n,ni
! !        integer :: ii,jj
! !        real(flyt) :: f0,sigma
! !        real(flyt) :: ei,eimin,eimax,eirange,foursigma
! !        !
! !        n=size(x,1)
! !        ni=size(xi,1)
! !        yi=0.0_flyt
! !        sigma=(maxval(xi)-minval(xi))*3.0_flyt/ni
! !        !
! !        eimin=minval(x)
! !        eimax=maxval(x)
! !        eirange=1.0_flyt/(eimax-eimin)
! !        foursigma=4*sigma
! !        do i=1,ni
! !            ! find lower and upper bound
! !            ei=xi(i)
! !            ii=floor(n*(ei-foursigma-eimin)*eirange)
! !            ii=max(ii,1)
! !            ii=min(ii,n)
! !            !
! !            jj=ceiling(n*(ei+foursigma-eimin)*eirange)
! !            jj=min(jj,n)
! !            jj=max(jj,1)
! !            do j=ii,jj
! !                f0=lo_gauss(xi(i),x(j),sigma)
! !                yi(i)=yi(i)+f0*y(j)
! !            enddo
! !        enddo
! !    end subroutine
!
!     !> get the trajectory for one pair
!     subroutine extract_pair_trajectory(pm,a1,i1,sim,ss,trajectory)
!         !> pair mapping
!         type(lo_pairmapping), intent(in) :: pm
!         !> supercell atom
!         integer, intent(in) :: a1
!         !> supercell pair
!         integer, intent(in) :: i1
!         !> the simulation data
!         type(lo_mdsim), intent(in) :: sim
!         !> supercell
!         type(lo_crystalstructure), intent(in) :: ss
!         !> the trajectory
!         real(flyt), dimension(:,:), intent(out) :: trajectory
!
!         real(flyt), dimension(3) :: r0,r1,v1,v2
!         integer :: ii,jj,t,i,j
!
!         ! lattice vector added to all the pair vectors
!         r0=pm%ssatom(a1)%pair(i1)%lv
!         ! indices to the supercell
!         ii=pm%ssatom(a1)%pair(i1)%i1
!         jj=pm%ssatom(a1)%pair(i1)%i2
!         ! the equilibrium pair distance
!         r0=r0+ss%rcart(:,jj)-ss%rcart(:,ii)
!
!         ! this pair is supposed to be lv + ss(jj) - ss(ii)
!         do t=1,sim%nt
!             trajectory(:,t)=r0+sim%u(:,jj,t)-sim%u(:,ii,t)
!         enddo
!     end subroutine
!
! end subroutine
!

end module
