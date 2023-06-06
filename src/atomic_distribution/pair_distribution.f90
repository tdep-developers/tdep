#include "precompilerdefinitions"
module pair_distribution
use konstanter, only: r8, lo_huge, lo_hugeint, lo_pi, lo_sqtol, lo_tol, lo_status, lo_bohr_to_A, lo_time_au_to_fs
use gottochblandat, only: walltime, lo_progressbar_init, lo_progressbar, lo_trapezoid_integration, lo_gauss, tochar, lo_linspace, &
                          lo_return_unique
use type_crystalstructure, only: lo_crystalstructure
use type_mdsim, only: lo_mdsim
use pairmapping, only: lo_pairmapping
use lo_memtracker, only: lo_mem_helper
use mpi_wrappers, only: lo_mpi_helper, lo_stop_gracefully
use hdf5_wrappers, only: lo_hdf5_helper

implicit none
private
public :: lo_pair_distribution

!> pair distribution per shell
type lo_pair_distribution_shell
    !> smallest distance
    real(r8) :: rmin = -lo_huge
    !> largest distance
    real(r8) :: rmax = -lo_huge
    !> shell-resolved histogram
    real(r8), dimension(:), allocatable :: x
    real(r8), dimension(:), allocatable :: y
    real(r8), dimension(:, :), allocatable :: z
end type

!> symmetry decomposed pair distribution function
type lo_pair_distribution
    !> number of bins in the histogram
    integer :: nbin = -lo_hugeint

    !> shell-resolved pdf
    type(lo_pair_distribution_shell), dimension(:), allocatable :: sh
    !> smearing parameter for binning
    real(r8) :: sigma = -lo_huge
contains
    !> Bin the timesteps
    procedure :: bin
    !> Write it to file
    procedure :: write_to_hdf5
end type

contains

!> Bin the pair distribution function
subroutine bin(pdf, uc, ss, pm, sim, nbin, nodiffusion, mw, mem, verbosity)
    !> the histogram stuff
    class(lo_pair_distribution), intent(out) :: pdf
    !> unitcell
    type(lo_crystalstructure), intent(in) :: uc
    !> supercell
    type(lo_crystalstructure), intent(in) :: ss
    !> pair mapping
    type(lo_pairmapping), intent(in) :: pm
    !> the simulation data
    type(lo_mdsim), intent(in) :: sim
    !> how many bins in the histogram
    integer, intent(in) :: nbin
    !> always assume no diffusion?
    logical, intent(in) :: nodiffusion
    !> mpi helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem
    !> talk a lot?
    integer, intent(in) :: verbosity

    real(r8), parameter :: buf = 0.1_r8 ! buffer the histograms a little bit
    real(r8) :: timer, t0, t1

    ! start timers
    timer = walltime()
    t0 = timer
    t1 = timer

    ! Set some basics
    init: block
        integer :: ish

        if (verbosity .gt. 0) then
            write (*, *) ''
            write (*, *) 'Calculating pair distribution function'
        end if

        ! How many bins?
        pdf%nbin = nbin

        ! Space for shells
        allocate (pdf%sh(pm%n_shell))
        do ish = 1, pm%n_shell
            pdf%sh(ish)%rmin = lo_huge
            pdf%sh(ish)%rmax = -lo_huge
            allocate (pdf%sh(ish)%x(pdf%nbin))
            allocate (pdf%sh(ish)%y(pdf%nbin))
            allocate (pdf%sh(ish)%z(pdf%nbin, sim%nt))
            pdf%sh(ish)%x = 0.0_r8
            pdf%sh(ish)%y = 0.0_r8
            pdf%sh(ish)%z = 0.0_r8
        end do
    end block init

    ! Do a first sweep
    firstpass: block
        real(r8), dimension(:, :), allocatable :: dr0
        real(r8) :: f0, f1
        integer :: ish, ctr, ipair, t, ctrtot

        ctrtot = 0
        do ish = 1, pm%n_shell
            ctrtot = ctrtot + pm%sh(ish)%n_ss_pair
        end do

        call mem%allocate(dr0, [3, sim%nt], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        dr0 = 0.0_r8

        if (verbosity .gt. 0) call lo_progressbar_init()
        ctr = 0
        do ish = 1, pm%n_shell
        do ipair = 1, pm%sh(ish)%n_ss_pair
            ctr = ctr + 1
            if (mod(ctr, mw%n) .ne. mw%r) cycle
            ! fetch the trajectory?
            call extract_trajectory(pm, ish, ipair, sim, dr0)
            do t = 1, sim%nt
                f0 = norm2(dr0(:, t))
                pdf%sh(ish)%rmin = min(pdf%sh(ish)%rmin, f0)
                pdf%sh(ish)%rmax = max(pdf%sh(ish)%rmax, f0)
            end do

            if (verbosity .gt. 0 .and. ctr .lt. ctrtot) then
                call lo_progressbar(' ... measuring span', ctr, ctrtot, walltime() - t0)
            end if
        end do
        end do

        ! sync and build axis
        f0 = 0.0_r8
        ctr = 0
        do ish = 1, pm%n_shell
            call mw%allreduce('min', pdf%sh(ish)%rmin)
            call mw%allreduce('max', pdf%sh(ish)%rmax)
            if (norm2(pm%sh(ish)%r) .gt. lo_tol) then
                ctr = ctr + 1
                f1 = (pdf%sh(ish)%rmax - pdf%sh(ish)%rmin)/real(pdf%nbin, r8)
                f0 = max(f0, f1)
            end if
        end do
        pdf%sigma = f0*2

        do ish = 1, pm%n_shell
            pdf%sh(ish)%rmin = pdf%sh(ish)%rmin - 7*pdf%sigma
            pdf%sh(ish)%rmax = pdf%sh(ish)%rmax + 7*pdf%sigma
            call lo_linspace(pdf%sh(ish)%rmin, pdf%sh(ish)%rmax, pdf%sh(ish)%x)
        end do

        if (verbosity .gt. 0) then
            t1 = walltime()
            call lo_progressbar(' ... measuring span', ctrtot, ctrtot, t1 - t0)
            t0 = t1
        end if

        call mem%deallocate(dr0, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    end block firstpass

    ! Do a second sweep
    secondpass: block
        real(r8), dimension(:, :), allocatable :: dr0
        real(r8) :: f0, f1, f2, f3, invf
        integer :: ish, ctr, ipair, t, ctrtot
        integer :: i, j, k, l, ilo, ihi

        ctrtot = 0
        do ish = 1, pm%n_shell
            ctrtot = ctrtot + pm%sh(ish)%n_ss_pair
        end do

        call mem%allocate(dr0, [3, sim%nt], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        dr0 = 0.0_r8

        if (verbosity .gt. 0) call lo_progressbar_init()
        ctr = 0
        do ish = 1, pm%n_shell
            invf = pdf%nbin/(pdf%sh(ish)%rmax - pdf%sh(ish)%rmin)
            do ipair = 1, pm%sh(ish)%n_ss_pair
                ctr = ctr + 1
                if (mod(ctr, mw%n) .ne. mw%r) cycle
                ! fetch the trajectory?
                call extract_trajectory(pm, ish, ipair, sim, dr0)
                do t = 1, sim%nt
                    f0 = norm2(dr0(:, t))  ! pair distance
                    f1 = (f0 - pdf%sh(ish)%rmin - 4*pdf%sigma)*invf
                    f2 = (f0 - pdf%sh(ish)%rmin + 4*pdf%sigma)*invf
                    ilo = floor(f1) + 1
                    ihi = ceiling(f2) + 1
                    do i = ilo, ihi
                        f3 = lo_gauss(pdf%sh(ish)%x(i), f0, pdf%sigma)
                        pdf%sh(ish)%z(i, t) = pdf%sh(ish)%z(i, t) + f3
                    end do
                end do

                if (verbosity .gt. 0 .and. ctr .lt. ctrtot - 1) then
                    call lo_progressbar(' ... binning', ctr, ctrtot, walltime() - t0)
                end if
            end do
        end do

        ! sync and build axis
        do ish = 1, pm%n_shell
            call mw%allreduce('sum', pdf%sh(ish)%z)
        end do

        if (verbosity .gt. 0) then
            t1 = walltime()
            call lo_progressbar(' ... binning', ctrtot, ctrtot, t1 - t0)
            t0 = t1
        end if

        call mem%deallocate(dr0, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    end block secondpass

    ! Normalize this sensibly.
    finalize: block
        real(r8) :: f0
        integer :: ish, t

        do ish = 1, pm%n_shell
            pdf%sh(ish)%y = 0.0_r8
            do t = 1, sim%nt
                pdf%sh(ish)%y = pdf%sh(ish)%y + pdf%sh(ish)%z(:, t)
            end do
            f0 = lo_trapezoid_integration(pdf%sh(ish)%x, pdf%sh(ish)%y)
            pdf%sh(ish)%y = pdf%sh(ish)%y*pm%sh(ish)%n_ss_pair/f0
        end do
    end block finalize

    ! integer, dimension(:,:), allocatable :: uniquepairs
    ! real(r8) :: timer,cutoff
    !
    ! init: block
    !     integer, dimension(:,:), allocatable :: di1
    !     integer :: i,ii,jj,ctr,a1,sh
    !
    !     timer=walltime()
    !     write(*,*) ''
    !     write(*,*) 'Building pair distribution function'
    !
    !     ! First fill out my structure thing
    !     ! make some space
    !     pdf%na=pm%na
    !     pdf%nbin=nbin
    !     lo_allocate(pdf%atom(pdf%na))
    !     do a1=1,pdf%na
    !         pdf%atom(a1)%nshell=pm%atom(a1)%nshell
    !         lo_allocate(pdf%atom(a1)%shell( pdf%atom(a1)%nshell ))
    !         do sh=1,pdf%atom(a1)%nshell
    !             pdf%atom(a1)%shell(sh)%rmin=lo_huge
    !             pdf%atom(a1)%shell(sh)%rmax=-lo_huge
    !             pdf%atom(a1)%shell(sh)%multiplicity=0.0_r8
    !
    !             ii=pm%atom(a1)%shell(sh)%i1
    !             jj=pm%atom(a1)%shell(sh)%i2
    !             pdf%atom(a1)%label=trim(uc%atomic_symbol(uc%species(ii)))
    !             pdf%atom(a1)%shell(sh)%label=trim(uc%atomic_symbol(uc%species(jj)))
    !             pdf%atom(a1)%atnum=uc%atomic_number(ii)
    !             pdf%atom(a1)%shell(sh)%atnum=uc%atomic_number(jj)
    !         enddo
    !     enddo
    !
    !     ! Figure out which kind of pairs there are
    !     ctr=0
    !     do a1=1,pdf%na
    !         ctr=ctr+pdf%atom(a1)%nshell
    !     enddo
    !     lo_allocate(di1(2,ctr))
    !     di1=0
    !     ctr=0
    !     do a1=1,pdf%na
    !     do sh=1,pdf%atom(a1)%nshell
    !         ctr=ctr+1
    !         if ( pdf%atom(a1)%atnum < pdf%atom(a1)%shell(sh)%atnum ) then
    !             ii=pdf%atom(a1)%atnum
    !             jj=pdf%atom(a1)%shell(sh)%atnum
    !         else
    !             ii=pdf%atom(a1)%shell(sh)%atnum
    !             jj=pdf%atom(a1)%atnum
    !         endif
    !         di1(:,ctr)=[ii,jj]
    !     enddo
    !     enddo
    !     call lo_return_unique(di1,uniquepairs)
    !     pdf%npair=size(uniquepairs,2)
    !     ! And finally make sure I know what kind of pair each pair is
    !     do a1=1,pdf%na
    !     do sh=1,pdf%atom(a1)%nshell
    !         ctr=ctr+1
    !         if ( pdf%atom(a1)%atnum < pdf%atom(a1)%shell(sh)%atnum ) then
    !             ii=pdf%atom(a1)%atnum
    !             jj=pdf%atom(a1)%shell(sh)%atnum
    !         else
    !             ii=pdf%atom(a1)%shell(sh)%atnum
    !             jj=pdf%atom(a1)%atnum
    !         endif
    !         do i=1,pdf%npair
    !             if ( sum(abs([ii,jj]-uniquepairs(:,i))) .eq. 0 ) pdf%atom(a1)%shell(sh)%pairind=i
    !         enddo
    !     enddo
    !     enddo
    !
    !     ! And some space for the species-resolved pairs
    !     lo_allocate(pdf%pair(pdf%npair))
    ! end block init
    !
    ! ! Do a first pass to see if things are diffusing or not
    ! firstpass: block
    !     real(r8), dimension(:,:), allocatable :: trajectory
    !     real(r8) :: t0,npair,f0,rmin,rmax
    !     integer :: a1,sh,i1,una,unsh,t,step
    !
    !     ! we will loop over everything twice, just to be on the safe side. The
    !     ! first pass is used to determine the limits of the histogram.
    !     t0=walltime()
    !     call lo_progressbar_init()
    !     cutoff=0.0_r8
    !     npair=0.0_r8
    !     lo_allocate(trajectory(3,sim%nt))
    !     do a1=1,pm%nass
    !         npair=npair+pm%ssatom(a1)%npair
    !         do i1=1,pm%ssatom(a1)%npair
    !             ! fetch the trajectory
    !             call extract_pair_trajectory(pm,a1,i1,sim,ss,trajectory)
    !             ! where should it be binned?
    !             una=pm%ssatom(a1)%unique_atom
    !             unsh=pm%ssatom(a1)%pair(i1)%shell
    !             step=floor(sim%nt/100.0_r8)+1
    !             rmin=lo_huge
    !             rmax=-lo_huge
    !             do t=1,sim%nt,step
    !                 f0=norm2(trajectory(:,t))
    !                 rmin=min(rmin,f0)
    !                 rmax=max(rmax,f0)
    !             enddo
    !             pdf%atom(una)%shell(unsh)%rmin=min(pdf%atom(una)%shell(unsh)%rmin,rmin)
    !             pdf%atom(una)%shell(unsh)%rmax=max(pdf%atom(una)%shell(unsh)%rmax,rmax)
    !             cutoff=max(cutoff,rmax)
    !             ! count how common this type of pair is
    !             pdf%atom(una)%shell(unsh)%multiplicity=pdf%atom(una)%shell(unsh)%multiplicity+1.0_r8
    !         enddo
    !         call lo_progressbar(' ... finding limits',a1,pm%nass,walltime()-t0)
    !     enddo
    !
    !     ! try to figure out wether the system is diffusing or not. If the 'cutoff' number is
    !     ! significantly larger than the input cutoff, there has been some diffusion, and we
    !     ! have to deal with it in a slightly different fashion.
    !     if ( nodiffusion ) then
    !         write(*,*) '... you have decided this is not diffusing'
    !         pdf%diffusion=.false.
    !     else
    !         if ( cutoff .gt. pm%nncutoff+pm%nndist ) then
    !             write(*,*) '... this might be diffusing'
    !             pdf%diffusion=.true.
    !         else
    !             write(*,*) '... probably not diffusing'
    !             pdf%diffusion=.false.
    !         endif
    !     endif
    !
    !     if ( pdf%diffusion ) then
    !         cutoff=pm%nncutoff+0.05*pm%nndist
    !         pdf%nbin=ceiling(nbin*cutoff/pm%nndist)
    !     else
    !         cutoff=cutoff+0.5_r8*pm%nndist
    !         pdf%nbin=ceiling(nbin*cutoff/pm%nndist)
    !         ! Make some space to bin things projected
    !         do a1=1,pdf%na
    !             do sh=1,pdf%atom(a1)%nshell
    !                 lo_allocate(pdf%atom(a1)%shell(sh)%x( nbin ))
    !                 lo_allocate(pdf%atom(a1)%shell(sh)%y( nbin ))
    !                 pdf%atom(a1)%shell(sh)%x=0.0_r8
    !                 pdf%atom(a1)%shell(sh)%y=0.0_r8
    !                 f0=abs(pdf%atom(a1)%shell(sh)%rmax-pdf%atom(a1)%shell(sh)%rmin)
    !                 pdf%atom(a1)%shell(sh)%rmin=pdf%atom(a1)%shell(sh)%rmin-f0*0.25_r8
    !                 pdf%atom(a1)%shell(sh)%rmax=pdf%atom(a1)%shell(sh)%rmax+f0*0.25_r8
    !                 call lo_linspace(pdf%atom(a1)%shell(sh)%rmin,pdf%atom(a1)%shell(sh)%rmax,pdf%atom(a1)%shell(sh)%x)
    !             enddo
    !         enddo
    !     endif
    !
    !     ! Some space for the species-pair projected
    !     do i1=1,pdf%npair
    !         lo_allocate(pdf%pair(i1)%y(pdf%nbin))
    !         pdf%pair(i1)%y=0.0_r8
    !         pdf%pair(i1)%paircounter=0
    !     enddo
    ! end block firstpass
    !
    ! if ( pdf%diffusion ) then
    ! diffusionbin: block
    !     ! Non-diffusing case
    !     real(r8), dimension(:), allocatable :: dum
    !     real(r8), dimension(3) :: v0
    !     real(r8) :: t0,f0,f1,ridiff
    !     integer :: a1,a2,i,ii,jj,pairind,t,ctr
    !
    !     lo_allocate(pdf%x(pdf%nbin))
    !
    !     lo_allocate(dum(pdf%nbin))
    !     call lo_linspace(0.0_r8,cutoff,pdf%x)
    !     ridiff=pdf%nbin/(cutoff)
    !
    !     t0=walltime()
    !     call lo_progressbar_init()
    !     do a1=1,pm%nass
    !         do a2=1,pm%nass
    !             ! What kind of pair is this.
    !             pairind=-1
    !             do i=1,pdf%npair
    !                 if ( ss%atomic_number(a1) < ss%atomic_number(a2) ) then
    !                     ii=ss%atomic_number(a1)
    !                     jj=ss%atomic_number(a2)
    !                 else
    !                     jj=ss%atomic_number(a1)
    !                     ii=ss%atomic_number(a2)
    !                 endif
    !                 if ( sum(abs([ii,jj]-uniquepairs(:,i))) .eq. 0 ) pairind=i
    !             enddo
    !
    !             ! get a partial histogram
    !             dum=0.0_r8
    !             ctr=0
    !             do t=1,sim%nt
    !                 v0=ss%displacement_fractional_to_cartesian( sim%r(:,a2,t)-sim%r(:,a1,t) )
    !                 f0=norm2(v0)
    !                 ii=floor( (f0)*ridiff )+1
    !                 if ( ii .gt. 3 .and. ii .le. pdf%nbin ) then
    !                     ctr=ctr+1
    !                     dum(ii)=dum(ii)+1.0_r8
    !                 endif
    !             enddo
    !             ! add it to the total histogram
    !             if ( ctr .gt. 0 ) then
    !                 pdf%pair(pairind)%y=pdf%pair(pairind)%y+dum/real(ctr,r8)
    !                 pdf%pair(pairind)%paircounter=pdf%pair(pairind)%paircounter+1
    !             endif
    !         enddo
    !         call lo_progressbar(' ... binning',a1,pm%nass,walltime()-t0)
    !     enddo
    !     lo_deallocate(dum)
    !
    !     ! Figure out the normalization
    !     lo_allocate(pdf%y(pdf%nbin))
    !     pdf%y=0.0_r8
    !     jj=0
    !     do ii=1,pdf%npair
    !         jj=jj+pdf%pair(ii)%paircounter
    !         pdf%y=pdf%y+pdf%pair(ii)%y
    !     enddo
    !     f0=lo_trapezoid_integration( pdf%x,pdf%y )
    !     f1=jj/f0/ss%na
    !
    !     do ii=1,pdf%npair
    !         ! Normalize to correct number of pairs
    !         pdf%pair(ii)%y=pdf%pair(ii)%y*f1
    !         ! Then the radial shell thing
    !         do i=1,pdf%nbin
    !             ! volume of radial shell
    !             f0=4*lo_pi*(pdf%x(i)**2)
    !             f0=f0*ss%na/ss%volume
    !             if ( f0 .gt. lo_sqtol ) then
    !                 pdf%pair(ii)%y(i)=pdf%pair(ii)%y(i)/(f0)
    !             else
    !                 pdf%pair(ii)%y(i)=0.0_r8
    !             endif
    !         enddo
    !     enddo
    !     ! Add things up to the total
    !     pdf%y=0.0_r8
    !     do ii=1,pdf%npair
    !         pdf%y=pdf%y+pdf%pair(ii)%y
    !     enddo
    !
    ! end block diffusionbin
    ! else
    ! normalbin: block
    !     ! Diffusion-free case, a little messy
    !     real(r8), dimension(:,:), allocatable :: trajectory
    !     real(r8), dimension(:), allocatable :: dum
    !     real(r8) :: t0,f0,f1,rmax,rmin,foursigma,ridiff,sigma
    !     integer :: a1,i1,una,unsh,nb,ctr,i,ii,jj,kk,t,sh
    !
    !     lo_allocate(trajectory(3,sim%nt))
    !     trajectory=0.0_r8
    !     t0=walltime()
    !     call lo_progressbar_init()
    !     do a1=1,pm%nass
    !     do i1=1,pm%ssatom(a1)%npair
    !         ! fetch the trajectory
    !         call extract_pair_trajectory(pm,a1,i1,sim,ss,trajectory)
    !         ! where should it be binned?
    !         una=pm%ssatom(a1)%unique_atom
    !         unsh=pm%ssatom(a1)%pair(i1)%shell
    !         ! parameters binning
    !         nb=size(pdf%atom(una)%shell(unsh)%x,1)
    !         lo_allocate(dum(nb))
    !         dum=0.0_r8
    !         sigma=(pdf%atom(una)%shell(unsh)%x(2)-pdf%atom(una)%shell(unsh)%x(1))*2
    !         foursigma=4*sigma
    !
    !         rmin=pdf%atom(una)%shell(unsh)%rmin
    !         rmax=pdf%atom(una)%shell(unsh)%rmax
    !         if ( abs(rmax-rmin) .lt. lo_sqtol ) then
    !             ! makes no sense?
    !             lo_deallocate(dum)
    !             cycle
    !         else
    !             ridiff=nb/(rmax-rmin)
    !         endif
    !         ctr=0
    !         do t=1,sim%nt
    !             ! the distance
    !             f0=norm2(trajectory(:,t))
    !             ! gaussian binning
    !             ii=floor( (f0-rmin)*ridiff )+1
    !             if ( ii .gt. 0 .and. ii .le. nb ) then
    !                 ctr=ctr+1
    !                 jj=floor( (f0-foursigma-rmin)*ridiff )+1
    !                 kk=floor( (f0+foursigma-rmin)*ridiff )+1
    !                 jj=max(1,jj)
    !                 kk=min(nb,kk)
    !                 do i=jj,kk
    !                     dum(i)=dum(i)+lo_gauss(0.0_r8,f0-pdf%atom(una)%shell(unsh)%x(i),sigma)
    !                 enddo
    !             endif
    !         enddo
    !
    !         if ( ctr .gt. 0 ) then
    !             dum=ctr*dum/lo_trapezoid_integration(pdf%atom(una)%shell(unsh)%x,dum)
    !             pdf%atom(una)%shell(unsh)%y=pdf%atom(una)%shell(unsh)%y+dum/(ctr*1.0_r8)
    !         endif
    !         lo_deallocate(dum)
    !     enddo
    !     call lo_progressbar(' ... binning',a1,pm%nass,walltime()-t0)
    !     enddo
    !
    !     ! Now that it is binned, normalize it properly
    !     do a1=1,pdf%na
    !         ! how many neighbours should this atom have, excluding itself?
    !         ii=pm%atom(a1)%npair-1
    !         ! count all the probabilities, or what you want to call it
    !         f0=0.0_r8
    !         do sh=1,pdf%atom(a1)%nshell
    !             f0=f0+lo_trapezoid_integration(pdf%atom(a1)%shell(sh)%x,pdf%atom(a1)%shell(sh)%y)
    !         enddo
    !         ! magic renormalizing factor, so that the number of pairs become correct
    !         f1=ii/f0
    !         do sh=1,pdf%atom(a1)%nshell
    !             pdf%atom(a1)%shell(sh)%y=pdf%atom(a1)%shell(sh)%y*f1
    !         enddo
    !         ! at this stage, int n(r) dr = N, the correct number of atoms.
    !         ! We just need to adjust with the density, and the volume of a radial shell.
    !         f1=ss%na/ss%volume ! particle density
    !         do sh=1,pm%atom(a1)%nshell
    !             nb=size(pdf%atom(a1)%shell(sh)%x,1)
    !             do i=1,nb
    !                 ! volume of radial shell
    !                 f0=4*lo_pi*(pdf%atom(a1)%shell(sh)%x(i)**2)
    !                 if ( f0 .gt. lo_sqtol ) then
    !                     pdf%atom(a1)%shell(sh)%y(i)=pdf%atom(a1)%shell(sh)%y(i)/(f0*f1)
    !                 else
    !                     pdf%atom(a1)%shell(sh)%y(i)=0.0_r8
    !                 endif
    !             enddo
    !         enddo
    !     enddo
    !
    !     ! Add things up to the total
    !     lo_allocate(pdf%x(pdf%nbin))
    !     lo_allocate(pdf%y(pdf%nbin))
    !     lo_allocate(dum(pdf%nbin))
    !     call lo_linspace(0.0_r8,cutoff,pdf%x)
    !     pdf%y=0.0_r8
    !     do a1=1,pdf%na
    !     do sh=1,pdf%atom(a1)%nshell
    !         if ( norm2(pm%atom(a1)%shell(sh)%prototype_vector) .lt. lo_tol ) cycle
    !         ii=pdf%atom(a1)%shell(sh)%pairind
    !         f0=lo_trapezoid_integration(pdf%atom(a1)%shell(sh)%x,pdf%atom(a1)%shell(sh)%y)
    !         call put_shape_on_new_mesh(pdf%atom(a1)%shell(sh)%x,pdf%atom(a1)%shell(sh)%y,pdf%x,dum)
    !         dum=dum*f0/lo_trapezoid_integration(pdf%x,dum)
    !         pdf%y=pdf%y+dum
    !         pdf%pair(ii)%y=pdf%pair(ii)%y+dum
    !     enddo
    !     enddo
    ! end block normalbin
    ! endif
end subroutine

!> Dump everything to file
subroutine write_to_hdf5(pdf, pm, uc, sim, filename)
    !> pair distribution function
    class(lo_pair_distribution), intent(in) :: pdf
    !> pair mapping
    type(lo_pairmapping), intent(in) :: pm
    !> unitcell
    type(lo_crystalstructure), intent(in) :: uc
    !> simulation
    type(lo_mdsim), intent(in) :: sim
    !> filename
    character(len=*), intent(in) :: filename

    type(lo_hdf5_helper) :: h5
    real(r8), dimension(:), allocatable :: dr0
    integer :: ish, i, nshell, ii

    call h5%init(__FILE__, __LINE__)
    call h5%open_file('write', trim(filename))

    ! Count number of relevant shells?
    nshell = 0
    do i = 1, pm%n_shell
        if (norm2(pm%sh(i)%r) .gt. lo_tol) then
            nshell = nshell + 1
        end if
    end do

    ! Store number of shells?
    call h5%store_attribute(nshell, h5%file_id, 'n_shell')

    ! Store time axis?
    allocate (dr0(sim%nt))
    dr0 = 0.0_r8
    do i = 1, sim%nt
        dr0(i) = (i - 1)*sim%timestep*lo_time_au_to_fs
    end do
    call h5%store_data(dr0, h5%file_id, 'time_axis', enhet='fs')
    deallocate (dr0)

    allocate (dr0(pdf%nbin))
    dr0 = 0.0_r8

    nshell = 0
    do ish = 1, pm%n_shell
        if (norm2(pm%sh(ish)%r) .lt. lo_tol) cycle
        nshell = nshell + 1
        ! store things per shell
        call h5%open_group('write', 'shell_'//tochar(nshell))

        ! store the prototype vector
        call h5%store_data(pm%sh(ish)%r*lo_bohr_to_A, h5%group_id, 'pair_vector', enhet='A')

        ! store histogram
        call h5%store_data(pdf%sh(ish)%x*lo_bohr_to_A, h5%group_id, 'x', enhet='A')
        call h5%store_data(pdf%sh(ish)%y/lo_bohr_to_A, h5%group_id, 'y', enhet='1/A')
        call h5%store_data(pdf%sh(ish)%z/lo_bohr_to_A, h5%group_id, 'z', enhet='1/A')

        ! create ideal peak at the average position
        dr0 = 0.0_r8
        do i = 1, pdf%nbin
            dr0(i) = lo_gauss(pdf%sh(ish)%x(i), norm2(pm%sh(ish)%r), pdf%sigma)
        end do
        dr0 = dr0/lo_bohr_to_A
        call h5%store_data(dr0, h5%group_id, 'y_ideal', enhet='1/A')

        call h5%close_group()
    end do

    deallocate (dr0)

    call h5%close_file()
    call h5%destroy(__FILE__, __LINE__)

    ! real(r8), dimension(:,:), allocatable :: dumx,dumy
    ! integer(HID_T) :: file_id
    ! integer :: a1,sh,nb,i
    ! character(len=1000) :: filename,dname
    ! character(len=4000) :: dstr

    ! filename='outfile.pair_distribution_function.hdf5'
    !
    ! call h5open_f(lo_status)
    ! call h5fcreate_f(trim(filename), H5F_ACC_TRUNC_F, file_id, lo_status)
    !
    ! ! store the total
    ! call lo_h5_store_data(pdf%x*lo_bohr_to_A,file_id,'r_axis',enhet='Angstrom')
    ! call lo_h5_store_data(pdf%y,file_id,'radial_pair_distribution_function')
    ! ! store the pair-resolved
    ! lo_allocate(dumy(pdf%nbin,pdf%npair))
    ! dumy=0.0_r8
    ! do i=1,pdf%npair
    !     dumy(:,i)=pdf%pair(i)%y
    ! enddo
    ! call lo_h5_store_data(dumy,file_id,'radial_pair_distribution_function_species_projected')
    ! lo_deallocate(dumy)
    ! ! Get a neat string for the legend.
    ! dstr=''
    ! ploop: do i=1,pdf%npair
    !     do a1=1,pdf%na
    !     do sh=1,pdf%atom(a1)%nshell
    !         if ( pdf%atom(a1)%shell(sh)%pairind .eq. i ) then
    !             dstr=trim(adjustl(dstr))//' '//trim(pdf%atom(a1)%label)//'-'//trim(pdf%atom(a1)%shell(sh)%label)
    !             cycle ploop
    !         endif
    !     enddo
    !     enddo
    ! enddo ploop
    ! call lo_h5_store_attribute(dstr,file_id,'projected_species_labels',lo_status)
    !
    ! ! store the projected
    ! if ( pdf%diffusion .eqv. .false. ) then
    !     call lo_h5_store_attribute(pdf%na,file_id,'number_unique_atoms',lo_status)
    !     do a1=1,pdf%na
    !         ! pack the data
    !         nb=size(pdf%atom(a1)%shell(1)%x,1)
    !         lo_allocate(dumx(pdf%atom(a1)%nshell-1,nb))
    !         lo_allocate(dumy(pdf%atom(a1)%nshell-1,nb))
    !         do sh=2,pdf%atom(a1)%nshell
    !             dumx(sh-1,:)=pdf%atom(a1)%shell(sh)%x*lo_bohr_to_A
    !             dumy(sh-1,:)=pdf%atom(a1)%shell(sh)%y
    !         enddo
    !
    !         dname='projected_r_axis_'//tochar(a1)
    !         call lo_h5_store_data(dumx,file_id,trim(dname),enhet='Angstrom')
    !         dname='projected_pair_distribution_function_atom_'//tochar(a1)
    !         call lo_h5_store_data(dumy,file_id,trim(dname))
    !
    !         ! Build a string that says what all the atoms are.
    !         dstr=''
    !         do i=1,pdf%atom(a1)%nshell
    !             dstr=trim(adjustl(dstr))//' '//trim(pdf%atom(a1)%label)//'-'//trim(pdf%atom(a1)%shell(i)%label)
    !         enddo
    !         dname='projected_pair_distribution_function_atomspec_'//tochar(a1)
    !         call lo_h5_store_attribute(trim(adjustl(dstr)),file_id,trim(dname))
    !
    !         lo_deallocate(dumx)
    !         lo_deallocate(dumy)
    !     enddo
    ! endif
    !
    ! call h5fclose_f(file_id, lo_status)
    ! call h5close_f(lo_status)

end subroutine

!> take a histogram from one set of x-y axes to another
subroutine put_shape_on_new_mesh(x, y, xi, yi)
    real(r8), dimension(:), intent(in) :: x, y
    real(r8), dimension(:), intent(in) :: xi
    real(r8), dimension(:), intent(out) :: yi
    !
    integer :: i, j, n, ni
    integer :: ii, jj
    real(r8) :: f0, sigma
    real(r8) :: ei, eimin, eimax, eirange, foursigma
    !
    n = size(x, 1)
    ni = size(xi, 1)
    yi = 0.0_r8
    sigma = (maxval(xi) - minval(xi))*3.0_r8/ni
    !
    eimin = minval(x)
    eimax = maxval(x)
    eirange = 1.0_r8/(eimax - eimin)
    foursigma = 4*sigma
    do i = 1, ni
        ! find lower and upper bound
        ei = xi(i)
        ii = floor(n*(ei - foursigma - eimin)*eirange)
        ii = max(ii, 1)
        ii = min(ii, n)
        !
        jj = ceiling(n*(ei + foursigma - eimin)*eirange)
        jj = min(jj, n)
        jj = max(jj, 1)
        do j = ii, jj
            f0 = lo_gauss(xi(i), x(j), sigma)
            yi(i) = yi(i) + f0*y(j)
        end do
    end do
end subroutine

!> get the trajectory for one pair
subroutine extract_trajectory(pm, ish, ipair, sim, trajectory)
    !> pair mapping
    type(lo_pairmapping), intent(in) :: pm
    !> coordination shell
    integer, intent(in) :: ish
    !> supercell pair
    integer, intent(in) :: ipair
    !> simulation data
    type(lo_mdsim), intent(in) :: sim
    !> trajectory
    real(r8), dimension(:, :), intent(out) :: trajectory

    real(r8), dimension(3) :: v0, v1
    integer :: a1, a2, t

    trajectory = 0.0_r8
    a1 = pm%sh(ish)%pairind(1, ipair)
    a2 = pm%sh(ish)%pairind(2, ipair)
    v0 = pm%sh(ish)%pairvec(:, ipair)
    do t = 1, sim%nt
        v1 = v0 + sim%u(:, a2, t) - sim%u(:, a1, t)
        trajectory(:, t) = v1
    end do

end subroutine

! !> get the trajectory for one pair
! subroutine extract_pair_trajectory(pm,a1,i1,sim,ss,trajectory)
!     !> pair mapping
!     type(lo_pairmapping), intent(in) :: pm
!     !> supercell atom
!     integer, intent(in) :: a1
!     !> supercell pair
!     integer, intent(in) :: i1
!     !> the simulation data
!     type(lo_mdsim), intent(in) :: sim
!     !> supercell
!     type(lo_crystalstructure), intent(in) :: ss
!     !> the trajectory
!     real(r8), dimension(:,:), intent(out) :: trajectory
!
!     real(r8), dimension(3) :: r0
!     integer :: ii,jj,t
!
!     ! lattice vector added to all the pair vectors
!     r0=pm%ssatom(a1)%pair(i1)%lv
!     ! indices to the supercell
!     ii=pm%ssatom(a1)%pair(i1)%i1
!     jj=pm%ssatom(a1)%pair(i1)%i2
!     ! the equilibrium pair distance
!     r0=r0+ss%rcart(:,jj)-ss%rcart(:,ii)
!     ! this pair is supposed to be lv + ss(jj) - ss(ii)
!     do t=1,sim%nt
!         trajectory(:,t)=r0+sim%u(:,jj,t)-sim%u(:,ii,t)
!     enddo
! end subroutine

end module
