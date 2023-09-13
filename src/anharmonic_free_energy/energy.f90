module energy
!! get the anharmonic free energy
use konstanter, only: r8,i8,lo_twopi,lo_freqtol,lo_imag,lo_sqtol,lo_status
use gottochblandat, only: lo_stop_gracefully,tochar,walltime,lo_chop,lo_trueNtimes,&
                   lo_progressbar_init,lo_progressbar,lo_planck,lo_trapezoid_integration,open_file,&
                   lo_flattentensor,lo_sqnorm,lo_linspace,lo_mean,lo_clean_fractional_coordinates
use mpi_wrappers, only: lo_mpi_helper,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_DOUBLE_COMPLEX,MPI_IN_PLACE
use lo_memtracker, only: lo_mem_helper
use type_crystalstructure, only: lo_crystalstructure
use type_forceconstant_secondorder, only: lo_forceconstant_secondorder
use type_forceconstant_thirdorder, only: lo_forceconstant_thirdorder
use type_forceconstant_fourthorder, only: lo_forceconstant_fourthorder
use type_qpointmesh, only: lo_fft_mesh
use type_phonon_dispersions, only: lo_phonon_dispersions

implicit none
private

public :: perturbative_anharmonic_free_energy

contains

!> Calculates the anharmonic contributions to the free energy
subroutine perturbative_anharmonic_free_energy(p,fct,fcf,qp,dr,temperature,free_energy_thirdorder,free_energy_fourthorder,mw,mem,verbosity)
    !> crystal structure
    type(lo_crystalstructure), intent(in) :: p
    !> third order force constant
    type(lo_forceconstant_thirdorder), intent(in) :: fct
    !> fourth order force constant
    type(lo_forceconstant_fourthorder), intent(in) :: fcf
    !> q-point mesh
    type(lo_fft_mesh), intent(in) :: qp
    !> dispersions
    type(lo_phonon_dispersions), intent(in) :: dr
    !> temperature
    real(r8) :: temperature
    !> free energies
    real(r8), intent(out) :: free_energy_thirdorder,free_energy_fourthorder
    !> mpi helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem
    !> talk a lot?
    integer, intent(in) :: verbosity

    complex(r8), dimension(:,:,:), allocatable :: nuvec1,nuvec2
    real(r8), dimension(:,:), allocatable :: en3,en4
    real(r8) :: t0,t1,timer

    ! Start timers
    timer=walltime()
    t0=timer
    t1=timer

    ! Set basic things
    init: block
        real(r8) :: f0,f1
        integer :: iq,imode,ctr,iatom,ix,ialpha

        if ( verbosity .gt. 0 ) then
            write(*,*) ''
            write(*,*) 'CALCULATING ANHARMONIC FREE ENERGY'
        endif

        ! Get the scaled eigenvectors
        call mem%allocate(nuvec1,[dr%n_mode,dr%n_mode,qp%n_irr_point],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%allocate(nuvec2,[dr%n_mode,dr%n_mode,qp%n_full_point],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        nuvec1=0.0_r8
        nuvec2=0.0_r8

        ctr=0
        do iq=1,qp%n_irr_point
        do imode=1,dr%n_mode
            ctr=ctr+1
            if ( mod(ctr,mw%n) .ne. mw%r ) cycle
            if ( dr%iq(iq)%omega(imode) .gt. lo_freqtol ) then
                f0=1.0_r8/sqrt(dr%iq(iq)%omega(imode))
            else
                f0=0.0_r8
            endif
            do iatom=1,p%na
                f1=p%invsqrtmass(iatom)
                do ix=1,3
                    ialpha=(iatom-1)*3+ix
                    nuvec1(ialpha,imode,iq)=dr%iq(iq)%egv(ialpha,imode)*f0*f1
                enddo
            enddo
        enddo
        enddo
        do iq=1,qp%n_full_point
        do imode=1,dr%n_mode
            ctr=ctr+1
            if ( mod(ctr,mw%n) .ne. mw%r ) cycle
            if ( dr%aq(iq)%omega(imode) .gt. lo_freqtol ) then
                f0=1.0_r8/sqrt(dr%aq(iq)%omega(imode))
            else
                f0=0.0_r8
            endif
            do iatom=1,p%na
                f1=p%invsqrtmass(iatom)
                do ix=1,3
                    ialpha=(iatom-1)*3+ix
                    nuvec2(ialpha,imode,iq)=dr%aq(iq)%egv(ialpha,imode)*f0*f1
                enddo
            enddo
        enddo
        enddo
        call mw%allreduce('sum',nuvec1)
        call mw%allreduce('sum',nuvec2)

        ! Space for free energies per mode
        call mem%allocate(en3,[dr%n_mode,qp%n_irr_point],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%allocate(en4,[dr%n_mode,qp%n_irr_point],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        en3=0.0_r8
        en4=0.0_r8

        if ( verbosity .gt. 0 ) then
            t1=walltime()
            write(*,*) '... scaled vectors and made space (',tochar(t1-t0),'s)'
            t0=t1
        endif
    end block init

    ! Integrate out the free energy
    third: block
        real(r8), external :: zdotu
        real(r8), parameter :: threephonon_prefactor=-1.0_r8/48.0_r8
        real(r8), parameter :: fourphonon_prefactor=1.0_r8/32.0_r8
        real(r8), parameter :: onethird=1.0_r8/3.0_r8
        complex(r8), dimension(:,:), allocatable :: buf_ev1,buf_ev2,buf_ev3
        complex(r8), dimension(:), allocatable :: evp1,evp2,evp3,ptf
        complex(r8) :: c0
        real(r8), dimension(3) :: qv1,qv2,qv3
        real(r8) :: psisq,sigma,om1,om2,om3,n1,n2,n3,f1,f2,prefactor,s1,s2,s3
        integer :: ctr,b1,b2,b3,q1,q2,q3

        ! Space for intermediate products
        call mem%allocate(evp1,dr%n_mode**2,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%allocate(evp2,dr%n_mode**3,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%allocate(ptf ,dr%n_mode**3,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%allocate(buf_ev1,[dr%n_mode,dr%n_mode],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%allocate(buf_ev2,[dr%n_mode,dr%n_mode],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%allocate(buf_ev3,[dr%n_mode,dr%n_mode],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        !call mem%allocate(buf_ev4,[dr%n_mode,dr%n_mode],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        evp1=0.0_r8
        evp2=0.0_r8
        ptf =0.0_r8
        buf_ev1=0.0_r8
        buf_ev2=0.0_r8
        buf_ev3=0.0_r8
        !buf_ev4=0.0_r8

        ! Optimizations:
        ! 1) check the prefactors, can maybe do those before in a neater way
        ! 2) only have to go over half the q-vectors I think.

        if ( verbosity .gt. 0 ) call lo_progressbar_init()
        en3=0.0_r8
        ctr=0
        do q1=1,qp%n_irr_point
        do q2=1,qp%n_full_point
            ! make it parallel
            ctr=ctr+1
            if ( mod(ctr,mw%n) .ne. mw%r ) cycle
            ! locate third q-vector
            q3=fft_third_grid_index( qp%ip(q1)%full_index, q2, qp%griddensity)

            qv1=qp%ip(q1)%r
            qv2=qp%ap(q2)%r
            qv3=qp%ap(q3)%r

            prefactor=threephonon_prefactor*qp%ap(q2)%integration_weight

            ! pre-transform the matrix element
            call pretransform_phi3(fct,qv2,qv3,ptf)
            ! pre-fetch eigenvectors. Maybe a speedup, dunno.
            buf_ev1=nuvec1(:,:,q1)
            buf_ev2=nuvec2(:,:,q2)
            buf_ev3=nuvec2(:,:,q3)

            do b1=1,dr%n_mode
                do b2=1,dr%n_mode
                    evp1=0.0_r8
                    call zgeru(dr%n_mode, dr%n_mode, (1.0_r8,0.0_r8), buf_ev2(:,b2), 1, buf_ev1(:,b1), 1, evp1, dr%n_mode)
                    do b3=1,dr%n_mode
                        evp2=0.0_r8
                        call zgeru(dr%n_mode, dr%n_mode*dr%n_mode, (1.0_r8,0.0_r8), buf_ev3(:,b3), 1, evp1, 1, evp2, dr%n_mode)
                        evp2=conjg(evp2)
                        c0=dot_product(evp2,ptf)
                        !c0=zdotu(dr%n_mode**3,evp2,1,ptf,1)
                        ! And now we have the matrix element.
                        psisq=real(conjg(c0)*c0,r8)
                        ! Get a smearing parameter
                        s1=dr%default_smearing(b1)
                        s2=dr%default_smearing(b2)
                        s3=dr%default_smearing(b3)
                        sigma=sqrt(s1**2+s2**2+s3**2)
                        ! Same sigma as the usual case for debugging.
                        !sigma=lo_mean(dr%default_smearing)
                        om1=dr%iq(q1)%omega(b1)
                        om2=dr%aq(q2)%omega(b2)
                        om3=dr%aq(q3)%omega(b3)

                        n1=lo_planck(temperature,om1)
                        n2=lo_planck(temperature,om2)
                        n3=lo_planck(temperature,om3)

                        ! This is the Wallace expression
                        ! f1=n1*n2+n1+onethird
                        ! f1=3*f1*real(1.0_r8/( om1+om2+om3+lo_imag*sigma ))
                        ! f2=2*n1*n3-n1*n2+n3
                        ! f2=3*f2*real(1.0_r8/( om1+om2-om3+lo_imag*sigma ))
                        ! en3(b1,q1)=en3(b1,q1)+( (f1+f2)*psisq )*prefactor

                        ! Try the Cowley expression instead?
                        f1=(n1+1)*(n1+n3+1)+n2*n3
                        f1=f1*real(1.0_r8/( om1+om2+om3+lo_imag*sigma ))
                        f2=n1*n2+n1*n3-n2*n3+1
                        f2=3*f2*real(1.0_r8/( om1+om2-om3+lo_imag*sigma ))
                        en3(b1,q1)=en3(b1,q1)+( (f1+f2)*psisq )*prefactor
                    enddo
                enddo
            enddo
        enddo
        if ( verbosity .gt. 0 .and. q1 .lt. qp%n_irr_point ) then
            call lo_progressbar(' ... third order',q1,qp%n_irr_point,walltime()-t0)
        endif
        enddo

        ! Some cleanup
        call mem%deallocate(evp1,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%deallocate(evp2,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%deallocate(ptf ,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)

        if ( verbosity .gt. 0 ) then
            t1=walltime()
            call lo_progressbar(' ... third order',qp%n_irr_point,qp%n_irr_point,t1-t0)
            t0=t1
        endif

        ! Do the fourth order? Might as well.
        call mem%allocate(evp1,dr%n_mode**2,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%allocate(evp2,dr%n_mode**2,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%allocate(evp3,dr%n_mode**4,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%allocate(ptf ,dr%n_mode**4,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        evp1=0.0_r8
        evp2=0.0_r8
        evp3=0.0_r8
        ptf =0.0_r8

        if ( verbosity .gt. 0 ) call lo_progressbar_init()
        en4=0.0_r8
        ctr=0
        do q1=1,qp%n_irr_point
        do q2=1,qp%n_full_point
            ! make it parallel
            ctr=ctr+1
            if ( mod(ctr,mw%n) .ne. mw%r ) cycle
            qv1=qp%ip(q1)%r
            qv2=qp%ap(q2)%r
            prefactor=fourphonon_prefactor*qp%ap(q2)%integration_weight

            ! pre-transform the matrix element
            call pretransform_phi4(fcf,qv1,qv2,ptf)
            ! pre-fetch eigenvectors. Maybe a speedup, dunno.
            buf_ev1=nuvec1(:,:,q1)
            buf_ev2=nuvec2(:,:,q2)
            do b1=1,dr%n_mode
            do b2=1,dr%n_mode
                om1=dr%iq(q1)%omega(b1)
                om2=dr%aq(q2)%omega(b2)
                if ( om1 .lt. lo_freqtol ) cycle
                if ( om2 .lt. lo_freqtol ) cycle
                n1=lo_planck(temperature,om1)
                n2=lo_planck(temperature,om2)
                ! Now to get
                evp1=0.0_r8
                evp2=0.0_r8
                evp3=0.0_r8
                call zgerc(dr%n_mode, dr%n_mode, (1.0_r8,0.0_r8), buf_ev1(:,b1), 1, buf_ev1(:,b1), 1, evp1, dr%n_mode)
                call zgerc(dr%n_mode, dr%n_mode, (1.0_r8,0.0_r8), buf_ev2(:,b2), 1, buf_ev2(:,b2), 1, evp2, dr%n_mode)
                call zgeru(dr%n_mode**2, dr%n_mode**2, (1.0_r8,0.0_r8), evp2, 1, evp1, 1, evp3, dr%n_mode**2)
                evp3=conjg(evp3)
                psisq=real(dot_product(evp3,ptf),r8)
                f1=(2*n1+1)*(2*n2+1)*psisq*prefactor

                en4(b1,q1)=en4(b1,q1)+f1
            enddo
            enddo
        enddo
        if ( verbosity .gt. 0 .and. q1 .lt. qp%n_irr_point ) then
            call lo_progressbar(' ... fourth order',q1,qp%n_irr_point,walltime()-t0)
        endif
        enddo

        call mem%deallocate(evp1,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%deallocate(evp2,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%deallocate(evp3,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%deallocate(ptf ,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%deallocate(buf_ev1,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%deallocate(buf_ev2,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%deallocate(buf_ev3,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)

        if ( verbosity .gt. 0 ) then
            t1=walltime()
            call lo_progressbar(' ... fourth order',qp%n_irr_point,qp%n_irr_point,t1-t0)
            t0=t1
        endif
    end block third

    ! add things up
    finalize: block
        integer :: q1

        ! Add it up
        call mw%allreduce('sum',en4)
        call mw%allreduce('sum',en3)
        free_energy_thirdorder=0.0_r8
        free_energy_fourthorder=0.0_r8
        do q1=1,qp%n_irr_point
            free_energy_thirdorder=free_energy_thirdorder+sum(en3(:,q1))*qp%ip(q1)%integration_weight
            free_energy_fourthorder=free_energy_fourthorder+sum(en4(:,q1))*qp%ip(q1)%integration_weight
        enddo
        ! And normalize it to be per atom
        free_energy_thirdorder=free_energy_thirdorder/p%na
        free_energy_fourthorder=free_energy_fourthorder/p%na

        ! and cleanup
        call mem%deallocate(nuvec1,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%deallocate(nuvec2,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%deallocate(en3   ,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%deallocate(en4   ,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    end block finalize

    if ( verbosity .gt. 0 ) then
        t1=walltime()
        write(*,*) '... got anharmonic free energy (',tochar(t1-timer),'s)'
    endif

    ! !> and the fourth order things
    ! fourth: block
    !     real(r8), parameter :: fourphonon_prefactor=1.0_r8/32.0_r8
    !     complex(r8) :: c0
    !     !real(r8), dimension(:,:,:), allocatable :: psibuf
    !     real(r8), dimension(4) :: omega
    !     real(r8), dimension(3) :: qv1,qv2,qv3,qv4
    !     real(r8) :: psi,prefactor,n1,n2,f1
    !     integer :: q1,q2,b1,b2,b3,i,j,l
    !     ! buffers for speed
    !     complex(r8), dimension(:,:), allocatable :: egv,flat_cwrk
    !     real(r8), dimension(:,:), allocatable :: flat_m,flat_r2,flat_r3,flat_r4
    !     real(r8), dimension(dr%n_mode) :: dumom1,dumom3,dump1,dump3
    !     real(r8) :: omthres
    !     integer(i8) :: ctrlong
    !     integer, dimension(:,:), allocatable :: flat_ai
    !     integer :: nquartet
    !
    !
    !     allocate(en4(dr%n_mode,qp%n_irr_point))
    !     en4=0.0_r8
    !
    !     omthres=dr%omega_min*0.5_r8
    !     ! create some buffers
    !     nquartet=0
    !     do i=1,fcf%na
    !         nquartet=nquartet+fcf%atom(i)%n
    !     enddo
    !     allocate(flat_m(81,nquartet))
    !     allocate(flat_r2(3,nquartet))
    !     allocate(flat_r3(3,nquartet))
    !     allocate(flat_r4(3,nquartet))
    !     allocate(flat_cwrk(81,fcf%na**4))
    !     allocate(flat_ai(4,nquartet))
    !     l=0
    !     do i=1,fcf%na
    !     do j=1,fcf%atom(i)%n
    !         l=l+1
    !         flat_r2(:,l)=fcf%atom(i)%quartet(j)%lv2
    !         flat_r3(:,l)=fcf%atom(i)%quartet(j)%lv3
    !         flat_r4(:,l)=fcf%atom(i)%quartet(j)%lv4
    !         flat_m(:,l)=lo_flattentensor( fcf%atom(i)%quartet(j)%mwm )
    !         flat_ai(:,l)=[fcf%atom(i)%quartet(j)%i1,fcf%atom(i)%quartet(j)%i2,fcf%atom(i)%quartet(j)%i3,fcf%atom(i)%quartet(j)%i4]
    !     enddo
    !     enddo
    !     allocate(egv(dr%n_mode,4))
    !     egv=0.0_r8
    !     ! the unit
    !     !prefactor=1.0_r8/32.0_r8
    !
    !     ! start integrating
    !     ctrlong=0
    !     do q1=1,qp%n_irr_point
    !     do q2=1,qp%n_full_point
    !         ! q-vectors
    !         qv1= qp%ip(q1)%r*lo_twopi
    !         qv2=-qp%ip(q1)%r*lo_twopi
    !         qv3= qp%ap(q2)%r*lo_twopi
    !         qv4=-qp%ap(q2)%r*lo_twopi
    !         ! store omega and planck-factors in a buffer
    !         do b1=1,dr%n_mode
    !             dumom1(b1)=dr%iq(q1)%omega(b1)
    !             dumom3(b1)=dr%aq(q2)%omega(b1)
    !             dump1(b1)=lo_planck(temperature,dumom1(b1))
    !             dump3(b1)=lo_planck(temperature,dumom3(b1))
    !         enddo
    !
    !         prefactor=fourphonon_prefactor*qp%ap(q2)%integration_weight
    !
    !         ! technically four band-loops, but b1=b2 and b3=b4
    !         do b1=1,dr%n_mode
    !         ctrlong=ctrlong+1; if ( mod(ctrlong,int(mw%n,r8)) .ne. mw%r ) cycle
    !         do b3=1,dr%n_mode
    !             omega=[dumom1(b1),dumom1(b1),dumom3(b3),dumom3(b3)]
    !             if ( minval(omega) .lt. omthres ) cycle
    !             egv(:,1)=dr%iq(q1)%egv(:,b1)
    !             egv(:,2)=conjg(dr%iq(q1)%egv(:,b1))
    !             egv(:,3)=dr%aq(q2)%egv(:,b3)
    !             egv(:,4)=conjg(dr%aq(q2)%egv(:,b3))
    !             call loc_4phpsi(fcf%na,nquartet,flat_r2,flat_r3,flat_r4,flat_m,flat_cwrk,flat_ai,omega,egv,qv2,qv3,qv4,c0)
    !             psi=real(c0,r8)
    !
    !             n1=dump1(b1)
    !             n2=dump3(b3)
    !             f1=(2*n1+1)*(2*n2+1)*psi*prefactor
    !
    !
    !             en4(b1,q1)=en4(b1,q1)+f1
    !         enddo
    !         enddo
    !     enddo
    !     if ( mw%talk ) write(*,*) 'did',q1,qp%n_irr_point
    !     enddo
    !
    !     if ( mw%talk ) then
    !         write(*,*) 'fourth new:',free_energy_fourthorder
    !     endif
    !
    !     call mw%allreduce('sum',en4)
    !     free_energy_fourthorder=0.0_r8
    !     do q1=1,qp%n_irr_point
    !         free_energy_fourthorder=free_energy_fourthorder+sum(en4(:,q1))*qp%ip(q1)%integration_weight
    !     enddo
    !     free_energy_fourthorder=free_energy_fourthorder/p%na
    !
    !     if ( mw%talk ) then
    !         write(*,*) 'fourth old:',free_energy_fourthorder
    !     endif
    ! end block fourth
end subroutine

!> pre-transform to get half of the third order matrix element
subroutine pretransform_phi3(fct,q2,q3,ptf)
    !> third order forceconstant
    type(lo_forceconstant_thirdorder), intent(in) :: fct
    !> q-vectors
    real(r8), dimension(3), intent(in) :: q2,q3
    !> flattened, pretransformed matrix element
    complex(r8), dimension(:), intent(out) :: ptf

    integer :: i,j,k,l

    complex(r8) :: expiqr
    real(r8), dimension(3) :: rv2,rv3
    real(r8) :: iqr
    integer :: a1,a2,a3,ia,ib,ic,t,nb

    nb=fct%na*3
    ptf=0.0_r8
    do a1=1,fct%na
    do t=1,fct%atom(a1)%n
        a2=fct%atom(a1)%triplet(t)%i2
        a3=fct%atom(a1)%triplet(t)%i3

        rv2=fct%atom(a1)%triplet(t)%lv2
        rv3=fct%atom(a1)%triplet(t)%lv3

        iqr=dot_product(q2,rv2)+dot_product(q3,rv3)
        iqr=-iqr*lo_twopi
        expiqr=cmplx(cos(iqr),sin(iqr),r8)
        do i=1,3
        do j=1,3
        do k=1,3
            ia=(a1-1)*3+i
            ib=(a2-1)*3+j
            ic=(a3-1)*3+k
            ! Now for the grand flattening scheme, consistent with the zgeru operations above.
            l=(ia-1)*nb*nb + (ib-1)*nb + ic
            ptf(l)=ptf(l)+fct%atom(a1)%triplet(t)%m(i,j,k)*expiqr
        enddo
        enddo
        enddo
    enddo
    enddo
end subroutine

!> pre-transform to get half of the fourth order matrix element
subroutine pretransform_phi4(fcf,q1,q2,ptf)
    !> third order forceconstant
    type(lo_forceconstant_fourthorder), intent(in) :: fcf
    !> q-vectors
    real(r8), dimension(3), intent(in) :: q1,q2
    !> flattened, pretransformed matrix element
    complex(r8), dimension(:), intent(out) :: ptf

    integer :: i,j,k,l,m

    complex(r8) :: expiqr
    real(r8), dimension(3) :: rv2,rv3,rv4
    real(r8) :: iqr
    integer :: a1,a2,a3,a4,ia,ib,ic,id,q,nb

    nb=fcf%na*3
    ptf=0.0_r8
    do a1=1,fcf%na
    do q=1,fcf%atom(a1)%n
        a2=fcf%atom(a1)%quartet(q)%i2
        a3=fcf%atom(a1)%quartet(q)%i3
        a4=fcf%atom(a1)%quartet(q)%i4

        rv2=fcf%atom(a1)%quartet(q)%lv2
        rv3=fcf%atom(a1)%quartet(q)%lv3
        rv4=fcf%atom(a1)%quartet(q)%lv4

        iqr=-dot_product(q1,rv2)+dot_product(q2,rv3)-dot_product(q2,rv4)
        iqr=-iqr*lo_twopi
        expiqr=cmplx(cos(iqr),sin(iqr),r8)
        do l=1,3
        do k=1,3
        do j=1,3
        do i=1,3
            ia=(a1-1)*3+i
            ib=(a2-1)*3+j
            ic=(a3-1)*3+k
            id=(a4-1)*3+l
            ! Now for the grand flattening scheme, consistent with the zgeru operations above.
            m=(ia-1)*nb*nb*nb + (ib-1)*nb*nb + (ic-1)*nb + id
            ptf(m)=ptf(m) + fcf%atom(a1)%quartet(q)%m(i,j,k,l)*expiqr
        enddo
        enddo
        enddo
        enddo
    enddo
    enddo
end subroutine

!> returns the index on the grid that gives q3=-q1-q2
pure function fft_third_grid_index(i1,i2,dims) result(i3)
    !> index to q1
    integer, intent(in) :: i1
    !> index to q2
    integer, intent(in) :: i2
    !> dimensions of the grid
    integer, dimension(3), intent(in) :: dims
    !> index to q3
    integer :: i3

    integer, dimension(3) :: gi1,gi2,gi3
    integer :: l,k

    ! Convert triplet to singlet
    gi1=singlet_to_triplet(i1,dims(2),dims(3))
    gi2=singlet_to_triplet(i2,dims(2),dims(3))
    do l=1,3
        gi3(l)=3-gi1(l)-gi2(l)
    enddo
    do k=1,3
    do l=1,3
        if ( gi3(l) .lt. 1 ) gi3(l)=gi3(l)+dims(l)
        if ( gi3(l) .gt. dims(l) ) gi3(l)=gi3(l)-dims(l)
    enddo
    enddo
    ! convert it back to a singlet
    i3=triplet_to_singlet(gi3,dims(2),dims(3))

    contains
    !> convert a linear index to a triplet
    pure function singlet_to_triplet(l,ny,nz) result(gi)
        !> linear index
        integer, intent(in) :: l
        !> second dimension
        integer, intent(in) :: ny
        !> third dimension
        integer, intent(in) :: nz
        !> grid-index
        integer, dimension(3) :: gi

        integer :: i,j,k

        k=mod(l,nz)
        if ( k .eq. 0 ) k=nz
        j=mod( (l-k)/nz , ny )+1
        i=(l-k-(j-1)*nz)/(nz*ny)+1
        gi=[i,j,k]
    end function
    !> convert a triplet index to a singlet
    pure function triplet_to_singlet(gi,ny,nz) result(l)
        !> grid-index
        integer, dimension(3), intent(in) :: gi
        !> second dimension
        integer, intent(in) :: ny
        !> third dimension
        integer, intent(in) :: nz
        !> linear index
        integer :: l

        l=(gi(1)-1)*ny*nz+(gi(2)-1)*nz+gi(3)
    end function
end function

! subroutine loc_4phpsi(na,nquartet,flv2,flv3,flv4,fm,fcw,fai,omega,egv,q2,q3,q4,psi)
!     !> temporary work arrays, that actually is the force constant
!     integer, intent(in) :: na,nquartet
!     real(r8), dimension(3,nquartet), intent(in) :: flv2,flv3,flv4
!     real(r8), dimension(81,nquartet), intent(in) :: fm
!     complex(r8), dimension(81,na*na*na*na), intent(inout) :: fcw
!     integer, dimension(4,nquartet), intent(in) :: fai
!     !> four frequencies
!     real(r8), dimension(4), intent(in) :: omega
!     !> four eigenvectors
!     complex(r8), dimension(na*3,4), intent(in) :: egv
!     !> three q-points
!     real(r8), dimension(3), intent(in) :: q2,q3,q4
!     !> the matrix element
!     complex(r8), intent(out) :: psi
!
!     complex(r8) :: expiqr,c0
!     real(r8), dimension(3) :: rv2,rv3,rv4
!     real(r8) :: iqr,omegaprod
!     real(r8), parameter :: enhet=1.0_r8 !5.8104841317859603E074_r8
!     integer :: atom1,atom2,atom3,atom4,i,j,k,l,k1,k2,k3,k4,ii,jj
!     integer :: na3,na2
!
!     na3=na*na*na
!     na2=na*na
!     ! Eigenvector product, again
!     do atom1=1,na
!     do atom2=1,na
!     do atom3=1,na
!     do atom4=1,na
!         ii=(atom1-1)*na3+(atom2-1)*na2+(atom3-1)*na+atom4
!         do i=1,3
!         do j=1,3
!         do k=1,3
!         do l=1,3
!             jj=(i-1)*27+(j-1)*9+(k-1)*3+l
!             k1=(atom1-1)*3+i
!             k2=(atom2-1)*3+j
!             k3=(atom3-1)*3+k
!             k4=(atom4-1)*3+l
!             fcw(jj,ii)=egv(k1,1)*egv(k2,2)*egv(k3,3)*egv(k4,4)
!         enddo
!         enddo
!         enddo
!         enddo
!     enddo
!     enddo
!     enddo
!     enddo
!
!     ! Frequency product
!     omegaprod=omega(1)*omega(2)*omega(3)*omega(4)
!     omegaprod=sqrt(omegaprod)
!     ! The actual matrix element
!     c0=0.0_r8
!     do i=1,nquartet
!         atom1=fai(1,i)
!         atom2=fai(2,i)
!         atom3=fai(3,i)
!         atom4=fai(4,i)
!         ii=(atom1-1)*na3+(atom2-1)*na2+(atom3-1)*na+atom4
!         rv2=flv2(:,i)
!         rv3=flv3(:,i)
!         rv4=flv4(:,i)
!         iqr=dot_product(q2,rv2)+dot_product(q3,rv3)+dot_product(q4,rv4)
!         expiqr=cmplx(cos(iqr),sin(iqr),r8)
!         c0=c0+sum(fcw(:,ii)*fm(:,i))*expiqr
!     enddo
!     psi=enhet*c0/omegaprod
! end subroutine

! !> impossibly convoluted way of calculating the anharmonic free energy. Really parallel though!
! subroutine impossible_free_energy(uc,fc,fct,fcf,qp,dr,temperature,free_energy_thirdorder,free_energy_fourthorder,mw)
!     !> crystal structure
!     type(lo_crystalstructure), intent(in) :: uc
!     !> second order force constant
!     type(lo_forceconstant_secondorder), intent(inout) :: fc
!     !> third order force constant
!     type(lo_forceconstant_thirdorder), intent(in) :: fct
!     !> fourth order force constant
!     type(lo_forceconstant_fourthorder), intent(in) :: fcf
!     !> q-point mesh
!     type(lo_fft_mesh), intent(in) :: qp
!     !> dispersions
!     type(lo_phonon_dispersions), intent(in) :: dr
!     !> temperature
!     real(r8) :: temperature
!     !> free energies
!     real(r8), intent(out) :: free_energy_thirdorder,free_energy_fourthorder
!     !> mpi helper
!     type(lo_mpi_helper) :: mw
!
!     real(r8), dimension(:,:), allocatable :: en3,en4
!     real(r8) :: omthres,sigma
!     integer :: index_gamma
!     real(r8) :: t0,t1
!
!     t0=walltime()
!     t1=t0
!
!     init: block
!         integer :: i
!         sigma=lo_mean(dr%default_smearing)
!         omthres=dr%omega_min*0.5_r8
!         index_gamma=-1
!         do i=1,qp%n_irr_point
!             if ( lo_sqnorm(qp%ip(i)%r) .lt. lo_sqtol ) then
!                 index_gamma=i
!                 exit
!             endif
!         enddo
!         lo_allocate(en3(dr%n_mode,qp%n_irr_point))
!         lo_allocate(en4(dr%n_mode,qp%n_irr_point))
!         en3=0.0_r8
!         en4=0.0_r8
!
!         if ( mw%talk ) then
!             t1=walltime()
!             write(*,*) 'init:',t1-t0
!             t0=t1
!         endif
!     end block init
!
!     ! Start with the third order
!     third: block
!         complex(r8) :: c0
!         real(r8), dimension(3) :: v0,qv1,qv2,qv3,omega
!         real(r8) :: psisquare,doublepsi,prefactor
!         real(r8) :: n1,n2,n3,f0,f1,f2,f3
!         integer, dimension(3) :: gi
!         integer :: q1,q2,q3,iq2,b1,b2,b3,ctr
!         ! buffers for speed
!         complex(r8), dimension(:,:,:,:,:,:), allocatable :: cwrk,cdwrk
!         real(r8), dimension(dr%n_mode) :: pbuf1,pbuf2,pbuf3,ombuf1,ombuf2,ombuf3
!         real(r8), parameter :: onethird=1.0_r8/3.0_r8
!
!         ! work-array for the scattering amplitude
!         allocate(cwrk(3,3,3,fct%na,fct%na,fct%na))  ! work-array for scattering amplitude
!         allocate(cdwrk(3,3,3,fct%na,fct%na,fct%na)) ! work-array for scattering amplitude
!
!         ctr=0
!         prefactor=-1.0_r8/48.0_r8
!         do q1=1,qp%n_irr_point
!         do q2=1,qp%n_full_point
!             ! locate the third q-vector
!             v0=uc%cartesian_to_fractional(-qp%ip(q1)%r-qp%ap(q2)%r,reciprocal=.true.)
!             gi=qp%index_from_coordinate(v0)
!             q3=qp%gridind2ind(gi(1),gi(2),gi(3))
!             ! buffer things
!             do b1=1,dr%n_mode
!                 ombuf1(b1)=dr%iq(q1)%omega(b1)
!                 ombuf2(b1)=dr%aq(q2)%omega(b1)
!                 ombuf3(b1)=dr%aq(q3)%omega(b1)
!                 pbuf1(b1)=lo_planck(temperature,ombuf1(b1))
!                 pbuf2(b1)=lo_planck(temperature,ombuf2(b1))
!                 pbuf3(b1)=lo_planck(temperature,ombuf3(b1))
!             enddo
!             ! get the vectors, and add 2*pi
!             qv1=-qp%ip(q1)%r*lo_twopi
!             qv2=-qp%ap(q2)%r*lo_twopi
!             qv3=-qp%ap(q3)%r*lo_twopi
!             ! do the MPI-division here. Should make it really parallel.
!             ! should also mean minimal overhead.
!             do b1=1,dr%n_mode
!             ctr=ctr+1; if ( mod(ctr,mw%n) .ne. mw%r ) cycle
!             do b2=1,dr%n_mode
!             do b3=1,dr%n_mode
!                 omega(1)=ombuf1(b1)
!                 omega(2)=ombuf2(b2)
!                 omega(3)=ombuf3(b3)
!                 if ( minval(omega) .gt. omthres ) then
!                     ! normal three-phonon thing, matrix element
!                     call loc_3phpsi(fct,omega,dr%iq(q1)%egv(:,b1),dr%aq(q2)%egv(:,b2),dr%aq(q3)%egv(:,b3),qv2,qv3,cwrk,c0)
!                     psisquare=abs(conjg(c0)*c0)
!                     n1=pbuf1(b1)
!                     n2=pbuf2(b2)
!                     n3=pbuf3(b3)
!                     ! prefactors for the other stuff
!                     f1=n1*n2+n1+onethird
!                     f1=3*f1*real(1.0_r8/( omega(1)+omega(2)+omega(3)+lo_imag*sigma ))
!                     f2=2*n1*n3-n1*n2+n3
!                     f2=3*f2*real(1.0_r8/( omega(1)+omega(2)-omega(3)+lo_imag*sigma ))
!                 else
!                     f1=0.0_r8
!                     f2=0.0_r8
!                     psisquare=0.0_r8
!                 endif
!                 ! get the double-psi term
!                 !omega(3)=dr%iq(index_gamma)%omega(b3)
!                 !if ( minval(omega) .gt. omthres ) then
!                 !    call loc_3phdbl(fct,omega,dr%iq(q1)%egv(:,b1),dr%aq(q2)%egv(:,b2),dr%aq(index_gamma)%egv(:,b3),qv1,qv2,cwrk,cdwrk,c0)
!                 !    doublepsi=real(c0)
!                 !    f3=6*(n1*(n2+1)+0.25_r8)
!                 !    f3=f3*real(1.0_r8/( omega(3)+lo_imag*sigma ))
!                 !else
!                 !    doublepsi=0.0_r8
!                 !    f3=0.0_r8
!                 !endif
!                 doublepsi=0.0_r8
!                 f3=0.0_r8
!                 ! Add it up
!                 !psisquare=1.0_r8
!                 f0=( (f1+f2)*psisquare+f3*doublepsi )*prefactor
!                 en3(b1,q1)=en3(b1,q1)+f0
!             enddo
!             enddo
!             enddo
!         enddo
!         enddo
!
!         if ( mw%talk ) then
!             t1=walltime()
!             write(*,*) 'third:',t1-t0
!             t0=t1
!         endif
!
!     end block third
!
!     !> and the fourth order things
!     fourth: block
!         complex(r8) :: c0
!         real(r8), dimension(:,:,:), allocatable :: psibuf
!         real(r8), dimension(4) :: omega
!         real(r8), dimension(3) :: qv1,qv2,qv3,qv4
!         real(r8) :: psi,prefactor,n1,n2,f1
!         integer :: q1,q2,b1,b2,b3,b4,i,j,l
!         ! buffers for speed
!         complex(r8), dimension(:,:), allocatable :: egv,flat_cwrk
!         real(r8), dimension(:,:), allocatable :: flat_m,flat_r2,flat_r3,flat_r4
!         real(r8), dimension(dr%n_mode) :: dumom1,dumom3,dump1,dump3
!         integer(i8) :: ctrlong
!         integer, dimension(:,:), allocatable :: flat_ai
!         integer :: nquartet
!
!         ! create some buffers
!         nquartet=0
!         do i=1,fcf%na
!             nquartet=nquartet+fcf%atom(i)%n
!         enddo
!         allocate(flat_m(81,nquartet))
!         allocate(flat_r2(3,nquartet))
!         allocate(flat_r3(3,nquartet))
!         allocate(flat_r4(3,nquartet))
!         allocate(flat_cwrk(81,fcf%na**4))
!         allocate(flat_ai(4,nquartet))
!         l=0
!         do i=1,fcf%na
!         do j=1,fcf%atom(i)%n
!             l=l+1
!             flat_r2(:,l)=fcf%atom(i)%quartet(j)%lv2
!             flat_r3(:,l)=fcf%atom(i)%quartet(j)%lv3
!             flat_r4(:,l)=fcf%atom(i)%quartet(j)%lv4
!             flat_m(:,l)=lo_flattentensor( fcf%atom(i)%quartet(j)%mwm )
!             flat_ai(:,l)=[fcf%atom(i)%quartet(j)%i1,fcf%atom(i)%quartet(j)%i2,fcf%atom(i)%quartet(j)%i3,fcf%atom(i)%quartet(j)%i4]
!         enddo
!         enddo
!         lo_allocate(egv(dr%n_mode,4))
!         egv=0.0_r8
!         ! the unit
!         prefactor=1.0_r8/32.0_r8
!
!         ! start integrating
!         ctrlong=0
!         do q1=1,qp%n_irr_point
!         do q2=1,qp%n_full_point
!             ! q-vectors
!             qv1= qp%ip(q1)%r*lo_twopi
!             qv2=-qp%ip(q1)%r*lo_twopi
!             qv3= qp%ap(q2)%r*lo_twopi
!             qv4=-qp%ap(q2)%r*lo_twopi
!             ! store omega and planck-factors in a buffer
!             do b1=1,dr%n_mode
!                 dumom1(b1)=dr%iq(q1)%omega(b1)
!                 dumom3(b1)=dr%aq(q2)%omega(b1)
!                 dump1(b1)=lo_planck(temperature,dumom1(b1))
!                 dump3(b1)=lo_planck(temperature,dumom3(b1))
!             enddo
!             ! technically four band-loops, but b1=b2 and b3=b4
!             do b1=1,dr%n_mode
!             ctrlong=ctrlong+1; if ( mod(ctrlong,int(mw%n,r8)) .ne. mw%r ) cycle
!             do b3=1,dr%n_mode
!                 omega=[dumom1(b1),dumom1(b1),dumom3(b3),dumom3(b3)]
!                 if ( minval(omega) .lt. omthres ) cycle
!                 egv(:,1)=dr%iq(q1)%egv(:,b1)
!                 egv(:,2)=conjg(dr%iq(q1)%egv(:,b1))
!                 egv(:,3)=dr%aq(q2)%egv(:,b3)
!                 egv(:,4)=conjg(dr%aq(q2)%egv(:,b3))
!                 call loc_4phpsi(fcf%na,nquartet,flat_r2,flat_r3,flat_r4,flat_m,flat_cwrk,flat_ai,omega,egv,qv2,qv3,qv4,c0)
!                 psi=real(c0)
!                 n1=dump1(b1)
!                 n2=dump3(b3)
!                 f1=(2*n1+1)*(2*n2+1)*psi*prefactor
!                 en4(b1,q1)=en4(b1,q1)+f1
!             enddo
!             enddo
!         enddo
!         enddo
!
!         if ( mw%talk ) then
!             t1=walltime()
!             write(*,*) 'fourth:',t1-t0
!             t0=t1
!         endif
!
!     end block fourth
!
!     !> add things up
!     finalize: block
!         integer :: i,j
!         call mpi_allreduce(MPI_IN_PLACE,en3,dr%n_mode*qp%n_irr_point,MPI_DOUBLE_PRECISION,MPI_SUM,mw%comm,mw%error )
!         call mpi_allreduce(MPI_IN_PLACE,en4,dr%n_mode*qp%n_irr_point,MPI_DOUBLE_PRECISION,MPI_SUM,mw%comm,mw%error )
!         en3=en3/real(qp%n_full_point,r8)/real(fc%na,r8)
!         en4=en4/real(qp%n_full_point,r8)/real(fc%na,r8)
!         free_energy_thirdorder=0.0_r8
!         free_energy_fourthorder=0.0_r8
!         do i=1,qp%n_irr_point
!             free_energy_thirdorder=free_energy_thirdorder+sum(en3(:,i))*qp%ip(i)%integration_weight
!             free_energy_fourthorder=free_energy_fourthorder+sum(en4(:,i))*qp%ip(i)%integration_weight
!         enddo
!     end block finalize
!
!     contains
!     ! specialized four-phonon matrix element, where q2=-q1, q4=-q3
!     subroutine loc_4phpsi(na,nquartet,flv2,flv3,flv4,fm,fcw,fai,omega,egv,q2,q3,q4,psi)
!         !> temporary work arrays, that actually is the force constant
!         integer, intent(in) :: na,nquartet
!         real(r8), dimension(3,nquartet), intent(in) :: flv2,flv3,flv4
!         real(r8), dimension(81,nquartet), intent(in) :: fm
!         complex(r8), dimension(81,na*na*na*na), intent(inout) :: fcw
!         integer, dimension(4,nquartet), intent(in) :: fai
!         !> four frequencies
!         real(r8), dimension(4), intent(in) :: omega
!         !> four eigenvectors
!         complex(r8), dimension(na*3,4), intent(in) :: egv
!         !> three q-points
!         real(r8), dimension(3), intent(in) :: q2,q3,q4
!         !> the matrix element
!         complex(r8), intent(out) :: psi
!
!         complex(r8) :: expiqr,c0
!         real(r8), dimension(3) :: rv2,rv3,rv4
!         real(r8) :: iqr,omegaprod
!         real(r8), parameter :: enhet=1.0_r8 !5.8104841317859603E074_r8
!         integer :: atom1,atom2,atom3,atom4,i,j,k,l,k1,k2,k3,k4,ii,jj
!         integer :: na3,na2
!
!         na3=na*na*na
!         na2=na*na
!         ! Eigenvector product, again
!         do atom1=1,na
!         do atom2=1,na
!         do atom3=1,na
!         do atom4=1,na
!             ii=(atom1-1)*na3+(atom2-1)*na2+(atom3-1)*na+atom4
!             do i=1,3
!             do j=1,3
!             do k=1,3
!             do l=1,3
!                 jj=(i-1)*27+(j-1)*9+(k-1)*3+l
!                 k1=(atom1-1)*3+i
!                 k2=(atom2-1)*3+j
!                 k3=(atom3-1)*3+k
!                 k4=(atom4-1)*3+l
!                 fcw(jj,ii)=egv(k1,1)*egv(k2,2)*egv(k3,3)*egv(k4,4)
!             enddo
!             enddo
!             enddo
!             enddo
!         enddo
!         enddo
!         enddo
!         enddo
!
!         ! Frequency product
!         omegaprod=omega(1)*omega(2)*omega(3)*omega(4)
!         omegaprod=sqrt(omegaprod)
!         ! The actual matrix element
!         c0=0.0_r8
!         do i=1,nquartet
!             atom1=fai(1,i)
!             atom2=fai(2,i)
!             atom3=fai(3,i)
!             atom4=fai(4,i)
!             ii=(atom1-1)*na3+(atom2-1)*na2+(atom3-1)*na+atom4
!             rv2=flv2(:,i)
!             rv3=flv3(:,i)
!             rv4=flv4(:,i)
!             iqr=dot_product(q2,rv2)+dot_product(q3,rv3)+dot_product(q4,rv4)
!             expiqr=cmplx(cos(iqr),sin(iqr),r8)
!             c0=c0+sum(fcw(:,ii)*fm(:,i))*expiqr
!         enddo
!         psi=enhet*c0/omegaprod
!     end subroutine
!     !> The three-phonon matrix element
!     pure subroutine loc_3phpsi(fc,omega,egv1,egv2,egv3,q2,q3,egvprod,psi)
!         !> the third order force constant
!         class(lo_forceconstant_thirdorder), intent(in) :: fc
!         !> the frequencies in question
!         real(r8), dimension(3), intent(in) :: omega
!         !> the eigenvectors
!         complex(r8), dimension(fc%na*3), intent(in) :: egv1,egv2,egv3
!         !> the two q-vectors that matter
!         real(r8), dimension(3), intent(in) :: q2,q3
!         !> work array
!         complex(r8), dimension(3,3,3,fc%na,fc%na,fc%na), intent(inout) :: egvprod
!         ! Matrix element
!         complex(r8), intent(out) :: psi
!
!         complex(r8) :: expiqr,c0
!         real(r8), dimension(3) :: rv2,rv3
!         real(r8) :: iqr,omegaprod
!         real(r8), parameter :: enhet=1.0_r8 !2.367755374064161E51_r8
!         integer :: atom1,atom2,atom3,i,j,k,k1,k2,k3,t
!
!         ! Eigenvector product
!         do atom3=1,fc%na
!         do atom2=1,fc%na
!         do atom1=1,fc%na
!             do k=1,3
!             do j=1,3
!             do i=1,3
!                 k1=(atom1-1)*3+i
!                 k2=(atom2-1)*3+j
!                 k3=(atom3-1)*3+k
!                 egvprod(i,j,k,atom1,atom2,atom3)=egv1(k1)*egv2(k2)*egv3(k3)
!             enddo
!             enddo
!             enddo
!         enddo
!         enddo
!         enddo
!         ! Frequency product
!         omegaprod=omega(1)*omega(2)*omega(3)
!         omegaprod=sqrt(omegaprod)
!         ! The actual matrix element
!         c0=0.0_r8
!         do atom1=1,fc%na
!         do t=1,fc%atom(atom1)%n
!             atom2=fc%atom(atom1)%triplet(t)%i2
!             atom3=fc%atom(atom1)%triplet(t)%i3
!
!             rv2=fc%atom(atom1)%triplet(t)%lv2
!             rv3=fc%atom(atom1)%triplet(t)%lv3
!
!             iqr=dot_product(q2,rv2)+dot_product(q3,rv3)
!             expiqr=cmplx(cos(iqr),sin(iqr),r8)
!             ! Note that the masses are premultiplied into the forceconstants, much faster.
!             c0=c0+sum(fc%atom(atom1)%triplet(t)%mwm*egvprod(:,:,:,atom1,atom2,atom3))*expiqr
!         enddo
!         enddo
!         psi=enhet*c0/omegaprod
!     end subroutine
!     !> The three-phonon matrix element, but product of two complex thingies
!     subroutine loc_3phdbl(fc,omega,egv1,egv2,egv3,q1,q2,egvprod1,egvprod2,doublepsi)
!         !> the third order force constant
!         class(lo_forceconstant_thirdorder), intent(in) :: fc
!         !> the frequencies in question
!         real(r8), dimension(3), intent(in) :: omega
!         !> the eigenvectors
!         complex(r8), dimension(fc%na*3), intent(in) :: egv1,egv2,egv3
!         !> the two q-vectors that matter
!         real(r8), dimension(3), intent(in) :: q1,q2
!         !> work arrays
!         complex(r8), dimension(3,3,3,fc%na,fc%na,fc%na), intent(inout) :: egvprod1,egvprod2
!         !> double psi
!         complex(r8), intent(out) :: doublepsi
!         !
!         integer :: atom1,atom2,atom3,i,j,k,k1,k2,k3,t
!         complex(r8) :: expiqr,c0,c1
!         real(r8), dimension(3) :: rv2,rv3
!         real(r8), parameter :: enhet=2.367755374064161E51_r8
!         real(r8) :: iqr,omegaprod1,omegaprod2
!
!         ! Not enough room on the stack for large unitcells, have to allocate space
!         do atom3=1,fc%na
!         do atom2=1,fc%na
!         do atom1=1,fc%na
!             do k=1,3
!             do j=1,3
!             do i=1,3
!                 k1=(atom1-1)*3+i
!                 k2=(atom2-1)*3+j
!                 k3=(atom3-1)*3+k
!                 egvprod1(i,j,k,atom1,atom2,atom3)=egv1(k1)*conjg(egv1(k2))*egv3(k3)
!                 egvprod2(i,j,k,atom1,atom2,atom3)=egv2(k1)*conjg(egv2(k2))*conjg(egv3(k3))
!             enddo
!             enddo
!             enddo
!         enddo
!         enddo
!         enddo
!
!         ! Frequency product
!         omegaprod1=omega(1)*omega(1)*omega(3)
!         omegaprod2=omega(2)*omega(2)*omega(3)
!         omegaprod1=sqrt(omegaprod1)
!         omegaprod2=sqrt(omegaprod2)
!         ! The actual matrix element
!         c0=0.0_r8
!         c1=0.0_r8
!         do atom1=1,fc%na
!         do t=1,fc%atom(atom1)%n
!             atom2=fc%atom(atom1)%triplet(t)%i2
!             atom3=fc%atom(atom1)%triplet(t)%i3
!             rv2=fc%atom(atom1)%triplet(t)%lv2
!             rv3=fc%atom(atom1)%triplet(t)%lv3
!             !
!             iqr=dot_product(-q1,rv2)
!             expiqr=cmplx(cos(iqr),sin(iqr),r8)
!             c0=c0+sum(fc%atom(atom1)%triplet(t)%mwm*egvprod1(:,:,:,atom1,atom2,atom3))*expiqr
!             !
!             iqr=dot_product(-q2,rv2)
!             expiqr=cmplx(cos(iqr),sin(iqr),r8)
!             c1=c1+sum(fc%atom(atom1)%triplet(t)%mwm*egvprod2(:,:,:,atom1,atom2,atom3))*expiqr
!         enddo
!         enddo
!         ! I predid the multiplication
!         c0=enhet*c0/omegaprod1
!         c1=enhet*c1/omegaprod2
!         doublepsi=c0*c1
!     end subroutine
! end subroutine

end module
