
!!> calculate things as a density of states
!subroutine get_free_energy_on_grid(qpd,drd,uc,fc,fct,fcf,qp,dr,loto,opts,mw,lsmpi)
!    !> q-point mesh
!    type(lo_monkhorst_pack_mesh), intent(in) :: qpd
!    !> harmonic properties on this mesh
!    type(lo_phonon_dispersions), intent(in) :: drd
!    !> crystal structure
!    type(lo_crystalstructure), intent(inout) :: uc
!    !> second order force constant
!    type(lo_forceconstant_secondorder), intent(in) :: fc
!    !> third order force constant
!    type(lo_forceconstant_thirdorder), intent(in) :: fct
!    !> fourth order force constant
!    type(lo_forceconstant_fourthorder), intent(in) :: fcf
!    !> q-point mesh
!    type(lo_monkhorst_pack_mesh), intent(in) :: qp
!    !> harmonic properties on this mesh
!    type(lo_phonon_dispersions), intent(in) :: dr
!    !> electrostatic corrections
!    type(lo_loto), intent(in) :: loto
!    !> all settings
!    type(lo_opts), intent(in) :: opts
!    !> mpi communicator
!    type(lo_mpiinfo), intent(in) :: mw
!    !> mpi helper
!    type(lo_lsmpi), intent(in) :: lsmpi
!
!    ! stuff for the DOS mesh
!    type(lo_phonon_selfenergy) :: se
!    type(lo_phonon_dos) :: pd
!    complex(flyt), dimension(3) :: cv0
!    real(flyt), dimension(:,:,:), allocatable :: dumdos,dumdosbuf
!    real(flyt), dimension(:,:), allocatable :: dum_site,dum_mode
!    real(flyt), dimension(:), allocatable :: dum
!    real(flyt) :: f0,f1,t0,sigma,weight
!    integer :: q1,lqp,mode,i
!
!    t0=mpi_wtime()
!    ! switch off openMP
!    call omp_set_num_threads(1)
!
!    if ( mw%talk ) call lo_progressbar_init()
!    lo_allocate(dumdos(opts%nf,dr%nb,qpd%nq_irr))
!    lo_allocate(dumdosbuf(opts%nf,dr%nb,qpd%nq_irr))
!    dumdos=0.0_flyt
!    dumdosbuf=0.0_flyt
!    do q1=1,lsmpi%nq
!        lqp=lsmpi%ind(q1)
!        ! get the actual self-energy
!        call se%generate(qpd%ip(lqp),drd%iq(lqp),uc,fc,fct,fcf,qp,dr,loto,opts)
!        ! and the intensity
!        do mode=1,drd%nb
!            call getintensity(se%faxis,se%im_3ph(:,mode)+se%im_iso(:,mode),&
!                 se%re_3ph(:,mode)+se%re_4ph(:,mode),se%p%omega(mode),se%faxis,&
!                 se%intensity(:,mode))
!            f0=lo_trapezoid_integration(se%faxis,se%intensity(:,mode))
!            ! if it's tiny, add a really narry gaussian instead.
!            if ( f0 .lt. lo_tol ) then
!                if ( drd%iq(q1)%omega(mode) .gt. drd%omega_min*0.5_flyt ) then
!                    do i=1,se%nf
!                        se%intensity(i,mode)=lo_gauss(se%faxis(i),drd%iq(lqp)%omega(mode),abs(se%faxis(3)-se%faxis(1)))
!                    enddo
!                endif
!            endif
!            ! store the intensities
!            dumdosbuf(:,mode,lqp)=dumdosbuf(:,mode,lqp)+se%intensity(:,mode)
!        enddo
!        !
!        if ( mw%talk ) call lo_progressbar(' ... calculating phonon dos',q1,lsmpi%nq,mpi_wtime()-t0)
!        !
!    enddo
!
!    ! Sum it up
!    call mpi_allreduce(dumdosbuf,dumdos,qpd%nq_irr*opts%nf*dr%nb,MPI_DOUBLE_PRECISION,MPI_SUM,mw%comm,mw%error)
!
!    ! Build the phonon dos structure
!    lo_allocate(pd%omega(opts%nf))
!    lo_allocate(pd%pdos_site(opts%nf,uc%na))
!    lo_allocate(pd%pdos_mode(opts%nf,dr%nb))
!    lo_allocate(pd%dos(opts%nf))
!    pd%na=uc%na
!    pd%nb=dr%nb
!    pd%ndos=opts%nf
!    pd%enhet=opts%enhet
!    pd%omega=se%intensityaxis
!    pd%pdos_site=0.0_flyt
!    pd%pdos_mode=0.0_flyt
!    pd%dos=0.0_flyt
!    pd%dosprojection=.true.
!    pd%dossmear=drd%default_smearing()
!    pd%dosmax=maxval(pd%omega) 
!
!    t0=mpi_wtime()
!    ! Add it up, with a little smearing
!    lo_allocate(dum_site(pd%ndos,pd%na))
!    lo_allocate(dum_mode(pd%ndos,pd%nb))
!    lo_allocate(dum(pd%ndos))
!    dum_site=0.0_flyt
!    dum_mode=0.0_flyt
!    if ( mw%talk ) call lo_progressbar_init()
!    do q1=1,lsmpi%nq
!        lqp=lsmpi%ind(q1)
!        weight=qpd%ip(lqp)%weight
!        !
!        do mode=1,drd%nb
!            ! skip acoustic at gamma
!            if ( drd%iq(q1)%omega(mode) .lt. drd%omega_min*0.5_flyt ) cycle
!            ! ok, smear a little
!            sigma=qpd%smearingparameter(drd%iq(lqp)%vel(:,mode),pd%dossmear,opts%sigma)
!            call lo_put_function_on_new_axis(se%faxis,dumdos(:,mode,lqp),pd%omega,dum,sigma)
!            ! make sure it's normalized after convolution
!            dum=dum/lo_trapezoid_integration(pd%omega,dum)
!            dum=dum*weight
!            ! add it to the mode-projected
!            dum_mode(:,mode)=dum_mode(:,mode)+dum
!            ! get the site-projection thing
!            do i=1,uc%na
!                cv0=drd%iq(lqp)%egv((i-1)*3+1:i*3,mode)
!                f0=abs(dot_product(cv0,conjg(cv0)))
!                dum_site(:,i)=dum_site(:,i)+f0*dum
!            enddo
!        enddo
!        !
!        if ( mw%talk ) call lo_progressbar(' ... integrating phonon dos',q1,lsmpi%nq,mpi_wtime()-t0)
!        !
!    enddo    
!
!    ! Add things together
!    call mpi_allreduce(dum_site,pd%pdos_site,opts%nf*uc%na,MPI_DOUBLE_PRECISION,MPI_SUM,mw%comm,mw%error)
!    call mpi_allreduce(dum_mode,pd%pdos_mode,opts%nf*dr%nb,MPI_DOUBLE_PRECISION,MPI_SUM,mw%comm,mw%error)
!    
!    if ( mw%r .eq. 0 ) then
!        ! Start normalizing, first the whole thing
!        do mode=1,drd%nb
!            pd%dos=pd%dos+pd%pdos_mode(:,mode)
!            f0=lo_trapezoid_integration(pd%omega,pd%pdos_mode(:,mode))
!            pd%pdos_mode=pd%pdos_mode/f0
!        enddo
!        ! remove stupid stuff
!        f1=pd%dos(1)
!        do i=1,pd%ndos
!            f0=(pd%ndos-i)*1.0_flyt/( (pd%ndos-1)*1.0_flyt )
!            f0=f0**2
!            pd%dos(i)=pd%dos(i)-f0*f1
!            if ( pd%dos(i) .lt. 0.0_flyt ) pd%dos(i)=0.0_flyt
!        enddo
!        pd%dos(1)=0.0_flyt
!        !
!        f0=lo_trapezoid_integration(pd%omega,pd%dos)
!        pd%dos=pd%dos*dr%nb/f0
!        ! and the site-projected
!        do i=1,uc%na
!            f0=lo_trapezoid_integration(pd%omega,pd%pdos_site(:,i))
!            pd%pdos_site(:,i)=pd%pdos_site(:,i)/f0
!        enddo
!        ! make sure the pojected adds up to the total
!        do i=1,pd%ndos
!            if ( pd%dos(i) .gt. lo_tol/lo_twopi/1E12_flyt ) then
!                f0=sum(pd%pdos_mode(i,:))
!                pd%pdos_mode(i,:)=pd%pdos_mode(i,:)*pd%dos(i)/f0
!                f0=sum(pd%pdos_site(i,:))
!                pd%pdos_site(i,:)=pd%pdos_site(i,:)*pd%dos(i)/f0
!            endif
!        enddo
!        ! Write it to file
!        call pd%write_to_file(uc,opts%enhet,pdfplot=.false.,modeproj=.false.,&
!                              filename='outfile.phonon_dos_lineshape')
!        !
!    endif
!
!end subroutine


