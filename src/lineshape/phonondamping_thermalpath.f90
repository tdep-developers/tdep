
!> calculate the intensity along a path
subroutine get_thermally_broadened_intensity_along_path(bs, uc, fc, fct, fcf, qp, dr, loto, opts, mw, lsmpi)
    !> the bandstructure
    type(lo_phonon_bandstructure), intent(inout) :: bs
    !> crystal structure
    type(lo_crystalstructure), intent(inout) :: uc
    !> second order force constant
    type(lo_forceconstant_secondorder), intent(in) :: fc
    !> third order force constant
    type(lo_forceconstant_thirdorder), intent(in) :: fct
    !> fourth order force constant
    type(lo_forceconstant_fourthorder), intent(in) :: fcf
    !> q-point mesh
    type(lo_monkhorst_pack_mesh), intent(in) :: qp
    !> harmonic properties on this mesh
    type(lo_phonon_dispersions), intent(in) :: dr
    !> electrostatic corrections
    type(lo_loto), intent(in) :: loto
    !> all settings
    type(lo_opts), intent(in) :: opts
    !> mpi communicator
    type(lo_mpiinfo), intent(in) :: mw
    !> mpi helper
    type(lo_lsmpi), intent(in) :: lsmpi

    type(lo_phonon_selfenergy) :: se
    !real(flyt), dimension(:,:,:), allocatable :: rebuf,imbuf,dumbuf
    !real(flyt), dimension(:,:), allocatable :: intbuf,thintbuf,lwbuf,shbuf,kernel
    !real(flyt), dimension(:), allocatable :: dum
    !real(flyt), dimension(3) :: qv1,qv2
    !real(flyt), dimension(2) :: lsintx,lsinty
    !real(flyt) :: t0,f0,f1,f2
    !complex(flyt) :: c0,c1,c2
    !integer :: q1,nb,lqp,i,j,path,ii,jj,k,nkern,iii,jjj
    type(lo_phonon_dispersions_qpoint) :: ompoint
    type(lo_qpoint) :: qpoint
    real(flyt), dimension(:, :), allocatable :: q_in_sphere, intbuf1, intbuf2
    real(flyt), dimension(:, :), allocatable :: spf1
    real(flyt), dimension(3) :: v0
    real(flyt) :: sigma_q, sigma_e, f0, f1, t0
    integer :: q1, q2, b1, i, j, k, l, lqp
    integer :: nspherepoints
    logical :: onlyharmonic

    call lo_seed_random_numbers()
    sigma_q = 0.015_flyt
    sigma_e = 0.1*lo_twopi*1E12_flyt !dr%default_smearing()
    onlyharmonic = .true.
    t0 = walltime()
    ! Get some points on a sphere
    nspherepoints = 1000 !400
    getpointsonsphere: block
        real(flyt), dimension(:, :), allocatable :: dum
        real(flyt), dimension(3) :: v0
        integer :: ctr
        !
        lo_allocate(q_in_sphere(3, nspherepoints))
        lo_allocate(dum(3, nspherepoints))
        q_in_sphere = 0.0_flyt
        dum = 0.0_flyt
        if (mw%r .eq. 0) then
            q_in_sphere = 0.0_flyt
            ctr = 0
            do
                call random_number(v0)
                v0 = 2.0_flyt*(v0 - 0.5_flyt)
                if (lo_sqnorm(v0) .gt. 1.0_flyt) cycle
                ctr = ctr + 1
                dum(:, ctr) = v0
                if (ctr .eq. nspherepoints) exit
            end do
            ! make the sphere 3 sigmas large
            dum = dum*3*sigma_q
            dum(:, 1) = 0.0_flyt

!do i=1,size(dum,2)
!    write(*,*) i,dum(:,i)
!enddo

        end if
        ! make sure all ranks have this
        call mpi_allreduce(dum, q_in_sphere, 3*nspherepoints, MPI_DOUBLE_PRECISION, MPI_SUM, mw%comm, mw%error)
    end block getpointsonsphere

    if (onlyharmonic) then
        se%n_energy = opts%nf
        lo_allocate(se%energy_axis_spectralfunction(se%n_energy))
        call lo_linspace(0.0_flyt, opts%maxf*dr%omega_max, se%energy_axis_spectralfunction)
    end if

    ! Make some space for stuff
    lo_allocate(spf1(se%n_energy, dr%nb))
    lo_allocate(intbuf1(bs%nptot, se%n_energy))
    lo_allocate(intbuf2(bs%nptot, se%n_energy))
    intbuf1 = 0.0_flyt
    intbuf2 = 0.0_flyt

    ! Do the entire brutal calculation:
    if (mw%talk) call lo_progressbar_init()
    do lqp = 1, lsmpi%nq
        q1 = lsmpi%ind(lqp)
        do q2 = 1, nspherepoints

            ! Get the harmonic stuff
            v0 = bs%q(q1)%v + q_in_sphere(:, q2)
            call harmonic_things_at_single_q(v0, uc, fc, loto, qpoint, ompoint)

            ! Get the raw lineshapes
            if (onlyharmonic) then
                ! Get a fake lineshape, just from the experimental energyresolution or something
                do b1 = 1, dr%nb
                    do i = 1, se%n_energy
                        spf1(i, b1) = lo_gauss(ompoint%omega(b1), se%energy_axis_spectralfunction(i), sigma_e)
                    end do
                end do
            else
                ! Get the actual lineshape stuff
                write (*, *) 'Actual lineshape not done yet'
                stop
            end if

            ! Add it together
            f0 = lo_gauss(0.0_flyt, norm2(q_in_sphere(:, q2)), sigma_q)
            do b1 = 1, dr%nb
                intbuf1(q1, :) = intbuf1(q1, :) + f0*spf1(:, b1)
                do i = 1, se%n_energy
                    f1 = se%energy_axis_spectralfunction(i) !max(dr%omega_max/50.0_flyt,se%energy_axis_spectralfunction(i))
                    spf1(i, b1) = spf1(i, b1)*(1 + lo_planck(opts%temperature, f1))
                end do
                intbuf2(q1, :) = intbuf2(q1, :) + f0*spf1(:, b1)*ompoint%thermal_prefactor(b1)
            end do

        end do
        if (mw%talk) then
            call lo_progressbar(' ... broadened spectral function', lqp, lsmpi%nq, walltime() - t0)
        end if
    end do

    lo_allocate(bs%intensity(bs%nptot, se%n_energy))
    lo_allocate(bs%intensity_with_prefactor(bs%nptot, se%n_energy))
    lo_allocate(bs%faxis(se%n_energy))
    bs%intensity = 0.0_flyt
    bs%intensity_with_prefactor = 0.0_flyt
    bs%faxis = se%energy_axis_spectralfunction

    call mpi_allreduce(intbuf1, bs%intensity, bs%nptot*opts%nf, MPI_DOUBLE_PRECISION, MPI_SUM, mw%comm, mw%error)
    call mpi_allreduce(intbuf2, bs%intensity_with_prefactor, bs%nptot*opts%nf, MPI_DOUBLE_PRECISION, MPI_SUM, mw%comm, mw%error)
    lo_deallocate(intbuf1)
    lo_deallocate(intbuf2)
    if (mw%r .eq. 0) then
        call bs%write_intensity(opts%enhet, logscale=.true.)
    end if

    !stop

!    t0=mpi_wtime()
!    ! Make space for linewidths
!    do q1=1,bs%nptot
!        allocate(bs%p(q1)%linewidth(bs%nb))
!        allocate(bs%p(q1)%shift(bs%nb))
!        allocate(bs%p(q1)%threephononphasespace(bs%nb))
!        allocate(bs%p(q1)%thermal_prefactor(bs%nb))
!    enddo
!    ! Space for intensity
!    allocate(bs%intensity(bs%nptot,opts%nf))
!    allocate(bs%intensity_with_prefactor(bs%nptot,opts%nf))
!    allocate(bs%selfenergy_real(bs%nptot,opts%nf,dr%nb))
!    allocate(bs%selfenergy_imag(bs%nptot,opts%nf,dr%nb))
!
!    ! Calculate the prefactors
!    do q1=1,bs%nptot
!    do i=1,bs%nb
!        f0=0.0_flyt
!        do j=1,uc%na
!            f1=dot_product(lo_twopi*bs%q(q1)%v,uc%rcart(:,j))
!            c0=dcmplx( cos(f1), sin(f1) )
!            ii=(j-1)*3+1
!            jj=j*3
!            c1=dot_product(lo_twopi*bs%q(q1)%v,bs%p(q1)%egv(ii:jj,i))
!            c2=c0*c1*uc%inelastic_neutron_cross_section(j)
!            f0=f0+abs(conjg(c2)*c2)
!        enddo
!        bs%p(q1)%thermal_prefactor(i)=f0
!    enddo
!    enddo
!
!    bs%intensity=0.0_flyt
!    bs%intensity_with_prefactor=0.0_flyt
!    bs%selfenergy_real=0.0_flyt
!    bs%selfenergy_imag=0.0_flyt
!
!    ! Dump some general info
!    if ( mw%talk ) then
!        write(*,*) '        isotope:',opts%isotopescattering
!        write(*,*) '    threephonon:',opts%thirdorder
!        write(*,*) '     fourphonon:',opts%fourthorder
!        write(*,*) '           loto:',opts%loto
!        write(*,*) 'integrationtype:',opts%integrationtype
!    endif
!
!    ! Turn off openmp
!    call omp_set_num_threads(1)
!
!    ! Calculate self energy
!    if ( mw%talk ) call lo_progressbar_init()
!    allocate(rebuf(bs%nptot,opts%nf,dr%nb))
!    allocate(imbuf(bs%nptot,opts%nf,dr%nb))
!    allocate(lwbuf(dr%nb,bs%nptot))
!    allocate(shbuf(dr%nb,bs%nptot))
!    rebuf=0.0_flyt
!    imbuf=0.0_flyt
!    lwbuf=0.0_flyt
!    shbuf=0.0_flyt
!    do q1=1,lsmpi%nq
!        ! global q-point index
!        lqp=lsmpi%ind(q1)
!        ! get the actual self-energy
!        call se%generate(bs%q(lqp),bs%p(lqp),uc,fc,fct,fcf,qp,dr,loto,opts)
!        ! Add it to the intensity
!        do j=1,bs%nb
!            do i=2,se%n_energy
!                imbuf(lqp,i,j)=se%im_3ph(i,j)+se%im_iso(i,j)
!                rebuf(lqp,i,j)=se%re_3ph(i,j)+se%re_4ph(i,j)
!            enddo
!            ! get the linewidth exactly at the harmonic frequency
!            lwbuf(j,lqp)=lo_linear_interpolation(se%energy_axis_selfenergy,se%im_3ph(:,j)+se%im_iso(:,j),bs%p(lqp)%omega(j))*2.0_flyt
!            ! get the anharmonic shift at the harmonic frequency
!            shbuf(j,lqp)=lo_linear_interpolation(se%energy_axis_selfenergy,se%re_3ph(:,j)+se%re_4ph(:,j),bs%p(lqp)%omega(j))
!        enddo
!        if ( mw%talk ) call lo_progressbar(' ... lineshape on path',q1,lsmpi%nq,mpi_wtime()-t0)
!    enddo
!
!    ! Add these up!
!    call mpi_allreduce(rebuf,bs%selfenergy_real,bs%nptot*opts%nf*bs%nb,MPI_DOUBLE_PRECISION,MPI_SUM,mw%comm,mw%error)
!    call mpi_allreduce(imbuf,bs%selfenergy_imag,bs%nptot*opts%nf*bs%nb,MPI_DOUBLE_PRECISION,MPI_SUM,mw%comm,mw%error)
!    allocate(intbuf(dr%nb,bs%nptot))
!    intbuf=0.0_flyt
!    call mpi_allreduce(lwbuf,intbuf,bs%nptot*bs%nb,MPI_DOUBLE_PRECISION,MPI_SUM,mw%comm,mw%error)
!    do i=1,bs%nptot
!    do j=1,bs%nb
!        bs%p(i)%linewidth(j)=intbuf(j,i)
!    enddo
!    enddo
!    intbuf=0.0_flyt
!    call mpi_allreduce(shbuf,intbuf,bs%nptot*bs%nb,MPI_DOUBLE_PRECISION,MPI_SUM,mw%comm,mw%error)
!    do i=1,bs%nptot
!    do j=1,bs%nb
!        bs%p(i)%shift(j)=intbuf(j,i)
!    enddo
!    enddo
!    deallocate(intbuf)
!    deallocate(lwbuf)
!    deallocate(shbuf)
!    ! Dump the shifts and widths
!    if ( mw%r .eq. 0 ) then
!        call bs%write_dispersive_property(opts%enhet,'shift','outfile.dispersion_shifts',.false.)
!        call bs%write_dispersive_property(opts%enhet,'linewidth','outfile.dispersion_linewidths',.false.)
!    endif
!    if ( mw%talk ) write(*,*) '... dumped some intermediate stuff'
!
!    rebuf=0.0_flyt
!    imbuf=0.0_flyt
!    t0=mpi_wtime()
!    ! Figure out some neat interpolation of self-energy for really small q
!    if ( mw%talk ) call lo_progressbar_init()
!    do q1=1,lsmpi%nq
!        lqp=lsmpi%ind(q1)
!        ! what path am I on?
!        path=bs%q(lqp)%path
!        ! The start and end-points
!        qv1=bs%segment(path)%r1-uc%bz%gshift( bs%segment(path)%r1 + lo_degenvector )
!        qv2=bs%segment(path)%r2-uc%bz%gshift( bs%segment(path)%r2 + lo_degenvector )
!        ! does it contain gamma?
!
!        if ( norm2(qv1) .gt. lo_tol .and. norm2(qv2) .gt. lo_tol ) cycle
!        ! seems it does, have to fix this, maybe.
!        ! Good small number to use
!        f0=(se%energy_axis_spectralfunction(2)-se%energy_axis_spectralfunction(1))*0.25_flyt ! smallest selfenergy
!        ! Is it in the beginning or the end?
!        if ( norm2(qv1) .lt. lo_tol ) then
!            ! Index of gamma
!            ii=(path-1)*bs%npts+1
!            ! Fix the acoustic branches
!            do j=1,3
!                ! Find index of point that is ok
!                do i=ii,ii+bs%npts-1
!                    if ( bs%p(i)%omega(j) .gt. dr%omega_min*0.5_flyt ) then
!                        jj=i
!                        exit
!                    endif
!                enddo
!                ! now I know that things are zero at ii, and ok at jj
!                bs%selfenergy_imag(ii,:,j)=f0       ! set imaginary at gamma
!                bs%selfenergy_real(ii,:,j)=0.0_flyt ! set real at gamma
!                lsintx(1)=bs%q_axis(ii)-lo_sqtol
!                lsintx(2)=bs%q_axis(jj)+lo_sqtol
!                ! Interpolate the missing self-energies at this point
!                do k=1,se%n_energy
!                    ! y-axis for interpolation, imaginary part
!                    lsinty(1)=bs%selfenergy_imag(ii,k,j)
!                    lsinty(2)=bs%selfenergy_imag(jj,k,j)
!                    imbuf(lqp,k,j)=lo_linear_interpolation(lsintx,lsinty,bs%q_axis(lqp))
!                    ! y-axis for interpolation, real part
!                    lsinty(1)=bs%selfenergy_real(ii,k,j)
!                    lsinty(2)=bs%selfenergy_real(jj,k,j)
!                    rebuf(lqp,k,j)=lo_linear_interpolation(lsintx,lsinty,bs%q_axis(lqp))
!                enddo
!            enddo
!        else
!            ! Same thing again, but this time gamma is at the end.
!            ii=path*bs%npts
!            ! loop over the three lowest branches
!            do j=1,3
!                jj=0
!                do i=ii,(path-1)*bs%npts+1,-1
!                    if ( bs%p(i)%omega(j) .gt. dr%omega_min*0.5_flyt ) then
!                        jj=i
!                        exit
!                    endif
!                enddo
!                ! Interpolate this, somehow
!                bs%selfenergy_imag(ii,:,j)=f0       ! set imaginary at gamma
!                bs%selfenergy_real(ii,:,j)=0.0_flyt ! set real at gamma
!                ! x-axis for interpolation
!                lsintx(2)=bs%q_axis(ii)-lo_sqtol
!                lsintx(1)=bs%q_axis(jj)+lo_sqtol
!                ! interpolate to missing points
!                do k=1,se%n_energy
!                    ! y-axis for interpolation
!                    lsinty(2)=bs%selfenergy_imag(ii,k,j)
!                    lsinty(1)=bs%selfenergy_imag(jj,k,j)
!                    imbuf(lqp,k,j)=lo_linear_interpolation(lsintx,lsinty,bs%q_axis(lqp))
!                    ! y-axis for interpolation
!                    lsinty(2)=bs%selfenergy_real(ii,k,j)
!                    lsinty(1)=bs%selfenergy_real(jj,k,j)
!                    rebuf(lqp,k,j)=lo_linear_interpolation(lsintx,lsinty,bs%q_axis(lqp))
!                enddo
!            enddo
!        endif
!        !
!        if ( mw%talk ) call lo_progressbar(' ... fixing tiny q',q1,lsmpi%nq,mpi_wtime()-t0)
!        !
!    enddo
!
!    ! Add this together, and add it to the self energy
!    allocate(dumbuf(bs%nptot,opts%nf,dr%nb))
!    dumbuf=0.0_flyt
!    call mpi_allreduce(rebuf,dumbuf,bs%nptot*opts%nf*bs%nb,MPI_DOUBLE_PRECISION,MPI_SUM,mw%comm,mw%error)
!    bs%selfenergy_real=bs%selfenergy_real+dumbuf
!    dumbuf=0.0_flyt
!    call mpi_allreduce(imbuf,dumbuf,bs%nptot*opts%nf*bs%nb,MPI_DOUBLE_PRECISION,MPI_SUM,mw%comm,mw%error)
!    bs%selfenergy_imag=bs%selfenergy_imag+dumbuf
!    deallocate(dumbuf)
!
!    ! Can't have anything negative in the imaginary selfenergy
!    do i=1,bs%nptot
!    do j=1,se%n_energy
!    do k=1,se%n_mode
!        bs%selfenergy_imag(i,j,k)=max(bs%selfenergy_imag(i,j,k),0.0_flyt)
!    enddo
!    enddo
!    enddo
!
!    if ( mw%talk ) then
!        write(*,*) '... selfenergies probably ok'
!        write(*,*) 'Im min,max,sum',minval(bs%selfenergy_imag/lo_twopi/1E12_flyt),&
!                                 maxval(bs%selfenergy_imag/lo_twopi/1E12_flyt),&
!                                 sum(bs%selfenergy_imag/lo_twopi/1E12_flyt)
!        write(*,*) 'Re min,max,sum',minval(bs%selfenergy_real/lo_twopi/1E12_flyt),&
!                                 maxval(bs%selfenergy_real/lo_twopi/1E12_flyt),&
!                                 sum(bs%selfenergy_real/lo_twopi/1E12_flyt)
!    endif

!    ! Probably a neat idea to ever so slightly smear the self-energies
!    rebuf=0.0_flyt
!    imbuf=0.0_flyt
!    ! A kernel to smear with
!    nkern=3
!    lo_allocate(kernel(2*nkern+1,2*nkern+1))
!    kernel=0.0_flyt
!    do i=-nkern,nkern
!    do j=-nkern,nkern
!        ii=i+nkern+1
!        jj=j+nkern+1
!        f0=norm2([i,j]*1.0_flyt)
!        kernel(ii,jj)=lo_gauss(0.0_flyt,f0,nkern*0.5_flyt)
!    enddo
!    enddo
!    f0=sum(bs%selfenergy_real)
!    f1=sum(bs%selfenergy_imag)
!    kernel=kernel/sum(kernel)
!    do q1=1,lsmpi%nq
!        lqp=lsmpi%ind(q1)
!        do j=1,bs%nb
!        do i=1,se%n_energy
!            do ii=1,2*nkern+1
!            do jj=1,2*nkern+1
!                iii=min(max(q1+ii-1-nkern,1),bs%nptot)
!                jjj=min(max(i+jj-1-nkern,1),se%n_energy)
!                rebuf(lqp,i,j)=rebuf(lqp,i,j)+kernel(ii,jj)*bs%selfenergy_real(iii,jjj,j)
!                imbuf(lqp,i,j)=imbuf(lqp,i,j)+kernel(ii,jj)*bs%selfenergy_imag(iii,jjj,j)
!            enddo
!            enddo
!        enddo
!        enddo
!    enddo
!    bs%selfenergy_real=0.0_flyt
!    bs%selfenergy_imag=0.0_flyt
!    call mpi_allreduce(rebuf,bs%selfenergy_real,bs%nptot*opts%nf*bs%nb,MPI_DOUBLE_PRECISION,MPI_SUM,mw%comm,mw%error)
!    call mpi_allreduce(imbuf,bs%selfenergy_imag,bs%nptot*opts%nf*bs%nb,MPI_DOUBLE_PRECISION,MPI_SUM,mw%comm,mw%error)
!    bs%selfenergy_real=bs%selfenergy_real*f0/sum(bs%selfenergy_real)
!    bs%selfenergy_imag=bs%selfenergy_imag*f1/sum(bs%selfenergy_imag)

!    ! Now all the self-energies are nice, time to get the intensities
!    allocate(intbuf(bs%nptot,opts%nf))
!    allocate(thintbuf(bs%nptot,opts%nf))
!    allocate(dum(opts%nf))
!    intbuf=0.0_flyt
!    thintbuf=0.0_flyt
!    t0=mpi_wtime()
!    ! Figure out some neat interpolation of self-energy for really small q
!    if ( mw%talk ) call lo_progressbar_init()
!    do q1=1,lsmpi%nq
!        lqp=lsmpi%ind(q1)
!        do j=1,bs%nb
!            ! Get the lineshape
!            dum=0.0_flyt
!            if ( bs%p(lqp)%omega(j) .gt. lo_freqtol ) then
!                call getintensity(se%energy_axis_selfenergy,bs%selfenergy_imag(lqp,:,j),bs%selfenergy_real(lqp,:,j),&
!                bs%p(lqp)%omega(j),se%energy_axis_spectralfunction,dum)
!            else
!                ! acoustic branch at Gamma. Add a gaussian at 0 to no make it disappear.
!                do i=1,se%n_energy
!                    dum(i)=lo_gauss(se%energy_axis_spectralfunction(i),0.0_flyt,se%energy_axis_spectralfunction(2)-se%energy_axis_spectralfunction(1))
!                enddo
!            endif
!            ! Add it to the intensity
!            intbuf(lqp,:)=intbuf(lqp,:)+dum
!            ! And the one with thermal factors
!            !do i=2,se%n_energy
!            !    f0=max(dr%omega_min*0.5_flyt,se%energy_axis_spectralfunction(i))
!            !    dum(i)=dum(i)*(lo_planck(opts%temperature,f0)+1.0_flyt)
!            !enddo
!            !thintbuf(lqp,:)=thintbuf(lqp,:)+dum*bs%p(lqp)%thermal_prefactor(j)
!        enddo
!        !
!        if ( mw%talk ) call lo_progressbar(' ... intensities',q1,lsmpi%nq,mpi_wtime()-t0)
!    enddo
!    ! add them up
!    call mpi_allreduce(intbuf,bs%intensity,bs%nptot*opts%nf,MPI_DOUBLE_PRECISION,MPI_SUM,mw%comm,mw%error)
!    call mpi_allreduce(thintbuf,bs%intensity_with_prefactor,bs%nptot*opts%nf,MPI_DOUBLE_PRECISION,MPI_SUM,mw%comm,mw%error)
!    deallocate(intbuf)
!    deallocate(thintbuf)
!    deallocate(dum)
!
!    ! Dump to file
!    if ( mw%r .eq. 0 ) then
!        write(*,*) 'Writing intensity to file'
!        lo_allocate(bs%faxis(se%n_energy))
!        bs%faxis=se%energy_axis_spectralfunction
!        call bs%write_intensity(opts%enhet,logscale=.true.)
!    endif
    ! And it is done!
end subroutine
