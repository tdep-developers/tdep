
!> Get the three-phonon self-energy via convolutions of the spectral function.
subroutine convolution_imaginary_selfenergy(se, wp, qp, dr, sr, ise, isf, p, temperature, mw, mem, verbosity)
    !> self-energy
    type(lo_phonon_selfenergy), intent(inout) :: se
    !> harmonic properties at this q-point
    type(lo_phonon_dispersions_qpoint), intent(in) :: wp
    !> qpoint mesh
    class(lo_qpoint_mesh), intent(in) :: qp
    !> harmonic dispersions
    type(lo_phonon_dispersions), intent(in) :: dr
    !> scattering rates
    type(lo_listofscatteringrates), intent(in) :: sr
    !> interpolated selfenergy
    type(lo_interpolated_selfenergy), intent(in) :: ise
    !> tabulated spectral functions
    type(lo_spectralfunction_helper), intent(in) :: isf
    !> structure
    type(lo_crystalstructure), intent(in) :: p
    !> temperature
    real(r8), intent(in) :: temperature
    !> MPI communicator
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem
    !> talk a lot?
    integer, intent(in) :: verbosity

    real(r8) :: timer, t0, t1

    timer = walltime()
    t0 = timer
    t1 = timer

    ! Might add some more pre-processing here.
    init: block
        if (verbosity .gt. 0) then
            write (lo_iou, *) '... convolution based integration'
        end if
        ! Reset values
        se%im_iso = 0.0_r8
        se%im_3ph = 0.0_r8
    end block init

    ! Try this again, this time with more feeling. The next attempt will
    ! make do with much less Fourier transforms, this one is to get it right.
    newattempt: block
        type(lo_convolution_helper) :: ch
        real(r8), dimension(:), allocatable :: buf_sfun, buf_sigma2, buf_sigma3
        real(r8), dimension(3) :: v0
        real(r8) :: pref, sigma2, sigma3, sigma, psisq
        integer :: imode, iq, jq, kq, mode1, mode2, mode3

        ! Set up the helper guy for the convolutions
        call ch%generate(se%energy_axis, temperature, dr%n_mode)

        ! Some temporary buffers
        allocate (buf_sigma2(dr%n_mode))
        allocate (buf_sigma3(dr%n_mode))
        buf_sigma2 = 0.0_r8
        buf_sigma3 = 0.0_r8
        allocate (buf_sfun(se%n_energy))
        buf_sfun = 0.0_r8

        if (sr%atgamma) then
            ! Use symmetry
            if (verbosity .gt. 0) call lo_progressbar_init()
            qploopirr: do iq = 1, qp%n_irr_point
                ! Make it MPI parallel
                if (mod(iq, mw%n) .ne. mw%r) cycle

                ! Pre-fetch and massage spectral functions before the mode loop
                select case (se%integrationtype)
                case (5)
                    ! This is the closed-grid or Gamma-point case, where everything has been prepared beforehand.
                    ! Get the smearing parameters per mode?
                    do imode = 1, dr%n_mode
                        buf_sigma2(imode) = qp%adaptive_sigma(qp%ip(iq)%radius, dr%iq(iq)%vel(:, imode), dr%default_smearing(imode), se%smearing_prefactor)
                        buf_sigma3(imode) = buf_sigma2(imode)
                    end do
                    ! Pre-fetch the spectral functions
                    call ch%buffer_spectral_functions(isf%spectralfunction(:, :, iq), isf%spectralfunction(:, :, iq))

                    ! do imode=1,dr%n_mode
                    !     do ie=1,n
                    !         buf_sf2_planck(ie,imode)=isf%spectralfunction(ie+1,imode,iq)
                    !         buf_sf2_planck(-ie,imode)=-buf_sf2_planck(ie,imode)
                    !     enddo
                    !     buf_sf2_planck(:,imode)=buf_sf2_planck(:,imode)*buf_planck_plus_one
                    !     ! Insert FFT here
                    ! enddo
                    !
                    ! ! Set them to sharp values instead, for debugging
                    ! buf_sf2_planck=0.0_r8
                    ! do imode=1,dr%n_mode
                    !     if ( dr%iq(iq)%omega(imode) .lt. lo_freqtol ) cycle
                    !     pref=lo_huge
                    !     j=-1
                    !     do ie=1,n
                    !         if ( abs(buf_x(ie)-dr%iq(iq)%omega(imode)) .lt. pref ) then
                    !             pref=abs(buf_x(ie)-dr%iq(iq)%omega(imode))
                    !             j=ie
                    !         endif
                    !     enddo
                    !     buf_sf2_planck(j,imode)=1.0_r8/deltax
                    !     buf_sf2_planck(-j,imode)=-1.0_r8/deltax
                    ! enddo
                    !
                    ! call ch%buffer_spectral_functions( buf_sf2_planck(0:n,:),buf_sf2_planck(0:n,:) )
                    !
                    ! do imode=1,dr%n_mode
                    !     buf_sf2_planck(:,imode)=buf_sf2_planck(:,imode)*buf_planck_plus_one
                    ! enddo
                case default
                    call lo_stop_gracefully(['Bad integration type'], lo_exitcode_param, __FILE__, __LINE__)
                end select

                ! Prefactor for this q-point
                pref = threephonon_prefactor*qp%ip(iq)%integration_weight

                ! Excellent. Now start convoluting every which way.
                mode2loop1: do mode2 = 1, dr%n_mode
                mode3loop1: do mode3 = mode2, dr%n_mode
                    ! Skip acoustic modes?
                    if (dr%iq(iq)%omega(mode2) .lt. lo_freqtol) cycle
                    if (dr%iq(iq)%omega(mode3) .lt. lo_freqtol) cycle
                    ! Smear?
                    sigma2 = qp%adaptive_sigma(qp%ip(iq)%radius, dr%iq(iq)%vel(:, mode2), dr%default_smearing(mode2), se%smearing_prefactor)
                    sigma3 = qp%adaptive_sigma(qp%ip(iq)%radius, dr%iq(iq)%vel(:, mode3), dr%default_smearing(mode3), se%smearing_prefactor)
                    sigma = sqrt(sigma2**2 + sigma3**2)
                    call ch%convolute_to_sfun(mode2, mode3, sigma, buf_sfun)

                    ! Excellent, everything was convoluted properly. Start accumulating self-energy
                    mode1loop1: do mode1 = 1, dr%n_mode
                        if (mode2 .eq. mode3) then
                            psisq = abs(sr%psi_3ph(mode1, mode2, mode3, iq)*conjg(sr%psi_3ph(mode1, mode2, mode3, iq)))
                            se%im_3ph(:, mode1) = se%im_3ph(:, mode1) + buf_sfun*pref*psisq
                        else
                            psisq = abs(sr%psi_3ph(mode1, mode2, mode3, iq)*conjg(sr%psi_3ph(mode1, mode2, mode3, iq))) + &
                                    abs(sr%psi_3ph(mode1, mode3, mode2, iq)*conjg(sr%psi_3ph(mode1, mode3, mode2, iq)))
                            se%im_3ph(:, mode1) = se%im_3ph(:, mode1) + buf_sfun*pref*psisq
                        end if
                    end do mode1loop1
                end do mode3loop1
                end do mode2loop1

                ! Since I have the spectral functions available I might as well sort
                ! out the isotope scattering here.
                pref = isotope_prefactor*qp%ip(iq)%integration_weight
                do mode2 = 1, dr%n_mode
                    buf_sfun = isf%spectralfunction(:, mode2, iq)
                    sigma2 = qp%adaptive_sigma(qp%ip(iq)%radius, dr%iq(iq)%vel(:, mode2), dr%default_smearing(mode2), se%smearing_prefactor)
                    call gaussian_smear_spectral_function(se%energy_axis, sigma2, buf_sfun)
                    do mode1 = 1, dr%n_mode
                        se%im_iso(:, mode1) = se%im_iso(:, mode1) + buf_sfun*sr%psi_iso(mode1, mode2, iq)*pref
                    end do
                end do

                if (verbosity .gt. 0 .and. iq .lt. qp%n_irr_point) then
                    t1 = walltime()
                    call lo_progressbar(' ... imaginary selfenergy', iq, qp%n_irr_point, t1 - t0)
                end if
            end do qploopirr
        else
            ! Not using symmetry, at a general q-point
            if (verbosity .gt. 0) call lo_progressbar_init()
            qploopfull: do iq = 1, qp%n_full_point
                ! Make it MPI parallel
                if (mod(iq, mw%n) .ne. mw%r) cycle

                ! Pre-fetch and massage spectral functions before the mode loop
                select case (se%integrationtype)
                case (5)
                    ! This is the closed-grid or Gamma-point case, where everything has been prepared beforehand.
                    ! First sort out the q-indices. This is q'
                    jq = qp%ap(iq)%irreducible_index
                    ! Sort out q''
                    v0 = matmul(p%inv_reciprocal_latticevectors, sr%qvec3(:, iq))
                    kq = index_on_grid(qp, v0)
                    kq = qp%ap(kq)%irreducible_index

                    ! Pre-fetch the spectral functions
                    call ch%buffer_spectral_functions(isf%spectralfunction(:, :, jq), isf%spectralfunction(:, :, kq))

                    ! ! Pre-fetch the spectral functions
                    ! do imode=1,dr%n_mode
                    !     do ie=1,n
                    !         buf_sf2_planck(ie,imode)=isf%spectralfunction(ie+1,imode,jq)
                    !         buf_sf2_planck(-ie,imode)=-buf_sf2_planck(ie,imode)
                    !         buf_sf3_planck(ie,imode)=isf%spectralfunction(ie+1,imode,kq)
                    !         buf_sf3_planck(-ie,imode)=-buf_sf3_planck(ie,imode)
                    !     enddo
                    !     buf_sf2_planck(:,imode)=buf_sf2_planck(:,imode)*buf_planck_plus_one
                    !     buf_sf3_planck(:,imode)=buf_sf3_planck(:,imode)*buf_planck_plus_one
                    !     ! Insert FFT here
                    ! enddo

                    ! ! Set them to sharp values instead, for debugging
                    ! buf_sf2_planck=0.0_r8
                    ! do imode=1,dr%n_mode
                    !     if ( dr%iq(iq)%omega(imode) .lt. lo_freqtol ) cycle
                    !     pref=lo_huge
                    !     j=-1
                    !     do ie=1,n
                    !         if ( abs(buf_x(ie)-dr%iq(iq)%omega(imode)) .lt. pref ) then
                    !             pref=abs(buf_x(ie)-dr%iq(iq)%omega(imode))
                    !             j=ie
                    !         endif
                    !     enddo
                    !     buf_sf2_planck(j,imode)=1.0_r8/deltax
                    !     buf_sf2_planck(-j,imode)=-1.0_r8/deltax
                    !     buf_sf2_planck(:,imode)=buf_sf2_planck(:,imode)*buf_planck_plus_one
                    ! enddo
                case default
                    call lo_stop_gracefully(['Bad integration type'], lo_exitcode_param, __FILE__, __LINE__)
                end select

                ! Prefactor for this q-point
                pref = threephonon_prefactor*qp%ap(iq)%integration_weight

                ! Excellent. Now start convoluting every which way.
                mode2loop2: do mode2 = 1, dr%n_mode
                mode3loop2: do mode3 = 1, dr%n_mode
                    ! Skip acoustic modes?
                    if (dr%iq(jq)%omega(mode2) .lt. lo_freqtol) cycle
                    if (dr%iq(kq)%omega(mode3) .lt. lo_freqtol) cycle
                    ! Smear?
                    sigma2 = qp%adaptive_sigma(qp%ip(jq)%radius, dr%iq(jq)%vel(:, mode2), dr%default_smearing(mode2), se%smearing_prefactor)
                    sigma3 = qp%adaptive_sigma(qp%ip(kq)%radius, dr%iq(kq)%vel(:, mode3), dr%default_smearing(mode3), se%smearing_prefactor)
                    sigma = sqrt(sigma2**2 + sigma3**2)

                    ! ! Do the convolution over modes
                    ! call lo_convolution(buf_sf2_planck(:,mode2),buf_sf2_planck(:,mode3),deltax,.false.,buf_conv)
                    ! ! Adjust with planck factors
                    ! buf_conv=buf_conv*buf_inv_planck_plus_one
                    ! buf_sfun=0.0_r8
                    ! do ie=1,n
                    !     buf_sfun(ie+1)=( buf_conv(ie)-buf_conv(-ie) )*0.5_r8
                    ! enddo
                    ! call gaussian_smear_spectral_function(se%energy_axis,sigma,buf_sfun)

                    call ch%convolute_to_sfun(mode2, mode3, sigma, buf_sfun)

                    ! Excellent, everything was convoluted properly. Start accumulating self-energy
                    mode1loop2: do mode1 = 1, dr%n_mode
                        psisq = abs(sr%psi_3ph(mode1, mode2, mode3, iq)*conjg(sr%psi_3ph(mode1, mode2, mode3, iq)))
                        se%im_3ph(:, mode1) = se%im_3ph(:, mode1) + buf_sfun*pref*psisq
                    end do mode1loop2
                end do mode3loop2
                end do mode2loop2

                ! Since I have the spectral functions available I might as well sort
                ! out the isotope scattering here.
                pref = isotope_prefactor*qp%ap(iq)%integration_weight
                do mode2 = 1, dr%n_mode
                    buf_sfun = isf%spectralfunction(:, mode2, jq)
                    sigma2 = qp%adaptive_sigma(qp%ip(jq)%radius, dr%iq(jq)%vel(:, mode2), dr%default_smearing(mode2), se%smearing_prefactor)
                    call gaussian_smear_spectral_function(se%energy_axis, sigma2, buf_sfun)
                    do mode1 = 1, dr%n_mode
                        se%im_iso(:, mode1) = se%im_iso(:, mode1) + buf_sfun*sr%psi_iso(mode1, mode2, iq)*pref
                    end do
                end do

                if (verbosity .gt. 0 .and. iq .lt. qp%n_full_point) then
                    t1 = walltime()
                    call lo_progressbar(' ... imaginary selfenergy', iq, qp%n_full_point, t1 - t0)
                end if
            end do qploopfull
        end if
        ! Cleanup
        call ch%destroy()

        ! Add together
        call mw%allreduce('sum', se%im_3ph)
        call mw%allreduce('sum', se%im_iso)

        if (verbosity .gt. 0) then
            t1 = walltime()
            call lo_progressbar(' ... imaginary selfenergy', qp%n_full_point, qp%n_full_point, t1 - t0)
            t0 = t1
        end if

        ! if ( mw%talk ) then
        !         call h5%init(__FILE__,__LINE__)
        !         call h5%open_file('write','sbuf.hdf5')
        !
        !         call h5%store_data(allsbuf,h5%file_id,'sfun')
        !
        !         call h5%close_file()
        !         call h5%destroy(__FILE__,__LINE__)
        ! endif
    end block newattempt

    ! actualcalculation: block
    !     real(r8), dimension(:,:), allocatable :: bufRe,bufIm,buf_spectral2,buf_spectral3
    !     real(r8), dimension(:), allocatable :: sfun,taperfn
    !     real(r8), dimension(:), allocatable :: buf_fpre,buf_om,buf_conv
    !     real(r8), dimension(:), allocatable :: buf_odd_sf2a,buf_odd_sf2b
    !     real(r8), dimension(:), allocatable :: buf_odd_sf3a,buf_odd_sf3b
    !     real(r8), dimension(3) :: v0
    !     real(r8) :: psisq,pref,minIm,sigma
    !     integer :: mode1,mode2,mode3,iq,jq,ctr
    !     integer :: ie,n
    !
    !     ! Oddly allocated buffers:
    !     n=se%n_energy
    !     allocate(buf_fpre(-n+1:n))
    !     allocate(buf_om(-n+1:n))
    !     allocate(buf_odd_sf2a(-n+1:n))
    !     allocate(buf_odd_sf2b(-n+1:n))
    !     allocate(buf_odd_sf3a(-n+1:n))
    !     allocate(buf_odd_sf3b(-n+1:n))
    !     allocate(buf_conv(-n+1:n))
    !     buf_fpre=0.0_r8
    !     buf_om=0.0_r8
    !     buf_odd_sf2a=0.0_r8
    !     buf_odd_sf2b=0.0_r8
    !     buf_odd_sf3a=0.0_r8
    !     buf_odd_sf3b=0.0_r8
    !     buf_conv=0.0_r8
    !
    !     ! get the planck prefactor
    !     buf_fpre=0.0_r8
    !     do ie=2,se%n_energy
    !         buf_fpre(ie)=lo_planck(temperature,se%energy_axis(ie))+0.5_r8
    !         buf_fpre(1-ie)=-buf_fpre(ie)
    !     enddo
    !
    !     ! Normally allocated buffers:
    !     call mem%allocate(sfun,se%n_energy,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    !     call mem%allocate(taperfn,se%n_energy,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    !     call mem%allocate(bufRe,[se%n_energy,dr%n_mode],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    !     call mem%allocate(bufIm,[se%n_energy,dr%n_mode],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    !     call mem%allocate(buf_spectral2,[se%n_energy,dr%n_mode],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    !     call mem%allocate(buf_spectral3,[se%n_energy,dr%n_mode],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    !     sfun=0.0_r8
    !     taperfn=0.0_r8
    !     bufRe=0.0_r8
    !     bufIm=0.0_r8
    !     buf_spectral2=0.0_r8
    !     buf_spectral3=0.0_r8
    !     minIm=( se%energy_axis(2)-se%energy_axis(1) )
    !     call taperfn_im(se%energy_axis,dr%omega_max,dr%omega_min,taperfn)
    !     ! Stop the init timer
    !
    !     if ( sr%atgamma ) then
    !         ! Use symmetry
    !         if ( verbosity .gt. 0 ) call lo_progressbar_init()
    !         qploopirr: do iq=1,qp%n_irr_point
    !             ! Make it MPI parallel
    !             if ( mod(iq,mw%n) .ne. mw%r ) cycle
    !
    !             ! First thing, I need the spectral functions at q'
    !
    !             ! We need the spectral functions at q', which is also the one at q''. Two different ways to go about it.
    !             select case(se%integrationtype)
    !             case(4)
    !                 ! Here I calculate the spectral functions from the interpolated self-energy, and then smear them appropriately.
    !
    !                 ! First thing, I need the spectral functions at q'
    !                 call ise%diagonal_selfenergy(p,qp%ip(iq)%r,dr%iq(iq)%omega,dr%iq(iq)%egv,bufRe,bufIm)
    !                 bufIm=max(bufIm,minIm)
    !                 do mode2=1,se%n_mode
    !                     if ( dr%aq(iq)%omega(mode2) .gt. lo_freqtol ) then
    !                         bufIm(:,mode2)=bufIm(:,mode2)*taperfn
    !                         call evaluate_spectral_function(ise%energy,bufIm(:,mode2),bufRe(:,mode2),dr%iq(iq)%omega(mode2),buf_spectral2(:,mode2))
    !                         sigma=qp%adaptive_sigma(qp%ip(iq)%radius,dr%iq(iq)%vel(:,mode2),dr%default_smearing(mode2),se%smearing_prefactor)
    !                         call gaussian_smear_spectral_function(ise%energy,sigma,buf_spectral2(:,mode2))
    !                         buf_spectral2(:,mode2)=buf_spectral2(:,mode2)/lo_trapezoid_integration(ise%energy,buf_spectral2(:,mode2))
    !                     else
    !                         buf_spectral2(:,mode2)=0.0_r8
    !                     endif
    !                 enddo
    !             case(5)
    !                 ! This is the closed-grid or Gamma-point case, where everything has been prepared beforehand.
    !                 ! Sort out q' first
    !                 buf_spectral2=isf%spectralfunction(:,:,iq)
    !             end select
    !
    !             ! Prefactor for this q-point
    !             pref=threephonon_prefactor*qp%ip(iq)%integration_weight
    !
    !             ! Excellent. Now start convoluting every which way.
    !             mode2loop1: do mode2=1,dr%n_mode
    !             mode3loop1: do mode3=mode2,dr%n_mode
    !                 ! Prepare for convolution
    !                 buf_odd_sf2a=0.0_r8
    !                 buf_odd_sf2b=0.0_r8
    !                 buf_odd_sf3a=0.0_r8
    !                 buf_odd_sf3b=0.0_r8
    !                 do ie=2,se%n_energy
    !                     buf_odd_sf2a(ie)=buf_spectral2(ie,mode2)
    !                     buf_odd_sf3a(ie)=buf_spectral2(ie,mode3)
    !                     buf_odd_sf2a(1-ie)=-buf_odd_sf2a(ie)
    !                     buf_odd_sf3a(1-ie)=-buf_odd_sf3a(ie)
    !                 enddo
    !                 buf_odd_sf2b=buf_odd_sf2a*buf_fpre
    !                 buf_odd_sf3b=buf_odd_sf3a*buf_fpre
    !
    !                 call lo_abcd_convolution(buf_odd_sf2a,buf_odd_sf3b,buf_odd_sf3a,buf_odd_sf2b,se%energy_axis(2)-se%energy_axis(1),buf_conv)
    !                 ! Post-convolution the values are FFT-shifted oddly so collect to the right place.
    !                 sfun=0.0_r8
    !                 do ie=1,se%n_energy
    !                     sfun(ie)=buf_conv(se%n_energy+1-ie)
    !                 enddo
    !
    !                 ! Excellent, everything was convoluted properly. Start accumulating self-energy
    !                 mode1loop1: do mode1=1,dr%n_mode
    !                     if ( mode2 .eq. mode3 ) then
    !                         psisq=abs( sr%psi_3ph(mode1,mode2,mode3,iq) * conjg( sr%psi_3ph(mode1,mode2,mode3,iq) ) )
    !                         se%im_3ph(:,mode1)=se%im_3ph(:,mode1)-sfun*pref*psisq
    !                     else
    !                         psisq=abs( sr%psi_3ph(mode1,mode2,mode3,iq) * conjg( sr%psi_3ph(mode1,mode2,mode3,iq) ) )+&
    !                               abs( sr%psi_3ph(mode1,mode3,mode2,iq) * conjg( sr%psi_3ph(mode1,mode3,mode2,iq) ) )
    !                         se%im_3ph(:,mode1)=se%im_3ph(:,mode1)-sfun*pref*psisq
    !                     endif
    !                 enddo mode1loop1
    !
    !             enddo mode3loop1
    !             enddo mode2loop1
    !
    !             ! Since I have the spectral functions available I might as well sort
    !             ! out the isotope scattering here.
    !             pref=isotope_prefactor*qp%ip(iq)%integration_weight
    !             do mode2=1,dr%n_mode
    !                 sfun=buf_spectral2(:,mode2)*pref
    !                 do mode1=1,dr%n_mode
    !                     se%im_iso(:,mode1)=se%im_iso(:,mode1)+sfun*sr%psi_iso(mode1,mode2,iq)
    !                 enddo
    !             enddo
    !
    !             if ( verbosity .gt. 0 .and. iq .lt. qp%n_irr_point ) then
    !                 t1=walltime()
    !                 call lo_progressbar(' ... imaginary selfenergy',iq,qp%n_irr_point,t1-t0)
    !             endif
    !         enddo qploopirr
    !     else
    !         ! Do not use symmetry
    !         if ( verbosity .gt. 0 ) call lo_progressbar_init()
    !         ctr=0
    !         qploopfull: do iq=1,qp%n_full_point
    !             ! Make it MPI parallel
    !             if ( mod(iq,mw%n) .ne. mw%r ) cycle
    !
    !             ! We need the spectral functions at q' and q''. Two different ways to go about it.
    !             select case(se%integrationtype)
    !             case(4)
    !                 ! Here I calculate the spectral functions from the interpolated self-energy, and then smear them appropriately.
    !
    !                 ! First thing, I need the spectral functions at q'
    !                 call ise%diagonal_selfenergy(p,qp%ap(iq)%r,dr%aq(iq)%omega,dr%aq(iq)%egv,bufRe,bufIm)
    !                 bufIm=max(bufIm,minIm)
    !                 do mode2=1,se%n_mode
    !                     if ( dr%aq(iq)%omega(mode2) .gt. lo_freqtol ) then
    !                         bufIm(:,mode2)=bufIm(:,mode2)*taperfn
    !                         call evaluate_spectral_function(ise%energy,bufIm(:,mode2),bufRe(:,mode2),dr%aq(iq)%omega(mode2),buf_spectral2(:,mode2))
    !                         sigma=qp%adaptive_sigma(qp%ap(iq)%radius,dr%aq(iq)%vel(:,mode2),dr%default_smearing(mode2),se%smearing_prefactor)
    !                         call gaussian_smear_spectral_function(ise%energy,sigma,buf_spectral2(:,mode2))
    !                         buf_spectral2(:,mode2)=buf_spectral2(:,mode2)/lo_trapezoid_integration(ise%energy,buf_spectral2(:,mode2))
    !                     else
    !                         buf_spectral2(:,mode2)=0.0_r8
    !                     endif
    !                 enddo
    !
    !                 ! Then the spectral functions at q''
    !                 call ise%diagonal_selfenergy(p,sr%qvec3(:,iq),sr%omega3(:,iq),sr%egv3(:,:,iq),bufRe,bufIm)
    !                 bufIm=max(bufIm,minIm)
    !                 do mode3=1,se%n_mode
    !                     if ( sr%omega3(mode3,iq) .gt. lo_freqtol ) then
    !                         bufIm(:,mode3)=bufIm(:,mode3)*taperfn
    !                         call evaluate_spectral_function(ise%energy,bufIm(:,mode3),bufRe(:,mode3),sr%omega3(mode3,iq),buf_spectral3(:,mode3))
    !                         sigma=qp%adaptive_sigma(qp%ap(iq)%radius,sr%vel3(:,mode3,iq),dr%default_smearing(mode3),se%smearing_prefactor)
    !                         call gaussian_smear_spectral_function(ise%energy,sigma,buf_spectral3(:,mode3))
    !                         buf_spectral3(:,mode3)=buf_spectral3(:,mode3)/lo_trapezoid_integration(ise%energy,buf_spectral3(:,mode3))
    !                     else
    !                         buf_spectral3(:,mode3)=0.0_r8
    !                     endif
    !                 enddo
    !             case(5)
    !                 ! This is the closed-grid case, where everything has been prepared beforehand.
    !                 ! Sort out q' first
    !                 jq=qp%ap(iq)%irreducible_index
    !                 buf_spectral2=isf%spectralfunction(:,:,jq)
    !
    !                 ! Sort out q''
    !                 v0=matmul(p%inv_reciprocal_latticevectors,sr%qvec3(:,iq))
    !                 jq=index_on_grid(qp,v0)
    !                 jq=qp%ap(jq)%irreducible_index
    !                 buf_spectral3=isf%spectralfunction(:,:,jq)
    !             end select
    !
    !             ! Prefactor for this q-point
    !             pref=threephonon_prefactor*qp%ap(iq)%integration_weight
    !             ! Excellent. Now start convoluting every which way.
    !             mode2loop2: do mode2=1,dr%n_mode
    !             mode3loop2: do mode3=1,dr%n_mode
    !                 ! ! Prepare for convolution
    !                 buf_odd_sf2a=0.0_r8
    !                 buf_odd_sf2b=0.0_r8
    !                 buf_odd_sf3a=0.0_r8
    !                 buf_odd_sf3b=0.0_r8
    !                 do ie=2,se%n_energy
    !                     buf_odd_sf2a(ie)=buf_spectral2(ie,mode2)
    !                     buf_odd_sf3a(ie)=buf_spectral3(ie,mode3)
    !                     buf_odd_sf2a(1-ie)=-buf_odd_sf2a(ie)
    !                     buf_odd_sf3a(1-ie)=-buf_odd_sf3a(ie)
    !                 enddo
    !                 buf_odd_sf2b=buf_odd_sf2a*buf_fpre
    !                 buf_odd_sf3b=buf_odd_sf3a*buf_fpre
    !
    !                 call lo_abcd_convolution(buf_odd_sf2a,buf_odd_sf3b,buf_odd_sf3a,buf_odd_sf2b,se%energy_axis(2)-se%energy_axis(1),buf_conv)
    !                 ! Post-convolution the values are FFT-shifted oddly so collect to the right place.
    !                 sfun=0.0_r8
    !                 do ie=1,se%n_energy
    !                     sfun(ie)=buf_conv(se%n_energy+1-ie)
    !                 enddo
    !
    !                 ! Excellent, everything was convoluted properly. Start accumulating self-energy
    !                 mode1loop2: do mode1=1,dr%n_mode
    !                     psisq=abs( sr%psi_3ph(mode1,mode2,mode3,iq) * conjg( sr%psi_3ph(mode1,mode2,mode3,iq) ) )
    !                     se%im_3ph(:,mode1)=se%im_3ph(:,mode1)-sfun*pref*psisq
    !                 enddo mode1loop2
    !
    !             enddo mode3loop2
    !             enddo mode2loop2
    !
    !             ! Since I have the spectral functions available I might as well sort
    !             ! out the isotope scattering here.
    !             pref=isotope_prefactor*qp%ap(iq)%integration_weight
    !             do mode2=1,dr%n_mode
    !                 sfun=buf_spectral2(:,mode2)*pref
    !                 do mode1=1,dr%n_mode
    !                     se%im_iso(:,mode1)=se%im_iso(:,mode1)+sfun*sr%psi_iso(mode1,mode2,iq)
    !                 enddo
    !             enddo
    !
    !             if ( verbosity .gt. 0 .and. iq .lt. qp%n_full_point ) then
    !                 t1=walltime()
    !                 call lo_progressbar(' ... imaginary selfenergy',iq,qp%n_full_point,t1-t0)
    !             endif
    !         enddo qploopfull
    !     endif ! use symmetry
    !
    !     call mem%deallocate(sfun,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    !     call mem%deallocate(taperfn,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    !     call mem%deallocate(bufRe,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    !     call mem%deallocate(bufIm,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    !     call mem%deallocate(buf_spectral2,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    !     call mem%deallocate(buf_spectral3,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    !
    !     ! Add together
    !     call mw%allreduce('sum',se%im_3ph)
    !     call mw%allreduce('sum',se%im_iso)
    !
    !     if ( verbosity .gt. 0 ) then
    !         t1=walltime()
    !         call lo_progressbar(' ... imaginary selfenergy',qp%n_full_point,qp%n_full_point,t1-t0)
    !         t0=t1
    !     endif
    ! end block actualcalculation

    ! Make sure the whole thing makes sense in the end.
    massage: block
        real(r8), dimension(:), allocatable :: taper, buf1, buf2
        integer :: mode1, mode2, i

        call mem%allocate(taper, se%n_energy, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(buf1, se%n_energy, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(buf2, se%n_energy, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        taper = 0.0_r8
        buf1 = 0.0_r8
        buf2 = 0.0_r8

        ! Sort out degeneracies
        do mode1 = 1, dr%n_mode
            buf1 = 0.0_r8
            buf2 = 0.0_r8
            do i = 1, wp%degeneracy(mode1)
                mode2 = wp%degenmode(i, mode1)
                buf1 = buf1 + se%im_iso(:, mode2)
                buf2 = buf2 + se%im_3ph(:, mode2)
            end do
            buf1 = buf1/real(wp%degeneracy(mode1), r8)
            buf2 = buf2/real(wp%degeneracy(mode1), r8)
            do i = 1, wp%degeneracy(mode1)
                mode2 = wp%degenmode(i, mode1)
                se%im_iso(:, mode2) = buf1
                se%im_3ph(:, mode2) = buf2
            end do
        end do
        ! Get the tapering function
        call taperfn_im(se%energy_axis, dr%omega_max, dr%omega_min, taper)
        ! Taper the self-energies so that they are zero where they should be.
        do mode1 = 1, dr%n_mode
            se%im_3ph(:, mode1) = max(se%im_3ph(:, mode1), 0.0_r8)
            se%im_iso(:, mode1) = max(se%im_iso(:, mode1), 0.0_r8)
            se%im_3ph(:, mode1) = se%im_3ph(:, mode1)*taper
            se%im_iso(:, mode1) = se%im_iso(:, mode1)*taper
        end do

        call mem%deallocate(taper, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(buf1, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(buf2, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    end block massage
end subroutine
