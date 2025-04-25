submodule(lo_evaluate_phonon_self_energy) lo_evaluate_phonon_self_energy_fourphonon
implicit none
contains

!> The fourth order self-energy
module subroutine fourphonon_real_selfenergy(qpoint, wp, gp, qp, uc, temperature, dr, fcf, delta, skipsym, mw, mem, verbosity)
    !> qpoint for q
    type(lo_qpoint), intent(in) :: qpoint
    !> harmonic properties for q
    type(lo_phonon_dispersions_qpoint), intent(in) :: wp
    !> harmonic properties for Gamma, with proper direction
    type(lo_phonon_dispersions_qpoint), intent(in) :: gp
    !> grid for q',q'',q'''
    class(lo_qpoint_mesh), intent(in) :: qp
    !> unit cell
    type(lo_crystalstructure), intent(in) :: uc
    !> harmonic properties for q',q'',q''
    type(lo_phonon_dispersions), intent(in) :: dr
    !> third order force constant
    type(lo_forceconstant_fourthorder), intent(in) :: fcf
    !> temperature
    real(r8), intent(in) :: temperature
    !> real four-phonon self-energy
    real(r8), dimension(:), intent(out) :: delta
    !> skip symmetry for some reason?
    logical, intent(in) :: skipsym
    !> mpi helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem
    !> how much to talk
    integer, intent(in) :: verbosity

    logical :: atgamma
    complex(r8), dimension(:, :, :), allocatable :: nuvec2
    complex(r8), dimension(:, :), allocatable :: nuvec1
    integer :: ind_gamma_irr, ind_gamma_full
    real(r8) :: timer, t0, t1

    ! Some basic things
    init: block
        real(r8) :: f0, f1
        integer :: iq, imode, ialpha, iatom, ix, ctr

        timer = walltime()
        t0 = timer
        t1 = timer

        ! Check if we are at gamma?
        if (norm2(qpoint%r) .lt. lo_sqtol) then
            atgamma = .true.
        else
            atgamma = .false.
        end if
        ! Skip symmetry if on path.
        if (skipsym) atgamma = .false.

        ! Locate the index for Gamma on the q-grid, always useful
        ! to make sure we are consistent in everything, and that
        ! whenever I reference Gamma I pick a Gamma with the correct
        ! q-direction attached?
        ind_gamma_irr = -1
        ind_gamma_full = -1
        do iq = 1, qp%n_irr_point
            if (lo_sqnorm(qp%ip(iq)%r) .lt. lo_sqtol) then
                ! found Gamma
                ind_gamma_irr = iq
                ind_gamma_full = qp%ip(iq)%full_index
                exit
            end if
        end do

        ! Some temporary space for scaled eigenvectors
        if (atgamma) then
            call mem%allocate(nuvec1, [dr%n_mode, dr%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%allocate(nuvec2, [dr%n_mode, dr%n_mode, qp%n_irr_point], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            nuvec1 = 0.0_r8
            nuvec2 = 0.0_r8
        else
            call mem%allocate(nuvec1, [dr%n_mode, dr%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%allocate(nuvec2, [dr%n_mode, dr%n_mode, qp%n_full_point], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            nuvec1 = 0.0_r8
            nuvec2 = 0.0_r8
        end if

        ! First fix the one at the main q-point, always the same
        ! just pre-multiply with masses and frequencies
        do imode = 1, dr%n_mode
        do iatom = 1, uc%na
        do ix = 1, 3
            if (atgamma) then
                ! Get consistent Gamma
                if (gp%omega(imode) .gt. lo_freqtol) then
                    f0 = 1.0_r8/sqrt(gp%omega(imode))
                else
                    f0 = 0.0_r8
                end if
                ialpha = (iatom - 1)*3 + ix
                f1 = uc%invsqrtmass(iatom)
                nuvec1(ialpha, imode) = gp%egv(ialpha, imode)*f0*f1
            else
                ! Just normal copy
                if (wp%omega(imode) .gt. lo_freqtol) then
                    f0 = 1.0_r8/sqrt(wp%omega(imode))
                else
                    f0 = 0.0_r8
                end if
                ialpha = (iatom - 1)*3 + ix
                f1 = uc%invsqrtmass(iatom)
                nuvec1(ialpha, imode) = wp%egv(ialpha, imode)*f0*f1
            end if
        end do
        end do
        end do

        ! Pre-multiply eigenvectors with masses and frequencies
        if (atgamma) then
            ! Case that will use symmetry eventually
            ctr = 0
            do iq = 1, qp%n_irr_point
            do imode = 1, dr%n_mode
                ctr = ctr + 1
                if (mod(ctr, mw%n) .ne. mw%r) cycle
                if (iq .eq. ind_gamma_irr) then
                    if (gp%omega(imode) .gt. lo_freqtol) then
                        f0 = 1.0_r8/sqrt(gp%omega(imode))
                    else
                        f0 = 0.0_r8
                    end if

                    do iatom = 1, uc%na
                        f1 = uc%invsqrtmass(iatom)
                        do ix = 1, 3
                            ialpha = (iatom - 1)*3 + ix
                            nuvec2(ialpha, imode, iq) = gp%egv(ialpha, imode)*f0*f1
                        end do
                    end do
                else
                    if (dr%iq(iq)%omega(imode) .gt. lo_freqtol) then
                        f0 = 1.0_r8/sqrt(dr%iq(iq)%omega(imode))
                    else
                        f0 = 0.0_r8
                    end if

                    do iatom = 1, uc%na
                        f1 = uc%invsqrtmass(iatom)
                        do ix = 1, 3
                            ialpha = (iatom - 1)*3 + ix
                            nuvec2(ialpha, imode, iq) = dr%iq(iq)%egv(ialpha, imode)*f0*f1
                        end do
                    end do
                end if
            end do
            end do
            call mw%allreduce('sum', nuvec2)
        else
            ! Case that will not use symmetry
            ctr = 0
            do iq = 1, qp%n_full_point
            do imode = 1, dr%n_mode
                ctr = ctr + 1
                if (mod(ctr, mw%n) .ne. mw%r) cycle

                if (iq .eq. ind_gamma_full) then
                    ! pick consistent Gamma
                    if (gp%omega(imode) .gt. lo_freqtol) then
                        f0 = 1.0_r8/sqrt(gp%omega(imode))
                    else
                        f0 = 0.0_r8
                    end if
                    do iatom = 1, uc%na
                        f1 = uc%invsqrtmass(iatom)
                        do ix = 1, 3
                            ialpha = (iatom - 1)*3 + ix
                            nuvec2(ialpha, imode, iq) = gp%egv(ialpha, imode)*f0*f1
                        end do
                    end do
                else
                    ! Not at Gamma, don't have to care
                    if (dr%aq(iq)%omega(imode) .gt. lo_freqtol) then
                        f0 = 1.0_r8/sqrt(dr%aq(iq)%omega(imode))
                    else
                        f0 = 0.0_r8
                    end if
                    do iatom = 1, uc%na
                        f1 = uc%invsqrtmass(iatom)
                        do ix = 1, 3
                            ialpha = (iatom - 1)*3 + ix
                            nuvec2(ialpha, imode, iq) = dr%aq(iq)%egv(ialpha, imode)*f0*f1
                        end do
                    end do
                end if
            end do
            end do
            call mw%allreduce('sum', nuvec2)
        end if

        if (verbosity .gt. 0) then
            t1 = walltime()
            write (*, *) "... pretransformed and made temporary space (", tochar(t1 - t0), "s)"
            t0 = t1
        end if
    end block init

    ! New, fun way of doing things
    new: block
        complex(r8), dimension(:, :), allocatable :: bf1, bf2
        complex(r8), dimension(:), allocatable :: evp3, ptf
        real(r8), dimension(dr%n_mode) :: bfom2
        real(r8) :: prefactor, psi, omegathres, f0
        integer :: iq, b1, b2, ctr, i

        ! Init some things
        delta = 0.0_r8
        omegathres = dr%omega_min*0.5_r8

        call mem%allocate(bf1, [dr%n_mode**2, dr%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(bf2, [dr%n_mode**2, dr%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(evp3, dr%n_mode**4, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(ptf, dr%n_mode**4, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        bf1 = 0.0_r8
        bf2 = 0.0_r8
        evp3 = 0.0_r8
        ptf = 0.0_r8
        ctr = 0

        if (verbosity .gt. 0) call lo_progressbar_init()

        ! I can always pre-transform the first guy
        bf1 = 0.0_r8
        do b1 = 1, dr%n_mode
            if (wp%omega(b1) .gt. omegathres) then
                call zgerc(dr%n_mode, dr%n_mode, (1.0_r8, 0.0_r8), nuvec1(:, b1), 1, nuvec1(:, b1), 1, bf1(:, b1), dr%n_mode)
            end if
        end do

        if (atgamma) then
            ! Can not use symmetry
            do iq = 1, qp%n_irr_point
                ! Make it parallel
                ctr = ctr + 1
                if (mod(ctr, mw%n) .ne. mw%r) cycle

                prefactor = qp%ip(iq)%integration_weight*fourphonon_prefactor
                ! pre-transform the matrix element
                call pretransform_phi4(fcf, qpoint%r, qp%ip(iq)%r, ptf)

                ! pre-fetch frequencies
                bfom2 = dr%iq(iq)%omega

                ! first outer product
                bf2 = 0.0_r8
                do b1 = 1, dr%n_mode
                    if (bfom2(b1) .gt. omegathres) then
                        call zgerc(dr%n_mode, dr%n_mode, (1.0_r8, 0.0_r8), nuvec2(:, b1, iq), 1, nuvec2(:, b1, iq), 1, bf2(:, b1), dr%n_mode)
                        ! and multiply in the prefactor here?
                        bf2(:, b1) = bf2(:, b1)*(1 + 2*lo_planck(temperature, bfom2(b1)))*prefactor
                    end if
                end do

                ! second larger outer product to get matrix elements
                do b1 = 1, dr%n_mode
                    if (wp%omega(b1) .lt. omegathres) cycle
                    do b2 = 1, dr%n_mode
                        if (bfom2(b2) .lt. omegathres) cycle
                        evp3 = 0.0_r8
                        call zgeru(dr%n_mode**2, dr%n_mode**2, (1.0_r8, 0.0_r8), bf2(:, b2), 1, bf1(:, b1), 1, evp3, dr%n_mode**2)
                        psi = real(dot_product(evp3, ptf), r8)
                        delta(b1) = delta(b1) + psi
                    end do
                end do

                if (verbosity .gt. 0) then
                    if (lo_trueNtimes(iq, 127, qp%n_irr_point)) call lo_progressbar(' ... fourphonon self-energy', iq, qp%n_irr_point, walltime() - t0)
                end if
            end do
        else
            ! Can not use symmetry
            do iq = 1, qp%n_full_point
                ! Make it parallel
                ctr = ctr + 1
                if (mod(ctr, mw%n) .ne. mw%r) cycle

                prefactor = qp%ap(iq)%integration_weight*fourphonon_prefactor
                ! pre-transform the matrix element
                call pretransform_phi4(fcf, qpoint%r, qp%ap(iq)%r, ptf)

                ! pre-fetch frequencies
                bfom2 = dr%aq(iq)%omega

                ! first outer product
                bf2 = 0.0_r8
                do b1 = 1, dr%n_mode
                    if (bfom2(b1) .gt. omegathres) then
                        call zgerc(dr%n_mode, dr%n_mode, (1.0_r8, 0.0_r8), nuvec2(:, b1, iq), 1, nuvec2(:, b1, iq), 1, bf2(:, b1), dr%n_mode)
                        ! and multiply in the prefactor here?
                        bf2(:, b1) = bf2(:, b1)*(1 + 2*lo_planck(temperature, bfom2(b1)))*prefactor
                    end if
                end do

                ! second larger outer product to get matrix elements
                do b1 = 1, dr%n_mode
                    if (wp%omega(b1) .lt. omegathres) cycle
                    do b2 = 1, dr%n_mode
                        if (bfom2(b2) .lt. omegathres) cycle
                        evp3 = 0.0_r8
                        call zgeru(dr%n_mode**2, dr%n_mode**2, (1.0_r8, 0.0_r8), bf2(:, b2), 1, bf1(:, b1), 1, evp3, dr%n_mode**2)
                        psi = real(dot_product(evp3, ptf), r8)
                        delta(b1) = delta(b1) + psi
                    end do
                end do

                if (verbosity .gt. 0) then
                    if (lo_trueNtimes(iq, 127, qp%n_full_point)) call lo_progressbar(' ... fourphonon self-energy', iq, qp%n_full_point, walltime() - t0)
                end if
            end do
        end if

        call mw%allreduce('sum', delta)
        ! Fix degeneracies
        do b1 = 1, dr%n_mode
            f0 = 0.0_r8
            do i = 1, wp%degeneracy(b1)
                b2 = wp%degenmode(i, b1)
                f0 = f0 + delta(b2)
            end do
            f0 = f0/real(wp%degeneracy(b1), r8)
            do i = 1, wp%degeneracy(b1)
                b2 = wp%degenmode(i, b1)
                delta(b2) = f0
            end do
        end do

        call mem%deallocate(bf1, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(bf2, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(evp3, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(ptf, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)

        if (verbosity .gt. 0) then
            call lo_progressbar(' ... fourphonon self-energy', qp%n_full_point, qp%n_full_point, walltime() - t0)
        end if
    end block new

    ! dumdum: block
    !     integer, parameter :: nb=4
    !     real(r8), dimension(nb) :: v0,v1,v2,v3
    !     real(r8), dimension(nb**2) :: pr0,pr1
    !     real(r8), dimension(nb**4) :: fl0,fl1
    !     integer :: i,j,k,l,m,ctr,ii,jj
    !
    !     ctr=0
    !     do i=1,nb
    !         v0(i)=i
    !         v1(i)=(i+nb)
    !         v2(i)=(i+2*nb)
    !         v3(i)=(i+3*nb)
    !     enddo
    !     pr0=0.0_r8
    !     pr1=0.0_r8
    !     fl0=0.0_r8
    !     fl1=0.0_r8
    !     call dger(nb, nb, 1.0_r8, v3, 1, v2, 1, pr0, nb)
    !     call dger(nb, nb, 1.0_r8, v1, 1, v0, 1, pr1, nb)
    !     call dger(nb**2, nb**2, 1.0_r8, pr0, 1, pr1, 1, fl0, nb**2)
    !
    !     do i=1,nb
    !     do j=1,nb
    !     do k=1,nb
    !     do l=1,nb
    !         m=(i-1)*nb*nb*nb + (j-1)*nb*nb + (k-1)*nb + l
    !         fl1(m)=v0(i)*v1(j)*v2(k)*v3(l)
    !         !,v0(i)*v1(j),v2(k)*v3(l)
    !     enddo
    !     enddo
    !     write(*,*) tochar([i,j,k,l,m]),fl1(m),fl0(m),fl1(m)-fl0(m)
    !     enddo
    !     enddo
    !
    !     write(*,*) ''
    !     do i=1,nb
    !     do j=1,nb
    !         k=(i-1)*nb+j
    !         write(*,*) i,j,k,pr0(k),v0(i)*v1(j)
    !     enddo
    !     enddo
    !     write(*,*) ''
    ! end block dumdum

    call mem%deallocate(nuvec1, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(nuvec2, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)

    ! ! Old, solid way of doing it
    ! old: block
    !     complex(r8), dimension(:,:), allocatable :: egv
    !     real(r8), dimension(4) :: omega
    !     real(r8), dimension(3) :: qv1,qv2,qv3,qv4
    !     real(r8) :: f0,t0,omegathres
    !     integer :: i,q,b1,b2,b3,b4,ctr
    !
    !     ! Init some things
    !     delta=0.0_r8
    !     t0=walltime()
    !     omegathres=dr%omega_min*0.5_r8
    !     call mem%allocate(egv,[dr%n_mode,4],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    !     egv=0.0_r8
    !
    !     if ( verbosity .gt. 0 ) call lo_progressbar_init()
    !     ctr=0
    !     ! Get harmonic things at negative q1
    !     do q=1,qp%n_full_point
    !     do b1=1,dr%n_mode
    !         ! MPI thing
    !         ctr=ctr+1
    !         if ( mod(ctr,mw%n) .ne. mw%r ) cycle
    !         ! Get the q-vectors
    !         qv1=qpoint%r*lo_twopi
    !         qv2=-qv1
    !         qv3=qp%ap(q)%r*lo_twopi
    !         qv4=-qv3
    !         ! and the self-energy, eventually
    !         do b3=1,dr%n_mode
    !             b2=b1
    !             b4=b3
    !             ! fetch frequencies
    !             omega(1)=wp%omega(b1)
    !             omega(2)=wp%omega(b2)
    !             omega(3)=dr%aq(q)%omega(b3)
    !             omega(4)=dr%aq(q)%omega(b4)
    !             ! fetch eigenvectors
    !             egv(:,1)=wp%egv(:,b1)
    !             egv(:,2)=conjg(wp%egv(:,b2))
    !             egv(:,3)=dr%aq(q)%egv(:,b3)
    !             egv(:,4)=conjg(dr%aq(q)%egv(:,b4))
    !             ! scatteringrate + selfenergy in one go!
    !             if ( minval(omega) .gt. omegathres ) then
    !                 f0=real(fcf%scatteringamplitude(omega,egv,-qv2,-qv3,-qv4),r8)
    !                 f0=f0*(1.0_r8+lo_planck(temperature,omega(3))+lo_planck(temperature,omega(4)) )
    !                 f0=f0*qp%ap(q)%integration_weight*fourphonon_prefactor
    !             else
    !                 f0=0.0_r8
    !             endif
    !             delta(b1)=delta(b1)+f0
    !         enddo
    !
    !         if ( verbosity .gt. 0 ) then
    !             if ( lo_trueNtimes(ctr,127,qp%n_full_point*dr%n_mode) ) call lo_progressbar(' ... fourphonon self-energy',ctr,dr%n_mode*qp%n_full_point,walltime()-t0)
    !         endif
    !     enddo
    !     enddo
    !     ! Intermediate cleanup
    !     call mem%deallocate(egv,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    !
    !     ! Add it together
    !     call mw%allreduce('sum',delta)
    !     ! Fix degeneracies
    !     do b1=1,dr%n_mode
    !         f0=0.0_r8
    !         do i=1,wp%degeneracy(b1)
    !             b2=wp%degenmode(i,b1)
    !             f0=f0+delta(b2)
    !         enddo
    !         f0=f0/real(wp%degeneracy(b1),r8)
    !         do i=1,wp%degeneracy(b1)
    !             b2=wp%degenmode(i,b1)
    !             delta(b2)=f0
    !         enddo
    !     enddo
    !
    !     if ( verbosity .gt. 0 ) call lo_progressbar(' ... fourphonon self-energy',dr%n_mode*qp%n_full_point,dr%n_mode*qp%n_full_point,walltime()-t0)
    !
    !     if ( mw%talk ) then
    !         do b1=1,dr%n_mode
    !             write(*,*) 'old',b1,delta(b1)
    !         enddo
    !     endif
    ! end block old
end subroutine

!> pre-transform to get half of the fourth order matrix element
subroutine pretransform_phi4(fcf, q1, q2, ptf)
    !> third order forceconstant
    type(lo_forceconstant_fourthorder), intent(in) :: fcf
    !> q-vectors
    real(r8), dimension(3), intent(in) :: q1, q2
    !> flattened, pretransformed matrix element
    complex(r8), dimension(:), intent(out) :: ptf

    integer :: i, j, k, l, m

    complex(r8) :: expiqr
    real(r8), dimension(3) :: rv2, rv3, rv4
    real(r8) :: iqr
    integer :: a1, a2, a3, a4, ia, ib, ic, id, q, nb

    nb = fcf%na*3
    ptf = 0.0_r8
    do a1 = 1, fcf%na
    do q = 1, fcf%atom(a1)%n
        a2 = fcf%atom(a1)%quartet(q)%i2
        a3 = fcf%atom(a1)%quartet(q)%i3
        a4 = fcf%atom(a1)%quartet(q)%i4

        rv2 = fcf%atom(a1)%quartet(q)%lv2
        rv3 = fcf%atom(a1)%quartet(q)%lv3
        rv4 = fcf%atom(a1)%quartet(q)%lv4

        iqr = -dot_product(q1, rv2) + dot_product(q2, rv3) - dot_product(q2, rv4)
        iqr = -iqr*lo_twopi
        expiqr = cmplx(cos(iqr), sin(iqr), r8)
        do l = 1, 3
        do k = 1, 3
        do j = 1, 3
        do i = 1, 3
            ia = (a1 - 1)*3 + i
            ib = (a2 - 1)*3 + j
            ic = (a3 - 1)*3 + k
            id = (a4 - 1)*3 + l
            ! Now for the grand flattening scheme, consistent with the zgeru operations above.
            m = (ia - 1)*nb*nb*nb + (ib - 1)*nb*nb + (ic - 1)*nb + id
            ptf(m) = ptf(m) + fcf%atom(a1)%quartet(q)%m(i, j, k, l)*expiqr
        end do
        end do
        end do
        end do
    end do
    end do
end subroutine

end submodule