module lo_epot
    !! Deal with many kinds of potential energy differences
    use konstanter, only: r8, lo_pi, lo_twopi, lo_tol, lo_sqtol, lo_status, lo_Hartree_to_eV, lo_kb_hartree
    use gottochblandat, only: tochar, walltime, lo_chop, lo_trueNtimes, lo_progressbar_init, &
                              lo_progressbar, lo_frobnorm, open_file, lo_flattentensor, lo_sqnorm, lo_outerproduct, lo_mean, &
                              lo_points_on_sphere, lo_mean, lo_stddev
    use mpi_wrappers, only: lo_mpi_helper, lo_stop_gracefully
    use lo_memtracker, only: lo_mem_helper
    use type_crystalstructure, only: lo_crystalstructure
    use type_forceconstant_secondorder, only: lo_forceconstant_secondorder
    use type_forceconstant_thirdorder, only: lo_forceconstant_thirdorder
    use type_forceconstant_fourthorder, only: lo_forceconstant_fourthorder
    use type_mdsim, only: lo_mdsim
    implicit none
    
    private
    public :: lo_energy_differences
    
    type lo_energy_differences
        ! forceconstants needed for evaluation
        type(lo_forceconstant_secondorder) :: fc2
        type(lo_forceconstant_thirdorder) :: fc3
        type(lo_forceconstant_fourthorder) :: fc4
        real(r8), dimension(:, :, :, :), allocatable :: fcp
        logical :: polar = .false.
    
    contains
        procedure :: setup => setup_potential_energy_differences
        procedure :: statistical_sampling
        procedure :: energies_and_forces
    end type
    
    contains
    
    !> statistically sample
    subroutine statistical_sampling(pot, uc, ss, fc, nstep, temperature, quantum, ebuf, mw, mem, verbosity)
        !> container for potential energy differences
        class(lo_energy_differences), intent(inout) :: pot
        !> unitcell
        type(lo_crystalstructure), intent(in) :: uc
        !> supercell
        type(lo_crystalstructure), intent(in) :: ss
        !> second order forceconstant
        type(lo_forceconstant_secondorder), intent(in) :: fc
        !> number of configs to generate
        integer, intent(in) :: nstep
        !> temperature
        real(r8), intent(in) :: temperature
        !> quantum configurations?
        logical :: quantum
        !> storage for potential energies
        real(r8), dimension(:, :), intent(out) :: ebuf
        !> MPI helper
        type(lo_mpi_helper), intent(inout) :: mw
        !> memory tracker
        type(lo_mem_helper), intent(inout) :: mem
        !> talk a lot?
        integer, intent(in) :: verbosity
    
        type(lo_crystalstructure) :: p
        integer :: ctr, i
        real(r8), dimension(:, :), allocatable :: f2, f3, f4, fp
        real(r8) :: e2, e3, e4, ep
    
        ! Copy of structure to work with
        p = ss
    
        ebuf = 0.0_r8
    
        ! Dummy space for force
        call mem%allocate(f2, [3, ss%na], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(f3, [3, ss%na], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(f4, [3, ss%na], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(fp, [3, ss%na], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        f2 = 0.0_r8
        f3 = 0.0_r8
        f4 = 0.0_r8
        fp = 0.0_r8
        
        do i = 1, nstep

            if (mod(i, mw%n) .ne. mw%r) cycle

            ! Reset structure
            p%u = 0.0_r8
            p%v = 0.0_r8
            p%r = ss%r
            ! Get new structure
            call pot%fc2%initialize_cell(p, uc, fc, temperature, quantum, .false., -1.0_r8, mw, nosync=.true.)
            ! Calculate the energy
            call pot%energies_and_forces(p%u, e2, e3, e4, ep, f2, f3, f4, fp)

            ! ek = p%kinetic_energy()/(p%na)

            ebuf(i, 1) = e2
            ebuf(i, 2) = e3
            ebuf(i, 3) = e4
            ebuf(i, 4) = ep            
        end do

        call mw%allreduce('sum', ebuf)
            
        call mem%deallocate(f2, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(f3, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(f4, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(fp, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    end subroutine
    
    ! !> calculate energies and forces for a given configuration
    subroutine energies_and_forces(pot, u, e2, e3, e4, ep, f2, f3, f4, fp)
        !> container for potential energy differences
        class(lo_energy_differences), intent(in) :: pot
        !> displacements
        real(r8), dimension(:, :), intent(in) :: u
        !> energies
        real(r8), intent(out) :: e2, e3, e4, ep
        !> forces
        real(r8), dimension(:, :), intent(out) :: f2, f3, f4, fp
    
        real(r8), dimension(3, 3, 3, 3) :: m4
        real(r8), dimension(3, 3, 3) :: m3
        real(r8), dimension(3, 3) :: m2
        real(r8), dimension(3) :: v0, u2, u3, u4
        integer :: a1, a2, a3, a4, i1, i2, i3, i4, i
    
        e2 = 0.0_r8
        e3 = 0.0_r8
        e4 = 0.0_r8
        ep = 0.0_r8
        f2 = 0.0_r8
        f3 = 0.0_r8
        f4 = 0.0_r8
        fp = 0.0_r8
    
        ! pair term
        do a1 = 1, size(u, 2)
            v0 = 0.0_r8
            do i1 = 1, pot%fc2%atom(a1)%n
                a2 = pot%fc2%atom(a1)%pair(i1)%i2
                m2 = pot%fc2%atom(a1)%pair(i1)%m
                v0 = v0 - matmul(m2, u(:, a2))
            end do
            e2 = e2 - dot_product(u(:, a1), v0)*0.5_r8
            f2(:, a1) = v0
        end do
    
        ! polar term
        if (pot%polar) then
            ep = 0.0_r8
            do a1 = 1, size(u, 2)
                v0 = 0.0_r8
                do a2 = 1, size(u, 2)
                    v0 = v0 - matmul(pot%fcp(:, :, a1, a2), u(:, a2))
                end do
                ep = ep - dot_product(u(:, a1), v0)*0.5_r8
                fp(:, a1) = v0
            end do
        end if
    
        ! triplet term
        if (pot%fc3%na .eq. size(u, 2)) then
            do a1 = 1, size(u, 2)
                v0 = 0.0_r8
                do i = 1, pot%fc3%atom(a1)%n
                    m3 = pot%fc3%atom(a1)%triplet(i)%m
                    a2 = pot%fc3%atom(a1)%triplet(i)%i2
                    a3 = pot%fc3%atom(a1)%triplet(i)%i3
                    u2 = u(:, a2)
                    u3 = u(:, a3)
                    do i1 = 1, 3
                    do i2 = 1, 3
                    do i3 = 1, 3
                        v0(i1) = v0(i1) - m3(i1, i2, i3)*u2(i2)*u3(i3)
                    end do
                    end do
                    end do
                end do
                v0 = v0*0.5_r8
                f3(:, a1) = v0
                e3 = e3 - dot_product(v0, u(:, a1))/3.0_r8
            end do
        end if
    
        ! quartet term
        if (pot%fc4%na .eq. size(u, 2)) then
            do a1 = 1, size(u, 2)
                v0 = 0.0_r8
                do i = 1, pot%fc4%atom(a1)%n
                    m4 = pot%fc4%atom(a1)%quartet(i)%m
                    a2 = pot%fc4%atom(a1)%quartet(i)%i2
                    a3 = pot%fc4%atom(a1)%quartet(i)%i3
                    a4 = pot%fc4%atom(a1)%quartet(i)%i4
                    u2 = u(:, a2)
                    u3 = u(:, a3)
                    u4 = u(:, a4)
                    do i1 = 1, 3
                    do i2 = 1, 3
                    do i3 = 1, 3
                    do i4 = 1, 3
                        v0(i1) = v0(i1) - m4(i1, i2, i3, i4)*u2(i2)*u3(i3)*u4(i4)
                    end do
                    end do
                    end do
                    end do
                end do
                v0 = v0/6.0_r8
                f4(:, a1) = v0
                e4 = e4 - dot_product(v0, u(:, a1))/4.0_r8
            end do
        end if
    end subroutine
    
    !> Calculate potential energy differences in several ways
    subroutine setup_potential_energy_differences(pot, uc, ss, fc2, fc3, fc4, mw, verbosity)
        !> container for potential energy differences
        class(lo_energy_differences), intent(out) :: pot
        !> unitcell
        type(lo_crystalstructure), intent(inout) :: uc
        !> supercell
        type(lo_crystalstructure), intent(inout) :: ss
        !> second order forceconstant
        type(lo_forceconstant_secondorder), intent(in) :: fc2
        !> third order forceconstant
        type(lo_forceconstant_thirdorder), intent(in) :: fc3
        !> fourth order forceconstant
        type(lo_forceconstant_fourthorder), intent(in) :: fc4
        !> mpi things
        type(lo_mpi_helper), intent(inout) :: mw
        !> how much to talk
        integer, intent(in) :: verbosity
    
        real(r8) :: timer, t0, t1
    
        timer = walltime()
        t0 = timer
        t1 = timer
    
        if (verbosity .gt. 0) then
            write (*, *) ''
            write (*, *) 'PREPARING POTENTIAL ENERGY DIFFERENCES'
        end if
    
        call fc2%remap(uc, ss, pot%fc2)
        if (verbosity .gt. 0) then
            t1 = walltime()
            write (*, *) '... remapped second order (', tochar(t1 - t0), ')'
            t0 = t1
        end if
    
        if (fc3%na .gt. 0) then
            call fc3%remap(uc, ss, pot%fc3)
            if (verbosity .gt. 0) then
                t1 = walltime()
                write (*, *) '... remapped third order (', tochar(t1 - t0), ')'
                t0 = t1
            end if
        end if
    
        if (fc4%na .gt. 0) then
            call fc4%remap(uc, ss, pot%fc4)
            if (verbosity .gt. 0) then
                t1 = walltime()
                write (*, *) '... remapped fourth order (', tochar(t1 - t0), ')'
                t0 = t1
            end if
        end if
    
        if (fc2%polar) then
            allocate (pot%fcp(3, 3, ss%na, ss%na))
            pot%fcp = 0.0_r8
            call fc2%supercell_longrange_dynamical_matrix_at_gamma(ss, pot%fcp, 1E-15_r8)
            pot%polar = .true.
        else
            pot%polar = .false.
        end if
        if (verbosity .gt. 0) then
            t1 = walltime()
            write (*, *) '... built polar forceconstant (', tochar(t1 - t0), ')'
            t0 = t1
        end if
    end subroutine
    
    ! !> remove longrange interactions
    ! subroutine potential_energy_difference(sim,uc,ss,fc,fct,fcf,Uref,mw,verbosity)
    !     !> force-displacement data
    !     type(lo_mdsim), intent(inout) :: sim
    !     !> unitcell
    !     type(lo_crystalstructure), intent(inout) :: uc
    !     !> supercell
    !     type(lo_crystalstructure), intent(inout) :: ss
    !     !> second order forceconstant
    !     type(lo_forceconstant_secondorder), intent(in) :: fc
    !     !> third order forceconstant
    !     type(lo_forceconstant_thirdorder), intent(in) :: fct
    !     !> fourth order forceconstant
    !     type(lo_forceconstant_fourthorder), intent(in) :: fcf
    !     !> reference energy
    !     real(flyt), intent(in) :: Uref
    !     !> mpi things
    !     type(lo_mpi_helper), intent(inout) :: mw
    !     !> how much to talk
    !     integer, intent(in) :: verbosity
    !
    !     type(lo_forceconstant_secondorder) :: fcss
    !     type(lo_forceconstant_thirdorder) :: fctss
    !     type(lo_forceconstant_fourthorder) :: fcfss
    !     real(flyt), dimension(:,:,:,:), allocatable :: forceconstant
    !     real(flyt), dimension(:,:), allocatable :: fp,f2,f3,f4,de
    ! !    real(flyt), dimension(:,:), allocatable :: f,de
    ! !    real(flyt), dimension(3,3) :: m0,m1
    ! !    real(flyt), dimension(3) :: v0
    ! !    real(flyt) :: epot,epot_2,epot_3,epot_4,epot_dd
    ! !    integer :: a1,a2,t,u
    !
    !     init: block
    !         real(flyt), dimension(3,3) :: m0,m1
    !         real(flyt) :: t0,f0
    !         integer :: a1,a2
    !         if ( mw%talk ) then
    !             write(*,*) ''
    !             write(*,*) 'CALCULATING POTENTIAL ENERGY DIFFERENCES'
    !         endif
    !         ! Remap the forceconstants to the supercell
    !         call ss%classify('supercell',uc)
    !         call fc%remap(uc,ss,fcss)
    !         call fct%remap(uc,ss,fctss)
    !         call fcf%remap(uc,ss,fcfss)
    !         if ( mw%talk ) write(*,*) '... remapped forceconstants'
    !         ! Get the dipole-dipole forceconstant
    !         if ( fc%polar ) then
    !             lo_allocate(forceconstant(3,3,ss%na,ss%na))
    !             call fc%supercell_longrange_dynamical_matrix_at_gamma(ss,forceconstant,1E-15_flyt)
    !             ! small sanity-check, just in case
    !             f0=0.0_flyt
    !             do a1=1,ss%na
    !             do a2=a1,ss%na
    !                 m0=forceconstant(:,:,a1,a2)
    !                 m1=transpose(forceconstant(:,:,a2,a1))
    !                 f0=f0+lo_frobnorm(m0-m1)
    !             enddo
    !             enddo
    !             f0=f0/(ss%na**2)
    !             if ( f0 .gt. lo_tol ) then
    !                 write(*,*) 'ERROR, non-Hermitian electrostatic forceconstant: ',f0
    !                 write(*,*) 'If this keeps happening I should not hard-code the Ewald tolerance.'
    !                 stop
    !             endif
    !         endif
    !
    !         ! some space for energies and forces
    !         lo_allocate(de(sim%nt,5))
    !         lo_allocate(fp(3,ss%na))
    !         lo_allocate(f2(3,ss%na))
    !         lo_allocate(f3(3,ss%na))
    !         lo_allocate(f4(3,ss%na))
    !         de=0.0_flyt
    !         fp=0.0_flyt
    !         f2=0.0_flyt
    !         f3=0.0_flyt
    !         f4=0.0_flyt
    !     end block init
    !
    !     forceenergy: block
    !         real(flyt), dimension(3,3,3,3) :: m4
    !         real(flyt), dimension(3,3,3) :: m3
    !         real(flyt), dimension(3,3) :: m2
    !         real(flyt), dimension(3) :: u1,u2,u3,u4,v0
    !         real(flyt) :: energy
    !         integer :: a1,a2,a3,a4,i1,i2,i3,i4,i,j,k,l,t
    !
    !         do t=1,sim%nt
    !             ! get the raw energy
    !             de(t,1)=sim%stat%potential_energy(t)
    !
    !             ! secondorder energy, first the polar term (if applicable)
    !             energy=0.0_flyt
    !             if ( fc%polar ) then
    !                 do a1=1,ss%na
    !                     v0=0.0_flyt
    !                     do a2=1,ss%na
    !                         v0=v0+matmul(forceconstant(:,:,a2,a1),u(:,a2,t))
    !                     enddo
    !                     energy=energy+dot_product(u(:,a1,t),v0)*0.5_flyt
    !                 enddo
    !             endif
    !             ! then the pair term, always there
    !             do a1=1,ss%na
    !                 v0=0.0_flyt
    !                 do i1=1,fcss%atom(a1)%n
    !                     a2=fcss%atom(a1)%pair(i1)%i2
    !                     m2=fcss%atom(a1)%pair(i1)%m
    !                     v0=v0+matmul(m2,u(:,a2,t))
    !                 enddo
    !                 energy=energy+dot_product(u(:,a1,t),v0)*0.5_flyt
    !             enddo
    !             de(t,2)=energy
    !
    !             ! triplet term
    !             f3=0.0_flyt
    !             energy=0.0_flyt
    !             do a1=1,fcss%na
    !                 v0=0.0_flyt
    !                 do i=1,fctss%atom(a1)%n
    !                     m3=fctss%atom(a1)%triplet(i)%m
    !                     a2=fctss%atom(a1)%triplet(i)%i2
    !                     a3=fctss%atom(a1)%triplet(i)%i3
    !                     u2=u(:,a2,t)
    !                     u3=u(:,a3,t)
    !                     do i1=1,3
    !                     do i2=1,3
    !                     do i3=1,3
    !                         v0(i1)=v0(i1)-m3(i1,i2,i3)*u2(i2)*u3(i3)
    !                     enddo
    !                     enddo
    !                     enddo
    !                 enddo
    !                 f3(:,a1)=v0*0.5_flyt
    !                 energy=energy-dot_product(f3(:,a1),u(:,a1,t))/3.0_flyt
    !             enddo
    !             de(t,3)=energy
    !
    !             ! quartet term
    !             f4=0.0_flyt
    !             energy=0.0_flyt
    !             do a1=1,fcfss%na
    !                 v0=0.0_flyt
    !                 do i=1,fcfss%atom(a1)%n
    !                     m4=fcfss%atom(a1)%quartet(i)%m
    !                     a2=fcfss%atom(a1)%quartet(i)%i2
    !                     a3=fcfss%atom(a1)%quartet(i)%i3
    !                     a4=fcfss%atom(a1)%quartet(i)%i4
    !                     u2=u(:,a2,t)
    !                     u3=u(:,a3,t)
    !                     u4=u(:,a4,t)
    !                     do i1=1,3
    !                     do i2=1,3
    !                     do i3=1,3
    !                     do i4=1,3
    !                         v0(i1)=v0(i1)-m4(i1,i2,i3,i4)*u2(i2)*u3(i3)*u4(i4)
    !                     enddo
    !                     enddo
    !                     enddo
    !                     enddo
    !                 enddo
    !                 f4(:,a1)=v0/6.0_flyt
    !                 energy=energy-dot_product(f4(:,a1),u(:,a1,t))/4.0_flyt
    !             enddo
    !             de(t,4)=energy
    !         enddo
    !     end block forceenergy
    !
    !     if ( mw%talk ) then
    !         write(*,*) '    raw avg potential energy:',lo_mean(de(:,1))
    !         write(*,*) '                        E-E2:',lo_mean(de(:,1)-de(:,2))
    !         write(*,*) '                     E-E2-E3:',lo_mean(de(:,1)-de(:,2)-de(:,3))
    !         write(*,*) '                  E-E2-E3-E4:',lo_mean(de(:,1)-de(:,2)-de(:,3)-de(:,4))
    !     endif
    !
    !
    ! !! Energy per atom
    ! !de=1000*de/sim%na
    ! !
    ! !if ( mw%talk ) then
    ! !write(*,*) 'RAW:',lo_mean(de(:,1)),lo_stddev(de(:,1))
    ! !write(*,*) 'RAW:',lo_mean(de(:,1)-de(:,3)),lo_stddev(de(:,1)-de(:,3))
    ! !write(*,*) 'RAW:',lo_mean(de(:,1)-de(:,2)-de(:,3)),lo_stddev(de(:,1)-de(:,2)-de(:,3))
    ! !write(*,*) 'RAW:',lo_mean(de(:,1)-de(:,2)-de(:,3)-de(:,4)),lo_stddev(de(:,1)-de(:,2)-de(:,3)-de(:,4))
    ! !write(*,*) 'RAW:',lo_mean(de(:,1)-de(:,2)-de(:,3)-de(:,4)-de(:,5)),lo_stddev(de(:,1)-de(:,2)-de(:,3)-de(:,4)-de(:,5))
    ! !u=open_file('out','outfile.U0')
    ! !    do t=1,sim%nt
    ! !        write(u,"(6(1X,E18.12))") de(t,:)/1000*sim%na,(de(t,1)-de(t,2)-de(t,3)-de(t,4)-de(t,5))/1000*sim%na
    ! !    enddo
    ! !close(u)
    ! !endif
    !
    ! !    t0=walltime()
    ! !    call lo_progressbar_init()
    ! !    ! Subtract forces
    ! !    lo_allocate(f(3,ss%na))
    ! !    do t=1,sim%nt
    ! !        ! calculate the long-range forces
    ! !        f=0.0_flyt
    ! !        do a1=1,ss%na
    ! !            do a2=1,ss%na
    ! !                v0=matmul(forceconstant(:,:,a1,a2),u(:,a2,t))
    ! !                f(:,a1)=f(:,a1)-v0
    ! !            enddo
    ! !        enddo
    ! !        ! sanity check that they add up to zero
    ! !        if ( abs(sum(f)) .gt. lo_sqtol ) then
    ! !            write(*,*) ''
    ! !            write(*,*) 'ERROR: dipole-dipole forces do not add up to zero!'
    ! !            stop
    ! !        endif
    ! !        ! subtract them
    ! !        sim%f(:,:,t)=sim%f(:,:,t)-f
    ! !        if( lo_trueNtimes(t,50,sim%nt) ) call lo_progressbar(' ... subtracting dipole-dipole forces',t,sim%nt)
    ! !    enddo
    ! !    call lo_progressbar(' ... subtracting dipole-dipole forces',sim%nt,sim%nt,walltime()-t0)
    ! end subroutine
    
    end module
    