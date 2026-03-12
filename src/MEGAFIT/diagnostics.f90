module diagnostics
!! Calculate R^2 and related diagnostics on the entire grid
use konstanter, only: r8,i8,lo_huge,lo_hugeint,lo_pi,lo_twopi,lo_imag,lo_status,lo_exitcode_param,lo_exitcode_symmetry,&
                      lo_tol,lo_sqtol,lo_pressure_HartreeBohr_to_GPa,lo_pressure_GPa_to_HartreeBohr,lo_Hartree_to_eV,&
                      lo_volume_bohr_to_A,lo_volume_A_to_bohr,lo_freqtol
use gottochblandat, only: walltime,tochar,lo_sqnorm,lo_mean,lo_rsquare,lo_stddev,open_file
use mpi_wrappers, only: lo_mpi_helper,lo_stop_gracefully,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_IN_PLACE
use lo_memtracker, only: lo_mem_helper
use type_crystalstructure, only: lo_crystalstructure
use type_forceconstant_firstorder, only: lo_forceconstant_firstorder
use type_forceconstant_secondorder, only: lo_forceconstant_secondorder
use type_forceconstant_thirdorder, only: lo_forceconstant_thirdorder
use type_forceconstant_fourthorder, only: lo_forceconstant_fourthorder
use lo_dielectric_interaction, only: lo_dielectric_tensor
use type_forcemap, only: lo_forcemap,lo_secondorder_rot_herm_huang
use type_mdsim, only: lo_mdsim
use type_gridsim, only: lo_gridsim
use type_polynomial_interpolation, only: lo_grid_interpolation

implicit none

private
public :: get_diagnostics

type lo_grid_diagnostics
    !> R^2 from fits
    real(r8), dimension(:), allocatable :: rsq_fp
    real(r8), dimension(:), allocatable :: rsq_f2
    real(r8), dimension(:), allocatable :: rsq_f3
    real(r8), dimension(:), allocatable :: rsq_f4
    real(r8), dimension(:), allocatable :: rsq_ep
    real(r8), dimension(:), allocatable :: rsq_e2
    real(r8), dimension(:), allocatable :: rsq_e3
    real(r8), dimension(:), allocatable :: rsq_e4
    real(r8), dimension(:), allocatable :: rsq_eps0
    real(r8), dimension(:), allocatable :: rsq_eps1
    real(r8), dimension(:), allocatable :: rsq_eps2

    !> residual standard deviation
    real(r8), dimension(:), allocatable :: dev_f0
    real(r8), dimension(:), allocatable :: dev_fp
    real(r8), dimension(:), allocatable :: dev_f2
    real(r8), dimension(:), allocatable :: dev_f3
    real(r8), dimension(:), allocatable :: dev_f4
    real(r8), dimension(:), allocatable :: dev_e0
    real(r8), dimension(:), allocatable :: dev_ep
    real(r8), dimension(:), allocatable :: dev_e2
    real(r8), dimension(:), allocatable :: dev_e3
    real(r8), dimension(:), allocatable :: dev_e4
    real(r8), dimension(:), allocatable :: dev_epso
    real(r8), dimension(:), allocatable :: dev_eps0
    real(r8), dimension(:), allocatable :: dev_eps1
    real(r8), dimension(:), allocatable :: dev_eps2

    ! Evaluated energies?
    real(r8), dimension(:), allocatable :: delta_U0
end type

contains

!> evaluate the irreducible representation at a certain point
subroutine get_diagnostics(gs,map,mw,mem,verbosity)
    !> simulation grid
    type(lo_gridsim), intent(inout) :: gs
    !> forcemap
    type(lo_forcemap), intent(inout) :: map
    !> MPI helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem
    !> talk a lot?
    integer, intent(in) :: verbosity

    type(lo_grid_diagnostics) :: diag
    character(len=5000), dimension(:), allocatable :: filenames

    ! Set some basic things
    init: block
        real(r8), dimension(gs%ndim) :: d1
        character(len=200), dimension(40) :: prtbuf
        character(len=1000) :: opf
        character(len=1) :: dum
        integer :: i,u,ctr

        if ( verbosity .gt. 0 ) then
            write(*,*) ''
            write(*,*) 'CALCULATING DIAGNOSTICS'
        endif

        allocate(filenames(gs%nsim))
        ! grab simulations again. Don't have the energy to sort again.
        u=open_file('in','infile.simulations')
            read(u,*) dum
            read(u,*) dum
            read(u,*) dum
            read(u,*) dum
            read(u,*) dum
            if ( allocated(gs%eos) ) then
                read(u,*) dum
            endif
            do i=1,gs%nsim
                read(u,*) d1,filenames(i)
            enddo
        close(u)

        allocate(diag%rsq_fp(gs%nsim))
        allocate(diag%rsq_f2(gs%nsim))
        allocate(diag%rsq_f3(gs%nsim))
        allocate(diag%rsq_f4(gs%nsim))
        allocate(diag%rsq_ep(gs%nsim))
        allocate(diag%rsq_e2(gs%nsim))
        allocate(diag%rsq_e3(gs%nsim))
        allocate(diag%rsq_e4(gs%nsim))
        allocate(diag%rsq_eps0(gs%nsim))
        allocate(diag%rsq_eps1(gs%nsim))
        allocate(diag%rsq_eps2(gs%nsim))
        diag%rsq_fp=0.0_r8
        diag%rsq_f2=0.0_r8
        diag%rsq_f3=0.0_r8
        diag%rsq_f4=0.0_r8
        diag%rsq_ep=0.0_r8
        diag%rsq_e2=0.0_r8
        diag%rsq_e3=0.0_r8
        diag%rsq_e4=0.0_r8
        diag%rsq_eps0=0.0_r8
        diag%rsq_eps1=0.0_r8
        diag%rsq_eps2=0.0_r8
        allocate(diag%dev_f0(gs%nsim))
        allocate(diag%dev_fp(gs%nsim))
        allocate(diag%dev_f2(gs%nsim))
        allocate(diag%dev_f3(gs%nsim))
        allocate(diag%dev_f4(gs%nsim))
        allocate(diag%dev_e0(gs%nsim))
        allocate(diag%dev_ep(gs%nsim))
        allocate(diag%dev_e2(gs%nsim))
        allocate(diag%dev_e3(gs%nsim))
        allocate(diag%dev_e4(gs%nsim))
        allocate(diag%dev_epso(gs%nsim))
        allocate(diag%dev_eps0(gs%nsim))
        allocate(diag%dev_eps1(gs%nsim))
        allocate(diag%dev_eps2(gs%nsim))
        diag%dev_f0=0.0_r8
        diag%dev_fp=0.0_r8
        diag%dev_f2=0.0_r8
        diag%dev_f3=0.0_r8
        diag%dev_f4=0.0_r8
        diag%dev_e0=0.0_r8
        diag%dev_ep=0.0_r8
        diag%dev_e2=0.0_r8
        diag%dev_e3=0.0_r8
        diag%dev_e4=0.0_r8
        diag%dev_epso=0.0_r8
        diag%dev_eps0=0.0_r8
        diag%dev_eps1=0.0_r8
        diag%dev_eps2=0.0_r8

        if ( mw%talk ) then
            ctr=0

            !ctr=ctr+1
            !prtbuf(ctr)='2nd'


            if ( map%polar .gt. 0 ) then
                ctr=ctr+1
                prtbuf(ctr)='polar'
            endif
            if ( map%have_fc_pair ) then
                ctr=ctr+1
                prtbuf(ctr)='2nd'
            endif
            if ( map%have_fc_triplet ) then
                ctr=ctr+1
                prtbuf(ctr)='3rd'
            endif
            if ( map%have_fc_quartet ) then
                ctr=ctr+1
                prtbuf(ctr)='4th'
            endif
            if ( map%polar .gt. 1 ) then
                ctr=ctr+1
                prtbuf(ctr)='eps0'
            endif
            if ( map%have_eps_singlet ) then
                ctr=ctr+1
                prtbuf(ctr)='eps1'
            endif
            if ( map%have_eps_pair ) then
                ctr=ctr+1
                prtbuf(ctr)='eps2'
            endif
            opf="("//tochar(gs%ndim*19)//"X,"//tochar(ctr)//"(4X,A10))"
            write(*,opf) prtbuf(1:ctr)
        endif

    end block init

    ! One way to do it is to read each simulation from file, fresh,
    ! and work from there? Maybe that makes sense.
    getforcesandenergies: block
        type(lo_mdsim) :: sim
        type(lo_crystalstructure) :: uc,ss
        type(lo_forceconstant_secondorder) :: fc,fc2_ss
        type(lo_forceconstant_thirdorder) :: fct,fc3_ss
        type(lo_forceconstant_fourthorder) :: fcf,fc4_ss
        type(lo_dielectric_tensor) :: di,diss
        real(r8), dimension(:,:,:,:), allocatable :: polar_fc
        real(r8), dimension(:,:), allocatable :: pairconstraints
        real(r8), dimension(:,:,:), allocatable :: f0,fp,f2,f3,f4
        real(r8), dimension(:,:,:), allocatable :: epso,eps0,eps1,eps2
        real(r8), dimension(:), allocatable :: e0,ep,e2,e3,e4,ebuf
        real(r8), dimension(3,3,3,3) :: m4
        real(r8), dimension(3,3,3) :: m3
        real(r8), dimension(3,3) :: m2,epsavg0,epsavgo
        real(r8), dimension(3) :: v0,u2,u3,u4
        real(r8) :: energy,baseline
        integer :: nconstr,isim,t,i1,i2,i3,i4
        integer :: a1,a2,a3,a4,i,ctr
        real(r8), dimension(gs%ndim) :: prtcoord
        real(r8), dimension(20) :: prtbuf,devbuf
        character(len=1000) :: opf

        do isim=1,gs%nsim
            ! Grab the simulation and related cells
            call sim%read_from_hdf5(trim(filenames(isim)),verbosity=-1,mw=mw)
            call uc%generate(sim%extra%unitcell_latticevectors, sim%extra%unitcell_positions,sim%extra%unitcell_atomic_numbers, 2 )
            call ss%generate(sim%extra%supercell_latticevectors,sim%extra%supercell_positions,sim%extra%supercell_atomic_numbers, enhet=2 )
            call ss%classify('supercell',uc)
            !@TODO think about what happens with reference positions. Later problem for now.

            ! forceconstants, first the constraints
            !call lo_secondorder_rot_herm_huang( map,uc,pairconstraints,nconstr,.true.,.true.,.true. )
nconstr=0
            if ( nconstr .gt. 0 ) then
                call gs%eval(map,gs%grid_coordinates(:,isim),pairconstraints)
            else
                call gs%eval(map,gs%grid_coordinates(:,isim))
            endif
            if ( map%have_fc_pair ) then
                call map%get_secondorder_forceconstant(uc,fc,mem,-1)
                call fc%remap(uc,ss,fc2_ss)
            endif
            if ( map%polar .gt. 0 .and. map%xuc%nx_Z_singlet .gt. 0 ) then
                call mem%allocate(polar_fc,[3,3,ss%na,ss%na],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
                polar_fc=0.0_r8
                call fc%supercell_longrange_dynamical_matrix_at_gamma(ss,polar_fc,1E-12_r8)
            else
                call mem%allocate(polar_fc,[1,1,1,1],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
                polar_fc=-lo_huge
            endif
            if ( map%have_fc_triplet ) then
                call map%get_thirdorder_forceconstant(uc,fct)
                call fct%remap(uc,ss,fc3_ss)
            endif
            if ( map%have_fc_quartet ) then
                call map%get_fourthorder_forceconstant(uc,fcf)
                call fcf%remap(uc,ss,fc4_ss)
            endif
            if ( map%polar .gt. 1 ) then
                call map%get_dielectric_tensors(uc,di)
                call di%remap(diss,uc,ss,mw,mem,-1)
            endif

            ! Space for forces and energies
            call mem%allocate(f0,[3,ss%na,sim%nt],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
            call mem%allocate(fp,[3,ss%na,sim%nt],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
            call mem%allocate(f2,[3,ss%na,sim%nt],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
            call mem%allocate(f3,[3,ss%na,sim%nt],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
            call mem%allocate(f4,[3,ss%na,sim%nt],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
            call mem%allocate(e0,sim%nt,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
            call mem%allocate(ep,sim%nt,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
            call mem%allocate(e2,sim%nt,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
            call mem%allocate(e3,sim%nt,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
            call mem%allocate(e4,sim%nt,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
            call mem%allocate(ebuf,sim%nt,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
            call mem%allocate(epso,[3,3,sim%nt],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
            call mem%allocate(eps0,[3,3,sim%nt],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
            call mem%allocate(eps1,[3,3,sim%nt],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
            call mem%allocate(eps2,[3,3,sim%nt],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)

            f0=0.0_r8
            fp=0.0_r8
            f2=0.0_r8
            f3=0.0_r8
            f4=0.0_r8
            e0=0.0_r8
            ep=0.0_r8
            e2=0.0_r8
            e3=0.0_r8
            e4=0.0_r8
            ebuf=0.0_r8
            epso=0.0_r8
            eps0=0.0_r8
            eps1=0.0_r8
            eps2=0.0_r8

            ! Start calculating forces and energies
            ! Calculate energies and stuff
            do t=1,sim%nt
                ! make it parallel to not confuse anyone
                if ( mod(t,mw%n) .ne. mw%r ) cycle
                ! Copy of DFT force and energy
                e0(t)=sim%stat%potential_energy(t)
                f0(:,:,t)=sim%f(:,:,t)
                ! then the pair term
                if ( map%have_fc_pair ) then
                    energy=0.0_r8
                    do a1=1,ss%na
                        v0=0.0_r8
                        do i1=1,fc2_ss%atom(a1)%n
                            a2=fc2_ss%atom(a1)%pair(i1)%i2
                            m2=fc2_ss%atom(a1)%pair(i1)%m
                            v0=v0-matmul(m2,sim%u(:,a2,t))
                        enddo
                        energy=energy-dot_product(sim%u(:,a1,t),v0)*0.5_r8
                        f2(:,a1,t)=v0
                    enddo
                    e2(t)=energy
                endif
                ! Possible polar term?
                if ( map%polar .gt. 0 .and. map%xuc%nx_Z_singlet .gt. 0 ) then
                    energy=0.0_r8
                    do a1=1,ss%na
                        v0=0.0_r8
                        do a2=1,ss%na
                            v0=v0-matmul(polar_fc(:,:,a1,a2),sim%u(:,a2,t))
                        enddo
                        energy=energy-dot_product(sim%u(:,a1,t),v0)*0.5_r8
                        fp(:,a1,t)=v0
                    enddo
                    ep(t)=energy
                endif
                ! triplet term
                if ( map%have_fc_triplet ) then
                    energy=0.0_r8
                    do a1=1,fc3_ss%na
                        v0=0.0_r8
                        do i=1,fc3_ss%atom(a1)%n
                            m3=fc3_ss%atom(a1)%triplet(i)%m
                            a2=fc3_ss%atom(a1)%triplet(i)%i2
                            a3=fc3_ss%atom(a1)%triplet(i)%i3
                            u2=sim%u(:,a2,t)
                            u3=sim%u(:,a3,t)
                            do i1=1,3
                            do i2=1,3
                            do i3=1,3
                                v0(i1)=v0(i1)-m3(i1,i2,i3)*u2(i2)*u3(i3)
                            enddo
                            enddo
                            enddo
                        enddo
                        v0=v0*0.5_r8
                        f3(:,a1,t)=v0
                        energy=energy-dot_product(v0,sim%u(:,a1,t))/3.0_r8
                    enddo
                    e3(t)=energy
                endif
                ! quartet term
                if ( map%have_fc_quartet ) then
                    energy=0.0_r8
                    do a1=1,fc4_ss%na
                        v0=0.0_r8
                        do i=1,fc4_ss%atom(a1)%n
                            m4=fc4_ss%atom(a1)%quartet(i)%m
                            a2=fc4_ss%atom(a1)%quartet(i)%i2
                            a3=fc4_ss%atom(a1)%quartet(i)%i3
                            a4=fc4_ss%atom(a1)%quartet(i)%i4
                            u2=sim%u(:,a2,t)
                            u3=sim%u(:,a3,t)
                            u4=sim%u(:,a4,t)
                            do i1=1,3
                            do i2=1,3
                            do i3=1,3
                            do i4=1,3
                                v0(i1)=v0(i1)-m4(i1,i2,i3,i4)*u2(i2)*u3(i3)*u4(i4)
                            enddo
                            enddo
                            enddo
                            enddo
                        enddo
                        v0=v0/6.0_r8
                        f4(:,a1,t)=v0
                        energy=energy-dot_product(v0,sim%u(:,a1,t))/4.0_r8
                    enddo
                    e4(t)=energy
                endif

                ! Also possibly the dielectric stuff?
                if ( sim%have_dielectric ) then
                    ! raw eps
                    epso(:,:,t)=sim%eps(:,:,t)
                    ! average eps
                    eps0(:,:,t)=diss%eps_inf

                    if ( map%have_eps_singlet ) then
                        m2=0.0_r8
                        do a1=1,diss%n_eps_singlet
                            do i1=1,3
                                m2=m2+diss%eps_singlet(a1)%m(:,:,i1)*sim%u(i1,a1,t)
                            enddo
                        enddo
                        ! convert to dielectric constant?
                        m2=-m2*4*lo_pi/ss%volume
                        ! do i=1,3
                        !     m2(i,i)=1.0_r8+m2(i,i)
                        ! enddo
                        eps1(:,:,t)=m2
                    endif

                    if ( map%have_eps_pair ) then
                        m2=0.0_r8
                        do i=1,diss%n_eps_pair
                            a1=diss%eps_pair(i)%a1
                            a2=diss%eps_pair(i)%a2
                            do i1=1,3
                            do i2=1,3
                                m2=m2+diss%eps_pair(i)%m(:,:,i1,i2)*sim%u(i1,a1,t)*sim%u(i2,a2,t)
                            enddo
                            enddo
                        enddo
                        ! convert to dielectric constant?
                        m2=-m2*4*lo_pi/ss%volume
                        ! do i=1,3
                        !     m2(i,i)=1.0_r8+m2(i,i)
                        ! enddo
                        eps2(:,:,t)=m2
                    endif
                endif
            enddo
            ! sync across ranks
            call mw%allreduce('sum',f0)
            call mw%allreduce('sum',fp)
            call mw%allreduce('sum',f2)
            call mw%allreduce('sum',f3)
            call mw%allreduce('sum',f4)
            call mw%allreduce('sum',e0)
            call mw%allreduce('sum',ep)
            call mw%allreduce('sum',e2)
            call mw%allreduce('sum',e3)
            call mw%allreduce('sum',e4)

            call mw%allreduce('sum',epso)
            call mw%allreduce('sum',eps0)
            call mw%allreduce('sum',eps1)
            call mw%allreduce('sum',eps2)

            ! Calculate the actual diagnostics
            ebuf=e0-e2-e3-e4-ep
            baseline=lo_mean(ebuf)
            e0=e0-baseline

            ! Calculate some R^2 thingies
            diag%rsq_fp(isim)=lo_rsquare(f0,fp)
            diag%rsq_f2(isim)=lo_rsquare(f0,fp+f2)
            diag%rsq_f3(isim)=lo_rsquare(f0,fp+f2+f3)
            diag%rsq_f4(isim)=lo_rsquare(f0,fp+f2+f3+f4)
            diag%rsq_ep(isim)=lo_rsquare(e0,ep)
            diag%rsq_e2(isim)=lo_rsquare(e0,ep+e2)
            diag%rsq_e3(isim)=lo_rsquare(e0,ep+e2+e3)
            diag%rsq_e4(isim)=lo_rsquare(e0,ep+e2+e3+e4)
            diag%rsq_eps0(isim)=lo_rsquare(epso,eps0)
            diag%rsq_eps1(isim)=lo_rsquare(epso,eps0+eps1)
            diag%rsq_eps2(isim)=lo_rsquare(epso,eps0+eps1+eps2)

            diag%dev_fp(isim)=lo_stddev(f0)
            diag%dev_fp(isim)=lo_stddev(f0-fp)
            diag%dev_f2(isim)=lo_stddev(f0-fp-f2)
            diag%dev_f3(isim)=lo_stddev(f0-fp-f2-f3)
            diag%dev_f4(isim)=lo_stddev(f0-fp-f2-f3-f4)
            diag%dev_ep(isim)=lo_stddev(e0)
            diag%dev_ep(isim)=lo_stddev(e0-ep)
            diag%dev_e2(isim)=lo_stddev(e0-ep-e2)
            diag%dev_e3(isim)=lo_stddev(e0-ep-e2-e3)
            diag%dev_e4(isim)=lo_stddev(e0-ep-e2-e3-e4)
            m2=lo_stddev_33matrix(epso)
            diag%dev_epso(isim)=sum(m2)
            m2=lo_stddev_33matrix(epso-eps0)
            diag%dev_eps0(isim)=sum(m2)
            m2=lo_stddev_33matrix(epso-eps0-eps1)
            diag%dev_eps1(isim)=sum(m2)
            m2=lo_stddev_33matrix(epso-eps0-eps1-eps2)
            diag%dev_eps2(isim)=sum(m2)

            if ( mw%talk ) then
                ctr=0
                if ( map%polar .gt. 0 ) then
                    ctr=ctr+1
                    prtbuf(ctr)=diag%rsq_fp(isim)
                endif
                if ( map%have_fc_pair ) then
                    ctr=ctr+1
                    prtbuf(ctr)=diag%rsq_f2(isim)
                endif
                if ( map%have_fc_triplet ) then
                    ctr=ctr+1
                    prtbuf(ctr)=diag%rsq_f3(isim)
                endif
                if ( map%have_fc_quartet ) then
                    ctr=ctr+1
                    prtbuf(ctr)=diag%rsq_f4(isim)
                endif
                if ( map%have_eps_singlet .or. map%have_eps_pair ) then
                    ctr=ctr+1
                    prtbuf(ctr)=diag%rsq_eps0(isim)
                endif
                if ( map%have_eps_singlet ) then
                    ctr=ctr+1
                    prtbuf(ctr)=diag%rsq_eps1(isim)
                endif
                if ( map%have_eps_pair ) then
                    ctr=ctr+1
                    prtbuf(ctr)=diag%rsq_eps2(isim)
                endif

                prtcoord=gs%grid_coordinates(:,isim)
                if ( gs%info%dim_volume .gt. 0 ) then
                    prtcoord(gs%info%dim_volume)=prtcoord(gs%info%dim_volume)*lo_volume_bohr_to_A
                endif

                opf="(1X,"//tochar(gs%ndim)//"(1X,F18.11),3X,"//tochar(ctr)//"(1X,F10.6))"
                write(*,opf) prtcoord,prtbuf(1:ctr)

                ! m2=0.0_r8
                ! do i=1,sim%nt
                !     m2=m2+epso(:,:,i)/sim%nt
                ! enddo
                ! do i=1,3
                !     write(*,*) 'epso',m2(:,i)
                ! enddo
                !
                ! m2=0.0_r8
                ! do i=1,sim%nt
                !     m2=m2+eps0(:,:,i)/sim%nt
                ! enddo
                ! do i=1,3
                !     write(*,*) 'eps0',m2(:,i)
                ! enddo
                !
                ! m2=0.0_r8
                ! do i=1,sim%nt
                !     m2=m2+eps1(:,:,i)/sim%nt
                ! enddo
                ! do i=1,3
                !     write(*,*) 'eps1',m2(:,i)
                ! enddo
                !
                ! m2=0.0_r8
                ! do i=1,sim%nt
                !     m2=m2+eps2(:,:,i)/sim%nt
                ! enddo
                ! do i=1,3
                !     write(*,*) 'eps2',m2(:,i)
                ! enddo
            endif

            ! Some cleanup
            call mem%deallocate(polar_fc,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
            call mem%deallocate(f0,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
            call mem%deallocate(fp,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
            call mem%deallocate(f2,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
            call mem%deallocate(f3,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
            call mem%deallocate(f4,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
            call mem%deallocate(e0,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
            call mem%deallocate(ep,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
            call mem%deallocate(e2,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
            call mem%deallocate(e3,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
            call mem%deallocate(e4,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
            call mem%deallocate(ebuf,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
            call mem%deallocate(epso,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
            call mem%deallocate(eps0,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
            call mem%deallocate(eps1,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
            call mem%deallocate(eps2,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        enddo

    end block getforcesandenergies

    report: block
        real(r8), dimension(20) :: prtbuf
        integer :: isim,ctr
        !character(len=1000), dimension(20) :: opf
        character(len=1000) :: opf

        do isim=1,gs%nsim
        if ( mw%talk ) then
            ctr=0

            ctr=ctr+1
            prtbuf(ctr)=diag%dev_f0(isim)

            if ( map%polar .gt. 0 ) then
                ctr=ctr+1
                prtbuf(ctr)=diag%dev_fp(isim)
            endif
            if ( map%have_fc_pair ) then
                ctr=ctr+1
                prtbuf(ctr)=diag%dev_f2(isim)
            endif
            if ( map%have_fc_triplet ) then
                ctr=ctr+1
                prtbuf(ctr)=diag%dev_f3(isim)
            endif
            if ( map%have_fc_quartet ) then
                ctr=ctr+1
                prtbuf(ctr)=diag%dev_f4(isim)
            endif
            if ( map%polar .gt. 1 ) then
                ctr=ctr+1
                prtbuf(ctr)=diag%dev_epso(isim)
                ctr=ctr+1
                prtbuf(ctr)=diag%dev_eps0(isim)
            endif
            if ( map%have_eps_singlet ) then
                ctr=ctr+1
                prtbuf(ctr)=diag%dev_eps1(isim)
            endif
            if ( map%have_eps_pair ) then
                ctr=ctr+1
                prtbuf(ctr)=diag%dev_eps2(isim)
            endif
            opf="("//tochar(ctr)//"(1X,F10.6))"
            write(*,opf) prtbuf(1:ctr)
        endif
        enddo

    end block report

end subroutine

!> standard deviation for 3x3-matrices
function lo_stddev_33matrix(x) result(s)
    !> values
    real(r8), dimension(:,:,:), intent(in) :: x
    !> standard deviation
    real(r8), dimension(3,3) :: s

    real(r8), dimension(3,3) :: mean
    integer :: i

    mean=0.0_r8
    do i=1,size(x,3)
        mean=mean+x(:,:,i)
    enddo
    mean=mean/real(size(x,3),r8)

    s=0.0_r8
    do i=1,size(x,3)
        s=s+(mean-x(:,:,i))**2
    enddo
    s=sqrt(s/real(size(x,3),r8))
end function

end module
