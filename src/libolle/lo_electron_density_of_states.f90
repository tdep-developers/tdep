module lo_electron_density_of_states
!!
!! Calculate the electronic density of states
!!
use konstanter, only: r8,i8,lo_iou,lo_tol,lo_sqtol,lo_huge,lo_hugeint,lo_exitcode_param,lo_eV_to_Hartree,lo_Hartree_to_eV,lo_gnuplotterminal,lo_twopi
use lo_memtracker, only: lo_mem_helper
use gottochblandat, only: open_file,lo_gauss,lo_trueNtimes,lo_trapezoid_integration,tochar,walltime,lo_fermi,&
                          lo_chop,lo_progressbar_init,lo_progressbar,lo_linspace,lo_untangle_one_tetrahedron,&
                          lo_clean_fractional_coordinates
use mpi_wrappers, only: lo_mpi_helper,lo_stop_gracefully
use hdf5_wrappers, only: lo_hdf5_helper
use type_crystalstructure, only: lo_crystalstructure
use lo_electron_dispersion_relations, only: lo_electron_dispersions
use type_qpointmesh, only: lo_qpoint_mesh,lo_integration_weights_for_one_tetrahedron,lo_LV_tetrahedron_weights

implicit none
private
public :: lo_electron_dos

!> Electron density of states
type lo_electron_dos
    !> energy axis
    real(r8), dimension(:), allocatable :: energy
    !> density of states
    real(r8), dimension(:,:), allocatable :: dos
    !> band projected density of states
    real(r8), dimension(:,:,:), allocatable :: pdos_band

    !> number of dos-points
    integer :: n_dos_point=-lo_hugeint
    !> number atoms
    integer :: n_atom=-lo_hugeint
    !> number bands
    integer :: n_band=-lo_hugeint
    !> number of spin channels
    integer :: n_spin=-lo_hugeint
    !> number of electrons
    real(r8) :: n_electron=-lo_huge

    !> smearing
    real(r8) :: dossmear=-lo_huge
    !> maximum energy
    real(r8) :: dosmax=-lo_huge
    !> minimum energy
    real(r8) :: dosmin=-lo_huge
    !> Fermi energy
    real(r8) :: efermi=-lo_huge
    !> band gap
    real(r8) :: bandgap=-lo_huge

    !> how many bands are far below the fermi level?
    integer :: n_below=-lo_hugeint
    integer, dimension(:), allocatable :: below
    !> how many are far above the conduction band?
    integer :: n_above=-lo_hugeint
    integer, dimension(:), allocatable :: above
    !> how many have to be treated carefully?
    integer :: n_near=-lo_hugeint
    integer, dimension(:), allocatable :: near
    !> how many bands are really really far below the fermi level?
    integer :: n_core=-lo_hugeint
    integer, dimension(:), allocatable :: core
    !> how many bands are really really far above the fermi level
    integer :: n_sky=-lo_hugeint
    integer, dimension(:), allocatable :: sky

    contains
        !> Calculate the DOS
        procedure :: generate
        !> dump it to file
        procedure :: write_to_file
        !> dump to hdf5
        procedure :: write_to_hdf5
        !> destroy?
        procedure :: destroy
        !> measure size in memory
        procedure :: size_in_mem
end type

contains

!> Calculate the electronic density of states
subroutine generate(ed,edr,kp,p,integrationtype,smearing,mw,mem,verbosity,ndos)
    !> electronic dos
    class(lo_electron_dos), intent(out) :: ed
    !> dispersion relations
    type(lo_electron_dispersions), intent(in) :: edr
    !> k-point grid
    class(lo_qpoint_mesh), intent(inout) :: kp
    !> crystalstructure
    type(lo_crystalstructure), intent(in) :: p
    !> integrationtype
    integer, intent(in) :: integrationtype
    !> how much to smear, of applicable
    real(r8), intent(in) :: smearing
    !> mpi helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem
    !> how much to talk
    integer, intent(in) :: verbosity
    !> change the default number of points on the x-axis
    integer, intent(in), optional :: ndos

    ! Some magic parameters, first, what is considered close in energy? Should be on the safe side, for sure.
    real(r8), parameter :: efilter=4.0_r8*lo_eV_to_Hartree
    real(r8) :: timer,t0,t1

    ! Start timers
    timer=walltime()
    t0=timer
    t1=timer

    ! Set some basics
    init: block

        if ( verbosity .gt. 0 ) then
            write(lo_iou,*) ''
            write(lo_iou,*) 'Calculating electronic density of states over ',tochar(mw%n),' MPI ranks'
        endif

        ! number of dospoints
        if ( present(ndos) ) then
            ed%n_dos_point=ndos
        else
            ed%n_dos_point=10000
        endif

        ! Store other sensible information
        ed%n_atom=p%na
        ed%n_band=edr%n_band
        ed%n_spin=edr%n_spin
        ed%bandgap=edr%bandgap
        ed%efermi=edr%efermi
    end block init

    ! Filter bands slightly
    filterbands: block
        integer :: i
        ! filter bands, since they are supposed to be treated slightly different!
        ed%n_core=0
        ed%n_sky=0
        ed%n_below=0
        ed%n_above=0
        ed%n_near=0
        do i=1,ed%n_band
            if ( edr%bandmax(i)-edr%evalence .lt. 0.0_r8 ) then
                ed%n_core=ed%n_core+1
            elseif ( edr%bandmin(i)-edr%efermi-edr%bandgap .gt. efilter*2 ) then
                ed%n_sky=ed%n_sky+1
            elseif ( edr%bandmax(i)-edr%efermi .lt. -efilter ) then
                ed%n_below=ed%n_below+1
            elseif ( edr%bandmin(i)-edr%efermi-edr%bandgap .gt. efilter ) then
                ed%n_above=ed%n_above+1
            else
                ed%n_near=ed%n_near+1
            endif
        enddo
        if ( ed%n_core  .gt. 0 ) allocate(ed%core( ed%n_core ))
        if ( ed%n_sky   .gt. 0 ) allocate(ed%sky(  ed%n_sky  ))
        if ( ed%n_below .gt. 0 ) allocate(ed%below(ed%n_below))
        if ( ed%n_above .gt. 0 ) allocate(ed%above(ed%n_above))
        if ( ed%n_near  .gt. 0 ) allocate(ed%near( ed%n_near ))
        ed%n_core=0
        ed%n_sky=0
        ed%n_below=0
        ed%n_above=0
        ed%n_near=0
        do i=1,ed%n_band
            if ( edr%bandmax(i)-edr%evalence .lt. 0.0_r8 ) then
                ed%n_core=ed%n_core+1
                ed%core(ed%n_core)=i
            elseif ( edr%bandmin(i)-edr%efermi-edr%bandgap .gt. efilter*2 ) then
                ed%n_sky=ed%n_sky+1
                ed%sky(ed%n_sky)=i
            elseif ( edr%bandmax(i)-edr%efermi .lt. -efilter ) then
                ed%n_below=ed%n_below+1
                ed%below(ed%n_below)=i
            elseif ( edr%bandmin(i)-edr%efermi-edr%bandgap .gt. efilter ) then
                ed%n_above=ed%n_above+1
                ed%above(ed%n_above)=i
            else
                ed%n_near=ed%n_near+1
                ed%near(ed%n_near)=i
            endif
        enddo

        ! Adjust the number of electrons in the DOS
        ed%n_electron=edr%n_electron-ed%n_core*2.0_r8
    end block filterbands

    ! Get the relevant energy ranges and guess smearing parameters
    rangesandsmearing: block
        integer :: i,j
        ! Smearing parameter?
        ed%dossmear=0.0_r8
        if ( ed%n_below .gt. 0 ) ed%dossmear=max(ed%dossmear,maxval(edr%default_smearing(ed%below,:)) )
        if ( ed%n_near  .gt. 0 ) ed%dossmear=max(ed%dossmear,maxval(edr%default_smearing(ed%near, :)) )
        if ( ed%n_above .gt. 0 ) ed%dossmear=max(ed%dossmear,maxval(edr%default_smearing(ed%above,:)) )

        ! Get the energy range
        ed%dosmin=lo_huge
        ed%dosmax=-lo_huge
        do j=1,ed%n_below
            i=ed%below(j)
            ed%dosmin=min(ed%dosmin,edr%bandmin(i))
            ed%dosmax=max(ed%dosmax,edr%bandmax(i))
        enddo
        do j=1,ed%n_near
            i=ed%near(j)
            ed%dosmin=min(ed%dosmin,edr%bandmin(i))
            ed%dosmax=max(ed%dosmax,edr%bandmax(i))
        enddo
        do j=1,ed%n_above
            i=ed%above(j)
            ed%dosmin=min(ed%dosmin,edr%bandmin(i))
            ed%dosmax=max(ed%dosmax,edr%bandmax(i))
        enddo
        ed%dosmin=ed%dosmin-ed%dossmear*10.0_r8
        ed%dosmax=ed%dosmax+ed%dossmear*10.0_r8

        ! use these to get energy axis
        allocate(ed%energy(ed%n_dos_point))
        call lo_linspace(ed%dosmin,ed%dosmax,ed%energy)
        ! and to make space for the actual dos
        allocate(ed%dos(ed%n_dos_point,ed%n_spin))
        allocate(ed%pdos_band(ed%n_dos_point,ed%n_band,ed%n_spin))
        ed%dos=0.0_r8
        ed%pdos_band=0.0_r8
    end block rangesandsmearing

    ! Intermediate report
    if ( verbosity .gt. 0 ) then
        write(lo_iou,*) '... got smearing parameter: ',tochar(ed%dossmear*lo_Hartree_to_eV),' eV'
        write(lo_iou,*) '... starting with rough Fermi level: ',tochar(edr%efermi*lo_Hartree_to_eV)
        if ( ed%n_core  .gt. 0 ) write(lo_iou,*) '... ',tochar(ed%n_core), ' bands designated core'
        if ( ed%n_below .gt. 0 ) write(lo_iou,*) '... ',tochar(ed%n_below),' bands far below eF'
        if ( ed%n_near  .gt. 0 ) write(lo_iou,*) '... ',tochar(ed%n_near), ' bands near eF'
        if ( ed%n_above .gt. 0 ) write(lo_iou,*) '... ',tochar(ed%n_above),' bands far above eF'
        if ( ed%n_sky   .gt. 0 ) write(lo_iou,*) '... ',tochar(ed%n_sky),  ' bands up in the sky'
    endif

    ! Perform operation actual operation
    select case(integrationtype)
    case(1:2)
        ! This is the default hybrid gaussian/tetrahedron thing.
        if ( verbosity .gt. 0 ) then
            write(lo_iou,*) '... hybrid tetrahedron/gaussian integration'
        endif
        call get_electron_dos_hybrid(ed,p,kp,edr,integrationtype,mw,mem,verbosity)
    case(3:6)
        ! Gaussian family of integrations?
        if ( verbosity .gt. 0 ) then
            write(lo_iou,*) '... hybrid tetrahedron/gaussian integration'
        endif
        call get_electron_dos_smearing(ed,p,kp,edr,integrationtype,smearing,mw,mem,verbosity)
    case default
        call lo_stop_gracefully(['Unknown integration type'],lo_exitcode_param,__FILE__,__LINE__,mw%comm)
    end select

    ! Make sure the DOS is clean and neat.
    normalizeandclean: block
        real(r8) :: f0
        integer :: ispin,iband

        ! Sum up contributions and normalize things
        ed%dos=0.0_r8
        do ispin=1,ed%n_spin
        do iband=1,ed%n_band
            f0=lo_trapezoid_integration(ed%energy,ed%pdos_band(:,iband,ispin))
            ed%dos(:,ispin)=ed%dos(:,ispin)+lo_chop(ed%pdos_band(:,iband,ispin),lo_sqtol)
            ed%pdos_band(:,iband,ispin)=ed%pdos_band(:,iband,ispin)/f0
        enddo
        enddo

        ! ! Add a tiny tiny smearing to make things just a little bit smoother.
        ! lo_allocate(dum(ed%ndos))
        ! f0=deltae*6
        ! emin=minval(ed%energy)
        ! erange=1.0_r8/(maxval(ed%energy)-minval(ed%energy))
        ! do k=1,ed%nspin
        !     dum=0.0_r8
        !     do i=1,ed%ndos
        !         e=ed%energy(i)
        !         ii=floor( ed%ndos*(e-4*f0-emin)*erange )
        !         ii=max(ii,1)
        !         ii=min(ii,ed%ndos)
        !         jj=ceiling( ed%ndos*(e+4*f0-emin)*erange )
        !         jj=max(jj,1)
        !         jj=min(jj,ed%ndos)
        !         do j=ii,jj
        !             dum(i)=dum(i)+ed%dos(j,k)*lo_gauss( ed%energy(i),ed%energy(j),f0 )
        !         enddo
        !     enddo
        !     ed%dos(:,k)=dum
        ! enddo

        ! Normalize the total
        f0=0.0_r8
        do ispin=1,ed%n_spin
            f0=f0+lo_trapezoid_integration(ed%energy,ed%dos(:,ispin))
        enddo

        ! not sure about factor 2 here
        select case(ed%n_spin)
        case(1)
            ed%dos=ed%dos*2.0_r8*(ed%n_below+ed%n_near+ed%n_above)/f0
        case(2)
            ed%dos=ed%dos*(ed%n_below+ed%n_near+ed%n_above)/f0
        end select
    end block normalizeandclean

    if ( verbosity .gt. 0 ) then
        write(lo_iou,*) '... eF: ',tochar(ed%efermi*lo_Hartree_to_eV),' gap: ',tochar(ed%bandgap*lo_Hartree_to_eV)
        write(lo_iou,*) 'Done integration density of states ('//tochar(walltime()-timer)//'s)'
    endif
end subroutine

!> Integrate the electron dos with the tetrahedron method
subroutine get_electron_dos_smearing(ed,p,kp,edr,integrationtype,smearing_adjustment,mw,mem,verbosity)
    !> electronic dos
    type(lo_electron_dos), intent(inout) :: ed
    !> structure
    type(lo_crystalstructure), intent(in) :: p
    !> k-point grid
    class(lo_qpoint_mesh), intent(in) :: kp
    !> dispersion relations
    type(lo_electron_dispersions), intent(in) :: edr
    !> how to integrate
    integer, intent(in) :: integrationtype
    !> adjustment for the smearing parameter
    real(r8), intent(in) :: smearing_adjustment
    !> mpi helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem
    !> talk a lot?
    integer, intent(in) :: verbosity

    ! Tetrahedron energies
    real(r8) :: t0
    integer, dimension(:), allocatable :: relevant_bands
    integer :: n_band,n_kp

    ! Start integration timer
    t0=walltime()

    ! Set basic things
    init: block
        integer :: i,j

        ! Select the relevant bands
        n_band=ed%n_below+ed%n_near+ed%n_above
        call mem%allocate(relevant_bands,n_band,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        j=0
        do i=1,ed%n_below
            j=j+1
            relevant_bands(j)=ed%below(i)
        enddo
        do i=1,ed%n_near
            j=j+1
            relevant_bands(j)=ed%near(i)
        enddo
        do i=1,ed%n_above
            j=j+1
            relevant_bands(j)=ed%above(i)
        enddo

        ! Also count number of k-points per rank
        n_kp=0
        do i=1,kp%n_irr_point
            if ( mod(i,mw%n) .eq. mw%r ) n_kp=n_kp+1
        enddo
    end block init

    ! Do the actual integration.
    actualintegration: block
        real(r8) :: invf,avgsigma,totmine
        integer :: i,ispin,iband,ctr,ikp

        ! some misc things to make it a little faster
        ed%pdos_band=0.0_r8
        invf=real(ed%n_dos_point,r8)/(maxval(ed%energy)-minval(ed%energy))
        totmine=minval(ed%energy)

        ! the average smearing for the relevant bands
        avgsigma=sum(edr%default_smearing(relevant_bands,:))/size(relevant_bands)/edr%n_spin

        t0=walltime()

        ! Also gaussian integrate the DOS for the bands far from the Fermi level.
        if ( verbosity .gt. 0 ) call lo_progressbar_init()
        ctr=0
        do ispin=1,ed%n_spin
        do ikp=1,kp%n_irr_point
            if ( mod(ikp,mw%n) .ne. mw%r ) cycle
            ctr=ctr+1
            do i=1,n_band
                iband=relevant_bands(i)
                ! Now it will depend...
                select case(integrationtype)
                case(3:5)
                gaussfamily: block
                    real(r8) :: f0,f1,sigma,fivesigma
                    integer :: ie,ii,jj

                    ! First decide on a smearing
                    select case(integrationtype)
                    case(3)
                        ! Fix Gaussian smearing
                        sigma=avgsigma*smearing_adjustment
                    case(4)
                        ! Semi-adaptive Gaussian smearing
                        sigma=edr%default_smearing(iband,ispin)*smearing_adjustment*2
                    case(5)
                        ! Even more adaptive gaussian smearing thing
                        f0=norm2(edr%ik(ikp)%groupvelocity(:,iband,ispin))
                        f0=f0*kp%ip(ikp)%radius
                        sigma=f0*smearing_adjustment
                        ! Don't make the smearing parameter too agressive.
                        sigma=max(sigma,avgsigma/10)
                        sigma=min(sigma,avgsigma*10)
                    end select
                    fivesigma=5*sigma

                    ! Then add to the histogram
                    ! Add to histogram
                    f0=edr%ik(ikp)%eigenvalue(iband,ispin)
                    ii=max(floor( (f0-totmine-fivesigma)*invf ),1)
                    jj=min(ceiling( (f0-totmine+fivesigma)*invf ),ed%n_dos_point)
                    do ie=ii,jj
                        f1=lo_gauss(ed%energy(ie),f0,sigma)*kp%ip(ikp)%integration_weight
                        ed%pdos_band(ie,iband,ispin)=ed%pdos_band(ie,iband,ispin)+f1
                    enddo
                end block gaussfamily
                case(6)
                spherefamily: block
                    real(r8) :: f0,f1,grad,igrad,r,bigrsq,rsq,deltae,prefactor
                    integer :: ie,ii,jj

                    ! Linear sphere thing, not sure if meaningful.
                    grad=norm2(edr%ik(ikp)%groupvelocity(:,iband,ispin))
                    ! insert tolerance here. If no gradient, revert to adaptive Gaussian.
                    if ( grad .lt. 1E-10_r8 ) cycle

                    igrad=1.0_r8/grad
                    r=kp%ip(ikp)%radius
                    deltae=grad*r
                    bigrsq=r**2
                    prefactor=3.0_r8*igrad/(2*r**3)

                    f0=edr%ik(ikp)%eigenvalue(iband,ispin)
                    ii=max(floor( (f0-totmine-deltae)*invf ),1)
                    jj=min(ceiling( (f0-totmine+deltae)*invf ),ed%n_dos_point)
                    do ie=ii,jj
                        rsq=bigrsq-((ed%energy(ie)-f0)*igrad)**2
                        if ( rsq .lt. 0.0_r8 ) cycle
                        f1=prefactor*rsq*kp%ip(ikp)%integration_weight
                        ed%pdos_band(ie,iband,ispin)=ed%pdos_band(ie,iband,ispin)+f1
                    enddo

                end block spherefamily
                end select
            enddo

            if ( verbosity .gt. 0 ) then
            if ( lo_trueNtimes(ctr,111,n_kp*ed%n_spin) ) then
                call lo_progressbar(' ... gaussian integration',ctr,ed%n_spin*n_kp,walltime()-t0)
            endif
            endif
        enddo
        enddo

        ! Collect things over ranks
        call mw%allreduce('sum',ed%pdos_band)

        if ( verbosity .gt. 0 ) then
            call lo_progressbar(' ... gaussian integration',ed%n_spin*n_kp,ed%n_spin*n_kp,walltime()-t0)
        endif
    end block actualintegration
end subroutine

!> Integrate the electron dos with the tetrahedron method
subroutine get_electron_dos_hybrid(ed,p,kp,edr,integrationtype,mw,mem,verbosity)
    !> electronic dos
    type(lo_electron_dos), intent(inout) :: ed
    !> structure
    type(lo_crystalstructure), intent(in) :: p
    !> k-point grid
    class(lo_qpoint_mesh), intent(in) :: kp
    !> dispersion relations
    type(lo_electron_dispersions), intent(in) :: edr
    !> how to integration
    integer, intent(in) :: integrationtype
    !> mpi helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem
    !> talk a lot?
    integer, intent(in) :: verbosity

    ! Tetrahedron energies
    real(r8), dimension(:,:,:,:), allocatable :: tetenergies
    real(r8) :: deltae,t0
    integer, dimension(:), allocatable :: ind_tet
    integer :: n_tet,n_kp

    ! Start integration timer
    t0=walltime()

    ! Set basic things
    init: block
        integer :: itet,ikp
        ! Decide on a tiny smearing parameter!
        deltae=(ed%energy(2)-ed%energy(1))*0.5_r8

        ! Divide tetrahedrons across ranks?
        n_tet=0
        do itet=1,kp%n_irr_tet
            if ( mod(itet,mw%n) .eq. mw%r ) n_tet=n_tet+1
        enddo

        if ( n_tet .gt. 0 ) then
            call mem%allocate(ind_tet,n_tet,persistent=.false.,scalable=.true.,file=__FILE__,line=__LINE__)
            ind_tet=0
        endif

        n_tet=0
        do itet=1,kp%n_irr_tet
            if ( mod(itet,mw%n) .eq. mw%r ) then
                n_tet=n_tet+1
                ind_tet(n_tet)=itet
            endif
        enddo

        ! Also count number of k-points per rank
        n_kp=0
        do ikp=1,kp%n_irr_point
            if ( mod(ikp,mw%n) .eq. mw%r ) n_kp=n_kp+1
        enddo
    end block init

    ! Prefetch energies
    select case(integrationtype)
    case(1)
    ! This is just a straight prefetch of tetrahedron energies, nothing magical.
    prefetchstupid: block
        integer :: i,itet,icrn,iband,ispin,ikp
        ! Make space for the tetrahedron energies
        if ( n_tet .gt. 0 ) then
            call mem%allocate(tetenergies,[4,ed%n_band,n_tet,ed%n_spin],persistent=.false.,scalable=.true.,file=__FILE__,line=__LINE__)
            tetenergies=0.0_r8
        endif
        ! Fetch the tetrahedron energies
        do ispin=1,ed%n_spin
        do i=1,n_tet
            itet=ind_tet(i)
            do icrn=1,4
                ikp=kp%it(itet)%irreducible_index(icrn)
                do iband=1,ed%n_band
                    tetenergies(icrn,iband,i,ispin)=edr%ik(ikp)%eigenvalue(iband,ispin)
                enddo
            enddo
        enddo
        enddo
    end block prefetchstupid
    case(2)
    prefetchsmarter: block
        real(r8), parameter :: degentol=1E-6_r8*lo_eV_to_Hartree
        real(r8), dimension(:,:,:), allocatable :: gradient
        real(r8), dimension(:,:), allocatable :: energy
        real(r8), dimension(3,4) :: corner
        real(r8), dimension(3) :: v0,v1
        integer, dimension(:,:), allocatable :: permutation
        integer :: i,itet,icrn,iband,jband,ispin,ikp,jkp,iop,success

        ! Anyhow, make space
        if ( n_tet .gt. 0 ) then
            call mem%allocate(tetenergies,[4,ed%n_band,n_tet,ed%n_spin],persistent=.false.,scalable=.true.,file=__FILE__,line=__LINE__)
            tetenergies=0.0_r8
        endif
        call mem%allocate(gradient,[3,ed%n_band,4],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%allocate(energy,[ed%n_band,4],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%allocate(permutation,[ed%n_band,4],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        gradient=0.0_r8
        energy=0.0_r8
        permutation=0

        if ( verbosity .gt. 0 ) call lo_progressbar_init()

        do ispin=1,ed%n_spin
        do i=1,n_tet
            itet=ind_tet(i)
            ! Fetch absolute coordinates of tetrahedron
            v0=kp%ap( kp%it(itet)%full_index(1) )%r
            do icrn=1,4
                jkp=kp%it(itet)%full_index(icrn)
                v1=kp%ap( jkp )%r
                v1=matmul(p%inv_reciprocal_latticevectors,v1-v0)
                v1=lo_clean_fractional_coordinates(v1+0.5_r8)-0.5_r8
                corner(:,icrn)=matmul(p%reciprocal_latticevectors,v1)+v0
            enddo

            ! Fetch things for tetrahedron
            do icrn=1,4
                jkp=kp%it(itet)%full_index(icrn)
                ikp=kp%it(itet)%irreducible_index(icrn)
                iop=kp%ap(jkp)%operation_from_irreducible
                do iband=1,ed%n_band
                    energy(iband,icrn)=edr%ik( ikp )%eigenvalue(iband,ispin)
                    if ( iop .gt. 0 ) then
                        v0=matmul(p%sym%op(iop)%m,edr%ik( ikp )%groupvelocity(:,iband,ispin))
                    else
                        v0=-matmul(p%sym%op(-iop)%m,edr%ik( ikp )%groupvelocity(:,iband,ispin))
                    endif
                    gradient(:,iband,icrn)=v0
                enddo
            enddo

            ! Figure out the permutation
            call lo_untangle_one_tetrahedron(corner,energy,gradient,degentol,permutation,success)
            ! Store energies using the permutation
            do icrn=1,4
                do iband=1,ed%n_band
                    jband=permutation(iband,icrn)
                    tetenergies(icrn,iband,i,ispin)=energy(jband,icrn)
                enddo
            enddo

            ! Report how it is going
            if ( verbosity .gt. 0 ) then
            if ( lo_trueNtimes(i,111,n_tet) ) then
                call lo_progressbar(' ... untangling tetrahedrons',i,n_tet,walltime()-t0)
            endif
            endif

        enddo
        enddo

        ! Local cleanup
        call mem%deallocate(gradient,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%deallocate(energy,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%deallocate(permutation,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)

        ! And possibly make a final report.
        if ( verbosity .gt. 0 ) then
            call lo_progressbar(' ... untangling tetrahedrons',n_tet,n_tet,walltime()-t0)
        endif
        t0=walltime()
    end block prefetchsmarter
    end select

    ! Do the actual integration.
    actualintegration: block
        real(r8), dimension(4) :: tete,w1,w2
        real(r8) :: e,f0,f1,sigma,fivesigma
        real(r8) :: mine,maxe,invf,totmine,totmaxe
        integer :: i,j,ii,jj,ispin,itet,iband,ie,ctr,ikp

        ! some misc things to make it a little faster
        ed%pdos_band=0.0_r8
        totmine=minval(ed%energy)
        totmaxe=maxval(ed%energy)
        invf=real(ed%n_dos_point,r8)/(totmaxe-totmine)
        sigma=ed%dossmear
        fivesigma=sigma*5.0_r8

        ! Calculate the normal DOS via tetrahedron integration
        if ( verbosity .gt. 0 ) call lo_progressbar_init()
        ctr=0
        do ispin=1,ed%n_spin
        do i=1,n_tet
            ctr=ctr+1
            itet=ind_tet(i)
            do j=1,ed%n_near
                iband=ed%near(j)
                tete=tetenergies(:,iband,i,ispin)
                mine=minval(tete)-totmine
                maxe=maxval(tete)-totmine
                ii=max(floor( mine*invf ),1)
                jj=min(ceiling( maxe*invf ),ed%n_dos_point)
                do ie=ii,jj
                    e=ed%energy(ie)
                    w1=lo_LV_tetrahedron_weights(tete,e-deltae*0.25_r8,1E-10_r8,sigma)*kp%it(itet)%integration_weight
                    w2=lo_LV_tetrahedron_weights(tete,e+deltae*0.25_r8,1E-10_r8,sigma)*kp%it(itet)%integration_weight
                    ed%pdos_band(ie,iband,ispin)=ed%pdos_band(ie,iband,ispin)+sum(w1+w2)*0.5_r8
                enddo
            enddo

            if ( verbosity .gt. 0 ) then
            if ( lo_trueNtimes(i,111,n_tet*ed%n_spin) ) then
                call lo_progressbar(' ... tetrahedron integration',ctr,ed%n_spin*n_tet,walltime()-t0)
            endif
            endif
        enddo
        enddo

        if ( verbosity .gt. 0 ) then
            call lo_progressbar(' ... tetrahedron integration',ed%n_spin*n_tet,ed%n_spin*n_tet,walltime()-t0)
        endif

        t0=walltime()

        ! Also gaussian integrate the DOS for the bands far from the Fermi level.
        if ( verbosity .gt. 0 ) call lo_progressbar_init()
        ctr=0
        do ispin=1,ed%n_spin
        do ikp=1,kp%n_irr_point
            if ( mod(ikp,mw%n) .ne. mw%r ) cycle
            ctr=ctr+1
            do j=1,ed%n_below
                iband=ed%below(j)
                f0=edr%ik(ikp)%eigenvalue(iband,ispin)
                ii=max(floor( (f0-totmine-fivesigma)*invf ),1)
                jj=min(ceiling( (f0-totmine+fivesigma)*invf ),ed%n_dos_point)
                do ie=ii,jj
                    f1=lo_gauss(ed%energy(ie),f0,sigma)*kp%ip(ikp)%integration_weight
                    ed%pdos_band(ie,iband,ispin)=ed%pdos_band(ie,iband,ispin)+f1
                enddo
            enddo
            do j=1,ed%n_above
                iband=ed%above(j)
                f0=edr%ik(ikp)%eigenvalue(iband,ispin)
                ii=max(floor( (f0-totmine-fivesigma)*invf ),1)
                jj=min(ceiling( (f0-totmine+fivesigma)*invf ),ed%n_dos_point)
                do ie=ii,jj
                    f1=lo_gauss(ed%energy(ie),f0,sigma)*kp%ip(ikp)%integration_weight
                    ed%pdos_band(ie,iband,ispin)=ed%pdos_band(ie,iband,ispin)+f1
                enddo
            enddo

            if ( verbosity .gt. 0 ) then
            if ( lo_trueNtimes(ctr,111,n_kp*ed%n_spin) ) then
                call lo_progressbar(' ... gaussian integration',ctr,ed%n_spin*n_kp,walltime()-t0)
            endif
            endif

        enddo
        enddo

        ! Collect things over ranks
        call mw%allreduce('sum',ed%pdos_band)

        ! And now we can clean things no longer in use.
        if ( n_tet .gt. 0 ) then
            call mem%deallocate(tetenergies,persistent=.false.,scalable=.true.,file=__FILE__,line=__LINE__)
            call mem%deallocate(ind_tet,persistent=.false.,scalable=.true.,file=__FILE__,line=__LINE__)
        endif

        if ( verbosity .gt. 0 ) then
            call lo_progressbar(' ... gaussian integration',ed%n_spin*n_kp,ed%n_spin*n_kp,walltime()-t0)
        endif
    end block actualintegration
end subroutine

!> Write electron dos to file
subroutine write_to_file(ed,filename)
    !> electron dos
    class(lo_electron_dos), intent(in) :: ed
    !> optionally a different filename
    character(len=*), intent(in), optional :: filename

    real(r8), parameter :: narrowrange=1.0_r8

    integer :: i,u
    character(len=100) :: opf,fn
    real(r8) :: emin,emax

    ! Maybe another filename
    if ( present(filename) ) then
        fn=trim(filename)
    else
        fn='outfile.electron_dos'
    endif

    ! Print the raw data
    u=open_file('out',trim(fn))
        opf="(2(1X,E18.12))"
        do i=1,ed%n_dos_point
            write(u,opf) ed%energy(i)*lo_Hartree_to_eV,ed%dos(i,:)/lo_Hartree_to_eV
        enddo
    close(u)

    ! select a range for the detailed window
    emin=-narrowrange
    emax=narrowrange+ed%bandgap*lo_hartree_to_eV

    ! Get the gnuplot file
    u=open_file('out',trim(fn)//'.gnuplot')
        ! Choose terminal
        write(u,*) 'set terminal '//lo_gnuplotterminal//' size 700,350 enhanced font "CMU Serif,8"'
        write(u,*) 'set multiplot'
        !
        write(u,*) 'set origin 0.0, 0.0'
        write(u,*) 'set size 0.5, 1.0'
        write(u,*) 'set border lw 0.5'
        write(u,*) 'set ytics scale 0.5'
        write(u,*) 'set xtics scale 0.5'
        write(u,*) 'set mytics 10'
        write(u,*) 'set mxtics 10'
        write(u,*) ' set xlabel "Energy (eV)" '
        write(u,*) ' set ylabel "DOS" '
        write(u,*) 'unset key'
        write(u,'(A)') ' plot "'//trim(fn)//'" u ($1-'//tochar(ed%efermi*lo_Hartree_to_eV)//'):($2) w line lc rgb "#318712"'
        !
        write(u,*) 'set origin 0.5, 0.0'
        write(u,*) 'set size 0.5, 1.0'
        write(u,*) 'set border lw 0.5'
        write(u,*) 'set ytics scale 0.5'
        write(u,*) 'set xtics scale 0.5'
        write(u,*) 'set mytics 10'
        write(u,*) 'set mxtics 10'
        write(u,*) ' set xlabel "Energy-Efermi (eV)" '
        write(u,*) ' set ylabel "DOS" '
        write(u,*) 'unset key'
        write(u,*) 'set xrange ['//tochar(emin)//':'//tochar(emax)//']'
        write(u,*) 'set yrange [0:]'
        write(u,'(A)') ' plot "'//trim(fn)//'" u ($1-'//tochar(ed%efermi*lo_Hartree_to_eV)//'):($2) w line lc rgb "#318712"'
    close(u)
end subroutine

!> Write electron dos to hdf5
subroutine write_to_hdf5(ed,filename,hdftag)
    !> electron dos
    class(lo_electron_dos), intent(in) :: ed
    !> filename
    character(len=*), intent(in) :: filename
    !> optionally, write to an already open hdf5 file.
    integer(i8), intent(in), optional :: hdftag

    type(lo_hdf5_helper) :: h5

    ! create the file, or use a provided tag
    if ( present(hdftag) ) then
        ! Write to the tag that was provided
        h5%file_id=hdftag
    else
        ! Create a new file.
        call h5%init(__FILE__,__LINE__)
        call h5%open_file('write',trim(filename))
    endif

    ! store all the data
    call h5%store_data(ed%energy*lo_Hartree_to_eV   ,h5%file_id,'energies',enhet='eV',dimensions='energy')
    call h5%store_data(ed%dos/lo_Hartree_to_eV      ,h5%file_id,'dos',enhet='states/eV',dimensions='dos')
    call h5%store_data(ed%pdos_band/lo_Hartree_to_eV,h5%file_id,'dos_per_band',enhet='states/eV',dimensions='band,dos')
    ! store some metadata
    call h5%store_attribute(ed%efermi*lo_Hartree_to_eV,h5%file_id,'fermi_level')
    ! close the file and hdf5, if relevant.
    if ( present(hdftag) .eqv. .false. ) then
        call h5%close_file()
        call h5%destroy(__FILE__,__LINE__)
    endif
end subroutine

!> Destroy and deallocate
subroutine destroy(pd)
    !> electron density of states
    class(lo_electron_dos), intent(inout) :: pd

    ! Deallocate everything
    if ( allocated(pd%energy   ) ) deallocate(pd%energy   )
    if ( allocated(pd%dos      ) ) deallocate(pd%dos      )
    if ( allocated(pd%pdos_band) ) deallocate(pd%pdos_band)
    if ( allocated(pd%below    ) ) deallocate(pd%below    )
    if ( allocated(pd%above    ) ) deallocate(pd%above    )
    if ( allocated(pd%near     ) ) deallocate(pd%near     )
    if ( allocated(pd%core     ) ) deallocate(pd%core     )
    if ( allocated(pd%sky      ) ) deallocate(pd%sky      )
end subroutine

!> Measure size in memory, in bytes
pure function size_in_mem(pd) result(mem)
    !> electron density of states
    class(lo_electron_dos), intent(in) :: pd
    !> memory in bytes
    integer :: mem

    mem=0
    mem=mem+storage_size(pd)
    if ( allocated(pd%energy   ) ) mem=mem+storage_size(pd%energy   )*size(pd%energy   )
    if ( allocated(pd%dos      ) ) mem=mem+storage_size(pd%dos      )*size(pd%dos      )
    if ( allocated(pd%pdos_band) ) mem=mem+storage_size(pd%pdos_band)*size(pd%pdos_band)
    if ( allocated(pd%below    ) ) mem=mem+storage_size(pd%below    )*size(pd%below    )
    if ( allocated(pd%above    ) ) mem=mem+storage_size(pd%above    )*size(pd%above    )
    if ( allocated(pd%near     ) ) mem=mem+storage_size(pd%near     )*size(pd%near     )
    if ( allocated(pd%core     ) ) mem=mem+storage_size(pd%core     )*size(pd%core     )
    if ( allocated(pd%sky      ) ) mem=mem+storage_size(pd%sky      )*size(pd%sky      )
    mem=mem/8
end function

end module
