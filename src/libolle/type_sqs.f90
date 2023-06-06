#include "precompilerdefinitions"
module type_sqs
    ! I use the formalism in
    ! A. van de Walle, Calphad 33, 266 (2009) doi:10.1016/j.calphad.2008.12.005
    ! to generate SQS structures
    use konstanter, only: flyt,lo_huge,lo_hugeint,lo_pi,lo_tol,lo_sqtol,lo_status,lo_exitcode_symmetry
    use gottochblandat, only: walltime,tochar,lo_stddev,lo_mean,lo_sqnorm,qsort,lo_progressbar_init,lo_progressbar
    use mpi_wrappers, only: lo_mpi_helper,lo_stop_gracefully
    use type_crystalstructure, only: lo_crystalstructure
    use type_symmetrylist, only: lo_symlist
    use lo_randomnumbers, only: lo_mersennetwister

    implicit none
    private
    public :: lo_sqs

    !> a symmetry irreducible shell
    type, abstract :: lo_sqs_shell
        !> Identifier of this shell in the grand scheme of things
        integer :: tag=-lo_hugeint
        !> weighting for this shell
        real(flyt) :: weight=-lo_huge
        !> reference cluster function
        real(flyt) :: ideal_clusterfunction=-lo_huge
        !> cluster function
        real(flyt) :: clusterfunction=-lo_huge
    end type

    !> a relevant coordination shell
    type, extends(lo_sqs_shell) :: lo_sqs_pairshell
        !> How many pairs in the supercell map to this shell?
        integer :: npair=lo_hugeint
        !> pair distance of this shell
        real(flyt) :: r=-lo_huge
        !> Species
        integer :: sp1=-lo_hugeint,sp2=-lo_hugeint
        !> Number of components
        integer :: m1=-lo_hugeint,m2=-lo_hugeint
        !> indices to the first atom in the pair
        integer, dimension(:), allocatable :: i1
        !> and the second atom
        integer, dimension(:), allocatable :: i2

        real(flyt), dimension(:,:), allocatable :: ideal_correlationfunction
        real(flyt), dimension(:,:), allocatable :: correlationfunction
    end type

    !> collected information about each species
    type lo_sqs_species
        !> how many shells involve this species
        integer :: npairshell=-lo_hugeint
        !> coordination shell
        type(lo_sqs_pairshell), dimension(:), allocatable :: pairshell

        !> Number of components
        integer :: ncomponent=-lo_hugeint
        !> Possible components
        integer, dimension(:), allocatable :: component
        !> Desired concentration
        real(flyt), dimension(:), allocatable :: concentration

        !> Number of sites for this speies
        integer :: nsites
        !> Which sites
        integer, dimension(:), allocatable :: sites
        !> What is the current occupation
        integer, dimension(:), allocatable :: occupation
    end type

    !> container that handles special quasi-random structures
    type lo_sqs
        !> how many species
        integer :: nsp=-lo_hugeint
        !> how many atoms in the supercell
        integer :: nss=-lo_hugeint
        !> max number of pair correlation thingies
        integer :: max_n_paircf=-lo_hugeint
        !> One container per species
        type(lo_sqs_species), dimension(:), allocatable :: sp
        !> Perhaps I want more than one SQS?
        integer, dimension(:,:), allocatable :: candidate_occupation
        contains
            !> Set up everything
            procedure :: generate
            !> get the correlation function
            procedure :: correlation_function_distance
            !> take the current configuration and return a crystalstructure
            procedure :: returnstructure
    end type

contains

!> return a proper poscar
subroutine returnstructure(sqs,ss,p,cnf,sort,verbosity)
    !> disorder generator thingy
    class(lo_sqs), intent(inout) :: sqs
    !> the perfect alloy supercell thingy
    type(lo_crystalstructure), intent(in) :: ss
    !> normal crystal structure thingy
    type(lo_crystalstructure), intent(out) :: p
    !> which configuration to generate
    integer, intent(in) :: cnf
    !> sort it or keep order
    logical, intent(in) :: sort
    !> how much to talk
    integer, intent(in), optional :: verbosity

    integer :: verb

    ! Talk?
    if ( present(verbosity) ) then
        verb=verbosity
    else
        verb=0
    endif

    build: block
        real(flyt), dimension(3,3) :: latticevectors
        real(flyt), dimension(:,:), allocatable :: positions
        integer, dimension(:), allocatable :: atomic_numbers,di

        latticevectors=ss%latticevectors
        lo_allocate(positions(3,ss%na))
        lo_allocate(atomic_numbers(ss%na))
        lo_allocate(di(ss%na))
        positions=ss%r

        atomic_numbers=sqs%candidate_occupation(:,cnf)
        ! sort it properly, maybe
        if ( sort ) then
            call qsort(atomic_numbers,di)
            positions=positions(:,di)
        endif
        ! Build a structure from this
        call p%generate(latticevectors=ss%latticevectors,positions=positions,atomic_numbers=atomic_numbers,enhet=2)
    end block build
end subroutine

!> Set everything up!
subroutine generate(sqs,uc,ss,sl,ncandidate,verbosity,mw)
    !> disorder generator thingy
    class(lo_sqs), intent(out) :: sqs
    !> unitcell
    type(lo_crystalstructure), intent(in) :: uc
    !> supercell
    type(lo_crystalstructure), intent(in) :: ss
    !> symmetries
    type(lo_symlist), intent(in) :: sl
    !> how many SQSs to generate
    integer, intent(in) :: ncandidate
    !> how much to talk
    integer, intent(in) :: verbosity
    !> MPI helper
    type(lo_mpi_helper), intent(inout) :: mw

    type(lo_mersennetwister) :: tw
    real(flyt) :: timer

    timer=walltime()
    if ( verbosity .gt. 0 ) then
        write(*,*) ''
        write(*,*) 'Building SQS'
    endif

    ! Start by seeding the random numbers
    call tw%init(iseed=mw%r,rseed=walltime())

    ! rearrange shells into a more convenient form
    buildshells: block
        real(flyt), dimension(3) :: v
        integer :: i,j,l,ll,sh,m1,m2
        integer :: a1,a2,sp1,sp2,ssp1,ssp2

        ! First step is to get the number of sublattices
        sqs%nss=ss%na
        sqs%nsp=uc%nelements
        lo_allocate(sqs%sp( sqs%nsp ))

        ! Count pair shells that start and end from this species, and store them
        do i=1,sqs%nsp

            ! count shells
            sqs%sp(i)%npairshell=0
            do sh=1,sl%npairshells
                sp1=uc%species( sl%pairshell(sh)%protpair%ui1 )
                sp2=uc%species( sl%pairshell(sh)%protpair%ui2 )
                m1=uc%alloyspecies( sp1 )%n
                m2=uc%alloyspecies( sp2 )%n
                v=sl%pairshell(sh)%protpair%v
                if ( sp1 .eq. i .or. sp2 .eq. i ) then
                if ( m1 .gt. 1 .and. m2 .gt. 1 ) then
                if ( lo_sqnorm(v) .gt. lo_sqtol ) then
                    sqs%sp(i)%npairshell=sqs%sp(i)%npairshell+1
                endif
                endif
                endif
            enddo
            lo_allocate(sqs%sp(i)%pairshell( sqs%sp(i)%npairshell ))

            ! identify shells
            l=0
            do sh=1,sl%npairshells
                sp1=uc%species( sl%pairshell(sh)%protpair%ui1 )
                sp2=uc%species( sl%pairshell(sh)%protpair%ui2 )
                m1=uc%alloyspecies(sp1)%n
                m2=uc%alloyspecies(sp2)%n
                v=sl%pairshell(sh)%protpair%v
                if ( sp1 .eq. i .or. sp2 .eq. i ) then
                if ( m1 .gt. 1 .and. m2 .gt. 1 ) then
                if ( lo_sqnorm(v) .gt. lo_sqtol ) then
                    l=l+1
                    ! Store some extra things to be on the safe side
                    sqs%sp(i)%pairshell(l)%r=norm2(v)
                    sqs%sp(i)%pairshell(l)%weight=1.0_flyt/(norm2(v)**4)
                    sqs%sp(i)%pairshell(l)%sp1=sp1
                    sqs%sp(i)%pairshell(l)%sp2=sp2
                    sqs%sp(i)%pairshell(l)%m1=m1
                    sqs%sp(i)%pairshell(l)%m2=m2
                    ! store the unique identifier for this shell
                    sqs%sp(i)%pairshell(l)%tag=sh
                    ! Count the number of supercell pairs that match this shell.
                    ! I take care to only count half the pairs, because transpose blablabla
                    ll=0
                    do a1=1,sl%nss
                    do j=1,sl%ss(a1)%npair
                        if ( sl%ss(a1)%pair(j)%unique_shell .ne. sqs%sp(i)%pairshell(l)%tag ) cycle
                        a2=sl%ss(a1)%pair(j)%j2
                        ssp1=ss%species(a1)
                        ssp2=ss%species(a2)
                        if ( sp1 .ne. ssp1 ) cycle
                        if ( sp2 .ne. ssp2 ) cycle
                        if ( sp1 .eq. sp2 .and. a2 .gt. a1 ) cycle
                        ll=ll+1 ! yup, good pair
                    enddo
                    enddo
                    sqs%sp(i)%pairshell(l)%npair=ll

                    ! Store the supercell pairs
                    lo_allocate(sqs%sp(i)%pairshell(l)%i1( ll ))
                    lo_allocate(sqs%sp(i)%pairshell(l)%i2( ll ))
                    ll=0
                    do a1=1,sl%nss
                    do j=1,sl%ss(a1)%npair
                        if ( sl%ss(a1)%pair(j)%unique_shell .ne. sqs%sp(i)%pairshell(l)%tag ) cycle
                        a2=sl%ss(a1)%pair(j)%j2
                        ssp1=ss%species(a1)
                        ssp2=ss%species(a2)
                        if ( sp1 .ne. ssp1 ) cycle
                        if ( sp2 .ne. ssp2 ) cycle
                        if ( sp1 .eq. sp2 .and. a2 .gt. a1 ) cycle
                        ll=ll+1 ! yup, good pair
                        sqs%sp(i)%pairshell(l)%i1( ll )=a1
                        sqs%sp(i)%pairshell(l)%i2( ll )=a2
                    enddo
                    enddo

                endif
                endif
                endif
            enddo

            if ( verbosity .gt. 0 ) then
                write(*,*) 'species',i,'nshell',sqs%sp(i)%npairshell
                do sh=1,sqs%sp(i)%npairshell
                    write(*,*) '    shell',tochar(sh,-3),' npair',tochar(sqs%sp(i)%pairshell(sh)%npair,-5),sqs%sp(i)%pairshell(sh)%r
                enddo
            endif
        enddo

        if ( verbosity .gt. 0 ) write(*,*) '... rearranged shells (',tochar(walltime()-timer),'s)'

    end block buildshells

    !> Setup the ideally random correlation functions
    setupcorr: block
        integer, dimension(:), allocatable :: di
        integer :: i,j,l,s,m1,m2,sh

        do s=1,sqs%nsp
            ! Which sites in the supercell are of this species?
            l=0
            do i=1,ss%na
                j=ss%species(i)
                m1=uc%alloyspecies(s)%n
                m2=ss%alloyspecies(j)%n
                if ( m1 .ne. m2 ) cycle
                if ( sum(ss%alloyspecies(j)%atomic_number-uc%alloyspecies(s)%atomic_number ) .ne. 0 ) cycle
                if ( sum(abs(ss%alloyspecies(j)%concentration-uc%alloyspecies(s)%concentration )) .gt. lo_tol ) cycle
                l=l+1
            enddo
            sqs%sp(s)%nsites=l
            lo_allocate(sqs%sp(s)%sites(l))
            l=0
            do i=1,ss%na
                j=ss%species(i)
                m1=uc%alloyspecies(s)%n
                m2=ss%alloyspecies(j)%n
                if ( m1 .ne. m2 ) cycle
                if ( sum(ss%alloyspecies(j)%atomic_number-uc%alloyspecies(s)%atomic_number ) .ne. 0 ) cycle
                if ( sum(abs(ss%alloyspecies(j)%concentration-uc%alloyspecies(s)%concentration )) .gt. lo_tol ) cycle
                l=l+1
                sqs%sp(s)%sites(l)=i
            enddo
            ! Some space for the occupation
            lo_allocate(sqs%sp(s)%occupation( sqs%sp(s)%nsites ))
            sqs%sp(s)%occupation=0

            ! Small sanity check
            if ( sqs%sp(s)%nsites .ne. ss%element_counter(s) ) then
                write(*,*) 'Inconsistency with sites'
                stop
            endif

            ! store the alloy stuff
            sqs%sp(s)%ncomponent=uc%alloyspecies(s)%n
            lo_allocate(sqs%sp(s)%component( sqs%sp(s)%ncomponent ))
            lo_allocate(sqs%sp(s)%concentration( sqs%sp(s)%ncomponent ))
            sqs%sp(s)%component=uc%alloyspecies(s)%atomic_number
            sqs%sp(s)%concentration=uc%alloyspecies(s)%concentration

            ! Adjust the concentration slightly, since with a certain supercell you can't really get what you want
            if ( sqs%sp(s)%ncomponent .gt. 1 ) then
                lo_allocate(di(sqs%sp(s)%ncomponent))
                do i=1,sqs%sp(s)%ncomponent-1
                    j=ceiling( sqs%sp(s)%concentration(i)*sqs%sp(s)%nsites )
                    j=max(j,1)
                    j=min(j,sqs%sp(s)%nsites-sqs%sp(s)%ncomponent+1)
                    di(i)=j
                enddo
                di(sqs%sp(s)%ncomponent)=sqs%sp(s)%nsites-sum(di(1:sqs%sp(s)%ncomponent-1))

                do i=1,sqs%sp(s)%ncomponent
                    sqs%sp(s)%concentration(i)=(1.0_flyt*di(i))/(1.0_flyt*sqs%sp(s)%nsites)
                enddo
                lo_deallocate(di)

                if ( verbosity .gt. 0 ) then
                    write(*,*) 'Sublattice ',tochar(s)
                    write(*,*) '     desired concentrations: ',tochar(uc%alloyspecies(s)%concentration)
                    write(*,*) '    possible concentrations: ',tochar(sqs%sp(s)%concentration)
                endif
            endif

            ! Initialize the correlation functions
            do sh=1,sqs%sp(s)%npairshell
                m1=sqs%sp(s)%pairshell(sh)%m1-1
                m2=sqs%sp(s)%pairshell(sh)%m2-1
                sqs%max_n_paircf=max(sqs%max_n_paircf,m1+m2)
                ! Space for the correlation functions
                lo_allocate( sqs%sp(s)%pairshell(sh)%ideal_correlationfunction(m1,m2) )
                lo_allocate( sqs%sp(s)%pairshell(sh)%correlationfunction(m1,m2) )
                sqs%sp(s)%pairshell(sh)%ideal_correlationfunction=0.0_flyt
                sqs%sp(s)%pairshell(sh)%correlationfunction=0.0_flyt
            enddo
        enddo
    end block setupcorr

    ! Get the numbers for the ideal correlation functions, and an initial guess for the occupation
    idealcorr: block
        real(flyt), dimension(:), allocatable :: dcf
        real(flyt) :: f0,f1
        integer, dimension(:), allocatable :: dj
        integer :: i,j,k,l,s,al1,al2,m1,m2,sh,sp1,sp2,nc,ns,ali

        lo_allocate(dcf(sqs%max_n_paircf))
        do s=1,sqs%nsp
            ! Initialize the correlation functions
            do sh=1,sqs%sp(s)%npairshell
                ! Set the ideal value
                m1=sqs%sp(s)%pairshell(sh)%m1
                m2=sqs%sp(s)%pairshell(sh)%m2
                sp1=sqs%sp(s)%pairshell(sh)%sp1
                sp2=sqs%sp(s)%pairshell(sh)%sp2
                do al1=1,m1-1
                do al2=1,m2-1
                    f0=0.0_flyt
                    f1=0.0_flyt
                    do i=1,m1
                        f0=f0+smallgamma(al1,i-1,m1)*sqs%sp(sp1)%concentration(i)
                    enddo
                    do i=1,m2
                        f1=f1+smallgamma(al2,i-1,m2)*sqs%sp(sp2)%concentration(i)
                    enddo
                    sqs%sp(s)%pairshell(sh)%ideal_correlationfunction(al1,al2)=f0*f1
                enddo
                enddo
            enddo

            do sh=1,sqs%sp(s)%npairshell
                ! Set the ideal value
                sp1=sqs%sp(s)%pairshell(sh)%sp1
                sp2=sqs%sp(s)%pairshell(sh)%sp2
                m1=sqs%sp(s)%pairshell(sh)%m1
                m2=sqs%sp(s)%pairshell(sh)%m2
                ! number of things in the correlation function
                ns=(m1-1)+(m2-1)
                dcf=0.0_flyt
                l=0
                do ali=1,m1-1
                    l=l+1
                    do i=1,m1
                        dcf(l)=dcf(l)+smallgamma(ali,i-1,m1)*sqs%sp(sp1)%concentration(i)
                    enddo
                enddo
                do ali=1,m2-1
                    l=l+1
                    do i=1,m2
                        dcf(l)=dcf(l)+smallgamma(ali,i-1,m2)*sqs%sp(sp2)%concentration(i)
                    enddo
                enddo
                sqs%sp(s)%pairshell(sh)%ideal_clusterfunction=product(dcf)
            enddo
        enddo
        if ( verbosity .gt. 0 ) write(*,*) '... got ideal correlationfunctions (',tochar(walltime()-timer),'s)'

        ! First set an initial configuration
        do s=1,sqs%nsp
            nc=sqs%sp(s)%ncomponent
            ns=sqs%sp(s)%nsites
            lo_allocate(dj(ns))
            l=0
            do i=1,nc
                j=int(anint(sqs%sp(s)%concentration(i)*ns))
                do k=1,j
                    l=l+1
                    dj(l)=i
                enddo
            enddo
            if ( l .ne. ns ) then
                call lo_stop_gracefully(['Something went wrong when generating SQS'],lo_exitcode_symmetry,__FILE__,__LINE__)
            endif
            ! shuffle it
            call tw%shuffle_int_array(dj)
            sqs%sp(s)%occupation=dj
            lo_deallocate(dj)
        enddo
        if ( verbosity .gt. 0 ) write(*,*) '... got initial guess (',tochar(walltime()-timer),'s)'
    end block idealcorr

    !> Now optimize the structure
    optimize: block
        ! Magic parameters for the minimization
        integer, parameter :: nouter=1000
        integer, parameter :: nconv=15
        real(flyt), dimension(:), allocatable :: bdr
        real(flyt), dimension(10) :: convarr
        real(flyt), dimension(ncandidate) :: cand_corr
        real(flyt) :: dist_curr,dist_new,temperature,f0,f1,f2,thres,scalefactor
        integer, dimension(:,:), allocatable :: bdi
        integer, dimension(ncandidate) :: cand_di
        integer, dimension(:), allocatable :: oldocc,dj
        integer :: outiter,initer,ninner,inflip,cctr,candctr
        integer :: i,j,k,s
        logical :: newcorr

        lo_allocate(oldocc(sqs%nss))
        lo_allocate(sqs%candidate_occupation(sqs%nss,ncandidate))
        sqs%candidate_occupation=-1
        cand_corr=1E10_flyt !lo_huge

        ! All random numbers, subject to change. Seems close enough though.
        convarr=0.0_flyt
        temperature=10000.0_flyt
        ninner=ceiling(sqrt(real(ss%na,flyt)))*sqs%max_n_paircf
        thres=1E-3_flyt
        cctr=0
        candctr=0

        if ( verbosity .ge. 0 ) call lo_progressbar_init()
        do outiter=1,nouter
            ! Current distance
            dist_curr=sqs%correlation_function_distance()
            inflip=0
            if ( outiter .eq. 1 ) scalefactor=0.0_flyt
            do initer=1,ninner
                ! Fetch the current occupation
                do s=1,sqs%nsp
                do i=1,sqs%sp(s)%nsites
                    j=sqs%sp(s)%sites(i)
                    oldocc(j)=sqs%sp(s)%occupation(i)
                enddo
                enddo

                ! Flip it around. If very warm, flip all atoms to end up with a resonable thing to start from
                if ( temperature .gt. 10.0_flyt ) then
                    call CHANGEPLACES(sqs,ss,tw,changeall=.true.)
                else
                    call CHANGEPLACES(sqs,ss,tw,changeall=.false.)
                endif
                dist_new=sqs%correlation_function_distance()

                ! This might be a good one, store it!
                if ( dist_new .lt. cand_corr(ncandidate) ) then
                    newcorr=.true.
                    do i=1,candctr
                        if ( abs(cand_corr(i)-dist_new) .lt. 1E-12_flyt ) then
                            newcorr=.false.
                            exit
                        endif
                    enddo
                    if ( newcorr ) then
                        ! Store this configuration!
                        if ( candctr .lt. ncandidate ) candctr=candctr+1
                        ! Store this configuration
                        do s=1,sqs%nsp
                        do i=1,sqs%sp(s)%nsites
                            j=sqs%sp(s)%sites(i)
                            k=sqs%sp(s)%occupation(i)
                            sqs%candidate_occupation(j,ncandidate)=sqs%sp(s)%component(k)
                        enddo
                        enddo
                        ! And store the correlation distance
                        cand_corr(ncandidate)=dist_new
                        ! Now sort it
                        call qsort(cand_corr,cand_di)
                        sqs%candidate_occupation=sqs%candidate_occupation(:,cand_di)
                    endif
                endif

                f0=dist_curr-dist_new
                f0=exp(f0/temperature)
                f1=tw%rnd_real()
                if ( f0 .gt. f1 ) then
                    inflip=inflip+1
                    dist_curr=dist_new
                else
                    ! reset the occupation
                    do s=1,sqs%nsp
                    do i=1,sqs%sp(s)%nsites
                        j=sqs%sp(s)%sites(i)
                        sqs%sp(s)%occupation(i)=oldocc(j)
                    enddo
                    enddo
                endif
                if ( outiter .eq. 1 ) scalefactor=scalefactor+dist_new/ninner
            enddo

            if ( inflip .gt. ninner/2 ) then
                temperature=temperature*0.5_flyt
            elseif ( inflip .lt. max(ninner/50,1) ) then
                temperature=temperature*1.001_flyt
            else
                temperature=temperature*0.99_flyt
            endif

            ! Figure out if we are converged
            if ( outiter .gt. 10 ) then
                convarr(2:10)=convarr(1:9)
                convarr(1)=cand_corr(1)/(cand_corr(candctr)-cand_corr(1))
                f0=lo_mean(convarr)
                f1=lo_stddev(convarr)
            else
                i=outiter
                convarr(i)=cand_corr(1)/(cand_corr(candctr)-cand_corr(1))
                f0=lo_mean(convarr(1:i))
                f1=lo_stddev(convarr(1:i))
            endif

            if ( f1 .lt. 1E-10_flyt ) then
                cctr=cctr+1
            else
                cctr=0
            endif

            ! Check for convergence
            if ( cctr .ge. nconv ) exit

            if ( outiter .lt. nouter .and. verbosity .ge. 0 ) then
                call lo_progressbar(' ... optimizing correlationfunction',outiter,nouter,walltime()-timer)
            endif
        enddo
        if ( verbosity .ge. 0 ) call lo_progressbar(' ... optimizing correlationfunction',nouter,nouter,walltime()-timer)

        ! Now pick the best candidates from all ranks.
        lo_allocate(bdr(mw%n*ncandidate))
        lo_allocate(dj(mw%n*ncandidate))
        lo_allocate(bdi(ss%na,mw%n*ncandidate))
        bdr=0.0_flyt
        bdi=0
        dj=0
        do i=1,ncandidate
            j=mw%r*ncandidate+i
            bdr(j)=cand_corr(i)
            bdi(:,j)=sqs%candidate_occupation(:,i)
        enddo
        call mw%allreduce('sum',bdr)
        call mw%allreduce('sum',bdi)
        call qsort(bdr,dj)
        sqs%candidate_occupation(:,1)=bdi(:,dj(1))
        j=1
        do i=2,mw%n*ncandidate
            if ( abs(bdr(i)-bdr(i-1)) .lt. lo_sqtol ) cycle
            j=j+1
            sqs%candidate_occupation(:,j)=bdi(:,dj(i))
            if ( j .eq. ncandidate ) exit
        enddo

        ! Now set this occupation? Maybe?
        do s=1,sqs%nsp
        do i=1,sqs%sp(s)%nsites
            j=sqs%sp(s)%sites(i)
            do k=1,sqs%sp(s)%ncomponent
                if ( sqs%candidate_occupation(j,1) .eq. sqs%sp(s)%component(k) ) then
                    sqs%sp(s)%occupation(i)=k
                endif
            enddo
        enddo
        enddo
        dist_new=sqs%correlation_function_distance()

        ! Short report how the optimization went
        if ( verbosity .gt. 0 ) then
            write(*,*) ''
            write(*,*) 'Final correlation functions:'
            do s=1,sqs%nsp
            do i=1,sqs%sp(s)%npairshell
                f0=(abs(sqs%sp(s)%pairshell(i)%clusterfunction))
                f1=(abs(sqs%sp(s)%pairshell(i)%ideal_clusterfunction))
                f2=(abs(sqs%sp(s)%pairshell(i)%clusterfunction-sqs%sp(s)%pairshell(i)%ideal_clusterfunction))
                write(*,*) s,i,f0,f1,f2
            enddo
            enddo
        endif

!        ! Ok this is a decent configuration. Now generate a
!        ! couple of new ones close to this?
!        ! Fetch the current occupation
!        do s=1,sqs%nsp
!        do i=1,sqs%sp(s)%nsites
!            j=sqs%sp(s)%sites(i)
!            oldocc(j)=sqs%sp(s)%occupation(i)
!        enddo
!        enddo
!        sqs%candidate_occupation(:,1)=oldocc
!        cand_corr(1)=sqs%correlation_function_distance()
!
!        ! How many configurations do I have?
!        cctr=1
!        do outiter=1,500 !100
!            ! Fetch the current occupation
!            do s=1,sqs%nsp
!            do i=1,sqs%sp(s)%nsites
!                j=sqs%sp(s)%sites(i)
!                oldocc(j)=sqs%sp(s)%occupation(i)
!            enddo
!            enddo
!            ! How many guys will I flip?
!            j=tw%rnd_int(max(ss%na/20,1))
!            do i=1,j
!                call CHANGEPLACES(sqs,ss,tw,changeall=.false.)
!            enddo
!            f0=sqs%correlation_function_distance()
!            ! Is this a new correlation distance?
!            newcorr=.true.
!            do i=1,cctr
!                if ( abs(f0-cand_corr(i)) .lt. 1E-12_flyt ) then
!                    newcorr=.false.
!                    exit
!                endif
!            enddo
!
!            ! If it's new, try to add it
!            if ( newcorr ) then
!            if ( f0 .lt. cand_corr(ncandidate) ) then
!                if ( cctr .lt. ncandidate ) cctr=cctr+1
!                ! Store this configuration
!                do s=1,sqs%nsp
!                do i=1,sqs%sp(s)%nsites
!                    j=sqs%sp(s)%sites(i)
!                    sqs%candidate_occupation(j,ncandidate)=sqs%sp(s)%occupation(i)
!                enddo
!                enddo
!                ! And store the correlation distance
!                cand_corr(ncandidate)=f0
!                ! Now sort it
!                call qsort(cand_corr,cand_di)
!                !sqs%candidate_occupation(:,cand_di)=sqs%candidate_occupation
!                sqs%candidate_occupation=sqs%candidate_occupation(:,cand_di) !=sqs%candidate_occupation
!            endif
!            endif
!
!            ! Reset to the original occupation
!            do s=1,sqs%nsp
!            do i=1,sqs%sp(s)%nsites
!                j=sqs%sp(s)%sites(i)
!                sqs%sp(s)%occupation(i)=oldocc(j)
!            enddo
!            enddo
!        enddo

    end block optimize

end subroutine

!> Switch place of two atoms
subroutine CHANGEPLACES(sqs,ss,tw,changeall)
    !> disorder generator thingy
    class(lo_sqs), intent(inout) :: sqs
    !> supercell
    type(lo_crystalstructure), intent(in) :: ss
    !> mersenne twister
    type(lo_mersennetwister), intent(inout) :: tw
    !> shuffle all atoms
    logical, intent(in) :: changeall

    integer :: s,i,j,k,l

    if ( changeall ) then
        do s=1,sqs%nsp
            call tw%shuffle_int_array(sqs%sp(s)%occupation)
        enddo
    else
        do
            i=tw%rnd_int(ss%na)
            s=ss%species(i)
            if ( ss%alloyspecies(s)%n .gt. 1 ) exit
        enddo
        ! Now I have a sublattice, switch to atoms
        do
            i=tw%rnd_int( sqs%sp(s)%nsites )
            j=tw%rnd_int( sqs%sp(s)%nsites )
            if ( sqs%sp(s)%occupation(i) .ne. sqs%sp(s)%occupation(j) ) exit
        enddo
        k=sqs%sp(s)%occupation(i)
        l=sqs%sp(s)%occupation(j)
        sqs%sp(s)%occupation(i)=l
        sqs%sp(s)%occupation(j)=k
    endif
end subroutine

!> calculate the distance from the ideal correlation function
function correlation_function_distance(sqs) result(distance)
    !> disorder generator thingy
    class(lo_sqs), intent(inout) :: sqs
    !> distance from the ideal
    real(flyt) :: distance

    real(flyt), dimension(sqs%max_n_paircf) :: dcf
    real(flyt) :: f0,pref
    integer, dimension(:), allocatable :: occupation
    integer :: s,p,sh,i,j,l,m1,m2,sig1,sig2
    integer :: ali,ns

    lo_allocate(occupation(sqs%nss))
    occupation=-1
    do s=1,sqs%nsp
        do i=1,sqs%sp(s)%nsites
            j=sqs%sp(s)%sites(i)
            occupation(j)=sqs%sp(s)%occupation(i)
        enddo
    enddo

    distance=0.0_flyt
    do s=1,sqs%nsp
    do sh=1,sqs%sp(s)%npairshell
        pref=1.0_flyt/( sqs%sp(s)%pairshell(sh)%npair )
        m1=sqs%sp(s)%pairshell(sh)%m1
        m2=sqs%sp(s)%pairshell(sh)%m2
        ns=(m1-1)+(m2-1)
        sqs%sp(s)%pairshell(sh)%clusterfunction=0.0_flyt
        do p=1,sqs%sp(s)%pairshell(sh)%npair
            sig1=occupation( sqs%sp(s)%pairshell(sh)%i1(p) )-1
            sig2=occupation( sqs%sp(s)%pairshell(sh)%i2(p) )-1
            dcf=0.0_flyt
            l=0
            do ali=1,m1-1
                l=l+1
                dcf(l)=smallgamma(ali,sig1,m1)
            enddo
            do ali=1,m2-1
                l=l+1
                dcf(l)=smallgamma(ali,sig2,m2)
            enddo
            sqs%sp(s)%pairshell(sh)%clusterfunction=sqs%sp(s)%pairshell(sh)%clusterfunction+product(dcf(1:ns))*pref
        enddo
        f0=abs(sqs%sp(s)%pairshell(sh)%ideal_clusterfunction-sqs%sp(s)%pairshell(sh)%clusterfunction)
        distance=distance+f0*sqs%sp(s)%pairshell(sh)%weight
    enddo
    enddo

    lo_deallocate(occupation)
end function

!> helper to create the correlation functions
pure function smallgamma(alphai,sigmai,mi) result(sg)
    !> which sublattice subcomponent, 1 .. mi-1
    integer, intent(in) :: alphai
    !> number of possible components
    integer, intent(in) :: mi
    !> 0 .. mi-1, the current occupation
    integer, intent(in) :: sigmai
    !> gamma
    real(flyt) :: sg

    integer :: i
    real(flyt) :: f0

    select case(mi)
    case(1)
        sg=0.0_flyt
    case(2)
        ! in this case alphai has to be 1
        sg=-(-1)**sigmai
    case default
        i=mod(alphai,2)
        f0=2.0_flyt*lo_pi*ceiling( alphai*0.5_flyt )*sigmai/(1.0_flyt*mi)
        select case(i)
            case(1) ! it's odd
                sg=-cos(f0)
            case(0) ! it's even
                sg=-sin(f0)
        end select
    end select
end function

end module
