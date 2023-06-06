#include "precompilerdefinitions"
module magneticdisorder
use konstanter, only: flyt, lo_huge, lo_hugeint, lo_sqtol, lo_tol
use gottochblandat, only: tochar, lo_mean, open_file, lo_sqnorm, lo_chop, walltime
use geometryfunctions, only: lo_rotation_matrix_from_vector_a_to_b
use type_crystalstructure, only: lo_crystalstructure
use type_symmetrylist, only: lo_symlist
use type_forcemap, only: lo_forcemap
use type_blas_lapack_wrappers, only: lo_dgesvd
use lo_randomnumbers, only: lo_mersennetwister
implicit none
private
public :: lo_magdisorder

type lo_magdisorder_shell
    ! How many pairs in the supercell map to this shell?
    integer :: npair = -lo_hugeint
    ! indices to the first atom in the pair
    integer, dimension(:), allocatable :: i1
    ! and the second atom
    integer, dimension(:), allocatable :: i2
end type

type lo_magdisorder
    !> ferromagnetic?
    logical :: ferromagnetic = .false.
    !> collinear or noncollinear
    logical :: coll = .false.
    !> How many bins of different levels of disorder do we want?
    integer :: nbin = -lo_hugeint
    !> How many configurations per bin
    integer :: nconf = -lo_hugeint
    !> number of magnetic coordination shells
    integer :: nshell = -lo_hugeint
    !> info about the coordination shells
    type(lo_magdisorder_shell), dimension(:), allocatable :: sh
    !> history of configurations
    integer, dimension(:, :, :), allocatable :: collhistory
    real(flyt), dimension(:, :, :, :), allocatable :: noncollhistory
    !> initial configuration
    integer, dimension(:), allocatable :: initial_collinear_configuration
    real(flyt), dimension(:, :), allocatable :: initial_noncollinear_configuration
    !> which sites are switchable?
    integer, dimension(:), allocatable :: sites
contains
    !> create the structure
    procedure :: generate
    !> get the correlation function
    procedure :: correlation_function
    !> generate magnetic sqs
    procedure :: optimize
    !> dump to file
    procedure :: dump_configurations
end type

contains

!> make sure an AFM configuration has net magnetic moment of 0
subroutine zerosum(x, rel)
    real(flyt), dimension(:, :), intent(inout) :: x
    logical, dimension(:), intent(in) :: rel
    !
    real(flyt), dimension(3) :: v
    integer :: i, j, n

    v = 0.0_flyt
    n = size(x, 2)
    j = 0
    do i = 1, n
        if (rel(i)) then
            v = v + x(:, i)
            j = j + 1
        end if
    end do
    v = v/(j*1.0_flyt)
    do i = 1, n
        if (rel(i)) then
            x(:, i) = x(:, i) - v
            x(:, i) = x(:, i)/norm2(x(:, i))
        end if
    end do
end subroutine

!> Print a lot of configurations
subroutine dump_configurations(mag, ss)
    !> shells and stuff
    class(lo_magdisorder), intent(in) :: mag
    type(lo_crystalstructure), intent(in) :: ss
    !
    real(flyt), dimension(ss%na) :: moment
    real(flyt) :: f0
    integer :: i, ii, k, l, m, ncf, u
    character(len=10000) :: dum, fn

    printbasic: block
        do i = 1, ss%na
            if (mag%coll) then
                moment(i) = abs(ss%mag%collinear_moment(i))
            else
                moment(i) = norm2(ss%mag%noncollinear_moment(:, i))
            end if
        end do

        ! This is just the normal, random dumps where each configuration is completely uncorrelated
        ncf = mag%nbin + 1
        do i = 1, ncf
            fn = 'outfile.magmom_'//tochar(i)
            u = open_file('out', trim(fn))
            do k = 1, mag%nconf
                dum = "MAGMOM = "
                if (mag%coll) then
                    do l = 1, size(mag%collhistory, 1)
                        if (i .eq. 1) then
                            ii = mag%initial_collinear_configuration(l)
                        else
                            ii = mag%collhistory(l, k, i - 1)
                        end if
                        ii = ii*int(anint(moment(l)))
                        dum = trim(dum)//" "//tochar(ii)
                    end do
                else
                    do l = 1, size(mag%noncollhistory, 2)
                    do m = 1, 3
                        if (i .eq. 1) then
                            f0 = mag%initial_noncollinear_configuration(m, l)
                        else
                            f0 = mag%noncollhistory(m, l, k, i - 1)
                        end if
                        f0 = f0*moment(l)
                        dum = trim(dum)//" "//tochar(f0)
                    end do
                    end do
                end if
                write (u, '(1X,A)') trim(dum)
            end do
            close (u)
        end do
    end block printbasic

!    if ( mag%coll .eqv. .false. ) then
!    printfancy: block
!        integer, parameter :: nflip=12
!        real(flyt), dimension(:,:), allocatable :: baseline,coeffM
!        real(flyt), dimension(:,:,:), allocatable :: newconf,oldconf
!        real(flyt), dimension(:), allocatable :: Cline,singval
!        real(flyt), dimension(3) :: v
!        real(flyt) :: cn1,cn2
!        integer, dimension(:), allocatable :: magatoms
!        integer, dimension(nflip) :: flipatoms
!        integer :: iter
!        ! Now do it slightly more sophisticated: In each bin I have decent configurations. Pick one of these
!        ! And then flip some spins 180 and 90 degrees to yield a better-conditioned problem to start from?
!        ! I only care about the exchange parameters of the magnetic ions, I think. So perhaps filter them a
!        ! little.
!        lo_allocate(di(map%ntheta_magpair))
!        di=0
!        do a1=1,map%nuc
!        do i=1,map%uc(a1)%nmagpair
!            a2=map%uc(a1)%magpair(i)%i2
!            if ( mag%coll ) then
!                f0=abs(uc%mag%collinear_moment(a1))
!                f1=abs(uc%mag%collinear_moment(a2))
!            else
!                f0=norm2(uc%mag%noncollinear_moment(:,a1))
!                f1=norm2(uc%mag%noncollinear_moment(:,a2))
!            endif
!            sh=map%uc(a1)%magpair(i)%irreducible_shell
!            if ( f0 .gt. lo_tol .and. f1 .gt. lo_tol ) then
!                do j=1,map%magpairshell(sh)%ntheta
!                    k=map%magpairshell(sh)%thetaind(j)
!                    if ( k .gt. 0 ) di(k)=1
!                enddo
!            endif
!        enddo
!        enddo
!        ntheta=sum(di)
!        lo_allocate(relirr(ntheta))
!        j=0
!        do i=1,size(di)
!            if ( di(i) .gt. 0 ) then
!                j=j+1
!                relirr(j)=i
!            endif
!        enddo
!
!        ! Get a list of the magnetic atoms
!        j=0
!        do i=1,ss%na
!            if ( norm2(ss%mag%noncollinear_moment(:,i)) .gt. lo_tol ) j=j+1
!        enddo
!        lo_allocate(magatoms(j))
!        j=0
!        do i=1,ss%na
!            if ( norm2(ss%mag%noncollinear_moment(:,i)) .gt. lo_tol ) then
!                j=j+1
!                magatoms(j)=i
!            endif
!        enddo
!
!        lo_allocate(baseline(3,ss%na))
!        lo_allocate(newconf(3,ss%na,nsubconf))
!        lo_allocate(oldconf(3,ss%na,nsubconf))
!        lo_allocate(coeffM(nsubconf-1,ntheta))
!        lo_allocate(Cline(map%ntheta_magpair))
!        lo_allocate(singval(ntheta))
!
!        do bin=1,1 !ncf
!        do conf=1,1 !mag%nconf
!            if ( bin .eq. 1 ) then
!                baseline=mag%initial_noncollinear_configuration
!            else
!                baseline=mag%noncollhistory(:,:,conf,bin-1)
!            endif
!
!            ! Grab the baseline distribution for this configuration
!            do i=1,nsubconf
!                newconf(:,:,i)=baseline
!            enddo
!            oldconf=newconf
!            cn1=lo_huge
!            cn2=0.0_flyt
!            do iter=1,5000
!
!                do i=2,nsubconf
!                    ! Pick some atoms to flip
!                    newconf(:,:,i)=baseline
!                    call randomsubset(magatoms,flipatoms)
!                    do j=1,nflip
!                        newconf(:,flipatoms(j),i)=random_unit_vector()
!                    enddo
!                    call lo_coeffmatrix_magpair_red(newconf(:,:,i),baseline,Cline,map)
!                    coeffM(i-1,:)=Cline(relirr)
!                enddo
!                ! SVD this matrix
!                call lo_dgesvd(coeffM,singval)
!                cn2=abs(maxval(singval)/minval(singval))
!write(*,*) iter,singval,cn1
!                if ( cn2 .lt. cn1 ) then
!                    oldconf=newconf
!                    cn1=cn2
!                endif
!                !
!!write(*,*) 'iter:',iter,'ss',tochar(flipatoms)
!            enddo
!
!        enddo
!        enddo
!    end block printfancy
!    endif
!
!    contains
!
!    subroutine randomsubset(original,subset)
!    !> list to choose from
!    integer, dimension(:), intent(in) :: original
!    !> list to be returned
!    integer, dimension(:), intent(out) :: subset
!    !
!    real(flyt) :: f0,f1
!    integer :: i,j,n_original,n_remaining,n_needed
!    !
!    n_original=size(original,1)
!    n_needed=size(subset,1)
!    ! sanity check
!    n_needed=max(n_needed,1)
!    n_needed=min(n_original,n_needed)
!    !
!    n_remaining=n_original
!    j=0
!    do i=1,n_original
!        f0=(1.0_flyt*n_needed)/(1.0_flyt*n_remaining)
!        call random_number(f1)
!        if ( f1 .le. f0 ) then
!            j=j+1
!            subset(j)=original(i)
!            n_needed=n_needed-1
!        endif
!        n_remaining=n_remaining-1
!        if ( n_needed .eq. 0 ) exit
!    enddo
!    end subroutine

end subroutine

!> calculate the correlation function
subroutine correlation_function(mag, collconf, noncollconf, cf)
    !> list of shells and stuff
    class(lo_magdisorder), intent(in) :: mag
    !> current collinear magnetic configuration
    integer, dimension(:), intent(in), optional :: collconf
    !> current noncollinear magnetic configuration
    real(flyt), dimension(:, :), intent(in), optional :: noncollconf
    !> the correlation function per shell
    real(flyt), dimension(:), intent(out) :: cf
    !
    integer :: i, j, l, i1, i2
    real(flyt) :: f0
    !
    if (mag%coll) then
        cf = 0.0_flyt
        do i = 1, mag%nshell
            l = 0
            do j = 1, mag%sh(i)%npair
                i1 = mag%sh(i)%i1(j)
                i2 = mag%sh(i)%i2(j)
                l = l + collconf(i1)*collconf(i2)
            end do
            cf(i) = (l*1.0_flyt)/(mag%sh(i)%npair*1.0_flyt)
        end do
    else
        cf = 0.0_flyt
        do i = 1, mag%nshell
            f0 = 0.0_flyt
            do j = 1, mag%sh(i)%npair
                i1 = mag%sh(i)%i1(j)
                i2 = mag%sh(i)%i2(j)
                f0 = f0 + dot_product(noncollconf(:, i1), noncollconf(:, i2))
            end do
            cf(i) = (f0*1.0_flyt)/(mag%sh(i)%npair*1.0_flyt)
        end do
    end if
end subroutine

!> generate optimized configuration
subroutine optimize(mag, ss, nconf, nbin)
    !> shells and stuff
    class(lo_magdisorder), intent(inout) :: mag
    !> crystal structure
    type(lo_crystalstructure), intent(in) :: ss
    !> how many configurations do I want in each bin?
    integer, intent(in) :: nconf
    !> how many bins?
    integer, intent(in) :: nbin

    type(lo_mersennetwister) :: tw
    real(flyt), dimension(mag%nshell) :: cf
    real(flyt), dimension(:, :), allocatable :: cftargets
    real(flyt), dimension(:), allocatable :: cffactor
    integer :: bin, i, j

    ! Seed random numbers
    call tw%init(iseed=0, rseed=walltime())

    ! Space for the history
    mag%nbin = nbin
    mag%nconf = nconf
    if (mag%coll) then
        lo_allocate(mag%collhistory(ss%na, nconf, nbin))
        mag%collhistory = 0
    else
        lo_allocate(mag%noncollhistory(3, ss%na, nconf, nbin))
        mag%noncollhistory = 0.0_flyt
    end if

    ! Set the target correlation functions
    lo_allocate(cftargets(mag%nshell, mag%nbin))
    lo_allocate(cffactor(mag%nbin))
    call mag%correlation_function(mag%initial_collinear_configuration, mag%initial_noncollinear_configuration, cf)
    cftargets = 0.0_flyt
    do i = 1, size(cftargets, 2)
        cffactor(i) = lo_chop(abs((1.0_flyt - (i)/(1.0_flyt*mag%nbin))), lo_sqtol)
        cftargets(:, i) = cf*cffactor(i)
        cffactor(i) = (cffactor(i) + 2.0_flyt)/3.0_flyt
    end do

    ! I have a series of correlation function targets, one minimization for each:
    do bin = 1, mag%nbin
        findonetarget: block
            real(flyt), dimension(:), allocatable :: convcheck
            real(flyt), dimension(3, 3) :: m0
            real(flyt), dimension(3) :: v0
            real(flyt) :: cf0, cf1, f0, f1
            real(flyt) :: temperature, tempinc, tempdec, breaktol
            real(flyt), dimension(:, :, :), allocatable :: noncollhist
            real(flyt), dimension(3, ss%na) :: ncconf0, ncconf1
            integer, dimension(:, :), allocatable :: collhist
            integer, dimension(ss%na) :: conf0, conf1
            integer :: nouter, ninner, na, histcounter
            integer :: outiter, initer, nflip

            cf = 0.0_flyt
            ! Initial configuration
            if (mag%coll) then
                conf0 = mag%initial_collinear_configuration
                conf1 = 0
                call mag%correlation_function(collconf=conf0, cf=cf)
            else
                ncconf0 = mag%initial_noncollinear_configuration
                ncconf1 = 0
                call mag%correlation_function(noncollconf=ncconf0, cf=cf)
            end if
            cf0 = sum(abs(cf - cftargets(:, bin)))/mag%nshell
            cf1 = 0.0_flyt
            na = size(mag%sites, 1)
            ! Some counters for the minimization
            nouter = 100
            ninner = mag%nconf*ss%na
            temperature = 0.3_flyt
            tempinc = 1.5_flyt
            tempdec = 0.5_flyt
            lo_allocate(convcheck(nouter))
            convcheck = 0.0_flyt
            breaktol = 1E-2_flyt

            write (*, *) 'Simulated annealing, bin #'//tochar(bin)

            ! Start minimizing
            outerloop1: do outiter = 1, nouter
                nflip = 0
                do initer = 1, ninner
                    ! Flip a spin, get new correlation function
                    if (mag%coll) then
                        conf1 = conf0
                        if (mag%ferromagnetic) then
                            ! Just flip a random spin
                            i = tw%rnd_int(na)
                            conf1(mag%sites(i)) = -1*conf1(mag%sites(i))
                            call mag%correlation_function(collconf=conf1, cf=cf)
                        else
                            ! Change place of two spins so that the total moment is preserved
                            do
                                i = tw%rnd_int(na)
                                j = tw%rnd_int(na)
                                if (conf1(mag%sites(i))*conf1(mag%sites(j)) .lt. 0) exit
                            end do
                            conf1(mag%sites(i)) = -1*conf1(mag%sites(i))
                            conf1(mag%sites(j)) = -1*conf1(mag%sites(j))
                            call mag%correlation_function(collconf=conf1, cf=cf)
                        end if
                    else
                        ncconf1 = ncconf0
                        i = tw%rnd_int(na)
                        if (mag%ferromagnetic) then
                            v0 = tw%rnd_unitvector()*(1.0_flyt - cffactor(bin)) + cffactor(bin)*ncconf1(:, mag%sites(i))
                            v0 = v0/norm2(v0)
                            ncconf1(:, mag%sites(i)) = v0 !random_unit_vector()
                        else
                            v0 = tw%rnd_unitvector()*(1.0_flyt - cffactor(bin)) + cffactor(bin)*ncconf1(:, mag%sites(i))
                            v0 = v0/norm2(v0)
                            ncconf1(:, mag%sites(i)) = v0 !random_unit_vector()
                            call zerosum(ncconf1, ss%mag%atom_has_moment)
                        end if
                        call mag%correlation_function(noncollconf=ncconf1, cf=cf)
                    end if
                    cf1 = sum(abs(cf - cftargets(:, bin)))/mag%nshell
                    ! Add a small bias to ferromagnetic?

                    ! MC compare thingy
                    f0 = exp(-(cf1 - cf0)/temperature)
                    f1 = tw%rnd_real()
                    ! keep?
                    if (f0 .gt. f1) then
                        if (mag%coll) then
                            conf0 = conf1
                        else
                            ncconf0 = ncconf1
                        end if
                        cf0 = cf1
                        nflip = nflip + 1
                    end if
                end do

                if (mag%coll .eqv. .false. .and. mag%ferromagnetic) then
                    v0 = 0.0_flyt
                    do i = 1, na
                        v0 = v0 + ncconf0(:, i)
                    end do
                    v0 = v0/norm2(v0)
                    m0 = lo_rotation_matrix_from_vector_a_to_b(v0, [1.0_flyt, 0.0_flyt, 0.0_flyt])
                    do i = 1, na
                        ncconf0(:, i) = matmul(m0, ncconf0(:, i))
                    end do
                end if

                ! check how many flips there were, and maybe adjust the temperature
                f0 = (1.0_flyt*nflip)/(1.0_flyt*ninner)
                if (f0 .lt. 0.02_flyt) then
                    temperature = temperature*tempinc
                elseif (f0 .gt. 0.10_flyt) then
                    temperature = temperature*tempdec
                end if
                ! check for convergence
                convcheck(outiter) = cf0
                if (outiter .ge. 5) then
                    f1 = lo_mean(convcheck(outiter - 4:outiter))
                else
                    f1 = 1.0_flyt
                end if
                if (f1 .lt. breaktol) exit outerloop1
                !
                if (mag%coll) then
                    write (*, '(1X,I8,1X,F10.7,4(1X,F10.7))') outiter, cf0, temperature, real(f0), f1, sum(conf0*1.0_flyt/na)
                else
                    write (*, '(1X,I8,1X,F10.7,6(1X,F10.7))') outiter, cf0, temperature, real(f0), f1, &
                        sum(ncconf0(1, :)*1.0_flyt/na), sum(ncconf0(2, :)*1.0_flyt/na), sum(ncconf0(3, :)*1.0_flyt/na)
                end if
            end do outerloop1

            convcheck = 0.0_flyt
            histcounter = 0
            if (mag%coll) then
                lo_allocate(collhist(ss%na, mag%nconf))
                collhist = 0
            else
                lo_allocate(noncollhist(3, ss%na, mag%nconf))
                noncollhist = 0.0_flyt
            end if
            ! Now gather some statistics for the history
            outerloop2: do outiter = 1, nouter
                nflip = 0
                do initer = 1, ninner
                    ! Flip a spin, get new correlation function
                    if (mag%coll) then
                        conf1 = conf0
                        if (mag%ferromagnetic) then
                            ! Just flip a random spin
                            i = tw%rnd_int(na)
                            conf1(mag%sites(i)) = -1*conf1(mag%sites(i))
                            call mag%correlation_function(collconf=conf1, cf=cf)
                        else
                            ! Change place of two spins so that the total moment is preserved
                            do
                                i = tw%rnd_int(na)
                                j = tw%rnd_int(na)
                                if (conf1(mag%sites(i))*conf1(mag%sites(j)) .lt. 0) exit
                            end do
                            conf1(mag%sites(i)) = -1*conf1(mag%sites(i))
                            conf1(mag%sites(j)) = -1*conf1(mag%sites(j))
                            call mag%correlation_function(collconf=conf1, cf=cf)
                        end if
                    else
                        ncconf1 = ncconf0
                        i = tw%rnd_int(na)
                        if (mag%ferromagnetic) then
                            v0 = tw%rnd_unitvector()*(1.0_flyt - cffactor(bin)) + cffactor(bin)*ncconf1(:, mag%sites(i))
                            v0 = v0/norm2(v0)
                            ncconf1(:, mag%sites(i)) = v0 !random_unit_vector()
                        else
                            v0 = tw%rnd_unitvector()*(1.0_flyt - cffactor(bin)) + cffactor(bin)*ncconf1(:, mag%sites(i))
                            v0 = v0/norm2(v0)
                            ncconf1(:, mag%sites(i)) = v0 !random_unit_vector()
                            call zerosum(ncconf1, ss%mag%atom_has_moment)
                        end if
                        call mag%correlation_function(noncollconf=ncconf1, cf=cf)
                    end if
                    cf1 = sum(abs(cf - cftargets(:, bin)))/mag%nshell
                    ! MC compare thingy
                    f0 = exp(-(cf1 - cf0)/temperature)
                    f1 = tw%rnd_real()
                    ! keep?
                    if (f0 .gt. f1) then
                        if (mag%coll) then
                            conf0 = conf1
                            ! Is this a new configuration?
                            j = 0
                            do i = 1, histcounter
                                if (sum(abs(conf1 - collhist(:, i))) .le. 2) then
                                    j = j + 1
                                    exit
                                end if
                            end do
                            ! If it is new, keep it
                            if (j .eq. 0) then
                                histcounter = histcounter + 1
                                collhist(:, histcounter) = conf1
                                ! We might be done now!
                                if (histcounter .eq. mag%nconf) exit outerloop2
                            end if
                        else
                            ncconf0 = ncconf1
                            ! Is this a new configuration?
                            j = 0
                            do i = 1, histcounter
                                if (sum(abs(ncconf1 - noncollhist(:, :, i))) .le. 0.2_flyt) then
                                    j = j + 1
                                    exit
                                end if
                            end do
                            ! If it is new, keep it!
                            if (j .eq. 0) then
                                histcounter = histcounter + 1
                                noncollhist(:, :, histcounter) = ncconf1
                                ! We might be done now!
                                if (histcounter .eq. mag%nconf) exit outerloop2
                            end if
                        end if
                        !
                        cf0 = cf1
                        nflip = nflip + 1
                        !
                    end if
                end do
                ! check how many flips there were, and maybe adjust the temperature
                f0 = (1.0_flyt*nflip)/(1.0_flyt*ninner)
                if (f0 .lt. 0.02_flyt) then
                    temperature = temperature*tempinc
                elseif (f0 .gt. 0.10_flyt) then
                    temperature = temperature*tempdec
                end if
            end do outerloop2
            ! And store the history
            if (mag%coll) then
                mag%collhistory(:, :, bin) = collhist
            else
                ! Rotate all spins such that the average points in the x-direction?
                do j = 1, mag%nconf
                    v0 = 0.0_flyt
                    do i = 1, na
                        v0 = v0 + noncollhist(:, i, j)
                    end do
                    v0 = v0/norm2(v0)
                    m0 = lo_rotation_matrix_from_vector_a_to_b(v0, [1.0_flyt, 0.0_flyt, 0.0_flyt])
                    do i = 1, na
                        noncollhist(:, i, j) = matmul(m0, noncollhist(:, i, j))
                    end do
                end do
                mag%noncollhistory(:, :, :, bin) = noncollhist
            end if
        end block findonetarget
    end do

end subroutine

!> set up all the coordination shells and stuff
subroutine generate(mag, sl, uc, ss)
    !> to keep track of all the coordination shells
    class(lo_magdisorder), intent(out) :: mag
    !> a symmetry table, useful for a bunch of stuff
    type(lo_symlist), intent(in) :: sl
    !> unit cell
    type(lo_crystalstructure), intent(in) :: uc
    !> supercell
    type(lo_crystalstructure), intent(in) :: ss

    ! Figure out wether it is ferro or antoferromagnetic!
    fmorafm: block
        integer :: i
        real(flyt), dimension(3) :: v0, v1
        real(flyt) :: f0

        write (*, *) '... figuring out if it is ferro- or antiferromagnetic.'
        f0 = 0.0_flyt
        v0 = 0.0_flyt
        v1 = 0.0_flyt
        do i = 1, uc%na
            ! Only bother with atoms with moments
            if (uc%mag%atom_has_moment(i)) then
                if (uc%info%collmag) then
                    v0 = v0 + [1, 1, 1]*uc%mag%collinear_moment(i)
                    v1 = v1 + [1, 1, 1]*abs(uc%mag%collinear_moment(i))
                else
                    v0 = v0 + uc%mag%noncollinear_moment(:, i)
                    v1 = v1 + abs(uc%mag%noncollinear_moment(:, i))
                end if
            end if
        end do
        f0 = sum(abs(abs(v0) - abs(v1)))
        if (f0 .lt. lo_tol) then
            mag%ferromagnetic = .true.
            write (*, *) '... suppose it is ferromagnetic'
        else
            mag%ferromagnetic = .false.
            write (*, *) '... suppose it is antiferromagnetic'
        end if
    end block fmorafm

    ! Count unique atoms with magnetic moments
    findrelevantshells: block
        integer, dimension(:), allocatable :: shelli
        integer :: i, j, k, l
        logical, dimension(:), allocatable :: shellmag

        ! Count number of relevant coordination shells, and create an
        ! index that maps from shell index to magnetic shell index
        lo_allocate(shelli(sl%npairshells))
        lo_allocate(shellmag(sl%npairshells))
        shelli = -1
        shellmag = .false.
        l = 0
        do i = 1, sl%npairshells
            if (lo_sqnorm(sl%pairshell(i)%protpair%v) .lt. lo_sqtol) cycle
            if (uc%mag%atom_has_moment(sl%pairshell(i)%protpair%i1) .and. &
                uc%mag%atom_has_moment(sl%pairshell(i)%protpair%i2)) then
                l = l + 1
                shelli(i) = l
                shellmag(i) = .true.
            end if
        end do

        ! initialize the shells
        mag%nshell = l
        lo_allocate(mag%sh(mag%nshell))
        do i = 1, mag%nshell
            mag%sh(i)%npair = 0
        end do

        ! Count pairs per shell
        do i = 1, sl%nss
        do j = 1, sl%ss(i)%npair
            k = sl%ss(i)%pair(j)%unique_shell
            if (shellmag(k) .eqv. .false.) cycle
            if (sl%ss(i)%pair(j)%j1 < sl%ss(i)%pair(j)%j2) cycle
            l = shelli(k)
            mag%sh(l)%npair = mag%sh(l)%npair + 1
        end do
        end do

        ! make some space
        do i = 1, mag%nshell
            lo_allocate(mag%sh(i)%i1(mag%sh(i)%npair))
            lo_allocate(mag%sh(i)%i2(mag%sh(i)%npair))
            mag%sh(i)%i1 = 0
            mag%sh(i)%i2 = 0
            mag%sh(i)%npair = 0
        end do

        ! now store all the pairs for each shell
        do i = 1, sl%nss
        do j = 1, sl%ss(i)%npair
            k = sl%ss(i)%pair(j)%unique_shell
            if (shellmag(k) .eqv. .false.) cycle
            if (sl%ss(i)%pair(j)%j1 < sl%ss(i)%pair(j)%j2) cycle
            l = shelli(k)
            mag%sh(l)%npair = mag%sh(l)%npair + 1
            mag%sh(l)%i1(mag%sh(l)%npair) = sl%ss(i)%pair(j)%j1
            mag%sh(l)%i2(mag%sh(l)%npair) = sl%ss(i)%pair(j)%j2
        end do
        end do

    end block findrelevantshells
    write (*, *) '... identified magnetic coordination shells'

    ! Now set up the starting configuration
    setupconf: block
        real(flyt), dimension(:), allocatable :: cf
        integer :: i, j

        ! collinear or noncollinear?
        if (uc%info%collmag) then
            mag%coll = .true.
        else
            mag%coll = .false.
        end if

        ! Get a list of switchable sites
        j = 0
        do i = 1, ss%na
            if (ss%mag%atom_has_moment(i)) then
                j = j + 1
            end if
        end do
        lo_allocate(mag%sites(j))
        j = 0
        do i = 1, ss%na
            if (ss%mag%atom_has_moment(i)) then
                j = j + 1
                mag%sites(j) = i
            end if
        end do

        ! Set up the initial configuration
        if (mag%coll) then
            lo_allocate(mag%initial_collinear_configuration(ss%na))
            mag%initial_collinear_configuration = 0
            do i = 1, ss%na
                if (ss%mag%atom_has_moment(i)) then
                    mag%initial_collinear_configuration(i) = int(anint(ss%mag%collinear_moment(i)/abs(ss%mag%collinear_moment(i))))
                end if
            end do
        else
            lo_allocate(mag%initial_noncollinear_configuration(3, ss%na))
            mag%initial_noncollinear_configuration = 0.0_flyt
            do i = 1, ss%na
                if (ss%mag%atom_has_moment(i)) then
                    mag%initial_noncollinear_configuration(:, i) = ss%mag%noncollinear_moment(:, i)/norm2(ss%mag%noncollinear_moment(:, i))
                end if
            end do
        end if

        ! Get the reference correlation function
        lo_allocate(cf(mag%nshell))
        if (mag%coll) then
            call mag%correlation_function(collconf=mag%initial_collinear_configuration, cf=cf)
        else
            call mag%correlation_function(noncollconf=mag%initial_noncollinear_configuration, cf=cf)
        end if

    end block setupconf
    write (*, *) '... generated initial magnetic configuration'

end subroutine

end module
