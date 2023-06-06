module lo_dielectric_interaction
!!
!! Contains the higher order interactions between electric fields and
!! atomic positions, i.e. the higher order Born charges and dielectric
!! constants.
!!
use konstanter, only: r8,i8,lo_iou,lo_huge,lo_hugeint,lo_bohr_to_A,lo_exitcode_param,lo_sqtol
use gottochblandat, only: open_file,tochar,lo_chop,lo_sqnorm,lo_clean_fractional_coordinates,walltime
use type_crystalstructure, only: lo_crystalstructure
use mpi_wrappers, only: lo_mpi_helper,lo_stop_gracefully
use lo_memtracker, only: lo_mem_helper
use lo_verletboxes, only: lo_verletbox
implicit none

private
public :: lo_dielectric_tensor

!> dielectric singlet
type eps_singlet
    !> interaction tensor (E,E,u)
    real(r8), dimension(3,3,3) :: m=-lo_huge
end type

!> dielectric pair
type eps_pair
    !> index to first atom
    integer :: a1=-lo_hugeint
    !> index to second atom
    integer :: a2=-lo_hugeint
    !> lattice vector between atoms
    real(r8), dimension(3) :: lv=-lo_huge
    !> absolute vector between atoms
    real(r8), dimension(3) :: v=-lo_huge
    !> interaction tensor (E,E,u,u)
    real(r8), dimension(3,3,3,3) :: m=-lo_huge
end type

!> second order born charge
type Z_pair
    !> index to first atom
    integer :: a1=-lo_hugeint
    !> index to second atom
    integer :: a2=-lo_hugeint
    !> lattice vector between atoms
    real(r8), dimension(3) :: lv=-lo_huge
    !> absolute vector between atoms
    real(r8), dimension(3) :: v=-lo_huge
    !> interaction tensor (E,u,u)
    real(r8), dimension(3,3,3) :: m=-lo_huge
end type

!> second order born charge
type Z_triplet
    !> index to first atom
    integer :: a1=-lo_hugeint
    !> index to second atom
    integer :: a2=-lo_hugeint
    !> index to third atom
    integer :: a3=-lo_hugeint
    !> lattice vector between atoms
    real(r8), dimension(3) :: lv2=-lo_huge
    real(r8), dimension(3) :: lv3=-lo_huge
    !> absolute vector between atoms
    real(r8), dimension(3) :: v2=-lo_huge
    real(r8), dimension(3) :: v3=-lo_huge
    !> interaction tensor (E,u,u,u)
    real(r8), dimension(3,3,3,3) :: m=-lo_huge
end type

!> dielectric interactions (to many orders)
type lo_dielectric_tensor
    !> number of atoms
    integer :: n_atom=-lo_hugeint

    !> number of born charges
    integer :: n_Z_singlet=-lo_hugeint
    integer :: n_Z_pair=-lo_hugeint
    integer :: n_Z_triplet=-lo_hugeint
    !> born charges
    real(r8), dimension(:,:,:), allocatable :: Z_singlet
    type(Z_pair), dimension(:), allocatable :: Z_pair
    type(Z_triplet), dimension(:), allocatable :: Z_triplet

    !> number of epsilon singlets
    integer :: n_eps_singlet=-lo_hugeint
    integer :: n_eps_pair=-lo_hugeint
    !> dielectric singlets
    real(r8), dimension(3,3) :: eps_inf=-lo_huge
    type(eps_singlet), dimension(:), allocatable :: eps_singlet
    type(eps_pair), dimension(:), allocatable :: eps_pair
    contains
        !> read from file
        procedure :: readfromfile
        !> write to file
        procedure :: writetofile
        !> remap to another unitcell
        procedure :: remap
        !> size in memory, in bytes
        procedure :: size_in_mem
        !> destroy
        procedure :: destroy
        !> matrix elements, various ways
        procedure :: epsilon_matrixelement_firstorder
        procedure :: epsilon_matrixelement_secondorder
        procedure :: polarization_matrixelement_firstorder
        procedure :: polarization_matrixelement_secondorder
        procedure :: pretransform_Z_secondorder
        procedure :: pretransform_epsilon_secondorder
end type

!> interfaces to matrix elements
interface
    module function epsilon_matrixelement_firstorder(di,nuvector) result(M)
        class(lo_dielectric_tensor), intent(in) :: di
        complex(r8), dimension(:), intent(in) :: nuvector
        complex(r8), dimension(3,3) :: M
    end function
    module function epsilon_matrixelement_secondorder(di,q,nuvector1,nuvector2) result(M)
        class(lo_dielectric_tensor), intent(in) :: di
        real(r8), dimension(3) :: q
        complex(r8), dimension(:), intent(in) :: nuvector1,nuvector2
        complex(r8), dimension(3,3) :: M
    end function
    module function polarization_matrixelement_firstorder(di,nuvector) result(M)
        class(lo_dielectric_tensor), intent(in) :: di
        complex(r8), dimension(:), intent(in) :: nuvector
        complex(r8), dimension(3) :: M
    end function
    module function polarization_matrixelement_secondorder(di,q,nuvector1,nuvector2) result(M)
        class(lo_dielectric_tensor), intent(in) :: di
        real(r8), dimension(3) :: q
        complex(r8), dimension(:), intent(in) :: nuvector1,nuvector2
        complex(r8), dimension(3) :: M
    end function
    module function polarization_matrixelement_thirdorder(di,q2,q3,nuvector1,nuvector2,nuvector3) result(M)
        type(lo_dielectric_tensor), intent(in) :: di
        real(r8), dimension(3) :: q2,q3
        complex(r8), dimension(:), intent(in) :: nuvector1,nuvector2,nuvector3
        complex(r8), dimension(3) :: M
    end function
    module subroutine pretransform_Z_secondorder(di,q,ptf)
        class(lo_dielectric_tensor), intent(in) :: di
        real(r8), dimension(3) :: q
        complex(r8), dimension(:,:), intent(out) :: ptf
    end subroutine
    module subroutine pretransform_epsilon_secondorder(di,q,ptf)
        class(lo_dielectric_tensor), intent(in) :: di
        real(r8), dimension(3) :: q
        complex(r8), dimension(:,:), intent(out) :: ptf
    end subroutine
end interface

contains

!> rearrange to a different unit cell
subroutine remap(di,dj,uc,ss,mw,mem,verbosity)
    !> unitcell interactions
    class(lo_dielectric_tensor), intent(in) :: di
    !> supercell interactions
    class(lo_dielectric_tensor), intent(out) :: dj
    !> unitcell
    type(lo_crystalstructure), intent(in) :: uc
    !> supercell
    type(lo_crystalstructure), intent(in) :: ss
    !> mpi helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem
    !> talk a lot?
    integer, intent(in) :: verbosity

    type(lo_verletbox) :: vb
    integer :: ncell
    real(r8) :: timer,t0,t1

    timer=walltime()
    t0=timer
    t1=timer

    init: block
        ! Some heuristics
        if ( ss%info%supercell .eqv. .false. ) then
            call lo_stop_gracefully(['The supercell needs to be related to the unitcell'],lo_exitcode_param,__FILE__,__LINE__)
        endif

        if ( verbosity .gt. 0 ) then
            write(lo_iou,*) ''
            write(lo_iou,*) 'REMAPPING DIELECTRIC RESPONSE'
        endif

        ! How many unitcells is the supercell?
        ncell=ss%na/uc%na

        ! How I can figure out the number of tuplets of each kind
        dj%n_atom        = ss%na
        dj%n_Z_singlet   = di%n_Z_singlet   * ncell
        dj%n_Z_pair      = di%n_Z_pair      * ncell
        dj%n_Z_triplet   = di%n_Z_triplet   * ncell
        dj%n_eps_singlet = di%n_eps_singlet * ncell
        dj%n_eps_pair    = di%n_eps_pair    * ncell

        ! Copy dielectric tensor
        dj%eps_inf       = di%eps_inf

        ! Generate Verlet boxes for fast lookup of where atoms are
        call vb%generate(ss%r,[11,11,11],mem)

        if ( verbosity .gt. 0 ) then
            write(*,*) 'dj%n_atom       ',dj%n_atom
            write(*,*) 'dj%n_Z_singlet  ',dj%n_Z_singlet
            write(*,*) 'dj%n_Z_pair     ',dj%n_Z_pair
            write(*,*) 'dj%n_Z_triplet  ',dj%n_Z_triplet
            write(*,*) 'dj%n_eps_singlet',dj%n_eps_singlet
            write(*,*) 'dj%n_eps_pair   ',dj%n_eps_pair
            t1=walltime()
            write(*,*) '... init remapping:',t1-t0
            t0=t1
        endif
    end block init

    if ( dj%n_Z_singlet .gt. 0 ) then
    Zsing: block
        integer :: is,js

        allocate(dj%Z_singlet(3,3,ss%na))
        dj%Z_singlet=0.0_r8
        do js=1,ss%na
            is=ss%info%index_in_unitcell(js)
            dj%Z_singlet(:,:,js)=di%Z_singlet(:,:,is)
        enddo
    end block Zsing
    endif

    if ( dj%n_Z_pair .gt. 0 ) then
    Zpair: block
        real(r8), dimension(3) :: v1,v2
        integer, dimension(:,:), allocatable :: pair_per_uca,ind
        integer, dimension(:), allocatable :: npair_per_uca,offset_per_ssa
        integer :: ip,jp,a1,a2,a3

        ! Count pairs per atom
        call mem%allocate(npair_per_uca,uc%na,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        npair_per_uca=0
        do ip=1,di%n_Z_pair
            if ( mod(ip,mw%n) .ne. mw%r ) cycle
            a1=di%Z_pair(ip)%a1
            npair_per_uca(a1)=npair_per_uca(a1)+1
        enddo
        call mw%allreduce('sum',npair_per_uca)

        ! Store pairs per atom
        call mem%allocate(pair_per_uca,[maxval(npair_per_uca),uc%na],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        pair_per_uca=0
        npair_per_uca=0
        do ip=1,di%n_Z_pair
            a1=di%Z_pair(ip)%a1
            npair_per_uca(a1)=npair_per_uca(a1)+1
            pair_per_uca(npair_per_uca(a1),a1)=ip
        enddo

        ! Now start mapping for reals
        call mem%allocate(offset_per_ssa,ss%na,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%allocate(ind,[3,dj%n_Z_pair],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        offset_per_ssa=0
        ind=0

        ! Well not just yet, I need the offset in index per supercell pair to
        ! make it parallel and neat. Also gets the pairs sorted somewhat.
        ip=0
        do a1=1,ss%na
            offset_per_ssa(a1)=ip
            a2=ss%info%index_in_unitcell(a1)
            ip=ip+npair_per_uca(a2)
        enddo

        ! Start locating pairs?
        do a1=1,ss%na
            a2=ss%info%index_in_unitcell(a1)
            v1=ss%r(:,a1)
            do ip=1,npair_per_uca(a2)

                ! make it parallel?
                if ( mod(ip+offset_per_ssa(a1),mw%n) .ne. mw%r ) cycle
                jp=pair_per_uca(ip,a2)  ! index to unitcell pair

                v2=v1+matmul(ss%inv_latticevectors,di%Z_pair(jp)%v)
                v2=lo_clean_fractional_coordinates(v2)
                v2=lo_chop(v2,1E-11_r8)
                v2=lo_clean_fractional_coordinates(v2)
                ! Now I have the location of the second atom in fractional coordinates.
                a3=vb%locate(ss%r,v2)
                if ( a3 .gt. 0 ) then
                    ind(1,ip+offset_per_ssa(a1))=a1
                    ind(2,ip+offset_per_ssa(a1))=a3
                    ind(3,ip+offset_per_ssa(a1))=jp
                else
                    call lo_stop_gracefully(['Could not locate pair in supercell'],lo_exitcode_param,__FILE__,__LINE__)
                endif
            enddo
        enddo
        call mw%allreduce('max',ind)

        if ( verbosity .gt. 0 ) then
            t1=walltime()
            write(*,*) '... sorted out Z pairs (',tochar(t1-t0),'s)'
            t0=t1
        endif

        ! Now store tensors and things properly
        allocate(dj%Z_pair(dj%n_Z_pair))
        do ip=1,dj%n_Z_pair
            a1=ind(1,ip)
            a2=ind(2,ip)
            jp=ind(3,ip)
            dj%Z_pair(ip)%a1=a1
            dj%Z_pair(ip)%a2=a2
            dj%Z_pair(ip)%v=di%Z_pair(jp)%v
            ! v  = a2 + lv - a1
            ! lv = v  + a1 - a2
            v1=dj%Z_pair(ip)%v+ss%rcart(:,a1)-ss%rcart(:,a2)
            v1=matmul(ss%inv_latticevectors,v1)
            if ( sum(abs(v1-anint(v1))) .gt. lo_sqtol ) then
                call lo_stop_gracefully(['I do not understand vectors'],lo_exitcode_param,__FILE__,__LINE__)
            else
                v1=anint(v1)
                dj%Z_pair(ip)%lv=matmul(ss%latticevectors,v1)
            endif
            ! And the actual tensor
            dj%Z_pair(ip)%m=di%Z_pair(jp)%m
        enddo

        call mem%deallocate(npair_per_uca ,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%deallocate(pair_per_uca  ,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%deallocate(offset_per_ssa,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%deallocate(ind           ,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)

        if ( verbosity .gt. 0 ) then
            t1=walltime()
            write(*,*) '... stored Z pairs (',tochar(t1-t0),'s)'
            t0=t1
        endif
    end block Zpair
    endif

    if ( dj%n_Z_triplet .gt. 0 ) then
    Ztriplet: block
        real(r8), dimension(3) :: v1,v2,v3
        integer, dimension(:,:), allocatable :: triplet_per_uca,ind
        integer, dimension(:), allocatable :: ntriplet_per_uca,offset_per_ssa
        integer :: ip,jp,a1,a2,a3,uca

        ! Count triplets per atom
        call mem%allocate(ntriplet_per_uca,uc%na,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        ntriplet_per_uca=0
        do ip=1,di%n_Z_triplet
            if ( mod(ip,mw%n) .ne. mw%r ) cycle
            a1=di%Z_triplet(ip)%a1
            ntriplet_per_uca(a1)=ntriplet_per_uca(a1)+1
        enddo
        call mw%allreduce('sum',ntriplet_per_uca)

        ! Store triplets per atom
        call mem%allocate(triplet_per_uca,[maxval(ntriplet_per_uca),uc%na],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        triplet_per_uca=0
        ntriplet_per_uca=0
        do ip=1,di%n_Z_triplet
            a1=di%Z_triplet(ip)%a1
            ntriplet_per_uca(a1)=ntriplet_per_uca(a1)+1
            triplet_per_uca(ntriplet_per_uca(a1),a1)=ip
        enddo

        ! Now start mapping for reals
        call mem%allocate(offset_per_ssa,ss%na,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%allocate(ind,[4,dj%n_Z_triplet],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        offset_per_ssa=0
        ind=0

        ! Well not just yet, I need the offset in index per supercell triplet to
        ! make it parallel and neat. Also gets the triplets sorted somewhat.
        ip=0
        do a1=1,ss%na
            offset_per_ssa(a1)=ip
            a2=ss%info%index_in_unitcell(a1)
            ip=ip+ntriplet_per_uca(a2)
        enddo

        ! Start locating triplets?
        do a1=1,ss%na
            uca=ss%info%index_in_unitcell(a1)
            v1=ss%r(:,a1)
            do ip=1,ntriplet_per_uca(uca)

                ! make it parallel?
                if ( mod(ip+offset_per_ssa(a1),mw%n) .ne. mw%r ) cycle
                jp=triplet_per_uca(ip,uca)  ! index to unitcell triplet

                v2=v1+matmul(ss%inv_latticevectors,di%Z_triplet(jp)%v2)
                v2=lo_clean_fractional_coordinates(v2)
                v2=lo_chop(v2,1E-11_r8)
                v2=lo_clean_fractional_coordinates(v2)
                v3=v1+matmul(ss%inv_latticevectors,di%Z_triplet(jp)%v3)
                v3=lo_clean_fractional_coordinates(v3)
                v3=lo_chop(v3,1E-11_r8)
                v3=lo_clean_fractional_coordinates(v3)


                ! Now I have the location of the second atom in fractional coordinates.
                a2=vb%locate(ss%r,v2)
                a3=vb%locate(ss%r,v3)

                if ( a2 .gt. 0 .and. a3 .gt. 0 ) then
                    ind(1,ip+offset_per_ssa(a1))=a1
                    ind(2,ip+offset_per_ssa(a1))=a2
                    ind(3,ip+offset_per_ssa(a1))=a3
                    ind(4,ip+offset_per_ssa(a1))=jp
                else
                    call lo_stop_gracefully(['Could not locate triplet in supercell'],lo_exitcode_param,__FILE__,__LINE__)
                endif
            enddo
        enddo
        call mw%allreduce('max',ind)

        if ( verbosity .gt. 0 ) then
            t1=walltime()
            write(*,*) '... sorted out Z triplets (',tochar(t1-t0),'s)'
            t0=t1
        endif

        ! Now store tensors and things properly
        allocate(dj%Z_triplet(dj%n_Z_triplet))
        do ip=1,dj%n_Z_triplet
            a1=ind(1,ip)
            a2=ind(2,ip)
            a3=ind(3,ip)
            jp=ind(4,ip)
            dj%Z_triplet(ip)%a1=a1
            dj%Z_triplet(ip)%a2=a2
            dj%Z_triplet(ip)%a3=a3
            dj%Z_triplet(ip)%v2=di%Z_triplet(jp)%v2
            dj%Z_triplet(ip)%v3=di%Z_triplet(jp)%v3
            ! v  = a2 + lv - a1
            ! lv = v  + a1 - a2
            v1=dj%Z_triplet(ip)%v2+ss%rcart(:,a1)-ss%rcart(:,a2)
            v1=matmul(ss%inv_latticevectors,v1)
            if ( sum(abs(v1-anint(v1))) .gt. lo_sqtol ) then
                call lo_stop_gracefully(['I do not understand vectors'],lo_exitcode_param,__FILE__,__LINE__)
            else
                v1=anint(v1)
                dj%Z_triplet(ip)%lv2=matmul(ss%latticevectors,v1)
            endif
            v1=dj%Z_triplet(ip)%v3+ss%rcart(:,a1)-ss%rcart(:,a3)
            v1=matmul(ss%inv_latticevectors,v1)
            if ( sum(abs(v1-anint(v1))) .gt. lo_sqtol ) then
                call lo_stop_gracefully(['I do not understand vectors'],lo_exitcode_param,__FILE__,__LINE__)
            else
                v1=anint(v1)
                dj%Z_triplet(ip)%lv3=matmul(ss%latticevectors,v1)
            endif
            ! And the actual tensor
            dj%Z_triplet(ip)%m=di%Z_triplet(jp)%m
        enddo

        call mem%deallocate(ntriplet_per_uca ,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%deallocate(triplet_per_uca  ,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%deallocate(offset_per_ssa,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%deallocate(ind           ,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)

        if ( verbosity .gt. 0 ) then
            t1=walltime()
            write(*,*) '... stored Z triplets (',tochar(t1-t0),'s)'
            t0=t1
        endif
    end block Ztriplet
    endif

    if ( dj%n_eps_singlet .gt. 0 ) then
    epssing: block
        integer :: is,js

        allocate(dj%eps_singlet(dj%n_eps_singlet))
        do is=1,dj%n_eps_singlet
            js=ss%info%index_in_unitcell(is)
            dj%eps_singlet(is)%m=di%eps_singlet(js)%m
        enddo
    end block epssing
    endif

    if ( dj%n_eps_pair .gt. 0 ) then
    epspair: block
        real(r8), dimension(3) :: v1,v2
        integer, dimension(:,:), allocatable :: pair_per_uca,ind
        integer, dimension(:), allocatable :: npair_per_uca,offset_per_ssa
        integer :: ip,jp,a1,a2,a3

        ! Count pairs per atom
        call mem%allocate(npair_per_uca,uc%na,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        npair_per_uca=0
        do ip=1,di%n_eps_pair
            if ( mod(ip,mw%n) .ne. mw%r ) cycle
            a1=di%eps_pair(ip)%a1
            npair_per_uca(a1)=npair_per_uca(a1)+1
        enddo
        call mw%allreduce('sum',npair_per_uca)

        ! Store pairs per atom
        call mem%allocate(pair_per_uca,[maxval(npair_per_uca),uc%na],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        pair_per_uca=0
        npair_per_uca=0
        do ip=1,di%n_eps_pair
            a1=di%eps_pair(ip)%a1
            npair_per_uca(a1)=npair_per_uca(a1)+1
            pair_per_uca(npair_per_uca(a1),a1)=ip
        enddo

        ! Now start mapping for reals
        call mem%allocate(offset_per_ssa,ss%na,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%allocate(ind,[3,dj%n_eps_pair],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        offset_per_ssa=0
        ind=0

        ! Well not just yet, I need the offset in index per supercell pair to
        ! make it parallel and neat. Also gets the pairs sorted somewhat.
        ip=0
        do a1=1,ss%na
            offset_per_ssa(a1)=ip
            a2=ss%info%index_in_unitcell(a1)
            ip=ip+npair_per_uca(a2)
        enddo

        ! Start locating pairs?
        do a1=1,ss%na
            a2=ss%info%index_in_unitcell(a1)
            v1=ss%r(:,a1)
            do ip=1,npair_per_uca(a2)

                ! make it parallel?
                if ( mod(ip+offset_per_ssa(a1),mw%n) .ne. mw%r ) cycle
                jp=pair_per_uca(ip,a2)  ! index to unitcell pair

                v2=v1+matmul(ss%inv_latticevectors,di%eps_pair(jp)%v)
                v2=lo_clean_fractional_coordinates(v2)
                v2=lo_chop(v2,1E-11_r8)
                v2=lo_clean_fractional_coordinates(v2)
                ! Now I have the location of the second atom in fractional coordinates.
                a3=vb%locate(ss%r,v2)
                if ( a3 .gt. 0 ) then
                    ind(1,ip+offset_per_ssa(a1))=a1
                    ind(2,ip+offset_per_ssa(a1))=a3
                    ind(3,ip+offset_per_ssa(a1))=jp
                else
                    call lo_stop_gracefully(['Could not locate pair in supercell'],lo_exitcode_param,__FILE__,__LINE__)
                endif
            enddo
        enddo
        call mw%allreduce('max',ind)

        if ( verbosity .gt. 0 ) then
            t1=walltime()
            write(*,*) '... sorted out eps pairs (',tochar(t1-t0),'s)'
            t0=t1
        endif

        ! Now store tensors and things properly
        allocate(dj%eps_pair(dj%n_eps_pair))
        do ip=1,dj%n_eps_pair
            a1=ind(1,ip)
            a2=ind(2,ip)
            jp=ind(3,ip)
            dj%eps_pair(ip)%a1=a1
            dj%eps_pair(ip)%a2=a2
            dj%eps_pair(ip)%v=di%eps_pair(jp)%v
            ! v  = a2 + lv - a1
            ! lv = v  + a1 - a2
            v1=dj%eps_pair(ip)%v+ss%rcart(:,a1)-ss%rcart(:,a2)
            v1=matmul(ss%inv_latticevectors,v1)
            if ( sum(abs(v1-anint(v1))) .gt. lo_sqtol ) then
                call lo_stop_gracefully(['I do not understand vectors'],lo_exitcode_param,__FILE__,__LINE__)
            else
                v1=anint(v1)
                dj%eps_pair(ip)%lv=matmul(ss%latticevectors,v1)
            endif
            ! And the actual tensor
            dj%eps_pair(ip)%m=di%eps_pair(jp)%m
        enddo

        call mem%deallocate(npair_per_uca ,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%deallocate(pair_per_uca  ,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%deallocate(offset_per_ssa,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%deallocate(ind           ,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)

        if ( verbosity .gt. 0 ) then
            t1=walltime()
            write(*,*) '... stored eps pairs (',tochar(t1-t0),'s)'
            t0=t1
        endif
    end block epspair
    endif

    ! Destroy the Verlet boxes
    call vb%destroy(mem)
end subroutine

!> write to plaintext file
subroutine writetofile(di,p,fn)
    !> dielectric interactions
    class(lo_dielectric_tensor), intent(in) :: di
    !> structure
    type(lo_crystalstructure), intent(in) :: p
    !> filename
    character(len=*), intent(in) :: fn

    real(r8), dimension(3) :: v0
    integer :: u,i,ii,jj,kk

    u=open_file('out',trim(fn))
        write(u,'(1X,I6,1X,A)') di%n_atom,'! number of atoms'
        write(u,'(1X,I6,1X,A)') di%n_Z_singlet,'! number of born charge singlets'
        write(u,'(1X,I6,1X,A)') di%n_Z_pair,'! number of born charge pairs'
        write(u,'(1X,I6,1X,A)') di%n_Z_triplet,'! number of born charge triplets'
        write(u,'(1X,I6,1X,A)') di%n_eps_singlet,'! number of eps singlets'
        write(u,'(1X,I6,1X,A)') di%n_eps_pair,'! number of eps singlets'


        write(u,"(1X,3(1X,E19.12),1X,A)") di%eps_inf(1,:),'! eps xx xy xz'
        write(u,"(1X,3(1X,E19.12),1X,A)") di%eps_inf(2,:),'! eps yx yy yz'
        write(u,"(1X,3(1X,E19.12),1X,A)") di%eps_inf(3,:),'! eps zx zy zz'

        do i=1,di%n_Z_singlet
            write(u,'(1X,I6,1X,A)') i,'! index to atom'
            do ii=1,3
                v0=lo_chop(di%Z_singlet(ii,:,i)/lo_bohr_to_A,1E-13_r8)
                write(u,*) v0
            enddo
        enddo

        do i=1,di%n_Z_pair
            write(u,'(1X,I6,1X,I6,1X,A)') di%Z_pair(i)%a1,di%Z_pair(i)%a2,'! indices to atoms in pair '//tochar(i)
            v0=matmul(p%inv_latticevectors,di%Z_pair(i)%lv)
            v0=anint(v0)
            write(u,'(3(1X,E21.14),1X,A)') v0,'! lattice vector'
            do ii=1,3
            do jj=1,3
                v0=lo_chop(di%Z_pair(i)%m(ii,jj,:)/lo_bohr_to_A**2,1E-13_r8)
                write(u,*) v0
            enddo
            enddo
        enddo
        do i=1,di%n_Z_triplet
            write(u,'(1X,I6,1X,I6,1X,I6,1X,A)') di%Z_triplet(i)%a1,di%Z_triplet(i)%a2,di%Z_triplet(i)%a3,'! indices to atoms in triplet '//tochar(i)
            v0=matmul(p%inv_latticevectors,di%Z_triplet(i)%lv2)
            v0=anint(v0)
            write(u,'(3(1X,E21.14),1X,A)') v0,'! lattice vector 2'
            v0=matmul(p%inv_latticevectors,di%Z_triplet(i)%lv3)
            v0=anint(v0)
            write(u,'(3(1X,E21.14),1X,A)') v0,'! lattice vector 3'
            do ii=1,3
            do jj=1,3
            do kk=1,3
                v0=lo_chop(di%Z_triplet(i)%m(ii,jj,kk,:)/lo_bohr_to_A**3,1E-13_r8)
                write(u,*) v0
            enddo
            enddo
            enddo
        enddo

        do i=1,di%n_eps_singlet
            do ii=1,3
            do jj=1,3
                v0=lo_chop(di%eps_singlet(i)%m(ii,jj,:)/lo_bohr_to_A,1E-13_r8)
                write(u,*) v0
            enddo
            enddo
        enddo
        do i=1,di%n_eps_pair
            write(u,'(1X,I6,1X,I6,1X,A)') di%eps_pair(i)%a1,di%eps_pair(i)%a2,'! indices to atoms in pair '//tochar(i)
            v0=matmul(p%inv_latticevectors,di%eps_pair(i)%lv)
            v0=anint(v0)
            write(u,'(3(1X,E21.14),1X,A)') v0,'! lattice vector'
            do ii=1,3
            do jj=1,3
            do kk=1,3
                v0=lo_chop(di%eps_pair(i)%m(ii,jj,kk,:)/lo_bohr_to_A**2,1E-13_r8)
                write(u,*) v0
            enddo
            enddo
            enddo
        enddo
    close(u)
end subroutine

!> read from plaintext
subroutine readfromfile(di,p,fn)
    !> dielectric interactions
    class(lo_dielectric_tensor), intent(out) :: di
    !> structure
    type(lo_crystalstructure), intent(in) :: p
    !> filename
    character(len=*), intent(in) :: fn

    real(r8), dimension(3) :: v0
    integer :: u,i,ii,jj,kk

    u=open_file('in',trim(fn))
        read(u,*) di%n_atom
        read(u,*) di%n_Z_singlet
        read(u,*) di%n_Z_pair
        read(u,*) di%n_Z_triplet
        read(u,*) di%n_eps_singlet
        read(u,*) di%n_eps_pair

        do i=1,3
            read(u,*) di%eps_inf(i,:)
        enddo

        if ( di%n_Z_singlet .gt. 0 ) then
            allocate(di%Z_singlet(3,3,di%n_Z_singlet))
            di%Z_singlet=0.0_r8
        endif
        if ( di%n_Z_pair .gt. 0 ) then
            allocate(di%Z_pair(di%n_Z_pair))
        endif
        if ( di%n_Z_triplet .gt. 0 ) then
            allocate(di%Z_triplet(di%n_Z_triplet))
        endif
        if ( di%n_eps_singlet .gt. 0 ) then
            allocate(di%eps_singlet(di%n_eps_singlet))
        endif
        if ( di%n_eps_pair .gt. 0 ) then
            allocate(di%eps_pair(di%n_eps_pair))
        endif

        do i=1,di%n_Z_singlet
            read(u,*) kk
            do ii=1,3
                read(u,*) v0
                di%Z_singlet(ii,:,i)=v0*lo_bohr_to_A
            enddo
        enddo

        do i=1,di%n_Z_pair
            read(u,*) di%Z_pair(i)%a1,di%Z_pair(i)%a2
            read(u,*) v0
            di%Z_pair(i)%lv=v0
            di%Z_pair(i)%v=p%r(:,di%Z_pair(i)%a2)+v0-p%r(:,di%Z_pair(i)%a1)
            di%Z_pair(i)%lv=lo_chop(matmul(p%latticevectors,di%Z_pair(i)%lv),1E-12_r8)
            di%Z_pair(i)%v=lo_chop(matmul(p%latticevectors,di%Z_pair(i)%v),1E-12_r8)
            do ii=1,3
            do jj=1,3
                read(u,*) v0
                di%Z_pair(i)%m(ii,jj,:)=v0*lo_bohr_to_A**2
            enddo
            enddo
        enddo

        do i=1,di%n_Z_triplet
            read(u,*) di%Z_triplet(i)%a1,di%Z_triplet(i)%a2,di%Z_triplet(i)%a3
            read(u,*) v0
            di%Z_triplet(i)%lv2=v0
            di%Z_triplet(i)%v2=p%r(:,di%Z_triplet(i)%a2)+v0-p%r(:,di%Z_triplet(i)%a1)
            di%Z_triplet(i)%lv2=lo_chop(matmul(p%latticevectors,di%Z_triplet(i)%lv2),1E-12_r8)
            di%Z_triplet(i)%v2=lo_chop(matmul(p%latticevectors,di%Z_triplet(i)%v2),1E-12_r8)
            read(u,*) v0
            di%Z_triplet(i)%lv3=v0
            di%Z_triplet(i)%v3=p%r(:,di%Z_triplet(i)%a3)+v0-p%r(:,di%Z_triplet(i)%a1)
            di%Z_triplet(i)%lv3=lo_chop(matmul(p%latticevectors,di%Z_triplet(i)%lv3),1E-12_r8)
            di%Z_triplet(i)%v3=lo_chop(matmul(p%latticevectors,di%Z_triplet(i)%v3),1E-12_r8)
            do ii=1,3
            do jj=1,3
            do kk=1,3
                read(u,*) v0
                di%Z_triplet(i)%m(ii,jj,kk,:)=v0*lo_bohr_to_A**3
            enddo
            enddo
            enddo
        enddo

        do i=1,di%n_eps_singlet
            do ii=1,3
            do jj=1,3
                read(u,*) v0
                di%eps_singlet(i)%m(ii,jj,:)=v0*lo_bohr_to_A
            enddo
            enddo
        enddo

        do i=1,di%n_eps_pair
            read(u,*) di%eps_pair(i)%a1,di%eps_pair(i)%a2
            read(u,*) v0
            di%eps_pair(i)%lv=v0
            di%eps_pair(i)%v =p%r(:,di%eps_pair(i)%a2)+v0-p%r(:,di%eps_pair(i)%a1)
            di%eps_pair(i)%lv=lo_chop(matmul(p%latticevectors,di%eps_pair(i)%lv),1E-12_r8)
            di%eps_pair(i)%v =lo_chop(matmul(p%latticevectors,di%eps_pair(i)%v) ,1E-12_r8)
            do ii=1,3
            do jj=1,3
            do kk=1,3
                read(u,*) v0
                di%eps_pair(i)%m(ii,jj,kk,:)=v0*lo_bohr_to_A**2
            enddo
            enddo
            enddo
        enddo
    close(u)
end subroutine

!> size in memory
function size_in_mem(di) result(mem)
    !> dielectric interactions
    class(lo_dielectric_tensor), intent(in) :: di
    !> memory, in bytes
    integer(i8) :: mem

    mem=0
    mem=mem+storage_size(di)
    if ( allocated(di%Z_pair      ) ) mem=mem+storage_size(di%Z_pair     )*size(di%Z_pair     )
    if ( allocated(di%Z_triplet   ) ) mem=mem+storage_size(di%Z_triplet  )*size(di%Z_triplet  )
    if ( allocated(di%eps_singlet ) ) mem=mem+storage_size(di%eps_singlet)*size(di%eps_singlet)
    if ( allocated(di%eps_pair    ) ) mem=mem+storage_size(di%eps_pair   )*size(di%eps_pair   )
    mem=mem/8
end function

!> destroy
subroutine destroy(di)
    !> dielectric interactions
    class(lo_dielectric_tensor), intent(inout) :: di

    if ( allocated(di%Z_pair      ) ) deallocate(di%Z_pair     )
    if ( allocated(di%Z_triplet   ) ) deallocate(di%Z_triplet  )
    if ( allocated(di%eps_singlet ) ) deallocate(di%eps_singlet)
    if ( allocated(di%eps_pair    ) ) deallocate(di%eps_pair   )
    di%n_atom=-lo_hugeint
    di%n_Z_pair=-lo_hugeint
    di%n_Z_triplet=-lo_hugeint
    di%n_eps_singlet=-lo_hugeint
    di%n_eps_pair=-lo_hugeint
end subroutine

end module
