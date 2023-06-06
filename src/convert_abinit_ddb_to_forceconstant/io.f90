#include "precompilerdefinitions"
module io
! Centralized place where you put routines that read data from files
use konstanter, only: flyt, lo_huge, lo_status, lo_tol, lo_sqtol, lo_twopi, lo_tiny
use gottochblandat, only: open_file, lo_clean_fractional_coordinates, tochar, lo_chop, lo_progressbar_init, lo_progressbar, &
                          lo_nullspace_coefficient_matrix, lo_flattentensor, lo_linear_least_squares, lo_identitymatrix, &
                          lo_unflatten_2tensor
use type_crystalstructure, only: lo_crystalstructure
use type_qpointmesh, only: lo_qpoint
use type_symmetryoperation, only: lo_expandoperation_pair, lo_operate_on_secondorder_tensor
implicit none

private 
public :: read_ddb_file

contains

!> get raw dynamical matrices from ddb files
subroutine read_ddb_file(filelist, uc, dynmat, qvecs, born_effective_charges, dielectric_tensor, verbosity)
    !> list of ddb files
    character(len=*), dimension(:), intent(in) :: filelist
    !> crystal structure
    type(lo_crystalstructure), intent(out) :: uc
    !> dynamical matrices
    complex(flyt), dimension(:, :, :), allocatable, intent(out) :: dynmat
    !> q-vectors
    real(flyt), dimension(:, :), allocatable, intent(out) :: qvecs
    !> born effective charges
    real(flyt), dimension(:, :, :), allocatable, intent(out) :: born_effective_charges
    !> dielectric tensor
    real(flyt), dimension(3, 3), intent(out) :: dielectric_tensor
    !> talk a lot?
    integer, intent(in) :: verbosity
    integer, dimension(:),allocatable :: prtorder
    integer, dimension(:),allocatable :: inverse_prtorder
    integer :: qctr

    ! grab the structure from the first file
    readstructure: block
        real(flyt), dimension(:,:,:), allocatable :: sym_rot
        real(flyt), dimension(:,:), allocatable :: positions,sym_tr,positions_ordered
        real(flyt), dimension(:), allocatable :: d1,znucl
        integer, dimension(:), allocatable :: atomic_numbers,typat,nb_kpt,atomic_numbers_ordered
        real(flyt), dimension(3,3) :: latticevectors
        real(flyt), dimension(3) :: v0
        integer :: u, i, j
        integer :: na, nk, nb, nq, usepaw
        integer :: nsppol, nsym, ntypat, occopt, nspden, nspinor
        character(len=500) :: ds, ds1, ds2, ds3 !,ds4,ds5,ds6,ds7
        character(len=500) :: longline

        u = open_file('in', trim(filelist(1)))
        do i = 1, 6 ! skip stupid lines
            read (u, *)
        end do
        ! Grab the number of atoms, k-points and bands
        read (u, *) ds, usepaw !usepaw
        call ddb_chkname(ds, 'usepaw')
        read (u, *) ds, na
        call ddb_chkname(ds, 'natom')
        read (u, *) ds, nk
        call ddb_chkname(ds, 'nkpt')
        read (u, *) ds, nsppol
        call ddb_chkname(ds, 'nsppol')
        read (u, *) ds, nsym
        call ddb_chkname(ds, 'nsym')
        read (u, *) ds, ntypat
        call ddb_chkname(ds, 'ntypat')
        read (u, *) ds, occopt
        call ddb_chkname(ds, 'occopt')
        !smth special for nband, nband can be different for each kpt (when
        !occpot==2)
        if (occopt == 2) then
            lo_allocate(nb_kpt(nk))
            read (u, *) ds, nb_kpt
            call ddb_chkname(ds, 'nband')
        else
            read (u, *) ds, nb !nband
            call ddb_chkname(ds, 'nband')
        end if
        read (u, *) ds, v0  ! acell
        call ddb_chkname(ds, 'acell')
        lo_allocate(d1(ntypat))
        read (u, *) ds, d1
        call ddb_chkname(ds, 'amu')
        lo_deallocate(d1)
        !if paw is used, pawecutdg is in the DDB file so we have to skip one more line
        if (usepaw == 0) then
            j = 6 + nk + 2
        else
            j = 7 + nk + 2
        end if
        do i = 1, j ! pointless lines
            read (u, *)
        end do
        read (u, *) ds, nspden
        call ddb_chkname(ds, 'nspden')
        read (u, *) ds, nspinor
        call ddb_chkname(ds, 'nspinor')
        lo_allocate(d1(nb))
        read (u, *) ds, d1
        call ddb_chkname(ds, 'occ')
        lo_deallocate(d1)
        read (u, *) ds, latticevectors
        call ddb_chkname(ds, 'rprim')
        lo_allocate(d1(1))
        read (u, *) ds, d1 ! sciss
        ! No check for this one, because the expected name can be either sciss or
        ! dfpt_sciss
        lo_deallocate(d1)
        lo_allocate(d1(na*3))
        read (u, *) ds, d1 ! spinat
        call ddb_chkname(ds, 'spinat')
        lo_deallocate(d1)
        lo_allocate(d1(nsym))
        read (u, *) ds, d1 ! symafm
        call ddb_chkname(ds, 'symafm')
        lo_deallocate(d1)
        lo_allocate(sym_rot(3, 3, nsym)) ! symops rotation
        read (u, *) ds, sym_rot
        call ddb_chkname(ds, 'symrel')
        !lo_deallocate(d3)
        !lo_allocate(d2(3,nsym))   ! symops translation
        lo_allocate(sym_tr(3, nsym))   ! symops translation
        !@todo Here I should  really keep the rotation and translation operations.
        read (u, *) ds, sym_tr
        call ddb_chkname(ds, 'tnons')
        !lo_deallocate(d2)
        do i = 1, 3
            read (u, *)
        end do
        lo_allocate(typat(na))  !typat
        read (u, *) ds, typat
        call ddb_chkname(ds, 'typat')
        lo_allocate(d1(nk)) ! k-weights
        read (u, *) ds, d1
        call ddb_chkname(ds, 'wtk')
        lo_deallocate(d1)
        lo_allocate(positions(3, na))
        read (u, *) ds, positions  !xred
        call ddb_chkname(ds, 'xred')
        lo_allocate(znucl(ntypat)) !znucl
        read (u, *) ds, znucl
        call ddb_chkname(ds, 'znucl')
        close(u)

        ! Now I should have everything to create the structure!
        do i=1,3
            latticevectors(:,i)=latticevectors(:,i)*v0(i)
        enddo
        lo_allocate(atomic_numbers_ordered(na))
        lo_allocate(atomic_numbers(na))
        do i=1,na
            atomic_numbers(i)=int(anint(znucl(typat(i))))
        enddo
        !In abinit, the atoms may not be in vasp order, we have to check that
        !and to rearrange the order of the atoms and related quantities.

        call sort_atomic_numbers(atomic_numbers,atomic_numbers_ordered,prtorder,inverse_prtorder)

        lo_allocate(positions_ordered(3,na))
        do i=1,na
                positions_ordered(:,inverse_prtorder(i)) = positions(:,i)
        enddo
        
        call uc%generate(latticevectors,positions_ordered,atomic_numbers_ordered,enhet=2,header='Structure from Abinit DDB file') !,verbosity=3)
        call uc%classify('wedge',timereversal=.true.)
        !@todo Insert test to see if my symmetry and Abinits symmetry agrees
        if (verbosity .gt. 0) then
            write (*, *) '... grabbed structure from ddb file:'
            call uc%writetofile('stdout', 1)
            write (*, *) ''
        end if
    end block readstructure

    ! Now I have the structure, that is mildly useful. Now start over and read everything again.
    readdynmat: block
        real(flyt), parameter :: bignum = 1E50_flyt
        ! somewhat useful
        real(flyt), dimension(:, :, :, :, :), allocatable :: dmr, dmi
        real(flyt), dimension(:, :, :), allocatable :: dmza, dmzb
        real(flyt), dimension(:, :), allocatable :: qvbuf
        ! dummy things
        real(flyt), dimension(:, :, :), allocatable :: d3
        real(flyt), dimension(:, :), allocatable :: d2, positions
        real(flyt), dimension(:), allocatable :: d1, znucl, zion
        integer, dimension(:), allocatable :: atomic_numbers, typat, block_type, nb_kpt
        real(flyt), dimension(3, 3) :: latticevectors, m0, diel, m1
        real(flyt), dimension(3) :: v0
        real(flyt) :: f0, f1
        integer :: u, i, j, fi, ia, ib, ix, iy, ii, jj
        integer :: na, nk, nb, nq, nfiles, nqtot, nelem, usepaw
        integer :: nsppol, nsym, ntypat, occopt, nspden, nspinor
        integer :: ctr_eps, ctr_Z, ctr_dynmat
        character(len=500) :: ds, ds1, ds2, ds3, longline
        character(len=32) :: name

        ! Create empty buffers
        nqtot = 1000 ! or some other random number
        lo_allocate(dmr(3, 3, uc%na, uc%na, nqtot))
        lo_allocate(dmi(3, 3, uc%na, uc%na, nqtot))
        lo_allocate(dmza(3, 3, uc%na))
        lo_allocate(dmzb(3, 3, uc%na))
        lo_allocate(qvbuf(3, nqtot))
        lo_allocate(zion(uc%na))
        dmr = bignum
        dmi = bignum
        dmza = bignum
        dmzb = bignum
        qvbuf = bignum
        zion = bignum
        diel = bignum
        qctr = 0
        ctr_eps = 0
        ctr_Z = 0
        ctr_dynmat = 0
        ! Read all files
        if (verbosity .gt. 0) call lo_progressbar_init()
        nfiles = size(filelist, 1)
        do fi = 1, nfiles
            u = open_file('in', trim(filelist(fi)))
            do i = 1, 6 ! skip stupid lines
                read (u, *)
            end do
            read (u, *) ds, usepaw
            read (u, *) ds, na
            read (u, *) ds, nk
            read (u, *) ds, nsppol
            read (u, *) ds, nsym
            read (u, *) ds, ntypat
            read (u, *) ds, occopt
            !smth special for nband, nband can be different for each kpt (when
            !occpot==2)
            if (occopt == 2) then
                lo_allocate(nb_kpt(nk))
                read (u, *) ds, nb_kpt
            else
                read (u, *) ds, nb !nband
            end if
            read (u, *) ds, v0  ! acell
            lo_allocate(d1(ntypat))
            read (u, *) ds, d1
            lo_deallocate(d1)
            if (usepaw == 0) then !if paw is used, pawecutdg is in the DDB file
                !so we have to skip one more line
                j = 6 + nk + 2
            else
                j = 7 + nk + 2
            end if
            do i = 1, j ! pointless lines
                read (u, *)
            end do
            read (u, *) ds, nspden
            read (u, *) ds, nspinor
            lo_allocate(d1(nb)) ! occ
            read (u, *) ds, d1
            lo_deallocate(d1)
            read (u, *) ds, latticevectors
            lo_allocate(d1(1))
            read (u, *) ds, d1 ! sciss
            lo_deallocate(d1)
            lo_allocate(d1(na*3))
            read (u, *) ds, d1 ! spinat
            lo_deallocate(d1)
            lo_allocate(d1(nsym))
            read (u, *) ds, d1 ! symafm
            lo_deallocate(d1)
            lo_allocate(d3(3, 3, nsym)) ! symops rotation
            read (u, *) ds, d3
            lo_deallocate(d3)
            lo_allocate(d2(3, nsym))   ! symops translation
            read (u, *) ds, d2
            lo_deallocate(d2)
            do i = 1, 3  ! nothing important
                read (u, *)
            end do
            lo_allocate(typat(na))
            read (u, *) ds, typat
            lo_allocate(d1(nk)) ! k-weights
            read (u, *) ds, d1
            lo_deallocate(d1) !xred
            lo_allocate(positions(3, na))
            read (u, *) ds, positions
            lo_allocate(znucl(ntypat)) !znucl
            read (u, *) ds, znucl
            lo_allocate(d1(ntypat)) !zion
            read (u, *) ds, d1
            lo_allocate(atomic_numbers(na))
            do i=1,na
                atomic_numbers(i)=int(anint(znucl(typat(i))))
                zion(i)=anint(d1(typat(i)))
            enddo
            lo_deallocate(d1) 
            do i=1,3
                latticevectors(:,i)=latticevectors(:,i)*v0(i)
            enddo    ! Sanity test the structure perhaps? That has to be the same
            ! Sanity test the structure perhaps? That has to be the same
            i=0
            if ( sum(abs(latticevectors-uc%latticevectors)) .gt. lo_tol ) i=i+1
            !if ( sum(abs(atomic_numbers-uc%atomic_number)) .gt. lo_tol ) i=i+1
            if ( i .gt. 0 ) then
                write(*,*) 'ERROR: mismatch in structures'
                stop
            end if

            lo_deallocate(atomic_numbers)
            lo_deallocate(typat)
            lo_deallocate(positions)
            lo_deallocate(znucl)

            ! read infos over pseudopotentials not important and the size changes
            do
                read (u, '(A)') longline
                ds = adjustl(longline)
                if (ds(1:4) .eq. '****') then
                    read (longline, *) ds, ds1
                    if (trim(ds1) .eq. 'Database') exit
                end if
            end do
            ! how many block is there?
            read (u, *) ds, ds1, ds2, ds3, nq
            lo_allocate(block_type(nq))
            ! Now the important part starts, with dynamical matrices. Maybe.
            ! first check if I have to increment my read buffers
            i = qctr + nq
            if (i .gt. size(qvbuf, 2)) then
                do
                    call extend_matrices(dmr, dmi, qvbuf, qctr)
                    if (size(qvbuf, 2) .gt. i) exit
                end do
            end if

            ! Read stuff for the q-points
            write(*,*) "inverse_prtorder" , inverse_prtorder
            do i=1,nq
                read(u,*)
                read(u, '(a32,12x,i8)' )name,nelem
                if(name==' 2nd derivatives (non-stat.)  - ' .or. name==' 2rd derivatives(non-stat.)  - ')then
                    block_type(i)=1
                else if(name==' 2nd derivatives (stationary) - ' .or. name==' 2rd derivatives(stationary) - ')then
                    block_type(i)=2
                else if(name==' 3rd derivatives              - ')then
                    block_type(i)=3
                else if(name==' Total energy                 - ')then
                    block_type(i)=0
                else if(name==' 1st derivatives              - ')then
                    block_type(i)=4
                else if(name==' 2nd eigenvalue derivatives   - ' .or. name==' 2rd eigenvalue derivatives   - ')then
                    block_type(i)=5
                else
                    write (*, '(a,a,a,a)')&
                    &   'The following string appears in the DDB in place of',&
                    &   ' the block type description :', trim(name),&
                    &   'Action: check your DDB.'
                end if

                if (block_type(i) == 1 .or. block_type(i) == 2) then
                    qctr = qctr + 1 !only store the 2nd derivative blocks
                    read (u, *) ds, qvbuf(:, qctr)
                    qvbuf(:, qctr) = lo_clean_fractional_coordinates(qvbuf(:, qctr))
                    do j = 1, nelem
                        read (u, *) ix, ia, iy, ib, f0, f1
                        ! This should mask the dynamical matrices
                        if(ia .le. uc%na .and. ib .le. uc%na ) then
                            ctr_dynmat=ctr_dynmat+1
                            dmr(ix,iy,inverse_prtorder(ia),inverse_prtorder(ib),qctr)=f0
                            dmi(ix,iy,inverse_prtorder(ia),inverse_prtorder(ib),qctr)=f1
                        end if
                        ! This should get the Born effective charges
                        if (ia .le. uc%na .and. ib .eq. uc%na + 2) then
                            dmza(ix, iy, ia) = f0
                            ctr_Z = ctr_Z + 1
                        end if
                        if (ib .le. uc%na .and. ia .eq. uc%na + 2) then
                        ! NB: this is thrown out!!
                            dmzb(ix, iy, ib) = f0
                            ctr_Z = ctr_Z + 1
                        end if
                        ! And this should be the dielectric tensor
                        if (ib .eq. uc%na + 2 .and. ia .eq. uc%na + 2) then
                            diel(ix, iy) = f0
                            ctr_eps = ctr_eps + 1
                        end if
                    end do
                elseif (block_type(i) == 3) then
                    do j = 1, 3 + nelem ! skip the perturbation wavevectors + all the elements
                        read (u, *)
                    end do
                elseif (block_type(i) == 0) then
                    read (u, *)  !skip total energy
                elseif (block_type(i) == 4) then
                    do j = 1, nelem ! skip all the elements
                        read (u, *)
                    end do
                elseif (block_type(i) == 5) then
                    write (*, '(a)') "a block containing the 2nd derivative eigenvalues has been detected. Could you remove it from the DDB file please."
                end if
            end do

            ! Print a summary of what I found
            !write(*,*) '=== info ==='
            !write(*,*) '      na',na
            !write(*,*) '      nk',nk
            !write(*,*) '  nsppol',nsppol
            !write(*,*) '    nsym',nsym
            !write(*,*) '  ntypat',ntypat
            !write(*,*) '  occopt',occopt
            !write(*,*) '      nb',nb
            !write(*,*) '   acell',v0
            !write(*,*) '  nspden',nspden
            !write(*,*) ' nspinor',nspinor
            !write(*,*) '      nq',nspinor
            !write(*,*) '==='

            close (u)
            if (verbosity .gt. 0) call lo_progressbar(' ... reading files', fi, nfiles)
        end do

        if (verbosity .gt. 0) write (*, *) '... read ', tochar(qctr), ' q-points'

        ! Now store this, properly. Start by multiplying in the masses, and convert to Cartesian
        ! coordinates. Also add the energy unit.
        do i = 1, qctr
            do ia = 1, uc%na
            do ib = 1, uc%na
                m0 = dmr(:, :, ia, ib, i)
                m0 = matmul(uc%reciprocal_latticevectors, matmul(m0, transpose(uc%reciprocal_latticevectors)))
                dmr(:, :, ia, ib, i) = m0

                m0 = dmi(:, :, ia, ib, i)
                m0 = matmul(uc%reciprocal_latticevectors, matmul(m0, transpose(uc%reciprocal_latticevectors)))
                dmi(:, :, ia, ib, i) = m0
            end do
            end do
        end do

        ! Store the dynamical matrices and the q-vectors
        lo_allocate(dynmat(3*uc%na, 3*uc%na, qctr))
        lo_allocate(qvecs(3, qctr))
        dynmat = 0.0_flyt
        qvecs = 0.0_flyt
        do i = 1, qctr
            qvecs(:, i) = qvbuf(:, i)
            do ix = 1, 3
            do ia = 1, uc%na
            do iy = 1, 3
            do ib = 1, uc%na
                ii = (ia - 1)*3 + ix
                jj = (ib - 1)*3 + iy
                dynmat(jj, ii, i) = cmplx(dmr(ix, iy, ia, ib, i), dmi(ix, iy, ia, ib, i), flyt)
            end do
            end do
            end do
            end do
        end do

        ! Fill in the blanks, if necessary
        if (ctr_Z .gt. 0) call fill_in_blanks(dmza, diel, uc)

        ! It might have failed and the provided Born charges are stupid
        if (count(abs(dmza) > 1E10_flyt) .gt. 0) ctr_Z = 0
        ! Same with dielectric tensor
        if (count(abs(diel) > 1E10_flyt) .gt. 0) ctr_eps = 0

        ! Store the Born effective charges, whatever units they might be in.
        m1 = 0.0_flyt
        do i = 1, 3
            m1(i, i) = 1.0_flyt
        end do
        lo_allocate(born_effective_charges(3, 3, uc%na))
        if (ctr_Z .gt. 0) then
            do i = 1, uc%na
                m0 = dmza(:, :, i)
                ! From Antoine. Skeptical
                !m0=matmul(transpose(uc%reciprocal_latticevectors),matmul(m0,uc%latticevectors))/lo_twopi
                ! My Guess:
                m0=matmul(uc%reciprocal_latticevectors,matmul(m0,transpose(uc%latticevectors)))/lo_twopi
                m0=m0+m1*zion(i)
! enforce diagonal BEC
                m0(1,2) = 0.0_flyt; m0(2,1) = 0.0_flyt
                m0(1,3) = 0.0_flyt; m0(3,1) = 0.0_flyt
                m0(3,2) = 0.0_flyt; m0(2,3) = 0.0_flyt
                born_effective_charges(:,:,inverse_prtorder(i))=m0
            enddo
        else
            born_effective_charges = 0.0_flyt
        end if
        if (ctr_eps .gt. 0) then
            m0 = diel
            m0 = matmul(uc%latticevectors, matmul(m0, transpose(uc%latticevectors)))
            m0 = m0/(lo_twopi**2)
            m0 = m1 - 2*lo_twopi*m0/uc%volume
            dielectric_tensor = lo_chop(m0, lo_sqtol)
        else
            dielectric_tensor = 0.0_flyt
        end if
    end block readdynmat
end subroutine

! Recover the Born effective charges and dielectric tensor in case not all are there
subroutine fill_in_blanks(Z, epsilon, uc)
    !> Raw Born effective charges
    real(flyt), dimension(:, :, :), intent(inout) :: Z
    !> Raw dielectric tensor
    real(flyt), dimension(3, 3), intent(inout) :: epsilon
    !> unitcell
    type(lo_crystalstructure), intent(in) :: uc

    real(flyt), dimension(:, :, :), allocatable :: dz
    real(flyt), dimension(3, 3) :: m0, sa, sb
    integer :: i, j, a1, a2, o

    ! On entry, there might be elements in Z and epsilon set to 1E50, that is they
    ! were not read from file. This routine should fix that. Not that in this routine
    ! the Born charges and dielectric tensor are in bizarro coordinates/units directly
    ! from Abinit, they are just the raw derivatives.

    lo_allocate(dz(3, 3, uc%na))
    do o = 1, uc%sym%n
        dz = 0.0_flyt
        ! Now transform this guys to the best of our abilities
        do a1 = 1, uc%na
            ! Transform this Born charge
            a2 = uc%sym%op(o)%fmap(a1)
            sa = uc%sym%op(o)%rfm
            sb = transpose(uc%sym%op(o)%fm)
            dz(:, :, a2) = matmul(matmul(sa, Z(:, :, a1)), sb)
        end do
        ! See if I can use the transformed to fill in some blanks in the original
        do a1 = 1, uc%na
            do i = 1, 3
            do j = 1, 3
                if (abs(Z(i, j, a1)) .gt. 1E10_flyt) then
                if (abs(dz(i, j, a1)) .lt. 1E10_flyt .and. abs(dz(i, j, a1)) .gt. lo_tiny) then
                    Z(i, j, a1) = dz(i, j, a1)
                end if
                end if
            end do
            end do
        end do
        ! Maybe we are done?
        if (count(abs(Z) > 1E10_flyt) .eq. 0) then
            write (*, *) 'Filled in blanks of Born charges'
            exit
        end if
    end do
    lo_deallocate(dz)

    do o = 1, uc%sym%n
        ! Transform the dielectric tensor
        sa = transpose(uc%sym%op(o)%rfm)
        sb = uc%sym%op(o)%rfm
        m0 = matmul(matmul(sa, epsilon), sb)
        ! See if I can use the transformed to fill in some blanks in the original
        do j = 1, 3
        do i = 1, 3
            if (abs(epsilon(i, j)) .gt. 1E10_flyt) then
            if (abs(m0(i, j)) .lt. 1E10_flyt .and. abs(m0(i, j)) .gt. lo_tiny) then
                epsilon(i, j) = m0(i, j)
            end if
            end if
        end do
        end do
        ! Maybe we are done?
        if (count(abs(epsilon) > 1E10_flyt) .eq. 0) exit
    end do
end subroutine

! ! Recover the Born effective charges in case not all are there
! subroutine quickmassage_borncharge(zba,diel,uc)
!     !> Raw Born effective charges
!     real(flyt), dimension(:,:,:), intent(inout) :: zba
!     !> Raw dielectric tensor
!     real(flyt), dimension(3,3), intent(inout) :: diel
!     !> unitcell
!     type(lo_crystalstructure), intent(in) :: uc
!
!     real(flyt), dimension(:,:,:), allocatable :: dz0,dz1
!     real(flyt), dimension(3,3,8) :: possops
!     real(flyt), dimension(3,3) :: m0,m1,m2,sa,sb
!
!     real(flyt), dimension(3) :: v0
!     integer :: i,j,k,l,a1,a2,o
!
!     lo_allocate(dz0(3,3,uc%na))
!     lo_allocate(dz1(3,3,uc%na))
!     dz0=0.0_flyt
!     dz1=0.0_flyt
!
!
!     do o=5,5
!         possops=0
!         ! Get all the possible operation thingies
!         possops(:,:,1)=uc%sym%op(o)%fm
!         possops(:,:,2)=uc%sym%op(o)%rfm
!         possops(:,:,3)=transpose(uc%sym%op(o)%fm)
!         possops(:,:,4)=transpose(uc%sym%op(o)%rfm)
!         m0=diel
!         do i=1,4
!         do j=1,4
!             sa=possops(:,:,i)
!             sb=possops(:,:,j)
!             m1=matmul(matmul(sa,m0),sb)
!             if ( count(abs(m1)>1E10) .ne. count(abs(m0)>1E10) ) cycle
!             write(*,*) 'op',o,i,j
!             do k=1,3
!                 write(*,"(1X,6(1X,F12.6))") m0(:,k),m1(:,k)
!             enddo
!             write(*,*) count(abs(m1)>1E10)
!         enddo
!         enddo
!     enddo
!
!
! stop
!
!     a1=1
!     do o=5,5
!         possops=0
!         ! Get all the possible operation thingies
!         possops(:,:,1)=uc%sym%op(o)%fm
!         possops(:,:,2)=uc%sym%op(o)%rfm
!         possops(:,:,3)=transpose(uc%sym%op(o)%fm)
!         possops(:,:,4)=transpose(uc%sym%op(o)%rfm)
!
!         ! possops(:,:,3)=uc%sym%op(o)%rfm
!         ! possops(:,:,4)=uc%sym%op(o)%irfm
!         ! do i=1,4
!         !     possops(:,:,4+i)=transpose(possops(:,:,i))
!         ! enddo
!         m0=zba(:,:,a1)
!         do i=1,4
!         do j=1,4
!             sa=possops(:,:,i)
!             sb=possops(:,:,j)
!             m1=matmul(matmul(sa,m0),sb)
!             if ( count(abs(m1)>1E10) .ne. count(abs(m0)>1E10) ) cycle
!             write(*,*) 'op',o,i,j
!             do k=1,3
!                 write(*,"(1X,6(1X,F12.6))") m0(:,k),m1(:,k)
!             enddo
!             write(*,*) count(abs(m1)>1E10)
!         enddo
!         enddo
!         ! sa=uc%sym%op(o)%fm
!         ! sb=uc%sym%op(o)%fm
!         !
!         ! m1=matmul(matmul(sa,m0),sb)
!         ! write(*,*) 'op',o
!         ! do i=1,3
!         !     write(*,"(1X,6(1X,F12.6))") m0(:,i),m1(:,i)
!         ! enddo
!     enddo
!
! end subroutine

!> make matrices larger, if necessary
subroutine extend_matrices(dynbr, dynbi, qvb, ctr)
    real(flyt), dimension(:, :, :, :, :), allocatable, intent(inout) :: dynbr, dynbi
    real(flyt), dimension(:, :), allocatable, intent(inout) :: qvb
    integer, intent(in) :: ctr

    integer, parameter :: increment = 100 ! how much larger to make the matrices
    real(flyt), dimension(:, :, :, :, :), allocatable :: dbr, dbi
    real(flyt), dimension(:, :), allocatable :: qv
    integer :: n, nn

    allocate (dbr, source=dynbr)
    allocate (dbi, source=dynbi)
    allocate (qv, source=qvb)
    ! Store the current
    !dbr=dynbr
    !dbi=dynbi
    !qv=qvb
    ! make larger arrays
    n = size(qvb, 2)
    nn = n + increment
    deallocate (dynbr, dynbi, qvb)
    allocate (dynbr(size(dbr, 1), size(dbr, 2), size(dbr, 3), size(dbr, 4), nn))
    allocate (dynbi(size(dbr, 1), size(dbr, 2), size(dbr, 3), size(dbr, 4), nn))
    allocate (qvb(3, nn))
    dynbr = lo_huge
    dynbi = lo_huge
    qvb = lo_huge
    ! restore the data
    if (ctr .gt. 0) then
        dynbr(:, :, :, :, 1:ctr) = dbr(:, :, :, :, 1:ctr)
        dynbi(:, :, :, :, 1:ctr) = dbi(:, :, :, :, 1:ctr)
        qvb(:, 1:ctr) = qv(:, 1:ctr)
    end if
    deallocate (dbr, dbi, qv)
end subroutine

!> check name of a string, routine probably stolen from Abinit.
subroutine ddb_chkname(nmfound, nmxpct)
    character(len=*), intent(in) :: nmfound, nmxpct

    if (trim(adjustl(nmxpct)) .ne. trim(adjustl(nmfound))) then
        write (*, '(a,a,a,a,a,a,a)')&
        &   'Reading DDB, expected name was "', trim(nmxpct), '" ',&
        &   'and name found is "', trim(nmfound), '".  ',&
        &   'Likely your DDB is incorrect.'
        stop
    end if
end subroutine

subroutine sort_atomic_numbers(atomic_numbers, atomic_numbers_ordered, prtorder, inverse_prtorder)
    integer, Dimension(:), intent(in) :: atomic_numbers
    integer, Dimension(:), allocatable, intent(out):: atomic_numbers_ordered, prtorder, inverse_prtorder
        

    integer :: i,j,temp,na
    !In abinit, the atoms may not be in vasp order, we have to check that
    !and to rearrange the order of the atoms and related quantities.

    ! bubble sort of the atomic_numbers array 
    allocate(atomic_numbers_ordered,source=atomic_numbers)
    allocate(prtorder,source=atomic_numbers)
    atomic_numbers_ordered = atomic_numbers
    na = size(atomic_numbers)
    do i=1,na
        prtorder(i)=i
    enddo
    do i=na-1,1,-1
        do j=1,i
            if (atomic_numbers_ordered(j) > atomic_numbers_ordered(j+1)) then
                temp = atomic_numbers_ordered(j+1) 
                atomic_numbers_ordered(j+1) = atomic_numbers_ordered(j)
                atomic_numbers_ordered(j) = temp
                                
                temp = prtorder(j+1)
                prtorder(j+1) = prtorder(j)
                prtorder(j) = temp
            endif
        end do
    end do
    allocate(inverse_prtorder,source=atomic_numbers)
    do i=1,na
        inverse_prtorder(prtorder(i))=i
    enddo
end subroutine
end module
