#include "precompilerdefinitions"
submodule(type_forcemap) type_forcemap_io
implicit none
contains

!> store the forcemap on disk in a hdf5 file
module subroutine write_to_hdf5(map, filename, verbosity)
    !> forcemap
    class(lo_forcemap), intent(in) :: map
    !> filename
    character(len=*), intent(in) :: filename
    !> talk?
    integer, intent(in) :: verbosity

    ! Number of digits in my enumeration of tuplets
    integer, parameter :: ndig = 6
    ! Other stuffs
    integer(HID_T) :: file_id
    real(r8) :: t0

    t0 = walltime()

    ! Dump all the constants and basic attributes
    init: block
        ! open hdf5 interface and a file
        call h5open_f(lo_status)
        call h5fcreate_f(trim(filename), H5F_ACC_TRUNC_F, file_id, lo_status)

        if (verbosity .gt. 0) then
            write (*, *) ''
            write (*, *) 'WRITING FORCEMAP TO FILE'
        end if

        ! store basic attributes
        ! call lo_h5_store_attribute(map%n_atom_uc,                  file_id,'number_of_unitcell_atoms' )
        ! call lo_h5_store_attribute(map%n_atom_ss,                  file_id,'number_of_supercell_atoms')
        ! call lo_h5_store_attribute(map%have_fc_singlet,            file_id,'contains_firstorder_forceconstant' )
        ! call lo_h5_store_attribute(map%have_fc_pair,               file_id,'contains_secondorder_forceconstant')
        ! call lo_h5_store_attribute(map%thirdorder,                 file_id,'contains_thirdorder_forceconstant' )
        ! call lo_h5_store_attribute(map%fourthorder,                file_id,'contains_fourthorder_forceconstant')
        ! call lo_h5_store_attribute(map%magnetic_pair_interactions, file_id,'contains_magnetic_jij')
        ! call lo_h5_store_attribute(map%polar,                      file_id,'contains_polar_information')
        ! call lo_h5_store_attribute(map%polarcorrectiontype,        file_id,'type_of_polar_correction')
        ! !call lo_h5_store_attribute(map%nfc_singlet,                file_id,'number_of_firstorder_ifc' )
        ! !call lo_h5_store_attribute(map%nx_fc_pair,                 file_id,'number_of_secondorder_ifc')
        ! call lo_h5_store_attribute(map%nfc_triplet,                file_id,'number_of_thirdorder_ifc' )
        ! call lo_h5_store_attribute(map%nfc_quartet,                file_id,'number_of_fourthorder_ifc')
        ! ! call lo_h5_store_attribute(map%nx_Z,                   file_id,'number_of_borncharge_theta')
        ! ! call lo_h5_store_attribute(map%nx_eps,                 file_id,'number_of_dielectric_theta')
        ! ! call lo_h5_store_attribute(map%nx_magpair_jij,         file_id,'number_of_magnetic_pair_theta_jij')
        ! ! call lo_h5_store_attribute(map%nx_magpair_tij,         file_id,'number_of_magnetic_pair_theta_tij')
        ! !call lo_h5_store_attribute(map%nsingletshells,             file_id,'n_firstorder_shells' )
        ! call lo_h5_store_attribute(map%n_fc_pair_shell,                file_id,'n_secondorder_shells')
        ! call lo_h5_store_attribute(map%ntripletshells,             file_id,'n_thirdorder_shells' )
        ! call lo_h5_store_attribute(map%nquartetshells,             file_id,'n_fourthorder_shells')
        ! call lo_h5_store_attribute(map%nmagpairshells,             file_id,'n_magpair_shells')
        !
        ! ! call lo_h5_store_attribute(map%neq_rot1,       file_id,'n_firstorder_rotational_constraints')
        ! ! call lo_h5_store_attribute(map%neq_rot2,       file_id,'n_secondorder_rotational_constraints')
        ! ! call lo_h5_store_attribute(map%neq_rot3,       file_id,'n_thirdorder_rotational_constraints')
        ! ! call lo_h5_store_attribute(map%neq_asr3,       file_id,'n_thirdorder_sumrule_constraints')
        ! ! call lo_h5_store_attribute(map%neq_asr4,       file_id,'n_fourthorder_sumrule_constraints')

        if (verbosity .gt. 0) write (*, *) '... wrote initial stats'
    end block init

    ! Store the coordination shells and such
    storeshells: block
        integer(HID_T) :: group_id
        integer :: sh

        ! if ( map%firstorder ) then
        !     do sh=1,map%nsingletshells
        !         call h5gcreate_f(file_id,'singletshell_'//tochar(sh,ndig),group_id,lo_status)
        !         call lo_h5_store_attribute(map%singletshell(sh)%nx,         group_id,'ntheta')
        !         if ( map%singletshell(sh)%nx .gt. 0 ) then
        !             call lo_h5_store_data( map%singletshell(sh)%ind_global,       group_id,'thetaind')
        !             call lo_h5_store_data( map%singletshell(sh)%ind_local,    group_id,'thetalocind')
        !             call lo_h5_store_data( map%singletshell(sh)%coeff,      group_id,'redcoeffM')
        !         endif
        !         call h5gclose_f(group_id,lo_status)
        !     enddo
        ! endif

        if (map%have_fc_pair) then
            do sh = 1, map%n_fc_pair_shell
                call h5gcreate_f(file_id, 'pairshell_'//tochar(sh, ndig), group_id, lo_status)
                call lo_h5_store_attribute(map%fc_pair_shell(sh)%nx, group_id, 'ntheta')
                if (map%fc_pair_shell(sh)%nx .gt. 0) then
                    call lo_h5_store_data(map%fc_pair_shell(sh)%ind_global, group_id, 'thetaind')
                    call lo_h5_store_data(map%fc_pair_shell(sh)%ind_local, group_id, 'thetalocind')
                    call lo_h5_store_data(map%fc_pair_shell(sh)%coeff, group_id, 'redcoeffM')
                end if
                call h5gclose_f(group_id, lo_status)
            end do
        end if

        ! if ( map%thirdorder ) then
        !     do sh=1,map%ntripletshells
        !         call h5gcreate_f(file_id,'tripletshell_'//tochar(sh,ndig),group_id,lo_status)
        !         call lo_h5_store_attribute(map%tripletshell(sh)%nx,         group_id,'ntheta')
        !         if ( map%tripletshell(sh)%nx .gt. 0 ) then
        !             call lo_h5_store_data( map%tripletshell(sh)%ind_global,       group_id,'thetaind')
        !             call lo_h5_store_data( map%tripletshell(sh)%ind_local,    group_id,'thetalocind')
        !             call lo_h5_store_data( map%tripletshell(sh)%coeff,      group_id,'redcoeffM')
        !         endif
        !         call h5gclose_f(group_id,lo_status)
        !     enddo
        ! endif

        ! if ( map%fourthorder ) then
        !     do sh=1,map%nquartetshells
        !         call h5gcreate_f(file_id,'quartetshell_'//tochar(sh,ndig),group_id,lo_status)
        !         call lo_h5_store_attribute(map%quartetshell(sh)%nx,         group_id,'ntheta')
        !         if ( map%quartetshell(sh)%nx .gt. 0 ) then
        !             call lo_h5_store_data( map%quartetshell(sh)%ind_global,       group_id,'thetaind')
        !             call lo_h5_store_data( map%quartetshell(sh)%ind_local,    group_id,'thetalocind')
        !             call lo_h5_store_data( map%quartetshell(sh)%coeff,      group_id,'redcoeffM')
        !         endif
        !         call h5gclose_f(group_id,lo_status)
        !     enddo
        ! endif

        ! if ( map%magnetic_pair_interactions ) then
        !     do sh=1,map%nmagpairshells
        !         call h5gcreate_f(file_id,'magpairshell_'//tochar(sh,ndig),group_id,lo_status)
        !         call lo_h5_store_attribute(map%magpairshell(sh)%nx_jij,         group_id,'ntheta_jij')
        !         if ( map%magpairshell(sh)%nx_jij .gt. 0 ) then
        !             call lo_h5_store_data( map%magpairshell(sh)%ind_global_jij,       group_id,'thetaind_jij')
        !             call lo_h5_store_data( map%magpairshell(sh)%ind_local_jij,    group_id,'thetalocind_jij')
        !             call lo_h5_store_data( map%magpairshell(sh)%coeff_jij,      group_id,'redcoeffM_jij')
        !         endif
        !         call lo_h5_store_attribute(map%magpairshell(sh)%nx_tij,         group_id,'ntheta_tij')
        !         if ( map%magpairshell(sh)%nx_tij .gt. 0 ) then
        !             call lo_h5_store_data( map%magpairshell(sh)%ind_global_tij,       group_id,'thetaind_tij')
        !             call lo_h5_store_data( map%magpairshell(sh)%ind_local_tij,    group_id,'thetalocind_tij')
        !             call lo_h5_store_data( map%magpairshell(sh)%coeff_tij,      group_id,'redcoeffM_tij')
        !         endif
        !         call h5gclose_f(group_id,lo_status)
        !     enddo
        ! endif

        ! if ( map%polar ) then
        !     write(*,*) 'FIXME WRITE POLAR FORCEMAP'
        !     stop
        !     ! do sh=1,map%nZshells
        !     !    call h5gcreate_f(file_id,'Zshell_'//tochar(sh,ndig),group_id,lo_status)
        !     !    call lo_h5_store_attribute(map%Zshell(sh)%nx,         group_id,'ntheta')
        !     !    if ( map%Zshell(sh)%nx .gt. 0 ) then
        !     !        call lo_h5_store_data( map%Zshell(sh)%ind_global,       group_id,'thetaind')
        !     !        call lo_h5_store_data( map%Zshell(sh)%ind_local,    group_id,'thetalocind')
        !     !        call lo_h5_store_data( map%Zshell(sh)%coeff,      group_id,'redcoeffM')
        !     !    endif
        !     !    call h5gclose_f(group_id,lo_status)
        !     ! enddo
        !     !
        !     ! call h5gcreate_f(file_id,'epsshell_'//tochar(sh,ndig),group_id,lo_status)
        !     ! call lo_h5_store_attribute(map%epsshell%nx,         group_id,'ntheta')
        !     ! if ( map%epsshell%nx .gt. 0 ) then
        !     !    call lo_h5_store_data( map%epsshell%coeff,      group_id,'redcoeffM')
        !     ! endif
        !     ! call h5gclose_f(group_id,lo_status)
        ! endif
        if (verbosity .gt. 0) write (*, *) '... wrote shell information'
    end block storeshells

    ! dump the forceconstant constraints
    storeconstraints: block
        integer(HID_T) :: group_id

        call h5gcreate_f(file_id, 'constraints', group_id, lo_status)

        call lo_h5_store_attribute(map%constraints%nconstr_tot, group_id, 'nconstr_tot')
        call lo_h5_store_attribute(map%constraints%neq1, group_id, 'neq1')
        call lo_h5_store_attribute(map%constraints%neq2, group_id, 'neq2')
        call lo_h5_store_attribute(map%constraints%neq3, group_id, 'neq3')
        call lo_h5_store_attribute(map%constraints%neq4, group_id, 'neq4')
        !call lo_h5_store_attribute(map%constraints%neqZ,                group_id,'neqZ')
        if (allocated(map%constraints%eq1)) call lo_h5_store_data(map%constraints%eq1, group_id, 'eq1')
        if (allocated(map%constraints%eq2)) call lo_h5_store_data(map%constraints%eq2, group_id, 'eq2')
        if (allocated(map%constraints%eq3)) call lo_h5_store_data(map%constraints%eq3, group_id, 'eq3')
        if (allocated(map%constraints%eq4)) call lo_h5_store_data(map%constraints%eq4, group_id, 'eq4')
        !if ( allocated(map%constraints%eqZ) ) call lo_h5_store_data( map%constraints%eqZ, group_id, 'eqZ')

        call h5gclose_f(group_id, lo_status)
    end block storeconstraints

    ! ! store the symmetry operations
    ! storeops: block
    !     integer(HID_T) :: group_id
    !     integer, dimension(:,:), allocatable :: di
    !     real(r8), dimension(:,:,:), allocatable :: dr1,dr2,dr3
    !     integer :: i,nsinglet,npair,ntriplet,nquartet
    !
    !     ! if ( allocated(map%singletop) ) then
    !     !     nsinglet=size(map%singletop)
    !     ! else
    !     !     nsinglet=0
    !     ! endif
    !     if ( allocated(map%op_pair) ) then
    !         npair=size(map%op_pair)
    !     else
    !         npair=0
    !     endif
    !     if ( allocated(map%tripletop) ) then
    !         ntriplet=size(map%tripletop)
    !     else
    !         ntriplet=0
    !     endif
    !     if ( allocated(map%quartetop) ) then
    !         nquartet=size(map%quartetop)
    !     else
    !         nquartet=0
    !     endif
    !
    !     call h5gcreate_f(file_id,'operations', group_id, lo_status)
    !     call lo_h5_store_attribute(nsinglet  , group_id, 'n_singlet_operations')
    !     call lo_h5_store_attribute(npair     , group_id, 'n_pair_operations')
    !     call lo_h5_store_attribute(ntriplet  , group_id, 'n_triplet_operations')
    !     call lo_h5_store_attribute(nquartet  , group_id, 'n_quartet_operations')
    !
    !     if ( nsinglet .gt. 0 ) then
    !         allocate(di(2,nsinglet))
    !         allocate(dr1(3,3,nsinglet))
    !         allocate(dr2(3,3,nsinglet))
    !         do i=1,nsinglet
    !             di(1,i)=map%singletop(i)%permind
    !             di(2,i)=map%singletop(i)%opind
    !             dr1(:,:,i)=map%singletop(i)%m3
    !             dr2(:,:,i)=map%singletop(i)%im3
    !         enddo
    !         call lo_h5_store_data(di,  group_id,'singletind')
    !         call lo_h5_store_data(dr1, group_id,'singlet_m3')
    !         call lo_h5_store_data(dr2, group_id,'singlet_im3')
    !         deallocate(di)
    !         deallocate(dr1)
    !         deallocate(dr2)
    !     endif
    !     if ( npair .gt. 0 ) then
    !         allocate(di(4,npair))
    !         allocate(dr1(3,3,npair))
    !         allocate(dr2(3,3,npair))
    !         allocate(dr3(9,9,npair))
    !         do i=1,npair
    !             di(1,i)=map%op_pair(i)%permind
    !             di(2,i)=map%op_pair(i)%opind
    !             di(3:4,i)=map%op_pair(i)%perm
    !             dr1(:,:,i)=map%op_pair(i)%m3
    !             dr2(:,:,i)=map%op_pair(i)%im3
    !             dr3(:,:,i)=map%op_pair(i)%sotr
    !         enddo
    !         call lo_h5_store_data(di,   group_id,'pairind')
    !         call lo_h5_store_data(dr1,  group_id,'pair_m3')
    !         call lo_h5_store_data(dr2,  group_id,'pair_im3')
    !         call lo_h5_store_data(dr3,  group_id,'pair_sotr')
    !         deallocate(dr1)
    !         deallocate(dr2)
    !         deallocate(dr3)
    !         deallocate(di)
    !     endif
    !     if ( ntriplet .gt. 0 ) then
    !         allocate(di(5,ntriplet))
    !         allocate(dr1(3,3,ntriplet))
    !         allocate(dr2(3,3,ntriplet))
    !         allocate(dr3(27,27,ntriplet))
    !         do i=1,ntriplet
    !             di(1,i)=map%tripletop(i)%permind
    !             di(2,i)=map%tripletop(i)%opind
    !             di(3:5,i)=map%tripletop(i)%perm
    !             dr1(:,:,i)=map%tripletop(i)%m3
    !             dr2(:,:,i)=map%tripletop(i)%im3
    !             dr3(:,:,i)=map%tripletop(i)%sotr
    !         enddo
    !         call lo_h5_store_data(di,   group_id,'tripletind')
    !         call lo_h5_store_data(dr1,  group_id,'triplet_m3')
    !         call lo_h5_store_data(dr2,  group_id,'triplet_im3')
    !         call lo_h5_store_data(dr3,  group_id,'triplet_sotr')
    !         deallocate(di)
    !         deallocate(dr1)
    !         deallocate(dr2)
    !         deallocate(dr3)
    !     endif
    !     if ( nquartet .gt. 0 ) then
    !         allocate(di(6,nquartet))
    !         allocate(dr1(3,3,nquartet))
    !         allocate(dr2(3,3,nquartet))
    !         allocate(dr3(81,81,nquartet))
    !         do i=1,nquartet
    !             di(1,i)=map%quartetop(i)%permind
    !             di(2,i)=map%quartetop(i)%opind
    !             di(3:6,i)=map%quartetop(i)%perm
    !             dr1(:,:,i)=map%quartetop(i)%m3
    !             dr2(:,:,i)=map%quartetop(i)%im3
    !             dr3(:,:,i)=map%quartetop(i)%sotr
    !         enddo
    !         call lo_h5_store_data(di,   group_id,'quartetind')
    !         call lo_h5_store_data(dr1,  group_id,'quartet_m3')
    !         call lo_h5_store_data(dr2,  group_id,'quartet_im3')
    !         call lo_h5_store_data(dr3,  group_id,'quartet_sotr')
    !         deallocate(di)
    !         deallocate(dr1)
    !         deallocate(dr2)
    !         deallocate(dr3)
    !     endif
    !     call h5gclose_f(group_id,lo_status)
    !     if ( verbosity .gt. 0 ) write(*,*) '... wrote symmetry operations'
    ! end block storeops

    ! write the unitcell info to file
    storeuc: block
        integer, dimension(:, :), allocatable :: dumi
        integer(HID_T) :: group_id
        integer :: a1, i1

        ! Now dumpt the unitcell information, one group per atom
        do a1 = 1, map%n_atom_uc
            call h5gcreate_f(file_id, 'unitcell_atom_'//tochar(a1, ndig), group_id, lo_status)
            ! store some metadata for this atom
            !call lo_h5_store_attribute(map%uc(a1)%npair,    group_id,'number_of_pairs')
            ! call lo_h5_store_attribute(map%uc(a1)%ntriplet, group_id,'number_of_triplets')
            ! call lo_h5_store_attribute(map%uc(a1)%nquartet, group_id,'number_of_quartets')
            ! if ( map%firstorder ) then
            !     call lo_h5_store_attribute(map%uc(a1)%singlet%irreducible_shell,    group_id,'singlet_irreducible_shell')
            !     call lo_h5_store_attribute(map%uc(a1)%singlet%operation_from_shell, group_id,'singlet_operation')
            ! else
            !     call lo_h5_store_attribute(0,                                       group_id,'singlet_irreducible_shell')
            !     call lo_h5_store_attribute(0,                                       group_id,'singlet_operation')
            ! endif
            !@todo Switch this to Jii
            !if ( map%polar ) then
            !    call lo_h5_store_attribute(map%uc(a1)%Z%irreducible_shell,          group_id,'Z_irreducible_shell')
            !    call lo_h5_store_attribute(map%uc(a1)%Z%operation_from_shell,       group_id,'Z_operation')
            !else
            !    call lo_h5_store_attribute(0,                                       group_id,'Z_irreducible_shell')
            !    call lo_h5_store_attribute(0,                                       group_id,'Z_operation')
            !endif
            ! Acoustic sum rules
            ! call lo_h5_store_attribute(map%uc(a1)%ngroup_triplet, group_id,'ngroup_triplet')
            ! if ( map%uc(a1)%ngroup_triplet .gt. 0 ) then
            !     call lo_h5_store_attribute(map%uc(a1)%group_triplet_selfterm, group_id,'group_triplet_selfterm')
            !     call lo_h5_store_data(map%uc(a1)%group_triplet_ctr,group_id,'group_triplet_ctr')
            !     call lo_h5_store_data(map%uc(a1)%group_triplet_ind,group_id,'group_triplet_ind')
            ! endif
            ! call lo_h5_store_attribute(map%uc(a1)%ngroup_quartet, group_id,'ngroup_quartet')
            ! if ( map%uc(a1)%ngroup_quartet .gt. 0 ) then
            !     call lo_h5_store_attribute(map%uc(a1)%group_quartet_selfterm, group_id,'group_quartet_selfterm')
            !     call lo_h5_store_data(map%uc(a1)%group_quartet_ctr,group_id,'group_quartet_ctr')
            !     call lo_h5_store_data(map%uc(a1)%group_quartet_ind,group_id,'group_quartet_ind')
            ! endif

            ! pack the pairs
            if (map%have_fc_pair) then
                ! ! pack data
                ! allocate(dumi(7,map%uc(a1)%npair))
                ! dumi=0
                ! do i1=1,map%uc(a1)%npair
                !     ! pack the metadata to ints
                !     dumi(1,i1)=map%uc(a1)%pair(i1)%irreducible_shell
                !     dumi(2,i1)=map%uc(a1)%pair(i1)%operation_from_shell
                !     if ( map%uc(a1)%pair(i1)%selfterm ) then
                !         dumi(3,i1)=1
                !     else
                !         dumi(3,i1)=0
                !     endif
                !     dumi(4,i1)=map%uc(a1)%pair(i1)%i2
                !     dumi(5:7,i1)=int(anint(map%uc(a1)%pair(i1)%flv))
                ! enddo
                ! ! store the packed versions
                ! call lo_h5_store_data(dumi,group_id,'pair_int')
                ! deallocate(dumi)
            end if
            ! if ( map%thirdorder ) then
            !     ! pack data
            !     allocate(dumi(11,map%uc(a1)%ntriplet))
            !     dumi=0
            !     do i1=1,map%uc(a1)%ntriplet
            !         ! pack the metadata to ints
            !         dumi(1,i1)=map%uc(a1)%triplet(i1)%irreducible_shell
            !         dumi(2,i1)=map%uc(a1)%triplet(i1)%operation_from_shell
            !         if ( map%uc(a1)%triplet(i1)%selfterm ) then
            !             dumi(3,i1)=1
            !         else
            !             dumi(3,i1)=0
            !         endif
            !         dumi(4,i1)=map%uc(a1)%triplet(i1)%i2
            !         dumi(5,i1)=map%uc(a1)%triplet(i1)%i3
            !         dumi(6:8,i1)=int(anint(map%uc(a1)%triplet(i1)%flv2))
            !         dumi(9:11,i1)=int(anint(map%uc(a1)%triplet(i1)%flv3))
            !     enddo
            !     ! store the packed versions
            !     call lo_h5_store_data(dumi,group_id,'triplet_int')
            !     deallocate(dumi)
            ! endif
            ! if ( map%fourthorder ) then
            !     ! pack data
            !     allocate(dumi(15,map%uc(a1)%nquartet))
            !     dumi=0
            !     do i1=1,map%uc(a1)%nquartet
            !         ! pack the metadata to ints
            !         dumi(1,i1)=map%uc(a1)%quartet(i1)%irreducible_shell
            !         dumi(2,i1)=map%uc(a1)%quartet(i1)%operation_from_shell
            !         if ( map%uc(a1)%quartet(i1)%selfterm ) then
            !             dumi(3,i1)=1
            !         else
            !             dumi(3,i1)=0
            !         endif
            !         dumi(4,i1)=map%uc(a1)%quartet(i1)%i2
            !         dumi(5,i1)=map%uc(a1)%quartet(i1)%i3
            !         dumi(6,i1)=map%uc(a1)%quartet(i1)%i4
            !         dumi(7:9,i1)=int(anint(map%uc(a1)%quartet(i1)%flv2))
            !         dumi(10:12,i1)=int(anint(map%uc(a1)%quartet(i1)%flv3))
            !         dumi(13:15,i1)=int(anint(map%uc(a1)%quartet(i1)%flv4))
            !     enddo
            !     ! store the packed versions
            !     call lo_h5_store_data(dumi,group_id,'quartet_int')
            !     deallocate(dumi)
            ! endif
            ! close the group
            call h5gclose_f(group_id, lo_status)
        end do
        if (verbosity .gt. 0) write (*, *) '... wrote unitcell'

        ! ! Same thing but the supercell
        ! do a1=1,map%nss
        !     call h5gcreate_f(file_id,'supercell_atom_'//tochar(a1,ndig),group_id,lo_status)
        !     ! store some metadata for this atom
        !     call lo_h5_store_attribute(map%ss(a1)%npair,    group_id,'number_of_pairs')
        !     call lo_h5_store_attribute(map%ss(a1)%ntriplet, group_id,'number_of_triplets')
        !     call lo_h5_store_attribute(map%ss(a1)%nquartet, group_id,'number_of_quartets')
        !     call lo_h5_store_attribute(map%ss(a1)%nmagpair, group_id,'number_of_magnetic_pairs')
        !     if ( map%firstorder ) then
        !         call lo_h5_store_attribute(map%ss(a1)%singlet%irreducible_shell,    group_id,'singlet_irreducible_shell')
        !         call lo_h5_store_attribute(map%ss(a1)%singlet%operation_from_shell, group_id,'singlet_operation')
        !     else
        !         call lo_h5_store_attribute(0,                                       group_id,'singlet_irreducible_shell')
        !         call lo_h5_store_attribute(0,                                       group_id,'singlet_operation')
        !     endif
        !
        !     if ( allocated(map%ss(a1)%pairind)    ) call lo_h5_store_data(map%ss(a1)%pairind,group_id,'pairind')
        !     if ( allocated(map%ss(a1)%tripletind) ) call lo_h5_store_data(map%ss(a1)%tripletind,group_id,'tripletind')
        !     if ( allocated(map%ss(a1)%quartetind) ) call lo_h5_store_data(map%ss(a1)%quartetind,group_id,'quartetind')
        !     if ( allocated(map%ss(a1)%magpairind) ) call lo_h5_store_data(map%ss(a1)%magpairind,group_id,'magpairind')
        !     ! close the group
        !     call h5gclose_f(group_id,lo_status)
        ! enddo
        if (verbosity .gt. 0) write (*, *) '... wrote supercell'
    end block storeuc

    ! close file and interface
    call h5fclose_f(file_id, lo_status)
    call h5close_f(lo_status)
    if (verbosity .gt. 0) write (*, *) '... done writing forcemap (', tochar(walltime() - t0), 's)'
end subroutine

!> read the forcemap from hdf5 file
module subroutine read_from_hdf5(map, uc, filename, verbosity)
    !> forcemap
    class(lo_forcemap), intent(out) :: map
    !> unitcell
    type(lo_crystalstructure), intent(in) :: uc
    !> filename
    character(len=*), intent(in) :: filename
    !> talk?
    integer, intent(in) :: verbosity
    !
    integer, parameter :: ndig = 6
    real(r8) :: timer
    integer(HID_T) :: file_id

    ! open file and make some space
    init: block
        if (verbosity .gt. 0) then
            write (*, *) ''
            write (*, *) 'Reading forcemap from file'
        end if
        timer = walltime()
        call h5open_f(lo_status)
        call h5fopen_f(trim(filename), H5F_ACC_RDONLY_F, file_id, lo_status)

        ! Grab constants and sizes of all arrays
        call lo_h5_read_attribute(map%n_atom_uc, file_id, 'number_of_unitcell_atoms')
        call lo_h5_read_attribute(map%n_atom_ss, file_id, 'number_of_supercell_atoms')
        !call lo_h5_read_attribute(map%firstorder,                 file_id,'contains_firstorder_forceconstant' )
        call lo_h5_read_attribute(map%have_fc_pair, file_id, 'contains_secondorder_forceconstant')
        ! call lo_h5_read_attribute(map%thirdorder,                 file_id,'contains_thirdorder_forceconstant' )
        ! call lo_h5_read_attribute(map%fourthorder,                file_id,'contains_fourthorder_forceconstant')
        ! call lo_h5_read_attribute(map%magnetic_pair_interactions, file_id,'contains_magnetic_jij')
        ! call lo_h5_read_attribute(map%polar,                      file_id,'contains_polar_information')
        ! call lo_h5_read_attribute(map%polarcorrectiontype,        file_id,'type_of_polar_correction')
        !call lo_h5_read_attribute(map%nfc_singlet,                file_id,'number_of_firstorder_ifc' )
        !call lo_h5_read_attribute(map%nx_fc_pair,                   file_id,'number_of_secondorder_ifc')
        ! call lo_h5_read_attribute(map%nfc_triplet,                file_id,'number_of_thirdorder_ifc' )
        ! call lo_h5_read_attribute(map%nfc_quartet,                file_id,'number_of_fourthorder_ifc')
        ! call lo_h5_read_attribute(map%nx_Z,                   file_id,'number_of_borncharge_theta')
        ! call lo_h5_read_attribute(map%nx_eps,                 file_id,'number_of_dielectric_theta')
        ! call lo_h5_read_attribute(map%nx_magpair_jij,         file_id,'number_of_magnetic_pair_theta_jij')
        ! call lo_h5_read_attribute(map%nx_magpair_tij,         file_id,'number_of_magnetic_pair_theta_tij')
        !call lo_h5_read_attribute(map%nsingletshells,             file_id,'n_firstorder_shells' )
        ! call lo_h5_read_attribute(map%n_fc_pair_shell,                file_id,'n_secondorder_shells')
        ! call lo_h5_read_attribute(map%ntripletshells,             file_id,'n_thirdorder_shells' )
        ! call lo_h5_read_attribute(map%nquartetshells,             file_id,'n_fourthorder_shells')
        ! call lo_h5_read_attribute(map%nmagpairshells,             file_id,'n_magpair_shells')

        ! ! Start making some space
        ! if ( map%n_atom_uc .gt. 0 ) allocate(map%uc(map%n_atom_uc))
        ! if ( map%n_atom_ss .gt. 0 ) allocate(map%ss(map%n_atom_ss))
        ! ! space for shells
        ! !if ( map%nsingletshells .gt. 0 )     allocate(map%singletshell(map%nsingletshells))
        ! if ( map%n_fc_pair_shell .gt. 0 )        allocate(map%fc_pair_shell(map%n_fc_pair_shell))
        ! if ( map%ntripletshells .gt. 0 )     allocate(map%tripletshell(map%ntripletshells))
        ! if ( map%nquartetshells .gt. 0 )     allocate(map%quartetshell(map%nquartetshells))
        ! if ( map%nmagpairshells .gt. 0 )     allocate(map%magpairshell(map%nmagpairshells))
        !
        ! ! space for forceconstants
        ! !if ( map%nfc_singlet .gt. 0 )        allocate(map%ifc_singlet(map%nfc_singlet))
        ! !if ( map%nx_fc_pair .gt. 0 )           allocate(map%ifc_pair(map%nx_fc_pair))
        ! if ( map%nfc_triplet .gt. 0 )        allocate(map%ifc_triplet(map%nfc_triplet))
        ! if ( map%nfc_quartet .gt. 0 )        allocate(map%ifc_quartet(map%nfc_quartet))
        ! ! if ( map%nx_eps .gt. 0 )         allocate(map%theta_eps(map%nx_eps))
        ! ! if ( map%nx_Z .gt. 0 )           allocate(map%theta_Z(map%nx_Z))
        ! ! if ( map%nx_magpair_jij .gt. 0 ) allocate(map%theta_magpair_jij(map%nx_magpair_jij))
        ! ! if ( map%nx_magpair_tij .gt. 0 ) allocate(map%theta_magpair_tij(map%nx_magpair_tij))
        ! !if ( map%nfc_singlet .gt. 0 )        map%ifc_singlet=0.0_r8
        ! !if ( map%nx_fc_pair .gt. 0 )           map%ifc_pair=0.0_r8
        ! if ( map%nfc_triplet .gt. 0 )        map%ifc_triplet=0.0_r8
        ! if ( map%nfc_quartet .gt. 0 )        map%ifc_quartet=0.0_r8
        ! ! if ( map%nx_eps .gt. 0 )         map%theta_eps=0.0_r8
        ! ! if ( map%nx_Z .gt. 0 )           map%theta_Z=0.0_r8
        ! ! if ( map%nx_magpair_jij .gt. 0 ) map%theta_magpair_jij=0.0_r8
        ! ! if ( map%nx_magpair_tij .gt. 0 ) map%theta_magpair_tij=0.0_r8
    end block init

    ! Store the coordination shells and such
    readshells: block
        integer(HID_T) :: group_id
        integer :: sh

        ! if ( map%firstorder ) then
        !     do sh=1,map%nsingletshells
        !         call h5gopen_f(file_id,'singletshell_'//tochar(sh,ndig),group_id,lo_status)
        !         call lo_h5_read_attribute(map%singletshell(sh)%nx,         group_id,'ntheta')
        !         if ( map%singletshell(sh)%nx .gt. 0 ) then
        !             call lo_h5_read_data( map%singletshell(sh)%ind_global,       group_id,'thetaind')
        !             call lo_h5_read_data( map%singletshell(sh)%ind_local,    group_id,'thetalocind')
        !             call lo_h5_read_data( map%singletshell(sh)%coeff,      group_id,'redcoeffM')
        !         endif
        !         call h5gclose_f(group_id,lo_status)
        !     enddo
        ! endif

        if (map%have_fc_pair) then
            ! do sh=1,map%n_fc_pair_shell
            !     call h5gopen_f(file_id,'pairshell_'//tochar(sh,ndig),group_id,lo_status)
            !     call lo_h5_read_attribute(map%fc_pair_shell(sh)%nx,         group_id,'ntheta')
            !     if ( map%fc_pair_shell(sh)%nx .gt. 0 ) then
            !         call lo_h5_read_data( map%fc_pair_shell(sh)%ind_global, group_id,'thetaind')
            !         call lo_h5_read_data( map%fc_pair_shell(sh)%ind_local,  group_id,'thetalocind')
            !         call lo_h5_read_data( map%fc_pair_shell(sh)%coeff,      group_id,'redcoeffM')
            !     endif
            !     call h5gclose_f(group_id,lo_status)
            ! enddo
        end if

        ! if ( map%thirdorder ) then
        !     do sh=1,map%ntripletshells
        !         call h5gopen_f(file_id,'tripletshell_'//tochar(sh,ndig),group_id,lo_status)
        !         call lo_h5_read_attribute(map%tripletshell(sh)%nx,         group_id,'ntheta')
        !         if ( map%tripletshell(sh)%nx .gt. 0 ) then
        !             call lo_h5_read_data( map%tripletshell(sh)%ind_global,       group_id,'thetaind')
        !             call lo_h5_read_data( map%tripletshell(sh)%ind_local,    group_id,'thetalocind')
        !             call lo_h5_read_data( map%tripletshell(sh)%coeff,      group_id,'redcoeffM')
        !         endif
        !         call h5gclose_f(group_id,lo_status)
        !     enddo
        ! endif
        !
        ! if ( map%fourthorder ) then
        !     do sh=1,map%nquartetshells
        !         call h5gopen_f(file_id,'quartetshell_'//tochar(sh,ndig),group_id,lo_status)
        !         call lo_h5_read_attribute(map%quartetshell(sh)%nx,         group_id,'ntheta')
        !         if ( map%quartetshell(sh)%nx .gt. 0 ) then
        !             call lo_h5_read_data( map%quartetshell(sh)%ind_global,       group_id,'thetaind')
        !             call lo_h5_read_data( map%quartetshell(sh)%ind_local,    group_id,'thetalocind')
        !             call lo_h5_read_data( map%quartetshell(sh)%coeff,      group_id,'redcoeffM')
        !         endif
        !         call h5gclose_f(group_id,lo_status)
        !     enddo
        ! endif

        ! if ( map%magnetic_pair_interactions ) then
        !     do sh=1,map%nmagpairshells
        !         call h5gopen_f(file_id,'magpairshell_'//tochar(sh,ndig),group_id,lo_status)
        !         call lo_h5_read_attribute(map%magpairshell(sh)%nx_jij,         group_id,'ntheta_jij')
        !         if ( map%magpairshell(sh)%nx_jij .gt. 0 ) then
        !             call lo_h5_read_data( map%magpairshell(sh)%ind_global_jij,       group_id,'thetaind_jij')
        !             call lo_h5_read_data( map%magpairshell(sh)%ind_local_jij,    group_id,'thetalocind_jij')
        !             call lo_h5_read_data( map%magpairshell(sh)%coeff_jij,      group_id,'redcoeffM_jij')
        !         endif
        !         call lo_h5_read_attribute(map%magpairshell(sh)%nx_tij,         group_id,'ntheta_tij')
        !         if ( map%magpairshell(sh)%nx_tij .gt. 0 ) then
        !             call lo_h5_read_data( map%magpairshell(sh)%ind_global_tij,       group_id,'thetaind_tij')
        !             call lo_h5_read_data( map%magpairshell(sh)%ind_local_tij,    group_id,'thetalocind_tij')
        !             call lo_h5_read_data( map%magpairshell(sh)%coeff_tij,      group_id,'redcoeffM_tij')
        !         endif
        !         call h5gclose_f(group_id,lo_status)
        !     enddo
        ! endif

        ! if ( map%polar ) then
        !     write(*,*) 'FIXME READ POLAR FORCEMAP'
        !     stop
        !     ! do sh=1,map%nZshells
        !     !    call h5gopen_f(file_id,'Zshell_'//tochar(sh,ndig),group_id,lo_status)
        !     !    call lo_h5_read_attribute(map%Zshell(sh)%nx,         group_id,'ntheta')
        !     !    if ( map%Zshell(sh)%nx .gt. 0 ) then
        !     !        call lo_h5_read_data( map%Zshell(sh)%ind_global,       group_id,'thetaind')
        !     !        call lo_h5_read_data( map%Zshell(sh)%ind_local,    group_id,'thetalocind')
        !     !        call lo_h5_read_data( map%Zshell(sh)%coeff,      group_id,'redcoeffM')
        !     !    endif
        !     !    call h5gclose_f(group_id,lo_status)
        !     ! enddo
        !     !
        !     ! call h5gopen_f(file_id,'epsshell_'//tochar(sh,ndig),group_id,lo_status)
        !     ! call lo_h5_read_attribute(map%epsshell%nx,         group_id,'ntheta')
        !     ! if ( map%epsshell%nx .gt. 0 ) then
        !     !    call lo_h5_read_data( map%epsshell%coeff,      group_id,'redcoeffM')
        !     ! endif
        !     ! call h5gclose_f(group_id,lo_status)
        ! endif
        if (verbosity .gt. 0) write (*, *) '... read shell information'
    end block readshells

    ! ! read the symmetry operations
    ! readops: block
    !     integer(HID_T) :: group_id
    !     integer, dimension(:,:), allocatable :: di
    !     real(r8), dimension(:,:,:), allocatable :: dr1,dr2,dr3
    !     integer :: i,nsinglet,npair,ntriplet,nquartet
    !
    !     call h5gopen_f(file_id,'operations', group_id, lo_status)
    !     call lo_h5_read_attribute(nsinglet , group_id, 'n_singlet_operations')
    !     call lo_h5_read_attribute(npair    , group_id, 'n_pair_operations')
    !     call lo_h5_read_attribute(ntriplet , group_id, 'n_triplet_operations')
    !     call lo_h5_read_attribute(nquartet , group_id, 'n_quartet_operations')
    !
    !     if ( nsinglet .gt. 0 ) then
    !         allocate(map%singletop(nsinglet))
    !         call lo_h5_read_data(di,  group_id,'singletind')
    !         call lo_h5_read_data(dr1, group_id,'singlet_m3')
    !         call lo_h5_read_data(dr2, group_id,'singlet_im3')
    !         do i=1,nsinglet
    !             map%singletop(i)%permind=di(1,i)
    !             map%singletop(i)%opind=di(2,i)
    !             map%singletop(i)%m3=dr1(:,:,i)
    !             map%singletop(i)%im3=dr2(:,:,i)
    !         enddo
    !         deallocate(di)
    !         deallocate(dr1)
    !         deallocate(dr2)
    !     endif
    !     if ( npair .gt. 0 ) then
    !         allocate(map%op_pair(npair))
    !         call lo_h5_read_data(di,   group_id,'pairind')
    !         call lo_h5_read_data(dr1,  group_id,'pair_m3')
    !         call lo_h5_read_data(dr2,  group_id,'pair_im3')
    !         call lo_h5_read_data(dr3,  group_id,'pair_sotr')
    !         do i=1,npair
    !             map%op_pair(i)%permind=di(1,i)
    !             map%op_pair(i)%opind=di(2,i)
    !             map%op_pair(i)%perm=di(3:4,i)
    !             map%op_pair(i)%m3=dr1(:,:,i)
    !             map%op_pair(i)%im3=dr2(:,:,i)
    !             map%op_pair(i)%sotr=dr3(:,:,i)
    !         enddo
    !         deallocate(dr1)
    !         deallocate(dr2)
    !         deallocate(dr3)
    !         deallocate(di)
    !     endif
    !     if ( ntriplet .gt. 0 ) then
    !         allocate(map%tripletop(ntriplet))
    !         call lo_h5_read_data(di,   group_id,'tripletind')
    !         call lo_h5_read_data(dr1,  group_id,'triplet_m3')
    !         call lo_h5_read_data(dr2,  group_id,'triplet_im3')
    !         call lo_h5_read_data(dr3,  group_id,'triplet_sotr')
    !         do i=1,ntriplet
    !             map%tripletop(i)%permind=di(1,i)
    !             map%tripletop(i)%opind=di(2,i)
    !             map%tripletop(i)%perm=di(3:5,i)
    !             map%tripletop(i)%m3=dr1(:,:,i)
    !             map%tripletop(i)%im3=dr2(:,:,i)
    !             map%tripletop(i)%sotr=dr3(:,:,i)
    !         enddo
    !         deallocate(di)
    !         deallocate(dr1)
    !         deallocate(dr2)
    !         deallocate(dr3)
    !     endif
    !     if ( nquartet .gt. 0 ) then
    !         allocate(map%quartetop(nquartet))
    !         call lo_h5_read_data(di,   group_id,'quartetind')
    !         call lo_h5_read_data(dr1,  group_id,'quartet_m3')
    !         call lo_h5_read_data(dr2,  group_id,'quartet_im3')
    !         call lo_h5_read_data(dr3,  group_id,'quartet_sotr')
    !         do i=1,nquartet
    !             map%quartetop(i)%permind=di(1,i)
    !             map%quartetop(i)%opind=di(2,i)
    !             map%quartetop(i)%perm=di(3:6,i)
    !             map%quartetop(i)%m3=dr1(:,:,i)
    !             map%quartetop(i)%im3=dr2(:,:,i)
    !             map%quartetop(i)%sotr=dr3(:,:,i)
    !         enddo
    !         deallocate(di)
    !         deallocate(dr1)
    !         deallocate(dr2)
    !         deallocate(dr3)
    !     endif
    !     call h5gclose_f(group_id,lo_status)
    !     if ( verbosity .gt. 0 ) write(*,*) '... read symmetry operations'
    ! end block readops

    ! dump the forceconstant constraints
    readconstraints: block
        integer(HID_T) :: group_id

        call h5gopen_f(file_id, 'constraints', group_id, lo_status)

        call lo_h5_read_attribute(map%constraints%nconstr_tot, group_id, 'nconstr_tot')
        call lo_h5_read_attribute(map%constraints%neq1, group_id, 'neq1')
        call lo_h5_read_attribute(map%constraints%neq2, group_id, 'neq2')
        call lo_h5_read_attribute(map%constraints%neq3, group_id, 'neq3')
        call lo_h5_read_attribute(map%constraints%neq4, group_id, 'neq4')
        !call lo_h5_read_attribute(map%constraints%neqZ,          group_id,'neqZ')
        if (map%constraints%neq1 .gt. 0) call lo_h5_read_data(map%constraints%eq1, group_id, 'eq1')
        if (map%constraints%neq2 .gt. 0) call lo_h5_read_data(map%constraints%eq2, group_id, 'eq2')
        if (map%constraints%neq3 .gt. 0) call lo_h5_read_data(map%constraints%eq3, group_id, 'eq3')
        if (map%constraints%neq4 .gt. 0) call lo_h5_read_data(map%constraints%eq4, group_id, 'eq4')
        !if ( map%constraints%neqZ .gt. 0 ) call lo_h5_read_data( map%constraints%eqZ, group_id, 'eqZ')

        call h5gclose_f(group_id, lo_status)
    end block readconstraints

    ! read unitcell information
    readuc: block
        real(r8), dimension(3) :: lv
        integer, dimension(:, :), allocatable :: dumi
        integer(HID_T) :: group_id
        integer :: a1, i1

        ! Now dumpt the unitcell information, one group per atom
        ! do a1=1,map%nuc
        !     call h5gopen_f(file_id,'unitcell_atom_'//tochar(a1,ndig),group_id,lo_status)
        !     call lo_h5_read_attribute(map%uc(a1)%npair,    group_id,'number_of_pairs')
        !     call lo_h5_read_attribute(map%uc(a1)%ntriplet, group_id,'number_of_triplets')
        !     call lo_h5_read_attribute(map%uc(a1)%nquartet, group_id,'number_of_quartets')
        !
        !     if ( map%uc(a1)%npair .gt. 0 )    allocate(map%uc(a1)%pair( map%uc(a1)%npair ))
        !     if ( map%uc(a1)%ntriplet .gt. 0 ) allocate(map%uc(a1)%triplet( map%uc(a1)%ntriplet ))
        !     if ( map%uc(a1)%nquartet .gt. 0 ) allocate(map%uc(a1)%quartet( map%uc(a1)%nquartet ))
        !
        !     if ( map%firstorder ) then
        !         call lo_h5_read_attribute(map%uc(a1)%singlet%irreducible_shell,    group_id,'singlet_irreducible_shell')
        !         call lo_h5_read_attribute(map%uc(a1)%singlet%operation_from_shell, group_id,'singlet_operation')
        !     endif
        !     !@todo Switch this to Jii
        !     !if ( map%polar ) then
        !     !    call lo_h5_store_attribute(map%uc(a1)%Z%irreducible_shell,          group_id,'Z_irreducible_shell')
        !     !    call lo_h5_store_attribute(map%uc(a1)%Z%operation_from_shell,       group_id,'Z_operation')
        !     !else
        !     !    call lo_h5_store_attribute(0,                                       group_id,'Z_irreducible_shell')
        !     !    call lo_h5_store_attribute(0,                                       group_id,'Z_operation')
        !     !endif
        !     ! Acoustic sum rules
        !     call lo_h5_read_attribute(map%uc(a1)%ngroup_triplet, group_id,'ngroup_triplet')
        !     if ( map%uc(a1)%ngroup_triplet .gt. 0 ) then
        !         call lo_h5_read_attribute(map%uc(a1)%group_triplet_selfterm, group_id,'group_triplet_selfterm')
        !         call lo_h5_read_data(map%uc(a1)%group_triplet_ctr,group_id,'group_triplet_ctr')
        !         call lo_h5_read_data(map%uc(a1)%group_triplet_ind,group_id,'group_triplet_ind')
        !     endif
        !     call lo_h5_read_attribute(map%uc(a1)%ngroup_quartet, group_id,'ngroup_quartet')
        !     if ( map%uc(a1)%ngroup_quartet .gt. 0 ) then
        !         call lo_h5_read_attribute(map%uc(a1)%group_quartet_selfterm, group_id,'group_quartet_selfterm')
        !         call lo_h5_read_data(map%uc(a1)%group_quartet_ctr,group_id,'group_quartet_ctr')
        !         call lo_h5_read_data(map%uc(a1)%group_quartet_ind,group_id,'group_quartet_ind')
        !     endif
        !
        !     ! pack the pairs
        !     if ( map%have_fc_pair ) then
        !         call lo_h5_read_data(dumi,group_id,'pair_int')
        !         map%uc(a1)%npair=size(dumi,2)
        !         do i1=1,map%uc(a1)%npair
        !             map%uc(a1)%pair(i1)%irreducible_shell=dumi(1,i1)
        !             map%uc(a1)%pair(i1)%operation_from_shell=dumi(2,i1)
        !             if ( dumi(3,i1) .eq. 1 ) then
        !                 map%uc(a1)%pair(i1)%selfterm=.true.
        !             else
        !                 map%uc(a1)%pair(i1)%selfterm=.false.
        !             endif
        !             map%uc(a1)%pair(i1)%i2=dumi(4,i1)
        !             lv=dumi(5:7,i1)
        !             map%uc(a1)%pair(i1)%flv=lv
        !             map%uc(a1)%pair(i1)%lv=uc%fractional_to_cartesian(lv)
        !             map%uc(a1)%pair(i1)%r=lo_chop( uc%fractional_to_cartesian(lv)+uc%rcart(:,map%uc(a1)%pair(i1)%i2)-uc%rcart(:,a1), lo_sqtol )
        !         enddo
        !         deallocate(dumi)
        !     endif
        !     if ( map%thirdorder ) then
        !         call lo_h5_read_data(dumi,group_id,'triplet_int')
        !         map%uc(a1)%ntriplet=size(dumi,2)
        !         do i1=1,map%uc(a1)%ntriplet
        !             map%uc(a1)%triplet(i1)%irreducible_shell=dumi(1,i1)
        !             map%uc(a1)%triplet(i1)%operation_from_shell=dumi(2,i1)
        !             if ( dumi(3,i1) .eq. 1 ) then
        !                 map%uc(a1)%triplet(i1)%selfterm=.true.
        !             else
        !                 map%uc(a1)%triplet(i1)%selfterm=.false.
        !             endif
        !             map%uc(a1)%triplet(i1)%i2=dumi(4,i1)
        !             map%uc(a1)%triplet(i1)%i3=dumi(5,i1)
        !             lv=dumi(6:8,i1)
        !             map%uc(a1)%triplet(i1)%flv2=lv
        !             map%uc(a1)%triplet(i1)%lv2=uc%fractional_to_cartesian(lv)
        !             map%uc(a1)%triplet(i1)%v2=lo_chop( uc%fractional_to_cartesian(lv)+uc%rcart(:,map%uc(a1)%triplet(i1)%i2)-uc%rcart(:,a1), lo_sqtol )
        !             lv=dumi(9:11,i1)
        !             map%uc(a1)%triplet(i1)%flv3=lv
        !             map%uc(a1)%triplet(i1)%lv3=uc%fractional_to_cartesian(lv)
        !             map%uc(a1)%triplet(i1)%v3=lo_chop( uc%fractional_to_cartesian(lv)+uc%rcart(:,map%uc(a1)%triplet(i1)%i3)-uc%rcart(:,a1), lo_sqtol )
        !         enddo
        !         deallocate(dumi)
        !     endif
        !     if ( map%fourthorder ) then
        !         call lo_h5_read_data(dumi,group_id,'quartet_int')
        !         map%uc(a1)%nquartet=size(dumi,2)
        !         do i1=1,map%uc(a1)%nquartet
        !             map%uc(a1)%quartet(i1)%irreducible_shell=dumi(1,i1)
        !             map%uc(a1)%quartet(i1)%operation_from_shell=dumi(2,i1)
        !             if ( dumi(3,i1) .eq. 1 ) then
        !                 map%uc(a1)%quartet(i1)%selfterm=.true.
        !             else
        !                 map%uc(a1)%quartet(i1)%selfterm=.false.
        !             endif
        !             map%uc(a1)%quartet(i1)%i2=dumi(4,i1)
        !             map%uc(a1)%quartet(i1)%i3=dumi(5,i1)
        !             map%uc(a1)%quartet(i1)%i4=dumi(6,i1)
        !             lv=dumi(7:9,i1)
        !             map%uc(a1)%quartet(i1)%flv2=lv
        !             map%uc(a1)%quartet(i1)%lv2=uc%fractional_to_cartesian(lv)
        !             map%uc(a1)%quartet(i1)%v2=lo_chop( uc%fractional_to_cartesian(lv)+uc%rcart(:,map%uc(a1)%quartet(i1)%i2)-uc%rcart(:,a1), lo_sqtol )
        !             lv=dumi(10:12,i1)
        !             map%uc(a1)%quartet(i1)%flv3=lv
        !             map%uc(a1)%quartet(i1)%lv3=uc%fractional_to_cartesian(lv)
        !             map%uc(a1)%quartet(i1)%v3=lo_chop( uc%fractional_to_cartesian(lv)+uc%rcart(:,map%uc(a1)%quartet(i1)%i3)-uc%rcart(:,a1), lo_sqtol )
        !             lv=dumi(13:15,i1)
        !             map%uc(a1)%quartet(i1)%flv4=lv
        !             map%uc(a1)%quartet(i1)%lv4=uc%fractional_to_cartesian(lv)
        !             map%uc(a1)%quartet(i1)%v4=lo_chop( uc%fractional_to_cartesian(lv)+uc%rcart(:,map%uc(a1)%quartet(i1)%i4)-uc%rcart(:,a1), lo_sqtol )
        !         enddo
        !         deallocate(dumi)
        !     endif
        !     ! close the group
        !     call h5gclose_f(group_id,lo_status)
        ! enddo
        if (verbosity .gt. 0) write (*, *) '... read unitcell'

        ! ! Same thing but the supercell
        ! do a1=1,map%nss
        !     call h5gopen_f(file_id,'supercell_atom_'//tochar(a1,ndig),group_id,lo_status)
        !     ! store some metadata for this atom
        !     call lo_h5_read_attribute(map%ss(a1)%npair,    group_id,'number_of_pairs')
        !     call lo_h5_read_attribute(map%ss(a1)%ntriplet, group_id,'number_of_triplets')
        !     call lo_h5_read_attribute(map%ss(a1)%nquartet, group_id,'number_of_quartets')
        !     call lo_h5_read_attribute(map%ss(a1)%nmagpair, group_id,'number_of_magnetic_pairs')
        !     if ( map%firstorder ) then
        !         call lo_h5_read_attribute(map%ss(a1)%singlet%irreducible_shell,    group_id,'singlet_irreducible_shell')
        !         call lo_h5_read_attribute(map%ss(a1)%singlet%operation_from_shell, group_id,'singlet_operation')
        !     endif
        !     if ( map%ss(a1)%npair .gt. 0 )    call lo_h5_read_data(map%ss(a1)%pairind,group_id,'pairind')
        !     if ( map%ss(a1)%ntriplet .gt. 0 ) call lo_h5_read_data(map%ss(a1)%tripletind,group_id,'tripletind')
        !     if ( map%ss(a1)%nquartet .gt. 0 ) call lo_h5_read_data(map%ss(a1)%quartetind,group_id,'quartetind')
        !     if ( map%ss(a1)%nmagpair .gt. 0 ) call lo_h5_read_data(map%ss(a1)%magpairind,group_id,'magpairind')
        !     ! close the group
        !     call h5gclose_f(group_id,lo_status)
        ! enddo
        if (verbosity .gt. 0) write (*, *) '... read supercell'
    end block readuc

    ! close the file
    call h5fclose_f(file_id, lo_status)
    call h5close_f(lo_status)
end subroutine

end submodule
