#include "precompilerdefinitions"
module hdf5_wrappers
!!
!! Wrapper to read/write hdf5 files. Adding stuff to the interfaces as they are needed. Mostly because I can never remember the syntax.
!!
use konstanter, only: flyt,lo_status,lo_exitcode_io,lo_hugeint
use gottochblandat, only: lo_stop_gracefully,tochar,lo_does_file_exist
use hdf5
use iso_c_binding

implicit none
private
public :: lo_hdf5_helper
public :: lo_h5_store_attribute
public :: lo_h5_store_data
public :: lo_h5_read_attribute
public :: lo_h5_read_data
public :: lo_h5_does_dataset_exist
public :: HID_T,H5F_ACC_TRUNC_F,H5F_ACC_RDONLY_F,H5F_SCOPE_GLOBAL_F
public :: h5close_f,h5open_f,h5fclose_f,h5fopen_f,h5fcreate_f,h5fflush_f
public :: h5gclose_f,h5gopen_f,h5gcreate_f
public :: h5dclose_f,h5dopen_f,h5dcreate_f

!> Container to store all the named constants to pass around, so that I don't need one billion use statements
type lo_hdf5_helper
    !> Store the error codes somewhere
    integer :: errcode=-lo_hugeint
    !> place to store file id
    integer(HID_T) :: file_id=-huge(HID_T)
    !> place to store dataset id
    integer(HID_T) :: dset_id=-huge(HID_T)
    !> place to store group id
    integer(HID_T) :: group_id=-huge(HID_T)
    !> place to store subgroup id
    integer(HID_T) :: subgroup_id=-huge(HID_T)
    contains
        ! open-close-create-destroy
        procedure :: init
        procedure :: destroy
        procedure :: open_file
        procedure :: close_file
        procedure :: open_group
        procedure :: close_group
        procedure :: open_subgroup
        procedure :: close_subgroup
        ! read data
        generic :: read_data=>read_int_1D_array_as_data,&
                              read_int_2D_array_as_data,&
                              read_int_3D_array_as_data,&
                              read_double_1D_array_as_data,&
                              read_double_2D_array_as_data,&
                              read_double_3D_array_as_data,&
                              read_double_4D_array_as_data,&
                              read_double_5D_array_as_data
        procedure, nopass, private :: read_int_1D_array_as_data
        procedure, nopass, private :: read_int_2D_array_as_data
        procedure, nopass, private :: read_int_3D_array_as_data
        procedure, nopass, private :: read_double_1D_array_as_data
        procedure, nopass, private :: read_double_2D_array_as_data
        procedure, nopass, private :: read_double_3D_array_as_data
        procedure, nopass, private :: read_double_4D_array_as_data
        procedure, nopass, private :: read_double_5D_array_as_data

        ! store data
        generic :: store_data=>store_int_1D_array_as_data,&
                               store_int_2D_array_as_data,&
                               store_int_3D_array_as_data,&
                               store_double_1D_array_as_data,&
                               store_double_2D_array_as_data,&
                               store_double_3D_array_as_data,&
                               store_double_4D_array_as_data,&
                               store_double_5D_array_as_data
        procedure, nopass, private :: store_int_1D_array_as_data
        procedure, nopass, private :: store_int_2D_array_as_data
        procedure, nopass, private :: store_int_3D_array_as_data
        procedure, nopass, private :: store_double_1D_array_as_data
        procedure, nopass, private :: store_double_2D_array_as_data
        procedure, nopass, private :: store_double_3D_array_as_data
        procedure, nopass, private :: store_double_4D_array_as_data
        procedure, nopass, private :: store_double_5D_array_as_data

        ! read attributes
        generic :: read_attribute=>read_int_as_attribute,&
                                   read_logical_as_attribute,&
                                   read_float_as_attribute
        procedure, nopass, private :: read_int_as_attribute
        procedure, nopass, private :: read_logical_as_attribute
        procedure, nopass, private :: read_float_as_attribute

        ! store attributes
        generic :: store_attribute=>store_int_as_attribute,&
                                    store_float_as_attribute,&
                                    store_float_1D_as_attribute,&
                                    store_float_2D_as_attribute,&
                                    store_logical_as_attribute,&
                                    store_char_as_attribute
        procedure, nopass, private :: store_int_as_attribute
        procedure, nopass, private :: store_float_as_attribute
        procedure, nopass, private :: store_float_1D_as_attribute
        procedure, nopass, private :: store_float_2D_as_attribute
        procedure, nopass, private :: store_logical_as_attribute
        procedure, nopass, private :: store_char_as_attribute
end type

!> interface to write simple things as attributes
interface lo_h5_store_attribute
    module procedure store_int_as_attribute
    module procedure store_float_as_attribute
    module procedure store_float_1D_as_attribute
    module procedure store_float_2D_as_attribute
    module procedure store_logical_as_attribute
    module procedure store_char_as_attribute
end interface

!> interface to store normal arrays
interface lo_h5_store_data
    module procedure store_int_1D_array_as_data
    module procedure store_int_2D_array_as_data
    module procedure store_int_3D_array_as_data
    module procedure store_double_1D_array_as_data
    module procedure store_double_2D_array_as_data
    module procedure store_double_3D_array_as_data
    module procedure store_double_4D_array_as_data
    module procedure store_double_5D_array_as_data
end interface

!> interface to read normal arrays
interface lo_h5_read_data
    module procedure read_int_1D_array_as_data
    module procedure read_int_2D_array_as_data
    module procedure read_int_3D_array_as_data
    module procedure read_double_1D_array_as_data
    module procedure read_double_2D_array_as_data
    module procedure read_double_3D_array_as_data
    module procedure read_double_4D_array_as_data
    module procedure read_double_5D_array_as_data
end interface

!> interface to read attributes
interface lo_h5_read_attribute
    module procedure read_int_as_attribute
    module procedure read_logical_as_attribute
    module procedure read_float_as_attribute
end interface

contains

!> Initialize HDF5
subroutine init(h5,filename,line)
    !> hdf5 helper
    class(lo_hdf5_helper), intent(out) :: h5
    !> perhaps say where we called from, for help with debugging
    character(len=*), intent(in), optional :: filename
    integer, intent(in), optional :: line

    h5%errcode=0
    ! Start hdf5
    call h5open_f(h5%errcode)
    ! Did it go ok?
    if ( h5%errcode .ne. 0 ) then
        if ( present(filename) .and. present(line) ) then
            call lo_stop_gracefully(['hdf5 nonzero exitcode: '//tochar(h5%errcode)],lo_exitcode_io,trim(filename),line)
        else
            call lo_stop_gracefully(['hdf5 nonzero exitcode: '//tochar(h5%errcode)],lo_exitcode_io)
        endif
    endif

    ! Silence weird errors that confuse me.
    call h5eset_auto_f(0,h5%errcode)
    if ( h5%errcode .ne. 0 ) then
        if ( present(filename) .and. present(line) ) then
            call lo_stop_gracefully(['hdf5 nonzero exitcode: '//tochar(h5%errcode)],lo_exitcode_io,trim(filename),line)
        else
            call lo_stop_gracefully(['hdf5 nonzero exitcode: '//tochar(h5%errcode)],lo_exitcode_io)
        endif
    endif
end subroutine

!> Close hdf5
subroutine destroy(h5,filename,line)
    !> hdf5 helper
    class(lo_hdf5_helper), intent(inout) :: h5
    !> perhaps say where we called from, for help with debugging
    character(len=*), intent(in), optional :: filename
    integer, intent(in), optional :: line

    !@todo insert check that hdf5 is running
    call h5close_f(h5%errcode)
    if ( h5%errcode .ne. 0 ) then
        if ( present(filename) .and. present(line) ) then
            call lo_stop_gracefully(['hdf5 nonzero exitcode: '//tochar(h5%errcode)],lo_exitcode_io,trim(filename),line)
        else
            call lo_stop_gracefully(['hdf5 nonzero exitcode: '//tochar(h5%errcode)],lo_exitcode_io)
        endif
    endif
end subroutine

!> Open a file and return a handle. Will overwrite any existing file.
subroutine open_file(h5,acc,filename)
    !> hdf5 helper
    class(lo_hdf5_helper), intent(inout) :: h5
    !> read or write
    character(len=*), intent(in) :: acc
    !> name of file
    character(len=*), intent(in) :: filename

    integer(HID_T) :: fapl

    select case(trim(adjustl(acc)))
    case('write')
        ! Create the property (?) that makes hdf close the file in a much angrier
        ! way, I hope. Stupid clusters. Maybe I should learn hdf5 better.
        call h5pcreate_f(H5P_FILE_ACCESS_F, fapl, h5%errcode)
        if ( h5%errcode .ne. 0 ) then
            call lo_stop_gracefully(['hdf5 nonzero exitcode: '//tochar(h5%errcode)],lo_exitcode_io)
        endif
        call h5pset_fclose_degree_f(fapl, H5F_CLOSE_STRONG_F, h5%errcode)
        if ( h5%errcode .ne. 0 ) then
            call lo_stop_gracefully(['hdf5 nonzero exitcode: '//tochar(h5%errcode)],lo_exitcode_io)
        endif
        ! Actually create the file. Destroys any file already there.
        call h5fcreate_f(trim(filename), H5F_ACC_TRUNC_F, h5%file_id, h5%errcode,access_prp=fapl)
        if ( h5%errcode .ne. 0 ) then
            call lo_stop_gracefully(['hdf5 nonzero exitcode: '//tochar(h5%errcode)],lo_exitcode_io)
        endif
    case('read')
        if ( lo_does_file_exist(trim(filename)) ) then
            call h5fopen_f(trim(filename), H5F_ACC_RDONLY_F, h5%file_id, h5%errcode)
            if ( h5%errcode .ne. 0 ) then
                call lo_stop_gracefully(['hdf5 nonzero exitcode: '//tochar(h5%errcode)],lo_exitcode_io)
            endif
        else
            call lo_stop_gracefully(['could not find "'//trim(filename)//'"'],lo_exitcode_io)
        endif
    case default
        call lo_stop_gracefully(['Choose either "read" or "write" for hdf5'],lo_exitcode_io)
    end select
end subroutine

!> Open a group
subroutine open_group(h5,acc,groupname)
    !> hdf5 helper
    class(lo_hdf5_helper), intent(inout) :: h5
    !> read or write
    character(len=*), intent(in) :: acc
    !> name of file
    character(len=*), intent(in) :: groupname

    select case(trim(adjustl(acc)))
    case('write')
        ! In this case we have to create the group
        call h5gcreate_f(h5%file_id,trim(groupname),h5%group_id,h5%errcode)
        if ( h5%errcode .ne. 0 ) then
            call lo_stop_gracefully(['hdf5 nonzero exitcode: '//tochar(h5%errcode)],lo_exitcode_io)
        endif
    case('read')
        call h5gopen_f(h5%file_id,trim(groupname),h5%group_id,h5%errcode)
        if ( h5%errcode .ne. 0 ) then
            call lo_stop_gracefully(['hdf5 nonzero exitcode: '//tochar(h5%errcode)],lo_exitcode_io)
        endif
    case default
        call lo_stop_gracefully(['Choose either "read" or "write"'],lo_exitcode_io)
    end select
end subroutine

!> Open a subgroup
subroutine open_subgroup(h5,acc,groupname)
    !> hdf5 helper
    class(lo_hdf5_helper), intent(inout) :: h5
    !> read or write
    character(len=*), intent(in) :: acc
    !> name of file
    character(len=*), intent(in) :: groupname

    select case(trim(adjustl(acc)))
    case('write')
        ! In this case we have to create the group
        call h5gcreate_f(h5%group_id,trim(groupname),h5%subgroup_id,h5%errcode)
        if ( h5%errcode .ne. 0 ) then
            call lo_stop_gracefully(['hdf5 nonzero exitcode: '//tochar(h5%errcode)],lo_exitcode_io)
        endif
    case('read')
        call h5gopen_f(h5%group_id,trim(groupname),h5%subgroup_id,h5%errcode)
        if ( h5%errcode .ne. 0 ) then
            call lo_stop_gracefully(['hdf5 nonzero exitcode: '//tochar(h5%errcode)],lo_exitcode_io)
        endif
    case default
        call lo_stop_gracefully(['Choose either "read" or "write"'],lo_exitcode_io)
    end select
end subroutine

!> Close a file
subroutine close_file(h5)
    !> hdf5 helper
    class(lo_hdf5_helper), intent(inout) :: h5

    ! Flush the file first?

    call h5fclose_f(h5%file_id,h5%errcode)
    if ( h5%errcode .ne. 0 ) then
        call lo_stop_gracefully(['hdf5 nonzero exitcode: '//tochar(h5%errcode)],lo_exitcode_io)
    endif
    h5%file_id=-huge(HID_T)
end subroutine

!> Close a group
subroutine close_group(h5)
    !> hdf5 helper
    class(lo_hdf5_helper), intent(inout) :: h5

    call h5gclose_f(h5%group_id,h5%errcode)
    if ( h5%errcode .ne. 0 ) then
        call lo_stop_gracefully(['hdf5 nonzero exitcode: '//tochar(h5%errcode)],lo_exitcode_io)
    endif
    h5%group_id=-huge(HID_T)
end subroutine

!> Close a group
subroutine close_subgroup(h5)
    !> hdf5 helper
    class(lo_hdf5_helper), intent(inout) :: h5

    call h5gclose_f(h5%subgroup_id,h5%errcode)
    if ( h5%errcode .ne. 0 ) then
        call lo_stop_gracefully(['hdf5 nonzero exitcode: '//tochar(h5%errcode)],lo_exitcode_io)
    endif
    h5%subgroup_id=-huge(HID_T)
end subroutine

!> query if a dataset exist
function lo_h5_does_dataset_exist(obj_id,data_name) result(isthere)
    !> where to look
    integer(HID_T), intent(in) :: obj_id
    !> name thing to look for
    character(len=*), intent(in) :: data_name
    !> is it there?
    logical :: isthere

    integer :: err
    ! just wrap the subroutine to a function
    call h5lexists_f(obj_id,trim(data_name),isthere, err)
end function

subroutine read_int_1D_array_as_data(buf,obj_id,data_name,error)
    !> array to store
    integer, dimension(:), allocatable, intent(out) :: buf
    !> identifier to attach the attribute to
    integer(HID_T), intent(in) :: obj_id
    !> name of the attribute
    character(len=*), intent(in) :: data_name
    !> status
    integer, optional, intent(out) :: error
    !
    integer :: err,rank
    integer(HID_T) :: dspace_id,dset_id
    integer(HSIZE_T), dimension(1) :: dims,maxdims
    ! sife of the array
    call h5dopen_f(obj_id, data_name, dset_id, err)
    ! figure out size
    call h5dget_space_f(dset_id, dspace_id, err)
    call h5sget_simple_extent_ndims_f(dspace_id, rank, err)
    if ( rank .ne. 1 ) then
        write(*,*) 'I expected this data to be 1D'
        stop
    endif
    call h5sget_simple_extent_dims_f(dspace_id, dims, maxdims, err)
    lo_allocate(buf(maxdims(1)))
    call h5dread_f(dset_id, H5T_NATIVE_INTEGER, buf, dims, err)
    call h5dclose_f(dset_id,err)
    if ( present(error) ) error=err
end subroutine
subroutine read_int_2D_array_as_data(buf,obj_id,data_name,error)
    !> array to store
    integer, dimension(:,:), allocatable, intent(out) :: buf
    !> identifier to attach the attribute to
    integer(HID_T), intent(in) :: obj_id
    !> name of the attribute
    character(len=*), intent(in) :: data_name
    !> status
    integer, optional, intent(out) :: error
    !
    integer :: err,rank
    integer(HID_T) :: dspace_id,dset_id
    integer(HSIZE_T), dimension(2) :: dims,maxdims
    ! sife of the array
    call h5dopen_f(obj_id, data_name, dset_id, err)
    ! figure out size
    call h5dget_space_f(dset_id, dspace_id, err)
    call h5sget_simple_extent_ndims_f(dspace_id, rank, err)
    if ( rank .ne. 2 ) then
        write(*,*) 'I expected this data to be 2D'
        stop
    endif
    call h5sget_simple_extent_dims_f(dspace_id, dims, maxdims, err)
    lo_allocate(buf(maxdims(1),maxdims(2)))
    call h5dread_f(dset_id, H5T_NATIVE_INTEGER, buf, dims, err)
    call h5dclose_f(dset_id,err)
    if ( present(error) ) error=err
end subroutine
subroutine read_int_3D_array_as_data(buf,obj_id,data_name,error)
    !> array to store
    integer, dimension(:,:,:), allocatable, intent(out) :: buf
    !> identifier to attach the attribute to
    integer(HID_T), intent(in) :: obj_id
    !> name of the attribute
    character(len=*), intent(in) :: data_name
    !> status
    integer, optional, intent(out) :: error
    !
    integer :: err,rank
    integer(HID_T) :: dspace_id,dset_id
    integer(HSIZE_T), dimension(3) :: dims,maxdims
    ! sife of the array
    call h5dopen_f(obj_id, data_name, dset_id, err)
    ! figure out size
    call h5dget_space_f(dset_id, dspace_id, err)
    call h5sget_simple_extent_ndims_f(dspace_id, rank, err)
    if ( rank .ne. 3 ) then
        write(*,*) 'I expected this data to be 3D'
        stop
    endif
    call h5sget_simple_extent_dims_f(dspace_id, dims, maxdims, err)
    lo_allocate(buf(maxdims(1),maxdims(2),maxdims(3)))
    call h5dread_f(dset_id, H5T_NATIVE_INTEGER, buf, dims, err)
    call h5dclose_f(dset_id,err)
    if ( present(error) ) error=err
end subroutine
subroutine read_double_1D_array_as_data(buf,obj_id,data_name,error)
    !> double to store
    real(flyt), dimension(:), allocatable, intent(out) :: buf
    !> identifier to attach the attribute to
    integer(HID_T), intent(in) :: obj_id
    !> name of the attribute
    character(len=*), intent(in) :: data_name
    !> status
    integer, optional, intent(out) :: error
    !
    integer :: err,rank
    integer(HID_T) :: dspace_id,dset_id
    integer(HSIZE_T), dimension(1) :: dims,maxdims
    ! sife of the array
    call h5dopen_f(obj_id, data_name, dset_id, err)
    ! figure out size
    call h5dget_space_f(dset_id, dspace_id, err)
    call h5sget_simple_extent_ndims_f(dspace_id, rank, err)
    if ( rank .ne. 1 ) then
        write(*,*) 'I expected this data to be 1D'
        stop
    endif
    call h5sget_simple_extent_dims_f(dspace_id, dims, maxdims, err)
    lo_allocate(buf(maxdims(1)))
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, buf, dims, err)
    call h5dclose_f(dset_id,err)
    if ( present(error) ) error=err
end subroutine
subroutine read_double_2D_array_as_data(buf,obj_id,data_name,error)
    !> double to store
    real(flyt), dimension(:,:), allocatable, intent(out) :: buf
    !> identifier to attach the attribute to
    integer(HID_T), intent(in) :: obj_id
    !> name of the attribute
    character(len=*), intent(in) :: data_name
    !> status
    integer, optional, intent(out) :: error
    !
    integer :: err,rank
    integer(HID_T) :: dspace_id,dset_id
    integer(HSIZE_T), dimension(2) :: dims,maxdims
    ! sife of the array
    call h5dopen_f(obj_id, data_name, dset_id, err)
    ! figure out size
    call h5dget_space_f(dset_id, dspace_id, err)
    call h5sget_simple_extent_ndims_f(dspace_id, rank, err)
    if ( rank .ne. 2 ) then
        write(*,*) 'I expected this data to be 2D'
        stop
    endif
    call h5sget_simple_extent_dims_f(dspace_id, dims, maxdims, err)
    lo_allocate(buf(maxdims(1),maxdims(2)))
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, buf, dims, err)
    call h5dclose_f(dset_id,err)
    if ( present(error) ) error=err
end subroutine
subroutine read_double_3D_array_as_data(buf,obj_id,data_name,error)
    !> double to store
    real(flyt), dimension(:,:,:), allocatable, intent(out) :: buf
    !> identifier to attach the attribute to
    integer(HID_T), intent(in) :: obj_id
    !> name of the attribute
    character(len=*), intent(in) :: data_name
    !> status
    integer, optional, intent(out) :: error
    !
    integer :: err,rank
    integer(HID_T) :: dspace_id,dset_id
    integer(HSIZE_T), dimension(3) :: dims,maxdims
    ! sife of the array
    call h5dopen_f(obj_id, data_name, dset_id, err)
    ! figure out size
    call h5dget_space_f(dset_id, dspace_id, err)
    call h5sget_simple_extent_ndims_f(dspace_id, rank, err)
    if ( rank .ne. 3 ) then
        write(*,*) 'I expected this data to be 3D'
        stop
    endif
    call h5sget_simple_extent_dims_f(dspace_id, dims, maxdims, err)
    lo_allocate(buf(maxdims(1),maxdims(2),maxdims(3)))
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, buf, dims, err)
    call h5dclose_f(dset_id,err)
    if ( present(error) ) error=err
end subroutine
subroutine read_double_4D_array_as_data(buf,obj_id,data_name,error)
    !> double to store
    real(flyt), dimension(:,:,:,:), allocatable, intent(out) :: buf
    !> identifier to attach the attribute to
    integer(HID_T), intent(in) :: obj_id
    !> name of the attribute
    character(len=*), intent(in) :: data_name
    !> status
    integer, optional, intent(out) :: error
    !
    integer :: err,rank
    integer(HID_T) :: dspace_id,dset_id
    integer(HSIZE_T), dimension(4) :: dims,maxdims
    ! sife of the array
    call h5dopen_f(obj_id, data_name, dset_id, err)
    ! figure out size
    call h5dget_space_f(dset_id, dspace_id, err)
    call h5sget_simple_extent_ndims_f(dspace_id, rank, err)
    if ( rank .ne. 4 ) then
        write(*,*) 'I expected this data to be 4D'
        stop
    endif
    call h5sget_simple_extent_dims_f(dspace_id, dims, maxdims, err)
    lo_allocate(buf(maxdims(1),maxdims(2),maxdims(3),maxdims(4)))
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, buf, dims, err)
    call h5dclose_f(dset_id,err)
    if ( present(error) ) error=err
end subroutine
subroutine read_double_5D_array_as_data(buf,obj_id,data_name,error)
    !> double to store
    real(flyt), dimension(:,:,:,:,:), allocatable, intent(out) :: buf
    !> identifier to attach the attribute to
    integer(HID_T), intent(in) :: obj_id
    !> name of the attribute
    character(len=*), intent(in) :: data_name
    !> status
    integer, optional, intent(out) :: error
    !
    integer :: err,rank
    integer(HID_T) :: dspace_id,dset_id
    integer(HSIZE_T), dimension(5) :: dims,maxdims
    ! sife of the array
    call h5dopen_f(obj_id, data_name, dset_id, err)
    ! figure out size
    call h5dget_space_f(dset_id, dspace_id, err)
    call h5sget_simple_extent_ndims_f(dspace_id, rank, err)
    if ( rank .ne. 5 ) then
        write(*,*) 'I expected this data to be 5D'
        stop
    endif
    call h5sget_simple_extent_dims_f(dspace_id, dims, maxdims, err)
    lo_allocate(buf(maxdims(1),maxdims(2),maxdims(3),maxdims(4),maxdims(5)))
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, buf, dims, err)
    call h5dclose_f(dset_id,err)
    if ( present(error) ) error=err
end subroutine

subroutine store_int_1D_array_as_data(buf,obj_id,data_name,error,enhet,dimensions)
        !> 1D array of int to store
        integer, dimension(:), intent(in) :: buf
        !> identifier to attach the attribute to
        integer(HID_T), intent(in) :: obj_id
        !> name of the attribute
        character(len=*), intent(in) :: data_name
        !> status
        integer, optional, intent(out) :: error
        !> set the unit right away
        character(len=*), optional, intent(in) :: enhet
        !> specify what the dimensions are
        character(len=*), optional, intent(in) :: dimensions
        !
        integer :: err,rank
        integer(HID_T) :: dspace_id,dset_id
        integer(HSIZE_T), dimension(1) :: dims
        ! sife of the array
        rank=1
        dims(1)=size(buf,1)
        call h5screate_simple_f(rank,dims,dspace_id,err)
            call h5dcreate_f(obj_id,trim(data_name),H5T_NATIVE_INTEGER,dspace_id,dset_id,err)
                call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, buf, dims, err)
                if ( present(enhet) ) call lo_h5_store_attribute(trim(enhet),dset_id,'unit',lo_status)
                if ( present(dimensions) ) call lo_h5_store_attribute(trim(dimensions),dset_id,'dimensions',lo_status)
            call h5dclose_f(dset_id,err)
        call h5sclose_f(dspace_id,err)
        if ( present(error) ) error=err
end subroutine
subroutine store_int_2D_array_as_data(buf,obj_id,data_name,error,enhet,dimensions)
        !> double to store
        integer, dimension(:,:), intent(in) :: buf
        !> identifier to attach the attribute to
        integer(HID_T), intent(in) :: obj_id
        !> name of the attribute
        character(len=*), intent(in) :: data_name
        !> status
        integer, optional, intent(out) :: error
        !> set the unit right away
        character(len=*), optional, intent(in) :: enhet
        !> specify what the dimensions are
        character(len=*), optional, intent(in) :: dimensions
        !
        integer :: err,rank
        integer(HID_T) :: dspace_id,dset_id,plist_id
        integer(HSIZE_T), dimension(2) :: dims,chunks
        ! sife of the array
        rank=2
        dims(1)=size(buf,1)
        dims(2)=size(buf,2)
        chunks=dims ! supposed to be bad, but ok for what I use it for, so far
        call h5screate_simple_f(rank,dims,dspace_id,err)      ! create dataset
            call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, err) ! create properties of dataset
            call h5pset_chunk_f(plist_id, 2, chunks, err)         ! chunk it, whatever that means.
            call h5pset_deflate_f(plist_id, 6, err)               ! szip compression lvl 6, whatever that means
                call h5dcreate_f(obj_id,trim(data_name),H5T_NATIVE_INTEGER,dspace_id,dset_id,err,dcpl_id=plist_id)
                    call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, buf, dims, err) ! write to it
                    if ( present(enhet) ) call lo_h5_store_attribute(trim(enhet),dset_id,'unit',lo_status)
                    if ( present(dimensions) ) call lo_h5_store_attribute(trim(dimensions),dset_id,'dimensions',lo_status)
                call h5dclose_f(dset_id,err)
            call h5pclose_f(plist_id,err)
        call h5sclose_f(dspace_id,err)
        if ( present(error) ) error=err
end subroutine
subroutine store_int_3D_array_as_data(buf,obj_id,data_name,error,enhet,dimensions)
        !> array to store
        integer, dimension(:,:,:), intent(in) :: buf
        !> identifier to attach the attribute to
        integer(HID_T), intent(in) :: obj_id
        !> name of the attribute
        character(len=*), intent(in) :: data_name
        !> status
        integer, optional, intent(out) :: error
        !> set the unit right away
        character(len=*), optional, intent(in) :: enhet
        !> specify what the dimensions are
        character(len=*), optional, intent(in) :: dimensions
        !
        integer :: err,rank
        integer(HID_T) :: dspace_id,dset_id,plist_id
        integer(HSIZE_T), dimension(3) :: dims,chunks
        ! sife of the array
        rank=3
        dims(1)=size(buf,1)
        dims(2)=size(buf,2)
        dims(3)=size(buf,3)
        chunks=dims ! supposed to be bad, but ok for what I use it for, so far
        call h5screate_simple_f(rank,dims,dspace_id,err)      ! create dataset
            call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, err) ! create properties of dataset
            call h5pset_chunk_f(plist_id, 3, chunks, err)         ! chunk it, whatever that means.
            call h5pset_deflate_f(plist_id, 6, err)               ! szip compression lvl 6, whatever that means
                call h5dcreate_f(obj_id,trim(data_name),H5T_NATIVE_INTEGER,dspace_id,dset_id,err,dcpl_id=plist_id)
                    call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, buf, dims, err) ! write to it
                    if ( present(enhet) ) call lo_h5_store_attribute(trim(enhet),dset_id,'unit',lo_status)
                    if ( present(dimensions) ) call lo_h5_store_attribute(trim(dimensions),dset_id,'dimensions',lo_status)
                call h5dclose_f(dset_id,err)
            call h5pclose_f(plist_id,err)
        call h5sclose_f(dspace_id,err)
        if ( present(error) ) error=err
end subroutine

subroutine store_double_1D_array_as_data(buf,obj_id,data_name,error,enhet,dimensions)
        !> double to store
        real(flyt), dimension(:), intent(in) :: buf
        !> identifier to attach the attribute to
        integer(HID_T), intent(in) :: obj_id
        !> name of the attribute
        character(len=*), intent(in) :: data_name
        !> status
        integer, optional, intent(out) :: error
        !> set the unit right away
        character(len=*), optional, intent(in) :: enhet
        !> specify what the dimensions are
        character(len=*), optional, intent(in) :: dimensions
        !
        integer :: err,rank
        integer(HID_T) :: dspace_id,dset_id,plist_id
        integer(HSIZE_T), dimension(1) :: dims,chunks
        ! sife of the array
        rank=1
        dims(1)=size(buf,1)
        chunks=dims
        call h5screate_simple_f(rank,dims,dspace_id,err)
            call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, err) ! create properties of dataset
            call h5pset_chunk_f(plist_id, 1, chunks, err)         ! chunk it, whatever that means.
            call h5pset_deflate_f(plist_id, 6, err)               ! szip compression lvl 6, whatever that means
                call h5dcreate_f(obj_id,trim(data_name),H5T_NATIVE_DOUBLE,dspace_id,dset_id,err)
                    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, buf, dims, err)
                    if ( present(enhet) ) call lo_h5_store_attribute(trim(enhet),dset_id,'unit',lo_status)
                    if ( present(dimensions) ) call lo_h5_store_attribute(trim(dimensions),dset_id,'dimensions',lo_status)
                call h5dclose_f(dset_id,err)
            call h5pclose_f(plist_id,err)
        call h5sclose_f(dspace_id,err)
        if ( present(error) ) error=err
end subroutine
subroutine store_double_2D_array_as_data(buf,obj_id,data_name,error,enhet,dimensions)
        !> double to store
        real(flyt), dimension(:,:), intent(in) :: buf
        !> identifier to attach the attribute to
        integer(HID_T), intent(in) :: obj_id
        !> name of the attribute
        character(len=*), intent(in) :: data_name
        !> status
        integer, optional, intent(out) :: error
        !> set the unit right away
        character(len=*), optional, intent(in) :: enhet
        !> specify what the dimensions are
        character(len=*), optional, intent(in) :: dimensions
        !
        integer :: err,rank
        integer(HID_T) :: dspace_id,dset_id,plist_id
        integer(HSIZE_T), dimension(2) :: dims,chunks
        ! sife of the array
        rank=2
        dims(1)=size(buf,1)
        dims(2)=size(buf,2)
        chunks=dims ! supposed to be bad, but ok for what I use it for, so far
        call h5screate_simple_f(rank,dims,dspace_id,err)      ! create dataset
            call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, err) ! create properties of dataset
            call h5pset_chunk_f(plist_id, 2, chunks, err)         ! chunk it, whatever that means.
            call h5pset_deflate_f(plist_id, 6, err)               ! szip compression lvl 6, whatever that means
                call h5dcreate_f(obj_id,trim(data_name),H5T_NATIVE_DOUBLE,dspace_id,dset_id,err,dcpl_id=plist_id)
                    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, buf, dims, err) ! write to it
                    if ( present(enhet) ) call lo_h5_store_attribute(trim(enhet),dset_id,'unit',lo_status)
                    if ( present(dimensions) ) call lo_h5_store_attribute(trim(dimensions),dset_id,'dimensions',lo_status)
                call h5dclose_f(dset_id,err)
            call h5pclose_f(plist_id,err)
        call h5sclose_f(dspace_id,err)
        if ( present(error) ) error=err
end subroutine
subroutine store_double_3D_array_as_data(buf,obj_id,data_name,error,enhet,dimensions)
        !> double to store
        real(flyt), dimension(:,:,:), intent(in) :: buf
        !> identifier to attach the attribute to
        integer(HID_T), intent(in) :: obj_id
        !> name of the attribute
        character(len=*), intent(in) :: data_name
        !> status
        integer, optional, intent(out) :: error
        !> set the unit right away
        character(len=*), optional, intent(in) :: enhet
        !> specify what the dimensions are
        character(len=*), optional, intent(in) :: dimensions
        !
        integer :: err,rank
        integer(HID_T) :: dspace_id,dset_id,plist_id
        integer(HSIZE_T), dimension(3) :: dims,chunks
        ! sife of the array
        rank=3
        dims(1)=size(buf,1)
        dims(2)=size(buf,2)
        dims(3)=size(buf,3)
        chunks=dims ! supposed to be bad, but ok for what I use it for, so far
        call h5screate_simple_f(rank,dims,dspace_id,err)      ! create dataset
            call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, err) ! create properties of dataset
            call h5pset_chunk_f(plist_id, 3, chunks, err)         ! chunk it, whatever that means.
            call h5pset_deflate_f(plist_id, 6, err)               ! szip compression lvl 6, whatever that means
                call h5dcreate_f(obj_id,trim(data_name),H5T_NATIVE_DOUBLE,dspace_id,dset_id,err,dcpl_id=plist_id) ! create the dataset
                    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, buf, dims, err) ! write to it
                    if ( present(enhet) ) call lo_h5_store_attribute(trim(enhet),dset_id,'unit',lo_status)
                    if ( present(dimensions) ) call lo_h5_store_attribute(trim(dimensions),dset_id,'dimensions',lo_status)
                call h5dclose_f(dset_id,err)
            call h5pclose_f(plist_id,err)
        call h5sclose_f(dspace_id,err)
        if ( present(error) ) error=err
end subroutine
subroutine store_double_4D_array_as_data(buf,obj_id,data_name,error,enhet,dimensions)
        !> double to store
        real(flyt), dimension(:,:,:,:), intent(in) :: buf
        !> identifier to attach the attribute to
        integer(HID_T), intent(in) :: obj_id
        !> name of the attribute
        character(len=*), intent(in) :: data_name
        !> status
        integer, optional, intent(out) :: error
        !> set the unit right away
        character(len=*), optional, intent(in) :: enhet
        !> specify what the dimensions are
        character(len=*), optional, intent(in) :: dimensions
        !
        integer :: err,rank
        integer(HID_T) :: dspace_id,dset_id,plist_id
        integer(HSIZE_T), dimension(5) :: dims,chunks
        ! sife of the array
        rank=4
        dims(1)=size(buf,1)
        dims(2)=size(buf,2)
        dims(3)=size(buf,3)
        dims(4)=size(buf,4)
        chunks=dims ! supposed to be bad, but ok for what I use it for, so far
        call h5screate_simple_f(rank,dims,dspace_id,err)      ! create dataset
            call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, err) ! create properties of dataset
            call h5pset_chunk_f(plist_id, 4, chunks, err)         ! chunk it, whatever that means.
            call h5pset_deflate_f(plist_id, 6, err)               ! szip compression lvl 6, whatever that means
                call h5dcreate_f(obj_id,trim(data_name),H5T_NATIVE_DOUBLE,dspace_id,dset_id,err,dcpl_id=plist_id) ! create the dataset
                    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, buf, dims, err) ! write to it
                    if ( present(enhet) ) call lo_h5_store_attribute(trim(enhet),dset_id,'unit',lo_status)
                    if ( present(dimensions) ) call lo_h5_store_attribute(trim(dimensions),dset_id,'dimensions',lo_status)
                call h5dclose_f(dset_id,err)
            call h5pclose_f(plist_id,err)
        call h5sclose_f(dspace_id,err)
        if ( present(error) ) error=err
end subroutine
subroutine store_double_5D_array_as_data(buf,obj_id,data_name,error,enhet,dimensions)
        !> double to store
        real(flyt), dimension(:,:,:,:,:), intent(in) :: buf
        !> identifier to attach the attribute to
        integer(HID_T), intent(in) :: obj_id
        !> name of the attribute
        character(len=*), intent(in) :: data_name
        !> status
        integer, optional, intent(out) :: error
        !> set the unit right away
        character(len=*), optional, intent(in) :: enhet
        !> specify what the dimensions are
        character(len=*), optional, intent(in) :: dimensions
        !
        integer :: err,rank
        integer(HID_T) :: dspace_id,dset_id,plist_id
        integer(HSIZE_T), dimension(5) :: dims,chunks
        ! sife of the array
        rank=5
        dims(1)=size(buf,1)
        dims(2)=size(buf,2)
        dims(3)=size(buf,3)
        dims(4)=size(buf,4)
        dims(5)=size(buf,5)
        chunks=dims ! supposed to be bad, but ok for what I use it for, so far
        call h5screate_simple_f(rank,dims,dspace_id,err)          ! create dataset
            call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, err) ! create properties of dataset
            call h5pset_chunk_f(plist_id, 5, chunks, err)         ! chunk it, whatever that means.
            call h5pset_deflate_f(plist_id, 6, err)               ! szip compression lvl 6, whatever that means
                call h5dcreate_f(obj_id,trim(data_name),H5T_NATIVE_DOUBLE,dspace_id,dset_id,err,dcpl_id=plist_id) ! create the dataset
                    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, buf, dims, err) ! write to it
                    if ( present(enhet) ) call lo_h5_store_attribute(trim(enhet),dset_id,'unit',lo_status)
                    if ( present(dimensions) ) call lo_h5_store_attribute(trim(dimensions),dset_id,'dimensions',lo_status)
                call h5dclose_f(dset_id,err)
            call h5pclose_f(plist_id,err)
        call h5sclose_f(dspace_id,err)
        if ( present(error) ) error=err
end subroutine

subroutine store_char_as_attribute(buf,obj_id,data_name,error)
    !> characters to store
    character(len=*), intent(in) :: buf
    !> identifier to attach the attribute to
    integer(HID_T), intent(in) :: obj_id
    !> name of the attribute
    character(len=*), intent(in) :: data_name
    !> status
    integer, optional, intent(out) :: error
    !
    integer :: err,rank
    integer(HID_T) :: dspace_id,aset_id,dtype_id
    integer(HSIZE_T), dimension(1) :: dims
    integer(HSIZE_T) :: strlen
    rank=0
    dims=0
    strlen=len(buf)
    !
    call H5Tcopy_f( H5T_FORTRAN_S1, dtype_id, err)
    call H5Tset_size_f(dtype_id, strlen, err)
    !
    call h5screate_simple_f(rank, dims, dspace_id, err)
    call h5acreate_f(obj_id, trim(data_name), dtype_id, dspace_id, aset_id, err)
    call h5awrite_f(aset_id,dtype_id,buf,dims,err)

    call h5aclose_f(aset_id, err)
    call h5sclose_f(dspace_id, err)
    call h5tclose_f(dtype_id, err)
    if ( present(error) ) error=err
end subroutine
subroutine store_float_as_attribute(buf,obj_id,attribute_name,error)
    !> double to store
    real(flyt), intent(in) :: buf
    !> identifier to attach the attribute to
    integer(HID_T), intent(in) :: obj_id
    !> name of the attribute
    character(len=*), intent(in) :: attribute_name
    !> status
    integer, optional, intent(out) :: error
    !
    integer :: err,rank
    integer(HID_T) :: attr_ds_id,attr_id
    integer(HSIZE_T), dimension(1) :: dims

    rank=1
    dims=1
    call h5screate_simple_f(rank,dims, attr_ds_id, err)
    call h5acreate_f(obj_id,trim(attribute_name),H5T_NATIVE_DOUBLE,attr_ds_id,attr_id,err)
    call h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, [buf], dims, err)
    call h5sclose_f(attr_ds_id,err)
    if ( present(error) ) error=err
end subroutine
subroutine store_float_1D_as_attribute(buf,obj_id,attribute_name,error)
    !> double to store
    real(flyt), intent(in), dimension(:) :: buf
    !> identifier to attach the attribute to
    integer(HID_T), intent(in) :: obj_id
    !> name of the attribute
    character(len=*), intent(in) :: attribute_name
    !> status
    integer, optional, intent(out) :: error
    !
    integer :: err,rank
    integer(HID_T) :: attr_ds_id,attr_id
    integer(HSIZE_T), dimension(1) :: dims

    rank=1
    dims(1)=size(buf,1)
    call h5screate_simple_f(rank,dims, attr_ds_id, err)
    call h5acreate_f(obj_id,trim(attribute_name),H5T_NATIVE_DOUBLE,attr_ds_id,attr_id,err)
    call h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, buf, dims, err)
    call h5sclose_f(attr_ds_id,err)
    if ( present(error) ) error=err
end subroutine
subroutine store_float_2D_as_attribute(buf,obj_id,attribute_name,error)
    !> double to store
    real(flyt), intent(in), dimension(:,:) :: buf
    !> identifier to attach the attribute to
    integer(HID_T), intent(in) :: obj_id
    !> name of the attribute
    character(len=*), intent(in) :: attribute_name
    !> status
    integer, optional, intent(out) :: error
    !
    integer :: err,rank
    integer(HID_T) :: attr_ds_id,attr_id
    integer(HSIZE_T), dimension(2) :: dims

    rank=2
    dims(1)=size(buf,1)
    dims(2)=size(buf,1)
    call h5screate_simple_f(rank,dims, attr_ds_id, err)
    call h5acreate_f(obj_id,trim(attribute_name),H5T_NATIVE_DOUBLE,attr_ds_id,attr_id,err)
    call h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, buf, dims, err)
    call h5sclose_f(attr_ds_id,err)
    if ( present(error) ) error=err
end subroutine

subroutine store_int_as_attribute(buf,obj_id,attribute_name,error)
    !> integer to store
    integer, intent(in) :: buf
    !> identifier to attach the attribute to
    integer(HID_T), intent(in) :: obj_id
    !> name of the attribute
    character(len=*), intent(in) :: attribute_name
    !> status
    integer, optional, intent(out) :: error
    !
    integer :: err,rank
    integer(HID_T) :: attr_ds_id,attr_id
    integer(HSIZE_T), dimension(1) :: dims
    integer, dimension(1) :: dumbuf
    dumbuf(1)=buf

    rank=1
    dims=1
    call h5screate_simple_f(rank,dims, attr_ds_id, err)
    call h5acreate_f(obj_id,trim(attribute_name),H5T_NATIVE_INTEGER,attr_ds_id,attr_id,err)
    call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, dumbuf, dims, err)
    call h5sclose_f(attr_ds_id,err)
    if ( present(error) ) error=err
end subroutine
subroutine store_logical_as_attribute(buf,obj_id,attribute_name,error)
    !> logical to store
    logical, intent(in) :: buf
    !> identifier to attach the attribute to
    integer(HID_T), intent(in) :: obj_id
    !> name of the attribute
    character(len=*), intent(in) :: attribute_name
    !> status
    integer, optional, intent(out) :: error
    !
    integer :: err,i,rank
    integer(HID_T) :: attr_ds_id,attr_id
    integer(HSIZE_T), dimension(1) :: dims
    ! hdf5 is stupid and can not store logicals, so I store it as an int, and when
    ! it's read it's converted back again, or something.
    if ( buf ) then
        i=1
    else
        i=0
    endif
    rank=1
    dims=1
    call h5screate_simple_f(rank,dims, attr_ds_id, err)
    call h5acreate_f(obj_id,trim(attribute_name),H5T_NATIVE_INTEGER,attr_ds_id,attr_id,err)
    call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, [i], dims, err)
    call h5sclose_f(attr_ds_id,err)
    if ( present(error) ) error=err
end subroutine

subroutine read_int_as_attribute(buf,obj_id,attribute_name,error)
    !> integer to store
    integer, intent(out) :: buf
    !> identifier to attach the attribute to
    integer(HID_T), intent(in) :: obj_id
    !> name of the attribute
    character(len=*), intent(in) :: attribute_name
    !> status
    integer, optional, intent(out) :: error
    !
    integer :: err
    integer(HID_T) :: attr_id
    integer(HSIZE_T), dimension(1) :: dims
    dims=0
    call h5aopen_f(obj_id,trim(attribute_name),attr_id,err)
    call h5aread_f(attr_id, H5T_NATIVE_INTEGER, buf, dims, err)
    call h5aclose_f(attr_id,err)
    if ( present(error) ) error=err
end subroutine
subroutine read_logical_as_attribute(buf,obj_id,attribute_name,error)
    !> integer to store
    logical, intent(out) :: buf
    !> identifier to attach the attribute to
    integer(HID_T), intent(in) :: obj_id
    !> name of the attribute
    character(len=*), intent(in) :: attribute_name
    !> status
    integer, optional, intent(out) :: error
    !
    integer :: err,i
    integer(HID_T) :: attr_id
    integer(HSIZE_T), dimension(1) :: dims
    dims=0
    call h5aopen_f(obj_id,trim(attribute_name),attr_id,err)
    call h5aread_f(attr_id, H5T_NATIVE_INTEGER, i, dims, err)
    call h5aclose_f(attr_id,err)
    if ( i .eq. 1 ) then
        buf=.true.
    else
        buf=.false.
    endif
    if ( present(error) ) error=err
end subroutine
subroutine read_float_as_attribute(buf,obj_id,attribute_name,error)
    !> integer to store
    real(flyt), intent(out) :: buf
    !> identifier to attach the attribute to
    integer(HID_T), intent(in) :: obj_id
    !> name of the attribute
    character(len=*), intent(in) :: attribute_name
    !> status
    integer, optional, intent(out) :: error
    !
    integer :: err
    integer(HID_T) :: attr_id
    integer(HSIZE_T), dimension(1) :: dims
    dims=0
    call h5aopen_f(obj_id,trim(attribute_name),attr_id,err)
    call h5aread_f(attr_id, H5T_NATIVE_DOUBLE, buf, dims, err)
    call h5aclose_f(attr_id,err)
    if ( present(error) ) error=err
end subroutine

end module
