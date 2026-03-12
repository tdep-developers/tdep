
!> report memory usage to stdout
subroutine memreport(gs,mw)
    !> grid with stuff
    class(lo_gridsim), intent(in) :: gs
    !> mpi helper
    type(lo_mpi_helper), intent(in) :: mw

    integer, dimension(mw%n) :: bs

write(*,*) 'REVISE ME GRIDSIM MEM'
stop
!    if ( mw%talk ) then
!        write(*,*) ''
!        write(*,*) 'GRIDSIM MEMORY REPORT:'
!    endif
!    bs=0
!    bs(mw%r+1)=gs%raw%size_in_mem()
!    call mpi_allreduce(MPI_IN_PLACE,bs,mw%n,MPI_INTEGER,MPI_SUM,mw%comm,mw%error)
!    if ( mw%talk ) write(*,"(1X,a,3(3X,F14.6))") '        raw:',lo_mean(bs/1024.0_flyt**2),minval(bs/1024.0_flyt**2),maxval(bs/1024.0_flyt**2)
!
!    bs=0
!    bs(mw%r+1)=gs%ref(1)%size_in_mem()*size(gs%ref)
!    call mpi_allreduce(MPI_IN_PLACE,bs,mw%n,MPI_INTEGER,MPI_SUM,mw%comm,mw%error)
!    if ( mw%talk ) write(*,"(1X,a,3(3X,F14.6))") '        ref:',lo_mean(bs/1024.0_flyt**2),minval(bs/1024.0_flyt**2),maxval(bs/1024.0_flyt**2)
!
!    bs=0
!    bs(mw%r+1)=gs%structure%size_in_mem()
!    call mpi_allreduce(MPI_IN_PLACE,bs,mw%n,MPI_INTEGER,MPI_SUM,mw%comm,mw%error)
!    if ( mw%talk ) write(*,"(1X,a,3(3X,F14.6))") '  structure:',lo_mean(bs/1024.0_flyt**2),minval(bs/1024.0_flyt**2),maxval(bs/1024.0_flyt**2)
!
!    bs=0
!    bs(mw%r+1)=gs%energy%size_in_mem()
!    call mpi_allreduce(MPI_IN_PLACE,bs,mw%n,MPI_INTEGER,MPI_SUM,mw%comm,mw%error)
!    if ( mw%talk ) write(*,"(1X,a,3(3X,F14.6))") '     energy:',lo_mean(bs/1024.0_flyt**2),minval(bs/1024.0_flyt**2),maxval(bs/1024.0_flyt**2)
!
!    bs=0
!    bs(mw%r+1)=gs%polar%size_in_mem()+gs%pair%size_in_mem()+gs%triplet%size_in_mem()+gs%quartet%size_in_mem()
!    call mpi_allreduce(MPI_IN_PLACE,bs,mw%n,MPI_INTEGER,MPI_SUM,mw%comm,mw%error)
!    if ( mw%talk ) write(*,"(1X,a,3(3X,F14.6))") '      coeff:',lo_mean(bs/1024.0_flyt**2),minval(bs/1024.0_flyt**2),maxval(bs/1024.0_flyt**2)
end subroutine

!> size in memory, in bytes
function memsize_raw(rw) result(mem)
    class(lo_gridsim_rawdata), intent(in) :: rw
    integer :: mem
    mem=0
!    mem=mem+storage_size(rw)
!    if ( allocated( rw%u ) ) mem=mem+size(rw%u)*storage_size(rw%u)
!    if ( allocated( rw%f ) ) mem=mem+size(rw%f)*storage_size(rw%f)
!    if ( allocated( rw%f0 ) ) mem=mem+size(rw%f0)*storage_size(rw%f0)
!    if ( allocated( rw%m ) ) mem=mem+size(rw%m)*storage_size(rw%m)
!    if ( allocated( rw%e ) ) mem=mem+size(rw%e)*storage_size(rw%e)
!    if ( allocated( rw%e_polar ) ) mem=mem+size(rw%e_polar)*storage_size(rw%e_polar)
!    if ( allocated( rw%e_secondorder ) ) mem=mem+size(rw%e_secondorder)*storage_size(rw%e_secondorder)
!    if ( allocated( rw%e_thirdorder ) ) mem=mem+size(rw%e_thirdorder)*storage_size(rw%e_thirdorder)
!    if ( allocated( rw%e_fourthorder ) ) mem=mem+size(rw%e_fourthorder)*storage_size(rw%e_fourthorder)
!    if ( allocated( rw%gridind ) ) mem=mem+size(rw%gridind)*storage_size(rw%gridind)
!    mem=mem/8
end function

!> size in memory, in bytes
function memsize_ref(rf) result(mem)
    class(lo_gridsim_refdata), intent(in) :: rf
    integer :: mem
    mem=0
!    mem=mem+storage_size(rf)
!    if ( allocated( rf%unitcell_positions ) ) mem=mem+size(rf%unitcell_positions)*storage_size(rf%unitcell_positions)
!    if ( allocated( rf%supercell_positions ) ) mem=mem+size(rf%supercell_positions)*storage_size(rf%supercell_positions)
!    if ( allocated( rf%supercell_atomic_numbers ) ) mem=mem+size(rf%supercell_atomic_numbers)*storage_size(rf%supercell_atomic_numbers)
!    if ( allocated( rf%unitcell_atomic_numbers ) ) mem=mem+size(rf%unitcell_atomic_numbers)*storage_size(rf%unitcell_atomic_numbers)
!    if ( allocated( rf%born_effective_charges ) ) mem=mem+size(rf%born_effective_charges)*storage_size(rf%born_effective_charges)
!    if ( allocated( rf%dipole_forceconstant ) ) mem=mem+size(rf%dipole_forceconstant)*storage_size(rf%dipole_forceconstant)
!    if ( allocated( rf%constr_pair ) ) mem=mem+size(rf%constr_pair)*storage_size(rf%constr_pair)
!    mem=mem/8
end function

!> size in memory, in bytes
function memsize_coeff(cf) result(mem)
    class(lo_gridsim_coeff), intent(in) :: cf
    integer :: mem
    mem=0
!    mem=mem+storage_size(cf)
!    if ( allocated( cf%coeff ) ) mem=mem+size(cf%coeff)*storage_size(cf%coeff)
!    if ( allocated( cf%constraints ) ) mem=mem+size(cf%constraints)*storage_size(cf%constraints)
!    if ( allocated( cf%rsquare ) ) mem=mem+size(cf%rsquare)*storage_size(cf%rsquare)
!    mem=mem/8
end function

!> size in memory, in bytes
function memsize_structure(st) result(mem)
    class(lo_gridsim_structure), intent(in) :: st
    integer :: mem
    mem=0
!    mem=mem+storage_size(st)
!!    if ( allocated( st%xi ) ) mem=mem+size(st%xi)*storage_size(st%xi)
!!    if ( allocated( st%vol ) ) mem=mem+size(st%vol)*storage_size(st%vol)
!!    if ( allocated( st%latticevectors ) ) mem=mem+size(st%lv)*storage_size(st%ilv)
!!    if ( allocated( st%latticevectors ) ) mem=mem+size(st%lv0)*storage_size(st%ilv0)
!    if ( allocated( st%atomic_numbers ) ) mem=mem+size(st%atomic_numbers)*storage_size(st%atomic_numbers)
!    if ( allocated( st%r0 ) ) mem=mem+size(st%r0)*storage_size(st%r0)
!    if ( allocated( st%coeffM ) ) mem=mem+size(st%coeffM)*storage_size(st%coeffM)
!    if ( allocated( st%wV ) ) mem=mem+size(st%wV)*storage_size(st%wV)
!    if ( allocated( st%wU ) ) mem=mem+size(st%wU)*storage_size(st%wU)
!    if ( allocated( st%wM ) ) mem=mem+size(st%wM)*storage_size(st%wM)    
!    mem=mem/8 ! to bytes
!    mem=mem+st%ipint%size_in_mem()
!    mem=mem+st%iplv%size_in_mem()
!    mem=mem+st%ipvol%size_in_mem()
end function

!> size in memory, in bytes
function memsize_energy(en) result(mem)
    class(lo_gridsim_energy), intent(in) :: en
    integer :: mem
    mem=0
!    mem=mem+storage_size(en)
!    if ( allocated( en%xi ) ) mem=mem+size(en%xi)*storage_size(en%xi)
!    if ( allocated( en%x0 ) ) mem=mem+size(en%x0)*storage_size(en%x0)
!    if ( allocated( en%ixs ) ) mem=mem+size(en%ixs)*storage_size(en%ixs)
!    if ( allocated( en%energy ) ) mem=mem+size(en%energy)*storage_size(en%energy)
!    mem=mem/8+en%pl%size_in_mem()
end function

!> size in memory, in bytes
function memsize_polar(plr) result(mem)
    class(lo_gridsim_polar), intent(in) :: plr
    integer :: mem
    mem=0
!    mem=mem+storage_size(plr)
!    if ( allocated(plr%constraints) ) mem=mem+size(plr%constraints)*storage_size(plr%constraints)
!    if ( allocated(plr%rsquare) ) mem=mem+size(plr%rsquare)*storage_size(plr%rsquare)
!    mem=mem/8+plr%ipZ%size_in_mem()+plr%ipeps%size_in_mem()
end function

!> size in memory, in bytes
function memsize_pair(pair) result(mem)
    class(lo_gridsim_pair), intent(in) :: pair
    integer :: mem
    mem=0
!    mem=mem+storage_size(pair)
!    if ( allocated(pair%rsquare) ) mem=mem+size(pair%rsquare)*storage_size(pair%rsquare)
!    mem=mem/8+pair%ip%size_in_mem()
end function
