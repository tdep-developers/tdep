
!> size in memory, in bytes
function memsize_polynomial(pl) result(mem)
    class(lo_polynomial), intent(in) :: pl
    integer :: mem

    mem=0
    mem=mem+storage_size(pl)
    if ( allocated( pl%exponents ) ) mem=mem+size(pl%exponents)*storage_size(pl%exponents)
    if ( allocated( pl%coeffM ) ) mem=mem+size(pl%coeffM)*storage_size(pl%coeffM)
    if ( allocated( pl%variablename ) ) mem=mem+size(pl%variablename)*storage_size(pl%variablename)
    mem=mem/8
end function

!> size in memory, in bytes
function memsize_interpolation(ip) result(mem)
    class(lo_grid_interpolation), intent(in) :: ip
    integer :: mem
    mem=0
    mem=mem+storage_size(ip)
    if ( allocated( ip%gridcoord )) mem=mem+size(ip%gridcoord)*storage_size(ip%gridcoord)
    if ( allocated( ip%gridvals )) mem=mem+size(ip%gridvals)*storage_size(ip%gridvals)
    if ( allocated( ip%coord_shift )) mem=mem+size(ip%coord_shift)*storage_size(ip%coord_shift)
    if ( allocated( ip%coord_scale )) mem=mem+size(ip%coord_scale)*storage_size(ip%coord_scale)
    if ( allocated( ip%wM )) mem=mem+size(ip%wM)*storage_size(ip%wM)
    if ( allocated( ip%wV )) mem=mem+size(ip%wV)*storage_size(ip%wV)
    if ( allocated( ip%weights )) mem=mem+size(ip%weights)*storage_size(ip%weights)
    mem=mem/8+ip%pl%size_in_mem()
end function
