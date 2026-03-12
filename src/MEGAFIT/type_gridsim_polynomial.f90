
!> transform coordinates for a better fit
subroutine coordinate_transformation(gs,xin,xout)
    !> grid simulation thing
    class(lo_gridsim), intent(in) :: gs
    !> raw coordinate, such as V,T,P,eta
    real(flyt), dimension(gs%ndim), intent(in) :: xin
    !> scaled coordinate, something magical.
    real(flyt), dimension(gs%ndim), intent(out) :: xout

    real(flyt) :: temperature,volume

    ! Per default, just pass it as-is
    xout=xin
    ! Scale the temperature with strange function?
    if ( gs%info%temperature_scale .gt. 0.0_flyt .and. gs%info%dim_temperature .gt. 0 ) then
        temperature=xin( gs%info%dim_temperature )
        xout( gs%info%dim_temperature )=tempscaler( temperature, gs%info%temperature_scale )
    endif

    ! Scale temperature to 0-1 interval
    if ( gs%info%dim_temperature .gt. 0 ) then
        xout( gs%info%dim_temperature )=xout( gs%info%dim_temperature )/maxval(gs%grid_coordinates( gs%info%dim_temperature,: ))
    endif

    ! Change volume to that x-thingy a Birch-Murnaghan uses?
    if ( gs%info%dim_volume .gt. 0 ) then
        volume=xin( gs%info%dim_volume )
        xout( gs%info%dim_volume )=(1.0_flyt/volume)**(2.0_flyt/3.0_flyt)
    endif

    ! Change volume to pressure instead?
    !select type(gs%eos)
    !class is(lo_eos_1d)
    !class is(lo_eos_2d)
    !end select
end subroutine

! rescales temperatures in a magical way
pure function tempscaler(T,Tref) result(x)
    real(flyt), intent(in) :: T,Tref

    real(flyt) :: x

    if ( T .lt. Tref*lo_tol ) then
        x=Tref*0.5_flyt
    else
        x= Tref*0.5_flyt + T*2.0_flyt/(exp(Tref/T)+1.0_flyt)
    endif
end function

! ! rescales temperatures in a magical way
! pure function magtempscaler(T,Tref) result(x)
!     real(flyt), intent(in) :: T,Tref
!
!     real(flyt) :: x
!
!     if ( T .lt. Tref*lo_tol ) then
!         x=Tref*0.5_flyt
!     else
!         x= Tref*0.5_flyt + T*2.0_flyt/(exp(Tref/T)+1.0_flyt)
!     endif
!     x=x*lo_kB_Hartree
! end function
