
!> write as a POV-Ray scene
subroutine write_to_povray(ps, uc, filename, theta, phi, modes_plusminus_clr_alpha, quality)
    !> surface
    class(lo_phasespacesurface), intent(in) :: ps
    !> unitcell
    type(lo_crystalstructure), intent(in) :: uc
    !> filename
    character(len=*), intent(in) :: filename
    !> theta angle
    real(flyt), intent(in) :: theta
    !> phi angle
    real(flyt), intent(in) :: phi
    !> which surfaces to look at, how to color them, alpha channel
    integer, dimension(:, :), intent(in) :: modes_plusminus_clr_alpha
    !> how nice to render
    integer, intent(in) :: quality

    integer :: i, j, k, l, u, ii, jj, kk, ntri, npts, nsurf, srf
    integer :: b1, b2, b3
    real(flyt) :: f0, cylwidth
    real(flyt), dimension(3) :: v0, v1, v2, g0, g1, g2
    real(flyt), dimension(3) :: v_cam, v_up, v_light1, v_light2
    real(flyt), dimension(:, :), allocatable :: grad
    logical :: plus

    write (*, *) ''
    write (*, *) 'Writing surfaces to POV-Ray'

    ! Set the rendering options
    u = open_file('out', trim(filename)//'.ini')
    write (u, *) 'Input_File_Name='//trim(filename)//'.pov'
    write (u, *) 'Quality = ', tochar(quality)
    write (u, *) 'Antialias = 1'
    write (u, *) 'Antialias_Depth = 6'
    write (u, *) 'Height = 1080'
    write (u, *) 'Width = 1920'
    close (u)

    u = open_file('out', trim(filename)//'.pov')
    ! First I write the povray header, setting up the scene.
    ! set some colors
    write (u, *) '#include "colors.inc"'
    write (u, *) '#declare BZCOLOR = rgbt<0.6,0.6,0.6,0.98>;'
    write (u, *) 'background { color White }'
    ! figure out where to put the camera
    v0(1) = sin(theta*lo_pi/180.0_flyt)*cos(phi*lo_pi/180.0_flyt)
    v0(2) = sin(theta*lo_pi/180.0_flyt)
    v0(3) = cos(theta*lo_pi/180.0_flyt)
    ! scale the location of the camera such at the distance from origin is 2.5 times the radius of the
    ! bounding sphere of the BZ, seems like a robust choice.
    v0 = v0*uc%bz%rmax/norm2(v0)
    v_cam = v0*2.5
    write (u, *) 'camera {'
    write (u, *) '  location <', tochar(v_cam(1)), ',', tochar(v_cam(2)), ',', tochar(v_cam(3)), '>'
    write (u, *) '  up <0,0.9,0>'
    write (u, *) '  right <1.6,0,0>  // right,up -> 16:9'
    write (u, *) '  look_at <0, 0, 0>'
    write (u, *) '}'
    v_up = [0, 0, 1]*1.0_flyt
    v0 = lo_cross(v_up, v_cam)
    v0 = v0/norm2(v0) + 0.5_flyt*v_cam/norm2(v0) + v_up
    v_light1 = v0
    v_light1 = v_light1*10.0_flyt/norm2(v_light1)
    ! Set the light source
    write (u, *) 'light_source {'
    write (u, *) '  <', tochar(v_light1(1)), ',', tochar(v_light1(2)), ',', tochar(v_light1(3)), '>'
    write (u, *) '  color White'
    write (u, *) '  area_light <.5, 0, 0>, <0, 0, .5>, 5, 5'
    write (u, *) '  adaptive 1'
    write (u, *) '  jitter'
    write (u, *) '}'
    write (u, *) 'light_source {'
    write (u, *) '  <10, 20, 10>'
    write (u, *) '  color Blue'
    write (u, *) '  area_light <5, 0, 0>, <0, 0, 5>, 5, 5'
    write (u, *) '  adaptive 1'
    write (u, *) '  jitter'
    write (u, *) '}'
    write (u, *) ''

    ! Now the header is done, with cameras and lights. Next up is defining the BZ, which
    ! I do as both a polyhedron that is very transparent. Then I put tiny cylinders along the
    ! edges. These cylinder reflect no light, and cast no shadow, so it looksjust like black lines.

    ! define the BZ
!        do i=1,uc%bz%nfaces
!            write(u,*) 'polygon {'
!            j=uc%bz%face(i)%n+1
!            write(u,*) '    ',tochar(j),','
!            write(u,"(4X)",advance='no')
!            do j=1,uc%bz%face(i)%n
!                write(u,"('<',F11.6,',',F11.6,',',F11.6,'>,')",advance='no') uc%bz%r(:,uc%bz%face(i)%ind(j))
!            enddo
!            write(u,"('<',F11.6,',',F11.6,',',F11.6,'>')") uc%bz%r(:,uc%bz%face(i)%ind(1))
!            write(u,*) '    texture { pigment { color BZCOLOR }}'
!            write(u,*) '}'
!        enddo
    ! draw some cylinders along the edges
    cylwidth = uc%bz%rmax/500.0_flyt
    do i = 1, uc%bz%nedges
        write (u, *) 'cylinder {'
        v0 = uc%bz%r(:, uc%bz%edge(i)%i1)
        v1 = uc%bz%r(:, uc%bz%edge(i)%i2)
        write (u, "(4X,'<',F11.6,',',F11.6,',',F11.6,'>,<',F11.6,',',F11.6,',',F11.6,'>,',F11.6)") v0, v1, cylwidth
        write (u, *) '    texture { pigment { color Black } }'
        write (u, *) '    no_shadow'
        write (u, *) '}'
    end do

    ! Now the idea is to dump triangles defining the phasespace surface. There are nb^3*2 phasespace surfaces for
    ! each q, so we have to define which one to look for. That's done in the input with 6 integers per surface.
    nsurf = size(modes_plusminus_clr_alpha, 2)
    do srf = 1, nsurf
        ! get the specification for this surface
        if (modes_plusminus_clr_alpha(4, srf) .eq. 0) then
            plus = .true.
        else
            plus = .false.
        end if
        b1 = modes_plusminus_clr_alpha(1, srf)
        b2 = modes_plusminus_clr_alpha(2, srf)
        b3 = modes_plusminus_clr_alpha(3, srf)
        ! define the color
        i = modes_plusminus_clr_alpha(5, srf)
        f0 = (100 - modes_plusminus_clr_alpha(6, srf))/100.0_flyt
        select case (i)
        case (1)
            write (u, *) '#declare SFCOLOR'//tochar(srf)//' = rgbft<0.8,0.1,0.1,0,'//tochar(f0)//'>;'
        case (2)
            write (u, *) '#declare SFCOLOR'//tochar(srf)//' = rgbft<0.1,0.8,0.1,0,'//tochar(f0)//'>;'
        case (3)
            write (u, *) '#declare SFCOLOR'//tochar(srf)//' = rgbft<0.3,0.3,0.9,0,'//tochar(f0)//'>;'
        end select

        ! Dump the triangles that make up the surface
        if (plus) then
            ntri = ps%plus(b1, b2, b3)%ntri
            npts = ps%plus(b1, b2, b3)%npts
            if (ntri .gt. 0) then
                ! print the mesh
                write (u, *) 'mesh {'
                do i = 1, ntri
                    ii = ps%plus(b1, b2, b3)%tri(1, i)
                    jj = ps%plus(b1, b2, b3)%tri(2, i)
                    kk = ps%plus(b1, b2, b3)%tri(3, i)
                    v0 = ps%plus(b1, b2, b3)%pts(:, ii)
                    v1 = ps%plus(b1, b2, b3)%pts(:, jj)
                    v2 = ps%plus(b1, b2, b3)%pts(:, kk)
                    g0 = ps%plus(b1, b2, b3)%grad(:, ii)
                    g1 = ps%plus(b1, b2, b3)%grad(:, jj)
                    g2 = ps%plus(b1, b2, b3)%grad(:, kk)
                    !write(u,"(1X,A)",advance='no') 'triangle { '
                    !write(u,"('<',F11.6,',',F11.6,',',F11.6,'>,')",advance='no') v0
                    !write(u,"('<',F11.6,',',F11.6,',',F11.6,'>,')",advance='no') v1
                    !write(u,"('<',F11.6,',',F11.6,',',F11.6,'>')",advance='no') v2
                    !write(u,*) '}'
                    write (u, "(1X,A)", advance='no') 'smooth_triangle { '
                    write (u, "('<',F11.6,',',F11.6,',',F11.6,'>,')", advance='no') v0
                    write (u, "('<',F11.6,',',F11.6,',',F11.6,'>,')", advance='no') g0
                    write (u, "('<',F11.6,',',F11.6,',',F11.6,'>')", advance='no') v1
                    write (u, "('<',F11.6,',',F11.6,',',F11.6,'>,')", advance='no') g1
                    write (u, "('<',F11.6,',',F11.6,',',F11.6,'>,')", advance='no') v2
                    write (u, "('<',F11.6,',',F11.6,',',F11.6,'>')", advance='no') g2
                    write (u, *) '}'
                end do
                write (u, *) '    texture { pigment { color SFCOLOR'//tochar(srf)//' }}'
                write (u, *) '    finish { ambient .4 diffuse .6 phong .2 }'
                write (u, *) '}'
            end if
        else
            ntri = ps%minus(b1, b2, b3)%ntri
            npts = ps%minus(b1, b2, b3)%npts
            if (ntri .gt. 0) then
                ! print the mesh
                write (u, *) 'mesh {'
                do i = 1, ntri
                    ii = ps%minus(b1, b2, b3)%tri(1, i)
                    jj = ps%minus(b1, b2, b3)%tri(2, i)
                    kk = ps%minus(b1, b2, b3)%tri(3, i)
                    v0 = ps%minus(b1, b2, b3)%pts(:, ii)
                    v1 = ps%minus(b1, b2, b3)%pts(:, jj)
                    v2 = ps%minus(b1, b2, b3)%pts(:, kk)
                    g0 = ps%minus(b1, b2, b3)%grad(:, ii)
                    g1 = ps%minus(b1, b2, b3)%grad(:, jj)
                    g2 = ps%minus(b1, b2, b3)%grad(:, kk)
                    !write(u,"(1X,A)",advance='no') 'triangle { '
                    !write(u,"('<',F11.6,',',F11.6,',',F11.6,'>,')",advance='no') v0
                    !write(u,"('<',F11.6,',',F11.6,',',F11.6,'>,')",advance='no') v1
                    !write(u,"('<',F11.6,',',F11.6,',',F11.6,'>')",advance='no') v2
                    !write(u,*) '}'
                    write (u, "(1X,A)", advance='no') 'smooth_triangle { '
                    write (u, "('<',F11.6,',',F11.6,',',F11.6,'>,')", advance='no') v0
                    write (u, "('<',F11.6,',',F11.6,',',F11.6,'>,')", advance='no') g0
                    write (u, "('<',F11.6,',',F11.6,',',F11.6,'>')", advance='no') v1
                    write (u, "('<',F11.6,',',F11.6,',',F11.6,'>,')", advance='no') g1
                    write (u, "('<',F11.6,',',F11.6,',',F11.6,'>,')", advance='no') v2
                    write (u, "('<',F11.6,',',F11.6,',',F11.6,'>')", advance='no') g2
                    write (u, *) '}'
                end do
                write (u, *) '    texture { pigment { color SFCOLOR'//tochar(srf)//' }}'
                write (u, *) '    finish { ambient .4 diffuse .6 phong .2 }'
                write (u, *) '}'
            end if
        end if
        !
    end do

    close (u)

end subroutine
