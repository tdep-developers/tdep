!> Tools to symmetrize things, mostly from
!> https://github.com/ollehellman/aims_wrappers/blob/master/src/relax_with_symmetry_constraints/type_irredcell.f90

module helper
use konstanter, only: r8, lo_sqtol, lo_exitcode_symmetry, lo_status
use gottochblandat, only: open_file

implicit none

contains

subroutine write_stress_3x3(stress, stress_err, header)
    real(r8), dimension(3, 3), intent(in) :: stress, stress_err
    character(len=*), intent(in) :: header
    character(len=1000) :: opf

    write (*, *) header
    opf = "(1X,'  average xx xy xz (GPa): ',3(F12.5,' '),' +/- ',3(F12.5,' ') )"
    write (*, opf) stress(:, 1), stress_err(:, 1)
    opf = "(1X,'  average yx yy yz (GPa): ',3(F12.5,' '),' +/- ',3(F12.5,' ') )"
    write (*, opf) stress(:, 2), stress_err(:, 2)
    opf = "(1X,'  average zx zy zz (GPa): ',3(F12.5,' '),' +/- ',3(F12.5,' ') )"
    write (*, opf) stress(:, 3), stress_err(:, 3)
    write (*, *)
end subroutine

subroutine write_stress_voigt(stress, stress_std, stress_err, header)
    real(r8), dimension(3, 3), intent(in) :: stress, stress_std, stress_err
    character(len=*), intent(in) :: header
    character(len=1000) :: opf

    write (*, *) header
    opf = "(1X,'  average xx yy zz xz yz xy (GPa): ',6(F12.5,' '),' +/- ')"
    associate (s => stress)
        write (*, opf) s(1, 1), s(2, 2), s(3, 3), s(3, 2), s(3, 1), s(2, 1)
    end associate

    opf = "(1X,'   stddev xx yy zz xz yz xy (GPa): ',6(F12.5,' '))"
    associate (s => stress_std)
        write (*, opf) s(1, 1), s(2, 2), s(3, 3), s(3, 2), s(3, 1), s(2, 1)
    end associate

    opf = "(1X,'   stderr xx yy zz xz yz xy (GPa): ',6(F12.5,' '))"
    associate (s => stress_err)
        write (*, opf) s(1, 1), s(2, 2), s(3, 3), s(3, 2), s(3, 1), s(2, 1)
    end associate

    write (*, *)
end subroutine

subroutine write_file_stress_voigt_time(stress, file, header)
    real(r8), dimension(:, :, :), intent(in) :: stress
    character(len=*), intent(in) :: header, file
    integer :: u, i

    u = open_file('out', file)
    write (u, '(a)') header
    ! write (u, '(a)') 'timestep,stress_xx,stress_yy,stress_zz,stress_yz,stress_xz,stress_xy'
    write (u, '(a)') 'timestep,xx,yy,zz,yz,xz,xy'
    associate (s => stress)
        do i = 1, size(s, dim=3)
            write (u, '(I5,6(",",ES15.8))') i, &
                s(1, 1, i), s(2, 2, i), s(3, 3, i), &
                s(2, 3, i), s(1, 3, i), s(1, 2, i)
        end do
    end associate
    close (u)

end subroutine

subroutine write_file_stress_voigt(stress, stress_std, stress_err, file, header)
    real(r8), dimension(3, 3), intent(in) :: stress, stress_std, stress_err
    character(len=*), intent(in) :: header, file
    integer :: u

    u = open_file('out', file)
    write (u, '(a)') header
    write (u, '(a)') '# stress (Voigt xx, yy, zz, yz, xz, xy)'
    associate (s => stress)
        write (u, '(1X,6(1X,ES15.8))') s(1, 1), s(2, 2), s(3, 3), s(3, 2), s(3, 1), s(2, 1)
    end associate
    write (u, '(a)') '# std. dev. (Voigt xx, yy, zz, yz, xz, xy)'
    associate (s => stress_std)
        write (u, '(1X,6(1X,ES15.8))') s(1, 1), s(2, 2), s(3, 3), s(3, 2), s(3, 1), s(2, 1)
    end associate
    write (u, '(a)') '# std. err. (Voigt xx, yy, zz, yz, xz, xy)'
    associate (s => stress_err)
        write (u, '(1X,6(1X,ES15.8))') s(1, 1), s(2, 2), s(3, 3), s(3, 2), s(3, 1), s(2, 1)
    end associate
    close (u)

end subroutine

end module
