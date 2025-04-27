
subroutine compute_isotope_scattering(il, sr, qp, dr, uc, temperature, &
                                      g0, integrationtype, smearing, mw, mem)
    !> The local point
    integer, intent(in) :: il
    !> The scattering amplitudes
    type(lo_scattering_rates), intent(inout) :: sr
    !> The q-point mesh
    class(lo_qpoint_mesh), intent(in) :: qp
    !> phonon dispersions
    type(lo_phonon_dispersions), intent(inout) :: dr
    !> structure
    type(lo_crystalstructure), intent(in) :: uc
    !> The temperature
    real(r8), intent(in) :: temperature
    !> The linewidth for this mode
    real(r8), intent(inout) :: g0
    !> what kind of integration are we doing
    integer, intent(in) :: integrationtype
    !> The smearing width
    real(r8), intent(in) :: smearing
    !> MPI helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem

    ! Eigenvectors
    complex(r8), dimension(uc%na*3, 2) :: egviso
    !> For the broadening calculation
    real(r8), dimension(3) :: allsig
    ! prefactor and phonon buffers
    real(r8) :: om1, om2, sigma, psisq, prefactor, f0
    ! Integers for do loops
    integer :: q1, b1, q2, b2, i, niso

    q1 = sr%my_qpoints(il)
    b1 = sr%my_modes(il)
    om1 = dr%iq(q1)%omega(b1)
    egviso(:, 1) = dr%iq(q1)%egv(:, b1)

    do q2 = 1, qp%n_full_point
        prefactor = isotope_prefactor*qp%ap(q2)%integration_weight
        do b2 = 1, dr%n_mode
            om2 = dr%aq(q2)%omega(b2)
            if (om2 .lt. lo_freqtol) cycle

            select case (integrationtype)
            case (1)
                sigma = lo_frequency_THz_to_Hartree*smearing
            case (2)
                sigma = sqrt(sr%sigsq(q1, b1) + &
                             sr%sigsq(qp%ap(q2)%irreducible_index, b2))
            case (6)
                allsig = matmul(dr%aq(q2)%vel(:, b2), sr%reclat)**2
                sigma = sqrt(maxval(allsig) * 0.5_r8)
                sigma = max(sigma, dr%default_smearing(b2) * 0.25_r8)
!               sigma = qp%smearingparameter(dr%aq(q2)%vel(:, b2), &
!                                            dr%default_smearing(b2), smearing)
            end select

            i = (q2 - 1)*dr%n_mode + b2

            egviso(:, 2) = dr%aq(q2)%egv(:, b2)

            psisq = isotope_scattering_strength(uc, egviso)*prefactor

            f0 = psisq*om1*om2*lo_gauss(om1, om2, sigma)
            g0 = g0 + f0
            sr%Xi(il, i) = sr%Xi(il, i) + f0*om2/om1
        end do
    end do
end subroutine

real(r8) function isotope_scattering_strength(uc, egv)
    type(lo_crystalstructure), intent(in) :: uc
    complex(r8), dimension(:, :), intent(in) :: egv
    !
    integer :: i, j
    real(r8) :: f0, f1
    complex(r8), dimension(3) :: cv0, cv1

    f1 = 0.0_r8
    do i = 1, uc%na
        cv0 = egv((i - 1)*3 + 1:(i*3), 1)
        cv1 = egv((i - 1)*3 + 1:(i*3), 2)
        f0 = 0.0_r8
        do j = 1, 3
            f0 = f0 + abs(conjg(cv0(j))*cv1(j))
        end do
        f0 = f0**2
        f1 = f1 + f0*uc%isotope(i)%disorderparameter
    end do
    isotope_scattering_strength = f1
end function
