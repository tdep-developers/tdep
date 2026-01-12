module bubble
use konstanter, only: r8,lo_iou,lo_frequency_THz_to_Hartree,lo_kb_Hartree,lo_freqtol,lo_pi,lo_kappa_au_to_SI,lo_imag
use gottochblandat, only: lo_linspace,walltime,tochar,lo_array_of_array,qsort,lo_trapezoid_integration,lo_outerproduct
use type_qpointmesh, only: lo_qpoint_mesh
use lo_distributed_phonon_dispersion_relations, only: lo_distributed_phonon_dispersions
use type_crystalstructure, only: lo_crystalstructure
use lo_selfenergy_interpolation, only: lo_interpolated_selfenergy_grid
use lo_memtracker, only: lo_mem_helper
use mpi_wrappers, only: lo_mpi_helper
use lo_spectralfunction_helpers, only: lo_tapering_function,lo_integrate_spectral_function_and_return_nodes,&
    lo_find_spectral_function_max_and_fwhm,lo_interpolated_spectral_function,lo_evaluate_spectral_function,lo_gaussian_smear_spectral_function
use lo_accumulate_thermal_transport, only: lo_thermal_conductivity
use quadratures_stencils, only: lo_gaussianquadrature
implicit none


private
public :: bubble_only_transport

contains

subroutine bubble_only_transport(qp,dr,p,ise,smearing_prefactor,temperature,mw,mem,verbosity)
    !> q-point grid
    class(lo_qpoint_mesh), intent(in) :: qp
    !> harmonic (distributed) dispersions
    type(lo_distributed_phonon_dispersions), intent(in) :: dr
    !> structure
    type(lo_crystalstructure), intent(in) :: p
    !> self-energy interpolation
    type(lo_interpolated_selfenergy_grid), intent(inout) :: ise
    !> smearing prefactor
    real(r8), intent(in) :: smearing_prefactor
    !> temperature
    real(r8), intent(in) :: temperature
    !> mpi helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem
    !> talk a lot?
    integer, intent(in) :: verbosity

    type(lo_thermal_conductivity) :: tc
    type(lo_array_of_array), dimension(:), allocatable :: integration_node
    real(r8), parameter :: integraltol=1E-13_r8
    integer, parameter :: n_omegaT=20
    real(r8), dimension(:, :, :, :), allocatable :: buf_velsq
    real(r8), dimension(:, :), allocatable :: buf_sigmaIm, buf_sigmaRe, buf_smeared_Ja, buf_smeared_Jb
    real(r8), dimension(:), allocatable :: normalizationfactor
    real(r8), dimension(:), allocatable :: probing_omegaT
    integer :: ipt,iq
    real(r8) :: timer, t0, t1

    init: block
        real(r8) :: f0
        integer :: i
        ! Start timers
        timer = walltime()
        t0 = timer
        t1 = timer

        call mem%allocate(buf_sigmaIm,    [ise%n_energy, dr%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(buf_sigmaRe,    [ise%n_energy, dr%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(buf_smeared_Ja, [ise%n_energy, dr%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(buf_smeared_Jb, [ise%n_energy, dr%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(normalizationfactor, dr%n_mode, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(probing_omegaT, n_omegaT, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(buf_velsq, [3,3,dr%n_mode,dr%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        buf_sigmaIm=0.0_r8
        buf_sigmaRe=0.0_r8
        buf_smeared_Ja=0.0_r8
        buf_smeared_Jb=0.0_r8
        normalizationfactor=0.0_r8
        probing_omegaT=0.0_r8
        buf_velsq=0.0_r8

        ! Space for integration nodes
        allocate(integration_node(dr%n_mode))

        probing_omegaT=0.0_r8
        f0=0.1*lo_frequency_THz_to_Hartree
        do i=n_omegaT,2,-1
            probing_omegaT(i)=f0
            f0=f0*0.5_r8
        enddo

        ! Decide on a range of frequencies somehow
        !call lo_linspace(0.0_r8,0.002_r8*lo_frequency_THz_to_Hartree,probing_omegaT)

        ! Initialize a container for thermal transport information
        call tc%initialize(dr,ise%omega,probing_omegaT,temperature,mw)
    end block init

    qploop: do ipt=1,dr%n_irr_qpoint_local

        if (mw%talk) then
            t1 = walltime()
            write (lo_iou, *) '... q-point '//tochar(ipt)//' out of '//tochar(dr%n_irr_qpoint_local)
            t0 = t1
        end if

        fetchselfenergy: block
            real(r8), dimension(:), allocatable :: buf_taper, xmid, xlo, xhi, buf_dummyim
            real(r8) :: f0,f1,f2,omega,sigma
            integer :: imode

            call mem%allocate(xmid, dr%n_mode, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%allocate(xlo, dr%n_mode, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%allocate(xhi, dr%n_mode, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%allocate(buf_taper, ise%n_energy, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%allocate(buf_dummyim, ise%n_energy, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            xmid=0.0_r8
            xlo=0.0_r8
            xhi=0.0_r8
            buf_taper=0.0_r8
            buf_dummyim=0.0_r8
            call lo_tapering_function(ise%omega, buf_taper)

            buf_sigmaIm=0.0_r8
            buf_sigmaRe=0.0_r8
            buf_smeared_Ja=0.0_r8
            buf_smeared_Jb=0.0_r8
            ! interpolate self-energy to this q
            iq=dr%iq(ipt)%global_irreducible_index
            call ise%evaluate(p,&
                qp%ip(iq)%r,&
                dr%iq(ipt)%omega,&
                dr%iq(ipt)%egv,&
                buf_sigmaRe,buf_sigmaIm,mem)

            ! Figure out how to integrate over spectral functions, it can be a little wonky.
            do imode=1,dr%n_mode
                omega=dr%iq(ipt)%omega(imode)
                if ( omega .lt. lo_freqtol ) cycle
                ! Make sure we taper appropriately
                buf_sigmaIm(:, imode) = buf_sigmaIm(:, imode)*buf_taper
                ! Figure out extrema for integration
                call lo_find_spectral_function_max_and_fwhm(&
                    omega, &
                    ise%omega, &
                    buf_sigmaIm(:,imode), &
                    buf_sigmaRe(:,imode), &
                    xmid(imode),&
                    xlo(imode),&
                    xhi(imode))
                ! Get the integration nodes
                call integration_node(imode)%destroy()
                call lo_integrate_spectral_function_and_return_nodes(&
                    ise%omega, &
                    omega, &
                    buf_sigmaIm(:,imode), &
                    buf_sigmaRe(:,imode), &
                    xmid(imode),xlo(imode),xhi(imode),1.0_r8,temperature, integraltol, f0, f1, f2, &
                    integration_node(imode)%buf_1d_r8)
                normalizationfactor(imode)=1.0_r8/f0

                ! We likely want the smeared spectral functions as well. These are
                ! only used for plots, not for evaluating anything serious. I try to
                ! make sure the imaginary part of the self-energy is not so small that
                ! the peak might disappear between two frequency grid points.
                f0=ise%omega(2)-ise%omega(1)
                buf_dummyim=max(buf_sigmaIm(:,imode),f0)
                sigma=dr%iq(ipt)%sigma(imode)
                call lo_evaluate_spectral_function(ise%omega, buf_dummyim, buf_sigmaRe(:, imode), omega, buf_smeared_Ja(:, imode))
                call lo_gaussian_smear_spectral_function(ise%omega, sigma, buf_smeared_Ja(:, imode))
                f0 = lo_trapezoid_integration(ise%omega, buf_smeared_Ja(:, imode))
                buf_smeared_Ja(:,imode)=buf_smeared_Ja(:,imode)/f0
            enddo

            call mem%deallocate(xmid,                persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%deallocate(xlo,                 persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%deallocate(xhi,                 persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%deallocate(buf_taper,           persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%deallocate(buf_dummyim,         persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        end block fetchselfenergy

        prep_velsquare: block
            real(r8), dimension(3) :: v0,v1
            integer :: imode,jmode
            integer :: k,iop

            buf_velsq=0.0_r8
            do imode = 1, dr%n_mode
            do jmode = 1, dr%n_mode
                v0 = dr%iq(ipt)%genvel(:, imode, jmode)
                do k = 1, qp%ip(iq)%n_full_point
                    iop = qp%ip(iq)%operation_full_point(k)
                    v1 = matmul(p%sym%op(abs(iop))%m, v0)
                    buf_velsq(:, :, imode, jmode) = buf_velsq(:, :, imode, jmode) + lo_outerproduct(v1, v1)
                end do
            end do
            end do
            buf_velsq = buf_velsq/real(qp%ip(iq)%n_full_point, r8)
        end block prep_velsquare

        frequencyintegral: block
            real(r8) :: omegaA,omegaB,f0,f1,kappapref
            integer :: imode,jmode,i

            kappapref=0.5_r8*lo_pi*qp%ip(iq)%integration_weight/p%volume/temperature

            do imode=1,dr%n_mode
            do jmode=1,dr%n_mode
                omegaA=dr%iq(ipt)%omega(imode)
                omegaB=dr%iq(ipt)%omega(jmode)
                if ( omegaA .lt. lo_freqtol ) cycle
                if ( omegaB .lt. lo_freqtol ) cycle

                do i=1,n_omegaT
                    call integrate_over_two_spectral(&
                        integration_node(imode)%buf_1d_r8,&
                        integration_node(jmode)%buf_1d_r8,&
                        probing_omegaT(i),&
                        temperature,&
                        omegaA,&
                        omegaB,&
                        ise%omega,&
                        buf_sigmaIm(:,imode),&
                        buf_sigmaIm(:,jmode),&
                        buf_sigmaRe(:,imode),&
                        buf_sigmaRe(:,jmode),&
                        normalizationfactor(imode),&
                        normalizationfactor(jmode),&
                        f0,f1)
                    ! Accumulate
                    !tc%AC_kappa(:,:,i)=tc%AC_kappa(:,:,i) + f1*buf_velsq(:,:,imode,jmode)*kappapref
                    if ( imode .eq. jmode ) then
                        tc%AC_kappa(:,:,i)=tc%AC_kappa(:,:,i) + f1*buf_velsq(:,:,imode,jmode)*kappapref
                    else
                        tc%AC_kappa_od(:,:,i)=tc%AC_kappa_od(:,:,i) + f1*buf_velsq(:,:,imode,jmode)*kappapref
                    endif
                enddo
            enddo
            enddo
        end block frequencyintegral
    enddo qploop

    syncstuff: block
        real(r8) :: f0,f1,f2
        integer :: i

        call mw%allreduce('sum',tc%AC_kappa)

        if ( mw%talk ) then
            do i=1,n_omegaT
                f0=probing_omegaT(i)/lo_frequency_THz_to_Hartree
                f1=tc%AC_kappa(1,1,i)*lo_kappa_au_to_SI
                f2=tc%AC_kappa_od(1,1,i)*lo_kappa_au_to_SI
                write(*,*) i,f0,f1,f2
            enddo
        endif
    end block syncstuff

    ! Cleanup
    call mem%deallocate(buf_sigmaIm,         persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(buf_sigmaRe,         persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(buf_smeared_Ja,      persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(buf_smeared_Jb,      persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)

    call mem%deallocate(normalizationfactor, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(probing_omegaT,      persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(buf_velsq,           persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)

    call mem%assertzero(__FILE__,__LINE__)
end subroutine

subroutine integrate_over_two_spectral(nodeA,nodeB,omegaT,temperature,omegaA,omegaB,Omega,imA,imB,reA,reB,scaleA,scaleB,integrated_spectral,integrated_bubble)
    !> integration nodes for spectral function A
    real(r8), dimension(:), intent(in) :: nodeA
    !> integration nodes for spectral function B
    real(r8), dimension(:), intent(in) :: nodeB
    !> probing frequency
    real(r8), intent(in) :: omegaT
    !> actual temperature
    real(r8), intent(in) :: temperature
    !> harmonic frequencies
    real(r8), intent(in) :: omegaA,omegaB
    !> frequency axis
    real(r8), dimension(:), intent(in) :: Omega
    !> imaginary self energy
    real(r8), dimension(:), intent(in) :: imA,imB
    !> real self energy
    real(r8), dimension(:), intent(in) :: reA,reB
    !> normalization factors
    real(r8), intent(in) :: scaleA,scaleB
    !> integral of spectral function
    real(r8) :: integrated_spectral
    !> integral of bubble contribution
    real(r8) :: integrated_bubble

    real(r8), parameter :: small_omegaT_tol=1E-2_r8*lo_frequency_THz_to_Hartree
    real(r8), parameter :: mergetol=1E-5_r8*lo_frequency_THz_to_Hartree
    integer, parameter :: quadorder=8
    real(r8), dimension(2, quadorder) :: quad
    real(r8), dimension(:), allocatable :: nodes,weight
    real(r8) :: kappa_00_taylor
    real(r8) :: kappa_00_full

    ! Figure out the integrate nodes/weights
    mergenodes: block
        real(r8), dimension(:), allocatable :: dr0,dr1
        integer :: n0
        integer :: i,l,ii,jj

        n0=( size(nodeA)+size(nodeB) )*2

        ! Merge nodes, and duplicate them positive/negative
        allocate(dr0(n0))
        allocate(dr1(n0))
        dr0=0.0_r8
        dr1=0.0_r8
        l=0
        do i=1,size(nodeA)
            l=l+1
            dr0(l)=nodeA(i)
            l=l+1
            dr0(l)=-nodeA(i)
        enddo
        do i=1,size(nodeB)
            l=l+1
            dr0(l)=(nodeB(i)-omegaT)
            l=l+1
            dr0(l)=-(nodeB(i)-omegaT)
        enddo

        ! Merge all the integration nodes into one array, removing duplicates
        call qsort(dr0)

        l=1
        dr1(l)=dr0(1)
        do i=2,size(dr0)
            if ( abs(dr0(i)-dr1(l)) .gt. mergetol ) then
                l=l+1
                dr1(l)=dr0(i)
            endif
        enddo

        ! Length of integration array:
        i=(l-1)*quadorder
        ! Populate Gaussian quadrature
        allocate(nodes(i))
        allocate(weight(i))
        nodes=0.0_r8
        weight=0.0_r8
        do i=1,l-1
            ii=(i-1)*quadorder+1
            jj=(i)*quadorder
            call lo_gaussianquadrature(quadorder,dr1(i),dr1(i + 1), quad)
            nodes(ii:jj)=quad(1,:)
            weight(ii:jj)=quad(2,:)
        enddo
    end block mergenodes

    ! Now for the actual integral
    actualintegral: block
        real(r8) :: nA,nB,jA,jB,x
        real(r8) :: n_pref,maxomega
        real(r8) :: f0,f1,f2
        integer :: i

        n_pref=1.0_r8/(lo_kb_Hartree*temperature)
        maxomega=Omega(size(Omega)-1)

        kappa_00_full=0.0_r8
        kappa_00_taylor=0.0_r8

        f0=0.0_r8
        f1=0.0_r8
        do i=1,size(nodes)
            ! evaluate the two spectral functions. I only have them defined for positive
            ! but I know they are odd so it's simple enough.
            x=abs(nodes(i))
            if ( x < maxomega ) then
                jA=lo_interpolated_spectral_function(Omega, imA,reA,omegaA,scaleA, x)
            else
                jA=0.0_r8
            endif
            if ( nodes(i) .lt. 0.0_r8 ) then
                jA=-jA
            endif

            x=abs(nodes(i)-omegaT)
            if ( x < maxomega ) then
                jB=lo_interpolated_spectral_function(Omega, imB,reB,omegaB,scaleB, x)
            else
                jB=0.0_r8
            endif
            if ( nodes(i)-omegaT .lt. 0.0_r8 ) then
                jB=-jB
            endif

            ! evaluate bose-einstein factors
            x=nodes(i)
            if ( abs(x) .lt. lo_freqtol ) then
                nA=0.0_r8
            elseif ( x*n_pref .gt. 50.0_r8 ) then
                nA=0.0_r8
            elseif ( x*n_pref .lt. -50.0_r8 ) then
                nA=-1.0_r8
            else
                nA=1.0_r8/(exp(x*n_pref)-1.0_r8)
            endif

            x=(nodes(i)-omegaT)
            if ( abs(x) .lt. lo_freqtol ) then
                nB=0.0_r8
            elseif ( x*n_pref .gt. 50.0_r8 ) then
                nB=0.0_r8
            elseif ( x*n_pref .lt. -50.0_r8 ) then
                nB=-1.0_r8
            else
                nB=1.0_r8/(exp(x*n_pref)-1.0_r8)
            endif

            ! evaluate scalar conductivity vertices


            ! Accumulate integral. First is just spectral function.
            f0=f0 + abs(jA)*weight(i)

            ! The the bubble contribution
            if ( omegaT .lt. small_omegaT_tol ) then
                ! We use the Taylor expansion
                f2=(nA*(1.0_r8 + nA))*n_pref
                f2=f2 + 0.5_r8*(nA*(1.0_r8 + nA)*(1.0_r8 + 2*nA))*n_pref*n_pref*omegaT
                f2=f2*jA*jB*nodes(i)*nodes(i)
                f1=f1 + f2*weight(i)
            else
                ! Full expression
                f2=(nB-nA)/omegaT
                f2=f2*jA*jB*nodes(i)*nodes(i)
                f1=f1+f2*weight(i)
            endif
        enddo

        ! kappa_00_taylor
        ! kappa_00_full
        ! kappa_ij_taylor
        ! kappa_ij_full

        integrated_spectral=f0
        integrated_bubble=f1

    end block actualintegral
end subroutine

end module