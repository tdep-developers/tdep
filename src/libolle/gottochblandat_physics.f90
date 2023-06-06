#include "precompilerdefinitions"
submodule (gottochblandat) gottochblandat_physics
implicit none
contains

!> Calculate the reciprocal basis
#ifdef AGRESSIVE_SANITY
module function lo_reciprocal_basis(a) result(b)
#else
pure module function lo_reciprocal_basis(a) result(b)
#endif
    !> the basis
    real(flyt), dimension(3,3), intent(in) :: a
    !> reciprocal basis
    real(flyt), dimension(3,3) :: b

    real(flyt), dimension(3) :: a1,a2,a3,b1,b2,b3
    real(flyt) :: vol
    a1=a(:,1)
    a2=a(:,2)
    a3=a(:,3)
    vol=abs(lo_determ(a))
#ifdef AGRESSIVE_SANITY
    if ( vol .lt. lo_sqtol ) then
        call lo_stop_gracefully(['Trying to calculate reciprocal basis of a cell with zero volume'],&
                                lo_exitcode_param,__FILE__,__LINE__)
    endif
#endif
    b1=lo_cross(a2,a3)
    b2=lo_cross(a3,a1)
    b3=lo_cross(a1,a2)
    b(:,1)=b1/vol
    b(:,2)=b2/vol
    b(:,3)=b3/vol
end function

!> Return a sensible-ish k-mesh density for a given cell from a linear density
pure module function lo_kmesh_density(reciprocal_latticevectors,na,nq) result(density)
    !> reciprocal latticevectors
    real(flyt), dimension(3,3), intent(in) :: reciprocal_latticevectors
    !> number of atoms in the cell
    integer, intent(in) :: na
    !> linear density
    integer, intent(in) :: nq
    !> resulting q-density
    integer, dimension(3) :: density

    real(flyt) :: f0,a,b,c,avg

    ! So, what do I want. I want that product(density)*natoms = constant
    ! So this is approximately the total number of points.
    f0=real(nq**3,flyt)/real(na,flyt)
    ! And this would be the sensible number per direction
    f0=f0**(1.0_flyt/3.0_flyt)

    ! Then, It should be proportional to the length of the reciprocal lattice vectors
    a=norm2(reciprocal_latticevectors(:,1))
    b=norm2(reciprocal_latticevectors(:,2))
    c=norm2(reciprocal_latticevectors(:,3))
    avg=(a+b+c)/3.0_flyt
    ! This should be somewhat ok, I guess
    density(1)=ceiling(f0*a/avg)
    density(2)=ceiling(f0*b/avg)
    density(3)=ceiling(f0*c/avg)
end function

!> Calculate a,b,c,alpha,beta,gamma from a set of lattice vectors
module pure subroutine lo_get_axis_angles(basis,a,b,c,al,be,gm)
    !> the basis
    real(flyt), dimension(3,3), intent(in) :: basis
    !> lattice parameter a
    real(flyt), intent(out) :: a
    !> lattice parameter a
    real(flyt), intent(out) :: b
    !> lattice parameter a
    real(flyt), intent(out) :: c
    !> angle alpha
    real(flyt), intent(out) :: al
    !> angle beta
    real(flyt), intent(out) :: be
    !> angle gamma
    real(flyt), intent(out) :: gm
    !
    real(flyt), dimension(3) :: a1,a2,a3
    !
    a1=basis(:,1)
    a2=basis(:,2)
    a3=basis(:,3)
    a=norm2(a1)
    b=norm2(a2)
    c=norm2(a3)
    al=acos( dot_product( a2,a3 )/(b*c) )
    be=acos( dot_product( a1,a3 )/(a*c) )
    gm=acos( dot_product( a1,a2 )/(a*b) )
end subroutine

!> Fermi-Dirac smearing. The derivative of the Fermi-Dirac occupation function: $$ \frac{\partial f}{\partial \epsilon} $$. Can perhaps be used in place of a Gaussian or Lorentzian for smearing.
module elemental function lo_fermidirac(x,mu,sigma) result(f)
    !> point to evaluate
    real(flyt), intent(in) :: x
    !> mean
    real(flyt), intent(in) :: mu
    !> width
    real(flyt), intent(in) :: sigma
    !> the intensity
    real(flyt) :: f

    real(flyt) :: f0,f1

    ! Avoid overflows
    if ( abs(x-mu) .gt. 20*sigma ) then
        f=0.0_flyt
        return
    else
        f0=exp( (x-mu)/sigma )
        f1=sigma*((1.0_flyt+f0)**2)
        f=f0/f1
    endif
end function lo_fermidirac

!> Fermi-Dirac occupation function $$ f=\frac{1}{\exp\left[(\epsilon-\epsilon_F)/k_BT \right]+1 } $$ Temperature tolerance it is replaced with a step function, same thing when the energy is very far from the Fermi level.
module elemental function lo_fermi(energy,efermi,temperature) result(f)
    !> energy
    real(flyt), intent(in) :: energy
    !> Fermi level
    real(flyt), intent(in) :: efermi
    !> Temperature
    real(flyt), intent(in) :: temperature
    !> occupation
    real(flyt) :: f

    ! defaults to 1E-3 K
    if ( temperature .lt. lo_temperaturetol ) then
        ! it's a step function. Technically it should be a Heaviside step function, such that if
        ! energy = efermi the value is 1/2. This gets a little weird with floating point precision.
        ! Probably better to handle that with a tiny temperature. Something future Olle will figure
        ! out when something breaks.
        if ( energy .gt. efermi ) then
            f=0.0_flyt
        else
            f=1.0_flyt
        endif
    else
        ! 36 sigma makes makes the error from using a step function disappear in the finite precision.
        if ( abs(energy-efermi) .lt. 36*lo_kb_hartree*temperature ) then
            ! use the proper distribution
            f=1.0_flyt/(exp( (energy-efermi)/(lo_kb_hartree*temperature) )+1.0_flyt)
        else
            ! replace with step function, avoids overflows.
            if ( energy .lt. efermi ) then
                f=1.0_flyt
            else
                f=0.0_flyt
            endif
        endif
    endif
end function lo_fermi

!> Planck distribution function. $$ f=\frac{1}{\exp\left( \hbar\omega/k_BT \right)-1 } $$ Assumes frequencies in Hz and temperature in K. Below a frequency or temperature tolerance it returns zero for numerical stability.
module elemental function lo_planck(temperature,omega) result(n)
    !> Temperature in K
    real(flyt), intent(in) :: temperature
    !> Frequency in Hartree
    real(flyt), intent(in) :: omega
    !> the occupation
    real(flyt) :: n

    real(flyt) :: x
    ! Get the T=0 limit correct
    if ( temperature .lt. lo_temperaturetol ) then
        n=0.0_flyt
        return
    endif
    ! And the omega = 0 limit
    if ( omega .lt. lo_freqtol ) then
        n=0.0_flyt
        return
    endif

    x=omega/(lo_kb_Hartree*temperature)
    if ( x .gt. 1E2_flyt ) then
        ! corresponds to 1E-44 or something like that.
        n=0.0_flyt
    else
        ! Use the real distribution
        n=1.0_flyt/(exp(x)-1.0_flyt)
    endif
end function lo_planck

!> Temperature derivative of the Planck distribution, in 1/K
module elemental function lo_planck_deriv(T,omega) result(dndt)
    !> Temperature in K
    real(flyt), intent(in) :: T
    !> Frequency in Hartree
    real(flyt), intent(in) :: omega
    !> The derivative
    real(flyt) :: dndt
    !
    real(flyt) :: x,n
    if ( T .lt. lo_temperaturetol ) then
        dndt=0.0_flyt
        return
    endif
    if ( omega .lt. lo_freqtol ) then
        dndt=0.0_flyt
        return
    endif
    ! ok, finite derivative
    x=(omega/lo_kb_hartree/T)
    n=lo_planck(T,omega)
    dndt=n*n*exp(x)*x/T
end function lo_planck_deriv

!> Second temperature derivative of the Planck distribution, in 1/K^2
module elemental function lo_planck_secondderiv(T,omega) result(ddnddt)
    !> Temperature in K
    real(flyt), intent(in) :: T
    !> Frequency in Hartree
    real(flyt), intent(in) :: omega
    !> Second derivative
    real(flyt) :: ddnddt
    !
    real(flyt) :: x,ex,n
    if ( T .lt. lo_temperaturetol ) then
        ddnddt=0.0_flyt
        return
    endif
    if ( omega .lt. lo_freqtol ) then
        ddnddt=0.0_flyt
        return
    endif
    x=(omega/lo_kb_hartree/T)
    ex=exp(x)
    n=lo_planck(T,omega)
    ddnddt=(ex*x*(2 + ex*(-2 + x) + x))/((-1 + ex)**3*T**2)
end function lo_planck_secondderiv

!> The free energy for a single quantum harmonic oscillator in eV $$ F = \frac{\hbar\omega}{2} - k_b T \log(n+1) $$ Assumes temperatures in K and frequency in Hz.
module elemental function lo_harmonic_oscillator_free_energy(temp,omega) result(f)
    !> Temperature in K
    real(flyt), intent(in) :: temp
    !> Angular frequency in Hartree
    real(flyt), intent(in) :: omega
    !> the free energy in Hartree
    real(flyt) :: f

    real(flyt) :: n
    ! return 0 for omega <= 0. Should be undefined, but that should be taken care of elsewhere.
    if ( omega .lt. lo_freqtol ) then
        f=0.0_flyt
        return
    endif
    ! correct T=0 limit without numerical instability
    if ( temp .lt. lo_temperaturetol ) then
        f=omega*0.5_flyt
        return
    endif
    ! all normal temperatures/frequencies
    n=lo_planck(temp,omega)
    f=omega*0.5_flyt-lo_kb_hartree*temp*log(1.0_flyt+n)
end function

!> Harmonic oscillator internal energy
module elemental function lo_harmonic_oscillator_internal_energy(temp,omega) result(u)
    !> Temperature in K
    real(flyt), intent(in) :: temp
    !> Angular frequency in Hartree
    real(flyt), intent(in) :: omega
    !> Internal energy in Hartree
    real(flyt) :: u

    u=omega*0.5_flyt+lo_planck(temp,omega)*omega
end function

!> Free energy of classical harmonic oscillator
module elemental function lo_classical_harmonic_oscillator_free_energy(temp,omega) result(f)
    !> Temperature in K
    real(flyt), intent(in) :: temp
    !> Angular frequency in Hartree
    real(flyt), intent(in) :: omega
    !> the free energy in Hartree
    real(flyt) :: f

    real(flyt) :: x
    ! undefined for tiny omega
    if ( omega .lt. lo_freqtol ) then
        f=0.0_flyt
        return
    endif
    ! nothing for tiny temperatures
    if ( temp .lt. lo_temperaturetol ) then
        f=0.0_flyt
        return
    endif
    ! normal T and omega
    x=(omega/lo_kb_hartree/temp)
    f=lo_kb_hartree*temp*log(x)
end function

!> The entropy for a single harmonic oscillator in eV/K
module elemental function lo_harmonic_oscillator_entropy(temp,omega) result(s)
    !> Temperature in K
    real(flyt), intent(in) :: temp
    !> Angular frequency in hartree
    real(flyt), intent(in) :: omega
    !> The entropy in Hartreee/K
    real(flyt) :: s
    !
    real(flyt) :: n

    if ( omega .lt. lo_freqtol ) then
        s=0.0_flyt
        return
    endif

    if ( temp .lt. lo_temperaturetol ) then
        s=0.0_flyt
        return
    endif

    n=lo_planck(temp,omega)
    ! make sure we don't take log zero, use the Taylor expansion below that
    if ( n .gt. 1E-8_flyt ) then
        s=lo_kb_hartree*( (1.0_flyt+n)*log(1.0_flyt+n)-n*log(n) )
    else
        s=n*lo_kb_hartree
    endif
end function

!> Heat capacity (Cv) for a single harmonic oscillator in eV/K
module elemental function lo_harmonic_oscillator_cv(temp,omega) result(cv)
    !> Temperature in K
    real(flyt), intent(in) :: temp
    !> Angular frequency in Hartree
    real(flyt), intent(in) :: omega
    !> The heat capacity in Hartree/K
    real(flyt) :: cv

    real(flyt) :: x

    if ( omega .lt. lo_freqtol ) then
        cv=0.0_flyt
        return
    endif

    if ( temp .lt. lo_temperaturetol ) then
        cv=0.0_flyt
        return
    endif

    x=omega/temp/lo_kb_hartree
    ! have to make sure we avoid overflows. Replace with high/low temperature Taylor expansions
    ! where applicable. This should still be ok to machine precision.
    if ( x .gt. 200.0_flyt ) then
        ! very low temperature
        cv=0.0_flyt
    elseif ( x .gt. 21.0_flyt ) then
        ! pretty low temperature
        cv=lo_kb_hartree*x*x*exp(-x)
    elseif ( x .lt. 2E-2_flyt ) then
        ! really high temperature
        cv=lo_kb_hartree*(1.0_flyt-x*x/12_flyt+x*x*x*x/240.0_flyt)
    else
        ! in the middle
        cv=lo_kb_hartree*x*x*exp(x)/(exp(x)-1)**2
    endif
end function

!> Untangle the bands the best I can for a single tetrahedron.
module subroutine lo_untangle_one_tetrahedron(corner,energy,groupvelocity,degentol,permutation,success,error_interp_energy,error_interp_gradient)
    !> corners of tetrahedron
    real(r8), dimension(3,4), intent(in) :: corner
    !> energies at corners
    real(r8), dimension(:,:), intent(in) :: energy
    !> group velocities at corners
    real(r8), dimension(:,:,:), intent(in) :: groupvelocity
    !> tolerance for determining degeneracy
    real(r8), intent(in) :: degentol
    !> resulting permutation of bands
    integer, dimension(:,:), intent(out) :: permutation
    !> exit status, 0 if all fine, -1 if all bands degenerate
    integer, intent(out) :: success
    !> perhaps return a measure of the interpolation error in energy?
    real(r8), dimension(:), intent(out), optional :: error_interp_energy
    !> perhaps return the error in interpolated group velocities?
    real(r8), dimension(:), intent(out), optional :: error_interp_gradient

    integer :: nband,refcorner
    logical :: checkerr

    ! First decide on some simple things
    init: block
        integer :: i,j,icrn

        ! Number of bands we are working with
        nband=size(energy,1)

        ! Should I measure the error?
        checkerr=.false.
        if ( present(error_interp_energy) ) checkerr=.true.
        if ( present(error_interp_gradient) ) checkerr=.true.
        ! Either we measure both errors or none.
        if ( present(error_interp_energy) .neqv. present(error_interp_gradient) ) then
            !call lo_stop_gracefully(['Either measure both errors or none, no half-assing.'],lo_exitcode_param,__FILE__,__LINE__)
        endif

        ! First we need to find a corner that is not
        ! degenerate to start untangling from:
        refcorner=0
        do icrn=1,4
            j=0
            do i=1,nband-1
                if ( abs(energy(i,icrn)-energy(i+1,icrn)) .lt. degentol ) j=j+1
            enddo
            if ( j .eq. 0 ) then
                ! yup, found a corner that is not degenerate
                refcorner=icrn
                exit
            endif
        enddo

        ! Not clear what to do if all corners are degenerate.
        ! The sensible thing to do is to inject a new point in the
        ! middle of the tetrahedron to sort things out. Anyway, that
        ! is not something that should be handled in this routine,
        ! just return a flag explaining what happened.
        if ( refcorner .eq. 0 ) then
            success=-1
            ! But still return a valid permutation. Just pick the normal one.
            do i=1,nband
                permutation(i,:)=i
            enddo
            ! If we want the errors, not sure what to do. Set the error to a huge number.
            if ( checkerr ) then
                error_interp_energy=1000.0_r8
                error_interp_gradient=1000.0_r8
            endif
            ! And return early.
            return
        else
            ! This went fine! Go ahead like we planned!
            success=0
            ! Set the permutation to nothing
            permutation=0
            ! And set the reference permutation for the reference corner
            do i=1,nband
                permutation(i,refcorner)=i
            enddo
        endif
    end block init

    ! With the basics taken care of, untangle it the best we can.
    untangle: block
        real(flyt), dimension(3) :: v0
        real(flyt) :: e0
        real(r8) :: f0,f1
        integer :: icrn,iband,jband,k

        cornerloop: do icrn=1,4
            ! the reference corner we already know is ok.
            if ( icrn .eq. refcorner ) cycle cornerloop
            ! Vector from the reference corner to this corner
            v0=corner(:,icrn)-corner(:,refcorner)
            ! Now untangle each band as best we can.
            do iband=1,nband
                ! Energy at reference corner
                e0=energy(iband,refcorner)
                ! Assuming a linear approximation based on the group velocity
                ! I guess what the energy is supposed to be at the next corner
                e0=e0+dot_product(groupvelocity(:,iband,refcorner),v0)

                ! Find the closest match at this corner
                f0=lo_huge
                k=0
                do jband=1,nband
                    f1=abs( energy(jband,icrn)-e0 )
                    if ( f1 .lt. f0 ) then
                        k=jband
                        f0=f1
                    endif
                enddo
                ! And we know the permutation, sort of.
                permutation(iband,icrn)=k
            enddo
        enddo cornerloop
    end block untangle

    if ( checkerr ) then
    testinterpolation: block
        real(r8), dimension(3,4) :: dvel
        real(r8), dimension(4) :: tete
        real(r8), dimension(3) :: v0,v1
        real(r8) :: e0,e1,e2
        integer :: icrn,iband

        ! Measure error in energy you get by using the group velocity at
        ! the reference tetrahedron as a linear interpolation.
        error_interp_energy=0.0_r8
        do icrn=1,4
            ! the reference corner we already know is ok.
            if ( icrn .eq. refcorner ) cycle
            ! Vector from the reference corner to this corner
            v0=corner(:,icrn)-corner(:,refcorner)
            e2=0.0_r8
            do iband=1,nband
                e0=energy(iband,refcorner)
                e0=e0+dot_product(groupvelocity(:,iband,refcorner),v0)
                e1=energy(iband,icrn)
                error_interp_energy(iband)=error_interp_energy(iband)+( (e1-e0)**2 )*0.25_r8
            enddo
        enddo
        error_interp_energy=sqrt(error_interp_energy)

        error_interp_gradient=0.0_flyt
        ! Measure the absolute error in group velocity. Or relative? No absolute I think.
        do iband=1,nband
            v0=groupvelocity(:,iband,refcorner)
            do icrn=1,4
                tete(icrn)=energy(iband,icrn)
                dvel(:,icrn)=groupvelocity(:,iband,icrn)
            enddo
            v1=lo_linear_gradient_in_tetrahedron(corner,tete,1E-12_r8)
            do icrn=1,4
                error_interp_gradient(iband)=error_interp_gradient(iband)+lo_sqnorm(v1-dvel(:,icrn))*0.25_flyt
            enddo
        enddo
        error_interp_gradient=sqrt(error_interp_gradient)
    end block testinterpolation

    ! testinterpolation: block
    !     real(r8), dimension(3,4) :: dvel
    !     real(r8), dimension(4) :: tete
    !     real(r8), dimension(3) :: v0,v1
    !     real(r8) :: e0,e1,e2
    !     integer :: icrn,iband,jband
    !
    !     ! Measure error in energy you get by using the group velocity at
    !     ! the reference tetrahedron as a linear interpolation.
    !     error_interp_energy=0.0_r8
    !     do icrn=1,4
    !         ! the reference corner we already know is ok.
    !         if ( icrn .eq. refcorner ) cycle
    !         ! Vector from the reference corner to this corner
    !         v0=corner(:,icrn)-corner(:,refcorner)
    !         e2=0.0_r8
    !         do iband=1,nband
    !             jband=permutation(iband,icrn)
    !             e0=energy(iband,refcorner)
    !             e0=e0+dot_product(groupvelocity(:,iband,refcorner),v0)
    !             e1=energy(jband,icrn)
    !             error_interp_energy(iband)=error_interp_energy(iband)+( (e1-e0)**2 )*0.25_r8
    !         enddo
    !     enddo
    !     error_interp_energy=sqrt(error_interp_energy)
    !
    !     error_interp_gradient=0.0_flyt
    !     ! Measure the absolute error in group velocity. Or relative? No absolute I think.
    !     do iband=1,nband
    !         v0=groupvelocity(:,iband,refcorner)
    !         do icrn=1,4
    !             jband=permutation(iband,icrn)
    !             tete(icrn)=energy(jband,icrn)
    !             dvel(:,icrn)=groupvelocity(:,jband,icrn)
    !         enddo
    !         v1=lo_linear_gradient_in_tetrahedron(corner,tete,1E-12_r8)
    !         do icrn=1,4
    !             error_interp_gradient(iband)=error_interp_gradient(iband)+lo_sqnorm(v1-dvel(:,icrn))*0.25_flyt
    !         enddo
    !     enddo
    !     error_interp_gradient=sqrt(error_interp_gradient)
    !
    ! end block testinterpolation
    endif
end subroutine

!> Finds the rotation matrix such that lattice'=rotation*lattice represents a diagnoal distortion of ref_lattice.
module subroutine lo_find_rotation_that_makes_strain_diagonal( ref_lattice,lattice,rotation,strain,guess )
    !> reference lattice
    real(flyt), dimension(3,3), intent(in) :: ref_lattice
    !> current lattice
    real(flyt), dimension(3,3), intent(in) :: lattice
    !> rotation matrix
    real(flyt), dimension(3,3), intent(out) :: rotation
    !> strain matrix
    real(flyt), dimension(3,3), intent(out) :: strain
    !> guess for the symmetric strain?
    real(flyt), dimension(3,3), intent(in), optional :: guess

    integer, parameter :: nouter=1000
    integer, parameter :: ninner=1000
    integer, parameter :: maxctr=20
    real(flyt), parameter :: defstep=1E-5_flyt
    real(flyt), parameter :: errtol=1E-20_flyt

    real(flyt), dimension(6,6,6) :: coeffM
    real(flyt), dimension(6) :: bvec,eta0,eta1,grad
    real(flyt), dimension(3,3) :: m0
    real(flyt) :: e0,e1,f0,step
    integer :: initer,outiter,ctr,i

    ! So, what do I solve here. The idea is that any transformation of a lattice takes
    ! a,b,c,alpha,beta,gamma to new values. I need six numbers to define that transformation
    ! uniquely. But you often come across these transformations as 3x3-matrices, so that they
    ! might contain an extra rotation. This routine finds the diagonal form of a general
    ! transformation between lattice vectors, since that can always be done.
    !
    ! So, in practice, I write the transformed lattice vectors
    !
    ! a'=eta*a
    ! b'=eta*b
    ! c'=eta*c
    !
    ! And find the symmetric matrix eta that produces the correct values for
    !
    !  a.a b.b c.c a.b a.c b.c
    !
    ! when applied to the original lattice. I end up with a bunch of matrix
    ! equations on the form X^T A X = b, that I solve with a steepest descent.
    ! Seems to be a robust procedure.

    ! Get the coefficient matrices, A in XAX=b
    call cfvars(coeffM(:,:,1),coeffM(:,:,2),coeffM(:,:,3),coeffM(:,:,4),coeffM(:,:,5),coeffM(:,:,6),ref_lattice)
    ! Calculate the target quantities. All combinations of scalar products of
    ! the lattice vectors, b in XAX=b
    bvec(1)=dot_product(lattice(:,1),lattice(:,1))
    bvec(2)=dot_product(lattice(:,2),lattice(:,2))
    bvec(3)=dot_product(lattice(:,3),lattice(:,3))
    bvec(4)=dot_product(lattice(:,1),lattice(:,2))
    bvec(5)=dot_product(lattice(:,2),lattice(:,3))
    bvec(6)=dot_product(lattice(:,3),lattice(:,1))
    ! Do we have a guess for the strain?
    if ( present(guess) ) then
        eta0(1)=guess(1,1)
        eta0(2)=guess(2,2)
        eta0(3)=guess(3,3)
        eta0(4)=guess(1,2)
        eta0(5)=guess(1,3)
        eta0(6)=guess(2,3)
    else
        ! Start from something not completely crazy
        m0=matmul(ref_lattice,lo_invert3x3matrix(lattice))
        eta0(1)=m0(1,1)
        eta0(2)=m0(2,2)
        eta0(3)=m0(3,3)
        eta0(4)=m0(1,2)
        eta0(5)=m0(1,3)
        eta0(6)=m0(2,3)
    endif
    eta1=0.0_flyt
    grad=0.0_flyt

    outloop: do outiter=1,nouter
        ! Get the error and gradient
        e0=0.0_flyt
        grad=0.0_flyt
        do i=1,6
            f0=dot_product(eta0,matmul(coeffM(:,:,i),eta0))-bvec(i)
            e0=e0+f0**2
            grad=grad+matmul( coeffM(:,:,i)+transpose(coeffM(:,:,i)),eta0)*2*f0
        enddo

        ! Line search in this direction
        ctr=0
        step=defstep
        inloop: do initer=1,ninner
            ! Guess new point
            eta1=eta0-step*grad
            ! Get error at new point
            e1=0.0_flyt
            do i=1,6
                e1=e1+( dot_product(eta1,matmul(coeffM(:,:,i),eta1))-bvec(i) )**2
            enddo

            ! See how it went:
            if ( e1 .lt. errtol ) then
                ! We are done
                eta0=eta1
                e0=e1
                exit outloop
            elseif ( e1 .lt. e0 ) then
                ! Good step, keep going
                eta0=eta1
                e0=e1
                step=step*1.3_flyt
            else
                ! Bad step, decrease and try again
                step=step*0.1_flyt
                ctr=ctr+1
            endif

            ! If we have decreased too many times, start over with new gradient
            if ( ctr .ge. maxctr ) exit inloop
        enddo inloop
    enddo outloop

    ! Convert to a strain matrix
    strain(1,1)=eta0(1)
    strain(2,2)=eta0(2)
    strain(3,3)=eta0(3)
    strain(1,2)=eta0(4)
    strain(2,1)=eta0(4)
    strain(1,3)=eta0(5)
    strain(3,1)=eta0(5)
    strain(2,3)=eta0(6)
    strain(3,2)=eta0(6)
    ! Get the rotation
    m0=matmul(strain,ref_lattice)
    rotation=matmul(m0,lo_invert3x3matrix(lattice))
    rotation=rotation/lo_determ(rotation)**(1.0_flyt/3.0_flyt)

    contains

    ! These are magic matrices that define the quadratic form. Thanks Mathematica.
    subroutine cfvars(cm1,cm2,cm3,cm4,cm5,cm6,basis)
        real(flyt), dimension(6,6), intent(out) :: cm1,cm2,cm3,cm4,cm5,cm6
        real(flyt), dimension(3,3), intent(in) :: basis
        real(flyt) :: ax,ay,az,bx,by,bz,cx,cy,cz
        ax=basis(1,1)
        ay=basis(2,1)
        az=basis(3,1)
        bx=basis(1,2)
        by=basis(2,2)
        bz=basis(3,2)
        cx=basis(1,3)
        cy=basis(2,3)
        cz=basis(3,3)

        cm1(1,1)=ax**2
        cm1(1,2)=0
        cm1(1,3)=0
        cm1(1,4)=2*ax*ay
        cm1(1,5)=2*ax*az
        cm1(1,6)=0
        cm1(2,1)=0
        cm1(2,2)=ay**2
        cm1(2,3)=0
        cm1(2,4)=2*ax*ay
        cm1(2,5)=0
        cm1(2,6)=2*ay*az
        cm1(3,1)=0
        cm1(3,2)=0
        cm1(3,3)=az**2
        cm1(3,4)=0
        cm1(3,5)=2*ax*az
        cm1(3,6)=2*ay*az
        cm1(4,1)=0
        cm1(4,2)=0
        cm1(4,3)=0
        cm1(4,4)=ax**2 + ay**2
        cm1(4,5)=2*ay*az
        cm1(4,6)=2*ax*az
        cm1(5,1)=0
        cm1(5,2)=0
        cm1(5,3)=0
        cm1(5,4)=0
        cm1(5,5)=ax**2 + az**2
        cm1(5,6)=2*ax*ay
        cm1(6,1)=0
        cm1(6,2)=0
        cm1(6,3)=0
        cm1(6,4)=0
        cm1(6,5)=0
        cm1(6,6)=ay**2 + az**2

        cm2(1,1)=bx**2
        cm2(1,2)=0
        cm2(1,3)=0
        cm2(1,4)=2*bx*by
        cm2(1,5)=2*bx*bz
        cm2(1,6)=0
        cm2(2,1)=0
        cm2(2,2)=by**2
        cm2(2,3)=0
        cm2(2,4)=2*bx*by
        cm2(2,5)=0
        cm2(2,6)=2*by*bz
        cm2(3,1)=0
        cm2(3,2)=0
        cm2(3,3)=bz**2
        cm2(3,4)=0
        cm2(3,5)=2*bx*bz
        cm2(3,6)=2*by*bz
        cm2(4,1)=0
        cm2(4,2)=0
        cm2(4,3)=0
        cm2(4,4)=bx**2 + by**2
        cm2(4,5)=2*by*bz
        cm2(4,6)=2*bx*bz
        cm2(5,1)=0
        cm2(5,2)=0
        cm2(5,3)=0
        cm2(5,4)=0
        cm2(5,5)=bx**2 + bz**2
        cm2(5,6)=2*bx*by
        cm2(6,1)=0
        cm2(6,2)=0
        cm2(6,3)=0
        cm2(6,4)=0
        cm2(6,5)=0
        cm2(6,6)=by**2 + bz**2

        cm3(1,1)=cx**2
        cm3(1,2)=0
        cm3(1,3)=0
        cm3(1,4)=2*cx*cy
        cm3(1,5)=2*cx*cz
        cm3(1,6)=0
        cm3(2,1)=0
        cm3(2,2)=cy**2
        cm3(2,3)=0
        cm3(2,4)=2*cx*cy
        cm3(2,5)=0
        cm3(2,6)=2*cy*cz
        cm3(3,1)=0
        cm3(3,2)=0
        cm3(3,3)=cz**2
        cm3(3,4)=0
        cm3(3,5)=2*cx*cz
        cm3(3,6)=2*cy*cz
        cm3(4,1)=0
        cm3(4,2)=0
        cm3(4,3)=0
        cm3(4,4)=cx**2 + cy**2
        cm3(4,5)=2*cy*cz
        cm3(4,6)=2*cx*cz
        cm3(5,1)=0
        cm3(5,2)=0
        cm3(5,3)=0
        cm3(5,4)=0
        cm3(5,5)=cx**2 + cz**2
        cm3(5,6)=2*cx*cy
        cm3(6,1)=0
        cm3(6,2)=0
        cm3(6,3)=0
        cm3(6,4)=0
        cm3(6,5)=0
        cm3(6,6)=cy**2 + cz**2

        cm4(1,1)=ax*bx
        cm4(1,2)=0
        cm4(1,3)=0
        cm4(1,4)=ay*bx + ax*by
        cm4(1,5)=az*bx + ax*bz
        cm4(1,6)=0
        cm4(2,1)=0
        cm4(2,2)=ay*by
        cm4(2,3)=0
        cm4(2,4)=ay*bx + ax*by
        cm4(2,5)=0
        cm4(2,6)=az*by + ay*bz
        cm4(3,1)=0
        cm4(3,2)=0
        cm4(3,3)=az*bz
        cm4(3,4)=0
        cm4(3,5)=az*bx + ax*bz
        cm4(3,6)=az*by + ay*bz
        cm4(4,1)=0
        cm4(4,2)=0
        cm4(4,3)=0
        cm4(4,4)=ax*bx + ay*by
        cm4(4,5)=az*by + ay*bz
        cm4(4,6)=az*bx + ax*bz
        cm4(5,1)=0
        cm4(5,2)=0
        cm4(5,3)=0
        cm4(5,4)=0
        cm4(5,5)=ax*bx + az*bz
        cm4(5,6)=ay*bx + ax*by
        cm4(6,1)=0
        cm4(6,2)=0
        cm4(6,3)=0
        cm4(6,4)=0
        cm4(6,5)=0
        cm4(6,6)=ay*by + az*bz

        cm5(1,1)=bx*cx
        cm5(1,2)=0
        cm5(1,3)=0
        cm5(1,4)=by*cx + bx*cy
        cm5(1,5)=bz*cx + bx*cz
        cm5(1,6)=0
        cm5(2,1)=0
        cm5(2,2)=by*cy
        cm5(2,3)=0
        cm5(2,4)=by*cx + bx*cy
        cm5(2,5)=0
        cm5(2,6)=bz*cy + by*cz
        cm5(3,1)=0
        cm5(3,2)=0
        cm5(3,3)=bz*cz
        cm5(3,4)=0
        cm5(3,5)=bz*cx + bx*cz
        cm5(3,6)=bz*cy + by*cz
        cm5(4,1)=0
        cm5(4,2)=0
        cm5(4,3)=0
        cm5(4,4)=bx*cx + by*cy
        cm5(4,5)=bz*cy + by*cz
        cm5(4,6)=bz*cx + bx*cz
        cm5(5,1)=0
        cm5(5,2)=0
        cm5(5,3)=0
        cm5(5,4)=0
        cm5(5,5)=bx*cx + bz*cz
        cm5(5,6)=by*cx + bx*cy
        cm5(6,1)=0
        cm5(6,2)=0
        cm5(6,3)=0
        cm5(6,4)=0
        cm5(6,5)=0
        cm5(6,6)=by*cy + bz*cz

        cm6(1,1)=ax*cx
        cm6(1,2)=0
        cm6(1,3)=0
        cm6(1,4)=ay*cx + ax*cy
        cm6(1,5)=az*cx + ax*cz
        cm6(1,6)=0
        cm6(2,1)=0
        cm6(2,2)=ay*cy
        cm6(2,3)=0
        cm6(2,4)=ay*cx + ax*cy
        cm6(2,5)=0
        cm6(2,6)=az*cy + ay*cz
        cm6(3,1)=0
        cm6(3,2)=0
        cm6(3,3)=az*cz
        cm6(3,4)=0
        cm6(3,5)=az*cx + ax*cz
        cm6(3,6)=az*cy + ay*cz
        cm6(4,1)=0
        cm6(4,2)=0
        cm6(4,3)=0
        cm6(4,4)=ax*cx + ay*cy
        cm6(4,5)=az*cy + ay*cz
        cm6(4,6)=az*cx + ax*cz
        cm6(5,1)=0
        cm6(5,2)=0
        cm6(5,3)=0
        cm6(5,4)=0
        cm6(5,5)=ax*cx + az*cz
        cm6(5,6)=ay*cx + ax*cy
        cm6(6,1)=0
        cm6(6,2)=0
        cm6(6,3)=0
        cm6(6,4)=0
        cm6(6,5)=0
        cm6(6,6)=ay*cy + az*cz
    end subroutine
end subroutine

!> Atomic masses from atomic number
module elemental function lo_mass_from_Z(Z) result(mass)
    !> atomic number
    integer, intent(in) :: Z
    !> atomic mass
    real(flyt) :: mass

    real(flyt), parameter :: mass_of_z(103)=[&
        1837.3635941738078_flyt,7296.2969734687958_flyt,12650.914674013131_flyt,16428.203160097353_flyt,19707.298616949192_flyt,&
        21894.232166580099_flyt,25532.658002279986_flyt,29165.131013836126_flyt,34631.970469439431_flyt,36785.974179179480_flyt,&
        41907.785702600464_flyt,44305.398719494158_flyt,49184.336085177551_flyt,51196.732196748162_flyt,56461.713422756598_flyt,&
        58450.531435178993_flyt,64626.751682478913_flyt,72820.159896487196_flyt,71271.842915840840_flyt,73057.766007199767_flyt,&
        81949.614150126494_flyt,87255.746029294125_flyt,92860.614527125974_flyt,94783.152256241054_flyt,100145.92981946046_flyt,&
        101799.47291883679_flyt,107428.64256934724_flyt,106991.43472875118_flyt,115837.34411593614_flyt,119232.81727182020_flyt,&
        127097.37348457018_flyt,132414.33382488810_flyt,136573.71556256805_flyt,143934.15963176650_flyt,145655.22107549547_flyt,&
        152754.42102390452_flyt,155798.01983481817_flyt,159715.37166228655_flyt,162065.44714477652_flyt,166290.53681230778_flyt,&
        169357.96685242554_flyt,174883.53175668753_flyt,176648.83155961873_flyt,184230.12531840615_flyt,187585.25831892420_flyt,&
        193983.27749829018_flyt,196631.61010862543_flyt,204913.72412466665_flyt,209300.56637053887_flyt,216395.28954064651_flyt,&
        221954.51211759393_flyt,232606.27244494593_flyt,231332.70256887272_flyt,239332.06149861717_flyt,242271.81796656057_flyt,&
        250331.60985197860_flyt,253209.15261516630_flyt,255415.34392933775_flyt,256858.93778433156_flyt,262926.38304727047_flyt,&
        264154.77040951385_flyt,274101.08530984761_flyt,277014.10111818346_flyt,286653.08288060426_flyt,289703.18470917374_flyt,&
        296218.40489202098_flyt,300649.58503824903_flyt,304894.64416939113_flyt,307948.23220325267_flyt,315428.42281563085_flyt,&
        318944.82173007610_flyt,325358.20052182802_flyt,329847.79912501748_flyt,335123.06042591197_flyt,339434.06137212575_flyt,&
        346758.71246973320_flyt,350388.43593453628_flyt,355605.08972890303_flyt,359048.09007865738_flyt,365669.91099485726_flyt,&
        372568.16925712326_flyt,377733.31494506379_flyt,380947.96245039790_flyt,380947.23566475883_flyt,382788.35303494456_flyt,&
        404717.70151349297_flyt,406540.58999882534_flyt,412027.48433967575_flyt,413850.37282500812_flyt,422979.49916528584_flyt,&
        421152.65264218533_flyt,433900.15969218750_flyt,432115.71544803848_flyt,444894.16373021837_flyt,443071.27524488600_flyt,&
        450381.05807106884_flyt,450381.05807106884_flyt,457690.84089725168_flyt,459513.72938258405_flyt,468664.62957895256_flyt,&
        470487.51806428493_flyt,472310.40654961730_flyt,477797.30089046771_flyt]
    mass=mass_of_Z(Z)
end function

end submodule
