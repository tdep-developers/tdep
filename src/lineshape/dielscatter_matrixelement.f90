
! !> lowest order dipole matrix element
! function polarization_matrixelement_firstorder(di,p,Z,omega,eigenvector) result(M)
!     !> interaction tensors
!     type(lo_dielectric_tensor), intent(in) :: di
!     !> structure
!     type(lo_crystalstructure), intent(in) :: p
!     !> Born effective charges
!     real(r8), dimension(:,:,:), intent(in) :: Z
!     !> phonon frequency
!     real(r8), intent(in) :: omega
!     !> phonon eigenvector
!     complex(r8), dimension(:), intent(in) :: eigenvector
!     !> matrix element
!     complex(r8), dimension(3) :: M
!
!     real(r8), dimension(3) :: v0
!     integer :: iatom
!
!     ! need a nonzero frequency
!     if ( omega .lt. lo_freqtol ) then
!         M=0.0_r8
!         return
!     endif
!
!     M=0.0_r8
!     do iatom=1,p%na
!         v0=real(eigenvector( (iatom-1)*3+1:iatom*3 ),r8)
!         v0=v0*p%invsqrtmass(iatom)/sqrt(2.0_r8*omega)
!         M=M+matmul(Z(:,:,iatom),v0)
!     enddo
! end function

! !> second order dipole matrix element
! function polarization_matrixelement_secondorder(di,p,q,omega1,omega2,eigenvector1,eigenvector2) result(M)
!     !> interaction tensors
!     type(lo_dielectric_tensor), intent(in) :: di
!     !> structure
!     type(lo_crystalstructure), intent(in) :: p
!     !> q-vector (without 2*pi)
!     real(r8), dimension(3) :: q
!     !> phonon frequencies
!     real(r8), intent(in) :: omega1,omega2
!     !> phonon eigenvectors
!     complex(r8), dimension(:), intent(in) :: eigenvector1,eigenvector2
!     !> matrix element
!     complex(r8), dimension(3) :: M
!
!     complex(r8), dimension(3) :: cv0
!     real(r8) :: qdotr,pref,f0
!     complex(r8) :: expiqr
!     integer :: ipair,a1,a2,i,j,k,ii,jj
!
!     ! need nonzero frequencies
!     if ( omega1 .lt. lo_freqtol ) then
!         M=0.0_r8
!         return
!     endif
!     if ( omega2 .lt. lo_freqtol ) then
!         M=0.0_r8
!         return
!     endif
!
!     M=0.0_r8
!     do ipair=1,di%n_Z_pair
!         ! phase factor
!         qdotr=dot_product(di%Z_pair(ipair)%lv,q)*lo_twopi
!         expiqr=cmplx(cos(qdotr),sin(qdotr),r8)
!         a1=di%Z_pair(ipair)%a1
!         a2=di%Z_pair(ipair)%a2
!         ! Get the prefactor
!         f0=p%invsqrtmass(a1)*p%invsqrtmass(a2)
!         cv0=0.0_r8
!         do i=1,3
!         do j=1,3
!             ii=(a1-1)*3+i
!             jj=(a2-1)*3+j
!             do k=1,3
!                 cv0(k)=cv0(k)+eigenvector1(ii)*conjg(eigenvector2(jj))*di%Z_pair(ipair)%m(k,i,j)
!             enddo
!         enddo
!         enddo
!         M=M+cv0*expiqr*f0
!     enddo
!     ! and the common prefactor
!     pref=1.0_r8/sqrt(4*omega1*omega2)
!     M=M*pref
! end function

! function compound_ir_I(di,p,Z,omega,egv,mode) result(M)
!     !> interaction tensors
!     type(lo_dielectric_tensor), intent(in) :: di
!     !> structure
!     type(lo_crystalstructure), intent(in) :: p
!     !> Born effective charges
!     real(r8), dimension(:,:,:), intent(in) :: Z
!     !> frequencies at this q-vector (mode)
!     real(r8), dimension(:), intent(in) :: omega
!     !> eigenvectors at q (ialpha,mode)
!     complex(r8), dimension(:,:), intent(in) :: egv
!     !> mode 1
!     integer, intent(in) :: mode
!     !> compound matrix element
!     real(r8), dimension(3,3) :: M
!
!     complex(r8), dimension(3) :: M1
!     integer :: i,j
!
!     M1=polarization_matrixelement_firstorder(di,p,Z,omega(mode),egv(:,mode))
!     do j=1,3
!     do i=1,3
!         M(i,j)=real(conjg(M1(i))*M1(j),r8)
!     enddo
!     enddo
!     M=lo_chop(M,1E-13_r8)
! end function

! function compound_ir_II(di,p,q,omega,egv,mode1,mode2) result(M)
!     !> interaction tensors
!     type(lo_dielectric_tensor), intent(in) :: di
!     !> structure
!     type(lo_crystalstructure), intent(in) :: p
!     !> q-vector (without 2*pi)
!     real(r8), dimension(3) :: q
!     !> frequencies at this q-vector (mode)
!     real(r8), dimension(:), intent(in) :: omega
!     !> eigenvectors at q (ialpha,mode)
!     complex(r8), dimension(:,:), intent(in) :: egv
!     !> mode 1
!     integer, intent(in) :: mode1
!     !> mode 2
!     integer, intent(in) :: mode2
!     !> compound matrix element
!     real(r8), dimension(3,3) :: M
!
!     real(r8) :: prefactor=2.0_r8
!     complex(r8), dimension(3) :: M1
!     integer :: i,j,k,l
!
!     M1=polarization_matrixelement_secondorder(di,p,q,omega(mode1),omega(mode2),egv(:,mode1),egv(:,mode2))
!     do j=1,3
!     do i=1,3
!         M(i,j)=real( M1(i)*conjg(M1(j)),r8)
!     enddo
!     enddo
!     M=lo_chop(M,1E-13_r8)
!     ! And add the prefactor? Unsure about volume.
!     M=M*prefactor
! end function

! function compound_ir_III(di,p,fc,fct,q,omg,omq,egvg,egvq,mode1,mode2,mode3) result(M)
!     !> interaction tensors
!     type(lo_dielectric_tensor), intent(in) :: di
!     !> structure
!     type(lo_crystalstructure), intent(in) :: p
!     !> second order forceconstant
!     type(lo_forceconstant_secondorder), intent(in) :: fc
!     !> third order forceconstant
!     type(lo_forceconstant_thirdorder), intent(in) :: fct
!     !> q-vector (without 2*pi)
!     real(r8), dimension(3), intent(in) :: q
!     !> frequencies it gamma
!     real(r8), dimension(:), intent(in) :: omg
!     !> frequencies at q-vector (mode)
!     real(r8), dimension(:), intent(in) :: omq
!     !> eigenvectors at gamma
!     complex(r8), dimension(:,:), intent(in) :: egvg
!     !> eigenvectors at q (ialpha,mode)
!     complex(r8), dimension(:,:), intent(in) :: egvq
!     !> mode 1
!     integer, intent(in) :: mode1
!     !> mode 2
!     integer, intent(in) :: mode2
!     !> mode 3
!     integer, intent(in) :: mode3
!     !> compound matrix element
!     real(r8), dimension(3,3) :: M
!
!     real(r8) :: prefactor=-6.0_r8
!     complex(r8) :: v3
!     complex(r8), dimension(fc%na*3,3) :: egv
!     complex(r8), dimension(3,3) :: M3
!     complex(r8), dimension(3) :: M1,M2
!     real(r8), dimension(3) :: omega
!     real(r8), dimension(3) :: qv1,qv2,qv3
!     integer :: i,j
!
!     ! Now it's a little tricky. We will need many different guys here.
!     ! firstorder IR
!     M1=polarization_matrixelement_firstorder(di,p,fc%loto%born_effective_charges,omg(mode3),egvg(:,mode3))
!     ! secondorder IR
!     M2=polarization_matrixelement_secondorder(di,p,q,omq(mode1),omq(mode2),egvq(:,mode1),egvq(:,mode2))
!     ! thirdorder phonon
!     qv1=-q*lo_twopi
!     qv2=-qv1
!     qv3=0.0_r8
!     omega(1)=omq(mode1)
!     omega(2)=omq(mode2)
!     omega(3)=omg(mode3)
!     egv(:,1)=egvq(:,mode1)
!     egv(:,2)=conjg(egvq(:,mode2))
!     egv(:,3)=egvg(:,mode3)
!
!     if ( minval(omega) .lt. lo_freqtol ) then
!         v3=0.0_r8
!     else
!         v3=fct%scatteringamplitude(omega,egv,qv2,qv3)
!     endif
!     ! put it together
!     do i=1,3
!     do j=1,3
!         M3(i,j)=M1(i)*M2(j)*v3+M1(j)*M2(i)*v3
!     enddo
!     enddo
!
!     M=real(M3,r8)
!     M=min(M,0.0_r8)
!     M=lo_chop(M,1E-13_r8)
!     M=M*prefactor
! end function

! function compound_raman_I(di,p,omega,egv,mode) result(M)
!     !> interaction tensors
!     type(lo_dielectric_tensor), intent(in) :: di
!     !> structure
!     type(lo_crystalstructure), intent(in) :: p
!     !> frequencies at this q-vector (mode)
!     real(r8), dimension(:), intent(in) :: omega
!     !> eigenvectors at q (ialpha,mode)
!     complex(r8), dimension(:,:), intent(in) :: egv
!     !> mode 1
!     integer, intent(in) :: mode
!     !> compound matrix element
!     real(r8), dimension(3,3,3,3) :: M
!
!     complex(r8), dimension(3,3) :: M1
!     integer :: i,j,k,l
!
!     M1=raman_matrixelement_firstorder(di,p,omega(mode),egv(:,mode))
!     do l=1,3
!     do k=1,3
!     do j=1,3
!     do i=1,3
!         !M(i,j,k,l)=real( M1(i,j)*conjg(M1(k,l)),r8)
!         !M(i,j,k,l)=real( M1(i,j)*M1(k,l),r8)
!         M(i,j,k,l)=real( conjg(M1(i,j))*M1(k,l),r8)
!     enddo
!     enddo
!     enddo
!     enddo
!     M=lo_chop(M,1E-13_r8)
! end function

! function compound_raman_II(di,p,q,omega,egv,mode1,mode2) result(M)
!     !> interaction tensors
!     type(lo_dielectric_tensor), intent(in) :: di
!     !> structure
!     type(lo_crystalstructure), intent(in) :: p
!     !> q-vector (without 2*pi)
!     real(r8), dimension(3) :: q
!     !> frequencies at this q-vector (mode)
!     real(r8), dimension(:), intent(in) :: omega
!     !> eigenvectors at q (ialpha,mode)
!     complex(r8), dimension(:,:), intent(in) :: egv
!     !> mode 1
!     integer, intent(in) :: mode1
!     !> mode 2
!     integer, intent(in) :: mode2
!     !> compound matrix element
!     real(r8), dimension(3,3,3,3) :: M
!
!     complex(r8), dimension(3,3) :: M1
!     integer :: i,j,k,l
!
!     M1=raman_matrixelement_secondorder(di,p,q,omega(mode1),omega(mode2),egv(:,mode1),egv(:,mode2))
!     !write(*,*) 're',sum(abs(real(M1))),sum(abs(aimag(M1)))
!     do l=1,3
!     do k=1,3
!     do j=1,3
!     do i=1,3
!         M(i,j,k,l)=real( M1(i,j)*conjg(M1(k,l)) ,r8 )
!     enddo
!     enddo
!     enddo
!     enddo
!     M=lo_chop(M,1E-13_r8)
! end function

! function compound_raman_II_flat(di,p,q,omega,egv,mode1,mode2) result(M)
!     !> interaction tensors
!     type(lo_dielectric_tensor), intent(in) :: di
!     !> structure
!     type(lo_crystalstructure), intent(in) :: p
!     !> q-vector (without 2*pi)
!     real(r8), dimension(3) :: q
!     !> frequencies at this q-vector (mode)
!     real(r8), dimension(:), intent(in) :: omega
!     !> eigenvectors at q (ialpha,mode)
!     complex(r8), dimension(:,:), intent(in) :: egv
!     !> mode 1
!     integer, intent(in) :: mode1
!     !> mode 2
!     integer, intent(in) :: mode2
!     !> compound matrix element
!     real(r8), dimension(81) :: M
!
!     complex(r8), dimension(3,3) :: M1
!     integer :: i,j,k,l,ii
!
!     M1=raman_matrixelement_secondorder(di,p,q,omega(mode1),omega(mode2),egv(:,mode1),egv(:,mode2))
!     do l=1,3
!     do k=1,3
!     do j=1,3
!     do i=1,3
!         ii=l+(k-1)*3+(j-1)*9+(i-1)*27
!         M(ii)=real( M1(i,j)*conjg(M1(k,l)) ,r8 )
!     enddo
!     enddo
!     enddo
!     enddo
!     M=lo_chop(M,1E-13_r8)
! end function

! function raman_matrixelement_firstorder(di,p,omega,eigenvector) result(M)
!     !> interaction tensors
!     type(lo_dielectric_tensor), intent(in) :: di
!     !> structure
!     type(lo_crystalstructure), intent(in) :: p
!     !> phonon frequencies
!     real(r8), intent(in) :: omega
!     !> phonon eigenvectors
!     complex(r8), dimension(:), intent(in) :: eigenvector
!     !> matrix element
!     complex(r8), dimension(3,3) :: M
!
!     real(r8), dimension(3) :: v0
!     integer :: iatom,i
!
!     ! need a nonzero frequency
!     if ( omega .lt. lo_freqtol ) then
!         M=0.0_r8
!         return
!     endif
!
!     M=0.0_r8
!     do iatom=1,p%na
!         v0=real(eigenvector( (iatom-1)*3+1:iatom*3 ),r8)
!         v0=v0*p%invsqrtmass(iatom)/sqrt(2.0_r8*omega)
!         do i=1,3
!             M=M + di%eps_singlet(iatom)%m(:,:,i) * v0(i)
!         enddo
!     enddo
! end function

! !> second order dipole matrix element
! function raman_matrixelement_secondorder(di,p,q,omega1,omega2,eigenvector1,eigenvector2) result(M)
!     !> interaction tensors
!     type(lo_dielectric_tensor), intent(in) :: di
!     !> structure
!     type(lo_crystalstructure), intent(in) :: p
!     !> q-vector (without 2*pi)
!     real(r8), dimension(3) :: q
!     !> phonon frequencies
!     real(r8), intent(in) :: omega1,omega2
!     !> phonon eigenvectors
!     complex(r8), dimension(:), intent(in) :: eigenvector1,eigenvector2
!     !> matrix element
!     complex(r8), dimension(3,3) :: M
!
!     complex(r8), dimension(3,3) :: cv0
!     real(r8) :: qdotr,pref,f0
!     complex(r8) :: expiqr
!     integer :: ipair,a1,a2,i,j,k,l,ii,jj
!
!     ! need nonzero frequencies
!     if ( omega1 .lt. lo_freqtol ) then
!         M=0.0_r8
!         return
!     endif
!     if ( omega2 .lt. lo_freqtol ) then
!         M=0.0_r8
!         return
!     endif
!
!     M=0.0_r8
!     do ipair=1,di%n_eps_pair
!         ! phase factor (lv = r_j - r_i, I think)
!         qdotr=dot_product(di%eps_pair(ipair)%lv,q)*lo_twopi
!         expiqr=cmplx(cos(qdotr),sin(qdotr),r8)
!         a1=di%eps_pair(ipair)%a1
!         a2=di%eps_pair(ipair)%a2
!
!         ! Get the prefactor
!         f0=p%invsqrtmass(a1)*p%invsqrtmass(a2)
!         cv0=0.0_r8
!         do i=1,3
!         do j=1,3
!             ii=(a1-1)*3+i
!             jj=(a2-1)*3+j
!             do k=1,3
!             do l=1,3
!                 cv0(k,l)=cv0(k,l)+eigenvector1(ii)*conjg(eigenvector2(jj))*di%eps_pair(ipair)%m(k,l,i,j)
!             enddo
!             enddo
!         enddo
!         enddo
!         M=M+cv0*expiqr*f0
!     enddo
!     ! and the common prefactor
!     pref=1.0_r8/sqrt(4*omega1*omega2)
!     M=M*pref
! end function

! !> second order dipole matrix element
! function raman_matrixelement_ugv_secondorder(di,p,q,omega1,omega2,eigenvector1,eigenvector2) result(M)
!     !> interaction tensors
!     type(lo_dielectric_tensor), intent(in) :: di
!     !> structure
!     type(lo_crystalstructure), intent(in) :: p
!     !> q-vector (without 2*pi)
!     real(r8), dimension(3) :: q
!     !> phonon frequencies
!     real(r8), intent(in) :: omega1,omega2
!     !> phonon eigenvectors
!     complex(r8), dimension(:), intent(in) :: eigenvector1,eigenvector2
!     !> matrix element
!     complex(r8), dimension(3,3) :: M
!
!     complex(r8), dimension(3,3) :: cv0
!     real(r8) :: qdotr,pref,f0
!     complex(r8) :: expiqr
!     integer :: ipair,a1,a2,i,j,k,l,ii,jj
!
!     ! need nonzero frequencies
!     if ( omega1 .lt. lo_freqtol ) then
!         M=0.0_r8
!         return
!     endif
!     if ( omega2 .lt. lo_freqtol ) then
!         M=0.0_r8
!         return
!     endif
!
!     M=0.0_r8
!     do ipair=1,di%n_eps_pair
!         ! phase factor (lv = r_j - r_i, I think)
!         qdotr=dot_product(di%eps_pair(ipair)%lv,q)*lo_twopi
!         expiqr=cmplx(cos(qdotr),sin(qdotr),r8)
!         a1=di%eps_pair(ipair)%a1
!         a2=di%eps_pair(ipair)%a2
!
!         ! Get the prefactor
!         !f0=p%invsqrtmass(a1)*p%invsqrtmass(a2)
!         cv0=0.0_r8
!         do i=1,3
!         do j=1,3
!             ii=(a1-1)*3+i
!             jj=(a2-1)*3+j
!             do k=1,3
!             do l=1,3
!                 cv0(k,l)=cv0(k,l)+eigenvector1(ii)*conjg(eigenvector2(jj))*di%eps_pair(ipair)%m(k,l,i,j)
!             enddo
!             enddo
!         enddo
!         enddo
!         M=M+cv0*expiqr !*f0
!     enddo
!     ! and the common prefactor
!     !pref=1.0_r8/sqrt(4*omega1*omega2)
!     M=M !*pref
! end function

! function threephonon_matrixelement(fct) result(M)
!     complex(r8) :: M
!
!
! end function
