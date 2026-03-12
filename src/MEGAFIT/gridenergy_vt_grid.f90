!> finalize V-T grids
subroutine VTfinalize( gr,gs,map,filename,qgrid,quasiharmonic,pressurestep,dumpgrid,mw,mem )
    !> grid
    type(lo_gridenergy_vol_temp), intent(in) :: gr
    !> gridsim
    type(lo_gridsim), intent(inout) :: gs
    !> forcemap
    type(lo_forcemap), intent(inout) :: map
    !> filename
    character(len=*), intent(in) :: filename
    !> q-grid for density of states?
    integer, dimension(3) :: qgrid
    !> print quasiharmonic things
    logical, intent(in) :: quasiharmonic
    !> step in pressure
    real(r8), intent(in) :: pressurestep
    !> dump interpolated data on the grids
    logical, intent(in) :: dumpgrid
    !> mpi helper
    type(lo_mpi_helper) :: mw
    !> memory tracker
    type(lo_mem_helper) :: mem

    integer, parameter :: iporder=2            !@todo make this input
    integer, parameter :: iprexp=2             !@todo make this input
    real(r8), parameter :: ipscale=0.15_r8 !@todo make this input

    class(lo_eos_1d), allocatable, dimension(:) :: eos,qh_eos
    type(lo_hdf5_helper) :: h5
    type(lo_grid_interpolation) :: ip
    real(r8), dimension(:), allocatable :: pressure

    call mem%tick()

    ! open files for writing and things like that
    init: block
        if ( mw%talk ) then
            ! Open file for writing and start dumping stuff
            call h5%init(__FILE__,__LINE__)
            call h5%open_file('write',trim(filename))
            !call h5open_f(lo_status)
            !call h5fcreate_f(trim(filename), H5F_ACC_TRUNC_F, file_id, lo_status)
        endif
        if ( mw%talk ) then
            write(*,*) '... dumping interpolated data, qgrid: ',tochar(qgrid)
            if ( quasiharmonic ) write(*,*) '... also dumping quasiharmonic stuff'
        endif
    end block init

    ! Dump NVT-ensemble grid
    dumpgridNVT: block
        !real(r8), dimension(:,:,:), allocatable :: gm,gj,gtt,g0
        real(r8), dimension(:,:), allocatable :: gv,gt,ge,gridcoord,gridval
        !real(r8), dimension(2) :: gc,gcs
        integer :: i,j,k,l

        ! Start with NVT ensemble
        call mem%allocate(gv,[gr%nv,gr%nt],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%allocate(gt,[gr%nv,gr%nt],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%allocate(ge,[gr%nv,gr%nt],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        gv=0.0_r8
        gt=0.0_r8
        ge=0.0_r8
        do i=1,gr%nv
        do j=1,gr%nt
            gv(i,j)=gr%volume(i,j)
            gt(i,j)=gr%temperature(j)
        enddo
        enddo

        ! Also evaluate the magnetic things
        if ( gs%info%magnetic_pair_interactions ) then
        !    lo_allocate(gm(map%nmagsingletshells,gr%nv,gr%nt))
        !    lo_allocate(gj(map%ntheta_magpair,gr%nv,gr%nt))
        !    lo_allocate(gtt(map%ntheta_magpair_T,gr%nv,gr%nt))
        !    lo_allocate(g0(map%ntheta_magpair,gr%nv,gr%nt))
        !    gm=0.0_r8
        !    gj=0.0_r8
        !    gtt=0.0_r8
        !    g0=0.0_r8
        !    do i=1,gr%nv
        !    do j=1,gr%nt
        !        gc( gs%info%dim_volume )=gr%volume(i,j)
        !        gc( gs%info%dim_temperature )=gr%temperature(j)
        !        call gs%coordinate_transformation(gc,gcs)
        !        do k=1,map%nmagsingletshells
        !            gm(k,i,j)=gs%magpair%ipM0%eval(k,gcs,rexp=2)
        !        enddo
        !        do k=1,map%ntheta_magpair
        !            gj(k,i,j)=gs%magpair%ipJ%eval(k,gcs,rexp=2)
        !            g0(k,i,j)=gs%poly%eval( gcs,gs%magpair%coeff(:,k))
        !        enddo
        !        do k=1,map%ntheta_magpair_T
        !            gtt(k,i,j)=gs%magpair%ipT%eval(k,gcs,rexp=2)
        !        enddo
        !    enddo
        !    enddo
        endif

        ! Store this in a group
        if ( mw%talk ) then
            call h5%open_group('write','grid_NVT')

            call h5%store_data(gv*lo_volume_Bohr_to_A,  h5%group_id,'volumes',enhet='A^3/atom')
            call h5%store_data(gt,                      h5%group_id,'temperatures',enhet='K')
            call h5%store_data(gr%U*lo_Hartree_to_eV,   h5%group_id,'static_internal_energy',enhet='eV/atom')
            call h5%store_data(gr%U0*lo_Hartree_to_eV,  h5%group_id,'delta_U0',enhet='eV/atom')
            call h5%store_data(gr%fph*lo_Hartree_to_eV, h5%group_id,'phonon_free_energy',enhet='eV/atom')
            call h5%store_data(gr%ah3*lo_Hartree_to_eV, h5%group_id,'anharmonic_free_energy_thirdorder',enhet='eV/atom')
            call h5%store_data(gr%ah4*lo_Hartree_to_eV, h5%group_id,'anharmonic_free_energy_fourthorder',enhet='eV/atom')
            ge=gr%U+gr%U0+gr%fph+gr%ah3+gr%ah4
            call h5%store_data(ge*lo_Hartree_to_eV,     h5%group_id,'Helmholtz_free_energy',enhet='eV/atom')

            if ( gs%info%magnetic_pair_interactions ) then
                !call lo_h5_store_data(gm,                  group_id,'mean_magnetic_moments',enhet='Bohr_magneton')
                !call lo_h5_store_data(gj*lo_Hartree_to_eV, group_id,'irreducible_Jij',enhet='eV/bohr')
                !call lo_h5_store_data(g0*lo_Hartree_to_eV, group_id,'irreducible_Jij_poly',enhet='eV/bohr')
                !call lo_h5_store_data(gtt*lo_Hartree_to_eV, group_id,'irreducible_Tij',enhet='eV/bohr')
            endif

            ! And the same thing, quasiharmonic
            if ( quasiharmonic ) then
                call h5%store_data(gr%qh_fph*lo_Hartree_to_eV,h5%group_id,'quasiharmonic_phonon_free_energy',enhet='eV/atom')
                call h5%store_data(gr%qh_ah3*lo_Hartree_to_eV,h5%group_id,'quasiharmonic_anharmonic_free_energy_thirdorder',enhet='eV/atom')
                call h5%store_data(gr%qh_ah4*lo_Hartree_to_eV,h5%group_id,'quasiharmonic_anharmonic_free_energy_fourthorder',enhet='eV/atom')
            endif

            call h5%close_group()
        endif

        ! create interpolations for later
        call mem%allocate(gridcoord,[2,gr%nv*gr%nt],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%allocate(gridval  ,[gr%nv*gr%nt,4],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        gridcoord=0.0_r8
        gridval=0.0_r8
        l=0
        do i=1,gr%nv
        do j=1,gr%nt
            l=l+1
            gridcoord(:,l)=[gr%volume(i,j),gr%temperature(j)]
            gridval(l,1)=gr%U0(i,j)
            gridval(l,2)=gr%fph(i,j)
            gridval(l,3)=gr%ah3(i,j)
            gridval(l,4)=gr%ah4(i,j)
        enddo
        enddo
        call ip%generate(gridcoord,gridval,iporder,ipscale,[iporder,iporder])

        ! And a little cleanup
        call mem%deallocate(gv,       persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%deallocate(gt,       persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%deallocate(ge,       persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%deallocate(gridcoord,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%deallocate(gridval  ,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    end block dumpgridNVT

    ! get the equation of state at each temperature
    geteos: block
        real(r8), dimension(:,:), allocatable :: energy_part
        real(r8), dimension(gr%nv) :: dumv,dume,energy
        real(r8), dimension(gr%nt) :: dr
        real(r8) :: dist,f1,f2,f3,f4,f5,dl1,dl2,dl3,y
        integer :: i,j,l

        if ( mw%talk ) write(*,*) '... calculating thermal equation of state'
        ! do some EOS fitting? Which one? Same is in the input I assume for now.
        !@todo This should be determined automagically with some N-fold cross-validation
        select type(a=>gs%eos)
        type is(lo_eos_birch_murnaghan)
            allocate(lo_eos_birch_murnaghan::eos(gr%nt))
            allocate(lo_eos_birch_murnaghan::qh_eos(gr%nt))
        type is(lo_eos_vinet)
            allocate(lo_eos_vinet::eos(gr%nt))
            allocate(lo_eos_vinet::qh_eos(gr%nt))
        end select
        if ( mw%talk ) then
            write(*,"(1X,4(1X,A18))") 'Temperature (K)','F (eV/atom)','V (A^3/atom)','B (Gpa)'
        endif
        do i=1,gr%nt
            energy=0.0_r8
            ! get the fit, first add up the energies
            energy=gr%U(:,i)+gr%U0(:,i)+gr%fph(:,i)+gr%ah3(:,i)+gr%ah4(:,i)
            ! filter out the stable parts
            dumv=0.0_r8
            dume=0.0_r8
            l=0
            do j=1,gr%nv
                if ( gr%fph(j,i) .lt. 123456788_r8 ) then
                if ( gr%ah3(j,i) .lt. 123456788_r8 ) then
                if ( gr%ah4(j,i) .lt. 123456788_r8 ) then
                    l=l+1
                    dumv(l)=gr%volume(j,i)
                    dume(l)=energy(j)
                endif
                endif
                endif
            enddo
            if ( l .ge. 4 ) then
                ! we can fit
                select type(a=>eos(i))
                type is(lo_eos_birch_murnaghan)
                    call a%fit( dumv(1:l),dume(1:l),verbosity=0 )
                type is(lo_eos_vinet)
                    call a%fit( dumv(1:l),dume(1:l),verbosity=0 )
                end select
            else
                if ( mw%talk ) write(*,*) 'WARNING: few points for EOS fit at T=',tochar(gr%temperature(i))
                ! nothing stable. Use the 0K thing.
                select type(a=>eos(i))
                type is(lo_eos_birch_murnaghan)
                    call a%generate(  gs%eos%E0*lo_Hartree_to_eV,gs%eos%V0*lo_volume_bohr_to_A,gs%eos%B0*lo_pressure_hartreebohr_to_GPa,gs%eos%B0p )
                type is(lo_eos_vinet)
                    call a%generate(  gs%eos%E0*lo_Hartree_to_eV,gs%eos%V0*lo_volume_bohr_to_A,gs%eos%B0*lo_pressure_hartreebohr_to_GPa,gs%eos%B0p )
                end select
            endif

            if ( quasiharmonic ) then
                select type(a=>qh_eos(i)); type is(lo_eos_birch_murnaghan)
                    energy=gr%U(:,i)+gr%qh_fph(:,i)
                    call a%fit( gr%volume(:,i),energy,verbosity=0 )
                end select
            endif
            ! report a little
            if ( mw%talk ) write(*,"(1X,4(1X,F18.12))") gr%temperature(i),eos(i)%E0*lo_Hartree_to_eV,&
                eos(i)%V0*lo_volume_bohr_to_A,eos(i)%B0*lo_pressure_Hartreebohr_to_GPa
        enddo
        ! for the unstable points -- replace the EOS with that of the closest stable point.
        do i=1,gr%nt
            if ( abs(eos(i)%E0-gs%eos%E0) .lt. lo_sqtol ) then
                ! find closest stable point in temperature
                dist=lo_huge
                l=0
                do j=i,1,-1
                    if ( abs(eos(j)%E0-gs%eos%E0) .gt. lo_sqtol ) then
                    if ( abs(gr%temperature(j)-gr%temperature(i)) .lt. dist ) then
                        l=j
                        dist=abs(gr%temperature(j)-gr%temperature(i))
                    endif
                    endif
                enddo
                ! if so, replace the EOS
                if ( l .gt. 0 ) then
                select type(a=>eos(i)); type is(lo_eos_birch_murnaghan)
                    call a%generate(  eos(j)%E0*lo_Hartree_to_eV,eos(j)%V0*lo_volume_bohr_to_A,eos(j)%B0*lo_pressure_hartreebohr_to_GPa,eos(j)%B0p )
                end select
                endif
            endif
        enddo

        ! Write this to file.
        if ( mw%talk ) then
            ! start with the general stuff, like F and V as a function of T.
            call h5%store_data(gr%temperature,h5%file_id,'temperature',enhet='K')
            do i=1,gr%nt
                dr(i)=eos(i)%V0*lo_volume_bohr_to_A
            enddo
            call h5%store_data(dr,h5%file_id,'volume',enhet='A^3/atom')
            do i=1,gr%nt
                dr(i)=eos(i)%E0*lo_Hartree_to_eV
            enddo
            call h5%store_data(dr,h5%file_id,'Gibbs_free_energy',enhet='eV/atom')
            do i=1,gr%nt
                dr(i)=eos(i)%B0*lo_pressure_hartreebohr_to_GPa
            enddo
            call h5%store_data(dr,h5%file_id,'bulk_modulus',enhet='GPa')

            ! Get the energy chopped into parts
            call mem%allocate(energy_part,[gr%nt,5],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
            energy_part=0.0_r8
            do i=1,gr%nt
                !select type(a=>gs%eos); class is(lo_eos_1d)
                !    f1=a%energy_from_volume( eos(i)%V0 )
                !end select
                f1=gs%eos%energy_from_volume( eos(i)%V0 )
                f2=ip%eval(1,[ eos(i)%V0, gr%temperature(i) ],iprexp)
                f3=ip%eval(2,[ eos(i)%V0, gr%temperature(i) ],iprexp)
                f4=ip%eval(3,[ eos(i)%V0, gr%temperature(i) ],iprexp)
                f5=ip%eval(4,[ eos(i)%V0, gr%temperature(i) ],iprexp)
                ! use lagrange multipliers to make sure the parts add up to the total.
                y=eos(i)%E0-f1-f2-f3-f4-f5  ! total error
                dl1=f3-f3**2*( f3+f4+f5-y )/(f3**2+f4**2+f5**2)
                dl2=f4-f4**2*( f3+f4+f5-y )/(f3**2+f4**2+f5**2)
                dl3=f5-f5**2*( f3+f4+f5-y )/(f3**2+f4**2+f5**2)
                energy_part(i,1)=f1
                energy_part(i,2)=f2
                energy_part(i,3)=f3+dl1
                energy_part(i,4)=f4+dl2
                energy_part(i,5)=f5+dl3
            enddo
            ! store the partial energies
            call h5%store_data(energy_part(:,1)*lo_Hartree_to_eV,h5%file_id,'partial_energy_U',enhet='eV/atom')
            call h5%store_data(energy_part(:,2)*lo_Hartree_to_eV,h5%file_id,'partial_energy_U0',enhet='eV/atom')
            call h5%store_data(energy_part(:,3)*lo_Hartree_to_eV,h5%file_id,'partial_energy_fph',enhet='eV/atom')
            call h5%store_data(energy_part(:,4)*lo_Hartree_to_eV,h5%file_id,'partial_energy_ah3',enhet='eV/atom')
            call h5%store_data(energy_part(:,5)*lo_Hartree_to_eV,h5%file_id,'partial_energy_ah4',enhet='eV/atom')
            call mem%deallocate(energy_part,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)

            if ( quasiharmonic ) then
                ! store the same thing quasiharmonically
                do i=1,gr%nt
                    dr(i)=qh_eos(i)%V0*lo_volume_bohr_to_A
                enddo
                call h5%store_data(dr,h5%file_id,'quasiharmonic_volume',enhet='A^3/atom')
                do i=1,gr%nt
                    dr(i)=qh_eos(i)%E0*lo_Hartree_to_eV
                enddo
                call h5%store_data(dr,h5%file_id,'quasiharmonic_Gibbs_free_energy',enhet='eV/atom')
                do i=1,gr%nt
                    dr(i)=qh_eos(i)%B0*lo_pressure_hartreebohr_to_GPa
                enddo
                call h5%store_data(dr,h5%file_id,'quasiharmonic_bulk_modulus',enhet='GPa')
            endif
        endif
    end block geteos

    ! Dump NPT-ensemble grid
    dumpgridNPT: block
        real(r8), dimension(gr%nv) :: dumpr

        real(r8), dimension(:,:), allocatable :: gv,gt,ge,gp
        real(r8), dimension(:,:), allocatable :: g_U,g_U0,g_fph,g_ah3,g_ah4,g_pv
        real(r8) :: minP,maxP,dl1,dl2,dl3,f1,f2,f3,f4,f5,y
        !integer(HID_T) :: group_id
        integer :: i,j,l,t,npress

        ! Find a reasonable pressure range, and figure out the dimensions of the pressure grid
        minP=-lo_huge
        maxP=lo_huge
        do t=1,gr%nt
            dumpr=eos(t)%pressure_from_volume( gr%volume(:,t) )
            minP=max(minP,minval(dumpr))
            maxP=min(maxP,maxval(dumpr))
        enddo
        ! Count the number of pressures I get
        l=0
        do i=ceiling(minP/pressurestep),floor(maxP/pressurestep)
            l=l+1
        enddo
        npress=l
        ! Store actual pressures available
        call mem%allocate(pressure,npress,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        pressure=0.0_r8
        l=0
        do i=ceiling(minP/pressurestep),floor(maxP/pressurestep)
            l=l+1
            pressure(l)=i*pressurestep
        enddo

        if ( npress .lt. 3 .and. mw%talk ) then
            write(*,*) ''
            write(*,*) 'WARNING: found only ',tochar(npress),' points within the specified pressure interval.'
            write(*,*) '    perhaps lower the pressure step a little?'
            write(*,*) ''
        endif

        call mem%allocate(gv   ,[npress,gr%nt],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%allocate(gt   ,[npress,gr%nt],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%allocate(ge   ,[npress,gr%nt],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%allocate(gp   ,[npress,gr%nt],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%allocate(g_u  ,[npress,gr%nt],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%allocate(g_u0 ,[npress,gr%nt],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%allocate(g_fph,[npress,gr%nt],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%allocate(g_ah3,[npress,gr%nt],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%allocate(g_ah4,[npress,gr%nt],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%allocate(g_pv ,[npress,gr%nt],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        gv   =0.0_r8
        gt   =0.0_r8
        ge   =0.0_r8
        gp   =0.0_r8
        g_u  =0.0_r8
        g_u0 =0.0_r8
        g_fph=0.0_r8
        g_ah3=0.0_r8
        g_ah4=0.0_r8
        g_pv =0.0_r8

        ! invert the equation of state to get some stuff at these pressure points.
        do j=1,gr%nt
            minP=minval( eos(j)%pressure_from_volume(gr%volume(:,j)) )
            maxP=maxval( eos(j)%pressure_from_volume(gr%volume(:,j)) )
            do i=1,npress
                gp(i,j)=pressure(i)
                if ( pressure(i) .gt. minP .and. pressure(i) .lt. maxP ) then
                    ge(i,j)=eos(j)%energy_from_pressure( pressure(i) )
                else
                    ge(i,j)=-100.0_r8
                endif
                gv(i,j)=eos(j)%volume_from_pressure( pressure(i) )
                gt(i,j)=gr%temperature(j)
            enddo
        enddo

        ! Get the energy chopped into parts
        do i=1,npress
        do j=1,gr%nt
            select type(a=>gs%eos); class is(lo_eos_1d)
                f1=a%energy_from_volume( gv(i,j) )
            end select
            !f1=gs%eos%energy_from_volume( eos(i)%V0 )
            f2=ip%eval(1,[ gv(i,j), gt(i,j) ],iprexp)
            f3=ip%eval(2,[ gv(i,j), gt(i,j) ],iprexp)
            f4=ip%eval(3,[ gv(i,j), gt(i,j) ],iprexp)
            f5=ip%eval(4,[ gv(i,j), gt(i,j) ],iprexp)
            ! use lagrange multiplier to ensure it adds up to the total
            y=ge(i,j)-f1-f2-f3-f4-f5  ! total error
            dl1=f3-f3**2*( f3+f4+f5-y )/(f3**2+f4**2+f5**2)
            dl2=f4-f4**2*( f3+f4+f5-y )/(f3**2+f4**2+f5**2)
            dl3=f5-f5**2*( f3+f4+f5-y )/(f3**2+f4**2+f5**2)
            g_u  (i,j)=f1
            g_u0 (i,j)=f2
            g_fph(i,j)=f3+dl1
            g_ah3(i,j)=f4+dl2
            g_ah4(i,j)=f5+dl3
            g_pv (i,j)=gp(i,j)*gv(i,j)
        enddo
        enddo
        ! now add the PV term to the total
        ge=ge+g_pv

        ! Store NPT grid to a group
        if ( mw%talk ) then
            call h5%open_group('write','grid_NPT')
            ! store real data
            call h5%store_data(gv*lo_volume_Bohr_to_A            ,h5%group_id,'volumes',enhet='A^3/atom')
            call h5%store_data(gt                                ,h5%group_id,'temperatures',enhet='K')
            call h5%store_data(gp*lo_pressure_HartreeBohr_to_GPa ,h5%group_id,'pressures',enhet='GPa')
            call h5%store_data(ge*lo_Hartree_to_eV               ,h5%group_id,'Gibbs_free_energy',enhet='eV/atom')
            ! store partial data
            call h5%store_data(g_u*lo_Hartree_to_eV              ,h5%group_id,'partial_U'  ,enhet='eV/atom')
            call h5%store_data(g_u0*lo_Hartree_to_eV             ,h5%group_id,'partial_U0' ,enhet='eV/atom')
            call h5%store_data(g_fph*lo_Hartree_to_eV            ,h5%group_id,'partial_fph',enhet='eV/atom')
            call h5%store_data(g_ah3*lo_Hartree_to_eV            ,h5%group_id,'partial_ah3',enhet='eV/atom')
            call h5%store_data(g_ah4*lo_Hartree_to_eV            ,h5%group_id,'partial_ah4',enhet='eV/atom')
            call h5%store_data(g_pv*lo_Hartree_to_eV             ,h5%group_id,'partial_pv' ,enhet='eV/atom')
            call h5%close_group()
        endif

        call mem%deallocate(gv   ,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%deallocate(gt   ,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%deallocate(ge   ,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%deallocate(gp   ,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%deallocate(g_u  ,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%deallocate(g_u0 ,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%deallocate(g_fph,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%deallocate(g_ah3,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%deallocate(g_ah4,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%deallocate(g_pv ,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    end block dumpgridNPT

    ! Dump phonon dispersions at each temperature, and pressure.
    ! but only if I have provided a q-point path to do so, so I can be sure it's consistent.
    if ( lo_does_file_exist('infile.qpoints_dispersion') ) then
    dispersions: block
        integer, parameter :: n_pstep=2 ! maximally, how many steps in pressure from zero.
        type(lo_crystalstructure) :: p
        type(lo_forceconstant_secondorder) :: fc
        type(lo_forceconstant_thirdorder) :: fct
        type(lo_forceconstant_fourthorder) :: fcf
        type(lo_phonon_bandstructure) :: bs
        real(r8), dimension(:,:), allocatable :: pairconstraints,gv,gt,gp
        real(r8), dimension(gs%ndim) :: gridcoord
        real(r8) :: t0,f0
        integer :: itemp,ipress,ctr,nconstr,npress
        integer :: i,l
        !integer(HID_T) :: group_id,subgroup_id
        character(len=1000) :: groupname

        ! Count the number of pressures I should evaluate the dispersions on.
        npress=0
        do i=-n_pstep,n_pstep
            f0=i*pressurestep
            if ( f0-minval(pressure) .lt. -lo_sqtol ) cycle
            if ( f0-maxval(pressure) .gt. lo_sqtol ) cycle
            npress=npress+1
        enddo

        if ( npress .gt. 0 ) then
            ! Ok, I have some sensible points. First determine the grid I will
            ! calculate the dispersions on.
            call mem%allocate(gv,[npress,gr%nt],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
            call mem%allocate(gp,[npress,gr%nt],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
            call mem%allocate(gt,[npress,gr%nt],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
            gv=0.0_r8
            gp=0.0_r8
            gt=0.0_r8
            npress=0
            do i=-n_pstep,n_pstep
                f0=i*pressurestep
                if ( f0-minval(pressure) .lt. -lo_sqtol ) cycle
                if ( f0-maxval(pressure) .gt. lo_sqtol ) cycle
                npress=npress+1
                do itemp=1,gr%nt
                    gp(npress,itemp)=f0
                    gv(npress,itemp)=eos(itemp)%volume_from_pressure( f0 )
                    gt(npress,itemp)=gr%temperature(itemp)
                enddo
            enddo

            if ( mw%talk ) then
                write(*,*) ''
                write(*,*) 'Calculating phonon bandstructures at constant pressure.'
                write(*,*) '... pressures range from ',tochar(minval(gp)*lo_pressure_HartreeBohr_to_GPa),' to ',tochar(maxval(gp)*lo_pressure_HartreeBohr_to_GPa),' in steps of ',tochar(pressurestep*lo_pressure_HartreeBohr_to_GPa),' GPa'
            endif

            if ( mw%talk ) then
                ! create a group for the dispersions
                call h5%open_group('write','NPT_dispersions')
                ! Store the relevant grids
                call h5%store_data(gv*lo_volume_Bohr_to_A            ,h5%group_id,'volumes',enhet='A^3/atom')
                call h5%store_data(gp*lo_pressure_HartreeBohr_to_GPa ,h5%group_id,'pressures',enhet='GPa')
                call h5%store_data(gt                                ,h5%group_id,'temperatures',enhet='K')
                ! Also start the progress bar since this could take a little while.
                call lo_progressbar_init()
                t0=walltime()
            endif

            ctr=0
            do ipress=1,npress
            do itemp=1,gr%nt
                ctr=ctr+1
                ! Interpolate structure and forceconstant here
                gridcoord( gs%info%dim_temperature ) = gt(ipress,itemp)
                gridcoord( gs%info%dim_volume ) = gv(ipress,itemp)
                call gs%structure%interpolate(gridcoord,p,gs%info%dim_volume)
                call p%classify('wedge',timereversal=.true.)
                !call lo_secondorder_rot_herm_huang( map,p,pairconstraints,nconstr,.true.,.true.,.true. )
nconstr=0
                if ( nconstr .gt. 0 ) then
                    call gs%eval(map,gridcoord,pairconstraints)
                else
                    call gs%eval(map,gridcoord)
                endif
                call map%get_secondorder_forceconstant(p,fc,mem,-1)
                call bs%generate(p,fc,timereversal=.true.,mw=mw,mem=mem,verbosity=-1,npts=100,readpathfromfile=.true.)
                !@TODO add destructors here?

                ! Store dispersions in their own subgroup
                if ( mw%talk ) then
                    groupname='pressure_'//tochar(ipress)//'_temperature_'//tochar(itemp)
                    call h5%open_subgroup('write',trim(groupname))
                    call bs%write_to_hdf5(p,enhet='thz',filename='null',mem=mem,hdftag=h5%subgroup_id)
                    call h5%close_subgroup()

                    if ( dumpgrid ) then
                        ! Possibly dump the explicit input files for further processing if necessary.
                        call p%writetofile('NPT_uc_'//tochar(ipress)//'_'//tochar(itemp),1)
                        call fc%writetofile(p,'NPT_fc2_'//tochar(ipress)//'_'//tochar(itemp))
                        if ( map%have_fc_triplet ) then
                            call map%get_thirdorder_forceconstant(p,fct)
                            call fct%writetofile(p,'NPT_fc3_'//tochar(ipress)//'_'//tochar(itemp))
                        endif
                        if ( map%have_fc_quartet ) then
                            call map%get_fourthorder_forceconstant(p,fcf)
                            call fcf%writetofile(p,'NPT_fc4_'//tochar(ipress)//'_'//tochar(itemp))
                        endif
                    endif
                endif

                ! Do things quasiharmonically? But with the anharmonic thermal expansion for now, way
                ! too annoying otherwise. Will revisit this in the future is someone cares.
                if ( quasiharmonic ) then
                    gridcoord( gs%info%dim_temperature ) = minval(gr%temperature)
                    gridcoord( gs%info%dim_volume ) = gv(ipress,itemp)
                    call gs%structure%interpolate(gridcoord,p,gs%info%dim_volume)
                    call p%classify('wedge',timereversal=.true.)
                    !call lo_secondorder_rot_herm_huang( map,p,pairconstraints,nconstr,.true.,.true.,.true. )
nconstr=0
                    if ( nconstr .gt. 0 ) then
                        call gs%eval(map,gridcoord,pairconstraints)
                    else
                        call gs%eval(map,gridcoord)
                    endif
                    call map%get_secondorder_forceconstant(p,fc,mem,-1)
                    call bs%generate(p,fc,timereversal=.true.,mw=mw,mem=mem,verbosity=-1,npts=100,readpathfromfile=.true.)

                    if ( mw%talk ) then
                        groupname='qh_pressure_'//tochar(ipress)//'_temperature_'//tochar(itemp)
                        call h5%open_subgroup('write',trim(groupname))
                        call bs%write_to_hdf5(p,enhet='thz',filename='null',mem=mem,hdftag=h5%subgroup_id)
                        call h5%close_subgroup()

                        if ( dumpgrid ) then
                            ! Possibly dump the explicit input files for further processing if necessary.
                            call p%writetofile('qh_NPT_uc_'//tochar(ipress)//'_'//tochar(itemp),1)
                            call fc%writetofile(p,'qh_NPT_fc2_'//tochar(ipress)//'_'//tochar(itemp))
                            if ( map%have_fc_triplet ) then
                                call map%get_thirdorder_forceconstant(p,fct)
                                call fct%writetofile(p,'qh_NPT_fc3_'//tochar(ipress)//'_'//tochar(itemp))
                            endif
                            if ( map%have_fc_quartet ) then
                                call map%get_fourthorder_forceconstant(p,fcf)
                                call fcf%writetofile(p,'qh_NPT_fc4_'//tochar(ipress)//'_'//tochar(itemp))
                            endif
                        endif
                    endif
                endif

                if ( mw%talk ) call lo_progressbar(' ... NPT phonon dispersions',ctr,size(gv),walltime()-t0)
            enddo
            enddo
            ! and close the dispersion group
            if ( mw%talk ) call h5%close_group()

            ! Now I might as well do it at constant volume as well. Will pick the volumes that
            ! correspond to the lowest temperature, suppose that makes sense.
            do ipress=1,npress
            do itemp=2,gr%nt
                gv(ipress,itemp)=gv(ipress,1)
                gp(ipress,itemp)=eos(itemp)%pressure_from_volume( gv(ipress,1) )
                gt(ipress,itemp)=gr%temperature(itemp)
            enddo
            enddo

            if ( mw%talk ) then
                write(*,*) ''
                write(*,*) 'Calculating phonon bandstructures at constant volume.'
                write(*,*) '... volumes range from ',tochar(minval(gv)*lo_volume_Bohr_to_A),' to ',tochar(maxval(gv)*lo_volume_Bohr_to_A),' A^3'
            endif

            if ( mw%talk ) then
                ! create a group for the dispersions
                call h5%open_group('write','NVT_dispersions')
                ! Store the relevant grids
                call h5%store_data(gv*lo_volume_Bohr_to_A            ,h5%group_id,'volumes',enhet='A^3/atom')
                call h5%store_data(gp*lo_pressure_HartreeBohr_to_GPa ,h5%group_id,'pressures',enhet='GPa')
                call h5%store_data(gt                                ,h5%group_id,'temperatures',enhet='K')
                ! Also start the progress bar since this could take a little while.
                call lo_progressbar_init()
                t0=walltime()
            endif

            ctr=0
            do ipress=1,npress
            do itemp=1,gr%nt
                ctr=ctr+1

                ! Interpolate structure and forceconstant here
                gridcoord( gs%info%dim_temperature ) = gt(ipress,itemp)
                gridcoord( gs%info%dim_volume ) = gv(ipress,itemp)
                call gs%structure%interpolate(gridcoord,p,gs%info%dim_volume)
                call p%classify('wedge',timereversal=.true.)
                !call lo_secondorder_rot_herm_huang( map,p,pairconstraints,nconstr,.true.,.true.,.true. )
nconstr=0
                if ( nconstr .gt. 0 ) then
                    call gs%eval(map,gridcoord,pairconstraints)
                else
                    call gs%eval(map,gridcoord)
                endif
                call map%get_secondorder_forceconstant(p,fc,mem,-1)
                call bs%generate(p,fc,timereversal=.true.,mw=mw,mem=mem,verbosity=-1,npts=100,readpathfromfile=.true.)
                !@TODO add destructors here?

                ! Store dispersions in their own subgroup
                if ( mw%talk ) then
                    groupname='volume_'//tochar(ipress)//'_temperature_'//tochar(itemp)
                    call h5%open_subgroup('write',trim(groupname))
                    call bs%write_to_hdf5(p,enhet='thz',filename='null',mem=mem,hdftag=h5%subgroup_id)
                    call h5%close_subgroup()

                    if ( dumpgrid ) then
                        ! Possibly dump the explicit input files for further processing if necessary.
                        call p%writetofile('NVT_uc_'//tochar(ipress)//'_'//tochar(itemp),1)
                        call fc%writetofile(p,'NVT_fc2_'//tochar(ipress)//'_'//tochar(itemp))
                        if ( map%have_fc_triplet ) then
                            call map%get_thirdorder_forceconstant(p,fct)
                            call fct%writetofile(p,'NVT_fc3_'//tochar(ipress)//'_'//tochar(itemp))
                        endif
                        if ( map%have_fc_quartet ) then
                            call map%get_fourthorder_forceconstant(p,fcf)
                            call fcf%writetofile(p,'NVT_fc4_'//tochar(ipress)//'_'//tochar(itemp))
                        endif
                    endif
                endif

                if ( mw%talk ) call lo_progressbar(' ... NVT phonon dispersions',ctr,size(gv),walltime()-t0)
            enddo
            enddo
            ! and close the dispersion group
            if ( mw%talk ) call h5%close_group()

            call mem%deallocate(gv,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
            call mem%deallocate(gp,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
            call mem%deallocate(gt,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        else
            ! Found not sensible pressures. Throw a warning.
            if ( mw%talk ) then
                write(*,*) ''
                write(*,*) 'WARNING: found no pressures valid for the entire grid around zero pressure.'
                write(*,*) 'Not that big of a deal, but you should look into it. Probably means that your'
                write(*,*) 'grid does not cover the equilibrium volume.'
                write(*,*) ''
            endif
        endif

        ! If I have done quasiharmonic things, repeat the procedure here.
        if ( quasiharmonic ) then
        endif
    end block dispersions
    endif

    ! And we are done with pressures
    call mem%deallocate(pressure,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)

    call mem%tock(__FILE__,__LINE__,mw%comm)

    ! Close files
    if ( mw%talk ) then
        call h5%close_file()
        call h5%destroy(__FILE__,__LINE__)
    endif

end subroutine
