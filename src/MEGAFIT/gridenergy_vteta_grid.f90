
!> finalize V-T grids
subroutine VTetafinalize( gr,gs,map,filename,qgrid,quasiharmonic,pressurestep,mw,mem )
    !> grid
    type(lo_gridenergy_vol_temp_eta), intent(in) :: gr
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
    real(flyt), intent(in) :: pressurestep
    !> mpi helper
    type(lo_mpi_helper) :: mw
    !> memory tracker
    type(lo_mem_helper) :: mem

    integer, parameter :: iporder=2 !@todo make this input
    integer, parameter :: iprexp=2 !@todo make this input
    real(flyt), parameter :: ipscale=0.15_flyt !@todo make this input

    class(lo_eos_2d), allocatable, dimension(:) :: eos,qh_eos
    integer(HID_T) :: file_id
    type(lo_grid_interpolation) :: ip

    ! open files for writing and things like that
    init: block
        if ( mw%talk ) then
            ! Open file for writing and start dumping stuff
            call h5open_f(lo_status)
            call h5fcreate_f(trim(filename), H5F_ACC_TRUNC_F, file_id, lo_status)
        endif
        if ( mw%talk ) then
            write(*,*) '... dumping interpolated data'
        endif
    end block init

    ! Dump NVT-ensemble grid
    dumpgridNVT: block
        real(flyt), dimension(:,:,:), allocatable :: gv,gt,ge,geta
        real(flyt), dimension(:,:), allocatable :: gridcoord,gridval
        integer(HID_T) :: group_id
        integer :: i,j,k,l

        ! Start with NVT ensemble
        lo_allocate(gv  (gr%nv,gr%nt,gr%neta))
        lo_allocate(gt  (gr%nv,gr%nt,gr%neta))
        lo_allocate(ge  (gr%nv,gr%nt,gr%neta))
        lo_allocate(geta(gr%nv,gr%nt,gr%neta))
        gv  =0.0_flyt
        gt  =0.0_flyt
        ge  =0.0_flyt
        geta=0.0_flyt
        do i=1,gr%nv
        do j=1,gr%nt
        do k=1,gr%neta
            gv(i,j,k)=gr%volume(i)
            gt(i,j,k)=gr%temperature(j)
            geta(i,j,k)=gr%eta(k)
        enddo
        enddo
        enddo
        ! write this to file
        if ( mw%talk ) then
            call h5gcreate_f(file_id,'grid_NVT',group_id,lo_status)
            ! store this as a group
            call lo_h5_store_data(gv*lo_volume_Bohr_to_A,  group_id,'volumes',enhet='A^3/atom')
            call lo_h5_store_data(gt,                      group_id,'temperatures',enhet='K')
            call lo_h5_store_data(geta,                    group_id,'eta',enhet='dimensionless')
            call lo_h5_store_data(gr%U*lo_Hartree_to_eV,   group_id,'static_internal_energy',enhet='eV/atom')
            call lo_h5_store_data(gr%U0*lo_Hartree_to_eV,  group_id,'delta_U0',enhet='eV/atom')
            call lo_h5_store_data(gr%fph*lo_Hartree_to_eV, group_id,'phonon_free_energy',enhet='eV/atom')
            call lo_h5_store_data(gr%ah3*lo_Hartree_to_eV, group_id,'anharmonic_free_energy_thirdorder',enhet='eV/atom')
            call lo_h5_store_data(gr%ah4*lo_Hartree_to_eV, group_id,'anharmonic_free_energy_fourthorder',enhet='eV/atom')
            ge=gr%U+gr%U0+gr%fph+gr%ah3+gr%ah4
            call lo_h5_store_data(ge*lo_Hartree_to_eV,     group_id,'Helmholtz_free_energy',enhet='eV/atom')

            if ( quasiharmonic ) then
                call lo_h5_store_data(gr%qh_fph*lo_Hartree_to_eV,group_id,'quasiharmonic_phonon_free_energy',enhet='eV/atom')
                call lo_h5_store_data(gr%qh_ah3*lo_Hartree_to_eV,group_id,'quasiharmonic_anharmonic_free_energy_thirdorder',enhet='eV/atom')
                call lo_h5_store_data(gr%qh_ah4*lo_Hartree_to_eV,group_id,'quasiharmonic_anharmonic_free_energy_fourthorder',enhet='eV/atom')
            endif

            call h5gclose_f(group_id,lo_status)
        endif

        ! create interpolations for later
        lo_allocate(gridcoord(2,gr%nv*gr%nt*gr%neta))
        lo_allocate(gridval(gr%nv*gr%nt*gr%neta,4))
        l=0
        do i=1,gr%nv
        do j=1,gr%nt
        do k=1,gr%neta
            l=l+1
            gridcoord(:,l)=[gr%volume(i),gr%temperature(j),gr%eta(k)]
            gridval(l,1)=gr%U0(i,j,k)
            gridval(l,2)=gr%fph(i,j,k)
            gridval(l,3)=gr%ah3(i,j,k)
            gridval(l,4)=gr%ah4(i,j,k)
        enddo
        enddo
        enddo
        call ip%generate(gridcoord,gridval,iporder,ipscale,[iporder,iporder,iporder])
        if ( mw%talk ) write(*,*) '... stored grid data and created interpolations'
    end block dumpgridNVT

    ! get the equation of state at each temperature
    geteos: block
        real(flyt), dimension(:,:), allocatable :: energy_part
        real(flyt), dimension(:), allocatable :: energy,volume,eta,dr
        real(flyt) :: f1,f2,f3,f4,f5,dl1,dl2,dl3,y
        integer :: i,j,l,t
        if ( mw%talk ) then
            write(*,*) '... calculating thermal equation of state'
        endif
        ! do some EOS fitting? Which one? Same is in the input I assume for now.
        select type(a=>gs%eos)
        type is(lo_eos_2d_birch_murnaghan)
            allocate(lo_eos_2d_birch_murnaghan::eos(gr%nt))
            allocate(lo_eos_2d_birch_murnaghan::qh_eos(gr%nt))
        end select
        lo_allocate(energy( gr%nv*gr%neta ))
        lo_allocate(volume( gr%nv*gr%neta ))
        lo_allocate(eta( gr%nv*gr%neta ))

        do t=1,gr%nt
            select type(a=>eos(t)); type is(lo_eos_2d_birch_murnaghan)
                ! add up the energy
                l=0
                do i=1,gr%nv
                do j=1,gr%neta
                    if ( gr%fph(i,t,j) .lt. 123456788_flyt ) then
                        l=l+1
                        volume(l)=gr%volume(i)
                        eta(l)=gr%eta(j)
                        energy(l)=gr%U(i,t,j)+gr%U0(i,t,j)+gr%fph(i,t,j)+gr%ah3(i,t,j)+gr%ah4(i,t,j)
                    endif
                enddo
                enddo

                if ( l .ge. 10 ) then
                    ! I need quite a few points to do a reliable fit
                    call a%fit( volume(1:l),energy(1:l),eta(1:l),verbosity=0 )
                else
                    ! just make a copy of the zero K one.
                    if ( mw%talk ) write(*,*) 'warning: insufficient stable points for temperature ',tochar(t)
                    select type(b=>gs%eos); class is(lo_eos_2d)
                    call a%generate_2d( b%E0*lo_hartree_to_eV,b%V0*lo_volume_bohr_to_A,b%eta0,&
                        b%B0*lo_pressure_Hartreebohr_to_GPa,b%B0p,b%C0,b%C1,b%C2,b%C3 )
                    end select
                endif
            end select
            ! if ( quasiharmonic ) then
            !     select type(a=>qh_eos(i)); type is(lo_eos_birch_murnaghan)
            !         energy=gr%U(:,i)+gr%qh_fph(:,i) !+gr%qh_ah3(:,i)+gr%qh_ah4(:,i)
            !         call a%fit( gr%volume,energy,verbosity=0 )
            !     end select
            ! endif
            ! report a little
            if ( mw%talk ) write(*,"(1X,5(1X,F18.12))") gr%temperature(t),eos(t)%E0*lo_Hartree_to_eV,eos(t)%V0*lo_volume_bohr_to_A,&
                eos(t)%eta0,eos(t)%B0*lo_pressure_hartreebohr_to_GPa
        enddo

        ! Write this to file.
        if ( mw%talk ) then
            ! start with the general stuff, like F and V as a function of T.
            lo_allocate(dr(gr%nt))
            call lo_h5_store_data(gr%temperature,file_id,'temperature',enhet='K')
            do i=1,gr%nt
                dr(i)=eos(i)%V0*lo_volume_bohr_to_A
            enddo
            call lo_h5_store_data(dr,file_id,'volume',enhet='A^3/atom')
            do i=1,gr%nt
                dr(i)=eos(i)%E0*lo_Hartree_to_eV
            enddo
            call lo_h5_store_data(dr,file_id,'Gibbs_free_energy',enhet='eV/atom')
            do i=1,gr%nt
                dr(i)=eos(i)%B0*lo_pressure_Hartreebohr_to_GPa
            enddo
            call lo_h5_store_data(dr,file_id,'bulk_modulus',enhet='GPa')
            do i=1,gr%nt
                dr(i)=eos(i)%eta0
            enddo
            call lo_h5_store_data(dr,file_id,'structural_parameter',enhet='dimensionless')

            ! Get the energy chopped into parts
            lo_allocate(energy_part(gr%nt,5))
            energy_part=0.0_flyt
            do i=1,gr%nt
                select type(a=>gs%eos); class is(lo_eos_2d)
                    f1=a%energy_from_volume_eta( eos(i)%V0, eos(i)%eta0 )
                end select
                f2=ip%eval (1,[ eos(i)%V0, gr%temperature(i), eos(i)%eta0 ],iprexp)
                f3=ip%eval (2,[ eos(i)%V0, gr%temperature(i), eos(i)%eta0 ],iprexp)
                f4=ip%eval (3,[ eos(i)%V0, gr%temperature(i), eos(i)%eta0 ],iprexp)
                f5=ip%eval (4,[ eos(i)%V0, gr%temperature(i), eos(i)%eta0 ],iprexp)
                ! use lagrange multiplier to fix
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
            call lo_h5_store_data(energy_part(:,1)*lo_Hartree_to_eV,file_id,'partial_energy_U',enhet='eV/atom')
            call lo_h5_store_data(energy_part(:,2)*lo_Hartree_to_eV,file_id,'partial_energy_U0',enhet='eV/atom')
            call lo_h5_store_data(energy_part(:,3)*lo_Hartree_to_eV,file_id,'partial_energy_fph',enhet='eV/atom')
            call lo_h5_store_data(energy_part(:,4)*lo_Hartree_to_eV,file_id,'partial_energy_ah3',enhet='eV/atom')
            call lo_h5_store_data(energy_part(:,5)*lo_Hartree_to_eV,file_id,'partial_energy_ah4',enhet='eV/atom')
            lo_deallocate(dr)
        endif
    end block geteos
    call mpi_barrier(mw%comm,mw%error)

    ! Dump NPT-ensemble grid
    dumpgridNPT: block
        real(flyt), dimension(gr%nv) :: dumpr
        real(flyt), dimension(:), allocatable :: pressure
        real(flyt), dimension(:,:), allocatable :: gv,gt,ge,gp,geta
        real(flyt), dimension(:,:), allocatable :: g_U,g_U0,g_fph,g_ah3,g_ah4,g_pv
        real(flyt) :: minP,maxP,dl1,dl2,dl3,f1,f2,f3,f4,f5,y
        integer(HID_T) :: group_id
        integer :: i,j,t,l,npress

        ! Find a reasonable pressure range, and figure out the dimensions of the pressure grid
        minP=lo_huge
        maxP=-lo_huge
        do t=1,gr%nt
            dumpr=eos(t)%pressure_from_volume( gr%volume )
            minP=min(minP,minval(dumpr))
            maxP=max(maxP,maxval(dumpr))
        enddo
        ! Count get pressures
        l=0
        do i=ceiling(minP/pressurestep),floor(maxP/pressurestep)
            l=l+1
        enddo
        npress=l
        lo_allocate(pressure(npress))
        l=0
        do i=ceiling(minP/pressurestep),floor(maxP/pressurestep)
            l=l+1
            pressure(l)=i*pressurestep
        enddo

        lo_allocate(gv(npress,gr%nt))
        lo_allocate(gt(npress,gr%nt))
        lo_allocate(ge(npress,gr%nt))
        lo_allocate(gp(npress,gr%nt))
        lo_allocate(geta(npress,gr%nt))
        gv=0.0_flyt
        gt=0.0_flyt
        ge=0.0_flyt
        gp=0.0_flyt
        geta=0.0_flyt
        ! invert the equation of state to get some stuff at these pressure points.
        do j=1,gr%nt
            minP=minval( eos(j)%pressure_from_volume(gr%volume) )
            maxP=maxval( eos(j)%pressure_from_volume(gr%volume) )
            do i=1,npress
                gp(i,j)=pressure(i)
                select type(a=>eos(j)); class is(lo_eos_2d)
                    if ( pressure(i) .gt. minP .and. pressure(i) .lt. maxP ) then
                        ge(i,j)=a%energy_from_pressure( pressure(i) )
                    else
                        ge(i,j)=-100.0_flyt
                    endif
                    gv(i,j)=a%volume_from_pressure( pressure(i) )
                    gt(i,j)=gr%temperature(j)
                    ge(i,j)=a%energy_from_pressure( pressure(i) )
                    geta(i,j)=a%eta_from_volume( gv(i,j) )
                end select
            enddo
        enddo
        ! Get the energy chopped into parts
        lo_allocate(g_u(  npress,gr%nt))
        lo_allocate(g_u0( npress,gr%nt))
        lo_allocate(g_fph(npress,gr%nt))
        lo_allocate(g_ah3(npress,gr%nt))
        lo_allocate(g_ah4(npress,gr%nt))
        lo_allocate(g_pv( npress,gr%nt))
        g_u=0.0_flyt
        g_u0=0.0_flyt
        g_fph=0.0_flyt
        g_ah3=0.0_flyt
        g_ah4=0.0_flyt
        g_pv=0.0_flyt

        do i=1,npress
        do j=1,gr%nt
            select type(a=>gs%eos); class is(lo_eos_2d)
                f1=a%energy_from_volume_eta( gv(i,j), geta(i,j) )
            end select
            f2=ip%eval(1,[ gv(i,j), gt(i,j), geta(i,j) ],iprexp)
            f3=ip%eval(2,[ gv(i,j), gt(i,j), geta(i,j) ],iprexp)
            f4=ip%eval(3,[ gv(i,j), gt(i,j), geta(i,j) ],iprexp)
            f5=ip%eval(4,[ gv(i,j), gt(i,j), geta(i,j) ],iprexp)
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

        ! Store this in a group
        if ( mw%talk ) then
            call h5gcreate_f(file_id,'grid_NPT',group_id,lo_status)
            ! store real data
            call lo_h5_store_data(gv*lo_volume_Bohr_to_A,             group_id,'volumes',enhet='A^3/atom')
            call lo_h5_store_data(gt,                                 group_id,'temperatures',enhet='K')
            call lo_h5_store_data(gp*lo_pressure_HartreeBohr_to_GPa,  group_id,'pressures',enhet='GPa')
            call lo_h5_store_data(ge*lo_Hartree_to_eV,                group_id,'Gibbs_free_energy',enhet='eV/atom')
            ! store partial data
            call lo_h5_store_data(g_u*lo_Hartree_to_eV  , group_id,'partial_U'  ,enhet='eV/atom')
            call lo_h5_store_data(g_u0*lo_Hartree_to_eV , group_id,'partial_U0' ,enhet='eV/atom')
            call lo_h5_store_data(g_fph*lo_Hartree_to_eV, group_id,'partial_fph',enhet='eV/atom')
            call lo_h5_store_data(g_ah3*lo_Hartree_to_eV, group_id,'partial_ah3',enhet='eV/atom')
            call lo_h5_store_data(g_ah4*lo_Hartree_to_eV, group_id,'partial_ah4',enhet='eV/atom')
            call lo_h5_store_data(g_pv*lo_Hartree_to_eV , group_id,'partial_pv' ,enhet='eV/atom')

            call h5gclose_f(group_id,lo_status)
        endif
    end block dumpgridNPT

    ! Dump phonon dispersions at each temperature, at the equilibrium volume.
    ! but only if I have provided a q-point path to do so, so I can be sure it's consistent.
    if ( lo_does_file_exist('infile.qpoints_dispersion') ) then
    dispersions: block
        type(lo_crystalstructure) :: p
        type(lo_forceconstant_secondorder) :: fc
        type(lo_phonon_bandstructure) :: bs
        real(flyt), dimension(:,:), allocatable :: pairconstraints
        real(flyt), dimension(gs%ndim) :: gridcoord
        real(flyt) :: t0
        integer :: t,nconstr
        integer(HID_T) :: group_id,subgroup_id
        character(len=1000) :: groupname

        if ( mw%talk ) then
            ! create a group for the dispersions
            call h5gcreate_f(file_id,'dispersions',group_id,lo_status)
            call lo_progressbar_init()
            t0=walltime()
        endif

        do t=1,gr%nt
            gridcoord( gs%info%dim_temperature ) = gr%temperature(t)
            gridcoord( gs%info%dim_volume ) = eos(t)%V0
            gridcoord( gs%info%dim_eta ) = eos(t)%eta0
            ! get a structure and a forceconstant
            call gs%structure%interpolate(gridcoord,p,gs%info%dim_volume)
            !call lo_secondorder_rot_herm_huang( map,p,pairconstraints,nconstr,.true.,.true.,.true. )
nconstr=0
            if ( nconstr .gt. 0 ) then
                call gs%eval(map,gridcoord,pairconstraints)
            else
                call gs%eval(map,gridcoord)
            endif
            call map%get_secondorder_forceconstant(p,fc,mem,-1)
            ! a phonon bandstructure
            write(*,*) 'FIXME GENERATE BANDSTRUCTURE',__FILE__
            ! call bs%generate(p,fc,timereversal=.true.,verbosity=0,npts=100,readpathfromfile=.true.,mpi_communicator=mw%comm)
            ! if ( mw%talk ) then
            !     ! dump this to file
            !     groupname='temperature_'//tochar(t)
            !     call h5gcreate_f(group_id,trim(groupname),subgroup_id,lo_status)
            !     call bs%write_to_hdf5( p, 'thz' , 'null', hdftag=subgroup_id )
            !     call h5gclose_f(subgroup_id,lo_status)
            ! endif
            if ( mw%talk ) call lo_progressbar(' ... dispersions at equilibrium',t,gr%nt,walltime()-t0)
        enddo

        ! and close the dispersion group
        if ( mw%talk ) call h5gclose_f(group_id,lo_status)
    end block dispersions
    endif

    ! Close files
    if ( mw%talk ) then
        call h5fclose_f(file_id, lo_status)
        call h5close_f(lo_status)
    endif

end subroutine
