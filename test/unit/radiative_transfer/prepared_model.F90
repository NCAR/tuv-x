! Copyright (C) 2026 University Corporation for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
! Exercise the complete prepared path against generic spectra/rates and the
! scalar transfer solver, including live profiles, grid changes and MPI rebuilds.
program test_prepared_model
  use musica_assert, only : assert
  use musica_config, only : config_t
  use musica_constants, only : dk => musica_dk
  use musica_mpi, only : musica_mpi_init, musica_mpi_finalize, MPI_COMM_WORLD
  use tuvx_core, only : core_t
  use tuvx_grid, only : grid_t
  use tuvx_grid_warehouse, only : grid_warehouse_t
  use tuvx_profile, only : profile_t
  use tuvx_profile_warehouse, only : profile_warehouse_t
  use tuvx_solver, only : radiation_field_t
  implicit none
    type profile_copy_t
      class(profile_t), pointer :: value => null()
    end type
  character(len=11), parameter :: names(*) = [character(len=11) :: 'temperature','air','O2','O3']
  call musica_mpi_init()
  call check_model('examples/tuv_5_4.json')
  call check_model('examples/ts1_tsmlt.json')
  call check_nonfinite_flux()
  call musica_mpi_finalize()
contains
  subroutine check_nonfinite_flux()
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_positive_inf, ieee_is_nan
    use, intrinsic :: ieee_exceptions, only : ieee_status_type, ieee_get_status, ieee_set_status, &
        ieee_set_halting_mode, ieee_invalid
    use tuvx_photolysis_rates, only : photolysis_rates_t
    use tuvx_la_sr_bands, only : la_sr_bands_t
    use tuvx_spherical_geometry, only : spherical_geometry_t
    type(config_t) :: config,part
    class(grid_warehouse_t), pointer :: grids
    class(profile_warehouse_t), pointer :: profiles
    type(photolysis_rates_t), pointer :: rates
    type(radiation_field_t), pointer :: field
    type(la_sr_bands_t) :: bands
    type(spherical_geometry_t) :: geometry
    type(ieee_status_type) :: saved_status
    real(dk) :: result(5,1)
    call config%from_file('test/data/solver_from_host.json')
    call config%get('grids',part,'nonfinite flux test')
    grids => grid_warehouse_t(part)
    call config%get('profiles',part,'nonfinite flux test')
    profiles => profile_warehouse_t(part,grids)
    part = '{"reactions":[{"name":"zero","cross section":{"type":"base",'// &
        '"data":{"default value":0},"apply O2 bands":false},'// &
        '"quantum yield":{"type":"base","constant value":1}}]}'
    rates => photolysis_rates_t(part,grids,profiles)
    field => radiation_field_t(5,6)
    field%fdr_ = 1
    field%fup_ = 0
    field%fdn_ = 0
    call rates%get(bands,geometry,grids,profiles,field,result,'')
    call assert(826130012,all(result == 0))
    ! A prepared zero term must not hide the generic 0*Inf -> NaN result.
    ! Save and restore flags/modes so this test does not change other tests.
    call ieee_get_status(saved_status)
    call ieee_set_halting_mode(ieee_invalid,.false.)
    field%fdr_(1,1) = ieee_value(0.0_dk,ieee_positive_inf)
    call rates%get(bands,geometry,grids,profiles,field,result,'')
    call assert(826130013,ieee_is_nan(result(1,1)))
    call assert(826130014,all(result(2:,1) == 0))
    call ieee_set_status(saved_status)
    deallocate(rates,field,grids,profiles)
  end subroutine

  subroutine check_model(path)
    character(*), intent(in) :: path
    type(config_t) :: config,reference_config,rt,solver,part
    type(config_t), allocatable :: items(:)
    class(core_t), pointer :: model,reference,unpacked
    type(profile_warehouse_t), pointer :: p,q
    type(grid_warehouse_t), pointer :: grids
    class(grid_t), pointer :: grid
    type(profile_copy_t) :: base(4)
    type(radiation_field_t) :: a,b
    real(dk), allocatable :: photo(:,:),heat(:,:),dose(:,:),ref_photo(:,:),ref_heat(:,:),ref_dose(:,:)
    real(dk), parameter :: angles(*) = [0.0_dk,45.0_dk,89.0_dk,90.0_dk,95.0_dk,110.0_dk,180.0_dk,52.0_dk]
    real(dk) :: scale,factor,offset
    integer :: nz,nw,np,nh,nd,column,i,j,pos
    character, allocatable :: buffer(:)
    logical :: found, diagnostic
    call config%from_file(path)
    call config%add('enable diagnostics',.false.,'prepared model test')
    call config%get('profiles',items,'prepared model test')
    do i = 1, size(items)
      call items(i)%get('enable diagnostics',diagnostic,'prepared model test',found=found)
      if( found ) call items(i)%add('enable diagnostics',.false.,'prepared model test')
    end do
    call config%add('profiles',items,'prepared model test')
    call config%get('radiative transfer',rt,'prepared model test')
    call rt%get('radiators',items,'prepared model test')
    do i = 1, size(items)
      call items(i)%add('enable diagnostics',.false.,'prepared model test')
    end do
    call rt%add('radiators',items,'prepared model test')
    call config%add('radiative transfer',rt,'prepared model test')
    call config%get('photolysis',part,'prepared model test')
    call part%add('enable diagnostics',.false.,'prepared model test')
    call config%add('photolysis',part,'prepared model test')
    call config%get('dose rates',part,'prepared model test',found=found)
    if( found ) then
      call part%add('enable diagnostics',.false.,'prepared model test')
      call config%add('dose rates',part,'prepared model test')
    end if
    reference_config = config
    call reference_config%add('column execution','generic','prepared model test')
    call rt%get('solver',solver,'prepared model test')
    call solver%add('execution','scalar','prepared model test')
    call rt%add('solver',solver,'prepared model test')
    call reference_config%add('radiative transfer',rt,'prepared model test')
    model => core_t(config)
    reference => core_t(reference_config)
    grid => model%get_grid('height','km')
    nz = grid%ncells_+1
    deallocate(grid)
    grid => model%get_grid('wavelength','nm')
    nw = grid%ncells_
    deallocate(grid)
    np = model%number_of_photolysis_reactions()
    nh = model%number_of_heating_rates()
    nd = model%number_of_dose_rates()
    allocate(photo(nz,np),heat(nz,nh),dose(nz,nd),ref_photo(nz,np),ref_heat(nz,nh),ref_dose(nz,nd))
    p => model%get_profile_warehouse()
    do i = 1, size(names)
      if( i == 1 ) then
        base(i)%value => p%get_profile(trim(names(i)),'K')
      else
        base(i)%value => p%get_profile(trim(names(i)),'molecule cm-3')
      end if
    end do
    do column = 1, 16
      ! Serialized models rebuild plans from their own objects, never pointers
      ! or cached values from the sending model. Exercise both execution modes.
      if( column == 3 .or. column == 6 ) then
        call model%set_prepared_execution(column /= 6)
        if( model%pack_size(MPI_COMM_WORLD) > 0 ) then
          allocate(buffer(model%pack_size(MPI_COMM_WORLD)))
          pos = 0
          call model%mpi_pack(buffer,pos,MPI_COMM_WORLD)
          allocate(unpacked)
          pos = 0
          call unpacked%mpi_unpack(buffer,pos,MPI_COMM_WORLD)
          deallocate(model,buffer)
          model => unpacked
          nullify(unpacked)
        end if
      end if
      if( column == 7 ) call model%set_prepared_execution(.true.)
      if( column == 10 ) call model%invalidate_prepared()
      if( column == 11 ) then
        grids => model%get_grid_warehouse()
        call shift_wavelengths(grids)
        grids => reference%get_grid_warehouse()
        call shift_wavelengths(grids)
      end if
      p => model%get_profile_warehouse()
      q => reference%get_profile_warehouse()
      factor = 1+.05_dk*sin(real(column,dk))
      offset = 10*sin(real(column,dk))
      call update_profiles(p,base,factor,offset,column == 14)
      call update_profiles(q,base,factor,offset,column == 14)
      scale = 1+.017_dk*sin(real(column,dk))
      if( column == 8 ) scale = 0
      call model%run(angles(mod(column-1,size(angles))+1),scale,photo,dose,heat)
      call reference%run(angles(mod(column-1,size(angles))+1),scale,ref_photo,ref_dose,ref_heat)
      a = model%get_radiation_field()
      b = reference%get_radiation_field()
      call compare(a%fdr_,b%fdr_,1.0_dk)
      call compare(a%fdn_,b%fdn_,1.0_dk)
      call compare(a%fup_,b%fup_,1.0_dk)
      call compare(a%edr_,b%edr_,1.0_dk)
      call compare(a%edn_,b%edn_,1.0_dk)
      call compare(a%eup_,b%eup_,1.0_dk)
      call compare(photo,ref_photo,tiny(1.0_dk))
      call compare(heat,ref_heat,tiny(1.0_dk))
      call compare(dose,ref_dose,tiny(1.0_dk))
    end do
    do i = 1, size(base)
      deallocate(base(i)%value)
    end do
    deallocate(model,reference)
    print *, 'Verified changing prepared model:',path
  end subroutine
    subroutine update_profiles(profiles,original,factor,offset,hot)
      type(profile_warehouse_t), intent(inout) :: profiles
      type(profile_copy_t), intent(in) :: original(:)
      real(dk), intent(in) :: factor,offset
      logical, intent(in) :: hot
      integer :: i,j
      do i = 1, size(original)
        do j = 1, size(profiles%profiles_)
          associate(p => profiles%profiles_(j)%val_,base => original(i)%value)
            if( p%handle_ /= names(i) ) cycle
            if( i == 1 ) then
              p%edge_val_ = base%edge_val_+offset
              p%mid_val_ = base%mid_val_+offset
              if( hot ) p%edge_val_(1) = 11000.0_dk
            else
              p%edge_val_ = base%edge_val_*factor
              p%mid_val_ = base%mid_val_*factor
              if( allocated(base%layer_dens_) ) p%layer_dens_ = base%layer_dens_*factor
              if( allocated(base%exo_layer_dens_) ) p%exo_layer_dens_ = base%exo_layer_dens_*factor
              if( allocated(base%burden_dens_) ) p%burden_dens_ = base%burden_dens_*factor
            end if
          end associate
        end do
      end do
    end subroutine
  subroutine shift_wavelengths(grids)
    type(grid_warehouse_t), intent(inout) :: grids
    integer :: i
    do i = 1, size(grids%grids_)
      associate(g => grids%grids_(i)%val_)
        if( g%handle_ /= 'wavelength' ) cycle
        g%mid_ = g%mid_+1.e-5_dk
        g%edge_ = g%edge_+1.e-5_dk
      end associate
    end do
  end subroutine
  subroutine compare(a,b,floor)
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    real(dk), intent(in) :: a(:,:),b(:,:),floor
    real(dk), allocatable :: bound(:,:)
    call assert(826130010,all(shape(a) == shape(b)))
    call assert(826130011,all(ieee_is_finite(a)) .and. all(ieee_is_finite(b)))
    if( size(a) == 0 ) return
    bound = spread(max(maxval(abs(a),dim=1),maxval(abs(b),dim=1),floor),1,size(a,1))
    if( any(abs(a-b) > 64*epsilon(1.0_dk)*bound) ) then
      print *, 'Model mismatch, maximum normalized eps:',maxval(abs(a-b)/bound/epsilon(1.0_dk))
      stop 1
    end if
  end subroutine
end program
