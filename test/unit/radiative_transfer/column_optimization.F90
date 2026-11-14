! Copyright (C) 2026 University Corporation for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
! Compare the blocked solver against the preserved scalar algorithm, across
! complete and partial blocks, thin/thick columns, twilight and night.
program test_column_optimization
  use musica_assert, only : assert
  use musica_config, only : config_t
  use musica_constants, only : dk => musica_dk
  use musica_mpi, only : musica_mpi_init, musica_mpi_finalize, MPI_COMM_WORLD
  use tuvx_grid_from_host, only : grid_from_host_t
  use tuvx_grid_warehouse, only : grid_warehouse_t
  use tuvx_profile_from_host, only : profile_from_host_t
  use tuvx_profile_warehouse, only : profile_warehouse_t
  use tuvx_radiator_warehouse, only : radiator_warehouse_t
  use tuvx_solver, only : radiation_field_t
  use tuvx_solver_delta_eddington, only : solver_delta_eddington_t
  use tuvx_spherical_geometry, only : spherical_geometry_t
  use tuvx_temperature_parameterization, only : temperature_parameterization_t
  use tuvx_temperature_parameterization_taylor_series, only : temperature_parameterization_taylor_series_t
  use tuvx_cross_section, only : cross_section_t
  use tuvx_quantum_yield, only : quantum_yield_t
  implicit none
  integer :: n, w, k
  integer, parameter :: layers(*) = [1,2,7,81]
  integer, parameter :: bins(*) = [1,15,16,17,31,32,33]
  real(dk), parameter :: angles(*) = [0.0_dk,45.0_dk,89.0_dk,90.0_dk,95.0_dk,110.0_dk,180.0_dk]

  call musica_mpi_init()
  call test_spectra()
  call test_parameterization()
  do n = 1, size(layers)
    do w = 1, size(bins)
      do k = 1, 5
        call test_solver(layers(n),bins(w),k)
      end do
    end do
  end do
  call musica_mpi_finalize()
contains
  subroutine compare(a,b,absolute_scale)
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    real(dk), intent(in) :: a(:,:), b(:,:)
    real(dk), optional, intent(in) :: absolute_scale
    real(dk) :: floor
    floor = tiny(1.0_dk)
    if( present(absolute_scale) ) floor = absolute_scale
    call assert(581964321, all(shape(a) == shape(b)))
    call assert(581964322, all(ieee_is_finite(a)) .and. all(ieee_is_finite(b)))
    ! Radiation is normalized to unit incident flux; cancellation at a
    ! zero-albedo boundary needs an absolute roundoff floor on that scale.
    ! The bound includes propagation through up to 162 elimination rows.
    call assert(581964323, all(abs(a-b) <= 1024*epsilon(1.0_dk)*max(abs(a),abs(b),floor)))
  end subroutine compare

  subroutine test_solver(n,nw,regime)
    integer, intent(in) :: n, nw, regime
    type(grid_warehouse_t) :: grids
    type(profile_warehouse_t) :: profiles
    type(radiator_warehouse_t) :: radiators
    type(grid_from_host_t), pointer :: height, wavelength
    type(profile_from_host_t), pointer :: albedo
    type(solver_delta_eddington_t), pointer :: solver, scalar_solver
    type(spherical_geometry_t), pointer :: geometry
    type(radiation_field_t), pointer :: field, reference
    type(config_t) :: config
    integer :: i, j, r, angle, surface, pos
    character, allocatable :: buffer(:)
    real(dk), parameter :: albedos(*) = [0.0_dk,0.1_dk,1.0_dk]
    height => grid_from_host_t('height','km',n)
    wavelength => grid_from_host_t('wavelength','nm',nw)
    height%edge_ = [(120.0_dk*real(i,dk)/n,i=0,n)]
    height%mid_ = (height%edge_(:n)+height%edge_(2:))*0.5_dk
    height%delta_ = height%edge_(2:)-height%edge_(:n)
    allocate(grids%grids_(2))
    grids%grids_(1)%val_ => height
    grids%grids_(2)%val_ => wavelength
    albedo => profile_from_host_t('surface albedo','none',nw)
    allocate(profiles%profiles_(1))
    profiles%profiles_(1)%val_ => albedo
    allocate(radiators%radiators_(2))
    do r = 1, 2
      allocate(radiators%radiators_(r)%val_)
      associate(state => radiators%radiators_(r)%val_%state_)
        allocate(state%layer_OD_(n,nw),state%layer_SSA_(n,nw),state%layer_G_(n,nw,1))
        do j = 1, nw
          do i = 1, n
            state%layer_OD_(i,j) = real(mod(13*i+7*j+3*r,101)+1,dk)*0.003_dk
            state%layer_SSA_(i,j) = real(mod(17*i+11*j+r,97),dk)/97.0_dk
            state%layer_G_(i,j,1) = real(mod(5*i+3*j+r,19)-9,dk)/10.0_dk
          end do
        end do
        select case(regime)
        case(1)
          state%layer_OD_ = 0.0_dk
        case(2)
          state%layer_OD_ = state%layer_OD_*1.e-10_dk
        case(3)
          state%layer_OD_ = state%layer_OD_*100.0_dk
        case(4)
          state%layer_SSA_ = 1.0_dk
          state%layer_G_ = 0.0_dk
        end select
      end associate
    end do
    call config%add('type','delta eddington','column optimization test')
    solver => solver_delta_eddington_t(config,grids,profiles)
    call config%add('execution','scalar','column optimization test')
    scalar_solver => solver_delta_eddington_t(config,grids,profiles)
    ! Exercise both execution settings after MPI serialization too.
    if( solver%pack_size(MPI_COMM_WORLD) > 0 ) then
      allocate(buffer(solver%pack_size(MPI_COMM_WORLD)))
      pos = 0
      call solver%mpi_pack(buffer,pos,MPI_COMM_WORLD)
      pos = 0
      call solver%mpi_unpack(buffer,pos,MPI_COMM_WORLD)
      deallocate(buffer)
      allocate(buffer(scalar_solver%pack_size(MPI_COMM_WORLD)))
      pos = 0
      call scalar_solver%mpi_pack(buffer,pos,MPI_COMM_WORLD)
      pos = 0
      call scalar_solver%mpi_unpack(buffer,pos,MPI_COMM_WORLD)
      deallocate(buffer)
    end if
    geometry => spherical_geometry_t(grids)
    do surface = 1, size(albedos)
      albedo%mid_val_ = albedos(surface)
      do angle = 1, size(angles)
        call geometry%set_parameters(angles(angle),grids)
        field => solver%update_radiation_field(angles(angle),n,geometry,grids,profiles,radiators)
        reference => scalar_solver%update_radiation_field(angles(angle),n,geometry,grids,profiles,radiators)
        call compare(field%fdr_,reference%fdr_,1.0_dk)
        call compare(field%fdn_,reference%fdn_,1.0_dk)
        call compare(field%fup_,reference%fup_,1.0_dk)
        call compare(field%edr_,reference%edr_,1.0_dk)
        call compare(field%edn_,reference%edn_,1.0_dk)
        call compare(field%eup_,reference%eup_,1.0_dk)
        deallocate(field,reference)
      end do
    end do
    deallocate(solver,scalar_solver,geometry)
  end subroutine test_solver

  subroutine test_parameterization()
    type(temperature_parameterization_t) :: param
    type(temperature_parameterization_taylor_series_t) :: taylor
    type(config_t) :: config
    real(dk), parameter :: temperatures(*) = [180.0_dk,200.0_dk,210.0_dk,250.0_dk,295.2_dk,300.0_dk,320.0_dk]
    real(dk) :: wavelengths(5), expected(5,size(temperatures)), actual(5,size(temperatures))
    integer :: i, base, inverted
    wavelengths = [302.0_dk,304.0_dk,306.0_dk,308.0_dk,310.0_dk]
    param%AA_ = [-20.0_dk,0.01_dk,-0.0001_dk]
    param%BB_ = [0.001_dk,-0.00001_dk,0.0_dk]
    param%lp_ = [0.0_dk,1.0_dk,2.0_dk]
    param%base_temperature_ = 273.0_dk
    param%base_wavelength_ = 300.0_dk
    param%min_wavelength_index_ = 2
    param%max_wavelength_index_ = 4
    allocate(param%ranges_(2))
    param%ranges_(1)%max_temperature_ = 210.0_dk
    param%ranges_(1)%is_fixed_ = .true.
    param%ranges_(1)%fixed_temperature_ = 210.0_dk
    param%ranges_(2)%min_temperature_ = 210.0_dk
    do base = 1, 2
      param%is_base_10_ = base == 1
      do inverted = 1, 2
        param%is_temperature_inverted_ = inverted == 1
        expected = 1.e-23_dk
        actual = expected
        do i = 1, size(temperatures)
          call param%calculate(temperatures(i),wavelengths,expected(:,i))
        end do
        call param%calculate_column(temperatures,wavelengths,actual)
        call compare(actual,expected)
      end do
    end do
    ! The column wrapper must still dispatch to derived parameterizations.
    call config%from_file('test/data/cross_sections/util/taylor.config.json')
    taylor = temperature_parameterization_taylor_series_t(config)
    wavelengths = taylor%merge_wavelength_grids([real(dk) ::])
    expected = 0.0_dk
    actual = expected
    do i = 1, size(temperatures)
      call taylor%calculate(temperatures(i),wavelengths,expected(:,i))
    end do
    call taylor%calculate_column(temperatures,wavelengths,actual)
    call compare(actual,expected)
  end subroutine test_parameterization

  subroutine test_spectra()
    type(cross_section_t) :: cs
    type(quantum_yield_t) :: qy
    type(grid_warehouse_t) :: grids
    type(profile_warehouse_t) :: profiles
    type(grid_from_host_t), pointer :: height
    real(dk), allocatable :: full(:,:), compact(:,:)
    height => grid_from_host_t('height','km',7)
    allocate(grids%grids_(1))
    grids%grids_(1)%val_ => height
    cs%height_grid_ = grids%get_ptr('height','km')
    qy%height_grid_ = cs%height_grid_
    allocate(cs%cross_section_parms(1),qy%quantum_yield_parms(1))
    allocate(cs%cross_section_parms(1)%array(4,1),qy%quantum_yield_parms(1)%array(4,1))
    cs%cross_section_parms(1)%array(:,1) = [0.0_dk,1.e-20_dk,2.e-20_dk,0.0_dk]
    qy%quantum_yield_parms(1)%array(:,1) = [0.0_dk,0.3_dk,0.4_dk,1.0_dk]
    allocate(cs%overrides_(1),qy%overrides_(1))
    cs%overrides_(1)%min_wavelength_index_ = 2
    cs%overrides_(1)%max_wavelength_index_ = 3
    cs%overrides_(1)%value_ = 3.e-20_dk
    qy%overrides_(1)%min_wavelength_index_ = 1
    qy%overrides_(1)%max_wavelength_index_ = 2
    qy%overrides_(1)%value_ = 0.5_dk
    full = cs%calculate(grids,profiles)
    compact = cs%calculate_column(grids,profiles)
    call assert(581964324,size(compact,1) == 1)
    call compare(full,spread(compact(1,:),1,8))
    full = qy%calculate(grids,profiles)
    compact = qy%calculate_column(grids,profiles)
    call assert(581964325,size(compact,1) == 1)
    call compare(full,spread(compact(1,:),1,8))
  end subroutine test_spectra
end program test_column_optimization
