! Copyright (C) 2024 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
program test_quantum_yield_h2so4_mills

  use musica_mpi
  use tuvx_quantum_yield_h2so4_mills

  implicit none

  call musica_mpi_init( )
  call test_quantum_yield_h2so4_mills_t( )
  call musica_mpi_finalize( )

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Test the temperature and pressure dependent H2SO4 quantum yield
  !! calculations against previously generated results from an older version
  !! of TUV
  subroutine test_quantum_yield_h2so4_mills_t( )

    use musica_assert,                 only : assert, die, almost_equal
    use musica_config,                 only : config_t
    use musica_constants,              only : dk => musica_dk
    use musica_io,                     only : io_t
    use musica_io_netcdf,              only : io_netcdf_t
    use musica_string,                 only : string_t
    use tuvx_constants,                only : gas_constant, Avogadro, pi
    use tuvx_cross_section,            only : cross_section_t
    use tuvx_grid,                     only : grid_t
    use tuvx_grid_from_host,           only : grid_from_host_t
    use tuvx_grid_warehouse,           only : grid_warehouse_t
    use tuvx_profile,                  only : profile_t
    use tuvx_profile_from_host,        only : profile_from_host_t,            &
                                              profile_updater_t
    use tuvx_profile_warehouse,        only : profile_warehouse_t
    use tuvx_quantum_yield,            only : quantum_yield_t
    use tuvx_quantum_yield_factory
    use tuvx_test_utils,               only : check_values

    character(len=*), parameter :: my_name = "h2so4 quantum yield test"
    class(grid_warehouse_t),    pointer :: grids
    class(profile_warehouse_t), pointer :: profiles
    class(quantum_yield_t),     pointer :: quantum_yield
    class(cross_section_t),     pointer :: cross_section
    type(config_t) :: config, qy_config, cs_config, grids_config,             &
                      profiles_config

    class(grid_from_host_t), pointer :: heights
    class(profile_from_host_t), pointer :: temperature, air
    type(profile_updater_t) :: temperature_updater, air_updater
    class(grid_t), pointer :: wavelength_grid
    class(profile_t), pointer :: temperature_profile, air_profile
    character, allocatable :: buffer(:)
    real(kind=dk), allocatable :: update_temperatures(:), update_air(:)
    real(kind=dk), allocatable :: file_temperatures(:), file_pressures(:),    &
                                  file_photo_rates(:,:,:),                    &
                                  file_wavelengths(:),                        &
                                  quantum_yields(:,:), cross_sections(:,:)
    class(io_t), pointer :: file
    type(string_t) :: type_name, file_name, var_name
    integer :: pos, pack_size, n_bins, i_temp, i_pres, i_wl, i_height, n_wl,  &
               i_file_offset
    integer, parameter :: comm = MPI_COMM_WORLD

    file_name = "test/data/quantum_yields/jh2so4.nc"
    file => io_netcdf_t( file_name )
    var_name = "temperature"
    call file%read( var_name, file_temperatures, my_name )
    var_name = "pressure"
    call file%read( var_name, file_pressures, my_name )
    file_pressures(:) = file_pressures(:) * 100.0_dk ! convert hPa to Pa
    var_name = "jh2so4"
    call file%read( var_name, file_photo_rates, my_name )
    var_name = "wavelength"
    call file%read( var_name, file_wavelengths, my_name )
    deallocate( file )

    call config%from_file( "test/data/quantum_yields/h2so4_mills.config.json" )
    call config%get( "grids",         grids_config,    my_name )
    call config%get( "quantum yield", qy_config,       my_name )
    call config%get( "cross section", cs_config,       my_name )

    n_bins = size( file_temperatures ) * size( file_pressures ) - 1
    heights => grid_from_host_t( "height", "km", n_bins )
    temperature => profile_from_host_t( "temperature", "K", n_bins )
    air => profile_from_host_t( "air", "molecule cm-3", n_bins )

    grids => grid_warehouse_t( grids_config )
    call grids%add( heights )
    call profiles_config%empty( )
    profiles => profile_warehouse_t( profiles_config, grids )
    call profiles%add( temperature )
    call profiles%add( air         )
    temperature_updater = profiles%get_updater( temperature )
    air_updater = profiles%get_updater( air )

    cross_section => cross_section_t( cs_config, grids, profiles )

    if( musica_mpi_rank( comm ) == 0 ) then
      quantum_yield => quantum_yield_builder( qy_config, grids, profiles )
      type_name = quantum_yield_type_name( quantum_yield )
      pack_size = type_name%pack_size( comm ) + quantum_yield%pack_size( comm )
      allocate( buffer( pack_size ) )
      pos = 0
      call type_name%mpi_pack(     buffer, pos, comm )
      call quantum_yield%mpi_pack( buffer, pos, comm )
      call assert( 837477723, pos <= pack_size )
    end if

    call musica_mpi_bcast( pack_size, comm )
    if( musica_mpi_rank( comm ) .ne. 0 ) allocate( buffer( pack_size ) )
    call musica_mpi_bcast( buffer, comm )

    if( musica_mpi_rank( comm ) .ne. 0 ) then
      pos = 0
      call type_name%mpi_unpack( buffer, pos, comm )
      quantum_yield => quantum_yield_allocate( type_name )
      call quantum_yield%mpi_unpack( buffer, pos, comm )
      call assert( 599405941, pos <= pack_size )
    end if
    deallocate( buffer )

    i_height = 0
    allocate( update_temperatures( n_bins + 1 ) )
    allocate( update_air(          n_bins + 1 ) )
    do i_temp = 1, size( file_temperatures )
      do i_pres = 1, size( file_pressures )
        i_height = i_height + 1
        update_temperatures( i_height ) = file_temperatures( i_temp )
        update_air( i_height ) = file_pressures( i_pres ) / gas_constant /    &
                                 update_temperatures( i_height ) *            &
                                 Avogadro * 1.0e-6_dk ! convert mol m-3 to molec cm-3
      end do
    end do
    call temperature_updater%update( edge_values = update_temperatures )
    call air_updater%update(         edge_values = update_air          )

    quantum_yields = quantum_yield%calculate( grids, profiles )
    cross_sections = cross_section%calculate( grids, profiles )

    wavelength_grid     => grids%get_grid( "wavelength", "nm" )
    temperature_profile => profiles%get_profile( "temperature", "K" )
    air_profile         => profiles%get_profile( "air", "molecule cm-3" )

    n_wl = size( quantum_yields, 2 )
    i_file_offset = 1
    do i_wl = 1, n_wl
      if ( almost_equal( file_wavelengths( 1 ),                               &
                         wavelength_grid%mid_( i_wl ) ) ) then
        exit
      end if
      i_file_offset = i_file_offset + 1
    end do
    call assert( 462064586, size( cross_sections, 2 ) .eq. n_wl )
    call assert( 234069123, size( file_photo_rates, 1 ) + i_file_offset - 1   &
                            .eq. n_wl )

    i_height = 0
    do i_temp = 1, size( file_temperatures )
      do i_pres = 1, size( file_pressures )
        i_height = i_height + 1
        do i_wl = 1, i_file_offset - 1
            call assert( 897976065,                                           &
                         almost_equal( quantum_yields( i_height, i_wl ),      &
                                       1.0_dk ) )
          if( i_wl .eq. 2 ) then
            call assert( 126374455,                                           &
                         almost_equal( cross_sections( i_height, i_wl ),      &
                                       6.3e-17_dk ) )
          else
            call assert( 291267052,                                           &
                         almost_equal( cross_sections( i_height, i_wl ),      &
                                       0.0_dk ) )
          end if
        end do
        do i_wl = i_file_offset, n_wl
          ! the top pressure level has different logic, but we're putting
          ! all pressure/temperature combos in one profile for this test,
          ! so skip the lowest pressure util we're on the last profile element
          if( i_pres .eq. size( file_pressures ) .and.                        &
              i_temp < size( file_temperatures ) ) cycle
          call assert( 342289277,                                             &
                       almost_equal( quantum_yields( i_height, i_wl ) *       &
                                     cross_sections( i_height, i_wl ),        &
                file_photo_rates( i_wl - i_file_offset + 1, i_temp, i_pres ), &
                relative_tolerance = 1.0e-3_dk ))
        end do
      end do
    end do

    ! clean up
    deallocate( grids               )
    deallocate( profiles            )
    deallocate( quantum_yield       )
    deallocate( cross_section       )
    deallocate( wavelength_grid     )
    deallocate( temperature_profile )
    deallocate( air_profile         )
    deallocate( heights             )
    deallocate( temperature         )
    deallocate( air                 )

  end subroutine test_quantum_yield_h2so4_mills_t

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program test_quantum_yield_h2so4_mills