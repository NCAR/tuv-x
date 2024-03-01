! Copyright (C) 2024 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
program test_heating_rates

  use musica_mpi,                      only : musica_mpi_init,                &
                                              musica_mpi_finalize
  use tuvx_heating_rates

  implicit none

  call musica_mpi_init( )
  call test_heating_rates_t( )
  call musica_mpi_finalize( )

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> @brief Test the heating rates
  subroutine test_heating_rates_t( )

    use musica_assert,                 only : assert
    use musica_config,                 only : config_t
    use musica_constants,              only : dk => musica_dk
    use musica_mpi,                    only : musica_mpi_bcast,               &
                                              musica_mpi_rank,                &
                                              MPI_COMM_WORLD
    use musica_string,                 only : string_t
    use tuvx_grid_warehouse,           only : grid_warehouse_t
    use tuvx_la_sr_bands,              only : la_sr_bands_t
    use tuvx_profile_warehouse,        only : profile_warehouse_t
    use tuvx_solver,                   only : radiation_field_t
    use tuvx_spherical_geometry,       only : spherical_geometry_t
    use tuvx_test_utils,               only : check_values

    type(heating_rates_t),      pointer :: heating_rates
    class(grid_warehouse_t),    pointer :: grids
    class(profile_warehouse_t), pointer :: profiles

    character(len=*), parameter :: Iam = "heating_rates_t tests"
    type(config_t) :: config, sub_config, reactions_config
    type(string_t), allocatable :: labels(:)
    character, allocatable :: buffer(:)
    integer :: pos, pack_size, i_height, i_wavelength
    integer, parameter :: comm = MPI_COMM_WORLD
    real(dk), parameter :: hc = 6.62608e-34_dk * 2.9979e8_dk / 1.e-9_dk
    real(dk) :: bde(6,2), actinic_flux(5,6), etfl(6)
    real(dk) :: wc(6) = (/ 425.0_dk, 475.0_dk, 525.0_dk, 575.0_dk, 625.0_dk,  &
                           675.0_dk /)
    type(radiation_field_t) :: radiation_field
    real(dk) :: calculated_rates(5,2), expected_rates(5,2)
    type(la_sr_bands_t) :: la_srb
    type(spherical_geometry_t) :: spherical_geometry

    call config%from_file( "test/data/heating_rates.json" )
    call config%get( "grids", sub_config, Iam )
    grids => grid_warehouse_t( sub_config )
    call config%get( "profiles", sub_config, Iam )
    profiles => profile_warehouse_t( sub_config, grids )

    if( musica_mpi_rank( comm ) == 0 ) then
      call config%get( "reactions", reactions_config, Iam )
      call sub_config%empty( )
      call sub_config%add( "reactions", reactions_config, Iam )
      heating_rates => heating_rates_t( sub_config, grids, profiles )
      pack_size = heating_rates%pack_size( comm )
      allocate( buffer( pack_size ) )
      pos = 0
      call heating_rates%mpi_pack( buffer, pos, comm )
      call assert( 534250649, pos <= pack_size )
    end if

    call musica_mpi_bcast( pack_size, comm )
    if( musica_mpi_rank( comm ) .ne. 0 ) allocate( buffer( pack_size ) )
    call musica_mpi_bcast( buffer, comm )

    if( musica_mpi_rank( comm ) .ne. 0 ) then
      pos = 0
      allocate( heating_rates )
      call heating_rates%mpi_unpack( buffer, pos, comm )
      call assert( 192483602, pos <= pack_size )
    end if
    deallocate( buffer )

    ! check labels
    labels = heating_rates%labels( )
    call assert( 152892147, size(labels) == 2 )
    call assert( 437272930, labels(1) == "jfoo" )
    call assert( 884640776, labels(2) == "jbaz" )

    ! check bond dissociation energies
    call assert( 613305591, size( heating_rates%heating_parameters_ ) == 2 )
    bde(:,1) = max( 0.0_dk, hc * ( 2.0_dk - wc(:) ) / ( 2.0_dk * wc(:) ) )
    call check_values( heating_rates%heating_parameters_(1)%energy_,          &
                       bde(:,1), 1.0e-4_dk )
    bde(:,2) = max( 0.0_dk, hc * ( 3000.0_dk - wc(:) ) / ( 3000.0_dk * wc(:) ) )
    call check_values( heating_rates%heating_parameters_(2)%energy_,          &
                       bde(:,2), 1.0e-4_dk )

    ! check calculated heating rates
    calculated_rates(:,:) = 0.0_dk
    allocate( radiation_field%fdr_(5,6), radiation_field%fdn_(5,6),           &
              radiation_field%fup_(5,6) )
    do i_wavelength = 1, 6
      etfl(i_wavelength) = 0.5_dk * ( 1.0e3_dk * 10.0_dk**i_wavelength +      &
                                      1.0e4_dk * 10.0_dk**i_wavelength )
    end do
    do i_height = 1, 5
      do i_wavelength = 1, 6
        radiation_field%fdr_( i_height, i_wavelength ) =                      &
            1.0_dk * i_height * i_wavelength
        radiation_field%fdn_( i_height, i_wavelength ) =                      &
            2.0_dk * i_height * i_wavelength
        radiation_field%fup_( i_height, i_wavelength ) =                      &
            3.0_dk * i_height * i_wavelength
        actinic_flux( i_height, i_wavelength ) =                              &
            6.0_dk * i_height * i_wavelength * etfl( i_wavelength )
        expected_rates( i_height, 1 ) =                                       &
          dot_product( actinic_flux(i_height,:), bde(:,1) * 12.3_dk * 0.75_dk )
        expected_rates( i_height, 2 ) = 1.1_dk *                              &
          dot_product( actinic_flux(i_height,:), bde(:,2) * 78.9_dk * 0.5_dk )
      end do
    end do
    call heating_rates%get( la_srb, spherical_geometry, grids, profiles,      &
                            radiation_field, calculated_rates )

    call check_values( calculated_rates, expected_rates, 1.0e-4_dk )

    deallocate( grids )
    deallocate( profiles )
    deallocate( heating_rates )

  end subroutine test_heating_rates_t

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program test_heating_rates
