! Copyright (C) 2021 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0

program test_profile_warehouse
  ! Test module for the profile warehouse

  use musica_mpi,                      only : musica_mpi_init,                &
                                              musica_mpi_finalize

  implicit none

  call musica_mpi_init( )
  call profile_warehouse_test( )
  call musica_mpi_finalize( )

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine profile_warehouse_test()
    use musica_assert,              only : assert, almost_equal
    use musica_config,              only : config_t
    use musica_constants,           only : dk => musica_dk
    use musica_mpi
    use musica_string,              only : string_t
    use tuvx_grid_warehouse,        only : grid_warehouse_t
    use tuvx_profile,               only : profile_t
    use tuvx_profile_from_host,     only : profile_from_host_t,               &
                                           profile_updater_t
    use tuvx_profile_warehouse,     only : profile_warehouse_t,               &
                                           profile_warehouse_ptr

    character(len=*), parameter :: grid_config = 'test/data/grid.test.config.json'
    character(len=*), parameter :: profile_config = 'test/data/profile.temperature.config.json'
    type(config_t) :: grid_tst_config
    type(config_t) :: profile_tst_config
    class(grid_warehouse_t), pointer :: grid_warehouse => null()
    class(profile_warehouse_t), pointer :: profile_warehouse => null()
    class(profile_t),           pointer :: profile_ptr => null()
    class(profile_warehouse_t), pointer :: host_profiles
    class(profile_from_host_t), pointer :: host_profile
    type(profile_updater_t) :: my_updater
    type(string_t) :: profile_name, units
    logical :: found
    character, allocatable :: buffer(:)
    integer :: pos, pack_size
    integer, parameter :: comm = MPI_COMM_WORLD
    type(profile_warehouse_ptr) :: ptr

    host_profile => profile_from_host_t( "foo", "bars", 3 )
    host_profiles => profile_warehouse_t( )
    call host_profiles%add( host_profile )

    call grid_tst_config%from_file( grid_config )
    grid_warehouse => grid_warehouse_t( grid_tst_config )

    if( musica_mpi_rank( comm ) == 0 ) then
      call profile_tst_config%from_file( profile_config )
      profile_warehouse =>                                                    &
          profile_warehouse_t( profile_tst_config, grid_warehouse )
      call profile_warehouse%add( host_profiles )
      ptr = profile_warehouse%get_ptr( "temperature", "K" )
      pack_size = profile_warehouse%pack_size( comm ) + ptr%pack_size( comm )
      allocate( buffer( pack_size ) )
      pos = 0
      call profile_warehouse%mpi_pack( buffer, pos, comm )
      call ptr%mpi_pack( buffer, pos, comm )
      call assert( 893633789, pos <= pack_size )
    end if
    deallocate( host_profiles )

    call musica_mpi_bcast( pack_size, comm )
    if( musica_mpi_rank( comm ) .ne. 0 ) allocate( buffer( pack_size ) )
    call musica_mpi_bcast( buffer, comm )

    if( musica_mpi_rank( comm ) .ne. 0 ) then
      pos = 0
      allocate( profile_warehouse )
      call profile_warehouse%mpi_unpack( buffer, pos, comm )
      call ptr%mpi_unpack( buffer, pos, comm )
      call assert( 436189624, pos <= pack_size )
    end if
    deallocate( buffer )

    profile_ptr => profile_warehouse%get_profile( "temperature", "K" )
    call assert( 418832741, associated(profile_ptr) )
    call assert( 811107812, profile_ptr%handle_ == "temperature" )
    call assert( 241345302, profile_ptr%units( ) == "K" )
    deallocate( profile_ptr )

    profile_ptr => profile_warehouse%get_profile( ptr )
    call assert( 418832741, associated(profile_ptr) )
    call assert( 811107812, profile_ptr%handle_ == "temperature" )
    call assert( 241345302, profile_ptr%units( ) == "K" )
    deallocate( profile_ptr )

    profile_name = "foo"
    units = "bars"
    my_updater = profile_warehouse%get_updater( host_profile, found )
    call assert( 289872557, found )
    call my_updater%update( mid_point_values = (/ 1.0_dk, 12.3_dk, 32.4_dk /),&
                         edge_values = (/ 0.5_dk, 9.8_dk, 15.4_dk, 45.0_dk /),&
                         layer_densities = (/ 94.3_dk, 0.52_dk, -12.3_dk /) )
    profile_ptr => profile_warehouse%get_profile( profile_name, units )
    call assert( 637694837, profile_ptr%size( ) == 3 )
    call assert( 233174517, profile_ptr%mid_val_(1) ==  1.0_dk )
    call assert( 345492862, profile_ptr%mid_val_(2) == 12.3_dk )
    call assert( 510385459, profile_ptr%mid_val_(3) == 32.4_dk )
    call assert( 957753305, profile_ptr%edge_val_(1) ==  0.5_dk )
    call assert( 505121152, profile_ptr%edge_val_(2) ==  9.8_dk )
    call assert( 399972648, profile_ptr%edge_val_(3) == 15.4_dk )
    call assert( 847340494, profile_ptr%edge_val_(4) == 45.0_dk )
    call assert( 394708341,                                                   &
                 almost_equal( profile_ptr%delta_val_(1),  9.8_dk -  0.5_dk ) )
    call assert( 559600938,                                                   &
                 almost_equal( profile_ptr%delta_val_(2), 15.4_dk -  9.8_dk ) )
    call assert( 106968785,                                                   &
                 almost_equal( profile_ptr%delta_val_(3), 45.0_dk - 15.4_dk ) )
    call assert( 619345031, profile_ptr%layer_dens_(1) == 94.3_dk )
    call assert( 166712878, profile_ptr%layer_dens_(2) == 0.52_dk )
    call assert( 331605475, profile_ptr%layer_dens_(3) == -12.3_dk )
    deallocate( profile_ptr )
    deallocate( host_profile )

    profile_name = "not there"
    units = "none"
    host_profile => profile_from_host_t( profile_name, units, 27 )
    my_updater = profile_warehouse%get_updater( host_profile, found )
    call assert( 783097109, .not. found )
    deallocate( host_profile )

    deallocate( grid_warehouse )
    deallocate( profile_warehouse )

  end subroutine profile_warehouse_test

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program test_profile_warehouse
