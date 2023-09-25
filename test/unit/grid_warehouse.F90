! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
program test_grid_warehouse

  use musica_mpi,                      only : musica_mpi_init,                &
                                              musica_mpi_finalize
  use tuvx_grid_warehouse

  implicit none

  call musica_mpi_init( )
  call test_grid_warehouse_t( )
  call musica_mpi_finalize( )

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine test_grid_warehouse_t( )

    use musica_assert,                 only : almost_equal, assert
    use musica_constants,              only : dk => musica_dk
    use musica_config,                 only : config_t
    use musica_mpi
    use musica_string,                 only : string_t
    use tuvx_grid,                     only : grid_t
    use tuvx_grid_from_host,           only : grid_from_host_t, grid_updater_t
#ifdef MUSICA_USE_MPI
    use mpi,                           only : MPI_COMM_WORLD
#endif

    class(grid_warehouse_t), pointer :: grids

    character(len=*), parameter :: my_name = "grid warehouse tests"
    type(config_t) :: config
    class(grid_t), pointer :: grid
    class(grid_warehouse_t), pointer :: host_grids
    class(grid_from_host_t), pointer :: host_grid
    type(grid_updater_t) :: grid_updater
    character, allocatable :: buffer(:)
    type(string_t) :: grid_name, units
    integer :: pos, pack_size
    integer, parameter :: comm = MPI_COMM_WORLD
    real(dk) :: grid_mids(3), grid_edges(4)
    logical :: found
    type(grid_warehouse_ptr) :: ptr

    host_grid => grid_from_host_t( "foo", "bars", 3 )
    host_grids => grid_warehouse_t( )
    call host_grids%add( host_grid )

    ! load test grids
    if( musica_mpi_rank( comm ) == 0 ) then
      call config%from_file( "test/data/grid.simple.config.json" )
      grids => grid_warehouse_t( config )
      call grids%add( host_grids )
      ptr = grids%get_ptr( "time", "hours" )
      pack_size = grids%pack_size( comm ) + ptr%pack_size( comm )
      allocate( buffer( pack_size ) )
      pos = 0
      call grids%mpi_pack( buffer, pos, comm )
      call ptr%mpi_pacK(   buffer, pos, comm )
      call assert( 441119696, pos <= pack_size )
    end if
    deallocate( host_grids )

    call musica_mpi_bcast( pack_size, comm )
    if( musica_mpi_rank( comm ) .ne. 0 ) allocate( buffer( pack_size ) )
    call musica_mpi_bcast( buffer, comm )

    if( musica_mpi_rank( comm ) .ne. 0 ) then
      pos = 0
      allocate( grids )
      call grids%mpi_unpack( buffer, pos, comm )
      call ptr%mpi_unpack(   buffer, pos, comm )
      call assert( 211670994, pos <= pack_size )
    end if

    grid_name = "time"
    units = "hours"
    grid => grids%get_grid( grid_name, units )
    call assert( 355697132, allocated( grid%mid_   ) )
    call assert( 115719816, allocated( grid%edge_  ) )
    call assert( 285876720, allocated( grid%delta_ ) )
    call assert( 571710742, almost_equal( grid%edge_(1), 5.0_dk  ) )
    call assert( 237400659, almost_equal( grid%edge_(2), 10.0_dk ) )
    call assert( 684768505, almost_equal( grid%edge_(3), 15.0_dk ) )
    deallocate( grid )

    grid => grids%get_grid( ptr )
    call assert( 700665311, allocated( grid%mid_   ) )
    call assert( 530508407, allocated( grid%edge_  ) )
    call assert( 977876253, allocated( grid%delta_ ) )
    call assert( 525244100, almost_equal( grid%edge_(1), 5.0_dk  ) )
    call assert( 420095596, almost_equal( grid%edge_(2), 10.0_dk ) )
    call assert( 867463442, almost_equal( grid%edge_(3), 15.0_dk ) )
    deallocate( grid )

    grid_name = "foo"
    units = "bars"
    grid_updater = grids%get_updater( host_grid, found )
    call assert( 429189037, found )
    grid_mids = (/ 12.5_dk, 20.0_dk, 30.0_dk /)
    grid_edges = (/ 10.0_dk, 15.0_dk, 25.0_dk, 35.0_dk /)
    call grid_updater%update( edges = grid_edges, mid_points = grid_mids )
    grid => grids%get_grid( grid_name, units )
    call assert( 522094162, allocated( grid%mid_   ) )
    call assert( 351937258, allocated( grid%edge_  ) )
    call assert( 181780354, allocated( grid%delta_ ) )
    call assert( 629148200, almost_equal( grid%edge_(1),  10.0_dk ) )
    call assert( 241524447, almost_equal( grid%edge_(2),  15.0_dk ) )
    call assert( 688892293, almost_equal( grid%edge_(3),  25.0_dk ) )
    call assert( 236260140, almost_equal( grid%edge_(4),  35.0_dk ) )
    call assert( 683627986, almost_equal( grid%mid_(1),   12.5_dk ) )
    call assert( 230995833, almost_equal( grid%mid_(2),   20.0_dk ) )
    call assert( 125847329, almost_equal( grid%mid_(3),   30.0_dk ) )
    call assert( 573215175, almost_equal( grid%delta_(1),  5.0_dk ) )
    call assert( 968008769, almost_equal( grid%delta_(2), 10.0_dk ) )
    call assert( 180327115, almost_equal( grid%delta_(3), 10.0_dk ) )
    deallocate( grid )
    deallocate( host_grid )

    grid_name = "not there"
    units     = "none"
    host_grid => grid_from_host_t( grid_name, units, 27 )
    grid_updater = grids%get_updater( host_grid, found )
    call assert( 190403434, .not. found )
    deallocate( host_grid )

    deallocate( grids )

  end subroutine test_grid_warehouse_t

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program test_grid_warehouse
