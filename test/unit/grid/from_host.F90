! Copyright (C) 2021 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
program test_grid_from_host

  use musica_mpi,                      only : musica_mpi_init,                &
                                              musica_mpi_finalize

  implicit none

  call musica_mpi_init( )
  call test_grid_from_host_t( )
  call musica_mpi_finalize( )

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine test_grid_from_host_t( )

    use musica_assert,                 only : assert, almost_equal, die
    use musica_constants,              only : dk => musica_dk
    use musica_string,                 only : string_t
    use musica_mpi
    use tuvx_grid,                     only : grid_t
    use tuvx_grid_from_host,           only : grid_from_host_t, grid_updater_t
    use tuvx_grid_factory,             only : grid_type_name, grid_allocate
    use tuvx_test_utils,               only : check_values

    class(grid_t), pointer :: my_grid
    type(grid_updater_t) :: my_updater
    character, allocatable :: buffer(:)
    integer :: pos, pack_size
    type(string_t) :: type_name
    integer, parameter :: comm = MPI_COMM_WORLD
    real(kind=dk) :: mids(3), edges(4), deltas(3)
    real(kind=dk) :: tol = 1.0e-10_dk

    if( musica_mpi_rank( comm ) == 0 ) then
      my_grid => grid_from_host_t( "foo", "bars", 3 )
      type_name = grid_type_name( my_grid )
      pack_size = type_name%pack_size( comm ) + my_grid%pack_size( comm )
      allocate( buffer( pack_size ) )
      pos = 0
      call type_name%mpi_pack( buffer, pos, comm )
      call my_grid%mpi_pack(   buffer, pos, comm )
      call assert( 581434754, pos <= pack_size )
    end if

    call musica_mpi_bcast( pack_size, comm )
    if( musica_mpi_rank( comm ) .ne. 0 ) allocate( buffer( pack_size ) )
    call musica_mpi_bcast( buffer, comm )

    if( musica_mpi_rank( comm ) .ne. 0 ) then
      pos = 0
      call type_name%mpi_unpack( buffer, pos, comm )
      my_grid => grid_allocate( type_name )
      call my_grid%mpi_unpack( buffer, pos, comm )
      call assert( 916936550, pos <= pack_size )
    end if
    deallocate( buffer )

    call assert( 340004199, my_grid%handle_ == "foo" )
    call assert( 282165640, my_grid%units( ) == "bars" )
    call assert( 119178577, size( my_grid%mid_ )   == 3 )
    call assert( 508707864, size( my_grid%edge_ )  == 4 )
    call assert( 956075710, size( my_grid%delta_ ) == 3 )

    select type( my_grid )
    class is( grid_from_host_t )
      my_updater = grid_updater_t( my_grid )
    class default
      call die( 945999391 )
    end select

    ! pass edges and mid points
    edges(:) = (/ 0.5_dk,  9.8_dk, 15.4_dk, 45.0_dk /)
    mids(:)  = (/ 1.0_dk, 12.3_dk, 32.4_dk /)
    ! calculated
    deltas(:) = (/ 9.8_dk - 0.5_dk, 15.4_dk - 9.8_dk, 45.0_dk - 15.4_dk /)
    call my_updater%update( mid_points = mids, edges = edges )
    call check_values( 238395985,   my_grid%mid_,   mids, tol )
    call check_values( 575351020,  my_grid%edge_,  edges, tol )
    call check_values( 793270164, my_grid%delta_, deltas, tol )

    ! pass edges only
    edges(:)  = (/ 12.0_dk, 20.0_dk, 30.0_dk, 60.0_dk /)
    ! calculated
    mids(:)   = (/ 16.0_dk, 25.0_dk, 45.0_dk /)
    deltas(:) = (/  8.0_dk, 10.0_dk, 30.0_dk /)
    call my_updater%update( edges = edges )
    call check_values( 453607393,   my_grid%mid_,   mids, tol )
    call check_values( 965983639,  my_grid%edge_,  edges, tol )
    call check_values( 513351486, my_grid%delta_, deltas, tol )

    deallocate( my_grid )

  end subroutine test_grid_from_host_t

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program test_grid_from_host
