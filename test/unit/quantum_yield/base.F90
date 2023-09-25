! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
program test_quantum_yield

  use musica_mpi,                      only : musica_mpi_init,                &
                                              musica_mpi_finalize
  use tuvx_quantum_yield

  implicit none

  call musica_mpi_init( )
  call test_quantum_yield_mpi( )
  call test_quantum_yield_t( )
  call musica_mpi_finalize( )

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine test_quantum_yield_mpi( )
    ! Test functionality of the :f:type:`~tuvx_quantum_yield/quantum_yield_t`
    ! class.
    !
    ! This test only checks the MPI functions currently. Additional test will
    ! be added as part of issue #177

    use musica_assert,                 only : assert
    use musica_constants,              only : dk => musica_dk
    use musica_mpi
    use musica_string,                 only : string_t
    use tuvx_quantum_yield_factory,    only : quantum_yield_type_name,        &
                                              quantum_yield_allocate
    use tuvx_test_utils

    class(quantum_yield_t), pointer :: quantum_yield
    character, allocatable :: buffer(:)
    type(string_t) :: type_name
    integer :: pos, pack_size
    integer, parameter :: comm = MPI_COMM_WORLD

    real(dk) :: temperature1(3), temperature2(4)
    real(dk) :: array1(3,2), array2(2,1)

    temperature1(:) = (/  12.5_dk,  16.3_dk, -290.4_dk /)
    array1(:,1)     = (/  42.3_dk, 132.4_dk,   13.4_dk /)
    array1(:,2)     = (/ 132.4_dk,  0.43_dk,   2.34_dk /)

    temperature2(:) = (/ -123.4_dk, 41.2_dk, 0.053_dk, 1.2e-7_dk /)
    array2(:,1)     = (/ 12.34_dk, -142.3_dk /)

    if( musica_mpi_rank( comm ) == 0 ) then
      allocate( quantum_yield )
      allocate( quantum_yield%quantum_yield_parms(2) )
      quantum_yield%quantum_yield_parms(1)%temperature = temperature1
      quantum_yield%quantum_yield_parms(1)%array = array1
      quantum_yield%quantum_yield_parms(2)%temperature = temperature2
      quantum_yield%quantum_yield_parms(2)%array = array2
      type_name = quantum_yield_type_name( quantum_yield )
      pack_size = type_name%pack_size( comm ) + quantum_yield%pack_size( comm )
      allocate( buffer( pack_size ) )
      pos = 0
      call type_name%mpi_pack(     buffer, pos , comm )
      call quantum_yield%mpi_pack( buffer, pos , comm )
      call assert( 209765802, pos <= pack_size )
    end if

    call musica_mpi_bcast( pack_size , comm )
    if( musica_mpi_rank( comm ) .ne. 0 ) allocate( buffer( pack_size ) )
    call musica_mpi_bcast( buffer , comm )

    if( musica_mpi_rank( comm ) .ne. 0 ) then
      pos = 0
      call type_name%mpi_unpack( buffer, pos , comm )
      quantum_yield => quantum_yield_allocate( type_name )
      call quantum_yield%mpi_unpack( buffer, pos , comm )
      call assert( 264697883, pos <= pack_size )
    end if
    deallocate( buffer )

    call assert( 962406487, associated( quantum_yield ) )
    call assert( 964312021, allocated( quantum_yield%quantum_yield_parms ) )
    call assert( 401267057, size( quantum_yield%quantum_yield_parms ) == 2 )
    associate( params => quantum_yield%quantum_yield_parms( 1 ) )
      call assert( 784078798, allocated( params%temperature ) )
      call assert( 891132836, allocated( params%array ) )
      call check_values_1D_no_code( params%temperature, temperature1, 1.0e-6_dk )
      call check_values_2D_no_code( params%array,       array1,       1.0e-6_dk )
    end associate
    associate( params => quantum_yield%quantum_yield_parms( 2 ) )
      call assert( 784078798, allocated( params%temperature ) )
      call assert( 891132836, allocated( params%array ) )
      call check_values_1D_no_code( params%temperature, temperature2, 1.0e-6_dk )
      call check_values_2D_no_code( params%array,       array2,       1.0e-6_dk )
    end associate

    deallocate( quantum_yield )

  end subroutine test_quantum_yield_mpi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine test_quantum_yield_t( )

    use musica_assert,                 only : assert
    use musica_constants,              only : dk => musica_dk
    use musica_config,                 only : config_t
    use musica_iterator,               only : iterator_t
    use musica_mpi
    use musica_string,                 only : string_t
    use tuvx_quantum_yield_factory,    only : quantum_yield_type_name,        &
                                              quantum_yield_allocate
    use tuvx_grid_warehouse,           only : grid_warehouse_t
    use tuvx_profile_warehouse,        only : profile_warehouse_t

    class(grid_warehouse_t),    pointer :: grids
    class(profile_warehouse_t), pointer :: profiles
    class(quantum_yield_t),     pointer :: quantum_yield

    character(len=*), parameter :: my_name = "base quantum yield tests"
    type(config_t) :: config, qy_set, qy_config
    class(iterator_t), pointer :: iter
    real(kind=dk), allocatable :: results(:,:)
    real(kind=dk), allocatable :: input(:), input_grid(:)
    real(kind=dk) :: input_base(4) =                                          &
      (/ 0.5_dk, 0.1_dk, 0.4_dk, 1.0_dk /)
    real(kind=dk) :: input_grid_base(4) =                                     &
      (/ 101.0_dk, 102.0_dk, 103.0_dk, 104.0_dk /)
    integer :: i_height
    character, allocatable :: buffer(:)
    type(string_t) :: type_name
    integer :: pos, pack_size
    integer, parameter :: comm = MPI_COMM_WORLD

    ! load test grids
    call config%from_file( "test/data/grid.la_sr.config.json" )
    grids => grid_warehouse_t( config )

    ! load test profiles
    call config%from_file( "test/data/profile.simple.config.json" )
    profiles => profile_warehouse_t( config, grids )

    ! get quantum yield config data
    call config%from_file( "test/data/quantum_yields/base.config.json" )
    call config%get( "quantum yields", qy_set, my_name )
    iter => qy_set%get_iterator( )

    ! load and test quantum yield w/o extrapolation
    call assert( 398844056, iter%next( ) )
    call qy_set%get( iter, qy_config, my_name )
    if( musica_mpi_rank( comm ) == 0 ) then
      quantum_yield => quantum_yield_t( qy_config, grids, profiles )
      type_name = quantum_yield_type_name( quantum_yield )
      pack_size = type_name%pack_size( comm ) + quantum_yield%pack_size( comm )
      allocate( buffer( pack_size ) )
      pos = 0
      call type_name%mpi_pack(     buffer, pos, comm )
      call quantum_yield%mpi_pack( buffer, pos, comm )
      call assert( 725722772, pos <= pack_size )
    end if

    call musica_mpi_bcast( pack_size, comm )
    if( musica_mpi_rank( comm ) .ne. 0 ) allocate( buffer( pack_size ) )
    call musica_mpi_bcast( buffer, comm )

    if( musica_mpi_rank( comm ) .ne. 0 ) then
      pos = 0
      call type_name%mpi_unpack( buffer, pos, comm )
      quantum_yield => quantum_yield_allocate( type_name )
      call quantum_yield%mpi_unpack( buffer, pos, comm )
      call assert( 705570134, pos <= pack_size )
    end if
    deallocate( buffer )

    results = quantum_yield%calculate( grids, profiles )
    input = input_base
    input_grid = input_grid_base
    call add_points( input, input_grid, 0.0_dk, 0.0_dk )
    call check_values( results(:,1:4), input, input_grid, 6 )
    deallocate( input )
    deallocate( input_grid )
    deallocate( quantum_yield )

    ! load and test quantum yield w/ band overrides
    call assert( 675215609, iter%next( ) )
    call qy_set%get( iter, qy_config, my_name )
    if( musica_mpi_rank( comm ) == 0 ) then
      quantum_yield => quantum_yield_t( qy_config, grids, profiles )
      type_name = quantum_yield_type_name( quantum_yield )
      pack_size = type_name%pack_size( comm ) + quantum_yield%pack_size( comm )
      allocate( buffer( pack_size ) )
      pos = 0
      call type_name%mpi_pack(     buffer, pos, comm )
      call quantum_yield%mpi_pack( buffer, pos, comm )
      call assert( 505058705, pos <= pack_size )
    end if

    call musica_mpi_bcast( pack_size, comm )
    if( musica_mpi_rank( comm ) .ne. 0 ) allocate( buffer( pack_size ) )
    call musica_mpi_bcast( buffer, comm )

    if( musica_mpi_rank( comm ) .ne. 0 ) then
      pos = 0
      call type_name%mpi_unpack( buffer, pos, comm )
      quantum_yield => quantum_yield_allocate( type_name )
      call quantum_yield%mpi_unpack( buffer, pos, comm )
      call assert( 617377050, pos <= pack_size )
    end if
    deallocate( buffer )

    results = quantum_yield%calculate( grids, profiles )
    input = input_base
    input_grid = input_grid_base
    do i_height = 1, 6
      call assert( 836749433, results( i_height, 6  ) == 0.423_dk )
      call assert( 384117280, results( i_height, 7  ) == 0.932_dk )
      call assert( 896493526, results( i_height, 8  ) == 0.932_dk )
      call assert( 443861373, results( i_height, 9  ) == 0.122_dk )
      call assert( 891229219, results( i_height, 10 ) == 0.122_dk )
      call assert( 438597066, results( i_height, 11 ) ==   0.0_dk )
    end do
    call add_points( input, input_grid, 0.0_dk, 0.0_dk )
    call check_values( results(:,1:4), input, input_grid, 6 )
    deallocate( input )
    deallocate( input_grid )
    deallocate( quantum_yield )

    ! load and test quantum yield w/ fixed lower extrapolation and no upper
    ! extrapolation
    call assert( 772193328, iter%next( ) )
    call qy_set%get( iter, qy_config, my_name )
    if( musica_mpi_rank( comm ) == 0 ) then
      quantum_yield => quantum_yield_t( qy_config, grids, profiles )
      type_name = quantum_yield_type_name( quantum_yield )
      pack_size = type_name%pack_size( comm ) + quantum_yield%pack_size( comm )
      allocate( buffer( pack_size ) )
      pos = 0
      call type_name%mpi_pack(     buffer, pos, comm )
      call quantum_yield%mpi_pack( buffer, pos, comm )
      call assert( 884511673, pos <= pack_size )
    end if

    call musica_mpi_bcast( pack_size, comm )
    if( musica_mpi_rank( comm ) .ne. 0 ) allocate( buffer( pack_size ) )
    call musica_mpi_bcast( buffer, comm )

    if( musica_mpi_rank( comm ) .ne. 0 ) then
      pos = 0
      call type_name%mpi_unpack( buffer, pos, comm )
      quantum_yield => quantum_yield_allocate( type_name )
      call quantum_yield%mpi_unpack( buffer, pos, comm )
      call assert( 996830018, pos <= pack_size )
    end if
    deallocate( buffer )

    results = quantum_yield%calculate( grids, profiles )
    input = input_base
    input_grid = input_grid_base
    call add_points( input, input_grid, 0.125_dk, 0.0_dk )
    call check_values( results(:,1:4), input, input_grid, 6 )
    deallocate( input )
    deallocate( input_grid )
    deallocate( quantum_yield )

    ! load and test quantum yield w/ extrapolation from lower boundary and
    ! fixed upper extrapolation
    call assert( 374493256, iter%next( ) )
    call qy_set%get( iter, qy_config, my_name )
    if( musica_mpi_rank( comm ) == 0 ) then
      quantum_yield => quantum_yield_t( qy_config, grids, profiles )
      type_name = quantum_yield_type_name( quantum_yield )
      pack_size = type_name%pack_size( comm ) + quantum_yield%pack_size( comm )
      allocate( buffer( pack_size ) )
      pos = 0
      call type_name%mpi_pack(     buffer, pos, comm )
      call quantum_yield%mpi_pack( buffer, pos, comm )
      call assert( 139327952, pos <= pack_size )
    end if

    call musica_mpi_bcast( pack_size, comm )
    if( musica_mpi_rank( comm ) .ne. 0 ) allocate( buffer( pack_size ) )
    call musica_mpi_bcast( buffer, comm )

    if( musica_mpi_rank( comm ) .ne. 0 ) then
      pos = 0
      call type_name%mpi_unpack( buffer, pos, comm )
      quantum_yield => quantum_yield_allocate( type_name )
      call quantum_yield%mpi_unpack( buffer, pos, comm )
      call assert( 534121546, pos <= pack_size )
    end if
    deallocate( buffer )

    results = quantum_yield%calculate( grids, profiles )
    input = input_base
    input_grid = input_grid_base
    call add_points( input, input_grid, 0.5_dk, 0.323_dk )
    call check_values( results(:,1:4), input, input_grid, 6 )
    deallocate( input )
    deallocate( input_grid )
    deallocate( quantum_yield )

    ! clean up
    deallocate( iter )
    deallocate( grids )
    deallocate( profiles )

  end subroutine test_quantum_yield_t

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Adds additional points input data
  subroutine add_points( values, grid, lower_val, upper_val )

    use musica_constants,              only : dk => musica_dk
    use tuvx_util,                     only : add_point

    real(kind=dk), allocatable, intent(inout) :: values(:)
    real(kind=dk), allocatable, intent(inout) :: grid(:)
    real(kind=dk),              intent(in)    :: lower_val
    real(kind=dk),              intent(in)    :: upper_val

    call add_point( x = grid, y = values,                                     &
                    xnew = ( 1.0_dk - 1.0e-5_dk ) * grid(1), ynew = lower_val )
    call add_point( x = grid, y = values,                                     &
                    xnew = 0.0_dk, ynew = lower_val )
    call add_point( x = grid, y = values,                                     &
                    xnew = ( 1.0_dk + 1.0e-5_dk ) * grid( size( grid ) ),     &
                    ynew = upper_val )
    call add_point( x = grid, y = values,                                     &
                    xnew = 1.0e38_dk, ynew = upper_val )

  end subroutine add_points

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Checks results against expectations
  subroutine check_values( results, input, input_grid, n_levels )

    use musica_assert,                 only : assert, almost_equal
    use musica_constants,              only : dk => musica_dk
    use tuvx_interpolate,              only : interpolator_conserving_t

    real(kind=dk), intent(in) :: results(:,:)
    real(kind=dk), intent(in) :: input(:)
    real(kind=dk), intent(in) :: input_grid(:)
    integer,       intent(in) :: n_levels

    integer :: i_level, i_wavelength
    real(kind=dk) :: output_grid(5) =                                         &
      (/ 100.0_dk, 101.5_dk, 102.0_dk, 103.0_dk, 104.5_dk /)
    real(kind=dk) :: output(4)
    type(interpolator_conserving_t) :: interpolator

    call assert( 896869030, size( results, dim = 1 ) == n_levels )
    call assert( 444236877, size( results, dim = 2 ) == size( output ) )
    do i_wavelength = 1, size( results, dim = 2 )
      do i_level = 2, size( results, dim = 1 )
        call assert( 956613123, results( i_level, i_wavelength ) ==           &
                                results( 1,       i_wavelength ) )
      end do
    end do

    ! x data are assumed to be at interfaces and y data at midpoints
    output = interpolator%interpolate( x_target = output_grid,                &
                                       x_source = input_grid,                 &
                                       y_source = input )
    call assert( 103923069, almost_equal( results( 1, 1 ), output( 1 ) ) )
    call assert( 898774564, almost_equal( results( 1, 2 ), output( 2 ) ) )
    call assert( 446142411, almost_equal( results( 1, 3 ), output( 3 ) ) )
    call assert( 893510257, almost_equal( results( 1, 4 ), output( 4 ) ) )

  end subroutine check_values

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program test_quantum_yield
