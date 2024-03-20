! Copyright (C) 2022 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> Tests for the musica_map module

!> Test module for the musica_map module
program test_util_map

  use musica_constants,                only : dk => musica_dk
  use musica_assert
  use musica_map
  use musica_mpi
#ifdef MUSICA_USE_OPENMP
  use omp_lib
#endif

  implicit none

  character(len=256) :: failure_test_type

  call musica_mpi_init( )

  if( command_argument_count( ) == 0 ) then
    call test_map_t( )
  else if( command_argument_count( ) == 1 ) then
    call get_command_argument( 1, failure_test_type )
    call failure_test( failure_test_type )
  else
    call die( 725972035 )
  end if

  call musica_mpi_finalize( )

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Test map_t functionality
  subroutine test_map_t( )

    use iso_fortran_env,               only : output_unit
    use musica_config,                 only : config_t
    use musica_string,                 only : string_t

    character(len=*), parameter :: my_name = "map_t tests"
    type(map_t) :: map
    type(string_t), allocatable :: from_labels(:), to_labels(:)
    real(kind=dk), allocatable :: from(:), to(:), omp_to(:,:)
    type(config_t) :: config
    character, allocatable :: buffer(:)
    integer :: pos, pack_size, i_thread, n_threads
    integer, parameter :: comm = MPI_COMM_WORLD

    config = '{'//                                                            &
             '  "pairs" : ['//                                                &
             '    {'//                                                        &
             '       "from": "foo",'//                                        &
             '       "to": "bar",'//                                          &
             '       "scale by": 2.5'//                                       &
             '    },'//                                                       &
             '    {'//                                                        &
             '       "from": "foo",'//                                        &
             '       "to": "baz"'//                                           &
             '    }'//                                                        &
             '  ]'//                                                          &
             '}'
    allocate( from_labels( 1 ) )
    allocate( to_labels(   2 ) )
    from_labels(1) = "foo"
    to_labels(1)   = "bar"
    to_labels(2)   = "baz"

    if( musica_mpi_rank( comm ) == 0 ) then
      map = map_t( config, from_labels, to_labels )
      pack_size = map%pack_size( comm )
      allocate( buffer( pack_size ) )
      pos = 0
      call map%mpi_pack( buffer, pos, comm )
      call assert( 796105167, pos == pack_size )
    end if

    call musica_mpi_bcast( pack_size, comm )
    if( musica_mpi_rank( comm ) .ne. 0 ) allocate( buffer( pack_size ) )
    call musica_mpi_bcast( buffer, comm )

    if( musica_mpi_rank( comm ) .ne. 0 ) then
      pos = 0
      call map%mpi_unpack( buffer, pos, comm )
      call assert( 124100631, pos == pack_size )
    end if

    from = (/ 10.0_dk /)
    allocate( to( 2 ) )

    call map%apply( from, to )

    call assert( 671804969, almost_equal( to(1), 25.0_dk ) )
    call assert( 338661002, almost_equal( to(2), 10.0_dk ) )
    deallocate( from_labels )
    deallocate( to_labels   )
    deallocate( from        )
    deallocate( to          )

    config = '{'//                                                            &
             '  "match full source": false,'//                                &
             '  "match full destination": false,'//                           &
             '  "sum multiple matches": true,'//                              &
             '  "pairs" : ['//                                                &
             '    {'//                                                        &
             '       "from": "foo",'//                                        &
             '       "to": "bar",'//                                          &
             '       "scale by": 2.5'//                                       &
             '    },'//                                                       &
             '    {'//                                                        &
             '       "from": "foo",'//                                        &
             '       "to": "baz"'//                                           &
             '    },'//                                                       &
             '    {'//                                                        &
             '       "from": "bar",'//                                        &
             '       "to": "bar"'//                                           &
             '    }'//                                                        &
             '  ]'//                                                          &
             '}'
    allocate( from_labels( 3 ) )
    allocate( to_labels(   3 ) )
    from_labels(1) = "foo"
    from_labels(2) = "bar"
    from_labels(3) = "baz"
    to_labels(1)   = "bar"
    to_labels(2)   = "baz"
    to_labels(3)   = "quz"

    from = (/ 10.0_dk, 20.0_dk, 30.0_dk /)
    
#ifdef MUSICA_USE_OPENMP
    n_threads = omp_get_num_threads( )
#else
    n_threads = 1
#endif

    allocate( omp_to( n_threads, 3 ) )

    map = map_t( config, from_labels, to_labels )

    !$omp parallel do private(i_thread)
    do i_thread = 1, n_threads
      call check_omp_case( map, from, omp_to( i_thread, : ) )
    end do
    !$omp end parallel do
    
    deallocate( from_labels )
    deallocate( to_labels   )
    deallocate( from        )
    deallocate( omp_to      )

    config = '{'//                                                            &
             '  "match full source": false,'//                                &
             '  "sum multiple matches": true,'//                              &
             '  "default matching": "backup",'//                              &
             '  "pairs" : ['//                                                &
             '    {'//                                                        &
             '       "from": "foo",'//                                        &
             '       "to": "bar",'//                                          &
             '       "scale by": 2.5'//                                       &
             '    },'//                                                       &
             '    {'//                                                        &
             '       "from": "foo",'//                                        &
             '       "to": "baz"'//                                           &
             '    },'//                                                       &
             '    {'//                                                        &
             '       "from": "bar",'//                                        &
             '       "to": "bar"'//                                           &
             '    }'//                                                        &
             '  ]'//                                                          &
             '}'
    allocate( from_labels( 3 ) )
    allocate( to_labels(   3 ) )
    from_labels(1) = "foo"
    from_labels(2) = "bar"
    from_labels(3) = "quz"
    to_labels(1)   = "bar"
    to_labels(2)   = "baz"
    to_labels(3)   = "quz"

    from = (/ 10.0_dk, 20.0_dk, 30.0_dk /)
    allocate( to( 3 ) )

    map = map_t( config, from_labels, to_labels )
    call map%apply( from, to )

    call assert( 393064014, almost_equal( to(1), 45.0_dk ) )
    call assert( 157898710, almost_equal( to(2), 10.0_dk ) )
    call assert( 952750205, almost_equal( to(3), 30.0_dk ) )
    deallocate( from_labels )
    deallocate( to_labels   )
    deallocate( from        )
    deallocate( to          )

    config = '{'//                                                            &
             '  "match full source": false,'//                                &
             '  "sum multiple matches": true,'//                              &
             '  "default matching": "always",'//                              &
             '  "pairs" : ['//                                                &
             '    {'//                                                        &
             '       "from": "foo",'//                                        &
             '       "to": "bar",'//                                          &
             '       "scale by": 2.5'//                                       &
             '    },'//                                                       &
             '    {'//                                                        &
             '       "from": "foo",'//                                        &
             '       "to": "baz"'//                                           &
             '    }'//                                                        &
             '  ]'//                                                          &
             '}'
    allocate( from_labels( 3 ) )
    allocate( to_labels(   3 ) )
    from_labels(1) = "foo"
    from_labels(2) = "bar"
    from_labels(3) = "quz"
    to_labels(1)   = "bar"
    to_labels(2)   = "baz"
    to_labels(3)   = "quz"

    from = (/ 10.0_dk, 20.0_dk, 30.0_dk /)
    allocate( to( 3 ) )

    map = map_t( config, from_labels, to_labels )
    call map%apply( from, to )
    if( musica_mpi_rank( comm ) == 0 ) then
      call map%print( from_labels, to_labels, output_unit )
    end if

    call assert( 884835327, almost_equal( to(1), 45.0_dk ) )
    call assert( 432203174, almost_equal( to(2), 10.0_dk ) )
    call assert( 597095771, almost_equal( to(3), 30.0_dk ) )
    deallocate( from_labels )
    deallocate( to_labels   )
    deallocate( from        )
    deallocate( to          )

  end subroutine test_map_t

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Checks a specific case on parallel OpenMP threads (if available)
  subroutine check_omp_case( map, from, to )

    type(map_t),   intent(in)    :: map
    real(kind=dk), intent(in)    :: from(:)
    real(kind=dk), intent(inout) :: to(:)

    call map%apply( from, to )

    call assert( 162210850, almost_equal( to(1), 45.0_dk ) )
    call assert( 495807112, almost_equal( to(2), 10.0_dk ) )
    call assert( 943174958, almost_equal( to(3),  0.0_dk ) )

  end subroutine check_omp_case

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Failure tests for map_t class
  subroutine failure_test( test_type )

    character(len=*), intent(in) :: test_type

    if( test_type .eq. "170733942" ) then
      call failure_test_170733942( )
    else if( test_type .eq. "764798475" ) then
      call failure_test_764798475( )
    else if( test_type .eq. "133386338" ) then
      call failure_test_133386338( )
    else if( test_type .eq. "956987954" ) then
      call failure_test_956987954( )
    else if( test_type .eq. "200274675" ) then
      call failure_test_200274675( )
    else if( test_type .eq. "240867074" ) then
      call failure_test_240867074( )
    else if( test_type .eq. "309595761" ) then
      call failure_test_309595761( )
    else if( test_type .eq. "122570601" ) then
      call failure_test_122570601( )
    else if( test_type .eq. "740547646" ) then
      call failure_test_740547646( )
    else if( test_type .eq. "548594113" ) then
      call failure_test_548594113( )
    else
      call die( 609154398 )
    end if

  end subroutine failure_test

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Test invalid map configuration
  subroutine failure_test_170733942( )

    use musica_config,                 only : config_t
    use musica_string,                 only : string_t

    type(map_t) :: map
    type(config_t) :: config
    type(string_t) :: from_labels(2), to_labels(1)

    from_labels(1) = "foo"
    from_labels(2) = "bar"
    to_labels(1)   = "baz"
    config = '{ "bad": "config" }'
    map = map_t( config, from_labels, to_labels )

  end subroutine failure_test_170733942

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Test invalid map pair configuration
  subroutine failure_test_309595761( )

    use musica_config,                 only : config_t
    use musica_string,                 only : string_t

    type(map_t) :: map
    type(config_t) :: config
    type(string_t) :: from_labels(2), to_labels(1)

    from_labels(1) = "foo"
    from_labels(2) = "bar"
    to_labels(1)   = "baz"
    config = '{ "pairs": [ { "bad": "config" } ] }'
    map = map_t( config, from_labels, to_labels )

  end subroutine failure_test_309595761

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Test invalid default matching options
  subroutine failure_test_548594113( )

    use musica_config,                 only : config_t
    use musica_string,                 only : string_t

    type(map_t) :: map
    type(config_t) :: config
    type(string_t) :: from_labels(2), to_labels(1)

    from_labels(1) = "foo"
    from_labels(2) = "bar"
    to_labels(1)   = "baz"
    config = '{ "match full destination": false, '//                          &
             '  "default matching": "always", '//                             &
             '  "pairs": [ { "from": "foo", "to": "baz" } ] }'
    map = map_t( config, from_labels, to_labels )

  end subroutine failure_test_548594113

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Test missing source element
  subroutine failure_test_122570601( )

    use musica_config,                 only : config_t
    use musica_string,                 only : string_t

    type(map_t) :: map
    type(config_t) :: config
    type(string_t) :: from_labels(2), to_labels(1)

    from_labels(1) = "foo"
    from_labels(2) = "bar"
    to_labels(1)   = "baz"
    config = '{ "pairs": [ { "from": "quz", "to": "baz" } ] }'
    map = map_t( config, from_labels, to_labels )

  end subroutine failure_test_122570601

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Test missing destination element
  subroutine failure_test_740547646( )

    use musica_config,                 only : config_t
    use musica_string,                 only : string_t

    type(map_t) :: map
    type(config_t) :: config
    type(string_t) :: from_labels(2), to_labels(1)

    from_labels(1) = "foo"
    from_labels(2) = "bar"
    to_labels(1)   = "baz"
    config = '{ "pairs": [ { "from": "foo", "to": "bar" } ] }'
    map = map_t( config, from_labels, to_labels )

  end subroutine failure_test_740547646

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Test unmatched source element
  subroutine failure_test_956987954( )

    use musica_config,                 only : config_t
    use musica_string,                 only : string_t

    type(map_t) :: map
    type(config_t) :: config
    type(string_t) :: from_labels(2), to_labels(1)

    from_labels(1) = "foo"
    from_labels(2) = "bar"
    to_labels(1)   = "baz"
    config = '{ "pairs": [ { "from": "foo", "to": "baz" } ] }'
    map = map_t( config, from_labels, to_labels )

  end subroutine failure_test_956987954

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Test unmatched destination element
  subroutine failure_test_200274675( )

    use musica_config,                 only : config_t
    use musica_string,                 only : string_t

    type(map_t) :: map
    type(config_t) :: config
    type(string_t) :: from_labels(1), to_labels(2)

    from_labels(1) = "foo"
    to_labels(1)   = "baz"
    to_labels(2)   = "quz"
    config = '{ "pairs": [ { "from": "foo", "to": "quz" } ] }'
    map = map_t( config, from_labels, to_labels )

  end subroutine failure_test_200274675

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Test multiple destination element matches
  subroutine failure_test_240867074( )

    use musica_config,                 only : config_t
    use musica_string,                 only : string_t

    type(map_t) :: map
    type(config_t) :: config
    type(string_t) :: from_labels(2), to_labels(1)

    from_labels(1) = "foo"
    from_labels(2) = "bar"
    to_labels(1)   = "baz"
    config = '{ "pairs": [ { "from": "foo", "to": "baz" }, '//                &
                         ' { "from": "bar", "to": "baz" } ] }'
    map = map_t( config, from_labels, to_labels )

  end subroutine failure_test_240867074

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Test wrong source array size
  subroutine failure_test_764798475( )

    use musica_config,                 only : config_t
    use musica_string,                 only : string_t

    type(map_t) :: map
    type(config_t) :: config
    type(string_t) :: from_labels(1), to_labels(2)
    real(kind=dk) :: from(3), to(2)

    from_labels(1) = "foo"
    to_labels(1)   = "baz"
    to_labels(2)   = "quz"
    config = '{ "pairs": [ { "from": "foo", "to": "quz" }, '//                &
                          '{ "from": "foo", "to": "baz" } ] }'
    map = map_t( config, from_labels, to_labels )
    call map%apply( from, to )

  end subroutine failure_test_764798475

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Test wrong destination array size
  subroutine failure_test_133386338( )

    use musica_config,                 only : config_t
    use musica_string,                 only : string_t

    type(map_t) :: map
    type(config_t) :: config
    type(string_t) :: from_labels(1), to_labels(2)
    real(kind=dk) :: from(1), to(1)

    from_labels(1) = "foo"
    to_labels(1)   = "baz"
    to_labels(2)   = "quz"
    config = '{ "pairs": [ { "from": "foo", "to": "quz" }, '//                &
                          '{ "from": "foo", "to": "baz" } ] }'
    map = map_t( config, from_labels, to_labels )
    call map%apply( from, to )

  end subroutine failure_test_133386338

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program test_util_map
