! Copyright (C) 2007-2021 Barcelona Supercomputing Center and University of
! Illinois at Urbana-Champaign
! SPDX-License-Identifier: MIT
program test_mpi
  ! Tests for MPI wrapper functions.
  !
  ! This module was adapted from CAMP (https://github.com/open-atmos/camp).

  use musica_assert
  use musica_mpi

  implicit none

  call test_mpi_wrappers( )

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine test_mpi_wrappers( )
    ! test MPI wrapper functions

    use musica_constants,              only : dp => musica_dk
    use musica_string,                 only : to_char

#ifdef MUSICA_USE_MPI
    integer, parameter :: comm = MPI_COMM_WORLD
    integer, parameter :: dc = dp
    real(kind=dp), parameter :: test_real = 2.718281828459d0
    complex(kind=dc), parameter :: test_complex &
         = (0.707106781187d0, 1.4142135624d0)
    logical, parameter :: test_logical = .true.
    character(len=100), parameter :: test_string &
         = "a truth universally acknowledged"
    integer, parameter :: test_integer = 314159

    character, allocatable :: buffer(:) ! memory buffer
    integer :: buffer_size, max_buffer_size, position
    real(kind=dp) :: send_real, recv_real
    complex(kind=dc) :: send_complex, recv_complex
    logical :: send_logical, recv_logical
    character(len=100) :: send_string, recv_string
    integer :: send_integer, recv_integer
    integer :: test_integer_array(2) = (/ 4, 2 /)
    integer, allocatable :: send_integer_array(:)
    integer, allocatable :: recv_integer_array(:)
    real(kind=dp), allocatable :: send_real_array(:)
    real(kind=dp), allocatable :: recv_real_array(:)
    character(len=5) :: test_string_array(2) = (/ "forty", "two  " /)
    character(len=5), allocatable :: send_string_array(:)
    character(len=5), allocatable :: recv_string_array(:)
    real(kind=dp) :: test_real_array_2d(2,2)
    real(kind=dp), allocatable :: send_real_array_2d(:,:)
    real(kind=dp), allocatable :: recv_real_array_2d(:,:)
    real(kind=dp) :: test_real_array_3d(2,2,2)
    real(kind=dp), allocatable :: send_real_array_3d(:,:,:)
    real(kind=dp), allocatable :: recv_real_array_3d(:,:,:)

    test_real_array_2d(1,1) = 42.0_dp
    test_real_array_2d(2,1) = 4.2_dp
    test_real_array_2d(1,2) = 0.42_dp
    test_real_array_2d(2,2) = 0.042_dp

    test_real_array_3d(1,1,1) = 412.3_dp
    test_real_array_3d(2,1,1) = 312.0_dp
    test_real_array_3d(1,2,1) = 212.9_dp
    test_real_array_3d(2,2,1) = 132.8_dp
    test_real_array_3d(1,1,2) = 312.7_dp
    test_real_array_3d(2,1,2) = 712.6_dp
    test_real_array_3d(1,2,2) = 452.2_dp
    test_real_array_3d(2,2,2) = 912.3_dp

    call assert( 357761664, musica_mpi_support( ) )
    call musica_mpi_init( )

    call assert( 455191678, musica_mpi_size( comm ) > 1 )

    call musica_mpi_barrier( comm )

    send_integer = 0
    if( musica_mpi_rank( comm ) == 0 ) send_integer = 42
    call musica_mpi_bcast( send_integer, comm )
    call assert( 353714667, send_integer == 42 )

    send_string = ""
    if( musica_mpi_rank( comm ) == 0 ) send_string = "forty two"
    call musica_mpi_bcast( send_string, comm )
    call assert( 904777778, trim( send_string ) == "forty two" )

    buffer = (/ 'x', 'x' /)
    if( musica_mpi_rank( comm ) == 0 ) buffer(:) = (/ '4', '2' /)
    call musica_mpi_bcast( buffer, comm )
    call assert( 954445552, buffer(1) == '4' )
    call assert( 496549092, buffer(2) == '2' )
    deallocate( buffer )

    if( musica_mpi_rank( comm ) == 0 ) then
       send_real = test_real
       send_complex = test_complex
       send_logical = test_logical
       send_string = test_string
       send_integer = test_integer
       allocate( send_real_array(2) )
       send_real_array(1) = real( test_complex )
       send_real_array(2) = aimag( test_complex )
       allocate( send_integer_array(2) )
       send_integer_array(:) = test_integer_array(:)
       allocate( send_string_array(2) )
       send_string_array(:) = test_string_array(:)
       allocate( send_real_array_2d(2,2) )
       send_real_array_2d(:,:) = test_real_array_2d(:,:)
       allocate( send_real_array_3d(2,2,2) )
       send_real_array_3d(:,:,:) = test_real_array_3d(:,:,:)

       max_buffer_size = 0
       max_buffer_size = max_buffer_size                                      &
            + musica_mpi_pack_size( send_integer, comm )
       max_buffer_size = max_buffer_size                                      &
            + musica_mpi_pack_size( send_real, comm )
       max_buffer_size = max_buffer_size                                      &
            + musica_mpi_pack_size( send_complex, comm )
       max_buffer_size = max_buffer_size                                      &
            + musica_mpi_pack_size( send_logical, comm )
       max_buffer_size = max_buffer_size                                      &
            + musica_mpi_pack_size( send_string, comm )
       max_buffer_size = max_buffer_size                                      &
            + musica_mpi_pack_size( send_real_array, comm )
       max_buffer_size = max_buffer_size                                      &
            + musica_mpi_pack_size( send_integer_array, comm )
       max_buffer_size = max_buffer_size                                      &
            + musica_mpi_pack_size( send_string_array, comm )
       max_buffer_size = max_buffer_size                                      &
            + musica_mpi_pack_size( send_real_array_2d, comm )
       max_buffer_size = max_buffer_size                                      &
            + musica_mpi_pack_size( send_real_array_3d, comm )

       allocate( buffer( max_buffer_size ) )

       position = 0
       call musica_mpi_pack( buffer, position, send_real         , comm )
       call musica_mpi_pack( buffer, position, send_complex      , comm )
       call musica_mpi_pack( buffer, position, send_logical      , comm )
       call musica_mpi_pack( buffer, position, send_string       , comm )
       call musica_mpi_pack( buffer, position, send_integer      , comm )
       call musica_mpi_pack( buffer, position, send_real_array   , comm )
       call musica_mpi_pack( buffer, position, send_integer_array, comm )
       call musica_mpi_pack( buffer, position, send_string_array , comm )
       call musica_mpi_pack( buffer, position, send_real_array_2d, comm )
       call musica_mpi_pack( buffer, position, send_real_array_3d, comm )
       call assert_msg( 350740830, position <= max_buffer_size,               &
                        "MPI test failure: pack position "                    &
                        // trim( to_char( position ) )                        &
                        // " greater than max_buffer_size "                   &
                        // trim( to_char( max_buffer_size ) ) )
       buffer_size = position ! might be less than we allocated
    end if

    call musica_mpi_bcast( buffer_size, comm )

    if( musica_mpi_rank( comm ) /= 0 ) then
       allocate( buffer( buffer_size ) )
    end if

    call musica_mpi_bcast( buffer, comm )

    if( musica_mpi_rank( comm ) /= 0 ) then
       position = 0
       call musica_mpi_unpack( buffer, position, recv_real         , comm )
       call musica_mpi_unpack( buffer, position, recv_complex      , comm )
       call musica_mpi_unpack( buffer, position, recv_logical      , comm )
       call musica_mpi_unpack( buffer, position, recv_string       , comm )
       call musica_mpi_unpack( buffer, position, recv_integer      , comm )
       call musica_mpi_unpack( buffer, position, recv_real_array   , comm )
       call musica_mpi_unpack( buffer, position, recv_integer_array, comm )
       call musica_mpi_unpack( buffer, position, recv_string_array , comm )
       call musica_mpi_unpack( buffer, position, recv_real_array_2d, comm )
       call musica_mpi_unpack( buffer, position, recv_real_array_3d, comm )

       call assert_msg( 787677020, position == buffer_size,                   &
                        "MPI test failure: unpack position "                  &
                        // trim( to_char( position ) )                        &
                        // " not equal to buffer_size "                       &
                        // trim( to_char( buffer_size ) ) )
    end if

    deallocate( buffer )

    if( musica_mpi_rank( comm ) /= 0 ) then
       call assert_msg( 567548916, recv_real == test_real,                    &
                        "MPI test failure: real recv "                        &
                        // trim( to_char( recv_real ) )                       &
                        // " not equal to "                                   &
                        // trim( to_char( test_real ) ) )
       call assert_msg( 653908509, recv_complex == test_complex,              &
                        "MPI test failure: complex recv "                     &
                        // trim( to_char( recv_complex ) )                    &
                        // " not equal to "                                   &
                        // trim( to_char( test_complex ) ) )
       call assert_msg( 307746296, recv_logical .eqv. test_logical,           &
                        "MPI test failure: logical recv "                     &
                        // trim( to_char( recv_logical ) )                    &
                        // " not equal to "                                   &
                        // trim( to_char( test_logical ) ) )
       call assert_msg( 155693492, recv_string == test_string,                &
                        "MPI test failure: string recv '"                     &
                        // trim( recv_string )                                &
                        // "' not equal to '"                                 &
                        // trim( test_string ) // "'" )
       call assert_msg( 875699427, recv_integer == test_integer,              &
                        "MPI test failure: integer recv "                     &
                        // trim( to_char( recv_integer ) )                    &
                        // " not equal to "                                   &
                        // trim( to_char( test_integer ) ) )
       call assert_msg( 326982363, size( recv_real_array ) == 2,              &
                        "MPI test failure: real array recv size " //          &
                        trim( to_char( size( recv_real_array ) ) )            &
                        // " not equal to 2" )
       call assert_msg( 744394323,                                            &
                        recv_real_array(1) == real( test_complex ),           &
                        "MPI test failure: real array recv index 1 "          &
                        // trim( to_char( recv_real_array(1) ) )              &
                        // " not equal to "                                   &
                        // trim( to_char( real( test_complex ) ) ) )
       call assert_msg( 858902527,                                            &
                        recv_real_array(2) == aimag( test_complex ),          &
                        "MPI test failure: real array recv index 2 "          &
                        // trim( to_char( recv_real_array(2) ) )              &
                        // " not equal to "                                   &
                        // trim( to_char( aimag( test_complex ) ) ) )
       call assert_msg( 785767484, size( recv_integer_array ) == 2,           &
                        "MPI test failure: integer array recv size " //       &
                        trim( to_char( size( recv_integer_array ) ) )         &
                        // " not equal to 2" )
       call assert_msg( 874548821,                                            &
                        recv_integer_array(1) == test_integer_array(1),       &
                        "MPI test failure: integer array recv index 1 "       &
                        // trim( to_char( recv_integer_array(1) ) )           &
                        // " not equal to "                                   &
                        // trim( to_char( test_integer_array(1) ) ) )
       call assert_msg( 422368963,                                            &
                        recv_integer_array(2) == test_integer_array(2),       &
                        "MPI test failure: integer array recv index 2 "       &
                        // trim( to_char( recv_integer_array(2) ) )           &
                        // " not equal to "                                   &
                        // trim( to_char( test_integer_array(2) ) ) )
       call assert_msg( 858519971, size( recv_string_array ) == 2,            &
                        "MPI test failure: string array recv size " //        &
                        trim( to_char( size( recv_string_array ) ) )          &
                        // " not equal to 2" )
       call assert_msg( 742842853, size( recv_real_array_2d, 1 ) == 2,        &
                        "MPI test failure: 2d real array recv size " //       &
                        trim( to_char( size( recv_real_array_2d ) ) )         &
                        // " not equal to 2 for dimension 1" )
       call assert_msg( 848443652,                                            &
                        recv_string_array(1) == test_string_array(1),         &
                        "MPI test failure: string array recv index 1 "        &
                        // trim( recv_string_array(1) )                       &
                        // " not equal to "                                   &
                        // trim( test_string_array(1) ) )
       call assert_msg( 553986550,                                            &
                        recv_string_array(2) == test_string_array(2),         &
                        "MPI test failure: string array recv index 2 "        &
                        // trim( recv_string_array(2) )                       &
                        // " not equal to "                                   &
                        // trim( test_string_array(2) ) )
       call assert_msg( 346596020, size( recv_real_array_2d, 2 ) == 2,        &
                        "MPI test failure: 2d real array recv size " //       &
                        trim( to_char( size( recv_real_array_2d ) ) )         &
                        // " not equal to 2 for dimension 2" )
       call assert_msg( 833103026,                                            &
                        recv_real_array_2d(1,1) == test_real_array_2d(1,1),   &
                        "MPI test failure: 2d real array recv index 1,1 "     &
                        // trim( to_char( recv_real_array_2d(1,1) ) )         &
                        // " not equal to "                                   &
                        // trim( to_char( test_real_array_2d(1,1) ) ) )
       call assert_msg( 145757864,                                            &
                        recv_real_array_2d(2,1) == test_real_array_2d(2,1),   &
                        "MPI test failure: 2d real array recv index 2,1 "     &
                        // trim( to_char( recv_real_array_2d(2,1) ) )         &
                        // " not equal to "                                   &
                        // trim( to_char( test_real_array_2d(2,1) ) ) )
       call assert_msg( 940609359,                                            &
                        recv_real_array_2d(1,2) == test_real_array_2d(1,2),   &
                        "MPI test failure: 2d real array recv index 1,2 "     &
                        // trim( to_char( recv_real_array_2d(1,2) ) )         &
                        // " not equal to "                                   &
                        // trim( to_char( test_real_array_2d(1,2) ) ) )
       call assert_msg( 770452455,                                            &
                        recv_real_array_2d(2,2) == test_real_array_2d(2,2),   &
                        "MPI test failure: 2d real array recv index 2,2 "     &
                        // trim( to_char( recv_real_array_2d(2,2) ) )         &
                        // " not equal to "                                   &
                        // trim( to_char( test_real_array_2d(2,2) ) ) )
       call assert_msg( 399792135,                                            &
                        recv_real_array_3d(1,1,1) ==                          &
                        test_real_array_3d(1,1,1),                            &
                        "MPI test failure: 3d real array recv index 1,1,1 "   &
                        // trim( to_char( recv_real_array_3d(1,1,1) ) )       &
                        // " not equal to "                                   &
                        // trim( to_char( test_real_array_3d(1,1,1) ) ) )
       call assert_msg( 229635231,                                            &
                        recv_real_array_3d(2,1,1) ==                          &
                        test_real_array_3d(2,1,1),                            &
                        "MPI test failure: 3d real array recv index 2,1,1 "   &
                        // trim( to_char( recv_real_array_3d(2,1,1) ) )       &
                        // " not equal to "                                   &
                        // trim( to_char( test_real_array_3d(2,1,1) ) ) )
       call assert_msg( 341953576,                                            &
                        recv_real_array_3d(1,2,1) ==                          &
                        test_real_array_3d(1,2,1),                            &
                        "MPI test failure: 3d real array recv index 1,2,1 "   &
                        // trim( to_char( recv_real_array_3d(1,2,1) ) )       &
                        // " not equal to "                                   &
                        // trim( to_char( test_real_array_3d(1,2,1) ) ) )
       call assert_msg( 171796672,                                            &
                        recv_real_array_3d(2,2,1) ==                          &
                        test_real_array_3d(2,2,1),                            &
                        "MPI test failure: 3d real array recv index 2,2,1 "   &
                        // trim( to_char( recv_real_array_3d(2,2,1) ) )       &
                        // " not equal to "                                   &
                        // trim( to_char( test_real_array_3d(2,2,1) ) ) )
       call assert_msg( 901639767,                                            &
                        recv_real_array_3d(1,1,2) ==                          &
                        test_real_array_3d(1,1,2),                            &
                        "MPI test failure: 3d real array recv index 1,1,2 "   &
                        // trim( to_char( recv_real_array_3d(1,1,2) ) )       &
                        // " not equal to "                                   &
                        // trim( to_char( test_real_array_3d(1,1,2) ) ) )
       call assert_msg( 449007614,                                            &
                        recv_real_array_3d(2,1,2) ==                          &
                        test_real_array_3d(2,1,2),                            &
                        "MPI test failure: 3d real array recv index 2,1,2 "   &
                        // trim( to_char( recv_real_array_3d(2,1,2) ) )       &
                        // " not equal to "                                   &
                        // trim( to_char( test_real_array_3d(2,1,2) ) ) )
       call assert_msg( 561325959,                                            &
                        recv_real_array_3d(1,2,2) ==                          &
                        test_real_array_3d(1,2,2),                            &
                        "MPI test failure: 3d real array recv index 1,2,2 "   &
                        // trim( to_char( recv_real_array_3d(1,2,2) ) )       &
                        // " not equal to "                                   &
                        // trim( to_char( test_real_array_3d(1,2,2) ) ) )
       call assert_msg( 391169055,                                            &
                        recv_real_array_3d(2,2,2) ==                          &
                        test_real_array_3d(2,2,2),                            &
                        "MPI test failure: 3d real array recv index 2,2,2 "   &
                        // trim( to_char( recv_real_array_3d(2,2,2) ) )       &
                        // " not equal to "                                   &
                        // trim( to_char( test_real_array_3d(2,2,2) ) ) )
    end if

    call musica_mpi_finalize( )

#else

    call assert( 242084546, .not. musica_mpi_support( ) )

#endif

  end subroutine test_mpi_wrappers

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program test_mpi
