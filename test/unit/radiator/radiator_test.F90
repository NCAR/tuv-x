! Copyright (C) 2021 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> Tests the photo_decomp radiator module

!> Test module for the radiator_core_t type
program radiator_test

  use musica_string,                   only : string_t
  use musica_mpi,                      only : musica_mpi_init,                &
                                              musica_mpi_finalize
  use radiator_core,                   only : radiator_core_t
#ifdef MUSICA_USE_OPENMP
  use omp_lib
#endif

  implicit none

  type :: core_ptr
    class(radiator_core_t), pointer :: core_
  end type core_ptr
  type(core_ptr), allocatable :: ptrs(:)

  ! Command-line options
  character(len=256) :: argument
  type(string_t)     :: configFileSpec
  integer            :: i_thread

  call musica_mpi_init( )

  ! Get the model configuration file and options from the command line
  argument = 'test/data/radiator.test.config.json'

  configFileSpec = argument

  ! instatiate and initialize radiator core object
#ifdef MUSICA_USE_OPENMP
  write(*,*) "Testing with ", omp_get_max_threads( ), " threads"
  allocate( ptrs( omp_get_max_threads( ) ) )
#else
  write(*,*) "Testing without OpenMP support"
  allocate( ptrs( 1 ) )
#endif

  do i_thread = 1, size( ptrs )
    ptrs( i_thread )%core_ => radiator_core_t( configFileSpec )
  end do

  ! set radiator cross sections
  ! \todo redo core so it can be used to test run functions with OpenMP
  associate( core => ptrs( 1 )%core_ )
    call core%test( )
  end associate

  do i_thread = 1, size( ptrs )
    deallocate( ptrs( i_thread )%core_ )
  end do
  deallocate( ptrs )

  ! test radiator_state_t functions
  !$omp parallel
  call test_state( )
  !$omp end parallel

  call musica_mpi_finalize( )

  contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine test_state( )
    ! Perform tests on the radiator_state_t functions

    use musica_assert,                 only : almost_equal, assert
    use musica_constants,              only : dk => musica_dk
    use tuvx_radiator,                 only : radiator_ptr, radiator_state_t

    type(radiator_state_t) :: a_state, b_state, total_state
    type(radiator_ptr)     :: radiators(2), one_radiator(1)

    allocate( radiators(1)%val_ )
    allocate( radiators(2)%val_ )
    allocate( one_radiator(1)%val_ )

    a_state%layer_OD_  = reshape( (/ 0.1_dk, 0.2_dk, 0.3_dk, 0.4_dk /),       &
                                  (/ 2, 2 /) )
    b_state%layer_OD_  = reshape( (/ 0.4_dk, 0.3_dk, 0.3_dk, 0.4_dk /),       &
                                  (/ 2, 2 /) )
    a_state%layer_SSA_ = reshape( (/ 0.05_dk, 0.03_dk, 0.02_dk, 0.07_dk /),   &
                                  (/ 2, 2 /) )
    b_state%layer_SSA_ = reshape( (/ 0.04_dk, 0.01_dk, 0.02_dk, 0.03_dk /),   &
                                  (/ 2, 2 /) )
    a_state%layer_G_   = reshape( (/ 0.2_dk, 0.4_dk, 0.2_dk, 0.3_dk /),       &
                                  (/ 2, 2, 1 /) )
    b_state%layer_G_   = reshape( (/ 0.4_dk, 0.2_dk, 0.1_dk, 0.5_dk /),       &
                                  (/ 2, 2, 1 /) )

    ! test accumulating a single radiator
    one_radiator(1)%val_%state_ = a_state

    allocate( total_state%layer_G_( 2, 2, 1 ) )
    call total_state%accumulate( one_radiator )

    call assert( 413535517, almost_equal( a_state%layer_OD_(1,1),             &
                                          total_state%layer_OD_(1,1) ) )
    call assert( 509059997, almost_equal( a_state%layer_OD_(1,2),             &
                                          total_state%layer_OD_(1,2) ) )
    call assert( 338903093, almost_equal( a_state%layer_OD_(2,1),             &
                                          total_state%layer_OD_(2,1) ) )
    call assert( 233754589, almost_equal( a_state%layer_OD_(2,2),             &
                                          total_state%layer_OD_(2,2) ) )
    call assert( 681122435, almost_equal( a_state%layer_SSA_(1,1),            &
                                          total_state%layer_SSA_(1,1) ) )
    call assert( 228490282, almost_equal( a_state%layer_SSA_(1,2),            &
                                          total_state%layer_SSA_(1,2) ) )
    call assert( 675858128, almost_equal( a_state%layer_SSA_(2,1),            &
                                          total_state%layer_SSA_(2,1) ) )
    call assert( 288234375, almost_equal( a_state%layer_SSA_(2,2),            &
                                          total_state%layer_SSA_(2,2) ) )
    call assert( 735602221, almost_equal( a_state%layer_G_(1,1,1),            &
                                          total_state%layer_G_(1,1,1) ) )
    call assert( 282970068, almost_equal( a_state%layer_G_(1,2,1),            &
                                          total_state%layer_G_(1,2,1) ) )
    call assert( 447862665, almost_equal( a_state%layer_G_(2,1,1),            &
                                          total_state%layer_G_(2,1,1) ) )
    call assert( 277705761, almost_equal( a_state%layer_G_(2,2,1),            &
                                          total_state%layer_G_(2,2,1) ) )

    radiators(1)%val_%state_ = a_state
    radiators(2)%val_%state_ = b_state

    call total_state%accumulate( radiators )

    call assert( 330018487, almost_equal( total_state%layer_OD_(2,1),         &
                 ( 0.2_dk * 0.03_dk + 0.3_dk * 0.01_dk ) +                    &
                 ( 0.2_dk * ( 1.0_dk - 0.03_dk ) +                            &
                   0.3_dk * ( 1.0_dk - 0.01_dk ) ) ) )

    call assert( 728457977, almost_equal( total_state%layer_SSA_(2,1),        &
                 1.0_dk / ( 1.0_dk +                                          &
                           ( 0.2_dk * ( 1.0_dk - 0.03_dk ) +                  &
                             0.3_dk * ( 1.0_dk - 0.01_dk ) ) /                &
                           ( 0.2_dk * 0.03_dk + 0.3_dk * 0.01_dk ) ) ) )

    call assert( 802376580, almost_equal( total_state%layer_G_(2,1,1),        &
                 ( 0.2_dk * 0.03_dk * 0.4_dk + 0.3_dk * 0.01_dk * 0.2_dk ) /  &
                 ( 0.2_dk * 0.03_dk + 0.3_dk * 0.01_dk ) ) )

    deallocate( radiators(1)%val_    )
    deallocate( radiators(2)%val_    )
    deallocate( one_radiator(1)%val_ )

  end subroutine test_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Fail run and print usage info
  subroutine fail_run( )

    write(*,*) "Usage: ./radiator_test configuration_file.json"
    stop 3

  end subroutine fail_run

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program radiator_test
