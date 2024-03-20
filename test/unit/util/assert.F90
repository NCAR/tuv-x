! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> Tests for the musica_assert module

!> Test module for the musica_assert module
program test_util_assert

  use musica_assert
#ifdef MUSICA_USE_OPENMP
  use omp_lib
#endif

  implicit none

  character(len=256) :: failure_test_type

  if( command_argument_count( ) .eq. 0 ) then
#ifdef MUSICA_USE_OPENMP
    write(*,*) "Testing with ", omp_get_max_threads( ), " threads"
#else
    write(*,*) "Testing without OpenMP support"
#endif
    !$omp parallel
    call test_assert( )
    !$omp end parallel
  else if( command_argument_count( ) .eq. 1 ) then
    call get_command_argument( 1, failure_test_type )
    call failure_test( failure_test_type )
  else
    call die( 233227610 )
  end if

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Test assert functions
  subroutine test_assert( )

    use musica_constants,              only : rk => musica_rk, dk => musica_dk
    use musica_string,                 only : string_t

    type(string_t) :: str
    real(kind=dk) :: a1d(3), b1d(3), c1d(2), a2d(3,3), b2d(3,3), c2d(2,3)
    real(kind=rk) :: a1r(3), b1r(3), c1r(2), a2r(3,3), b2r(3,3), c2r(2,3)

    str = "foo"

    call assert_msg( 449241220, .true., "foo" )
    call assert_msg( 612680578, .true., str )
    call assert( 549577712, .true. )

    ! test almost_equal( )
    ! for real
    call assert( 126460695, almost_equal( 12.5_rk, 12.5_rk ) )
    call assert( 740626672, .not. almost_equal( 12.5_rk, 12.6_rk ) )
    call assert( 172317401, almost_equal( 12.5_rk, 12.6_rk,                   &
                                          relative_tolerance = 0.11_rk ) )
    call assert( 955187043, almost_equal( 12.5_rk, 12.6_rk,                   &
                                          absolute_tolerance = 0.11_rk ) )
    call assert( 293998244, .not.                                             &
                 almost_equal( 12.5e34_rk,                                    &
                               12.5e34_rk + 12.5e34_rk * 1.0e-5_rk ) )
    call assert( 881294037,                                                   &
                 almost_equal( 12.5e-34_rk, 12.5e-34_rk + 1.0e-32_rk ) )
    call assert( 151450942, .not.                                             &
                 almost_equal( 12.5e34_rk,                                    &
                               12.5e34_rk - 12.5e34_rk * 1.0e-5_rk ) )
    call assert( 328325392,                                                   &
                 almost_equal( 12.5e-34_rk, 12.5e-34_rk - 1.0e-32_rk ) )
    call assert( 597365549,                                                   &
                 almost_equal( 12.5e34_rk,                                    &
                               12.5e34_rk + 12.5e34_rk * 1.0e-5_rk,           &
                               relative_tolerance = 1.0e-4_rk ) )
    call assert( 709683894, .not.                                             &
                 almost_equal( 12.5e-34_rk, 12.5e-34_rk + 1.0e-32_rk,         &
                               absolute_tolerance = 1.0e-33_rk ) )
    call assert( 539526990,                                                   &
                 almost_equal( 12.5e34_rk,                                    &
                               12.5e34_rk - 12.5e34_rk * 1.0e-5_rk,           &
                               relative_tolerance = 1.0e-4_rk ) )
    call assert( 986894836, .not.                                             &
                 almost_equal( 12.5e-34_rk, 12.5e-34_rk - 1.0e-32_rk,         &
                               absolute_tolerance = 1.0e-33_rk ) )

    ! for double
    call assert( 799568563, almost_equal( 12.5_dk, 12.5_dk ) )
    call assert( 794304256, .not. almost_equal( 12.5_dk, 12.6_dk ) )
    call assert( 341672103, almost_equal( 12.5_dk, 12.6_dk,                   &
                                          relative_tolerance = 0.11_dk ) )
    call assert( 236523599, almost_equal( 12.5_dk, 12.6_dk,                   &
                                          absolute_tolerance = 0.11_dk ) )
    call assert( 966366694, .not.                                             &
                 almost_equal( 12.5e94_dk,                                    &
                               12.5e94_dk + 12.5e94_dk * 1.0e-5_dk ) )
    call assert( 796209790,                                                   &
                 almost_equal( 12.5e-94_dk, 12.5e-94_dk + 1.0e-92_dk ) )
    call assert( 343577637, .not.                                             &
                 almost_equal( 12.5e94_dk,                                    &
                               12.5e94_dk - 12.5e94_dk * 1.0e-5_dk ) )
    call assert( 173420733,                                                   &
                 almost_equal( 12.5e-94_dk, 12.5e-94_dk - 1.0e-92_dk ) )
    call assert( 903263828,                                                   &
                 almost_equal( 12.5e94_dk,                                    &
                               12.5e94_dk + 12.5e94_dk * 1.0e-5_dk,           &
                               relative_tolerance = 1.0e-4_dk ) )
    call assert( 733106924, .not.                                             &
                 almost_equal( 12.5e-94_dk, 12.5e-94_dk + 1.0e-92_dk,         &
                               absolute_tolerance = 1.0e-93_dk ) )
    call assert( 562950020,                                                   &
                 almost_equal( 12.5e94_dk,                                    &
                               12.5e94_dk - 12.5e94_dk * 1.0e-5_dk,           &
                               relative_tolerance = 1.0e-4_dk ) )
    call assert( 675268365, .not.                                             &
                 almost_equal( 12.5e-94_dk, 12.5e-94_dk - 1.0e-92_dk,         &
                               absolute_tolerance = 1.0e-93_dk ) )

    ! for cmplx real
    call assert( 677913317, almost_equal( cmplx( 12.5_rk, 3.2_rk, kind=rk ),  &
                                          cmplx( 12.5_rk, 3.2_rk, kind=rk ) ) )
    call assert( 837993902, .not.                                             &
                            almost_equal( cmplx( 12.5_rk, 3.2_rk, kind=rk ),  &
                                          cmplx( 12.6_rk, 3.2_rk, kind=rk ) ) )
    call assert( 264420324, .not.                                             &
                            almost_equal( cmplx( 12.5_rk, 3.2_rk, kind=rk ),  &
                                          cmplx( 12.6_rk, 3.3_rk, kind=rk ) ) )
    call assert( 827917583, almost_equal( cmplx( 12.5_rk, 3.2_rk, kind=rk ),  &
                                          cmplx( 12.6_rk, 3.2_rk, kind=rk ),  &
                                          relative_tolerance = 0.11_rk ) )
    call assert( 538724788, almost_equal( cmplx( 12.5_rk, 3.2_rk, kind=rk ),  &
                                          cmplx( 12.5_rk, 3.3_rk, kind=rk ),  &
                                          relative_tolerance = 0.11_rk ) )
    call assert( 754738398, almost_equal( cmplx( 12.5_rk, 3.2_rk, kind=rk ),  &
                                          cmplx( 12.6_rk, 3.2_rk, kind=rk ),  &
                                          absolute_tolerance = 0.11_rk ) )
    call assert( 584581494, almost_equal( cmplx( 12.5_rk, 3.2_rk, kind=rk ),  &
                                          cmplx( 12.5_rk, 3.3_rk, kind=rk ),  &
                                          absolute_tolerance = 0.11_rk ) )

    ! for cmplx double
    call assert( 556258071, almost_equal( cmplx( 12.5_dk, 3.2_dk, kind=dk ),  &
                                          cmplx( 12.5_dk, 3.2_dk, kind=dk ) ) )
    call assert( 268518515, .not.                                             &
                            almost_equal( cmplx( 12.5_dk, 3.2_dk, kind=dk ),  &
                                          cmplx( 12.6_dk, 3.2_dk, kind=dk ) ) )
    call assert( 163370011, .not.                                             &
                            almost_equal( cmplx( 12.5_dk, 3.2_dk, kind=dk ),  &
                                          cmplx( 12.6_dk, 3.3_dk, kind=dk ) ) )
    call assert( 610737857, almost_equal( cmplx( 12.5_dk, 3.2_dk, kind=dk ),  &
                                          cmplx( 12.6_dk, 3.2_dk, kind=dk ),  &
                                          relative_tolerance = 0.11_dk ) )
    call assert( 440580953, almost_equal( cmplx( 12.5_dk, 3.2_dk, kind=dk ),  &
                                          cmplx( 12.5_dk, 3.3_dk, kind=dk ),  &
                                          relative_tolerance = 0.11_dk ) )
    call assert( 270424049, almost_equal( cmplx( 12.5_dk, 3.2_dk, kind=dk ),  &
                                          cmplx( 12.6_dk, 3.2_dk, kind=dk ),  &
                                          absolute_tolerance = 0.11_dk ) )
    call assert( 100267145, almost_equal( cmplx( 12.5_dk, 3.2_dk, kind=dk ),  &
                                          cmplx( 12.5_dk, 3.3_dk, kind=dk ),  &
                                          absolute_tolerance = 0.11_dk ) )

    ! test are_equal( )
    ! for 1d real arrays
    a1r = [ 2.3_rk, 4.2_rk, 5.2_rk ]
    b1r = [ 2.3_rk, 4.2_rk, 5.2_rk ]
    c1r = [ 2.3_rk, 4.2_rk ]
    call assert( 197244864, are_equal( a1r, b1r ) )
    call assert( 316733050, .not. are_equal( a1r, c1r ) )
    b1r(3) = 42.5_rk
    call assert( 478266874, .not. are_equal( a1r, b1r ) )

    ! for 1d double arrays
    a1d = [ 2.3_dk, 4.2_dk, 5.2_dk ]
    b1d = [ 2.3_dk, 4.2_dk, 5.2_dk ]
    c1d = [ 2.3_dk, 4.2_dk ]
    call assert( 197244864, are_equal( a1d, b1d ) )
    call assert( 316733050, .not. are_equal( a1d, c1d ) )
    b1d(3) = 42.5_dk
    call assert( 478266874, .not. are_equal( a1d, b1d ) )

    ! for 2d real arrays
    a2r(1,:) = [ 2.3_rk, 4.2_rk,   5.2_rk ]
    a2r(2,:) = [ 5.2_rk, 3.2_rk, -42.3_rk ]
    a2r(3,:) = [ 7.3_rk, 1.2_rk, 423.1_rk ]
    b2r(1,:) = [ 2.3_rk, 4.2_rk,   5.2_rk ]
    b2r(2,:) = [ 5.2_rk, 3.2_rk, -42.3_rk ]
    b2r(3,:) = [ 7.3_rk, 1.2_rk, 423.1_rk ]
    c2r(1,:) = [ 2.3_rk, 4.2_rk,   5.2_rk ]
    c2r(2,:) = [ 5.2_rk, 3.2_rk, -42.3_rk ]
    call assert( 787185609, are_equal( a2r, b2r ) )
    call assert( 341723297, .not. are_equal( a2r, c2r ) )
    b2r(3,3) = 94.2_rk
    call assert( 613669932, .not. are_equal( a2r, b2r ) )

    ! for 2d double arrays
    a2d(1,:) = [ 2.3_dk, 4.2_dk,   5.2_dk ]
    a2d(2,:) = [ 5.2_dk, 3.2_dk, -42.3_dk ]
    a2d(3,:) = [ 7.3_dk, 1.2_dk, 423.1_dk ]
    b2d(1,:) = [ 2.3_dk, 4.2_dk,   5.2_dk ]
    b2d(2,:) = [ 5.2_dk, 3.2_dk, -42.3_dk ]
    b2d(3,:) = [ 7.3_dk, 1.2_dk, 423.1_dk ]
    c2d(1,:) = [ 2.3_dk, 4.2_dk,   5.2_dk ]
    c2d(2,:) = [ 5.2_dk, 3.2_dk, -42.3_dk ]
    call assert( 787185609, are_equal( a2d, b2d ) )
    call assert( 341723297, .not. are_equal( a2d, c2d ) )
    b2d(3,3) = 94.2_dk
    call assert( 613669932, .not. are_equal( a2d, b2d ) )

  end subroutine test_assert

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Failure tests for assert functions
  subroutine failure_test( test_type )

    character(len=*), intent(in) :: test_type

    if( test_type .eq. "903602145" ) then
      call failure_test_903602145( )
    else if( test_type .eq. "151700878" ) then
      call failure_test_151700878( )
    else
      call die( 634624772 )
    end if

  end subroutine failure_test

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Test failure of assert_msg with string
  subroutine failure_test_903602145( )

    use musica_string,                 only : string_t

    type(string_t) :: msg

    msg = "foo"
    call assert_msg( 903602145, .false., msg )

  end subroutine failure_test_903602145

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Test failure of assert_msg with char array
  subroutine failure_test_151700878( )

    call assert_msg( 151700878, .false., "bar" )

  end subroutine failure_test_151700878

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program test_util_assert
