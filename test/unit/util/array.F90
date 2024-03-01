! Copyright (C) 2021 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> Tests for the musica_array module

!> Tests for the musica_array module
program test_util_array

  use musica_assert,                   only : assert, almost_equal
  use musica_array
  use musica_constants,                only : musica_ik, musica_rk, musica_dk
  use musica_string,                   only : string_t
#ifdef MUSICA_USE_OPENMP
  use omp_lib
#endif

  implicit none

#ifdef MUSICA_USE_OPENMP
  write(*,*) "Testing with ", omp_get_max_threads( ), " threads"
#else
  write(*,*) "Testing without OpenMP support"
#endif

  !$omp parallel
  call test_array_functions( )
  !$omp end parallel

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Tests array functions
  subroutine test_array_functions( )

    type(string_t) :: str_array(3)
    real(kind=musica_dk), allocatable :: dbl_array(:), dbl_array_2(:)
    real(kind=musica_dk), allocatable :: merged_array(:)
    real(kind=musica_rk), allocatable :: flt_array(:)
    integer(kind=musica_ik), allocatable :: int_array(:)
    logical, allocatable :: bool_array(:)
    type(string_t) :: str
    integer(kind=musica_ik) :: idx

    allocate( dbl_array(  0 ) )
    allocate( flt_array(  0 ) )
    allocate( int_array(  0 ) )
    allocate( bool_array( 0 ) )

    ! test find_string_in_array( )
    str_array(1) = "foo"
    str_array(2) = "bar"
    str_array(3) = "foObar"

    call assert( 301097835, size( str_array ) .eq. 3 )
    call assert( 184681299, find_string_in_array( str_array, "foo", idx ) )
    call assert( 360841928, idx .eq. 1 )
    call assert( 520470218, find_string_in_array( str_array, "foObar", idx,   &
                                                  case_sensitive = .true. ) )
    call assert( 745106908, idx .eq. 3 )
    call assert( 239900503, .not. find_string_in_array( str_array, "fooBar",  &
                                              idx, case_sensitive = .true. ) )
    call assert( 234636196, .not. find_string_in_array( str_array,            &
                                                        "not there", idx ) )
    str = "bar"
    call assert( 911905039, find_string_in_array( str_array, str, idx ) )
    call assert( 689173883, idx .eq. 2 )
    str = "Bar"
    call assert( 183967478, .not. find_string_in_array( str_array, str, idx,  &
                                                   case_sensitive = .true. ) )
    str = "not there"
    call assert( 231277423, .not. find_string_in_array( str_array, str, idx ) )

    ! test find_string_in_split_array( )
    str_array( 1 ) = "foo.BaR"
    str_array( 2 ) = "Bar.foO"
    str_array( 3 ) = "justfoo"

    call assert( 100527721, find_string_in_split_array( str_array, "foo", ".",&
                                                        1, idx ) )
    call assert( 253438465, idx .eq. 1 )
    call assert( 192693428, find_string_in_split_array( str_array, "foo", ".",&
                                                        2, idx ) )
    call assert( 522478622, idx .eq. 2 )
    call assert( 634796967, .not. find_string_in_split_array( str_array,      &
                              "foo", ".", 2, idx, case_sensitive = .true. ) )
    call assert( 747115312, find_string_in_split_array( str_array, "BaR", ".",&
                              2, idx, case_sensitive = .true. ) )
    call assert( 859433657, idx .eq. 1 )
    str = "foo"
    call assert( 929884076, find_string_in_split_array( str_array, str, ".",  &
                                                        1, idx ) )
    call assert( 477251923, idx .eq. 1 )
    call assert( 924619769, find_string_in_split_array( str_array, str, ".",  &
                                                        2, idx ) )
    call assert( 471987616, idx .eq. 2 )
    call assert( 366839112, .not. find_string_in_split_array( str_array,      &
                              str, ".", 2, idx, case_sensitive = .true. ) )
    str = "BaR"
    call assert( 196682208, find_string_in_split_array( str_array, str, ".",  &
                              2, idx, case_sensitive = .true. ) )
    call assert( 926525303, idx .eq. 1 )

    ! test merge_series( )
    if( allocated( dbl_array   ) ) deallocate( dbl_array   )
    if( allocated( dbl_array_2 ) ) deallocate( dbl_array_2 )
    dbl_array   = [ 3.5_musica_dk, 5.0_musica_dk, 12.3_musica_dk ]
    dbl_array_2 = [ 1.0_musica_dk, 4.2_musica_dk, 5.0_musica_dk,              &
                    12.3_musica_dk, 24.3_musica_dk ]
    merged_array = merge_series( dbl_array, dbl_array_2 )
    call assert( 182507698, size( merged_array ) .eq. 6 )
    call assert( 105969740, merged_array( 1 ) .eq.  1.0_musica_dk )
    call assert( 835812835, merged_array( 2 ) .eq.  3.5_musica_dk )
    call assert( 948131180, merged_array( 3 ) .eq.  4.2_musica_dk )
    call assert( 777974276, merged_array( 4 ) .eq.  5.0_musica_dk )
    call assert( 607817372, merged_array( 5 ) .eq. 12.3_musica_dk )
    call assert( 155185219, merged_array( 6 ) .eq. 24.3_musica_dk )
    merged_array = merge_series( dbl_array, dbl_array_2,                      &
                                 with_bounds_from = dbl_array )
    call assert( 267503564, size( merged_array ) .eq. 4 )
    call assert( 162355060, merged_array( 1 ) .eq.  3.5_musica_dk )
    call assert( 609722906, merged_array( 2 ) .eq.  4.2_musica_dk )
    call assert( 439566002, merged_array( 3 ) .eq.  5.0_musica_dk )
    call assert( 886933848, merged_array( 4 ) .eq. 12.3_musica_dk )

    ! test calculate_linear_array( )
    if( allocated( dbl_array ) ) deallocate( dbl_array )
    dbl_array = calculate_linear_array( 1.0_musica_dk, 5.0_musica_dk, 5 )
    call assert( 781682679, size( dbl_array ) .eq. 5 )
    call assert( 106319370, dbl_array( 1 ) .eq. 1.0_musica_dk )
    call assert( 824180612, almost_equal( dbl_array( 2 ), 2.0_musica_dk ) )
    call assert( 654023708, almost_equal( dbl_array( 3 ), 3.0_musica_dk ) )
    call assert( 201391555, almost_equal( dbl_array( 4 ), 4.0_musica_dk ) )
    call assert( 996243050, dbl_array( 5 ) .eq. 5.0_musica_dk )

    ! test calculate_logarithmic_array( )
    if( allocated( dbl_array ) ) deallocate( dbl_array )
    dbl_array =                                                               &
        calculate_logarithmic_array( 1.0_musica_dk, 10000.0_musica_dk, 5 )
    call assert( 764888814, size( dbl_array ) .eq. 5 )
    call assert( 312256661, dbl_array( 1 ) .eq. 1.0_musica_dk )
    call assert( 142099757, almost_equal( dbl_array( 2 ), 10.0_musica_dk ) )
    call assert( 589467603, almost_equal( dbl_array( 3 ), 100.0_musica_dk ) )
    call assert( 136835450, almost_equal( dbl_array( 4 ), 1000.0_musica_dk ) )
    call assert( 931686945, dbl_array( 5 ) .eq. 10000.0_musica_dk )

  end subroutine test_array_functions

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program test_util_array
