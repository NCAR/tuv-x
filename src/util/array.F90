! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> The musica_array module

!> Functions for working with allocatable arrays
module musica_array

  use musica_constants,                only : musica_ik, musica_dk

  implicit none
  private

  public :: find_string_in_array, find_string_in_split_array,                 &
            merge_series, calculate_linear_array, calculate_logarithmic_array

  ! Find a string in an array of strings
  interface find_string_in_array
    module procedure :: find_string_in_array_string
    module procedure :: find_string_in_array_char
  end interface find_string_in_array

  ! Find a string in an array of split strings
  interface find_string_in_split_array
    module procedure :: find_string_in_split_array_string
    module procedure :: find_string_in_split_array_char
  end interface find_string_in_split_array

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Finds a string in a string array (case insensitive by default)
  logical function find_string_in_array_char( array, string, id,              &
      case_sensitive )

    use musica_string,                 only : string_t

    !> Array to search
    type(string_t), intent(in) :: array(:)
    !> String to search for
    character(len=*), intent(in) :: string
    !> Index of located string
    integer(kind=musica_ik), intent(out) :: id
    !> Do a case sensitive search
    logical, intent(in), optional :: case_sensitive

    type(string_t) :: temp_string, array_string
    integer :: i_str
    logical :: is_case_sensitive

    is_case_sensitive = .false.
    if( present( case_sensitive ) ) then
      is_case_sensitive = case_sensitive
    end if
    id = 0
    find_string_in_array_char = .false.
    temp_string = trim( string )
    if( .not. is_case_sensitive ) temp_string = temp_string%to_lower( )
    do i_str = 1, size( array )
      array_string = array( i_str )
      if( .not. is_case_sensitive ) array_string = array_string%to_lower( )
      if( temp_string .eq. array_string ) then
        id = i_str
        find_string_in_array_char = .true.
        exit
      end if
    end do

  end function find_string_in_array_char

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Finds a string in an array ( case insensitive by default)
  logical function find_string_in_array_string( array, string, id,            &
      case_sensitive )

    use musica_string,                 only : string_t

    !> Array to search
    type(string_t), intent(in) :: array(:)
    !> String to search for
    type(string_t), intent(in) :: string
    !> Index of located string
    integer(kind=musica_ik), intent(out) :: id
    !> Do a case sensitive search
    logical, intent(in), optional :: case_sensitive

    find_string_in_array_string = find_string_in_array_char( array,           &
        string%to_char( ), id, case_sensitive )

  end function find_string_in_array_string

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Find a string in an array of strings after splitting the array elements
  !!
  !! Case insensitive by default
  logical function find_string_in_split_array_char( array, string, splitter,  &
      element_id, id, case_sensitive )

    use musica_string,                 only : string_t

    !> Array to search
    type(string_t), intent(in) :: array(:)
    !> String to search for
    character(len=*), intent(in) :: string
    !> Splitting characters
    character(len=*), intent(in) :: splitter
    !> Element to compare in split strings
    integer(kind=musica_ik), intent(in) :: element_id
    !> Index of located string
    integer(kind=musica_ik), intent(out) :: id
    !> Do a case sensitive search
    logical, intent(in), optional :: case_sensitive

    type(string_t) :: temp_string, array_string
    type(string_t), allocatable :: split_string(:)
    integer :: i_str
    logical :: is_case_sensitive

    is_case_sensitive = .false.
    if( present( case_sensitive ) ) then
      is_case_sensitive = case_sensitive
    end if
    id = 0
    find_string_in_split_array_char = .false.
    temp_string = trim( string )
    if( .not. is_case_sensitive ) temp_string = temp_string%to_lower( )
    do i_str = 1, size( array )
      array_string = array( i_str )
      if( .not. is_case_sensitive ) array_string = array_string%to_lower( )
      split_string = array_string%split( splitter )
      if( size( split_string ) .ge. element_id ) then
        array_string = split_string( element_id )
      else
        cycle
      end if
      if( temp_string .eq. array_string ) then
        id = i_str
        find_string_in_split_array_char = .true.
        exit
      end if
    end do

  end function find_string_in_split_array_char

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Find a string in an array of strings after splitting the array elements
  !!
  !! Case insensitive by default
  logical function find_string_in_split_array_string( array, string, splitter, &
      element_id, id, case_sensitive )

    use musica_string,                 only : string_t

    !> Array to search
    type(string_t), intent(in) :: array(:)
    !> String to search for
    type(string_t), intent(in) :: string
    !> Splitting characters
    character(len=*), intent(in) :: splitter
    !> Element to compare in split strings
    integer(kind=musica_ik), intent(in) :: element_id
    !> Index of located string
    integer(kind=musica_ik), intent(out) :: id
    !> Do a case sensitive search
    logical, intent(in), optional :: case_sensitive

    find_string_in_split_array_string =                                       &
        find_string_in_split_array_char( array, string%to_char( ), splitter,  &
                                         element_id, id, case_sensitive )

  end function find_string_in_split_array_string

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Merge two sets of values into a single set without duplicates
  !!
  !! Both sets must be arranged in increasing order
  !!
  function merge_series( a, b, with_bounds_from ) result( new_set )

    !> New series
    real(kind=musica_dk), allocatable :: new_set(:)
    !> First series
    real(kind=musica_dk), intent(in) :: a(:)
    !> Second series
    real(kind=musica_dk), intent(in) :: b(:)
    !> Restricts series to bounds in this array
    real(kind=musica_dk), intent(in), optional :: with_bounds_from(:)

    real(kind=musica_dk) :: curr_val, val_a, val_b, min_val, max_val
    integer :: n_total, i_a, i_b, n_a, n_b

    if( present( with_bounds_from ) ) then
      min_val = with_bounds_from( 1 )
      max_val = with_bounds_from( size( with_bounds_from ) )
    else
      min_val = -huge( 0.0_musica_dk )
      max_val =  huge( 0.0_musica_dk )
    endif

    n_a = size( a )
    n_b = size( b )
    if( n_a + n_b .eq. 0 ) then
      allocate( new_set( 0 ) )
      return
    end if

    curr_val = huge( 1.0_musica_dk )
    if( n_a .gt. 0 ) curr_val = a( 1 )
    if( n_b .gt. 0 ) then
      if( b( 1 ) .lt. curr_val ) curr_val = b( 1 )
    end if
    if( curr_val .lt. min_val ) curr_val = min_val
    if( curr_val .gt. max_val ) curr_val = max_val

    i_a = 1
    i_b = 1
    n_total = 0
    do while( i_a .le. n_a )
      if( a( i_a ) .ge. min_val ) exit
      i_a = i_a + 1
    end do
    do while( i_b .le. n_b )
      if( b( i_b ) .ge. min_val ) exit
      i_b = i_b + 1
    end do
    do while( i_a .le. n_a .or. i_b .le. n_b )
      if( i_a .le. n_a ) then
        val_a = a( i_a )
        if( val_a .gt. max_val ) then
          i_a = n_a + 1
          cycle
        end if
      else
        val_a = huge( 1.0_musica_dk )
      end if
      if( i_b .le. n_b ) then
        val_b = b( i_b )
        if( val_b .gt. max_val ) then
          i_b = n_b + 1
          cycle
        end if
      else
        val_b = huge( 1.0_musica_dk )
      end if
      curr_val = min( val_a, val_b )
      n_total = n_total + 1
      if( val_a .le. curr_val ) i_a = i_a + 1
      if( val_b .le. curr_val ) i_b = i_b + 1
    end do

    allocate( new_set( n_total ) )

    i_a = 1
    i_b = 1
    n_total = 0
    do while( i_a .le. n_a )
      if( a( i_a ) .ge. min_val ) exit
      i_a = i_a + 1
    end do
    do while( i_b .le. n_b )
      if( b( i_b ) .ge. min_val ) exit
      i_b = i_b + 1
    end do
    do while( i_a .le. n_a .or. i_b .le. n_b )
      if( i_a .le. n_a ) then
        val_a = a( i_a )
        if( val_a .gt. max_val ) then
          i_a = n_a + 1
          cycle
        end if
      else
        val_a = huge( 1.0_musica_dk )
      end if
      if( i_b .le. n_b ) then
        val_b = b( i_b )
        if( val_b .gt. max_val ) then
          i_b = n_b + 1
          cycle
        end if
      else
        val_b = huge( 1.0_musica_dk )
      end if
      curr_val = min( val_a, val_b )
      n_total = n_total + 1
      new_set( n_total ) = curr_val
      if( val_a .le. curr_val ) i_a = i_a + 1
      if( val_b .le. curr_val ) i_b = i_b + 1
    end do

  end function merge_series

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Allocates and calculates an array of linearly increasing value with
  !! specified minimum and maximum values and number of elements
  function calculate_linear_array( minimum, maximum, number_of_elements )     &
      result( new_array )

    use musica_assert,                 only : assert

    !> Calculated array
    real(kind=musica_dk), allocatable :: new_array(:)
    !> Minimum array value
    real(kind=musica_dk), intent(in) :: minimum
    !> Maximum array value
    real(kind=musica_dk), intent(in) :: maximum
    !> Number of array elements
    integer(kind=musica_ik), intent(in) :: number_of_elements

    integer(kind=musica_ik) :: i_elem
    real(kind=musica_dk) :: space

    call assert( 167917803, maximum .gt. minimum )
    call assert( 211868975, number_of_elements .ge. 1 )
    allocate( new_array( number_of_elements ) )
    space = ( maximum - minimum ) /                                           &
            real( number_of_elements - 1, kind=musica_dk )
    new_array( 1 ) = minimum
    do i_elem = 2, number_of_elements - 1
      new_array( i_elem ) = new_array( i_elem - 1 ) + space
    end do
    new_array( number_of_elements ) = maximum

  end function calculate_linear_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Allocates and calculates an array of logarithmically increasing value
  !! with specified minimum and maximum values and number of elements
  function calculate_logarithmic_array( minimum, maximum, number_of_elements )&
      result( new_array )

    use musica_assert,                 only : assert

    !> Calculated array
    real(kind=musica_dk), allocatable :: new_array(:)
    !> Minimum array value
    real(kind=musica_dk), intent(in) :: minimum
    !> Maximum array value
    real(kind=musica_dk), intent(in) :: maximum
    !> Number of array elements
    integer(kind=musica_ik), intent(in) :: number_of_elements

    integer(kind=musica_ik) :: i_elem
    real(kind=musica_dk) :: space

    call assert( 527530853, maximum .gt. minimum )
    call assert( 752167543, number_of_elements .gt. 1 )
    allocate( new_array( number_of_elements ) )
    space = ( log( maximum ) - log( minimum ) ) /                             &
            real( number_of_elements - 1, kind=musica_dk )
    new_array( 1 ) = minimum
    do i_elem = 2, number_of_elements - 1
      new_array( i_elem ) = exp( log( new_array( i_elem - 1 ) ) + space )
    end do
    new_array( number_of_elements ) = maximum

  end function calculate_logarithmic_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module musica_array
