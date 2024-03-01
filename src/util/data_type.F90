! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> The musica_data_type module

!> The data_type_t type and related functions
module musica_data_type

  implicit none
  private

  public :: data_type_t

  !> Enum type for referencing specific primitive data-types
  type :: data_type_t
    private
    integer :: id_ = 0
  contains
    !> Equality comparisons
    !! @{
    procedure, private :: equals_data_type
    generic :: operator(==) => equals_data_type
    procedure, private :: not_equals_data_type
    generic :: operator(/=) => not_equals_data_type
    !> @}
    !> Name of the data type
    procedure :: name
  end type data_type_t

  !> Constructor of data_type_t objects
  interface data_type_t
    module procedure :: constructor_string
    module procedure :: constructor_char
  end interface

  !> @name MUSICA data types
  !! @{
  type(data_type_t), parameter, public :: kTypeUnknown = data_type_t(0)
  type(data_type_t), parameter, public :: kInteger     = data_type_t(1)
  type(data_type_t), parameter, public :: kFloat       = data_type_t(2)
  type(data_type_t), parameter, public :: kDouble      = data_type_t(3)
  type(data_type_t), parameter, public :: kBoolean     = data_type_t(4)
  !> @}

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Constructor of data_type_t objects
  function constructor_string( name ) result( new_obj )

    use musica_string,                 only : string_t

    !> New data type object
    type(data_type_t) :: new_obj
    !> Name of the data type
    type(string_t), intent(in) :: name

    new_obj = constructor_char( name%to_char( ) )

  end function constructor_string

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Constructor of data_type_t objects
  function constructor_char( name ) result( new_obj )

    use musica_string,                 only : string_t

    !> New data type object
    type(data_type_t) :: new_obj
    !> Name of the data type
    character(len=*), intent(in) :: name

    type(string_t) :: l_name

    l_name = name
    l_name = l_name%to_lower( )
    if( l_name .eq. "integer" ) then
      new_obj = kInteger
    else if( l_name .eq. "float" ) then
      new_obj = kFloat
    else if( l_name .eq. "double" ) then
      new_obj = kDouble
    else if( l_name .eq. "boolean" ) then
      new_obj = kBoolean
    end if

  end function constructor_char

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compares two data types for equality
  logical function equals_data_type( a, b )

    !> First data type
    class(data_type_t), intent(in) :: a
    !> Second data type
    class(data_type_t), intent(in) :: b

    equals_data_type = a%id_ .eq. b%id_

  end function equals_data_type

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compares two data types for inequality
  logical function not_equals_data_type( a, b )

    !> First data type
    class(data_type_t), intent(in) :: a
    !> Second data type
    class(data_type_t), intent(in) :: b

    not_equals_data_type = .not. a .eq. b

  end function not_equals_data_type

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the name of a data type
  type(string_t) function name( this )

    use musica_assert,                 only : die
    use musica_string,                 only : string_t

    !> Data type
    class(data_type_t), intent(in) :: this

    select case( this%id_ )
    case( kInteger%id_ )
      name = "integer"
    case( kFloat%id_ )
      name = "float"
    case( kDouble%id_ )
      name = "double"
    case( kBoolean%id_ )
      name = "boolean"
    case( kTypeUnknown%id_ )
      name = "unknown"
    case default
      call die( 451420098 )
    end select

  end function name

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module musica_data_type
