! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> The musica_target module

!> The target_t type and related functions
module musica_target

  implicit none
  private

  public :: target_t

  !> Target subset
  !!
  !! Model components can extend this type to define specific targets for
  !! property_t, iterator_t, and other types
  !!
  type, abstract :: target_t
  contains
    !> @name Equality comparisons
    !! @{
    procedure(equals_target), deferred :: equals_target
    generic :: operator(==) => equals_target
    procedure, private :: not_equals_target
    generic :: operator(/=) => not_equals_target
    !> @}
    !> Name of the target
    procedure(target_name), deferred :: name
  end type target_t

interface

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compares two targets for equality
  logical function equals_target( a, b )
    import target_t
    !> First target
    class(target_t), intent(in) :: a
    !> Second target
    class(target_t), intent(in) :: b
  end function equals_target

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the name of the target
  function target_name( this )
    use musica_string,                 only : string_t
    import target_t
    !> Target name
    type(string_t) :: target_name
    !> Target
    class(target_t), intent(in) :: this
  end function target_name

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end interface

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compares two targets for inequality
  logical function not_equals_target( a, b )

    !> First target
    class(target_t), intent(in) :: a
    !> Second target
    class(target_t), intent(in) :: b

    not_equals_target = .not. a .eq. b

  end function not_equals_target

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module musica_target
