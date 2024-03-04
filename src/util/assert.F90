! Portions Copyright (C) 2005-2016 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.
!
! Portions Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> The musica_assert module.

!> Assertion functions
module musica_assert

  implicit none

  !> Unit for error output files
  integer, parameter :: kErrorFileId = 10
  !> Error output id
  integer, parameter :: kErrorId = 0

  interface assert_msg
    procedure :: assert_msg_string
    procedure :: assert_msg_char
  end interface

  interface assert_warn_msg
    procedure :: assert_warn_msg_string
    procedure :: assert_warn_msg_char
  end interface

  interface die_msg
    procedure :: die_msg_string
    procedure :: die_msg_char
  end interface

  interface almost_equal
    procedure :: almost_equal_complex_real
    procedure :: almost_equal_complex_double
    procedure :: almost_equal_real
    procedure :: almost_equal_double
  end interface

  interface are_equal
    procedure :: compare_arrays_1D_real
    procedure :: compare_arrays_2D_real
    procedure :: compare_arrays_1D_double
    procedure :: compare_arrays_2D_double
  end interface

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Asserts condition to be true or fails with provided message
  subroutine assert_msg_string( code, condition, error_message )

    use musica_string,                 only : string_t

    !> Unique code for the assertion
    integer, intent(in) :: code
    !> Condition to evaluate
    logical, intent(in) :: condition
    !> Message to display on failure
    type(string_t), intent(in) :: error_message

    character(len=50) :: str_code

    if( .not. condition ) then
      write(str_code,'(i30)') code
      write(kErrorId,*) "ERROR (Musica-"//trim( adjustl( str_code ) )//"): "  &
                        //error_message%val_
      open( unit = kErrorFileId, file = "error.json", action = "WRITE" )
      write(kErrorFileId,'(A)') '{'
      write(kErrorFileId,'(A)') '  "code" : "'//trim( adjustl( str_code ) )//'",'
      write(kErrorFileId,'(A)') '  "message" : "'//error_message%val_//'"'
      write(kErrorFileId,'(A)') '}'
      close(kErrorFileId)
      stop 3
    end if

  end subroutine assert_msg_string

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Asserts condition to be true or fails with provided message
  subroutine assert_msg_char( code, condition, error_message )

    !> Unique code for the assertion
    integer, intent(in) :: code
    !> Condition to evaluate
    logical, intent(in) :: condition
    !> Message to display on failure
    character(len=*), intent(in) :: error_message

    character(len=50) :: str_code

    if( .not. condition ) then
      write(str_code,'(i30)') code
      write(kErrorId,*) "ERROR (Musica-"//trim( adjustl( str_code ) )//"): "  &
                        //error_message
      open( unit = kErrorFileId, file = "error.json", action = "WRITE" )
      write(kErrorFileId,'(A)') '{'
      write(kErrorFileId,'(A)') '  "code" : "'//trim( adjustl( str_code ) )//'",'
      write(kErrorFileId,'(A)') '  "message" : "'//error_message//'"'
      write(kErrorFileId,'(A)') '}'
      close(kErrorFileId)
      stop 3
    end if

  end subroutine assert_msg_char

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Asserts condition to be true or fails
  subroutine assert( code, condition )

    !> Unique code for the assertion
    integer, intent(in) :: code
    !> Condition to evaluate
    logical, intent(in) :: condition

    call assert_msg( code, condition, 'assertion failed' )

  end subroutine assert

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Asserts condition to be true or prints a provided warning message
  subroutine assert_warn_msg_string( code, condition, warning_message )

    use musica_string,                 only : string_t

    !> Unique code for the assertion
    integer, intent(in) :: code
    !> Condition to evaluate
    logical, intent(in) :: condition
    !> Message to display on failure
    type(string_t), intent(in) :: warning_message

    character(len=50) :: str_code

    if( .not. condition ) then
      write(str_code,'(i30)') code
      write(kErrorId,*) "WARNING (Musica-"//trim( adjustl( str_code ) )//     &
                        "): "//warning_message%val_
    end if

  end subroutine assert_warn_msg_string

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Asserts condition to be true or prints a provided warning message
  subroutine assert_warn_msg_char( code, condition, warning_message )

    !> Unique code for the assertion
    integer, intent(in) :: code
    !> Condition to evaluate
    logical, intent(in) :: condition
    !> Message to display on failure
    character(len=*), intent(in) :: warning_message

    character(len=50) :: str_code

    if( .not. condition ) then
      write(str_code,'(i30)') code
      write(kErrorId,*) "WARNING (Musica-"//trim( adjustl( str_code ) )//     &
                        "): "//warning_message
    end if

  end subroutine assert_warn_msg_char

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Errors immediately and prints a provided message
  subroutine die_msg_string( code, error_message )

    use musica_string,                 only : string_t

    !> Unique code for the failure
    integer, intent(in) :: code
    !> Message to display with failure
    type(string_t), intent(in) :: error_message

    call assert_msg( code, .false., error_message )

  end subroutine die_msg_string

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Errors immediately and prints a provided message
  subroutine die_msg_char( code, error_message )

    !> Unique code for the failure
    integer, intent(in) :: code
    !> Message to display with failure
    character(len=*), intent(in) :: error_message

    call assert_msg( code, .false., error_message )

  end subroutine die_msg_char

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Errors immediately
  subroutine die( code )

    !> Unique code for the failure
    integer, intent(in) :: code

    call die_msg( code, "Internal error" )

  end subroutine die

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determines whether two real numbers are equal within a provided or
  !! standard tolerance
  logical function almost_equal_real( a, b, relative_tolerance,               &
      absolute_tolerance ) result( almost_equal )

    use musica_constants,              only : musica_rk

    !> First number to compare
    real(kind=musica_rk), intent(in) :: a
    !> Second number to compare
    real(kind=musica_rk), intent(in) :: b
    !> Relative tolerance
    real(kind=musica_rk), intent(in), optional :: relative_tolerance
    !> Absolute tolerance
    real(kind=musica_rk), intent(in), optional :: absolute_tolerance

    real(kind=musica_rk) :: rel_tol, abs_tol

    rel_tol = 1.0e-10_musica_rk
    abs_tol = 1.0e-30_musica_rk
    if( present( relative_tolerance ) ) rel_tol = relative_tolerance
    if( present( absolute_tolerance ) ) abs_tol = absolute_tolerance

    almost_equal = .false.
    if( a .eq. b ) then
      almost_equal = .true.
    else
      if( 2.0_musica_rk * abs( a - b ) / ( abs( a ) + abs( b ) )              &
          .lt. rel_tol ) then
        almost_equal = .true.
      else if( abs( a - b ) .le. abs_tol ) then
        almost_equal = .true.
      end if
    end if

  end function almost_equal_real

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determines whether two real numbers are equal within a provided or
  !! standard tolerance
  logical function almost_equal_double( a, b, relative_tolerance,             &
      absolute_tolerance ) result( almost_equal )

    use musica_constants,              only : musica_dk

    !> First number to compare
    real(kind=musica_dk), intent(in) :: a
    !> Second number to compare
    real(kind=musica_dk), intent(in) :: b
    !> Relative tolerance
    real(kind=musica_dk), intent(in), optional :: relative_tolerance
    !> Absolute tolerance
    real(kind=musica_dk), intent(in), optional :: absolute_tolerance

    real(kind=musica_dk) :: rel_tol, abs_tol

    rel_tol = 1.0e-10_musica_dk
    abs_tol = 1.0e-30_musica_dk
    if( present( relative_tolerance ) ) rel_tol = relative_tolerance
    if( present( absolute_tolerance ) ) abs_tol = absolute_tolerance

    almost_equal = .false.
    if( a .eq. b ) then
      almost_equal = .true.
    else
      if( 2.0_musica_dk * dabs( a - b ) / ( dabs( a ) + dabs( b ) )           &
          .lt. rel_tol ) then
        almost_equal = .true.
      else if( dabs( a - b ) .le. abs_tol ) then
        almost_equal = .true.
      end if
    end if

  end function almost_equal_double

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determines whether two complex numbers are equal within a provided or
  !! standard tolerance
  logical function almost_equal_complex_real( a, b, relative_tolerance,       &
      absolute_tolerance ) result( almost_equal )

    use musica_constants,              only : musica_rk

    !> First number to compare
    complex(kind=musica_rk), intent(in) :: a
    !> Second number to compare
    complex(kind=musica_rk), intent(in) :: b
    !> Relative tolerance
    real(kind=musica_rk), intent(in), optional :: relative_tolerance
    !> Absolute tolerance
    real(kind=musica_rk), intent(in), optional :: absolute_tolerance

    real(kind=musica_rk) :: rel_tol, abs_tol

    rel_tol = 1.0e-10_musica_rk
    abs_tol = 1.0e-30_musica_rk
    if( present( relative_tolerance ) ) rel_tol = relative_tolerance
    if( present( absolute_tolerance ) ) abs_tol = absolute_tolerance

    almost_equal = .false.
    if( a .eq. b ) then
      almost_equal = .true.
    else
    associate( ra => real( a ), ia => aimag( a ),                             &
               rb => real( b ), ib => aimag( b ) )
      if( 2.0_musica_rk * abs( ra - rb ) / ( abs( ra ) + abs( rb ) )          &
          .lt. rel_tol .and.                                                  &
          2.0_musica_rk * abs( ia - ib ) / ( abs( ia ) + abs( ib ) )          &
          .lt. rel_tol ) then
        almost_equal = .true.
      else if( abs( ra - rb ) .le. abs_tol .and.                              &
               abs( ia - ib ) .le. abs_tol ) then
        almost_equal = .true.
      end if
    end associate
    end if

  end function almost_equal_complex_real

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determines whether two complex numbers are equal within a provided or
  !! standard tolerance
  logical function almost_equal_complex_double( a, b, relative_tolerance,     &
      absolute_tolerance ) result( almost_equal )

    use musica_constants,              only : musica_dk

    !> First number to compare
    complex(kind=musica_dk), intent(in) :: a
    !> Second number to compare
    complex(kind=musica_dk), intent(in) :: b
    !> Relative tolerance
    real(kind=musica_dk), intent(in), optional :: relative_tolerance
    !> Absolute tolerance
    real(kind=musica_dk), intent(in), optional :: absolute_tolerance

    real(kind=musica_dk) :: rel_tol, abs_tol

    rel_tol = 1.0e-10_musica_dk
    abs_tol = 1.0e-30_musica_dk
    if( present( relative_tolerance ) ) rel_tol = relative_tolerance
    if( present( absolute_tolerance ) ) abs_tol = absolute_tolerance

    almost_equal = .false.
    if( a .eq. b ) then
      almost_equal = .true.
    else
    associate( ra => real( a, kind=musica_dk ), ia => aimag( a ),             &
               rb => real( b, kind=musica_dk ), ib => aimag( b ) )
      if( 2.0_musica_dk * dabs( ra - rb ) / ( dabs( ra ) + dabs( rb ) )       &
          .lt. rel_tol .and.                                                  &
          2.0_musica_dk * dabs( ia - ib ) / ( dabs( ia ) + dabs( ib ) )       &
          .lt. rel_tol ) then
        almost_equal = .true.
      else if( dabs( ra - rb ) .le. abs_tol .and.                             &
               dabs( ia - ib ) .le. abs_tol ) then
        almost_equal = .true.
      end if
    end associate
    end if

  end function almost_equal_complex_double

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compares two 1D arrays for equality
  logical pure function compare_arrays_1D_real( a, b ) result( equal )

    use musica_constants,              only : musica_rk

    !> First array to compare
    real(kind=musica_rk), intent(in) :: a(:)
    !> Second array to compare
    real(kind=musica_rk), intent(in) :: b(:)

    integer :: i_elem

    equal = .false.
    if( size( a ) .ne. size( b ) ) return
    do i_elem = 1, size( a )
      if( a( i_elem ) .ne. b( i_elem ) ) return
    end do
    equal = .true.

  end function compare_arrays_1D_real

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compares two 2D arrays for equality
  logical pure function compare_arrays_2D_real( a, b ) result( equal )

    use musica_constants,              only : musica_rk

    !> First array to compare
    real(kind=musica_rk), intent(in) :: a(:,:)
    !> Second array to compare
    real(kind=musica_rk), intent(in) :: b(:,:)

    integer :: i_elem

    equal = .false.
    if( size( a, 1 ) .ne. size( b, 1 ) ) return
    if( size( a, 2 ) .ne. size( b, 2 ) ) return
    do i_elem = 1, size( a, 1 )
      if( .not. compare_arrays_1D_real( a(:,i_elem), b(:,i_elem) ) ) return
    end do
    equal = .true.

  end function compare_arrays_2D_real

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compares two 1D arrays for equality
  logical pure function compare_arrays_1D_double( a, b ) result( equal )

    use musica_constants,              only : musica_dk

    !> First array to compare
    real(kind=musica_dk), intent(in) :: a(:)
    !> Second array to compare
    real(kind=musica_dk), intent(in) :: b(:)

    integer :: i_elem

    equal = .false.
    if( size( a ) .ne. size( b ) ) return
    do i_elem = 1, size( a )
      if( a( i_elem ) .ne. b( i_elem ) ) return
    end do
    equal = .true.

  end function compare_arrays_1D_double

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compares two 2D arrays for equality
  logical pure function compare_arrays_2D_double( a, b ) result( equal )

    use musica_constants,              only : musica_dk

    !> First array to compare
    real(kind=musica_dk), intent(in) :: a(:,:)
    !> Second array to compare
    real(kind=musica_dk), intent(in) :: b(:,:)

    integer :: i_elem

    equal = .false.
    if( size( a, 1 ) .ne. size( b, 1 ) ) return
    if( size( a, 2 ) .ne. size( b, 2 ) ) return
    do i_elem = 1, size( a, 1 )
      if( .not. compare_arrays_1D_double( a(:,i_elem), b(:,i_elem) ) ) return
    end do
    equal = .true.

  end function compare_arrays_2D_double

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module musica_assert
