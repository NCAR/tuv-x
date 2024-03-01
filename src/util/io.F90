! Copyright (C) 2021 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> The musica_io module

!> The io_t type and related functions
module musica_io

  implicit none
  private

  public :: io_t

  !> General input/output class
  type, abstract :: io_t
  contains
    !> @name Data read functions
    !! @{
    procedure(read_0D_double), deferred :: read_0D_double
    procedure(read_1D_double), deferred :: read_1D_double
    procedure(read_2D_double), deferred :: read_2D_double
    procedure(read_3D_double), deferred :: read_3D_double
    procedure(read_4D_double), deferred :: read_4D_double
    procedure(read_0D_int),    deferred :: read_0D_int
    procedure(read_1D_int),    deferred :: read_1D_int
    generic :: read => read_0D_double, read_1D_double, read_2D_double,        &
                       read_3D_double, read_4D_double, read_0D_int,           &
                       read_1D_int
    !> @}
    !> @name Data write functions
    !! @{
    procedure(write_0D_double), deferred :: write_0D_double
    procedure(write_1D_double), deferred :: write_1D_double
    procedure(write_2D_double), deferred :: write_2D_double
    procedure(write_3D_double), deferred :: write_3D_double
    procedure(write_4D_double), deferred :: write_4D_double
    procedure(write_0D_int),    deferred :: write_0D_int
    procedure(write_1D_int),    deferred :: write_1D_int
    generic :: write => write_0D_double, write_1D_double, write_2D_double,    &
                        write_3D_double, write_4D_double, write_0D_int,       &
                        write_1D_int
    !> @}
    !> @name Data append functions
    !! @{
    procedure(append_0D_double), deferred :: append_0D_double
    procedure(append_1D_double), deferred :: append_1D_double
    procedure(append_2D_double), deferred :: append_2D_double
    procedure(append_3D_double), deferred :: append_3D_double
    procedure(append_0D_int),    deferred :: append_0D_int
    generic :: append => append_0D_double, append_1D_double, append_2D_double,&
                         append_3D_double, append_0D_int
    !> @}
    !> Returns whether a variable exists in the file
    !! @{
    procedure(exists_char),         deferred :: exists_char
    procedure(exists_string),       deferred :: exists_string
    generic :: exists => exists_char, exists_string
    !> @}
    !> Returns the dimension names for a given variable
    procedure(variable_dimensions), deferred :: variable_dimensions
    !> Returns the units for a given variable
    procedure(variable_units),      deferred :: variable_units
    !> Sets the units for a given variable
    procedure(set_variable_units),  deferred :: set_variable_units
  end type io_t

interface

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Reads 0D double-precision floating-point data
  subroutine read_0D_double( this, variable_name, container, requestor_name )
    use musica_constants,              only : musica_dk
    use musica_string,                 only : string_t
    import io_t
    class(io_t),          intent(inout) :: this
    class(string_t),      intent(in)    :: variable_name
    real(kind=musica_dk), intent(out)   :: container
    character(len=*),     intent(in)    :: requestor_name
  end subroutine read_0D_double

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Reads 1D double-precision floating-point data
  !!
  !! If \c container is unallocated, it will be allocated to the dimensions
  !! of the read variable. Otherwise, its dimensions must match those of the
  !! read variable.
  !!
  subroutine read_1D_double( this, variable_name, container, requestor_name )
    use musica_constants,              only : musica_dk
    use musica_string,                 only : string_t
    import io_t
    class(io_t),                       intent(inout) :: this
    class(string_t),                   intent(in)    :: variable_name
    real(kind=musica_dk), allocatable, intent(inout) :: container(:)
    character(len=*),                  intent(in)    :: requestor_name
  end subroutine read_1D_double

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Reads 2D double-precision floating-point data
  !!
  !! If \c container is unallocated, it will be allocated to the dimensions
  !! of the read variable. Otherwise, its dimensions must match those of the
  !! read variable.
  !!
  subroutine read_2D_double( this, variable_name, container, requestor_name )
    use musica_constants,              only : musica_dk
    use musica_string,                 only : string_t
    import io_t
    class(io_t),                       intent(inout) :: this
    class(string_t),                   intent(in)    :: variable_name
    real(kind=musica_dk), allocatable, intent(inout) :: container(:,:)
    character(len=*),                  intent(in)    :: requestor_name
  end subroutine read_2D_double

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Reads 3D double-precision floating-point data
  !!
  !! If \c container is unallocated, it will be allocated to the dimensions
  !! of the read variable. Otherwise, its dimensions must match those of the
  !! read variable.
  !!
  subroutine read_3D_double( this, variable_name, container, requestor_name )
    use musica_constants,              only : musica_dk
    use musica_string,                 only : string_t
    import io_t
    class(io_t),                       intent(inout) :: this
    class(string_t),                   intent(in)    :: variable_name
    real(kind=musica_dk), allocatable, intent(inout) :: container(:,:,:)
    character(len=*),                  intent(in)    :: requestor_name
  end subroutine read_3D_double

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Reads 4D double-precision floating-point data
  !!
  !! If \c container is unallocated, it will be allocated to the dimensions
  !! of the read variable. Otherwise, its dimensions must match those of the
  !! read variable.
  !!
  subroutine read_4D_double( this, variable_name, container, requestor_name )
    use musica_constants,              only : musica_dk
    use musica_string,                 only : string_t
    import io_t
    class(io_t),                       intent(inout) :: this
    class(string_t),                   intent(in)    :: variable_name
    real(kind=musica_dk), allocatable, intent(inout) :: container(:,:,:,:)
    character(len=*),                  intent(in)    :: requestor_name
  end subroutine read_4D_double

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Reads 0D integer data
  subroutine read_0D_int( this, variable_name, container, requestor_name )
    use musica_string,                 only : string_t
    import io_t
    class(io_t),      intent(inout) :: this
    class(string_t),  intent(in)    :: variable_name
    integer,          intent(out)   :: container
    character(len=*), intent(in)    :: requestor_name
  end subroutine read_0D_int

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Reads 1D integer data
  !!
  !! If \c container is unallocated, it will be allocated to the dimensions
  !! of the read variable. Otherwise, its dimensions must match those of the
  !! read variable.
  !!
  subroutine read_1D_int( this, variable_name, container, requestor_name )
    use musica_string,                 only : string_t
    import io_t
    class(io_t),                   intent(inout) :: this
    class(string_t),               intent(in)    :: variable_name
    integer,          allocatable, intent(inout) :: container(:)
    character(len=*),              intent(in)    :: requestor_name
  end subroutine read_1D_int

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Writes 0D double data
  subroutine write_0D_double( this, variable_name, variable_data,             &
      requestor_name )
    use musica_constants,              only : musica_dk
    use musica_string,                 only : string_t
    import io_t
    class(io_t),          intent(inout) :: this
    type(string_t),       intent(in)    :: variable_name
    real(kind=musica_dk), intent(in)    :: variable_data
    character(len=*),     intent(in)    :: requestor_name
  end subroutine write_0D_double

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Writes 1D double data
  !!
  !! If the provided dimensions do not exist, they will be created based on
  !! the shape of the given data. If they do exist, they must be compatible
  !! with the shape of the given data.
  subroutine write_1D_double( this, variable_name, dimensions, variable_data, &
      requestor_name )
    use musica_constants,              only : musica_dk
    use musica_string,                 only : string_t
    import io_t
    class(io_t),          intent(inout) :: this
    type(string_t),       intent(in)    :: variable_name
    type(string_t),       intent(in)    :: dimensions
    real(kind=musica_dk), intent(in)    :: variable_data(:)
    character(len=*),     intent(in)    :: requestor_name
  end subroutine write_1D_double

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Writes 2D double data
  !!
  !! If the provided dimensions do not exist, they will be created based on
  !! the shape of the given data. If they do exist, they must be compatible
  !! with the shape of the given data.
  subroutine write_2D_double( this, variable_name, dimensions, variable_data, &
      requestor_name )
    use musica_constants,              only : musica_dk
    use musica_string,                 only : string_t
    import io_t
    class(io_t),          intent(inout) :: this
    type(string_t),       intent(in)    :: variable_name
    type(string_t),       intent(in)    :: dimensions(2)
    real(kind=musica_dk), intent(in)    :: variable_data(:,:)
    character(len=*),     intent(in)    :: requestor_name
  end subroutine write_2D_double

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Writes 3D double data
  !!
  !! If the provided dimensions do not exist, they will be created based on
  !! the shape of the given data. If they do exist, they must be compatible
  !! with the shape of the given data.
  subroutine write_3D_double( this, variable_name, dimensions, variable_data, &
      requestor_name )
    use musica_constants,              only : musica_dk
    use musica_string,                 only : string_t
    import io_t
    class(io_t),          intent(inout) :: this
    type(string_t),       intent(in)    :: variable_name
    type(string_t),       intent(in)    :: dimensions(3)
    real(kind=musica_dk), intent(in)    :: variable_data(:,:,:)
    character(len=*),     intent(in)    :: requestor_name
  end subroutine write_3D_double

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Writes 4D double data
  !!
  !! If the provided dimensions do not exist, they will be created based on
  !! the shape of the given data. If they do exist, they must be compatible
  !! with the shape of the given data.
  subroutine write_4D_double( this, variable_name, dimensions, variable_data, &
      requestor_name )
    use musica_constants,              only : musica_dk
    use musica_string,                 only : string_t
    import io_t
    class(io_t),          intent(inout) :: this
    type(string_t),       intent(in)    :: variable_name
    type(string_t),       intent(in)    :: dimensions(4)
    real(kind=musica_dk), intent(in)    :: variable_data(:,:,:,:)
    character(len=*),     intent(in)    :: requestor_name
  end subroutine write_4D_double

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Writes 0D int data
  subroutine write_0D_int( this, variable_name, variable_data,             &
      requestor_name )
    use musica_string,                 only : string_t
    import io_t
    class(io_t),      intent(inout) :: this
    type(string_t),   intent(in)    :: variable_name
    integer,          intent(in)    :: variable_data
    character(len=*), intent(in)    :: requestor_name
  end subroutine write_0D_int

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Writes 1D int data
  !!
  !! If the provided dimensions do not exist, they will be created based on
  !! the shape of the given data. If they do exist, they must be compatible
  !! with the shape of the given data.
  subroutine write_1D_int( this, variable_name, dimensions, variable_data, &
      requestor_name )
    use musica_string,                 only : string_t
    import io_t
    class(io_t),      intent(inout) :: this
    type(string_t),   intent(in)    :: variable_name
    type(string_t),   intent(in)    :: dimensions
    integer,          intent(in)    :: variable_data(:)
    character(len=*), intent(in)    :: requestor_name
  end subroutine write_1D_int

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Writes 0D double data to append 1D double data
  !!
  !! If the provided dimensions do not exist, they will be created based on
  !! the shape of the given data. If they do exist, they must be compatible
  !! with the shape of the given data.
  subroutine append_0D_double( this, variable_name, variable_units,           &
      append_dimension, append_index, variable_data, requestor_name )
    use musica_constants,              only : musica_dk
    use musica_string,                 only : string_t
    import io_t
    class(io_t),          intent(inout) :: this
    type(string_t),       intent(in)    :: variable_name
    type(string_t),       intent(in)    :: variable_units
    type(string_t),       intent(in)    :: append_dimension
    integer,              intent(in)    :: append_index
    real(kind=musica_dk), intent(in)    :: variable_data
    character(len=*),     intent(in)    :: requestor_name
  end subroutine append_0D_double

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Writes 1D double data to append 2D data
  !!
  !! If the provided dimensions do not exist, they will be created based on
  !! the shape of the given data. If they do exist, they must be compatible
  !! with the shape of the given data.
  subroutine append_1D_double( this, variable_name, variable_units,           &
      append_dimension, append_index, dimensions, variable_data,              &
      requestor_name )
    use musica_constants,              only : musica_dk
    use musica_string,                 only : string_t
    import io_t
    class(io_t),          intent(inout) :: this
    type(string_t),       intent(in)    :: variable_name
    type(string_t),       intent(in)    :: variable_units
    type(string_t),       intent(in)    :: append_dimension
    integer,              intent(in)    :: append_index
    type(string_t),       intent(in)    :: dimensions
    real(kind=musica_dk), intent(in)    :: variable_data(:)
    character(len=*),     intent(in)    :: requestor_name
  end subroutine append_1D_double

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Writes 2D double data to append 3D double data
  !!
  !! If the provided dimensions do not exist, they will be created based on
  !! the shape of the given data. If they do exist, they must be compatible
  !! with the shape of the given data.
  subroutine append_2D_double( this, variable_name, variable_units,           &
      append_dimension, append_index, dimensions, variable_data,              &
      requestor_name )
    use musica_constants,              only : musica_dk
    use musica_string,                 only : string_t
    import io_t
    class(io_t),          intent(inout) :: this
    type(string_t),       intent(in)    :: variable_name
    type(string_t),       intent(in)    :: variable_units
    type(string_t),       intent(in)    :: append_dimension
    integer,              intent(in)    :: append_index
    type(string_t),       intent(in)    :: dimensions(2)
    real(kind=musica_dk), intent(in)    :: variable_data(:,:)
    character(len=*),     intent(in)    :: requestor_name
  end subroutine append_2D_double

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Writes 3D double data to append 4D double data
  !!
  !! If the provided dimensions do not exist, they will be created based on
  !! the shape of the given data. If they do exist, they must be compatible
  !! with the shape of the given data.
  subroutine append_3D_double( this, variable_name, variable_units,           &
      append_dimension, append_index, dimensions, variable_data,              &
      requestor_name )
    use musica_constants,              only : musica_dk
    use musica_string,                 only : string_t
    import io_t
    class(io_t),          intent(inout) :: this
    type(string_t),       intent(in)    :: variable_name
    type(string_t),       intent(in)    :: variable_units
    type(string_t),       intent(in)    :: append_dimension
    integer,              intent(in)    :: append_index
    type(string_t),       intent(in)    :: dimensions(3)
    real(kind=musica_dk), intent(in)    :: variable_data(:,:,:)
    character(len=*),     intent(in)    :: requestor_name
  end subroutine append_3D_double

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Writes 0D int data to append 1D int data
  !!
  !! If the provided dimensions do not exist, they will be created based on
  !! the shape of the given data. If they do exist, they must be compatible
  !! with the shape of the given data.
  subroutine append_0D_int( this, variable_name, variable_units,              &
      append_dimension, append_index, variable_data, requestor_name )
    use musica_string,                 only : string_t
    import io_t
    class(io_t),      intent(inout) :: this
    type(string_t),   intent(in)    :: variable_name
    type(string_t),   intent(in)    :: variable_units
    type(string_t),   intent(in)    :: append_dimension
    integer,          intent(in)    :: append_index
    integer,          intent(in)    :: variable_data
    character(len=*), intent(in)    :: requestor_name
  end subroutine append_0D_int

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns whether a variable exists in the file
  logical function exists_char( this, variable_name, requestor_name )         &
      result( exists )
    import io_t
    class(io_t), intent(in) :: this
    character(len=*), intent(in) :: variable_name
    character(len=*), intent(in) :: requestor_name
  end function exists_char

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns whether a variable exists in the file
  logical function exists_string( this, variable_name, requestor_name )       &
      result( exists )
    use musica_string,                 only : string_t
    import io_t
    class(io_t),      intent(in) :: this
    type(string_t),   intent(in) :: variable_name
    character(len=*), intent(in) :: requestor_name
  end function exists_string

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the dimension names for a given variable
  function variable_dimensions( this, variable_name, requestor_name )         &
      result( dimensions )
    use musica_string,                 only : string_t
    import io_t
    type(string_t), allocatable  :: dimensions(:)
    class(io_t),      intent(in) :: this
    class(string_t),  intent(in) :: variable_name
    character(len=*), intent(in) :: requestor_name
  end function variable_dimensions

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the units for a given variable
  function variable_units( this, variable_name, requestor_name )
    use musica_string,                 only : string_t
    import io_t
    type(string_t)               :: variable_units
    class(io_t),      intent(in) :: this
    class(string_t),  intent(in) :: variable_name
    character(len=*), intent(in) :: requestor_name
  end function variable_units

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Sets the units for a given variable
  subroutine set_variable_units( this, variable_name, units, requestor_name )
    use musica_string,                 only : string_t
    import io_t
    class(io_t),      intent(inout) :: this
    class(string_t),  intent(in)    :: variable_name
    class(string_t),  intent(in)    :: units
    character(len=*), intent(in)    :: requestor_name
  end subroutine set_variable_units

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end interface

end module musica_io
