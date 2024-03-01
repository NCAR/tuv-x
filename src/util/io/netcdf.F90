! Copyright (C) 2021 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> The musica_io_netcdf module

!> The io_netcdf_t type and related functions
module musica_io_netcdf

  use musica_io,                       only : io_t
  use musica_string,                   only : string_t

  implicit none
  private

  public :: io_netcdf_t

  integer, parameter :: kUnknownFileId = -9999

  !> NetCDF file reader
  type, extends(io_t) :: io_netcdf_t
    integer        :: file_id_ = kUnknownFileId
    type(string_t) :: file_name_
  contains
    !> @name Data read functions
    !! @{
    procedure :: read_0D_double
    procedure :: read_1D_double
    procedure :: read_2D_double
    procedure :: read_3D_double
    procedure :: read_4D_double
    procedure :: read_0D_int
    procedure :: read_1D_int
    !> @}
    !> @name Data write functions
    !! @{
    procedure :: write_0D_double
    procedure :: write_1D_double
    procedure :: write_2D_double
    procedure :: write_3D_double
    procedure :: write_4D_double
    procedure :: write_0D_int
    procedure :: write_1D_int
    !> @}
    !> @name Data append functions
    !! @{
    procedure :: append_0D_double
    procedure :: append_1D_double
    procedure :: append_2D_double
    procedure :: append_3D_double
    procedure :: append_0D_int
    !! @}
    !> @name Returns whether a variable exists in the file
    !! @{
    procedure :: exists_char
    procedure :: exists_string
    !> @}
    !> Returns the dimension names for a given variable
    procedure :: variable_dimensions
    !> Returns the units for a given variable
    procedure :: variable_units
    !> Sets the units for a given variable
    procedure :: set_variable_units
    procedure, private :: is_open
    procedure, private :: variable_id
    procedure, private :: dimension_sizes
    procedure, private :: check_add_dimension
    procedure, private :: check_add_variable
    final :: finalize
  end type io_netcdf_t

  interface io_netcdf_t
    procedure :: constructor
  end interface

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Constructor for NetCDF file readers
  function constructor( file_name, read_only ) result( new_io )

    use musica_string,                 only : string_t
    use netcdf,                        only : nf90_create, nf90_open,         &
                                              NF90_NETCDF4, NF90_WRITE,       &
                                              NF90_NOWRITE

    type(io_netcdf_t), pointer    :: new_io
    type(string_t),    intent(in) :: file_name
    logical, optional, intent(in) :: read_only

    logical :: file_exists

    allocate( new_io )
    new_io%file_name_ = file_name
    if( present( read_only ) ) then
      if( read_only ) then
        call check_status( 233000996,                                         &
            nf90_open( file_name%to_char( ), NF90_NOWRITE, new_io%file_id_ ), &
            "Error openning file '"//file_name%to_char( )//"'" )
        return
      end if
    end if
    inquire( file = file_name%to_char( ), exist = file_exists )
    if( file_exists ) then
      call check_status( 126279520,                                           &
          nf90_open( file_name%to_char( ), NF90_WRITE, new_io%file_id_ ),     &
          "Error openning file '"//file_name%to_char( )//"'" )
    else
      call check_status( 427923808,                                           &
          nf90_create( file_name%to_char( ), NF90_NETCDF4, new_io%file_id_ ), &
          "Error creating file '"//file_name%to_char( )//"'" )
    end if

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Reads 0D double-precision floating-pointer data
  subroutine read_0D_double( this, variable_name, container, requestor_name )

    use musica_assert,                 only : assert_msg
    use musica_constants,              only : musica_dk
    use musica_string,                 only : to_char
    use netcdf,                        only : nf90_get_var

    class(io_netcdf_t),   intent(inout) :: this
    class(string_t),      intent(in)    :: variable_name
    real(kind=musica_dk), intent(out)   :: container
    character(len=*),     intent(in)    :: requestor_name

    integer :: var_id
    integer, allocatable :: dim_sizes(:)
    type(string_t) :: id_str

    id_str = "variable '"//variable_name//"' in file '"//this%file_name_//"'"
    call assert_msg( 879207328, this%is_open( ),                              &
                     "Trying to read from unopen file: "//id_str )
    var_id    = this%variable_id( variable_name )
    dim_sizes = this%dimension_sizes( variable_name )
    call assert_msg( 712409197, size( dim_sizes ) .eq. 0,                     &
                     "Wrong number of dimensions for "//                      &
                     trim( id_str%to_char( ) )//": Expected 0 got "//         &
                     trim( to_char( size( dim_sizes ) ) ) )
    call check_status( 190408927,                                             &
                       nf90_get_var( this%file_id_, var_id, container ),      &
                       "Error getting value for "//trim( id_str%to_char( ) ) )

  end subroutine read_0D_double

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Reads 1D double-precision floating-pointer data
  subroutine read_1D_double( this, variable_name, container, requestor_name )

    use musica_assert,                 only : assert_msg
    use musica_constants,              only : musica_dk
    use musica_string,                 only : to_char
    use netcdf,                        only : nf90_get_var

    class(io_netcdf_t),                intent(inout) :: this
    class(string_t),                   intent(in)    :: variable_name
    real(kind=musica_dk), allocatable, intent(inout) :: container(:)
    character(len=*),                  intent(in)    :: requestor_name

    integer :: var_id
    integer, allocatable :: dim_sizes(:)
    type(string_t) :: id_str

    id_str = "variable '"//variable_name//"' in file '"//this%file_name_//"'"
    call assert_msg( 163123652, this%is_open( ),                              &
                     "Trying to read from unopen file: "//id_str )
    var_id    = this%variable_id( variable_name )
    dim_sizes = this%dimension_sizes( variable_name )
    call assert_msg( 275441997, size( dim_sizes ) .eq. 1,                     &
                     "Wrong number of dimensions for "//                      &
                     trim( id_str%to_char( ) )//": Expected 1 got "//         &
                     trim( to_char( size( dim_sizes ) ) ) )
    if( allocated( container ) ) then
      call assert_msg( 976961669, size( container ) .eq. dim_sizes(1),        &
                       "Wrong size container for "//trim( id_str%to_char( ) ) &
                       //": Expected "//trim( to_char( dim_sizes(1) ) )//     &
                       " got "//trim( to_char( size( container ) ) ) )
    else
      allocate( container( dim_sizes(1) ) )
    end if
    call check_status( 722809843,                                             &
                       nf90_get_var( this%file_id_, var_id, container ),      &
                       "Error getting value for "//trim( id_str%to_char( ) ) )

  end subroutine read_1D_double

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Reads 2D double-precision floating-pointer data
  subroutine read_2D_double( this, variable_name, container, requestor_name )

    use musica_assert,                 only : assert_msg
    use musica_constants,              only : musica_dk
    use musica_string,                 only : to_char
    use netcdf,                        only : nf90_get_var

    class(io_netcdf_t),                intent(inout) :: this
    class(string_t),                   intent(in)    :: variable_name
    real(kind=musica_dk), allocatable, intent(inout) :: container(:,:)
    character(len=*),                  intent(in)    :: requestor_name

    integer :: var_id, i_dim
    integer, allocatable :: dim_sizes(:)
    type(string_t) :: id_str

    id_str = "variable '"//variable_name//"' in file '"//this%file_name_//"'"
    call assert_msg( 675787021, this%is_open( ),                              &
                     "Trying to read from unopen file: "//id_str )
    var_id    = this%variable_id( variable_name )
    dim_sizes = this%dimension_sizes( variable_name )
    call assert_msg( 400481613, size( dim_sizes ) .eq. 2,                     &
                     "Wrong number of dimensions for "//                      &
                     trim( id_str%to_char( ) )//": Expected 2 got "//         &
                     trim( to_char( size( dim_sizes ) ) ) )
    if( allocated( container ) ) then
      do i_dim = 1, 2
        call assert_msg( 230324709, size( container, i_dim ) .eq.             &
                                    dim_sizes( i_dim ),                       &
                         "Wrong size container for "//                        &
                         trim( id_str%to_char( ) )//": Expected "//           &
                         trim( to_char( dim_sizes( i_dim ) ) )//     &
                         " got "//trim( to_char( size( container, i_dim ) ) ) &
                         //" for dimension "//trim( to_char( i_dim ) ) )
      end do
    else
      allocate( container( dim_sizes(1), dim_sizes(2) ) )
    end if
    call check_status( 960167804,                                             &
                       nf90_get_var( this%file_id_, var_id, container ),      &
                       "Error getting value for "//trim( id_str%to_char( ) ) )

  end subroutine read_2D_double

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Reads 3D double-precision floating-pointer data
  subroutine read_3D_double( this, variable_name, container, requestor_name )

    use musica_assert,                 only : assert_msg
    use musica_constants,              only : musica_dk
    use musica_string,                 only : to_char
    use netcdf,                        only : nf90_get_var

    class(io_netcdf_t),                intent(inout) :: this
    class(string_t),                   intent(in)    :: variable_name
    real(kind=musica_dk), allocatable, intent(inout) :: container(:,:,:)
    character(len=*),                  intent(in)    :: requestor_name

    integer :: var_id, i_dim
    integer, allocatable :: dim_sizes(:)
    type(string_t) :: id_str

    id_str = "variable '"//variable_name//"' in file '"//this%file_name_//"'"
    call assert_msg( 539957265, this%is_open( ),                              &
                     "Trying to read from unopen file: "//id_str )
    var_id    = this%variable_id( variable_name )
    dim_sizes = this%dimension_sizes( variable_name )
    call assert_msg( 603060131, size( dim_sizes ) .eq. 3,                     &
                     "Wrong number of dimensions for "//                      &
                     trim( id_str%to_char( ) )//": Expected 3 got "//         &
                     trim( to_char( size( dim_sizes ) ) ) )
    if( allocated( container ) ) then
      do i_dim = 1, 3
        call assert_msg( 715378476, size( container, i_dim ) .eq.             &
                                    dim_sizes( i_dim ),                       &
                         "Wrong size container for "//                        &
                         trim( id_str%to_char( ) )//": Expected "//           &
                         trim( to_char( dim_sizes( i_dim ) ) )//     &
                         " got "//trim( to_char( size( container, i_dim ) ) ) &
                         //" for dimension "//trim( to_char( i_dim ) ) )
      end do
    else
      allocate( container( dim_sizes(1), dim_sizes(2), dim_sizes(3) ) )
    end if
    call check_status( 210172071,                                             &
                       nf90_get_var( this%file_id_, var_id, container ),      &
                       "Error getting value for "//trim( id_str%to_char( ) ) )

  end subroutine read_3D_double

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Reads 4D double-precision floating-pointer data
  subroutine read_4D_double( this, variable_name, container, requestor_name )

    use musica_assert,                 only : assert_msg
    use musica_constants,              only : musica_dk
    use musica_string,                 only : to_char
    use netcdf,                        only : nf90_get_var

    class(io_netcdf_t),                intent(inout) :: this
    class(string_t),                   intent(in)    :: variable_name
    real(kind=musica_dk), allocatable, intent(inout) :: container(:,:,:,:)
    character(len=*),                  intent(in)    :: requestor_name

    integer :: var_id, i_dim
    integer, allocatable :: dim_sizes(:)
    type(string_t) :: id_str

    id_str = "variable '"//variable_name//"' in file '"//this%file_name_//"'"
    call assert_msg( 198190218, this%is_open( ),                              &
                     "Trying to read from unopen file: "//id_str )
    var_id    = this%variable_id( variable_name )
    dim_sizes = this%dimension_sizes( variable_name )
    call assert_msg( 650822371, size( dim_sizes ) .eq. 4,                     &
                     "Wrong number of dimensions for "//                      &
                     trim( id_str%to_char( ) )//": Expected 4 got "//         &
                     trim( to_char( size( dim_sizes ) ) ) )
    if( allocated( container ) ) then
      do i_dim = 1, 4
        call assert_msg( 820979275, size( container, i_dim ) .eq.             &
                                    dim_sizes( i_dim ),                       &
                         "Wrong size container for "//                        &
                         trim( id_str%to_char( ) )//": Expected "//           &
                         trim( to_char( dim_sizes( i_dim ) ) )//     &
                         " got "//trim( to_char( size( container, i_dim ) ) ) &
                         //" for dimension "//trim( to_char( i_dim ) ) )
      end do
    else
      allocate( container( dim_sizes(1), dim_sizes(2), dim_sizes(3),          &
                           dim_sizes(4) ) )
    end if
    call check_status( 708660930,                                             &
                       nf90_get_var( this%file_id_, var_id, container ),      &
                       "Error getting value for "//trim( id_str%to_char( ) ) )

  end subroutine read_4D_double

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Reads 0D integer data
  subroutine read_0D_int( this, variable_name, container, requestor_name )

    use musica_assert,                 only : assert_msg
    use musica_string,                 only : to_char
    use netcdf,                        only : nf90_get_var

    class(io_netcdf_t), intent(inout) :: this
    class(string_t),    intent(in)    :: variable_name
    integer,            intent(out)   :: container
    character(len=*),   intent(in)    :: requestor_name

    integer :: var_id
    integer, allocatable :: dim_sizes(:)
    type(string_t) :: id_str

    id_str = "variable '"//variable_name//"' in file '"//this%file_name_//"'"
    call assert_msg( 418014896, this%is_open( ),                              &
                     "Trying to read from unopen file: "//id_str )
    var_id    = this%variable_id( variable_name )
    dim_sizes = this%dimension_sizes( variable_name )
    call assert_msg( 747800090, size( dim_sizes ) .eq. 0,                     &
                     "Wrong number of dimensions for "//                      &
                     trim( id_str%to_char( ) )//": Expected 0 got "//         &
                     trim( to_char( size( dim_sizes ) ) ) )
    call check_status( 860118435,                                             &
                       nf90_get_var( this%file_id_, var_id, container ),      &
                       "Error getting value for "//trim( id_str%to_char( ) ) )

  end subroutine read_0D_int

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Reads 1D integer data
  subroutine read_1D_int( this, variable_name, container, requestor_name )

    use musica_assert,                 only : assert_msg
    use musica_string,                 only : to_char
    use netcdf,                        only : nf90_get_var

    class(io_netcdf_t),   intent(inout) :: this
    class(string_t),      intent(in)    :: variable_name
    integer, allocatable, intent(inout) :: container(:)
    character(len=*),     intent(in)    :: requestor_name

    integer :: var_id
    integer, allocatable :: dim_sizes(:)
    type(string_t) :: id_str

    id_str = "variable '"//variable_name//"' in file '"//this%file_name_//"'"
    call assert_msg( 121652260, this%is_open( ),                              &
                     "Trying to read from unopen file: "//id_str )
    var_id    = this%variable_id( variable_name )
    dim_sizes = this%dimension_sizes( variable_name )
    call assert_msg( 798921103, size( dim_sizes ) .eq. 1,                     &
                     "Wrong number of dimensions for "//                      &
                     trim( id_str%to_char( ) )//": Expected 1 got "//         &
                     trim( to_char( size( dim_sizes ) ) ) )
    if( allocated( container ) ) then
      call assert_msg( 346288950, size( container ) .eq. dim_sizes(1),        &
                       "Wrong size container for "//trim( id_str%to_char( ) ) &
                       //": Expected "//trim( to_char( dim_sizes(1) ) )//     &
                       " got "//trim( to_char( size( container ) ) ) )
    else
      allocate( container( dim_sizes(1) ) )
    end if
    call check_status( 458607295,                                             &
                       nf90_get_var( this%file_id_, var_id, container ),      &
                       "Error getting value for "//trim( id_str%to_char( ) ) )

  end subroutine read_1D_int

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Writes 0D double data
  subroutine write_0D_double( this, variable_name, variable_data,             &
      requestor_name )

    use musica_assert,                 only : assert_msg
    use musica_constants,              only : musica_dk
    use musica_string,                 only : string_t
    use netcdf,                        only : nf90_def_var, nf90_put_var,     &
                                              NF90_DOUBLE

    class(io_netcdf_t),   intent(inout) :: this
    type(string_t),       intent(in)    :: variable_name
    real(kind=musica_dk), intent(in)    :: variable_data
    character(len=*),     intent(in)    :: requestor_name

    integer :: var_id
    type(string_t) :: id_str

    id_str = "variable '"//variable_name//"' in file '"//this%file_name_//"'"
    call assert_msg( 576950310, this%is_open( ),                              &
                     "Trying to write to an unopen file: "//id_str )
    call check_status( 550080126, nf90_def_var( this%file_id_,                &
                                                variable_name%to_char( ),     &
                                                NF90_DOUBLE,                  &
                                                varid = var_id  ),            &
                       "Error creating "//trim( id_str%to_char( ) ) )
    call check_status( 540003807,                                             &
                       nf90_put_var( this%file_id_, var_id, variable_data ),  &
                       "Error writing to "//trim( id_str%to_char( ) ) )

  end subroutine write_0D_double

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Writes 1D double data
  !!
  !! If the provided dimensions do not exist, they will be created based on
  !! the shape of the given data. If they do exist, they must be compatible
  !! with the shape of the given data.
  subroutine write_1D_double( this, variable_name, dimensions, variable_data, &
      requestor_name )

    use musica_assert,                 only : assert_msg
    use musica_constants,              only : musica_dk
    use musica_string,                 only : string_t
    use netcdf,                        only : nf90_def_var, nf90_put_var,     &
                                              NF90_DOUBLE

    class(io_netcdf_t),   intent(inout) :: this
    type(string_t),       intent(in)    :: variable_name
    type(string_t),       intent(in)    :: dimensions
    real(kind=musica_dk), intent(in)    :: variable_data(:)
    character(len=*),     intent(in)    :: requestor_name

    integer :: var_id, dimids(1)
    type(string_t) :: id_str

    id_str = "variable '"//variable_name//"' in file '"//this%file_name_//"'"
    call assert_msg( 616828888, this%is_open( ),                              &
                     "Trying to write to an unopen file: "//id_str )
    dimids(1) = this%check_add_dimension( dimensions%to_char( ),              &
                                          size( variable_data ) )
    call check_status( 111622483, nf90_def_var( this%file_id_,                &
                                                variable_name%to_char( ),     &
                                                NF90_DOUBLE, dimids = dimids, &
                                                varid = var_id  ),            &
                       "Error creating "//trim( id_str%to_char( ) ) )
    call check_status( 841465578,                                             &
                       nf90_put_var( this%file_id_, var_id, variable_data ),  &
                       "Error writing to "//trim( id_str%to_char( ) ) )

  end subroutine write_1D_double

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Writes 2D double data
  !!
  !! If the provided dimensions do not exist, they will be created based on
  !! the shape of the given data. If they do exist, they must be compatible
  !! with the shape of the given data.
  subroutine write_2D_double( this, variable_name, dimensions, variable_data, &
      requestor_name )

    use musica_assert,                 only : assert_msg
    use musica_constants,              only : musica_dk
    use musica_string,                 only : string_t
    use netcdf,                        only : nf90_def_var, nf90_put_var,     &
                                              NF90_DOUBLE

    class(io_netcdf_t),   intent(inout) :: this
    type(string_t),       intent(in)    :: variable_name
    type(string_t),       intent(in)    :: dimensions(2)
    real(kind=musica_dk), intent(in)    :: variable_data(:,:)
    character(len=*),     intent(in)    :: requestor_name

    integer :: var_id, dimids(2)
    type(string_t) :: id_str

    id_str = "variable '"//variable_name//"' in file '"//this%file_name_//"'"
    call assert_msg( 186994325, this%is_open( ),                              &
                     "Trying to write to an unopen file: "//id_str )
    dimids(1) = this%check_add_dimension( dimensions(1)%to_char( ),           &
                                          size( variable_data, 1 ) )
    dimids(2) = this%check_add_dimension( dimensions(2)%to_char( ),           &
                                          size( variable_data, 2 ) )
    call check_status( 916837420, nf90_def_var( this%file_id_,                &
                                                variable_name%to_char( ),     &
                                                NF90_DOUBLE, dimids = dimids, &
                                                varid = var_id  ),            &
                       "Error creating "//trim( id_str%to_char( ) ) )
    call check_status( 464205267,                                             &
                       nf90_put_var( this%file_id_, var_id, variable_data ),  &
                       "Error writing to "//trim( id_str%to_char( ) ) )

  end subroutine write_2D_double

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Writes 3D double data
  !!
  !! If the provided dimensions do not exist, they will be created based on
  !! the shape of the given data. If they do exist, they must be compatible
  !! with the shape of the given data.
  subroutine write_3D_double( this, variable_name, dimensions, variable_data, &
      requestor_name )

    use musica_assert,                 only : assert_msg
    use musica_constants,              only : musica_dk
    use musica_string,                 only : string_t
    use netcdf,                        only : nf90_def_var, nf90_put_var,     &
                                              NF90_DOUBLE

    class(io_netcdf_t),   intent(inout) :: this
    type(string_t),       intent(in)    :: variable_name
    type(string_t),       intent(in)    :: dimensions(3)
    real(kind=musica_dk), intent(in)    :: variable_data(:,:,:)
    character(len=*),     intent(in)    :: requestor_name

    integer :: var_id, dimids(3)
    type(string_t) :: id_str

    id_str = "variable '"//variable_name//"' in file '"//this%file_name_//"'"
    call assert_msg( 232851031, this%is_open( ),                              &
                     "Trying to write to an unopen file: "//id_str )
    dimids(1) = this%check_add_dimension( dimensions(1)%to_char( ),           &
                                          size( variable_data, 1 ) )
    dimids(2) = this%check_add_dimension( dimensions(2)%to_char( ),           &
                                          size( variable_data, 2 ) )
    dimids(3) = this%check_add_dimension( dimensions(3)%to_char( ),           &
                                          size( variable_data, 3 ) )
    call check_status( 403007935, nf90_def_var( this%file_id_,                &
                                                variable_name%to_char( ),     &
                                                NF90_DOUBLE, dimids = dimids, &
                                                varid = var_id  ),            &
                       "Error creating "//trim( id_str%to_char( ) ) )
    call check_status( 573164839,                                             &
                       nf90_put_var( this%file_id_, var_id, variable_data ),  &
                       "Error writing to "//trim( id_str%to_char( ) ) )

  end subroutine write_3D_double

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Writes 4D double data
  !!
  !! If the provided dimensions do not exist, they will be created based on
  !! the shape of the given data. If they do exist, they must be compatible
  !! with the shape of the given data.
  subroutine write_4D_double( this, variable_name, dimensions, variable_data, &
      requestor_name )

    use musica_assert,                 only : assert_msg
    use musica_constants,              only : musica_dk
    use musica_string,                 only : string_t
    use netcdf,                        only : nf90_def_var, nf90_put_var,     &
                                              NF90_DOUBLE

    class(io_netcdf_t),   intent(inout) :: this
    type(string_t),       intent(in)    :: variable_name
    type(string_t),       intent(in)    :: dimensions(4)
    real(kind=musica_dk), intent(in)    :: variable_data(:,:,:,:)
    character(len=*),     intent(in)    :: requestor_name

    integer :: var_id, dimids(4)
    type(string_t) :: id_str

    id_str = "variable '"//variable_name//"' in file '"//this%file_name_//"'"
    call assert_msg( 338451830, this%is_open( ),                              &
                     "Trying to write to an unopen file: "//id_str )
    dimids(1) = this%check_add_dimension( dimensions(1)%to_char( ),           &
                                          size( variable_data, 1 ) )
    dimids(2) = this%check_add_dimension( dimensions(2)%to_char( ),           &
                                          size( variable_data, 2 ) )
    dimids(3) = this%check_add_dimension( dimensions(3)%to_char( ),           &
                                          size( variable_data, 3 ) )
    dimids(4) = this%check_add_dimension( dimensions(4)%to_char( ),           &
                                          size( variable_data, 4 ) )
    call check_status( 233303326, nf90_def_var( this%file_id_,                &
                                                variable_name%to_char( ),     &
                                                NF90_DOUBLE, dimids = dimids, &
                                                varid = var_id  ),            &
                       "Error creating "//trim( id_str%to_char( ) ) )
    call check_status( 680671172,                                             &
                       nf90_put_var( this%file_id_, var_id, variable_data ),  &
                       "Error writing to "//trim( id_str%to_char( ) ) )

  end subroutine write_4D_double

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Writes 0D int data
  subroutine write_0D_int( this, variable_name, variable_data,                &
      requestor_name )

    use musica_assert,                 only : assert_msg
    use musica_string,                 only : string_t
    use netcdf,                        only : nf90_def_var, nf90_put_var,     &
                                              NF90_INT

    class(io_netcdf_t), intent(inout) :: this
    type(string_t),     intent(in)    :: variable_name
    integer,            intent(in)    :: variable_data
    character(len=*),   intent(in)    :: requestor_name

    integer :: var_id
    type(string_t) :: id_str

    id_str = "variable '"//variable_name//"' in file '"//this%file_name_//"'"
    call assert_msg( 834034211, this%is_open( ),                              &
                     "Trying to write to an unopen file: "//id_str )
    call check_status( 998926808, nf90_def_var( this%file_id_,                &
                                                variable_name%to_char( ),     &
                                                NF90_INT,                     &
                                                varid = var_id  ),            &
                       "Error creating "//trim( id_str%to_char( ) ) )
    call check_status( 546294655,                                             &
                       nf90_put_var( this%file_id_, var_id, variable_data ),  &
                       "Error writing to "//trim( id_str%to_char( ) ) )

  end subroutine write_0D_int

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Writes 1D int data
  !!
  !! If the provided dimensions do not exist, they will be created based on
  !! the shape of the given data. If they do exist, they must be compatible
  !! with the shape of the given data.
  subroutine write_1D_int( this, variable_name, dimensions, variable_data,    &
      requestor_name )

    use musica_assert,                 only : assert_msg
    use musica_string,                 only : string_t
    use netcdf,                        only : nf90_def_var, nf90_put_var,     &
                                              NF90_INT

    class(io_netcdf_t), intent(inout) :: this
    type(string_t),     intent(in)    :: variable_name
    type(string_t),     intent(in)    :: dimensions
    integer,            intent(in)    :: variable_data(:)
    character(len=*),   intent(in)    :: requestor_name

    integer :: var_id, dimids(1)
    type(string_t) :: id_str

    id_str = "variable '"//variable_name//"' in file '"//this%file_name_//"'"
    call assert_msg( 769478106, this%is_open( ),                              &
                     "Trying to write to an unopen file: "//id_str )
    dimids(1) = this%check_add_dimension( dimensions%to_char( ),              &
                                          size( variable_data ) )
    call check_status( 257101860, nf90_def_var( this%file_id_,                &
                                                variable_name%to_char( ),     &
                                                NF90_INT, dimids = dimids,    &
                                                varid = var_id  ),            &
                       "Error creating "//trim( id_str%to_char( ) ) )
    call check_status( 427258764,                                             &
                       nf90_put_var( this%file_id_, var_id, variable_data ),  &
                       "Error writing to "//trim( id_str%to_char( ) ) )

  end subroutine write_1D_int

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Writes 0D double data to append 1D double data
  !!
  !! If the provided dimensions do not exist, they will be created based on
  !! the shape of the given data. If they do exist, they must be compatible
  !! with the shape of the given data.
  subroutine append_0D_double( this, variable_name, variable_units,           &
      append_dimension, append_index, variable_data, requestor_name )

    use musica_assert,                 only : assert_msg
    use musica_constants,              only : musica_dk
    use musica_string,                 only : string_t
    use netcdf,                        only : nf90_def_var, nf90_put_var,     &
                                              NF90_DOUBLE, NF90_UNLIMITED

    class(io_netcdf_t),   intent(inout) :: this
    type(string_t),       intent(in)    :: variable_name
    type(string_t),       intent(in)    :: variable_units
    type(string_t),       intent(in)    :: append_dimension
    integer,              intent(in)    :: append_index
    real(kind=musica_dk), intent(in)    :: variable_data
    character(len=*),     intent(in)    :: requestor_name

    integer :: varid, dimids(1), start_ids(1), dim_sizes(0)
    type(string_t) :: id_str, dimensions(0)

    id_str = "variable '"//variable_name//"' in file '"//this%file_name_//"'"
    call assert_msg( 660803774, this%is_open( ),                              &
                     "Trying to write to an unopen file: "//id_str )
    call this%check_add_variable( variable_name, variable_units, NF90_DOUBLE, &
                                  append_dimension, dimensions,               &
                                  dim_sizes, varid, dimids, start_ids )
    start_ids(1) = append_index
    call check_status( 320489966,                                             &
                       nf90_put_var( this%file_id_, varid,                    &
                                     (/ variable_data /), start = start_ids,  &
                                     count = (/ 1 /) ),                       &
                       "Error writing to "//trim( id_str%to_char( ) ) )

  end subroutine append_0D_double

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Writes 1D double data to append 2D double data
  !!
  !! If the provided dimensions do not exist, they will be created based on
  !! the shape of the given data. If they do exist, they must be compatible
  !! with the shape of the given data.
  subroutine append_1D_double( this, variable_name, variable_units,           &
      append_dimension, append_index, dimensions, variable_data,              &
      requestor_name )

    use musica_assert,                 only : assert_msg
    use musica_constants,              only : musica_dk
    use musica_string,                 only : string_t
    use netcdf,                        only : nf90_def_var, nf90_put_var,     &
                                              NF90_DOUBLE, NF90_UNLIMITED

    class(io_netcdf_t),   intent(inout) :: this
    type(string_t),       intent(in)    :: variable_name
    type(string_t),       intent(in)    :: variable_units
    type(string_t),       intent(in)    :: append_dimension
    integer,              intent(in)    :: append_index
    type(string_t),       intent(in)    :: dimensions
    real(kind=musica_dk), intent(in)    :: variable_data(:)
    character(len=*),     intent(in)    :: requestor_name

    integer :: varid, dim_sizes(1), dimids(2), start_ids(2)
    type(string_t) :: id_str

    id_str = "variable '"//variable_name//"' in file '"//this%file_name_//"'"
    call assert_msg( 246721328, this%is_open( ),                              &
                     "Trying to write to an unopen file: "//id_str )
    dim_sizes(1) = size( variable_data )
    call this%check_add_variable( variable_name, variable_units, NF90_DOUBLE, &
                                  append_dimension, (/ dimensions /),         &
                                  dim_sizes, varid, dimids, start_ids )
    start_ids(1) = append_index
    call check_status( 641514922,                                             &
                       nf90_put_var( this%file_id_, varid,                    &
                                     (/ variable_data /), start = start_ids,  &
                                     count = (/ 1, size( variable_data ) /) ),&
                       "Error writing to "//trim( id_str%to_char( ) ) )

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

    use musica_assert,                 only : assert_msg
    use musica_constants,              only : musica_dk
    use musica_string,                 only : string_t
    use netcdf,                        only : nf90_def_var, nf90_put_var,     &
                                              NF90_DOUBLE, NF90_UNLIMITED

    class(io_netcdf_t),   intent(inout) :: this
    type(string_t),       intent(in)    :: variable_name
    type(string_t),       intent(in)    :: variable_units
    type(string_t),       intent(in)    :: append_dimension
    integer,              intent(in)    :: append_index
    type(string_t),       intent(in)    :: dimensions(2)
    real(kind=musica_dk), intent(in)    :: variable_data(:,:)
    character(len=*),     intent(in)    :: requestor_name

    integer :: varid, dim_sizes(2), dimids(3), start_ids(3)
    type(string_t) :: id_str

    id_str = "variable '"//variable_name//"' in file '"//this%file_name_//"'"
    call assert_msg( 264592928, this%is_open( ),                              &
                     "Trying to write to an unopen file: "//id_str )
    dim_sizes(1) = size( variable_data, 1 )
    dim_sizes(2) = size( variable_data, 2 )
    call this%check_add_variable( variable_name, variable_units, NF90_DOUBLE, &
                                  append_dimension, dimensions,               &
                                  dim_sizes, varid, dimids, start_ids )
    start_ids(1) = append_index
    call check_status( 889287519,                                             &
                       nf90_put_var( this%file_id_, varid,                    &
                                     variable_data, start = start_ids,        &
                                     count = (/ 1, size( variable_data, 1 ),  &
                                                size( variable_data, 2 ) /) ),&
                       "Error writing to "//trim( id_str%to_char( ) ) )

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

    use musica_assert,                 only : assert_msg
    use musica_constants,              only : musica_dk
    use musica_string,                 only : string_t
    use netcdf,                        only : nf90_def_var, nf90_put_var,     &
                                              NF90_DOUBLE, NF90_UNLIMITED

    class(io_netcdf_t),   intent(inout) :: this
    type(string_t),       intent(in)    :: variable_name
    type(string_t),       intent(in)    :: variable_units
    type(string_t),       intent(in)    :: append_dimension
    integer,              intent(in)    :: append_index
    type(string_t),       intent(in)    :: dimensions(3)
    real(kind=musica_dk), intent(in)    :: variable_data(:,:,:)
    character(len=*),     intent(in)    :: requestor_name

    integer :: varid, dim_sizes(3), dimids(4), start_ids(4)
    type(string_t) :: id_str

    id_str = "variable '"//variable_name//"' in file '"//this%file_name_//"'"
    call assert_msg( 351946623, this%is_open( ),                              &
                     "Trying to write to an unopen file: "//id_str )
    dim_sizes(1) = size( variable_data, 1 )
    dim_sizes(2) = size( variable_data, 2 )
    dim_sizes(3) = size( variable_data, 3 )
    call this%check_add_variable( variable_name, variable_units, NF90_DOUBLE, &
                                  append_dimension, dimensions,               &
                                  dim_sizes, varid, dimids, start_ids )
    start_ids(1) = append_index
    call check_status( 181789719,                                             &
                       nf90_put_var( this%file_id_, varid,                    &
                                     variable_data, start = start_ids,        &
                                     count = (/ 1, size( variable_data, 1 ),  &
                                                size( variable_data, 2 ),     &
                                                size( variable_data, 3 ) /) ),&
                       "Error writing to "//trim( id_str%to_char( ) ) )

  end subroutine append_3D_double

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Writes 0D int data to append 1D int data
  !!
  !! If the provided dimensions do not exist, they will be created based on
  !! the shape of the given data. If they do exist, they must be compatible
  !! with the shape of the given data.
  subroutine append_0D_int( this, variable_name, variable_units,              &
      append_dimension, append_index, variable_data, requestor_name )

    use musica_assert,                 only : assert_msg
    use musica_string,                 only : string_t
    use netcdf,                        only : nf90_def_var, nf90_put_var,     &
                                              NF90_INT, NF90_UNLIMITED

    class(io_netcdf_t), intent(inout) :: this
    type(string_t),     intent(in)    :: variable_name
    type(string_t),     intent(in)    :: variable_units
    type(string_t),     intent(in)    :: append_dimension
    integer,            intent(in)    :: append_index
    integer,            intent(in)    :: variable_data
    character(len=*),   intent(in)    :: requestor_name

    integer :: varid, dimids(1), start_ids(1), dim_sizes(0)
    type(string_t) :: id_str, dimensions(0)

    id_str = "variable '"//variable_name//"' in file '"//this%file_name_//"'"
    call assert_msg( 896317785, this%is_open( ),                              &
                     "Trying to write to an unopen file: "//id_str )
    call this%check_add_variable( variable_name, variable_units, NF90_INT,    &
                                  append_dimension, dimensions,               &
                                  dim_sizes, varid, dimids, start_ids )
    start_ids(1) = append_index
    call check_status( 108636131,                                             &
                       nf90_put_var( this%file_id_, varid,                    &
                                     (/ variable_data /), start = start_ids,  &
                                     count = (/ 1 /) ),                       &
                       "Error writing to "//trim( id_str%to_char( ) ) )

  end subroutine append_0D_int

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns whether a variable exists in the file
  logical function exists_char( this, variable_name, requestor_name )         &
      result( exists )

    use netcdf,                        only : nf90_inq_varid,                 &
                                              NF90_ENOTVAR

    class(io_netcdf_t), intent(in) :: this
    character(len=*),   intent(in) :: variable_name
    character(len=*),   intent(in) :: requestor_name

    integer :: var_id, err_id

    err_id = nf90_inq_varid( this%file_id_, variable_name, var_id )

    exists = .false.
    if( err_id == NF90_ENOTVAR ) return
    call check_status( 855364555, err_id, "Error trying to find variable '"// &
                       variable_name//"' in NetCDF file '"//                  &
                       trim( this%file_name_%to_char( ) )//"' for "//         &
                       requestor_name )
    exists = .true.

  end function exists_char

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns whether a variable exists in the file
  logical function exists_string( this, variable_name, requestor_name )       &
      result( exists )

    use musica_string,                 only : string_t

    class(io_netcdf_t), intent(in) :: this
    type(string_t),     intent(in) :: variable_name
    character(len=*),   intent(in) :: requestor_name

    exists = this%exists_char( variable_name%to_char( ), requestor_name )

  end function exists_string

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the dimension names for a given variable
  function variable_dimensions( this, variable_name, requestor_name )         &
      result( dimensions )

    use musica_string,                 only : to_char
    use netcdf,                        only : NF90_MAX_NAME,                  &
                                              nf90_inquire_variable,          &
                                              nf90_inquire_dimension

    type(string_t),     allocatable :: dimensions(:)
    class(io_netcdf_t), intent(in)  :: this
    class(string_t),    intent(in)  :: variable_name
    character(len=*),   intent(in)  :: requestor_name

    integer :: var_id, i_dim, n_dims
    integer, allocatable :: dimids(:)
    type(string_t) :: id_str
    character(len=NF90_MAX_NAME) :: dim_name

    id_str = "variable '"//variable_name//"' in file '"//this%file_name_//"'"
    var_id = this%variable_id( variable_name )
    call check_status( 744311319,                                             &
        nf90_inquire_variable( this%file_id_, var_id, ndims = n_dims ),       &
        "Error getting number of dimensions for "//id_str%to_char( ) )
    allocate( dimids( n_dims ) )
    call check_status( 104014576,                                             &
        nf90_inquire_variable( this%file_id_, var_id, dimids = dimids ),      &
        "Error getting dimesions for "//id_str%to_char( ) )
    allocate( dimensions( n_dims ) )
    do i_dim = 1, n_dims
      call check_status( 788714786,                                           &
          nf90_inquire_dimension( this%file_id_, dimids( i_dim ),             &
                                  name = dim_name ),&
          "Error getting dimesion size "//trim( to_char( i_dim ) )//" for "// &
          id_str%to_char( ) )
      dimensions( i_dim ) = trim( dim_name )
    end do

  end function variable_dimensions

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the units for a given variable
  function variable_units( this, variable_name, requestor_name )

    use netcdf,                        only : NF90_MAX_NAME,                  &
                                              nf90_get_att

    type(string_t)                 :: variable_units
    class(io_netcdf_t), intent(in) :: this
    class(string_t),    intent(in) :: variable_name
    character(len=*),   intent(in) :: requestor_name

    integer :: var_id
    type(string_t) :: id_str
    character(len=NF90_MAX_NAME) :: units

    id_str = "variable '"//variable_name//"' in file '"//this%file_name_//"'"
    var_id = this%variable_id( variable_name )
    call check_status( 301987512,                                             &
                       nf90_get_att( this%file_id_, var_id, "units", units ), &
                       "Error getting units for "//trim( id_str%to_char( ) ) )
    variable_units = trim( units )

  end function variable_units

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Sets the units for a given variable
  subroutine set_variable_units( this, variable_name, units, requestor_name )

    use musica_string,                 only : string_t
    use netcdf,                        only : nf90_put_att

    class(io_netcdf_t), intent(inout) :: this
    class(string_t),    intent(in)    :: variable_name
    class(string_t),    intent(in)    :: units
    character(len=*),   intent(in)    :: requestor_name

    integer :: var_id
    type(string_t) :: id_str

    id_str = "variable '"//variable_name//"' in file '"//this%file_name_//"'"
    var_id = this%variable_id( variable_name )
    call check_status( 235495983,                                             &
                       nf90_put_att( this%file_id_, var_id, "units",          &
                                     units%to_char( ) ),                      &
                       "Error setting units for "//trim( id_str%to_char( ) ) )

  end subroutine set_variable_units

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns whether a file is open or not
  logical function is_open( this )

    class(io_netcdf_t), intent(in) :: this

    is_open = this%file_id_ .ne. kUnknownFileId

  end function is_open

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns a variable's id in the NetCDF file
  integer function variable_id( this, variable_name )

    use musica_assert,                 only : assert_msg
    use netcdf,                        only : nf90_inq_varid

    class(io_netcdf_t), intent(in) :: this
    class(string_t),    intent(in) :: variable_name

    call assert_msg( 249726322, this%is_open( ),                              &
                     "Trying to read from unopen file: '"//                   &
                     this%file_name_//"'" )
    call check_status( 153462424,                                             &
                       nf90_inq_varid( this%file_id_,                         &
                                       variable_name%to_char( ),              &
                                       variable_id ), &
                       "Cannot find variable '"//                             &
                       trim( variable_name%to_char( ) )//"' in file '"//      &
                       trim( this%file_name_%to_char( ) )//"'" )

  end function variable_id

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the dimensions for variable in the NetCDF file
  function dimension_sizes( this, variable_name ) result( dim_sizes )

    use musica_assert,                 only : assert_msg
    use musica_string,                 only : to_char
    use netcdf,                        only : nf90_inquire_variable,          &
                                              nf90_inquire_dimension

    integer,            allocatable :: dim_sizes(:)
    class(io_netcdf_t), intent(in)  :: this
    class(string_t),    intent(in)  :: variable_name

    integer :: var_id, n_dims, i_dim
    integer, allocatable :: dimids(:)
    type(string_t) :: id_str

    id_str = "variable '"//variable_name//"' in file '"//this%file_name_//"'"
    call assert_msg( 191887763, this%is_open( ),                              &
                     "Trying to read from unopen file: "//id_str )
    var_id = this%variable_id( variable_name )
    call check_status( 516121527,                                             &
        nf90_inquire_variable( this%file_id_, var_id, ndims = n_dims ),       &
        "Error getting number of dimensions for "//trim( id_str%to_char( ) ) )
    allocate( dimids( n_dims ) )
    call check_status( 269878960,                                             &
        nf90_inquire_variable( this%file_id_, var_id, dimids = dimids ),      &
        "Error getting dimensions for "//trim( id_str%to_char( ) ) )
    allocate( dim_sizes( n_dims ) )
    do i_dim = 1, n_dims
      call check_status( 770273353,                                           &
          nf90_inquire_dimension( this%file_id_, dimids( i_dim ),             &
                                  len = dim_sizes( i_dim ) ),                 &
          "Error getting dimension size "//trim( to_char( i_dim ) )//" for "//&
          trim( id_str%to_char( ) ) )
    end do

  end function dimension_sizes

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Checks if a dimension exists and verifies its size
  !!
  !! If the dimension does not exist, it is created. The dimension id is
  !! returned.
  function check_add_dimension( this, dim_name, dim_size ) result( dimid )

    use musica_assert,                 only : assert_msg
    use musica_string,                 only : string_t, to_char
    use netcdf,                        only : nf90_inq_dimid,                 &
                                              nf90_inquire_dimension,         &
                                              nf90_def_dim,                   &
                                              NF90_NOERR, NF90_UNLIMITED

    integer                           :: dimid
    class(io_netcdf_t), intent(inout) :: this
    character(len=*),   intent(in)    :: dim_name
    integer,            intent(in)    :: dim_size

    integer :: ierr, curr_size
    type(string_t) :: id_str

    id_str = "dimension '"//dim_name//"' in file '"//this%file_name_//"'"

    ierr = nf90_inq_dimid( this%file_id_, dim_name, dimid )
    if( ierr == NF90_NOERR ) then
      ! dimension exists, check its size, unless it's unlimited
      if( dim_size .ne. NF90_UNLIMITED ) then
        call check_status( 737744716,                                         &
                           nf90_inquire_dimension( this%file_id_, dimid,      &
                                                   len = curr_size ),         &
                           "NetCDF file error for "//                         &
                           trim( id_str%to_char( ) ) )
        call assert_msg( 343403417, curr_size == dim_size,                    &
                         "Dimension mismatch for "//trim( id_str%to_char( ) ) &
                         //"; Expected "//trim( to_char( curr_size ) )//      &
                         ", got "//trim( to_char( dim_size ) ) )
      end if
    else
      call check_status( 947493075,                                           &
                         nf90_def_dim( this%file_id_, dim_name, dim_size,     &
                                       dimid ),                               &
                         "NetCDF file error for "//trim( id_str%to_char( ) ) )
    end if

  end function check_add_dimension

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Checks for an appendable variable in the file and adds it if it does not
  !! exist yet
  subroutine check_add_variable( this, variable_name, variable_units,         &
      variable_type, append_dimension, dimensions, dimension_sizes, varid,    &
      dimids, start_ids )

    use musica_assert,                 only : assert_msg
    use musica_string,                 only : string_t, to_char
    use netcdf,                        only : nf90_inq_varid,                 &
                                              nf90_inquire_dimension,         &
                                              nf90_inquire_variable,          &
                                              nf90_def_var,                   &
                                              nf90_put_att,                   &
                                              NF90_NOERR, NF90_UNLIMITED

    class(io_netcdf_t), intent(inout) :: this
    type(string_t),     intent(in)    :: variable_name
    type(string_t),     intent(in)    :: variable_units
    integer,            intent(in)    :: variable_type
    type(string_t),     intent(in)    :: append_dimension
    type(string_t),     intent(in)    :: dimensions(:)
    integer,            intent(in)    :: dimension_sizes(:)
    integer,            intent(out)   :: varid
    integer,            intent(out)   :: dimids(size(dimensions)+1)
    integer,            intent(out)   :: start_ids(size(dimensions)+1)

    integer :: ierr, i_dim, ndims, ldimids(size(dimensions)+1)
    type(string_t) :: id_str

    id_str = "variable '"//variable_name//"' in file '"//this%file_name_//"'"

    dimids = this%check_add_dimension( trim( append_dimension%to_char( ) ),   &
                                       NF90_UNLIMITED )
    do i_dim = 1, size( dimensions )
      dimids( i_dim + 1 ) =                                                   &
          this%check_add_dimension( trim( dimensions( i_dim )%to_char( ) ),   &
                                    dimension_sizes( i_dim ) )
    end do
    start_ids = 1
    ierr = nf90_inq_varid( this%file_id_, variable_name%to_char( ), varid )
    if( ierr == NF90_NOERR ) then
      ! Check the dimension ids and units
      call check_status( 372537549,                                           &
                         nf90_inquire_variable( this%file_id_, varid,         &
                                                ndims = ndims,                &
                                                dimids = ldimids ),           &
                        "NetCDF file error for "//trim( id_str%to_char( ) ) )
      call assert_msg( 192756621, ndims == size( dimensions ) + 1,            &
                       "Dimension mismatch for "//trim( id_str%to_char( ) ) )
      do i_dim = 1, size( dimids )
        call assert_msg( 900541544, dimids( i_dim ) == ldimids( i_dim ),      &
                         "Dimension "//trim( to_char( i_dim ) )//             &
                         " mismatch for "//trim( id_str%to_char( ) ) )
      end do
      call check_status( 302003316,                                           &
                         nf90_inquire_dimension( this%file_id_, dimids(1),    &
                                                 len = start_ids(1) ),        &
                         "NetCDF file error for "//trim( id_str%to_char( ) ) )
      start_ids(1) = start_ids(1) + 1
    else
      call check_status( 497577165,                                           &
                         nf90_def_var( this%file_id_,                         &
                                       variable_name%to_char( ),              &
                                       variable_type, dimids, varid ),        &
                        "NetCDF file error for "//trim( id_str%to_char( ) ) )
      call check_status( 757618738,                                           &
                         nf90_put_att( this%file_id_, varid, "units",         &
                                       variable_units%to_char( ) ),           &
                         "Error setting units for "//                         &
                         trim( id_str%to_char( ) ) )
    end if

  end subroutine check_add_variable

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Finalizes a NetCDF file reader
  subroutine finalize( this )

    use netcdf,                        only : nf90_close

    type(io_netcdf_t), intent(inout) :: this

    if( this%file_id_ .ne. kUnknownFileId ) then
      call check_status( 708311006, nf90_close( this%file_id_ ),              &
                         "Error closing file" )
    end if
    this%file_id_   = kUnknownFileId
    this%file_name_ = ""

  end subroutine finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! @name Private NetCDF support functions
!! @{
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Checks a NetCDF status code and fail with a message if an error occurred
  subroutine check_status( code, status, error_message )

    use musica_assert,                 only : die_msg
    use netcdf,                        only : NF90_NOERR, nf90_strerror

    !> Unique code to associate with any failure
    integer,          intent(in) :: code
    !> NetCDF status code
    integer,          intent(in) :: status
    !> Error message to display on failure
    character(len=*), intent(in) :: error_message

    if( status .eq. NF90_NOERR ) return
    call die_msg( 330311277, "NetCDF error: "//trim( error_message )//": "//  &
                  trim( nf90_strerror( status ) ) )

  end subroutine check_status

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> @}

end module musica_io_netcdf
