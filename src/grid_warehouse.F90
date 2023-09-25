! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0

module tuvx_grid_warehouse
! A warehouse to hold and distribute grids.

  use tuvx_grid, only : grid_ptr

  implicit none

  private
  public :: grid_warehouse_t, grid_warehouse_ptr

  !> Grid warehouse type
  type :: grid_warehouse_t
    private
    type(grid_ptr), allocatable :: grids_(:) ! grid objects
  contains
    !> get a copy of a grid object
    procedure, private :: get_grid_char, get_grid_string, get_grid_ptr
    generic :: get_grid => get_grid_char, get_grid_string, get_grid_ptr
    !> returns a pointer to a grid object
    procedure, private :: get_ptr_char, get_ptr_string
    generic :: get_ptr => get_ptr_char, get_ptr_string
    !> checks if a grid is present in the warehouse
    procedure, private :: exists_char, exists_string
    generic :: exists => exists_char, exists_string
    !> Adds a grid or set of grids to the warehouse
    procedure, private :: add_grid
    procedure, private :: add_grids
    generic :: add => add_grid, add_grids
    !> Returns a updater for a `grid_from_host_t` grid
    procedure :: get_updater
    !> Returns the number of bytes required to pack the warehouse onto a buffer
    procedure :: pack_size
    !> Packs the warehouse onto a character buffer
    procedure :: mpi_pack
    !> Unpacks a warehouse from a character buffer into the object
    procedure :: mpi_unpack
    !> Finalize the object
    final :: finalize
  end type grid_warehouse_t

  !> Grid warehouse_t constructor
  interface grid_warehouse_t
    module procedure :: constructor_empty
    module procedure :: constructor
  end interface

  !> Pointer to a grid in the warehouse
  type :: grid_warehouse_ptr
    private
    integer :: index_ = 0
  contains
    !> Returns the number of bytes required to pack the pointer onto a buffer
    procedure :: pack_size => ptr_pack_size
    !> Packs the pointer onto a character buffer
    procedure :: mpi_pack => ptr_mpi_pack
    !> Unpacks a pointer from a character buffer into the object
    procedure :: mpi_unpack => ptr_mpi_unpack
  end type grid_warehouse_ptr

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function constructor_empty( ) result( grid_warehouse )

    class(grid_warehouse_t), pointer :: grid_warehouse ! Empty grid warehouse

    allocate( grid_warehouse )
    allocate( grid_warehouse%grids_(0) )

  end function constructor_empty

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function constructor( config ) result( grid_warehouse )
    ! Grid warehouse constructor

    use musica_assert,                 only : assert_msg
    use musica_config,                 only : config_t
    use musica_iterator,               only : iterator_t
    use musica_string,                 only : string_t
    use tuvx_grid_factory,             only : grid_builder

    type(config_t), intent(inout) :: config ! grid configuration data
    class(grid_warehouse_t), pointer :: grid_warehouse ! This :f:type:`~tuvx_grid_warehouse/grid_warehouse_t`

    !> local variables
    character(len=*), parameter :: Iam = "Grid warehouse constructor: "
    type(config_t)              :: grid_config
    class(iterator_t), pointer  :: iter
    type(grid_ptr), allocatable :: temp_grids(:)
    type(grid_ptr)              :: grid_obj
    character(len=32)           :: keychar
    type(string_t)              :: aswkey

    allocate( grid_warehouse )
    allocate( grid_warehouse%grids_(0) )

    ! iterate over grids
    iter => config%get_iterator()
    do while( iter%next() )
      call config%get( iter, grid_config, Iam )

      ! Build grid objects
      grid_obj%val_ => grid_builder( grid_config )
      call assert_msg( 101630104,                                             &
                       .not. grid_warehouse%exists( grid_obj%val_%handle_,    &
                                                    grid_obj%val_%units( ) ), &
                       "Grid '"//grid_obj%val_%handle_//                      &
                       "' duplicated in grid warehouse." )
      temp_grids = grid_warehouse%grids_
      deallocate( grid_warehouse%grids_ )
      allocate( grid_warehouse%grids_( size( temp_grids ) + 1 ) )
      grid_warehouse%grids_( 1 : size( temp_grids ) ) = temp_grids(:)
      grid_warehouse%grids_( size( temp_grids ) + 1 ) = grid_obj
      deallocate( temp_grids )
      nullify( grid_obj%val_ )
    end do
    deallocate( iter )

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function get_grid_char( this, name, units ) result( a_grid_ptr )
    ! Get copy of a grid object

    use musica_assert,                 only : assert_msg, die_msg
    use musica_string,                 only : string_t
    use tuvx_grid,                     only : grid_t

    class(grid_warehouse_t), intent(in) :: this     ! This :f:type:`~tuvx_grid_warehouse/grid_warehouse_t`
    character(len=*),        intent(in) :: name     ! The name of a grid, see :ref:`configuration-grids` for grid names
    character(len=*),        intent(in) :: units    ! The units of the grid
    class(grid_t), pointer              :: a_grid_ptr ! The :f:type:`~tuvx_grid/grid_t` which matches the name passed in

    type(grid_warehouse_ptr) :: ptr

    ptr = this%get_ptr_char( name, units )
    allocate( a_grid_ptr, source = this%grids_( ptr%index_ )%val_ )

  end function get_grid_char

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function get_grid_string( this, name, units ) result( a_grid_ptr )
    ! Get a copy of a grid object

    use musica_string,                 only : string_t
    use tuvx_grid,                     only : grid_t

    class(grid_warehouse_t), intent(inout) :: this ! This :f:type:`~tuvx_grid_warehouse/grid_warehouse_t`
    type(string_t),          intent(in)    :: name ! The name of a grid, see :ref:`configuration-grids` for grid names
    type(string_t),          intent(in)    :: units ! The units of the grid
    class(grid_t), pointer                 :: a_grid_ptr ! The :f:type:`~tuvx_grid/grid_t` which matches the name passed in

    a_grid_ptr => this%get_grid_char( name%to_char( ), units%to_char( ) )

  end function get_grid_string

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function get_grid_ptr( this, ptr ) result( grid )
    ! Returns a copy of a grid from a grid pointer

    use musica_assert,                 only : assert_msg
    use tuvx_grid,                     only : grid_t

    class(grid_warehouse_t),  intent(inout) :: this ! This grid warehouse
    type(grid_warehouse_ptr), intent(in)    :: ptr  ! Pointer to a grid in the warehouse
    class(grid_t),            pointer       :: grid

    call assert_msg( 870082797, ptr%index_ > 0, "Invalid grid pointer" )
    allocate( grid, source = this%grids_( ptr%index_ )%val_ )

  end function get_grid_ptr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function get_ptr_char( this, name, units ) result( ptr )
    ! Returns a pointer to a grid object in the warehouse

    use musica_assert,                 only : assert_msg, die_msg
    use musica_string,                 only : string_t

    class(grid_warehouse_t), intent(in) :: this  ! This :f:type:`~tuvx_grid_warehouse/grid_warehouse_t`
    character(len=*),        intent(in) :: name  ! The name of a grid, see :ref:`configuration-grids` for grid names
    character(len=*),        intent(in) :: units ! The units of the grid
    type(grid_warehouse_ptr)            :: ptr   ! Pointer to the grid in the warehouse

    integer :: ndx
    logical :: found
    type(string_t) :: this_units

    found = .false.
    do ndx = 1, size( this%grids_ )
      if( name .eq. this%grids_( ndx )%val_%handle_ ) then
        found = .true.
        exit
      endif
    end do
    call assert_msg( 345804219, found, "Invalid grid name: '"//name//"'" )
    this_units = this%grids_( ndx )%val_%units( )
    if( units .ne. this_units ) then
      call die_msg( 509243577,                                                &
                    "Grid '"//name//"' has units of '"//this_units//          &
                    "' not '"//units//"' as requested." )
    end if
    ptr%index_ = ndx

  end function get_ptr_char

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function get_ptr_string( this, name, units ) result( ptr )
    ! Returns a pointer to a grid object in the warehouse

    use musica_string,                 only : string_t
    use tuvx_grid,                     only : grid_t

    class(grid_warehouse_t), intent(inout) :: this  ! This :f:type:`~tuvx_grid_warehouse/grid_warehouse_t`
    type(string_t),          intent(in)    :: name  ! The name of a grid, see :ref:`configuration-grids` for grid names
    type(string_t),          intent(in)    :: units ! The units of the grid
    type(grid_warehouse_ptr)               :: ptr   ! Pointer to the grid in the warehouse

    ptr = this%get_ptr_char( name%to_char( ), units%to_char( ) )

  end function get_ptr_string

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  logical function exists_char( this, name, units ) result( exists )
    ! checks if a grid exists in the warehouse

    use musica_assert,                 only : assert_msg

    class(grid_warehouse_t), intent(inout) :: this ! This :f:type:`~tuvx_grid_warehouse/grid_warehouse_t`
    character(len=*),        intent(in)    :: name ! The name of a grid, see :ref:`configuration-grids` for grid names
    character(len=*),        intent(in)    :: units ! The units of the grid

    integer :: ndx

    exists = .false.
    do ndx = 1, size( this%grids_ )
      if( name .eq. this%grids_( ndx )%val_%handle_ ) then
        call assert_msg( 287595102, this%grids_( ndx )%val_%units( ) == units,&
                         "Units mismatch for grid '"//name//"': '"//units//   &
                         "' != '"//this%grids_( ndx )%val_%units( ) )
        exists = .true.
        exit
      endif
    end do

  end function exists_char

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  logical function exists_string( this, name, units ) result( exists )
    ! checks if a grid exists in the warehouse

    use musica_string,                 only : string_t

    class(grid_warehouse_t), intent(inout) :: this ! This :f:type:`~tuvx_grid_warehouse/grid_warehouse_t`
    type(string_t),          intent(in)    :: name ! The name of a grid, see :ref:`configuration-grids` for grid names
    type(string_t),          intent(in)    :: units ! The units of the grid

    exists = this%exists_char( name%to_char( ), units%to_char( ) )

  end function exists_string

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine add_grid( this, grid )
    ! adds a grid to the warehouse

    use musica_assert,                 only : assert, assert_msg
    use tuvx_grid,                     only : grid_t

    class(grid_warehouse_t), intent(inout) :: this ! This :f:type:`~tuvx_grid_warehouse/grid_warehouse_t`
    class(grid_t),           intent(in)    :: grid ! Grid to add

    type(grid_ptr) :: ptr

    call assert( 900933280, allocated( this%grids_  ) )
    call assert_msg( 244177406,                                               &
                     .not. this%exists( grid%handle_, grid%units( ) ),        &
                     "Grid '"//grid%handle_//"' already exists." )
    allocate( ptr%val_, source = grid )
    this%grids_ = [ this%grids_, ptr ]

  end subroutine add_grid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine add_grids( this, grids )
    ! adds a set of grids to the warehouse

    use musica_assert,                 only : assert

    class(grid_warehouse_t), intent(inout) :: this ! This :f:type:`~tuvx_grid_warehouse/grid_warehouse_t`
    class(grid_warehouse_t), intent(in)    :: grids ! Set of grids to add

    integer :: i_grid

    call assert( 169184651, allocated( this%grids_  ) )
    call assert( 176354492, allocated( grids%grids_ ) )
    do i_grid = 1, size( grids%grids_ )
      call this%add_grid( grids%grids_( i_grid )%val_ )
    end do

  end subroutine add_grids

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  type(grid_updater_t) function get_updater( this, grid, found )              &
      result( updater )
    ! returns an updater for a `grid_from_host_t` grid
    !
    ! If the optional `found` flag is omitted, an error is returned if the
    ! grid does not exist in the warehouse.

    use musica_assert,                 only : assert_msg, die_msg
    use tuvx_grid,                     only : grid_t
    use tuvx_grid_from_host,           only : grid_from_host_t, grid_updater_t

    class(grid_warehouse_t), intent(in)  :: this ! This :f:type:`~tuvx_grid_warehouse/grid_warehouse_t`
    class(grid_t),           intent(in)  :: grid ! The grid to find in the warehouse
    logical, optional,       intent(out) :: found ! Flag indicating whether the grid was found

    integer :: i_grid
    logical :: l_found

    l_found = .false.
    do i_grid = 1, size( this%grids_ )
      if( grid%handle_ == this%grids_( i_grid )%val_%handle_ ) then
        call assert_msg( 782725188,                                           &
                         this%grids_( i_grid )%val_%units( ) == grid%units( ),&
                         "Units mismatch for grid '"//grid%handle_//"': '"//  &
                         grid%units( )//"' != '"//                            &
                         this%grids_( i_grid )%val_%units( ) )
        l_found = .true.
        exit
      end if
    end do

    if( present( found ) ) then
      found = l_found
      if( .not. found ) return
    end if

    call assert_msg( 311845931, l_found,                                      &
                     "Cannot find grid '"//grid%handle_//"'" )

    select type( w_grid => this%grids_( i_grid )%val_ )
    class is( grid_from_host_t )
      updater = grid_updater_t( w_grid )
    class default
      call die_msg( 953621510, "Cannot update grid '"//w_grid%handle_//       &
                               "' from a host application" )
    end select

  end function get_updater

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function pack_size( this, comm )
    ! Returns the number of bytes required to pack the warehouse onto a buffer

    use musica_assert,                 only : assert
    use musica_mpi,                    only : musica_mpi_pack_size
    use musica_string,                 only : string_t
    use tuvx_grid_factory,             only : grid_type_name

    class(grid_warehouse_t), intent(in) :: this ! warehouse to be packed
    integer,                 intent(in) :: comm ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: i_grid
    type(string_t) :: type_name

    call assert( 880887280, allocated( this%grids_ ) )
    pack_size = musica_mpi_pack_size( size( this%grids_ ), comm )
    do i_grid = 1, size( this%grids_ )
    associate( grid => this%grids_( i_grid )%val_ )
      type_name = grid_type_name( grid )
      pack_size = pack_size +                                                 &
                  type_name%pack_size( comm ) +                               &
                  grid%pack_size( comm )
    end associate
    end do
#else
    pack_size = 0
#endif

  end function pack_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine mpi_pack( this, buffer, position, comm )
    ! Packs the warehouse onto a character buffer

    use musica_assert,                 only : assert
    use musica_mpi,                    only : musica_mpi_pack
    use musica_string,                 only : string_t
    use tuvx_grid_factory,             only : grid_type_name

    class(grid_warehouse_t), intent(in)    :: this      ! warehouse to be packed
    character,               intent(inout) :: buffer(:) ! memory buffer
    integer,                 intent(inout) :: position  ! current buffer position
    integer,                 intent(in)    :: comm      ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: prev_pos, i_grid
    type(string_t) :: type_name

    prev_pos = position
    call assert( 607322234, allocated( this%grids_ ) )
    call musica_mpi_pack( buffer, position, size( this%grids_ ), comm )
    do i_grid = 1, size( this%grids_ )
    associate( grid => this%grids_( i_grid )%val_ )
      type_name = grid_type_name( grid )
      call type_name%mpi_pack( buffer, position, comm )
      call grid%mpi_pack(      buffer, position, comm )
    end associate
    end do
    call assert( 379779066, position - prev_pos <= this%pack_size( comm ) )
#endif

  end subroutine mpi_pack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine mpi_unpack( this, buffer, position, comm )
    ! Unpacks a warehouse from a character buffer into the object

    use musica_assert,                 only : assert
    use musica_mpi,                    only : musica_mpi_unpack
    use musica_string,                 only : string_t
    use tuvx_grid_factory,             only : grid_allocate

    class(grid_warehouse_t), intent(out)   :: this      ! warehouse to be unpacked
    character,               intent(inout) :: buffer(:) ! memory buffer
    integer,                 intent(inout) :: position  ! current buffer position
    integer,                 intent(in)    :: comm      ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: prev_pos, i_grid, n_grids
    type(string_t) :: type_name

    prev_pos = position
    call musica_mpi_unpack( buffer, position, n_grids, comm )
    if( allocated( this%grids_ ) ) deallocate( this%grids_ )
    allocate( this%grids_( n_grids ) )
    do i_grid = 1, n_grids
    associate( grid => this%grids_( i_grid ) )
      call type_name%mpi_unpack( buffer, position, comm )
      grid%val_ => grid_allocate( type_name )
      call grid%val_%mpi_unpack( buffer, position, comm )
    end associate
    end do
    call assert( 459962920, position - prev_pos <= this%pack_size( comm ) )
#endif

  end subroutine mpi_unpack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine finalize( this )
    ! Finalize grid warehouse

    !> Arguments
    type(grid_warehouse_t), intent(inout) :: this ! This :f:type:`~tuvx_grid_warehouse/grid_warehouse_t`

    integer :: ndx

    if( allocated( this%grids_ ) ) then
      do ndx = 1, size( this%grids_ )
        if( associated( this%grids_( ndx )%val_ ) ) then
          deallocate( this%grids_( ndx )%val_ )
        end if
      end do
      deallocate( this%grids_ )
    endif

  end subroutine finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function ptr_pack_size( this, comm ) result( pack_size )
    ! Returns the number of bytes required to pack the pointer onto a buffer

    use musica_mpi,                    only : musica_mpi_pack_size

    class(grid_warehouse_ptr), intent(in) :: this ! This grid pointer
    integer,                   intent(in) :: comm ! MPI communicator

#ifdef MUSICA_USE_MPI
    pack_size = musica_mpi_pack_size( this%index_, comm )
#else
    pack_size = 0
#endif

  end function ptr_pack_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ptr_mpi_pack( this, buffer, position, comm )
    ! Packs a pointer onto a character buffer

    use musica_assert,                 only : assert
    use musica_mpi,                    only : musica_mpi_pack

    class(grid_warehouse_ptr), intent(in)    :: this      ! pointer to be packed
    character,                 intent(inout) :: buffer(:) ! memory buffer
    integer,                   intent(inout) :: position  ! current buffer position
    integer,                   intent(in)    :: comm      ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: prev_pos

    prev_pos = position
    call musica_mpi_pack( buffer, position, this%index_, comm )
    call assert( 282334709, position - prev_pos <= this%pack_size( comm ) )
#endif

  end subroutine ptr_mpi_pack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ptr_mpi_unpack( this, buffer, position, comm )
    ! Unpacks a pointer from a character buffer

    use musica_assert,                 only : assert
    use musica_mpi,                    only : musica_mpi_unpack

    class(grid_warehouse_ptr), intent(out)   :: this      ! pointer to be packed
    character,                 intent(inout) :: buffer(:) ! memory buffer
    integer,                   intent(inout) :: position  ! current buffer position
    integer,                   intent(in)    :: comm      ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: prev_pos

    prev_pos = position
    call musica_mpi_unpack( buffer, position, this%index_, comm )
    call assert( 212966592, position - prev_pos <= this%pack_size( comm ) )
#endif

  end subroutine ptr_mpi_unpack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_grid_warehouse
