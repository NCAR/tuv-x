! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0

module tuvx_radiator_warehouse
! A class that holds and gives out
! :f:type:`~tuvx_radiator/radiator_t`'s built by the
! :f:mod:`~tuvx_radiator_factory`.

  use musica_string,                   only : string_t
  use tuvx_radiator,                   only : radiator_ptr

  implicit none

  private
  public :: radiator_warehouse_t, warehouse_iterator_t, radiator_warehouse_ptr

  type radiator_warehouse_t
    ! Radiator warehouse

    private
    type(radiator_ptr), allocatable :: radiators_(:) ! Radiators
    type(string_t),     allocatable :: handle_(:) ! Radiator "handles"
  contains
    !> @name Returns a pointer to a requested radiator
    !! @{
    procedure, private :: get_radiator_char, get_radiator_string,             &
                          get_radiator_ptr, get_radiator_iterator
    generic   :: get_radiator => get_radiator_char, get_radiator_string,      &
                                 get_radiator_ptr, get_radiator_iterator
    !> @}
    !> Returns a pointer to a radiator in the warehouse
    procedure, private :: get_ptr_char, get_ptr_string
    generic :: get_ptr => get_ptr_char, get_ptr_string
    !> Returns whether a radiator exists in the warehouse
    procedure, private :: exists_char, exists_string
    generic :: exists => exists_char, exists_string
    !> Adds a radiator or set of radiators to the warehouse
    procedure, private :: add_radiator, add_radiators
    generic :: add => add_radiator, add_radiators
    !> Returns an updater for a `radiator_from_host_t` radiator
    procedure :: get_updater
    !> Returns the name for a radiator from an iterator
    procedure :: name => get_name
    !> Gets an iterator for the warehouse
    procedure :: get_iterator
    !> Accumulates the state of all radiators in the warehouse
    procedure :: accumulate_states
    !> Returns the number of bytes required to pack the warehouse onto a binary buffer
    procedure :: pack_size
    !> Packs the warehouse onto a character buffer
    procedure :: mpi_pack
    !> Unpacks a warehouse from a character buffer
    procedure :: mpi_unpack
    !> Cleans up memory
    final     :: finalize
  end type radiator_warehouse_t

  type warehouse_iterator_t
    !  Radiator warehouse iterator

    type(radiator_warehouse_t), pointer :: warehouse_ => null( ) ! Pointer to the radiator warehouse
    integer :: id_ = 0 ! Current index in the data set
  contains
    !> Advances to the next key-value pair
    procedure :: next => iterator_next
    !> Resets the iterator
    procedure :: reset => iterator_reset
  end type warehouse_iterator_t

  interface radiator_warehouse_t
    ! constructor of an empty radiator_warehouse_t
    module procedure :: constructor_empty
    ! radiator_warehouse_t constructor
    module procedure :: constructor
  end interface

  !> Pointer to a radiatorin the warehouse
  type :: radiator_warehouse_ptr
    private
    integer :: index_ = 0
  contains
    !> Returns the number of bytes required to pack the pointer onto a buffer
    procedure :: pack_size => ptr_pack_size
    !> Packs the pointer onto a character buffer
    procedure :: mpi_pack => ptr_mpi_pack
    !> Unpacks a pointer from a character buffer into the object
    procedure :: mpi_unpack => ptr_mpi_unpack
  end type radiator_warehouse_ptr

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function constructor_empty( ) result( radiator_warehouse )
    ! Constructs an empty radiator warehouse

    class(radiator_warehouse_t), pointer :: radiator_warehouse

    allocate( radiator_warehouse )
    allocate( radiator_warehouse%radiators_(0) )
    allocate( radiator_warehouse%handle_(0) )

  end function constructor_empty

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function constructor( config, grid_warehouse, profile_warehouse,            &
     cross_section_warehouse ) result( radiator_warehouse )
    ! Constructs radiator_warehouse_t abjects

    use musica_config,                 only : config_t
    use musica_iterator,               only : iterator_t
    use musica_string,                 only : string_t
    use tuvx_cross_section_warehouse,  only : cross_section_warehouse_t
    use tuvx_grid_warehouse,           only : grid_warehouse_t
    use tuvx_profile_warehouse,        only : profile_warehouse_t
    use tuvx_radiator_factory,         only : radiator_builder

    type(config_t),                  intent(inout) :: config ! Radiator configuration
    type(grid_warehouse_t),          intent(inout) :: grid_warehouse ! A :f:type:`~tuvx_grid_warehouse/grid_warehouse_t`
    type(profile_warehouse_t),       intent(inout) :: profile_warehouse ! profile warehouse
    type(cross_section_warehouse_t), intent(inout) :: cross_section_warehouse ! cross section warehouse
    class(radiator_warehouse_t),     pointer       :: radiator_warehouse ! A :f:type:`~tuvx_radiator_warehouse/radiator_warehouse_t`

    ! local variables
    character(len=*), parameter     :: Iam = "Radiator warehouse constructor"
    type(config_t)                  :: radiator_config
    class(iterator_t), pointer      :: iter
    type(radiator_ptr), allocatable :: temp_rads(:)
    type(radiator_ptr)              :: aRadiator
    type(string_t), allocatable     :: temp_names(:)

    allocate( radiator_warehouse )
    allocate( radiator_warehouse%radiators_(0) )
    allocate( radiator_warehouse%handle_(0) )

    iter => config%get_iterator()
    do while( iter%next() )
      call config%get( iter, radiator_config, Iam )

      ! build and store the radiator
      aRadiator%val_ => radiator_builder( radiator_config, grid_warehouse,    &
                                          profile_warehouse,                  &
                                          cross_section_warehouse )

      temp_names = radiator_warehouse%handle_
      deallocate( radiator_warehouse%handle_ )
      allocate( radiator_warehouse%handle_( size( temp_names ) + 1 ) )
      radiator_warehouse%handle_( 1 : size( temp_names ) ) = temp_names(:)
      radiator_warehouse%handle_( size( temp_names ) + 1 ) =                  &
          aRadiator%val_%handle_
      deallocate( temp_names )

      temp_rads = radiator_warehouse%radiators_
      deallocate( radiator_warehouse%radiators_ )
      allocate( radiator_warehouse%radiators_( size( temp_rads ) + 1 ) )
      radiator_warehouse%radiators_( 1 : size( temp_rads ) ) = temp_rads(:)
      radiator_warehouse%radiators_( size( temp_rads ) + 1 ) = aRadiator
      deallocate( temp_rads )
      nullify( aRadiator%val_ )
    end do
    deallocate( iter )

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function get_radiator_char( this, name ) result( radiator )
    ! Returns a pointer to a requested radiator

    use tuvx_radiator,                 only : radiator_t

    class(radiator_warehouse_t), intent(in) :: this ! This :f:type:`~tuvx_radiator_warehouse/radiator_warehouse_t`
    character(len=*),            intent(in) :: name ! Name associated with requested radiator
    class(radiator_t),           pointer    :: radiator ! Pointer to the requested radiator of type :f:type:`~tuvx_radiator/radiator_t`

    type(radiator_warehouse_ptr) :: ptr

    ptr = this%get_ptr_char( name )
    radiator => this%radiators_( ptr%index_ )%val_

  end function get_radiator_char

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function get_radiator_string( this, name ) result( radiator )
    ! Returns a pointer to a radiator

    use musica_string,                 only : string_t
    use tuvx_radiator,                 only : radiator_t

    class(radiator_warehouse_t), intent(in) :: this ! this radiator warehouse
    type(string_t),              intent(in) :: name ! name of the radiator
    class(radiator_t),           pointer    :: radiator

    radiator => this%get_radiator_char( name%to_char( ) )

  end function get_radiator_string

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function get_radiator_ptr( this, ptr ) result( radiator )
    ! Returns a pointer to a radiator

    use musica_assert,                 only : assert_msg
    use tuvx_radiator,                 only : radiator_t

    class(radiator_warehouse_t),  intent(in) :: this     ! this radiator warehouse
    type(radiator_warehouse_ptr), intent(in) :: ptr      ! pointer to the raditor in the warehouse
    class(radiator_t),            pointer    :: radiator ! pointer to the radiator

    call assert_msg( 432301788, ptr%index_ > 0, "Invalid radiator pointer" )
    radiator => this%radiators_( ptr%index_ )%val_

  end function get_radiator_ptr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function get_radiator_iterator( this, iterator ) result( radiator )
    ! Returns a pointer to a radiator in the warehouse from an iterator

    use tuvx_radiator,                 only : radiator_t

    class(radiator_warehouse_t), intent(inout) :: this ! This :f:type:`~tuvx_radiator_warehouse/radiator_warehouse_t`
    type(warehouse_iterator_t),  intent(in)    :: iterator ! ! A :f:type:`~tuvx_radiator_warehouse/warehouse_iterator_t`
    class(radiator_t),           pointer       :: radiator ! Pointer to requested :f:type:`~tuvx_radiator/radiator_t`

    ! Local variables
    character(len=*), parameter :: Iam =                                      &
        'radiator warehouse get radiator from iterator'

    radiator => this%radiators_( iterator%id_ )%val_

  end function get_radiator_iterator

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function get_ptr_char( this, name ) result( ptr )
    ! Returns a pointer to a radiator in the warehouse

    use musica_assert,                 only : assert_msg

    class(radiator_warehouse_t), intent(in) :: this
    character(len=*),            intent(in) :: name
    type(radiator_warehouse_ptr)            :: ptr

    integer :: ndx
    logical :: found

    found = .false.
    do ndx = 1, size( this%handle_ )
      if( name .eq. this%handle_( ndx ) ) then
        found = .true.
        exit
      endif
    end do
    call assert_msg( 200578546, found,                                        &
                     "Invalid radiator handle: '"//name//"'" )
    ptr%index_ = ndx

  end function get_ptr_char

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function get_ptr_string( this, name ) result( ptr )

    use musica_string,                 only : string_t

    class(radiator_warehouse_t), intent(in) :: this
    type(string_t),              intent(in) :: name
    type(radiator_warehouse_ptr)            :: ptr

    ptr = this%get_ptr_char( name%to_char( ) )

  end function get_ptr_string

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function exists_char( this, name )
    ! Returns whether a radiator exists in the warehouse

    class(radiator_warehouse_t), intent(inout) :: this ! This :f:type:`~tuvx_radiator_warehouse/radiator_warehouse_t`
    character(len=*),            intent(in)    :: name ! Name associated with requested radiator
    logical                                    :: exists_char ! Flag indicating whether the radiator exists in the warehouse

    ! Local variables
    character(len=*), parameter :: Iam = 'radiator in warehouse'
    integer :: ndx

    exists_char = .false.
    do ndx = 1, size( this%handle_ )
      if( name == this%handle_( ndx ) ) then
        exists_char = .true.
        exit
      endif
    end do

  end function exists_char

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function exists_string( this, name )
    ! Returns whether a radiator exists in the warehouse

    use musica_string,      only : string_t

    class(radiator_warehouse_t), intent(inout) :: this
    type(string_t),              intent(in)    :: name
    logical                                    :: exists_string

    exists_string = this%exists_char( name%to_char( ) )

  end function exists_string

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine add_radiator( this, radiator )
    ! adds a radiator to the warehouse

    use musica_assert,                 only : assert, assert_msg
    use tuvx_radiator,                 only : radiator_t

    class(radiator_warehouse_t), intent(inout) :: this
    class(radiator_t),           intent(in)    :: radiator

    type(radiator_ptr) :: ptr

    call assert( 439108556, allocated( this%radiators_ ) )
    call assert_msg( 781327898,                                               &
                     .not. this%exists( radiator%handle_ ),                   &
                     "Radiator '"//radiator%handle_//"' already exists." )
    allocate( ptr%val_, source = radiator )
    this%radiators_ = [ this%radiators_, ptr ]
    this%handle_    = [ this%handle_, radiator%handle_ ]

  end subroutine add_radiator

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine add_radiators( this, radiators )
    ! adds a set of radiators to the warehouse

    use musica_assert,                 only : assert

    class(radiator_warehouse_t), intent(inout) :: this
    class(radiator_warehouse_t), intent(in)    :: radiators

    integer :: i_radiator

    call assert( 330601279, allocated( this%radiators_ ) )
    call assert( 777969125, allocated( radiators%radiators_ ) )
    do i_radiator = 1, size( radiators%radiators_ )
      call this%add_radiator( radiators%radiators_( i_radiator )%val_ )
    end do

  end subroutine add_radiators

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  type(radiator_updater_t) function get_updater( this, radiator, found )      &
      result( updater )
    ! returns an updater for a `radiator_from_host_t` radiator
    !
    ! If the optional `found` argument is omitted, an error is returned if the
    ! radiator does not exist in the warehouse.

    use musica_assert,                 only : assert, assert_msg, die_msg
    use tuvx_radiator,                 only : radiator_t
    use tuvx_radiator_from_host,       only : radiator_from_host_t,           &
                                              radiator_updater_t

    class(radiator_warehouse_t), intent(in)  :: this     ! radiator warehouse
    class(radiator_t),           intent(in)  :: radiator ! the radiator to find in the warehouse
    logical, optional,           intent(out) :: found    ! flag indicating whether the
                                                         ! radiator was found

    integer :: i_radiator
    logical :: l_found

    l_found = .false.
    do i_radiator = 1, size( this%radiators_ )
      if( radiator%handle_ == this%radiators_( i_radiator )%val_%handle_ ) then
        l_found = .true.
        exit
      end if
    end do

    if( present( found ) ) then
      found = l_found
      if( .not. found ) return
    end if

    call assert_msg( 995888269, l_found,                                      &
                     "Cannot find radiator '"//radiator%handle_//"'" )

    select type( w_radiator => this%radiators_( i_radiator )%val_ )
    class is( radiator_from_host_t )
      updater = radiator_updater_t( w_radiator )
    class default
      call die_msg( 208206615, "Cannot update radiator '"//w_radiator%handle_ &
                               //"' from a host application." )
    end select

  end function get_updater

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  type(string_t) function get_name( this, iterator )
    ! Returns the name of a radiator from an iterator

    class(radiator_warehouse_t), intent(in) :: this
    class(warehouse_iterator_t), intent(in) :: iterator

    get_name = this%handle_( iterator%id_ )

  end function get_name

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function get_iterator( this )
    ! Gets an iterator for the radiator warehouse

    use musica_assert,                 only : assert

    class(warehouse_iterator_t), pointer            :: get_iterator ! Pointer to the :f:type:`~tuvx_radiator_warehouse/warehouse_iterator_t`
    class(radiator_warehouse_t), intent(in), target :: this ! This :f:type:`~tuvx_radiator_warehouse/radiator_warehouse_t`

    call assert( 753334333, allocated( this%radiators_ ) )
    allocate( warehouse_iterator_t :: get_iterator )
    select type( iter => get_iterator )
      type is( warehouse_iterator_t )
        iter%warehouse_ => this
    end select

  end function get_iterator

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine accumulate_states( this, state )
    ! Accumulates the states off all radiators in the warehouse into a
    ! single representative state.

    use tuvx_radiator,                 only : radiator_state_t

    class(radiator_warehouse_t), intent(in)    :: this
    class(radiator_state_t),     intent(inout) :: state

    call state%accumulate( this%radiators_ )

  end subroutine accumulate_states

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function pack_size( this, comm )
    ! Returns the size of a character buffer required to pack the warehouse

    use musica_assert,                 only : assert
    use musica_mpi,                    only : musica_mpi_pack_size
    use tuvx_radiator_factory,         only : radiator_type_name

    class(radiator_warehouse_t), intent(in) :: this ! warehouse to be packed
    integer,                     intent(in) :: comm ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: i_radiator
    type(string_t) :: type_name

    call assert( 208683985, allocated( this%radiators_ ) )
    call assert( 263163771, allocated( this%handle_ ) )
    call assert( 710531617,                                                   &
                 size( this%radiators_ ) .eq. size( this%handle_ ) )
    pack_size = musica_mpi_pack_size( size( this%radiators_ ), comm )
    do i_radiator = 1, size( this%radiators_ )
    associate( radiator => this%radiators_( i_radiator )%val_ )
      type_name = radiator_type_name( radiator )
      pack_size = pack_size +                                                 &
                  type_name%pack_size( comm ) +                               &
                  radiator%pack_size( comm ) +                                &
                  this%handle_( i_radiator )%pack_size( comm )
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
    use tuvx_radiator_factory,         only : radiator_type_name

    class(radiator_warehouse_t), intent(in)    :: this      ! warehouse to be packed
    character,                   intent(inout) :: buffer(:) ! memory buffer
    integer,                     intent(inout) :: position  ! currently buffer position
    integer,                     intent(in)    :: comm      ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: prev_pos, i_radiator
    type(string_t) :: type_name

    prev_pos = position
    call assert( 629181647, allocated( this%radiators_ ) )
    call assert( 459024743, allocated( this%handle_ ) )
    call assert( 906392589,                                                   &
                 size( this%radiators_ ) .eq. size( this%handle_ ) )
    call musica_mpi_pack( buffer, position, size( this%radiators_ ), comm )
    do i_radiator = 1, size( this%radiators_ )
    associate( radiator => this%radiators_( i_radiator )%val_ )
      type_name = radiator_type_name( radiator )
      call type_name%mpi_pack( buffer, position, comm )
      call radiator%mpi_pack( buffer, position, comm )
      call this%handle_( i_radiator )%mpi_pack( buffer, position, comm )
    end associate
    end do
    call assert( 463435654, position - prev_pos <= this%pack_size( comm ) )
#endif

  end subroutine mpi_pack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine mpi_unpack( this, buffer, position, comm )
    ! Unpacks a warehouse from a character buffer

    use musica_assert,                 only : assert
    use musica_mpi,                    only : musica_mpi_unpack
    use tuvx_radiator_factory,         only : radiator_allocate

    class(radiator_warehouse_t), intent(out)   :: this      ! warehouse to unpack
    character,                   intent(inout) :: buffer(:) ! memory buffer
    integer,                     intent(inout) :: position  ! current buffer position
    integer,                     intent(in)    :: comm      ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: prev_pos, i_radiator, n_radiators
    type(string_t) :: type_name

    prev_pos = position
    call musica_mpi_unpack( buffer, position, n_radiators, comm )
    if( allocated( this%radiators_ ) ) deallocate( this%radiators_ )
    if( allocated( this%handle_   ) ) deallocate( this%handle_   )
    allocate( this%radiators_( n_radiators ) )
    allocate( this%handle_(   n_radiators ) )
    do i_radiator = 1, n_radiators
    associate( radiator => this%radiators_( i_radiator ) )
      call type_name%mpi_unpack( buffer, position, comm )
      radiator%val_ => radiator_allocate( type_name )
      call radiator%val_%mpi_unpack( buffer, position, comm )
      call this%handle_( i_radiator )%mpi_unpack( buffer, position, comm )
    end associate
    end do
    call assert( 928789078, position - prev_pos <= this%pack_size( comm ) )
#endif

  end subroutine mpi_unpack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function iterator_next( this ) result( continue )
    ! Advances the iterator
    !
    ! Returns false if the end of the collection has been reached

    class(warehouse_iterator_t), intent(inout) :: this ! A :f:type:`~tuvx_radiator_warehouse/warehouse_iterator_t`

    logical :: continue

    this%id_ = this%id_ + 1
    continue = this%id_ <= size( this%warehouse_%radiators_ )

  end function iterator_next

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine iterator_reset( this )
    ! Resets the iterator

    class(warehouse_iterator_t), intent(inout) :: this ! This :f:type:`~tuvx_radiator_warehouse/warehouse_iterator_t`

    this%id_ = 0

  end subroutine iterator_reset

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine finalize( this )
    ! Finalize the radiator warehouse

    type(radiator_warehouse_t), intent(inout) :: this ! This :f:type:`~tuvx_radiator_warehouse/radiator_warehouse_t`

    integer :: ndx
    character(len=*), parameter :: Iam = 'radiator_warehouse finalize: '

    if( allocated( this%radiators_ ) ) then
      do ndx = 1, size( this%radiators_ )
        if( associated( this%radiators_( ndx )%val_ ) ) then
          deallocate( this%radiators_( ndx )%val_ )
        end if
      end do
      deallocate( this%radiators_ )
    endif

  end subroutine finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function ptr_pack_size( this, comm ) result( pack_size )
    ! Returns the number of bytes required to pack the pointer onto a buffer

    use musica_mpi,                    only : musica_mpi_pack_size

    class(radiator_warehouse_ptr), intent(in) :: this ! This radiator pointer
    integer,                       intent(in) :: comm ! MPI communicator

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

    class(radiator_warehouse_ptr), intent(in)    :: this      ! pointer to be packed
    character,                     intent(inout) :: buffer(:) ! memory buffer
    integer,                       intent(inout) :: position  ! current buffer position
    integer,                       intent(in)    :: comm      ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: prev_pos

    prev_pos = position
    call musica_mpi_pack( buffer, position, this%index_, comm )
    call assert( 776713787, position - prev_pos <= this%pack_size( comm ) )
#endif

  end subroutine ptr_mpi_pack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ptr_mpi_unpack( this, buffer, position, comm )
    ! Unpacks a pointer from a character buffer

    use musica_assert,                 only : assert
    use musica_mpi,                    only : musica_mpi_unpack

    class(radiator_warehouse_ptr), intent(out)   :: this      ! pointer to be packed
    character,                     intent(inout) :: buffer(:) ! memory buffer
    integer,                       intent(inout) :: position  ! current buffer position
    integer,                       intent(in)    :: comm      ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: prev_pos

    prev_pos = position
    call musica_mpi_unpack( buffer, position, this%index_, comm )
    call assert( 381920193, position - prev_pos <= this%pack_size( comm ) )
#endif

  end subroutine ptr_mpi_unpack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_radiator_warehouse
