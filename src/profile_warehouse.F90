! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0

module tuvx_profile_warehouse
  ! The profile warehouse type and related functions. Holds one or more
  ! :f:type:`~tuvx_profile/profile_t` s created by the
  ! :f:mod:`tuvx_profile_factory`

  use tuvx_profile, only : profile_ptr

  implicit none

  private
  public :: profile_warehouse_t, profile_warehouse_ptr

  type profile_warehouse_t
    private
    type(profile_ptr), allocatable :: profiles_(:)
  contains
    ! returns a copy of a profile object
    procedure, private :: get_profile_char, get_profile_string, get_profile_ptr
    generic :: get_profile => get_profile_char, get_profile_string,           &
                              get_profile_ptr
    ! returns a pointer to a profile object
    procedure, private :: get_ptr_char, get_ptr_string
    generic :: get_ptr => get_ptr_char, get_ptr_string
    ! checks if a profile is present in the warehouse
    procedure :: exists_char, exists_string
    generic :: exists => exists_char, exists_string
    ! adds a profile or set of profiles to the warehouse
    procedure, private :: add_profile
    procedure, private :: add_profiles
    generic :: add => add_profile, add_profiles
    ! returns an updater for a `profile_from_host_t` profile
    procedure :: get_updater
    ! returns the number of bytes required to pack the warehouse onto a buffer
    procedure :: pack_size
    ! packs the warehouse onto a character buffer
    procedure :: mpi_pack
    ! unpacks a warehouse from a character buffer into the object
    procedure :: mpi_unpack
    final :: finalize
  end type profile_warehouse_t

  interface profile_warehouse_t
    ! empty constructor for profile_warehouse_t
    module procedure :: constructor_empty
    ! profile_warehouse_t constructor
    module procedure :: constructor
  end interface

  !> Pointer to a profile in the warehouse
  type :: profile_warehouse_ptr
    private
    integer :: index_ = 0
  contains
    !> Returns the number of bytes required to pack the pointer onto a buffer
    procedure :: pack_size => ptr_pack_size
    !> Packs the pointer onto a character buffer
    procedure :: mpi_pack => ptr_mpi_pack
    !> Unpacks a pointer from a character buffer into the object
    procedure :: mpi_unpack => ptr_mpi_unpack
  end type profile_warehouse_ptr

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function constructor_empty( ) result( profile_warehouse )

    class(profile_warehouse_t), pointer :: profile_warehouse ! Empty profile warehouse

    allocate( profile_warehouse )
    allocate( profile_warehouse%profiles_(0) )

  end function constructor_empty

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function constructor( config, grid_warehouse ) result( profile_warehouse )
    ! profile warehouse constructor

    use musica_assert,        only : assert_msg
    use musica_config,        only : config_t
    use musica_iterator,      only : iterator_t
    use musica_string,        only : string_t
    use tuvx_grid_warehouse,  only : grid_warehouse_t
    use tuvx_profile_factory, only : profile_builder

    type(config_t),             intent(inout) :: config ! profile configuration data
    type(grid_warehouse_t),     intent(inout) :: grid_warehouse ! A :f:type:`~tuvx_grid_warehouse/grid_warehouse_t`
    class(profile_warehouse_t), pointer       :: profile_warehouse ! A :f:type:`~tuvx_profile_warehouse/profile_warehouse_t`

    ! local variables
    type(config_t)                 :: profile_config
    class(iterator_t), pointer     :: iter
    type(profile_ptr)              :: profile_obj
    character(len=32)              :: keychar
    type(string_t)                 :: aswkey
    character(len=*), parameter    :: Iam = "profile warehouse constructor: "
    integer                        :: profile_idx

    allocate( profile_warehouse )
    allocate( profile_warehouse%profiles_(config%number_of_children()) )

    ! iterate over profiles
    profile_idx = 1
    iter => config%get_iterator( )
    do while( iter%next( ) )
      call config%get( iter, profile_config, Iam )

      ! Build profile objects
      profile_obj%val_ => profile_builder( profile_config, grid_warehouse )

      call assert_msg( 178168062,                                         &
                      .not. profile_warehouse%exists(                     &
                                             profile_obj%val_%handle_,    &
                                             profile_obj%val_%units( ) ), &
                      "Profile '"//profile_obj%val_%handle_//             &
                      "' duplicated in profile warehouse." )

      profile_warehouse%profiles_( profile_idx ) = profile_obj
      nullify( profile_obj%val_ )

      profile_idx = profile_idx + 1
    end do

    deallocate( iter )

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function get_profile_char( this, name, units ) result( a_profile_ptr )
    ! Get copy of a profile object

    use musica_assert,     only : assert_msg
    use tuvx_profile,      only : profile_t

    class(profile_warehouse_t), intent(in) :: this
    character(len=*),           intent(in) :: name
    character(len=*),           intent(in) :: units
    class(profile_t),           pointer    :: a_profile_ptr

    type(profile_warehouse_ptr) :: ptr

    ptr = this%get_ptr( name, units )
    allocate( a_profile_ptr, source = this%profiles_( ptr%index_ )%val_ )

  end function get_profile_char

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function get_profile_string( this, name, units ) result( a_profile_ptr )
    ! Get a copy of a profile object

    use musica_string,                 only : string_t
    use tuvx_profile,                  only : profile_t

    class(profile_warehouse_t), intent(in) :: this
    type(string_t),             intent(in) :: name
    type(string_t),             intent(in) :: units
    class(profile_t), pointer              :: a_profile_ptr

    a_profile_ptr => this%get_profile_char( name%to_char( ), units%to_char( ) )

  end function get_profile_string

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function get_profile_ptr( this, ptr ) result( profile )
    ! Returns a copy of a profile from the warehouse

    use musica_assert,                 only : assert_msg
    use tuvx_profile,                  only : profile_t

    class(profile_warehouse_t),  intent(in) :: this    ! The profile warehouse
    type(profile_warehouse_ptr), intent(in) :: ptr     ! Pointer to a profile in the warehouse
    class(profile_t),            pointer    :: profile ! Copy of the profile

    call assert_msg( 419668898, ptr%index_ > 0, "Invalid profile pointer" )
    allocate( profile, source = this%profiles_( ptr%index_ )%val_ )

  end function get_profile_ptr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function get_ptr_char( this, name, units ) result( ptr )
    ! Returns a pointer to a profile in the warehouse

    use musica_assert,                 only : assert_msg

    class(profile_warehouse_t), intent(in) :: this
    character(len=*),           intent(in) :: name
    character(len=*),           intent(in) :: units
    type(profile_warehouse_ptr)            :: ptr

    integer :: ndx
    logical :: found

    found = .false.
    do ndx = 1, size( this%profiles_ )
      if( name .eq. this%profiles_( ndx )%val_%handle_ ) then
        found = .true.
        exit
      endif
    end do

    call assert_msg( 460768214, found, "Invalid profile handle: '"//name//"'" )
    call assert_msg( 465479557,                                               &
                     units .eq. this%profiles_( ndx )%val_%units( ),          &
                     "Profile '"//name//"' has units of '"//                  &
                     this%profiles_( ndx )%val_%units( )//"' not '"//         &
                     units//"' as requested." )
    ptr%index_ = ndx

  end function get_ptr_char

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function get_ptr_string( this, name, units ) result( ptr )
    ! Returns a pointer to a profile object in the warehouse

    use musica_string,                 only : string_t

    class(profile_warehouse_t), intent(in) :: this
    type(string_t),             intent(in) :: name
    type(string_t),             intent(in) :: units
    type(profile_warehouse_ptr)            :: ptr

    ptr = this%get_ptr_char( name%to_char( ), units%to_char( ) )

  end function get_ptr_string

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  logical function exists_char( this, name, units ) result( exists )
    ! checks if a profile exists in the warehouse

    use musica_assert,                 only : assert_msg

    class(profile_warehouse_t), intent(in) :: this
    character(len=*),           intent(in) :: name
    character(len=*),           intent(in) :: units

    integer :: ndx

    exists = .false.
    do ndx = 1, size( this%profiles_ )
      if ( .not. associated(this%profiles_(ndx)%val_)) cycle

      if( name .eq. this%profiles_( ndx )%val_%handle_ ) then
        call assert_msg( 496262126,                                           &
                         this%profiles_( ndx )%val_%units( ) == units,        &
                         "Units mismatch for profile '"//name//"': '"//units//&
                         "' != '"//this%profiles_( ndx )%val_%units( ) )
        exists  = .true.
        exit
      endif
    end do

  end function exists_char

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  logical function exists_string( this, name, units ) result( exists )
    ! checks if a profile exists in the warehouse

    use musica_string,                 only : string_t

    class(profile_warehouse_t), intent(inout) :: this
    type(string_t),             intent(in)    :: name
    type(string_t),             intent(in)    :: units

    exists = this%exists_char( name%to_char( ), units%to_char( ) )

  end function exists_string

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine add_profile( this, profile )
    ! adds a profile to the warehouse

    use musica_assert,                 only : assert, assert_msg
    use tuvx_profile,                  only : profile_t

    class(profile_warehouse_t), intent(inout) :: this
    class(profile_t),           intent(in)    :: profile

    type(profile_ptr) :: ptr

    call assert( 809705750, allocated( this%profiles_ ) )
    call assert_msg( 490846671,                                               &
                     .not. this%exists( profile%handle_, profile%units( ) ),  &
                     "Profile '"//profile%handle_//"' already exists." )
    allocate( ptr%val_, source = profile )
    this%profiles_ = [ this%profiles_, ptr ]

  end subroutine add_profile

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine add_profiles( this, profiles )
    ! adds a set of profiles to the warehouse

    use musica_assert,                 only : assert

    class(profile_warehouse_t), intent(inout) :: this
    class(profile_warehouse_t), intent(in)    :: profiles

    integer :: i_profile

    call assert( 220530020, allocated( this%profiles_ ) )
    call assert( 162691461, allocated( profiles%profiles_ ) )
    do i_profile = 1, size( profiles%profiles_ )
      call this%add_profile( profiles%profiles_( i_profile )%val_ )
    end do

  end subroutine add_profiles

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  type(profile_updater_t) function get_updater( this, profile, found )        &
      result( updater )
    ! returns an updater for a `profile_from_host_t` profile
    !
    ! If the optional `found` argument is omitted, an error is returned if the
    ! profile does not exist in the warehouse.

    use musica_assert,                 only : assert, assert_msg, die_msg
    use tuvx_profile,                  only : profile_t
    use tuvx_profile_from_host,        only : profile_from_host_t,            &
                                              profile_updater_t

    class(profile_warehouse_t), intent(in)  :: this    ! profile warehouse
    class(profile_t),           intent(in)  :: profile ! the profile to find in the warehouse
    logical, optional,          intent(out) :: found   ! flag indicating whether the profile was found

    integer :: i_profile
    logical :: l_found

    l_found = .false.
    do i_profile = 1, size( this%profiles_ )
      if( profile%handle_ == this%profiles_( i_profile )%val_%handle_ ) then
        call assert_msg( 585094657,                                           &
                         this%profiles_( i_profile )%val_%units( )            &
                         == profile%units( ),                                 &
                         "Units mismatch for profile '"//profile%handle_//    &
                         "': '"//profile%units( )//"' != '"//                 &
                         this%profiles_( i_profile )%val_%units( ) )
        l_found = .true.
        exit
      end if
    end do

    if( present( found ) ) then
      found = l_found
      if( .not. found ) return
    end if

    call assert_msg( 348928409, l_found,                                      &
                     "Cannot find profile '"//profile%handle_//"'" )

    select type( w_profile => this%profiles_( i_profile )%val_ )
    class is( profile_from_host_t )
      updater = profile_updater_t( w_profile )
    class default
      call die_msg( 166789652, "Cannot update profile '"//w_profile%handle_// &
                               "' from a host application." )
    end select

  end function get_updater

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function pack_size( this, comm )
    ! returns the number of bytes required to pack the warehouse onto a buffer

    use musica_assert,                 only : assert
    use musica_mpi,                    only : musica_mpi_pack_size
    use musica_string,                 only : string_t
    use tuvx_profile_factory,          only : profile_type_name

    class(profile_warehouse_t), intent(in) :: this ! warehouse to be packed
    integer,                    intent(in) :: comm ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: i_profile
    type(string_t) :: type_name

    call assert( 821308000, allocated( this%profiles_ ) )
    pack_size = musica_mpi_pack_size( size( this%profiles_ ), comm )
    do i_profile = 1, size( this%profiles_ )
    associate( profile => this%profiles_( i_profile )%val_ )
      type_name = profile_type_name( profile )
      pack_size = pack_size                                                   &
                  + type_name%pack_size( comm )                               &
                  + profile%pack_size(   comm )
    end associate
    end do
#else
    pack_size = 0
#endif

  end function pack_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine mpi_pack( this, buffer, position, comm )
    ! packs the warehouse onto a character buffer

    use musica_assert,                 only : assert
    use musica_mpi,                    only : musica_mpi_pack
    use musica_string,                 only : string_t
    use tuvx_profile_factory,          only : profile_type_name

    class(profile_warehouse_t), intent(in)    :: this      ! warehouse to be packed
    character,                  intent(inout) :: buffer(:) ! memory buffer
    integer,                    intent(inout) :: position  ! current buffer position
    integer,                    intent(in)    :: comm      ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: prev_pos, i_profile
    type(string_t) :: type_name

    prev_pos = position
    call assert( 333182583, allocated( this%profiles_ ) )
    call musica_mpi_pack( buffer, position, size( this%profiles_ ), comm )
    do i_profile = 1, size( this%profiles_ )
    associate( profile => this%profiles_( i_profile )%val_ )
      type_name = profile_type_name( profile )
      call type_name%mpi_pack( buffer, position, comm )
      call profile%mpi_pack(   buffer, position, comm )
    end associate
    end do
    call assert( 777643951, position - prev_pos <= this%pack_size( comm ) )
#endif

  end subroutine mpi_pack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine mpi_unpack( this, buffer, position, comm )
    ! unpacks a warehouse from a character buffer

    use musica_assert,                 only : assert
    use musica_mpi,                    only : musica_mpi_unpack
    use musica_string,                 only : string_t
    use tuvx_profile_factory,          only : profile_allocate

    class(profile_warehouse_t), intent(out)   :: this      ! warehouse to be unpacked
    character,                  intent(inout) :: buffer(:) ! memory buffer
    integer,                    intent(inout) :: position  ! current buffer position
    integer,                    intent(in)    :: comm      ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: prev_pos, i_profile, n_profiles
    type(string_t) :: type_name

    prev_pos = position
    call musica_mpi_unpack( buffer, position, n_profiles, comm )
    if( allocated( this%profiles_ ) ) deallocate( this%profiles_ )
    allocate( this%profiles_( n_profiles ) )
    do i_profile = 1, n_profiles
    associate( profile => this%profiles_( i_profile ) )
      call type_name%mpi_unpack( buffer, position, comm )
      profile%val_ => profile_allocate( type_name )
      call profile%val_%mpi_unpack( buffer, position, comm )
    end associate
    end do
    call assert( 294782841, position - prev_pos <= this%pack_size( comm ) )
#endif

  end subroutine mpi_unpack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine finalize( this )
    ! Finalize profile warehouse

    use musica_constants, only : ik => musica_ik

    type(profile_warehouse_t), intent(inout) :: this

    ! Local variables
    integer(kind=ik) :: ndx
    character(len=*), parameter :: Iam = 'profile warehouse finalize: '

    if( allocated( this%profiles_ ) ) then
      do ndx = 1, size( this%profiles_ )
        if( associated( this%profiles_( ndx )%val_ ) ) then
          deallocate( this%profiles_( ndx )%val_ )
        end if
      end do
      deallocate( this%profiles_ )
    endif

  end subroutine finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function ptr_pack_size( this, comm ) result( pack_size )
    ! Returns the number of bytes required to pack the pointer onto a buffer

    use musica_mpi,                    only : musica_mpi_pack_size

    class(profile_warehouse_ptr), intent(in) :: this ! This profile pointer
    integer,                      intent(in) :: comm ! MPI communicator

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

    class(profile_warehouse_ptr), intent(in)    :: this      ! pointer to be packed
    character,                    intent(inout) :: buffer(:) ! memory buffer
    integer,                      intent(inout) :: position  ! current buffer position
    integer,                      intent(in)    :: comm      ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: prev_pos

    prev_pos = position
    call musica_mpi_pack( buffer, position, this%index_, comm )
    call assert( 153924730, position - prev_pos <= this%pack_size( comm ) )
#endif

  end subroutine ptr_mpi_pack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ptr_mpi_unpack( this, buffer, position, comm )
    ! Unpacks a pointer from a character buffer

    use musica_assert,                 only : assert
    use musica_mpi,                    only : musica_mpi_unpack

    class(profile_warehouse_ptr), intent(out)   :: this      ! pointer to be packed
    character,                    intent(inout) :: buffer(:) ! memory buffer
    integer,                      intent(inout) :: position  ! current buffer position
    integer,                      intent(in)    :: comm      ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: prev_pos

    prev_pos = position
    call musica_mpi_unpack( buffer, position, this%index_, comm )
    call assert( 266243075, position - prev_pos <= this%pack_size( comm ) )
#endif

  end subroutine ptr_mpi_unpack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_profile_warehouse
