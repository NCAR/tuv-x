! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0

module tuvx_cross_section_warehouse
! A type to store all cross sections needed for a particular run
!

  use musica_string,                   only : string_t
  use tuvx_cross_section,              only : cross_section_ptr

  implicit none

  private
  public :: cross_section_warehouse_t, cross_section_warehouse_ptr

  !> Radiative tranfser cross section type
  type cross_section_warehouse_t
    private
    type(cross_section_ptr), allocatable :: cross_sections_(:) ! A:f:type:`~tuvx_cross_section/cross_section_ptr`
    type(string_t), allocatable          :: handles_(:) ! cross section "handle"
  contains
    !> Get a copy of a specific cross section
    procedure, private :: get_copy_char, get_copy_string
    procedure, private :: get_copy_ptr
    generic :: get => get_copy_char, get_copy_string, get_copy_ptr
    !> Returns a pointer to a cross section in the warehouse
    procedure, private :: get_ptr_char, get_ptr_string
    generic :: get_ptr => get_ptr_char, get_ptr_string
    !> Returns the number of bytes required to pack the warehouse onto a buffer
    procedure :: pack_size
    !> Packs the warehouse onto a character buffer
    procedure :: mpi_pack
    !> Unpacks the warehouse from a character buffer into the object
    procedure :: mpi_unpack
    !> Finalize the object
    final :: finalize
  end type cross_section_warehouse_t

  !> cross_section_warehouse_t constructor
  interface cross_section_warehouse_t
    module procedure :: constructor
  end interface

  !> Pointer to a cross section in the warehouse
  type :: cross_section_warehouse_ptr
    private
    integer :: index_ = 0
  contains
    !> Returns the number of bytes required to pack the pointer onto a buffer
    procedure :: pack_size => ptr_pack_size
    !> Packs the pointer onto a character buffer
    procedure :: mpi_pack => ptr_mpi_pack
    !> Unpacks a pointer from a character buffer into the object
    procedure :: mpi_unpack => ptr_mpi_unpack
  end type cross_section_warehouse_ptr

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function constructor( config, grid_warehouse, profile_warehouse )           &
      result( new_obj )
    ! Constructor of cross_section_warehouse_t objects

    use musica_config,                 only : config_t
    use musica_iterator,               only : iterator_t
    use tuvx_cross_section_factory,    only : cross_section_builder
    use tuvx_grid_warehouse,           only : grid_warehouse_t
    use tuvx_profile_warehouse,        only : Profile_warehouse_t

    class(cross_section_warehouse_t), pointer :: new_obj ! New radiative :f:type:`~tuvx_cross_section_warehouse/cross_section_warehouse_t` object
    type(config_t),             intent(inout) :: config ! Cross section configuration data
    type(grid_warehouse_t),     intent(inout) :: grid_warehouse ! A :f:type:`~tuvx_grid_warehouse/grid_warehouse_t`
    type(Profile_warehouse_t),  intent(inout) :: profile_warehouse ! A :f:type:`~tuvx_profile_warehouse/profile_warehouse_t`

    ! local variables
    character(len=*), parameter :: Iam = "Cross section constructor"
    type(config_t) :: cross_section_config
    class(iterator_t), pointer :: iter
    type(cross_section_ptr), allocatable :: temp_cross_sections(:)
    type(cross_section_ptr) :: cross_section
    type(string_t)              :: cross_section_name
    type(string_t), allocatable :: temp_names(:)

    allocate( new_obj )
    allocate( string_t :: new_obj%handles_(0) )
    allocate( new_obj%cross_sections_(0) )

    ! iterate over cross sections
    iter => config%get_iterator( )
    do while( iter%next( ) )
      call config%get( iter, cross_section_config, Iam )

      ! save the cross section name for lookups
      call cross_section_config%get( "name", cross_section_name, Iam )
      temp_names = new_obj%handles_
      deallocate( new_obj%handles_ )
      allocate( new_obj%handles_( size( temp_names ) + 1 ) )
      new_obj%handles_( 1 : size( temp_names ) ) = temp_names(:)
      new_obj%handles_( size( temp_names ) + 1 ) = cross_section_name
      deallocate( temp_names )

      ! build and store cross section object
      cross_section%val_ => cross_section_builder( cross_section_config,      &
                                           grid_warehouse, profile_warehouse )
      temp_cross_sections = new_obj%cross_sections_
      deallocate( new_obj%cross_sections_ )
      allocate( new_obj%cross_sections_( size( temp_cross_sections ) + 1 ) )
      new_obj%cross_sections_( 1 : size( temp_cross_sections ) ) =            &
          temp_cross_sections(:)
      new_obj%cross_sections_( size( temp_cross_sections ) + 1 ) =            &
          cross_section
      deallocate( temp_cross_sections )
      nullify( cross_section%val_ )
    end do
    deallocate( iter )

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function get_copy_char( this, name ) result( cross_section )
    ! Get a copy of a specific cross section object

    use tuvx_cross_section,            only : cross_section_t

    class(cross_section_warehouse_t), intent(in) :: this          ! cross section warehouse
    character(len=*),                 intent(in) :: name          ! name of the cross section to copy
    class(cross_section_t), pointer              :: cross_section ! copy of the cross section

    type(cross_section_warehouse_ptr) :: ptr

    ptr = this%get_ptr( name )
    allocate( cross_section, source = this%cross_sections_( ptr%index_ )%val_ )

  end function get_copy_char

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function get_copy_string( this, name ) result( cross_section )
    ! Get a copy of a specific radiative transfer cross section object

    use musica_assert,                 only : die_msg
    use musica_string,                 only : string_t
    use tuvx_cross_section,            only : cross_section_t

    class(cross_section_t),           pointer       :: cross_section ! Pointer to a copy of the requested :f:type:`~tuvx_cross_section/cross_section_t`
    class(cross_section_warehouse_t), intent(inout) :: this ! A :f:type:`~tuvx_cross_section_warehouse/cross_section_warehouse_t`
    type(string_t),                   intent(in)    :: name ! Name of the cross section to find

    cross_section => this%get_copy_char( name%to_char( ) )

  end function get_copy_string

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function get_copy_ptr( this, ptr ) result( cross_section )
    ! Returns a copy of a cross section from the warehouse

    use musica_assert,                 only : assert_msg
    use tuvx_cross_section,            only : cross_section_t

    class(cross_section_warehouse_t),  intent(in) :: this ! cross section warehouse
    type(cross_section_warehouse_ptr), intent(in) :: ptr ! cross section pointer
    class(cross_section_t), pointer               :: cross_section ! copy of the cross section

    call assert_msg( 829999477, ptr%index_ > 0,                               &
                     "Invalid cross section pointer" )
    allocate( cross_section, source = this%cross_sections_( ptr%index_ )%val_ )

  end function get_copy_ptr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function get_ptr_char( this, name ) result( ptr )
    ! Returns a pointer to a cross section in the warehouse

    use musica_assert,                 only : die_msg

    class(cross_section_warehouse_t), intent(in) :: this ! cross section warehouse
    character(len=*),                 intent(in) :: name ! name of the cross section to find
    type(cross_section_warehouse_ptr)            :: ptr  ! pointer to the cross section

    integer :: ndx
    logical :: found

    found = .false.
    do ndx = 1, size( this%handles_ )
      if( name .eq. this%handles_( ndx ) ) then
        found = .true.
        exit
      endif
    end do

    if( .not. found ) then
      call die_msg( 980636382, "Invalid cross_section_name: '"//name//"'" )
    endif
    ptr%index_ = ndx

  end function get_ptr_char

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function get_ptr_string( this, name ) result( ptr )
    ! Returns a pointer to a cross section in the warehouse

    use musica_string,                 only : string_t

    class(cross_section_warehouse_t), intent(in) :: this ! cross section warehouse
    type(string_t),                   intent(in) :: name ! name of the cross section to find
    type(cross_section_warehouse_ptr)            :: ptr  ! pointer to the cross section

    ptr = this%get_ptr_char( name%to_char( ) )

  end function get_ptr_string

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function pack_size( this, comm )
    ! Returns the number of bytes required to pack the warehouse onto a buffer

    use musica_assert,                 only : assert
    use musica_mpi,                    only : musica_mpi_pack_size
    use tuvx_cross_section_factory,    only : cross_section_type_name

    class(cross_section_warehouse_t), intent(in) :: this ! warehouse to be packed
    integer,                          intent(in) :: comm ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: i_cross_section
    type(string_t) :: type_name

    call assert( 161062144, allocated( this%cross_sections_ ) )
    call assert( 385698834, allocated( this%handles_ ) )
    call assert( 259945397,                                                   &
                 size( this%cross_sections_ ) .eq. size( this%handles_ ) )
    pack_size = musica_mpi_pack_size( size( this%cross_sections_ ), comm )
    do i_cross_section = 1, size( this%cross_sections_ )
    associate( cross_section => this%cross_sections_( i_cross_section )%val_ )
      type_name = cross_section_type_name( cross_section )
      pack_size = pack_size +                                                 &
                  type_name%pack_size( comm ) +                               &
                  cross_section%pack_size( comm ) +                           &
                  this%handles_( i_cross_section )%pack_size( comm )
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
    use tuvx_cross_section_factory,    only : cross_section_type_name

    class(cross_section_warehouse_t), intent(in) :: this ! warehouse to be packed
    character,                        intent(inout) :: buffer(:) ! memory buffer
    integer,                          intent(inout) :: position  ! current buffer position
    integer,                          intent(in)    :: comm      ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: prev_pos, i_cross_section
    type(string_t) :: type_name

    prev_pos = position
    call assert( 442823572, allocated( this%cross_sections_ ) )
    call assert( 386890547, allocated( this%handles_        ) )
    call assert( 264043588,                                                   &
                 size( this%cross_sections_ ) .eq. size( this%handles_ ) )
    call musica_mpi_pack( buffer, position, size( this%cross_sections_ ),     &
                          comm )
    do i_cross_section = 1, size( this%cross_sections_ )
    associate( cross_section => this%cross_sections_( i_cross_section )%val_ )
      type_name = cross_section_type_name( cross_section )
      call type_name%mpi_pack( buffer, position, comm )
      call cross_section%mpi_pack( buffer, position, comm )
      call this%handles_( i_cross_section )%mpi_pack( buffer, position, comm )
    end associate
    end do
    call assert( 359568068, position - prev_pos <= this%pack_size( comm ) )
#endif

  end subroutine mpi_pack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine mpi_unpack( this, buffer, position, comm )
    ! Unpacks the warehouse from a character buffer

    use musica_assert,                 only : assert
    use musica_mpi,                    only : musica_mpi_unpack
    use tuvx_cross_section_factory,    only : cross_section_allocate

    class(cross_section_warehouse_t), intent(out)   :: this      ! warehouse to be unpacked
    character,                        intent(inout) :: buffer(:) ! memory buffer
    integer,                          intent(inout) :: position  ! current buffer position
    integer,                          intent(in)    :: comm      ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: prev_pos, i_cross_section, n_cross_sections
    type(string_t) :: type_name

    prev_pos = position
    call musica_mpi_unpack( buffer, position, n_cross_sections, comm )
    if( allocated( this%cross_sections_ ) ) deallocate( this%cross_sections_ )
    if( allocated( this%handles_        ) ) deallocate( this%handles_        )
    allocate( this%cross_sections_( n_cross_sections ) )
    allocate( this%handles_(        n_cross_sections ) )
    do i_cross_section = 1, n_cross_sections
    associate( cross_section => this%cross_sections_( i_cross_section ) )
      call type_name%mpi_unpack( buffer, position, comm )
      cross_section%val_ => cross_section_allocate( type_name )
      call cross_section%val_%mpi_unpack( buffer, position, comm )
      call this%handles_( i_cross_section )%mpi_unpack( buffer, position,     &
                                                        comm )
    end associate
    end do
    call assert( 659576600, position - prev_pos <= this%pack_size( comm ) )
#endif

  end subroutine mpi_unpack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine finalize( this )
    ! Finalize the cross section warehouse

    type(cross_section_warehouse_t), intent(inout) :: this ! This :f:type:`~tuvx_cross_section_warehouse/cross_section_warehouse_t`

    integer :: ndx

    if( allocated( this%cross_sections_ ) ) then
      do ndx = 1,size( this%cross_sections_ )
        if( associated( this%cross_sections_( ndx )%val_ ) ) then
          deallocate( this%cross_sections_( ndx )%val_ )
        endif
      enddo
      deallocate( this%cross_sections_ )
    end if

    if( allocated( this%handles_ ) ) then
      deallocate( this%handles_ )
    end if

  end subroutine finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function ptr_pack_size( this, comm ) result( pack_size )
    ! Returns the number of bytes required to pack the pointer onto a buffer

    use musica_mpi,                    only : musica_mpi_pack_size

    class(cross_section_warehouse_ptr), intent(in) :: this ! This cross section pointer
    integer,                            intent(in) :: comm ! MPI communicator

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

    class(cross_section_warehouse_ptr), intent(in)    :: this      ! pointer to be packed
    character,                          intent(inout) :: buffer(:) ! memory buffer
    integer,                            intent(inout) :: position  ! current buffer position
    integer,                            intent(in)    :: comm      ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: prev_pos

    prev_pos = position
    call musica_mpi_pack( buffer, position, this%index_, comm )
    call assert( 344231889, position - prev_pos <= this%pack_size( comm ) )
#endif

  end subroutine ptr_mpi_pack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ptr_mpi_unpack( this, buffer, position, comm )
    ! Unpacks a pointer from a character buffer

    use musica_assert,                 only : assert
    use musica_mpi,                    only : musica_mpi_unpack

    class(cross_section_warehouse_ptr), intent(out)   :: this      ! pointer to be packed
    character,                          intent(inout) :: buffer(:) ! memory buffer
    integer,                            intent(inout) :: position  ! current buffer position
    integer,                            intent(in)    :: comm      ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: prev_pos

    prev_pos = position
    call musica_mpi_unpack( buffer, position, this%index_, comm )
    call assert( 286393330, position - prev_pos <= this%pack_size( comm ) )
#endif

  end subroutine ptr_mpi_unpack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_cross_section_warehouse
