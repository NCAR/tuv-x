! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> The musica_map module

!> Utility for mapping among arrays
module musica_map

  use musica_constants,                only : dk => musica_dk

  implicit none
  private

  public :: map_t

  !> Matched pair
  type :: pair_t
    private
    !> Index in source array
    integer :: from_index_
    !> Index in destination array
    integer :: to_index_
    !> Scaling factor applied to source data
    real(kind=dk) :: scale_factor_ = 1.0
  contains
    !> Returns the size of a binary buffer required to pack the pair
    procedure :: pack_size => pair_pack_size
    !> Packs the pair onto a characcter buffer
    procedure :: mpi_pack => pair_mpi_pack
    !> Unpacks a pair from a character buffer
    procedure :: mpi_unpack => pair_mpi_unpack
  end type pair_t

  !> Constructor of pair_t objects
  interface pair_t
    module procedure :: pair_constructor
  end interface pair_t

  !> Map between arrays
  !!
  !! Maps can be used to transfer data from a source to a destination array
  !! with optional scaling.
  !!
  !! The mapped elements are identified by name according to the passed
  !! configuration. The configuration format for a map is:
  !! \code{json}
  !!   {
  !!     "match full source": false,
  !!     "match full destination": false,
  !!     "sum multiple matches": true,
  !!     "default matching": "backup",
  !!     "pairs": [
  !!       {
  !!         "from": "foo",
  !!         "to": "bar"
  !!       },
  !!       {
  !!         "from": "baz",
  !!         "to": "quz",
  !!         "scale by": 1.2
  !!       }
  !!     ]
  !!   }
  !! \endcode
  !!
  !! The "match full source" and "match full destination" terms are optional
  !! and default to \c true.
  !! When these are \c true unmatched source/destination array elements will
  !! trigger an error.
  !! If unmatched destination elements are allowed, they will be set to
  !! zero when the map is applied to transfer data.
  !! The "sum multiple matches" term is optional and defaults to \c false.
  !! When this is \c true, multiple matches to a single destination array
  !! element will be summed when the map is applied to transfer data.
  !! When this is \c false, the second match to a destination array
  !! element will trigger an error.
  !! The "default matching" term is optional and indicates how matching
  !! names that appear in both source and destination label arrays
  !! should be treated.
  !! The three options for default matching are "always", "backup",
  !! and "never"; the default option is "never".
  !! Default matching "always" indicates that every time a name appears
  !! in both the source and destination label arrays, a set of paired
  !! elements should be created in the map with a scaling factor of 1.0.
  !! Default matching "backup" indicates that such a pair is only
  !! created when no explicit entries for the destination element
  !! exist in the configuration.
  !! Default mapping "never" means that no such pairs are created.
  !! If the default mapping is set to "always" or "backup", the
  !! "match full destination" term must be \c true.
  !!
  !! The "pairs" term is required and is an array that
  !! describes each matched pair of elements. The matched pair terms
  !! must include "from" and "to" terms.
  !! The "scale by" term is optional and defaults to 1.0.
  !! This scaling factor will be applied to the source array element
  !! before additon to the destination array element.
  !!
  !! The \c map_t constructor accepts an array of source element labels
  !! and an array of destination element labels that are used to
  !! identify the mapped array indices.
  !!
  type :: map_t
    private
    !> Mapped pairs of array elements
    type(pair_t), allocatable :: pairs_(:)
    !> Source array size
    integer :: from_size_
    !> Destination array size
    integer :: to_size_
  contains
    !> Transfers data from source to destination arrays
    procedure :: apply
    !> Returns the size of a character buffer required to pack the map
    procedure :: pack_size
    !> Packs the map onto a character buffer
    procedure :: mpi_pack
    !> Unpacks the map from a character buffer
    procedure :: mpi_unpack
    !> Prints the map
    procedure :: print => print_map
    !> Adds default matches by name to the map
    procedure, private :: add_default_matches
    !> Validates the matches based on user-selected options
    procedure, private :: validate
  end type map_t

  !> Constructor of map_t objects
  interface map_t
    module procedure :: constructor
  end interface map_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Constructs a map_t object
  type(map_t) function constructor( config, from_labels, to_labels )          &
      result( this )

    use musica_assert,                 only : assert_msg
    use musica_config,                 only : config_t
    use musica_iterator,               only : iterator_t
    use musica_string,                 only : string_t

    !> Map configuration
    type(config_t), intent(inout) :: config
    !> Source array element labels
    type(string_t), intent(in) :: from_labels(:)
    !> Destination array element labels
    type(string_t), intent(in) :: to_labels(:)

    character(len=*), parameter :: my_name = "Map constructor"
    type(config_t) :: pairs, pair
    class(iterator_t), pointer :: iter
    integer :: i_pair
    integer, allocatable :: source_match(:), dest_match(:)
    type(string_t) :: default_matching
    type(string_t) :: required_keys(1), optional_keys(4)

    required_keys(1) = "pairs"
    optional_keys(1) = "match full source"
    optional_keys(2) = "match full destination"
    optional_keys(3) = "sum multiple matches"
    optional_keys(4) = "default matching"

    call assert_msg( 170733942,                                               &
                     config%validate( required_keys, optional_keys ),         &
                     "Bad configuration format for map." )
    call config%get( "default matching", default_matching, my_name,           &
                     default = "never" )
    call config%get( "pairs", pairs, my_name )

    this%from_size_ = size( from_labels )
    this%to_size_   = size( to_labels   )

    ! Get all matched pairs
    allocate( this%pairs_( pairs%number_of_children( ) ) )
    iter => pairs%get_iterator( )
    i_pair = 0
    do while( iter%next( ) )
      call pairs%get( iter, pair, my_name )
      i_pair = i_pair + 1
      this%pairs_( i_pair ) = pair_t( pair, from_labels, to_labels )
    end do
    deallocate( iter )

    if( default_matching == "always" ) then
      call this%add_default_matches( from_labels, to_labels, always = .true. )
    else if( default_matching == "backup" ) then
      call this%add_default_matches( from_labels, to_labels, always = .false. )
    else
      call assert_msg( 135980113, default_matching == "never",                &
                       "Invalid default matching option for map creation: '"//&
                       default_matching//"'" )
    end if

    call this%validate( config, from_labels, to_labels )

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Transfers data from source to destination array based on map
  subroutine apply( this, from, to )

    use musica_assert,                 only : assert_msg

    !> Map
    class(map_t),  intent(in)  :: this
    !> Source array
    real(kind=dk), intent(in)  :: from(:)
    !> Destination array
    real(kind=dk), intent(out) :: to(:)

    integer :: i_elem

    call assert_msg( 764798475, size( from ) .eq. this%from_size_,            &
                     "Wrong size for mapped source array." )
    call assert_msg( 133386338, size( to ) .eq. this%to_size_,                &
                     "Wrong size for mapped destination array." )
    to(:) = 0.0_dk
    do i_elem = 1, size( this%pairs_ )
    associate( pair => this%pairs_( i_elem ) )
      to( pair%to_index_ ) = to( pair%to_index_ ) +                           &
                             from( pair%from_index_ ) *                       &
                             pair%scale_factor_
    end associate
    end do

  end subroutine apply

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the size of a binary buffer required to pack the map
  integer function pack_size( this, comm )

    use musica_mpi

    !> Map to pack
    class(map_t),   intent(in) :: this
    !> MPI communicator
    integer,        intent(in) :: comm

#ifdef MUSICA_USE_MPI
    integer :: i_pair

    pack_size = musica_mpi_pack_size( allocated( this%pairs_ ), comm )
    if( allocated( this%pairs_ ) ) then
      pack_size = pack_size +                                                 &
                  musica_mpi_pack_size( size( this%pairs_ ), comm )
      do i_pair = 1, size( this%pairs_ )
        pack_size = pack_size + this%pairs_( i_pair )%pack_size( comm )
      end do
    end if
    pack_size = pack_size + musica_mpi_pack_size( this%from_size_, comm ) +   &
                            musica_mpi_pack_size( this%to_size_,   comm )
#else
    pack_size = 0
#endif

  end function pack_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Packs the map onto a character buffer
  subroutine mpi_pack( this, buffer, position, comm )

    use musica_assert,                 only : assert
    use musica_mpi

    !> Map to pack
    class(map_t),   intent(in)    :: this
    !> Memory buffer
    character,      intent(inout) :: buffer(:)
    !> Current buffer position
    integer,        intent(inout) :: position
    !> MPI communicator
    integer,        intent(in)    :: comm

#ifdef MUSICA_USE_MPI
    integer :: i_pair, prev_position

    prev_position = position
    call musica_mpi_pack( buffer, position, allocated( this%pairs_ ), comm )
    if( allocated( this%pairs_ ) ) then
      call musica_mpi_pack( buffer, position, size( this%pairs_ ), comm )
      do i_pair = 1, size( this%pairs_ )
        call this%pairs_( i_pair )%mpi_pack( buffer, position, comm )
      end do
    end if
    call musica_mpi_pack( buffer, position, this%from_size_, comm )
    call musica_mpi_pack( buffer, position, this%to_size_,   comm )
    call assert( 419959778,                                                   &
                 position - prev_position <= this%pack_size( comm ) )
#endif

  end subroutine mpi_pack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Unpacks a map from a character buffer
  subroutine mpi_unpack( this, buffer, position, comm )

    use musica_assert,                 only : assert
    use musica_mpi

    !> Map to unpack
    class(map_t),   intent(out)   :: this
    !> Memory buffer
    character,      intent(inout) :: buffer(:)
    !> Current buffer position
    integer,        intent(inout) :: position
    !> MPI communicator
    integer,        intent(in)    :: comm

#ifdef MUSICA_USE_MPI
    logical :: alloced
    integer :: i_pair, n_pairs, prev_position

    prev_position = position
    call musica_mpi_unpack( buffer, position, alloced, comm )
    if( alloced ) then
      call musica_mpi_unpack( buffer, position, n_pairs, comm )
      allocate( this%pairs_( n_pairs ) )
      do i_pair = 1, size( this%pairs_ )
        call this%pairs_( i_pair )%mpi_unpack( buffer, position, comm )
      end do
    end if
    call musica_mpi_unpack( buffer, position, this%from_size_, comm )
    call musica_mpi_unpack( buffer, position, this%to_size_,   comm )
    call assert( 576681590,                                                   &
                 position - prev_position <= this%pack_size( comm ) )
#endif

  end subroutine mpi_unpack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Prints the map details to a specified output unit
  subroutine print_map( this, from_labels, to_labels, out_unit )

    use musica_assert,                 only : assert_msg
    use musica_string,                 only : string_t, output_table

    !> Map
    class(map_t),   intent(in) :: this
    !> Source array element labels
    type(string_t), intent(in)    :: from_labels(:)
    !> Destination array element labels
    type(string_t), intent(in)    :: to_labels(:)
    !> Output unit
    integer,        intent(in) :: out_unit

    type(string_t) :: header(3)
    type(string_t), allocatable :: table(:,:)
    integer :: i_pair

    call assert_msg( 727878410, size( from_labels ) .eq. this%from_size_,     &
                     "Wrong size for map source label array." )
    call assert_msg( 161474673, size( to_labels ) .eq. this%to_size_,         &
                     "Wrong size for map destination label array." )
    if( .not. allocated( this%pairs_ ) ) then
      write(out_unit,*) "Map not initialized"
      return
    end if
    header(1) = "from"
    header(2) = "to"
    header(3) = "scaling factor"
    allocate( table( 3, size( this%pairs_ ) ) )
    do i_pair = 1, size( this%pairs_ )
    associate( pair => this%pairs_( i_pair ) )
      table( 1, i_pair ) = from_labels( pair%from_index_ )
      table( 2, i_pair ) = to_labels(   pair%to_index_   )
      table( 3, i_pair ) = pair%scale_factor_
    end associate
    end do
    call output_table( header, table, out_unit )

  end subroutine print_map

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Adds default matches by name to the map
  !!
  !! If the \c always option is set to \c false, only unmatched source
  !! elements are included in the default matching
  subroutine add_default_matches( this, from_labels, to_labels, always )

    use musica_array,                  only : find_string_in_array
    use musica_string,                 only : string_t

    !> Map
    class(map_t),   intent(inout) :: this
    !> Source array element labels
    type(string_t), intent(in)    :: from_labels(:)
    !> Destination array element labels
    type(string_t), intent(in)    :: to_labels(:)
    !> Flag indicating whether to always add default matches, or only do so
    !! for unmatched source elements
    logical,        intent(in)    :: always

    integer :: matches( size( to_labels ) )
    integer :: i_to, i_from, i_pair
    type(pair_t) :: pair

    matches(:) = 0
    if( .not. always ) then
      do i_pair = 1, size( this%pairs_ )
        i_to = this%pairs_( i_pair )%to_index_
        matches( i_to ) = matches( i_to ) + 1
      end do
    end if
    do i_to = 1, size( to_labels )
      if( matches( i_to ) > 0 ) cycle
      if( find_string_in_array( from_labels, to_labels( i_to ), i_from,       &
                                case_sensitive = .true. ) ) then
        pair%to_index_ = i_to
        pair%from_index_ = i_from
        pair%scale_factor_ = 1.0_dk
        this%pairs_ = [ this%pairs_, pair ]
      end if
    end do

  end subroutine add_default_matches

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Validates the map based on user-selected options
  subroutine validate( this, config, from_labels, to_labels )

    use musica_assert,                 only : assert_msg
    use musica_config,                 only : config_t
    use musica_string,                 only : string_t

    !> Map
    class(map_t),   intent(in)    :: this
    !> Map configuration
    type(config_t), intent(inout) :: config
    !> Source array element labels
    type(string_t), intent(in)    :: from_labels(:)
    !> Destination array element labels
    type(string_t), intent(in)    :: to_labels(:)

    character(len=*), parameter :: my_name = "Map validation"
    integer, allocatable :: match(:)
    type(string_t) :: default_matching
    integer :: i_pair, i_elem
    logical :: match_source
    logical :: match_dest
    logical :: allow_sum

    call config%get( "match full source", match_source, my_name,              &
                     default = .true. )
    call config%get( "match full destination", match_dest, my_name,           &
                     default = .true. )
    call config%get( "sum multiple matches", allow_sum, my_name,              &
                     default = .false. )
    call config%get( "default matching", default_matching, my_name,           &
                     default = "never" )

    call assert_msg( 548594113, match_dest .or. default_matching == "never",  &
                     "Default matching is only possible when matching the "// &
                     "full destination array for maps." )

    if( match_source ) then
      allocate( match( this%from_size_ ) )
      match(:) = 0
      do i_pair = 1, size( this%pairs_ )
      associate( match_elem => match( this%pairs_( i_pair )%from_index_ ) )
        match_elem = match_elem + 1
      end associate
      end do
      do i_elem = 1, size( match )
        call assert_msg( 956987954, match( i_elem ) > 0,                      &
                         "Unmatched element '"//from_labels( i_elem )//       &
                         "' in source array of map." )
      end do
      deallocate( match )
    end if

    if( match_dest .or. .not. allow_sum ) then
      allocate( match( this%to_size_ ) )
      match(:) = 0
      do i_pair = 1, size( this%pairs_ )
      associate( match_elem => match( this%pairs_( i_pair )%to_index_ ) )
        match_elem = match_elem + 1
      end associate
      end do
      do i_elem = 1, size( match )
        call assert_msg( 200274675,                                           &
                         match( i_elem ) > 0 .or. .not. match_dest,           &
                         "Unmatched element '"//to_labels( i_elem )//         &
                         "' in destination array of map." )
        call assert_msg( 240867074,                                           &
                         match( i_elem ) < 2 .or. allow_sum,                  &
                         "Multiple matches found for element '"//             &
                         to_labels( i_elem )//                                &
                         "' in destination array of map." )
      end do
      deallocate( match )
    end if

  end subroutine validate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Constructor of pair_t objects
  type(pair_t) function pair_constructor( config, from_labels, to_labels )    &
      result( this )

    use musica_array,                  only : find_string_in_array
    use musica_assert,                 only : assert_msg
    use musica_config,                 only : config_t
    use musica_string,                 only : string_t

    !> Matched pair configuration
    type(config_t), intent(inout) :: config
    !> Source array element labels
    type(string_t), intent(in) :: from_labels(:)
    !> Destination array element labels
    type(string_t), intent(in) :: to_labels(:)

    character(len=*), parameter :: my_name = "Map pair constructor"
    type(string_t) :: label
    type(string_t) :: required_keys(2), optional_keys(1)

    required_keys(1) = "from"
    required_keys(2) = "to"
    optional_keys(1) = "scale by"

    call assert_msg( 309595761,                                               &
                     config%validate( required_keys, optional_keys ),         &
                     "Bad configuration format for map pair." )

    call config%get( "from", label, my_name )
    call assert_msg( 122570601,                                               &
                     find_string_in_array( from_labels, label,                &
                         this%from_index_, case_sensitive = .true. ),         &
                     "Cannot find source label '"//label//"' building map." )
    call config%get( "to", label, my_name )
    call assert_msg( 740547646,                                               &
                     find_string_in_array( to_labels, label,                  &
                         this%to_index_, case_sensitive = .true. ),           &
                     "Cannot find destination label '"//label//               &
                     "' building map." )
    call config%get( "scale by", this%scale_factor_, my_name,                 &
                     default = 1.0_dk )

  end function pair_constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the size of a binary buffer required to pack the pair
  integer function pair_pack_size( this, comm ) result( pack_size )

    use musica_mpi

    !> Pair to pack
    class(pair_t),   intent(in) :: this
    !> MPI communicator
    integer,         intent(in) :: comm

#ifdef MUSICA_USE_MPI
    pack_size = musica_mpi_pack_size( this%from_index_,   comm ) +            &
                musica_mpi_pack_size( this%to_index_,     comm ) +            &
                musica_mpi_pack_size( this%scale_factor_, comm )
#else
    pack_size = 0
#endif

  end function pair_pack_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Packs the pair onto a character buffer
  subroutine pair_mpi_pack( this, buffer, position, comm )

    use musica_assert,                 only : assert
    use musica_mpi

    !> Pair to pack
    class(pair_t),  intent(in)    :: this
    !> Memory buffer
    character,      intent(inout) :: buffer(:)
    !> Current buffer position
    integer,        intent(inout) :: position
    !> MPI communicator
    integer,        intent(in)    :: comm

#ifdef MUSICA_USE_MPI
    integer :: prev_position

    prev_position = position
    call musica_mpi_pack( buffer, position, this%from_index_,   comm )
    call musica_mpi_pack( buffer, position, this%to_index_,     comm )
    call musica_mpi_pack( buffer, position, this%scale_factor_, comm )
    call assert( 995726013,                                                   &
                 position - prev_position <= this%pack_size( comm ) )
#endif

  end subroutine pair_mpi_pack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Unpacks a pair from a character buffer
  subroutine pair_mpi_unpack( this, buffer, position, comm )

    use musica_assert,                 only : assert
    use musica_mpi

    !> Pair to unpack
    class(pair_t),  intent(out)   :: this
    !> Memory buffer
    character,      intent(inout) :: buffer(:)
    !> Current buffer position
    integer,        intent(inout) :: position
    !> MPI communicator
    integer,        intent(in)    :: comm

#ifdef MUSICA_USE_MPI
    integer :: prev_position

    prev_position = position
    call musica_mpi_unpack( buffer, position, this%from_index_,   comm )
    call musica_mpi_unpack( buffer, position, this%to_index_,     comm )
    call musica_mpi_unpack( buffer, position, this%scale_factor_, comm )
    call assert( 143488254,                                                   &
                 position - prev_position <= this%pack_size( comm ) )
#endif

  end subroutine pair_mpi_unpack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module musica_map
