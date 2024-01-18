! Copyright (C) 2020-4 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0

module tuvx_temperature_range
! Defines a temperature range for use in temperature-based cross
! section parameterizations

  ! Including musica_config at the module level to avoid an ICE
  ! with Intel 2022.1 compiler
  use musica_config,                   only : config_t
  use musica_constants,                only : dk => musica_dk

  implicit none

  private
  public :: temperature_range_t


  !> Range for temperature-based calculations
  type :: temperature_range_t
    !> Minimum temperature [K] for inclusion in range
    real(kind=dk) :: min_temperature_ = 0.0_dk
    !> Maximum temperature [K] for include in range
    real(kind=dk) :: max_temperature_ = huge(1.0_dk)
    !> Indicates whether to use a fixed temperature for the
    !! parameterization calculation. If FALSE, the actual
    !! temperature is used.
    logical :: is_fixed_ = .false.
    !> Fixed temperature [K] to use in paramterization calculation
    !!
    !! Is only used if is_fixed == TRUE
    real(kind=dk) :: fixed_temperature_ = 0.0_dk
  contains
    !> Returns the number of bytes required to pack the range onto a
    !! character buffer
    procedure :: pack_size => pack_size
    !> Packs the range onto a character buffer
    procedure :: mpi_pack => mpi_pack
    !> Unpacks a range from a character buffer
    procedure :: mpi_unpack => mpi_unpack
  end type temperature_range_t

  !> Constructor for temperature_range_t
  interface temperature_range_t
    module procedure :: constructor
  end interface temperature_range_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function constructor( config ) result( this )
    ! Constructs temperature range objects

    use musica_assert,                 only : assert_msg
    use musica_config,                 only : config_t
    use musica_string,                 only : string_t

    type(temperature_range_t)               :: this
    type(config_t),           intent(inout) :: config

    character(len=*), parameter :: my_name = "temperature range constructor"
    type(string_t) :: required_keys(0), optional_keys(3)
    logical :: found

    optional_keys(1) = "minimum"
    optional_keys(2) = "maximum"
    optional_keys(3) = "fixed value"
    call assert_msg( 355912601,                                               &
                     config%validate( required_keys, optional_keys ),         &
                     "Bad configuration for temperature range" )

    call config%get( "minimum", this%min_temperature_, my_name,               &
                     default = 0.0_dk )
    call config%get( "maximum", this%max_temperature_, my_name,               &
                     default = huge(1.0_dk) )
    call config%get( "fixed value", this%fixed_temperature_, my_name,         &
                     found = found )
    this%is_fixed_ = found

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function pack_size( this, comm )
    ! Returns the size of a character buffer required to pack the range

    use musica_mpi,                    only : musica_mpi_pack_size

    class(temperature_range_t), intent(in) :: this ! temperature range to be packed
    integer,                    intent(in) :: comm ! MPI communicator

#ifdef MUSICA_USE_MPI
    pack_size = musica_mpi_pack_size( this%min_temperature_,   comm ) +       &
                musica_mpi_pack_size( this%max_temperature_,   comm ) +       &
                musica_mpi_pack_size( this%is_fixed_,          comm ) +       &
                musica_mpi_pack_size( this%fixed_temperature_, comm )
#else
    pack_size = 0
#endif

  end function pack_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine mpi_pack( this, buffer, position, comm )
    ! Packs the temperature range onto a character buffer

    use musica_assert,                 only : assert
    use musica_mpi,                    only : musica_mpi_pack

    class(temperature_range_t), intent(in)    :: this      ! temperature range to be packed
    character,                  intent(inout) :: buffer(:) ! memory buffer
    integer,                    intent(inout) :: position  ! current buffer position
    integer,                    intent(in)    :: comm      ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: prev_pos

    prev_pos = position
    call musica_mpi_pack( buffer, position, this%min_temperature_,   comm )
    call musica_mpi_pack( buffer, position, this%max_temperature_,   comm )
    call musica_mpi_pack( buffer, position, this%is_fixed_,          comm )
    call musica_mpi_pack( buffer, position, this%fixed_temperature_, comm )
    call assert( 409699380, position - prev_pos <= this%pack_size( comm ) )
#endif

  end subroutine mpi_pack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine mpi_unpack( this, buffer, position, comm )
    ! Unpacks a temperature range from a character buffer

    use musica_assert,                 only : assert
    use musica_mpi,                    only : musica_mpi_unpack

    class(temperature_range_t), intent(out)   :: this      ! temperature range to be unpacked
    character,                  intent(inout) :: buffer(:) ! memory buffer
    integer,                    intent(inout) :: position  ! current buffer position
    integer,                    intent(in)    :: comm      ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: prev_pos

    prev_pos = position
    call musica_mpi_unpack( buffer, position, this%min_temperature_,   comm )
    call musica_mpi_unpack( buffer, position, this%max_temperature_,   comm )
    call musica_mpi_unpack( buffer, position, this%is_fixed_,          comm )
    call musica_mpi_unpack( buffer, position, this%fixed_temperature_, comm )
    call assert( 164457757, position - prev_pos <= this%pack_size( comm ) )
#endif

  end subroutine mpi_unpack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_temperature_range