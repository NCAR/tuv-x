! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0

module tuvx_radiator
! Represents an atmospheric constituent that affects radiative transfer calculations by absorbing or scattering radiation

  use musica_assert,                 only : assert, assert_msg
  use musica_config,                 only : config_t
  use musica_constants,              only : dk => musica_dk
  use musica_mpi,                    only : musica_mpi_pack, musica_mpi_pack_size, musica_mpi_unpack
  use musica_string,                 only : string_t
  use tuvx_constants,                only : largest, precis
  use tuvx_cross_section,            only : cross_section_t
  use tuvx_cross_section_warehouse,  only : cross_section_warehouse_ptr
  use tuvx_cross_section_warehouse,  only : cross_section_warehouse_t
  use tuvx_diagnostic_util,          only : diagout
  use tuvx_grid,                     only : grid_t
  use tuvx_grid_warehouse,           only : grid_warehouse_ptr
  use tuvx_grid_warehouse,           only : grid_warehouse_t
  use tuvx_profile,                  only : profile_t
  use tuvx_profile_warehouse,        only : profile_warehouse_ptr
  use tuvx_profile_warehouse,        only : profile_warehouse_t

  implicit none

  private
  public :: radiator_t, radiator_ptr, radiator_state_t, base_constructor

  type radiator_state_t
    ! Optical properties for a radiator

    real(kind=dk), allocatable :: layer_OD_(:,:)  ! layer optical depth (vertical layer, wavelength bin)
    real(kind=dk), allocatable :: layer_SSA_(:,:) ! layer single scattering albedo (vertical layer, wavelength bin)
    real(kind=dk), allocatable :: layer_G_(:,:,:) ! layer asymmetry factor (vertical layer, wavelength bin, stream)
  contains
    ! Accumulates a net radiator state for a set of radiators
    procedure :: accumulate
    ! Returns the number of bytes needed to pack the object onto a buffer
    procedure :: pack_size => state_pack_size
    ! Packs the object onto a character buffer
    procedure :: mpi_pack => state_mpi_pack
    ! Unpacks data from a character buffer into the object
    procedure :: mpi_unpack => state_mpi_unpack
    final :: finalize
  end type radiator_state_t

  type radiator_t
    ! Optically active species

    type(string_t)         :: handle_
    type(string_t)         :: type_
    type(string_t)         :: vertical_profile_name_ ! Name of the vertical profile to use
    type(string_t)         :: vertical_profile_units_ ! Units for the vertical profile
    type(string_t)         :: cross_section_name_ ! Name of the absorption cross-section to use
    type(radiator_state_t) :: state_ ! Optical properties, a :f:type:`~tuvx_radiator/radiator_state_t`
    logical                :: enable_diagnostics_ ! determines if diagnostic output is written or not
    logical                :: is_air_ = .false. ! Indicates whether the radiator should be treated as "air" in optical property calculations
    type(grid_warehouse_ptr)          :: height_grid_ ! pointer to the height grid in the grid warehouse
    type(grid_warehouse_ptr)          :: wavelength_grid_ ! pointer to the wavelength grid in the grid warehouse
    type(profile_warehouse_ptr)       :: radiator_profile_ ! pointer to the radiator profile in the profile warehouse
    type(cross_section_warehouse_ptr) :: cross_section_ ! pointer to the cross section for this radiator
  contains
    ! Update radiator for new environmental conditions
    procedure :: update_state
    ! Outputs radiation diagnostics if enabled
    procedure :: output_diagnostics
    ! Returns the number of bytes needed to pack the object onto a buffer
    procedure :: pack_size
    ! Packs the object onto a character buffer
    procedure :: mpi_pack
    ! Unpacks data from a character buffer into the object
    procedure :: mpi_unpack
  end type radiator_t

  interface radiator_t
    module procedure :: constructor
  end interface

  type radiator_ptr
    ! Pointer type for building sets of radiator objects

    class(radiator_t), pointer :: val_ => null( )
  end type radiator_ptr

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function constructor( config, grid_warehouse, profile_warehouse,            &
     cross_section_warehouse ) result( new_radiator )
    ! Constructs a base_radiator_t object

    class(radiator_t),               pointer       :: new_radiator ! New :f:type:`~tuvx_radiator/radiator_t` object
    type(config_t),                  intent(inout) :: config ! Radiator configuration
    type(grid_warehouse_t),          intent(inout) :: grid_warehouse ! A :f:type:`~tuvx_grid_warehouse/grid_warehouse_t`
    type(profile_warehouse_t),       intent(inout) :: profile_warehouse ! profile warehouse
    type(cross_section_warehouse_t), intent(inout) :: cross_section_warehouse ! cross section warehouse
    allocate( new_radiator )
    call base_constructor( new_radiator, config, grid_warehouse,              &
                           profile_warehouse, cross_section_warehouse )

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine base_constructor( this, config, grid_warehouse,                  &
      profile_warehouse, cross_section_warehouse )
    ! Initializes a radiator_t object
    !
    ! This should only be called by subclasses of radiator_t

    class(radiator_t),               intent(inout) :: this ! New :f:type:`~tuvx_radiator/radiator_t` object
    type(config_t),                  intent(inout) :: config ! Radiator configuration
    type(grid_warehouse_t),          intent(inout) :: grid_warehouse ! A :f:type:`~tuvx_grid_warehouse/grid_warehouse_t`
    type(profile_warehouse_t),       intent(inout) :: profile_warehouse ! profile warehouse
    type(cross_section_warehouse_t), intent(inout) :: cross_section_warehouse ! cross section warehouse

    ! local variables
    character(len=*), parameter   :: Iam = "Base radiator constructor"
    class(grid_t),    pointer     :: z_grid, lambda_grid
    type(string_t)                :: required_keys(5), optional_keys(2)

    required_keys(1) = "name"
    required_keys(2) = "type"
    required_keys(3) = "cross section"
    required_keys(4) = "vertical profile"
    required_keys(5) = "vertical profile units"
    optional_keys(1) = "enable diagnostics"
    optional_keys(2) = "treat as air"

    call assert_msg( 691711954,                                               &
                     config%validate( required_keys, optional_keys ),         &
                     "Bad configuration data format for "//                   &
                     "base radiator." )

    call config%get( 'name',             this%handle_,                Iam )
    call config%get( 'type',             this%type_,                Iam )
    call config%get( 'vertical profile', this%vertical_profile_name_, Iam )
    call config%get( 'vertical profile units', this%vertical_profile_units_,  &
                     Iam )
    call config%get( 'cross section',    this%cross_section_name_,    Iam )

    call config%get( 'enable diagnostics', this%enable_diagnostics_, Iam,     &
      default=.false. )
    call config%get( 'treat as air', this%is_air_, Iam, default = .false. )

    this%height_grid_ = grid_warehouse%get_ptr( "height", "km" )
    this%wavelength_grid_ = grid_warehouse%get_ptr( "wavelength", "nm" )
    this%radiator_profile_ =                                                  &
        profile_warehouse%get_ptr( this%vertical_profile_name_,               &
                                   this%vertical_profile_units_ )
    this%cross_section_ =                                                     &
        cross_section_warehouse%get_ptr( this%cross_section_name_ )
    z_grid => grid_warehouse%get_grid( this%height_grid_ )
    lambda_grid => grid_warehouse%get_grid( this%wavelength_grid_ )

    ! allocate radiator state variables
    allocate( this%state_%layer_OD_(  z_grid%ncells_, lambda_grid%ncells_ ) )
    allocate( this%state_%layer_SSA_( z_grid%ncells_, lambda_grid%ncells_ ) )
    allocate( this%state_%layer_G_(   z_grid%ncells_, lambda_grid%ncells_, 1 ) )

    this%state_%layer_OD_( :,:) = 0.0_dk
    this%state_%layer_SSA_(:,:) = 0.0_dk
    this%state_%layer_G_(:,:,:) = 0.0_dk

    deallocate( z_grid      )
    deallocate( lambda_grid )

  end subroutine base_constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine update_state( this, grid_warehouse, profile_warehouse,           &
      cross_section_warehouse )
    ! Update radiator state

    class(radiator_t),               intent(inout) :: this ! A :f:type:`~tuvx_radiator/radiator_state_t`
    type(grid_warehouse_t),          intent(inout) :: grid_warehouse ! A :f:type:`~tuvx_grid_warehouse/grid_warehouse_t`
    type(profile_warehouse_t),       intent(inout) :: profile_warehouse ! A :f:type:`~tuvx_profile_warehouse/profile_warehouse_t`
    type(cross_section_warehouse_t), intent(inout) :: cross_section_warehouse ! A :f:type:`~tuvx_cross_section_warehouse/cross_section_warehouse_t`

    ! Local variables
    character(len=*), parameter     :: Iam = 'Base radiator update state'
    real(dk) ,        parameter     :: km2cm = 1.e5_dk
    integer                         :: w_index
    real(dk),         allocatable   :: cross_section(:,:)
    class(grid_t),          pointer :: z_grid
    class(grid_t),          pointer :: lambda_grid
    class(profile_t),       pointer :: radiator_profile
    class(cross_section_t), pointer :: radiator_cross_section

    ! get specific grids and profiles
    z_grid => grid_warehouse%get_grid( this%height_grid_ )
    lambda_grid => grid_warehouse%get_grid( this%wavelength_grid_ )

    radiator_profile =>                                                       &
      profile_warehouse%get_profile( this%radiator_profile_ )

    radiator_cross_section =>                                                 &
      cross_section_warehouse%get( this%cross_section_ )

    ! check radiator state type allocation
    call assert_msg( 345645215, allocated( this%state_%layer_OD_ ),           &
                     "Radiator state not allocated" )

    ! check that the profile has a set of layer densities avaialble
    call assert_msg( 108830786, allocated( radiator_profile%layer_dens_ ),    &
                     "Radiator profiles must provide layer densities" )

    ! set radiator state members
    cross_section = radiator_cross_section%calculate( grid_warehouse,         &
                                                      profile_warehouse,      &
                                                      at_mid_point = .true. )
    call diagout( 'o2xs.new', cross_section, this%enable_diagnostics_ )
    do w_index = 1,lambda_grid%ncells_
      this%state_%layer_OD_(:,w_index) = radiator_profile%layer_dens_         &
                                         * cross_section(:,w_index)
    enddo

    ! Settings for a gas phase radiator
    if( this%is_air_ ) then
      this%state_%layer_SSA_ = 1._dk
      this%state_%layer_G_   = 0._dk
    else
      this%state_%layer_SSA_ = 0._dk
      this%state_%layer_G_   = 0._dk
    endif

    deallocate( z_grid )
    deallocate( lambda_grid )
    deallocate( radiator_profile )
    deallocate( radiator_cross_section )

  end subroutine update_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine output_diagnostics( this )

    class(radiator_t), intent(in) :: this ! A :f:type:`~tuvx_radiator/radiator_state_t`
    character(len=:), allocatable :: filename

    select case( this%handle_%to_char( ) )
      case( 'air' )
        filename = 'dtrl.new'
      case( 'Aerosols' )
        filename = 'dtaer.new'
      case( 'O3' )
        filename = 'dto3.new'
      case( 'O2' )
        filename = 'dto2.new'
    case default
        filename = this%type_%to_char( )
    end select

    call diagout( filename, this%state_%layer_OD_, this%enable_diagnostics_ )

  end subroutine output_diagnostics

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function pack_size( this, comm )
    ! Returns the size of a character buffer required to pack the radiator


    class(radiator_t), intent(in) :: this ! radiator to be packed
    integer,           intent(in) :: comm ! MPI communicator

#ifdef MUSICA_USE_MPI
    pack_size = this%handle_%pack_size( comm ) +                              &
                this%type_%pack_size( comm ) +                                &
                this%vertical_profile_name_%pack_size( comm ) +               &
                this%vertical_profile_units_%pack_size( comm ) +              &
                this%cross_section_name_%pack_size( comm ) +                  &
                this%state_%pack_size( comm ) +                               &
                musica_mpi_pack_size( this%enable_diagnostics_, comm ) +      &
                musica_mpi_pack_size( this%is_air_, comm ) +                  &
                this%height_grid_%pack_size( comm ) +                         &
                this%wavelength_grid_%pack_size( comm ) +                     &
                this%radiator_profile_%pack_size( comm ) +                    &
                this%cross_section_%pack_size( comm )
#else
    pack_size = 0
#endif

  end function pack_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine mpi_pack( this, buffer, position, comm )
    ! Packs the radiator onto a character buffer

    class(radiator_t), intent(in)    :: this      ! radiator to be packed
    character,         intent(inout) :: buffer(:) ! memory buffer
    integer,           intent(inout) :: position  ! current buffer position
    integer,           intent(in)    :: comm      ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: prev_pos

    prev_pos = position
    call this%handle_%mpi_pack(                 buffer, position, comm )
    call this%type_%mpi_pack(                   buffer, position, comm )
    call this%vertical_profile_name_%mpi_pack(  buffer, position, comm )
    call this%vertical_profile_units_%mpi_pack( buffer, position, comm )
    call this%cross_section_name_%mpi_pack(     buffer, position, comm )
    call this%state_%mpi_pack(                  buffer, position, comm )
    call musica_mpi_pack( buffer, position, this%enable_diagnostics_, comm )
    call musica_mpi_pack( buffer, position, this%is_air_,             comm )
    call this%height_grid_%mpi_pack(            buffer, position, comm )
    call this%wavelength_grid_%mpi_pack(        buffer, position, comm )
    call this%radiator_profile_%mpi_pack(       buffer, position, comm )
    call this%cross_section_%mpi_pack(          buffer, position, comm )
    call assert( 449676235, position - prev_pos <= this%pack_size( comm ) )
#endif

  end subroutine mpi_pack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine mpi_unpack( this, buffer, position, comm )
    ! Unpacks a radiator from a character buffer

    class(radiator_t), intent(out)   :: this      ! radiator to be unpacked
    character,         intent(inout) :: buffer(:) ! memory buffer
    integer,           intent(inout) :: position  ! current buffer position
    integer,           intent(in)    :: comm      ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: prev_pos

    prev_pos = position
    call this%handle_%mpi_unpack(                 buffer, position, comm )
    call this%type_%mpi_unpack(                   buffer, position, comm )
    call this%vertical_profile_name_%mpi_unpack(  buffer, position, comm )
    call this%vertical_profile_units_%mpi_unpack( buffer, position, comm )
    call this%cross_section_name_%mpi_unpack(     buffer, position, comm )
    call this%state_%mpi_unpack(                  buffer, position, comm )
    call musica_mpi_unpack( buffer, position, this%enable_diagnostics_, comm )
    call musica_mpi_unpack( buffer, position, this%is_air_,             comm )
    call this%height_grid_%mpi_unpack(            buffer, position, comm )
    call this%wavelength_grid_%mpi_unpack(        buffer, position, comm )
    call this%radiator_profile_%mpi_unpack(       buffer, position, comm )
    call this%cross_section_%mpi_unpack(          buffer, position, comm )
    call assert( 216868760, position - prev_pos <= this%pack_size( comm ) )
#endif

  end subroutine mpi_unpack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine accumulate( this, radiators )
    ! Create a single radiator state that corresponds to the cumulative
    ! state of a set of radiators, such that the optical properties of the
    ! accumulated state can be used to solve radiative transfer equations for
    ! the set of radiators
    !
    ! Optical properties for radiators configured to 'treat as air' are
    ! unique.

    class(radiator_state_t), intent(inout) :: this
    class(radiator_ptr),     intent(in)    :: radiators(:)

    real(dk), parameter :: kfloor = 1.0_dk / largest ! smallest value for radiative properties
    real(dk), parameter :: kair_asym_factor = 0.1_dk

    integer :: i_radiator, i_stream, n_streams
    real(dk), allocatable :: dscat(:,:)
    real(dk), allocatable :: dscat_accum(:,:)
    real(dk), allocatable :: dabs_accum(:,:)
    real(dk), allocatable :: asym_accum(:,:,:)

    allocate( dscat,       mold = radiators(1)%val_%state_%layer_OD_ )
    allocate( dscat_accum, mold = dscat )
    allocate( dabs_accum,  mold = dscat )
    allocate( asym_accum,  mold = this%layer_G_ )

    dscat_accum = 0.0_dk
    dabs_accum  = 0.0_dk
    asym_accum  = 0.0_dk

    n_streams = size( this%layer_G_,3 )

    ! iterate over radiators accumulating radiative properties
    do i_radiator = 1, size( radiators )
      associate( OD     => radiators( i_radiator )%val_%state_%layer_OD_,     &
                 SSA    => radiators( i_radiator )%val_%state_%layer_SSA_,    &
                 G      => radiators( i_radiator )%val_%state_%layer_G_,      &
                 is_air => radiators( i_radiator )%val_%is_air_ )
        dscat       = OD * SSA
        dscat_accum = dscat_accum + dscat
        dabs_accum  = dabs_accum + OD * ( 1.0_dk - SSA )
        if( .not. is_air ) then
          do i_stream = 1, n_streams
            asym_accum(:,:,i_stream)  = asym_accum(:,:,i_stream) &
                                      + G(:,:,1)**i_stream * dscat
          end do
        else
        ! rayleigh scattering
          if( n_streams >= 2 ) then
            asym_accum(:,:,2)  = asym_accum(:,:,2) + kair_asym_factor * dscat
          else
            asym_accum(:,:,1)  = asym_accum(:,:,1) + G(:,:,1) * dscat
          end if
        end if
      end associate
    end do

    ! set atmosphere radiative properties
    dscat_accum = max( dscat_accum, kfloor )
    dabs_accum  = max( dabs_accum, kfloor )

    this%layer_OD_ = dscat_accum + dabs_accum
    if( .not. allocated( this%layer_SSA_ ) ) then
      allocate( this%layer_SSA_, mold = this%layer_OD_ )
    end if
    where( dscat_accum == kfloor )
      this%layer_SSA_ = kfloor
    elsewhere
      this%layer_SSA_ = dscat_accum / this%layer_OD_
    endwhere

    if( n_streams >= 2 ) then
      this%layer_SSA_ = max( min( this%layer_SSA_,1._dk - precis ),precis )
      this%layer_OD_  = max( this%layer_OD_,precis )
    end if

    do i_stream = 1, n_streams
      this%layer_G_(:,:,i_stream) = asym_accum(:,:,i_stream) / dscat_accum
    end do

  end subroutine accumulate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function state_pack_size( this, comm ) result( pack_size )
    ! Returns the size of a character buffer required to pack the radiator
    ! state

    use musica_mpi,                    only : musica_mpi_pack_size

    class(radiator_state_t), intent(in) :: this ! radiator state to be packed
    integer,                 intent(in) :: comm ! MPI communicator

#ifdef MUSICA_USE_MPI
    pack_size = musica_mpi_pack_size( this%layer_OD_,  comm ) +               &
                musica_mpi_pack_size( this%layer_SSA_, comm ) +               &
                musica_mpi_pack_size( this%layer_G_,   comm )
#else
    pack_size = 0
#endif

  end function state_pack_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine state_mpi_pack( this, buffer, position, comm )
    ! Packs the radiator state onto a character buffer

    class(radiator_state_t), intent(in)    :: this      ! radiator state to be packed
    character,               intent(inout) :: buffer(:) ! memory buffer
    integer,                 intent(inout) :: position  ! current buffer position
    integer,                 intent(in)    :: comm      ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: prev_pos

    prev_pos = position
    call musica_mpi_pack( buffer, position, this%layer_OD_,  comm )
    call musica_mpi_pack( buffer, position, this%layer_SSA_, comm )
    call musica_mpi_pack( buffer, position, this%layer_G_,   comm )
    call assert( 942613664, position - prev_pos <= this%pack_size( comm ) )
#endif

  end subroutine state_mpi_pack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine state_mpi_unpack( this, buffer, position, comm )
    ! Unpacks a radiator state from a character buffer

    class(radiator_state_t), intent(out)   :: this      ! radiator state to be unpacked
    character,               intent(inout) :: buffer(:) ! memory buffer
    integer,                 intent(inout) :: position  ! current buffer position
    integer,                 intent(in)    :: comm      ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: prev_pos

    prev_pos = position
    call musica_mpi_unpack( buffer, position, this%layer_OD_,  comm )
    call musica_mpi_unpack( buffer, position, this%layer_SSA_, comm )
    call musica_mpi_unpack( buffer, position, this%layer_G_,   comm )
    call assert( 709806189, position - prev_pos <= this%pack_size( comm ) )
#endif

  end subroutine state_mpi_unpack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine finalize( this )
    ! Finalizes a radiator state object

    type(radiator_state_t) :: this ! A :f:type:`~tuvx_radiator/radiator_state_t`

    if( allocated( this%layer_OD_ ) ) then
      deallocate( this%layer_OD_ )
    endif
    if( allocated( this%layer_SSA_ ) ) then
      deallocate( this%layer_SSA_ )
    endif
    if( allocated( this%layer_G_ ) ) then
      deallocate( this%layer_G_ )
    endif

  end subroutine finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_radiator
