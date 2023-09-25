! Copyright (C) 2021 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0

module tuvx_radiative_transfer
! A calculator for atmospheric radiation

  use musica_config,                   only : config_t
  use musica_constants,                only : dk => musica_dk
  use tuvx_cross_section_warehouse,    only : cross_section_warehouse_t
  use tuvx_profile_warehouse,          only : profile_warehouse_ptr
  use tuvx_radiator_warehouse,         only : radiator_warehouse_t,           &
                                              radiator_warehouse_ptr
  use tuvx_solver,                     only : solver_t

  implicit none
  private

  public :: radiative_transfer_t

  type radiative_transfer_t
    ! radXfer component core
    !
    ! Calculates the atmospheric radiation field
    private
    class(solver_t), pointer                    :: solver_ => null()
    type(cross_section_warehouse_t),    pointer :: cross_section_warehouse_   &
                                                       => null( ) ! Copy of original radiative transfer :f:type:`~tuvx_cross_section_warehouse/cross_section_warehouse_t`
    type(radiator_warehouse_t), public, pointer :: radiator_warehouse_        &
                                                       => null( ) ! Copy of original :f:type:`~tuvx_radiator/radiator_warehouse_t`
    logical                                     :: O2_exists_   ! indicates whether O2 exists as a profile
    type(radiator_warehouse_ptr)                :: O2_radiator_ ! pointer to the O2 radiator in the radiator warehouse
    type(profile_warehouse_ptr)                 :: air_profile_ ! pointer to the air profile in the profile warehouse
  contains
    procedure :: name => component_name
    procedure :: description
    procedure :: calculate
    ! Returns an updater for a radiator in the warehouse
    procedure :: get_radiator_updater
    ! Returns the number of bytes needed to pack the object onto a buffer
    procedure :: pack_size
    ! Packs the object onto a character buffer
    procedure :: mpi_pack
    ! Unpacks data from a character buffer into the object
    procedure :: mpi_unpack
    final :: finalize
  end type radiative_transfer_t

  interface radiative_transfer_t
    ! Constructor
    module procedure constructor
  end interface radiative_transfer_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function constructor( config, grid_warehouse, profile_warehouse, radiators )&
      result( this )
    ! Initializes the components necessary to solve radiative transfer

    use musica_assert,                 only : assert_msg, die_msg
    use musica_string,                 only : string_t
    use tuvx_grid_warehouse,           only : grid_warehouse_t
    use tuvx_profile_warehouse,        only : profile_warehouse_t
    use tuvx_solver_factory,           only : solver_builder

    type(radiative_transfer_t), pointer :: this ! New :f:type:`~tuvx_radiative_transfer/radxfer_component_core_t`
    type(config_t),                        intent(inout) :: config            ! radXfer configuration data
    type(grid_warehouse_t),                intent(inout) :: grid_warehouse    ! A :f:type:`~tuvx_grid_warehouse/grid_warehouse_t`
    type(profile_warehouse_t),             intent(inout) :: profile_warehouse ! A :f:type:`~tuvx_profile_warehouse/profile_warehouse_t`
    class(radiator_warehouse_t), optional, intent(in)    :: radiators         ! Radiators to include in the configuration

    character(len=*), parameter :: Iam = 'radiative transfer constructor: '
    type(config_t) :: solver_config
    type(config_t) :: child_config
    type(string_t) :: required_keys(2), optional_keys(1)

    required_keys(1) = "cross sections"
    required_keys(2) = "radiators"
    optional_keys(1) = "solver"

    call assert_msg( 817033232,                                               &
                     config%validate( required_keys, optional_keys ),         &
                     "Bad configuration data format for "//                   &
                     "radiative transfer." )

    allocate( this )

    ! instantiate and initialize the radXfer cross section warehouse
    call config%get( "cross sections", child_config, Iam )
    this%cross_section_warehouse_ =>                        &
        cross_section_warehouse_t( child_config, grid_warehouse,              &
                                   profile_warehouse )

    ! instantiate and initialize the radiator warehouse
    call config%get( "radiators", child_config, Iam )
    this%radiator_warehouse_ =>                             &
        radiator_warehouse_t( child_config, grid_warehouse, profile_warehouse,&
                              this%cross_section_warehouse_ )
    if( present( radiators ) ) call this%radiator_warehouse_%add( radiators )

    ! set up pointers to radiators and profiles
    this%O2_exists_ = this%radiator_warehouse_%exists( "O2" )
    if( this%O2_exists_ ) then
      this%O2_radiator_ = this%radiator_warehouse_%get_ptr( "O2" )
      this%air_profile_ =                                     &
          profile_warehouse%get_ptr( "air", "molecule cm-3" )
    end if

    ! get the radiative transfer solver
    call config%get( "solver", solver_config, Iam )
    this%solver_ => solver_builder( solver_config, grid_warehouse,            &
                                    profile_warehouse )

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  type(string_t) function component_name( this )
    ! Model component name

    use musica_string,                 only : string_t

    class(radiative_transfer_t), intent(in) :: this ! A :f:type:`~tuvx_radiative_transfer/radxfer_component_core_t`

    component_name = "radiative transfer"

  end function component_name

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  type(string_t) function description( this )
    ! Model component description

    use musica_string,                 only : string_t

    class(radiative_transfer_t), intent(in) :: this ! A :f:type:`~tuvx_radiative_transfer/radxfer_component_core_t`

    description = "calculates the radiation field"

  end function description

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine calculate( this, la_srb, spherical_geometry, grid_warehouse,     &
      profile_warehouse, radiation_field )
    ! Calculate the radiation field

    use musica_assert,                 only : die_msg
    use tuvx_grid_warehouse,           only : grid_warehouse_t
    use tuvx_radiator_warehouse,       only : warehouse_iterator_t
    use tuvx_radiator,                 only : radiator_t
    use tuvx_radiator,                 only : radiator_state_t
    use tuvx_profile,                  only : profile_t
    use tuvx_profile_warehouse,        only : profile_warehouse_t
    use tuvx_spherical_geometry,       only : spherical_geometry_t
    use tuvx_la_sr_bands,              only : la_sr_bands_t
    use tuvx_solver,                   only : radiation_field_t

    class(radiative_transfer_t),       intent(inout) :: this               ! A :f:type:`~tuvx_radiative_transfer/radxfer_component_core_t`
    type(grid_warehouse_t),            intent(inout) :: grid_warehouse     ! :f:type:`~tuvx_grid_warehouse/grid_warehouse_t`
    type(profile_warehouse_t),         intent(inout) :: profile_warehouse  ! A :f:type:`~tuvx_profile_warehouse/profile_warehouse_t`
    type(spherical_geometry_t),        intent(inout) :: spherical_geometry ! A :f:type:`~tuvx_spherical_geometry/spherical_geometry_t`
    type(la_sr_bands_t),               intent(inout) :: la_srb             ! A :f:type:`~tuvx_la_sr_bands/la_sr_bands_t`

    type(radiation_field_t), pointer, intent(out)   :: radiation_field

    ! Local variables
    character(len=*), parameter          :: Iam = 'radXfer component calculate: '

    integer                              :: nlyr
    integer                              :: radNdx
    real(dk)                             :: zenithAngle
    real(dk), allocatable                :: airVcol(:), airScol(:)
    type(warehouse_iterator_t), pointer  :: iter
    class(radiator_t),          pointer  :: aRadiator
    class(profile_t),           pointer  :: airprofile

    ! iterate over radiators
    iter => this%radiator_warehouse_%get_iterator( )
    do while( iter%next( ) )
      aRadiator => this%radiator_warehouse_%get_radiator( iter )
      call aRadiator%update_state( grid_warehouse, profile_warehouse,         &
                                   this%cross_section_warehouse_ )
    enddo
    deallocate( iter )

    ! look for O2 radiator; Lyman Alpha and SR bands
    if( this%O2_exists_ ) then
      aRadiator => this%radiator_warehouse_%get_radiator( this%O2_radiator_ )
      airprofile => profile_warehouse%get_profile( this%air_profile_ )
      allocate( airVcol( airprofile%ncells_ ),                                &
                airScol( airprofile%ncells_ + 1 ) )
      call spherical_geometry%air_mass( airprofile%exo_layer_dens_, airVcol,  &
                                        airScol )
      call la_srb%optical_depth( grid_warehouse, profile_warehouse, airVcol,  &
                                 airScol, aRadiator%state_%layer_OD_,         &
                                 spherical_geometry )
      deallocate( airVcol, airScol )
      deallocate( airprofile )
    endif

    nlyr = size( aRadiator%state_%layer_OD_, dim = 1 )

    zenithAngle = spherical_geometry%solar_zenith_angle_
    associate( theSolver => this%solver_ )
    radiation_field => theSolver%update_radiation_field(                      &
                     zenithAngle, nlyr, spherical_geometry,                   &
                     grid_warehouse, profile_warehouse,                       &
                     this%radiator_warehouse_ )
    end associate

  end subroutine calculate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function get_radiator_updater( this, radiator, found ) result( updater )
    ! Returns an updater for a radiator that a host application can use to
    ! update TUV-x state data at runtime
    !
    ! If the optional `found` flag is omitted, an error is returned if the
    ! radiator does not exist in TUV-x

    use musica_assert,                 only : assert_msg
    use tuvx_radiator,                 only : radiator_t
    use tuvx_radiator_from_host,       only : radiator_updater_t

    class(radiative_transfer_t), intent(in)  :: this     ! Radiative transfer calculator
    class(radiator_t),           intent(in)  :: radiator ! The radiator to get an updater for
    logical, optional,           intent(out) :: found    ! Flag indicating whether the
                                                         ! radiator was found
    type(radiator_updater_t)       :: updater

    call assert_msg( 101866398, associated( this%radiator_warehouse_ ),       &
                     "Radiators not available" )
    updater = this%radiator_warehouse_%get_updater( radiator, found )

  end function get_radiator_updater

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function pack_size( this, comm )
    ! Returns the size of a character buffer required to pack the radiative
    ! transfer calculator

    use musica_mpi,                    only : musica_mpi_pack_size
    use musica_string,                 only : string_t
    use tuvx_solver_factory,           only : solver_type_name

    class(radiative_transfer_t), intent(inout) :: this ! radiative transfer to be packed
    integer,                     intent(in)    :: comm ! MPI communicator

#ifdef MUSICA_USE_MPI
    type(string_t) :: solver_type

    solver_type = solver_type_name( this%solver_ )
    pack_size = solver_type%pack_size( comm ) +                               &
                this%solver_%pack_size( comm ) +                              &
                this%cross_section_warehouse_%pack_size( comm ) +             &
                this%radiator_warehouse_%pack_size( comm ) +                  &
                musica_mpi_pack_size( this%O2_exists_, comm ) +               &
                this%O2_radiator_%pack_size( comm ) +                         &
                this%air_profile_%pack_size( comm )
#else
    pack_size = 0
#endif

  end function pack_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine mpi_pack( this, buffer, position, comm )
    ! Packs the radiative transfer calculator onto a character buffer

    use musica_assert,                 only : assert
    use musica_mpi,                    only : musica_mpi_pack
    use musica_string,                 only : string_t
    use tuvx_solver_factory,           only : solver_type_name

    class(radiative_transfer_t), intent(inout)    :: this   ! radiative transfer to be packed
    character,                   intent(inout) :: buffer(:) ! memory buffer
    integer,                     intent(inout) :: position  ! current buffer position
    integer,                     intent(in)    :: comm      ! MPI communicator

#ifdef MUSICA_USE_MPI
    type(string_t) :: solver_type
    integer :: prev_pos
    prev_pos = position

    solver_type = solver_type_name( this%solver_ )
    call solver_type%mpi_pack( buffer, position , comm )
    call this%solver_%mpi_pack( buffer, position, comm )
    call this%cross_section_warehouse_%mpi_pack( buffer, position, comm )
    call this%radiator_warehouse_%mpi_pack( buffer, position, comm )
    call musica_mpi_pack( buffer, position, this%O2_exists_, comm )
    call this%O2_radiator_%mpi_pack( buffer, position, comm )
    call this%air_profile_%mpi_pack( buffer, position, comm )

    call assert( 742641642, position - prev_pos <= this%pack_size( comm ) )
#endif

  end subroutine mpi_pack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine mpi_unpack( this, buffer, position, comm )
    ! Unpacks a radiative transfer calculator from a character buffer

    use musica_assert,                 only : assert, die_msg
    use musica_string,                 only : string_t
    use musica_mpi,                    only : musica_mpi_unpack
    use tuvx_solver_factory,           only : solver_allocate

    class(radiative_transfer_t), intent(out)   :: this      ! radiative transfer to be packed
    character,                   intent(inout) :: buffer(:) ! memory buffer
    integer,                     intent(inout) :: position  ! current buffer position
    integer,                     intent(in)    :: comm      ! MPI communicator

#ifdef MUSICA_USE_MPI
    type(string_t) :: solver_type
    integer :: prev_pos

    prev_pos = position

    allocate( this%cross_section_warehouse_ )
    allocate( this%radiator_warehouse_ )

    call solver_type%mpi_unpack( buffer, position, comm )
    this%solver_ => solver_allocate( solver_type )
    call this%solver_%mpi_unpack( buffer, position, comm )
    call this%cross_section_warehouse_%mpi_unpack( buffer, position, comm )
    call this%radiator_warehouse_%mpi_unpack( buffer, position, comm )
    call musica_mpi_unpack( buffer, position, this%O2_exists_, comm )
    call this%O2_radiator_%mpi_unpack( buffer, position, comm )
    call this%air_profile_%mpi_unpack( buffer, position, comm )

    call assert( 559826176, position - prev_pos <= this%pack_size( comm ) )
#endif

  end subroutine mpi_unpack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine finalize( this )
    ! Frees the warehouse and solver data associated with this class

    type(radiative_transfer_t), intent(inout) :: this ! This :f:type:`~tuvx_radiative_transfer/radxfer_component_core_t`

    if( associated( this%cross_section_warehouse_ ) )                         &
        deallocate( this%cross_section_warehouse_ )
    if( associated( this%radiator_warehouse_ ) )                              &
        deallocate( this%radiator_warehouse_ )
    if( associated( this%solver_ ) ) deallocate( this%solver_ )

  end subroutine finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_radiative_transfer
