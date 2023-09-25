! Copyright (C) 2021 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
module tuvx_core
  ! Top-level TUV-x interface

  use musica_config,                   only : config_t
  use musica_string,                   only : string_t
  use musica_constants,                only : dk => musica_dk
  use tuvx_dose_rates,                 only : dose_rates_t
  use tuvx_grid_warehouse,             only : grid_warehouse_t
  use tuvx_la_sr_bands,                only : la_sr_bands_t
  use tuvx_photolysis_rates,           only : photolysis_rates_t
  use tuvx_profile_warehouse,          only : profile_warehouse_t
  use tuvx_radiative_transfer,         only : radiative_transfer_t
  use tuvx_solver,                     only : radiation_field_t
  use tuvx_spherical_geometry,         only : spherical_geometry_t

  implicit none

  private
  public :: core_t

  type :: core_t
    private
    ! The TUV-x core_t class defines the API for interactions
    ! with a host application
    type(grid_warehouse_t),      pointer :: grid_warehouse_ => null()
    type(profile_warehouse_t),   pointer :: profile_warehouse_ => null()
    type(spherical_geometry_t),  pointer :: spherical_geometry_ => null()
    type(la_sr_bands_t),         pointer :: la_sr_bands_ => null()
    type(radiative_transfer_t),  pointer :: radiative_transfer_ => null()
    type(photolysis_rates_t),    pointer :: photolysis_rates_ => null()
    type(dose_rates_t),          pointer :: dose_rates_ => null()
    type(radiation_field_t),     pointer :: radiation_field_ => null()
    logical                              :: enable_diagnostics_ ! determines if diagnostic output is written or not
  contains
    ! Calculate photolysis rate constants and dose rates
    procedure :: run
    ! Returns a grid from the warehouse
    procedure :: get_grid
    ! Returns a profile from the warehouse
    procedure :: get_profile
    ! Returns an updater for use TUV-x data
    procedure, private :: get_grid_updater, get_profile_updater,              &
                          get_radiator_updater
    generic :: get_updater => get_grid_updater, get_profile_updater,          &
                              get_radiator_updater
    ! Returns the number of photolysis reactions
    procedure :: number_of_photolysis_reactions
    ! Returns the number of dose rates
    procedure :: number_of_dose_rates
    ! Returns the set of photolysis reaction labels
    procedure :: photolysis_reaction_labels
    ! Returns the set of dose rate labels
    procedure :: dose_rate_labels
    ! Returns the photolysis reaction cross section for the current conditions
    procedure :: get_photolysis_cross_section
    ! Returns the photolysis reaction quantum yield for the current conditions
    procedure :: get_photolysis_quantum_yield
    ! Returns the radiation field for the current conditions
    procedure :: get_radiation_field
    ! Returns the number of bytes required to pack the core onto a buffer
    procedure :: pack_size
    ! Packs the core onto a character buffer
    procedure :: mpi_pack
    ! Unpacks a core from a character buffer
    procedure :: mpi_unpack
    final     :: finalize
  end type core_t

  interface core_t
    module procedure constructor
  end interface core_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function constructor( config, grids, profiles, radiators ) result( new_core )
    ! Constructor of TUV-x core objects

    use musica_assert,                 only : assert_msg
    use musica_string,                 only : string_t
    use tuvx_diagnostic_util,          only : diagout
    use tuvx_profile,                  only : profile_t
    use tuvx_radiator_warehouse,       only : radiator_warehouse_t

    type(string_t),                        intent(in) :: config    ! Full TUV-x configuration data
    class(grid_warehouse_t),     optional, intent(in) :: grids     ! Set of grids to include in the configuration
    class(profile_warehouse_t),  optional, intent(in) :: profiles  ! Set of profiles to include in the configuration
    class(radiator_warehouse_t), optional, intent(in) :: radiators ! Set of radiators to include in the configuration
    class(core_t),                         pointer    :: new_core

    ! Local variables
    character(len=*), parameter :: Iam = 'Photolysis core constructor: '
    logical                     :: found
    type(config_t)              :: core_config, child_config
    class(profile_t),  pointer  :: aprofile
    type(string_t)              :: required_keys(4), optional_keys(3)

    call core_config%from_file( config%to_char() )

    ! Check json configuration file for basic structure, integrity
    required_keys(1) = "radiative transfer"
    required_keys(2) = "grids"
    required_keys(3) = "profiles"
    required_keys(4) = "O2 absorption"
    optional_keys(1) = "photolysis"
    optional_keys(2) = "dose rates"
    optional_keys(3) = "enable diagnostics"
    call assert_msg( 255400232,                                               &
                     core_config%validate( required_keys, optional_keys ),    &
                     "Bad configuration data format for tuv-x core." )

    ! Instantiate photolysis core
    allocate( new_core )

    call core_config%get( 'enable diagnostics', new_core%enable_diagnostics_,  &
      Iam, default=.false. )

    ! Instantiate and initialize grid warehouse
    call core_config%get( "grids", child_config, Iam )
    new_core%grid_warehouse_ => grid_warehouse_t( child_config )
    if( present( grids ) ) call new_core%grid_warehouse_%add( grids )

    ! Instantiate and initialize profile warehouse
    call core_config%get( "profiles", child_config, Iam )
    new_core%profile_warehouse_ =>                                            &
       profile_warehouse_t( child_config, new_core%grid_warehouse_ )
     if( present( profiles ) ) call new_core%profile_warehouse_%add( profiles )

    aprofile => new_core%profile_warehouse_%get_profile( "temperature", "K" )
    call diagout( 'vptmp.new', aprofile%edge_val_,                            &
      new_core%enable_diagnostics_ )
    deallocate( aprofile )

    aprofile => new_core%profile_warehouse_%get_profile( "air",               &
                                                         "molecule cm-3" )
    call diagout( 'vpair.new', aprofile%edge_val_,                            &
      new_core%enable_diagnostics_  )
    deallocate( aprofile )

    if( new_core%profile_warehouse_%exists( "O3", "molecule cm-3" ) ) then
      aprofile => new_core%profile_warehouse_%get_profile( "O3",              &
                                                           "molecule cm-3" )
      call diagout( 'vpco3.new', aprofile%layer_dens_,                        &
        new_core%enable_diagnostics_  )
      deallocate( aprofile )
    end if

    ! Set up radiative transfer calculator
    call core_config%get( "radiative transfer", child_config, Iam )
    new_core%radiative_transfer_ => &
        radiative_transfer_t( child_config,                                   &
                              new_core%grid_warehouse_,                       &
                              new_core%profile_warehouse_,                    &
                              radiators )

    ! photolysis rate constants
    call core_config%get( "photolysis", child_config, Iam,          &
                          found = found )
    if( found ) then
      new_core%photolysis_rates_ => &
          photolysis_rates_t( child_config,                                   &
                              new_core%grid_warehouse_,                       &
                              new_core%profile_warehouse_ )
    end if

    ! dose rates
    call core_config%get( "dose rates", child_config, Iam, found = found )
    if( found ) then
      new_core%dose_rates_ => &
          dose_rates_t( child_config, new_core%grid_warehouse_,               &
                        new_core%profile_warehouse_ )
    end if

    ! instantiate and initialize spherical geometry type
    new_core%spherical_geometry_ =>                                           &
        spherical_geometry_t( new_core%grid_warehouse_ )

    ! instantiate and initialize lyman alpha, srb type
    call core_config%get( "O2 absorption", child_config, Iam )
    new_core%la_sr_bands_ => la_sr_bands_t( child_config,                     &
                                            new_core%grid_warehouse_,         &
                                            new_core%profile_warehouse_ )


  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine run( this, solar_zenith_angle, earth_sun_distance,               &
      photolysis_rate_constants, dose_rates, diagnostic_label )
    ! Performs calculations for specified photolysis and dose rates for a
    ! given set of conditions

    use tuvx_profile,                    only : profile_t
    use tuvx_radiator,                   only : radiator_t
    use tuvx_diagnostic_util,            only : diagout
    use tuvx_radiator_warehouse,         only : warehouse_iterator_t

    class(core_t),              intent(inout) :: this ! TUV-x core
    real(dk),                   intent(in)    :: solar_zenith_angle             ! [degrees]
    real(dk),                   intent(in)    :: earth_sun_distance             ! [AU]
    real(dk),         optional, intent(out)   :: photolysis_rate_constants(:,:) ! (vertical level, reaction) [s-1]
    real(dk),         optional, intent(out)   :: dose_rates(:,:)                ! (vertical level, reaction) [s-1]
    character(len=*), optional, intent(in)    :: diagnostic_label               ! label used in diagnostic file names

    ! Local variables
    character(len=*), parameter         :: Iam = 'Photolysis core run: '
    character(len=2)                    :: number
    class(radiator_t),          pointer :: radiator
    type(warehouse_iterator_t), pointer :: warehouse_iter
    character(len=:), allocatable       :: diag_label

    if( present( diagnostic_label ) ) then
      diag_label = diagnostic_label
    else
      diag_label = ""
    end if

    ! calculate the radiation field
    if( associated( this%radiation_field_ ) )                                 &
        deallocate( this%radiation_field_ )
    call this%spherical_geometry_%set_parameters( solar_zenith_angle,         &
                                                  this%grid_warehouse_ )
    call this%radiative_transfer_%calculate( this%la_sr_bands_,               &
                                             this%spherical_geometry_,        &
                                             this%grid_warehouse_,            &
                                             this%profile_warehouse_,         &
                                             this%radiation_field_ )
    if( this%enable_diagnostics_ ) then
      call diagout( 'radField.' // diag_label // '.new',                      &
                    this%radiation_field_%fdr_ + this%radiation_field_%fup_ + &
                    this%radiation_field_%fdn_, this%enable_diagnostics_  )
    end if
    ! scale the radiation field by the Earth-Sun distance
    call this%radiation_field_%apply_scale_factor( earth_sun_distance )
    if( associated( this%photolysis_rates_ ) .and.                            &
        present( photolysis_rate_constants ) ) then
      call this%photolysis_rates_%get( this%la_sr_bands_,                     &
                                       this%spherical_geometry_,              &
                                       this%grid_warehouse_,                  &
                                       this%profile_warehouse_,               &
                                       this%radiation_field_,                 &
                                       photolysis_rate_constants,             &
                                       diag_label )
    end if
    if( associated( this%dose_rates_ ) .and. present( dose_rates ) ) then
      call this%dose_rates_%get( this%grid_warehouse_,                        &
                                 this%profile_warehouse_,                     &
                                 this%radiation_field_,                       &
                                 dose_rates,                                  &
                                 diag_label )
    endif

    ! diagnostic output
    if( this%enable_diagnostics_ ) then
      warehouse_iter =>                                                       &
          this%radiative_transfer_%radiator_warehouse_%get_iterator( )
      do while( warehouse_iter%next( ) )
        radiator => this%radiative_transfer_%                                 &
          radiator_warehouse_%get_radiator( warehouse_iter )
        call radiator%output_diagnostics()
      enddo
      deallocate( warehouse_iter )
    end if

  end subroutine run

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function get_grid( this, grid_name, units ) result( grid )
    ! Returns a copy of a grid from the warehouse

    use musica_assert,                 only : assert_msg
    use tuvx_grid,                     only : grid_t

    class(core_t),    intent(in) :: this
    character(len=*), intent(in) :: grid_name
    character(len=*), intent(in) :: units
    class(grid_t),    pointer    :: grid

    call assert_msg( 285057977, associated( this%grid_warehouse_ ),           &
                     "Grids not available" )
    grid => this%grid_warehouse_%get_grid( grid_name, units )

  end function get_grid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function get_profile( this, profile_name, units ) result( profile )
    ! Returns a copy of a profile from the warehouse

    use musica_assert,                 only : assert_msg
    use tuvx_profile,                  only : profile_t

    class(core_t),    intent(in) :: this
    character(len=*), intent(in) :: profile_name
    character(len=*), intent(in) :: units
    class(profile_t), pointer    :: profile

    call assert_msg( 780188063, associated( this%profile_warehouse_ ),        &
                     "Profiles not available" )
    profile => this%profile_warehouse_%get_profile( profile_name, units )

  end function get_profile

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function get_grid_updater( this, grid, found ) result( updater )
    ! Returns an updater for a grid that a host application can use to update
    ! TUV-x state data at runtime
    !
    ! If the optional `found` flag is omitted, an error is returned if the
    ! grid does not exist in the TUV-x

    use musica_assert,                 only : assert_msg
    use tuvx_grid,                     only : grid_t
    use tuvx_grid_from_host,           only : grid_updater_t

    class(core_t),     intent(in)  :: this  ! TUV-x core
    class(grid_t),     intent(in)  :: grid  ! The grid to get an updater for
    logical, optional, intent(out) :: found ! Flag indicating whether the grid
                                            ! was found
    type(grid_updater_t)           :: updater

    call assert_msg( 938167254, associated( this%grid_warehouse_ ),           &
                     "Grids not available" )
    updater = this%grid_warehouse_%get_updater( grid, found )

  end function get_grid_updater

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function get_profile_updater( this, profile, found ) result( updater )
    ! Returns an updater for a profile that a host application can use to
    ! update TUV-x state data at runtime
    !
    ! If the optional `found` flag is omitted, an error is returned if the
    ! profile does not exist in TUV-x

    use musica_assert,                 only : assert_msg
    use tuvx_profile,                  only : profile_t
    use tuvx_profile_from_host,        only : profile_updater_t

    class(core_t),     intent(in)  :: this    ! TUV-x core
    class(profile_t),  intent(in)  :: profile ! The profile to get an updater for
    logical, optional, intent(out) :: found   ! Flag indicating whether the
                                              ! profile was found
    type(profile_updater_t)        :: updater

    call assert_msg( 804243032, associated( this%profile_warehouse_ ),        &
                     "Profiles not available" )
    updater = this%profile_warehouse_%get_updater( profile, found )

  end function get_profile_updater

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

    class(core_t),     intent(in)  :: this     ! TUV-x core
    class(radiator_t), intent(in)  :: radiator ! The radiator to get an updater for
    logical, optional, intent(out) :: found    ! Flag indicating whether the
                                               ! radiator was found
    type(radiator_updater_t)       :: updater

    call assert_msg( 676440920, associated( this%radiative_transfer_ ),       &
                     "Radiators not available" )
    updater = this%radiative_transfer_%get_radiator_updater( radiator, found )

  end function get_radiator_updater

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function number_of_photolysis_reactions( this )
    ! Returns the number of photolysis reactions

    class(core_t), intent(in) :: this

    number_of_photolysis_reactions = 0
    if( associated( this%photolysis_rates_ ) ) then
      number_of_photolysis_reactions = this%photolysis_rates_%size( )
    end if

  end function number_of_photolysis_reactions

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function number_of_dose_rates( this )
    ! Returns the number of dose rates

    class(core_t), intent(in) :: this

    number_of_dose_rates = 0
    if( associated( this%dose_rates_ ) ) then
      number_of_dose_rates = this%dose_rates_%size( )
    end if

  end function number_of_dose_rates

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function photolysis_reaction_labels( this ) result( labels )
    ! Returns the set of photolysis reaction labels

    class(core_t),  intent(in)  :: this
    type(string_t), allocatable :: labels(:)

    if( associated( this%photolysis_rates_ ) ) then
      labels = this%photolysis_rates_%labels( )
    else
      allocate( labels( 0 ) )
    end if

  end function photolysis_reaction_labels

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function dose_rate_labels( this ) result( labels )
    ! Returns the set of dose rate labels

    class(core_t),  intent(in)  :: this
    type(string_t), allocatable :: labels(:)

    if( associated( this%dose_rates_ ) ) then
      labels = this%dose_rates_%labels( )
    else
      allocate( labels( 0 ) )
    end if

  end function dose_rate_labels

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function get_photolysis_cross_section( this, reaction_label, found )        &
      result( cross_section )

    use musica_assert,                 only : assert_msg
    use tuvx_cross_section,            only : cross_section_t

    class(core_t),          intent(in)  :: this
    type(string_t),         intent(in)  :: reaction_label      ! Reaction to find
    logical,      optional, intent(out) :: found               ! Flag indicating if reaction exists
    real(kind=dk), allocatable          :: cross_section(:,:)  ! cross section (height,wavelength)

    logical :: l_found
    class(cross_section_t), pointer :: cs

    cs => this%photolysis_rates_%get_cross_section( reaction_label, l_found )
    if( present( found ) ) found = l_found
    call assert_msg( 135956210, present( found ) .or. l_found,                &
                     "Photolysis cross section not found for '"//             &
                     reaction_label//"'" )
    if( l_found ) then
      cross_section = cs%calculate( this%grid_warehouse_,                     &
                                    this%profile_warehouse_ )
      deallocate( cs )
    end if

  end function get_photolysis_cross_section

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function get_photolysis_quantum_yield( this, reaction_label, found )        &
      result( quantum_yield )

    use musica_assert,                 only : assert_msg
    use tuvx_quantum_yield,            only : quantum_yield_t

    class(core_t),          intent(in)  :: this
    type(string_t),         intent(in)  :: reaction_label      ! Reaction to find
    logical,      optional, intent(out) :: found               ! Flag indicating if reaction exists
    real(kind=dk), allocatable          :: quantum_yield(:,:)  ! quantum yield (height,wavelength)

    logical :: l_found
    class(quantum_yield_t), pointer :: qy

    qy => this%photolysis_rates_%get_quantum_yield( reaction_label, l_found )
    if( present( found ) ) found = l_found
    call assert_msg( 454211846, present( found ) .or. l_found,                &
                     "Photolysis quantum yield not found for '"//             &
                     reaction_label//"'" )
    if( l_found ) then
      quantum_yield = qy%calculate( this%grid_warehouse_,                     &
                                    this%profile_warehouse_ )
      deallocate( qy )
    end if

  end function get_photolysis_quantum_yield

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  type(radiation_field_t) function get_radiation_field( this ) result( field )

    class(core_t), intent(in) :: this

    field = this%radiation_field_

  end function get_radiation_field

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function pack_size( this, comm )
    ! Returns the number of bytes required to pack the core onto a buffer

    use musica_mpi,                    only : musica_mpi_pack_size

    class(core_t),     intent(in) :: this ! core to be packed
    integer,           intent(in) :: comm ! MPI communicator

#ifdef MUSICA_USE_MPI
    pack_size =                                                               &
        musica_mpi_pack_size( associated( this%grid_warehouse_ ), comm )
    if( associated( this%grid_warehouse_ ) ) then
      pack_size = pack_size + this%grid_warehouse_%pack_size( comm )
    end if
    pack_size = pack_size +                                                   &
        musica_mpi_pack_size( associated( this%profile_warehouse_ ), comm )
    if( associated( this%profile_warehouse_ ) ) then
      pack_size = pack_size + this%profile_warehouse_%pack_size( comm )
    end if
    pack_size = pack_size +                                                   &
        musica_mpi_pack_size( associated( this%spherical_geometry_ ), comm )
    if( associated( this%spherical_geometry_ ) ) then
      pack_size = pack_size + this%spherical_geometry_%pack_size( comm )
    end if
    pack_size = pack_size +                                                   &
        musica_mpi_pack_size( associated( this%la_sr_bands_ ), comm )
    if( associated( this%la_sr_bands_ ) ) then
      pack_size = pack_size + this%la_sr_bands_%pack_size( comm )
    end if
    pack_size = pack_size +                                                   &
        musica_mpi_pack_size( this%enable_diagnostics_ , comm )
    pack_size = pack_size +                                                   &
        musica_mpi_pack_size( associated( this%radiative_transfer_ ), comm )
    if( associated( this%radiative_transfer_ ) ) then
      pack_size = pack_size + this%radiative_transfer_%pack_size( comm )
    end if
    pack_size = pack_size +                                                   &
        musica_mpi_pack_size( associated( this%photolysis_rates_ ), comm )
    if( associated( this%photolysis_rates_ ) ) then
      pack_size = pack_size + this%photolysis_rates_%pack_size( comm )
    end if
    pack_size = pack_size +                                                   &
      musica_mpi_pack_size( associated( this%dose_rates_ ), comm )
    if( associated( this%dose_rates_ ) ) then
      pack_size = pack_size + this%dose_rates_%pack_size( comm )
    end if
#else
    pack_size = 0
#endif

  end function pack_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine mpi_pack( this, buffer, position, comm )
    ! Packs the core onto a character buffer

    use musica_assert,                 only : assert
    use musica_mpi,                    only : musica_mpi_pack

    class(core_t),     intent(in)    :: this      ! core to be packed
    character,         intent(inout) :: buffer(:) ! memory buffer
    integer,           intent(inout) :: position  ! current buffer position
    integer,           intent(in)    :: comm      ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: prev_pos

    prev_pos = position
    call musica_mpi_pack( buffer, position,                                   &
                          associated( this%grid_warehouse_ ), comm )
    if( associated( this%grid_warehouse_ ) ) then
      call this%grid_warehouse_%mpi_pack( buffer, position, comm )
    end if
    call musica_mpi_pack( buffer, position,                                   &
                          associated( this%profile_warehouse_ ), comm )
    if( associated( this%profile_warehouse_ ) ) then
      call this%profile_warehouse_%mpi_pack( buffer, position, comm )
    end if
    call musica_mpi_pack( buffer, position,                                   &
                          associated( this%spherical_geometry_ ), comm )
    if( associated( this%spherical_geometry_ ) ) then
      call this%spherical_geometry_%mpi_pack( buffer, position, comm )
    end if
    call musica_mpi_pack( buffer, position,                                   &
                          associated( this%la_sr_bands_ ), comm )
    if( associated( this%la_sr_bands_ ) ) then
      call this%la_sr_bands_%mpi_pack( buffer, position, comm )
    end if
    call musica_mpi_pack( buffer, position, this%enable_diagnostics_ , comm )
    call musica_mpi_pack( buffer, position,                                   &
                          associated( this%radiative_transfer_ ), comm )
    if( associated( this%radiative_transfer_ ) ) then
      call this%radiative_transfer_%mpi_pack( buffer, position, comm )
    end if
    call musica_mpi_pack( buffer, position,                                   &
                          associated( this%photolysis_rates_ ), comm )
    if( associated( this%photolysis_rates_ ) ) then
      call this%photolysis_rates_%mpi_pack( buffer, position, comm )
    end if
    call musica_mpi_pack( buffer, position,                                   &
                          associated( this%dose_rates_ ), comm )
    if( associated( this%dose_rates_ ) ) then
      call this%dose_rates_%mpi_pack( buffer, position, comm )
    end if
    call assert( 332208077, position - prev_pos <= this%pack_size( comm ) )
#endif

  end subroutine mpi_pack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine mpi_unpack( this, buffer, position, comm )
    ! Unpacks a core from a character buffer

    use musica_assert,                 only : assert
    use musica_mpi,                    only : musica_mpi_unpack

    class(core_t),     intent(out)   :: this      ! core to be unpacked
    character,         intent(inout) :: buffer(:) ! memory buffer
    integer,           intent(inout) :: position  ! current buffer position
    integer,           intent(in)    :: comm      ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: prev_pos
    logical :: alloced

    prev_pos = position
    call musica_mpi_unpack( buffer, position, alloced, comm )
    if( alloced ) then
      allocate( this%grid_warehouse_ )
      call this%grid_warehouse_%mpi_unpack( buffer, position, comm )
    end if
    call musica_mpi_unpack( buffer, position, alloced, comm )
    if( alloced ) then
      allocate( this%profile_warehouse_ )
      call this%profile_warehouse_%mpi_unpack( buffer, position, comm )
    end if
    call musica_mpi_unpack( buffer, position, alloced, comm )
    if( alloced ) then
      allocate( this%spherical_geometry_ )
      call this%spherical_geometry_%mpi_unpack( buffer, position, comm )
    end if
    call musica_mpi_unpack( buffer, position, alloced, comm )
    if( alloced ) then
      allocate( this%la_sr_bands_ )
      call this%la_sr_bands_%mpi_unpack( buffer, position, comm )
    end if
    call musica_mpi_unpack( buffer, position, this%enable_diagnostics_, comm )
    call musica_mpi_unpack( buffer, position, alloced, comm )
    if( alloced ) then
      allocate( this%radiative_transfer_ )
      call this%radiative_transfer_%mpi_unpack( buffer, position, comm )
    end if
    call musica_mpi_unpack( buffer, position, alloced, comm )
    if( alloced ) then
      allocate( this%photolysis_rates_ )
      call this%photolysis_rates_%mpi_unpack( buffer, position, comm )
    end if
    call musica_mpi_unpack( buffer, position, alloced, comm )
    if( alloced ) then
      allocate( this%dose_rates_ )
      call this%dose_rates_%mpi_unpack( buffer, position, comm )
    end if
    call assert( 332208077, position - prev_pos <= this%pack_size( comm ) )
#endif

  end subroutine mpi_unpack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine finalize( this )
    ! Finalizes the core

    !> Photolysis core
    type(core_t), intent(inout) :: this

    if( associated( this%grid_warehouse_ ) ) then
      deallocate( this%grid_warehouse_ )
    end if
    if( associated( this%profile_warehouse_ ) ) then
      deallocate( this%profile_warehouse_ )
    end if
    if( associated( this%spherical_geometry_ ) ) then
      deallocate( this%spherical_geometry_ )
    end if
    if( associated( this%la_sr_bands_ ) ) then
      deallocate( this%la_sr_bands_ )
    end if
    if( associated( this%radiative_transfer_ ) ) then
      deallocate( this%radiative_transfer_ )
    end if
    if( associated( this%photolysis_rates_ ) ) then
      deallocate( this%photolysis_rates_ )
    end if
    if( associated( this%dose_rates_ ) ) then
      deallocate( this%dose_rates_ )
    end if
    if( associated( this%radiation_field_ ) ) then
      deallocate( this%radiation_field_ )
    end if

  end subroutine finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_core
