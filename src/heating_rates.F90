! Copyright (C) 2024-2026 University Corporation for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
module tuvx_heating_rates
  ! The chemical potential heating rates type heating_rates_t and related functions

  use musica_assert,                 only : assert, assert_msg
  use musica_config,                 only : config_t
  use musica_constants,              only : dk => musica_dk
  use musica_iterator,               only : iterator_t
  use musica_mpi,                    only : musica_mpi_pack, musica_mpi_pack_size, musica_mpi_unpack
  use musica_string,                 only : string_t
  use tuvx_constants,                only : hc
  use tuvx_cross_section,            only : cross_section_ptr
  use tuvx_cross_section_factory,    only : cross_section_allocate, cross_section_builder, cross_section_type_name
  use tuvx_grid,                     only : grid_t
  use tuvx_grid_warehouse,           only : grid_warehouse_ptr, grid_warehouse_t
  use tuvx_la_sr_bands,              only : la_sr_bands_t
  use tuvx_profile,                  only : profile_t
  use tuvx_profile_warehouse,        only : profile_warehouse_ptr, profile_warehouse_t
  use tuvx_quantum_yield,            only : quantum_yield_ptr
  use tuvx_quantum_yield_factory,    only : quantum_yield_allocate, quantum_yield_builder, quantum_yield_type_name
  use tuvx_solver,                   only : radiation_field_t
  use tuvx_spherical_geometry,       only : spherical_geometry_t

  use tuvx_prepared_spectra, only : spectrum_cache_t, spectrum_plan_t, rate_plan_t, prepare_rate, integrate_rate

  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  implicit none

  private
  public :: heating_rates_t

  type :: heating_parameters_t
    ! Heating parameters for a single photolyzing species
    type(string_t)             :: label_          ! label for the heating rate
    type(cross_section_ptr)    :: cross_section_  ! cross section
    type(quantum_yield_ptr)    :: quantum_yield_  ! quantum yield
    real(kind=dk)              :: scaling_factor_ ! scaling factor for the heating rate
    real(kind=dk), allocatable :: energy_(:)      ! wavelength resolved bond-dissociation energy [J]
  contains
    !> Returns the size of a character buffer needed to pack the heating parameters
    procedure :: pack_size => heating_parameters_pack_size
    !> Packs the heating parameters into a character buffer
    procedure :: mpi_pack => heating_parameters_mpi_pack
    !> Unpacks the heating parameters from a character buffer
    procedure :: mpi_unpack => heating_parameters_mpi_unpack
  end type heating_parameters_t

  !> heating_parameters_t constructor
  interface heating_parameters_t
    module procedure :: heating_parameters_constructor
  end interface heating_parameters_t

  type :: heating_rates_t
    type(heating_parameters_t), allocatable :: heating_parameters_(:) ! heating parameters for each photolyzing species
    type(grid_warehouse_ptr) :: height_grid_     ! height grid
    type(grid_warehouse_ptr) :: wavelength_grid_ ! wavelength grid
    type(profile_warehouse_ptr) :: etfl_profile_ ! Extraterrestrial flux profile
    type(profile_warehouse_ptr) :: air_profile_  ! Air profile
    integer, allocatable :: o2_rate_indices_(:)  ! indices in the heating rates array where O2
                                                 ! corrections to the cross-section in the
                                                 ! Lyman-Alpha and Schumann-Runge bands should
                                                 ! be applied
    type(spectrum_cache_t), pointer :: spectrum_cache_ => null()
    logical :: own_spectrum_cache_ = .false.
    type(rate_plan_t), allocatable :: prepared_rates_(:)
  contains
    procedure :: set_spectrum_cache
    !> Calulates the heating rates
    procedure :: get
    !> Returns the names of each photolysis reaction with a heating rate
    procedure :: labels
    !> Returns the number of heating rates
    procedure :: size => get_number
    !> Returns the size of a character buffer needed to pack the heating rates
    procedure :: pack_size
    !> Packs the heating rates into a character buffer
    procedure :: mpi_pack
    !> Unpacks the heating rates from a character buffer
    procedure :: mpi_unpack
    !> Cleans up memory
    final :: destructor
  end type heating_rates_t

  !> heating_rates_t constructor
  interface heating_rates_t
    module procedure :: constructor
  end interface heating_rates_t

contains

  subroutine set_spectrum_cache(this,cache)
    class(heating_rates_t), intent(inout) :: this
    type(spectrum_cache_t), pointer, intent(in) :: cache
    if( this%own_spectrum_cache_ .and. associated(this%spectrum_cache_) ) deallocate(this%spectrum_cache_)
    this%spectrum_cache_ => cache
    this%own_spectrum_cache_ = .false.
    if( allocated(this%prepared_rates_) ) deallocate(this%prepared_rates_)
  end subroutine set_spectrum_cache
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> heating_rates_t constructor
  function constructor( config, grids, profiles ) result( this )


    !> Heating rate collection
    type(heating_rates_t),     pointer       :: this
    !> Configuration
    type(config_t),            intent(inout) :: config
    !> Grids
    type(grid_warehouse_t),    intent(inout) :: grids
    !> Profiles
    type(profile_warehouse_t), intent(inout) :: profiles

    character(len=*), parameter :: Iam = 'heating rates constructor'
    type(config_t) :: reaction_set, reaction_config, heating_config
    type(config_t) :: cross_section_config
    class(iterator_t), pointer :: iter
    type(string_t) :: required_keys(1), optional_keys(1)
    logical :: found, do_apply_bands
    integer :: n_hr, i_hr, n_O2, i_O2

    required_keys(1) = "reactions"
    optional_keys(1) = "enable diagnostics"

    call assert_msg( 310567326,                                               &
                     config%validate( required_keys, optional_keys ),         &
                      "Invalid configuration for heating rates" )

    allocate( this )
    this%height_grid_ = grids%get_ptr( "height", "km" )
    this%wavelength_grid_ = grids%get_ptr( "wavelength", "nm" )
    this%etfl_profile_ = profiles%get_ptr( "extraterrestrial flux",           &
                                           "photon cm-2 s-1" )
    this%air_profile_ = profiles%get_ptr( "air", "molecule cm-3" )

    ! iterate over photolysis reactions looking for those with
    ! heating rate parameters
    call config%get( "reactions", reaction_set, Iam )
    iter => reaction_set%get_iterator( )
    n_hr = 0
    n_O2 = 0
    do while( iter%next( ) )
      call reaction_set%get( iter, reaction_config, Iam )
      call reaction_config%get( "heating", heating_config, Iam, found = found )
      if( found ) then
        n_hr = n_hr + 1
        call reaction_config%get( "cross section", cross_section_config, Iam )
        call cross_section_config%get( "apply O2 bands", do_apply_bands, Iam, &
                                  default = .false. )
        if( do_apply_bands ) n_O2 = n_O2 + 1
      end if
    end do
    allocate( this%heating_parameters_( n_hr ) )
    allocate( this%o2_rate_indices_( n_O2 ) )
    call iter%reset( )
    i_hr = 0
    i_O2 = 0
    do while( iter%next( ) )
      call reaction_set%get( iter, reaction_config, Iam )
      call reaction_config%get( "heating", heating_config, Iam, found = found )
      if( found ) then
        i_hr = i_hr + 1
        this%heating_parameters_( i_hr ) =                                    &
          heating_parameters_constructor( reaction_config, grids, profiles )
        call reaction_config%get( "cross section", cross_section_config, Iam )
        call cross_section_config%get( "apply O2 bands", do_apply_bands, Iam, &
                                  default = .false. )
        if( do_apply_bands ) then
          i_O2 = i_O2 + 1
          this%o2_rate_indices_( i_O2 ) = i_hr
        end if
      end if
    end do
    call assert( 357615745, i_hr .eq. n_hr )
    call assert( 336635308, i_O2 .eq. n_O2 )
    deallocate( iter )

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> heating_parameters_t constructor
  function heating_parameters_constructor( config, grids, profiles )          &
      result( this )


    !> Heating parameters for a single photolyzing species
    type(heating_parameters_t)               :: this
    !> Configuration for the photolysis reaction
    type(config_t),            intent(inout) :: config
    !> Grids
    type(grid_warehouse_t),    intent(inout) :: grids
    !> Profiles
    type(profile_warehouse_t), intent(inout) :: profiles

    character(len=*), parameter :: Iam = 'heating parameters constructor'
    type(config_t) :: heating_config, cs_config, qy_config
    class(grid_t), pointer :: wavelengths
    real(kind=dk) :: energy_term
    type(string_t) :: required_keys(4), optional_keys(1)
    type(string_t) :: heating_required_keys(1)
    type(string_t), allocatable :: heating_optional_keys(:)

    allocate(heating_optional_keys(0))

    required_keys(1) = "name"
    required_keys(2) = "cross section"
    required_keys(3) = "quantum yield"
    required_keys(4) = "heating"
    optional_keys(1) = "scaling factor"

    call assert_msg( 316144353,                                               &
                     config%validate( required_keys, optional_keys ),         &
                     "Invalid configuration for photolysis reactions with "// &
                     "heating parameters" )

    call config%get( "heating", heating_config, Iam )

    heating_required_keys(1) = "energy term"

    call assert_msg( 316144354,                                               &
                     heating_config%validate( heating_required_keys,          &
                                              heating_optional_keys ), &
                     "Invalid configuration for heating parameters" )

    call config%get( "name", this%label_, Iam )
    call config%get( "cross section", cs_config, Iam )
    this%cross_section_%val_ => cross_section_builder( cs_config, grids,      &
                                                       profiles )
    call config%get( "quantum yield", qy_config, Iam )
    this%quantum_yield_%val_ => quantum_yield_builder( qy_config, grids,      &
                                                       profiles )
    call config%get( "scaling factor", this%scaling_factor_, Iam,             &
                     default = 1.0_dk )
    call heating_config%get( "energy term", energy_term, Iam )
    wavelengths => grids%get_grid( "wavelength", "nm" )
    allocate( this%energy_( wavelengths%ncells_ ) )
    this%energy_(:) =                                                         &
        max( 0.0_dk, hc * 1.0e9_dk * ( energy_term - wavelengths%mid_(:) ) /  &
                                     ( energy_term * wavelengths%mid_(:) ) )
    deallocate( wavelengths )

  end function heating_parameters_constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate heating rates
  subroutine get( this, la_srb, spherical_geometry, grids, profiles,          &
                  radiation_field, heating_rates )


    !> Heating rate collection
    class(heating_rates_t),      intent(inout) :: this
    !> Lyman Alpha and Schumann-Runge bands
    class(la_sr_bands_t),        intent(inout) :: la_srb
    !> Spherical geometry
    class(spherical_geometry_t), intent(inout) :: spherical_geometry
    !> Grids
    class(grid_warehouse_t),     intent(inout) :: grids
    !> Profiles
    class(profile_warehouse_t),  intent(inout) :: profiles
    !> Radiation field
    class(radiation_field_t),    intent(in)    :: radiation_field
    !> Heating rates (vertical interface, reaction) [J s-1]
    real(kind=dk),               intent(inout) :: heating_rates(:,:)

    character(len=*), parameter :: Iam = 'heating rates get'
    class(grid_t), pointer :: heights, wavelengths
    class(profile_t), pointer :: etfl, air
    real(kind=dk), allocatable :: actinic_flux(:,:)
    real(kind=dk), allocatable :: cross_section(:,:), quantum_yield(:,:)
    real(kind=dk), allocatable :: air_vertical_column(:), air_slant_column(:)
    integer :: i_rate, n_rates, i_height, i_wavelength

    type(spectrum_plan_t), pointer :: cs_plan, qy_plan
    logical :: prepared, oxygen
    logical, allocatable :: oxygen_bands(:)

    if( .not. associated(this%spectrum_cache_) ) then
      allocate(this%spectrum_cache_)
      this%own_spectrum_cache_ = .true.
    end if
    if( this%own_spectrum_cache_ ) call this%spectrum_cache_%begin_column(grids,profiles)
    heights => grids%get_grid( this%height_grid_ )
    wavelengths => grids%get_grid( this%wavelength_grid_ )
    etfl => profiles%get_profile( this%etfl_profile_ )
    air => profiles%get_profile( this%air_profile_ )

    n_rates = size( this%heating_parameters_ )
    call assert( 966855732,                                                   &
                 size( heating_rates, 1 ) .eq. heights%ncells_ + 1 .and.      &
                 size( heating_rates, 2 ) .eq. n_rates )
    if( n_rates == 0 ) then
      deallocate( heights )
      deallocate( wavelengths )
      deallocate( etfl )
      deallocate( air )
      return
    end if

    actinic_flux = radiation_field%fdr_ + radiation_field%fup_ + radiation_field%fdn_
    do i_wavelength = 1, size( actinic_flux, 2 )
      actinic_flux(:,i_wavelength) = actinic_flux(:,i_wavelength) * etfl%mid_val_(i_wavelength)
    end do
    where( actinic_flux < 0.0_dk )
      actinic_flux = 0.0_dk
    end where

    allocate(oxygen_bands(wavelengths%ncells_))
    oxygen_bands = .false.
    if( size(this%o2_rate_indices_) > 0 ) then
      if( la_srb%has_la ) oxygen_bands(la_srb%ila) = .true.
      if( la_srb%has_srb ) oxygen_bands(la_srb%isrb:la_srb%isrb+16) = .true.
    end if
    prepared = this%spectrum_cache_%usable .and. all(ieee_is_finite(actinic_flux))
    if( allocated(this%prepared_rates_) ) then
      if( size(this%prepared_rates_) /= n_rates ) deallocate(this%prepared_rates_)
    end if
    if( .not. allocated(this%prepared_rates_) ) allocate(this%prepared_rates_(n_rates))
rate_loop: do i_rate = 1, n_rates
    associate( params => this%heating_parameters_( i_rate ) )
      if( prepared ) then
        cs_plan => this%spectrum_cache_%cross_section(params%cross_section_%val_,grids,profiles)
        qy_plan => this%spectrum_cache_%quantum_yield(params%quantum_yield_%val_,grids,profiles)
        oxygen = any(this%o2_rate_indices_ == i_rate)
        associate(r => this%prepared_rates_(i_rate))
          if( r%generation /= this%spectrum_cache_%generation ) then
            call prepare_rate(r,cs_plan,qy_plan,this%spectrum_cache_%generation,oxygen,params%energy_,oxygen_bands)
          else if( any(r%energy /= params%energy_) ) then
            call prepare_rate(r,cs_plan,qy_plan,this%spectrum_cache_%generation,oxygen,params%energy_,oxygen_bands)
          end if
          if( oxygen ) then
            if( cs_plan%oxygen_epoch /= this%spectrum_cache_%epoch ) then
              allocate(air_vertical_column(air%ncells_),air_slant_column(air%ncells_+1))
              call spherical_geometry%air_mass(air%exo_layer_dens_,air_vertical_column,air_slant_column)
              cs_plan%oxygen = cs_plan%materialize()
              call la_srb%cross_section(grids,profiles,air_vertical_column,air_slant_column, &
                                        cs_plan%oxygen,spherical_geometry)
              cs_plan%oxygen_epoch = this%spectrum_cache_%epoch
              deallocate(air_vertical_column,air_slant_column)
            end if
            call integrate_rate(r,actinic_flux,params%scaling_factor_,heating_rates(:,i_rate),cs_plan%oxygen)
          else
            call integrate_rate(r,actinic_flux,params%scaling_factor_,heating_rates(:,i_rate))
          end if
        end associate
        cycle rate_loop
      end if
      cross_section = params%cross_section_%val_%calculate_column( grids, profiles )
      quantum_yield = params%quantum_yield_%val_%calculate_column( grids, profiles )

      ! O2 photolysis can have special la & srb band handling
      if( any( this%o2_rate_indices_ == i_rate ) ) then
        if( size( cross_section, 1 ) == 1 ) &
          cross_section = spread( cross_section(1,:), 1, heights%ncells_ + 1 )
        allocate( air_vertical_column( air%ncells_ ),                         &
                  air_slant_column( air%ncells_ + 1 ) )
        call spherical_geometry%air_mass( air%exo_layer_dens_,                &
                                          air_vertical_column,                &
                                          air_slant_column )
        call la_srb%cross_section( grids, profiles, air_vertical_column,      &
                                   air_slant_column, cross_section,           &
                                   spherical_geometry )
        deallocate( air_vertical_column, air_slant_column )
      end if

      heating_rates(:,i_rate) = 0.0_dk
      do i_wavelength = 1, size( actinic_flux, 2 )
        do i_height = 1, heights%ncells_ + 1
          heating_rates(i_height,i_rate) = heating_rates(i_height,i_rate) + &
              actinic_flux(i_height,i_wavelength) * &
              ( params%energy_(i_wavelength) * &
                ( cross_section(min(i_height,size(cross_section,1)),i_wavelength) * &
                  quantum_yield(min(i_height,size(quantum_yield,1)),i_wavelength) ) )
        end do
      end do
      heating_rates(:,i_rate) = heating_rates(:,i_rate) * params%scaling_factor_
    end associate
    end do rate_loop

    deallocate( heights )
    deallocate( wavelengths )
    deallocate( etfl )
    deallocate( air )


  end subroutine get

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the names of each photolysis reaction with a heating rate
  function labels( this ) result(res)

    !> Heating rate collection
    class(heating_rates_t), intent(in) :: this
    !> Photolysis reaction labels
    type(string_t), allocatable :: res(:)
    integer :: i

    allocate( res( size( this%heating_parameters_ ) ) )

    do i = 1, size(res)
      res(i) = this%heating_parameters_(i)%label_
    end do

  end function labels

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the number of heating rates
  function get_number( this ) result( n_rates )

    !> Number of heating rates
    integer :: n_rates
    !> Heating rate collection
    class(heating_rates_t), intent(in) :: this

    n_rates = 0
    if( allocated( this%heating_parameters_ ) ) then
      n_rates = size( this%heating_parameters_ )
    end if

  end function get_number

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the size of a character buffer needed to pack the heating rates
  function pack_size( this, comm )


    !> Heating rate collection
    class(heating_rates_t), intent(in) :: this
    !> MPI communicator
    integer,                intent(in) :: comm
    !> Size of the character buffer
    integer                            :: pack_size

#ifdef MUSICA_USE_MPI
    integer :: i_elem

    pack_size = musica_mpi_pack_size( allocated( this%heating_parameters_ ),  &
                                      comm )
    if( allocated( this%heating_parameters_ ) ) then
      pack_size = pack_size +                                                 &
          musica_mpi_pack_size( size( this%heating_parameters_ ), comm )
      do i_elem = 1, size( this%heating_parameters_ )
        pack_size = pack_size +                                               &
            this%heating_parameters_( i_elem )%pack_size( comm )
      end do
    end if
    pack_size = pack_size +                                                   &
                this%height_grid_%pack_size( comm ) +                         &
                this%wavelength_grid_%pack_size( comm ) +                     &
                this%etfl_profile_%pack_size( comm ) +                        &
                this%air_profile_%pack_size( comm ) +                         &
                musica_mpi_pack_size( this%o2_rate_indices_, comm )
#else
    pack_size = 0
#endif

  end function pack_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Packs the heating rates into a character buffer
  subroutine mpi_pack( this, buffer, position, comm )

    !> Heating rate collection
    class(heating_rates_t), intent(in)    :: this
    !> Character buffer
    character,              intent(inout) :: buffer(:)
    !> Position in the buffer
    integer,                intent(inout) :: position
    !> MPI communicator
    integer,                intent(in) :: comm

#ifdef MUSICA_USE_MPI
    integer :: prev_pos, i_elem

    prev_pos = position
    call musica_mpi_pack( buffer, position,                                   &
                          allocated( this%heating_parameters_ ), comm )
    if( allocated( this%heating_parameters_ ) ) then
      call musica_mpi_pack( buffer, position,                                 &
                            size( this%heating_parameters_ ), comm )
      do i_elem = 1, size( this%heating_parameters_ )
        call this%heating_parameters_( i_elem )%mpi_pack( buffer, position,   &
                                                          comm )
      end do
    end if
    call this%height_grid_%mpi_pack( buffer, position, comm )
    call this%wavelength_grid_%mpi_pack( buffer, position, comm )
    call this%etfl_profile_%mpi_pack( buffer, position, comm )
    call this%air_profile_%mpi_pack( buffer, position, comm )
    call musica_mpi_pack( buffer, position, this%o2_rate_indices_, comm )
    call assert( 247051769, position - prev_pos <= this%pack_size( comm ) )
#endif

  end subroutine mpi_pack
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Unpacks the heating rates from a character buffer
  subroutine mpi_unpack( this, buffer, position, comm )

    !> Heating rate collection
    class(heating_rates_t), intent(out) :: this
    !> Character buffer
    character,              intent(inout) :: buffer(:)
    !> Position in the buffer
    integer,                intent(inout) :: position
    !> MPI communicator
    integer,                intent(in)    :: comm

#ifdef MUSICA_USE_MPI
    integer :: prev_pos, i_elem, n_elems
    logical :: is_allocated

    prev_pos = position
    call musica_mpi_unpack( buffer, position, is_allocated, comm )
    if( is_allocated ) then
      call musica_mpi_unpack( buffer, position, n_elems, comm )
      allocate( this%heating_parameters_( n_elems ) )
      do i_elem = 1, n_elems
        call this%heating_parameters_( i_elem )%mpi_unpack( buffer, position, &
                                                            comm )
      end do
    end if
    call this%height_grid_%mpi_unpack( buffer, position, comm )
    call this%wavelength_grid_%mpi_unpack( buffer, position, comm )
    call this%etfl_profile_%mpi_unpack( buffer, position, comm )
    call this%air_profile_%mpi_unpack( buffer, position, comm )
    call musica_mpi_unpack( buffer, position, this%o2_rate_indices_, comm )
    call assert( 631316749, position - prev_pos <= this%pack_size( comm ) )
#endif

  end subroutine mpi_unpack
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the size of a character buffer needed to pack the heating
  !! parameters
  function heating_parameters_pack_size( this, comm ) result( pack_size )

    !> Heating parameters for a single photolyzing species
    class(heating_parameters_t), intent(in) :: this
    !> MPI communicator
    integer,                     intent(in) :: comm
    !> Size of the character buffer
    integer                                 :: pack_size

#ifdef MUSICA_USE_MPI
    type(string_t) :: cs_type_name, qy_type_name

    cs_type_name = cross_section_type_name( this%cross_section_%val_ )
    qy_type_name = quantum_yield_type_name( this%quantum_yield_%val_ )
    pack_size = this%label_%pack_size(  comm ) +                              &
                cs_type_name%pack_size( comm ) +                              &
                this%cross_section_%val_%pack_size( comm ) +                  &
                qy_type_name%pack_size( comm ) +                              &
                this%quantum_yield_%val_%pack_size( comm ) +                  &
                musica_mpi_pack_size( this%scaling_factor_, comm ) +          &
                musica_mpi_pack_size( this%energy_, comm )
#else
    pack_size = 0
#endif

  end function heating_parameters_pack_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Packs the heating parameters into a character buffer
  subroutine heating_parameters_mpi_pack( this, buffer, position, comm )

    !> Heating parameters for a single photolyzing species
    class(heating_parameters_t), intent(in)    :: this
    !> Character buffer
    character,                   intent(inout) :: buffer(:)
    !> Position in the buffer
    integer,                     intent(inout) :: position
    !> MPI communicator
    integer,                     intent(in)    :: comm

#ifdef MUSICA_USE_MPI
    integer :: prev_pos
    type(string_t) :: cs_type_name, qy_type_name

    prev_pos = position
    cs_type_name = cross_section_type_name( this%cross_section_%val_ )
    qy_type_name = quantum_yield_type_name( this%quantum_yield_%val_ )
    call this%label_%mpi_pack( buffer, position, comm )
    call cs_type_name%mpi_pack( buffer, position, comm )
    call this%cross_section_%val_%mpi_pack( buffer, position, comm )
    call qy_type_name%mpi_pack( buffer, position, comm )
    call this%quantum_yield_%val_%mpi_pack( buffer, position, comm )
    call musica_mpi_pack( buffer, position, this%scaling_factor_, comm )
    call musica_mpi_pack( buffer, position, this%energy_, comm )
    call assert( 243240701, position - prev_pos <= this%pack_size( comm ) )
#endif

  end subroutine heating_parameters_mpi_pack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Unpacks the heating parameters from a character buffer
  subroutine heating_parameters_mpi_unpack( this, buffer, position, comm )

    !> Heating parameters for a single photolyzing species
    class(heating_parameters_t), intent(out) :: this
    !> Character buffer
    character,                   intent(inout) :: buffer(:)
    !> Position in the buffer
    integer,                     intent(inout) :: position
    !> MPI communicator
    integer,                     intent(in)    :: comm

#ifdef MUSICA_USE_MPI
    integer :: prev_pos
    type(string_t) :: cs_type_name, qy_type_name

    prev_pos = position
    call this%label_%mpi_unpack( buffer, position, comm )
    call cs_type_name%mpi_unpack( buffer, position, comm )
    this%cross_section_%val_ => cross_section_allocate( cs_type_name )
    call this%cross_section_%val_%mpi_unpack( buffer, position, comm )
    call qy_type_name%mpi_unpack( buffer, position, comm )
    this%quantum_yield_%val_ => quantum_yield_allocate( qy_type_name )
    call this%quantum_yield_%val_%mpi_unpack( buffer, position, comm )
    call musica_mpi_unpack( buffer, position, this%scaling_factor_, comm )
    call musica_mpi_unpack( buffer, position, this%energy_, comm )
    call assert( 243240702, position - prev_pos <= this%pack_size( comm ) )
#endif

  end subroutine heating_parameters_mpi_unpack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Cleans up memory
  subroutine destructor( this )

    !> Heating rates
    type(heating_rates_t), intent(inout) :: this

    integer :: i_rate

    if( this%own_spectrum_cache_ .and. associated(this%spectrum_cache_) ) deallocate(this%spectrum_cache_)

    if( allocated( this%heating_parameters_ ) ) then
      do i_rate = 1, size( this%heating_parameters_ )
      associate( params => this%heating_parameters_( i_rate ) )
        if( associated( params%cross_section_%val_ ) ) then
          deallocate( params%cross_section_%val_ )
          nullify( params%cross_section_%val_ )
        end if
        if( associated( params%quantum_yield_%val_ ) ) then
          deallocate( params%quantum_yield_%val_ )
          nullify( params%quantum_yield_%val_ )
        end if
      end associate
      end do
      deallocate( this%heating_parameters_ )
    end if

  end subroutine destructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_heating_rates
