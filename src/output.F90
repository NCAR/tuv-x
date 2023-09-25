! Copyright (C) 2023 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
module tuvx_output
  ! Writer of TUV-x data

  use musica_io,                       only : io_t
  use musica_string,                   only : string_t
  use tuvx_grid_warehouse,             only : grid_warehouse_ptr
  use tuvx_profile_warehouse,          only : profile_warehouse_ptr

  implicit none
  private

  public :: output_t

  !> Writer class for TUV-x data
  !!
  !! Instances of \c output_t can be used by host applications to write TUV-x
  !! data and diagnostics to a file.
  type :: output_t
    private
    class(io_t), pointer        :: file_ => null( )
    logical                     :: do_photo_ = .false.
    logical                     :: do_dose_ = .false.
    logical                     :: do_radiation_ = .false.
    type(string_t), allocatable :: photo_labels_(:)
    type(string_t), allocatable :: dose_labels_(:)
    type(string_t), allocatable :: photo_cross_sections_(:)
    type(string_t), allocatable :: photo_quantum_yields_(:)
  contains
    procedure :: output
    procedure, private :: add_grids
    procedure, private :: add_photolysis_diagnostics
    final :: finalize
  end type output_t

  interface output_t
    module procedure :: constructor
  end interface output_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Creates an output_t object for a given set of variables
  function constructor( config, core ) result( this )

    use musica_assert,                 only : assert_msg
    use musica_config,                 only : config_t
    use musica_io_netcdf,              only : io_netcdf_t
    use tuvx_core,                     only : core_t

    type(output_t), pointer       :: this
    type(config_t), intent(inout) :: config
    class(core_t),  intent(inout) :: core

    character(len=*), parameter :: Iam = "output writer"
    integer        :: stat
    type(string_t) :: file_path
    type(string_t) :: required_keys(2), optional_keys(3)
    type(config_t) :: tuvx_config, rad_config

    required_keys(1) = "file path"
    required_keys(2) = "tuv-x configuration"
    optional_keys(1) = "include photolysis"
    optional_keys(2) = "include dose rates"

    call assert_msg( 215370625,                                               &
                     config%validate( required_keys, optional_keys ),         &
                     "Bad configuration for TUV-x output writer" )

    allocate( this )

    ! Get the file path and overwrite any existing file with the same name
    ! NOTE: Could add option for other output file types
    call config%get( "file path", file_path, Iam )
    open( unit = 16, iostat = stat, file = file_path%to_char( ),              &
          status = 'old' )
    if( stat == 0 ) close( 16, status = 'delete' )
    this%file_ => io_netcdf_t( file_path )

    ! Add all grids as file dimensions
    call this%add_grids( core )

    ! Add photolysis and dose rates
    call config%get( "include photolysis", this%do_photo_, Iam,               &
                     default = .false. )
    if( this%do_photo_ )                                                      &
        this%photo_labels_ = core%photolysis_reaction_labels( )

    call config%get( "include dose rates", this%do_dose_, Iam,                &
                     default = .false. )
    if( this%do_dose_ ) this%dose_labels_ = core%dose_rate_labels( )

    ! Add custom diagnostics
    call config%get( "tuv-x configuration", tuvx_config, Iam )
    call this%add_photolysis_diagnostics( tuvx_config )
    call tuvx_config%get( "radiative transfer", rad_config, Iam )
    call rad_config%get( "__output", this%do_radiation_, Iam,                 &
                         default = .false. )

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Outputs results
  subroutine output( this, step, core, photolysis_rate_constants, dose_rates, &
     time, solar_zenith_angle, earth_sun_distance )

    use musica_assert,                 only : assert_msg
    use musica_constants,              only : dk => musica_dk
    use tuvx_core,                     only : core_t
    use tuvx_netcdf,                   only : clean_string
    use tuvx_profile,                  only : profile_t
    use tuvx_solver,                   only : radiation_field_t

    class(output_t), intent(inout) :: this
    !> Output step (time, sza, etc.)
    integer,            intent(in) :: step
    !> TUV-x core
    class(core_t),      intent(in) :: core
    !> Photolysis rate constants [s-1] (vertical level, reaction)
    real(dk), optional, intent(in) :: photolysis_rate_constants(:,:)
    !> Dose rates (vertical level, dose rate type)
    real(dk), optional, intent(in) :: dose_rates(:,:)
    !> Time [hours]
    real(dk), optional, intent(in) :: time
    !> Solar zenith angle [degrees]
    real(dk), optional, intent(in) :: solar_zenith_angle
    !> Earth-Sun distance [AU]
    real(dk), optional, intent(in) :: earth_sun_distance

    character(len=*), parameter :: Iam = "TUV-x results output"
    integer        :: i_rate, i_elem
    real(kind=dk), allocatable :: values_2D(:,:)
    type(radiation_field_t) :: rad_field
    class(profile_t), pointer :: profile
    type(string_t) :: var_name, append_dim, dim_names(2), units

    append_dim = "time"

    if( present( time ) ) then
      var_name = "time"
      units = "hours"
      call this%file_%append( var_name, units, append_dim, step, time, Iam )
    end if

    if( present( solar_zenith_angle ) ) then
      var_name = "solar zenith angle"
      units = "degrees"
      call this%file_%append( var_name, units, append_dim, step,              &
                              solar_zenith_angle, Iam )
    end if

    if( present( earth_sun_distance ) ) then
      var_name = "Earth-Sun distance"
      units = "AU"
      call this%file_%append( var_name, units, append_dim, step,              &
                              earth_sun_distance, Iam )
    end if

    profile => core%get_profile( "temperature", "K" )
    var_name = "temperature"
    units = "K"
    dim_names(1) = "vertical_level"
    call this%file_%append( var_name, units, append_dim, step, dim_names(1),  &
                            profile%edge_val_, Iam )
    deallocate( profile )

    if( this%do_radiation_ ) then
      rad_field = core%get_radiation_field( )
      units = "photon s-1 nm-1 cm-2"
      dim_names(1) = "vertical_level"
      dim_names(2) = "wavelength"
      var_name = "direct radiation"
      call this%file_%append( var_name, units, append_dim, step,              &
                              dim_names(1:2), rad_field%fdr_, Iam )
      var_name = "upward radiation"
      call this%file_%append( var_name, units, append_dim, step,              &
                              dim_names(1:2), rad_field%fup_, Iam )
      var_name = "downward radiation"
      call this%file_%append( var_name, units, append_dim, step,              &
                              dim_names(1:2), rad_field%fdn_, Iam )
    end if

    if( present( photolysis_rate_constants ) ) then
      call assert_msg( 316406965, this%do_photo_, "Photolysis rate constants "&
                       //"are not configured to be output" )
      dim_names(1) = "vertical_level"
      units = "s-1"
      do i_rate = 1, size( this%photo_labels_ )
        var_name = clean_string( this%photo_labels_( i_rate ) )
        call this%file_%append( var_name, units, append_dim, step,            &
                                dim_names(1),                                 &
                                photolysis_rate_constants( :, i_rate ), Iam )
      end do
    end if

    if( present( dose_rates ) ) then
      call assert_msg( 981693955, this%do_dose_, "Dose rates are not "        &
                       //"configured to be output" )
      dim_names(1) = "vertical_level"
      units = "various"
      do i_rate = 1, size( this%dose_labels_ )
        var_name = clean_string( this%dose_labels_( i_rate ) )
        call this%file_%append( var_name, units, append_dim, step,            &
                                dim_names(1), dose_rates( :, i_rate ), Iam )
      end do
    end if

    dim_names(1) = "vertical_level"
    dim_names(2) = "wavelength"
    units = "cm2 molecule-1"
    do i_elem = 1, size( this%photo_cross_sections_ )
    associate( label => this%photo_cross_sections_( i_elem ) )
      values_2D = core%get_photolysis_cross_section( label )
      var_name = "cross section "//label
      call this%file_%append( var_name, units, append_dim, step,              &
                              dim_names(1:2), values_2D, Iam )
    end associate
    end do

    dim_names(1) = "vertical_level"
    dim_names(2) = "wavelength"
    units = "unitless"
    do i_elem = 1, size( this%photo_quantum_yields_ )
    associate( label => this%photo_quantum_yields_( i_elem ) )
      values_2D = core%get_photolysis_quantum_yield( label )
      var_name = "quantum yield "//label
      call this%file_%append( var_name, units, append_dim, step,              &
                              dim_names(1:2), values_2D, Iam )
    end associate
    end do

  end subroutine output

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Adds grids to the output
  !!
  !! Grids are used as dimensions in the output
  subroutine add_grids( this, core )

    use tuvx_core,                     only : core_t
    use tuvx_grid,                     only : grid_t

    class(output_t), intent(inout) :: this
    class(core_t),   intent(inout) :: core

    character(len=*), parameter :: Iam = "Output grid dimension adder"
    type(string_t)              :: var_name, dim_name, units
    class(grid_t), pointer      :: grid

    var_name = "altitude"
    dim_name = "vertical_level"
    units    = "km"
    grid => core%get_grid( "height", "km" )
    call this%file_%write( var_name, dim_name, grid%edge_, Iam )
    call this%file_%set_variable_units( var_name, units, Iam )
    deallocate( grid )

    var_name = "wavelength"
    dim_name = "wavelength"
    units    = "nm"
    grid => core%get_grid( "wavelength", "nm" )
    call this%file_%write( var_name, dim_name, grid%mid_, Iam )
    call this%file_%set_variable_units( var_name, units, Iam )
    deallocate( grid )

  end subroutine add_grids

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Adds photolysis diagnostics to the output
  subroutine add_photolysis_diagnostics( this, config )

    use musica_config,                 only : config_t
    use musica_iterator,               only : iterator_t

    class(output_t), intent(inout) :: this
    type(config_t),  intent(inout) :: config

    character(len=*), parameter :: my_name = "Photolysis output diagnostics"
    logical :: found, do_output
    class(iterator_t), pointer :: iter
    type(config_t) :: photolysis, rxns, photo, cs, qy
    type(string_t) :: output_cs, output_qy, photo_label

    output_cs = ""
    output_qy = ""
    call config%get( "photolysis", photolysis, my_name, found = found )
    if( .not. found ) then
      allocate( this%photo_cross_sections_( 0 ) )
      allocate( this%photo_quantum_yields_( 0 ) )
      return
    end if
    call photolysis%get( "reactions", rxns, my_name )
    iter => rxns%get_iterator( )
    do while( iter%next( ) )
      call rxns%get( iter, photo, my_name )
      call photo%get( "name", photo_label, my_name )
      call photo%get( "cross section", cs, my_name )
      call cs%get( "__output", do_output, my_name, default = .false. )
      if( do_output ) output_cs = output_cs//photo_label//";"
      call photo%get( "quantum yield", qy, my_name )
      call qy%get( "__output", do_output, my_name, default = .false. )
      if( do_output ) output_qy = output_qy//photo_label//";"
    end do
    this%photo_cross_sections_ = output_cs%split( ";", compress = .true. )
    this%photo_quantum_yields_ = output_qy%split( ";", compress = .true. )
    deallocate( iter )

  end subroutine add_photolysis_diagnostics

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Cleans up memory
  subroutine finalize( this )

    type(output_t), intent(inout) :: this

    if( associated( this%file_ ) ) deallocate( this%file_ )

  end subroutine finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_output
