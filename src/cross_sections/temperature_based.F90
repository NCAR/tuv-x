! Copyright (C) 2020-4 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0

module tuvx_cross_section_temperature_based
! Calculate a cross section using a temperature-based parameterization

  ! Including musica_config at the module level to avoid an ICE
  ! with Intel 2022.1 compiler
  use musica_config,                   only : config_t
  use musica_constants,                only : dk => musica_dk
  use tuvx_cross_section,              only : cross_section_t
  use tuvx_interpolate,                only : interpolator_conserving_t
  use tuvx_temperature_parameterization,                                      &
      only : temperature_parameterization_t

  implicit none

  private
  public :: cross_section_temperature_based_t

  integer, parameter :: PARAM_BASE          = 1
  integer, parameter :: PARAM_TAYLOR_SERIES = 2
  integer, parameter :: PARAM_BURKHOLDER    = 3
  integer, parameter :: PARAM_HARWOOD       = 4

  !> Calculator for temperature-based cross sections
  type, extends(cross_section_t) :: cross_section_temperature_based_t
    real(kind=dk), allocatable :: raw_wavelengths_(:) ! [nm]
    real(kind=dk), allocatable :: raw_data_(:)
    class(temperature_parameterization_t), pointer :: parameterization_ =>    &
                                                                       null( )
    type(interpolator_conserving_t) :: interpolator_
  contains
    !> Calculate the cross section
    procedure :: calculate
    !> Returns the number of bytes required to pack the cross section onto
    !! a character buffer
    procedure :: pack_size
    !> Packs the cross section onto a character buffer
    procedure :: mpi_pack
    !> Unpacks a cross section from a character buffer
    procedure :: mpi_unpack
    !> Clean up memory
    final :: finalize
  end type cross_section_temperature_based_t

  !> Constructor
  interface cross_section_temperature_based_t
    module procedure :: constructor
  end interface cross_section_temperature_based_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function constructor( config, grid_warehouse, profile_warehouse )           &
      result( this )
    ! Constructs cross_section_temperature_based_t objects

    use musica_assert,                 only : assert, assert_msg, die_msg
    use musica_string,                 only : string_t
    use tuvx_cross_section,            only : base_constructor
    use tuvx_grid,                     only : grid_t
    use tuvx_grid_factory,             only : grid_builder
    use tuvx_grid_warehouse,           only : grid_warehouse_t
    use tuvx_netcdf,                   only : netcdf_t
    use tuvx_profile_warehouse,        only : profile_warehouse_t
    use tuvx_temperature_parameterization_burkholder,                         &
        only : temperature_parameterization_burkholder_t
    use tuvx_temperature_parameterization_harwood,                         &
        only : temperature_parameterization_harwood_t
    use tuvx_temperature_parameterization_taylor_series,                      &
        only : temperature_parameterization_taylor_series_t

    class(cross_section_t),    pointer       :: this
    type(config_t),            intent(inout) :: config
    type(grid_warehouse_t),    intent(inout) :: grid_warehouse
    type(profile_warehouse_t), intent(inout) :: profile_warehouse

    ! local variables
    character(len=*), parameter :: my_name =                                  &
        'Temperature-based cross section constructor'
    type(string_t) :: required_keys(2), optional_keys(3)
    class(grid_t), pointer :: wavelengths
    type(config_t) :: param_config, interpolator_config, grid_config
    type(string_t) :: file_path, param_type
    type(netcdf_t) :: netcdf
    real(kind=dk), allocatable :: file_data(:), file_wl(:)
    logical :: found
    integer :: i_file, i_wl

    required_keys(1) = "type"
    required_keys(2) = "parameterization"
    optional_keys(1) = "name"
    optional_keys(2) = "parameterization wavelength grid"
    optional_keys(3) = "netcdf file"
    call assert_msg( 483410000,                                               &
                     config%validate( required_keys, optional_keys ),         &
                     "Bad configuration for temperature-based cross section" )

    allocate( cross_section_temperature_based_t :: this )

    ! Get grid and profile pointers
    this%wavelength_grid_ = grid_warehouse%get_ptr( "wavelength", "nm" )
    this%height_grid_ = grid_warehouse%get_ptr( "height", "km" )
    this%temperature_profile_ = profile_warehouse%get_ptr( "temperature", "K" )

    ! Load NetCDF files
    call config%get( "netcdf file", file_path, my_name, found = found )
    if( found ) then
      call netcdf%read_netcdf_file( file_path = file_path%to_char( ),         &
                                    variable_name = "cross_section_" )
      call assert_msg( 793476078, size( netcdf%parameters, dim = 2 ) == 1,    &
                       "File: "//file_path//" should contain 1 parameter" )
      file_data = netcdf%parameters(:,1)
      file_wl   = netcdf%wavelength(:)
    else
      allocate( file_data(0) )
      allocate( file_wl(0) )
    end if

    ! Check for custom wavelength grid for parameterization
    call config%get( "parameterization wavelength grid", grid_config, my_name,&
                     found = found)
    if( found ) then
      wavelengths => grid_builder( grid_config )
      call assert_msg( 993335233, wavelengths%units( ) .eq. "nm",             &
                       "Invalid units for custom wavelength grid in "//       &
                       "temperature-based cross section. Expected 'nm' "//    &
                       "but got '"//wavelengths%units( )//"'" )
    else
      wavelengths => grid_warehouse%get_grid( this%wavelength_grid_ )
    end if

    ! Load parameters
    select type( this )
    type is( cross_section_temperature_based_t )
      call config%get( "parameterization", param_config, my_name )
      call param_config%get( "type", param_type, my_name, found = found )
      if( found ) then
        if( param_type == "TAYLOR_SERIES" ) then
          allocate( this%parameterization_, source =                          &
            temperature_parameterization_taylor_series_t( param_config ) )
        else if( param_type == "BURKHOLDER" ) then
          allocate( this%parameterization_, source =                          &
            temperature_parameterization_burkholder_t( param_config ) )
        else if( param_type == "HARWOOD" ) then
          allocate( this%parameterization_, source =                          &
            temperature_parameterization_harwood_t( param_config,             &
                                                    wavelengths ) )
        else
          call die_msg( 370773773, "Invalid temperature-based "//             &
                        "parameterization type: '"//param_type//"'" )
        end if
      else
        allocate( this%parameterization_, source =                            &
            temperature_parameterization_t( param_config, wavelengths ) )
      end if
      this%raw_wavelengths_ =                                                 &
          this%parameterization_%merge_wavelength_grids( file_wl )
      allocate( this%raw_data_( size( this%raw_wavelengths_ ) ) )
      i_file = 1
      do i_wl = 1, size( this%raw_wavelengths_ )
        if( i_wl >= this%parameterization_%min_wavelength_index_ .and.        &
            i_wl <= this%parameterization_%max_wavelength_index_ ) then
          this%raw_data_( i_wl ) = 0.0_dk
          cycle
        end if
        do while( file_wl( i_file ) < this%raw_wavelengths_( i_wl ) )
          i_file = i_file + 1
        end do
        this%raw_data_( i_wl ) = file_data( i_file )
        i_file = i_file + 1
      end do
      call assert( 950874524, i_file <= size( file_data ) + 1 )
    end select
    deallocate( wavelengths )

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function calculate( this, grid_warehouse, profile_warehouse, at_mid_point ) &
      result( cross_section )
    ! Calculates cross section by combining NetCDF data and temperature-based
    ! parameterization results

    use tuvx_grid,                     only : grid_t
    use tuvx_grid_warehouse,           only : grid_warehouse_t
    use tuvx_profile_warehouse,        only : profile_warehouse_t
    use tuvx_profile,                  only : profile_t
    use tuvx_util,                     only : add_point

    real(kind=dk), allocatable                           :: cross_section(:,:)
    class(cross_section_temperature_based_t), intent(in) :: this
    type(grid_warehouse_t),                intent(inout) :: grid_warehouse
    type(profile_warehouse_t),             intent(inout) :: profile_warehouse
    logical, optional,                     intent(in)    :: at_mid_point

    ! local variables
    character(len=*),  parameter :: Iam =                                      &
      'Temperature-based cross section calculate'
    real(kind=dk),   parameter :: deltax = 1.0e-5
    class(grid_t),     pointer :: heights
    class(grid_t),     pointer :: wavelengths
    class(profile_t),  pointer :: temperatures
    real(kind=dk)              :: temperature
    real(kind=dk), allocatable :: raw_data(:), raw_wl(:)
    logical                    :: l_at_mid_point
    integer                    :: i_wl, i_height

    ! Add temperature-based cross section values
    temperatures => profile_warehouse%get_profile( this%temperature_profile_ )
    wavelengths  => grid_warehouse%get_grid( this%wavelength_grid_ )
    if( present( at_mid_point ) ) then
      l_at_mid_point = at_mid_point
    else
      l_at_mid_point = .false.
    end if
    if( l_at_mid_point ) then
      allocate( cross_section( temperatures%size( ), wavelengths%size( ) ) )
    else
      allocate( cross_section( temperatures%size( ) + 1,                      &
                               wavelengths%size( ) ) )
    end if
    do i_height = 1, size( cross_section, 1 )
      if( l_at_mid_point ) then
        temperature = temperatures%mid_val_( i_height )
      else
        temperature = temperatures%edge_val_( i_height )
      end if
      raw_data = this%raw_data_
      raw_wl   = this%raw_wavelengths_
      call this%parameterization_%calculate( temperature, raw_wl, raw_data )
      call add_point( x = raw_wl, y = raw_data,                               &
                      xnew = ( 1.0_dk - deltax ) * raw_wl(1), ynew = 0.0_dk )
      call add_point( x = raw_wl, y = raw_data,                               &
                      xnew = 0.0_dk, ynew = 0.0_dk )
      call add_point( x = raw_wl, y = raw_data,                               &
                      xnew = ( 1.0_dk + deltax ) * raw_wl( size( raw_wl ) ),  &
                      ynew = 0.0_dk )
      call add_point( x = raw_wl, y = raw_data,                               &
                      xnew = 1.0e38_dk, ynew = 0.0_dk )
      cross_section( i_height, : ) =                                          &
          this%interpolator_%interpolate( x_target = wavelengths%edge_,       &
                                          x_source = raw_wl,                  &
                                          y_source = raw_data,                &
                                          requested_by =                      &
                           "temperature based cross section wavelength grid" )
    end do
    deallocate( temperatures )
    deallocate( wavelengths  )

  end function calculate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function pack_size( this, comm )
    ! Returns the size of a character buffer required to pack the cross
    ! section

    use musica_mpi,                    only : musica_mpi_pack_size

    class(cross_section_temperature_based_t), intent(in) :: this ! cross section to be packed
    integer,                                  intent(in) :: comm ! MPI communicator

#ifdef MUSICA_USE_MPI
    pack_size = this%cross_section_t%pack_size( comm ) +                      &
                musica_mpi_pack_size( this%raw_wavelengths_, comm ) +         &
                musica_mpi_pack_size( this%raw_data_,        comm ) +         &
                musica_mpi_pack_size( .false.,               comm )
    if( associated( this%parameterization_ ) ) then
      pack_size = pack_size +                                                 &
                  musica_mpi_pack_size( 1, comm ) +                           &
                  this%parameterization_%pack_size( comm )
    end if
#else
    pack_size = this%cross_section_t%pack_size( comm )
#endif

  end function pack_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine mpi_pack( this, buffer, position, comm )
    ! Packs the cross section onto a character buffer

    use musica_assert,                 only : assert, die
    use musica_mpi,                    only : musica_mpi_pack
    use tuvx_temperature_parameterization_burkholder,                         &
        only : temperature_parameterization_burkholder_t
    use tuvx_temperature_parameterization_harwood,                            &
        only : temperature_parameterization_harwood_t
    use tuvx_temperature_parameterization_taylor_series,                      &
        only : temperature_parameterization_taylor_series_t

    class(cross_section_temperature_based_t), intent(in)    :: this      ! cross section to be packed
    character,                                intent(inout) :: buffer(:) ! memory buffer
    integer,                                  intent(inout) :: position  ! current buffer position
    integer,                                  intent(in)    :: comm      ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: prev_pos, param_type
    logical :: is_alloced

    is_alloced = associated( this%parameterization_ )

    prev_pos = position
    call this%cross_section_t%mpi_pack( buffer, position, comm )
    call musica_mpi_pack( buffer, position, this%raw_wavelengths_, comm )
    call musica_mpi_pack( buffer, position, this%raw_data_,        comm )
    call musica_mpi_pack( buffer, position, is_alloced,            comm )
    if( is_alloced ) then
      select type( param => this%parameterization_ )
        type is( temperature_parameterization_t )
          param_type = PARAM_BASE
        type is( temperature_parameterization_taylor_series_t )
          param_type = PARAM_TAYLOR_SERIES
        type is( temperature_parameterization_burkholder_t )
          param_type = PARAM_BURKHOLDER
        type is( temperature_parameterization_harwood_t )
          param_type = PARAM_HARWOOD
        class default
          call die( 424852458 )
      end select
      call musica_mpi_pack( buffer, position, param_type, comm )
      call this%parameterization_%mpi_pack( buffer, position, comm )
    end if
    call assert( 322345685, position - prev_pos <= this%pack_size( comm ) )
#endif

  end subroutine mpi_pack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine mpi_unpack( this, buffer, position, comm )
    ! Unpacks a cross section from a character buffer

    use musica_assert,                 only : assert, die
    use musica_mpi,                    only : musica_mpi_unpack
    use tuvx_temperature_parameterization_burkholder,                         &
        only : temperature_parameterization_burkholder_t
    use tuvx_temperature_parameterization_harwood,                            &
        only : temperature_parameterization_harwood_t
    use tuvx_temperature_parameterization_taylor_series,                      &
        only : temperature_parameterization_taylor_series_t

    class(cross_section_temperature_based_t), intent(out)   :: this      ! cross section to be unpacked
    character,                                intent(inout) :: buffer(:) ! memory buffer
    integer,                                  intent(inout) :: position  ! current buffer position
    integer,                                  intent(in)    :: comm      ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: prev_pos, param_type
    logical :: is_alloced

    prev_pos = position
    call this%cross_section_t%mpi_unpack( buffer, position, comm )
    call musica_mpi_unpack( buffer, position, this%raw_wavelengths_, comm )
    call musica_mpi_unpack( buffer, position, this%raw_data_,        comm )
    call musica_mpi_unpack( buffer, position, is_alloced,            comm )
    if( is_alloced ) then
      call musica_mpi_unpack( buffer, position, param_type, comm )
      select case( param_type )
        case( PARAM_BASE )
          allocate( temperature_parameterization_t :: this%parameterization_ )
        case( PARAM_TAYLOR_SERIES )
          allocate( temperature_parameterization_taylor_series_t ::          &
                    this%parameterization_ )
        case( PARAM_BURKHOLDER )
          allocate( temperature_parameterization_burkholder_t ::             &
                    this%parameterization_ )
        case( PARAM_HARWOOD )
          allocate( temperature_parameterization_harwood_t ::                &
                    this%parameterization_ )
        case default
          call die( 324803089 )
      end select
      call this%parameterization_%mpi_unpack( buffer, position, comm )
    end if
    call assert( 820834544, position - prev_pos <= this%pack_size( comm ) )
#endif

  end subroutine mpi_unpack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Clean-up memory
  subroutine finalize( this )

    type(cross_section_temperature_based_t), intent(inout) :: this

    if( associated( this%parameterization_ ) ) then
      deallocate( this%parameterization_ )
    end if

  end subroutine finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_cross_section_temperature_based
