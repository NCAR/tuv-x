! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0

module tuvx_cross_section_temperature_based
! Calculate a cross section using a temperature-based parameterization

  ! Including musica_config at the module level to avoid an ICE
  ! with Intel 2022.1 compiler
  use musica_config,                   only : config_t
  use musica_constants,                only : dk => musica_dk
  use tuvx_cross_section,              only : cross_section_t
  use tuvx_interpolate,                only : interpolator_conserving_t

  implicit none

  private
  public :: cross_section_temperature_based_t

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
    procedure :: pack_size => temperature_range_pack_size
    !> Packs the range onto a character buffer
    procedure :: mpi_pack => temperature_range_mpi_pack
    !> Unpacks a range from a character buffer
    procedure :: mpi_unpack => temperature_range_mpi_unpack
  end type temperature_range_t

  !> Constructor for temperature_range_t
  interface temperature_range_t
    module procedure :: temperature_range_constructor
  end interface temperature_range_t

  !> Parameters for calculating cross section values based on
  !! temperature
  !!
  !! Cross section elements are calculated as:
  !!
  !! \f[
  !! 10^{\sum_i{(AA_i + (T-273)*BB_i)*\lambda^{lp_i}}}
  !! \f]
  !!
  !! where \f$\lambda\f$ is the wavelength [nm] and
  !! \f$T\f$ is the temperature [K].
  type :: temperature_parameterization_t
    integer :: n_sets_ = 0
    real(kind=dk), allocatable :: AA_(:)
    real(kind=dk), allocatable :: BB_(:)
    real(kind=dk), allocatable :: lp_(:)
    !> Minimum wavelength [nm] to calculate values for
    real(kind=dk) :: min_wavelength_
    !> Maximum wavelength [nm] to calculate values for
    real(kind=dk) :: max_wavelength_
    !> Index of minimum wavelength [nm] to calculate values for
    integer :: min_wavelength_index_
    !> Index of maximum wavelength to calculate values for
    integer :: max_wavelength_index_
    !> Temperature ranges used in parameterization
    type(temperature_range_t), allocatable :: ranges_(:)
  contains
    !> Merges NetCDF wavelength grid with parameterization grid
    procedure :: merge_wavelength_grids
    !> Calculate the cross section value for a specific temperature
    !! and wavelength
    procedure :: calculate => temperature_parameterization_calculate
    !> Returns the number of bytes required to pack the parameterization
    !! onto a character buffer
    procedure :: pack_size => temperature_parameterization_pack_size
    !> Packs the parameterization onto a character buffer
    procedure :: mpi_pack => temperature_parameterization_mpi_pack
    !> Unpacks the parameterization from a character buffer
    procedure :: mpi_unpack => temperature_parameterization_mpi_unpack
  end type temperature_parameterization_t

  !> Constructor for temperature_parameterization_t
  interface temperature_parameterization_t
    module procedure :: temperature_parameterization_constructor
  end interface temperature_parameterization_t

  !> Calculator for temperature-based cross sections
  type, extends(cross_section_t) :: cross_section_temperature_based_t
    real(kind=dk), allocatable :: raw_wavelengths_(:) ! [nm]
    real(kind=dk), allocatable :: raw_data_(:)
    type(temperature_parameterization_t) :: parameterization_
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

    use musica_assert,                 only : assert, assert_msg
    use musica_string,                 only : string_t
    use tuvx_cross_section,            only : base_constructor
    use tuvx_grid,                     only : grid_t
    use tuvx_grid_warehouse,           only : grid_warehouse_t
    use tuvx_netcdf,                   only : netcdf_t
    use tuvx_profile_warehouse,        only : profile_warehouse_t
    use tuvx_util,                     only : add_point

    class(cross_section_t),    pointer       :: this
    type(config_t),            intent(inout) :: config
    type(grid_warehouse_t),    intent(inout) :: grid_warehouse
    type(profile_warehouse_t), intent(inout) :: profile_warehouse

    ! local variables
    character(len=*), parameter :: my_name =                                  &
        'Temperature-based cross section constructor'
    real(kind=dk), parameter :: deltax = 1.0e-5
    type(string_t) :: required_keys(3), optional_keys(1)
    class(grid_t), pointer :: wavelengths
    type(config_t) :: param_config, interpolator_config
    type(string_t) :: file_path
    type(netcdf_t) :: netcdf
    real(kind=dk), allocatable :: file_data(:), file_wl(:)
    logical :: found
    integer :: i_file, i_wl

    required_keys(1) = "type"
    required_keys(2) = "parameterization"
    required_keys(3) = "netcdf file"
    optional_keys(1) = "name"
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
    call netcdf%read_netcdf_file( file_path = file_path%to_char( ),           &
                                  variable_name = "cross_section_" )
    call assert_msg( 793476078, size( netcdf%parameters, dim = 2 ) == 1,      &
                     "File: "//file_path//" should contain 1 parameter" )
    file_data = netcdf%parameters(:,1)
    file_wl   = netcdf%wavelength(:)
    call add_point( x = file_wl, y = file_data,                               &
                    xnew = ( 1.0_dk - deltax ) * file_wl(1), ynew = 0.0_dk )
    call add_point( x = file_wl, y = file_data,                               &
                    xnew = 0.0_dk, ynew = 0.0_dk )
    call add_point( x = file_wl, y = file_data,                               &
                    xnew = ( 1.0_dk + deltax ) * file_wl( size( file_wl ) ),  &
                    ynew = 0.0_dk )
    call add_point( x = file_wl, y = file_data,                               &
                    xnew = 1.0e38_dk, ynew = 0.0_dk )

    ! Load parameters
    select type( this )
    type is( cross_section_temperature_based_t )
      wavelengths => grid_warehouse%get_grid( this%wavelength_grid_ )
      call config%get( "parameterization", param_config, my_name )
      this%parameterization_ =                                                &
          temperature_parameterization_t( param_config, wavelengths )
      this%raw_wavelengths_ =                                                 &
          this%parameterization_%merge_wavelength_grids( file_wl, wavelengths )
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
      call assert( 950874524, i_file == size( file_data ) + 1 )
      deallocate( wavelengths )
    end select

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

    real(kind=dk), allocatable                           :: cross_section(:,:)
    class(cross_section_temperature_based_t), intent(in) :: this
    type(grid_warehouse_t),                intent(inout) :: grid_warehouse
    type(profile_warehouse_t),             intent(inout) :: profile_warehouse
    logical, optional,                     intent(in)    :: at_mid_point

    ! local variables
    character(len=*),  parameter :: Iam =                                      &
      'Temperature-based cross section calculate'
    class(grid_t),     pointer :: heights
    class(grid_t),     pointer :: wavelengths
    class(profile_t),  pointer :: temperatures
    real(kind=dk)              :: temperature
    real(kind=dk), allocatable :: raw_data(:)
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
      call this%parameterization_%calculate( temperature,                     &
                                             this%raw_wavelengths_, raw_data )
      cross_section( i_height, : ) =                                          &
          this%interpolator_%interpolate( x_target = wavelengths%edge_,       &
                                          x_source = this%raw_wavelengths_,   &
                                          y_source = raw_data )
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
                musica_mpi_pack_size( this%raw_data_, comm ) +                &
                this%parameterization_%pack_size( comm )
#else
    pack_size = this%cross_section_t%pack_size( comm )
#endif

  end function pack_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine mpi_pack( this, buffer, position, comm )
    ! Packs the cross section onto a character buffer

    use musica_assert,                 only : assert
    use musica_mpi,                    only : musica_mpi_pack

    class(cross_section_temperature_based_t), intent(in)    :: this      ! cross section to be packed
    character,                                intent(inout) :: buffer(:) ! memory buffer
    integer,                                  intent(inout) :: position  ! current buffer position
    integer,                                  intent(in)    :: comm      ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: prev_pos

    prev_pos = position
    call this%cross_section_t%mpi_pack( buffer, position, comm )
    call musica_mpi_pack( buffer, position, this%raw_wavelengths_, comm )
    call musica_mpi_pack( buffer, position, this%raw_data_, comm )
    call this%parameterization_%mpi_pack( buffer, position, comm )
    call assert( 322345685, position - prev_pos <= this%pack_size( comm ) )
#endif

  end subroutine mpi_pack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine mpi_unpack( this, buffer, position, comm )
    ! Unpacks a cross section from a character buffer

    use musica_assert,                 only : assert
    use musica_mpi,                    only : musica_mpi_unpack

    class(cross_section_temperature_based_t), intent(out)   :: this      ! cross section to be unpacked
    character,                                intent(inout) :: buffer(:) ! memory buffer
    integer,                                  intent(inout) :: position  ! current buffer position
    integer,                                  intent(in)    :: comm      ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: prev_pos

    prev_pos = position
    call this%cross_section_t%mpi_unpack( buffer, position, comm )
    call musica_mpi_unpack( buffer, position, this%raw_wavelengths_, comm )
    call musica_mpi_unpack( buffer, position, this%raw_data_, comm )
    call this%parameterization_%mpi_unpack( buffer, position, comm )
    call assert( 820834544, position - prev_pos <= this%pack_size( comm ) )
#endif

  end subroutine mpi_unpack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function temperature_range_constructor( config ) result( this )
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

  end function temperature_range_constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function temperature_range_pack_size( this, comm )                  &
      result( pack_size )
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

  end function temperature_range_pack_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine temperature_range_mpi_pack( this, buffer, position, comm )
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

  end subroutine temperature_range_mpi_pack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine temperature_range_mpi_unpack( this, buffer, position, comm )
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

  end subroutine temperature_range_mpi_unpack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function temperature_parameterization_constructor( config, wavelengths )    &
      result( this )
    ! Constructs temperature_parameterization_t objects

    use musica_assert,                 only : assert_msg
    use musica_config,                 only : config_t
    use musica_iterator,               only : iterator_t
    use musica_string,                 only : string_t
    use tuvx_grid,                     only : grid_t

    type(temperature_parameterization_t)               :: this
    type(config_t),                      intent(inout) :: config
    class(grid_t),                       intent(in)    :: wavelengths

    character(len=*), parameter :: my_name =                                  &
        "temperature parameterization constructor"
    type(string_t) :: required_keys(3), optional_keys(3)
    type(config_t) :: temp_ranges, temp_range
    class(iterator_t), pointer :: iter
    integer :: i_range
    logical :: found

    required_keys(1) = "AA"
    required_keys(2) = "BB"
    required_keys(3) = "lp"
    optional_keys(1) = "minimum wavelength"
    optional_keys(2) = "maximum wavelength"
    optional_keys(3) = "temperature ranges"
    call assert_msg( 256315527,                                               &
                     config%validate( required_keys, optional_keys ),         &
                     "Bad configuration for temperature parameterization." )

    call config%get( "AA", this%AA_, my_name )
    call config%get( "BB", this%BB_, my_name )
    call config%get( "lp", this%lp_, my_name )
    call assert_msg( 467090427, size( this%AA_ ) == size( this%BB_ ) .and.    &
                                size( this%AA_ ) == size( this%lp_ ),         &
                     "Arrays AA, BB, and lp must be the same size for "//     &
                     "temperature-based cross sections." )
    call config%get( "minimum wavelength", this%min_wavelength_, my_name,     &
                     default = 0.0_dk )
    call config%get( "maximum wavelength", this%max_wavelength_, my_name,     &
                     default = huge(1.0_dk) )
    this%min_wavelength_index_ = 1
    do while( wavelengths%mid_( this%min_wavelength_index_ )                  &
                < this%min_wavelength_                                        &
              .and. this%min_wavelength_index_ <= wavelengths%ncells_ )
      this%min_wavelength_index_ = this%min_wavelength_index_ + 1
    end do
    call assert_msg( 286143383,                                               &
             wavelengths%mid_( this%min_wavelength_index_ )                   &
               >= this%min_wavelength_,                                       &
             "Minimum wavelength for temperature-based cross section is "//   &
             "outside the bounds of the wavelength grid." )
    this%max_wavelength_index_ = wavelengths%ncells_
    do while( wavelengths%mid_( this%max_wavelength_index_ )                  &
                > this%max_wavelength_                                        &
              .and. this%max_wavelength_index_ >= 1 )
      this%max_wavelength_index_ = this%max_wavelength_index_ - 1
    end do
    call assert_msg( 490175140,                                               &
             wavelengths%mid_( this%max_wavelength_index_ )                   &
               <= this%max_wavelength_,                                       &
             "Maximum wavelength for temperature-based cross section is "//   &
             "outside the bounds of the wavelength grid." )
    call config%get( "temperature ranges", temp_ranges, my_name,              &
                     found = found )
    if( .not. found ) then
      allocate( this%ranges_( 1 ) )
      return
    end if
    allocate( this%ranges_( temp_ranges%number_of_children( ) ) )
    iter => temp_ranges%get_iterator( )
    i_range = 0
    do while( iter%next( ) )
      i_range = i_range + 1
      call temp_ranges%get( iter, temp_range, my_name )
      this%ranges_( i_range ) = temperature_range_t( temp_range )
    end do
    deallocate( iter )

  end function temperature_parameterization_constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function merge_wavelength_grids( this, input_grid, tuv_grid )               &
      result( merged_grid )
    ! Merges wavelength grid from NetCDF input data with parameterization
    ! grid (same as the TUV-x grid).
    ! Where they overlap, the parameterization is used.
    ! Updates the parameterization wavelength indices for new grid.
    ! Returns merged wavelength grid.
    !
    ! NOTE: Uses mid-points on the TUV-x wavelength grid

    use musica_assert,                 only : assert
    use tuvx_grid,                     only : grid_t

    class(temperature_parameterization_t), intent(inout) :: this
    real(kind=dk),                         intent(in)    :: input_grid(:)
    class(grid_t),                         intent(in)    :: tuv_grid
    real(kind=dk), allocatable                           :: merged_grid(:)

    logical :: found_min
    integer :: i_wl, n_wl, i_input_wl, i_tuv_wl, n_tuv_wl

    if( size( input_grid ) == 0 ) then
      merged_grid = tuv_grid%mid_
      return
    end if

    associate( wl_min_index => this%min_wavelength_index_,                    &
               wl_max_index => this%max_wavelength_index_,                    &
               min_wl       => this%min_wavelength_,                          &
               max_wl       => this%max_wavelength_ )
    n_wl = 0
    do i_input_wl = 1, size( input_grid(:) )
      if( min_wl > input_grid( i_input_wl ) .or.                              &
          max_wl < input_grid( i_input_wl ) ) n_wl = n_wl + 1
    end do
    i_tuv_wl   = wl_min_index
    n_tuv_wl   = wl_max_index
    n_wl       = n_wl + ( n_tuv_wl - i_tuv_wl + 1 )
    allocate( merged_grid( n_wl ) )
    i_input_wl = 1
    i_wl = 1
    found_min = .false.
    do
      if( i_wl > n_wl ) then
        ! end of merged grid
        exit
      else if( i_tuv_wl > n_tuv_wl .and.                                      &
               input_grid( i_input_wl ) <= max_wl ) then
        ! skipping input data wavelengths in parameterization range
        i_input_wl = i_input_wl + 1
      else if( .not. ( min_wl <= input_grid( i_input_wl ) .and.               &
                       max_wl >= input_grid( i_input_wl ) ) ) then
        ! adding input data wavelengths outside of parameterization range
        merged_grid( i_wl ) = input_grid( i_input_wl )
        i_input_wl = i_input_wl + 1
        i_wl = i_wl + 1
      else if( i_tuv_wl <= n_tuv_wl ) then
        ! adding TUV-x wavelengths in parameterization range
        !
        ! TODO This follows logic from original TUV, but perhaps should
        !      be modified to assign TUV-x wavelength edges
        merged_grid( i_wl ) = tuv_grid%mid_( i_tuv_wl )
        if( .not. found_min ) then
          found_min = .true.
          wl_min_index = i_wl
        end if
        wl_max_index = i_wl
        i_tuv_wl = i_tuv_wl + 1
        i_wl = i_wl + 1
      end if
    end do
    call assert( 265861594, i_tuv_wl   == n_tuv_wl + 1 )
    call assert( 537808229, i_input_wl == size( input_grid ) + 1 )
    call assert( 422870529, i_wl       == n_wl + 1 )
    end associate

  end function merge_wavelength_grids

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine temperature_parameterization_calculate( this, temperature,       &
      wavelengths, cross_section )

    use tuvx_profile,                  only : profile_t

    class(temperature_parameterization_t), intent(in)    :: this
    real(kind=dk),                         intent(in)    :: temperature
    real(kind=dk),                         intent(in)    :: wavelengths(:)
    real(kind=dk),                         intent(inout) :: cross_section(:)

    ! local variables
    real(kind=dk) :: temp, temp_xs( size( cross_section ) )
    integer :: i_lp, i_range, w_min, w_max

    w_min = this%min_wavelength_index_
    w_max = this%max_wavelength_index_
    do i_range = 1, size( this%ranges_ )
    associate( temp_range => this%ranges_( i_range ) )
      if( temperature < temp_range%min_temperature_ .or.       &
          temperature > temp_range%max_temperature_ ) cycle
      if( temp_range%is_fixed_ ) then
        temp = temp_range%fixed_temperature_
      else
        temp = temperature
      end if
      temp_xs(:) = 0.0_dk
      do i_lp = 1, size( this%lp_ )
        temp_xs( w_min:w_max ) = temp_xs( w_min:w_max ) +                     &
            ( this%AA_( i_lp ) + (temp - 273.0_dk) * this%BB_( i_lp ) ) *     &
              wavelengths( w_min:w_max )**this%lp_( i_lp )
      end do
      cross_section( w_min:w_max ) = cross_section( w_min:w_max )             &
                                     + 10**temp_xs( w_min:w_max )
    end associate
    end do

  end subroutine temperature_parameterization_calculate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function temperature_parameterization_pack_size( this, comm )       &
      result( pack_size )
    ! Returns the size of a character buffer required to pack the
    ! parameterization

    use musica_mpi,                    only : musica_mpi_pack_size

    class(temperature_parameterization_t), intent(in) :: this ! parameterization to be packed
    integer,                               intent(in) :: comm ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: i_range

    pack_size = musica_mpi_pack_size( this%AA_,                   comm ) +    &
                musica_mpi_pack_size( this%BB_,                   comm ) +    &
                musica_mpi_pack_size( this%lp_,                   comm ) +    &
                musica_mpi_pack_size( this%min_wavelength_,       comm ) +    &
                musica_mpi_pack_size( this%max_wavelength_,       comm ) +    &
                musica_mpi_pack_size( this%min_wavelength_index_, comm ) +    &
                musica_mpi_pack_size( this%max_wavelength_index_, comm ) +    &
                musica_mpi_pack_size( allocated( this%ranges_ ),  comm )
    if( allocated( this%ranges_ ) ) then
      pack_size = pack_size +                                                 &
                  musica_mpi_pack_size( size( this%ranges_ ), comm )
      do i_range = 1, size( this%ranges_ )
        pack_size = pack_size + this%ranges_( i_range )%pack_size( comm )
      end do
    end if
#else
    pack_size = 0
#endif

  end function temperature_parameterization_pack_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine temperature_parameterization_mpi_pack( this, buffer, position,   &
      comm )
    ! Packs the parameterization onto a character buffer

    use musica_assert,                 only : assert
    use musica_mpi,                    only : musica_mpi_pack

    class(temperature_parameterization_t), intent(in)    :: this      ! parameterization to be packed
    character,                             intent(inout) :: buffer(:) ! memory buffer
    integer,                               intent(inout) :: position  ! current buffer position
    integer,                               intent(in)    :: comm      ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: prev_pos, i_range

    prev_pos = position
    call musica_mpi_pack( buffer, position, this%AA_,                   comm )
    call musica_mpi_pack( buffer, position, this%BB_,                   comm )
    call musica_mpi_pack( buffer, position, this%lp_,                   comm )
    call musica_mpi_pack( buffer, position, this%min_wavelength_,       comm )
    call musica_mpi_pack( buffer, position, this%max_wavelength_,       comm )
    call musica_mpi_pack( buffer, position, this%min_wavelength_index_, comm )
    call musica_mpi_pack( buffer, position, this%max_wavelength_index_, comm )
    call musica_mpi_pack( buffer, position, allocated( this%ranges_ ),  comm )
    if( allocated( this%ranges_ ) ) then
      call musica_mpi_pack( buffer, position, size( this%ranges_ ), comm )
      do i_range = 1, size( this%ranges_ )
        call this%ranges_( i_range )%mpi_pack( buffer, position, comm )
      end do
    end if
    call assert( 267439201, position - prev_pos <= this%pack_size( comm ) )
#endif

  end subroutine temperature_parameterization_mpi_pack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine temperature_parameterization_mpi_unpack( this, buffer, position, &
      comm )
    ! Unpacks a parameterization from a character buffer

    use musica_assert,                 only : assert
    use musica_mpi,                    only : musica_mpi_unpack

    class(temperature_parameterization_t), intent(out)   :: this      ! parameterization to be unpacked
    character,                             intent(inout) :: buffer(:) ! memory buffer
    integer,                               intent(inout) :: position  ! current buffer position
    integer,                               intent(in)    :: comm      ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: prev_pos, i_range, n_ranges
    logical :: alloced

    prev_pos = position
    call musica_mpi_unpack( buffer, position, this%AA_,                  comm )
    call musica_mpi_unpack( buffer, position, this%BB_,                  comm )
    call musica_mpi_unpack( buffer, position, this%lp_,                  comm )
    call musica_mpi_unpack( buffer, position, this%min_wavelength_,      comm )
    call musica_mpi_unpack( buffer, position, this%max_wavelength_,      comm )
    call musica_mpi_unpack( buffer, position, this%min_wavelength_index_,comm )
    call musica_mpi_unpack( buffer, position, this%max_wavelength_index_,comm )
    call musica_mpi_unpack( buffer, position, alloced,                   comm )
    if( alloced ) then
      call musica_mpi_unpack( buffer, position, n_ranges, comm )
      allocate( this%ranges_( n_ranges ) )
      do i_range = 1, size( this%ranges_ )
        call this%ranges_( i_range )%mpi_unpack( buffer, position, comm )
      end do
    end if
    call assert( 483905106, position - prev_pos <= this%pack_size( comm ) )
#endif

  end subroutine temperature_parameterization_mpi_unpack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_cross_section_temperature_based
