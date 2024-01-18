! Copyright (C) 2020-4 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0

module tuvx_temperature_parameterization_taylor_series
! Calculates cross-section elements using a Taylor-series temperature-based
! parameterization

  ! Including musica_config at the module level to avoid an ICE
  ! with Intel 2022.1 compiler
  use musica_config,                   only : config_t
  use musica_constants,                only : dk => musica_dk
  use tuvx_temperature_parameterization,                                      &
      only : temperature_parameterization_t
  use tuvx_temperature_range,          only : temperature_range_t

  implicit none

  private
  public :: temperature_parameterization_taylor_series_t

  !> Parameters for calculating cross section values based on
  !! temperature using a Taylor series
  !!
  !! Cross section elements are calculated as:
  !!
  !! \f[
  !! \sigma(T_{base}) * \[ 1.0 + A_1 * (T - T_{base}) + A_2 * (T - T_{base})^2 \]
  !! \f]
  !!
  !! where \f$\sigma\f$ is a reference cross section at temperature [K]
  !! \f$T_{base}\f$, \f$A_1\f$ and \f$A_2\f$ are fitting parameters, and
  !! \f$T\f$ is temperature [K].
  type, extends(temperature_parameterization_t) :: temperature_parameterization_taylor_series_t
    !> Wavelength grid for temperature parameterization [nm]
    real(kind=dk), allocatable :: wavelengths_(:)
    !> Base cross section element
    real(kind=dk), allocatable :: sigma_(:)
    !> Taylor-series coefficients A_n (n,wavelength)
    real(kind=dk), allocatable :: A_(:,:)
  contains
    !> Calculate the cross section value for a specific temperature and wavelength
    procedure :: calculate
    !> Returns the number of bytes required to pack the parameterization
    !! onto a character buffer
    procedure :: pack_size => pack_size
    !> Packs the parameterization onto a character buffer
    procedure :: mpi_pack => mpi_pack
    !> Unpacks the parameterization from a character buffer
    procedure :: mpi_unpack => mpi_unpack
  end type temperature_parameterization_taylor_series_t

  !> Constructor for temperature_parameterization_taylor_series_t
  interface temperature_parameterization_taylor_series_t
    module procedure :: constructor
  end interface temperature_parameterization_taylor_series_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Constructs a Taylor-series temperature-based parameterization
  function constructor( config, wavelengths ) result ( this )

    use musica_assert,                 only : assert_msg
    use musica_config,                 only : config_t
    use musica_iterator,               only : iterator_t
    use musica_string,                 only : string_t
    use tuvx_grid,                     only : grid_t
    use tuvx_netcdf,                   only : netcdf_t

    type(temperature_parameterization_taylor_series_t) :: this
    type(config_t), intent(inout)                      :: config
    class(grid_t),  intent(in)                         :: wavelengths

    character(len=*), parameter :: my_name =                                  &
        "Taylor-series temperature parameterization constructor"
    type(string_t) :: required_keys(2), optional_keys(4), file_path
    type(config_t) :: temp_ranges, temp_range, netcdf_file
    class(iterator_t), pointer :: iter
    type(netcdf_t) :: netcdf
    integer :: i_range, i_param, n_param
    logical :: found

    required_keys(1) = "netcdf file"
    required_keys(2) = "base temperature"
    optional_keys(1) = "minimum wavelength"
    optional_keys(2) = "maximum wavelength"
    optional_keys(3) = "temperature ranges"
    optional_keys(4) = "type"
    call assert_msg( 235183546,                                               &
                     config%validate( required_keys, optional_keys ),         &
                     "Bad configuration for temperature parameterization." )

    ! Load NetCDF file
    call config%get( "netcdf file", netcdf_file, my_name )
    call netcdf_file%get( "file path", file_path, my_name )
    call netcdf%read_netcdf_file( file_path = file_path%to_char( ),           &
                                  variable_name = "temperature_" )
    n_param = size( netcdf%parameters, dim = 2 ) - 1
    call assert_msg( 164185428, n_param >= 1, "Taylor-series temperature "//  &
                     "parameterization must have at least one set of "//      &
                     "coefficients" )
    allocate( this%A_( n_param, size( netcdf%wavelength ) ) )
    this%wavelengths_ = netcdf%wavelength
    this%sigma_ = netcdf%parameters(:,1)
    do i_param = 1, n_param
      this%A_( i_param, : ) = netcdf%parameters( : , i_param + 1 )
    end do

    call config%get( "base temperature", this%base_temperature_, my_name )
    call config%get( "minimum wavelength", this%min_wavelength_, my_name,     &
                     default = 0.0_dk )
    call config%get( "maximum wavelength", this%max_wavelength_, my_name,     &
                     default = huge( 1.0_dk ) )
    this%min_wavelength_index_ = 1
    do while( wavelengths%mid_( this%min_wavelength_index_ )                  &
                < this%max_wavelength_                                        &
              .and. this%min_wavelength_index_ <= wavelengths%ncells_ )
      this%min_wavelength_index_ = this%min_wavelength_index_ + 1
    end do
    call assert_msg( 504874740,                                               &
             wavelengths%mid_( this%min_wavelength_index_ )                   &
               >= this%min_wavelength_,                                       &
             "Minimum wavelength for Taylor-series temperature-based cross "//&
             "section is outside the bounds of the wavelength grid." )
    this%max_wavelength_index_ = wavelengths%ncells_
    do while( wavelengths%mid_( this%max_wavelength_index_ )                  &
                > this%max_wavelength_                                        &
              .and. this%max_wavelength_index_ >= 1 )
      this%max_wavelength_index_ = this%max_wavelength_index_ - 1
    end do
    call assert_msg( 587703546,                                               &
             wavelengths%mid_( this%max_wavelength_index_ )                   &
               <= this%max_wavelength_,                                       &
             "Maximum wavelength for Taylor-series temperature-based cross "//&
             "section is outside the bounds of the wavelength grid." )
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

    ! initialize unused data members
    allocate( this%AA_(0) )
    allocate( this%BB_(0) )
    allocate( this%lp_(0) )

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine calculate( this, temperature, wavelengths, cross_section )

    use tuvx_profile,                  only : profile_t

    class(temperature_parameterization_taylor_series_t), intent(in) :: this
    real(kind=dk), intent(in)    :: temperature
    real(kind=dk), intent(in)    :: wavelengths(:)
    real(kind=dk), intent(inout) :: cross_section(:)

  end subroutine calculate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the size of a character buffer required to pack the
  !! parameterization
  integer function pack_size( this, comm )

    use musica_mpi,                    only : musica_mpi_pack_size

    !> Parameterization to be packed
    class(temperature_parameterization_taylor_series_t), intent(in) :: this
    !> MPI communicator
    integer,                                             intent(in) :: comm

#ifdef MUSICA_USE_MPI
    pack_size = this%temperature_parameterization_t%pack_size( comm ) +    &
                musica_mpi_pack_size( this%wavelengths_, comm ) +             &
                musica_mpi_pack_size( this%sigma_,       comm ) +             &
                musica_mpi_pack_size( this%A_,           comm )
#else
    pack_size = 0
#endif

  end function pack_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Packs the parameterization onto a character buffer
  subroutine mpi_pack( this, buffer, position, comm )

    use musica_assert,                 only : assert
    use musica_mpi,                    only : musica_mpi_pack

    !> Parameterization to be packed
    class(temperature_parameterization_taylor_series_t), intent(in) :: this
    !> Memory buffer
    character, intent(inout) :: buffer(:)
    !> Current buffer position
    integer,   intent(inout) :: position
    !> MPI communicator
    integer,   intent(in)    :: comm

#ifdef MUSICA_USE_MPI
    integer :: prev_pos

    prev_pos = position
    call this%temperature_parameterization_t%mpi_pack( buffer, position, comm )
    call musica_mpi_pack( buffer, position, this%wavelengths_, comm )
    call musica_mpi_pack( buffer, position, this%sigma_,       comm )
    call musica_mpi_pack( buffer, position, this%A_,           comm )
    call assert( 342538714, position - prev_pos <= this%pack_size( comm ) )
#endif

  end subroutine mpi_pack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Unpacks a parameterization from a character buffer
  subroutine mpi_unpack( this, buffer, position, comm )

    use musica_assert,                 only : assert
    use musica_mpi,                    only : musica_mpi_unpack

    !> The parameterization to be unpacked
    class(temperature_parameterization_taylor_series_t), intent(out) :: this
    !> Memory buffer
    character, intent(inout) :: buffer(:)
    !> Current buffer position
    integer,   intent(inout) :: position
    !> MPI communicator
    integer,   intent(in)    :: comm

#ifdef MUSICA_USE_MPI
    integer :: prev_pos

    prev_pos = position
    call this%temperature_parameterization_t%mpi_unpack( buffer, position, comm )
    call musica_mpi_unpack( buffer, position, this%wavelengths_, comm )
    call musica_mpi_unpack( buffer, position, this%sigma_,       comm )
    call musica_mpi_unpack( buffer, position, this%A_,           comm )
#endif

  end subroutine mpi_unpack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_temperature_parameterization_taylor_series