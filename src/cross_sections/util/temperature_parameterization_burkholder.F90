! Copyright (C) 2024 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0

module tuvx_temperature_parameterization_burkholder
! Calculates cross-section elements using a temperature-based
! parameterization from Burkholder et al. Phys. Chem. Chem. Phys. 4, 1432-1437 (2002).

  ! Including musica_config at the module level to avoid an ICE
  ! with Intel 2022.1 compiler
  use musica_config,                   only : config_t
  use musica_constants,                only : dk => musica_dk
  use tuvx_temperature_parameterization,                                      &
      only : temperature_parameterization_t
  use tuvx_temperature_range,          only : temperature_range_t

  implicit none

  private
  public :: temperature_parameterization_burkholder_t

  !> Parameters for calculating cross section values based on
  !! temperature using the algoritm in Burkholder et al.
  !! Phys. Chem. Chem. Phys. 4, 1432-1437 (2002).
  !!
  !! Cross section elements are calculated as:
  !!
  !! \f[
  !! Q(T) = 1 + e^{\frac{A}{B*T}}
  !! \sigma(T,\lambda) = \frac{aa(\lambda)}{Q(T)} + bb(\lambda)*\[1-\frac{1}{Q(T)}\]
  !! \f]
  !!
  !! where A, B, aa, and bb are constants, T is temperature [K] and \f$\lambda\f$ is
  !! wavelength [nm].
  type, extends(temperature_parameterization_t) :: temperature_parameterization_burkholder_t
    real(kind=dk) :: A_
    real(kind=dk) :: B_
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
  end type temperature_parameterization_burkholder_t

  !> Constructor for temperature_parameterization_burkholder_t
  interface temperature_parameterization_burkholder_t
    module procedure :: constructor
  end interface temperature_parameterization_burkholder_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Constructs a Burkholder (2002) temperature-based parameterization
  function constructor( config ) result ( this )

    use musica_assert,                 only : assert_msg
    use musica_config,                 only : config_t
    use musica_iterator,               only : iterator_t
    use musica_string,                 only : string_t
    use tuvx_grid,                     only : grid_t
    use tuvx_netcdf,                   only : netcdf_t

    type(temperature_parameterization_burkholder_t) :: this
    type(config_t), intent(inout)                      :: config

    character(len=*), parameter :: my_name =                                  &
        "Burkholder (2002) temperature parameterization constructor"
    type(string_t) :: required_keys(3), optional_keys(4), file_path
    type(config_t) :: temp_ranges, temp_range, netcdf_file
    class(iterator_t), pointer :: iter
    type(netcdf_t) :: netcdf
    integer :: i_range, i_param, n_param
    logical :: found

    required_keys(1) = "netcdf file"
    required_keys(2) = "A"
    required_keys(3) = "B"
    optional_keys(1) = "type"
    optional_keys(2) = "temperature ranges"
    optional_keys(3) = "minimum wavelength"
    optional_keys(4) = "maximum wavelength"
    call assert_msg( 235183546,                                               &
                     config%validate( required_keys, optional_keys ),         &
                     "Bad configuration for Burkholder (2002) temperature "// &
                     "parameterization." )

    ! Load NetCDF file
    call config%get( "netcdf file", netcdf_file, my_name )
    call netcdf_file%get( "file path", file_path, my_name )
    call netcdf%read_netcdf_file( file_path = file_path%to_char( ),           &
                                  variable_name = "temperature_" )
    n_param = size( netcdf%parameters, dim = 2 )
    call assert_msg( 164185428, n_param >= 2, "Burkholder (2002) "//          &
                     "parameterization must have at two sets of "//           &
                     "coefficients" )

    call config%get( "minimum wavelength", this%min_wavelength_, my_name,     &
                     default = netcdf%wavelength(1) )
    call config%get( "maximum wavelength", this%max_wavelength_, my_name,     &
                     default = netcdf%wavelength( size( netcdf%wavelength ) ) )
    this%min_wavelength_ = max( this%min_wavelength_, netcdf%wavelength(1) )
    this%max_wavelength_ = min( this%max_wavelength_,                         &
                               netcdf%wavelength( size( netcdf%wavelength ) ) )
    call assert_msg( 856954069, this%min_wavelength_ < this%max_wavelength_,  &
                     "Invalid wavelength range for Burkholder temperature "// &
                     "parameterization" )

    ! Load parameters
    call config%get( "A", this%A_, my_name )
    call config%get( "B", this%B_, my_name )
    this%wavelengths_ = netcdf%wavelength(:)
    this%AA_ = netcdf%parameters(:,1)
    this%BB_ = netcdf%parameters(:,2)
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
    allocate( this%lp_(0) )

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine calculate( this, temperature, wavelengths, cross_section )

    use tuvx_profile,                  only : profile_t

    class(temperature_parameterization_burkholder_t), intent(in) :: this
    real(kind=dk), intent(in)    :: temperature
    real(kind=dk), intent(in)    :: wavelengths(:)
    real(kind=dk), intent(inout) :: cross_section(:)

    real(kind=dk) :: temp, Q
    integer :: i_range, w_min, w_max

    w_min = this%min_wavelength_index_
    w_max = this%max_wavelength_index_
    do i_range = 1, size( this%ranges_ )
    associate( temp_range => this%ranges_( i_range ) )
      if( temperature < temp_range%min_temperature_ .or.                      &
          temperature > temp_range%max_temperature_ ) cycle
      if( temp_range%is_fixed_ ) then
        temp = temp_range%fixed_temperature_ - this%base_temperature_
      else
        temp = temperature - this%base_temperature_
      end if
      Q = 1.0_dk + exp( this%A_ / ( this%B_ * temp ) )
      cross_section( w_min:w_max ) = ( this%AA_(:) / Q +                      &
                                       this%BB_(:) * ( 1.0_dk - 1.0_dk / Q )  &
                                     ) * 1.0e-20_dk
    end associate
    end do

  end subroutine calculate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the size of a character buffer required to pack the
  !! parameterization
  integer function pack_size( this, comm )

    use musica_mpi,                    only : musica_mpi_pack_size

    !> Parameterization to be packed
    class(temperature_parameterization_burkholder_t), intent(in) :: this
    !> MPI communicator
    integer,                                             intent(in) :: comm

#ifdef MUSICA_USE_MPI
    pack_size = this%temperature_parameterization_t%pack_size( comm ) +       &
                musica_mpi_pack_size( this%A_, comm ) +                       &
                musica_mpi_pack_size( this%B_, comm )
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
    class(temperature_parameterization_burkholder_t), intent(in) :: this
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
    call musica_mpi_pack( buffer, position, this%A_, comm )
    call musica_mpi_pack( buffer, position, this%B_, comm )
    call assert( 190816083, position - prev_pos <= this%pack_size( comm ) )
#endif

  end subroutine mpi_pack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Unpacks a parameterization from a character buffer
  subroutine mpi_unpack( this, buffer, position, comm )

    use musica_assert,                 only : assert
    use musica_mpi,                    only : musica_mpi_unpack

    !> The parameterization to be unpacked
    class(temperature_parameterization_burkholder_t), intent(out) :: this
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
    call musica_mpi_unpack( buffer, position, this%A_, comm )
    call musica_mpi_unpack( buffer, position, this%B_, comm )
    call assert( 634825156, position - prev_pos <= this%pack_size( comm ) )
#endif

  end subroutine mpi_unpack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_temperature_parameterization_burkholder