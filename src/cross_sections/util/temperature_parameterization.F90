! Copyright (C) 2020-2025 University Corporation for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0

module tuvx_temperature_parameterization
! Calculates cross-section elements based on a temperature parameterization

  ! Including musica_config at the module level to avoid an ICE
  ! with Intel 2022.1 compiler
  use musica_config,                   only : config_t
  use musica_constants,                only : dk => musica_dk
  use tuvx_temperature_range,          only : temperature_range_t

  implicit none

  private
  public :: temperature_parameterization_t

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
    !> Wavelengths in parameterization range [nm]
    real(kind=dk), allocatable :: wavelengths_(:)
    !> Base temperature [K] to use in calculations
    real(kind=dk) :: base_temperature_ = 0.0_dk
    !> Base wavelength [nm] to use in calcuations
    real(kind=dk) :: base_wavelength_ = 0.0_dk
    !> Flag indicating whether cross section algorithm is base 10 (true)
    !! or base e (false)
    logical :: is_base_10_ = .true.
    !> Flad indicating whether to subtract base temperature from
    !! actual temperature (false) or to subtract actual temperature
    !! from base temperature (true)
    logical :: is_temperature_inverted_ = .false.
    !> Minimum wavelength [nm] to calculate values for
    real(kind=dk) :: min_wavelength_ = 0.0_dk
    !> Maximum wavelength [nm] to calculate values for
    real(kind=dk) :: max_wavelength_ = 0.0_dk
    !> Index of minimum wavelength [nm] to calculate values for
    integer :: min_wavelength_index_ = 0
    !> Index of maximum wavelength to calculate values for
    integer :: max_wavelength_index_ = 0
    !> Temperature ranges used in parameterization
    type(temperature_range_t), allocatable :: ranges_(:)
  contains
    !> Merges NetCDF wavelength grid with parameterization grid
    procedure :: merge_wavelength_grids
    !> Calculate the cross section value for a specific temperature
    !! and wavelength
    procedure :: calculate => calculate
    !> Returns the number of bytes required to pack the parameterization
    !! onto a character buffer
    procedure :: pack_size => pack_size
    !> Packs the parameterization onto a character buffer
    procedure :: mpi_pack => mpi_pack
    !> Unpacks the parameterization from a character buffer
    procedure :: mpi_unpack => mpi_unpack
  end type temperature_parameterization_t

  !> Constructor for temperature_parameterization_t
  interface temperature_parameterization_t
    module procedure :: constructor
  end interface temperature_parameterization_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function constructor( config, wavelengths ) result( this )
    ! Constructs temperature_parameterization_t objects

    use musica_assert,                 only : assert_msg, die_msg
    use musica_iterator,               only : iterator_t
    use musica_string,                 only : string_t
    use tuvx_grid,                     only : grid_t

    type(temperature_parameterization_t)               :: this
    type(config_t),                      intent(inout) :: config
    class(grid_t),                       intent(in)    :: wavelengths

    character(len=*), parameter :: my_name =                                  &
        "temperature parameterization constructor"
    type(string_t) :: required_keys(6), optional_keys(4), exp_base
    type(config_t) :: temp_ranges, temp_range
    class(iterator_t), pointer :: iter
    integer :: i_range
    logical :: found

    required_keys(1) = "AA"
    required_keys(2) = "BB"
    required_keys(3) = "lp"
    required_keys(4) = "base temperature"
    required_keys(5) = "base wavelength"
    required_keys(6) = "logarithm"
    optional_keys(1) = "minimum wavelength"
    optional_keys(2) = "maximum wavelength"
    optional_keys(3) = "temperature ranges"
    optional_keys(4) = "invert temperature offset"
    call assert_msg( 256315527,                                               &
                     config%validate( required_keys, optional_keys ),         &
                     "Bad configuration for temperature parameterization." )

    call config%get( "AA", this%AA_, my_name )
    call config%get( "BB", this%BB_, my_name )
    call config%get( "lp", this%lp_, my_name )
    call config%get( "base temperature", this%base_temperature_, my_name )
    call config%get( "base wavelength",  this%base_wavelength_,  my_name )
    call config%get( "logarithm", exp_base, my_name )
    call config%get( "invert temperature offset",                             &
                     this%is_temperature_inverted_, my_name, default = .false.)
    if( exp_base == "base 10" ) then
      this%is_base_10_ = .true.
    else if( exp_base == "natural" ) then
      this%is_base_10_ = .false.
    else
      call die_msg( 104603249, "Invalid logarithm type in temperature-based"//&
                               " cross section: '"//exp_base//"'" )
    end if
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
    ! TODO This follows logic from original TUV, but perhaps should
    !      be modified to assign TUV-x wavelength edges
    this%wavelengths_ = wavelengths%mid_( this%min_wavelength_index_ :        &
                                          this%max_wavelength_index_ )
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

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function merge_wavelength_grids( this, input_grid ) result( merged_grid )
    ! Merges wavelength grid from NetCDF input data with parameterization
    ! grid.
    ! Where they overlap, the parameterization is used.
    ! Updates the parameterization wavelength indices for new grid.
    ! Returns merged wavelength grid.

    use musica_assert,                 only : assert
    use tuvx_grid,                     only : grid_t

    class(temperature_parameterization_t), intent(inout) :: this
    real(kind=dk),                         intent(in)    :: input_grid(:)
    real(kind=dk), allocatable                           :: merged_grid(:)

    logical :: found_min
    integer :: i_wl, n_wl, i_input_wl, i_param_wl

    if( size( input_grid ) == 0 ) then
      merged_grid = this%wavelengths_
      this%min_wavelength_index_ = 1
      this%max_wavelength_index_ = size( merged_grid )
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
    n_wl = n_wl + size( this%wavelengths_ )
    allocate( merged_grid( n_wl ) )
    i_input_wl = 1
    i_param_wl = 1
    i_wl = 1
    found_min = .false.
    do
      if( i_wl > n_wl ) then
        ! end of merged grid
        exit
      else if( i_param_wl > size( this%wavelengths_ ) .and.                                      &
               input_grid( i_input_wl ) <= max_wl ) then
        ! skipping input data wavelengths in parameterization range
        i_input_wl = i_input_wl + 1
      else if( .not. ( min_wl <= input_grid( i_input_wl ) .and.               &
                       max_wl >= input_grid( i_input_wl ) ) ) then
        ! adding input data wavelengths outside of parameterization range
        merged_grid( i_wl ) = input_grid( i_input_wl )
        i_input_wl = i_input_wl + 1
        i_wl = i_wl + 1
      else if( i_param_wl <= size( this%wavelengths_ ) ) then
        ! adding TUV-x wavelengths in parameterization range
        merged_grid( i_wl ) = this%wavelengths_( i_param_wl )
        if( .not. found_min ) then
          found_min = .true.
          wl_min_index = i_wl
        end if
        wl_max_index = i_wl
        i_param_wl = i_param_wl + 1
        i_wl = i_wl + 1
      end if
    end do
    call assert( 265861594, i_param_wl == size( this%wavelengths_ ) + 1 )
    call assert( 537808229, i_input_wl <= size( input_grid ) + 1 )
    call assert( 422870529, i_wl       == n_wl + 1 )
    end associate

  end function merge_wavelength_grids

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine calculate( this, temperature, wavelengths, cross_section )

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
      if ( this%is_temperature_inverted_ ) then
        temp = this%base_temperature_ - temp
      else
        temp = temp - this%base_temperature_
      end if
      temp_xs(:) = 0.0_dk
      do i_lp = 1, size( this%lp_ )
        temp_xs( w_min:w_max ) = temp_xs( w_min:w_max ) +                     &
            ( this%AA_( i_lp ) + temp * this%BB_( i_lp ) ) *                  &
              ( wavelengths( w_min:w_max )                                    &
                - this%base_wavelength_ )**this%lp_( i_lp )
      end do
      if (this%is_base_10_) then
        cross_section( w_min:w_max ) = cross_section( w_min:w_max )           &
                                       + 10**temp_xs( w_min:w_max )
      else
        cross_section( w_min:w_max ) = cross_section( w_min:w_max )           &
                                       + exp( temp_xs( w_min:w_max ) )
      end if
    end associate
    end do

  end subroutine calculate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function pack_size( this, comm )
    ! Returns the size of a character buffer required to pack the
    ! parameterization

    use musica_mpi,                    only : musica_mpi_pack_size

    class(temperature_parameterization_t), intent(in) :: this ! parameterization to be packed
    integer,                               intent(in) :: comm ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: i_range

    pack_size = musica_mpi_pack_size( this%AA_,                      comm ) + &
                musica_mpi_pack_size( this%BB_,                      comm ) + &
                musica_mpi_pack_size( this%lp_,                      comm ) + &
                musica_mpi_pack_size( this%wavelengths_,             comm ) + &
                musica_mpi_pack_size( this%base_temperature_,        comm ) + &
                musica_mpi_pack_size( this%base_wavelength_,         comm ) + &
                musica_mpi_pack_size( this%is_base_10_,              comm ) + &
                musica_mpi_pack_size( this%is_temperature_inverted_, comm ) + &
                musica_mpi_pack_size( this%min_wavelength_,          comm ) + &
                musica_mpi_pack_size( this%max_wavelength_,          comm ) + &
                musica_mpi_pack_size( this%min_wavelength_index_,    comm ) + &
                musica_mpi_pack_size( this%max_wavelength_index_,    comm ) + &
                musica_mpi_pack_size( allocated( this%ranges_ ),     comm )
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

  end function pack_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine mpi_pack( this, buffer, position, comm )
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
    call musica_mpi_pack( buffer, position, this%wavelengths_,          comm )
    call musica_mpi_pack( buffer, position, this%base_temperature_,     comm )
    call musica_mpi_pack( buffer, position, this%base_wavelength_,      comm )
    call musica_mpi_pack( buffer, position, this%is_base_10_,           comm )
    call musica_mpi_pack( buffer, position, this%is_temperature_inverted_,    &
                          comm )
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

  end subroutine mpi_pack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine mpi_unpack( this, buffer, position, comm )
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
    call musica_mpi_unpack( buffer, position, this%wavelengths_,         comm )
    call musica_mpi_unpack( buffer, position, this%base_temperature_,    comm )
    call musica_mpi_unpack( buffer, position, this%base_wavelength_,     comm )
    call musica_mpi_unpack( buffer, position, this%is_base_10_,          comm )
    call musica_mpi_unpack( buffer, position, this%is_temperature_inverted_,  &
                            comm )
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

  end subroutine mpi_unpack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_temperature_parameterization