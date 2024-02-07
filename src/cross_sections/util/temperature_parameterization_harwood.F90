! Copyright (C) 2024 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0

module tuvx_temperature_parameterization_harwood
! Calculates cross-section elements using a temperature-based
! parameterization. TODO: need reference

  use musica_config,                   only : config_t
  use musica_constants,                only : dk => musica_dk
  use tuvx_temperature_parameterization,                                      &
      only : temperature_parameterization_t
  use tuvx_temperature_range,          only : temperature_range_t

  implicit none

  private
  public :: temperature_parameterization_harwood_t

  !> Parameterization for calculating cross section values
  !! TODO: need reference
  !!
  !! Cross section elements are calculated as:
  !!
  !! \f[
  !! \sigma(T,\lambda_i) = 10^{\(aa_i + bb_i / T\)}
  !! \f]
  !!
  !! where aa_i and bb_i are constants, T is temperature [K], and
  !! \f$\lambda_i\f$ is wavelength [nm]. The size of the aa and bb
  !! arrays must equal the number of wavelengths in the parameterization
  !! range.
  type, extends(temperature_parameterization_t) :: temperature_parameterization_harwood_t
  contains
    !> Calculate the cross section value for a specific temperature and wavelength
    procedure :: calculate
  end type temperature_parameterization_harwood_t

  !> Constructor for temperature_parameterization_harwood_t
  interface temperature_parameterization_harwood_t
    module procedure :: constructor
  end interface temperature_parameterization_harwood_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Constructs a Harwood (TODO: need ref) temperature-based parameterization
  function constructor( config, wavelengths ) result( this )

    use musica_assert,                 only : assert_msg, die_msg
    use musica_config,                 only : config_t
    use musica_iterator,               only : iterator_t
    use musica_string,                 only : string_t, to_char
    use tuvx_grid,                     only : grid_t

    type(temperature_parameterization_harwood_t) :: this
    type(config_t), intent(inout)                :: config
    class(grid_t),  intent(in)                   :: wavelengths

    character(len=*), parameter :: my_name =                                  &
        "Harwood temperature parameterization"
    type(string_t) :: required_keys(5), optional_keys(5), exp_base
    type(config_t) :: temp_ranges, temp_range
    class(iterator_t), pointer :: iter
    integer :: i_range, i_param, n_param, n_wl
    logical :: found

    required_keys(1) = "aa"
    required_keys(2) = "bb"
    required_keys(3) = "base temperature"
    required_keys(4) = "base wavelength"
    required_keys(5) = "logarithm"
    optional_keys(1) = "type"
    optional_keys(2) = "minimum wavelength"
    optional_keys(3) = "maximum wavelength"
    optional_keys(4) = "temperature ranges"
    optional_keys(5) = "invert temperature offset"
    call assert_msg( 581965121,                                               &
                     config%validate( required_keys, optional_keys ),         &
                     "Bad configuration for Harwood temperature "//           &
                     "parameterization." )
    call config%get( "aa", this%aa_, my_name )
    call config%get( "bb", this%bb_, my_name )
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
      call die_msg( 768514789, "Invalid logarithm type in Harwood "//         &
                               "temperature-based cross section: '"//         &
                               exp_base//"'" )
    end if
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
    call assert_msg( 654743205,                                               &
             wavelengths%mid_( this%min_wavelength_index_ )                   &
               >= this%min_wavelength_,                                       &
             "Minimum wavelength for Harawood temperature-based cross "//     &
             "section is outside the bounds of the wavelength grid." )
    this%max_wavelength_index_ = wavelengths%ncells_
    do while( wavelengths%mid_( this%max_wavelength_index_ )                  &
                > this%max_wavelength_                                        &
              .and. this%max_wavelength_index_ >= 1 )
      this%max_wavelength_index_ = this%max_wavelength_index_ - 1
    end do
    call assert_msg( 309165090,                                               &
             wavelengths%mid_( this%max_wavelength_index_ )                   &
               <= this%max_wavelength_,                                       &
             "Maximum wavelength for Harwood temperature-based cross "//      &
             "section is outside the bounds of the wavelength grid." )
    ! TODO This follows logic from original TUV, but perhaps should
    !      be modified to assign TUV-x wavelength edges
    this%wavelengths_ = wavelengths%mid_( this%min_wavelength_index_ :        &
                                          this%max_wavelength_index_ )
    call assert_msg( 760344004, size( this%aa_ ) .eq. size( this%bb_ ),       &
                     "Parameter arrays for Harwood temperature-based cross "//&
                     "section (aa and bb) must be of the same size." )
    n_param = size( this%aa_ )
    n_wl = this%max_wavelength_index_ - this%min_wavelength_index_ + 1
    call assert_msg( 641308113, n_param .eq. n_wl,                            &
                     "Parameter arrays for Harwood temperature-based cross "//&
                     "section (aa and bb) must be the same size as the "//    &
                     "parameterized portion of the wavelength grid. "//       &
                     "Expected "//trim( to_char( n_wl ) )//" but got "//      &
                     trim( to_char( n_param ) )//"." )
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
    allocate( this%lp_( 0 ) )

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine calculate( this, temperature, wavelengths, cross_section )

    use tuvx_profile,                  only : profile_t

    class(temperature_parameterization_harwood_t), intent(in) :: this
    real(kind=dk), intent(in)    :: temperature
    real(kind=dk), intent(in)    :: wavelengths(:)
    real(kind=dk), intent(inout) :: cross_section(:)

    real(kind=dk) :: temp
    integer :: i_range, i_wl, w_min, w_max

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
      do i_wl = 1, w_max - w_min + 1
        cross_section( w_min + i_wl - 1 ) =                                   &
            10.0d0**( this%aa_( i_wl ) + this%bb_( i_wl ) / temp )
      end do
    end associate
    end do

  end subroutine calculate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_temperature_parameterization_harwood