! Copyright (C) 2024 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0

module tuvx_quantum_yield_taylor_series
  ! A quantum yield calculator based on a Taylor series

  use musica_config,                   only : config_t
  use tuvx_quantum_yield,              only : quantum_yield_t, base_constructor

  implicit none

  private
  public :: quantum_yield_taylor_series_t

  type, extends(quantum_yield_t) :: quantum_yield_taylor_series_t
    ! Calculator of quantum yields using a Taylor series
  contains
  end type quantum_yield_taylor_series_t

  interface quantum_yield_taylor_series_t
    module procedure :: constructor
  end interface quantum_yield_taylor_series_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function constructor( config, grid_warehouse, profile_warehouse )           &
      result( this )
    ! Constructor

    use musica_assert,                 only : assert_msg
    use musica_constants,              only : dk => musica_dk
    use musica_string,                 only : string_t
    use tuvx_grid,                     only : grid_t
    use tuvx_grid_warehouse,           only : grid_warehouse_t
    use tuvx_profile_warehouse,        only : profile_warehouse_t

    class(quantum_yield_t),    pointer       :: this ! This :f:type:`~tuvx_quantum_yield/quantum_yield_t` calculator
    type(config_t),            intent(inout) :: config ! Quantum yield configuration data
    type(grid_warehouse_t),    intent(inout) :: grid_warehouse ! A :f:type:`~tuvx_grid_warehouse/grid_warehouse_t`
    type(profile_warehouse_t), intent(inout) :: profile_warehouse ! A :f:type:`~tuvx_profile_warehouse/profile_warehouse_t`

    character(len=*), parameter :: my_name = "Taylor-series quantum yield constructor"
    type(string_t) :: required_keys(1), optional_keys(9)
    real(kind=dk) :: min_wl, max_wl
    class(grid_t), pointer :: wavelengths
    real(kind=dk), allocatable :: coeff(:)
    integer :: i_wl, i_coeff

    required_keys(1) = "coefficients"
    optional_keys(1) = "type"
    optional_keys(2) = "netcdf files"
    optional_keys(3) = "lower extrapolation"
    optional_keys(4) = "upper extrapolation"
    optional_keys(5) = "name"
    optional_keys(6) = "constant value"
    optional_keys(7) = "override bands"
    optional_keys(8) = "minimum wavelength"
    optional_keys(9) = "maximum wavelength"
    
    allocate( quantum_yield_taylor_series_t :: this )

    call assert_msg( 268043622,                                               &
                     config%validate( required_keys, optional_keys ),         &
                     "Bad configuration for Taylor-series quantum yield." )
    call base_constructor( this, config, grid_warehouse, profile_warehouse )

    call config%get( "coefficients", coeff, my_name )
    call assert_msg( 589032591, size( coeff ) .ge. 1,                         &
                     "Taylor-series quantum yield must have at least one "//  &
                     "coefficient.")
    call config%get( "minimum wavelength", min_wl, my_name, default = 0.0_dk )
    call config%get( "maximum wavelength", max_wl, my_name,                   &
                     default = huge(1.0_dk) )
    wavelengths => grid_warehouse%get_grid( this%wavelength_grid_ )
    call assert_msg( 401342404, size( this%quantum_yield_parms ) .eq. 1,      &
                     "Taylor-series quantum yield cannot be used with "//     &
                     "multiple data files" )
    associate( params => this%quantum_yield_parms(1)%array )
    do i_wl = 1, wavelengths%ncells_
      if( wavelengths%mid_( i_wl ) .lt. min_wl .or.                           &
          wavelengths%mid_( i_wl ) .gt. max_wl ) cycle
      params( i_wl, 1 ) = coeff(1)
      do i_coeff = 2, size( coeff )
        params( i_wl, 1 ) = params( i_wl, 1 ) +                               &
                                coeff( i_coeff )                              &
                                * wavelengths%mid_( i_wl )**( i_coeff - 1 )
      end do
      params( i_wl, 1 ) = max( 0.0, min( 1.0, params( i_wl, 1 ) ) )
    end do
    end associate
    deallocate( wavelengths )

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_quantum_yield_taylor_series