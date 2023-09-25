! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0

module tuvx_spectral_weight_phytoplankton_boucher
  ! The phytoplankton boucher spectral weight type and related functions

  use tuvx_spectral_weight,            only : spectral_weight_t
  use musica_constants,                only : dk => musica_dk

  implicit none

  private
  public :: spectral_weight_phytoplankton_boucher_t

  type, extends(spectral_weight_t) :: spectral_weight_phytoplankton_boucher_t
    ! Calculator for phytoplankton boucher spectral weight
  contains
    procedure :: calculate => run
  end type spectral_weight_phytoplankton_boucher_t

  interface spectral_weight_phytoplankton_boucher_t
    module procedure constructor
  end interface spectral_weight_phytoplankton_boucher_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function constructor( config, grid_warehouse, profile_warehouse )           &
      result ( this )
    ! Initialize the phytoplankton boucher spectral weight

    use musica_assert,                 only : assert_msg
    use musica_config,                 only : config_t
    use musica_string,                 only : string_t
    use tuvx_grid_warehouse,           only : grid_warehouse_t
    use tuvx_profile_warehouse,        only : profile_warehouse_t

    class(spectral_weight_t),  pointer       :: this   ! This :f:type:`~tuvx_spectral_weight/spectral_weight_t`
    type(config_t),            intent(inout) :: config ! Spectral weight configuration data
    type(grid_warehouse_t),    intent(inout) :: grid_warehouse ! A :f:type:`~tuvx_grid_warehouse/grid_warehouse_t`
    type(profile_warehouse_t), intent(inout) :: profile_warehouse ! A :f:type:`~tuvx_profile_warehouse/profile_warehouse_t`

    ! Local variables
    type(string_t) :: required_keys(1), optional_keys(1)

    required_keys(1) = "type"
    optional_keys(1) = "name"
    call assert_msg( 399197234,                                               &
                     config%validate( required_keys, optional_keys ),         &
                     "Bad configuration data format for "//                   &
                     "phytoplankton boucher spectral wght." )

    allocate( spectral_weight_phytoplankton_boucher_t :: this )
    this%wavelength_grid_ = grid_warehouse%get_ptr( "wavelength", "nm" )

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function run( this, grid_warehouse, profile_warehouse )                     &
      result( spectral_weight )
    ! Calculate the phytoplankton boucher spectral weight

    use tuvx_grid,                     only : grid_t
    use tuvx_grid_warehouse,           only : grid_warehouse_t
    use tuvx_profile_warehouse,        only : profile_warehouse_t

    class(spectral_weight_phytoplankton_boucher_t),  intent(in)     :: this ! This :f:type:`~tuvx_spectral_weight_phytoplankton_boucher/spectral_weight_phytoplankton_boucher_t`
    type(grid_warehouse_t),    intent(inout) :: grid_warehouse ! A :f:type:`~tuvx_grid_warehouse/grid_warehouse_t`
    type(profile_warehouse_t), intent(inout) :: profile_warehouse ! A :f:type:`~tuvx_profile_warehouse/profile_warehouse_t`
    real(kind=dk), allocatable               :: spectral_weight(:) ! The calculated spectral weights (wavelength) [unitless]

    ! Local variables
    real(dk), parameter  :: em = -3.17e-6_dk
    real(dk), parameter  :: a  = 112.5_dk
    real(dk), parameter  :: b  = -.6223_dk
    real(dk), parameter  :: c  = 7.67e-4_dk
    class(grid_t), pointer      :: lambdaGrid

    lambdaGrid => grid_warehouse%get_grid( this%wavelength_grid_ )

    allocate( spectral_weight( lambdaGrid%ncells_ ) )

    where( lambdaGrid%mid_ > 290._dk .and. lambdaGrid%mid_ < 400._dk )
      spectral_weight = em + exp( a + lambdaGrid%mid_                         &
                                      * ( b + lambdaGrid%mid_ * c ) )
    elsewhere
      spectral_weight = 0.0_dk
    endwhere
    spectral_weight = max( 0.0_dk,spectral_weight )

    if( associated( lambdaGrid ) ) deallocate( lambdaGrid )

  end function run

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_spectral_weight_phytoplankton_boucher
