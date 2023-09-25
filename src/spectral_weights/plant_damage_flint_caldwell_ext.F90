! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0

module tuvx_spectral_weight_plant_damage_flint_caldwell_ext
  ! The Flint-Caldwell plant damage extension spectral weight type and
  ! related functions

  use tuvx_spectral_weight,            only : spectral_weight_t
  use musica_constants,                only : dk => musica_dk

  implicit none

  private
  public :: spectral_weight_plant_damage_flint_caldwell_ext_t

  type, extends(spectral_weight_t) ::                                         &
      spectral_weight_plant_damage_flint_caldwell_ext_t
    ! Calculator for Flint-Caldwell plant damage extension spectral weight
  contains
    procedure :: calculate => run
  end type spectral_weight_plant_damage_flint_caldwell_ext_t

  interface spectral_weight_plant_damage_flint_caldwell_ext_t
    module procedure constructor
  end interface spectral_weight_plant_damage_flint_caldwell_ext_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function constructor( config, grid_warehouse, profile_warehouse )           &
      result ( this )
    ! Initialize the Flint-Caldwell plant damage extension spectral weight

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
                     "plant damage flint caldwell ext spectral wght." )

    allocate( spectral_weight_plant_damage_flint_caldwell_ext_t :: this )
    this%wavelength_grid_ = grid_warehouse%get_ptr( "wavelength", "nm" )

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function run( this, grid_warehouse, profile_warehouse )                     &
      result( spectral_weight )
    ! Calculate the Flint-Caldwell plant damage extension spectral weight

    use tuvx_grid,                     only : grid_t
    use tuvx_grid_warehouse,           only : grid_warehouse_t
    use tuvx_profile_warehouse,        only : profile_warehouse_t

    class(spectral_weight_plant_damage_flint_caldwell_ext_t), intent(in) &
      :: this ! This :f:type:`~tuvx_spectral_weight_plant_damage_flint_caldwell_ext/spectral_weight_plant_damage_flint_caldwell_ext_t`
    type(grid_warehouse_t),    intent(inout) :: grid_warehouse ! A :f:type:`~tuvx_grid_warehouse/grid_warehouse_t`
    type(profile_warehouse_t), intent(inout) :: profile_warehouse ! A :f:type:`~tuvx_profile_warehouse/profile_warehouse_t`
    real(kind=dk), allocatable               :: spectral_weight(:) ! The calculated spectral weights (wavelength) [unitless]

    ! Local variables
    real(dk), parameter  :: a0 = 4.688272_dk
    real(dk), parameter  :: a1 = .1703411_dk
    real(dk), parameter  :: w1 = 307.867_dk
    real(dk), parameter  :: w2 = 390._dk
    class(grid_t), pointer      :: lambdaGrid

    lambdaGrid => grid_warehouse%get_grid( this%wavelength_grid_ )

    spectral_weight = exp( a0 * exp( -exp( a1                                 &
                                     * ( lambdaGrid%mid_ - w1 ) / 1.15_dk ) ) &
                   + ( ( w2 - lambdaGrid%mid_ ) / 121.7557_dk - 4.183832_dk ) )
    spectral_weight = spectral_weight * lambdaGrid%mid_ / 300._dk
    where( spectral_weight < 0.0_dk .or. lambdaGrid%mid_ > 390._dk )
      spectral_weight = 0.0_dk
    endwhere

    if( associated( lambdaGrid ) ) deallocate( lambdaGrid )

  end function run

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_spectral_weight_plant_damage_flint_caldwell_ext
