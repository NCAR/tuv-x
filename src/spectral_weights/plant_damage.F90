! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0

module tuvx_spectral_weight_plant_damage
  ! The plant damage spectral weight type and related functions

  use tuvx_spectral_weight,            only : spectral_weight_t
  use musica_constants,                only : dk => musica_dk

  implicit none

  private
  public :: spectral_weight_plant_damage_t

  type, extends(spectral_weight_t) :: spectral_weight_plant_damage_t
    ! Calculator for plant damage spectral weight
  contains
    procedure :: calculate => run
  end type spectral_weight_plant_damage_t

  interface spectral_weight_plant_damage_t
    module procedure constructor
  end interface spectral_weight_plant_damage_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function constructor( config, grid_warehouse, profile_warehouse )           &
      result ( this )
    ! Initialize the plant damage spectral weight

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
                     "plant damage spectral wght." )

    allocate( spectral_weight_plant_damage_t :: this )
    this%wavelength_grid_ = grid_warehouse%get_ptr( "wavelength", "nm" )

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function run( this, grid_warehouse, profile_warehouse )                     &
      result( spectral_weight )
    ! Calculate the plant damage spectral weight

    use tuvx_grid,                     only : grid_t
    use tuvx_grid_warehouse,           only : grid_warehouse_t
    use tuvx_profile_warehouse,        only : profile_warehouse_t

    class(spectral_weight_plant_damage_t),  intent(in) :: this ! This :f:type:`~tuvx_spectral_weight_plant_damage/spectral_weight_plant_damage_t`
    type(grid_warehouse_t),    intent(inout) :: grid_warehouse ! A :f:type:`~tuvx_grid_warehouse/grid_warehouse_t`
    type(profile_warehouse_t), intent(inout) :: profile_warehouse ! A :f:type:`~tuvx_profile_warehouse/profile_warehouse_t`
    real(kind=dk), allocatable               :: spectral_weight(:) ! The calculated spectral weights (wavelength) [unitless]

    ! Local variables
    real(dk), parameter  :: a0 = 570.25_dk
    real(dk), parameter  :: a1 = -4.70144_dk
    real(dk), parameter  :: a2 = .01274_dk
    real(dk), parameter  :: a3 = -1.13118e-5_dk

    class(grid_t), pointer      :: lambdaGrid

    lambdaGrid => grid_warehouse%get_grid( this%wavelength_grid_ )

    spectral_weight = a0 + lambdaGrid%mid_                                    &
                     * ( a1 + lambdaGrid%mid_ * ( a2 + lambdaGrid%mid_ * a3 ) )
    where( spectral_weight < 0.0_dk .or. lambdaGrid%mid_ > 313._dk )
      spectral_weight = 0.0_dk
    endwhere

    if( associated( lambdaGrid ) ) deallocate( lambdaGrid )

  end function run

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_spectral_weight_plant_damage
