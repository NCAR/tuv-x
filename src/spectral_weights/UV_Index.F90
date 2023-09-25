! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0

module tuvx_spectral_weight_uv_index
  ! The uv index type and related functions

  use tuvx_spectral_weight,            only : spectral_weight_t
  use musica_constants,                only : dk => musica_dk

  implicit none

  private
  public :: spectral_weight_uv_index_t

  type, extends(spectral_weight_t) :: spectral_weight_uv_index_t
    ! Calculator for uv_index_spectral_weight
  contains
    procedure :: calculate => run
  end type spectral_weight_uv_index_t

  !> Constructor
  interface spectral_weight_uv_index_t
    module procedure constructor
  end interface spectral_weight_uv_index_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function constructor( config, grid_warehouse, profile_warehouse )           &
      result ( this )
    ! Initialize the spectral wght

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
                     "uv index spectral wght." )

    allocate( spectral_weight_uv_index_t :: this )
    this%wavelength_grid_ = grid_warehouse%get_ptr( "wavelength", "nm" )

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function run( this, grid_warehouse, profile_warehouse )                     &
      result( spectral_weight )
    ! Calculate the UV Index spectral weight

    use tuvx_grid,              only  :  grid_t
    use tuvx_grid_warehouse,    only  :  grid_warehouse_t
    use tuvx_profile_warehouse, only  :  profile_warehouse_t

    class(spectral_weight_uv_index_t),  intent(in)     :: this ! This :f:type:`~tuvx_spectral_weight_uv_index/spectral_weight_uv_index_t`
    type(grid_warehouse_t),    intent(inout) :: grid_warehouse ! A :f:type:`~tuvx_grid_warehouse/grid_warehouse_t`
    type(profile_warehouse_t), intent(inout) :: profile_warehouse ! A :f:type:`~tuvx_profile_warehouse/profile_warehouse_t`
    real(kind=dk), allocatable               :: spectral_weight(:) ! The calculated spectral weights (wavelength) [unitless]

    ! Local variables
    class(grid_t), pointer      :: lambdaGrid

    lambdaGrid => grid_warehouse%get_grid( this%wavelength_grid_ )

    spectral_weight = 40._dk * sw_fery( lambdaGrid%mid_ )

    if( associated( lambdaGrid ) ) deallocate( lambdaGrid )

  end function run

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function sw_fery( w ) result( fery )
    ! Calculate the action spectrum value for erythema at a given wavelength
    ! Webb, A.R., H. Slaper, P. Koepke, and A. W. Schmalwieser,
    ! Know your standard: Clarifying the CIE erythema action spectrum,
    ! Photochem. Photobiol. 87, 483-486, 2011.
    !
    ! Value at 300 nm = 0.6486
    ! 
    ! `doi:10.1111/j.1751-1097.2010.00871.x
    ! <https://doi.org/10.1111/j.1751-1097.2010.00871.x>`_

    real(dk), intent(in)  :: w(:)    ! Wavelength [nm]
    real(dk), allocatable :: fery(:) ! Action spectrum

    allocate( fery( size( w ) ) )

    where( w <= 298._dk )
      fery = 1._dk
    elsewhere( w > 298._dk .and. w <= 328._dk )
      fery = 10._dk**( .094_dk * ( 298._dk - w ) )
    elsewhere( w > 328._dk .and. w <= 400._dk )
      fery = 10._dk**( .015_dk * ( 140._dk - w ) )
    elsewhere
      fery = 1.e-36_dk
    endwhere

  end function sw_fery

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_spectral_weight_uv_index
