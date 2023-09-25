! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
module tuvx_spectral_weight_scup_mice
  ! The scup mice spectral weight type and related functions

  use tuvx_spectral_weight,            only : spectral_weight_t
  use musica_constants,                only : dk => musica_dk

  implicit none

  private
  public :: spectral_weight_scup_mice_t

  type, extends(spectral_weight_t) :: spectral_weight_scup_mice_t
    ! Calculator for Scup mice spectral weight
  contains
    procedure :: calculate => run
  end type spectral_weight_scup_mice_t

  interface spectral_weight_scup_mice_t
    module procedure constructor
  end interface spectral_weight_scup_mice_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function constructor( config, grid_warehouse, profile_warehouse )           &
      result ( this )
    ! Initialize the Scup mice spectral weight

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
                     "scup mice spectral wght." )

    allocate( spectral_weight_scup_mice_t :: this )
    this%wavelength_grid_ = grid_warehouse%get_ptr( "wavelength", "nm" )

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function run( this, grid_warehouse, profile_warehouse )                     &
      result( spectral_weight )
    ! Calculate the Scup mice spectral weight

    use tuvx_grid,                     only : grid_t
    use tuvx_grid_warehouse,           only : grid_warehouse_t
    use tuvx_profile_warehouse,        only : profile_warehouse_t

    class(spectral_weight_scup_mice_t),  intent(in)     :: this ! This :f:type:`~tuvx_spectral_weight_scup_mice/spectral_weight_scup_mice_t`
    type(grid_warehouse_t),    intent(inout) :: grid_warehouse ! A :f:type:`~tuvx_grid_warehouse/grid_warehouse_t`
    type(profile_warehouse_t), intent(inout) :: profile_warehouse ! A :f:type:`~tuvx_profile_warehouse/profile_warehouse_t`
    real(kind=dk), allocatable               :: spectral_weight(:) ! The calculated spectral weights (wavelength) [unitless]

    ! Local variables
    real(dk), allocatable       :: factor(:)

    class(grid_t), pointer      :: lambdaGrid

    lambdaGrid => grid_warehouse%get_grid( this%wavelength_grid_ )
    factor = 1._dk / sw_futr( (/ 300._dk /) )
    spectral_weight = sw_futr( lambdaGrid%mid_ ) * factor(1)

    if( associated( lambdaGrid ) ) deallocate( lambdaGrid )

  end function run

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function sw_futr( w ) result( futr )
    ! Calculate the action spectrum value for skin cancer of albino hairless
    ! mice at a given wavelength according to:  deGRuijl, F.R., H.J.C.M.Steren-
    ! borg, P.D.Forbes, R.E.Davies, C.Colse, G.Kelfkens, H.vanWeelden,
    ! and J.C.van der Leun, Wavelength dependence of skin cancer induction by
    ! ultraviolet irradiation of albino hairless mice, Cancer Research, vol 53,
    ! pp. 53-60, 1993
    ! (Action spectrum for carcinomas)
    !
    ! `Link to the article
    ! <https://aacrjournals.org/cancerres/article/53/1/53/498742/Wavelength-Dependence-of-Skin-Cancer-Induction-by>`_

    real(dk), intent(in)  :: w(:)    ! Wavelength [nm]
    real(dk), allocatable :: futr(:) ! Calculated action spectrum

    real(dk), allocatable :: t1(:), t2(:), t3(:), t4(:), t5(:)
    real(dk), allocatable :: p(:)

    real(dk), parameter :: a1 = -10.91_dk
    real(dk), parameter :: a2 = - 0.86_dk
    real(dk), parameter :: a3 = - 8.60_dk
    real(dk), parameter :: a4 = - 9.36_dk
    real(dk), parameter :: a5 = -13.15_dk

    real(dk), parameter :: x1 = 270._dk
    real(dk), parameter :: x2 = 302._dk
    real(dk), parameter :: x3 = 334._dk
    real(dk), parameter :: x4 = 367._dk
    real(dk), parameter :: x5 = 400._dk

    real(dk), parameter :: b1 =                                               &
        ( x1 - x2 ) * ( x1 - x3 ) * ( x1 - x4 ) * ( x1 - x5 )
    real(dk), parameter :: b2 =                                               &
        ( x2 - x1 ) * ( x2 - x3 ) * ( x2 - x4 ) * ( x2 - x5 )
    real(dk), parameter :: b3 =                                               &
        ( x3 - x1 ) * ( x3 - x2 ) * ( x3 - x4 ) * ( x3 - x5 )
    real(dk), parameter :: b4 =                                               &
        ( x4 - x1 ) * ( x4 - x2 ) * ( x4 - x3 ) * ( x4 - x5 )
    real(dk), parameter :: b5 =                                               &
        ( x5 - x1 ) * ( x5 - x2 ) * ( x5 - x3 ) * ( x5 - x4 )

    real(dk), parameter :: w1 = a1 / b1
    real(dk), parameter :: w2 = a2 / b2
    real(dk), parameter :: w3 = a3 / b3
    real(dk), parameter :: w4 = a4 / b4
    real(dk), parameter :: w5 = a5 / b5

    t1 = ( w - x2 ) * ( w - x3 ) * ( w - x4 ) * ( w - x5 )
    t2 = ( w - x1 ) * ( w - x3 ) * ( w - x4 ) * ( w - x5 )
    t3 = ( w - x1 ) * ( w - x2 ) * ( w - x4 ) * ( w - x5 )
    t4 = ( w - x1 ) * ( w - x2 ) * ( w - x3 ) * ( w - x5 )
    t5 = ( w - x1 ) * ( w - x2 ) * ( w - x3 ) * ( w - x4 )

    p = w1 * t1 + w2 * t2 + w3 * t3 + w4 * t4 + w5 * t5

    futr  = exp( p )

  end function sw_futr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_spectral_weight_scup_mice
