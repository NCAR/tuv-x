! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0

module tuvx_cross_section_ccl4
! Calculates the cross section for CCl4

  use tuvx_cross_section,              only : cross_section_t,                &
                                              base_constructor

  implicit none

  private
  public :: cross_section_ccl4_t

  !> Calculator for ccl4 cross section
  type, extends(cross_section_t) :: cross_section_ccl4_t
  contains
    !> Calculate the cross section
    procedure :: calculate
  end type cross_section_ccl4_t

  !> Constructor
  interface cross_section_ccl4_t
    module procedure constructor
  end interface cross_section_ccl4_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function constructor( config, grid_warehouse, profile_warehouse )           &
      result ( this )
    ! Initialize the cross section

    use musica_assert,                 only : assert_msg
    use musica_config,                 only : config_t
    use musica_string,                 only : string_t
    use tuvx_grid_warehouse,           only : grid_warehouse_t
    use tuvx_profile_warehouse,        only : profile_warehouse_t

    class(cross_section_t),    pointer       :: this ! This :f:type:`~tuvx_cross_section/cross_section_t`
    type(config_t),            intent(inout) :: config ! Cross section configuration data
    type(grid_warehouse_t),    intent(inout) :: grid_warehouse ! A :f:type:`~tuvx_grid_warehouse/grid_warehouse_t`
    type(profile_warehouse_t), intent(inout) :: profile_warehouse ! A :f:type:`~tuvx_profile_warehouse/profile_warehouse_t`

    type(string_t) :: required_keys(2), optional_keys(3)

    required_keys(1) = "type"
    required_keys(2) = "netcdf files"
    optional_keys(1) = "lower extrapolation"
    optional_keys(2) = "upper extrapolation"
    optional_keys(3) = "name"
    call assert_msg( 399197573,                                               &
                     config%validate( required_keys, optional_keys ),         &
                     "Bad configuration data format for "//                   &
                     "ccl4 cross section." )
    allocate (cross_section_ccl4_t :: this )
    call base_constructor( this, config, grid_warehouse, profile_warehouse )

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function calculate( this, grid_warehouse, profile_warehouse, at_mid_point ) &
      result( cross_section )
    ! Calculate the cross section for a given set of environmental conditions

    use musica_constants,              only : dk => musica_dk
    use tuvx_grid,                     only : grid_t
    use tuvx_grid_warehouse,           only : grid_warehouse_t
    use tuvx_profile,                  only : profile_t
    use tuvx_profile_warehouse,        only : profile_warehouse_t

    real(kind=dk), allocatable                 :: cross_section(:,:) ! Calculated cross section
    class(cross_section_ccl4_t), intent(in)    :: this ! A :f:type:`~tuvx_cross_section_ccl4/cross_section_ccl4_t`
    type(grid_warehouse_t),    intent(inout) :: grid_warehouse ! A :f:type:`~tuvx_grid_warehouse/grid_warehouse_t`
    type(profile_warehouse_t), intent(inout) :: profile_warehouse ! A :f:type:`~tuvx_profile_warehouse/profile_warehouse_t`
    logical, optional,           intent(in)    :: at_mid_point ! Flag indicating whether cross-section data should be at mid-points on the wavelength grid.  If this is false or omitted, cross-section data are calculated at interfaces on the wavelength grid.

    ! Local variables
    character(len=*), parameter :: Iam = 'ccl4 cross section run: '
    real(dk), parameter :: rZERO = 0.0_dk
    real(dk), parameter :: b0 = 1.0739_dk
    real(dk), parameter :: b1 = -1.6275e-2_dk
    real(dk), parameter :: b2 = 8.8141e-5_dk
    real(dk), parameter :: b3 = -1.9811e-7_dk
    real(dk), parameter :: b4 = 1.5022e-10_dk

    real(dk)    :: Temp, Wpoly, w1
    integer :: lambdaNdx, vertNdx, nzdim
    real(dk),         allocatable :: modelTemp(:)
    class(grid_t),    pointer     :: zGrid
    class(grid_t),    pointer     :: lambdaGrid
    class(profile_t), pointer     :: mdlTemperature

    zGrid => grid_warehouse%get_grid( this%height_grid_ )
    lambdaGrid => grid_warehouse%get_grid( this%wavelength_grid_ )
    mdlTemperature =>                                                         &
        profile_warehouse%get_profile( this%temperature_profile_ )

    nzdim = zGrid%ncells_ + 1
    if( present( at_mid_point ) ) then
      if( at_mid_point ) then
        nzdim = nzdim - 1
        modelTemp = mdlTemperature%mid_val_
      else
        modelTemp = mdlTemperature%edge_val_
      endif
    else
      modelTemp = mdlTemperature%edge_val_
    endif

    allocate( cross_section( lambdaGrid%ncells_, nzdim ) )
    cross_section = rZERO

    do vertNdx = 1,nzdim
      Temp = max( min( 300._dk, modelTemp( vertNdx ) ), 210._dk )
      Temp = Temp - 295._dk
      do lambdaNdx = 1, lambdaGrid%ncells_
        w1 = lambdaGrid%mid_( lambdaNdx )
        if( w1 > 194._dk .and. w1 < 250._dk ) then
          Wpoly = b0 + w1 * ( b1 + w1 * ( b2 + w1 * ( b3 + b4 * w1 ) ) )
          cross_section( lambdaNdx, vertNdx ) =                               &
              this%cross_section_parms(1)%array( lambdaNdx, 1 )               &
              * 10._dk**( Wpoly * Temp )
        else
          cross_section( lambdaNdx, vertNdx ) =                               &
              this%cross_section_parms(1)%array(lambdaNdx,1)
        endif
      enddo
    enddo

    cross_section = transpose( cross_section )

    deallocate( zGrid )
    deallocate( lambdaGrid )
    deallocate( mdlTemperature )

  end function calculate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_cross_section_ccl4
