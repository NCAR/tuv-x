! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0

module tuvx_cross_section_hcfc
! Calculates the cross section for an HCFC

  use tuvx_cross_section,              only : cross_section_t,                &
                                              base_constructor

  implicit none

  private
  public :: cross_section_hcfc_t

  !> Calculator for acetone cross_section
  type, extends(cross_section_t) :: cross_section_hcfc_t
  contains
    !> Initialize the cross section
    procedure :: calculate
  end type cross_section_hcfc_t

  !> Constructor
  interface cross_section_hcfc_t
    module procedure constructor
  end interface cross_section_hcfc_t

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

    type(string_t) :: required_keys(2), optional_keys(1)

    required_keys(1) = "type"
    required_keys(2) = "netcdf files"
    optional_keys(1) = "name"
    call assert_msg( 369759228,                                               &
                     config%validate( required_keys, optional_keys ),         &
                     "Bad configuration data format for "//                   &
                     "cf3chcl2+hv->products cross section." )
    allocate( cross_section_hcfc_t :: this )
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
    class(cross_section_hcfc_t), intent(in)    :: this ! A :f:type:`~tuvx_cross_section_hcfc/cross_section_hcfc_t`
    type(grid_warehouse_t),    intent(inout) :: grid_warehouse ! A :f:type:`~tuvx_grid_warehouse/grid_warehouse_t`
    type(profile_warehouse_t), intent(inout) :: profile_warehouse ! A :f:type:`~tuvx_profile_warehouse/profile_warehouse_t`
    logical, optional,           intent(in)    :: at_mid_point ! Flag indicating whether cross-section data should be at mid-points on the wavelength grid.  If this is false or omitted, cross-section data are calculated at interfaces on the wavelength grid.

    ! Local variables
    character(len=*), parameter :: Iam =                                      &
        'cf3chcl2+hv->products cross section calculate'
    real(dk), parameter    :: rZERO = 0.0_dk
    real(dk), parameter    :: LBar  = 206.214_dk
    integer                   :: nzdim, vertNdx
    integer                   :: lambdaNdx, polyNdx
    real(dk)                      :: Tadj, sigma, uLambda
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

    uLambda = this%cross_section_parms(1)%temperature(2)
vert_loop:                                                                    &
    do vertNdx = 1, nzdim
      Tadj = min( 295._dk, max( 203._dk,modelTemp( vertNdx ) ) )              &
             - this%cross_section_parms(1)%temperature(1)
lambda_loop:                                                                  &
      do lambdaNdx = 1, lambdaGrid%ncells_
        if( lambdaGrid%mid_( lambdaNdx ) >= 190._dk                           &
            .and. lambdaGrid%mid_( lambdaNdx ) <= uLambda ) then
          sigma = rZERO
          associate( coefficient => this%cross_section_parms(1)%array )
          do polyNdx = 1, size( coefficient, dim = 1 )
            sigma = sigma                                                     &
                  + ( coefficient( polyNdx, 1 )                               &
                      + Tadj * ( coefficient( polyNdx, 2 )                    &
                                 + Tadj * coefficient( polyNdx, 3 ) ) )       &
                    * ( lambdaGrid%mid_( lambdaNdx ) - LBar )**( polyNdx - 1 )
          enddo
          end associate
          sigma = exp( sigma )
        else
          sigma = rZERO
        endif
        cross_section( lambdaNdx, vertNdx ) = sigma
      enddo lambda_loop
    enddo vert_loop

    cross_section = transpose( cross_section )

    deallocate( zGrid )
    deallocate( lambdaGrid )
    deallocate( mdlTemperature )

  end function calculate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_cross_section_hcfc
