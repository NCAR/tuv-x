! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0

module tuvx_cross_section_ch3coch3_ch3co_ch3
  ! Computes the cross section for acetone

  use tuvx_cross_section,              only : cross_section_t,                &
                                              base_constructor

  implicit none

  private
  public :: cross_section_ch3coch3_ch3co_ch3_t

  !> Calculator for acetone cross_section
  type, extends(cross_section_t) :: cross_section_ch3coch3_ch3co_ch3_t
  contains
    !> Initialize the cross section
    procedure :: calculate
  end type cross_section_ch3coch3_ch3co_ch3_t

  !> Constructor
  interface cross_section_ch3coch3_ch3co_ch3_t
    module procedure constructor
  end interface cross_section_ch3coch3_ch3co_ch3_t

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
    call assert_msg( 577238139,                                               &
                     config%validate( required_keys, optional_keys ),         &
                     "Bad configuration data format for "//                   &
                     "ch3coch3+hv->ch3co+ch3 cross section." )
    allocate ( cross_section_ch3coch3_ch3co_ch3_t :: this )
    call base_constructor( this, config, grid_warehouse, profile_warehouse )

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function calculate( this, grid_warehouse, profile_warehouse, at_mid_point ) &
      result( cross_section )
    ! Calculate the cross section for a given set of environmental conditions
    ! qyacet - q.y. for acetone, based on Blitz et al. (2004)
    ! Compute acetone quantum yields according to the parameterization of:
    ! Blitz, M. A., D. E. Heard, M. J. Pilling, S. R. Arnold,
    ! and M. P. Chipperfield (2004), Pressure and temperature-dependent 
    ! quantum yields for the photodissociation of acetone between 279 
    ! and 327.5 nm, Geophys. Res. Lett., 
    ! 31, L06111, 
    ! `doi:10.1029/2003GL018793 <https://doi.org/10.1029/2003GL018793>`_.

    use musica_assert,                 only : assert_msg
    use musica_constants,              only : dk => musica_dk
    use tuvx_grid,                     only : grid_t
    use tuvx_grid_warehouse,           only : grid_warehouse_t
    use tuvx_profile,                  only : profile_t
    use tuvx_profile_warehouse,        only : profile_warehouse_t

    real(kind=dk), allocatable               :: cross_section(:,:) ! Calculated cross section
    class(cross_section_ch3coch3_ch3co_ch3_t), intent(in) :: this ! This :f:type:`~tuvx_cross_section_ch3coch3_ch3co_ch3/cross_section_ch3coch3_ch3co_ch3_t`
    type(grid_warehouse_t),    intent(inout) :: grid_warehouse ! A :f:type:`~tuvx_grid_warehouse/grid_warehouse_t`
    type(profile_warehouse_t), intent(inout) :: profile_warehouse ! A :f:type:`~tuvx_profile_warehouse/profile_warehouse_t`
    logical, optional,         intent(in)    :: at_mid_point ! Flag indicating whether cross-section data should be at mid-points on the wavelength grid.  If this is false or omitted, cross-section data are calculated at interfaces on the wavelength grid.

    ! Local variables
    character(len=*), parameter :: Iam =                                      &
        'ch3coch3+hv->ch3co+ch3 cross section calculate'
    real(dk), parameter :: rZERO    = 0.0_dk
    real(dk), parameter :: rONE     = 1.0_dk
    integer                   :: nzdim, vertNdx
    real(dk)                      :: Tadj
    real(dk),         allocatable :: modelTemp(:)
    class(grid_t),    pointer     :: zGrid
    class(grid_t),    pointer     :: lambdaGrid
    class(profile_t), pointer     :: mdlTemperature

    zGrid => grid_warehouse%get_grid( this%height_grid_ )
    lambdaGrid => grid_warehouse%get_grid( this%wavelength_grid_ )
    mdlTemperature =>                                                         &
        profile_warehouse%get_profile( this%temperature_profile_ )

    nzdim = zGrid%ncells_ + 1
    if( present(at_mid_point) ) then
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

    call assert_msg(811958314, &
      size( this%cross_section_parms(1)%array, dim = 2 ) == 4, &
      Iam//' array must have 4 parameters')

    associate( coefficient => this%cross_section_parms(1)%array )
      do vertNdx = 1, nzdim
        Tadj = min( 298._dk, max( 235._dk, modelTemp( vertNdx ) ) )
        cross_section(:,vertNdx) = coefficient(:,1)                           &
                          * ( rONE + Tadj * ( coefficient(:,2)                &
                                              + Tadj*(coefficient(:,3)        &
                                              + Tadj*coefficient(:,4) ) ) )
      enddo
    end associate

    cross_section = transpose( cross_section )

    deallocate( zGrid )
    deallocate( lambdaGrid )
    deallocate( mdlTemperature )

  end function calculate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_cross_section_ch3coch3_ch3co_ch3
