! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0

module tuvx_cross_section_cl2_cl_cl
! Calculates the cross section for gaseous chlorine

  use tuvx_cross_section,              only : cross_section_t,                &
                                              base_constructor

  implicit none

  private
  public :: cross_section_cl2_cl_cl_t

  !> Calculator for base_cross_section
  type, extends(cross_section_t) :: cross_section_cl2_cl_cl_t
  contains
    !> Calculate the cross section
    procedure :: calculate
  end type cross_section_cl2_cl_cl_t

  !> Constructor
  interface cross_section_cl2_cl_cl_t
    module procedure constructor
  end interface cross_section_cl2_cl_cl_t

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

    type(string_t) :: required_keys(1), optional_keys(1)

    required_keys(1) = "type"
    optional_keys(1) = "name"
    call assert_msg( 864124299,                                               &
                     config%validate( required_keys, optional_keys ),         &
                     "Bad configuration data format for "//                   &
                     "cl2+hv->cl+cl cross section." )
    allocate( cross_section_cl2_cl_cl_t :: this )
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

    real(kind=dk), allocatable                      :: cross_section(:,:) ! Calculated cross section
    class(cross_section_cl2_cl_cl_t), intent(in)    :: this ! A :f:type:`~tuvx_cross_section_cl2_cl_cl/cross_section_cl2_cl_cl_t`
    type(grid_warehouse_t),    intent(inout) :: grid_warehouse ! A :f:type:`~tuvx_grid_warehouse/grid_warehouse_t`
    type(profile_warehouse_t), intent(inout) :: profile_warehouse ! A :f:type:`~tuvx_profile_warehouse/profile_warehouse_t`
    logical, optional,           intent(in)    :: at_mid_point ! Flag indicating whether cross-section data should be at mid-points on the wavelength grid.  If this is false or omitted, cross-section data are calculated at interfaces on the wavelength grid.

    ! Local variables
    character(len=*), parameter :: Iam =                                      &
        'cl2+hv->cl+cl cross section calculate'
    real(dk), parameter     :: rONE = 1.0_dk
    integer :: lambdaNdx, vertNdx, nzdim
    real(dk)    :: aa, bb, bbsq, alpha, ex1, ex2
    real(dk),         allocatable :: modelTemp(:)
    class(grid_t),    pointer     :: lambdaGrid
    class(grid_t),    pointer     :: zGrid
    class(profile_t), pointer     :: temperature

    zGrid => grid_warehouse%get_grid( this%height_grid_ )
    lambdaGrid => grid_warehouse%get_grid( this%wavelength_grid_ )
    temperature =>                                                            &
        profile_warehouse%get_profile( this%temperature_profile_ )

    nzdim = zGrid%ncells_ + 1
    if( present( at_mid_point ) ) then
      if( at_mid_point ) then
        nzdim = nzdim - 1
        modelTemp = temperature%mid_val_
      else
        modelTemp = temperature%edge_val_
      endif
    else
      modelTemp = temperature%edge_val_
    endif

    allocate( cross_section( lambdaGrid%ncells_, zGrid%ncells_ + 1 ) )

    associate( wc => lambdaGrid%mid_ )
    do vertNdx = 1, nzdim
      aa    = 402.7_dk / modelTemp( vertNdx )
      bb    = exp( aa )
      bbsq  = bb * bb
      alpha = ( bbsq - rONE ) / ( bbsq + rONE )
      do lambdaNdx = 1, lambdaGrid%ncells_
        ex1 = 27.3_dk * exp( -99.0_dk * alpha                                 &
                      * ( log( 329.5_dk / wc( lambdaNdx ) ) )**2 )
        ex2 =  .932_dk * exp( -91.5_dk * alpha                                &
                       * ( log( 406.5_dk / wc( lambdaNdx ) ) )**2 )
        cross_section( lambdaNdx, vertNdx ) =                                 &
            1.e-20_dk * sqrt( alpha ) * ( ex1 + ex2 )
      enddo
    enddo
    end associate

    cross_section = transpose( cross_section )

    deallocate( zGrid )
    deallocate( lambdaGrid )
    deallocate( temperature )

  end function calculate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_cross_section_cl2_cl_cl
