! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0

module tuvx_cross_section_h2o2_oh_oh
! Calculates the cross section for hydrogen peroxide

  use tuvx_cross_section,              only : cross_section_t,                &
                                              base_constructor

  implicit none

  private
  public :: cross_section_h2o2_oh_oh_t

  !> Calculator for base_cross_section
  type, extends(cross_section_t) :: cross_section_h2o2_oh_oh_t
  contains
    !> Calculate the cross section
    procedure :: calculate
  end type cross_section_h2o2_oh_oh_t

  !> Constructor
  interface cross_section_h2o2_oh_oh_t
    module procedure constructor
  end interface cross_section_h2o2_oh_oh_t

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
    call assert_msg( 969725098,                                               &
                     config%validate( required_keys, optional_keys ),         &
                     "Bad configuration data format for "//                   &
                     "h2o2+hv->oh+oh cross section." )
    allocate( cross_section_h2o2_oh_oh_t :: this )
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

    real(kind=dk), allocatable                       :: cross_section(:,:) ! Calculated cross section
    class(cross_section_h2o2_oh_oh_t), intent(in)    :: this ! A :f:type:`~tuvx_cross_section_h2o2_oh_oh/cross_section_h2o2_oh_oh_t`
    type(grid_warehouse_t),    intent(inout) :: grid_warehouse ! A :f:type:`~tuvx_grid_warehouse/grid_warehouse_t`
    type(profile_warehouse_t), intent(inout) :: profile_warehouse ! A :f:type:`~tuvx_profile_warehouse/profile_warehouse_t`
    logical, optional,           intent(in)    :: at_mid_point ! Flag indicating whether cross-section data should be at mid-points on the wavelength grid.  If this is false or omitted, cross-section data are calculated at interfaces on the wavelength grid.

    ! local variables
    real(dk), parameter ::  rONE = 1.0_dk
    real(dk), parameter ::  A0 = 6.4761E+04_dk
    real(dk), parameter ::  A1 = -9.2170972E+02_dk
    real(dk), parameter ::  A2 = 4.535649_dk
    real(dk), parameter ::  A3 = -4.4589016E-03_dk
    real(dk), parameter ::  A4 = -4.035101E-05_dk
    real(dk), parameter ::  A5 = 1.6878206E-07_dk
    real(dk), parameter ::  A6 = -2.652014E-10_dk
    real(dk), parameter ::  A7 = 1.5534675E-13_dk

    real(dk), parameter ::  B0 = 6.8123E+03_dk
    real(dk), parameter ::  B1 = -5.1351E+01_dk
    real(dk), parameter ::  B2 = 1.1522E-01_dk
    real(dk), parameter ::  B3 = -3.0493E-05_dk
    real(dk), parameter ::  B4 = -1.0924E-07_dk

    character(len=*), parameter :: Iam =                                      &
        'h2o2+hv->oh+oh cross section calculate'
    integer    :: vertNdx, wNdx
    real(dk)       :: lambda, sumA, sumB, t, chi
    class(grid_t),    pointer :: zGrid
    class(grid_t),    pointer :: lambdaGrid
    class(profile_t), pointer :: temperature

    zGrid => grid_warehouse%get_grid( this%height_grid_ )
    lambdaGrid => grid_warehouse%get_grid( this%wavelength_grid_ )
    temperature =>                                                            &
        profile_warehouse%get_profile( this%temperature_profile_ )

    allocate( cross_section( lambdaGrid%ncells_, zGrid%ncells_ + 1 ) )

    associate( wl => lambdaGrid%edge_, wc => lambdaGrid%mid_ )
    do vertNdx = 1, zGrid%ncells_ + 1
      do wNdx = 1, lambdaGrid%ncells_
        ! Parameterization (JPL94)
        ! Range 260-350 nm; 200-400 K
        if( wl( wNdx ) >= 260._dk .and. wl( wNdx ) < 350._dk ) then
           lambda = wc( wNdx )
           sumA = ( ( ( ( ( ( A7 * lambda + A6 ) * lambda + A5 ) * lambda     &
                  + A4 ) * lambda + A3 ) * lambda + A2 ) * lambda             &
                  + A1 ) * lambda + A0
           sumB = ( ( ( B4 * lambda + B3 ) * lambda + B2 ) * lambda + B1 )    &
                  * lambda + B0
           t = min( max( temperature%edge_val_( vertNdx ), 200._dk ), 400._dk )
           chi = rONE / ( rONE + exp( -1265._dk / t ) )
           cross_section( wNdx, vertNdx ) =                                   &
               ( chi * sumA + ( rONE - chi ) * sumB ) * 1.E-21_dk
         else
           cross_section( wNdx, vertNdx ) =                                   &
               this%cross_section_parms(1)%array( wNdx, 1 )
         endif
      enddo
    enddo
    end associate

    cross_section = transpose( cross_section )

    deallocate( zGrid )
    deallocate( lambdaGrid )
    deallocate( temperature )

  end function calculate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_cross_section_h2o2_oh_oh
