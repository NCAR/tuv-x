! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0

module tuvx_cross_section_hobr_oh_br
! Calculates the cross section for hypobromous acid

  use tuvx_cross_section,              only : cross_section_t,                &
                                              base_constructor

  implicit none

  private
  public :: cross_section_hobr_oh_br_t

  !> Calculator for hobr-oh_br cross section
  type, extends(cross_section_t) :: cross_section_hobr_oh_br_t
  contains
    !> Calculate the cross section
    procedure :: calculate
  end type cross_section_hobr_oh_br_t

  !> Constructor
  interface cross_section_hobr_oh_br_t
    module procedure constructor
  end interface cross_section_hobr_oh_br_t

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
    call assert_msg( 916010327,                                               &
                     config%validate( required_keys, optional_keys ),         &
                     "Bad configuration data format for "//                   &
                     "hobr_oh_br cross section." )
    allocate( cross_section_hobr_oh_br_t :: this )
    call base_constructor( this, config, grid_warehouse, profile_warehouse )

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function calculate( this, grid_warehouse, profile_warehouse, at_mid_point ) &
      result( cross_section )
    ! Calculate the cross section for a given set of environmental conditions

    use musica_constants,              only : dk => musica_dk
    use tuvx_grid,                     only : grid_t
    use tuvx_grid_warehouse,           only : grid_warehouse_t
    use tuvx_profile_warehouse,        only : profile_warehouse_t

    real(kind=dk), allocatable                       :: cross_section(:,:) ! Calculated cross section
    class(cross_section_hobr_oh_br_t), intent(in)    :: this ! A :f:type:`~tuvx_cross_section_hobr_oh_br/cross_section_hobr_oh_br_t`
    type(grid_warehouse_t),    intent(inout) :: grid_warehouse ! A :f:type:`~tuvx_grid_warehouse/grid_warehouse_t`
    type(profile_warehouse_t), intent(inout) :: profile_warehouse ! A :f:type:`~tuvx_profile_warehouse/profile_warehouse_t`
    logical, optional,           intent(in)    :: at_mid_point ! Flag indicating whether cross-section data should be at mid-points on the wavelength grid.  If this is false or omitted, cross-section data are calculated at interfaces on the wavelength grid.

    ! Local variables
    character(len=*), parameter :: Iam =                                      &
        'hobr_oh_br cross section calculate: '
    real(dk), parameter :: rZERO = 0.0_dk
    real(dk), parameter :: a = -2.359E-3_dk
    real(dk), parameter :: b = 1.2478_dk
    real(dk), parameter :: c = -210.4_dk

    integer                    :: nzdim
    integer                    :: vertNdx
    real(dk),      allocatable :: wrkCrossSection(:)
    class(grid_t), pointer     :: zGrid
    class(grid_t), pointer     :: lambdaGrid

    zGrid => grid_warehouse%get_grid( this%height_grid_ )
    lambdaGrid => grid_warehouse%get_grid( this%wavelength_grid_ )

    nzdim = zGrid%ncells_ + 1
    if( present( at_mid_point ) ) then
      if( at_mid_point ) then
        nzdim = nzdim - 1
      endif
    endif

    allocate( cross_section( lambdaGrid%ncells_,nzdim ) )
    allocate( wrkCrossSection( lambdaGrid%ncells_ ) )

    associate( wc => lambdaGrid%mid_ )
    where( wc >= 250._dk .and. wc <= 550._dk )
      wrkCrossSection = &
               24.77_dk * exp( -109.80_dk * ( log( 284.01_dk / wc) )**2 )     &
             + 12.22_dk * exp(  -93.63_dk * ( log( 350.57_dk / wc) )**2 )     &
             + 2.283_dk * exp(- 242.40_dk * ( log( 457.38_dk / wc) )**2 )
      wrkCrossSection = wrkCrossSection * 1.e-20_dk
    elsewhere
      wrkCrossSection = rZERO
    endwhere
    end associate

    do vertNdx = 1,nzdim
      cross_section( :, vertNdx ) = wrkCrossSection
    enddo

    cross_section = transpose( cross_section )

    deallocate( zGrid )
    deallocate( lambdaGrid )

  end function calculate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_cross_section_hobr_oh_br
