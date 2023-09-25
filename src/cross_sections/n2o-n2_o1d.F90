! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0

module tuvx_cross_section_n2o_n2_o1d
! Calculates the cross section for nitrogen oxide

  use tuvx_cross_section,              only : cross_section_t,                &
                                              base_constructor

  implicit none

  private
  public :: cross_section_n2o_n2_o1d_t

  !> Calculator for base_cross_section
  type, extends(cross_section_t) :: cross_section_n2o_n2_o1d_t
  contains
    !> Calculate the cross section
    procedure :: calculate
  end type cross_section_n2o_n2_o1d_t

  !> Constructor
  interface cross_section_n2o_n2_o1d_t
    module procedure constructor
  end interface cross_section_n2o_n2_o1d_t

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
    call assert_msg( 210014560,                                               &
                     config%validate( required_keys, optional_keys ),         &
                     "Bad configuration data format for "//                   &
                     "n2o_n2+o1d cross section." )
    allocate( cross_section_n2o_n2_o1d_t :: this )
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
    class(cross_section_n2o_n2_o1d_t), intent(in)    :: this ! A :f:type:`~tuvx_cross_section_n2o_n2_o1d/cross_section_n2o_n2_o1d_t`
    type(grid_warehouse_t),    intent(inout) :: grid_warehouse ! A :f:type:`~tuvx_grid_warehouse/grid_warehouse_t`
    type(profile_warehouse_t), intent(inout) :: profile_warehouse ! A :f:type:`~tuvx_profile_warehouse/profile_warehouse_t`
    logical, optional,           intent(in)    :: at_mid_point ! Flag indicating whether cross-section data should be at mid-points on the wavelength grid.  If this is false or omitted, cross-section data are calculated at interfaces on the wavelength grid.

    !> Local variables
    character(len=*), parameter :: Iam =                                      &
        'n2o_n2+o1d cross section calculate'

    real(dk), parameter :: rZERO = 0._dk
    real(dk), parameter :: A0 = 68.21023_dk
    real(dk), parameter :: A1 = -4.071805_dk
    real(dk), parameter :: A2 = 4.301146E-02_dk
    real(dk), parameter :: A3 = -1.777846E-04_dk
    real(dk), parameter :: A4 = 2.520672E-07_dk

    real(dk), parameter :: B0 = 123.4014_dk
    real(dk), parameter :: B1 = -2.116255_dk
    real(dk), parameter :: B2 = 1.111572E-02_dk
    real(dk), parameter :: B3 = -1.881058E-05_dk

    real(dk), parameter :: Tlower = 173._dk
    real(dk), parameter :: Tupper = 240._dk
    real(dk), parameter :: Tfloor = 194._dk
    real(dk), parameter :: Tceil  = 320._dk
    real(dk), parameter :: Thold  = 300._dk

    integer :: lambdaNdx, nzdim, vertNdx
    real(dk)    :: lambda, Tadj, A, B
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

    !*** quantum yield of N(4s) and NO(2Pi) is less than 1% (Greenblatt and
    !*** Ravishankara), so quantum yield of O(1D) is assumed to be unity
    do vertNdx = 1, nzdim
      Tadj = max( Tfloor, min( modelTemp( vertNdx ), Tceil ) )
      do lambdaNdx = 1, lambdaGrid%ncells_
        lambda = lambdaGrid%mid_( lambdaNdx )
        if( lambda >= Tlower .and. lambda <= Tupper) then
          A = ( ( ( A4 * lambda + A3 ) * lambda + A2 ) * lambda + A1 )        &
              * lambda+A0
          B = ( ( B3 * lambda + B2 ) * lambda + B1 ) * lambda + B0
          B = ( Tadj - Thold ) * exp( B )
          cross_section( lambdaNdx, vertNdx ) = exp( A + B )
        endif
      enddo
    enddo

    cross_section = transpose( cross_section )

    deallocate( zGrid )
    deallocate( lambdaGrid )
    deallocate( mdlTemperature )

  end function calculate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_cross_section_n2o_n2_o1d
