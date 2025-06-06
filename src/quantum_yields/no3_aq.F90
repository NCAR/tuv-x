! Copyright (C) 2020-2025 University Corporation for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0

module tuvx_quantum_yield_no3m_aq
  ! The no3-aq+hv->no2(aq)+o- quantum yield type and related functions

  ! Including musica_config at the module level to avoid an ICE
  ! with Intel 2022.1 compiler
  use musica_config,                   only : config_t
  use tuvx_quantum_yield,              only : quantum_yield_t, base_constructor

  implicit none

  private
  public :: quantum_yield_no3m_aq_t

  type, extends(quantum_yield_t) :: quantum_yield_no3m_aq_t
    ! Calculator for no3m(aq)+hv->no2(aq)+o- quantum yield
  contains
    procedure :: calculate => run
  end type quantum_yield_no3m_aq_t

  interface quantum_yield_no3m_aq_t
    module procedure constructor
  end interface quantum_yield_no3m_aq_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function constructor( config, grid_warehouse, profile_warehouse )           &
      result( this )
    ! Constructor

    use tuvx_grid_warehouse,           only : grid_warehouse_t
    use tuvx_profile_warehouse,        only : profile_warehouse_t

    class(quantum_yield_t), pointer :: this ! This :f:type:`~tuvx_quantum_yield/quantum_yield_t` calculator
    type(config_t),            intent(inout) :: config ! Quantum yield configuration data
    type(grid_warehouse_t),    intent(inout) :: grid_warehouse ! A :f:type:`~tuvx_grid_warehouse/grid_warehouse_t`
    type(profile_warehouse_t), intent(inout) :: profile_warehouse ! A :f:type:`~tuvx_profile_warehouse/profile_warehouse_t`

    allocate ( quantum_yield_no3m_aq_t :: this )

    call base_constructor( this, config, grid_warehouse, profile_warehouse )

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function run( this, grid_warehouse, profile_warehouse )                     &
      result( quantum_yield )
    ! Calculate the quantum yield for a given set of environmental conditions

    use musica_constants,              only : dk => musica_dk
    use tuvx_grid,                     only : grid_t
    use tuvx_grid_warehouse,           only : grid_warehouse_t
    use tuvx_profile,                  only : profile_t
    use tuvx_profile_warehouse,        only : profile_warehouse_t

    class(quantum_yield_no3m_aq_t), intent(in) :: this ! This :f:type:`~tuvx_quantum_yield_no3m_aq/quantum_yield_no3m_aq_t`
    type(grid_warehouse_t),    intent(inout) :: grid_warehouse ! A :f:type:`~tuvx_grid_warehouse/grid_warehouse_t`
    type(profile_warehouse_t), intent(inout) :: profile_warehouse ! A :f:type:`~tuvx_profile_warehouse/profile_warehouse_t`
    real(kind=dk), allocatable               :: quantum_yield(:,:) ! Calculated quantum_yield

    ! Local variables
    character(len=*), parameter :: Iam = 'no3-_(aq)+hv->products calculate'
    real(dk), parameter ::    rZERO = 0.0_dk
    real(dk), parameter ::    rONE  = 1.0_dk

    integer                       :: nzdim, vertNdx
    real(dk),         allocatable :: modelTemp(:)
    class(grid_t),    pointer     :: zGrid
    class(grid_t),    pointer     :: lambdaGrid
    class(profile_t), pointer     :: mdlTemperature

    zGrid => grid_warehouse%get_grid( this%height_grid_ )
    lambdaGrid => grid_warehouse%get_grid( this%wavelength_grid_ )
    mdlTemperature =>                                                         &
        profile_warehouse%get_profile( this%temperature_profile_ )

    nzdim = zGrid%ncells_ + 1
    modelTemp = mdlTemperature%edge_val_

    allocate( quantum_yield( lambdaGrid%ncells_, nzdim ) )
    quantum_yield = rZERO

    do vertNdx = 1, nzdim
      quantum_yield( :, vertNdx ) =                                           &
        exp( -2400._dk / modelTemp(vertNdx) + 3.6_dk ) ! Chu & Anastasio, 2003
    enddo

    quantum_yield = transpose( quantum_yield )

    deallocate( zGrid )
    deallocate( lambdaGrid )
    deallocate( mdlTemperature )

  end function run

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_quantum_yield_no3m_aq
