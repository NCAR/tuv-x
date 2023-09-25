! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0

module tuvx_quantum_yield_clono2_clo_no2
  ! The clono2+hv->clo+no2 quantum yield type and related functions

  ! Including musica_config at the module level to avoid an ICE
  ! with Intel 2022.1 compiler
  use musica_config,                   only : config_t
  use tuvx_quantum_yield,              only : quantum_yield_t, base_constructor

  implicit none

  private
  public :: quantum_yield_clono2_clo_no2_t

  type, extends(quantum_yield_t) :: quantum_yield_clono2_clo_no2_t
    ! Calculator for clono2+hv->clo+no2 quantum yield
  contains
    procedure :: calculate => run
  end type quantum_yield_clono2_clo_no2_t

  interface quantum_yield_clono2_clo_no2_t
    module procedure constructor
  end interface quantum_yield_clono2_clo_no2_t

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

    allocate ( quantum_yield_clono2_clo_no2_t :: this )

    call base_constructor( this, config, grid_warehouse, profile_warehouse )

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function run( this, grid_warehouse, profile_warehouse )                     &
      result( quantum_yield )
    ! Calculate the photorate quantum yield for a given set of environmental
    ! conditions

    use musica_constants,              only : dk => musica_dk
    use tuvx_grid,                     only : grid_t
    use tuvx_grid_warehouse,           only : grid_warehouse_t
    use tuvx_profile_warehouse,        only : profile_warehouse_t

    class(quantum_yield_clono2_clo_no2_t), intent(in) :: this ! This :f:type:`~tuvx_quantum_yield_clono2_clo_no2/quantum_yield_clono2_clo_no2_t`
    type(grid_warehouse_t),    intent(inout) :: grid_warehouse ! A :f:type:`~tuvx_grid_warehouse/grid_warehouse_t`
    type(profile_warehouse_t), intent(inout) :: profile_warehouse ! A :f:type:`~tuvx_profile_warehouse/profile_warehouse_t`
    real(kind=dk), allocatable               :: quantum_yield(:,:) ! Calculated quantum_yield

    ! Local variables
    character(len=*), parameter :: Iam = 'clono2+hv->clo+no2 calculate'
    real(dk), parameter :: rZERO = 0.0_dk
    real(dk), parameter :: rONE  = 1.0_dk

    integer                    :: vertNdx, lambdaNdx
    real(dk)                   :: lambda
    real(dk),      allocatable :: wrkQuantumYield(:)
    class(grid_t), pointer     :: lambdaGrid
    class(grid_t), pointer     :: zGrid

    zGrid => grid_warehouse%get_grid( this%height_grid_ )
    lambdaGrid => grid_warehouse%get_grid( this%wavelength_grid_ )

    allocate( wrkQuantumYield( lambdaGrid%ncells_ ) )
    allocate( quantum_yield( lambdaGrid%ncells_, zGrid%ncells_ + 1 ) )

    do lambdaNdx = 1, lambdaGrid%ncells_
      lambda = lambdaGrid%mid_( lambdaNdx )
      if( lambda < 308._dk ) then
        wrkQuantumYield( lambdaNdx ) = .6_dk
      elseif( 308._dk <= lambda .and. lambda <= 364._dk ) then
        wrkQuantumYield( lambdaNdx ) = 7.143e-3_dk * lambda - 1.6_dk
      else
        wrkQuantumYield( lambdaNdx ) = rONE
      endif
    enddo

    do vertNdx = 1, zGrid%ncells_ + 1
      quantum_yield( :, vertNdx ) = rONE - wrkQuantumYield
    enddo

    quantum_yield = transpose( quantum_yield )

    deallocate( zGrid )
    deallocate( lambdaGrid )

  end function run

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_quantum_yield_clono2_clo_no2
