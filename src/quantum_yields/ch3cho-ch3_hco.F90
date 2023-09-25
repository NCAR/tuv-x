! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0

module tuvx_quantum_yield_ch3cho_ch3_hco
  ! The ch3cho+hv->ch3+hco quantum yield type and related functions

  ! Including musica_config at the module level to avoid an ICE
  ! with Intel 2022.1 compiler
  use musica_config,                   only : config_t
  use tuvx_quantum_yield,              only : quantum_yield_t, base_constructor

  implicit none

  private
  public :: quantum_yield_ch3cho_ch3_hco_t

  type, extends(quantum_yield_t) :: quantum_yield_ch3cho_ch3_hco_t
    ! Calculator for ch3cho+hv->ch3+hco quantum yield
  contains
    procedure :: calculate => run
  end type quantum_yield_ch3cho_ch3_hco_t

  interface quantum_yield_ch3cho_ch3_hco_t
    module procedure constructor
  end interface quantum_yield_ch3cho_ch3_hco_t

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

    allocate ( quantum_yield_ch3cho_ch3_hco_t :: this )

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
    use tuvx_profile,                  only : profile_t
    use tuvx_profile_warehouse,        only : profile_warehouse_t

    class(quantum_yield_ch3cho_ch3_hco_t), intent(in) :: this ! This :f:type:`~tuvx_quantum_yield_ch3cho_ch3_hco/quantum_yield_ch3cho_ch3_hco_t`
    type(grid_warehouse_t),    intent(inout) :: grid_warehouse ! A :f:type:`~tuvx_grid_warehouse/grid_warehouse_t`
    type(profile_warehouse_t), intent(inout) :: profile_warehouse ! A :f:type:`~tuvx_profile_warehouse/profile_warehouse_t`
    real(kind=dk), allocatable               :: quantum_yield(:,:) ! Calculated quantum_yield

    ! Local variables
    character(len=*), parameter :: Iam =                                      &
        'ch3cho+hv->ch3_hco quantum yield calculate'
    real(dk), parameter ::    rZERO = 0.0_dk
    real(dk), parameter ::    rONE  = 1.0_dk
    integer                       :: nzdim, vertNdx
    real(dk)                      :: air_dens_factor
    real(dk),         allocatable :: quantum_yield_chnl1(:)
    real(dk),         allocatable :: quantum_yield_chnl2(:)
    real(dk),         allocatable :: quantum_yield_wrk(:)
    real(dk),         allocatable :: modelDens(:)
    class(grid_t),    pointer     :: zGrid
    class(grid_t),    pointer     :: lambdaGrid
    class(profile_t), pointer     :: mdlDensity

    zGrid => grid_warehouse%get_grid( this%height_grid_ )
    lambdaGrid => grid_warehouse%get_grid( this%wavelength_grid_ )
    mdlDensity => profile_warehouse%get_profile( this%air_profile_ )

    nzdim = zGrid%ncells_ + 1
    modelDens = mdlDensity%edge_val_

    allocate( quantum_yield( lambdaGrid%ncells_, nzdim ) )
    allocate( quantum_yield_wrk( lambdaGrid%ncells_ ) )
    quantum_yield = rZERO

    quantum_yield_chnl1 = this%quantum_yield_parms(1)%array(:,2)
    quantum_yield_chnl2 = rONE - this%quantum_yield_parms(1)%array(:,1)
    where( quantum_yield_chnl1 > rZERO )
      quantum_yield_wrk = quantum_yield_chnl2 / quantum_yield_chnl1 - rONE
    elsewhere
      quantum_yield_wrk = rZERO
    endwhere
    do vertNdx = 1,nzdim
      air_dens_factor = modelDens( vertNdx ) / 2.465e19_dk
      quantum_yield( :, vertNdx ) =                                           &
                      quantum_yield_chnl1 * ( rONE + quantum_yield_wrk )      &
                      / ( rONE + quantum_yield_wrk * air_dens_factor )
      quantum_yield( :, vertNdx ) =                                           &
          min( rONE, max( rZERO, quantum_yield( :, vertNdx ) ) )
    enddo

    quantum_yield = transpose( quantum_yield )

    deallocate( zGrid )
    deallocate( lambdaGrid )
    deallocate( mdlDensity )

  end function run

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_quantum_yield_ch3cho_ch3_hco
