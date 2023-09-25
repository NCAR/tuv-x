! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0

module tuvx_quantum_yield_ch2o_h2_co
  ! The ch2o+hv->h2+co quantum yield type and related functions

  ! Including musica_config at the module level to avoid an ICE
  ! with Intel 2022.1 compiler
  use musica_config,                   only : config_t
  use tuvx_quantum_yield,              only : quantum_yield_t, base_constructor

  implicit none

  private
  public :: quantum_yield_ch2o_h2_co_t


  type, extends(quantum_yield_t) :: quantum_yield_ch2o_h2_co_t
    ! Calculator for ch2o+hv->h2+co quantum yield
  contains
    procedure :: calculate => run
  end type quantum_yield_ch2o_h2_co_t

  interface quantum_yield_ch2o_h2_co_t
    module procedure constructor
  end interface quantum_yield_ch2o_h2_co_t

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

    allocate ( quantum_yield_ch2o_h2_co_t :: this )

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

    class(quantum_yield_ch2o_h2_co_t), intent(in) :: this ! This :f:type:`~tuvx_quantum_yield_ch2o_h2_co/quantum_yield_ch2o_h2_co_t`
    type(grid_warehouse_t),    intent(inout) :: grid_warehouse ! A :f:type:`~tuvx_grid_warehouse/grid_warehouse_t`
    type(profile_warehouse_t), intent(inout) :: profile_warehouse ! A :f:type:`~tuvx_profile_warehouse/profile_warehouse_t`
    real(kind=dk), allocatable               :: quantum_yield(:,:) ! Calculated quantum_yield

    ! Local variables
    character(len=*), parameter :: Iam =                                      &
        'ch2o+hv->h2_co quantum yield calculate'
    real(dk), parameter ::    rZERO = 0.0_dk
    real(dk), parameter ::    rONE  = 1.0_dk
    real(dk), parameter  :: lambdaL = 330._dk
    real(dk), parameter  :: lambdaU = 360._dk

    integer                       :: nzdim, vertNdx
    real(dk)                      :: Tfactor
    real(dk),         allocatable :: quantum_yield_tmp(:)
    real(dk),         allocatable :: quantum_yield_wrk(:)
    real(dk),         allocatable :: modelTemp(:), modelDens(:)
    class(grid_t),    pointer     :: zGrid
    class(grid_t),    pointer     :: lambdaGrid
    class(profile_t), pointer     :: mdlTemperature
    class(profile_t), pointer     :: mdlDensity

    zGrid => grid_warehouse%get_grid( this%height_grid_ )
    lambdaGrid => grid_warehouse%get_grid( this%wavelength_grid_ )
    mdlTemperature =>                                                         &
        profile_warehouse%get_profile( this%temperature_profile_ )
    mdlDensity => profile_warehouse%get_profile( this%air_profile_ )

    nzdim = zGrid%ncells_ + 1
    modelTemp = mdlTemperature%edge_val_
    modelDens = mdlDensity%edge_val_

    allocate( quantum_yield( lambdaGrid%ncells_, nzdim ) )
    quantum_yield = rZERO

    associate( quantum_yield_chnl1 => this%quantum_yield_parms(1)%array(:,1), &
               quantum_yield_chnl2 => this%quantum_yield_parms(1)%array(:,2) )
    quantum_yield_tmp   = rONE - quantum_yield_chnl1
    allocate( quantum_yield_wrk( lambdaGrid%ncells_ ) )
    do vertNdx = 1, nzdim
      Tfactor = ( 300._dk - modelTemp( vertNdx ) ) / 80._dk
      where( lambdaGrid%mid_ >= lambdaL .and. lambdaGrid%mid_ < lambdaU &
                                        .and. quantum_yield_chnl2 > rZERO )
        quantum_yield_wrk = ( rONE -                                           &
                               ( quantum_yield_chnl1 + quantum_yield_chnl2 ) ) &
                      / ( 2.45e19_dk * quantum_yield_chnl2 * quantum_yield_tmp )
        quantum_yield_wrk = quantum_yield_wrk * ( rONE                         &
                          + .05_dk * ( lambdaGrid%mid_ - 329._dk ) * Tfactor )
        quantum_yield(:,vertNdx) = rONE / ( rONE / quantum_yield_tmp +         &
                                   quantum_yield_wrk * modelDens( vertNdx ) )
      elsewhere
        quantum_yield( :, vertNdx ) = quantum_yield_chnl2
      endwhere
    enddo
    end associate

    quantum_yield = transpose( quantum_yield )

    deallocate( zGrid )
    deallocate( lambdaGrid )
    deallocate( mdlTemperature )
    deallocate( mdlDensity )

  end function run

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_quantum_yield_ch2o_h2_co
