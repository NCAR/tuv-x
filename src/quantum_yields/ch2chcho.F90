! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
module tuvx_quantum_yield_ch2chcho
  ! The ch2chcho+hv->prodcuts quantum yield type and related functions

  ! Including musica_config at the module level to avoid an ICE
  ! with Intel 2022.1 compiler
  use musica_config,                   only : config_t
  use tuvx_quantum_yield,              only : quantum_yield_t, base_constructor

  implicit none

  private
  public :: quantum_yield_ch2chcho_t

  type, extends(quantum_yield_t) :: quantum_yield_ch2chcho_t
    ! Calculator for ch2chcho+hv->oh+h quantum yield
  contains
    procedure :: calculate => run
  end type quantum_yield_ch2chcho_t

  interface quantum_yield_ch2chcho_t
    module procedure constructor
  end interface quantum_yield_ch2chcho_t

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

    allocate ( quantum_yield_ch2chcho_t :: this )

    call base_constructor( this, config, grid_warehouse, profile_warehouse )

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate the photorate quantum yield for a given set of environmental
  !! conditions
  function run( this, grid_warehouse, profile_warehouse )                     &
      result( quantum_yield )

    use musica_constants,              only : dk => musica_dk, ik => musica_ik
    use tuvx_grid_warehouse,           only : grid_warehouse_t
    use tuvx_grid,                     only : grid_t
    use tuvx_profile_warehouse,        only : profile_warehouse_t
    use tuvx_profile,                  only : profile_t

    !> Arguments
    class(quantum_yield_ch2chcho_t), intent(in) :: this ! This :f:type:`~tuvx_quantum_yield_ch2chcho/quantum_yield_ch2chcho_t`
    type(grid_warehouse_t),    intent(inout) :: grid_warehouse ! A :f:type:`~tuvx_grid_warehouse/grid_warehouse_t`
    type(profile_warehouse_t), intent(inout) :: profile_warehouse ! A :f:type:`~tuvx_profile_warehouse/profile_warehouse_t`
    real(kind=dk), allocatable               :: quantum_yield(:,:) ! Calculated quantum_yield

    !> Local variables
    character(len=*), parameter :: Iam = 'ch2chcho+hv->products calculate'
    real(dk), parameter :: phiL = .004_dk
    real(dk), parameter :: phiU = .086_dk
    integer    , parameter :: iONE = 1_ik
    real(dk), parameter ::    rZERO = 0.0_dk
    real(dk), parameter ::    rONE  = 1.0_dk

    integer                       :: nzdim, vertNdx
    real(dk),         allocatable :: phi0(:)
    real(dk),         allocatable :: modelDens(:)
    class(grid_t),    pointer     :: zGrid
    class(grid_t),    pointer     :: lambdaGrid
    class(profile_t), pointer     :: mdlDensity

    zGrid => grid_warehouse%get_grid( this%height_grid_ )
    lambdaGrid => grid_warehouse%get_grid( this%wavelength_grid_ )
    mdlDensity => profile_warehouse%get_profile( this%air_profile_ )

    nzdim = zGrid%ncells_ + iONE
    modelDens = mdlDensity%edge_val_

    allocate( quantum_yield(lambdaGrid%ncells_,nzdim) )
    allocate( phi0(lambdaGrid%ncells_) )
    quantum_yield = rZERO

    do vertNdx = iONE,nzdim
      associate( M => modelDens(vertNdx) )
      if( M > 2.6e19_dk ) then
        quantum_yield(:,vertNdx) = phiL
      else
        if( M <= 8.e17_dk ) then
          phi0 = phiU + 1.613e-17_dk*8.e17_dk
        else
          phi0 = phiU + 1.613e-17_dk*M
        endif
        quantum_yield(:,vertNdx) = phiL + rONE/phi0
      endif
      end associate
    enddo

    quantum_yield = transpose( quantum_yield )

    deallocate( zGrid )
    deallocate( lambdaGrid )
    deallocate( mdlDensity )

  end function run

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_quantum_yield_ch2chcho
