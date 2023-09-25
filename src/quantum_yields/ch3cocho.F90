! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0

module tuvx_quantum_yield_ch3cocho_ch3co_hco
  ! The ch3cocho+hv->ch3co+hco quantum yield type and related functions

  ! Including musica_config at the module level to avoid an ICE
  ! with Intel 2022.1 compiler
  use musica_config,                   only : config_t
  use tuvx_quantum_yield,              only : quantum_yield_t, base_constructor

  implicit none

  private
  public :: quantum_yield_ch3cocho_ch3co_hco_t

  type, extends(quantum_yield_t) :: quantum_yield_ch3cocho_ch3co_hco_t
    ! Calculator for ch3cocho+hv->ch3co+hco quantum yield
  contains
    procedure :: calculate => run
  end type quantum_yield_ch3cocho_ch3co_hco_t

  interface quantum_yield_ch3cocho_ch3co_hco_t
    module procedure constructor
  end interface quantum_yield_ch3cocho_ch3co_hco_t

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

    allocate ( quantum_yield_ch3cocho_ch3co_hco_t :: this )

    call base_constructor( this, config, grid_warehouse, profile_warehouse )

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate the photorate quantum yield for a given set of environmental
  !! conditions
  function run( this, grid_warehouse, profile_warehouse )                     &
      result( quantum_yield )

    use musica_constants,              only : dk => musica_dk
    use tuvx_grid,                     only : grid_t
    use tuvx_grid_warehouse,           only : grid_warehouse_t
    use tuvx_profile,                  only : profile_t
    use tuvx_profile_warehouse,        only : profile_warehouse_t

    class(quantum_yield_ch3cocho_ch3co_hco_t), intent(in) :: this ! This :f:type:`~tuvx_quantum_yield_ch3cocho_ch3co_hco/quantum_yield_ch3cocho_ch3co_hco_t`
    type(grid_warehouse_t),    intent(inout) :: grid_warehouse ! A :f:type:`~tuvx_grid_warehouse/grid_warehouse_t`
    type(profile_warehouse_t), intent(inout) :: profile_warehouse ! A :f:type:`~tuvx_profile_warehouse/profile_warehouse_t`
    real(kind=dk), allocatable               :: quantum_yield(:,:) ! Calculated quantum_yield

    ! Local variables
    character(len=*), parameter :: Iam =                                      &
        'ch3cocho+hv->ch3co_hco quantum yield calculate'
    real(dk), parameter ::    rZERO = 0.0_dk
    real(dk), parameter ::    rONE  = 1.0_dk
    real(dk), parameter ::  lambdaL = 380._dk
    real(dk), parameter ::  lambdaU = 440._dk

    integer                       :: nzdim, vertNdx
    integer                       :: lambdaNdx
    real(dk)                      :: phi0, kq, lambda, airfac, qy
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
    quantum_yield = rZERO

    ! zero pressure yield:
    ! 1.0 for wc < 380 nm
    ! 0.0 for wc > 440 nm
    ! linear in between:

    ! Pressure correction: quenching coefficient, torr-1
    ! in air, Koch and Moortgat:
    do vertNdx = 1, nzdim
      airfac = modelDens( vertNdx ) * 760._dk / 2.456E19_dk
      do lambdaNdx = 1, lambdaGrid%ncells_
        lambda = lambdaGrid%mid_( lambdaNdx )
        phi0 = rONE - ( lambda - 380._dk ) / 60._dk
        phi0 = max( min( phi0, rONE ), rZERO )
        kq = 1.36e8_dk * exp( -8793._dk / lambda )
        if( phi0 > rZERO ) then
          if( lambda >= lambdaL .and. lambda <= lambdaU ) then
            qy = phi0 / ( phi0 + kq * airfac )
          else
            qy = phi0
          endif
        else
          qy = rZERO
        endif
        quantum_yield( lambdaNdx, vertNdx ) = qy
      enddo
    enddo

    quantum_yield = transpose( quantum_yield )

    deallocate( zGrid )
    deallocate( lambdaGrid )
    deallocate( mdlDensity )

  end function run

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_quantum_yield_ch3cocho_ch3co_hco
