! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0

module tuvx_quantum_yield_ch3coch3_ch3co_ch3
  ! This calculates the quantum_yield for acetone

  ! Including musica_config at the module level to avoid an ICE
  ! with Intel 2022.1 compiler
  use musica_config,                   only : config_t
  use tuvx_quantum_yield,              only : quantum_yield_t, base_constructor

  implicit none

  private
  public :: quantum_yield_ch3coch3_ch3co_ch3_t

  type, extends(quantum_yield_t) :: quantum_yield_ch3coch3_ch3co_ch3_t
    ! Calculator for acetone quantum_yield
  contains
    !> Initialize the quantum_yield
    procedure :: calculate => run
  end type quantum_yield_ch3coch3_ch3co_ch3_t

  interface quantum_yield_ch3coch3_ch3co_ch3_t
    module procedure constructor
  end interface quantum_yield_ch3coch3_ch3co_ch3_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function constructor( config, grid_warehouse, profile_warehouse )           &
      result( this )
    ! Build the quantum yield

    use tuvx_grid_warehouse,           only : grid_warehouse_t
    use tuvx_profile_warehouse,        only : profile_warehouse_t

    class(quantum_yield_t), pointer :: this ! This :f:type:`~tuvx_quantum_yield/quantum_yield_t` calculator
    type(config_t),            intent(inout) :: config ! Quantum yield configuration data
    type(grid_warehouse_t),    intent(inout) :: grid_warehouse ! A :f:type:`~tuvx_grid_warehouse/grid_warehouse_t`
    type(profile_warehouse_t), intent(inout) :: profile_warehouse ! A :f:type:`~tuvx_profile_warehouse/profile_warehouse_t`

    allocate ( quantum_yield_ch3coch3_ch3co_ch3_t :: this )

    call base_constructor( this, config, grid_warehouse, profile_warehouse )

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function run( this, grid_warehouse, profile_warehouse )                     &
      result( quantum_yield )
    ! Calculate the photorate quantum_yield for a given set of environmental
    ! conditions
    ! qyacet - q.y. for acetone, based on Blitz et al. (2004)
    ! Compute acetone quantum yields according to the parameterization of:
    ! Blitz, M. A., D. E. Heard, M. J. Pilling, S. R. Arnold, and
    ! M. P. Chipperfield
    ! (2004), Pressure and temperature-dependent quantum yields for the
    ! photodissociation of acetone between 279 and 327.5 nm, Geophys.
    ! Res. Lett., 31, L06111, `doi:10.1029/2003GL018793. 
    ! <https://doi.org/10.1029/2003GL018793>`_

    use musica_constants,              only : dk => musica_dk
    use tuvx_grid,                     only : grid_t
    use tuvx_grid_warehouse,           only : grid_warehouse_t
    use tuvx_profile,                  only : profile_t
    use tuvx_profile_warehouse,        only : profile_warehouse_t

    class(quantum_yield_ch3coch3_ch3co_ch3_t), intent(in) :: this ! This :f:type:`~tuvx_quantum_yield_ch3coch3_ch3co_ch3/quantum_yield_ch3coch3_ch3co_ch3_t`
    type(grid_warehouse_t),    intent(inout) :: grid_warehouse ! A :f:type:`~tuvx_grid_warehouse/grid_warehouse_t`
    type(profile_warehouse_t), intent(inout) :: profile_warehouse ! A :f:type:`~tuvx_profile_warehouse/profile_warehouse_t`
    real(kind=dk), allocatable               :: quantum_yield(:,:) ! Calculated quantum_yield

    character(len=*), parameter :: Iam =                                      &
      'ch3coch3+hv->ch3co+ch3 quantum_yield calculate'
    real(dk), parameter :: rZERO = 0.0_dk
    real(dk), parameter :: rONE  = 1.0_dk

    integer                       :: lambdaNdx
    integer                       :: nzdim, vertNdx
    real(dk),         allocatable :: modelTemp(:), modelDens(:)
    class(grid_t),    pointer     :: zGrid
    class(grid_t),    pointer     :: lambdaGrid
    class(profile_t), pointer     :: mdlTemperature
    class(profile_t), pointer     :: mdlDensity

    ! w = wavelength, nm
    ! T = temperature, K
    ! M = air number density, molec. cm-3
    real(dk)    :: w, wadj, Tadj, M
    real(dk)    :: a0, a1, a2, a3, a4
    real(dk)    :: b0, b1, b2, b3, b4
    real(dk)    :: c3
    real(dk)    :: cA0, cA1, cA2, cA3, cA4
    real(dk)    :: dumexp
    real(dk)    :: fco, fac

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

vert_loop: &
    do vertNdx = 1, nzdim
      Tadj = modelTemp( vertNdx ) / 295._dk
      M    = modelDens( vertNdx )
lambda_loop: &
      do lambdaNdx = 1, lambdaGrid%ncells_
        w = lambdaGrid%mid_( lambdaNdx )
        if( w < 279._dk ) then
           fac = 0.95_dk
        elseif( w > 327._dk ) then
           fac = rZERO
        else
          ! CO (carbon monoxide) quantum yields:
          a0 = 0.350_dk * Tadj**( -1.28_dk )
          b0 = 0.068_dk * Tadj**( -2.65_dk )
          ! SM: prevent exponent overflow in rare cases:
          dumexp = b0 * ( w - 248._dk )
          if( dumexp > 80._dk ) then
            cA0 = 5.e34_dk
          else
            cA0 = exp( dumexp ) * a0 / ( rONE - a0 )
          endif

          fco = rONE / ( rONE + cA0 )
          ! CH3CO (acetyl radical) quantum yields:
          wadj = 1.e7_dk / w
          if( w >= 279._dk .and. w < 302._dk ) then
            a1 = 1.600E-19_dk * Tadj**( -2.38_dk )
            b1 = 0.55E-3_dk   * Tadj**( -3.19_dk )
            cA1 = a1 * EXP( -b1 * ( wadj - 33113._dk ) )
            fac = ( rONE - fco ) / ( rONE + cA1 * M )
          else if( w >= 302._dk .and. w <= 327._dk ) then
            a2 = 1.62E-17_dk * Tadj**( -10.03_dk )
            b2 = 1.79E-3_dk  * Tadj**( -1.364_dk )
            cA2 = a2 * EXP( -b2 * ( wadj - 30488._dk ) )

            a3 = 26.29_dk   * Tadj**( -6.59_dk )
            b3 = 5.72E-7_dk * Tadj**( -2.93_dk )
            c3 = 30006._dk  * Tadj**( -0.064_dk )
            ca3 = a3 * EXP( -b3 * ( wadj - c3 )**2 )

            a4 = 1.67E-15_dk * Tadj**( -7.25_dk )
            b4 = 2.08E-3_dk  * Tadj**( -1.16_dk )
            cA4 = a4 * EXP( -b4 * ( wadj - 30488._dk ) )

            fac = ( rONE - fco ) * ( rONE + cA3 + cA4 * M ) &
                  / ( ( rONE + cA3 + cA2 * M ) * ( rONE + cA4 * M ) )
          endif
        endif
        quantum_yield( lambdaNdx, vertNdx ) = fac
      enddo lambda_loop
    enddo vert_loop

    quantum_yield = transpose( quantum_yield )

    deallocate( zGrid )
    deallocate( lambdaGrid )
    deallocate( mdlTemperature )
    deallocate( mdlDensity )

  end function run

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_quantum_yield_ch3coch3_ch3co_ch3
