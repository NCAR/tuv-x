! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0

module tuvx_quantum_yield_o3_o2_o1d
  ! The o3+hv->o2+o1d quantum yield type and related functions

  ! Including musica_config at the module level to avoid an ICE
  ! with Intel 2022.1 compiler
  use musica_config,                   only : config_t
  use tuvx_quantum_yield,              only : quantum_yield_t, base_constructor

  implicit none

  private
  public :: quantum_yield_o3_o2_o1d_t

  type, extends(quantum_yield_t) :: quantum_yield_o3_o2_o1d_t
    ! Calculator for o3+hv->o2+o1d quantum yield
  contains
    procedure :: calculate => run
  end type quantum_yield_o3_o2_o1d_t

  interface quantum_yield_o3_o2_o1d_t
    module procedure constructor
  end interface quantum_yield_o3_o2_o1d_t

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

    allocate ( quantum_yield_o3_o2_o1d_t :: this )

    call base_constructor( this, config, grid_warehouse, profile_warehouse )

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function run( this, grid_warehouse, profile_warehouse )                     &
      result( quantum_yield )
    ! Calculate the quantum yield for a given set of environmental conditions
    !
    ! Function to calculate the quantum yield O3 + hv -> O(1D) + O2,
    ! according to:
    ! Matsumi, Y., F. J. Comes, G. Hancock, A. Hofzumanhays, A. J. Hynes,
    ! M. Kawasaki, and A. R. Ravishankara, QUantum yields for production of O(1D)
    ! in the ultraviolet photolysis of ozone:  Recommendation based on evaluation
    ! of laboratory data, J. Geophys. Res., 107, `10.1029/2001JD000510. 
    ! <https://doi.org/10.1029/2001JD000510>`_, 2002.

    use musica_constants,              only : dk => musica_dk
    use tuvx_grid,                     only : grid_t
    use tuvx_grid_warehouse,           only : grid_warehouse_t
    use tuvx_profile,                  only : profile_t
    use tuvx_profile_warehouse,        only : profile_warehouse_t

    class(quantum_yield_o3_o2_o1d_t), intent(in) :: this ! This :f:type:`~tuvx_quantum_yield_o3_o2_o1d/quantum_yield_o3_o2_o1d_t`
    type(grid_warehouse_t),    intent(inout) :: grid_warehouse ! A :f:type:`~tuvx_grid_warehouse/grid_warehouse_t`
    type(profile_warehouse_t), intent(inout) :: profile_warehouse ! A :f:type:`~tuvx_profile_warehouse/profile_warehouse_t`
    real(kind=dk), allocatable               :: quantum_yield(:,:) ! Calculated quantum_yield

    ! Local variables
    real(dk), parameter :: a(3)  = (/ 0.8036_dk, 8.9061_dk, 0.1192_dk /)
    real(dk), parameter :: x(3)  = (/ 304.225_dk, 314.957_dk, 310.737_dk /)
    real(dk), parameter :: om(3) = (/ 5.576_dk, 6.601_dk, 2.187_dk /)
    real(dk), parameter ::   rZERO = 0.0_dk
    real(dk), parameter ::   rONE  = 1.0_dk

    character(len=*), parameter :: Iam = 'o3+hv->o2+o1d quantum yield calculate'

    integer     :: wNdx, vertNdx
    real(dk)    :: kt, q1, q2, T300, lambda
    real(dk)    :: qfac1, qfac2
    class(grid_t),    pointer :: lambdaGrid
    class(grid_t),    pointer :: zGrid
    class(profile_t), pointer :: temperature

    zGrid => grid_warehouse%get_grid( this%height_grid_ )
    lambdaGrid => grid_warehouse%get_grid( this%wavelength_grid_ )
    temperature => profile_warehouse%get_profile( this%temperature_profile_ )

    allocate( quantum_yield( lambdaGrid%ncells_, zGrid%ncells_ + 1 ) )
    quantum_yield = rZERO

    associate( w => lambdaGrid%mid_, Temp => temperature%edge_val_ )

    do vertNdx = 1, zGrid%ncells_ + 1
      kt = 0.695_dk * Temp( vertNdx )
      q1 = rONE
      q2 = exp( -825.518_dk / kt )
      qfac1 = q1 / (q1 + q2)
      qfac2 = q2 / (q1 + q2)
      T300 = Temp( vertNdx ) / 300._dk

      where( w(:) <= 305._dk )
        quantum_yield( :, vertNdx ) = 0.90_dk
      elsewhere( w(:) > 328._dk .and. w(:) <= 340._dk )
        quantum_yield( :, vertNdx ) = 0.08_dk
      endwhere
      do wNdx = 1, size( w )
        lambda = w( wNdx )
        if( lambda > 305._dk .and. lambda <= 328._dk ) then
          quantum_yield( wNdx, vertNdx ) = 0.0765_dk                          &
            + a(1) * qfac1 * EXP( -( (x(1) - lambda ) / om(1) )**4 )          &
            + a(2) * T300 * T300 * qfac2 *                                    &
                                     EXP( -( ( x(2) - lambda ) / om(2) )**2 ) &
            + a(3) * T300**1.5_dk * EXP( -( ( x(3) - lambda ) / om(3) )**2 )
        endif
      enddo
    enddo

    end associate

    quantum_yield = transpose( quantum_yield )

    deallocate( zGrid )
    deallocate( lambdaGrid )
    deallocate( temperature )

  end function run

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_quantum_yield_o3_o2_o1d
