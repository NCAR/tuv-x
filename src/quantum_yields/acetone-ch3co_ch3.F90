! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0

module tuvx_quantum_yield_ch3coch3_ch3co_ch3
  ! This calculates the quantum_yield for acetone

  ! Including musica_config at the module level to avoid an ICE
  ! with Intel 2022.1 compiler
  use musica_constants,                only : dk => musica_dk
  use musica_config,                   only : config_t
  use tuvx_quantum_yield,              only : quantum_yield_t, base_constructor

  implicit none

  private
  public :: quantum_yield_ch3coch3_ch3co_ch3_t

  type, extends(quantum_yield_t) :: quantum_yield_ch3coch3_ch3co_ch3_t
    ! Calculator for acetone quantum_yield
    logical :: do_CO_ = .false.
    logical :: do_CH3CO_ = .false.
    real(kind=dk) :: low_wavelength_value_
    real(kind=dk) :: high_wavelength_value_
    real(kind=dk) :: minimum_temperature_
    real(kind=dk) :: maximum_temperature_
  contains
    !> Initialize the quantum_yield
    procedure :: calculate => run
    ! returns the number of bytes required to pack the object onto a buffer
    procedure :: pack_size
    ! packs the object onto a character buffer
    procedure :: mpi_pack
    ! unpacks an object from a character buffer
    procedure :: mpi_unpack
  end type quantum_yield_ch3coch3_ch3co_ch3_t

  interface quantum_yield_ch3coch3_ch3co_ch3_t
    module procedure constructor
  end interface quantum_yield_ch3coch3_ch3co_ch3_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function constructor( config, grid_warehouse, profile_warehouse )           &
      result( this )
    ! Build the quantum yield

    use musica_assert,                 only : assert_msg, die_msg
    use musica_string,                 only : string_t
    use tuvx_grid_warehouse,           only : grid_warehouse_t
    use tuvx_profile_warehouse,        only : profile_warehouse_t

    type(quantum_yield_ch3coch3_ch3co_ch3_t), pointer :: this ! This :f:type:`~tuvx_quantum_yield/quantum_yield_t` calculator
    type(config_t),                           intent(inout) :: config ! Quantum yield configuration data
    type(grid_warehouse_t),                   intent(inout) :: grid_warehouse ! A :f:type:`~tuvx_grid_warehouse/grid_warehouse_t`
    type(profile_warehouse_t),                intent(inout) :: profile_warehouse ! A :f:type:`~tuvx_profile_warehouse/profile_warehouse_t`

    character(len=*), parameter :: my_name =                                  &
        "Acetone quantum yield constructor"
    type(string_t) :: required_keys(1), optional_keys(6), branch

    required_keys(1) = "type"
    optional_keys(1) = "name"
    optional_keys(2) = "low wavelength value"
    optional_keys(3) = "high wavelength value"
    optional_keys(4) = "minimum temperature"
    optional_keys(5) = "maximum temperature"
    optional_keys(6) = "branch"
    call assert_msg( 253342443,                                               &
                     config%validate( required_keys, optional_keys ),         &
                     "Bad configuration for acetone quantum yield." )
    allocate ( this )
    call base_constructor( this, config, grid_warehouse, profile_warehouse )
    call config%get( "low wavelength value", this%low_wavelength_value_,      &
                     my_name, default = 0.95_dk )
    call config%get( "high wavelength value", this%high_wavelength_value_,    &
                     my_name, default = 0.0_dk )
    call config%get( "minimum temperature", this%minimum_temperature_,        &
                     my_name, default = 0.0_dk )
    call config%get( "maximum temperature", this%maximum_temperature_,        &
                     my_name, default = huge( 1.0_dk ) )
    call config%get( "branch", branch, my_name, default = "CH3CO" )
    if( branch .eq. "CO" ) then
      this%do_CO_ = .true.
      this%do_CH3CO_ = .false.
    else if( branch .eq. "CH3CO" ) then
      this%do_CO_ = .false.
      this%do_CH3CO_ = .true.
    else if( branch .eq. "CO+CH3CO" ) then
      this%do_CO_ = .true.
      this%do_CH3CO_ = .true.
    else
      call die_msg( 534162111, "Invalid branch for acetone quantum yield: '"//&
                               branch//"'." )
    end if

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
    real(dk)    :: fco, fac, qy

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
      Tadj = max( this%minimum_temperature_,                                  &
                  min( this%maximum_temperature_, modelTemp( vertNdx ) ) )    &
                / 295._dk
      M    = modelDens( vertNdx )
lambda_loop: &
      do lambdaNdx = 1, lambdaGrid%ncells_
        w = lambdaGrid%mid_( lambdaNdx )
        if( w < 279._dk ) then
           qy = this%low_wavelength_value_
        elseif( w > 327._dk ) then
           qy = this%high_wavelength_value_
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
          qy = 0.0_dk
          if( this%do_CO_ ) qy = qy + fco
          if( this%do_CH3CO_ ) qy = qy + fac
        endif
        quantum_yield( lambdaNdx, vertNdx ) = qy
      enddo lambda_loop
    enddo vert_loop

    quantum_yield = transpose( quantum_yield )

    deallocate( zGrid )
    deallocate( lambdaGrid )
    deallocate( mdlTemperature )
    deallocate( mdlDensity )

  end function run

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the number of bytes required to pack the object onto a buffer
  integer function pack_size( this, comm )

    use musica_mpi,                    only : musica_mpi_pack_size

    !> Quantum yield to be packed
    class(quantum_yield_ch3coch3_ch3co_ch3_t), intent(in) :: this
    !> MPI communicator
    integer,                                   intent(in) :: comm

#ifdef MUSICA_USE_MPI
    pack_size = this%quantum_yield_t%pack_size( comm ) +                      &
                musica_mpi_pack_size( this%do_CO_,                 comm ) +   &
                musica_mpi_pack_size( this%do_CH3CO_,              comm ) +   &
                musica_mpi_pack_size( this%low_wavelength_value_,  comm ) +   &
                musica_mpi_pack_size( this%high_wavelength_value_, comm ) +   &
                musica_mpi_pack_size( this%minimum_temperature_,   comm ) +   &
                musica_mpi_pack_size( this%maximum_temperature_,   comm )
#else
    pack_size = 0
#endif

  end function pack_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Packs the quantum yield onto a character buffer
  subroutine mpi_pack( this, buffer, position, comm )

    use musica_assert,                 only : assert
    use musica_mpi,                    only : musica_mpi_pack

    !> Quantum yield to pack
    class(quantum_yield_ch3coch3_ch3co_ch3_t), intent(in)    :: this
    !> Memory buffer
    character,                                 intent(inout) :: buffer(:)
    !> Current buffer position
    integer,                                   intent(inout) :: position
    !> MPI communicator
    integer,                                   intent(in)    :: comm

#ifdef MUSICA_USE_MPI
    integer :: prev_pos

    prev_pos = position
    call this%quantum_yield_t%mpi_pack( buffer, position, comm )
    call musica_mpi_pack( buffer, position, this%do_CO_,                 comm )
    call musica_mpi_pack( buffer, position, this%do_CH3CO_,              comm )
    call musica_mpi_pack( buffer, position, this%low_wavelength_value_,  comm )
    call musica_mpi_pack( buffer, position, this%high_wavelength_value_, comm )
    call musica_mpi_pack( buffer, position, this%minimum_temperature_,   comm )
    call musica_mpi_pack( buffer, position, this%maximum_temperature_,   comm )
    call assert( 985830490, position - prev_pos <= this%pack_size( comm ) )
#endif

  end subroutine mpi_pack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Unpacks a quantum yield calculator from a character buffer
  subroutine mpi_unpack( this, buffer, position, comm )

    use musica_assert,                 only : assert
    use musica_mpi,                    only : musica_mpi_unpack

    !> Quantum yield to unpack
    class(quantum_yield_ch3coch3_ch3co_ch3_t), intent(out)   :: this
    !> Memory buffer
    character,                                 intent(inout) :: buffer(:)
    !> Current buffer position
    integer,                                   intent(inout) :: position
    !> MPI communicator
    integer,                                   intent(in)    :: comm

#ifdef MUSICA_USE_MPI
    integer :: prev_pos

    prev_pos = position
    call this%quantum_yield_t%mpi_unpack( buffer, position, comm )
    call musica_mpi_unpack( buffer, position, this%do_CO_,                 comm )
    call musica_mpi_unpack( buffer, position, this%do_CH3CO_,              comm )
    call musica_mpi_unpack( buffer, position, this%low_wavelength_value_,  comm )
    call musica_mpi_unpack( buffer, position, this%high_wavelength_value_, comm )
    call musica_mpi_unpack( buffer, position, this%minimum_temperature_,   comm )
    call musica_mpi_unpack( buffer, position, this%maximum_temperature_,   comm )
    call assert( 301844101, position - prev_pos <= this%pack_size( comm ) )
#endif

  end subroutine mpi_unpack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_quantum_yield_ch3coch3_ch3co_ch3
