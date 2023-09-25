! Copyright (C) 2022 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
module tuvx_solver_discrete_ordinate

  use musica_constants,               only : dk => musica_dk
  use tuvx_constants,                 only : kPi => pi
  use tuvx_discrete_ordinate_util,    only : solver_constants_t
  use tuvx_grid_warehouse,            only : grid_warehouse_ptr
  use tuvx_profile_warehouse,         only : profile_warehouse_ptr
  use tuvx_solver,                    only : solver_t, radiation_field_t

  implicit none

  private
  public :: solver_discrete_ordinate_t

  type, extends(solver_t) :: solver_discrete_ordinate_t
     ! Radiative flux calculator that applies the discrete-ordinate method.
     integer :: n_streams_ ! number of streams
     type(grid_warehouse_ptr) :: height_grid_
     type(grid_warehouse_ptr) :: wavelength_grid_
     type(profile_warehouse_ptr) :: surface_albedo_profile_
     type(solver_constants_t) :: solver_constants_
   contains
    procedure :: update_radiation_field
    procedure :: pack_size
    procedure :: mpi_pack
    procedure :: mpi_unpack
  end type solver_discrete_ordinate_t

  interface solver_discrete_ordinate_t
    procedure :: constructor
  end interface solver_discrete_ordinate_t

  real(dk), parameter :: rZERO = 0.0_dk
  real(dk), parameter :: rONE  = 1.0_dk
  real(dk), parameter :: rTWO  = 2.0_dk
  real(dk), parameter :: kFOURPi  = 4.0_dk*kPi
  real(dk), parameter :: d2r   = kPi/180._dk

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Constructs a discrete ordinate solver
  function constructor( config, grid_warehouse, profile_warehouse )           &
      result( solver )

    use musica_assert,                 only : assert_msg
    use musica_config,                 only : config_t
    use musica_string,                 only : string_t
    use tuvx_grid_warehouse,           only : grid_warehouse_t
    use tuvx_profile_warehouse,        only : profile_warehouse_t

    type(config_t),                    intent(inout) :: config
    class(solver_discrete_ordinate_t), pointer       :: solver
    type(grid_warehouse_t),            intent(in)    :: grid_warehouse
    type(profile_warehouse_t),         intent(in)    :: profile_warehouse

    character(len=*), parameter :: Iam = "Discrete ordinate solver constrctor"
    type(string_t) :: required_keys(2), optional_keys(0)

    required_keys(1) = "type"
    required_keys(2) = "number of streams"

    call assert_msg( 108215931,                                               &
                     config%validate( required_keys, optional_keys ),         &
                     "Bad configuration format for discrete ordinate solver" )

    allocate( solver )

    call config%get( "number of streams", solver%n_streams_, Iam )

    call assert_msg( 326135075, solver%n_streams_ >= 2 .and.                  &
                                solver%n_streams_ <= 32,                      &
                     "Discrete ordinate solver must have between 2 and 32 "// &
                     "streams" )
    call assert_msg( 741820725, mod( solver%n_streams_, 2 ) == 0,             &
                     "Discrete ordinate solver must have an even number of "//&
                     "streams" )
    solver%height_grid_ = grid_warehouse%get_ptr( "height", "km" )
    solver%wavelength_grid_ = grid_warehouse%get_ptr( "wavelength", "nm" )
    solver%surface_albedo_profile_ =                                          &
        profile_warehouse%get_ptr( "surface albedo", "none" )
    solver%solver_constants_ = solver_constants_t( )

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function update_radiation_field( this, solar_zenith_angle, n_layers,        &
      spherical_geometry, grid_warehouse, profile_warehouse,                  &
      radiator_warehouse ) result( radiation_field )

    use musica_string,                 only : string_t
    use tuvx_diagnostic_util,          only : diagout
    use tuvx_grid,                     only : grid_t
    use tuvx_grid_warehouse,           only : grid_warehouse_t
    use tuvx_linear_algebra_linpack,   only : linear_algebra_linpack_t
    use tuvx_profile,                  only : profile_t
    use tuvx_profile_warehouse,        only : Profile_warehouse_t
    use tuvx_radiator,                 only : radiator_t, radiator_state_t
    use tuvx_radiator_warehouse,       only : radiator_warehouse_t
    use tuvx_radiator_warehouse,       only : warehouse_iterator_t
    use tuvx_solver,                   only : slant_optical_depth
    use tuvx_spherical_geometry,       only : spherical_geometry_t
    use tuvx_discrete_ordinate_util,   only : psndo

    class(solver_discrete_ordinate_t), intent(inout) :: this ! Discrete ordinate solver

    integer,                    intent(in)    :: n_layers  ! number of vertical layers
    real(dk),                   intent(in)    :: solar_zenith_angle ! solar zenith angle [degrees]
    type(grid_warehouse_t),     intent(inout) :: grid_warehouse
    type(Profile_warehouse_t),  intent(inout) :: profile_warehouse
    type(radiator_warehouse_t), intent(inout) :: radiator_warehouse
    type(spherical_geometry_t), intent(inout) :: spherical_geometry

    type(radiation_field_t),   pointer       :: radiation_field

    ! Local variables
    character(len=*), parameter :: Iam = 'Update radiation field: '

    ! Local variables
    real(dk), parameter                  :: largest = 1.e36_dk
    real(dk), parameter                  :: kfloor = rONE/largest
    real(dk), parameter                  :: precis = 1.e-7_dk
    real(dk), parameter                  :: eps    = 1.e-3_dk

    integer                              :: nlambda, lambdaNdx
    integer                              :: streamNdx
    real(dk)                             :: umu0
    real(dk), allocatable                :: pmom(:,:)
    type(radiator_state_t)               :: atmRadiatorState
    class(grid_t),    pointer            :: zGrid
    class(grid_t),    pointer            :: lambdaGrid
    class(profile_t), pointer            :: surfaceAlbedo

    zGrid => grid_warehouse%get_grid( this%height_grid_ )
    lambdaGrid => grid_warehouse%get_grid( this%wavelength_grid_ )
    surfaceAlbedo =>                                                          &
        profile_warehouse%get_profile( this%surface_albedo_profile_ )

    nlambda = lambdaGrid%ncells_
    radiation_field => radiation_field_t( n_layers + 1, nlambda )

    allocate( atmRadiatorState%layer_G_( n_layers, nlambda, this%n_streams_ ) )
    ! Create cumulative state from all radiators
    call radiator_warehouse%accumulate_states( atmRadiatorState )

    ! UMU0   = cosine of solar zenith angle
    ! ALBEDO = surface albedo
    ! DTAUC  =  unscaled optical depth of each layer
    ! SSALB  =  unscaled single scattering albedo
    ! GU     =  unscaled asymmetry factor
    ! N_LAYERS = number of layers in the atmosphere

    umu0 = cos( solar_zenith_angle*d2r )
    associate( nid  => spherical_geometry%nid_,                               &
               dsdh => spherical_geometry%dsdh_ )

    allocate( pmom( 0:this%n_streams_, n_layers ) )

    radiation_field%fdr_ = rZERO
    radiation_field%fup_ = rZERO
    radiation_field%fdn_ = rZERO
    wavelength_loop: do lambdaNdx = 1, nlambda
      associate( albedo => surfaceAlbedo%mid_val_( lambdaNdx ),                 &
             dtauc => atmRadiatorState%layer_OD_( n_layers:1:-1, lambdaNdx ),   &
             ssalb  => atmRadiatorState%layer_SSA_( n_layers:1:-1, lambdaNdx ) )
        pmom( 0, : ) = rONE
        do streamNdx = 1, this%n_streams_
          pmom( streamNdx, : ) = atmRadiatorState%layer_G_( n_layers:1:-1, lambdaNdx, streamNdx )
        end do
        associate( fdr => radiation_field%fdr_(: , lambdaNdx ),        &
                   fup => radiation_field%fup_(: , lambdaNdx ),        &
                   fdn => radiation_field%fdn_(: , lambdaNdx ),        &
                   edr => radiation_field%edr_(: , lambdaNdx ),        &
                   eup => radiation_field%eup_(: , lambdaNdx ),        &
                   edn => radiation_field%edn_(: , lambdaNdx ) )
          call psndo( dsdh, nid, n_layers,                                    &
                      dtauc, ssalb, pmom,                                     &
                      albedo, this%n_streams_, umu0,                          &
                      edr, edn, eup,                                          &
                      fdr, fup, fdn, this%solver_constants_ )
          fdr = fdr(n_layers+1:1:-1) * kFOURPi
          fup = fup(n_layers+1:1:-1) * kFOURPi
          fdn = fdn(n_layers+1:1:-1) * kFOURPi
          edr = edr(n_layers+1:1:-1)
          eup = eup(n_layers+1:1:-1)
          edn = edn(n_layers+1:1:-1)
        end associate
      end associate
    enddo wavelength_loop

    end associate

    deallocate( zGrid )
    deallocate( lambdaGrid )
    deallocate( surfaceAlbedo )
    deallocate( pmom )

  end function update_radiation_field

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function pack_size( this, comm )
    ! Returns the size of a character buffer required to pack the solver

    use musica_mpi,                    only : musica_mpi_pack_size

    class(solver_discrete_ordinate_t), intent(in) :: this ! solver to be packed
    integer,                           intent(in) :: comm ! MPI communicator

#ifdef MUSICA_USE_MPI
    pack_size = musica_mpi_pack_size( this%n_streams_, comm ) +               &
                this%height_grid_%pack_size(            comm ) +              &
                this%wavelength_grid_%pack_size(        comm ) +              &
                this%surface_albedo_profile_%pack_size( comm )
#else
    pack_size = 0
#endif

  end function pack_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine mpi_pack( this, buffer, position, comm )
    ! Packs the solver onto a character buffer

    use musica_assert,                 only : assert
    use musica_mpi,                    only : musica_mpi_pack

    class(solver_discrete_ordinate_t), intent(in)    :: this      ! solver to be packed
    character,                         intent(inout) :: buffer(:) ! memory buffer
    integer,                           intent(inout) :: position  ! current buffer position
    integer,                           intent(in)    :: comm      ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: prev_pos
    prev_pos = position

    call musica_mpi_pack( buffer, position, this%n_streams_, comm )
    call this%height_grid_%mpi_pack(            buffer, position, comm )
    call this%wavelength_grid_%mpi_pack(        buffer, position, comm )
    call this%surface_albedo_profile_%mpi_pack( buffer, position, comm )
    call assert( 807855666, position - prev_pos <= this%pack_size( comm ) )
#endif

  end subroutine mpi_pack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine mpi_unpack( this, buffer, position, comm )
    ! Unpacks a solver from a character buffer

    use musica_assert,                 only : assert
    use musica_mpi,                    only : musica_mpi_unpack

    class(solver_discrete_ordinate_t), intent(out)   :: this      ! solver to be packed
    character,                         intent(inout) :: buffer(:) ! memory buffer
    integer,                           intent(inout) :: position  ! current buffer position
    integer,                           intent(in)    :: comm      ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: prev_pos
    prev_pos = position

    call musica_mpi_unpack( buffer, position, this%n_streams_, comm )
    call this%height_grid_%mpi_unpack(            buffer, position, comm )
    call this%wavelength_grid_%mpi_unpack(        buffer, position, comm )
    call this%surface_albedo_profile_%mpi_unpack( buffer, position, comm )
    call assert( 120962799, position - prev_pos <= this%pack_size( comm ) )
    this%solver_constants_ = solver_constants_t( )
#endif

  end subroutine mpi_unpack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_solver_discrete_ordinate
