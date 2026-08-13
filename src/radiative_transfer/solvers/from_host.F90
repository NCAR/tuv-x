! Copyright (C) 2026 University Corporation for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
module tuvx_solver_from_host
  ! Radiation field solver whose radiation field the host application supplies
  ! at runtime.
  !
  ! This solver does no radiative transfer calculation. It returns the
  ! radiation field that the host application sets through a
  ! :f:type:`~tuvx_solver_from_host/radiation_field_updater_t`.
  !
  ! The host application must supply the same quantities that the internal
  ! solvers produce. TUV-x normalizes the radiation field by the
  ! extraterrestrial flux. The components are therefore dimensionless
  ! transmission factors, and TUV-x forms the photolysis rate constant as:
  !
  ! .. math::
  !
  !    J = \sum_\lambda ( f_{dr} + f_{up} + f_{dn} )_\lambda \,
  !        F^\infty_\lambda \, \sigma_\lambda \, \phi_\lambda
  !
  ! where :math:`F^\infty` is the ``extraterrestrial flux`` profile. The host
  ! application must not fold the extraterrestrial flux or the Earth-Sun
  ! distance factor into the values that it supplies.
  !
  ! The solver starts with a radiation field of zero. A host application that
  ! never calls
  ! :f:func:`~tuvx_solver_from_host/radiation_field_updater_t%update` gets
  ! rates of zero. See :ref:`configuration-solvers-from-host` for more
  ! information.

  ! Including musica_config at the module level to avoid an ICE
  ! with the Intel compiler
#ifdef MUSICA_IS_INTEL_COMPILER
  use musica_config,                   only : config_t
#endif
  use musica_constants,                only : dk => musica_dk
  use tuvx_solver,                     only : solver_t, radiation_field_t

  implicit none

  private
  public :: solver_from_host_t, radiation_field_updater_t

  type, extends(solver_t) :: solver_from_host_t
    ! solver whose radiation field the host application sets at runtime
    private
    integer  :: n_vertical_interfaces_ = 0   ! number of vertical interfaces
    integer  :: n_wavelength_bins_     = 0   ! number of wavelength bins
    real(dk), allocatable :: edr_(:,:) ! direct component of the spectral irradiance (vertical interface, wavelength bin)
    real(dk), allocatable :: eup_(:,:) ! diffuse upwelling component of the spectral irradiance (vertical interface, wavelength bin)
    real(dk), allocatable :: edn_(:,:) ! diffuse downwelling component of the spectral irradiance (vertical interface, wavelength bin)
    real(dk), allocatable :: fdr_(:,:) ! direct component of the actinic flux (vertical interface, wavelength bin)
    real(dk), allocatable :: fup_(:,:) ! diffuse upwelling component of the actinic flux (vertical interface, wavelength bin)
    real(dk), allocatable :: fdn_(:,:) ! diffuse downwelling component of the actinic flux (vertical interface, wavelength bin)
  contains
    procedure :: update_radiation_field
    ! Returns the number of vertical interfaces the host must supply
    procedure :: number_of_vertical_interfaces
    ! Returns the number of wavelength bins the host must supply
    procedure :: number_of_wavelength_bins
    procedure :: pack_size
    procedure :: mpi_pack
    procedure :: mpi_unpack
  end type solver_from_host_t

  interface solver_from_host_t
    ! Constructor
    module procedure :: constructor
  end interface solver_from_host_t

  type :: radiation_field_updater_t
    ! Updater for a `solver_from_host_t` solver
    !
    ! A host application gets an updater from
    ! :f:func:`~tuvx_core/core_t%get_radiation_field_updater` and calls
    ! :f:func:`~tuvx_solver_from_host/radiation_field_updater_t%update` before
    ! each call to :f:func:`~tuvx_core/core_t%run`.
    private
    class(solver_from_host_t), pointer :: solver_ => null( )
  contains
    procedure :: update
  end type radiation_field_updater_t

  interface radiation_field_updater_t
    ! Constructor
    module procedure :: updater_constructor
  end interface radiation_field_updater_t

  real(dk), parameter :: rZERO = 0.0_dk

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function constructor( config, grid_warehouse, profile_warehouse )           &
      result( solver )
    ! Constructs a solver that the host application updates

    ! avoid a GCC13 ICE when including musica_config at the module level
#ifndef MUSICA_IS_INTEL_COMPILER
    use musica_config,                 only : config_t
#endif
    use musica_assert,                 only : assert_msg
    use musica_string,                 only : string_t
    use tuvx_grid,                     only : grid_t
    use tuvx_grid_warehouse,           only : grid_warehouse_t
    use tuvx_profile_warehouse,        only : profile_warehouse_t

    type(config_t),              intent(inout) :: config            ! solver configuration
    type(grid_warehouse_t),      intent(in)    :: grid_warehouse    ! available grids
    type(profile_warehouse_t),   intent(in)    :: profile_warehouse ! available profiles
    class(solver_from_host_t),   pointer       :: solver            ! new solver

    class(grid_t), pointer :: height_grid, wavelength_grid
    type(string_t) :: required_keys(1), optional_keys(0)

    required_keys(1) = "type"

    call assert_msg( 271038467,                                               &
                     config%validate( required_keys, optional_keys ),         &
                     "Bad configuration format for from host solver" )

    allocate( solver )

    ! The grid cell counts are fixed at construction, so the solver keeps the
    ! counts rather than a pointer to either grid
    height_grid     => grid_warehouse%get_grid( "height",     "km" )
    wavelength_grid => grid_warehouse%get_grid( "wavelength", "nm" )
    solver%n_vertical_interfaces_ = height_grid%ncells_ + 1
    solver%n_wavelength_bins_     = wavelength_grid%ncells_
    deallocate( height_grid     )
    deallocate( wavelength_grid )

    allocate( solver%edr_( solver%n_vertical_interfaces_,                     &
                           solver%n_wavelength_bins_ ) )
    allocate( solver%eup_( solver%n_vertical_interfaces_,                     &
                           solver%n_wavelength_bins_ ) )
    allocate( solver%edn_( solver%n_vertical_interfaces_,                     &
                           solver%n_wavelength_bins_ ) )
    allocate( solver%fdr_( solver%n_vertical_interfaces_,                     &
                           solver%n_wavelength_bins_ ) )
    allocate( solver%fup_( solver%n_vertical_interfaces_,                     &
                           solver%n_wavelength_bins_ ) )
    allocate( solver%fdn_( solver%n_vertical_interfaces_,                     &
                           solver%n_wavelength_bins_ ) )
    solver%edr_(:,:) = rZERO
    solver%eup_(:,:) = rZERO
    solver%edn_(:,:) = rZERO
    solver%fdr_(:,:) = rZERO
    solver%fup_(:,:) = rZERO
    solver%fdn_(:,:) = rZERO

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function updater_constructor( solver ) result( this )
    ! Constructs an updater for a `solver_from_host_t` solver

    class(solver_from_host_t), target, intent(inout) :: solver ! solver to update
    type(radiation_field_updater_t)                  :: this   ! new updater

    this%solver_ => solver

  end function updater_constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine update( this, direct_actinic_flux, upward_actinic_flux,          &
      downward_actinic_flux, direct_irradiance, upward_irradiance,            &
      downward_irradiance )
    ! Sets the radiation field from the host application
    !
    ! Every argument has the shape (vertical interface, wavelength bin). The
    ! vertical index follows the ``height`` grid, so index 1 is the lowest
    ! altitude and the last index is the top of the atmosphere. Both internal
    ! solvers return the radiation field in this order.
    !
    ! The three actinic flux components are required. TUV-x uses their sum for
    ! photolysis rate constants and for heating rates.
    !
    ! The three irradiance components are optional, but TUV-x uses their sum
    ! for dose rates. A host application that requests dose rates must supply
    ! them. TUV-x sets an omitted component to zero, so a dose rate that has
    ! no irradiance behind it is zero.
    !
    ! All values are dimensionless. Do not fold the extraterrestrial flux or
    ! the Earth-Sun distance factor into them.

    use musica_assert,                 only : assert_msg

    class(radiation_field_updater_t), intent(inout) :: this ! updater
    real(dk),           intent(in) :: direct_actinic_flux(:,:)    ! direct component of the normalized actinic flux (vertical interface, wavelength bin)
    real(dk),           intent(in) :: upward_actinic_flux(:,:)    ! diffuse upwelling component of the normalized actinic flux (vertical interface, wavelength bin)
    real(dk),           intent(in) :: downward_actinic_flux(:,:)  ! diffuse downwelling component of the normalized actinic flux (vertical interface, wavelength bin)
    real(dk), optional, intent(in) :: direct_irradiance(:,:)      ! direct component of the normalized spectral irradiance (vertical interface, wavelength bin)
    real(dk), optional, intent(in) :: upward_irradiance(:,:)      ! diffuse upwelling component of the normalized spectral irradiance (vertical interface, wavelength bin)
    real(dk), optional, intent(in) :: downward_irradiance(:,:)    ! diffuse downwelling component of the normalized spectral irradiance (vertical interface, wavelength bin)

    call assert_msg( 419530284, associated( this%solver_ ),                   &
                     "Cannot update an unassociated radiation field" )
    associate( solver => this%solver_ )

    call check_shape( 254866372, "direct actinic flux",                       &
                      direct_actinic_flux,   solver )
    call check_shape( 431793118, "upward actinic flux",                       &
                      upward_actinic_flux,   solver )
    call check_shape( 608719864, "downward actinic flux",                     &
                      downward_actinic_flux, solver )

    solver%fdr_(:,:) = direct_actinic_flux(:,:)
    solver%fup_(:,:) = upward_actinic_flux(:,:)
    solver%fdn_(:,:) = downward_actinic_flux(:,:)

    if( present( direct_irradiance ) ) then
      call check_shape( 785646610, "direct irradiance", direct_irradiance,    &
                        solver )
      solver%edr_(:,:) = direct_irradiance(:,:)
    else
      solver%edr_(:,:) = rZERO
    end if
    if( present( upward_irradiance ) ) then
      call check_shape( 962573356, "upward irradiance", upward_irradiance,    &
                        solver )
      solver%eup_(:,:) = upward_irradiance(:,:)
    else
      solver%eup_(:,:) = rZERO
    end if
    if( present( downward_irradiance ) ) then
      call check_shape( 139500102, "downward irradiance",                     &
                        downward_irradiance, solver )
      solver%edn_(:,:) = downward_irradiance(:,:)
    else
      solver%edn_(:,:) = rZERO
    end if

    end associate

  end subroutine update

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine check_shape( code, name, values, solver )
    ! Checks that a host-supplied array matches the shape of the radiation
    ! field

    use musica_assert,                 only : assert_msg
    use musica_string,                 only : to_char

    integer,                   intent(in) :: code      ! assertion code
    character(len=*),          intent(in) :: name      ! name of the array
    real(dk),                  intent(in) :: values(:,:) ! host-supplied array
    class(solver_from_host_t), intent(in) :: solver    ! solver to check against

    call assert_msg( code,                                                    &
        size( values, 1 ) == solver%n_vertical_interfaces_ .and.              &
        size( values, 2 ) == solver%n_wavelength_bins_,                       &
        "Bad shape for host-supplied "//name//". Expected ("//                &
        trim( to_char( solver%n_vertical_interfaces_ ) )//","//               &
        trim( to_char( solver%n_wavelength_bins_ ) )//") but got ("//         &
        trim( to_char( size( values, 1 ) ) )//","//                           &
        trim( to_char( size( values, 2 ) ) )//")" )

  end subroutine check_shape

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function update_radiation_field( this, solar_zenith_angle, n_layers,        &
      spherical_geometry, grid_warehouse, profile_warehouse,                  &
      radiator_warehouse ) result( radiation_field )
    ! Returns a copy of the radiation field that the host application set
    !
    ! This solver ignores the solar zenith angle, the spherical geometry, and
    ! the radiators. The host application accounts for them when it computes
    ! the radiation field.

    use musica_assert,                 only : assert_msg
    use musica_string,                 only : to_char
    use tuvx_grid_warehouse,           only : grid_warehouse_t
    use tuvx_profile_warehouse,        only : profile_warehouse_t
    use tuvx_radiator_warehouse,       only : radiator_warehouse_t
    use tuvx_spherical_geometry,       only : spherical_geometry_t

    class(solver_from_host_t),  intent(inout) :: this               ! from host solver
    integer,                    intent(in)    :: n_layers           ! number of vertical layers
    real(dk),                   intent(in)    :: solar_zenith_angle ! solar zenith angle [degrees]
    type(grid_warehouse_t),     intent(inout) :: grid_warehouse     ! available grids
    type(profile_warehouse_t),  intent(inout) :: profile_warehouse  ! available profiles
    type(radiator_warehouse_t), intent(inout) :: radiator_warehouse ! set of radiators
    type(spherical_geometry_t), intent(inout) :: spherical_geometry ! spherical geometry calculator

    type(radiation_field_t), pointer :: radiation_field

    call assert_msg( 493353776,                                               &
                     n_layers + 1 == this%n_vertical_interfaces_,             &
                     "The from host solver was built for "//                  &
                     trim( to_char( this%n_vertical_interfaces_ ) )//         &
                     " vertical interfaces but the height grid has "//        &
                     trim( to_char( n_layers + 1 ) )//"." )

    ! Return a copy. The caller takes ownership and deallocates the field
    ! before the next call to this function.
    radiation_field => radiation_field_t( this%n_vertical_interfaces_,         &
                                          this%n_wavelength_bins_ )
    radiation_field%edr_(:,:) = this%edr_(:,:)
    radiation_field%eup_(:,:) = this%eup_(:,:)
    radiation_field%edn_(:,:) = this%edn_(:,:)
    radiation_field%fdr_(:,:) = this%fdr_(:,:)
    radiation_field%fup_(:,:) = this%fup_(:,:)
    radiation_field%fdn_(:,:) = this%fdn_(:,:)

  end function update_radiation_field

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function number_of_vertical_interfaces( this ) result( n_interfaces )
    ! Returns the number of vertical interfaces that the host must supply

    class(solver_from_host_t), intent(in) :: this ! from host solver

    n_interfaces = this%n_vertical_interfaces_

  end function number_of_vertical_interfaces

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function number_of_wavelength_bins( this ) result( n_bins )
    ! Returns the number of wavelength bins that the host must supply

    class(solver_from_host_t), intent(in) :: this ! from host solver

    n_bins = this%n_wavelength_bins_

  end function number_of_wavelength_bins

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function pack_size( this, comm )
    ! Returns the size of a character buffer required to pack the solver

    use musica_mpi,                    only : musica_mpi_pack_size

    class(solver_from_host_t), intent(in) :: this ! solver to be packed
    integer,                   intent(in) :: comm ! MPI communicator

#ifdef MUSICA_USE_MPI
    pack_size = musica_mpi_pack_size( this%n_vertical_interfaces_, comm ) +   &
                musica_mpi_pack_size( this%n_wavelength_bins_, comm ) +       &
                musica_mpi_pack_size( this%edr_, comm ) +                     &
                musica_mpi_pack_size( this%eup_, comm ) +                     &
                musica_mpi_pack_size( this%edn_, comm ) +                     &
                musica_mpi_pack_size( this%fdr_, comm ) +                     &
                musica_mpi_pack_size( this%fup_, comm ) +                     &
                musica_mpi_pack_size( this%fdn_, comm )
#else
    pack_size = 0
#endif

  end function pack_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine mpi_pack( this, buffer, position, comm )
    ! Packs the solver onto a character buffer

    use musica_assert,                 only : assert
    use musica_mpi,                    only : musica_mpi_pack

    class(solver_from_host_t), intent(in)    :: this      ! solver to be packed
    character,                 intent(inout) :: buffer(:) ! memory buffer
    integer,                   intent(inout) :: position  ! current buffer position
    integer,                   intent(in)    :: comm      ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: prev_pos

    prev_pos = position
    call musica_mpi_pack( buffer, position, this%n_vertical_interfaces_, comm )
    call musica_mpi_pack( buffer, position, this%n_wavelength_bins_, comm )
    call musica_mpi_pack( buffer, position, this%edr_, comm )
    call musica_mpi_pack( buffer, position, this%eup_, comm )
    call musica_mpi_pack( buffer, position, this%edn_, comm )
    call musica_mpi_pack( buffer, position, this%fdr_, comm )
    call musica_mpi_pack( buffer, position, this%fup_, comm )
    call musica_mpi_pack( buffer, position, this%fdn_, comm )
    call assert( 806640397, position - prev_pos <= this%pack_size( comm ) )
#endif

  end subroutine mpi_pack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine mpi_unpack( this, buffer, position, comm )
    ! Unpacks a solver from a character buffer

    use musica_assert,                 only : assert
    use musica_mpi,                    only : musica_mpi_unpack

    class(solver_from_host_t), intent(out)   :: this      ! solver to be unpacked
    character,                 intent(inout) :: buffer(:) ! memory buffer
    integer,                   intent(inout) :: position  ! current buffer position
    integer,                   intent(in)    :: comm      ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: prev_pos

    prev_pos = position
    call musica_mpi_unpack( buffer, position, this%n_vertical_interfaces_,     &
                            comm )
    call musica_mpi_unpack( buffer, position, this%n_wavelength_bins_, comm )
    call musica_mpi_unpack( buffer, position, this%edr_, comm )
    call musica_mpi_unpack( buffer, position, this%eup_, comm )
    call musica_mpi_unpack( buffer, position, this%edn_, comm )
    call musica_mpi_unpack( buffer, position, this%fdr_, comm )
    call musica_mpi_unpack( buffer, position, this%fup_, comm )
    call musica_mpi_unpack( buffer, position, this%fdn_, comm )
    call assert( 349040634, position - prev_pos <= this%pack_size( comm ) )
#endif

  end subroutine mpi_unpack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_solver_from_host
