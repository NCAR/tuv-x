! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
module tuvx_solver
! General interface for radiative transfer solvers

  use musica_constants,               only : dk => musica_dk

  implicit none

  private
  public :: solver_t, radiation_field_t, slant_optical_depth

  type :: radiation_field_t
    real(dk), allocatable :: edr_(:,:) ! Contribution of the direct component to the total spectral irradiance (vertical interface, wavelength)
    real(dk), allocatable :: eup_(:,:) ! Contribution of the diffuse upwelling component to the total spectral irradiance (vertical interface, wavelength)
    real(dk), allocatable :: edn_(:,:) ! Contribution of the diffuse downwelling component to the total spectral irradiance (vertical interface, wavelength)
    real(dk), allocatable :: fdr_(:,:) ! Contribution of the direct component to the total actinic flux (vertical interface, wavelength)
    real(dk), allocatable :: fup_(:,:) ! Contribution of the diffuse upwelling component to the total actinic flux (vertical interface, wavelength)
    real(dk), allocatable :: fdn_(:,:) ! Contribution of the diffuse downwelling component to the total actinic flux (vertical interface, wavelength)
  contains
    ! Scale the radiation field values
    procedure :: apply_scale_factor
    ! Returns the number of bytes needed to pack the object onto a buffer
    procedure :: pack_size => field_pack_size
    ! Packs the object onto a character buffer
    procedure :: mpi_pack => field_pack
    ! Unpacks data from a character buffer into the object
    procedure :: mpi_unpack => field_unpack
    final :: finalize
  end type radiation_field_t

  type, abstract :: solver_t
    contains
    procedure(update_radiation_field), deferred :: update_radiation_field
    ! Returns the number of bytes needed to pack the object onto a buffer
    procedure :: pack_size => solver_pack_size
    ! Packs the object onto a character buffer
    procedure :: mpi_pack => solver_pack
    ! Unpacks data from a character buffer into the object
    procedure :: mpi_unpack => solver_unpack
  end type solver_t

  interface radiation_field_t
    module procedure :: field_constructor
  end interface radiation_field_t

interface

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function update_radiation_field( this, solar_zenith_angle, n_layers,        &
      spherical_geometry, grid_warehouse, profile_warehouse,                  &
      radiator_warehouse ) result( radiation_field )
    ! Solves for the radiation field based on given conditions

     use musica_constants,             only : dk => musica_dk
     use tuvx_grid_warehouse,          only : grid_warehouse_t
     use tuvx_profile_warehouse,       only : profile_warehouse_t
     use tuvx_radiator_warehouse,      only : radiator_warehouse_t
     use tuvx_spherical_geometry,      only : spherical_geometry_t

     import solver_t, radiation_field_t

     class(solver_t), intent(inout) :: this
     integer, intent(in)                       :: n_layers           ! Number of vertical layers
     real(dk), intent(in)                      :: solar_zenith_angle ! Solar zenith angle [degrees]
     type(grid_warehouse_t), intent(inout)     :: grid_warehouse     ! Available grids
     type(profile_warehouse_t), intent(inout)  :: profile_warehouse  ! Available profiles
     type(radiator_warehouse_t), intent(inout) :: radiator_warehouse ! Set of radiators
     type(spherical_geometry_t), intent(inout) :: spherical_geometry ! Spherical geometry calculator

     type(radiation_field_t), pointer         :: radiation_field

  end function update_radiation_field

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end interface

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function solver_pack_size( this, comm )
    ! Returns the size of a character buffer required to pack the solver

    class(solver_t),   intent(in) :: this ! solver to be packed
    integer,           intent(in) :: comm ! MPI communicator

    solver_pack_size = 0

#ifdef MUSICA_USE_MPI
    ! nothing to do for now
#endif
  end function solver_pack_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine solver_pack( this, buffer, position, comm )
    ! Packs the solver onto a character buffer

    class(solver_t), intent(in)    :: this              ! solver to be packed
    character,               intent(inout) :: buffer(:) ! memory buffer
    integer,                 intent(inout) :: position  ! current buffer position
    integer,                 intent(in)    :: comm      ! MPI communicator

#ifdef MUSICA_USE_MPI
    ! nothing to do for now
#endif
  end subroutine solver_pack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine solver_unpack( this, buffer, position, comm )
    ! Unpacks a solver from a character buffer

    class(solver_t),         intent(out)   :: this      ! solver to be packed
    character,               intent(inout) :: buffer(:) ! memory buffer
    integer,                 intent(inout) :: position  ! current buffer position
    integer,                 intent(in)    :: comm      ! MPI communicator

#ifdef MUSICA_USE_MPI
    ! nothing to do for now
#endif
  end subroutine solver_unpack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function field_constructor( n_vertical_interfaces, n_wavelength_bins )      &
      result( field )
    ! Constructor of radiation field objects

    type(radiation_field_t), pointer :: field
    integer, intent(in) :: n_vertical_interfaces
    integer, intent(in) :: n_wavelength_bins

    allocate( field )
    allocate( field%edr_( n_vertical_interfaces, n_wavelength_bins ) )
    allocate( field%edn_( n_vertical_interfaces, n_wavelength_bins ) )
    allocate( field%eup_( n_vertical_interfaces, n_wavelength_bins ) )
    allocate( field%fdr_( n_vertical_interfaces, n_wavelength_bins ) )
    allocate( field%fdn_( n_vertical_interfaces, n_wavelength_bins ) )
    allocate( field%fup_( n_vertical_interfaces, n_wavelength_bins ) )

  end function field_constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine apply_scale_factor( this, scale_factor )
    ! Applies a scaling factor to the radiation field data members

    class(radiation_field_t), intent(inout) :: this
    real(dk),                 intent(in)    :: scale_factor

    this%edr_ = this%edr_ * scale_factor
    this%edn_ = this%edn_ * scale_factor
    this%eup_ = this%eup_ * scale_factor
    this%fdr_ = this%fdr_ * scale_factor
    this%fdn_ = this%fdn_ * scale_factor
    this%fup_ = this%fup_ * scale_factor

  end subroutine apply_scale_factor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function field_pack_size( this, comm ) result( pack_size )
    ! Returns the size of a character buffer required to pack the radiation
    ! field

    use musica_mpi,                    only : musica_mpi_pack_size

    class(radiation_field_t), intent(in) :: this ! radiation field state to be packed
    integer,                  intent(in) :: comm ! MPI communicator

    pack_size = 0

#ifdef MUSICA_USE_MPI
    pack_size = musica_mpi_pack_size( this%edr_, comm) +                      &
                musica_mpi_pack_size( this%eup_, comm) +                      &
                musica_mpi_pack_size( this%edn_, comm) +                      &
                musica_mpi_pack_size( this%fdr_, comm) +                      &
                musica_mpi_pack_size( this%fup_, comm) +                      &
                musica_mpi_pack_size( this%fdn_, comm)
#endif

  end function field_pack_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine field_pack( this, buffer, position, comm )
    ! Packs the radiation field onto a character buffer

    use musica_assert,                 only : assert
    use musica_mpi,                    only : musica_mpi_pack

    class(radiation_field_t), intent(in)    :: this     ! radiation field state to be packed
    character,               intent(inout) :: buffer(:) ! memory buffer
    integer,                 intent(inout) :: position  ! current buffer position
    integer,                 intent(in)    :: comm      ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: prev_pos

    prev_pos = position
    call musica_mpi_pack( buffer, position, this%edr_, comm)
    call musica_mpi_pack( buffer, position, this%eup_, comm)
    call musica_mpi_pack( buffer, position, this%edn_, comm)
    call musica_mpi_pack( buffer, position, this%fdr_, comm)
    call musica_mpi_pack( buffer, position, this%fup_, comm)
    call musica_mpi_pack( buffer, position, this%fdn_, comm)

    call assert( 942613664, position - prev_pos <= this%pack_size( comm ) )
#endif

  end subroutine field_pack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine field_unpack( this, buffer, position, comm )
    ! Unpacks a radiation field from a character buffer

    use musica_assert,                 only : assert
    use musica_mpi,                    only : musica_mpi_unpack

    class(radiation_field_t), intent(out)   :: this      ! radiation field state to be unpacked
    character,                intent(inout) :: buffer(:) ! memory buffer
    integer,                  intent(inout) :: position  ! current buffer position
    integer,                  intent(in)    :: comm      ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: prev_pos

    prev_pos = position
    call musica_mpi_unpack( buffer, position, this%edr_, comm)
    call musica_mpi_unpack( buffer, position, this%eup_, comm)
    call musica_mpi_unpack( buffer, position, this%edn_, comm)
    call musica_mpi_unpack( buffer, position, this%fdr_, comm)
    call musica_mpi_unpack( buffer, position, this%fup_, comm)
    call musica_mpi_unpack( buffer, position, this%fdn_, comm)
    call assert( 709806189, position - prev_pos <= this%pack_size( comm ) )
#endif

  end subroutine field_unpack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine finalize( this )
    ! Cleans up memory for a radiative transfer solver

    type(radiation_field_t), intent(inout) :: this

    if( allocated( this%edr_ ) ) deallocate( this%edr_ )
    if( allocated( this%eup_ ) ) deallocate( this%eup_ )
    if( allocated( this%edn_ ) ) deallocate( this%edn_ )
    if( allocated( this%fdr_ ) ) deallocate( this%fdr_ )
    if( allocated( this%fup_ ) ) deallocate( this%fup_ )
    if( allocated( this%fdn_ ) ) deallocate( this%fdn_ )

  end subroutine finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  pure real(dk) function slant_optical_depth( i_layer, n_layers_crossed,      &
      slant_path, optical_depth )
    ! Calculates the total slant-path optical depth at a given vertical layer

    integer, intent(in) :: i_layer ! Index of the interface to be calculated for (0 = top of the atmosphere, increasing with decreasing altitude)
    integer, intent(in) :: n_layers_crossed ! number of layers crossed by the direct beam when travelling from the top of the atmosphere to layer ``i_layer``
    real(dk), intent(in) :: slant_path(:) ! slant path of direct beam through each layer crossed when travelling from the top of the atmosphere to layer ``i_layer``
    real(dk), intent(in) :: optical_depth(:) ! scaled optical depth in layer ``i_layer``

    integer :: j_layer

    if( n_layers_crossed < 0 ) then
      slant_optical_depth = 1.0e36
      return
    end if
    slant_optical_depth = 0.0_dk
    do j_layer = 1, min( n_layers_crossed, i_layer )
      slant_optical_depth = slant_optical_depth +                             &
                            optical_depth( j_layer ) * slant_path( j_layer )
    end do
    do j_layer = min( n_layers_crossed, i_layer ) + 1, n_layers_crossed
      slant_optical_depth = slant_optical_depth +                             &
                   2.0_dk * optical_depth( j_layer ) * slant_path( j_layer )
    end do

  end function slant_optical_depth

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_solver
