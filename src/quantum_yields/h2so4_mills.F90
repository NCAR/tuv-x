! Copyright (C) 2024 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0

module tuvx_quantum_yield_h2so4_mills
  ! The temperature and pressure dependent quantum yield calculations
  ! used in WACCM simulations

  use musica_constants,                only : dk => musica_dk
  use tuvx_quantum_yield

  implicit none

  private
  public :: quantum_yield_h2so4_mills_t

  !> Quantum yield calculator for H2SO4
  !!
  !! See Miller et al. (GRL, 2007)
  !!
  !! Quantum yields qy_x are based on parameters r_x set in the
  !! configuration:
  !!
  !!  lambda = R T / ( 2^(1/2) pi d^2 Na P )
  !!  v = ( 8 R T / ( pi MW ) )^(1/2)
  !!  qy_x = r_x / ( r_x + v / lambda )
  !!
  !! where R is the universal gas constant (J mol-1 K-1 ), T is
  !! temperature (K), P is pressure (Pa), Na is  Avogadro's number (mol-1)
  !! and MW is the molecular weight of H2SO4 (kg mol-1), d is the
  !! molecular diameter (m), and x corresponds to the wavelength band
  !! being parameterized.
  type, extends(quantum_yield_t) :: quantum_yield_h2so4_mills_t
    !> Indices for wavelengths to update
    integer, allocatable :: wavelength_indices_(:)
    !> Collision rate [s-1]
    real(kind=dk), allocatable :: collision_rate_(:)
    !> Molecular diameter [m]
    real(kind=dk) :: molecular_diameter_
    !> Molecular weight [kg mol-1]
    real(kind=dk) :: molecular_weight_
  contains
    !> Calculate the quantum yields
    procedure :: calculate
    ! returns the number of bytes required to pack the object onto a buffer
    procedure :: pack_size
    ! packs the object onto a character buffer
    procedure :: mpi_pack
    ! unpacks an object from a character buffer
    procedure :: mpi_unpack
  end type quantum_yield_h2so4_mills_t

  interface quantum_yield_h2so4_mills_t
    module procedure :: constructor
  end interface quantum_yield_h2so4_mills_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Constructor of H2SO4 quantum yield calculators
  function constructor( config, grids, profiles ) result ( this )

    use musica_assert,                 only : assert_msg
    use musica_config,                 only : config_t
    use musica_string,                 only : string_t, to_char
    use tuvx_grid,                     only : grid_t
    use tuvx_grid_warehouse,           only : grid_warehouse_t
    use tuvx_profile_warehouse,        only : profile_warehouse_t

    class(quantum_yield_h2so4_mills_t), pointer       :: this
    type(config_t),                     intent(inout) :: config
    type(grid_warehouse_t),             intent(inout) :: grids
    type(profile_warehouse_t),          intent(inout) :: profiles

    character(len=*), parameter :: my_name = "H2SO4 quantum yield constructor"
    type(string_t) :: required_keys(5), optional_keys(6)
    real(kind=dk), allocatable :: param_wavelengths(:)
    type(grid_t), pointer :: wavelengths
    integer :: i_param, i_wl

    required_keys(1) = "type"
    required_keys(2) = "parameterized wavelengths"
    required_keys(3) = "collision interval s"
    required_keys(4) = "molecular diameter m"
    required_keys(5) = "molecular weight kg mol-1"
    optional_keys(1) = "netcdf files"
    optional_keys(2) = "lower extrapolation"
    optional_keys(3) = "upper extrapolation"
    optional_keys(4) = "name"
    optional_keys(5) = "constant value"
    optional_keys(6) = "override bands"
    call assert_msg( 157064056,                                               &
                     config%validate( required_keys, optional_keys ),         &
                     "Bad configration data format for H2SO4 quantum yield." )
    allocate( this )
    call base_constructor( this, config, grids, profiles )

    call config%get( "parameterized wavelengths", param_wavelengths, my_name )
    call config%get( "collision interval s", this%collision_rate_,       &
                     my_name )
    this%collision_rate_(:) = 1.0_dk / this%collision_rate_(:)
    call config%get( "molecular diameter m", this%molecular_diameter_,      &
                     my_name )
    call config%get( "molecular weight kg mol-1", this%molecular_weight_,   &
                     my_name )
    call assert_msg( 472700337, size( param_wavelengths ) .eq.                &
                                size( this%collision_rate_ ),                 &
                     "Size mismatch between parameterized wavelengths and "// &
                     "collision frequency in H2SO4 quantum yield calculator" )    
    wavelengths => grids%get_grid( "wavelength", "nm" )
    allocate( this%wavelength_indices_( size( param_wavelengths ) ) )
    this%wavelength_indices_(:) = 0
    do i_param = 1, size( param_wavelengths )
      do i_wl = 1, wavelengths%ncells_
        if( wavelengths%mid_( i_wl ) .eq. param_wavelengths( i_param ) ) then
          this%wavelength_indices_( i_param ) = i_wl
          exit
        end if
      end do
      call assert_msg( 170811868, this%wavelength_indices_( i_param ) > 0,    &
                       "Parameterized wavelength in H2SO4 quantum yield "//   &
                       "configuration not on wavelength grid: "//             &
                       trim( to_char( param_wavelengths( i_param ) ) ) )
    end do
    deallocate( wavelengths )

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculates the quantum yield
  function calculate( this, grid_warehouse, profile_warehouse )               &
      result( quantum_yield )

    use tuvx_constants,                only : gas_constant, Avogadro, pi
    use tuvx_grid,                     only : grid_t
    use tuvx_grid_warehouse,           only : grid_warehouse_t
    use tuvx_profile,                  only : profile_t
    use tuvx_profile_warehouse,        only : profile_warehouse_t

    class(quantum_yield_h2so4_mills_t), intent(in)    :: this
    type(grid_warehouse_t),             intent(inout) :: grid_warehouse
    type(profile_warehouse_t),          intent(inout) :: profile_warehouse
    real(dk), allocatable                             :: quantum_yield(:,:)

    class(profile_t), pointer :: temperature, air
    integer :: i_wl
    real(dk) :: lambda, velocity

    quantum_yield =                                                           &
        this%quantum_yield_t%calculate( grid_warehouse, profile_warehouse )
    temperature => profile_warehouse%get_profile( this%temperature_profile_ )
    air => profile_warehouse%get_profile( this%air_profile_ )

    ! Overwrite the quantum yields for the parameterized wavelengths
    do i_wl = 1, size( this%wavelength_indices_ )
      quantum_yield( :, this%wavelength_indices_( i_wl ) ) =                  &
          this%collision_rate_( i_wl ) /                                      &
            ( this%collision_rate_( i_wl ) +                                  &
              sqrt( 16.0_dk * gas_constant * temperature%edge_val_(:)         &
                    / ( pi * this%molecular_weight_ ) )                       &
              * ( pi * this%molecular_diameter_**2 *                          &
                  air%edge_val_(:) * 1.0e6_dk ) )
    end do
    ! The top layer has quantum yields set to 1.0
    quantum_yield( size( quantum_yield, dim=1 ),                              &
                   this%wavelength_indices_(:) ) = 1.0_dk

    deallocate( temperature )
    deallocate( air         )

  end function calculate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the number of bytes required to pack the object onto a buffer
  integer function pack_size( this, comm )

    use musica_mpi,                    only : musica_mpi_pack_size

    !> Quantum yield to be packed
    class(quantum_yield_h2so4_mills_t), intent(in) :: this
    !> MPI communicator
    integer,                            intent(in) :: comm

#ifdef MUSICA_USE_MPI
    pack_size = this%quantum_yield_t%pack_size( comm ) +                      &
                musica_mpi_pack_size( this%wavelength_indices_, comm ) +      &
                musica_mpi_pack_size( this%collision_rate_,     comm ) +      &
                musica_mpi_pack_size( this%molecular_diameter_, comm ) +      &
                musica_mpi_pack_size( this%molecular_weight_,   comm )
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
    class(quantum_yield_h2so4_mills_t), intent(in)    :: this
    !> Memory buffer
    character,                          intent(inout) :: buffer(:)
    !> Current buffer position
    integer,                            intent(inout) :: position
    !> MPI communicator
    integer,                            intent(in)    :: comm

#ifdef MUSICA_USE_MPI
    integer :: prev_pos

    prev_pos = position
    call this%quantum_yield_t%mpi_pack( buffer, position, comm )
    call musica_mpi_pack( buffer, position, this%wavelength_indices_, comm )
    call musica_mpi_pack( buffer, position, this%collision_rate_,     comm )
    call musica_mpi_pack( buffer, position, this%molecular_diameter_, comm )
    call musica_mpi_pack( buffer, position, this%molecular_weight_,   comm )
    call assert( 931898871, position - prev_pos <= this%pack_size( comm ) )
#endif

  end subroutine mpi_pack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Unpacks a quantum yield calculator from a character buffer
  subroutine mpi_unpack( this, buffer, position, comm )

    use musica_assert,                 only : assert
    use musica_mpi,                    only : musica_mpi_unpack

    !> Quantum yield to unpack
    class(quantum_yield_h2so4_mills_t), intent(out)   :: this
    !> Memory buffer
    character,                          intent(inout) :: buffer(:)
    !> Current buffer position
    integer,                            intent(inout) :: position
    !> MPI communicator
    integer,                            intent(in)    :: comm

#ifdef MUSICA_USE_MPI
    integer :: prev_pos

    prev_pos = position
    call this%quantum_yield_t%mpi_unpack( buffer, position, comm )
    call musica_mpi_unpack( buffer, position, this%wavelength_indices_, comm )
    call musica_mpi_unpack( buffer, position, this%collision_rate_,     comm )
    call musica_mpi_unpack( buffer, position, this%molecular_diameter_, comm )
    call musica_mpi_unpack( buffer, position, this%molecular_weight_,   comm )
    call assert( 237836163, position - prev_pos <= this%pack_size( comm ) )
#endif

  end subroutine mpi_unpack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_quantum_yield_h2so4_mills