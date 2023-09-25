! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
module tuvx_dose_rates
  ! The dose_rates_t type and related functions

  use musica_constants,                only : dk => musica_dk
  use musica_string,                   only : string_t
  use tuvx_spectral_weight,            only : spectral_weight_ptr
  use tuvx_grid_warehouse,             only : grid_warehouse_ptr
  use tuvx_grid,                       only : grid_t
  use tuvx_profile_warehouse,          only : profile_warehouse_ptr
  use tuvx_profile,                    only : profile_t

  implicit none

  private
  public :: dose_rates_t

  !> Photolysis rate constant calculator
  type :: dose_rates_t
    ! Spectral weights
    type(spectral_weight_ptr), allocatable :: spectral_weights_(:)
    ! Configuration label for the dose rate
    type(string_t),            allocatable :: handles_(:)
    ! Flag for enabling diagnostic output
    logical                                :: enable_diagnostics_
    ! Height grid
    type(grid_warehouse_ptr) :: height_grid_
    ! Wavelength grid
    type(grid_warehouse_ptr) :: wavelength_grid_
    ! Extraterrestrial flux profile
    type(profile_warehouse_ptr) :: etfl_profile_
  contains
    ! Returns the dose rates for a given set of conditions
    procedure :: get
    ! Returns the names of each dose rate
    procedure :: labels
    ! Returns the number of dose rates
    procedure :: size => get_number
    ! Returns the number of bytes required to pack the dose rates onto a
    ! buffer
    procedure :: pack_size
    ! Packs the dose rates onto a character buffer
    procedure :: mpi_pack
    ! Unpacks dose rates from a character buffer
    procedure :: mpi_unpack
    ! Finalize the object
    final :: finalize
  end type dose_rates_t

  !> dose_rates_t constructor
  interface dose_rates_t
    module procedure :: constructor
  end interface

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function constructor( config, grid_warehouse, profile_warehouse )           &
      result( dose_rates )
    ! Constructor of dose_rates_t objects

    use musica_assert,                 only : assert_msg
    use musica_config,                 only : config_t
    use musica_iterator,               only : iterator_t
    use tuvx_grid_warehouse,           only : grid_warehouse_t
    use tuvx_spectral_weight_factory,  only : spectral_weight_builder
    use tuvx_profile_warehouse,        only : profile_warehouse_t

    !> Dose rate configuration
    type(config_t),            intent(inout) :: config
    !> grid warehouse
    type(grid_warehouse_t),    intent(inout) :: grid_warehouse
    !> profile warehouse
    type(profile_warehouse_t), intent(inout) :: profile_warehouse
    !> New dose rates
    class(dose_rates_t),        pointer      :: dose_rates

    ! Local variables
    character(len=*), parameter :: Iam = "dose_rates_t constructor"
    type(config_t) :: wght_config, spectral_weight_config
    class(iterator_t), pointer  :: iter
    type(spectral_weight_ptr)   :: spectral_weight
    character(len=64)           :: keychar
    type(string_t)              :: wght_key
    type(string_t) :: required_keys(1), optional_keys(1)
    type(string_t) :: rate_required_keys(2), rate_optional_keys(0)
    type(config_t) :: rate_config

    required_keys(1) = "rates"
    optional_keys(1) = "enable diagnostics"

    rate_required_keys(1) = "name"
    rate_required_keys(2) = "weights"

    call assert_msg( 100983245,                                               &
                     config%validate( required_keys, optional_keys ),         &
                     "Bad configuration data format for "//                   &
                     "dose rates." )

    allocate( dose_rates )

    associate( rates => dose_rates )

    allocate( string_t :: rates%handles_(0) )
    allocate( rates%spectral_weights_(0) )

    call config%get( "enable diagnostics", rates%enable_diagnostics_, Iam,    &
                     default = .false. )

    rates%height_grid_ = grid_warehouse%get_ptr( "height", "km" )
    rates%wavelength_grid_ = grid_warehouse%get_ptr( "wavelength", "nm" )
    rates%etfl_profile_ = profile_warehouse%get_ptr( "extraterrestrial flux", &
                                                     "photon cm-2 s-1" )

    ! iterate over dose rates
    call config%get( "rates", rate_config, Iam )


    iter => rate_config%get_iterator( )
    do while( iter%next( ) )
      call rate_config%get( iter, wght_config, Iam )

      call assert_msg( 200983245,                                               &
        wght_config%validate( rate_required_keys, rate_optional_keys ),         &
        "Bad configuration data format for dose rates")

      ! get spectral wght
      call wght_config%get( "weights", spectral_weight_config, Iam )
      call wght_config%get( "name", wght_key, Iam )

      rates%handles_ = [ rates%handles_, wght_key ]
      call config%get( iter, wght_config, Iam )

      spectral_weight%val_ => &
         spectral_weight_builder( spectral_weight_config, grid_warehouse,     &
                                  profile_warehouse )
      rates%spectral_weights_ = [ rates%spectral_weights_,                    &
                                  spectral_weight ]
    end do
    deallocate( iter )

    end associate

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine get( this, grid_warehouse, profile_warehouse, radiation_field,   &
      dose_rates, file_tag )
    ! calculates dose rate constants

    use musica_assert,                 only : assert_msg
    use tuvx_diagnostic_util,          only : diagout
    use tuvx_grid_warehouse,           only : grid_warehouse_t
    use tuvx_solver,                   only : radiation_field_t
    use tuvx_profile_warehouse,        only : profile_warehouse_t

    !> Dose rate constant calculator
    class(dose_rates_t),       intent(inout) :: this
    !> Warehouses
    type(grid_warehouse_t),    intent(inout) :: grid_warehouse
    type(profile_warehouse_t), intent(inout) :: profile_warehouse
    !> Actinic flux
    type(radiation_field_t),   intent(in)    :: radiation_field
    !> Tag used in file name of output data
    character(len=*),          intent(in)    :: file_tag
    !> Calculated dose rate constants (vertical layer, dose rate type)
    real(dk),                  intent(inout) :: dose_rates(:,:)

    ! Local variables
    character(len=*), parameter :: Iam = "dose rates calculator:"
    integer               :: wavNdx, rateNdx, nRates
    real(dk), allocatable :: spectral_weight(:)
    real(dk), allocatable :: tmp_spectral_weight(:)
    real(dk), allocatable :: sirrad(:,:)
    class(grid_t),    pointer :: zGrid
    class(grid_t),    pointer :: lambdaGrid
    class(profile_t), pointer :: etfl

    zGrid => grid_warehouse%get_grid( this%height_grid_ )
    lambdaGrid => grid_warehouse%get_grid( this%wavelength_grid_ )
    etfl => profile_warehouse%get_profile( this%etfl_profile_ )

    if( this%enable_diagnostics_ ) then
      allocate( tmp_spectral_weight(0) )
    end if

    nRates = size( this%spectral_weights_ )
    call assert_msg( 116265931,                                               &
                     size( dose_rates, 1 ) == zGrid%ncells_ + 1 .and.         &
                     size( dose_rates, 2 ) == nRates ,                        &
                     "Bad shape for dose rates array" )

    !> spectral irradiance
    sirrad = radiation_field%edr_ + radiation_field%eup_ + radiation_field%edn_
    do wavNdx = 1, lambdaGrid%ncells_
      sirrad( :, wavNdx ) = sirrad( :, wavNdx ) * etfl%mid_val_( wavNdx )
    enddo
    where( sirrad < 0.0_dk )
      sirrad = 0.0_dk
    end where

rate_loop:                                                                    &
    do rateNdx = 1, nRates
      associate( calc_ftn => this%spectral_weights_( rateNdx )%val_ )
        spectral_weight = calc_ftn%calculate( grid_warehouse,                 &
                                              profile_warehouse )
      end associate

      if( this%enable_diagnostics_ ) then
        tmp_spectral_weight = [ tmp_spectral_weight, spectral_weight ]
      end if
      dose_rates( :, rateNdx ) = matmul( sirrad, spectral_weight )

      if( allocated( spectral_weight ) ) deallocate( spectral_weight )
    end do rate_loop

    if( this%enable_diagnostics_ ) then
      call diagout( 'annotatedslabels.new', this%handles_,                      &
        this%enable_diagnostics_ )
      call diagout( 'sw.'//file_tag//'.new', tmp_spectral_weight,               &
        this%enable_diagnostics_ )
    end if

    if( associated( zGrid ) ) deallocate( zGrid )
    if( associated( lambdaGrid ) ) deallocate( lambdaGrid )
    if( associated( etfl ) ) deallocate( etfl )

  end subroutine get

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function labels( this )
    ! Returns the names of each dose rate

    type(string_t), allocatable     :: labels(:)
    class(dose_rates_t), intent(in) :: this

    labels = this%handles_

  end function labels

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function get_number( this )
    ! Returns the number of dose rates

    use musica_assert,                 only : assert_msg

    class(dose_rates_t), intent(in) :: this

    call assert_msg( 952537624, allocated( this%handles_ ),                   &
                     "Dose rates not initialized" )
    get_number = size( this%handles_ )

  end function get_number

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function pack_size( this, comm )
    ! Returns the number of bytes required to pack the dose rates onto a
    ! buffer

    use musica_mpi,                    only : musica_mpi_pack_size
    use tuvx_spectral_weight_factory,  only : spectral_weight_type_name

    class(dose_rates_t), intent(in) :: this ! dose rates to be packed
    integer,             intent(in) :: comm ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: i_elem
    type(string_t) :: type_name

    pack_size = musica_mpi_pack_size( allocated( this%spectral_weights_ ),    &
                                      comm )
    if( allocated( this%spectral_weights_ ) ) then
      pack_size = pack_size +                                                 &
                  musica_mpi_pack_size( size( this%spectral_weights_ ), comm )
      do i_elem = 1, size( this%spectral_weights_ )
      associate( weight => this%spectral_weights_( i_elem )%val_ )
        type_name = spectral_weight_type_name( weight )
        pack_size = pack_size +                                               &
                    type_name%pack_size( comm ) +                             &
                    weight%pack_size( comm )
      end associate
      end do
    end if
    pack_size = pack_size +                                                   &
                musica_mpi_pack_size( allocated( this%handles_ ), comm )
    if( allocated( this%handles_ ) ) then
      pack_size = pack_size +                                                 &
                  musica_mpi_pack_size( size( this%handles_  ), comm )
      do i_elem = 1, size( this%handles_ )
        pack_size = pack_size + this%handles_( i_elem )%pack_size( comm )
      end do
    end if
    pack_size = pack_size +                                                   &
                musica_mpi_pack_size( this%enable_diagnostics_, comm ) +      &
                this%height_grid_%pack_size( comm ) +                         &
                this%wavelength_grid_%pack_size( comm ) +                     &
                this%etfl_profile_%pack_size( comm )
#else
    pack_size = 0
#endif

  end function pack_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine mpi_pack( this, buffer, position, comm )
    ! Packs the dose rates onto a character buffer

    use musica_assert,                 only : assert
    use musica_mpi,                    only : musica_mpi_pack
    use tuvx_spectral_weight_factory,  only : spectral_weight_type_name

    class(dose_rates_t), intent(in)    :: this      ! dose rates to be packed
    character,           intent(inout) :: buffer(:) ! memory buffer
    integer,             intent(inout) :: position  ! current buffer position
    integer,             intent(in)    :: comm      ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: prev_pos, i_elem
    type(string_t) :: type_name

    prev_pos = position
    call musica_mpi_pack( buffer, position,                                   &
                          allocated( this%spectral_weights_ ), comm )
    if( allocated( this%spectral_weights_ ) ) then
      call musica_mpi_pack( buffer, position, size( this%spectral_weights_ ), &
                            comm )
      do i_elem = 1, size( this%spectral_weights_ )
      associate( weight => this%spectral_weights_( i_elem )%val_ )
        type_name = spectral_weight_type_name( weight )
        call type_name%mpi_pack( buffer, position, comm )
        call weight%mpi_pack(    buffer, position, comm )
      end associate
      end do
    end if
    call musica_mpi_pack( buffer, position, allocated( this%handles_ ), comm )
    if( allocated( this%handles_ ) ) then
      call musica_mpi_pack( buffer, position, size( this%handles_ ), comm )
      do i_elem = 1, size( this%handles_ )
        call this%handles_( i_elem )%mpi_pack( buffer, position, comm )
      end do
    end if
    call musica_mpi_pack( buffer, position, this%enable_diagnostics_, comm )
    call this%height_grid_%mpi_pack(     buffer, position, comm )
    call this%wavelength_grid_%mpi_pack( buffer, position, comm )
    call this%etfl_profile_%mpi_pack(    buffer, position, comm )
    call assert( 258716172, position - prev_pos <= this%pack_size( comm ) )
#endif

  end subroutine mpi_pack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine mpi_unpack( this, buffer, position, comm )
    ! Unpacks the dose rates onto a character buffer

    use musica_assert,                 only : assert
    use musica_mpi,                    only : musica_mpi_unpack
    use tuvx_spectral_weight_factory,  only : spectral_weight_allocate

    class(dose_rates_t), intent(out)   :: this      ! dose rates to be unpacked
    character,           intent(inout) :: buffer(:) ! memory buffer
    integer,             intent(inout) :: position  ! current buffer position
    integer,             intent(in)    :: comm      ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: prev_pos, i_elem, n_elems
    logical :: alloced
    type(string_t) :: type_name

    prev_pos = position
    call musica_mpi_unpack( buffer, position, alloced, comm )
    if( alloced ) then
      call musica_mpi_unpack( buffer, position, n_elems, comm )
      allocate( this%spectral_weights_( n_elems ) )
      do i_elem = 1, n_elems
      associate( weight => this%spectral_weights_( i_elem ) )
        call type_name%mpi_unpack( buffer, position, comm )
        weight%val_ => spectral_weight_allocate( type_name )
        call weight%val_%mpi_unpack( buffer, position, comm )
      end associate
      end do
    end if
    call musica_mpi_unpack( buffer, position, alloced, comm )
    if( alloced ) then
      call musica_mpi_unpack( buffer, position, n_elems, comm )
      allocate( this%handles_( n_elems ) )
      do i_elem = 1, n_elems
        call this%handles_( i_elem )%mpi_unpack( buffer, position, comm )
      end do
    end if
    call musica_mpi_unpack( buffer, position, this%enable_diagnostics_, comm )
    call this%height_grid_%mpi_unpack(     buffer, position, comm )
    call this%wavelength_grid_%mpi_unpack( buffer, position, comm )
    call this%etfl_profile_%mpi_unpack(    buffer, position, comm )
    call assert( 143039054, position - prev_pos <= this%pack_size( comm ) )
#endif

  end subroutine mpi_unpack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine finalize( this )
    ! Finalize the dose rates object

    !> Dose rate object
    type(dose_rates_t), intent(inout) :: this

    integer :: ndx

    if( allocated( this%spectral_weights_ ) ) then
      do ndx = 1,size( this%spectral_weights_ )
        if( associated( this%spectral_weights_( ndx )%val_ ) ) then
          deallocate( this%spectral_weights_( ndx )%val_ )
        endif
      enddo
      deallocate( this%spectral_weights_ )
    end if

    if( allocated( this%handles_ ) ) deallocate( this%handles_ )

  end subroutine finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_dose_rates
