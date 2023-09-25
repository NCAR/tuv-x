! Copyright (C) 2022 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
module tuvx_radiator_from_host
  ! Radiator whose optical properties are provided by the host application
  ! at runtime. See :ref:`configuration-radiators-from-host` for more
  ! information.

  use musica_constants,                only : dk => musica_dk
  use tuvx_radiator,                   only : radiator_t

  implicit none

  public :: radiator_from_host_t, radiator_updater_t

  type, extends(radiator_t) :: radiator_from_host_t
    ! radiator that can be updated from a host application
  contains
    procedure :: update_state
  end type radiator_from_host_t

  ! radiator constructor
  interface radiator_from_host_t
    module procedure :: constructor_char
    module procedure :: constructor_string
  end interface radiator_from_host_t

  type :: radiator_updater_t
#ifndef MUSICA_IS_NAG_COMPILER
    private
#endif
    ! updater for `radiator_from_host_t` radiators
    class(radiator_from_host_t), pointer :: radiator_ => null( )
  contains
    procedure :: update
  end type radiator_updater_t

  ! radiator updater constructor
  interface radiator_updater_t
    module procedure :: updater_constructor
  end interface radiator_updater_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function constructor_char( name, height_grid, wavelength_grid )             &
      result( this )
    ! Constructs a radiator that can be updated by the host application

    use tuvx_grid,                     only : grid_t

    character(len=*),           intent(in) :: name            ! name of the radiator
    class(grid_t),              intent(in) :: height_grid     ! height grid
    class(grid_t),              intent(in) :: wavelength_grid ! wavelength grid
    type(radiator_from_host_t), pointer    :: this            ! new radiator

    allocate( this )

    this%handle_                 = name
    this%type_                   = "from host"
    this%vertical_profile_name_  = "none"
    this%vertical_profile_units_ = "none"
    this%cross_section_name_     = "none"
    this%enable_diagnostics_     = .false.

    allocate( this%state_%layer_OD_(  height_grid%size( ),                    &
                                      wavelength_grid%size( ) ) )
    allocate( this%state_%layer_SSA_( height_grid%size( ),                    &
                                      wavelength_grid%size( ) ) )
    allocate( this%state_%layer_G_(   height_grid%size( ),                    &
                                      wavelength_grid%size( ), 1 ) )

    this%state_%layer_OD_( :,:) = 0.0_dk
    this%state_%layer_SSA_(:,:) = 0.0_dk
    this%state_%layer_G_(:,:,:) = 0.0_dk

  end function constructor_char

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function constructor_string( name, height_grid, wavelength_grid )           &
      result( this )
    ! Constructs a radiator that can be updated by the host application

    use musica_string,                 only : string_t
    use tuvx_grid,                     only : grid_t

    type(string_t),             intent(in) :: name            ! name of the radiator
    class(grid_t),              intent(in) :: height_grid     ! height grid
    class(grid_t),              intent(in) :: wavelength_grid ! wavelength grid
    type(radiator_from_host_t), pointer    :: this            ! new radiator

    this => constructor_char( name%to_char( ), height_grid, wavelength_grid )

  end function constructor_string

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function updater_constructor( radiator ) result( this )
    ! Constructs an updater for a `radiator_from_host_t` radiator

    class(radiator_from_host_t), target, intent(inout) :: radiator ! radiator to be updated
    type(radiator_updater_t)                           :: this     ! new updater

    this%radiator_ => radiator

  end function updater_constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine update( this, optical_depths, single_scattering_albedos,         &
      asymmetry_factors )
    ! Updates the radiator optical properties
    !
    ! If single scattering albedos or asymmetry factors are omitted, they are
    ! assumed to be zero.

    use musica_assert,                 only : die_msg, assert_msg
    use musica_string,                 only : to_char

    class(radiator_updater_t), intent(inout) :: this ! radiator updater
    real(kind=dk),             intent(in)    :: optical_depths(:,:) ! optical depths (vertical layer, wavelength bin) [TODO units?]
    real(kind=dk), optional,   intent(in)    :: single_scattering_albedos(:,:) ! single scattering albedos (vertical layer, wavelength bin) [TODO units?]
    real(kind=dk), optional,   intent(in)    :: asymmetry_factors(:,:) ! asymetry factors (vertical level, wavelength bin) [TODO units?]

    integer :: n_levels_host, n_wavelengths_host
    integer :: n_levels, n_wavelengths

    call assert_msg( 513442250, associated( this%radiator_ ),                 &
                     "Cannot update an unspecified radiator" )
    associate( radiator => this%radiator_, state => this%radiator_%state_ )

      ! optical depths
      n_levels           = size( state%layer_OD_, 1 )
      n_wavelengths      = size( state%layer_OD_, 2 )
      n_levels_host      = size( optical_depths,  1 )
      n_wavelengths_host = size( optical_depths,  2 )
      if( n_levels .ne. n_levels_host ) then
        call die_msg( 366082936, "Size mismatch for radiator '"//             &
                      trim( radiator%handle_%to_char( ) )//                   &
                      "' optical depth levels. Expected "//                   &
                      trim( to_char( n_levels ) )//", got "//                 &
                      trim( to_char( n_levels_host ) ) )
      end if
      if( n_wavelengths .ne. n_wavelengths_host ) then
        call die_msg( 770637729, "Size mismatch for radiator '"//             &
                      trim( radiator%handle_%to_char( ) )//                   &
                      "' optical depth wavelengths. Expected "//              &
                      trim( to_char( n_wavelengths ) )//", got "//            &
                      trim( to_char( n_wavelengths_host ) ) )
      end if
      state%layer_OD_(:,:) = optical_depths(:,:)

      ! single scattering albedos
      if( present( single_scattering_albedos ) ) then
        n_levels           = size( state%layer_SSA_, 1 )
        n_wavelengths      = size( state%layer_SSA_, 2 )
        n_levels_host      = size( single_scattering_albedos, 1 )
        n_wavelengths_host = size( single_scattering_albedos, 2 )
        if( n_levels .ne. n_levels_host ) then
          call die_msg( 517225321, "Size mismatch for radiator '"//           &
                        trim( radiator%handle_%to_char( ) )//                 &
                        "' single scattering albedo levels. Expected "//      &
                        trim( to_char( n_levels ) )//", got "//               &
                        trim( to_char( n_levels_host ) ) )
        end if
        if( n_wavelengths .ne. n_wavelengths_host ) then
          call die_msg( 629543666, "Size mismatch for radiator '"//           &
                        trim( radiator%handle_%to_char( ) )//                 &
                        "' single scattering albedo wavelengths. Expected "// &
                        trim( to_char( n_wavelengths ) )//", got "//          &
                        trim( to_char( n_wavelengths_host ) ) )
        end if
        state%layer_SSA_(:,:) = single_scattering_albedos(:,:)
      else
        state%layer_SSA_(:,:) = 0.0_dk
      end if

      ! aymmetry factros
      if( present( asymmetry_factors ) ) then
        n_levels           = size( state%layer_G_, 1 )
        n_wavelengths      = size( state%layer_G_, 2 )
        n_levels_host      = size( asymmetry_factors, 1 )
        n_wavelengths_host = size( asymmetry_factors, 2 )
        if( n_levels .ne. n_levels_host ) then
          call die_msg( 789624251, "Size mismatch for radiator '"//           &
                        trim( radiator%handle_%to_char( ) )//                 &
                        "' asymmetry factor levels. Expected "//              &
                        trim( to_char( n_levels ) )//", got "//               &
                        trim( to_char( n_levels_host ) ) )
        end if
        if( n_wavelengths .ne. n_wavelengths_host ) then
          call die_msg( 336992098, "Size mismatch for radiator '"//           &
                        trim( radiator%handle_%to_char( ) )//                 &
                        "' asymmetry factor wavelengths. Expected "//         &
                        trim( to_char( n_wavelengths ) )//", got "//          &
                        trim( to_char( n_wavelengths_host ) ) )
        end if
        state%layer_G_(:,:,1) = asymmetry_factors(:,:)
      else
        state%layer_G_(:,:,:) = 0.0_dk
      end if

    end associate

  end subroutine update

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine update_state( this, grid_warehouse, profile_warehouse,           &
      cross_section_warehouse )
    ! No updates based on profile state are done for this radiation type

    use tuvx_cross_section_warehouse,  only : cross_section_warehouse_t
    use tuvx_grid_warehouse,           only : grid_warehouse_t
    use tuvx_profile_warehouse,        only : profile_warehouse_t

    class(radiator_from_host_t),     intent(inout) :: this
    type(grid_warehouse_t),          intent(inout) :: grid_warehouse
    type(profile_warehouse_t),       intent(inout) :: profile_warehouse
    type(cross_section_warehouse_t), intent(inout) :: cross_section_warehouse

  end subroutine update_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_radiator_from_host
