! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
module tuvx_profile_from_host
  ! Profile whose values will be provided by the host application at runtime
  ! See :ref:`configuration-profiles-from-host` for more information.

  use musica_constants,                only : dk => musica_dk
  use tuvx_profile,                    only : profile_t

  implicit none

  public :: profile_from_host_t, profile_updater_t

  type, extends(profile_t) :: profile_from_host_t
    ! profile that can be updated from a host application
  contains
  end type profile_from_host_t

  ! profile constructor
  interface profile_from_host_t
    module procedure :: constructor_char
    module procedure :: constructor_string
  end interface profile_from_host_t

  type :: profile_updater_t
#ifndef MUSICA_IS_NAG_COMPILER
    private
#endif
    ! updater for `profile_from_host_t` profiles
    class(profile_from_host_t), pointer :: profile_ => null( )
  contains
    procedure :: update
  end type profile_updater_t

  ! profile updater constructor
  interface profile_updater_t
    module procedure :: updater_constructor
  end interface profile_updater_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function constructor_char( name, units, number_of_cells ) result( this )
    ! Initializes the profile

    use musica_assert,                 only : assert_msg
    use musica_string,                 only : to_char

    character(len=*),          intent(in) :: name            ! name of the profile
    character(len=*),          intent(in) :: units           ! units for the profile
    integer,                   intent(in) :: number_of_cells ! number of discreet profile cells
    type(profile_from_host_t), pointer    :: this            ! constructor profile

    allocate( this )

    this%handle_ = name
    this%units_  = units
    this%ncells_ = number_of_cells
    this%hscale_ = 0.0_dk
    this%enable_diagnostics = .false.

    call assert_msg( 253742246, this%ncells_ >= 0,                            &
                     "Invalid profile size for profile from host: '"//        &
                     trim( to_char( number_of_cells ) ) )
    allocate( this%mid_val_(        this%ncells_     ) )
    allocate( this%edge_val_(       this%ncells_ + 1 ) )
    allocate( this%delta_val_(      this%ncells_     ) )
    allocate( this%layer_dens_(     this%ncells_     ) )
    allocate( this%exo_layer_dens_( this%ncells_ + 1 ) )
    allocate( this%burden_dens_(    this%ncells_     ) )
    this%mid_val_(:)        = 0.0_dk
    this%edge_val_(:)       = 0.0_dk
    this%delta_val_(:)      = 0.0_dk
    this%layer_dens_(:)     = 0.0_dk
    this%exo_layer_dens_(:) = 0.0_dk
    this%burden_dens_(:)    = 0.0_dk

  end function constructor_char

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function constructor_string( name, units, number_of_cells ) result( this )
    ! Initializes the profile

    use musica_string,                 only : string_t

    type(string_t),            intent(in) :: name            ! name of the profile
    type(string_t),            intent(in) :: units           ! units for the profile
    integer,                   intent(in) :: number_of_cells ! number of discreet profile cells
    type(profile_from_host_t), pointer    :: this            ! constructed profile

    this => constructor_char( name%to_char( ), units%to_char( ),              &
                              number_of_cells )

  end function constructor_string

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function updater_constructor( profile ) result( this )
    ! Constructs an updater for a `profile_from_host_t` profile

    class(profile_from_host_t), target, intent(inout) :: profile ! profile to be updated
    type(profile_updater_t)                           :: this    ! new updater

    this%profile_ => profile

  end function updater_constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine update( this, mid_point_values, edge_values, layer_densities,    &
     scale_height, exo_density )
    ! Updates the target profile
    !
    ! Profiles are updated according to which optional values are passed to
    ! this function as follows:
    !
    ! - edge_values must always be provided
    ! - if mid_point_values are provided they are used, otherwise they are
    !   calculated as the average of the cell edge values
    ! - if layer_densities are not provided, then layer_densities,
    !   exo_density and burden_densities are not changed
    ! - if layer_densities are provided, burden_densities are calculated as
    !   the cumulative layer density moving from the upper to the lower
    !   grid boundary. If an exo_density is not provided, the scale_height
    !   is used to calculate the exo_density.
    ! - if a scale_height is provided, layer_densities must be provided,
    !   and exo_density must NOT be provided, The exo density is calculated
    !   using the scale height. If a scale_height is not provided it is
    !   assumed to be zero.
    ! - if an exo_density is provided, layer_densities must be provided,
    !   and scale_height must NOT be provided.
    ! - any provided or calculated exo density is added to the topmost
    !   layer of the layer_densities array

    use musica_assert,                 only : die_msg, assert_msg
    use musica_string,                 only : to_char

    class(profile_updater_t), intent(inout) :: this                ! profile updater
    real(kind=dk), optional,  intent(in)    :: mid_point_values(:) ! new mid-point values
    real(kind=dk), optional,  intent(in)    :: edge_values(:)      ! new edge values
    real(kind=dk), optional,  intent(in)    :: layer_densities(:)  ! new layer densities
    real(kind=dk), optional,  intent(in)    :: scale_height        ! new scale height [km]
    real(kind=dk), optional,  intent(in)    :: exo_density         ! density beyond the upper grid boundary

    integer :: size_profile, size_host, i_elem
    real(kind=dk) :: accum

    call assert_msg( 314048995, associated( this%profile_ ),                  &
                     "Cannot update an unspecified profile" )
    associate( profile => this%profile_ )

      ! mid points provided
      if( present( mid_point_values ) ) then
        size_profile = size( profile%mid_val_ )
        size_host    = size( mid_point_values )
        if( size_profile .ne. size_host ) then
          call die_msg( 307331449,                                            &
                        "Size mismatch for profile mid-point values for "//   &
                        "profile '"//profile%handle_//"'. Expected "//        &
                        trim( to_char( size_profile ) )//", got "//           &
                        trim( to_char( size_host ) ) )
        end if
        profile%mid_val_(:) = mid_point_values(:)
      end if

      ! edge values povided
      call assert_msg( 829266883, present( edge_values ),                     &
                       "Edge values must always be provided to a profile "//  &
                       "updater." )
      size_profile = size( profile%edge_val_ )
      size_host    = size( edge_values )
      if( size_profile .ne. size_host ) then
        call die_msg( 573012833,                                              &
                      "Size mismatch for profile edge values for "//          &
                      "profile '"//profile%handle_//"'. Expected "//          &
                      trim( to_char( size_profile ) )//", got "//             &
                      trim( to_char( size_host ) ) )
      end if
      profile%edge_val_(:) = edge_values(:)
      if( profile%ncells_ > 0 ) then
        profile%delta_val_(:) =                                               &
            profile%edge_val_( 2 : profile%ncells_ + 1 ) -                    &
            profile%edge_val_( 1 : profile%ncells_ )
        if( .not. present( mid_point_values ) ) then
          profile%mid_val_(:) =                                               &
              ( profile%edge_val_( 1 : profile%ncells_ ) +                    &
                profile%edge_val_( 2 : profile%ncells_ + 1 ) ) * 0.5_dk
        end if
      end if

      ! scale height provided
      if( present( scale_height ) ) then
        call assert_msg( 700154673, present( layer_densities ),               &
                         "Layer densities must be provded if scale height "// &
                         "is specified." )
        call assert_msg( 191589495, .not. present( exo_density ),             &
                         "Exo denisty and scale height cannot both be "//     &
                         "provided." )
        profile%hscale_ = scale_height
      else
        profile%hscale_ = 0.0_dk
      end if

      ! layer densities provided
      if( present( layer_densities ) ) then
        size_profile = size( profile%layer_dens_ )
        size_host    = size( layer_densities )
        if( size_profile .ne. size_host ) then
          call die_msg( 229340252,                                            &
                        "Size mismatch for profile layer densities for "//    &
                        "profile '"//trim( profile%handle_%to_char( ) )//     &
                        "'. Expected "//trim( to_char( size_profile ) )//     &
                        ", got "//trim( to_char( size_host ) ) )
        end if
        profile%layer_dens_(:) = layer_densities(:)

        ! exo densities
        profile%exo_layer_dens_( 1 : profile%ncells_ ) = profile%layer_dens_(:)
        if( present( exo_density ) ) then
          profile%exo_layer_dens_( profile%ncells_ + 1 ) = exo_density
        else
          profile%exo_layer_dens_( profile%ncells_ + 1 ) =                    &
              profile%edge_val_( profile%ncells_ + 1 ) * profile%hscale_      &
              * 1.0e5_dk ! km to cm conversion
        end if
        profile%layer_dens_( profile%ncells_ ) =                              &
            profile%layer_dens_( profile%ncells_ ) +                          &
            profile%exo_layer_dens_( profile%ncells_ + 1 )

        ! burden densities
        accum = 0.0_dk
        do i_elem = profile%ncells_, 1, -1
          accum = accum + profile%layer_dens_( i_elem )
          profile%burden_dens_( i_elem ) = accum
        end do
      end if
    end associate

  end subroutine update

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_profile_from_host
