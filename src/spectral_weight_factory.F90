! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0

module tuvx_spectral_weight_factory
  ! Builds spectral weight calculators

  use tuvx_spectral_weight,            only : spectral_weight_t
  use tuvx_spectral_weight_notch_filter,                                      &
    only : spectral_weight_notch_filter_t
  use tuvx_spectral_weight_gaussian,   only : spectral_weight_gaussian_t
  use tuvx_spectral_weight_eppley,     only : spectral_weight_eppley_t
  use tuvx_spectral_weight_par,        only : spectral_weight_par_t
  use tuvx_spectral_weight_exp_decay,  only : spectral_weight_exp_decay_t
  use tuvx_spectral_weight_scup_mice,  only : spectral_weight_scup_mice_t
  use tuvx_spectral_weight_standard_human_erythema,                           &
    only : spectral_weight_standard_human_erythema_t
  use tuvx_spectral_weight_uv_index,   only : spectral_weight_uv_index_t
  use tuvx_spectral_weight_phytoplankton_boucher,                             &
    only : spectral_weight_phytoplankton_boucher_t
  use tuvx_spectral_weight_plant_damage,                                      &
    only : spectral_weight_plant_damage_t
  use tuvx_spectral_weight_plant_damage_flint_caldwell,                       &
    only : spectral_weight_plant_damage_flint_caldwell_t
  use tuvx_spectral_weight_plant_damage_flint_caldwell_ext,                   &
    only : spectral_weight_plant_damage_flint_caldwell_ext_t

  implicit none

  private
  public :: spectral_weight_builder, spectral_weight_type_name,               &
            spectral_weight_allocate

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function spectral_weight_builder( config, grid_warehouse,                   &
      profile_warehouse ) result( new_spectral_weight )

    use musica_assert,                 only : die_msg
    use musica_config,                 only : config_t
    use musica_string,                 only : string_t
    use tuvx_grid_warehouse,           only : grid_warehouse_t
    use tuvx_profile_warehouse,        only : profile_warehouse_t

    class(spectral_weight_t),  pointer       :: new_spectral_weight   ! This :f:type:`~tuvx_spectral_weight/spectral_weight_t`
    type(config_t),            intent(inout) :: config ! Spectral weight configuration data
    type(grid_warehouse_t),    intent(inout) :: grid_warehouse ! A :f:type:`~tuvx_grid_warehouse/grid_warehouse_t`
    type(profile_warehouse_t), intent(inout) :: profile_warehouse ! A :f:type:`~tuvx_profile_warehouse/profile_warehouse_t`

    type(string_t) :: spectral_weight_type
    character(len=*), parameter :: Iam = 'spectral weight builder: '

    new_spectral_weight => null()
    call config%get( 'type', spectral_weight_type, Iam )

    select case( spectral_weight_type%to_char() )
      case( 'base' )
        new_spectral_weight => spectral_weight_t( config, grid_warehouse,     &
                                                  profile_warehouse )
      case( 'Notch Filter' )
        new_spectral_weight => spectral_weight_notch_filter_t( config,        &
                                                           grid_warehouse,    &
                                                           profile_warehouse )
      case( 'Gaussian' )
        new_spectral_weight => spectral_weight_gaussian_t( config,            &
                                                           grid_warehouse,    &
                                                           profile_warehouse )
      case( 'Eppley UV Photometer' )
        new_spectral_weight => spectral_weight_eppley_t( config,              &
                                                         grid_warehouse,      &
                                                         profile_warehouse )
      case( 'PAR, 400-700 nm, umol m-2 s-1' )
        new_spectral_weight => spectral_weight_par_t( config,                 &
                                                      grid_warehouse,         &
                                                      profile_warehouse )
      case( 'Exponential decay, 14 nm/10' )
        new_spectral_weight => spectral_weight_exp_decay_t( config,           &
                                                            grid_warehouse,   &
                                                            profile_warehouse )
      case( 'SCUP-mice (de Gruijl et al., 1993)' )
        new_spectral_weight => spectral_weight_scup_mice_t( config,           &
                                                            grid_warehouse,   &
                                                            profile_warehouse )
      case( 'Standard human erythema (Webb et al., 2011)' )
        new_spectral_weight =>                                                &
            spectral_weight_standard_human_erythema_t( config, grid_warehouse,&
                                                       profile_warehouse )
      case( 'UV index (WMO, 1994; Webb et al., 2011)' )
        new_spectral_weight => spectral_weight_uv_index_t( config,            &
                                                           grid_warehouse,    &
                                                           profile_warehouse )
      case( 'Phytoplankton (Boucher et al., 1994)' )
        new_spectral_weight =>                                                &
            spectral_weight_phytoplankton_boucher_t( config,                  &
                                                     grid_warehouse,          &
                                                     profile_warehouse )
      case( 'Plant damage (Caldwell, 1971)' )
        new_spectral_weight => spectral_weight_plant_damage_t( config,        &
                                                            grid_warehouse,   &
                                                            profile_warehouse )
      case( 'Plant damage,Flint&Caldwell,2003,orig.' )
        new_spectral_weight =>                                                &
            spectral_weight_plant_damage_flint_caldwell_t( config,            &
                                                           grid_warehouse,    &
                                                           profile_warehouse )
      case( 'Plant damage,Flint&Caldwell,2003,ext390' )
        new_spectral_weight =>                                                &
            spectral_weight_plant_damage_flint_caldwell_ext_t( config,        &
                                                            grid_warehouse,   &
                                                            profile_warehouse )
      case default
        call die_msg( 450768215, "Invalid spectral weight type: '"//          &
                                 spectral_weight_type%to_char( )//"'" )
    end select

  end function spectral_weight_builder

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  type(string_t) function spectral_weight_type_name( spectral_weight )        &
      result( name )
    ! Returns the type of a spectral weight as a string

    use musica_assert,                 only : die
    use musica_string,                 only : string_t

    class(spectral_weight_t), intent(in) :: spectral_weight ! spectral weight to return type for

    select type( spectral_weight )
      type is( spectral_weight_t )
        name = "spectral_weight_t"
      type is( spectral_weight_notch_filter_t )
        name = "spectral_weight_notch_filter_t"
      type is( spectral_weight_gaussian_t )
        name = "spectral_weight_gaussian_t"
      type is( spectral_weight_eppley_t )
        name = "spectral_weight_eppley_t"
      type is( spectral_weight_par_t )
        name = "spectral_weight_par_t"
      type is( spectral_weight_exp_decay_t )
        name = "spectral_weight_exp_decay_t"
      type is( spectral_weight_scup_mice_t )
        name = "spectral_weight_scup_mice_t"
      type is( spectral_weight_standard_human_erythema_t )
        name = "spectral_weight_standard_human_erythema_t"
      type is( spectral_weight_uv_index_t )
        name = "spectral_weight_uv_index_t"
      type is( spectral_weight_phytoplankton_boucher_t )
        name = "spectral_weight_phytoplankton_boucher_t"
      type is( spectral_weight_plant_damage_t )
        name = "spectral_weight_plant_damage_t"
      type is( spectral_weight_plant_damage_flint_caldwell_t )
        name = "spectral_weight_plant_damage_flint_caldwell_t"
      type is( spectral_weight_plant_damage_flint_caldwell_ext_t )
        name = "spectral_weight_plant_damage_flint_caldwell_ext_t"
      class default
        call die( 312672872 )
    end select

  end function spectral_weight_type_name

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function spectral_weight_allocate( type_name ) result( spectral_weight )
    ! Allocates a spectral weight as the base class or a subclass by type

    use musica_assert,                 only : die_msg
    use musica_string,                 only : string_t

    type(string_t),           intent(in) :: type_name       ! Name of the type to allocate
    class(spectral_weight_t), pointer    :: spectral_weight ! Allocated spectral weight

    select case( type_name%to_char( ) )
      case( 'spectral_weight_t' )
        allocate( spectral_weight_t :: spectral_weight )
      case( 'spectral_weight_notch_filter_t' )
        allocate( spectral_weight_notch_filter_t :: spectral_weight )
      case( 'spectral_weight_gaussian_t' )
        allocate( spectral_weight_gaussian_t :: spectral_weight )
      case( 'spectral_weight_eppley_t' )
        allocate( spectral_weight_eppley_t :: spectral_weight )
      case( 'spectral_weight_par_t' )
        allocate( spectral_weight_par_t :: spectral_weight )
      case( 'spectral_weight_exp_decay_t' )
        allocate( spectral_weight_exp_decay_t :: spectral_weight )
      case( 'spectral_weight_scup_mice_t' )
        allocate( spectral_weight_scup_mice_t :: spectral_weight )
      case( 'spectral_weight_standard_human_erythema_t' )
        allocate( spectral_weight_standard_human_erythema_t ::                &
                  spectral_weight )
      case( 'spectral_weight_uv_index_t' )
        allocate( spectral_weight_uv_index_t :: spectral_weight )
      case( 'spectral_weight_phytoplankton_boucher_t' )
        allocate( spectral_weight_phytoplankton_boucher_t :: spectral_weight )
      case( 'spectral_weight_plant_damage_t' )
        allocate( spectral_weight_plant_damage_t :: spectral_weight )
      case( 'spectral_weight_plant_damage_flint_caldwell_t' )
        allocate( spectral_weight_plant_damage_flint_caldwell_t ::            &
                  spectral_weight )
      case( 'spectral_weight_plant_damage_flint_caldwell_ext_t' )
        allocate( spectral_weight_plant_damage_flint_caldwell_ext_t ::        &
                  spectral_weight )
      case default
        call die_msg( 299237780, "Invalid spectral weight type: '"//          &
                                 type_name//"'" )
    end select

  end function spectral_weight_allocate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_spectral_weight_factory
