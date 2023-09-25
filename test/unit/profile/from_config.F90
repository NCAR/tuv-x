! Copyright (C) 2023 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> Tests the photo_decomp Profile module

!> Test module for the Profile_t type
program test_profile_from_config

  use musica_string,                   only : string_t

  implicit none

  type(string_t) :: config_file

  config_file = "test/data/profile.from_config.json"
  call test_profile_from_config_t( config_file )

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Tests the profile_from_config_t type
  subroutine test_profile_from_config_t( config_file_path )

    use musica_config,                 only : config_t
    use musica_assert,                 only : assert, almost_equal
    use musica_constants,              only : dk => musica_dk
    use tuvx_grid,                     only : grid_t
    use tuvx_grid_warehouse,           only : grid_warehouse_t
    use tuvx_profile,                  only : profile_t
    use tuvx_profile_warehouse,        only : profile_warehouse_t
    use tuvx_test_utils,               only : check_values

    type(string_t), intent(in) :: config_file_path

    character(len=*), parameter :: Iam = "test profile_from_config_t"
    type(config_t)              :: config, child_config
    type(grid_warehouse_t),    pointer :: grids
    type(profile_warehouse_t), pointer :: profiles
    class(profile_t),          pointer :: profile

    call config%from_file( config_file_path%to_char( ) )

    call config%get( "grids", child_config, Iam )
    grids => grid_warehouse_t( child_config )

    call config%get( "profiles", child_config, Iam )
    profiles => profile_warehouse_t( child_config, grids )

    profile => profiles%get_profile( "by height", "mol cm-3" )

    call assert( 183062109, profile%hscale_ .eq. 0.0_dk )
    call check_values( 296546570, profile%edge_val_,                          &
                      [ 0.0_dk, 10.0_dk, 20.0_dk, 30.0_dk, 40.0_dk, 50.0_dk ],&
                      1.0e-8_dk )
    call check_values( 744653834, profile%mid_val_,                           &
                      [ 5.0_dk, 15.0_dk, 25.0_dk, 35.0_dk, 45.0_dk ],         &
                      1.0e-8_dk )
    call check_values( 678192195, profile%delta_val_,                         &
                      [ 10.0_dk, 10.0_dk, 10.0_dk, 10.0_dk, 10.0_dk ],        &
                      1.0e-8_dk )
    call check_values( 283850896, profile%layer_dens_,                        &
                     [ 5.0e5_dk, 15.0e5_dk, 25.0e5_dk, 35.0e5_dk, 45.0e5_dk ],&
                     1.0e-8_dk )
    call check_values( 209218472, profile%exo_layer_dens_(1:5),               &
                       profile%layer_dens_, 1.0e-8_dk )
    call assert( 990182580, almost_equal( profile%exo_layer_dens_(6),         &
                 0.0_dk ) )
    call assert( 216223141, size( profile%exo_layer_dens_ ) .eq. 6 )

    deallocate( profile )

    profile => profiles%get_profile( "by not height", "foos" )

    call assert( 808043823, profile%hscale_ .eq. 10.0_dk )
    call check_values( 697631012, profile%edge_val_,                          &
                       [ 12.5_dk, 12.5_dk, 12.5_dk ], 1.0e-8_dk )
    call check_values( 971483181, profile%mid_val_,                           &
                       [ 12.5_dk, 12.5_dk ], 1.0e-8_dk )
    call check_values( 120698661, profile%delta_val_,                         &
                       [ 0.0_dk, 0.0_dk ], 1.0e-8_dk )
    call check_values( 117339888, profile%layer_dens_,                        &
                       [ 25.0_dk, 25.0_dk + 125.0_dk ], 1.0e-8_dk )
    call check_values( 343882112, profile%exo_layer_dens_,                    &
                       [ 25.0_dk, 25.0_dk, 125.0_dk ], 1.0e-8_dk )

    deallocate( profile )

    profile => profiles%get_profile( "single value", "foos" )

    call assert( 711376262, profile%ncells_ .eq. 0 )
    call check_values( 529237505, profile%edge_val_, [ 12.0_dk ], 1.0e-8_dk )
    call assert( 990492738, size( profile%mid_val_        ) .eq. 0 )
    call assert( 468040173, size( profile%delta_val_      ) .eq. 0 )
    call assert( 527784266, size( profile%layer_dens_     ) .eq. 0 )
    call check_values( 305053110, profile%exo_layer_dens_, [ 0.0_dk ],        &
                       1.0e-8_dk )
    call assert( 752420956, size( profile%burden_dens_    ) .eq. 0 )

    deallocate( profile )

    deallocate( grids )
    deallocate( profiles )

  end subroutine test_profile_from_config_t

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program test_profile_from_config
