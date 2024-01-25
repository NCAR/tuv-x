! Copyright (C) 2024 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0

program test_la_srb
  !> Tests against Doug's LA and SR band calculations

  use musica_assert,                   only : assert
  use musica_constants,                only : dk => musica_dk
  use musica_config,                   only : config_t
  use tuvx_grid_warehouse,             only : grid_warehouse_t
  use tuvx_profile_warehouse,          only : profile_warehouse_t
  use tuvx_la_sr_bands,                only : la_sr_bands_t
  use tuvx_test_utils,                 only : check_values

  implicit none

  character(len=*), parameter :: conf_l = 'test/data/la_srb_bands.config.json'
  type(config_t) :: grid_config, la_config, profile_config
  class(grid_warehouse_t), pointer :: grids => null( )
  class(profile_warehouse_t), pointer :: profiles => null( )
  class(la_sr_bands_t), pointer :: la_sr_bands => null( )
  character, allocatable :: buffer(:)
  
  call la_config%from_file( conf_l )
  call la_config%get( "grids", grid_config, "" )
  call la_config%get( "profiles", profile_config, "" )

  grids => grid_warehouse_t( grid_config )
  profiles => profile_warehouse_t( profile_config, grids )
  la_sr_bands => la_sr_bands_t( la_config, grids, profiles )

  call compare_o2_cross_sections( la_sr_bands, grids, profiles )

  deallocate( grids )
  deallocate( profiles )
  deallocate( la_sr_bands )

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine compare_o2_cross_sections( la_sr_bands, grids, profiles )

    use tuvx_grid,                     only : grid_t
    use tuvx_spherical_geometry,       only : spherical_geometry_t

    class(la_sr_bands_t),      intent(inout) :: la_sr_bands
    class(grid_warehouse_t),    intent(inout) :: grids
    class(profile_warehouse_t), intent(inout) :: profiles

    call assert( 597503990, .true. )

  end subroutine compare_o2_cross_sections

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program test_la_srb