! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> Tests for the base cross_section_t type
program test_cross_section

  use tuvx_cross_section, only : cross_section_t
  use tuvx_cross_section_ccl4
  use tuvx_test_utils, only : check_values

  implicit none

  call test_cross_section_ccl4_t( )

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine test_cross_section_ccl4_t( )

    use musica_assert,                 only : assert
    use musica_constants,              only : dk => musica_dk
    use musica_config,                 only : config_t
    use musica_iterator,               only : iterator_t
    use tuvx_grid_warehouse,           only : grid_warehouse_t
    use tuvx_profile_warehouse,        only : profile_warehouse_t

    class(grid_warehouse_t),    pointer :: grids
    class(profile_warehouse_t), pointer :: profiles
    class(cross_section_t),     pointer :: cross_section

    character(len=*), parameter :: Iam = "ccl4 cross section tests"
    type(config_t) :: config, cs_set, cs_config
    class(iterator_t), pointer :: iter
    real(kind=dk), allocatable :: results(:,:)
    real(dk), allocatable :: no_extrap(:,:)
    real(dk), allocatable :: lower_extrap(:,:)
    real(dk), allocatable :: upper_extrap(:,:)
    allocate(no_extrap(20, 6))
    allocate(lower_extrap(20, 6))
    allocate(upper_extrap(20, 6))

    ! All of these values were produced by one run of the cross section.
    ! So, these tests are testing that any changes don't produce unexpected
    ! changes. The values here are meaningless.
    no_extrap = reshape([                                                     &
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0,                                           &
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0,                                           &
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0,                                           &
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0,                                           &
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0,                                           &
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0,                                           &
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0,                                           &
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0,                                           &
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0,                                           &
      3.29e-22, 3.29e-22, 3.29e-22, 3.29e-22, 3.29e-22, 3.29e-22,             &
      6.48e-19, 6.48e-19, 6.48e-19, 6.47e-19, 6.47e-19, 6.47e-19,             &
      6.22e-19, 6.22e-19, 6.21e-19, 6.20e-19, 6.19e-19, 6.18e-19,             &
      5.88e-19, 5.86e-19, 5.84e-19, 5.82e-19, 5.80e-19, 5.78e-19,             &
      5.44e-19, 5.41e-19, 5.38e-19, 5.35e-19, 5.33e-19, 5.30e-19,             &
      4.92e-19, 4.89e-19, 4.85e-19, 4.81e-19, 4.77e-19, 4.73e-19,             &
      2.43e-22, 2.40e-22, 2.38e-22, 2.35e-22, 2.33e-22, 2.30e-22,             &
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0,                                           &
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0,                                           &
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0,                                           &
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0],                                          &
      (/ size(no_extrap, 2), size(no_extrap, 1) /)                            &
    )

    lower_extrap = reshape([                                                  &
      180.0, 180.0, 180.0, 180.0, 180.0, 180.0,                               &
      180.0, 180.0, 180.0, 180.0, 180.0, 180.0,                               &
      180.0, 180.0, 180.0, 180.0, 180.0, 180.0,                               &
      180.0, 180.0, 180.0, 180.0, 180.0, 180.0,                               &
      180.0, 180.0, 180.0, 180.0, 180.0, 180.0,                               &
      180.0, 180.0, 180.0, 180.0, 180.0, 180.0,                               &
      180.0, 180.0, 180.0, 180.0, 180.0, 180.0,                               &
      179.77, 179.56, 179.35, 179.14, 178.93, 178.71,                         &
      179.95, 179.90, 179.86, 179.82, 179.77, 179.73,                         &
      179.90, 179.89, 179.89, 179.88, 179.88, 179.87,                         &
      6.48e-19, 6.48e-19, 6.48e-19, 6.47e-19, 6.47e-19, 6.47e-19,             &
      6.22e-19, 6.22e-19, 6.21e-19, 6.20e-19, 6.19e-19, 6.18e-19,             &
      5.88e-19, 5.86e-19, 5.84e-19, 5.82e-19, 5.80e-19, 5.78e-19,             &
      5.44e-19, 5.41e-19, 5.38e-19, 5.35e-19, 5.33e-19, 5.30e-19,             &
      4.92e-19, 4.89e-19, 4.85e-19, 4.81e-19, 4.77e-19, 4.73e-19,             &
      2.43e-22, 2.40e-22, 2.38e-22, 2.35e-22, 2.33e-22, 2.30e-22,             &
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0,                                           &
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0,                                           &
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0,                                           &
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0],                                          &
      (/ size(lower_extrap, 2), size(lower_extrap, 1) /)                      &
    )

    upper_extrap = reshape([                                                  &
      6.59e-19, 6.59e-19, 6.59e-19, 6.59e-19, 6.59e-19, 6.59e-19,             &
      6.59e-19, 6.59e-19, 6.59e-19, 6.59e-19, 6.59e-19, 6.59e-19,             &
      6.59e-19, 6.59e-19, 6.59e-19, 6.59e-19, 6.59e-19, 6.59e-19,             &
      6.59e-19, 6.59e-19, 6.59e-19, 6.59e-19, 6.59e-19, 6.59e-19,             &
      6.59e-19, 6.59e-19, 6.59e-19, 6.59e-19, 6.59e-19, 6.59e-19,             &
      6.59e-19, 6.59e-19, 6.59e-19, 6.59e-19, 6.59e-19, 6.59e-19,             &
      6.59e-19, 6.59e-19, 6.59e-19, 6.59e-19, 6.59e-19, 6.59e-19,             &
      6.59e-19, 6.58e-19, 6.57e-19, 6.56e-19, 6.56e-19, 6.55e-19,             &
      6.59e-19, 6.59e-19, 6.59e-19, 6.59e-19, 6.59e-19, 6.59e-19,             &
      6.59e-19, 6.59e-19, 6.59e-19, 6.59e-19, 6.59e-19, 6.59e-19,             &
      6.48e-19, 6.48e-19, 6.48e-19, 6.47e-19, 6.47e-19, 6.47e-19,             &
      6.22e-19, 6.22e-19, 6.21e-19, 6.20e-19, 6.19e-19, 6.18e-19,             &
      5.88e-19, 5.86e-19, 5.84e-19, 5.82e-19, 5.80e-19, 5.78e-19,             &
      5.44e-19, 5.41e-19, 5.38e-19, 5.35e-19, 5.33e-19, 5.30e-19,             &
      4.92e-19, 4.89e-19, 4.85e-19, 4.81e-19, 4.77e-19, 4.73e-19,             &
      217.39, 215.05, 212.74, 210.45, 208.19, 205.95,                         &
      216.74, 213.71, 210.71, 207.76, 204.85, 201.98,                         &
      215.92, 212.13, 208.41, 204.75, 201.16, 197.63,                         &
      215.05, 210.47, 205.98, 201.59, 197.29, 193.09,                         &
      214.14, 208.73, 203.46, 198.32, 193.32, 188.44],                        &
      (/ size(upper_extrap, 2), size(upper_extrap, 1) /)                      &
    )


    ! load test grids
    call config%from_file( "test/data/grid.180-220.config.json" )
    grids => grid_warehouse_t( config )

    ! load test profiles
    call config%from_file( "test/data/profile.temperature.config.json" )
    profiles => profile_warehouse_t( config, grids )

    ! get cross section config data
    call config%from_file( "test/data/cross_sections/cross_section.ccl4.config.json" )
    call config%get( "cross sections", cs_set, Iam )
    iter => cs_set%get_iterator( )

    ! load and test cross section w/o extrapolation
    call assert( 101264914, iter%next( ) )
    call cs_set%get( iter, cs_config, Iam )
    cross_section => cross_section_ccl4_t( cs_config, grids, profiles )
    results = cross_section%calculate( grids, profiles )
    call check_values( results, no_extrap, .01_dk )
    deallocate( cross_section )

    ! load and test cross section with lower extrapolation
    call assert( 101264914, iter%next( ) )
    call cs_set%get( iter, cs_config, Iam )
    cross_section => cross_section_ccl4_t( cs_config, grids, profiles )
    results = cross_section%calculate( grids, profiles )
    call check_values( results, lower_extrap, .01_dk )
    deallocate( cross_section )

    ! load and test cross section with lower and upper extrapolation
    call assert( 101264914, iter%next( ) )
    call cs_set%get( iter, cs_config, Iam )
    cross_section => cross_section_ccl4_t( cs_config, grids, profiles )
    results = cross_section%calculate( grids, profiles )
    call check_values( results, upper_extrap, .01_dk )
    deallocate( cross_section )

    ! clean up
    deallocate( iter )
    deallocate( grids )
    deallocate( profiles )
    deallocate( no_extrap )
    deallocate( lower_extrap )
    deallocate( upper_extrap )

  end subroutine test_cross_section_ccl4_t

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program test_cross_section
