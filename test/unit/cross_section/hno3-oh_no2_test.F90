! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> Tests for the base cross_section_t type
program test_cross_section

  use tuvx_cross_section, only : cross_section_t
  use tuvx_cross_section_hno3_oh_no2
  use tuvx_test_utils, only : check_values

  implicit none

  call test_cross_section_hno3_oh_no2_t( )

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine test_cross_section_hno3_oh_no2_t( )

    use musica_assert,                 only : assert
    use musica_constants,              only : dk => musica_dk
    use musica_config,                 only : config_t
    use musica_iterator,               only : iterator_t
    use tuvx_grid_warehouse,           only : grid_warehouse_t
    use tuvx_profile_warehouse,        only : profile_warehouse_t

    class(grid_warehouse_t),    pointer :: grids
    class(profile_warehouse_t), pointer :: profiles
    class(cross_section_t),     pointer :: cross_section

    character(len=*), parameter :: Iam = "hno3-oh_no2 cross section tests"
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
      6.35e-21, 6.28e-21, 6.21e-21, 6.14e-21, 6.07e-21, 6.01e-21,             &
      1.27e-17, 1.25e-17, 1.24e-17, 1.22e-17, 1.21e-17, 1.20e-17,             &
      1.14e-17, 1.12e-17, 1.11e-17, 1.10e-17, 1.09e-17, 1.07e-17,             &
      1.00e-17, 9.89e-18, 9.78e-18, 9.67e-18, 9.57e-18, 9.46e-18,             &
      8.41e-18, 8.31e-18, 8.22e-18, 8.13e-18, 8.05e-18, 7.96e-18,             &
      6.68e-18, 6.60e-18, 6.53e-18, 6.46e-18, 6.39e-18, 6.33e-18,             &
      5.09e-18, 5.03e-18, 4.98e-18, 4.92e-18, 4.87e-18, 4.82e-18,             &
      3.81e-18, 3.76e-18, 3.72e-18, 3.68e-18, 3.64e-18, 3.60e-18,             &
      2.74e-18, 2.71e-18, 2.68e-18, 2.65e-18, 2.62e-18, 2.59e-18,             &
      1.90e-18, 1.87e-18, 1.85e-18, 1.83e-18, 1.81e-18, 1.79e-18,             &
      1.27e-18, 1.26e-18, 1.24e-18, 1.22e-18, 1.21e-18, 1.19e-18,             &
      5.35e-22, 5.28e-22, 5.21e-22, 5.15e-22, 5.08e-22, 5.02e-22,             &
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0,                                           &
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0,                                           &
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0,                                           &
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0],                                          &
      (/ size(no_extrap, 2), size(no_extrap, 1) /)                            &
    )

    lower_extrap = reshape([ real::                                           &
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0,                                           &
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0,                                           &
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0,                                           &
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0,                                           &
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0,                                           &
      1.27e-17, 1.25e-17, 1.24e-17, 1.22e-17, 1.21e-17, 1.20e-17,             &
      1.14e-17, 1.12e-17, 1.11e-17, 1.10e-17, 1.09e-17, 1.07e-17,             &
      1.00e-17, 9.89e-18, 9.78e-18, 9.67e-18, 9.57e-18, 9.46e-18,             &
      8.41e-18, 8.31e-18, 8.22e-18, 8.13e-18, 8.05e-18, 7.96e-18,             &
      6.68e-18, 6.60e-18, 6.53e-18, 6.46e-18, 6.39e-18, 6.33e-18,             &
      5.09e-18, 5.03e-18, 4.98e-18, 4.92e-18, 4.87e-18, 4.82e-18,             &
      3.81e-18, 3.76e-18, 3.72e-18, 3.68e-18, 3.64e-18, 3.60e-18,             &
      2.74e-18, 2.71e-18, 2.68e-18, 2.65e-18, 2.62e-18, 2.59e-18,             &
      1.90e-18, 1.87e-18, 1.85e-18, 1.83e-18, 1.81e-18, 1.79e-18,             &
      1.27e-18, 1.26e-18, 1.24e-18, 1.22e-18, 1.21e-18, 1.19e-18,             &
      5.35e-22, 5.28e-22, 5.21e-22, 5.15e-22, 5.08e-22, 5.02e-22,             &
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0,                                           &
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0,                                           &
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0,                                           &
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0],                                          &
      (/ size(lower_extrap, 2), size(lower_extrap, 1) /)                      &
    )

    upper_extrap = reshape([ real::                                           &
      1.33e-17, 1.32e-17, 1.30e-17, 1.29e-17, 1.27e-17, 1.26e-17,             &
      1.33e-17, 1.32e-17, 1.30e-17, 1.29e-17, 1.27e-17, 1.26e-17,             &
      1.33e-17, 1.32e-17, 1.30e-17, 1.29e-17, 1.27e-17, 1.26e-17,             &
      1.33e-17, 1.32e-17, 1.30e-17, 1.29e-17, 1.27e-17, 1.26e-17,             &
      1.33e-17, 1.32e-17, 1.30e-17, 1.29e-17, 1.27e-17, 1.26e-17,             &
      1.27e-17, 1.25e-17, 1.24e-17, 1.22e-17, 1.21e-17, 1.20e-17,             &
      1.14e-17, 1.12e-17, 1.11e-17, 1.10e-17, 1.09e-17, 1.07e-17,             &
      1.00e-17, 9.89e-18, 9.78e-18, 9.67e-18, 9.57e-18, 9.46e-18,             &
      8.41e-18, 8.31e-18, 8.22e-18, 8.13e-18, 8.05e-18, 7.96e-18,             &
      6.68e-18, 6.60e-18, 6.53e-18, 6.46e-18, 6.39e-18, 6.33e-18,             &
      5.09e-18, 5.03e-18, 4.98e-18, 4.92e-18, 4.87e-18, 4.82e-18,             &
      3.81e-18, 3.76e-18, 3.72e-18, 3.68e-18, 3.64e-18, 3.60e-18,             &
      2.74e-18, 2.71e-18, 2.68e-18, 2.65e-18, 2.62e-18, 2.59e-18,             &
      1.90e-18, 1.87e-18, 1.85e-18, 1.83e-18, 1.81e-18, 1.79e-18,             &
      1.27e-18, 1.26e-18, 1.24e-18, 1.22e-18, 1.21e-18, 1.19e-18,             &
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0,                                           &
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0,                                           &
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0,                                           &
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0,                                           &
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0],                                          &
      (/ size(upper_extrap, 2), size(upper_extrap, 1) /)                      &
    )

    ! load test grids
    call config%from_file( "test/data/grid.180-220.config.json" )
    grids => grid_warehouse_t( config )

    ! load test profiles
    call config%from_file( "test/data/profile.temperature.config.json" )
    profiles => profile_warehouse_t( config, grids )

    ! get cross section config data
    call config%from_file( "test/data/cross_sections/cross_section.hno3-oh_no2.config.json" )
    call config%get( "cross sections", cs_set, Iam )
    iter => cs_set%get_iterator( )

    ! load and test cross section w/o extrapolation
    call assert( 101264914, iter%next( ) )
    call cs_set%get( iter, cs_config, Iam )
    cross_section => cross_section_hno3_oh_no2_t( cs_config, grids, profiles )
    results = cross_section%calculate( grids, profiles )
    call check_values( results, no_extrap, .01_dk )
    deallocate( cross_section )

    ! load and test cross section with lower extrapolation
    call assert( 101264915, iter%next( ) )
    call cs_set%get( iter, cs_config, Iam )
    cross_section => cross_section_hno3_oh_no2_t( cs_config, grids, profiles )
    results = cross_section%calculate( grids, profiles )
    call check_values( results, lower_extrap, .01_dk )
    deallocate( cross_section )

    ! load and test cross section with upper extrapolation
    call assert( 101264914, iter%next( ) )
    call cs_set%get( iter, cs_config, Iam )
    cross_section => cross_section_hno3_oh_no2_t( cs_config, grids, profiles )
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

  end subroutine test_cross_section_hno3_oh_no2_t

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program test_cross_section
