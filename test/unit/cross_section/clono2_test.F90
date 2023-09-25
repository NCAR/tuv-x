! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> Tests for the base cross_section_t type
program test_cross_section

  use tuvx_cross_section, only : cross_section_t
  use tuvx_cross_section_clono2
  use tuvx_test_utils, only : check_values

  implicit none

  call test_cross_section_clono2_t( )

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine test_cross_section_clono2_t( )

    use musica_assert,                 only : assert
    use musica_constants,              only : dk => musica_dk
    use musica_config,                 only : config_t
    use musica_iterator,               only : iterator_t
    use tuvx_grid_warehouse,           only : grid_warehouse_t
    use tuvx_profile_warehouse,        only : profile_warehouse_t

    class(grid_warehouse_t),    pointer :: grids
    class(profile_warehouse_t), pointer :: profiles
    class(cross_section_t),     pointer :: cross_section

    character(len=*), parameter :: Iam = "clono2 cross section tests"
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
      1.51e-21, 1.51e-21, 1.51e-21, 1.51e-21, 1.51e-21, 1.51e-21,             &
      3.01e-18, 3.01e-18, 3.00e-18, 2.99e-18, 2.98e-18, 2.96e-18,             &
      2.87e-18, 2.87e-18, 2.86e-18, 2.86e-18, 2.85e-18, 2.83e-18,             &
      2.79e-18, 2.79e-18, 2.78e-18, 2.78e-18, 2.77e-18, 2.76e-18,             &
      2.78e-18, 2.78e-18, 2.78e-18, 2.78e-18, 2.77e-18, 2.77e-18,             &
      2.84e-18, 2.84e-18, 2.85e-18, 2.84e-18, 2.84e-18, 2.84e-18,             &
      2.94e-18, 2.95e-18, 2.95e-18, 2.96e-18, 2.96e-18, 2.95e-18,             &
      3.08e-18, 3.09e-18, 3.09e-18, 3.10e-18, 3.10e-18, 3.10e-18,             &
      3.23e-18, 3.24e-18, 3.25e-18, 3.26e-18, 3.26e-18, 3.26e-18,             &
      3.36e-18, 3.37e-18, 3.38e-18, 3.39e-18, 3.40e-18, 3.40e-18,             &
      1.81e-21, 1.81e-21, 1.81e-21, 1.81e-21, 1.81e-21, 1.81e-21,             &
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0,                                           &
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0],                                          &
      (/ size(no_extrap, 2), size(no_extrap, 1) /)                            &
    )

    lower_extrap = reshape([                                                  &
      1900723.24, 6770147.44, 14622344.20,                                    &
      25454557.94, 39264033.97, 56045203.29,                                  &
      1900723.24, 6770147.44, 14622344.20,                                    &
      25454557.94, 39264033.97, 56045203.29,                                  &
      1900723.24, 6770147.44, 14622344.20,                                    &
      25454557.94, 39264033.97, 56045203.29,                                  &
      1900723.24, 6770147.44, 14622344.20,                                    &
      25454557.94, 39264033.97, 56045203.29,                                  &
      1900723.24, 6770147.44, 14622344.20,                                    &
      25454557.94, 39264033.97, 56045203.29,                                  &
      1900723.24, 6770147.44, 14622344.20,                                    &
      25454557.94, 39264033.97, 56045203.29,                                  &
      1900723.24, 6770147.44, 14622344.20,                                    &
      25454557.94, 39264033.97, 56045203.29,                                  &
      1898861.07, 6763514.41, 14608017.90,                                    &
      25429618.68, 39225564.73, 55990292.54,                                  &
      3.01e-18, 3.01e-18, 3.00e-18, 2.99e-18, 2.98e-18, 2.96e-18,             &
      2.87e-18, 2.87e-18, 2.86e-18, 2.86e-18, 2.85e-18, 2.83e-18,             &
      2.79e-18, 2.79e-18, 2.78e-18, 2.78e-18, 2.77e-18, 2.76e-18,             &
      2.78e-18, 2.78e-18, 2.78e-18, 2.78e-18, 2.77e-18, 2.77e-18,             &
      2.84e-18, 2.84e-18, 2.85e-18, 2.84e-18, 2.84e-18, 2.84e-18,             &
      2.94e-18, 2.95e-18, 2.95e-18, 2.96e-18, 2.96e-18, 2.95e-18,             &
      3.08e-18, 3.09e-18, 3.09e-18, 3.10e-18, 3.10e-18, 3.10e-18,             &
      3.23e-18, 3.24e-18, 3.25e-18, 3.26e-18, 3.26e-18, 3.26e-18,             &
      3.36e-18, 3.37e-18, 3.38e-18, 3.39e-18, 3.40e-18, 3.40e-18,             &
      1.81e-21, 1.81e-21, 1.81e-21, 1.81e-21, 1.81e-21, 1.81e-21,             &
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0,                                           &
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0],                                          &
      (/ size(lower_extrap, 2), size(lower_extrap, 1) /)                      &
    )

    upper_extrap = reshape([                                                  &
      3.09e-18, 3.09e-18, 3.08e-18, 3.07e-18, 3.05e-18, 3.04e-18,             &
      3.09e-18, 3.09e-18, 3.08e-18, 3.07e-18, 3.05e-18, 3.04e-18,             &
      3.09e-18, 3.09e-18, 3.08e-18, 3.07e-18, 3.05e-18, 3.04e-18,             &
      3.09e-18, 3.09e-18, 3.08e-18, 3.07e-18, 3.05e-18, 3.04e-18,             &
      3.09e-18, 3.09e-18, 3.08e-18, 3.07e-18, 3.05e-18, 3.04e-18,             &
      3.09e-18, 3.09e-18, 3.08e-18, 3.07e-18, 3.05e-18, 3.04e-18,             &
      3.09e-18, 3.09e-18, 3.08e-18, 3.07e-18, 3.05e-18, 3.04e-18,             &
      3.09e-18, 3.09e-18, 3.08e-18, 3.07e-18, 3.05e-18, 3.04e-18,             &
      3.01e-18, 3.01e-18, 3.00e-18, 2.99e-18, 2.98e-18, 2.96e-18,             &
      2.87e-18, 2.87e-18, 2.86e-18, 2.86e-18, 2.85e-18, 2.83e-18,             &
      2.79e-18, 2.79e-18, 2.78e-18, 2.78e-18, 2.77e-18, 2.76e-18,             &
      2.78e-18, 2.78e-18, 2.78e-18, 2.78e-18, 2.77e-18, 2.77e-18,             &
      2.84e-18, 2.84e-18, 2.85e-18, 2.84e-18, 2.84e-18, 2.84e-18,             &
      2.94e-18, 2.95e-18, 2.95e-18, 2.96e-18, 2.96e-18, 2.95e-18,             &
      3.08e-18, 3.09e-18, 3.09e-18, 3.10e-18, 3.10e-18, 3.10e-18,             &
      3.23e-18, 3.24e-18, 3.25e-18, 3.26e-18, 3.26e-18, 3.26e-18,             &
      3.36e-18, 3.37e-18, 3.38e-18, 3.39e-18, 3.40e-18, 3.40e-18,             &
      2414377.89, 8599792.48, 18574091.35, 32333774.23, 49875341.95, 71191720.46,&
      2416963.24, 8609001.53, 18593981.47, 32368399.04, 49928751.30, 71267956.66,&
      2416963.24, 8609001.53, 18593981.47, 32368399.04, 49928751.30, 71267956.66],&
      (/ size(upper_extrap, 2), size(upper_extrap, 1) /)                      &
    )

    ! load test grids
    call config%from_file( "test/data/grid.180-220.config.json" )
    grids => grid_warehouse_t( config )

    ! load test profiles
    call config%from_file( "test/data/profile.temperature.config.json" )
    profiles => profile_warehouse_t( config, grids )

    ! get cross section config data
    call config%from_file( "test/data/cross_sections/cross_section.clono2.config.json" )
    call config%get( "cross sections", cs_set, Iam )
    iter => cs_set%get_iterator( )

    ! load and test cross section w/o extrapolation
    call assert( 101264914, iter%next( ) )
    call cs_set%get( iter, cs_config, Iam )
    cross_section => cross_section_clono2_t( cs_config, grids, profiles )
    results = cross_section%calculate( grids, profiles )
    call check_values( results, no_extrap, .01_dk )
    deallocate( cross_section )

    ! load and test cross section with lower extrapolation
    call assert( 101264915, iter%next( ) )
    call cs_set%get( iter, cs_config, Iam )
    cross_section => cross_section_clono2_t( cs_config, grids, profiles )
    results = cross_section%calculate( grids, profiles )
    call check_values( results, lower_extrap, .01_dk )
    deallocate( cross_section )

    ! load and test cross section with lower and upper extrapolation
    call assert( 101264916, iter%next( ) )
    call cs_set%get( iter, cs_config, Iam )
    cross_section => cross_section_clono2_t( cs_config, grids, profiles )
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

  end subroutine test_cross_section_clono2_t

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program test_cross_section
