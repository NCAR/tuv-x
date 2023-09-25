! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> Tests for the base cross_section_t type
program test_cross_section

  use tuvx_cross_section, only : cross_section_t
  use tuvx_cross_section_chcl3
  use tuvx_test_utils, only : check_values

  implicit none

  call test_cross_section_chcl3_t( )

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine test_cross_section_chcl3_t( )

    use musica_assert,                 only : assert
    use musica_constants,              only : dk => musica_dk
    use musica_config,                 only : config_t
    use musica_iterator,               only : iterator_t
    use tuvx_grid_warehouse,           only : grid_warehouse_t
    use tuvx_profile_warehouse,        only : profile_warehouse_t

    class(grid_warehouse_t),    pointer :: grids
    class(profile_warehouse_t), pointer :: profiles
    class(cross_section_t),     pointer :: cross_section

    character(len=*), parameter :: Iam = "chcl3 cross section tests"
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
      5.36e-22, 5.36e-22, 5.36e-22, 5.36e-22, 5.36e-22, 5.36e-22,             &
      1.00e-18, 9.99e-19, 9.92e-19, 9.84e-19, 9.77e-19, 9.70e-19,             &
      8.21e-19, 8.13e-19, 8.05e-19, 7.97e-19, 7.89e-19, 7.82e-19,             &
      6.92e-19, 6.83e-19, 6.74e-19, 6.65e-19, 6.57e-19, 6.48e-19,             &
      5.75e-19, 5.66e-19, 5.57e-19, 5.48e-19, 5.39e-19, 5.30e-19,             &
      4.67e-19, 4.58e-19, 4.49e-19, 4.40e-19, 4.31e-19, 4.22e-19,             &
      3.75e-19, 3.66e-19, 3.57e-19, 3.48e-19, 3.40e-19, 3.31e-19,             &
      2.98e-19, 2.90e-19, 2.82e-19, 2.73e-19, 2.66e-19, 2.58e-19,             &
      2.31e-19, 2.23e-19, 2.16e-19, 2.09e-19, 2.02e-19, 1.95e-19,             &
      1.72e-19, 1.65e-19, 1.59e-19, 1.53e-19, 1.47e-19, 1.42e-19,             &
      1.23e-19, 1.18e-19, 1.13e-19, 1.08e-19, 1.04e-19, 9.99e-20,             &
      5.34e-23, 5.10e-23, 4.87e-23, 4.65e-23, 4.44e-23, 4.24e-23,             &
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0,                                           &
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0,                                           &
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0,                                           &
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0],                                          &
      (/ size(no_extrap, 2), size(no_extrap, 1) /)                            &
    )

    lower_extrap = reshape([                                                  &
      188.0, 188.0, 188.0, 188.0, 188.0, 188.0,                               &
      188.0, 188.0, 188.0, 188.0, 188.0, 188.0,                               &
      188.0, 188.0, 188.0, 188.0, 188.0, 188.0,                               &
      188.0, 188.0, 188.0, 188.0, 188.0, 188.0,                               &
      187.91, 187.91, 187.91, 187.91, 187.91, 187.91,                         &
      1.00e-18, 9.99e-19, 9.92e-19, 9.84e-19, 9.77e-19, 9.70e-19,             &
      8.21e-19, 8.13e-19, 8.05e-19, 7.97e-19, 7.89e-19, 7.82e-19,             &
      6.92e-19, 6.83e-19, 6.74e-19, 6.65e-19, 6.57e-19, 6.48e-19,             &
      5.75e-19, 5.66e-19, 5.57e-19, 5.48e-19, 5.39e-19, 5.30e-19,             &
      4.67e-19, 4.58e-19, 4.49e-19, 4.40e-19, 4.31e-19, 4.22e-19,             &
      3.75e-19, 3.66e-19, 3.57e-19, 3.48e-19, 3.40e-19, 3.31e-19,             &
      2.98e-19, 2.90e-19, 2.82e-19, 2.73e-19, 2.66e-19, 2.58e-19,             &
      2.31e-19, 2.23e-19, 2.16e-19, 2.09e-19, 2.02e-19, 1.95e-19,             &
      1.72e-19, 1.65e-19, 1.59e-19, 1.53e-19, 1.47e-19, 1.42e-19,             &
      1.23e-19, 1.18e-19, 1.13e-19, 1.08e-19, 1.04e-19, 9.99e-20,             &
      5.34e-23, 5.10e-23, 4.87e-23, 4.65e-23, 4.44e-23, 4.24e-23,             &
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0,                                           &
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0,                                           &
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0,                                           &
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0],                                          &
      (/ size(lower_extrap, 2), size(lower_extrap, 1) /)                      &
    )

    upper_extrap = reshape([                                                  &
      1.13e-18, 1.13e-18, 1.13e-18, 1.13e-18, 1.13e-18, 1.13e-18,             &
      1.13e-18, 1.13e-18, 1.13e-18, 1.13e-18, 1.13e-18, 1.13e-18,             &
      1.13e-18, 1.13e-18, 1.13e-18, 1.13e-18, 1.13e-18, 1.13e-18,             &
      1.13e-18, 1.13e-18, 1.13e-18, 1.13e-18, 1.13e-18, 1.13e-18,             &
      1.13e-18, 1.13e-18, 1.13e-18, 1.13e-18, 1.13e-18, 1.13e-18,             &
      1.00e-18, 9.99e-19, 9.92e-19, 9.84e-19, 9.77e-19, 9.70e-19,             &
      8.21e-19, 8.13e-19, 8.05e-19, 7.97e-19, 7.89e-19, 7.82e-19,             &
      6.92e-19, 6.83e-19, 6.74e-19, 6.65e-19, 6.57e-19, 6.48e-19,             &
      5.75e-19, 5.66e-19, 5.57e-19, 5.48e-19, 5.39e-19, 5.30e-19,             &
      4.67e-19, 4.58e-19, 4.49e-19, 4.40e-19, 4.31e-19, 4.22e-19,             &
      3.75e-19, 3.66e-19, 3.57e-19, 3.48e-19, 3.40e-19, 3.31e-19,             &
      2.98e-19, 2.90e-19, 2.82e-19, 2.73e-19, 2.66e-19, 2.58e-19,             &
      2.31e-19, 2.23e-19, 2.16e-19, 2.09e-19, 2.02e-19, 1.95e-19,             &
      1.72e-19, 1.65e-19, 1.59e-19, 1.53e-19, 1.47e-19, 1.42e-19,             &
      1.23e-19, 1.18e-19, 1.13e-19, 1.08e-19, 1.04e-19, 9.99e-20,             &
      201.75, 192.59, 183.84, 175.49, 167.53, 159.93,                         &
      200.96, 191.02, 181.57, 172.60, 164.07, 155.96,                         &
      200.09, 189.41, 179.31, 169.75, 160.70, 152.14,                         &
      199.26, 187.88, 177.15, 167.04, 157.51, 148.53,                         &
      198.46, 186.41, 175.10, 164.48, 154.51, 145.15],                        &
      (/ size(upper_extrap, 2), size(upper_extrap, 1) /)                      &
    )

    ! load test grids
    call config%from_file( "test/data/grid.180-220.config.json" )
    grids => grid_warehouse_t( config )

    ! load test profiles
    call config%from_file( "test/data/profile.temperature.config.json" )
    profiles => profile_warehouse_t( config, grids )

    ! get cross section config data
    call config%from_file( "test/data/cross_sections/cross_section.chcl3.config.json" )
    call config%get( "cross sections", cs_set, Iam )
    iter => cs_set%get_iterator( )

    ! load and test cross section w/o extrapolation
    call assert( 101264914, iter%next( ) )
    call cs_set%get( iter, cs_config, Iam )
    cross_section => cross_section_chcl3_t( cs_config, grids, profiles )
    results = cross_section%calculate( grids, profiles )
    call check_values( results, no_extrap, .01_dk )
    deallocate( cross_section )

    ! load and test cross section with lower extrapolation
    call assert( 101264915, iter%next( ) )
    call cs_set%get( iter, cs_config, Iam )
    cross_section => cross_section_chcl3_t( cs_config, grids, profiles )
    results = cross_section%calculate( grids, profiles )
    call check_values( results, lower_extrap, .01_dk )
    deallocate( cross_section )

    ! load and test cross section with lower and upper extrapolation
    call assert( 101264916, iter%next( ) )
    call cs_set%get( iter, cs_config, Iam )
    cross_section => cross_section_chcl3_t( cs_config, grids, profiles )
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

  end subroutine test_cross_section_chcl3_t

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program test_cross_section
