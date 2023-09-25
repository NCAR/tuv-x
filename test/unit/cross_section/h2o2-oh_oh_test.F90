! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> Tests for the base cross_section_t type
program test_cross_section

  use tuvx_cross_section, only : cross_section_t
  use tuvx_cross_section_h2o2_oh_oh
  use tuvx_test_utils, only : check_values

  implicit none

  call test_cross_section_h2o2_oh_oh_t( )

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine test_cross_section_h2o2_oh_oh_t( )

    use musica_assert,                 only : assert
    use musica_constants,              only : dk => musica_dk
    use musica_config,                 only : config_t
    use musica_iterator,               only : iterator_t
    use tuvx_grid_warehouse,           only : grid_warehouse_t
    use tuvx_profile_warehouse,        only : profile_warehouse_t

    class(grid_warehouse_t),    pointer :: grids
    class(profile_warehouse_t), pointer :: profiles
    class(cross_section_t),     pointer :: cross_section

    character(len=*), parameter :: Iam = "h2o2-oh_oh cross section tests"
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
      3.19e-22, 3.19e-22, 3.19e-22, 3.19e-22, 3.19e-22, 3.19e-22,             &
      6.50e-19, 6.50e-19, 6.50e-19, 6.50e-19, 6.50e-19, 6.50e-19,             &
      6.07e-19, 6.07e-19, 6.07e-19, 6.07e-19, 6.07e-19, 6.07e-19,             &
      5.64e-19, 5.64e-19, 5.64e-19, 5.64e-19, 5.64e-19, 5.64e-19,             &
      5.28e-19, 5.28e-19, 5.28e-19, 5.28e-19, 5.28e-19, 5.28e-19,             &
      4.92e-19, 4.92e-19, 4.92e-19, 4.92e-19, 4.92e-19, 4.92e-19,             &
      4.61e-19, 4.61e-19, 4.61e-19, 4.61e-19, 4.61e-19, 4.61e-19,             &
      4.34e-19, 4.34e-19, 4.34e-19, 4.34e-19, 4.34e-19, 4.34e-19,             &
      4.08e-19, 4.08e-19, 4.08e-19, 4.08e-19, 4.08e-19, 4.08e-19,             &
      3.87e-19, 3.87e-19, 3.87e-19, 3.87e-19, 3.87e-19, 3.87e-19,             &
      3.67e-19, 3.67e-19, 3.67e-19, 3.67e-19, 3.67e-19, 3.67e-19,             &
      1.87e-22, 1.87e-22, 1.87e-22, 1.87e-22, 1.87e-22, 1.87e-22,             &
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0,                                           &
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0,                                           &
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0,                                           &
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0],                                          &
      (/ size(no_extrap, 2), size(no_extrap, 1) /)                            &
    )

    lower_extrap = reshape([ real::                                           &
      188, 188, 188, 188, 188, 188,                                           &
      188, 188, 188, 188, 188, 188,                                           &
      188, 188, 188, 188, 188, 188,                                           &
      188, 188, 188, 188, 188, 188,                                           &
      187.91, 187.91, 187.91, 187.91, 187.91, 187.91,                         &
      6.50e-19, 6.50e-19, 6.50e-19, 6.50e-19, 6.50e-19, 6.50e-19,             &
      6.07e-19, 6.07e-19, 6.07e-19, 6.07e-19, 6.07e-19, 6.07e-19,             &
      5.64e-19, 5.64e-19, 5.64e-19, 5.64e-19, 5.64e-19, 5.64e-19,             &
      5.28e-19, 5.28e-19, 5.28e-19, 5.28e-19, 5.28e-19, 5.28e-19,             &
      4.92e-19, 4.92e-19, 4.92e-19, 4.92e-19, 4.92e-19, 4.92e-19,             &
      4.61e-19, 4.61e-19, 4.61e-19, 4.61e-19, 4.61e-19, 4.61e-19,             &
      4.34e-19, 4.34e-19, 4.34e-19, 4.34e-19, 4.34e-19, 4.34e-19,             &
      4.08e-19, 4.08e-19, 4.08e-19, 4.08e-19, 4.08e-19, 4.08e-19,             &
      3.87e-19, 3.87e-19, 3.87e-19, 3.87e-19, 3.87e-19, 3.87e-19,             &
      3.67e-19, 3.67e-19, 3.67e-19, 3.67e-19, 3.67e-19, 3.67e-19,             &
      1.87e-22, 1.87e-22, 1.87e-22, 1.87e-22, 1.87e-22, 1.87e-22,             &
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0,                                           &
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0,                                           &
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0,                                           &
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0],                                          &
      (/ size(lower_extrap, 2), size(lower_extrap, 1) /)                      &
    )

    upper_extrap = reshape([ real::                                           &
      6.71e-19, 6.71e-19, 6.71e-19, 6.71e-19, 6.71e-19, 6.71e-19,             &
      6.71e-19, 6.71e-19, 6.71e-19, 6.71e-19, 6.71e-19, 6.71e-19,             &
      6.71e-19, 6.71e-19, 6.71e-19, 6.71e-19, 6.71e-19, 6.71e-19,             &
      6.71e-19, 6.71e-19, 6.71e-19, 6.71e-19, 6.71e-19, 6.71e-19,             &
      6.71e-19, 6.71e-19, 6.71e-19, 6.71e-19, 6.71e-19, 6.71e-19,             &
      6.50e-19, 6.50e-19, 6.50e-19, 6.50e-19, 6.50e-19, 6.50e-19,             &
      6.07e-19, 6.07e-19, 6.07e-19, 6.07e-19, 6.07e-19, 6.07e-19,             &
      5.64e-19, 5.64e-19, 5.64e-19, 5.64e-19, 5.64e-19, 5.64e-19,             &
      5.28e-19, 5.28e-19, 5.28e-19, 5.28e-19, 5.28e-19, 5.28e-19,             &
      4.92e-19, 4.92e-19, 4.92e-19, 4.92e-19, 4.92e-19, 4.92e-19,             &
      4.61e-19, 4.61e-19, 4.61e-19, 4.61e-19, 4.61e-19, 4.61e-19,             &
      4.34e-19, 4.34e-19, 4.34e-19, 4.34e-19, 4.34e-19, 4.34e-19,             &
      4.08e-19, 4.08e-19, 4.08e-19, 4.08e-19, 4.08e-19, 4.08e-19,             &
      3.87e-19, 3.87e-19, 3.87e-19, 3.87e-19, 3.87e-19, 3.87e-19,             &
      3.67e-19, 3.67e-19, 3.67e-19, 3.67e-19, 3.67e-19, 3.67e-19,             &
      211.88, 211.88, 211.88, 211.88, 211.88, 211.88,                         &
      212, 212, 212, 212, 212, 212,                                           &
      212, 212, 212, 212, 212, 212,                                           &
      212, 212, 212, 212, 212, 212,                                           &
      212, 212, 212, 212, 212, 212],                                          &
      (/ size(upper_extrap, 2), size(upper_extrap, 1) /)                      &
    )

    ! load test grids
    call config%from_file( "test/data/grid.180-220.config.json" )
    grids => grid_warehouse_t( config )

    ! load test profiles
    call config%from_file( "test/data/profile.temperature.config.json" )
    profiles => profile_warehouse_t( config, grids )

    ! get cross section config data
    call config%from_file( "test/data/cross_sections/cross_section.h2o2-oh_oh.config.json" )
    call config%get( "cross sections", cs_set, Iam )
    iter => cs_set%get_iterator( )

    ! load and test cross section w/o extrapolation
    call assert( 101264914, iter%next( ) )
    call cs_set%get( iter, cs_config, Iam )
    cross_section => cross_section_h2o2_oh_oh_t( cs_config, grids, profiles )
    results = cross_section%calculate( grids, profiles )
    call check_values( results, no_extrap, .01_dk )
    deallocate( cross_section )

    ! load and test cross section with lower extrapolation
    call assert( 101264915, iter%next( ) )
    call cs_set%get( iter, cs_config, Iam )
    cross_section => cross_section_h2o2_oh_oh_t( cs_config, grids, profiles )
    results = cross_section%calculate( grids, profiles )
    call check_values( results, lower_extrap, .01_dk )
    deallocate( cross_section )

    ! load and test cross section with lower and upper extrapolation
    call assert( 101264916, iter%next( ) )
    call cs_set%get( iter, cs_config, Iam )
    cross_section => cross_section_h2o2_oh_oh_t( cs_config, grids, profiles )
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

  end subroutine test_cross_section_h2o2_oh_oh_t

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program test_cross_section
