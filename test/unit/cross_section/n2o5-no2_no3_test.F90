! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> Tests for the base cross_section_t type
program test_cross_section

  use tuvx_cross_section, only : cross_section_t
  use tuvx_cross_section_n2o5_no2_no3
  use tuvx_test_utils, only : check_values

  implicit none

  call test_cross_section_n2o5_no2_no3_t( )

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine test_cross_section_n2o5_no2_no3_t( )

    use musica_assert,                 only : assert
    use musica_constants,              only : dk => musica_dk
    use musica_config,                 only : config_t
    use musica_iterator,               only : iterator_t
    use tuvx_grid_warehouse,           only : grid_warehouse_t
    use tuvx_profile_warehouse,        only : profile_warehouse_t

    class(grid_warehouse_t),    pointer :: grids
    class(profile_warehouse_t), pointer :: profiles
    class(cross_section_t),     pointer :: cross_section

    character(len=*), parameter :: Iam = "n2o5-no2_no3 cross section tests"
    type(config_t) :: config, cs_set, cs_config
    class(iterator_t), pointer :: iter
    real(kind=dk), allocatable :: results(:,:)
    real(dk), allocatable :: no_extrap(:,:)
    real(dk), allocatable :: lower_extrap(:,:)
    real(dk), allocatable :: upper_extrap(:,:)
    allocate(no_extrap(15, 6))
    allocate(lower_extrap(15, 6))
    allocate(upper_extrap(15, 6))

    ! All of these values were produced by one run of the cross section.
    ! So, these tests are testing that any changes don't produce unexpected
    ! changes. The values here are meaningless.
    no_extrap = reshape([ real::                                              &
      1.87e-20, 1.87e-20, 1.87e-20, 1.87e-20, 1.87e-20, 1.87e-20,             &
      2.48e-20, 2.48e-20, 2.48e-20, 2.48e-20, 2.48e-20, 2.48e-20,             &
      1.84e-20, 1.84e-20, 1.84e-20, 1.84e-20, 1.84e-20, 1.84e-20,             &
      8.71e-21, 8.71e-21, 8.71e-21, 8.71e-21, 8.71e-21, 8.71e-21,             &
      0, 0, 0, 0, 0, 0,                                                       &
      0, 0, 0, 0, 0, 0,                                                       &
      0, 0, 0, 0, 0, 0,                                                       &
      0, 0, 0, 0, 0, 0,                                                       &
      0, 0, 0, 0, 0, 0,                                                       &
      0, 0, 0, 0, 0, 0,                                                       &
      0, 0, 0, 0, 0, 0,                                                       &
      0, 0, 0, 0, 0, 0,                                                       &
      0, 0, 0, 0, 0, 0,                                                       &
      0, 0, 0, 0, 0, 0,                                                       &
      0, 0, 0, 0, 0, 0],                                                      &
      (/ size(no_extrap, 2), size(no_extrap, 1) /)                            &
    )

    lower_extrap = reshape([ real::                                           &
      0.04125, 0.04201, 0.04283, 0.04371, 0.04465, 0.04565,                   &
      2.562e-20, 2.609e-20, 2.660e-20, 2.715e-20, 2.773e-20, 2.835e-20,       &
      1.907e-20, 1.943e-20, 1.981e-20, 2.021e-20, 2.065e-20, 2.111e-20,       &
      8.998e-21, 9.165e-21, 9.344e-21, 9.534e-21, 9.739e-21, 9.958e-21,       &
      0, 0, 0, 0, 0, 0,                                                       &
      0, 0, 0, 0, 0, 0,                                                       &
      0, 0, 0, 0, 0, 0,                                                       &
      0, 0, 0, 0, 0, 0,                                                       &
      0, 0, 0, 0, 0, 0,                                                       &
      0, 0, 0, 0, 0, 0,                                                       &
      0, 0, 0, 0, 0, 0,                                                       &
      0, 0, 0, 0, 0, 0,                                                       &
      0, 0, 0, 0, 0, 0,                                                       &
      0, 0, 0, 0, 0, 0,                                                       &
      0, 0, 0, 0, 0, 0],                                                      &
      (/ size(lower_extrap, 2), size(lower_extrap, 1) /)                      &
    )

    upper_extrap = reshape([ real::                                           &
        2.94e-20, 2.79e-20, 2.63e-20, 2.48e-20, 2.33e-20, 2.18e-20,           &
        2.26e-20, 2.14e-20, 2.02e-20, 1.90e-20, 1.79e-20, 1.67e-20,           &
        1.68e-20, 1.59e-20, 1.50e-20, 1.42e-20, 1.33e-20, 1.25e-20,           &
        0.36426, 0.34504, 0.32600, 0.30717, 0.28859, 0.27028,                 &
        0.90464, 0.85320, 0.80247, 0.75255, 0.70350, 0.65543,                 &
        0.89132, 0.83338, 0.77676, 0.72154, 0.66782, 0.61569,                 &
        1.37094, 1.64839, 1.99920, 2.44725, 3.02568, 3.78096,                 &
        1.37113, 1.64876, 1.99983, 2.44825, 3.02721, 3.78326,                 &
        1.37113, 1.64876, 1.99983, 2.44825, 3.02721, 3.78326,                 &
        1.37113, 1.64876, 1.99983, 2.44825, 3.02721, 3.78326,                 &
        1.37113, 1.64876, 1.99983, 2.44825, 3.02721, 3.78326,                 &
        1.37113, 1.64876, 1.99983, 2.44825, 3.02721, 3.78326,                 &
        1.37113, 1.64876, 1.99983, 2.44825, 3.02721, 3.78326,                 &
        1.37113, 1.64876, 1.99983, 2.44825, 3.02721, 3.78326,                 &
        1.37113, 1.64876, 1.99983, 2.44825, 3.02721, 3.78326],                &
      (/ size(upper_extrap, 2), size(upper_extrap, 1) /)                      &
    )

    ! load test grids
    call config%from_file( "test/data/grid.300-375.config.json" )
    grids => grid_warehouse_t( config )

    ! load test profiles
    call config%from_file( "test/data/profile.temperature.config.json" )
    profiles => profile_warehouse_t( config, grids )

    ! get cross section config data
    call config%from_file( "test/data/cross_sections/cross_section.n2o5-no2_no3.config.json" )
    call config%get( "cross sections", cs_set, Iam )
    iter => cs_set%get_iterator( )

    ! load and test cross section w/o extrapolation
    call assert( 101264914, iter%next( ) )
    call cs_set%get( iter, cs_config, Iam )
    cross_section => cross_section_n2o5_no2_no3_t( cs_config, grids, profiles )
    results = cross_section%calculate( grids, profiles )
    call check_values( results, no_extrap, .01_dk )
    deallocate( cross_section )

    ! load and test cross section with lower extrapolation
    call assert( 101264915, iter%next( ) )
    call cs_set%get( iter, cs_config, Iam )
    cross_section => cross_section_n2o5_no2_no3_t( cs_config, grids, profiles )
    results = cross_section%calculate( grids, profiles )
    call check_values( results, lower_extrap, .01_dk )
    deallocate( cross_section )

    ! load and test cross section with lower and upper extrapolation
    call assert( 101264916, iter%next( ) )
    call cs_set%get( iter, cs_config, Iam )
    cross_section => cross_section_n2o5_no2_no3_t( cs_config, grids, profiles )
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

  end subroutine test_cross_section_n2o5_no2_no3_t

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program test_cross_section
