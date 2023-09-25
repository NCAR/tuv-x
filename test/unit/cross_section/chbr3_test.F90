! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> Tests for the base cross_section_t type
program test_cross_section

  use tuvx_cross_section, only : cross_section_t
  use tuvx_cross_section_chbr3
  use tuvx_test_utils, only : check_values

  implicit none

  call test_cross_section_chbr3_t( )

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine test_cross_section_chbr3_t( )

    use musica_assert,                 only : assert
    use musica_constants,              only : dk => musica_dk
    use musica_config,                 only : config_t
    use musica_iterator,               only : iterator_t
    use tuvx_grid_warehouse,           only : grid_warehouse_t
    use tuvx_profile_warehouse,        only : profile_warehouse_t

    class(grid_warehouse_t),    pointer :: grids
    class(profile_warehouse_t), pointer :: profiles
    class(cross_section_t),     pointer :: cross_section

    character(len=*), parameter :: Iam = "chbr3 cross section tests"
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
      1.89e-21, 1.89e-21, 1.89e-21, 1.89e-21, 1.89e-21, 1.89e-21,             &
      3.79e-18, 3.79e-18, 3.79e-18, 3.79e-18, 3.79e-18, 3.79e-18,             &
      3.55e-18, 3.55e-18, 3.55e-18, 3.55e-18, 3.55e-18, 3.55e-18,             &
      3.58e-18, 3.58e-18, 3.58e-18, 3.58e-18, 3.58e-18, 3.58e-18,             &
      3.79e-18, 3.79e-18, 3.79e-18, 3.79e-18, 3.79e-18, 3.79e-18,             &
      4.04e-18, 4.04e-18, 4.04e-18, 4.04e-18, 4.04e-18, 4.04e-18,             &
      4.24e-18, 4.24e-18, 4.24e-18, 4.24e-18, 4.24e-18, 4.24e-18,             &
      4.36e-18, 4.36e-18, 4.36e-18, 4.36e-18, 4.36e-18, 4.36e-18,             &
      4.42e-18, 4.42e-18, 4.42e-18, 4.42e-18, 4.42e-18, 4.42e-18,             &
      4.48e-18, 4.48e-18, 4.48e-18, 4.48e-18, 4.48e-18, 4.48e-18,             &
      4.59e-18, 4.59e-18, 4.59e-18, 4.59e-18, 4.59e-18, 4.59e-18,             &
      2.45e-21, 2.45e-21, 2.45e-21, 2.45e-21, 2.45e-21, 2.45e-21,             &
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
      179.91, 179.91, 179.91, 179.91, 179.91, 179.91,                         &
      3.79e-18, 3.79e-18, 3.79e-18, 3.79e-18, 3.79e-18, 3.79e-18,             &
      3.55e-18, 3.55e-18, 3.55e-18, 3.55e-18, 3.55e-18, 3.55e-18,             &
      3.58e-18, 3.58e-18, 3.58e-18, 3.58e-18, 3.58e-18, 3.58e-18,             &
      3.79e-18, 3.79e-18, 3.79e-18, 3.79e-18, 3.79e-18, 3.79e-18,             &
      4.04e-18, 4.04e-18, 4.04e-18, 4.04e-18, 4.04e-18, 4.04e-18,             &
      4.24e-18, 4.24e-18, 4.24e-18, 4.24e-18, 4.24e-18, 4.24e-18,             &
      4.36e-18, 4.36e-18, 4.36e-18, 4.36e-18, 4.36e-18, 4.36e-18,             &
      4.42e-18, 4.42e-18, 4.42e-18, 4.42e-18, 4.42e-18, 4.42e-18,             &
      4.48e-18, 4.48e-18, 4.48e-18, 4.48e-18, 4.48e-18, 4.48e-18,             &
      4.59e-18, 4.59e-18, 4.59e-18, 4.59e-18, 4.59e-18, 4.59e-18,             &
      2.45e-21, 2.45e-21, 2.45e-21, 2.45e-21, 2.45e-21, 2.45e-21,             &
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0,                                           &
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0,                                           &
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0,                                           &
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0],                                          &
      (/ size(lower_extrap, 2), size(lower_extrap, 1) /)                      &
    )

    upper_extrap = reshape([                                                  &
      3.99e-18, 3.99e-18, 3.99e-18, 3.99e-18, 3.99e-18, 3.99e-18,             &
      3.99e-18, 3.99e-18, 3.99e-18, 3.99e-18, 3.99e-18, 3.99e-18,             &
      3.99e-18, 3.99e-18, 3.99e-18, 3.99e-18, 3.99e-18, 3.99e-18,             &
      3.99e-18, 3.99e-18, 3.99e-18, 3.99e-18, 3.99e-18, 3.99e-18,             &
      3.99e-18, 3.99e-18, 3.99e-18, 3.99e-18, 3.99e-18, 3.99e-18,             &
      3.79e-18, 3.79e-18, 3.79e-18, 3.79e-18, 3.79e-18, 3.79e-18,             &
      3.55e-18, 3.55e-18, 3.55e-18, 3.55e-18, 3.55e-18, 3.55e-18,             &
      3.58e-18, 3.58e-18, 3.58e-18, 3.58e-18, 3.58e-18, 3.58e-18,             &
      3.79e-18, 3.79e-18, 3.79e-18, 3.79e-18, 3.79e-18, 3.79e-18,             &
      4.04e-18, 4.04e-18, 4.04e-18, 4.04e-18, 4.04e-18, 4.04e-18,             &
      4.24e-18, 4.24e-18, 4.24e-18, 4.24e-18, 4.24e-18, 4.24e-18,             &
      4.36e-18, 4.36e-18, 4.36e-18, 4.36e-18, 4.36e-18, 4.36e-18,             &
      4.42e-18, 4.42e-18, 4.42e-18, 4.42e-18, 4.42e-18, 4.42e-18,             &
      4.48e-18, 4.48e-18, 4.48e-18, 4.48e-18, 4.48e-18, 4.48e-18,             &
      4.59e-18, 4.59e-18, 4.59e-18, 4.59e-18, 4.59e-18, 4.59e-18,             &
      219.88, 219.88, 219.88, 219.88, 219.88, 219.88,                         &
      220.0, 220.0, 220.0, 220.0, 220.0, 220.0,                               &
      220.0, 220.0, 220.0, 220.0, 220.0, 220.0,                               &
      220.0, 220.0, 220.0, 220.0, 220.0, 220.0,                               &
      220.0, 220.0, 220.0, 220.0, 220.0, 220.0],                              &
      (/ size(upper_extrap, 2), size(upper_extrap, 1) /)                      &
    )

    ! load test grids
    call config%from_file( "test/data/grid.180-220.config.json" )
    grids => grid_warehouse_t( config )

    ! load test profiles
    call config%from_file( "test/data/profile.temperature.config.json" )
    profiles => profile_warehouse_t( config, grids )

    ! get cross section config data
    call config%from_file( "test/data/cross_sections/cross_section.chbr3.config.json" )
    call config%get( "cross sections", cs_set, Iam )
    iter => cs_set%get_iterator( )

    ! load and test cross section w/o extrapolation
    call assert( 101264914, iter%next( ) )
    call cs_set%get( iter, cs_config, Iam )
    cross_section => cross_section_chbr3_t( cs_config, grids, profiles )
    results = cross_section%calculate( grids, profiles )
    call check_values( results, no_extrap, .01_dk )
    deallocate( cross_section )

    ! load and test cross section with lower extrapolation
    call assert( 101264914, iter%next( ) )
    call cs_set%get( iter, cs_config, Iam )
    cross_section => cross_section_chbr3_t( cs_config, grids, profiles )
    results = cross_section%calculate( grids, profiles )
    call check_values( results, lower_extrap, .01_dk )
    deallocate( cross_section )

    ! load and test cross section with lower and upper extrapolation
    call assert( 101264914, iter%next( ) )
    call cs_set%get( iter, cs_config, Iam )
    cross_section => cross_section_chbr3_t( cs_config, grids, profiles )
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

  end subroutine test_cross_section_chbr3_t

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program test_cross_section
