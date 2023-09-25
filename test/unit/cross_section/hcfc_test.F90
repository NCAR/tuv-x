! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> Tests for the base cross_section_t type
program test_cross_section

  use tuvx_cross_section, only : cross_section_t
  use tuvx_cross_section_hcfc
  use tuvx_test_utils, only : check_values

  implicit none

  call test_cross_section_hcfc_t( )

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine test_cross_section_hcfc_t( )

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
    real(dk), allocatable :: expected(:,:)
    allocate(expected(20, 6))

    ! All of these values were produced by one run of the cross section.
    ! So, these tests are testing that any changes don't produce unexpected
    ! changes. The values here are meaningless.
    expected = reshape([                                                     &
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0,                                           &
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0,                                           &
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0,                                           &
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0,                                           &
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0,                                           &
      4.58e-19, 4.56e-19, 4.54e-19, 4.51e-19, 4.48e-19, 4.45e-19,             &
      3.42e-19, 3.39e-19, 3.37e-19, 3.35e-19, 3.32e-19, 3.29e-19,             &
      2.50e-19, 2.48e-19, 2.45e-19, 2.43e-19, 2.41e-19, 2.38e-19,             &
      1.79e-19, 1.77e-19, 1.75e-19, 1.73e-19, 1.71e-19, 1.69e-19,             &
      1.27e-19, 1.25e-19, 1.23e-19, 1.22e-19, 1.20e-19, 1.18e-19,             &
      8.87e-20, 8.73e-20, 8.59e-20, 8.45e-20, 8.31e-20, 8.16e-20,             &
      6.13e-20, 6.02e-20, 5.90e-20, 5.79e-20, 5.68e-20, 5.56e-20,             &
      4.20e-20, 4.11e-20, 4.02e-20, 3.93e-20, 3.85e-20, 3.77e-20,             &
      2.85e-20, 2.79e-20, 2.72e-20, 2.66e-20, 2.60e-20, 2.54e-20,             &
      1.93e-20, 1.88e-20, 1.84e-20, 1.79e-20, 1.75e-20, 1.70e-20,             &
      1.31e-20, 1.27e-20, 1.24e-20, 1.20e-20, 1.17e-20, 1.15e-20,             &
      8.86e-21, 8.61e-21, 8.38e-21, 8.16e-21, 7.96e-21, 7.77e-21,             &
      6.00e-21, 5.84e-21, 5.68e-21, 5.54e-21, 5.41e-21, 5.28e-21,             &
      4.09e-21, 3.98e-21, 3.88e-21, 3.79e-21, 3.70e-21, 3.63e-21,             &
      2.80e-21, 2.73e-21, 2.67e-21, 2.61e-21, 2.56e-21, 2.52e-21],            &
      (/ size(expected, 2), size(expected, 1) /)                            &
    )

    ! load test grids
    call config%from_file( "test/data/grid.180-220.config.json" )
    grids => grid_warehouse_t( config )

    ! load test profiles
    call config%from_file( "test/data/profile.temperature.config.json" )
    profiles => profile_warehouse_t( config, grids )

    ! get cross section config data
    call config%from_file( "test/data/cross_sections/cross_section.hcfc.config.json" )
    call config%get( "cross sections", cs_set, Iam )
    iter => cs_set%get_iterator( )

    ! load and test cross section w/o extrapolation
    call assert( 101264914, iter%next( ) )
    call cs_set%get( iter, cs_config, Iam )
    cross_section => cross_section_hcfc_t( cs_config, grids, profiles )
    results = cross_section%calculate( grids, profiles )
    call check_values( results, expected, .01_dk )
    deallocate( cross_section )

    ! clean up
    deallocate( iter )
    deallocate( grids )
    deallocate( profiles )
    deallocate( expected )

  end subroutine test_cross_section_hcfc_t

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program test_cross_section
