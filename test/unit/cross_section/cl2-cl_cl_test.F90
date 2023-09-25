! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> Tests for the base cross_section_t type
program test_cross_section

  use tuvx_cross_section, only : cross_section_t
  use tuvx_cross_section_cl2_cl_cl
  use tuvx_test_utils, only : check_values

  implicit none

  call test_cross_section_cl2_cl_cl_t( )

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine test_cross_section_cl2_cl_cl_t( )

    use musica_assert,                 only : assert, assert_msg
    use musica_constants,              only : dk => musica_dk
    use musica_config,                 only : config_t
    use musica_iterator,               only : iterator_t
    use tuvx_grid_warehouse,           only : grid_warehouse_t
    use tuvx_profile_warehouse,        only : profile_warehouse_t

    class(grid_warehouse_t),    pointer :: grids
    class(profile_warehouse_t), pointer :: profiles
    class(cross_section_t),     pointer :: cross_section

    character(len=*), parameter :: Iam = "cl2_cl_cl cross section tests"
    type(config_t) :: config, cs_set, cs_config
    class(iterator_t), pointer :: iter
    real(kind=dk), allocatable :: results(:,:)
    real(dk), allocatable :: expected(:,:)
    allocate(expected(20, 6))

    ! All of these values were produced by one run of the cross section.
    ! So, these tests are testing that any changes don't produce unexpected
    ! changes. The values here are meaningless.
    expected = reshape([                                                      &
      5.70e-33, 4.49e-33, 3.55e-33, 2.82e-33, 2.24e-33, 1.79e-33,             &
      1.78e-32, 1.42e-32, 1.13e-32, 9.07e-33, 7.29e-33, 5.88e-33,             &
      5.42e-32, 4.35e-32, 3.49e-32, 2.82e-32, 2.28e-32, 1.85e-32,             &
      1.59e-31, 1.28e-31, 1.04e-31, 8.49e-32, 6.93e-32, 5.68e-32,             &
      4.53e-31, 3.69e-31, 3.02e-31, 2.47e-31, 2.03e-31, 1.68e-31,             &
      1.25e-30, 1.02e-30, 8.46e-31, 6.99e-31, 5.79e-31, 4.82e-31,             &
      3.35e-30, 2.77e-30, 2.30e-30, 1.91e-30, 1.60e-30, 1.34e-30,             &
      8.72e-30, 7.27e-30, 6.08e-30, 5.09e-30, 4.28e-30, 3.61e-30,             &
      2.20e-29, 1.85e-29, 1.56e-29, 1.31e-29, 1.11e-29, 9.47e-30,             &
      5.44e-29, 4.60e-29, 3.90e-29, 3.31e-29, 2.82e-29, 2.41e-29,             &
      1.30e-28, 1.11e-28, 9.48e-29, 8.10e-29, 6.95e-29, 5.98e-29,             &
      3.04e-28, 2.61e-28, 2.24e-28, 1.93e-28, 1.66e-28, 1.44e-28,             &
      6.94e-28, 5.99e-28, 5.17e-28, 4.48e-28, 3.89e-28, 3.38e-28,             &
      1.54e-27, 1.34e-27, 1.16e-27, 1.01e-27, 8.86e-28, 7.76e-28,             &
      3.35e-27, 2.92e-27, 2.55e-27, 2.24e-27, 1.97e-27, 1.73e-27,             &
      7.10e-27, 6.24e-27, 5.49e-27, 4.83e-27, 4.27e-27, 3.78e-27,             &
      1.47e-26, 1.30e-26, 1.15e-26, 1.01e-26, 9.05e-27, 8.05e-27,             &
      2.98e-26, 2.65e-26, 2.35e-26, 2.10e-26, 1.87e-26, 1.67e-26,             &
      5.92e-26, 5.28e-26, 4.72e-26, 4.23e-26, 3.79e-26, 3.41e-26,             &
      1.15e-25, 1.03e-25, 9.27e-26, 8.34e-26, 7.52e-26, 6.79e-26],            &
      (/ size(expected, 2), size(expected, 1) /)                              &
    )

    ! load test grids
    call config%from_file( "test/data/grid.180-220.config.json" )
    grids => grid_warehouse_t( config )

    ! load test profiles
    call config%from_file( "test/data/profile.temperature.config.json" )
    profiles => profile_warehouse_t( config, grids )

    ! get cross section config data
    call config%from_file( "test/data/cross_sections/cross_section.cl2_cl_cl.config.json" )
    call config%get( "cross sections", cs_set, Iam )
    iter => cs_set%get_iterator( )

    ! load and test cross section w/o extrapolation
    call assert( 101268914, iter%next( ) )
    call cs_set%get( iter, cs_config, Iam )
    cross_section => cross_section_cl2_cl_cl_t( cs_config, grids, profiles )
    results = cross_section%calculate( grids, profiles )
    call check_values( results, expected, .01_dk )
    deallocate( cross_section )

    ! clean up
    deallocate( iter )
    deallocate( grids )
    deallocate( profiles )
    deallocate( expected )

  end subroutine test_cross_section_cl2_cl_cl_t

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program test_cross_section
