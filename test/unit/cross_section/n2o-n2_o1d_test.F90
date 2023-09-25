! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> Tests for the base cross_section_t type
program test_cross_section

  use tuvx_cross_section, only : cross_section_t
  use tuvx_cross_section_n2o_n2_o1d
  use tuvx_test_utils, only : check_values

  implicit none

  call test_cross_section_n2o_n2_o1d_t( )

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine test_cross_section_n2o_n2_o1d_t( )

    use musica_assert,                 only : assert
    use musica_constants,              only : dk => musica_dk
    use musica_config,                 only : config_t
    use musica_iterator,               only : iterator_t
    use tuvx_grid_warehouse,           only : grid_warehouse_t
    use tuvx_profile_warehouse,        only : profile_warehouse_t

    class(grid_warehouse_t),    pointer :: grids
    class(profile_warehouse_t), pointer :: profiles
    class(cross_section_t),     pointer :: cross_section

    character(len=*), parameter :: Iam = "n2o-n2_o1d cross section tests"
    type(config_t) :: config, cs_set, cs_config
    class(iterator_t), pointer :: iter
    real(kind=dk), allocatable :: results(:,:)
    real(dk), allocatable :: expected(:,:)
    allocate(expected(20, 6))

    ! All of these values were produced by one run of the cross section.
    ! So, these tests are testing that any changes don't produce unexpected
    ! changes. The values here are meaningless.
    expected = reshape([                                                      &
      1.45e-19, 1.44e-19, 1.43e-19, 1.42e-19, 1.41e-19, 1.41e-19,             &
      1.43e-19, 1.42e-19, 1.42e-19, 1.41e-19, 1.40e-19, 1.39e-19,             &
      1.38e-19, 1.37e-19, 1.36e-19, 1.35e-19, 1.34e-19, 1.33e-19,             &
      1.28e-19, 1.27e-19, 1.26e-19, 1.25e-19, 1.24e-19, 1.23e-19,             &
      1.15e-19, 1.14e-19, 1.13e-19, 1.12e-19, 1.11e-19, 1.10e-19,             &
      1.01e-19, 1.00e-19, 9.89e-20, 9.78e-20, 9.67e-20, 9.56e-20,             &
      8.58e-20, 8.47e-20, 8.36e-20, 8.26e-20, 8.15e-20, 8.04e-20,             &
      7.07e-20, 6.97e-20, 6.86e-20, 6.76e-20, 6.66e-20, 6.56e-20,             &
      5.66e-20, 5.57e-20, 5.47e-20, 5.38e-20, 5.28e-20, 5.19e-20,             &
      4.41e-20, 4.32e-20, 4.24e-20, 4.15e-20, 4.07e-20, 3.99e-20,             &
      3.35e-20, 3.27e-20, 3.19e-20, 3.12e-20, 3.05e-20, 2.98e-20,             &
      2.47e-20, 2.41e-20, 2.34e-20, 2.28e-20, 2.22e-20, 2.16e-20,             &
      1.78e-20, 1.73e-20, 1.68e-20, 1.63e-20, 1.58e-20, 1.53e-20,             &
      1.26e-20, 1.21e-20, 1.17e-20, 1.13e-20, 1.09e-20, 1.05e-20,             &
      8.71e-21, 8.36e-21, 8.03e-21, 7.71e-21, 7.40e-21, 7.11e-21,             &
      5.89e-21, 5.62e-21, 5.37e-21, 5.13e-21, 4.90e-21, 4.68e-21,             &
      3.91e-21, 3.71e-21, 3.53e-21, 3.35e-21, 3.18e-21, 3.02e-21,             &
      2.55e-21, 2.41e-21, 2.27e-21, 2.14e-21, 2.02e-21, 1.91e-21,             &
      1.64e-21, 1.54e-21, 1.44e-21, 1.35e-21, 1.27e-21, 1.19e-21,             &
      1.04e-21, 9.77e-22, 9.10e-22, 8.47e-22, 7.89e-22, 7.35e-22],            &
      (/ size(expected, 2), size(expected, 1) /)                              &
    )

    ! load test grids
    call config%from_file( "test/data/grid.180-220.config.json" )
    grids => grid_warehouse_t( config )

    ! load test profiles
    call config%from_file( "test/data/profile.temperature.config.json" )
    profiles => profile_warehouse_t( config, grids )

    ! get cross section config data
    call config%from_file( "test/data/cross_sections/cross_section.hobr-oh_br.config.json" )
    call config%get( "cross sections", cs_set, Iam )
    iter => cs_set%get_iterator( )

    ! load and test cross section w/o extrapolation
    call assert( 101264914, iter%next( ) )
    call cs_set%get( iter, cs_config, Iam )
    cross_section => cross_section_n2o_n2_o1d_t( cs_config, grids, profiles )
    results = cross_section%calculate( grids, profiles )
    call check_values( results, expected, .01_dk )
    deallocate( cross_section )

    ! clean up
    deallocate( iter )
    deallocate( grids )
    deallocate( profiles )
    deallocate( expected )

  end subroutine test_cross_section_n2o_n2_o1d_t

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program test_cross_section
