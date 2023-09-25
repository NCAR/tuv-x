! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> Tests for the base cross_section_t type
program test_cross_section

  use tuvx_cross_section, only : cross_section_t
  use tuvx_cross_section_no2_tint
  use tuvx_test_utils, only : check_values

  implicit none

  call test_cross_section_no2_tint_t( )

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine test_cross_section_no2_tint_t( )

    use musica_assert,                 only : assert
    use musica_constants,              only : dk => musica_dk
    use musica_config,                 only : config_t
    use musica_iterator,               only : iterator_t
    use tuvx_grid_warehouse,           only : grid_warehouse_t
    use tuvx_profile_warehouse,        only : profile_warehouse_t

    class(grid_warehouse_t),    pointer :: grids
    class(profile_warehouse_t), pointer :: profiles
    class(cross_section_t),     pointer :: cross_section

    character(len=*), parameter :: Iam = "no2_tint cross section tests"
    type(config_t) :: config, cs_set, cs_config
    class(iterator_t), pointer :: iter
    real(kind=dk), allocatable :: results(:,:)
    real(dk), allocatable :: no_extrap(:,:)
    real(dk), allocatable :: lower_extrap(:,:)
    real(dk), allocatable :: upper_extrap(:,:)
    allocate(no_extrap(5, 6))
    allocate(lower_extrap(5, 6))
    allocate(upper_extrap(5, 6))

    ! All of these values were produced by one run of the cross section.
    ! So, these tests are testing that any changes don't produce unexpected
    ! changes. The values here are meaningless.
    no_extrap = reshape([ real::                                              &
      8.37e-20, 8.37e-20, 8.37e-20, 8.37e-20, 8.37e-20, 8.37e-20,             &                             
      1.46e-19, 1.46e-19, 1.45e-19, 1.45e-19, 1.45e-19, 1.45e-19,             &                             
      1.06e-19, 1.06e-19, 1.06e-19, 1.06e-19, 1.06e-19, 1.06e-19,             &                             
      0, 0, 0, 0, 0, 0,                                                       &
      0, 0, 0, 0, 0, 0],                                                      &
      (/ size(no_extrap, 2), size(no_extrap, 1) /)                            &
    )

    lower_extrap = reshape([ real::                                           &
      115.04, 115.04, 115.04, 115.04, 115.04, 115.04,                         &                 
      1.46e-19, 1.46e-19, 1.45e-19, 1.45e-19, 1.45e-19, 1.45e-19,             &                             
      1.06e-19, 1.06e-19, 1.06e-19, 1.06e-19, 1.06e-19, 1.06e-19,             &                             
      0, 0, 0, 0, 0, 0,                                                       &
      0, 0, 0, 0, 0, 0],                                                      &
      (/ size(lower_extrap, 2), size(lower_extrap, 1) /)                      &
    )

    upper_extrap = reshape([ real::                                           &
      1.34e-19, 1.34e-19, 1.34e-19, 1.34e-19, 1.34e-19, 1.34e-19,             &                             
      1.46e-19, 1.46e-19, 1.45e-19, 1.45e-19, 1.45e-19, 1.45e-19,             &                             
      98.80, 98.80, 98.80, 98.80, 98.80, 98.80,                               &           
      310, 310, 310, 310, 310, 310,                                           &
      310, 310, 310, 310, 310, 310],                                          &
      (/ size(upper_extrap, 2), size(upper_extrap, 1) /)                      &
    )

    ! load test grids
    call config%from_file( "test/data/grid.300-310.config.json" )
    grids => grid_warehouse_t( config )

    ! load test profiles
    call config%from_file( "test/data/profile.temperature.config.json" )
    profiles => profile_warehouse_t( config, grids )

    ! get cross section config data
    call config%from_file( "test/data/cross_sections/cross_section.no2_tint.config.json" )
    call config%get( "cross sections", cs_set, Iam )
    iter => cs_set%get_iterator( )

    ! load and test cross section w/o extrapolation
    call assert( 101264914, iter%next( ) )
    call cs_set%get( iter, cs_config, Iam )
    cross_section => cross_section_no2_tint_t( cs_config, grids, profiles )
    results = cross_section%calculate( grids, profiles )
    call check_values( results, no_extrap, .01_dk )
    deallocate( cross_section )

    ! load and test cross section with lower extrapolation
    call assert( 101264915, iter%next( ) )
    call cs_set%get( iter, cs_config, Iam )
    cross_section => cross_section_no2_tint_t( cs_config, grids, profiles )
    results = cross_section%calculate( grids, profiles )
    call check_values( results, lower_extrap, .01_dk )
    deallocate( cross_section )

    ! load and test cross section with upper extrapolation
    call assert( 101264914, iter%next( ) )
    call cs_set%get( iter, cs_config, Iam )
    cross_section => cross_section_no2_tint_t( cs_config, grids, profiles )
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

  end subroutine test_cross_section_no2_tint_t

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program test_cross_section
