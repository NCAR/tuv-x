! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> Tests for the base cross_section_t type
program test_cross_section

  use tuvx_cross_section, only : cross_section_t
  use tuvx_cross_section_cfc11
  use tuvx_test_utils, only : check_values

  implicit none

  call test_cross_section_cfc11_t( )

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine test_cross_section_cfc11_t( )

    use musica_assert,                 only : assert
    use musica_constants,              only : dk => musica_dk
    use musica_config,                 only : config_t
    use musica_iterator,               only : iterator_t
    use tuvx_grid_warehouse,           only : grid_warehouse_t
    use tuvx_profile_warehouse,        only : profile_warehouse_t

    class(grid_warehouse_t),    pointer :: grids
    class(profile_warehouse_t), pointer :: profiles
    class(cross_section_t),     pointer :: cross_section

    character(len=*), parameter :: Iam = "cfc-11 cross section tests"
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
      8.420e-22, 8.398e-22, 8.376e-22, 8.353e-22, 8.331e-22, 8.309e-22,       &
      1.625e-18, 1.618e-18, 1.612e-18, 1.605e-18, 1.599e-18, 1.593e-18,       &
      1.349e-18, 1.342e-18, 1.335e-18, 1.328e-18, 1.321e-18, 1.314e-18,       &
      1.099e-18, 1.091e-18, 1.084e-18, 1.077e-18, 1.070e-18, 1.063e-18,       &
      8.848e-19, 8.779e-19, 8.710e-19, 8.642e-19, 8.575e-19, 8.507e-19,       &
      7.140e-19, 7.075e-19, 7.010e-19, 6.946e-19, 6.883e-19, 6.820e-19,       &
      5.684e-19, 5.624e-19, 5.566e-19, 5.508e-19, 5.451e-19, 5.394e-19,       &
      4.400e-19, 4.349e-19, 4.298e-19, 4.248e-19, 4.198e-19, 4.149e-19,       &
      3.338e-19, 3.294e-19, 3.252e-19, 3.209e-19, 3.168e-19, 3.127e-19,       &
      2.470e-19, 2.435e-19, 2.400e-19, 2.366e-19, 2.332e-19, 2.299e-19,       &
      1.787e-19, 1.759e-19, 1.731e-19, 1.705e-19, 1.678e-19, 1.652e-19,       &
      7.879e-23, 7.747e-23, 7.616e-23, 7.488e-23, 7.363e-23, 7.239e-23,       &
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0,                                           &
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0,                                           &
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0,                                           &
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0],                                          &
      (/ size(no_extrap, 2), size(no_extrap, 1) /)                            &
    )

    lower_extrap = reshape([                                                  &
      180.692, 181.151, 181.610, 182.071, 182.533, 182.995,                   &
      180.337, 180.560, 180.783, 181.006, 181.229, 181.453,                   &
      179.982, 179.970, 179.958, 179.947, 179.935, 179.923,                   &
      179.628, 179.383, 179.138, 178.894, 178.650, 178.407,                   &
      179.189, 178.712, 178.237, 177.763, 177.290, 176.819,                   &
      1.625e-18, 1.618e-18, 1.612e-18, 1.605e-18, 1.599e-18, 1.593e-18,       &
      1.349e-18, 1.342e-18, 1.335e-18, 1.328e-18, 1.321e-18, 1.314e-18,       &
      1.099e-18, 1.091e-18, 1.084e-18, 1.077e-18, 1.070e-18, 1.063e-18,       &
      8.848e-19, 8.779e-19, 8.710e-19, 8.642e-19, 8.575e-19, 8.507e-19,       &
      7.140e-19, 7.075e-19, 7.010e-19, 6.946e-19, 6.883e-19, 6.820e-19,       &
      5.684e-19, 5.624e-19, 5.566e-19, 5.508e-19, 5.451e-19, 5.394e-19,       &
      4.400e-19, 4.349e-19, 4.298e-19, 4.248e-19, 4.198e-19, 4.149e-19,       &
      3.338e-19, 3.294e-19, 3.252e-19, 3.209e-19, 3.168e-19, 3.127e-19,       &
      2.470e-19, 2.435e-19, 2.400e-19, 2.366e-19, 2.332e-19, 2.299e-19,       &
      1.787e-19, 1.759e-19, 1.731e-19, 1.705e-19, 1.678e-19, 1.652e-19,       &
      7.879e-23, 7.747e-23, 7.616e-23, 7.488e-23, 7.363e-23, 7.239e-23,       &
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0,                                           &
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0,                                           &
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0,                                           &
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0],                                          &
      (/ size(lower_extrap, 2), size(lower_extrap, 1) /)                      &
    )

    upper_extrap = reshape([                                                  &
      1.786e-18, 1.791e-18, 1.795e-18, 1.800e-18, 1.805e-18, 1.809e-18,       &
      1.783e-18, 1.785e-18, 1.787e-18, 1.789e-18, 1.792e-18, 1.794e-18,       &
      1.779e-18, 1.779e-18, 1.779e-18, 1.779e-18, 1.779e-18, 1.779e-18,       &
      1.776e-18, 1.773e-18, 1.771e-18, 1.769e-18, 1.766e-18, 1.764e-18,       &
      1.772e-18, 1.768e-18, 1.763e-18, 1.758e-18, 1.754e-18, 1.749e-18,       &
      1.625e-18, 1.618e-18, 1.612e-18, 1.605e-18, 1.599e-18, 1.593e-18,       &
      1.349e-18, 1.342e-18, 1.335e-18, 1.328e-18, 1.321e-18, 1.314e-18,       &
      1.099e-18, 1.091e-18, 1.084e-18, 1.077e-18, 1.070e-18, 1.063e-18,       &
      8.848e-19, 8.779e-19, 8.710e-19, 8.642e-19, 8.575e-19, 8.507e-19,       &
      7.140e-19, 7.075e-19, 7.010e-19, 6.946e-19, 6.883e-19, 6.820e-19,       &
      5.684e-19, 5.624e-19, 5.566e-19, 5.508e-19, 5.451e-19, 5.394e-19,       &
      4.400e-19, 4.349e-19, 4.298e-19, 4.248e-19, 4.198e-19, 4.149e-19,       &
      3.338e-19, 3.294e-19, 3.252e-19, 3.209e-19, 3.168e-19, 3.127e-19,       &
      2.470e-19, 2.435e-19, 2.400e-19, 2.366e-19, 2.332e-19, 2.299e-19,       &
      1.787e-19, 1.759e-19, 1.731e-19, 1.705e-19, 1.678e-19, 1.652e-19,       &
      214.303, 210.699, 207.156, 203.674, 200.251, 196.8883,                  &
      213.994, 210.121, 206.320, 202.589, 198.926, 195.3313,                  &
      213.573, 209.435, 205.379, 201.403, 197.505, 193.6848,                  &
      213.152, 208.752, 204.443, 200.225, 196.095, 192.0522,                  &
      212.733, 208.070, 203.511, 199.053, 194.694, 190.4334],                 &
      (/ size(upper_extrap, 2), size(upper_extrap, 1) /)                      &
    )

    ! load test grids
    call config%from_file( "test/data/grid.180-220.config.json" )
    grids => grid_warehouse_t( config )

    ! load test profiles
    call config%from_file( "test/data/profile.temperature.config.json" )
    profiles => profile_warehouse_t( config, grids )

    ! get cross section config data
    call config%from_file( "test/data/cross_sections/cross_section.cfc-11.config.json" )
    call config%get( "cross sections", cs_set, Iam )
    iter => cs_set%get_iterator( )

    ! load and test cross section w/o extrapolation
    call assert( 101264914, iter%next( ) )
    call cs_set%get( iter, cs_config, Iam )
    cross_section => cross_section_cfc11_t( cs_config, grids, profiles )
    results = cross_section%calculate( grids, profiles )
    call check_values( results, no_extrap, .01_dk )
    deallocate( cross_section )

    ! load and test cross section with lower extrapolation
    call assert( 101264914, iter%next( ) )
    call cs_set%get( iter, cs_config, Iam )
    cross_section => cross_section_cfc11_t( cs_config, grids, profiles )
    results = cross_section%calculate( grids, profiles )
    call check_values( results, lower_extrap, .01_dk )
    deallocate( cross_section )

    ! load and test cross section with lower and upper extrapolation
    call assert( 101264914, iter%next( ) )
    call cs_set%get( iter, cs_config, Iam )
    cross_section => cross_section_cfc11_t( cs_config, grids, profiles )
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

  end subroutine test_cross_section_cfc11_t

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program test_cross_section
