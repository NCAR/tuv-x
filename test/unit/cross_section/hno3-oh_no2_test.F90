! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> Tests for the base cross_section_t type
program test_cross_section

  use tuvx_cross_section, only : cross_section_t
  use tuvx_cross_section_hno3_oh_no2
  use tuvx_test_utils, only : check_values

  implicit none

  call test_cross_section_hno3_oh_no2_t( )

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine test_cross_section_hno3_oh_no2_t( )

    use musica_assert,                 only : assert
    use musica_constants,              only : dk => musica_dk
    use musica_config,                 only : config_t
    use musica_iterator,               only : iterator_t
    use tuvx_grid_warehouse,           only : grid_warehouse_t
    use tuvx_profile_warehouse,        only : profile_warehouse_t

    class(grid_warehouse_t),    pointer :: grids
    class(profile_warehouse_t), pointer :: profiles
    class(cross_section_t),     pointer :: cross_section

    character(len=*), parameter :: Iam = "hno3-oh_no2 cross section tests"
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
        0.00000000E+00, 0.00000000E+00, 0.00000000E+00,                       &
        0.00000000E+00, 0.00000000E+00, 0.00000000E+00,                       &
        0.00000000E+00, 0.00000000E+00, 0.00000000E+00,                       &
        0.00000000E+00, 0.00000000E+00, 0.00000000E+00,                       &
        0.00000000E+00, 0.00000000E+00, 0.00000000E+00,                       &
        0.00000000E+00, 0.00000000E+00, 0.00000000E+00,                       &
        0.00000000E+00, 0.00000000E+00, 0.00000000E+00,                       &
        0.00000000E+00, 0.00000000E+00, 0.00000000E+00,                       &
        6.35272794E-21, 6.28292740E-21, 6.21391491E-21,                       &
        6.14568136E-21, 6.07821774E-21, 6.01152535E-21,                       &
        1.27103729E-17, 1.25707177E-17, 1.24326394E-17,                       &
        1.22961194E-17, 1.21611400E-17, 1.20277036E-17,                       &
        1.14073753E-17, 1.12820368E-17, 1.11581135E-17,                       &
        1.10355888E-17, 1.09144467E-17, 1.07946895E-17,                       &
        1.00060382E-17, 9.89609694E-18, 9.78739694E-18,                       &
        9.67992382E-18, 9.57366339E-18, 9.46861771E-18,                       &
        8.41009301E-18, 8.31903871E-18, 8.22899780E-18,                       &
        8.13995871E-18, 8.05191002E-18, 7.96485376E-18,                       &
        6.68020845E-18, 6.60874219E-18, 6.53806213E-18,                       &
        6.46815940E-18, 6.39902523E-18, 6.33066143E-18,                       &
        5.09031945E-18, 5.03520764E-18, 4.98070919E-18,                       &
        4.92681712E-18, 4.87352449E-18, 4.82083254E-18,                       &
        3.81009039E-18, 3.76785972E-18, 3.72610992E-18,                       &
        3.68483537E-18, 3.64403051E-18, 3.60369606E-18,                       &
        2.74709876E-18, 2.71594406E-18, 2.68515211E-18,                       &
        2.65471858E-18, 2.62463919E-18, 2.59491429E-18,                       &
        1.90080756E-18, 1.87857899E-18, 1.85661710E-18,                       &
        1.83491859E-18, 1.81348025E-18, 1.79230211E-18,                       &
        1.27577108E-18, 1.26003270E-18, 1.24449323E-18,                       &
        1.22915010E-18, 1.21400077E-18, 1.19904502E-18,                       &
        5.35507275E-22, 5.28694862E-22, 5.21971170E-22,                       &
        5.15335016E-22, 5.08785237E-22, 5.02321672E-22,                       &
        0.00000000E+00, 0.00000000E+00, 0.00000000E+00,                       &
        0.00000000E+00, 0.00000000E+00, 0.00000000E+00,                       &
        0.00000000E+00, 0.00000000E+00, 0.00000000E+00,                       &
        0.00000000E+00, 0.00000000E+00, 0.00000000E+00,                       &
        0.00000000E+00, 0.00000000E+00, 0.00000000E+00,                       &
        0.00000000E+00, 0.00000000E+00, 0.00000000E+00,                       &
        0.00000000E+00, 0.00000000E+00, 0.00000000E+00,                       &
        0.00000000E+00, 0.00000000E+00, 0.00000000E+00],                      &
       (/ size(no_extrap, 2), size(no_extrap, 1) /)                           &
    )

    lower_extrap = reshape([ real::                                           &
        1.84878151E+02, 1.82846803E+02, 1.80838391E+02,                       &
        1.78852646E+02, 1.76889309E+02, 1.74948416E+02,                       &
        1.84878151E+02, 1.82846803E+02, 1.80838391E+02,                       &
        1.78852646E+02, 1.76889309E+02, 1.74948416E+02,                       &
        1.84878151E+02, 1.82846803E+02, 1.80838391E+02,                       &
        1.78852646E+02, 1.76889309E+02, 1.74948416E+02,                       &
        1.84878151E+02, 1.82846803E+02, 1.80838391E+02,                       &
        1.78852646E+02, 1.76889309E+02, 1.74948416E+02,                       &
        1.84790334E+02, 1.82759951E+02, 1.80752492E+02,                       &
        1.78767691E+02, 1.76805286E+02, 1.74865315E+02,                       &
        1.27103729E-17, 1.25707177E-17, 1.24326394E-17,                       &
        1.22961194E-17, 1.21611400E-17, 1.20277036E-17,                       &
        1.14073753E-17, 1.12820368E-17, 1.11581135E-17,                       &
        1.10355888E-17, 1.09144467E-17, 1.07946895E-17,                       &
        1.00060382E-17, 9.89609694E-18, 9.78739694E-18,                       &
        9.67992382E-18, 9.57366339E-18, 9.46861771E-18,                       &
        8.41009301E-18, 8.31903871E-18, 8.22899780E-18,                       &
        8.13995871E-18, 8.05191002E-18, 7.96485376E-18,                       &
        6.68020845E-18, 6.60874219E-18, 6.53806213E-18,                       &
        6.46815940E-18, 6.39902523E-18, 6.33066143E-18,                       &
        5.09031945E-18, 5.03520764E-18, 4.98070919E-18,                       &
        4.92681712E-18, 4.87352449E-18, 4.82083254E-18,                       &
        3.81009039E-18, 3.76785972E-18, 3.72610992E-18,                       &
        3.68483537E-18, 3.64403051E-18, 3.60369606E-18,                       &
        2.74709876E-18, 2.71594406E-18, 2.68515211E-18,                       &
        2.65471858E-18, 2.62463919E-18, 2.59491429E-18,                       &
        1.90080756E-18, 1.87857899E-18, 1.85661710E-18,                       &
        1.83491859E-18, 1.81348025E-18, 1.79230211E-18,                       &
        1.27577108E-18, 1.26003270E-18, 1.24449323E-18,                       &
        1.22915010E-18, 1.21400077E-18, 1.19904502E-18,                       &
        5.35507275E-22, 5.28694862E-22, 5.21971170E-22,                       &
        5.15335016E-22, 5.08785237E-22, 5.02321672E-22,                       &
        0.00000000E+00, 0.00000000E+00, 0.00000000E+00,                       &
        0.00000000E+00, 0.00000000E+00, 0.00000000E+00,                       &
        0.00000000E+00, 0.00000000E+00, 0.00000000E+00,                       &
        0.00000000E+00, 0.00000000E+00, 0.00000000E+00,                       &
        0.00000000E+00, 0.00000000E+00, 0.00000000E+00,                       &
        0.00000000E+00, 0.00000000E+00, 0.00000000E+00,                       &
        0.00000000E+00, 0.00000000E+00, 0.00000000E+00,                       &
        0.00000000E+00, 0.00000000E+00, 0.00000000E+00],                      &
       (/ size(lower_extrap, 2), size(lower_extrap, 1) /)                     &
    )

    upper_extrap = reshape([ real::                                           &
        1.33741641E-17, 1.32272156E-17, 1.30819261E-17,                       &
        1.29382766E-17, 1.27962479E-17, 1.26558428E-17,                       &
        1.33741641E-17, 1.32272156E-17, 1.30819261E-17,                       &
        1.29382766E-17, 1.27962479E-17, 1.26558428E-17,                       &
        1.33741641E-17, 1.32272156E-17, 1.30819261E-17,                       &
        1.29382766E-17, 1.27962479E-17, 1.26558428E-17,                       &
        1.33741641E-17, 1.32272156E-17, 1.30819261E-17,                       &
        1.29382766E-17, 1.27962479E-17, 1.26558428E-17,                       &
        1.33741641E-17, 1.32272156E-17, 1.30819261E-17,                       &
        1.29382766E-17, 1.27962479E-17, 1.26558428E-17,                       &
        1.27103729E-17, 1.25707177E-17, 1.24326394E-17,                       &
        1.22961194E-17, 1.21611400E-17, 1.20277036E-17,                       &
        1.14073753E-17, 1.12820368E-17, 1.11581135E-17,                       &
        1.10355888E-17, 1.09144467E-17, 1.07946895E-17,                       &
        1.00060382E-17, 9.89609694E-18, 9.78739694E-18,                       &
        9.67992382E-18, 9.57366339E-18, 9.46861771E-18,                       &
        8.41009301E-18, 8.31903871E-18, 8.22899780E-18,                       &
        8.13995871E-18, 8.05191002E-18, 7.96485376E-18,                       &
        6.68020845E-18, 6.60874219E-18, 6.53806213E-18,                       &
        6.46815940E-18, 6.39902523E-18, 6.33066143E-18,                       &
        5.09031945E-18, 5.03520764E-18, 4.98070919E-18,                       &
        4.92681712E-18, 4.87352449E-18, 4.82083254E-18,                       &
        3.81009039E-18, 3.76785972E-18, 3.72610992E-18,                       &
        3.68483537E-18, 3.64403051E-18, 3.60369606E-18,                       &
        2.74709876E-18, 2.71594406E-18, 2.68515211E-18,                       &
        2.65471858E-18, 2.62463919E-18, 2.59491429E-18,                       &
        1.90080756E-18, 1.87857899E-18, 1.85661710E-18,                       &
        1.83491859E-18, 1.81348025E-18, 1.79230211E-18,                       &
        1.27577108E-18, 1.26003270E-18, 1.24449323E-18,                       &
        1.22915010E-18, 1.21400077E-18, 1.19904502E-18,                       &
        2.07816741E+02, 2.05173017E+02, 2.02563723E+02,                       &
        1.99988400E+02, 1.97446598E+02, 1.94938253E+02,                       &
        2.07925902E+02, 2.05280789E+02, 2.02670124E+02,                       &
        2.00093449E+02, 1.97550312E+02, 1.95040649E+02,                       &
        2.07925902E+02, 2.05280789E+02, 2.02670124E+02,                       &
        2.00093449E+02, 1.97550312E+02, 1.95040649E+02,                       &
        2.07925902E+02, 2.05280789E+02, 2.02670124E+02,                       &
        2.00093449E+02, 1.97550312E+02, 1.95040649E+02,                       &
        2.07925902E+02, 2.05280789E+02, 2.02670124E+02,                       &
        2.00093449E+02, 1.97550312E+02, 1.95040649E+02],                      &
       (/ size(upper_extrap, 2), size(upper_extrap, 1) /)                     &
    )

    ! load test grids
    call config%from_file( "test/data/grid.180-220.config.json" )
    grids => grid_warehouse_t( config )

    ! load test profiles
    call config%from_file( "test/data/profile.temperature.config.json" )
    profiles => profile_warehouse_t( config, grids )

    ! get cross section config data
    call config%from_file( "test/data/cross_sections/cross_section.hno3-oh_no2.config.json" )
    call config%get( "cross sections", cs_set, Iam )
    iter => cs_set%get_iterator( )

    ! load and test cross section w/o extrapolation
    call assert( 101264914, iter%next( ) )
    call cs_set%get( iter, cs_config, Iam )
    cross_section => cross_section_hno3_oh_no2_t( cs_config, grids, profiles )
    results = cross_section%calculate( grids, profiles )
    call check_values( results, no_extrap, .01_dk )
    deallocate( cross_section )

    ! load and test cross section with lower extrapolation
    call assert( 101264915, iter%next( ) )
    call cs_set%get( iter, cs_config, Iam )
    cross_section => cross_section_hno3_oh_no2_t( cs_config, grids, profiles )
    results = cross_section%calculate( grids, profiles )
    
    call check_values( results, lower_extrap, .01_dk )
    deallocate( cross_section )

    ! load and test cross section with upper extrapolation
    call assert( 101264914, iter%next( ) )
    call cs_set%get( iter, cs_config, Iam )
    cross_section => cross_section_hno3_oh_no2_t( cs_config, grids, profiles )
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

  end subroutine test_cross_section_hno3_oh_no2_t

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program test_cross_section
