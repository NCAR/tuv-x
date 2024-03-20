! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> Tests for the base cross_section_t type
program test_cross_section

  use tuvx_cross_section, only : cross_section_t
  use tuvx_cross_section_rono2
  use tuvx_test_utils, only : check_values

  implicit none

  call test_cross_section_rono2_t( )

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine test_cross_section_rono2_t( )

    use musica_assert,                 only : assert
    use musica_constants,              only : dk => musica_dk
    use musica_config,                 only : config_t
    use musica_iterator,               only : iterator_t
    use tuvx_grid_warehouse,           only : grid_warehouse_t
    use tuvx_profile_warehouse,        only : profile_warehouse_t

    class(grid_warehouse_t),    pointer :: grids
    class(profile_warehouse_t), pointer :: profiles
    class(cross_section_t),     pointer :: cross_section

    character(len=*), parameter :: Iam = "rono2 cross section tests"
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
  3.3302044212254328E-024, 3.2161952124512202E-024, 3.1061223947672814E-024,  &
  2.9998489273903812E-024, 2.8972425790838159E-024, 2.7981907562917913E-024,  &
  4.0122656353068610E-021, 3.8717595318338420E-021, 3.7362148177227141E-021,  &
  3.6054548800794594E-021, 3.4793094427501407E-021, 3.3576327537719567E-021,  &
  3.2269583072774368E-021, 3.1085946499675142E-021, 2.9946069663498917E-021,  &
  2.8848322221146335E-021, 2.7791135029650159E-021, 2.6773151777054408E-021,  &
  2.5576809099900323E-021, 2.4589071748754612E-021, 2.3639765819585908E-021,  &
  2.2727385030188113E-021, 2.1850482559046444E-021, 2.1007795988510739E-021,  &
  1.7524893730671802E-024, 1.6829505050933260E-024, 1.6161910805734185E-024,  &
  1.5520992136486835E-024, 1.4905675538255535E-024, 1.4315020191061525E-024], &
      (/ size(no_extrap, 2), size(no_extrap, 1) /)                            &
    )

    lower_extrap = reshape([ real::                                           &
  284.35719862387975, 274.62225892536395, 265.22343707474488,                 &
  256.14903152812633, 247.38775141585484, 238.92998267703231,                 &
  4.0122656353068610E-021, 3.8717595318338420E-021, 3.7362148177227141E-021,  &
  3.6054548800794594E-021, 3.4793094427501407E-021, 3.3576327537719567E-021,  &
  3.2269583072774368E-021, 3.1085946499675142E-021, 2.9946069663498917E-021,  &
  2.8848322221146335E-021, 2.7791135029650159E-021, 2.6773151777054408E-021,  &
  2.5576809099900323E-021, 2.4589071748754612E-021, 2.3639765819585908E-021,  &
  2.2727385030188113E-021, 2.1850482559046444E-021, 2.1007795988510739E-021,  &
  1.7524893730671802E-024, 1.6829505050933260E-024, 1.6161910805734185E-024,  &
  1.5520992136486835E-024, 1.4905675538255535E-024, 1.4315020191061525E-024], &
      (/ size(lower_extrap, 2), size(lower_extrap, 1) /)                      &
    )

    upper_extrap = reshape([ real::                                           &
  4.4108667830913617E-021, 4.2598612085555927E-021, 4.1140693970533126E-021,  &
  3.9733098376133487E-021, 3.8374073895248339E-021, 3.7062129222502928E-021,  &
  4.0122656353068610E-021, 3.8717595318338420E-021, 3.7362148177227141E-021,  &
  3.6054548800794594E-021, 3.4793094427501407E-021, 3.3576327537719567E-021,  &
  3.2269583072774368E-021, 3.1085946499675142E-021, 2.9946069663498917E-021,  &
  2.8848322221146335E-021, 2.7791135029650159E-021, 2.6773151777054408E-021,  &
  2.5576809099900323E-021, 2.4589071748754612E-021, 2.3639765819585908E-021,  &
  2.2727385030188113E-021, 2.1850482559046444E-021, 2.1007795988510739E-021,  &
  291.32413139180113, 279.76437495516359, 268.66665781098493,                 &
  258.01238067351579, 247.78369819101036, 237.96497069299082],                &
      (/ size(upper_extrap, 2), size(upper_extrap, 1) /)                      &
    )
    ! load test grids
    call config%from_file( "test/data/grid.300-310.config.json" )
    grids => grid_warehouse_t( config )

    ! load test profiles
    call config%from_file( "test/data/profile.temperature.config.json" )
    profiles => profile_warehouse_t( config, grids )

    ! get cross section config data
    call config%from_file( "test/data/cross_sections/cross_section.rono2.config.json" )
    call config%get( "cross sections", cs_set, Iam )
    iter => cs_set%get_iterator( )

    ! load and test cross section w/o extrapolation
    call assert( 101264914, iter%next( ) )
    call cs_set%get( iter, cs_config, Iam )
    cross_section => cross_section_rono2_t( cs_config, grids, profiles )
    results = cross_section%calculate( grids, profiles )
    call check_values( results, no_extrap, .01_dk )
    deallocate( cross_section )

    ! load and test cross section with lower extrapolation
    call assert( 101264915, iter%next( ) )
    call cs_set%get( iter, cs_config, Iam )
    cross_section => cross_section_rono2_t( cs_config, grids, profiles )
    results = cross_section%calculate( grids, profiles )
    call check_values( results, lower_extrap, .01_dk )
    deallocate( cross_section )

    ! load and test cross section with upper extrapolation
    call assert( 101264914, iter%next( ) )
    call cs_set%get( iter, cs_config, Iam )
    cross_section => cross_section_rono2_t( cs_config, grids, profiles )
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

  end subroutine test_cross_section_rono2_t

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program test_cross_section
