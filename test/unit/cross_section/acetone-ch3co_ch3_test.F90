! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> Tests for the base cross_section_t type
program test_cross_section

  use musica_mpi,                      only : musica_mpi_init,                &
                                              musica_mpi_finalize
  use tuvx_cross_section,              only : cross_section_t
  use tuvx_cross_section_ch3coch3_ch3co_ch3
  use tuvx_test_utils,                 only : check_values

  implicit none

  call musica_mpi_init( )
  call test_cross_section_ch3coch3_ch3co_ch3_t( )
  call musica_mpi_finalize( )

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine test_cross_section_ch3coch3_ch3co_ch3_t( )

    use musica_assert,                 only : assert
    use musica_constants,              only : dk => musica_dk
    use musica_config,                 only : config_t
    use musica_iterator,               only : iterator_t
    use musica_mpi
    use musica_string,                 only : string_t
    use tuvx_cross_section_factory,    only : cross_section_type_name,        &
                                              cross_section_allocate
    use tuvx_grid_warehouse,           only : grid_warehouse_t
    use tuvx_profile_warehouse,        only : profile_warehouse_t

    class(grid_warehouse_t),    pointer :: grids
    class(profile_warehouse_t), pointer :: profiles
    class(cross_section_t),     pointer :: cross_section

    character(len=*), parameter :: Iam = "acetone cross section tests"
    type(config_t) :: config, cs_set, cs_config
    class(iterator_t), pointer :: iter
    real(kind=dk), allocatable :: results(:,:)
    real(dk), allocatable :: acetone_no_extrap(:,:)
    real(dk), allocatable :: acetone_lower_extrap(:,:)
    real(dk), allocatable :: acetone_upper_extrap(:,:)
    character, allocatable :: buffer(:)
    type(string_t) :: type_name
    integer :: pos, pack_size
    integer, parameter :: comm = MPI_COMM_WORLD

    allocate(acetone_no_extrap(4, 6))
    allocate(acetone_lower_extrap(4, 5))
    allocate(acetone_upper_extrap(4, 6))

    ! All of these values were produced by one run of the cross section.
    ! So, these tests are testing that any changes don't produce unexpected
    ! changes. The values here are meaningless.
    acetone_no_extrap = reshape([                                             &
      7.0301862E+08,6.5655223E+08,6.1219315E+08,                              &
      5.6989122E+08,5.2959639E+08,4.9126447E+08,                              &
      1.1707670E+09,1.0933860E+09,1.0195143E+09,                              &
      9.4906831E+08,8.8196471E+08,8.1812988E+08,                              &
      1.5530895E+09,1.4504406E+09,1.3524467E+09,                              &
      1.2589971E+09,1.1699813E+09,1.0853017E+09,                              &
      9.3653653E+08,8.7463863E+08,8.1554773E+08,                              &
      7.5919704E+08,7.0551989E+08,6.5445740E+08],                             &
      (/ size(acetone_no_extrap, 2), size(acetone_no_extrap, 1) /)            &
    )

    acetone_lower_extrap = reshape([                                          &
      679520663.87, 634114120.66, 590789797.00,                               &
      549497594.58, 510190457.65, 1131635623.82,                              &
      1056019580.15, 983871012.82, 915106495.77,                              &
      849647669.83, 1501180166.35, 1400872429.35,                             &
      1305164378.50, 1213945348.72, 1127111396.33,                            &
      905234910.36, 844748758.10, 787036198.63,                               &
      732030502.40, 679668992.93],                                            &
      (/ size(acetone_lower_extrap, 2), size(acetone_lower_extrap, 1) /)      &
    )

    acetone_upper_extrap = reshape([                                    &
      703018619.71846068, 656552233.31681216, 612193147.37404847,             &
      569891218.73771667, 529596393.89548576, 491264467.90404767,             &
      1170767016.0155232, 1093386038.4216070, 1019514303.0582997,             &
      949068309.21753716, 881964705.46405852, 818129880.08553731,             &
      1553089515.1841559, 1450440553.7819972, 1352446681.2051401,             &
      1258997133.3324850, 1169981344.0539927, 1085301667.4253523,             &
      7844654649.0400229, 7326251638.3173409, 6831354471.3212347,             &
      6359403928.4437160, 5909841789.7788944, 5482175087.4377661],            &
      (/                                                                      &
        size(acetone_upper_extrap, 2),                                  & 
        size(acetone_upper_extrap, 1)                                   &
      /)                                                                      &
    )


    ! load test grids
    call config%from_file( "test/data/grid.simple.config.json" )
    grids => grid_warehouse_t( config )

    ! load test profiles
    call config%from_file( "test/data/profile.temperature.config.json" )
    profiles => profile_warehouse_t( config, grids )

    ! get cross section config data
    call config%from_file( "test/data/cross_sections/cross_section.acetone.config.json" )
    call config%get( "cross sections", cs_set, Iam )
    iter => cs_set%get_iterator( )

    ! load and test cross section w/o extrapolation
    call assert( 560066370, iter%next( ) )
    call cs_set%get( iter, cs_config, Iam )

    if( musica_mpi_rank( comm ) == 0 ) then
      cross_section =>                                                          &
        cross_section_ch3coch3_ch3co_ch3_t( cs_config, grids, profiles )
      type_name = cross_section_type_name( cross_section )
      pack_size = type_name%pack_size( comm ) + cross_section%pack_size( comm )
      allocate( buffer( pack_size ) )
      pos = 0
      call type_name%mpi_pack(     buffer, pos , comm )
      call cross_section%mpi_pack( buffer, pos , comm )
      call assert( 779377506, pos <= pack_size )
    end if

    call musica_mpi_bcast( pack_size , comm )
    if( musica_mpi_rank( comm ) .ne. 0 ) allocate( buffer( pack_size ) )
    call musica_mpi_bcast( buffer , comm )

    if( musica_mpi_rank( comm ) .ne. 0 ) then
      pos = 0
      call type_name%mpi_unpack( buffer, pos , comm )
      cross_section => cross_section_allocate( type_name )
      call cross_section%mpi_unpack( buffer, pos , comm )
      call assert( 589067964, pos <= pack_size )
    end if
    deallocate( buffer )

    results = cross_section%calculate( grids, profiles )
    call check_values( results, acetone_no_extrap, 0.01_dk )
    deallocate( cross_section )

    ! load and test cross section w/ fixed lower extrapolation and no upper
    ! extrapolation
    call assert( 102622205, iter%next( ) )
    call cs_set%get( iter, cs_config, Iam )

    if( musica_mpi_rank( comm ) == 0 ) then
      cross_section =>                                                        &
        cross_section_ch3coch3_ch3co_ch3_t( cs_config, grids, profiles )
      type_name = cross_section_type_name( cross_section )
      pack_size = type_name%pack_size( comm ) + cross_section%pack_size( comm )
      allocate( buffer( pack_size ) )
      pos = 0
      call type_name%mpi_pack(     buffer, pos , comm )
      call cross_section%mpi_pack( buffer, pos , comm )
      call assert( 973785239, pos <= pack_size )
    end if

    call musica_mpi_bcast( pack_size , comm )
    if( musica_mpi_rank( comm ) .ne. 0 ) allocate( buffer( pack_size ) )
    call musica_mpi_bcast( buffer , comm )

    if( musica_mpi_rank( comm ) .ne. 0 ) then
      pos = 0
      call type_name%mpi_unpack( buffer, pos , comm )
      cross_section => cross_section_allocate( type_name )
      call cross_section%mpi_unpack( buffer, pos , comm )
      call assert( 403570434, pos <= pack_size )
    end if
    deallocate( buffer )

    results = cross_section%calculate( grids, profiles, at_mid_point = .true. )
    call check_values( results, acetone_lower_extrap, 0.01_dk )
    deallocate( cross_section )

    ! load and test cross section w/ extrpolation from lower boundary and
    ! fixed upper extrpolation
    call assert( 101168966, iter%next( ) )
    call cs_set%get( iter, cs_config, Iam )

    if( musica_mpi_rank( comm ) == 0 ) then
      cross_section =>                                                        &
        cross_section_ch3coch3_ch3co_ch3_t( cs_config, grids, profiles )
      type_name = cross_section_type_name( cross_section )
      pack_size = type_name%pack_size( comm ) + cross_section%pack_size( comm )
      allocate( buffer( pack_size ) )
      pos = 0
      call type_name%mpi_pack(     buffer, pos , comm )
      call cross_section%mpi_pack( buffer, pos , comm )
      call assert( 973785239, pos <= pack_size )
    end if

    call musica_mpi_bcast( pack_size , comm )
    if( musica_mpi_rank( comm ) .ne. 0 ) allocate( buffer( pack_size ) )
    call musica_mpi_bcast( buffer , comm )

    if( musica_mpi_rank( comm ) .ne. 0 ) then
      pos = 0
      call type_name%mpi_unpack( buffer, pos , comm )
      cross_section => cross_section_allocate( type_name )
      call cross_section%mpi_unpack( buffer, pos , comm )
      call assert( 403570434, pos <= pack_size )
    end if
    deallocate( buffer )

    results = cross_section%calculate( grids, profiles, at_mid_point = .false.)
    call check_values( results, acetone_upper_extrap, 0.01_dk )
    deallocate( cross_section )

    ! clean up
    deallocate( iter )
    deallocate( grids )
    deallocate( profiles )
    deallocate( acetone_no_extrap )
    deallocate( acetone_lower_extrap )
    deallocate( acetone_upper_extrap )

  end subroutine test_cross_section_ch3coch3_ch3co_ch3_t

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program test_cross_section
