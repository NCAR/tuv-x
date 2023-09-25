! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> Tests for the base cross_section_t type
program test_cross_section

  use tuvx_cross_section, only : cross_section_t, cross_section_ptr
  use tuvx_cross_section_rayliegh
  use tuvx_test_utils, only : check_values
#ifdef MUSICA_USE_OPENMP
  use omp_lib
#endif

  implicit none

  call test_cross_section_rayliegh_t( )

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine test_cross_section_rayliegh_t( )

    use musica_assert,                 only : assert
    use musica_constants,              only : dk => musica_dk
    use musica_config,                 only : config_t
    use musica_iterator,               only : iterator_t
    use tuvx_grid_warehouse,           only : grid_warehouse_t
    use tuvx_profile_warehouse,        only : profile_warehouse_t

    class(grid_warehouse_t),      pointer :: grids
    class(profile_warehouse_t),   pointer :: profiles
    class(cross_section_ptr), allocatable :: cross_sections(:)

    character(len=*), parameter :: Iam = "rayliegh cross section tests"
    type(config_t) :: config, cs_set, cs_config
    class(iterator_t), pointer :: iter
    real(kind=dk), allocatable :: results(:,:)
    real(dk), allocatable :: expected(:,:)
    integer :: i_thread

    allocate(expected(5, 6))

    ! All of these values were produced by one run of the cross section.
    ! So, these tests are testing that any changes don't produce unexpected
    ! changes. The values here are meaningless.
    expected = reshape([real::                                                &
      5.57e-26, 5.57e-26, 5.57e-26, 5.57e-26, 5.57e-26, 5.57e-26,             &
      5.41e-26, 5.41e-26, 5.41e-26, 5.41e-26, 5.41e-26, 5.41e-26,             &
      5.26e-26, 5.26e-26, 5.26e-26, 5.26e-26, 5.26e-26, 5.26e-26,             &
      5.11e-26, 5.11e-26, 5.11e-26, 5.11e-26, 5.11e-26, 5.11e-26,             &
      4.97e-26, 4.97e-26, 4.97e-26, 4.97e-26, 4.97e-26, 4.97e-26],            &
      (/ size(expected, 2), size(expected, 1) /)                              &
    )

    ! load test grids
    call config%from_file( "test/data/grid.300-310.config.json" )
    grids => grid_warehouse_t( config )

    ! load test profiles
    call config%from_file( "test/data/profile.temperature.config.json" )
    profiles => profile_warehouse_t( config, grids )

    ! get cross section config data
    call config%from_file( "test/data/cross_sections/cross_section.rayliegh.config.json" )
    call config%get( "cross sections", cs_set, Iam )
    iter => cs_set%get_iterator( )

    ! load and test cross section w/o extrapolation
    call assert( 101264914, iter%next( ) )
    call cs_set%get( iter, cs_config, Iam )

#ifdef MUSICA_USE_OPENMP
    write(*,*) "Testing with ", omp_get_max_threads( ), " threads"
    allocate( cross_sections( omp_get_max_threads( ) ) )
#else
    write(*,*) "Testing without OpenMP support"
    allocate( cross_sections( 1 ) )
#endif
    do i_thread = 1, size( cross_sections )
      cross_sections( i_thread )%val_ =>                                      &
          cross_section_rayliegh_t( cs_config, grids, profiles )
    end do

    !$omp parallel private( results )
#ifdef MUSICA_OPENMP
    associate( cross_section =>                                               &
                 cross_sections( omp_get_thread_num( ) + 1 )%val_ )
#else
    associate( cross_section => cross_sections( 1 )%val_ )
#endif
      results = cross_section%calculate( grids, profiles )
      call check_values( results, expected, .01_dk )
    end associate
    !$omp end parallel

    do i_thread = 1, size( cross_sections )
      deallocate( cross_sections( i_thread )%val_ )
    end do

    ! clean up
    deallocate( iter )
    deallocate( grids )
    deallocate( profiles )
    deallocate( expected )

  end subroutine test_cross_section_rayliegh_t

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program test_cross_section
