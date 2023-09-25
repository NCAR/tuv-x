! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
program test_cross_section_warehouse

  use musica_mpi,                      only : musica_mpi_init,                &
                                              musica_mpi_finalize
  use tuvx_cross_section_warehouse

  implicit none

  call musica_mpi_init( )
  call test_cross_section_warehouse_t( )
  call musica_mpi_finalize( )

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine test_cross_section_warehouse_t( )

    use musica_assert,                 only : almost_equal, assert
    use musica_constants,              only : dk => musica_dk
    use musica_config,                 only : config_t
    use musica_mpi
    use musica_string,                 only : string_t
    use tuvx_cross_section,            only : cross_section_t
    use tuvx_grid_warehouse,           only : grid_warehouse_t
    use tuvx_profile_warehouse,        only : profile_warehouse_t

    class(grid_warehouse_t),          pointer :: grids
    class(profile_warehouse_t),       pointer :: profiles
    class(cross_section_warehouse_t), pointer :: cross_sections
    class(cross_section_t),           pointer :: cross_section

    character(len=*), parameter :: my_name = "cross section warehouse tests"
    type(config_t) :: config, cs_config
    character, allocatable :: buffer(:)
    type(string_t) :: cs_name
    integer :: pos, pack_size
    integer, parameter :: comm = MPI_COMM_WORLD

    ! load test grids
    call config%from_file( "test/data/grid.simple.config.json" )
    grids => grid_warehouse_t( config )

    ! load test profiles
    call config%from_file( "test/data/profile.temperature.config.json" )
    profiles => profile_warehouse_t( config, grids )

    ! load test cross sections on the primary MPI task
    if( musica_mpi_rank( comm ) == 0 ) then
      call config%from_file( "test/data/cross_section_warehouse.config.json" )
      call config%get( "cross sections", cs_config, my_name )
      cross_sections => cross_section_warehouse_t( cs_config, grids, profiles )
      pack_size = cross_sections%pack_size( comm )
      allocate( buffer( pack_size ) )
      pos = 0
      call cross_sections%mpi_pack( buffer, pos , comm )
      call assert( 996567639, pos <= pack_size )
    end if

    call musica_mpi_bcast( pack_size , comm )
    if( musica_mpi_rank( comm ) .ne. 0 ) allocate( buffer( pack_size ) )
    call musica_mpi_bcast( buffer , comm )

    if( musica_mpi_rank( comm ) .ne. 0 ) then
      pos = 0
      allocate( cross_sections )
      call cross_sections%mpi_unpack( buffer, pos , comm )
      call assert( 416276515, pos <= pack_size )
    end if

    cs_name = "base-extend-fix"
    cross_section => cross_sections%get( cs_name )
    call assert( 911858896, allocated( cross_section%cross_section_parms ) )
    call assert( 720007734, size( cross_section%cross_section_parms ) == 1 )
    associate( params => cross_section%cross_section_parms(1) )
      call assert( 622577720, .not. allocated( params%temperature ) )
      call assert( 452420816, .not. allocated( params%deltaT ) )
      call assert( 564739161, allocated( params%array ) )
      call assert( 394582257, size( params%array, 2 ) == 1 )
    end associate
    deallocate( cross_section )

    cs_name = "acetone-fix-zero"
    cross_section => cross_sections%get( cs_name )
    call assert( 401752098, allocated( cross_section%cross_section_parms ) )
    call assert( 231595194, size( cross_section%cross_section_parms ) == 1 )
    associate( params => cross_section%cross_section_parms(1) )
      call assert( 678963040, .not. allocated( params%temperature ) )
      call assert( 508806136, .not. allocated( params%deltaT ) )
      call assert( 956173982, allocated( params%array ) )
      call assert( 221066580, size( params%array, 2 ) == 4 )
    end associate
    deallocate( cross_section )

    cs_name = "o3"
    cross_section => cross_sections%get( cs_name )
    call assert( 675604267, allocated( cross_section%cross_section_parms ) )
    call assert( 505447363, size( cross_section%cross_section_parms ) == 4 )
    associate( params => cross_section%cross_section_parms(2) )
      call assert( 617765708, allocated( params%temperature ) )
      call assert( 447608804, allocated( params%deltaT ) )
      call assert( 894976650, allocated( params%array ) )
      call assert( 107294996, size( params%array, 2 ) == 4 )
      end associate
    deallocate( cross_section )

    deallocate( grids          )
    deallocate( profiles       )
    deallocate( cross_sections )

  end subroutine test_cross_section_warehouse_t

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program test_cross_section_warehouse

