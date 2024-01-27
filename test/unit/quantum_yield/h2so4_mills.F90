! Copyright (C) 2024 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
program test_quantum_yield_h2so4_mills

  use musica_mpi
  use tuvx_quantum_yield_h2so4_mills

  implicit none

  call musica_mpi_init( )
  call test_quantum_yield_h2so4_mills_t( )
  call musica_mpi_finalize( )

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Test the temperature and pressure dependent H2SO4 quantum yield
  !! calculations against previously generated results from an older version
  !! of TUV
  subroutine test_quantum_yield_h2so4_mills_t( )

    use musica_assert,                 only : assert, die
    use musica_config,                 only : config_t
    use musica_constants,              only : dk => musica_dk
    use musica_string,                 only : string_t
    use tuvx_grid_warehouse,           only : grid_warehouse_t
    use tuvx_profile_warehouse,        only : profile_warehouse_t
    use tuvx_quantum_yield,            only : quantum_yield_t
    use tuvx_quantum_yield_factory
    use tuvx_test_utils,               only : check_values

    character(len=*), parameter :: my_name = "h2so4 quantum yield test"
    class(grid_warehouse_t),    pointer :: grids
    class(profile_warehouse_t), pointer :: profiles
    class(quantum_yield_t),     pointer :: quantum_yield
    type(config_t) :: config, qy_config, grids_config, profiles_config

    character, allocatable :: buffer(:)
    type(string_t) :: type_name
    integer :: pos, pack_size
    integer, parameter :: comm = MPI_COMM_WORLD

    call config%from_file( "test/data/quantum_yields/h2so4_mills.config.json" )
    call config%get( "grids",         grids_config,    my_name )
    call config%get( "profiles",      profiles_config, my_name )
    call config%get( "quantum yield", qy_config,       my_name )

    grids => grid_warehouse_t( grids_config )
    profiles => profile_warehouse_t( profiles_config, grids )

    if( musica_mpi_rank( comm ) == 0 ) then
      quantum_yield => quantum_yield_builder( qy_config, grids, profiles )
      type_name = quantum_yield_type_name( quantum_yield )
      pack_size = type_name%pack_size( comm ) + quantum_yield%pack_size( comm )
      allocate( buffer( pack_size ) )
      pos = 0
      call type_name%mpi_pack(     buffer, pos, comm )
      call quantum_yield%mpi_pack( buffer, pos, comm )
      call assert( 837477723, pos <= pack_size )
    end if

    call musica_mpi_bcast( pack_size, comm )
    if( musica_mpi_rank( comm ) .ne. 0 ) allocate( buffer( pack_size ) )
    call musica_mpi_bcast( buffer, comm )

    if( musica_mpi_rank( comm ) .ne. 0 ) then
      pos = 0
      call type_name%mpi_unpack( buffer, pos, comm )
      quantum_yield => quantum_yield_allocate( type_name )
      call quantum_yield%mpi_unpack( buffer, pos, comm )
      call assert( 599405941, pos <= pack_size )
    end if
    deallocate( buffer )

    ! clean up
    deallocate( grids )
    deallocate( profiles )
    deallocate( quantum_yield )
      
  end subroutine test_quantum_yield_h2so4_mills_t

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program test_quantum_yield_h2so4_mills