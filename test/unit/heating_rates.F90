! Copyright (C) 2024 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
program test_heating_rates

  use musica_mpi,                      only : musica_mpi_init,                &
                                              musica_mpi_finalize
  use tuvx_heating_rates

  implicit none

  call musica_mpi_init( )
  call test_heating_rates_t( )
  call musica_mpi_finalize( )

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> @brief Test the heating rates
  subroutine test_heating_rates_t( )

    use musica_assert,                 only : assert
    use musica_config,                 only : config_t
    use musica_constants,              only : dk => musica_dk
    use musica_mpi,                    only : musica_mpi_bcast,               &
                                              musica_mpi_rank,                &
                                              MPI_COMM_WORLD
    use musica_string,                 only : string_t
    use tuvx_grid_warehouse,           only : grid_warehouse_t
    use tuvx_profile_warehouse,        only : profile_warehouse_t
    use tuvx_test_utils,               only : check_values

    type(heating_rates_t),      pointer :: heating_rates
    class(grid_warehouse_t),    pointer :: grids
    class(profile_warehouse_t), pointer :: profiles

    character(len=*), parameter :: Iam = "heating_rates_t tests"
    type(config_t) :: config, sub_config, reactions_config
    type(string_t), allocatable :: labels(:)
    character, allocatable :: buffer(:)
    integer :: pos, pack_size
    integer, parameter :: comm = MPI_COMM_WORLD

    call config%from_file( "test/data/heating_rates.json" )
    call config%get( "grids", sub_config, Iam )
    grids => grid_warehouse_t( sub_config )
    call config%get( "profiles", sub_config, Iam )
    profiles => profile_warehouse_t( sub_config, grids )

    if( musica_mpi_rank( comm ) == 0 ) then
      call config%get( "reactions", reactions_config, Iam )
      call sub_config%empty( )
      call sub_config%add( "reactions", reactions_config, Iam )
      heating_rates => heating_rates_t( sub_config, grids, profiles )
      pack_size = heating_rates%pack_size( comm )
      allocate( buffer( pack_size ) )
      pos = 0
      call heating_rates%mpi_pack( buffer, pos, comm )
      call assert( 534250649, pos <= pack_size )
    end if

    call musica_mpi_bcast( pack_size, comm )
    if( musica_mpi_rank( comm ) .ne. 0 ) allocate( buffer( pack_size ) )
    call musica_mpi_bcast( buffer, comm )

    if( musica_mpi_rank( comm ) .ne. 0 ) then
      pos = 0
      allocate( heating_rates )
      call heating_rates%mpi_unpack( buffer, pos, comm )
      call assert( 192483602, pos <= pack_size )
    end if
    deallocate( buffer )

    ! check labels
    labels = heating_rates%labels( )
    call assert( 152892147, size(labels) == 2 )
    call assert( 437272930, labels(1) == "jfoo" )
    call assert( 884640776, labels(2) == "jbaz" )

    deallocate( grids )
    deallocate( profiles )
    deallocate( heating_rates )

  end subroutine test_heating_rates_t

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program test_heating_rates
