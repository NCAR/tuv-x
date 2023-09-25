! Copyright (C) 2021 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0

program test_radiative_transfer_mpi

  use musica_mpi,              only : musica_mpi_init,                        &
                                      musica_mpi_finalize, MPI_COMM_WORLD
  use musica_config,           only : config_t
  use tuvx_grid_warehouse,     only : grid_warehouse_t
  use tuvx_profile_warehouse,  only : profile_warehouse_t
  use tuvx_radiative_transfer, only : radiative_transfer_t

  character(len=*), parameter :: conf = 'test/data/radiative_transfer.config.json'
  type(config_t)                  :: grid_config, config, profile_config
  type(config_t)                  :: rad_config
  class(grid_warehouse_t), pointer :: grid_warehouse => null()
  class(profile_warehouse_t), pointer :: profile_warehouse => null()
  class(radiative_transfer_t), pointer :: radiative_transfer => null()

  call config%from_file( conf )

  call config%get( "grids", grid_config, "" )
  call config%get( "profiles", profile_config, "" )
  call config%get( "radiative transfer", rad_config, "" )

  call musica_mpi_init( )

  grid_warehouse => grid_warehouse_t( grid_config )
  profile_warehouse => profile_warehouse_t( profile_config, grid_warehouse )
  radiative_transfer => radiative_transfer_t( rad_config, grid_warehouse,       &
    profile_warehouse)

  call test_mpi( radiative_transfer )

  deallocate( profile_warehouse )
  deallocate( grid_warehouse )
  deallocate( radiative_transfer )

  call musica_mpi_finalize( )

contains

  subroutine test_mpi( the_radiative_transfer )
    use musica_mpi,    only : musica_mpi_bcast, musica_mpi_rank
    use musica_assert, only : assert

    class(radiative_transfer_t), pointer, intent(in) :: the_radiative_transfer
    class(radiative_transfer_t), pointer  :: unpacked

    character, allocatable :: buffer(:)
    integer :: pos, pack_size
    integer, parameter :: comm = MPI_COMM_WORLD

    ! test MPI functions
    if( musica_mpi_rank( comm ) == 0 ) then
      pack_size = the_radiative_transfer%pack_size( comm )

      allocate( buffer( pack_size ) )
      pos = 0

      call the_radiative_transfer%mpi_pack( buffer, pos , comm )

      call assert( 793189081, pos <= pack_size )
    end if

    call musica_mpi_bcast( pack_size , comm )
    if( musica_mpi_rank( comm ) .ne. 0 ) allocate( buffer( pack_size ) )
    call musica_mpi_bcast( buffer , comm )

    if( musica_mpi_rank( comm ) .ne. 0 ) then
      pos = 0
      allocate( unpacked )
      call unpacked%mpi_unpack( buffer, pos , comm )

      call assert( 871384211, pos > 0 )
      call assert( 132313383, pos <= pack_size )
      deallocate( unpacked )
    end if
  end subroutine test_mpi

end program test_radiative_transfer_mpi
