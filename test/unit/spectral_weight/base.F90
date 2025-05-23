! Copyright (C) 2021-2025 University Corporation for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
program test_spectral_weight_base

  use musica_constants,       only : dk => musica_dk
  use musica_config,          only : config_t
  use tuvx_grid_warehouse,    only : grid_warehouse_t
  use tuvx_profile_warehouse, only : profile_warehouse_t
  use tuvx_spectral_weight,   only : spectral_weight_t
  use musica_assert,          only : assert
  use musica_mpi

  implicit none

  character(len=*), parameter :: conf = 'test/data/spectral_weights/base.json'
  type(config_t)                  :: grid_config, config, profile_config
  type(config_t)                  :: weights_config
  class(grid_warehouse_t), pointer :: grid_warehouse => null()
  class(spectral_weight_t), pointer :: spectral_weight => null()
  class(profile_warehouse_t), pointer :: profile_warehouse => null()

  call config%from_file( conf )

  call config%get( "grids", grid_config, "" )
  call config%get( "profiles", profile_config, "" )
  call config%get( "weights", weights_config, "" )

  grid_warehouse => grid_warehouse_t( grid_config )
  profile_warehouse => profile_warehouse_t( profile_config, grid_warehouse )
  spectral_weight => spectral_weight_t( weights_config, grid_warehouse,       &
    profile_warehouse)

  call musica_mpi_init( )

  call test_spectral_weight_t( spectral_weight )

  deallocate( spectral_weight )
  deallocate( profile_warehouse )
  deallocate( grid_warehouse )

  call musica_mpi_finalize( )

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine test_spectral_weight_t( the_weight )
    ! Currently only tests MPI functions

    use musica_string,                 only : string_t
    use musica_constants,              only : dk => musica_dk
    use tuvx_spectral_weight_factory,  only : spectral_weight_type_name,      &
                                              spectral_weight_allocate
    use tuvx_test_utils,               only : check_values

    class(spectral_weight_t), pointer :: the_weight
    class(spectral_weight_t), pointer :: unpacked

    character, allocatable :: buffer(:)
    type(string_t) :: type_name
    integer :: pos, pack_size
    integer, parameter :: comm = MPI_COMM_WORLD

    ! Get copy of the rayliegh radiator and test MPI functions
    if( musica_mpi_rank( comm ) == 0 ) then
      type_name = spectral_weight_type_name( the_weight )
      pack_size = type_name%pack_size( comm ) + the_weight%pack_size( comm )

      allocate( buffer( pack_size ) )
      pos = 0
      call type_name%mpi_pack(  buffer, pos , comm )
      call the_weight%mpi_pack( buffer, pos , comm )
      call assert( 192131787, pos <= pack_size )
    end if

    call musica_mpi_bcast( pack_size , comm )
    if( musica_mpi_rank( comm ) .ne. 0 ) allocate( buffer( pack_size ) )
    call musica_mpi_bcast( buffer , comm )

    if( musica_mpi_rank( comm ) .ne. 0 ) then
      pos = 0
      call type_name%mpi_unpack( buffer, pos , comm )
      unpacked => spectral_weight_allocate( type_name )
      call unpacked%mpi_unpack( buffer, pos , comm )
      call assert( 862314213, pos > 0 )
      call check_values(                                                      &
        unpacked%spectral_weight_parms(1)%array,                              &
        the_weight%spectral_weight_parms(1)%array,                            &
        0.01_dk )
      deallocate( unpacked )
    end if

  end subroutine test_spectral_weight_t

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program test_spectral_weight_base
