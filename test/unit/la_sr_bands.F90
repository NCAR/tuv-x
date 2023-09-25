! Copyright (C) 2021 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0

program test_la_sr_bands
  !> Tests for the ls sr bands

  use musica_assert,          only : assert
  use musica_constants,       only : dk => musica_dk
  use musica_config,          only : config_t
  use musica_mpi
  use tuvx_grid_warehouse,    only : grid_warehouse_t
  use tuvx_profile_warehouse, only : profile_warehouse_t
  use tuvx_la_sr_bands,       only : la_sr_bands_t
  use tuvx_test_utils,        only : check_values

  implicit none

  character(len=*), parameter :: conf_l = 'test/data/la_srb_bands.config.json'
  type(config_t)                  :: grid_config, la_config, profile_config
  class(grid_warehouse_t), pointer :: grid_warehouse => null()
  class(la_sr_bands_t), pointer    :: la_sr_bands_ => null()
  class(profile_warehouse_t), pointer :: profile_warehouse => null()
  character, allocatable :: buffer(:)
  integer :: pos, pack_size
  integer, parameter :: comm = MPI_COMM_WORLD

  call musica_mpi_init( )

  call la_config%from_file( conf_l )
  call la_config%get( "grids", grid_config, "" )
  call la_config%get( "profiles", profile_config, "" )

  grid_warehouse => grid_warehouse_t( grid_config )
  profile_warehouse => profile_warehouse_t( profile_config, grid_warehouse )

  if( musica_mpi_rank( comm ) == 0 ) then
    la_sr_bands_ =>                                                           &
        la_sr_bands_t( la_config, grid_warehouse, profile_warehouse )
    pack_size = la_sr_bands_%pack_size( comm )
    allocate( buffer( pack_size ) )
    pos = 0
    call la_sr_bands_%mpi_pack( buffer, pos , comm )
    call assert( 328389002, pos <= pack_size )
  end if

  call musica_mpi_bcast( pack_size , comm )
  if( musica_mpi_rank( comm ) .ne. 0 ) allocate( buffer( pack_size ) )
  call musica_mpi_bcast( buffer , comm )

  if( musica_mpi_rank( comm ) .ne. 0 ) then
    pos = 0
    allocate( la_sr_bands_ )
    call la_sr_bands_%mpi_unpack( buffer, pos , comm )
    call assert( 301066523, pos <= pack_size )
  end if
  deallocate( buffer )

  call test_optical_depth( )

  deallocate( grid_warehouse )
  deallocate( profile_warehouse )
  deallocate( la_sr_bands_ )
  call musica_mpi_finalize( )

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine test_optical_depth( )
    use tuvx_grid, only : grid_t
    use tuvx_spherical_geometry, only : spherical_geometry_t

    real(dk), allocatable :: air_vertical_column(:), air_slant_column(:)
    real(dk), allocatable :: o2_optical_depth(:,:)
    class(grid_t), pointer :: height_grid => null( ) ! specified altitude working grid [km]
    type(spherical_geometry_t) :: spherical_geometry

    height_grid => grid_warehouse%get_grid( "height", "km" )
    allocate( air_vertical_column( height_grid%ncells_ ),                      &
              air_slant_column( height_grid%ncells_ + 1 ) )

    allocate( o2_optical_depth(120, 38) )

    air_vertical_column(:) = 1
    air_slant_column(:) = 3
    o2_optical_depth(:,:) = 0

    ! just checking that it runs. This method apparently requires at least
    ! 18 columns of output if you are in the Schumann-Runge wavelegth
    ! running is good enough
    call la_sr_bands_%optical_depth( grid_warehouse, profile_warehouse,       &
      air_vertical_column, air_slant_column, o2_optical_depth,                &
      spherical_geometry )

    deallocate( height_grid )
    deallocate( o2_optical_depth )
    deallocate( air_vertical_column )
    deallocate( air_slant_column )

  end subroutine test_optical_depth

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program test_la_sr_bands
