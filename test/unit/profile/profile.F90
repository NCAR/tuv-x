! Copyright (C) 2021 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> Tests the photo_decomp Profile module

!> Test module for the Profile_t type
program test_Profile

  use musica_mpi,                      only : musica_mpi_init,                &
                                              musica_mpi_finalize
  use musica_string,                   only : string_t

  implicit none

  type(string_t)     :: configFileSpec

  call musica_mpi_init( )
  configFileSpec = 'test/data/profile.test.config.json'
  call test_Profile_t( configFileSpec )
  call musica_mpi_finalize( )

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Tests the Profile_t type
  subroutine test_Profile_t( config_flsp )

    use musica_config,    only : config_t
    use musica_string,    only : string_t
    use musica_assert,    only : assert, almost_equal
    use musica_mpi
    use musica_constants, only : ik => musica_ik, dk => musica_dk
    use tuvx_grid_warehouse, only : grid_warehouse_t
    use tuvx_grid,    only : grid_t
    use tuvx_profile_warehouse, only : Profile_warehouse_t
    use tuvx_profile,           only : profile_t

    !> Arguments
    type(string_t), intent(in) :: config_flsp
    !> Local variables
    character(len=*), parameter :: Iam = 'test_Profile: '
    type(config_t)              :: tst_config, child_config
    type(grid_warehouse_t), pointer :: theGridWarehouse
    class(grid_t), pointer   :: zGrid, lambdaGrid
    type(Profile_warehouse_t), pointer :: theProfileWarehouse
    class(profile_t), pointer      :: aProfile
    type(string_t)                  :: Handle
    character, allocatable :: buffer(:)
    integer :: pos, pack_size
    integer, parameter :: comm = MPI_COMM_WORLD

    !> master configuration -> config type
    call tst_config%from_file( config_flsp%to_char() )

    !> Initialize grid warehouse
    call tst_config%get( "grids", child_config, Iam )
    theGridWarehouse => grid_warehouse_t( child_config )

    !> Get copy of grid
    zGrid => theGridWarehouse%get_grid( "height", "km" )
    call assert( 412238768, zGrid%ncells_ .eq. 120_ik )
    call assert( 412238769, all( zGrid%delta_ .eq. 1._dk ) )

    !> Get copy of wavelength grid
    lambdaGrid => theGridWarehouse%get_grid( "wavelength", "nm" )

    !> Initialize profile warehouse
    if( musica_mpi_rank( comm ) == 0 ) then
      call tst_config%get( "profiles", child_config, Iam )
      theProfileWarehouse =>                                                  &
          Profile_warehouse_t( child_config, theGridWareHouse )
      pack_size = theProfileWarehouse%pack_size( comm )
      allocate( buffer( pack_size ) )
      pos = 0
      call theProfileWarehouse%mpi_pack( buffer, pos , comm )
      call assert( 423920648, pos <= pack_size )
    end if

    call musica_mpi_bcast( pack_size , comm )
    if( musica_mpi_rank( comm ) .ne. 0 ) allocate( buffer( pack_size ) )
    call musica_mpi_bcast( buffer , comm )

    if( musica_mpi_rank( comm ) .ne. 0 ) then
      pos = 0
      allocate( theProfileWarehouse )
      call theProfileWarehouse%mpi_unpack( buffer, pos , comm )
      call assert( 743629523, pos <= pack_size )
    end if
    deallocate( buffer )

    !> Get copy of the Air Profile
        aProfile => theProfileWarehouse%get_profile( "air", "molecule cm-3" )
    call assert( 422238771, all( aProfile%delta_val_ < 0._dk ) )

    deallocate( aProfile )

    !> Get copy of the temperature Profile
    aProfile => theProfileWarehouse%get_profile( "temperature", "K" )
    call assert( 412238778, all( aProfile%edge_val_ < 400._dk ) )
    call assert( 412238772, all( aProfile%edge_val_ > 150._dk ) )
    call assert( 412238773, all( abs(aProfile%delta_val_) < 20._dk ) )

    deallocate( aProfile )

    !> Get copy of the Unit test
    aProfile => theProfileWarehouse%get_profile( "FromConfigValues", "bars" )
    call assert( 412238770, aProfile%mid_val_(1) .eq. 22.35_dk )
    call assert( 412238771, aProfile%mid_val_(2) .eq. 67.8_dk )

    deallocate( aProfile )

    aProfile => theProfileWarehouse%get_profile( "FromConfigUniformValue", "bars" )
    call assert( 412238780, all( aProfile%mid_val_ .eq. 12.3_dk ) )

    deallocate( aProfile )

    aProfile => theProfileWarehouse%get_profile( "Earth-Sun distance", "AU" )

    call assert( 412238374,                                                   &
      almost_equal(aProfile%mid_val_(1), 0.9962_dk, 0.01_dk ))
    call assert( 412238375,                                                   &
      almost_equal(aProfile%edge_val_(1), 0.996229_dk, 0.01_dk ))
    call assert( 412238376,                                                   &
      almost_equal(aProfile%edge_val_(2), 0.996252_dk, 0.01_dk ))

    deallocate( aProfile )

    aProfile => theProfileWarehouse%get_profile( "solar zenith angle",        &
      "degrees" )
    call assert( 412239371,                                                   &
      almost_equal(aProfile%mid_val_(1), 20.127_dk, 0.01_dk ))
    call assert( 412239372,                                                   &
      almost_equal(aProfile%edge_val_(1), 21.523_dk, 0.01_dk ))
    call assert( 412239373,                                                   &
      almost_equal(aProfile%edge_val_(2), 18.732_dk, 0.01_dk ))

    deallocate( aProfile )

    !> Get copy of O2 profile
    aProfile => theProfileWarehouse%get_profile( "O2", "molecule cm-3" )

    deallocate( aProfile )

    !> Get copy of ozone profile
    aProfile => theProfileWarehouse%get_profile( "O3", "molecule cm-3" )

    ! print the ozone profile
    call aProfile%output( zGrid, 6 )

    deallocate( aProfile )

    deallocate( zGrid )
    deallocate( lambdaGrid )
    deallocate( theGridWarehouse )
    deallocate( theProfileWarehouse )

  end subroutine test_Profile_t

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program test_Profile
