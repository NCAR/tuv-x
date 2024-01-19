! Copyright (C) 2024 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> Tests for the temperature_parameterization_taylor_series_t type
program test_temperature_parameterization_taylor_series

  use musica_mpi,                      only : musica_mpi_init,                &
                                              musica_mpi_finalize
  use tuvx_temperature_parameterization_taylor_series

  implicit none

  call musica_mpi_init( )
  call test_taylor_series_t( )
  call musica_mpi_finalize( )

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine test_taylor_series_t( )

    use musica_assert,                 only : assert
    use musica_constants,              only : dk => musica_dk
    use musica_config,                 only : config_t
    use musica_mpi
    use musica_string,                 only : string_t
    use tuvx_test_utils,               only : check_values

    type(temperature_parameterization_taylor_series_t) :: taylor_param
    type(config_t) :: config
    character, allocatable :: buffer(:)
    integer :: pack_size, pos
    integer, parameter :: comm = MPI_COMM_WORLD

    call config%from_file( "test/data/cross_sections/util/taylor.config.json" )

    if( musica_mpi_rank( comm ) == 0 ) then
      taylor_param = temperature_parameterization_taylor_series_t( config )
      pack_size = taylor_param%pack_size( comm )
      allocate( buffer( pack_size ) )
      pos = 0
      call taylor_param%mpi_pack( buffer, pos, comm )
      call assert( 857895829, pos <= pack_size )
    end if     

    call musica_mpi_bcast( pack_size, comm )
    if( musica_mpi_rank( comm ) .ne. 0 ) allocate( buffer( pack_size ) )
    call musica_mpi_bcast( buffer, comm )

    if( musica_mpi_rank( comm ) .ne. 0 ) then
      pos = 0
      call taylor_param%mpi_unpack( buffer, pos, comm )
      call assert( 960137855, pos <= pack_size )
    end if
    deallocate( buffer )

    ! Check temperature parameterization data members
    call check_values( 405965907, taylor_param%wavelengths_,                  &
                       (/ 302.0_dk, 304.0_dk, 306.0_dk, 308.0_dk, 310.0_dk /),&
                       1.0e-6_dk )
    call check_values( 973109062, taylor_param%sigma_,                        &
                       (/ 13.3_dk, 14.4_dk, 15.5_dk, 16.6_dk, 17.7_dk /),     &
                       1.0e-6_dk )
    call assert( 614095855, size( taylor_param%A_, dim = 1) == 2 )
    call check_values( 178709862, taylor_param%A_(1,:),                       &
                       (/ 21.4_dk, 22.3_dk, 23.2_dk, 24.1_dk, 25.0_dk /),     &
                       1.0e-6_dk )
    call check_values( 842091318, taylor_param%A_(2,:),                       &
                       (/ 6.0_dk, 7.0_dk, 8.0_dk, 9.0_dk, 10.0_dk /),         &
                       1.0e-6_dk )
    call assert( 161915997, taylor_param%base_temperature_ == 295.2_dk )
    call assert( 940974571, taylor_param%min_wavelength_ == 280.5_dk )
    call assert( 992095584, taylor_param%max_wavelength_ == 540.2_dk )
    call assert( 483530406, size( taylor_param%ranges_ ) == 3 )
    call assert( 815221134, taylor_param%ranges_(1)%min_temperature_ ==       &
                            0.0_dk )
    call assert( 182355758, taylor_param%ranges_(1)%max_temperature_ ==       &
                            209.999999999999_dk )
    call assert( 977207253, taylor_param%ranges_(1)%is_fixed_ .eqv. .true. )
    call assert( 242099851, taylor_param%ranges_(1)%fixed_temperature_ ==     &
                            210.0_dk )
    call assert( 689467697, taylor_param%ranges_(2)%min_temperature_ ==       &
                            210.0_dk )
    call assert( 301843944, taylor_param%ranges_(2)%max_temperature_ ==       &
                            300.0_dk )
    call assert( 466736541, taylor_param%ranges_(2)%is_fixed_ .eqv. .false. )
    call assert( 914104387, taylor_param%ranges_(2)%fixed_temperature_ ==     &
                            0.0_dk )
    call assert( 178996985, taylor_param%ranges_(3)%min_temperature_ ==       &
                            300.00000000001_dk )
    call assert( 691373231, taylor_param%ranges_(3)%max_temperature_ ==       &
                            huge(1.0_dk) )
    call assert( 856265828, taylor_param%ranges_(3)%is_fixed_ .eqv. .true. )
    call assert( 403633675, taylor_param%ranges_(3)%fixed_temperature_ ==     &
                            300.0_dk )

  end subroutine test_taylor_series_t

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program test_temperature_parameterization_taylor_series