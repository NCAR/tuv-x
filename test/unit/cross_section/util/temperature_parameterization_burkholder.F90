! Copyright (C) 2024 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> Tests for the temperature_parameterization_burkholder_t type
program test_temperature_parameterization_burkholder

  use musica_mpi,                      only : musica_mpi_init,                &
                                              musica_mpi_finalize
  use tuvx_temperature_parameterization_burkholder

  implicit none

  call musica_mpi_init( )
  call test_burkholder_t( )
  call musica_mpi_finalize( )

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine test_burkholder_t( )

    use musica_assert,                 only : assert
    use musica_constants,              only : dk => musica_dk
    use musica_config,                 only : config_t
    use musica_mpi
    use musica_string,                 only : string_t
    use tuvx_test_utils,               only : check_values

    type(temperature_parameterization_burkholder_t) :: burkholder_param
    type(config_t) :: config
    character, allocatable :: buffer(:)
    integer :: pack_size, pos
    integer, parameter :: comm = MPI_COMM_WORLD

    call config%from_file(                                                    &
        "test/data/cross_sections/util/burkholder.config.json" )

    if( musica_mpi_rank( comm ) == 0 ) then
      burkholder_param = temperature_parameterization_burkholder_t( config )
      pack_size = burkholder_param%pack_size( comm )
      allocate( buffer( pack_size ) )
      pos = 0
      call burkholder_param%mpi_pack( buffer, pos, comm )
      call assert( 984049167, pos <= pack_size )
    end if     

    call musica_mpi_bcast( pack_size, comm )
    if( musica_mpi_rank( comm ) .ne. 0 ) allocate( buffer( pack_size ) )
    call musica_mpi_bcast( buffer, comm )

    if( musica_mpi_rank( comm ) .ne. 0 ) then
      pos = 0
      call burkholder_param%mpi_unpack( buffer, pos, comm )
      call assert( 761318011, pos <= pack_size )
    end if
    deallocate( buffer )

    ! Check temperature parameterization data members
    call check_values( 756053704, burkholder_param%wavelengths_,              &
                       (/ 302.0_dk, 304.0_dk, 306.0_dk, 308.0_dk, 310.0_dk /),&
                       1.0e-6_dk )
    call check_values( 973109062, burkholder_param%AA_,                       &
                       (/ 13.3_dk, 14.4_dk, 15.5_dk, 16.6_dk, 17.7_dk /),     &
                       1.0e-6_dk )
    call check_values( 187744433, burkholder_param%BB_,                       &
                       (/ 21.4_dk, 22.3_dk, 23.2_dk, 24.1_dk, 25.0_dk /),     &
                       1.0e-6_dk )
    call assert( 301968312, burkholder_param%A_ == 12.5_dk )
    call assert( 521340695, burkholder_param%B_ == 202.3_dk )
    call assert( 405663577, size( burkholder_param%ranges_ ) == 3 )
    call assert( 800457171, burkholder_param%ranges_(1)%min_temperature_ ==   &
                            0.0_dk )
    call assert( 412833418, burkholder_param%ranges_(1)%max_temperature_ ==   &
                            209.999999999999_dk )
    call assert( 860201264, burkholder_param%ranges_(1)%is_fixed_ .eqv.       &
                            .true. )
    call assert( 125093862, burkholder_param%ranges_(1)%fixed_temperature_ == &
                            210.0_dk )
    call assert( 289986459, burkholder_param%ranges_(2)%min_temperature_ ==   &
                            210.0_dk )
    call assert( 802362705, burkholder_param%ranges_(2)%max_temperature_ ==   &
                            300.0_dk )
    call assert( 967255302, burkholder_param%ranges_(2)%is_fixed_ .eqv.       &
                            .false. )
    call assert( 514623149, burkholder_param%ranges_(2)%fixed_temperature_ == &
                            0.0_dk )
    call assert( 126999396, burkholder_param%ranges_(3)%min_temperature_ ==   &
                            300.00000000001_dk )
    call assert( 574367242, burkholder_param%ranges_(3)%max_temperature_ ==   &
                            huge(1.0_dk) )
    call assert( 739259839, burkholder_param%ranges_(3)%is_fixed_ .eqv.       &
                            .true. )
    call assert( 351636086, burkholder_param%ranges_(3)%fixed_temperature_ == &
                            300.0_dk )

  end subroutine test_burkholder_t

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program test_temperature_parameterization_burkholder