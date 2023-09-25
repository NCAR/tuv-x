! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
program test_quantum_yield_no2_tint

  use musica_mpi,                      only : musica_mpi_init,                &
                                              musica_mpi_finalize
  use tuvx_quantum_yield_no2_tint

  implicit none

  call musica_mpi_init( )
  call test_quantum_yield_no2_tint_t( )
  call musica_mpi_finalize( )

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine test_quantum_yield_no2_tint_t( )
    ! Test functionality of the :f:type:`~tuvx_quantum_yield_no2_tint/quantum_yield_no2_tint_t`
    ! class.
    !
    ! This test only checks the MPI functions currently. Additional test will
    ! be added as part of issue #177

    use musica_assert,                 only : assert, die
    use musica_constants,              only : dk => musica_dk
    use musica_mpi
    use musica_string,                 only : string_t
    use tuvx_quantum_yield,            only : quantum_yield_t
    use tuvx_quantum_yield_factory,    only : quantum_yield_type_name,        &
                                              quantum_yield_allocate
    use tuvx_test_utils,               only : check_values

    class(quantum_yield_t), pointer :: quantum_yield
    character, allocatable :: buffer(:)
    type(string_t) :: type_name
    integer :: pos, pack_size
    integer, parameter :: comm = MPI_COMM_WORLD

    real(dk) :: temperature1(3), temperature2(4)
    real(dk) :: deltaT1(2), deltaT2(3)
    real(dk) :: array1(3,2), array2(2,1)

    temperature1(:) = (/  12.5_dk,  16.3_dk, -290.4_dk /)
    deltaT1(:)      = (/ 412.3_dk, 0.032_dk /)
    array1(:,1)     = (/  42.3_dk, 132.4_dk,   13.4_dk /)
    array1(:,2)     = (/ 132.4_dk,  0.43_dk,   2.34_dk /)

    temperature2(:) = (/ -123.4_dk,  41.2_dk, 0.053_dk, 1.2e-7_dk /)
    deltaT2(:)      = (/ 132.45_dk,  13.4_dk, 1.324_dk /)
    array2(:,1)     = (/ 12.34_dk, -142.3_dk /)

    if( musica_mpi_rank( comm ) == 0 ) then
      allocate( quantum_yield_no2_tint_t :: quantum_yield )
      select type( quantum_yield )
      class is( quantum_yield_no2_tint_t )
        allocate( quantum_yield%parameters(2) )
        quantum_yield%parameters(1)%temperature = temperature1
        quantum_yield%parameters(1)%deltaT = deltaT1
        quantum_yield%parameters(1)%array = array1
        quantum_yield%parameters(2)%temperature = temperature2
        quantum_yield%parameters(2)%deltaT = deltaT2
        quantum_yield%parameters(2)%array = array2
      class default
        call die( 637746025 )
      end select
      type_name = quantum_yield_type_name( quantum_yield )
      pack_size = type_name%pack_size( comm ) + quantum_yield%pack_size( comm )
      allocate( buffer( pack_size ) )
      pos = 0
      call type_name%mpi_pack(     buffer, pos , comm )
      call quantum_yield%mpi_pack( buffer, pos , comm )
      call assert( 209765802, pos <= pack_size )
    end if

    call musica_mpi_bcast( pack_size , comm )
    if( musica_mpi_rank( comm ) .ne. 0 ) allocate( buffer( pack_size ) )
    call musica_mpi_bcast( buffer , comm )

    if( musica_mpi_rank( comm ) .ne. 0 ) then
      pos = 0
      call type_name%mpi_unpack( buffer, pos , comm )
      quantum_yield => quantum_yield_allocate( type_name )
      call quantum_yield%mpi_unpack( buffer, pos , comm )
      call assert( 264697883, pos <= pack_size )
    end if
    deallocate( buffer )

    call assert( 962406487, associated( quantum_yield ) )
    select type( quantum_yield )
    class is( quantum_yield_no2_tint_t )
      call assert( 964312021, allocated( quantum_yield%parameters) )
      call assert( 401267057, size( quantum_yield%parameters ) == 2 )
      associate( params => quantum_yield%parameters( 1 ) )
        call assert( 784078798, allocated( params%temperature ) )
        call assert( 160863043, allocated( params%deltaT ) )
        call assert( 891132836, allocated( params%array ) )
        call check_values( params%temperature, temperature1, 1.0e-6_dk )
        call check_values( params%deltaT,      deltaT1,      1.0e-6_dk )
        call check_values( params%array,       array1,       1.0e-6_dk )
      end associate
      associate( params => quantum_yield%parameters( 2 ) )
        call assert( 784078798, allocated( params%temperature ) )
        call assert( 104930018, allocated( params%deltaT ) )
        call assert( 891132836, allocated( params%array ) )
        call check_values( params%temperature, temperature2, 1.0e-6_dk )
        call check_values( params%deltaT,      deltaT2,      1.0e-6_dk )
        call check_values( params%array,       array2,       1.0e-6_dk )
      end associate
    class default
      call die( 464230348 )
    end select
    deallocate( quantum_yield )

  end subroutine test_quantum_yield_no2_tint_t

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program test_quantum_yield_no2_tint
