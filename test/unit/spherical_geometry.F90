! Copyright (C) 2022 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
program test_spherical_geometry
  ! Tests components of the :f:type:`~tuvx_spherical_geometry` module

  use musica_mpi,                      only : musica_mpi_init,                &
                                              musica_mpi_finalize
  use tuvx_spherical_geometry

  implicit none

  call musica_mpi_init( )
  call test_spherical_geometry_t( )
  call musica_mpi_finalize( )

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine test_spherical_geometry_t( )
    ! Tests functionality of the :f:type:`~tuvx_spherical_geometry/spherical_geometry_t`
    ! class.
    !
    ! Currently only mpi functions are tested. Remaining tests will be added
    ! with issue #35

    use musica_assert,                 only : assert
    use musica_constants,              only : dk => musica_dk
    use musica_mpi
    use tuvx_test_utils,               only : check_values

    class(spherical_geometry_t), pointer :: calculator => null( )
    character, allocatable :: buffer(:)
    integer :: pos, pack_size
    integer, parameter :: comm = MPI_COMM_WORLD

    integer :: nid(0:2)
    real(dk) :: dsdh(0:2,3), sza

    nid = (/ 14, 53, 12 /)
    dsdh(:,1) = (/ 41.25_dk, 0.241_dk, -412.3_dk /)
    dsdh(:,2) = (/ 12.25_dk, 3.123_dk, -424.2_dk /)
    dsdh(:,3) = (/ 63.63_dk, 9.421_dk, -923.5_dk /)
    sza = 53.3_dk

    if( musica_mpi_rank( comm ) == 0 ) then
      allocate( calculator )
      calculator%nid_ = nid
      calculator%solar_zenith_angle_ = sza
      calculator%dsdh_ = dsdh
      pack_size = calculator%pack_size( comm )
      allocate( buffer( pack_size ) )
      pos = 0
      call calculator%mpi_pack( buffer, pos , comm )
      call assert( 181405192, pos <= pack_size )
    end if

    call musica_mpi_bcast( pack_size , comm )
    if( musica_mpi_rank( comm ) .ne. 0 ) allocate( buffer( pack_size ) )
    call musica_mpi_bcast( buffer , comm )

    if( musica_mpi_rank( comm ) .ne. 0 ) then
      pos = 0
      allocate( calculator )
      call calculator%mpi_unpack( buffer, pos , comm )
      call assert( 174687646, pos <= pack_size )
    end if
    deallocate( buffer )

    call assert( 846692182, associated( calculator ) )
    call assert( 731015064, allocated( calculator%nid_ ) )
    call assert( 108226007, allocated( calculator%dsdh_ ) )
    call assert( 615337946, calculator%solar_zenith_angle_ == sza )
    call assert( 885117521, size( calculator%nid_ ) == size( nid ) )
    call assert( 830637735, calculator%nid_(0) == nid(0) )
    call assert( 992171559, calculator%nid_(1) == nid(1) )
    call assert( 539539406, calculator%nid_(2) == nid(2) )
    call check_values( 473236945, calculator%dsdh_, dsdh, 1.0e-6_dk )

    deallocate( calculator )

  end subroutine test_spherical_geometry_t

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program test_spherical_geometry
