! Copyright (C) 2021 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
program test_profile_from_host

  use musica_mpi,                      only : musica_mpi_init,                &
                                              musica_mpi_finalize
  use tuvx_test_utils,                 only : check_values

  implicit none

  call musica_mpi_init( )
  call test_profile_from_host_t( )
  call musica_mpi_finalize( )

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine test_profile_from_host_t( )

    use musica_assert,                 only : assert, almost_equal, die
    use musica_constants,              only : dk => musica_dk
    use musica_string,                 only : string_t
    use musica_mpi
    use tuvx_profile,                  only : profile_t
    use tuvx_profile_from_host,        only : profile_from_host_t,            &
                                              profile_updater_t
    use tuvx_profile_factory,          only : profile_type_name,              &
                                              profile_allocate

    class(profile_t), pointer :: my_profile
    type(profile_updater_t) :: my_updater
    character, allocatable :: buffer(:)
    integer :: pos, pack_size
    type(string_t) :: type_name
    integer, parameter :: comm = MPI_COMM_WORLD
    real(kind=dk), parameter :: tol = 1.0e-10_dk
    real(kind=dk) :: edges(4)
    real(kind=dk) :: mids(3)
    real(kind=dk) :: deltas(3)
    real(kind=dk) :: dens(3)
    real(kind=dk) :: exos(4)
    real(kind=dk) :: burden(3)

    if( musica_mpi_rank( comm ) == 0 ) then
      my_profile => profile_from_host_t( "foo", "bars", 3 )
      type_name = profile_type_name( my_profile )
      pack_size = type_name%pack_size( comm ) + my_profile%pack_size( comm )
      allocate( buffer( pack_size ) )
      pos = 0
      call type_name%mpi_pack( buffer, pos, comm )
      call my_profile%mpi_pack( buffer, pos, comm )
      call assert( 102131176, pos <= pack_size )
    end if

    call musica_mpi_bcast( pack_size, comm )
    if( musica_mpi_rank( comm ) .ne. 0 ) allocate( buffer( pack_size ) )
    call musica_mpi_bcast( buffer, comm )

    if( musica_mpi_rank( comm ) .ne. 0 ) then
      pos = 0
      call type_name%mpi_unpack( buffer, pos, comm )
      my_profile => profile_allocate( type_name )
      call my_profile%mpi_unpack( buffer, pos, comm )
      call assert( 997319163, pos <= pack_size )
    end if
    deallocate( buffer )

    call assert( 204373202, my_profile%handle_ == "foo" )
    call assert( 941386138, my_profile%units( ) == "bars" )
    call assert( 995865924, size( my_profile%mid_val_    ) == 3 )
    call assert( 992507151, size( my_profile%edge_val_   ) == 4 )
    call assert( 199561190, size( my_profile%delta_val_  ) == 3 )
    call assert( 371623628, size( my_profile%layer_dens_ ) == 3 )

    select type( my_profile )
    class is( profile_from_host_t )
      my_updater = profile_updater_t( my_profile )
    class default
      call die( 870112487 )
    end select

    ! specify edges only
    edges  = (/ 1.0_dk, 2.0_dk, 4.0_dk, 10.0_dk /)
    ! calculated
    mids   = (/ 1.5_dk, 3.0_dk, 7.0_dk /)
    deltas = (/ 1.0_dk, 2.0_dk, 6.0_dk /)
    dens   = (/ 0.0_dk, 0.0_dk, 0.0_dk /)
    exos   = (/ 0.0_dk, 0.0_dk, 0.0_dk, 0.0_dk /)
    burden = (/ 0.0_dk, 0.0_dk, 0.0_dk /)
    call my_updater%update( edge_values = edges )
    call check_values( 275069338, my_profile%mid_val_,          mids, tol )
    call check_values( 169920834, my_profile%edge_val_,        edges, tol )
    call check_values( 617288680, my_profile%delta_val_,      deltas, tol )
    call check_values( 447131776, my_profile%layer_dens_,       dens, tol )
    call check_values( 612024373, my_profile%exo_layer_dens_,   exos, tol )
    call check_values( 159392220, my_profile%burden_dens_,    burden, tol )

    ! specify edges, mids, dens
    edges  = (/ 0.5_dk, 9.8_dk, 15.4_dk, 45.0_dk /)
    mids   = (/ 1.0_dk, 12.3_dk, 32.4_dk /)
    dens   = (/ 94.3_dk, 0.52_dk, 12.3_dk /)
    ! calculated
    deltas = (/ 9.8_dk - 0.5_dk, 15.4_dk - 9.8_dk, 45.0_dk - 15.4_dk /)
    exos   = (/ 94.3_dk, 0.52_dk, 12.3_dk, 0.0_dk /)
    burden = (/ 94.3_dk + 0.52_dk + 12.3_dk, 0.52_dk + 12.3_dk, 12.3_dk /)
    call my_updater%update( mid_point_values = mids,                          &
                            edge_values      = edges,                         &
                            layer_densities  = dens )
    call check_values( 792121290, my_profile%mid_val_,          mids, tol )
    call check_values( 399233230, my_profile%edge_val_,        edges, tol )
    call check_values( 106229367, my_profile%delta_val_,      deltas, tol )
    call check_values( 613341306, my_profile%layer_dens_,       dens, tol )
    call check_values( 466544996, my_profile%exo_layer_dens_,   exos, tol )
    call check_values( 796782485, my_profile%burden_dens_,    burden, tol )

    ! specify edges, dens, scale height
    edges  = (/ 1.0_dk, 2.0_dk, 4.0_dk, 10.0_dk /)
    dens   = (/ 10.0_dk, 2.5_dk, 5.0_dk /)
    ! calculated
    mids   = (/ 1.5_dk, 3.0_dk, 7.0_dk /)
    deltas = (/ 1.0_dk, 2.0_dk, 6.0_dk /)
    exos   = (/ 10.0_dk, 2.5_dk, 5.0_dk, 10.0_dk * 2.0_dk * 1.0e5_dk /)
    burden = (/ 10.0_dk + 2.5_dk + 5.0_dk, 2.5_dk + 5.0_dk, 5.0_dk /)
    call my_updater%update( edge_values = edges,                             &
                            layer_densities = dens,                          &
                            scale_height = 2.0_dk )
    dens(3) = 5.0_dk + 10.0_dk * 2.0_dk * 1.0e5_dk
    burden(:) = burden(:) + 10.0_dk * 2.0_dk * 1.0e5_dk
    call check_values( 420261592, my_profile%mid_val_,          mids, tol )
    call check_values( 585154189, my_profile%edge_val_,        edges, tol )
    call check_values( 132522036, my_profile%delta_val_,      deltas, tol )
    call check_values( 927373531, my_profile%layer_dens_,       dens, tol )
    call check_values( 192266129, my_profile%exo_layer_dens_,   exos, tol )
    call check_values( 639633975, my_profile%burden_dens_,    burden, tol )

    ! specify edges, mids, dens, exo density
    edges  = (/ 0.5_dk, 9.8_dk, 15.4_dk, 45.0_dk /)
    mids   = (/ 1.0_dk, 12.3_dk, 32.4_dk /)
    dens   = (/ 94.3_dk, 0.52_dk, 12.3_dk /)
    ! calculated
    deltas = (/ 9.8_dk - 0.5_dk, 15.4_dk - 9.8_dk, 45.0_dk - 15.4_dk /)
    exos   = (/ 94.3_dk, 0.52_dk, 12.3_dk, 12.0_dk /)
    burden = (/ 94.3_dk + 0.52_dk + 12.3_dk, 0.52_dk + 12.3_dk, 12.3_dk /)
    call my_updater%update( mid_point_values = mids,                          &
                            edge_values      = edges,                         &
                            layer_densities  = dens,                          &
                            exo_density      = 12.0_dk )
    dens(3) = dens(3) + 12.0_dk
    burden(:) = burden(:) + 12.0_dk
    call check_values( 234764062, my_profile%mid_val_,          mids, tol )
    call check_values( 747140308, my_profile%edge_val_,        edges, tol )
    call check_values( 294508155, my_profile%delta_val_,      deltas, tol )
    call check_values( 124351251, my_profile%layer_dens_,       dens, tol )
    call check_values( 571719097, my_profile%exo_layer_dens_,   exos, tol )
    call check_values( 184095344, my_profile%burden_dens_,    burden, tol )

    deallocate( my_profile )

    my_profile => profile_from_host_t( "baz", "qux", 0 )
    select type( my_profile )
    class is( profile_from_host_t )
      my_updater = profile_updater_t( my_profile )
    class default
      call die( 174774784 )
    end select
    call my_updater%update( edge_values = (/ 12.0_dk /) )
    call check_values( 918052971, my_profile%edge_val_, (/ 12.0_dk /), tol )
    call assert( 281541699, size( my_profile%mid_val_        ) .eq. 0 )
    call assert( 387142498, size( my_profile%delta_val_      ) .eq. 0 )
    call assert( 216985594, size( my_profile%layer_dens_     ) .eq. 0 )
    call check_values( 946828689, my_profile%exo_layer_dens_, (/ 0.0_dk /),   &
                       tol )
    call assert( 776671785, size( my_profile%burden_dens_    ) .eq. 0 )

    deallocate( my_profile )

  end subroutine test_profile_from_host_t

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program test_profile_from_host
