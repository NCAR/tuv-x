! Copyright (C) 2022 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
program test_radiator_from_host

  use musica_mpi,                      only : musica_mpi_init,                &
                                              musica_mpi_finalize
  use tuvx_test_utils,                 only : check_values

  implicit none

  call musica_mpi_init( )
  call test_radiator_from_host_t( )
  call musica_mpi_finalize( )

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine test_radiator_from_host_t( )

    use musica_assert,                 only : assert, almost_equal, die
    use musica_constants,              only : dk => musica_dk
    use musica_mpi
    use musica_string,                 only : string_t
    use tuvx_cross_section_warehouse,  only : cross_section_warehouse_t
    use tuvx_grid,                     only : grid_t
    use tuvx_grid_from_host,           only : grid_from_host_t
    use tuvx_grid_warehouse,           only : grid_warehouse_t
    use tuvx_profile_warehouse,        only : profile_warehouse_t
    use tuvx_radiator,                 only : radiator_t
    use tuvx_radiator_from_host,       only : radiator_from_host_t,           &
                                              radiator_updater_t
    use tuvx_radiator_factory,         only : radiator_type_name,             &
                                              radiator_allocate

    class(radiator_t), pointer :: radiator
    type(radiator_updater_t) :: radiator_updater
    class(grid_t), pointer :: height, wavelength
    type(grid_warehouse_t) :: grids
    type(profile_warehouse_t) :: profiles
    type(cross_section_warehouse_t) :: cross_sections
    character, allocatable :: buffer(:)
    integer :: pos, pack_size
    type(string_t) :: type_name
    integer, parameter :: comm = MPI_COMM_WORLD

    real(kind=dk), parameter :: tol = 1.0e-10_dk
    real(kind=dk) :: od(3,2)
    real(kind=dk) :: ssa(3,2)
    real(kind=dk) :: asym(3,2)

    if( musica_mpi_rank( comm ) == 0 ) then
      height => grid_from_host_t( "height", "km", 3 )
      wavelength => grid_from_host_t( "wavelength", "nm", 2 )
      radiator => radiator_from_host_t( "foo", height, wavelength )
      deallocate( height )
      deallocate( wavelength )
      type_name = radiator_type_name( radiator )
      pack_size = type_name%pack_size( comm ) + radiator%pack_size( comm )
      allocate( buffer( pack_size ) )
      pos = 0
      call type_name%mpi_pack( buffer, pos, comm )
      call radiator%mpi_pack(  buffer, pos, comm )
      call assert( 134439177, pos <= pack_size )
    end if

    call musica_mpi_bcast( pack_size, comm )
    if( musica_mpi_rank( comm ) .ne. 0 ) allocate( buffer( pack_size ) )
    call musica_mpi_bcast( buffer, comm )

    if( musica_mpi_rank( comm ) .ne. 0 ) then
      pos = 0
      call type_name%mpi_unpack( buffer, pos, comm )
      radiator => radiator_allocate( type_name )
      call radiator%mpi_unpack( buffer, pos, comm )
      call assert( 202578394, pos <= pack_size )
    end if
    deallocate( buffer )

    call assert( 769434426, radiator%handle_ == "foo" )
    call assert( 876488464, size( radiator%state_%layer_OD_,  1 ) == 3 )
    call assert( 310084727, size( radiator%state_%layer_OD_,  2 ) == 2 )
    call assert( 704878321, size( radiator%state_%layer_SSA_, 1 ) == 3 )
    call assert( 534721417, size( radiator%state_%layer_SSA_, 2 ) == 2 )
    call assert( 699614014, size( radiator%state_%layer_G_,   1 ) == 3 )
    call assert( 311990261, size( radiator%state_%layer_G_,   2 ) == 2 )

    select type( radiator )
    type is( radiator_from_host_t )
      radiator_updater = radiator_updater_t( radiator )
    class default
      call die( 298102874 )
    end select

    ! specify optical depths only
    od(:,1) = (/ 12.5_dk, 42.3_dk,  0.4_dk /)
    od(:,2) = (/ 49.2_dk, 12.5_dk, 92.1_dk /)
    ! calculated
    ssa(:,:)  = 0.0_dk
    asym(:,:) = 0.0_dk
    call radiator_updater%update( optical_depths = od )
    call check_values( 109698866, radiator%state_%layer_OD_,         od, tol )
    call check_values( 164178652, radiator%state_%layer_SSA_,       ssa, tol )
    call check_values( 388815342, radiator%state_%layer_G_(:,:,1), asym, tol )

    ! specify all optical properties
    od(:,1)   = (/ 0.32_dk, -2.3_dk,  9.2_dk /)
    od(:,2)   = (/ 61.2_dk, 8.32_dk, 0.42_dk /)
    ssa(:,1)  = (/ 12.3_dk, 0.23_dk, 12.4_dk /)
    ssa(:,2)  = (/ 92.3_dk, 12.3_dk, 31.2_dk /)
    asym(:,1) = (/ 9.53_dk, 62.3_dk,  6.4_dk /)
    asym(:,2) = (/  3.4_dk, 82.3_dk, 9.52_dk /)
    call radiator_updater%update( optical_depths = od,                        &
                                  single_scattering_albedos = ssa,            &
                                  asymmetry_factors = asym )
    call check_values( 712335285, radiator%state_%layer_OD_,         od, tol )
    call check_values( 877227882, radiator%state_%layer_SSA_,       ssa, tol )
    call check_values( 707070978, radiator%state_%layer_G_(:,:,1), asym, tol )

    ! update state should do nothing
    call radiator%update_state( grids, profiles, cross_sections )

    deallocate( radiator )

  end subroutine test_radiator_from_host_t

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program test_radiator_from_host
