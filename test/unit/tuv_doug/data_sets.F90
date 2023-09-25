! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> Tests for the H2O -> H + OH rate constant data sets
program doug_data_set

  use musica_mpi,                      only : musica_mpi_init,                &
                                              musica_mpi_finalize

  implicit none

  call musica_mpi_init( )
  call test_data_set( )
  call musica_mpi_finalize( )

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine test_data_set( )

    use musica_assert,                 only : assert
    use musica_config,                 only : config_t
    use musica_constants,              only : dk => musica_dk
    use musica_iterator,               only : iterator_t
    use musica_mpi
    use musica_string,                 only : string_t
    use tuv_doug,                      only : get_grids, calculate, wc, wl
    use tuvx_cross_section,            only : cross_section_t
    use tuvx_cross_section_factory,    only : cross_section_builder,          &
                                              cross_section_type_name,        &
                                              cross_section_allocate
    use tuvx_grid,                     only : grid_t
    use tuvx_grid_warehouse,           only : grid_warehouse_t
    use tuvx_profile,                  only : profile_t
    use tuvx_profile_warehouse,        only : profile_warehouse_t
    use tuvx_quantum_yield,            only : quantum_yield_t
    use tuvx_quantum_yield_factory,    only : quantum_yield_type_name,        &
                                              quantum_yield_allocate
    use tuvx_test_utils,               only : check_values

    class(grid_warehouse_t),    pointer :: grids
    class(profile_warehouse_t), pointer :: profiles
    class(cross_section_t),     pointer :: cross_section
    class(quantum_yield_t),     pointer :: quantum_yield

    character(len=*), parameter :: Iam = "H2O cross section test"
    type(config_t) :: config, config_pair, cs_config, qy_config
    class(iterator_t), pointer :: iter
    type(string_t) :: cs_type_name, qy_type_name, label
    character, allocatable :: buffer(:)
    integer :: pos, pack_size
    integer, parameter :: comm = MPI_COMM_WORLD
    real(kind=dk), allocatable :: cross_section_data(:,:)
    real(kind=dk), allocatable :: quantum_yield_data(:,:)
    real(kind=dk), allocatable :: tuvx_xsqy(:,:)
    real, allocatable :: doug_xsqy(:,:)
    class(profile_t), pointer :: air, temperature
    class(grid_t), pointer :: wavelength
    real(kind=dk) :: tolerance
    integer :: i

    ! Load grids based on Doug's TUV
    grids => get_grids( )

    ! Load test profiles
    call config%from_file( "test/data/profile.doug.config.json" )
    profiles    => profile_warehouse_t( config, grids )
    air         => profiles%get_profile( "air", "molecule cm-3" )
    temperature => profiles%get_profile( "temperature", "K" )

    ! Load cross section and quantum yield configurations for Doug's
    ! original TUV code
    call config%from_file( "test/data/xsqy.doug.config.json" )
    iter => config%get_iterator( )

    ! Check each cross section/quantum yield configuration and data
    do while( iter%next( ) )

      call config%get( iter, config_pair, Iam )
      call config_pair%get( "cross section", cs_config, Iam )
      call config_pair%get( "quantum yield", qy_config, Iam )
      call config_pair%get( "label",         label,     Iam )
      call config_pair%get( "tolerance",     tolerance, Iam,                  &
                            default = 1.0e-6_dk )

      ! Load and test cross section
      if( musica_mpi_rank( comm ) == 0 ) then
        cross_section => cross_section_builder( cs_config, grids, profiles )
        cs_type_name = cross_section_type_name( cross_section )
        quantum_yield => quantum_yield_t( qy_config, grids, profiles )
        qy_type_name = quantum_yield_type_name( quantum_yield )
        pack_size = cs_type_name%pack_size( comm ) +                          &
                    cross_section%pack_size( comm ) +                         &
                    qy_type_name%pack_size( comm ) +                          &
                    quantum_yield%pack_size( comm )
        allocate( buffer( pack_size ) )
        pos = 0
        call cs_type_name%mpi_pack(  buffer, pos, comm )
        call cross_section%mpi_pack( buffer, pos, comm )
        call qy_type_name%mpi_pack(  buffer, pos, comm )
        call quantum_yield%mpi_pack( buffer, pos, comm )
        call assert( 108537023, pos <= pack_size )
      end if

      call musica_mpi_bcast( pack_size, comm )
      if( musica_mpi_rank( comm ) .ne. 0 ) allocate( buffer( pack_size ) )
      call musica_mpi_bcast( buffer, comm )

      if( musica_mpi_rank( comm ) .ne. 0 ) then
        pos = 0
        call cs_type_name%mpi_unpack( buffer, pos, comm )
        cross_section => cross_section_allocate( cs_type_name )
        call cross_section%mpi_unpack( buffer, pos, comm )
        call qy_type_name%mpi_unpack( buffer, pos, comm )
        quantum_yield => quantum_yield_allocate( qy_type_name )
        call quantum_yield%mpi_unpack( buffer, pos, comm )
        call assert( 598402802, pos <= pack_size )
      end if
      deallocate( buffer )

      cross_section_data = cross_section%calculate( grids, profiles )
      quantum_yield_data = quantum_yield%calculate( grids, profiles )
      allocate( tuvx_xsqy, mold = cross_section_data )
      tuvx_xsqy(:,:) = cross_section_data(:,:) * quantum_yield_data(:,:)

      call calculate( label%val_,                                             &
                      real( temperature%edge_val_(:temperature%ncells_+1) ),  &
                      real( air%mid_val_ ), doug_xsqy )

      wavelength => grids%get_grid( "wavelength", "nm" )
      write(*,*) label%val_
      do i = 1, size( tuvx_xsqy, dim=2 )
        write(*,*) i, wavelength%edge_(i), wavelength%mid_(i),                &
                   cross_section_data(62,i),             &
                   quantum_yield_data(62,i), tuvx_xsqy(62,i), wl(i),          &
                   real( doug_xsqy(62,i), kind=dk )
      end do
      write(*,*) size( tuvx_xsqy, dim=2 ) + 1,                                &
                 wavelength%edge_(wavelength%ncells_+1)
      deallocate( wavelength )

      ! Skip first two bins because Lyman-Alpha bins are different in
      ! Doug's version of TUV-x. Data sets were adapted to have Lyman-Alpha
      ! specific data go into the TUV-x Lyman-Alpha bin 121.4-121.9 nm
      call check_values( 377150482, tuvx_xsqy(:,3:),                          &
                         real( doug_xsqy(:,3:), kind=dk ), tolerance )

      deallocate( cross_section      )
      deallocate( quantum_yield      )
      deallocate( cross_section_data )
      deallocate( quantum_yield_data )
      deallocate( tuvx_xsqy          )

    end do

    deallocate( iter        )
    deallocate( air         )
    deallocate( temperature )
    deallocate( grids       )
    deallocate( profiles    )

  end subroutine test_data_set

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program doug_data_set

