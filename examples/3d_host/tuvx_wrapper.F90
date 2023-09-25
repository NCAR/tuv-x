! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
module tuvx_wrapper
  ! Wrapper for TUV-x functionality
  !
  ! This module handles parallel TUV-x calculations and conversion
  ! to and from host model data structures

  use tuvx_core,                       only : core_t
  use tuvx_grid_updater,               only : grid_updater_t
  use tuvx_profile_updater,            only : profile_updater_t

  implicit none
  private

  public :: tuvx_init, tuvx_run, tuvx_finalize

  type :: tuvx_wrapper_t
    private
    type(core_t),             pointer :: core_
    class(grid_updater_t),    pointer :: height_
    class(profile_updater_t), pointer :: air_
    class(profile_updater_t), pointer :: temperature_
  contains
    final :: finalize
  end type tuvx_wrapper_t

  integer :: tuvx_comm                                  ! MPI communicator
  type(tuvx_wrapper_t), allocatable :: tuvx_wrappers(:) ! wrappers ( OMP thread )

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine tuvx_init( mpi_comm, omp_threads, tuvx_config_path,              &
      n_vertical_levels, n_photolysis_reactions )
    ! Initializes TUV-x cores for each process/thread and sets up grids
    ! and profiles that will be updated based on 3D model state data at
    ! run time

    use musica_mpi
    use musica_string,                 only : string_t
    use tuvx_grid_from_host,           only : grid_from_host_t
    use tuvx_grid_warehouse,           only : grid_warehouse_t
    use tuvx_profile_from_host,        only : profile_from_host_t
    use tuvx_profile_warehouse,        only : profile_warehouse_t

    integer,          intent(in) :: mpi_comm           ! MPI communicator for processes that will call tuvx_run( )
    integer,          intent(in) :: omp_threads        ! Number of OMP threads that will call tuvx_run( )
    character(len=*), intent(in) :: tuvx_config_path   ! Path to the TUV-x configuration data file
    integer,          intent(in) :: n_vertical_levels  ! Number of vertical levels per column
    integer,          intent(in) :: n_photolysis_reactions  ! Number of photolysis reactions

    class(core_t), pointer :: core
    character, allocatable :: buffer(:)
    integer :: pack_size, pos
    type(grid_from_host_t),    pointer :: height
    type(profile_from_host_t), pointer :: temperature, air
    type(grid_warehouse_t),    pointer :: grids
    type(profile_warehouse_t), pointer :: profiles
    type(string_t), allocatable :: photo_labels(:)

    ! allocate a wrapper for each OMP thread on every MPI process
    allocate( tuvx_wrappers( omp_threads ) )

    ! set up the grids that will be updated at each time step
    height => grid_from_host_t( "height", "km", n_vertical_levels )
    grids => grid_warehouse_t( )
    call grids%add( height )

    ! set up the profiles that will be updated at each time step
    air => profile_from_host_t( "air", "molecule cm-3", n_vertical_levels )
    temperature => profile_from_host_t( "temperature", "K", n_vertical_levels )
    profiles => profile_warehouse_t( )
    call profiles%add( air )
    call profiles%add( temperature )

    ! save the MPI communicator for use at run time
    tuvx_comm = mpi_comm

    ! pack the core on the primary MPI process
    if( musica_mpi_rank( ) == 0 ) then
      core => core_t( tuvx_config_path, grids, profiles )

      ! this could be used to dynamically set the number of photolysis
      ! reaction rate constants for the 3D model and map to chemistry
      ! packages
      photo_labels = core%photolysis_reaction_labels( )
      call assert_msg( 287966507, size( photo_labels ) .eq.                   &
                                  n_photolysis_reactions,                     &
                                  "Photolysis reaction mismatch" )
      pack_size = core%pack_size( tuvx_comm )
      allocate( buffer( pack_size ) )
      pos = 0
      call core%mpi_pack( buffer, pos, tuvx_comm )
      deallocate( core )
    end if

    ! broadcast the core data to all MPI processes
    call musica_mpi_bcast( pack_size )
    if( musica_mpi_rank( ) .ne. 0 ) allocate( buffer( pack_size ) )
    call musica_mpi_bcast( buffer )

    ! unpack the core for each OMP thread on every MPI process
    allocate( tuvx_wrappers )
    do i_thread = 1, omp_threads
    associate( wrapper => tuvx_wrappers( i_thread ) )
      allocate( wrapper%core_ )
      pos = 0
      call wrapper%core_%mpi_unpack( buffer, pos, tuvx_comm )
      wrapper%height_      => wrapper%core_%get_updater( height      )
      wrapper%temperature_ => wrapper%core_%get_updater( temperature )
      wrapper%air_         => wrapper%core_%get_updater( air         )
    end associate
    end do

    deallocate( grids    )
    deallocate( profiles )

  end subroutine tuvx_init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine tuvx_run( omp_thread, height, air_density, temperature,          &
      photolysis_rate_constants )
    ! Calculates photolysis rate constants for the current model conditions

    real(kind=dk), intent(in)  :: height(:,:) ! height above sea level [km] (layer, column)
    real(kind=dk), intent(in)  :: air_density(:,:) ! number density of dry (?) air [molecule cm-3] (layer, column)
    real(kind=dk), intent(in)  :: temperature(:,:) ! temperature [K] (layer, column)
    real(kind=dk), intent(out) :: photolysis_rate_constants(:,:,:) ! photolysis rate constants [s-1] (reaction, layer, column)

    integer :: i_col

    associate( wrapper => tuvx_wrappers( omp_thread ) )
      do i_col = 1, size( height, 2 )
        call wrapper%height_%update( height( :, i_col ) )
        call wrapper%air_%update( air_density( :, i_col ) )
        call wrapper%temperature_%update( temperature( :, i_col ) )
        call core%calculate( photolysis_rate_constants                        &
                             = photolysis_rate_constants( :, :, i_col ) )
      end do
    end associate

  end subroutine tuvx_run

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine tuvx_finalize( )
    ! Cleans up memory associated with TUV-x

    deallocate( tuvx_wrappers )

  end subroutine tuvx_finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine finalize( this )
    ! Cleans up memory associated with a TUV-x wrapper

    type(tuvx_wrapper_t), intent(inout) :: this

    if( associated( this%core_        ) ) deallocate( this%core_        )
    if( associated( this%height_      ) ) deallocate( this%height_      )
    if( associated( this%air_         ) ) deallocate( this%air_         )
    if( associated( this%temperature_ ) ) deallocate( this%temperature_ )

  end subroutine finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_wrapper
