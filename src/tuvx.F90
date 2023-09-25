! Copyright (C) 2021 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
program tuvx
  ! Driver for stand-alone TUV-x and regression tests
  !
  ! If MPI support is compiled in, core driver constructs the core
  ! on the primary process, passes it to the other processes, and
  ! performs calculations on process 1. This allows MPI functions to
  ! be tested as part of the regression tests.

  use musica_assert,                   only : assert
  use musica_config,                   only : config_t
  use musica_constants,                only : dk => musica_dk
  use musica_string,                   only : string_t, to_char
  use musica_mpi
#ifdef MUSICA_USE_OPENMP
  use omp_lib
#endif
  use tuvx_core,                       only : core_t
  use tuvx_version,                    only: get_tuvx_version

  implicit none

#if MUSICA_USE_OPENMP
  type :: thread_t
    type(core_t), pointer :: core_ => null( )
  end type thread_t
  type(thread_t), allocatable :: threads(:)
#endif

  type(core_t), pointer :: core

  ! Command-line options
  character(len=*), parameter :: Iam = "TUV-x driver"
  character(len=256) :: argument
  type(string_t)     :: config_file_path
  character, allocatable :: buffer(:)
  integer :: pos, pack_size
  integer, parameter :: comm = MPI_COMM_WORLD
  integer :: i_thread
  logical :: found
  type(config_t) :: tuvx_config

  call musica_mpi_init( )

#ifdef MUSICA_USE_OPENMP
  write(*,*) "Running TUV-x(" // get_tuvx_version() // ") on ", omp_get_max_threads( ), " threads"
  allocate( threads( omp_get_max_threads( ) ) )
#else
  write(*,*) "Running TUV-x(" // get_tuvx_version() // ") without OpenMP support"
#endif

  ! Get the model configuration file and options from the command line
  if( command_argument_count() /= 1 ) then
    call fail_run( )
  endif
  call get_command_argument( 1, argument )

  config_file_path = argument

  ! instatiate and initialize photolysis core object on the
  ! primary MPI process
  if( musica_mpi_rank( comm ) == 0 ) then
    call tuvx_config%from_file( config_file_path%to_char( ) )
    core => core_t( config_file_path )
    pack_size = core%pack_size( comm ) + tuvx_config%pack_size( comm )
    allocate( buffer( pack_size ) )
    pos = 0
    call core%mpi_pack( buffer, pos, comm )
    call tuvx_config%mpi_pack( buffer, pos, comm )
    call assert( 115315474, pos <= pack_size )
  end if

  ! broadcast the core data to other MPI processes
  call musica_mpi_bcast( pack_size, comm )
  if( musica_mpi_rank( comm ) .ne. 0 ) allocate( buffer( pack_size ) )
  call musica_mpi_bcast( buffer, comm )

  ! Unpack the core on other MPI processes
  if( musica_mpi_rank( comm ) .ne. 0 ) then
    pos = 0
    allocate( core )
    call core%mpi_unpack( buffer, pos, comm )
    call tuvx_config%mpi_unpack( buffer, pos, comm )
    call assert( 501938283, pos <= pack_size )
  end if

#ifdef MUSICA_USE_OPENMP
  ! Create a core for each OpenMP thread
  do i_thread = 1, size( threads )
  associate( thread => threads( i_thread ) )
    pos = 0
    allocate( thread%core_ )
    call thread%core_%mpi_unpack( buffer, pos, comm )
    call assert( 736677220, pos <= pack_size )
  end associate
  end do
#endif

  ! Perform photolysis calculations on MPI process 1
  ! if MPI support is compiled in
  if( musica_mpi_size( comm ) == 1 .or. musica_mpi_rank( comm ) == 1 ) then
    call run_tuvx( )
  end if

  deallocate( buffer )
  deallocate( core )
#ifdef MUSICA_USE_OPENMP
  do i_thread = 1, size( threads )
  associate( thread => threads( i_thread ) )
    deallocate( thread%core_ )
  end associate
  end do
  deallocate( threads )
#endif

  call musica_mpi_finalize( )

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine run_tuvx( )
    ! Runs TUV-x under SZA conditions provided in the configuration

    use musica_assert,                 only : assert_msg
    use tuvx_grid,                     only : grid_t
    use tuvx_output,                   only : output_t
    use tuvx_profile,                  only : profile_t

    integer                   :: i_sza, i_thread, i_photo, i_level, i_dose
    class(grid_t),    pointer :: height
    class(grid_t),    pointer :: time               ! [hours]
    class(profile_t), pointer :: sza                ! [degrees]
    class(profile_t), pointer :: earth_sun_distance ! [AU]
    real(dk), allocatable     :: photo_rates(:,:,:) ! (time, vertical level, reaction) [s-1]
    real(dk), allocatable     :: dose_rates(:,:,:)  ! (time, vertical level, dose rate) [?]
    real(dk), allocatable     :: thread_photo_rates(:,:,:) ! (vertical level, reaction, thread) [s-1]
    real(dk), allocatable     :: thread_dose_rates(:,:,:)  ! (vertical level, dose rate, thread) [?]
    type(string_t)            :: file_path
    character(len=2)          :: diagnostic_label
    class(output_t), pointer  :: photo_output, dose_output
    type(config_t)            :: config

    height => core%get_grid( "height", "km" )
    time => core%get_grid( "time", "hours" )
    sza => core%get_profile( "solar zenith angle", "degrees" )
    earth_sun_distance => core%get_profile( "Earth-Sun distance", "AU" )

    call assert_msg( 891821582, sza%ncells_ == earth_sun_distance%ncells_,    &
                     "Solar zenith angle and Earth-Sun distance must be on "//&
                     "the same grid." )

    allocate( photo_rates( sza%ncells_ + 1,                                   &
                           height%ncells_ + 1,                                &
                           core%number_of_photolysis_reactions( ) ) )
    allocate( dose_rates(           sza%ncells_ + 1,                          &
                                    height%ncells_ + 1,                       &
                                    core%number_of_dose_rates( ) ) )

    ! set up output files
    nullify( photo_output )
    nullify( dose_output  )
    if( core%number_of_photolysis_reactions( ) > 0 ) then
      call config%empty( )
      call config%add( "file path", "photolysis_rate_constants.nc", Iam )
      call config%add( "include photolysis", .true., Iam )
      call config%add( "tuv-x configuration", tuvx_config, Iam )
      photo_output => output_t( config, core )
    end if
    if( core%number_of_dose_rates( ) > 0 ) then
      call config%empty( )
      call config%add( "file path", "dose_rates.nc", Iam )
      call config%add( "include dose rates", .true., Iam )
      call config%add( "tuv-x configuration", tuvx_config, Iam )
      dose_output => output_t( config, core )
    end if

    ! calculate photolysis and dose rates
    do i_sza = 1, sza%ncells_ + 1
      write(diagnostic_label,'(i2.2)') i_sza
      call core%run( sza%edge_val_( i_sza ),                                  &
                     earth_sun_distance%edge_val_( i_sza ),                   &
                     photolysis_rate_constants = photo_rates( i_sza, :, : ),  &
                     dose_rates = dose_rates( i_sza, :, : ),                  &
                     diagnostic_label = diagnostic_label )

      ! output results
      if( associated( photo_output ) ) then
        call photo_output%output( i_sza, core,                               &
            photolysis_rate_constants = photo_rates( i_sza, : , : ),         &
            time = time%edge_( i_sza ),                                      &
            solar_zenith_angle = sza%edge_val_( i_sza ),                     &
            earth_sun_distance = earth_sun_distance%edge_val_( i_sza ) )
      end if
      if( associated( dose_output ) ) then
        call dose_output%output( i_sza, core,                                &
            dose_rates= dose_rates( i_sza, : , : ),                          &
            time = time%edge_( i_sza ),                                      &
            solar_zenith_angle = sza%edge_val_( i_sza ),                     &
            earth_sun_distance = earth_sun_distance%edge_val_( i_sza ) )
      end if
    end do

    deallocate( height             )
    deallocate( time               )
    deallocate( sza                )
    deallocate( earth_sun_distance )
    if( associated( photo_output ) ) deallocate( photo_output )
    if( associated( dose_output  ) ) deallocate( dose_output  )

#if MUSICA_USE_OPENMP
    ! Compare results from threads for fixed solar zenith angle
    allocate( thread_photo_rates( size( photo_rates, 2 ),                     &
                                  size( photo_rates, 3 ),                     &
                                  omp_get_max_threads( ) ) )
    allocate( thread_dose_rates(  size( dose_rates,  2 ),                     &
                                  size( dose_rates,  3 ),                     &
                                  omp_get_max_threads( ) ) )
    !$omp parallel &
    !$omp   shared( threads, thread_photo_rates, thread_dose_rates )
    associate( thread => threads( omp_get_thread_num( ) + 1 ),                &
               photos => thread_photo_rates(:,:, omp_get_thread_num( ) + 1 ), &
               doses  => thread_dose_rates( :,:, omp_get_thread_num( ) + 1 ) )
      call thread%core_%run( 40.0_dk,                                         &
                             1.0_dk,                                          &
                             photolysis_rate_constants = photos,              &
                             dose_rates = doses )
    end associate
    !$omp end parallel

    do i_thread = 2, size( threads )
      do i_level = 1, size( thread_photo_rates, 1 )
        do i_photo = 1, size( thread_photo_rates, 2 )
          call assert_msg( 365420634,                                         &
                           thread_photo_rates( i_level, i_photo, i_thread )   &
                           .eq. thread_photo_rates( i_level, i_photo, 1 ),    &
                           "Thread result mismatch for thread "//             &
                           to_char( i_thread ) )
        end do
        do i_dose = 1, size( thread_dose_rates, 2 )
          call assert_msg( 401201021,                                         &
                           thread_dose_rates( i_level, i_dose, i_thread )     &
                           .eq. thread_dose_rates( i_level, i_dose, 1 ),      &
                           "Thread result mismatch for thread "//             &
                           to_char( i_thread ) )
        end do
      end do
    end do
#endif

  end subroutine run_tuvx

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine fail_run( )
    ! Fail run and print usage info

    write(*,*) "Usage: ./tuv-x configuration_file.json"
    stop 3

  end subroutine fail_run

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program tuvx
