! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
program atmos_model
  ! Mock 3D atmospheric model with example interface to TUV-x

  use tuvx_wrapper,                    only : tuvx_init, tuvx_run,            &
                                              tuvx_finalize
  use omp_lib
  use mpi

  implicit none

  integer, parameter :: kColumns = 16              ! Columns per MPI process
  integer, parameter :: kVerticalLevels = 24       ! Number of vertical levels per column
  integer, parameter :: kPhotolysisReactions = 4   ! Number of photolysis reactions
  integer, parameter :: kTimeSteps = 10            ! Model time steps
  integer, parameter :: kDouble = kind(0.0d0)      ! Double precision floating point kind
  integer, parameter :: mpi_comm = MPI_COMM_WORLD  ! MPI communicator

  real(kDouble) :: height(      kVerticalLevels, kColumns ) ! height above sea level [km]
  real(kDouble) :: air_density( kVerticalLevels, kColumns ) ! number density of dry (?) air [molecule cm-3]
  real(kDouble) :: temperature( kVerticalLevels, kColumns ) ! temperature [K]
  real(kDouble) :: photolysis_rate_constants( kPhotolysisReactions,           &
                                              kVerticalLevels,                &
                                              kColumns ) ! rate constants [s-1]

  character(len=:), allocatable :: tuvx_config_path
  integer :: omp_threads
  integer :: ierr

  ! initialize MPI, OpenMP
  call mpi_init( ierr )
  call check_mpi_status( ierr )
  omp_threads = omp_get_max_threads( )

  ! initialize model
  call model_init( )

  do i_time = 1, kTimeSteps
    call model_run( )
  end do

  call model_finalize( )
  call mpi_finalize( ierr )
  call check_mpi_status( ierr )

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine model_init( )
    ! Initialize the model data structures

    tuvx_config_path = "path/to/tuvx/config.json"
    call tuvx_init( mpi_comm, omp_threads, tuvx_config_path, kVerticalLevels, &
                   kPhotolysisReactions )

  end subroutine model_init( )

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine model_run( )
    ! Advance the model state in time, including calculation of photolysis
    ! rate constants and dose rates

    integer :: start_col, end_col
    call update_model_state( )

    !$omp parallel &
    !$omp shared ( height, air_density, temperature, photolysis_rate_constants ) &
    !$omp private ( start_col, end_col )

    ! this should modified for cases where the number of columns is not a multiple
    ! of the number for OMP threads
    start_col = kColumns / omp_get_max_threads( ) * omp_get_thread_num( )
    end_col = kColumns / omp_get_max_threads( ) * ( omp_get_thread_num( ) + 1 )
    call tuvx_run( omp_get_thread_num( ),                                     &
                   height( :, start_col : end_col ),                          &
                   air_density( :, start_col : end_col ),                     &
                   temperature( :, start_col : end_col ),                     &
                   photolysis_rate_constants( :, :, start_col : end_col ) )
    !$omp end parallel

    ! here you might merge / output the model state, etc.

  end subroutine model_run

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine model_finalize( )
    ! Clean up memory

    call tuvx_finalize( )

  end subroutine model_finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine update_model_state( )
    ! Randomizes the model state as a stand-in for advancing the state in a
    ! real 3D model

    integer :: i_cell, i_col
    real(kDouble) :: random

    do i_col = 1, kColumns
      do i_cell = 1, kVerticalLevels
        call random_number( random )
        height(      i_cell, i_col ) = ( i_cell - 1 ) * 10.0d0 + 1.0d0 * random
        air_density( i_cell, i_col ) = ( kVerticalLevels - i_cell - 1 )       &
                                       * 2.54d19 + 100.0d0 * random
        temperature( i_cell, i_col ) = 250.0d0 + 100.0d0 * random
      end do
    end do

  end subroutine update_model_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine check_mpi_status( ierr )

    integer, intent(in) :: ierr

    if( ierr /= MPI_SUCCESS ) then
      write(*,*) "MPI error: ", ierr
      stop 3
    end if

  end subroutine check_mpi_status

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program atmos_model
