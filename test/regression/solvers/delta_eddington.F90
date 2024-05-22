! Copyright (C) 2023-2024 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
program test_cpp_delta_eddington_solver

  use, intrinsic :: iso_c_binding
  use musica_constants,                only : dk => musica_dk
  use musica_mpi,                      only : musica_mpi_init,                &
                                              musica_mpi_finalize
  use musica_string,                   only : string_t

  implicit none

  ! Struct to hold input conditions for the C++ Delta-Eddington solver
  type, bind(c) :: solver_input_t_c
    integer(c_int) :: n_wavelengths_
    integer(c_int) :: n_levels_
    integer(c_int) :: n_columns_
    type(c_ptr) :: solar_zenith_angles_   ! (columns)
    type(c_ptr) :: earth_sun_distances_   ! (columns)
    type(c_ptr) :: altitude_mid_points_   ! (columns, levels)
    type(c_ptr) :: altitude_edges_        ! (columns, levels+1)
    type(c_ptr) :: wavelength_mid_points_ ! (wavelengths)
    type(c_ptr) :: wavelength_edges_      ! (wavelengths+1)
  end type solver_input_t_c

  ! Struct to hold output radiation fields from the C++ Delta-Eddington solver
  type, bind(c) :: solver_output_t_c
    integer(c_int) :: n_wavelengths_
    integer(c_int) :: n_levels_
    integer(c_int) :: n_columns_
    type(c_ptr) :: flux_direct_  ! (columns, levels, wavelengths)
    type(c_ptr) :: flux_up_      ! (columns, levels, wavelengths)
    type(c_ptr) :: flux_down_    ! (columns, levels, wavelengths)
    type(c_ptr) :: irrad_direct_ ! (columns, levels, wavelengths)
    type(c_ptr) :: irrad_up_     ! (columns, levels, wavelengths)
    type(c_ptr) :: irrad_down_   ! (columns, levels, wavelengths)
  end type solver_output_t_c

  ! Function to calculate the radiation fields using the C++ Delta-Eddington solver
  interface
    function run_delta_eddington_solver_c( input ) bind(c, name='RunDeltaEddingtonSolver')
      import :: solver_input_t_c
      import :: solver_output_t_c
      type(solver_input_t_c), value :: input
      type(solver_output_t_c) run_delta_eddington_solver_c
    end function run_delta_eddington_solver_c
    subroutine free_output_c( output ) bind(c, name='FreeOutput')
      import :: solver_output_t_c
      type(solver_output_t_c), value :: output
    end subroutine free_output_c
  end interface

  type(string_t) :: config_file_path

  call musica_mpi_init( )

  ! Run the TUV 5.4 test
  config_file_path = 'examples/tuv_5_4.json'
  call test_cpp_delta_eddington_solver_t(config_file_path)

  call musica_mpi_finalize( )

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Runs the Fortran-based Delta-Eddington solver under prescribed conditions
  ! and then runs the C++-based Delta-Eddington solver under the same
  ! conditions, comparing the results.
  subroutine test_cpp_delta_eddington_solver_t(config_file_path)

    use tuvx_core,                     only: core_t
    use tuvx_grid,                     only: grid_t
    use tuvx_profile,                  only: profile_t
    use tuvx_solver,                   only: radiation_field_t
    use musica_string,                 only: string_t

    type(string_t), intent(in) :: config_file_path

    type(core_t),     pointer :: core
    class(grid_t),    pointer :: columns
    class(profile_t), pointer :: solar_zenith_angle ! [degrees]
    class(profile_t), pointer :: earth_sun_distance ! [AU]
    type(radiation_field_t), allocatable :: f90_radiation_fields(:),          &
                                            cpp_radiation_fields(:)
    integer :: i_column

    core => core_t(config_file_path)
    columns => core%get_grid( "time", "hours" )
    solar_zenith_angle => core%get_profile( "solar zenith angle", "degrees" )
    earth_sun_distance => core%get_profile( "Earth-Sun distance", "AU" )

    ! Run the solver for each set of conditions
    allocate( f90_radiation_fields( columns%ncells_ ) )
    do i_column = 1, columns%ncells_
      call core%run( solar_zenith_angle%edge_val_( i_column ),                &
                     earth_sun_distance%edge_val_( i_column ) )

      f90_radiation_fields( i_column ) = core%get_radiation_field( )
    end do
    cpp_radiation_fields =                                                    &
        calculate_cpp_radiation_fields( core,                                 &
                                        solar_zenith_angle%edge_val_,         &
                                        earth_sun_distance%edge_val_ )
    call compare_radiation_fields( f90_radiation_fields, cpp_radiation_fields )

    ! Clean up
    deallocate( earth_sun_distance )
    deallocate( solar_zenith_angle )
    deallocate( columns )
    deallocate( core )
    
  end subroutine test_cpp_delta_eddington_solver_t

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Calculates the radiation field using the C++ Delta-Eddington solver
  function calculate_cpp_radiation_fields( tuvx_core, solar_zenith_angles,      &
      earth_sun_distances ) result( radiation_fields )

    use tuvx_constants,                only: pi
    use tuvx_core,                     only: core_t
    use tuvx_grid,                     only: grid_t
    use tuvx_profile,                  only: profile_t
    use tuvx_solver,                   only: radiation_field_t

    type(core_t),     intent(in) :: tuvx_core
    real(dk),         intent(in) :: solar_zenith_angles(:)
    real(dk),         intent(in) :: earth_sun_distances(:)
    type(radiation_field_t), allocatable :: radiation_fields(:)

    class(grid_t), pointer :: heights
    class(grid_t), pointer :: wavelengths
    real(kind=c_double), allocatable, target :: solar_zenith_angles_c(:)
    real(kind=c_double), allocatable, target :: earth_sun_distances_c(:)
    real(kind=c_double), allocatable, target :: altitude_mid_points_c(:,:)
    real(kind=c_double), allocatable, target :: altitude_edges_c(:,:)
    real(kind=c_double), allocatable, target :: wavelength_mid_points_c(:)
    real(kind=c_double), allocatable, target :: wavelength_edges_c(:)
    type(solver_input_t_c) :: input
    type(solver_output_t_c) :: output
    integer :: i_column, n_lev, n_wl

    ! Set up the input struct
    heights => tuvx_core%get_grid( "height", "km" )
    wavelengths => tuvx_core%get_grid( "wavelength", "nm" )
    input%n_wavelengths_ = wavelengths%ncells_
    input%n_levels_      = heights%ncells_
    input%n_columns_     = size( solar_zenith_angles )
    solar_zenith_angles_c   = real( solar_zenith_angles(:), kind=c_double )   &
                                  * real(pi, kind=c_double) / 180.0_c_double ! degrees -> radians
    earth_sun_distances_c   = real( earth_sun_distances(:), kind=c_double )
    allocate( altitude_mid_points_c( size( solar_zenith_angles ),             &
                                     heights%ncells_ ) )
    allocate( altitude_edges_c(      size( solar_zenith_angles ),             &
                                     heights%ncells_+1 ) )
    do i_column = 1, size( solar_zenith_angles )
      altitude_mid_points_c(i_column,:) =                                     &
          real( heights%mid_(:), kind=c_double ) * 1.0e3_c_double ! km -> m
      altitude_edges_c(i_column,:)      =                                     &
          real( heights%edge_(:), kind=c_double ) * 1.0e3_c_double ! km -> m
    end do
    wavelength_mid_points_c = real( wavelengths%mid_(:), kind=c_double )      &
                                    * 1.0e-9_c_double ! nm -> m
    wavelength_edges_c      = real( wavelengths%edge_(:), kind=c_double )     &
                                    * 1.0e-9_c_double ! nm -> m
    input%solar_zenith_angles_   = c_loc( solar_zenith_angles_c )
    input%earth_sun_distances_   = c_loc( earth_sun_distances_c )
    input%altitude_mid_points_   = c_loc( altitude_mid_points_c )
    input%altitude_edges_        = c_loc( altitude_edges_c )
    input%wavelength_mid_points_ = c_loc( wavelength_mid_points_c )
    input%wavelength_edges_      = c_loc( wavelength_edges_c )

    ! run the C++ Delta-Eddington solver
    output = run_delta_eddington_solver_c( input )

    ! copy output to radiation_fields
    allocate( radiation_fields( size( solar_zenith_angles ) ) )
    do i_column = 1, size( solar_zenith_angles )
      n_lev = heights%ncells_
      n_wl = wavelengths%ncells_
      allocate( radiation_fields(i_column)%edr_( n_lev+1, n_wl ) )
      allocate( radiation_fields(i_column)%eup_( n_lev+1, n_wl ) )
      allocate( radiation_fields(i_column)%edn_( n_lev+1, n_wl ) )
      allocate( radiation_fields(i_column)%fdr_( n_lev+1, n_wl ) )
      allocate( radiation_fields(i_column)%fup_( n_lev+1, n_wl ) )
      allocate( radiation_fields(i_column)%fdn_( n_lev+1, n_wl ) )
      call copy_c_array_to_fortran( output%irrad_direct_,                      &
                               radiation_fields(i_column)%edr_,                &
                               i_column, heights%ncells_, wavelengths%ncells_ )
      call copy_c_array_to_fortran( output%irrad_up_,                          &
                               radiation_fields(i_column)%eup_,                &
                               i_column, heights%ncells_, wavelengths%ncells_ )
      call copy_c_array_to_fortran( output%irrad_down_,                        &
                               radiation_fields(i_column)%edn_,                &
                               i_column, heights%ncells_, wavelengths%ncells_ )
      call copy_c_array_to_fortran( output%flux_direct_,                       &
                               radiation_fields(i_column)%fdr_,                &
                               i_column, heights%ncells_, wavelengths%ncells_ )
      call copy_c_array_to_fortran( output%flux_up_,                           &
                               radiation_fields(i_column)%fup_,                &
                               i_column, heights%ncells_, wavelengths%ncells_ )
      call copy_c_array_to_fortran( output%flux_down_,                         &
                               radiation_fields(i_column)%fdn_,                &
                               i_column, heights%ncells_, wavelengths%ncells_ )
    end do

    call free_output_c( output )
    deallocate( heights )
    deallocate( wavelengths )

  end function calculate_cpp_radiation_fields

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Compares the radiation fields calculated by the Fortran and C++ solvers
  subroutine compare_radiation_fields( f90_radiation_fields, cpp_radiation_fields )

    use tuvx_solver, only: radiation_field_t

    type(radiation_field_t), intent(in) :: f90_radiation_fields(:)
    type(radiation_field_t), intent(in) :: cpp_radiation_fields(:)

    ! Compare the radiation fields
    ! ...

  end subroutine compare_radiation_fields

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Copies a 3D C array pointer to a Fortran array
  subroutine copy_c_array_to_fortran( c_array_ptr, f_array, i_col, n_lev, n_wl )

    type(c_ptr), intent(in) :: c_array_ptr
    real(dk), intent(inout) :: f_array(:,:)
    integer, intent(in) :: i_col, n_lev, n_wl

    real(kind=c_double), pointer :: c_array(:,:,:)

    call c_f_pointer( c_array_ptr, c_array, [i_col, n_lev+1, n_wl] )
    f_array(:,:) = real( c_array(i_col,:,:), kind=dk )

  end subroutine copy_c_array_to_fortran

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program test_cpp_delta_eddington_solver