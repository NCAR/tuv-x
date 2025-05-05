! Copyright (C) 2023-2024 National Science Foundation-National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
program test_cpp_delta_eddington_solver

#define ASSERT(x) call assert(x, __FILE__, __LINE__)

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
    type(c_ptr) :: layer_OD_              ! (columns, levels, wavelengths)
    type(c_ptr) :: layer_SSA_             ! (columns, levels, wavelengths)
    type(c_ptr) :: layer_G_               ! (columns, levels, wavelengths)
  end type solver_input_t_c

  ! Struct to hold output radiation fields from the C++ Delta-Eddington solver
  type, bind(c) :: solver_output_t_c
    integer(c_int) :: n_wavelengths_
    integer(c_int) :: n_levels_
    integer(c_int) :: n_columns_
    type(c_ptr) :: flux_direct_  ! (columns, levels+1, wavelengths)
    type(c_ptr) :: flux_up_      ! (columns, levels+1, wavelengths)
    type(c_ptr) :: flux_down_    ! (columns, levels+1, wavelengths)
    type(c_ptr) :: irrad_direct_ ! (columns, levels+1, wavelengths)
    type(c_ptr) :: irrad_up_     ! (columns, levels+1, wavelengths)
    type(c_ptr) :: irrad_down_   ! (columns, levels+1, wavelengths)
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
    use tuvx_radiator,                 only: radiator_state_t
    use tuvx_solver,                   only: radiation_field_t
    use musica_string,                 only: string_t

    type(string_t), intent(in) :: config_file_path

    type(core_t),     pointer :: core
    class(grid_t),    pointer :: columns, heights, wavelengths
    class(profile_t), pointer :: solar_zenith_angle ! [degrees]
    class(profile_t), pointer :: earth_sun_distance ! [AU]
    type(radiation_field_t), allocatable :: f90_radiation_fields(:),          &
                                            cpp_radiation_fields(:)
    type(radiator_state_t), allocatable :: radiator_states(:)
    integer :: i_col

    core => core_t(config_file_path)
    columns => core%get_grid( "time", "hours" )
    heights => core%get_grid( "height", "km" )
    wavelengths => core%get_grid( "wavelength", "nm" )
    solar_zenith_angle => core%get_profile( "solar zenith angle", "degrees" )
    earth_sun_distance => core%get_profile( "Earth-Sun distance", "AU" )

    ! Run the solver for each set of conditions
    allocate( f90_radiation_fields( columns%ncells_ + 1 ) )
    allocate( radiator_states( columns%ncells_ + 1 ) )
    do i_col = 1, columns%ncells_ + 1
      call core%run( solar_zenith_angle%edge_val_( i_col ),                   &
                     earth_sun_distance%edge_val_( i_col ) )

      f90_radiation_fields( i_col ) = core%get_radiation_field( )
      allocate( radiator_states( i_col )%layer_G_( heights%ncells_,           &
                                                      wavelengths%ncells_, 1 ) )
      call core%radiative_transfer_%radiator_warehouse_%accumulate_states(    &
                                                  radiator_states( i_col ) )
    end do
    cpp_radiation_fields =                                                    &
        calculate_cpp_radiation_fields( core,                                 &
                                        solar_zenith_angle%edge_val_,         &
                                        earth_sun_distance%edge_val_,         &
                                        radiator_states )
    call compare_radiation_fields( f90_radiation_fields, cpp_radiation_fields )

    ! Clean up
    deallocate( earth_sun_distance )
    deallocate( solar_zenith_angle )
    deallocate( columns )
    deallocate( heights )
    deallocate( wavelengths )
    deallocate( core )
    
  end subroutine test_cpp_delta_eddington_solver_t

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Calculates the radiation field using the C++ Delta-Eddington solver
  function calculate_cpp_radiation_fields( tuvx_core, solar_zenith_angles,      &
      earth_sun_distances, radiator_states ) result( radiation_fields )

    use tuvx_constants,                only: pi
    use tuvx_core,                     only: core_t
    use tuvx_grid,                     only: grid_t
    use tuvx_profile,                  only: profile_t
    use tuvx_radiator,                 only: radiator_state_t
    use tuvx_solver,                   only: radiation_field_t

    type(core_t),             intent(in) :: tuvx_core
    real(dk),                 intent(in) :: solar_zenith_angles(:)
    real(dk),                 intent(in) :: earth_sun_distances(:)
    type(radiator_state_t),   intent(in) :: radiator_states(:)
    type(radiation_field_t), allocatable :: radiation_fields(:)

    class(grid_t), pointer :: heights
    class(grid_t), pointer :: wavelengths
    real(kind=c_double), allocatable, target :: solar_zenith_angles_c(:)
    real(kind=c_double), allocatable, target :: earth_sun_distances_c(:)
    real(kind=c_double), allocatable, target :: altitude_mid_points_c(:,:)
    real(kind=c_double), allocatable, target :: altitude_edges_c(:,:)
    real(kind=c_double), allocatable, target :: wavelength_mid_points_c(:)
    real(kind=c_double), allocatable, target :: wavelength_edges_c(:)
    real(kind=c_double), allocatable, target :: layer_OD_c(:,:,:)
    real(kind=c_double), allocatable, target :: layer_SSA_c(:,:,:)
    real(kind=c_double), allocatable, target :: layer_G_c(:,:,:)
    type(solver_input_t_c) :: input
    type(solver_output_t_c) :: output
    integer :: i_col, n_col, n_lev, n_wl

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
    allocate( layer_OD_c(            size( solar_zenith_angles ),             &
                                     heights%ncells_, wavelengths%ncells_ ) )
    allocate( layer_SSA_c(           size( solar_zenith_angles ),             &
                                     heights%ncells_, wavelengths%ncells_ ) )
    allocate( layer_G_c(             size( solar_zenith_angles ),             &
                                     heights%ncells_, wavelengths%ncells_ ) )
    do i_col = 1, size( solar_zenith_angles )
      altitude_mid_points_c(i_col,:) =                                        &
          real( heights%mid_(:), kind=c_double ) * 1.0e3_c_double ! km -> m
      altitude_edges_c(i_col,:)      =                                        &
          real( heights%edge_(:), kind=c_double ) * 1.0e3_c_double ! km -> m
      layer_OD_c(i_col,:,:)          =                                        &
          real( radiator_states(i_col)%layer_OD_(:,:), kind=c_double )
      layer_SSA_c(i_col,:,:)         =                                        &
          real( radiator_states(i_col)%layer_SSA_(:,:), kind=c_double )
      layer_G_c(i_col,:,:)           =                                        &
          real( radiator_states(i_col)%layer_G_(:,:,1), kind=c_double )
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
    input%layer_OD_              = c_loc( layer_OD_c )
    input%layer_SSA_             = c_loc( layer_SSA_c )
    input%layer_G_               = c_loc( layer_G_c )

    ! run the C++ Delta-Eddington solver
    output = run_delta_eddington_solver_c( input )

    ! copy output to radiation_fields
    allocate( radiation_fields( output%n_columns_ ) )
    n_col = output%n_columns_
    do i_col = 1, n_col
      n_lev = heights%ncells_
      n_wl = wavelengths%ncells_
      allocate( radiation_fields(i_col)%edr_( n_lev+1, n_wl ) )
      allocate( radiation_fields(i_col)%eup_( n_lev+1, n_wl ) )
      allocate( radiation_fields(i_col)%edn_( n_lev+1, n_wl ) )
      allocate( radiation_fields(i_col)%fdr_( n_lev+1, n_wl ) )
      allocate( radiation_fields(i_col)%fup_( n_lev+1, n_wl ) )
      allocate( radiation_fields(i_col)%fdn_( n_lev+1, n_wl ) )
      call copy_c_array_to_fortran( output%irrad_direct_,                      &
                               radiation_fields(i_col)%edr_,                   &
                               i_col, n_col, n_lev, n_wl )
      call copy_c_array_to_fortran( output%irrad_up_,                          &
                               radiation_fields(i_col)%eup_,                   &
                               i_col, n_col, n_lev, n_wl )
      call copy_c_array_to_fortran( output%irrad_down_,                        &
                               radiation_fields(i_col)%edn_,                   &
                               i_col, n_col, n_lev, n_wl )
      call copy_c_array_to_fortran( output%flux_direct_,                       &
                               radiation_fields(i_col)%fdr_,                   &
                               i_col, n_col, n_lev, n_wl )
      call copy_c_array_to_fortran( output%flux_up_,                           &
                               radiation_fields(i_col)%fup_,                   &
                               i_col, n_col, n_lev, n_wl )
      call copy_c_array_to_fortran( output%flux_down_,                         &
                               radiation_fields(i_col)%fdn_,                   &
                               i_col, n_col, n_lev, n_wl )
    end do

    call free_output_c( output )
    deallocate( heights )
    deallocate( wavelengths )

  end function calculate_cpp_radiation_fields

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Compares two double values for closeness
  function is_close( a, b )

    real(dk), intent(in) :: a, b
    logical :: is_close

    is_close = abs( a - b ) <= (abs( a ) + abs( b )) * 1.0e-6

  end function is_close

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Asserts that a condition is true
  subroutine assert( condition, file, line )

    logical, intent(in) :: condition
    character(len=*), intent(in) :: file
    integer, intent(in) :: line

    if ( .not. condition ) then
      write(*,*) 'Assertion failed at line ', line, ' in ', file
      stop 3
    end if

  end subroutine assert

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Compares the radiation fields calculated by the Fortran and C++ solvers
  subroutine compare_radiation_fields( f90_radiation_fields, cpp_radiation_fields )

    use tuvx_solver, only: radiation_field_t

    type(radiation_field_t), intent(in) :: f90_radiation_fields(:)
    type(radiation_field_t), intent(in) :: cpp_radiation_fields(:)

    integer :: i_col, i_lev, i_wl
    integer :: n_col, n_lev, n_wl
    logical :: is_val_close

    ! [DEV NOTES] Temporarily check for the fixed output of the C++ solver
    n_col = size( f90_radiation_fields )
    n_lev = size( f90_radiation_fields(1)%edr_, 1 )
    n_wl = size( f90_radiation_fields(1)%edr_, 2 )
    do i_wl = 1, n_wl
      do i_lev = 1, n_lev
        do i_col = 1, n_col
          is_val_close =                                                     &
              is_close( cpp_radiation_fields(i_col)%edr_(i_lev,i_wl),        &
                   real( (i_wl-1) * n_lev * n_col + (i_lev-1) * n_col +      &
                         (i_col-1) + 42, kind=dk ))
          ASSERT( is_val_close )
          is_val_close =                                                     &
              is_close( cpp_radiation_fields(i_col)%eup_(i_lev,i_wl),        &
                   real( (i_wl-1) * n_lev * n_col + (i_lev-1) * n_col +      &
                         (i_col-1) + 93, kind=dk ))
          ASSERT( is_val_close )
          is_val_close =                                                     &
              is_close( cpp_radiation_fields(i_col)%edn_(i_lev,i_wl),        &
                   real( (i_wl-1) * n_lev * n_col + (i_lev-1) * n_col +      &
                         (i_col-1) + 52, kind=dk ))
          ASSERT( is_val_close )
          is_val_close =                                                     &
              is_close( cpp_radiation_fields(i_col)%fdr_(i_lev,i_wl),        &
                   real( (i_wl-1) * n_lev * n_col + (i_lev-1) * n_col +      &
                         (i_col-1) + 5, kind=dk ))
          ASSERT( is_val_close )
          is_val_close =                                                     &
              is_close( cpp_radiation_fields(i_col)%fup_(i_lev,i_wl),        &
                   real( (i_wl-1) * n_lev * n_col + (i_lev-1) * n_col +      &
                         (i_col-1) + 24, kind=dk ))
          ASSERT( is_val_close )
          is_val_close =                                                     &
              is_close( cpp_radiation_fields(i_col)%fdn_(i_lev,i_wl),        &
                   real( (i_wl-1) * n_lev * n_col + (i_lev-1) * n_col +      &
                         (i_col-1) + 97, kind=dk ))
          ASSERT( is_val_close )
        end do
      end do
    end do

  end subroutine compare_radiation_fields

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Copies a 3D C array pointer to a Fortran array
  subroutine copy_c_array_to_fortran( c_array_ptr, f_array, i_col, n_col,     &
      n_lev, n_wl )

    type(c_ptr), intent(in) :: c_array_ptr
    real(dk), intent(inout) :: f_array(:,:)
    integer, intent(in) :: i_col, n_col, n_lev, n_wl

    real(kind=c_double), pointer :: c_array(:,:,:)

    call c_f_pointer( c_array_ptr, c_array, [n_col, n_lev+1, n_wl] )
    f_array(:,:) = real( c_array(i_col,:,:), kind=dk )

  end subroutine copy_c_array_to_fortran

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program test_cpp_delta_eddington_solver