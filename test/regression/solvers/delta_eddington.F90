! Copyright (C) 2023-2024 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
program test_cpp_delta_eddington_solver

  use musica_constants,                only: dk => musica_dk
  use musica_string,                   only: string_t

  implicit none

  type(string_t) :: config_file_path

  ! Run the TUV 5.4 test
  config_file_path = 'examples/tuv_5_4.json'
  call test_cpp_delta_eddington_solver_t(config_file_path)

  ! Run the TS1/TSMLT test
  config_file_path = 'examples/ts1_tsmlt.json'
  call test_cpp_delta_eddington_solver_t(config_file_path)

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
    class(grid_t),    pointer :: time
    class(profile_t), pointer :: solar_zenith_angle ! [degrees]
    class(profile_t), pointer :: earth_sun_distance ! [AU]
    type(radiation_field_t), allocatable :: f90_radiation_fields(:),          &
                                            cpp_radiation_fields(:)
    integer :: i_time

    core => core_t(config_file_path)
    time => core%get_grid( "time", "hours" )
    solar_zenith_angle => core%get_profile( "solar zenith angle", "degrees" )
    earth_sun_distance => core%get_profile( "Earth-Sun distance", "AU" )

    ! Run the solver for each set of conditions
    allocate( f90_radiation_fields( time%ncells_ ) )
    do i_time = 1, time%ncells_
      call core%run( solar_zenith_angle%edge_val_( i_time ),                  &
                     earth_sun_distance%edge_val_( i_time ) )

      f90_radiation_fields( i_time ) = core%get_radiation_field( )
    end do
    cpp_radiation_fields =                                                    &
        calculate_cpp_radiation_fields( core,                                 &
                                        solar_zenith_angle%edge_val_,         &
                                        earth_sun_distance%edge_val_ )
    call compare_radiation_fields( f90_radiation_fields, cpp_radiation_fields )

    ! Clean up
    deallocate( earth_sun_distance )
    deallocate( solar_zenith_angle )
    deallocate( time )
    deallocate( core )
    
  end subroutine test_cpp_delta_eddington_solver_t

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Calculates the radiation field using the C++ Delta-Eddington solver
  function calculate_cpp_radiation_fields( tuvx_core, solar_zenith_angle,      &
      earth_sun_distance) result( radiation_fields )

    use tuvx_core,                     only: core_t
    use tuvx_profile,                  only: profile_t
    use tuvx_solver,                   only: radiation_field_t

    type(core_t),     intent(in) :: tuvx_core
    real(dk),         intent(in) :: solar_zenith_angle(:)
    real(dk),         intent(in) :: earth_sun_distance(:)
    type(radiation_field_t), allocatable :: radiation_fields(:)

    ! run the C++ Delta-Eddington solver
    allocate( radiation_fields( size( solar_zenith_angle ) ) )

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

end program test_cpp_delta_eddington_solver