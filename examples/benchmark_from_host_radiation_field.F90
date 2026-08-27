! Copyright (C) 2026 University Corporation for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
program benchmark_from_host_radiation_field
  ! Times the TS1-TSMLT photolysis rate constant calculation with and
  ! without an internal radiative transfer solve, and separates the two
  ! costs.
  !
  ! This program runs the benchmark once for each of two radiative
  ! transfer solvers: delta-Eddington, in ``examples/ts1_tsmlt.json``, and
  ! the discrete ordinate solver (4 streams), in
  ! ``examples/ts1_tsmlt_discrete_ordinate.json``. The two solvers differ
  ! greatly in cost, so running both shows how much of a full call the
  ! radiative transfer solve is, and how much the ``from host`` solver in
  ! ``examples/ts1_tsmlt_from_host.json`` saves, for each.
  !
  ! For each solver, this program calls core_t%run() three ways, each
  ! averaged over several iterations to reduce timing noise:
  !
  ! 1. The radiative transfer solve alone, with no photolysis rate
  !    constants requested.
  ! 2. The same solve, together with the photolysis rate constant
  !    calculation.
  ! 3. The photolysis rate constant calculation alone, with the ``from
  !    host`` solver, fed the radiation field that run 2 calculated.
  !    TUV-x does no radiative transfer solve in this run.
  !
  ! Run 2 minus run 1 is an independent estimate of the photolysis rate
  ! constant calculation cost, and should be close to run 3. Comparing run
  ! 1 to run 2 shows what share of a full call the radiative transfer solve
  ! actually is. See ``docs/source/host_radiation_field.rst`` for the full
  ! description of the ``from host`` solver.
  !
  ! Run this program from a directory that contains the ``examples`` and
  ! ``data`` directories, for example the CMake build directory.

  use musica_constants,              only : dk => musica_dk
  use musica_string,                 only : string_t
  use tuvx_core,                     only : core_t
  use tuvx_grid,                     only : grid_t
  use tuvx_solver,                   only : radiation_field_t
  use tuvx_solver_from_host,         only : radiation_field_updater_t

  implicit none

  character(len=*), parameter :: kReuseConfigPath =                          &
      "examples/ts1_tsmlt_from_host.json"
  integer,          parameter :: kIterations = 20 ! calls averaged per timing

  ! The solar zenith angle and the Earth-Sun distance are the same for
  ! every run. The Earth-Sun distance must be 1.0 AU. TUV-x applies the
  ! Earth-Sun distance factor to the radiation field on every call to
  ! core_t%run(). Run 2 applies the factor once, while it builds the
  ! radiation field. Run 3 applies the factor again, to the field that this
  ! program copied from run 2. A factor of 1.0 leaves the field unchanged,
  ! so run 3 stays consistent with run 2.
  real(dk), parameter :: kSolarZenithAngle = 40.0_dk ! [degrees]
  real(dk), parameter :: kEarthSunDistance = 1.0_dk  ! [AU]

  integer(8) :: count_rate

  call system_clock( count_rate = count_rate )

  call run_benchmark( "Delta-Eddington", "examples/ts1_tsmlt.json" )
  call run_benchmark( "Discrete ordinate (4 streams)",                       &
                      "examples/ts1_tsmlt_discrete_ordinate.json" )

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine run_benchmark( label, solve_config_path )
    ! Times one radiative transfer solver against the ``from host`` solver

    character(len=*), intent(in) :: label             ! solver name, for the report
    character(len=*), intent(in) :: solve_config_path  ! path to its TUV-x configuration

    type(string_t)                    :: config_path
    class(core_t),           pointer  :: solve_core, reuse_core
    class(grid_t),           pointer  :: height
    type(radiation_field_t)           :: field
    type(radiation_field_updater_t)   :: updater
    logical                           :: found
    integer                           :: n_levels, n_reactions, i_iteration
    real(dk), allocatable             :: solve_rates(:,:), reuse_rates(:,:)
    integer(8) :: start_count, end_count
    real(dk)   :: solve_only_seconds, solve_and_photolysis_seconds,          &
                 reuse_seconds, implied_photolysis_seconds
    real(dk)   :: max_relative_difference

    write(*,'(A)') ""
    write(*,'(A)') "=== " // trim( label ) // " ==="

    ! Set up the solver run
    config_path = solve_config_path
    solve_core => core_t( config_path )

    height => solve_core%get_grid( "height", "km" )
    n_levels = height%ncells_ + 1
    deallocate( height )
    n_reactions = solve_core%number_of_photolysis_reactions( )
    allocate( solve_rates( n_levels, n_reactions ) )
    allocate( reuse_rates( n_levels, n_reactions ) )

    ! Run 1: the radiative transfer solve alone. photolysis_rate_constants
    ! is absent, so TUV-x skips the photolysis rate constant calculation.
    call system_clock( start_count )
    do i_iteration = 1, kIterations
      call solve_core%run( solar_zenith_angle = kSolarZenithAngle,           &
                           earth_sun_distance = kEarthSunDistance )
    end do
    call system_clock( end_count )
    solve_only_seconds = real( end_count - start_count, dk ) /               &
        real( count_rate, dk ) / real( kIterations, dk )

    ! Run 2: the radiative transfer solve together with the photolysis rate
    ! constant calculation
    call system_clock( start_count )
    do i_iteration = 1, kIterations
      call solve_core%run( solar_zenith_angle = kSolarZenithAngle,           &
                           earth_sun_distance = kEarthSunDistance,           &
                           photolysis_rate_constants = solve_rates )
    end do
    call system_clock( end_count )
    solve_and_photolysis_seconds = real( end_count - start_count, dk ) /     &
        real( count_rate, dk ) / real( kIterations, dk )

    ! Take a copy of the radiation field from the last iteration of run 2
    field = solve_core%get_radiation_field( )
    deallocate( solve_core )

    ! Set up the ``from host`` run
    config_path = kReuseConfigPath
    reuse_core => core_t( config_path )

    updater = reuse_core%get_radiation_field_updater( found )
    if( .not. found ) then
      write(*,*) "examples/ts1_tsmlt_from_host.json must set the "//        &
                 "radiative transfer solver type to 'from host'."
      stop 3
    end if

    ! Run 3: the photolysis rate constant calculation alone, with a reused
    ! radiation field. A host application must call update() before every
    ! call to run(), so that call is inside the timing loop too, even
    ! though it is only an array copy.
    call system_clock( start_count )
    do i_iteration = 1, kIterations
      call updater%update( direct_actinic_flux   = field%fdr_,               &
                           upward_actinic_flux   = field%fup_,               &
                           downward_actinic_flux = field%fdn_ )
      call reuse_core%run( solar_zenith_angle = kSolarZenithAngle,           &
                           earth_sun_distance = kEarthSunDistance,           &
                           photolysis_rate_constants = reuse_rates )
    end do
    call system_clock( end_count )
    reuse_seconds = real( end_count - start_count, dk ) /                    &
        real( count_rate, dk ) / real( kIterations, dk )
    deallocate( reuse_core )

    implied_photolysis_seconds = solve_and_photolysis_seconds -              &
        solve_only_seconds

    ! The two sets of rate constants should agree. This is a check on the
    ! ``from host`` setup, not a timing result.
    max_relative_difference = maxval( abs( reuse_rates - solve_rates ) /     &
        max( maxval( abs( solve_rates ) ), tiny( 1.0_dk ) ) )

    write(*,'(A,I0,A,I0,A,I0,A)') "Photolysis rate constants: ", n_reactions,&
        " reactions, ", n_levels, " vertical levels, averaged over ",       &
        kIterations, " calls"
    write(*,'(A,F12.6,A)') "1. Radiative transfer solve alone:            ", &
        solve_only_seconds, " s"
    write(*,'(A,F12.6,A)') "2. Radiative transfer solve + photolysis:     ", &
        solve_and_photolysis_seconds, " s"
    write(*,'(A,F12.6,A)') "3. Photolysis alone (reused radiation field): ", &
        reuse_seconds, " s"
    write(*,'(A,F12.6,A)') "   (2 - 1), for comparison with 3:            ", &
        implied_photolysis_seconds, " s"
    write(*,'(A,ES12.4)') "Maximum relative difference between the two "//   &
        "sets of rate constants: ", max_relative_difference

  end subroutine run_benchmark

end program benchmark_from_host_radiation_field
